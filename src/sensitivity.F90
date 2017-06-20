module sensitivity

#ifdef MPI
  use message_passing
#endif

  use algorithm,   only: binary_search
  use constants,   only: ONE
  use error,       only: fatal_error, warning
  use global
  use mesh,        only: count_bank_sites_bin, get_mesh_bin, &
                         mesh_indices_to_bin, mesh_intersects_2d, &
                         mesh_intersects_3d, get_mesh_indices, &
                         bin_to_mesh_indices
  use message_passing
  use particle_header,  only: Particle
  use sensitivity_header
  use fissionmatrix,    only: collect_fm, fm_eigenvalue, free_fm_memory, &
       count_source_for_fm
  use tally_filter

  implicit none
  integer :: position(N_FILTER_TYPES - 3) = 0 ! Tally map positioning array

contains

!===============================================================================
! SENSITIVITY_INITIALIZE_BATCH
!===============================================================================

  subroutine sensitivity_initialize_batch()

    if (adjointmethod == 1) then
      if (current_batch > n_inactive) call ifptally_reset_batch()
    end if
    if (adjointmethod == 2) then
      if (current_batch <= n_inactive) call ifptally_reset_batch()
      if (current_batch > n_inactive) call clutch_reset_batch()
    end if
    if (adjointmethod == 3) then   ! fission matrix approach
      if (current_batch > 10 .and. current_batch <= n_inactive) call count_source_for_fm()
      if (current_batch > n_inactive) call clutch_reset_batch()
    end if
    if (adjointmethod == 4) then
      if (current_batch > n_inactive) call ifptally_reset_batch()
    end if
    if (adjointmethod == 5) then
      if (current_batch <= n_inactive) call ifptally_reset_batch()
      if (current_batch > n_inactive) call clutch_reset_batch()
    end if
    if (adjointmethod == 6) then   ! fission matrix approach
      if (current_batch > 10 .and. current_batch <= n_inactive) call count_source_for_fm()
      if (current_batch > n_inactive) call clutch_reset_batch()
    end if

  end subroutine sensitivity_initialize_batch

!===============================================================================
! SENSITIVITY_FINALIZE_BATCH
!===============================================================================

  subroutine sensitivity_finalize_batch()

    if (adjointmethod == 1) then
      if (asymptotic) call sensitivity_calc()
    end if
    if (adjointmethod == 2) then
      if (asymptotic) call importance_calc()
      if (clutch_second) call sensitivity_calc_clutch()
    end if
    if (adjointmethod == 3) then
      if (current_batch == n_inactive) call importance_calc_fm()
      if (clutch_second) call sensitivity_calc_clutch()
    end if
    if (adjointmethod == 4) then
      if (current_batch == n_inactive) call reaction_rates()
      if (asymptotic) call sensitivity_gpt_calc()
    end if
    if (adjointmethod == 5) then
      if (asymptotic) call importance_gpt_calc()
      if (clutch_second) call sensitivity_calc_gclutch()
    end if
    if (adjointmethod == 6) then
      if (current_batch == n_inactive) call reaction_rates()
      if (current_batch == n_inactive) call importance_gpt_fm()
      if (clutch_second) call sensitivity_calc_gclutch()
    end if
    if (current_batch == n_max_batches) call sen_statistics()

  end subroutine sensitivity_finalize_batch

!===============================================================================
! SENSITIVITY_INITIALIZE_HISTORY
!===============================================================================

  subroutine sensitivity_initialize_history()

    if (adjointmethod == 1 .and. original) call tally_reset_history()
    if (adjointmethod == 2 .and. clutch_first) call tally_reset_history()
    if (adjointmethod == 3 .and. clutch_first) call tally_reset_history()
    if (adjointmethod == 4 .and. original) call tally_reset_history()
    if (adjointmethod == 5 .and. clutch_first) call tally_reset_history()
    if (adjointmethod == 6 .and. clutch_first) call tally_reset_history()

  end subroutine sensitivity_initialize_history

!===============================================================================
! SCORE_TRACKLENGTH_SENSITIVITY
!===============================================================================

  subroutine score_tracklength_sensitivity(p, distance)

    type(Particle), intent(in) :: p
    real(8),        intent(in) :: distance

    integer :: next_bin_mesh        ! index for next mesh
    integer :: next_bin_energy           ! index for energy bin

    integer :: i
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: i_nuclide            ! index in nuclides array (from bins)
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: atom_density         ! atom density of single nuclide in atom/b-cm
    real(8) :: mesh_weight          ! weight of mesh filter
    type(SensitivityObject), pointer :: t
    type(Material),    pointer :: mat

    ! Determine track-length estimate of flux
    flux = distance  ! w0 = 1, initial weight of particle

    atom_density = ZERO

    ! A loop over all sensitivities is necessary

    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      call get_next_bin_sen(t, p, &
           NO_BIN_FOUND, next_bin_mesh, mesh_weight)  ! get mesh bin

      if (next_bin_mesh == NO_BIN_FOUND) cycle SENSITIVITY_LOOP

      next_bin_energy = binary_search(t % energystructure, &   ! get energy bin
           t % n_energy_bins + 1, p % E)

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      MESH_LOOP: do

        ! ======================================================================
        ! Nuclide logic

        NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins
          ! Get index of nuclide in nuclides array
          i_nuclide = t % nuclide_bins(k)

          if (i_nuclide > 0) then
            if (p % material /= MATERIAL_VOID) then
              ! Get pointer to current material
              mat => materials(p % material)

              ! Determine if nuclide is actually in material
              NUCLIDE_MAT_LOOP: do j = 1, mat % n_nuclides
                ! If index of nuclide matches the j-th nuclide listed in the
                ! material, break out of the loop
                if (i_nuclide == mat % nuclide(j)) exit

                ! If we've reached the last nuclide in the material, it means
                ! the specified nuclide to be tallied is not in this material
                if (j == mat % n_nuclides) then
                  cycle NUCLIDE_BIN_LOOP
                end if
              end do NUCLIDE_MAT_LOOP

              atom_density = mat % atom_density(j)
            else
              atom_density = ZERO
            end if
          end if

          ! Determine score for each bin
          call score_track_general(p, t, k, next_bin_mesh, next_bin_energy, &
               i_nuclide, atom_density, flux * mesh_weight)

        end do NUCLIDE_BIN_LOOP

        ! ======================================================================
        ! Mesh logic

        call get_next_bin_sen(t, p, &
             next_bin_mesh, next_bin_mesh, mesh_weight)
        if (next_bin_mesh == NO_BIN_FOUND) exit MESH_LOOP

      end do MESH_LOOP

    end do SENSITIVITY_LOOP

  end subroutine score_tracklength_sensitivity

!===============================================================================
! SCORE_SCATTERING_SENSITIVITY
!===============================================================================

  subroutine score_scattering_sensitivity(p, i_nuclide, mt_number)

    type(Particle), intent(in) :: p
    integer, intent(in) :: i_nuclide
    integer, intent(in) :: mt_number
    integer :: i_nuclide_sen
    integer :: i_score_sen

    integer :: next_bin_mesh        ! index for next mesh
    integer :: next_bin_energy      ! index for energy bin
    integer :: next_bin_nuclide     ! index for nuclide bin
    integer :: next_bin_score       ! index for score bin

    integer :: i_mesh
    integer :: i
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    type(SensitivityObject), pointer :: t
    type(Material),    pointer :: mat
    type(RegularMesh), pointer :: m


    ! A loop over all sensitivities is necessary

    SENSITIVITY_LOOP: do i = 1, n_sens

      ! clear the information of last sensitivity
      next_bin_mesh = NO_BIN_FOUND
      next_bin_energy = 0
      next_bin_nuclide = 0
      next_bin_score = 0

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      if (senmesh_dict % has_key(t % meshid)) then
        i_mesh = senmesh_dict % get_key(t % meshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % meshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if

      ! ======================================================================
      ! Mesh logic
      call get_mesh_bin(m, p % coord(1) % xyz, next_bin_mesh)
      if (next_bin_mesh == NO_BIN_FOUND) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Energy logic
      next_bin_energy = binary_search(t % energystructure, &   ! get energy bin
           t % n_energy_bins + 1, p % last_E)
      if (next_bin_energy == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Nuclide logic,  to check if the collision nuclide is in the
      ! sensitivity tally list
      NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins

        ! Get index of nuclide in nuclides array
        i_nuclide_sen = t % nuclide_bins(k)

        ! the collision nuclide is in the sensitivity tally list
        if (i_nuclide_sen == i_nuclide) then

          next_bin_nuclide = k

          exit NUCLIDE_BIN_LOOP

        end if

      end do NUCLIDE_BIN_LOOP

      if (next_bin_nuclide == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Score logic,  to check if the score is in the
      ! sensitivity tally list
      SCORE_BIN_LOOP: do j = 1, t % n_score_bins

        ! determine what type of score bin
        i_score_sen = t % score_bins(j)

        ! the reaction type is in the sensitivity tally list
        if (i_score_sen == mt_number) then

          next_bin_score = j

          exit SCORE_BIN_LOOP

        end if

      end do SCORE_BIN_LOOP

      if (next_bin_score == 0) cycle SENSITIVITY_LOOP

      t % cumtally(next_bin_nuclide, next_bin_score, &
           next_bin_mesh, next_bin_energy) = &
           t % cumtally(next_bin_nuclide, next_bin_score, &
           next_bin_mesh, next_bin_energy) + 1

    end do SENSITIVITY_LOOP

  end subroutine score_scattering_sensitivity

!===============================================================================
! SCORE_FISSION_SENSITIVITY
!===============================================================================

  subroutine score_fission_sensitivity(p, b, i_nuclide, mt_number)

    ! progenitor is a global parameter

    type(Particle), intent(in) :: p
    type(Bank), intent(in) :: b
    integer, intent(in) :: i_nuclide
    integer, intent(in) :: mt_number
    integer :: i_nuclide_sen
    integer :: i_score_sen

    integer :: next_bin_mesh        ! index for next mesh
    integer :: next_bin_energy      ! index for energy bin
    integer :: next_bin_nuclide     ! index for nuclide bin
    integer :: next_bin_score       ! index for score bin

    integer :: i_mesh
    integer :: i
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    type(SensitivityObject), pointer :: t
    type(Material),    pointer :: mat
    type(RegularMesh), pointer :: m


    ! A loop over all sensitivities is necessary

    SENSITIVITY_LOOP: do i = 1, n_sens

      ! clear the information of last sensitivity
      next_bin_mesh = NO_BIN_FOUND
      next_bin_energy = 0
      next_bin_nuclide = 0
      next_bin_score = 0

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      if (senmesh_dict % has_key(t % meshid)) then
        i_mesh = senmesh_dict % get_key(t % meshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % meshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if

      ! ======================================================================
      ! Mesh logic
      call get_mesh_bin(m, p % coord(1) % xyz, next_bin_mesh)
      if (next_bin_mesh == NO_BIN_FOUND) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Energy logic
      if (mt_number == FISSION_CHI) then
        next_bin_energy = binary_search(t % energystructure, & ! get energybin
             t % n_energy_bins + 1, b % E)
      else
        next_bin_energy = binary_search(t % energystructure, & ! get energybin
             t % n_energy_bins + 1, p % last_E)
      end if
      if (next_bin_energy == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Nuclide logic,  to check if the collision nuclide is in the
      ! sensitivity tally list
      NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins

        ! Get index of nuclide in nuclides array
        i_nuclide_sen = t % nuclide_bins(k)

        ! the collision nuclide is in the sensitivity tally list
        if (i_nuclide_sen == i_nuclide) then

          next_bin_nuclide = k

          exit NUCLIDE_BIN_LOOP

        end if

      end do NUCLIDE_BIN_LOOP

      if (next_bin_nuclide == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Score logic,  to check if the score is in the
      ! sensitivity tally list
      SCORE_BIN_LOOP: do j = 1, t % n_score_bins

        ! determine what type of score bin
        i_score_sen = t % score_bins(j)

        ! the reaction type is in the sensitivity tally list
        if (i_score_sen == mt_number) then

          next_bin_score = j

          exit SCORE_BIN_LOOP

        end if

      end do SCORE_BIN_LOOP

      if (next_bin_score == 0) cycle SENSITIVITY_LOOP

      t % neutrontally(progenitornum, next_bin_nuclide, next_bin_score, &
           next_bin_mesh, next_bin_energy) = &
           t % neutrontally(progenitornum, next_bin_nuclide, next_bin_score, &
           next_bin_mesh, next_bin_energy) + 1

    end do SENSITIVITY_LOOP

  end subroutine score_fission_sensitivity

!===============================================================================
! SCORE_IMPORTANCE_FUNCTION (IFP-->CLUTCH)
!===============================================================================

  subroutine score_importance_dis(p)

    ! progenitor is a global parameter

    type(Particle), intent(in) :: p

    integer :: next_bin_mesh        ! index for next mesh
    integer :: i_mesh
    integer :: i
    type(SensitivityObject), pointer :: t
    type(RegularMesh), pointer :: m


    ! A loop over all sensitivities is necessary

    SENSITIVITY_LOOP: do i = 1, n_sens

      ! clear the information of last sensitivity
      next_bin_mesh = NO_BIN_FOUND

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      if (senmesh_dict % has_key(t % impmeshid)) then
        i_mesh = senmesh_dict % get_key(t % impmeshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % impmeshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if

      ! ======================================================================
      ! Mesh logic
      call get_mesh_bin(m, p % coord(1) % xyz, next_bin_mesh)
      if (next_bin_mesh == NO_BIN_FOUND) cycle SENSITIVITY_LOOP

      t % neutrontally(progenitornum, 1, 1, &
           next_bin_mesh, 1) = &
           t % neutrontally(progenitornum, 1, 1, &
           next_bin_mesh, 1) + 1

    end do SENSITIVITY_LOOP

  end subroutine score_importance_dis

!===============================================================================
! SCORE_DENOM_SENSITIVITY calculates the denominators of sensitivities
!===============================================================================

  subroutine score_denom_sensitivity()

    integer :: i

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      associate (t => sensitivities(i))
        t % neutronfission(progenitornum) = t % neutronfission(progenitornum) + 1
      end associate
    end do SENSITIVITY_LOOP

  end subroutine score_denom_sensitivity

!===============================================================================
! SCORE_FISSION_SITES calculates number of fission neutrons (IFP-->CLUTCH)
!===============================================================================

  subroutine score_fission_sites(p)

    type(Particle), intent(in) :: p
    integer :: next_bin_mesh        ! index for next mesh
    integer :: i_mesh
    integer :: i
    type(SensitivityObject), pointer :: t
    type(RegularMesh), pointer :: m

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens

      ! clear the information of last sensitivity
      next_bin_mesh = NO_BIN_FOUND

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      if (senmesh_dict % has_key(t % impmeshid)) then
        i_mesh = senmesh_dict % get_key(t % impmeshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % impmeshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if

      ! ======================================================================
      ! Mesh logic
      call get_mesh_bin(m, p % coord(1) % xyz, next_bin_mesh)
      if (next_bin_mesh == NO_BIN_FOUND) cycle SENSITIVITY_LOOP

      t % neutronfission(next_bin_mesh) = t % neutronfission(next_bin_mesh) + 1
    end do SENSITIVITY_LOOP

  end subroutine score_fission_sites

!===============================================================================
! SCORE_NEUTRON_VALUE calculates the number of progenies for a given ifp_id
!===============================================================================

  subroutine score_neutron_value(p, nu_born)

    type(Particle), intent(in) :: p
    integer, intent(in) :: nu_born  ! the number of fission neutrons

    type(SensitivityObject), pointer :: t

    integer :: i
    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      t % neutronvalue(p % ifp_id) = t % neutronvalue(p % ifp_id) + nu_born
    end do SENSITIVITY_LOOP

  end subroutine score_neutron_value

!===============================================================================
! ADD_BRANCH_SENSITIVITY
!===============================================================================

  subroutine add_branch_sensitivity()

    integer :: i

    progenitornum = progenitornum + 1 ! number of progenitor has been added

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      associate (t => sensitivities(i))
        if (t % method == 1) then
          t % neutrontally(progenitornum,:,:,:,:) = t % cumtally(:,:,:,:)
        end if
      end associate
    end do SENSITIVITY_LOOP

  end subroutine add_branch_sensitivity

!===============================================================================
! TALLY_CUMTOSECONDARY transfer the tally information from cumulative tally to
! secondary array.
!===============================================================================

  subroutine tally_cumtosecondary(p)

    type(Particle), intent(in) :: p
    integer :: i

    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      associate (t => sensitivities(i))
        t % secondtally(p % n_secondary,:,:,:,:) = t % cumtally(:,:,:,:)
      end associate
    end do SENSITIVITY_LOOP
  end subroutine tally_cumtosecondary

!===============================================================================
! TALLY_SECONDARYTOCUM transfer the tally information from secondary array to
! cumulative tally.
!===============================================================================

  subroutine tally_secondarytocum(p)

    type(Particle), intent(in) :: p
    type(SensitivityObject), pointer :: t
    integer :: i

    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      t % cumtally(:,:,:,:) = t % secondtally(p % n_secondary,:,:,:,:)
    end do SENSITIVITY_LOOP
  end subroutine tally_secondarytocum

!===============================================================================
! IFPTALLY_RESET_BATCH reset the neutrontally variables into zero at the beginning
! and the end of each block calculation
!===============================================================================

  subroutine ifptally_reset_batch()

    integer :: i
    integer :: j
    integer :: n_discard

    original = .false.
    asymptotic = .false.

    if (adjointmethod == 1 .or. adjointmethod ==4) n_discard = n_inactive
    if (adjointmethod == 2) n_discard = 0
    if (adjointmethod == 5) n_discard = 0!100

    if (mod(current_batch - n_discard, ifp_block) == 1) then
      original = .true.          !  original generation
      progenitornum = 3 * work_index(rank) ! the initial index
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        associate (t => sensitivities(i))
          t % neutrontally   = 0
          t % neutronvalue   = 0
          if (adjointmethod /= 4) then
            t % neutronfission = 0
          end if
          t % results(3,:,:,:,:) = 0 ! used as direct term in gpt ifp
          if (adjointmethod == 5) then ! calculate importance
            t % gptvaluenumer  = 0
            t % gptvaluedenom  = 0
          end if
        end associate
      end do SENSITIVITY_LOOP
    end if
    if (mod(current_batch - n_discard, ifp_block) == 0) then
      asymptotic = .true.          ! asymptotic generation
    end if

  end subroutine ifptally_reset_batch

!===============================================================================
! CLUTCH_RESET_BATCH reset some variables related to CLUTCH calculation
!===============================================================================

  subroutine clutch_reset_batch()

    integer :: i
    integer :: n_discard

    original = .false.
    asymptotic = .false.
    clutch_first = .false.   ! first batch in CLUTCH calculation
    clutch_second = .false.  ! second batch in CLUTCH calculation
    fismatrix_on = .false.   ! turn off the fission matrix calculation

    if (mod(current_batch - n_inactive, 2) == 1) then
      clutch_first = .true.          !  first batch in one CLUTCH block
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        sensitivities(i) % clutchsen = 0   ! every process has an array
        sensitivities(i) % denom = 0       ! every process has one value
      end do SENSITIVITY_LOOP
    end if
    if (mod(current_batch - n_inactive, 2) == 0) then
      clutch_second = .true.         !  second batch in one CLUTCH block
    end if

  end subroutine clutch_reset_batch

!===============================================================================
! TALLY_RESET_HISTORY reset the neutrontally variables into zero at the
! beginning of each particle history
!===============================================================================

  subroutine tally_reset_history()

    type(SensitivityObject), pointer :: t
    integer :: i

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      t % cumtally = 0
      t % secondtally = 0
      !t % secondtally(1:maxsecondnum,:,:,:,:) = 0
    end do SENSITIVITY_LOOP
    maxsecondnum  = 0

  end subroutine tally_reset_history

#ifdef MPI
!===============================================================================
! COLLECT_IFPCAL means collecting sensitivity tallies from different cores
! (focused on IFP method)
!===============================================================================
  subroutine collect_ifpcal()

    type(SensitivityObject), pointer :: t
    integer :: i         ! sensitivity loop index
    integer :: n1        ! total size of count variable neutrontally
    integer :: n2        ! total size of count variable neutronfission
    integer :: n3        ! total size of count variable neutronvalue
    real(8) :: dummy     ! temporary receive buffer for non-root reductions

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      if (t % method == 1) then
        n1 = 3 * n_particles * t % n_nuclide_bins * t % n_score_bins * &
             t % n_mesh_bins * t % n_energy_bins
        n2 = 3 * n_particles
        n3 = 3 * n_particles
      end if
      if (t % method == 2) then
        n1 = 3 * n_particles * t % imp_mesh_bins
        n2 = t % imp_mesh_bins
        n3 = 3 * n_particles
      end if
      ! collect values from all processors
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, t % neutrontally, n1, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(MPI_IN_PLACE, t % neutronfission, n2, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(MPI_IN_PLACE, t % neutronvalue, n3, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      else
        ! Receive buffer not significant at other processors
        call MPI_REDUCE(t % neutrontally, dummy, n1, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(t % neutronfission, dummy, n2, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(t % neutronvalue, dummy, n3, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
      end if
    end do SENSITIVITY_LOOP

  end subroutine collect_ifpcal
#endif

!===============================================================================
! SENSITIVITY_CALC combines the ifp tallies and neutron importances to calculate
! real sensitivities
!===============================================================================

  subroutine sensitivity_calc()

    integer :: i, j, k, l, m, n
    real(8) :: val

#ifdef MPI
    call collect_ifpcal()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        associate (t => sensitivities(i))
          t % n_realizations = t % n_realizations + 1
          t % denom = 0
          do j = 1, 3 * n_particles
            t % denom = t % denom + t % neutronfission(j) * t % neutronvalue(j)
          end do

          do n = 1, t % n_energy_bins
            do m = 1, t % n_mesh_bins
              do l = 1, t % n_score_bins
                do k = 1, t % n_nuclide_bins
                  val = sum(t % neutrontally(:,k,l,m,n) * t % neutronvalue(:)) / &
                       t % denom
                  t % results(1,k,l,m,n) = t % results(1,k,l,m,n) + val
                  t % results(2,k,l,m,n) = t % results(2,k,l,m,n) + val * val
                end do
              end do
            end do
          end do
        end associate
      end do SENSITIVITY_LOOP
    end if

  end subroutine sensitivity_calc

!===============================================================================
! IMPORTANCE_CALC calculates the importance function via IFP method
!===============================================================================

  subroutine importance_calc()

    type(SensitivityObject), pointer :: t
    integer :: i
    integer :: j
    integer :: k
    integer :: n

#ifdef MPI
    call collect_ifpcal()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        t % imp_realizations = t % imp_realizations + 1

        do k = 1, t % imp_mesh_bins
          if (t % neutronfission(k) /= 0) then
            t % importance(k) = t % importance(k) + &
                 sum(t % neutrontally(:,1,1,k,1) * t % neutronvalue(:))/&
                 t % neutronfission(k)
          else
            t % importance(k) = t % importance(k)
          end if
          if (current_batch == n_inactive) then
            t % importance(k) = t % importance(k) / t % imp_realizations
            !print *, t % importance(k)
          end if
        end do

      end do SENSITIVITY_LOOP
    end if

    if (current_batch ==  n_inactive) then
      ! distribute this importance distribution to different process
#ifdef MPI
      ! A loop over all sensitivities is necessary
      do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        n = t % imp_mesh_bins
        call MPI_BCAST(t % importance, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
      end do
#endif
    end if

  end subroutine importance_calc

!===============================================================================
! IMPORTANCE_CALC_FM calculates the importance function via fission matrix method
!===============================================================================

  subroutine importance_calc_fm()

    type(SensitivityObject), pointer :: t
    integer :: i
    integer :: j
    integer :: k
    integer :: n

#ifdef MPI
    call collect_fm()
#endif

    if (master) then
      call fm_eigenvalue(1)
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        t % importance = fismatrix % eigenvector
        !t % importance = ONE
        !do j = 1, t % imp_mesh_bins
        !print *, t % importance(j)
        !end do
      end do SENSITIVITY_LOOP
    end if

    call free_fm_memory()

    ! distribute this importance distribution to different process
#ifdef MPI
    ! A loop over all sensitivities is necessary
    do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      n = t % imp_mesh_bins
      call MPI_BCAST(t % importance, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    end do
#endif

  end subroutine importance_calc_fm

!===============================================================================
! SENSITIVITY_CLUTCH_SCACOL CALCULATES SCATTERING AND COLLISION TERMS IN CLUTCH
! METHOD
!===============================================================================

  subroutine sensitivity_clutch_scacol(p)

    type(Particle), intent(in) :: p
    type(RegularMesh), pointer :: m
    type(Material),    pointer :: mat
    type(SensitivityObject), pointer :: t

    integer :: i
    integer :: i_mesh
    integer :: imp_mesh_bin

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      if (senmesh_dict % has_key(t % impmeshid)) then
        i_mesh = senmesh_dict % get_key(t % impmeshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % impmeshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if

      ! Mesh logic to determine the importance function
      call get_mesh_bin(m, p % coord(1) % xyz, imp_mesh_bin)
      if (imp_mesh_bin == NO_BIN_FOUND) cycle SENSITIVITY_LOOP

      ! cumulate sensitivity at each collision point
      if (p % material /= MATERIAL_VOID) then
        mat => materials(p % material)
        t % clutchsen(:,:,:,:) = t % clutchsen(:,:,:,:) + &
             t % cumtally(:,:,:,:) * p % wgt * t % importance(imp_mesh_bin) * &
             material_xs % nu_fission / material_xs % total
        ! t % clutchsen(:,:,:,:) = t % clutchsen(:,:,:,:) + &
        ! t % cumtally(:,:,:,:) * p % wgt * 1 * &
        !      material_xs % nu_fission / material_xs % total
      end if

    end do SENSITIVITY_LOOP

  end subroutine sensitivity_clutch_scacol

!===============================================================================
! SENSITIVITY_CLUTCH_FISSION calculates fission related sensitivities
!===============================================================================

  subroutine sensitivity_clutch_fission(p, mt_number)

    ! progenitor is a global parameter

    type(Particle), intent(in) :: p
    integer, intent(in) :: mt_number
    integer :: i_nuclide_sen
    integer :: i_score_sen
    integer :: i_mesh

    integer :: energy_bin      ! index for energy bin
    integer :: score_bin       ! index for score bin
    integer :: nuclide_bin     ! index for nuclide bin
    integer :: imp_mesh_bin    ! index for importance mesh bin

    integer :: i
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    type(SensitivityObject), pointer :: t
    type(Material),    pointer :: mat
    type(RegularMesh), pointer :: m


    ! A loop over all sensitivities is necessary

    SENSITIVITY_LOOP: do i = 1, n_sens

      ! clear the information of last sensitivity
      score_bin = 0
      energy_bin = 0
      nuclide_bin = 0
      imp_mesh_bin = 0

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      ! ======================================================================
      ! Mesh logic to determine the importance function
      if (senmesh_dict % has_key(t % impmeshid)) then
        i_mesh = senmesh_dict % get_key(t % impmeshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % impmeshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if
      call get_mesh_bin(m, p % coord(1) % xyz, imp_mesh_bin)
      if (imp_mesh_bin == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Energy logic
      if (mt_number == FISSION_CHI) then
        energy_bin = p % energy_born
      else
        energy_bin = p % energy_fission
      end if
      if (energy_bin == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Nuclide logic,  to check if the collision nuclide is in the
      ! sensitivity tally list
      NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins

        ! Get index of nuclide in nuclides array
        i_nuclide_sen = t % nuclide_bins(k)

        ! the collision nuclide is in the sensitivity tally list
        if (i_nuclide_sen == p % nuclide_born) then

          nuclide_bin = k

          exit NUCLIDE_BIN_LOOP

        end if

      end do NUCLIDE_BIN_LOOP
      if (nuclide_bin == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Score logic,  to check if the score is in the
      ! sensitivity tally list
      SCORE_BIN_LOOP: do j = 1, t % n_score_bins

        ! determine what type of score bin
        i_score_sen = t % score_bins(j)

        ! the reaction type is in the sensitivity tally list
        if (i_score_sen == mt_number) then

          score_bin = j

          exit SCORE_BIN_LOOP

        end if

      end do SCORE_BIN_LOOP

      if (score_bin == 0) cycle SENSITIVITY_LOOP

      t % clutchsen(nuclide_bin,score_bin,p % mesh_born,energy_bin) = &
           t % clutchsen(nuclide_bin,score_bin,p % mesh_born,energy_bin) + &
           p % wgt * t % importance(imp_mesh_bin) * material_xs % nu_fission / &
           material_xs % total

    end do SENSITIVITY_LOOP

  end subroutine sensitivity_clutch_fission

!===============================================================================
! transfer_clutch_info
!===============================================================================

  subroutine transfer_clutch_info(p, b, i_nuclide, mt_number)

    type(Particle), intent(in) :: p
    type(Bank), intent(inout) :: b
    integer, intent(in) :: i_nuclide
    integer, intent(in) :: mt_number

    integer :: i_mesh
    integer :: i
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    type(SensitivityObject), pointer :: t
    type(RegularMesh), pointer :: m


    ! A loop over all sensitivities isn't necessary, need to be modified

    SENSITIVITY_LOOP: do i = 1, n_sens

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      ! ======================================================================
      ! Nuclide_born logic
      b % nuclide_born = i_nuclide

      ! ======================================================================
      ! Energy_fission logic
      b % energy_fission = binary_search(t % energystructure, &
           t % n_energy_bins + 1, p % last_E)

      ! ======================================================================
      ! Energy_born logic
      b % energy_born = binary_search(t % energystructure, &
           t % n_energy_bins + 1, b % E)

      ! ======================================================================
      ! MT number logic
      b % mtnum_born = mt_number

      ! ======================================================================
      ! Mesh_born logic
      if (senmesh_dict % has_key(t % meshid)) then
        i_mesh = senmesh_dict % get_key(t % meshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % meshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if

      call get_mesh_bin(m, p % coord(1) % xyz, b % mesh_born)

    end do SENSITIVITY_LOOP

  end subroutine transfer_clutch_info

!===============================================================================
! SCORE_DENOM_CLUTCH calculates denominators in clutch calculations
!===============================================================================

  subroutine score_denom_clutch(p)

    type(Particle), intent(in) :: p
    type(RegularMesh), pointer :: m
    type(SensitivityObject), pointer :: t

    integer :: i
    integer :: i_mesh
    integer :: imp_mesh_bin    ! index for importance mesh bin

    ! A loop over all sensitivities isn't necessary, need to be modified
    SENSITIVITY_LOOP: do i = 1, n_sens

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      ! ======================================================================
      ! Mesh logic to determine the importance function
      if (senmesh_dict % has_key(t % impmeshid)) then
        i_mesh = senmesh_dict % get_key(t % impmeshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % impmeshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if
      call get_mesh_bin(m, p % coord(1) % xyz, imp_mesh_bin)

      if (imp_mesh_bin == 0) cycle SENSITIVITY_LOOP

      t % denom = t % denom + &
           p % wgt * t % importance(imp_mesh_bin) * material_xs % nu_fission / &
           material_xs % total
    end do SENSITIVITY_LOOP

  end subroutine score_denom_clutch

!===============================================================================
! SENSITIVITY_CLUTCH calculates sensitivities in clutch
!===============================================================================

  subroutine sensitivity_clutch(p)
    type(Particle), intent(in) :: p

    if (clutch_first) call sensitivity_clutch_scacol(p)

    if (clutch_second) then
      call sensitivity_clutch_fission(p, SCORE_FISSION)
      call sensitivity_clutch_fission(p, FISSION_NUBAR)
      call sensitivity_clutch_fission(p, FISSION_CHI)
      call sensitivity_clutch_fission(p, p % mtnum_born)
      call score_denom_clutch(p)
    end if

  end subroutine sensitivity_clutch

#ifdef MPI
!===============================================================================
! COLLECT_CLUTCHCAL means collecting sensitivity tallies from different cores
! (focused on CLUTCH method)
!===============================================================================

  subroutine collect_clutchcal()

    type(SensitivityObject), pointer :: t
    integer :: i         ! sensitivity loop index
    integer :: n         ! total size of count variable neutrontally
    real(8) :: dummy     ! temporary receive buffer for non-root reductions

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      n = t % n_nuclide_bins * t % n_score_bins * &
           t % n_mesh_bins * t % n_energy_bins
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, t % clutchsen, n, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(MPI_IN_PLACE, t % denom, 1, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      else
        ! Receive buffer not significant at other processors
        call MPI_REDUCE(t % clutchsen, dummy, n, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(t % denom, dummy, 1, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
      end if
    end do SENSITIVITY_LOOP

  end subroutine collect_clutchcal
#endif

!===============================================================================
! SENSITIVITY_CALC_CLUTCH calculates sensitivities at every second clutch batches
!===============================================================================

  subroutine sensitivity_calc_clutch()

    type(SensitivityObject), pointer :: t
    integer :: i
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    real(8) :: value

#ifdef MPI
    call collect_clutchcal()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        t % n_realizations = t % n_realizations + 1

        do k = 1, t % n_nuclide_bins
          do l = 1, t % n_score_bins
            do m = 1, t % n_mesh_bins
              do n = 1, t % n_energy_bins
                value = t%clutchsen(k,l,m,n)/t % denom
                t%results(1,k,l,m,n) = t%results(1,k,l,m,n) + value
                t%results(2,k,l,m,n) = t%results(2,k,l,m,n) + value *value
              end do
            end do
          end do
        end do

      end do SENSITIVITY_LOOP
    end if

  end subroutine sensitivity_calc_clutch

!===============================================================================
! SEN_STATISTICS computes the mean and standard deviation of the mean of each
! tally and stores them in the val and val_sq attributes of the TallyResults
! respectively
!===============================================================================

  subroutine sen_statistics()

    integer :: i    ! index in tallies array
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    type(SensitivityObject), pointer :: t

    ! Calculate statistics for user-defined tallies
    do i = 1, n_sens
      t => sensitivities(i)
      do k = 1, t % n_nuclide_bins
        do l = 1, t % n_score_bins
          do m = 1, t % n_mesh_bins
            do n = 1, t % n_energy_bins
              t%results(1,k,l,m,n) = t%results(1,k,l,m,n)/t%n_realizations
              t%results(2,k,l,m,n) = sqrt((t%results(2,k,l,m,n)/t%n_realizations - &
                   t%results(1,k,l,m,n)*t%results(1,k,l,m,n))/(t%n_realizations-1))
            end do
          end do
        end do
      end do
    end do

  end subroutine sen_statistics

!===============================================================================
! GET_NEXT_BIN_SEN gives the index for the next valid filter bin and a weight
! that will be applied to the flux, this subroutine will be called when using
! tracklength estimator to tally in sensitivity calculations
!===============================================================================

  subroutine get_next_bin_sen(t, p, current_bin, next_bin, weight)
    type(SensitivityObject),    intent(in)  :: t
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: current_bin
    integer,           intent(out) :: next_bin
    real(8),           intent(out) :: weight

    integer, parameter :: MAX_SEARCH_ITER = 100 ! Maximum number of times we can
    !  can loop while trying to find
    !  the first intersection.

    integer :: i_mesh
    integer :: j                    ! loop index for direction
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: search_iter          ! loop count for intersection search
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross            ! coordinates of next boundary
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: total_distance       ! distance of entire particle track
    real(8) :: distance             ! distance traveled in mesh cell
    logical :: start_in_mesh        ! starting coordinates inside mesh?
    logical :: end_in_mesh          ! ending coordinates inside mesh?
    type(RegularMesh), pointer :: m


    if (senmesh_dict % has_key(t % meshid)) then
      i_mesh = senmesh_dict % get_key(t % meshid)
      m => sen_meshes(i_mesh)
    else
      call fatal_error("Could not find mesh " // trim(to_str(t % meshid)) &
           // " specified on sensitivity " // trim(to_str(t % id)))
    end if

    ! A track can span multiple mesh bins so we need to handle a lot of
    ! intersection logic for tracklength tallies.

    ! Copy the starting and ending coordinates of the particle.  Offset these
    ! just a bit for the purposes of determining if there was an intersection
    ! in case the mesh surfaces coincide with lattice/geometric surfaces which
    ! might produce finite-precision errors.
    xyz0 = p % last_xyz + TINY_BIT * p % coord(1) % uvw
    xyz1 = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw

    ! Determine indices for starting and ending location.
    call get_mesh_indices(m, xyz0, ijk0(:m % n_dimension), start_in_mesh)
    call get_mesh_indices(m, xyz1, ijk1(:m % n_dimension), end_in_mesh)

    ! If this is the first iteration of the filter loop, check if the track
    ! intersects any part of the mesh.
    if (current_bin == NO_BIN_FOUND) then
      if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
        if (m % n_dimension == 2) then
          if (.not. mesh_intersects_2d(m, xyz0, xyz1)) then
            next_bin = NO_BIN_FOUND
            return
          end if
        else
          if (.not. mesh_intersects_3d(m, xyz0, xyz1)) then
            next_bin = NO_BIN_FOUND
            return
          end if
        end if
      end if
    end if

    ! Copy the un-modified coordinates the particle direction.
    xyz0 = p % last_xyz
    xyz1 = p % coord(1) % xyz
    uvw = p % coord(1) % uvw

    ! Compute the length of the entire track.
    total_distance = sqrt(sum((xyz1 - xyz0)**2))

    ! If we're looking for the first valid bin, check to see if the particle
    ! starts inside the mesh.
    if (current_bin == NO_BIN_FOUND) then
      if (any(ijk0(:m % n_dimension) < 1) &
           .or. any(ijk0(:m % n_dimension) > m % dimension)) then

        ! The particle does not start in the mesh so keep iterating the ijk0
        ! indices to cross the nearest mesh surface until we've found a valid
        ! bin.  MAX_SEARCH_ITER prevents an infinite loop.
        search_iter = 0
        do while (any(ijk0(:m % n_dimension) < 1) &
             .or. any(ijk0(:m % n_dimension) > m % dimension))
          if (search_iter == MAX_SEARCH_ITER) then
            call warning("Failed to find a mesh intersection on a tally mesh &
                 &filter.")
            next_bin = NO_BIN_FOUND
            return
          end if

          do j = 1, m % n_dimension
            if (abs(uvw(j)) < FP_PRECISION) then
              d(j) = INFINITY
            else if (uvw(j) > 0) then
              xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
              d(j) = (xyz_cross - xyz0(j)) / uvw(j)
            else
              xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
              d(j) = (xyz_cross - xyz0(j)) / uvw(j)
            end if
          end do
          j = minloc(d(:m % n_dimension), 1)
          if (uvw(j) > ZERO) then
            ijk0(j) = ijk0(j) + 1
          else
            ijk0(j) = ijk0(j) - 1
          end if
          search_iter = search_iter + 1
        end do
        distance = d(j)
        xyz0 = xyz0 + distance * uvw

      end if
    end if

    ! ========================================================================
    ! If we've already scored some mesh bins, figure out which mesh cell is
    ! next and where the particle enters that cell.

    if (current_bin /= NO_BIN_FOUND) then
      ! Get the indices to the last bin.
      call bin_to_mesh_indices(m, current_bin, ijk0(:m % n_dimension))

      ! If the particle track ends in that bin, then we are done.
      if (all(ijk0(:m % n_dimension) == ijk1(:m % n_dimension))) then
        next_bin = NO_BIN_FOUND
        return
      end if

      ! Figure out which face of the previous mesh cell our track exits, i.e.
      ! the closest surface of that cell for which
      ! dot(p % uvw, face_normal) > 0.
      do j = 1, m % n_dimension
        if (abs(uvw(j)) < FP_PRECISION) then
          d(j) = INFINITY
        else if (uvw(j) > 0) then
          xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
          d(j) = (xyz_cross - xyz0(j)) / uvw(j)
        else
          xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
          d(j) = (xyz_cross - xyz0(j)) / uvw(j)
        end if
      end do
      j = minloc(d(:m % n_dimension), 1)

      ! Translate the starting coordintes by the distance to that face. This
      ! should be the xyz that we computed the distance to in the last
      ! iteration of the filter loop.
      distance = d(j)
      xyz0 = xyz0 + distance * uvw

      ! Increment the indices into the next mesh cell.
      if (uvw(j) > ZERO) then
        ijk0(j) = ijk0(j) + 1
      else
        ijk0(j) = ijk0(j) - 1
      end if

      ! If the next indices are invalid, then the track has left the mesh and
      ! we are done.
      if (any(ijk0(:m % n_dimension) < 1) &
           .or. any(ijk0(:m % n_dimension) > m % dimension)) then
        next_bin = NO_BIN_FOUND
        return
      end if
    end if

    ! Compute the length of the track segment in this mesh cell.
    if (all(ijk0(:m % n_dimension) == ijk1(:m % n_dimension))) then
      ! The track ends in this cell.  Use the particle end location rather
      ! than the mesh surface.
      distance = sqrt(sum((xyz1 - xyz0)**2))
    else
      ! The track exits this cell.  Use the distance to the mesh surface.
      do j = 1, m % n_dimension
        if (abs(uvw(j)) < FP_PRECISION) then
          d(j) = INFINITY
        else if (uvw(j) > 0) then
          xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
          d(j) = (xyz_cross - xyz0(j)) / uvw(j)
        else
          xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
          d(j) = (xyz_cross - xyz0(j)) / uvw(j)
        end if
      end do
      distance = minval(d(:m % n_dimension))
    end if

    ! Assign the next tally bin and the score
    next_bin = mesh_indices_to_bin(m, ijk0(:m % n_dimension))
    weight = distance / total_distance

  end subroutine get_next_bin_sen

!===============================================================================
! SCORE_TRACK_GENERAL adds scores to the sensitivity tally
!===============================================================================

  subroutine score_track_general(p, t, bin_nuclide, bin_mesh, &
       bin_energy, i_nuclide, atom_density, flux)
    type(Particle),    intent(in)    :: p
    type(SensitivityObject), intent(inout) :: t
    integer,           intent(in)    :: bin_nuclide
    integer,           intent(in)    :: bin_mesh
    integer,           intent(in)    :: bin_energy
    integer,           intent(in)    :: i_nuclide
    real(8),           intent(in)    :: flux           ! flux estimate
    real(8),           intent(in)    :: atom_density   ! atom/b-cm

    integer :: i                    ! loop index for scoring bins
    integer :: l                    ! loop index for nuclides in material
    integer :: m                    ! loop index for reactions
    integer :: q                    ! loop index for scoring bins
    integer :: i_temp               ! temperature index
    integer :: i_nuc                ! index in nuclides array (from material)
    integer :: i_energy             ! index in nuclide energy grid
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    integer :: d                    ! delayed neutron index
    integer :: g                    ! delayed neutron index
    integer :: k                    ! loop index for bank sites
    real(8) :: f                    ! interpolation factor
    real(8) :: score                ! analog tally score
    real(8) :: E                    ! particle energy

    i = 0
    SCORE_LOOP: do q = 1, t % n_score_bins
      i = i + 1

      ! determine what type of score bin
      score_bin = t % score_bins(i)

      !#########################################################################
      ! Determine appropirate scoring value.

      select case(score_bin)

      case (SCORE_TOTAL)
        score = micro_xs(i_nuclide) % total * atom_density * flux

      case (SCORE_SCATTER)
        score = (micro_xs(i_nuclide) % total &
             - micro_xs(i_nuclide) % absorption) * atom_density * flux

      case (SCORE_ABSORPTION)
        score = micro_xs(i_nuclide) % absorption * atom_density * flux

      case (SCORE_FISSION)
        score = micro_xs(i_nuclide) % fission * atom_density * flux

      case (SCORE_CAPTURE)
        score = (micro_xs(i_nuclide) % absorption &
             - micro_xs(i_nuclide) % fission) * atom_density * flux

      case (ELASTIC)
        score = micro_xs(i_nuclide) % elastic * atom_density * flux

      case (FISSION_CHI)
        score = ZERO

      case (FISSION_NUBAR)
        score = ZERO

      case default
        ! Any other cross section has to be calculated on-the-fly. For
        ! cross sections that are used often (e.g. n2n, ngamma, etc. for
        ! depletion), it might make sense to optimize this section or
        ! pre-calculate cross sections
        if (score_bin > 1) then
          ! Set default score
          score = ZERO

          if (i_nuclide > 0) then
            if (nuclides(i_nuclide)%reaction_index%has_key(score_bin)) then
              m = nuclides(i_nuclide)%reaction_index%get_key(score_bin)

              ! Retrieve temperature and energy grid index and interpolation
              ! factor
              i_temp = micro_xs(i_nuclide) % index_temp
              i_energy = micro_xs(i_nuclide) % index_grid
              f = micro_xs(i_nuclide) % interp_factor

              associate (xs => nuclides(i_nuclide) % reactions(m) % xs(i_temp))
                if (i_energy >= xs % threshold) then
                  score = ((ONE - f) * xs % value(i_energy - &
                       xs % threshold + 1) + f * xs % value(i_energy - &
                       xs % threshold + 2)) * atom_density * flux
                end if
              end associate
            end if
          end if

        else
          call fatal_error("Invalid score type on sensitivity " &
               // to_str(t % id) // ".")
        end if

      end select

      t % cumtally(bin_nuclide, i, bin_mesh, bin_energy) = &
           t % cumtally(bin_nuclide, i, bin_mesh, bin_energy) - score

    end do SCORE_LOOP

  end subroutine score_track_general

!===============================================================================
!============================GPT CALCULATION BELOW==============================
!===============================================================================


!===============================================================================
! REACTION_RATES calculates the averaged reaction rates in inactive batches
!===============================================================================

  subroutine reaction_rates()

    integer :: i_sen       ! index in sensitivities array
    integer :: n
    integer :: j
    type(SensitivityObject), pointer :: s
    real(8) :: dummy     ! temporary receive buffer for non-root reductions


    ! need parallel collection
    do i_sen = 1, n_sens
      s => sensitivities(i_sen)

#ifdef MPI
      ! collect values from all processors
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, s % respnumer, 1, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(MPI_IN_PLACE, s % respdenom, 1, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        if (adjointmethod == 6) then
          n = s % imp_mesh_bins
          call MPI_REDUCE(MPI_IN_PLACE, s % tallynumer, n, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(MPI_IN_PLACE, s % tallydenom, n, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        end if
      else
        ! Receive buffer not significant at other processors
        call MPI_REDUCE(s % respnumer, dummy, 1, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(s % respdenom, dummy, 1, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
        if (adjointmethod == 6) then
          n = s % imp_mesh_bins
          call MPI_REDUCE(s % tallynumer, dummy, n, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(s % tallydenom, dummy, n, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
        end if
      end if
#endif

      if (master) then
        s % respnumer = s % respnumer / n_inactive
        s % respdenom = s % respdenom / n_inactive
        if (adjointmethod == 6) then
          ! mesh_born_fm is not correct in first generation, so discard it
          s % tallynumer = s % tallynumer / (n_inactive-1)
          s % tallydenom = s % tallydenom / (n_inactive-1)
        end if
      end if

#ifdef MPI
      call MPI_BCAST(s % respnumer, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
      call MPI_BCAST(s % respdenom, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
      if (adjointmethod == 6) then
        n = s % imp_mesh_bins
        call MPI_BCAST(s % tallynumer, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(s % tallydenom, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
      end if
#endif

    end do

  end subroutine reaction_rates

!===============================================================================
! SCORE_INACTIVE calculates different tallies in inactive batches <Sigma1*Flux>,
! <Sigma2*FLux> and so on.
!===============================================================================

  subroutine score_inactive(p)

    type(Particle), intent(in) :: p
    integer :: i_sen       ! index in sensitivities array
    integer :: i_response  ! index in responses array
    integer :: i_tally     ! index in tallies array

    type(SensitivityObject), pointer :: s
    type(ResponseObject), pointer :: r
    type(TallyObject), pointer :: t

    call score_tally(p)  ! get all tallies at one collision point

    do i_sen = 1, n_sens
      s => sensitivities(i_sen)
      ! decide whether sensitivity has a corresponding response
      if (response_dict % has_key(s % response)) then
        i_response = response_dict % get_key(s % response)
        r => responses(i_response)
        i_tally = resptally_dict % get_key(r % numerid)
        s % respnumer = s % respnumer + resptalresult(i_tally)
        ! mesh_born_fm is not correct in first generation, so discard it
        if (adjointmethod == 6 .and. current_batch > 1) then
          s % tallynumer(p % mesh_born_fm) = s % tallynumer(p % mesh_born_fm)+&
               resptalresult(i_tally)
        end if
        i_tally = resptally_dict % get_key(r % denomid)
        s % respdenom = s % respdenom + resptalresult(i_tally)
        if (adjointmethod == 6 .and. current_batch > 1) then
          s % tallydenom(p % mesh_born_fm) = s % tallydenom(p % mesh_born_fm)+&
               resptalresult(i_tally)
        end if
      end if
    end do

  end subroutine score_inactive

!===============================================================================
! SCORE_DIRECTANDINTRA calculates direct sensitivity and intragenerational
! sensitivity term.
!===============================================================================

  subroutine score_directandintra(p)

    type(Particle), intent(in) :: p
    integer :: i_sen         ! index in sensitivities array
    integer :: i_response    ! index in responses array
    integer :: i_tally       ! index in tallies array
    integer :: i_mesh        ! index in meshes array
    real(8) :: gptimportance ! for each response, calculate one importance

    type(SensitivityObject), pointer :: s
    type(ResponseObject), pointer :: r
    type(TallyObject), pointer :: t
    type(RegularMesh), pointer :: m

    integer :: k
    integer :: j
    integer :: bin_mesh
    integer :: bin_energy
    integer :: bin_nuclide_numer
    integer :: bin_score_numer
    integer :: bin_nuclide_denom
    integer :: bin_score_denom

    call score_tally(p)  ! get all tally values at one collision point

    SENSITIVITY_LOOP: do i_sen = 1, n_sens
      s => sensitivities(i_sen)
      i_response = response_dict % get_key(s % response)
      r => responses(i_response)

      ! ======================================================================
      ! calculate the intragenerational term
      i_tally = resptally_dict % get_key(r % numerid)
      gptimportance = resptalresult(i_tally) / s % respnumer
      i_tally = resptally_dict % get_key(r % denomid)
      gptimportance = gptimportance - resptalresult(i_tally) / s % respdenom
      s % results(3,:,:,:,:) = s % results(3,:,:,:,:) + &
           s % cumtally(:,:,:,:) * gptimportance

      ! ======================================================================
      ! Mesh logic
      bin_mesh = 0
      i_mesh = senmesh_dict % get_key(s % meshid)
      m => sen_meshes(i_mesh)
      call get_mesh_bin(m, p % coord(1) % xyz, bin_mesh)

      ! ======================================================================
      ! Energy logic
      bin_energy = 0
      bin_energy = binary_search(s % energystructure, &   ! get energy bin
           s % n_energy_bins + 1, p % E)

      ! ======================================================================
      ! decide whether there is a direct sensitivity term
      i_tally = resptally_dict % get_key(r % numerid)
      t => resp_tallies(i_tally)
      call tally_in_senlist(s,t,bin_nuclide_numer,bin_score_numer)

      if (bin_nuclide_numer /=0 .and. bin_score_numer /=0 .and. &
           bin_mesh /=0 .and. bin_energy /=0) then
        s % results(3,bin_nuclide_numer,bin_score_numer,bin_mesh, bin_energy) = &
             s % results(3,bin_nuclide_numer,bin_score_numer,bin_mesh, bin_energy) + &
             resptalresult(i_tally) / s % respnumer
      end if


      i_tally = resptally_dict % get_key(r % denomid)
      t => resp_tallies(i_tally)
      call tally_in_senlist(s,t,bin_nuclide_denom,bin_score_denom)

      if (bin_nuclide_denom /=0 .and. bin_score_denom /=0 .and. &
           bin_mesh /=0 .and. bin_energy /=0) then
        s % results(3,bin_nuclide_denom,bin_score_denom,bin_mesh, bin_energy) = &
             s % results(3,bin_nuclide_denom,bin_score_denom,bin_mesh, bin_energy) - &
             resptalresult(i_tally) / s % respdenom
      end if

    end do SENSITIVITY_LOOP

  end subroutine score_directandintra

!===============================================================================
! SCORE_GPT_VALUE calculates generalized importance at every collision point.
!===============================================================================

  subroutine score_gpt_value(p)

    type(Particle), intent(in) :: p
    type(SensitivityObject), pointer :: s
    type(ResponseObject), pointer :: r
    integer :: i_sen
    integer :: i_response
    integer :: i_tally
    real(8) :: gptimportance

    call score_tally(p)  ! get all tally values at one collision point

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i_sen = 1, n_sens

      ! Get index of tally and pointer to tally
      s => sensitivities(i_sen)
      i_response = response_dict % get_key(s % response)
      r => responses(i_response)

      ! ======================================================================
      ! calculate the intergenerational importance
      if (adjointmethod == 4) then
        i_tally = resptally_dict % get_key(r % numerid)
        gptimportance = resptalresult(i_tally) / s % respnumer
        i_tally = resptally_dict % get_key(r % denomid)
        gptimportance = gptimportance - resptalresult(i_tally) / s % respdenom
        s % neutronvalue(p % ifp_id) = s % neutronvalue(p % ifp_id) + gptimportance
      else if (adjointmethod == 5) then
        ! need this gptvaluenumer to calculate s % respdenom in inactive batches
        i_tally = resptally_dict % get_key(r % numerid)
        s % gptvaluenumer(p % ifp_id) = s % gptvaluenumer(p % ifp_id) + &
             resptalresult(i_tally)
        i_tally = resptally_dict % get_key(r % denomid)
        s % gptvaluedenom(p % ifp_id) = s % gptvaluedenom(p % ifp_id) + &
             resptalresult(i_tally)
      end if

    end do SENSITIVITY_LOOP

  end subroutine score_gpt_value

#ifdef MPI
!===============================================================================
! COLLECT_GPTIFPCAL means collecting sensitivity tallies from different cores
! (focused on IFP method)
!===============================================================================

  subroutine collect_gptifpcal()

    type(SensitivityObject), pointer :: s
    integer :: i         ! sensitivity loop index
    integer :: n1        ! total size of count variable neutrontally
    integer :: n2        ! total size of count variable neutronfission
    integer :: n3        ! total size of count variable neutronvalue
    real(8) :: dummy     ! temporary receive buffer for non-root reductions

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      s => sensitivities(i)
      if (s % method == 4) then
        n1 = 3 * n_particles * s % n_nuclide_bins * s % n_score_bins * &
             s % n_mesh_bins * s % n_energy_bins
        n2 = 3 * n_particles
        n3 = s % n_nuclide_bins * s % n_score_bins * &
             s % n_mesh_bins * s % n_energy_bins
        ! collect values from all processors
        if (master) then
          call MPI_REDUCE(MPI_IN_PLACE, s % neutrontally, n1, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(MPI_IN_PLACE, s % neutronvalue, n2, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(MPI_IN_PLACE, s % results(3,:,:,:,:), n3, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(s % neutrontally, dummy, n1, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(s % neutronvalue, dummy, n2, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(s % results(3,:,:,:,:), dummy, n3, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
        end if
      end if
      if (s % method == 5) then
        n1 = 3 * n_particles * s % imp_mesh_bins
        n2 = 3 * n_particles
        n3 = s % imp_mesh_bins
        ! collect values from all processors
        if (master) then
          call MPI_REDUCE(MPI_IN_PLACE, s % neutrontally, n1, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(MPI_IN_PLACE, s % gptvaluenumer, n2, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(MPI_IN_PLACE, s % gptvaluedenom, n2, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(MPI_IN_PLACE, s % neutronfission, n3, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(s % neutrontally, dummy, n1, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(s % gptvaluenumer, dummy, n2, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(s % gptvaluedenom, dummy, n2, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
          call MPI_REDUCE(s % neutronfission, dummy, n3, MPI_REAL8, MPI_SUM, 0, &
               MPI_COMM_WORLD, mpi_err)
        end if
      end if

    end do SENSITIVITY_LOOP

  end subroutine collect_gptifpcal
#endif

!===============================================================================
! SENSITIVITY_GPT_CALC combines the ifp tallies and neutron importances to calculate
! real sensitivities
!===============================================================================

  subroutine sensitivity_gpt_calc()

    type(SensitivityObject), pointer :: s
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    real(8) :: value

#ifdef MPI
    call collect_gptifpcal()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        s => sensitivities(i)
        s % n_realizations = s % n_realizations + 1

        do k = 1, s % n_nuclide_bins
          do l = 1, s % n_score_bins
            do m = 1, s % n_mesh_bins
              do n = 1, s % n_energy_bins
                value = sum(s%neutrontally(:,k,l,m,n)*s%neutronvalue(:))
                value = value + s%results(3,k,l,m,n) ! add direct and intra
                !value = s%results(3,k,l,m,n)
                s%results(1,k,l,m,n) = s%results(1,k,l,m,n) + value
                s%results(2,k,l,m,n) = s%results(2,k,l,m,n) + value * value
              end do
            end do
          end do
        end do
      end do SENSITIVITY_LOOP
    end if

  end subroutine sensitivity_gpt_calc

!===============================================================================
! IMPORTANCE_GPT_CALC calculates the GPT importance function via IFP method
!===============================================================================

  subroutine importance_gpt_calc()

    type(SensitivityObject), pointer :: t
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    real(8) :: tallynumer
    real(8) :: tallydenom

    asymptotic = .false.

#ifdef MPI
    call collect_gptifpcal()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        t % imp_realizations = t % imp_realizations + 1
        tallynumer = sum(t % gptvaluenumer) / (ifp_block-1) ! cumulative tallies
        tallydenom = sum(t % gptvaluedenom) / (ifp_block-1) ! negative values
        t % respnumer = t % respnumer + tallynumer
        t % respdenom = t % respdenom + tallydenom

        do k = 1, t % imp_mesh_bins
          if (t % neutronfission(k) /= 0) then
            t % importance(k) = t % importance(k) + &
                 sum(t % neutrontally(:,1,1,k,1) * &
                 (t%gptvaluenumer(:)/tallynumer-t%gptvaluedenom(:)/tallydenom))/&
                 t % neutronfission(k)
          else
            t % importance(k) = t % importance(k)
          end if
          if (current_batch == n_inactive) then
            t % importance(k) = t % importance(k) / t % imp_realizations
          end if
        end do

        if (current_batch == n_inactive) then
          t % respnumer = t % respnumer / t % imp_realizations
          t % respdenom = t % respdenom / t % imp_realizations
        end if

      end do SENSITIVITY_LOOP
    end if

    if (current_batch == n_inactive) then
      ! distribute this importance distribution to different process
#ifdef MPI
      ! A loop over all sensitivities is necessary
      do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        n = t % imp_mesh_bins
        call MPI_BCAST(t % importance, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(t % respnumer,  1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(t % respdenom,  1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
      end do
#endif
    end if

  end subroutine importance_gpt_calc

!===============================================================================
! IMPORTANCE_GPT_FM calculates the importance function via fission matrix method
!===============================================================================

  subroutine importance_gpt_fm()

    type(SensitivityObject), pointer :: t
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    real(8), allocatable :: gptsource(:)
    allocate(gptsource(fismatrix % fm_dimension))

#ifdef MPI
    call collect_fm()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        gptsource = t%tallynumer/t%respnumer-t%tallydenom/t%respdenom
        call fm_fixedsource(gptsource,t % importance)
        !do j = 1, t % imp_mesh_bins
        !   print *, t % importance(j)
        !end do
      end do SENSITIVITY_LOOP
    end if

    call free_fm_memory()


    ! distribute this importance distribution to different process
#ifdef MPI
    ! A loop over all sensitivities is necessary
    do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      n = t % imp_mesh_bins
      call MPI_BCAST(t % importance, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    end do
#endif

    if(allocated(gptsource)) deallocate (gptsource)

  end subroutine importance_gpt_fm

!===============================================================================
! FM_FIXEDSOURCE means calculation of the generalized function
!===============================================================================

  subroutine fm_fixedsource(gptsource,gptflux)

    real(8), intent(in) :: gptsource(fismatrix % fm_dimension)
    real(8), intent(out) :: gptflux(fismatrix % fm_dimension)

    real(8), allocatable :: tmp1(:)
    real(8), allocatable :: forwardflux(:)
    real(8), allocatable :: adjointflux(:)
    real(8) :: eigenvalue


    real(8) :: normfa       ! norm of adjointT*forward
    real(8) :: norm         ! norm of gptfluxT*forward
    integer :: i            ! loop index i
    integer :: j            ! loop index j
    integer :: k            ! loop index k
    integer :: si           ! loop index for source iteration

    allocate(tmp1(fismatrix % fm_dimension))
    allocate(forwardflux(fismatrix % fm_dimension))
    allocate(adjointflux(fismatrix % fm_dimension))

    call fm_eigenvalue(0)
    forwardflux = fismatrix % eigenvector
    eigenvalue = fismatrix % fm_keff
    call fm_eigenvalue(1)
    adjointflux = fismatrix % eigenvector

    normfa = 0
    do k = 1, fismatrix % fm_dimension
      normfa = normfa + forwardflux(k) * adjointflux(k)
    end do
    gptflux = ZERO
    tmp1 = ZERO

    do si = 1, 20  ! source iteration
      do i = 1, fismatrix % fm_dimension
        tmp1(i) = 0
        do j = 1, fismatrix % fm_dimension
          if (fismatrix % source(i) /=0) then
            tmp1(i) = tmp1(i) + fismatrix % greenterm(j,i) * &
                 gptflux(j) / fismatrix % source(i) / eigenvalue
          end if
        end do
        if (fismatrix % source(i) /=0) then
          ! here the fismatrix % source is a cumulative value
          tmp1(i) = tmp1(i) + gptsource(i) * n_inactive / fismatrix % source(i)
        end if
      end do
      norm = 0
      do k = 1, fismatrix % fm_dimension
        norm = norm + tmp1(k) * forwardflux(k)
      end do
      gptflux = tmp1 - norm * adjointflux / normfa
    end do

    if(allocated(tmp1)) deallocate (tmp1)
    if(allocated(forwardflux)) deallocate (forwardflux)
    if(allocated(adjointflux)) deallocate (adjointflux)

  end subroutine fm_fixedsource

!===============================================================================
! SCORE_GCLUTCH_SCACOL calculates all term is GPT-CLUTCH method
!===============================================================================

  subroutine sensitivity_gclutch_scacol(p)

    type(Particle), intent(in) :: p
    integer :: i_sen         ! index in sensitivities array
    integer :: i_response    ! index in responses array
    integer :: i_tally       ! index in tallies array
    integer :: i_mesh        ! index in meshes array
    integer :: imp_mesh      ! index in importance meshes array
    real(8) :: gptimportance ! for each response, calculate one importance

    type(SensitivityObject), pointer :: s
    type(ResponseObject), pointer :: r
    type(TallyObject), pointer :: t
    type(RegularMesh), pointer :: m
    type(RegularMesh), pointer :: impm

    integer :: k
    integer :: j
    integer :: bin_mesh
    integer :: bin_impmesh
    integer :: bin_energy
    integer :: bin_nuclide_numer
    integer :: bin_score_numer
    integer :: bin_nuclide_denom
    integer :: bin_score_denom

    call score_tally(p)  ! get all tally values at one collision point

    SENSITIVITY_LOOP: do i_sen = 1, n_sens
      s => sensitivities(i_sen)
      imp_mesh = senmesh_dict % get_key(s % impmeshid)
      impm => sen_meshes(imp_mesh)
      i_response = response_dict % get_key(s % response)
      r => responses(i_response)

      ! ======================================================================
      ! calculate the intragenerational term
      i_tally = resptally_dict % get_key(r % numerid)
      gptimportance = resptalresult(i_tally) / s % respnumer
      i_tally = resptally_dict % get_key(r % denomid)
      gptimportance = gptimportance - resptalresult(i_tally) / s % respdenom
      ! calculate the intergenerational term
      call get_mesh_bin(impm, p % coord(1) % xyz, bin_impmesh)
      gptimportance = gptimportance + p % wgt * s % importance(bin_impmesh) * &
           material_xs % nu_fission / material_xs % total / keff
      s % clutchsen(:,:,:,:) = s % clutchsen(:,:,:,:) + &
           s % cumtally(:,:,:,:) * gptimportance

      ! ======================================================================
      ! Mesh logic
      bin_mesh = 0
      i_mesh = senmesh_dict % get_key(s % meshid)
      m => sen_meshes(i_mesh)
      call get_mesh_bin(m, p % coord(1) % xyz, bin_mesh)

      ! ======================================================================
      ! Energy logic
      bin_energy = 0
      bin_energy = binary_search(s % energystructure, &   ! get energy bin
           s % n_energy_bins + 1, p % E)

      ! ======================================================================
      ! decide whether there is a direct sensitivity term
      i_tally = resptally_dict % get_key(r % numerid)
      t => resp_tallies(i_tally)
      call tally_in_senlist(s,t,bin_nuclide_numer,bin_score_numer)

      if (bin_nuclide_numer /=0 .and. bin_score_numer /=0 .and. &
           bin_mesh /=0 .and. bin_energy /=0) then
        s % clutchsen(bin_nuclide_numer,bin_score_numer,bin_mesh, bin_energy) = &
             s % clutchsen(bin_nuclide_numer,bin_score_numer,bin_mesh, bin_energy) + &
             resptalresult(i_tally) / s % respnumer
      end if


      i_tally = resptally_dict % get_key(r % denomid)
      t => resp_tallies(i_tally)
      call tally_in_senlist(s,t,bin_nuclide_denom,bin_score_denom)

      if (bin_nuclide_denom /=0 .and. bin_score_denom /=0 .and. &
           bin_mesh /=0 .and. bin_energy /=0) then
        s % clutchsen(bin_nuclide_denom,bin_score_denom,bin_mesh, bin_energy) = &
             s % clutchsen(bin_nuclide_denom,bin_score_denom,bin_mesh, bin_energy) - &
             resptalresult(i_tally) / s % respdenom
      end if

    end do SENSITIVITY_LOOP

  end subroutine sensitivity_gclutch_scacol

!===============================================================================
! SENSITIVITY_GCLUTCH_FISSION calculates fission related sensitivities
!===============================================================================

  subroutine sensitivity_gclutch_fission(p, mt_number)

    ! progenitor is a global parameter

    type(Particle), intent(in) :: p
    integer, intent(in) :: mt_number
    integer :: i_nuclide_sen
    integer :: i_score_sen
    integer :: i_mesh
    integer :: i_tally
    integer :: i_response

    integer :: energy_bin      ! index for energy bin
    integer :: score_bin       ! index for score bin
    integer :: nuclide_bin     ! index for nuclide bin
    integer :: imp_mesh_bin    ! index for importance mesh bin

    integer :: i
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    type(SensitivityObject), pointer :: t
    type(ResponseObject), pointer :: r
    type(Material),    pointer :: mat
    type(RegularMesh), pointer :: m
    real(8) :: gptimportance

    call score_tally(p)  ! get all tally values at one collision point
    ! A loop over all sensitivities is necessary

    SENSITIVITY_LOOP: do i = 1, n_sens

      ! clear the information of last sensitivity
      score_bin = 0
      energy_bin = 0
      nuclide_bin = 0
      imp_mesh_bin = 0

      ! Get index of tally and pointer to tally
      t => sensitivities(i)

      ! ======================================================================
      ! Mesh logic to determine the importance function
      if (senmesh_dict % has_key(t % impmeshid)) then
        i_mesh = senmesh_dict % get_key(t % impmeshid)
        m => sen_meshes(i_mesh)
      else
        call fatal_error("Could not find mesh " // trim(to_str(t % impmeshid)) &
             // " specified on sensitivity " // trim(to_str(t % id)))
      end if
      call get_mesh_bin(m, p % coord(1) % xyz, imp_mesh_bin)
      if (imp_mesh_bin == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Energy logic
      if (mt_number == FISSION_CHI) then
        energy_bin = p % energy_born
      else
        energy_bin = p % energy_fission
      end if
      if (energy_bin == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Nuclide logic,  to check if the collision nuclide is in the
      ! sensitivity tally list
      NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins

        ! Get index of nuclide in nuclides array
        i_nuclide_sen = t % nuclide_bins(k)

        ! the collision nuclide is in the sensitivity tally list
        if (i_nuclide_sen == p % nuclide_born) then

          nuclide_bin = k

          exit NUCLIDE_BIN_LOOP

        end if

      end do NUCLIDE_BIN_LOOP
      if (nuclide_bin == 0) cycle SENSITIVITY_LOOP

      ! ======================================================================
      ! Score logic,  to check if the score is in the
      ! sensitivity tally list
      SCORE_BIN_LOOP: do j = 1, t % n_score_bins

        ! determine what type of score bin
        i_score_sen = t % score_bins(j)

        ! the reaction type is in the sensitivity tally list
        if (i_score_sen == mt_number) then

          score_bin = j

          exit SCORE_BIN_LOOP

        end if

      end do SCORE_BIN_LOOP

      if (score_bin == 0) cycle SENSITIVITY_LOOP

      ! calculate the intragenerational term
      i_response = response_dict % get_key(t % response)
      r => responses(i_response)
      i_tally = resptally_dict % get_key(r % numerid)
      gptimportance = resptalresult(i_tally) / t % respnumer
      i_tally = resptally_dict % get_key(r % denomid)
      gptimportance = gptimportance - resptalresult(i_tally) / t % respdenom

      ! calculate the intergenerational term
      gptimportance = gptimportance + p % wgt * t % importance(imp_mesh_bin) * &
           material_xs % nu_fission / material_xs % total / keff


      t % clutchsen(nuclide_bin,score_bin,p % mesh_born,energy_bin) = &
           t % clutchsen(nuclide_bin,score_bin,p % mesh_born,energy_bin) + &
           gptimportance

    end do SENSITIVITY_LOOP

  end subroutine sensitivity_gclutch_fission

!===============================================================================
! SENSITIVITY_CLUTCH calculates sensitivities in clutch
!===============================================================================

  subroutine sensitivity_gptclutch(p)

    type(Particle), intent(in) :: p

    if (clutch_first) call sensitivity_gclutch_scacol(p)

    if (clutch_second) then
      call sensitivity_gclutch_fission(p, SCORE_FISSION)
      call sensitivity_gclutch_fission(p, FISSION_NUBAR)
      call sensitivity_gclutch_fission(p, FISSION_CHI)
      call sensitivity_gclutch_fission(p, p % mtnum_born)
    end if

  end subroutine sensitivity_gptclutch

!===============================================================================
! COLLECT_CLUTCHCAL means collecting sensitivity tallies from different cores
! (focused on CLUTCH method)
!===============================================================================

#ifdef MPI
  subroutine collect_gptclutchcal()

    type(SensitivityObject), pointer :: t
    integer :: i         ! sensitivity loop index
    integer :: n         ! total size of count variable neutrontally
    real(8) :: dummy     ! temporary receive buffer for non-root reductions

    ! A loop over all sensitivities is necessary
    SENSITIVITY_LOOP: do i = 1, n_sens
      ! Get index of tally and pointer to tally
      t => sensitivities(i)
      n = t % n_nuclide_bins * t % n_score_bins * &
           t % n_mesh_bins * t % n_energy_bins
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, t % clutchsen, n, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      else
        ! Receive buffer not significant at other processors
        call MPI_REDUCE(t % clutchsen, dummy, n, MPI_REAL8, MPI_SUM, 0, &
             MPI_COMM_WORLD, mpi_err)
      end if
    end do SENSITIVITY_LOOP

  end subroutine collect_gptclutchcal
#endif

!===============================================================================
! SENSITIVITY_CALC_GCLUTCH calculates sensitivities at every second clutch batches
!===============================================================================

  subroutine sensitivity_calc_gclutch()

    type(SensitivityObject), pointer :: t
    integer :: i
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    real(8) :: value

#ifdef MPI
    call collect_gptclutchcal()
#endif

    if (master) then
      ! A loop over all sensitivities is necessary
      SENSITIVITY_LOOP: do i = 1, n_sens
        ! Get index of tally and pointer to tally
        t => sensitivities(i)
        t % n_realizations = t % n_realizations + 1

        do k = 1, t % n_nuclide_bins
          do l = 1, t % n_score_bins
            do m = 1, t % n_mesh_bins
              do n = 1, t % n_energy_bins
                value = t%clutchsen(k,l,m,n)
                t%results(1,k,l,m,n) = t%results(1,k,l,m,n) + value
                t%results(2,k,l,m,n) = t%results(2,k,l,m,n) + value *value
              end do
            end do
          end do
        end do

      end do SENSITIVITY_LOOP
    end if

  end subroutine sensitivity_calc_gclutch

!===============================================================================
! TALLY_IN_SENLIST decides whether one tally is in the sensitivity list
!===============================================================================

  subroutine tally_in_senlist(s,t,nuclide,score)

    type(SensitivityObject), intent(in) :: s
    type(TallyObject), intent(in) :: t
    integer, intent(inout) :: nuclide
    integer, intent(inout) :: score
    integer :: i
    integer :: i_nuclide_sen
    integer :: i_score_sen

    nuclide = 0
    score   = 0

    ! ======================================================================
    ! Nuclide logic,  to check if the response tally's nuclide is in the
    ! sensitivity tally list
    NUCLIDE_BIN_LOOP: do i = 1, s % n_nuclide_bins
      ! Get index of nuclide in nuclides array
      i_nuclide_sen = s % nuclide_bins(i)
      ! the collision nuclide is in the sensitivity tally list
      if (s % nuclide_bins(i)  == t % nuclide_bins(1)) then
        nuclide = i
        exit NUCLIDE_BIN_LOOP
      end if
    end do NUCLIDE_BIN_LOOP

    ! ======================================================================
    ! Score logic,  to check if the response tally's score is in the
    ! sensitivity tally list
    SCORE_BIN_LOOP: do i = 1, s % n_score_bins
      ! determine what type of score bin
      i_score_sen = s % score_bins(i)
      ! the reaction type is in the sensitivity tally list
      if (s % score_bins(i) == t % score_bins(1)) then
        score = i
        exit SCORE_BIN_LOOP
      end if
    end do SCORE_BIN_LOOP

  end subroutine tally_in_senlist

!===============================================================================
! SCORE_TALLY calculates different tallies, update the tally values, and store
! in resptalresult array, and this subroutine can be used in many occasions.
!===============================================================================

  subroutine score_tally(p)

    type(Particle), intent(in) :: p

    integer :: i
    integer :: i_tally
    integer :: i_filt
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array (from bins)
    integer :: matching_bins        ! matching bin for filter
    real(8) :: filter_weights       ! weight for filter
    real(8) :: flux                 ! collision estimate of flux
    real(8) :: atom_density         ! atom density of single nuclide
    !   in atom/b-cm
    real(8) :: filter_weight        ! combined weight of all filters
    type(TallyObject), pointer :: t
    type(Material),    pointer :: mat

    ! Determine collision estimate of flux
    !if (survival_biasing) then
    !  ! We need to account for the fact that some weight was already absorbed
    !  flux = (p % last_wgt + p % absorb_wgt) / material_xs % total
    !else
    !  flux = p % last_wgt / material_xs % total
    !end if
    flux = p % wgt / material_xs % total

    resptalresult = 0

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, n_resptallies
      ! Get index of tally and pointer to tally
      i_tally = i
      t => resp_tallies(i_tally)
      ! Find bin in each filter, only exist up to one bin in each filter
      do i_filt = 1, size(t % filters)
        call t % filters(i_filt) % obj % get_next_bin(p, ESTIMATOR_COLLISION, &
             NO_BIN_FOUND, respmatching_bins(i_filt), respfilter_weights(i_filt))
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (respmatching_bins(i_filt) == NO_BIN_FOUND) cycle TALLY_LOOP
      end do
      ! ======================================================================
      ! Nuclide logic

      if (t % all_nuclides) then
        if (p % material /= MATERIAL_VOID) then
          i_nuclide = -1
          atom_density = ZERO
        end if
      else
        i_nuclide = t % nuclide_bins(1)
        if (i_nuclide > 0) then
          if (p % material /= MATERIAL_VOID) then
            ! Get pointer to current material
            mat => materials(p % material)

            ! Determine if nuclide is actually in material
            NUCLIDE_MAT_LOOP: do j = 1, mat % n_nuclides
              ! If index of nuclide matches the j-th nuclide listed in the
              ! material, break out of the loop
              if (i_nuclide == mat % nuclide(j)) exit

              ! If we've reached the last nuclide in the material, it means
              ! the specified nuclide to be tallied is not in this material,
              ! and there isn't contribution of this collision to i_tally
              if (j == mat % n_nuclides) then
                flux = 0
              end if
            end do NUCLIDE_MAT_LOOP

            atom_density = mat % atom_density(j)
          else
            atom_density = ZERO
          end if
        end if

      end if

      ! Determine score for each bin
      call score_tally_general(p, t, 1, i_nuclide, &
           atom_density, flux, resptalresult(i_tally))

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

  end subroutine score_tally

!===============================================================================
! SCORE_TALLY_GENERAL tallies score for each response tally
!===============================================================================

  subroutine score_tally_general(p, t, filter_index, i_nuclide, &
       atom_density, flux, score)
    type(Particle),    intent(in)    :: p
    type(TallyObject), intent(inout) :: t
    integer,           intent(in)    :: i_nuclide
    integer,           intent(in)    :: filter_index   ! for % results
    real(8),           intent(in)    :: flux           ! flux estimate
    real(8),           intent(in)    :: atom_density   ! atom/b-cm
    real(8),           intent(out)   :: score          ! score of this collision

    integer :: l                    ! loop index for nuclides in material
    integer :: m                    ! loop index for reactions
    integer :: q                    ! loop index for scoring bins
    integer :: i_temp               ! temperature index
    integer :: i_nuc                ! index in nuclides array (from material)
    integer :: i_energy             ! index in nuclide energy grid
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    integer :: d                    ! delayed neutron index
    integer :: g                    ! delayed neutron index
    integer :: k                    ! loop index for bank sites
    integer :: d_bin                ! delayed group bin index
    integer :: dg_filter            ! index of delayed group filter
    real(8) :: yield                ! delayed neutron yield
    real(8) :: atom_density_        ! atom/b-cm
    real(8) :: f                    ! interpolation factor
    real(8) :: E                    ! particle energy

    ! determine what type of score bin
    score_bin = t % score_bins(1)

    !#########################################################################
    ! Determine appropirate scoring value.

    select case(score_bin)

    case (SCORE_FLUX)
      ! For flux, we need no cross section
      score = flux

    case (SCORE_TOTAL)
      if (i_nuclide > 0) then
        score = micro_xs(i_nuclide) % total * atom_density * flux
      else
        score = material_xs % total * flux
      end if

    case (SCORE_SCATTER)
      if (i_nuclide > 0) then
        score = (micro_xs(i_nuclide) % total &
             - micro_xs(i_nuclide) % absorption) * atom_density * flux
      else
        score = (material_xs % total - material_xs % absorption) * flux
      end if

    case (SCORE_ABSORPTION)
      if (i_nuclide > 0) then
        score = micro_xs(i_nuclide) % absorption * atom_density * flux
      else
        score = material_xs % absorption * flux
      end if

    case (SCORE_CAPTURE)
      if (i_nuclide > 0) then
        score = (micro_xs(i_nuclide) % absorption &
             - micro_xs(i_nuclide) % fission) * atom_density * flux
      else
        score = (material_xs % absorption - material_xs % fission) * flux
      end if

    case (SCORE_FISSION)
      if (i_nuclide > 0) then
        score = micro_xs(i_nuclide) % fission * atom_density * flux
      else
        score = material_xs % fission * flux
      end if

    case (SCORE_NU_FISSION)
      if (i_nuclide > 0) then
        score = micro_xs(i_nuclide) % nu_fission * atom_density * flux
      else
        score = material_xs % nu_fission * flux
      end if

    case (ELASTIC)
      if (i_nuclide > 0) then
        score = micro_xs(i_nuclide) % elastic * atom_density * flux
      else
        score = material_xs % elastic * flux
      end if

    case default
      ! Any other cross section has to be calculated on-the-fly. For
      ! cross sections that are used often (e.g. n2n, ngamma, etc. for
      ! depletion), it might make sense to optimize this section or
      ! pre-calculate cross sections
      if (score_bin > 1) then
        ! Set default score
        score = ZERO

        if (i_nuclide > 0) then
          if (nuclides(i_nuclide)%reaction_index%has_key(score_bin)) then
            m = nuclides(i_nuclide)%reaction_index%get_key(score_bin)

            ! Retrieve temperature and energy grid index and interpolation
            ! factor
            i_temp = micro_xs(i_nuclide) % index_temp
            i_energy = micro_xs(i_nuclide) % index_grid
            f = micro_xs(i_nuclide) % interp_factor

            associate (xs => nuclides(i_nuclide) % reactions(m) % xs(i_temp))
              if (i_energy >= xs % threshold) then
                score = ((ONE - f) * xs % value(i_energy - &
                     xs % threshold + 1) + f * xs % value(i_energy - &
                     xs % threshold + 2)) * atom_density * flux
              end if
            end associate
          end if
        end if

      else
        call fatal_error("Invalid score type on response tally " &
             // to_str(t % id) // ".")
      end if

    end select

  end subroutine score_tally_general

end module sensitivity
