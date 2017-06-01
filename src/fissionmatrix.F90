module fissionmatrix

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
  use particle_header,  only: Particle

  implicit none

contains

!===============================================================================
! COUNT_SOURCE_FOR_FM determines the number of source particles in each FM mesh
! cell. The 'source' will be used as the denominotor of the fission matrix
!===============================================================================

  subroutine count_source_for_fm()

    logical :: sites_outside ! were there sites outside the fission matrix mesh?

    call count_bank_sites_bin(fismatrix % fm_mesh, source_bank, fismatrix % source, &
         sites_outside=sites_outside, size_bank=work)

    ! Check for sites outside of the mesh
    if (master .and. sites_outside) then
      call fatal_error("Source sites outside of the FM mesh!")
    end if

  end subroutine count_source_for_fm

!===============================================================================
! COUNT_GREENTERM_ANALOG determines the numerator of the fission matrix by
! the analog estimator
!===============================================================================

  subroutine count_greenterm_analog(p, current_mesh)

    type(Particle), intent(in) :: p
    integer, intent(in) :: current_mesh
    integer :: i                    ! index for current mesh

    call get_mesh_bin(fismatrix % fm_mesh, p % coord(1) % xyz, i)

    fismatrix % greenterm(current_mesh,p % mesh_born_fm) = &
         fismatrix % greenterm(current_mesh,p % mesh_born_fm) + keff ! FIJ (J->I)

  end subroutine count_greenterm_analog

!===============================================================================
! COUNT_GREENTERM_COLLISION determines the numerator of the fission matrix by
! the collision estimator
!===============================================================================

  subroutine count_greenterm_collision(p)

    type(Particle), intent(in) :: p
    real(8) :: flux                 ! collision estimate of flux
    integer :: i                    ! index for current mesh

    ! Determine collision estimate of flux
    if (survival_biasing) then
      ! We need to account for the fact that some weight was already absorbed
      flux = (p % last_wgt + p % absorb_wgt) / material_xs % total
    else
      flux = p % last_wgt / material_xs % total
    end if

    call get_mesh_bin(fismatrix % fm_mesh, p % coord(1) % xyz, i)

    fismatrix % greenterm(i,p % mesh_born_fm) = fismatrix % greenterm(i,p % mesh_born_fm) + &
         material_xs % nu_fission * flux  ! FIJ (J->I)


  end subroutine count_greenterm_collision

!===============================================================================
! COUNT_GREENTERM_COLLISION determines the numerator of the fission matrix by
! the collision estimator
!===============================================================================

  subroutine count_greenterm_tracklength(p,distance)

    type(Particle), intent(in) :: p
    real(8),        intent(in) :: distance

    real(8) :: flux                 ! collision estimate of flux
    integer :: next_bin_mesh        ! index for next mesh
    real(8) :: mesh_weight          ! weight for corresponding mesh

    flux = p % wgt * distance

    call get_next_bin_fm(fismatrix, p, &
         NO_BIN_FOUND, next_bin_mesh, mesh_weight)

    if (next_bin_mesh /= NO_BIN_FOUND) then

      fismatrix % greenterm(next_bin_mesh,p % mesh_born_fm) = &
           fismatrix % greenterm(next_bin_mesh,p % mesh_born_fm) + &
           material_xs % nu_fission * flux * mesh_weight  ! FIJ (J->I)

    end if

    do while (next_bin_mesh /= NO_BIN_FOUND)

      call get_next_bin_fm(fismatrix, p, &
           next_bin_mesh, next_bin_mesh, mesh_weight)

      if (next_bin_mesh /= NO_BIN_FOUND) then

        fismatrix % greenterm(next_bin_mesh,p % mesh_born_fm) = &
             fismatrix % greenterm(next_bin_mesh,p % mesh_born_fm) + &
             material_xs % nu_fission * flux * mesh_weight  ! FIJ (J->I)

      end if

    end do

  end subroutine count_greenterm_tracklength

!===============================================================================
! GET_NEXT_BIN_FM gives the index for the next valid filter bin and a weight that
! will be applied to the flux, this subroutine will be called when using
! tracklength estimator to tally the greenterm of fission matrix
!===============================================================================

  subroutine get_next_bin_fm(this, p, current_bin, next_bin, weight)
    type(FissionmatrixObject), intent(in) :: this
    type(Particle),    intent(in)  :: p
    integer, value,    intent(in)  :: current_bin
    integer,           intent(out) :: next_bin
    real(8),           intent(out) :: weight

    integer, parameter :: MAX_SEARCH_ITER = 100 ! Maximum number of times we can
    !  can loop while trying to find
    !  the first intersection.

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

    ! Get a pointer to the mesh.
    m => this % fm_mesh

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

  end subroutine get_next_bin_fm

!===============================================================================
! FREE_FM_MEMORY free the memory of fission matrix
!===============================================================================

  subroutine free_fm_memory()

    if (associated(fismatrix)) then
      if (associated(fismatrix % fm_mesh)) then
        if (allocated(fismatrix % fm_mesh % lower_left)) &
             deallocate(fismatrix % fm_mesh % lower_left)
        if (allocated(fismatrix % fm_mesh % upper_right)) &
             deallocate(fismatrix % fm_mesh % upper_right)
        if (allocated(fismatrix % fm_mesh % width)) &
             deallocate(fismatrix % fm_mesh % width)
        if (allocated(fismatrix % fm_mesh % dimension)) &
             deallocate(fismatrix % fm_mesh % dimension)
        deallocate(fismatrix % fm_mesh)
      end if
      if(allocated(fismatrix % greenterm)) deallocate (fismatrix % greenterm)
      if(allocated(fismatrix % source)) deallocate (fismatrix % source)
      if(allocated(fismatrix % eigenvector)) deallocate (fismatrix % eigenvector)
      deallocate(fismatrix)
    end if

  end subroutine free_fm_memory

!===============================================================================
! COLLECT_FM means collecting fission matrix values from different cores
!===============================================================================

  subroutine collect_fm()

    integer :: n1        ! total size of count variable
    integer :: n2        ! total size of count variable
    real(8) :: dummy    ! temporary receive buffer for non-root reductions

    n1 = fismatrix % fm_dimension * fismatrix % fm_dimension
    n2 = fismatrix % fm_dimension

#ifdef MPI
    ! collect values from all processors
    if (master) then
      call MPI_REDUCE(MPI_IN_PLACE, fismatrix % greenterm, n1, MPI_REAL8, MPI_SUM, 0, &
           MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(MPI_IN_PLACE, fismatrix % source, n2, MPI_REAL8, MPI_SUM, 0, &
           MPI_COMM_WORLD, mpi_err)
    else
      ! Receive buffer not significant at other processors
      call MPI_REDUCE(fismatrix % greenterm, dummy, n1, MPI_REAL8, MPI_SUM, 0, &
           MPI_COMM_WORLD, mpi_err)
      call MPI_REDUCE(fismatrix % source, dummy, n2, MPI_REAL8, MPI_SUM, 0, &
           MPI_COMM_WORLD, mpi_err)
    end if
#endif

  end subroutine collect_fm

!===============================================================================
! FM_EIGENVALUE means the power iteration of fission matrix eigenvalue calculation
!===============================================================================

  subroutine fm_eigenvalue(run_adjoint)

    integer, intent(in) :: run_adjoint  ! run adjoint mode or not
    real(8), allocatable :: tmp(:)
    real(8) :: norm         ! norm of the eigenvector bT*b
    real(8) :: crossnorm    ! norm of the bT*A*b
    integer :: i            ! loop index i
    integer :: j            ! loop index j
    integer :: k            ! loop index k
    integer :: pi           ! loop index for power iteration

    !open(1023,file='fismatrix.dat')

    allocate(tmp(fismatrix % fm_dimension))

    fismatrix % eigenvector = ONE
    norm = 0
    do k = 1, fismatrix % fm_dimension
      norm = norm + fismatrix % eigenvector(k) * fismatrix % eigenvector(k)
    end do
    norm = sqrt(norm)
    fismatrix % eigenvector = fismatrix % eigenvector / norm


    do pi = 1, 1000
      do i = 1, fismatrix % fm_dimension
        tmp(i) = 0
        do j = 1, fismatrix % fm_dimension
          ! forward mode calculation
          if (run_adjoint == 0 .and. fismatrix % source(j) /=0) then
            tmp(i) = tmp(i) + fismatrix % greenterm(i,j) * &
                 fismatrix % eigenvector(j) / fismatrix % source(j)
          end if
          ! adjoint mode calculation
          if (run_adjoint == 1 .and. fismatrix % source(i) /=0) then
            tmp(i) = tmp(i) + fismatrix % greenterm(j,i) * &
                 fismatrix % eigenvector(j) / fismatrix % source(i)
          end if
        end do
      end do

      norm = 0
      do k = 1, fismatrix % fm_dimension
        norm = norm + tmp(k) * tmp(k)
      end do
      norm = sqrt(norm)
      fismatrix % eigenvector = tmp / norm
    end do

    !do i = 1, fismatrix % fm_dimension
    !   write(1023,*) fismatrix % eigenvector(i)
    !end do

    norm = 0
    crossnorm = 0
    do i = 1, fismatrix % fm_dimension
      norm = norm + fismatrix % eigenvector(i) * fismatrix % eigenvector(i)
      do j = 1, fismatrix % fm_dimension
        ! forward mode calculation
        if (run_adjoint == 0 .and. fismatrix % source(j) /=0) then
          crossnorm = crossnorm + fismatrix % eigenvector(i) * &
               fismatrix % greenterm(i,j) * fismatrix % eigenvector(j) / &
               fismatrix % source(j)
        end if
        ! adjoint mode calculation
        if (run_adjoint == 1 .and. fismatrix % source(i) /=0) then
          crossnorm = crossnorm + fismatrix % eigenvector(i) * &
               fismatrix % greenterm(j,i) * fismatrix % eigenvector(j) / &
               fismatrix % source(i)
        end if
      end do
    end do

    fismatrix % fm_keff = crossnorm / norm

    if (allocated(tmp)) deallocate (tmp)

  end subroutine fm_eigenvalue

end module fissionmatrix
