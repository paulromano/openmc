module tally

  use constants
  use cross_section, only: get_macro_xs
  use error,         only: fatal_error
  use global
  use output,        only: message
  use search,        only: binary_search
  use string,        only: int_to_str
  use tally_header,  only: TallyScore, TallyMapItem, TallyMapElement

#ifdef MPI
  use mpi
#endif

  implicit none

  integer, allocatable :: position(:)

contains

!===============================================================================
! CALCULATE_KEFF
!===============================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer(8)              :: total_bank ! total number of source sites
    integer                 :: n          ! active cycle number
    real(8)                 :: kcoll      ! keff collision estimator         
    real(8), save           :: k1 = 0.    ! accumulated keff
    real(8), save           :: k2 = 0.    ! accumulated keff**2
    real(8)                 :: std        ! stdev of keff over active cycles
    character(MAX_LINE_LEN) :: msg        ! output/error message
#ifdef MPI
    integer :: ierr
#endif

    msg = "Calculate cycle keff..."
    call message(msg, 8)

    ! set k1 and k2 at beginning of run
    if (i_cycle == 1) then
       k1 = ZERO
       k2 = ZERO
    end if

#ifdef MPI
    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         & MPI_COMM_WORLD, ierr)
#else
    total_bank = n_bank
#endif

    ! Collect statistics and print output
    if (master) then
       kcoll = real(total_bank)/real(n_particles)*keff
       if (i_cycle > n_inactive) then
          n = i_cycle - n_inactive
          k1 = k1 + kcoll
          k2 = k2 + kcoll**2
          keff = k1/n
          std  = sqrt((k2/n-keff**2)/n)
          if (i_cycle > n_inactive+1) then
             write(6,101) i_cycle, kcoll, keff, std
          else
             write(6,100) i_cycle, kcoll
          end if
       else
          write(6,100) i_cycle, kcoll
          keff = kcoll
       end if
    end if

#ifdef MPI
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
#endif

100 format (2X,I4,2X,F8.5)
101 format (2X,I4,2X,F8.5,9X,F8.5,1X,F8.5)

  end subroutine calculate_keff

!===============================================================================
! CREATE_TALLY_MAP creates a map that allows a quick determination of which
! tallies and bins need to be scored to when a particle makes a collision.
!===============================================================================

  subroutine create_tally_map()

    integer :: i
    integer :: j
    integer :: index
    integer :: n
    integer :: filter_bins
    integer :: score_bins
    character(MAX_LINE_LEN) :: msg        ! output/error message
    type(TallyObject), pointer :: t => null()

    ! allocate tally map array -- note that we don't need a tally map for the
    ! energy_in and energy_out filters
    allocate(tally_maps(TALLY_TYPES - 2))

    ! allocate list of items for each different filter type
    allocate(tally_maps(T_CELL)     % items(n_cells))
    allocate(tally_maps(T_SURFACE)  % items(n_surfaces))
    allocate(tally_maps(T_UNIVERSE) % items(n_universes))
    allocate(tally_maps(T_MATERIAL) % items(n_materials))
    allocate(tally_maps(T_MESH)     % items(100)) ! TODO: Change this
    allocate(tally_maps(T_CELLBORN)   % items(n_cells))

    ! Allocate and initialize tally map positioning for finding bins
    allocate(position(TALLY_TYPES))
    position = 0

    do i = 1, n_tallies
       t => tallies(i)

       ! initialize number of scoring bins
       filter_bins = 1

       ! determine if there are subdivisions for incoming or outgoing energy to
       ! adjust the number of filter bins appropriately
       n = t % n_bins(T_ENERGYOUT)
       t % stride(T_ENERGYOUT) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       n = t % n_bins(T_ENERGYIN)
       t % stride(T_ENERGYIN) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       ! TODO: Determine size of mesh to increase number of scoring bins

       ! Add map elements for surface bins
       n = t % n_bins(T_SURFACE)
       t % stride(T_SURFACE) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % surface_bins(j) % scalar
             call add_map_element(tally_maps(T_SURFACE) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cellborn bins
       n = t % n_bins(T_CELLBORN)
       t % stride(T_CELLBORN) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % cellborn_bins(j) % scalar
             call add_map_element(tally_maps(T_CELLBORN) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cell bins
       n = t % n_bins(T_CELL)
       t % stride(T_CELL) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % cell_bins(j) % scalar
             call add_map_element(tally_maps(T_CELL) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for material bins
       n = t % n_bins(T_MATERIAL)
       t % stride(T_MATERIAL) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % material_bins(j) % scalar
             call add_map_element(tally_maps(T_MATERIAL) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for universe bins
       n = t % n_bins(T_UNIVERSE)
       t % stride(T_UNIVERSE) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % universe_bins(j) % scalar
             call add_map_element(tally_maps(T_UNIVERSE) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Finally add scoring bins for the macro tallies
       n = t % n_macro_bins
       if (n > 0) then
          score_bins = n
       else
          msg = "Must have macro tally bins!"
          call fatal_error(msg)
       end if

       ! Allocate scores for tally
       t % n_total_bins = filter_bins
       allocate(t % scores(filter_bins, score_bins))

    end do

  end subroutine create_tally_map

!===============================================================================
! ADD_MAP_ELEMENT adds a pair of tally and bin indices to the list for a given
! cell/surface/etc.
!===============================================================================

  subroutine add_map_element(item, index_tally, index_bin)

    type(TallyMapItem), intent(inout) :: item
    integer, intent(in) :: index_tally ! index in tallies array
    integer, intent(in) :: index_bin   ! index in bins array

    integer :: n
    type(TallyMapElement), allocatable :: temp(:)

    if (.not. allocated(item % elements)) then
       allocate(item % elements(1))
       item % elements(1) % index_tally = index_tally
       item % elements(1) % index_bin   = index_bin
    else
       ! determine size of elements array
       n = size(item % elements)

       ! allocate temporary storage and copy elements
       allocate(temp(n+1))
       temp(1:n) = item % elements

       ! move allocation back to main array
       call move_alloc(FROM=temp, TO=item%elements)

       ! set new element
       item % elements(n+1) % index_tally = index_tally
       item % elements(n+1) % index_bin   = index_bin
    end if

  end subroutine add_map_element

!===============================================================================
! SCORE_TALLY contains the main logic for scoring user-specified tallies
!===============================================================================

  subroutine score_tally(p)

    type(Particle), pointer :: p     ! particle

    integer :: i
    integer :: j
    integer :: n
    integer :: bins(TALLY_TYPES)

    integer :: score_index
    real(8) :: score
    type(TallyObject), pointer :: t

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_tallies
       t => tallies(i)
       
       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       ! determine next cell bin
       if (t % n_bins(T_CELL) > 0) then
          bins(T_CELL) = get_next_bin(T_CELL, p % cell, i)
          if (bins(T_CELL) == NO_BIN_FOUND) cycle
       else
          bins(T_CELL) = 1
       end if

       ! determine next surface bin
       if (t % n_bins(T_SURFACE) > 0) then
          bins(T_SURFACE) = get_next_bin(T_SURFACE, p % surface, i)
          if (bins(T_SURFACE) == NO_BIN_FOUND) cycle
       else
          bins(T_SURFACE) = 1
       end if

       ! determine next universe bin
       if (t % n_bins(T_UNIVERSE) > 0) then
          bins(T_UNIVERSE) = get_next_bin(T_UNIVERSE, p % universe, i)
          if (bins(T_UNIVERSE) == NO_BIN_FOUND) cycle
       else
          bins(T_UNIVERSE) = 1
       end if

       ! determine next material bin
       if (t % n_bins(T_MATERIAL) > 0) then
          bins(T_MATERIAL) = get_next_bin(T_MATERIAL, p % material, i)
          if (bins(T_MATERIAL) == NO_BIN_FOUND) cycle
       else
          bins(T_MATERIAL) = 1
       end if

       ! determine next cellborn bin
       if (t % n_bins(T_CELLBORN) > 0) then
          bins(T_CELLBORN) = get_next_bin(T_CELLBORN, p % cell_born, i)
          if (bins(T_CELLBORN) == NO_BIN_FOUND) cycle
       else
          bins(T_CELLBORN) = 1
       end if

       ! determine incoming energy bin
       n = t % n_bins(T_ENERGYIN)
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % E < t % energy_in(1) .or. p % E > t % energy_in(n)) cycle

          ! search to find incoming energy bin
          bins(T_ENERGYIN) = binary_search(t % energy_in, n, p % E)
       else
          bins(T_ENERGYIN) = 1
       end if

       ! determine outgoing energy bin
       n = t % n_bins(T_ENERGYOUT)
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % E < t % energy_out(1) .or. p % E > t % energy_out(n)) cycle

          ! search to find incoming energy bin
          bins(T_ENERGYOUT) = binary_search(t % energy_out, n, p % E)
       else
          bins(T_ENERGYOUT) = 1
       end if

       ! =======================================================================
       ! DETERMINE INDEX IN SCORES ARRAY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       score_index = sum((bins - 1) * t % stride) + 1

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       do j = 1, t % n_macro_bins
          ! Determine score
          select case(t % macro_bins(j) % scalar)
          case (MACRO_FLUX)
             score = p % wgt / material_xs % total
          case (MACRO_TOTAL)
             score = p % wgt
          case (MACRO_SCATTER)
             score = p % wgt * (material_xs % total - material_xs % absorption) &
                  / material_xs % total
          case (MACRO_ABSORPTION)
             score = p % wgt * material_xs % absorption / material_xs % total
          case (MACRO_FISSION)
             score = p % wgt * material_xs % fission / material_xs % total
          case (MACRO_NU_FISSION)
             score = p % wgt * material_xs % nu_fission / material_xs % total
          end select

          ! Add score to tally
          call add_to_score(t % scores(score_index, j), score)

       end do

    end do

    ! Reset tally map positioning
    position = 0

  end subroutine score_tally

!===============================================================================
! GET_NEXT_BIN determines the next scoring bin for a particular filter variable
!===============================================================================

  function get_next_bin(i_map, i_cell, i_tally) result(bin)

    integer, intent(in) :: i_map
    integer, intent(in) :: i_cell
    integer, intent(in) :: i_tally
    integer             :: bin

    integer :: index_tally
    integer :: index_bin
    integer :: n

    n = size(tally_maps(i_map) % items(i_cell) % elements)

    do
       ! Increment position in elements
       position(i_map) = position(i_map) + 1

       ! If we've reached the end of the array, there is no more bin to score to
       if (position(i_map) > n) then
          position(i_map) = 0
          bin = NO_BIN_FOUND
          return
       end if

       index_tally = tally_maps(i_map) % items(i_cell) % &
            elements(position(i_map)) % index_tally
       index_bin = tally_maps(i_map) % items(i_cell) % &
            elements(position(i_map)) % index_bin

       if (index_tally > i_tally) then
          ! Since the index being checked against is greater than the index we
          ! need (and the tally indices were added to elements sequentially), we
          ! know that no more bins will be scoring bins for this tally
          position(i_map) = 0
          bin = NO_BIN_FOUND
          return
       elseif (index_tally == i_tally) then
          ! Found a match
          bin = index_bin
          return
       end if

    end do

  end function get_next_bin

!===============================================================================
! ADD_TO_SCORE accumulates a scoring contribution to a specific tally bin and
! specific response function. Note that we don't need to add the square of the
! contribution since that is done at the cycle level, not the history level
!===============================================================================

  subroutine add_to_score(score, val)

    type(TallyScore), intent(inout) :: score
    real(8),          intent(in)    :: val
    
    score % n_events    = score % n_events    + 1
    score % val_history = score % val_history + val
    
  end subroutine add_to_score

!===============================================================================
! SYNCHRONIZE_TALLIES accumulates the sum of the contributions from each history
! within the cycle to a new random variable
!===============================================================================

  subroutine synchronize_tallies()

    integer :: i
    integer :: j
    integer :: k
    real(8) :: val
    type(TallyObject), pointer :: t


    do i = 1, n_tallies
       t => tallies(i)

       do j = 1, t % n_total_bins
          do k = 1, t % n_macro_bins
             ! Add the sum and square of the sum of contributions from each
             ! history within a cycle to the variables val and val_sq. This will
             ! later allow us to calculate a variance on the tallies

             val = t % scores(j,k) % val_history / n_particles
             t % scores(j,k) % val    = t % scores(j,k) % val    + val
             t % scores(j,k) % val_sq = t % scores(j,k) % val_sq + val*val

             ! Reset the within-cycle accumulation variable

             t % scores(j,k) % val_history = ZERO
          end do
       end do

    end do

  end subroutine synchronize_tallies

end module tally
