module sensitivity_header

  use constants,           only: NONE, N_FILTER_TYPES
  use tally_filter_header, only: TallyFilterContainer
  use mesh_header,         only: RegularMesh

  use, intrinsic :: ISO_C_BINDING

  implicit none

!===============================================================================
! SENSITIVITYOBJECT describes a user-specified sensitivity. The region of phase
! space to tally in is given by the TallyFilters and the results are stored in a
! SensitivityResult array.
!===============================================================================

  type SensitivityObject
    ! Basic data

    integer :: id                   ! user-defined identifier
    integer :: meshid               ! user-defined mesh identifier
    integer :: impmeshid            ! importance mesh identifier (for CLUTCH)
    integer :: blocklen = 0         ! user-defined block length (IFP and CLUTCH(IFP))
    character(len=104) :: name = "" ! user-defined name
    ! If method == 1, it is IFP method
    ! If method == 2, it is CLUTCH method using IFP to calculate IMP
    ! If method == 3, it is CLUTCH method using FM to calculate IMP, blocklen = 0

    ! response numerator and denominator for a GPT sensitivity calculation
    integer              :: response = 0  ! response id
    real(8)              :: respnumer = 0
    real(8)              :: respdenom = 0

    ! Individual nuclides to tally
    integer              :: n_nuclide_bins = 0
    integer, allocatable :: nuclide_bins(:)

    ! Values to tally, for different MT reaction for different nuclide
    integer              :: n_score_bins = 0
    integer, allocatable :: score_bins(:)

    ! mesh bins and energy bins
    integer              :: n_mesh_bins = 0
    integer              :: n_energy_bins = 0
    real(8), allocatable :: energystructure(:)

    ! importance mesh bins
    integer              :: imp_mesh_bins = 0
    real(8), allocatable :: importance(:)

    ! Results for each bin -- the first dimension of the array is for scores
    ! (e.g. different reaction for different nuclide) and the
    ! second dimension of the array is for the combination of filters
    ! (e.g. mesh discretization and specific energy group, etc.)

    ! sensitivity for different mesh, different energy, different nuclide
    ! and reaction
    real(8), allocatable :: results(:,:,:,:,:)

    ! cumulative tally for a single tracking
    real(8), allocatable :: cumtally(:,:,:,:)

    ! clutch sensitivity for a process
    real(8), allocatable :: clutchsen(:,:,:,:)

    ! clutch denon for a process
    real(8)              :: denom

    ! cumulative tally for different secondary neutrons
    real(8), allocatable :: secondtally(:,:,:,:,:)

    ! tally result of different progenitor
    real(8), allocatable :: neutrontally(:,:,:,:,:)

    ! neutron importance of different progenitor
    real(8), allocatable :: neutronvalue(:)

    ! gpt importance of different progenitor
    real(8), allocatable :: gptvaluenumer(:)
    real(8), allocatable :: gptvaluedenom(:)

    ! gpt importance tally used in fission matrix calculation
    real(8), allocatable :: tallynumer(:)
    real(8), allocatable :: tallydenom(:)

    ! used in the calculation of denominator of sensitivity
    real(8), allocatable :: neutronfission(:)

    ! reset property - allows a sen_tally to be reset after every batch
    logical :: reset = .false.

    ! Number of realizations of tally random variables
    integer :: n_realizations = 0

    ! Number of realizations of tally random variables
    integer :: imp_realizations = 0

  end type SensitivityObject

end module sensitivity_header
