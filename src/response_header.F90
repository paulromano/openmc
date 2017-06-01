module response_header

  use constants,           only: NONE, N_FILTER_TYPES
  use tally_header,        only: TallyObject
  use mesh_header,         only: RegularMesh

  use, intrinsic :: ISO_C_BINDING

  implicit none

!===============================================================================
! SENSITIVITYOBJECT describes a user-specified sensitivity. The region of phase
! space to tally in is given by the TallyFilters and the results are stored in a
! SensitivityResult array.
!===============================================================================

  type ResponseObject
    ! Basic data

    integer :: id                                ! user-defined identifier
    integer :: numerid                           ! numerator tally identifier
    integer :: denomid                           ! denominator tally identifier

  end type ResponseObject

end module response_header
