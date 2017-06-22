module mesh_header

  use hdf5

  use hdf5_interface
  use string, only: to_str

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tessellation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type RegularMesh
    integer :: id                          ! user-specified id
    integer :: type                        ! rectangular, hexagonal
    integer :: n_dimension                 ! rank of mesh
    real(8) :: volume_frac                 ! volume fraction of each cell
    integer, allocatable :: dimension(:)   ! number of cells in each direction
    real(8), allocatable :: lower_left(:)  ! lower-left corner of mesh
    real(8), allocatable :: upper_right(:) ! upper-right corner of mesh
    real(8), allocatable :: width(:)       ! width of each mesh cell
  contains
    procedure :: to_hdf5 => regularmesh_to_hdf5
  end type RegularMesh

contains

  subroutine regularmesh_to_hdf5(this, group)
    class(RegularMesh), intent(in) :: this
    integer(HID_T), intent(in) :: group

    integer(HID_T) :: mesh_group

    mesh_group = create_group(group, "mesh " // trim(to_str(this % id)))

    call write_dataset(mesh_group, "type", "regular")
    call write_dataset(mesh_group, "dimension", this % dimension)
    call write_dataset(mesh_group, "lower_left", this % lower_left)
    call write_dataset(mesh_group, "upper_right", this % upper_right)
    call write_dataset(mesh_group, "width", this % width)

    call close_group(mesh_group)
  end subroutine regularmesh_to_hdf5


end module mesh_header
