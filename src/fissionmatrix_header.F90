module fissionmatrix_header

  use constants, only: NONE, ZERO

  use mesh_header,      only: RegularMesh

  implicit none

!===============================================================================
! FISSIONMATRIXOBJECT stores the numerator, denominator and mesh
! for user-specified fissionmatrix.
!===============================================================================

  type FissionmatrixObject

    ! Fission matrix mesh
    type(RegularMesh), pointer :: fm_mesh

    ! Fission matrix numerator (greenterm) and denominator (source)
    real(8), allocatable :: greenterm(:,:)
    real(8), allocatable :: source(:)

    ! The eigenvector of the fission matrix problem
    real(8), allocatable :: eigenvector(:)

    ! The size of the fission matrix
    integer :: fm_dimension

    real(8) :: fm_keff

  end type FissionmatrixObject

end module fissionmatrix_header
