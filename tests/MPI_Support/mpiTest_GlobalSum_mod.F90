! *** MPI Test Cases ***
! This module contains tests of the GlobalSum_mod module. 
module mpiTest_GlobalSum_mod
  use pFUnit, only: assertEqual, amRoot, mpiCommunicator, numProcesses, &
    processRank, TestInfo_type
  use dist_grid_mod, only: setCommunicator, dist_grid, init_grid, destroy_grid, &
    getDomainBounds
  use GlobalSum_mod, only: GlobalSum

  implicit none
  private

  public :: setUp
  public :: tearDown
  public :: fixture

! This data structure is used internally to hold a copy of the distributed grid 
! and its domain bounds
  type fixture
    type (DIST_GRID) :: distGrid
    integer :: ibeg, iend, jbeg, jend      ! local dimensions
    integer :: ibegh, iendh, jbegh, jendh  ! local dimensions + halo
  end type fixture

! output a bit-reproducible global-hemisphere-zonal sum for an array
  public :: test_GlobalSumJ
  public :: nproc_test_GlobalSumJ
  integer, parameter :: NPROC_TEST_GLOBALSUMJ(3) = [ 1, 3, 5 ]

  integer, parameter :: IM = 5 
  integer, parameter :: JM = 10
  integer, parameter :: LM = 2

contains

! ----------------------------------------------------------------------
  subroutine setUp(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info

    call setCommunicator(mpiCommunicator(info))
    call init_grid(this%distGrid, IM, JM, LM)

    CALL getDomainBounds( this%distgrid, &
      I_STRT=this%ibeg, &
      I_STOP=this%iend, &
      J_STRT=this%jbeg, &
      J_STOP=this%jend)

  end subroutine setUp

! ----------------------------------------------------------------------
  subroutine tearDown(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    call destroy_grid(this%distGrid)
  end subroutine tearDown

! ----------------------------------------------------------------------
  subroutine test_GlobalSumJ(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:) :: distArray
    real(8) :: gsum
    real(8) :: arithSum
    integer :: j

    allocate (distArray(this%jbeg:this%jend))
    do j = this%jbeg, this%jend
      distArray(j) = j
    end do

    call globalSum(this%distGrid, distArray, gsum)

    if (amRoot(info)) then
      arithSum = JM * ( JM + 1) / 2
      call assertEqual(arithSum, gsum)
    end if
    deallocate(distArray)

  end subroutine test_GlobalSumJ

end module mpiTest_GlobalSum_mod
