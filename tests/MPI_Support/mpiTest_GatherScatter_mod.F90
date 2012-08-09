! *** MPI Test Cases ***
module mpiTest_GatherScatter_mod
  use pFUnit, only: assertEqual, amRoot, mpiCommunicator, numProcesses, &
    processRank, TestInfo_type
  use dist_grid_mod, only: setCommunicator, dist_grid, init_grid, destroy_grid, &
    getDomainBounds
  use GatherScatter_mod, only: pack_data, unpack_data
  implicit none
  private

  public :: setUp
  public :: tearDown
  public :: fixture

  type fixture
    type (DIST_GRID) :: distGrid
    integer :: ibeg, iend, jbeg, jend ! local dimensions
  end type fixture

  public :: test_GatherReal
  public :: nproc_test_GatherReal
  integer, parameter :: NPROC_TEST_GATHERREAL(3) = [ 1, 3, 5 ]

  public :: test_ScatterReal
  public :: nproc_test_ScatterReal
  integer, parameter :: NPROC_TEST_SCATTERREAL(3) = [ 1, 3, 5 ]

  public :: test_GatherInteger
  public :: nproc_test_GatherInteger
  integer, parameter :: NPROC_TEST_GATHERINTEGER(3) = [ 1, 3, 5 ]

  public :: test_ScatterInteger
  public :: nproc_test_ScatterInteger
  integer, parameter :: NPROC_TEST_SCATTERINTEGER(3) = [ 1, 3, 5 ]

  integer, parameter :: IM = 4 
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

    call getDomainBounds( this%distgrid, &
      I_STRT=this%ibeg, &
      I_STOP=this%iend, &
      J_STRT=this%jbeg, &
      J_STOP=this%jend  )

  end subroutine setUp

! ----------------------------------------------------------------------
  subroutine tearDown(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    call destroy_grid(this%distGrid)
  end subroutine tearDown

! ----------------------------------------------------------------------
  subroutine test_GatherReal(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:) :: distArray
    real(8), allocatable, dimension(:,:) :: globalArray
    integer :: j

    allocate (distArray(this%ibeg:this%iend, this%jbeg:this%jend))    
    allocate (globalArray(IM, JM)) 
    do j = this%jbeg, this%jend   
      distArray(:,j) = j
    end do

    ! test gather (pack_data)
    call pack_data(this%distgrid, distArray, globalArray)

    if ( amRoot(info) ) then
      do j = 1, JM
        call assertEqual(j, globalArray(:,j) )
      end do
    end if

    deallocate(distArray, globalArray)

  end subroutine test_GatherReal

! ----------------------------------------------------------------------
  subroutine test_GatherInteger(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    integer, allocatable, dimension(:,:) :: distArray
    integer, allocatable, dimension(:,:) :: globalArray
    integer :: j

    allocate (distArray(this%ibeg:this%iend, this%jbeg:this%jend))    
    allocate (globalArray(IM, JM)) 
    do j = this%jbeg, this%jend   
      distArray(:,j) = j
    end do

    ! test gather (pack_data)
    call pack_data(this%distgrid, distArray, globalArray)

    if ( amRoot(info) ) then
      do j = 1, JM
        call assertEqual(j, globalArray(:,j) )
      end do
    end if

    deallocate(distArray, globalArray)

  end subroutine test_GatherInteger

! ----------------------------------------------------------------------
  subroutine test_ScatterReal(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:) :: distArray
    real(8), allocatable, dimension(:,:) :: globalArray
    integer :: j

    allocate (distArray(this%ibeg:this%iend, this%jbeg:this%jend))    
    allocate(globalArray(IM, JM))

    if ( amRoot(info) ) then
      do j=1,JM
        globalArray(:,j) = j
      end do
    end if

    ! test scatter (unpack_data)
    call unpack_data(this%distGrid, globalArray, distArray)

    do j = this%jbeg, this%jend   
      call assertEqual(j, distArray(:,j) )
    end do

    deallocate(distArray, globalArray)

  end subroutine test_ScatterReal

! ----------------------------------------------------------------------
  subroutine test_ScatterInteger(this, info)
    ! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    integer, allocatable, dimension(:,:) :: distArray
    integer, allocatable, dimension(:,:) :: globalArray
    integer :: j

    allocate (distArray(this%ibeg:this%iend, this%jbeg:this%jend))    
    allocate (globalArray(IM, JM)) 

    if ( amRoot(info) ) then
      do j=1,JM
        globalArray(:,j) = j
      end do
    end if

    ! test scatter (unpack_data)
    call unpack_data(this%distGrid, globalArray, distArray)

    do j = this%jbeg, this%jend   
      call assertEqual(j, distArray(:,j) )
    end do

    deallocate(distArray, globalArray)

  end subroutine test_ScatterInteger

end module mpiTest_GatherScatter_mod
