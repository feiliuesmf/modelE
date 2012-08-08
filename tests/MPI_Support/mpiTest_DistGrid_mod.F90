! *** MPI Test Cases ***
! This module contains tests of the dist_grid_mod module. 
module mpiTest_DistGrid_mod
  use pFUnit, only: assertEqual, assertTrue, amRoot, mpiCommunicator, numProcesses, &
    processRank, TestInfo_type
  use dist_grid_mod, only: setCommunicator, dist_grid, init_grid, destroy_grid, &
    getDomainBounds, AxisIndex, getAxisIndex, sumxpe, globalmin, globalmax, &
    broadcast, transp

  implicit none
  private

  public :: setUp
  public :: tearDown
  public :: fixture

  ! This data structure is used internally to hold a copy of the distGrid 
  type fixture
    type (DIST_GRID) :: distGrid
  end type fixture

  ! compute axes indices of a distributed grid
  public :: test_getAxisIndexNoOverlap
  public :: nproc_test_getAxisIndexNoOverlap
  integer, parameter :: nproc_test_getAxisIndexNoOverlap(3) = [ 1, 3, 5 ]

  public :: test_getAxisIndexNoGaps
  public :: nproc_test_getAxisIndexNoGaps
  integer, parameter :: nproc_test_getAxisIndexNoGaps(3) = [ 1, 3, 5 ]

  ! sum an array over processors without reducing its rank
  public :: test_sumxpeMasterPresent
  public :: nproc_test_sumxpeMasterPresent
  integer, parameter :: nproc_test_sumxpeMasterPresent(3) = [ 1, 3, 5 ]
  public :: test_sumxpeIncrPresent
  public :: nproc_test_sumxpeIncrPresent
  integer, parameter :: nproc_test_sumxpeIncrPresent(3) = [ 1, 3, 5 ]

  ! sum an array over processors without reducing its rank (Real(8))
  public :: test_sumxpe
  public :: nproc_test_sumxpe
  integer, parameter :: nproc_test_sumxpe(3) = [ 1, 3, 5 ]

  ! sum an array over processors without reducing its rank (integer)
  public :: test_sumxpeInt
  public :: nproc_test_sumxpeInt
  integer, parameter :: nproc_test_sumxpeint(3) = [ 1, 3, 5 ]

  ! determine min, max value across pes
  public :: test_globalmax1d
  public :: nproc_test_globalmax1d
  integer, parameter :: nproc_test_globalmax1d(3) = [ 1, 3, 5 ]
  ! determine min value across pes
  public :: test_globalminmax
  public :: nproc_test_globalminmax
  integer, parameter :: nproc_test_globalminmax(3) = [ 1, 3, 5 ]

  ! broadcast data to all PEs. (real(8))
  public :: test_broadcast
  public :: nproc_test_broadcast
  integer, parameter :: nproc_test_broadcast(3) = [ 1, 3, 5 ]
  ! broadcast data to all PEs. (integer)
  public :: test_broadcastInt
  public :: nproc_test_broadcastInt
  integer, parameter :: nproc_test_broadcastInt(3) = [ 1, 3, 5 ]

  ! TODO : need transp tests

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

  end subroutine setUp

  ! ----------------------------------------------------------------------
  subroutine tearDown(this, info)
    ! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    call destroy_grid(this%distGrid)
  end subroutine tearDown

! ----------------------------------------------------------------------
  subroutine test_getAxisIndexNoOverlap(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    Type (AxisIndex), Pointer :: ai(:,:)
    integer :: p1, p2

    allocate(ai(0:numProcesses(info)-1,3))

    ! test getAxisIndex
    call getAxisIndex(this%distGrid, ai)

    ! NOTE : This will not work for 2D decomposition (CS)
    do p1 = 0, numProcesses(info) - 1
      do p2 = p1+1, numProcesses(info) - 1
        call assertTrue(ai(p2,2)%min > ai(p1,2)%max)
      end do
    end do

    deallocate(ai)

  end subroutine test_getAxisIndexNoOverlap

! ----------------------------------------------------------------------
  subroutine test_getAxisIndexNoGaps(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    Type (AxisIndex), Pointer :: ai(:,:)
    integer :: p, j
    logical :: found

    allocate(ai(0:numProcesses(info)-1,3))

    ! test getAxisIndex
    call getAxisIndex(this%distGrid, ai)

    ! NOTE : This will not work for 2D decomposition (CS)
    do j = 1, JM
      found = .false.
      do p = 0, numProcesses(info) - 1
        if (ai(p,2)%min <= j .and. j <= ai(p,2)%max) then
          found = .true.
          exit
        end if
      end do
      call assertTrue(found)
    end do

    deallocate(ai)

  end subroutine test_getAxisIndexNoGaps

! ----------------------------------------------------------------------
  subroutine test_sumxpe(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    real(8), allocatable, dimension(:) :: distArray
    real(8) :: arithSum 

    allocate (distArray(IM))
    distArray = processRank(info) 

    ! test sumxpe
    call sumxpe(distArray)

    arithSum = numProcesses(info) * (numProcesses(info) - 1) /2
    if (amRoot(info)) then
      call assertEqual(arithSum, distArray)
    endif

    deallocate(distArray)

  end subroutine test_sumxpe

! ----------------------------------------------------------------------
  subroutine test_sumxpeMasterPresent(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    real(8), allocatable, dimension(:) :: distArray
    real(8), allocatable, dimension(:) :: globalArray
    real(8) :: arithSum 

    allocate (distArray(IM))
    distArray = processRank(info) 
    if (amRoot(info)) then
      allocate(globalArray(IM))
    else
      allocate(globalArray(1))
    endif

    ! test sumxpe
    call sumxpe(distArray, globalArray)

    arithSum = numProcesses(info) * (numProcesses(info) - 1) /2
    if (amRoot(info)) then
      call assertEqual(arithSum, globalArray)
    endif

    deallocate(distArray, globalArray)

  end subroutine test_sumxpeMasterPresent

! ----------------------------------------------------------------------
  subroutine test_sumxpeIncrPresent(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    real(8), allocatable, dimension(:) :: distArray
    real(8) :: arithSum 

    allocate (distArray(IM))
    distArray = processRank(info) + 1

    ! test sumxpe
    call sumxpe(distArray, increment=.true.)

    arithSum = numProcesses(info) * (numProcesses(info) + 1) /2
    if (amRoot(info)) then
      call assertEqual(arithSum, distArray)
    endif

    deallocate(distArray)

  end subroutine test_sumxpeIncrPresent

! ----------------------------------------------------------------------
  subroutine test_sumxpeInt(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    integer, allocatable, dimension(:) :: distArray
    integer, allocatable, dimension(:) :: globalArray
    integer :: arithSum

    allocate (distArray(IM))
    distArray = processRank(info) 
    if (amRoot(info)) then
      allocate(globalArray(IM))
    else
      allocate(globalArray(1))
    endif

    ! test sumxpe
    call sumxpe(distArray, globalArray)

    arithSum = numProcesses(info) * (numProcesses(info) - 1) /2
    if (amRoot(info)) then
      call assertEqual(arithSum, globalArray)
    endif

    deallocate(distArray, globalArray)

  end subroutine test_sumxpeInt

! ----------------------------------------------------------------------
  subroutine test_globalminmax(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    real(8) :: expectedMin, expectedMax
    real(8) :: actualMin, actualMax, distValue

    distValue = processRank(info) + 1

    ! test globalmin
    call globalmin(this%distGrid, distValue, actualMin) 

    expectedMin = 1
    call assertEqual(expectedMin, actualMin)

    ! test globalmax
    call globalmax(this%distGrid, distValue, actualMax) 

    expectedMax = numProcesses(info)
    call assertEqual(expectedMax, actualMax)

  end subroutine test_globalminmax

! ----------------------------------------------------------------------
  subroutine test_globalmax1d(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    integer, allocatable, dimension(:) :: distArray
    integer, allocatable, dimension(:) :: max_array
    integer :: my_max

    allocate (distArray(IM))
    allocate (max_array(IM))

    distArray = processRank(info) + 1

    ! test globalmax
    call globalmax(this%distGrid, distArray, max_array) 

    my_max = numProcesses(info)
    call assertEqual(my_max, max_array)

    deallocate(distArray, max_array)

  end subroutine test_globalmax1d

! ----------------------------------------------------------------------
  subroutine test_broadcast(this, info)
! ----------------------------------------------------------------------
    type (TestInfo_type), intent(in) :: info
    type (fixture) :: this
    real(8), allocatable, dimension(:) :: array1d
    real(8), allocatable, dimension(:,:) :: array2d
    real(8) :: expected

    allocate (array1d(IM))
    allocate (array2d(IM, JM))

    expected = 1.0
    if (amRoot(info)) then
      array1d = expected
      array2d = expected
    else
      array1d = -999
      array2d = -999
    end if

    ! test mpi_bcast
    call broadcast(this%distGrid, array1d)

    call assertEqual(expected, array1d)

    ! test mpi_bcast
    call broadcast(this%distGrid, array2d)

    call assertEqual(expected, array2d)

    deallocate(array1d, array2d)

  end subroutine test_broadcast

! ----------------------------------------------------------------------
  subroutine test_broadcastInt(this, info)
! ----------------------------------------------------------------------
    type (TestInfo_type), intent(in) :: info
    type (fixture) :: this
    integer, allocatable, dimension(:) :: array1d
    integer, allocatable, dimension(:,:) :: array2d
    integer :: expected

    allocate (array1d(IM))
    allocate (array2d(IM, JM))

    expected = 1.0
    if (amRoot(info)) then
      array1d = expected
      array2d = expected
    else
      array1d = -999
      array2d = -999
    end if

    ! test mpi_bcast
    call broadcast(this%distGrid, array1d)

    call assertEqual(expected, array1d)

    ! test mpi_bcast
    call broadcast(this%distGrid, array2d)

    call assertEqual(expected, real(array2d))

    deallocate(array1d, array2d)

  end subroutine test_broadcastInt

end module mpiTest_DistGrid_mod
