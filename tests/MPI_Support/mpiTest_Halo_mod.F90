! *** MPI Test Cases ***
! This module contains tests of the Halo_mod module. 
module mpiTest_Halo_mod
  use pFUnit, only: assertEqual, assertTrue, amRoot, mpiCommunicator, numProcesses, &
    processRank, TestInfo_type
  use dist_grid_mod, only: setCommunicator, dist_grid, init_grid, &
    destroy_grid, getDomainBounds, NORTH, SOUTH
  use Halo_mod
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

  public :: test_HaloUpdateSouthPole1D
  public :: nproc_test_HaloUpdateSouthPole1D
  integer, parameter :: NPROC_TEST_HALOUPDATESOUTHPOLE1D(4) = [1, 2 , 9, 10]

  public :: test_HaloUpdateNorthPole1D
  public :: nproc_test_HaloUpdateNorthPole1D
  integer, parameter :: NPROC_TEST_HALOUPDATENORTHPOLE1D(4) = [1, 2 , 9, 10]

  public :: test_HaloUpdateSouthPole2D
  public :: nproc_test_HaloUpdateSouthPole2D
  integer, parameter :: NPROC_TEST_HALOUPDATESOUTHPOLE2D(4) = [1, 2, 9, 10]

  public :: test_HaloUpdateNorthPole2D
  public :: nproc_test_HaloUpdateNorthPole2D
  integer, parameter :: NPROC_TEST_HALOUPDATENORTHPOLE2D(4) = [1, 2 , 9, 10]

  public :: test_HaloUpdateSouthPole3D
  public :: nproc_test_HaloUpdateSouthPole3D
  integer, parameter :: NPROC_TEST_HALOUPDATESOUTHPOLE3D(4) = [1, 2 , 9, 10]

  public :: test_HaloUpdateNorthPole3D
  public :: nproc_test_HaloUpdateNorthPole3D
  integer, parameter :: NPROC_TEST_HALOUPDATENORTHPOLE3D(4) = [1, 2 , 9, 10]

  public :: test_HaloUpdate1d ! j
  public :: nproc_test_HaloUpdate1d
  integer, parameter :: NPROC_TEST_HALOUPDATE1D(3) = [ 1, 3, 5 ]

  public :: test_HaloUpdate2d ! i,j
  public :: nproc_test_HaloUpdate2d
  integer, parameter :: NPROC_TEST_HALOUPDATE2D(3) = [ 1, 3, 5 ]

  public :: test_HaloUpdate2dint ! i,j
  public :: nproc_test_HaloUpdate2dint
  integer, parameter :: NPROC_TEST_HALOUPDATE2DINT(3) = [ 1, 3, 5 ]

  public :: test_HaloUpdateJ2d ! j,LM
  public :: nproc_test_HaloUpdateJ2d
  integer, parameter :: NPROC_TEST_HALOUPDATEJ2D(3) = [ 1, 3, 5 ]

  public :: test_HaloUpdateColumn ! k,i,j
  public :: nproc_test_HaloUpdateColumn
  integer, parameter :: NPROC_TEST_HALOUPDATECOLUMN(3) = [ 1, 3, 5 ]

  public :: test_HaloUpdateColumnInt ! k,i,j
  public :: nproc_test_HaloUpdateColumnInt
  integer, parameter :: NPROC_TEST_HALOUPDATECOLUMNINT(3) = [ 1, 3, 5 ]

  public :: test_HaloUpdateBlock ! k,LM,i,j
  public :: nproc_test_HaloUpdateBlock
  integer, parameter :: NPROC_TEST_HALOUPDATEBLOCK(3) = [ 1, 3, 5 ]

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

    CALL getDomainBounds( this%distgrid, &
      I_STRT=this%ibeg, &
      I_STOP=this%iend, &
      J_STRT=this%jbeg, &
      J_STOP=this%jend, &
      I_STRT_HALO=this%ibegh, &
      I_STOP_HALO=this%iendh, &
      J_STRT_HALO=this%jbegh, &
      J_STOP_HALO=this%jendh  )

  end subroutine setUp

! ----------------------------------------------------------------------
  subroutine tearDown(this, info)
! ---------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type) :: info
    call destroy_grid(this%distGrid)
  end subroutine tearDown

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateSouthPole1D(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:) :: distArray
    integer :: rank

    rank = processRank(info)
    ! real array with dimensions (J)
    allocate (distArray(this%jbegh:this%jendh))     
    distArray = rank

    call halo_update_south_pole(this%distGrid, distArray)

    if (this%jend == 1) &
      call assertEqual(1, distArray(this%jendh))

    if (this%jbeg == 2) &
      call assertEqual(0, distArray(this%jbegh))

! Do not update halos except at SP
    if (this%jend > 1) &
       call assertEqual(rank, distArray(this%jendh))

    ! also test that interior does not change
    call assertEqual(rank, distArray(this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdateSouthPole1D

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateNorthPole1D(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:) :: distArray
    integer :: rank

    rank = processRank(info)
    ! real array with dimensions (J)
    allocate (distArray(this%jbegh:this%jendh))     
    distArray = rank

    call halo_update_north_pole(this%distGrid, distArray)

    if (this%jbeg == JM) then 
      call assertEqual(rank-1, distArray(this%jbegh))
    else if (this%jend == JM-1) then
      call assertEqual(rank+1, distArray(this%jendh))
    end if

! Do not update halos except at NP
    if (this%jbeg < JM) &
       call assertEqual(rank, distArray(this%jbegh))

    ! also test that interior does not change
    call assertEqual(rank, distArray(this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdateNorthPole1D

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateSouthPole2D(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:) :: distArray
    integer :: rank

    rank = processRank(info)
    ! real array with dimensions (I,J)
    allocate (distArray(this%ibegh:this%iendh, this%jbegh:this%jendh))    
    distArray = rank

    call halo_update_south_pole(this%distGrid, distArray)

    if (this%jend == 1) &
      call assertEqual(1, distArray(:,this%jendh))

    if (this%jbeg == 2) &
      call assertEqual(0, distArray(:,this%jbegh))

! Do not update halos except at SP
    if (this%jend > 1) &
       call assertEqual(rank, distArray(:,this%jendh))

    ! also test that interior does not change
    call assertEqual(rank, distArray(this%ibeg:this%iend, this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdateSouthPole2D

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateNorthPole2D(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:) :: distArray
    integer :: rank

    rank = processRank(info)
    ! real array with dimensions (I,J)
    allocate (distArray(this%ibegh:this%iendh, this%jbegh:this%jendh))    
    distArray = rank

    call halo_update_north_pole(this%distGrid, distArray)

    if (this%jbeg == JM) then 
      call assertEqual(rank-1, distArray(:,this%jbegh))
    else if (this%jend == JM-1) then
      call assertEqual(rank+1, distArray(:,this%jendh))
    end if

! Do not update halos except at NP
    if (this%jbeg < JM) &
       call assertEqual(rank, distArray(:,this%jbegh))

    ! also test that interior does not change
    call assertEqual(rank, distArray(this%ibeg:this%iend, this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdateNorthPole2D

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateSouthPole3D(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:,:) :: distArray
    integer :: rank

    rank = processRank(info)
    ! real array with dimensions (I,J,K)
    allocate (distArray(this%ibegh:this%iendh, this%jbegh:this%jendh, LM))    
    distArray = rank

    call halo_update_south_pole(this%distGrid, distArray)

    if (this%jend == 1) &
      call assertEqual(1, distArray(:,this%jendh,:))

    if (this%jbeg == 2) &
      call assertEqual(0, distArray(:,this%jbegh,:))

! Do not update halos except at SP
    if (this%jend > 1) &
       call assertEqual(rank, distArray(:,this%jendh,:))

    ! also test that interior does not change
    call assertEqual(rank, distArray(this%ibeg:this%iend, this%jbeg:this%jend,:))

    deallocate(distArray)

  end subroutine test_HaloUpdateSouthPole3D

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateNorthPole3D(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:,:) :: distArray
    integer :: rank

    rank = processRank(info)
    ! real array with dimensions (I,J,K)
    allocate (distArray(this%ibegh:this%iendh, this%jbegh:this%jendh, LM))    
    distArray = rank

    call halo_update_north_pole(this%distGrid, distArray)

    if (this%jbeg == JM) then 
      call assertEqual(rank-1, distArray(:,this%jbegh,:))
    else if (this%jend == JM-1) then
      call assertEqual(rank+1, distArray(:,this%jendh,:))
    end if

! Do not update halos except at NP
    if (this%jbeg < JM) &
       call assertEqual(rank, distArray(:,this%jbegh,:))

    ! also test that interior does not change
    call assertEqual(rank, distArray(this%ibeg:this%iend, this%jbeg:this%jend,:))

    deallocate(distArray)

  end subroutine test_HaloUpdateNorthPole3D

! ----------------------------------------------------------------------
  subroutine test_HaloUpdate1d(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:) :: distArray
    integer :: peSouth, peNorth, rank

    rank = processRank(info)
    ! real array with dimensions (J)
    allocate (distArray(this%jbegh:this%jendh))    
    distArray = rank

    call halo_update(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1 

    if (rank > 0) call assertEqual(peSouth, distArray(this%jbegh))
    if (rank < numProcesses(info)-1) call assertEqual(peNorth, distArray(this%jendh))
    ! also test that interior does not change
    call assertEqual(rank, distArray(this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdate1d

! ----------------------------------------------------------------------
  subroutine test_HaloUpdate2d(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:) :: distArray
    integer :: peSouth, peNorth, rank

    rank = processRank(info)
    ! real array with dimensions (J)
    allocate (distArray(this%ibegh:this%iendh, this%jbegh:this%jendh))    
    distArray = rank

    call halo_update(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1 

    if (rank > 0) call assertEqual(peSouth, distArray(:,this%jbegh))
    if (rank < numProcesses(info)-1) call assertEqual(peNorth, distArray(:,this%jendh))
    ! also test that interior does not change
    call assertEqual(rank, distArray(this%ibeg:this%iend, this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdate2d

! ----------------------------------------------------------------------
  subroutine test_HaloUpdate2dint(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    integer, allocatable, dimension(:,:) :: distArray
    integer :: peSouth, peNorth, rank

    rank = processRank(info)
    ! real array with dimensions (J)
    allocate (distArray(this%ibegh:this%iendh, this%jbegh:this%jendh))    
    distArray = rank

    call halo_update(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1 

    if (rank > 0) call assertEqual(peSouth, distArray(:,this%jbegh))
    if (rank < numProcesses(info)-1) call assertEqual(peNorth, distArray(:,this%jendh))
    ! also test that interior does not change
    call assertEqual(real(rank), real(distArray(this%ibeg:this%iend, this%jbeg:this%jend)))

    deallocate(distArray)

  end subroutine test_HaloUpdate2dint

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateJ2d(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:) :: distArray
    integer :: peSouth, peNorth, rank
   
    rank = processRank(info)
    ! array (column) with dimensions (J,L)
    allocate (distArray(this%jbegh:this%jendh,LM))    
    distArray = rank

    call halo_updatej(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1 

    if (rank > 0) call assertEqual(peSouth, distArray(this%jbegh, :))
    if (rank < numProcesses(info)-1) call assertEqual(peNorth, distArray(this%jendh,:))
    ! also test that interior does not change
    call assertEqual(real(rank), real(distArray(this%jbeg:this%jend, :)))

    deallocate(distArray)

  end subroutine test_HaloUpdateJ2d

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateColumn(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:,:) :: distArray
    integer :: peSouth, peNorth, rank
   
    rank = processRank(info)
    ! real array with dimensions (L,I,J)
    allocate (distArray(LM, this%ibegh:this%iendh,this%jbegh:this%jendh))    
    distArray = rank

    call halo_update_column(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1

    if (rank > 0) call assertEqual(peSouth, distArray(:,:,this%jbegh))
    if (rank < numProcesses(info)-1) call assertEqual(peNorth, distArray(:,:,this%jendh))
    ! also test that interior does not change
    call assertEqual(rank, distArray(:, this%ibeg:this%iend,this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdateColumn

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateColumnInt(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    integer, allocatable, dimension(:,:,:) :: distArray
    integer :: peSouth, peNorth, rank
   
    rank = processRank(info)
    ! integer array with dimensions (L,I,J)
    allocate (distArray(LM, this%ibegh:this%iendh,this%jbegh:this%jendh))    
    distArray = rank

    call halo_update_column(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1

    if (rank > 0) call assertEqual(real(peSouth), real(distArray(:,:,this%jbegh)))
    if (rank < numProcesses(info)-1) call assertEqual(real(peNorth), real(distArray(:,:,this%jendh)))
    ! also test that interior does not change
    call assertEqual(real(rank), real(distArray(:, this%ibeg:this%iend,this%jbeg:this%jend)))

    deallocate(distArray)

  end subroutine test_HaloUpdateColumnInt

! ----------------------------------------------------------------------
  subroutine test_HaloUpdateBlock(this, info)
! ----------------------------------------------------------------------
    type (fixture) :: this
    type (TestInfo_type), intent(in) :: info
    real(8), allocatable, dimension(:,:,:,:) :: distArray
    integer :: peSouth, peNorth, rank
   
    rank = processRank(info)
    ! array  with dimensions (K,L,I,J)
    allocate (distArray(LM, LM, this%ibegh:this%iendh,this%jbegh:this%jendh))    
    distArray = rank

    call halo_update_block(this%distGrid, distArray)

    peSouth = rank - 1 
    peNorth = rank + 1

    if (rank > 0) call assertEqual(peSouth, distArray(:,:,:,this%jbegh))
    if (rank < numProcesses(info)-1) call assertEqual(peNorth, distArray(:,:,:,this%jendh))
    ! also test that interior does not change
    call assertEqual(rank, distArray(:, :, this%ibeg:this%iend,this%jbeg:this%jend))

    deallocate(distArray)

  end subroutine test_HaloUpdateBlock

end module mpiTest_Halo_mod
