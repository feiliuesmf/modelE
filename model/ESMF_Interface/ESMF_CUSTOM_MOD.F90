#include "assert.h"
Module ESMF_CUSTOM_MOD

#ifdef USE_ESMF
  Use ESMF_MOD, Only: ESMF_MAXSTR
  Use ESMF_MOD, Only: ESMF_DistGrid
  Use ESMF_MOD, Only: ESMF_DistGridGet
  Use ESMF_MOD, Only: ESMF_GridGet
  Use ESMF_MOD, Only: ESMF_DELayoutGet

  Use ESMF_MOD, Only: tGrid => ESMF_Grid
  Use ESMF_MOD, Only: Field => ESMF_Field
  Use ESMF_MOD, Only: ESMF_Grid
  Use ESMF_MOD, Only: ESMF_Field
  Use ESMF_MOD, Only: KIND_R4 => ESMF_KIND_R4
  Use ESMF_MOD, Only: KIND_R8 => ESMF_KIND_R8
  Use ESMF_MOD, Only: ESMF_SUCCESS

  Use ESMF_MOD_private, Only: Field_GetDataPointer
  Use ESMF_MOD_private, Only: Field_SetDataPointer
  Use ESMF_MOD, Only: VM => ESMF_VM
  Use ESMF_MOD, Only: ESMF_DELayout

  Implicit None
  Private

  Public :: ESMF_GridGetAxisIndex

  Public :: tGrid
  Public :: Field

  Public :: Field_GetDataPointer
  Public :: Field_SetDataPointer

  Public :: KIND_R4
  Public :: KIND_R8
  Public :: ESMF_SUCCESS

  ! To avoid the need to modify modelE interfaces for passing
  ! the ESMF grid throughout the source code, this one object is
  ! mode publicly available.
  Public :: modelE_grid
  Public :: modelE_vm
  Public :: modelE_default_layout

  Logical :: initQ = .false.
  Type (ESMF_Grid), Target, Save :: modelE_grid
  Type (VM), Target, Save :: modelE_vm
  Type (ESMF_DELayout), Target, Save :: modelE_default_layout

  type ESMF_AxisIndex
     sequence
     integer :: min
     integer :: max
     integer :: stride
  end type ESMF_AxisIndex
  Public :: ESMF_AxisIndex

  Integer, Parameter :: N_DIMENSIONS = 3

Contains

  subroutine ESMF_GridGetAxisIndex(grid, axisIndex, my_pet)

    type (ESMF_Grid), intent(IN) :: grid
    type (ESMF_AxisIndex) :: axisIndex(0:,:)
    integer, intent(in) :: my_pet

! local vars
    integer                               :: status
    character(len=ESMF_MAXSTR)            :: IAm='ESMF_GridGetAxisIndex'

    type (ESMF_DistGrid)                  :: distGrid
    type(ESMF_DELayout)                   :: LAYOUT
    integer,               allocatable    :: AL(:,:)
    integer,               allocatable    :: AU(:,:)
    integer                               :: nDEs
    integer                               :: gridRank
    integer :: p, npes_end, deId, I1,IN,J1,JN
    integer                               :: deList(1)

    call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=STATUS)
    call ESMF_DistGridGet(distGRID, delayout=layout, rc=STATUS)
    call ESMF_DELayoutGet(layout, deCount =nDEs, localDeList=deList, rc=status)
    deId = deList(1)
!    print *,'deId = ',deId

    allocate (AL(gridRank,0:nDEs-1),  stat=status)
    allocate (AU(gridRank,0:nDEs-1),  stat=status)

    call ESMF_DistGridGet(distgrid, &
         minIndexPDimPDe=AL, maxIndexPDimPDe=AU, rc=status)

    I1 = AL(1, deId)
    IN = AU(1, deId)
!    ASSERT_(gridRank > 1) !ALT: tilegrid is 1d (without RC this only for info)
    J1 = AL(2, deId)
    JN = AU(2, deId)
!    write(*,23)my_pet,' : I1,IN,J1,JN = ',I1,IN,J1,JN
!23 format(i3,a,4i3)

    do p = 0, nDEs - 1
      axisIndex(p,1)%min = AL(1, p) ! I1
      axisIndex(p,1)%max = AU(1, p) ! IN
      axisIndex(p,2)%min = AL(2, p) ! J1
      axisIndex(p,2)%max = AU(2, p) ! JN
    end do

    npes_end = size(axisIndex,1)
    axisIndex(nDEs:npes_end-1,2)%min = axisIndex(nDEs-1,2)%max + 1
    axisIndex(nDEs:npes_end-1,2)%max = axisIndex(nDEs-1,2)%max

    deallocate(AU, AL)

  end subroutine ESMF_GridGetAxisIndex

#endif

End Module ESMF_CUSTOM_MOD
