#include "rundeck_opts.h"
#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#ifndef CUBED_SPHERE
#define DOMAIN_DECOMP_ATM_IS_1D
#endif

#ifdef NEW_IO
#define USE_DD2D_UTILS
#endif

MODULE dist_grid_mod
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@+    for the message passing (ESMF) implementation.
!@auth NCCS ASTG

#if ( defined USE_ESMF )  || ( defined USE_MPP )
#define USE_MPI
#endif

#ifdef USE_ESMF
   use ESMF_Mod
   use MpiSupport_mod, only: createDecompMpiType
   use ESMF_CUSTOM_mod, only: ESMF_AxisIndex, ESMF_GridGetAxisIndex
#endif
   use MpiSupport_mod, only: am_i_root
   use Domain_mod
   use Hidden_mod


! retaining for now, but disabling, the MPP+FVCUBED coding in this file
#undef USE_MPP

#ifdef CUBED_SPHERE
#define USE_DD2D_UTILS
#undef CUBED_SPHERE
#endif

#ifdef USE_DD2D_UTILS
   use dd2d_utils, only : dist_grid,init_dist_grid
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif
   SAVE
   PRIVATE ! Except for

!aoo since DIST_GRID is public ESMF_GRID has to be public
!aoo (SGI compiler complains)
   PUBLIC :: ESMF_GRID

#ifdef USE_ESMF
   public :: load_cap_config
   TYPE(ESMF_GridComp)  :: compmodelE
#endif


!@var DIST_GRID derived type to provide ESMF decomposition info
!@+   public components are used to minimize overhead for accessing
!@+   routine components
   PUBLIC :: DIST_GRID
   public :: setMpiCommunicator
!@var  grid Default decomposition; globally accessible for convenience.
   PUBLIC :: grid
!@var INIT_APP set some parameters and initialize ESMF
   PUBLIC :: INIT_APP
   PUBLIC :: INIT_GRID
   PUBLIC :: DESTROY_GRID
!@var FINISH_APP Cleans up at the end of the run (closes debugging file)
   PUBLIC :: FINISH_APP
!@var GLOBALMIN determine max value across pes
   PUBLIC :: GLOBALMIN
!@var GLOBALMAX determine max value across pes
   PUBLIC :: GLOBALMAX
!@var SUMXPE sum an array over processors without reducing its rank
   PUBLIC :: SUMXPE
!@var GET - extracts bounds information from DIST_GRID object
   PUBLIC :: GET
   PUBLIC :: HERE
   PUBLIC :: LOG_PARALLEL
  
   PUBLIC :: TRANSP
   PUBLIC :: TRANSPOSE_COLUMN
!@var GLOBALMIN Generic wrapper for Real
   INTERFACE GLOBALMIN
      MODULE PROCEDURE GLOBALMIN_R
   END INTERFACE

!@var GLOBALMAX Generic wrapper for Real/integer
   INTERFACE GLOBALMAX
      MODULE PROCEDURE GLOBALMAX_R
      MODULE PROCEDURE GLOBALMAX_I
      MODULE PROCEDURE GLOBALMAX_I_1D
   END INTERFACE

   INTERFACE TRANSP
      MODULE PROCEDURE TRANSPOSE_ij
      MODULE PROCEDURE TRANSPOSE_ijk
   END INTERFACE
  
   INTERFACE SUMXPE
      MODULE PROCEDURE SUMXPE_1D
      MODULE PROCEDURE SUMXPE_1D_I
      MODULE PROCEDURE SUMXPE_2D
      MODULE PROCEDURE SUMXPE_3D
      MODULE PROCEDURE SUMXPE_4D
!        MODULE PROCEDURE SUMXPE_5D
   END INTERFACE

#ifndef USE_ESMF
      ! Place holders for the real things
   TYPE ESMF_DELayout
      Integer :: i
   END TYPE ESMF_DELayout
   TYPE ESMF_GridComp
      Integer :: i
   END TYPE ESMF_GridComp
   TYPE ESMF_GRID
      Integer :: i
   END TYPE ESMF_GRID
   INTEGER, PARAMETER :: ESMF_MAXSTR = 40
   INTEGER, PARAMETER :: ESMF_KIND_R8 = Selected_Real_Kind(15)
#endif

   type (ESMF_DELayout) :: ESMF_LAYOUT
   public :: ESMF_LAYOUT
   public :: amRootESMF


!@var ESMF_BCAST Generic routine to broadcast data to all PEs.
   PUBLIC :: ESMF_BCAST
   INTERFACE ESMF_BCAST
      MODULE PROCEDURE ESMF_BCAST_0D
      MODULE PROCEDURE ESMF_BCAST_1D
      MODULE PROCEDURE ESMF_BCAST_2D
      MODULE PROCEDURE ESMF_BCAST_3D
      MODULE PROCEDURE ESMF_BCAST_4D
      MODULE PROCEDURE ESMF_IBCAST_0D
      MODULE PROCEDURE ESMF_IBCAST_1D
      MODULE PROCEDURE ESMF_IBCAST_2D
      MODULE PROCEDURE ESMF_IBCAST_3D
      MODULE PROCEDURE ESMF_IBCAST_4D
   END INTERFACE

!@var BAND_PACK Procedure in which each PE receives data from other PEs
!@+             to fill a contiguous pre-requested range of J indices
   public :: band_pack,band_pack_column
   interface band_pack
      module procedure band_pack_ij
      module procedure band_pack_ijl
   end interface
!@var BAND_PACK_TYPE a data structure needed by BAND_PACK, initialized
!@var via INIT_BAND_PACK_TYPE
   type band_pack_type
!        integer :: im_world
      integer :: j_strt,j_stop
      integer :: j_strt_halo,j_stop_halo
      integer :: jband_strt,jband_stop
      integer, dimension(:), pointer :: scnts,sdspl,sdspl_inplace
      integer, dimension(:), pointer :: rcnts,rdspl,rdspl_inplace
      integer, dimension(:), pointer :: j0_send,j1_send
      integer, dimension(:), pointer :: j0_recv,j1_recv
      integer :: mpi_comm
      integer :: npes_comm
   end type band_pack_type
!@var INIT_BAND_PACK_TYPE initialization routine during which each PE
!@+   requests a range of J indices and sets up the necessary send/receive
!@+   information for the BAND_PACK procedure
   public :: band_pack_type,init_band_pack_type

   PUBLIC SEND_TO_J
   interface SEND_TO_J
      module procedure SEND_TO_J_1D
      module procedure ISEND_TO_J_0D
   end interface

   PUBLIC RECV_FROM_J
   interface RECV_FROM_J
      module procedure RECV_FROM_J_1D
      module procedure IRECV_FROM_J_0D
   end interface

   INTEGER, PARAMETER :: HALO_WIDTH = 1
   integer ::  root


#ifndef USE_DD2D_UTILS
      ! Local grid information
   TYPE DIST_GRID

      type (Hidden_type) :: private
!!$$        type (Domain_type) :: localSubdomain

      TYPE (ESMF_Grid) :: ESMF_GRID
      ! Parameters for Global domain
      INTEGER :: IM_WORLD        ! Number of Longitudes
      INTEGER :: JM_WORLD        ! Number of latitudes
      ! Parameters for local domain
      LOGICAL :: HAVE_DOMAIN = .false.     ! Whether this PE has any of the domain
      INTEGER :: I_STRT          ! Begin local domain longitude index
      INTEGER :: I_STOP          ! End   local domain longitude index
      INTEGER :: J_STRT          ! Begin local domain latitude  index
      INTEGER :: J_STOP          ! End   local domain latitude  index
      INTEGER :: J_STRT_SKP      ! Begin local domain exclusive of S pole
      INTEGER :: J_STOP_SKP      ! End   local domain exclusive of N pole
      INTEGER :: ni_loc ! for transpose
      ! Parameters for halo of local domain
      INTEGER :: I_STRT_HALO     ! Begin halo longitude index
      INTEGER :: I_STOP_HALO     ! End   halo longitude index
      INTEGER :: J_STRT_HALO     ! Begin halo latitude  index
      INTEGER :: J_STOP_HALO     ! End   halo latitude  index
      ! Parameters for staggered "B" grid
      ! Note that global staggered grid begins at "2".
      INTEGER :: J_STRT_STGR     ! Begin local staggered domain
      INTEGER :: J_STOP_STGR     ! End   local staggered domain
      
      INTEGER, DIMENSION(:), POINTER :: DJ_MAP
      INTEGER :: DJ
      
   END TYPE DIST_GRID
#endif
   type (DIST_GRID) :: GRID
!not used      TYPE (DIST_GRID) :: GRID_TRANS

   public :: haveLatitude

! Remaining variables are private to the module.

!@var NPES_WORLD number of total processes
   INTEGER :: NPES_WORLD
!@var NP_LON number of azimuthal processes.
   INTEGER :: NP_LON
!@var NP_LAT number of meridional     processes.
   INTEGER :: NP_LAT
!@var MY_PET index of _this_ PET (analagous to MPI rank)
   INTEGER :: my_pet
!@var RANK_LON index of _this_ process in azimuthal set.
   INTEGER :: RANK_LON
!@var RANK_LAT_RANK index of _this_ process in meridional set.
   INTEGER :: RANK_LAT

   INTEGER, PUBLIC :: CHECKSUM_UNIT

   ! temporary public during refactoring
   public :: isPeriodic
   public :: my_pet, npes_world
#ifdef USE_ESMF
   public :: esmf_grid_my_pe_loc
   public :: esmf_grid_pe_layout
#endif
   public :: am_i_root

   ! Direction bits
   public :: NORTH, SOUTH
   Integer, Parameter :: NORTH = 2**0
   Integer, Parameter :: SOUTH = 2**1

   public :: getMpiCommunicator
   public :: getNumProcesses
   public :: getNumAllProcesses
   public :: getMpiTag
   public :: incrementMpiTag
   public :: hasSouthPole
   public :: hasNorthPole
   public :: hasPeriodicBC
   
   public :: ESMF_DELayout
   public :: ESMF_MAXSTR

   public :: getLogUnit

 CONTAINS

! ----------------------------------------------------------------------
   ! This routine initializes the quantities described above.
   ! The initialization should proceed prior to any grid computations.
   subroutine init_app
! ----------------------------------------------------------------------
#ifdef USE_ESMF
     USE ESMF_CUSTOM_MOD, Only: modelE_vm
#endif
     use FILEMANAGER, only: openunit
     integer :: rc
     character(len=10) :: fileName
     character(len=10) :: hostName
     character(len=80) :: command
!$    integer :: omp_get_max_threads, omp_get_thread_num
!$    external omp_get_max_threads
!$    external omp_get_thread_num
!$    integer :: numThreads
     integer :: myUnit
     NPES_WORLD = 1                  ! default NPES = 1 for serial run

#ifdef USE_ESMF
     Call ESMF_Initialize(vm=modelE_vm, rc=rc)
     Call ESMF_VMGet(modelE_vm, localPET=my_pet, petCount=NPES_WORLD, &
    &     rc=rc)
!$    call omp_set_dynamic(.false.)
!$    numThreads = omp_get_max_threads()
!$    call mpi_bcast(numThreads, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, rc)
!$    call omp_set_num_threads(numThreads)
     if(my_pet == 0) then
        write(*,*)'Num MPI Processes: ', NPES_WORLD
!$       write(*,*)'   with ',numThreads, 'threads per process'
     end if
!CC why create gridded component if not used?
     compmodelE  = ESMF_GridCompCreate(name="ModelE ESMF", rc=rc)
!CC ditto
     ESMF_Layout = ESMF_DELayoutCreate(modelE_vm, &
    &     petList = (/ 1, NPES_WORLD /))
#endif

     NP_LON   = 1
     RANK_LON = 0
     NP_LAT   = NPES_WORLD
     RANK_LAT = my_pet
 
     return
   end subroutine init_app

! ----------------------------------------------------------------------
#ifdef USE_ESMF
   subroutine INIT_GRID(distGrid,IM,JM, LM,width,vm,J_SCM,bc_periodic, &
  &                     CREATE_CAP,npes_max,excess_on_pe0)
     USE ESMF_CUSTOM_MOD, Only : modelE_vm
#else
   SUBROUTINE INIT_GRID(distGrid,IM,JM,LM,width,J_SCM,bc_periodic, &
  &                     CREATE_CAP,npes_max,excess_on_pe0)
#endif
! ----------------------------------------------------------------------
     USE FILEMANAGER, Only : openunit
     IMPLICIT NONE
     TYPE (DIST_GRID), INTENT(INOUT) :: distGrid
     INTEGER, INTENT(IN) :: IM, JM,LM
     INTEGER, OPTIONAL, INTENT(IN) :: J_SCM ! single column model
     INTEGER, OPTIONAL :: width
     LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic
     LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP
     INTEGER, OPTIONAL, INTENT(IN) :: npes_max
     LOGICAL, OPTIONAL, INTENT(IN) :: excess_on_pe0
     integer, parameter :: numDims=2
#ifdef USE_ESMF
     TYPE (ESMF_VM), INTENT(IN), Target, Optional :: vm
#endif
     integer, dimension(numDims) :: grid_size
     integer             :: rc
     real(ESMF_KIND_R8), dimension(numDims) :: range_min,range_max

     INTEGER :: RANK_LON, RANK_LAT
     INTEGER :: J_EQUATOR
     INTEGER :: I0_DUM, I1_DUM
     INTEGER :: J0_DUM, J1_DUM
     INTEGER :: width_
     INTEGER :: pet
     INTEGER :: NTILES
#ifdef USE_ESMF
     TYPE(ESMF_Grid), external :: AppGridCreateF
     TYPE(ESMF_Config) :: cf
     TYPE(ESMF_VM), Pointer :: vm_
     Type (ESMF_DELayout)::layout
     REAL*8 :: deltaZ
     INTEGER :: L
     integer, allocatable            :: IMS(:), JMS(:)
#endif
#ifdef USE_MPI
     Type (ESMF_Axisindex), Pointer :: AI(:,:)
#endif
     INTEGER :: p, pindex
     integer :: npes_used
     integer, dimension(:), allocatable :: pelist
     integer :: group_world,group_used,ierr
     integer :: newCommunicator
     
     logical :: useCubedSphere
     integer :: minIndex(3), maxIndex(3), indexArray(2,3)
     type(ESMF_DistGrid) :: tmpDistGrid
     integer :: lbnd(3) 
     integer :: ubnd(3) 
     real(ESMF_KIND_R8), pointer :: xcoord(:)
     real(ESMF_KIND_R8), pointer :: zcoord(:)
     real(ESMF_KIND_R8),  pointer :: ycoord(:)


#ifdef CUBED_SPHERE
     useCubedSphere = .true.
#else
     useCubedSphere = .false.
#endif

     if (useCubedSphere) then
        grid_size = (/IM, JM*6/)
!!$$        distGrid%localSubdomain = newDomain(IM, 6*JM)
     else
        grid_size = (/IM, JM/)
!!$$        distGrid%localSubdomain = newDomain(IM, JM)
     end if

     range_min(1)=0.;   range_min(2)=-90.
     range_max(1)=360.; range_max(2)=90.

     npes_used = 1

#ifdef USE_ESMF

     Allocate(distGrid%dj_map(0:npes_world-1))

     vm_ => modelE_vm
     If (Present(vm)) vm_ => vm

!CC In the old code but still here: A new layout is created - why?

     layout = ESMF_DELayoutCreate(vm=vm, &
    &     petList = (/ 1, NPES_WORLD /), rc=rc)

!     call ESMF_DELayoutPrint(layout)

    ! Distribute the horizontal dimensions
    !-------------------------------------

     allocate(ims(0:0), jms(0:npes_world-1))
     npes_used = min(npes_world, jm-2) ! jm-2 is latlon-specific
     if (present(npes_max)) npes_used = min(npes_max, npes_used)
     jms(0:npes_used-1) = getLatitudeDistribution(jm, npes_used)
     if(npes_used<npes_world) jms(npes_used:npes_world-1) = 0
     ims(0) = im

!     print *,'ims: ',ims
!     print *,'jms: ',jms
!     print *,' npes_world = ',npes_world,' --- npes_used = ',npes_used

     minIndex = (/1,1,1/)
     maxIndex = (/IM,JM,LM/)

     distGrid%ESMF_GRID = ESMF_GridCreateShapeTile(    &
            name="modelE grid",            &
            countsPerDEDim1=ims,           &
            countsPerDEDim2=jms,           &
            indexFlag = ESMF_INDEX_USER,   &
            gridMemLBound = (/1,1/),       &
            gridEdgeLWidth = (/0,0/),      &
            gridEdgeUWidth = (/0,0/),      &
            coordDep1 = (/1,2/),           &
            coordDep2 = (/1,2/),           &
            rc=rc)

! Allocate coords at default stagger location
!    call ESMF_GridAddCoord(distGrid%ESMF_GRID, rc=rc)

!    call ESMF_AttributeSet(distGrid%ESMF_GRID, name='GRID_LM', &
!         value=LM, rc=rc)

     RANK_LON=0
     RANK_LAT=my_pet
     Call ESMF_GRID_BOUNDS(distGrid, RANK_LON, RANK_LAT, &
    &        I0_DUM, I1_DUM, J0_DUM, J1_DUM)
     write(*,*)'esmf-bounds',my_pet,I0_DUM, I1_DUM, J0_DUM, J1_DUM

#else
     RANK_LON = 0
     RANK_LAT = 0
     I0_DUM = 1
     I1_DUM = IM
     J0_DUM = 1
     J1_DUM = JM

     if (present(J_SCM)) then
        J0_DUM = J_SCM
        J1_DUM = J_SCM
     end if
#endif

     width_ = HALO_WIDTH
     If (Present(width)) width_=width

     distGrid%IM_WORLD      = IM
     distGrid%JM_WORLD      = JM

     ! Wrapped ESMF grid
     distGrid%I_STRT        = I0_DUM
     distGrid%I_STOP        = I1_DUM
     distGrid%I_STRT_HALO   = MAX( 1, I0_DUM-width_)
     distGrid%I_STOP_HALO   = MIN(IM, I1_DUM+width_)
     distGrid%ni_loc = (RANK_LAT+1)*IM/NPES_used - RANK_LAT*IM/NPES_used
     
     distGrid%j_strt        = J0_DUM
     distGrid%J_STOP        = J1_DUM
#ifdef USE_ESMF
!CC original code call routine but does not use layout - 
!it simply call MPI_Barrier
!     call ESMF_DELayoutBarrier(layout, rc)
!CC We can replace by the following:
     call ESMF_VMBarrier(modelE_vm, rc=rc)
#endif
     distGrid%HAVE_DOMAIN   = J0_DUM <= JM

     distGrid%J_STRT_SKP = max (   2, J0_DUM)
     distGrid%J_STOP_SKP = min (JM-1, J1_DUM)

#ifdef USE_MPI
     distGrid%J_STRT_HALO   = J0_DUM - width_
     distGrid%J_STOP_HALO   = J1_DUM + width_
     distGrid%private%numProcesses = npes_used
     distGrid%private%numAllProcesses = npes_world
     distGrid%private%mpi_tag = 10  ! initial value

! Create a new MPI communicator including all the PEs with a nonzero
! domain size.  Even when NPES_USED == NPES_WORLD, this is convenient for
! avoiding collisions of MPI tag sequences.

     call mpi_comm_group(MPI_COMM_WORLD,group_world,ierr)
     allocate(pelist(0:npes_used-1))
     do p=0,npes_used-1
        pelist(p) = p
     enddo
     call mpi_group_incl(group_world,npes_used,pelist,group_used,ierr)
     deallocate(pelist)
     call mpi_comm_create(MPI_COMM_WORLD,group_used, newCommunicator, ierr)
     if(.not. distGrid%HAVE_DOMAIN) newCommunicator = MPI_COMM_NULL
     call setMpiCommunicator(distGrid, newCommunicator)
#else
     ! I guess we don't need HALO in SCM mode...
     !distGrid%J_STRT_HALO = MAX(1,  distGrid % J_STRT - 1)
     !distGrid%J_STOP_HALO = MIN(JM, distGrid % J_STOP + 1)
     distGrid%J_STRT_HALO = MAX(1,  distGrid % J_STRT)
     distGrid%J_STOP_HALO = MIN(JM, distGrid % J_STOP)
#endif
     
     distGrid%J_STRT_STGR   = max(2,J0_DUM)
     distGrid%J_STOP_STGR   = J1_DUM
     
     distGrid%private%hasSouthPole = J0_DUM == 1  !(RANK_LAT == 0)
     distGrid%private%hasNorthPole = J1_DUM == JM .and. J0_DUM <= JM !(RANK_LAT == NP_LAT - 1) &

     J_EQUATOR = JM/2
     distGrid%private%hasEquator =  (J0_DUM <= J_EQUATOR) .AND. (J1_DUM >= J_EQUATOR)

#ifdef USE_DD2D_UTILS
! need to initialize the dd2d version of dist_grid for I/O
     call init_dist_grid( &
    &     distGrid%IM_WORLD,distGrid%JM_WORLD,1,  &
    &     distGrid%I_STRT,distGrid%I_STOP, &
    &     distGrid%j_strt,distGrid%J_STOP, &
    &     distGrid%I_STRT_HALO,distGrid%I_STOP_HALO, &
    &     distGrid%J_STRT_HALO,distGrid%J_STOP_HALO,distGrid)
#endif

     if (present(J_SCM)) then
        ! assume J_SCM is in "general position"
        distGrid%private%hasSouthPole = .false.
        distGrid%private%hasNorthPole = .false.
        distGrid%private%hasEquator    = .false.
     endif

     distGrid % private%PERIODICBC = isPeriodic(bc_periodic)

     ! set lookup table PET(J)
     Allocate(distGrid%private%lookup_pet(1:JM))
     distGrid%private%lookup_pet(:) = 0

#ifdef USE_MPI

     Allocate(AI(NPES_WORLD,3))
     call ESMF_GridGetAxisIndex(distGrid%ESMF_GRID, AI, my_pet)

      write(*,23), " DE : ",my_pet, &
           "  (",AI(my_pet+1,1)%min, " , ",AI(my_pet+1,1)%max, &
         ") - (",AI(my_pet+1,2)%min, " , ",AI(my_pet+1,2)%max, &
         ") - (",AI(my_pet+1,3)%min, " , ",AI(my_pet+1,3)%max
23  format(7(a,i3))

! original:
     Do p = 1, npes_world
       distGrid%private%lookup_pet( AI(p,2)%min : AI(p,2)%max ) = p-1
     end do

     Do p = 0, npes_world-1
       distGrid%dj_map(p) = AI(p+1,2)%max - AI(p+1,2)%min + 1
     end do
     distGrid%dj=distGrid%dj_map(my_pet)

     Deallocate(AI)

#endif

   contains

! ----------------------------------------------------------------------
     function getLatitudeDistribution(jm, numProcesses) result(latsPerProcess)
! ----------------------------------------------------------------------
! Contstraint: assumes jm >=4.
       integer, intent(in) :: jm
       integer :: numProcesses
       integer :: latsPerProcess(0:numProcesses-1)
       
       integer :: excess, npes_used, p
       integer :: localAdjustment

       latsPerProcess = 0

       ! Set minimum requirements per processor
       ! Currently this is 1 lat/proc away from poles
       ! and 2 lat/proc at poles
       select case (npes_world)
       case (1)
          latsPerProcess = jm
          return
       case (2)
          latsPerProcess(0) = jm/2
          latsPerProcess(1) = jm - (jm/2)
          return
       case (3:)
          npes_used = min(numProcesses, jm-2)
          
          ! 1st cut - round down
          latsPerProcess(0:numProcesses-1) = JM/npes_used  ! round down
          
          ! Fix at poles
          latsPerProcess(0) = max(2, latsPerProcess(0))
          latsPerProcess(numProcesses-1) = max(2, latsPerProcess(numProcesses-1))
          
          ! redistribute excess
          excess = JM - sum(latsPerProcess(0:numProcesses-1))
          if(present(excess_on_pe0)) then
             if(excess_on_pe0) then
                latsPerProcess(0) = latsPerProcess(0) + excess
                excess = 0
             endif
          endif
          ! redistribute any remaining excess among interior processors
          do p = 1, numProcesses - 2
             localAdjustment = (p+1)*excess/(numProcesses-2) - (p*excess)/(numProcesses-2)
             latsPerProcess(p) = latsPerProcess(p) + localAdjustment
          end do
       end select
       
     end function getLatitudeDistribution

   END SUBROUTINE INIT_GRID

! ----------------------------------------------------------------------
   SUBROUTINE DESTROY_GRID(distGrid)
! ----------------------------------------------------------------------
     TYPE (DIST_GRID), INTENT(INOUT) :: distGrid
#ifdef USE_ESMF
     Call ESMF_GridDestroy(distGrid%ESMF_Grid)
#endif
   END SUBROUTINE DESTROY_GRID

! ----------------------------------------------------------------------
   SUBROUTINE GET(distGrid, I_STRT, I_STOP, &
  &                        I_STRT_HALO, I_STOP_HALO, &
  &                        J_STRT, J_STOP, J_STRT_HALO, J_STOP_HALO, &
  &                        J_STRT_SKP, J_STOP_SKP, &
  &                        J_STRT_STGR, J_STOP_STGR, &
  &                        have_south_pole, have_north_pole)
! ----------------------------------------------------------------------
     TYPE (DIST_GRID), INTENT(IN) :: distGrid
     INTEGER, OPTIONAL :: I_STRT, I_STOP
     INTEGER, OPTIONAL :: I_STRT_HALO, I_STOP_HALO
     INTEGER, OPTIONAL :: J_STRT, J_STOP
     INTEGER, OPTIONAL :: J_STRT_HALO, J_STOP_HALO
     INTEGER, OPTIONAL :: J_STRT_SKP, J_STOP_SKP
     INTEGER, OPTIONAL :: J_STRT_STGR, J_STOP_STGR
     LOGICAL, OPTIONAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

     IF (PRESENT(I_STRT)) I_STRT = distGrid%I_STRT
     IF (PRESENT(I_STOP)) I_STOP = distGrid%I_STOP
     
     IF (PRESENT(I_STRT_HALO)) I_STRT_HALO = distGrid%I_STRT_HALO
     IF (PRESENT(I_STOP_HALO)) I_STOP_HALO = distGrid%I_STOP_HALO
     
     IF (PRESENT(J_STRT)) J_STRT = distGrid%j_strt
     IF (PRESENT(J_STOP)) J_STOP = distGrid%J_STOP

     IF (PRESENT(J_STRT_HALO)) J_STRT_HALO = distGrid%J_STRT_HALO
     IF (PRESENT(J_STOP_HALO)) J_STOP_HALO = distGrid%J_STOP_HALO
     
     IF (PRESENT(J_STRT_SKP)) J_STRT_SKP = distGrid%J_STRT_SKP
     IF (PRESENT(J_STOP_SKP)) J_STOP_SKP = distGrid%J_STOP_SKP
     
     IF (PRESENT(J_STRT_STGR)) J_STRT_STGR = distGrid%J_STRT_STGR
     IF (PRESENT(J_STOP_STGR)) J_STOP_STGR = distGrid%J_STOP_STGR
     
     IF (PRESENT(HAVE_SOUTH_POLE)) &
    &             HAVE_SOUTH_POLE= distGrid%private%hasSouthPole
     IF (PRESENT(HAVE_NORTH_POLE)) &
    &             HAVE_NORTH_POLE= distGrid%private%hasNorthPole

   END SUBROUTINE GET


! ----------------------------------------------------------------------
   SUBROUTINE FINISH_APP()
! ----------------------------------------------------------------------
     USE FILEMANAGER, ONLY : closeunit
     IMPLICIT NONE
     
     INTEGER :: rc

#ifdef DEBUG_DECOMP
     CALL closeunit(CHECKSUM_UNIT)
     CALL closeunit(grid%private%log_unit)
#endif

#ifdef USE_ESMF
     CALL ESMF_FINALIZE(rc=rc)
#endif

   END SUBROUTINE FINISH_APP

#ifdef USE_ESMF

! ----------------------------------------------------------------------
   subroutine ESMF_GRID_BOUNDS(GRID, X_LOC, Y_LOC, I1, IN, J1, JN)
! ----------------------------------------------------------------------
     type (DIST_Grid), intent(INOUT) :: grid
     integer, intent(IN)          :: X_LOC, Y_LOC
     integer, intent(OUT)         :: I1, IN, J1, JN
     integer :: i
     integer :: rc

     type(ESMF_AxisIndex), dimension(:,:), pointer :: AI

     Allocate(AI(0:NPES_WORLD-1,3))
     call ESMF_GridGetAxisIndex(grid%ESMF_GRID, AI, my_pet)

     I1 = AI(my_pet,1)%min
     IN = AI(my_pet,1)%max
     J1 = AI(my_pet,2)%min
     JN = AI(my_pet,2)%max

     do i=0, npes_world-1
        grid%dj_map(i)=AI(i,2)%max - AI(i,2)%min + 1
     end do
     grid%dj=grid%dj_map(my_pet)

     Deallocate(AI)

   end subroutine ESMF_GRID_BOUNDS

! ----------------------------------------------------------------------
   function load_cap_config(config_file,IM,JM,LM,NP_X,NP_Y)  &
  &         result( config )
! ----------------------------------------------------------------------
     use ESMF_mod
     use FILEMANAGER    
     character(len=*), parameter :: Iam= &
    &               "DOMAIN_DECOMP::load_cap_config"
     character(len=*), intent(in) :: config_file
     integer,          intent(in) :: IM,JM,LM,NP_X,NP_Y
     type (esmf_config)           :: config
    
     integer :: rc, iunit
     type (ESMF_VM) :: vm
    
     config = esmf_configcreate(rc=rc)
   
   !     print*, 'load_cap_config: ', IM, JM, LM, NP_X, NP_Y
 
     call openunit(config_file, iunit, qbin=.false., qold=.false.)
     if(am_i_root()) then
        write(iunit,*)'IM:  ', IM
        write(iunit,*)'JM:  ', JM
        write(iunit,*)'LM:  ', LM
        write(iunit,*)'NX:  ', NP_X
        write(iunit,*)'NY:  ', NP_Y
     endif
     call closeUnit(iunit)
     
     Call ESMF_VMGetGlobal(vm, rc)
     call esmF_VMbarrier(vm, rc)
     call esmf_configloadfile(config, config_file, rc=rc)
     
   end function load_cap_config

#endif
      
! ----------------------------------------------------------------------
   function amRootESMF(layout) result(R)
! ----------------------------------------------------------------------
     type (ESMF_DELayout) :: layout
     logical              :: R
     
#ifdef USE_MPI
     R = .false.
     if (my_pet == root) R = .true.
#else
     R = .true.
#endif
     
   end function amRootESMF

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_1D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:) :: arr
      REAL*8, DIMENSION(:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_ESMF
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + arr_tmp
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_1D_I(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, DIMENSION(:) :: arr
      INTEGER, DIMENSION(:), optional :: arr_master
      logical, intent(in), optional :: increment
      INTEGER, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_ESMF
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) allocate(arr_tmp(arr_size))
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_INTEGER, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + arr_tmp
          deallocate(arr_tmp)
        endif
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_INTEGER, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_INTEGER,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D_I

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_2D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:) :: arr
      REAL*8, DIMENSION(:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_ESMF
         arr_size = size(arr)
         if(increment_) then
            if(am_i_root()) then
               allocate(arr_tmp(arr_size))
            else
               allocate(arr_tmp(1))
            end if
            call mpi_reduce(arr,arr_tmp,arr_size, &
    &           MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &           MPI_COMM_WORLD, ierr)
            if(am_i_root()) then
              arr_master = arr_master + reshape(arr_tmp,shape(arr))
            endif
            deallocate(arr_tmp)
         else
            call mpi_reduce(arr,arr_master,arr_size, &
    &           MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &           MPI_COMM_WORLD, ierr)
         endif
#else 
         if(increment_) then
            arr_master = arr_master + arr
         else
            arr_master = arr
         endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_2D

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_3D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:,:) :: arr
      REAL*8, DIMENSION(:,:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_ESMF
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr  is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_3D

! ----------------------------------------------------------------------
      SUBROUTINE SUMXPE_4D(arr, arr_master, increment)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:,:,:) :: arr
      REAL*8, DIMENSION(:,:,:,:), optional :: arr_master
      logical, intent(in), optional :: increment
      REAL*8, DIMENSION(:), ALLOCATABLE :: arr_tmp
      logical :: increment_
      logical :: loc_
      integer :: ierr,arr_size
      if(present(increment)) then
        increment_ = increment
      else
        increment_ = .false.
      endif
      if (present(arr_master)) then
         loc_ = .true.
      else
         loc_ = .false.
         increment_ = .false.
      endif
      if (loc_) then
#ifdef USE_ESMF
      arr_size = size(arr)
      if(increment_) then
        if(am_i_root()) then
           allocate(arr_tmp(arr_size))
        else
           allocate(arr_tmp(1))
        end if
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION, &
    &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
      endif
#else
      if(increment_) then
        arr_master = arr_master + arr
      else
        arr_master = arr
      endif
#endif
      else  
!**** arr plays both roles of local and global array
!**** arr  is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size, &
    &        MPI_DOUBLE_PRECISION,MPI_SUM,root, &
    &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_4D

#ifdef USE_MPI

      subroutine ESMF_GRID_PE_LAYOUT  (GRID, NX, NY)
        type (Dist_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: NX, NY

        NX = 1
        NY = getNumProcesses(grid)

      end subroutine ESMF_GRID_PE_LAYOUT

      subroutine ESMF_GRID_MY_PE_LOC  (GRID,NX0,NY0)
        type (Dist_Grid), intent(IN) :: grid
        integer, intent(OUT)          :: NX0, NY0

        NX0 = 0
        NY0 = my_pet

      end subroutine ESMF_GRID_MY_PE_LOC

#endif

! ----------------------------------------------------------------------
      subroutine ESMF_DELayoutBarrier(layout, rc)
! ----------------------------------------------------------------------
      type (ESMF_DELayout) :: layout
      integer, optional  :: rc

      integer :: ierr

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (present(rc)) then
         rc = ESMF_SUCCESS
      endif
#endif

      end subroutine ESMF_DELayoutBarrier

! ----------------------------------------------------------------------
      SUBROUTINE HERE(file, line)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      Character(Len=*) :: file
      Integer :: line

      INTEGER :: ierr
      Integer, Allocatable :: lines(:)

#ifdef DEBUG_DECOMP
      CALL LOG_PARALLEL(grid, file, line)
      If (AM_I_ROOT()) Then
         WRITE(CHECKSUM_UNIT,*)'HERE: ',file, line
         CALL SYS_FLUSH(CHECKSUM_UNIT)
       End If
#ifdef USE_MPI
       ALLOCATE(lines(npes_world))
       Call MPI_Allgather(line, 1, MPI_INTEGER, lines, 1, MPI_INTEGER, &
    &      MPI_COMM_WORLD, ierr)
       If (Any(lines /= line)) &
    &      call stop_model('HERE: synchronization error -severe.',255)
       Deallocate(lines)
#endif
#endif
#ifdef USE_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE HERE

! ----------------------------------------------------------------------
      SUBROUTINE LOG_PARALLEL(distGrid, file, line, i0, i1, x0, x1)
! ----------------------------------------------------------------------
      Use FILEMANAGER, only : nameunit
      IMPLICIT NONE

      TYPE(DIST_GRID), INTENT(IN) :: distGrid
      CHARACTER(Len=*) :: file
      INTEGER          :: line

      INTEGER, OPTIONAL :: i0, i1(:)
      REAL*8, OPTIONAL  :: x0, x1(:)

      INTEGER :: iu
      INTEGER :: n

#ifdef DEBUG_DECOMP
      iu = distGrid%private%log_unit
      WRITE(iu, *) file, line
      If (PRESENT(i0)) WRITE(iu, *) '   i0=',i0
      If (PRESENT(x0)) WRITE(iu, *) '   x0=',x0
      IF (PRESENT(i1)) THEN
         DO n = 1, Size(i1)
            WRITE(iu, '(10x,i4,1x,a,i6)')n, '   i1=',i1(n)
         END DO
      END IF
      IF (PRESENT(x1)) THEN
         DO n = 1, Size(x1)
            WRITE(iu, '(10x,i4,1x,a,e22.16)')n, '   x1=',x1(n)
         END DO
      END IF
      CALL SYS_FLUSH(iu)
#endif

      END SUBROUTINE LOG_PARALLEL

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMIN_R(distGrid, val, val_min)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_min

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_min, 1, MPI_DOUBLE_PRECISION,MPI_MIN, &
    &     getMpiCommunicator(DISTGRID), ierr)
#else
      val_min = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMAX_R(distGrid, val, val_max)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_DOUBLE_PRECISION,MPI_MAX, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMAX_I(distGrid, val, val_max)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_INTEGER, MPI_MAX, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

! ----------------------------------------------------------------------
      SUBROUTINE GLOBALMAX_I_1D(distGrid, val, val_max)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: distGrid
      INTEGER,            INTENT(IN)  :: val(:)
      INTEGER,            INTENT(OUT) :: val_max(:)

      INTEGER  :: n,ierr

#ifdef USE_MPI
      n = size(val)
      CALL MPI_Allreduce(val, val_max, n, MPI_INTEGER, MPI_MAX, &
    &     getMpiCommunicator(distGrid), ierr)
#else
      val_max(:) = val(:)
#endif

      END SUBROUTINE


! ----------------------------------------------------------------------
      SUBROUTINE INIT_BAND_PACK_TYPE(grd_src, grd_dst, band_j0,band_j1, &
    &     bandpack)
! initialize the bandpack derived type with the information needed
! for the band_pack procedure to fill output arrays with data from
! J indices band_j0 to band_j1
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_src,grd_dst
      INTEGER, INTENT(IN) :: band_j0,band_j1
      TYPE (BAND_PACK_TYPE), intent(OUT) :: bandpack
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: &
    &     j0_have,j1_have,j0_requested,j1_requested
      integer :: p, ierr, im,jm, j0send,j1send,j0recv,j1recv,npes

!
! Store some general MPI info (NOTE: for now we assume that the process
! set of grd_src is a subset of that in grd_dst, or vice versa, which
! allows us to take this info from the larger set).
!
      if(getNumProcesses(grd_src) > getNumProcesses(grd_dst)) then
        npes = getNumProcesses(grd_src)
        bandpack%mpi_comm = getMpiCommunicator(grd_src)
      else
        npes = getNumProcesses(grd_dst)
        bandpack%mpi_comm = getMpiCommunicator(grd_dst)
      endif
      bandpack%npes_comm = npes
      allocate(bandpack%j0_send(0:npes-1), &
    &         bandpack%j1_send(0:npes-1), &
    &         bandpack%j0_recv(0:npes-1), &
    &         bandpack%j1_recv(0:npes-1), &
    &         bandpack%scnts(0:npes-1), &
    &         bandpack%sdspl(0:npes-1), &
    &         bandpack%rcnts(0:npes-1), &
    &         bandpack%rdspl(0:npes-1), &
    &         bandpack%sdspl_inplace(0:npes-1), &
    &         bandpack%rdspl_inplace(0:npes-1) &
    &     )
#endif
!      bandpack%im_world = grd_src%im_world
      bandpack%j_strt = grd_src%j_strt
      bandpack%j_stop = grd_src%j_stop
      bandpack%j_strt_halo = grd_src%j_strt_halo
      bandpack%j_stop_halo = grd_src%j_stop_halo
      bandpack%jband_strt = band_j0
      bandpack%jband_stop = band_j1
#ifdef USE_MPI
      im = 1!grd_src%im_world
      jm = grd_src%jm_world
!
! Set up the MPI send/receive information
!
      call mpi_allgather(grd_src%j_strt,1,MPI_INTEGER,j0_have,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(grd_src%j_stop,1,MPI_INTEGER,j1_have,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(band_j0,1,MPI_INTEGER,j0_requested,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      call mpi_allgather(band_j1,1,MPI_INTEGER,j1_requested,1, &
    &     MPI_INTEGER,bandpack%MPI_COMM,ierr)
      do p=0,npes-1
        j0send = max(grd_src%j_strt,j0_requested(p))
        j1send = min(grd_src%j_stop,j1_requested(p))
        bandpack%j0_send(p) = j0send
        bandpack%j1_send(p) = j1send
        if(j0send <= j1send) then
          bandpack%scnts(p) = im*(j1send-j0send+1)
          bandpack%sdspl_inplace(p) = im*(j0send-grd_src%j_strt_halo)
        else
          bandpack%scnts(p) = 0
          bandpack%sdspl_inplace(p) = 0
        endif
        j0recv = max(j0_have(p),band_j0)
        j1recv = min(j1_have(p),band_j1)
        bandpack%j0_recv(p) = j0recv
        bandpack%j1_recv(p) = j1recv
        if(j0recv <= j1recv) then
          bandpack%rcnts(p) = im*(j1recv-j0recv+1)
          bandpack%rdspl_inplace(p) = im*(j0recv-band_j0)
        else
          bandpack%rcnts(p) = 0
          bandpack%rdspl_inplace(p) = 0
        endif
      enddo
      bandpack%rdspl(0) = 0
      bandpack%sdspl(0) = 0
      do p=1,npes-1
        bandpack%sdspl(p) = bandpack%sdspl(p-1)+bandpack%scnts(p-1)
        bandpack%rdspl(p) = bandpack%rdspl(p-1)+bandpack%rcnts(p-1)
      enddo
#endif
      RETURN
      END SUBROUTINE INIT_BAND_PACK_TYPE

! ----------------------------------------------------------------------
      SUBROUTINE BAND_PACK_ij(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr
      integer :: im,npes
      npes = bandpack%npes_comm
      im = size(arr,1)
      scnts(0:npes-1) = im*bandpack%scnts
      sdspl(0:npes-1) = im*bandpack%sdspl_inplace
      rcnts(0:npes-1) = im*bandpack%rcnts
      rdspl(0:npes-1) = im*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision, &
    &                   arr_band, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
#else
      arr_band(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP) = &
    &     arr(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ij

! ----------------------------------------------------------------------
      SUBROUTINE BAND_PACK_ijl(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:,:)
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,im,lm,i,j,l,n,p
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer :: npes
      npes = bandpack%npes_comm
      im = size(arr,1)
      lm = size(arr,3)
      scnts = im*lm*bandpack%scnts
      sdspl = im*lm*bandpack%sdspl
      rcnts = im*lm*bandpack%rcnts
      rdspl = im*lm*bandpack%rdspl
      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_send(p),bandpack%j1_send(p)
            do i=1,im
              n = n + 1
              bufsend(n) = arr(i,j,l)
            enddo
          enddo
        enddo
      enddo
      call mpi_alltoallv(bufsend, scnts, sdspl, mpi_double_precision, &
    &                   bufrecv, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_recv(p),bandpack%j1_recv(p)
            do i=1,im
              n = n + 1
              arr_band(i,j,l) = bufrecv(n)
            enddo
          enddo
        enddo
      enddo
      deallocate(bufsend,bufrecv)
#else
      arr_band(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP,:) = &
    &     arr(:,bandpack%JBAND_STRT:bandpack%JBAND_STOP,:)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ijl

! ----------------------------------------------------------------------
      SUBROUTINE BAND_PACK_COLUMN(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, dimension(0:npes_world-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,im,lm
      integer :: npes
      npes = bandpack%npes_comm
      im = size(arr,2)
      lm = size(arr,1)
      scnts(0:npes-1) = im*lm*bandpack%scnts
      sdspl(0:npes-1) = im*lm*bandpack%sdspl_inplace
      rcnts(0:npes-1) = im*lm*bandpack%rcnts
      rdspl(0:npes-1) = im*lm*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision, &
    &                   arr_band, rcnts, rdspl, mpi_double_precision, &
    &                   bandpack%mpi_comm, ierr)
#else
      arr_band(:,:,bandpack%JBAND_STRT:bandpack%JBAND_STOP) = &
    &     arr(:,:,bandpack%JBAND_STRT:bandpack%JBAND_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_COLUMN

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_BCAST_0D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_0D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_BCAST_1D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_1D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_BCAST_2D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_2D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_BCAST_3D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_3D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_BCAST_4D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(InOut) :: arr(:,:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif

      END SUBROUTINE ESMF_BCAST_4D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_IBCAST_0D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_INTEGER,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_0D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_IBCAST_1D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_1D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_IBCAST_2D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_2D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_IBCAST_3D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_3D

! ----------------------------------------------------------------------
      SUBROUTINE ESMF_IBCAST_4D(distGrid, arr)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(InOut) :: arr(:,:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root, &
    &     getMpiCommunicator(distGrid), ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_4D


! ----------------------------------------------------------------------
      SUBROUTINE TRANSPOSE_ijk(grid, x_in, x_out, reverse)
! ----------------------------------------------------------------------
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_out(:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES_WORLD-1), I1(0:NPES_WORLD-1)
      INTEGER :: J0(0:NPES_WORLD-1), J1(0:NPES_WORLD-1)
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,nk,k
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER, DIMENSION(0:NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X_IN(:,1:grid%JM_WORLD,:) = X_OUT(:,1:grid%JM_WORLD,:)
      Else
         X_OUT(:,1:grid%JM_WORLD,:) = X_IN(:,1:grid%JM_WORLD,:)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes_world-1,3))
      call ESMF_GridGetAxisIndex(grid%ESMF_GRID, AI, my_pet)
      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(my_pet) - I0(my_pet) + 1
      nj_loc = J1(my_pet) - J0(my_pet) + 1

      nk = SIZE(X_IN,3)

      ALLOCATE(rbuf(grid%JM_WORLD * ni_loc * nk))
      ALLOCATE(sbuf(grid%IM_WORLD * nj_loc * nk))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(icnt) = X_out(i,j,k)
               END DO
            END DO
         ELSE
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(icnt) = X_in(i,j,k)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = nj_loc * nip * nk
         rcnts(p) = ni_loc * njp * nk
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
       END DO
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X_in(i,j,k) = sbuf(icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_out(i,j,k) = rbuf(icnt)
               END DO
            END DO
         End If
       END DO
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_ijk

! ----------------------------------------------------------------------
      SUBROUTINE TRANSPOSE_ij(grid, x_in, x_out, reverse)
! ----------------------------------------------------------------------
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:)
      REAL*8 :: x_out(:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES_WORLD-1), I1(0:NPES_WORLD-1)
      INTEGER :: J0(0:NPES_WORLD-1), J1(0:NPES_WORLD-1)
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER, DIMENSION(0:NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X_IN(:,1:grid%JM_WORLD) = X_OUT(:,1:grid%JM_WORLD)
      Else
         X_OUT(:,1:grid%JM_WORLD) = X_IN(:,1:grid%JM_WORLD)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes_world-1,3))
      call ESMF_GridGetAxisIndex(grid%ESMF_GRID, AI, my_pet)
      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(my_pet) - I0(my_pet) + 1
      nj_loc = J1(my_pet) - J0(my_pet) + 1


      ALLOCATE(rbuf(grid%JM_WORLD * ni_loc ))
      ALLOCATE(sbuf(grid%IM_WORLD * nj_loc ))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(icnt) = X_out(i,j)
               END DO
            END DO
         ELSE
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(icnt) = X_in(i,j)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = nj_loc * nip
         rcnts(p) = ni_loc * njp
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
         If (reverse_) Then
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X_in(i,j) = sbuf(icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_out(i,j) = rbuf(icnt)
               END DO
            END DO
         End If
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_ij

! ----------------------------------------------------------------------
      SUBROUTINE TRANSPOSE_COLUMN(grid, x, x_tr, reverse)
! ----------------------------------------------------------------------
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x(:,:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_tr(:,:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES_WORLD-1), I1(0:NPES_WORLD-1)
      INTEGER :: J0(0:NPES_WORLD-1), J1(0:NPES_WORLD-1)
      REAL*8, ALLOCATABLE :: sbuf(:,:), rbuf(:,:)
#ifdef USE_MPI
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,k
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt, npes
      INTEGER, DIMENSION(0:NPES_WORLD-1) :: scnts, rcnts, sdspl, rdspl
      INTEGER :: n, nk
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_MPI
      If (reverse_) Then
         X(:,:,1:grid%JM_WORLD,:) = X_TR(:,:,1:grid%JM_WORLD,:)
      Else
         X_TR(:,:,1:grid%JM_WORLD,:) = X(:,:,1:grid%JM_WORLD,:)
      End If
#else
      npes = getNumProcesses(grid)

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes_world-1,3))

      call ESMF_GridGetAxisIndex(grid%ESMF_GRID, AI, my_pet)
      DO p = 0, npes - 1
         J0(p) = AI(p,2)%min
         J1(p) = AI(p,2)%max
      END DO
      DEALLOCATE(AI)

      ni_loc = I1(my_pet) - I0(my_pet) + 1
      nj_loc = J1(my_pet) - J0(my_pet) + 1

      n  = SIZE(X, 1)
      nk = SIZE(X,4)

      ALLOCATE(rbuf(n, grid%JM_WORLD * ni_loc * nk))
      ALLOCATE(sbuf(n, grid%IM_WORLD * nj_loc * nk))

      sdspl(0) = 0
      rdspl(0) = 0
      icnt = 0
      DO p = 0, npes -1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  rbuf(:,icnt) = X_tr(:,i,j,k)
               END DO
            END DO
         ELSE
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  sbuf(:,icnt) = X(:,i,j,k)
               End Do
            END DO
         END IF
         nip = I1(p) - I0(p) + 1
         njp = J1(p) - J0(p) + 1
         scnts(p) = n * nj_loc * nip * nk
         rcnts(p) = n * ni_loc * njp * nk
         If (p > 0) sdspl(p) = sdspl(p-1) + scnts(p-1)
         If (p > 0) rdspl(p) = rdspl(p-1) + rcnts(p-1)
       END DO
      END DO

      If (reverse_) Then
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION, &
    &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION, &
    &        getMpiCommunicator(grid), ierr)
      End If

      icnt = 0
      DO p = 0, npes - 1
        Do k = 1, nk
         If (reverse_) Then
            DO j = J0(my_pet), J1(my_pet)
               DO i = I0(p), I1(p)
                  icnt = icnt + 1
                  X(:,i,j,k) = sbuf(:,icnt)
               End Do
            END DO
         Else
            DO j = J0(p), J1(p)
               DO i = 1, ni_loc
                  icnt = icnt + 1
                  X_tr(:,i,j,k) = rbuf(:,icnt)
               END DO
            END DO
         End If
       END DO
      END DO

      DEALLOCATE(sbuf)
      DEALLOCATE(rbuf)
#endif

      END SUBROUTINE TRANSPOSE_COLUMN

! ----------------------------------------------------------------------
      logical function haveLatitude(distGrid, j)
! ----------------------------------------------------------------------
      type (DIST_GRID), intent(in) :: distGrid
      integer, intent(in) :: j

      haveLatitude = (j >= distGrid%j_strt .and. j <= distGrid%J_STOP)

      end function haveLatitude

! ----------------------------------------------------------------------
      subroutine SEND_TO_J_1D(distGrid, arr, j_dest, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, Size(arr), MPI_DOUBLE_PRECISION, &
    &     distGrid%private%lookup_pet(j_dest), tag, MPI_COMM_WORLD, ierr)
#endif
      end subroutine SEND_TO_J_1D

! ----------------------------------------------------------------------
      subroutine ISEND_TO_J_0D(distGrid, arr, j_dest, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, 1, MPI_INTEGER, &
    &     distGrid%private%lookup_pet(j_dest), tag, MPI_COMM_WORLD, ierr)
#endif
      end subroutine ISEND_TO_J_0D

! ----------------------------------------------------------------------
      subroutine RECV_FROM_J_1D(distGrid, arr, j_src, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, Size(arr), MPI_DOUBLE_PRECISION, &
    &     distGrid%private%lookup_pet(j_src), tag, MPI_COMM_WORLD, status, ierr)
#endif
      end subroutine RECV_FROM_J_1D

! ----------------------------------------------------------------------
      subroutine IRECV_FROM_J_0D(distGrid, arr, j_src, tag)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: distGrid
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, 1, MPI_INTEGER, &
    &     distGrid%private%lookup_pet(j_src), tag, MPI_COMM_WORLD, status, ierr)
#endif
      end subroutine IRECV_FROM_J_0D

! ----------------------------------------------------------------------
! Helper function to handle optional arguments related to periodic boundaries
      logical function isPeriodic(override)
! ----------------------------------------------------------------------
        logical, optional, intent(in) :: override

        isPeriodic = .false.
        if (present(override)) isPeriodic = override

      end function isPeriodic

! ----------------------------------------------------------------------
      subroutine setMpiCommunicator(this, comm)
! ----------------------------------------------------------------------
        type (dist_grid), intent(inout) :: this
        integer, intent(in) :: comm
        this%private%MPI_COMM = comm
      end subroutine setMpiCommunicator

      integer function getMpiCommunicator(this) result(comm)
        type (dist_grid), intent(in) :: this
        comm = this%private%mpi_comm
      end function getMpiCommunicator

! ----------------------------------------------------------------------
      integer function getNumProcesses(this) result(numProcesses)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        numProcesses = this%private%numProcesses
      end function getNumProcesses

! ----------------------------------------------------------------------
      integer function getNumAllProcesses(this) result(numAllProcesses)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        numAllProcesses = this%private%numAllProcesses
      end function getNumAllProcesses


! ----------------------------------------------------------------------
      subroutine incrementMpiTag(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(inout) :: this
        integer, parameter :: MIN_TAG = 10
        integer, parameter :: MAX_TAG = 128

        integer :: tag
        tag = this%private%mpi_tag
        this%private%mpi_tag = max(mod(tag,MAX_TAG),MIN_TAG) + 1
      end subroutine incrementMpiTag

! ----------------------------------------------------------------------
      integer function getMpiTag(this) result(mpiTag)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        mpiTag = this%private%mpi_tag
      end function getMpiTag

! ----------------------------------------------------------------------
      logical function hasSouthPole(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        hasSouthPole = this%private%hasSouthPole
      end function hasSouthPole

! ----------------------------------------------------------------------
      logical function hasNorthPole(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        hasNorthPole = this%private%hasNorthPole
      end function hasNorthPole

! ----------------------------------------------------------------------
      logical function hasPeriodicBC(this)
! ----------------------------------------------------------------------
        type (dist_grid), intent(in) :: this
        hasPeriodicBC = this%private%periodicBC
      end function hasPeriodicBC

! ----------------------------------------------------------------------
      integer function getLogUnit()
! ----------------------------------------------------------------------
        use FileManager, only: openUnit
        integer :: unit
        character(len=40) :: logFileName

        integer, parameter :: UNINITIALIZED = -1
        integer, save :: logUnit = UNINITIALIZED

        if (logUnit == UNINITIALIZED) then
          write(logFileName,'(a,i4.4)') 'debug.', my_pet
          call openUnit(logFileName, logUnit, qbin=.false., qold=.false.)
        end if

        getLogUnit = logUnit

      end function getLogUnit

END MODULE dist_grid_mod

