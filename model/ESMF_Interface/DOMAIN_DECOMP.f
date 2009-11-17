#include "rundeck_opts.h"
#ifdef MPI_DEFS_HACK
#include "mpi_defs.h"
#endif

#ifndef USE_FVCUBED
#define DOMAIN_DECOMP_ATM_IS_1D
#endif

#ifdef NEW_IO
#define USE_DD2D_UTILS
#endif

      MODULE DOMAIN_DECOMP_1D
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@+    for the message passing (ESMF) implementation.
!@auth NCCS ASTG

#define FILL(N) IAND(USABLE_FROM,N)==N

#if ( defined USE_ESMF )  || ( defined USE_MPP )
#define USE_MPI
#endif


#ifdef USE_ESMF
      use ESMF_Mod
#endif
      use ESMF_CUSTOM_MOD, Only: NORTH, SOUTH

c retaining for now, but disabling, the MPP+FVCUBED coding in this file
#undef USE_MPP
#ifdef USE_FVCUBED
#define USE_DD2D_UTILS
#undef USE_FVCUBED
#endif

#ifdef USE_MPP
      use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe
      use mpp_mod,         only : mpp_error, NOTE, FATAL
      use mpp_domains_mod, only : mpp_domains_init, MPP_DOMAIN_TIME
      use mpp_domains_mod, only : mpp_domains_set_stack_size
      use mpp_domains_mod, only : mpp_define_layout, mpp_define_mosaic
      use mpp_domains_mod, only : domain2D
      use mpp_domains_mod, only : mpp_get_compute_domain
      use mpp_domains_mod, only : mpp_get_data_domain
      use mpp_domains_mod, only : mpp_update_domains
      use mpp_parameter_mod, only : WUPDATE, EUPDATE, SUPDATE, NUPDATE
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
#ifdef USE_MPP
      type ESMF_AXISINDEX
         sequence
         integer :: min
         integer :: max
         integer :: stride
      end type ESMF_AXISINDEX
      integer, parameter :: ESMF_SUCCESS = 1
#endif
#endif

#ifdef USE_MPP
      public init_domain, UPDATE_HALO_2D
#endif

!aoo since DIST_GRID is public ESMF_GRID has to be public
!aoo (SGI compiler complains)
      PUBLIC :: ESMF_GRID

#ifdef USE_ESMF
      public :: load_cap_config
#endif

      TYPE(ESMF_GridComp)  :: compmodelE
      TYPE (ESMF_DELayout) :: ESMF_LAYOUT

!@var DIST_GRID derived type to provide ESMF decomposition info
!@+   public components are used to minimize overhead for accessing
!@+   routine components
      PUBLIC :: DIST_GRID
!@var  grid Default decomposition; globally accessible for convenience.
      PUBLIC :: grid
!@var INIT_APP set some parameters and initialize ESMF
      PUBLIC :: INIT_APP
      PUBLIC :: INIT_GRID
      PUBLIC :: DESTROY_GRID
!@var FINISH_APP Cleans up at the end of the run (closes debugging file)
      PUBLIC :: FINISH_APP
!@var HALO_UPDATE Update data in halo for local domain using data from
!@+   neighbouring processes
      PUBLIC :: HALO_UPDATE ! Communicate overlapping portions of subdomains
      PUBLIC :: HALO_UPDATEj ! jx
      PUBLIC :: HALO_UPDATE_COLUMN ! K, I, J
      PUBLIC :: HALO_UPDATE_BLOCK ! K, L, I, J
      PUBLIC :: HALO_UPDATE_MASK ! K, L, I, J
!@var CHECKSUM output a bit-reproducible checksum for an array
      PUBLIC :: CHECKSUM ! Communicate overlapping portions of subdomains
      PUBLIC :: CHECKSUMj! Communicate overlapping portions of subdomains
      PUBLIC :: CHECKSUM_COLUMN ! K, I, J
!@var GLOBALSUM output a bit-reproducible global-hemisphere-zonal sum for an array
      PUBLIC :: GLOBALSUM
!@var GLOBALMIN determine max value across pes
      PUBLIC :: GLOBALMIN
!@var GLOBALMAX determine max value across pes
      PUBLIC :: GLOBALMAX
!@var SUMXPE sum an array over processors without reducing its rank
      PUBLIC :: SUMXPE
!@var ARRAYSCATTER scatter a global array to a decomposed array
      PUBLIC :: ARRAYSCATTER
!@var ARRAYGATHER gather a decomposed array to a global array
      PUBLIC :: ARRAYGATHER
#ifdef USE_MPI
      PUBLIC :: xESMF_ARRAYGATHER
      PUBLIC :: xESMF_ARRAYSCATTER
#endif
!@var GET - extracts bounds information from DIST_GRID object
      PUBLIC :: GET
      PUBLIC :: HERE
      PUBLIC :: LOG_PARALLEL

      PUBLIC :: BACKSPACE_PARALLEL
      PUBLIC :: REWIND_PARALLEL
      PUBLIC :: SKIP_PARALLEL
      PUBLIC :: DREAD_PARALLEL,DREAD8_PARALLEL,DWRITE8_PARALLEL
      PUBLIC :: MREAD_PARALLEL
      PUBLIC :: READT_PARALLEL,WRITET_PARALLEL,READT_PARALLEL_COLUMN
      PUBLIC :: READT8_PARALLEL,READT8_COLUMN,WRITET8_COLUMN
      PUBLIC :: READ_PARALLEL
      PUBLIC :: WRITE_PARALLEL
      PUBLIC :: WRITEI_PARALLEL
      PUBLIC :: WRITEI8_PARALLEL
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

!@var HALO_UPDATE Generic wrapper for 2D and 3D routines
      INTERFACE HALO_UPDATE
        MODULE PROCEDURE HALO_UPDATE_1D  ! J
        MODULE PROCEDURE HALO_UPDATE_2D  ! I,J
        MODULE PROCEDURE HALO_UPDATE_3D  ! I,J,K
        MODULE PROCEDURE HALO_UPDATE_2Dint  ! I,J,K
      END INTERFACE

      INTERFACE HALO_UPDATEj
      MODULE PROCEDURE HALO_UPDATEj_2d
      ENd INTERFACE HALO_UPDATEj

      INTERFACE HALO_UPDATE_COLUMN
        MODULE PROCEDURE HALO_UPDATE_COLUMN_2D  ! M,J
        MODULE PROCEDURE HALO_UPDATE_COLUMN_3D  ! M,I,J
        MODULE PROCEDURE INT_HALO_UPDATE_COLUMN_3D  ! M,I,J
        MODULE PROCEDURE HALO_UPDATE_COLUMN_4D  ! M,I,J,K
        MODULE PROCEDURE HALO_UPDATE_COLUMN_7D  ! M1,M2,M3,M4,I,J
      END INTERFACE

      INTERFACE HALO_UPDATE_BLOCK
        MODULE PROCEDURE HALO_UPDATE_BLOCK_4D  ! K,L,I,J
      END INTERFACE

      INTERFACE CHECKSUM
        MODULE PROCEDURE CHECKSUM_1D
        MODULE PROCEDURE CHECKSUM_2D
        MODULE PROCEDURE CHECKSUM_3D
      END INTERFACE

      INTERFACE TRANSP
        MODULE PROCEDURE TRANSPOSE_ij
        MODULE PROCEDURE TRANSPOSE_ijk
      END INTERFACE

      INTERFACE CHECKSUMj
        MODULE PROCEDURE CHECKSUMj_2D
        MODULE PROCEDURE CHECKSUMj_3D
      END INTERFACE

      INTERFACE CHECKSUM_COLUMN
        MODULE PROCEDURE CHECKSUM_COLUMN_2D
        MODULE PROCEDURE CHECKSUM_COLUMN_3D
        MODULE PROCEDURE INT_CHECKSUM_COLUMN_3D
        MODULE PROCEDURE CHECKSUM_COLUMN_4D
      END INTERFACE

      INTERFACE GLOBALSUM
        MODULE PROCEDURE GLOBALSUM_INT_REDUCE
        MODULE PROCEDURE GLOBALSUM_J
        MODULE PROCEDURE GLOBALSUM_IJ
        MODULE PROCEDURE GLOBALSUM_IJK
        MODULE PROCEDURE GLOBALSUM_IJK_IK
        MODULE PROCEDURE GLOBALSUM_OTHER_IJK
        MODULE PROCEDURE GLOBALSUM_OTHER_IJK_IK
        MODULE PROCEDURE GLOBALSUM_JK
        MODULE PROCEDURE GLOBALSUM_XXXJ_XXX
        MODULE PROCEDURE GLOBALSUM_XXXIJ_XXX
      END INTERFACE

      INTERFACE SUMXPE
        MODULE PROCEDURE SUMXPE_1D
        MODULE PROCEDURE SUMXPE_1D_I
        MODULE PROCEDURE SUMXPE_2D
        MODULE PROCEDURE SUMXPE_3D
c        MODULE PROCEDURE SUMXPE_4D
c        MODULE PROCEDURE SUMXPE_5D
      END INTERFACE

      INTERFACE ARRAYSCATTER
        MODULE PROCEDURE ARRAYSCATTER_J
        MODULE PROCEDURE ARRAYSCATTER_IJ
        MODULE PROCEDURE IARRAYSCATTER_IJ
        MODULE PROCEDURE LARRAYSCATTER_IJ
      END INTERFACE

#ifdef USE_MPI
      INTERFACE xESMF_ARRAYGATHER
        MODULE PROCEDURE ESMF_ARRAYGATHER_J_int
        MODULE PROCEDURE ESMF_IARRAYGATHER_IJ
        MODULE PROCEDURE ESMF_LARRAYGATHER_IJ
      END INTERFACE
      INTERFACE xESMF_ARRAYSCATTER
        MODULE PROCEDURE ESMF_IARRAYSCATTER_IJ
        MODULE PROCEDURE ESMF_LARRAYSCATTER_IJ
      END INTERFACE
#endif

      INTERFACE ARRAYGATHER
        MODULE PROCEDURE ARRAYGATHER_J
        MODULE PROCEDURE ARRAYGATHER_J_INT
        MODULE PROCEDURE ARRAYGATHER_IJ
        MODULE PROCEDURE IARRAYGATHER_IJ
        MODULE PROCEDURE LARRAYGATHER_IJ
      END INTERFACE

      INTERFACE DREAD_PARALLEL
        MODULE PROCEDURE DREAD_PARALLEL_2D
        MODULE PROCEDURE DREAD_PARALLEL_3D
      END INTERFACE

      INTERFACE DREAD8_PARALLEL
c        MODULE PROCEDURE DREAD_PARALLEL_2D
        MODULE PROCEDURE DREAD8_PARALLEL_3D
      END INTERFACE

      interface DWRITE8_PARALLEL
        module procedure DWRITE8_PARALLEL_3D
      end interface

      INTERFACE MREAD_PARALLEL
        MODULE PROCEDURE MREAD_PARALLEL_2D
        MODULE PROCEDURE MREAD_PARALLEL_3D
      END INTERFACE

      INTERFACE READT_PARALLEL
        MODULE PROCEDURE READT_PARALLEL_2D
        MODULE PROCEDURE READT_PARALLEL_3D
      END INTERFACE

      INTERFACE WRITET_PARALLEL
        MODULE PROCEDURE WRITET_PARALLEL_2D
c        MODULE PROCEDURE READT_PARALLEL_3D
      END INTERFACE

      INTERFACE READT_PARALLEL_COLUMN
        MODULE PROCEDURE READT_PARALLEL_COLUMN_3D
      END INTERFACE

      INTERFACE READT8_PARALLEL
c        MODULE PROCEDURE READT8_PARALLEL_2D
        MODULE PROCEDURE READT8_PARALLEL_3D
      END INTERFACE

      INTERFACE READT8_COLUMN
        MODULE PROCEDURE READT8_COLUMN_3D
      END INTERFACE

      INTERFACE WRITET8_COLUMN
        MODULE PROCEDURE WRITET8_COLUMN_3D
      END INTERFACE

      interface READ_PARALLEL
         module procedure READ_PARALLEL_INTEGER_0
         module procedure READ_PARALLEL_INTEGER_1
         module procedure READ_PARALLEL_REAL8_1
      end interface

      interface WRITE_PARALLEL
         module procedure WRITE_PARALLEL_INTEGER_0
         module procedure WRITE_PARALLEL_INTEGER_1
         module procedure WRITE_PARALLEL_REAL8_0
         module procedure WRITE_PARALLEL_REAL8_1
         module procedure WRITE_PARALLEL_STRING_0
         module procedure WRITE_PARALLEL_STRING_1

      end interface

      interface WRITEI_PARALLEL
        module procedure WRITEI_PARALLEL_2D
      end interface

      interface WRITEI8_PARALLEL
        module procedure WRITEI8_PARALLEL_3D
      end interface


      PUBLIC :: AM_I_ROOT
      interface AM_I_ROOT
         module procedure mpi_am_i_root
      end interface

!@var PACK Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
      PUBLIC :: PACK_DATA
      interface PACK_DATA
         module procedure PACK_1D       ! (i)
         module procedure IPACK_1D      ! (i)
         module procedure PACK_2D       ! (i,j)
         module procedure IPACK_2D      ! (i,j)
         module procedure LPACK_2D      ! (i,j)
         module procedure PACK_3D       ! (i,j,l)
         module procedure IPACK_3D      ! (i,j,l)
         module procedure PACK_4D       ! (i,j,l,m)
         module procedure PACK_5D       ! (i,j,l,m,n)
      end interface

      PUBLIC :: PACK_DATAj
      interface PACK_DATAj
         module procedure PACKj_2D     ! (j,k)
         module procedure IPACKj_2D    ! (j,k)
         module procedure PACKj_3D     ! (j,k,l)
         module procedure PACKj_4D     ! (j,k,l,m)
      end interface

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

!@var UNPACK Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
      PUBLIC :: UNPACK_DATA
      interface UNPACK_DATA
         module procedure UNPACK_1D      ! (i)
         module procedure UNPACK_2D      ! (i,j)
         module procedure IUNPACK_2D     ! (i,j)
         module procedure LUNPACK_2D     ! (i,j)
         module procedure UNPACK_3D      ! (i,j,l)
         module procedure IUNPACK_3D     ! (i,j,l)
         module procedure UNPACK_4D      ! (i,j,l,m)
         module procedure UNPACK_5D      ! (i,j,l,m,n)
      end interface

      PUBLIC :: UNPACK_DATAj
      interface UNPACK_DATAj
         module procedure UNPACKj_2D     ! (j,k)
         module procedure UNPACKj_3D     ! (j,k,l)
         module procedure UNPACKj_4D     ! (j,k,l,m)
      end interface

!@var PACK_COLUMN Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
      PUBLIC :: PACK_COLUMN
      interface PACK_COLUMN
         module procedure PACK_COLUMN_1D  ! (k,  j  )
         module procedure PACK_COLUMN_2D  ! (k,i,j  )
         module procedure PACK_COLUMN_i2D ! (k,i,j  )
         module procedure PACK_COLUMN_3D  ! (k,i,j,l)
      end interface

!@var UNPACK_COLUMN Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
      PUBLIC :: UNPACK_COLUMN
      interface UNPACK_COLUMN
         module procedure UNPACK_COLUMN_1D  ! (k,  j  )
         module procedure UNPACK_COLUMN_2D  ! (k,i,j  )
         module procedure IUNPACK_COLUMN_2D ! (k,i,j  )
         module procedure UNPACK_COLUMN_3D  ! (k,i,j,l)
      end interface

!@var PACK_BLOCK  Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
      PUBLIC :: PACK_BLOCK
      interface PACK_BLOCK
         module procedure IPACK_BLOCK_2D    ! (k,l,i,j  )
         module procedure  PACK_BLOCK_2D    ! (k,l,i,j  )
         module procedure  PACK_BLOCK_3D    ! (k,l,m,i,j)
      end interface

!@var UNPACK_BLOCK  Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
      PUBLIC :: UNPACK_BLOCK
      interface UNPACK_BLOCK
         module procedure IUNPACK_BLOCK_2D    ! (k,l,i,j  )
         module procedure  UNPACK_BLOCK_2D    ! (k,l,i,j  )
         module procedure  UNPACK_BLOCK_3D    ! (k,l,m,i,j)
      end interface

!@var PACK_J Generic routine to pack  a global array
!@+   with the data from the corresponding distributed array.
      PUBLIC :: PACK_J
      interface PACK_J
         module procedure PACK_J_2D
         module procedure PACK_J_3D
         module procedure PACK_J_4D
      end interface
!@var UNPACK Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
      PUBLIC :: UNPACK_J
      interface UNPACK_J
         module procedure UNPACK_J_2D
         module procedure UNPACK_J_3D
         module procedure UNPACK_J_4D
      end interface

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
        integer :: im_world
        integer :: j_strt,j_stop
        integer :: j_strt_halo,j_stop_halo
        integer :: jband_strt,jband_stop
        integer, dimension(:), pointer :: scnts,sdspl,sdspl_inplace
        integer, dimension(:), pointer :: rcnts,rdspl,rdspl_inplace
        integer, dimension(:), pointer :: j0_send,j1_send
        integer, dimension(:), pointer :: j0_recv,j1_recv
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

      ! Direction bits
      PUBLIC :: NORTH, SOUTH

!, EAST, WEST
      ! Use powers of two such that addition can be used
      ! to create combination
c***      INTEGER, PARAMETER :: NORTH = 2**0, SOUTH = 2**1
c***      INTEGER, PARAMETER :: EAST  = 2**2, WEST  = 2**3
      INTEGER, PARAMETER :: ALL = NORTH + SOUTH ! no east/west for now


      INTEGER, PARAMETER :: HALO_WIDTH = 1
      integer ::  root

#ifndef USE_DD2D_UTILS
      ! Local grid information
      TYPE DIST_GRID

         TYPE (ESMF_Grid) :: ESMF_GRID
#ifdef USE_MPP
         TYPE (domain2D ) :: domain
#endif
         ! Parameters for Global domain
         INTEGER :: IM_WORLD        ! Number of Longitudes
         INTEGER :: JM_WORLD        ! Number of latitudes
         ! Parameters for local domain
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
         ! Controls for special cases
         LOGICAL :: HAVE_SOUTH_POLE ! South pole is in local domain
         LOGICAL :: HAVE_NORTH_POLE ! North pole is in local domain
         LOGICAL :: HAVE_EQUATOR    ! Equator (JM+1)/2 is in local domain

         INTEGER, DIMENSION(:), POINTER :: DJ_MAP
         INTEGER :: DJ
         INTEGER :: log_unit ! for debugging
         !@var lookup_pet index of PET for a given J
         INTEGER, DIMENSION(:), POINTER :: lookup_pet
         LOGICAL :: BC_PERIODIC         
      END TYPE DIST_GRID
#endif

      TYPE (DIST_GRID) :: GRID, GRID_TRANS

      public :: haveLatitude

! Remaining variables are private to the module.

!@var NPES number of processes upon which work is to be distributed
      INTEGER :: NPES
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

      integer, parameter :: ROOT_ID=0
      INTEGER, PUBLIC :: CHECKSUM_UNIT

      Integer :: tag = 10

      CONTAINS


      ! This routine initializes the quantities described above.
      ! The initialization should proceed prior to any grid computations.
      subroutine init_app
#ifdef USE_ESMF
      USE ESMF_CUSTOM_MOD, Only: modelE_vm
#endif
      integer :: rc
      character(len=10) :: fileName
      NPES = 1                  ! default NPES = 1 for serial run

#ifdef USE_ESMF
      Call ESMF_Initialize(vm=modelE_vm, rc=rc)
      Call ESMF_VMGet(modelE_vm, localPET=my_pet, petCount=NPES, rc=rc)
      compmodelE  = ESMF_GridCompCreate(modelE_vm,"ModelE ESMF", rc=rc)
      ESMF_Layout = ESMF_DELayoutCreate(modelE_vm,
     &     deCountList = (/ 1, NPES /))
c along with ESMF_GridCompCreate, move this somewhere into init_grid?
c Separate atm/ocn/seaice components
c      Call ESMF_GridCompSet(compmodelE, grid=grd_dum%ESMF_GRID, rc=rc)
c$$$      write(fileName,'(a,i4.4)')'trace.',my_pet
c$$$      open(78,file=fileName,form='formatted')
#endif

#ifdef USE_MPP
c fms_init() has already been called.  Move that call here?
#ifndef USE_ESMF
      MY_PET = mpp_pe()
      NPES   = mpp_npes()
#endif
#endif

      NP_LON   = 1
      RANK_LON = 0
      NP_LAT   = NPES
      RANK_LAT = my_pet
      return
      end subroutine init_app

!      ! This routine initializes the quantities described above.
!      ! The initialization should proceed prior to any grid computations.
!      SUBROUTINE INIT_APP(grd_dum,IM,JM,LM, J_SCM)
!      USE FILEMANAGER, Only : openunit
!
!#ifdef USE_ESMF
!      USE ESMF_CUSTOM_MOD, Only: Initialize_App
!      USE ESMF_CUSTOM_MOD, Only: vm => modelE_vm
!#endif
!!AOO      USE ESMF_CUSTOM_MOD, Only: modelE_grid
!
!      IMPLICIT NONE
!      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
!      INTEGER, INTENT(IN) :: IM, JM, LM
!      INTEGER, OPTIONAL, INTENT(IN) :: J_SCM ! single column model
!      INTEGER             :: rc
!      INTEGER             :: pet
!      CHARACTER(LEN=20) :: buffer
!
!#ifndef USE_MPP
!#ifdef USE_ESMF
!      ! Initialize ESMF
!      Call Initialize_App(IM, JM, LM,rc=rc)
!
!      Call ESMF_VMGet(vm, localPET = my_pet, petCount = NPES, rc=rc)
!      root = ROOT_ID
!      compmodelE  = ESMF_GridCompCreate(vm,"ModelE ESMF", rc=rc)
!
!      ! The default layout is not what we want - it splits in the "I" direction.
!      ESMF_Layout = ESMF_DELayoutCreate(vm, deCountList = (/ 1, NPES /))
!
!      NP_LON = 1
!      NP_LAT = NPES
!      RANK_LAT = my_pet
!      RANK_LON = 0
!#else
!      MY_PET = 0
!      NPES = 1
!      NP_LON = 1
!      NP_LAT = 1
!      RANK_LON = 0
!      RANK_LAT = 0
!#endif
!
!#ifdef USE_ESMF
!      !! write(*,*) "INIT_GRID 1 USE ESMF"
!      call INIT_GRID(grd_dum,IM,JM,LM,vm=vm,CREATE_CAP=.true.)
!      Call ESMF_GridCompSet(compmodelE, grid=grd_dum%ESMF_GRID, rc=rc)
!      WRITE(*,*)'Domain Decomposition for rank: ',MY_PET,RANK_LAT,
!     &     RANK_LON
!#else
!      !! write(*,*) "INIT_GRID 2"
!      call INIT_GRID(grd_dum,IM,JM,LM, J_SCM=J_SCM)
!#endif
!#endif
!
!#ifdef USE_MPP
!      NP_LON = 1
!      NP_LAT = mpp_npes()
!      RANK_LAT = mpp_pe()
!      RANK_LON = 0
!      MY_PET=mpp_pe()
!      NPES=mpp_npes()
!      root=mpp_root_pe()
!
!#ifdef USE_ESMF
!      ! Initialize ESMF
!      Call Initialize_App(IM, JM, LM,rc=rc)
!      compmodelE  = ESMF_GridCompCreate(vm,"ModelE ESMF", rc=rc)
!
!      ESMF_Layout = ESMF_DELayoutCreate(vm, deCountList = (/ 1, NPES /))
!#endif
!      !! write(*,*) "INIT_GRID 3"
!      call INIT_GRID(grd_dum,IM,JM,LM, J_SCM=J_SCM,CREATE_CAP=.true.)
!#endif
!
!#ifdef DEBUG_DECOMP
!      IF (AM_I_ROOT()) CALL openunit('CHKSUM_DECOMP', CHECKSUM_UNIT)
!      WRITE(buffer,'(a,i3.3)') 'LOG_',my_pet
!      CALL openunit(TRIM(buffer), grd_dum%log_unit)
!#endif
!
!      END SUBROUTINE INIT_APP

      SUBROUTINE DESTROY_GRID(grd_dum)
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
#ifdef USE_ESMF
      Call ESMF_GridDestroy(grd_dum%ESMF_Grid)
#endif
      END SUBROUTINE DESTROY_GRID

#ifdef USE_ESMF
      SUBROUTINE INIT_GRID(grd_dum,IM,JM, LM,width,vm,J_SCM,bc_periodic,
     &                     CREATE_CAP)
      USE ESMF_CUSTOM_MOD, Only : modelE_vm
#else
      SUBROUTINE INIT_GRID(grd_dum,IM,JM,LM,width,J_SCM,bc_periodic,
     &                     CREATE_CAP)
#endif
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM,LM
      INTEGER, OPTIONAL, INTENT(IN) :: J_SCM ! single column model
      INTEGER, OPTIONAL :: width
      LOGICAL, OPTIONAL, INTENT(IN) :: bc_periodic
      LOGICAL, OPTIONAL, INTENT(IN) :: CREATE_CAP
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
      INTEGER :: p
#ifdef USE_MPP
      integer :: npx, npy, ng
      integer :: isd, ied , jsd, jed
      integer :: capnx, capny
#endif
      
      !! write(*,*) "BEGIN INIT_GRID"

#ifdef USE_FVCUBED
      grid_size(1)=IM;   grid_size(2)=JM*6
#else
      grid_size(1)=IM;   grid_size(2)=JM
#endif
      range_min(1)=0.;   range_min(2)=-90.
      range_max(1)=360.; range_max(2)=90.

#ifndef USE_MPP
#ifdef USE_ESMF
      grd_dum%ESMF_GRID = ESMF_GridCreateHorzLatLonUni(counts=grid_size,
     &     minGlobalCoordPerDim=range_min,
     &     maxGlobalCoordPerDim=range_max,
     &     horzStagger=ESMF_GRID_HORZ_STAGGER_A,
     &     name="source grid", rc=rc)

      Allocate(grd_dum%dj_map(0:npes-1))

      if (LM > 1) then
         deltaZ = 1.0d0
         call ESMF_GridAddVertHeight(grd_dum%ESMF_GRID,
     &         delta=(/(deltaZ, L=1,LM) /),
     &    vertStagger=ESMF_GRID_VERT_STAGGER_TOP,
     &    rc=rc)
         if (rc /= ESMF_SUCCESS)
     &        call stop_model('Failure when adding vert grid cooords.',
     &        255)
      end if

      vm_ => modelE_vm
      If (Present(vm)) vm_ => vm

      Call ESMF_VMGet(vm_, localPET = my_pet, petCount = NPES, rc=rc)
      ! The default layout is not what we want - it splits in the "I" direction.
      layout = ESMF_DELayoutCreate(vm_, deCountList = (/ 1, NPES /))

      if(npes > jm-2)
     &     call stop_model('init_grid: jm too large',255)

      if(npes == jm-2) then ! maximum for decomposition in latitude
        allocate(ims(0:0), jms(0:npes-1))
        ims(0) = im
        jms(:) = 1
        jms(0)      = 2 ! poles need 2 lats
        jms(npes-1) = 2
        Call ESMF_GridDistribute(grid=grd_dum%ESMF_GRID,
     &       countsPerDEDim1=ims, countsPerDEDim2=jms,
     &       delayout = layout, rc=rc)
        deallocate(ims, jms)
      else ! ask ESMF to make the decomposition
        Call ESMF_GridDistribute(grid=grd_dum%ESMF_GRID,
     &       delayout = layout, rc=rc)
      endif
      call ESMF_GridGet(grd_dum%esmf_grid, delayout=layout, rc=rc)
      RANK_LON=0
      RANK_LAT=my_pet
      Call ESMF_GRID_BOUNDS(grd_dum, RANK_LON, RANK_LAT,
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
#endif

      width_ = HALO_WIDTH
      If (Present(width)) width_=width

#ifdef USE_MPP
      if (present(bc_periodic)) then
         grd_dum % BC_PERIODIC = isPeriodic(bc_periodic)
      else
         grd_dum % BC_PERIODIC = .false.
      endif
#ifdef USE_FVCUBED
      npx=IM+1; npy=JM+1; ng = width_;
      ntiles=6
#else
      npx=IM+1; npy=JM+1; ng = width_;
      ntiles=1
#endif
      call init_domain(grd_dum%domain,npx,npy,ntiles,ng,
     &                 grd_dum%bc_periodic)
      RANK_LON=0
      RANK_LAT=mpp_pe()
      my_pet = mpp_pe()
      NPES = mpp_npes()
      root   = mpp_root_pe()
      call mpp_get_compute_domain( grd_dum%domain,
     &                         I0_DUM, I1_DUM, J0_DUM, J1_DUM)
      !following includes halo regions
      call mpp_get_data_domain( grd_dum%domain, isd, ied , jsd, jed)
      write(*,*)'mpp-bounds',my_pet,I0_DUM, I1_DUM, J0_DUM, J1_DUM
#endif

      grd_dum%IM_WORLD      = IM
      grd_dum%JM_WORLD      = JM

      ! Wrapped ESMF grid
      grd_dum%I_STRT        = I0_DUM
      grd_dum%I_STOP        = I1_DUM
      grd_dum%I_STRT_HALO   = MAX( 1, I0_DUM-width_)
      grd_dum%I_STOP_HALO   = MIN(IM, I1_DUM+width_)
      grd_dum%ni_loc = (RANK_LAT+1)*IM/NPES - RANK_LAT*IM/NPES

      grd_dum%J_STRT        = J0_DUM
      grd_dum%J_STOP        = J1_DUM
#ifdef USE_ESMF
      call ESMF_DElayoutBarrier(layout, rc)
#endif

cddd      IF (RANK_LAT > 0) THEN
cddd        grd_dum%J_STRT_SKP = J0_DUM
cddd      ELSE
cddd        grd_dum%J_STRT_SKP    = 2
cddd      END IF
cddd      IF (RANK_LAT < NP_LAT - 1) THEN
cddd        grd_dum%J_STOP_SKP    = J1_DUM
cddd      ELSE
cddd        grd_dum%J_STOP_SKP    = JM-1
cddd      END IF

ccc I think the following will do the same and will be compatible with SCM
      grd_dum%J_STRT_SKP = max (   2, J0_DUM)
      grd_dum%J_STOP_SKP = min (JM-1, J1_DUM)

#ifdef USE_MPI
      grd_dum%J_STRT_HALO   = J0_DUM - width_
      grd_dum%J_STOP_HALO   = J1_DUM + width_
#else
      ! I guess we don't need HALO in SCM mode...
      !grd_dum%J_STRT_HALO = MAX(1,  grd_dum % J_STRT - 1)
      !grd_dum%J_STOP_HALO = MIN(JM, grd_dum % J_STOP + 1)
      grd_dum%J_STRT_HALO = MAX(1,  grd_dum % J_STRT)
      grd_dum%J_STOP_HALO = MIN(JM, grd_dum % J_STOP)
#endif

cddd      IF (RANK_LAT > 0) THEN
cddd        grd_dum%J_STRT_STGR = J0_DUM
cddd      ELSE
cddd        grd_dum%J_STRT_STGR = 2
cddd      ENDIF
      grd_dum%J_STRT_STGR   = max(2,J0_DUM)
      grd_dum%J_STOP_STGR   = J1_DUM

      grd_dum%HAVE_SOUTH_POLE = (RANK_LAT == 0)
      grd_dum%HAVE_NORTH_POLE = (RANK_LAT == NP_LAT - 1)

      J_EQUATOR = JM/2
      grd_dum%HAVE_EQUATOR    =
     &      (J0_DUM <= J_EQUATOR) .AND. (J1_DUM >= J_EQUATOR)

#ifdef USE_DD2D_UTILS
c need to initialize the dd2d version of dist_grid for I/O
           call init_dist_grid(
     &     grd_dum%IM_WORLD,grd_dum%JM_WORLD,1, 
     &     grd_dum%I_STRT,grd_dum%I_STOP,
     &     grd_dum%J_STRT,grd_dum%J_STOP,
     &     grd_dum%I_STRT_HALO,grd_dum%I_STOP_HALO,
     &     grd_dum%J_STRT_HALO,grd_dum%J_STOP_HALO,grd_dum)
#endif

      if (present(J_SCM)) then
        ! assume J_SCM is in "general position"
        grd_dum%HAVE_SOUTH_POLE = .false.
        grd_dum%HAVE_NORTH_POLE = .false.
        grd_dum%HAVE_EQUATOR    = .false.
      endif

      grd_dum % BC_PERIODIC = isPeriodic(bc_periodic)

      ! set lookup table PET(J)
      Allocate(grd_dum%lookup_pet(1:JM))
      grd_dum%lookup_pet(:) = 0

#ifdef USE_MPI
      Allocate(AI(NPES,3))

#ifdef USE_MPP
      Allocate(grd_dum%dj_map(0:npes-1))
      Call GridGetAllAxisIndex(grd_dum%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grd_dum%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif
      Do p = 1, npes
        grd_dum%lookup_pet( AI(p,2)%min : AI(p,2)%max ) = p-1
      End Do

      do p=0, npes-1
         grd_dum%dj_map(p)=AI(p+1,2)%max - AI(p+1,2)%min + 1
      end do
      grd_dum%dj=grd_dum%dj_map(my_pet)

      Deallocate(AI)
#endif

#if (defined USE_MPP) && (defined  USE_ESMF)
# ifdef USE_FVCUBED

      if (present(CREATE_CAP)) then
         if (CREATE_CAP) then

            capnx = int(floor(sqrt(real(NPES/6))))
            capny = NPES / capnx

            cf = load_cap_config('cap.rc',IM,JM*6,LM,capnx,capny)
            vm_ => modelE_vm
            If (Present(vm)) vm_ => vm
            print*, 'Started AppGridCreateF'
            grd_dum%ESMF_GRID = AppGridCreateF(cf, vm_, rc)
            call ESMF_GridGet(grd_dum%ESMF_GRID, delayout=layout, rc=rc)
            print*, 'Finished AppGridCreateF'
         endif
      endif

# else

      grd_dum%ESMF_GRID = ESMF_GridCreateHorzLatLonUni(counts=grid_size,
     &     minGlobalCoordPerDim=range_min,
     &     maxGlobalCoordPerDim=range_max,
     &     horzStagger=ESMF_GRID_HORZ_STAGGER_A,
     &     name="source grid", rc=rc)

      if (LM > 1) then
         deltaZ = 1.0d0
         call ESMF_GridAddVertHeight(grd_dum%ESMF_GRID,
     &         delta=(/(deltaZ, L=1,LM) /),
     &    vertStagger=ESMF_GRID_VERT_STAGGER_TOP,
     &    rc=rc)
         if (rc /= ESMF_SUCCESS)
     &        call stop_model('Failure when adding vert grid cooords.',
     &        255)
      end if

      allocate( IMS(0:0) ); allocate( JMS(0:NPES-1) );
      IMS(0)=im;
      do p=0, npes-1
        JMS(p) = grd_dum%dj_map(p)
      end do
      write(*,*)'mkbhat: IMS are ',IMS
      write(*,*)'mkbhat: JMS are ',JMS

      vm_ => modelE_vm
      If (Present(vm)) vm_ => vm

      ! The default layout is not what we want - it splits in the "I" direction.
      layout = ESMF_DELayoutCreate(vm_, deCountList = (/ 1, NPES /))
      Call ESMF_GridDistribute(grid=grd_dum%ESMF_GRID,
     &     delayout = layout, 
     &     countsPerDEDim1=ims, 
     &     countsPerDEDim2=jms, rc=rc)
      call ESMF_GridValidate(grd_dum%ESMF_GRID, rc=rc)
      deallocate(ims); deallocate(jms);

      call ESMF_GridGet(grd_dum%esmf_grid, delayout=layout, rc=rc)
!     Call ESMF_GRID_BOUNDS(grd_dum, RANK_LON, RANK_LAT,
!    &        I0_DUM, I1_DUM, J0_DUM, J1_DUM)
!     write(*,*)'esmf-bounds',my_pet,I0_DUM, I1_DUM, J0_DUM, J1_DUM

# endif
#endif

      END SUBROUTINE INIT_GRID

#ifdef USE_MPP
      SUBROUTINE GridGetAllAxisIndex(domain, AI)

      implicit none
      TYPE (domain2D ) :: domain
      Type (ESMF_Axisindex) :: AI(npes,3)

      integer :: is, ie , js, je, myRank, ierr
      integer :: Axistype, oldtypes(1), blockcounts(1)
      INTEGER(KIND=MPI_ADDRESS_KIND):: offsets

      call mpp_get_compute_domain( domain, is, ie , js, je)
      myRank = mpp_pe()

      AI(myRank+1,1)%min=is; AI(myRank+1,1)%max=ie;
      AI(myRank+1,2)%min=js; AI(myRank+1,2)%max=je;

      offsets=0
      oldtypes(1)=MPI_INTEGER; blockcounts(1)=3;

      call MPI_Type_create_struct(1,blockcounts, offsets, oldtypes,
     &     Axistype, ierr)
      call MPI_TYPE_COMMIT(Axistype, ierr)

      call MPI_Allgather(AI(myRank+1,1),1, Axistype,
     &          AI(1,1),1, Axistype,MPI_COMM_WORLD,ierr )
      call MPI_Allgather(AI(myRank+1,2),1, Axistype,
     &          AI(1,2),1, Axistype,MPI_COMM_WORLD,ierr )

      END SUBROUTINE GridGetAllAxisIndex
#endif

      SUBROUTINE GET(grd_dum, I_STRT, I_STOP,
     &                        I_STRT_HALO, I_STOP_HALO,
     &                        J_STRT, J_STOP, J_STRT_HALO, J_STOP_HALO,
     &                        J_STRT_SKP, J_STOP_SKP,
     &                        J_STRT_STGR, J_STOP_STGR,
     &                        HAVE_SOUTH_POLE, HAVE_NORTH_POLE)
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER, OPTIONAL :: I_STRT, I_STOP
      INTEGER, OPTIONAL :: I_STRT_HALO, I_STOP_HALO
      INTEGER, OPTIONAL :: J_STRT, J_STOP
      INTEGER, OPTIONAL :: J_STRT_HALO, J_STOP_HALO
      INTEGER, OPTIONAL :: J_STRT_SKP, J_STOP_SKP
      INTEGER, OPTIONAL :: J_STRT_STGR, J_STOP_STGR
      LOGICAL, OPTIONAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      IF (PRESENT(I_STRT)) I_STRT = grd_dum%I_STRT
      IF (PRESENT(I_STOP)) I_STOP = grd_dum%I_STOP

      IF (PRESENT(I_STRT_HALO)) I_STRT_HALO = grd_dum%I_STRT_HALO
      IF (PRESENT(I_STOP_HALO)) I_STOP_HALO = grd_dum%I_STOP_HALO

      IF (PRESENT(J_STRT)) J_STRT = grd_dum%J_STRT
      IF (PRESENT(J_STOP)) J_STOP = grd_dum%J_STOP

      IF (PRESENT(J_STRT_HALO)) J_STRT_HALO = grd_dum%J_STRT_HALO
      IF (PRESENT(J_STOP_HALO)) J_STOP_HALO = grd_dum%J_STOP_HALO

      IF (PRESENT(J_STRT_SKP)) J_STRT_SKP = grd_dum%J_STRT_SKP
      IF (PRESENT(J_STOP_SKP)) J_STOP_SKP = grd_dum%J_STOP_SKP

      IF (PRESENT(J_STRT_STGR)) J_STRT_STGR = grd_dum%J_STRT_STGR
      IF (PRESENT(J_STOP_STGR)) J_STOP_STGR = grd_dum%J_STOP_STGR

      IF (PRESENT(HAVE_SOUTH_POLE))
     &             HAVE_SOUTH_POLE= grd_dum%HAVE_SOUTH_POLE
      IF (PRESENT(HAVE_NORTH_POLE))
     &             HAVE_NORTH_POLE= grd_dum%HAVE_NORTH_POLE

      END SUBROUTINE GET

      SUBROUTINE HALO_UPDATE_1D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                        arr(grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 1, from
     &     ,grd_dum%BC_PERIODIC)
#endif

      END SUBROUTINE HALO_UPDATE_1D

      SUBROUTINE HALO_UPDATE_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                    arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 2, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATE_2D

      SUBROUTINE HALO_UPDATE_2Dint(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      INTEGER,            INTENT(INOUT) ::
     &                    arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_MPI
      Call sendrecv_int(grd_dum%ESMF_GRID, arr, shape(arr), 2, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATE_2Dint

      SUBROUTINE HALO_UPDATEj_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                    arr(grd_dum%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 1, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATEj_2D

      SUBROUTINE HALO_UPDATE_3D(grd_dum, arr, from, jdim)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: arr(:,:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from
      INTEGER, OPTIONAL, INTENT(IN)    :: jdim

      INTEGER :: jd

#ifdef USE_MPI
      if(present(jdim)) then
        jd = jdim
      else
        jd = 2
      endif
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), jd, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATE_3D

      SUBROUTINE HALO_UPDATE_mask(grd_dum, sBufS, sBufN, rBufS, rBufN)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(IN) :: sBufN(:,:)
      REAL*8,            INTENT(IN) :: sBufS(:,:)
      REAL*8,            INTENT(OUT) :: rBufN(:,:)
      REAL*8,            INTENT(OUT) :: rBufS(:,:)

#ifdef USE_MPI
      integer :: numSendSouth, numSendNorth
      integer :: numRecvSouth, numRecvNorth
      integer :: pe_south, pe_north
      integer, save :: tag = 1
      integer, parameter :: NUM_TAGS = 100
      integer :: status(MPI_STATUS_SIZE)
      integer :: ier
#endif

#ifdef USE_MPI
      numSendSouth = size(sBufS)
      numSendNorth = size(sBufN)

      numRecvSouth = size(rBufS)
      numRecvNorth = size(rBufN)

      call getNeighbors(my_pet, npes, pe_south, pe_north, .false.)

      tag = 1 + mod(tag - 1, NUM_TAGS)
      call MPI_sendrecv(sBufS, numSendSouth, MPI_DOUBLE_PRECISION,
     &     pe_south, tag,
     &     rBufN, numRecvNorth, MPI_DOUBLE_PRECISION, pe_north, tag,
     &     MPI_COMM_WORLD, status, ier)

      tag = 1 + mod(tag - 1, NUM_TAGS)

      call MPI_sendrecv(sBufN, numSendNorth, MPI_DOUBLE_PRECISION,
     &     pe_north, tag,
     &     rBufS, numRecvSouth, MPI_DOUBLE_PRECISION, pe_south, tag,
     &     MPI_COMM_WORLD, status, ier)

#endif
      END SUBROUTINE HALO_UPDATE_mask

      SUBROUTINE HALO_UPDATE_COLUMN_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                  arr(:,grd_dum%j_strt_halo:)

      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 2, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_2D

      SUBROUTINE HALO_UPDATE_COLUMN_3D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                  arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 3, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_3D

      SUBROUTINE INT_HALO_UPDATE_COLUMN_3D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      INTEGER,           INTENT(INOUT) ::
     &                  arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      REAL*8 :: foo(size(arr,1),size(arr,2),size(arr,3))

      ! cheat
      foo = arr
      CALL HALO_UPDATE_COLUMN(grd_dum, foo, from)
      arr = nint(foo)

      END SUBROUTINE INT_HALO_UPDATE_COLUMN_3D


      SUBROUTINE HALO_UPDATE_COLUMN_4D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &               arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 3, from)
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_4D

      SUBROUTINE HALO_UPDATE_BLOCK_4D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &               arr(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 4, from)
#endif
      END SUBROUTINE HALO_UPDATE_BLOCK_4D

      SUBROUTINE HALO_UPDATE_COLUMN_7D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &     arr(:,:,:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_MPI
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 6, from
     &     ,grd_dum%BC_PERIODIC)
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_7D

      SUBROUTINE CHECKSUM_1D(grd_dum, arr, line, file, unit, STGR, SKIP)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr
      LOGICAL, OPTIONAL, INTENT(IN) :: skip


      INTEGER :: unit_
      REAL*8  :: asum, L1norm
      REAL*8  :: t_arr(grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      INTEGER :: J_0, J_1
      INTEGER :: stgr_, skip_

#ifdef DEBUG_DECOMP
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      skip_ = 0
      If (Present(skip)) THEN
         If (skip) skip_=1
      End If

      t_arr = arr
      Call GLOBALSUM(grd_dum, t_arr, asum,istag=stgr_,iskip=skip_)
      t_arr(J_0:J_1) = ABS(t_arr(J_0:J_1))
      Call GLOBALSUM(grd_dum, t_arr, L1norm,istag=stgr_,iskip=skip_)

      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &     file,line, asum, L1norm

#endif

      END SUBROUTINE CHECKSUM_1D

      SUBROUTINE CHECKSUM_2D(grd_dum, arr, line, file, unit, stgr,
     &     skip)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr
      LOGICAL, OPTIONAL, INTENT(IN) :: skip

      INTEGER :: unit_
      REAL*8  :: asum, L1norm
      REAL*8  :: asum_glob, L1norm_glob
      REAL*8 ::
     &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      INTEGER :: J_0, J_1, I,J
      INTEGER :: stgr_,skip_

      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      skip_ = 0
      If (Present(skip)) THEN
         If (skip) skip_=1
      End If

      t_arr(:,J_0:J_1) = arr(:,J_0:J_1)
      Call GLOBALSUM(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call GLOBALSUM(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)

      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &     file,line, asum, L1norm


      END SUBROUTINE CHECKSUM_2D

      SUBROUTINE CHECKSUM_3D(grd_dum, arr, line, file, unit, stgr, skip)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr
      LOGICAL, OPTIONAL, INTENT(IN) :: skip


      INTEGER :: unit_
      INTEGER :: k
      REAL*8, DIMENSION(Size(arr,3))  :: asum, L1norm

      REAL*8 ::
     &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      INTEGER :: J_0, J_1
      Integer :: stgr_,skip_


      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      skip_ = 0
      If (Present(skip)) THEN
         If (skip) skip_=1
      End If

      Do k = 1, Size(arr, 3)
        t_arr(:,J_0:J_1) = arr(:,J_0:J_1,k)
        Call GLOBALSUM(grd_dum, t_arr, asum(k), istag=stgr_,iskip=skip_)
        t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
        Call GLOBALSUM(grd_dum, t_arr,L1norm(k),istag=stgr_,iskip=skip_)
      End Do
      If (AM_I_ROOT()) Then
        Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &       file,line, Sum(asum), Sum(L1norm)
      End If

      END SUBROUTINE CHECKSUM_3D

      SUBROUTINE CHECKSUMj_2D(grd_dum, arr, line, file, unit, stgr,
     &     skip)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(grd_dum%j_strt_halo:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr
      LOGICAL, OPTIONAL, INTENT(IN) :: skip

      INTEGER :: unit_
      REAL*8  :: asum, L1norm
      REAL*8  :: asum_glob, L1norm_glob
      REAL*8 ::
     &  t_arr(size(arr,2),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      INTEGER :: J_0, J_1
      INTEGER :: stgr_,skip_

      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      skip_ = 0
      If (Present(skip)) THEN
         If (skip) skip_=1
      End If

      t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:))
      Call GLOBALSUM(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call GLOBALSUM(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)

      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &     file,line, asum, L1norm


      END SUBROUTINE CHECKSUMj_2D

      SUBROUTINE CHECKSUMj_3D(grd_dum, arr, line, file, unit, stgr)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(grd_dum%j_strt_halo:,:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr


      INTEGER :: unit_
      INTEGER :: k
      REAL*8, DIMENSION(Size(arr,3))  :: asum, L1norm

      REAL*8 ::
     &  t_arr(size(arr,2),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      INTEGER :: J_0, J_1
      Integer :: stgr_


      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      Do k = 1, Size(arr, 3)
         t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:,k))
         Call GLOBALSUM(grd_dum, t_arr, asum(k), istag=stgr_)
         t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
         Call GLOBALSUM(grd_dum, t_arr, L1norm(k), istag=stgr_)
      End Do
      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &     file,line, Sum(asum), Sum(L1norm)


      END SUBROUTINE CHECKSUMj_3D


      SUBROUTINE CHECKSUM_COLUMN_2D(grd_dum,arr,line,file,unit,stgr)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(:,grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr


      INTEGER :: unit_
      INTEGER :: k
      REAL*8, DIMENSION(SIZE(arr,1))  :: asum, L1norm
      INTEGER :: stgr_

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit
      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      Do k = 1, Size(arr, 1)
         Call GLOBALSUM(grd_dum,      arr(k,:), asum(k), istag=stgr_)
         Call GLOBALSUM(grd_dum, Abs(arr(k,:)), L1norm(k), istag=stgr_)
      END DO

      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &     file,line, Sum(asum), Sum(L1norm)


      END SUBROUTINE CHECKSUM_COLUMN_2D

      SUBROUTINE CHECKSUM_COLUMN_3D(grd_dum,arr,line,file,unit,stgr)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr


      INTEGER :: unit_
      INTEGER :: k
      REAL*8, DIMENSION(Size(arr,1))  :: asum, L1norm
      INTEGER :: stgr_

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit
      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      Do k = 1, Size(arr, 1)
        Call GLOBALSUM(grd_dum,      arr(k,:,:), asum(k), istag=stgr_)
        Call GLOBALSUM(grd_dum, Abs(arr(k,:,:)), L1norm(k), istag=stgr_)
      END DO

      If (AM_I_ROOT()) Then
        Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &       file,line, Sum(asum), Sum(L1norm)
        CALL SYS_FLUSH(unit_)
      End If

      END SUBROUTINE CHECKSUM_COLUMN_3D

      SUBROUTINE INT_CHECKSUM_COLUMN_3D(grd_dum, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      INTEGER,           INTENT(IN) ::
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: k
      REAL*8, DIMENSION(Size(arr,1))  :: asum, L1norm

      Call CHECKSUM_COLUMN(grd_dum, Real(arr,KIND=KIND(1.0D+0)), line,
     &     file, unit)

      END SUBROUTINE INT_CHECKSUM_COLUMN_3D

      SUBROUTINE CHECKSUM_COLUMN_4D(grd_dum,arr,line,file,unit,stgr,
     &                              skip)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) ::
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr
      LOGICAL, OPTIONAL, INTENT(IN) :: skip


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1, k
      REAL*8, DIMENSION(Size(arr,1),Size(arr,4)) :: asum, L1norm
      Integer :: stgr_,skip_

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grd_dum%i_strt
      i_1 = grd_dum%i_stop
      j_0 = grd_dum%j_strt
      j_1 = grd_dum%j_stop

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      skip_ = 0
      If (Present(skip)) THEN
         If (skip) skip_=1
      End If

      Do k = 1, Size(arr,1)
        CALL GLOBALSUM(grd_dum,    arr(k,:,:,:),   asum(k,:),
     &        istag=stgr_, iskip=skip_)
        CALL GLOBALSUM(grd_dum,Abs(arr(k,:,:,:)),L1norm(k,:),
     &       istag=stgr_, iskip=skip_)
      End Do

      If (AM_I_ROOT()) Then
        Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))')
     &       file,line, SUM(asum), SUM(L1norm)
      End If

      END SUBROUTINE CHECKSUM_COLUMN_4D


      SUBROUTINE FINISH_APP()
      USE FILEMANAGER, ONLY : closeunit
      IMPLICIT NONE

      INTEGER :: rc

#ifdef DEBUG_DECOMP
      CALL closeunit(CHECKSUM_UNIT)
      CALL closeunit(grid%log_unit)
#endif

#ifdef USE_ESMF
      CALL ESMF_FINALIZE(rc=rc)
#endif

      END SUBROUTINE FINISH_APP

      SUBROUTINE SUMXPE_1D(arr, arr_master, increment)
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
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION,
     &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + arr_tmp
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION,
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
c**** arr plays both roles of local and global array
c**** arr is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,root,
     &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D

      SUBROUTINE SUMXPE_1D_I(arr, arr_master, increment)
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
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_INTEGER,
     &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + arr_tmp
          deallocate(arr_tmp)
        endif
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_INTEGER,
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
c**** arr plays both roles of local and global array
c**** arr is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size,
     &        MPI_INTEGER,MPI_SUM,root,
     &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_1D_I

      SUBROUTINE SUMXPE_2D(arr, arr_master, increment)
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
            call mpi_reduce(arr,arr_tmp,arr_size,
     &           MPI_DOUBLE_PRECISION,MPI_SUM,root,
     &           MPI_COMM_WORLD, ierr)
            if(am_i_root()) then
              arr_master = arr_master + reshape(arr_tmp,shape(arr))
            endif
            deallocate(arr_tmp)
         else
            call mpi_reduce(arr,arr_master,arr_size,
     &           MPI_DOUBLE_PRECISION,MPI_SUM,root,
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
c**** arr plays both roles of local and global array
c**** arr is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,root,
     &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_2D

      SUBROUTINE SUMXPE_3D(arr, arr_master, increment)
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
        call mpi_reduce(arr,arr_tmp,arr_size,MPI_DOUBLE_PRECISION,
     &       MPI_SUM,root,MPI_COMM_WORLD, ierr)
        if(am_i_root()) then
          arr_master = arr_master + reshape(arr_tmp,shape(arr))
        endif
        deallocate(arr_tmp)
      else
        call mpi_reduce(arr,arr_master,arr_size,MPI_DOUBLE_PRECISION,
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
c**** arr plays both roles of local and global array
c**** arr  is overwritten by itself after reduction
#ifdef USE_ESMF
         arr_size = size(arr)
         allocate(arr_tmp(arr_size))
         call mpi_reduce(arr,arr_tmp,arr_size,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,root,
     &        MPI_COMM_WORLD, ierr)
         arr=reshape(arr_tmp,shape(arr))
         deallocate(arr_tmp)
#endif
      endif
      END SUBROUTINE SUMXPE_3D



!---------------------------
      SUBROUTINE GLOBALSUM_INT_REDUCE(grd_dum, ivar, isum, all)
      IMPLICIT NONE
      TYPE(DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER, INTENT(IN)  :: ivar
      INTEGER, INTENT(OUT) :: isum
      LOGICAL, OPTIONAL, INTENT(IN) :: all

      LOGICAL :: all_
      INTEGER :: ierr

      all_ = .false.
      If (Present(all)) all_ = all

#ifdef USE_MPI
      If (all_) Then
         call MPI_Allreduce(ivar, isum, 1, MPI_INTEGER, MPI_SUM,
     &        MPI_COMM_WORLD, ierr)
      Else
         call MPI_Reduce(ivar, isum, 1, MPI_INTEGER, MPI_SUM, root,
     &        MPI_COMM_WORLD, ierr)
      End If
#else
      isum = ivar
#endif

      END SUBROUTINE GLOBALSUM_INT_REDUCE

      SUBROUTINE GLOBALSUM_J(grd_dum, arr, gsum,
     &                       hsum, istag, iskip, all, jband)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: arr(grd_dum%j_strt_halo:)
      REAL*8,            INTENT(OUT):: gsum
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2)
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      INTEGER,OPTIONAL,  INTENT(IN) :: iskip
      LOGICAL,OPTIONAL,  INTENT(IN) :: all
      INTEGER, OPTIONAL, INTENT(IN) :: jband(2)

      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, J, ierr
      REAL*8  :: garr(grd_dum%jm_world)
      LOGICAL :: istag_, iskip_

      INTEGER :: JSTRT, JSTOP
    ! now local

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      iskip_ = .false.
      If (Present(iskip)) Then
        If (iskip == 1) iskip_ = .true.
      End If

      JSTRT=1
      JSTOP=JM

      IF (PRESENT(jband)) THEN
        JSTRT=jband(1)
        JSTOP=jband(2)
      END IF

#ifdef DEBUG_DECOMP
      If (Size(arr) /= grd_dum%j_stop_halo-grd_dum%j_strt_halo+1) Then
#ifdef USE_MPI
        CALL MPI_ABORT(ierr)
#else
        STOP
#endif
      End If
#endif

#ifdef USE_MPI
      Call gather(grd_dum, arr, garr, shape(arr), 1)
#else
      garr = 0.d0
      garr(j_0:j_1) = arr(j_0:j_1)
#endif

      If (AM_I_ROOT()) then
         If (istag_) then
           gsum = sum(garr(2:JM),1)
         ElseIf (iskip_) then
           gsum = sum(garr(2:JM-1),1)
         Else
           gsum = sum(garr(JSTRT:JSTOP),1)
         EndIf
         If (Present(hsum)) then
            If (istag_) then
               hsum(1)   = Sum( garr(2     :JM/2),1   )
               hsum(2)   = Sum( garr(2+JM/2:JM  ),1   )
               hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
               hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
            Else
               hsum(1)   = Sum( garr(1     :JM/2),1   )
               hsum(2)   = Sum( garr(1+JM/2:JM  ),1   )
            EndIf
         EndIf
      EndIf

#ifdef USE_MPI
      If (Present(all)) Then
         If (all) THEN
            Call MPI_BCAST(gsum,1,MPI_DOUBLE_PRECISION,root,
     &           MPI_COMM_WORLD, ierr)
            If (Present(hsum)) Call MPI_BCAST(hsum,2,
     &           MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, ierr)
         End If
      End If
#endif
      END SUBROUTINE GLOBALSUM_J

      SUBROUTINE GLOBALSUM_IJ(grd_dum, arr, gsum, hsum, zsum, istag,
     &                        all, iskip)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: arr(grd_dum%i_strt_halo:,
     &                                     grd_dum%j_strt_halo:)
      REAL*8,            INTENT(OUT):: gsum
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2)
      REAL*8, OPTIONAL,  INTENT(OUT):: zsum(grd_dum%j_strt:
     &                                      grd_dum%j_stop)
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      INTEGER,OPTIONAL,  INTENT(IN) :: iskip

      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, J,J_0STG, ierr, J_0S,J_1S
      REAL*8  :: zon(grd_dum%j_strt:grd_dum%j_stop)
      REAL*8  :: garr(grd_dum%jm_world)
      LOGICAL :: istag_,iskip_
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

    ! now local
#ifdef USE_ESMF
!     type (ESMF_Grid)                           :: GRID
#endif

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_0STG  = grd_dum%j_strt_stgr
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      j_0S = grd_dum%j_strt_skp
      j_1S  = grd_dum%j_stop_skp

      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      iskip_ = .false.
      If (Present(iskip)) Then
        If (iskip == 1) iskip_ = .true.
      End If

      If (istag_) Then
        zon(J_0STG:J_1)  = sum(arr(:,j_0STG:j_1),1)
      Else If (iskip_) Then
        zon(J_0S:J_1S)  = sum(arr(:,j_0S:j_1S),1)
      Else
        zon  = sum(arr(:,j_0:j_1),1)
      End If

#ifdef USE_MPI
      Call gather(grd_dum, zon, garr, shape(zon), 1)
#else
      garr = 0.d0
      garr(j_0:j_1) = zon(j_0:j_1)
#endif
      If (AM_I_ROOT()) then
        If (istag_) then
          gsum = sum(garr(2:JM),1)
        Else if(iskip_) then
          gsum = sum(garr(2:JM-1),1)
        Else
          gsum = sum(garr(1:JM),1)
        EndIf
        If (Present(hsum)) then
          If (istag_) then
            hsum(1)   = Sum( garr(2     :JM/2),1   )
            hsum(2)   = Sum( garr(2+JM/2:JM  ),1   )
            hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
            hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
          Else
            hsum(1)   = Sum( garr(1     :JM/2),1   )
            hsum(2)   = Sum( garr(1+JM/2:JM  ),1   )
          EndIf
        EndIf
      EndIf
      If (Present(zsum)) zsum = zon

#ifdef USE_MPI
      If (Present(all)) Then
         If (all) Then
            Call MPI_BCAST(gsum,1,MPI_DOUBLE_PRECISION,root,
     &           MPI_COMM_WORLD, ierr)
            If (Present(hsum))     Call MPI_BCAST(hsum,2,
     &           MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, ierr)
         End If
      End If
#endif

      END SUBROUTINE GLOBALSUM_IJ

      SUBROUTINE GLOBALSUM_IJK(grd_dum, arr, gsum, hsum, zsum, istag,
     &                         iskip)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,             INTENT(IN) :: arr(grd_dum%i_strt_halo:,
     &                                      grd_dum%j_strt_halo:,:)
      REAL*8,             INTENT(OUT):: gsum(size(arr,3))
      REAL*8, OPTIONAL,   INTENT(OUT):: hsum(2,size(arr,3))
      REAL*8, OPTIONAL,   INTENT(OUT):: zsum(grd_dum%j_strt:
     &                                       grd_dum%j_stop,
     &                                       size(arr,3))
      INTEGER,OPTIONAL,   INTENT(IN) :: istag,iskip

      INTEGER :: k
      INTEGER :: i_0, i_1, j_0, j_1, j_0S, j_1S, IM, JM, j_0STG
      REAL*8  :: zon(grd_dum%j_strt:grd_dum%j_stop,size(arr,3))
      REAL*8  :: garr(grd_dum%jm_world,size(arr,3))
      LOGICAL :: istag_,iskip_
    ! now local
#ifdef USE_ESMF
!     type (ESMF_Grid)                           :: GRID
#endif

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_0S = grd_dum%j_strt_skp
      j_0STG  = grd_dum%j_strt_stgr
      j_1  = grd_dum%j_stop
      j_1S = grd_dum%j_stop_skp
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      iskip_ = .false.
      If (Present(iskip)) Then
        If (iskip == 1) iskip_ = .true.
      End If

      If (istag_) Then
         zon(j_0STG:j_1,:)  = sum(arr(:,j_0STG:j_1,:),1)
      Else if (iskip_) Then
         zon(j_0S:j_1S,:)  = sum(arr(:,j_0S:j_1S,:),1)
      Else
         zon  = sum(arr(:,j_0:j_1,:),1)
      End If

#ifdef USE_MPI
      Call gather(grd_dum, zon, garr, shape(zon), 1)
#else
      garr(:,:) = 0.d0
      garr(j_0:j_1,:)  = zon(j_0:j_1,:)
#endif
      If (AM_I_ROOT()) then
         If (istag_) then
           gsum = sum(garr(2:JM,:),1)
         Else if(iskip_) then
           gsum = sum(garr(2:JM-1,:),1)
         Else
           gsum = sum(garr(1:JM,:),1)
         EndIf
         If (Present(hsum)) then
            If (istag_) then
               hsum(1,:)   = Sum( garr(2     :JM/2,:),1   )
               hsum(2,:)   = Sum( garr(2+JM/2:JM  ,:),1   )
               hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
               hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
            Else
               hsum(1,:)   = Sum( garr(1     :JM/2,:),1   )
               hsum(2,:)   = Sum( garr(1+JM/2:JM  ,:),1   )
            EndIf
         EndIf
      EndIf
      If (Present(zsum)) zsum = zon

      END SUBROUTINE GLOBALSUM_IJK

      SUBROUTINE GLOBALSUM_OTHER_IJK(grd_dum, arr, gsum, jband, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,             INTENT(IN) :: arr(:,grd_dum%j_strt_halo:)
      REAL*8,             INTENT(OUT):: gsum(size(arr,1))
      INTEGER,           INTENT(IN) :: jband(2)
      Logical,OPTIONAL,   INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
      Logical :: all_
      INTEGER :: ierr
#ifdef USE_MPI
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world)
#endif
    ! now local
#ifdef USE_ESMF
!     type (ESMF_Grid)                           :: GRID
#endif

      all_ = .false.
      If (Present(all)) all_ = all

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM = SIZE(arr,1)
      JM   = grd_dum%JM_WORLD

      jb1 = jband(1)
      jb2 = jband(2)

#ifdef USE_MPI
      Call gather(grd_dum, arr, garr, shape(arr), 2)
      IF (AM_I_ROOT()) gsum = Sum(garr(:,jb1:jb2),2)
      If (all_) Then
         call MPI_BCAST(gsum, Size(gsum), MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, ierr)
      End If
#else
      gsum = Sum(arr(:, max(jb1,j_0):min(jb2,j_1) ),2)
#endif
      END SUBROUTINE GLOBALSUM_OTHER_IJK

      SUBROUTINE GLOBALSUM_OTHER_IJK_IK(grd_dum, arr, gsum, jband, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,             INTENT(IN) :: arr(:,grd_dum%j_strt_halo:,:)
      REAL*8,             INTENT(OUT):: gsum(size(arr,1), size(arr,3))
      INTEGER,            INTENT(IN) :: jband(2)
      Logical,OPTIONAL,   INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
      Logical :: all_
      INTEGER :: ierr
#ifdef USE_MPI
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world,size(arr,3))
#endif
    ! now local
#ifdef USE_ESMF
!     type (ESMF_Grid)                           :: GRID
#endif

      all_ = .false.
      If (Present(all)) all_ = all

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM = SIZE(arr,1)
      JM   = grd_dum%JM_WORLD

      jb1 = jband(1)
      jb2 = jband(2)

#ifdef USE_MPI
      Call gather(grd_dum, arr, garr, shape(arr), 2)
      IF (AM_I_ROOT()) gsum = Sum(garr(:,jb1:jb2,:),2)
      If (all_) Then
         call MPI_BCAST(gsum, Size(gsum), MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, ierr)
      End If
#else
      gsum = Sum(arr(:, max(jb1,j_0):min(jb2,j_1) ,:),2)
#endif
      END SUBROUTINE GLOBALSUM_OTHER_IJK_IK

      SUBROUTINE GLOBALSUM_IJK_IK(grd_dum, arr, gsum, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,             INTENT(IN) :: arr(:,grd_dum%j_strt_halo:,:)
      REAL*8,             INTENT(OUT):: gsum(size(arr,1), size(arr,3))
      Logical,OPTIONAL,   INTENT(IN) :: all

      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, LM
      Logical :: all_
      INTEGER :: ierr
#ifdef USE_MPI
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world,size(arr,3))
#endif
    ! now local
#ifdef USE_MPI
!     type (ESMF_Grid)                           :: GRID
      Integer :: scnts(0:npes-1), sdspl(0:npes-1)
      Integer :: rcnts(0:npes-1), rdspl(0:npes-1)
      Integer ::  dik_map(0:npes-1), dik, dik_sum
      Integer :: nik, i,k,j,p, ik, ijk, iremain
      Real*8, Allocatable :: tsum(:)
      Real*8, Allocatable :: send_buf(:)
      Real*8, Allocatable :: recv_buf(:,:)
#endif

      all_ = .false.
      If (Present(all)) all_ = all

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM = SIZE(arr,1)
      JM   = grd_dum%JM_WORLD
      LM =  SIZE(arr,3)

#ifdef USE_MPI
! Number of sums each processor computes is dik
      Do p = 0, npes-1
         dik_map(p)=(IM*LM)/npes
         iremain=mod(IM*LM, npes)
         if (iremain > 0 .and. iremain > p)  dik_map(p)= dik_map(p)+1
      end do

      dik=dik_map(my_pet)

      Allocate(tsum(dik))
      Allocate(send_buf(maxval(dik_map) *(grd_dum%dj)*npes))
      Allocate(recv_buf(dik, JM))

! ugly packing for transpose
      ijk = 0
      Do j = 1, (grd_dum%dj)
         p = 0
         dik_sum=0
         nik = 0
         ijk = dik_map(p)*(j-1)
         Do k = 1, lm
            do i = 1, im
               nik = nik+1
               ijk=ijk+1
               send_buf(ijk+dik_sum) = arr(i,j+grd_dum%j_strt-1,k)
               If (nik == dik_map(p)) Then
                  dik_sum=dik_sum+dik_map(p)*(grd_dum%dj)
                  p = p+1
                  if (p == npes) exit
                  nik= 0
                  ijk = dik_map(p)*(j-1)
               End If
            end do
         end do
      end do

      scnts=dik_map*(grd_dum%dj)
      sdspl(0)=0
      Do p = 1, npes-1
         sdspl(p)=sdspl(p-1)+(grd_dum%dj)*dik_map(p-1)
      End Do

      rcnts=dik*(grd_dum%dj_map)
      rdspl(0)=0
      Do p = 1, npes-1
         rdspl(p)=rdspl(p-1)+rcnts(p-1)
      End Do

      Call MPI_ALLTOALLV(send_buf, scnts, sdspl, mpi_double_precision,
     &             recv_buf, rcnts, rdspl, mpi_double_precision,
     &             mpi_comm_world, ierr)

      tsum=sum(recv_buf,2)

      rcnts=dik_map
      rdspl(0)=0
      Do p = 1, npes-1
         rdspl(p)=rdspl(p-1)+rcnts(p-1)
      End Do

      Call MPI_GatherV(tsum, dik, mpi_double_precision,
     & gsum, dik_map, rdspl, mpi_double_precision,
     & root, mpi_comm_world, ierr)

      Deallocate(recv_buf)
      Deallocate(send_buf)
      Deallocate(tsum)

      if (all_) Then
         call MPI_BCAST(gsum, Size(gsum), MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, ierr)
      End If
#else
      gsum = Sum(arr(:,j_0:j_1,:),2)
#endif
      END SUBROUTINE GLOBALSUM_IJK_IK

      SUBROUTINE GLOBALSUM_JK(grd_dum, arr, gsum, hsum, istag, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: arr(grd_dum%j_strt_halo:,:)
      REAL*8,            INTENT(OUT):: gsum(size(arr,2))
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2,size(arr,2))
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: ierr
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM
      REAL*8  :: garr(grd_dum%jm_world,size(arr,2))
      LOGICAL :: istag_

    ! now local

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If


#ifdef USE_MPI
      Call gather(grd_dum, arr, garr, shape(arr), 1)
#else
      garr(:,:) = 0.d0
      garr(j_0:j_1,:) = arr(j_0:j_1,:)
#endif
      If (AM_I_ROOT()) then
         If (istag_) then
           gsum = sum(garr(2:JM,:),1)
         Else
           gsum = sum(garr(1:JM,:),1)
         EndIf
         If (Present(hsum)) then
            If (istag_) then
               hsum(1,:)   = Sum( garr(2     :JM/2,:),1   )
               hsum(2,:)   = Sum( garr(2+JM/2:JM  ,:),1   )
               hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
               hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
            Else
               hsum(1,:)   = Sum( garr(1     :JM/2,:),1   )
               hsum(2,:)   = Sum( garr(1+JM/2:JM  ,:),1   )
            EndIf
         EndIf
      EndIf

#ifdef USE_MPI
      If (Present(all)) Then
         If (all) Then
            Call MPI_BCAST(gsum,Size(gsum),MPI_DOUBLE_PRECISION,root,
     &           MPI_COMM_WORLD, ierr)
            If (Present(hsum)) Call MPI_BCAST(hsum,size(hsum),
     &           MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
         End If
      End If
#endif
      END SUBROUTINE GLOBALSUM_JK

      SUBROUTINE GLOBALSUM_XXXJ_XXX(grd_dum, arr, gsum, all)
      Type (DIST_GRID), INTENT(IN) :: grd_dum
      Real*8, INTENT(IN) :: arr(:,:,:,grd_dum%j_strt_halo:)
      Real*8, INTENT(Out) :: gsum(:,:,:)
      Logical, Optional, INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM
      REAL*8  :: garr(size(arr,1),size(arr,2),size(arr,3),
     &     grd_dum%jm_world)
      LOGICAL :: istag_,all_

    ! now local
#ifdef USE_ESMF
!     type (ESMF_Grid)                           :: GRID
#endif

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

#ifdef USE_MPI
      Call gather(grd_dum, arr, garr, shape(arr), 4, all=all)
      all_=.false.
      If (Present(all)) all_=all
      If (AM_I_ROOT() .or. all_)  gsum = SUM(garr,4)
#else
      gsum = SUM(arr(:,:,:,J_0:J_1),4)
#endif

      END SUBROUTINE GLOBALSUM_XXXJ_XXX

      SUBROUTINE GLOBALSUM_XXXIJ_XXX(grd_dum, arr, gsum, all)
      Type (DIST_GRID), INTENT(IN) :: grd_dum
      Real*8, INTENT(IN) ::
     &     arr(:,:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      Real*8, INTENT(Out) :: gsum(:,:,:)
      Logical, Optional, INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM
      REAL*8  :: garr(size(arr,1),size(arr,2),size(arr,3),
     &     grd_dum%jm_world)
      REAL*8  :: larr(size(arr,1),size(arr,2),size(arr,3),
     &     grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      LOGICAL :: istag_,all_

    ! now local
#ifdef USE_ESMF
!     type (ESMF_Grid)                           :: GRID
#endif

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

#ifdef USE_MPI
      larr = sum(arr,4)
      Call gather(grd_dum, larr, garr, shape(larr), 4, all=all)
      all_=.false.
      If (Present(all)) all_=all
      If (AM_I_ROOT() .or. all_)  gsum = SUM(garr,4)
#else
      gsum = SUM(SUM(arr(:,:,:,I_0:I_1,J_0:J_1),4),4)
#endif

      END SUBROUTINE GLOBALSUM_XXXIJ_XXX

      SUBROUTINE BACKSPACE_PARALLEL(IUNIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT

      If (AM_I_ROOT()) backspace IUNIT

      END SUBROUTINE BACKSPACE_PARALLEL

      SUBROUTINE REWIND_PARALLEL(IUNIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT

      If (AM_I_ROOT()) rewind IUNIT

      END SUBROUTINE REWIND_PARALLEL

      SUBROUTINE SKIP_PARALLEL(IUNIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT

      If (AM_I_ROOT()) read(IUNIT)

      END SUBROUTINE SKIP_PARALLEL

      SUBROUTINE DREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,AVAR,
     &     recs_to_skip)
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME      !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:) !@var  AOUT real*8 array
      INTEGER, INTENT(IN), OPTIONAL :: recs_to_skip
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array

      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      INTEGER :: rc
    ! now local

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
c      write(*,*) "DD dread parallel IM,JM",IM,JM
      If (AM_I_ROOT()) then
         if(present(recs_to_skip)) then
           do n=1,recs_to_skip
             READ (IUNIT,IOSTAT=rc)
           enddo
         endif
         READ (IUNIT,IOSTAT=rc) AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     AOUT(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif
      if (AM_I_ROOT()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ': IOSTAT=',rc
            call stop_model('DREAD_PARALLEL_2D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE DREAD_PARALLEL_2D

      SUBROUTINE DREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR,
     &     recs_to_skip)
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME        !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var  AOUT real*8 array
      INTEGER, INTENT(IN), OPTIONAL :: recs_to_skip
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      INTEGER :: rc
    ! now local
#ifdef USE_ESMF
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)
c      write(*,*) "DD dread parallel IM,JM,NM",IM,JM,NM

      If (AM_I_ROOT()) then
         if(present(recs_to_skip)) then
           do n=1,recs_to_skip
             READ (IUNIT,IOSTAT=rc)
           enddo
         endif
         READ (IUNIT,IOSTAT=rc) AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=
     &     AOUT(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      if (AM_I_ROOT()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ': IOSTAT=',rc
            call stop_model('DREAD_PARALLEL_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE DREAD_PARALLEL_3D

      SUBROUTINE DREAD8_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR,
     &     recs_to_skip)
!@sum DREAD_PARALLEL  read an array real*8 avar(im,jm,:)
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME        !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var AVAR real*8 array
      INTEGER, INTENT(IN), OPTIONAL :: recs_to_skip
      REAL*8 :: AGLOB(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var AGLOB global array
      INTEGER :: N,rc

      If (AM_I_ROOT()) then
        if(present(recs_to_skip)) then
          do n=1,recs_to_skip
            READ (IUNIT,IOSTAT=rc)
          enddo
        endif
        READ (IUNIT,IOSTAT=rc) AGLOB
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AGLOB, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=
     &     AGLOB(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      if (AM_I_ROOT()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ': IOSTAT=',rc
            call stop_model('DREAD8_PARALLEL_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE DREAD8_PARALLEL_3D

      SUBROUTINE MREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,M,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME      !@var  NAME  name of file being read
      INTEGER,      INTENT(OUT) :: M         !@var  M      initial integer
      REAL*8,      INTENT(OUT)  :: AVAR(:,grd_dum%J_STRT_HALO:) !@var  AOUT real*8 array
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      INTEGER :: rc
    ! now local
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=rc) M, AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AOUT, AVAR, shape(AVAR), 2)
      CALL ESMF_BCAST(grd_dum, M   )
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     AOUT(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      if (AM_I_ROOT()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ': IOSTAT=',rc
            call stop_model('MREAD_PARALLEL_2D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE MREAD_PARALLEL_2D

      SUBROUTINE MREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,M,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME        !@var  NAME  name of file being read
      INTEGER,      INTENT(OUT) :: M           !@var  M      initial integer
      REAL*8,      INTENT(OUT)  :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var  AOUT real*8 array
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      INTEGER :: rc 
    ! now local

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=rc) M, AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AOUT, AVAR, shape(AVAR), 2)
      CALL ESMF_BCAST(grd_dum, M   )
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=
     &     AOUT(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      if (AM_I_ROOT()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ': IOSTAT=',rc
            call stop_model('MREAD_PARALLEL_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE MREAD_PARALLEL_3D

      SUBROUTINE READT_PARALLEL_2D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS       !@var  IPOS  no. of recs. to advance
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: rc
    ! now local

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=rc)
         END DO
         READ (IUNIT, IOSTAT=rc) TITLE, AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     AOUT(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      if (AM_I_ROOT()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ': ',
     &           TRIM(TITLE),' IOSTAT=',rc
            call stop_model('READT_PARALLEL_2D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT_PARALLEL_2D

      SUBROUTINE READT_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: rc
    ! now local
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=rc)
         END DO
         READ (IUNIT, IOSTAT=rc) TITLE, AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=
     &     AOUT(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      if (am_i_root()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ':
     &           ',TRIM(TITLE),' IOSTAT=',rc
            call stop_model('READT_PARALLEL_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT_PARALLEL_3D

      SUBROUTINE READT_PARALLEL_COLUMN_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT_PARALLEL_COLUMN read in real*4 (:,im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,:,grd_dum%J_STRT_HALO:)  !@var  AVAR real*8 array
      INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
      REAL*4 :: AGLOB4(size(AVAR,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var AGLOB global array
      REAL*8 :: AGLOB (size(AVAR,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var AGLOB global array
      INTEGER :: N                        !@var  N loop variable
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: rc

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=rc)
         END DO
         READ (IUNIT, IOSTAT=rc) TITLE, AGLOB4
C****  convert from real*4 to real*8
         AGLOB = AGLOB4
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AGLOB, AVAR, shape(AVAR), 3)
#else
      AVAR(:,:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     AGLOB(:,:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      if (am_i_root()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ':
     &           ',TRIM(TITLE),' IOSTAT=',rc
            call stop_model('READT_PARALLEL_COLUMN_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT_PARALLEL_COLUMN_3D

      SUBROUTINE READT8_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT8_PARALLEL read in real*8 (im,jm,:) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
      REAL*8 :: AGLOB(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var AGLOB global array
      INTEGER :: N                        !@var  N loop variable
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: rc

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=rc)
         END DO
         READ (IUNIT, IOSTAT=rc) TITLE, AGLOB
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AGLOB, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=
     &     AGLOB(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      if (am_i_root()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ':
     &           ',TRIM(TITLE),' IOSTAT=',rc
            call stop_model('READT8_PARALLEL_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT8_PARALLEL_3D

      SUBROUTINE READT8_COLUMN_3D (grd_dum,IUNIT,NAME,AVAR,IPOS)
!@sum READT8_COLUMN read in real*8 (:,im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME         !@var  NAME  name of file being read
      REAL*8,       INTENT(OUT) :: AVAR(:,:,grd_dum%J_STRT_HALO:)  !@var  AVAR real*8 array
      INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
      REAL*8 :: AGLOB(size(AVAR,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var AGLOB global array
      INTEGER :: N                        !@var  N loop variable
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: rc

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=rc)
         END DO
         READ (IUNIT, IOSTAT=rc) TITLE, AGLOB
      EndIf

#ifdef USE_MPI
      Call scatter(grd_dum, AGLOB, AVAR, shape(AVAR), 3)
#else
      AVAR(:,:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     AGLOB(:,:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      if (am_i_root()) then
         If (rc==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME), ':
     &           ',TRIM(TITLE),' IOSTAT=',rc
            call stop_model('READT8_COLUMN_3D: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT8_COLUMN_3D

      SUBROUTINE WRITET8_COLUMN_3D (grd_dum,IUNIT,NAME,buf,title)
!@sum WRITET8_COLUMN write title*80, real*8 buf(:,im,jm)
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being read
      REAL*8,       INTENT(IN) :: buf(:,:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
      CHARACTER*80,  INTENT(IN)  :: title
      REAL*8 :: buf_glob(size(buf,1),grd_dum%IM_WORLD,grd_dum%JM_WORLD)
      INTEGER :: rc

#ifdef USE_MPI
      Call gather(grd_dum, buf, buf_glob, shape(buf), 3)
#else
      buf_glob(:,:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     buf(:,:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      If (AM_I_ROOT()) then
        WRITE (IUNIT, IOSTAT=rc) title, buf_glob
         If (rc==0) Then
            WRITE(6,*) "Wrote to file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'WRITE ERROR ON FILE ', NAME, ' IOSTAT=',rc
            call stop_model('WRITET8_COLUMN: WRITE ERROR',255)
         EndIf
      end if

      END SUBROUTINE WRITET8_COLUMN_3D

      SUBROUTINE WRITEI_PARALLEL_2D (grd_dum,IUNIT,NAME,buf,it)
!@sum WRITEI_PARALLEL  Parallel version of UTILDBL.f:WRITEI for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of record being read
      REAL*4,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
      INTEGER,      INTENT(IN)  :: it       !@var  it iteration
      REAL*8 :: buf8(grd_dum%IM_WORLD,
     &     grd_dum%J_STRT_HALO:grd_dum%J_STOP_HALO)
!@var buf_glob real*4 array
      REAL*4 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD)
      REAL*8 :: buf_glob8(grd_dum%IM_WORLD,grd_dum%JM_WORLD)
      INTEGER :: rc

!!! not sure if it is implemented for real*4 ...
#ifdef USE_MPI
c***      Call gather(grd_dum%ESMF_GRID, buf, buf_glob, shape(buf), 2)
      !buf_glob = -99999999. ! not implemented
      buf8 = buf
      Call gather(grd_dum, buf8, buf_glob8, shape(buf), 2)
      buf_glob = buf_glob8
#else
      buf_glob(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     buf(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      If (AM_I_ROOT()) then
        WRITE (IUNIT, IOSTAT=rc) it, buf_glob, it
         If (rc==0) Then
            WRITE(6,*) "Wrote to file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'WRITE ERROR ON FILE ', NAME, ' IOSTAT=',rc
            call stop_model('WRITEI_PARALLEL: WRITE ERROR',255)
         EndIf
      end if

      END SUBROUTINE WRITEI_PARALLEL_2D

      SUBROUTINE WRITET_PARALLEL_2D (grd_dum,IUNIT,NAME,buf,title)
!@sum WRITET_PARALLEL  write character*80 title, real(buf(1:im,1:jm),kind=4)
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being written
      REAL*8,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
      CHARACTER*80,  INTENT(IN)  :: title
      REAL*8 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD)
      INTEGER :: rc

#ifdef USE_MPI
      Call gather(grd_dum, buf, buf_glob, shape(buf), 2)
#else
      buf_glob(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     buf(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      If (AM_I_ROOT()) then
        WRITE (IUNIT, IOSTAT=rc) title, real(buf_glob,kind=4)
         If (rc==0) Then
            WRITE(6,*) "Wrote to file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'WRITE ERROR ON FILE ', NAME, ' IOSTAT=',rc
            call stop_model('WRITET_PARALLEL: WRITE ERROR',255)
         EndIf
      end if

      END SUBROUTINE WRITET_PARALLEL_2D

      SUBROUTINE WRITEI8_PARALLEL_3D (grd_dum,IUNIT,NAME,buf,it)
!@sum WRITEI8_PARALLEL  Parallel version of UTILDBL.f:WRITEI8 for (im,jm,:) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being written
      REAL*8,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:,:)  !@var  buf real*8 array
      INTEGER,      INTENT(IN)  :: it       !@var  it iteration
      REAL*8 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(buf,3))
      INTEGER :: rc 

#ifdef USE_MPI
      Call gather(grd_dum, buf, buf_glob, shape(buf), 2)
#else
      buf_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:) =
     &     buf(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      If (AM_I_ROOT()) then
        WRITE (IUNIT, IOSTAT=rc) it, buf_glob, it
         If (rc==0) Then
            WRITE(6,*) "Wrote to file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'WRITE ERROR ON FILE ', NAME, ' IOSTAT=',rc
            call stop_model('WRITEI8_PARALLEL: WRITE ERROR',255)
         EndIf
      end if

      END SUBROUTINE WRITEI8_PARALLEL_3D

      SUBROUTINE DWRITE8_PARALLEL_3D (grd_dum,IUNIT,NAME,buf)
!@sum DWRITE8_PARALLEL  Writes real*8 (im,jm,:) arrays as real*8 on disk
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*(*), INTENT(IN)  :: NAME       !@var  NAME  name of file being written
      REAL*8,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:,:)  !@var  buf real*8 array
      REAL*8, dimension (:,:,:), allocatable :: buf_glob ! global array written to disk
      INTEGER :: rc

      if(am_i_root()) then
         allocate(
     &     buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(buf,3)))
      else
         allocate(buf_glob(1,1,1))
      end if
#ifdef USE_MPI
      Call gather(grd_dum, buf, buf_glob, shape(buf), 2)
#else
      buf_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:) =
     &     buf(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      If (AM_I_ROOT()) then
        WRITE (IUNIT, IOSTAT=rc) buf_glob
        deallocate(buf_glob)
         If (rc==0) Then
            WRITE(6,*) "Wrote to file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'WRITE ERROR ON FILE ', NAME, ' IOSTAT=',rc
            call stop_model('DWRITE8_PARALLEL: WRITE ERROR',255)
         EndIf
      else
         deallocate(buf_glob)
      end if

      END SUBROUTINE DWRITE8_PARALLEL_3D


      subroutine ESMF_IArrayScatter_IJ(egrid, local_array, global_array)
      integer      , dimension (:,:) :: local_array, global_array
!     type (ESMF_Grid)      :: egrid
      type (DIST_GRID)      :: egrid

#ifdef USE_MPI
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            ::
     &     sendcounts, displs
      integer                                       :: nDEs
      integer                                       :: rc
      integer                                       :: ierr
      integer                                       :: recvcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      integer      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(egrid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(egrid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

      allocate (sendcounts(NPES), displs(0:NPES), stat=rc)

      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=rc)
      else
        allocate(VAR(1))
      endif

      displs(0) = 0
      do I = 1,NPES
        J = I - 1
        I1 = AI(I,1)%min
        IN = AI(I,1)%max
        J1 = AI(I,2)%min
        JN = AI(I,2)%max

        sendcounts(I) = (IN - I1 + 1) * (JN - J1 + 1)
        if (J == MY_PET) then
          recvcount = sendcounts(I)
        endif
        displs(I) = displs(J) + sendcounts(I)
        if (AM_I_ROOT()) then

          var(displs(J):displs(I)-1) =
     &         RESHAPE(global_array(I1:IN,J1:JN),
     &         shape=(/sendcounts(I)/))
        endif
      enddo

      call MPI_ScatterV(var, sendcounts, displs,
     &     MPI_INTEGER         , local_array, recvcount,
     &     MPI_INTEGER         , root, MPI_COMM_WORLD, ierr)

        deallocate(VAR, stat=rc)

      deallocate(sendcounts, displs, stat=rc)
      deallocate(AI)
#else
      local_array=global_array
#endif

      end subroutine ESMF_IArrayScatter_IJ

!-----------------------------------------------------------
      subroutine ESMF_LArrayScatter_IJ(egrid, local_array, global_array)
      logical      , dimension (:,:) :: local_array, global_array
!     type (ESMF_Grid)      :: egrid
      type (DIST_GRID)      :: egrid

#ifdef USE_MPI
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            ::
     &     sendcounts, displs
      integer                                       :: nDEs
      integer                                       :: rc
      integer                                       :: ierr
      integer                                       :: recvcount

      integer                                       :: I, J, II, JJ, III
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      Logical      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(egrid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(egrid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif
      allocate (sendcounts(NPES), displs(0:NPES), stat=rc)

      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=rc)
      else
        allocate(VAR(1))
      endif

      displs(0) = 0
      do I = 1,NPES
        J = I - 1
        I1 = AI(I,1)%min
        IN = AI(I,1)%max
        J1 = AI(I,2)%min
        JN = AI(I,2)%max

        sendcounts(I) = (IN - I1 + 1) * (JN - J1 + 1)
        if (J == MY_PET) then
          recvcount = sendcounts(I)
        endif
        displs(I) = displs(J) + sendcounts(I)
        if (AM_I_ROOT()) then


          var(displs(J):displs(I)-1) =
     &         RESHAPE(global_array(I1:IN,J1:JN),
     &         shape=(/sendcounts(I)/))


        endif
      enddo

      call MPI_ScatterV(var, sendcounts, displs,
     &     MPI_LOGICAL         , local_array, recvcount,
     &     MPI_LOGICAL         , root, MPI_COMM_WORLD, ierr)

        deallocate(VAR, stat=rc)

      deallocate(sendcounts, displs, stat=rc)
      deallocate(AI)
#else
      local_array=global_array
#endif
      end subroutine ESMF_LArrayScatter_IJ

#ifdef USE_MPI
!--------------------------------
      subroutine Esmf_IArrayGather_IJ(grid, local_array, global_array)
!     type (ESMF_Grid)      :: grid
      type (DIST_GRID)      :: grid
      integer      , dimension (:,:) :: local_array, global_array

      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      integer, allocatable, dimension(:)            ::
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: rc
      integer                                       :: ierr
      integer                                       :: sendcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      integer      , allocatable                    :: var(:)


      Allocate(AI(1:NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

      allocate (recvcounts(NPES), displs(0:NPES), stat=rc)

      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=rc)

      displs(0) = 0
      do I = 1,NPES
        J = I - 1
        I1 = AI(I,1)%min
        IN = AI(I,1)%max
        J1 = AI(I,2)%min
        JN = AI(I,2)%max

        recvcounts(I) = (IN - I1 + 1) * (JN - J1 + 1)
        if (J == MY_PET) then
          sendcount = recvcounts(I)
        endif
        displs(I) = displs(J) + recvcounts(I)
      enddo


      call MPI_GatherV(local_array, sendcount, MPI_INTEGER         ,
     &     var, recvcounts, displs,
     &     MPI_INTEGER         , root, MPI_COMM_WORLD, ierr)

      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          I1 = AI(I,1)%min
          IN = AI(I,1)%max
          J1 = AI(I,2)%min
          JN = AI(I,2)%max

          global_array(I1:IN,J1:JN) =
     &         RESHAPE(var(displs(J):displs(I)-1),
     &         shape=(/IN-I1+1,
     &         JN-J1+1/))
        enddo
      endif
      deallocate(VAR, stat=rc)
      deallocate(recvcounts, displs, stat=rc)
      DEALLOCATE(AI)


      end subroutine Esmf_IArrayGather_IJ

!--------------------------------
      subroutine Esmf_LArrayGather_IJ(e_grid, local_array, global_array)
!     type (ESMF_Grid)      :: e_grid
      type (DIST_GRID)      :: e_grid
      logical      , dimension (:,:) :: local_array, global_array

      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      integer, allocatable, dimension(:)            ::
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: rc
      integer                                       :: ierr
      integer                                       :: sendcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN

      logical      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(e_grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(e_grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif
      allocate (recvcounts(NPES), displs(0:NPES), stat=rc)

      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=rc)


      displs(0) = 0
      do I = 1,NPES
        J = I - 1
        I1 = AI(I,1)%min
        IN = AI(I,1)%max
        J1 = AI(I,2)%min
        JN = AI(I,2)%max

        recvcounts(I) = (IN - I1 + 1) * (JN - J1 + 1)
        if (J == MY_PET) then
          sendcount = recvcounts(I)
        endif
        displs(I) = displs(J) + recvcounts(I)
      enddo

      call MPI_GatherV(local_array, sendcount, MPI_LOGICAL,
     &     var, recvcounts, displs,
     &     MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)

      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          I1 = AI(I,1)%min
          IN = AI(I,1)%max
          J1 = AI(I,2)%min
          JN = AI(I,2)%max

          global_array(I1:IN,J1:JN) =
     &         RESHAPE(var(displs(J):displs(I)-1),
     &         shape=(/IN-I1+1,
     &         JN-J1+1/))
        enddo
      endif

      deallocate(VAR, stat=rc)
      deallocate(recvcounts, displs, stat=rc)
      Deallocate(AI)
      end subroutine Esmf_LArrayGather_IJ
!------------------------------------------------------------

      subroutine Esmf_ArrayGather_J_int(grid, local_array, global_array)
!     type (ESMF_Grid)      :: grid
      type (DIST_GRID)      :: grid
      INTEGER, dimension (:) :: local_array, global_array

      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI

      integer, allocatable, dimension(:)            ::
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: rc 
      integer                                       :: ierr
      integer                                       :: sendcount

      integer                                       :: I,J
      integer                                       :: NY
      integer                                      :: J1, JN


      integer, allocatable                    :: var(:)


      Allocate(AI(1:NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

      allocate (recvcounts(NPES), displs(0:NPES), stat=rc)

      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=rc)

      displs(0) = 0
      do I = 1,NPES
        J = I - 1
        J1 = AI(I,2)%min
        JN = AI(I,2)%max

        recvcounts(I) = (JN - J1 + 1)
        if (J == MY_PET) then
          sendcount = recvcounts(I)
        endif
        displs(I) = displs(J) + recvcounts(I)
      enddo


      call MPI_GatherV(local_array,sendcount,MPI_INTEGER,
     &     var, recvcounts, displs,
     &     MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          J1 = AI(I,2)%min
          JN = AI(I,2)%max

          global_array(J1:JN) = var(displs(J):displs(I)-1)
        enddo
      endif

      deallocate(VAR, stat=rc)
      deallocate(recvcounts, displs, stat=rc)
      Deallocate(AI)

      end subroutine Esmf_ArrayGather_J_int

#endif
!---------------------------

      subroutine ArrayGather_J(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:)               :: global_array

#ifdef USE_MPI
      Call Gather(grd_dum, local_array, global_array,
     &     shape(local_array), 1)
#else
      global_array(grd_dum%J_STRT:grd_dum%J_STOP) =
     &     local_array(grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine ArrayGather_J

!---------------------------

      subroutine ArrayGather_J_int(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      INTEGER, dimension (grd_dum%J_STRT_HALO:) :: local_array
      INTEGER, dimension (:)               :: global_array

#ifdef USE_MPI
      call Esmf_ArrayGather_J_int(grd_dum,
     &     local_array(grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array(grd_dum%J_STRT:grd_dum%J_STOP) =
     &     local_array(grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine ArrayGather_J_int

!---------------------------

      subroutine ArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:,:)               :: global_array

#ifdef USE_MPI
      Call gather(grd_dum, local_array, global_array,
     &     shape(local_array), 2)

#else
      global_array(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine ArrayGather_IJ

!---------------------------

      subroutine LArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      logical, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      logical, dimension (:,:)                    :: global_array

#ifdef USE_MPI
          call Esmf_LArrayGather_IJ(grd_dum,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &         local_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine LArrayGather_IJ

!---------------------------


      subroutine IArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      integer, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      integer, dimension (:,:)               :: global_array

#ifdef USE_MPI
      call Esmf_IArrayGather_IJ(grd_dum,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine IArrayGather_IJ


    !---------------------------
!---------------------------

      subroutine ArrayScatter_J(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:)               :: global_array

#ifdef USE_MPI
      Call scatter(grd_dum, global_array, local_array,
     &     shape(local_array), 1)
#else
      local_array(grd_dum%J_STRT:grd_dum%J_STOP) =
     &     global_array(grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine ArrayScatter_J

!---------------------------

      subroutine ArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:,:)               :: global_array

#ifdef USE_MPI
      Call scatter(grd_dum, global_array, local_array,
     &     shape(local_array), 2)
#else
      local_array(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     global_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine ArrayScatter_IJ

!---------------------------

      subroutine IArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      integer, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      integer, dimension (:,:)               :: global_array

#ifdef USE_MPI
      call Esmf_IArrayScatter_IJ(grd_dum,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     & global_array)
#else
      local_array(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     global_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine IArrayScatter_IJ

!---------------------------

      subroutine LArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      logical, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      logical, dimension (:,:)               :: global_array

#ifdef USE_MPI
      call Esmf_LArrayScatter_IJ(grd_dum,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     & global_array)
#else
      local_array(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     global_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine LArrayScatter_IJ


    !---------------------------
#ifdef USE_MPI
      subroutine ESMF_GRID_PE_LAYOUT  (GRID, NX, NY)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: NX, NY

        NX = 1
        NY = NPES

      end subroutine ESMF_GRID_PE_LAYOUT

    !---------------------------

      subroutine ESMF_GRID_MY_PE_LOC  (GRID,NX0,NY0)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)          :: NX0, NY0

    ! local vars
        type (ESMF_DELayout) :: layout

        integer :: rc

#ifdef USE_ESMF
        call ESMF_GridGet(grid, delayout=layout, rc=rc)
#endif
        NX0 = 0
        NY0 = my_pet

      end subroutine ESMF_GRID_MY_PE_LOC

    !---------------------------
#endif

      subroutine WRITE_PARALLEL_INTEGER_0 ( data, UNIT, format, CRIT)

        INTEGER, intent(in   )            :: data
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        logical,            intent(in   ), optional  :: CRIT

        character(len=ESMF_MAXSTR) :: FORMATTED

        logical :: crit_,crit2_ ! local

        crit_ = .false.
        crit2_= .false.
        if (present(CRIT)) then
          crit_ = crit
        else
          crit2_ = AM_I_ROOT()
        endif

        if (crit_ .or. crit2_) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(format)) then;    write(UNIT, format) data
                 else;    write(UNIT, *     ) data
                 end if
              elseif(FORMATTED == "NO") then; write(UNIT        ) data
              end if
           else
              if (present(format)) then; write(*, format) data
              else; write(*,      *) data
              end if
           end if
        end if

        return

      end subroutine WRITE_PARALLEL_INTEGER_0

    !---------------------------
      subroutine WRITE_PARALLEL_INTEGER_1 ( data, UNIT, format, CRIT)

        INTEGER, intent(in   )            :: data (:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        logical,            intent(in   ),  optional :: CRIT

        character(len=ESMF_MAXSTR) :: FORMATTED
        logical :: crit_,crit2_ ! local

        crit_ = .false.
        crit2_= .false.
        if (present(CRIT)) then
          crit_ = crit
        else
          crit2_ = AM_I_ROOT()
        endif

        if (crit_ .or. crit2_) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(format)) then;    write(UNIT, format) data
                 else;    write(UNIT, *     ) data
                 end if
              elseif(FORMATTED == "NO") then; write(UNIT        ) data
              end if
           else
              if (present(format)) then; write(*, format) data
              else; write(*,      *) data
              end if
           end if
        end if

        return

      end subroutine WRITE_PARALLEL_INTEGER_1

    !---------------------------

      subroutine WRITE_PARALLEL_REAL8_0 ( data, UNIT, format, CRIT)

        REAL (KIND=8), intent(in   )            :: data
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        logical,            intent(in   ),  optional :: CRIT

        character(len=ESMF_MAXSTR) :: FORMATTED

        logical :: crit_,crit2_ ! local

        crit_ = .false.
        crit2_= .false.
        if (present(CRIT)) then
          crit_ = crit
        else
          crit2_ = AM_I_ROOT()
        endif

        if (crit_ .or. crit2_) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(format)) then;    write(UNIT, format) data
                 else;    write(UNIT, *     ) data
                 end if
              elseif(FORMATTED == "NO") then; write(UNIT        ) data
              end if
           else
              if (present(format)) then; write(*, format) data
              else; write(*,      *) data
              end if
           end if
        end if

        return

      end subroutine WRITE_PARALLEL_REAL8_0

    !---------------------------
      subroutine WRITE_PARALLEL_REAL8_1 ( data, UNIT, format, CRIT)

        REAL (KIND=8), intent(in   )            :: data(:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        logical,            intent(in   ), optional  :: CRIT

        character(len=ESMF_MAXSTR) :: FORMATTED

        logical :: crit_,crit2_ ! local

        crit_ = .false.
        crit2_= .false.
        if (present(CRIT)) then
          crit_ = crit
        else
          crit2_ = AM_I_ROOT()
        endif

        if (crit_ .or. crit2_) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(format)) then;    write(UNIT, format) data
                 else;    write(UNIT, *     ) data
                 end if
              elseif(FORMATTED == "NO") then; write(UNIT        ) data
              end if
           else
              if (present(format)) then; write(*, format) data
              else; write(*,      *) data
              end if
           end if
        end if

        return

      end subroutine WRITE_PARALLEL_REAL8_1

    !---------------------------

      subroutine WRITE_PARALLEL_STRING_0 (STRING, FORMAT, UNIT, CRIT)
        character(LEN=*), intent(in   )            :: STRING
        character(LEN=*), intent(in   ), optional  :: FORMAT
        integer,          intent(in   ),  optional :: UNIT
        logical,          intent(in   ), optional  :: CRIT
        logical :: crit_,crit2_ ! local

        crit_ = .false.
        crit2_= .false.
        if (present(CRIT)) then
          crit_ = crit
        else
          crit2_ = AM_I_ROOT()
        endif

        if (crit_ .or. crit2_) then
          if (present(UNIT)) then
            if (present(FORMAT)) then
              write(UNIT, FMT=FORMAT) STRING
            else
              write(UNIT, FMT='(A)') STRING
            endif
          else
            if (present(FORMAT)) then
              write(*, FMT=FORMAT) STRING
            else
              write(*, FMT='(A)') STRING
            endif
          end if
        end if

      end subroutine WRITE_PARALLEL_STRING_0

    !---------------------------

      subroutine WRITE_PARALLEL_STRING_1 (STRING, FORMAT, UNIT, CRIT)
        character(LEN=*), intent(in   )            :: STRING(:)
        character(LEN=*), intent(in   ), optional  :: FORMAT
        integer,          intent(in   ),  optional :: UNIT
        logical,          intent(in   ), optional  :: CRIT
        logical :: crit_,crit2_ ! local

        crit_ = .false.
        crit2_= .false.
        if (present(CRIT)) then
          crit_ = crit
        else
          crit2_ = AM_I_ROOT()
        endif

        if (crit_ .or. crit2_) then
          if (present(UNIT)) then
            if (present(FORMAT)) then
              write(UNIT, FMT=FORMAT) STRING
            else
              write(UNIT, FMT='(A)') STRING
            endif
          else
            if (present(FORMAT)) then
              write(*, FMT=FORMAT) STRING
            else
              write(*, FMT='(A)') STRING
            endif
          end if
        end if

      end subroutine WRITE_PARALLEL_STRING_1

    !---------------------------

      subroutine READ_PARALLEL_INTEGER_0 ( DATA, UNIT, FORMAT)
    ! Wrapper for ESMF_READ_PARALLEL_INTEGER_0

        integer, intent(out  )                       :: DATA
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: FORMAT

   ! now local
        type (ESMF_DELayout)                         :: layout

        layout=ESMF_LAYOUT
        call ESMF_READ_PARALLEL_INTEGER_0( layout, DATA, UNIT, FORMAT)

      END SUBROUTINE READ_PARALLEL_INTEGER_0

    !---------------------------
      subroutine ESMF_READ_PARALLEL_INTEGER_0(layout,DATA,UNIT,FORMAT)

        type (ESMF_DELayout)                         :: layout
        integer, intent(out  )                       :: DATA
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: FORMAT

        character(len=ESMF_MAXSTR) :: FORMATTED
        character(LEN=ESMF_MAXSTR) :: FILENAME
        logical                :: IS_NAMED
        integer                :: IOSTAT
        integer                :: ierr

        if (ESMF_AM_I_ROOT(layout)) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(FORMAT)) then
                    read(UNIT, FORMAT, IOSTAT=IOSTAT) DATA
                 else
                    read(UNIT, *, IOSTAT=IOSTAT) DATA
                 end if
              elseif(FORMATTED == "NO") then
                 read(UNIT,       IOSTAT=IOSTAT) DATA
              end if

              if (IOSTAT < 0) then
                 inquire(unit=UNIT, NAMED=IS_NAMED, NAME=FILENAME)
                 if (IS_NAMED) then
                    print *, "Premature end of file ",FILENAME
                    goto 9999
                 else
                    print *, "Premature end of file UNNAMED"
                    goto 9999
                 end if
              end if
           else
              if (present(FORMAT)) then; read(*, FORMAT ) DATA
              else; read          *, DATA
              end if
           end if
        end if

#ifdef USE_MPI
        CALL MPI_BCAST(data, 1, MPI_INTEGER, root,
     &       MPI_COMM_WORLD, ierr)
#endif

        return
 9999 continue
        call stop_model('ESMF_READ_PARALLEL: READ ERROR',255)
      END SUBROUTINE ESMF_READ_PARALLEL_INTEGER_0

    !---------------------------

      subroutine READ_PARALLEL_INTEGER_1 ( DATA, UNIT, FORMAT)
    ! Wrapper for ESMF_READ_PARALLEL_INTEGER_1

        integer, intent(out  )                       :: DATA(:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: FORMAT

   ! now local
        type (ESMF_DELayout)                         :: layout

        layout = ESMF_LAYOUT
        call ESMF_READ_PARALLEL_INTEGER_1 ( layout, DATA, UNIT, FORMAT)

      END SUBROUTINE READ_PARALLEL_INTEGER_1

    !---------------------------
      subroutine ESMF_READ_PARALLEL_INTEGER_1(layout,DATA,UNIT,FORMAT)

        type (ESMF_DELayout)                         :: layout
        integer, intent(out  )                       :: DATA(:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: FORMAT

        character(len=ESMF_MAXSTR) :: FORMATTED
        character(LEN=ESMF_MAXSTR) :: FILENAME
        logical                :: IS_NAMED
        integer                :: IOSTAT
        integer                :: ierr

        if (ESMF_AM_I_ROOT(layout)) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(FORMAT)) then
                    read(UNIT, FORMAT, IOSTAT=IOSTAT) DATA
                 else
                    read(UNIT, *, IOSTAT=IOSTAT) DATA
                 end if
              elseif(FORMATTED == "NO") then
                 read(UNIT,       IOSTAT=IOSTAT) DATA
              end if

              if (IOSTAT < 0) then
                 inquire(unit=UNIT, NAMED=IS_NAMED, NAME=FILENAME)
                 if (IS_NAMED) then
                    print *, "Premature end of file ",FILENAME
                    goto 9999
                 else
                    print *, "Premature end of file UNNAMED"
                    goto 9999
                 end if
              end if
           else
              if (present(FORMAT)) then; read(*, FORMAT ) DATA
              else; read          *, DATA
              end if
           end if
        end if

#ifdef USE_MPI
        CALL MPI_BCAST(data, size(data), MPI_INTEGER, root,
     &       MPI_COMM_WORLD, ierr)
#endif

        return
 9999 continue
        call stop_model('ESMF_READ_PARALLEL: READ ERROR',255)
      END SUBROUTINE ESMF_READ_PARALLEL_INTEGER_1

    !---------------------------

      subroutine READ_PARALLEL_REAL8_1 ( DATA, UNIT, FORMAT)
    ! Wrapper for READ_PARALLEL_REAL8_1

        real(kind=8),       intent(out  )            :: DATA(:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: FORMAT

  ! now local
        type (ESMF_DELayout)                         :: layout

        layout = ESMF_LAYOUT
        call ESMF_READ_PARALLEL_REAL8_1 ( layout, DATA, UNIT, FORMAT)

      END SUBROUTINE READ_PARALLEL_REAL8_1

    !---------------------------

      subroutine ESMF_READ_PARALLEL_REAL8_1 (layout,DATA,UNIT,FORMAT)

        type (ESMF_DELayout)                         :: layout
        real(kind=8),       intent(out  )            :: DATA(:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: FORMAT

        character(len=ESMF_MAXSTR) :: FORMATTED
        character(LEN=ESMF_MAXSTR) :: FILENAME
        logical                :: IS_NAMED
        integer                :: IOSTAT
        integer                :: ierr

        if (ESMF_AM_I_ROOT(layout)) then
           if (present(UNIT)) then
              inquire(unit=UNIT, formatted=FORMATTED)
              if   (FORMATTED == "YES") then
                 if (present(FORMAT)) then
                    read(UNIT, FORMAT, IOSTAT=IOSTAT) DATA
                 else
                    read(UNIT, *, IOSTAT=IOSTAT) DATA
                 end if
              elseif(FORMATTED == "NO") then
                 read(UNIT,       IOSTAT=IOSTAT) DATA
              end if

              if (IOSTAT < 0) then
                 inquire(unit=UNIT, NAMED=IS_NAMED, NAME=FILENAME)
                 if (IS_NAMED) then
                    print *, "Premature end of file ",FILENAME
                    goto 9999
                 else
                    print *, "Premature end of file UNNAMED"
                    goto 9999
                 end if
              end if
           else
              if (present(FORMAT)) then; read(*, FORMAT ) DATA
              else; read          *, DATA
              end if
           end if
        end if

#ifdef USE_MPI
        CALL MPI_BCAST(data, Size(data), MPI_DOUBLE_PRECISION, root,
     &       MPI_COMM_WORLD, ierr)
#endif

        return
 9999 continue
        call stop_model('ESMF_READ_PARALLEL: READ ERROR',255)
      END SUBROUTINE ESMF_READ_PARALLEL_REAL8_1

    !---------------------------

#ifdef USE_ESMF
      subroutine ESMF_GRID_BOUNDS(GRID, X_LOC, Y_LOC, I1, IN, J1, JN)
        type (DIST_Grid), intent(INOUT) :: grid
        integer, intent(IN)          :: X_LOC, Y_LOC
        integer, intent(OUT)         :: I1, IN, J1, JN
        integer :: i
    ! local vars
        integer :: deId
        integer :: rc

        type(ESMF_AxisIndex), dimension(:,:), pointer :: AI

        Allocate(AI(1:NPES,3))
        call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, AI,
     &       horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
        deId = ESMF_GRID_PE_NUM_FM_PE_LOC(grid%ESMF_GRID, X_LOC, Y_LOC)


    ! AI uses 1-based index for deId
        deId = deId + 1

        I1 = AI(deId,1)%min
        IN = AI(deId,1)%max
        J1 = AI(deId,2)%min
        JN = AI(deId,2)%max

        do i=0, npes-1
          grid%dj_map(i)=AI(i+1,2)%max - AI(i+1,2)%min + 1
        end do
        grid%dj=grid%dj_map(deid-1)

        Deallocate(AI)
      end subroutine ESMF_GRID_BOUNDS
#endif
    !---------------------------

      function ESMF_AM_I_ROOT(layout) result(R)
        type (ESMF_DELayout) :: layout
        logical              :: R

#ifdef USE_MPI
        R = .false.
        if (my_pet == root) R = .true.
#else
        R = .true.
#endif

      end function ESMF_AM_I_ROOT

    !---------------------------
      function MPI_AM_I_ROOT() result(R)
        logical              :: R

        integer :: rank
        integer :: ierr

#ifdef USE_MPI
        R = .false.
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        if (rank == root) R = .true.
#else
      R = .true.
#endif

      end function MPI_AM_I_ROOT

    !---------------------------

      function ESMF_GRID_PE_NUM_FM_PE_LOC (
     &  GRID                              ,
     &  X_LOC                             ,
     &  Y_LOC                             ) result (deId)
        type (ESMF_GRID), intent(in ) :: GRID
        integer           , intent(in ) :: X_LOC
        integer           , intent(in ) :: Y_LOC
        integer                         :: deId

        type (ESMF_DELayout) :: layout
        integer :: nx, ny
        integer :: rc

#ifdef USE_MPI
#ifdef USE_ESMF
        call ESMF_GridGet(grid, delayout=layout, rc=rc)
#endif
        deid = y_loc
#else
        deid = 1
#endif


      end function ESMF_GRID_PE_NUM_FM_PE_LOC


    !---------------------------

    !---------------------------

#ifdef USE_ESMF
      subroutine ESMF_GRID_WORLD(GRID,IM_WORLD,JM_WORLD)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: IM_WORLD, JM_WORLD

! local vars
        integer rc
        integer, dimension(ESMF_MAXGRIDDIM) :: dims

        call ESMF_GridGet(grid, horzRelLoc=ESMF_CELL_CENTER,
     &     globalcellcountperdim=dims, rc=rc)

        IM_WORLD = dims(1)
        JM_WORLD = dims(2)

      end subroutine ESMF_GRID_WORLD
#endif

!---------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ESMF_DELayoutBarrier(layout, rc)
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

#ifndef USE_MPI
      Subroutine ESMF_GRID_WORLD(grd, im, jm)
      Type (ESMF_GRID) :: grd
      Integer :: im, jm

      im = GRID%im_world
      jm = GRID%jm_world

      End subroutine ESMF_GRID_WORLD
#endif

      SUBROUTINE HERE(file, line)
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
       ALLOCATE(lines(npes))
       Call MPI_Allgather(line, 1, MPI_INTEGER, lines, 1, MPI_INTEGER,
     &      MPI_COMM_WORLD, ierr)
       If (Any(lines /= line))
     &      call stop_model('HERE: synchronization error -severe.',255)
       Deallocate(lines)
#endif
#endif
#ifdef USE_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE HERE

      SUBROUTINE LOG_PARALLEL(grd_dum, file, line, i0, i1, x0, x1)
      Use FILEMANAGER, only : nameunit
      IMPLICIT NONE

      TYPE(DIST_GRID), INTENT(IN) :: grd_dum
      CHARACTER(Len=*) :: file
      INTEGER          :: line

      INTEGER, OPTIONAL :: i0, i1(:)
      REAL*8, OPTIONAL  :: x0, x1(:)

      INTEGER :: iu
      INTEGER :: n

#ifdef DEBUG_DECOMP
      iu = grd_dum%log_unit
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

      SUBROUTINE GLOBALMIN_R(grd_dum, val, val_min)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_min

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_min, 1, MPI_DOUBLE_PRECISION,MPI_MIN,
     &     MPI_COMM_WORLD, ierr)
#else
      val_min = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_R(grd_dum, val, val_max)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_DOUBLE_PRECISION,MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I(grd_dum, val, val_max)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_max

      INTEGER  :: ierr

#ifdef USE_MPI
      CALL MPI_Allreduce(val, val_max, 1, MPI_INTEGER, MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
#else
      val_max = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I_1D(grd_dum, val, val_max)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val(:)
      INTEGER,            INTENT(OUT) :: val_max(:)

      INTEGER  :: n,ierr

#ifdef USE_MPI
      n = size(val)
      CALL MPI_Allreduce(val, val_max, n, MPI_INTEGER, MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
#else
      val_max(:) = val(:)
#endif

      END SUBROUTINE

      SUBROUTINE PACK_1D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(grd_dum%JM_WORLD)

#ifdef USE_MPI
      Call Gather(grd_dum, arr, arr_glob, shape(arr), 1)
#else
      arr_glob(grd_dum%J_STRT:grd_dum%J_STOP) =
     &     arr(grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE PACK_1D

      SUBROUTINE IPACK_1D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:)
      INTEGER, INTENT(OUT) :: ARR_GLOB(grd_dum%JM_WORLD)

      CALL ARRAYGATHER(grd_dum,ARR,ARR_GLOB)

      RETURN
      END SUBROUTINE IPACK_1D

      SUBROUTINE PACK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:)

#ifdef USE_MPI
      Call Gather(grd_dum, arr, arr_glob, shape(arr), 2)
#else
      arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP) =
     & arr(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE PACK_2D

      SUBROUTINE IPACK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, INTENT(INOUT) :: ARR_GLOB(:,:)

      Call arrayGather(grd_dum, arr, arr_glob)

      RETURN
      END SUBROUTINE IPACK_2D

      SUBROUTINE LPACK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      LOGICAL, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      LOGICAL, INTENT(INOUT) :: ARR_GLOB(:,:)

      CALL ARRAYGATHER(grd_dum,ARR,ARR_GLOB)

      RETURN
      END SUBROUTINE LPACK_2D


      SUBROUTINE PACK_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:,:)
      
      INTEGER :: l

#ifdef USE_MPI
      Call Gather(grd_dum, arr, arr_glob, shape(arr), 2)
#else
      arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:) =
     & arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      RETURN
      END SUBROUTINE PACK_3D



      SUBROUTINE IPACK_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER, INTENT(OUT) :: ARR_GLOB(:,:,:)
      INTEGER :: l

      do l=1,size(arr_glob,3)
       CALL ARRAYGATHER(grd_dum,ARR(:,:,l),ARR_GLOB(:,:,l))
      end do

      RETURN
      END SUBROUTINE IPACK_3D



      SUBROUTINE PACK_4D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: l,m

#ifdef USE_MPI
      Call Gather(grd_dum, arr, arr_glob, shape(arr), 2)
#else
      arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:) =
     & arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:)
#endif

      RETURN
      END SUBROUTINE PACK_4D

      SUBROUTINE PACK_5D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:,:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:,:,:,:)
      INTEGER :: l,m

#ifdef USE_MPI
      Call Gather(grd_dum, arr, arr_glob, shape(arr), 2)
#else
      arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:,:) =
     & arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)
#endif

      RETURN
      END SUBROUTINE PACK_5D

      SUBROUTINE PACKj_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:)

      INTEGER :: k

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 1)
#else
      arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:)=
     & arr(grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      RETURN
      END SUBROUTINE PACKj_2D

      SUBROUTINE IPACKj_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      INTEGER, INTENT(INOUT) :: ARR_GLOB(:,:)

      INTEGER :: k

#ifdef USE_MPI
      do k=1,size(arr,2)
        call arraygather(grd_dum,arr(:,k),arr_glob(:,k))
      end do
#else
      arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:)=
     & arr(grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      RETURN
      END SUBROUTINE IPACKj_2D


      SUBROUTINE PACKj_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)

      INTEGER :: k

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 1)
#else
      arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:,:)=
     & arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:)
#endif

      RETURN
      END SUBROUTINE PACKj_3D

      SUBROUTINE PACKj_4D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:,:)

      INTEGER :: k

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 1)
#else
      arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)=
     & arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)
#endif

      RETURN
      END SUBROUTINE PACKj_4D

      SUBROUTINE UNPACK_1D(grd_dum,ARR_GLOB,ARR)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(grd_dum%JM_WORLD)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP)=
     & arr_glob(grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE UNPACK_1D


      SUBROUTINE UNPACK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1
      LOGICAL, OPTIONAL :: local

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 2)
#else
      arr(:,grd_dum%J_STRT:grd_dum%J_STOP)=
     & arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE UNPACK_2D

      SUBROUTINE IUNPACK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) :: ARR_GLOB(:,:)
      INTEGER, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,J_0:J_1)=ARR_GLOB(:,J_0:J_1)
        else
            call arrayscatter(grd_dum, arr(:,:), arr_glob(:,:))
        end if
      else
          call arrayscatter(grd_dum, arr(:,:), arr_glob(:,:))
      end if

      RETURN
      END SUBROUTINE IUNPACK_2D

      SUBROUTINE LUNPACK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      LOGICAL, INTENT(IN) :: ARR_GLOB(:,:)
      LOGICAL, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,J_0:J_1)=ARR_GLOB(:,J_0:J_1)
        else
          call arrayscatter(grd_dum, arr, arr_glob)
        end if
      else
        call arrayscatter(grd_dum, arr, arr_glob)
      end if

      RETURN
      END SUBROUTINE LUNPACK_2D

      SUBROUTINE UNPACK_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0, J_1, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt_halo
          J_1=grd_dum%j_stop_halo
          ARR(:,J_0:J_1,:)=ARR_GLOB(:,J_0:J_1,:)
        else
          do l=1,size(arr,3)
            call arrayscatter(grd_dum, arr(:,:,l), arr_glob(:,:,l))
          end do
        end if
      else
        do l=1,size(arr,3)
          call arrayscatter(grd_dum, arr(:,:,l), arr_glob(:,:,l))
        end do
      end if

      RETURN
      END SUBROUTINE UNPACK_3D


      SUBROUTINE IUNPACK_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) :: ARR_GLOB(:,:,:)
      INTEGER, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0, J_1, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,J_0:J_1,:)=ARR_GLOB(:,J_0:J_1,:)
        else
          do l=1,size(arr,3)
            call arrayscatter(grd_dum, arr(:,:,l), arr_glob(:,:,l))
          end do
        end if
      else
        do l=1,size(arr,3)
          call arrayscatter(grd_dum, arr(:,:,l), arr_glob(:,:,l))
        end do
      end if

      RETURN
      END SUBROUTINE IUNPACK_3D

      SUBROUTINE UNPACK_4D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:,:,:)
      INTEGER :: J_0, J_1, L,M
      LOGICAL, OPTIONAL :: local

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 2)
#else
      arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:)=
     & arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:)
#endif

      RETURN
      END SUBROUTINE UNPACK_4D

      SUBROUTINE UNPACK_5D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:,:,:,:)
      INTEGER :: J_0, J_1, L,M
      LOGICAL, OPTIONAL :: local

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 2)
#else
      arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)=
     & arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)
#endif

      RETURN
      END SUBROUTINE UNPACK_5D

      SUBROUTINE UNPACKj_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      INTEGER :: J_0, J_1
      LOGICAL, OPTIONAL :: local

      LOGICAL :: local_
      INTEGER :: k

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP,:)=
     & arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      RETURN
      END SUBROUTINE UNPACKj_2D

      SUBROUTINE UNPACKj_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      LOGICAL, OPTIONAL :: local

      INTEGER :: k

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:)=
     &     arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:,:)
#endif

      RETURN
      END SUBROUTINE UNPACKj_3D

      SUBROUTINE UNPACKj_4D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:,:,:)
      LOGICAL, OPTIONAL :: local

      INTEGER :: k

#ifdef USE_MPI
      Call scatter(grd_dum, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)=
     &     arr_glob(grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)
#endif

      RETURN
      END SUBROUTINE UNPACKj_4D

C*** Overloaded PACK_COLUMN routines:
C------------------------------------
      SUBROUTINE PACK_COLUMN_1D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(:,grd_dum%j_strt_halo:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:)
      integer :: K

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 2)
#else
      arr_glob(:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     arr(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE PACK_COLUMN_1D


      SUBROUTINE PACK_COLUMN_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(:,:,grd_dum%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)
      INTEGER :: K

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 3)
#else
      arr_glob(:,:,grd_dum%J_STRT:grd_dum%J_STOP)=
     &     arr(:,:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE PACK_COLUMN_2D

      SUBROUTINE PACK_COLUMN_i2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(:,:,grd_dum%j_strt_halo:)
      INTEGER, INTENT(INOUT) :: ARR_GLOB(:,:,:)
      INTEGER :: I, K

      DO K=1,SIZE(ARR,1)
        DO I = 1, SIZE(ARR,2)
          CALL ARRAYGATHER(grd_dum,ARR(K,I,:),ARR_GLOB(K,I,:))
        END DO
      END DO

      RETURN
      END SUBROUTINE PACK_COLUMN_i2D

      SUBROUTINE PACK_COLUMN_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(:,:,grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: l,K

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 3)
#else
      arr_glob(:,:,grd_dum%J_STRT:grd_dum%J_STOP,:)=
     & arr(:,:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      RETURN
      END SUBROUTINE PACK_COLUMN_3D

      SUBROUTINE PACK_BLOCK_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(:,:,:,:,grd_dum%j_strt_halo:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:,:,:,:)
      INTEGER :: l,K

#ifdef USE_MPI
      CALL gather(grd_dum, arr, arr_glob, shape(arr), 5)
      !mkb5
#else
      arr_glob(:,:,:,:,grd_dum%J_STRT:grd_dum%J_STOP)=
     & arr(:,:,:,:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      RETURN
      END SUBROUTINE PACK_BLOCK_3D

C*** Overloaded UNPACK_COLUMN routines:
C--------------------------------------
      SUBROUTINE UNPACK_COLUMN_1D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1, K
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,J_0:J_1)=ARR_GLOB(:,J_0:J_1)
        else
          do k=1,size(arr,1)
            call arrayscatter(grd_dum, arr(k,:), arr_glob(k,:))
          end do
        end if
      else
        do k=1,size(arr,1)
          call arrayscatter(grd_dum, arr(k,:), arr_glob(k,:))
        end do
      end if

      RETURN
      END SUBROUTINE UNPACK_COLUMN_1D

      SUBROUTINE UNPACK_COLUMN_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1
      LOGICAL, OPTIONAL, intent(in) :: local

      logical :: local_


#ifdef USE_MPI
      !mkb5
      local_ = .false.
      if (present(local)) local_ = local
#else
      local_ = .true.
#endif

      if (local_)then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,J_0:J_1)=ARR_GLOB(:,:,J_0:J_1)
#ifdef USE_MPI
      else
         call scatter(grd_dum , arr_glob, arr,
     &        shape(arr), dist_idx = 3)
#endif
      end if

      RETURN
      END SUBROUTINE UNPACK_COLUMN_2D

c      SUBROUTINE UNPACK_COLUMN_2D(grd_dum,ARR_GLOB,ARR,local)
c      IMPLICIT NONE
c      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
c
c      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
c      REAL*8, INTENT(OUT) ::
c     &        ARR(:,:,grd_dum%j_strt_halo:)
c      INTEGER :: J_0, J_1, K
c      LOGICAL, OPTIONAL :: local
c
c      if (present(local)) then
c        if (local) then
c          J_0=grd_dum%j_strt
c          J_1=grd_dum%j_stop
c          ARR(:,:,J_0:J_1)=ARR_GLOB(:,:,J_0:J_1)
c        else
c          do k=1,size(arr,1)
c            call arrayscatter(grd_dum, arr(k,:,:), arr_glob(k,:,:))
c          end do
c        end if
c      else
c        do k=1,size(arr,1)
c          call arrayscatter(grd_dum, arr(k,:,:), arr_glob(k,:,:))
c        end do
c      end if
c
c      RETURN
c      END SUBROUTINE UNPACK_COLUMN_2D

      SUBROUTINE IUNPACK_COLUMN_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) :: ARR_GLOB(:,:,:)
      INTEGER, INTENT(OUT) ::
     &        ARR(:,:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1, K
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,J_0:J_1)=ARR_GLOB(:,:,J_0:J_1)
        else
          do k=1,size(arr,1)
            call arrayscatter(grd_dum, arr(k,:,:), arr_glob(k,:,:))
          end do
        end if
      else
        do k=1,size(arr,1)
          call arrayscatter(grd_dum, arr(k,:,:), arr_glob(k,:,:))
        end do
      end if

      RETURN
      END SUBROUTINE IUNPACK_COLUMN_2D

      SUBROUTINE UNPACK_COLUMN_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0, J_1, K, L
      LOGICAL, OPTIONAL, intent(in) :: local

      logical :: local_


#ifdef USE_MPI
      local_ = .false.
      if (present(local)) local_ = local
#else
      local_ = .true.
#endif

      if (local_)then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,J_0:J_1,:)=ARR_GLOB(:,:,J_0:J_1,:)
#ifdef USE_MPI
      else
         call scatter(grd_dum , arr_glob, arr,
     &        shape(arr), dist_idx = 3)
#endif
      end if

      RETURN
      END SUBROUTINE UNPACK_COLUMN_3D

      SUBROUTINE UNPACK_BLOCK_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,:,:,:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1, K, L
      LOGICAL, OPTIONAL, intent(in) :: local

      logical :: local_


#ifdef USE_MPI
      !mkb5
      local_ = .false.
      if (present(local)) local_ = local
#else
      local_ = .true.
#endif

      if (local_)then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,:,:,J_0:J_1)=ARR_GLOB(:,:,:,:,J_0:J_1)
#ifdef USE_MPI
      else
         call scatter(grd_dum , arr_glob, arr,
     &        shape(arr), dist_idx = 5)
#endif
      end if

      RETURN
      END SUBROUTINE UNPACK_BLOCK_3D


      SUBROUTINE IPACK_BLOCK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(:,:,:,grd_dum%j_strt_halo:)
      INTEGER, INTENT(INOUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: K,L

      DO K=1,SIZE(ARR,1)
        DO L=1,SIZE(ARR,2)
          CALL ARRAYGATHER(grd_dum,ARR(K,L,:,:),ARR_GLOB(K,L,:,:))
        END DO
      END DO

      RETURN
      END SUBROUTINE IPACK_BLOCK_2D

      SUBROUTINE PACK_BLOCK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8 , INTENT(IN) ::
     &        ARR(:,:,:,grd_dum%j_strt_halo:)
      REAL*8 , INTENT(INOUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: K,L

#ifdef USE_MPI
      Call Gather(grd_dum, arr, arr_glob, shape(arr), 4)
#else
      arr_glob(:,:,:,grd_dum%J_STRT:grd_dum%J_STOP) =
     &     arr(:,:,:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif
c      DO K=1,SIZE(ARR,1)
c        DO L=1,SIZE(ARR,2)
c          CALL ARRAYGATHER(grd_dum,ARR(K,L,:,:),ARR_GLOB(K,L,:,:))
c        END DO
c      END DO

      RETURN
      END SUBROUTINE PACK_BLOCK_2D


c      SUBROUTINE PACK_BLOCK_3D(grd_dum,ARR,ARR_GLOB)  mkbhat5
c      IMPLICIT NONE
c      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
c
c      REAL*8 , INTENT(IN) ::
c     &        ARR(:,:,:,grd_dum%j_strt_halo:,:)
c      REAL*8 , INTENT(INOUT) :: ARR_GLOB(:,:,:,:,:)
c      INTEGER :: K,L,M
c
c      DO M=1,SIZE(ARR,5)
c        DO L=1,SIZE(ARR,2)
c          DO K=1,SIZE(ARR,1)
c            CALL ARRAYGATHER(grd_dum,ARR(K,L,:,:,M),ARR_GLOB(K,L,:,:,M))
c          END DO
c        END DO
c      END DO
c
c      RETURN
c      END SUBROUTINE PACK_BLOCK_3D


      SUBROUTINE IUNPACK_BLOCK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      INTEGER, INTENT(OUT) ::
     &        ARR(:,:,:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1, K, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,:,J_0:J_1)=ARR_GLOB(:,:,:,J_0:J_1)
        else
          do k=1,size(arr,1)
            do l=1,size(arr,2)
              call arrayscatter(grd_dum,arr(k,l,:,:),arr_glob(k,l,:,:))
            end do
          end do
        end if
      else
        do k=1,size(arr,1)
          do l=1,size(arr,2)
            call arrayscatter(grd_dum, arr(k,l,:,:), arr_glob(k,l,:,:))
          end do
        end do
      end if
      END SUBROUTINE IUNPACK_BLOCK_2D

      SUBROUTINE UNPACK_BLOCK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,:,:,grd_dum%j_strt_halo:)
      INTEGER :: J_0, J_1
      LOGICAL, OPTIONAL, intent(in) :: local

      logical :: local_


#ifdef USE_MPI
      !mkb5
      local_ = .false.
      if (present(local)) local_ = local
#else
      local_ = .true.
#endif

      if (local_)then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,:,J_0:J_1)=ARR_GLOB(:,:,:,J_0:J_1)
#ifdef USE_MPI
      else
         call scatter(grd_dum , arr_glob, arr,
     &        shape(arr), dist_idx = 4)
#endif
      end if

      RETURN
      END SUBROUTINE UNPACK_BLOCK_2D

c      SUBROUTINE UNPACK_BLOCK_2D(grd_dum,ARR_GLOB,ARR,local)
c      IMPLICIT NONE
c      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
c
c      REAL*8 , INTENT(IN) :: ARR_GLOB(:,:,:,:)
c      REAL*8 , INTENT(OUT) ::
c     &        ARR(:,:,:,grd_dum%j_strt_halo:)
c      INTEGER :: J_0, J_1, K, L
c      LOGICAL, OPTIONAL :: local
c
c      if (present(local)) then
c        if (local) then
c          J_0=grd_dum%j_strt
c          J_1=grd_dum%j_stop
c          ARR(:,:,:,J_0:J_1)=ARR_GLOB(:,:,:,J_0:J_1)
c        else
c          do k=1,size(arr,1)
c            do l=1,size(arr,2)
c              call arrayscatter(grd_dum,arr(k,l,:,:),arr_glob(k,l,:,:))
c            end do
c          end do
c        end if
c      else
c        do k=1,size(arr,1)
c          do l=1,size(arr,2)
c            call arrayscatter(grd_dum, arr(k,l,:,:), arr_glob(k,l,:,:))
c          end do
c        end do
c      end if
c      END SUBROUTINE UNPACK_BLOCK_2D

c      SUBROUTINE UNPACK_BLOCK_3D(grd_dum,ARR_GLOB,ARR,local)   mkbhat5
c      IMPLICIT NONE
c      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
c
c      REAL*8 , INTENT(IN) :: ARR_GLOB(:,:,:,:,:)
c      REAL*8 , INTENT(OUT) ::
c     &        ARR(:,:,:,grd_dum%j_strt_halo:,:)
c      INTEGER :: J_0, J_1, K, L, M
c      LOGICAL, OPTIONAL :: local
c
c      if (present(local)) then
c        if (local) then
c          J_0=grd_dum%j_strt
c          J_1=grd_dum%j_stop
c          ARR(:,:,:,J_0:J_1,:)=ARR_GLOB(:,:,:,J_0:J_1,:)
c        else
c          do m=1,size(arr,5)
c            do l=1,size(arr,2)
c              do k=1,size(arr,1)
c                call arrayscatter(grd_dum,arr(k,l,:,:,m)
c     &                              ,arr_glob(k,l,:,:,m))
c              end do
c            end do
c          end do
c        end if
c      else
c        do m=1,size(arr,5)
c          do l=1,size(arr,2)
c            do k=1,size(arr,1)
c              call arrayscatter(grd_dum, arr(k,l,:,:,m)
c     &                            , arr_glob(k,l,:,:,m))
c            end do
c          end do
c        end do
c      end if
c      END SUBROUTINE UNPACK_BLOCK_3D


C*** Overloaded PACK_J routines:
C-------------------------------
      SUBROUTINE PACK_J_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:)
      integer :: K

      DO K=1,SIZE(ARR,2)
       CALL ARRAYGATHER(grd_dum,ARR(:,k),ARR_GLOB(:,k) )
      END DO

      RETURN
      END SUBROUTINE PACK_J_2D


      SUBROUTINE PACK_J_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)
      INTEGER :: K,L

      DO L=1,SIZE(ARR,3)
        DO K=1,SIZE(ARR,2)
          CALL ARRAYGATHER(grd_dum,ARR(:,K,L),ARR_GLOB(:,K,L))
        END DO
      END DO

      RETURN
      END SUBROUTINE PACK_J_3D

      SUBROUTINE PACK_J_4D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: K,L,M

      DO M=1,SIZE(ARR,4)
        DO L=1,SIZE(ARR,3)
          DO K=1,SIZE(ARR,2)
            CALL ARRAYGATHER( grd_dum,ARR(:,K,L,M), ARR_GLOB(:,K,L,M) )
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE PACK_J_4D


C*** Overloaded UNPACK_J routines:
C--------------------------------
      SUBROUTINE UNPACK_J_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      INTEGER :: J_0H, J_1H, K
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=Max(1,grd_dum%j_strt_halo)
          J_1H=Min(grd_dum%JM_WORLD,grd_dum%j_stop_halo)
          ARR(J_0H:J_1H,:)=ARR_GLOB(J_0H:J_1H,:)
        else
          do k=1,size(arr,2)
            call arrayscatter(grd_dum, arr(:,k), arr_glob(:,k))
          end do
        end if
      else
        do k=1,size(arr,2)
          call arrayscatter(grd_dum, arr(:,k), arr_glob(:,k))
        end do
      end if

      RETURN
      END SUBROUTINE UNPACK_J_2D


      SUBROUTINE UNPACK_J_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      INTEGER :: J_0H, J_1H, K,L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=Max(1,grd_dum%j_strt_halo)
          J_1H=Min(grd_dum%JM_WORLD,grd_dum%j_stop_halo)
          ARR(J_0H:J_1H,:,:)=ARR_GLOB(J_0H:J_1H,:,:)
        else
          do l=1,size(arr,3)
            do k=1,size(arr,2)
              call arrayscatter(grd_dum, arr(:,k,l), arr_glob(:,k,l))
            end do
          end do
        end if
      else
          do l=1,size(arr,3)
            do k=1,size(arr,2)
              call arrayscatter(grd_dum, arr(:,k,l), arr_glob(:,k,l))
            end do
          end do
      end if

      RETURN
      END SUBROUTINE UNPACK_J_3D


      SUBROUTINE UNPACK_J_4D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:,:,:)
      INTEGER :: J_0H, J_1H, K, L, M
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(J_0H:J_1H,:,:,:)=ARR_GLOB(J_0H:J_1H,:,:,:)
        else
          do m=1,size(arr,4)
            do l=1,size(arr,3)
              do k=1,size(arr,2)
                call arrayscatter( grd_dum, arr(:,k,l,m),
     &                                 arr_glob(:,k,l,m) )
              end do
            end do
          end do
        end if
      else
        do m=1,size(arr,4)
          do l=1,size(arr,3)
            do k=1,size(arr,2)
              call arrayscatter( grd_dum, arr(:,k,l,m),
     &                               arr_glob(:,k,l,m) )
            end do
          end do
        end do
      end if

      RETURN
      END SUBROUTINE UNPACK_J_4D

      SUBROUTINE INIT_BAND_PACK_TYPE(grd_dum, band_j0,band_j1, bandpack)
c initialize the bandpack derived type with the information needed
c for the band_pack procedure to fill output arrays with data from
c J indices band_j0 to band_j1
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER, INTENT(IN) :: band_j0,band_j1
      TYPE (BAND_PACK_TYPE), intent(OUT) :: bandpack
#ifdef USE_MPI
      integer, dimension(0:npes-1) ::
     &     j0_have,j1_have,j0_requested,j1_requested
      integer :: p, ierr, im,jm, j0send,j1send,j0recv,j1recv

      allocate(bandpack%j0_send(0:npes-1), bandpack%j1_send(0:npes-1))
      allocate(bandpack%j0_recv(0:npes-1), bandpack%j1_recv(0:npes-1))
      allocate(bandpack%scnts(0:npes-1), bandpack%sdspl(0:npes-1))
      allocate(bandpack%rcnts(0:npes-1), bandpack%rdspl(0:npes-1))
      allocate(bandpack%sdspl_inplace(0:npes-1))
      allocate(bandpack%rdspl_inplace(0:npes-1))
#endif
      bandpack%im_world = grd_dum%im_world
      bandpack%j_strt = grd_dum%j_strt
      bandpack%j_stop = grd_dum%j_stop
      bandpack%j_strt_halo = grd_dum%j_strt_halo
      bandpack%j_stop_halo = grd_dum%j_stop_halo
      bandpack%jband_strt = band_j0
      bandpack%jband_stop = band_j1
#ifdef USE_MPI
      im = grd_dum%im_world
      jm = grd_dum%jm_world
c
c Set up the MPI send/receive information
c
      call mpi_allgather(grd_dum%j_strt,1,MPI_INTEGER,j0_have,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(grd_dum%j_stop,1,MPI_INTEGER,j1_have,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(band_j0,1,MPI_INTEGER,j0_requested,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(band_j1,1,MPI_INTEGER,j1_requested,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      do p=0,npes-1
        j0send = max(grd_dum%j_strt,j0_requested(p))
        j1send = min(grd_dum%j_stop,j1_requested(p))
        bandpack%j0_send(p) = j0send
        bandpack%j1_send(p) = j1send
        if(j0send <= j1send) then
          bandpack%scnts(p) = im*(j1send-j0send+1)
          bandpack%sdspl_inplace(p) = im*(j0send-grd_dum%j_strt_halo)
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

      SUBROUTINE BAND_PACK_ij(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, dimension(0:npes-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr
      scnts = bandpack%scnts
      sdspl = bandpack%sdspl_inplace
      rcnts = bandpack%rcnts
      rdspl = bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision,
     &                   arr_band, rcnts, rdspl, mpi_double_precision,
     &                   mpi_comm_world, ierr)
#else
      arr_band(:,bandpack%J_STRT:bandpack%J_STOP) =
     &     arr(:,bandpack%J_STRT:bandpack%J_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ij

      SUBROUTINE BAND_PACK_ijl(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,bandpack%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,bandpack%jband_strt:,:)
#ifdef USE_MPI
      integer, dimension(0:npes-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,lm,i,j,l,n,p
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      lm = size(arr,3)
      scnts = lm*bandpack%scnts
      sdspl = lm*bandpack%sdspl
      rcnts = lm*bandpack%rcnts
      rdspl = lm*bandpack%rdspl
      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_send(p),bandpack%j1_send(p)
            do i=1,bandpack%im_world
              n = n + 1
              bufsend(n) = arr(i,j,l)
            enddo
          enddo
        enddo
      enddo
      call mpi_alltoallv(bufsend, scnts, sdspl, mpi_double_precision,
     &                   bufrecv, rcnts, rdspl, mpi_double_precision,
     &                   mpi_comm_world, ierr)
      n = 0
      do p=0,npes-1
        do l=1,lm
          do j=bandpack%j0_recv(p),bandpack%j1_recv(p)
            do i=1,bandpack%im_world
              n = n + 1
              arr_band(i,j,l) = bufrecv(n)
            enddo
          enddo
        enddo
      enddo
      deallocate(bufsend,bufrecv)
#else
      arr_band(:,bandpack%J_STRT:bandpack%J_STOP,:) =
     &     arr(:,bandpack%J_STRT:bandpack%J_STOP,:)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_ijl

      SUBROUTINE BAND_PACK_COLUMN(bandpack,ARR,ARR_band)
!@var bandpack (input) instance of the band_pack_type structure
!@var ARR      (input) local domain-decomposed array on this PE
!@var ARR_band (output) array dimensioned and filled over the
!@+   J range requested during the call to
!@+   init_band_pack_type that initialized bandpack
      IMPLICIT NONE
      TYPE (BAND_PACK_TYPE),  INTENT(IN) :: bandpack
      REAL*8, INTENT(IN) :: ARR(:,:,bandpack%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_band(:,:,bandpack%jband_strt:)
#ifdef USE_MPI
      integer, dimension(0:npes-1) :: scnts,sdspl, rcnts,rdspl
      integer :: ierr,lm
      lm = size(arr,1)
      scnts = lm*bandpack%scnts
      sdspl = lm*bandpack%sdspl_inplace
      rcnts = lm*bandpack%rcnts
      rdspl = lm*bandpack%rdspl_inplace
      call mpi_alltoallv(arr, scnts, sdspl, mpi_double_precision,
     &                   arr_band, rcnts, rdspl, mpi_double_precision,
     &                   mpi_comm_world, ierr)
#else
      arr_band(:,:,bandpack%J_STRT:bandpack%J_STOP) =
     &     arr(:,:,bandpack%J_STRT:bandpack%J_STOP)
#endif
      RETURN
      END SUBROUTINE BAND_PACK_COLUMN

      SUBROUTINE ESMF_BCAST_0D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ierr)
#endif

      END SUBROUTINE ESMF_BCAST_0D

      SUBROUTINE ESMF_BCAST_1D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ierr)
#endif

      END SUBROUTINE ESMF_BCAST_1D

      SUBROUTINE ESMF_BCAST_2D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ierr)
#endif

      END SUBROUTINE ESMF_BCAST_2D

      SUBROUTINE ESMF_BCAST_3D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ierr)
#endif

      END SUBROUTINE ESMF_BCAST_3D

      SUBROUTINE ESMF_BCAST_4D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:,:,:)

      INTEGER :: ierr

#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ierr)
#endif

      END SUBROUTINE ESMF_BCAST_4D

      SUBROUTINE ESMF_IBCAST_0D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,1,MPI_INTEGER,root,
     &     MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_0D

      SUBROUTINE ESMF_IBCAST_1D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root,
     &     MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_1D

      SUBROUTINE ESMF_IBCAST_2D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root,
     &     MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_2D

      SUBROUTINE ESMF_IBCAST_3D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root,
     &     MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_3D

      SUBROUTINE ESMF_IBCAST_4D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:,:,:)
      INTEGER :: ierr
#ifdef USE_MPI
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root,
     &     MPI_COMM_WORLD, ierr)
#endif
      END SUBROUTINE ESMF_IBCAST_4D


      SUBROUTINE TRANSPOSE_ijk(grid, x_in, x_out, reverse)
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_out(:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES-1), I1(0:NPES-1)
      INTEGER :: J0(0:NPES-1), J1(0:NPES-1)
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,nk,k
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt
      INTEGER, DIMENSION(0:NPES-1) :: scnts, rcnts, sdspl, rdspl
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
      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes-1,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

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
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ierr)
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

      SUBROUTINE TRANSPOSE_ij(grid, x_in, x_out, reverse)
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x_in(:,grid%J_STRT_HALO:)
      REAL*8 :: x_out(:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES-1), I1(0:NPES-1)
      INTEGER :: J0(0:NPES-1), J1(0:NPES-1)
      REAL*8, ALLOCATABLE :: sbuf(:), rbuf(:)
#ifdef USE_MPI
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt
      INTEGER, DIMENSION(0:NPES-1) :: scnts, rcnts, sdspl, rdspl
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
      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes-1,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

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
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ierr)
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




      SUBROUTINE TRANSPOSE_COLUMN(grid, x, x_tr, reverse)
      TYPE (DIST_GRID), INTENT(IN) :: grid
      REAL*8 :: x(:,:,grid%J_STRT_HALO:,:)
      REAL*8 :: x_tr(:,:,:,:)
      Logical, Optional, INTENT(IN) :: reverse

      INTEGER :: I0(0:NPES-1), I1(0:NPES-1)
      INTEGER :: J0(0:NPES-1), J1(0:NPES-1)
      REAL*8, ALLOCATABLE :: sbuf(:,:), rbuf(:,:)
#ifdef USE_MPI
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,k
      INTEGER :: ierr, p, rc
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt
      INTEGER, DIMENSION(0:NPES-1) :: scnts, rcnts, sdspl, rdspl
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

      DO p = 0, npes - 1
         I0(p) = 1 + p * grid%IM_WORLD / NPES
         I1(p) = (p+1) * grid%IM_WORLD / NPES
      END DO

      ALLOCATE(AI(0:npes-1,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

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
         CALL MPI_ALLTOALLV(rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ierr)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ierr)
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

#ifdef USE_MPI
      Function CreateDist_MPI_Type(base_type, counts, dist_idx)
     &     Result(new_type)
      Integer, Intent(In) :: base_type
      Integer, Intent(In) :: counts(:)
      Integer, Intent(In) :: dist_idx
      Integer :: new_type

      Integer :: stride, n_blocks, blocklen
      Integer :: new_len
      INTEGER(KIND=MPI_ADDRESS_KIND) disp(2), ext_lb, base_byte_len
      Integer :: vtype1, vtype2
      Integer :: ierr

#ifdef MPITYPE_LOOKUP_HACK
      type mpi_types_str
      integer n_blocks,blocklen,stride,base_type,vtype
      end type mpi_types_str

      type(mpi_types_str), save :: mt(1024)
      integer, save :: mt_count = 0
      integer m


      !write(444,*) base_type,dist_idx,counts

#endif

      n_blocks = Product(counts(dist_idx+1:))
      blocklen = Product(counts(:dist_idx-1))
      stride = counts(dist_idx) * blocklen
      ext_lb = 0

#ifdef MPITYPE_LOOKUP_HACK
      do m=1,mt_count
        if ( n_blocks == mt(m)%n_blocks .and.
     &       blocklen == mt(m)%blocklen .and.
     &       stride   == mt(m)%stride   .and.
     &       base_type   == mt(m)%base_type ) then
          new_type = mt(m)%vtype 
          return
        endif
      enddo
#endif

      Call MPI_Type_vector(n_blocks, blocklen, stride, base_type,
     &     vtype1, ierr)
      !Call MPI_Type_extent(base_type, base_byte_len, ierr)
      Call MPI_Type_get_extent(base_type, ext_lb, base_byte_len, ierr)
      new_len = base_byte_len * blocklen
      !Call MPI_Type_struct(2, (/ 1, 1 /), (/ 0, new_len /),
      disp(1)=0; disp(2)=new_len;
      Call MPI_Type_create_struct(2, (/ 1, 1 /), disp,
     &     (/ vtype1, MPI_UB /), vtype2, ierr)
      Call MPI_Type_Free(vtype1, ierr)

      Call MPI_Type_Commit(vtype2, ierr)

#ifdef MPITYPE_LOOKUP_HACK
      mt_count = mt_count + 1
      print *,"NEW MPI_Type: ", mt_count
      if ( mt_count > 1024 ) call stop_model("mt_count > 1024",255)
      mt(mt_count)%n_blocks = n_blocks
      mt(mt_count)%blocklen = blocklen
      mt(mt_count)%stride = stride
      mt(mt_count)%base_type = base_type
      mt(mt_count)%vtype = vtype2
#endif

      new_type = vtype2

      End Function CreateDist_MPI_Type

      Subroutine GetNeighbors(rank, npes, pe_south,pe_north,bc_periodic)
      Integer, Intent(In)  :: rank
      Integer, Intent(In)  :: npes
      Integer, Intent(Out) :: pe_south
      Integer, Intent(Out) :: pe_north
      Logical, Intent(In)  :: bc_periodic

cddd      If (rank > 0) Then
cddd        pe_south = rank - 1
cddd      Else
cddd        pe_south = MPI_PROC_NULL
cddd      End If
cddd
cddd      If (rank < npes-1) Then
cddd        pe_north = rank + 1
cddd      Else
cddd        pe_north = MPI_PROC_NULL
cddd      End If

      pe_south = mod( npes + rank - 1, npes)
      pe_north = mod(        rank + 1, npes)

      if ( .not. bc_periodic ) then
        if ( rank == 0      ) pe_south = MPI_PROC_NULL
        if ( rank == npes-1 ) pe_north = MPI_PROC_NULL
      endif

      End Subroutine GetNeighbors

      Subroutine sendrecv(grid, arr, shp, dist_idx, from, bc_periodic_)
      Type (Esmf_Grid) :: grid
      Real(Kind=8) :: arr(*)
      Integer :: shp(:)
      Integer :: dist_idx
      Integer, optional :: from
      Logical, optional :: bc_periodic_

      Integer :: new_type
      Integer :: npy, npx, px, py, pe_south, pe_north
      Integer :: off(4)
      Integer :: USABLE_FROM
      Integer :: status(MPI_STATUS_SIZE), ierr
      Integer :: n, sz
      Logical :: bc_periodic

      integer :: requests(2)
      integer :: tagS, tagN
      integer :: nSendMessages, nRecvMessages
      integer :: index, i

#ifdef USE_SYSUSAGE
      call sysusage(3,1)
#endif
 
      USABLE_FROM = usableFrom(from)
      bc_periodic = isPeriodic(bc_periodic_)

      ! create a new mpi type for use in communication
      !-------------------------------
      new_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION, shp,dist_idx)

      ! Determine neigboring processes
      !-------------------------------
      call ESMF_GRID_MY_PE_LOC(grid,  px,  py)
      call ESMF_GRID_PE_LAYOUT(grid, npx, npy)
      Call GetNeighbors(py, npy, pe_south, pe_north, bc_periodic)

      sz = Product(shp(1:dist_idx-1))
      n  = shp(dist_idx)
      off(1) = 1
      off(2) = 1+sz*1
      off(3) = 1+sz*(n-2)
      off(4) = 1+sz*(n-1)

#ifdef USE_SYSUSAGE
      call sysusage(4,1)
#endif
      nSendMessages = 0
      nRecvMessages = 0
      tag = max(MOD(tag,128),10) + 1
      IF(FILL(NORTH)) THEN
        if (pe_north /= MPI_PROC_NULL) nRecvMessages = nRecvMessages + 1
        if (pe_south /= MPI_PROC_NULL) then
           nSendMessages = nSendMessages + 1
           call MPI_Isend(arr(off(2)), 1, new_type, pe_south, tag, 
     &          MPI_COMM_WORLD, requests(nSendMessages), ierr)
        end if
      end if

      IF(FILL(SOUTH)) THEN
        if (pe_south /= MPI_PROC_NULL) nRecvMessages = nRecvMessages + 1
        if (pe_north /= MPI_PROC_NULL) then
           nSendMessages = nSendMessages + 1
           call MPI_Isend(arr(off(3)), 1, new_type, pe_north, tag, 
     &          MPI_COMM_WORLD, requests(nSendMessages), ierr)
        end if
      end if

      do i = 1, nRecvMessages
         call MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, 
     &        status, ierr)
         if (status(MPI_SOURCE) == pe_north) then
            call MPI_Recv( arr(off(4)), 1, new_type, pe_north, tag,
     &           MPI_COMM_WORLD, status, ierr)
         else
            call MPI_Recv( arr(off(1)), 1, new_type, pe_south, tag,
     &           MPI_COMM_WORLD, status, ierr)
         end If
      end do

      call MPI_WaitAll(nSendMessages, requests,MPI_STATUSES_IGNORE,ierr)


#ifdef USE_SYSUSAGE
      call sysusage(4,2)
#endif

#ifndef MPITYPE_LOOKUP_HACK
      Call MPI_Type_Free(new_type, ierr)
#endif
#ifdef USE_SYSUSAGE
      call sysusage(3,2)
#endif

      End SUBROUTINE SendRecv

      Subroutine sendrecv_int(grid, arr, shp, dist_idx, from,
     &     bc_periodic_)
      Type (Esmf_Grid) :: grid
      Integer :: arr(*)
      Integer :: shp(:)
      Integer :: dist_idx
      Integer, optional :: from
      Logical, optional :: bc_periodic_

      Integer :: new_type
      Integer :: npy, npx, px, py, pe_south, pe_north
      Integer :: off(4)
      Integer :: USABLE_FROM
      Integer :: status(MPI_STATUS_SIZE), ierr
      Integer :: n, sz
      Logical :: bc_periodic

      USABLE_FROM = usableFrom(from)
      bc_periodic = isPeriodic(bc_periodic_)
      ! create a new mpi type for use in communication
      !-------------------------------
      new_type = CreateDist_MPI_Type(MPI_INTEGER, shp,dist_idx)

      ! Determine neigboring processes
      !-------------------------------
      call ESMF_GRID_MY_PE_LOC(grid,  px,  py)
      call ESMF_GRID_PE_LAYOUT(grid, npx, npy)
      Call GetNeighbors(py, npy, pe_south, pe_north, bc_periodic)

      sz = Product(shp(1:dist_idx-1))
      n  = shp(dist_idx)
      off(1) = 1
      off(2) = 1+sz*1
      off(3) = 1+sz*(n-2)
      off(4) = 1+sz*(n-1)

      IF(FILL(NORTH)) THEN
        tag = max(MOD(tag,128),10) + 1
        Call MPI_SendRecv(arr(off(2)), 1, new_type, pe_south, tag,
     &                    arr(off(4)), 1, new_type, pe_north, tag,
     &                    MPI_COMM_WORLD, status, ierr)
      End If

      IF(FILL(SOUTH)) THEN
        tag = max(MOD(tag,128),10) + 1
        Call MPI_SendRecv(arr(off(3)), 1, new_type, pe_north, tag,
     &                    arr(off(1)), 1, new_type, pe_south, tag,
     &                    MPI_COMM_WORLD, status, ierr)
      End If

#ifndef MPITYPE_LOOKUP_HACK
      Call MPI_Type_Free(new_type, ierr)
#endif

      End SUBROUTINE SendRecv_int

      Subroutine gather(grid, arr_loc, arr_glob, shp, dist_idx, all)
      !Type (Esmf_Grid) :: grid
      Type (DIST_GRID) :: grid
      Real(Kind=8) :: arr_loc(*)
      Real(Kind=8) :: arr_glob(*)
      Integer :: shp(:),shp_glob(size(shp))
      Integer :: dist_idx
      Logical, Optional :: all

      Logical :: all_
      Integer :: new_type, orig_type
      Integer :: ierr, rc
      Integer :: p, scount, offset
      Integer, Allocatable :: rcounts(:), displs(:)
      Type (ESMF_Axisindex), Pointer :: AI(:,:)
      Integer :: i, n_its = 1
      ! create a new mpi type for use in communication
      !-------------------------------
      orig_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION,shp,dist_idx)


      Allocate(AI(NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

      allocate (rcounts(NPES), displs(NPES), stat=rc)
      Do p = 1, npes
        rcounts(p) = AI(p,2)%max - AI(p,2)%min + 1
        displs(p) = AI(p,2)%min - 1
      End Do

      shp_glob = shp
      shp_glob(dist_idx) = AI(npes,2)%max-AI(1,2)%min + 1
      new_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION,
     &     shp_glob,dist_idx)
      Deallocate(AI)

      scount = rcounts(my_pet+1)
      If (scount == shp(dist_idx)) Then ! no halo
        offset = 1
      Else
        offset = 1 + Product(shp(:dist_idx-1))
      End If

      all_=.false.
      If (Present(all)) all_=all

      Do i = 1, n_its
      If (all_) Then
        Call MPI_AllGatherV(arr_loc(offset), scount, orig_type,
     &       arr_glob(1), rcounts, displs, new_type,
     &       MPI_COMM_WORLD, ierr)
      Else
        Call MPI_GatherV(arr_loc(offset), scount, orig_type,
     &       arr_glob(1), rcounts, displs, new_type,
     &       root, MPI_COMM_WORLD, ierr)
      End If

      End Do

#ifndef MPITYPE_LOOKUP_HACK
      Call MPI_Type_Free(new_type, ierr)
#endif

      Deallocate(rcounts,displs)
#ifndef MPITYPE_LOOKUP_HACK
      Call MPI_Type_Free(orig_type, ierr)
#endif

      End SUBROUTINE gather

      Subroutine scatter(grid, arr_glob, arr_loc, shp, dist_idx)
      !Type (Esmf_Grid) :: grid
      Type (DIST_GRID) :: grid
      Real(Kind=8) :: arr_loc(*)
      Real(Kind=8) :: arr_glob(*)
      Integer :: shp(:),shp_glob(size(shp))
      Integer :: dist_idx

      Integer :: new_type, orig_type
      Integer :: rc, ierr
      Integer :: p, rcount, offset
      Integer, Allocatable :: scounts(:), displs(:)
      Type (ESMF_Axisindex), Pointer :: AI(:,:)

      ! create a new mpi type for use in communication
      !-------------------------------
      new_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION,shp,dist_idx)


      Allocate(AI(NPES,3))
#ifdef USE_MPP
      Call GridGetAllAxisIndex(grid%domain, AI)
#else
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &   horzRelLoc=ESMF_CELL_CENTER, vertRelLoc=ESMF_CELL_CELL, rc=rc)
#endif

      allocate (scounts(NPES), displs(NPES), stat=rc)
      Do p = 1, npes
        scounts(p) = AI(p,2)%max - AI(p,2)%min + 1
        displs(p) = AI(p,2)%min - 1
      End Do

      shp_glob = shp
      shp_glob(dist_idx) = AI(npes,2)%max-AI(1,2)%min + 1
      orig_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION,
     &     shp_glob,dist_idx)

      rcount = scounts(my_pet+1)
      If (rcount == shp(dist_idx)) Then ! no halo
        offset = 1
      Else
        offset = 1 + Product(shp(:dist_idx-1))
      End If

      Call MPI_ScatterV(arr_glob(1), scounts, displs, orig_type,
     &     arr_loc(offset), rcount, new_type,
     &     root, MPI_COMM_WORLD, ierr)

#ifndef MPITYPE_LOOKUP_HACK
      Call MPI_Type_Free(new_type, ierr)
      Call MPI_Type_Free(orig_type, ierr)
#endif

      Deallocate(scounts,displs)
      Deallocate(AI)

      End SUBROUTINE scatter
#endif

      logical function haveLatitude(grd_dum, j)
      type (DIST_GRID), intent(in) :: grd_dum
      integer, intent(in) :: j

      haveLatitude = (j >= grd_dum%J_STRT .and. j <= grd_dum%J_STOP)

      end function haveLatitude


      subroutine SEND_TO_J_1D(grd_dum, arr, j_dest, tag)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, Size(arr), MPI_DOUBLE_PRECISION,
     &     grd_dum%lookup_pet(j_dest), tag, MPI_COMM_WORLD, ierr)
#endif
      end subroutine SEND_TO_J_1D

      subroutine ISEND_TO_J_0D(grd_dum, arr, j_dest, tag)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_dest, tag
      INTEGER :: ierr
#ifdef USE_MPI
      call MPI_Send(arr, 1, MPI_INTEGER,
     &     grd_dum%lookup_pet(j_dest), tag, MPI_COMM_WORLD, ierr)
#endif
      end subroutine ISEND_TO_J_0D

      subroutine RECV_FROM_J_1D(grd_dum, arr, j_src, tag)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(In) :: arr(:)
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, Size(arr), MPI_DOUBLE_PRECISION,
     &     grd_dum%lookup_pet(j_src), tag, MPI_COMM_WORLD, status, ierr)
#endif
      end subroutine RECV_FROM_J_1D

      subroutine IRECV_FROM_J_0D(grd_dum, arr, j_src, tag)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(In) :: arr
      Integer, Intent(In) :: j_src, tag
#ifdef USE_MPI
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      call MPI_Recv(arr, 1, MPI_INTEGER,
     &     grd_dum%lookup_pet(j_src), tag, MPI_COMM_WORLD, status, ierr)
#endif
      end subroutine IRECV_FROM_J_0D

! Helper function to handle optional arguments related to periodic boundaries
      logical function isPeriodic(override)
        logical, optional, intent(in) :: override

        isPeriodic = .false.
        if (present(override)) isPeriodic = override

      end function isPeriodic

! Helper function to handle optional arguments related to halo directions
      integer function usableFrom(fromDirection)
        integer, optional, intent(in) :: fromDirection
        usableFrom = ALL
        if (present(fromDirection)) usableFrom = fromDirection
      end function usableFrom

#ifdef USE_ESMF

  !----------------------------
      function load_cap_config(config_file,IM,JM,LM,NP_X,NP_Y) 
     &         result( config )
         use ESMF_mod
         use FILEMANAGER    
         character(len=*), parameter :: Iam=
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
  !----------------------------

#endif

#ifdef USE_MPP
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: init_domain --- MPP SPMD parallel decompostion/communication module
!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
!     init_domain :: Setup MPP domains
!
      subroutine init_domain(domain,npx,npy,ntiles,ng,bc_periodic)

      type(domain2D) :: domain
      integer, intent(IN)  :: npx,npy,ntiles,ng
      logical :: bc_periodic

      integer :: npes, grid_type
      integer :: layout(2)
      integer, allocatable :: pe_start(:), pe_end(:)

      integer :: ios,nx,ny,n,num_alloc
      character(len=32) :: type
      integer :: num_contact, npes_per_tile, tile
      integer, allocatable, dimension(:)  :: npes_tile, tile1, tile2
      integer, allocatable, dimension(:)  :: istart1,iend1,jstart1,jend1
      integer, allocatable, dimension(:)  :: istart2,iend2,jstart2,jend2
      integer, allocatable, dimension(:,:) :: layout2D, global_indices
      integer :: is, ie, js, je, isd, ied, jsd, jed

      npes = mpp_npes()
      nx = npx-1
      ny = npy-1

      call mpp_domains_init(MPP_DOMAIN_TIME)

      call mpp_domains_set_stack_size(1500000)

      select case(ntiles)
      case ( 1 )  ! Lat-Lon "cyclic"

       if ( bc_periodic ) then
         grid_type = 4
       else
         grid_type = 7
       end if

       select case (grid_type)
        case (4)   ! Cartesian, double periodic
          num_contact = 2
          npes_per_tile = npes/ntiles
          call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, 
     &      layout )
          layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
        case (7)   ! Cartesian, channel, non-periodic
          num_contact = 1
          npes_per_tile = npes/ntiles
          call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile,
     &     layout )
          layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
       end select

      case ( 6 )  ! Cubed-Sphere

         num_contact = 12
         !--- cubic grid always have six tiles, so npes should be multiple of 6
         if( mod(npes,ntiles) .NE. 0 .OR. npx-1 .NE. npy-1) then
               call stop_model('NPES not divisible by ntiles', 255)
               return
         end if
         npes_per_tile = npes/ntiles
         call  mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, 
     &      layout )

         if ( (npx/layout(1) < ng) .or. (npy/layout(2) < ng) ) then
               write(*,310) layout(1), layout(2),
     &                      npx/layout(1), npy/layout(2)
 310           format('Invalid layout, NPES_X:',i4.4,'NPES_Y:',i4.4,
     &                'ncells_X:',i4.4,'ncells_Y:',i4.4)
               call mpp_error(FATAL, 'init_domain:Invalid layout')
         endif

         layout = (/layout(1),layout(2)/)

      case default

            call mpp_error(FATAL, 'init_domain: no such test: '//type)

      end select


      allocate( layout2D(2,ntiles), global_indices(4,ntiles) )
      allocate( npes_tile(ntiles) )
      allocate(pe_start(ntiles),pe_end(ntiles))
      npes_tile = npes_per_tile
      do n = 1, ntiles
            global_indices(:,n) = (/1,npx-1,1,npy-1/)
            layout2D(:,n)         = layout
            pe_start(n) = (n-1)*layout(1)*layout(2)
            pe_end(n)   = pe_start(n) + layout(1)*layout(2) -1
      end do
      num_alloc=max(1,num_contact)
      allocate( tile1(num_alloc), tile2(num_alloc) )
      allocate( istart1(num_alloc), iend1(num_alloc) )
      allocate( jstart1(num_alloc), jend1(num_alloc) )
      allocate( istart2(num_alloc), iend2(num_alloc) )
      allocate( jstart2(num_alloc), jend2(num_alloc) )

      select case(ntiles)
      case ( 1 )  ! Lat-Lon "cyclic"

       type ='Lat-Lon'
       select case (grid_type)
        case (4)   ! Cartesian, double periodic
          !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
          tile1(1) = 1; tile2(1) = 1
          istart1(1) = nx; iend1(1) = nx;
          jstart1(1) = 1;  jend1(1) = ny;
          istart2(1) = 1;  iend2(1) = 1;
          jstart2(1) = 1;  jend2(1) = ny;
          !--- Contact line 2, between tile 1 (SOUTH) and
          !    tile 1 (NORTH)  --- cyclic
          tile1(2) = 1; tile2(2) = 1
          istart1(2) = 1;  iend1(2) = nx; 
          jstart1(2) = 1;   jend1(2) = 1;
          istart2(2) = 1;  iend2(2) = nx;
          jstart2(2) = ny;  jend2(2) = ny;
          call mpp_define_mosaic(global_indices, layout2D,
     &         domain, ntiles, num_contact, tile1, tile2,
     &         istart1, iend1, jstart1, jend1, istart2, iend2,
     &         jstart2, jend2,
     &         pe_start=pe_start, pe_end=pe_end, symmetry=.true.,
     &      shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, name = type)

         case (7)   ! Cartesian, channel
          !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
          tile1(1) = 1; tile2(1) = 1
          istart1(1) = nx; iend1(1) = nx;
          jstart1(1) = 1;  jend1(1) = ny;
          istart2(1) = 1;  iend2(1) = 1;
          jstart2(1) = 1;  jend2(1) = ny;
          call mpp_define_mosaic(global_indices, layout2D,
     &         domain, ntiles, num_contact, tile1, tile2,
     &         istart1, iend1, jstart1, jend1, istart2, iend2,
     &         jstart2, jend2,
     &         pe_start=pe_start, pe_end=pe_end, symmetry=.true.,
     &      shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, name = type)
         end select

        case ( 6 )  ! Cubed-Sphere
         type="Cubic-Grid"
         !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
         tile1(1) = 1; tile2(1) = 2
         istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
         istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
         !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
         tile1(2) = 1; tile2(2) = 3
         istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
         istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
         !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
         tile1(3) = 1; tile2(3) = 5
         istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
         istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
         !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
         tile1(4) = 1; tile2(4) = 6
         istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
         istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
         !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
         tile1(5) = 2; tile2(5) = 3
         istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
         istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
         !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
         tile1(6) = 2; tile2(6) = 4
         istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
         istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
         !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
         tile1(7) = 2; tile2(7) = 6
         istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
         istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
         !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
         tile1(8) = 3; tile2(8) = 4
         istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
         istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
         !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
         tile1(9) = 3; tile2(9) = 5
         istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
         istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
         !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
         tile1(10) = 4; tile2(10) = 5
         istart1(10) = 1;  iend1(10) = nx 
         jstart1(10) = ny; jend1(10) = ny
         istart2(10) = 1;  iend2(10) = nx 
         jstart2(10) = 1;  jend2(10) = 1
         !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
         tile1(11) = 4; tile2(11) = 6
         istart1(11) = nx; iend1(11) = nx
         jstart1(11) = 1;  jend1(11) = ny
         istart2(11) = nx; iend2(11) = 1  
         jstart2(11) = 1;  jend2(11) = 1
         !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
         tile1(12) = 5; tile2(12) = 6
         istart1(12) = nx; iend1(12) = nx
         jstart1(12) = 1;  jend1(12) = ny
         istart2(12) = 1;  iend2(12) = 1
         jstart2(12) = 1;  jend2(12) = ny
 
         call mpp_define_mosaic(global_indices, layout2D, domain, 
     &                        ntiles, num_contact, tile1, tile2, 
     &                        istart1, iend1, jstart1, jend1,
     &                        istart2, iend2, jstart2, jend2, 
     &                        pe_start=pe_start, pe_end=pe_end,
     &                        symmetry=.true., shalo = ng, nhalo = ng,
     &                        whalo = ng, ehalo = ng, name = type)

      end select

      deallocate(npes_tile, tile1, tile2)
      deallocate(istart1,iend1,jstart1,jend1)
      deallocate(istart2,iend2,jstart2,jend2)
      deallocate(layout2D, global_indices)
      deallocate(pe_start,pe_end)

      !--- find the tile number
      tile = mpp_pe()/npes_per_tile+1
      call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

!     tile=1
!     write(*,200) tile, is, ie, js, je
 200  format(i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ')

      end subroutine init_domain

      subroutine UPDATE_HALO_2D(array, domain, flags1)
      real array(:,:)
      type(domain2D) :: domain
      integer, optional :: flags1

      call mpp_update_domains ( array, domain, flags=SUPDATE+NUPDATE)
      !call mpp_update_domains ( array, domain, flags=NUPDATE)
      !call mpp_update_domains ( array, domain, flags=SUPDATE)

      end subroutine UPDATE_HALO_2D
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------
#endif


      END MODULE DOMAIN_DECOMP_1D

#ifdef DOMAIN_DECOMP_ATM_IS_1D
c If the atmosphere has a 1D domain decomposition, pass along the contents
c of DOMAIN_DECOMP_1D
      MODULE DOMAIN_DECOMP_ATM
      USE DOMAIN_DECOMP_1D
      END MODULE DOMAIN_DECOMP_ATM
#endif
