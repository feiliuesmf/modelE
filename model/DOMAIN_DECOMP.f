      MODULE DOMAIN_DECOMP
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@+    for the message passing (ESMF) implementation.
!@auth NCCS ASTG

#ifdef USE_ESMF
      use ESMF_Mod
#endif

      IMPLICIT NONE
#ifdef USE_ESMF
#include "mpif.h"
#endif
      SAVE
      PRIVATE ! Except for

#ifndef USE_ESMF
      ! Place holders for the real things
      TYPE ESMF_VM
         Integer :: i
      END TYPE ESMF_VM
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

!aoo since DIST_GRID is public ESMF_GRID has to be public
!aoo (SGI compiler complains)
      PUBLIC :: ESMF_GRID

      TYPE(ESMF_GridComp)  :: compmodelE
      TYPE (ESMF_DELayout) :: ESMF_LAYOUT

!@var DIST_GRID derived type to provide ESMF decomposition info
!@+   public components are used to minimize overhead for accessing
!@+   routine components
      PUBLIC :: DIST_GRID 
!@var  grid Default decomposition; globally accessible for convenience.
      PUBLIC :: grid
!@var INIT_APP Initialize default decomposition
      PUBLIC :: INIT_APP
      PUBLIC :: INIT_GRID
!@var FINISH_APP Cleans up at the end of the run (closes debugging file)
      PUBLIC :: FINISH_APP
!@var HALO_UPDATE Update data in halo for local domain using data from
!@+   neighbouring processes       
      PUBLIC :: HALO_UPDATE ! Communicate overlapping portions of subdomains
      PUBLIC :: HALO_UPDATE_COLUMN ! K, I, J
!@var CHECKSUM output a bit-reproducible checksum for an array
      PUBLIC :: CHECKSUM ! Communicate overlapping portions of subdomains
      PUBLIC :: CHECKSUMj! Communicate overlapping portions of subdomains
      PUBLIC :: CHECKSUM_COLUMN ! K, I, J
!@var GLOBALSUM output a bit-reproducible global-hemisphere-zonal sum for an array
      PUBLIC :: GLOBALSUM
!@var GLOBALMAX determine max value across pes
      PUBLIC :: GLOBALMAX
!@var ARRAYSCATTER scatter a global array to a decomposed array
      PUBLIC :: ARRAYSCATTER
!@var ARRAYGATHER gather a decomposed array to a global array
      PUBLIC :: ARRAYGATHER
#ifdef USE_ESMF
      PUBLIC :: ESMF_ARRAYGATHER
      PUBLIC :: ESMF_ARRAYSCATTER
#endif
!@var GET - extracts bounds information from DIST_GRID object
      PUBLIC :: GET
      PUBLIC :: HERE
      PUBLIC :: LOG_PARALLEL

      PUBLIC :: BACKSPACE_PARALLEL
      PUBLIC :: REWIND_PARALLEL
      PUBLIC :: SKIP_PARALLEL
      PUBLIC :: DREAD_PARALLEL
      PUBLIC :: MREAD_PARALLEL
      PUBLIC :: READT_PARALLEL
      PUBLIC :: READ_PARALLEL
      PUBLIC :: WRITE_PARALLEL
      PUBLIC :: READ_GRID_VAR
      PUBLIC :: WRITE_GRID_VAR

!@var HALO_UPDATE Generic wrapper for 2D and 3D routines
      INTERFACE HALO_UPDATE
        MODULE PROCEDURE HALO_UPDATE_1D  ! J
        MODULE PROCEDURE HALO_UPDATE_2D  ! I,J
        MODULE PROCEDURE HALO_UPDATE_3D  ! I,J,K
      END INTERFACE

      INTERFACE HALO_UPDATE_COLUMN
        MODULE PROCEDURE HALO_UPDATE_COLUMN_2D  ! M,J
        MODULE PROCEDURE HALO_UPDATE_COLUMN_3D  ! M,I,J
        MODULE PROCEDURE INT_HALO_UPDATE_COLUMN_3D  ! M,I,J
        MODULE PROCEDURE HALO_UPDATE_COLUMN_4D  ! M,I,J,K
      END INTERFACE

      INTERFACE CHECKSUM
        MODULE PROCEDURE CHECKSUM_1D
        MODULE PROCEDURE CHECKSUM_2D
        MODULE PROCEDURE CHECKSUM_3D
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

#ifdef USE_ESMF
      INTERFACE GHOST
        MODULE PROCEDURE GHOST_R8_1D
        MODULE PROCEDURE GHOST_R8_2D
        MODULE PROCEDURE GHOST_R8_3D
        MODULE PROCEDURE GHOST_R8_4D
      END INTERFACE

      INTERFACE GHOST_PAYLOAD
        MODULE PROCEDURE GHOST_PAY_R8_2D
      END INTERFACE

      interface COMM_BSR
        module procedure COMM_BSR_R8_0D
        module procedure COMM_BSR_R8_1D
        module procedure COMM_BSR_R8_2D
      end interface
#endif

      INTERFACE GLOBALSUM
        MODULE PROCEDURE GLOBALSUM_INT_REDUCE
        MODULE PROCEDURE GLOBALSUM_J
        MODULE PROCEDURE GLOBALSUM_IJ
        MODULE PROCEDURE GLOBALSUM_IJK
        MODULE PROCEDURE GLOBALSUM_IJK_IK
        MODULE PROCEDURE GLOBALSUM_JK
      END INTERFACE
 
      INTERFACE ARRAYSCATTER
        MODULE PROCEDURE ARRAYSCATTER_J
        MODULE PROCEDURE ARRAYSCATTER_IJ
        MODULE PROCEDURE IARRAYSCATTER_IJ
        MODULE PROCEDURE LARRAYSCATTER_IJ
      END INTERFACE   

#ifdef USE_ESMF
      INTERFACE ESMF_ARRAYGATHER
        MODULE PROCEDURE ESMF_ARRAYGATHER_J
        MODULE PROCEDURE ESMF_ARRAYGATHER_J_int
        MODULE PROCEDURE ESMF_ARRAYGATHER_IJ
        MODULE PROCEDURE ESMF_IARRAYGATHER_IJ
        MODULE PROCEDURE ESMF_LARRAYGATHER_IJ
      END INTERFACE   
      INTERFACE ESMF_ARRAYSCATTER
        MODULE PROCEDURE ESMF_ARRAYSCATTER_J
        MODULE PROCEDURE ESMF_ARRAYSCATTER_IJ
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

      INTERFACE MREAD_PARALLEL
        MODULE PROCEDURE MREAD_PARALLEL_2D
        MODULE PROCEDURE MREAD_PARALLEL_3D
      END INTERFACE

      INTERFACE READT_PARALLEL
        MODULE PROCEDURE READT_PARALLEL_2D
        MODULE PROCEDURE READT_PARALLEL_3D
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

      interface READ_GRID_VAR
         module procedure READ_GRID_VAR_8_2D
         module procedure READ_GRID_VAR_8_3D
         module procedure READ_GRID_VAR_8_4D
      end interface

      interface WRITE_GRID_VAR
         module procedure WRITE_GRID_VAR_8_2D
         module procedure WRITE_GRID_VAR_8_3D
         module procedure WRITE_GRID_VAR_8_4D
      end interface

      interface ESMF_READ_GRID_VAR
         module procedure ESMF_READ_GRID_VAR_8_2D
         module procedure ESMF_READ_GRID_VAR_8_3D
         module procedure ESMF_READ_GRID_VAR_8_4D
      end interface

      interface ESMF_WRITE_GRID_VAR
         module procedure ESMF_WRITE_GRID_VAR_8_2D
         module procedure ESMF_WRITE_GRID_VAR_8_3D
         module procedure ESMF_WRITE_GRID_VAR_8_4D
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
         module procedure PACK_2D       ! (i,j)
         module procedure LPACK_2D      ! (i,j)
         module procedure PACK_3D       ! (i,j,l)
         module procedure IPACK_3D      ! (i,j,l)
         module procedure PACK_4D       ! (i,j,l,m)
      end interface

      PUBLIC :: PACK_DATAj
      interface PACK_DATAj
         module procedure PACKj_2D     ! (j,k)
         module procedure PACKj_3D     ! (j,k,l)
         module procedure PACKj_4D     ! (j,k,l,m)
      end interface

!@var ESMF_BCAST Generic routine to broadcast data to all PEs.
      PUBLIC :: ESMF_BCAST
      INTERFACE ESMF_BCAST
        MODULE PROCEDURE ESMF_BCAST_1D
      END INTERFACE

!@var UNPACK Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
      PUBLIC :: UNPACK_DATA
      interface UNPACK_DATA
         module procedure UNPACK_1D      ! (i)
         module procedure UNPACK_2D      ! (i,j)
         module procedure LUNPACK_2D     ! (i,j)
         module procedure UNPACK_3D      ! (i,j,l)
         module procedure IUNPACK_3D     ! (i,j,l)
         module procedure UNPACK_4D      ! (i,j,l,m)
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
         module procedure  PACK_BLOCK_3D    ! (k,l,i,j,m)
      end interface

!@var UNPACK_BLOCK  Generic routine to unpack into a distributed
!@+   array the data from the corresponding global array.
      PUBLIC :: UNPACK_BLOCK 
      interface UNPACK_BLOCK 
         module procedure IUNPACK_BLOCK_2D    ! (k,l,i,j  )
         module procedure  UNPACK_BLOCK_2D    ! (k,l,i,j  )
         module procedure  UNPACK_BLOCK_3D    ! (k,l,i,j,m)
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


      ! Direction bits
      PUBLIC :: NORTH, SOUTH, EAST, WEST
      ! Use powers of two such that addition can be used
      ! to create combination
      INTEGER, PARAMETER :: NORTH = 2**0, SOUTH = 2**1
      INTEGER, PARAMETER :: EAST  = 2**2, WEST  = 2**3
      INTEGER, PARAMETER :: ALL = NORTH + SOUTH ! no east/west for now

      INTEGER, PARAMETER :: HALO_WIDTH = 1
      integer, parameter :: root=0

      ! Local grid information
      TYPE DIST_GRID

         TYPE (ESMF_Grid) :: ESMF_GRID
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

#ifdef DEBUG_DECOMP
         INTEGER :: log_unit ! for debugging
#endif
      END TYPE DIST_GRID

      TYPE (DIST_GRID) :: GRID

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

      TYPE (ESMF_DELayout) :: ESMF_LAYOUT_def
      Integer :: pe
      integer, parameter :: ROOT_ID=0
      INTEGER :: CHECKSUM_UNIT

      Integer :: tag = 10

      CONTAINS

      ! This routine initializes the quantities described above.
      ! The initialization should proceed prior to any grid computations.
      SUBROUTINE INIT_APP(grd_dum,IM,JM)
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM
      INTEGER             :: rc
      INTEGER             :: pet
      CHARACTER(LEN=20) :: buffer
#ifdef USE_ESMF
      Type (ESMF_VM) :: vm
#endif

#ifdef USE_ESMF
      ! Initialize ESMF
      Call ESMF_Initialize(vm=vm, rc=rc)
      Call ESMF_VMGet(vm, localPET = my_pet, petCount = NPES, rc=rc)
      compmodelE  = ESMF_GridCompCreate(vm,"ModelE ESMF", rc=rc)

      ! The default layout is not what we want - it splits in the "I" direction.
      ESMF_Layout = ESMF_DELayoutCreate(vm, 
     & dePetList=(/ (pet,pet=0,NPES) /), deCountList = (/ 1, NPES /))

      NP_LON = 1
      NP_LAT = NPES
      RANK_LAT = my_pet
      RANK_LON = 0
#else
      NP_LON = 1
      NP_LAT = 1
      RANK_LON = 0
      RANK_LAT = 0
#endif
      
      call INIT_GRID(grd_dum,IM,JM)
      WRITE(*,*)'INIT_APP: ', IM, JM, NP_LON, NP_LAT, RANK_LON, RANK_LAT

#ifdef DEBUG_DECOMP
      IF (AM_I_ROOT()) CALL openunit('CHKSUM_DECOMP', CHECKSUM_UNIT)
      WRITE(buffer,'(a,i3.3)') 'LOG_',my_pet
      CALL openunit(TRIM(buffer), grd_dum%log_unit)
#endif

      END SUBROUTINE INIT_APP

      SUBROUTINE INIT_GRID(grd_dum,IM,JM)
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM
      integer, parameter :: numDims=2
      integer, dimension(numDims) :: grid_size
      integer             :: rc
      real(ESMF_KIND_R8), dimension(numDims) :: range_min,range_max

      INTEGER :: J_EQUATOR
      INTEGER :: I0_DUM, I1_DUM
      INTEGER :: J0_DUM, J1_DUM

      grid_size(1)=IM;   grid_size(2)=JM
      range_min(1)=0.;   range_min(2)=-90.
      range_max(1)=360.; range_max(2)=90.

#ifdef USE_ESMF
      grd_dum%ESMF_GRID = ESMF_GridCreateHorzXYUni(counts=grid_size,
     &     minGlobalCoordPerDim=range_min, 
     &     maxGlobalCoordPerDim=range_max, 
     &     horzStagger=ESMF_GRID_HORZ_STAGGER_C_NE,
     &     name="source grid", rc=rc)

      Call ESMF_GridDistribute(grid=grd_dum%ESMF_GRID, 
     &     delayout = ESMF_Layout, rc=rc)
      Call ESMF_GridCompSet(compmodelE, grid=grd_dum%ESMF_GRID, rc=rc)
      Call ESMF_GRID_BOUNDS(grd_dum%ESMF_GRID, RANK_LON, RANK_LAT,
     &        I0_DUM, I1_DUM, J0_DUM, J1_DUM)

#else
      RANK_LON = 0
      RANK_LAT = 0
      I0_DUM = 1
      I1_DUM = IM
      J0_DUM = 1
      J1_DUM = JM
#endif

      grd_dum%IM_WORLD      = IM
      grd_dum%JM_WORLD      = JM

      ! Wrapped ESMF grid
      grd_dum%I_STRT        = I0_DUM
      grd_dum%I_STOP        = I1_DUM
      grd_dum%I_STRT_HALO   = MAX( 1, I0_DUM-HALO_WIDTH)
      grd_dum%I_STOP_HALO   = MIN(IM, I1_DUM+HALO_WIDTH)

      grd_dum%J_STRT        = J0_DUM
      grd_dum%J_STOP        = J1_DUM
      IF (RANK_LAT > 0) THEN
        grd_dum%J_STRT_SKP = J0_DUM
      ELSE
        grd_dum%J_STRT_SKP    = 2
      END IF
      IF (RANK_LAT < NP_LAT - 1) THEN
        grd_dum%J_STOP_SKP    = J1_DUM
      ELSE
        grd_dum%J_STOP_SKP    = JM-1
      END IF

      grd_dum%J_STRT_HALO   = MAX( 1, J0_DUM - HALO_WIDTH)
      grd_dum%J_STOP_HALO   = MIN(JM, J1_DUM + HALO_WIDTH)

      IF (RANK_LAT > 0) THEN
        grd_dum%J_STRT_STGR = J0_DUM
      ELSE
        grd_dum%J_STRT_STGR = 2
      ENDIF
      grd_dum%J_STOP_STGR   = J1_DUM

      grd_dum%HAVE_SOUTH_POLE = (RANK_LAT == 0)
      grd_dum%HAVE_NORTH_POLE = (RANK_LAT == NP_LAT - 1)

      J_EQUATOR = JM/2
      grd_dum%HAVE_EQUATOR    = 
     &      (J0_DUM <= J_EQUATOR) .AND. (J1_DUM >= J_EQUATOR)

      WRITE(*,*)'init_grid: ',my_pet,grd_dum%j_strt,grd_dum%j_stop
      END SUBROUTINE INIT_GRID

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

#ifdef USE_ESMF
      Call GHOST(grd_dum%ESMF_GRID, arr, from)
#endif

      END SUBROUTINE HALO_UPDATE_1D

      SUBROUTINE HALO_UPDATE_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                    arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_ESMF
      Call GHOST(grd_dum%ESMF_GRID, arr, from)
#endif
      END SUBROUTINE HALO_UPDATE_2D

      SUBROUTINE HALO_UPDATE_3D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                 arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L
#ifdef USE_ESMF
      Call GHOST(grd_dum%ESMF_GRID, arr, from)
#endif
      END SUBROUTINE HALO_UPDATE_3D

      SUBROUTINE HALO_UPDATE_COLUMN_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                  arr(:,grd_dum%j_strt_halo:)

      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_ESMF
      ! This case can be handled by GHOST_2D
      Call GHOST(grd_dum%ESMF_GRID, arr, IAND(from,NORTH+SOUTH))
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_2D

      SUBROUTINE HALO_UPDATE_COLUMN_3D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                  arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_ESMF
      call GHOST_PAYLOAD(grd_dum%ESMF_GRID, arr, FROM)
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

#ifdef USE_ESMF
      do L = 1, size(arr,4)
        call GHOST_PAYLOAD(grd_dum%ESMF_GRID, arr(:,:,:,L), FROM)
      end do
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_4D

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

      t_arr(:,J_0:J_1) = arr(:,J_0:J_1)
      Call GLOBALSUM(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call GLOBALSUM(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)

      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &     file,line, asum, L1norm

#endif

      END SUBROUTINE CHECKSUM_2D

      SUBROUTINE CHECKSUM_3D(grd_dum, arr, line, file, unit, stgr)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: 
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit
      LOGICAL, OPTIONAL, INTENT(IN) :: stgr


      INTEGER :: unit_
      INTEGER :: k
      REAL*8, DIMENSION(Size(arr,3))  :: asum, L1norm

      REAL*8 :: 
     &  t_arr(size(arr,1),grd_dum%j_strt_halo:grd_dum%j_stop_halo)
      INTEGER :: J_0, J_1
      Integer :: stgr_

#ifdef DEBUG_DECOMP

      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      stgr_ = 0
      If (Present(stgr)) THEN
         If (stgr) stgr_=1
      End If

      Do k = 1, Size(arr, 3)
         t_arr(:,J_0:J_1) = arr(:,J_0:J_1,k)
         Call GLOBALSUM(grd_dum, t_arr, asum(k), istag=stgr_)
         t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
         Call GLOBALSUM(grd_dum, t_arr, L1norm(k), istag=stgr_)
      End Do
      If (AM_I_ROOT()) Then
        Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &       file,line, Sum(asum), Sum(L1norm)
      End If

#endif

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

      t_arr(:,J_0:J_1) = Transpose(arr(J_0:J_1,:))
      Call GLOBALSUM(grd_dum,      t_arr,asum, istag=stgr_,iskip=skip_)
      t_arr(:,J_0:J_1) = ABS(t_arr(:,J_0:J_1))
      Call GLOBALSUM(grd_dum,      t_arr,L1norm,istag=stgr_,iskip=skip_)
 
      If (AM_I_ROOT()) Write(unit_,'(a20,1x,i6,1x,2(e24.17,1x))') 
     &     file,line, asum, L1norm

#endif

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

#ifdef DEBUG_DECOMP

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

#endif

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

#ifdef DEBUG_DECOMP
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

#endif

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

#ifdef DEBUG_DECOMP
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

#endif

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

#ifdef DEBUG_DECOMP
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

#endif

      END SUBROUTINE CHECKSUM_COLUMN_4D


      SUBROUTINE FINISH_APP()
      USE FILEMANAGER, ONLY : closeunit
      IMPLICIT NONE

      INTEGER :: ier

#ifdef DEBUG_DECOMP
      CALL closeunit(CHECKSUM_UNIT)
      CALL HERE(__FILE__,__LINE__)
      CALL closeunit(grid%log_unit)
#endif

#ifdef USE_ESMF      
      CALL MPI_Finalize(ier)
#endif

      END SUBROUTINE FINISH_APP

!---------------------------
      SUBROUTINE GLOBALSUM_INT_REDUCE(grd_dum, ivar, isum, all)
      IMPLICIT NONE
      TYPE(DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER, INTENT(IN)  :: ivar
      INTEGER, INTENT(OUT) :: isum
      LOGICAL, OPTIONAL, INTENT(IN) :: all

      LOGICAL :: all_
      INTEGER :: status

      all_ = .false.
      If (Present(all)) all_ = all

#ifdef USE_ESMF      
      If (all_) Then
         call MPI_Allreduce(ivar, isum, 1, MPI_INTEGER, MPI_SUM, 
     &        MPI_COMM_WORLD, status)
      Else
         call MPI_Reduce(ivar, isum, 1, MPI_INTEGER, MPI_SUM, root,
     &        MPI_COMM_WORLD, status)
      End If
#else
      isum = ivar
#endif
      
      END SUBROUTINE GLOBALSUM_INT_REDUCE

      SUBROUTINE GLOBALSUM_J(grd_dum, arr, gsum,
     &                       hsum, istag, iskip, polefirst,all, jband)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: arr(grd_dum%j_strt_halo:)
      REAL*8,            INTENT(OUT):: gsum
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2)
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      INTEGER,OPTIONAL,  INTENT(IN) :: iskip
      LOGICAL,OPTIONAL,  INTENT(IN) :: polefirst
      LOGICAL,OPTIONAL,  INTENT(IN) :: all
      INTEGER, OPTIONAL, INTENT(IN) :: jband(2)

      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, J, ier
      REAL*8  :: garr(grd_dum%jm_world)
      LOGICAL :: istag_, iskip_

      INTEGER :: JSTRT, JSTOP
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
#endif

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
        CALL HERE(__FILE__//' should never happen',__LINE__)
#ifdef USE_ESMF
        CALL MPI_ABORT(ier)
#else
        STOP
#endif
      End If
#endif

#ifdef USE_ESMF
      GRID = grd_dum%ESMF_GRID
      call Esmf_ArrayGather(GRID, arr(j_0:j_1), garr)
#else
      garr = arr(j_0:j_1)
#endif
    
      If (AM_I_ROOT()) then
         If (Present(polefirst)) Then
            If (polefirst) Then
               gsum = garr(1) + garr(JM)
               DO J = 2, JM-1
                  gsum = gsum + garr(J)
               END DO
            End IF
         Else
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
         Endif
      EndIf

#ifdef USE_ESMF
      If (Present(all)) Then
         If (all) THEN
            Call MPI_BCAST(gsum,1,MPI_REAL8,root,
     &           MPI_COMM_WORLD, ier)
            If (Present(hsum)) Call MPI_BCAST(hsum,2,MPI_REAL8,root,
     &              MPI_COMM_WORLD, ier)
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

      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, J,J_0STG, ier, J_0S,J_1S
      REAL*8  :: zon(grd_dum%j_strt:grd_dum%j_stop)
      REAL*8  :: garr(grd_dum%jm_world)
      LOGICAL :: istag_,iskip_
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
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

#ifdef USE_ESMF
      GRID = grd_dum%ESMF_GRID
      call Esmf_ArrayGather(GRID, zon(j_0:j_1), garr)
#else
      garr = zon(j_0:j_1)
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

#ifdef USE_ESMF
      If (Present(all)) Then
         If (all) Then
            Call MPI_BCAST(gsum,1,MPI_REAL8,root,
     &           MPI_COMM_WORLD, ier)
            If (Present(hsum))     Call MPI_BCAST(hsum,2,MPI_REAL8,root,
     &           MPI_COMM_WORLD, ier)
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
      type (ESMF_Grid)                           :: GRID
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

#ifdef USE_ESMF
      GRID = grd_dum%ESMF_GRID
      Do k=1,size(arr,3)
         call Esmf_ArrayGather(GRID, zon(j_0:j_1,k), garr(:,k))
      End Do
#else
      garr(:,:) = zon(j_0:j_1,:)
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

      SUBROUTINE GLOBALSUM_IJK_IK(grd_dum, arr, gsum, jband)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN) :: grd_dum
      REAL*8,             INTENT(IN) :: arr(grd_dum%i_strt_halo:,
     &                                      grd_dum%j_strt_halo:,:)
      REAL*8,             INTENT(OUT):: gsum(size(arr,1), size(arr,3))
      INTEGER,OPTIONAL,   INTENT(IN) :: jband(2)

      INTEGER :: k
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM, jb1, jb2
#ifdef USE_ESMF
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world)
#endif
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
#endif

      i_0  = grd_dum%i_strt
      i_1  = grd_dum%i_stop
      j_0  = grd_dum%j_strt
      j_1  = grd_dum%j_stop
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      jb1 = 1
      jb2 = JM
      If (Present(jband)) Then
         jb1 = jband(1)
         jb2 = jband(2)
      End If

#ifdef USE_ESMF
      GRID = grd_dum%ESMF_GRID
      Do k=1,size(arr,3)
         call Esmf_ArrayGather(GRID, arr(:,j_0:j_1,k), garr(:,:))
         If (AM_I_ROOT()) gsum(:,k) = Sum(garr(:,jb1:jb2),2)
      End Do
#else
      gsum = Sum(arr(:,jb1:jb2,:),2)
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
      INTEGER :: ier
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM
      REAL*8  :: garr(grd_dum%jm_world,size(arr,2))
      LOGICAL :: istag_

    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
#endif

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


#ifdef USE_ESMF
      GRID = grd_dum%ESMF_GRID
      Do k=1,size(arr,2)
         call Esmf_ArrayGather(GRID, arr(j_0:j_1,k), garr(:,k))
      End Do
#else
      garr(:,:) = arr(j_0:j_1,:)
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

#ifdef USE_ESMF
      If (Present(all)) Then
         If (all) Then
            Call MPI_BCAST(gsum,Size(gsum),MPI_REAL8,root,
     &           MPI_COMM_WORLD, ier)
            If (Present(hsum)) Call MPI_BCAST(hsum,size(hsum),MPI_REAL8,
     &           root, MPI_COMM_WORLD, ier)
         End If
      End If
#endif

      END SUBROUTINE GLOBALSUM_JK

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

      SUBROUTINE DREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,AVAR)
!@sum	DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth	NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME      !@var  NAME  name of record being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:) !@var  AOUT real*8 array
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array

      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      INTEGER :: IERR
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
      GRID = grd_dum%ESMF_GRID
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      call ESMF_ArrayScatter(GRID, AVAR(:,grd_dum%J_STRT:), AOUT)
#else
      AVAR=AOUT
#endif

      If (IERR==0) Then
         WRITE(6,*) "Read from file ",TRIM(NAME)
         RETURN
      Else
         WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
         call stop_model('DREAD_PARALLEL: READ ERROR',255)
      EndIf

      END SUBROUTINE DREAD_PARALLEL_2D

      SUBROUTINE DREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR)
!@sum	DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth	NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME        !@var  NAME  name of record being read
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var  AOUT real*8 array
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      INTEGER :: IERR
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
      GRID = grd_dum%ESMF_GRID
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      DO N=1,NM
         call ESMF_ArrayScatter(GRID, 
     &        AVAR(:,grd_dum%J_STRT:,N), 
     &        AOUT(:,:,N))
      END DO
#else
      AVAR=AOUT
#endif

      If (IERR==0) Then
         WRITE(6,*) "Read from file ",TRIM(NAME)
         RETURN
      Else
         WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
         call stop_model('DREAD_PARALLEL: READ ERROR',255)
      EndIf

      END SUBROUTINE DREAD_PARALLEL_3D

      SUBROUTINE MREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,M,NSKIP,AVAR)
!@sum	MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth	NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME      !@var  NAME  name of record being read
      INTEGER,      INTENT(OUT) :: M         !@var  M      initial integer
      INTEGER,      INTENT(IN)  :: NSKIP     !@var  NSKIP no. of R*4's to skip
      REAL*8,      INTENT(OUT)  :: AVAR(:,:) !@var  AOUT real*8 array
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      INTEGER :: IERR
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
      GRID = grd_dum%ESMF_GRID
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) M,(X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      call ESMF_ArrayScatter(GRID, AVAR, AOUT)
#else
      AVAR=AOUT
#endif

      If (IERR==0) Then
         WRITE(6,*) "Read from file ",TRIM(NAME)
         RETURN
      Else
         WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
         call stop_model('MREAD_PARALLEL: READ ERROR',255)
      EndIf

      END SUBROUTINE MREAD_PARALLEL_2D

      SUBROUTINE MREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,M,NSKIP,AVAR)
!@sum	MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth	NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME        !@var  NAME  name of record being read
      INTEGER,      INTENT(OUT) :: M           !@var  M      initial integer
      INTEGER,      INTENT(IN)  :: NSKIP       !@var  NSKIP no. of R*4's to skip
      REAL*8,      INTENT(OUT)  :: AVAR(:,:,:) !@var  AOUT real*8 array
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      INTEGER :: IERR
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
      GRID = grd_dum%ESMF_GRID
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) M,(X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      DO N=1,NM
         call ESMF_ArrayScatter(GRID, AVAR(:,:,N), AOUT(:,:,N))
      END DO
#else
      AVAR=AOUT
#endif

      If (IERR==0) Then
         WRITE(6,*) "Read from file ",TRIM(NAME)
         RETURN
      Else
         WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
         call stop_model('MREAD_PARALLEL: READ ERROR',255)
      EndIf

      END SUBROUTINE MREAD_PARALLEL_3D

      SUBROUTINE READT_PARALLEL_2D (grd_dum,IUNIT,NAME,NSKIP,AVAR,IPOS)
!@sum	READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth	NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME       !@var  NAME  name of record being read
      INTEGER,      INTENT(IN)  :: NSKIP      !@var  NSKIP no. of R*4's to skip
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS       !@var  IPOS  no. of recs. to advance
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: IERR
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
      GRID = grd_dum%ESMF_GRID
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=IERR)
         END DO
         READ (IUNIT, IOSTAT=IERR) TITLE, (X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      call ESMF_ArrayScatter(GRID, 
     &   AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP), AOUT)
#else
      AVAR=AOUT
#endif

      Call CHECKSUM(grd_dum, AVAR, __LINE__,__FILE__)

      If (IERR==0) Then
         WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
         RETURN
      Else
         WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': ',
     &                      TRIM(TITLE),' IOSTAT=',IERR
         call stop_model('READT_PARALLEL: READ ERROR',255)
      EndIf

      END SUBROUTINE READT_PARALLEL_2D

      SUBROUTINE READT_PARALLEL_3D (grd_dum,IUNIT,NAME,NSKIP,AVAR,IPOS)
!@sum	READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth	NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME         !@var  NAME  name of record being read
      INTEGER,      INTENT(IN)  :: NSKIP        !@var  NSKIP no. of R*4's to skip
      REAL*8,       INTENT(OUT) :: AVAR(:,grd_dum%J_STRT_HALO:,:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: IERR
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
      GRID = grd_dum%ESMF_GRID
#endif

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)

      If (AM_I_ROOT()) then
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=IERR)
         END DO
         READ (IUNIT, IOSTAT=IERR) TITLE, (X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
         CALL LOG_PARALLEL(grd_dum, __FILE__, __LINE__,
     &        x0 = SUM(AOUT**2),x1=Sum(Sum(AOUT**2,1),1))
      EndIf

#ifdef USE_ESMF
      DO N=1,NM
         call ESMF_ArrayScatter(GRID, 
     &        AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,N), AOUT(:,:,N))
      END DO
#else
      AVAR=AOUT
#endif

      Call CHECKSUM(grd_dum, AVAR, __LINE__,__FILE__)

      If (IERR==0) Then
         WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
         RETURN
      Else
         WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': 
     &                      ',TRIM(TITLE),' IOSTAT=',IERR
         call stop_model('READT_PARALLEL: READ ERROR',255)
      EndIf

      END SUBROUTINE READT_PARALLEL_3D

      subroutine ESMF_ArrayScatter_IJ(egrid, local_array, global_array)
      real (kind=8), dimension (:,:) :: local_array, global_array
      type (ESMF_Grid)      :: egrid
      
#ifdef USE_ESMF
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            :: 
     &     sendcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: recvcount
      
      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN
      
      
      real (kind=8), allocatable                    :: var(:)
      

      allocate (sendcounts(NPES), displs(0:NPES), stat=status)
      
      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
      else
         allocate(var(0:0), stat=status)
      endif
      
      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(egrid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)
      
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
        if (MY_PET == root) then
          
          var(displs(J):displs(I)-1) = 
     &         RESHAPE(global_array(I1:IN,J1:JN), 
     &         shape=(/sendcounts(I)/))
        endif
      enddo

      call MPI_ScatterV(
     &       var, sendcounts, displs, MPI_DOUBLE_PRECISION,
     &       local_array, recvcount,  MPI_DOUBLE_PRECISION,
     &       root, MPI_COMM_WORLD, status)
      
        deallocate(VAR, stat=status)
      
      deallocate(sendcounts, displs, stat=status)
      Deallocate(AI)
#else
      local_array=global_array
#endif
      
      
      end subroutine ESMF_ArrayScatter_IJ

      subroutine ESMF_IArrayScatter_IJ(egrid, local_array, global_array)
      integer      , dimension (:,:) :: local_array, global_array
      type (ESMF_Grid)      :: egrid

#ifdef USE_ESMF
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            ::
     &     sendcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: recvcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      integer      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(egrid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)

      allocate (sendcounts(NPES), displs(0:NPES), stat=status)

      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
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
     &     MPI_INTEGER         , root, MPI_COMM_WORLD, status)

      if (AM_I_ROOT()) then
        deallocate(VAR, stat=status)
      endif

      deallocate(sendcounts, displs, stat=status)
      deallocate(AI)
#else
      local_array=global_array
#endif

      end subroutine ESMF_IArrayScatter_IJ

!-----------------------------------------------------------
      subroutine ESMF_LArrayScatter_IJ(egrid, local_array, global_array)
      logical      , dimension (:,:) :: local_array, global_array
      type (ESMF_Grid)      :: egrid

#ifdef USE_ESMF
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            ::
     &     sendcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: recvcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      integer      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(egrid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)
      allocate (sendcounts(NPES), displs(0:NPES), stat=status)

      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
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
     &     MPI_LOGICAL         , root, MPI_COMM_WORLD, status)

      if (AM_I_ROOT()) then
        deallocate(VAR, stat=status)
      endif

      deallocate(sendcounts, displs, stat=status)
      deallocate(AI)
#else
      local_array=global_array
#endif
      end subroutine ESMF_LArrayScatter_IJ



      subroutine ESMF_ArrayScatter_J(grid, local_array, global_array)
      real (kind=8), dimension (:) :: local_array, global_array
      type (ESMF_Grid)      :: grid
      
#ifdef USE_ESMF
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            :: 
     &     sendcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: recvcount
      
      integer                                       :: I, J
      integer                                       :: NY
      integer                                       :: J1, JN
      
      
      real (kind=8), allocatable                    :: var(:)

      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)
      allocate (sendcounts(NPES), displs(0:NPES), stat=status)
      
      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
      endif
      
      displs(0) = 0
      do I = 1,NPES
        J = I - 1
        J1 = AI(I,2)%min
        JN = AI(I,2)%max
        
        sendcounts(I) = (JN - J1 + 1)
        if (J == MY_PET) then
          recvcount = sendcounts(I)
        endif
        displs(I) = displs(J) + sendcounts(I)
        if (AM_I_ROOT()) then
          
          var(displs(J):displs(I)-1) = global_array(J1:JN)
        endif
      enddo
      
      
      call MPI_ScatterV(var, sendcounts, displs, 
     &     MPI_DOUBLE_PRECISION, local_array, recvcount,
     &     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, status)
      
      
      if (AM_I_ROOT()) then
        deallocate(VAR, stat=status)
      endif
      
      deallocate(sendcounts, displs, stat=status)
      deallocate(AI)
#else
      local_array=global_array
#endif

      end subroutine ESMF_ArrayScatter_J

#ifdef USE_ESMF
      subroutine Esmf_ArrayGather_IJ(grid, local_array, global_array)
      type (ESMF_Grid)      :: grid
      real (kind=8), dimension (:,:) :: local_array, global_array
      
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      type (ESMF_DELayout)                          :: layout
      integer, allocatable, dimension(:)            :: 
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: sendcount
      
      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN
      
      
      real (kind=8), allocatable                    :: var(:)
      

      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)

      allocate (recvcounts(NPES), displs(0:NPES), stat=status)
      
      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
      
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
      
      
      call MPI_GatherV(local_array, sendcount, MPI_DOUBLE_PRECISION, 
     &     var, recvcounts, displs, 
     &     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, status)
      
      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          I1 = AI(I,1)%min
          IN = AI(I,1)%max
          J1 = AI(I,2)%min
          JN = AI(I,2)%max
          
          global_array(I1:IN,J1:JN) = 
     &         RESHAPE(var(displs(J):displs(I)-1), 
     &         shape=(/size(local_array,1),
     &         size(local_array,2)/))
        enddo
      endif
      
      deallocate(VAR, stat=status)
      deallocate(recvcounts, displs, stat=status)
      deallocate(AI)
      end subroutine Esmf_ArrayGather_IJ

!--------------------------------
      subroutine Esmf_IArrayGather_IJ(grid, local_array, global_array)
      type (ESMF_Grid)      :: grid
      integer      , dimension (:,:) :: local_array, global_array

      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      integer, allocatable, dimension(:)            ::
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: sendcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      integer      , allocatable                    :: var(:)


      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)

      allocate (recvcounts(NPES), displs(0:NPES), stat=status)

      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)

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
     &     MPI_INTEGER         , root, MPI_COMM_WORLD, status)

      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          I1 = AI(I,1)%min
          IN = AI(I,1)%max
          J1 = AI(I,2)%min
          JN = AI(I,2)%max

          global_array(I1:IN,J1:JN) =
     &         RESHAPE(var(displs(J):displs(I)-1),
     &         shape=(/size(local_array,1),
     &         size(local_array,2)/))
        enddo
      endif
      deallocate(VAR, stat=status)
      deallocate(recvcounts, displs, stat=status)
      DEALLOCATE(AI)


      end subroutine Esmf_IArrayGather_IJ

!--------------------------------
      subroutine Esmf_LArrayGather_IJ(e_grid, local_array, global_array)
      type (ESMF_Grid)      :: e_grid
      logical      , dimension (:,:) :: local_array, global_array

      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      integer, allocatable, dimension(:)            ::
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: sendcount

      integer                                       :: I, J
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN

      logical      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(e_grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)
      allocate (recvcounts(NPES), displs(0:NPES), stat=status)

      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)


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
     &     MPI_LOGICAL, root, MPI_COMM_WORLD, status)

      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          I1 = AI(I,1)%min
          IN = AI(I,1)%max
          J1 = AI(I,2)%min
          JN = AI(I,2)%max

          global_array(I1:IN,J1:JN) =
     &         RESHAPE(var(displs(J):displs(I)-1),
     &         shape=(/size(local_array,1),
     &         size(local_array,2)/))
        enddo
      endif

      deallocate(VAR, stat=status)
      deallocate(recvcounts, displs, stat=status)
      Deallocate(AI)
      end subroutine Esmf_LArrayGather_IJ
!------------------------------------------------------------



      subroutine Esmf_ArrayGather_J(grid, local_array, global_array)
      type (ESMF_Grid)      :: grid
      real (kind=8), dimension (:) :: local_array, global_array
      
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
      integer, allocatable, dimension(:)            :: 
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: sendcount
      
      integer                                       :: I,J
      integer                                       :: NY
      integer                                      :: J1, JN
      
      
      real (kind=8), allocatable                    :: var(:)
      
      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)

      allocate (recvcounts(NPES), displs(0:NPES), stat=status)
      
      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
      
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
      
      
      call MPI_GatherV(local_array,sendcount,MPI_DOUBLE_PRECISION,
     &     var, recvcounts, displs, 
     &     MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, status)
      
      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          J1 = AI(I,2)%min
          JN = AI(I,2)%max
          
          global_array(J1:JN) = var(displs(J):displs(I)-1)
        enddo
      endif
      
      deallocate(VAR, stat=status)
      deallocate(recvcounts, displs, stat=status)
      deallocate(AI)
      end subroutine Esmf_ArrayGather_J

!---------------------------

      subroutine Esmf_ArrayGather_J_int(grid, local_array, global_array)
      type (ESMF_Grid)      :: grid
      INTEGER, dimension (:) :: local_array, global_array
      
      type(ESMF_AxisIndex), dimension(:,:), pointer :: AI

      integer, allocatable, dimension(:)            :: 
     &     recvcounts, displs
      integer                                       :: nDEs
      integer                                       :: status
      integer                                       :: sendcount
      
      integer                                       :: I,J
      integer                                       :: NY
      integer                                      :: J1, JN
      
      
      integer, allocatable                    :: var(:)
      
      
      Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)

      allocate (recvcounts(NPES), displs(0:NPES), stat=status)
      
      allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
      
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
     &     MPI_INTEGER, root, MPI_COMM_WORLD, status)
      
      if (AM_I_ROOT()) then
        do I = 1,NPES
          J = I - 1
          J1 = AI(I,2)%min
          JN = AI(I,2)%max
          
          global_array(J1:JN) = var(displs(J):displs(I)-1)
        enddo
      endif
      
      deallocate(VAR, stat=status)
      deallocate(recvcounts, displs, stat=status)
      Deallocate(AI)
      
      end subroutine Esmf_ArrayGather_J_int

#endif
!---------------------------

      subroutine ArrayGather_J(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:)               :: global_array
      
#ifdef USE_ESMF
      call Esmf_ArrayGather_J(grd_dum%ESMF_GRID, 
     &     local_array(grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array = local_array
#endif

      end subroutine ArrayGather_J

!---------------------------

      subroutine ArrayGather_J_int(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      INTEGER, dimension (grd_dum%J_STRT_HALO:) :: local_array
      INTEGER, dimension (:)               :: global_array
      
#ifdef USE_ESMF
      call Esmf_ArrayGather_J_int(grd_dum%ESMF_GRID, 
     &     local_array(grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array = local_array
#endif

      end subroutine ArrayGather_J_int

!---------------------------

      subroutine ArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:,:)               :: global_array
      
#ifdef USE_ESMF
      call Esmf_ArrayGather_IJ(grd_dum%ESMF_GRID,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array = local_array
#endif

      end subroutine ArrayGather_IJ

!---------------------------

      subroutine LArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      logical, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      logical, dimension (:,:)                    :: global_array

#ifdef USE_ESMF
          call Esmf_LArrayGather_IJ(grd_dum%ESMF_GRID,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array = local_array
#endif

      end subroutine LArrayGather_IJ

!---------------------------


      subroutine IArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      integer, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      integer, dimension (:,:)               :: global_array

#ifdef USE_ESMF
      call Esmf_IArrayGather_IJ(grd_dum%ESMF_GRID,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      global_array = local_array
#endif

      end subroutine IArrayGather_IJ


    !---------------------------
!---------------------------

      subroutine ArrayScatter_J(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:)               :: global_array
      
#ifdef USE_ESMF
      call Esmf_ArrayScatter_J(grd_dum%ESMF_GRID,
     &     local_array(grd_dum%J_STRT:grd_dum%J_STOP),
     &     global_array)
#else
      local_array=global_array
#endif

      end subroutine ArrayScatter_J

!---------------------------

      subroutine ArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:,:)               :: global_array
      
#ifdef USE_ESMF
      call Esmf_ArrayScatter_IJ(grd_dum%ESMF_GRID, 
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     & global_array)
#else
      local_array = global_array
#endif

      end subroutine ArrayScatter_IJ

!---------------------------

      subroutine IArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      integer, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      integer, dimension (:,:)               :: global_array

#ifdef USE_ESMF
      call Esmf_IArrayScatter_IJ(grd_dum%ESMF_GRID,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     & global_array)
#else
      local_array = global_array
#endif

      end subroutine IArrayScatter_IJ

!---------------------------

      subroutine LArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      logical, dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      logical, dimension (:,:)               :: global_array

#ifdef USE_ESMF
      call Esmf_LArrayScatter_IJ(grd_dum%ESMF_GRID,
     &     local_array(:,grd_dum%J_STRT:grd_dum%J_STOP),
     & global_array)
#else
      local_array = global_array
#endif

      end subroutine LArrayScatter_IJ


    !---------------------------
#ifdef USE_ESMF
      subroutine ESMF_GRID_PE_LAYOUT  (GRID, NX, NY)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: NX, NY
   

    ! local vars
   
        integer :: counts(2)
        integer :: status
   
        NX = 1
        NY = NPES

      end subroutine ESMF_GRID_PE_LAYOUT

    !---------------------------

      subroutine ESMF_GRID_MY_PE_LOC  (GRID,NX0,NY0)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)          :: NX0, NY0
   
    ! local vars
        type (ESMF_DELayout) :: layout
   
        integer :: status
   
        call ESMF_GridGet(grid, delayout=layout, rc=status)
        NX0 = 0
        NY0 = my_pet

      end subroutine ESMF_GRID_MY_PE_LOC

    !---------------------------
#endif

      subroutine WRITE_PARALLEL_INTEGER_0 ( data, UNIT, format)
        
        INTEGER, intent(in   )            :: data
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        
        character(len=ESMF_MAXSTR) :: FORMATTED
        

        if (AM_I_ROOT()) then
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
      subroutine WRITE_PARALLEL_INTEGER_1 ( data, UNIT, format)
        
        INTEGER, intent(in   )            :: data (:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        
        character(len=ESMF_MAXSTR) :: FORMATTED
        
        
        if (AM_I_ROOT()) then
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
    
      subroutine WRITE_PARALLEL_REAL8_0 ( data, UNIT, format)
        
        REAL (KIND=8), intent(in   )            :: data 
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        
        character(len=ESMF_MAXSTR) :: FORMATTED
        
        
        if (AM_I_ROOT()) then
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
      subroutine WRITE_PARALLEL_REAL8_1 ( data, UNIT, format)
        
        REAL (KIND=8), intent(in   )            :: data(:)
        integer,            intent(in   ),  optional :: UNIT
        character(len=*),   intent(in   ),  optional :: format
        
        character(len=ESMF_MAXSTR) :: FORMATTED
        
        
        if (AM_I_ROOT()) then
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
    
      subroutine WRITE_PARALLEL_STRING_0 (STRING, FORMAT, UNIT)
        character(LEN=*), intent(in   )            :: STRING
        character(LEN=*), intent(in   ), optional  :: FORMAT
        integer,          intent(in   ),  optional :: UNIT
    
        if (AM_I_ROOT()) then
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
    
      subroutine WRITE_PARALLEL_STRING_1 (STRING, FORMAT, UNIT)
        character(LEN=*), intent(in   )            :: STRING(:)
        character(LEN=*), intent(in   ), optional  :: FORMAT
        integer,          intent(in   ),  optional :: UNIT
    
        if (AM_I_ROOT()) then
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
        integer                :: status
    
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
      
#ifdef USE_ESMF
        CALL MPI_BCAST(data, 1, MPI_INTEGER, root, 
     &       MPI_COMM_WORLD, status)
#endif
    
        return
 9999 continue
        call exit_rc(99)
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
        integer                :: status
    
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
      
#ifdef USE_ESMF
        CALL MPI_BCAST(data, size(data), MPI_INTEGER, root, 
     &       MPI_COMM_WORLD, status)
#endif
    
        return
 9999 continue
        call exit_rc(99)
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
        integer                :: status
    
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
      
#ifdef USE_ESMF
        CALL MPI_BCAST(data, Size(data), MPI_DOUBLE_PRECISION, root, 
     &       MPI_COMM_WORLD, status)
#endif
    
        return
 9999 continue
        call exit_rc(99)
      END SUBROUTINE ESMF_READ_PARALLEL_REAL8_1

    !---------------------------

#ifdef USE_ESMF
      subroutine ESMF_GRID_BOUNDS(GRID, X_LOC, Y_LOC, I1, IN, J1, JN)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(IN)          :: X_LOC, Y_LOC
        integer, intent(OUT)         :: I1, IN, J1, JN
   
    ! local vars
        integer :: deId
        integer :: status
   
        type(ESMF_AxisIndex), dimension(:,:), pointer :: AI

        Allocate(AI(1:NPES,2))
        call ESMF_GridGetAllAxisIndex(grid, AI, 
     &       horzRelLoc=ESMF_CELL_CENTER, rc=status)
        deId = ESMF_GRID_PE_NUM_FM_PE_LOC(grid, X_LOC, Y_LOC)
   

    ! AI uses 1-based index for deId
        deId = deId + 1
   
        I1 = AI(deId,1)%min
        IN = AI(deId,1)%max
        J1 = AI(deId,2)%min
        JN = AI(deId,2)%max
   
        Deallocate(AI)
      end subroutine ESMF_GRID_BOUNDS
#endif
    !---------------------------

      function ESMF_AM_I_ROOT(layout) result(R)
        type (ESMF_DELayout) :: layout
        logical              :: R

        integer :: status

#ifdef USE_ESMF
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
        integer :: status

#ifdef USE_ESMF
        R = .false.
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, status)
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
        integer :: status

#ifdef USE_ESMF
        call ESMF_GridGet(grid, delayout=layout, rc=status)
        deid = y_loc
#else
        deid = 1
#endif


      end function ESMF_GRID_PE_NUM_FM_PE_LOC


    !---------------------------

      subroutine READ_GRID_VAR_8_2D(UNIT,grd_dum,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
        integer                     , intent(IN)   :: UNIT
        type (DIST_GRID)            , intent(IN)   :: grd_dum
        real (kind=8)               , intent(OUT)  :: A(:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM

    ! now local
        type (ESMF_Grid)                           :: GRID

        GRID = grd_dum%ESMF_GRID
        call ESMF_READ_GRID_VAR(UNIT,GRID,A,CHECK_SUM_STRING,CHECK_SUM)

      END SUBROUTINE READ_GRID_VAR_8_2D

    !---------------------------

      subroutine ESMF_READ_GRID_VAR_8_2D(UNIT,GRID,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
    
        integer                     , intent(IN)   :: UNIT
        type (ESMF_Grid)            , intent(IN)   :: GRID
        real (kind=8)               , intent(OUT)  :: A(:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM
    
    ! Local variables
        real (kind=8),  allocatable                :: VAR(:,:)
        integer                                    :: IM_WORLD
        integer                                    :: JM_WORLD
        integer                                    :: status
        type (ESMF_DELayout)                       :: layout
    
    

        call ESMF_GRID_WORLD(GRID,IM_WORLD,JM_WORLD)
#ifdef USE ESMF
        call ESMF_GridGet(grid, delayout=layout, rc=status)
#endif

        allocate(VAR(IM_WORLD,JM_WORLD), stat=status)
    
    
        if (ESMF_am_i_root(layout)) then
           read (UNIT) VAR
           if (present(CHECK_SUM_STRING)) then
              print *, "In READ_GRID_VAR_2D_8 reading ", 
     &                              IM_WORLD,JM_WORLD, 
     &             CHECK_SUM_STRING
           else
              print *, "In READ_GRID_VAR_2D_8 reading ", 
     &                              IM_WORLD,JM_WORLD
           end if
#ifdef READY
           if (present(CHECK_SUM_STRING)) then
              call PRINT_CHECKSUM(CHECK_SUM_STRING, sum(VAR))
           end if
#endif
           if (present(CHECK_SUM)) then
              CHECK_SUM = sum(VAR)
           end if
        end if
#ifdef USE_ESMF
        call ArrayScatter(grid, A, VAR)
#else
      A = VAR
#endif
        
        deallocate(VAR)
        
        return
#ifdef USE_ESMF
        contains

          subroutine ArrayScatter(grid, local_array, global_array)
            real (kind=8), dimension (:,:) :: 
     &                              local_array, global_array
            type (ESMF_Grid)      :: grid
    
            type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
            type (ESMF_DELayout)                          :: layout
            integer, allocatable, dimension(:)            :: 
     &                                sendcounts, displs
            integer                                       :: nDEs
            integer                                       :: status
            integer                                       :: recvcount
    
            integer                                       :: I,J
            integer                                       :: NX, NY
            integer                                       :: I1, IN
            integer                                       :: J1, JN
    
    
            real (kind=8), allocatable                    :: var(:)
    
    
        Allocate(AI(1:NPES,2))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI, 
     &     horzRelLoc=ESMF_CELL_CENTER, rc=status)
    
            allocate (sendcounts(NPES), displs(0:NPES), stat=status)
    
            if (AM_I_ROOT()) then
               allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
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
     &               RESHAPE(global_array(I1:IN,J1:JN), 
     &               shape=(/sendcounts(I)/))
               endif
            enddo
               
    
            call MPI_ScatterV(var, sendcounts, displs, 
     &       MPI_DOUBLE_PRECISION, local_array,recvcount, 
     &       MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, status)
    
    
            if (AM_I_ROOT()) then
               deallocate(VAR, stat=status)
            endif
    
            deallocate(sendcounts, displs, stat=status)
            deallocate(AI)
          end subroutine ArrayScatter
#endif    
      end subroutine ESMF_READ_GRID_VAR_8_2D
    
    !---------------------------

      subroutine READ_GRID_VAR_8_3D(UNIT,grd_dum,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
        integer                     , intent(IN)   :: UNIT
        type (DIST_GRID)            , intent(IN)   :: grd_dum
        real (kind=8)               , intent(OUT)  :: A(:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM

    ! now local
        type (ESMF_Grid)                           :: GRID

        GRID = grd_dum%ESMF_GRID
        call ESMF_READ_GRID_VAR(UNIT,GRID,A,CHECK_SUM_STRING,CHECK_SUM)

      END SUBROUTINE READ_GRID_VAR_8_3D

    !---------------------------
      subroutine ESMF_READ_GRID_VAR_8_3D(UNIT,GRID,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
    
        integer                     , intent(IN)   :: UNIT
        type (ESMF_Grid)            , intent(IN)   :: GRID
        real (kind=8)               , intent(OUT)  :: A(:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM
    
        integer :: L
    
        do L = 1, size(A,3)
           call ESMF_READ_GRID_VAR(UNIT, GRID, A(:,:,L), 
     &          check_sum_string=check_sum_string, check_sum=check_sum)
        end do
    
        return
      end subroutine ESMF_READ_GRID_VAR_8_3D
      
    !---------------------------

      subroutine READ_GRID_VAR_8_4D(UNIT,grd_dum,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
        integer                     , intent(IN)   :: UNIT
        type (DIST_GRID)            , intent(IN)   :: grd_dum
        real (kind=8)               , intent(OUT)  :: A(:,:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM

    ! now local
        type (ESMF_Grid)                           :: GRID

        GRID = grd_dum%ESMF_GRID
        call ESMF_READ_GRID_VAR(UNIT,GRID,A,CHECK_SUM_STRING,CHECK_SUM)

      END SUBROUTINE READ_GRID_VAR_8_4D

    !---------------------------

      subroutine ESMF_READ_GRID_VAR_8_4D(UNIT,GRID,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
    
        integer                     , intent(IN)   :: UNIT
        type (ESMF_Grid)            , intent(IN)   :: GRID
        real (kind=8)               , intent(OUT)  :: A(:,:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM
    
        integer :: L
    
        do L = 1, size(A,4)
           call ESMF_READ_GRID_VAR(UNIT, GRID, A(:,:,:,L), 
     &          check_sum_string=check_sum_string, check_sum=check_sum)
        end do
    
        return
      end subroutine ESMF_READ_GRID_VAR_8_4D
    
    !---------------------------

      subroutine WRITE_GRID_VAR_8_2D(UNIT,grd_dum,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
        integer                     , intent(IN)   :: UNIT
        type (DIST_GRID)            , intent(IN)   :: grd_dum
        real (kind=8)               , intent(OUT)  :: A(:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM

    ! now local
        type (ESMF_Grid)                           :: GRID

        GRID = grd_dum%ESMF_GRID
        call ESMF_WRITE_GRID_VAR(UNIT,GRID,A,CHECK_SUM_STRING,CHECK_SUM)

      end subroutine WRITE_GRID_VAR_8_2D

    !---------------------------

      subroutine ESMF_WRITE_GRID_VAR_8_2D(UNIT,GRID,A,CHECK_SUM_STRING,
     &                              CHECK_SUM)
    
    
        integer                     , intent(IN)   :: UNIT
        type (ESMF_Grid)            , intent(IN)   :: GRID
        real (kind=8)               , intent(IN)   :: A(:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM
    
    ! Local variables
        real (kind=8),  allocatable                :: VAR(:,:)
        integer                                    :: IM_WORLD
        integer                                    :: JM_WORLD
        integer                                    :: status
        integer                                    :: I
        type (ESMF_DELayout)                       :: layout
#ifdef USE_ESMF
        type(ESMF_AxisIndex),dimension(2)         :: AI_local,AI_global
        integer, parameter,  dimension(2)         :: decompids=(/1,2/)
#endif
    
        call ESMF_GRID_WORLD(GRID,IM_WORLD,JM_WORLD)
    
        allocate(VAR(IM_WORLD,JM_WORLD), stat=status)
        var = -9.0 ! dummy init
    
#ifdef USE_ESMF
        call ESMF_GridGet(grid, delayout=layout, rc=status)
    
        call Esmf_ArrayGather(grid, A, VAR)
#else
      VAR = A
#endif
#ifdef LAYOUTGATHERARRAY_WORKS
    !    call ESMF_GridGetDE(grid, lcellexc_index=AI_local, rc=status)
    
        do I = 1, 2
           AI_local(I)%r  = AI_local(I)%r - AI_local(I)%l
           AI_local(I)%l  = 0
           AI_global(I)%l = 0
        enddo
        AI_global(1)%r = IM_WORLD - 1
        AI_global(2)%r = JM_WORLD - 1 
    
#ifdef USE_ESMF
    !    call ESMF_DELayoutGatherArrayR(layout, A, decompids, 
    !&                                  AI_local, AI_global, VAR, rc=status)
        call c_ESMC_DELayoutGatherArrayR(layout, A, decompids, 2, 
     &                                   AI_local, AI_global, 
     &                                   VAR, status)
#else
      VAR = A
#endif
#endif
        if (ESMF_am_i_root(layout)) then
           write (UNIT) VAR
    
#ifdef READY
           if (present(CHECK_SUM_STRING)) then
              call PRINT_CHECKSUM(CHECK_SUM_STRING, sum(VAR))
           end if
#endif
           if (present(CHECK_SUM)) then
              CHECK_SUM = sum(VAR)
           end if
        end if
        
        deallocate(VAR)
        
        return
      end subroutine ESMF_WRITE_GRID_VAR_8_2D
    
    !---------------------------

      subroutine WRITE_GRID_VAR_8_3D(UNIT,grd_dum,A,CHECK_SUM_STRING,
     &                               CHECK_SUM)
    
        integer                     , intent(IN)   :: UNIT
        type (DIST_GRID)            , intent(IN)   :: grd_dum
        real (kind=8)               , intent(OUT)  :: A(:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM

    ! now local
        type (ESMF_Grid)                           :: GRID

        GRID = grd_dum%ESMF_GRID
        call ESMF_WRITE_GRID_VAR(UNIT,GRID,A,CHECK_SUM_STRING,CHECK_SUM)

      end subroutine WRITE_GRID_VAR_8_3D
    
    !---------------------------

      subroutine ESMF_WRITE_GRID_VAR_8_3D(UNIT,GRID,A,CHECK_SUM_STRING,
     &                                    CHECK_SUM)
    
    
        integer                     , intent(IN)   :: UNIT
        type (ESMF_Grid)            , intent(IN)   :: GRID
        real (kind=8)               , intent(IN)   :: A(:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM
    
        integer :: L
    
        do L = 1, size(A,3)
           call ESMF_WRITE_GRID_VAR(UNIT, GRID, A(:,:,L),
     &           check_sum=check_sum)
        end do
    
        return
      end subroutine ESMF_WRITE_GRID_VAR_8_3D
      
    !---------------------------

      subroutine WRITE_GRID_VAR_8_4D(UNIT,grd_dum,A,CHECK_SUM_STRING,
     &                               CHECK_SUM)
    
        integer                     , intent(IN)   :: UNIT
        type (DIST_GRID)            , intent(IN)   :: grd_dum
        real (kind=8)               , intent(OUT)  :: A(:,:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM

    ! now local
        type (ESMF_Grid)                           :: GRID

        GRID = grd_dum%ESMF_GRID
        call ESMF_WRITE_GRID_VAR(UNIT,GRID,A,CHECK_SUM_STRING,CHECK_SUM)

      end subroutine WRITE_GRID_VAR_8_4D
    
    !---------------------------
      subroutine ESMF_WRITE_GRID_VAR_8_4D(UNIT,GRID,A,CHECK_SUM_STRING,
     &                                    CHECK_SUM)
    
    
        integer                     , intent(IN)   :: UNIT
        type (ESMF_Grid)            , intent(IN)   :: GRID
        real (kind=8)               , intent(IN)   :: A(:,:,:,:)
        character (LEN=*), optional , intent(IN)   :: CHECK_SUM_STRING
        real (kind=8),     optional , intent(OUT)  :: CHECK_SUM
    
        integer :: L
    
        do L = 1, size(A,4)
           call ESMF_WRITE_GRID_VAR
     &          (UNIT, GRID, A(:,:,:,L),check_sum=check_sum)
        end do
    
        return
      end subroutine ESMF_WRITE_GRID_VAR_8_4D
    
    !---------------------------

#ifdef USE_ESMF
      subroutine ESMF_GRID_WORLD(GRID,IM_WORLD,JM_WORLD)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: IM_WORLD, JM_WORLD
    
! local vars
        integer status
        integer, dimension(ESMF_MAXGRIDDIM) :: dims

        call ESMF_GridGet(grid, globalcellcountperdim=dims, rc=status)
   
        IM_WORLD = dims(1)
        JM_WORLD = dims(2)

      end subroutine ESMF_GRID_WORLD
#endif

!---------------------------

#ifdef USE_ESMF

      SUBROUTINE GHOST_R8_1D(GRID,A,FROM)
      TYPE(ESMF_Grid),     INTENT(IN)    :: GRID
      REAL*8,              INTENT(INOUT) :: A(:)
      INTEGER, OPTIONAL,   INTENT(IN)  :: FROM
      
      
      INTEGER :: LAST, LEN1, LEN2, STR
      INTEGER :: USABLE_FROM

      integer :: nn_north, nn_south, nn_east, nn_west
      integer :: j_north, j_south
      integer                 :: npx, npy
      integer                 :: px, py
      type (ESMF_DELayout)           :: layout
      integer                 :: status

      REAL*8 :: junk(1) = 0.
      
      call ESMF_GridGet(grid, delayout=layout, rc=status)
      call ESMF_GRID_PE_LAYOUT(grid, npx, npy)
      call ESMF_GRID_MY_PE_LOC(grid,  px,  py) 
      
      j_north = MOD(py  + 1, npy)
      j_south = MOD(py+npy-1,npy)
      
      nn_north = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, px, j_north)
      nn_south = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, px, j_south)
      
!       The default for FROM is to do all four sides.
      
      IF(.NOT.PRESENT(FROM)) THEN
        USABLE_FROM = ALL
      ELSE
        USABLE_FROM = FROM
      ENDIF


#define FILL(N) IAND(USABLE_FROM,N)==N

!       Fill NORTHERN ghost region

      IF(FILL(NORTH) .OR.
     &  FILL(SOUTH)) THEN
        LAST = SIZE(A,1)-1 
        LEN2 = 1
        STR = LEN2*SIZE(A,1)
        LEN1 = SIZE(A)/STR

        IF(FILL(NORTH)) THEN
           IF (py == npy - 1) THEN ! No northern neigbor
              CALL COMM_BSR( layout, junk(1),
     &             A(2    ), NN_SOUTH, LEN1,LEN2,STR,STR)
           ELSE
              CALL COMM_BSR( layout, A(last+1  ),
     &             A(2    ), NN_SOUTH, LEN1,LEN2,STR,STR)
           END IF
        ENDIF
          
!       Fill SOUTHERN ghost region

        IF(FILL(SOUTH)) THEN
           IF (py == 0) THEN ! No southern neigbbor
              CALL COMM_BSR( layout, junk(1),
     &             A(last   ), NN_NORTH, LEN1,LEN2,STR,STR)
           ELSE
              CALL COMM_BSR( layout, A(1    ),
     &             A(last   ), NN_NORTH, LEN1,LEN2,STR,STR)
           END IF
        ENDIF

      ENDIF

      RETURN

#undef FILL

      end SUBROUTINE GHOST_R8_1D

      SUBROUTINE GHOST_R8_2D(GRID,A,FROM)
      TYPE(ESMF_Grid),     INTENT(IN)    :: GRID
      REAL*8,              INTENT(INOUT) :: A(:,:)
      INTEGER, OPTIONAL,   INTENT(IN)  :: FROM
      
      
      INTEGER :: LAST, LEN1, LEN2, STR
      INTEGER :: USABLE_FROM

      integer :: nn_north, nn_south, nn_east, nn_west
      integer :: j_north, j_south, i_east, i_west
      integer                 :: npx, npy
      integer                 :: px, py
      type (ESMF_DELayout)           :: layout
      integer                 :: status
      REAL*8 :: junk(size(A,1))
      
      call ESMF_GridGet(grid, delayout=layout, rc=status)
      call ESMF_GRID_PE_LAYOUT(grid, npx, npy)
      call ESMF_GRID_MY_PE_LOC(grid,  px,  py) 
      
      j_north = MOD(py  +1,npy)
      j_south = MOD(py+npy-1,npy)
      i_east = MOD(px  +1,npx)
      i_west = MOD(px+npx-1,npx)
      
      nn_north = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, px, j_north)
      nn_south = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, px, j_south)
      
      nn_west = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, i_west, py)
      nn_east = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, i_east, py)
      
!       The default for FROM is to do all four sides.
      
      IF(.NOT.PRESENT(FROM)) THEN
        USABLE_FROM = ALL
      ELSE
        USABLE_FROM = FROM
      ENDIF


#define FILL(N) IAND(USABLE_FROM,N)==N

!       Fill NORTHERN ghost region

      IF(FILL(NORTH) .OR.
     &  FILL(SOUTH)) THEN
        LAST = SIZE(A,2)-1 
        LEN2 = SIZE(A,1)
        STR = LEN2*SIZE(A,2)
        LEN1 = SIZE(A)/STR
 
        IF(FILL(NORTH)) THEN
           IF (py == npy-1) THEN  ! No northern neigbor
              CALL COMM_BSR( layout, junk(:),
     &             A(:,2    ), NN_SOUTH, LEN1,LEN2,STR,STR)
           ELSE
              CALL COMM_BSR( layout, A(:,last+1  ),
     &             A(:,2    ), NN_SOUTH, LEN1,LEN2,STR,STR)
           END IF
        ENDIF
          
!       Fill SOUTHERN ghost region

        IF(FILL(SOUTH)) THEN
           IF (py == 0) THEN ! No southern neigbbor
              CALL COMM_BSR( layout, junk(:),
     &             A(:,last   ), NN_NORTH, LEN1,LEN2,STR,STR)
           ELSE
              CALL COMM_BSR( layout, A(:,1    ),
     &             A(:,last   ), NN_NORTH, LEN1,LEN2,STR,STR)
           END IF
        END IF

      END IF
cgsfc      IF(FILL(EAST) .OR. FILL(WEST)) THEN
cgsfc        LAST = SIZE(A,1)-1
cgsfc        LEN2 = 1
cgsfc        STR = SIZE(A,1)
cgsfc        LEN1 = SIZE(A)/STR
cgsfc        
cgsfc!       Fill EASTERN ghost region
cgsfc        
cgsfc        IF(FILL(EAST)) THEN
cgsfc          CALL COMM_BSR(layout, A(last+1,:  ),
cgsfc     &         A(2   ,:  ), NN_WEST, LEN1,LEN2,STR,STR)
cgsfc        ENDIF
cgsfc
cgsfc! Fill WESTERN ghost region
cgsfc
cgsfc        IF(FILL(WEST)) THEN
cgsfc          CALL COMM_BSR( layout, A(1   ,:  ), 
cgsfc     &         A(last ,:  ), NN_EAST, LEN1,LEN2,STR,STR)
cgsfc        ENDIF
cgsfc
cgsfc      ENDIF

      RETURN

#undef FILL

      end SUBROUTINE GHOST_R8_2D

      SUBROUTINE GHOST_R8_3D(GRID,A,FROM)
      TYPE(ESMF_Grid),     INTENT(IN)  :: GRID
      REAL*8,         INTENT(INOUT) :: A(:,:,:)
      INTEGER, OPTIONAL,    INTENT(IN)  :: FROM

      integer :: L

      do L = 1, size(A,3)
        call GHOST(GRID,A(:,:,L),FROM)
      end do
      
      end SUBROUTINE GHOST_R8_3D

      SUBROUTINE GHOST_R8_4D(GRID,A,FROM)
      TYPE(ESMF_Grid),     INTENT(IN)  :: GRID
      REAL*8,          INTENT(INOUT) :: A(:,:,:,:)
      INTEGER, OPTIONAL,    INTENT(IN)  :: FROM

      integer :: L

      do L = 1, size(A,4)
        call GHOST(GRID,A(:,:,:,L),FROM)
      end do

      end SUBROUTINE GHOST_R8_4D

      SUBROUTINE GHOST_PAY_R8_2D(GRID,A,FROM)
      TYPE(ESMF_Grid),     INTENT(IN)    :: GRID
      REAL*8,              INTENT(INOUT) :: A(:,:,:)
      INTEGER, OPTIONAL,   INTENT(IN)    :: FROM
      
      INTEGER :: LAST, LEN1, LEN2, STR
      INTEGER :: USABLE_FROM

      integer :: nn_north, nn_south, nn_east, nn_west
      integer :: j_north, j_south, i_east, i_west
      integer                 :: npx, npy
      integer                 :: px, py
      type (ESMF_DELayout)           :: layout
      integer                 :: status
      REAL*8 :: junk(size(A,1),size(A,2))
      
      call ESMF_GridGet(grid, delayout=layout, rc=status)
      call ESMF_GRID_PE_LAYOUT(grid, npx, npy)
      call ESMF_GRID_MY_PE_LOC(grid,  px,  py) 
      
      j_north = MOD(py  +1,npy)
      j_south = MOD(py+npy-1,npy)
      i_east = MOD(px  +1,npx)
      i_west = MOD(px+npx-1,npx)
      
      nn_north = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, px, j_north)
      nn_south = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, px, j_south)
      
      nn_west = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, i_west, py)
      nn_east = ESMF_GRID_PE_NUM_FM_PE_LOC ( grid, i_east, py)
      
!       The default for FROM is to do all four sides.
      
      IF(.NOT.PRESENT(FROM)) THEN
        USABLE_FROM = ALL
      ELSE
        USABLE_FROM = FROM
      ENDIF


#define FILL(N) IAND(USABLE_FROM,N)==N

!       Fill NORTHERN ghost region

      IF(FILL(NORTH) .OR.
     &  FILL(SOUTH)) THEN
        LAST = SIZE(A,3)-1 
        LEN2 = SIZE(A,1)*SIZE(A,2)
        STR = LEN2*SIZE(A,1)*SIZE(A,2)
        LEN1 = SIZE(A)/STR
 
        IF(FILL(NORTH)) THEN
           IF (py == npy-1) THEN  ! No northern neigbor
              CALL COMM_BSR( layout, junk(:,:),
     &             A(:,:,2    ), NN_SOUTH, LEN1,LEN2,STR,STR)
           ELSE
              CALL COMM_BSR( layout, A(:,:,last+1  ),
     &             A(:,:,2    ), NN_SOUTH, LEN1,LEN2,STR,STR)
           END IF
        ENDIF
          
!       Fill SOUTHERN ghost region
        IF(FILL(SOUTH)) THEN
           IF (py == 0) THEN ! No southern neigbbor
              CALL COMM_BSR( layout, junk(:,:),
     &             A(:,:,last   ), NN_NORTH, LEN1,LEN2,STR,STR)
           ELSE
              CALL COMM_BSR( layout, A(:,:,1 ),
     &             A(:,:,last   ), NN_NORTH, LEN1,LEN2,STR,STR)
           END IF
        END IF

      END IF
cgsfc      IF(FILL(EAST) .OR. FILL(WEST)) THEN
cgsfc        LAST = SIZE(A,1)-1
cgsfc        LEN2 = 1
cgsfc        STR = SIZE(A,1)
cgsfc        LEN1 = SIZE(A)/STR
cgsfc        
cgsfc!       Fill EASTERN ghost region
cgsfc        
cgsfc        IF(FILL(EAST)) THEN
cgsfc          CALL COMM_BSR(layout, A(:,last+1,:  ),
cgsfc     &         A(:,2   ,:  ), NN_WEST, LEN1,LEN2,STR,STR)
cgsfc        ENDIF
cgsfc
cgsfc! Fill WESTERN ghost region
cgsfc
cgsfc        IF(FILL(WEST)) THEN
cgsfc          CALL COMM_BSR( layout, A(:,1   ,:  ), 
cgsfc     &         A(:, last ,:  ), NN_EAST, LEN1,LEN2,STR,STR)
cgsfc        ENDIF
cgsfc
cgsfc      ENDIF

      RETURN

#undef FILL

      end SUBROUTINE GHOST_PAY_R8_2D

      subroutine COMM_BSR_R8_0D ( layout, RBUF, SBUF, DEST, 
     &      NUM_BLOCKS, BLK_LEN, RSTR, SSTR )

      IMPLICIT NONE

      type (ESMF_DELayout)       :: layout
      real*8, INTENT(IN)        :: SBUF
!       real, INTENT(OUT)        :: RBUF
      real*8               :: RBUF
      INTEGER,  INTENT(IN)      :: DEST
      INTEGER,  OPTIONAL, INTENT(IN) :: NUM_BLOCKS, BLK_LEN, RSTR, SSTR

!       local vars  
      integer :: status
      integer :: source = MPI_ANY_SOURCE
      integer :: errcode
      
      
      integer :: SEND_BUF_SIZE
      integer :: RECV_BUF_SIZE
      
      SEND_BUF_SIZE = 1
      RECV_BUF_SIZE = 1
      
      call ESMF_DELayoutBarrier(layout, status)
#if 0
      call ESMF_DELayoutSendRecv(layout, sbuf, rbuf, send_buf_size, 
     &       recv_buf_size, source, dest, rc=status)
      call ESMF_DELayoutBarrier(layout, status)
#else
      
      tag = max(MOD(tag,128),10) + 1

      call MPI_SendRecv                  (
     &      sbuf, send_buf_size, MPI_REAL8,
     &      dest, tag,
     &      rbuf, recv_buf_size, MPI_REAL8,
     &      source, tag,
     &  MPI_COMM_WORLD, status,errcode   )
#endif
      
      RETURN
      
      end subroutine COMM_BSR_R8_0D

      subroutine COMM_BSR_R8_1D ( layout, RBUF, SBUF, DEST, 
     &      NUM_BLOCKS, BLK_LEN, RSTR, SSTR )

      IMPLICIT NONE

      type (ESMF_DELayout)       :: layout
      real*8, INTENT(IN)        :: SBUF(:)
!       real, INTENT(OUT)        :: RBUF(:)
      real*8               :: RBUF(:)
      INTEGER,  INTENT(IN)      :: DEST
      INTEGER,  OPTIONAL, INTENT(IN) :: NUM_BLOCKS, BLK_LEN, RSTR, SSTR

!       local vars  
      integer :: status
      integer :: source = MPI_ANY_SOURCE
      integer :: errcode
      
      
      integer :: SEND_BUF_SIZE
      integer :: RECV_BUF_SIZE
      
      
      SEND_BUF_SIZE = SIZE(SBUF)
      RECV_BUF_SIZE = SIZE(RBUF)
      
      call ESMF_DELayoutBarrier(layout, status)
#if 0
      call ESMF_DELayoutSendRecv(layout, sbuf, rbuf, send_buf_size, 
     &recv_buf_size, source, dest, rc=status)
      call ESMF_DELayoutBarrier(layout, status)
#else
      
      tag = max(MOD(tag,128),10) + 1
      
      call MPI_SendRecv                  (
     &      sbuf, send_buf_size, MPI_REAL8,
     &      dest, tag,
     &      rbuf, recv_buf_size, MPI_REAL8,
     &      source, tag,
     &      MPI_COMM_WORLD, status,errcode   )
#endif
      
      RETURN
      
      end subroutine COMM_BSR_R8_1D
!********************************************************
      
      subroutine COMM_BSR_R8_2D ( layout, RBUF, SBUF, DEST,
     &  NUM_BLOCKS, BLK_LEN, RSTR, SSTR )
      
      IMPLICIT NONE
      
      type (ESMF_DELayout)       :: layout
      real*8, INTENT(IN)        :: SBUF(:,:)
      real*8, INTENT(OUT)        :: RBUF(:,:)
      INTEGER,  INTENT(IN)      :: DEST
      INTEGER,  OPTIONAL, INTENT(IN) :: NUM_BLOCKS, BLK_LEN, RSTR, SSTR
      
!       local vars  
      integer :: status
      integer :: source = MPI_ANY_SOURCE
      integer :: errcode
      
      integer :: SEND_BUF_SIZE
      integer :: RECV_BUF_SIZE
      
      
      SEND_BUF_SIZE = SIZE(SBUF)
      RECV_BUF_SIZE = SIZE(RBUF)
      
      call ESMF_DELayoutBarrier(layout, status)
#if 0
!  call ESMF_DELayoutSendRecv(layout, sbuf, rbuf, send_buf_size, &
!       recv_buf_size, source, dest, rc=status)
!       Fortran interface takes only 1d real arrays - call C wrapper directly  
      
      call c_ESMC_DELayoutSendRecv(layout, sbuf, rbuf, send_buf_size,
     &recv_buf_size, source, dest, ESMF_FLOAT, &
      status)
      call ESMF_DELayoutBarrier(layout, status)
#else
      
      tag = max(MOD(tag,128),10) + 1
      
      call MPI_SendRecv                  (
     &      sbuf, send_buf_size, MPI_REAL8,
     &      dest, tag,
     &      rbuf, recv_buf_size, MPI_REAL8,
     &      source, tag,
     &      MPI_COMM_WORLD, status,errcode   )
#endif
      
      RETURN
      
      end subroutine COMM_BSR_R8_2D

#endif 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ESMF_DELayoutBarrier(layout, rc)
      type (ESMF_DELayout) :: layout
      integer, optional  :: rc
      
      integer :: status

#ifdef USE_ESMF
      call MPI_BARRIER(MPI_COMM_WORLD, status)
      if (present(rc)) then
         rc = ESMF_SUCCESS
      endif
#endif
      
      end subroutine ESMF_DELayoutBarrier

#ifndef USE_ESMF
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

      INTEGER :: ier
      Integer, Allocatable :: lines(:)

#ifdef DEBUG_DECOMP
      CALL LOG_PARALLEL(grid, file, line)
      If (AM_I_ROOT()) Then
         WRITE(CHECKSUM_UNIT,*)'HERE: ',file, line
         CALL SYS_FLUSH(CHECKSUM_UNIT)
       End If
#ifdef USE_ESMF       
       ALLOCATE(lines(npes))
       Call MPI_Allgather(line, 1, MPI_INTEGER, lines, 1, MPI_INTEGER, 
     &      MPI_COMM_WORLD, ier)
       If (Any(lines /= line)) 
     &      call stop_model('HERE: synchronization error -severe.',255)
       Deallocate(lines)
#endif
#endif
#ifdef USE_ESMF
      CALL MPI_Barrier(MPI_COMM_WORLD, ier)
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


      SUBROUTINE GLOBALMAX(grd_dum, val, val_max)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ier

#ifdef USE_ESMF
      CALL MPI_Allreduce(val, val_max, 1, MPI_REAL8, MPI_MAX, 
     &     MPI_COMM_WORLD, ier)
#else
      val_max = val
#endif
        
      END SUBROUTINE


      SUBROUTINE PACK_1D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(grd_dum%JM_WORLD)

      CALL ARRAYGATHER(grd_dum,ARR,ARR_GLOB)
 
      RETURN
      END SUBROUTINE PACK_1D
 

      SUBROUTINE PACK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:)

      CALL ARRAYGATHER(grd_dum,ARR,ARR_GLOB)

      RETURN
      END SUBROUTINE PACK_2D

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

      do l=1,size(arr_glob,3)
       CALL ARRAYGATHER(grd_dum,ARR(:,:,l),ARR_GLOB(:,:,l))
      end do

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

      do m=1,size(arr_glob,4)
        do l=1,size(arr_glob,3)
          CALL ARRAYGATHER(grd_dum,ARR(:,:,l,m),ARR_GLOB(:,:,l,m))
        end do
      end do

      RETURN
      END SUBROUTINE PACK_4D

      SUBROUTINE PACKj_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:)

      INTEGER :: k

      DO k = 1, SIZE(ARR,2)
        CALL ARRAYGATHER(grd_dum,ARR(:,k),ARR_GLOB(:,k))
      END DO

      RETURN
      END SUBROUTINE PACKj_2D

      SUBROUTINE PACKj_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)

      INTEGER :: k

      DO k = 1, SIZE(ARR,3)
        CALL PACKj_2D(grd_dum,ARR(:,:,k),ARR_GLOB(:,:,k))
      END DO

      RETURN
      END SUBROUTINE PACKj_3D

      SUBROUTINE PACKj_4D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:,:)

      INTEGER :: k

      DO k = 1, SIZE(ARR,4)
        CALL PACKj_3D(grd_dum,ARR(:,:,:,k),ARR_GLOB(:,:,:,k))
      END DO

      RETURN
      END SUBROUTINE PACKj_4D

      SUBROUTINE UNPACK_1D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(grd_dum%JM_WORLD)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H
      LOGICAL, OPTIONAL :: local
 
      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(J_0H:J_1H)=ARR_GLOB(J_0H:J_1H)
        else
          call ARRAYSCATTER(grd_dum, ARR, ARR_GLOB)
        end if
      else
        call ARRAYSCATTER(grd_dum, ARR, ARR_GLOB)
      end if

      RETURN
      END SUBROUTINE UNPACK_1D
      

      SUBROUTINE UNPACK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,J_0H:J_1H)=ARR_GLOB(:,J_0H:J_1H)
        else
          call arrayscatter(grd_dum, arr, arr_glob)
        end if
      else
        call arrayscatter(grd_dum, arr, arr_glob)
      end if

      RETURN 
      END SUBROUTINE UNPACK_2D

      SUBROUTINE LUNPACK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      LOGICAL, INTENT(IN) :: ARR_GLOB(:,:)
      LOGICAL, INTENT(OUT) ::
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,J_0H:J_1H)=ARR_GLOB(:,J_0H:J_1H)
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
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0H, J_1H, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,J_0H:J_1H,:)=ARR_GLOB(:,J_0H:J_1H,:)
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
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0H, J_1H, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,J_0H:J_1H,:)=ARR_GLOB(:,J_0H:J_1H,:)
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
     &        ARR(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:,:)
      INTEGER :: J_0H, J_1H, L,M
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,J_0H:J_1H,:,:)=ARR_GLOB(:,J_0H:J_1H,:,:)
        else
          do m=1,size(arr,4)
            do l=1,size(arr,3)
              call arrayscatter(grd_dum,arr(:,:,l,m), arr_glob(:,:,l,m))
            end do
          end do
        end if
      else
        do m=1,size(arr,4)
          do l=1,size(arr,3)
            call arrayscatter(grd_dum, arr(:,:,l,m), arr_glob(:,:,l,m))
          end do
        end do
      end if

      RETURN
      END SUBROUTINE UNPACK_4D

      SUBROUTINE UNPACKj_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      INTEGER :: J_0H, J_1H
      LOGICAL, OPTIONAL :: local
      
      LOGICAL :: local_
      INTEGER :: k
      local_ = .true.
      if (present(local)) local_ = local

      if (local_) then
        J_0H=grd_dum%j_strt_halo
        J_1H=grd_dum%j_stop_halo
        ARR(J_0H:J_1H,:)=ARR_GLOB(J_0H:J_1H,:)
      else
        DO k = 1, SIZE(arr,2)
          call arrayscatter(grd_dum, arr(:,k), arr_glob(:,k))
        END DO
      end if

      RETURN 
      END SUBROUTINE UNPACKj_2D

      SUBROUTINE UNPACKj_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      INTEGER :: J_0H, J_1H
      LOGICAL, OPTIONAL :: local
      
      INTEGER :: k

      DO k = 1, SIZE(arr,3)
        call UNPACKj_2D(grd_dum, arr_glob(:,:,k), arr(:,:,k), local)
      END DO

      RETURN 
      END SUBROUTINE UNPACKj_3D

      SUBROUTINE UNPACKj_4D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(grd_dum%j_strt_halo:,:,:,:)
      INTEGER :: J_0H, J_1H
      LOGICAL, OPTIONAL :: local
      
      INTEGER :: k

      DO k = 1, SIZE(arr,4)
        call UNPACKj_3D(grd_dum, arr_glob(:,:,:,k), arr(:,:,:,k), local)
      END DO

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

      DO K=1,SIZE(ARR,1)
       CALL ARRAYGATHER(grd_dum,ARR(k,:),ARR_GLOB(k,:) )
      END DO
 
      RETURN
      END SUBROUTINE PACK_COLUMN_1D
 

      SUBROUTINE PACK_COLUMN_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)
      INTEGER :: K

      DO K=1,SIZE(ARR,1)
      CALL ARRAYGATHER(grd_dum,ARR(K,:,:),ARR_GLOB(K,:,:))
      END DO

      RETURN
      END SUBROUTINE PACK_COLUMN_2D

      SUBROUTINE PACK_COLUMN_i2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
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
     &        ARR(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(OUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: l,K

      DO K=1,SIZE(ARR,1)
        do l=1,size(arr,4)
          CALL ARRAYGATHER(grd_dum,ARR(K,:,:,l),ARR_GLOB(K,:,:,l))
        end do
      END DO
      
      RETURN
      END SUBROUTINE PACK_COLUMN_3D

C*** Overloaded UNPACK_COLUMN routines:
C--------------------------------------
      SUBROUTINE UNPACK_COLUMN_1D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
      REAL*8, INTENT(OUT) ::
     &        ARR(:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H, K
      LOGICAL, OPTIONAL :: local
 
      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,J_0H:J_1H)=ARR_GLOB(:,J_0H:J_1H)
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
     &        ARR(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H, K
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,:,J_0H:J_1H)=ARR_GLOB(:,:,J_0H:J_1H)
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
      END SUBROUTINE UNPACK_COLUMN_2D

      SUBROUTINE IUNPACK_COLUMN_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) :: ARR_GLOB(:,:,:)
      INTEGER, INTENT(OUT) ::
     &        ARR(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H, K
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,:,J_0H:J_1H)=ARR_GLOB(:,:,J_0H:J_1H)
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
     &        ARR(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0H, J_1H, K, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,:,J_0H:J_1H,:)=ARR_GLOB(:,:,J_0H:J_1H,:)
        else
          do k=1,size(arr,1)
            do l=1,size(arr,4)
             call arrayscatter(grd_dum, arr(k,:,:,l),arr_glob(k,:,:,l))
            end do
          end do
        end if
      else
        do k=1,size(arr,1)
          do l=1,size(arr,4)
            call arrayscatter(grd_dum, arr(k,:,:,l), arr_glob(k,:,:,l))
          end do
        end do
      end if

      RETURN
      END SUBROUTINE UNPACK_COLUMN_3D


      SUBROUTINE IPACK_BLOCK_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) ::
     &        ARR(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
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
     &        ARR(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      REAL*8 , INTENT(INOUT) :: ARR_GLOB(:,:,:,:)
      INTEGER :: K,L

      DO K=1,SIZE(ARR,1)
        DO L=1,SIZE(ARR,2)
          CALL ARRAYGATHER(grd_dum,ARR(K,L,:,:),ARR_GLOB(K,L,:,:))
        END DO
      END DO

      RETURN
      END SUBROUTINE PACK_BLOCK_2D


      SUBROUTINE PACK_BLOCK_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8 , INTENT(IN) ::
     &        ARR(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      REAL*8 , INTENT(INOUT) :: ARR_GLOB(:,:,:,:,:)
      INTEGER :: K,L,M

      DO M=1,SIZE(ARR,5) 
        DO L=1,SIZE(ARR,2)
          DO K=1,SIZE(ARR,1)
            CALL ARRAYGATHER(grd_dum,ARR(K,L,:,:,M),ARR_GLOB(K,L,:,:,M))
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE PACK_BLOCK_3D


      SUBROUTINE IUNPACK_BLOCK_2D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      INTEGER, INTENT(IN) :: ARR_GLOB(:,:,:,:)
      INTEGER, INTENT(OUT) ::
     &        ARR(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H, K, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,:,:,J_0H:J_1H)=ARR_GLOB(:,:,:,J_0H:J_1H)
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

      REAL*8 , INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8 , INTENT(OUT) ::
     &        ARR(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER :: J_0H, J_1H, K, L
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,:,:,J_0H:J_1H)=ARR_GLOB(:,:,:,J_0H:J_1H)
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
      END SUBROUTINE UNPACK_BLOCK_2D

      SUBROUTINE UNPACK_BLOCK_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8 , INTENT(IN) :: ARR_GLOB(:,:,:,:,:)
      REAL*8 , INTENT(OUT) ::
     &        ARR(:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0H, J_1H, K, L, M
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
          ARR(:,:,:,J_0H:J_1H,:)=ARR_GLOB(:,:,:,J_0H:J_1H,:)
        else
          do m=1,size(arr,5)
            do l=1,size(arr,2)
              do k=1,size(arr,1)
                call arrayscatter(grd_dum,arr(k,l,:,:,m)
     &                              ,arr_glob(k,l,:,:,m))
              end do
            end do
          end do
        end if
      else
        do m=1,size(arr,5)
          do l=1,size(arr,2)
            do k=1,size(arr,1)
              call arrayscatter(grd_dum, arr(k,l,:,:,m)
     &                            , arr_glob(k,l,:,:,m))
            end do
          end do
        end do
      end if
      END SUBROUTINE UNPACK_BLOCK_3D


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
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
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
          J_0H=grd_dum%j_strt_halo
          J_1H=grd_dum%j_stop_halo
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

      SUBROUTINE ESMF_BCAST_1D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:)

      INTEGER :: ier

#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_REAL8,root,
     &     MPI_COMM_WORLD, ier)
#endif

      END SUBROUTINE ESMF_BCAST_1D


      END MODULE DOMAIN_DECOMP
