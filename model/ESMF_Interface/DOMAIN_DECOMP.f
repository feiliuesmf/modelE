#include "mpi_defs.h"

      MODULE DOMAIN_DECOMP
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@+    for the message passing (ESMF) implementation.
!@auth NCCS ASTG

#define FILL(N) IAND(USABLE_FROM,N)==N

#ifdef USE_ESMF
      use ESMF_Mod
#endif
      use ESMF_CUSTOM_MOD, Only: NORTH, SOUTH
      IMPLICIT NONE
#ifdef USE_ESMF
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
      PUBLIC :: grid, grid_TRANS
!@var INIT_APP Initialize default decomposition
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
      PUBLIC :: DREAD_PARALLEL
      PUBLIC :: MREAD_PARALLEL
      PUBLIC :: READT_PARALLEL
      PUBLIC :: READ_PARALLEL
      PUBLIC :: WRITE_PARALLEL
      PUBLIC :: WRITEI_PARALLEL
      PUBLIC :: TRANSP
      PUBLIC :: TRANSPOSE_COLUMN
!@var GLOBALMAX Generic wrapper for Real/integer
      INTERFACE GLOBALMAX
        MODULE PROCEDURE GLOBALMAX_R
        MODULE PROCEDURE GLOBALMAX_I
      END INTERFACE

!@var HALO_UPDATE Generic wrapper for 2D and 3D routines
      INTERFACE HALO_UPDATE
        MODULE PROCEDURE HALO_UPDATE_1D  ! J
        MODULE PROCEDURE HALO_UPDATE_2D  ! I,J
        MODULE PROCEDURE HALO_UPDATE_3D  ! I,J,K
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
      END INTERFACE

      INTERFACE ARRAYSCATTER
        MODULE PROCEDURE ARRAYSCATTER_J
        MODULE PROCEDURE ARRAYSCATTER_IJ
        MODULE PROCEDURE IARRAYSCATTER_IJ
        MODULE PROCEDURE LARRAYSCATTER_IJ
      END INTERFACE

#ifdef USE_ESMF
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

      interface WRITEI_PARALLEL
        module procedure WRITEI_PARALLEL_2D
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
      PUBLIC :: NORTH, SOUTH

!, EAST, WEST
      ! Use powers of two such that addition can be used
      ! to create combination
c***      INTEGER, PARAMETER :: NORTH = 2**0, SOUTH = 2**1
c***      INTEGER, PARAMETER :: EAST  = 2**2, WEST  = 2**3
      INTEGER, PARAMETER :: ALL = NORTH + SOUTH ! no east/west for now


      INTEGER, PARAMETER :: HALO_WIDTH = 1
      integer ::  root

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
#ifdef DEBUG_DECOMP
         INTEGER :: log_unit ! for debugging
#endif
      END TYPE DIST_GRID

      TYPE (DIST_GRID) :: GRID, GRID_TRANS

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
      INTEGER, PUBLIC :: CHECKSUM_UNIT

      Integer :: tag = 10

      CONTAINS

      ! This routine initializes the quantities described above.
      ! The initialization should proceed prior to any grid computations.
      SUBROUTINE INIT_APP(grd_dum,IM,JM,LM)
      USE FILEMANAGER, Only : openunit
      USE ESMF_CUSTOM_MOD, Only: Initialize_App
#ifdef USE_ESMF
      USE ESMF_CUSTOM_MOD, Only: vm => modelE_vm
#endif
!AOO      USE ESMF_CUSTOM_MOD, Only: modelE_grid
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM, LM
      INTEGER             :: rc
      INTEGER             :: pet
      CHARACTER(LEN=20) :: buffer

#ifdef USE_ESMF
#endif

#ifdef USE_ESMF
      ! Initialize ESMF
      Call Initialize_App(IM, JM, LM,rc=rc)

      Call ESMF_VMGet(vm, localPET = my_pet, petCount = NPES, rc=rc)
      root = ROOT_ID
      compmodelE  = ESMF_GridCompCreate(vm,"ModelE ESMF", rc=rc)

      ! The default layout is not what we want - it splits in the "I" direction.
      ESMF_Layout = ESMF_DELayoutCreate(vm, deCountList = (/ 1, NPES /))

      NP_LON = 1
      NP_LAT = NPES
      RANK_LAT = my_pet
      RANK_LON = 0
#else
      NPES = 1
      NP_LON = 1
      NP_LAT = 1
      RANK_LON = 0
      RANK_LAT = 0
#endif

#ifdef USE_ESMF
      call INIT_GRID(grd_dum,IM,JM,LM,vm=vm)
      Call ESMF_GridCompSet(compmodelE, grid=grd_dum%ESMF_GRID, rc=rc)
      call INIT_GRID(grid_TRANS,JM,IM,LM,width=0,vm=vm)
#else
      call INIT_GRID(grd_dum,IM,JM,LM)
      call INIT_GRID(grid_TRANS,JM,IM,LM,width=0)
#endif

      WRITE(*,*)'Domain Decomposition for rank: ',MY_PET,RANK_LAT,
     &     RANK_LON

#ifdef DEBUG_DECOMP
      IF (AM_I_ROOT()) CALL openunit('CHKSUM_DECOMP', CHECKSUM_UNIT)
      WRITE(buffer,'(a,i3.3)') 'LOG_',my_pet
      CALL openunit(TRIM(buffer), grd_dum%log_unit)
#endif

      END SUBROUTINE INIT_APP

      SUBROUTINE DESTROY_GRID(grd_dum)
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
#ifdef USE_ESMF
      Call ESMF_GridDestroy(grd_dum%ESMF_Grid)
#endif
      END SUBROUTINE DESTROY_GRID

#ifdef USE_ESMF
      SUBROUTINE INIT_GRID(grd_dum,IM,JM, LM,width,vm)
      USE ESMF_CUSTOM_MOD, Only : modelE_vm
#else
      SUBROUTINE INIT_GRID(grd_dum,IM,JM,LM,width)
#endif
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM,LM
      INTEGER, OPTIONAL :: width
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
#ifdef USE_ESMF
      TYPE(ESMF_VM), Pointer :: vm_
      Type (ESMF_DELayout)::layout
      REAL*8 :: deltaZ
      INTEGER :: L
#endif
      grid_size(1)=IM;   grid_size(2)=JM
      range_min(1)=0.;   range_min(2)=-90.
      range_max(1)=360.; range_max(2)=90.

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
      Call ESMF_GridDistribute(grid=grd_dum%ESMF_GRID,
     &     delayout = layout, rc=rc)
      call ESMF_GridGet(grd_dum%esmf_grid, delayout=layout, rc=rc)
      RANK_LON=0
      RANK_LAT=my_pet
      Call ESMF_GRID_BOUNDS(grd_dum, RANK_LON, RANK_LAT,
     &        I0_DUM, I1_DUM, J0_DUM, J1_DUM)

#else
      RANK_LON = 0
      RANK_LAT = 0
      I0_DUM = 1
      I1_DUM = IM
      J0_DUM = 1
      J1_DUM = JM
#endif

      width_ = HALO_WIDTH
      If (Present(width)) width_=width
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

#ifdef USE_ESMF
      grd_dum%J_STRT_HALO   = J0_DUM - width_
      grd_dum%J_STOP_HALO   = J1_DUM + width_
#else
      grd_dum%J_STRT_HALO   = 1
      grd_dum%J_STOP_HALO   = JM
#endif

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
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 1, from)
#endif

      END SUBROUTINE HALO_UPDATE_1D

      SUBROUTINE HALO_UPDATE_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                    arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_ESMF
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 2, from)
#endif
      END SUBROUTINE HALO_UPDATE_2D

      SUBROUTINE HALO_UPDATEj_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                    arr(grd_dum%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_ESMF
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 1, from)
#endif
      END SUBROUTINE HALO_UPDATEj_2D

      SUBROUTINE HALO_UPDATE_3D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                 arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_ESMF

      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 2, from)
#endif
      END SUBROUTINE HALO_UPDATE_3D

      SUBROUTINE HALO_UPDATE_COLUMN_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &                  arr(:,grd_dum%j_strt_halo:)

      INTEGER, OPTIONAL, INTENT(IN)    :: from

#ifdef USE_ESMF
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 2, from)
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
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 3, from)
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
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 3, from)
#endif
      END SUBROUTINE HALO_UPDATE_COLUMN_4D

      SUBROUTINE HALO_UPDATE_COLUMN_7D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DIST_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) ::
     &     arr(:,:,:,:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      INTEGER :: L

#ifdef USE_ESMF
      Call sendrecv(grd_dum%ESMF_GRID, arr, shape(arr), 6, from)
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
      CALL closeunit(grid%log_unit)
#endif

#ifdef USE_ESMF
      CALL ESMF_FINALIZE(rc=ier)
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
#ifdef USE_ESMF
        CALL MPI_ABORT(ier)
#else
        STOP
#endif
      End If
#endif

#ifdef USE_ESMF
      Call gather(grd_dum%ESMF_GRID, arr, garr, shape(arr), 1)
#else
      garr = arr(j_0:j_1)
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

#ifdef USE_ESMF
      If (Present(all)) Then
         If (all) THEN
            Call MPI_BCAST(gsum,1,MPI_DOUBLE_PRECISION,root,
     &           MPI_COMM_WORLD, ier)
            If (Present(hsum)) Call MPI_BCAST(hsum,2,
     &           MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, ier)
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
      Call gather(grd_dum%ESMF_GRID, zon, garr, shape(zon), 1)
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
            Call MPI_BCAST(gsum,1,MPI_DOUBLE_PRECISION,root,
     &           MPI_COMM_WORLD, ier)
            If (Present(hsum))     Call MPI_BCAST(hsum,2,
     &           MPI_DOUBLE_PRECISION,root, MPI_COMM_WORLD, ier)
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
      Call gather(grd_dum%ESMF_GRID, zon, garr, shape(zon), 1)
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
      INTEGER :: ier
#ifdef USE_ESMF
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world)
#endif
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
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

#ifdef USE_ESMF
      Call gather(grd_dum%ESMF_GRID, arr, garr, shape(arr), 2)
      IF (AM_I_ROOT()) gsum = Sum(garr(:,jb1:jb2),2)
      If (all_) Then
         call MPI_BCAST(gsum, Size(gsum), MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, ier)
      End If
#else
      gsum = Sum(arr(:,jb1:jb2),2)
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
      INTEGER :: ier
#ifdef USE_ESMF
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world,size(arr,3))
#endif
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
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

#ifdef USE_ESMF
      Call gather(grd_dum%ESMF_GRID, arr, garr, shape(arr), 2)
      IF (AM_I_ROOT()) gsum = Sum(garr(:,jb1:jb2,:),2)
      If (all_) Then
         call MPI_BCAST(gsum, Size(gsum), MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, ier)
      End If
#else
      gsum = Sum(arr(:,jb1:jb2,:),2)
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
      INTEGER :: ier
#ifdef USE_ESMF
      REAL*8  :: garr(size(arr,1),grd_dum%jm_world,size(arr,3))
#endif
    ! now local
#ifdef USE_ESMF
      type (ESMF_Grid)                           :: GRID
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

#ifdef USE_ESMF
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
     &             mpi_comm_world, ier)

      tsum=sum(recv_buf,2)

      rcnts=dik_map
      rdspl(0)=0
      Do p = 1, npes-1
         rdspl(p)=rdspl(p-1)+rcnts(p-1)
      End Do

      Call MPI_GatherV(tsum, dik, mpi_double_precision,
     & gsum, dik_map, rdspl, mpi_double_precision,
     & root, mpi_comm_world, ier)

      Deallocate(recv_buf)
      Deallocate(send_buf)
      Deallocate(tsum)

      if (all_) Then
         call MPI_BCAST(gsum, Size(gsum), MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, ier)
      End If
#else
      gsum = Sum(arr(:,1:JM,:),2)
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
      Call gather(grd_dum%ESMF_GRID, arr, garr, shape(arr), 1)
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
            Call MPI_BCAST(gsum,Size(gsum),MPI_DOUBLE_PRECISION,root,
     &           MPI_COMM_WORLD, ier)
            If (Present(hsum)) Call MPI_BCAST(hsum,size(hsum),
     &           MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ier)
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
      INTEGER :: ier
      INTEGER :: i_0, i_1, j_0, j_1, IM, JM
      REAL*8  :: garr(size(arr,1),size(arr,2),size(arr,3),
     &     grd_dum%jm_world)
      LOGICAL :: istag_,all_

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

#ifdef USE_ESMF
      Call gather(grd_dum%ESMF_GRID, arr, garr, shape(arr), 4, all=all)
      all_=.false.
      If (Present(all)) all_=all
      If (AM_I_ROOT() .or. all_)  gsum = SUM(garr,4)
#else
      gsum = SUM(arr(:,:,:,J_0:J_1),4)
#endif

      END SUBROUTINE GLOBALSUM_XXXJ_XXX

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
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
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

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,1:grd_dum%JM_WORLD)=AOUT
#endif
      if (AM_I_ROOT()) then
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('DREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE DREAD_PARALLEL_2D

      SUBROUTINE DREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR)
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
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
      Call scatter(grd_dum%ESMF_GRID, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=AOUT
#endif

      if (AM_I_ROOT()) then
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('DREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE DREAD_PARALLEL_3D

      SUBROUTINE MREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,M,NSKIP,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME      !@var  NAME  name of record being read
      INTEGER,      INTENT(OUT) :: M         !@var  M      initial integer
      INTEGER,      INTENT(IN)  :: NSKIP     !@var  NSKIP no. of R*4's to skip
      REAL*8,      INTENT(OUT)  :: AVAR(:,grd_dum%J_STRT_HALO:) !@var  AOUT real*8 array
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD)  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM
      INTEGER :: IERR
    ! now local
      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) M,(X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, AOUT, AVAR, shape(AVAR), 2)
      CALL ESMF_BCAST(grd_dum, M   )
#else
      AVAR=AOUT
#endif

      if (AM_I_ROOT()) then
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('MREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE MREAD_PARALLEL_2D

      SUBROUTINE MREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,M,NSKIP,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME        !@var  NAME  name of record being read
      INTEGER,      INTENT(OUT) :: M           !@var  M      initial integer
      INTEGER,      INTENT(IN)  :: NSKIP       !@var  NSKIP no. of R*4's to skip
      REAL*8,      INTENT(OUT)  :: AVAR(:,grd_dum%J_STRT_HALO:,:) !@var  AOUT real*8 array
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4 :: AIN(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3))  !@var  AIN  real*4 array
      REAL*8 :: AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3)) !@var  AOUT real*8 array
      INTEGER :: N                        !@var  N loop variable
      INTEGER :: IM,JM,NM
      INTEGER :: IERR
    ! now local

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      NM   = size(AVAR,3)

      If (AM_I_ROOT()) then
         READ (IUNIT,IOSTAT=IERR) M,(X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, AOUT, AVAR, shape(AVAR), 2)
      CALL ESMF_BCAST(grd_dum, M   )
#else
      AVAR=AOUT
#endif

      if (AM_I_ROOT()) then
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('MREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE MREAD_PARALLEL_3D

      SUBROUTINE READT_PARALLEL_2D (grd_dum,IUNIT,NAME,NSKIP,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
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
      Call scatter(grd_dum%ESMF_GRID, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP)=AOUT
#endif

      if (AM_I_ROOT()) then
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': ',
     &           TRIM(TITLE),' IOSTAT=',IERR
            call stop_model('READT_PARALLEL: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT_PARALLEL_2D

      SUBROUTINE READT_PARALLEL_3D (grd_dum,IUNIT,NAME,NSKIP,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
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
      EndIf

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, AOUT, AVAR, shape(AVAR), 2)
#else
      AVAR(:,grd_dum%J_STRT:grd_dum%J_STOP,:)=AOUT
#endif

      if (am_i_root()) then
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ':
     &           ',TRIM(TITLE),' IOSTAT=',IERR
            call stop_model('READT_PARALLEL: READ ERROR',255)
         EndIf
      end if

      END SUBROUTINE READT_PARALLEL_3D


      SUBROUTINE WRITEI_PARALLEL_2D (grd_dum,IUNIT,NAME,buf,it)
!@sum WRITEI_PARALLEL  Parallel version of UTILDBL.f:WRITEI for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER*16, INTENT(IN)  :: NAME       !@var  NAME  name of record being read
      REAL*4,       INTENT(IN) :: buf(:,grd_dum%J_STRT_HALO:)  !@var  buf real*8 array
      INTEGER,      INTENT(IN)  :: it       !@var  it iteration
!@var buf_glob real*4 array
      REAL*4 :: buf_glob(grd_dum%IM_WORLD,grd_dum%JM_WORLD) 
      INTEGER :: IERR

!!! not sure if it is implemented for real*4 ...
#ifdef USE_ESMF
c***      Call gather(grd_dum%ESMF_GRID, buf, buf_glob, shape(buf), 2)
      buf_glob = -99999999. ! not implemented
#else
      buf_glob = buf(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      If (AM_I_ROOT()) then
        WRITE (IUNIT, IOSTAT=IERR) it, buf_glob, it
         If (IERR==0) Then
            WRITE(6,*) "Wrote to file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'WRITE ERROR ON FILE ', NAME, ' IOSTAT=',IERR
            call stop_model('WRITEI_PARALLEL: WRITE ERROR',255)
         EndIf
      end if

      END SUBROUTINE WRITEI_PARALLEL_2D


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

      Allocate(AI(1:NPES,3))
      Call ESMF_GridGetAllAxisIndex(egrid, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)

      allocate (sendcounts(NPES), displs(0:NPES), stat=status)

      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
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
     &     MPI_INTEGER         , root, MPI_COMM_WORLD, status)

        deallocate(VAR, stat=status)

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

      integer                                       :: I, J, II, JJ, III
      integer                                       :: NX, NY
      integer                                       :: I1, IN
      integer                                       :: J1, JN


      Logical      , allocatable                    :: var(:)

      Allocate(AI(1:NPES,3))
      Call ESMF_GridGetAllAxisIndex(egrid, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)
      allocate (sendcounts(NPES), displs(0:NPES), stat=status)

      if (AM_I_ROOT()) then
        allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
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
     &     MPI_LOGICAL         , root, MPI_COMM_WORLD, status)

        deallocate(VAR, stat=status)

      deallocate(sendcounts, displs, stat=status)
      deallocate(AI)
#else
      local_array=global_array
#endif
      end subroutine ESMF_LArrayScatter_IJ

#ifdef USE_ESMF
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


      Allocate(AI(1:NPES,3))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)

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
     &         shape=(/IN-I1+1,
     &         JN-J1+1/))
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

      Allocate(AI(1:NPES,3))
      Call ESMF_GridGetAllAxisIndex(e_grid, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)
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
     &         shape=(/IN-I1+1,
     &         JN-J1+1/))
        enddo
      endif

      deallocate(VAR, stat=status)
      deallocate(recvcounts, displs, stat=status)
      Deallocate(AI)
      end subroutine Esmf_LArrayGather_IJ
!------------------------------------------------------------

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


      Allocate(AI(1:NPES,3))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)

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
      Call Gather(grd_dum%ESMF_GRID, local_array, global_array,
     &     shape(local_array), 1)
#else
      global_array = local_array(grd_dum%J_STRT:grd_dum%J_STOP)
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
      global_array = local_array(grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine ArrayGather_J_int

!---------------------------

      subroutine ArrayGather_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:,:)               :: global_array

#ifdef USE_ESMF
      Call gather(grd_dum%ESMF_GRID, local_array, global_array,
     &     shape(local_array), 2)

#else
      global_array = local_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
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
      global_array = local_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
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
      global_array = local_array(:,grd_dum%J_STRT:grd_dum%J_STOP)
#endif

      end subroutine IArrayGather_IJ


    !---------------------------
!---------------------------

      subroutine ArrayScatter_J(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:)               :: global_array

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, global_array, local_array,
     &     shape(local_array), 1)
#else
      local_array(grd_dum%J_STRT:grd_dum%J_STOP) = global_array
#endif

      end subroutine ArrayScatter_J

!---------------------------

      subroutine ArrayScatter_IJ(grd_dum, local_array, global_array)
      type (DIST_GRID)      :: grd_dum
      real (kind=8), dimension (:,grd_dum%J_STRT_HALO:) :: local_array
      real (kind=8), dimension (:,:)               :: global_array

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, global_array, local_array,
     &     shape(local_array), 2)
#else
      local_array(:,grd_dum%J_STRT:grd_dum%J_STOP) = global_array
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
      local_array(:,grd_dum%J_STRT:grd_dum%J_STOP) = global_array
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
      local_array(:,grd_dum%J_STRT:grd_dum%J_STOP) = global_array
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
        integer :: status

        type(ESMF_AxisIndex), dimension(:,:), pointer :: AI

        Allocate(AI(1:NPES,3))
        call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, AI,
     &       horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)
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

    !---------------------------

#ifdef USE_ESMF
      subroutine ESMF_GRID_WORLD(GRID,IM_WORLD,JM_WORLD)
        type (ESMF_Grid), intent(IN) :: grid
        integer, intent(OUT)         :: IM_WORLD, JM_WORLD

! local vars
        integer status
        integer, dimension(ESMF_MAXGRIDDIM) :: dims

        call ESMF_GridGet(grid, horzRelLoc=ESMF_CELL_CENTER,
     &     globalcellcountperdim=dims, rc=status)

        IM_WORLD = dims(1)
        JM_WORLD = dims(2)

      end subroutine ESMF_GRID_WORLD
#endif

!---------------------------

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
      CALL MPI_BARRIER(MPI_COMM_WORLD, ier)
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


      SUBROUTINE GLOBALMAX_R(grd_dum, val, val_max)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      REAL*8,            INTENT(IN)  :: val
      REAL*8,            INTENT(OUT) :: val_max

      INTEGER  :: ier

#ifdef USE_ESMF
      CALL MPI_Allreduce(val, val_max, 1, MPI_DOUBLE_PRECISION,MPI_MAX,
     &     MPI_COMM_WORLD, ier)
#else
      val_max = val
#endif

      END SUBROUTINE

      SUBROUTINE GLOBALMAX_I(grd_dum, val, val_max)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN)  :: grd_dum
      INTEGER,            INTENT(IN)  :: val
      INTEGER,            INTENT(OUT) :: val_max

      INTEGER  :: ier

#ifdef USE_ESMF
      CALL MPI_Allreduce(val, val_max, 1, MPI_INTEGER, MPI_MAX,
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

#ifdef USE_ESMF
      Call Gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 1)
#else
      arr_glob = arr(grd_dum%J_STRT:grd_dum%J_STOP)
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

#ifdef USE_ESMF
      Call Gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 2)
#else
      arr_glob = arr(:,grd_dum%J_STRT:grd_dum%J_STOP)
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

#ifdef USE_ESMF
      Call Gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 2)
#else
      arr_glob = arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:)
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

#ifdef USE_ESMF
      Call Gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 2)
#else
      arr_glob = arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:)
#endif

      RETURN
      END SUBROUTINE PACK_4D

      SUBROUTINE PACKj_2D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:)

      INTEGER :: k

#ifdef USE_ESMF
      CALL gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 1)
#else
      arr_glob=arr(grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

      RETURN
      END SUBROUTINE PACKj_2D

      SUBROUTINE PACKj_3D(grd_dum,ARR,ARR_GLOB)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8, INTENT(IN) ::
     &        ARR(grd_dum%j_strt_halo:,:,:)
      REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)

      INTEGER :: k

#ifdef USE_ESMF
      CALL gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 1)
#else
      arr_glob=arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:)
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

#ifdef USE_ESMF
      CALL gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 1)
#else
      arr_glob=arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)
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

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP)=arr_glob
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

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, arr_glob, arr, shape(arr), 2)
#else
      arr(:,grd_dum%J_STRT:grd_dum%J_STOP)=arr_glob
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

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, arr_glob, arr, shape(arr), 2)
#else
      arr(:,grd_dum%J_STRT:grd_dum%J_STOP,:,:)=arr_glob
#endif

      RETURN
      END SUBROUTINE UNPACK_4D

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

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP,:)=arr_glob
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

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:)=arr_glob
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

#ifdef USE_ESMF
      Call scatter(grd_dum%ESMF_GRID, arr_glob, arr, shape(arr), 1)
#else
      arr(grd_dum%J_STRT:grd_dum%J_STOP,:,:,:)=arr_glob
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

#ifdef USE_ESMF
      CALL gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 2)
#else
      arr_glob=arr(:,grd_dum%J_STRT:grd_dum%J_STOP)
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

#ifdef USE_ESMF
      CALL gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 3)
#else
      arr_glob=arr(:,:,grd_dum%J_STRT:grd_dum%J_STOP)
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

#ifdef USE_ESMF
      CALL gather(grd_dum%ESMF_GRID, arr, arr_glob, shape(arr), 3)
#else
      arr_glob=arr(:,:,grd_dum%J_STRT:grd_dum%J_STOP,:)
#endif

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
      END SUBROUTINE UNPACK_COLUMN_2D

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
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,J_0:J_1,:)=ARR_GLOB(:,:,J_0:J_1,:)
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
     &        ARR(:,:,:,grd_dum%j_strt_halo:,:)
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

      REAL*8 , INTENT(IN) :: ARR_GLOB(:,:,:,:)
      REAL*8 , INTENT(OUT) ::
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
      END SUBROUTINE UNPACK_BLOCK_2D

      SUBROUTINE UNPACK_BLOCK_3D(grd_dum,ARR_GLOB,ARR,local)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum

      REAL*8 , INTENT(IN) :: ARR_GLOB(:,:,:,:,:)
      REAL*8 , INTENT(OUT) ::
     &        ARR(:,:,:,grd_dum%j_strt_halo:,:)
      INTEGER :: J_0, J_1, K, L, M
      LOGICAL, OPTIONAL :: local

      if (present(local)) then
        if (local) then
          J_0=grd_dum%j_strt
          J_1=grd_dum%j_stop
          ARR(:,:,:,J_0:J_1,:)=ARR_GLOB(:,:,:,J_0:J_1,:)
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

      SUBROUTINE ESMF_BCAST_0D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr

      INTEGER :: ier

#ifdef USE_ESMF
      Call MPI_BCAST(arr,1,MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ier)
#endif

      END SUBROUTINE ESMF_BCAST_0D

      SUBROUTINE ESMF_BCAST_1D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:)

      INTEGER :: ier

#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ier)
#endif

      END SUBROUTINE ESMF_BCAST_1D

      SUBROUTINE ESMF_BCAST_2D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:)

      INTEGER :: ier

#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ier)
#endif

      END SUBROUTINE ESMF_BCAST_2D

      SUBROUTINE ESMF_BCAST_3D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:,:)

      INTEGER :: ier

#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ier)
#endif

      END SUBROUTINE ESMF_BCAST_3D

      SUBROUTINE ESMF_BCAST_4D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Real*8, Intent(InOut) :: arr(:,:,:,:)

      INTEGER :: ier

#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_DOUBLE_PRECISION,root,
     &     MPI_COMM_WORLD, ier)
#endif

      END SUBROUTINE ESMF_BCAST_4D

      SUBROUTINE ESMF_IBCAST_0D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr
      INTEGER :: ier
#ifdef USE_ESMF
      Call MPI_BCAST(arr,1,MPI_INTEGER,root,
     &     MPI_COMM_WORLD, ier)
#endif
      END SUBROUTINE ESMF_IBCAST_0D

      SUBROUTINE ESMF_IBCAST_1D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:)
      INTEGER :: ier
#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root,
     &     MPI_COMM_WORLD, ier)
#endif
      END SUBROUTINE ESMF_IBCAST_1D

      SUBROUTINE ESMF_IBCAST_2D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:)
      INTEGER :: ier
#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root,
     &     MPI_COMM_WORLD, ier)
#endif
      END SUBROUTINE ESMF_IBCAST_2D

      SUBROUTINE ESMF_IBCAST_3D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:,:)
      INTEGER :: ier
#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER, root,
     &     MPI_COMM_WORLD, ier)
#endif
      END SUBROUTINE ESMF_IBCAST_3D

      SUBROUTINE ESMF_IBCAST_4D(grd_dum, arr)
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(In) :: grd_dum
      Integer, Intent(InOut) :: arr(:,:,:,:)
      INTEGER :: ier
#ifdef USE_ESMF
      Call MPI_BCAST(arr,Size(arr),MPI_INTEGER ,root,
     &     MPI_COMM_WORLD, ier)
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
#ifdef USE_ESMF
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,nk,k
      INTEGER :: status, ier, p
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt
      INTEGER, DIMENSION(0:NPES-1) :: scnts, rcnts, sdspl, rdspl
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_ESMF
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
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)

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
     &        MPI_COMM_WORLD, ier)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ier)
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
#ifdef USE_ESMF
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ
      INTEGER :: status, ier, p
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt
      INTEGER, DIMENSION(0:NPES-1) :: scnts, rcnts, sdspl, rdspl
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_ESMF
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
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)

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
     &        MPI_COMM_WORLD, ier)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ier)
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
#ifdef USE_ESMF
      TYPE (ESMF_AXISINDEX), Pointer :: AI(:,:)
#endif
      INTEGER :: I,J, II,JJ,k
      INTEGER :: status, ier, p
      INTEGER :: ni_loc, nj_loc, nip, njp, icnt
      INTEGER, DIMENSION(0:NPES-1) :: scnts, rcnts, sdspl, rdspl
      INTEGER :: n, nk
      LOGICAL :: reverse_

      reverse_=.false.
      If (PRESENT(reverse)) reverse_=reverse

#ifndef USE_ESMF
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
      Call ESMF_GridGetAllAxisIndex(grid%ESMF_GRID, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=status)

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
     &        MPI_COMM_WORLD, ier)
      Else
         CALL MPI_ALLTOALLV(sbuf, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &        rbuf, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD, ier)
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

#ifdef USE_ESMF
      Function CreateDist_MPI_Type(base_type, counts, dist_idx)
     &     Result(new_type)
      Integer, Intent(In) :: base_type
      Integer, Intent(In) :: counts(:)
      Integer, Intent(In) :: dist_idx
      Integer :: new_type

      Integer :: stride, n_blocks, blocklen
      Integer :: ext_lb
      Integer :: base_byte_len, new_len
      Integer :: vtype1, vtype2
      Integer :: ier

      n_blocks = Product(counts(dist_idx+1:))
      blocklen = Product(counts(:dist_idx-1))
      stride = counts(dist_idx) * blocklen

      Call MPI_Type_vector(n_blocks, blocklen, stride, base_type,
     &     vtype1, ier)
      Call MPI_Type_extent(base_type, base_byte_len, ier)
      new_len = base_byte_len * blocklen
      Call MPI_Type_struct(2, (/ 1, 1 /), (/ 0, new_len /),
     &     (/ vtype1, MPI_UB /), vtype2, ier)
      Call MPI_Type_Free(vtype1, ier)

      Call MPI_Type_Commit(vtype2, ier)
      new_type = vtype2

      End Function CreateDist_MPI_Type

      Subroutine GetNeighbors(rank, npes, pe_south, pe_north)
      Integer, Intent(In)  :: rank
      Integer, Intent(In)  :: npes
      Integer, Intent(Out) :: pe_south
      Integer, Intent(Out) :: pe_north

      If (rank > 0) Then
        pe_south = rank - 1
      Else
        pe_south = MPI_PROC_NULL
      End If

      If (rank < npes-1) Then
        pe_north = rank + 1
      Else
        pe_north = MPI_PROC_NULL
      End If

      End Subroutine GetNeighbors

      Subroutine sendrecv(grid, arr, shp, dist_idx, from)
      Type (Esmf_Grid) :: grid
      Real(Kind=8) :: arr(*)
      Integer :: shp(:)
      Integer :: dist_idx
      Integer, optional :: from

      Integer :: new_type
      Integer :: npy, npx, px, py, pe_south, pe_north
      Integer :: off(4)
      Integer :: USABLE_FROM
      Integer :: status(MPI_STATUS_SIZE), ier
      Integer :: n, sz

      IF(.NOT.PRESENT(FROM)) THEN
        USABLE_FROM = ALL
      ELSE
        USABLE_FROM = FROM
      ENDIF

      ! create a new mpi type for use in communication
      !-------------------------------
      new_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION, shp,dist_idx)

      ! Determine neigboring processes
      !-------------------------------
      call ESMF_GRID_MY_PE_LOC(grid,  px,  py)
      call ESMF_GRID_PE_LAYOUT(grid, npx, npy)
      Call GetNeighbors(py, npy, pe_south, pe_north)

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
     &                    MPI_COMM_WORLD, status, ier)
      End If

      IF(FILL(SOUTH)) THEN
        tag = max(MOD(tag,128),10) + 1
        Call MPI_SendRecv(arr(off(3)), 1, new_type, pe_north, tag,
     &                    arr(off(1)), 1, new_type, pe_south, tag,
     &                    MPI_COMM_WORLD, status, ier)
      End If

      Call MPI_Type_Free(new_type, ier)

      End SUBROUTINE SendRecv

      Subroutine gather(grid, arr_loc, arr_glob, shp, dist_idx, all)
      Type (Esmf_Grid) :: grid
      Real(Kind=8) :: arr_loc(*)
      Real(Kind=8) :: arr_glob(*)
      Integer :: shp(:),shp_glob(size(shp))
      Integer :: dist_idx
      Logical, Optional :: all

      Logical :: all_
      Integer :: new_type, orig_type
      Integer :: ier
      Integer :: p, scount, offset
      Integer, Allocatable :: rcounts(:), displs(:)
      Type (ESMF_Axisindex), Pointer :: AI(:,:)
      Integer :: i, n_its = 1
      ! create a new mpi type for use in communication
      !-------------------------------
      orig_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION,shp,dist_idx)


      Allocate(AI(NPES,3))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI,
     &     horzRelLoc=ESMF_CELL_CENTER,
     &       vertRelLoc=ESMF_CELL_CELL, rc=ier)

      allocate (rcounts(NPES), displs(NPES), stat=ier)
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
     &       MPI_COMM_WORLD, ier)
      Else
        Call MPI_GatherV(arr_loc(offset), scount, orig_type,
     &       arr_glob(1), rcounts, displs, new_type,
     &       root, MPI_COMM_WORLD, ier)
      End If

      End Do

      Call MPI_Type_Free(new_type, ier)

      Deallocate(rcounts,displs)
      Call MPI_Type_Free(orig_type, ier)

      End SUBROUTINE gather

      Subroutine scatter(grid, arr_glob, arr_loc, shp, dist_idx)
      Type (Esmf_Grid) :: grid
      Real(Kind=8) :: arr_loc(*)
      Real(Kind=8) :: arr_glob(*)
      Integer :: shp(:),shp_glob(size(shp))
      Integer :: dist_idx

      Integer :: new_type, orig_type
      Integer :: ier
      Integer :: p, rcount, offset
      Integer, Allocatable :: scounts(:), displs(:)
      Type (ESMF_Axisindex), Pointer :: AI(:,:)

      ! create a new mpi type for use in communication
      !-------------------------------
      new_type = CreateDist_MPI_Type(MPI_DOUBLE_PRECISION,shp,dist_idx)


      Allocate(AI(NPES,3))
      Call ESMF_GridGetAllAxisIndex(grid, globalAI=AI,
     &   horzRelLoc=ESMF_CELL_CENTER, vertRelLoc=ESMF_CELL_CELL, rc=ier)

      allocate (scounts(NPES), displs(NPES), stat=ier)
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
     &     root, MPI_COMM_WORLD, ier)

      Call MPI_Type_Free(new_type, ier)
      Call MPI_Type_Free(orig_type, ier)

      Deallocate(scounts,displs)
      Deallocate(AI)

      End SUBROUTINE scatter
#endif

      END MODULE DOMAIN_DECOMP

