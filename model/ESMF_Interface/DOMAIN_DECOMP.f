      MODULE DOMAIN_DECOMP
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@+    for the message passing (ESMF) implementation.
!@auth NCCS ASTG

      IMPLICIT NONE
      SAVE
      PRIVATE ! Except for

!@var DYN_GRID derived type to provide ESMF decomposition info
!@+   public components are used to minimize overhead for accessing
!@+   routine components
      PUBLIC :: DYN_GRID 
!@var  grid Default decomposition; globally accessible for convenience.
      PUBLIC :: grid
!@var INIT_DECOMP Initialize default decomposition
      PUBLIC :: INIT_DECOMP
!@var FINISH_DECOMP Cleans up at the end of the run (closes debugging file)
      PUBLIC :: FINISH_DECOMP
!@var HALO_UPDATE Update data in halo for local domain using data from
!@+   neighbouring processes       
      PUBLIC :: HALO_UPDATE ! Communicate overlapping portions of subdomains
!@var CHECKSUM output a bit-reproducible checksum for an array
      PUBLIC :: CHECKSUM ! Communicate overlapping portions of subdomains
!@var GET - extracts bounds information from DYN_GRID object
      PUBLIC :: GET

!ia since DYN_GRID is public ESMF_GRID_TYPE has to be public
!ia (SGI compiler complains)
      PUBLIC :: ESMF_GRID_TYPE

!@var HALO_UPDATE Generic wrapper for 2D and 3D routines
      INTERFACE HALO_UPDATE
        MODULE PROCEDURE HALO_UPDATE_1D  ! J
        MODULE PROCEDURE HALO_UPDATE_2D  ! I,J
        MODULE PROCEDURE HALO_UPDATE_3D  ! I,J,K
      END INTERFACE

      INTERFACE CHECKSUM
        MODULE PROCEDURE CHECKSUM_1D
        MODULE PROCEDURE CHECKSUM_2D
        MODULE PROCEDURE CHECKSUM_3D
      END INTERFACE

      PUBLIC :: HALO_UPDATE_COLUMN ! K, I, J
      PUBLIC :: CHECKSUM_COLUMN ! K, I, J

      ! Direction bits
      PUBLIC :: NORTH, SOUTH, EAST, WEST
      INTEGER, PARAMETER :: NORTH = 1
      INTEGER, PARAMETER :: SOUTH = 2
      INTEGER, PARAMETER :: EAST  = 4
      INTEGER, PARAMETER :: WEST  = 8

      ! Place holder for actual ESMF type
      TYPE ESMF_GRID_TYPE
         PRIVATE
         INTEGER :: id ! stub - cannot have empty derived types
      END TYPE ESMF_GRID_TYPE

      ! Local grid information
      TYPE DYN_GRID
         TYPE (ESMF_GRID_TYPE), POINTER :: ESMF_GRID
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
      END TYPE DYN_GRID

      TYPE (DYN_GRID) :: GRID

! Remaining variables are private to the module.

!@var NPES number of processes upon which work is to be distributed
      INTEGER :: NPES
!@var NP_LON number of azimuthal processes.
      INTEGER :: NP_LON
!@var NP_LAT number of meridional     processes.
      INTEGER :: NP_LAT
!@var MY_RANK index of _this_ process among 2D process topology.
      INTEGER :: MY_RANK
!@var RANK_LON index of _this_ process in azimuthal set.
      INTEGER :: RANK_LON
!@var RANK_LAT_RANK index of _this_ process in meridional set.
      INTEGER :: RANK_LAT

      INTEGER :: CHECKSUM_UNIT

      CONTAINS

      ! This routine initializes the quantities described above.
      ! The initialization should proceed prior to any grid computations.
      SUBROUTINE INIT_DECOMP(grd_dum, IM,JM)
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      TYPE (DYN_GRID), INTENT(INOUT) :: grd_dum
      INTEGER, INTENT(IN) :: IM, JM


      ! Trivial topology for now
      NPES     = 1
      NP_LON   = 1
      NP_LAT   = 1
      MY_RANK  = 0
      RANK_LON = 0
      RANK_LAT = 0

      grd_dum%I_STRT        = 1
      grd_dum%I_STOP        = IM
      grd_dum%I_STRT_HALO   = 1
      grd_dum%I_STOP_HALO   = IM

      grd_dum%J_STRT        = 1
      grd_dum%J_STOP        = JM
      grd_dum%J_STRT_SKP    = 2
      grd_dum%J_STOP_SKP    = JM-1
      grd_dum%J_STRT_HALO   = 1
      grd_dum%J_STOP_HALO   = JM

      grd_dum%J_STRT_STGR   = 2
      grd_dum%J_STOP_STGR   = JM

      grd_dum%HAVE_SOUTH_POLE = .TRUE.
      grd_dum%HAVE_NORTH_POLE = .TRUE.
      grd_dum%HAVE_EQUATOR    = .TRUE.


#ifdef DEBUG_DECOMP
      CALL openunit('CHKSUM_DECOMP', CHECKSUM_UNIT)
#endif

      END SUBROUTINE INIT_DECOMP

      SUBROUTINE GET(grd_dum, I_STRT, I_STOP, I_STRT_HALO, I_STOP_HALO,
     &                        J_STRT, J_STOP, J_STRT_HALO, J_STOP_HALO,
     &                        J_STRT_SKP, J_STOP_SKP,      
     &                        J_STRT_STGR, J_STOP_STGR,      
     &                        HAVE_SOUTH_POLE, HAVE_NORTH_POLE)
      TYPE (DYN_GRID), INTENT(IN) :: grd_dum
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
     &             HAVE_SOUTH_POLE=grd_dum%HAVE_SOUTH_POLE
      IF (PRESENT(HAVE_NORTH_POLE)) 
     &             HAVE_NORTH_POLE=grd_dum%HAVE_NORTH_POLE

      END SUBROUTINE GET

      SUBROUTINE HALO_UPDATE_1D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                        arr(grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_1D

      SUBROUTINE HALO_UPDATE_2D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                    arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_2D

      SUBROUTINE HALO_UPDATE_3D(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                 arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_3D

      SUBROUTINE HALO_UPDATE_COLUMN(grd_dum, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd_dum
      REAL*8,            INTENT(INOUT) :: 
     &                  arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_COLUMN

      SUBROUTINE CHECKSUM_1D(grd_dum, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: 
     &                arr(grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      j_0 = grd_dum%j_strt
      j_1 = grd_dum%j_stop

      asum   = Sum(    arr(j_0:j_1)   )
      L1norm = Sum(Abs(arr(j_0:j_1))  )
      L2norm = Sum(    arr(j_0:j_1)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_1D

      SUBROUTINE CHECKSUM_2D(grd_dum, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: 
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grd_dum%i_strt
      i_1 = grd_dum%i_stop
      j_0 = grd_dum%j_strt
      j_1 = grd_dum%j_stop

      asum   = Sum(    arr(i_0:i_1,j_0:j_1)   )
      L1norm = Sum(Abs(arr(i_0:i_1,j_0:j_1))  )
      L2norm = Sum(    arr(i_0:i_1,j_0:j_1)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_2D

      SUBROUTINE CHECKSUM_3D(grd_dum, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: 
     &                arr(grd_dum%i_strt_halo:,grd_dum%j_strt_halo:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grd_dum%i_strt
      i_1 = grd_dum%i_stop
      j_0 = grd_dum%j_strt
      j_1 = grd_dum%j_stop

      asum   = Sum(    arr(i_0:i_1,j_0:j_1,:)   )
      L1norm = Sum(Abs(arr(i_0:i_1,j_0:j_1,:))  )
      L2norm = Sum(    arr(i_0:i_1,j_0:j_1,:)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_3D

      SUBROUTINE CHECKSUM_COLUMN(grd_dum, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: 
     &                arr(:,grd_dum%i_strt_halo:,grd_dum%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grd_dum%i_strt
      i_1 = grd_dum%i_stop
      j_0 = grd_dum%j_strt
      j_1 = grd_dum%j_stop

      asum   = Sum(    arr(:,i_0:i_1,j_0:j_1)   )
      L1norm = Sum(Abs(arr(:,i_0:i_1,j_0:j_1))  )
      L2norm = Sum(    arr(:,i_0:i_1,j_0:j_1)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_COLUMN

      SUBROUTINE FINISH_DECOMP()
      USE FILEMANAGER, ONLY : closeunit
      IMPLICIT NONE

#ifdef DEBUG_DECOMP
      CALL closeunit(CHECKSUM_UNIT)
#endif

      END SUBROUTINE FINISH_DECOMP

      END MODULE DOMAIN_DECOMP
