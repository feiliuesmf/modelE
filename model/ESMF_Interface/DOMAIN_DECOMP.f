#define DEBUG_DECOMP
      MODULE DOMAIN_DECOMP
!@sum  DOMAIN_DECOMP encapsulates lat-lon decomposition information
!@sum  for the message passing (ESMF) implementation.
!@auth NCCS ASTG

      IMPLICIT NONE
      SAVE
      PRIVATE ! Except for

!@var DYN_GRID derived type to provide ESMF decomposition info
!@var public components are used to minimize overhead for accessing
!@var routine components
      PUBLIC :: DYN_GRID 
!@var  grid Default decomposition; globally accessible for convenience.
      PUBLIC :: grid
!@var INIT_DECOMP Initialize default decomposition
      PUBLIC :: INIT_DECOMP
!@var INIT_DECOMP Cleans up at the end of the run (closes debugging file)
      PUBLIC :: FINISH_DECOMP
!@var HALO_UPDATE Update data in halo for local domain using data from
!@var neihboring processes       
      PUBLIC :: HALO_UPDATE ! Communicate overlapping portions of subdamains
!@var CHECKSUM output a bit-reproducible checksum for an array
      PUBLIC :: CHECKSUM ! Communicate overlapping portions of subdamains

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
         LOGICAL :: HAVE_NORTH_POLE ! North pole is in local domain
         LOGICAL :: HAVE_SOUTH_POLE ! South pole is in local domain
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

      LOGICAL :: init = .false.
      INTEGER :: CHECKSUM_UNIT

      CONTAINS

      ! This routine initializes the quantities described above.
      ! The initialization should proceed prior to any grid computations.
      SUBROUTINE INIT_DECOMP(IM,JM)
      USE FILEMANAGER, Only : openunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IM, JM

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      ! Trivial topology for now
      NPES     = 1
      NP_LON   = 1
      NP_LAT   = 1
      MY_RANK  = 0
      RANK_LON = 0
      RANK_LAT = 0

      GRID%I_STRT        = 1
      GRID%I_STOP        = IM
      GRID%I_STRT_HALO   = 1
      GRID%I_STOP_HALO   = IM

      GRID%J_STRT        = 1
      GRID%J_STOP        = JM
      GRID%J_STRT_SKP    = 2
      GRID%J_STOP_SKP    = JM-1
      GRID%J_STRT_HALO   = 1
      GRID%J_STOP_HALO   = JM

      GRID%J_STRT_STGR   = 2
      GRID%J_STOP_STGR   = JM

      GRID%HAVE_NORTH_POLE = .TRUE.
      GRID%HAVE_SOUTH_POLE = .TRUE.
      GRID%HAVE_EQUATOR    = .TRUE.


#ifdef DEBUG_DECOMP
      CALL openunit('CHKSUM_DECOMP', CHECKSUM_UNIT)
#endif

      END SUBROUTINE INIT_DECOMP

      SUBROUTINE HALO_UPDATE_1D(grd, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd
      REAL*8,            INTENT(INOUT) :: 
     &                        arr(grd%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_1D

      SUBROUTINE HALO_UPDATE_2D(grd, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd
      REAL*8,            INTENT(INOUT) :: 
     &                        arr(grd%i_strt_halo:,grd%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_2D

      SUBROUTINE HALO_UPDATE_3D(grd, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd
      REAL*8,            INTENT(INOUT) :: 
     &                       arr(grd%i_strt_halo:,grid%j_strt_halo:,:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_3D

      SUBROUTINE HALO_UPDATE_COLUMN(grd, arr, from)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN)    :: grd
      REAL*8,            INTENT(INOUT) :: 
     &                       arr(:,grd%i_strt_halo:,grid%j_strt_halo:)
      INTEGER, OPTIONAL, INTENT(IN)    :: from

      END SUBROUTINE HALO_UPDATE_COLUMN

      SUBROUTINE CHECKSUM_1D(grd, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd
      REAL*8,            INTENT(IN) :: 
     &                arr(grid%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      j_0 = grid%j_strt
      j_1 = grid%j_stop

      asum   = Sum(    arr(j_0:j_1)   )
      L1norm = Sum(Abs(arr(j_0:j_1))  )
      L2norm = Sum(    arr(j_0:j_1)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_1D

      SUBROUTINE CHECKSUM_2D(grd, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd
      REAL*8,            INTENT(IN) :: 
     &                arr(grd%i_strt_halo:,grid%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      asum   = Sum(    arr(i_0:i_1,j_0:j_1)   )
      L1norm = Sum(Abs(arr(i_0:i_1,j_0:j_1))  )
      L2norm = Sum(    arr(i_0:i_1,j_0:j_1)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_2D

      SUBROUTINE CHECKSUM_3D(grd, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd
      REAL*8,            INTENT(IN) :: 
     &                arr(grd%i_strt_halo:,grid%j_strt_halo:,:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      asum   = Sum(    arr(i_0:i_1,j_0:j_1,:)   )
      L1norm = Sum(Abs(arr(i_0:i_1,j_0:j_1,:))  )
      L2norm = Sum(    arr(i_0:i_1,j_0:j_1,:)**2)
      
      Write(unit_,'(a20,1x,i5,1x,3(e22.17,1x))') 
     &        file,line, asum, L1norm, L2norm

#endif

      END SUBROUTINE CHECKSUM_3D

      SUBROUTINE CHECKSUM_COLUMN(grd, arr, line, file, unit)
      IMPLICIT NONE
      TYPE (DYN_GRID),   INTENT(IN) :: grd
      REAL*8,            INTENT(IN) :: 
     &                arr(:,grd%i_strt_halo:,grid%j_strt_halo:)
      INTEGER,           INTENT(IN) :: line
      CHARACTER(LEN=*),  INTENT(IN) :: file
      INTEGER, OPTIONAL, INTENT(IN) :: unit


      INTEGER :: unit_
      INTEGER :: i_0, i_1, j_0, j_1
      REAL*8  :: asum, L1norm, L2norm

#ifdef DEBUG_DECOMP
      unit_ = CHECKSUM_UNIT ! default
      If (Present(unit)) unit_ = unit

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

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
