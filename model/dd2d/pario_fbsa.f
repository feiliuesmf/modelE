!@auth M. Kelley
!@ver 1.0
!@sum  pario_fbsa contains multitile-enabled versions of routines
!@+    that fill real*8 domain-decomposed arrays with the contents of
!@+    files written using Fortran Binary Sequential Access.
!@+    The first two indices of arrays are assumed to correspond to
!@+    the two horizontal dimensions.   Fortran records have the
!@+    following structure in each routine, where im/jm are the
!@+    number of points in i/j on each tile and ntiles is the number
!@+    of tiles (1 for a lat-lon grid, 6 for a cubed sphere grid).
!@+
!@+    READT_PARALLEL:  title*80, real*4 arr(im,jm,ntiles)
!@+                or   title*80, real*4 arr(im,jm,:,ntiles)
!@+    DREAD_PARALLEL:  real*4 arr(im,jm,ntiles)
!@+                or   real*4 arr(im,jm,:,ntiles)
!@+    MREAD_PARALLEL:  integer m, real*4 arr(im,jm,ntiles)
!@+                or   integer m, real*4 arr(im,jm,:,ntiles)
!@+
!@+    BACKSPACE, REWIND, and SKIP are provided for
!@+    plug compatibility with the DOMAIN_DECOMP module.
!@+    
!@+    Although intended for cubed sphere migration purposes, the
!@+    routines in this module will also work for a single-tile
!@+    latitude-longitude grid with 1D domain decomposition provided
!@+    that the dd2d component of the dist_grid input argument is
!@+    properly initialized.
!@+
      module pario_fbsa
      use filemanager, only : openunit,closeunit
      use domain_decomp, only : dist_grid,am_i_root,esmf_bcast
      use dd2d_utils, only : unpack_data
      private

      public :: openunit,closeunit
      PUBLIC :: BACKSPACE_PARALLEL
      PUBLIC :: REWIND_PARALLEL
      PUBLIC :: SKIP_PARALLEL
      PUBLIC :: DREAD_PARALLEL
      PUBLIC :: MREAD_PARALLEL
      PUBLIC :: READT_PARALLEL
      PUBLIC :: READ_PARALLEL

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
      end interface

      contains

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
      CHARACTER(len=*), INTENT(IN)  :: NAME      !@var  NAME  name of record being read
      REAL*8,INTENT(OUT) :: AVAR(:,:)
      REAL*4,allocatable :: AIN(:,:,:)  !@var  AIN  real*4 array for reading
      REAL*8,allocatable :: AOUT(:,:,:) !@var  AOUT real*8 array for scatter
      INTEGER :: IERR
    ! now local

      If (AM_I_ROOT()) then
        allocate(
     &       AIN (grd_dum%IM_WORLD,grd_dum%JM_WORLD,
     &       grd_dum%dd2d%ntiles),
     &       AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,
     &       grd_dum%dd2d%ntiles)   )
         READ (IUNIT,IOSTAT=IERR) AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

      call unpack_data(grd_dum%dd2d,aout,avar)

      if (AM_I_ROOT()) then
         deallocate (ain,aout)
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('DREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if
      return
      END SUBROUTINE DREAD_PARALLEL_2D

      SUBROUTINE DREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,AVAR)
!@sum DREAD_PARALLEL  Parallel version of UTILDBL.f:DREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER(len=*), INTENT(IN)  :: NAME        !@var  NAME  name of record being read
      REAL*8, INTENT(OUT) :: AVAR(:,:,:) !@var  AOUT real*8 array
      real*4, allocatable :: ain(:,:,:,:)  !@var  AIN  real*4 array for reading
      real*8, allocatable :: aout(:,:,:,:) !@var  AOUT real*8 array for scatter
      INTEGER :: IERR

      If (AM_I_ROOT()) then
         allocate(
     &       AIN (grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3),
     &       grd_dum%dd2d%ntiles),
     &       AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3),
     &       grd_dum%dd2d%ntiles)   )
         READ (IUNIT,IOSTAT=IERR) AIN
C**** convert from real*4 to real*8
         AOUT=AIN
      EndIf

      call unpack_data(grd_dum%dd2d,aout,avar)
      
      if (AM_I_ROOT()) then
         deallocate (ain,aout)
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('DREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if
      return
      END SUBROUTINE DREAD_PARALLEL_3D

      SUBROUTINE MREAD_PARALLEL_2D (grd_dum,IUNIT,NAME,M,NSKIP,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT     !@var  IUNIT file unit number
      CHARACTER(len=*), INTENT(IN)  :: NAME      !@var  NAME  name of record being read
      INTEGER,      INTENT(OUT) :: M         !@var  M      initial integer
      INTEGER,      INTENT(IN)  :: NSKIP     !@var  NSKIP no. of R*4's to skip
      REAL*8, INTENT(OUT)  :: AVAR(:,:) !@var  AOUT real*8 array
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4,allocatable :: AIN(:,:,:)  !@var  AIN  real*4 array for reading
      REAL*8,allocatable :: AOUT(:,:,:) !@var  AOUT real*8 array for scatter
      INTEGER :: N,IERR

      If (AM_I_ROOT()) then
        allocate(
     &       AIN (grd_dum%IM_WORLD,grd_dum%JM_WORLD,
     &       grd_dum%dd2d%ntiles),
     &       AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,
     &       grd_dum%dd2d%ntiles)   )
         READ (IUNIT,IOSTAT=IERR) M,(X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

      call unpack_data(grd_dum%dd2d,aout,avar)
      CALL ESMF_BCAST(grd_dum, M   )

      if (AM_I_ROOT()) then
         deallocate (ain,aout)
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('MREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if
      return
      END SUBROUTINE MREAD_PARALLEL_2D

      SUBROUTINE MREAD_PARALLEL_3D (grd_dum,IUNIT,NAME,M,NSKIP,AVAR)
!@sum MREAD_PARALLEL  Parallel version of UTILDBL.f:MREAD for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT       !@var  IUNIT file unit number
      CHARACTER(len=*), INTENT(IN)  :: NAME        !@var  NAME  name of record being read
      INTEGER,      INTENT(OUT) :: M           !@var  M      initial integer
      INTEGER,      INTENT(IN)  :: NSKIP       !@var  NSKIP no. of R*4's to skip
      REAL*8,      INTENT(OUT)  :: AVAR(:,:,:) !@var  AOUT real*8 array
      REAL*4 :: X                         !@var  X dummy variable
      real*4, allocatable :: ain(:,:,:,:) !@var  AIN  real*4 array for reading
      real*8, allocatable :: aout(:,:,:,:)!@var  AOUT real*8 array for scatter
      INTEGER :: N,IERR

      If (AM_I_ROOT()) then
         allocate(
     &       AIN (grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3),
     &       grd_dum%dd2d%ntiles),
     &       AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3),
     &       grd_dum%dd2d%ntiles)   )
         READ (IUNIT,IOSTAT=IERR) M,(X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

      call unpack_data(grd_dum%dd2d,aout,avar)
      CALL ESMF_BCAST(grd_dum, M   )

      if (AM_I_ROOT()) then
         deallocate(ain,aout)
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': IOSTAT=',IERR
            call stop_model('MREAD_PARALLEL: READ ERROR',255)
         EndIf
      end if
      return
      END SUBROUTINE MREAD_PARALLEL_3D

      SUBROUTINE READT_PARALLEL_2D (grd_dum,IUNIT,NAME,NSKIP,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT      !@var  IUNIT file unit number
      CHARACTER(len=*), INTENT(IN)  :: NAME       !@var  NAME  name of record being read
      INTEGER,      INTENT(IN)  :: NSKIP      !@var  NSKIP no. of R*4's to skip
      REAL*8,       INTENT(OUT) :: AVAR(:,:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS       !@var  IPOS  no. of recs. to advance
      REAL*4 :: X                         !@var  X dummy variable
      REAL*4,allocatable :: AIN(:,:,:)  !@var  AIN  real*4 array for reading
      REAL*8,allocatable :: AOUT(:,:,:) !@var  AOUT real*8 array for scatter
      INTEGER :: N
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: IERR

      If (AM_I_ROOT()) then
        allocate(
     &       AIN (grd_dum%IM_WORLD,grd_dum%JM_WORLD,
     &       grd_dum%dd2d%ntiles),
     &       AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,
     &       grd_dum%dd2d%ntiles)   )
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=IERR)
         END DO
         READ (IUNIT, IOSTAT=IERR) TITLE, (X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

      call unpack_data(grd_dum%dd2d,aout,avar)

      if (AM_I_ROOT()) then
         deallocate(ain,aout)
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ': ',
     &           TRIM(TITLE),' IOSTAT=',IERR
            call stop_model('READT_PARALLEL: READ ERROR',255)
         EndIf
      end if
      return
      END SUBROUTINE READT_PARALLEL_2D

      SUBROUTINE READT_PARALLEL_3D (grd_dum,IUNIT,NAME,NSKIP,AVAR,IPOS)
!@sum READT_PARALLEL  Parallel version of UTILDBL.f:READT for (im,jm) arrays
!@auth NCCS-ESMF Development Team
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      INTEGER,      INTENT(IN)  :: IUNIT        !@var  IUNIT file unit number
      CHARACTER(len=*), INTENT(IN)  :: NAME         !@var  NAME  name of record being read
      INTEGER,      INTENT(IN)  :: NSKIP        !@var  NSKIP no. of R*4's to skip
      REAL*8,       INTENT(OUT) :: AVAR(:,:,:)  !@var  AOUT real*8 array
      INTEGER,      INTENT(IN)  :: IPOS         !@var  IPOS  no. of recs. to advance
      REAL*4 :: X                         !@var  X dummy variable
      real*4, allocatable :: ain(:,:,:,:) !@var  AIN  real*4 array for reading
      real*8, allocatable :: aout(:,:,:,:)!@var  AOUT real*8 array for scatter
      INTEGER :: N                        !@var  N loop variable
      CHARACTER*80 :: TITLE               !@var  TITLE title of file record
      INTEGER :: IERR

      If (AM_I_ROOT()) then
         allocate(
     &       AIN (grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3),
     &       grd_dum%dd2d%ntiles),
     &       AOUT(grd_dum%IM_WORLD,grd_dum%JM_WORLD,size(AVAR,3),
     &       grd_dum%dd2d%ntiles)   )
         DO N=1,IPOS-1
            READ (IUNIT,IOSTAT=IERR)
         END DO
         READ (IUNIT, IOSTAT=IERR) TITLE, (X,N=1,NSKIP), AIN
C****  convert from real*4 to real*8
         AOUT=AIN
      EndIf

      call unpack_data(grd_dum%dd2d,aout,avar)

      if (am_i_root()) then
         deallocate(ain,aout)
         If (IERR==0) Then
            WRITE(6,*) "Read from file ",TRIM(NAME),": ",TRIM(TITLE)
            RETURN
         Else
            WRITE(6,*) 'READ ERROR ON FILE ',NAME, ':
     &           ',TRIM(TITLE),' IOSTAT=',IERR
            call stop_model('READT_PARALLEL: READ ERROR',255)
         EndIf
      end if
      return
      END SUBROUTINE READT_PARALLEL_3D

      subroutine read_parallel_integer_0 (data_int, unit)
      integer, intent(out) :: data_int
      integer, intent(in) :: unit
      type(dist_grid) :: grd_dum ! only used to satisfy esmf_bcast interface
      if(am_i_root()) read(unit) data_int
      CALL ESMF_BCAST(grd_dum, data_int)
      end subroutine read_parallel_integer_0

      end module pario_fbsa
