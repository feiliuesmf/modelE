C****
C**** ICEDYN_DRV.f    Sea ICE DYNamics    2006/12/21
C****
#include "rundeck_opts.h"

      MODULE ICEDYN_COM
!@sum  ICEDYN_COM holds global variables for dynamic sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE ICEDYN, only : imic,jmic
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE

C**** variables used in advection (on ICE grid)
!@var RSIX,RSIY first order moments for seaice concentration
!@var USI,VSI east-west, and north-south sea ice velocities (m/s)
!     REAL*8, DIMENSION(IMIC,JMIC) :: RSIX,RSIY,USI,VSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RSIX,RSIY,USI,VSI

!@var USIDT,VSIDT sea ice fluxes, saved for advection (m)
!     REAL*8, DIMENSION(IMIC,JMIC) :: USIDT,VSIDT
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: USIDT,VSIDT

C**** Needed for ADVSI (on ATM grid)
!@var RSISAVE saved value of sea ice concentration before DYNSI
!     REAL*8, DIMENSION(IM,JM) :: RSISAVE
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RSISAVE

C**** Ice advection diagnostics
      INTEGER, PARAMETER :: KICIJ=12
!@var IJ_xxx Names for ICIJ diagnostics
      INTEGER IJ_USI,IJ_VSI,IJ_DMUI,IJ_DMVI,IJ_PICE,IJ_MUSI,IJ_MVSI
     *     ,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI,IJ_RSI
!@var ICIJ lat-lon ice dynamic diagnostics (on atm grid)
!     REAL*8, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: ICIJ
!@var lname_icij Long names for ICIJ diagnostics
      CHARACTER*50, DIMENSION(KICIJ) :: LNAME_ICIJ
!@var sname_icij Short names for ICIJ diagnostics
      CHARACTER*30, DIMENSION(KICIJ) :: SNAME_ICIJ
!@var units_icij Units for ICIJ diagnostics
      CHARACTER*50, DIMENSION(KICIJ) :: UNITS_ICIJ
!@var ia_icij IDACC numbers for ICIJ diagnostics
      INTEGER, DIMENSION(KICIJ) :: IA_ICIJ
!@var scale_icij scales for ICIJ diagnostics
      REAL*8, DIMENSION(KICIJ) :: SCALE_ICIJ
!@var ijgrid_icij Grid descriptor for ICIJ diagnostics
       INTEGER, DIMENSION(KICIJ) :: IJGRID_ICIJ
#ifdef TRACERS_WATER
!@var KTICIJ number of lat/lon ice dynamic tracer diagnostics
      INTEGER, PARAMETER :: KTICIJ=2
!@var TICIJ  lat/lon ice dynamic tracer diagnostics
!     REAL*8, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)  :: TICIJ
!@var ticij_xxx indices for TICIJ diags
      INTEGER :: tICIJ_tusi,tICIJ_tvsi
!@var lname_ticij Long names for TICIJ diagnostics
      CHARACTER*50, DIMENSION(KTICIJ) :: LNAME_TICIJ
!@var sname_ticij Short names for TICIJ diagnostics
      CHARACTER*30, DIMENSION(KTICIJ) :: SNAME_TICIJ
!@var units_ticij Units for TICIJ diagnostics
      CHARACTER*50, DIMENSION(KTICIJ) :: UNITS_TICIJ
!@var ia_ticij IDACC numbers for TICIJ diagnostics
      INTEGER, DIMENSION(KTICIJ) :: IA_TICIJ
!@var scale_ticij scales for TICIJ diagnostics
      REAL*8, DIMENSION(KTICIJ) :: SCALE_TICIJ
!@var ijgrid_ticij Grid descriptor for TICIJ diagnostics
      INTEGER, DIMENSION(KTICIJ) :: IJGRID_TICIJ
#endif

      END MODULE ICEDYN_COM

      SUBROUTINE ALLOC_ICEDYN_COM(grid)
!@sum ALLOC_ICEDYN_COM allocates arrays defined in the ICEDYN_COM module
!@auth Rosalinda de Fainchtein

      USE DOMAIN_DECOMP, only : GET
      USE DOMAIN_DECOMP, only : DIST_GRID
      USE MODEL_COM, only : im
      USE ICEDYN, only : imic
      USE ICEDYN, only : grid_MIC
      USE ICEDYN_COM, only : KICIJ
      USE ICEDYN_COM, only : RSIX,RSIY,USI,VSI,USIDT,VSIDT,
     &                       RSISAVE,ICIJ
#ifdef TRACERS_WATER
      USE ICEDYN_COM, only : TICIJ,KTICIJ,NTM
#endif
      IMPLICIT NONE

      LOGICAL, SAVE :: init=.false.
      INTEGER :: J_1H    , J_0H
      INTEGER :: J_1H_MIC, J_0H_MIC
      INTEGER :: IER
      TYPE(DIST_GRID) :: grid

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

!*** For now set grid_MIC to be the same as grid
!    This is consistent with the current status of the code for
!    parallelization along latitude (j)
      grid_MIC =grid

C**** Get dimensioning parameters for arrays defined in the grid  and 
C**** grid_MIC stencils.

      CALL GET(grid    , J_STRT_HALO=J_0H    , J_STOP_HALO=J_1H    )
      CALL GET(grid_MIC, J_STRT_HALO=J_0H_MIC, J_STOP_HALO=J_1H_MIC)

      ALLOCATE( RSIX(IMIC, J_0H_MIC:J_1H_MIC),
     &          RSIY(IMIC, J_0H_MIC:J_1H_MIC),
     &           USI(IMIC, J_0H_MIC:J_1H_MIC),
     &           VSI(IMIC, J_0H_MIC:J_1H_MIC),
     &     STAT = IER)

      ALLOCATE( USIDT(IMIC, J_0H_MIC:J_1H_MIC),
     &          VSIDT(IMIC, J_0H_MIC:J_1H_MIC),
     &     STAT = IER)

      ALLOCATE( RSISAVE(IM, J_0H:J_1H),
     &     STAT = IER)

      ALLOCATE(  ICIJ(IMIC, J_0H_MIC:J_1H_MIC, KICIJ),
     &     STAT = IER)

#ifdef TRACERS_WATER
      ALLOCATE( TICIJ(IMIC, J_0H_MIC:J_1H_MIC, KTICIJ, NTM),
     &     STAT = IER)
#endif

      return
      END SUBROUTINE ALLOC_ICEDYN_COM


      SUBROUTINE io_icedyn(kunit,iaction,ioerr)
!@sum  io_icedyn reads and writes dynamic ice arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irsficno,irsficnt
     *     ,irerun,lhead
      USE DOMAIN_DECOMP, only: AM_I_ROOT, PACK_DATA, UNPACK_DATA
      USE ICEDYN, only : grid_MIC
      USE ICEDYN_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICEDYN01"
      REAL*8, DIMENSION(IMIC, JMIC) :: RSIX_glob,RSIY_glob,USI_glob
     &                                ,VSI_glob

      write(MODULE_HEADER(lhead+1:80),'(a7,i3,a1,i3,a)')
     *     'R8 dim(',imic,',',jmic,'):RSIX,RSIY,USI,VSI'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_DATA(grid_MIC, RSIX, RSIX_GLOB)
        CALL PACK_DATA(grid_MIC, RSIY, RSIY_GLOB)
        CALL PACK_DATA(grid_MIC,  USI,  USI_GLOB)
        CALL PACK_DATA(grid_MIC,  VSI,  VSI_GLOB)
        IF (AM_I_ROOT())
     &   WRITE (kunit,err=10) MODULE_HEADER,RSIX_glob,RSIY_glob
     &                                     , USI_glob, VSI_glob
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFICNO)           ! initial conditions (no ocean)
        CASE (ioread,irerun,irsfic,irsficnt)    ! restarts
          if (AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,       RSIX_glob,RSIY_glob
     &                                         , USI_glob, VSI_glob
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER,
     &             MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_DATA(grid_MIC, RSIX_GLOB, RSIX)
          CALL UNPACK_DATA(grid_MIC, RSIY_GLOB, RSIY)
          CALL UNPACK_DATA(grid_MIC,  USI_GLOB,  USI)
          CALL UNPACK_DATA(grid_MIC,  VSI_GLOB,  VSI)
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icedyn

      SUBROUTINE io_icdiag(kunit,it,iaction,ioerr)
!@sum  io_icdiag reads and writes ice dynamic diagnostic arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irsfic,irsficnt,irerun,ioread_single,lhead
      USE DOMAIN_DECOMP, only : GET, AM_I_ROOT
      USE DOMAIN_DECOMP, only : PACK_DATA, UNPACK_DATA
      USE DOMAIN_DECOMP, only : ESMF_BCAST
      USE ICEDYN, only : grid_MIC
      USE ICEDYN_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICDIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var ICIJ4 dummy arrays for reading diag. files
      REAL*8, DIMENSION(:,:,:), allocatable  :: ICIJ4
      REAL*4, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ4_GLOB
      REAL*8, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ4_GLOB8
      REAL*8, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ_GLOB
      INTEGER :: J_0H_MIC, J_1H_MIC

#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), allocatable  :: TICIJ4
      REAL*4, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ4_GLOB
      REAL*8, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ4_GLOB8
      REAL*8, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ_GLOB
!@var TR_HEADER Character string label for individual tracer records
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRICDIAG01"

      write(TR_MODULE_HEADER(lhead+1:80),'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')
     *     'R8 Ticij(',imic,',',jmic,',',kticij,',',ntm,'),it'
#endif
      write(MODULE_HEADER(lhead+1:80),'(a8,i3,a1,i3,a1,i2,a4)')
     *     'R8 ICij(',imic,',',jmic,',',kicij,'),it'

      CALL GET(grid_MIC, J_STRT_HALO=J_0H_MIC, J_STOP_HALO=J_1H_MIC)

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON)  ! output to standard restart file
        CALL PACK_DATA(grid_mic, icij, icij_glob)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) MODULE_HEADER,ICIJ_glob,it
#ifdef TRACERS_WATER
        CALL PACK_DATA(grid_mic, TICIJ, TICIJ_GLOB)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) TR_MODULE_HEADER,TICIJ_GLOB,it
#endif
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        CALL PACK_DATA(grid_mic, icij, icij_glob)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) MODULE_HEADER,REAL(ICIJ_GLOB,KIND=4),it
#ifdef TRACERS_WATER
        TR_MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        CALL PACK_DATA(grid_mic, TICIJ, TICIJ_GLOB)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) TR_MODULE_HEADER
     &          ,REAL(TICIJ_GLOB,KIND=4),it
#endif
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,ICIJ4_GLOB,it
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
          endif
C**** accumulate diagnostics
          allocate( ICIJ4(IMIC, J_0H_MIC:J_1H_MIC, KICIJ) )
          ICIJ4 = 0.d0 ! should do "halo_ipdate" instead?
          ICIJ4_GLOB8 = ICIJ4_GLOB ! convert to real*8
          CALL UNPACK_DATA(grid_MIC, ICIJ4_GLOB8, ICIJ4)
          call ESMF_BCAST(grid_MIC, it)
          ICIJ(:,J_0H_MIC:J_1H_MIC,:)=ICIJ(:,J_0H_MIC:J_1H_MIC,:)
     &                            +ICIJ4(:,J_0H_MIC:J_1H_MIC,:)
          deallocate( ICIJ4 )
#ifdef TRACERS_WATER
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TR_HEADER,TICIJ4_GLOB,it
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TR_HEADER
     *             ,TR_MODULE_HEADER
              GO TO 10
            END IF
          endif
C**** accumulate diagnostics
          allocate( TICIJ4(IMIC, J_0H_MIC:J_1H_MIC, KTICIJ, NTM) )
          TICIJ4 = 0.d0 ! should do "halo_ipdate" instead?
          TICIJ4_GLOB8 = TICIJ4_GLOB ! convert to real*8
          CALL UNPACK_DATA(grid_MIC, TICIJ4_GLOB8, TICIJ4)
          call ESMF_BCAST(grid_MIC, it)
          TICIJ(:,J_0H_MIC:J_1H_MIC,:,:)=TICIJ(:,J_0H_MIC:J_1H_MIC,:,:)
     &                                 +TICIJ4(:,J_0H_MIC:J_1H_MIC,:,:)
          deallocate( TICIJ4 )
#endif
        CASE (ioread,irerun)    ! restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,ICIJ_GLOB,it
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_DATA(grid_MIC, ICIJ_GLOB, ICIJ)
          call ESMF_BCAST(grid_MIC, it)
#ifdef TRACERS_WATER
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TR_HEADER,TICIJ_GLOB,it
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TR_HEADER
     *             ,TR_MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_DATA(grid_MIC, TICIJ_GLOB, TICIJ)
          call ESMF_BCAST(grid_MIC, it)
#endif
        CASE (IRSFIC)  ! initial conditions
          READ (kunit)
#ifdef TRACERS_WATER
          READ (kunit)
#endif
        CASE (IRSFICNT)  ! initial conditions (with no tracers)
          READ (kunit)
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icdiag

      SUBROUTINE reset_icdiag
!@sum reset_icdiag resets ice dynamic diagnostic arrays
!@auth Gavin Schmidt
      USE ICEDYN_COM
      IMPLICIT NONE

      ICIJ=0.
#ifdef TRACERS_WATER
      TICIJ=0.
#endif
      RETURN
      END SUBROUTINE reset_icdiag

      SUBROUTINE DYNSI
!@sum  DYNSI calculate ice velocites
!@+    Note that the ice velocities are calculated on the ice grid
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      USE CONSTANT, only : rhoi,grav,omega,rhows
      USE MODEL_COM, only : im,jm,p,ptop,dts=>dtsrc,focean
      USE DOMAIN_DECOMP, only : grid, DIST_GRID, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE, NORTH, SOUTH
      USE GEOM, only : dxyn,dxys,dxyv,dxyp,bydxyp,dxp,dyv,imaxj
      USE ICEDYN, only : imic,jmic,nx1,ny1,press,heffm,uvm,dwatn,cor
     *     ,sinwat,coswat,bydts,sinen,uice,vice,heff,area,gairx,gairy
     *     ,gwatx,gwaty,pgfub,pgfvb,amass,uicec,vicec,uib,vib,dmu,dmv
      USE ICEDYN_COM, only : usi,vsi,usidt,vsidt,rsisave,icij,ij_usi
     *     ,ij_vsi,ij_dmui,ij_dmvi,ij_pice,ij_rsi
      USE FLUXES, only : dmua,dmva,dmui,dmvi,UI2rho,ogeoza,uosurf,vosurf
     *     ,apress,uisurf,visurf
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      IMPLICIT NONE
      SAVE
C**** intermediate calculation for pressure gradient terms
      REAL*8, DIMENSION(IM, grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                            PGFU,PGFV
C****
      REAL*8, PARAMETER :: BYRHOI=1D0/RHOI
      REAL*8 :: hemi
      INTEGER I,J,n,k,ip1,im1,l
      REAL*8 USINP,DMUINP,duA,dvA

C**** Declare a new grid data type grid_NXY (to handle do j=1,ny1 type
C     Loops. Set grid_NXY=grid, consistent with curr. state of the code.

      TYPE(DIST_GRID) :: grid_NXY
      INTEGER :: J_1NXY, J_0NXY
      INTEGER :: J_1NXYS
      INTEGER :: J_1   , J_0
      INTEGER :: J_1S  , J_0S
      INTEGER :: J_1STG,J_0STG
      grid_NXY=grid

C**** Get loop indices  corresponding to grid and grid_NXY structures
      CALL GET(grid_NXY, J_STRT=J_0NXY   , J_STOP=J_1NXY
     &                                   , J_STOP_SKP=J_1NXYS)
      call GET(grid    , J_STRT=J_0      , J_STOP=J_1     )
      call GET(grid    , J_STRT_SKP=J_0S , J_STOP_SKP=J_1S)
      call GET(grid    , J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)

C**** Start main loop
C**** Replicate polar boxes
      if (grid%HAVE_NORTH_POLE) then
        RSI(2:IM,JM)=RSI(1,JM)
        MSI(2:IM,JM)=MSI(1,JM)
        DMUA(2:IM,JM,2) = DMUA(1,JM,2)
        DMVA(2:IM,JM,2) = DMVA(1,JM,2)
      end if

C**** save current value of sea ice concentration for ADVSI
C**** RSISAVE is on atmospheric grid
      RSISAVE(:,:)=RSI(:,:)

C**** Pressure anomaly at surface APRESS: calculated by sea ice routines
C**** APRESS is on atmospheric grid. We are now no longer using this as
C**** a forcing in the sea ice dynamics because it is partially
C**** included already in the internal ice pressure gradient. The
C**** atmospheric pressure gradient term does not in general produce a
C**** horizontal force in a solid (such as ice).

C**** calculate sea surface tilt on atmospheric C grid
C**** (using OGEOZA on atmospheric grid plus displacement of free
C**** surface due to presence of ice). This is ignored in favour of
C**** geostrophy if osurf_tilt=0.
C**** PGF is an accelaration
      if (grid%HAVE_NORTH_POLE) PGFU(1:IM,JM)=0
      if (grid%HAVE_SOUTH_POLE) PGFU(1:IM, 1)=0  !RKF
      DO J=J_0S,J_1S
        I=IM
        DO IP1=1,IM
          IF(FOCEAN(I,J).gt.0 .and. FOCEAN(IP1,J).gt.0. .and.
     *         RSI(I,J)+RSI(IP1,J).gt.0.) THEN
c            PGFU(I,J)=-((APRESS(IP1,J)-APRESS(I,J))*BYRHOI
c     *           +OGEOZA(IP1,J)-OGEOZA(I,J))/DXP(J)
            PGFU(I,J)=-(OGEOZA(IP1,J)-OGEOZA(I,J)+
     *           (RSI(IP1,J)*(MSI(IP1,J)+SNOWI(IP1,J)+ACE1I)
     *           -RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I))*GRAV/RHOWS )
     *           /DXP(J)
          ELSE
            PGFU(I,J)=0.
          END IF
          I=IP1
        END DO
      END DO
C**** Fill halos for arrays FOCEAN, RSI,OGEOZA,DYV,MSI,SNOWI
C**** Commented Halo fill for array APRESS supports commented statement
      CALL HALO_UPDATE(grid, FOCEAN, from=NORTH )
      CALL HALO_UPDATE(grid, RSI   , from=NORTH )
      CALL HALO_UPDATE(grid, OGEOZA, from=NORTH )
      CALL HALO_UPDATE(grid, MSI   , from=NORTH )
      CALL HALO_UPDATE(grid, SNOWI , from=NORTH )
      DO J=J_0,J_1S
        DO I=1,IM
          IF(FOCEAN(I,J+1).gt.0 .and. FOCEAN(I,J).gt.0. .and.
     *         RSI(I,J)+RSI(I,J+1).gt.0.) THEN
c            PGFV(I,J)=-((APRESS(I,J+1)-APRESS(I,J))*BYRHOI
c     *           +OGEOZA(I,J+1)-OGEOZA(I,J))/DYV(J+1)
            PGFV(I,J)=-(OGEOZA(I,J+1)-OGEOZA(I,J)+
     *           (RSI(I,J+1)*(MSI(I,J+1)+SNOWI(I,J+1)+ACE1I) -
     *           RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I))*GRAV/RHOWS  )
     *           /DYV(J+1)
          ELSE
            PGFV(I,J)=0.
          END IF
        END DO
      END DO

C**** DMUA is defined over the whole box (not just over ptype)
C**** Convert to stress over ice fraction only (on atmospheric grid)
      DO J=J_0,J_1
        DO I=1,IM
          IF (FOCEAN(I,J)*RSI(I,J).gt.0) THEN
            DMUA(I,J,2) = DMUA(I,J,2)/(FOCEAN(I,J)*RSI(I,J))
            DMVA(I,J,2) = DMVA(I,J,2)/(FOCEAN(I,J)*RSI(I,J))
          ELSE
            DMUA(I,J,2) = 0.
            DMVA(I,J,2) = 0.
          END IF
        END DO
      END DO

C**** Set up ice grid variables
C**** HEFF,AREA on primary (tracer) grid for ice

      do j=J_0NXY,J_1NXY
        do i=2,NX1-1
          HEFF(I,J)=RSI(I-1,J)*(ACE1I+MSI(I-1,J))*BYRHOI
          AREA(I,J)=RSI(I-1,J)
c if heffm(i,j).eq.1 .and. focean(i-1,j).eq.0 interpolate.....
c          if (heffm(i,j).eq.1 .and. focean(i-1,j).eq.0) then
c     rsifix=(rsi(i-2,j)*focean(i-2,j)+rsi(i,j)*focean(i,j)+rsi(i-1,j+1)
c    *  *focean(i-1,j+1)+rsi(i-1,j-1)*focean(i-1,j-1))/(focean(i-2,j)
c    *  +focean(i,j)+focean(i-1,j+1)+focean(i-1,j-1))
c     msifix=(msi(i-2,j)*focean(i-2,j)+msi(i,j)*focean(i,j)+msi(i-1,j+1)
c    *  *focean(i-1,j+1)+msi(i-1,j-1)*focean(i-1,j-1))/(focean(i-2,j)
c    *  +focean(i,j)+focean(i-1,j+1)+focean(i-1,j-1))
c      HEFF(I,J)=rsifix*(ACE1I+msifix)*BYRHOI
c      AREA(I,J)=rsifix
c          end if
          HEFF(I,J)=HEFF(I,J)*HEFFM(I,J)
        enddo
      enddo
C**** fill in overlap regions
      DO J=J_0NXY,J_1NXY
        HEFF(1,J)=HEFF(NX1-1,J)
        AREA(1,J)=AREA(NX1-1,J)
        HEFF(NX1,J)=HEFF(2,J)
        AREA(NX1,J)=AREA(2,J)
      END DO

C****
C**** Set up mass per unit area and coriolis term (on ice grid B)
C****
C**** Update halo for HEFF
      CALL HALO_UPDATE(grid, HEFF, from=NORTH    )

      DO J=J_0NXY,J_1NXYS
      DO I=1,NX1-1
        AMASS(I,J)=RHOI*0.25*(HEFF(I,J)
     *       +HEFF(I+1,J)+HEFF(I,J+1)+HEFF(I+1,J+1))
        COR(I,J)=AMASS(I,J)*2.0*OMEGA*SINEN(I,J)
      END DO
      END DO
c**** set north pole
      if (grid%HAVE_NORTH_POLE) then
        do i=1,nx1
          AMASS(i,jm)= RHOI*HEFF(1,JM)
          COR  (i,jm)= AMASS(i,jm)*2.0*OMEGA*SINEN(1,JM)
        end do
      end if                    !end NORTH_POLE block if

c**** interpolate air, current and ice velocity from C grid to B grid
C**** This should be more generally from ocean grid to ice grid
C**** NOTE: UOSURF, VOSURF are expected to be on the C-grid

C**** Update halo for USI,UOSURF,PGFU
      CALL HALO_UPDATE(grid, USI   , from=NORTH    )
      CALL HALO_UPDATE(grid, UOSURF, from=NORTH    )
      CALL HALO_UPDATE(grid, PGFU  , from=NORTH    )

      do j=j_0,j_1s
        im1=im
        do i=1,im
          UIB  (i,j)=0.5*(USI (im1,j)  +USI (im1,j+1))   ! iceC--> iceB
          GWATX(i,j)=0.5*(UOSURF(im1,j)+UOSURF(im1,j+1)) ! ocnC--> iceB
          PGFUB(i,j)=0.5*(PGFU(im1,j)  +PGFU(im1,j+1))   ! atmC--> iceB
          VIB  (i,j)=0.5*(VSI (im1,j)  +VSI (i,j))
          GWATY(i,j)=0.5*(VOSURF(im1,j)+VOSURF(i,j))
          PGFVB(i,j)=0.5*(PGFV(im1,j)  +PGFV(i,j))
          im1=i
        enddo
      enddo
c**** set north pole
      if (grid%HAVE_NORTH_POLE) then
        do i=1,im
          UIB  (i,jm)=USI(1,jm)
          GWATX(i,jm)=0.
          PGFUB(i,jm)=0.
          VIB  (i,jm)=0.
          GWATY(i,jm)=0.
          PGFVB(i,jm)=0.
        enddo
      end if                    !end NORTH_POLE block if
      DO J=J_0NXY,J_1NXY
        UIB(nx1-1,J)=UIB(1,J)
        VIB(nx1-1,J)=VIB(1,J)
        GWATX(nx1-1,J)=GWATX(1,J)
        GWATY(nx1-1,J)=GWATY(1,J)
        PGFUB(nx1-1,J)=PGFUB(1,J)
        PGFVB(nx1-1,J)=PGFVB(1,J)
        UIB(nx1,J)=UIB(2,J)
        VIB(nx1,J)=VIB(2,J)
        GWATX(nx1,J)=GWATX(2,J)
        GWATY(nx1,J)=GWATY(2,J)
        PGFUB(nx1,J)=PGFUB(2,J)
        PGFVB(nx1,J)=PGFVB(2,J)
      END DO

c**** interpolate air stress from A grid in atmos, to B grid in ice
C**** change of unit from change of momentum, to flux
C**** Update halo for USI,UOSURF,PGFU
      CALL HALO_UPDATE(grid, DMUA  , from=NORTH    )
      CALL HALO_UPDATE(grid, DMVA  , from=NORTH    )

      do j=j_0,j_1s
        im1=im
        do i=1,im
          GAIRX(i,j)=0.25*(dmua(i,j,2)+dmua(im1,j,2)+dmua(im1,j+1,2)
     &         +dmua(i,j+1,2))*bydts
          GAIRY(i,j)=0.25*(dmva(i,j,2)+dmva(im1,j,2)+dmva(im1,j+1,2)
     &         +dmva(i,j+1,2))*bydts
          im1=i
        enddo
      enddo
      IF (grid%HAVE_NORTH_POLE) THEN
        GAIRX(1:nx1,jm)=dmua(1,jm,2)*bydts
        GAIRY(1:nx1,jm)=dmva(1,jm,2)*bydts
      END IF
      do j=J_0NXY,J_1NXYS
       GAIRX(nx1-1,j)=GAIRX(1,j)
       GAIRY(nx1-1,j)=GAIRY(1,j)
       GAIRX(nx1,j)=GAIRX(2,j)
       GAIRY(nx1,j)=GAIRY(2,j)
      enddo

c**** read in sea ice velocity
      DO J=J_0NXY,J_1NXY
      DO I=1,NX1
       UICE(I,J,1)=UIB(I,J)
       VICE(I,J,1)=VIB(I,J)
       UICE(I,J,2)=0.
       VICE(I,J,2)=0.
       UICE(I,J,3)=0.
       VICE(I,J,3)=0.
      END DO
      END DO

C**** do the looping over pseudo-timesteps
      CALL VPICEDYN

C**** Calculate stress on ice velocity grid (B grid)
      DO J=J_0NXY,J_1NXY
        hemi=-1.
        if (J.gt.NY1/2) hemi=1.
        DO I=1,NX1-1
          DMU(i,j)=DTS*dwatn(i,j)*(COSWAT*(UICE(i,j,1)-GWATX(i,j))-
     *         HEMI*SINWAT*(VICE(i,j,1)-GWATY(i,j)))
          DMV(i,j)=DTS*dwatn(i,j)*(HEMI*SINWAT*(UICE(i,j,1)-GWATX(i,j))
     *         +COSWAT*(VICE(i,j,1)-GWATY(i,j)))
        END DO
      END DO
      DO J=J_0NXY,J_1NXY
        DMU(1,J)=DMU(NX1-1,J)
        DMU(NX1,J)=DMU(2,J)
      END DO

C**** interpolate ice velocity and stress from B grid to C grid in atm
C**** Update halos for UICE and DMU
      CALL HALO_UPDATE(grid,  UICE, from=SOUTH     )
      CALL HALO_UPDATE(grid,   DMU, from=SOUTH     )
      do j=J_0STG,J_1STG
        i=im
        do ip1=1,im
          usi(i,j)=0.5*(uice(i+1,j-1,1)+uice(i+1,j,1))
          IF (abs(USI(I,J)).lt.1d-10) USI(I,J)=0
          DMUI(I,J) = 0.5*(dmu(i+1,j-1)+dmu(i+1,j))
C**** Rescale DMUI to be net momentum into ocean
          DMUI(I,J) = 0.5*DMUI(I,J)*(FOCEAN(I,J)*RSI(I,J)+FOCEAN(IP1,J)
     *         *RSI(IP1,J))
          i=ip1
        enddo
      enddo
C**** Update halos for FOCEAN, DXYS, DXYV
      CALL HALO_UPDATE(grid,  FOCEAN, from=NORTH     )

      do j=j_0,j_1s
        do i=1,im
          vsi(i,j)=0.5*(vice(i,j,1)+vice(i+1,j,1))
          IF (abs(VSI(I,J)).lt.1d-10) VSI(I,J)=0
          DMVI(I,J) = 0.5*(dmv(i,j)+dmv(i+1,j))
C**** Rescale DMVI to be net momentum into ocean
          IF (J.lt.JM-1) DMVI(I,J) = 0.5*DMVI(I,J)*(FOCEAN(I,J)*RSI(I,J)
     *         *DXYN(J)+FOCEAN(I,J+1)*RSI(I,J+1)*DXYS(J+1))/DXYV(J+1)
          IF (J.eq.JM-1) DMVI(I,JM-1) = 0.5*DMVI(I,JM-1)*(FOCEAN(I,JM-1)
     *         *RSI(I,JM-1)*DXYN(JM-1)+FOCEAN(1,JM)*RSI(1,JM)*DXYS(JM))
     *         /DXYV(JM)
        enddo
      enddo
      IF (grid%HAVE_SOUTH_POLE) then
        usi(1:im,1)=0.
        dmui(1:im,1)=0.
      END IF
      IF (grid%HAVE_NORTH_POLE) then
        vsi(1:im,jm)=0.
        dmvi(1:im,jm)=0.
      END IF

C**** Calculate ustar*2*rho for ice-ocean fluxes on atmosphere grid
C**** UI2rho = | tau |
C**** Update halos for DMU, and DMV
      CALL HALO_UPDATE(grid,  DMU   , from=SOUTH     )
      CALL HALO_UPDATE(grid,  DMV   , from=SOUTH     )

      do j=J_0,J_1
        do i=1,imaxj(j)
          UI2rho(i,j)=0
          if (FOCEAN(I,J)*RSI(i,j).gt.0) THEN
C**** calculate 4 point average of B grid values of stresses
            duA = 0.5*(DXYN(J)*(dmu(i+1,j)+dmu(i,j))+DXYS(j)*(dmu(i+1
     *           ,j-1)+dmu(i,j-1)))*BYDXYP(J)
            dvA = 0.5*(DXYN(J)*(dmv(i+1,j)+dmv(i,j))+DXYS(j)*(dmv(i+1
     *           ,j-1)+dmv(i,j-1)))*BYDXYP(J)
            UI2rho(i,j)= sqrt (duA**2 + dvA**2) * bydts
          end if
        end do
      end do

C**** set north pole in C grid (atm)
      IF (grid%HAVE_NORTH_POLE) THEN
        USINP=0.
        DMUINP=0.
        do i=1,im
          USINP = USINP + USI(i,jm)
          DMUINP = DMUINP + DMUI(i,jm)
        enddo
        USINP=USINP/IM
        DMUINP=DMUINP/IM
        USI(1:im,jm)=USINP
        DMUI(1:im,jm)=DMUINP
        VSI(1:im,jm)=0.
        DMVI(1:im,jm)=0.
      END IF

C**** calculate mass fluxes for the ice advection
C**** Update halos for FOCEAN, and RSI
      CALL HALO_UPDATE(grid,  FOCEAN, from=NORTH     )
      CALL HALO_UPDATE(grid,  RSI   , from=NORTH     )

      DO J=J_0,J_1S
        I=IM
        DO IP1=1,IM
          USIDT(I,J)=0.
          IF (FOCEAN(I,J).gt.0 .and. FOCEAN(IP1,J).gt.0. .and.
     *         RSI(I,J)+RSI(IP1,J).gt.1d-4) THEN
            USIDT(I,J)=USI(I,J)*DTS
            ICIJ(I,J,IJ_USI) =ICIJ(I,J,IJ_USI) +(RSI(I,J)+RSI(IP1,J))
     *           *USI(i,j)
            ICIJ(I,J,IJ_DMUI)=ICIJ(I,J,IJ_DMUI)+DMUI(i,j)
          END IF
          VSIDT(I,J)=0.
          IF (FOCEAN(I,J+1).gt.0 .and. FOCEAN(I,J).gt.0. .and.
     *         RSI(I,J)+RSI(I,J+1).gt.1d-4) THEN
            VSIDT(I,J)=VSI(I,J)*DTS
            ICIJ(I,J,IJ_VSI) =ICIJ(I,J,IJ_VSI) +(RSI(I,J)+RSI(I,J+1))
     *           *VSI(i,j)
            ICIJ(I,J,IJ_DMVI)=ICIJ(I,J,IJ_DMVI)+DMVI(i,j)
          END IF
          ICIJ(I,J,IJ_PICE)=ICIJ(I,J,IJ_PICE)+ RSI(I,J)*press(i+1,j)
          ICIJ(I,J,IJ_RSI) =ICIJ(I,J,IJ_RSI) + RSI(I,J)
          I=IP1
        END DO
      END DO
      IF (grid%HAVE_NORTH_POLE) THEN
        VSIDT(1:IM,JM)=0.
        USIDT(1:IM,JM)=USI(1,JM)*DTS
        ICIJ(1,JM,IJ_USI) =ICIJ(1,JM,IJ_USI) +RSI(1,JM)*USI(1,JM)
        ICIJ(1,JM,IJ_DMUI)=ICIJ(1,JM,IJ_DMUI)+DMUI(1,JM)
        ICIJ(1,JM,IJ_RSI) =ICIJ(1,JM,IJ_RSI) +RSI(1,JM)
        ICIJ(1,JM,IJ_PICE)=ICIJ(1,JM,IJ_PICE)+RSI(1,JM)*press(1,JM)
      END IF
C**** Set uisurf,visurf for use in atmospheric drag calculations
      uisurf=usi
      visurf=vsi
C****
      END SUBROUTINE DYNSI

      SUBROUTINE ADVSI
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
      USE CONSTANT, only : byshi,lhm,grav
      USE MODEL_COM, only : im,jm,focean,p,ptop,kocean
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE, SOUTH, NORTH
      USE DOMAIN_DECOMP, only : HALO_UPDATE_COLUMN
      USE GEOM, only : dxyp,dyp,dxp,dxv,bydxyp,imaxj
c      USE ICEGEOM, only : dxyp,dyp,dxp,dxv,bydxyp ?????
      USE ICEDYN_COM, only : usidt,vsidt,rsix,rsiy,rsisave,icij,ij_musi
     *     ,ij_mvsi,ij_husi,ij_hvsi,ij_susi,ij_svsi
#ifdef TRACERS_WATER
     *     ,ticij,ticij_tusi,ticij_tvsi
#endif
      USE ICEDYN, only : grid_MIC
      USE SEAICE, only : ace1i,xsi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,lmi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,apress,msicnv,fwsim
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DIAG_COM, only : oa
      IMPLICIT NONE
!      REAL*8, DIMENSION(IM) :: FAW,FASI,FXSI,FYSI
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     FASI, FXSI, FYSI, FAW
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=2+2*LMI
#else
      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
      INTEGER ITR
      REAL*8 TRSNOW(NTM), TRICE(NTM)
#endif
      REAL*8 FMSI(NTRICE,IM),SFMSI(NTRICE),AMSI(NTRICE)
      REAL*8 FMSJ(IM,NTRICE,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 BYFOA(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cmsi1(grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cmsi2(grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cssi1(grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cssi2(grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER I,J,L,IM1,IP1,K
      REAL*8 SFASI,DMHSI,ASI,YRSI,XRSI,FRSI,SICE,TMP
!@var MHS mass/heat/salt content of sea ice
      REAL*8 MHS(NTRICE,IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
C****
C**** FLUXCB  USIDT  U compon of time integrated sea ice velocity (m)
C****         VSIDT  V compon of time integrated sea ice velocity (m)
C****
C****         FAW    flux of surface water area (m^2) = USIDT*DYP
C****         FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C****         FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)

      INTEGER J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_NORTH_POLE
C**** Get grid parameters
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S ,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO =J_1H ,
     &          HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Regularise ice concentration gradients to prevent advection errors
      DO J=J_0S, J_1S
      DO I=1,IM
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)

          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
      END DO
      END DO

C**** update RSISAVE for diagnostics
      RSISAVE(:,:)=RSI(:,:)

C**** set up local MHS array to contain all advected quantities
C**** MHS(1:2) = MASS, MHS(3:2+LMI) = HEAT, MHS(3+LMI:2+2*LMI)=SALT
C**** Currently this is on atmospheric grid
      MHS(1,:,J_0:J_1) = ACE1I + SNOWI(:,J_0:J_1)
      MHS(2,:,J_0:J_1) = MSI(:,J_0:J_1)
      DO L=1,LMI
        MHS(L+2,:,J_0:J_1) = HSI(L,:,J_0:J_1)
        MHS(L+2+LMI,:,J_0:J_1) = SSI(L,:,J_0:J_1)
      END DO
#ifdef TRACERS_WATER
C**** add tracers to advected arrays
      DO J=J_0, J_1
        DO I=1,IM
          DO ITR=1,NTM
          IF (SNOWI(I,J)*XSI(2).gt.XSI(1)*ACE1I) THEN ! layer 1:all snow
            SICE=SSI(1,I,J)+SSI(2,I,J)
            TRSNOW(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J)*MAX(1.
     *           -(ACE1I-SICE)/(XSI(2)*(ACE1I+SNOWI(I,J))-SICE),0d0)
          ELSE                  ! first layer is snow and some ice
            TRSNOW(ITR) = TRSI(ITR,1,I,J)*MIN(SNOWI(I,J)/(XSI(1)*(ACE1I
     *           +SNOWI(I,J))-SSI(1,I,J)),1d0)
          END IF
          TRICE(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J) - TRSNOW(ITR)
          MHS(1+2+(1+ITR)*LMI,I,J)=TRSNOW(ITR)
          MHS(2+2+(1+ITR)*LMI,I,J)=TRICE(ITR)
          DO L=3,LMI
            MHS(L+2+(1+ITR)*LMI,I,J)=TRSI(ITR,L,I,J)
          END DO
          END DO
        END DO
      END DO
#endif
C**** define inverse area array
      DO J=J_0, J_1
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          BYFOA(I,J)=BYDXYP(J)/FOCEAN(I,J)
        ELSE
          BYFOA(I,J)=0.
        END IF
      END DO
      END DO
C****
C**** North-South Advection of Sea Ice
C****
      SFASI  = 0.
      SFMSI(1:NTRICE) = 0.
C**** Update halo of DXV,RSIY,RSI,RSIX,FOCEAN,BYDXYP,and MHS
      CALL HALO_UPDATE(grid, RSI  , FROM=NORTH)
      CALL HALO_UPDATE(grid,FOCEAN, FROM=NORTH+SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, MHS  , FROM=NORTH)

      CALL HALO_UPDATE(grid, VSIDT, FROM=SOUTH)

      CALL HALO_UPDATE(grid_MIC, RSIY, FROM=NORTH)
      CALL HALO_UPDATE(grid_MIC, RSIX, FROM=NORTH)

      DO I=1,IM
C****
C**** Calculate south-north sea ice fluxes at grid box edges
C****
      DO 120 J=J_0S,MIN(JM-2,J_1)
      IF(VSIDT(I,J).eq.0.)  GO TO 120
      FAW(I,J) = VSIDT(I,J)*DXV(J+1) ! be careful with atm. grid index
      IF(VSIDT(I,J).le.0.) THEN
C**** Sea ice velocity is southward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J+1)-
     *       (1d0+FAW(I,J)*BYDXYP(J+1))*RSIY(I,J+1))*FOCEAN(I,J+1)
        FXSI(I,J)=FAW(I,J)*RSIX(I,J+1)*FOCEAN(I,J+1)
        FYSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J+1)*
     *       FAW(I,J)*RSIY(I,J+1)*FOCEAN(I,J+1) - 3d0*FASI(I,J))
        FMSJ(I,1:NTRICE,J) = FASI(I,J)*MHS(1:NTRICE,I,J+1)
      ELSE
C**** Sea ice velocity is northward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J)+(1d0-FAW(I,J)*BYDXYP(J))*RSIY(I,J))
     *       *FOCEAN(I,J)
        FXSI(I,J)=FAW(I,J)*RSIX(I,J)*FOCEAN(I,J)
        FYSI(I,J)=FAW(I,J)*
     *       (FAW(I,J)*BYDXYP(J)*FAW(I,J)*RSIY(I,J)*FOCEAN(I,J)
     *       -3d0*FASI(I,J))
        FMSJ(I,1:NTRICE,J) = FASI(I,J)*MHS(1:NTRICE,I,J)
      END IF
        ICIJ(I,J,IJ_MVSI)=ICIJ(I,J,IJ_MVSI)+SUM(FMSJ(I,1:2,J))
        ICIJ(I,J,IJ_HVSI)=ICIJ(I,J,IJ_HVSI)+SUM(FMSJ(I,3:2+LMI,J))
        ICIJ(I,J,IJ_SVSI)=ICIJ(I,J,IJ_SVSI)+SUM(FMSJ(I,3+LMI:2+2*LMI,J))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,J,TICIJ_TVSI,ITR)=TICIJ(I,J,TICIJ_TVSI,ITR)+
     *          SUM(FMSJ(I,3+(1+ITR)*LMI:2+(2+ITR)*LMI,J))
         END DO
#endif
  120 CONTINUE
C****
C**** Calculate south-north sea ice fluxes near North Pole
C****
      IF (HAVE_NORTH_POLE) THEN
        IF(VSIDT(I,JM-1).eq.0.)  EXIT
        FAW(I,JM-1) = VSIDT(I,JM-1)*DXV(JM) ! careful with atm.grid index!
        IF(VSIDT(I,JM-1).le.0.) THEN
C**** Sea ice velocity is southward from North Pole box
          FASI(I,JM-1) = FAW(I,JM-1)*RSI(1,JM)*FOCEAN(1,JM)
          FXSI(I,JM-1) = 0.
          FYSI(I,JM-1) = -FAW(I,JM-1)*FASI(I,JM-1)
          FMSJ(I,1:NTRICE,JM-1) = FASI(I,JM-1)*MHS(1:NTRICE,1,JM)
        ELSE
C**** Sea ice velocity is northward into North Pole box
          FASI(I,JM-1) = FAW(I,JM-1)*FOCEAN(I,JM-1)*
     *       (RSI(I,JM-1)+(1d0-FAW(I,JM-1)*BYDXYP(JM-1))*RSIY(I,JM-1))
          FXSI(I,JM-1) = FAW(I,JM-1)*RSIX(I,JM-1)*FOCEAN(I,JM-1)
          FYSI(I,JM-1) = FAW(I,JM-1)*(FAW(I,JM-1)*BYDXYP(JM-1)*
     *       FAW(I,JM-1)*RSIY(I,JM-1)*FOCEAN(I,JM-1)-3d0*FASI(I,JM-1))
          FMSJ(I,1:NTRICE,JM-1) = FASI(I,JM-1)*MHS(1:NTRICE,I,JM-1)
        END IF
C**** Accumulate sea ice leaving and entering North Pole box
        SFASI = SFASI + FASI(I,JM-1)
        SFMSI(1:NTRICE) = SFMSI(1:NTRICE) + FMSJ(I,1:NTRICE,JM-1)
         ICIJ(I,JM-1,IJ_MVSI)=ICIJ(I,JM-1,IJ_MVSI)+SUM(FMSJ(I,1:2,JM-1))
         ICIJ(I,JM-1,IJ_HVSI)=ICIJ(I,JM-1,IJ_HVSI)
     &                        +SUM(FMSJ(I,3:2+LMI,JM-1))
         ICIJ(I,JM-1,IJ_SVSI)=ICIJ(I,JM-1,IJ_SVSI)+
     *      SUM(FMSJ(I,3+LMI:2+2*LMI,JM-1))
#ifdef TRACERS_WATER
           DO ITR=1,NTM
             TICIJ(I,JM-1,TICIJ_TVSI,ITR)=TICIJ(I,JM-1,TICIJ_TVSI,ITR)+
     *            SUM(FMSJ(I,3+(1+ITR)*LMI:2+(2+ITR)*LMI,JM-1))
           END DO
#endif
      ENDIF

      END DO ! I loop
C****
C**** Update sea ice variables due to south-north fluxes
C****

C****Update halo of VSIDT, FASI, FAW, FOCEAN, FMSI, and FXSI
      CALL HALO_UPDATE_COLUMN(grid, FMSJ, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FASI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FXSI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FYSI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FAW,  FROM=SOUTH)

      DO I=1,IM
      DO 330 J=J_0S, J_1S
      IF(VSIDT(I,J-1)) 240,210,280
C**** VSIDT(J-1)=0.
  210 IF(VSIDT(I,J))  220,330,230
C**** VSIDT(J-1)=0, VSIDT(J)<0.
  220 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I,J)
      DO 225 K=1,NTRICE
  225 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSJ(I,K,J)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FYSI(I,J)
     *  + 3d0*(FAW(I,J)*ASI-DXYP(J)*FASI(I,J))) / (DXYP(J)-FAW(I,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) - FXSI(I,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)=0, VSIDT(J)>0.
  230 RSI(I,J)  =  RSI(I,J) -  FASI(I,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I,J)*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I,J)*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)<0.
  240 IF(VSIDT(I,J))  260,250,270
C**** VSIDT(J-1)<0, VSIDT(J)=0.
  250 RSI(I,J)  =  RSI(I,J) +  FASI(I,J-1)*BYFOA(I,J)
      TMP = (1d0+FAW(I,J-1)*FOCEAN(I,J-1)*BYFOA(I,J))
      RSIX(I,J) = RSIX(I,J)*TMP
      RSIY(I,J) = RSIY(I,J)*TMP**2
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)<0  or  VSIDT(J-1)>0, VSIDT(J) not 0.
  260 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + ( FASI(I,J-1)- FASI(I,J))
      DO 265 K=1,NTRICE
  265 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSJ(I,K,J-1)-FMSJ(I,K,J))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+
     *     (FYSI(I,J-1)-FYSI(I,J)) + 3d0*((FAW(I,J-1)+
     *     FAW(I,J))*ASI-DXYP(J)*(FASI(I,J-1)+FASI(I,J))))
     *    / (DXYP(J) + (FAW(I,J-1)-FAW(I,J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + (FXSI(I,J-1)-FXSI(I,J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)>0.
  270 RSI(I,J)  =  RSI(I,J) + (FASI(I,J-1)-FASI(I,J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(I,J-1)*FOCEAN(I,J-1)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(I,J-1)*FOCEAN(I,J-1)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))**2
      GO TO 310
C**** VSIDT(J-1)>0.
  280 IF(VSIDT(I,J).ne.0.)  GO TO 260
C**** VSIDT(J-1)>0, VSIDT(J)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(I,J-1)
      DO 285 K=1,NTRICE
  285 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSJ(I,K,J-1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FYSI(I,J-1)
     *    + 3d0*(FAW(I,J-1)*ASI-DXYP(J)*FASI(I,J-1))) /
     *     (DXYP(J)+FAW(I,J-1))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + FXSI(I,J-1)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  310 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 330
C**** Sea ice crunches into itself and completely covers grid box
  320 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYFOA(I,J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(I,J) -1d0 / ASI )
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C**** End of loop over J
  330 CONTINUE
C**** End of loop over I
      END DO
C****
C**** Advection of Sea Ice leaving and entering North Pole box
C****
      IF (HAVE_NORTH_POLE) THEN
        ASI = RSI(1,JM)*DXYP(JM)*FOCEAN(1,JM) + SFASI/IM
        DO 345 K=1,NTRICE
  345   AMSI(K) = RSI(1,JM)*DXYP(JM)*MHS(K,1,JM)*FOCEAN(1,JM)+
     &            SFMSI(K)/IM
        IF(ASI.gt.DXYP(JM)*FOCEAN(1,JM))  GO TO 350
        RSI(1,JM)   = ASI*BYFOA(1,JM)
        IF (ASI.gt.0) MHS(1:NTRICE,1,JM) = AMSI(1:NTRICE)/ASI
        GO TO 400
C**** Sea ice crunches into itself at North Pole box
  350   RSI(1,JM)   = 1d0
        MHS(1,1,JM) = AMSI(1)/ASI
        MHS(2,1,JM) =(AMSI(1)+AMSI(2))*BYFOA(1,JM)-MHS(1,1,JM)
        DO K=1,(NTRICE-2)/LMI
          MHS(3+LMI*(K-1),1,JM) = AMSI(3+LMI*(K-1)) / ASI
          MHS(4+LMI*(K-1),1,JM) = AMSI(4+LMI*(K-1)) / ASI
          DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+
     &              AMSI(5+LMI*(K-1))
     *         +AMSI(6+LMI*(K-1)))*(BYFOA(1,JM) -1d0/ ASI)
          MHS(5+LMI*(K-1),1,JM) = AMSI(5+LMI*(K-1)) / ASI +
     *         XSI(3)*DMHSI
          MHS(6+LMI*(K-1),1,JM) = AMSI(6+LMI*(K-1)) / ASI +
     *         XSI(4)*DMHSI
        END DO
      END IF   !HAVE_NORTH_POLE
C****
C**** East-West Advection of Sea Ice
C****
  400 DO 640 J=J_0S, J_1S
C****
C**** Calculate west-east sea ice fluxes at grid box edges
C****
      I=IM
      DO IP1=1,IM
      IF(USIDT(I,J).eq.0.)  GO TO 420
      FAW(I,J) = USIDT(I,J)*DYP(J)
      IF(USIDT(I,J).le.0.) THEN
C**** Sea ice velocity is westward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(IP1,J)-
     *       (1d0+FAW(I,J)*BYDXYP(J))*RSIX(IP1,J))
     *       *FOCEAN(IP1,J)
        FXSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J)*
     *       FAW(I,J)*RSIX(IP1,J)*FOCEAN(IP1,J)-3d0*FASI(I,J))
        FYSI(I,J)=FAW(I,J)*RSIY(IP1,J)*FOCEAN(IP1,J)
        FMSI(1:NTRICE,I) = FASI(I,J)*MHS(1:NTRICE,IP1,J)
      ELSE
C**** Sea ice velocity is eastward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J)+(1d0-FAW(I,J)*BYDXYP(J))*RSIX(I,J))
     *       *FOCEAN(I,J)
        FXSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J)*FAW(I,J)*
     *       RSIX(I,J)*FOCEAN(I,J)-3d0*FASI(I,J))
        FYSI(I,J)=FAW(I,J)*RSIY(I,J)*FOCEAN(I,J)
        FMSI(1:NTRICE,I) = FASI(I,J)*MHS(1:NTRICE,I,J)
      END IF
         ICIJ(I,J,IJ_MUSI)=ICIJ(I,J,IJ_MUSI)+SUM(FMSI(1:2,I))
         ICIJ(I,J,IJ_HUSI)=ICIJ(I,J,IJ_HUSI)+SUM(FMSI(3:2+LMI,I))
         ICIJ(I,J,IJ_SUSI)=ICIJ(I,J,IJ_SUSI)+SUM(FMSI(3+LMI:2+2*LMI,I))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,J,TICIJ_TUSI,ITR)=TICIJ(I,J,TICIJ_TUSI,ITR)+
     *          SUM(FMSI(3+(1+ITR)*LMI:2+(2+ITR)*LMI,I))
         END DO
#endif
  420 I=IP1
      END DO
C****
C**** Update sea ice variables due to west-east fluxes
C****
      IM1=IM
      DO 630 I=1,IM
      IF(USIDT(IM1,J)) 540,510,580
C**** USIDT(IM1)=0.
  510 IF(USIDT(I,J))  520,630,530
C**** USIDT(IM1)=0, USIDT(I)<0.
  520 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I,J)
      DO 525 K=1,NTRICE
  525 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSI(K,I)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FXSI(I,J)
     *    + 3d0*(FAW(I,J)*ASI-DXYP(J)*FASI(I,J))) / (DXYP(J)-FAW(I,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) - FYSI(I,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)=0, USIDT(I)>0.
  530 RSI(I,J)  =  RSI(I,J) -  FASI(I,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I,J)*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I,J)*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)<0.
  540 IF(USIDT(I,J))  560,550,570
C**** USIDT(IM1)<0, USIDT(I)=0.
  550 RSI(I,J)  =  RSI(I,J) +  FASI(IM1,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(IM1,J)*FOCEAN(IM1,J)*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(IM1,J)*FOCEAN(IM1,J)*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)<0  or  USIDT(IM1)>0, USIDT(I) not 0.
  560 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + (FASI(IM1,J)- FASI(I,J))
      DO 565 K=1,NTRICE
  565 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSI(K,IM1)-FMSI(K,I))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+
     *     (FXSI(IM1,J)-FXSI(I,J)) + 3d0*((FAW(IM1,J)+FAW(I,J))*ASI-
     *     DXYP(J)*(FASI(IM1,J)+FASI(I,J))))
     *    / (DXYP(J) + (FAW(IM1,J)-FAW(I,J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + (FYSI(IM1,J)-FYSI(I,J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)>0.
  570 RSI(I,J)  =  RSI(I,J) + (FASI(IM1,J)-FASI(I,J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(IM1,J)*FOCEAN(IM1,J)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(IM1,J)*FOCEAN(IM1,J)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)>0.
  580 IF(USIDT(I,J).ne.0.)  GO TO 560
C**** USIDT(IM1)>0, USIDT(I)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(IM1,J)
      DO 585 K=1,NTRICE
  585 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSI(K,IM1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FXSI(IM1,J)
     *    + 3d0*(FAW(IM1,J)*ASI-DXYP(J)*FASI(IM1,J))) /
     *     (DXYP(J)+FAW(IM1,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + FYSI(IM1,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  610 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 630
C**** Sea ice crunches into itself and completely covers grid box
  620 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYFOA(I,J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(I,J) -1d0/ ASI)
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C**** End of loop over I
  630 IM1=I
C**** End of loop over J
  640 CONTINUE

      IF (KOCEAN.ge.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=J_0, J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
C**** Fresh water sea ice mass convergence (needed for qflux model)
            MSICNV(I,J) = RSI(I,J)*(MHS(1,I,J)+MHS(2,I,J)-SUM(MHS(3
     *           +LMI:2*LMI+2,I,J))) - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     *           +MSI(I,J)-SUM(SSI(1:LMI,I,J)))
C**** sea ice prognostic variables
            SNOWI(I,J)= MAX(0d0,MHS(1,I,J) - ACE1I)
            MSI(I,J)  = MHS(2,I,J)
            DO L=1,LMI
              HSI(L,I,J) = MHS(L+2,I,J)
            END DO
C**** ensure that salinity is only associated with ice
            SICE=MHS(1+2+LMI,I,J)+MHS(2+2+LMI,I,J)
            IF (SNOWI(I,J).gt.XSI(2)*(ACE1I+SNOWI(I,J))) THEN
              SSI(1,I,J)=0.
            ELSE
              SSI(1,I,J)=SICE*(XSI(1)*ACE1I-XSI(2)*SNOWI(I,J))/ACE1I
            END IF
            SSI(2,I,J)=SICE-SSI(1,I,J)
C**** correction of heat energy to compensate for salinity fix 
            HSI(1,I,J)=HSI(1,I,J)-(MHS(1+2+LMI,I,J)-SSI(1,I,J))*LHM
            HSI(2,I,J)=HSI(2,I,J)+(MHS(1+2+LMI,I,J)-SSI(1,I,J))*LHM
            DO L=3,LMI
               SSI(L,I,J) = MHS(L+2+LMI,I,J)
            END DO
#ifdef TRACERS_WATER
C**** reconstruct tracer arrays
            DO ITR=1,NTM
              IF (ACE1I.gt.XSI(2)*(SNOWI(I,J)+ACE1I)) THEN
                TRSI(ITR,1,I,J)= MHS(2+2+(1+ITR)*LMI,I,J) *(ACE1I
     *               -XSI(2)*(SNOWI(I,J)+ACE1I))/ACE1I +MHS(1+2+(1+ITR)
     *               *LMI,I,J)
              ELSE
                TRSI(ITR,1,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)*XSI(1)*(ACE1I
     *               +SNOWI(I,J))/SNOWI(I,J)
              END IF
              TRSI(ITR,2,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)+MHS(2+2+(1+ITR)
     *             *LMI,I,J)-TRSI(ITR,1,I,J)
              DO L=3,LMI
                TRSI(ITR,L,I,J)=MHS(L+2+(1+ITR)*LMI,I,J)
              END DO
            END DO
#endif
            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     *           SUM(SSI(1:LMI,I,J)))
            END IF
          END DO
        END DO
C**** Set atmospheric arrays
        DO J=J_0, J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
C**** set total atmopsheric pressure anomaly in case needed by ocean
              APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+RSI(I,J)
     *             *(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV
              GTEMP(1:2,2,I,J)=((HSI(1:2,I,J)-SSI(1:2,I,J)*LHM)/
     *             (XSI(1:2)*(SNOWI(I,J)+ACE1I))+LHM)*BYSHI
#ifdef TRACERS_WATER
              GTRACER(:,2,I,J)=TRSI(:,1,I,J)/(XSI(1)*MHS(1,I,J)
     *             -SSI(1,I,J))
#endif
            END IF
          END DO
        END DO
        IF (HAVE_NORTH_POLE) THEN
          IF (FOCEAN(1,JM).gt.0) THEN
            DO I=2,IM           ! North pole
              APRESS(I,JM)=APRESS(1,JM)
              GTEMP(1:2,2,I,JM)= GTEMP(1:2,2,1,JM)
#ifdef TRACERS_WATER
              GTRACER(:,2,I,JM)=GTRACER(:,2,1,JM)
#endif
            END DO
          END IF
        END IF
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=J_0, J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
              OA(I,J,13)=OA(I,J,13)+(RSI(I,J)*SUM(MHS(3:2+LMI,I,J))
     *             -RSISAVE(I,J)*SUM(HSI(1:LMI,I,J)))
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****
      RETURN
      END SUBROUTINE ADVSI

      SUBROUTINE init_icedyn(iniOCEAN)
!@sum  init_icedyn initializes ice dynamics variables
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,dtsrc,foceanA=>focean
      USE DIAG_COM, only : ia_src
      USE ICEDYN_COM
      USE ICEDYN, only : setup_icedyn_grid,focean,osurf_tilt
      USE FLUXES, only : uisurf,visurf
      USE PARAM
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER k,i,j

C**** setup ice dynamics grid
C**** Currently using ATM grid:
C**** -------------------------
C**** If a different grid (but still lat/lon) is required, edit
C**** definition of focean here and the definition of imic,jmic in
C**** ICEDYN.f, and obviously the code in DYNSI.
C**** -------------------------
C**** Set land masks for ice dynamics
      focean(:,:)=foceanA(:,:)         ! EDIT FOR GRID CHANGE

      call setup_icedyn_grid

C**** Initiallise ice dynamics if ocean model needs initialising
      if (iniOCEAN) THEN
        RSIX=0.
        RSIY=0.
        USI=0.
        VSI=0.
      end if

C**** set uisurf,visurf for atmopsherice drag calculations
      uisurf=usi
      visurf=vsi

C**** set properties for ICIJ diagnostics
      k=0

      k=k+1
      IJ_USI=k
      lname_icij(k)="Sea ice EW velocity x POICEU"
      sname_icij(k)="icij_usi"
      units_icij(k)="m/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1.
      ijgrid_icij(k)=2

      k=k+1
      IJ_VSI=k
      lname_icij(k)="Sea ice NS velocity x POICEV"
      sname_icij(k)="icij_vsi"
      units_icij(k)="m/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1.
      ijgrid_icij(k)=2

      k=k+1
      IJ_DMUI=k
      lname_icij(k)="Ice-ocean EW stress"
      sname_icij(k)="icij_dmui"
      units_icij(k)="kg/m s^2"
      ia_icij(k)=ia_src
      scale_icij(k)=1./dtsrc
      ijgrid_icij(k)=1

      k=k+1
      IJ_DMVI=k
      lname_icij(k)="Ice-ocean NS stress"
      sname_icij(k)="icij_dmvi"
      units_icij(k)="kg/m s^2"
      ia_icij(k)=ia_src
      scale_icij(k)=1./dtsrc
      ijgrid_icij(k)=1

      k=k+1
      IJ_PICE=k
      lname_icij(k)="Sea ice internal pressure x POICE"
      sname_icij(k)="icij_psi"
      units_icij(k)="10^3 kg/m s^2"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-3
      ijgrid_icij(k)=1

      k=k+1
      IJ_MUSI=k
      lname_icij(k)="Sea ice NS mass flux"
      sname_icij(k)="icij_musi"
      units_icij(k)="10^7 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-7/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_MVSI=k
      lname_icij(k)="Sea ice EW mass flux"
      sname_icij(k)="icij_mvsi"
      units_icij(k)="10^7 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-7/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_HUSI=k
      lname_icij(k)="Sea ice NS heat flux"
      sname_icij(k)="icij_husi"
      units_icij(k)="10^12 W"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-12/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_HVSI=k
      lname_icij(k)="Sea ice EW heat flux"
      sname_icij(k)="icij_hvsi"
      units_icij(k)="10^12 W"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-12/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_SUSI=k
      lname_icij(k)="Sea ice NS salt flux"
      sname_icij(k)="icij_susi"
      units_icij(k)="10^3 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-3/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_SVSI=k
      lname_icij(k)="Sea ice EW salt flux"
      sname_icij(k)="icij_svsi"
      units_icij(k)="10^3 kg/s"
      ia_icij(k)=ia_src
      scale_icij(k)=1d-3/dtsrc
      ijgrid_icij(k)=2

      k=k+1
      IJ_RSI=k
      lname_icij(k)="Ocean ice fraction (ice dynamic grid)"
      sname_icij(k)="icij_rsi"
      units_icij(k)="%"
      ia_icij(k)=ia_src
      scale_icij(k)=100.
      ijgrid_icij(k)=1

      if (k.gt.KICIJ) then
        write(6,*) "Too many ICIJ diags: increase KICIJ to at least"
     *       ,k
        call stop_model("ICIJ diagnostic error",255)
      end if

#ifdef TRACERS_WATER
C**** simple tracer diags same description for all tracers
C**** set properties for TICIJ diagnostics
      k=0

      k=k+1
      TICIJ_TUSI=k
      lname_ticij(k)="Sea ice NS tracer flux"
      sname_ticij(k)="ticij_tusi"
      units_ticij(k)="kg/s"
      ia_ticij(k)=ia_src
      scale_ticij(k)=1./dtsrc
      ijgrid_ticij(k)=2

      k=k+1
      TICIJ_TVSI=k
      lname_ticij(k)="Sea ice EW tracer flux"
      sname_ticij(k)="ticij_tvsi"
      units_ticij(k)="kg/s"
      ia_ticij(k)=ia_src
      scale_ticij(k)=1./dtsrc
      ijgrid_ticij(k)=2

      if (k.gt.KTICIJ) then
        write(6,*) "Too many TICIJ diags: increase KTICIJ to at least"
     *       ,k
        call stop_model("TICIJ diagnostic error",255)
      end if
#endif

      RETURN
      END SUBROUTINE init_icedyn

      SUBROUTINE diag_ICEDYN
!@sum  diag_ICEDYN prints out diagnostics for ice dynamics
!@&    ESMF: It should only be called from a serial region.
!@$          It is NOT parallelized
!@auth Gavin Schmidt
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm,trname,ntrocn
#endif
      USE DIAG_COM, only : qdiag,acc_period
      USE ICEDYN_COM
      USE ICEDYN, only : focean
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE
      REAL*8, DIMENSION(IMIC,JMIC) :: Q,ADENOM
      INTEGER I,J,L,N,KXLB,ijgrid,IP1,k1,k
      CHARACTER XLB*30
      CHARACTER TITLE*80,lname*50,sname*30,units*50
      character*50 :: unit_string
      REAL*8 QJ(JM),QSUM
      REAL*8 byiacc

      IF (.not. QDIAG) RETURN
C**** determine label to be added to all titles
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)

C**** Open output files
      call open_ij(trim(acc_period)//'.icij'//
     *     XLABEL(1:LRUNID),imic,jmic)
C****
C**** Simple scaled ICIJ diagnostics
C****
      DO K=1,KICIJ
        byiacc=1./(IDACC(IA_ICIJ(K))+teeny)
        adenom=1.
        lname=lname_icij(k)
        k1 = index(lname,' x ')
        if (k1 .gt. 0) then
          if (index(lname,' x POICE ') .gt. 0) then
            do j=1,jmic
            do i=1,imic
              adenom(i,j)=icij(i,j,ij_rsi) * byiacc
            end do
            end do
          else if (index(lname,' x POICEU') .gt. 0) then
            do j=1,jmic
            i=imic
            do ip1=1,imic
              adenom(i,j)=0.5*(icij(i,j,ij_rsi)+icij(ip1,j,ij_rsi))
     *             * byiacc
              i=ip1
            end do
            end do
          else if (index(lname,' x POICEV') .gt. 0) then
            do j=1,jmic-1
            do i=1,imic
              adenom(i,j)=0.5*(icij(i,j,ij_rsi)+icij(i,j+1,ij_rsi))
     *             * byiacc
            end do
            end do
          end if
          lname(k1:50) = ' '
        end if

        Q=UNDEF ; QJ=UNDEF ; QSUM=UNDEF
        DO J=1,JMIC
          DO I=1,IMIC
            IF (ADENOM(I,J).gt.0 .and. FOCEAN(I,J).gt.0.5)
     *           Q(I,J)=SCALE_ICIJ(K)*ICIJ(I,J,K)*byiacc/adenom(i,j)
          END DO
        END DO
        Q(2:IMIC,JMIC)=Q(1,JMIC)
        Q(2:IMIC,1)=Q(1,1)
        TITLE=trim(LNAME)//" ("//trim(UNITS_ICIJ(K))//") "
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME_ICIJ(K),LNAME_ICIJ(K),UNITS_ICIJ(K),Q
     *       ,QJ,QSUM,IJGRID_ICIJ(K),IJGRID_ICIJ(K)) ! assume igrid=jgrid

      END DO

#ifdef TRACERS_WATER
C**** simple tracer diags (no need for weighting)
C**** Name and scale are tracer dependent
      DO K=1,KTICIJ
        byiacc=1./(IDACC(IA_TICIJ(K))+teeny)

        DO N=1,NTM
        lname=trim(trname(n))//" "//lname_ticij(k)
        sname=trim(trname(n))//"_"//sname_ticij(k)
        Q=UNDEF
        DO J=1,JMIC
          DO I=1,IMIC
            IF (FOCEAN(I,J).gt.0.5) Q(I,J)=10**(-ntrocn(n))*
     *           SCALE_TICIJ(K)*TICIJ(I,J,K,N)*byiacc
          END DO
        END DO
        Q(2:IMIC,JMIC)=Q(1,JMIC)
        Q(2:IMIC,1)=Q(1,1)
        UNITS=unit_string(ntrocn(n),UNITS_TICIJ(K))
        TITLE=trim(LNAME)//" ("//trim(UNITS)//")"
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,
     &       IJGRID_TICIJ(K),IJGRID_TICIJ(K)) ! assume igrid=jgrid
      END DO
      END DO
#endif
      call close_ij
C****
      RETURN
      END SUBROUTINE diag_ICEDYN

