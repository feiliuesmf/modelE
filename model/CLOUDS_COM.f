#include "rundeck_opts.h"
      MODULE CLOUDS_COM
!@sum  CLOUDS_COM model variables for moist convction and
!@+    large-scale condensation
!@+    cloud droplet number added to list of saves for rsf files
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE MODEL_COM, only : IM,JM,LM
      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      SAVE
!@var TTOLD,QTOLD previous potential temperature, humidity
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RHSAV
C**** some arrays here for compatility with new clouds
!@var CLDSAV, CLDSAV1 previous cloud cover area (percent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLDSAV,CLDSAV1
!@var ULS,VLS,UMC,VMC velocity work arrays
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ULS,VLS,UMC,VMC
!@var TLS,QLS,TMC,QMC temperature and humidity work arrays
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TLS,QLS,TMC,QMC
!@var FSS grid fraction for large-scale clouds
!@+   initialised as 1. for compatibility with previous clouds
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FSS
#ifdef CLD_AER_CDNC
!@var OLDNO, OLDNL old CDNC for ocean and land ,SMFPM:CTEI parameter
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: OLDNO,OLDNL,SMFPM
#endif

C**** variables saved for radiation calculations
!@var TAUSS optical depth from super-saturated clouds
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAUSS
!@var TAUMC optical depth from moist-convective clouds
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAUMC
!@var CLDSS super-saturated cloud cover area (percent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLDMC
!@var CSIZMC,CSIZSS mc,ss effective cloud droplet radius (microns)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CSIZMC,CSIZSS

C**** variables saved for surface wind spectrum calculations
!@var DDM1 downdraft mass flux / rho at lowest level (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DDM1

C**** variables used (and saved) for gravity wave drag calculations
!@var AIRX, AIRMX*DXYP(J) convective mass flux (kg/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: AIRX
!@var LMC max layer of mc convective mass flux. 
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: LMC

!@var LLOW,LMID,LHI max levels for low, mid and high clouds
      INTEGER LLOW,LMID,LHI

      END MODULE CLOUDS_COM

      SUBROUTINE ALLOC_CLOUDS_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE MODEL_COM, ONLY : IM,LM
      USE CLOUDS_COM, ONLY : TTOLD,QTOLD,SVLHX,SVLAT,RHSAV,CLDSAV,
     *                       CLDSAV1,FSS,
#ifdef CLD_AER_CDNC
     *                       OLDNO,OLDNL,SMFPM,
#endif
     *                       TAUSS,TAUMC, CLDSS,CLDMC,CSIZMC,CSIZSS,
     *                       ULS,VLS,UMC,VMC,TLS,QLS,
     *                       TMC,QMC,DDM1,AIRX,LMC
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE(     TTOLD(LM,IM,J_0H:J_1H),
     *              QTOLD(LM,IM,J_0H:J_1H),
     *              SVLHX(LM,IM,J_0H:J_1H),
     *              SVLAT(LM,IM,J_0H:J_1H),
     *              RHSAV(LM,IM,J_0H:J_1H),
     *             CLDSAV(LM,IM,J_0H:J_1H),
     *            CLDSAV1(LM,IM,J_0H:J_1H),
     *                FSS(LM,IM,J_0H:J_1H),
#ifdef CLD_AER_CDNC
     *             OLDNO(LM,IM,J_0H:J_1H),
     *             OLDNL(LM,IM,J_0H:J_1H),
     *             SMFPM(LM,IM,J_0H:J_1H),
#endif
     *              TAUSS(LM,IM,J_0H:J_1H),
     *              TAUMC(LM,IM,J_0H:J_1H),
     *              CLDSS(LM,IM,J_0H:J_1H),
     *              CLDMC(LM,IM,J_0H:J_1H),
     *             CSIZMC(LM,IM,J_0H:J_1H),
     *             CSIZSS(LM,IM,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     ULS(IM,J_0H:J_1H,LM),
     *              VLS(IM,J_0H:J_1H,LM),
     *              UMC(IM,J_0H:J_1H,LM),
     *              VMC(IM,J_0H:J_1H,LM),
     *              TLS(IM,J_0H:J_1H,LM),
     *              QLS(IM,J_0H:J_1H,LM),
     *              TMC(IM,J_0H:J_1H,LM),
     *              QMC(IM,J_0H:J_1H,LM),
     *         STAT=IER)

!@var FSS initialized to 1.
      FSS = 1. 
#ifdef CLD_AER_CDNC
!@var OLDNO and OLDNL are initialized to 10.
      OLDNO = 10. 
      OLDNL = 10. 
!@var SMFPM is initialised to 0.5 (proxy for cloud top turbulence)
      SMFPM = 0.5
#endif 

      ALLOCATE(     DDM1(IM,J_0H:J_1H),
     *              AIRX(IM,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     LMC(2,IM,J_0H:J_1H),
     *         STAT=IER)

C**** Initialise some output used in dynamics
      LMC(:,:,J_0H:J_1H)=0
      AIRX(:,J_0H:J_1H)=0.

      END SUBROUTINE ALLOC_CLOUDS_COM

      SUBROUTINE io_clouds(kunit,iaction,ioerr)
!@sum  io_clouds reads and writes cloud arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE CLOUDS_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "CLD01"

      write(MODULE_HEADER(lhead+1:80),'(a)')
     *'R8 dim(im,jm,lm):potT,Hum,LatHeat,RHum,CldCv,NO,NL,SM (all old)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,
     *     TTOLD,QTOLD,SVLHX,RHSAV,CLDSAV
#ifdef CLD_AER_CDNC
     *     ,OLDNO,OLDNL,SMFPM
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,
     *     TTOLD,QTOLD,SVLHX,RHSAV,CLDSAV
#ifdef CLD_AER_CDNC
     *     ,OLDNO,OLDNL,SMFPM
#endif
        IF (HEADER(1:15).NE.MODULE_HEADER(1:15)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_clouds
