#include "rundeck_opts.h"

      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm
#endif
      USE SOCPBL, only : npbl=>n
      IMPLICIT NONE
      SAVE
!@var uabl boundary layer profile for zonal wind
!@var vabl boundary layer profile for meridional wind
!@var tabl boundary layer profile for temperature
!@var qabl boundary layer profile for humidity
      real*8, allocatable, dimension(:,:,:,:) :: uabl,vabl,tabl,qabl
!@var eabl boundary layer profile for turbulent KE (calc. on sec. grid)
      real*8, allocatable, dimension(:,:,:,:) :: eabl
#ifdef TRACERS_ON
!@var trabl boundary layer profile for tracers
      real*8, allocatable, dimension(:,:,:,:,:) :: trabl
#endif

!@var cmgs drag coefficient (dimensionless surface momentum flux)
!@var chgs Stanton number   (dimensionless surface heat flux)
!@var cqgs Dalton number    (dimensionless surface moisture flux)
      real*8, allocatable, dimension(:,:,:) :: cmgs,chgs,cqgs

!@var ipbl flag for whether pbl properties were found at last timestep
      integer, allocatable, dimension(:,:,:) :: ipbl

!@var ROUGHL log10(zgs/roughness length), prescribed with zgs=30 m.
      REAL*8, allocatable, dimension(:,:) :: roughl
!@var WSAVG     COMPOSITE SURFACE WIND MAGNITUDE (M/S)
!@var TSAVG     COMPOSITE SURFACE AIR TEMPERATURE (K)
!@var QSAVG     COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
!@var DCLEV     LAYER TO WHICH DRY CONVECTION MIXES (1)
!@var USAVG     COMPOSITE SURFACE U WIND
!@var VSAVG     COMPOSITE SURFACE V WIND
!@var TAUAVG    COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
!@var w2_l1     COMPOSITE vertical component of t.k.e. at gcm layer 1
!@var USTAR_pbl friction velocity (sqrt of srfc mom flux) (m/s)
      REAL*8, allocatable, dimension(:,:) ::
     &     wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,tgvavg,qgavg
     &    ,w2_l1
      REAL*8, allocatable, dimension(:,:,:) :: ustar_pbl

!@var egcm  3-d turbulent kinetic energy in the whole atmosphere
!@var w2gcm vertical component of egcm
!@var t2gcm 3-d turbulent temperature variance in the whole atmosphere
      real*8, allocatable, dimension(:,:,:) :: egcm,w2gcm,t2gcm

!@var uflux surface turbulent u-flux (=-<uw>)
!@var vflux surface turbulent v-flux (=-<vw>)
!@var tflux surface turbulent t-flux (=-<tw>)
!@var qflux surface turbulent q-flux (=-<qw>)
      real*8, allocatable, dimension(:,:) :: uflux,vflux,tflux,qflux

      END MODULE PBLCOM

      SUBROUTINE io_pbl(kunit,iaction,ioerr)
!@sum  io_pbl reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,irsfic,irerun,iowrite,irsficno,lhead
      USE PBLCOM
      USE DOMAIN_DECOMP, only : grid, GET, CHECKSUM, AM_I_ROOT
      USE DOMAIN_DECOMP, only : pack_column, pack_data
      USE DOMAIN_DECOMP, only : pack_block , unpack_block
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
      REAL*8, DIMENSION(npbl,IM,JM,4) :: uabl_glob,vabl_glob,tabl_glob,
     &     qabl_glob,eabl_glob
      REAL*8, DIMENSION(IM,JM,4) :: cmgs_glob, chgs_glob, cqgs_glob
      INTEGER, DIMENSION(IM,JM,4) :: ipbl_glob
      INTEGER :: J_0, J_1
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "PBL01"
#ifdef TRACERS_ON
!@var TR_HEADER Character string label for tracer record
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRPBL01"
      REAL*8, DIMENSION(npbl,ntm,im,jm,4) :: trabl_glob
      write (TR_MODULE_HEADER(lhead+1:80),'(a7,i2,a,i2,a)') 'R8 dim(',
     *     npbl,',',ntm,',ijm,4):TRt'
#endif
      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)') 'R8 dim(',npbl,
     *  ',ijm,4):Ut,Vt,Tt,Qt,Et dim(ijm,4,3):Cmhq, I:Ipb(ijm,4)'

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_COLUMN(grid, uabl, uabl_glob)
        CALL PACK_COLUMN(grid, vabl, vabl_glob)
        CALL PACK_COLUMN(grid, tabl, tabl_glob)
        CALL PACK_COLUMN(grid, qabl, qabl_glob)
        CALL PACK_COLUMN(grid, eabl, eabl_glob)

        CALL PACK_DATA(grid, cmgs, cmgs_glob)
        CALL PACK_DATA(grid, chgs, chgs_glob)
        CALL PACK_DATA(grid, cqgs, cqgs_glob)
        CALL PACK_DATA(grid, ipbl, ipbl_glob)

#ifdef TRACERS_ON
        CALL PACK_BLOCK(grid, trabl, trabl_glob)
#endif
        IF (AM_I_ROOT()) THEN
          WRITE (KUNIT,ERR=10) MODULE_HEADER,UABL_GLOB,VABL_GLOB
     *       ,TABL_GLOB,QABL_GLOB,EABL_GLOB,CMGS_GLOB
     *       ,CHGS_GLOB,CQGS_GLOB,IPBL_GLOB
#ifdef TRACERS_ON
          WRITE (KUNIT,ERR=10) TR_MODULE_HEADER,TRABL_GLOB
#endif
        END IF

      CASE (IOREAD:)            ! input from restart file or restart
        READ (KUNIT,ERR=10) HEADER,UABL_glob,VABL_glob,TABL_glob,
     &       QABL_glob,EABL_glob,CMGS_glob,CHGS_glob,CQGS_glob,IPBL_glob

        uabl(:,:,J_0:J_1,:) = uabl_glob(:,:,J_0:J_1,:)
        vabl(:,:,J_0:J_1,:) = vabl_glob(:,:,J_0:J_1,:)
        tabl(:,:,J_0:J_1,:) = tabl_glob(:,:,J_0:J_1,:)
        qabl(:,:,J_0:J_1,:) = qabl_glob(:,:,J_0:J_1,:)
        eabl(:,:,J_0:J_1,:) = eabl_glob(:,:,J_0:J_1,:)

        cmgs(:,J_0:J_1,:) = cmgs_glob(:,J_0:J_1,:)
        chgs(:,J_0:J_1,:) = chgs_glob(:,J_0:J_1,:)
        cqgs(:,J_0:J_1,:) = cqgs_glob(:,J_0:J_1,:)
        ipbl(:,J_0:J_1,:) = ipbl_glob(:,J_0:J_1,:)
        CALL CHECKSUM(grid,cmgs,__LINE__,__FILE__//'::cmgs')

        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_ON
        SELECT CASE (IACTION)
        CASE (IOREAD,IRERUN,IRSFIC,IRSFICNO)    ! restarts
          READ (KUNIT,ERR=10) TR_HEADER,TRABL_GLOB
          IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in tracer module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
            GO TO 10
          END IF
          CALL UNPACK_BLOCK(grid, TRABL_GLOB, TRABL, local=.true.)
        END SELECT
#endif
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_pbl

      SUBROUTINE io_bldat(kunit,iaction,ioerr)
!@sum  io_bldat reads and writes boundary layer data to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE DOMAIN_DECOMP, only : GET, grid, ARRAYGATHER, AM_I_ROOT
      USE DOMAIN_DECOMP, only : CHECKSUM, CHECKSUM_COLUMN
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "BLD02"
      INTEGER :: J_0, J_1, J_0H, J_1H, L,K
      REAL*8, DIMENSION(IM,JM) :: wsavg_glob,tsavg_glob,
     *       qsavg_glob,dclev_glob,usavg_glob,
     *       vsavg_glob,tauavg_glob,qgavg_glob,tgvavg_glob
      REAL*8, DIMENSION(LM,IM,JM) :: egcm_glob, w2gcm_glob
      REAL*8, DIMENSION(IM,JM,4) :: ustar_pbl_glob

      MODULE_HEADER(lhead+1:80) = 'R8 dim(ijm):ws,ts,qs,'//
     *  'LvlDC,us,vs,tau,u*(.,4),ke;w2(lijm),tgv,qg'

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL ARRAYGATHER(grid, wsavg, wsavg_glob)
        CALL ARRAYGATHER(grid, tsavg, tsavg_glob)
        CALL ARRAYGATHER(grid, qsavg, qsavg_glob)
        CALL ARRAYGATHER(grid, dclev, dclev_glob)
        CALL ARRAYGATHER(grid, usavg, usavg_glob)
        CALL ARRAYGATHER(grid, vsavg, vsavg_glob)
        CALL ARRAYGATHER(grid, tauavg, tauavg_glob)
        CALL ARRAYGATHER(grid, tgvavg, tgvavg_glob)
        CALL ARRAYGATHER(grid, qgavg(:,:), qgavg_glob(:,:))

        DO K=1,4
          CALL ARRAYGATHER(grid,ustar_pbl(:,:,K),ustar_pbl_glob(:,:,K))
        END DO
        DO L = 1, LM
          CALL ARRAYGATHER(grid, egcm(L,:,:), egcm_glob(L,:,:))
          CALL ARRAYGATHER(grid, w2gcm(L,:,:), w2gcm_glob(L,:,:))
        END DO
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,wsavg_glob,tsavg_glob
     *         ,qsavg_glob,dclev_glob,usavg_glob,vsavg_glob,tauavg_glob
     *         ,ustar_pbl_glob,egcm_glob,w2gcm_glob,tgvavg_glob
     *         ,qgavg_glob
        END IF

      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,wsavg_glob,tsavg_glob,
     *       qsavg_glob,dclev_glob,usavg_glob,
     *       vsavg_glob,tauavg_glob,ustar_pbl_glob,egcm_glob,
     *       w2gcm_glob,tgvavg_glob,qgavg_glob
        wsavg(:,J_0:J_1) = wsavg_glob(:,J_0:J_1)
        tsavg(:,J_0:J_1) = tsavg_glob(:,J_0:J_1)
        qsavg(:,J_0:J_1) = qsavg_glob(:,J_0:J_1)
        dclev(:,J_0:J_1) = dclev_glob(:,J_0:J_1)
        usavg(:,J_0:J_1) = usavg_glob(:,J_0:J_1)
        vsavg(:,J_0:J_1) = vsavg_glob(:,J_0:J_1)
        tauavg(:,J_0:J_1)= tauavg_glob(:,J_0:J_1)
        tgvavg(:,J_0:J_1)= tgvavg_glob(:,J_0:J_1)
        qgavg(:,J_0:J_1) = qgavg_glob(:,J_0:J_1)

        CALL CHECKSUM(grid, wsavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, tsavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, qsavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, dclev, __LINE__, __FILE__)
        CALL CHECKSUM(grid, usavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, vsavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, tauavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, tgvavg, __LINE__, __FILE__)
        CALL CHECKSUM(grid, qgavg, __LINE__, __FILE__)

        ustar_pbl(:,J_0:J_1,:) = ustar_pbl_glob(:,J_0:J_1,:)
        egcm(:,:,J_0:J_1)= egcm_glob(:,:,J_0:J_1)
        w2gcm(:,:,J_0:J_1)=w2gcm_glob(:,:,J_0:J_1)
        CALL CHECKSUM       (grid, ustar_pbl, __LINE__, __FILE__)
        CALL CHECKSUM_COLUMN(grid, egcm, __LINE__, __FILE__)
        CALL CHECKSUM_COLUMN(grid, w2gcm, __LINE__, __FILE__)

        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_bldat

      SUBROUTINE ALLOC_PBL_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE PBLCOM
      USE DOMAIN_DECOMP, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(    uabl(npbl,    im,J_0H:J_1H,4),
     *             vabl(npbl,    im,J_0H:J_1H,4),
     *             tabl(npbl,    im,J_0H:J_1H,4),
     *             qabl(npbl,    im,J_0H:J_1H,4),
     *             eabl(npbl,    im,J_0H:J_1H,4),
     *         STAT=IER)
      qabl=0.  ! initialise to make life easier

#ifdef TRACERS_ON
      ALLOCATE(    trabl(npbl,ntm,im,J_0H:J_1H,4),
     *         STAT=IER)
#endif

      ALLOCATE(    cmgs(im,J_0H:J_1H,4),
     *             chgs(im,J_0H:J_1H,4),
     *             cqgs(im,J_0H:J_1H,4),
     *             ipbl(im,J_0H:J_1H,4),
     *         STAT=IER)
      ipbl(:,J_0H:J_1H,:) = 0.

      ALLOCATE(    roughl(im,J_0H:J_1H),
     *              wsavg(im,J_0H:J_1H),
     *              tsavg(im,J_0H:J_1H),
     *              qsavg(im,J_0H:J_1H),
     *              dclev(im,J_0H:J_1H),
     *              usavg(im,J_0H:J_1H),
     *              vsavg(im,J_0H:J_1H),
     *             tauavg(im,J_0H:J_1H),
     *             tgvavg(im,J_0H:J_1H),
     *              qgavg(im,J_0H:J_1H),
     *              w2_l1(im,J_0H:J_1H),
     *         STAT=IER)
      w2_l1(:,J_0H:J_1H) =0.

      ALLOCATE(    ustar_pbl(im,J_0H:J_1H,4),
     *         STAT=IER)

      ALLOCATE(    egcm(lm,im,J_0H:J_1H),
     *            w2gcm(lm,im,J_0H:J_1H),
     *            t2gcm(lm,im,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(   uflux(im,J_0H:J_1H),
     *            vflux(im,J_0H:J_1H),
     *            tflux(im,J_0H:J_1H),
     *            qflux(im,J_0H:J_1H),
     *         STAT=IER)
      uflux(:,J_0H:J_1H) = 0. ! define defaults
      vflux(:,J_0H:J_1H) = 0. 
      tflux(:,J_0H:J_1H) = 0. 
      qflux(:,J_0H:J_1H) = 0. 

      END SUBROUTINE ALLOC_PBL_COM

