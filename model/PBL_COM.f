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
     &    ,w2_l1,gustiwind
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
      USE DOMAIN_DECOMP, only : unpack_column, unpack_data
      USE DOMAIN_DECOMP, only : pack_block , unpack_block
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
      REAL*8, DIMENSION(npbl,4,IM,JM) :: uabl_glob,vabl_glob,tabl_glob,
     &     qabl_glob,eabl_glob
      REAL*8, DIMENSION(4,IM,JM) :: cmgs_glob, chgs_glob, cqgs_glob
      INTEGER, DIMENSION(4,IM,JM) :: ipbl_glob
      INTEGER :: J_0, J_1
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "PBL01"
#ifdef TRACERS_ON
!@var TR_HEADER Character string label for tracer record
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRPBL01"
      REAL*8, DIMENSION(npbl,ntm,4,im,jm) :: trabl_glob
      write (TR_MODULE_HEADER(lhead+1:80),'(a7,i2,a,i2,a)') 'R8 dim(',
     *     npbl,',',ntm,',4,ijm):TRt'
#endif
      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)') 'R8 dim(',npbl,
     *  ',4,ijm):Ut,Vt,Tt,Qt,Et dim(4,ijm,3):Cmhq, I:Ipb(4,ijm)'

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_BLOCK(grid, uabl, uabl_glob)
        CALL PACK_BLOCK(grid, vabl, vabl_glob)
        CALL PACK_BLOCK(grid, tabl, tabl_glob)
        CALL PACK_BLOCK(grid, qabl, qabl_glob)
        CALL PACK_BLOCK(grid, eabl, eabl_glob)

        CALL PACK_COLUMN(grid, cmgs, cmgs_glob)
        CALL PACK_COLUMN(grid, chgs, chgs_glob)
        CALL PACK_COLUMN(grid, cqgs, cqgs_glob)
        CALL PACK_COLUMN(grid, ipbl, ipbl_glob)

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
        if ( AM_I_ROOT() ) then
          READ (KUNIT,ERR=10) HEADER,UABL_glob,VABL_glob,TABL_glob,
     &         QABL_glob,EABL_glob,CMGS_glob,CHGS_glob,CQGS_glob,
     &         IPBL_glob
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if

        call UNPACK_BLOCK(grid, uabl_glob, uabl)
        call UNPACK_BLOCK(grid, vabl_glob, vabl)
        call UNPACK_BLOCK(grid, tabl_glob, tabl)
        call UNPACK_BLOCK(grid, qabl_glob, qabl)
        call UNPACK_BLOCK(grid, eabl_glob, eabl)

        call UNPACK_COLUMN(grid, cmgs_glob, cmgs)
        call UNPACK_COLUMN(grid, chgs_glob, chgs)
        call UNPACK_COLUMN(grid, cqgs_glob, cqgs)
        call UNPACK_COLUMN(grid, ipbl_glob, ipbl)
#ifdef TRACERS_ON
        SELECT CASE (IACTION)
        CASE (IOREAD,IRERUN,IRSFIC,IRSFICNO)    ! restarts
          if ( AM_I_ROOT() ) then
            READ (KUNIT,ERR=10) TR_HEADER,TRABL_GLOB
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in tracer module version ",TR_HEADER
     *             ,TR_MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_BLOCK(grid, TRABL_GLOB, TRABL)
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
      USE DOMAIN_DECOMP, only : GET, grid, AM_I_ROOT
      USE DOMAIN_DECOMP, only : CHECKSUM, CHECKSUM_COLUMN
      USE DOMAIN_DECOMP, only : UNPACK_DATA, UNPACK_COLUMN
      USE DOMAIN_DECOMP, only : PACK_DATA, PACK_COLUMN
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "BLD02"
      INTEGER :: J_0, J_1
      REAL*8, DIMENSION(IM,JM) :: wsavg_glob,tsavg_glob,
     *       qsavg_glob,dclev_glob,usavg_glob,
     *       vsavg_glob,tauavg_glob,qgavg_glob,tgvavg_glob
      REAL*8, DIMENSION(LM,IM,JM) :: egcm_glob, w2gcm_glob
      REAL*8, DIMENSION(4,IM,JM) :: ustar_pbl_glob

      MODULE_HEADER(lhead+1:80) = 'R8 dim(ijm):ws,ts,qs,'//
     *  'LvlDC,us,vs,tau,u*(4,.),ke;w2(lijm),tgv,qg'

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_DATA(grid, wsavg, wsavg_glob)
        CALL PACK_DATA(grid, tsavg, tsavg_glob)
        CALL PACK_DATA(grid, qsavg, qsavg_glob)
        CALL PACK_DATA(grid, dclev, dclev_glob)
        CALL PACK_DATA(grid, usavg, usavg_glob)
        CALL PACK_DATA(grid, vsavg, vsavg_glob)
        CALL PACK_DATA(grid, tauavg, tauavg_glob)
        CALL PACK_DATA(grid, tgvavg, tgvavg_glob)
        CALL PACK_DATA(grid, qgavg, qgavg_glob)

        CALL PACK_COLUMN(grid, ustar_pbl, ustar_pbl_glob)
        CALL PACK_COLUMN(grid, egcm, egcm_glob)
        CALL PACK_COLUMN(grid, w2gcm, w2gcm_glob)

        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,wsavg_glob,tsavg_glob
     *         ,qsavg_glob,dclev_glob,usavg_glob,vsavg_glob,tauavg_glob
     *         ,ustar_pbl_glob,egcm_glob,w2gcm_glob,tgvavg_glob
     *         ,qgavg_glob
        END IF

      CASE (IOREAD:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER,wsavg_glob,tsavg_glob,
     *       qsavg_glob,dclev_glob,usavg_glob,
     *       vsavg_glob,tauavg_glob,ustar_pbl_glob,egcm_glob,
     *       w2gcm_glob,tgvavg_glob,qgavg_glob
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if

        call UNPACK_DATA(grid, wsavg_glob, wsavg)
        call UNPACK_DATA(grid, tsavg_glob, tsavg)
        call UNPACK_DATA(grid, qsavg_glob, qsavg)
        call UNPACK_DATA(grid, dclev_glob, dclev)
        call UNPACK_DATA(grid, usavg_glob, usavg)
        call UNPACK_DATA(grid, vsavg_glob, vsavg)
        call UNPACK_DATA(grid, tauavg_glob, tauavg)
        call UNPACK_DATA(grid, tgvavg_glob, tgvavg)
        call UNPACK_DATA(grid, qgavg_glob, qgavg)

        call UNPACK_COLUMN(grid, ustar_pbl_glob, ustar_pbl)
        call UNPACK_COLUMN(grid, egcm_glob, egcm)
        call UNPACK_COLUMN(grid, w2gcm_glob, w2gcm)

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

      ALLOCATE(    uabl(npbl,4,    im,J_0H:J_1H),
     *             vabl(npbl,4,    im,J_0H:J_1H),
     *             tabl(npbl,4,    im,J_0H:J_1H),
     *             qabl(npbl,4,    im,J_0H:J_1H),
     *             eabl(npbl,4,    im,J_0H:J_1H),
     *         STAT=IER)
      qabl=0.  ! initialise to make life easier

#ifdef TRACERS_ON
      ALLOCATE(    trabl(npbl,ntm,4,im,J_0H:J_1H),
     *         STAT=IER)
#endif

      ALLOCATE(    cmgs(4,im,J_0H:J_1H),
     *             chgs(4,im,J_0H:J_1H),
     *             cqgs(4,im,J_0H:J_1H),
     *             ipbl(4,im,J_0H:J_1H),
     *         STAT=IER)
      ipbl(:,:,J_0H:J_1H) = 0

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

      ALLOCATE(    ustar_pbl(4,im,J_0H:J_1H),
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

      ALLOCATE( gustiwind(im,J_0H:J_1H),STAT=IER)
      gustiwind(:,J_0H:J_1H) = 0.

      END SUBROUTINE ALLOC_PBL_COM

