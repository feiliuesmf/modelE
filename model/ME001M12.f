C**** ME001M12 E001M12 SOMTQ M_f90 MB399M12   07/00
C****
C**** Second order scheme for momentum advection, with FLTRUV
C****
C**** Snow ages each day independent of temperature
C****
C**** No extra leads - works with both: KOCEAN=0 and KOCEAN=1
C****
C**** Basic model II (OA,PALMER omitted) .5 box longitude shift
C**** Pressure replaces Sigma above LS1 as the vertical coordinate
C**** Modified for using new MC codes, radiation, 11 veg.types
C**** Quadratic upstream scheme + 4th order scheme, Cor.term=0 at poles
C**** Routines included: MAIN,INPUT,DAILY,CHECKT,CHECK3
C**** f90 changes
*****
      USE CONSTANT, only : sday
      USE E001M12_COM
      USE GEOM
      USE RANDOM

      USE DAGCOM, only : aj,kacc,aij,tsfrez,keynr,kdiag,ij_tgo2
      USE DYNAMICS, only : FILTER,CALC_AMPK
      USE OCEAN, only : ODATA,OA

      IMPLICIT NONE

      INTEGER I,J,L,K,LLAB1,NDAILY,MSRCE,KSS6,NSTEP0,I24,MNOW,MINC
     *     ,MSUM,MSTART,ITAU,IDTHR,INTFX,MODD5D,iowrite,ioerr
      REAL*8 DTIME,PELSE,PDIAG,PSURF,PRAD,PCDNS,PDYN,TOTALT
     *     ,RUNON,PERCNT,DTHR,XTAU
!@param HR12,HR24 useful constants
      REAL*8, PARAMETER :: HR12=12., HR24=24.

      CHARACTER CYEAR*4,CMMND*80
      CHARACTER*8 :: LABSSW,OFFSSW = 'XXXXXXXX'
      CHARACTER*80 TITLE_T
      LOGICAL EVENT
C**** STATEMENT FUNCTION CONVERTING HOURS TO INTERNAL TIME UNITS
      INTFX(XTAU)=INT(XTAU*XINT+.5)
C**** STATEMENT FUNCTION EVENT IS TRUE IF TAU IS A MULTIPLE OF XTAU
      EVENT(XTAU)=MOD(ITAU,INTFX(XTAU)).LT.IDTHR

      iowrite=-1
      CALL TIMER (MNOW,MINC,MSUM)
      CALL INPUT
      WRITE (3) OFFSSW
      CLOSE (3)
      MSTART=MNOW-MDYN-MCNDS-MRAD-MSURF-MDIAG-MELSE
C**** INITIALIZE TIME PARAMETERS
      DTHR=DT/3600.
      IDTHR=INTFX(DTHR)
      I24=INTFX(HR24)
      NSTEP0=.5+TAUI/DTHR
      NSTEP=INT(.5+TAU/DTHR)-NSTEP0
      ITAU=(NSTEP+NSTEP0)*IDTHR
      TAU=DFLOAT(ITAU)/XINT
      IDAY=1+ITAU/I24
      TOFDAY=(ITAU-(IDAY-1)*I24)/XINT
         MODD5K=1000
      CALL DAILY0
      CALL CALC_AMPK(LM)
         CALL CHECKT ('INPUT ')

      CALL TIMER (MNOW,MINC,MELSE)
      PERCNT=100.*MELSE/(MNOW-MSTART+.001)
      WRITE (6,'(A,13X,A,I6,A,F6.2,I6,A5,I27,I7,F7.1,A,F11.2)')
     *  '0CLIMATE MODEL STARTED UP','DAY',IDAY,', HR',TOFDAY,
     *   JDATE,JMONTH,MINC,MELSE,PERCNT,' TAU',TAU

      RUNON=1.
      IF (TAU.GE.TAUE) RUNON=-1.
C****
C**** MAIN LOOP
C****
      DO WHILE (TAU.LT.TAUE)

C**** EVERY TAUT HOURS (OR IF IT IS FIRST TIME-STEP)
C**** WRITE RESTART INFORMATION ONTO DISK
      IF (EVENT(TAUT).or.TAU.LE.TAUI+DTHR*(NDYN+.5)) THEN
         CALL RFINAL (IRAND)
         call io_rsf(KDISK,TAU,iowrite,ioerr)
         CALL TIMER (MNOW,MINC,MELSE)
         PERCNT=100.*MELSE/(MNOW-MSTART+1.D-5)
         WRITE (6,'(A,I3,55X,2I7,F7.1,A,F11.2)')
     *        ' OUTPUT RECORD WRITTEN ON UNIT',KDISK,MINC,MELSE,PERCNT
     *        ,' TAU',TAU
         KDISK=3-KDISK
      END IF

C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY DAY
      IF (EVENT(HR24)) THEN
C**** CHECK FOR BEGINNING OF EACH MONTH => RESET DIAGNOSTICS
         DO K=1,13
            IF (JDAY.EQ.NDZERO(K)) call reset_DIAG
         END DO
C**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
         call daily_DIAG
      END IF
C****
C**** INTEGRATE DYNAMIC TERMS (DIAGA AND DIAGB ARE CALLED FROM DYNAM)
C****
         MODD5D=MOD(NSTEP,NDA5D)
         IF (MODD5D.EQ.0) CALL DIAG5A (2,0)
         IF (MODD5D.EQ.0) CALL DIAG9A (1)
      CALL DYNAM

      CALL CALC_AMPK(LS1-1)

      CALL CHECKT ('DYNAM ')
      CALL TIMER (MNOW,MINC,MDYN)
      PERCNT=100.*MDYN/(MNOW-MSTART)

         IF (MODD5D.EQ.0) CALL DIAG5A (7,NDYN)
         IF (MODD5D.EQ.0) CALL DIAG9A (2)
         IF (EVENT(HR12)) CALL DIAG7A
C****
C**** INTEGRATE SOURCE TERMS
C****
         MODRD=MOD(NSTEP,NRAD)
         MODD5S=MOD(NSTEP,NDA5S)
         IF (MODD5S.EQ.0) IDACC(8)=IDACC(8)+1
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG5A (1,0)
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG9A (1)
C**** CONDENSATION, SUPER SATURATION AND MOIST CONVECTION
c      CALL MSTCNV
c      CALL CONDSE
      CALL MC_COND
      CALL CHECKT ('CONDSE ')
      CALL TIMER (MNOW,MINC,MCNDS)
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NCNDS)
         IF (MODD5S.EQ.0) CALL DIAG9A (3)
C**** RADIATION, SOLAR AND THERMAL
      CALL RADIA
      CALL CHECKT ('RADIA ')
      CALL TIMER (MNOW,MINC,MRAD)
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NCNDS)
         IF (MODD5S.EQ.0) CALL DIAG9A (4)
C**** SURFACE INTERACTION AND GROUND CALCULATION
      CALL PRECIP
      CALL CHECKT ('PRECIP')
      CALL SURFCE
      CALL CHECKT ('SURFCE')
      CALL GROUND
      CALL CHECKT ('GROUND')
      CALL DRYCNV
      CALL CHECKT ('DRYCNV')
      CALL TIMER (MNOW,MINC,MSURF)
         IF (MODD5S.EQ.0) CALL DIAG9A (5)
C**** STRATOSPHERIC MOMENTUM DRAG
      CALL SDRAG
      CALL CHECKT ('SDRAG ')
      CALL TIMER (MNOW,MINC,MSURF)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NCNDS)
         IF (MODD5S.EQ.0) CALL DIAG9A (6)
      MSRCE=MCNDS+MRAD+MSURF
      PERCNT=100.*MSRCE/(MNOW-MSTART)
C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(NSTEP,NFILTR).EQ.0) THEN
         IDACC(10)=IDACC(10)+1
         IF (MODD5S.NE.0) CALL DIAG5A (1,0)
         CALL DIAG9A (1)
      CALL FILTER
      CALL CHECKT ('FILTER')
      CALL TIMER (MNOW,MINC,MDYN)
         CALL DIAG5A (14,NFILTR)
         CALL DIAG9A (7)
      END IF
C****
C**** UPDATE MODEL TIME AND CALL DAILY IF REQUIRED
C****
      NSTEP=NSTEP+NDYN
      ITAU=(NSTEP+NSTEP0)*IDTHR
      TAU=DFLOAT(ITAU)/XINT
      IDAY=1+ITAU/I24
      TOFDAY=(ITAU-(IDAY-1)*I24)/XINT

      IF (EVENT(HR24)) THEN
         CALL DIAG5A (1,0)
         CALL DIAG9A (1)
      CALL DAILY
      CALL TIMER (MNOW,MINC,MELSE)
         NDAILY=SDAY/DT
         CALL DIAG5A (16,NDAILY)
         CALL DIAG9A (8)
      call daily_SNOW   
      call daily_OCEAN
         CALL CHECKT ('DAILY ')
      CALL TIMER (MNOW,MINC,MSURF)
      END IF
C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY USET HOURS
C****
      IF (USET.GT.0.) THEN
         IF (EVENT(USET)) THEN
            WRITE (20) TAU,OA
            ENDFILE 20
            BACKSPACE 20
            WRITE (6,'(A,78X,A,F11.2)')
     *           ' INFORMATION WRITTEN ON UNIT 20',' TAU',TAU
C**** ZERO OUT INTEGRATED QUANTITIES
            OA(:,:,4:12)=0.
         ELSEIF (EVENT(USET/2.)) THEN
            call uset_OCEAN
         END IF
         CALL TIMER (MNOW,MINC,MELSE)
         PERCNT=100.*MELSE/(MNOW-MSTART)
      END IF
C****
C**** CALL DIAGNOSTIC ROUTINES
C****
         IF (MOD(NSTEP+NDA4,NDA4).EQ.0) CALL DIAG4A
         IF (USESLP.NE.0.) THEN
         IF (USESLP.LT.0..AND.EVENT(-USESLP)) USESLP=-USESLP
         IF (USESLP.GT.0..AND.EVENT( USESLP)) CALL DIAG10(0)
         END IF
         IF (NDPRNT(1).LT.0) THEN
C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
         IF (KDIAG(1).LT.9) CALL DIAGJ
         IF (KDIAG(2).LT.9) CALL DIAGJK
         IF (KDIAG(2).LT.9) CALL DIAGJL
         IF (KDIAG(10).EQ.0) CALL DIAGIL
         IF (KDIAG(7).LT.9) CALL DIAG7P
         IF (KDIAG(3).LT.9) CALL DIAGIJ
         IF (KDIAG(9).LT.9) CALL DIAG9P
         IF (KDIAG(5).LT.9) CALL DIAG5P
         IF (KDIAG(4).LT.9) CALL DIAG4
         IF (TAU.LE.TAUI+DTHR*(NDYN+.5)) THEN
            CALL DIAGKN
         ELSE ! RESET THE UNUSED KEYNUMBERS TO ZERO
            KEYNR(1:42,KEYCT)=0
         END IF
         NDPRNT(1)=NDPRNT(1)+1
         END IF

C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES ON NDPRNT-TH DAY OF RUN
      IF (EVENT(HR24)) THEN
      DO K=1,13
      IF (JDAY.EQ.NDPRNT(K)) THEN
         WRITE (6,'("1"/64(1X/))')
         IF (KDIAG(1).LT.9) CALL DIAGJ
         IF (KDIAG(2).LT.9) CALL DIAGJK
         IF (KDIAG(2).LT.9) CALL DIAGJL
         IF (KDIAG(10).EQ.0) CALL DIAGIL
         IF (KDIAG(7).LT.9) CALL DIAG7P
         IF (KDIAG(3).LT.9) CALL DIAGIJ
         IF (KDIAG(9).LT.9) CALL DIAG9P
         IF (KDIAG(5).LT.9) CALL DIAG5P
         IF (KDIAG(6).LT.9) CALL DIAG6
         IF (KDIAG(4).LT.9) CALL DIAG4
      END IF
      END DO

C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
C****   (NORMALLY DONE AT THE END OF A MONTH)
      DO K=1,13
      IF (JDAY.EQ.NDZERO(K)) THEN
C**** PRINT THE KEY DIAGNOSTICS
         CALL DIAGKN
      IF (KCOPY.GT.0) THEN
C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
         WRITE (CYEAR,'(I4)') JYEAR0
         LLAB1 = INDEX(LABEL1,'(') -1
         IF(LLAB1.LT.1) LLAB1=16
         IF (KCOPY.GT.1) THEN
C**** KCOPY > 1 : SAVE THE RESTART INFORMATION
            CALL RFINAL (IRAND)
            OPEN(30,FILE=JMNTH0(1:3)//CYEAR//'.rsf'//LABEL1(1:LLAB1),
     *           FORM='UNFORMATTED')

            call io_rsf(30,TAU,iowrite,ioerr)
            
            CLOSE (30)
         END IF
         IF (KCOPY.GT.2) THEN
C**** KCOPY > 2 : SAVE THE OCEAN DATA FOR INITIALIZING DEEP OCEAN RUNS
            OPEN (30,FILE=JMNTH0(1:3)//CYEAR//'.oda'//LABEL1(1:LLAB1),
     *           FORM='UNFORMATTED')
            WRITE (30) TAU,ODATA,((AIJ(I,J,IJ_TGO2),I=1,IM),J=1,JM)
            CLOSE (30)
         END IF
C**** SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
         OPEN (30,FILE=JMNTH0(1:3)//CYEAR//'.acc'//LABEL1(1:LLAB1),
     *        FORM='UNFORMATTED')
         WRITE (30) SNGL(TAU),JC,CLABEL,(SNGL(RC(I)),I=1,161),KEYNR,
     *        (SNGL(TSFREZ(I,1,1)),I=1,IM*JM*2),
     *        (SNGL(AJ(I,1)),I=1,KACC),SNGL(TAU)
         CLOSE (30)
      END IF

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
      CALL TIMER (MNOW,MINC,MDIAG)
      TOTALT=.01*(MNOW-MSTART)
      PDYN  = MDYN/TOTALT
      PCDNS = MCNDS/TOTALT
      PRAD  = MRAD/TOTALT
      PSURF = MSURF/TOTALT
      PDIAG = MDIAG/TOTALT
      PELSE = MELSE/TOTALT
      DTIME = 24.*TOTALT/(60.*(TAU-TAU0))
      WRITE (6,'(/A,F7.2,A,F5.1,A,F5.1,A,F5.1,A,F5.1,A,F5.1,A,F5.1//)')
     *  '0TIME',DTIME,'(MINUTES)    DYNAMICS',PDYN,
     *  '    CONDENSATION',PCDNS,'    RADIATION',PRAD,'    SURFACE',
     *  PSURF,'    DIAGNOSTICS',PDIAG,'    OTHER',PELSE
      MDYN  = 0
      MCNDS = 0
      MRAD  = 0
      MSURF = 0
      MDIAG = 0
      MELSE = 0
      MSTART= MNOW
      END IF
      END DO
      END IF
C**** TIME FOR CALLING DIAGNOSTICS
      CALL TIMER (MNOW,MINC,MDIAG)

C**** TEST FOR TERMINATION OF RUN
      READ (3,END=210) LABSSW
 210  CLOSE (3)
      IF (LABSSW.EQ.LABEL1(1:8)) THEN
C**** RUN TERMINATED BECAUSE SENSE SWITCH 6 WAS TURNED ON
         KSS6=1
         WRITE (6,'("0SENSE SWITCH 6 HAS BEEN TURNED ON.")')
         EXIT
      END IF

      END DO
C****
C**** END OF MAIN LOOP
C****

C**** PRINT OUT DIAGNOSTICS IF NO TIME-STEP WAS DONE 
      IF (RUNON.eq.-1) THEN
         WRITE (6,'("1"/64(1X/))')
         IF (KDIAG(1).LT.9) CALL DIAGJ
         IF (KDIAG(2).LT.9) CALL DIAGJK
         IF (KDIAG(2).LT.9) CALL DIAGJL
         IF (KDIAG(10).EQ.0) CALL DIAGIL
         IF (KDIAG(7).LT.9) CALL DIAG7P
         IF (KDIAG(3).LT.9) CALL DIAGIJ
         IF (KDIAG(9).LT.9) CALL DIAG9P
         IF (KDIAG(5).LT.9) CALL DIAG5P
         IF (KDIAG(6).LT.9) CALL DIAG6
         IF (KDIAG(4).LT.9) CALL DIAG4
C**** PRINT THE KEY DIAGNOSTICS IF END OF MONTH
         DO K=1,13
            IF (JDAY.EQ.NDZERO(K)) CALL DIAGKN
         END DO
      END IF

C**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      call io_rsf(KDISK,TAU,iowrite,ioerr)
      WRITE (6,'(A,I3,77X,A,F11.2)')
     *  ' OUTPUT RECORD WRITTEN ON UNIT',KDISK,'TAU',TAU

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)
      WRITE (6,'(/////4(1X,33("****")/)//,A,F11.2,I6,F7.2
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY.TAU,IDAY,TOFDAY=',TAU,IDAY,TOFDAY
      IF (KSS6.EQ.1) STOP 12
      STOP 13
      END

      BLOCK DATA BDINP
C****
C**** DEFAULT PARAMETERS FOR MODEL COMMON BLOCK
C****
      USE E001M12_COM
     &     , only : im,jm,lm

      IMPLICIT NONE 

      INTEGER, DIMENSION(3) :: IDUM
      INTEGER, DIMENSION(13) :: NDZERO,NDPRNT
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC
      INTEGER :: IM0,JM0,LM0,JMM1x,LMM1,LS1,LTMx,LBLMx,LMCMx,LSSMx,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRAx,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0
      COMMON /IPARMB/IM0,JM0,LM0,JMM1x,LMM1,LS1,LTMx,LBLMx,LMCMx,LSSMx,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRAx,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0,IDUM   ,NDZERO    ,NDPRNT    ,  IJD6     ,IDACC

      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(60+2*(36-LM)) :: RDM2
      DOUBLE PRECISION ::  TAU,TAU0,TOFDAY,TOFDY0,DT,TAUP,TAUI,TAUE,TAUT
     *     ,TAUO,PTOP,PSF,PSDRAG,PTRUNC,XINT,
     *     SKIPSE,USESLP,USEP,USET,RSDIST,SIND,COSD,PSFMPT,PSTRAT
       COMMON /RPARMB/
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  PTOP,PSF,PSDRAG,PTRUNC,
     *         XINT,                 SKIPSE,USESLP,USEP,USET,
     *  RSDIST,SIND,COSD,SIG,SIGE,   TAUTR0,PSFMPT,PSTRAT,RDM2


      CHARACTER*4 NAMD6,JMONTH,JMNTH0
      CHARACTER*132 XLABEL
      COMMON /TEXT/ XLABEL,NAMD6(4),JMONTH,JMNTH0

      DATA IM0,JM0,LM0,  LS1/
     *     IM ,JM ,LM ,   12/,
     *  KOCEAN,KDISK,KEYCT,KCOPY,     IRAND,MFILTR,NDYN/
     *       1,    1,    1,    2, 123456789,     1,   8/,
     *       NSURF,NGRND,  IYEAR/
     *           2,    1,   1976/,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,  MDIAG,MELSE,MODRD,MODD5K,MODD5S/
     *      0,   0,    0,   0,    0,      0,    0,    0,     0,     0/
      DATA  DT,  TAUP,TAUI,TAUE,TAUT/
     *    450.,   -1., -1.,  1., 24./,
     *  PTOP, PSF, PSDRAG,PTRUNC/
     *  150.,984.,  500.,    0./,
     *                      XINT, SKIPSE,USESLP,USEP,USET/
     *                      120.,     0.,    0.,  0.,  0./
      DATA SIGE /1.0000000,LM*0./
      DATA NAMD6 /'AUSD','MWST','SAHL','EPAC'/,
     *  NDZERO/ 0,1,32,60,91,121,152,182,213,244,274,305,335/,
     *  NDPRNT/-1,1,32,60,91,121,152,182,213,244,274,305,335/,
     *  IJD6/63,17, 17,34, 37,27, 13,23/

      LOGICAL Q_GISS,Q_HDF,Q_PRT,Q_NETCDF
      COMMON /Q_PP/Q_GISS,Q_HDF,Q_PRT,Q_NETCDF
      DATA Q_GISS,Q_HDF,Q_PRT,Q_NETCDF/.FALSE.,.FALSE.,.FALSE.,.FALSE./
      END

      SUBROUTINE INPUT
C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,rhow
      USE E001M12_COM
      USE SOMTQ_COM
      USE GEOM
      USE GHYCOM
     &  , only : ghdata
      USE RANDOM
      USE RADNCB, only : RQT,SRHR,TRHR,FSF,S0X,CO2
      USE CLD01_COM_E001, only : TTOLD,QTOLD,SVLHX,RHSAV,CLDSAV,
     *     U00wtr,U00ice,LMCM
      USE PBLCOM
     &     , only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq=>cqgs
     &  ,ipbl,bldata,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
      USE DAGCOM, only : aj,kacc,tsfrez,tdiurn,kdiag,keynr,jreg
     &  ,TITREG,NAMREG,iwrite,jwrite,itwrite,qcheck
      USE OCEAN, only : odata,OA,T50
      USE FILEMANAGER, only : getunit

      IMPLICIT NONE
!@var iu_AIC,iu_TOPO,iu_GIC,iu_REG unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG

      INTEGER I,J,L,K,KLAST,IUNIT2,IUNIT1,IUNIT,KDISK0,ITYPE,IM1,KTACC
     *     ,IR,IREC,NOFF,NDIFS,ISTART,ioread,ioerr
      REAL*8 TIJL,X,TAU2,TAUX,TAUY,TAU1,RDSIG,CDM,SNOAGE,TEMP

      REAL*8 EXPBYK
      INTEGER JC1(100)
      REAL*8 RC1(161)
      CHARACTER*4 RUNID
      CHARACTER*156 CLABEL1

      COMMON/WORK1/NLREC(256),SNOAGE(IM,JM,2)

      LOGICAL :: redoGH = .FALSE.,iniPBL = .FALSE.

      CHARACTER NLREC*80
      REAL*4 TAU4,TAUY4,XX4
      NAMELIST/INPUTZ/IM0,JM0,LM0,LS1,KOCEAN,ISTART,
     *  KDISK,TAUP,TAUI,TAUE,TAUT,TAUO,NDYN,NCNDS,NRAD,NSURF,NGRND,
     *  NFILTR,NDAA,NDA5D,NDA5K,NDA5S,NDA4,NDASF,DT,TAU,XINT,IYEAR,
     *     U00wtr,U00ice,S0X,CO2,
     *     PTOP,PSF,PTRUNC,IRAND,MFILTR
     *     ,NDIFS,KACC0,KEYCT,SKIPSE,USESLP,USEP,USET,KCOPY
     *     ,KDIAG,NDZERO,NDPRNT,IJD6,NAMD6,SIG,SIGE,KTACC0
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK
     *     ,Q_GISS,Q_HDF,Q_PRT,Q_NETCDF

      ioread=1
      ISTART=10
C**** READ SPECIAL REGIONS FROM UNIT 29 
      call getunit("REG",iu_REG,.TRUE.)
      READ(iu_REG) TITREG,JREG,NAMREG
      WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG

      KDIAG(1:12)=0
C**** INITIALIZE ADVECTION TERMS FOR SECOND ORDER MOMENTS
      DO I=1,IM*JM*LM*9
         TX(I,1,1) = 0.
      END DO
      DO I=1,IM*JM*LM*9
         QX(I,1,1) = 0.
      END DO
      WRITE (6,'(A,40X,A/)') '0','GISS N LAYER WEATHER MODEL'
C**** Read in Label and Namelist parameters from rundeck
      READ(5,'(A80)') NLREC(1),NLREC(2)
      NOFF=0
      IF (NLREC(1)(73:80).EQ.'        ') NOFF=8
      XLABEL(1:80-NOFF)=NLREC(1)(1:80-NOFF)
      XLABEL(81-NOFF:128)=NLREC(2)(1:48+NOFF)
      RUNID=XLABEL(1:4)
      XLABEL(129:132) = ' '
      WRITE (6,'(A,A/)') '0',XLABEL
C**** COPY INPUTZ NAMELIST ONTO CORE TAPE AND TITLE PAGE
      IREC=0
   60 IREC=IREC+1
      READ  (5,'    (A80)') NLREC(IREC)
      WRITE (6,'(35X,A80)') NLREC(IREC)
      IF (NLREC(IREC)(1:5).NE.' &END')  GO TO 60
CSGI  READ (NLREC,INPUTZ)
      REWIND 8
      WRITE(8,'(A)') (NLREC(IR),IR=1,IREC)
      REWIND 8
      READ (8,NML=INPUTZ)
      REWIND 8

      IF (ISTART.GE.10) GO TO 90
C**** get unit for atmospheric initial conditions
      call getunit("AIC",iu_AIC,.TRUE.)

C**** SET STRICTLY DEPENDENT QUANTITIES
      LMM1=LM-1
      DO L=1,LM
         IF(SIGE(L).EQ.0.) LS1=L
         SIG(L)=.5*(SIGE(L)+SIGE(L+1))
      END DO
      PSFMPT = PSF-PTOP
      PSTRAT = (PSF-PTOP)*(SIGE(LS1)-SIGE(LM+1))

      TAUO=TAUP
      IF (TAUP.LT.0.) TAUP=TAUI
      NCNDS=NDYN
      NRAD=5*NDYN
      NFILTR=2*NDYN
         NDAA=7*NDYN+2
         NDA5D=7*NDYN
         NDA5K=NDAA
         NDA5S=7*NDYN
         NDA4=24*NDYN
         NDASF=2*NSURF-1
         KACC0=KACC
         KTACC0 = KTACC
   90 GO TO (100,200,300,320,322,325,890,890,890,400,440,440,400),ISTART
      GO TO 400
C****
C**** FUNCTIONAL DETERMINATION OF PROGNOSTIC QUANTITIES, ISTART=1
C****
  100 TEMP=250.
      IF(TAUI.LT.0.) TAUI=0.
      P(:,:)=PSF-PTOP
      TSAVG(:,:)=TEMP
      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=TEMP
      Q(:,:,:)=3.D-6
      GO TO 210
CALT  DEFINE GDATA(1->14),GHDATA(1->4*NGM+5) AND GO TO 220
C****
C**** INITIALIZE A RUN FROM ATMOSPHERIC CONDITIONS AND
C**** GROUND CONDITIONS, ISTART=2
C****
  200 CALL READT (iu_AIC,0,P,IM*JM,P,1)             ! Psurf
      DO I=1,IM*JM
      P(I,1)=P(I,1)-PTOP                               ! Psurf -> P
      END DO
      DO L=1,LM
      CALL READT (iu_AIC,0,U(1,1,L),IM*JM,U(1,1,L),1) ! U
      END DO
      DO L=1,LM
      CALL READT (iu_AIC,0,V(1,1,L),IM*JM,V(1,1,L),1) ! V
      END DO
      DO L=1,LM
      CALL READT (iu_AIC,0,T(1,1,L),IM*JM,T(1,1,L),1) ! Temperature
      END DO
      DO L=1,LM
      CALL READT (iu_AIC,0,Q(1,1,L),IM*JM,Q(1,1,L),1) ! Q
      END DO
      CALL READT (iu_AIC,0,TSAVG(1,1),IM*JM,TSAVG(1,1),1)  ! Tsurf
      CLOSE (iu_AIC)
C**** new: GDATA(8) UNUSED,GDATA(9-11) SNOW AGE OVER OCN.ICE,L.ICE,EARTH
  210 call getunit("GIC",iu_GIC,.TRUE.)
      READ(iu_GIC,ERR=830) GDATA,GHDATA,(ODATA(I,1,1),I=1,IM*JM*2)
      CLOSE (iu_GIC)
C**** Check whether a proper TAUI is given - initialize TAU=model time
  220 IF(TAUI.LT.0.) THEN
         WRITE(6,*) 'PLEASE SET TAUI IN THE RUNDECK'
         STOP 'ERROR: TAUI UNDEFINED'
      END IF
      TAU=TAUI
      TAUX=TAUI
      WSAVG(1,1)=SQRT(U(1,2,1)*U(1,2,1)+V(1,2,1)*V(1,2,1))
      USAVG(1,1)=U(1,2,1)
      VSAVG(1,1)=V(1,2,1)
      WSAVG(1,JM)=SQRT(U(1,JM,1)*U(1,JM,1)+V(1,JM,1)*V(1,JM,1))
      USAVG(1,JM)=U(1,JM,1)
      VSAVG(1,JM)=V(1,JM,1)
      DO J=2,JM-1
         IM1=IM
         DO I=1,IM
            WSAVG(I,J)=.25*SQRT(
     *           (U(IM1,J,1)+U(I,J,1)+U(IM1,J+1,1)+U(I,J+1,1))**2
     *           +(V(IM1,J,1)+V(I,J,1)+V(IM1,J+1,1)+V(I,J+1,1))**2)
            USAVG(I,J)=.25*(U(IM1,J,1)+U(I,J,1)+U(IM1,J+1,1)+U(I,J+1,1))
            VSAVG(I,J)=.25*(V(IM1,J,1)+V(I,J,1)+V(IM1,J+1,1)+V(I,J+1,1))
            IM1=I
         END DO
      END DO
      CDM=.001
      DO J=1,JM
      DO I=1,IM
C**** SET SURFACE MOMENTUM TRANSFER TAU0
      TAUAVG(I,J)=CDM*WSAVG(I,J)**2
C**** SET LAYER THROUGH WHICH DRY CONVECTION MIXES TO 1
      DCLEV(I,J)=1.
C**** SET SURFACE SPECIFIC HUMIDITY FROM FIRST LAYER HUMIDITY
      QSAVG(I,J)=Q(I,J,1)
C**** SET RADIATION EQUILIBRIUM TEMPERATURES FROM LAYER LM TEMPERATURE
      DO K=1,3
         RQT(I,J,K)=T(I,J,LM)
      END DO
C**** REPLACE TEMPERATURE BY POTENTIAL TEMPERATURE
      DO L=1,LS1-1
         RHSAV(I,J,L)=.85
         CLDSAV(I,J,L)=0.
         SVLHX(I,J,L)=0.
         T(I,J,L)=T(I,J,L)/EXPBYK(SIG(L)*P(I,J)+PTOP)
      END DO
      DO L=LS1,LM
         RHSAV(I,J,L)=.85
         CLDSAV(I,J,L)=0.
         SVLHX(I,J,L)=0.
         T(I,J,L)=T(I,J,L)/((SIG(L)*(PSF-PTOP)+PTOP)**KAPA)
      END DO
      DO L=1,LM
         TTOLD(I,J,L)=T(I,J,L)
         QTOLD(I,J,L)=Q(I,J,L)
         WM(I,J,L)=0.
      END DO
      IF (LS1.LE.LM) THEN
C**** SET STRATOSPHERIC SPECIFIC HUMIDITY TO 3.D-6
         DO L=LS1,LM
            Q(I,J,L)=3.D-6
         END DO
      END IF
      END DO
      END DO
C**** INITIALIZE TSFREZ
      DO J=1,JM
         DO I=1,IM
            TSFREZ(I,J,1)=365.
            TSFREZ(I,J,2)=365.
         END DO
      END DO
      DO ITYPE=1,4
         DO J=1,JM
            DO I=1,IM
               USTAR(I,J,ITYPE)=WSAVG(I,J)*SQRT(CDM)
            END DO
         END DO
      END DO
CALT  GO TO 327       ! possibility to make tracer slopes more realistic
      GO TO 350
C****
C**** INITIALIZE RUN FROM PREVIOUS MODEL OUTPUT, ISTART=3-8
C**** ISTART=3-4  IC-file looks like restart file (same # of prog.var)
C****     ISTART=3: C ARRAY IS BUILT UP FROM DEFAULTS AND NAMELIST
 300  READ (iu_AIC,ERR=800,END=810) TAUX,JC1,CLABEL1,RC1,KEYNR,
     *     U,V,T,P,Q,
     2  ODATA,GDATA,GHDATA,BLDATA,
     *   uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     3  TTOLD,QTOLD,SVLHX,RHSAV,WM,CLDSAV,
     4  TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     5  QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ,
     6  RQT,T50
      CLOSE (iu_AIC)
      GO TO 398
C****     ISTART=4: C ARRAY IS COPIED FROM INPUT DATA EXCEPT FOR XLABEL
C**** USE ISTART=4 TO START FROM THIS RUN'S .rsf FILE
C****                                                               ****
C**** Use ISTART=4 for reruns and extensions, NOT to start new runs ****
C***********************************************************************
 320  READ(iu_AIC,ERR=800)TAUX,JC,CLABEL1,RC,KEYNR,U,V,T,P,Q,ODATA,
     *  GDATA,GHDATA,BLDATA,
     *  uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     2  TTOLD,QTOLD,SVLHX,RHSAV,WM,CLDSAV,
     *  TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     *  QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ,
     *  T50,RQT,SRHR,TRHR,FSF,TSFREZ
      CLABEL(133:156)=CLABEL1(133:156)
      JC(32:37)=0
      CLOSE (iu_AIC)
      GO TO 399
C**** ISTART=5-9; ICfile looks different from restart file
C**** ISTART=5; initial start from run B140 (snow ages not yet in GDATA)
 322  READ (iu_AIC,ERR=800,END=810) TAUX,JC1,CLABEL1,RC1,KEYNR,
     2  U,V,T,P,Q,ODATA,GDATA,GHDATA,BLDATA,
     3  TTOLD,QTOLD,SVLHX,RHSAV,WM,CLDSAV,
     4  TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     5  QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ,
     6  RQT,SRHR,TRHR,TSFREZ,SNOAGE
      close (iu_AIC)
      redoGH=.TRUE.
      GO TO 350
C**** ISTART=6 ; start from run B120 (no QUS, only 1 snow age)
 325  READ(iu_AIC,ERR=800)TAUX,JC1,CLABEL1,RC1,KEYNR,U,V,T,P,Q,ODATA,
     *  GDATA,GHDATA,BLDATA,
     2  TTOLD,QTOLD,SVLHX,RHSAV,WM,CLDSAV,
     *  RQT,SRHR,TRHR,TSFREZ
      CLOSE (iu_AIC)
      DO J=1,JM
         DO I=1,IM
            SNOAGE(I,J,1)=GDATA(I,J,11)
            SNOAGE(I,J,2)=GDATA(I,J,11)
         END DO
      END DO
      redoGH=.TRUE.
C**** INITIALIZE VERTICAL SLOPES OF T,Q
  327 DO J=1,JM
         DO I=1,IM
            RDSIG=(SIG(1)-SIGE(2))/(SIG(1)-SIG(2))
            TZ(I,J,1)=(T(I,J,2)-T(I,J,1))*RDSIG
            QZ(I,J,1)=(Q(I,J,2)-Q(I,J,1))*RDSIG
            IF (Q(I,J,1)+QZ(I,J,1).LT.0.) QZ(I,J,1) = -Q(I,J,1)
            DO L=2,LM-1
               RDSIG=(SIG(L)-SIGE(L+1))/(SIG(L-1)-SIG(L+1))
               TZ(I,J,L)=(T(I,J,L+1)-T(I,J,L-1))*RDSIG
               QZ(I,J,L)=(Q(I,J,L+1)-Q(I,J,L-1))*RDSIG
               IF (Q(I,J,L)+QZ(I,J,L).LT.0.) QZ(I,J,L) = -Q(I,J,L)
            END DO
            RDSIG=(SIG(LM)-SIGE(LM+1))/(SIG(LM-1)-SIG(LM))
            TZ(I,J,LM)=(T(I,J,LM)-T(I,J,LM-1))*RDSIG
            QZ(I,J,LM)=(Q(I,J,LM)-Q(I,J,LM-1))*RDSIG
            IF (Q(I,J,LM)+QZ(I,J,LM).LT.0.) QZ(I,J,LM) = -Q(I,J,LM)
         END DO
      END DO
 350  iniPBL=.TRUE.
C**** Set TAU to TAUI for initial starts
  398 IF (TAUI.LT.0.) TAUI=TAUX
      TAU=TAUI
  399 REWIND 9
      WRITE (6,'(A,I4,F11.2,3X,A/)')
     *  '0ATMOSPHERIC I.C. ISTART,TAUX=',ISTART,TAUX,CLABEL1(1:80)
C**** Check consistency of TAU and TAUX (from IC)
      IF (MOD(TAUX-TAU,8760.D0).NE.0.) THEN
         WRITE (6,'(A,2F11.2)')
     *     '0 *******  I.C. and TAU are inconsistent', TAUX,TAU
         IF (TAUX.NE.TAUP) STOP 'ERROR: TAU I.C. and TAU inconsistent'
      ENDIF
      GO TO 500
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10-13
C****
C**** CHOOSE DATA SET TO RESTART ON
  400 TAU1=-1.
      READ (1,ERR=410) TAU1
  410 REWIND 1
      TAU2=-1.
      READ (2,ERR=420) TAU2
  420 REWIND 2
      KDISK=1
      IF (TAU1+TAU2.LE.-2.) GO TO 850
      IF (TAU2.GT.TAU1) KDISK=2
      IF (ISTART.GE.13) KDISK=3-KDISK
      GO TO 450
  440 KDISK=ISTART-10
C**** RESTART ON UNIT KDISK
  450 KDISK0=KDISK

 460  call io_rsf(KDISK0,TAUX,ioread,ioerr)
      if (ioerr.eq.1) goto 840

      KDISK=KDISK0
      IF (RUNID.NE.XLABEL(1:4)) THEN
         WRITE (6,'(4A)')
     *     ' THIS RESTART FILE IS FOR RUN',XLABEL(1:4),' NOT RUN',RUNID
         STOP 'ERROR: WRONG RESTART FILES, MISMATCHED LABELS'
      ENDIF
c      IF (TAUX.NE.TAUY) GO TO 860
      WRITE (6,'(A,I2,A,F11.2,A,A/)') '0RESTART DISK READ, UNIT',
     *  KDISK,', TAUX=',TAUX,' ',CLABEL(1:80)
      IF (ISTART.GT.10) KDISK=3-KDISK
      TAU=TAUX
      TAUP=TAUX
C**** UPDATE C ARRAY FROM INPUTZ
CC500 READ (NLREC,INPUTZ)
  500 REWIND 8
      WRITE(8,'(A)') (NLREC(IR),IR=1,IREC)
      REWIND 8
      READ (8,NML=INPUTZ)
      REWIND 8
      IF (TAU.GT.TAUE) STOP 13
      IF (IM0.LT.IM.OR.JM0.LT.JM.OR.LM0.LT.LM) THEN
      WRITE (6,'('' ARRAY-DIMENSIONS IM,JM,LM '',3I3,
     *  '' ARE INSUFFICIENT FOR IM,JM,LM='',3I3)') IM0,JM0,LM0,IM,JM,LM
      STOP ' ERROR IN GRID SIZE DIMENSIONS '
      END IF
      IF (ISTART.GE.10 .AND. TAU.LT.TAUP) GO TO 900
         IF (USESLP.LE.0.) GO TO 515
C****    REPOSITION THE SEA LEVEL PRESSURE HISTORY DATA SET (UNIT 16)
         REWIND 16
  510    READ (16,ERR=870,END=880) TAU4,((XX4,I=1,IM),J=1,JM),TAUY4
         IF (TAU4.NE.TAUY4) GO TO 860
         IF(TAU.LT.TAU4-.5*USESLP) REWIND 16
         IF (TAU.GE.TAU4+USESLP) GO TO 510
         WRITE (6,'(A,F11.2/)')
     *   '0SLP HISTORY REPOSITIONED.  LAST TAU READ WAS',TAU4
  515    CONTINUE
      IF (USET.LE.0.) GO TO 600
C**** ZERO OUT ARRAYS
      OA = 0.
C**** REPOSITION THE OUTPUT TAPE ON UNIT 20 FOR RESTARTING
      REWIND 20
      IF (TAU.LT.TAUO+USET) GO TO 600
  520 READ (20,ERR=870,END=880) TAUX
      IF (TAU.GE.TAUX+USET) GO TO 520
      WRITE (6,'(A,F11.2/)')
     *  '0OUTPUT TAPE REPOSITIONED.  LAST TAU READ WAS',TAUX
C****
C**** CONSTANT ARRAYS TO BE CALCULATED OR READ IN EACH RUN
C****
  600 IF (KEYCT.LE.1) KEYNR=0
C**** CALCULATE SPHERICAL GEOMETRY
      CALL GEOM_B
C**** CALCULATE DSIG AND DSIGO
      DO L=1,LM
         DSIG(L)=SIGE(L)-SIGE(L+1)
         BYDSIG(L)=1./DSIG(L)
      END DO
      DO L=1,LM-1
         DSIGO(L)=SIG(L)-SIG(L+1)
      END DO

C***  READ IN LANDMASKS AND TOPOGRAPHIC DATA
      call getunit("TOPO",iu_TOPO,.TRUE.)
      
      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      CALL READT (iu_TOPO,0,FLAKE,IM*JM,FLAKE,1)   ! Lake fraction
      CALL READT (iu_TOPO,0,FEARTH,IM*JM,FEARTH,1) ! Earth fraction (no LI)
      CALL READT (iu_TOPO,0,FLICE,IM*JM,FLICE,1)   ! Land ice fraction  
      FLAND = FEARTH + FLICE                       ! Land fraction
C**** DON'T  adjust Land ice fraction to be fraction only over land
c      FLICE = (FLICE/(FLAND+1.D-20))
      CALL READT (iu_TOPO,0,ZATMO,IM*JM,ZATMO,1)   ! Topography
      ZATMO = ZATMO*GRAV                           ! Geopotential
      REWIND iu_TOPO

C**** Init pbl (and read in file containing roughness length data)
      call pblini(iniPBL)

C**** read in ocean heat transports and max. mixed layer depths 
      IF(KOCEAN.EQ.1) CALL OHT_INIT

C**** READ IN VEGETATION DATA SET: VDATA AND VADATA
      call getunit("VEG",iu_VEG,.TRUE.)
      DO K=1,11
         CALL READT (iu_VEG,0,VDATA(1,1,K),IM*JM,VDATA(1,1,K),1)
      END DO
      CLOSE (iu_VEG)
C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS
C**** Recompute GHDATA if redoGH (new soils data)
C****
      CALL GHINIT (DT*NDYN/NSURF,redoGH)
      IF (redoGH) THEN
        WRITE (*,*) 'GHDATA WAS MADE FROM GDATA'
C****   Copy Snow age info into GDATA array
        GDATA(:,:, 9)=SNOAGE(:,:,1)
        GDATA(:,:,10)=SNOAGE(:,:,2)
      END IF

      IF(IRAND.LT.0.AND.TAU.EQ.TAUI) THEN
C****   Perturb the Initial Temperatures by at most 1 degree C
        IRAND=-IRAND
        CALL RINIT (IRAND)
        DO L=1,LM
           DO J=1,JM
              DO I=1,IM
                 TIJL=T(I,J,L)*EXPBYK(P(I,J)*SIG(L)+PTOP)-1.+2*RANDU(X)
                 T(I,J,L)=TIJL/EXPBYK(P(I,J)*SIG(L)+PTOP)
              END DO
           END DO
        END DO
        WRITE(6,*) 'Initial conditions were perturbed !!',IRAND
        IRAND=123456789
      END IF
      CALL RINIT (IRAND)
      CALL FFT0 (IM)
      CALL init_CLD
      CALL init_DIAG
C**** MAKE NRAD A MULTIPLE OF NCNDS
      NRAD=(MAX(NRAD,NCNDS)/NCNDS)*NCNDS
      IF (KDIAG(2).EQ.9.AND.SKIPSE.EQ.0..AND.KDIAG(3).LT.9) KDIAG(2)=8
         KACC0=KACC
      WRITE (6,INPUTZ)
      WRITE (6,'(A6,12I6)') "IDACC=",(IDACC(I),I=1,12)
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  800 WRITE (6,'(A,I4)')
     *  '0ERROR ENCOUNTERED READING AIC ISTART=', ISTART
      STOP 'READ ERROR FOR AIC'
  810 WRITE (6,'(A,2F11.2)')
     *  '0EOF ON AIC.  LATER I.C. NEEDED. TAUP,TAUX=', TAUP,TAUX
      STOP 'ERROR: ALL TAUS<TAUP ON I.C. FILE FOR AIC'
  830 WRITE(6,*) 'READ ERROR FOR GIC: GDATA,GHDATA'
      STOP 'READ ERROR FOR GIC'
  840 IF (3-KDISK.EQ.KLAST) GO TO 850
      REWIND KDISK
      KLAST=KDISK
      KDISK=3-KDISK
      WRITE (6,'(A,I3/,A,I1)')
     *  '0ERROR ENCOUNTERED READING RESTART TAPE ON UNIT',KLAST,
     *  '  TRY TO RESTART THE JOB WITH ISTART=3,KDISK=',KDISK
      GO TO 450
  850 WRITE (6,'(A)')
     *  '0ERRORS ON BOTH RESTART DATA SETS. TERMINATE THIS JOB'
      STOP 'ERRORS ON BOTH RESTART FILES'
  860 WRITE (6,'(A,2F11.2/,A)')
     *  '0TAUX,TAUY=',TAUX,TAUY,'0DISK RESTART FILE DESTROYED'
  861 KDISK0=3-KDISK0
      IF (ISTART.EQ.57) GO TO 850
      ISTART=57
      GO TO 460
  870 WRITE (6,'(A,2F11.2)') '0ERROR ENCOUNTERED REPOSITIONING UNIT 16 O
     *R 20.    TAUX,TAU=', TAUX,TAU
      IF (ISTART.GE.10) GO TO 861
      STOP 'READ ERROR ON OUTPUT FILE ON UNIT 16 OR 20'
  880 WRITE (6,'(A,2F11.2)') '0EOF ON UNIT 16 OR 20 WHILE REPOSITIONING.
     *  TAUX,TAU=', TAUX,TAU
      IF (ISTART.GE.10) GO TO 861
      STOP 'POSITIONING ERROR: EOF REACHED ON UNIT 16 OR 20'
  890 WRITE (6,'(A,I5)') '0INCORRECT VALUE OF ISTART',ISTART
      STOP 'ERROR: ISTART-SPECIFICATION INVALID'
  900 WRITE (6,'(A,F11.2,A,F11.2,A)')
     *  '0PREVIOUS TAUE=',TAUP,' WAS NOT YET REACHED. TAU=',TAU,
     *  ' RESUBMIT JOB WITH AN EARLIER TAUE CARD'
      STOP 'ERROR: TAUE TOO LARGE, SINCE TAU<TAUP'
      END

      SUBROUTINE DAILY
C****
C**** THIS SUBROUTINE PERFORMS THOSE FUNCTIONS OF THE PROGRAM WHICH
C**** TAKE PLACE AT THE BEGINNING OF A NEW DAY.
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,rhow,edpery,orbit
      USE E001M12_COM
      USE GEOM
      USE SLE001
     &  , only : cosday=>cost, sinday=>sint
      USE PBLCOM, only : bldata
      USE DYNAMICS, only : CALC_AMPK
      USE OCEAN, only : OCLIM0,OCLIM,KLAKE

      IMPLICIT NONE

      REAL*8 DOZ1O,DELTAP,PBAR,SPRESS,PEARTH,WFC1,SMASS,TSAVG
      INTEGER I,J,IMAX

      CHARACTER*4 :: AMONTH(12) = (/
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)
      INTEGER :: JDOFM(13) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      INTEGER, PARAMETER :: JDPERY = 365, JMPERY = 12
      REAL*8 LAM

C**** ORBITAL PARAMETERS FOR EARTH FOR YEAR 2000 A.D.
      REAL*8, PARAMETER :: OMEGT = 282.9,OBLIQ=23.44,ECCN=.0167

C****
C**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS
C****
C**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE
  100 SMASS=0.
      DO J=1,JM
         SPRESS=0.
         DO I=1,IM
            SPRESS=SPRESS+P(I,J)
         END DO
         SMASS=SMASS+SPRESS*DXYP(J)
      END DO
      PBAR=SMASS/AREAG+PTOP
C**** CORRECT PRESSURE FIELD FOR ANY LOSS OF MASS BY TRUNCATION ERROR
      DELTAP=PSF-PBAR
      P=P+DELTAP

      CALL CALC_AMPK(LS1-1)

      WRITE (6,'(A25,F10.6/)') '0PRESSURE ADDED IN GMP IS',DELTAP
      DOZ1O=1.
      KLAKE=0
C****
C**** CALCULATE THE DAILY CALENDAR
C****
  200 JYEAR=IYEAR+(IDAY-1)/JDPERY
      JDAY=IDAY-(JYEAR-IYEAR)*JDPERY
      DO MONTH=1,JMPERY
      IF (JDAY.LE.JDOFM(MONTH+1)) GO TO 220
      END DO
  220 JDATE=JDAY-JDOFM(MONTH)
      JMONTH=AMONTH(MONTH)
C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
      CALL ORBIT (OBLIQ,ECCN,OMEGT,DFLOAT(JDAY)-.5,RSDIST,SIND,COSD,LAM)
C****
C**** FIND LEAF-AREA INDEX & WATER FIELD CAPACITY FOR GROUND LAYER 1
C****
      COSDAY=COS(TWOPI/EDPERY*JDAY)
      SINDAY=SIN(TWOPI/EDPERY*JDAY)
      DO J=1,JM
         DO I=1,IM
            PEARTH=FEARTH(I,J)
            WFCS(I,J)=24.
            IF (PEARTH.GT.0.) THEN
               CALL GHINIJ(I,J,WFC1)
               WFCS(I,J)=RHOW*WFC1 ! canopy part changes
            END IF
         END DO
      END DO
   
C      IF (KOCEAN.EQ.1) GO TO 500 ! Lake Ice is prescribed for KOCEAN=1
      CALL OCLIM(DOZ1O)

      RETURN
C****
      ENTRY DAILY0
C**** DO Z1O COMPUTATION ONLY IF BLDATA(I,J,5) DOES NOT CONTAIN Z1O
      DOZ1O=0.
      IF (BLDATA(IM,1,5).NE.-9999.) DOZ1O=1.

      CALL OCLIM0
      KLAKE=1
      IF (TAU.GT.TAUI+DT/7200.) GO TO 200
      KLAKE=0
      GO TO 100
C*****
      END

      SUBROUTINE CHECKT (SUBR)
!@sum  CHECKT Checks arrays for NaN/INF and reasonablness
!@auth Original Development Team
!@ver  1.0 

C**** CHECKT IS TURNED ON BY SETTING QCHECK=.TRUE. IN NAMELIST
C**** REMEMBER TO SET QCHECK BACK TO .FALSE. AFTER THE ERRORS ARE
C**** CORRECTED.

      USE E001M12_COM
      USE DAGCOM, only : QCHECK
      IMPLICIT NONE
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR 

      IF (QCHECK) THEN
C**** Check all prog. arrays for Non-numbers
         CALL CHECK3(U,IM,JM,LM,SUBR,'u ')
         CALL CHECK3(V,IM,JM,LM,SUBR,'v ')
         CALL CHECK3(T,IM,JM,LM,SUBR,'t ')
         CALL CHECK3(Q,IM,JM,LM,SUBR,'q ')
         CALL CHECK3(P,IM,JM,1,SUBR,'p ')
         CALL CHECK3(WM,IM,JM,LM,SUBR,'wm')
         CALL CHECK3(GDATA,IM,JM,16,SUBR,'gd')
C**** Check PBL arrays      
         CALL CHECKPBL(SUBR)
C**** Check Ocean arrays      
         CALL CHECKO(SUBR)
C**** Check Earth arrays      
         CALL CHECKE(SUBR)
      END IF
      RETURN
      END

      SUBROUTINE CHECK3(A,IN,JN,LN,SUBR,FIELD)
!@sum  CHECK3 Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: IN,JN,LN
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*2, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(IN,JN,LN),INTENT(IN) :: A

      INTEGER I,J,L !@var I,J,L loop variables

      DO L=1,LN
         DO J=1,JN
            DO I=1,IN
               IF(.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.)) THEN
                  WRITE (6,*) FIELD,': ',I,J,L,A(I,J,L),'after '
     *                 ,SUBR
                  IF (J.LT.JN.AND.J.GT.1) STOP 'CHECK3'
               END IF
            END DO
         END DO
      END DO
      RETURN
      END

C**** Temporary io_rsf

      SUBROUTINE io_rsf(kunit,TAU,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
      USE E001M12_COM, only : JC,CLABEL,RC,U,V,T,P,Q,WM,GDATA
      USE SOMTQ_COM
      USE GHYCOM, only : wbare,wvege,htbare,htvege,snowbv
      USE RADNCB, only : S0X,CO2,RQT,SRHR,TRHR,FSF
      USE CLD01_COM_E001, only : U00wtr,U00ice,LMCM,TTOLD,QTOLD,SVLHX
     *     ,RHSAV,CLDSAV
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,
     *     cq=>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,mld,usavg,vsavg
     *        ,tauavg,ustar
      USE DAGCOM 
      USE OCEAN, only : ODATA,OA,T50

      IMPLICIT NONE
!@var iaction flag for reading or writing rsf file
!@var kunit Fortran unit number of file i/o
      INTEGER, INTENT(IN) :: iaction,kunit
!@param IOWRITE,IOREAD Values for writing or reading file
      INTEGER, PARAMETER :: IOWRITE=-1, IOREAD=1
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(OUT) :: IOERR
!@var TAU1,TAU2 hours for correct reading check
      REAL*8 TAU1,TAU2
!@var TAU input/ouput value of hour
      REAL*8, INTENT(INOUT) :: TAU 
!@var HEADER Character string label for individual restart files
      CHARACTER HEADER*7
      ioerr=1
      rewind (kunit)

C**** headers are introduced so that individual modules will be able to
C**** tell whether the version of the model variables is current 
      select case (iaction) 
      case (:iowrite)
         WRITE (kunit,err=10) TAU,JC,CLABEL,RC
C**** need a blank line to fool 'qrsfnt' etc.
         WRITE (kunit,err=10) 
         WRITE (kunit,err=10) "E001M12",U,V,T,P,Q,WM
         WRITE (kunit,err=10) "OCN01  ",ODATA,OA,T50
         WRITE (kunit,err=10) "ERT01  ",GDATA
         WRITE (kunit,err=10) "SOL01  ",wbare,wvege,htbare,htvege,snowbv
         WRITE (kunit,err=10) "BLD01  ",wsavg,tsavg,qsavg,dclev,mld
     *        ,usavg,vsavg,tauavg,ustar
         WRITE (kunit,err=10) "PBL01  ",uabl,vabl,tabl,qabl,eabl,cm,ch
     *        ,cq,ipbl
         WRITE (kunit,err=10) "CLD01  ",U00wtr,U00ice,LMCM,TTOLD,QTOLD
     *        ,SVLHX,RHSAV,CLDSAV
         WRITE (kunit,err=10) "QUS01  ",TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ
     *        ,QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
         WRITE (kunit,err=10) "RAD01  ",S0X,CO2,RQT,SRHR,TRHR,FSF
         WRITE (kunit,err=10) "DAG01  ",TSFREZ,AJ,BJ,CJ,RJ,APJ,AJL,ASJL
     *        ,AIJ,AIL,AIJG,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,AJK
     *        ,AIJK,AIJL,AJLSP,TDIURN,KEYNR  
         WRITE (kunit,err=10) TAU
         ioerr=-1
      case (ioread:)
         READ (kunit,err=10) TAU1,JC,CLABEL,RC
         READ (kunit,err=10) 
         READ (kunit,err=10) HEADER,U,V,T,P,Q,WM
         READ (kunit,err=10) HEADER,ODATA,OA,T50
         READ (kunit,err=10) HEADER,GDATA
         READ (kunit,err=10) HEADER,wbare,wvege,htbare,htvege,snowbv
         READ (kunit,err=10) HEADER,wsavg,tsavg,qsavg,dclev,mld,usavg
     *        ,vsavg,tauavg,ustar 
         READ (kunit,err=10) HEADER,uabl,vabl,tabl,qabl,eabl,cm,ch,cq
     *        ,ipbl
         READ (kunit,err=10) HEADER,U00wtr,U00ice,LMCM,TTOLD,QTOLD,SVLHX
     *        ,RHSAV,CLDSAV
         READ (kunit,err=10) HEADER,TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     *        QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
         READ (kunit,err=10) HEADER,S0X,CO2,RQT,SRHR,TRHR,FSF
         READ (kunit,err=10) HEADER,TSFREZ,AJ,BJ,CJ,RJ,APJ,AJL,ASJL,AIJ
     *        ,AIL,AIJG,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,AJK,AIJK
     *        ,AIJL,AJLSP,TDIURN,KEYNR 
         READ (kunit,err=10) TAU2
         IF (TAU1.ne.TAU2) THEN
            STOP "PROBLEM READING RSF FILE"
         ELSE
            TAU=TAU1
            ioerr=-1
         END IF
      end select

 10   REWIND kunit
      RETURN
      END

C**** What io_rsf should look like....

C      SUBROUTINE io_rsf(kunit,iaction,ioerr)
C!@sum   io_rsf controls the reading and writing of the restart files
C!@auth  Gavin Schmidt
C!@ver   1.0
C!@calls io_model,io_ocean,io_seaice,io_lakes,io_ground,io_soils,io_bndry,
C!@calls io_pbl,io_clouds,io_radiation,io_diags
C
C      IMPLICIT NONE
C!@var iaction flag for reading or writing rsf file
C!@var kunit Fortran unit number of file i/o
C      INTEGER, INTENT(IN) :: iaction,kunit
C!@var TAU hour of model run
C      REAL*8, INTENT(INOUT) :: TAU
C!@param IOWRITE,IOREAD Values for writing or reading file
C      INTEGER, PARAMETER :: IOWRITE=-1, IOREAD=1
C!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
C      INTEGER, INTENT(OUT) :: IOERR
C!@var TAU1,TAU2 hours for correct reading check
C      REAL*8 TAU1,TAU2
C
C      rewind kunit
C
CC**** For all iaction < 0  ==> WRITE
CC**** For all iaction > 0  ==> READ
CC**** Particular values may produce variations in individual i/o routines
C
C      select case (iaction) 
C      case (:iowrite)
C         WRITE (kunit) TAU
C      case (ioread:)
C         READ (kunit) TAU1
C      end select
C
CC**** Calls to individual i/o routines
C      call io_model    (kunit,iaction,ioerr)
C      call io_ocean    (kunit,iaction,ioerr)
C      call io_seaice   (kunit,iaction,ioerr)
C      call io_lakes    (kunit,iaction,ioerr)
C      call io_ground   (kunit,iaction,ioerr)
C      call io_soils    (kunit,iaction,ioerr)
C      call io_bndry    (kunit,iaction,ioerr)
C      call io_pbl      (kunit,iaction,ioerr)
C      call io_clouds   (kunit,iaction,ioerr)
C      call io_radiation(kunit,iaction,ioerr)
C      call io_diags    (kunit,iaction,ioerr)
C
C      select case (iaction) 
C      case (:iowrite)
C         WRITE (kunit) TAU
C      case (ioread:)
C         READ (kunit) TAU2
C         IF (TAU1.ne.TAU2) STOP "PROBLEM READING RSF FILE"
C         TAU=TAU1
C      end select
C
C      RETURN
C      END
