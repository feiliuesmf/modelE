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

      INTEGER I,J,L,K,LLAB1,KSS6,MSTART,MNOW,MINC,MSUM
     *     ,MODD5D,iowrite,ioerr
      REAL*8 DTIME,PELSE,PDIAG,PSURF,PRAD,PCDNS,PDYN,TOTALT

      CHARACTER CYEAR*4,CMMND*80
      CHARACTER*8 :: LABSSW,OFFSSW = 'XXXXXXXX'

C****
C**** INITIALIZATIONS
C****
      CALL TIMER (MNOW,MINC,MSUM)
      CALL INPUT
      WRITE (3) OFFSSW
      CLOSE (3)
      MSTART=MNOW-MDYN-MCNDS-MRAD-MSURF-MDIAG-MELSE
C**** INITIALIZE TIME PARAMETERS
      NSTEP=(Itime-ItimeI)*NIdyn
         MODD5K=1000

      CALL DAILY(0)
      CALL daily_EARTH(0)
      CALL daily_LAKE(0)
      CALL daily_OCEAN(0)
      CALL CALC_AMPK(LS1-1)
         CALL CHECKT ('INPUT ')

      WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)')
     *   '0NASA/GISS Climate Model (re)started',
     *   'Year',JYEAR,AMON,JDATE,', Hr',JHOUR,
     *   'Internal clock: DTsrc-steps since 1/1/',IYEAR0,ITIME
      CALL TIMER (MNOW,MINC,MELSE)

C****
C**** If run is already done, just produce diagnostic printout
C****
      IF(Itime.GE.ItimeE) then
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
         CALL DIAGKN
         STOP 13          ! no output files are affected
      END IF

C****
C**** MAIN LOOP
C****
      DO WHILE (Itime.lt.ItimeE)

C**** Every Ndisk Time Steps (DTsrc), starting with the first one,
C**** write restart information alternatingly onto 2 disk files
      IF (MOD(Itime-ItimeI,Ndisk).eq.0) THEN
         CALL RFINAL (IRAND)
         iowrite=-1
         call io_rsf(KDISK, dfloat(Itime) ,iowrite,ioerr) ! no dfloat?
         WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *     ' Restart file written on fort.',KDISK,'Year',
     *     JYEAR,AMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
         KDISK=3-KDISK
         CALL TIMER (MNOW,MINC,MELSE)
      END IF

C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY DAY
      IF (MOD(Itime,NDAY).eq.0) THEN
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
         MODD5D=MOD(Itime-ItimeI,NDA5D)
         IF (MODD5D.EQ.0) CALL DIAG5A (2,0)
         IF (MODD5D.EQ.0) CALL DIAG9A (1)
      CALL DYNAM

      CALL CALC_AMPK(LS1-1)

      CALL CHECKT ('DYNAM ')
      CALL TIMER (MNOW,MINC,MDYN)

         IF (MODD5D.EQ.0) CALL DIAG5A (7,NIdyn)           ! ?
         IF (MODD5D.EQ.0) CALL DIAG9A (2)
         IF (MOD(Itime,NDAY/2).eq.0) CALL DIAG7A
C****
C**** INTEGRATE SOURCE TERMS
C****
         IDACC(1)=IDACC(1)+1
         MODRD=MOD(Itime-ItimeI,NRAD)
         MODD5S=MOD(Itime-ItimeI,NDA5S)
         IF (MODD5S.EQ.0) IDACC(8)=IDACC(8)+1
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG5A (1,0)
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG9A (1)
C**** CONDENSATION, SUPER SATURATION AND MOIST CONVECTION
      CALL CONDSE
      CALL CHECKT ('CONDSE ')
      CALL TIMER (MNOW,MINC,MCNDS)
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)           ! ?
         IF (MODD5S.EQ.0) CALL DIAG9A (3)
C**** RADIATION, SOLAR AND THERMAL
      CALL RADIA
      CALL CHECKT ('RADIA ')
      CALL TIMER (MNOW,MINC,MRAD)
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)          ! ?
         IF (MODD5S.EQ.0) CALL DIAG9A (4)
C**** SURFACE INTERACTION AND GROUND CALCULATION
      CALL PRECIP
      CALL CHECKT ('PRECIP')
      CALL SURFCE
      CALL CHECKT ('SURFCE')
      CALL GROUND
      CALL CHECKT ('GROUND')
C**** Calculate river runoff from lake mass
      CALL RIVERF
      CALL CHECKT ('RIVERF')
      CALL DRYCNV
      CALL CHECKT ('DRYCNV')
      CALL TIMER (MNOW,MINC,MSURF)
         IF (MODD5S.EQ.0) CALL DIAG9A (5)
C**** STRATOSPHERIC MOMENTUM DRAG
      CALL SDRAG
      CALL CHECKT ('SDRAG ')
      CALL TIMER (MNOW,MINC,MSURF)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)          ! ?
         IF (MODD5S.EQ.0) CALL DIAG9A (6)
C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(Itime-ItimeI,NFILTR).EQ.0) THEN
         IDACC(10)=IDACC(10)+1
         IF (MODD5S.NE.0) CALL DIAG5A (1,0)
         CALL DIAG9A (1)
      CALL FILTER
      CALL CHECKT ('FILTER')
      CALL TIMER (MNOW,MINC,MDYN)
         CALL DIAG5A (14,NFILTR*NIdyn)
         CALL DIAG9A (7)
      END IF
C****
C**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
C****
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear0
      Jyear=Itime/(Nday*JDperY) + Iyear0  !  calendar year (A.D.)
      Jhour=MOD(Itime*24/NDAY,24)         ! Hour (0-23)
      Nstep=Nstep+NIdyn                   ! counts DT(dyn)-steps

      IF (MOD(Itime,NDAY).eq.0) THEN
         CALL DIAG5A (1,0)
         CALL DIAG9A (1)
      CALL DAILY(1)
      CALL TIMER (MNOW,MINC,MELSE)
         CALL DIAG5A (16,NDAY*NIdyn)                      ! ?
         CALL DIAG9A (8)
      call daily_EARTH(1)
      CALL daily_LAKE(1)
      call daily_OCEAN(1)
         CALL CHECKT ('DAILY ')
      CALL TIMER (MNOW,MINC,MSURF)
      END IF
C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY 24 HOURS
C****
      IF (NOHT.NE.0.) THEN
         IF (MOD(Itime,NDAY).eq.0) THEN
            WRITE (20) ITIME,OA
            ENDFILE 20
            BACKSPACE 20
            IF(NOHT.LT.0) NOHT=-NOHT
            WRITE (6,'(A,78X,A,I8)')
     *           ' oht-info WRITTEN ON UNIT 20',' Time',ITIME
C**** ZERO OUT INTEGRATED QUANTITIES
            OA(:,:,4:12)=0.
         ELSEIF (MOD(Itime,NDAY/2).eq.0) THEN
            call uset_OCEAN
         END IF
         CALL TIMER (MNOW,MINC,MELSE)
      END IF
C****
C**** WRITE SEA LEVEL PRESSURES EVERY NSLP DTSRC-TIME STEPS
C****
      IF (NSLP.NE.0) THEN
         IF (MOD(ITIME, NSLP).eq.0) THEN
            CALL DIAG10(0)
            IF (NSLP.LT.0) NSLP=-NSLP
            WRITE (6,'(A,78X,A,I8)')
     *           ' Sea Level Pressure written to unit 16',' Time',ITIME
         END IF
      END IF
C****
C**** CALL DIAGNOSTIC ROUTINES
C****
      IF (MOD(Itime-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 ?
C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
      IF (NDPRNT(1).LT.0) THEN
         IF (KDIAG(1).LT.9) CALL DIAGJ
         IF (KDIAG(2).LT.9) CALL DIAGJK
         IF (KDIAG(2).LT.9) CALL DIAGJL
         IF (KDIAG(10).EQ.0) CALL DIAGIL
         IF (KDIAG(7).LT.9) CALL DIAG7P
         IF (KDIAG(3).LT.9) CALL DIAGIJ
         IF (KDIAG(9).LT.9) CALL DIAG9P
         IF (KDIAG(5).LT.9) CALL DIAG5P
         IF (KDIAG(4).LT.9) CALL DIAG4
         IF (Itime.LE.ItimeI+1) THEN
            CALL DIAGKN
         ELSE ! RESET THE UNUSED KEYNUMBERS TO ZERO
            KEYNR(1:42,KEYCT)=0
         END IF
         NDPRNT(1)=NDPRNT(1)+1
      END IF

C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES ON NDPRNT-TH DAY OF RUN
      IF (MOD(Itime,NDAY).eq.0) THEN
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
C****    (NORMALLY DONE AT THE END OF A MONTH)
         DO K=1,13
         IF (JDAY.EQ.NDZERO(K)) THEN
C**** PRINT THE KEY DIAGNOSTICS
           CALL DIAGKN
C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
           IF (KCOPY.GT.0) THEN
             WRITE (CYEAR,'(I4)') JYEAR0
             LLAB1 = INDEX(LABEL1,'(') -1
             IF(LLAB1.LT.1) LLAB1=16
C**** KCOPY > 1 : SAVE THE RESTART INFORMATION
             IF (KCOPY.GT.1) THEN
               CALL RFINAL (IRAND)
               OPEN(30,FILE=AMON0(1:3)//CYEAR//'.rsf'//LABEL1(1:LLAB1),
     *              FORM='UNFORMATTED')
               iowrite=-1   ! ? should be different
               call io_rsf(30,Dfloat(Itime),iowrite,ioerr)  ! ? temp

               CLOSE (30)
             END IF
C**** KCOPY > 2 : SAVE THE OCEAN DATA FOR INITIALIZING DEEP OCEAN RUNS
             IF (KCOPY.GT.2) THEN
               OPEN (30,FILE=AMON0(1:3)//CYEAR//'.oda'//LABEL1(1:LLAB1),
     *              FORM='UNFORMATTED')
               WRITE (30) Itime,ODATA,((AIJ(I,J,IJ_TGO2),I=1,IM),J=1,JM)
               CLOSE (30)
             END IF
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
             OPEN (30,FILE=AMON0(1:3)//CYEAR//'.acc'//LABEL1(1:LLAB1),
     *            FORM='UNFORMATTED')
             WRITE (30) Itime,JC,CLABEL,(SNGL(RC(I)),I=1,161),KEYNR,
     *            (SNGL(TSFREZ(I,1,1)),I=1,IM*JM*2),
     *            (SNGL(AJ(I,1)),I=1,KACC),Itime
             CLOSE (30)
           END IF

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
           CALL TIMER (MNOW,MINC,MDIAG)
           TOTALT=.01*(MNOW-MSTART)      ! in seconds
           PDYN  = MDYN/TOTALT
           PCDNS = MCNDS/TOTALT
           PRAD  = MRAD/TOTALT
           PSURF = MSURF/TOTALT
           PDIAG = MDIAG/TOTALT
           PELSE = MELSE/TOTALT
           DTIME = NDAY*TOTALT/(60.*(Itime-Itime0))  ! minutes/day
           WRITE (6,'(/A,F7.2,6(A,F5.1)//)')
     *      '0TIME',DTIME,'(MINUTES)    DYNAMICS',PDYN,
     *      '    CONDENSATION',PCDNS,'    RADIATION',PRAD,'    SURFACE',
     *      PSURF,'    DIAGNOSTICS',PDIAG,'    OTHER',PELSE
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
      If(MOD(Itime,Nssw).eq.0) READ (3,END=210) LABSSW
  210 CLOSE (3)
      IF (LABSSW.EQ.LABEL1(1:8)) THEN
C**** FLAG TO TERMINATE RUN WAS TURNED ON (GOOD OLE SENSE SWITCH 6)
         KSS6=1
         WRITE (6,'("0SENSE SWITCH 6 HAS BEEN TURNED ON.")')
         EXIT
      END IF

      END DO
C****
C**** END OF MAIN LOOP
C****

C**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      iowrite=-1
      call io_rsf(KDISK,Dfloat(Itime),iowrite,ioerr)   ! temporary ?
      WRITE (6,'(A,I1,55X,A4,I5,A5,I3,A4,I3,A,I8)')
     *  ' Restart file written on fort.',KDISK,'Year',JYEAR,
     *        AMON,JDATE,', Hr',JHOUR,'  Internal clock',ITIME

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)
      WRITE (6,'(/////4(1X,33("****")/)//,A,I8
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME
      IF (Itime.ge.ItimeE) STOP 13
      STOP 12                       ! stopped because of SSW6
      END

      BLOCK DATA BDINP
C****
C**** DEFAULT PARAMETERS FOR MODEL COMMON BLOCK
C****
      USE E001M12_COM
     &     , only : im,jm,lm
      USE DAGCOM, only : kacc

      IMPLICIT NONE

      INTEGER ::
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Noht,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH
      INTEGER                  IDUM
      INTEGER, DIMENSION(13) :: NDZERO,NDPRNT
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC

      COMMON /IPARMB/
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Noht,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH, IDUM   ,NDZERO,NDPRNT    ,  IJD6     ,IDACC

      DOUBLE PRECISION ::
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,  PTRUNC,SKIPSE,
     *  RSDIST,SIND,COSD 
      DOUBLE PRECISION, DIMENSION(4+12) :: TAUTR0  ! to keep sig fixed?
      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(161-29-2*LM) :: RDM2 
       COMMON /RPARMB/   
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,  PTRUNC,SKIPSE,
     *  RSDIST,SIND,COSD,        TAUTR0,SIG,SIGE,  RDM2   ! S0, ??


      CHARACTER*4 NAMD6,AMON,AMON0
      CHARACTER*132 XLABEL
      COMMON /TEXT/ XLABEL,NAMD6(4),AMON,AMON0

      DATA IM0,JM0,LM0, KACC0/         ! KTACC0 should be here too ???
     *     IM ,JM ,LM , KACC /,
     *  KOCEAN,KDISK,KEYCT,KCOPY,     IRAND,MFILTR,Ndisk,Noht,Nslp/
     *       1,    1,    1,    2, 123456789,     1,   24,   0,   0/,
     *  Nrad, Nfiltr, NIsurf, Nssw,                         IYEAR0/
     *     5,      2,      2,    1,                           1976/,
     *  NDAa,   NDA5d, NDA5k, NDA5s, NDA4, NDAsf/
     *     7,       7,     7,     7,   24,     1/,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,  MDIAG,MELSE,MODRD,MODD5K,MODD5S/
     *      0,   0,    0,   0,    0,      0,    0,    0,     0,     0/
      DATA  DT, DTsrc/
     *    450., 3600./,
C****
C**** Note:           DT = DTdyn and NIdyn = DTsrc/DTdyn (set in INPUT)
C**** In general      DTxxx = Nxxx*DTsrc  and  DTxxx = DTsrc/NIxxx
C**** except that the time steps related to NDAa, NDA5k, NDAsf are
C**** slightly larger:     NDAa:   NDAa*DTsrc + 2*DT(dyn),
C****                      NDA5k: NDA5k*DTsrc + 2*DT(dyn),
C****                      NDAsf: NDAsf*DTsrc + DTsrc/NIsurf
C****
     *  PTOP, PSF, PSDRAG,PTRUNC,SKIPSE/
     *  150.,984.,  500.,    0.,     0./
      DATA SIGE /1.0000000,LM*0./                    ! Define in rundeck
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
      USE E001M12_COM, only : im,jm,lm,wm,u,v,t,p,q,gdata,fearth,fland
     *     ,focean,flake,flice,hlake,zatmo,sig,dsig,sige,dsigo
     *     ,bydsig,xlabel,jc,rc,clabel,NAMD6,IJD6,NDPRNT,NDZERO
     *     ,SKIPSE,KEYCT,MFILTR,IRAND,PTRUNC,PSF,PTOP
     *     ,XCDLM,NDASF,NDA4,NDA5S,NDA5K,NDA5D,NDAA,NFILTR
     *     ,NIsurf,NRAD,NIdyn,NDAY,dt,dtsrc,KDISK
     *     ,Iyear0,Itime,ItimeI,ItimeE,  Noht,Nslp,Ndisk,Nssw,Kcopy
     *     ,KOCEAN,LS1,PSFMPT,PSTRAT,KACC0,KTACC0,IDACC,IM0,JM0,LM0
     *     ,VDATA,amonth,JDendOfM,JDperY  ,amon,amon0
      USE SOMTQ_COM
      USE GEOM, only : geom_b
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
     &  ,TITREG,NAMREG,HR_in_DAY,iwrite,jwrite,itwrite,qcheck
      USE DYNAMICS, only : FILTER,CALC_AMPK
      USE OCEAN, only : odata,OA
      USE LAKES_COM, only : T50
      USE FILEMANAGER, only : getunit

      IMPLICIT NONE
!@var iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG

      INTEGER I,J,L,K,KLAST,KDISK0,ITYPE,IM1,KTACC     ! ? ktacc ?
     *     ,IR,IREC,NOFF,ioread,ioerr,itime1,itime2
      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
     *  ISTART=10, HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHRE=-1
      REAL*8 TIJL,X,RDSIG,CDM,TEMP   ! ,SNOAGE ? obsolete
      REAL*4 XX4
      REAL*8 TAU2,TAUX,TAUY,TAU1   ! ? temporary for compatibility only

      INTEGER JC1(100)
      REAL*8 RC1(161)
      INTEGER :: LRUNID=4                       ! RUNID longer than 4?
      CHARACTER*156 CLABEL1

cc?   COMMON/WORK1/SNOAGE(IM,JM,2)

      LOGICAL :: redoGH = .FALSE.,iniPBL = .FALSE., inilake = .FALSE.,
     &           iniSNOW = .FALSE.  ! true = restart from "no snow" rsf

      CHARACTER NLREC*80
C****    List of parameters that CANNOT be changed during a run:
      NAMELIST/INPUTZ/ KOCEAN,PTOP,PSF,SIGE,SIG,LS1,DTsrc
C****    List of parameters that COULD be changed during a run:
     *     ,ISTART,DT,  NRAD,NIsurf,NFILTR, MFILTR
     *     ,SKIPSE,  NDAA,NDA5D,NDA5K,NDA5S,NDA4,NDASF
     *     ,U00wtr,U00ice,LMCM, S0X,CO2,XCDLM,IRAND ! ,PTRUNC obsolete
     *     ,Nslp,Noht,KCOPY,Ndisk, Nssw        ,keyct ! keyct??
     *     ,KDIAG,NDZERO,NDPRNT,IJD6,NAMD6
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK
     *     ,HOURE,DATEE,MONTHE,YEARE, IYEAR0
C****    List of parameters that are disregarded at restarts
     *     ,HOURI,DATEI,MONTHI,YEARI
C****
C**** More default settings
C****
      TEMP=250.
      P(:,:)=PSF-PTOP
      TSAVG(:,:)=TEMP
      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=TEMP  ! will be changed to pot.temp later
      Q(:,:,:)=3.D-6
C**** ADVECTION TERMS FOR SECOND ORDER MOMENTS
      DO I=1,IM*JM*LM*9
         TX(I,1,1) = 0.
      END DO
      DO I=1,IM*JM*LM*9
         QX(I,1,1) = 0.
      END DO
C**** Auxiliary clouds arrays
      RHSAV (:,:,:)=.85
      CLDSAV(:,:,:)=0.
      SVLHX (:,:,:)=0.
      WM    (:,:,:)=0.
C****    ocean info saved for ocean heat transport calculations
         OA = 0.
C**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
C****
C**** Print Header and Label (2 lines) from rundeck
C****
      WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(5,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      WRITE (6,'(A,A/)') '0',XLABEL
      CLABEL1=XLABEL                            ! save for istart=9 case
C****
C**** Print and Copy Namelist parameter changes to disk so they may be
C**** read in repeatedly. Then read them in to overwrite the defaults
C****
      DO WHILE (NLREC(1:5).NE.' &END')
         READ  (5,'    (A80)') NLREC
         WRITE (6,'(35X,A80)') NLREC
         WRITE (8,'(A)') NLREC
      END DO
      REWIND 8
      READ (8,NML=INPUTZ)
      REWIND 8

      IF(ISTART.GE.9) GO TO 400
C***********************************************************************
C****                                                               ****
C****                  INITIAL STARTS - ISTART< 0 - 8               ****
C****                                                               ****
C****   Current settings: 1 - from defaults                         ****
C****                     2 - from observed data                    ****
C****                     3   from first E-model M-file, old format ****
C****         not done  4-7 - from old model M-file - not all data  ****
C****                     8 - from current model M-file             ****
C****                                                               ****
C***********************************************************************
C**** get unit for atmospheric initial conditions if needed
      if(ISTART.gt.1) call getunit("AIC",iu_AIC,.TRUE.)
C****
C**** Set the derived quantities NDAY, Itime.., vert. layering, etc
C****
      KTACC0 = KTACC      ! ? should be set in BDINP
      NDAY = 2*NINT(.5*SDAY/DTsrc)  !  # of time steps in a day: even
C**** Correct the time step
c     DTsrc = SDAY/NDAY   ! currently 1 hour
C**** Get Start Time; at least YearI has to be specified
      IhrI = (yearI-Iyear0)*JDperY +
     *       (JDendofM(monthI-1) + dateI-1)*HR_IN_DAY + HourI
      ITimeI = IhrI*NDAY/24  !  internal clock counts DTsrc-steps
      IF (IhrI.lt.0) then
         WRITE(6,*) 'Please set a proper start time; current values:',
     *     'yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
           STOP 'ERROR: No proper start date'
      END IF
C**** Derived quantities related to vertical layering (sige->ls1,sig)
      DO L=1,LM                           ! alt.:pbot,ls1->sige,...??
         IF(SIGE(L).EQ.0.) LS1=L
         SIG(L)=.5*(SIGE(L)+SIGE(L+1))
      END DO
      PSFMPT = PSF-PTOP
      PSTRAT = (PSF-PTOP)*(SIGE(LS1)-SIGE(LM+1))
      JC1=JC          ! save from being overwritten by reading I.C.s
      RC1=RC          ! save from being overwritten by reading I.C.s
C****
C**** Get Ground conditions from a separate file - ISTART=1,2
C****
      IF(ISTART.LE.2) THEN
C**** Set flag to initialise pbl and snow variables
         iniPBL=.TRUE.
         iniSNOW = .TRUE.  ! extract snow data from first soil layer
C**** new: GDATA(8) UNUSED,GDATA(9-11) SNOW AGE OVER OCN.ICE,L.ICE,EARTH
         call getunit("GIC",iu_GIC,.TRUE.)
         READ(iu_GIC,ERR=830) GDATA,GHDATA,(ODATA(I,1,1),I=1,IM*JM*2)
         CLOSE (iu_GIC)
      END IF
C****
C**** Get primary Atmospheric data from NMC tapes - ISTART=2
C****
      IF(ISTART.EQ.2) THEN
C**** Use title of first record to get the date and make sure  ???
C**** it is consistent with IHRI (at least equal mod 8760)     ???
C****            not yet implemented but could easily be done  ???
         CALL READT (iu_AIC,0,P,IM*JM,P,1)             ! Psurf
         DO I=1,IM*JM
         P(I,1)=P(I,1)-PTOP                            ! Psurf -> P
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
         DO L=1,LM  ! alternatively, only read in L=1,LS1 ; skip rest
         CALL READT (iu_AIC,0,Q(1,1,L),IM*JM,Q(1,1,L),1) ! Q
         END DO
         CALL READT (iu_AIC,0,TSAVG(1,1),IM*JM,TSAVG(1,1),1)  ! Tsurf
         CLOSE (iu_AIC)
      END IF
C****
C**** Derive other data from primary data if necessary - ISTART=1,2
C****                                                    currently
      IF(ISTART.LE.2) THEN
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
     *         (U(IM1,J,1)+U(I,J,1)+U(IM1,J+1,1)+U(I,J+1,1))**2
     *        +(V(IM1,J,1)+V(I,J,1)+V(IM1,J+1,1)+V(I,J+1,1))**2)
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
              T(I,J,L)=T(I,J,L)/(SIG(L)*P(I,J)+PTOP)**KAPA
            END DO
            DO L=LS1,LM
              T(I,J,L)=T(I,J,L)/((SIG(L)*(PSF-PTOP)+PTOP)**KAPA)
            END DO
            DO L=1,LM
              TTOLD(I,J,L)=T(I,J,L)
              QTOLD(I,J,L)=Q(I,J,L)
            END DO
         END DO
         END DO
C**** INITIALIZE TSFREZ
         TSFREZ(:,:,1:2)=365.
C**** Initialize surface friction velocity
         DO ITYPE=1,4
         DO J=1,JM
         DO I=1,IM
            USTAR(I,J,ITYPE)=WSAVG(I,J)*SQRT(CDM)
         END DO
         END DO
         END DO
C**** Initiallise T50
         T50=TSAVG
C**** INITIALIZE VERTICAL SLOPES OF T,Q
C***?    or leave them =0 by removing this section
         DO J=1,JM
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
      END IF
C****
C**** I.C FROM OLDER INCOMPLETE MODEL OUTPUT, ISTART=4-7    just hints
C****
C**** Read what's there and substitute rest as needed (as above)
C**** To be implemented as needed. Sometimes it is safer to
C**** combine the ground layers into 2 layers (top 10cm and rest) and
C**** set   redoGH  to .true.  (after major changes in the GH code or
C**** after changing to a new horizontal grid)
C     redoGH=.TRUE.
C**** Set flag to initialise pbl/snow variables if they are not in I.C.
C     iniPBL=.TRUE.  ; iniSNOW = .TRUE.
      SELECT CASE (ISTART)
      CASE (4:7)
         go to 890   !  not available
C****
C**** I.C FROM RESTART FILE WITH COMPLETE DATA        ISTART=3,8
C****
      CASE (3)                  ! old format
         READ (iu_AIC,ERR=800,END=810) TAUX,JC1,CLABEL1,RC1,KEYNR,
     *     U,V,T,P,Q,
     2     ODATA,GDATA,GHDATA,BLDATA,
     *     uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     3     TTOLD,QTOLD,SVLHX,RHSAV,WM,CLDSAV,
     4     TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     5     QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ,
     6     RQT,T50
         CLOSE (iu_AIC)
C**** New additions
      iniSNOW = .TRUE.  ! extract snow data from first soil layer
C****
C****   Data from current type of RESTART FILE     ISTART=8
C****
      CASE (8)
        ioread=3 ! ??? no need to read SRHR,TRHR,FSF.TSFREZ,diag.arrays
        call io_rsf(iu_AIC,TAUX,ioread,ioerr)
        if (ioerr.eq.1) goto 800
        JC=JC1 ; CLABEL=CLABEL1 ; RC=RC1
      END SELECT
C**** Check consistency of starting date
      IF (ISTART.ge.3.and.(MOD(IHRI-nint(TAUX),8760).ne.0)) THEN
          WRITE (6,*) ' Difference in hours between ',
     *       'Starting date and Data date:',MOD(IHRI-nint(TAUX),8760)
          WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
          STOP 'ERROR: start date inconsistent with data'
      END IF
C**** Set flag to initialise lake variables if they are not in I.C.
      IF (ISTART.lt.8) inilake=.TRUE.
C****
C**** Use IRAND<0 to perturb initial temperatures to create ensembles
C****                              perturbation is at most 1 degree C
      IF(IRAND.LT.0) THEN
        IRAND=-IRAND ! in old Random#gen. all seeds were >0 (RANDIBM)
        CALL RINIT (IRAND)
        DO L=1,LM
        DO J=1,JM
        DO I=1,IM
           TIJL=T(I,J,L)*(P(I,J)*SIG(L)+PTOP)**KAPA-1.+2*RANDU(X)
           T(I,J,L)=TIJL/(P(I,J)*SIG(L)+PTOP)**KAPA
        END DO
        END DO
        END DO
        WRITE(6,*) 'Initial conditions were perturbed !!',IRAND
        IRAND=123456789  ! old Rand#gen: all seeds were >0 (RANDIBM)
      END IF
      write(6,*) 'before setting amon/0 dtsrc,dt:',dtsrc,dt
      write(6,*) 'after  setting amon/0 dtsrc,dt:',dtsrc,dt
      WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
     *  '0Model started on',datei,amonth(monthi),yeari,' Hr',houri,
     *  'ISTART =',ISTART,CLABEL1(1:80)

      GO TO 600
C***********************************************************************
C****                                                               ****
C****                  RESTARTS: ISTART > 8                         ****
C****                                                               ****
C****   Current settings: 9 - from own model M-file                 ****
C****                    10 - from later of fort.1 or fort.2        ****
C****                    11 - from fort.1                           ****
C****                    12 - from fort.2                           ****
C****               13 & up - from earlier of fort.1 or fort.2      ****
C****                                                               ****
C***********************************************************************
C****
C****   DATA FROM end-of-month RESTART FILE     ISTART=9
C****                          used for REPEATS and delayed EXTENSIONS
  400 SELECT CASE (ISTART)
      CASE (9)
        ioread=2 ! ??? no need to read diag.arrays
        call io_rsf(iu_AIC,TAUX,ioread,ioerr)
        if (ioerr.eq.1) goto 800
         WRITE (6,'(A,I2,A,F11.2,A,A/)') '0Model restarted; ISTART=',
     *     ISTART,', TAU=',TAUX,' ',CLABEL(1:80)
        JC(32:37)=0
        CLABEL(1:132)=CLABEL1(1:132)
        GO TO 500
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
C****
C**** CHOOSE DATA SET TO RESTART ON
      CASE (10,13:)
         TAU1=-1.                    ! temporarily tau still real*8 ???
         READ (1,ERR=410) TAU1
  410    REWIND 1
         TAU2=-1.
         READ (2,ERR=420) TAU2
  420    REWIND 2
         KDISK=1
         IF (TAU1+TAU2.LE.-2.) GO TO 850
         IF (TAU2.GT.TAU1) KDISK=2
         IF (ISTART.GE.13) KDISK=3-KDISK
      CASE (11,12)
         KDISK=ISTART-10
      END SELECT
C**** RESTART ON UNIT KDISK IN CASE 10-99
      KDISK0=KDISK
C**** After reading (which may change KDISK), reset KDISK such that
C**** the next write action overwrites that SAME file ONLY if ISTART=10
C**** was used WITHOUT PROBLEMS (since then - in case of trouble - we
C**** can go back to the earlier file). In all other cases we want to
C**** first overwrite the other (potentially bad) file. (The most likely
C**** reason not to use ISTART=10 is trouble with the other file.)
      ioread=1
      call io_rsf(KDISK0,TAUX,ioread,ioerr)
      if (ioerr.eq.1) then    ! try the other restart file
         rewind kdisk0
         KDISK=3-KDISK0
         WRITE (6,'(A,I1,A,I1)')
*          ' Read Error on fort.',kdisk0,' trying fort.',kdisk
         KDISK0=KDISK
         call io_rsf(KDISK0,TAUX,ioread,ioerr)
         if (ioerr.eq.1) go to 850
         if(istart.eq.10) KDISK0=3-KDISK
      end if
      KDISK=KDISK0
      IF (istart.gt.10) KDISK=3-KDISK
c     IF (CLABEL1(1:LRUNID).NE.XLABEL(1:LRUNID)) THEN ! used to be vital
c        WRITE (6,*) 'THIS RESTART FILE IS FOR RUN ', !   when rsfs were
c    *   XLABEL(1:LRUNID),' NOT RUN ',CLABEL1(1:LRUNID)   ! shared among
c        STOP 'ERROR: WRONG RESTART FILES, MISMATCHED LABELS'    !  runs
c     ENDIF
      WRITE (6,'(A,I2,A,F11.2,A,A/)') '0RESTART DISK READ, UNIT',
     *   KDISK,', TAUX=',TAUX,' ',CLABEL(1:80)
  500 CONTINUE
C**** UPDATE C ARRAY FROM INPUTZ
      REWIND 8
      READ (8,NML=INPUTZ)
      REWIND 8
c     IF (IM0.NE.IM.OR.JM0.NE.JM.OR.LM0.NE.LM) THEN   ! obsolete
c     WRITE (6,'('' ARRAY-DIMENSIONS IM,JM,LM '',3I3,
c    *  '' ARE INSUFFICIENT FOR IM,JM,LM='',3I3)') IM0,JM0,LM0,IM,JM,LM
c     STOP ' ERROR IN GRID SIZE DIMENSIONS '
c     END IF
      IF (NSLP.GT.0) THEN
C****    REPOSITION THE SEA LEVEL PRESSURE HISTORY DATA SET (UNIT 16)
         REWIND 16
  510    READ (16,ERR=870,END=880) Itime1,((XX4,I=1,IM),J=1,JM),Itime2
         IF (itime1.NE.itime2) THEN
            write(6,*) 'slp history record destroyed; time tags',
     *        itime1,itime2,' inconsistent. Try ISTART=99'
            stop 'slp record bad'
         end if
         IF (itime.LT.itime1) REWIND 16     ! for some false starts
         IF (itime.GE.itime1+NSLP) GO TO 510
         WRITE (6,'(A,I8/)')
     *      '0SLP HISTORY REPOSITIONED.  LAST RECORD:',Itime1
      END IF
      IF (Noht.GT.0) THEN
C****    REPOSITION THE OUTPUT TAPE ON UNIT 20 FOR RESTARTING
         REWIND 20
  520    READ (20,ERR=870,END=880) itime1
         IF (itime.LT.itime1) REWIND 20
         IF (itime.GE.itime1+Noht) GO TO 520
         WRITE (6,'(A,I8/)')
     *        '0VFLX-file for OHT repositioned, last record:',Itime1
      END IF
C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************
  600 CONTINUE
      IF (KEYCT.LE.1) KEYNR=0
C****
C**** Update  ItimE only if YearE is specified in the rundeck
C****
      if(yearE.ge.0) ItimeE = ((yearE-Iyear0)*JDperY +
     *   (JDendofM(monthE-1)+dateE-1)*HR_IN_DAY + HourE)*NDAY/24
C****
C**** Make sure  dtsrc/dt(dyn)  is a multiple of 2
C****
      NIdyn = 2*nint(.5*dtsrc/dt)
      DT = DTsrc/NIdyn
      write(6,*) 'after  redef. of dt dtsrc,dt:',dtsrc,dt
C****
C**** COMPUTE GRID RELATED VARIABLES AND READ IN TIME-INDEPENDENT ARRAYS
C****
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
      CALL CALC_AMPK(LM)

C**** READ SPECIAL REGIONS FROM UNIT 29
      call getunit("REG",iu_REG,.TRUE.)
      READ(iu_REG) TITREG,JREG,NAMREG
      WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG

C***  READ IN LANDMASKS AND TOPOGRAPHIC DATA
      call getunit("TOPO",iu_TOPO,.TRUE.)

      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      CALL READT (iu_TOPO,0,FLAKE,IM*JM,FLAKE,1)   ! Lake fraction
      CALL READT (iu_TOPO,0,FEARTH,IM*JM,FEARTH,1) ! Earth frac. (no LI)
      CALL READT (iu_TOPO,0,FLICE,IM*JM,FLICE,1)   ! Land ice fraction
      FLAND = FEARTH + FLICE                       ! Land fraction
C**** DON'T  adjust Land ice fraction to be fraction only over land
c      FLICE = (FLICE/(FLAND+1.D-20))
      CALL READT (iu_TOPO,0,ZATMO,IM*JM,ZATMO,1)   ! Topography
      ZATMO = ZATMO*GRAV                           ! Geopotential
      CALL READT (iu_TOPO,0,HLAKE,IM*JM,HLAKE,2)   ! Lake Depths
      REWIND iu_TOPO

C**** Initialize pbl (and read in file containing roughness length data)
      CALL init_pbl(iniPBL)
C**** Initialise lake variables (including river directions)
      CALL init_LAKES(inilake)

C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN

C**** READ IN VEGETATION DATA SET: VDATA
      call getunit("VEG",iu_VEG,.TRUE.)
      DO K=1,11
         CALL READT (iu_VEG,0,VDATA(1,1,K),IM*JM,VDATA(1,1,K),1)
      END DO
      CLOSE (iu_VEG)
C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS
C**** Recompute GHDATA if redoGH (new soils data)
C****
      CALL init_GH(DTsrc/NIsurf,redoGH,iniSNOW)
      IF (redoGH) THEN
        WRITE (*,*) 'GHDATA WAS MADE FROM GDATA'
cc?     Copy Snow age info into GDATA array
cc?     GDATA(:,:, 9)=SNOAGE(:,:,1)
cc?     GDATA(:,:,10)=SNOAGE(:,:,2)
      END IF
      CALL RINIT (IRAND)
      CALL FFT0 (IM)
      CALL init_CLD
      CALL init_DIAG
      IF (KDIAG(2).EQ.9.AND.SKIPSE.EQ.0..AND.KDIAG(3).LT.9) KDIAG(2)=8
      WRITE (6,INPUTZ)
      WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
      WRITE (6,'(A14,2I8)') "KACC0,KTACC0=",KACC0,KTACC0
      WRITE (6,'(A14,3I4)') "IM,JM,LM=",IM,JM,LM
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  800 WRITE (6,'(A,I4/" ",A)')
     *  '0ERROR ENCOUNTERED READING AIC ISTART=', ISTART,CLABEL1(1:80)
      STOP 'READ ERROR FOR AIC'
  810 WRITE (6,'(A,2F11.2)')
     *  '0EOF ON AIC.  ISTART=', ISTART
      STOP 'ERROR: End-of-file on AIC'
  830 WRITE(6,*) 'READ ERROR FOR GIC: GDATA,GHDATA'
      STOP 'READ ERROR FOR GIC'
  850 WRITE (6,'(A)')
     *  '0ERRORS ON BOTH RESTART DATA SETS. TERMINATE THIS JOB'
      STOP 'ERRORS ON BOTH RESTART FILES'
  870 WRITE (6,'(A,2i8)') '0READ ERROR WHILE POSITIONING UNIT 16 OR 20;
     *   MODEL TIME, LAST TIME READ ON OUTPUT FILE:', itime,itime1
      STOP 'READ ERROR ON OUTPUT FILE ON UNIT 16 OR 20'
  880 WRITE (6,'(A,2i8)') '0EOF ON UNIT 16 OR 20 WHILE REPOSITIONING.
     *   MODEL TIME, LAST TIME READ ON OUTPUT FILE:', itime,itime1
      STOP 'POSITIONING ERROR: EOF REACHED ON UNIT 16 OR 20'
  890 WRITE (6,'(A,I5)') '0INCORRECT VALUE OF ISTART',ISTART
      STOP 'ERROR: ISTART-SPECIFICATION INVALID'
      END

      SUBROUTINE DAILY(IEND)
!@sum  DAILY performs daily model-related tasks and at start
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : orbit
      USE E001M12_COM
      USE GEOM, only : areag,dxyp
      USE DYNAMICS, only : calc_ampk

      IMPLICIT NONE

      REAL*8 DELTAP,PBAR,SPRESS,SMASS
      INTEGER I,J,IEND,IDOZ1O

      REAL*8 LAM

C**** ORBITAL PARAMETERS FOR EARTH FOR YEAR 2000 A.D.
      REAL*8, PARAMETER :: OMEGT = 282.9,OBLIQ=23.44,ECCN=.0167

      IF (IEND.eq.0.and.Itime.gt.ItimeI) GO TO 200
C****
C**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS
C****
C**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE
      SMASS=0.
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
C****
C**** CALCULATE THE DAILY CALENDAR
C****
  200 JYEAR=IYEAR0+Itime/(Nday*JDperY)
      JDAY=1+Itime/Nday-(JYEAR-IYEAR0)*JDperY
         JMON=1
      DO WHILE (JDAY.GT.JDendOfM(JMON))
         JMON=JMON+1
      END DO
      JDATE=JDAY-JDendOfM(JMON-1)
      AMON=AMONTH(JMON)

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
      CALL ORBIT (OBLIQ,ECCN,OMEGT,DFLOAT(JDAY)-.5,RSDIST,SIND,COSD,LAM)

      RETURN
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
C**** Check Lake arrays
         CALL CHECKL(SUBR)
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

      SUBROUTINE io_rsf(kunit,TAU,iaction,ioerr)     ! ?? TAU->Itime
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
      USE E001M12_COM, only : JC,CLABEL,RC,U,V,T,P,Q,WM,GDATA
      USE SOMTQ_COM
      USE GHYCOM, only : wbare,wvege,htbare,htvege,snowbv,
     &     NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ,HSN_IJ,FR_SNOW_IJ
      USE RADNCB, only : S0X,CO2,RQT,SRHR,TRHR,FSF
      USE CLD01_COM_E001, only : U00wtr,U00ice,LMCM,TTOLD,QTOLD,SVLHX
     *     ,RHSAV,CLDSAV
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,
     *     cq=>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg
     *        ,tauavg,ustar
      USE DAGCOM
      USE OCEAN, only : ODATA,OA,Z1O
      USE LAKES_COM, only : MWL,TLAKE,GML,T50

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
         WRITE (kunit,err=10) "OCN01  ",ODATA,OA,T50,MWL,TLAKE,GML
         WRITE (kunit,err=10) "ERT01  ",GDATA
         WRITE (kunit,err=10) "SOL01  ",wbare,wvege,htbare,htvege,snowbv
         WRITE (kunit,err=10) "SNOW01 ",
     *         NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ,HSN_IJ,FR_SNOW_IJ
         WRITE (kunit,err=10) "BLD01  ",wsavg,tsavg,qsavg,dclev,Z1O
     *        ,usavg,vsavg,tauavg,ustar
         WRITE (kunit,err=10) "PBL01  ",uabl,vabl,tabl,qabl,eabl,cm,ch
     *        ,cq,ipbl
         WRITE (kunit,err=10) "CLD01  ",U00wtr,U00ice,LMCM,TTOLD,QTOLD
     *        ,SVLHX,RHSAV,CLDSAV
         WRITE (kunit,err=10) "QUS01  ",TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ
     *        ,QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
         WRITE (kunit,err=10) "RAD01  ",S0X,CO2,RQT,SRHR,TRHR,FSF
         WRITE (kunit,err=10) "DAG01  ",TSFREZ,AJ,BJ,CJ,AREG,APJ,AJL
     *        ,ASJL,AIJ,AIL,AIJG,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE
     *        ,AJK,AIJK,AIJL,AJLSP,TDIURN,KEYNR
         WRITE (kunit,err=10) TAU
         ioerr=-1
      case (ioread:)
         READ (kunit,err=10) TAU1,JC,CLABEL,RC
         READ (kunit,err=10)
         READ (kunit,err=10) HEADER,U,V,T,P,Q,WM
         READ (kunit,err=10) HEADER,ODATA,OA,T50,MWL,TLAKE,GML
         READ (kunit,err=10) HEADER,GDATA
         READ (kunit,err=10) HEADER,wbare,wvege,htbare,htvege,snowbv
         READ (kunit,err=10) HEADER,
     *         NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ,HSN_IJ,FR_SNOW_IJ
         READ (kunit,err=10) HEADER,wsavg,tsavg,qsavg,dclev,Z1O,usavg
     *        ,vsavg,tauavg,ustar
         READ (kunit,err=10) HEADER,uabl,vabl,tabl,qabl,eabl,cm,ch,cq
     *        ,ipbl
         READ (kunit,err=10) HEADER,U00wtr,U00ice,LMCM,TTOLD,QTOLD,SVLHX
     *        ,RHSAV,CLDSAV
         READ (kunit,err=10) HEADER,TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     *        QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
         READ (kunit,err=10) HEADER,S0X,CO2,RQT,SRHR,TRHR,FSF
         READ (kunit,err=10) HEADER,TSFREZ,AJ,BJ,CJ,AREG,APJ,AJL,ASJL
     *        ,AIJ,AIL,AIJG,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,AJK
     *        ,AIJK,AIJL,AJLSP,TDIURN,KEYNR
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
C!@calls io_model,io_ocean,io_seaice,io_lakes,io_ground,io_soils,
C!@calls io_bndry,io_pbl,io_clouds,io_radiation,io_diags
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
CC**** Particular values may produce variations in indiv. i/o routines
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
C      call io_lakes    (kunit,iaction,ioerr)
C      call io_ice      (kunit,iaction,ioerr)
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
