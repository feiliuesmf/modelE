      PROGRAM GISS_modelE
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  1.0 (Based originally on B399)
      USE CONSTANT, only : sday
      USE E001M12_COM
      USE GEOM
      USE RANDOM
      USE DAGCOM, only : aij,keynr,kdiag,ij_tgo2
      USE DYNAMICS, only : filter,calc_ampk
      USE OCEAN, only : tocean,oa
      USE SEAICE_COM, only : rsi,msi
      IMPLICIT NONE

      INTEGER I,J,L,K,LLAB1,KSS6,MSTART,MNOW,MINC,MSUM,MODD5D,ioerr
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
         call io_rsf(KDISK,Itime,iowrite,ioerr)
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
               call io_rsf(30,Itime,iowrite_mon,ioerr)

               CLOSE (30)
             END IF
C**** KCOPY > 2 : SAVE THE OCEAN DATA FOR INITIALIZING DEEP OCEAN RUNS
C**** Note minor change of format to accomodate replacement of ODATA
             IF (KCOPY.GT.2) THEN
               OPEN (30,FILE=AMON0(1:3)//CYEAR//'.oda'//LABEL1(1:LLAB1),
     *              FORM='UNFORMATTED')
               WRITE (30) Itime,TOCEAN,RSI,MSI,((AIJ(I,J,IJ_TGO2),I=1,IM
     *              ),J=1,JM)
               CLOSE (30)
             END IF
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
             OPEN (30,FILE=AMON0(1:3)//CYEAR//'.acc'//LABEL1(1:LLAB1),
     *            FORM='UNFORMATTED')
             call io_model(30,Itime,iowrite_single,ioerr)
             call io_diags(30,Itime,iowrite_single,ioerr)
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
      call io_rsf(KDISK,Itime,iowrite,ioerr)
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
      USE E001M12_COM, only : im,jm,lm
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
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,SKIPSE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(161-13-2*LM) :: RDM2 
       COMMON /RPARMB/   
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,  SKIPSE,
     *  TAUTR0,SIG,SIGE,  RDM2


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
     *  PTOP, PSF, PSDRAG,SKIPSE/
     *  150.,984.,  500.,     0./
      DATA SIGE /1.0000000,LM*0./                    ! Define in rundeck
      DATA NAMD6 /'AUSD','MWST','SAHL','EPAC'/,
     *  NDZERO/ 0,1,32,60,91,121,152,182,213,244,274,305,335/,
     *  NDPRNT/-1,1,32,60,91,121,152,182,213,244,274,305,335/,
     *  IJD6/63,17, 17,34, 37,27, 13,23/

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
     *     ,bydsig,xlabel,jc,rc,clabel,namd6,ijd6,ndprnt,ndzero
     *     ,skipse,keyct,mfiltr,irand,psf,ptop
     *     ,xcdlm,ndasf,nda4,nda5s,nda5k,nda5d,ndaa,nfiltr
     *     ,nisurf,nrad,nidyn,nday,dt,dtsrc,kdisk
     *     ,iyear0,itime,itimei,itimee,noht,nslp,ndisk,nssw,kcopy
     *     ,kocean,ls1,psfmpt,pstrat,kacc0,ktacc0,idacc,im0,jm0,lm0
     *     ,vdata,amonth,jdendofm,jdpery,amon,amon0,irestart,irerun
     *     ,irsfic,mdyn,mcnds,mrad,msurf,mdiag,melse
      USE SOMTQ_COM
      USE GEOM, only : geom_b
      USE GHYCOM
     &  , only : ghdata
      USE RANDOM
      USE RADNCB, only : rqt,s0x,co2,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav,
     *     U00wtr,U00ice,lmcm
      USE PBLCOM
     &     , only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq=>cqgs
     &     ,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
      USE DAGCOM, only : aj,kacc,tsfrez,tdiurn,kdiag,keynr,jreg
     &  ,titreg,namreg,hr_in_day,iwrite,jwrite,itwrite,qcheck
      USE DYNAMICS, only : filter,calc_ampk
      USE OCEAN, only : tocean,oa,z1o
      USE SEAICE_COM, only : rsi,msi
      USE LAKES_COM, only : t50
      USE FILEMANAGER, only : getunit

      IMPLICIT NONE
!@var iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG

      INTEGER I,J,L,K,KLAST,KDISK0,ITYPE,IM1,KTACC     ! ? ktacc ?
     *     ,IR,IREC,NOFF,ioerr
      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
     *  ISTART=10, HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHRE=-1
      REAL*8 TIJL,X,RDSIG,CDM,TEMP   ! ,SNOAGE ? obsolete
      REAL*4 XX4
      REAL*8 TAUX   ! ? temporary for compatibility only
      INTEGER Itime1,Itime2,ItimeX
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
     *     ,U00wtr,U00ice,LMCM, S0X,CO2,XCDLM,IRAND
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
      TX=0; TY=0; TZ=0; TXX=0; TYY=0; TZZ=0; TXY=0; TZX=0; TYZ=0
      QX=0; QY=0; QZ=0; QXX=0; QYY=0; QZZ=0; QXY=0; QZX=0; QYZ=0
C**** Auxiliary clouds arrays
      RHSAV (:,:,:)=.85
      CLDSAV(:,:,:)=0.
      SVLHX (:,:,:)=0.
      WM    (:,:,:)=0.
C**** Ocean info saved for ocean heat transport calculations
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
C****
C**** Get Ground conditions from a separate file - ISTART=1,2
C****
      IF(ISTART.LE.2) THEN
C**** Set flag to initialise pbl and snow variables
         iniPBL=.TRUE.
         iniSNOW = .TRUE.  ! extract snow data from first soil layer
C**** new: GDATA(8) UNUSED,GDATA(9-11) SNOW AGE OVER OCN.ICE,L.ICE,EARTH
         call getunit("GIC",iu_GIC,.TRUE.)
         READ(iu_GIC,ERR=830) GDATA,GHDATA,((TOCEAN(1,I,J),I=1,IM),J=1
     *        ,JM),RSI
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
            DO K=1,LM_REQ
              RQT(K,I,J)=T(I,J,LM)
            END DO
C**** REPLACE TEMPERATURE BY POTENTIAL TEMPERATURE
            DO L=1,LS1-1
              T(I,J,L)=T(I,J,L)/(SIG(L)*P(I,J)+PTOP)**KAPA
            END DO
            DO L=LS1,LM
              T(I,J,L)=T(I,J,L)/((SIG(L)*(PSF-PTOP)+PTOP)**KAPA)
            END DO
            DO L=1,LM
              TTOLD(L,I,J)=T(I,J,L)
              QTOLD(L,I,J)=Q(I,J,L)
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
C**** Changes to array index orders for greater efficiency (?)
        READ (iu_AIC,ERR=800,END=810) TAUX,JC1,CLABEL1,RC1,KEYNR,
     *     U,V,T,P,Q,
     2     ((TOCEAN(1,I,J),I=1,IM),J=1,JM),RSI,MSI,
     *     (((TOCEAN(L,I,J),I=1,IM),J=1,JM),L=2,3),
     *     GDATA,GHDATA,
     *     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg,ustar,
     *     uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     A     (((TTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     B     (((QTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     C     (((SVLHX(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     D     (((RHSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),WM,
     E     (((CLDSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     4     TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     5     QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ,
     6     (((RQT(L,I,J),I=1,IM),J=1,JM),L=1,LM_REQ),T50
        CLOSE (iu_AIC)
C**** New additions
        iniSNOW = .TRUE.        ! extract snow data from first soil layer
C****
C****   Data from current type of RESTART FILE     ISTART=8
C****
      CASE (8)  ! no need to read SRHR,TRHR,FSF.TSFREZ,diag.arrays
        call io_rsf(iu_AIC,ItimeX,irsfic,ioerr)
        if (ioerr.eq.1) goto 800
      END SELECT
C**** Check consistency of starting date
      IF (ISTART.ge.3.and.(MOD(IHRI-ItimeX,8760).ne.0)) THEN
        WRITE (6,*) ' Difference in hours between ',
     *       'Starting date and Data date:',MOD(IHRI-ItimeX,8760)
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
      CASE (9)    ! no need to read diag.arrays
        call io_rsf(iu_AIC,ItimeX,irerun,ioerr)
        if (ioerr.eq.1) goto 800
         WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *     ISTART,', HOUR=',ItimeX,' ',CLABEL(1:80)
        MDYN=0 ; MCNDS=0 ; MRAD=0 ; MSURF=0 ; MDIAG=0 ; MELSE=0   
c        CLABEL(1:132)=CLABEL1(1:132)
        GO TO 500
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
C****
C**** CHOOSE DATA SET TO RESTART ON
      CASE (10,13:)
         Itime1=-1
         READ (1,ERR=410) Itime1
  410    REWIND 1
         Itime2=-1
         READ (2,ERR=420) Itime2
  420    REWIND 2
         KDISK=1
         IF (Itime1+Itime2.LE.-2.) GO TO 850
         IF (Itime2.GT.Itime1) KDISK=2
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
      call io_rsf(KDISK0,ItimeX,irestart,ioerr)
      if (ioerr.eq.1) then    ! try the other restart file
         rewind kdisk0
         KDISK=3-KDISK0
         WRITE (6,'(A,I1,A,I1)')
*          ' Read Error on fort.',kdisk0,' trying fort.',kdisk
         KDISK0=KDISK
         call io_rsf(KDISK0,ItimeX,irestart,ioerr)
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
      WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
     *   KDISK,', HOUR=',ItimeX,' ',CLABEL(1:80)
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
         IF (Itime1.NE.Itime2) THEN
            write(6,*) 'slp history record destroyed; time tags',
     *        Itime1,Itime2,' inconsistent. Try ISTART=99'
            stop 'slp record bad'
         end if
         IF (Itime.LT.Itime1) REWIND 16     ! for some false starts
         IF (Itime.GE.Itime1+NSLP) GO TO 510
         WRITE (6,'(A,I8/)')
     *      '0SLP HISTORY REPOSITIONED.  LAST RECORD:',Itime1
      END IF
      IF (Noht.GT.0) THEN
C****    REPOSITION THE OUTPUT TAPE ON UNIT 20 FOR RESTARTING
         REWIND 20
  520    READ (20,ERR=870,END=880) Itime1
         IF (Itime.LT.Itime1) REWIND 20
         IF (Itime.GE.Itime1+Noht) GO TO 520
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
C**** Update  ItimeE only if YearE is specified in the rundeck
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
     *   MODEL TIME, LAST TIME READ ON OUTPUT FILE:', Itime,Itime1
      STOP 'READ ERROR ON OUTPUT FILE ON UNIT 16 OR 20'
  880 WRITE (6,'(A,2i8)') '0EOF ON UNIT 16 OR 20 WHILE REPOSITIONING.
     *   MODEL TIME, LAST TIME READ ON OUTPUT FILE:', Itime,Itime1
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
      USE RADNCB, only : RSDIST,COSD,SIND

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
C**** Check Ice arrays
         CALL CHECKI(SUBR)
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

      SUBROUTINE io_rsf(kunit,it,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
!@calls io_model,io_ocean,io_lakes,io_seaice,io_earth,io_soils,io_snow
!@calls io_landice,io_bldat,io_pbl,io_clouds,io_somtq,io_rad,io_diags

      IMPLICIT NONE
!@var iaction flag for reading or writing rsf file
!@var kunit Fortran unit number of file i/o
      INTEGER, INTENT(IN) :: iaction,kunit
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var IT1 hour for correct reading check
      INTEGER IT1

      ioerr=-1
      rewind kunit

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
C**** Particular values may produce variations in indiv. i/o routines

C**** Calls to individual i/o routines
      call io_model  (kunit,it,iaction,ioerr)
      it1=it
      call io_ocean  (kunit,iaction,ioerr)
      call io_lakes  (kunit,iaction,ioerr)
      call io_seaice (kunit,iaction,ioerr)
      call io_earth  (kunit,iaction,ioerr)
      call io_soils  (kunit,iaction,ioerr)
      call io_snow   (kunit,iaction,ioerr)
c      call io_landice(kunit,iaction,ioerr)
      call io_bldat  (kunit,iaction,ioerr)
      call io_pbl    (kunit,iaction,ioerr)
      call io_clouds (kunit,iaction,ioerr)
      call io_somtq  (kunit,iaction,ioerr)
      call io_rad    (kunit,iaction,ioerr)
      call io_diags  (kunit,it,iaction,ioerr)

      if (it1.ne.it) THEN
        WRITE(6,*) "TIMES DO NOT MATCH READING IN RSF FILE",it,it1
        ioerr=1
      END IF
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN RESTART FILE: KUNIT="
     *     ,kunit

      RETURN
      END SUBROUTINE io_rsf
