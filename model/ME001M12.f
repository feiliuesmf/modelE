      PROGRAM GISS_modelE
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  1.0 (Based originally on B399)
      USE MODEL_COM
      USE RANDOM
      USE DAGCOM, only : keynr,kdiag,oa,monacc
      USE FILEMANAGER, only : getunit
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      IMPLICIT NONE

      INTEGER I,J,L,K,M,LLAB1,KSS6,MSTART,MNOW,MODD5D,months,ioerr,Ldate
      INTEGER iu_VFLXO,iu_SLP
      INTEGER :: MDUM = 0
      REAL*8, DIMENSION(NTIMEMAX) :: PERCENT
      REAL*8 DTIME,PELSE,PDIAG,PSURF,PRAD,PCDNS,PDYN,TOTALT

      CHARACTER aDATE*14
      CHARACTER*8 :: LABSSW,OFFSSW = 'XXXXXXXX'
C****
C**** INITIALIZATIONS
C****
         CALL TIMER (MNOW,MDUM)

      CALL INPUT
C****
C**** If run is already done, just produce diagnostic printout
C****
      IF (Itime.GE.ItimeE) then
         WRITE (6,'("1"/64(1X/))')
         IF (KDIAG(1).LT.9) CALL DIAGJ
         IF (KDIAG(2).LT.9) CALL DIAGJK
         IF (KDIAG(2).LT.9) CALL DIAGJL
         IF (KDIAG(10).LT.9) CALL DIAGIL
         IF (KDIAG(7).LT.9) CALL DIAG7P
         IF (KDIAG(3).LT.9) CALL DIAGIJ
         IF (KDIAG(9).LT.9) CALL DIAG9P
         IF (KDIAG(5).LT.9) CALL DIAG5P
         IF (KDIAG(6).LT.9) CALL DIAG6
         IF (KDIAG(4).LT.9) CALL DIAG4
         IF (KDIAG(11).LT.9) CALL diag_RIVER
         CALL DIAGKN
         CALL exit_rc (13)  ! no output files are affected
      END IF
      WRITE (3) OFFSSW
      CLOSE (3)
         MSTART=MNOW
         DO M=1,NTIMEACC
           MSTART= MSTART-TIMING(M)
         END DO
C**** INITIALIZE TIME PARAMETERS
      NSTEP=(Itime-ItimeI)*NIdyn
         MODD5K=1000

      CALL DAILY(0)
      if (itime.eq.itimei) call reset_diag(0)
      CALL daily_EARTH(0)
      CALL daily_OCEAN(0)
      CALL CALC_AMPK(LS1-1)
         CALL CHECKT ('INPUT ')

      WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)')
     *   '0NASA/GISS Climate Model (re)started',
     *   'Year',JYEAR,aMON,JDATE,', Hr',JHOUR,
     *   'Internal clock: DTsrc-steps since 1/1/',IYEAR0,ITIME
         CALL TIMER (MNOW,MELSE)
C****
C**** Open and position output history files if needed
C****
      if (Kvflxo.ne.0) then
        write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
        call getunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
        call IO_pos(iu_VFLXO,Itime,2*im*jm*12,Nday)
      end if
      if (Nslp.ne.0) then
        write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
        call getunit('SLP'//aDATE(1:7),iu_SLP,.true.,.false.)
        call IO_pos(iu_SLP,Itime,im*jm,Nslp)
      end if
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
     *     JYEAR,aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
         KDISK=3-KDISK
         CALL TIMER (MNOW,MELSE)
      END IF

C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY DAY
      IF (MOD(Itime,NDAY).eq.0) THEN
C**** CHECK FOR BEGINNING OF EACH MONTH => RESET DIAGNOSTICS
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
        IF ( months.ge.NMONAV .and. JDAY.eq.1+JDendOfM(Jmon-1) ) then
          call reset_DIAG(0)
          if (Kvflxo.ne.0) then
            close (iu_VFLXO)
            write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
            open (iu_VFLXO,file='VFLXO'//aDATE(1:7),form='unformatted')
          end if
          if (Nslp.ne.0) then
            close (iu_SLP)
            write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
            open (iu_SLP,file='SLP'//aDATE(1:7),form='unformatted')
          end if
        end if
C**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
        call daily_DIAG
      END IF
C****
C**** INTEGRATE DYNAMIC TERMS (DIAGA AND DIAGB ARE CALLED FROM DYNAM)
C****
         MODD5D=MOD(Itime-ItimeI,NDA5D)
         IF (MODD5D.EQ.0) IDACC(7)=IDACC(7)+1
         IF (MODD5D.EQ.0) CALL DIAG5A (2,0)
         IF (MODD5D.EQ.0) CALL DIAG9A (1)
      CALL DYNAM

      CALL CALC_AMPK(LS1-1)
C**** calculate some dynamic variables for the PBL
      CALL PGRAD_PBL

         CALL CHECKT ('DYNAM ')
         CALL TIMER (MNOW,MDYN)

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
         CALL TIMER (MNOW,MCNDS)
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)           ! ?
         IF (MODD5S.EQ.0) CALL DIAG9A (3)
C**** RADIATION, SOLAR AND THERMAL
      CALL RADIA
         CALL CHECKT ('RADIA ')
         CALL TIMER (MNOW,MRAD)
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)          ! ?
         IF (MODD5S.EQ.0) CALL DIAG9A (4)
C****
C**** SURFACE INTERACTION AND GROUND CALCULATION
C****
C**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
C**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
C**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      CALL PRECIP_SI
      CALL PRECIP_LI
C**** APPLY PRECIPITATION AND RUNOFF TO LAKES/OCEANS
      CALL PRECIP_LK
      CALL PRECIP_OC
         IF (MODD5S.EQ.0) CALL DIAG9A (5)
         CALL CHECKT ('PRECIP')
C**** CALCULATE SURFACE FLUXES AND EARTH
      CALL SURFCE
         CALL CHECKT ('SURFCE')
         CALL CHECKT ('GRDNSI')
C**** APPLY SURFACE FLUXES TO SEA/LAKE/LAND ICE
      CALL GROUND_SI
      CALL GROUND_LI
C**** APPLY FLUXES TO LAKES AND DETERMINE ICE FORMATION
      CALL GROUND_LK
         IF (MODD5S.EQ.0) CALL DIAG9A (6)
C**** CALCULATE RIVER RUNOFF FROM LAKE MASS
      CALL RIVERF
      CALL GROUND_E    ! diagnostic only - should be merged with EARTH
C**** APPLY FLUXES TO OCEAN AND DETERMINE ICE FORMATION
      CALL GROUND_OC
         CALL CHECKT ('GRNDOC')
C**** CALCULATE ICE DYNAMICS
      CALL DYNSI
C**** APPLY ICE FORMED IN THE OCEAN/LAKES TO ICE VARIABLES
      CALL FORM_SI
         CALL CHECKT ('FORMSI')
C**** ADVECT ICE
      CALL ADVSI
         CALL CHECKT ('ADVSI ')
C**** if aturb_on=.false., CALCULATE DRY CONVECTION ABOVE PBL
      if(.not. aturb_on) CALL DRYCNV (2,LM-1)
         IF (MODD5S.EQ.0) CALL DIAG9A (7)
         CALL CHECKT ('DRYCNV')
         CALL TIMER (MNOW,MSURF)

C**** CALL OCEAN DYNAMIC ROUTINES
      CALL ODYNAM

C**** STRATOSPHERIC MOMENTUM DRAG
      CALL SDRAG(LM,DTSRC)
         CALL CHECKT ('SDRAG ')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAG9A (9)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)          ! ?
C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(Itime-ItimeI,NFILTR).EQ.0) THEN
           IDACC(10)=IDACC(10)+1
           IF (MODD5S.NE.0) CALL DIAG5A (1,0)
           CALL DIAG9A (1)
        CALL FILTER
           CALL CHECKT ('FILTER')
           CALL TIMER (MNOW,MDYN)
           CALL DIAG5A (14,NFILTR*NIdyn)
           CALL DIAG9A (8)
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
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
           CALL TIMER (MNOW,MELSE)
        call daily_EARTH(1)
        CALL daily_LAKE(1)
        call daily_OCEAN(1)
           CALL CHECKT ('DAILY ')
           CALL TIMER (MNOW,MSURF)
           CALL DIAG5A (16,NDAY*NIdyn)                      ! ?
           CALL DIAG9A (10)
        call flush(6)
      END IF
C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY 24 HOURS
C****
      IF (Kvflxo.NE.0.) THEN
         IF (MOD(Itime,NDAY).eq.0) THEN
            call WRITEI (iu_vflxo,Itime,OA,2*im*jm*12)
C**** ZERO OUT INTEGRATED QUANTITIES
            OA(:,:,4:12)=0.
         ELSEIF (MOD(Itime,NDAY/2).eq.0) THEN
            call vflx_OCEAN
         END IF
         CALL TIMER (MNOW,MELSE)
      END IF
C****
C**** WRITE SEA LEVEL PRESSURES EVERY NSLP DTSRC-TIME STEPS
C****
      IF (NSLP.NE.0) THEN
         IF (MOD(ITIME, NSLP).eq.0) CALL get_SLP(iu_SLP)
      END IF
C****
C**** CALL DIAGNOSTIC ROUTINES
C****
      IF (MOD(Itime-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 ?
C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
      IF (NIPRNT.GT.0) THEN
         IF (KDIAG(1).LT.9) CALL DIAGJ
         IF (KDIAG(2).LT.9) CALL DIAGJK
         IF (KDIAG(2).LT.9) CALL DIAGJL
         IF (KDIAG(10).EQ.0) CALL DIAGIL
         IF (KDIAG(7).LT.9) CALL DIAG7P
         IF (KDIAG(3).LT.9) CALL DIAGIJ
         IF (KDIAG(9).LT.9) CALL DIAG9P
         IF (KDIAG(5).LT.9) CALL DIAG5P
         IF (KDIAG(4).LT.9) CALL DIAG4
         IF (KDIAG(11).LT.9) CALL diag_RIVER
         IF (Itime.LE.ItimeI+1) THEN
            CALL DIAGKN
         ELSE ! RESET THE UNUSED KEYNUMBERS TO ZERO
            KEYNR(1:42,KEYCT)=0
         END IF
         NIPRNT=NIPRNT-1
      END IF

C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
C****    done at the end of (selected) months
      IF (months.ge.NMONAV .and.   ! next 2 conditions are rarely needed
     *    JDAY.eq.1+JDendOfM(JMON-1) .and. MOD(Itime,NDAY).eq.0) THEN
C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
c       WRITE (6,'("1"/64(1X/))')
        IF (KDIAG(1).LT.9) CALL DIAGJ
        IF (KDIAG(2).LT.9) CALL DIAGJK
        IF (KDIAG(2).LT.9) CALL DIAGJL
        IF (KDIAG(10).LT.9) CALL DIAGIL
        IF (KDIAG(7).LT.9) CALL DIAG7P
        IF (KDIAG(3).LT.9) CALL DIAGIJ
        IF (KDIAG(9).LT.9) CALL DIAG9P
        IF (KDIAG(5).LT.9) CALL DIAG5P
        IF (KDIAG(6).LT.9) CALL DIAG6
        IF (KDIAG(4).LT.9) CALL DIAG4
        IF (KDIAG(11).LT.9) CALL diag_RIVER
        CALL DIAGKN

C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
          LLAB1 = INDEX(XLABEL(1:16),'(') -1
          IF (LLAB1.LT.1) LLAB1=16
          if (index(XLABEL(1:LLAB1),' ').gt.0)
     *      LLAB1=index(XLABEL(1:LLAB1),' ')-1
          call aPERIOD (JMON0,JYEAR0,months,1,aDATE(1:12),Ldate)
          WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),JYEAR
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.12) m = m-12
            monacc(m) = 1
          end do
          OPEN (30,FILE=aDATE(1:7)//'.acc'//XLABEL(1:LLAB1),
     *         FORM='UNFORMATTED')
          call io_rsf (30,Itime,iowrite_single,ioerr)
          CLOSE (30)
C**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
          IF (KCOPY.GT.1) THEN
            CALL RFINAL (IRAND)
            OPEN(30,FORM='UNFORMATTED',
     *              FILE='1'//aDATE(8:14)//'.rsf'//XLABEL(1:LLAB1))
            call io_rsf(30,Itime,iowrite_mon,ioerr)
            CLOSE (30)
          END IF
C**** KCOPY > 2 : ALSO SAVE THE OCEAN DATA TO INITIALIZE DEEP OCEAN RUNS
          IF (KCOPY.GT.2) THEN
            OPEN (30,FILE=aDATE(1:7)//'.oda'//XLABEL(1:LLAB1),
     *           FORM='UNFORMATTED')
            call io_oda(30,Itime,iowrite,ioerr)
            CLOSE (30)
          END IF
        END IF

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (MNOW,MDIAG)
        TOTALT=.01*(MNOW-MSTART)      ! in seconds
        DO M=1,NTIMEACC
          PERCENT(M) = TIMING(M)/TOTALT
        END DO
        DTIME = NDAY*TOTALT/(60.*(Itime-Itime0))  ! minutes/day
        WRITE (6,'(/A,F7.2,A,10(A13,F5.1)//)')
     *   '0TIME',DTIME,'(MINUTES) ',(TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        TIMING = 0
        MSTART= MNOW
      END IF

C**** CPU TIME FOR CALLING DIAGNOSTICS
      CALL TIMER (MNOW,MDIAG)
C**** TEST FOR TERMINATION OF RUN
      IF (MOD(Itime,Nssw).eq.0) READ (3,END=210) LABSSW
  210 CLOSE (3)
      IF (LABSSW.EQ.XLABEL(1:8)) THEN
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
      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *  ' Restart file written on fort.',KDISK,'Year',JYEAR,
     *     aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)
      WRITE (6,'(/////4(1X,33("****")/)//,A,I8
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME
      IF (Itime.ge.ItimeE) CALL exit_rc (13)
      CALL exit_rc (12)             ! stopped because of SSW6
      END

      BLOCK DATA BDINP
C****
C**** DEFAULT PARAMETERS FOR MODEL COMMON BLOCK
C****
      USE MODEL_COM, only : im,jm,lm
      USE DAGCOM, only : kacc

      IMPLICIT NONE

      INTEGER ::
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Kvflxo,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH,NIPRNT,NMONAV
      INTEGER, DIMENSION(32) :: IDUM
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC

      COMMON /IPARMB/
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Kvflxo,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH,NIPRNT,NMONAV,  IDUM ,  IJD6     ,IDACC

      DOUBLE PRECISION ::
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,SKIPSE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(161-13-2*LM) :: RDM2
       COMMON /RPARMB/
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,  SKIPSE,
     *  TAUTR0,SIG,SIGE,  RDM2


      CHARACTER*4 NAMD6,aMON,aMON0
      CHARACTER*132 XLABEL
      COMMON /TEXT/ XLABEL,NAMD6(4),aMON,aMON0

      DATA IM0,JM0,LM0, KACC0/         ! KTACC0 should be here too ???
     *     IM ,JM ,LM , KACC /,
     *  KOCEAN,KDISK,KEYCT,KCOPY,     IRAND,MFILTR,Ndisk,Kvflxo,Nslp/
     *       1,    1,    1,    2, 123456789,     1,   24,   0,   0/,
     *  Nrad, Nfiltr, NIsurf, Nssw, NIPRNT, NMONAV,         IYEAR0/
     *     5,      2,      2,    1,      1,      1,           1976/,
     *  NDAa,   NDA5d, NDA5k, NDA5s, NDA4, NDAsf/
     *     7,       7,     7,     7,   24,     1/,
     *  MODRD,MODD5K,MODD5S/
     *      0,     0,     0/
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
     *  IJD6/63,17, 17,34, 37,27, 13,23/

      END

      SUBROUTINE INPUT
C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE CONSTANT, only : grav,kapa,sday,shi,lhm
      USE MODEL_COM, only : im,jm,lm,wm,u,v,t,p,q,fearth,fland
     *     ,focean,flake0,flice,hlake,zatmo,sig,dsig,sige,dsigo
     *     ,bydsig,xlabel,namd6,ijd6,niprnt,nmonav
     *     ,skipse,keyct,mfiltr,irand,psf,ptop
     *     ,xcdlm,ndasf,nda4,nda5s,nda5k,nda5d,ndaa,nfiltr
     *     ,nisurf,nrad,nidyn,nday,dt,dtsrc,kdisk,jmon0,jyear0
     *     ,iyear0,itime,itimei,itimee,Kvflxo,nslp,ndisk,nssw,kcopy
     *     ,kocean,ls1,psfmpt,pstrat,kacc0,ktacc0,idacc,im0,jm0,lm0
     *     ,vdata,aMONTH,jdendofm,jdpery,aMON,aMON0,ioread,irerun
     *     ,ioread_single,irsfic,iowrite_single,ftype,itearth,itlandi
     *     ,mdyn,mcnds,mrad,msurf,mdiag,melse
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : geom_b
      USE RANDOM
      USE RADNCB, only : rqt,s0x,co2,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav,
     *     U00wtr,U00ice,lmcm
      USE PBLCOM
     &     , only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
      USE DAGCOM, only : acc_period,monacc,kacc,tsfrez,kdiag,keynr,jreg
     &  ,titreg,namreg,hr_in_day,iwrite,jwrite,itwrite,qcheck,oa
     &  ,iu_ij,iu_jl,iu_il
      USE LAKES_COM, only : flake
      USE FILEMANAGER, only : getunit,closeunits
      USE TIMINGS, only : timing,ntimeacc

      IMPLICIT NONE
!@var iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG

      INTEGER I,J,L,K,KLAST,KDISK0,ITYPE,IM1,KTACC     ! ? ktacc ?
     *     ,IR,IREC,NOFF,ioerr,Ldate,months,years,mswitch
      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
     *  ISTART=10, HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHOURE=-1
      REAL*8 TIJL,CDM,TEMP,PLTOP(LM),X
      REAL*4 XX4
      INTEGER Itime1,Itime2,ItimeX,IhrX,llab1,iargc
      INTEGER :: LRUNID=4                       ! RUNID longer than 4?

      LOGICAL :: redoGH = .FALSE.,iniPBL = .FALSE., inilake = .FALSE.,
     &           iniSNOW = .FALSE.  ! true = restart from "no snow" rsf
     &           ,iniOCEAN = .FALSE.

      CHARACTER NLREC*80,filenm*100
C****    List of parameters that CANNOT be changed during a run:
      NAMELIST/INPUTZ/ KOCEAN,PTOP,PSF,LS1,DTsrc
C****    List of parameters that COULD be changed during a run:
     *     ,ISTART,DT,  NRAD,NIsurf,NFILTR, MFILTR
     *     ,SKIPSE,  NDAA,NDA5D,NDA5K,NDA5S,NDA4,NDASF
     *     ,U00wtr,U00ice,LMCM, S0X,CO2,XCDLM,IRAND
     *     ,Nslp,Kvflxo,KCOPY,Ndisk, Nssw
     *     ,KDIAG,NIPRNT,NMONAV,IJD6,NAMD6
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK
     *     ,IYEAR0,IHOURE, HOURE,DATEE,MONTHE,YEARE
C****    List of parameters that are disregarded at restarts
     *     ,PLTOP,         HOURI,DATEI,MONTHI,YEARI
C****
C**** More default settings
C****
      LMCM=-1  ! if not set in rundeck it will default to LS1-1
      TEMP=250.
      TSAVG(:,:)=TEMP
      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=TEMP  ! will be changed to pot.temp later
      Q(:,:,:)=3.D-6
C**** Advection terms for first and second order moments
      TMOM(:,:,:,:)=0.
      QMOM(:,:,:,:)=0.
C**** Auxiliary clouds arrays
      RHSAV (:,:,:)=.85
      CLDSAV(:,:,:)=0.
      SVLHX (:,:,:)=0.
      WM    (:,:,:)=0.
C****    Ocean info saved for ocean heat transport calculations
         OA = 0.
C**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
C**** Set global default timing descriptions
C**** Other speciality descriptions can be added/used locally
      NTIMEACC = 0
      CALL SET_TIMER("ATMOS. DYNAM",MDYN)
      CALL SET_TIMER("CONDENSATION",MCNDS)
      CALL SET_TIMER("   RADIATION",MRAD)
      CALL SET_TIMER("     SURFACE",MSURF)
      CALL SET_TIMER(" DIAGNOSTICS",MDIAG)
      CALL SET_TIMER("       OTHER",MELSE)
C****
C**** Vertical coord system: L=1->LS1-1 sigma ; L=LS1->LM const. press
C****                        PSF->PTOP          PTOP->PMTOP
C**** Specify in rundeck:    PLTOP 1->LM (mb) and LS1 (or PTOP)
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
      IF (ISTART.LT.0) GO TO 600  !  just do the diagnostics
      IF (ISTART.GE.9) GO TO 400
C***********************************************************************
C****                                                               ****
C****                  INITIAL STARTS - ISTART: 1 to 8              ****
C****                                                               ****
C****   Current settings: 1 - from defaults                         ****
C****                     2 - from observed data                    ****
C****         not done  3-6 - from old model M-file - not all data  ****
C****                     7 - from converted M-file - no snow model ****
C****                     8 - from current model M-file             ****
C****                                                               ****
C***********************************************************************
C**** get unit for atmospheric initial conditions if needed
      IF (ISTART.gt.1) call getunit("AIC",iu_AIC,.true.,.true.)
C****
C**** Set quantities that are derived from the namelist parameters
C****
      KTACC0 = KTACC      ! ????? not needed
!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc later if necessary
      NDAY = 2*NINT(.5*SDAY/DTsrc)
C**** Get Start Time; at least YearI HAS to be specified in the rundeck
      IhrI = ((yearI-Iyear0)*JDperY +
     *        JDendofM(monthI-1) + dateI-1)*HR_IN_DAY + HourI
      ITimeI = IhrI*NDAY/24  !  internal clock counts DTsrc-steps
      Itime=ItimeI
      IF (IhrI.lt.0) then
        WRITE(6,*) 'Please set a proper start time; current values:',
     *     'yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
        STOP 'INPUT: No proper start date - set IYEARI'
      END IF
C**** The vertical layering
      IF (PLTOP(1).lt.0.) then
         write(6,*) 'Please specify PLTOP(1->12) (mb) and if PTOP=',
     *      PTOP,' is not ok also specify either PTOP or LS1'
         STOP 'INPUT: Vertical layering not defined - PLTOP'
      END IF
      IF (LS1.gt.0) PTOP=PLTOP(LS1-1)
      DO L=1,LM    !  SIGE(1)=1.
        SIGE(L+1)=(PLTOP(L)-PTOP)/(PSF-PTOP)
        IF (SIGE(L).EQ.0.) LS1=L
        SIG(L)=.5*(SIGE(L)+SIGE(L+1))
      END DO
      PSFMPT = PSF-PTOP
      PSTRAT = PSFMPT*(SIGE(LS1)-SIGE(LM+1))
      P(:,:)=PSFMPT
C****
C**** Get Ground conditions from a separate file - ISTART=1,2
C****
      IF (ISTART.LE.2) THEN
C**** Set flag to initialise pbl and snow variables
        iniPBL=.TRUE.
        iniSNOW = .TRUE.  ! extract snow data from first soil layer
        iniOCEAN = .TRUE. ! read in ocean ic
C**** Read in ground initial conditions
        call getunit("GIC",iu_GIC,.true.,.true.)
        ioerr=-1
        call io_ocean  (iu_GIC,ioread,ioerr)
        call io_seaice (iu_GIC,ioread,ioerr)
        call io_earth  (iu_GIC,ioread,ioerr)
        call io_soils  (iu_GIC,ioread,ioerr)
        call io_landice(iu_GIC,ioread,ioerr)
        if (ioerr.eq.1) then
          WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
          STOP "INPUT: GIC READ IN ERROR"
        end if
        close (iu_GIC)
      END IF
C****
C**** Get primary Atmospheric data from NMC tapes - ISTART=2
C****
      IF (ISTART.EQ.2) THEN
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
      IF (ISTART.LE.2) THEN
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
     *         +(V(IM1,J,1)+V(I,J,1)+V(IM1,J+1,1)+V(I,J+1,1))**2)
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
C**** INITIALIZE VERTICAL SLOPES OF T,Q
        call tq_zmom_init(t,q)
      END IF
C****
C**** I.C FROM OLDER INCOMPLETE MODEL OUTPUT, ISTART=3-6    just hints
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
      CASE (3:6)
         go to 890   !  not available
C****
C**** I.C FROM RESTART FILE WITH almost COMPLETE DATA    ISTART=7
C****
      CASE (7)             ! converted model II' (B399) format (no snow)
        call io_rsf(iu_AIC,IhrX,irsfic,ioerr)
        CLOSE (iu_AIC)
        if (ioerr.eq.1) goto 800
        iniSNOW = .TRUE.      ! extract snow data from first soil layer
        iniOCEAN = .TRUE. ! read in ocean ic
C****
C****   Data from current type of RESTART FILE           ISTART=8
C****
      CASE (8)  ! no need to read SRHR,TRHR,FSF,TSFREZ,diag.arrays
        call io_rsf(iu_AIC,IhrX,irsfic,ioerr)
        CLOSE (iu_AIC)
        if (ioerr.eq.1) goto 800
        iniOCEAN = .TRUE. ! read in ocean ic
      END SELECT
C**** Check consistency of starting time
      IF (ISTART.ge.3.and.(MOD(IHRI-IHRX,8760).ne.0)) THEN
        WRITE (6,*) ' Difference in hours between ',
     *       'Starting date and Data date:',MOD(IHRI-IHRX,8760)
        WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
        STOP 'INPUT: start date inconsistent with data'
      END IF
C**** Set flag to initialise lake variables if they are not in I.C.
      IF (ISTART.lt.8) inilake=.TRUE.
C****
C**** Use IRAND<0 to perturb initial temperatures to create ensembles
C****                              perturbation is at most 1 degree C
      IF (IRAND.LT.0) THEN
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
      WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
     *  '0Model started on',datei,aMONTH(monthi),yeari,' Hr',houri,
     *  'ISTART =',ISTART,XLABEL(1:80)

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
        call getunit("AIC",iu_AIC,.true.,.true.)
        call io_rsf(iu_AIC,ItimeX,irerun,ioerr)
        if (ioerr.eq.1) goto 800
        WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *    ISTART,', HOUR=',ItimeX,' ',XLABEL(1:80)
        TIMING = 0
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
      call io_rsf(KDISK0,ItimeX,ioread,ioerr)
      if (ioerr.eq.1) then    ! try the other restart file
         rewind kdisk0
         KDISK=3-KDISK0
         WRITE (6,'(A,I1,A,I1)')
     *     ' Read Error on fort.',kdisk0,' trying fort.',kdisk
         KDISK0=KDISK
         call io_rsf(KDISK0,ItimeX,ioread,ioerr)
         if (ioerr.eq.1) go to 850
         IF (istart.eq.10) KDISK0=3-KDISK
      end if
      KDISK=KDISK0
      IF (istart.gt.10) KDISK=3-KDISK

      WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
     *   KDISK,', HOUR=',ItimeX,' ',XLABEL(1:80)
  500 CONTINUE
C**** UPDATE C ARRAY FROM INPUTZ
      REWIND 8
      READ (8,NML=INPUTZ)
      REWIND 8
C**** For documentation purposes only, find PLTOP (appears on printout)
      DO L=1,LM
      PLTOP(L)=PTOP+PSFMPT*SIGE(L+1)
      END DO

C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************
  600 CONTINUE
      IF (KEYCT.LE.1) KEYNR=0
      IF (KEYCT.LE.1) KEYCT=1
C****
C**** Update ItimeE only if YearE or IhourE is specified in the rundeck
C****
      IF (yearE.ge.0) ItimeE = (( (yearE-Iyear0)*JDperY +
     *    JDendofM(monthE-1)+dateE-1 )*HR_IN_DAY + HourE )*NDAY/24
C**** Alternate (old) way of specifying end time
      if(IHOURE.gt.0) ItimeE=IHOURE*NDAY/24

      IF(ISTART.LT.0) THEN
        call reset_diag(1)
        monacc = 0
        do k=1,iargc()
          call getarg(k,filenm)
          call getunit(filenm,iu_AIC,.true.,.true.)
          call io_rsf(iu_AIC,itime,ioread_single,ioerr)
          write(6,*) 'read: ',filenm(1:70)
          call closeunits
        end do
        months=0 ; years=monacc(jmon0) ; mswitch=0
        do k=1,12
          if (monacc(k).eq.years) then
            months=months+1
          else if (monacc(k).ne.0) then
            write(6,*) 'uneven period:',monacc
            stop 'uneven period'
          end if
          if(k.gt.1.and.monacc(k).ne.monacc(k-1)) mswitch = mswitch+1
        end do
        if (mswitch.gt.2) then
          write(6,*) 'non-consecutive period:',monacc
          stop 'non-consecutive period'
        end if
        call aPERIOD (JMON0,JYEAR0,months,years,acc_period,Ldate)
        LLAB1 = INDEX(XLABEL(1:17),'(') -1
        IF (LLAB1.LT.1) LLAB1=16
        if (index(XLABEL(1:LLAB1),' ').gt.0)
     *    LLAB1=index(XLABEL(1:LLAB1),' ')-1
        if (iargc().gt.1) then    ! save the summed acc-file
          write(6,*) iargc(),' files are summed up'
          keyct=1 ; KEYNR=0
          XLABEL(128:132)='     '
          XLABEL(120:132)=acc_period(1:3)//' '//acc_period(4:Ldate)
          OPEN (30,FILE=acc_period(1:Ldate)//'.acc'//XLABEL(1:LLAB1),
     *         FORM='UNFORMATTED')
          call io_rsf (30,Itime,iowrite_single,ioerr)
          CLOSE (30)
        end if
        if(qcheck) then          ! open the giss output files
          call getunit(acc_period(1:Ldate)//'.ij'//XLABEL(1:LLAB1),
     *                 iu_ij,.true.,.false.)
          call getunit(acc_period(1:Ldate)//'.jk'//XLABEL(1:LLAB1),
     *                 iu_jl,.true.,.false.)
          call getunit(acc_period(1:Ldate)//'.il'//XLABEL(1:LLAB1),
     *                 iu_il,.true.,.false.)
        end if
        ItimeE = -1
        close (6)
        open(6,file=acc_period(1:Ldate)//'.'//XLABEL(1:LLAB1)//'.PRT',
     *     FORM='FORMATTED')
      END IF
C****
C**** Recompute dtsrc,dt making NIdyn=dtsrc/dt(dyn) a multiple of 2
C****
      DTsrc = SDAY/NDAY   ! currently 1 hour
      NIdyn = 2*nint(.5*dtsrc/dt)
      DT = DTsrc/NIdyn
C**** Restrict NMONAV to 1(default),2,3,4,6,12, i.e. a factor of 12
      if (NMONAV.lt. 1) NMONAV=1
      if (NMONAV.gt.12) NMONAV=12
      NMONAV = 12/nint(12./nmonav)
      write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'
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
      call getunit("REG",iu_REG,.true.,.true.)
      READ(iu_REG) TITREG,JREG,NAMREG
      WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
C**** Note that FLAKE0 is read in only to provide initial values
C**** Actual array is set from restart file.
      call getunit("TOPO",iu_TOPO,.true.,.true.)

      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      CALL READT (iu_TOPO,0,FLAKE0,IM*JM,FLAKE0,1) ! Orig. Lake fraction
      CALL READT (iu_TOPO,0,FEARTH,IM*JM,FEARTH,1) ! Earth frac. (no LI)
      CALL READT (iu_TOPO,0,FLICE,IM*JM,FLICE,1)   ! Land ice fraction
C****
      CALL READT (iu_TOPO,0,ZATMO,IM*JM,ZATMO,1)   ! Topography
      ZATMO = ZATMO*GRAV                           ! Geopotential
      CALL READT (iu_TOPO,0,HLAKE,IM*JM,HLAKE,2)   ! Lake Depths
      CLOSE (iu_TOPO)

C**** Initialise some modules before finalising Land/Ocean/Lake/LI mask
C**** Initialize ice
      CALL init_ice              !???? needed if istart<0 ????
C**** Initialise lake variables (including river directions)
      CALL init_LAKES(inilake)
C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN(iniOCEAN)  !???? needed if istart<0 ????
C**** Initialize land ice (must come after oceans)
      CALL init_LI               !???? needed if istart<0 ????

C**** Make sure that constraints are satisfied by defining FLAND/FEARTH
C**** as residual terms. (deals with SP=>DP problem)
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          FLAND(I,J)=1.-FOCEAN(I,J) ! Land fraction
          IF (FLAKE(I,J).gt.0) THEN
            WRITE(6,*) "Ocean and lake cannot co-exist in same grid box"
     *           ,i,j,FOCEAN(I,J),FLAKE(I,J)
            FLAKE(I,J)=0
          END IF
        ELSEIF (FLAKE(I,J).gt.0) THEN
          FLAND(I,J)=1.-FLAKE(I,J)
        ELSE
          FLAND(I,J)=1.
        END IF
        FEARTH(I,J)=FLAND(I,J)-FLICE(I,J) ! Earth fraction
      END DO
      END DO
C**** set land components of FTYPE array. Summation is necessary for
C**** cases where Earth and Land Ice are lumped together
      FTYPE(ITLANDI,:,:)=0.
      FTYPE(ITEARTH,:,:)=FEARTH
      FTYPE(ITLANDI,:,:)=FTYPE(ITLANDI,:,:)+FLICE

C**** READ IN VEGETATION DATA SET: VDATA
c     if(istart.gt.0) then
        call getunit("VEG",iu_VEG,.true.,.true.)
        DO K=1,11
          CALL READT (iu_VEG,0,VDATA(1,1,K),IM*JM,VDATA(1,1,K),1)
        END DO
        CLOSE (iu_VEG)
c     end if
C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS
C**** Recompute GHDATA if redoGH (new soils data)
C****
      CALL init_GH(DTsrc/NIsurf,redoGH,iniSNOW)
C**** Initialize pbl (and read in file containing roughness length data)
      if(istart.gt.0) CALL init_pbl(iniPBL)
C****
      if(istart.gt.0) CALL RINIT (IRAND)
      CALL FFT0 (IM)
      if(istart.gt.0) CALL init_CLD
      CALL init_DIAG
      if(istart.gt.0) CALL init_QUS(im,jm,lm)
      IF (KDIAG(2).EQ.9.AND.SKIPSE.EQ.0..AND.KDIAG(3).LT.9) KDIAG(2)=8
      if(istart.gt.0) WRITE (6,INPUTZ)
      WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
      WRITE (6,'(A14,2I8)') "KACC0,KTACC0=",KACC0,KTACC0
      WRITE (6,'(A14,3I4)') "IM,JM,LM=",IM,JM,LM
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  800 WRITE (6,'(A,I4/" ",A)')
     *  '0ERROR ENCOUNTERED READING AIC ISTART=', ISTART,XLABEL(1:80)
      STOP 'INPUT: READ ERROR FOR AIC'
  830 WRITE(6,*) 'READ ERROR FOR GIC'
      STOP 'INPUT: READ ERROR FOR GIC'
  850 WRITE (6,'(A)')
     *  '0ERRORS ON BOTH RESTART DATA SETS. TERMINATE THIS JOB'
      STOP 'INPUT: ERRORS ON BOTH RESTART FILES'
  890 WRITE (6,'(A,I5)') '0INCORRECT VALUE OF ISTART',ISTART
      STOP 'INPUT: ISTART-SPECIFICATION INVALID'
      END SUBROUTINE INPUT

      SUBROUTINE DAILY(IEND)
!@sum  DAILY performs daily model-related tasks and at start
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : orbit
      USE MODEL_COM, only : im,jm,p,itime,itimei,ptop,psf,ls1,jday
     *     ,iyear0,nday,jdpery,jyear,jmon,jdendofm,jdate,aMON,aMONTH
      USE GEOM, only : areag,dxyp
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

      IF (ABS(DELTAP).gt.1d-6)
     *     WRITE (6,'(A25,F10.6/)') '0PRESSURE ADDED IN GMP IS',DELTAP
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
      aMON=aMONTH(JMON)

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
      CALL ORBIT (OBLIQ,ECCN,OMEGT,DFLOAT(JDAY)-.5,RSDIST,SIND,COSD,LAM)

      RETURN
      END SUBROUTINE DAILY

      SUBROUTINE CHECKT (SUBR)
!@sum  CHECKT Checks arrays for NaN/INF and reasonablness
!@auth Original Development Team
!@ver  1.0

C**** CHECKT IS TURNED ON BY SETTING QCHECK=.TRUE. IN NAMELIST
C**** REMEMBER TO SET QCHECK BACK TO .FALSE. AFTER THE ERRORS ARE
C**** CORRECTED.

      USE MODEL_COM
      USE DAGCOM, only : QCHECK
      IMPLICIT NONE
      INTEGER I,J
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

        DO J=1,JM
        DO I=1,IM
          IF (Q(I,J,1).gt.1d-1)print*,SUBR," Q BIG ",i,j,Q(I,J,1:LS1)
          IF (T(I,J,1).gt.50.) print*,SUBR," T BIG ",i,j,T(I,J,1:LS1)
        END DO
        END DO

C**** Check PBL arrays
         CALL CHECKPBL(SUBR)
C**** Check Ocean arrays
         CALL CHECKO(SUBR)
C**** Check Ice arrays
         CALL CHECKI(SUBR)
C**** Check Lake arrays
         CALL CHECKL(SUBR)
C**** Check Earth arrays
c        CALL CHECKE(SUBR)
      END IF
      RETURN
      END SUBROUTINE CHECKT

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
      LOGICAL QCHECK3
      INTEGER I,J,L !@var I,J,L loop variables

      DO L=1,LN
      DO J=1,JN
      DO I=1,IN
        IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.)) THEN
          WRITE (6,*) FIELD,': ',I,J,L,A(I,J,L),'after ',SUBR
          IF (J.LT.JN.AND.J.GT.1) QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
      IF (QCHECK3) STOP 'CHECK3'
      RETURN
      END SUBROUTINE CHECK3

      SUBROUTINE io_rsf(kunit,it,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
!@calls io_model,io_ocean,io_lakes,io_seaice,io_earth,io_soils,io_snow
!@calls io_landice,io_bldat,io_pbl,io_clouds,io_somtq,io_rad,io_diags
!@calls io_ocdiag
      USE MODEL_COM, only : ioread_single,iowrite_single

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
      call io_label  (kunit,it,iaction,ioerr)
      it1=it
      if(iaction.ne.ioread_single.and.iaction.ne.iowrite_single) then
        call io_model  (kunit,iaction,ioerr)
        call io_ocean  (kunit,iaction,ioerr)
        call io_lakes  (kunit,iaction,ioerr)
        call io_seaice (kunit,iaction,ioerr)
        call io_earth  (kunit,iaction,ioerr)
        call io_soils  (kunit,iaction,ioerr)
        call io_snow   (kunit,iaction,ioerr)
        call io_landice(kunit,iaction,ioerr)
        call io_bldat  (kunit,iaction,ioerr)
        call io_pbl    (kunit,iaction,ioerr)
        call io_clouds (kunit,iaction,ioerr)
        call io_somtq  (kunit,iaction,ioerr)
        call io_rad    (kunit,iaction,ioerr)
      end if
      call io_diags  (kunit,it,iaction,ioerr)
      call io_ocdiag (kunit,it,iaction,ioerr)

      if (it1.ne.it) THEN
        WRITE(6,*) "TIMES DO NOT MATCH READING IN RSF FILE",it,it1
        ioerr=1
      END IF
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN RESTART FILE: KUNIT="
     *     ,kunit
      close (kunit)

      RETURN
      END SUBROUTINE io_rsf
