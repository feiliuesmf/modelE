#include "rundeck_opts.h"

      PROGRAM GISS_modelE
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  1.0 (Based originally on B399)
      USE MODEL_COM
      USE RANDOM
      USE DAGCOM, only : keyct,keynr,kdiag,oa,monacc
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE PARAM
      USE SOIL_DRV, only: daily_earth, ground_e
      IMPLICIT NONE

      INTEGER I,J,L,K,M,KSS6,MSTART,MNOW,MODD5D,months,ioerr,Ldate
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
      IF (Itime.GE.ItimeE) then     ! includes ISTART<1 case
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
     *   'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME
         CALL TIMER (MNOW,MELSE)
C****
C**** Open and position output history files if needed
C****
      if (Kvflxo.ne.0) then
        write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
        call openunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
        call IO_pos(iu_VFLXO,Itime,2*im*jm*12,Nday)
      end if
      if (Nslp.ne.0) then
        write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
        call openunit('SLP'//aDATE(1:7),iu_SLP,.true.,.false.)
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
         call set_param( "IRAND", IRAND, 'o' )
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
            write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
            call closeunit( iu_VFLXO )
            call openunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
          end if
          if (Nslp.ne.0) then
            write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
            call closeunit( iu_SLP )
            call openunit('SLP'//aDATE(1:7),iu_SLP,.true.,.false.)
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

C**** calculate some dynamic variables for the PBL
      CALL PGRAD_PBL

         CALL CHECKT ('DYNAM ')
         CALL TIMER (MNOW,MDYN)

         IF (MODD5D.EQ.0) CALL DIAG5A (7,NIdyn)             
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
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)             
         IF (MODD5S.EQ.0) CALL DIAG9A (3)
C**** RADIATION, SOLAR AND THERMAL
      CALL RADIA
         CALL CHECKT ('RADIA ')
         CALL TIMER (MNOW,MRAD)
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)            
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
C**** IF ATURB is used in rundeck then this is a dummy call
C**** CALCULATE DRY CONVECTION ABOVE PBL
      CALL DIFFUS (2,LM-1,dtsrc)
         IF (MODD5S.EQ.0) CALL DIAG9A (7)
         CALL CHECKT ('DRYCNV')
         CALL TIMER (MNOW,MSURF)

C**** CALL OCEAN DYNAMIC ROUTINES
      CALL ODYNAM
         CALL CHECKT ('OCEAN ')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAG9A (9)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)            
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
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1
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
           CALL DIAG5A (16,NDAY*NIdyn)                        
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
      IF (MOD(Itime-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 E-history
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
         call set_param( "NIPRNT", NIPRNT, 'o' )
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
          call aPERIOD (JMON0,JYEAR0,months,1,aDATE(1:12),Ldate)
          WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),JYEAR
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.12) m = m-12
            monacc(m) = 1
          end do
          OPEN (30,FILE=aDATE(1:7)//'.acc'//XLABEL(1:LRUNID),
     *         FORM='UNFORMATTED')
          call io_rsf (30,Itime,iowrite_single,ioerr)
          CLOSE (30)
C**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
          IF (KCOPY.GT.1) THEN
            CALL RFINAL (IRAND)
            call set_param( "IRAND", IRAND, 'o' )
            OPEN(30,FORM='UNFORMATTED',
     *              FILE='1'//aDATE(8:14)//'.rsf'//XLABEL(1:LRUNID))
            call io_rsf(30,Itime,iowrite_mon,ioerr)
            CLOSE (30)
          END IF
C**** KCOPY > 2 : ALSO SAVE THE OCEAN DATA TO INITIALIZE DEEP OCEAN RUNS
          IF (KCOPY.GT.2) THEN
            OPEN (30,FILE=aDATE(1:7)//'.oda'//XLABEL(1:LRUNID),
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
      IF (LABSSW.ne.offssw) THEN
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
      call set_param( "IRAND", IRAND, 'o' )
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


      subroutine init_Model
!@sum This program reads most of parameters from the database (DB)
C**** get_param( "A", X ) reads parameter A into variable X
C**** if "A" is not in the database, it will generate an error
C**** message and stop
C**** sync_param( "B", Y ) reads parameter B into variable Y
C**** if "B" is not in the database, then Y is unchanged and its
C**** value is saved in the database as "B" (here sync = synchronize)
      USE MODEL_COM, only : LM,NAMD6,IJD6,NIPRNT,MFILTR,XCDLM,NDASF
     *     ,NDA4,NDA5S,NDA5K,NDA5D,NDAA,NFILTR,NRAD,Kvflxo,Nslp
     *     ,NMONAV,Ndisk,Nssw,KCOPY,KOCEAN,PSF,NIsurf,iyear1
     $     ,PTOP,LS1,IRAND,LSDRAG
     $     ,ItimeI,PSFMPT,PSTRAT,SIG,SIGE
      USE PARAM
      implicit none
      INTEGER L

C**** Rundeck parameters:
      call sync_param( "NAMD6", NAMD6, 4 )
      call sync_param( "IJD6", IJD6(1:8,1), 8)
      call sync_param( "NMONAV", NMONAV )
      call sync_param( "NIPRNT", NIPRNT )
      call sync_param( "MFILTR", MFILTR )
      call sync_param( "XCDLM", XCDLM, 2 )
      call sync_param( "NDASF", NDASF )
      call sync_param( "NDA4", NDA4 ) !!
      call sync_param( "NDA5S", NDA5S ) !!
      call sync_param( "NDA5K", NDA5K ) !!
      call sync_param( "NDA5D", NDA5D ) !!
      call sync_param( "NDAA", NDAA ) !!
      call sync_param( "NFILTR", NFILTR ) !!
      call sync_param( "NRAD", NRAD ) !!
      call sync_param( "Kvflxo", Kvflxo ) !!
      call sync_param( "Nslp", Nslp )
      call sync_param( "Ndisk", Ndisk )
      call sync_param( "Nssw", Nssw )
      call sync_param( "KCOPY", KCOPY )
      call sync_param( "KOCEAN", KOCEAN )
      call sync_param( "NIsurf", NIsurf )
      call sync_param( "IRAND", IRAND )

C**** Non-Rundeck parameters

C**** Calculate level for application of SDRAG
C**** All levels above and including 1 mb, or LM for trop model
      DO L=1,LM
        IF (PTOP+PSFMPT*SIGE(L+1)+1d-5.lt.1d0 .and.
     *      PTOP+PSFMPT*SIGE(L)+1d-5.gt.1d0) EXIT
      END DO
      LSDRAG=MIN(L,LM)
      WRITE(6,*) "Level for SDRAG = ",LSDRAG
      RETURN
C****
      end subroutine init_Model


      SUBROUTINE INPUT
C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE CONSTANT, only : grav,kapa,sday,shi,lhm
      USE MODEL_COM, only : im,jm,lm,wm,u,v,t,p,q,fearth,fland
     *     ,focean,flake0,flice,hlake,zatmo,sig,dsig,sige
     *     ,bydsig,xlabel,lrunid,nmonav
     *     ,irand,psf,ptop
     *     ,nisurf,nidyn,nday,dt,dtsrc,kdisk,jmon0,jyear0
     *     ,iyear1,itime,itimei,itimee
     *     ,ls1,psfmpt,pstrat,idacc,jyear,jmon,jday,jdate,jhour
     *     ,vdata,aMONTH,jdendofm,jdpery,aMON,aMON0,ioread,irerun
     *     ,ioread_single,irsfic,iowrite_single,ftype,itearth,itlandi
     *     ,mdyn,mcnds,mrad,msurf,mdiag,melse,Itime0,Jdate0,Jhour0
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : geom_b
      USE RANDOM
      USE RADNCB, only : rqt,lm_req
      USE CLOUDS_COM, only : ttold,qtold,svlhx,rhsav,cldsav
      USE PBLCOM
     &     , only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
      USE DAGCOM, only : acc_period,monacc,kacc,tsfrez,kdiag,jreg
     &  ,titreg,namreg,hr_in_day,iwrite,jwrite,itwrite,qcheck,oa
     &  ,iu_ij,iu_jl,iu_il,iu_j
      USE LAKES_COM, only : flake
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : timing,ntimeacc
      USE PARAM
      USE PARSER
      USE SOIL_DRV, only: init_gh


      IMPLICIT NONE
!@var iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_VEG

      INTEGER I,J,L,K,KLAST,KDISK0,ITYPE,IM1
     *     ,IR,IREC,NOFF,ioerr
!@nlparam HOURI,DATEI,MONTHI,YEARI        start of model run
!@nlparam HOURE,DATEE,MONTHE,YEARE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
     *             HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHOURE=-1,
!@nlparam ISTART  postprocessing(-1)/start(1-8)/restart(>8)  option
!@nlparam IRANDI  random number seed to perturb init.state (if>0)
     *             ISTART=10, IRANDI=0
      REAL*8 TIJL,CDM,TEMP,X
      REAL*4 XX4
      INTEGER Itime1,Itime2,ItimeX,IhrX,iargc

      LOGICAL :: redoGH = .FALSE.,iniPBL = .FALSE., inilake = .FALSE.,
     &           iniSNOW = .FALSE.  ! true = restart from "no snow" rsf
     &           ,iniOCEAN = .FALSE.

      CHARACTER NLREC*80,filenm*100,RLABEL*132
      NAMELIST/INPUTZ/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,KDIAG
     *     ,IHOURE, HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI

C****
C**** default settings for prog. variables etc
C****
      TEMP=250.
      TSAVG(:,:)=TEMP
      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=TEMP  ! will be changed to pot.temp later
      Q(:,:,:)=3.D-6
      P(:,:)=PSFMPT
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
C**** Print Header and Label (2 lines) from rundeck
C****
      WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(5,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label
C****
C**** Print preprocessing options (if any are defined)
C****
#ifdef TRACERS_ON
      write(6,*) 'This program includes tracers code'
#endif
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
C****
C**** Read parameters from the rundeck to the database
C****
      call parse_params( 8 )
C**** Get those parameters which are needed in this subroutine
      if(is_set_param("DTsrc"))  call get_param( "DTsrc", DTsrc )
      if(is_set_param("DT"))     call get_param( "DT", DT )
      if(is_set_param("NIsurf")) call get_param( "NIsurf", NIsurf ) !
      if(is_set_param("IRAND"))  call get_param( "IRAND", IRAND )
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      call reset_diag(1)

C***********************************************************************
C****                                                               ****
C****        Post-process one or more ACC-files : ISTART < 1        ****
C****                                                               ****
C***********************************************************************
      if (istart.le.0) then
        monacc = 0
        do k=1,iargc()
          call getarg(k,filenm)
          call openunit(filenm,iu_AIC,.true.,.true.)
          call io_rsf(iu_AIC,itime,ioread_single,ioerr)
          call closeunit(iu_AIC)
        end do
        GO TO 500 
      end if

      if (istart.ge.9) go to 400
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
      IF (ISTART.gt.1) call openunit("AIC",iu_AIC,.true.,.true.)
C****
C**** Set quantities that are derived from the namelist parameters
C****
!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc later if necessary
      NDAY = 2*NINT(.5*SDAY/DTsrc)

C**** Get Start Time; at least YearI HAS to be specified in the rundeck
      IF (YearI.lt.0) then
        WRITE(6,*) 'Please choose a proper start year yearI, not',yearI
        STOP 'INPUT: yearI not provided'
      END IF
      IF (Iyear1.lt.0) Iyear1 = yearI
      IhrI = HourI +
     +  HR_IN_DAY*(dateI-1 + JDendofM(monthI-1) + JDperY*(yearI-Iyear1))
      ITimeI = IhrI*NDAY/24  !  internal clock counts DTsrc-steps
      Itime=ItimeI
      IF (IhrI.lt.0) then
        WRITE(6,*) 'Improper start time OR Iyear1=',Iyear1,' > yearI;',
     *     ' yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
        STOP 'INPUT: Improper start date or base year Iyear1'
      END IF
C**** Check the vertical layering defined in RES_ (is sige(ls1)=0 ?)
      IF (SIGE(LS1).ne.0.) then
        write(6,*) 'bad vertical layering: ls1,sige(ls1)',ls1,sige(ls1)
        STOP 'INPUT: ls1 incorrectly set in RES_'
      END IF
C****
C**** Get Ground conditions from a separate file - ISTART=1,2
C****
      IF (ISTART.LE.2) THEN
C**** Set flag to initialise pbl and snow variables
        iniPBL=.TRUE.
        iniSNOW = .TRUE.  ! extract snow data from first soil layer
        iniOCEAN = .TRUE. ! read in ocean ic
C**** Read in ground initial conditions
        call openunit("GIC",iu_GIC,.true.,.true.)
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
        call closeunit (iu_GIC)
      END IF
C****
C**** Get primary Atmospheric data from NMC tapes - ISTART=2
C****
      IF (ISTART.EQ.2) THEN
C**** Use title of first record to get the date and make sure  ???
C**** it is consistent with IHRI (at least equal mod 8760)     ???
C****            not yet implemented but could easily be done  ???
        XLABEL(1:80)='Observed atmospheric data from NMC tape'
Csoon   READ (iu_AIC) XLABEL(1:80)
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
        if (ioerr.eq.1) goto 800
        iniSNOW = .TRUE.      ! extract snow data from first soil layer
        iniOCEAN = .TRUE. ! read in ocean ic
C****
C****   Data from current type of RESTART FILE           ISTART=8
C****
      CASE (8)  ! no need to read SRHR,TRHR,FSF,TSFREZ,diag.arrays
        call io_rsf(iu_AIC,IhrX,irsfic,ioerr)
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
!**** IRANDI seed for random perturbation of initial conditions (if/=0)
C****        perturbation is at most 1 degree C
      IF (IRANDI.NE.0) THEN
        CALL RINIT (IRANDI)
        DO L=1,LM
        DO J=1,JM
        DO I=1,IM
           TIJL=T(I,J,L)*(P(I,J)*SIG(L)+PTOP)**KAPA-1.+2*RANDU(X)
           T(I,J,L)=TIJL/(P(I,J)*SIG(L)+PTOP)**KAPA
        END DO
        END DO
        END DO
        WRITE(6,*) 'Initial conditions were perturbed !!',IRANDI
      END IF
C**** Close "AIC" here if it was opened
      IF (ISTART.gt.1) call closeunit(iu_AIC)
C**** Sending parameters which had just been set to the DB
      ! the following lines overwrite rundeck parameters
      ! the following are NON-rundeck parameters

      WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
     *  '0Model started on',datei,aMONTH(monthi),yeari,' Hr',houri,
     *  'ISTART =',ISTART,XLABEL(1:80)    ! report input file label
      XLABEL = RLABEL                     ! switch to rundeck label

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
  400 SELECT CASE (ISTART)
C****
C****   DATA FROM end-of-month RESTART FILE     ISTART=9
C****        mainly used for REPEATS and delayed EXTENSIONS
      CASE (9)    ! no need to read diag.arrays
        call openunit("AIC",iu_AIC,.true.,.true.)
        call io_rsf(iu_AIC,Itime,irerun,ioerr)
        call closeunit(iu_AIC)
        if (ioerr.eq.1) goto 800
        WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *    ISTART,', TIME=',Itime,' ',XLABEL(1:80) ! sho input file label
        XLABEL = RLABEL                        ! switch to rundeck label
C****
!**** IRANDI seed for random perturbation of current state (if/=0)
C****        perturbation is at most 1 degree C
        IF (IRANDI.NE.0) THEN
          CALL RINIT (IRANDI)
          DO L=1,LM
          DO J=1,JM
          DO I=1,IM
             TIJL=T(I,J,L)*(P(I,J)*SIG(L)+PTOP)**KAPA-1.+2*RANDU(X)
             T(I,J,L)=TIJL/(P(I,J)*SIG(L)+PTOP)**KAPA
          END DO
          END DO
          END DO
          WRITE(6,*) 'Current temperatures were perturbed !!',IRANDI
        END IF
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
         IF (Itime1+Itime2.LE.-2.) GO TO 850
                               KDISK=1
         IF (Itime2.GT.Itime1) KDISK=2
         IF (ISTART.GE.13)     KDISK=3-KDISK
      CASE (11,12)
                               KDISK=ISTART-10
      END SELECT
  430 call io_rsf(KDISK,Itime,ioread,ioerr)
      if (ioerr.eq.1) then
         if (istart.gt.10) go to 850  ! no 2nd chance if istart/=10
         KDISK=3-KDISK                ! try the earlier restart file
         WRITE (6,'(A,I1,A,I1)')
     *     ' Read Error on fort.',3-kdisk,' trying fort.',kdisk
         ISTART=110
         go to 430
      end if
      WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
     *   KDISK,', Time=',Itime,' ',XLABEL(1:80)

C**** Switch KDISK if the other file is (or may be) bad (istart>10)
C****     so both files will be fine after the next write execution
      IF (istart.gt.10) KDISK=3-KDISK
C**** Keep KDISK after reading from the later restart file, so that
C****     the same file is overwritten first; in case of trouble,
C****     the earlier restart file will still be available

  500 CONTINUE
C**** Get parameters we just read from rsf file. Only those
C**** parameters which we need in "INPUT" should be extracted here.
      if(is_set_param("DTsrc"))  call get_param( "DTsrc", DTsrc )
      if(is_set_param("DT"))     call get_param( "DT", DT )
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )

C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************
  600 CONTINUE

C**** initialize Lrunid (length of the identifying part of XLABEL)
C****
      LRUNID = INDEX(XLABEL(1:16),'(') -1
      IF (LRUNID.LT.1) LRUNID=16
      if (index(XLABEL(1:LRUNID),' ').gt.0)
     *     LRUNID=index(XLABEL(1:LRUNID),' ')-1

C**** Update ItimeE only if YearE or IhourE is specified in the rundeck
C****
      IF (yearE.ge.0) ItimeE = (( (yearE-iyear1)*JDperY +
     *    JDendofM(monthE-1)+dateE-1 )*HR_IN_DAY + HourE )*NDAY/24
C**** Alternate (old) way of specifying end time
      if(IHOURE.gt.0) ItimeE=IHOURE*NDAY/24

C**** Recompute dtsrc,dt making NIdyn=dtsrc/dt(dyn) a multiple of 2
C****
      if (is_set_param("DTsrc") .and. nint(sday/DTsrc).ne.NDAY) then
        write(6,*) 'DTsrc=',DTsrc,' has to stay at/be set to',SDAY/NDAY
        stop 'INPUT: DTsrc inappropriately set'
      end if
      DTsrc = SDAY/NDAY   ! currently 1 hour
      call set_param( "DTsrc", DTsrc, 'o' )   ! copy DTsrc into DB

      NIdyn = 2*nint(.5*dtsrc/dt)
      if (is_set_param("DT") .and. nint(DTsrc/dt).ne.NIdyn) then
        write(6,*) 'DT=',DT,' has to be changed to',DTsrc/NIdyn
        stop 'INPUT: DT inappropriately set'
      end if
      DT = DTsrc/NIdyn
      call set_param( "DT", DT, 'o' )         ! copy DT into DB

C**** NMONAV has to be 1(default),2,3,4,6,12, i.e. a factor of 12
      if (NMONAV.lt.1 .or. MOD(12,NMONAV).ne.0) then
        write (6,*) 'NMONAV has to be 1,2,3,4,6 or 12, not',NMONAV
        stop 'INPUT: nmonav inappropriately set'
      end if
      write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'

C**** Updating Parameters: If any of them changed beyond this line
C**** use set_param(.., .., 'o') to update them in the database (DB)

C**** Get the rest of parameters from DB or put defaults to DB
      call init_Model

C**** Set julian date information
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      call getdte(Itime0,Nday,iyear1,Jyear0,Jmon0,J,Jdate0,Jhour0,amon0)

C****
C**** COMPUTE GRID RELATED VARIABLES AND READ IN TIME-INDEPENDENT ARRAYS
C****
C**** CALCULATE SPHERICAL GEOMETRY
      CALL GEOM_B
      CALL CALC_AMPK(LM)

C**** READ SPECIAL REGIONS FROM UNIT 29
      call openunit("REG",iu_REG,.true.,.true.)
      READ(iu_REG) TITREG,JREG,NAMREG
      WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG
      call closeunit(iu_REG)

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
C**** Note that FLAKE0 is read in only to provide initial values
C**** Actual array is set from restart file.
      call openunit("TOPO",iu_TOPO,.true.,.true.)

      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      CALL READT (iu_TOPO,0,FLAKE0,IM*JM,FLAKE0,1) ! Orig. Lake fraction
      CALL READT (iu_TOPO,0,FEARTH,IM*JM,FEARTH,1) ! Earth frac. (no LI)
      CALL READT (iu_TOPO,0,FLICE,IM*JM,FLICE,1)   ! Land ice fraction
C****
      CALL READT (iu_TOPO,0,ZATMO,IM*JM,ZATMO,1)   ! Topography
      ZATMO = ZATMO*GRAV                           ! Geopotential
      CALL READT (iu_TOPO,0,HLAKE,IM*JM,HLAKE,2)   ! Lake Depths
      call closeunit(iu_TOPO)

C**** Initialise some modules before finalising Land/Ocean/Lake/LI mask
C**** Initialize ice
      CALL init_ice        
C**** Initialise lake variables (including river directions)
      CALL init_LAKES(inilake)
C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN(iniOCEAN)
C**** Initialize land ice (must come after oceans)
      CALL init_LI             

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
        call openunit("VEG",iu_VEG,.true.,.true.)
        DO K=1,11
          CALL READT (iu_VEG,0,VDATA(1,1,K),IM*JM,VDATA(1,1,K),1)
        END DO
        call closeunit(iu_VEG)
c     end if
C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS
C**** Recompute Ground hydrology data if redoGH (new soils data)
C****
      CALL init_GH(DTsrc/NIsurf,redoGH,iniSNOW)
C**** Initialize pbl (and read in file containing roughness length data)
      if(istart.gt.0) CALL init_pbl(iniPBL)
C****
      if(istart.gt.0) CALL RINIT (IRAND)
      CALL FFT0 (IM)
      if(istart.gt.0) CALL init_CLD
      CALL init_DIAG(ISTART)
      if(istart.gt.0) CALL init_QUS(im,jm,lm)
      if(istart.gt.0) CALL init_MOM
      CALL init_RAD
      IF (KDIAG(2).EQ.9.AND.KDIAG(3).LT.9) KDIAG(2)=8
      WRITE (6,INPUTZ)
      call print_param( 6 )
      WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
      WRITE (6,'(A14,2I8)') "KACC=",KACC
      WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
      WRITE (6,*) "PLbot=",PTOP+PSFMPT*SIGE
C****
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
     *     ,iyear1,nday,jdpery,jyear,jmon,jdendofm,jdate,aMON,aMONTH
     *     ,jhour
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
  200 call getdte(Itime,Nday,iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

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
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
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
        call io_strat  (kunit,iaction,ioerr)
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
