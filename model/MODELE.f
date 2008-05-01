#include "rundeck_opts.h"
CAOO   Just to test CVS
      PROGRAM GISS_modelE
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  1.0 (Based originally on B399)
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE PARAM
      USE MODEL_COM
      USE DOMAIN_DECOMP, ONLY : init_app,grid,AM_I_ROOT,pack_data
      USE DOMAIN_DECOMP, ONLY : ESMF_BCAST
      USE DYNAMICS
      USE RAD_COM, only : dimrad_sv
      USE RANDOM
      USE GETTIME_MOD
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: mtrace
#endif
      USE DIAG_COM, only : oa,monacc,koa,acc_period
      USE SOIL_DRV, only: daily_earth, ground_e
      USE SUBDAILY, only : nsubdd,init_subdd,get_subdd,reset_subdd
      USE DIAG_SERIAL, only : print_diags
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: fv_core
      USE FV_INTERFACE_MOD, only: Initialize
      USE FV_INTERFACE_MOD, only: Run
      USE FV_INTERFACE_MOD, only: Checkpoint
      USE FV_INTERFACE_MOD, only: Finalize
      USE FV_INTERFACE_MOD, only: Compute_Tendencies
      USE FV_INTERFACE_MOD, only: init_app_clock
c$$$      USE MODEL_COM, only: clock
      USE ESMF_MOD, only: ESMF_Clock
      USE ESMF_CUSTOM_MOD, Only: vm => modelE_vm
#endif
      USE ATMDYN, only : DYNAM,CALC_TROP,PGRAD_PBL,SDRAG
     &     ,DISSIP,FILTER,CALC_AMPK, COMPUTE_DYNAM_AIJ_DIAGNOSTICS
     &     ,COMPUTE_WSAVE, getTotalEnergy, addEnergyAsDiffuseHeat
#ifdef TRACERS_ON
     &     ,trdynam
#endif
      USE ATMDYN_QDYNAM, only : QDYNAM
#ifdef SCM
      USE SCMCOM , only : SG_CONV,SCM_SAVE_T,SCM_SAVE_Q,SCM_DEL_T,
     &    SCM_DEL_Q,iu_scm_prt,iu_scm_diag
#endif

      !use soil_drv, only : conserv_wtg, conserv_htg
      IMPLICIT NONE

      INTEGER K,M,MSTART,MNOW,MODD5D,months,ioerr,Ldate,istart
      INTEGER iu_VFLXO,iu_ACC,iu_RSF,iu_ODA
      INTEGER :: MDUM = 0


      REAL*8, DIMENSION(NTIMEMAX) :: PERCENT
      REAL*8 DTIME,TOTALT , oa_glob(im,jm,koa)

      CHARACTER aDATE*14
      CHARACTER*8 :: flg_go='___GO___'      ! green light
      integer :: iflag=1
      external sig_stop_model

C**** Command line options
      LOGICAL :: qcrestart=.false.
      CHARACTER*32 :: ifile
      real :: lat_min=-90.,lat_max=90.,longt_min=0.,longt_max=360.
      integer :: tloopbegin, tloopend
#ifdef USE_FVCORE
      Character(Len=*), Parameter :: fv_config = 'fv_config.rc'
      Type (FV_CORE) :: fv
      Type (ESMF_CLOCK) :: clock
      character(len=28) :: fv_fname, fv_dfname
#endif
      integer :: L
      real*8 :: initialTotalEnergy, finalTotalEnergy

#ifdef SCM
c     Hard Code J - cannot get syncparam to work here ???? 
      J_TARG = 39
      call init_app(grid,im,jm,lm,J_TARG)
#else
C****
C****
      call init_app(grid,im,jm,lm)
#endif

      call alloc_drv()
C****
C**** Processing command line options
C****
      call read_options( qcrestart, ifile )
      if ( qcrestart ) then
        call print_restart_info
        call stop_model("Terminated normally: printed restart info",-1)
      endif
C****
C**** INITIALIZATIONS
C****
         CALL TIMER (MNOW,MDUM)


#ifdef USE_FVCORE
      CALL INPUT (istart,ifile,clock)
#else
      CALL INPUT (istart,ifile)
#endif

C**** Set run_status to "run in progress"
      if(istart > 0) call write_run_status("Run in progress...",1)

C****
C**** If run is already done, just produce diagnostic printout
C****
      IF (Itime.GE.ItimeE.and.Kradia.le.0) then ! includes ISTART<1 case
        if (ItimeE.gt.0) then
          months=(Jyear-Jyear0)*JMperY + JMON-JMON0
          call aPERIOD (JMON0,JYEAR0,months,1,0, acc_period,Ldate)
        end if
        call print_diags(0)
        if(istart < 1) then
          call stop_model ('Finished post-processing',-1)
        else
          call stop_model ('The run has already completed',13)
        end if
        ! no output files are affected
      END IF

      IF (AM_I_ROOT()) Then
         open(3,file='flagGoStop',form='FORMATTED',status='REPLACE')
         write (3,'(A8)') flg_go
         close (3)
      END IF
      call sys_signal( 15, sig_stop_model )  ! works only on single CPU
         MSTART=MNOW
         DO M=1,NTIMEACC
           MSTART= MSTART-TIMING(M)
         END DO
C**** INITIALIZE TIME PARAMETERS
      NSTEP=(Itime-ItimeI)*NIdyn
#ifdef SCM
      NSTEPSCM=ITIME-ITIMEI
#endif
         MODD5K=1000
      CALL DAILY(.false.)                  ! not end_of_day
      CALL daily_RAD(.false.)
      if (istart.le.9) call reset_diag(0)
      if(Kradia==10) call daily_OCEAN(.false.) ! to test OCLIM
      if (Kradia.le.0) then
        CALL daily_EARTH(.false.)          ! not end_of_day
        CALL daily_OCEAN(.false.)          ! not end_of_day
        CALL CALC_AMPK(LS1-1)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
        CALL daily_tracer(0)
#endif
           if (kradia.le.0) CALL CHECKT ('INPUT ')
      end if
      CALL UPDTYPE

C****
C**** Initialize FV dynamical core (ESMF component) if requested
C****
#ifdef USE_FVCORE
      Call Initialize(fv, istart, vm, grid%esmf_grid, clock,fv_config)
#endif

      if (AM_I_ROOT())
     *   WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)')
     *   '0NASA/GISS Climate Model (re)started',
     *   'Year',JYEAR,aMON,JDATE,', Hr',JHOUR,
     *   'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME
         CALL TIMER (MNOW,MELSE)
C****
C**** Open and position output history files if needed
C****
C**** Monthly files
      if (Kradia.ne.0 .and. Kradia<10) then
        write(aDATE(1:7),'(a3,I4.4)') aMON(1:3),Jyear
        if (Kradia.gt.0) aDATE(4:7)='    '
        call openunit(trim('RAD'//aDATE(1:7)),iu_RAD,.true.,.false.)
        if (Kradia.lt.0) call io_POS(iu_RAD,Itime-1,2*dimrad_sv,Nrad)
      end if
C**** Files for an accumulation period (1-12 months)
      write(aDATE(1:7),'(a3,I4.4)') aMON0(1:3),Jyear0
      if (Kvflxo.ne.0) then
        call openunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
        call io_POS(iu_VFLXO,Itime,2*im*jm*koa,Nday) ! real*8-dim -> 2*
      end if
C**** Initiallise file for sub-daily diagnostics, controlled by
C**** space-separated string segments in SUBDD & SUBDD1 in the rundeck
      call init_subdd(aDATE)
      call sys_flush(6)

C****
C**** MAIN LOOP
C****
      call gettime(tloopbegin)


      DO WHILE (Itime.lt.ItimeE)

c$$$         call test_save(__LINE__, itime)
C**** Every Ndisk Time Steps (DTsrc), starting with the first one,
C**** write restart information alternately onto 2 disk files
      IF (MOD(Itime-ItimeI,Ndisk).eq.0) THEN
         CALL RFINAL (IRAND)
         call set_param( "IRAND", IRAND, 'o' )
         iu_RSF=-1
         IF (AM_I_ROOT())
     *        call openunit(rsf_file_name(KDISK),iu_RSF,.true.,.false.)
         call io_rsf(iu_RSF,Itime,iowrite,ioerr)
         IF (AM_I_ROOT()) call closeunit(iu_RSF)
#ifdef USE_FVCORE
         fv_fname='fv.'   ; write(fv_fname(4:4),'(i1)') kdisk
         fv_dfname='dfv.' ; write(fv_dfname(5:5),'(i1)') kdisk
         call Checkpoint(fv, clock, fv_fname, fv_dfname)
#endif
         if (AM_I_ROOT())
     *        WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *     '0Restart file written on fort.',KDISK,'Year',
     *     JYEAR,aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
         KDISK=3-KDISK
         CALL TIMER (MNOW,MELSE)
      END IF
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY DAY
      IF (MOD(Itime,NDAY).eq.0) THEN
C**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
        if (kradia.le.0) call daily_DIAG
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY MONTH
        IF ( JDAY.eq.1+JDendOfM(Jmon-1) ) then
          write(aDATE(1:7),'(a3,I4.4)') aMON(1:3),Jyear
          if (Kradia.ne.0 .and. Kradia<10) then
            if (Kradia.gt.0) aDATE(4:7)='    '
            call closeunit( iu_RAD )
            call openunit(trim('RAD'//aDATE(1:7)),iu_RAD,.true.,.false.)
          end if
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY ACC.PERIOD
          months=(Jyear-Jyear0)*JMperY + JMON-JMON0
          if ( months.ge.NMONAV ) then
            call reset_DIAG(0)
            if (Kvflxo.ne.0) then
              call closeunit( iu_VFLXO )
              call openunit('VFLXO'//aDATE(1:7),iu_VFLXO,.true.,.false.)
            end if
C**** reset sub-daily diag files
            call reset_subdd(aDATE)
          end if   !  beginning of acc.period
        END IF     !  beginning of month
      END IF       !  beginning of day
C****
C**** INTEGRATE DYNAMIC TERMS (DIAGA AND DIAGB ARE CALLED FROM DYNAM)
C****
      if(Kradia>9) go to 100 ! to test daily/monthly procedures fast
      CALL CHECKT ('DYNAM0')
      if (kradia.le.0) then                   ! full model,kradia le 0
         MODD5D=MOD(Itime-ItimeI,NDA5D)
         IF (MODD5D.EQ.0) IDACC(7)=IDACC(7)+1
         IF (MODD5D.EQ.0) CALL DIAG5A (2,0)
         IF (MODD5D.EQ.0) CALL DIAGCA (1)


      PTOLD = P ! save for clouds
C**** Initialize pressure for mass fluxes used by tracers and Q
      PS (:,:)   = P(:,:)

C**** Initialise total energy (J/m^2)
      initialTotalEnergy = getTotalEnergy()
#ifdef SCM
      NSTEPSCM = ITIME-ITIMEI
      write(0,*) 'NSTEPSCM ',NSTEPSCM
      do L=1,LM
         SCM_SAVE_T(L) = T(I_TARG,J_TARG,L)
         SCM_SAVE_Q(L) = Q(I_TARG,J_TARG,L)
      enddo
c     do L=1,LM
c        write(iu_scm_prt,'(a13,i3,4(f9.3))')
c    &              'before dynam ',
c    &               L,T(I_TARG,J_TARG,L)*PK(L,I_TARG,J_TARG),
c    &               Q(I_TARG,J_TARG,L)*1000.0,
c    &               U(I_TARG,J_TARG,L),V(I_TARG,J_TARG,L)
c     enddo 
#endif


#ifndef USE_FVCORE
      CALL DYNAM()
#else
      ! Using FV instead
         IF (MOD(Itime-ItimeI,NDAA).eq.0) CALL DIAGA0
      call Run(fv, clock)

      CALL SDRAG (DTsrc)
         if (MOD(Itime-ItimeI,NDAA).eq.0) THEN
           call DIAGA
           call DIAGB
           call EPFLUX (U,V,T,P)
         endif

#endif
C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
C**** Currently energy is put in uniformly weighted by mass
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)
      call COMPUTE_DYNAM_AIJ_DIAGNOSTICS(PUA, PVA, DT)
#ifdef SCM
       do L=1,LM
          CONV(I_TARG,J_TARG,L) = SG_CONV(L)
       enddo
#endif
      SD_CLOUDS(:,:,:) = CONV(:,:,:)
      call COMPUTE_WSAVE(wsave, sda, T, PK, PEDN, NIdyn)
C**** Scale WM mixing ratios to conserve liquid water
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LS1-1
        WM(:,:,L)=WM(:,:,L)* (PTOLD/P)
      END DO
!$OMP  END PARALLEL DO

      CALL QDYNAM  ! Advection of Q by integrated fluxes
         CALL TIMER (MNOW,MDYN)
#ifdef TRACERS_ON
      CALL TrDYNAM   ! tracer dynamics
         CALL TIMER (MNOW,MTRACE)
#endif
C****
C**** Calculate tropopause level and pressure
C****
      CALL CALC_TROP
C**** calculate some dynamic variables for the PBL
      CALL PGRAD_PBL
C**** calculate zenith angle for current time step
      CALL CALC_ZENITH_ANGLE

         CALL CHECKT ('DYNAM ')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5D.EQ.0) CALL DIAG5A (7,NIdyn)
         IF (MODD5D.EQ.0) CALL DIAGCA (2)
         IF (MOD(Itime,NDAY/2).eq.0) CALL DIAG7A
C****
C**** INTEGRATE SOURCE TERMS
C****
         IDACC(1)=IDACC(1)+1
         MODD5S=MOD(Itime-ItimeI,NDA5S)
         IF (MODD5S.EQ.0) IDACC(8)=IDACC(8)+1
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG5A (1,0)
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAGCA (1)

C**** FIRST CALL MELT_SI SO THAT TOO SMALL ICE FRACTIONS ARE REMOVED
C**** AND ICE FRACTION CAN THEN STAY CONSTANT UNTIL END OF TIMESTEP
      CALL MELT_SI
         CALL UPDTYPE
         CALL TIMER (MNOW,MSURF)
C**** CONDENSATION, SUPER SATURATION AND MOIST CONVECTION
      CALL CONDSE
         CALL CHECKT ('CONDSE')
         CALL TIMER (MNOW,MCNDS)
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)
         IF (MODD5S.EQ.0) CALL DIAGCA (3)
      end if                                  ! full model,kradia le 0
C**** RADIATION, SOLAR AND THERMAL
      MODRD=MOD(Itime-ItimeI,NRAD)
      if (kradia.le.0. or. MODRD.eq.0) then
         CALL RADIA
         if (kradia.le.0) CALL CHECKT ('RADIA ')
      end if
         CALL TIMER (MNOW,MRAD)
      if (kradia.le.0) then                    ! full model,kradia le 0
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)
         IF (MODD5S.EQ.0) CALL DIAGCA (4)
C****
C**** SURFACE INTERACTION AND GROUND CALCULATION
C****
C**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
C**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
C****
C**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      CALL PRECIP_SI
      CALL PRECIP_LI
C**** APPLY PRECIPITATION AND RUNOFF TO LAKES/OCEANS
      CALL PRECIP_LK
      CALL PRECIP_OC
         CALL TIMER (MNOW,MSURF)
         CALL CHECKT ('PRECIP')
#ifdef TRACERS_ON
C**** Calculate non-interactive tracer surface sources and sinks
         call set_tracer_2Dsource
         CALL TIMER (MNOW,MTRACE)
#endif
C**** CALCULATE SURFACE FLUXES AND EARTH
      CALL SURFCE
         CALL CHECKT ('SURFCE')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (5)
C**** CALCULATE ICE DYNAMICS
      CALL DYNSI
C**** CALCULATE BASE ICE-OCEAN/LAKE FLUXES
      CALL UNDERICE
C**** APPLY SURFACE/BASE FLUXES TO SEA/LAKE ICE
      CALL GROUND_SI
C**** APPLY SURFACE FLUXES TO LAND ICE
      CALL GROUND_LI
         CALL CHECKT ('GRNDSI')
C**** APPLY FLUXES TO LAKES AND DETERMINE ICE FORMATION
      CALL GROUND_LK
         CALL CHECKT ('GRNDLK')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (6)
C**** CALCULATE RIVER RUNOFF FROM LAKE MASS
      CALL RIVERF
      CALL GROUND_E    ! diagnostic only - should be merged with EARTH
C**** APPLY FLUXES TO OCEAN, DO OCEAN DYNAMICS AND CALC. ICE FORMATION
      CALL OCEANS
         CALL CHECKT ('OCEANS')
C**** APPLY ICE FORMED IN THE OCEAN/LAKES TO ICE VARIABLES
      CALL FORM_SI
         CALL CHECKT ('FORMSI')
C**** IF ATURB is used in rundeck then this is a dummy call
C**** CALCULATE DRY CONVECTION ABOVE PBL
      CALL ATM_DIFFUS (2,LM-1,dtsrc)
         CALL CHECKT ('DRYCNV')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (9)
C**** ADVECT ICE
      CALL ADVSI
      CALL ADVSI_DIAG ! needed to update qflux model, dummy otherwise
         CALL CHECKT ('ADVSI ')
C**** UPDATE DIAGNOSTIC TYPES
      CALL UPDTYPE
C**** ADD DISSIPATED KE FROM SURFACE CALCULATION BACK AS LOCAL HEAT
      CALL DISSIP
         CALL CHECKT ('DISSIP')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (7)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)

C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(Itime-ItimeI,NFILTR).EQ.0) THEN
           IDACC(10)=IDACC(10)+1
           IF (MODD5S.NE.0) CALL DIAG5A (1,0)
           CALL DIAGCA (1)
           CALL FILTER
           CALL CHECKT ('FILTER')
           CALL TIMER (MNOW,MDYN)
           CALL DIAG5A (14,NFILTR*NIdyn)
           CALL DIAGCA (8)
      END IF
#ifdef TRACERS_ON
C**** 3D Tracer sources and sinks
C**** Tracer gravitational settling for aerosols
      CALL TRGRAV
C**** Tracer radioactive decay (and possible source)
      CALL TDECAY
C**** Calculate 3D tracers sources and sinks

      call tracer_3Dsource
C**** Accumulate tracer distribution diagnostics
      CALL TRACEA
         CALL TIMER (MNOW,MTRACE)
         CALL CHECKT ('T3DSRC')
#endif
      end if                                  ! full model,kradia le 0
C****
C**** WRITE SUB-DAILY DIAGNOSTICS EVERY NSUBDD hours
C****
      if (Nsubdd.ne.0) then
        if (mod(Itime+1,Nsubdd).eq.0) call get_subdd
      end if
#ifdef TRACERS_DUST
      call ahourly
#endif
C****
C**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
C****
  100 Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1
      Jhour=MOD(Itime*24/NDAY,24)         ! Hour (0-23)
      Nstep=Nstep+NIdyn                   ! counts DT(dyn)-steps

      IF (MOD(Itime,NDAY).eq.0) THEN      ! NEW DAY
      if (kradia.gt.0) then               ! radiative forcing run
        CALL DAILY(.false.)
        if(Kradia<10)  CALL daily_RAD(.true.)
        if(Kradia==10) CALL daily_OCEAN(.true.) ! to test OCLIM
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
      else                                ! full model, kradia le 0
           CALL DIAG5A (1,0)
           CALL DIAGCA (1)
        CALL DAILY(.true.)                 ! end_of_day
        CALL daily_RAD(.true.)
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
           CALL TIMER (MNOW,MELSE)

        call daily_LAKE
        call daily_EARTH(.true.)           ! end_of_day

        call daily_OCEAN(.true.)           ! end_of_day
        call daily_ICE
        call daily_LI
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
        call daily_tracer(1)
           CALL TIMER (MNOW,MTRACE)
#endif
           CALL CHECKT ('DAILY ')
           CALL TIMER (MNOW,MSURF)
           CALL DIAG5A (16,NDAY*NIdyn)
           CALL DIAGCA (10)
        call sys_flush(6)
      end if   ! kradia: full model (or rad.forcing run)
      CALL UPDTYPE
      END IF   !  NEW DAY

#ifdef USE_FVCORE
      Call Compute_Tendencies(fv)
#endif

      if (kradia.le.0) then   ! full model
C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY 24 HOURS
C****
      IF (Kvflxo.EQ.0.) OA(:,:,4:KOA)=0. ! to prepare for future saves
      IF (Kvflxo.NE.0.) THEN
         IF (MOD(Itime,NDAY).eq.0) THEN
            call pack_data (grid, OA, OA_glob)
            if (am_I_root())
     *         call WRITEI8 (iu_vflxo,Itime,OA_glob,im*jm*koa)
C**** ZERO OUT INTEGRATED QUANTITIES
            OA(:,:,4:KOA)=0.
         ELSEIF (MOD(Itime,NDAY/2).eq.0) THEN
            call vflx_OCEAN
         END IF
         CALL TIMER (MNOW,MELSE)
      END IF

C****
C**** CALL DIAGNOSTIC ROUTINES
C****
#ifdef SCM
c*****call scm diagnostics every time step
      call scm_diag
#endif

      IF (MOD(Itime-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 E-history
C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
      IF (NIPRNT.GT.0) THEN
        acc_period='PARTIAL      '
        call print_diags(1)
        NIPRNT=NIPRNT-1
        call set_param( "NIPRNT", NIPRNT, 'o' )
      END IF

      end if   ! full model ; kradia le 0

C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
C****    done at the end of (selected) months
      IF (months.ge.NMONAV .and.   ! next 2 conditions are rarely needed
     *    JDAY.eq.1+JDendOfM(JMON-1) .and. MOD(Itime,NDAY).eq.0) THEN

C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
      call aPERIOD (JMON0,JYEAR0,months,1,0, aDATE(1:12),Ldate)
      acc_period=aDATE(1:12)
      WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),JYEAR
      if (kradia.le.0) call print_diags(0)
C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.12) m = m-12
            monacc(m) = 1
          end do
          If (AM_I_ROOT())
     *       call openunit(aDATE(1:7)//'.acc'//XLABEL(1:LRUNID),iu_ACC,
     *        .true.,.false.)
          call io_rsf (iu_ACC,Itime,iowrite_single,ioerr)
          If (AM_I_ROOT()) call closeunit(iu_ACC)
C**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
          IF (KCOPY.GT.1) THEN
            CALL RFINAL (IRAND)
            call set_param( "IRAND", IRAND, 'o' )
            IF (AM_I_ROOT())
     *          call openunit('1'//aDATE(8:14)//'.rsf'//XLABEL(1:LRUNID)
     *           ,iu_RSF,.true.,.false.)
            call io_rsf(iu_RSF,Itime,iowrite_mon,ioerr)
            IF (AM_I_ROOT()) call closeunit(iu_RSF)
#ifdef USE_FVCORE
            fv_fname  = '1'//aDATE(8:14)//'.fv'//XLABEL(1:LRUNID)
            fv_dfname = '1'//aDATE(8:14)//'.dfv'//XLABEL(1:LRUNID)
            call Checkpoint(fv, clock, fv_fname, fv_dfname)
#endif
          END IF
C**** KCOPY > 2 : ALSO SAVE THE OCEAN DATA TO INITIALIZE DEEP OCEAN RUNS
          IF (KCOPY.GT.2) THEN
            If (AM_I_ROOT())
     *           call openunit(aDATE(1:7)//'.oda'//XLABEL(1:LRUNID)
     *           ,iu_ODA,.true.,.false.)
            call io_oda(iu_ODA,Itime,iowrite,ioerr)
            IF (AM_I_ROOT()) call closeunit(iu_ODA)
          END IF
        END IF

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (MNOW,MDIAG)
        TOTALT=SUM(TIMING(1:NTIMEACC))
        DO M=1,NTIMEACC
          PERCENT(M) = 100d0*TIMING(M)/(TOTALT+.00001)
        END DO
c        TOTALT=(MNOW-MSTART)/(60.*100.) ! wrong when clock rolls over
        TOTALT=TOTALT/(60.*100.)           ! in minutes
        DTIME = NDAY*TOTALT/(Itime-Itime0)        ! minutes/day
        WRITE (6,'(/A,F7.2,A,/(8(A13,F5.1/))//)')
     *   '0TIME',DTIME,'(MINUTES) ',(TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        TIMING = 0
        MSTART= MNOW
      END IF

C**** CPU TIME FOR CALLING DIAGNOSTICS
      CALL TIMER (MNOW,MDIAG)
C**** TEST FOR TERMINATION OF RUN
ccc
      IF (MOD(Itime,Nssw).eq.0) then
       IF (AM_I_ROOT()) then
        flg_go = '__STOP__'     ! stop if flagGoStop if missing
        iflag=0
        open(3,file='flagGoStop',form='FORMATTED',status='OLD',err=210)
        read (3,'(A8)',end=210) flg_go
        close (3)
 210    continue
        IF (flg_go .eq. '___GO___') iflag=1
        call ESMF_BCAST( grid, iflag)
       else
        call ESMF_BCAST( grid, iflag)
        if (iflag .eq. 1) flg_go = '___GO___'
        if (iflag .eq. 0) flg_go = '__STOP__'
       end if
      endif
      IF (flg_go.ne.'___GO___' .or. stop_on) THEN
C**** Flag to continue run has been turned off
         WRITE (6,'("0Flag to continue run has been turned off.")')
         EXIT
      END IF

c$$$      call test_save(__LINE__, itime-1)
      END DO

      call gettime(tloopend)
      if (AM_I_ROOT())
     *     write(*,*) "Time spent in the main loop in seconds:",
     *     .01*(tloopend-tloopbegin)
C****
C**** END OF MAIN LOOP
C****

C**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      IF (AM_I_ROOT())
     *     call openunit(rsf_file_name(KDISK),iu_RSF,.true.,.false.)
      call io_rsf(iu_RSF,Itime,iowrite,ioerr)
      IF (AM_I_ROOT()) call closeunit(iu_RSF)
#ifdef USE_FVCORE
         fv_fname='fv.' ; write(fv_fname(4:4),'(i1)') kdisk
         fv_dfname='dfv.' ; write(fv_dfname(5:5),'(i1)') kdisk
         call Finalize(fv, clock, fv_fname, fv_dfname)
#endif
      if (AM_I_ROOT()) then
      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *  '0Restart file written on fort.',KDISK,'Year',JYEAR,
     *     aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
      end if

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)
#ifdef SCM
      call closeunit(iu_scm_prt)
      call closeunit(iu_scm_diag)
#endif
      IF (AM_I_ROOT())
     *   WRITE (6,'(/////4(1X,33("****")/)//,A,I8
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME

      IF (Itime.ge.ItimeE) CALL stop_model (
     &     'Terminated normally (reached maximum time)',13)
      CALL stop_model ('Run stopped with sswE',12)  ! voluntary stop

      END


      subroutine sig_stop_model
      USE MODEL_COM, only : stop_on
      implicit none
      stop_on = .true.
      end subroutine sig_stop_model


      subroutine init_Model
!@sum This program reads most of parameters from the database (DB)
!@+   get_param( "A", X ) reads parameter A into variable X
!@+   if "A" is not in the database, it will generate an error
!@+   message and stop
!@+   sync_param( "B", Y ) reads parameter B into variable Y
!@+   if "B" is not in the database, then Y is unchanged and its
!@+   value is saved in the database as "B" (here sync = synchronize)
      USE MODEL_COM, only : JM,LM,NIPRNT,MFILTR,NFILTR,NRAD
     *     ,NDASF,NDA4,NDA5S,NDA5K,NDA5D,NDAA,Kvflxo,kradia
     *     ,NMONAV,Ndisk,Nssw,KCOPY,KOCEAN,NIsurf,iyear1
     $     ,LS1,IRAND,ItimeI,PSTRAT,UOdrag
     $     ,X_SDRAG,C_SDRAG,LSDRAG,P_SDRAG,LPSDRAG,PP_SDRAG,ang_sdrag
     $     ,P_CSDRAG,CSDRAGL,Wc_Jdrag,wmax,COUPLED_CHEM,dt
     *     ,DT_XUfilter,DT_XVfilter,DT_YVfilter,DT_YUfilter,QUVfilter
     &     ,do_polefix,pednl00,pmidl00
      USE DOMAIN_DECOMP, only: AM_I_ROOT
      USE PARAM
      implicit none
      INTEGER L,LCSDRAG

C**** Rundeck parameters:
      call sync_param( "NMONAV", NMONAV )
      call sync_param( "NIPRNT", NIPRNT )
      call sync_param( "DT_XVfilter", DT_XVfilter )
      call sync_param( "DT_XUfilter", DT_XUfilter )
      call sync_param( "DT_YVfilter", DT_YVfilter )
      call sync_param( "DT_YUfilter", DT_YUfilter )
      call sync_param( "MFILTR", MFILTR )
      call sync_param( "X_SDRAG", X_SDRAG, 2 )
      call sync_param( "C_SDRAG", C_SDRAG )
      call sync_param( "P_CSDRAG", P_CSDRAG )
      call sync_param( "P_SDRAG", P_SDRAG )
      call sync_param( "PP_SDRAG", PP_SDRAG )
      call sync_param( "ANG_SDRAG", ANG_SDRAG )
      call sync_param( "Wc_Jdrag", Wc_Jdrag )
      call sync_param( "wmax", wmax )
      call sync_param( "do_polefix", do_polefix )
      call sync_param( "NDASF", NDASF )
      call sync_param( "NDA4", NDA4 ) !!
      call sync_param( "NDA5S", NDA5S ) !!
      call sync_param( "NDA5K", NDA5K ) !!
      call sync_param( "NDA5D", NDA5D ) !!
      call sync_param( "NDAA", NDAA ) !!
      call sync_param( "NFILTR", NFILTR ) !!
      call sync_param( "NRAD", NRAD ) !!
      call sync_param( "Kvflxo", Kvflxo ) !!
      call sync_param( "Ndisk", Ndisk )
      call sync_param( "Nssw", Nssw )
      call sync_param( "KCOPY", KCOPY )
      call sync_param( "KOCEAN", KOCEAN )
      call sync_param( "KRADIA", KRADIA )
      call sync_param( "NIsurf", NIsurf )
      call sync_param( "UOdrag", UOdrag )
      call sync_param( "IRAND", IRAND )
      call sync_param( "COUPLED_CHEM", COUPLED_CHEM )

C**** Non-Rundeck parameters

C**** Calculate levels for application of SDRAG: LSDRAG,LPSDRAG->LM i.e.
C**** all levels above and including P_SDRAG mb (PP_SDRAG near poles)
C**** If P is the edge between 2 levels, take the higher level.
C**** Also find CSDRAGL, the coefficients of C_Sdrag as a function of L

      LSDRAG=LM ; LPSDRAG=LM ; LCSDRAG=LM ; CSDRAGL=C_SDRAG
      DO L=1,LM
        IF (PEDNL00(L+1)-1d-5.lt.P_SDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.P_SDRAG)         LSDRAG=L
        IF (PEDNL00(L+1)-1d-5.lt.PP_SDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.PP_SDRAG)        LPSDRAG=L
        IF (PEDNL00(L+1)-1d-5.lt.P_CSDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.P_CSDRAG)        LCSDRAG=L
      END DO
      DO L=LCSDRAG,LSDRAG-1
         CSDRAGL(L) = C_SDRAG + max( 0.d0 , (X_SDRAG(1)-C_SDRAG) *
     *     LOG(P_CSDRAG/(PMIDL00(L))) / LOG(P_CSDRAG/P_SDRAG) )
      END DO
      if (AM_I_ROOT()) then
         WRITE(6,*) "Levels for  LSDRAG =",LSDRAG ,"->",LM
         WRITE(6,*) "Levels for LPSDRAG =",LPSDRAG,"->",LM," near poles"
         WRITE(6,*) "C_SDRAG coefficients:",CSDRAGL(LS1:LSDRAG-1)
      end if

C**** Determine if FLTRUV is called.
      QUVfilter = .false.
      if (DT_XUfilter>0. .or. DT_XVfilter>0. .or.
     *    DT_YUfilter>0. .or. DT_YVfilter>0.)  QUVfilter = .true.
      if (QUVfilter) then
         if (DT_XUfilter > 0. .and. DT_XUfilter < DT) then
             DT_XUfilter = DT
             WRITE(6,*) "DT_XUfilter too small; reset to :",DT_XUfilter
         end if
         if (DT_XVfilter > 0. .and. DT_XVfilter < DT) then
             DT_XVfilter = DT
             WRITE(6,*) "DT_XVfilter too small; reset to :",DT_XVfilter
         end if
         if (DT_YUfilter > 0. .and. DT_YUfilter < DT) then
             DT_YUfilter = DT
             WRITE(6,*) "DT_YUfilter too small; reset to :",DT_YUfilter
         end if
         if (DT_YVfilter > 0. .and. DT_YVfilter < DT) then
             DT_YVfilter = DT
             WRITE(6,*) "DT_YVfilter too small; reset to :",DT_YVfilter
         end if
      end if
c Warn if polar fixes requested for a model not having a half polar box
c     if(do_polefix.eq.1 .and. jm.ne.46) then
c        do_polefix = 0
c        write(6,*) 'Polar fixes are currently applicable only to'//
c    &           'models having a half polar box; no fixes applied'
c     endif
      RETURN
C****
      end subroutine init_Model

#ifdef USE_FVCORE
      SUBROUTINE INPUT (istart,ifile,clock)
#else
      SUBROUTINE INPUT (istart,ifile)
#endif
C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE FILEMANAGER, only : openunit,closeunit,nameunit
      USE TIMINGS, only : timing,ntimeacc
      USE PARAM
      USE PARSER
      USE CONSTANT, only : grav,kapa,sday,by3
      USE MODEL_COM, only : im,jm,lm,wm,u,v,t,p,q,fearth0,fland
     *     ,focean,flake0,flice,hlake,zatmo,plbot,sig,dsig,sige,kradia
     *     ,bydsig,xlabel,lrunid,nmonav,qcheck,irand,ptop
     *     ,nisurf,nidyn,nday,dt,dtsrc,kdisk,jmon0,jyear0
     *     ,iyear1,itime,itimei,itimee
     *     ,ls1,psfmpt,pstrat,idacc,jyear,jmon,jday,jdate,jhour
     *     ,aMONTH,jdendofm,jdpery,aMON,aMON0,ioread,irerun
     *     ,ioread_single,irsfic,irsficnt,iowrite_single,ioreadnt
     *     ,irsficno,mdyn,mcnds,mrad,msurf,mdiag,melse,Itime0,Jdate0
     *     ,Jhour0,rsf_file_name,lm_req
     *     ,pl00,aml00,pednl00,pdsigl00,pmidl00,byaml00
#ifdef SCM
     *     ,I_TARG,J_TARG
#endif
#ifdef SCM
      USE SCMCOM, only : iu_scm_prt
#endif
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : geom_b,imaxj
      USE RANDOM
      USE RAD_COM, only : rqt,cloud_rad_forc
      USE DYNAMICS, only : pk,pmid,pedn
      USE CLOUDS_COM, only : ttold,qtold,svlhx,rhsav,cldsav
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM,only: MTRACE,NTM,TRNAME
#ifdef TRACERS_SPECIAL_Shindell
     *     ,mchem
#endif
#endif
#ifdef TRACERS_AMP
      USE AERO_CONFIG
      USE AERO_COAG
      USE AERO_INIT
      USE AERO_SETUP
      USE AERO_SUBS
      USE AERO_NPF
      USE AERO_DIAM
      USE AMP_AEROSOL
#endif
      USE DIAG_COM, only : acc_period,monacc,jreg,titreg,namreg
     &  ,hr_in_day,iwrite,jwrite,itwrite,kdiag,qdiag,qdiag_ratios,oa
      USE PBLCOM
     &     , only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar_pbl
     &  ,egcm,w2gcm,tgvavg,qgavg
      USE LAKES_COM, only : flake
      USE GHY_COM, only : fearth
      USE SOIL_DRV, only: init_gh
      USE DOMAIN_DECOMP, only : grid, GET, READT_PARALLEL, AM_I_ROOT
      USE DOMAIN_DECOMP, only : HALO_UPDATE, NORTH, HERE
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: init_app_clock
      USE CONSTANT, only : hrday
      USE ESMF_MOD, only: ESMF_Clock
#endif
      USE ATMDYN, only : init_ATMDYN,CALC_AMPK
#ifdef USE_ENT
      USE ENT_DRV, only : init_module_ent
#endif

      IMPLICIT NONE
      CHARACTER(*) :: ifile
!@var iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_RSF unit numbers for input files
      INTEGER iu_AIC,iu_TOPO,iu_GIC,iu_REG,iu_RSF,iu_IFILE
!@var num_acc_files number of acc files for diag postprocessing
      INTEGER I,J,L,K,LID1,LID2,ITYPE,IM1,NOFF,ioerr,num_acc_files
!@nlparam HOURI,DATEI,MONTHI,YEARI        start of model run
!@nlparam TIMEE,HOURE,DATEE,MONTHE,YEARE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
     *    TIMEE=-1,HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHOURE=-1,
!@nlparam ISTART  postprocessing(-1)/start(1-8)/restart(>8)  option
!@nlparam IRANDI  random number seed to perturb init.state (if>0)
     *             ISTART, IRANDI=0
      REAL*8 TIJL,CDM,TEMP,X
      INTEGER Itime1,Itime2,ItimeX,IhrX, LMR

!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5d0

      LOGICAL :: redoGH = .FALSE.,iniPBL = .FALSE., inilake = .FALSE.,
     &           iniSNOW = .FALSE.  ! true = restart from "no snow" rsf
     &           ,iniOCEAN = .FALSE.
#ifdef USE_ENT
     &     ,iniENT = .FALSE.
#endif
#ifdef USE_FVCORE
      type (ESMF_Clock) :: clock
      integer :: minti,minte
      character(len=1) :: suffix
#endif
      CHARACTER NLREC*80,filenm*100,RLABEL*132
      NAMELIST/INPUTZ/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,QDIAG,KDIAG,QDIAG_RATIOS
     *     ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI
      integer ISTART_kradia, nl_soil

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
CCCC      INTEGER :: stdin ! used to read 'I' file
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP  = J_1S,
     &               J_STRT_HALO= J_0H, J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****
C**** Default setting for ISTART : restart from latest save-file (10)
C****
      ISTART=10
      num_acc_files=0
C****
C**** Set dependent vertical resolution variables
C****
      SIGE(:) = (PLbot(:)-PTOP)/PSFMPT
      SIG(:)  = (sige(1:lm)+sige(2:lm+1))*0.5d0
      DSIG(:) =  sige(1:lm)-sige(2:lm+1)
      byDSIG  =  1./DSIG
C**** CALCULATE SPHERICAL GEOMETRY
      CALL GEOM_B
C**** Calculate default vertical arrays (including rad. eq. layers)
      LMR=LM+LM_REQ
      CALL CALC_VERT_AMP(PSFMPT,LMR,PL00,AML00,PDSIGL00,PEDNL00,PMIDL00)
      BYAML00(:)=1./AML00(:)
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
      RHSAV (:,:,:)=.85d0
      CLDSAV(:,:,:)=0.
      SVLHX (:,:,:)=0.
      WM    (:,:,:)=0.
C****    Ocean info saved for ocean heat transport calculations
         OA = 0.
C**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
      KDIAG(13)=9
C**** Set global default timing descriptions
C**** Other speciality descriptions can be added/used locally
      NTIMEACC = 0
      CALL SET_TIMER("ATMOS. DYNAM",MDYN)
      CALL SET_TIMER("CONDENSATION",MCNDS)
      CALL SET_TIMER("   RADIATION",MRAD)
      CALL SET_TIMER("     SURFACE",MSURF)
      CALL SET_TIMER(" DIAGNOSTICS",MDIAG)
      CALL SET_TIMER("       OTHER",MELSE)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      CALL SET_TIMER("     TRACERS",MTRACE)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      CALL SET_TIMER("   CHEMISTRY",MCHEM)
#endif
C****
C**** Set some documentary parameters in the database
C****
      call set_param("IM",IM)
      call set_param("JM",JM)
      call set_param("LM",LM)
      call set_param("LS1",LS1)
      call set_param("PLBOT",Plbot,LM+1)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call set_param("NTM",NTM)
      call set_param("TRNAME",TRNAME,ntm)
#endif
C****
C**** Print Header and Label (2 lines) from rundeck
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      if (AM_I_ROOT()) WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(iu_IFILE,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      if (AM_I_ROOT()) WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label
      IF(AM_I_ROOT()) THEN
C****
C**** Print preprocessing options (if any are defined)
C****
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      write(6,*) 'This program includes tracer code'
#endif
#ifdef TRACERS_WATER
      write(6,*) '...and water tracer code'
#ifndef TRACERS_ON
      call stop_model(
     &' Water tracers need TRACERS_ON as well as TRACERS_WATER',255)
#endif
#endif
#ifdef TRACERS_OCEAN
      write(6,*) '...and ocean tracer code'
#endif
#ifdef TRACERS_SPECIAL_O18
      write(6,*) '...and water isotope code'
#ifndef TRACERS_WATER
      call stop_model('Water isotope tracers need TRACERS_WATER '//
     *'as well as TRACERS_SPECIAL_O18',255)
#endif
#endif
#ifdef TRACERS_SPECIAL_Lerner
      write(6,*) '...and Jean/David tracers and chemistry'
#endif
#ifdef TRACERS_GASEXCH_Natassa
      write(6,*) '          '
      write(6,*) '...and Natassa Romanou air-sea GAS EXCHANGE'
#ifdef TRACERS_GASEXCH_CFC_Natassa
      write(6,*) '****CFC flux across air/sea interface****'
#endif
#ifdef TRACERS_GASEXCH_CO2_Natassa
      write(6,*) '****CO2 flux across air/sea interface****'
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Drew Shindell tracers and chemistry'
#endif
#ifdef TRACERS_AEROSOLS_Koch
      write(6,*) '...and Dorothy Koch aerosols'
#endif
#ifdef TRACERS_DRYDEP
      write(6,*) '...and tracer dry deposition'
#endif
#ifdef EDGAR_HYDE_SOURCES
      write(6,*) '...and EDGAR HYDE sources instead of GISS'
#endif
#ifdef SHINDELL_STRAT_CHEM
      write(6,*) '...and Drew Shindell stratospheric chemistry'
#endif
#ifdef SHINDELL_STRAT_EXTRA
      write(6,*) '...and Drew Shindell extra strat tracers'
#endif
#ifdef regional_Ox_tracers
      write(6,*) '...and regional Ox tracers'
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      write(6,*) '...and interactive CH4 wetlands emissions'
#endif
#ifdef NUDGE_ON
      write(6,*) '...and nudging of meteorology'
#endif
      ENDIF ! AM_I_ROOT()
C****
C**** Read parameters from the rundeck to database and namelist
C****
      call parse_params(iu_IFILE)
      READ (iu_IFILE,NML=INPUTZ,ERR=900)
      call closeunit(iu_IFILE)

C**** Get those parameters which are needed in this subroutine
      if(is_set_param("DTsrc"))  call get_param( "DTsrc", DTsrc )
      if(is_set_param("DT"))     call get_param( "DT", DT )
      if(is_set_param("NIsurf")) call get_param( "NIsurf", NIsurf ) !
      if(is_set_param("IRAND"))  call get_param( "IRAND", IRAND )
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      if(is_set_param("Kradia")) call get_param( "Kradia", Kradia )

C***********************************************************************
C****                                                               ****
C****        Post-process one or more ACC-files : ISTART < 1        ****
C****                                                               ****
C***********************************************************************
      if (istart.le.0) then
        call reset_diag(1)
        monacc = 0
        do
          call nextarg(filenm, 0)
          if ( filenm == "" ) exit ! end of args
          call openunit(filenm,iu_AIC,.true.,.true.)
          call io_rsf(iu_AIC,itime,ioread_single,ioerr)
          call closeunit(iu_AIC)
          num_acc_files = num_acc_files + 1
        end do
        GO TO 500
      end if
      if (istart.ge.9 .or. Kradia.gt.0) go to 400
C***********************************************************************
C****                                                               ****
C****                  INITIAL STARTS - ISTART: 1 to 8              ****
C****                                                               ****
C****   Current settings: 1 - from defaults                         ****
C****                     2 - from observed data                    ****
C****                     3 - so far unused                         ****
C****                     4 - from coupled model M-file - reset ocn ****
C****                     5 - tracer run from M-file w/o tracers    ****
C****                     6 - pred.ocn run from M-file w/o ocn data ****
C****                     7 - from mod. II' M-file - reset snow/ocn ****
C****                     8 - from current model M-file - no resets ****
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
        IF (AM_I_ROOT())
     *   WRITE(6,*) 'Please choose a proper start year yearI, not',yearI
        call stop_model('INPUT: yearI not provided',255)
      END IF
      IF (Iyear1.lt.0) Iyear1 = yearI
      IhrI = HourI +
     +  HR_IN_DAY*(dateI-1 + JDendofM(monthI-1) + JDperY*(yearI-Iyear1))
      ITimeI = IhrI*NDAY/HR_IN_DAY ! internal clock counts DTsrc-steps
      Itime=ItimeI
      IF (IhrI.lt.0) then
        IF (AM_I_ROOT())
     *  WRITE(6,*) 'Improper start time OR Iyear1=',Iyear1,' > yearI;',
     *  ' yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
        call stop_model(
     &       'INPUT: Improper start date or base year Iyear1',255)
      END IF
C**** Check the vertical layering defined in RES_ (is sige(ls1)=0 ?)
      IF (SIGE(LS1).ne.0.) then
        if (AM_I_ROOT())
     *       write(6,*) 'bad vertical layering: ls1,sige(ls1)',
     &       ls1,sige(ls1)
        call stop_model('INPUT: ls1 incorrectly set in RES_',255)
      END IF
C****
C**** Get Ground conditions from a separate file - ISTART=1,2
C****
      IF (ISTART.LE.2) THEN

C**** Set flag to initialise pbl and snow variables
        iniPBL=.TRUE.
        iniSNOW = .TRUE.  ! extract snow data from first soil layer
        iniOCEAN = .TRUE. ! read in ocean ic
#ifdef USE_ENT
        iniENT = .TRUE.
#endif
        if (istart.eq.1) redogh=.true.
C**** Read in ground initial conditions
        call openunit("GIC",iu_GIC,.true.,.true.)
        ioerr=-1
        read(iu_GIC)  ! ignore first line (ocean ic done in init_OCEAN)
        call io_seaice (iu_GIC,ioreadnt,ioerr)
        call io_earth  (iu_GIC,ioreadnt,ioerr)
        call io_soils  (iu_GIC,ioreadnt,ioerr)
        call io_landice(iu_GIC,ioreadnt,ioerr)
        if (ioerr.eq.1) then
          IF (AM_I_ROOT())
     *          WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
          call stop_model("INPUT: GIC READ IN ERROR",255)
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
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),0,P,1) ! Psurf
        DO J=J_0,J_1
          DO I=1,IM
            P(I,J)=P(I,J)-PTOP                        ! Psurf -> P
          END DO
        END DO
        DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),0,U(:,:,L),1) ! U
        END DO
        DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),0,V(:,:,L),1) ! V
        END DO
        DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),0,T(:,:,L),1) ! Temperature
        END DO
        DO L=1,LM  ! alternatively, only read in L=1,LS1 ; skip rest
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),0,Q(:,:,L),1) ! Q
        END DO
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),0,TSAVG,1)  ! Tsurf
      END IF
C****
C**** Derive other data from primary data if necessary - ISTART=1,2
C****                                                    currently
      IF (ISTART.LE.2) THEN
         If (HAVE_SOUTH_POLE) THEN
           WSAVG(1:im,1)=SQRT(U(1,2,1)*U(1,2,1)+V(1,2,1)*V(1,2,1))
           USAVG(1:im,1)=U(1,2,1)
           VSAVG(1:im,1)=V(1,2,1)
         End If
         If (HAVE_NORTH_POLE) THEN
           WSAVG(1:im,JM)=SQRT(U(1,JM,1)*U(1,JM,1)+V(1,JM,1)*V(1,JM,1))
           USAVG(1:im,JM)=U(1,JM,1)
           VSAVG(1:im,JM)=V(1,JM,1)
         End If

         CALL HALO_UPDATE(grid, U, FROM=NORTH)
         CALL HALO_UPDATE(grid, V, FROM=NORTH)
#ifdef SCM
        DO J=J_0S,J_1S
        DO I=1,IM
          WSAVG(I,J)=SQRT(U(I,J,1)**2+V(I,J,1)**2)
          USAVG(I,J)=U(I,J,1)
          VSAVG(I,J)=V(I,J,1)
        END DO
        END DO
#else 
        DO J=J_0S,J_1S
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
#endif
        CDM=.001d0

#ifdef SCM
c      enter SCM part of INPUT 
       call openunit("scm.prt",iu_scm_prt,.false.,.false.)
       call sync_param( "I_TARG",I_TARG)
       call sync_param( "J_TARG",J_TARG)
       write(0,*) 'I/J Targets set ',I_TARG,J_TARG
       write(iu_scm_prt,*) 'I/J Targets set ',I_TARG,J_TARG

c      write(iu_scm_prt,*) 'before scm inputs L u v '
c      do L=1,LM
c         write(iu_scm_prt,'(a6,i5,2(f9.3))') 'l u v ',l,
c    &         u(i_targ,j_targ,l),v(i_targ,j_targ,l)
c      enddo

       if ((I_TARG.lt.1 .or. I_TARG.gt. 144) .or.
     &      (J_TARG.lt.2 .or. J_TARG.gt.89)) then
             write(0,*) 'Invalid grid coordinates for selected box ',
     &             I_TARG,J_TARG
             STOP 100
       endif
!      read scm data and initialize model 
!      note:  usavg,vsavg and wsavg filled from here
       call init_scmdata
       write(0,*) 'return from init_scmdata'
c      write(iu_scm_prt,*) 'return from init_scmdata'
c      do L=1,LM
c         write(iu_scm_prt,'(a6,i5,2(f10.4))') 'l u v ',l,
c    &         u(i_targ,j_targ,l),v(i_targ,j_targ,l)
c      enddo
#endif

        CALL CALC_AMPK(LM)

        DO J=J_0,J_1
        DO I=1,IM
C**** SET SURFACE MOMENTUM TRANSFER TAU0
          TAUAVG(I,J)=CDM*WSAVG(I,J)**2
C**** SET LAYER THROUGH WHICH DRY CONVECTION MIXES TO 1
          DCLEV(I,J)=1.
C**** SET SURFACE SPECIFIC HUMIDITY FROM FIRST LAYER HUMIDITY
          QSAVG(I,J)=Q(I,J,1)
          QGAVG(I,J)=Q(I,J,1)
          TGVAVG(I,J)=T(I,J,1)
C**** SET RADIATION EQUILIBRIUM TEMPERATURES FROM LAYER LM TEMPERATURE
          DO K=1,LM_REQ
            RQT(K,I,J)=T(I,J,LM)
          END DO
C**** REPLACE TEMPERATURE BY POTENTIAL TEMPERATURE
          DO L=1,LM
            T(I,J,L)=T(I,J,L)/PK(L,I,J)
            TTOLD(L,I,J)=T(I,J,L)
            QTOLD(L,I,J)=Q(I,J,L)
          END DO
C**** initialize egcm to be used in ATURB.f
          DO L=1,LM
            egcm(l,i,j)=egcm_init_max/(float(l)**2)
            w2gcm(l,i,j)=egcm(l,i,j)*2.*by3
          END DO
        END DO
        END DO
C**** Initialize surface friction velocity
        DO ITYPE=1,4
        DO J=J_0,J_1
        DO I=1,IM
          USTAR_pbl(I,J,ITYPE)=WSAVG(I,J)*SQRT(CDM)
        END DO
        END DO
        END DO
C**** INITIALIZE VERTICAL SLOPES OF T,Q
        call tq_zmom_init(t,q,PMID,PEDN)
      END IF
C****
C**** I.C from possibly older/incomplete MODEL OUTPUT, ISTART=3-8
C****
      SELECT CASE (ISTART)
      CASE (3)               ! just general hints - not to be used as is
C**** Read what's there and substitute rest as needed (as above)
C**** To be implemented as needed. Sometimes it is safer to
C**** combine the ground layers into 2 layers (top 10cm and rest) and
C**** set   redoGH  to .true.  (after major changes in the GH code or
C**** after changing to a new horizontal grid)
C     redoGH=.TRUE.
C**** Set flag to initialise pbl/snow variables if obsolete or missing
C     iniPBL=.TRUE.  ; iniSNOW = .TRUE.
        go to 890            !  not implemented; stop model
C****
C**** I.C FROM FULL MODEL RESTART FILE (but re-initialise ocean)
C****
      CASE (4)
        call io_rsf(iu_AIC,IhrX,irsficno,ioerr)
        if (ioerr.eq.1) goto 800
        iniOCEAN = .TRUE. ! read in ocean ic
C****
C**** I.C FROM FULL MODEL RESTART FILE (but no tracers)
C****
      CASE (5)             ! this model's rsf file, no tracers
        call io_rsf(iu_AIC,IhrX,irsficnt,ioerr)
        if (ioerr.eq.1) goto 800
C****
C**** I.C FROM RESTART FILE that may not match land-ocean mask  ISTART=6
C****
      CASE (6)             ! converted model II' (B399) format (no snow)
        call io_rsf(iu_AIC,IhrX,irsficno,ioerr)
        if (ioerr.eq.1) goto 800
        iniSNOW = .TRUE.      ! extract snow data from first soil layer
        inipbl  = .TRUE.      ! initialise pbl profiles
        iniOCEAN = .TRUE.     ! read in ocean ic
        redoGH=.TRUE.
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
        !iniSNOW = .TRUE.      ! extract snow data from first soil layer
        iniPBL=.TRUE.
        if (ioerr.eq.1) goto 800
      END SELECT
C**** Check consistency of starting time
      IF (ISTART.ge.3) THEN
        IF( (MOD(IHRI-IHRX,8760).ne.0) ) THEN
         WRITE (6,*) ' Difference in hours between ',
     *       'Starting date and Data date:',MOD(IHRI-IHRX,8760)
         WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
         call stop_model('INPUT: start date inconsistent with data',255)
        ENDIF
      END IF
C**** Set flag to initialise lake variables if they are not in I.C.
      IF (ISTART.lt.8) inilake=.TRUE.

      CALL CALC_AMPK(LM)
C****
!**** IRANDI seed for random perturbation of initial conditions (if/=0):
C****        tropospheric temperatures changed by at most 1 degree C
      IF (IRANDI.NE.0) THEN
        CALL RINIT (IRANDI)
        DO L=1,LS1-1
        call burn_random(n=im*(j_0-1))
        DO J=J_0,J_1
        DO I=1,IM
           TIJL=T(I,J,L)*PK(L,I,J)-1.+2*RANDU(X)
           T(I,J,L)=TIJL/PK(L,I,J)
        END DO
        END DO
        END DO
        IF (AM_I_ROOT())
     *       WRITE(6,*) 'Initial conditions were perturbed !!',IRANDI
      END IF
C**** Close "AIC" here if it was opened
      IF (ISTART.gt.1) call closeunit(iu_AIC)

      IF (AM_I_ROOT())
     *     WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
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
      CASE (1:9)                      !  diag.arrays are not read in
        call openunit("AIC",iu_AIC,.true.,.true.)
        if(istart.eq.9) call io_rsf(iu_AIC,Itime,irerun,ioerr)
#ifdef USE_FVCORE
        call system('cp AICfv  fv_internal_restart.dat')
        call system('cp AICdfv tendencies_checkpoint')
#endif
        if(istart.le.8) then         !  initial start of rad.forcing run
          call io_label(iu_AIC,Itime,ItimeX,irerun,ioerr)
          if (Kradia.gt.0) call io_rad (iu_AIC,irsfic,ioerr)
        end if
        call closeunit(iu_AIC)
        if (ioerr.eq.1) goto 800
        WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *    ISTART,', TIME=',Itime,' ',XLABEL(1:80) ! sho input file label
        XLABEL = RLABEL                        ! switch to rundeck label

        CALL CALC_AMPK(LM)
C****
!**** IRANDI seed for random perturbation of current state (if/=0)
C****        tropospheric temperatures are changed by at most 1 degree C
        IF (IRANDI.ne.0 .and. Kradia.le.0) THEN
          CALL RINIT (IRANDI)
          DO L=1,LS1-1
          DO J=J_0,J_1
          DO I=1,IM
             TIJL=T(I,J,L)*PK(L,I,J)-1.+2*RANDU(X)
             T(I,J,L)=TIJL/PK(L,I,J)
          END DO
          END DO
          END DO
          IF (AM_I_ROOT())
     *         WRITE(6,*) 'Current temperatures were perturbed !!',IRANDI
        END IF
        TIMING = 0
        GO TO 500
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
C****
C**** CHOOSE DATA SET TO RESTART ON
      CASE (10,13:)
         Itime1=-1
         call openunit(rsf_file_name(1),iu_RSF,.true.,.true.)
         READ (iu_RSF,ERR=410) Itime1
         call closeunit(iu_RSF)
  410    continue !REWIND 1
         Itime2=-1
         call openunit(rsf_file_name(2),iu_RSF,.true.,.true.)
         READ (iu_RSF,ERR=420) Itime2
         call closeunit(iu_RSF)
  420    continue !REWIND 2
         IF (Itime1+Itime2.LE.-2.) GO TO 850
                               KDISK=1
         IF (Itime2.GT.Itime1) KDISK=2
         IF (ISTART.GE.13)     KDISK=3-KDISK
      CASE (11,12)
                               KDISK=ISTART-10
      END SELECT
  430 continue
      call openunit(rsf_file_name(KDISK),iu_RSF,.true.,.true.)
#ifdef USE_FVCORE
        write(suffix,'(i1)') kdisk
        call system('cp  fv.'// suffix // ' fv_internal_restart.dat')
        call system('cp dfv.'// suffix // ' tendencies_checkpoint')
#endif
      CALL HERE(__FILE__//'::io_rsf',__LINE__ + 10000*KDISK)
      call io_rsf(iu_RSF,Itime,ioread,ioerr)
      call closeunit(iu_RSF)
      if (ioerr.eq.1) then
         if (istart.gt.10) go to 850  ! no 2nd chance if istart/=10
         KDISK=3-KDISK                ! try the earlier restart file
         WRITE (6,'(A,I1,A,I1)')
     *     ' Read Error on fort.',3-kdisk,' trying fort.',kdisk
         ISTART=110
         go to 430
      end if
      if (AM_I_ROOT())
     * WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
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
      if(is_set_param("Kradia")) call get_param( "Kradia", Kradia )

C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************
  600 CONTINUE

C**** initialize Lrunid (length of the identifying part of XLABEL)
C****
      lid1 = INDEX(XLABEL,'(') -1
      if (lid1.lt.1) lid1=17
      lid2 = INDEX(XLABEL,' ') -1
      if (lid2.lt.1) lid2=17
      LRUNID = min(lid1,lid2)
      IF (LRUNID.gt.16) call stop_model
     *     ('INPUT: Rundeck name too long. Shorten to 16 char or less'
     *     ,255)

C**** Update ItimeE only if YearE or IhourE is specified in the rundeck
C****
      if(timee.lt.0) timee=houre*nday/HR_IN_DAY
      IF(yearE.ge.0) ItimeE = (( (yearE-iyear1)*JDperY +
     *  JDendofM(monthE-1)+dateE-1 )*HR_IN_DAY )*NDAY/HR_IN_DAY + TIMEE
C**** Alternate (old) way of specifying end time
      if(IHOURE.gt.0) ItimeE=IHOURE*NDAY/HR_IN_DAY

C**** Check consistency of DTsrc (with NDAY) and dt (with NIdyn)
      if (is_set_param("DTsrc") .and. nint(sday/DTsrc).ne.NDAY) then
        if (AM_I_ROOT())
     *        write(6,*) 'DTsrc=',DTsrc,' has to stay at/be set to',SDAY/NDAY
        call stop_model('INPUT: DTsrc inappropriately set',255)
      end if
      DTsrc = SDAY/NDAY
      call set_param( "DTsrc", DTsrc, 'o' )   ! copy DTsrc into DB

      NIdyn=nint(dtsrc/dt)
#ifndef USE_FVCORE
C**** NIdyn=dtsrc/dt(dyn) has to be a multiple of 2
C****
      if(istart>0) then
        NIdyn = 2*nint(.5*dtsrc/dt)
        if (is_set_param("DT") .and. nint(DTsrc/dt).ne.NIdyn) then
          if (AM_I_ROOT())
     *        write(6,*) 'DT=',DT,' has to be changed to',DTsrc/NIdyn
          call stop_model('INPUT: DT inappropriately set',255)
        end if
      end if
#else
C**** need a clock to satisfy ESMF interfaces
      call getdte(itimei,nday,iyear1,YEARI,MONTHI,jday,DATEI,HOURI,amon)
      MINTI = nint(mod( mod(Itimei*hrday/Nday,hrday) * 60d0, 60d0))
      call getdte(itimee,nday,iyear1,YEARE,MONTHE,jday,DATEE,HOURE,amon)
      MINTE = nint(mod( mod(Itimee*hrday/Nday,hrday) * 60d0, 60d0))
      clock = init_app_clock( (/ YEARI, MONTHI, DATEI, HOURI, MINTI,0/),
     &                        (/ YEARE, MONTHE, DATEE, HOURE, MINTE,0/),
     &             interval = int(dt) )
#endif
      DT = DTsrc/NIdyn
      call set_param( "DT", DT, 'o' )         ! copy DT into DB

C**** NMONAV has to be 1(default),2,3,4,6,12, i.e. a factor of 12
      if (NMONAV.lt.1 .or. MOD(12,NMONAV).ne.0) then
        write (6,*) 'NMONAV has to be 1,2,3,4,6 or 12, not',NMONAV
        call stop_model('INPUT: nmonav inappropriately set',255)
      end if
      if (AM_I_ROOT())
     *     write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'

C**** Updating Parameters: If any of them changed beyond this line
C**** use set_param(.., .., 'o') to update them in the database (DB)

C**** Get the rest of parameters from DB or put defaults to DB
      call init_Model

C**** Set julian date information
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      call getdte(Itime0,Nday,iyear1,Jyear0,Jmon0,J,Jdate0,Jhour0,amon0)

C****
C**** READ IN TIME-INDEPENDENT ARRAYS
C****
      if (Kradia.le.0) then   !  full model
        CALL CALC_AMPK(LM)

C****   READ SPECIAL REGIONS FROM UNIT 29
        call openunit("REG",iu_REG,.true.,.true.)
        READ(iu_REG) TITREG,JREG,NAMREG
        IF (AM_I_ROOT())
     &       WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG
        call closeunit(iu_REG)
      end if  ! full model: Kradia le 0

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** Initialise tracer parameters and diagnostics
       call init_tracer
#endif
C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
C**** Note that FLAKE0 is read in only to provide initial values
C**** Actual array is set from restart file.
      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),0,FOCEAN,1) ! Ocean fraction
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),0,FLAKE0,1) ! Orig. Lake fraction
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),0,FEARTH0,1) ! Earth frac. (no LI)
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),0,FLICE ,1) ! Land ice fraction
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),0,ZATMO ,1) ! Topography
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),0,HLAKE ,2) ! Lake Depths
      ZATMO(:,J_0:J_1) = ZATMO(:,J_0:J_1)*GRAV                  ! Geopotential
      call closeunit(iu_TOPO)

C**** Check polar uniformity
      if (AM_I_ROOT()) then
        do j=1,jm,jm-1  ! only polar boxes
          do i=2,im
            if (zatmo(i,j).ne.zatmo(1,j)) then
              print*,"Polar topography not uniform, corrected",i,j
     *             ,zatmo(i,j),zatmo(1,j)
              zatmo(i,j)=zatmo(1,j)
            end if
          end do
        end do
      end if

C**** Initialise some modules before finalising Land/Ocean/Lake/LI mask
C**** Initialize ice
      CALL init_ice(iniOCEAN,istart)
C**** Initialize lake variables (including river directions)
      CALL init_LAKES(inilake,istart)
C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN(iniOCEAN,istart)
C**** Initialize ice dynamics code (if required)
      CALL init_icedyn(iniOCEAN)
C**** Initialize land ice (must come after oceans)
      CALL init_LI(istart)

C**** Make sure that constraints are satisfied by defining FLAND/FEARTH
C**** as residual terms. (deals with SP=>DP problem)
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          FLAND(I,J)=1.-FOCEAN(I,J) ! Land fraction
          IF (FLAKE(I,J).gt.0) THEN
           IF (AM_I_ROOT()) WRITE(6,*)
     *            "Ocean and lake cannot co-exist in same grid box"
     *       ,i,j,FOCEAN(I,J),FLAKE(I,J)
            FLAKE(I,J)=0
          END IF
        ELSEIF (FLAKE(I,J).gt.0) THEN
          FLAND(I,J)=1.-FLAKE(I,J)
        ELSE
          FLAND(I,J)=1.
        END IF
C**** Ensure that no round off error effects land with ice and earth
        IF (FLICE(I,J)-FLAND(I,J).gt.-1d-4 .and. FLICE(I,J).gt.0) THEN
          FLICE(I,J)=FLAND(I,J)
          FEARTH(I,J)=0.
        ELSE
          FEARTH(I,J)=FLAND(I,J)-FLICE(I,J) ! Earth fraction
        END IF
      END DO
      END DO
      
      If (HAVE_SOUTH_POLE) Then
         FLAND(2:IM,1)=FLAND(1,1)
         FEARTH(2:IM,1)=FEARTH(1,1)
         FLICE(2:IM,1)=FLICE(1,1)
      End If
      If (HAVE_NORTH_POLE) Then
         FLAND(2:IM,JM)=FLAND(1,JM)
         FEARTH(2:IM,JM)=FEARTH(1,JM)
         FLICE(2:IM,JM)=FLICE(1,JM)
      End If
C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS (INCL. VEGETATION)
C**** Recompute Ground hydrology data if redoGH (new soils data)
C****
!!! hack: make sure that ISTART_kradia==0 if Kradia>0
!!! do we need it ? I.A.
      ISTART_kradia = ISTART
      if ( Kradia.gt.0 ) ISTART_kradia = 0
      CALL init_GH(DTsrc/NIsurf,redoGH,iniSNOW,ISTART_kradia, nl_soil)
#ifdef USE_ENT
      CALL init_module_ent(iniENT, Jday, Jyear, FOCEAN) !!! FEARTH)
#endif
      if (Kradia.gt.0) then   !  radiative forcing run
        !CALL init_GH(DTsrc/NIsurf,redoGH,iniSNOW,0)
        !CALL init_module_ent(iniENT,grid,jday,dxyp)
        CALL init_RAD(istart)
        if(istart.lt.0) CALL init_DIAG(0,num_acc_files) !post-processing
        if (AM_I_ROOT()) Then
           WRITE (6,INPUTZ)
           call print_param( 6 )
           WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
           WRITE (6,*) "PLbot=",PLbot
        end if
        if(istart.lt.0)
     &       CALL stop_model ('Terminated normally, istart<0',13)
        return
      end if                  !  Kradia>0; radiative forcing run
!      CALL init_GH(DTsrc/NIsurf,redoGH,iniSNOW,ISTART)
!      CALL init_module_ent(iniENT,grid,jday,dxyp)
C**** Initialize pbl (and read in file containing roughness length data)
      if(istart.gt.0) CALL init_pbl(iniPBL)
C****
C**** Initialize the use of gravity wave drag diagnostics
C****
      CALL init_GWDRAG
C
C**** Initialize nudging
#ifdef NUDGE_ON
      if (istart.gt.0) CALL NUDGE_INIT
#endif
#ifdef TRACERS_AMP
      CALL SETUP_CONFIG
      CALL SETUP_SPECIES_MAPS
      CALL SETUP_DP0
      CALL SETUP_AERO_MASS_MAP
      CALL SETUP_COAG_TENSORS
      CALL SETUP_DP0
      CALL SETUP_KIJ
      CALL SETUP_EMIS
      CALL SETUP_KCI
      CALL SETUP_NPFMASS
      CALL SETUP_DIAM
#endif
C****
      if(istart.gt.0) CALL RINIT (IRAND)
      CALL FFT0 (IM)
      CALL init_CLD
      CALL init_DIAG(istart,num_acc_files) ! initialize for accumulation
      CALL UPDTYPE
      if(istart.gt.0) CALL init_QUS(grid,im,jm,lm)
      if(istart.gt.0) CALL init_ATMDYN
      CALL init_RAD(istart)
      if (AM_I_ROOT()) then
         WRITE (6,INPUTZ)
         call print_param( 6 )
         WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
         WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
         WRITE (6,*) "PLbot=",PLbot
      end if
C****
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  800 WRITE (6,'(A,I4/" ",A)')
     *  '0ERROR ENCOUNTERED READING AIC ISTART=', ISTART,XLABEL(1:80)
      call stop_model('INPUT: READ ERROR FOR AIC',255)
  830 WRITE(6,*) 'READ ERROR FOR GIC'
      call stop_model('INPUT: READ ERROR FOR GIC',255)
  850 WRITE (6,'(A)')
     *  '0ERRORS ON BOTH RESTART DATA SETS. TERMINATE THIS JOB'
      call stop_model('INPUT: ERRORS ON BOTH RESTART FILES',255)
  890 WRITE (6,'(A,I5)') '0INCORRECT VALUE OF ISTART',ISTART
      call stop_model('INPUT: ISTART-SPECIFICATION INVALID',255)
  900 write (6,*) 'Error in NAMELIST parameters'
      call stop_model('Error in NAMELIST parameters',255)
      END SUBROUTINE INPUT

      SUBROUTINE DAILY(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@ver  1.0
!@calls constant:orbit, calc_ampk, getdte
      USE MODEL_COM, only : im,jm,lm,ls1,ptop,psf,p,q
     *     ,itime,itimei,iyear1,nday,jdpery,jdendofm
     *     ,jyear,jmon,jday,jdate,jhour,aMON,aMONTH,ftype
      USE GEOM, only : areag,dxyp,imaxj
      USE DYNAMICS, only : byAM
      USE RADPAR, only : ghgam,ghgyr2,ghgyr1
      USE RAD_COM, only : RSDIST,COSD,SIND, dh2o,H2ObyCH4,ghg_yr,
     *     omegt,obliq,eccn
#ifdef TRACERS_WATER
      USE TRACER_COM, only: trm,tr_wd_type,nwater,tr_H2ObyCH4,itime_tr0
     *     ,ntm
#endif
      USE DIAG_COM, only : aj=>aj_loc,j_h2och4
      USE DOMAIN_DECOMP, only : grid, GET, GLOBALSUM, AM_I_ROOT
      USE ATMDYN, only : CALC_AMPK
      IMPLICIT NONE
      REAL*8 DELTAP,PBAR,SMASS,LAM,xCH4,EDPY,VEDAY
      REAL*8 :: SPRESS(grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER i,j,l,iy
      LOGICAL, INTENT(IN) :: end_of_day
#ifdef TRACERS_WATER
      INTEGER n
#endif
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)


C**** Tasks to be done at end of day and at each start or restart
C****
C**** CALCULATE THE DAILY CALENDAR
C****
      call getdte(Itime,Nday,iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
C**** This is for noon (GMT) for new day.

C**** The orbital calculation will need to vary depending on the kind
C**** of calendar adopted (i.e. a generic 365 day year, or a transient
C**** calendar including leap years etc.).  For transient calendars the
C**** JDAY passed to orbit needs to be adjusted to represent the number
C**** of days from Jan 1 2000AD.
c      EDPY=365.2425d0, VEDAY=79.3125d0  ! YR 2000AD
c      JDAY => JDAY + 365 * (JYEAR-2000) + appropriate number of leaps
C**** Default calculation (no leap, VE=Mar 21 hr 0)
      EDPY=365d0 ; VEDAY=79d0           ! Generic year
      CALL ORBIT (OBLIQ,ECCN,OMEGT,VEDAY,EDPY,REAL(JDAY,KIND=8)-.5
     *     ,RSDIST,SIND,COSD,LAM)

      IF (.not.(end_of_day.or.itime.eq.itimei)) RETURN

C**** Tasks to be done at end of day and at initial starts only
C****
C**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS
C****
C**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE
#ifndef SCM
      SMASS=0.
      DO J=J_0,J_1
        SPRESS(J)=0.
        DO I=1,IM
          SPRESS(J)=SPRESS(J)+P(I,J)
        END DO
        SPRESS(J) = SPRESS(J) * DXYP(J)
      END DO
      CALL GLOBALSUM(grid, SPRESS, SMASS, ALL=.TRUE.)
      PBAR=SMASS/AREAG+PTOP
C**** CORRECT PRESSURE FIELD FOR ANY LOSS OF MASS BY TRUNCATION ERROR
C****   except if it was just done (restart from itime=itimei)
      DELTAP=PSF-PBAR
      if(itime.eq.itimei .and. abs(deltap).lt.1.d-10) return
      P=P+DELTAP

      CALL CALC_AMPK(LS1-1)

      if (AM_I_ROOT()) then
         IF (ABS(DELTAP).gt.1d-6)
     *      WRITE (6,'(A25,F10.6/)') '0PRESSURE ADDED IN GMP IS',DELTAP
      end if
#endif

      IF (.not.end_of_day) RETURN

C**** Tasks to be done at end of day only
      if (H2ObyCH4.gt.0) then
C****   Add obs. H2O generated by CH4(*H2ObyCH4) using a 2 year lag
        iy = jyear - 2 - ghgyr1 + 1
        if (ghg_yr.gt.0) iy = ghg_yr - 2 - ghgyr1 + 1
        if (iy.lt.1) iy=1
        if (iy.gt.ghgyr2-ghgyr1+1) iy=ghgyr2-ghgyr1+1
        xCH4=ghgam(3,iy)*H2ObyCH4
c        If (AM_I_ROOT())
c     &    write(6,*) 'add in stratosphere: H2O gen. by CH4(ppm)=',xCH4

        do l=1,lm
        do j=J_0,J_1
        do i=1,imaxj(j)
          q(i,j,l)=q(i,j,l)+xCH4*dH2O(j,l,jmon)*byAM(l,i,j)
#ifdef TRACERS_WATER
C**** Add water to relevant tracers as well
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              select case (tr_wd_type(n))
              case (nWater)    ! water: add CH4-sourced water to tracers
                trm(i,j,l,n) = trm(i,j,l,n) +
     +                tr_H2ObyCH4(n)*xCH4*dH2O(j,l,jmon)*dxyp(j)
              end select
            end if
          end do
#endif
          aj(j,j_h2och4,:)=aj(j,j_h2och4,:)+
     +                       xCH4*dH2O(j,l,jmon)*ftype(:,i,j)
        end do
        end do
        If (HAVE_NORTH_POLE) q(2:im,jm,l)=q(1,jm,l)
        If (HAVE_SOUTH_POLE) q(2:im, 1,l)=q(1, 1,l)
#ifdef TRACERS_WATER
        do n=1,ntm
          If (HAVE_SOUTH_POLE) trm(2:im, 1,l,n)=trm(1, 1,l,n)
          If (HAVE_NORTH_POLE) trm(2:im,jm,l,n)=trm(1,jm,l,n)
        end do
#endif
        end do
      end if

      RETURN
      END SUBROUTINE DAILY

      SUBROUTINE CHECKT (SUBR)
!@sum  CHECKT Checks arrays for NaN/INF and reasonablness
!@auth Original Development Team
!@ver  1.0

C**** CHECKT IS TURNED ON BY SETTING QCHECK=.TRUE. IN NAMELIST
C**** REMEMBER TO SET QCHECK BACK TO .FALSE. AFTER THE ERRORS ARE
C**** CORRECTED.
      USE CONSTANT, only : tf
      USE MODEL_COM
      USE DYNAMICS, only : pk
      USE DOMAIN_DECOMP, only : grid, GET, AM_I_ROOT
      USE soil_drv, only : checke
      IMPLICIT NONE
      INTEGER I,J,L
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)


      IF (QCHECK) THEN
C**** Check all prog. arrays for Non-numbers
        CALL CHECK3(U,IM,JM,LM,SUBR,'u     ')
        CALL CHECK3(V,IM,JM,LM,SUBR,'v     ')
        CALL CHECK3(T,IM,JM,LM,SUBR,'t     ')
        CALL CHECK3(Q,IM,JM,LM,SUBR,'q     ')
        CALL CHECK3(P,IM,JM,1,SUBR,'p     ')
        CALL CHECK3(WM,IM,JM,LM,SUBR,'wm    ')

        DO J=J_0,J_1
        DO I=1,IM
          IF (Q(I,J,1).gt.1d-1)print*,SUBR," Q BIG ",i,j,Q(I,J,1:LS1)
          IF (T(I,J,1)*PK(1,I,J)-TF.gt.50.) print*,SUBR," T BIG ",i,j
     *         ,T(I,J,1:LS1)*PK(1:LS1,I,J)-TF
        END DO
        END DO
        DO L=1,LM
        DO J=J_0,J_1
        DO I=1,IM
          IF (Q(I,J,L).lt.0.) then
            print*,"After ",SUBR," Q < 0 ",i,j,Q(I,J,L)
            call stop_model('Q<0 in CHECKT',255)
          END IF
          IF (WM(I,J,L).lt.0.) then
            print*,"After ",SUBR," WM < 0 ",i,j,WM(I,J,L)
            call stop_model('WM<0 in CHECKT',255)
          END IF
        END DO
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
       CALL CHECKE(SUBR)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** check tracers
        CALL CHECKTR(SUBR)
#endif
      END IF

      RETURN
      END SUBROUTINE CHECKT


      subroutine read_options( qcrestart, ifile )
!@sum reads options from the command line (for now only one option)
!@auth I. Aleinov
!@ver 1.0
      implicit none
!@var qcrestart true if "-r" is present
      logical, intent(inout) :: qcrestart
      character(*),intent(inout)  :: ifile
      character*80 arg,arg1

      do
        call nextarg( arg, 1 )
        if ( arg == "" ) exit          ! end of args
        select case (arg)
        case ("-r")
          qcrestart = .true.
        case ("-i")
          call nextarg( arg1, 0 )
          ifile=arg1
        ! new options can be included here
        case default
          print *,'Unknown option specified: ', arg
          print *,'Aborting...'
          call stop_model("Unknown option on a command line",255)
        end select
      enddo

      return
      end subroutine read_options


      subroutine print_restart_info
!@sum prints timing information needed to restart the model
!@auth I. Aleinov
!@ver 1.0
      USE MODEL_COM
      USE FILEMANAGER, only : openunit,closeunit
      use DOMAIN_DECOMP, only: AM_I_ROOT
      implicit none
      integer :: ItimeMax=-1, Itime1, Itime2, itm, ioerr1=-1, ioerr2=-1
      integer :: iu_rsf

      call openunit(rsf_file_name(1),iu_rsf,.true.)
      call io_label(iu_rsf,Itime1,itm,ioread,ioerr1)
      call closeunit(iu_rsf)
      call openunit(rsf_file_name(2),iu_rsf,.true.)
      call io_label(iu_rsf,Itime2,itm,ioread,ioerr2)
      call closeunit(iu_rsf)

      if ( ioerr1==-1 ) ItimeMax = Itime1
      if ( ioerr2==-1 ) ItimeMax = max( ItimeMax, Itime2 )

      if ( Itime < 0 )
     $     call stop_model("Could not read fort.1, fort.2",255)

      call getdte(
     &     ItimeMax,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      if (AM_I_ROOT())
     &     write(6,"('QCRESTART_DATA: ',I10,1X,I2,'-',I2.2,'-',I4.4)")
     &     ItimeMax*24/Nday, Jmon, Jdate, Jyear

      return
      end subroutine print_restart_info

