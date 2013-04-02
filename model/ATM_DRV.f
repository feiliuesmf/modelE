#include "rundeck_opts.h"

      subroutine atm_phase1
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE Dictionary_mod
      use resolution, only : im,jm,lm,ls1,ptop
      USE MODEL_COM
      USE ATM_COM, only : p,wm
      USE ATM_COM, only : pua,pva,sd_clouds,ptold,ps,kea
      USE DYNAMICS, only : nstep,nidyn,nfiltr,mfiltr,dt,conv
      USE DOMAIN_DECOMP_ATM, only: grid
      use domain_decomp_atm, only: writei8_parallel
      USE RANDOM
      USE GETTIME_MOD
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: mtrace
#ifdef TRAC_ADV_CPU
      USE TRACER_COM, only: mtradv
#endif
#endif
      USE DIAG_COM, only : ia_src,ia_d5s,ia_d5d,ia_filt
     &     ,oa,koa
     &     ,MODD5S,NDAa, NDA5d,NDA5s
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: Run,fvstate
#endif
#ifndef CUBED_SPHERE
      USE ATMDYN, only : DYNAM,SDRAG
     &     ,FILTER, COMPUTE_DYNAM_AIJ_DIAGNOSTICS
#endif
#ifdef SCM
      USE ATM_COM, only : t,p,q
      USE SCMCOM , only : SG_CONV,SCM_SAVE_T,SCM_SAVE_Q,
     &    iu_scm_prt,iu_scm_diag,I_TARG,J_TARG,nstepscm
#endif
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only : NBINS, IDTNUMD,IDTSO4,IDTECIL, IDTECOB,
     &     IDTOCIL, IDTOCOB,IDTDUST,IDTH2O,IDTNA
#endif
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use SystemTimers_mod
      use RAD_COM, only : nrad,modrd
      use seaice_com, only : si_atm,si_ocn,iceocn ! temporary until
      use lakes_com, only : icelak                ! melt_si calls
      use fluxes, only : atmocn,atmice            ! are moved
      implicit none

      INTEGER K,M,MSTART,MNOW,MODD5D,months,ioerr,Ldate,istart
      INTEGER :: MDUM = 0

      REAL*8 start,now, DTIME,TOTALT
#ifdef TRACERS_TOMAS
      integer :: n
#endif
      integer :: I,J,L,I_0,I_1,J_0,J_1
      real*8 :: initialTotalEnergy, finalTotalEnergy
      real*8 :: gettotalenergy ! external for now

      I_0 = GRID%I_STRT; I_1 = GRID%I_STOP
      J_0 = GRID%J_STRT; J_1 = GRID%J_STOP

C**** INITIALIZE TIME PARAMETERS
      NSTEP=(Itime-ItimeI)*NIdyn

C****
C**** INTEGRATE DYNAMIC TERMS (DIAGA AND DIAGB ARE CALLED FROM DYNAM)
C****
      CALL CHECKT ('DYNAM0')
         MODD5D=MOD(Itime-ItimeI,NDA5D)

         IF (MODD5D.EQ.0) IDACC(ia_d5d)=IDACC(ia_d5d)+1
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

      call startTimer('Atm. Dynamics')

#ifndef USE_FVCORE
      CALL DYNAM()
#else

      ! Using FV instead
        IF (MOD(Itime-ItimeI,NDAA).eq.0) CALL DIAGA0

      call Run(fvstate)

#ifndef CUBED_SPHERE
      CALL SDRAG (DTsrc)
#endif
        if (MOD(Itime-ItimeI,NDAA).eq.0) THEN
          call DIAGA
          call DIAGB
#ifdef CUBED_SPHERE
          call EPFLUX
#endif
        endif
#endif /* USE_FVCORE */

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
C**** Currently energy is put in uniformly weighted by mass
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)
#ifndef CUBED_SPHERE
      call COMPUTE_DYNAM_AIJ_DIAGNOSTICS(PUA, PVA, DT)
#endif
#ifdef SCM
       do L=1,LM
          CONV(I_TARG,J_TARG,L) = SG_CONV(L)
       enddo
#endif
      SD_CLOUDS(:,:,:) = CONV(:,:,:)
      call COMPUTE_WSAVE
C**** Scale WM mixing ratios to conserve liquid water
      DO L=1,LS1-1
      DO J=J_0,J_1
      DO I=I_0,I_1
        WM(I,J,L)=WM(I,J,L)* (PTOLD(I,J)/P(I,J))
      END DO
      END DO
      END DO
      CALL QDYNAM  ! Advection of Q by integrated fluxes
         CALL TIMER (NOW,MDYN)
#ifdef TRACERS_ON
      CALL TrDYNAM   ! tracer dynamics
#ifdef TRACERS_TOMAS
!TOMAS- This next section of code ratios the higher order moments of
!       aerosol mass to those of aerosol number so the distributions of
!       aerosol mass and number within a grid cell are consistent
      do n=1,NBINS
         call momentfix(IDTNUMD-1+n, IDTSO4-1+n)  !sulfate mass
         call momentfix(IDTNUMD-1+n, IDTNA -1+n)  !na+ mass
         call momentfix(IDTNUMD-1+n, IDTECOB-1+n) !hydrophobic EC
         call momentfix(IDTNUMD-1+n, IDTECIL-1+n)
         call momentfix(IDTNUMD-1+n, IDTOCOB-1+n)
         call momentfix(IDTNUMD-1+n, IDTOCIL-1+n)
         call momentfix(IDTNUMD-1+n, IDTDUST-1+n)
         call momentfix(IDTNUMD-1+n, IDTH2O-1+n)  !water mass
      enddo
#endif
#ifdef TRAC_ADV_CPU
         CALL TIMER (NOW,MTRADV)
#else
         CALL TIMER (NOW,MTRACE)
#endif
#endif

      call stopTimer('Atm. Dynamics')

C****
C**** Calculate tropopause level and pressure
C****
      CALL CALC_TROP
C**** calculate some dynamic variables for the PBL
#ifndef SCM
      CALL PGRAD_PBL
#endif
C**** calculate zenith angle for current time step
      CALL CALC_ZENITH_ANGLE

         CALL CHECKT ('DYNAM ')
         CALL TIMER (NOW,MSURF)
         IF (MODD5D.EQ.0) CALL DIAG5A (7,NIdyn)
         IF (MODD5D.EQ.0) CALL DIAGCA (2)
         IF (MOD(Itime,NDAY/2).eq.0) CALL DIAG7A
C****
C**** INTEGRATE SOURCE TERMS
C****

c calculate KE before atmospheric column physics
         call calc_kea_3d(kea)

#ifdef CUBED_SPHERE
c GWDRAG, SDRAG considered as column physics so that their KE
c dissipation gets included in the KE->PE adjustment
      CALL GWDRAG
      CALL SDRAG (DTsrc)
#endif

         IDACC(ia_src)=IDACC(ia_src)+1
         MODD5S=MOD(Itime-ItimeI,NDA5S)
         atmocn%MODD5S = MODD5S
         IF (MODD5S.EQ.0) IDACC(ia_d5s)=IDACC(ia_d5s)+1
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG5A (1,0)
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAGCA (1)

C**** FIRST CALL MELT_SI SO THAT TOO SMALL ICE FRACTIONS ARE REMOVED
C**** AND ICE FRACTION CAN THEN STAY CONSTANT UNTIL END OF TIMESTEP
! todo: move melt_si(ocean) to the end of the ocean driver, and
! possibly unite melt_si(lakes) with the rest of the lakes calls
      CALL MELT_SI(si_ocn,iceocn,atmocn,atmice)
      CALL MELT_SI(si_atm,icelak,atmocn,atmice)
      call seaice_to_atmgrid(atmice)
         CALL UPDTYPE
         CALL TIMER (NOW,MSURF)
C**** CONDENSATION, SUPER SATURATION AND MOIST CONVECTION
      CALL CONDSE
         CALL CHECKT ('CONDSE')
         CALL TIMER (NOW,MCNDS)
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)
         IF (MODD5S.EQ.0) CALL DIAGCA (3)

C**** RADIATION, SOLAR AND THERMAL
      MODRD=MOD(Itime-ItimeI,NRAD)
      CALL RADIA
         CALL CHECKT ('RADIA ')
         CALL TIMER (NOW,MRAD)
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)
         IF (MODD5S.EQ.0) CALL DIAGCA (4)

#ifdef TRACERS_ON
C**** Calculate non-interactive tracer surface sources and sinks
         call set_tracer_2Dsource
         CALL TIMER (NOW,MTRACE)

C****
C**** Add up the non-interactive tracer surface sources.
C****
      call sum_prescribed_tracer_2Dsources(dtsrc)
#endif

      call atm_phase1_exports

      return
      end subroutine atm_phase1

      subroutine atm_phase1_exports
! Copies fields calculated by the atmosphere into the data structures
! seen by physics of the surface components (ocean, ice, land).
! Some fields are already type-classified (e.g. radiative fluxes).
! A per-type breakdown will soon be applied to other fields as
! appropriate (e.g. fields which vary with height will have
! different values over the ocean and the other portions of a gridbox).
! Some fields have already been stored in atmsrf%xxx and are not
! referenced here in Step 1.  For temporary convenience,
! fields depending on surface pressure are being copied into the
! per-surface-type structures by CALC_AMPK because subroutine FILTER
! is currently being called after the main surface physics, but
! before special "daily" surface coding which sometimes requires
! surface pressure - once FILTER is absorbed into DYNAM this hack
! can be eliminated.
      use fluxes, only : atmsrf,asflx4,prec,eprec
#ifdef TRACERS_ON
#ifdef TRACERS_WATER
      use fluxes, only : trprec
#endif
#endif
      use rad_com, only : trhr,fsf,trsurf,cosz1
#ifdef OBIO_RAD_coupling
      use fluxes, only : atmocn
      use rad_com, only : dirvis,fsrdif,dirnir,difnir
#endif
      implicit none
      integer :: it

      ! Step 1: copy fields not already stored in atmsrf%xxx
      atmsrf%prec = prec
      atmsrf%eprec = eprec
#ifdef TRACERS_WATER
      atmsrf%trprec = trprec
#endif
      atmsrf%cosz1 = cosz1
      atmsrf%flong = trhr(0,:,:)


      ! Copy from gridbox-mean data structure to per-type structures
      do it=1,4
        asflx4(it)%atm_exports_phase1 = atmsrf%atm_exports_phase1
#ifdef TRACERS_WATER
        asflx4(it)%tratm_exports_phase1 = atmsrf%tratm_exports_phase1
#endif
      enddo

      ! fields already available per type
      do it=1,4
        asflx4(it)%fshort = fsf(it,:,:)
        asflx4(it)%trup_in_rad = trsurf(it,:,:)
      enddo


      ! miscellaneous fields only defined for some types
#ifdef OBIO_RAD_coupling
      atmocn % DIRVIS = DIRVIS
      atmocn % DIFVIS = FSRDIF
      atmocn % DIRNIR = DIRNIR
      atmocn % DIFNIR = DIFNIR
#endif

      return
      end subroutine atm_phase1_exports

      subroutine atm_exports_phasesrf
! Copies fields calculated by the atmosphere into the data structures
! seen by the physics of the surface components (ocean, ice, land).
! Per-component exports will be introduced as appropriate (e.g. fields
! which vary with height will have different values over the ocean and
! the other portions of a gridbox).
      use fluxes, only : atmsrf,asflx4
      implicit none

      integer :: it

      call get_atm_layer1

      do it=1,4
        asflx4(it)%atm_exports_phasesrf = atmsrf%atm_exports_phasesrf
#ifdef TRACERS_ON
        asflx4(it)%tratm_exports_phasesrf =
     &       atmsrf%tratm_exports_phasesrf
#endif
      enddo

      return
      end subroutine atm_exports_phasesrf

      subroutine get_atm_layer1
C**** Copies first-layer atm. conditions into the 2D arrays
C**** in the atm-surf. coupling data structure.
      use fluxes, only : atmsrf
      use domain_decomp_atm, only : grid, getDomainBounds
      use atm_com, only : t,q,ualij,valij
      use geom, only : imaxj,byaxyp
#ifdef TRACERS_ON
      use tracer_com, only : ntm,trm
#endif
      implicit none
      integer :: n,i,j,i_0,i_1,j_0,j_1
c
      call getDomainBounds(grid, i_strt=i_0,i_stop=i_1,
     &  j_strt=j_0,j_stop=j_1)
c
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        atmsrf%temp1(i,j) = t(i,j,1)
        atmsrf%q1(i,j) = q(i,j,1)
        atmsrf%u1(i,j) = ualij(1,i,j)
        atmsrf%v1(i,j) = valij(1,i,j)
      enddo
      enddo
#ifdef TRACERS_ON
      do n=1,ntm
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        atmsrf%trm1(n,i,j) = trm(i,j,1,n)*byaxyp(i,j)
      enddo
      enddo
      enddo
#endif
      end subroutine get_atm_layer1

      subroutine atm_phase2
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      use resolution, only : lm
      USE MODEL_COM
      USE DYNAMICS, only : nidyn,nfiltr,mfiltr
      USE GETTIME_MOD
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: mtrace
#endif
      USE DIAG_COM, only : kvflxo,oa,koa,ia_filt
     &     ,MODD5S,NDAa, NDA5d,NDA5s,NDA4
      USE SUBDAILY, only : nsubdd,get_subdd,accSubdd
#ifndef CUBED_SPHERE
      USE ATMDYN, only : FILTER
      USE ATM_COM, only : P
      USE RESOLUTION, only : PTOP
#endif
#ifdef SCM
      USE SCMCOM , only : SG_CONV,SCM_SAVE_T,SCM_SAVE_Q,
     &    iu_scm_prt,iu_scm_diag
#endif
      USE FLUXES, only : atmice
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use SystemTimers_mod
      implicit none

      REAL*8 start,now

      call seaice_to_atmgrid(atmice)
      CALL ADVSI_DIAG ! needed to update qflux model, dummy otherwise
C**** SAVE some noon GMT ice quantities
      IF (MOD(Itime+1,NDAY).ne.0 .and. MOD(Itime+1,NDAY/2).eq.0)
     &        call vflx_OCEAN

C**** IF ATURB is used in rundeck then this is a dummy call
C**** CALCULATE DRY CONVECTION ABOVE PBL
      CALL ATM_DIFFUS (2,LM-1,dtsrc)
         CALL CHECKT ('DRYCNV')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (9)

C**** UPDATE DIAGNOSTIC TYPES
         CALL UPDTYPE
C**** ADD DISSIPATED KE FROM COLUMN PHYSICS CALCULATION BACK AS LOCAL HEAT
      CALL DISSIP ! uses kea calculated before column physics
         CALL CHECKT ('DISSIP')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (7)
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)

#ifdef CUBED_SPHERE
      IDACC(ia_filt)=IDACC(ia_filt)+1 ! prevent /0
#else
C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(Itime-ItimeI,NFILTR).EQ.0) THEN
           IDACC(ia_filt)=IDACC(ia_filt)+1
           IF (MODD5S.NE.0) CALL DIAG5A (1,0)
           CALL DIAGCA (1)
           CALL FILTER
           CALL CHECKT ('FILTER')
           CALL TIMER (NOW,MDYN)
           CALL DIAG5A (14,NFILTR*NIdyn)
           CALL DIAGCA (8)
      END IF
#endif
#ifdef TRACERS_ON
#ifdef CUBED_SPHERE
! Reinitialize instantaneous consrv qtys (every timestep since
! DIAGTCA is called every timestep for 3D sources)
      CALL DIAGCA (1) ! was not called w/ SLP filter
#endif
C**** 3D Tracer sources and sinks
C**** Tracer gravitational settling for aerosols
      CALL TRGRAV
C**** Tracer radioactive decay (and possible source)
      CALL TDECAY
C**** Calculate 3D tracers sources and sinks

      call tracer_3Dsource
C**** Accumulate tracer distribution diagnostics
      CALL TRACEA
         CALL TIMER (NOW,MTRACE)
         CALL CHECKT ('T3DSRC')
#endif
C****
C**** WRITE SUB-DAILY DIAGNOSTICS EVERY NSUBDD hours
C****
      if (Nsubdd.ne.0) then
        call accSubdd
        if (mod(Itime+1,Nsubdd).eq.0) call get_subdd
      end if
#ifdef TRACERS_DUST
      call ahourly
#endif

#ifdef SCM
c*****call scm diagnostics every time step
      call scm_diag
#endif

      IF (MOD(Itime+1-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 E-history

      IF (Kvflxo.EQ.0.) OA(:,:,4:KOA)=0. ! to prepare for future saves

      return
      end subroutine atm_phase2

      SUBROUTINE INPUT_atm (istart,istart_fixup,do_IC_fixups,
     &     is_coldstart,KDISK_restart,IRANDI)

C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE Dictionary_mod
      USE CONSTANT, only : grav
      USE FLUXES, only : nisurf,atmocn,atmice
      USE RESOLUTION, only : ls1,plbot
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only :
     *      irand,idacc ,nday,dtsrc ,iyear1,itime,itimei,itimee
     *     ,mdyn,mcnds,mrad,msurf,mdiag
      USE DIAG_ZONAL, only : imlon
      USE RANDOM
      USE DYNAMICS, only : USE_UNR_DRAG
      USE ATM_COM, only : ij_debug
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE SOMTQ_COM, only : mz,tmom
      USE ATM_COM, only : p,t
      USE TRACER_COM,only: MTRACE,daily_z
#ifdef TRACERS_SPECIAL_Shindell
     *     ,mchem
#endif
#ifdef TRAC_ADV_CPU
      USE TRACER_COM,only: MTRADV
#endif
#endif
      USE SOIL_DRV, only: init_LSM,daily_earth
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, AM_I_ROOT
#ifndef CUBED_SPHERE
      USE ATMDYN, only : init_ATMDYN
#endif
#ifdef IRRIGATION_ON
      use irrigate_crop, only : init_irrigate
#endif
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: fvstate,initialize
#endif
      IMPLICIT NONE
!@var istart start(1-8)/restart(>8)  option
      integer :: istart,istart_fixup,do_IC_fixups
      LOGICAL :: is_coldstart
      INTEGER :: KDISK_restart
      INTEGER :: IRANDI

      INTEGER I,J,L,K

!@nlparam HOURI,DATEI,MONTHI,YEARI        start of model run
!@nlparam TIMEE,HOURE,DATEE,MONTHE,YEARE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
c      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
c     *    TIMEE=-1,HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHOURE=-1

      LOGICAL :: redoGH, iniPBL, inilake, iniSNOW
      INTEGER :: J_0, J_1

#ifdef USE_ESMF
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      write(6,*) 'mpi-zone',J_0,' - ',J_1
#endif

      CALL SET_TIMER("ATMOS. DYNAM",MDYN)
      CALL SET_TIMER("CONDENSATION",MCNDS)
      CALL SET_TIMER("   RADIATION",MRAD)
      CALL SET_TIMER("     SURFACE",MSURF)
      CALL SET_TIMER(" DIAGNOSTICS",MDIAG)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      CALL SET_TIMER("     TRACERS",MTRACE)
#ifdef TRAC_ADV_CPU
      CALL SET_TIMER(" TRACER ADV.",MTRADV)
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      CALL SET_TIMER("   CHEMISTRY",MCHEM)
#endif


      call sync_param( "ij_debug",ij_debug , 2)

C****
C**** Set some documentary parameters in the database
C****
      call set_param("IM",IM,'o')
      call set_param("JM",JM,'o')
      call set_param("LM",LM,'o')
      call set_param("LS1",LS1,'o')
      call set_param("PLBOT",Plbot,LM+1,'o')

      if(istart.eq.2) call read_nmc()

C****
C**** IRANDI seed for random perturbation of current state (if/=0)
C****        tropospheric temperatures are changed by at most 1 degree C
      IF (ISTART.LT.10 .AND. IRANDI.NE.0) THEN
        CALL RINIT (IRANDI)
        CALL PERTURB_TEMPS
        IF (AM_I_ROOT())
     *       WRITE(6,*) 'Initial conditions were perturbed !!',IRANDI
      END IF

      CALL CALC_AMPK(LM)

#ifdef TRACERS_ON
      if(istart.le.2) then
        call COMPUTE_GZ(p,t,tmom(mz,:,:,:),daily_z)
        daily_z = daily_z/grav
      endif
#endif

C****
      CALL RINIT (IRAND)
c Note on FFT initialization: IMLON is defined by the diag_zonal module,
c not by the resolution module.  IMLON==IM for a latlon grid.
      CALL FFT0 (IMLON)  ! CALL FFT0(IM)
      CALL init_QUS(grid,im,jm,lm)
#ifndef CUBED_SPHERE
      CALL init_ATMDYN
#endif
      call init_sdrag
C**** Initialize nudging
#ifdef NUDGE_ON
      CALL NUDGE_INIT
#endif
C****
C**** Initialize the gravity wave drag scheme
C****
      CALL init_GWDRAG
      call sync_param( "USE_UNR_DRAG", USE_UNR_DRAG )
#ifndef SCM
#ifndef CUBED_SPHERE
      if (USE_UNR_DRAG==1) CALL init_UNRDRAG
#endif
#endif

#ifdef SCM
!      read scm data and initialize model
!      note:  usavg,vsavg and wsavg filled from here
       call init_scmdata
#endif

      CALL init_CLD(istart)
      CALL init_RAD(istart)
      CALL daily_orbit(.false.)             ! not end_of_day
      CALL daily_ch4ox(.false.)             ! not end_of_day
      CALL daily_RAD(.false.)

      call atm_phase1_exports

C**** Initialize lake variables (including river directions)
      if(istart.eq.2) call read_agrice_ic
      iniLAKE = is_coldstart
      CALL init_lakeice(inilake,do_IC_fixups)
c      call stop_model('set_noice_defaults prob that fwater==flake'//
c     &     ' not initialized yet?',255)
      call seaice_to_atmgrid(atmice) ! set gtemp etc.
      CALL init_LAKES(inilake,istart_fixup)

C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS (INCL. VEGETATION)
C**** Recompute Ground hydrology data if redoGH (new soils data)
C****
      if(istart.eq.2) call read_landsurf_ic
      iniSNOW = is_coldstart       ! extract snow data from first soil layer
      redoGH = .false.
      CALL init_LSM(DTsrc/NIsurf,redoGH,iniSNOW,inilake,ISTART)
#ifdef IRRIGATION_ON
      call init_irrigate()
#endif
      CALL daily_EARTH(.false.)            ! not end_of_day

#ifdef CALCULATE_FLAMMABILITY
      CALL init_flammability
#endif

C**** Initialize land ice (must come after oceans)
      if(istart.eq.2) call read_landice_ic
      CALL init_LI(istart_fixup)

C**** Initialize pbl (and read in file containing roughness length data)
      iniPBL = is_coldstart
#ifndef CUBED_SPHERE /* until a better solution is found */
      if (iniPBL) call recalc_agrid_uv   ! PBL needs A-grid winds
#endif
      CALL init_pbl(iniPBL,istart)

! note there is an all-component call to reset_diag in init_diag
      CALL init_DIAG
      CALL UPDTYPE   ! for atm-grid diags

! init tasks transplanted from main program.  cannot disperse to
! corresponding components w/o changing results.
! DAILY_atmdyn adjusts the global mean mass when ISTART=2, so
! moving it before init of other atm components will change results
      CALL DAILY_atmdyn(.false.)           ! not end_of_day
#ifdef USE_FVCORE
C****
C**** Initialize FV dynamical core (ESMF component) if requested
C**** For restarts/continuations, FV state and import files are
C**** copied to the appropriate names by this procedure, and for
C**** cold starts the required IC files are generated.
C****
      Call Initialize(fvstate, istart, kdisk_restart)
#endif

! this daily_OCEAN call belongs in ocean init, but
! ISTART=2 result differs if daily_OCEAN comes before init_pbl,
! since daily_OCEAN replaces the GIC values of gtemp
! with the values from the prescribed SST file
      CALL daily_OCEAN(.false.,atmocn)            ! not end_of_day
! need to do this after prescribed-ice daily_OCEAN changes si frac.
! not needed once daily_OCEAN is moved before atm init.
      call seaice_to_atmgrid(atmice) ! debug
      CALL UPDTYPE

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      CALL daily_tracer(.false.)
#endif
      CALL CHECKT ('INPUT ')

      if (AM_I_ROOT()) then
         WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
         WRITE (6,*) "PLbot=",PLbot
      end if

      return

      end subroutine INPUT_atm

      subroutine alloc_drv_atm()
#ifdef SCM
      USE SCMCOM, only : I_TARG,J_TARG
      use Dictionary_mod, only : sync_param
#endif
c Driver to allocate arrays that become dynamic as a result of
c set-up for MPI implementation
      USE DOMAIN_DECOMP_ATM, ONLY : grid,init_grid
#ifdef GLINT2
      USE DOMAIN_DECOMP_ATM, ONLY : glint2
      use MpiSupport_mod, only: ROOT_PROCESS
      Use glint2_modele
#endif
      USE RESOLUTION, only : im,jm,lm
#ifndef CUBED_SPHERE
      USE MOMENTS, only : initMoments
#endif
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      use TRACER_COM, only: initTracerCom, alloc_tracer_com
      use ghy_tracers, only: initGhyTracers
#endif
      IMPLICIT NONE
      include 'mpif.h'      ! Needed for GLINT2

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call initTracerCom
      call initGhyTracers
#endif

#ifdef SCM
!TODO push init_grid SCM option down into INIT_GRID.
      call sync_param( "J_TARG", J_TARG )
      call init_grid(grid, im, jm, lm, j_scm=j_targ)
#else
c initialize the atmospheric domain decomposition
c for now, CREATE_CAP is only relevant to the cubed sphere grid
      call init_grid(grid, im, jm, lm, CREATE_CAP=.true.)
#endif

#ifdef GLINT2
      glint2 = glint2_modele_new('GLINT2', 6, 'm', 1,
     &    im, jm,
     &    grid%i_strt_halo, grid%i_stop_halo,
     &    grid%j_strt_halo, grid%j_stop_halo,
     &    grid%i_strt, grid%i_stop, grid%j_strt, grid%j_stop,
     &    grid%j_strt_skp, grid%j_stop_skp,
     &    MPI_COMM_WORLD, ROOT_PROCESS)
#endif  ! GLINT2

      call alloc_dynamics(grid)
      call alloc_atm_com(grid)
      call alloc_smomtq(grid)

      call alloc_fluxes !(grid)
      call alloc_clouds_com(grid)
      call alloc_ghy_com(grid)
      call alloc_pbl_com(grid)
      call alloc_diag_com(grid)
      call alloc_diag_loc(grid)
      call alloc_strat_com(grid)
      call alloc_rad_com(grid)
#ifdef RAD_O3_GCM_HRES
      call alloc_RAD_native_O3(grid)
#endif
      call alloc_lakes(grid)
      call alloc_lakes_com(grid)
      call alloc_landice_com(grid)

#ifdef CALCULATE_FLAMMABILITY
      call alloc_flammability(grid)
#endif
#ifdef BIOGENIC_EMISSIONS
      call alloc_biogenic_emis(grid)
#endif
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      call alloc_tracer_com(grid)
#ifdef TRACERS_DRYDEP
      call alloc_trdrydep(grid)
#endif
#ifdef TRACERS_SPECIAL_Lerner
      call alloc_tracer_special_lerner_com(grid)
      call alloc_linoz_chem_com(grid)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call alloc_trchem_shindell_com(grid)
      call alloc_tracer_sources(grid)
      call alloc_lightning(grid)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      call alloc_aerosol_sources(grid)
#endif
#ifdef TRACERS_AMP
      call alloc_tracer_amp_com(grid)
#endif
#ifdef TRACERS_TOMAS
      call alloc_tracer_tomas_com(grid)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      CALL alloc_dust(grid)
#endif
#endif
      call alloc_tracer_adv(grid)
#ifdef USE_ENT
!!! should be done in init_module_ent
      call alloc_ent_com(grid)
#else
      call alloc_veg_com(grid)
#endif
#ifdef TRACERS_ON
      call alloc_trdiag_com
#ifdef TRACERS_SPECIAL_Shindell
      call interpolateAltitude()
#endif
#endif
#ifdef NUDGE_ON
      call alloc_nudge(grid)
#endif

#ifdef SCM 
      call alloc_scm_com()
#endif
#ifndef CUBED_SPHERE
      call initMoments
#endif
      end subroutine alloc_drv_atm

      subroutine def_rsf_atmvars(fid)
!@sum  def_rsf_atmvars defines atm prognostic array structure in rsf
!@auth M. Kelley
!@ver  beta
      implicit none
      integer :: fid
      call def_rsf_atm    (fid)
      call def_rsf_lakes  (fid)
      call def_rsf_agrice (fid)
      call def_rsf_icedyn (fid)
      call def_rsf_earth  (fid)
      call def_rsf_soils  (fid)
      call def_rsf_vegetation(fid)
#ifdef USE_ENT
      call def_rsf_veg_related(fid)
#endif
      call def_rsf_snow   (fid)
      call def_rsf_landice(fid)
      call def_rsf_bldat  (fid)
      call def_rsf_pbl    (fid)
      call def_rsf_clouds (fid)
      call def_rsf_somtq  (fid)
      call def_rsf_rad    (fid)
#ifdef CALCULATE_FLAMMABILITY
      call def_rsf_flammability(fid)
#endif
#ifdef TRACERS_ON
      call def_rsf_tracer (fid)
#endif
      call def_rsf_subdd  (fid)
      call def_rsf_fluxes (fid)
      return
      end subroutine def_rsf_atmvars

      subroutine new_io_atmvars(fid,iorw)
      implicit none
      integer, intent(in) :: fid,iorw
      call new_io_atm    (fid,iorw)
      call new_io_lakes  (fid,iorw)
      call new_io_agrice (fid,iorw)
      call new_io_earth  (fid,iorw)
      call new_io_soils  (fid,iorw)
      call new_io_vegetation  (fid,iorw)
#ifdef USE_ENT
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
      call new_io_veg_related(fid,iorw)
        !call io_ent    (kunit,iaction,ioerr) ! io_vegetation handles ent
#endif
      call new_io_snow   (fid,iorw)
      call new_io_landice(fid,iorw)
      call new_io_bldat  (fid,iorw)
      call new_io_pbl    (fid,iorw)
      call new_io_clouds (fid,iorw)
      call new_io_somtq  (fid,iorw)
      call new_io_rad    (fid,iorw)
      call new_io_icedyn (fid,iorw)
#ifdef CALCULATE_FLAMMABILITY
      call new_io_flammability(fid,iorw)
#endif
#ifdef TRACERS_ON
      call new_io_tracer (fid,iorw)
#endif
      call new_io_subdd  (fid,iorw)
      call new_io_fluxes (fid,iorw)
      return
      end subroutine new_io_atmvars

      subroutine daily_atm(end_of_day)
      use filemanager
      use MODEL_COM, only: nday,itime
      use DYNAMICS, only : nidyn
      USE SOIL_DRV, only: daily_earth
      use diag_com, only : kvflxo,iu_vflxo,oa,koa
      use domain_decomp_atm, only: grid,writei8_parallel
      implicit none
      logical, intent(in) :: end_of_day ! not used yet

      call DIAG5A (1,0)
      call DIAGCA (1)
      CALL daily_atmdyn(.true.)  ! end_of_day
      CALL daily_orbit(.true.)   ! end_of_day
      CALL daily_ch4ox(.true.)   ! end_of_day
      call daily_RAD(.true.)
        
      call daily_LAKE
      call daily_EARTH(.true.)  ! end_of_day
        
      call daily_LI
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call daily_tracer(.true.)
#endif
      call CHECKT ('DAILY ')
      call DIAG5A (16,NDAY*NIdyn)
      call DIAGCA (10)
      call sys_flush(6)
      call UPDTYPE

C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY 24 HOURS
C****
      IF (Kvflxo.NE.0.) THEN
        call writei8_parallel(grid,iu_vflxo,nameunit(iu_vflxo),oa,Itime)
C**** ZERO OUT INTEGRATED QUANTITIES
        OA(:,:,4:KOA)=0.
      END IF

      return
      end subroutine daily_atm

      subroutine finalize_atm
      USE SUBDAILY, only : close_subdd
#ifdef USE_FVCORE
      USE MODEL_COM, only : kdisk
      USE FV_INTERFACE_MOD, only: fvstate
      USE FV_INTERFACE_MOD, only: Finalize
#endif
#ifdef SCM
      USE FILEMANAGER, only : closeunit
      USE SCMCOM , only : iu_scm_prt,iu_scm_diag
#endif
      implicit none
#ifdef USE_FVCORE
         call Finalize(fvstate, kdisk)
#endif

C**** CLOSE SUBDAILY OUTPUT FILES
      CALL CLOSE_SUBDD

#ifdef SCM
      call closeunit(iu_scm_prt)
      call closeunit(iu_scm_diag)
#endif
      return
      end subroutine finalize_atm

      SUBROUTINE CHECKT (SUBR)
!@sum  CHECKT Checks arrays for NaN/INF and reasonablness
!@auth Original Development Team

C**** CHECKT IS TURNED ON BY SETTING QCHECK=.TRUE. IN NAMELIST
C**** REMEMBER TO SET QCHECK BACK TO .FALSE. AFTER THE ERRORS ARE
C**** CORRECTED.
      USE CONSTANT, only : tf
      USE RESOLUTION, only : ls1
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v,t,p,q,wm,pk
#ifdef BLK_2MOM
      USE ATM_COM, only : wmice
#endif
      USE MODEL_COM
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, AM_I_ROOT
      USE soil_drv, only : checke
      IMPLICIT NONE
      INTEGER I,J,L
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H, J_1H, I_0,I_1, I_0H,I_1H, njpol
      INTEGER :: I_0STG,I_1STG,J_0STG,J_1STG
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     *     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
c      I_0STG = grid%I_STRT_STGR
c      I_1STG = grid%I_STOP_STGR
      I_0STG = I_0
      I_1STG = I_1
      J_0STG = grid%J_STRT_STGR
      J_1STG = grid%J_STOP_STGR
      njpol = grid%J_STRT_SKP-grid%J_STRT

      IF (QCHECK) THEN
C**** Check all prog. arrays for Non-numbers
        CALL CHECK3B(U(I_0STG:I_1STG,J_0STG:J_1STG,:),
     &       I_0STG,I_1STG,J_0STG,J_1STG,0,LM,SUBR,'u     ')
        CALL CHECK3B(V(I_0STG:I_1STG,J_0STG:J_1STG,:),
     &       I_0STG,I_1STG,J_0STG,J_1STG,0,LM,SUBR,'v     ')
        CALL CHECK3B(T(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'t     ')
        CALL CHECK3B(Q(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'q     ')
        CALL CHECK3B(P(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &       SUBR,'p     ')
        CALL CHECK3B(WM(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'wm    ')
#ifdef BLK_2MOM
        CALL CHECK3B(WMICE(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'wmice    ')
#endif

        DO J=J_0,J_1
        DO I=I_0,I_1
          IF (Q(I,J,1).gt.1d-1)print*,SUBR," Q BIG ",i,j,Q(I,J,1:LS1)
          IF (T(I,J,1)*PK(1,I,J)-TF.gt.50.) print*,SUBR," T BIG ",i,j
     *         ,T(I,J,1:LS1)*PK(1:LS1,I,J)-TF
        END DO
        END DO
        DO L=1,LM
        DO J=J_0,J_1
        DO I=I_0,I_1
          IF (Q(I,J,L).lt.0.) then
            print*,"After ",SUBR," Q < 0 ",i,j,Q(I,J,L)
            call stop_model('Q<0 in CHECKT',255)
          END IF
          IF (WM(I,J,L).lt.0.) then
            print*,"After ",SUBR," WM < 0 ",i,j,WM(I,J,L)
            call stop_model('WM<0 in CHECKT',255)
          END IF
#ifdef BLK_2MOM
          IF (WMICE(I,J,L).lt.0.) then
            print*,"After ",SUBR," WMICE < 0 ",i,j,WMICE(I,J,L)
            call stop_model('WMICE<0 in CHECKT',255)
          END IF
#endif
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
C**** Check Land Ice arrays
        CALL CHECKLI(SUBR)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** check tracers
        CALL CHECKTR(SUBR)
#endif
      END IF

      RETURN
      END SUBROUTINE CHECKT
