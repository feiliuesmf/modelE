
C****   
C**** SURFACE.f    SURFACE fluxes    2006/12/21
C****
#include "rundeck_opts.h"

! Overall structure:
! 
! 1. Init
! 2. Loop over grid points (some diagnostics here)
! 3. Diagnostics (in separate loop over gridpoints; composite diagnostics)
! 4. Now that all the surface fluxes have been calculated... call atm_diffus()
!    to update atmospheric state
! 5. Other diagnostics
! 6. Things that were moved from the main loop to the end
!    (Each one has its own loop inside).
! 
! 



      SUBROUTINE SURFACE
!@sum SURFACE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@vers 2013/03/27
!@auth Nobody will claim responsibilty

      USE CONSTANT, only : rgas,lhm,lhe,lhs
     *     ,sha,tf,rhow,shv,shi,stbo,bygrav,by6
     *     ,deltx,teeny,rhows,grav
      USE ATM_COM, only : t,q
      USE MODEL_COM, only : modelEclock
      USE MODEL_COM, only : dtsrc,idacc,nday,itime,qcheck
      use TimeConstants_mod, only: INT_HOURS_PER_DAY
#ifdef SCM
      USE SCMDIAG, only : EVPFLX,SHFLX
      USE SCMCOM, only : iu_scm_prt, ALH, ASH, SCM_SURFACE_FLAG
     &     ,I_TARG,J_TARG
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, GLOBALSUM
      USE GEOM, only : axyp,imaxj,byaxyp
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE ATM_COM, only : pmid,pk,pedn,pek,MA,byMA
      USE RAD_COM, only : trhr,fsf,cosz1,trsurf
#ifdef TRACERS_ON
      use OldTracer_mod, only: itime_tr0, needtrs
      USE TRACER_COM, only : NTM,trm,trmom,
     *     n_Be7, n_Be10
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: dodrydep
#endif
#ifdef TRACERS_WATER
      use OldTracer_mod, only: nWATER,tr_wd_TYPE
#endif
#endif
      USE SOCPBL, only : npbl=>n
      USE PBL_DRV, only : alloc_pbl_args, dealloc_pbl_args
      USE PBL_DRV, only : pbl, t_pbl_args, xdelt
      USE DIAG_COM, only : MODD5S
      USE DIAG_COM, only : ndasf,ia_srf
     &     ,aij=>aij_loc,ij_dskin  ! temporarily still here
      USE SEAICE, only : xsi,ace1i,alami0,rhoi,byrls,solar_ice_frac
     *     ,tfrez,dEidTi,alami
      USE SEAICE_COM, only : si_atm
      USE LAKES_COM, only : mwl,gml,flake,icelak
      USE LAKES, only : minmld
      USE GHY_COM, only : FEARTH
      USE EXCHANGE_TYPES, only :
     &     avg_patches_pbl_exports
     &    ,avg_patches_srfflx_exports
     &    ,avg_patches_srfstate_exports
     &    ,avg_patches_srfflx_exports_gla
      USE FLUXES, only :
     *      uflux1,vflux1,tflux1,qflux1
     *     ,UOdrag
     &     ,after_atm_phase1,during_srfflx
     &     ,nisurf,fland,flice,focean
     &     ,atmocn,atmice,atmgla,atmlnd,asflx,atmsrf
     &     ,atmglas
#ifdef GLINT2
      USE FLUXES, only : atmglas_hp
      use hp2hc
#endif

#ifdef TRACERS_ON
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     &     ,dust_flux_glob,dust_flux2_glob
#ifdef TRACERS_DRYDEP
     &     ,depo_turb_glob,depo_grav_glob
#endif
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7D_acc
#endif
#ifndef SKIP_TRACER_DIAGS
#ifdef TRACERS_DRYDEP
      USE TRDIAG_COM, only : taijn=>taijn_loc,
     *      tij_drydep, tij_gsdep, itcon_dd
#endif
#endif /*SKIP_TRACER_DIAGS*/
#ifdef TRACERS_ON
      use trdiag_com, only: trcsurf,trcSurfByVol
#endif
#endif

      USE SOIL_DRV, only: earth,ground_e

      USE Timer_mod, only: Timer_type
      USE TimerList_mod, only: startTimer => start
      USE TimerList_mod, only: stopTimer => stop
      USE itype_enum     ! Surface Type enumeration: ITYPE_OCEAN, etc.

      IMPLICIT NONE

      ! ----------- For Debugging
      LOGICAL :: debug       ! Set to true for certain conditions where you want printout (eg, i=17 and j=34)
      ! ----------- End Debugging Variables

      INTEGER I,J,K,KR,JR,NS,NSTEPS,MODDSF,MODDD,ITYPE,IH,IHM
     *     ,ii,itype4,ipatch
#ifdef GLINT2
      INTEGER ihp
#endif
      REAL*8 PLAND,PLICE,POICE,POCEAN,PS,P1K
     *     ,ELHX,MSI2,CDTERM,CDENOM,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG
     *     ,dT2,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER
     *     ,PSK,Q1,THV1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,RCDHDWS,RCDQDWS
     *     ,SHEAT,TRHEAT,T2DEN,T2CON,T2MUL,FQEVAP,Z1BY6L,EVAPLIM,F2
     *     ,FSRI(2),HTLIM,dlwdt,byNIsurf,TGO,SRHDT
     *     ,PEARTH

      REAL*8 MA1, MSI1

!@var XXX_sv temporary arrays to store the values of fluxes for
!@+   the current surface timestep for a few flux diagnostics that
!@+   are not averages over the full physics timestep.
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     SHDT_sv,EVHDT_sv,SRHDT_sv,TRHDT_sv,EVAP_sv,TRHDT_sv2
      REAL*8, DIMENSION(2,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                    GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     TGRND,TGRN2,TGR4,E1 ! 2 for ocean,seaice

      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8 QSAT,DQSATDT,TR4
c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 qg_sat,dtsurf,uocean,vocean,qsrf,us,vs,ws,ws0,
     &     dmua_ij,dmva_ij
c
      logical :: lim_lake_evap,lim_dew ! for tracer convenience
#ifdef TRACERS_ON
      integer n,nx
      integer ntix(ntm), ntx
      real*8, dimension(ntm) :: trgrnd,trgrnd2
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: tevaplim
#endif
#ifdef TRACERS_DRYDEP
      real*8 tdryd, tdd, td1, rtsdt, rts, depvel, gsvel
#endif
#endif
C****
      INTEGER :: J_0, J_1, I_0,I_1

      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,SSS
      REAL*8, DIMENSION(:,:,:), POINTER :: SSI
#ifdef TRACERS_GASEXCH_ocean
      REAL*8, DIMENSION(:,:,:), POINTER :: TRGASEX
#endif

      type (Timer_type), pointer :: aTimer

      RSI => SI_ATM%RSI
      MSI => SI_ATM%MSI
      SNOWI => SI_ATM%SNOWI
      SSI => SI_ATM%SSI
      SSS => atmocn%SSS
#ifdef TRACERS_GASEXCH_ocean
      TRGASEX => atmocn%TRGASEX
#endif

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL getDomainBounds(grid, J_STRT=J_0,        J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C****

      call startTimer('SURFACE()')

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      byNIsurf=1.d0/real(NIsurf)
      IH=modelEclock%hour()+1
      IHM = IH+(modelEclock%date()-1)*INT_HOURS_PER_DAY

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        PLAND=FLAND(I,J)
        PWATER=1.-PLAND
      ! RSI = RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
        POICE=RSI(I,J)*PWATER
        POCEAN=PWATER-POICE
        atmocn%ftype(i,j) = pocean
        atmice%ftype(i,j) = poice
        atmgla%ftype(i,j) = flice(I,J)
        atmlnd%ftype(i,j) = fearth(i,j)
      enddo
      enddo

      call downscale_pressure_li

      do ipatch=1,size(asflx)
        ! solar,trheat,evapor,sensht,latht,e0,dmua,dmva,runo,eruno
        asflx(ipatch)%srfflx_exports(:,:,:) = 0.
#ifdef TRACERS_ON
        ! trdrydep,trevapor,trsrfflx,truno
        asflx(ipatch)%trsrfflx_exports(:,:,:,:) = 0.
#endif
      enddo

      CALL PRECIP_SI(si_atm,icelak,atmice)

#ifdef GLINT2
      do ihp=1,ubound(atmglas_hp,1)
        CALL PRECIP_LI(atmglas_hp(ihp),ihp)
      enddo
      call bundle_hp_to_hc(bundle_precip_li,
     &       atmglas_hp(1:), atmglas(1:))
#else
      do ipatch=1,ubound(atmglas,1)
        CALL PRECIP_LI(atmglas(ipatch),ipatch)
      enddo
#endif

c The following landice section will migrate once precip_li AND precip_lk
c are folded into their respective top-level drivers in a "single-entry-point"
c design (as already exists for the land surface).
c There is a transfer of precip-induced runoff from land ice to lakes,
c and the lakes need the aggregate over patches.
      call avg_patches_srfflx_exports(grid,
c     &     atmglas,atmgla, ! gfortran prob. if passed as class() args
     &     atmglas(:)%atmsrf_xchng_vars,atmgla%atmsrf_xchng_vars,
     &     rel=.true.)
      call avg_patches_srfflx_exports_gla(grid,atmglas,atmgla,
     &     rel=.true.)
      call surface_diag_post_precip_li


#ifdef IRRIGATION_ON
C**** CHECK FOR IRRIGATION POSSIBILITY
      CALL IRRIG_LK
#endif
      CALL PRECIP_LK
         CALL CHECKT ('PRECIP')

      call seaice_to_atmgrid(atmice)

c avoid uninitialized variable problems when the first gridpoint
c in the domain is ocean
      SNOW = 0.

C**** INITIALIZE TGRND: THIS IS USED TO UPDATE T OVER SURFACE STEPS
      DO J=J_0,J_1
      DO I=I_0,I_1
        ! GTEMP temperature of surface (C)
        ! GTEMP2 "ground" temperature of "second" layer (C)
        ! GTEMPR radiative ground temperature over surface type (K)
        TGRND(2,I,J)=atmice%GTEMP(I,J)
        TGRN2(2,I,J)=atmice%GTEMP2(I,J)
        TGR4(2,I,J)=atmice%GTEMPR(I,J)**4
      END DO 
      END DO 

C**** Zero out fluxes summed over type and surface time step

      ! itype == 1, ocean; 2, ocean ice; 3, land ice; 4, land

      E1=0.
#ifdef SCM
      EVPFLX= 0.0d0
      SHFLX = 0.0d0
#endif
#ifdef TRACERS_GASEXCH_ocean
      TRGASEX = 0.0d0
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) 
      dust_flux_glob = 0.d0
      dust_flux2_glob = 0.d0
      depo_turb_glob = 0.d0
      depo_grav_glob = 0.d0
#endif

      SHDT_sv = 0.
      EVHDT_sv = 0.
      SRHDT_sv = 0.
      TRHDT_sv = 0.
      EVAP_sv = 0.
      TRHDT_sv2 = 0.

      call alloc_pbl_args(pbl_args)

C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NIsurf TIMES EVERY HOUR
C****
      DO NS=1,NIsurf
         MODDSF=MOD(NSTEPS+NS-1,NDASF*NIsurf+1)
         IF(MODDSF.EQ.0) IDACC(ia_srf)=IDACC(ia_srf)+1
         MODDD=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
         pbl_args % moddd = moddd
         pbl_args % ih = ih
         pbl_args % ihm = ihm
#endif

C**** ZERO OUT FLUXES ACCUMULATED OVER SURFACE TYPES
#ifdef TRACERS_ON
         if(NS==1) then
           trcsurf = 0.
           trcSurfByVol = 0.
         end if

         do ipatch=1,size(asflx)
           asflx(ipatch)%trsrfflx(:,:,:) = 0.
         enddo
#endif
      ! Part of PBL_DRV.f, uses atm/srf exchange variables
      ! Initializes boundary layer calc each surface time step
      call loadbl

#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          nx=nx+1
          ntix(nx) = n
        end if
      end do 
      ntx = nx
      pbl_args%ntix(1:ntm) = ntix(1:ntm)
      pbl_args%ntx = ntx

#endif

      call recalc_agrid_uv

C**** Copy first-layer atm. conditions into the 2D arrays in
C**** atm-surf. coupling data structures.
      call atm_exports_phasesrf

      call get_dbl

      call surface_diag0(moddd,ih,ihm)

C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

      EVAPLIM = 0. ; HTLIM=0.  ! need initialisation
#ifdef TRACERS_WATER
      tevaplim = 0.
#endif

C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      ! RSI = RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
      POICE=RSI(I,J)*PWATER
      POCEAN=PWATER-POICE

      ! ./model/ATMDYN_COM.f:!@var PEDN edge pressure (top of box) (mb)
      PS=PEDN(1,I,J)
      !@var  PEK  PEDN**KAPA
      PSK=PEK(1,I,J)
      !@var  PK   PMID**KAPA
      P1K=PK(1,I,J)
      !@var Q specific humidity (kg water vapor/kg air)
      Q1=Q(I,J,1)

      ! T=Temperature
      THV1=T(I,J,1)*(1.+Q1*xdelt)

      MA1=MA(1,I,J) !@var MA1 mass of lowest atmospheric layer (kg/m^2)

C****
      DO ITYPE=ITYPE_MIN,ITYPE_OCEANICE ! no earth or landice type
  !    ipbl(itype,i,j)=0
C****
C**** OPEN OCEAN/LAKE
C****
      if ( ITYPE == ITYPE_OCEAN ) then

      PTYPE=POCEAN
      IF (PTYPE.gt.0) THEN
      TG1=atmocn%GTEMP(I,J)
      TG2=atmocn%GTEMP2(I,J)   ! diagnostic only
      TR4=atmocn%GTEMPR(I,J)**4
      IF (FLAKE(I,J).gt.0) THEN
C**** limit evap/cooling if between MINMLD and 40cm, no evap below 40cm
        IF (MWL(I,J).lt.MINMLD*RHOW*FLAKE(I,J)*AXYP(I,J)) THEN
          EVAPLIM=MAX(0.5*(MWL(I,J)/(FLAKE(I,J)*AXYP(I,J))-0.4d0*RHOW),
     *         0d0)
        ELSE
          EVAPLIM=MWL(I,J)/(FLAKE(I,J)*AXYP(I,J))
     &         -(0.5*MINMLD+0.2d0)*RHOW
        END IF
#ifdef TRACERS_WATER
C**** limit on tracer evporation from lake
        TEVAPLIM(1:NTX)=EVAPLIM*atmocn%GTRACER(NTIX(1:NTX),I,J)
#endif
        HTLIM = GML(I,J)/(FLAKE(I,J)*AXYP(I,J)) + 0.5*LHM*EVAPLIM
        uocean = 0. ; vocean = 0. ! no velocities for lakes
      ELSE
        if (UOdrag.eq.1) then ! use uocean for drag calculation
C**** UOSURF,VOSURF are on atm A grid
            uocean = atmocn%uosurf(i,j)
            vocean = atmocn%vosurf(i,j)
        else
          uocean=0. ; vocean=0.
        end if
      END IF
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      atmocn%SOLAR(I,J)=atmocn%SOLAR(I,J)+DTSURF*SRHEAT
      ELHX=LHE

C**** pass salinity (zero for lakes)
      pbl_args%sss_loc=sss(i,j)
c**** sanity check (to prevent rare anomalies that will be dealt with by
C**** addice next time)
      TG1=max(TG1,tfrez(sss(i,j)))

      END IF
C****
C**** OCEAN/LAKE ICE
C****
      else if ( ITYPE == ITYPE_OCEANICE ) then

      PTYPE=POICE
      IF (PTYPE.gt.0) THEN
      IF (FLAKE(I,J).gt.0) THEN
        uocean = 0. ; vocean = 0. ! no dynamic ice for lakes
      ELSE
        if (UOdrag.eq.1) then ! use ice velocities in drag calculation
            uocean = atmice%uisurf(i,j)
            vocean = atmice%visurf(i,j)
        else
          uocean = 0. ; vocean =0.
        end if
      END IF
      TG1=TGRND(2,I,J)
      TG2=TGRN2(2,I,J)
      TR4=TGR4(2,I,J)
      SNOW=SNOWI(I,J)
      MSI1=SNOW+ACE1I ! snow and first layer ice mass (kg/m^2)
      MSI2=MSI(I,J)   ! second (physical) layer ice mass (kg/m^2)
C**** determine heat capacity etc for top ice layers
      dF1dTG = 2./(ACE1I/(RHOI*alami(TG1,1d3*((SSI(1,I,J)+SSI(2,I,J))
     *     /ACE1I)))+SNOW*BYRLS)
      HCG1 = dEidTi(TG1,1d3*(SSI(1,I,J)/(XSI(1)*MSI1)))*XSI(1)*MSI1
      HCG2 = dEidTi(TG2,1d3*(SSI(2,I,J)/(XSI(2)*MSI1)))*XSI(2)*MSI1

      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      atmice%SOLAR(I,J)=atmice%SOLAR(I,J)+DTSURF*SRHEAT
C**** fraction of solar radiation leaving layer 1 and 2
      IF (SRHEAT.gt.0) THEN ! only bother if there is sun
        call solar_ice_frac(SNOW,MSI2,si_atm%FLAG_DSWS(I,J),FSRI,2)
      ELSE
        FSRI(1:2) = 0
      END IF
      ELHX=LHS

C**** pass salinity in underlying water (zero for lakes)
      pbl_args%sss_loc=sss(i,j)
C**** Underlying ocean temperature with sanity check 
C**** (to prevent rare anomalies that will be dealt with by
C**** addice next time)
      TGO=max(atmocn%GTEMP(I,J),tfrez(sss(i,j)))
      
      END IF

      endif ! itype check
C****
      IF (PTYPE.gt.0) THEN
C****
C**** BOUNDARY LAYER INTERACTION
C****
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.

      TG=TG1+TF
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG_SAT=0.98d0*QG_SAT
      pbl_args%TG=TG   ! actual ground temperature
      pbl_args%TR4=TR4 ! radiative temperature K^4
      !pbl_args%ELHX=ELHX   ! relevant latent heat
      !pbl_args%QSOL=SRHEAT   ! solar heating
      pbl_args%TGV=TG*(1.+QG_SAT*xdelt)

C =====================================================================
      pbl_args%dtsurf = dtsurf
      !pbl_args%TKV=THV1*PSK     ! TKV is referenced to the surface pressure
      !pbl_args%ZS1=.5d-2*RGAS*pbl_args%TKV*MA1/PMID(1,I,J)
      pbl_args%qg_sat = qg_sat
      !pbl_args%qg_aver = qg_sat   ! QG_AVER=QG_SAT
      !pbl_args%hemi = sign(1d0,lat2d(i,j))
c      pbl_args%pole = pole
c      pbl_args%evap_max = 1.
c      pbl_args%fr_sat = 1. ! entire surface is saturated
      pbl_args%uocean = uocean
      pbl_args%vocean = vocean
c      pbl_args%psurf = PS
c      pbl_args%trhr0 = TRHR(0,I,J)
      pbl_args%ocean = (ITYPE.eq.1 .and. FOCEAN(I,J).gt.0)
      pbl_args%snow = SNOW

!     Calculate drag coefficients, wind speed, air density, etc.
!     PBL = "Planetary Boundary Layer"
C**** Call pbl to calculate near surface profile
      if(itype == 1) then
        CALL PBL(I,J,1,ITYPE,PTYPE,pbl_args,atmocn)
      else
        CALL PBL(I,J,1,ITYPE,PTYPE,pbl_args,atmice)
      endif
c#ifdef TRACERS_ON
c      trs(1:ntm) = pbl_args%trs(1:ntm)
c#endif

      us = pbl_args%us
      vs = pbl_args%vs
      ws = pbl_args%ws
      ws0 = pbl_args%ws0
      qsrf = pbl_args%qsrf
      CM = pbl_args%cm
      CH = pbl_args%ch
      CQ = pbl_args%cq
      TS=pbl_args%TSV/(1.+QSRF*xdelt)
C =====================================================================

C**** Adjust ground variables to account for skin effects
      TG = TG + pbl_args%dskin
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (pbl_args%ocean) QG_SAT=0.98d0*QG_SAT
      TG1 = TG - TF
      TR4=(sqrt(sqrt(TR4))+pbl_args%dskin)**4

C**** CALCULATE RHOSRF*CM*WS AND RHOSRF*CH*WS
      RHOSRF=100.*PS/(RGAS*pbl_args%TSV)
      RCDMWS=CM*WS*RHOSRF
      RCDHWS=CH*WS*RHOSRF
      RCDQWS=CQ*WS*RHOSRF
      RCDHDWS=CH*(WS-WS0)*RHOSRF
      RCDQDWS=CQ*(WS-WS0)*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2) (positive down)
      ! Including gustiness in the sensible heat flux:
      SHEAT=SHA*(RCDHWS*(TS-TG)+RCDHDWS*pbl_args%tprime)
      ! Including gustiness in the latent heat flux:
      EVHEAT=(LHE+TG1*SHV)*(RCDQWS*(QSRF-QG_SAT)+
     *                      RCDQDWS*pbl_args%qprime)
      TRHEAT=TRHR(0,I,J)-STBO*TR4

C**** CASE (1) ! FLUXES USING EXPLICIT TIME STEP FOR OCEAN POINTS
      if ( ITYPE == 1) then
        SHDT = DTSURF*SHEAT
        EVHDT=DTSURF*EVHEAT              ! latent heat flux
        TRHDT=DTSURF*TRHEAT

C**** CASE (2) ! FLUXES USING IMPLICIT TIME STEP FOR ICE POINTS
      else if ( ITYPE == 2 ) then

! heat flux on first/second/third layers (W/m^2)
        F1 = (TG1-TG2)*dF1dTG + SRHEAT*FSRI(1)
        F2 = SRHEAT*FSRI(2)
        ! Including gustiness in the latent heat flux:
        EVHEAT=LHE*(RCDQWS*(QSRF-QG_SAT)+RCDQDWS*pbl_args%qprime) ! why is this different to above?
        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        dSNdTG=-RCDHWS*SHA
        dQGdTG=QG_SAT*DQSATDT(TG,ELHX) ! d(QG)/dTG
        dEVdTG = -dQGdTG*LHE*RCDQWS ! d(EVHEAT)/dTG
        dTRdTG = -4*STBO*sqrt(sqrt(TR4))**3 ! d(TRHEAT)/dTG
        dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG

        T2DEN = HCG2+DTSURF*dF1dTG
        T2CON = DTSURF*(F1-F2)/T2DEN
        T2MUL = DTSURF*dF1dTG/T2DEN

        DFDTG=DF0DTG-dF1dTG
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)

        IF (TG1+dTG .GT. 0.) dTG = -TG1
        dT2 = T2CON+T2MUL*dTG

        SHDT  = DTSURF*(SHEAT +dTG*dSNdTG) ! sensible
        EVHDT = DTSURF*(EVHEAT+dTG*dEVdTG) ! latent
        TRHDT = DTSURF*(TRHEAT+dTG*dTRdTG) ! thermal flux (J/m^2)
        F1DT = DTSURF*(F1+(dTG*dF1dTG-dT2*dF1dTG))
        TG1 = TG1+dTG          ! first layer sea ice temperature (degC)
        TG2 = TG2+dT2          ! second layer sea ice temperature (degC)
        TGRN2(ITYPE,I,J) = TG2

      endif

C**** CALCULATE EVAPORATION
      DQ1X =EVHDT/((LHE+TG1*SHV)*MA1)
      EVHDT0=EVHDT
C**** Limit evaporation if lake mass is at minimum
      lim_lake_evap=.false.
      lim_dew=.false.
      IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and.
     *     (atmocn%EVAPOR(I,J)-DQ1X*MA1).gt.EVAPLIM) THEN
        if (QCHECK) WRITE(99,*) "Lake EVAP limited: I,J,EVAP,MWL",I,J
     *     ,atmocn%EVAPOR(I,J)-DQ1X*MA1,
     *     MWL(I,J)/(RHOW*FLAKE(I,J)*AXYP(I,J))
        DQ1X=(atmocn%EVAPOR(I,J)-EVAPLIM)*byMA(1,I,J)
        lim_lake_evap=.true.
      ELSEIF (DQ1X.GT.Q1) THEN
        DQ1X=Q1
        lim_dew=.true.
      ELSE
        GO TO 3720
      END IF
      EVHDT=DQ1X*(LHE+TG1*SHV)*MA1
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*MA1

#ifdef TRACERS_WATER
      DO NX=1,NTX
        N=NTIX(NX)
        if (tr_wd_TYPE(n).eq.nWATER) THEN
          call water_tracer_evap(
     &         itype, i,j, n,
     &         tg1, rcdqws, rcdqdws, evap, snow, qg_sat, qsrf,
     &         lim_lake_evap, flake(i,j), tevaplim(nx), ! arguments for lakes only
     &         lim_dew, dtsurf,
     &         asflx(itype)%TRM1(n,I,J), pbl_args%trs(nx),
     &         asflx(itype)%gtracer(n,i,j), trgrnd2(nx),
     &         pbl_args%trprime(nx),
     &         asflx(itype)%trsrfflx(n,i,j),asflx(itype)%trevapor(n,i,j)
     &     )
        END IF
      END DO
#endif

#ifdef TRACERS_ON

#ifdef TRACERS_GASEXCH_ocean

      call calc_gasexch(i,j,itype,ns,moddsf,
     &     ptype,pocean,rsi(i,j),rhosrf,tgo,
     &     dtsurf,pbl_args)

#endif /* TRACERS_GASEXCH_ocean */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)

      call collect_ocean_emissions(i,j,ptype,dtsurf,pbl_args)

#endif

#endif /* TRACERS_ON */

C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
      asflx(itype)%latht(i,j) = asflx(itype)%latht(i,j) + evhdt
      asflx(itype)%trheat(i,j) = asflx(itype)%trheat(i,j) + trhdt
C**** Limit heat fluxes out of lakes if near minimum depth
      IF (ITYPE.eq.1 .and. FLAKE(I,J).gt.0 .and.
     *     atmocn%E0(I,J)+F0DT+HTLIM.lt.0 .and.
     *     atmocn%E0(I,J)+F0DT.lt.0) THEN
        if (QCHECK.and.HTLIM.le.0) write(6,*) "NEW case:"
        if (QCHECK) write(6,*) "Limiting heat flux from lake",i,j,SHDT
     *       ,F0DT,atmocn%E0(I,J),DTSURF*SRHEAT,TRHDT,EVHDT,HTLIM
        SHDT =-(max(0d0,HTLIM)+atmocn%E0(I,J)+DTSURF*SRHEAT+TRHDT+EVHDT)
        F0DT = -atmocn%E0(I,J)-max(0d0,HTLIM)
        if (QCHECK) write(6,*) "New SHDT,F0DT",i,j,SHDT,F0DT
      END IF
      !E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      asflx(itype)%E0(I,J)=asflx(itype)%E0(I,J)+F0DT
      E1(ITYPE,I,J)=E1(ITYPE,I,J)+F1DT
      !EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP
      asflx(itype)%EVAPOR(I,J)=asflx(itype)%EVAPOR(I,J)+EVAP
#ifdef SCM
      if (J.eq.J_TARG.and.I.eq.I_TARG) then
          if (SCM_SURFACE_FLAG.eq.0.or.SCM_SURFACE_FLAG.eq.2) then
              EVPFLX = EVPFLX -(DQ1X*MA1)*(PTYPE/DTSURF)*LHE
              SHFLX = SHFLX - SHDT*PTYPE/DTSURF
c             write(iu_scm_prt,*) 'srf  evpflx shflx ptype ',
c    *                   EVPFLX,SHFLX,ptype
          endif
      endif
#endif
      TGRND(ITYPE,I,J)=TG1  ! includes skin effects
      TGR4(ITYPE,I,J) =TR4
C**** calculate correction for different TG in radiation and surface
      !dLWDT = DTSURF*(TRSURF(ITYPE,I,J)-TRHR(0,I,J))+TRHDT
      dLWDT = DTSURF*(asflx(itype)%TRUP_in_rad(I,J)-TRHR(0,I,J))+TRHDT
C**** final fluxes
#ifdef SCM
cccccc for SCM use ARM provided fluxes for designated box
      if ((I.eq.I_TARG.and.J.eq.J_TARG).and.SCM_SURFACE_FLAG.eq.1) then
           DTH1(I,J)=DTH1(I,J)
     &              +ash*DTSURF*ptype/(SHA*MA1)
           DQ1(I,J)=DQ1(I,J) + ALH*DTSURF*ptype/(MA1*LHE)
           SHFLX = SHFLX + ASH*ptype
           EVPFLX = EVPFLX + ALH*ptype
           write(iu_scm_prt,980) I,PTYPE,DTH1(I,J),DQ1(I,J),
     &           EVPFLX,SHFLX
 980       format(1x,'SURFACE ARM   I PTYPE DTH1 DQ1 evpflx shflx',
     &            i5,f9.4,f9.5,f9.6,f9.5,f9.5)
      else
#endif
      asflx(itype)%DTH1(I,J)=-(SHDT+dLWDT)/(SHA*MA1) ! +ve up
      asflx(itype)%sensht(i,j) = asflx(itype)%sensht(i,j)+SHDT
      asflx(itype)%DQ1(I,J) = -DQ1X
#ifdef SCM
      if (i.eq.I_TARG.and.j.eq.J_TARG) then
          write(iu_scm_prt,988) I,PTYPE,DTH1(I,J),DQ1(I,J),SHDT,dLWDT
 988      format(1x,'988 SURFACE GCM  I PTYPE DTH1 DQ1 SHDT dLWDT ',
     &           i5,f9.4,f9.5,f9.6,f12.4,f10.4)
      endif
      endif
#endif
      DMUA_IJ=RCDMWS*(US-UOCEAN)
      DMVA_IJ=RCDMWS*(VS-VOCEAN)
      asflx(itype)%DMUA(I,J) = asflx(itype)%DMUA(I,J) + DMUA_IJ*DTSURF
      asflx(itype)%DMVA(I,J) = asflx(itype)%DMVA(I,J) + DMVA_IJ*DTSURF
      asflx(itype)%uflux1(i,j) = DMUA_IJ
      asflx(itype)%vflux1(i,j) = DMVA_IJ

C****
C**** SAVE SOME TYPE DEPENDENT FLUXES/DIAGNOSTICS
C****
!!!      CASE (1)  ! ocean
      if ( ITYPE == 1 ) then
        ! todo: use the (gtemps-gtemp) form outside this loop
        ! after figuring why it doesn't always give the same result
        IF(MODDSF.EQ.0)
     &       AIJ(I,J,IJ_DSKIN)=AIJ(I,J,IJ_DSKIN)+pbl_args%dskin
      endif
C****
      END IF
      END DO   ! end of itype loop
      END DO   ! end of I loop
      END DO   ! end of J loop

      atmice%gtemp2(:,:) = tgrn2(2,:,:)

C****
C**** dynamic vegetation time step
C****
!!! probably don't need this call unless something can be done
!   separately from ground hydrology on i,j grid
!      call step_dveg(dtsurf)


C****
C**** LAND ICE
C****
      call downscale_temperature_li
#ifdef GLINT2
      do ihp=1,ubound(atmglas_hp,1)
        CALL SURFACE_LANDICE(NS==1,MODDD,DTSURF,atmglas_hp(ihp),ihp)
      enddo

      ! Convert variables set in surface_landice from height point
      ! to height class space.
      call bundle_hp_to_hc(bundle_surface_landice,
     &    atmglas_hp(1:), atmglas(1:))
#else
      do ipatch=1,ubound(atmglas,1)
!        print *,'SURFACE_LANDICE',ipatch
        CALL SURFACE_LANDICE(NS==1,MODDD,DTSURF,atmglas(ipatch),ipatch)
      enddo

#endif


C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODDD)
      IF(NS==Nisurf) CALL GROUND_E ! diagnostic only - should be merged with EARTH

      atmlnd%gtemps(:,:) = atmlnd%gtemp(:,:) ! fill in for diags


#ifdef TRACERS_DRYDEP

C****
C**** Calculate Tracer Dry Deposition (including gravitational settling)
C****
      do j=j_0,j_1
      do i=i_0,imaxj(j)

      do ipatch=1,size(asflx)

        ptype = asflx(ipatch)%ftype(i,j)
        if(ptype.le.0.) cycle

      do n=1,ntm
        if(dodrydep(n))then
          depvel = asflx(ipatch)%dep_vel(n,i,j)
          gsvel = asflx(ipatch)%gs_vel(n,i,j)
          tdryd=asflx(ipatch)%drydflx(n,i,j)         ! kg/m2
          rtsdt = -tdryd/(depvel+gsvel+teeny)
          rts = rtsdt/dtsurf                             ! kg*s/m^3
          tdd = tdryd                                    ! kg/m2
          td1 = (asflx(ipatch)%trsrfflx(n,i,j)
     &          +asflx(ipatch)%trflux_prescr(n,i,j)
     &         )*dtsurf         ! kg/m2
          if (trm(i,j,1,n)+td1+tdd.le.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd surface",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,pbl_args%trs(nx),pbl_args%trtop(nx)
            tdd= -max(trm(i,j,1,n)+td1,0d0)*byaxyp(i,j)
            tdryd=tdd
          end if
          asflx(ipatch)%trsrfflx(n,i,j)=
     &         asflx(ipatch)%trsrfflx(n,i,j)+tdd/dtsurf

! trdrydep downward flux by surface type (kg/m^2)
          asflx(ipatch)%trdrydep(n,i,j)=
     &         asflx(ipatch)%trdrydep(n,i,j) - tdryd
! diagnose turbulent and settling fluxes separately
          taijn(i,j,tij_drydep,n)=taijn(i,j,tij_drydep,n) +
     &         ptype*rtsdt*depvel
          taijn(i,j,tij_gsdep ,n)=taijn(i,j,tij_gsdep ,n) +
     &         ptype*rtsdt* gsvel
#ifdef TRACERS_COSMO
          if (n .eq. n_Be7) BE7D_acc(i,j)=BE7D_acc(i,j)+ptype*rtsdt
     *         *depvel+ptype*rtsdt* gsvel
#endif

          if (itcon_dd(n,1).gt.0) call inc_diagtcb(i,j,-
     &     ptype*rtsdt*axyp(i,j)*depvel,itcon_dd(n,1),n)
          if (itcon_dd(n,2).gt.0) call inc_diagtcb(i,j,-
     &     ptype*rtsdt*axyp(i,j)*gsvel,itcon_dd(n,2),n)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
c**** for subdaily diagnostics
          depo_turb_glob( i, j, n ) = depo_turb_glob( i, j, n )
     &         + ptype * rts * depvel / nisurf
          depo_grav_glob( i, j, n ) = depo_grav_glob( i, j, n )
     &         + ptype * rts * gsvel / nisurf
#endif

        end if

      enddo ! tracer
      enddo ! ipatch
      enddo ! i
      enddo ! j
#endif

      atmice%gtemp(:,:) = tgrnd(2,:,:)
      atmice%gtemps(:,:) = tgrnd(2,:,:)

      atmocn%gtemps(:,:) = tgrnd(1,:,:)

C****
C**** Create gridbox composite values
C****
      call avg_patches_pbl_exports(grid,asflx,atmsrf)
      call avg_patches_srfflx_exports(grid,asflx,atmsrf)
      call avg_patches_srfstate_exports(grid,asflx,atmsrf)

C****
C**** UPDATE FIRST LAYER QUANTITIES
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        FTEVAP=0
        IF (atmsrf%DTH1(I,J)*T(I,J,1).lt.0)
     &       FTEVAP=-atmsrf%DTH1(I,J)/(T(I,J,1)*PK(1,I,J))
        FQEVAP=0
        IF (atmsrf%DQ1(I,J).lt.0.and.Q(I,J,1).gt.0)
     &       FQEVAP=-atmsrf%DQ1(I,J)/Q(I,J,1)
! Z-moments should be set from PBL
        TMOM(:,I,J,1) = TMOM(:,I,J,1)*(1.-FTEVAP)
        QMOM(:,I,J,1) = QMOM(:,I,J,1)*(1.-FQEVAP)
#ifdef WATER_PROPORTIONAL
        if(fqevap.gt.0.) then
          do n=1,ntm
            trmom(:,i,j,1,n) = trmom(:,i,j,1,n)*(1.-fqevap)
          enddo
        endif
#endif
        IF ( Q(I,J,1)+atmsrf%DQ1(I,J) .LT. qmin ) THEN
          QMOM(:,I,J,1)=0.
#ifdef WATER_PROPORTIONAL
          do n=1,ntm
            trmom(:,i,j,1,n) = 0.
          enddo
#endif
        ENDIF
c****   retrieve fluxes
        uflux1(i,j)=atmsrf%uflux1(i,j)
        vflux1(i,j)=atmsrf%vflux1(i,j)
        tflux1(i,j) = -atmsrf%dth1(i,j)*MA(1,I,J) / dtsurf
        qflux1(i,j) = -atmsrf%dq1(i,j) *MA(1,I,J) / dtsurf
      END DO 
      END DO

c create land ice composite values for a few diagnostics that
c happen to be sampled within this loop over surface sub-timesteps,
c in surface_diag1()
      if(moddsf==0) then
        ! just to get tsavg.
        call avg_patches_pbl_exports(grid,
c     &     atmglas,atmgla, ! gfortran prob. if passed as class() args
     &       atmglas(:)%atmsrf_xchng_vars,atmgla%atmsrf_xchng_vars,
     &       rel=.true.)
      endif
      ! to get gtemp
      call avg_patches_srfstate_exports(grid,
c     &     atmglas,atmgla, ! gfortran prob. if passed as class() args
     &     atmglas(:)%atmsrf_xchng_vars,atmgla%atmsrf_xchng_vars,
     &     rel=.true.)

      call surface_diag1(dtsurf,moddsf,trhdt_sv2)

#ifdef TRACERS_ON
#ifdef TRACERS_TOMAS
C**** Apply subgrid coagulation for freshly emitted particles
      call subgridcoag_drv_2D(dtsurf)

#endif
C****
C**** Apply tracer surface sources and sinks
C****
      call apply_tracer_2Dsource(dtsurf)
#endif

c****
c**** apply surface fluxes to the first layer of the atmosphere
c****  (replaced with dummy sub when ATURB is used)
c****
      call apply_fluxes_to_atm(dtsurf)

C**** Call dry convection or aturb depending on rundeck
      CALL ATM_DIFFUS(1,1,dtsurf)

      call surface_diag2(moddd,ih,ihm,
     &     srhdt_sv,trhdt_sv,evhdt_sv,shdt_sv,evap_sv)

C****
      END DO   ! end of surface time step

      call dealloc_pbl_args(pbl_args)

#ifdef mjo_subdd
      call surface_diag_mjo
#endif

#ifdef CALCULATE_FLAMMABILITY
      call flammability_drv
#endif

      atmice%e1(:,:) = e1(2,:,:)

C**** APPLY SURFACE FLUXES TO LAND ICE (and modify fluxes as well)
#ifdef GLINT2
      do ipatch=1,ubound(atmglas_hp,1)
        CALL GROUND_LI(atmglas_hp(ipatch),ipatch)
      enddo

      ! Convert variables set in surface_landice from height point
      ! to height class space.
      call bundle_hp_to_hc(bundle_ground_li,
     &     atmglas_hp(1:), atmglas(1:))
#else
      do ipatch=1,ubound(atmglas,1)
        CALL GROUND_LI(atmglas(ipatch),ipatch)
      enddo
#endif

c create composite values for land ice for diagnostics
      call avg_patches_srfflx_exports(grid,
c     &     atmglas,atmgla, ! gfortran prob. if passed as class() args
     &     atmglas(:)%atmsrf_xchng_vars,atmgla%atmsrf_xchng_vars,
     &     rel=.true.)
      call avg_patches_srfflx_exports_gla(grid,atmglas,atmgla,
     &     rel=.true.)
      call avg_patches_srfstate_exports(grid,
c     &     atmglas,atmgla, ! gfortran prob. if passed as class() args
     &     atmglas(:)%atmsrf_xchng_vars,atmgla%atmsrf_xchng_vars,
     &     rel=.true.)

c create grid-composite values
      call avg_patches_srfstate_exports(grid,asflx,atmsrf)

      call surface_diag1a
      call surface_diag3

C****
C**** LAKES
C****
      CALL UNDERICE(si_atm,icelak,atmocn)
      CALL GROUND_SI(si_atm,icelak,atmice,atmocn)
C**** APPLY FLUXES TO LAKES AND DETERMINE ICE FORMATION
      CALL GROUND_LK
         CALL CHECKT ('GRNDLK')
C**** CALCULATE RIVER RUNOFF FROM LAKE MASS
      CALL RIVERF
      CALL FORM_SI(si_atm,icelak,atmice)
      CALL SI_diags(si_atm,icelak,atmice)

c calculate global integral of heat of river discharge
      if(atmocn%need_eflow_gl) then
        do j=J_0,J_1
        do i=I_0,I_1
          atmocn%work1(i,j) = 0.
          if (focean(i,j).eq.0.) cycle
          atmocn%work1(i,j) = atmocn%eflowo(i,j)*focean(i,j)*axyp(i,j)
        enddo
        enddo
        call globalsum(grid, atmocn%work1, atmocn%eflow_gl, all=.true.)
      endif

         CALL CHECKT ('SURFACE')
         IF (MODD5S.EQ.0) CALL DIAGCA (5)

      call stopTimer('SURFACE()')

      RETURN
C****
      END SUBROUTINE SURFACE

      subroutine surface_diag_post_precip_li
      ! the contents of this routine will be merged with other
      ! landice diag bits once the "precip" entry point to landice
      ! physics is merged with the others.
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE DIAG_COM, only : aij=>aij_loc
     &     ,itlandi,jreg,ij_runli,j_run,ij_f1li
      USE FLUXES, only : flice,atmgla
      USE GEOM, only : imaxj
      implicit none
      INTEGER :: JR,I,J,I_0,I_1,J_0,J_1
      REAL*8 :: PLICE
      CALL getDomainBounds(grid, I_STRT=I_0,I_STOP=I_1,
     &  J_STRT=J_0,J_STOP=J_1)
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        JR=JREG(I,J)
        PLICE = FLICE(I,J)
        IF(PLICE.LE.0.) CYCLE
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+atmgla%RUNO(I,J)*PLICE
        CALL INC_AJ(I,J,ITLANDI,J_RUN ,  atmgla%RUNO(I,J)*PLICE)
        CALL INC_AREG(I,J,JR,J_RUN ,  atmgla%RUNO(I,J)*PLICE)
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI)+atmgla%E1(I,J)*PLICE
      ENDDO
      ENDDO
      return
      end subroutine surface_diag_post_precip_li

#ifdef mjo_subdd
      subroutine surface_diag_mjo
      USE CONSTANT, only : undef
      USE MODEL_COM, only : lm
      USE ATM_COM, only : phi,MWs
      USE DIAG_COM, only :
     *     ,PW_acc, E_acc,sst_avg,p_avg,lwu_avg
     *     ,u_avg,v_avg,w_avg,t_avg,q_avg,r_avg,z_avg
      implicit none
C**** Accumulate subdaily precipitable water (kg/m^2) PW_acc ***
C****   longwave upward flux lwu_avg,surface pres p_avg, sst sst_avg
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        PW_acc(I,J) = PW_acc(I,J) + Sum(Q(I,J,:)*MA(:,I,J))
        p_avg(I,J)=p_avg(I,J)+P(I,J)
        if (FOCEAN(I,J).gt.0) then
          sst_avg(i,j)=sst_avg(i,j)+atmocn%GTEMP(i,j)
        else
          sst_avg(i,j)=undef
        end if
        PLAND=FLAND(I,J)
        PWATER=1.-PLAND
        POICE=RSI(I,J)*PWATER
        POCEAN=PWATER-POICE
        PEARTH=FEARTH(I,J)
        PLICE=FLICE(I,J)
        lwu_avg(i,j)=lwu_avg(i,j)+STBO*(
     &        POCEAN*atmocn%GTEMPR(I,J)**4
     &       +POICE *atmice%GTEMPR(I,J)**4
     &       +PLICE *atmgla%GTEMPR(I,J)**4
     &       +PEARTH*atmlnd%GTEMPR(I,J)**4)
C**** Accumulate 3D subdaily quantities
       DO K=1,LM
        u_avg(I,J,K)=u_avg(I,J,K)+u(I,J,K)
        v_avg(I,J,K)=v_avg(I,J,K)+v(I,J,K)
        IF (K < LM) THEN
          w_avg(I,J,K)=w_avg(I,J,K)+MWs(I,J,K)*byaxyp(I,J)
        END IF
        t_avg(I,J,K)=t_avg(I,J,K)+t(I,J,K)*pk(K,I,J)
        q_avg(I,J,K)=q_avg(I,J,K)+q(I,J,K)
        r_avg(I,J,K)=r_avg(I,J,K)+q(I,J,K)/
     &      qsat(t(I,J,K)*pk(K,I,J),lhe,pmid(K,I,J))
        z_avg(I,J,K)=z_avg(I,J,K)+phi(I,J,K)*bygrav
       END DO
      END DO
      END DO
      return
      end subroutine surface_diag_mjo
#endif

      subroutine surface_diag0(moddd,ih,ihm)
      USE ATM_COM, only : t,q,pek,pedn
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE DIAG_COM, ONLY : idd_spr,idd_pt5,idd_q5
     &     ,ndiuvar,ndiupt,ijdd,adiurn=>adiurn_loc
#ifndef NO_HDIURN
     *     ,hdiurn=>hdiurn_loc
#endif
      implicit none
      integer, intent(in) :: ih,ihm
      integer, intent(in) :: moddd

C**** some shorthand indices and arrays for diurn diags
      INTEGER, PARAMETER :: n_idx1 = 11
      INTEGER :: idx1(n_idx1)
      REAL*8 :: tmp(NDIUVAR)

      integer :: kr,ii,i,j, j_0, j_1, i_0,i_1

C**** Initialise constant indices
      idx1 = (/ IDD_SPR, (IDD_PT5+ii-1,ii=1,5), (IDD_Q5+ii-1,ii=1,5) /)

      CALL getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
      IF(MODDD.EQ.0) THEN
        DO KR=1,NDIUPT
C**** For distributed implementation - ensure point is on local process.
          I = IJDD(1,KR)
          J = IJDD(2,KR)
          IF ((J >= J_0) .AND. (J <= J_1) .AND.
     &        (I >= I_0) .AND. (I <= I_1)) THEN
            tmp(IDD_SPR)=PEDN(1,I,J)
            do ii=1,5
              tmp(IDD_PT5+ii-1)=PEK(1,I,J)*T(I,J,ii)
              tmp(IDD_Q5+ii-1) =Q(I,J,ii)
            end do 
            ADIURN(idx1(:),kr,ih)=ADIURN(idx1(:),kr,ih)+tmp(idx1(:))
#ifndef NO_HDIURN
            HDIURN(idx1(:),kr,ihm)=HDIURN(idx1(:),kr,ihm)+tmp(idx1(:))
#endif
          END IF
        END DO 
      END IF
      return
      end subroutine surface_diag0

      subroutine surface_diag1(dtsurf,moddsf,trhdt_sv2)
      USE CONSTANT, only : tf,bygrav,rhows
      USE MODEL_COM, only : dtsrc,itime
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE GEOM, only : axyp,imaxj,byaxyp
      USE SOMTQ_COM, only : mz
      USE ATM_COM, only : byMA
      USE RAD_COM, only : trhr
#ifdef TRACERS_ON
      use OldTracer_mod, only: itime_tr0, needtrs
      USE TRACER_COM, only : NTM,trm,trmom
#ifdef TRACERS_WATER
      use OldTracer_mod, only: nWATER,tr_wd_TYPE

#endif
#endif
      USE DIAG_COM, only : aij=>aij_loc,aijmm
     &     ,itocean,itoice,itlake,itlkice,itlandi,itearth
     *     ,tdiurn,jreg
     *     ,ij_tsli,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,ij_tg1,ij_tgo,ij_rhs
     *     ,j_tsrf,j_type,j_tg1,j_tg2
     *     ,ij_pblht,ij_dskin
     *     ,ij_gusti,ij_mccon,ij_sss,ij_trsup,ij_trsdn,ij_ssh
     *     ,ij_silwu, ij_silwd
     *     ,ij_sish,ij_evap
     *     ,ij_tsurfmin,ij_tsurfmax
#if (defined mjo_subdd) || (defined etc_subdd)
     *     ,qsen_avg,qlat_avg,pblht_acc
#endif
#ifdef mjo_subdd
     *     ,E_acc
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
#endif
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : si_atm
      USE FLUXES, only :
     *      uflux1,vflux1,tflux1,qflux1
     &     ,nisurf,fland,flice,focean
     &     ,atmocn,atmice,atmgla,atmlnd,asflx4,atmsrf
#ifdef TRACERS_ON
#ifndef SKIP_TRACER_DIAGS
      USE TRDIAG_COM, only : taijn=>taijn_loc,
     *      taijs=>taijs_loc,jls_isrc, tij_surf,
     *      tij_surfbv, tij_evap, tij_grnd
#endif /*SKIP_TRACER_DIAGS*/
#ifdef TRACERS_ON
      use trdiag_com, only: trcsurf,trcSurfByVol
#endif
#endif

!@var DDMS downdraft mass flux in kg/(m^2 s), (i,j)
      USE CLOUDS_COM, only : DDMS

      implicit none
      real*8, intent(in) :: dtsurf
      integer, intent(in) :: moddsf
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: trhdt_sv2

      integer :: itype,idtype,jr
      integer :: n,i,j, j_0, j_1, i_0,i_1
      real*8 :: ptype,byNIsurf
      real*8 :: tg1,tg2,ts

      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,SSS
      REAL*8, DIMENSION(:,:,:), POINTER :: SSI

      byNIsurf = 1d0/real(nisurf,kind=8)

      RSI => SI_ATM%RSI
      MSI => SI_ATM%MSI
      SNOWI => SI_ATM%SNOWI
      SSI => SI_ATM%SSI
      SSS => atmocn%SSS

      CALL getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=j_0,j_1
      do i=i_0,imaxj(j)
        if(moddsf.eq.0) then
          trhdt_sv2(i,j) = atmsrf%trheat(i,j)-trhdt_sv2(i,j)
        else
          trhdt_sv2(i,j) = atmsrf%trheat(i,j)
        endif
      enddo
      enddo

      do j=j_0,j_1
      do i=i_0,imaxj(j)
#if (defined mjo_subdd) || (defined etc_subdd)
C**** SUBDD qsen_avg,qlat_avg for sensible/latent heat flux ***
        qlat_avg(i,j)=qlat_avg(i,j)+qflux1(i,j)
        qsen_avg(i,j)=qsen_avg(i,j)+tflux1(i,j)
#endif
C**** Diurnal cycle of temperature diagnostics
        tdiurn(i,j,5)=tdiurn(i,j,5)+(atmsrf%tsavg(i,j)-tf)
        if(atmsrf%tsavg(i,j).gt.tdiurn(i,j,6))
     &       tdiurn(i,j,6)=atmsrf%tsavg(i,j)
        if(atmsrf%tsavg(i,j).lt.tdiurn(i,j,9))
     &       tdiurn(i,j,9)=atmsrf%tsavg(i,j)
C*** min/max tsurf
        aijmm(i,j,ij_tsurfmin) =
     &       max( -(atmsrf%tsavg(i,j)-tf), aijmm(i,j,ij_tsurfmin) )
        aijmm(i,j,ij_tsurfmax) =
     &       max(  (atmsrf%tsavg(i,j)-tf), aijmm(i,j,ij_tsurfmax) )
        aij(i,j,ij_evap)=aij(i,j,ij_evap)-dtsurf*qflux1(i,j)
#ifdef mjo_subdd
C**** SUBDD E_acc for evaporation *** 
        E_acc(I,J)=E_acc(I,J)-dtsurf*qflux1(i,j)
#endif
      enddo
      enddo

      if(moddsf.eq.0) then
        do j=j_0,j_1
        do i=i_0,imaxj(j)
          aij(i,j,ij_tg1)=aij(i,j,ij_tg1)+atmsrf%gtemps(i,j)
          do itype=1,4
            ptype = asflx4(itype)%ftype(i,j)
            if(ptype.le.0.) cycle
            tg1 = asflx4(itype)%gtemps(i,j) !tgrnd(itype,i,j)
C****
C**** ACCUMULATE DIAGNOSTICS FOR EACH SURFACE TIME STEP AND ITYPE AND REGION
C****
            JR=JREG(I,J)
            if(itype==1) then
              if(focean(i,j)>0.) then
                idtype = itocean
              else
                idtype = itlake
              endif
            elseif(itype==2) then
              if(focean(i,j)>0.) then
                idtype = itoice
              else
                idtype = itlkice
              endif
            elseif(itype==3) then
              idtype = itlandi
            else
              idtype = itearth
            endif

            ts = asflx4(itype)%tsavg(i,j)
            CALL INC_AJ(I,J,IDTYPE,J_TSRF,(TS-TF)*PTYPE)
            CALL INC_AREG(I,J,JR,J_TSRF,(TS-TF)*PTYPE)
            CALL INC_AJ(I,J,IDTYPE,J_TYPE,        PTYPE)
            CALL INC_AJ(I,J,IDTYPE,J_TG1 ,    TG1*PTYPE)
            CALL INC_AREG(I,J,JR,J_TG1 ,    TG1*PTYPE)
            tg2 = asflx4(itype)%gtemp2(i,j)
            CALL INC_AJ(I,J,IDTYPE,J_TG2 ,    TG2*PTYPE)
            CALL INC_AREG(I,J,JR,J_TG2 ,    TG2*PTYPE)
          enddo
          aij(i,j,ij_us)=aij(i,j,ij_us)+atmsrf%usavg(i,j)
          aij(i,j,ij_vs)=aij(i,j,ij_vs)+atmsrf%vsavg(i,j)
          aij(i,j,ij_ws)=aij(i,j,ij_ws)+atmsrf%wsavg(i,j)
          aij(i,j,ij_ts)=aij(i,j,ij_ts)+(atmsrf%tsavg(i,j)-tf)
          aij(i,j,ij_qs)=aij(i,j,ij_qs)+atmsrf%qsavg(i,j)
          aij(i,j,ij_rhs)=aij(i,j,ij_rhs)+atmsrf%rsavg(i,j)

          aij(i,j,ij_taus)=aij(i,j,ij_taus)+atmsrf%tauavg(i,j)
          aij(i,j,ij_tauus)=aij(i,j,ij_tauus)+atmsrf%uflux1(i,j)
          aij(i,j,ij_tauvs)=aij(i,j,ij_tauvs)+atmsrf%vflux1(i,j)

          if(DDMS(I,J).lt.0.) ! ddms < 0 for down draft
     &         AIJ(I,J,ij_mccon)=AIJ(I,J,ij_mccon)+1.

          aij(i,j,ij_trsdn)=aij(i,j,ij_trsdn)+trhr(0,i,j)


          aij(i,j,ij_trsup)=aij(i,j,ij_trsup)+
     &         (trhr(0,i,j)-trhdt_sv2(i,j)/dtsurf)

          aij(i,j,ij_gusti)=aij(i,j,ij_gusti)+atmsrf%gustiwind(i,j)

          aij(i,j,ij_pblht)=aij(i,j,ij_pblht)+atmsrf%dblavg(i,j)

#if (defined mjo_subdd) || (defined etc_subdd)
C**** SUBDD qblht_acc for PBL height *** YH Chen ***
          pblht_acc(I,J)=pblht_acc(I,J)+dblavg(i,j)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          aij(i,j,ij_wsgcm)=aij(i,j,ij_wsgcm)+atmsrf%wsgcm(i,j)
          aij(i,j,ij_wspdf)=aij(i,j,ij_wspdf)+atmsrf%wspdf(i,j)
          aij(i,j,ij_wdry)=aij(i,j,ij_wdry)+atmsrf%wsubwd(i,j)
          aij(i,j,ij_wtke)=aij(i,j,ij_wtke)+atmsrf%wsubtke(i,j)
          aij(i,j,ij_wmoist)=aij(i,j,ij_wmoist)+atmsrf%wsubwm(i,j)
#endif

        enddo
        enddo
      endif

C****
C****
C**** SAVE SOME TYPE DEPENDENT FLUXES/DIAGNOSTICS
C****
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        JR=JREG(I,J)
        ptype = focean(i,j)
        if(ptype.gt.0.) then
          IF(MODDSF.eq.0) THEN
           AIJ(I,J,IJ_TGO)=AIJ(I,J,IJ_TGO)+atmocn%GTEMP(I,J)*FOCEAN(I,J)
            AIJ(I,J,IJ_SSS)=AIJ(I,J,IJ_SSS)+SSS(I,J)*FOCEAN(I,J)
            AIJ(I,J,IJ_SSH)=AIJ(I,J,IJ_SSH)+(atmocn%OGEOZA(I,J)*BYGRAV+
     *           RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS)*FOCEAN(I,J)
          END IF
        endif
        ! open water
        ptype = asflx4(1)%ftype(i,j)
        !if(ptype.gt.0.) then
        !endif
        ! floating ice
        ptype = asflx4(2)%ftype(i,j)
        if(ptype.gt.0.) then
          tg1 = atmice%gtemps(i,j) !tgrnd(2,i,j)
          IF (TG1.GT.TDIURN(I,J,7)) TDIURN(I,J,7) = TG1
        endif
        ! land ice
        ptype = asflx4(3)%ftype(i,j)
        if(ptype.gt.0.) then
          tg1 = atmgla%gtemps(i,j) !tgrnd(3,i,j)
          IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
          IF (MODDSF.eq.0) then
          AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(atmgla%TSAVG(I,J)-TF)*PTYPE
          ENDIF
        endif
      enddo
      enddo

#ifdef TRACERS_ON
#ifndef SKIP_TRACER_DIAGS
      do j=j_0,j_1
      do i=i_0,imaxj(j)
C**** Save surface tracer concentration whether calculated or not
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if(.not. needtrs(n)) then
            atmsrf%travg(n,i,j) = byMA(1,i,j)*byaxyp(i,j)*
     &           max(trm(i,j,1,n)-trmom(mz,i,j,1,n),0d0)
            atmsrf%travg_byvol(n,i,j) =
     &           atmsrf%travg(n,i,j)*atmsrf%rhoavg(i,j)
          endif
          taijn(i,j,tij_surf  ,n) = taijn(i,j,tij_surf  ,n)
     &         +atmsrf%travg(n,i,j)
          taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)
     &         +atmsrf%travg_byvol(n,i,j)
          trcsurf(i,j,n)=trcsurf(i,j,n)
     &         +atmsrf%travg(n,i,j)*byNIsurf
          trcSurfByVol(i,j,n)=trcSurfByVol(i,j,n)
     &         +atmsrf%travg_byvol(n,i,j)*byNIsurf
#ifdef TRACERS_WATER
          if (tr_wd_type(n).eq.nWater) then

            taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     &           atmsrf%trevapor(n,i,j)
            if (jls_isrc(1,n)>0) call inc_tajls2(i,j,1,jls_isrc(1,n),
     &           atmsrf%trevapor(n,i,j))

            do itype=1,3
              ptype = asflx4(itype)%ftype(i,j)
! ==== DEBUGGING: Remove these lines to get exact match in regression tests
                if (focean(i,j)>0 .and. jls_isrc(2,n)>0) call inc_tajls2
     &        (i,j,1,jls_isrc(2,n),asflx4(itype)%trevapor(n,i,j)*ptype)
! ===== END DEBUGGING
            enddo
c            if (focean(i,j)>0 .and. jls_isrc(2,n)>0) call inc_tajls2
c     &        (i,j,1,jls_isrc(2,n),asflx4(1)%trevapor(n,i,j)*focean(i,j))
          end if
          taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     &         atmsrf%gtracer(n,i,j)
#endif
        end if
      end do ! tracer n
      enddo
      enddo
#endif /*SKIP_TRACER_DIAGS*/
#endif

      return
      end subroutine surface_diag1

      subroutine surface_diag1a
      USE MODEL_COM, only : dtsrc
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE GEOM, only : axyp,imaxj,byaxyp
      USE RAD_COM, only : trhr,trsurf
#ifdef TRACERS_ON
      USE TRACER_COM, only : NTM
#endif
      USE DIAG_COM, only : oa,aij=>aij_loc,aijmm
     &     ,itocean,itoice,itlake,itlkice,itlandi,itearth
     *     ,tdiurn,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_popocn
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,ij_tg1,ij_evap,ij_evapo,ij_tgo,ij_f0oc,ij_rhs
     *     ,ij_evapi,ij_f0li,ij_evapli,j_evap,j_evhdt,j_lwcorr
     *     ,j_tsrf,j_shdt,j_trhdt,j_type,j_tg1,j_tg2
     *     ,ij_pblht,ij_dskin
     *     ,ij_gusti,ij_mccon,ij_sss,ij_trsup,ij_trsdn,ij_fwoc,ij_ssh
     *     ,ij_silwu, ij_silwd
     *     ,ij_sish ,ij_popwat
     *     ,ij_tsurfmin,ij_tsurfmax
     *     ,j_rvrd,j_ervr,ij_micb,ij_eicb,ij_rsit,ij_li
      USE FLUXES, only :
     &      focean
     &     ,atmocn,atmice,atmgla,atmlnd,asflx4,atmsrf
#ifdef TRACERS_ON
#ifndef SKIP_TRACER_DIAGS
      USE TRDIAG_COM, only : taijn=>taijn_loc, tij_evap
#endif /*SKIP_TRACER_DIAGS*/
#endif

#ifdef TRACERS_WATER
#ifdef TRACERS_OCEAN
      USE TRDIAG_COM, only : tij_icb
#endif
#endif

      implicit none
      integer :: itype,idtype,jr
      integer :: n,i,j, j_0, j_1, i_0,i_1
      real*8 :: ptype
      real*8 :: evap,shdt,evhdt,trhdt,dlwdt


      CALL getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=j_0,j_1
      do i=i_0,imaxj(j)
        aij(i,j,ij_neth)=aij(i,j,ij_neth)+atmsrf%e0(i,j)
        aij(i,j,ij_shdt)=aij(i,j,ij_shdt)+atmsrf%sensht(i,j)
        aij(i,j,ij_srtr)=aij(i,j,ij_srtr)+
     &         (atmsrf%solar(i,j)+atmsrf%trheat(i,j))

        do itype=1,4
          ptype = asflx4(itype)%ftype(i,j)
          if(ptype.le.0.) cycle


C****
C**** ACCUMULATE DIAGNOSTICS FOR EACH SURFACE TIME STEP AND ITYPE AND REGION
C****
          JR=JREG(I,J)
          if(itype==1) then
            if(focean(i,j)>0.) then
              idtype = itocean
            else
              idtype = itlake
            endif
          elseif(itype==2) then
            if(focean(i,j)>0.) then
              idtype = itoice
            else
              idtype = itlkice
            endif
          elseif(itype==3) then
            idtype = itlandi
          else
            idtype = itearth
          endif
          evap = asflx4(itype)%EVAPOR(I,J)
          shdt = asflx4(itype)%SENSHT(I,J)
          evhdt = asflx4(itype)%LATHT(I,J)
          trhdt = asflx4(itype)%TRHEAT(I,J)
          dLWDT = DTsrc*(TRSURF(ITYPE,I,J)-TRHR(0,I,J))+TRHDT
          CALL INC_AJ(I,J,IDTYPE,J_EVAP , EVAP*PTYPE)
          CALL INC_AREG(I,J,JR,J_EVAP ,EVAP *PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_EVHDT,EVHDT*PTYPE)
          CALL INC_AREG(I,J,JR,J_EVHDT,EVHDT*PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_SHDT , SHDT*PTYPE)
          CALL INC_AREG(I,J,JR,J_SHDT,  SHDT*PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_TRHDT,TRHDT*PTYPE)
          CALL INC_AREG(I,J,JR,J_TRHDT,TRHDT*PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_LWCORR,dLWDT*PTYPE)
          CALL INC_AREG(I,J,JR,J_LWCORR,dLWDT*PTYPE)

        enddo
      enddo
      enddo

C****
C****
C**** SAVE SOME TYPE DEPENDENT FLUXES/DIAGNOSTICS
C****
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        JR=JREG(I,J)
        ptype = focean(i,j)
        if(ptype.gt.0.) then
            OA(I,J,5)=OA(I,J,5)+atmocn%solar(i,j)
            OA(I,J,6)=OA(I,J,6)+atmocn%trheat(i,j)
            OA(I,J,7)=OA(I,J,7)+atmocn%sensht(i,j)
            OA(I,J,8)=OA(I,J,8)+atmocn%latht(i,j)
            AIJ(I,J,IJ_FWOC)= AIJ(I,J,IJ_FWOC)+ atmocn%GMELT(I,J)*ptype
            AIJ(I,J,IJ_F0OC)= AIJ(I,J,IJ_F0OC)+atmocn%EGMELT(I,J)*ptype
            AIJ(I,J,IJ_MICB)=AIJ(I,J,IJ_MICB) +
     &           atmocn%GMELT(I,J)*ptype*AXYP(I,J)
            AIJ(I,J,IJ_EICB)=AIJ(I,J,IJ_EICB) +
     &           atmocn%EGMELT(I,J)*ptype*AXYP(I,J)
            CALL INC_AREG(I,J,JR,J_RVRD, atmocn%GMELT(I,J)*ptype)
            CALL INC_AREG(I,J,JR,J_ERVR,atmocn%EGMELT(I,J)*ptype)
#ifdef TRACERS_WATER
#ifdef TRACERS_OCEAN
            TAIJN(I,J,TIJ_ICB,:)=TAIJN(I,J,TIJ_ICB,:)+
     &           atmocn%TRGMELT(:,I,J)*ptype
#endif
#endif  /* TNL: inserted */
        endif
        ! open water
        ptype = asflx4(1)%ftype(i,j)
        if(ptype.gt.0.) then
            AIJ(I,J,IJ_POPWAT)=AIJ(I,J,IJ_POPWAT)+PTYPE
            if(focean(i,j).gt.0.)
     &           aij(i,j,ij_popocn) = aij(i,j,ij_popocn) + ptype
            AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+atmocn%evapor(i,j)*ptype
            IF (FOCEAN(I,J).gt.0) THEN
              AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)
     &             +atmocn%e0(i,j)*ptype
              AIJ(I,J,IJ_FWOC)=AIJ(I,J,IJ_FWOC)
     &             -atmocn%evapor(i,j)*ptype
            CALL INC_AJ(I,J,ITOCEAN,J_RVRD, atmocn%GMELT(I,J)*ptype)
            CALL INC_AJ(I,J,ITOCEAN,J_ERVR,atmocn%EGMELT(I,J)*ptype)
            END IF
        endif
        ! floating ice
        ptype = asflx4(2)%ftype(i,j)
        if(ptype.gt.0.) then
            OA(I,J,9)=OA(I,J,9)+atmice%trheat(i,j)
            OA(I,J,10)=OA(I,J,10)+atmice%sensht(i,j)
            OA(I,J,11)=OA(I,J,11)+atmice%latht(i,j)
            OA(I,J,12)=OA(I,J,12)+atmice%solar(i,j)
            AIJ(I,J,atmice%IJ_F0OI) =AIJ(I,J,atmice%IJ_F0OI)
     &           +atmice%e0(i,j)*ptype
            AIJ(I,J,IJ_EVAPI)=AIJ(I,J,IJ_EVAPI)+atmice%evapor(i,j)*ptype
            AIJ(I,J,IJ_SILWU)=AIJ(I,J,IJ_SILWU)+
     &           (DTsrc*TRHR(0,I,J)-atmice%trheat(i,j))*PTYPE
            AIJ(I,J,IJ_SILWD)=AIJ(I,J,IJ_SILWD)+DTsrc*TRHR(0,I,J)*PTYPE
            AIJ(I,J,IJ_SISH) =AIJ(I,J,IJ_SISH) -atmice%sensht(i,j)*ptype ! postive up
            IF (FOCEAN(I,J).gt.0) THEN
              CALL INC_AJ(I,J,ITOICE,J_RVRD,atmocn%GMELT(I,J)*ptype)
              CALL INC_AJ(I,J,ITOICE,J_ERVR,atmocn%EGMELT(I,J)*ptype)
            ENDIF
        endif
        ! land ice
        ptype = asflx4(3)%ftype(i,j)
        if(ptype.gt.0.) then
            AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)
     &           + atmgla%SENSHT(I,J)*PTYPE
            AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)
     &           + atmgla%LATHT(I,J)*PTYPE
            AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)
     &           + atmgla%TRHEAT(I,J)*PTYPE
            AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)
     &           + atmgla%E0(I,J)*PTYPE
            AIJ(I,J,IJ_EVAPLI)=AIJ(I,J,IJ_EVAPLI)
     &           + atmgla%EVAPOR(I,J)*PTYPE
            ! do the same extraction for sea/lake ice?
            AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+ptype
            AIJ(I,J,IJ_LI)  =AIJ(I,J,IJ_LI)  +ptype
            AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+atmgla%EPREC(I,J)*ptype
        endif
      enddo
      enddo

      return
      end subroutine surface_diag1a

      subroutine surface_diag2(moddd,ih,ihm,
     &     srhdt_sv,trhdt_sv,evhdt_sv,shdt_sv,evap_sv)
      USE GEOM, only : imaxj
      USE CONSTANT, only : tf
      USE FLUXES, only : atmsrf
      USE PBLCOM, only : dclev,ugeo,vgeo,bldep
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE DIAG_COM, ONLY :
     &      ndiuvar,ndiupt,ijdd,adiurn=>adiurn_loc
#ifndef NO_HDIURN
     *     ,hdiurn=>hdiurn_loc
#endif
     *     ,idd_ts,idd_qs,idd_us,idd_vs,idd_ws,idd_dbl,idd_cm,idd_ch
     *     ,idd_cq
     *     ,idd_tg1,idd_qg,idd_ug,idd_vg,idd_wg,idd_cia,idd_eds
     *     ,idd_swg,idd_lwg,idd_lh,idd_sh,idd_ev,idd_hz0,idd_dcf,idd_ldc
      implicit none
      integer, intent(in) :: ih,ihm
      integer, intent(in) :: moddd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     srhdt_sv,trhdt_sv,evhdt_sv,shdt_sv,evap_sv
c
      integer :: kr,ii,i,j, j_0, j_1, i_0,i_1
C**** some shorthand indices and arrays for diurn diags
      INTEGER :: idx6(24)
      REAL*8 :: tmp(NDIUVAR)
      real*8 :: srhdt,trhdt,shdt,evhdt,evap

C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
C**** Initialise constant indices
      idx6 = (/ idd_ts, idd_qs, idd_us, idd_vs, idd_ws, idd_dbl,
     &     idd_cm, idd_ch, idd_cq, idd_tg1, idd_qg,
     &     idd_ug, idd_vg, idd_wg, idd_cia, idd_eds,
     &     idd_swg, idd_lwg, idd_lh, idd_sh, idd_ev, idd_hz0,
     &     idd_dcf, idd_ldc
     &     /)


      CALL getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=j_0,j_1
      do i=i_0,imaxj(j)
        srhdt = atmsrf%solar(i,j)
        trhdt = atmsrf%trheat(i,j)
        shdt = atmsrf%sensht(i,j)
        evhdt = atmsrf%latht(i,j)
        evap = atmsrf%evapor(i,j)
        if(moddd.eq.0) then
          shdt_sv(i,j) = shdt-shdt_sv(i,j)
          evhdt_sv(i,j) = evhdt-evhdt_sv(i,j)
          srhdt_sv(i,j) = srhdt-srhdt_sv(i,j)
          trhdt_sv(i,j) = trhdt-trhdt_sv(i,j)
          evap_sv(i,j) = evap-evap_sv(i,j)
        else
          shdt_sv(i,j) = shdt
          evhdt_sv(i,j) = evhdt
          srhdt_sv(i,j) = srhdt
          trhdt_sv(i,j) = trhdt
          evap_sv(i,j) = evap
        endif
      enddo
      enddo


      IF(MODDD.EQ.0) THEN
        DO KR=1,NDIUPT
C**** For distributed implementation - ensure point is on local process.
          I = IJDD(1,KR)
          J = IJDD(2,KR)
          IF ((J >= J_0) .AND. (J <= J_1) .AND.
     &        (I >= I_0) .AND. (I <= I_1)) THEN
            tmp(:) = 0.
            tmp(idd_swg) = srhdt_sv(i,j)
            tmp(idd_lwg) = trhdt_sv(i,j)
            tmp(idd_lh) = evhdt_sv(i,j)
            tmp(idd_sh) = shdt_sv(i,j)
            tmp(idd_ev) = evap_sv(i,j)
            tmp(idd_hz0) =
     &           srhdt_sv(i,j)+trhdt_sv(i,j)+shdt_sv(i,j)+evhdt_sv(i,j)
            tmp(idd_ts) = atmsrf%tsavg(i,j)
            tmp(idd_qs) = atmsrf%qsavg(i,j)
            tmp(idd_qg) = atmsrf%qgavg(i,j)
            tmp(idd_us) = atmsrf%usavg(i,j)
            tmp(idd_vs) = atmsrf%vsavg(i,j)
            tmp(idd_ws) = atmsrf%wsavg(i,j)
            tmp(idd_dbl) = bldep(i,j) !dblavg(i,j)
            tmp(idd_ug) = ugeo(i,j)
            tmp(idd_vg) = vgeo(i,j)
            tmp(idd_wg) = sqrt(ugeo(i,j)**2+vgeo(i,j)**2)
            tmp(idd_cia) = atmsrf%ciaavg(i,j)
            tmp(idd_eds) = atmsrf%khsavg(i,j)

            tmp(idd_cm) = atmsrf%cmgs(i,j)
            tmp(idd_ch) = atmsrf%chgs(i,j)
            tmp(idd_cq) = atmsrf%cqgs(i,j)

            tmp(idd_tg1) = (atmsrf%gtemps(i,j)+tf)

            IF(DCLEV(I,J).GT.1.) THEN ! CHECK IF DRY CONV HAS HAPPENED
              tmp(idd_dcf)=1.
              tmp(idd_ldc)=DCLEV(I,J)
            endif
            ADIURN(idx6(:),kr,ih)=ADIURN(idx6(:),kr,ih)+tmp(idx6(:))
#ifndef NO_HDIURN
            HDIURN(idx6(:),kr,ihm)=HDIURN(idx6(:),kr,ihm)+tmp(idx6(:))
#endif
          END IF
        END DO 
      END IF
      end subroutine surface_diag2

      subroutine surface_diag3
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE GEOM, only : axyp,imaxj
      USE DIAG_COM, only : jreg,itlandi,itearth,aij=>aij_loc
     &     ,IJ_RUNLI,J_RUN,J_ACE1,J_ACE2,J_IMPLH,J_IMPLM,IJ_IMPMLI
     &     ,IJ_IMPHLI,IJ_F1LI,J_RSNOW,J_SNOW,IJ_RSIT,IJ_RSNW,IJ_SNOW
     &     ,IJ_ZSNOW
      USE LANDICE, only : ace1li,ace2li
      USE FLUXES, only : flice,atmgla,atmlnd,asflx4,atmsrf
      USE LANDICE_COM, only : mdwnimp,edwnimp
#ifdef TRACERS_WATER
     &     ,trdwnimp
#endif
      implicit none
      integer :: jr
      integer :: i,j, j_0, j_1, i_0,i_1
      real*8 :: plice,pearth

      CALL getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        JR=JREG(I,J)
        PLICE = FLICE(I,J)
        IF(PLICE.LE.0.) CYCLE
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+atmgla%RUNO(I,J)*PLICE
        CALL INC_AJ(I,J,ITLANDI,J_RUN ,  atmgla%RUNO(I,J)*PLICE)
        CALL INC_AREG(I,J,JR,J_RUN ,  atmgla%RUNO(I,J)*PLICE)

        CALL INC_AJ(I,J,ITLANDI,J_ACE1,ACE1LI*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_ACE2,ACE2LI*PLICE)
        CALL INC_AREG(I,J,JR,J_ACE1 ,ACE1LI*PLICE)
        CALL INC_AREG(I,J,JR,J_ACE2 ,ACE2LI*PLICE)

        CALL INC_AJ(I,J,ITLANDI,J_IMPLH,atmgla%IMPLH(I,J)*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_IMPLM,atmgla%IMPLM(I,J)*PLICE)
        CALL INC_AREG(I,J,JR,J_IMPLH, atmgla%IMPLH(I,J)*PLICE)
        CALL INC_AREG(I,J,JR,J_IMPLM, atmgla%IMPLM(I,J)*PLICE)
        AIJ(I,J,IJ_IMPMLI)=AIJ(I,J,IJ_IMPMLI)+atmgla%IMPLM(I,J)*PLICE
        AIJ(I,J,IJ_IMPHLI)=AIJ(I,J,IJ_IMPHLI)+atmgla%IMPLH(I,J)*PLICE

        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI)+atmgla%E1(I,J)*PLICE

C**** accumulate implicit fluxes for setting ocean balance
        MDWNIMP(I,J)=MDWNIMP(I,J)+atmgla%IMPLM(I,J)*PLICE*AXYP(I,J)
        EDWNIMP(I,J)=EDWNIMP(I,J)+atmgla%IMPLH(I,J)*PLICE*AXYP(I,J)
#ifdef TRACERS_WATER
        TRDWNIMP(:,I,J)=TRDWNIMP(:,I,J)+
     &       atmgla%IMPLT(:,I,J)*PLICE*AXYP(I,J)
#endif

      ENDDO
      ENDDO

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        JR=JREG(I,J)
        PLICE = asflx4(3)%ftype(i,j)
        IF(PLICE.GT.0.) THEN
          CALL INC_AJ(I,J,ITLANDI,J_RSNOW,PLICE*atmgla%SNOWFR(I,J))
          CALL INC_AJ(I,J,ITLANDI,J_SNOW,PLICE*atmgla%SNOW(I,J))
        ENDIF
        PEARTH = asflx4(4)%ftype(i,j)
        IF(PEARTH.GT.0.) THEN
          call inc_aj(i,j,itearth,j_rsnow,atmlnd%snowfr(i,j)*pearth)
          call inc_aj(i,j,itearth,j_snow, atmlnd%snow(i,j)*pearth)
          aij(i,j,ij_rsit)=aij(i,j,ij_rsit)+atmlnd%snowfr(i,j)*pearth
        ENDIF
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+atmsrf%SNOWFR(I,J)
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+atmsrf%SNOW(I,J)
        AIJ(I,J,IJ_ZSNOW)=AIJ(I,J,IJ_ZSNOW)+atmsrf%SNOWDP(I,J)
        CALL INC_AREG(I,J,JR,J_RSNOW, atmsrf%SNOWFR(I,J))
        CALL INC_AREG(I,J,JR,J_SNOW , atmsrf%SNOW(I,J))
      ENDDO
      ENDDO
      return
      end subroutine surface_diag3

#ifdef TRACERS_WATER
      subroutine water_tracer_evap(
     &     itype, i,j, n,
     &     tg1, rcdqws, rcdqdws, evap, snow, qg_sat, qsrf,
     &     lim_lake_evap, flake, tevaplim,
     &     lim_dew, dtsurf,
     &     trm1, trs, trgrnd, trgrnd2, trprime,
     &     trsrfflx, trevapor
     &     )
      use constant, only : teeny
      use landice, only : snmin
      implicit none
      integer :: itype,i,j,n
      real*8 :: tg1
      real*8 :: rcdqws
      real*8 :: rcdqdws
      real*8 :: evap
      real*8 :: snow
      real*8 :: qg_sat
      real*8 :: qsrf
      real*8 :: dtsurf
      real*8 :: flake
      real*8, intent(in) :: trm1	! comes from atmgla%
      real*8 :: trs
      real*8, intent(in) :: trgrnd	! comes from atmgla%
      real*8 :: trgrnd2
      real*8 :: trprime
      real*8 :: tevaplim
      real*8, intent(out) :: trsrfflx	! comes from atmgla%
      real*8, intent(out) :: trevapor	! comes from atmgla%


      logical :: lim_lake_evap, lim_dew
c
      real*8 ::
     &     tev,tevap,tdp
#ifdef TRACERS_SPECIAL_O18
     *     ,FRACVL,FRACVS,frac
#endif

C****
C**** Calculate Water Tracer Evaporation
C****
      IF (ITYPE.EQ.1) THEN      ! OCEAN
#ifdef TRACERS_SPECIAL_O18
        TEV=-(RCDQWS*(trs-trgrnd*QG_SAT*fracvl(tg1,n))
     *       +RCDQDWS*trprime)*pbl_args%frack(nx)
#else
        TEV=-(RCDQWS*(trs-trgrnd*QG_SAT)
     *       +RCDQDWS*trprime)
#endif
        TEVAP=DTSURF*TEV
C**** Limit evaporation if lake mass is at minimum
        IF (FLAKE.GT.0) THEN
#ifdef WATER_PROPORTIONAL
          if(lim_lake_evap) then
#else
          if( TREVAPOR+TEVAP.gt.TEVAPLIM ) THEN
#endif
c            IF(QCHECK) 
c     &           WRITE(99,*) "Lake TEVAP limited: I,J,TEVAP,TMWL"
c     *           ,N,TREVAPOR+TEVAP,TEVAPLIM
            TEVAP= TEVAPLIM-TREVAPOR
          endif
        END IF
      ELSE                      ! ICE AND LAND ICE
C**** tracer flux is set by source tracer concentration
        IF (EVAP.GE.0) THEN     ! EVAPORATION
          IF (EVAP.le.SNOW .or. SNOW.lt.SNMIN .or. ITYPE.ne.3) THEN
            TEVAP=EVAP*trgrnd
          ELSE                  ! special treatment for landice when EVAP>SNOW>SNMIN
            TEVAP=SNOW*(trgrnd-trgrnd2)+EVAP*trgrnd2
          END IF
        ELSE                    ! DEW (fractionates)
#ifdef TRACERS_SPECIAL_O18
          IF (TG1.gt.0) THEN
            frac=FRACVL(TG1,n)
          ELSE
            frac=FRACVS(TG1,n)
          END IF
          TEVAP=EVAP*trs/(QSRF*frac+teeny)
#else
          TEVAP=EVAP*trs/(QSRF+teeny)
#endif
        END IF
      END IF
      TDP = TEVAP
#ifdef WATER_PROPORTIONAL
      if(lim_dew) then
#else
      IF (TRM1+TDP.lt.0..and.tdp.lt.0) THEN
#endif
c        IF (QCHECK)
c     &       WRITE(99,*) "LIMITING TRDEW",I,J,N,TDP,TRM(I,J,1,n),TDT1
        TEVAP = -TRM1
        TDP = -TRM1
      END IF
      trsrfflx = TDP/DTSURF
      TREVAPOR = TREVAPOR + TEVAP
      return
      end subroutine water_tracer_evap
#endif

#ifdef TRACERS_ON
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      subroutine collect_ocean_emissions(i,j,ptype,dtsurf,pbl_args)
      use OldTracer_mod, only : trname
      use fluxes, only : atmocn
      use geom, only : axyp
      use trdiag_com, only : taijs=>taijs_loc,ijts_isrc,jls_isrc
      use pbl_drv, only : t_pbl_args
      use trdiag_com, only : itcon_surf
#ifdef TRACERS_TOMAS
      use tracer_com, only : xk,nbins
#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: ptype,dtsurf
      type (t_pbl_args) :: pbl_args
c
      integer :: n,nx
      real*8 :: trc_flux
#ifdef TRACERS_TOMAS
      INTEGER ss_bin,num_bin
      real*8 ss_num(nbins)
c$$$      real*8, parameter :: scalesizeSS(nbins)=(/!0.0,0.0,0.0,
c$$$     *     6.4614E-08,5.0110E-07,2.7243E-06,1.1172E-05,
c$$$     *     3.7192E-05,1.2231E-04,4.4986E-04,1.4821E-03,
c$$$     *     3.7403E-03,7.9307E-03,1.8918E-01,7.9705E-01/)
#endif
c
#ifdef TRACERS_TOMAS
      ss_bin=0
      num_bin=0
#endif
C**** Loop over tracers
      DO NX=1,pbl_args%NTX
        N=pbl_args%NTIX(NX)

C****
C**** Calculate Aersosol Exchange
C****
        select case (trname(n))
        case ('DMS')
          trc_flux=pbl_args%DMS_flux
        case ('seasalt1', 'M_SSA_SS')
          trc_flux=pbl_args%ss1_flux
        case ('seasalt2', 'M_SSC_SS')
          trc_flux=pbl_args%ss2_flux
        case ('M_SSS_SS')
          trc_flux=(pbl_args%ss1_flux+pbl_args%ss2_flux)
#ifdef TRACERS_AEROSOLS_OCEAN
        case ('OCocean')
          trc_flux=pbl_args%OCocean_flux
#endif  /* TRACERS_AEROSOLS_OCEAN */
        case default
          trc_flux=0
#ifdef TRACERS_TOMAS
        case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04', 
     &         'ANACL_05','ANACL_06','ANACL_07','ANACL_08',
     &         'ANACL_09','ANACL_10','ANACL_11','ANACL_12'
#ifdef TOMAS_12_3NM 
     *    ,'ANACL_13','ANACL_14','ANACL_15'
#endif
     &         )
          ss_bin=ss_bin+1
          trc_flux=pbl_args%tomas_ss_flux(ss_bin)

          ss_num(ss_bin)=(pbl_args%tomas_ss_flux(ss_bin))
     &         /sqrt(xk(ss_bin)*xk(ss_bin+1))
! No subgrid coagulation for sea-salt
!        TOMAS_EMIS(I,J,ss_bin,1)= trc_flux*axyp(i,j)*ptype

        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04',
     &         'ANUM__05','ANUM__06','ANUM__07','ANUM__08',
     &         'ANUM__09','ANUM__10','ANUM__11','ANUM__12'
#ifdef TOMAS_12_3NM 
     *    ,'ANUM__13','ANUM__14','ANUM__15'
#endif
     &         )
           num_bin=num_bin+1
           trc_flux=ss_num(num_bin)
      
#endif
        end select

        atmocn%trsrfflx(n,i,j)=atmocn%trsrfflx(n,i,j)+trc_flux
        if (ijts_isrc(1,n)>0) then
           taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &       trc_flux*axyp(i,j)*ptype*dtsurf
        end if

#ifdef TRACERS_AMP
            select case (trname(n))
              case ('DMS','M_SSA_SS','M_SSC_SS','M_SSS_SS')
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)
            end select
#else
#ifndef TRACERS_TOMAS

        if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols?
#endif
#endif

      
#ifdef TRACERS_TOMAS

        select case (trname(n))

            case ('DMS')              
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *              trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)

        case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04', 
     &       'ANACL_05','ANACL_06','ANACL_07','ANACL_08',
     &       'ANACL_09','ANACL_10','ANACL_11','ANACL_12'
#ifdef TOMAS_12_3NM 
     *    ,'ANACL_13','ANACL_14','ANACL_15'
#endif
     &         )
        
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)
        
        if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols?
        
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04',
     &       'ANUM__05','ANUM__06','ANUM__07','ANUM__08',
     &       'ANUM__09','ANUM__10','ANUM__11','ANUM__12'
#ifdef TOMAS_12_3NM 
     *    ,'ANUM__13','ANUM__14','ANUM__15'
#endif
     &         )
        
!TOMAS - itcon_surf (1,3) is for SO4/EC/OC.
        if (itcon_surf(4,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(4,n),n)

        if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols? 

            end select
#endif

      enddo ! tracer loop
      return
      end subroutine collect_ocean_emissions
#endif
#endif

#ifdef TRACERS_GASEXCH_ocean
      subroutine calc_gasexch(i,j,itype,ns,moddsf,
     &     ptype,pocean,rsi,rhosrf,tgo,
     &     dtsurf,pbl_args)
      use TimeConstants_mod, only: SECONDS_PER_YEAR
      use fluxes, only : focean,atmocn,nisurf
      use geom, only : axyp
      use trdiag_com, only : taijs=>taijs_loc,
     &     ijts_isrc,ijts_gasex,jls_isrc
      use pbl_drv, only : t_pbl_args
      use model_com, only : itime,dtsrc
      use diag_com, only : aij=>aij_loc,ij_gasx,ij_kw,ij_alpha
      use OldTracer_mod, only: vol2mass, tr_mm, itime_tr0
      USE TRACER_COM, only: n_co2n,n_cfcn
#ifdef OBIO_ON_GARYocean
      USE MODEL_COM, only: nstep=>itime
#else
      USE HYCOM_SCALARS, only: nstep
#endif
      implicit none
      integer, intent(in) :: i,j,itype,ns,moddsf
      real*8, intent(in) :: ptype,pocean,rsi,rhosrf,tgo,dtsurf
      type (t_pbl_args) :: pbl_args
c
      real*8 :: trgrnd,alpha_gas2,trs
      real*8, dimension(:,:,:), pointer :: TRGASEX
      integer :: n,nx


      TRGASEX => atmocn%TRGASEX

C****
C**** Calculate Tracer Gas Exchange
C****
      IF (ITYPE.EQ.1 .and. focean(i,j).gt.0.) THEN ! OCEAN
        DO NX=1,pbl_args%NTX
          N=pbl_args%NTIX(NX)
          trgrnd=atmocn%gtracer(n,i,j)
          trs=pbl_args%trs(nx)
#ifdef TRACERS_GASEXCH_ocean_CFC
          TRGASEX(n,I,J) = TRGASEX(n,I,J) +
     .        pbl_args%Kw_gas * (pbl_args%beta_gas*trs-trgrnd)
          atmocn%trsrfflx(n,i,j) = atmocn%trsrfflx(n,i,j)
     .         -pbl_args%Kw_gas * (pbl_args%beta_gas*trs-trgrnd)
          taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n))
     .         -pbl_args%Kw_gas * (pbl_args%beta_gas*trs-trgrnd)
     .               * axyp(i,j)*ptype*dtsurf
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
! TRGASEX is the gas exchange flux btw ocean and atmosphere.
! Its sign is positive for flux entering the ocean (positive down)
! because obio_carbon needs such. Units mol,CO2/m2/s (accumulated over itype)

          TRGASEX(n,I,J) = TRGASEX(n,I,J) +
     .          pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs 
     .         - pbl_args%alpha_gas * trgrnd ) 
     .         * 1d6/vol2mass(n)
     .         * dtsurf/dtsrc      !in order to accumulate properly over time
     .         * (1.d0-RSI)   !units mol,co2/m2/s

          if (nstep.eq.48)
     .         write(*,'(a,3i5,6e12.4)')'SURFACE: test carbon flux: ',
     .         nstep,i,j,
     .         pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs
     .         - pbl_args%alpha_gas * trgrnd )
     .         * 1d6/vol2mass(n),
     .         dtsurf,dtsrc,1.d0-RSI,TRGASEX(n,I,J),
     .         focean(i,j)

! trsrfflx is positive up 
! units are kg,CO2/s
          atmocn%trsrfflx(n,i,j)=atmocn%trsrfflx(n,i,j)
     .         - pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs 
     .         - pbl_args%alpha_gas * trgrnd )
     .         * 1.0d6/vol2mass(n) 
     .         * tr_mm(n)*1.0d-3        !units kg,co2/m2/s

!units are kg,co2
          taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n))-
     $         pbl_args%Kw_gas * (pbl_args%beta_gas* trs-
     $         pbl_args%alpha_gas * trgrnd )* 1.0d6/vol2mass(n)
     $         * tr_mm(n)*1.0d-3* ptype* axyp(i,j) * dtsurf

          if(i.eq.1 .and. j.eq.45) then
            write(*,'(a,3i5,11e12.4)')'SURFACE, trgasex:',
!     write(*,'(a,3i5,11e12.4)')'22222222222222222',
     .           nstep,i,j,pbl_args%Kw_gas,pbl_args%beta_gas,trs,
     .           pbl_args%alpha_gas,trgrnd,TRGASEX(n,I,J),
     .      pbl_args%Kw_gas * pbl_args%beta_gas*trs*1.d6/vol2mass(n)
     .           *ptype,
     .           pbl_args%Kw_gas * pbl_args%alpha_gas * trgrnd 
     .           * 1.0d6/vol2mass(n) * ptype,
     .           atmocn%trsrfflx(n,i,j),rhosrf,taijs(i,j,ijts_isrc(1,n))
          endif
#endif
        END DO
      END IF

      do nx=1,pbl_args%ntx
        n=pbl_args%ntix(nx)
        trgrnd=atmocn%gtracer(n,i,j)
        trs=pbl_args%trs(nx)
c      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if (n.eq.n_CO2n .or. n.eq.n_CFCn) then
          if (focean(i,j).gt.0) then

! original versions
             if (POCEAN.gt.0) then ! original coding
               if (MODDSF.EQ.0) THEN
                 AIJ(i,j,ij_kw) = AIJ(i,j,ij_kw)  
     .            + pbl_args%Kw_gas * focean(i,j) ! m/s
     .            * (1.d0 - RSI) ! only over open water
                 AIJ(i,j,ij_alpha) = AIJ(i,j,ij_alpha) 
     .           + pbl_args%alpha_gas * focean(i,j) ! mol,CO2/m3/uatm
     .           * (1.d0 - RSI) ! only over open water
               endif
               if(NS==NIsurf .and. itype==1) then
                 AIJ(i,j,ij_gasx) = AIJ(i,j,ij_gasx) 
     .                + TRGASEX(n,I,J) * focean(i,j)
     .                * SECONDS_PER_YEAR    ! mol,CO2/m2/yr
     .                * (1.d0 - RSI) ! only over open water
               endif
            end if
! tracer diag versions
            if (ITYPE.eq.1) then 
! gas exchange
            taijs(i,j,ijts_gasex(3,n)) = taijs(i,j,ijts_gasex(3,n)) 
     .           + pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs 
     .         - pbl_args%alpha_gas * trgrnd ) 
     .         * 1d6/vol2mass(n)
     .         * dtsurf/dtsrc   !in order to accumulate properly over time
     .         * ptype * SECONDS_PER_YEAR        ! mol/m2/yr

! zonal mean diag accumulates kgCO2
            if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *           - pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs 
     .         - pbl_args%alpha_gas * trgrnd ) 
     .         * 1d6/vol2mass(n) * dtsurf  
     .         * ptype*tr_mm(n)*1d-3*axyp(i,j))

              if (MODDSF.EQ.0) THEN
! piston velocity
                taijs(i,j,ijts_gasex(1,n)) = taijs(i,j,ijts_gasex(1,n)) 
     .                + pbl_args%Kw_gas * pocean ! m/s only over open water
! solubility mol/m3/uatm
                taijs(i,j,ijts_gasex(2,n)) = taijs(i,j,ijts_gasex(2,n)) 
     .               + pbl_args%alpha_gas * focean(i,j) 
              endif

            elseif (POCEAN.eq.0) then  ! ITYPE=2 and all ice covered
! solubility mol/m3/uatm ice covered area
               if (MODDSF.EQ.0) taijs(i,j,ijts_gasex(2,n)) = taijs(i,j
     $              ,ijts_gasex(2,n))+ alpha_gas2(tgo,pbl_args%sss_loc)
     $              * focean(i,j) 
            endif                ! itype
          endif                  !focean
          endif                  !gasexch tracers
        end if
      end do 
      return
      end subroutine calc_gasexch
#endif /* TRACERS_GASEXCH_ocean */
