! Lakes = itype == 1
! Look in GHY_DRV to see which tracers in SURFACE.f are needed for LANDICE surface type

C****   
C**** SURFACE.f    SURFACE fluxes    2006/12/21
C****
#include "rundeck_opts.h"



      SUBROUTINE SURFACE_LANDICE (do_init,moddd,dtsurf,atmgla,ihc)

!@sum SURFACE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : rgas,lhm,lhe,lhs
     *     ,sha,tf,rhow,shv,shi,stbo,bygrav,by6
     *     ,deltx,teeny,grav
      USE MODEL_COM, only : modelEclock
      USE MODEL_COM, only : itime
#ifdef SCM
      USE SCMDIAG, only : EVPFLX,SHFLX
      USE SCMCOM, only : iu_scm_prt, ALH, ASH, SCM_SURFACE_FLAG
     &     ,I_TARG,J_TARG
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE GEOM, only : imaxj
#ifdef TRACERS_ON
      use OldTracer_mod, only: itime_tr0, needtrs
      USE TRACER_COM, only : NTM
#ifdef TRACERS_WATER
      use OldTracer_mod, only: nWATER,tr_wd_TYPE
#endif
#endif
      USE PBL_DRV, only : alloc_pbl_args, dealloc_pbl_args
      USE PBL_DRV, only : pbl, t_pbl_args, xdelt

      USE LANDICE, only : z1e,z2li,hc1li,hc2li,ace1li,ace2li,snmin
      USE LANDICE_COM, only : snowli
#ifdef TRACERS_WATER
     *     ,trlndi
#endif
      USE LANDICE_COM, only : ijhc,ijhc_tsurf,ijhc_tsli
      USE SEAICE, only : xsi,ace1i,alami0,rhoi,byrls,alami
      USE EXCHANGE_TYPES
      USE Timer_mod, only: Timer_type
      USE TimerList_mod, only: startTimer => start
      USE TimerList_mod, only: stopTimer => stop
      USE itype_enum

      IMPLICIT NONE
! ================ Parameter Declarations =====================
      logical, intent(in) :: do_init
      integer, intent(in) :: moddd,ihc
      real*8, intent(in) :: dtsurf
      type(atmgla_xchng_vars) :: atmgla
! ================ VARIABLE DECLARATIONS ======================
      INTEGER I,J,K,ITYPE,IH,IHM
      REAL*8 PLICE
      REAL*8 PS
     *     ,ELHX,MSI2,CDTERM,CDENOM,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG
     *     ,dT2,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER
     *     ,Q1,THV1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,RCDHDWS,RCDQDWS
     *     ,SHEAT,TRHEAT,T2DEN,T2CON,T2MUL,FQEVAP,Z1BY6L,F2
     *     ,FSRI(2),dlwdt,TGO

      REAL*8 MA1
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: Z2LI3L=Z2LI/(3.*ALAMI0), Z1LIBYL=Z1E/ALAMI0
      REAL*8, EXTERNAL :: QSAT
      REAL*8 DQSATDT,TR4
c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 qg_sat,qsrf,us,vs,ws,ws0

      ! This is for making a correction to the surface wind stress
      ! (which is proportional to |u_air - u_ocean|, or  | u_air - u_seaice |)
      ! For land ice, the difference in magnitudes makes the correction negligible.
      ! Hence, we will have always uocean == vocean == 0
c      real*8 uocean, vocean
c      logical pole
c
      logical :: lim_dew ! for tracer convenience
#ifdef TRACERS_ON
      integer n,nx
c      real*8, dimension(ntm) :: trs!,trsfac,trconstflx
      integer ntix(ntm), ntx
      real*8, dimension(ntm) :: trgrnd,trgrnd2
#endif
C****
      INTEGER :: J_0, J_1, J_0H, J_1H, I_0,I_1
      LOGICAL :: debug

      real*8, dimension(:,:), pointer :: trdn,srdn,e1,tgrnd,tgr4

      type (Timer_type), pointer :: aTimer
! ======================= MAIN BODY ======================================

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,        J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C****

      call startTimer('SURFACE_LANDICE()')

      trdn => atmgla%flong
      srdn => atmgla%fshort
      e1 => atmgla%e1
      tgrnd => atmgla%tgrnd
      tgr4 => atmgla%tgr4

! inits before looping over surface timesteps.
      if(do_init) then
        do j=j_0,j_1
        do i=i_0,i_1
          e1(i,j) = 0.
          tgrnd(i,j) = atmgla%gtemp(i,j)
          tgr4(i,j) = atmgla%gtempr(i,j)**4
        enddo
        enddo
      endif

      IH=modelEclock%hour()+1
      IHM = IH+(modelEclock%date()-1)*24
c avoid uninitialized variable problems when the first gridpoint
c in the domain is ocean
      SNOW = 0.

C**** Zero out fluxes summed over type and surface time step
c ---- Don't need this section, it's asflx is a global variable
c and it's already zeroed out before outer loop

!      E1=0.    ! I don't think local var E1 is really used.  It is set, but not read.
! /------------ All global variables, zeroed in SURFACE.f ------------\
!#ifdef SCM
!      EVPFLX= 0.0d0
!      SHFLX = 0.0d0
!#endif
! \--------------------------------------------------------------------/
! Next in SURFACE.f comes outside loop over timesteps.

! ==============================================================

      call alloc_pbl_args(pbl_args)

      ITYPE=3

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
      pbl_args%ntx = ntx
      pbl_args%ntix(1:ntm) = ntix(1:ntm)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      pbl_args % moddd = moddd
      pbl_args % ih = ih
      pbl_args % ihm = ihm
#endif

C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

      ! Model may expect that some things are zero
      ! call stop_model('Please double-check that model might expect some things to be defined for some boxes, even if there is no ice in that box.', -17)
      PLICE=atmgla%ftype(I,J)
      PTYPE = PLICE
      IF (PTYPE <= 0) CYCLE

C****
C**** DETERMINE SURFACE CONDITIONS
C****

      ! PS = Pressure at Surface
      ! PEDN = Pressure at lower edge of box (incl. surface) (mb)  (ATM_COM.f)
      !PS=PEDN(1,I,J)
      PS=atmgla%srfp(I,J)

      ! PEK = PEDN**KAPA
      ! kapa = ideal gas law exponent for dry air (.2862)
      ! kapa = (g-1)/g where g = 1.401 = c_p/c_v
      !     [or: kapa = (c_p - c_v) / c_p ]
      ! g (or srat) = c_p/c_v = ratio of specific heats at const. press. and vol. (=1.401)
      !PSK=atmgla%srfpk(i,j)

      ! Q = specific humidity (kg water vapor/kg air)
      !     (This is the same as mixing ratio)
      ! Q1 = Specific humidity of the bottom atmosphere layer
      ! TODO: This needs to be adjusted for downscaling
      Q1=atmgla%Q1(I,J) !Q(I,J,1)

      ! T(I,J,Z) = Dry Potential Temperature (theta) referenced to 1 millibar
      ! Actual Temperature temp = T*pk = T * [(p/p0)**(R/cp)]
      !      (p = pressure, p0 = 1 millibar,
      !       R is gas constant of air (NOT the Universal Gas Constant),
      !       cp is specific heat capacity at constant pressure.) 
      ! See: http://en.wikipedia.org/wiki/Potential_temperature
      !      http://en.wikipedia.org/wiki/Gas_constant
      ! If T is a temperature, then virtual temperature VT can be approximated with:
      !      VT = T * (1 + Q*epsilon)
      !      where epsilon = R_d/R_v = M_v/M_d ~= .622 in Earth's Atmosphere
      !      R_d = gas constant for air; R_v = gas constant for water vapor
      ! THV1 = Virtual Potential Temperature.
      ! deltx = epsilon = coeff. of humidity in virtual temperature defn. (0.6078)
      ! xdelt = deltx (or 0, if we want to define PBL code to take and receive
      !     actual temperatures to the external workd).
      ! Note from Max Kelley:
      !     Virtual potential temperature is for computing the static stability of
      !     the atmosphere (usually for turbulence parameterizations).  My
      !     understanding is that the new code to be introduced for downscaling
      !     probably does not need to re-compute any turbulence-related quantities
      !     differently than the model is already doing.  So just compute the
      !     virtual temperature passed to the PBL using co-located T and Q,
      !     getting Q for the height classes with the same interpolation used for
      !     T I guess.
      !
      !     But keep in mind that virtual temperature effects are negligible in
      !     "cold" regions of the globe, and that the parameter xdelt is actually
      !     set to zero these days (virtual temp. effects on buoyancy are
      !     considered only _within_ the PBL code).
      ! TODO: THV1 needs to get downscaled for topography
      !THV1=T(I,J,1)*(1.+Q1*xdelt)
      !THV1=atmgla%TEMP1(I,J)*(1.+Q1*xdelt)

      ! MA1 = mass of lowest atmospheric layer (kg/m^2)
      MA1=atmgla%AM1(I,J)

C****

      ! This is a good way to make sure you dont forget to do something
      ! before you run the model.
      ! call stop_model('Please double-check something or another.', -17)

! BEGIN ---------------------------------------------------------
      PTYPE=PLICE
      ! snow amount on land ice (kg/m^2)
      SNOW=SNOWLI(I,J,IHC)
      ! Temperature of top ice layer (C)
      TG1=TGRND(I,J)

      ! TR4 = TGR4(3,i,j) = atmgla%GTEMPR**4
      ! (Needed for Stefan-Boltzmann Law)
      ! GTEMPR radiative ground temperature over surface type (K)
      TR4=TGR4(I,J)

      ! SRHEAT = Solar Heating
      ! FSF = Solar Forcing over each type (W/m^2)
      ! FSF = net absorption (for cosZ = 1)
      ! COSZ1 = Mean Solar Zenith angle for curr. physics(not rad) time step
      SRHEAT=SRDN(I,J)*atmgla%COSZ1(I,J)

      ! LHS = latent heat of sublimation at 0 C (J/kg)
      ELHX=LHS

c      uocean = 0. ; vocean = 0. ! no land ice velocity
#ifdef TRACERS_WATER
      do nx=1,ntx
        trgrnd2(nx)=TRLNDI(ntix(nx),I,J,IHC)/(ACE1LI+ACE2LI)
      end do
      pbl_args%trgrnd2(1:ntm) = trgrnd2(1:ntm)
#endif
! END ---------------------------------------------------------

! ----------------------------------------------------------
! This code is replicated when factorut of SURFACE.f
! Bob should not have to touch this code, it is not LANDICE-specific.

C****
C****
C**** BOUNDARY LAYER INTERACTION
C****
      ! TG = Temperature of top ice layer (K)
      ! TF = freezing point of water at 1 atm (273.16 K)
      TG=TG1+TF

      ! LHS = latent heat of sublimation at 0 C (J/kg)
      ! QG_SAT = Saturation vapor mixing ratio (kg vapor / kg air in a given volume)
      ! PS = Pressure (TODO: Must change with height classes)
      QG_SAT=QSAT(TG,ELHX,PS)

!      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG_SAT=0.98d0*QG_SAT
      pbl_args%TG=TG   ! actual ground temperature
      pbl_args%TR4=TR4 ! radiative temperature K^4
      !pbl_args%ELHX=ELHX   ! relevant latent heat
      !pbl_args%QSOL=SRHEAT   ! solar heating
      pbl_args%TGV=TG*(1.+QG_SAT*xdelt)  ! Virtual temperature of the ground

C =====================================================================
      pbl_args%dtsurf = dtsurf

      ! TKV = Virtual temperature at the surface
      ! TKV = ttop = virtual potential temperature at the top of the
      !              Planetary Boundary layer (PBL).
      !              (if xdelt=0, ttop is the actual temperature)
      !pbl_args%TKV=THV1*PSK     ! TKV is referenced to the surface pressure

      !pbl_args%ZS1=.5d-2*RGAS*pbl_args%TKV*MA1/atmgla%p1(i,j) !PMID(1,I,J)
      ! TODO: qg_sat changes with height (this will come out automatically)
      pbl_args%qg_sat = qg_sat
      !pbl_args%qg_aver = qg_sat   ! QG_AVER=QG_SAT
      !pbl_args%hemi = sign(1d0,atmgla%lat(i,j))
c      pbl_args%pole = pole
c      pbl_args%evap_max = 1.
c      pbl_args%fr_sat = 1. ! entire surface is saturated
c      pbl_args%uocean = uocean
c      pbl_args%vocean = vocean
      ! TODO: This will change with height
c      pbl_args%psurf = PS
c      pbl_args%trhr0 = TRDN(I,J)
      pbl_args%ocean = .false.
      pbl_args%snow = SNOW

! Calculate drag coefficients, wind speed, air density, etc.
! PBL = "Planetary Boundary Layer"
C**** Call pbl to calculate near surface profile
      CALL PBL(I,J,ITYPE,PTYPE,pbl_args,atmgla)

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
      ! TG comes from the landice model.
      TG = TG + pbl_args%dskin
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (pbl_args%ocean) QG_SAT=0.98d0*QG_SAT
      TG1 = TG - TF
      TR4=(sqrt(sqrt(TR4))+pbl_args%dskin)**4
! %dskin trapped here!!!

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
      ! TS = Temperature @ surface air (10m)
      ! TG = Temperature of ground
      SHEAT=SHA*(RCDHWS*(TS-TG)+RCDHDWS*pbl_args%tprime)
      ! Including gustiness in the latent heat flux:
      EVHEAT=(LHE+TG1*SHV)*(RCDQWS*(QSRF-QG_SAT)+
     *                      RCDQDWS*pbl_args%qprime)
      TRHEAT=TRDN(I,J)-STBO*TR4

! BEGIN ------------------------------------------------------------
! Setting up implicit timestep aspects to surface flux calculations
! Setting up derivitives of fluxes w.r.t temperature, etc.
! Might need to change these if I decide to use layering, differnt
! variable names, etc.
C**** CASE (3) ! FLUXES USING IMPLICIT TIME STEP OVER LANDICE

      ! Z1E = Thickness of top landice layer (.1m const)
      ! ALAMI0 = Lambda coefficient for ice J/(m*degC*sec) = W/(m K) (2.11 const)
      !     (Lambda is thermal conductivity)
      ! Z1LIBYL = "Z1 Land Ice by Lambda" = Z1E / ALAMI0
      ! SNOW = snow amount on land ice (kg/m^2)
      !                                          +---- m -------+
      !                   m     J/(m*degC*sec)   kg/m^2    kg/m^3   J/(m*degC*sec)
      ! Z1BY6L = 1/6 * [(Z1E / ALAMI0)          + ((SNOW / RHOS)   /   ALAMS)]
      ! Z1BY6L units = (m^2 * degC * sec) / J
      ! Z1BY6L = "Z1 by 6 Lambda"
      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6

      ! Temperature of second ice layer
      TG2=atmgla%GTEMP2(I,J)

      ! TG2 = Temperature of second ice layer
      CDTERM=TG2

      ! Z2LI = Thickness of second layer of ice (2.9m)
      ! Z2LI3L = "Z2LI by 3 Lambda" = 1/3 (Z2LI / ALAMI0)
      ! Z2LI3L=Z2LI/(3.*ALAMI0)
      ! 2 * Z1BY6L = "Z1 by 3 Lambda"
      ! CDENOM = 3 Lambda(ice) / (Z1 + Z2 = 3m)
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)

      ! HC1LI = heat capacity of first layer land ice (J/m^2)
      ! SNOW = snow amount on land ice (kg/m^2)
      ! SHI = heat capacity of pure ice (at 0 C) (2060 J/kg C)
      ! SNOW*SHI = Heat capacity of snow layer (J/m^2)
      HCG1=HC1LI+SNOW*SHI

c      SHDT=0.
c      EVHDT=0.
c      TRHDT=0.
c      F1DT=0.

        ! F0 = Energy flux between atmosphere and surface
        ! SRHEAT = Solar Heating
        ! TRHEAT = Thermal-Radiation (longwave; down?)
        ! SHEAT = Sensible Heat
        ! EVHEAT = Latent Heat
        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT

        ! F1 = Energy flux between first and second ice layers
        ! TG1 = Temperature of top ice layer
        ! CDTERM = TG2 = Temperature of second ice layer
        ! F0 = Energy flux between atmosphere and surface
        ! Z1BY6L =
        ! CDENOM =
        F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM

        DSHDTG=-RCDHWS*SHA
        DQGDTG=QG_SAT*DQSATDT(TG,ELHX)
        DEVDTG=-RCDQWS*LHE*DQGDTG
        DTRDTG=-4.*STBO*sqrt(sqrt(TR4))**3
        DF0DTG=DSHDTG+DEVDTG+DTRDTG
        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)
        SHDT=DTSURF*(SHEAT+DTG*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+DTG*DEVDTG)
        TRHDT=DTSURF*(TRHEAT+DTG*DTRDTG)
        F1DT=DTSURF*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
        TG1=TG1+DTG

! END ------------------------------------------------------------

C**** CALCULATE EVAPORATION
      DQ1X =EVHDT/((LHE+TG1*SHV)*MA1)
      EVHDT0=EVHDT
      lim_dew=.false.


      IF (DQ1X.GT.Q1) THEN
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
     &         3, i,j, n,
     &         tg1, rcdqws, rcdqdws, evap, snow, qg_sat, qsrf,
     &         .false., 0d0, 0d0, ! arguments for lakes only
     &         lim_dew, dtsurf,
     &         atmgla%TRM1(n,I,J), pbl_args%trs(nx),
     &         atmgla%gtracer(n,i,j), trgrnd2(nx),
     &         pbl_args%trprime(nx),
     &         atmgla%trsrfflx(n,i,j), atmgla%trevapor(n,i,j)
     &     )
        END IF
      END DO
#endif


C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
      atmgla%latht(i,j) = atmgla%latht(i,j) + evhdt
      atmgla%trheat(i,j) = atmgla%trheat(i,j) + trhdt
      atmgla%solar(i,j) = atmgla%solar(i,j) + dtsurf*srheat

      ! atmgla%e0 = net energy flux at surface (J/m^2)
      atmgla%E0(I,J)=atmgla%E0(I,J)+F0DT

      ! atmgla%e1 <= e1 (set below in SURFACE.f)
      ! E1 = net energy flux at layer 1
      ! (I believe this means between the top 10cm layer and the 2.9m layer below it)
      E1(I,J)=E1(I,J)+F1DT

      atmgla%EVAPOR(I,J)=atmgla%EVAPOR(I,J)+EVAP
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
      TGRND(I,J)=TG1  ! includes skin effects
      TGR4(I,J) =TR4
C**** calculate correction for different TG in radiation and surface
      dLWDT = DTSURF*(atmgla%TRUP_in_rad(I,J)-TRDN(I,J))+TRHDT
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
      atmgla%DTH1(I,J)=-(SHDT+dLWDT)/(SHA*MA1) ! +ve up
      atmgla%sensht(i,j) = atmgla%sensht(i,j)+SHDT
      atmgla%DQ1(I,J) = -DQ1X
#ifdef SCM
      if (i.eq.I_TARG.and.j.eq.J_TARG) then
          write(iu_scm_prt,988) I,PTYPE,DTH1(I,J),DQ1(I,J),SHDT,dLWDT
 988      format(1x,'988 SURFACE GCM  I PTYPE DTH1 DQ1 SHDT dLWDT ',
     &           i5,f9.4,f9.5,f9.6,f12.4,f10.4)
      endif
      endif
#endif
!unused      DMUA_IJ=RCDMWS*US
!unused      DMVA_IJ=RCDMWS*VS
!unused      atmgla%DMUA(I,J) = atmgla%DMUA(I,J) + DMUA_IJ*DTSURF
!unused      atmgla%DMVA(I,J) = atmgla%DMVA(I,J) + DMVA_IJ*DTSURF
      atmgla%uflux1(i,j) = RCDMWS*US
      atmgla%vflux1(i,j) = RCDMWS*VS


      ! demo diagnostic: a PBL output
      ijhc(i,j,ihc,ijhc_tsurf) = ijhc(i,j,ihc,ijhc_tsurf)
     &     +ts*dtsurf  ! scale factor = dtsurf/dtsrc
      ijhc(i,j,ihc,ijhc_tsli) = ijhc(i,j,ihc,ijhc_tsli)
     &     +ts*dtsurf  ! scale factor = dtsurf/dtsrc

C****

      END DO  ! end of I loop
      END DO  ! end of J loop
! ============================================================
! Now we're outside the loop over grid points

      call dealloc_pbl_args(pbl_args)

      atmgla%gtemp(:,:) = tgrnd(:,:)
      atmgla%gtemps(:,:) = tgrnd(:,:)

      call stopTimer('SURFACE_LANDICE()')

      RETURN
C****
      END SUBROUTINE SURFACE_LANDICE

! -----------------------------------------------------------------------
      subroutine patchify_landice_inputs(phase)
! create patch-specific values of inputs to the land ice model from
! grid-mean values and patch-specific physical parameters (elevation etc.)
! This initial version is simply copying grid-mean values into all patches.
      use fluxes, only : atmglas,atmgla
     &     ,after_atm_phase1,during_srfflx
      use domain_decomp_atm, only : grid
      use exchange_types, only : atmgla_xchng_vars

      ! Stuff needed for downscaling
      use landice_com, only : elevhc, HC_T_LAPSE_RATE
      use constant, only : kapa, Grav
      use atm_com, only : zatmo
      USE GEOM, only : imaxj
      USE RESOLUTION, only : ptop
      USE DYNAMICS, only : dsig
      USE CONSTANT, only : LHS

      implicit none
      real*8, external :: QSAT      ! Import function from Utilities.F90
!@var phase indicates from which point this routine was called
      integer, intent(in) :: phase
c
      integer :: ipatch
      type(atmgla_xchng_vars), pointer :: igla
      integer :: i,j,i_0,i_1,j_0,j_1
      real*8, parameter :: H = 6800d0  ! See eq. 3.7: http://paoc.mit.edu/labweb/notes/chap3.pdf
      real*8 :: iglaT_K, atmglaT_K, P, zdiff
      real*8 :: AM1_hPa    ! AM1 in units of hecto-Pascals

      if(ubound(atmglas,1)==1) return ! only one patch. nothing to do

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop


      ! ======= atm_exports_phase1 is composed of:
      ! real*8, dimension(:,:), pointer ::
      !       SRFP  ! SRFP actual surface pressure (hecto-Pascals)
      !      ,SRFPK ! srfp**kapa
      !      ,AM1   ! first-layer air mass (kg/m2)
      !      ,BYAM1 ! 1/AM1
      !      ,P1    ! center pressure of first layer (mb)
      ! !@var PREC precipitation (kg/m^2)
      ! !@var EPREC energy of preciptiation (J/m^2)
      !      ,PREC,EPREC
      ! !@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
      !      ,COSZ1
      ! !@var FLONG, FSHORT downwelling longwave, shortwave radiation at surface
      !      ,FLONG,FSHORT
      !      ,TRUP_in_rad ! LW emission by surface during rad. timestep.
      if(phase==after_atm_phase1) then
        ! quantities that are not updated during surface sub-timesteps
        do ipatch=1+lbound(atmglas,1),ubound(atmglas,1)
          igla => atmglas(ipatch)
          igla%atm_exports_phase1 = atmgla%atm_exports_phase1
#ifdef TRACERS_WATER
          igla%tratm_exports_phase1 = atmgla%tratm_exports_phase1
#endif

          DO J=J_0,J_1
          DO I=I_0,IMAXJ(J)
            ! Downscale Pressure
            ! See eq. 3.7: http://paoc.mit.edu/labweb/notes/chap3.pdf
            zdiff = elevhc(i,j,ipatch) - zatmo(i,j)/Grav
            igla%SRFP(i,j) = atmgla%SRFP(i,j) * exp(-zdiff/H)

            ! Set other things that rely on downscaling
            igla%SRFPK(i,j) = igla%SRFP(i,j) ** kapa
            AM1_hPa = (igla%SRFP(i,j) - ptop) * dsig(1)

            ! Pa = kg / (m s^2)
            ! Grav = m/s^2
            igla%AM1(i,j) = AM1_hPA * (100d0 / Grav)      ! kg/m^2

            ! Center pressure of first layer
            igla%P1(i,j) = igla%SRFP(i,j) - .5*AM1_hPA

            ! Consider setting igla%FLONG.  This depends on the number
            ! of molecules above you radiating down.  The form is
            ! non-trivial.  For now, we are ignoring it.
            ! Consider looking at differences in FLOG between edge and
            ! center of an ice sheet, to see how much this might matter
            ! igla%FLONG(i,j) = ...

          enddo
          enddo
        enddo	! ipatch

      ! ======= atm_exports_phasesrf is composed of:
      !    real*8, dimension(:,:,:), pointer ::
      !         atm_exports_phasesrf=>null()
      !    real*8, dimension(:,:), pointer ::
      !@var TEMP1 pot. temp. of first layer w.r.t. 1 mb (K)
      !@var Q1 specific humidity of first layer
      !@var U1,V1 wind components of first layer (A-grid)
      !       TEMP1
      !      ,Q1
      !      ,U1,V1
      ! The "layer 1" arrays TEMP1,Q1,U1,V1,P1,AM1 will be redefined
      ! to correspond to the full boundary layer depth.
      elseif(phase==during_srfflx) then
        ! quantities that are updated during surface sub-timesteps
        do ipatch=1+lbound(atmglas,1),ubound(atmglas,1)
          igla => atmglas(ipatch)
          igla%atm_exports_phasesrf = atmgla%atm_exports_phasesrf
#ifdef TRACERS_ON
          igla%tratm_exports_phasesrf = atmgla%tratm_exports_phasesrf
#endif

        ! Downscale Temperature
        ! ***** But isn't this supposed to be POTENTIAL Temperature?
        ! Use a simple 8 K/km lapse rate
         DO J=J_0,J_1
         DO I=I_0,IMAXJ(J)
           ! ----------------------------------
           ! Convert from potential temperature to real temperature
           ! (reference pressure = 1mb = 1hPa)
           atmglaT_K = atmgla%TEMP1(i,j) * atmgla%SRFPK(i,j)

           ! Downscale temperature, 8K/km (.008 K/m)
           zdiff = elevhc(i,j,ipatch) - zatmo(i,j)/Grav
           iglaT_K = atmglaT_K - zdiff * HC_T_LAPSE_RATE

           ! Convert back to potential temperature
           igla%TEMP1(i,j) = iglaT_K / igla%SRFPK(i,j)
           ! ----------------------------------

           ! Scale Q1 (kg/kg mixing ratio) assuming constant
           ! humidity as we move up-slope.  That is appropriate,
           ! since we are still just as close to the saturated
           ! ice surface
           igla%Q1(i,j) = atmgla%Q1(i,j) * (
     &          QSAT(iglaT_K,   LHS,   igla%SRFP(i,j)) /
     &          QSAT(atmglaT_K, LHS, atmgla%SRFP(i,j)))

         enddo
         enddo
        enddo		! ipatch
      endif
      return
      end subroutine patchify_landice_inputs
