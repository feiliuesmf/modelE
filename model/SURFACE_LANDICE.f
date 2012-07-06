! Lakes = itype == 1
! Look in GHY_DRV to see which tracers in SURFACE.f are needed for LANDICE surface type

C****   
C**** SURFACE.f    SURFACE fluxes    2006/12/21
C****
#include "rundeck_opts.h"



      SUBROUTINE SURFACE_LANDICE (ns,moddsf,moddd,TGRND,TGRN2,TGR4,e1)
!      integer, intent(in) :: ns,moddsf,moddd
! TGRND...e1 are passed in from SURFACE.f (it's a local variable there)

!@sum SURFACE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : rgas,lhm,lhe,lhs
     *     ,sha,tf,rhow,shv,shi,stbo,bygrav,by6
     *     ,deltx,teeny,grav
      USE MODEL_COM, only : modelEclock
      USE MODEL_COM, only : dtsrc,nday,itime,qcheck
#ifdef SCM
      USE SCMDIAG, only : EVPFLX,SHFLX
      USE SCMCOM, only : iu_scm_prt, ALH, ASH, SCM_SURFACE_FLAG
     &     ,I_TARG,J_TARG
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID, GET
      USE GEOM, only : imaxj
#ifdef TRACERS_ON
      USE TRACER_COM, only : NTM,itime_tr0,needtrs
#ifdef TRACERS_WATER
     *     ,nWATER,tr_wd_TYPE
#endif
#endif
      USE PBL_DRV, only : alloc_pbl_args, dealloc_pbl_args
      USE PBL_DRV, only : pbl, t_pbl_args, xdelt

      USE LANDICE, only : z1e,z2li,hc1li,hc2li,ace1li,ace2li,snmin
      USE LANDICE_COM, only : snowli
#ifdef TRACERS_WATER
     *     ,trlndi
#endif
      USE SEAICE, only : xsi,ace1i,alami0,rhoi,byrls,alami
      USE FLUXES, only : nstype,nisurf,flice,atmgla

      USE Timer_mod, only: Timer_type
      USE TimerList_mod, only: startTimer => start
      USE TimerList_mod, only: stopTimer => stop
      USE itype_enum

      IMPLICIT NONE
! ================ Parameter Declarations =====================
      integer, intent(in) :: ns,moddsf,moddd
      REAL*8, intent(inout), DIMENSION(
     *       NSTYPE,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &       GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *     TGRND,TGRN2,TGR4

    !@var E1 net energy flux at layer 1 (from FLUXES.f)
      REAL*8, intent(inout), DIMENSION(
     &       GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &       GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSTYPE) :: E1

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
     *     ,FSRI(2),dlwdt,byNIsurf,TGO

      REAL*8 MA1
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: Z2LI3L=Z2LI/(3.*ALAMI0), Z1LIBYL=Z1E/ALAMI0
      REAL*8, EXTERNAL :: QSAT
      REAL*8 DQSATDT,TR4
c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 qg_sat,dtsurf,qsrf,us,vs,ws,ws0,
     &     dmua_ij,dmva_ij

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

      real*8, dimension(:,:), pointer :: trdn,srdn

      type (Timer_type), pointer :: aTimer
! ======================= MAIN BODY ======================================

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,        J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C****

      call startTimer('SURFACE_LANDICE()')

      trdn => atmgla%flong
      srdn => atmgla%fshort

      DTSURF=DTsrc/NIsurf
      byNIsurf=1.d0/real(NIsurf)
      IH=modelEclock%hour()+1
      IHM = IH+(modelEclock%date()-1)*24
c avoid uninitialized variable problems when the first gridpoint
c in the domain is ocean
      SNOW = 0.

!! /============================================================\
!! THIS CODE NEEDS to go in SURFACE_LANDICE.f!!!!
!C**** INITIALIZE TGRND: THIS IS USED TO UPDATE T OVER SURFACE STEPS
!! Because we're doing landice, just initialize TGRND(3,...)
!      DO J=J_0,J_1
!      DO I=I_0,I_1
!!        TGRND(2,I,J)=atmice%GTEMP(I,J)
!        TGRND(3,I,J)=atmgla%GTEMP(I,J)
!!        TGRN2(2,I,J)=atmice%GTEMP2(I,J)
!        TGRN2(3,I,J)=atmgla%GTEMP2(I,J)
!!        TGR4(2,I,J)=atmice%GTEMPR(I,J)**4
!        TGR4(3,I,J)=atmgla%GTEMPR(I,J)**4
!      END DO 
!      END DO 


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
      pbl_args % ih = 1+modelEclock%hour()
      pbl_args % ihm = pbl_args%ih+(modelEclock%date()-1)*24
#endif

C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

	  ! Model may expect that some things are zero
      ! call stop_model('Please double-check that model might expect some things to be defined for some boxes, even if there is no ice in that box.', -17)
      PLICE=FLICE(I,J)
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
      ! 	Virtual potential temperature is for computing the static stability of
      ! 	the atmosphere (usually for turbulence parameterizations).  My
      ! 	understanding is that the new code to be introduced for downscaling
      ! 	probably does not need to re-compute any turbulence-related quantities
      ! 	differently than the model is already doing.  So just compute the
      ! 	virtual temperature passed to the PBL using co-located T and Q,
      ! 	getting Q for the height classes with the same interpolation used for
      ! 	T I guess.
      !
      ! 	But keep in mind that virtual temperature effects are negligible in
      ! 	"cold" regions of the globe, and that the parameter xdelt is actually
      ! 	set to zero these days (virtual temp. effects on buoyancy are
      ! 	considered only _within_ the PBL code).
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
      SNOW=SNOWLI(I,J)
          ! Temperature of top ice layer (C)
      TG1=TGRND(3,I,J)

          ! TR4 = TGR4(3,i,j) = atmgla%GTEMPR**4
          ! (Needed for Stefan-Boltzmann Law)
          ! GTEMPR radiative ground temperature over surface type (K)
      TR4=TGR4(3,I,J)

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
        trgrnd2(nx)=TRLNDI(ntix(nx),I,J)/(ACE1LI+ACE2LI)
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
     &         atmgla%TRM1(I,J,n), pbl_args%trs(nx),
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

      ! atmgla%e0 = net energy flux at surface (J/m^2)
      atmgla%E0(I,J)=atmgla%E0(I,J)+F0DT

      ! atmgla%e1 <= e1 (set below in SURFACE.f)
      ! E1 = net energy flux at layer 1
      ! (I believe this means between the top 10cm layer and the 2.9m layer below it)
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT

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
      TGRND(ITYPE,I,J)=TG1  ! includes skin effects
      TGR4(ITYPE,I,J) =TR4
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
!unused      DMUA_IJ=PTYPE*DTSURF*RCDMWS*US
!unused      DMVA_IJ=PTYPE*DTSURF*RCDMWS*VS
!unused      atmgla%DMUA(I,J) = atmgla%DMUA(I,J) + DMUA_IJ
!unused      atmgla%DMVA(I,J) = atmgla%DMVA(I,J) + DMVA_IJ
      atmgla%uflux1(i,j) = RCDMWS*US
      atmgla%vflux1(i,j) = RCDMWS*VS

C****

      END DO  ! end of I loop
      END DO  ! end of J loop
! ============================================================
! Now we're outside the loop over grid points

      call dealloc_pbl_args(pbl_args)

      call stopTimer('SURFACE_LANDICE()')

      RETURN
C****
      END SUBROUTINE SURFACE_LANDICE

