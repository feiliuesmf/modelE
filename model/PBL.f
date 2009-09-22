#include "rundeck_opts.h"


      MODULE SOCPBL
!@sum  SOCPBL deals with boundary layer physics
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from PBLB336E)
!@cont pbl,advanc,stars,getl,dflux,simil,griddr,tfix
!@cont ccoeff0,getk,e_eqn,t_eqn,q_eqn,uv_eqn
!@cont t_eqn_sta,q_eqn_sta,uv_eqn_sta
!@cont inits,tcheck,ucheck,check1,output,rtsafe

      USE CONSTANT, only : grav,pi,radian,bygrav,teeny,deltx,tf
     &     ,by3,lhe,rgas,rhows,mair,byrhows,sha,shv,shw,stbo
      USE SEAICE, only : tfrez
      USE LANDICE, only : snmin
#ifdef TRACERS_ON
#ifdef RUNTIME_NTM
      USE TRACER_COM, only : maxntm
#endif
#ifdef WATER_PROPORTIONAL
      USE TRACER_COM, only : force_limit=>force_limit_pbl
#endif
      USE TRACER_COM, only : ntm,trname,tr_wd_TYPE, nWATER
#ifdef TRACERS_SPECIAL_O18
     *     ,n_water
#endif
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,Ntm_dust
#endif
#ifdef TRACERS_DUST
     &     ,n_clay
#else
#ifdef TRACERS_MINERALS
     &     ,n_clayilli
#else
#ifdef TRACERS_QUARZHEM
     &     ,n_sil1quhe
#endif
#endif
#endif
#endif
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      USE tracers_dust,ONLY : Mtrac
#endif
      USE TRIDIAG_MOD, only :  TRIDIAG
      IMPLICIT NONE

      private

      ! public parameters
      public n,zgs,XCDpbl,kappa,emax,skin_effect,ustar_min

      ! public interfaces
      public advanc,inits,ccoeff0

      ! model coefficients (actually a hack, but leave it for now)
      public rimax,ghmin,ghmax,gmmax0,gm_at_rimax,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,s7,s8,c1,c2,c3,c4,c5,b1,b123,b2,prt
     &     ,g0,d1_3,d2_3,d3_3,d4_3,d5_3
     *     ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *     ,g1,g2,g3,g4,g5,g6,g7,g8

      public t_pbl_args

      integer, parameter :: n=8  !@param n  no of pbl. layers
      integer, parameter :: npbl=n

#ifdef TRACERS_ON
#ifndef RUNTIME_NTM
      integer, parameter :: maxntm=ntm
#endif
#endif

c**** t_pbl_args is a derived type structure which contains all
c**** input/output arguments for PBL
c**** Please, use this structure to pass all your arguments to PBL
c**** Do not use global variables for that purpose !
      type t_pbl_args
        ! input:
        real*8 dtsurf,zs1,tgv,tkv,qg_sat,qg_aver,hemi,tr4
        real*8 evap_max,fr_sat,uocean,vocean,psurf,trhr0
        real*8 tg,elhx,qsol,sss_loc
        logical :: ocean,ddml_eq_1
        ! inout:
        real*8 gusti,tdns,qdns,tprime,qprime,snow
        ! output:
        real*8 us,vs,ws,tsv,qsrf,cm,ch,cq,dskin,ws0
        ! the following args needed for diagnostics
        real*8 psi,dbl,khs,ug,vg,wg,ustar,zgs
#ifdef TRACERS_ON
c**** Attention: tracer arrays in this structure have dim 1:ntm
c**** while local arrays in PBL code have dim 1:ntx
c**** Tracer input/output
!@var trtop,trs tracer mass ratio in level 1/surface
!@var trsfac, trconstflx factors in surface flux boundary cond.
!@var trgrnd2 factor for correction to landice evap
!@var ntx number of tracers that need pbl calculation
!@var ntix index array to map local tracer number to global
!@var trprime anomalous tracer concentration in downdraft
        real*8, dimension(maxntm) :: trtop,trs,trsfac,trconstflx
        real*8, dimension(maxntm) :: trdn1,trprime,trgrnd2
        integer ntx
        integer, dimension(maxntm) :: ntix
#ifdef TRACERS_SPECIAL_O18
        real*8, dimension(maxntm) :: frack
#endif
#ifdef BIOGENIC_EMISSIONS
        real*8 :: emisop 
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
c**** input
!@var pbl_args%wearth earth water of first layer [kg/m^2]
!@var pbl_args%aiearth earth ice of first layer [kg/m^2]
!@var pbl_args%wfcs water field capacity of first ground layer [kg/m^2]
        REAL*8 :: wearth,aiearth,wfcs
!@var pbl_args%ers_data ERS data
!@var pbl_args%src_fnct distribution of preferred sources
        REAL*8 :: ers_data,src_fnct
!@var pbl_args%frclay fraction of clay
!@var pbl_args%frsilt fraction of silt
!@var pbl_args%dryhr number of hours with evaporation-precipitation greater
!@+                  Zero to allow dust emission
!@var pbl_args%vtrsh thresh. wind speed above which dust emis. is allowed [m/s]
        REAL*8 :: frclay,frsilt,dryhr,vtrsh
!@var pbl_args%pprec precipitation at previous time step [kg/m^2]
!@var pbl_args%pevap evaporation at previous time step [kg/m^2]
        REAL*8 :: pprec,pevap
!@var pbl_args%d_dust prescribed daily dust emissions [kg/m^2/s] (e.g. AEROCOM)
        REAL*8 :: d_dust(Ntm_dust)
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
!@var pbl_args%minfr distribution of tracer fractions in grid box
        REAL*8 :: minfr(Mtrac)
#endif
c**** output
!@var pbl_args%pdfint integral of dust emission probability density function
        REAL*8 :: pdfint
!@var pbl_args%wsgcm magnitude of the GCM surface wind - ocean currents [m/s]
!@var pbl_args%wspdf mean surface wind calculated from PDF of wind speed [m/s]
        REAL*8 ::wsgcm,wspdf
!@var pbl_args%wsubtke turbulent kinetic energy velocity scale [m/s]
!@var pbl_args%wsubwd dry convective velocity scale [m/s]
!@var pbl_args%wsubwm moist convective velocity scale [m/s]
!@var pbl_args%wtrsh velocity threshold for dust emission (depends on soil
!@+                  moisture) [m/s]
!@var pbl_args%dust_event1 number of dust events [1]
!@var pbl_args%dust_event2 number of dust events above velocity threshold
!@+                        of cubic emission scheme (diagnositcs only) [1]
!@var pbl_args%dust_flux1 dust flux [kg/m^2/s]
!@var pbl_args%dust_flux2 dust flux from cubic emission scheme (diagnostics
!@+                       only) [kg/m^2/s]
        REAL*8 :: dust_flux(Ntm_dust),dust_flux2(Ntm_dust),wsubtke
     *       ,wsubwd,wsubwm,dust_event1,dust_event2,wtrsh
        REAL*8 :: z(npbl),km(npbl-1),gh(npbl-1),gm(npbl-1),zhat(npbl-1),
     &       lmonin
!@var pbl_args%hbaij accumulated precipitation - evaporation balance  [kg/m^2]
!@var pbl_args%ricntd no. of hours with negative precipitation - evaporation
!@+                   balance [1]
        REAL*8 :: hbaij,ricntd
!@var pbl_args%qdust flag whether conditions for dust emission are fulfilled
        LOGICAL :: qdust
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        real*8 :: DMS_flux,ss1_flux,ss2_flux
#endif
#ifdef TRACERS_DRYDEP
!@var dep_vel turbulent deposition velocity = 1/bulk sfc. res. (m/s)
!@var gs_vel gravitational settling velocity (m/s)
        real*8, dimension(maxntm) :: dep_vel,gs_vel
#endif

#ifdef TRACERS_WATER
!@var tr_evap_max maximum amount of tracer available in ground reservoir
        real*8, dimension(maxntm) :: tr_evap_max
#endif


#ifdef TRACERS_GASEXCH_ocean
        real*8  :: alati,Kw_gas,alpha_gas,beta_gas
#endif
#endif


      end type t_pbl_args


!@dbparam XCDpbl factor for momentum drag coefficient
      real*8 :: XCDpbl=1d0

      real*8 :: rimax,ghmin,ghmax,gmmax0,gm_at_rimax,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,s7,s8,c1,c2,c3,c4,c5,b1,b123,b2,prt

      !for level 3 model only:
      real*8 :: g0,d1_3,d2_3,d3_3,d4_3,d5_3
     *         ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *         ,g1,g2,g3,g4,g5,g6,g7,g8

C**** boundary layer parameters
      real*8, parameter :: kappa=0.40d0 !@var kappa  Von Karman constant
!@var  zgs  height of the surface layer (nominally 10 m)
      real*8, parameter :: zgs=10. !@var zgs height of surface layer (m)

C**** parameters for surface fluxes
      !Hogstrom 1988:
      real*8, parameter :: sigma=0.95d0,sigma1=1.-sigma
      real*8, parameter :: gamamu=19.3d0,gamahu=11.6d0,gamams=6.d0,
     *     gamahs=7.8d0/sigma

      ! Businger 1971:
ccc   real*8, parameter :: sigma=0.74d0,sigma1=1.-sigma
ccc   real*8, parameter :: gamamu=15.0d0,gamahu=9.d0,gamams=4.7d0,
ccc  *     gamahs=4.7d0/sigma


CCC !@var bgrid log-linear gridding parameter
CCC      real*8 :: bgrid

!@var smax,smin,cmax,cmin limits on drag coeffs.
!@var emax limit on turbulent kinetic energy
!@var ustar_min limit on surface friction speed 
      real*8, parameter :: smax=0.25d0,smin=0.005d0,cmax=smax*smax,
     *     cmin=smin*smin,emax=1.d5,ustar_min=5d-2

!@param twoby3 2/3 constant
      real*8, parameter :: twoby3 = 2d0/3d0

!@dbparam skin_effect sets whether skin effects are used or not
      integer :: skin_effect=1  ! Used by default

      CONTAINS

      subroutine advanc(pbl_args,coriol,utop,vtop,qtop,ztop,ts_guess,mdf
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0,ilong,jlat,itype
     &     ,kms,kqs,z0m,z0h,z0q,w2_1,ufluxs,vfluxs,tfluxs,qfluxs
     &     ,u,v,t,q,e
#if defined(TRACERS_ON)
     &     ,tr,ptype,trnradius,trndens,trnmm
#endif
     &     )
!@sum  advanc  time steps the solutions for the boundary layer variables
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
c    input:
!@var  coriol  2.*omega*sin(latitude), the coriolis factor
!@var  utop  x component of wind at the top of the layer
!@var  vtop  y component of wind at the top of the layer
!@var  ttop  virtual potential temperature at the top of the layer
!@var  qtop  moisture at the top of the layer
!@var  tgrnd0 bulk virt. pot. temp. of the ground, at the roughness height
!@var  tg  actual ground temp
!@var  elhx relevant latent heat for qg calc
!@var  qgrnd0 initial moisture at the ground, at the roughness height
!@var  qgrnd_sat saturated moisture at the ground, at the roughness height
!@var  evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
!@var  ztop height of the first model layer, approx 200 m if lm=9
!@var  dtime  time step
!@var  ilong  longitude identifier
!@var  jlat   latitude identifier
!@var  itype  1, ocean; 2, ocean ice; 3, land ice; 4, land
!@var  psurf surface pressure
!@var  trhr0 incident long wave radiation
!@var  sss_loc SSS at i,j
c   output:
!@var  us  x component of the surface wind (i.e., due east)
!@var  vs  y component of the surface wind (i.e., due north)
!@var  tsv  virtual potential surface temperature
!@var  qsrf  surface specific moisture
!@var  kms  surface value of km
!@var  khs  surface value of kh
!@var  kqs  surface value of kq
!@var  ws  magnitude of surface wind modified by buoyancy flux (m/s)
!@var  ustar  friction speed
!@var  cm  dimensionless momentum flux at surface (drag coeff.)
!@var  ch  dimensionless heat flux at surface (stanton number)
!@var  cq  dimensionless moisture flux at surface (dalton number)
!@var  z0m  roughness height for momentum (if itype=1 or 2)
!@var  z0h  roughness height for heat
!@var  z0q  roughness height for moisture
!@var  dskin skin-bulk SST difference
!  more output (may be duplicated)
!@var US     = x component of surface wind, positive eastward (m/s)
!@var VS     = y component of surface wind, positive northward (m/s)
!@var WSGCM  = magnitude of the GCM surface wind - ocean currents (m/s)
!@var WSPDF  = mean surface wind calculated from PDF of wind speed (m/s)
!@var WS     = magn of GCM surf wind - ocean curr + buoyancy + gust (m/s)
!@var TSV    = virtual potential temperature of the surface (K)
!@var QS     = surface value of the specific moisture
!@var DBL    = boundary layer height (m)
!@var KMS    = momentum transport coefficient at ZGS (m**2/s)
!@var KHS    = heat transport coefficient at ZGS (m**2/s)
!@var KHQ    = moist transport coefficient at ZGS (m**2/s)
!@var PPBL   = pressure at DBL (mb)
!@var USTAR  = friction speed (square root of momentum flux) (m/s)
!@var CM     = drag coefficient (dimensionless surface momentum flux)
!@var CH     = Stanton number   (dimensionless surface heat flux)
!@var CQ     = Dalton number    (dimensionless surface moisture flux)
!@var z0m   = roughness length for momentum,
!@+           prescribed for itype=3,4 but computed for itype=1,2 (m)
!@var z0h   = roughness length for temperature (m)
!@var z0q   = roughness length for water vapor (m)
!@var UG     = eastward component of the geostrophic wind (m/s)
!@var VG     = northward component of the geostrophic wind (m/s)
!@var MDF    = downdraft mass flux (m/s)
!@var WINT   = integrated surface wind speed over sgs wind distribution
!@var  u  local due east component of wind
!@var  v  local due north component of wind
!@var  t  local virtual potential temperature
!@var  q  local specific humidity (a passive scalar)
!@var  e  local turbulent kinetic energy
#if defined(TRACERS_ON)
!@var  trtop  tracer conc. at the top of the layer
!@var  trs  surface tracer conc.
!@var  trsfac  factor for pbl_args%trs in surface boundary condition
!@var  trconstflx  constant component of surface tracer flux
!@var  ntx  number of tracers to loop over
!@var  ntix index of tracers used in pbl
#endif
#if defined(TRACERS_ON) && defined(TRACERS_GASEXCH_ocean)
!@var  alati SSS at i,j
!@var  Kw_gas  gas exchange transfer velocity at i,j only over ocean
!@var  alpha_gas  solubility of gas
!@var  beta_gas  conversion term  that includes solubility
#endif
#if defined(TRACERS_ON) && defined(TRACERS_WATER)
!@var  tr_evap_max max amount of possible tracer evaporation
#endif
c  internals:
!@var  n     number of the local, vertical grid points
!@var  lscale turbulence length scale. computed on secondary grid.
!@var  z     altitude of primary vertical grid points
!@var  zhat  altitude of secondary vertical grid points
!@var  dzh   dz evaluated at zhat(i)
!@var  dz    dz evaluated at z(i)
!@var  dxi   (z(n)-z(1))/(n-1)
!@var  km    turbulent momentum tranport coefficient.
!@var  kh    turbulent thermometric conductivity. computed
!@var  ke    transport coefficient for the turbulent kinetic energy.

#ifdef TRACERS_GASEXCH_ocean
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE MODEL_COM, only : nstep=>itime
#endif
#endif

!@var tdns downdraft temperature in K, (i,j)
!@var qdns downdraft humidity in kg/kg, (i,j)

      implicit none

      !-- in/out structure
      type (t_pbl_args), intent(inout) :: pbl_args
      !-- input:
      real*8, intent(in) :: coriol,utop,vtop,qtop,ztop,ts_guess
      real*8, intent(in) :: mdf
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      integer, intent(in) :: ilong,jlat,itype
      !-- output:
      real*8, intent(out) :: kms,kqs,z0h,z0q,w2_1
      real*8, intent(out) :: ufluxs,vfluxs,tfluxs,qfluxs
      !-- inout:
      real*8, intent(inout) :: z0m
      real*8, dimension(n),   intent(inout) :: u,v,t,q
      real*8, dimension(n-1), intent(inout) :: e
#if defined(TRACERS_ON)
!@var  tr local tracer profile (passive scalars)
      real*8, intent(in) :: ptype
      real*8, dimension(ntm), intent(in) :: trnradius,trndens,trnmm
      real*8, dimension(n,ntm), intent(inout) :: tr
#endif

c**** local vars for input from pbl_args
      real*8 :: evap_max,fr_sat,uocean,vocean,psurf,trhr0,tg,elhx,qsol
      real*8 :: dtime,sss_loc,dbl,ug,vg,tgrnd0,ttop,qgrnd_sat,qgrnd0
      real*8 :: tdns,qdns,tprime,qprime,tr4
      logical :: ocean,ddml_eq_1
      real*8 :: gusti
c**** local vars for output to pbl_args
      real*8 :: us,vs,ws,tsv,qsrf,khs,dskin,ustar,cm,ch,cq,wsgcm
      real*8 :: ws0
c**** other local vars
      real*8 :: qsat,deltaSST,tgskin,qnet,ts,rhosrf,qgrnd,tg1
      real*8 :: tstar,qstar,ustar0,test,wstar3,wstar2h,tgrnd,ustar_oc
      real*8 :: bgrid,an2,as2,dudz,dvdz,tau,tgr4skin
      real*8 :: ws02
      real*8, parameter ::  tol=1d-3,w=.5d0
      integer, parameter ::  itmax=50
      integer, parameter :: iprint=0,jprint=41  ! set iprint>0 to debug
      real*8, dimension(n) :: dz,xi,usave,vsave,tsave,qsave
     *       ,usave1,vsave1,tsave1,qsave1
      real*8, dimension(n-1) :: lscale,dzh,xihat,kh,kq,ke,esave,esave1
      integer :: i,iter,ierr  !@var i,iter loop variable
C****
      REAL*8,DIMENSION(n) :: z
      REAL*8,DIMENSION(n-1) :: zhat,km,gm,gh
      REAL*8 :: lmonin,snow
#ifdef TRACERS_ON
      real*8, dimension(n,ntm) :: trsave
      real*8 trcnst,trsf,cqsave,byrho,rh1,evap
      real*8, dimension(n-1) :: kqsave
      integer itr
#ifdef TRACERS_WATER
      real*8 :: trc2         ! could be passed out....
#ifdef TRACERS_SPECIAL_O18
      real*8 :: trc1,trs1    ! could be passed out....
      real*8 :: fac_cq_tr(ntm),fk
#endif
#endif
#ifdef TRACERS_DRYDEP
      real*8 vgs
      real*8 :: tr_dens, tr_radius ! variable tracer density and size
      logical hydrate
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      INTEGER :: n1
      REAL*8 :: dsrcflx,dsrcflx2
      real*8 wspdf,delt
#endif

c**** get input from pbl_args structure
      dtime = pbl_args%dtsurf
      tgrnd0 = pbl_args%tgv
      ttop = pbl_args%tkv
      qgrnd_sat = pbl_args%qg_sat
      qgrnd0 = pbl_args%qg_aver
      evap_max = pbl_args%evap_max
      fr_sat = pbl_args%fr_sat
      uocean = pbl_args%uocean
      vocean = pbl_args%vocean
      psurf = pbl_args%psurf
      trhr0 = pbl_args%trhr0
      tg = pbl_args%tg
      tr4 =  pbl_args%tr4
      elhx = pbl_args%elhx
      qsol = pbl_args%qsol
      sss_loc = pbl_args%sss_loc
      ocean = pbl_args%ocean
      cm = pbl_args%cm
      ch = pbl_args%ch
      cq = pbl_args%cq
      dskin = pbl_args%dskin
      dbl = pbl_args%dbl
      khs = pbl_args%khs
      ug = pbl_args%ug
      vg = pbl_args%vg
      ddml_eq_1 = pbl_args%ddml_eq_1

      call griddr(z,zhat,xi,xihat,dz,dzh,zgs,ztop,bgrid,n,ierr)
      if (ierr.gt.0) then
        print*,"advanc: i,j,itype=",ilong,jlat,itype,u(1),v(1),t(1),q(1)
        call stop_model("PBL error in advanc",255)
      end if

      usave(:)=u(:)
      vsave(:)=v(:)
      tsave(:)=t(:)
      qsave(:)=q(:)
      esave(:)=e(:)

#ifdef TRACERS_ON
      trsave(:,1:pbl_args%ntx)=tr(:,1:pbl_args%ntx)
#endif
      do i=1,n-1
        usave1(i)=usave(i)
        vsave1(i)=vsave(i)
        tsave1(i)=tsave(i)
        qsave1(i)=qsave(i)
        esave1(i)=esave(i)
      end do
      ustar0=0.

      tgrnd=tgrnd0              ! use initial bulk ground temp
      qgrnd=qgrnd0              ! use initial sat humidity
      tgskin=tg                 ! initially assume no skin/bulk difference
      tgr4skin=tr4              ! initially assume no skin/bulk difference
      dskin=0
      ts=t(1)/(1+q(1)*deltx)

      call getl1(e,zhat,dzh,lscale,n)

      if(ddml_eq_1) then
        tdns=pbl_args%tdns
        qdns=pbl_args%qdns
        tprime=tdns-t(1)/(1.+deltx*q(1))
        qprime=qdns-q(1)
      else ! either ddml(ilong,jlat).ne.1 or USE_PBL_E1
        tdns=0.d0
        qdns=0.d0
        tprime=0.d0
        qprime=0.d0
#ifdef TRACERS_ON
        pbl_args%trprime(:)=0
#endif
      endif

      do iter=1,itmax

        if(iter.gt.1) then
          call getl(e,u,v,t,zhat,dzh,lmonin,ustar,lscale,dbl,n)
C**** adjust tgrnd/qgrnd for skin effects over the ocean & lakes
          if (itype.eq.1 .and. skin_effect.gt.0) then
c estimate net flux and ustar_oc from current tg,qg etc.
            ts=t(1)/(1+q(1)*deltx)
            rhosrf=100.*psurf/(rgas*t(1)) ! surface air density
            Qnet= (lhe+tgskin*shv)*cq*rhosrf*(ws*(q(1)-qgrnd)
     &           +(ws-ws0)*qprime)        ! Latent
     &           + sha*ch*rhosrf*(ws*(ts-tgskin)+(ws-ws0)*tprime) ! Sensible
     &           +trhr0-stbo*tgr4skin     ! LW

            ustar_oc=ustar*sqrt(rhosrf*byrhows)
            dskin=deltaSST(Qnet,Qsol,ustar_oc)
            tgskin=0.5*(tgskin+(tg+dskin))   ! smooth changes in iteration
            tgskin=max(tgskin,tf+tfrez(sss_loc))  ! prevent unphysical values
            dskin=tgskin-tg ! net dskin diagnostic
            tgr4skin=(sqrt(sqrt(tr4))+dskin)**4
            qgrnd=qsat(tgskin,elhx,psurf)
            if (ocean) qgrnd=0.98d0*qgrnd  ! use ocean adjustment
            tgrnd=tgskin*(1.+qgrnd*deltx)
          endif
        endif

        call getk(km,kh,kq,ke,gm,gh,u,v,t,e,lscale,dzh,n)
        call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,ts,
     2             u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
#ifdef TRACERS_SPECIAL_O18
     *             fac_cq_tr,
#endif
     3             km,kh,kq,dzh,itype,n)
#ifdef TRACERS_ON
        kqsave=kq
        cqsave=cq
#endif

#ifdef USE_PBL_E1
        call e_eqn(esave,e,u,v,t,km,kh,ke,lscale,dz,dzh,
     2                 ustar,dtime,n)
#endif

        !@var wstar the convection-induced wind according to
        !@+ M.J.Miller et al. 1992, J. Climate, 5(5), 418-434, Eqs(6-7),
        !@+ for heat and mositure

        if(t(2).lt.t(1)) then !convective
          wstar3=-dbl*grav*2.*(t(2)-t(1))*kh(1)/((t(2)+t(1))*dzh(1))
          wstar2h = wstar3**twoby3
#ifdef USE_PBL_E1
          gusti=0.
#else
          ! Redelsperger et al. 2000, eqn(13), J. Climate, 13, 402-421
          gusti=pbl_args%gusti ! defined in PBL_DRV.f
#endif
        else
          wstar3=0.
          wstar2h=0.
          gusti=0.
        endif

#ifndef USE_PBL_E1
        call e_eqn(esave,e,u,v,t,km,kh,ke,lscale,dz,dzh,
     2                 ustar,dtime,n)

        call e_les(tstar,ustar,wstar3,dbl,lmonin,zhat,lscale,e,n)
#endif

        ! Inclusion of gustiness in surface fluxes
        ! Redelsperger et al. 2000, eqn(13), J. Climate, 13, 402-421

        ws02=(u(1)-uocean)**2+(v(1)-vocean)**2+wstar2h
        ws0=sqrt(ws02)
        ws=sqrt(ws02+gusti*gusti)

        call q_eqn(qsave,q,kq,dz,dzh,cq,ws,qgrnd_sat,qtop,dtime,n
     &       ,evap_max,fr_sat,ws0,qprime,qdns,ddml_eq_1)

        call t_eqn(u,v,tsave,t,q,z,kh,kq,dz,dzh,ch,ws,tgrnd,ttop,dtime
     *       ,n,dpdxr,dpdyr,dpdxr0,dpdyr0,ws0,tprime,tdns,ddml_eq_1)

        call uv_eqn(usave,vsave,u,v,z,km,dz,dzh,ustar,cm,z0m,utop,vtop
     *       ,dtime,coriol,ug,vg,uocean,vocean,n,dpdxr,dpdyr,dpdxr0
     *       ,dpdyr0)

        if ((ttop.gt.tgrnd).and.(lmonin.lt.0.)) call tfix(t,z,ttop,tgrnd
     *       ,lmonin,tstar,ustar,kh(1),ts_guess,n)

        if(ddml_eq_1) then
          tprime=tdns-t(1)/(1.+deltx*q(1))
          qprime=qdns-q(1)
        endif

        test=abs(2.*(ustar-ustar0)/(ustar+ustar0))
        if (test.lt.tol) exit

        if (iter.lt.itmax) then
        do i=1,n-1
          u(i)=w*usave1(i)+(1.-w)*u(i)
          v(i)=w*vsave1(i)+(1.-w)*v(i)
          t(i)=w*tsave1(i)+(1.-w)*t(i)
          q(i)=w*qsave1(i)+(1.-w)*q(i)
          e(i)=w*esave1(i)+(1.-w)*e(i)
          usave1(i)=u(i)
          vsave1(i)=v(i)
          tsave1(i)=t(i)
          qsave1(i)=q(i)
          esave1(i)=e(i)
        end do
        end if
        ustar0=ustar
cc      WRITE(0,*) 'advanc iter t(1),q(1) ',iter,t(1),q(1)

      end do ! end of iter loop

c**** cannot update ws without taking care that ws used for tracers is
c**** the same as that used for q
      wsgcm=sqrt((u(1)-uocean)**2+(v(1)-vocean)**2)

C**** Preliminary coding for use of sub-gridscale wind distribution
C**** generic calculations for all tracers
C**** To use, uncomment next two lines and adapt the next chunk for
C**** your tracers. The integrated wind value is passed back to SURFACE
C**** and GHY_DRV. This may need to be tracer dependent?
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      delt = t(1)/(1.+q(1)*deltx) - tgrnd/(1.+qgrnd*deltx)
      CALL sig(e(1),mdf,dbl,delt,ch,wsgcm,t(1),pbl_args%wsubtke
     &     ,pbl_args%wsubwd,pbl_args%wsubwm)
      CALL get_wspdf(pbl_args%wsubtke,pbl_args%wsubwd,pbl_args%wsubwm
     &     ,wsgcm,wspdf)
#endif
csgs      sig0 = sig(e(1),mdf,dbl,delt,ch,ws,t(1))
csgsC**** possibly tracer specific coding
csgs      wt = 3.                 ! threshold velocity
csgs      wmin = wt               ! minimum wind velocity (usually wt?)
csgs      wmax = 50.              ! maxmimum wind velocity
csgs      icase=3                 ! icase=3 ==> w^3 dependency
csgs      call integrate_sgswind(sig0,wt,wmin,wmax,ws,icase,wint)

#ifdef TRACERS_ON
C**** tracer calculations are passive and therefore do not need to
C**** be inside the iteration. Use moisture diffusivity.
C**** First, define some useful quantities
      ts=t(1)/(1.+q(1)*deltx)   ! surface air temp (K)
      rhosrf=100.*psurf/(rgas*t(1)) ! surface air density
      byrho=1d0/rhosrf
      tg1 = tgskin-tf ! re-calculate ground T (C)
      rh1=q(1)/qsat(ts,lhe,psurf) ! rel. hum. at surface (wrt water)
      evap=cqsave*rhosrf*(ws*(qgrnd-q(1))-(ws-ws0)*qprime)  ! net evap

#ifdef TRACERS_DRYDEP
C**** Get tracer deposition velocity (= 1 / bulk sfc resistance)
C**** for all dry deposited tracers
      call get_dep_vel(ilong,jlat,itype,lmonin,dbl,ustar,ts
     &     ,pbl_args%dep_vel,trnmm)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      CALL dust_emission_constraints(itype,ptype,wsgcm,pbl_args)
#endif

#ifdef WATER_PROPORTIONAL
      force_limit = .false.
#endif
C**** loop over tracers
      do itr=1,pbl_args%ntx
C**** Define boundary conditions

C****   1) default air mass tracers
        trcnst=pbl_args%trconstflx(itr)*byrho   ! convert to (conc * m/s)
        trsf=pbl_args%trsfac(itr)

#ifdef TRACERS_WATER
C****   2) water mass tracers
C**** Water tracers need to multiply trsfac/trconstflx by cq*Usurf
C**** and qgrnd_sat (moved from driver routines to deal with skin effects)
        if (tr_wd_TYPE(pbl_args%ntix(itr)).eq.nWATER) then
          trcnst=cqsave*(pbl_args%trconstflx(itr)*ws*qgrnd_sat-
     *         (ws-ws0)*pbl_args%trdn1(itr))
          trc2=pbl_args%trconstflx(itr)
          if (ddml_eq_1) then   ! hmmm... gusti>0 even if ddml_eq_1=F
            trsf=pbl_args%trsfac(itr)*cqsave*ws0
          else
            trsf=pbl_args%trsfac(itr)*cqsave*ws
          end if
#ifdef TRACERS_SPECIAL_O18
          if (itype.eq.3) then  ! possible correction for large E over LI
            if (evap.gt.pbl_args%snow .and. pbl_args%snow.gt.snmin) then
              trc2 = (pbl_args%snow*pbl_args%trconstflx(itr)+(evap ! weighted mean tracer conc
     *             -pbl_args%snow)*pbl_args%trgrnd2(itr))/evap 
              trcnst=cqsave*(trc2*ws*qgrnd_sat-
     *             (ws-ws0)*pbl_args%trdn1(itr))
            end if
          end if
C**** get fractionation for isotopes
          call get_frac(itype,ws,tg1,trcnst,trsf,evap*byrho
     $         ,trc2,q(1),pbl_args%ntix(itr)
     $         ,fac_cq_tr(itr),trc1,trs1,fk)
          pbl_args%frack(itr)=fk  ! save for output
          trcnst=trc1
          trsf  =trs1
#endif
        end if
#endif

#ifdef TRACERS_DRYDEP
C****   3) dry deposited tracers (including gravitational settling)
C**** Tracer Dry Deposition boundary condition for dry dep tracers:
        if(dodrydep(pbl_args%ntix(itr))) then
C****   get settling velocity

c****   set tracer size and density (not constant anymore)
        tr_radius = trnradius(pbl_args%ntix(itr))
        tr_dens =   trndens(pbl_args%ntix(itr))

C**** need to hydrate the sea salt before determining settling
          hydrate = (trname(pbl_args%ntix(itr)).eq.'seasalt1'.or.
     *         trname(pbl_args%ntix(itr)).eq.'seasalt2'.or.
     *         trname(pbl_args%ntix(itr)).eq.'M_SSA_SS'.or.
     *         trname(pbl_args%ntix(itr)).eq.'M_SSC_SS'.or.
     *         trname(pbl_args%ntix(itr)).eq.'M_SSS_SS')

          if (trnradius(pbl_args%ntix(itr)).gt.0.) then
            pbl_args%gs_vel(pbl_args%ntix(itr))=vgs(rhosrf,rh1
     &           ,tr_radius,tr_dens,ts,hydrate)
          else
            pbl_args%gs_vel(pbl_args%ntix(itr))=0.
          end if
          trsf=pbl_args%trsfac(itr)*(pbl_args%dep_vel(pbl_args%ntix(itr)
     &         )+pbl_args%gs_vel(pbl_args%ntix(itr)))
        end if
#endif

C****   4) tracers with interactive sources
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        select case (trname(pbl_args%ntix(itr)))
        case ('DMS')
          call read_DMS_sources(ws,itype,ilong,jlat,pbl_args%DMS_flux)
          trcnst=pbl_args%DMS_flux*byrho
        case ('seasalt1', 'M_SSA_SS')
          call read_seasalt_sources(ws,itype,1,ilong,jlat
     &         ,pbl_args%ss1_flux)
          trcnst=pbl_args%ss1_flux*byrho
        case ('seasalt2', 'M_SSC_SS')
          call read_seasalt_sources(ws,itype,2,ilong,jlat
     &         ,pbl_args%ss2_flux)
          trcnst=pbl_args%ss2_flux *byrho
        case ('M_SSS_SS')
          call read_seasalt_sources(ws,itype,1,ilong,jlat
     &         ,pbl_args%ss1_flux)
          call read_seasalt_sources(ws,itype,2,ilong,jlat
     &         ,pbl_args%ss2_flux)
          trcnst=(pbl_args%ss2_flux + pbl_args%ss1_flux)*byrho
        end select
#endif

#ifdef BIOGENIC_EMISSIONS
! Nadine Unger test code for biogenic emissions
! emisop is in units kg C/m2/s
        select case (trname(pbl_args%ntix(itr)))
        case ('Isoprene')
          call isoprene_emission(ilong,jlat,itype,pbl_args%emisop)
          trcnst=pbl_args%emisop*byrho
        end select
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
ccc dust emission from earth
        SELECT CASE (trname(pbl_args%ntix(itr)))
#ifdef TRACERS_DUST
        CASE ('Clay','Silt1','Silt2','Silt3','Silt4')
          n1=pbl_args%ntix(itr)-n_clay+1
#else
#ifdef TRACERS_MINERALS
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &        'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &        'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &        'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &        'Sil1QuHe','Sil2QuHe','Sil3QuHe')
          n1=pbl_args%ntix(itr)-n_clayilli+1
#else
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          n1=pbl_args%ntix(itr)-n_sil1quhe+1
#endif
#endif
#endif
        END SELECT
        SELECT CASE (trname(pbl_args%ntix(itr)))
        CASE ('Clay','Silt1','Silt2','Silt3','Silt4',
     &        'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &        'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &        'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &        'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &        'Sil1QuHe','Sil2QuHe','Sil3QuHe')
          CALL local_dust_emission(pbl_args%ntix(itr),ptype,wsgcm,
     &       pbl_args,dsrcflx,dsrcflx2)
          trcnst=dsrcflx*byrho
          pbl_args%dust_flux(n1)=dsrcflx
          pbl_args%dust_flux2(n1)=dsrcflx2
        END SELECT
#endif

#ifdef TRACERS_AMP
        SELECT CASE (trname(pbl_args%ntix(itr)))
        CASE ('M_DD1_DU')
          DO n1=1,2
            CALL local_dust_emission(n1,ptype,wsgcm,pbl_args,dsrcflx,
     &           dsrcflx2)
            trcnst=dsrcflx*byrho
            pbl_args%dust_flux(n1)=dsrcflx
            pbl_args%dust_flux2(n1)=dsrcflx2
          END DO
        CASE ('M_DD2_DU')
          DO n1=3,4
            CALL local_dust_emission(n1,ptype,wsgcm,pbl_args,dsrcflx,
     &           dsrcflx2)
            trcnst=dsrcflx*byrho
            pbl_args%dust_flux(n1)=dsrcflx
            pbl_args%dust_flux2(n1)=dsrcflx2
          END DO
        CASE ('M_DDD_DU')
          DO n1=1,4
            CALL local_dust_emission(n1,ptype,wsgcm,pbl_args,dsrcflx,
     &           dsrcflx2)
            trcnst=dsrcflx*byrho
            pbl_args%dust_flux(n1)=dsrcflx
            pbl_args%dust_flux2(n1)=dsrcflx2
          END DO
        END SELECT
#endif

#ifdef TRACERS_GASEXCH_ocean
#ifdef TRACERS_GASEXCH_ocean_CO2
       IF (ocean) THEN  ! OCEAN only
          !note trcnst is already multiplied by byrho in TRACERS_GASEXCH_ocean_CO2_PBL

          call TRACERS_GASEXCH_ocean_CO2_PBL(tg1,ws,
     .          pbl_args%alati,psurf,itr,pbl_args%trconstflx(itr),
     .          byrho,pbl_args%Kw_gas,pbl_args%alpha_gas,
     .          pbl_args%beta_gas,trsf,trcnst,ilong,jlat)

!     write(*,'(a,2i5,4e12.4,i5,7e12.4)')'PBL:', 
!    .  ilong,jlat,tg1,ws,pbl_args%alati,psurf,itr,
!    .             pbl_args%trconstflx(itr),
!    .          byrho,pbl_args%Kw_gas,pbl_args%alpha_gas,
!    .          pbl_args%beta_gas,trsf,trcnst
       ELSE
        trsf = 0.d0
        trcnst = 0.d0
       ENDIF
#endif   /* TRACERS_GASEXCH_ocean_CO2 */

#ifdef TRACERS_GASEXCH_ocean_CFC
      IF (ocean) THEN  ! OCEAN only
       call TRACERS_GASEXCH_ocean_CFC_PBL(tg1,ws,
     .          pbl_args%alati,psurf,itr,pbl_args%trconstflx(itr),
     .          byrho,pbl_args%Kw_gas,pbl_args%alpha_gas,
     .          pbl_args%beta_gas,trsf,trcnst,ilong,jlat)
      ENDIF
#endif /* TRACERS_GASEXCH_ocean_CFC */
#endif /* TRACERS_GASEXCH_ocean */

C**** solve tracer transport equation
        call tr_eqn(trsave(1,itr),tr(1,itr),kqsave,dz,dzh,trsf
     *       ,trcnst,pbl_args%trtop(itr),
#ifdef TRACERS_WATER
     *       pbl_args%tr_evap_max(itr),fr_sat,
#endif
     *       dtime,n)

#ifdef TRACERS_DRYDEP
C**** put in a check to prevent unphysical solutions. If too much
C**** tracer is being taken out, replace profile with linear one
C**** with maximum allowed flux.
        if (dodrydep(pbl_args%ntix(itr))) then
          if ((trsf*tr(1,itr)-trcnst)*dtime
     &         .gt.pbl_args%trtop(itr)*ztop) then
            do i=1,n
              tr(i,itr)=(pbl_args%trtop(itr)*ztop/dtime+trcnst)/trsf
     &             +(i-1) *pbl_args%trtop(itr)/float(n-1)
            end do
          end if
        end if
#endif
      end do
#endif

      us  = u(1)
      vs  = v(1)
      tsv = t(1)
      qsrf  = q(1)

      kms = km(1)
      khs = kh(1)
      kqs = kq(1)

      ufluxs=km(1)*(u(2)-u(1))/dzh(1)
      vfluxs=km(1)*(v(2)-v(1))/dzh(1)
      tfluxs=kh(1)*(t(2)-t(1))/dzh(1)
      qfluxs=kq(1)*(q(2)-q(1))/dzh(1)

      an2=2.*grav*(t(n)-t(n-1))/((t(n)+t(n-1))*dzh(n-1))
      dudz=(u(n)-u(n-1))/dzh(n-1)
      dvdz=(v(n)-v(n-1))/dzh(n-1)
      as2=dudz*dudz+dvdz*dvdz
      tau=B1*lscale(n-1)/max(sqrt(2.*e(n-1)),teeny)
      !@var w2_1 the vertical component of 2*e at GCM layer 1
      w2_1=twoby3*e(n-1)-tau*by3*(s7*km(n-1)*as2+s8*kh(n-1)*an2)
      w2_1=max(0.24d0*e(n-1),w2_1) ! Mellor-Yamada 1982

c     call check1(ustar,1,ilong,jlat,2)

c Diagnostics printed at a selected point:

      if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,kq,ke,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif

c**** copy output to pbl_args
      pbl_args%us = us
      pbl_args%vs = vs
      pbl_args%ws = ws    ! wind speed including buoyancy + gusti
      pbl_args%ws0 = ws0  ! wind speed including buoyancy
      pbl_args%tsv = tsv
      pbl_args%qsrf = qsrf
      pbl_args%cm = cm
      pbl_args%ch = ch
      pbl_args%cq = cq
      pbl_args%dskin = dskin
      !!pbl_args%psi = psi   ! maybe compute it here ?
      pbl_args%dbl = dbl
      pbl_args%khs = khs
      pbl_args%ustar = ustar
      pbl_args%zgs = zgs
      pbl_args%gusti = gusti
      pbl_args%tprime=tprime
      pbl_args%qprime=qprime

C**** tracer code output
#ifdef TRACERS_ON
      pbl_args%trs(1:ntm) = tr(1,1:ntm)

      if (ddml_eq_1)
     &     pbl_args%trprime(1:ntm) = pbl_args%trdn1(1:ntm)-tr(1,1:ntm)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      pbl_args%wsgcm=wsgcm
      pbl_args%wspdf=wspdf
      pbl_args%z(:) = z(:)
      pbl_args%zhat(:) = zhat(:)
      pbl_args%km(:) = km(:)
      pbl_args%gm(:) = gm(:)
      pbl_args%gh(:) = gh(:)
      pbl_args%lmonin = lmonin
#endif
#endif

      return
      end subroutine advanc

      subroutine stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,ts,
     2                 u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
#ifdef TRACERS_SPECIAL_O18
     *                 fac_cq_tr,
#endif
     3                 km,kh,kq,dzh,itype,n)
!@sum computes USTAR,TSTAR and QSTAR
!@+   Momentum flux = USTAR*USTAR
!@+   Heat flux     = USTAR*TSTAR
!@+   MOISTURE flux = USTAR*QSTAR
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from PBLB336E)
!@var USTAR the friction speed
!@var TSTAR the virtual potential temperature scale
!@var QSTAR the moisture scale
!@var LMONIN the Monin-Obukhov length scale
      USE CONSTANT, only : teeny
      implicit none

      integer, intent(in) :: itype,n
      real*8, dimension(n), intent(in) :: u,v,t,q,z
      real*8, dimension(n-1), intent(in) :: km,kh,kq,dzh
      real*8, intent(in) :: tgrnd,qgrnd,ts
      real*8, intent(inout) :: z0m
      real*8, intent(out) :: ustar,tstar,qstar,lmonin
      real*8, intent(out) :: z0h,z0q,cm,ch,cq
#ifdef TRACERS_SPECIAL_O18
      real*8, intent(out) :: fac_cq_tr(ntm)
#endif


      real*8 dz,vel1,du1,dv1,dudz,dtdz,dqdz,zgs

      dz     = dzh(1)
      vel1   = sqrt(u(1)*u(1)+v(1)*v(1))
      du1=u(2)-u(1)
      dv1=v(2)-v(1)
      dudz=sqrt(du1*du1+dv1*dv1)/dz
      dtdz   = (t(2)-t(1))/dz
      dqdz   = (q(2)-q(1))/dz
      ustar  = sqrt(km(1)*dudz)
#ifdef USE_PBL_E1
      ustar  = max(ustar,teeny)
#else
      ustar  = max(ustar,ustar_min)
#endif
      tstar  = kh(1)*dtdz/ustar
      qstar  = kq(1)*dqdz/ustar
      zgs    = z(1)

      if (ustar.gt.smax*vel1) ustar=smax*vel1
      if (ustar.lt.smin*vel1) ustar=smin*vel1
      if (abs(tstar).gt.smax*abs(t(1)-tgrnd)) tstar=smax*(t(1)-tgrnd)
      if (abs(tstar).lt.smin*abs(t(1)-tgrnd)) tstar=smin*(t(1)-tgrnd)
#ifdef USE_PBL_E1
      ! do nothing
#else
      if (tstar.eq.0.) tstar=teeny
      ustar  = max(ustar,ustar_min)
#endif
      if (abs(qstar).gt.smax*abs(q(1)-qgrnd)) qstar=smax*(q(1)-qgrnd)
      if (abs(qstar).lt.smin*abs(q(1)-qgrnd)) qstar=smin*(q(1)-qgrnd)

      lmonin = ustar*ustar*tgrnd/(kappa*grav*tstar)
      if(abs(lmonin).lt.teeny) lmonin=sign(teeny,lmonin)
c**** To compute the drag coefficient,Stanton number and Dalton number
      call dflux(lmonin,ustar,vel1,ts,z0m,z0h,z0q,zgs,cm,ch,cq,
#ifdef TRACERS_SPECIAL_O18
     *     fac_cq_tr,
#endif
     *     itype)

      return
      end subroutine stars

      subroutine getl1(e,zhat,dzh,lscale,n)
!@sum getl1 estimates the master length scale of the turbulence model
!@+   on the secondary grid
!@auth  Ye Cheng/G. Hartke
      implicit none

      integer, intent(in) :: n     !@var n  array dimension
      real*8, dimension(n-1), intent(in) :: e,zhat,dzh
      real*8, dimension(n-1), intent(out) :: lscale
      real*8, parameter :: alpha=0.2d0

      real*8 :: sum1,sum2,l0,l1
      integer :: j

      sum1=0.
      sum2=0.
      do j=1,n-1
        sum1=sum1+sqrt(e(j))*zhat(j)*dzh(j)
        sum2=sum2+sqrt(e(j))*dzh(j)
      end do
      l0=alpha*sum1/sum2

      do j=1,n-1
        l1=kappa*zhat(j)
        lscale(j)=l0*l1/(l0+l1)
      end do

      return
      end subroutine getl1

      subroutine getl(e,u,v,t,zhat,dzh,lmonin,ustar,lscale,dbl,n)
!@sum   getl computes the master length scale of the turbulence model
!@+     on the secondary grid. l0 in this routine is 0.16*(pbl height)
!@+     according to the LES data (Moeng and Sullivan 1992)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var e z-profle of turbulent kinetic energy
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var lscale z-profile of the turbulent dissipation length scale
!@var z vertical grids (main, meter)
!@var zhat vertical grids (secondary, meter)
!@var dzh(j)  z(j+1)-z(j)
!@var dbl PBL height (m)
!@var n number of vertical subgrid main layers

      USE CONSTANT, only : by3
      implicit none

      integer, intent(in) :: n   !@var n  array dimension
      real*8, dimension(n-1), intent(in) :: e,zhat,dzh
      real*8, dimension(n), intent(in) :: u,v,t
      real*8, dimension(n-1), intent(out) :: lscale
      real*8, intent(in) :: dbl,lmonin,ustar

      integer :: i   !@var i  array dimension
      real*8 kz,l0,ls,lb,an2,an,qty,qturb,zeta

#ifdef USE_PBL_E1
      l0=.16d0*dbl ! Moeng and Sullivan 1994
#else
      l0=.3d0*dbl ! Moeng and Sullivan 1994
#endif

      if (l0.lt.zhat(1)) l0=zhat(1)

      kz=kappa*zhat(1)
      lscale(1)=l0*kz/(l0+kz)

      do i=2,n-1
        kz=kappa*zhat(i)
        zeta=zhat(i)/lmonin
        ! Nakanishi (2001)
        if(zeta.ge.1.) then
          ls=kz/3.7d0
        elseif(zeta.ge.0.) then
          ls=kz/(1.+2.7d0*zeta)
        else
          ls=kz*(1.-100.*zeta)**0.2d0
        endif
        if (t(i+1).gt.t(i)) then
          an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
          an=sqrt(an2)
          qturb=sqrt(2*e(i))
          if(zeta.ge.0.) then
             lb=qturb/an
          else
             qty=(ustar/((-kappa*lmonin*l0*l0)**by3*an))**0.5d0
             lb=qturb*(1.+5.*qty)/an
          endif
        else
          lb=1.d30
        endif
        lscale(i)=l0*ls*lb/(l0*ls+l0*lb+ls*lb)
      end do

      return
      end subroutine getl

      subroutine dflux(lmonin,ustar0,vsurf,ts,z0m,z0h,z0q,zgs,
     *                 cm,ch,cq,
#ifdef TRACERS_SPECIAL_O18
     *                 fac_cq_tr,
#endif
     *                 itype)
!@sum   dflux computes (dimensionless) surface fluxes of momemtun,
!@+     heat and moisture (drag coefficient, Stanton number,
!@+     and Dalton number)
!@+     Now with explicit Sc and Pr number dependence
!@+     and flexibility for water isotopes
!@auth  Ye Cheng/G. Hartke (mods by G. Schmidt)
!@ver   1.0
!@var lmonin = Monin-Obukhov length (m)
!@var ustar  = friction speed (sqrt of surface momentum flux) (m/sec)
!@var vsurf  = total surface wind speed, used to limit ustar -> cm
!@var ts     = surface air temperature (K)
!@var zgs    = height of the surface layer (m)
!@var itype  = integer identifying surface type
!@var z0m   = momentum roughness length, prescribed (itype=3,4) (m)
!@var z0m   = roughness length for momentum, computed (itype=1,2)
!@var cm    = drag coefficient for momentum
!@var ch    = Stanton number
!@var cq    = Dalton number
!@var z0h   = roughness length for temperature (m)
!@var z0q   = roughness length for water vapor (m)
!@var Sc    = Schmidt (no relation) number (visc_air_kin/diff)
!@var Pr    = Prandtl number (visc_air_kin/therm_diff)
!@var fac_cq_tr = ratio of cq for water isotopes = f(Sc_tr)
      USE CONSTANT, only : visc_air_kin
      implicit none

      real*8,  intent(in) :: lmonin,ustar0,vsurf,zgs,ts
      integer,  intent(in) :: itype
      real*8,  intent(inout) :: z0m
      real*8,  intent(out) :: cm,ch,cq,z0h,z0q
#ifdef TRACERS_SPECIAL_O18
      real*8, intent(out) :: fac_cq_tr(ntm)
      real*8 :: cq_tr(ntm),z0q_tr(ntm),Sc_tr,get_diff_rel
      integer :: itr
#endif

      real*8 :: nu,num,nuh,nuq
      real*8 dm,ustar,dum
      real*8, parameter :: Sc=0.595d0, Pr=0.71d0

C**** Kinematic viscosity
#ifdef PBL_E1
      nu=1.5d-5                 ! temperature independent (original E1 value)
#else
      nu=visc_air_kin(ts)       ! temperature dependent
#endif

      if ((itype.eq.1).or.(itype.eq.2)) then
c *********************************************************************
        ustar = max(ustar0,1.0125d-5)  ! make sure not too small
c Compute roughness lengths using smooth/rough surface formulation:

c       z0m=0.11d0*nu/ustar+0.011d0*ustar*ustar*bygrav ! COARE algorithm
        z0m=0.135d0*nu/ustar+0.018d0*ustar*ustar*bygrav ! Hartke and Rind (1996)

#ifdef USE_PBL_E1
        nuh=0.395d0*nu ;  nuq=0.624d0*nu
c**** more accurate nuh, nuq assuming Sc and Pr so that this is
c**** consistent with formula in getzhq. Note '0.624' is for Sc=0.6.
c       nuh=0.39522377113362589d0*nu ;  nuq=0.63943320118296587d0*nu

        z0h=nuh/ustar + 1.4d-5
        z0q=nuq/ustar + 1.3d-4
#else
        call getzhq(ustar,z0m,Pr,nu,1.4d-5,z0h)  ! heat
        call getzhq(ustar,z0m,Sc,nu,1.3d-4,z0q)  ! vapour
#endif

#ifdef TRACERS_SPECIAL_O18
C**** calculate different z0q for different diffusivities
          do itr=1,ntm
            if (tr_wd_TYPE(itr).eq.nWater) then
              Sc_tr=Sc*get_diff_rel(itr)
              call getzhq(ustar,z0m,Sc_tr,nu,1.3d-4,z0q_tr(itr))
            end if
          end do
#endif

      else
c *********************************************************************
c  For land and land ice, z0m is specified. For z0h and z0q,
c    empirical evidence suggests:
        z0h=z0m*.13533528d0    ! = exp(-2.)
        z0q=z0h
#ifdef TRACERS_SPECIAL_O18
        z0q_tr(:)=z0q
#endif
c *********************************************************************
      endif

      call getcm(zgs,z0m,lmonin,dm,cm)
      call getchq(zgs,z0m,lmonin,dm,z0h,dum,ch)
      call getchq(zgs,z0m,lmonin,dm,z0q,dum,cq)

#ifdef TRACERS_SPECIAL_O18
      do itr=1,ntm
        if (tr_wd_TYPE(itr).eq.nWater)
     *       call getchq(zgs,z0m,lmonin,dm,z0q_tr(itr),cq_tr(itr),dum)
      end do
      do itr=1,ntm
        if (tr_wd_TYPE(itr).eq.nWater)
     *       fac_cq_tr(itr)=cq_tr(itr)/cq_tr(n_Water)
      end do
#endif

      return
      end subroutine dflux

      subroutine getzhq(ustar,z0m,ScPr,nu,z0min,z0hq)
!@sum calculate z0hq heat/humidity roughness length
!@+   modified from eqs 5.24, 5.27 and 5.35 in Brutsaert (1982)
!!*** remove z0min for original HR97 code
      implicit none
      real*8, intent(in) :: ustar,z0m,ScPr,z0min,nu
      real*8, intent(out) :: z0hq
      real*8 r0q,beta,fac_smooth_ScPr,fac_rough_ScPr

C**** functional dependence on Sc,Pr for smooth, rough surfaces
      fac_smooth_ScPr = 30.*exp(-13.6d0*kappa*ScPr**twoby3)

      if (ustar.le.0.20d0) then ! smooth regime
        z0hq=nu*fac_smooth_ScPr/ustar + z0min
      else                      ! rough regime
        fac_rough_ScPr = -7.3d0*kappa*sqrt(ScPr)
        r0q=sqrt(sqrt(ustar*z0m/nu))
        z0hq=7.4d0*z0m*exp(fac_rough_ScPr*r0q)
        if (ustar.lt.0.2d0) then ! intermediate regime (lin. interp.)
          beta=(ustar-0.02d0)/0.18d0
          z0hq=(1.-beta)*(nu*fac_smooth_ScPr/ustar+z0min)+beta*z0hq
        endif
      endif

      return
      end subroutine getzhq

      subroutine getcm(zgs,z0m,lmonin,dm,cm)
!@sum calculate cm drag coefficient for momentum
!@+   Hartke and Rind (1997)
      implicit none
      real*8, intent(in) :: zgs,z0m,lmonin
      real*8, intent(out) :: dm,cm

      real*8 zgsbyl,z0mbyl,cmn,dpsim,xms,xm0,lzgsbyz0m

      lzgsbyz0m = log(zgs/z0m)
      cmn=kappa*kappa/(lzgsbyz0m**2)

      zgsbyl=zgs/lmonin
      z0mbyl=z0m/lmonin

c *********************************************************************
c  Now compute DPSI, which is the difference in the psi functions
c    computed at zgs and the relevant roughness height:
c *********************************************************************

      if (lmonin.gt.0.) then
c *********************************************************************
c  Here the atmosphere is stable with respect to the ground:
        dpsim=-gamams*(zgsbyl-z0mbyl)
c *********************************************************************
      else
c *********************************************************************
c  Here the atmosphere is unstable with respect to the ground:
        xms  =    (1.-gamamu*zgsbyl)**0.25d0
        xm0  =    (1.-gamamu*z0mbyl)**0.25d0
        dpsim=log((1.+xms)*(1.+xms)*(1.+xms*xms)/
     2           ((1.+xm0)*(1.+xm0)*(1.+xm0*xm0)))-
     3        2.*(atan(xms)-atan(xm0))
c *********************************************************************
      endif

      dm=      1./(1.-min(dpsim/lzgsbyz0m,.9d0))**2
      cm=XCDpbl*dm*cmn
      if (cm.gt.cmax) cm=cmax
      if (cm.lt.cmin) cm=cmin

      return
      end subroutine getcm

      subroutine getchq(zgs,z0m,lmonin,dm,z0hq,chq0,chq)
!@sum calculate chq drag coefficients for heat/water
!@+   Hartke and Rind (1997)
      implicit none
      real*8, intent(in) :: zgs,z0m,lmonin,dm,z0hq
      real*8, intent(out) :: chq,chq0  ! final and unlimited version

      real*8 zgsbyl,z0hqbyl,chqn,dpsihq,xhqs,xhq0,lzgsbyz0m,lzgsbyz0hq
     *     ,dhq

      lzgsbyz0m  = log(zgs/z0m)
      lzgsbyz0hq = log(zgs/z0hq)

      chqn=kappa*kappa/(lzgsbyz0m*lzgsbyz0hq)

      zgsbyl=zgs/lmonin
      z0hqbyl=z0hq/lmonin

c *********************************************************************
c  Now compute DPSI, which is the difference in the psi functions
c    computed at zgs and the relevant roughness height:
c *********************************************************************

      if (lmonin.gt.0.) then
c *********************************************************************
c  Here the atmosphere is stable with respect to the ground:
        dpsihq= sigma1*lzgsbyz0hq-sigma*gamahs*(zgsbyl-z0hqbyl)
c *********************************************************************
      else
c *********************************************************************
c  Here the atmosphere is unstable with respect to the ground:
        xhqs  =sqrt(1.-gamahu*zgsbyl)
        xhq0  =sqrt(1.-gamahu*z0hqbyl)
        dpsihq=sigma1*lzgsbyz0hq+2.*sigma*log((1.+xhqs)/(1.+xhq0))
c *********************************************************************
      endif

      dhq=sqrt(dm)/(1.-min(dpsihq/lzgsbyz0hq,.9d0))

      chq=dhq*chqn
      chq0=chq
      if (chq.gt.cmax) chq=cmax
      if (chq.lt.cmin) chq=cmin

      return
      end subroutine getchq

      subroutine simil(u,t,q,z,ustar,tstar,qstar,
     2                 z0m,z0h,z0q,lmonin,tg,qg)
!@sum   simil calculates the similarity solutions for wind speed,
!@+     virtual potential temperature, and moisture mixing ratio
!@+     at height z.
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var     z       height above ground at which solution is computed (m)
!@var     ustar   friction speed (m/sec)
!@var     tstar   temperature scale (K)
!@var     qstar   moisture scale
!@var     z0m     momentum roughness height (m)
!@var     z0h     temperature roughness height (m)
!@var     z0q     moisture roughness height (m)
!@var     lmonin  Monin-Obukhov length scale (m)
!@var     tg      ground temperature (K)
!@var     qg      ground moisture mixing ratio
!@var     u       computed similarity solution for wind speed (m/sec)
!@var     t       computed similarity solution for virtual potential
!@+               temperature (K)
!@var     q       computed similarity solution for moisture mixing ratio
      implicit none

      real*8,  intent(in) :: z,ustar,tstar,qstar,z0m,z0h,z0q
      real*8,  intent(in) :: lmonin,tg,qg
      real*8,  intent(out) :: u,t,q

      real*8 zbyl,z0mbyl,z0hbyl,z0qbyl,dpsim,dpsih,dpsiq,xm,xm0,xh,xh0
     *     ,xq,xq0,lzbyz0m,lzbyz0h,lzbyz0q

      zbyl  =z  /lmonin
      z0mbyl=z0m/lmonin
      z0hbyl=z0h/lmonin
      z0qbyl=z0q/lmonin
      lzbyz0m=log(z/z0m)
      lzbyz0h=log(z/z0h)
      lzbyz0q=log(z/z0q)
c *********************************************************************
c  Now compute DPSI, which is the difference in the psi functions
c    computed at zgs and the relevant roughness height:
c *********************************************************************

      if (lmonin.gt.0.) then
c *********************************************************************
c  Here the atmosphere is stable with respect to the ground:
        dpsim=-gamams*(zbyl-z0mbyl)
        dpsih= sigma1*lzbyz0h-sigma*gamahs*(zbyl-z0hbyl)
        dpsiq= sigma1*lzbyz0q-sigma*gamahs*(zbyl-z0qbyl)
c *********************************************************************
        else
c *********************************************************************
c  Here the atmosphere is unstable with respect to the ground:
        xm   =    (1.-gamamu*  zbyl)**0.25d0
        xm0  =    (1.-gamamu*z0mbyl)**0.25d0
        xh   =sqrt(1.-gamahu*  zbyl)
        xh0  =sqrt(1.-gamahu*z0hbyl)
        xq   =sqrt(1.-gamahu*  zbyl)
        xq0  =sqrt(1.-gamahu*z0qbyl)
        dpsim=log((1.+xm )*(1.+xm )*(1.+xm *xm )/
     2           ((1.+xm0)*(1.+xm0)*(1.+xm0*xm0)))-
     3        2.*(atan(xm )-atan(xm0))
        dpsih=sigma1*lzbyz0h+2.*sigma*log((1.+xh)/(1.+xh0))
        dpsiq=sigma1*lzbyz0q+2.*sigma*log((1.+xq)/(1.+xq0))
c *********************************************************************
      endif

      u=   (ustar/kappa)*(lzbyz0m-dpsim)
      t=tg+(tstar/kappa)*(lzbyz0h-dpsih)
      q=qg+(qstar/kappa)*(lzbyz0q-dpsiq)

      return
      end subroutine simil

      subroutine griddr(z,zhat,xi,xihat,dz,dzh,z1,zn,bgrid,n,ierr)
!@sum Computes altitudes on the vertical grid. The XI coordinates are
!@+   uniformly spaced and are mapped in a log-linear fashion onto the
!@+   Z grid. (The Z's are the physical coords.) Also computes the
!@+   altitudes on the secondary grid, ZHAT(I), and the derivatives
!@+   dxi/dz evaluated at both all Z(I) and ZHAT(I).
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
c     Grids:
c
c                n   - - - - - - - - - - - - -
c                    -------------------------  n-1
c                n-1 - - - - - - - - - - - - -
c                    -------------------------  j+1
c     (main,z)   j+1 - - - - - - - - - - - - -
c     (z)            -------------------------  j     (secondary, zhat)
c                j   - - - - - - - - - - - - -
c                    -------------------------  j-1
c                j-1 - - - - - - - - - - - - -
c                    -------------------------    2
c                2   - - - - - - - - - - - - -
c                    -------------------------    1
c                1   - - - - - - - - - - - - -
c
c     dz(j)==zhat(j)-zhat(j-1), dzh(j)==z(j+1)-z(j)
!@var bgrid determines how strongly non-linear the
!@+   mapping is. BGRID=0 gives linear mapping. Increasing BGRID
!@+   packs more points into the bottom of the layer.
!@+   Bgrid is calculated in the beginning of this subroutine every
!@+   time this subroutine is called. zs is the first grid height
!@+   and dzs is the grid separation near the surface. The values of
!@+   parameters byzs(=1/zs) and bydzs(=1/dzs) are obtained by
!@+   calling an old version of this subroutine with bgrid=0.2927
!@+   and has been tested as appropriate for pe(2)=934mb.
!@+   Now we impose that for other pe(2)s, zs and dzs be the same
!@+   as when pe(2)=934mb, so as to maintain the balance between
!@+   the accuracy and stability.
!@+   The new value of bgrid is then calculated below which
!@+   also depends on ztop (i.e., depends on pe(2)).
!@var  z       height of main grids (meter)
!@var  zhat    height of secondary grids (meter)
!@var  xi      an uniformly spaced coordinate mapped to z
!@var  xihat   an uniformly spaced coordinate mapped to zhat
!@var  dz   dxi/(dxi/dz)
!@var  dzh  dxi/(dxi/dzh)
!@var  dxi  (ztop - zbottom)/(n-1)
!@var  ierr Error reporting flag
      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(out) :: z,xi,dz
      real*8, dimension(n-1), intent(out) :: zhat,xihat,dzh
      integer, intent(out) :: ierr
      real*8, intent(in) :: z1,zn
      real*8, intent(out) :: bgrid

      real*8, parameter ::  tolz=1d-3
     &  ,byzs=1.d0/10.d0,bydzs=1.d0/4.7914d0
      real*8 z1pass,znpass,b,xipass,lznbyz1
      common /grids_99/z1pass,znpass,b,xipass,lznbyz1
!$OMP  THREADPRIVATE(/GRIDS_99/)
      external fgrid2
      real*8 rtsafe
      integer i,iter  !@var i,iter loop variable
      real*8 dxi,zmin,zmax,dxidz,dxidzh

      z1pass=z1
      znpass=zn
      dxi=(zn-z1)/float(n-1)
      bgrid=max((dxi*bydzs-1.)/((zn-z1)*byzs-log(zn/z1)),0.d0)

      b=bgrid
      zmin=z1
      zmax=zn

      do i=1,n-1
        xi(i)=z1+(zn-z1)*float(i-1)/float(n-1)
        xihat(i)=z1+(zn-z1)*(float(i)-0.5)/float(n-1)
      end do
      xi(n)=zn
      z(1)=z1

      lznbyz1 = log(zn/z1)
      dxidz=1.+bgrid*((zn-z1)/z1-lznbyz1)
      dz(1)=dxi/dxidz
      xipass=xihat(1)
      zhat(1)=rtsafe(fgrid2,zmin,zmax,tolz,ierr)
      if (ierr.gt.0) return
      dxidzh=1.+bgrid*((zn-z1)/zhat(1)-lznbyz1)
      dzh(1)=dxi/dxidzh

      do i=2,n-1
        xipass=xi(i)
        z(i)=rtsafe(fgrid2,zmin,zmax,tolz,ierr)
        if (ierr.gt.0) return
        xipass=xihat(i)
        zhat(i)=rtsafe(fgrid2,zmin,zmax,tolz,ierr)
        if (ierr.gt.0) return
        dxidz=1.+bgrid*((zn-z1)/z(i)-lznbyz1)
        dxidzh=1.+bgrid*((zn-z1)/zhat(i)-lznbyz1)
        dz(i)=dxi/dxidz
        dzh(i)=dxi/dxidzh
      end do
      z(n)=zn
      dxidz=1.+bgrid*((zn-z1)/zn-lznbyz1)
      dz(n)=dxi/dxidz

      return
      end subroutine griddr

      subroutine tfix(t,z,ttop,tgrnd,lmonin,tstar,ustar,khs,ts_guess,n)
!@sum   tfix
!@auth  Ye Cheng/G. Hartke
!@ver   1.0

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n),intent(in) :: z
      real*8, dimension(n),intent(inout) :: t
      real*8, intent(in) :: ttop,tgrnd,ustar,khs,ts_guess
      real*8, intent(inout) :: lmonin,tstar

      real*8 dtdz
      integer i  !@var i loop variable

      t(1)=ts_guess
      do i=2,n-1
        t(i)=t(1)+(z(i)-z(1))*(ttop-t(1))/(z(n)-z(1))
      end do

      dtdz   = (t(2)-t(1))/(z(2)-z(1))
      tstar  = khs*dtdz/ustar
      if (abs(tstar).gt.smax*abs(t(1)-tgrnd)) tstar=smax*(t(1)-tgrnd)
      if (abs(tstar).lt.smin*abs(t(1)-tgrnd)) tstar=smin*(t(1)-tgrnd)

      lmonin = ustar*ustar*tgrnd/(kappa*grav*tstar)
      if(abs(lmonin).lt.teeny) lmonin=sign(teeny,lmonin)

      return
      end subroutine tfix

      subroutine ccoeff0
!@sum   ccoeff0 sets/calculates model coefficients for the
!@+     Giss 2000 turbulence model (level2 2/2.5)
!@auth  Ye Cheng
!@ver   1.0
      implicit none

      ! temperary variable
#ifdef USE_PBL_E1
      real*8 :: del
#else
      real*8 :: del,aa,bb,cc,tmp
#endif

      prt=     0.82d0
      b1=     19.3d0
      b2=     15.8d0
      b123=b1**(2./3.)
      g1=       .1070d0
      g2=       .0032d0
      g3=       .0864d0
      g4=       .1000d0
      g5=     11.04d0
      g6=       .786d0
      g7=       .643d0
      g8=       .547d0
c
      d1=(7.*g4/3+g8)/g5
      d2=(g3**2-g2**2/3.)-1./(4.*g5**2)*(g6**2-g7**2)
      d3=g4/(3.*g5**2)*(4.*g4+3.*g8)
      d4=g4/(3.*g5**2)*(g2*g6-3.*g3*g7-g5*(g2**2-g3**2))
     &   +g8/g5*(g3**2-g2**2/3.)
      d5=-1./(4.*g5**2)*(g3**2-g2**2/3)*(g6**2-g7**2)
      s0=g1/2.
      s1=-g4/(3.*g5**2)*(g6+g7)+2.*g4/(3.*g5)*(g1-g2/3.-g3)
     &   +g1/(2.*g5)*g8
      s2=-g1/(8.*g5**2)*(g6**2-g7**2)
      s4=2./(3.*g5)
      s5=2.*g4/(3.*g5**2)
      s6=2./(3.*g5)*(g3**2-g2**2/3)-g1/(2.*g5)*(g3-g2/3.)
     &   +g1/(4*g5**2)*(g6-g7)

      s7=3.*g3-g2   !@ useful to calculate w2
      s8=4.*g4      !@ useful to calculate w2

c     find rimax:

      c1=s5+2*d3
      c2=s1-s6-2*d4
      c3=-s2+2.*d5
      c4=s4+2.*d1
      c5=-s0+2.*d2

      rimax=(c2+sqrt(c2**2-4.*c1*c3))/(2.*c1)
      rimax=int(rimax*1000.)/1000.

#ifdef USE_PBL_E1
      ! do nothing
#else
c     find gm_at_rimax

      aa=c1*rimax*rimax-c2*rimax+c3
      bb=c4*rimax+c5
      cc=2.d0
      if(abs(aa).lt.1d-8) then
         gm_at_rimax= -cc/bb
      else
         tmp=bb*bb-4.*aa*cc
         gm_at_rimax=(-bb-sqrt(tmp))/(2.*aa)
      endif
      ! rimax=.96,  gm_at_rimax=.1366285d6
#endif

c     find ghmin,ghmax,gmmax0:

      del=(s4+2.*d1)**2-8.*(s5+2.*d3)
      ghmin=(-s4-2.*d1+sqrt(del))/(2.*(s5+2.*d3))
      ghmin=int(ghmin*10000.)/10000.
      ghmax=(b1*0.53d0)**2
      gmmax0=(b1*0.34d0)**2  ! not in use yet

      ! for level 3 model only:
      g0=2./3.
      d1_3=(7.*g4/3.)/g5
      d2_3=d2
      d3_3=g4/(3.*g5**2)*(4.*g4)
      d4_3=g4/(3.*g5**2)*(g2*g6-3*g3*g7-g5*(g2**2-g3**2))
      d5_3=d5
      s0_3=s0
      s1_3=-g4/(3.*g5**2)*(g6+g7)+2.*g4/(3.*g5)*(g1-g2/3.-g3)
      s2_3=s2
      s3_3=g0*g4/g5*(g3+g2/3.+1./(2.*g5)*(g6+g7))
      s4_3=s4
      s5_3=s5
      s6_3=s6

      return
      end subroutine ccoeff0

      subroutine getk(km,kh,kq,ke,gma,gha,u,v,t,e,lscale,dzh,n)
!@sum   getk calculates eddy diffusivities Km, Kh and Ke
!@+     Giss 2000 turbulence model at level 2.5
!@auth  Ye Cheng
!@ver   1.0
c     Grids:
c
c                n   - - - - - - - - - - - - -
c                    -------------------------  n-1
c                n-1 - - - - - - - - - - - - -
c                    -------------------------  j+1
c     (main,z)   j+1 - - - - - - - - - - - - -
c     (z)            -------------------------  j     (secondary, zhat)
c                j   - - - - - - - - - - - - -
c                    -------------------------  j-1
c                j-1 - - - - - - - - - - - - -
c                    -------------------------    2
c                2   - - - - - - - - - - - - -
c                    -------------------------    1
c                1   - - - - - - - - - - - - -
c
c     dz(j)==zhat(j)-zhat(j-1), dzh(j)==z(j+1)-z(j)
c     at main: u,v,t,q,ke
c     at edge: e,lscale,km,kh,gm,gh
!@sum getk computes the turbulent viscosity, Km, and turbulent
!@+   conductivity, Kh, and turbulent diffusivity , Ke,
!@+   using the GISS second order closure model (2000)
!@+   at main: u,v,t,q,ke
!@+   at secondary: e,lscale,km,kh,gma,gha
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u,v,t,e,lscale,t_real z-profiles
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gma normalized velocity gradient, tau**2*as2
!@var gha normalized temperature gradient, tau**2*an2
!@var tau B1*lscale/sqrt(2*e)
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@var an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var se stability constant for e, adjustable

      USE CONSTANT, only : teeny

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t
      real*8, dimension(n-1), intent(in) :: e,lscale,dzh
      real*8, dimension(n-1), intent(out) :: km,kh,kq,ke,gma,gha

      real*8, parameter :: se=0.1d0,kmax=100.d0
     &  ,kmmin=1.5d-5,khmin=2.5d-5,kqmin=2.5d-5,kemin=1.5d-5
      real*8 :: an2,dudz,dvdz,as2,ell,den,qturb,tau,gh,gm,gmmax,sm,sh
     &  ,sq,taue
      integer :: i !@var i loop variable

      do i=1,n-1
        an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
        dudz=(u(i+1)-u(i))/dzh(i)
        dvdz=(v(i+1)-v(i))/dzh(i)
        as2=dudz*dudz+dvdz*dvdz
        ell=lscale(i)
        qturb=sqrt(2.*e(i))
        tau=B1*ell/max(qturb,teeny)
        gh=tau*tau*an2
        gm=tau*tau*as2
        if(gh.lt.ghmin) gh=ghmin
        if(gh.gt.ghmax) gh=ghmax
        gmmax=(1+d1*gh+d3*gh*gh)/(d2+d4*gh)
        if(gm.gt.gmmax) gm=gmmax
        den=1.+d1*gh+d2*gm+d3*gh*gh+d4*gh*gm+d5*gm*gm
        sm=(s0+s1*gh+s2*gm)/den
        sh=(s4+s5*gh+s6*gm)/den
        sq=sh
        taue=tau*e(i)
        km(i)=min(max(taue*sm,kmmin),kmax)
        kh(i)=min(max(taue*sh,khmin),kmax)
        kq(i)=min(max(taue*sq,kqmin),kmax)
        ke(i)=min(max(taue*se,kemin),kmax)
        gma(i)=gm
        gha(i)=gh
      end do
      return
      end subroutine getk

      subroutine e_eqn(esave,e,u,v,t,km,kh,ke,lscale,
     &                     dz,dzh,ustar,dtime,n)
!@sum e_eqn integrates differential eqn for e (tridiagonal method)
!@+   between the surface and the first GCM layer.
!@+   The boundary conditions at the bottom are:
!@+   e(1)=(1/2)*B1**(2/3)*ustar**2
!@+   at the top, dedz is continuous
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var ke z-profile of turbulent diffusion in eqn for e
!@var lscale z-profile of the turbulent dissipation length scale
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var dtime time step
!@var ustar friction velocity at the surface
!@var n number of vertical subgrid main layers
      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, intent(in) :: ustar, dtime
      real*8, dimension(n), intent(in) :: u,v,t,dz
      real*8, dimension(n-1), intent(in) :: esave,km,kh,ke,lscale,dzh
      real*8, dimension(n-1), intent(inout) :: e

      real*8 :: an2,dudz,dvdz,as2,qturb,ri,aa,bb,cc,gm,tmp
      integer :: j !@var j loop variable
c
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     from kirk:/u/acyxc/papers/2ndOrder/maple/phik.1,
c       sub(j)=-dtime*aj*P1a_jm1_k/dxi^2
c       sup(j)=-dtime*aj*P1a_jp1_k/dxi^2
c       dia(j)=1-(sub(j)+dia(j))+dtime*P3_j_k
c       rhs(j)=PHI_j_k + dtime*P4_j_k
c       where P1a == P1*a
c       aj=dxi/dz(j) if j refers to the primary grid
c       aj=dxi/dzh(j) if j refers to the secondary grid
c       now for e_eqn j refers to the secondary grid, therefore:
c       aj = dxi/dzh(j)
c       a(j-1/2) = dxi/dzh(j-1/2)=dxi/dz(j)
c       a(j+1/2) = dxi/dzh(j+1/2)=dxi/dz(j+1)
c       p1(j-1/2)=0.5*(ke(j)+ke(j-1))
c       p1(j+1/2)=0.5*(ke(j)+ke(j+1))
c       ke(j)=sq*lscale(j)*qturb, qturb=sqrt(2.*e(j))
c
      do j=2,n-2
          qturb=sqrt(2.*e(j))
          sub(j)=-dtime*0.5d0*(ke(j)+ke(j-1))/(dzh(j)*dz(j))
          sup(j)=-dtime*0.5d0*(ke(j)+ke(j+1))/(dzh(j)*dz(j+1))
          dia(j)=1.-(sub(j)+sup(j))+dtime*2*qturb/(b1*lscale(j))
          an2=2*grav*(t(j+1)-t(j))/((t(j+1)+t(j))*dzh(j))
          dudz=(u(j+1)-u(j))/dzh(j)
          dvdz=(v(j+1)-v(j))/dzh(j)
          as2=dudz*dudz+dvdz*dvdz
          rhs(j)=esave(j)+dtime*(km(j)*as2-kh(j)*an2)
       end do

      dia(1)=1.
      sup(1)=0.
      rhs(1)=0.5d0*b123*ustar*ustar

      j=n-1
      an2=2.*grav*(t(j+1)-t(j))/((t(j+1)+t(j))*dzh(j))
      dudz=(u(j+1)-u(j))/dzh(j)
      dvdz=(v(j+1)-v(j))/dzh(j)
      as2=max(dudz*dudz+dvdz*dvdz,teeny)
      ri=an2/as2
      if(ri.gt.rimax) ri=rimax
      aa=c1*ri*ri-c2*ri+c3
      bb=c4*ri+c5
      cc=2.d0
      if(abs(aa).lt.1d-8) then
        gm= -cc/bb
      else
        tmp=bb*bb-4.*aa*cc
        gm=(-bb-sqrt(tmp))/(2.*aa)
      endif
      sub(n-1)=0.
      dia(n-1)=1.
      rhs(n-1)=max(0.5d0*(B1*lscale(j))**2*as2/max(gm,teeny),teeny)

c     sub(n-1)=-1.
c     dia(n-1)=1.
c     rhs(n-1)=0.

      call TRIDIAG(sub,dia,sup,rhs,e,n-1)

      do j=1,n-1
         e(j)=min(max(e(j),teeny),emax)
      end do

      Return
      end subroutine e_eqn

      subroutine e_les(tstar,ustar,wstar3,dbl,lmonin,zhat,lscale,e,n)
!@sum e_gcm finds e according to the parameterization of les data
!@Ref Moeng and Sullivan 1994, J. Atmos. Sci., 51, 999-1022.
!@Ref Cheng et al. 2002, J. Atmos. Sci., 59, 1550-1565.
!@auth  Ye Cheng
!@ver   1.0
!@var (see subroutine k_gcm)
      USE CONSTANT, only : by3

      implicit none

      integer, intent(in) :: n   !@var n  array dimension
      real*8, intent(in) :: tstar,ustar,wstar3,dbl,lmonin
      real*8, dimension(n), intent(in) :: zhat,lscale
      real*8, dimension(n), intent(inout) :: e
      real*8, parameter :: emin=1.d-6
      integer :: j !@var j loop variable
      real*8 :: tvflx,ustar3,zj,kz,zeta,phi_m,eps,ej

      tvflx=ustar*tstar
      ustar3=ustar*ustar*ustar
      do j=1,n-1   ! Dyer 1974
        zj=zhat(j)
        kz=kappa*zj
        if(zj.le.dbl) then
          zeta=zj/lmonin
          if(zeta.ge.0.) then ! stable or neutral
            if(zeta.le.1.) then
              phi_m=1.+5.*zeta
            else
              phi_m=5.+zeta
            endif
          else                ! unstable
            phi_m=(1.-15.*zeta)**(-.25d0)
          endif
          eps=.4d0*wstar3/dbl+ustar3*(1.-zj/dbl)*phi_m/kz
          ej=.5d0*(24.d0*lscale(j)*eps)**(2.*by3)
          ej=min(max(ej,emin),emax)
        else
          ej=0.
        endif
        e(j)=max(e(j),ej)
      end do
      return
      end subroutine e_les

      subroutine t_eqn(u,v,t0,t,q,z,kh,kq,dz,dzh,ch,usurf,tgrnd
     &     ,ttop,dtime,n
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0,usurf0,tprime,tdns,ddml_eq_1)
!@sum t_eqn integrates differential eqn for t (tridiagonal method)
!@+   between the surface and the first GCM layer.
!@+   Boundary conditions at bottom: Mellor and Yamada 1982, Eq(72),
!@+   Redelsperger et al. 2000, J. Climate, 13, 402-421
!@+   Emanuel and Zivkovic 1999, JAS, 56, 1766-1782
!@+   including the effects on the surface flux
!@+   due to the moist convection wind gustiness and the
!@+   downdraft temperature perturbation
!@+   kh * dt/dz = ch * ( usurf*(t1 - tgrnd)
!@+                      +(1+deltx*q1)*(usurf-usurf0)*tprime )
!@+                + deltx * t1/(1+deltx*q1) * kq * dqdz
!@+   where tprime=tdns-t1/(1+deltx*q1), t1 is at surf
!@+   at the top, the virtual potential temperature is prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature ref. to the surface
!@var q z-profle of specific humidity
!@var t0 z-profle of t at previous time step
!@var kh z-profile of heat conductivity
!@var kq z-profile of moisture diffusivity
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var usurf effective surface velocity
!@var tgrnd virtual potential temperature at the ground
!@var ttop virtual potential temperature at the first GCM layer
!@var dtime time step
!@var n number of vertical subgrid main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: u,v,t0,q,z,dz
      real*8, dimension(n-1), intent(in) :: dzh,kh,kq
      real*8, dimension(n), intent(inout) :: t
      real*8, intent(in) :: ch,tgrnd
      real*8, intent(in) :: ttop,dtime,usurf
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      real*8, intent(in) ::  usurf0,tprime,tdns
      logical, intent(in) :: ddml_eq_1

      real*8 :: facth,factx,facty,rat
      integer :: i,j,iter  !@var i,j,iter loop variable

      do i=2,n-1
         sub(i)=-dtime/(dz(i)*dzh(i-1))*kh(i-1)
         sup(i)=-dtime/(dz(i)*dzh(i))*kh(i)
         dia(i)=1.-(sub(i)+sup(i))
      end do

cccc
c changes for using aturb with SCM
#ifdef SCM
c     factx=(dpdxr-dpdxr0)/(z(n)-z(1))
c     facty=(dpdyr-dpdyr0)/(z(n)-z(1))
      do i=2,n-1
c       rhs(i)=t0(i)-dtime*t(i)*bygrav*(v(i)*facty+u(i)*factx)
        rhs(i) = t0(i)
      end do
#else
      factx=(dpdxr-dpdxr0)/(z(n)-z(1))
      facty=(dpdyr-dpdyr0)/(z(n)-z(1))
      do i=2,n-1
        rhs(i)=t0(i)-dtime*t(i)*bygrav*(v(i)*facty+u(i)*factx)
      end do
#endif

      facth  = ch*usurf*dzh(1)/kh(1)
      sup(1) = -1.

      if(ddml_eq_1) then
         rat = usurf0/(usurf+teeny)
         dia(1) = 1+facth*rat
     &             +deltx*kq(1)*(q(2)-q(1))/(kh(1)*(1.+deltx*q(1)))
         rhs(1) = facth*(tgrnd-(1.d0-rat)*(1.+deltx*q(1))*tdns)
      else
         dia(1) = 1+facth
     &             +deltx*kq(1)*(q(2)-q(1))/(kh(1)*(1.+deltx*q(1)))
         rhs(1) = facth*tgrnd
      endif

      dia(n) = 1.
      sub(n) = 0.
      rhs(n) = ttop

      call TRIDIAG(sub,dia,sup,rhs,t,n)

      return
      end subroutine t_eqn

      subroutine q_eqn(q0,q,kq,dz,dzh,cq,usurf,qgrnd,qtop,dtime,n
     &     ,flux_max,fr_sat,usurf0,qprime,qdns,ddml_eq_1)
!@sum q_eqn integrates differential eqn q (tridiagonal method)
!@+   between the surface and the first GCM layer.
!@+   The boundary conditions at the bottom Ref:
!@+   Redelsperger et al. 2000, J. Climate, 13, 402-421
!@+   Emanuel and Zivkovic 1999, JAS, 56, 1766-1782
!@+   including the effects on the surface flux
!@+   due to the moist convection wind gustiness and the
!@+   downdraft specific humidity perturbation
!@+   kq * dq/dz = min ( cq * usurf * (q1 - qgrnd)
!@+                    + cq * (usurf-usurf0) * qprime ,
!@+           fr_sat * ( cq * usurf * (q1 - qgrnd)
!@+                    + cq * (usurf-usurf0) * qprime )
!@+       - ( 1 - fr_sat ) * flux_max )
!@+   where qprime=qdns-q1, q1 is q at surf
!@+   at the top, the moisture is prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var q z-profle of specific humidity
!@var q0 z-profle of q at previous time step
!@var kq z-profile of moisture diffusivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var usurf effective surface velocity
!@var qgrnd specific humidity at the ground
!@var qtop specific humidity at the first GCM layer
!@var dtime time step
!@var n number of vertical subgrid main layers
!@var flux_max maximal flux from the unsaturated soil
!@var fr_sat fraction of the saturated soil

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: q0,dz
      real*8, dimension(n-1), intent(in) :: dzh,kq
      real*8, dimension(n), intent(out) :: q
      real*8, intent(in) :: cq,qgrnd,qtop,dtime,usurf
      real*8, intent(in) :: flux_max,fr_sat
      real*8, intent(in) :: usurf0,qprime,qdns
      logical, intent(in) :: ddml_eq_1

      real*8 :: factq,rat
      integer :: i  !@var i loop variable

      do i=2,n-1
         sub(i)=-dtime/(dz(i)*dzh(i-1))*kq(i-1)
         sup(i)=-dtime/(dz(i)*dzh(i))*kq(i)
         dia(i)=1.-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=q0(i)
      end do

      factq  = cq*usurf*dzh(1)/kq(1)
      sup(1) = -1.

      if(ddml_eq_1) then
         rat = usurf0/(usurf+teeny)
         dia(1) = 1.+factq*rat
         rhs(1)= factq*(qgrnd-(1.-rat)*qdns)
      else
         dia(1) = 1.+factq
         rhs(1) = factq*qgrnd
      endif

      dia(n) = 1.
      sub(n) = 0.
      rhs(n) = qtop

      call TRIDIAG(sub,dia,sup,rhs,q,n)

c**** Now let us check if the computed flux doesn't exceed the maximum
c**** for unsaturated fraction

      if ( fr_sat .ge. 1. ) return   ! all soil is saturated
      if ( cq * usurf * (qgrnd - q(1)) + cq * (usurf0-usurf) * qprime
     &   .le. flux_max ) return

c**** Flux is too high, have to recompute with the following boundary
c**** conditions at the bottom:
c**** kq * dq/dz = fr_sat * ( cq * usurf * (q(1) - qgrnd)
c****                       + cq * (usurf-usurf0) * qprime )
c****             - ( 1 - fr_sat ) * flux_max

      if(ddml_eq_1) then
         rat = usurf0/(usurf+teeny)
         dia(1) = 1. + fr_sat*factq*rat
         rhs(1)= fr_sat*factq*(qgrnd-(1.-rat)*qdns)
     &            + (1.-fr_sat)*flux_max*dzh(1)/kq(1)
      else
         dia(1) = 1. + fr_sat*factq
         rhs(1)= fr_sat*factq*qgrnd
     &            + (1.-fr_sat)*flux_max*dzh(1)/kq(1)
      endif

      call TRIDIAG(sub,dia,sup,rhs,q,n)

      return
      end subroutine q_eqn

      subroutine tr_eqn(tr0,tr,kq,dz,dzh,sfac,constflx,trtop,
#ifdef TRACERS_WATER
     *     tr_evap_max,fr_sat,
#endif
     *     dtime,n)
!@sum tr_eqn integrates differential eqn for tracers (tridiag. method)
!@+   between the surface and the first GCM layer.
!@+   The boundary conditions at the bottom are:
!@+   kq * dtr/dz = sfac * trs - constflx
!@+   i.e. for moisture, sfac=cq*usurf, constflx=cq*usurf*qg
!@+        to get:  kq * dq/dz = cq * usurf * (qs - qg)
!@+   for new moisture (including downdraft effects)
!@+        sfac=cq*(usurf-dusurf), constflx=cq*(usurf*qg + dusurf*qdns)
!@+    or  sfac=cq*usurf0, constflx=cq*(usurf*(qg+qdns)-usurf0*qdns)
!@+        to get:  kq * dq/dz = cq*(usurf*(qs-qg) + dusurf*(qdns-qs))
!@+   This should be flexible enough to deal with most situations.
!@+   at the top, the tracer conc. is prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var tr z-profle of tracer concentration
!@var tr0 z-profle of tr at previous time step
!@var kq z-profile of moisture diffusivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var sfac factor multiplying surface tracer conc. in b.c.
!@var constflx the constant component of the surface tracer flux
!@var trtop tracer concentration at the first GCM layer
!@var dtime time step
!@var n number of vertical subgrid main layers
#ifdef TRACERS_WATER
!@var tr_evap_max  max amount of possible tracer evaporation
!@var fr_sat fraction of the saturated soil
#endif
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: tr0,dz
      real*8, dimension(n-1), intent(in) :: dzh,kq
      real*8, dimension(n), intent(out) :: tr
      real*8, intent(in) :: sfac,constflx,trtop,dtime
#ifdef TRACERS_WATER
      real*8, intent(in) :: tr_evap_max, fr_sat
#endif
      real*8 :: facttr
      integer :: i  !@var i loop variable

      do i=2,n-1
        sub(i)=-dtime/(dz(i)*dzh(i-1))*kq(i-1)
        sup(i)=-dtime/(dz(i)*dzh(i))*kq(i)
        dia(i)=1.-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=tr0(i)
      end do

      facttr  = dzh(1)/kq(1)

      dia(1) = 1.+sfac*facttr
      sup(1) = -1.
      rhs(1)= facttr*constflx

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = trtop

#ifdef WATER_PROPORTIONAL
      if( .not. force_limit ) then
#endif
      call TRIDIAG(sub,dia,sup,rhs,tr,n)
#ifdef WATER_PROPORTIONAL
      endif
#endif

#ifdef TRACERS_WATER
c**** Check as in q_eqn if flux is limited, and if it is, recalculate
c**** profile

      if ( fr_sat .ge. 1. ) return   ! all soil is saturated

#ifdef WATER_PROPORTIONAL
      if (.not.force_limit) then
#endif

      if ( constflx - sfac * tr(1) .le. tr_evap_max ) return

#ifdef WATER_PROPORTIONAL
      endif
      force_limit = .true.
#endif

c**** Flux is too high, have to recompute with the following boundary
c**** conditions at the bottom:
c**** kq * dq/dz = fr_sat * (sfac * trs - constflx)
c****              + ( 1 - fr_sat ) * tr_evap_max

      dia(1) = 1. + fr_sat*sfac*facttr
      sup(1) = -1.
      rhs(1)= fr_sat*facttr*constflx +
     *     (1.-fr_sat)*tr_evap_max*dzh(1)/kq(1)

      call TRIDIAG(sub,dia,sup,rhs,tr,n)
#endif

      return
      end subroutine tr_eqn

      subroutine uv_eqn(u0,v0,u,v,z,km,dz,dzh,
     2                  ustar,cm,z0m,utop,vtop,dtime,coriol,
     3                  ug,vg,uocean,vocean,n
     &                  ,dpdxr,dpdyr,dpdxr0,dpdyr0)
!@sum uv_eqn integrates differential eqns for u & v (tridiagonal method)
!@+   between the surface and the first GCM layer.
!@+   The boundary conditions at the bottom are:
!@+   km * du/dz = cm * usurf * u
!@+   km * dv/dz = cm * usurf * v
!@+    at the top, the winds are prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var u0 z-profle of u at previous time step
!@var v0 z-profle of v at previous time step
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ustar friction velocity at the surface
!@var cm  dimensionless  momentum flux at surface (drag coeff.)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var z0m  roughness height for momentum (if itype=1 or 2)
!@var utop due east component of the wind at the first GCM layer
!@var vtop due north component of the wind at the first GCM layer
!@var dtime time step
!@var coriol the Coriolis parameter
!@var ug due east component of the geostrophic wind
!@var vg due north component of the geostrophic wind
!@var n number of vertical subgrid main layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1

      real*8, dimension(n), intent(in) :: u0,v0,z,dz
      real*8, dimension(n), intent(inout) :: u,v
      real*8, dimension(n-1), intent(in) :: km,dzh
      real*8, intent(in) :: ustar,cm,z0m,utop,vtop,dtime,coriol,ug,vg
     &                     ,uocean,vocean
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0

      real*8 :: factx,facty,dpdx,dpdy,usurf,factor
      integer :: i,j,iter  !@var i,j,iter loop variable

      do i=2,n-1
         sub(i)=-dtime/(dz(i)*dzh(i-1))*km(i-1)
         sup(i)=-dtime/(dz(i)*dzh(i))*km(i)
         dia(i)=1.-(sub(i)+sup(i))
      end do

#ifndef SCM
      factx=(dpdxr-dpdxr0)/(z(n)-z(1))
      facty=(dpdyr-dpdyr0)/(z(n)-z(1))
#endif
      do i=2,n-1
ccc for SCM use ug and vg and not dpdx,dpdy
#ifdef SCM
        rhs(i)=u0(i)+dtime*coriol*(v(i)-vg)
        rhs1(i)=v0(i)-dtime*coriol*(u(i)-ug)
#else
        dpdx=factx*(z(i)-z(1))+dpdxr0
        dpdy=facty*(z(i)-z(1))+dpdyr0
        rhs(i)=u0(i)+dtime*(coriol*v(i)-dpdx)
        rhs1(i)=v0(i)-dtime*(coriol*u(i)+dpdy)
#endif
      end do

      usurf  = sqrt((u(1)-uocean)**2+(v(1)-vocean)**2)
      factor = cm*usurf*dzh(1)/km(1)

      dia(1) = 1.+factor
      sup(1) = -1.
      rhs(1)  = factor*uocean
      rhs1(1) = factor*vocean

      dia(n) = 1.
      sub(n) = 0.
      rhs(n)  = utop
      rhs1(n)  = vtop

      call TRIDIAG(sub,dia,sup,rhs,u,n)
      call TRIDIAG(sub,dia,sup,rhs1,v,n)

      return
      end subroutine uv_eqn

      subroutine t_eqn_sta(t,q,kh,kq,dz,dzh,ch,usurf,tgrnd,ttop,n)
!@sum  t_eqn_sta computes the static solutions of t
!@+    between the surface and the first GCM layer.
!@+    The boundary conditions at the bottom are:
!@+    kh * dt/dz = ch * usurf * (t - tg)
!@+    at the top, virtual potential temperature is prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var q z-profle of specific humidity
!@var kh z-profile of heat conductivity
!@var kq z-profile of specific humidity
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var usurf effective surface velocity
!@var tgrnd virtual potential temperature at the ground
!@var ttop virtual potential temperature at the first GCM layer
!@var n number of vertical main subgrid layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: q,dz
      real*8, dimension(n), intent(inout) :: t
      real*8, dimension(n-1), intent(in) :: kh,kq,dzh
      real*8, intent(in) :: ch,tgrnd,ttop,usurf

      real*8 :: facth
      integer :: i !@var i loop variable

      do i=2,n-1
         sub(i)=-1./(dz(i)*dzh(i-1))*kh(i-1)
         sup(i)=-1./(dz(i)*dzh(i))*kh(i)
         dia(i)=-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=0.
      end do

      facth  = ch*usurf*dzh(1)/kh(1)

      dia(1) = 1.+facth
     &        +deltx*kq(1)/kh(1)*(q(2)-q(1))/(1.+deltx*q(1))
      sup(1) = -1.
      rhs(1) = facth*tgrnd

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = ttop

      call TRIDIAG(sub,dia,sup,rhs,t,n)

      return
      end subroutine t_eqn_sta

      subroutine q_eqn_sta(q,kq,dz,dzh,cq,usurf,qgrnd,qtop,n)
!@sum  q_eqn_sta computes the static solutions of q
!@+    between the surface and the first GCM layer.
!@+    The boundary conditions at the bottom are:
!@+    kq * dq/dz = cq * usurf * (q - qg)
!@+    at the top, moisture is prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var q z-profle of specific humidity
!@var kq z-profile of moisture diffusivity
!@var z vertical grids (main, meter)
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var usurf effective surface velocity
!@var qgrnd specific humidity at the ground
!@var qtop specific humidity at the first GCM layer
!@var n number of vertical main subgrid layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs

      real*8, dimension(n), intent(in) :: dz
      real*8, dimension(n), intent(inout) :: q
      real*8, dimension(n-1), intent(in) :: kq,dzh
      real*8, intent(in) :: cq,qgrnd,qtop,usurf

      real*8 :: factq
      integer :: i  !@var i loop variable

      do i=2,n-1
         sub(i)=-1./(dz(i)*dzh(i-1))*kq(i-1)
         sup(i)=-1./(dz(i)*dzh(i))*kq(i)
         dia(i)=-(sub(i)+sup(i))
      end do

      do i=2,n-1
        rhs(i)=0.
      end do

      factq  = cq*usurf*dzh(1)/kq(1)

      dia(1) = 1.+factq
      sup(1) = -1.
      rhs(1)= factq*qgrnd

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = qtop

      call TRIDIAG(sub,dia,sup,rhs,q,n)

      return
      end subroutine q_eqn_sta

      subroutine uv_eqn_sta(u,v,z,km,dz,dzh,
     2            ustar,cm,utop,vtop,coriol,uocean,vocean,n
     &            ,dpdxr,dpdyr,dpdxr0,dpdyr0,ug,vg)
!@sum  uv_eqn_sta computes the static solutions of the u and v
!@+    between the surface and the first GCM layer.
!@+    The boundary conditions at the bottom are:
!@+    km * du/dz = cm * usurf * u
!@+    km * dv/dz = cm * usurf * v
!@+    at the top, the winds are prescribed.
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var z vertical height at the main grids (meter)
!@var zhat vertical height at the secondary grids (meter)
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var dz(j) zhat(j)-zhat(j-1)
!@var dzh(j)  z(j+1)-z(j)
!@var ustar friction velocity at the surface
!@var cm  dimensionless  momentum flux at surface (drag coeff.)
!@var ch  dimensionless heat flux at surface (stanton number)
!@var cq  dimensionless moisture flux at surface (dalton number)
!@var utop due east component of the wind at the first GCM layer
!@var vtop due north component of the wind at the first GCM layer
!@var coriol the Coriolis parameter
!@var n number of vertical subgrid main layers
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1

      real*8, dimension(n), intent(in) :: z,dz
      real*8, dimension(n), intent(inout) :: u,v
      real*8, dimension(n-1), intent(in) :: km,dzh
      real*8, intent(in) :: ustar,cm,utop,vtop,coriol,uocean,vocean
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0

      real*8 :: factx,facty,dpdx,dpdy,usurf,factor
      integer :: i,j,iter  !@var i,j,iter loop variable
c**** passed for SCM
      real*8 ug,vg

      do i=2,n-1
          sub(i)=-1./(dz(i)*dzh(i-1))*km(i-1)
          sup(i)=-1./(dz(i)*dzh(i))*km(i)
          dia(i)=-(sub(i)+sup(i))
       end do

      factx=(dpdxr-dpdxr0)/(z(n)-z(1))
      facty=(dpdyr-dpdyr0)/(z(n)-z(1))
c      write(99,*) factx,facty
      do i=2,n-1
ccc if running SCM then use ug and vg instead of dpdx,dpdy
#ifdef SCM
            rhs(i)=coriol*(v(i)-vg)
            rhs1(i)=-coriol*(u(i)-ug)
#else
            dpdx=factx*(z(i)-z(1))+dpdxr0
            dpdy=facty*(z(i)-z(1))+dpdyr0
            rhs(i)=(coriol*v(i)-dpdx)
            rhs1(i)=-(coriol*u(i)+dpdy)
#endif
      end do

      usurf  = sqrt((u(1)-uocean)**2+(v(1)-vocean)**2)
      factor = cm*usurf*dzh(1)/km(1)

      dia(1) = 1.+factor
      sup(1) = -1.
      rhs(1) = factor*uocean
      rhs1(1)= factor*vocean

      dia(n) = 1.
      sub(n) = 0.
      rhs(n)  = utop
      rhs1(n)  = vtop

      call TRIDIAG(sub,dia,sup,rhs,u,n)
      call TRIDIAG(sub,dia,sup,rhs1,v,n)

      return
      end subroutine uv_eqn_sta

      subroutine level2(e,u,v,t,lscale,dzh,n)
!@sum  level2 computes the turbulent kinetic energy e using the
!@+    GISS 2000 turbulence model at level 2
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var e z-profle of turbulent kinetic energy
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of virtual potential temperature
!@var lscale z-profile of the turbulent dissipation length scale
!@var z vertical grids (main, meter)
!@var dzh(j)  z(j+1)-z(j)
!@var n number of vertical subgrid main layers

      USE CONSTANT, only : teeny

      implicit none

      integer, intent(in) :: n     !@var n array dimension
      real*8, dimension(n), intent(in)   :: u,v,t
      real*8, dimension(n-1), intent(in) :: lscale,dzh
      real*8, dimension(n-1), intent(out) :: e
      real*8 :: dudz,dvdz,as2,an2,ri,gm
      real*8 :: aa,bb,cc,tmp
      integer :: i    !@var i loop variable

      do i=1,n-1
        an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
        dudz=(u(i+1)-u(i))/dzh(i)
        dvdz=(v(i+1)-v(i))/dzh(i)
        as2=max(dudz*dudz+dvdz*dvdz,teeny)
        ri=an2/as2
        if(ri.gt.rimax) ri=rimax
        aa=c1*ri*ri-c2*ri+c3
        bb=c4*ri+c5
        cc=2.d0
        if(abs(aa).lt.1d-8) then
          gm= -cc/bb
        else
          tmp=bb*bb-4.*aa*cc
          gm=(-bb-sqrt(tmp))/(2.*aa)
        endif
        e(i)=0.5d0*(B1*lscale(i))**2*as2/max(gm,teeny)
        e(i)=min(max(e(i),teeny),emax)
      end do

      return
      end subroutine level2

      subroutine inits(tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2          ttop,qtop,coriol,cm,ch,cq,ustar,uocean,vocean,
     3          ilong,jlat,itype
     &          ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &          ,u,v,t,q,e,ug,vg)
!@sum  inits initializes the winds, virtual potential temperature,
!@+    and humidity using static solutions of the GISS 2000
!@+    turbulence model at level 2
!@var  n number of sub-grid levels for the PBL
!@var  tgrnd virtual potential temperature of ground,at roughness height
!@var  qgrnd  moisture at the ground, at the roughness height
!@var  zgrnd
!@var  zgs  height of the surface layer (nominally 10 m)
!@var  ztop height of the first model layer, approx 200 m if lm=9
!@var  utop  x component of wind at the top of the layer
!@var  vtop  y component of wind at the top of the layer
!@var  ttop  virtual potential temperature at the top of the layer
!@var  qtop  moisture at the top of the layer
!@var  coriol the Coriolis parameter
!@var  cm  dimensionless  momentum flux at surface (drag coeff.)
!@var  ch  dimensionless heat flux at surface (stanton number)
!@var  cq  dimensionless moisture flux at surface (dalton number)
!@var  bgrid log-linear gridding parameter
!@var  ustar  friction speed
!@var  ilong  longitude identifier
!@var  jlat  latitude identifier
!@var  itype  surface type
!@var  ug,vg  passed for scm

!@var  iprint longitude for diagnostics
!@var  jprint latitude for diagnostics
      implicit none

      real*8, intent(in) :: tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,ttop
     *     ,qtop,coriol,uocean,vocean
      real*8, intent(out) :: cm,ch,cq,ustar
      integer, intent(in) :: ilong,jlat,itype
      real*8, intent(in) ::  dpdxr,dpdyr,dpdxr0,dpdyr0

      real*8, dimension(n-1) :: km,kh,kq,ke,gm,gh
      real*8, dimension(n) :: z,dz,xi,usave,vsave,tsave,qsave
      real*8, dimension(n-1) :: zhat,xihat,dzh,lscale,esave
      real*8 :: lmonin,bgrid,z0m,z0h,z0q,hemi,psi1,psi0,psi
     *     ,usurf,tstar,qstar,ustar0,dtime,test
     *     ,wstar3,wstar2h,usurfq,usurfh,ts

      integer, save :: iter_count=0
      integer, parameter ::  itmax=100
      integer, parameter ::  iprint=0,jprint=41 ! set iprint>0 to debug
      real*8, parameter ::  w=0.50,tol=1d-3
      integer :: i,j,iter,ierr  !@var i,j,iter loop variable
#ifdef TRACERS_SPECIAL_O18
      real*8 :: fac_cq_tr(ntm)   ! not used here
#endif

      real*8 dbl ! I hope it is really a local variable (was global before) I.A

c**** special threadprivate common block (compaq compiler stupidity)
      real*8, dimension(n), intent(out) :: u,v,t,q
      real*8, dimension(n-1), intent(out) :: e
c****  passed for scm
      real*8  ug,vg

!!      common/pbluvtq/u,v,t,q,e

!!!$OMP  THREADPRIVATE (/pbluvtq/)
C**** end special threadprivate common block


      dbl=1000.d0 !initial guess of dbl
      z0m=zgrnd
      z0h=z0m
      z0q=z0m

      call griddr(z,zhat,xi,xihat,dz,dzh,zgs,ztop,bgrid,n,ierr)
      if (ierr.gt.0) then
        print*,"In inits: i,j,itype =",ilong,jlat,itype,tgrnd,qgrnd
        call stop_model("PBL error in inits",255)
      end if

c Initialization for iteration:
      if (coriol.le.0.) then
        hemi=-1.
        else
        hemi= 1.
      endif
      if (tgrnd.gt.ttop) then
        psi1=hemi*5.*radian
        else
        psi1=hemi*15.*radian
      endif
      psi0=atan2(vtop,utop+teeny)
      if(psi0.lt.0.) psi0=psi0+2.*pi
      psi=psi0+psi1
      usurf=0.4d0*sqrt(utop*utop+vtop*vtop)
      u(1)=usurf*cos(psi)
      v(1)=usurf*sin(psi)
      t(1)=tgrnd-0.5d0*(tgrnd-ttop)
      q(1)=qgrnd-0.5d0*(qgrnd-qtop)
      e(1)=1.d-3
        usave(1)=u(1)
        vsave(1)=v(1)
        tsave(1)=t(1)
        qsave(1)=q(1)
        esave(1)=e(1)
      u(n)=utop
      v(n)=vtop
      t(n)=ttop
      q(n)=qtop
      e(n-1)=1.d-3

      do i=2,n-1
        u(i)=u(1)+(u(n)  -u(1))*((z(i)-z(1))/(z(n)  -z(1)))
        v(i)=v(1)+(v(n)  -v(1))*((z(i)-z(1))/(z(n)  -z(1)))
        t(i)=t(1)+(t(n)  -t(1))*((z(i)-z(1))/(z(n)  -z(1)))
        q(i)=q(1)+(q(n)  -q(1))*((z(i)-z(1))/(z(n)  -z(1)))
        e(i)=e(1)+(e(n-1)-e(1))*((z(i)-z(1))/(z(n-1)-z(1)))
        usave(i)=u(i)
        vsave(i)=v(i)
        tsave(i)=t(i)
        qsave(i)=q(i)
        esave(i)=e(i)
      end do

      ustar0=0.
      do iter=1,itmax

        if(iter.eq.1) then
          call getl1(e,zhat,dzh,lscale,n)
        else
          call getl(e,u,v,t,zhat,dzh,lmonin,ustar,lscale,dbl,n)
        endif

        call getk(km,kh,kq,ke,gm,gh,u,v,t,e,lscale,dzh,n)

        ts=t(1)/(1+q(1)*deltx)
        call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,ts,
     2             u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
#ifdef TRACERS_SPECIAL_O18
     *             fac_cq_tr,
#endif
     3             km,kh,kq,dzh,itype,n)
        !! dbl=.375d0*sqrt(ustar*abs(lmonin)/omega)
        !@+ M.J.Miller et al. 1992, J. Climate, 5(5), 418-434, Eqs(6-7),
        !@+ for heat and mositure

        if(t(2).lt.t(1)) then
          wstar3=-dbl*grav*2.*(t(2)-t(1))*kh(1)/((t(2)+t(1))*dzh(1))
          wstar2h = wstar3**twoby3
        else
          wstar3=0.
          wstar2h = 0.
        endif

        usurfh  = sqrt((u(1)-uocean)**2+(v(1)-vocean)**2+wstar2h)
        usurfq  = usurfh

        call t_eqn_sta(t,q,kh,kq,dz,dzh,ch,usurfh,tgrnd,ttop,n)

        call q_eqn_sta(q,kq,dz,dzh,cq,usurfq,qgrnd,qtop,n)

        call uv_eqn_sta(u,v,z,km,dz,dzh,ustar,cm,utop,vtop,coriol,
     &                  uocean,vocean,n
     &                  ,dpdxr,dpdyr,dpdxr0,dpdyr0,ug,vg)

        call tcheck(t,tgrnd,n)
        call tcheck(q,qgrnd,n)
        call ucheck(u,v,z,ustar,lmonin,z0m,hemi,psi0,psi1,n)

        test=abs(2.*(ustar-ustar0)/(ustar+ustar0))
        if (test.lt.tol) exit
        call level2(e,u,v,t,lscale,dzh,n)

        do i=1,n-1
          u(i)=w*usave(i)+(1.-w)*u(i)
          v(i)=w*vsave(i)+(1.-w)*v(i)
          t(i)=w*tsave(i)+(1.-w)*t(i)
          q(i)=w*qsave(i)+(1.-w)*q(i)
          e(i)=w*esave(i)+(1.-w)*e(i)
          usave(i)=u(i)
          vsave(i)=v(i)
          tsave(i)=t(i)
          qsave(i)=q(i)
          esave(i)=e(i)
        end do
        ustar0=ustar
      end do

c     iter_count=iter_count+min(iter,itmax)
c     write(96,*) "iter_count in inits =",iter_count, iter,dbl

c     call check1(ustar,1,ilong,jlat,1)

      if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,kq,ke,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif


      return
      end subroutine inits

      subroutine tcheck(t,tgrnd,n)
!@sum   tcheck checks for reasonable temperatures
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
c ----------------------------------------------------------------------
c This routine makes sure that the temperature remains within
c  reasonable bounds during the initialization process. (Sometimes the
c  the computed temperature iterated out in left field someplace,
c  *way* outside any reasonable range.) This routine keeps the temp
c  between the maximum and minimum of the boundary temperatures.
c ----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: n    !@var n array dimension
      real*8, dimension(n),intent(inout) ::  t !@var t temperature array
      real*8, intent(in) :: tgrnd  !@var tgrnd ground temperature

      integer, dimension(20) :: imax,imin
      integer nmin,nmax
      real*8 tmin,tmax,dt
      integer i                 !@var i  loop and dummy variables

      if (tgrnd.lt.t(n)) then
        tmin=tgrnd
        tmax=t(n)
      else
        tmin=t(n)
        tmax=tgrnd
      endif
      nmin=0
      nmax=0

      do i=1,n-1
        if (t(i).lt.tmin) then
          nmin=nmin+1
          imin(nmin)=i
        endif
        if (t(i).gt.tmax) then
          nmax=nmax+1
          imax(nmax)=i
        endif
      end do

      dt=0.01*(tmax-tmin)
      if (nmin.gt.0) then
        do i=1,nmin
          t(imin(i))=tmin+float(i)*dt
       end do
      endif

      if (nmax.gt.0) then
        do i=1,nmax
          t(imax(i))=tmax-float(i)*dt
       end do
      endif

      return
      end subroutine tcheck

      subroutine ucheck(u,v,z,ustar,lmonin,z0m,hemi,psi0,psi1,n)
!@sum  ucheck makes sure that the winds remain within reasonable
!@+    bounds during the initialization process. (Sometimes the computed
!@+    wind speed iterated out in left field someplace, *way* outside
!@+    any reasonable range.) Tests and corrects both direction and
!@+    magnitude of the wind rotation with altitude. Tests the total
!@+    wind speed via comparison to similarity theory. Note that it
!@+    works from the top down so that it can assume that at level (i),
!@+    level (i+1) displays reasonable behavior.
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
      implicit none
      real*8, parameter :: psistb=15.*radian, psiuns=5.*radian

      integer, intent(in) :: n  !@var n array dimension
      real*8, dimension(n),intent(in) :: z
      real*8, dimension(n),intent(inout) :: u,v

      real*8, intent(in) :: lmonin,z0m,hemi,psi0,psi1,ustar

      real*8 psilim,psirot,angle,utotal,x,x0,dpsim,psiu,zbyl,z0byl,utest
      integer i                 !@var i  loop and dummy variables

      if (lmonin.ge.0.) then
        psilim=psistb
        else
        psilim=psiuns
      endif

c First, check the rotation of the wind vector:

      do i=n-1,1,-1
        psiu=atan2(v(i),u(i)+teeny)
        if(psiu.lt.0.) psiu=psiu+2.*pi
        psirot=psiu-psi0
c --------------------------------------------------------------------
c Next, check the magnitude of the wind. If the gradient is incorrect,
c  set the wind magnitude to that given by similarity theory:

        utotal=sqrt(u(i)*u(i)+v(i)*v(i))
        utest =sqrt(u(i+1)*u(i+1)+v(i+1)*v(i+1))
        if (utotal.gt.utest) then
          zbyl=z(i)/lmonin
          z0byl=z0m/lmonin
          if (lmonin.gt.0.) then
            dpsim=-gamams*(zbyl-z0byl)
            else
            x  = (1.-gamamu* zbyl)**0.25d0
            x0 = (1.-gamamu*z0byl)**0.25d0
            dpsim=log((1.+x )*(1.+x )*(1.+x *x )/
     2               ((1.+x0)*(1.+x0)*(1.+x0*x0)))-
     3            2.*(atan(x)-atan(x0))
          endif
          utotal=(ustar/kappa)*(log(z(i)/z0m)-dpsim)
          if (utotal.gt.utest) utotal=0.95d0*utest
          u(i)=utotal*cos(psiu)
          v(i)=utotal*sin(psiu)
        endif

        if (hemi*psirot.lt.0.) then
          angle=psi0+psi1*float(n-i)/float(n-1)
          u(i)=utotal*cos(angle)
          v(i)=utotal*sin(angle)
          go to 100
        endif

        if (hemi*psirot.gt.psilim) then
          angle=psi0+hemi*psilim*float(n-i)/float(n-1)
          u(i)=utotal*cos(angle)
          v(i)=utotal*sin(angle)
        endif

 100  end do

      return
      end subroutine ucheck

      subroutine check1(a,n,ilong,jlat,id)
!@sum   check1 checks for NaN'S and INF'S in real 1-D arrays.
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
      implicit none

      integer, intent(in) :: n      !@var n  array dimension
      integer, intent(in) :: ilong  !@var ilong  longitude identifier
      integer, intent(in) :: jlat   !@var jlat  latitude identifier
      integer, intent(in) :: id     !@var n  integer id
      real*8, dimension(n), intent(in) ::  a !@var a  real array

      integer i,k !@var i,k loop and dummy variables
      character*16 str  !@var str  output string

      do i=1,n
        write(str,'(e16.8)')a(i)
        k=index(str,'N')+index(str,'n')
        if (k.ne.0) then
          write (99,1000) ilong,jlat,i,id,a(i)
          if (id.lt.100) call stop_model('check1',255)
        endif
      end do

      return
 1000 format (1x,'check1:  ilong = ',6x,i3,/,1x,
     2            '         jlat  = ',6x,i3,/,1x,
     3            '         i     = ',6x,i3,/,1x,
     4            '         id    = ',6x,i3,/,1x,
     5            '         value = ',1pe11.4,/)
      end subroutine check1

      subroutine output(u,v,t,q,e,lscale,z,zhat,dzh,
     2                  km,kh,kq,ke,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3                  ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4                  utop,vtop,ttop,qtop,
     5                  dtime,bgrid,ilong,jlat,iter,itype,n)
!@sum   output produces output for diagnostic purposes
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@calls simil
      implicit none
      real*8, parameter :: degree=1./radian

      integer, intent(in) :: n,itype,iter,jlat,ilong
      real*8,  intent(in) :: lscale(n-1),lmonin

      real*8, dimension(n), intent(in) :: u,v,t,q,z
      real*8, dimension(n-1), intent(in) :: e,zhat,dzh,km,kh,kq,ke,gm,gh

      real*8, intent(in) :: cm,ch,cq,z0m,z0h,z0q,ustar,tstar,qstar,tgrnd
     *     ,qgrnd,utop,vtop,ttop,qtop,dtime,bgrid

      integer i  !@var i loop variable

      real*8 :: psitop,psi,utotal,utot1,utot2,sign,shear,tgrad,qgrad
     *     ,egrad,uflux,hflux,qflux,bvfrq2,shear2,rich,dqdz,dtdz,dudz
     *     ,phim,phih,dudzs,dtdzs,dqdzs,uratio,tratio,qratio,dbydzh
     *     ,tgradl,prod,utest,ttest,qtest

      write (99,5000) ilong,jlat,itype,dtime,iter,ustar,tstar,qstar,
     2                lmonin,tgrnd,qgrnd,cm,ch,cq,z0m,z0h,z0q,
     3                utop,vtop,ttop,qtop,
     4                bgrid
      write (99,1000)
      psitop=atan2(v(n),u(n)+teeny)*degree
      if (psitop.lt.0.) psitop=psitop+360.
      do i=1,n
        psi=atan2(v(i),u(i)+teeny)*degree
        if (psi.lt.0.) psi=psi+360.
        psi=psi-psitop
        utotal=sqrt(u(i)*u(i)+v(i)*v(i))
        call simil(utest,ttest,qtest,z(i),ustar,tstar,qstar,
     2             z0m,z0h,z0q,lmonin,tgrnd,qgrnd)
        write (99,3000) i,z(i),u(i),v(i),psi,utotal,utest,
     2                    t(i),ttest,q(i),qtest
      end do
      write (99,9000)
      write (99,2000)
      do i=1,n-1
        utot1=u(i+1)*u(i+1)+v(i+1)*v(i+1)
        utot2=u(i)  *u(i)  +v(i)  *v(i)
        if (utot1.ge.utot2) then
          sign=1.
          else
          sign=-1.
        endif
        shear =sqrt((u(i+1)-u(i))*(u(i+1)-u(i))
     2             +(v(i+1)-v(i))*(v(i+1)-v(i)))/dzh(i)
        tgrad =(t(i+1)-t(i))/dzh(i)
        qgrad =(q(i+1)-q(i))/dzh(i)
        egrad =(e(i+1)-e(i))/dzh(i)
        uflux=km(i)*shear*sign
        hflux=kh(i)*tgrad
        qflux=kh(i)*qgrad
        bvfrq2=grav*log(t(i+1)/t(i))/dzh(i)+1d-12
        shear2=((u(i+1)-u(i))/dzh(i))**2+
     2         ((v(i+1)-v(i))/dzh(i))**2+1d-8
        rich=bvfrq2/shear2
        write (99,3000) i,zhat(i),e(i),lscale(i),km(i),kh(i),
     2                    gm(i),gh(i),uflux,hflux,qflux,rich
      end do
      write (99,9000)
      write (99,7000)
      do i=1,n-1
        utot1=sqrt(u(i+1)*u(i+1)+v(i+1)*v(i+1))
        utot2=sqrt(  u(i)*  u(i)+  v(i)*  v(i))
        dudz=(utot1-utot2)/dzh(i)
        dtdz=(t(i+1)-t(i))/dzh(i)
        dqdz=(q(i+1)-q(i))/dzh(i)
        if (lmonin.lt.0.) then
          phim = 1./((1.-gamamu*zhat(i)/lmonin)**0.25)
          phih = sigma/sqrt(1.-gamahu*zhat(i)/lmonin)
          else
          phim = 1.+gamams*zhat(i)/lmonin
          phih = sigma*(1.+gamahs*zhat(i)/lmonin)
        endif
        dudzs=ustar*phim/(kappa*zhat(i))
        dtdzs=tstar*phih/(kappa*zhat(i))
        dqdzs=qstar*phih/(kappa*zhat(i))
        uratio=dudzs/(dudz+teeny)
        tratio=dtdzs/(dtdz+teeny)
        qratio=dqdzs/(dqdz+teeny)

        dbydzh=1./dzh(i)
        shear2=dbydzh*dbydzh*((u(i+1)-u(i))**2+
     2                        (v(i+1)-v(i))**2)
        tgradl=dbydzh*log(t(i+1)/t(i))
        prod=km(i)*shear2-grav*kh(i)*tgradl

        write (99,3000) i,zhat(i),dudz,dudzs,dtdz,dtdzs,dqdz,dqdzs,
     2                            uratio,tratio,qratio,prod
      end do
      write (99,9000)
      write (99,9000)
c ----------------------------------------------------------------------
      return
1000  format (1x,'   i','      z     ','      U     ','      V     ',
     2                   '     psi    ','    Utotal  ','    Utest   ',
     2                   '      T     ','    Ttest   ','      Q     ',
     2                   '    Qtest   ',/)
2000  format (1x,'   i','     zhat   ','      E     ','      L     ',
     2                   '      Km    ','      Kh    ',
     3                   '      Gm    ','      Gh    ','    uflux   ',
     4                   '    hflux   ','    qflux   ','    rich #  ',/)
3000  format (1x,i4,11(1x,1pe11.4))
5000  format (1x,'i      = ',9x,i2,/,1x,
     2            'j      = ',9x,i2,/,1x,
     2            'itype  = ',9x,i2,/,1x,
     2            'dtime  = ',1pe11.4,/,1x,
     3            'iter   = ',7x,i4,/,1x,
     4            'ustar  = ',1pe11.4,/,1x,
     5            'tstar  = ',1pe11.4,/,1x,
     5            'qstar  = ',1pe11.4,/,1x,
     5            'lmonin = ',1pe11.4,/,1x,
     5            'tgrnd  = ',1pe11.4,/,1x,
     5            'qgrnd  = ',1pe11.4,/,1x,
     6            'cm     = ',1pe11.4,/,1x,
     6            'ch     = ',1pe11.4,/,1x,
     6            'cq     = ',1pe11.4,/,1x,
     6            'z0m    = ',1pe11.4,/,1x,
     7            'z0h    = ',1pe11.4,/,1x,
     8            'z0q    = ',1pe11.4,/,1x,
     6            'utop   = ',1pe11.4,/,1x,
     6            'vtop   = ',1pe11.4,/,1x,
     7            'ttop   = ',1pe11.4,/,1x,
     8            'qtop   = ',1pe11.4,/,1x,
     8            'bgrid  = ',1pe11.4,/)
7000  format (1x,'   i','     zhat   ','     dudz   ','   dudz sim ',
     2                   '     dtdz   ','   dtdz sim ','     dqdz   ',
     3                   '   dqdz sim ','    uratio  ','    tratio  ',
     4                   '    qratio  ','  production',/)
9000  format (1x)
      end subroutine output

      END MODULE SOCPBL

      subroutine fgrid2(z,f,df)
!@sum  fgrid2 computes functional relationship of z and xi + derivative
!@+    fgrid2 will be called in function rtsafe(fgrid2,x1,x2,xacc)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
      implicit none
      real*8, intent(in) :: z
      real*8, intent(out) :: f,df
      real*8 z1,zn,bgrid,xi,lznbyz1
      common /grids_99/z1,zn,bgrid,xi,lznbyz1
!$OMP  THREADPRIVATE(/GRIDS_99/)

      f=z+bgrid*((zn-z1)*log(z/z1)-(z-z1)*lznbyz1)-xi
      df=1.+bgrid*((zn-z1)/z-lznbyz1)

      return
      end subroutine fgrid2

      function rtsafe(funcd,x1,x2,xacc,ierr)
!@sum   rtsafe use Newton-Rapheson + safeguards to solve F(x)=0
!@auth  Numerical Recipes
!@ver   1.0
      integer,parameter :: maxit=100
      real*8, intent(in) :: x1,x2,xacc
      integer, intent(out) :: ierr
      real*8 rtsafe
      external funcd
      integer j
      real*8 df,dx,dxold,f,fh,fl,temp,xh,xl

      ierr = 0
      call funcd(x1,fl,df)
      call funcd(x2,fh,df)
      if(fl*fh.gt.0.) then
        ierr = 1
        print*, 'Error: root must be bracketed in rtsafe'
      end if
      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df)
      do j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or. abs(2.*
     *       f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return
        call funcd(rtsafe,f,df)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
      end do
      ierr = 1
      print*, 'Error: rtsafe exceeding maximum iterations'
      return
      END

C**** Functions and subroutines for sub-gridscale wind distribution calc

      subroutine sig(tke,mdf,dbl,delt,ch,ws,tsv,wtke,wd,wm)
!@sum calculate sub grid scale velocities
!@auth Reha Cakmur/Gavin Schmidt
      use constant, only : by3,grav
      implicit none
!@var tke local turbulent kinetic energy at surface (J/kg)
!@var mdf downdraft mass flux from moist convection (m/s)
!@var dbl boundary layer height (m)
!@var delt difference between surface and ground T (ts-tg) (K)
!@var tsv virtual surface temperature (K)
!@var ch local heat exchange coefficient
!@var ws grid box mean wind speed (m/s)
      real*8, intent(in):: tke,mdf,dbl,delt,tsv,ch,ws
      real*8, intent(out):: wtke,wd,wm

C**** TKE contribution  ~ sqrt( 2/3 * tke)
      wtke=sqrt(2d0*tke*by3)

C**** dry convection/turbulence contribution ~(Qsens*g*H/rho*cp*T)^(1/3)
      if (delt.lt.0d0) then
        wd=(-delt*ch*ws*grav*dbl/tsv)**by3
      else
        wd=0.d0
      endif
C**** moist convection contribution beta/frac_conv = 200.
      wm=200.d0*mdf
      return
      end SUBROUTINE sig

      subroutine integrate_sgswind(sig,wt,wmin,wmax,ws,icase,wint)
!@sum Integrate sgswind distribution for different cases
!@auth Reha Cakmur/Gavin Schmidt
      use constant, only : by3
      implicit none
!@var sig distribution parameter (m/s)
!@var wt threshold velocity (m/s)
!@var ws resolved wind speed (m/s)
!@var wmin,wmax min and max wind speed for integral
      real*8, intent(in):: sig,wt,ws,wmin,wmax
      real*8, intent(out) :: wint
      integer, intent(in) :: icase
      integer, parameter :: nstep = 100
      integer i
      real*8 x,sgsw,bysig2

C**** integrate distribution from wmin to wmax using Simpsons' rule
C**** depending on icase, integral is done on w, w^2 or w^3
C**** Integration maybe more efficient with a log transform (putting
C**** more points near wmin), but that's up to you...
C**** Use approximate value for small sig and unresolved delta function
      if ((wmax-wmin).lt.sig*dble(nstep)) then
        bysig2=1./(sig*sig)
        wint=0
        do i=1,nstep+1
          x=wmin+(wmax-wmin)*(i-1)/dble(nstep)
          if (i.eq.1.or.i.eq.nstep+1) then
            wint=wint+sgsw(x,ws,wt,bysig2,icase)*by3
          elseif (mod(i,2).eq.0) then
            wint=wint+4d0*sgsw(x,ws,wt,bysig2,icase)*by3
          else
            wint=wint+2d0*sgsw(x,ws,wt,bysig2,icase)*by3
          end if
        end do
        wint=wint*(wmax-wmin)/dble(nstep)
      else                      ! approximate delta function
        if (ws.ge.wmin.and.ws.le.wmax) then
          wint = ws**(icase-1)*(ws-wt)
        else
          wint = 0.
        end if
      end if

      return
      end

      real*8 function sgsw(x,ws,wt,bysig2,icase)
!@sum sgsw function to be integrated for sgs wind calc
!@auth Reha Cakmur/Gavin Schmidt
      use constant, only : pi
      implicit none
      real*8, intent(in) :: x,ws,wt,bysig2
      integer, intent(in) :: icase
      real*8 :: besf,exx,bessi0

      besf=x*ws*bysig2
      if (besf.lt.200d0 ) then
        exx=exp(-0.5d0*(x*x+ws*ws)*bysig2)
        sgsw=(x)**(icase)*(x-wt)*bysig2*BESSI0(besf)*exx
      else ! use bessel function expansion for large besf
        sgsw=(x)**(icase)*(x-wt)*bysig2*exp(-0.5d0*(x-ws)**2*bysig2)
     *       /sqrt(besf*2*pi)
      end if

      return
      end

      real*8 function bessi0(x)
!@sum Bessel's function I0
!@auth  Numerical Recipes
      implicit none
      real*8 ax,y
      real*8, parameter:: p1=1.0d0, p2=3.5156229d0, p3=3.0899424d0,
     *     p4=1.2067492d0, p5=0.2659732d0, p6=0.360768d-1,
     *     p7=0.45813d-2
      real*8, parameter:: q1=0.39894228d0,q2=0.1328592d-1,
     *     q3=0.225319d-2, q4=-0.157565d-2, q5=0.916281d-2,
     *     q6=-0.2057706d-1, q7=0.2635537d-1, q8=-0.1647633d-1,
     *     q9=0.392377d-2
      real*8, intent(in)::x

      if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2.d0
        bessi0=(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=3.75d0/ax
        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     *       +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      end if
      end function bessi0

      SUBROUTINE get_wspdf(wsubtke,wsubwd,wsubwm,wsgcm,wspdf)
!@sum calculates mean surface wind speed using the integral over the
!@sum probability density function of the wind speed from lookup table
!@auth Reha Cakmur/Jan Perlwitz

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)

      USE tracers_dust,ONLY : McFrac,kim,kjm,table1,x11,x21

      IMPLICIT NONE

c Input:
!@var wsubtke velocity scale of sub grid scale turbulence
!@var wsubwd velocity scale of sub grid scale dry convection
!@var wsubwm velocity scale of sub grid scale moist convection
!@var wsgcm GCM surface wind

      REAL*8,INTENT(IN) :: wsubtke,wsubwd,wsubwm,wsgcm

c Output:
!@var wspdf mean surface wind speed from integration over PDF

      REAL*8,INTENT(OUT) :: wspdf

!@var sigma standard deviation of sub grid fluctuations
      REAL*8 :: sigma,ans,dy

      REAL*8 :: work_wspdf1,work_wspdf2,wsgcm1
      CHARACTER(14) :: fname='WARNING_in_PBL'
      CHARACTER(9) :: subr='get_wspdf'
      CHARACTER(5) :: vname1='wsgcm',vname2='sigma'

      wsgcm1=wsgcm

c     This is the case when wsgcm is very small and we set it
c     equal to one of the smallest values in the table index

      IF (wsgcm1 < 0.0005) wsgcm1=0.0005D0

c     If sigma <= 0.0005:

      wspdf=wsgcm1

      CALL check_upper_limit(wsgcm1,x11(Kim),fname,subr,vname1)
c     There is no moist convection, sigma is composed of TKE and DRY
c     convective velocity scale
      IF (wsubwm == 0.) THEN
        sigma=wsubtke+wsubwd
        IF (sigma > 0.0005) THEN
          CALL check_upper_limit(sigma,x21(Kjm),fname,subr,vname2)
c     Linear Polynomial fit (Default)
          CALL polint2dlin(x11,x21,table1,kim,kjm,wsgcm1,sigma,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c          CALL polint2dcub(x11,x21,table1,kim,kjm,wsgcm1,sigma,ans,dy)
          wspdf=exp(ans)
        END IF
      ELSE

c     When there is moist convection, the sigma is the combination of
c     all three subgrid scale parameters (i.e. independent or dependent)
c     Takes into account that the moist convective velocity scale acts
c     only over 5% (mcfrac) of the area.

        work_wspdf1=0.D0
        sigma=wsubtke+wsubwd+wsubwm
        IF (sigma > 0.0005) THEN
          CALL check_upper_limit(sigma,x21(Kjm),fname,subr,vname2)
c     Linear Polynomial fit (Default)
          CALL polint2dlin(x11,x21,table1,kim,kjm,wsgcm1,sigma,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c          CALL polint2dcub(x11,x21,table1,kim,kjm,wsgcm1,sigma,ans,dy)
          work_wspdf1=exp(ans)*mcfrac
        END IF

        work_wspdf2=0.D0
        sigma=wsubtke+wsubwd
        IF (sigma > 0.0005) THEN
          CALL check_upper_limit(sigma,x21(Kjm),fname,subr,vname2)
c     Linear Polynomial fit (Default)
          CALL polint2dlin(x11,x21,table1,kim,kjm,wsgcm1,sigma,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c          CALL polint2dcub(x11,x21,table1,kim,kjm,wsgcm1,sigma,ans,dy)
          work_wspdf2=exp(ans)*(1.d0-mcfrac)
        END IF
        wspdf=work_wspdf1+work_wspdf2

      END IF
#endif

      RETURN
      END SUBROUTINE get_wspdf

      SUBROUTINE polint2dlin(x1a,x2a,ya,m,n,x1,x2,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: m,n
      REAL*8,INTENT(IN) :: x1a(m),x2a(n),ya(m,n),x1,x2
      REAL*8,INTENT(OUT) :: dy,y
      INTEGER, PARAMETER :: nmax=2
      INTEGER j,k,jjj,iii,xx,yy
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax)

      call locatepbl(x1a,m,x1,xx)
      call locatepbl(x2a,n,x2,yy)

      do k=1,nmax
        iii=k
        do j=1,nmax
          jjj=j
          yntmp(j)=ya(xx+jjj-1,yy+iii-1)
          x11(j)=x1a(xx+jjj-1)
        enddo
        if (yntmp(1).eq.-1000.) then
          ymtmp(k)=-1000.
          x22(k)=x2a(yy+iii-1)
        else
          call polintpbl(x11,yntmp,nmax,x1,ymtmp(k),dy)
          x22(k)=x2a(yy+iii-1)
        endif
      enddo
      if (ymtmp(1).eq.-1000.) then
        y=-1000.
      else
        call polintpbl(x22,ymtmp,nmax,x2,y,dy)
      endif
#endif

      return
      END SUBROUTINE polint2dlin
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint2dcub(x1a,x2a,ya,m,n,x1,x2,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: m,n
      REAL*8,INTENT(IN) :: x1a(m),x2a(n),ya(m,n),x1,x2
      REAL*8,INTENT(OUT) :: dy,y
      INTEGER, PARAMETER :: nmax=4
      INTEGER j,k,jjj,iii,xx,yy
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax)

      call locatepbl(x1a,m,x1,xx)
      call locatepbl(x2a,n,x2,yy)

      do k=1,nmax
        iii=k-1
        do j=1,nmax
          jjj=j-1
          yntmp(j)=ya(xx+jjj-1,yy+iii-1)
          x11(j)=x1a(xx+jjj-1)
        enddo
        if (yntmp(1).eq.-1000..or.yntmp(2).eq.-1000.) then
          ymtmp(k)=-1000.
          x22(k)=x2a(yy+iii-1)
        else
          call polintpbl(x11,yntmp,nmax,x1,ymtmp(k),dy)
          x22(k)=x2a(yy+iii-1)
        endif
      enddo
      if (ymtmp(1).eq.-1000..or.ymtmp(2).eq.-1000.) then
        y=-1000.
      else
        call polintpbl(x22,ymtmp,nmax,x2,y,dy)
      endif
#endif

      return
      END SUBROUTINE polint2dcub
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polintpbl(xa,ya,n,x,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(OUT) :: y,dy
      REAL*8, INTENT(IN) :: x,xa(n),ya(n)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(n),d(n)

      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
 11   continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)CALL stop_model('failure in polint in pbl',255)
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
 12     continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
 13   continue
#endif

      return
      END SUBROUTINE polintpbl
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE locatepbl(xx,n,x,j)
!@sum locates parameters of integration in lookup table

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      implicit none
      INTEGER, INTENT(IN):: n
      INTEGER, INTENT(OUT):: j
      REAL*8, INTENT(IN):: x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
#endif

      return
      END SUBROUTINE locatepbl

      SUBROUTINE check_upper_limit(var,varm,cwarn,subr,vname)
!@sum checks whether variable var exceeds maximum value varm. If it does
!@+   set it to maximum value and writes warning message to file cwarn
!@auth Jan Perlwitz
!@ver 1.0

      USE filemanager,ONLY : openunit,closeunit
      USE model_com,ONLY : itime

      IMPLICIT NONE

!@var var Variable to be checked and set to maximum value if larger (inout)
!@var varm Maximum value of variable to be checked (in)
!@var cwarn File name for warning message (in)
!@var subr Name of subroutine (in)
!@var vname Name of variable (in)
      REAL(KIND=8),INTENT(IN) :: varm
      CHARACTER(*),INTENT(IN) :: cwarn,subr,vname

      REAL(KIND=8),INTENT(INOUT) :: var

      INTEGER :: iu

      IF (var > varm) THEN
        CALL openunit(cwarn,iu)
        DO
          READ(iu,*,END=1)
          CYCLE
 1        EXIT
        END DO
        WRITE(iu,*) 'Warning in ',subr,':'
        WRITE(iu,*) ' Variable ',vname,' exceeds allowed maximum value.'
        WRITE(iu,*) ' ',vname,' set to maximum value.'
        WRITE(iu,*) ' itime, max of ',vname,', ',vname,':',itime,varm
     &       ,var
        CALL closeunit(iu)
        var=varm        
      END IF

      RETURN
      END SUBROUTINE check_upper_limit

      real*8 function deltaSST(Qnet,Qsol,ustar_oc)
!@sum deltaSST calculate skin-bulk SST difference (deg C)
!@var Qnet Net heat flux (not including solar, +ve dwn) (W/m2)
!@var Qsol Solar heat flux (W/m2)
!@var ustar_oc ocean u* (m/s)
      USE CONSTANT, only : rhows,visc_wtr_kin
      IMPLICIT NONE
      real*8, intent(in) :: Qnet,Qsol,ustar_oc
C**** k thermal cond. 0.596 W/mK (35 psu, 20 deg, 0 press)
      real*8, parameter :: byk= 1.677d0  ! 1/thermal conductivity (K/(W/m))
!@var lam coefficient (non-dim)
      real*8, parameter :: lam=2.4d0  ! Ward (2007)
!@var del micro-layer thickness (m)
      real*8 :: del
!@var fc fraction of solar absorbed in micro-layer
      real*8 :: fc

C**** calculate micro-layer thickness (m)
C**** Cap ustar so that it doesn't get too small in low wind conditions
C**** ustar_oc > 0.00098 corresponding to tau > 0.001 N/m2
      del = lam*visc_wtr_kin/max(ustar_oc,0.00098d0)

C**** fraction of solar absorbed (Fairall et al, 1996)
      if (del .lt. 1d-8) then
        fc = 0.0545d0 + 11.*del
      else
        fc = 0.137d0 + 11.*del - 6.6d-5*(1.-exp(-del*1250.))/del
      end if
      fc = min(0.24d0,fc)    ! assume max feasible del=1cm?

C**** delta SST = skin - bulk (+ve for flux going down)
C**** includes occurences of warm skin temperatures (generally < 0.01 C)
      deltaSST=(Qnet + fc*Qsol)*del*byk

      end function deltaSST
