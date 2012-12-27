#include "rundeck_opts.h"

!#define ROUGHL_HACK

      module PBL_DRV
      use SOCPBL, only : t_pbl_args, xdelt
      use SOCPBL, only : alloc_pbl_args, dealloc_pbl_args
      implicit none

      private

      public t_pbl_args, pbl, xdelt
      public alloc_pbl_args, dealloc_pbl_args
      contains

      SUBROUTINE PBL(I,J,ITYPE,PTYPE,pbl_args,atm)
!@sum  PBL calculate pbl profiles for each surface type
!@+        Contains code common for all surfaces
!@auth Greg. Hartke/Ye Cheng
!@var DDMS downdraft mass flux in kg/(m^2 s), (i,j)
!@var TDN1 downdraft temperature in K, (i,j)
!@var QDN1 downdraft humidity in kg/kg, (i,j)

      USE EXCHANGE_TYPES
      USE CONSTANT, only :  rgas,grav,omega2,deltx,teeny,lhe,lhs
      USE GEOM, only : sinlat2d
      USE ATM_COM, only : pk
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      USE CLOUDS_COM, only : ddm1
      USE CLOUDS_COM, only : DDMS,TDN1,QDN1,DDML
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,trdn1
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: trradius, trpdens, tr_mm
#endif
#endif
#ifdef TRACERS_AMP
      use TRACER_COM, only: AMP_MODES_MAP,AMP_NUMB_MAP,ntmAMPi,ntmAMPe
      USE AMP_AEROSOL, only : DIAM, AMP_dens,AMP_TR_MM
      USE AERO_SETUP,  only : CONV_DPAM_TO_DGN
#endif

      use SOCPBL, only : npbl=>n, zgs, advanc
      USE PBLCOM

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type
      type (t_pbl_args) :: pbl_args
      class (atmsrf_xchng_vars) :: atm
      REAL*8 Ts
      real*8 :: qsat ! external

#ifdef TRACERS_ON
      integer nx,n
#ifndef TRACERS_TOMAS
     *     ,nAMP
#endif
#endif

c
      REAL*8 ztop,coriol,rhosrf
      REAL*8 qtop,utop,vtop,ufluxs,vfluxs,tfluxs,qfluxs,psitop,psisrf
      INTEGER k
!@var uocean,vocean ocean/ice velocities for use in drag calulation
!@var evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
!@var ZS1    = height of the first model layer (m)
!@var TGV    = virtual potential temperature of the ground (K)
!@+            (if xdelt=0, TGV is the actual temperature)
!@var TKV    = virtual potential temperature of first model layer (K)
!@+            (if xdelt=0, TKV is the actual temperature)
!@var WS     = magnitude of the surface wind (m/s)
!@var PSI    = angular diff. btw geostrophic and surface winds (rads)
!@var WG     = magnitude of the geostrophic wind (m/s)
!@var HEMI   = 1 for northern hemisphere, -1 for southern hemisphere
!@var TG     = bulk ground temperature (K)
!@var ELHX   = latent heat for saturation humidity (J/kg)
!@var dskin  = skin-bulk SST difference (C)
!@VAR QSOL   = solar heating (W/m2)
      real*8 zs1,psi,hemi
!@var POLE   = .TRUE. if at the north or south pole, .FALSE. otherwise
c      logical pole

!**** the following is output from advance (mostly passed through pbl_args)
!@var US     = x component of surface wind, positive eastward (m/s)
!@var VS     = y component of surface wind, positive northward (m/s)
!@var WSGCM  = magnitude of the GCM surface wind - ocean currents (m/s)
!@var WSPDF  = mean surface wind calculated from PDF of wind speed (m/s)
!@var WS     = magn. of GCM surf wind - ocean curr + buoyancy + gust (m/s)
!@var TSV    = virtual potential temperature of the surface (K)
!@+            (if xdelt=0, TSV is the actual temperature)
!@var QS     = surface value of the specific moisture
!@var DBL    = boundary layer height (m)
!@var KMS    = momentum transport coefficient at ZGS (m**2/s)
!@var KHS    = heat transport coefficient at ZGS (m**2/s)
!@var KHQ    = moist transport coefficient at ZGS (m**2/s)
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
      real*8 :: dbl,kms,kqs,cm,ch,cq,z0m,z0h,z0q,ug,vg,w2_1,mdf
!@var dtdt_gcm temp. tendency from processes other than turbulence (K/s)
      real*8 ::  dpdxr,dpdyr,dpdxr0,dpdyr0,dtdt_gcm
      real*8 ::  mdn  ! ,mup
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
#if defined(TRACERS_ON)
!@var  tr local tracer profile (passive scalars)
      real*8, dimension(npbl,pbl_args%ntx) :: tr
      real*8, dimension(ntm) :: trnradius,trndens,trnmm
      real*8 :: rts,rtsdt
#endif

      pbl_args%psurf = atm%srfp(i,j)
      pbl_args%QSOL = atm%fshort(i,j)*atm%cosz1(i,j) ! solar heating

      !qtop=q(i,j,1)
      qtop = atm%q1(i,j)
      pbl_args%tkv = (atm%temp1(I,J)*(1.+qtop*xdelt))*atm%srfpk(i,j)
      pbl_args%ZS1=.5d-2*RGAS*pbl_args%TKV*atm%AM1(i,j)/atm%p1(i,j)
      pbl_args%hemi = sign(1d0,atm%lat(i,j))

      utop = atm%u1(i,j) !ua(1,i,j)
      vtop = atm%v1(i,j) !va(1,i,j)

      if(itype < 4) then
        pbl_args%evap_max = 1.
        pbl_args%fr_sat = 1.    ! entire surface is saturated
        pbl_args%qg_aver = pbl_args%qg_sat ! QG_AVER=QG_SAT
        if(itype==1) then
          pbl_args%elhx = lhe
        else
          pbl_args%elhx = lhs
        endif
      endif

      if(itype > 2) then
        pbl_args%uocean = 0.
        pbl_args%vocean = 0.
      endif

      pbl_args%trhr0 = atm%flong(i,j)

#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      do nx=1,pbl_args%ntx
        n=pbl_args%ntix(nx)
C**** Calculate first layer tracer concentration
          pbl_args%trtop(nx)=
     &         atm%trm1(n,i,j)*atm%byam1(i,j)
      end do
      call tracer_lower_bc(i,j,itype,pbl_args,atm)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) 
      call dust_emission_prep(i,j,itype,pbl_args)
#endif

#endif


ccc extract data needed in driver from the pbl_args structure
      zs1 = pbl_args%zs1
      hemi = pbl_args%hemi
c      pole = pbl_args%pole

      ! Redelsperger et al. 2000, eqn(13), J. Climate, 13, 402-421
      ! tprime,qprime are the pertubation of t and q due to gustiness

      ! pick up one of the following two expressions for gusti

      ! for down draft:
      mdn=max(DDMS(i,j), -0.07d0)
      pbl_args%gusti=log(1.-600.4d0*mdn-4375.*mdn*mdn)

      ! for up draft:
      ! mup=min(DDMS(i,j), 0.1d0)
      ! pbl_args%gusti=log(1.+386.6d0*mup-1850.*mup*mup)

C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces

      IF (ITYPE.GT.2) THEN
        Z0M=ROUGHL(I,J)           ! 30./(10.**ROUGHL(I,J))
      ENDIF
      ztop=zgs+zs1  ! zs1 is calculated before pbl is called
      IF (pbl_args%TKV.EQ.pbl_args%TGV)
     &     pbl_args%TGV = 1.0001d0*pbl_args%TGV

      dbl = bldep(i,j)
      ug   = ugeo(i,j) !ua(L,i,j)
      vg   = vgeo(i,j) !va(L,i,j)

      coriol=sinlat2d(i,j)*omega2

      upbl(:)=atm%uabl(:,i,j)
      vpbl(:)=atm%vabl(:,i,j)
      tpbl(:)=atm%tabl(:,i,j)
      qpbl(:)=atm%qabl(:,i,j)
      epbl(1:npbl-1)=atm%eabl(1:npbl-1,i,j)

#ifdef TRACERS_ON
      do nx=1,pbl_args%ntx
        tr(:,nx)=atm%trabl(:,pbl_args%ntix(nx),i,j)
      end do

      do n = 1,ntm
#ifdef TRACERS_DRYDEP
           trnradius(n) = trradius(n)
           trndens(n)   = trpdens(n)
           trnmm(n)     = tr_mm(n)
#endif
#ifdef TRACERS_AMP
       if (n.ge.ntmAMPi.and.n.le.ntmAMPe) then
         nAMP=n-ntmAMPi+1
        if(AMP_MODES_MAP(nAMP).gt.0) then
         if(DIAM(i,j,1,AMP_MODES_MAP(nAMP)).gt.0.) then
          if(AMP_NUMB_MAP(nAMP).eq. 0) then    ! Mass
        trnradius(n)=0.5*DIAM(i,j,1,AMP_MODES_MAP(nAMP))
          else                              ! Number
        trnradius(n)=0.5*DIAM(i,j,1,AMP_MODES_MAP(nAMP))
     +               *CONV_DPAM_TO_DGN(AMP_MODES_MAP(nAMP))
          endif

           call AMPtrdens(i,j,1,n)
           call AMPtrmass(i,j,1,n)

          trndens(n) =AMP_dens(i,j,1,AMP_MODES_MAP(nAMP))
          trnmm(n)   =AMP_TR_MM(i,j,1,AMP_MODES_MAP(nAMP))
        endif   
        endif   
       endif 
#endif  
      enddo
#endif

      cm=atm%cmgs(i,j)
      ch=atm%chgs(i,j)
      cq=atm%cqgs(i,j)
      dpdxr  = DPDX_BY_RHO(i,j)
      dpdyr  = DPDY_BY_RHO(i,j)
      dpdxr0 = DPDX_BY_RHO_0(i,j)
      dpdyr0 = DPDY_BY_RHO_0(i,j)

      mdf = ddm1(i,j)

!!! put some results from above to pbl_args
      pbl_args%dbl = dbl
      pbl_args%ug = ug
      pbl_args%vg = vg
      pbl_args%wg = sqrt(ug*ug+vg*vg)
      pbl_args%cm = cm
      pbl_args%ch = ch
      pbl_args%cq = cq

#ifdef USE_PBL_E1
      pbl_args%ddml_eq_1=.false.
#else
      pbl_args%ddml_eq_1=DDML(i,j).eq.1
#endif

      ! if ddml_eq_1=.false.,
      ! i.e., either USE_PBL_E1 or DDML(i,j) is not 1,
      ! then tdns,qdns,tprime,qprime are not in use

      if (pbl_args%ddml_eq_1) then
        pbl_args%tdns=TDN1(i,j)*atm%srfpk(i,j)/pk(1,i,j)
        pbl_args%qdns=QDN1(i,j)
#ifdef TRACERS_ON
        do nx=1,pbl_args%ntx
          pbl_args%trdn1(nx)=TRDN1(pbl_args%ntix(nx),i,j)
        end do
#endif
      else
        pbl_args%tdns=0.d0
        pbl_args%qdns=0.d0
#ifdef TRACERS_ON
        pbl_args%trdn1(:)=0.
#endif
      endif

      dtdt_gcm = (pbl_args%tkv - t1_after_aturb(i,j)*
     &     atm%srfpk(i,j))/pbl_args%dtsurf
      call advanc( pbl_args,coriol,utop,vtop,qtop,ztop,mdf
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &     ,dtdt_gcm,u1_after_aturb(i,j),v1_after_aturb(i,j)
     &     ,i,j,itype
     &     ,kms,kqs,z0m,z0h,z0q,w2_1,ufluxs,vfluxs,tfluxs,qfluxs
     &     ,upbl,vpbl,tpbl,qpbl,epbl
#if defined(TRACERS_ON)
     &     ,tr,trnradius,trndens,trnmm
#endif
     &     )

      atm%uabl(:,i,j)=upbl(:)
      atm%vabl(:,i,j)=vpbl(:)
      atm%tabl(:,i,j)=tpbl(:)
      atm%qabl(:,i,j)=qpbl(:)
      atm%eabl(1:npbl-1,i,j)=epbl(1:npbl-1)
      rhosrf=100.*pbl_args%psurf/(rgas*pbl_args%tsv)
#ifdef TRACERS_ON
      do nx=1,pbl_args%ntx
        n = pbl_args%ntix(nx)
        atm%trabl(:,n,i,j)=tr(:,nx)
        atm%travg(n,i,j) = pbl_args%trs(nx)
        atm%travg_byvol(n,i,j) = pbl_args%trs(nx)*rhosrf
      end do
#ifdef TRACERS_DRYDEP
      atm%dep_vel(:,i,j)=pbl_args%dep_vel(:)
      atm%gs_vel(:,i,j)=pbl_args%gs_vel(:)
      do nx=1,pbl_args%ntx
        n = pbl_args%ntix(nx)
        rts=rhosrf*pbl_args%trs(nx)
        rtsdt=rts*pbl_args%dtsurf        ! kg*s/m^3
        atm%drydflx(n,i,j)=-rtsdt*
     &       (atm%dep_vel(n,i,j)+atm%gs_vel(n,i,j)) ! kg/m2
      enddo
#endif
#endif

      atm%cmgs(i,j)=pbl_args%cm
      atm%chgs(i,j)=pbl_args%ch
      atm%cqgs(i,j)=pbl_args%cq
      atm%ipbl(i,j)=1  ! ipbl is used in subroutine init_pbl

      psitop=atan2(vg,ug+teeny)
      psisrf=atan2(pbl_args%vs,pbl_args%us+teeny)
      psi   =psisrf-psitop
      atm%ustar_pbl(i,j)=pbl_args%ustar
C ******************************************************************
      TS=pbl_args%TSV/(1.+pbl_args%QSRF*xdelt)
      if ( ts.lt.152d0 .or. ts.gt.423d0 ) then
        write(6,*) 'PBL: Ts bad at',i,j,' itype',itype,ts
        if (ts.gt.1d3) call stop_model("PBL: Ts out of range",255)
        if (ts.lt.50d0) call stop_model("PBL: Ts out of range",255)
      end if

      atm%tsavg(i,j) = ts
      atm%qsavg(i,j) = pbl_args%qsrf
      atm%usavg(i,j) = pbl_args%us
      atm%vsavg(i,j) = pbl_args%vs
      atm%wsavg(i,j) = pbl_args%ws
      atm%rsavg(i,j) = atm%qsavg(i,j)/
     &     qsat(ts,pbl_args%elhx,atm%srfp(i,j))


      atm%TAUAVG(I,J) = pbl_args%CM*pbl_args%WS*pbl_args%WS*rhosrf
      atm%tgvAVG(I,J) = pbl_args%tgv
      atm%qgAVG(I,J) = pbl_args%qg_aver
      atm%gustiwind(i,j) = pbl_args%gusti
      atm%dblavg(i,j) = dbl
      atm%rhoavg(i,j) = rhosrf
      atm%w2_l1(I,J) = w2_1
      atm%ciaavg(i,j)  =  psi
      atm%khsavg(i,j)  =  pbl_args%khs

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      atm%wsgcm(i,j) = pbl_args%wsgcm
      atm%wspdf(i,j) = pbl_args%wspdf
      atm%wsubwd(i,j) = pbl_args%wsubwd
      atm%wsubtke(i,j) = pbl_args%wsubtke
      atm%wsubwm(i,j) = pbl_args%wsubwm
#endif

ccc put drive output data to pbl_args structure
      pbl_args%psi = psi ! maybe also should be moved to ADVANC
                         ! or completely otside of PBL* ?

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      call PBL_adiurn_dust(I,J,ITYPE,PTYPE,pbl_args,atm)
#endif

#ifdef TRACERS_ON
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) 
      call save_dust_emission_vars(i,j,itype,pbl_args)
#endif
#endif

      RETURN
      END SUBROUTINE PBL

#ifdef TRACERS_ON
      subroutine tracer_lower_bc(i,j,itype,pbl_args,atm)
      use exchange_types
      use geom, only : byaxyp
      use OldTracer_mod, only :
     &     dodrydep,tr_wd_type,nWATER,nPART,nGAS
      implicit none
      integer, intent(in) :: i,j,itype  !@var itype surface type
      type (t_pbl_args) :: pbl_args
      class (atmsrf_xchng_vars) :: atm
c
      integer :: n,nx
c

c todo: try to merge the itype==4 and itype<4 cases
      if(itype.le.3) then ! ice/water surfaces
C 
C**** Set up b.c. for tracer PBL calculation if required
        do nx=1,pbl_args%ntx
          n=pbl_args%ntix(nx)
C**** set defaults
          pbl_args%trsfac(nx)=0.
          pbl_args%trconstflx(nx)=0.
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
#ifdef TRACERS_WATER
          pbl_args%tr_evap_max(nx)=1.
C**** This distinguishes water from gases or particle
          if ( tr_wd_TYPE(n) == nWATER ) then
c          trgrnd(nx)=atmgla%gtracer(n,i,j)
C**** trsfac and trconstflx are multiplied by cq*ws and QG_SAT in PBL
            pbl_args%trsfac(nx)=1.
            pbl_args%trconstflx(nx)=atm%gtracer(n,i,j)

          else if ( tr_wd_TYPE(n) == nGAS .or.
     &           tr_wd_TYPE(n) == nPART ) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
            if(dodrydep(n)) then
              pbl_args%trsfac(nx)=1. !then multiplied by deposition velocity in PBL
#ifdef TRACERS_WATER
              pbl_args%tr_evap_max(nx)=1.d30
c            trgrnd(nx)=0.
#endif
            end if
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
C**** Now send kg/m^2/s to PBL, and divided by rho there.
#ifndef SKIP_TRACER_SRCS
            pbl_args%trconstflx(nx)=atm%trflux_prescr(n,i,j) ! kg/m^2/s
#endif /*SKIP_TRACER_SRCS*/

#ifdef TRACERS_WATER
          endif
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
! transplanted from SURFACE.f:
          IF (pbl_args%ocean) THEN ! OCEAN
            !trgrnd(nx)=atm%gtracer(n,i,j)
            pbl_args%trsfac(nx)=1.
            pbl_args%trconstflx(nx)=atm%gtracer(n,i,j)
          END IF
        !need to redo this here because the previous line has changed trconstflx to zero.
        !because we have no sources. is there a better way to do this?
          !call stop_model('why 1/area in the following line?',255)
          pbl_args%trconstflx(nx)=atm%gtracer(n,i,j) * byaxyp(i,j) !kg,co2/kg,air/m2
#endif

        end do
      else
c LAND SURFACE
        do nx=1,pbl_args%ntx
          n=pbl_args%ntix(nx)
C**** set defaults
          pbl_args%trsfac(nx)=0.
          pbl_args%trconstflx(nx)=0.
#ifdef TRACERS_WATER
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
C**** The select is used to distinguish water from gases or particle
! select removed because of OMP compiler bug
!        select case (tr_wd_TYPE(n))
!        case (nWATER)
          if (tr_wd_TYPE(n) .eq. nWATER) then
C**** no fractionation from ground (yet)
C**** trsfac and trconstflx are multiplied by cq*ws and QG in PBL
            pbl_args%trsfac(nx)=1.
            pbl_args%trconstflx(nx)=atm%gtracer(n,i,j)
!        case (nGAS, nPART)
          elseif(tr_wd_TYPE(n).eq.nGAS .or. tr_wd_TYPE(n).eq.nPART) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
            if(dodrydep(n)) pbl_args%trsfac(nx) = 1.
          !then multiplied by deposition velocity in PBL
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
            pbl_args%trconstflx(nx)=atm%trflux_prescr(n,i,j) ! kg/m^2/s
#ifdef TRACERS_WATER
!        end select
          end if
#endif
        end do

#ifdef TRACERS_WATER
c**** water tracers are also flux limited
        do nx=1,pbl_args%ntx
          n=pbl_args%ntix(nx)
C       pbl_args%tr_evap_max(nx) = evap_max * trsoil_rat(nx)
          pbl_args%tr_evap_max(nx)=pbl_args%evap_max*atm%gtracer(n,i,j)
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) pbl_args%tr_evap_max(nx) = 1.d30
#endif
        end do
#endif
      endif ! itype check

      return
      end subroutine tracer_lower_bc
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) 
      subroutine dust_emission_prep(i,j,itype,pbl_args)
      use constant, only : by3
      use fluxes, only : pprec,pevap
      use tracers_dust, only : hbaij,ricntd
      use clouds_com, only : ddml,fss
      implicit none
      integer, intent(in) :: i,j  !@var i,j grid point
      integer, intent(in) :: itype  !@var itype surface type
      type (t_pbl_args) :: pbl_args
c
      pbl_args%hbaij=hbaij(i,j)
      pbl_args%ricntd=ricntd(i,j)
      pbl_args%pprec=pprec(i,j)
      pbl_args%pevap=pevap(i,j)
c**** fractional area of moist convection * fraction of downdrafts
c**** (=1/3), but only if downdrafts reach the lowest atmospheric
c**** layer. It's only needed to constrain dust emission due to
c**** downdrafts for soil type earth, but it's needed here to calculate
c**** wspdf in PBL.f for the other soil types.
      IF (ddml(i,j) == 1) THEN
        pbl_args%mcfrac=(1.D0-fss(1,i,j))*By3
      ELSE
        pbl_args%mcfrac=0.D0
      END IF
      return
      end subroutine dust_emission_prep

      subroutine save_dust_emission_vars(i,j,itype,pbl_args)
      use tracers_dust, only : hbaij,ricntd
      implicit none
      integer, intent(in) :: i,j  !@var i,j grid point
      integer, intent(in) :: itype  !@var itype surface type
      type (t_pbl_args) :: pbl_args
      hbaij(i,j)=pbl_args%hbaij
      ricntd(i,j)=pbl_args%ricntd
      return
      end subroutine save_dust_emission_vars
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      SUBROUTINE PBL_adiurn_dust(I,J,ITYPE,PTYPE,pbl_args,atm)
      use exchange_types
      USE CONSTANT, only :  rgas,grav,deltx,teeny
      use OldTracer_mod, only : dodrydep,trname
      USE TRACER_COM, only : ntm,ntm_dust
      use SOCPBL, only : npbl=>n
      USE PBLCOM
      USE DIAG_COM, only :
     *      adiurn=>adiurn_loc,ndiupt,ndiuvar,ijdd,adiurn_dust
#ifndef NO_HDIURN
     *     ,hdiurn=>hdiurn_loc
#endif
     *     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf,idd_wtrsh
#ifdef TRACERS_DUST
     *     ,idd_emis,idd_emis2
     *     ,idd_ws2,idd_ustar,idd_us3,idd_stress,idd_lmon
     *     ,idd_rifl,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,idd_tabl1
     *     ,idd_qabl1,idd_zhat1,idd_e1,idd_km1,idd_ri1,idd_grav,idd_turb
#endif

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type
      type (t_pbl_args) :: pbl_args
      class (atmsrf_xchng_vars) :: atm

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER,PARAMETER :: n_idxd=6
#else
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idxd=16+6*npbl+4*(npbl-1)
#endif
#endif

      REAL*8 ts,ws,rhosrf
      integer nx,n
      INTEGER :: kr,ii,idxd(n_idxd)
      REAL*8 :: tmp(NDIUVAR)


      IF (adiurn_dust == 1 .and. pbl_args%moddd==0) THEN
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
        DO KR=1,NDIUPT
          IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
            tmp(:) = 0.
            ii = i
            idxd=(/ idd_wsgcm,idd_wspdf,idd_wtke,idd_wd,idd_wm,idd_wtrsh
#ifdef TRACERS_DUST
     &        ,idd_emis,   idd_emis2,  idd_turb,   idd_grav,   idd_ws2
     &        ,idd_ustar,  idd_us3,    idd_stress, idd_lmon,   idd_rifl,
     &       (idd_zpbl1+ii-1,ii=1,npbl), (idd_uabl1+ii-1,ii=1,npbl),
     *       (idd_vabl1+ii-1,ii=1,npbl), (idd_uvabl1+ii-1,ii=1,npbl),
     *       (idd_tabl1+ii-1,ii=1,npbl), (idd_qabl1+ii-1,ii=1,npbl),
     *       (idd_zhat1+ii-1,ii=1,npbl-1), (idd_e1+ii-1,ii=1,npbl-1),
     *       (idd_km1+ii-1,ii=1,npbl-1), (idd_ri1+ii-1,ii=1,npbl-1)
#endif
     &      /)
            rhosrf=100.*pbl_args%psurf/(rgas*pbl_args%tsv)
            TS=pbl_args%TSV/(1.+pbl_args%QSRF*xdelt)
            ws = pbl_args%ws
            tmp(idd_wsgcm)=pbl_args%wsgcm*ptype
            tmp(idd_wspdf)=pbl_args%wspdf*ptype
            tmp(idd_wtke)=pbl_args%wsubtke*ptype
            tmp(idd_wd)=pbl_args%wsubwd*ptype
            tmp(idd_wm)=pbl_args%wsubwm*ptype
            if(itype==4) then
              tmp(idd_wtrsh)=pbl_args%wtrsh*ptype
#ifdef TRACERS_DUST
              tmp(idd_emis)=0.D0
              tmp(idd_emis2)=0.D0
              do n=1,Ntm_dust
                tmp(idd_emis)=tmp(idd_emis)
     &               +pbl_args%dust_flux(n)*ptype
                tmp(idd_emis2)=tmp(idd_emis2)
     &               +pbl_args%dust_flux2(n)*ptype
              enddo
#endif
            endif
#ifdef TRACERS_DUST
            tmp(idd_turb)=0.D0
            tmp(idd_grav)=0.D0
#ifdef TRACERS_DRYDEP
            do nx=1,pbl_args%ntx
              n = pbl_args%ntix(nx)
              if (dodrydep( n )) then
                select case(trname( n ))
                case('Clay', 'Silt1', 'Silt2', 'Silt3')
                  tmp( idd_turb ) = tmp( idd_turb ) + ptype * rhosrf
     &                 * pbl_args%trs(nx) * pbl_args%dep_vel(n)
                  tmp( idd_grav ) = tmp( idd_grav ) + ptype * rhosrf
     &                 * pbl_args%trs(nx) * pbl_args%gs_vel(n)
                end select
              end if
            end do
#endif
            tmp(idd_ws2)=ws*ws*ptype
            tmp(idd_ustar)=pbl_args%ustar*ptype
            tmp(idd_us3)=ptype*pbl_args%ustar**3
            tmp(idd_stress)=rhosrf*pbl_args%cm*(pbl_args%ws**2)*ptype
            tmp(idd_lmon)=pbl_args%lmonin*ptype
            tmp(idd_rifl)=
     &           +ptype*grav*(ts-pbl_args%tg)*pbl_args%zgs/
     &           (ws*ws*pbl_args%tg)
            tmp(idd_zpbl1:idd_zpbl1+npbl-1)=ptype*pbl_args%z(1:npbl)
            tmp(idd_uabl1:idd_uabl1+npbl-1)=
     *           ptype*atm%uabl(1:npbl,i,j)
            tmp(idd_vabl1:idd_vabl1+npbl-1)=
     *           ptype*atm%vabl(1:npbl,i,j)
            tmp(idd_uvabl1:idd_uvabl1+npbl-1)=ptype*sqrt(
     *           atm%uabl(1:npbl,i,j)*atm%uabl(1:npbl,i,j)+
     *           atm%vabl(1:npbl,i,j)*atm%vabl(1:npbl,i,j))
            tmp(idd_tabl1:idd_tabl1+npbl-1)=
     *           ptype*atm%tabl(1:npbl,i,j)
            tmp(idd_qabl1:idd_qabl1+npbl-1)=
     *           ptype*atm%qabl(1:npbl,i,j)
            tmp(idd_zhat1:idd_zhat1+npbl-2)=ptype
     *           *pbl_args%zhat(1:npbl-1)
            tmp(idd_e1:idd_e1+npbl-2)=atm%eabl(1:npbl-1,i,j)*ptype
            tmp(idd_km1:idd_km1+npbl-2)=ptype*pbl_args%km(1:npbl-1)
            tmp(idd_ri1:idd_ri1+npbl-2)=ptype*pbl_args%gh(1:npbl-1)
     *           /(pbl_args%gm(1:npbl-1)+1d-20)
#endif
            ADIURN(idxd(:),kr,pbl_args%ih)=
     &      ADIURN(idxd(:),kr,pbl_args%ih)  +tmp(idxd(:))
#ifndef NO_HDIURN
            HDIURN(idxd(:),kr,pbl_args%ihm)=
     &      HDIURN(idxd(:),kr,pbl_args%ihm) +tmp(idxd(:))
#endif
          END IF
        END DO
      END IF

      RETURN
      END SUBROUTINE PBL_adiurn_dust
#endif

      end module PBL_DRV
      
      subroutine read_pbl_tsurf_from_nmcfile
      USE FILEMANAGER
      USE DOMAIN_DECOMP_ATM, only : GRID, READT_PARALLEL
      use resolution, only : lm
      use fluxes, only : atmsrf
      implicit none
      integer :: iu_NMC,tsurf_record
      call openunit("AIC",iu_NMC,.true.,.true.)
      tsurf_record = 1+4*lm +1  ! skip over psrf,u,v,t,q
      CALL READT_PARALLEL(grid,iu_NMC,NAMEUNIT(iu_NMC),atmsrf%TSAVG,
     &     tsurf_record)
      call closeunit(iu_NMC)
      atmsrf%tgvavg(:,:) = atmsrf%tsavg(:,:) ! not used for init. set anyway.
      return
      end subroutine read_pbl_tsurf_from_nmcfile

      subroutine init_pbl(inipbl,istart)
c -------------------------------------------------------------
c These routines include the array ipbl which indicates if the
c  computation for a particular ITYPE was done last time step.
c Sets up the initialization of wind, temperature, and moisture
c  fields in the boundary layer. The initial values of these
c  fields are obtained by solving the static equations of the
c  Level 2 model. This is used when starting from a restart
c  file that does not have this data stored.
c -------------------------------------------------------------
      USE FILEMANAGER
      USE Dictionary_mod
      USE CONSTANT, only : lhe,lhs,tf,omega2,deltx
      USE ATM_COM, only : u,v,p,t,q
      USE GEOM, only : imaxj,sinlat2d
#ifdef TRACERS_ON
      use TRACER_COM, only: NTM
#endif
!      USE SOCPBL, only : dpdxr,dpdyr,dpdxr0,dpdyr0

      USE SOCPBL, only : npbl=>n,zgs,inits,XCDpbl,ccoeff0,skin_effect
     &     ,xdelt, maxNTM
      USE GHY_COM, only : fearth
      USE PBLCOM
      USE DOMAIN_DECOMP_ATM, only : GRID, READT_PARALLEL
      USE DOMAIN_DECOMP_1D, only : WRITET_PARALLEL, getDomainBounds
      USE ATM_COM, only : pmid,pk,pedn,pek
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
     &    ,ua=>ualij,va=>valij
      USE SEAICE_COM, only : si_atm
      USE FLUXES, only : atmocn,atmice,atmgla,atmlnd,flice,fland
     &     ,asflx,atmsrf
#ifdef USE_ENT
      use ent_mod, only: ent_get_exports
      use ent_com, only : entcells
#endif


      IMPLICIT NONE
C**** ignore ocean currents for initialisation.
      real*8, parameter :: uocean=0.,vocean=0.
!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
!@var istart what kind of (re)start is being done
      integer, intent(in) :: istart

!@var iu_CDN unit number for roughness length input file
      integer :: iu_CDN
      integer :: ilong  !@var ilong  longitude identifier
      integer :: jlat   !@var jlat  latitude identifier
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,4) ::
     *                                                      tgvdat

      integer :: ipatch,itype4  !@var itype surface type
      integer i,j,k,lpbl !@var i,j,k loop variable
      real*8 pland,pwater,plice,psoil,poice,pocean,
     *     ztop,elhx,coriol,tgrndv,pij,ps,psk,qgrnd
     *     ,utop,vtop,qtop,ttop,zgrnd,cm,ch,cq,ustar
      real*8 qsat
      real*8 ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
      real*8 ug,vg
      real*8, allocatable :: buf(:,:)
      real*8 :: canopy_height, fv
      integer, save :: roughl_from_file = 0

      real*8 :: cdm

      integer :: I_1, I_0, J_1, J_0
      integer :: I_1H, I_0H, J_1H, J_0H

       character*80 :: titrrr
       real*8 rrr(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &            grid%J_STRT_HALO:grid%J_STOP_HALO)


#ifdef TRACERS_ON
       maxNTM = NTM
#endif


        titrrr = "roughness length over land"
        rrr = 0.

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,       J_STOP=J_1)

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      if(istart==2) then ! replace with cold vs warm start logic

        call read_pbl_tsurf_from_nmcfile
        CDM=.001d0

        DO J=J_0,J_1
        DO I=I_0,I_1
#ifndef SCM /* scm set these already */
          atmsrf%usavg(i,j) = ua(1,i,j)
          atmsrf%vsavg(i,j) = va(1,i,j)
          atmsrf%wsavg(i,j) =
     &         sqrt(atmsrf%usavg(i,j)**2 + atmsrf%vsavg(i,j)**2)
#endif
C**** SET SURFACE MOMENTUM TRANSFER TAU0
          atmsrf%TAUAVG(I,J)=1.*CDM*atmsrf%WSAVG(I,J)**2  ! air density = 1 kg/m3
C**** Initialize surface friction velocity
          do ipatch=1,size(asflx)
            asflx(ipatch)%USTAR_pbl(I,J)=atmsrf%WSAVG(I,J)*SQRT(CDM)
          enddo
C**** SET SURFACE SPECIFIC HUMIDITY FROM FIRST LAYER HUMIDITY
          atmsrf%QSAVG(I,J)=Q(I,J,1)
          atmsrf%QGAVG(I,J)=Q(I,J,1)
        ENDDO
        ENDDO
      endif

C things to be done regardless of inipbl

      call sync_param( 'roughl_from_file', roughl_from_file )
!!#if ( ! defined ROUGHL_HACK ) || ( ! defined USE_ENT )
      if ( roughl_from_file .ne. 0 ) then
        allocate ( buf(I_0H:I_1H, J_0H:J_1H) )
        call openunit("CDN",iu_CDN,.TRUE.,.true.)
        CALL READT_PARALLEL(grid,iu_CDN,NAMEUNIT(iu_CDN),buf,1)
        call closeunit(iu_CDN)
        roughl(:,:)=30./(10.**buf(:,:))
        deallocate ( buf )
      endif
!!#endif

      call sync_param( 'XCDpbl', XCDpbl )
      call sync_param( 'skin_effect', skin_effect )

      do j=J_0,J_1
        do i=I_0,I_1
C**** fix roughness length for ocean ice that turned to land ice
          if (si_atm%snowi(i,j).lt.-1.and.flice(i,j).gt.0)
     &         roughl(i,j)=30./(10.**1.84d0)
          if (fland(i,j).gt.0.and.roughl(i,j) .gt. 29.d0) then
            print*,"Roughness length not defined for i,j",i,j
     *           ,roughl(i,j),fland(i,j),flice(i,j)
            print*,"Setting to .01"
            roughl(i,j)=1d-2
          end if
        end do
      end do

      call ccoeff0
      call getztop(zgs,ztop)

      if(.not.inipbl) return

      do j=J_0,J_1
      do i=I_0,I_1
        pland=fland(i,j)
        pwater=1.-pland
        plice=flice(i,j)
        psoil=fearth(i,j)
        poice=si_atm%rsi(i,j)*pwater
        pocean=pwater-poice
        if (pocean.le.0.) then
          tgvdat(i,j,1)=0.
        else
          tgvdat(i,j,1)=atmocn%gtemp(i,j)+TF
        end if
        if (poice.le.0.) then
          tgvdat(i,j,2)=0.
        else
          tgvdat(i,j,2)=atmice%gtemp(i,j)+TF
        end if
        if (plice.le.0.) then
          tgvdat(i,j,3)=0.
        else
          tgvdat(i,j,3)=atmgla%gtemp(i,j)+TF
        end if
        if (psoil.le.0.) then
          tgvdat(i,j,4)=0.
        else
          tgvdat(i,j,4)=atmlnd%gtemp(i,j)+TF
        end if
      end do
      end do

      do ipatch=1,size(asflx)
        itype4 = asflx(ipatch)%itype4
        if ((itype4.eq.1).or.(itype4.eq.4)) then
          elhx=lhe
        else
          elhx=lhs
        endif

        do j=J_0,J_1
          jlat=j
          do i=I_0,imaxj(j)
            coriol=sinlat2d(i,j)*omega2
            tgrndv=tgvdat(i,j,itype4)
            if (tgrndv.eq.0.) then
              asflx(ipatch)%ipbl(i,j)=0
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pedn(1,i,j)    !pij+ptop
            psk=pek(1,i,j)    !expbyk(ps)
            qgrnd=qsat(tgrndv,elhx,ps)

            utop = ua(1,i,j)
            vtop = va(1,i,j)
            qtop=q(i,j,1)
            ttop=t(i,j,1)*(1.+qtop*xdelt)*psk
            t1_after_aturb(i,j) = ttop/psk
            u1_after_aturb(i,j) = utop
            v1_after_aturb(i,j) = vtop

            zgrnd=.1d0 ! formal initialization
            if (itype4.gt.2) zgrnd=roughl(i,j) !         30./(10.**roughl(i,j))

            if (itype4.gt.2) rrr(i,j) = zgrnd

            dpdxr  = DPDX_BY_RHO(i,j)
            dpdyr  = DPDY_BY_RHO(i,j)
            dpdxr0 = DPDX_BY_RHO_0(i,j)
            dpdyr0 = DPDY_BY_RHO_0(i,j)
#ifdef SCM
            utop = u(i,j,1)
            vtop = v(i,j,1)
            ug = utop
            vg = vtop
#endif
            call inits(tgrndv,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,ustar,
     3                 uocean,vocean,ilong,jlat,itype4
     &                 ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &                 ,upbl,vpbl,tpbl,qpbl,epbl,ug,vg)
            asflx(ipatch)%cmgs(i,j)=cm
            asflx(ipatch)%chgs(i,j)=ch
            asflx(ipatch)%cqgs(i,j)=cq

            do lpbl=1,npbl
              asflx(ipatch)%uabl(lpbl,i,j)=upbl(lpbl)
              asflx(ipatch)%vabl(lpbl,i,j)=vpbl(lpbl)
              asflx(ipatch)%tabl(lpbl,i,j)=tpbl(lpbl)
              asflx(ipatch)%qabl(lpbl,i,j)=qpbl(lpbl)
            end do

            do lpbl=1,npbl-1
              asflx(ipatch)%eabl(lpbl,i,j)=epbl(lpbl)
            end do

            asflx(ipatch)%ipbl(i,j)=1
            asflx(ipatch)%ustar_pbl(i,j)=ustar

 200      end do
        end do
      end do

      !write(981) titrrr,rrr
#ifndef CUBED_SPHERE
      call WRITET_PARALLEL(grid,981,"fort.981",rrr,titrrr)
#endif

      return
 1000 format (1x,//,1x,'completed initialization, itype = ',i2,//)
      end subroutine init_pbl

      subroutine loadbl
!@sum loadbl initiallise boundary layer calc each surface time step
!@auth Ye Cheng
c ----------------------------------------------------------------------
c             This routine checks to see if ice has
c              melted or frozen out of a grid box.
c
c For ITYPE=1 (ocean; melted ocean ice since last time step):
c  If there was no computation made for ocean at the last time step,
c  this time step may start from ocean ice result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done. Also deals with newly created lake (from land)
c
c For ITYPE=2 (ocean ice; frozen from ocean since last time step):
c  If there was no computation made for ocean ice at the last time step,
c  this time step may start from ocean result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done.
c
c For ITYPE=3 (land ice; frozen on land since last time step):
c  If there was no computation made for land ice at the last time step,
c  this time step may start from land result. If there was no
c  land ice nor land computation at the last time step, nothing
c  need be done.
c
c For ITYPE=4 (land; melted land ice since last time step):
c  If there was no computation made for land at the last time step,
c  this time step may start from land ice result. If there was no
c  land nor land ice computation at the last time step, nothing
c  need be done. Also deal with newly created earth (from lake)
c
c In the current version of the GCM, there is no need to check the
c  land or land ice components of the grid box for ice formation and
c  melting because pland and plice are fixed. The source code to do
c  this is retained and deleted in the update deck in the event this
c  capability is added in future versions of the model.
c ----------------------------------------------------------------------
      USE EXCHANGE_TYPES
      USE MODEL_COM
      USE GEOM, only : imaxj
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE FLUXES, only : asflx,atmocns,atmices,atmglas,atmlnds
      IMPLICIT NONE
      integer i,j,ip  !@var i,j,ip loop variable
      type(atmsrf_xchng_vars), pointer :: ain,aout

      integer :: J_1, J_0, I_1, I_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

! NOTE ON SUBDIVISIONS OF SURFACE TYPES:
! when initializing profiles for a newly existent surface type,
! values are taken from the _first_ instance of the preexisting
! "donor" surface type (index 1).  This follows the coding before
! the introduction of subdivisions of surface types.  If an instance
! of the surface type already exists in the gridbox, new instances
! of that type could be initialized from it.

      do j=J_0,J_1
        do i=I_0,imaxj(j)

c ******* itype=1: Ocean

          do ip=1,ubound(atmocns,1)
            aout => atmocns(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmices(1)%ipbl(i,j).eq.1) then
                ain => atmices(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              elseif (atmlnds(1)%ipbl(i,j).eq.1) then ! initialise from land
                ain => atmlnds(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

c ******* itype=2: Ocean ice

          do ip=1,ubound(atmices,1)
            aout => atmices(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmocns(1)%ipbl(i,j).eq.1) then
                ain => atmocns(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

c ******* itype=3: Land ice

          do ip=1,ubound(atmglas,1)
            aout => atmglas(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmlnds(1)%ipbl(i,j).eq.1) then
                ain => atmlnds(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

c ******* itype=4: Land

          do ip=1,ubound(atmlnds,1)
            aout => atmlnds(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmglas(1)%ipbl(i,j).eq.1) then
                ain => atmglas(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              elseif (atmocns(1)%ipbl(i,j).eq.1) then
                ain => atmocns(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

C**** initialise some pbl common variables

          do ip=1,size(asflx)
            asflx(ip)%ipbl(i,j) = 0 ! - will be set to 1s when pbl is called
          enddo
        end do
      end do

      return
      end subroutine loadbl

      subroutine setbl(ain,aout,i,j)
!@sum setbl initiallise bl from another surface type for one grid box
!@auth Ye Cheng
      USE EXCHANGE_TYPES
      USE PBLCOM, only : npbl
      IMPLICIT NONE
      type (atmsrf_xchng_vars) :: ain,aout
      integer, INTENT(IN) :: i,j
      integer lpbl  !@var lpbl loop variable
      aout%uabl(:,i,j)=ain%uabl(:,i,j)
      aout%vabl(:,i,j)=ain%vabl(:,i,j)
      aout%tabl(:,i,j)=ain%tabl(:,i,j)
      aout%qabl(:,i,j)=ain%qabl(:,i,j)
      aout%eabl(:,i,j)=ain%eabl(:,i,j)
#ifdef TRACERS_ON
      aout%trabl(:,:,i,j)=ain%trabl(:,:,i,j)
#endif
      aout%cmgs(i,j)=ain%cmgs(i,j)
      aout%chgs(i,j)=ain%chgs(i,j)
      aout%cqgs(i,j)=ain%cqgs(i,j)
      aout%ustar_pbl(i,j)=ain%ustar_pbl(i,j)
      return
      end subroutine setbl

      subroutine getztop(zgs,ztop)
!@sum  getztop computes the value of ztop which is the height in meters
!@+  of the first GCM layer from the surface.
!@+  This subroutine only needs to be called when the BL fields require
!@+  initialization.
!@+  This form for z1 = zgs + zs1 (in terms of GCM parameters) yields an
!@+  average value for zs1. The quantity theta was computed on the
!@+  assumption of zs1=200 m from the original 9-layer model (actually
!@+  was misconstrued as z1 = 200m when it should have been zs1 = 200m)
!@+  and is then applied to all vertical resolutions.
!@auth Greg. Hartke/Ye Cheng
!@var zgs The height of the surface layer.
!@var ztop The height of the top of the BL simulation domain.
!@+   Corresponds to averaged height of the middle of first model layer.

      USE CONSTANT, only : rgas,grav
      USE RESOLUTION, only : psf
      use ATM_COM, only : pednl00 ! use plbot from res file instead
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: ZGS
      REAL*8, INTENT(OUT) :: ZTOP
      real*8, parameter :: theta=269.0727251d0

      ztop=zgs+0.5d0*(pednl00(1)-pednl00(2))*rgas*theta/(grav*psf)

      return
      end subroutine getztop

      subroutine get_dbl
      USE FLUXES, only : atmsrf
      USE CONSTANT, only :  rgas,grav,omega2,deltx,teeny
      USE ATM_COM, only : t,q,ua=>ualij,va=>valij
      USE ATM_COM, only : pmid,pk
      use SOCPBL, only : zgs
      USE PBLCOM
      use GEOM, only : imaxj
      use PBL_DRV
      use domain_decomp_atm, only : grid
      implicit none
      integer :: ldc
      integer :: i,j,l
      real*8 :: zpbl,zpbl1,tbar,tl,pl,tl1,pl1,dbl,ztop
      REAL*8, parameter :: dbl_max=3000., dbl_max_stable=500. ! meters
      real*8 :: thbar ! function

      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)

        ztop = zgs +
     &       .5d-2*RGAS*((atmsrf%temp1(I,J)*(1.+atmsrf%q1(i,j)*xdelt))*
     &       atmsrf%srfpk(i,j))*atmsrf%AM1(i,j)/atmsrf%p1(i,j)

      ! FIND THE PBL HEIGHT IN METERS (DBL) AND THE CORRESPONDING
      ! GCM LAYER (L) AT WHICH TO COMPUTE UG AND VG.
      ! LDC IS THE LAYER TO WHICH DRY CONVECTION/TURBULENCE MIXES

c       IF (TKV.GE.TGV) THEN
c         ! ATMOSPHERE IS STABLE WITH RESPECT TO THE GROUND
c         ! DETERMINE VERTICAL LEVEL CORRESPONDING TO HEIGHT OF PBL:
c         ! WHEN ATMOSPHERE IS STABLE, CAN COMPUTE DBL BUT DO NOT
c         ! KNOW THE INDEX OF THE LAYER.
c         ustar=ustar_pbl(itype,i,j)
c         DBL=min(0.3d0*USTAR/OMEGA2,dbl_max_stable)
c         if (dbl.le.ztop) then
c           dbl=ztop
c           L=1
c         else
c           ! FIND THE VERTICAL LEVEL NEXT HIGHER THAN DBL AND
c           ! COMPUTE Ug and Vg THERE:
c           zpbl=ztop
c           pl1=pmid(1,i,j)         ! pij*sig(1)+ptop
c           tl1=t(i,j,1)*(1.+xdelt*q(i,j,1))*pk(1,i,j)
c           do l=2,ls1
c             pl=pmid(l,i,j)        !pij*sig(l)+ptop
c             tl=t(i,j,l)*(1.+xdelt*q(i,j,l))*pk(l,i,j) !virtual,absolute
c             tbar=thbar(tl1,tl)
c             zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
c             if (zpbl.ge.dbl) exit
c             pl1=pl
c             tl1=tl
c           end do
c         endif

c     ELSE
        ! ATMOSPHERE IS UNSTABLE WITH RESPECT TO THE GROUND
        ! LDC IS THE LEVEL TO WHICH DRYCNV/ATURB MIXES.
        ! FIND DBL FROM LDC.  IF BOUNDARY
        ! LAYER HEIGHT IS LESS THAN DBL_MAX, ASSIGN LDC TO L, OTHERWISE
        ! MUST FIND INDEX FOR NEXT MODEL LAYER ABOVE 3 KM:

        LDC=max(int(DCLEV(I,J)+.5d0),1)
        IF (LDC.EQ.0) LDC=1
        if (ldc.eq.1) then
          dbl=ztop
          l=1
        else
          zpbl=ztop
          pl1=pmid(1,i,j)                             ! pij*sig(1)+ptop
          tl1=t(i,j,1)*(1.+xdelt*q(i,j,1))*pk(1,i,j)  ! expbyk(pl1)
          zpbl1=ztop
          do l=2,ldc
            pl=pmid(l,i,j)                            ! pij*sig(l)+ptop
            tl=t(i,j,l)*(1.+xdelt*q(i,j,l))*pk(l,i,j) ! expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
            if (zpbl.ge.dbl_max) then
              zpbl=zpbl1
              exit
            endif
            pl1=pl
            tl1=tl
            zpbl1=zpbl
          end do
          l=min(l,ldc)
          dbl=zpbl
        endif

c     ENDIF

      ugeo(i,j) = ua(L,i,j)
      vgeo(i,j) = va(L,i,j)
      bldep(i,j) = dbl

      enddo
      enddo
      return
      end subroutine get_dbl

      SUBROUTINE CHECKPBL(SUBR)
!@sum  CHECKPBL Checks whether PBL data are reasonable
!@auth Original Development Team
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE PBLCOM, only : dclev
      USE FLUXES, only : atmsrf
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      integer :: I_1, I_0, J_1, J_0, njpol
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1,
     *               J_STRT=J_0, J_STOP=J_1)
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in boundary layer data
      CALL CHECK3B(atmsrf%wsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'wsavg')
      CALL CHECK3B(atmsrf%tsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'tsavg')
      CALL CHECK3B(atmsrf%qsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'qsavg')
      CALL CHECK3B(dclev(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'dclev')
      CALL CHECK3B(atmsrf%usavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'usavg')
      CALL CHECK3B(atmsrf%vsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'vsavg')
      CALL CHECK3B(atmsrf%tauavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'tauavg')
c      CALL CHECK3C(ustar_pbl(:,I_0:I_1,J_0:J_1),4,I_0,I_1,J_0,J_1,NJPOL,
c     &     SUBR,'ustar')

      CALL CHECK3B(atmsrf%tgvavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'tgvavg')
      CALL CHECK3B(atmsrf%qgavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'qgavg')
      CALL CHECK3B(atmsrf%w2_l1(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'w2_l1')

      END SUBROUTINE CHECKPBL

