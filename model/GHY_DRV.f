#include "rundeck_opts.h"

      module soil_drv
!@sum soil_drv contains variables and routines for the ground
!@+   hydrology driver
!@auth I. Alienov/F. Abramopolous
      use model_com, only : im,jm
      use DOMAIN_DECOMP, only : GRID, GET
      use veg_drv, only : cosday,sinday
      implicit none
      private
      save

      public daily_earth, ground_e, init_gh, earth, conserv_wtg
     $     ,conserv_htg

      !real*8 cosday,sinday
      !real*8 cosdaym1, sindaym1               !nyk TEMPORARY for jday-1
      real*8 adlmass          ! accumulator for dleafmass in daily_earth

      real*8 spgsn !@var spgsn specific gravity of snow
!@dbparam snow_cover_coef coefficient for topography variance in
!@+       snow cover parameterisation for albedo
      real*8 :: snow_cover_coef = .15d0

      contains

      subroutine earth (ns,moddsf,moddd)
!@sum EARTH calculates surface fluxes of sensible heat,
!@+   evaporation, thermal radiation, and momentum drag over earth.
!@auth I. Alienov/F. Abramopolous
c****
      use constant, only : grav,rgas,lhe,lhs
     *     ,sha,tf,rhow,deltx
      use model_com, only : t,p,q,dtsrc,nisurf,dsig,qcheck,jdate
     *     ,jday,jhour,nday,itime,jeq,fearth,modrd,itearth
     *     ,u,v
      use DOMAIN_DECOMP, only : HALO_UPDATE, CHECKSUM, NORTH
      use geom, only : imaxj,dxyp,bydxyp
      use dynamics, only : pmid,pk,pek,pedn,pdsig,am,byam
      use somtq_com, only : mz
      use radncb, only : trhr,fsf, cosz1

      use surf_albedo, only: albvnh   ! added 5/23/03 from RADIATION.f
      !albvnh(9,6,2)=albvnh(sand+8veg,6bands,2hemi) - only need 1st band
      use sle001
     &    , only : advnc,evap_limits,
     &    pr,htpr,prs,htprs,w,ht,snowd,tp,fice,
     &    fv,fb,atrg,ashg,alhg,
     &    abetad,abetav,abetat,
     &    abetap,abetab,abeta,
     &    acna,acnc,agpp,
     &    aevap,aevapw,aevapd,aevapb,
     &    aruns,arunu,aeruns,aerunu,
     &    aepc,aepb,aepp,af0dt,af1dt,zw,tbcs,
     &    qm1,qs,
     &    pres,rho,ts,vsm,ch,srht,trht, !cd,snht,
     &    nsn,dzsn,wsn,hsn,fr_snow
     &     ,ghy_debug
#ifdef TRACERS_WATER
     &     ,tr_w,tr_wsn,trpr,tr_surf,ntg,ntgm,atr_evap,atr_rnff,atr_g
#endif
#ifdef TRACERS_SPECIAL_O18
     &     ,tr_name
#endif
      use veg_drv, only: veg_save_cell,veg_set_cell
      use vegetation, only: update_veg_locals

      use dagcom , only : aij,tsfrez,tdiurn,aj,areg,adiurn,jreg,hdiurn,
     *     ij_rune, ij_arunu, ij_pevap, ij_shdt, ij_beta, ij_trnfp0,
     *     ij_srtr, ij_neth, ij_ws, ij_ts, ij_us, ij_vs, ij_taus,
     *     ij_tauus, ij_tauvs, ij_qs, ij_tg1, ij_evap, j_trhdt, j_shdt,
     *     j_evhdt,j_evap,j_erun,j_run,j_tsrf,j_type,j_tg1,j_tg2,ij_g05
     *     ,ij_g06,ij_g11,ij_g12,ij_g13,ij_g14,ij_g15,ij_g16,ij_g17
     *     ,ij_gpp, ij_dleaf,ij_pblht
     *     ,ij_g18,ij_g19,ij_g20,ij_g21,ij_g22,ij_g23,ij_g24,ij_g25
     *     ,ij_g26,ij_g27,ijdd,idd_ts,idd_tg1,idd_qs,idd_qg,idd_swg
     *     ,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg,idd_wg,idd_us
     *     ,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq,idd_eds,idd_dbl
     *     ,idd_ev,tf_day1,tf_last,ndiupt
#ifdef TRACERS_ON
      use tracer_com, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE,trname
#endif
      use tracer_diag_com, only : taijn,tij_surf
#ifdef TRACERS_WATER
     *     ,tij_evap,tij_grnd,jls_source,tajls,tij_soil
#endif
#ifdef TRACERS_DRYDEP
     *     ,tij_drydep,tij_gsdep,itcon_dd
#endif
      use fluxes, only : trsource,trsrfflx
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer,trprec
#endif
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
      use tracers_DRYDEP, only : dtr_dd
#endif
#endif
      use fluxes, only : dth1,dq1,uflux1,vflux1,e0,e1,evapor,prec,eprec
     *     ,runoe,erunoe,gtemp,precss
      use ghycom, only : ngm,nlsn,
     &     gdeep,wbare,wvege,htbare,htvege,snowbv,
     &     nsn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij,
     *     canopy_temp_ij,snowe,tearth,wearth,aiearth,
     &     evap_max_ij, fr_sat_ij, qg_ij, fr_snow_rad_ij,top_dev_ij
#ifdef TRACERS_WATER
     &     ,tr_wbare,tr_wvege,tr_wsn_ij
#endif
!#ifdef TRACERS_WATER
!     *     ,trvege,trbare,trsnowbv
!#endif
      use vegetation, only :
     &    veg_srht=>srht,veg_pres=>pres,veg_ch=>ch,veg_vsm=>vsm !ia

      USE SOCPBL, only : dtsurf         ! zgs,     ! global
     &     ,zs1,tgv,tkv,qg_sat,qg_aver,hemi,pole     ! rest local
     &     ,us,vs,ws,wsm,wsh,tsv,qsrf,psi,dbl    ! ,edvisc=>kms
     &     ,khs,ug,vg,wg,wint   ! ,kq=>kqs ,ppbl
      use pblcom, only : ipbl,cmgs,chgs,cqgs,tsavg,qsavg
      use pbl_drv, only : pbl, evap_max,fr_sat,uocean,vocean,psurf,trhr0
#ifdef TRACERS_ON
     *     ,trtop,trs,trsfac,trconstflx,ntx,ntix
#ifdef TRACERS_WATER
     *     ,tr_evap_max
#endif
#ifdef TRACERS_DRYDEP
     *     ,dep_vel,gs_vel
#endif
#endif
      use snow_drvm, only : snow_cover_same_as_rad
      implicit none

      integer, intent(in) :: ns,moddsf,moddd
      integer i,j,kr,jr,itype,ih,ihm,ibv
      real*8 shdt,qsats,evap,evhdt,tg2av,ace2av,trhdt,rcdmws,rcdhws
     *     ,cdq,cdm,cdh,elhx,tg,srheat,tg1,ptype,trheat,wtr2av    !,dhgs
     *     ,rhosrf,ma1,tfs,th1,thv1,p1k,psk,ps,pij,psoil,pearth
     *     ,warmer,timez,spring,zs1co,q1

!@var rhosrf0 estimated surface air density
      real*8 rhosrf0

#ifdef TRACERS_ON
      integer n,nx,nsrc
      real*8 totflux(ntm)
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: trsoil_tot,tevapw,tevapd,
     *     tevapb,trruns,trrunu,trsoil_rat 
      real*8, dimension(ntm,0:ngm,2) :: trw
      real*8, dimension(ntm,2) :: trsnowd
      real*8  tdp, tdt1, wsoil_tot, frac, tevap
ccc new vars
      real*8 ntixw(ntm)
#ifdef TRACERS_SPECIAL_O18
      real*8 fracvl
#endif
#endif
#ifdef TRACERS_DRYDEP
      real*8 tdryd, tdd, td1, rtsdt
#endif
#endif
      real*8 qsat
      real*8 srhdt
      real*8 aregij(9,im,jm)
!@var qg rel. humidity at the ground, defined: total_evap = Cq V (qg-qs)
!@var qg_nsat rel. humidity at non-saturated fraction of soil
      real*8 qg, qg_nsat
c****
c**** fearth    soil covered land fraction (1)
c****
c**** snowe     earth snow amount (kg/m**2)
c**** tearth    earth temperature of first layer (c)
c**** wearth    earth water of first layer (kg/m**2)
c**** aiearth   earth ice of first layer (kg/m**2)
c****
c**** wbare  1-6 water of bare soil layer 1-6 (m)
c**** wvege   0  water of vegetation canopy (m)
c****        1-6 water of vegetated soil layer 1-6 (m)
c**** htbare  0  bare soil layer 0 is unused
c****        1-6 heat content of bare soil layer 1-6 (j m-2)
c**** htvege  0  heat content of vegetation canopy (j m-2)
c****        1-6 heat content of vegetated soil layer 1-6 (j m-2)
c**** snowbv  1  snow depth over bare soil (m)
c****         2  snow depth over vegetated soil (m)
c****

C****   define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      dtsurf=dtsrc/nisurf
      zs1co=.5*dsig(1)*rgas/grav

      spring=-1.
      if((jday.ge.32).and.(jday.le.212)) spring=1.
      ih=1+jhour
      ihm = ih+(jdate-1)*24
c****
c**** outside loop over time steps, executed nisurf times every hour
c****
      timez=jday+(mod(itime,nday)+(ns-1.)/nisurf)/nday ! -1 ??
      if(jday.le.31) timez=timez+365.

ccc set i,j - independent stuff for tracers
#ifdef TRACERS_ON
      ntx=0
#ifdef TRACERS_WATER
      ntg=0
#endif
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          ntx=ntx+1
          ntix(ntx) = n
#ifdef TRACERS_WATER
          if ( tr_wd_TYPE(n) == nWATER ) then
            ntg = ntg + 1
            if ( ntg > ntgm ) call stop_model("ghy_drv: ntg > ntgm",255)
            ntixw(ntg) = n
#ifdef TRACERS_SPECIAL_O18
            tr_name(ntg) = trname(n)
#endif
          endif
#endif
        end if
      end do
#endif

      aregij = 0.
c****
c**** outside loop over j and i, executed once for each grid point
c****
C**** halo update u and v for distributed parallelization
       call checksum   (grid, U, __LINE__, __FILE__)
       call halo_update(grid, U, from=NORTH)
       call checksum   (grid, V, __LINE__, __FILE__)
       call halo_update(grid, V, from=NORTH)

!$OMP  PARALLEL DO PRIVATE
!$OMP*  (ACE2AV, ELHX,EVAP,EVHDT, CDM,CDH,CDQ,
!$OMP*   I,ITYPE,ibv, J, KR, MA1,PIJ,PSK,PEARTH,PSOIL,PS,P1K,PTYPE, QG,
!$OMP*   QG_NSAT,QSATS, RHOSRF,RHOSRF0,RCDMWS,RCDHWS, SRHDT,SRHEAT,SHDT,
!$OMP*   TRHEAT, TH1,TFS,THV1,TG1,TG,TRHDT,TG2AV, WARMER,WTR2AV,q1
#if defined(TRACERS_ON)
!$OMP*   ,n,nx,totflux,nsrc
#if defined(TRACERS_WATER)
!$OMP*   ,trsoil_tot,tevapw,tevapd,tevapb,trruns,trrunu,
!$OMP*   trsoil_rat,trw,trsnowd,tdp, tdt1, wsoil_tot,frac,tevap
#endif
#if defined(TRACERS_DRYDEP)
!$OMP*   ,tdryd,tdd,td1,rtsdt
#endif
#endif
!$OMP*   )
!$OMP*   SCHEDULE(DYNAMIC,2)

      loop_j: do j=J_0,J_1
      hemi=1.
      if(j.le.jm/2) hemi=-1.
c**** conditions at the south/north pole
      pole= ( j.eq.1 .or. j.eq.jm )

ccc   if(j.lt.jeq) warmer=-spring
ccc   if(j.ge.jeq) warmer=spring
      if(j.lt.jeq)  then
         warmer=-spring
       else
         warmer=spring
      end if
      loop_i: do i=1,imaxj(j)
c****
c**** determine surface conditions
c****
      pearth=fearth(i,j)
      psoil=pearth
      pij=p(i,j)
      ps=pedn(1,i,j)
      psk=pek(1,i,j)
      p1k=pk(1,i,j)
      th1=t(i,j,1)
      q1=q(i,j,1)
      thv1=th1*(1.+q1*deltx)
      tkv=thv1*psk
      tfs=tf*psoil
      ma1=am(1,i,j)
      qm1=q1*ma1
c     rhosrf=100.*ps/(rgas*tsv)
c     rhosrf=100.*ps/(rgas*tkv)
ccc   jr=jreg(i,j)
c****
c**** earth
c****
      if (pearth.le.0.) then
        ipbl(i,j,4)=0
        cycle loop_i
      endif

#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      do nx=1,ntx
        n = ntix(nx)
C**** Calculate first layer tracer concentration
        trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
      end do
#endif
      itype=4
      ptype=pearth
      pr=prec(i,j)/(dtsrc*rhow)
C**** This variable was originally assoicated with super-saturated
C**** large-scale precip, but has not worked for many moons.
C**** If you want to reinstate it, uncomment this calculation.
c      prs=precss(i,j)/(dtsrc*rhow)
      prs=0.
      htpr=eprec(i,j)/dtsrc
!!! insert htprs here
      w(1:ngm,1) =  wbare(1:ngm,i,j)
      w(0:ngm,2) =  wvege(0:ngm,i,j)
      ht(0:ngm,1) = htbare(0:ngm,i,j)
      ht(0:ngm,2) = htvege(0:ngm,i,j)
      snowd(1:2)  = snowbv(1:2,i,j)
ccc extracting snow variables
      nsn(1:2)          = nsn_ij    (1:2, i, j)
      !isn(1:2)          = isn_ij    (1:2, i, j)
      dzsn(1:nlsn, 1:2) = dzsn_ij   (1:nlsn, 1:2, i, j)
      wsn(1:nlsn, 1:2)  = wsn_ij    (1:nlsn, 1:2, i, j)
      hsn(1:nlsn, 1:2)  = hsn_ij    (1:nlsn, 1:2, i, j)
      fr_snow(1:2)      = fr_snow_ij(1:2, i, j)
ccc tracers variables
#ifdef TRACERS_WATER
      do nx=1,ntg
        n = ntixw(nx)
        ! prognostic vars
        tr_w(nx,0,1) = 0.d0
        tr_w(nx,1:ngm,1) = tr_wbare(n,1:ngm,i,j)
        tr_w(nx,0:ngm,2) = tr_wvege(n,0:ngm,i,j)
        tr_wsn(nx,1:nlsn,1:2) = tr_wsn_ij(n,1:nlsn, 1:2, i, j)
        ! flux in
        trpr(nx)=(trprec(n,i,j)*bydxyp(j))/dtsrc ! kg/m^2 s
        ! concentration of tracers in atm. water at the surface
        if (qm1.gt.0) then  ! avoid very rare error
          tr_surf(nx) = trm(i,j,1,n)*bydxyp(j)*rhow/qm1 ! kg/m^3
        else
          tr_surf(nx) = 0.
        end if
      enddo
#endif
      tg1 = tearth(i,j)
      srheat=fsf(itype,i,j)*cosz1(i,j)
! need this
      srhdt=srheat*dtsurf
c****
c**** boundary layer interaction
c****
      zs1=zs1co*tkv*pij/pmid(1,i,j)
c**** loop over ground time steps
      tg=tg1+tf
      elhx=lhe
      if(tg1.lt.0.)  elhx=lhs
      qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step
      qg = qg_ij(i,j)
      ! if ( qg > 999.d0 ) qg = qg_sat
      qg_aver = qg
      tgv=tg*(1.+qg*deltx)
      psurf=ps
      trhr0 = TRHR(0,I,J)
      rhosrf0=100.*ps/(rgas*tgv) ! estimated surface density
C**** Obviously there are no ocean currents for earth points, but
C**** variables set for consistency with surfce
      uocean=0 ; vocean=0
#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
cddd#ifdef TRACERS_WATER_OLD
cdddC**** Quick and dirty calculation of water tracer amounts in
cdddC**** soils to ensure conservation. Should be replaced with proper
cdddC**** calculation at some point
cdddC**** Calculate mean tracer ratio
cddd      trsoil_tot = 0 ; wsoil_tot = 0
cddd      fb=afb(i,j) ; fv=1.-fb
cddd      do ibv=1,2
cddd        frac=fb
cddd        if (ibv.eq.2) frac=fv
cddd        wsoil_tot=wsoil_tot+snowd(ibv)*frac
cddd        do nx=1,ntx
cddd          n=ntix(nx)
cddd          if (tr_wd_TYPE(n).eq.nWATER) THEN
cddd            trsnowd(nx,ibv) = TRSNOWBV(n,ibv,i,j)
cddd            trsoil_tot(nx)=trsoil_tot(nx)+trsnowd(nx,ibv)*frac
cddd          end if
cddd        end do
cddd        do k= 2-ibv,ngm
cddd          wsoil_tot=wsoil_tot+w(k,ibv)*frac
cddd          do nx=1,ntx
cddd            n=ntix(nx)
cddd            if (tr_wd_TYPE(n).eq.nWATER) THEN
cddd              if (ibv.eq.1) then
cddd                trw(nx,k,ibv)= TRBARE(n,k,i,j)
cddd              else
cddd                trw(nx,k,ibv)= TRVEGE(n,k,i,j)
cddd              end if
cddd              trsoil_tot(nx)=trsoil_tot(nx)+trw(nx,k,ibv)*frac
cddd            end if
cddd          end do
cddd        end do
cddd      end do
cdddC**** calculate new tracer ratio after precip
cddd      wsoil_tot=wsoil_tot+dtsurf*pr
cddd      do nx=1,ntx
cddd        n=ntix(nx)
cddd        if (tr_wd_TYPE(n).eq.nWATER) THEN
cddd          trpr(nx)=(trprec(n,i,j)*bydxyp(j))/dtsrc ! kg/m^2 s
cddd          trsoil_tot(nx)=trsoil_tot(nx)+dtsurf*trpr(nx)
cddd          trsoil_rat(nx)=trsoil_tot(nx)/(rhow*wsoil_tot)
cddd        end if
cddd      end do
cddd#endif
C****
      do nx=1,ntx
        n=ntix(nx)
C**** set defaults
        trsfac(nx)=0.
        totflux(nx)=0.
        trconstflx(nx)=0.
#ifdef TRACERS_WATER
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
C**** The select is used to distinguish water from gases or particle
! select removed because of OMP compiler bug
!        select case (tr_wd_TYPE(n))
!        case (nWATER)
        if (tr_wd_TYPE(n) .eq. nWATER) then
C**** no fractionation from ground (yet)
C**** trsfac and trconstflx are multiplied by cq*wsh in PBL
          trsfac(nx)=1.
          trconstflx(nx)=gtracer(n,itype,i,j)*QG
!        case (nGAS, nPART)
        elseif (tr_wd_TYPE(n).eq.nGAS .or. tr_wd_TYPE(n).eq.nPART) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) trsfac(nx) = 1.
          !then multiplied by deposition velocity in PBL
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
          do nsrc=1,ntsurfsrc(n)
            totflux(nx) = totflux(nx)+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux(nx)*bydxyp(j)   ! kg/m^2/s
#ifdef TRACERS_WATER
!        end select
        end if
#endif
      end do
#endif
c***********************************************************************
c****
ccc actually PBL needs evap (kg/m^2*s) / rho_air
      evap_max = evap_max_ij(i,j) * 1000.d0 / rhosrf0
      fr_sat = fr_sat_ij(i,j)
#ifdef TRACERS_WATER
c**** water tracers are also flux limited
      do nx=1,ntx
        n=ntix(nx)
C       tr_evap_max(nx) = evap_max * trsoil_rat(nx)
        tr_evap_max(nx) = evap_max * gtracer(n,itype,i,j)
#ifdef TRACERS_DRYDEP
        if(dodrydep(n)) tr_evap_max(nx) = 1.d30
#endif
      end do
#endif
      
      call pbl(i,j,itype,ptype)
c****
      cdm = cmgs(i,j,itype)
      cdh = chgs(i,j,itype)
      cdq = cqgs(i,j,itype)
c***********************************************************************
c**** calculate qs
      qs=qsrf
      ts=tsv/(1.+qs*deltx)
c**** calculate rhosrf*cdm*ws
      rhosrf=100.*ps/(rgas*tsv)
      rcdmws=cdm*wsm*rhosrf
      rcdhws=cdh*wsh*rhosrf
c**** calculate fluxes of sensible heat, latent heat, thermal
c****   radiation, and conduction heat (watts/m**2)
c      snht=-sha*rcdhws*(ts-tg)  ! -not used
      trheat=trhr(0,i,j)
c***********************************************************************
c****
c  define extra variables to be passed in surfc:
      pres  =ps
      veg_pres = ps
      rho   =rhosrf
      vsm   =ws
      veg_vsm = ws
      ch    =cdh
      veg_ch = cdh
      srht  =srheat
      veg_srht = srheat
      trht  =trheat
  !    zs    =zgs  !!! will not need after qsbal is replaced
  !    eddy  =khs   !!! will not need after qsbal is replaced
c     tspass=ts
c***********************************************************************
c****
c**** calculate ground fluxes
c     call qsbal

      call ghinij (i,j)
      call veg_set_cell(i,j)
      !call init_localveg
      call advnc
      call evap_limits( .false., evap_max_ij(i,j), fr_sat_ij(i,j) )

      call veg_save_cell(i,j)

      tg1=tbcs
      !qg_ij(i,j) = qs  !!! - this seemed to work ok
      !! trying more precise value for qg :  qsat(tg1+tf,elhx,ps)
      qg_sat = qsat(tg1+tf,elhx,ps) ! saturated soil
      qg_nsat = qs              ! non-sat soil, no evap
      if ( rcdhws > 1.d-30 ) then   ! correction to non-sat, due to evap
        qg_nsat = qg_nsat + evap_max_ij(i,j)/(0.001*rcdhws)
      endif
      qg_nsat = min( qg_nsat, qg_sat )
      qg_ij(i,j) = fr_sat_ij(i,j) * qg_sat
     &     + (1.d0 -fr_sat_ij(i,j)) * qg_nsat

      wbare(1:ngm,i,j) = w(1:ngm,1)
      wvege(0:ngm,i,j) = w(0:ngm,2)
      htbare(0:ngm,i,j) = ht(0:ngm,1)
      htvege(0:ngm,i,j) = ht(0:ngm,2)
      snowbv(1:2,i,j)   = snowd(1:2)

!!! test test test
!      aevap = 0.d0

cddd#ifdef TRACERS_WATER_OLD
cdddC**** reset tracer variables
cdddc      wsoil_tot=wsoil_tot-(aevapw+aevapd+aevapb+aruns+arunu)/rhow
cddd      wsoil_tot=wsoil_tot-(aevap+aruns+arunu)/rhow
cddd      do nx=1,ntx
cddd        n=ntix(nx)
cddd        if (tr_wd_TYPE(n).eq.nWATER) THEN
cdddc**** fix outputs to mean ratio (TO BE REPLACED WITHIN SOIL TRACERS)
cddd        trruns(nx)=aruns * trsoil_rat(nx) ! kg/m^2
cddd        trrunu(nx)=arunu * trsoil_rat(nx)
cdddc#ifdef TRACERS_SPECIAL_O18
cdddc        tevapw(nx)=aevapw*trsoil_rat(nx)*fracvl(tp(1,1),trname(n))
cdddc#else
cdddc        tevapw(nx)=aevapw * trsoil_rat(nx)
cdddc#endif
cdddc        tevapd(nx)=aevapd * trsoil_rat(nx)
cdddc        tevapb(nx)=aevapb * trsoil_rat(nx)
cddd        tevapw(nx)=aevap * trsoil_rat(nx)
cdddc**** update ratio
cc       trsoil_tot(nx)=trsoil_tot(nx)-(tevapw(nx)+tevapd(nx)+tevapb(nx)
cdddc     *       +trruns(nx)+trrunu(nx))
c       trsoil_tot(nx)=trsoil_tot(nx)-(tevapw(nx)+trruns(nx)+trrunu(nx))
cddd        trsoil_rat(nx)=trsoil_tot(nx)/(rhow*wsoil_tot)
cddd        trbare(n,1:ngm,i,j) = trsoil_rat(nx)*w(1:ngm,1)*rhow
cddd        trvege(n,:,i,j) = trsoil_rat(nx)*w(:,2)*rhow
cddd        trsnowbv(n,:,i,j)=trsoil_rat(nx)*snowd(:)*rhow
cdddc       trbare(n,:,i,j) = trw(nx,:,1)
cdddc       trvege(n,:,i,j) = trw(nx,:,2)
cdddc       trsnowbv(n,:,i,j) = trsnowd(nx,:)
cddd        gtracer(n,itype,i,j)=trsoil_rat(nx)
cddd        end if
cddd      end do
cddd#endif
ccc copy snow variables back to storage
      nsn_ij    (1:2, i, j)         = nsn(1:2)
      !isn_ij    (1:2, i, j)         = isn(1:2)
      dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
      wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
      hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
      fr_snow_ij(1:2, i, j)         = fr_snow(1:2)
ccc Save canopy temperature.
      canopy_temp_ij(i,j) = tp(0,2)  !nyk
ccc tracers
#ifdef TRACERS_WATER
      do nx=1,ntg
        n = ntixw(nx)
        tr_wbare(n,1:ngm,i,j) = tr_w(nx,1:ngm,1)
        tr_wvege(n,0:ngm,i,j) = tr_w(nx,0:ngm,2)
        tr_wsn_ij(n,1:nlsn, 1:2, i, j) = tr_wsn(nx,1:nlsn,1:2)
      enddo
#endif

c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      if ( snow_cover_same_as_rad == 0 ) then
        do ibv=1,2
          call snow_cover(fr_snow_rad_ij(ibv,i,j),
     &         snowbv(ibv,i,j), top_dev_ij(i,j) )
          fr_snow_rad_ij(ibv,i,j) = min (
     &         fr_snow_rad_ij(ibv,i,j), fr_snow_ij(ibv, i, j) )
        enddo
      else
        do ibv=1,2
          fr_snow_rad_ij(ibv,i,j) = fr_snow_ij(ibv, i, j)
        enddo
      endif

      aij(i,j,ij_g18)=aij(i,j,ij_g18)+aevapb
      aij(i,j,ij_g19)=aij(i,j,ij_g19)+aevapd
      aij(i,j,ij_g20)=aij(i,j,ij_g20)+aevapw
      aij(i,j,ij_g05)=aij(i,j,ij_g05)+abetab/nisurf
      aij(i,j,ij_g06)=aij(i,j,ij_g06)+abetap/nisurf
      aij(i,j,ij_g11)=aij(i,j,ij_g11)+abeta/nisurf
      aij(i,j,ij_g12)=aij(i,j,ij_g12)+acna/nisurf
      aij(i,j,ij_g13)=aij(i,j,ij_g13)+acnc/nisurf
      aij(i,j,ij_gpp)=aij(i,j,ij_gpp)+agpp
      aij(i,j,ij_g26)=aij(i,j,ij_g26)+abetav/nisurf
      aij(i,j,ij_g27)=aij(i,j,ij_g27)+abetat/nisurf
      aij(i,j,ij_g14)=aij(i,j,ij_g14)+aepp
      if (moddsf.eq.0) then
        aij(i,j,ij_g15)=aij(i,j,ij_g15)+tp(1,1)
        aij(i,j,ij_g16)=aij(i,j,ij_g16)+tp(2,1)
        aij(i,j,ij_g17)=aij(i,j,ij_g17)+tp(3,1)
        aij(i,j,ij_g21)=aij(i,j,ij_g21)+tp(0,2)
        aij(i,j,ij_g22)=aij(i,j,ij_g22)+tp(1,2)
        aij(i,j,ij_g23)=aij(i,j,ij_g23)+tp(2,2)
        aij(i,j,ij_g24)=aij(i,j,ij_g24)+tp(3,2)
        aij(i,j,ij_g25)=aij(i,j,ij_g25)+fb*zw(1)+fv*zw(2)
      end if
      trhdt=trheat*dtsurf-atrg
c           for radiation find composite values over earth
c           for diagnostic purposes also compute gdeep 1 2 3
      snowe(i,j)=1000.*(snowd(1)*fb+snowd(2)*fv)
      tearth(i,j)=tg1
      wearth(i,j)=1000.*( fb*w(1,1)*(1.-fice(1,1)) +
     &     fv*(w(1,2)*(1.-fice(1,2))+w(0,2)*(1.-fice(0,2))) )
      aiearth(i,j)=1000.*( fb*w(1,1)*fice(1,1) +
     &     fv*(w(1,2)*fice(1,2)+w(0,2)*fice(0,2)) )
      call retp2 (tg2av,wtr2av,ace2av)
      gdeep(i,j,1)=tg2av
      gdeep(i,j,2)=wtr2av
      gdeep(i,j,3)=ace2av
      gtemp(1,4,i,j)=tearth(i,j)
c**** calculate fluxes using implicit time step for non-ocean points
      uflux1(i,j)=uflux1(i,j)+ptype*rcdmws*(us-uocean)
      vflux1(i,j)=vflux1(i,j)+ptype*rcdmws*(vs-vocean)
c**** accumulate surface fluxes and prognostic and diagnostic quantities
      !evap=aevapw+aevapd+aevapb
      evap = aevap
      evapor(i,j,4)=evapor(i,j,4)+evap
      evhdt=-alhg
C**** hack to correct energy flux
ccc      evhdt=-evap*lhe ! hopefully not needed any more
      shdt=-ashg
      dth1(i,j)=dth1(i,j)-shdt*ptype/(sha*ma1*p1k)
      dq1(i,j) =dq1(i,j)+evap*ptype/ma1
      qsavg(i,j)=qsavg(i,j)+qs*ptype
      qsats=qsat(ts,elhx,ps)
c**** save runoff for addition to lake mass/energy resevoirs
      runoe (i,j)=runoe (i,j)+ aruns+ arunu
      erunoe(i,j)=erunoe(i,j)+aeruns+aerunu
c****
      e0(i,j,4)=e0(i,j,4)+af0dt
      e1(i,j,4)=e1(i,j,4)+af1dt

#ifdef TRACERS_WATER
ccc accumulate tracer evaporation and runoff
      do nx=1,ntg
        n=ntixw(nx)
        trevapor(n,itype,i,j) = trevapor(n,itype,i,j) + atr_evap(nx)
        !trevapor(n,itype,i,j) = trevapor(n,itype,i,j) + aevap  !*rhow
        trunoe(n,i,j) = trunoe(n,i,j) + atr_rnff(nx)
        !trunoe(n,i,j) = trunoe(n,i,j) + (aruns+arunu)  !*rhow
        gtracer(n,itype,i,j) = atr_g(nx)   ! /dtsurf
        trsrfflx(i,j,n)=trsrfflx(i,j,n)+
     &       atr_evap(nx)/dtsurf *dxyp(j)*ptype
      enddo
#endif
#ifdef TRACERS_DRYDEP
      do nx=1,ntx
        n=ntix(nx)
ccc accumulate tracer dry deposition
        if(dodrydep(n)) then
          rtsdt=rhosrf*trs(nx)*dtsurf 
          tdryd=-rtsdt*(dep_vel(n)+gs_vel(n))          ! kg/m2
          tdd = tdryd*dxyp(j)*ptype                    ! kg
          td1 = (trsrfflx(i,j,n)+totflux(nx))*dtsurf   ! kg
          if (trm(i,j,1,n)+td1+tdd.lt.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd earth",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,trs(nx),trtop(nx)
            tdd= -max(trm(i,j,1,n)+td1,0d0)
            tdryd= tdd/(dxyp(j)*ptype)
            trsrfflx(i,j,n)= - trm(i,j,1,n)/dtsurf
          else
            trsrfflx(i,j,n)=trsrfflx(i,j,n)+tdd/dtsurf
          end if
! trdrydep downward flux by surface type (kg/m^2)
          trdrydep(n,itype,i,j)=trdrydep(n,itype,i,j) - tdryd
! diagnose turbulent and settling fluxes separately
          taijn(i,j,tij_drydep,n)=taijn(i,j,tij_drydep,n) +
     &         ptype*rtsdt*dep_vel(n)
          taijn(i,j,tij_gsdep ,n)=taijn(i,j,tij_gsdep ,n) +
     &         ptype*rtsdt* gs_vel(n)
          dtr_dd(j,n,1)=dtr_dd(j,n,1)-ptype*rtsdt*dxyp(j)*dep_vel(n)
          dtr_dd(j,n,2)=dtr_dd(j,n,2)-ptype*rtsdt*dxyp(j)* gs_vel(n)
        end if
      end do
#endif
c****
c**** accumulate diagnostics
c****
#ifdef TRACERS_ON
C**** Save surface tracer concentration whether calculated or not
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if (needtrs(n)) then
            nx=nx+1
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)+trs(nx)*ptype
          else
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *bydxyp(j),0d0)*ptype
          end if
        end if
      end do
#endif
ccc not sure about the code below. hopefully that''s what is meant above
#ifdef TRACERS_WATER
      do nx=1,ntg
        n=ntixw(nx)
        taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *       trevapor(n,itype,i,j)*ptype
        taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     *         gtracer(n,itype,i,j)*ptype
        taijn(i,j,tij_soil,n)=taijn(i,j,tij_soil,n) + (
     &       fb*(sum( tr_w(nx,1:ngm,1) ) + sum( tr_wsn(nx,1:nsn(1),1)))+
     &       fv*(sum( tr_w(nx,0:ngm,2) ) + sum( tr_wsn(nx,1:nsn(2),2) ))
     *       )
        tajls(j,1,jls_source(1,n))=tajls(j,1,jls_source(1,n))
     *       +trevapor(n,itype,i,j)*ptype
      enddo
#endif
!      print '(a,10(e12.4))','trevapor_trunoe',
!     &     trevapor(1,itype,i,j), trunoe(1,i,j)
!     &     ,aevap,atr_evap(1),aruns,arunu,atr_rnff(1),pr,trpr


      aij(i,j,ij_rune)=aij(i,j,ij_rune)+aruns
      aij(i,j,ij_arunu)=aij(i,j,ij_arunu)+arunu
      aij(i,j,ij_pevap)=aij(i,j,ij_pevap)+(aepc+aepb)

      if ( warmer >= 0 ) then
        if(ts.lt.tf) tsfrez(i,j,tf_day1)=timez
        tsfrez(i,j,tf_last)=timez
      else
        if ( tsfrez(i,j,tf_last)+.03 >= timez .and. ts >= tf )
     $       tsfrez(i,j,tf_last)=timez
      endif

      if(tg1.lt.tdiurn(i,j,1)) tdiurn(i,j,1)=tg1
      if(tg1.gt.tdiurn(i,j,2)) tdiurn(i,j,2)=tg1
      if(ts.lt.tdiurn(i,j,3)) tdiurn(i,j,3)=ts
      if(ts.gt.tdiurn(i,j,4)) tdiurn(i,j,4)=ts

c**** quantities accumulated for regions in diagj
ccc   areg(jr,j_trhdt)=areg(jr,j_trhdt)+trhdt*ptype*dxyp(j)
ccc   areg(jr,j_shdt )=areg(jr,j_shdt )+shdt*ptype*dxyp(j)
ccc   areg(jr,j_evhdt)=areg(jr,j_evhdt)+evhdt*ptype*dxyp(j)
ccc   areg(jr,j_evap )=areg(jr,j_evap )+evap*ptype*dxyp(j)
ccc   areg(jr,j_erun)=areg(jr,j_erun)+(aeruns+aerunu)*pearth*dxyp(j)
ccc   areg(jr,j_run )=areg(jr,j_run )+(aruns+arunu)*pearth*dxyp(j)
ccc   if ( moddsf == 0 )
ccc  $     areg(jr,j_tsrf )=areg(jr,j_tsrf )+(ts-tf)*ptype*dxyp(j)
      AREGIJ(1,I,J)  = trhdt*ptype*dxyp(j)
      AREGIJ(2,I,J)  = shdt*ptype*dxyp(j)
      AREGIJ(3,I,J)  = evhdt*ptype*dxyp(j)
      AREGIJ(4,I,J)  = evap*ptype*dxyp(j)
      AREGIJ(5,I,J)  = (aeruns+aerunu)*pearth*dxyp(j)
      AREGIJ(6,I,J)  = (aruns+arunu)*pearth*dxyp(j)
      if ( moddsf == 0 ) THEN
        AREGIJ(7,I,J)  = (ts-tf)*ptype*dxyp(j)
        AREGIJ(8,I,J)  = tg1    *ptype*dxyp(j)
        AREGIJ(9,I,J)  = tg2av  *ptype*dxyp(j)
      end if
c**** quantities accumulated for latitude-longitude maps in diagij
      aij(i,j,ij_shdt)=aij(i,j,ij_shdt)+shdt*ptype
      aij(i,j,ij_beta)=aij(i,j,ij_beta)+abetad/nisurf
      if(modrd.eq.0)aij(i,j,ij_trnfp0)=aij(i,j,ij_trnfp0)+trhdt*ptype
     *     /dtsrc
      aij(i,j,ij_srtr)=aij(i,j,ij_srtr)+(srhdt+trhdt)*ptype
      aij(i,j,ij_neth)=aij(i,j,ij_neth)+(srhdt+trhdt+shdt+evhdt)*ptype
      aij(i,j,ij_evap)=aij(i,j,ij_evap)+evap*ptype
      if ( moddsf == 0 ) then
        aij(i,j,ij_ws)=aij(i,j,ij_ws)+ws*ptype
        aij(i,j,ij_ts)=aij(i,j,ij_ts)+(ts-tf)*ptype
        aij(i,j,ij_us)=aij(i,j,ij_us)+us*ptype
        aij(i,j,ij_vs)=aij(i,j,ij_vs)+vs*ptype
        aij(i,j,ij_taus)=aij(i,j,ij_taus)+rcdmws*wsm*ptype
        aij(i,j,ij_tauus)=aij(i,j,ij_tauus)+rcdmws*(us-uocean)*ptype
        aij(i,j,ij_tauvs)=aij(i,j,ij_tauvs)+rcdmws*(vs-vocean)*ptype
        aij(i,j,ij_qs)=aij(i,j,ij_qs)+qs*ptype
        aij(i,j,ij_tg1)=aij(i,j,ij_tg1)+tg1*ptype
        aij(i,j,ij_pblht)=aij(i,j,ij_pblht)+dbl*ptype
chyd       aij(i,j,ij_arunu)=aij(i,j,ij_arunu)
chyd      *  +   (40.6*psoil+.72*(2.*(tss-tfs)-(qsatss-qss)*lhe/sha))
c**** quantities accumulated hourly for diagDD
      endif
      if ( moddd == 0 ) then
        do kr=1,ndiupt
          if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then
            adiurn(ih,idd_ts,kr)=adiurn(ih,idd_ts,kr)+ts*ptype
            adiurn(ih,idd_tg1,kr)=adiurn(ih,idd_tg1,kr)+(tg1+tf)*ptype
            adiurn(ih,idd_qs,kr)=adiurn(ih,idd_qs,kr)+qs*ptype
            adiurn(ih,idd_qg,kr)=adiurn(ih,idd_qg,kr)+qg*ptype
            adiurn(ih,idd_swg,kr)=adiurn(ih,idd_swg,kr)+srhdt*ptype
            adiurn(ih,idd_lwg,kr)=adiurn(ih,idd_lwg,kr)+trhdt*ptype
            adiurn(ih,idd_sh,kr)=adiurn(ih,idd_sh,kr)+shdt*ptype
            adiurn(ih,idd_lh,kr)=adiurn(ih,idd_lh,kr)+evhdt*ptype
            adiurn(ih,idd_hz0,kr)=adiurn(ih,idd_hz0,kr)
     *           +(srhdt+trhdt+shdt+evhdt)*ptype
            adiurn(ih,idd_ug,kr)=adiurn(ih,idd_ug,kr)+ug*ptype
            adiurn(ih,idd_vg,kr)=adiurn(ih,idd_vg,kr)+vg*ptype
            adiurn(ih,idd_wg,kr)=adiurn(ih,idd_wg,kr)+wg*ptype
            adiurn(ih,idd_us,kr)=adiurn(ih,idd_us,kr)+us*ptype
            adiurn(ih,idd_vs,kr)=adiurn(ih,idd_vs,kr)+vs*ptype
            adiurn(ih,idd_ws,kr)=adiurn(ih,idd_ws,kr)+ws*ptype
            adiurn(ih,idd_cia,kr)=adiurn(ih,idd_cia,kr)+psi*ptype
            adiurn(ih,idd_cm,kr)=adiurn(ih,idd_cm,kr)+cdm*ptype
            adiurn(ih,idd_ch,kr)=adiurn(ih,idd_ch,kr)+cdh*ptype
            adiurn(ih,idd_cq,kr)=adiurn(ih,idd_cq,kr)+cdq*ptype
            adiurn(ih,idd_eds,kr)=adiurn(ih,idd_eds,kr)+khs*ptype
            adiurn(ih,idd_dbl,kr)=adiurn(ih,idd_dbl,kr)+dbl*ptype
            adiurn(ih,idd_ev,kr)=adiurn(ih,idd_ev,kr)+evap*ptype
            hdiurn(ihm,idd_ts,kr)=hdiurn(ihm,idd_ts,kr)+ts*ptype
            hdiurn(ihm,idd_tg1,kr)=hdiurn(ihm,idd_tg1,kr)+(tg1+tf)*ptype
            hdiurn(ihm,idd_qs,kr)=hdiurn(ihm,idd_qs,kr)+qs*ptype
            hdiurn(ihm,idd_qg,kr)=hdiurn(ihm,idd_qg,kr)+qg*ptype
            hdiurn(ihm,idd_swg,kr)=hdiurn(ihm,idd_swg,kr)+srhdt*ptype
            hdiurn(ihm,idd_lwg,kr)=hdiurn(ihm,idd_lwg,kr)+trhdt*ptype
            hdiurn(ihm,idd_sh,kr)=hdiurn(ihm,idd_sh,kr)+shdt*ptype
            hdiurn(ihm,idd_lh,kr)=hdiurn(ihm,idd_lh,kr)+evhdt*ptype
            hdiurn(ihm,idd_hz0,kr)=hdiurn(ihm,idd_hz0,kr)
     *           +(srhdt+trhdt+shdt+evhdt)*ptype
            hdiurn(ihm,idd_ug,kr)=hdiurn(ihm,idd_ug,kr)+ug*ptype
            hdiurn(ihm,idd_vg,kr)=hdiurn(ihm,idd_vg,kr)+vg*ptype
            hdiurn(ihm,idd_wg,kr)=hdiurn(ihm,idd_wg,kr)+wg*ptype
            hdiurn(ihm,idd_us,kr)=hdiurn(ihm,idd_us,kr)+us*ptype
            hdiurn(ihm,idd_vs,kr)=hdiurn(ihm,idd_vs,kr)+vs*ptype
            hdiurn(ihm,idd_ws,kr)=hdiurn(ihm,idd_ws,kr)+ws*ptype
            hdiurn(ihm,idd_cia,kr)=hdiurn(ihm,idd_cia,kr)+psi*ptype
            hdiurn(ihm,idd_cm,kr)=hdiurn(ihm,idd_cm,kr)+cdm*ptype
            hdiurn(ihm,idd_ch,kr)=hdiurn(ihm,idd_ch,kr)+cdh*ptype
            hdiurn(ihm,idd_cq,kr)=hdiurn(ihm,idd_cq,kr)+cdq*ptype
            hdiurn(ihm,idd_eds,kr)=hdiurn(ihm,idd_eds,kr)+khs*ptype
            hdiurn(ihm,idd_dbl,kr)=hdiurn(ihm,idd_dbl,kr)+dbl*ptype
            hdiurn(ihm,idd_ev,kr)=hdiurn(ihm,idd_ev,kr)+evap*ptype
          end if
        end do
      endif
c**** quantities accumulated for surface type tables in diagj
      aj(j,j_evap ,itearth)=aj(j,j_evap ,itearth)+ evap*pearth
      aj(j,j_trhdt,itearth)=aj(j,j_trhdt,itearth)+trhdt*pearth
      aj(j,j_evhdt,itearth)=aj(j,j_evhdt,itearth)+evhdt*pearth
      aj(j,j_shdt ,itearth)=aj(j,j_shdt ,itearth)+ shdt*pearth
      aj(j,j_erun ,itearth)=aj(j,j_erun ,itearth)+(aeruns+aerunu)*pearth
      aj(j,j_run  ,itearth)=aj(j,j_run  ,itearth)+(aruns+arunu)*pearth
      if(moddsf.eq.0) then
        aj(j,j_tsrf,itearth)=aj(j,j_tsrf,itearth)+(ts-tf)*pearth
        aj(j,j_tg1 ,itearth)=aj(j,j_tg1 ,itearth)+    tg1*pearth
        aj(j,j_tg2 ,itearth)=aj(j,j_tg2 ,itearth)+  tg2av*pearth
        aj(j,j_type,itearth)=aj(j,j_type,itearth)+        pearth
      end if

      end do loop_i
      end do loop_j
!$OMP  END PARALLEL DO

      DO 825 J=J_0,J_1
      DO 825 I=1,IMAXJ(J)
         IF(FEARTH(I,J).LE.0.0)  GO TO 825
         JR=JREG(I,J)
         areg(jr,j_trhdt)=areg(jr,j_trhdt)+AREGIJ(1,I,J)
         areg(jr,j_shdt )=areg(jr,j_shdt )+AREGIJ(2,I,J)
         areg(jr,j_evhdt)=areg(jr,j_evhdt)+AREGIJ(3,I,J)
         areg(jr,j_evap )=areg(jr,j_evap )+AREGIJ(4,I,J)
         areg(jr,j_erun )=areg(jr,j_erun )+AREGIJ(5,I,J)
         areg(jr,j_run  )=areg(jr,j_run  )+AREGIJ(6,I,J)
         if( moddsf == 0 ) then
           areg(jr,j_tsrf)=areg(jr,j_tsrf)+AREGIJ(7,I,J)
           areg(jr,j_tg1 )=areg(jr,j_tg1 )+AREGIJ(8,I,J)
           areg(jr,j_tg2 )=areg(jr,j_tg2 )+AREGIJ(9,I,J)
         end if
  825 CONTINUE
C
      return
      end subroutine earth

      subroutine snow_cover( fract_snow, snow_water, top_dev )
!@sum computes snow cover from snow water eq. and topography
!@var fract_snow snow cover fraction (0-1)
!@var snow_water snow water equivalent (m)
!@var top_dev standard deviation of the surface elevation
      use constant, only : teeny
      real*8, intent(out) :: fract_snow
      real*8, intent(in) :: snow_water, top_dev

      ! using formula from the paper by A. Roesch et al
      ! (Climate Dynamics (2001), 17: 933-946)
      fract_snow =
ccc     $     .95d0 * tanh( 100.d0 * snow_water ) *
ccc                               currently using only topography part
     $     sqrt ( 1000.d0 * snow_water /
     $     (1000.d0 * snow_water + teeny + snow_cover_coef * top_dev) )

      end subroutine snow_cover


      subroutine init_gh(dtsurf,redogh,inisnow,istart)
c**** modifications needed for split of bare soils into 2 types
      use filemanager
      use param
      use constant, only : twopi,rhow,edpery,sha,lhe,tf
      use model_com, only : fearth,itime,nday,jeq,jyear
      use dagcom, only : npts,icon_wtg,icon_htg,conpt0
      use sle001
#ifdef TRACERS_WATER
      use tracer_com, only : ntm,tr_wd_TYPE,nwater,itime_tr0,needtrs
      use fluxes, only : gtracer
      use veg_com, only:  afb, avh
#endif
      use fluxes, only : gtemp
      use ghycom
      use dynamics, only : pedn
      use snow_drvm, only : snow_cover_coef2=>snow_cover_coef
     &     ,snow_cover_same_as_rad
      use veg_drv, only : init_vegetation
      use veg_com, only : vdata

      implicit none

      real*8, intent(in) :: dtsurf
      integer, intent(in) :: istart
      logical, intent(in) :: redogh, inisnow
      integer iu_soil,iu_top_index
      integer jday
      real*8 snowdp,wtr1,wtr2,ace1,ace2,tg1,tg2
      logical :: qcon(npts)
      integer i, j, ibv
      real*8 pearth
      logical ghy_data_missing
      character conpt(npts)*10
#ifdef TRACERS_WATER
      real*8 trsoil_tot,wsoil_tot,fm
#endif
c****
c**** contents of sdata(i,j,k):
c****       1 -   ngm   dz(ngm)
c****   ngm+1 - 6*ngm   q(is,ngm)
c**** 6*ngm+1 - 11*ngm   qk(is,ngm)
c**** 11*ngm+1           sl
      real*8, external :: qsat
!@dbparam ghy_default_data if == 1 reset all GHY data to defaults
!@+ (do not read it from files)
      integer :: ghy_default_data = 0

C****	define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

c**** set conservation diagnostics for ground water mass and energy
      conpt=conpt0
      conpt(4)="EARTH"
      qcon=(/ .false., .false., .false., .true., .false., .false.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,conpt,"GRND WTR","(kg/m^2)        ",
     *     "(10^-9 kg/s/m^2)",1d0,1d9,icon_wtg)
      qcon=(/ .false., .false., .false., .true., .false., .false.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,conpt,"GRND ENG","(10**6 J/m^2)   ",
     *     "(10^-3 W/m^2)   ",1d-6,1d3,icon_htg)

c**** read rundeck parameters
      call sync_param( "snow_cover_coef", snow_cover_coef )
! hack. snow_cover_coef should be moved to snow_drvm
      snow_cover_coef2 = snow_cover_coef
      call sync_param( "snow_cover_same_as_rad", snow_cover_same_as_rad)
      call sync_param( "snoage_def", snoage_def )
      call sync_param( "ghy_default_data", ghy_default_data )

c**** read land surface parameters or use defaults
      if ( ghy_default_data == 0 ) then ! read from files

        !!!if (istart.le.0) return ! avoid reading unneeded files
c**** read soils parameters
        call openunit("SOIL",iu_SOIL,.true.,.true.)
        call dread (iu_SOIL,dz_ij,im*jm*(11*ngm+1),dz_ij)
        call closeunit (iu_SOIL)
c**** read topmodel parameters
        call openunit("TOP_INDEX",iu_TOP_INDEX,.true.,.true.)
        call readt(iu_TOP_INDEX,0,top_index_ij,im*jm,top_index_ij,1)
        call readt(iu_TOP_INDEX,0,top_dev_ij,  im*jm,top_dev_ij,  1)
        call closeunit (iu_TOP_INDEX)
      else  ! reset to default data
        if ( istart>0 .and. istart<10 ) then ! reset all
          call reset_gh_to_defaults( .true. )
        else   ! do not reset ghy prognostic variables
          call reset_gh_to_defaults( .false. )
        endif
        !!!if (istart.le.0) return
      endif


      if (istart.le.0) return
c****
c**** initialize constants
c****
c**** time step for ground hydrology
      dt=dtsurf
c spgsn is the specific gravity of snow
      spgsn=.1d0

c**** check whether ground hydrology data exist at this point.
      ghy_data_missing = .false.
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).gt.0) then
            if ( top_index_ij(i,j).eq.-1. ) then
              print *,"No top_index data: i,j=",i,j,top_index_ij(i,j)
              ghy_data_missing = .true.
            end if
            if ( sum(dz_ij(i,j,1:ngm)).eq.0 ) then
              print *, "No soil data: i,j=",i,j,dz_ij(i,j,1:ngm)
              ghy_data_missing = .true.
            endif
            if (wbare(1,i,j) < 1.d-10 .and. wvege(1,i,j) < 1.d-10) then
              print*,"No gh data in restart file: i,j=",i,j,
     &             wbare(:,i,j),wvege(:,i,j)
              ghy_data_missing = .true.
            endif
          end if
        enddo
      enddo
      if ( ghy_data_missing ) then
        write(6,*) 'Ground Hydrology data is missing at some pts'
        write(6,*) 'If you have a non-standard land mask, please'
        write(6,*) 'consider using extended GH data and rfs file.'
        call stop_model(
     &       'Ground Hydrology data is missing at some cells',255)
      endif

ccc read and initialize vegetation here
      call init_vegetation(redogh,istart)

      ! no need to continue computations for postprocessing
      if (istart.le.0) return

      call hl0

c****
c      print *,' '
c      print *,'soils parameters'
c      sdstnc=100.
c      print *,'interstream distance (m) sdstnc:',sdstnc
c      c1=90.
c      print *,'canopy conductance related parameter c1:',c1
c      prfr=.1d0
c      print *,'fraction (by area) of precipitation prfr:',prfr
c      print *,' '
c****
c code transplanted from subroutine input
c**** recompute ground hydrology data if necessary (new soils data)
      if (redogh) then
        jday=1+mod(itime/nday,365)
        cosday=cos(twopi/edpery*jday)
        sinday=sin(twopi/edpery*jday)

        do j=J_0,J_1
        do i=1,im
          pearth=fearth(i,j)
          if(pearth.le.0.) then

            wbare(:,i,j)=0.
            wvege(:,i,j)=0.
            htbare(:,i,j)=0.
            htvege(:,i,j)=0.
            snowbv(:,i,j)=0.

          else
ccc   ??? remove next 5 lines? -check the old version
            w(1:ngm,1) =   wbare(1:ngm,i,j)
            w(0:ngm,2) =   wvege(0:ngm,i,j)
            ht(0:ngm,1) = htbare(0:ngm,i,j)
            ht(0:ngm,2) = htvege(0:ngm,i,j)
            snowd(1:2) =  snowbv(1:2,i,j)

c**** compute soil heat capacity and ground water saturation gws
            call ghinij (i,j)
c**** fill in soils common blocks
            snowdp=snowe(i,j)/rhow
            wtr1=wearth(i,j)
            ace1=aiearth(i,j)
            tg1 =tearth(i,j)
            wtr2=wtr1
            ace2=ace1
            tg2 =tg1
c           wtr2=gdata(i,j,9)   ! this cannot be right
c           ace2=gdata(i,j,10)
c           tg2 =gdata(i,j,8)
            call ghinht (snowdp, tg1,tg2, wtr1,wtr2, ace1,ace2)

c**** copy soils prognostic quantities to model variables
            wbare(1:ngm,i,j) = w(1:ngm,1)
            wvege(0:ngm,i,j) = w(0:ngm,2)
            htbare(0:ngm,i,j) = ht(0:ngm,1)
            htvege(0:ngm,i,j) = ht(0:ngm,2)
            snowbv(1:2,i,j)   = snowd(1:2)
          end if
        end do
        end do
        write (*,*) 'ground hydrology data was made from ground data'
      end if
c**** set gtemp array
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).gt.0) then
            gtemp(1,4,i,j)=tearth(i,j)
          end if
        end do
      end do

C GISS-ESMF EXCEPTIONAL CASE
C-BMP Global sum on evap_max_ij

ccc if not initialized yet, set evap_max_ij, fr_sat_ij, qg_ij
ccc to something more appropriate
      if ( sum(evap_max_ij(:,:)) > im*jm-1.d0 ) then ! old default
        do j=J_0,J_1
          do i=1,im
            if ( fearth(i,j) .le. 0.d0 ) cycle
            qg_ij(i,j) = qsat(tearth(i,j)+tf,lhe,pedn(1,i,j))
          enddo
        enddo
        fr_sat_ij(:,:) = 0.d0
        evap_max_ij(:,:) = 0.d0
      endif

ccc   init snow here
ccc hope this is the right place to split first layer into soil
ccc and snow  and to set snow arrays
ccc!!! this should be done only when restarting from an old
ccc!!! restart file (without snow model data)

      if (inisnow) then
        do j=J_0,J_1
        do i=1,im
          pearth=fearth(i,j)
          if(pearth.le.0.) then
            nsn_ij(:,i,j)     = 0
            !isn_ij(:,i,j)     = 0
            dzsn_ij(:,:,i,j)  = 0.
            wsn_ij(:,:,i,j)   = 0.
            hsn_ij(:,:,i,j)   = 0.
            fr_snow_ij(:,i,j) = 0.
          else
            jday=1+mod(itime/nday,365)
            cosday=cos(twopi/edpery*jday)
            sinday=sin(twopi/edpery*jday)

            w(1:ngm,1) =   wbare(1:ngm,i,j)
            w(0:ngm,2) =   wvege(0:ngm,i,j)
            ht(0:ngm,1) = htbare(0:ngm,i,j)
            ht(0:ngm,2) = htvege(0:ngm,i,j)
            snowd(1:2) =  snowbv(1:2,i,j)

            call ghinij (i,j)
            call set_snow

            nsn_ij    (1:2, i, j)         = nsn(1:2)
            !isn_ij    (1:2, i, j)         = isn(1:2)
            dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
            wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
            hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
            fr_snow_ij(1:2, i, j)         = fr_snow(1:2)

c****     copy soils prognostic quantities to model variables
             wbare(1:ngm,i,j) = w(1:ngm,1)
             wvege(0:ngm,i,j) = w(0:ngm,2)
            htbare(0:ngm,i,j) = ht(0:ngm,1)
            htvege(0:ngm,i,j) = ht(0:ngm,2)
            snowbv(1:2,i,j)   = snowd(1:2)

          end if
        end do
        end do
      end if

c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      fr_snow_rad_ij(:,:,:) = 0.d0
      do j=J_0,J_1
        do i=1,im
          if ( fearth(i,j) > 0.d0 ) then
            do ibv=1,2
              call snow_cover(fr_snow_rad_ij(ibv,i,j),
     &             snowbv(ibv,i,j), top_dev_ij(i,j) )
              fr_snow_rad_ij(ibv,i,j) = min (
     &             fr_snow_rad_ij(ibv,i,j), fr_snow_ij(ibv, i, j) )
            enddo
          endif
        enddo
      enddo

#ifdef TRACERS_WATER
ccc still not quite correct (assumes fw=1)
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).le.0.d0) cycle
          fb=afb(i,j) ; fv=1.-fb
          fm=1.d0-exp(-snowbv(2,i,j)/((avh(i,j)*spgsn) + 1d-12))
          if ( fm < 1.d-3 ) fm=0.d0
          wsoil_tot=fb*( wbare(1,i,j)*(1.d0-fr_snow_ij(1,i,j))
     &     + wsn_ij(1,1,i,j)*fr_snow_ij(1,i,j) )
     &     + fv*( wvege(0,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))   !*1.d0
     &     + wsn_ij(1,2,i,j)*fm*fr_snow_ij(2,i,j) )
          do n=1,ntm
            if (itime_tr0(n).gt.itime) cycle
            if ( .not. needtrs(n) ) cycle
            ! should also restrict to TYPE=nWATER ?
            if ( wsoil_tot > 1.d-30 ) then
            gtracer(n,4,i,j) = (
     &           fb*( tr_wbare(n,1,i,j)*(1.d0-fr_snow_ij(1,i,j))
     &           + tr_wsn_ij(n,1,1,i,j) )         !*fr_snow_ij(1,i,j)
     &           + fv*( tr_wvege(n,0,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))
     &           + tr_wsn_ij(n,1,2,i,j)*fm ) )    !*fr_snow_ij(2,i,j)
     &           /(rhow*wsoil_tot)
            else
              gtracer(n,4,i,j) = 0.
            end if
          enddo
        end do
      end do
#endif

      return
      end subroutine init_gh


      subroutine reset_gh_to_defaults( reset_prognostic )
      !use model_com, only: vdata
      use ghycom
      use veg_drv, only : reset_veg_to_defaults
      logical, intent(in) :: reset_prognostic
      integer i,j

C****	define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

ccc ugly, should fix later
      call reset_veg_to_defaults( reset_prognostic )

      do j=J_0,J_1
      do i=1,im

      dz_ij(i,j,1:ngm)= (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      q_ij(i,j,1:imt,1:ngm)=
     &     reshape( (/  0.33491623d+00,  0.52958947d+00,
     &     0.13549370d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.32995611d+00,  0.52192056d+00,  0.14812243d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.32145596d+00,
     &     0.48299056d+00,  0.19555295d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.47638881d+00,  0.40400982d+00,
     &     0.11959970d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.99985123d-01,  0.95771909d-01,  0.41175738d-01,
     &     0.00000000d+00,  0.76306665d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.10000000d+01 /), (/imt,ngm/) )
      qk_ij(i,j,1:imt,1:ngm)=
     &     reshape( (/  0.34238762d+00,  0.52882469d+00,
     &     0.12878728d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.32943058d+00,  0.52857041d+00,  0.14199871d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.30698991d+00,
     &     0.52528000d+00,  0.16772974d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.39890009d+00,  0.43742162d+00,
     &     0.16367787d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.46536058d+00,  0.39922065d+00,  0.13541836d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.10000000d+01 /), (/imt,ngm/) )
      sl_ij(i,j)= 0.22695422d+00
      top_index_ij(i,j)= 0.10832934d+02
      top_dev_ij(i,j)= 0.21665636d+03

      if ( .not. reset_prognostic ) cycle

      snowe(i,j)= 0.65458111d-01
      tearth(i,j)= -0.12476520d+00
      wearth(i,j)=  0.29203081d+02
      aiearth(i,j)=  0.93720329d-01
      wbare(:,i,j) = (/  0.17837750d-01,  0.40924843d-01,
     &     0.77932012d-01,  0.11919649d+00,  0.57237469d-01,
     &     0.10000000d-11 /)
      wvege(:,i,j) = (/  0.10000000d-11,  0.29362259d-01,
     &     0.50065177d-01,  0.82533140d-01,  0.10383620d+00,
     &     0.31552459d-01,  0.10000000d-11 /)
      htbare(:,i,j)= (/  0.00000000d+00, -0.15487181d+07,
     &     -0.50720067d+07,  0.18917623d+07,  0.77174974d+07,
     &     0.21716040d+08,  0.44723067d+08 /)
      htvege(:,i,j)= (/ -0.13991376d+05, -0.53165599d+05,
     &     0.65443775d+06,  0.29276050d+07,  0.81455096d+07,
     &     0.21575081d+08,  0.45952255d+08 /)
      snowbv(:,i,j)= (/  0.00000000d+00,  0.65458111d-04 /)

      enddo
      enddo

      end subroutine reset_gh_to_defaults


      subroutine ghinij (i0,j0)
c**** input:
c**** avh(i,j) - array of vegetation heights
c**** spgsn - specific gravity of snow
c**** output:
c**** vh - vegetation height
c**** snowm - snow masking depth
c**** wfcap - water field capacity of top soil layer, m
c****
      use snow_model, only : i_earth,j_earth
      use sle001, only : dz,qk,ng,zb,zc,q,sl,xklh0 !spgsn,
     *     ,fb,fv,prs,ijdebug,n
     *     ,thets,thetm,ws,thm,nth,shc,shw,htprs,pr !shcap,shtpr,
     *     ,htpr
     *     ,top_index,top_stdev
      use ghycom, only : ngm,imt,dz_ij,sl_ij,q_ij,qk_ij
     *     ,top_index_ij,top_dev_ij
      use veg_com, only: afb
!      use veg_drv, only : veg_set_cell


      implicit none
      integer i0,j0
!      real*8 wfcap
      integer k,ibv,i
      real*8 shtpr
!----------------------------------------------------------------------!
      real*8, parameter :: shcap(imt) = (/2d6,2d6,2d6,2.5d6,2.4d6/)


      ijdebug=i0*1000+j0
      i_earth = i0
      j_earth = j0

ccc setting vegetation
 !     call veg_set_cell(i0,j0)

ccc passing topmodel parameters
      top_index = top_index_ij(i0, j0)
      top_stdev = top_dev_ij(i0, j0)
c**** set up layers
      dz(1:ngm)=dz_ij(i0,j0,1:ngm)
      q(1:imt,1:ngm)=q_ij(i0,j0,1:imt,1:ngm)
      qk(1:imt,1:ngm)=qk_ij(i0,j0,1:imt,1:ngm)
      sl=sl_ij(i0,j0)


      do n=1,ngm
        if(dz(n).le.0.) go to 21
      end do
   21 n=n-1
   !  print *,'>>>>>> n= ', n
      !n = ngm
      if(n.le.0) then
         write (99,*) 'ghinij:  n <= 0:  i,j,n=',i0,j0,n,(dz(k),k=1,43)
         call stop_model('stopped in GHY_DRV.f',255)
      end if

c**** calculate the boundaries, based on the thicknesses.
      zb(1)=0.
      do k=1,n
        zb(k+1)=zb(k)-dz(k)
      end do
c**** calculate the layer centers, based on the boundaries.
      do k=1,n
        zc(k)=.5*(zb(k)+zb(k+1))
      end do
c**** fb,fv: bare, vegetated fraction (1=fb+fv)
      fb=afb(i0,j0)
      fv=1.-fb

c****
      do ibv=1,2
        do k=1,n
          thets(k,ibv)=0.
          thetm(k,ibv)=0.
          do i=1,imt-1
            thets(k,ibv)=thets(k,ibv)+q(i,k)*thm(0,i)
            thetm(k,ibv)=thetm(k,ibv)+q(i,k)*thm(nth,i)
          end do
          ws(k,ibv)=thets(k,ibv)*dz(k)
        end do
      end do
!veg      ws(0,2)=.0001d0*alai  ! from call veg_set_cell above
  !    wfcap=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
c****
      call xklh0
c****
      do ibv=1,2
        do k=1,n
          shc(k,ibv)=0.
          do i=1,imt
            shc(k,ibv)=shc(k,ibv)+q(i,k)*shcap(i)
          end do
          shc(k,ibv)=(1.-thets(k,ibv))*shc(k,ibv)*dz(k)
        end do
      end do
c****
c shc(0,2) is the heat capacity of the canopy
!veg      aa=ala(1,i0,j0)
!veg      shc(0,2)=(.010d0+.002d0*aa+.001d0*aa**2)*shw
c****
c htpr is the heat of precipitation.
c shtpr is the specific heat of precipitation.
      shtpr=0.
      if(pr.gt.0.)shtpr=htpr/pr
c htprs is the heat of large scale precipitation
      htprs=shtpr*prs
c****
      return
      end subroutine ghinij

      subroutine ghinht (snowdp,tg1,tg2,wtr1,wtr2,ace1,ace2)
c**** initializes new ground (w,ht,snw) from old (t,w,ice,snw)
c**** evaluates the heat in the soil layers based on the
c**** temperatures.
c**** input:
c**** w - water in soil layers, m
c**** tp - temperature of layers, c
c**** fice - fraction of ice of layers
c**** fsn - heat of fusion of water
c**** shc - specific heat capacity of soil
c**** shi - specific heat capacity of ice
c**** shw - specific heat capcity of water
c**** snowd - snow depth, equivalent water m
c**** output:
c**** ht - heat in soil layers
c**** add calculation of wfc2
c**** based on combination of layers 2-n, as in retp2
      use sle001
      implicit none

      real*8 snowdp,tg1,tg2,wtr1,wtr2,ace1,ace2
      real*8 wfc1, wfc2, wet1, wet2, wmin, fbv
      integer k, ibv, ll

      wfc1=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
      wfc2=0.
      fbv=fb
      do 30 ibv=1,2
      do 20 k=2,n
      wfc2=wfc2+fbv*ws(k,ibv)
   20 continue
      fbv=fv
   30 continue
      wfc1=1000.*wfc1
      wfc2=1000.*wfc2
      fice(0,2)=1.
      fice(1,1)=(ace1+snowdp*1000.)/(wtr1+ace1+snowdp*1000.+1.d-20)
      fice(1,2)=fice(1,1)
      tp(0,2)=tg1
c**** w = snow(if top layer) + wmin + (wmax-wmin)*(wtr+ice)/wfc
      w(0,2)=0.
      do ibv=1,2
        w(1,ibv)=snowdp
        wmin=thetm(1,ibv)*dz(1)
        wet1=(wtr1+ace1)/(wfc1+1.d-20)
        if(wet1.gt.1.) wet1=1.
        w(1,ibv)=w(1,ibv)+wmin+(ws(1,ibv)-wmin)*wet1
        snowd(ibv)=snowdp
        tp(1,ibv)=tg1
        do k=2,n
          fice(k,ibv)=ace2/(wtr2+ace2+1.d-20)
          wmin=thetm(k,ibv)*dz(k)
          wet2=(wtr2+ace2)/(wfc2+1.d-20)
          if(wet2.gt.1.) wet2=1.
          w(k,ibv)=wmin+(ws(k,ibv)-wmin)*wet2
          tp(k,ibv)=tg2
        end do
      end do
c****
      entry ghexht
c****
c**** compute ht (heat w/m+2)
      do ibv=1,2
        ll=2-ibv
        do k=ll,n
          if(tp(k,ibv)) 2,4,6
 2        ht(k,ibv)=tp(k,ibv)*(shc(k,ibv)+w(k,ibv)*shi)-w(k,ibv)*fsn
          cycle
 4        ht(k,ibv)=-fice(k,ibv)*w(k,ibv)*fsn
          cycle
 6        ht(k,ibv)=tp(k,ibv)*(shc(k,ibv)+w(k,ibv)*shw)
        end do
      end do
      if(ijdebug.eq.0)then
       write(99,*)'ghinht id check',ijdebug
       write(99,*)'tg1,tg2',tg1,tg2
       write(99,*)'tp',tp
       write(99,*)'ht',ht
       write(99,*)'w',w
       write(99,*)'wtr1,wtr2',wtr1,wtr2
       write(99,*)'ace1,ace2',ace1,ace2
       write(99,*)'wfc1,wfc2',wfc1,wfc2
       write(99,*)'shc',shc
       write(99,*)'fice',fice
      endif
      return
      end subroutine ghinht

      subroutine retp2 (tg2av,wtr2av,ace2av)
c**** evaluates the mean temperature in the soil layers 2-ngm
c**** as well as the water and ice content.
c**** input:
c**** w - water in soil layers, m
c**** ht - heat in soil layers
c**** fsn - heat of fusion of water
c**** shc - specific heat capacity of soil
c**** shi - specific heat capacity of ice
c**** shw - specific heat capcity of water
c**** output:
c**** tg2av - temperature of layers 2 to ngm, c
c**** ice2av - ice amount in layers 2 to ngm, kg/m+2
c**** wtr2av - water in layers 2 to ngm, kg/m+2
      use sle001
      implicit none
      real*8 tg2av,wtr2av,ace2av, wc,htc,shcc,tpc,ficec,ftp
      integer k, ibv
      tg2av=0.
      wtr2av=0.
      ace2av=0.
      do 3500 ibv=1,2
      wc=0.
      htc=0.
      shcc=0.
      do k=2,n
        wc=wc+w(k,ibv)
        htc=htc+ht(k,ibv)
        shcc=shcc+shc(k,ibv)
      end do
      tpc=0.
      ficec=0.
      if(wc.ne.0.)  ficec=-htc/(fsn*wc)
      if(fsn*wc+htc.ge.0.)go to 3430
      tpc=(htc+wc*fsn)/(shcc+wc*shi)
      ficec=1.
      go to 3440
 3430 if(htc.le.0.) go to 3440
      tpc=htc/(shcc+wc*shw)
      ficec=0.
 3440 continue
      ftp=fb
      if(ibv.eq.2) ftp=fv
      tg2av=tg2av+tpc*ftp
      wtr2av=wtr2av+wc*ftp*1000.*(1.-ficec)
      ace2av=ace2av+wc*ftp*1000.*ficec
 3500 continue
      return
      end subroutine retp2

      subroutine checke(subr)
!@sum  checke checks whether arrays are reasonable over earth
!@auth original development team
!@ver  1.0
      use model_com, only : fearth,itime,wfcs
      use geom, only : imaxj
      use ghycom, only : tearth,wearth,aiearth,snowe,wbare,wvege,htbare
     *     ,htvege,snowbv,ngm
      implicit none

      real*8 x,tgl,wtrl,acel
      integer i,j
!@var subr identifies where check was called from
      character*6, intent(in) :: subr

C****	define local grid
      integer I_0, I_1
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

c**** check for nan/inf in earth data
      call check3(wbare(1:ngm,I_0:I_1,J_0:J_1) ,ngm  ,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'wb')
      call check3(wvege(0:ngm,I_0:I_1,J_0:J_1) ,ngm+1,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'wv')
      call check3(htbare(0:ngm,I_0:I_1,J_0:J_1),ngm+1,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'hb')
      call check3(htvege(0:ngm,I_0:I_1,J_0:J_1),ngm+1,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'hv')
      call check3(snowbv(1:ngm,I_0:I_1,J_0:J_1),2    ,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'sn')

c**** check for reasonable temperatures over earth
      x=1.001
      do j=J_0,J_1
        do i=1,imaxj(j)
          if (fearth(i,j).gt.0.) then
            tgl=tearth(i,j)
            wtrl=wearth(i,j)
            acel=aiearth(i,j)
            if ((tgl+60.)*(60.-tgl).le.0.) write (6,901) subr,i,j,itime
     *           ,fearth(i,j),'tg1 ',snowe(i,j),tgl,wtrl,acel
            if (wtrl.lt.0..or.acel.lt.0..or.(wtrl+acel).gt.x*wfcs(i
     *           ,j)) write(6,901) subr,i,j,itime,fearth(i,j),'wtr '
     *           ,snowe(i,j),tgl,wtrl,acel,wfcs(i,j)
          end if
        end do
      end do

      return
 901  format ('0gdata off, subr,i,j,i-time,pearth,',a7,2i4,i10,f5.2,1x
     *     ,a4/' snw,x,tg1,wtr1,ice1, wfc1 ',6f12.4)

      end subroutine checke

      subroutine daily_earth(end_of_day)
!@sum  daily_earth performs daily tasks for earth related functions
!@auth original development team
!@ver  1.0
!@calls RDLAI
      use constant, only : rhow,twopi,edpery,tf
      use model_com, only : nday,nisurf,jday,jyear,fearth,wfcs
      use veg_com, only : vdata                 !nyk
      use geom, only : imaxj
      use dagcom, only : aij,tdiurn,ij_strngts,ij_dtgdts,ij_tmaxe
     *     ,ij_tdsl,ij_tmnmx,ij_tdcomp, ij_dleaf
      use ghycom, only : snoage, snoage_def
      use veg_com, only : almass,aalbveg       !nyk
      use vegetation, only: crops_yr,cond_scheme !nyk
      use surf_albedo, only: albvnh  !nyk
      use sle001, only : fb,fv,ws
      use veg_drv, only : veg_set_cell

      implicit none
      real*8 tsavg,wfc1
      real*8 aleafmass, aalbveg0, fvp, sfv  !nyk veg ! , aleafmasslast
      integer i,j,itype
      integer northsouth,iv  !nyk
      logical, intent(in) :: end_of_day

C****	define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C**** Update vegetation file if necessary  (i.e. if crops_yr=0)
      if(crops_yr.eq.0) call updveg(jyear,.true.)
      if(cond_scheme.eq.2) call updsur (0,jday)
c****
c**** find leaf-area index & water field capacity for ground layer 1
c****
      cosday=cos(twopi/edpery*jday)
      sinday=sin(twopi/edpery*jday)
      do j=J_0,J_1
        if(j.le.jm/2) then      !nyk added northsouth
          northsouth=1.d0       !southern hemisphere
        else
          northsouth=2.d0       !northern hemisphere
        end if
        do i=1,im
          wfcs(i,j)=24.
          if (fearth(i,j).gt.0.) then
            !-----------------------------------------------------------
            !nyk Update vegetation albedos.
            !WARNING:  ALBVNH is not saved in a restart file.
            !          Must run without restarting.
            !albvnh(9,6,2)=albvnh(1+8veg,6bands,2hemi), band 1 is VIS.
            !if RUNDECK selects new conductance scheme
            if (cond_scheme.eq.2) then
              aalbveg0 = 0.d0
              sfv=0.d0
              do iv=1,8
                fvp=vdata(i,j,iv+1)
                sfv=sfv+fvp
                aalbveg0 = aalbveg0 + fvp*(ALBVNH(iv+1,1,northsouth))
                !write (99,*) 'fvp',fvp
                !write (99,*) 'ALBVNH',ALBVNH(iv+1,1,northsouth)
              end do
              aalbveg(i,j) = 0.08D0
              if(sfv.gt.0.) aalbveg(i,j) = aalbveg0/sfv !nyk
             !write (99,*) 'daily aalbveg', aalbveg(i,j)
            end if
            !-----------------------------------------------------------

            call ghinij(i,j)
            call veg_set_cell(i,j,.true.)
            wfc1=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
            wfcs(i,j)=rhow*wfc1 ! canopy part changes

            !-----------------------------------------------------------
            !nyk - TEMPORARY calculate change in leaf mass per day
            !get aleafmass(i,j) at jday
            aleafmass=
     $           almass(1,i,j)+cosday*almass(2,i,j)+sinday*almass(3,i,j)

            !Calculate dlmass(i,j) increment from last jday
            !cosdaym1=cos(twopi/edpery*(jday-1))
            !sindaym1=sin(twopi/edpery*(jday-1))
            !aleafmasslast=almass(1,i,j)+cosdaym1*almass(2,i,j)+
!     $      !     sindaym1*almass(3,i,j)
            !accumulate dlmass
            !adlmass = aleafmass - aleafmasslast
            adlmass = aleafmass
            !aij(i,j,ij_dleaf)=aij(i,j,ij_dleaf)+adlmass
            aij(i,j,ij_dleaf)=adlmass  !accumulate just instant. value
            !PRINT '(F4.4)',adlmass                            !DEBUG
            !call stop_model('Just did adlmass',255)           !DEBUG
          end if
        end do
      end do

      if (end_of_day) then
        do j=J_0,J_1
        do i=1,imaxj(j)
c****
c**** increase snow age depending on snoage_def
c****
          if (snoage_def.eq.0) then ! update indep. of ts
            do itype=1,3
              snoage(itype,i,j)=1.+.98d0*snoage(itype,i,j)
            end do
          elseif (snoage_def.eq.1) then ! update if max T>0
            if (tdiurn(i,j,7).gt.0) snoage(1,i,j)=1.+.98d0
     *           *snoage(1,i,j) ! ocean ice (not currently used)
            if (tdiurn(i,j,8).gt.0) snoage(2,i,j)=1.+.98d0
     *           *snoage(2,i,j) ! land ice
            if (tdiurn(i,j,2).gt.0) snoage(3,i,j)=1.+.98d0
     *           *snoage(3,i,j) ! land
          else
            write(6,*) "This snoage_def is not defined: ",snoage_def
            write(6,*) "Please use: 0 (update indep of T)"
            write(6,*) "            1 (update if T>0)"
            call stop_model('stopped in GHY_DRV.f',255)
          end if
          tsavg=tdiurn(i,j,5)/(nday*nisurf)
          if(32.+1.8*tsavg.lt.65.)
     *         aij(i,j,ij_strngts)=aij(i,j,ij_strngts)+(33.-1.8*tsavg)
          aij(i,j,ij_dtgdts)=aij(i,j,ij_dtgdts)+18.*((tdiurn(i,j,2)-
     *         tdiurn(i,j,1))/(tdiurn(i,j,4)-tdiurn(i,j,3)+1.d-20)-1.)
          aij(i,j,ij_tdsl)=aij(i,j,ij_tdsl)+
     *         (tdiurn(i,j,4)-tdiurn(i,j,3))
          aij(i,j,ij_tdcomp)=aij(i,j,ij_tdcomp)+
     *         (tdiurn(i,j,6)-tdiurn(i,j,9))
          aij(i,j,ij_tmaxe)=aij(i,j,ij_tmaxe)+(tdiurn(i,j,4)-tf)
          if (tdiurn(i,j,6).lt.aij(i,j,ij_tmnmx))
     *         aij(i,j,ij_tmnmx)=tdiurn(i,j,6)
        end do
        end do
      end if
#ifdef TRACERS_DRYDEP
      CALL RDLAI ! read leaf area indices for tracer dry deposition
#endif

      return
      end subroutine daily_earth

      subroutine ground_e
!@sum  ground_e driver for applying surface fluxes to land fraction
!@auth original development team
!@ver  1.0
      use model_com, only : fearth,itearth
      use geom, only : imaxj,dxyp
      use ghycom, only : snowe, tearth,wearth,aiearth,wbare,wvege,snowbv
     *     ,fr_snow_ij,fr_snow_rad_ij, gdeep
      use veg_com, only : afb
      use dagcom, only : aj,areg,aij,jreg,ij_evap,ij_f0e,ij_evape
     *     ,ij_gwtr,ij_tg1,j_wtr1,j_ace1,j_wtr2,j_ace2
     *     ,j_snow,j_evap,j_type,ij_g01,ij_g07,ij_g28
     *     ,ij_g29,j_rsnow,ij_rsnw,ij_rsit,ij_snow
      use fluxes, only : e0,e1,evapor,eprec
      implicit none

      real*8 snow,tg1,tg2,f0dt,f1dt,evap,wtr1,wtr2,ace1,ace2
     *     ,pearth,enrgp,scove
      integer i,j,jr,k

C****	define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      do j=J_0,J_1
      do i=1,imaxj(j)
      pearth=fearth(i,j)
      jr=jreg(i,j)
      if (pearth.gt.0) then

        snow=snowe(i,j)
        tg1 = tearth(i,j)
        wtr1= wearth(i,j)
        ace1=aiearth(i,j)
        tg2=gdeep(i,j,1)
        wtr2=gdeep(i,j,2)
        ace2=gdeep(i,j,3)
        f0dt=e0(i,j,4)
        f1dt=e1(i,j,4)
        evap=evapor(i,j,4)
        enrgp=eprec(i,j)      ! including latent heat

c**** accumulate diagnostics
c**** the following is the actual snow cover of the snow model
c        scove = pearth *
c     *       ( afb(i,j)*fr_snow_ij(1,i,j)
c     *       + (1.-afb(i,j))*fr_snow_ij(2,i,j) )
c**** the following computes the snow cover as it is used in RAD_DRV.f
        scove = pearth *
     *       ( afb(i,j)*fr_snow_rad_ij(1,i,j)
     *       + (1.-afb(i,j))*fr_snow_rad_ij(2,i,j) )

        !if (snowe(i,j).gt.0.) scove=pearth
        aj(j,j_rsnow,itearth)=aj(j,j_rsnow,itearth)+scove
        areg(jr,j_rsnow)=areg(jr,j_rsnow)+scove*dxyp(j)
        aij(i,j,ij_rsnw)=aij(i,j,ij_rsnw)+scove
        aij(i,j,ij_snow)=aij(i,j,ij_snow)+snow*pearth
        aij(i,j,ij_rsit)=aij(i,j,ij_rsit)+scove

        aj(j,j_wtr1,itearth)=aj(j,j_wtr1,itearth)+wtr1*pearth
        aj(j,j_ace1,itearth)=aj(j,j_ace1,itearth)+ace1*pearth
        aj(j,j_wtr2,itearth)=aj(j,j_wtr2,itearth)+wtr2*pearth
        aj(j,j_ace2,itearth)=aj(j,j_ace2,itearth)+ace2*pearth
        aj(j,j_snow,itearth)=aj(j,j_snow,itearth)+snow*pearth
        areg(jr,j_snow)=areg(jr,j_snow)+snow*pearth*dxyp(j)
        areg(jr,j_wtr1)=areg(jr,j_wtr1)+wtr1*pearth*dxyp(j)
        areg(jr,j_ace1)=areg(jr,j_ace1)+ace1*pearth*dxyp(j)
        areg(jr,j_wtr2)=areg(jr,j_wtr2)+wtr2*pearth*dxyp(j)
        areg(jr,j_ace2)=areg(jr,j_ace2)+ace2*pearth*dxyp(j)

        aij(i,j,ij_f0e)  =aij(i,j,ij_f0e)  +f0dt+enrgp
        aij(i,j,ij_gwtr) =aij(i,j,ij_gwtr)+(wtr1+ace1+wtr2+ace2)
        aij(i,j,ij_evape)=aij(i,j,ij_evape)+evap
        do k=1,4
          aij(i,j,ij_g01+k-1)=aij(i,j,ij_g01+k-1)+wbare(k,i,j)
          aij(i,j,ij_g07+k-1)=aij(i,j,ij_g07+k-1)+wvege(k-1,i,j)
        end do
        aij(i,j,ij_g28)=aij(i,j,ij_g28)+snowbv(1,i,j)
        aij(i,j,ij_g29)=aij(i,j,ij_g29)+snowbv(2,i,j)
      end if
c****
      end do
      end do
      end subroutine ground_e

      subroutine conserv_wtg(waterg)
!@sum  conserv_wtg calculates zonal ground water incl snow
!@auth Gavin Schmidt
!@ver  1.0
      use constant, only : rhow
      use model_com, only : fim,fearth
      use geom, only : imaxj
      use ghycom, only : ngm,wbare,wvege,snowbv
      use veg_com, only : afb
      implicit none
!@var waterg zonal ground water (kg/m^2)
      real*8, dimension(grid%j_strt:grid%j_stop) :: waterg
      integer i,j,n
      real*8 wij,fb

C****	define local grid
      integer :: J_0, J_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      do j=J_0,J_1
        waterg(j)=0
        do i=1,imaxj(j)
          if (fearth(i,j).gt.0) then
            fb=afb(i,j)
            wij=fb*snowbv(1,i,j)+(1.-fb)*(wvege(0,i,j)+snowbv(2,i,j))
            do n=1,ngm
              wij=wij+fb*wbare(n,i,j)+(1.-fb)*wvege(n,i,j)
            end do
            waterg(j)=waterg(j)+fearth(i,j)*wij*rhow
          end if
        end do
      end do
      if (HAVE_SOUTH_POLE) waterg(1) =fim*waterg(1)
      if (HAVE_NORTH_POLE) waterg(jm)=fim*waterg(jm)
c****
      end subroutine conserv_wtg

      subroutine conserv_htg(heatg)
!@sum  conserv_htg calculates zonal ground energy incl. snow energy
!@auth Gavin Schmidt
!@ver  1.0
      use model_com, only : fim,fearth
      use geom, only : imaxj, dxyp
      use ghycom, only : ngm,htbare,htvege,fr_snow_ij,nsn_ij,hsn_ij
      use veg_com, only : afb
      implicit none
!@var heatg zonal ground heat (J/m^2)
      real*8, dimension(grid%j_strt:grid%j_stop) :: heatg
      integer i,j
      real*8 hij,fb,fv

C****	define local grid
      integer J_0, J_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      do j=J_0,J_1
        heatg(j)=0
        do i=1,imaxj(j)
          if (fearth(i,j).le.0) cycle
          fb=afb(i,j)
          fv=(1.d0-fb)
          hij=fb*sum( htbare(1:ngm,i,j) )
     &       +  fv*sum( htvege(0:ngm,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( hsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( hsn_ij(1:nsn_ij(2,i,j),2,i,j))
          heatg(j)=heatg(j)+fearth(i,j)*hij
        end do
      end do
      if (HAVE_SOUTH_POLE) heatg(1) =fim*heatg(1)
      if (HAVE_NORTH_POLE) heatg(jm)=fim*heatg(jm)
c****
ccc debugging ...
ccc      print *,'conserv_htg energy ',
ccc     &     sum(heatg(1:jm)*dxyp(1:jm))/(sum(dxyp(1:jm))*im)
      end subroutine conserv_htg


      end module soil_drv

      subroutine check_ghy_conservation( flag )
ccc debugging program: cam be put at the beginning and at the end
ccc of the 'surface' to check water conservation
      use constant, only : rhow
      use geom, only : imaxj
      use model_com, only : im,jm,fearth
      use DOMAIN_DECOMP, only : GRID, GET
      use fluxes, only : prec,evapor,runoe
      use ghycom, only : ngm,wbare,wvege,htbare,htvege,snowbv,dz_ij
      use veg_com, only : afb
      implicit none
      integer flag
      real*8 total_water(im,jm), error_water
      real*8, save :: old_total_water(im,jm)
!      real*8 total_energy(im,jm), error_energy
!      real*8, save :: old_total_energy(im,jm)
      integer i,j,n
      real*8 fb,fv
ccc enrgy check not implemented yet ...

C****	define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      do j=J_0,J_1
        do i=1,imaxj(j)
          if ( fearth(i,j) <= 0.d0 ) cycle

ccc just checking ...
          do n = 1,ngm
            if ( dz_ij(i,j,n) .le. 0.d0 )
     &           call stop_model('incompatible dz',255)
          enddo

          fb = afb(i,j)
          fv = 1.d0 - fb
          total_water(i,j) = fb*sum( wbare(1:ngm,i,j) )
     &         + fv*sum( wvege(0:ngm,i,j) )
     &         + fb*snowbv(1,i,j) + fv*snowbv(2,i,j)
        end do
      end do

      ! call stop_model('just testing...',255)

      if ( flag == 0 ) then
        old_total_water(:,:) = total_water(:,:)
        return
      endif

      do j=J_0,J_1
        do i=1,imaxj(j)

          !print *,'fearth = ', i, j, fearth(i,j)

          if ( fearth(i,j) <= 0.d0 ) cycle
          fb = afb(i,j)
          fv = 1.d0 - fb
          error_water = ( total_water(i,j) - old_total_water(i,j) )*rhow
     &         - prec(i,j) + evapor(i,j,4) + runoe(i,j)

          !print *, 'err H2O: ', i, j, error_water

  !        if ( abs( error_water ) > 1.d-9 ) print *, 'error'
          if ( abs( error_water ) > 1.d-9 ) call stop_model(  ! was -15
     &         'check_ghy_conservation: water conservation problem',255)

        end do
      end do

      end subroutine check_ghy_conservation

