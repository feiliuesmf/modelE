#include "rundeck_opts.h"

      module soil_drv
!@sum soil_drv contains variables and routines for the ground
!@+   hydrology driver
!@auth I. Alienov/F. Abramopolous
      use model_com, only : im,jm
      implicit none
      private
ccc   save

      public daily_earth, ground_e, init_gh, earth, conserv_wtg
     $     ,conserv_htg

      real*8 cosday,sinday
!@var gdeep keeps average (2:n) values of temperature, water and ice
      real*8, dimension(im,jm,3) :: gdeep

      real*8 spgsn !@var specific gravity of snow


      contains

      subroutine earth (ns,moddsf,moddd)
!@sum EARTH calculates surface fluxes of sensible heat,
!@+   evaporation, thermal radiation, and momentum drag over earth.
!@auth I. Alienov/F. Abramopolous
c****
      use constant, only : grav,rgas,lhe,lhs
     *     ,sha,tf,rhow,deltx
      use model_com, only : t,p,q,dtsrc,nisurf,dsig
     *     ,jday,jhour,nday,itime,jeq,fearth,modrd,itearth
      use geom, only : imaxj,dxyp,bydxyp
      use dynamics, only : pk,pek,pedn,pdsig,am,byam
      use somtq_com, only : mz
      use radncb, only : trhr,fsf,cosz1
      use sle001
     &    , only : reth,retp,advnc,evap_limits,
     &    ngm,
     &    pr,htpr,prs,htprs,w,ht,snowd,tp,fice,
     &    fv,fb,atrg,ashg,alhg,
     &    abetad,abetav,abetat,
     &    abetap,abetab,abeta,
     &    acna,acnc,
     &    aevapw,aevapd,aevapb,
     &    aruns,arunu,aeruns,aerunu,
     &    aepc,aepb,aepp,afhg,af0dt,af1dt,zw,tbcs,
     &    qm1,q1,qs,
     &    pres,rho,ts,vsm,ch,srht,trht,zs, !cd,snht,
     &           z1,
     &    eddy,
     &    nlsn,isn,nsn,dzsn,wsn,hsn,fr_snow

#ifdef TRACERS_ON
      use tracer_com, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE,trname
#endif
      use tracer_diag_com, only : taijn,tij_surf
#ifdef TRACERS_WATER
     *     ,tij_evap,tij_grnd,jls_source,tajls
#endif
      use fluxes, only : trsource
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer,trsrfflx,trprec
#endif
#endif
      use fluxes, only : dth1,dq1,uflux1,vflux1,e0,e1,evapor,prec,eprec
     *     ,runoe,erunoe,gtemp
      use ghycom, only : wbare,wvege,htbare,htvege,snowbv,
     &     nsn_ij,isn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij,
     *     snowe,tearth,wearth,aiearth,
     &     evap_max_ij, fr_sat_ij, qg_ij
#ifdef TRACERS_WATER
     *     ,trvege,trbare,trsnowbv
#endif
      USE SOCPBL, only : zgs,dtsurf              ! global
     &     ,zs1,tgv,tkv,qg_sat,hemi,pole     ! rest local
     &     ,us,vs,ws,wsh,tsv,qsrf,psi,dbl    ! ,edvisc=>kms
     &     ,khs,ppbl,ug,vg,wg,zmix   ! ,kq=>kqs
      use pblcom, only : ipbl,cmgs,chgs,cqgs,tsavg,qsavg
      use pbl_drv, only : pbl, evap_max,fr_sat
#ifdef TRACERS_WATER
     *     ,tr_evap_max
#endif
      use dagcom , only : aij,tsfrez,tdiurn,aj,areg,adiurn,jreg,
     *     ij_rune, ij_arunu, ij_pevap, ij_shdt, ij_beta, ij_trnfp0,
     *     ij_srtr, ij_neth, ij_ws, ij_ts, ij_us, ij_vs, ij_taus,
     *     ij_tauus, ij_tauvs, ij_qs, ij_tg1, ij_evap, j_trhdt, j_shdt,
     *     j_evhdt,j_evap,j_erun,j_run,j_tsrf,j_type,j_tg1,j_tg2,ij_g05
     *     ,ij_g06,ij_g11,ij_g12,ij_g13,ij_g14,ij_g15,ij_g16,ij_g17
     *     ,ij_g18,ij_g19,ij_g20,ij_g21,ij_g22,ij_g23,ij_g24,ij_g25
     *     ,ij_g26,ij_g27,ijdd,idd_ts,idd_tg1,idd_qs,idd_qg,idd_swg
     *     ,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg,idd_wg,idd_us
     *     ,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq,idd_eds,idd_dbl
     *     ,idd_ev

      implicit none

      integer, intent(in) :: ns,moddsf,moddd
      integer i,j,l,kr,jr,itype,ih
      real*8 shdt,qsats,evap,evhdt,tg2av,ace2av,trhdt,rcdmws,rcdhws,dhgs
     *     ,cdq,cdm,cdh,elhx,tg,srheat,tg1,ptype,trheat,wtr2av
     *     ,wfc1,rhosrf,ma1,tfs,th1,thv1,p1k,psk,ps,pij,psoil,pearth
     *     ,warmer,brun0,berun0,bts,bevhdt,brunu,berunu
     *     ,bshdt,btrhdt,timez,spring,zs1co

      real*8, dimension(im,jm) :: prcss
      common /workls/prcss
!@var rhosrf0 estimated surface air density
      real*8 rhosrf0

#ifdef TRACERS_ON
C**** Tracer input/output common block for PBL
!@var trsfac, trconstflx factors in surface flux boundary cond.
!@var ntx number of tracers that need pbl calculation
      real*8, dimension(ntm) :: trtop,trs,trsfac,trconstflx
      integer n,nx,ntx
      integer, dimension(ntm) :: ntix
      common /trspec/trtop,trs,trsfac,trconstflx,ntx,ntix
      real*8 totflux
      integer nsrc
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: trgrnd,trsoil_tot,tevapw,tevapd,
     *     tevapb,trruns,trrunu,trpr,trsoil_rat
      real*8, dimension(ntm,0:ngm,2) :: trw
      real*8, dimension(ntm,2) :: trsnowd
      real*8  tdp, tdt1, wsoil_tot, frac, tevap, fracvl
      integer ibv
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
c**** snowi     ocean ice snow amount (kg/m**2)
c**** snowe     earth snow amount (kg/m**2)
c**** tsi(1:2)  ocean ice temperature of first/second layer (c)
c****        4  earth temperature of first layer (c)
c****        5  earth water of first layer (kg/m**2)
c****        6  earth ice of first layer (kg/m**2)
c**** snowli      land ice snow amount (kg/m**2)
c**** tlandi(1:2) land ice temperature of first/second layer (c)
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

      dtsurf=dtsrc/nisurf
      zs1co=.5*dsig(1)*rgas/grav

      spring=-1.
      if((jday.ge.32).and.(jday.le.212)) spring=1.
      ih=1+jhour
c****
c**** outside loop over time steps, executed nisurf times every hour
c****
      timez=jday+(mod(itime,nday)+(ns-1.)/nisurf)/nday ! -1 ??
      if(jday.le.31) timez=timez+365.
c
      aregij = 0.
c****
c**** outside loop over j and i, executed once for each grid point
c****

C$OMP  PARALLEL DO PRIVATE
C$OMP*  (ACE2AV, ELHX,EVAP,EVHDT, CDM,CDH,CDQ,
C$OMP*   I,ITYPE, J, KR, L,MA1,PIJ,PSK,PEARTH,PSOIL,PS,P1K,PTYPE, QG,
C$OMP*   QG_NSAT,QSATS, RHOSRF,RHOSRF0,RCDMWS,RCDHWS, SRHDT,SRHEAT,SHDT,
C$OMP*   TRHEAT, TH1,TFS,THV1,TG1,TG,TRHDT,TG2AV, WARMER,WFC1,WTR2AV
#ifdef TRACERS_ON
C$OMP*   ,trtop,trs,trsfac,trconstflx,n,nx,ntx,ntix,totflux,nsrc
#ifdef TRACERS_WATER
C$OMP*   ,trgrnd,trsoil_tot,tevapw,tevapd,tevapb,trruns,trrunu,trpr,
C$OMP*   trsoil_rat,trw,trsnowd,tdp, tdt1, wsoil_tot,frac,tevap,ibv
#endif
#endif
C$OMP*   )
C$OMP*   SCHEDULE(DYNAMIC,2)

      loop_j: do j=1,jm
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
#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          nx=nx+1
          ntix(nx) = n
C**** Calculate first layer tracer concentration
          trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
        end if
      end do
      ntx = nx
#ifdef TRACERS_WATER
      TEVAPW=0. ; TEVAPD=0. ; TEVAPB=0. ; TRRUNS=0. ; TRRUNU=0.
#endif
#endif
c**** new quantities to be zeroed out over ground timesteps
         aruns=0.
         arunu=0.
         aeruns=0.
         aerunu=0.
      aevapw=0.
      aevapd=0.
      aevapb=0.
         alhg=0.
         aepc=0.
         aepb=0.
         afhg=0.
         atrg=0.
      ashg=0.
         af0dt=0.
         af1dt=0.
c****
c**** earth
c****
      if (pearth.le.0.) then
        ipbl(i,j,4)=0
        cycle loop_i
      endif

      itype=4
      ptype=pearth
      pr=prec(i,j)/(dtsrc*rhow)
      prs=prcss(i,j)/(dtsrc*rhow)
      htpr=eprec(i,j)/dtsrc
!!! insert htprs here
      w(1:ngm,1) =  wbare(1:ngm,i,j)
      w(0:ngm,2) =  wvege(0:ngm,i,j)
      ht(0:ngm,1) = htbare(0:ngm,i,j)
      ht(0:ngm,2) = htvege(0:ngm,i,j)
      snowd(1:2)  = snowbv(1:2,i,j)
ccc extracting snow variables
      nsn(1:2)          = nsn_ij    (1:2, i, j)
      isn(1:2)          = isn_ij    (1:2, i, j)
      dzsn(1:nlsn, 1:2) = dzsn_ij   (1:nlsn, 1:2, i, j)
      wsn(1:nlsn, 1:2)  = wsn_ij    (1:nlsn, 1:2, i, j)
      hsn(1:nlsn, 1:2)  = hsn_ij    (1:nlsn, 1:2, i, j)
      fr_snow(1:2)      = fr_snow_ij(1:2, i, j)
      ! call ghinij (i,j,wfc1) !! -seems that one doesn't need it here
      !call reth
      !call retp
c     call hydra
c??   snow = snowd(1)*fb + snowd(2)*fv
      !tg1=tbcs
      tg1 = tearth(i,j)
      srheat=fsf(itype,i,j)*cosz1(i,j)
      !srhdts=srhdts+srheat*dtsurf*ptype
! need this
      srhdt=srheat*dtsurf
c****
c**** boundary layer interaction
c****
      zs1=zs1co*tkv*pij/ps
c**** loop over ground time steps
      tg=tg1+tf
      elhx=lhe
      if(tg1.lt.0.)  elhx=lhs
      qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step
      qg = qg_ij(i,j)
      tgv=tg*(1.+qg*deltx)
      rhosrf0=100.*ps/(rgas*tgv) ! estimated surface density
#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
#ifdef TRACERS_WATER
C**** Quick and dirty calculation of water tracer amounts in
C**** soils to ensure conservation. Should be replaced with proper
C**** calculation at some point
C**** Calculate mean tracer ratio
      trsoil_tot = 0 ; wsoil_tot = 0
      do ibv=1,2
        frac=fb
        if (ibv.eq.2) frac=fv
        wsoil_tot=wsoil_tot+snowd(ibv)*frac
        do nx=1,ntx
          n=ntix(nx)
          if (tr_wd_TYPE(n).eq.nWATER) THEN
            trsnowd(nx,ibv) = TRSNOWBV(n,ibv,i,j)
            trsoil_tot(nx)=trsoil_tot(nx)+trsnowd(nx,ibv)*frac
          end if
        end do
        do l= 2-ibv,ngm
          wsoil_tot=wsoil_tot+w(l,ibv)*frac
          do nx=1,ntx
            n=ntix(nx)
            if (tr_wd_TYPE(n).eq.nWATER) THEN
              if (ibv.eq.1) then
                trw(nx,l,ibv)= TRBARE(n,l,i,j)
              else
                trw(nx,l,ibv)= TRVEGE(n,l,i,j)
              end if
              trsoil_tot(nx)=trsoil_tot(nx)+trw(nx,l,ibv)*frac
            end if
          end do
        end do
      end do
C**** calculate new tracer ratio after precip
      do nx=1,ntx
        n=ntix(nx)
        if (tr_wd_TYPE(n).eq.nWATER) THEN
          trpr(nx)=(trprec(n,i,j)*bydxyp(j))/dtsrc ! kg/m^2 s
          trsoil_rat(nx)=(trsoil_tot(nx)+dtsurf*trpr(nx))
     *         /(rhow*(wsoil_tot+dtsurf*pr))
        end if
      end do
#endif
C****
      do nx=1,ntx
        n=ntix(nx)
#ifdef TRACERS_WATER
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
C**** The select is used to distinguish water from gases or particle
        select case (tr_wd_TYPE(n))
        case (nWATER)
C**** no fractionation from ground (yet)
          trgrnd(nx)=gtracer(n,itype,i,j)*QG
C**** trsfac and trconstflx are multiplied by cq*wsh in PBL
          trsfac(nx)=1.
          trconstflx(nx)=trgrnd(nx)
        case (nGAS, nPART)
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
          trsfac(nx)=0.
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
          rhosrf0=100.*ps/(rgas*tgv) ! estimated surface density
          totflux=0.
          do nsrc=1,ntsurfsrc(n)
            totflux = totflux+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux/(dxyp(j)*rhosrf0)
#ifdef TRACERS_WATER
        end select
#endif
      end do
#endif
c***********************************************************************
c***
ccc just for test - setting all soil to saturated
      !fr_sat = 1.
      !evap_max = 1.
      !call evap_limits( .false., evap_max, fr_sat )
ccc actually PBL needs evap (kg/m^2*s) / rho_air
      evap_max = evap_max_ij(i,j) * 1000.d0 / rhosrf0
      fr_sat = fr_sat_ij(i,j)
      !print *, "evap lim. = ", evap_max, fr_sat
#ifdef TRACERS_WATER
c**** water tracers are also flux limited
      tr_evap_max(1:ntx) = evap_max * trsoil_rat(1:ntx)
#endif
      call pbl(i,j,itype,ptype)
c****
      cdm = cmgs(i,j,itype)
      cdh = chgs(i,j,itype)
      cdq = cqgs(i,j,itype)
c***********************************************************************
c**** calculate qs
c      dhgs=(zmix-zgs)*cdh*wsh
c      dqgs=(zmix-zgs)*cdq*wsh
      qs=qsrf
      ts=tsv/(1.+qs*deltx)
c**** calculate rhosrf*cdm*ws
      rhosrf=100.*ps/(rgas*tsv)
      rcdmws=cdm*ws*rhosrf
      rcdhws=cdh*wsh*rhosrf
c**** calculate fluxes of sensible heat, latent heat, thermal
c****   radiation, and conduction heat (watts/m**2)
c      snht=-sha*rcdhws*(ts-tg)  ! -not used
      trheat=trhr(0,i,j)
c **********************************************************************
c *****
c  define extra variables to be passed in surfc:
      pres  =ps
      rho   =rhosrf
      vsm   =ws
      ch    =cdh
      srht  =srheat
      trht  =trheat
      zs    =zgs  !!! will not need after qsbal is replaced
      z1    =zmix  !!! will not need after qsbal is replaced
      eddy  =khs   !!! will not need after qsbal is replaced
c     tspass=ts
c **********************************************************************
c *****
c**** calculate ground fluxes
c     call qsbal
      call ghinij (i,j,wfc1)
      call advnc
      call evap_limits( .false., evap_max_ij(i,j), fr_sat_ij(i,j) )
      tg1=tbcs
      !qg_ij(i,j) = qs  !!! - this seemed to work ok
      !! trying more precise value for qg :  qsat(tg1+tf,elhx,ps)
      qg_sat = qsat(tg1+tf,elhx,ps) ! saturated soil
      qg_nsat = qs                  ! non-sat soil, no evap
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

#ifdef TRACERS_WATER
C**** reset tracer variables
      do nx=1,ntx
        n=ntix(nx)
        if (tr_wd_TYPE(n).eq.nWATER) THEN
c**** fix outputs to mean ratio (TO BE REPLACED BY WITHIN SOIL TRACERS)
        trruns(nx)=aruns * trsoil_rat(nx) ! kg/m^2
        trrunu(nx)=arunu * trsoil_rat(nx)
#ifdef TRACERS_SPECIAL_O18
        tevapw(nx)=aevapw * trsoil_rat(nx)*fracvl(tp(1,1),trname(n))
#else
        tevapw(nx)=aevapw * trsoil_rat(nx)
#endif
        tevapd(nx)=aevapd * trsoil_rat(nx)
        tevapb(nx)=aevapb * trsoil_rat(nx)
c**** update ratio
        trsoil_rat(nx)=(trsoil_tot(nx)+dtsurf*trpr(nx)-
     *       (tevapw(nx)+tevapd(nx)+tevapb(nx)+trruns(nx)+trrunu(nx)))/
     *       (rhow*(wsoil_tot+dtsurf*pr)-(aevapw+aevapd+aevapb+aruns
     *       +arunu))
        trbare(n,1:ngm,i,j) = trsoil_rat(nx)*w(1:ngm,1)*rhow
        trvege(n,:,i,j) = trsoil_rat(nx)*w(:,2)*rhow
        trsnowbv(n,:,i,j)=trsoil_rat(nx)*snowd(:)*rhow
c       trbare(n,:,i,j) = trw(nx,:,1)
c       trvege(n,:,i,j) = trw(nx,:,2)
c       trsnowbv(n,:,i,j) = trsnowd(nx,:)
        gtracer(n,itype,i,j)=trsoil_rat(nx)
        end if
      end do
#endif
ccc copy snow variables back to storage
      nsn_ij    (1:2, i, j)         = nsn(1:2)
      isn_ij    (1:2, i, j)         = isn(1:2)
      dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
      wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
      hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
      fr_snow_ij(1:2, i, j)         = fr_snow(1:2)

      aij(i,j,ij_g18)=aij(i,j,ij_g18)+aevapb
      aij(i,j,ij_g19)=aij(i,j,ij_g19)+aevapd
      aij(i,j,ij_g20)=aij(i,j,ij_g20)+aevapw
      aij(i,j,ij_g05)=aij(i,j,ij_g05)+abetab/nisurf
      aij(i,j,ij_g06)=aij(i,j,ij_g06)+abetap/nisurf
      aij(i,j,ij_g11)=aij(i,j,ij_g11)+abeta/nisurf
      aij(i,j,ij_g12)=aij(i,j,ij_g12)+acna/nisurf
      aij(i,j,ij_g13)=aij(i,j,ij_g13)+acnc/nisurf
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
      uflux1(i,j)=uflux1(i,j)+ptype*rcdmws*us
      vflux1(i,j)=vflux1(i,j)+ptype*rcdmws*vs
c**** accumulate surface fluxes and prognostic and diagnostic quantities
      evap=aevapw+aevapd+aevapb
      evapor(i,j,4)=evapor(i,j,4)+evap
      evhdt=-alhg
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
C****
C**** Calculate Water Tracer Evaporation
C****
      do nx=1,ntx
        n=ntix(nx)
        if (tr_wd_TYPE(n).eq.nWATER) THEN
          tevap = tevapw(nx)+tevapd(nx)+tevapb(nx)
          tdp = tevap*dxyp(j)*ptype
          tdt1 = trsrfflx(i,j,n)*dtsurf
          if (trm(i,j,1,n)+tdt1+tdp.lt.0.and.tdp.lt.0) then
            write(99,*) "limiting trdew earth",i,j,n,tdp,trm(i,j,1,n)
            tevap=- (trm(i,j,1,n)+tdt1)/(dxyp(j)*ptype)
            trsrfflx(i,j,n)= - trm(i,j,1,n)/dtsurf
          else
            trsrfflx(i,j,n)=trsrfflx(i,j,n)+tdp/dtsurf
          end if
          trevapor(n,itype,i,j)=trevapor(n,itype,i,j)+tevap
          trunoe(n,i,j) = trunoe(n,i,j)+trruns(nx)+trrunu(nx)
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
#ifdef TRACERS_WATER
          taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *         trevapor(n,itype,i,j)*ptype
          taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     *         gtracer(n,1,i,j)*ptype
          tajls(j,1,jls_source(1,n))=tajls(j,1,jls_source(1,n))
     *         +trevapor(n,itype,i,j)*ptype
#endif
        end if
      end do
#endif
      aij(i,j,ij_rune)=aij(i,j,ij_rune)+aruns
      aij(i,j,ij_arunu)=aij(i,j,ij_arunu)+arunu
      aij(i,j,ij_pevap)=aij(i,j,ij_pevap)+(aepc+aepb)

      if ( warmer >= 0 ) then
        if(ts.lt.tf) tsfrez(i,j,1)=timez
        tsfrez(i,j,2)=timez
      else
        if ( tsfrez(i,j,2)+.03 >= timez .and. ts >= tf )
     $       tsfrez(i,j,2)=timez
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
        aij(i,j,ij_taus)=aij(i,j,ij_taus)+rcdmws*ws*ptype
        aij(i,j,ij_tauus)=aij(i,j,ij_tauus)+rcdmws*us*ptype
        aij(i,j,ij_tauvs)=aij(i,j,ij_tauvs)+rcdmws*vs*ptype
        aij(i,j,ij_qs)=aij(i,j,ij_qs)+qs*ptype
        aij(i,j,ij_tg1)=aij(i,j,ij_tg1)+tg1*ptype
chyd       aij(i,j,ij_arunu)=aij(i,j,ij_arunu)
chyd      *  +   (40.6*psoil+.72*(2.*(tss-tfs)-(qsatss-qss)*lhe/sha))
c**** quantities accumulated hourly for diagDD
      endif
      if ( moddd == 0 ) then
        do kr=1,4
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
C$OMP  END PARALLEL DO
C
C
      DO 825 J=1,JM
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

      subroutine init_gh(dtsurf,redogh,inisnow,istart)
c**** modifications needed for split of bare soils into 2 types
      use filemanager
      use constant, only : twopi,rhow,edpery,sha,shw_const=>shw,
     *     shi_const=>shi,lhe,lhm
      use model_com, only : fearth,vdata,itime,nday,jeq
      use dagcom, only : npts,icon_wtg,icon_htg
      use sle001
#ifdef TRACERS_WATER
      use tracer_com, only : ntm,tr_wd_TYPE,nwater,itime_tr0
      use fluxes, only : gtracer
#endif
      use fluxes, only : gtemp
      use ghycom
      implicit none

      real*8, intent(in) :: dtsurf
      integer, intent(in) :: istart
      logical, intent(in) :: redogh, inisnow
      integer iu_soil,iu_top_index,iu_veg
      integer jday
      real*8 snowdp,wtr1,wtr2,ace1,ace2,tg1,tg2
      logical :: qcon(npts)
      integer i, j, k
      real*8 one,wfc1
      real*8 dif,frdn,frup,pearth,phase,scs0,scsim,scsre,sfv,sla0
      real*8 slim,slre,svh,z
      integer iv, l
#ifdef TRACERS_WATER
      real*8 trsoil_tot,wsoil_tot
#endif

      real*8, parameter :: alamax(8) =
     $     (/ 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,10.0d0,8.0d0,4.5d0/)
      real*8, parameter :: alamin(8) =
     $     (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 8.0d0,6.0d0,1.0d0/)
      real*8, parameter :: aroot(8) =
     $     (/ 12.5d0, 0.9d0, 0.8d0,0.25d0,0.25d0,0.25d0,1.1d0,0.9d0/)
      real*8, parameter :: broot(8) =
     $     (/  1.0d0, 0.9d0, 0.4d0,2.00d0,2.00d0,2.00d0,0.4d0,0.9d0/)
      real*8, parameter :: rsar(8) =
     $     (/100d0, 100d0, 200d0, 200d0, 200d0, 300d0,250d0, 125d0/)
      real*8, parameter :: vhght(8) =
     $     (/0.1d0, 1.5d0,   5d0,  15d0,  20d0,  30d0, 25d0,1.75d0/)
      integer, parameter :: laday(8) =
     $     (/ 196,  196,  196,  196,  196,  196,  196,  196/)

c****             tundr grass shrub trees decid evrgr rainf crops
c****
c**** laday(veg type, lat belt) = day of peak lai
c old peak lai:  2nd line is for latitudes < 23.5 deg
c****    1  temperate latitudes
c****    2  non-temperate latitudes
c     data  laday/ 196,  196,  196,  196,  196,  196,  105,  196/
c     data  laday/ 196,  288,  288,  288,  288,  196,  105,  288/
c****
c**** contents of ala(k,i,j),  lai coefficients
c****   1  average leaf area index
c****   2  real amplitude of leaf area index
c****   3  imaginary amplitude of leaf area index
c****
c**** contents of acs(k,i,j),  cs coefficients
c****   1  average stomatal conductance
c****   2  real amplitude of stomatal conductance
c****   3  imaginary amplitude of stomatal conductance
c****
c**** contents of sdata(i,j,k):
c****       1 -   ngm   dz(ngm)
c****   ngm+1 - 6*ngm   q(is,ngm)
c**** 6*ngm+1 - 11*ngm   qk(is,ngm)
c**** 11*ngm+1           sl

c**** set conservation diagnostics for ground water mass and energy
      qcon=(/ .false., .false., .false., .false., .false., .true.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,"GRND WTR","(kg/m^2)        ",
     *     "(10^-9 kg/s/m^2)",1d0,1d9,icon_wtg)
      qcon=(/ .false., .false., .false., .false., .false., .true.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,"GRND ENG","(10**6 J/m^2)   ",
     *     "(10^-3 J/s/m^2) ",1d-6,1d3,icon_htg)

c read in vegetation data set: vdata
      call openunit("VEG",iu_VEG,.true.,.true.)
      do k=1,11
        call readt (iu_VEG,0,vdata(1,1,K),im*jm,vdata(1,1,k),1)
      end do
      call closeunit(iu_VEG)
      if (istart.le.0) return
c read soils parameters
      call openunit("SOIL",iu_SOIL,.true.,.true.)
      call dread (iu_SOIL,dz_ij,im*jm*(11*ngm+1),dz_ij)
      call closeunit (iu_SOIL)
ccc read topmodel parameters
      call openunit("TOP_INDEX",iu_TOP_INDEX,.true.,.true.)
      call readt(iu_TOP_INDEX,0,top_index_ij,im*jm,top_index_ij,1)
      call closeunit (iu_TOP_INDEX)
C
      one=1.
c****
c**** initialize constants
c****
c**** time step for ground hydrology
      dt=dtsurf
c**** units are mks
c**** water quantities are density times usual values in mks
c**** to get volumetric units
c**** 1m water = 1000 kg m-2; 1m3 water = 1000 kg
c fsn is the heat of fusion
      fsn= lhm * rhow
c elh is the heat of vaporization
      elh= lhe * rhow
c the sh's are the specific heat capacaties
      shw= shw_const * rhow
      shi= shi_const * rhow
c      sha= sha_const
c      shv=1911.
c the alam's are the heat conductivities
      alamw=.573345d0
      alami=2.1762d0
      alama=.025d0
ccc alamsn is not used
c      alamsn=0.088d0
      alambr=2.9d0
      alams(1)=8.8d0
      alams(2)=2.9d0
      alams(3)=2.9d0
      alams(4)=.25d0
c hw is the wilting point in meters
      hw=-100
c zhtb is depth for combining heat layers for stability
ccc looks like zhtb is not used
c      if(q(4,1).lt..01)then
c      zhtb=6.
c      else
c      zhtb=6.
cc     endif
c spgsn is the specifig gravity of snow
      spgsn=.1d0
c
c****
c**** initialize global arrays  ala, acs, afb, afr
c****
      ala(:,:,:)=0.
      acs(:,:,:)=0.
      afb(:,:)=0.
      afr(:,:,:)=0.
      acs(1,:,:)=.01d0

      do j=1,jm
        do i=1,im
          pearth=fearth(i,j)
c**** check whether data exist at this point. We should use a
c**** nearest neighbour approximation.
c**** In the meantime, use arbitrary point 10,40 (TEMPORARY FIX)
          if (pearth.gt.0) then
            if (sum(vdata(i,j,1:10)).eq.0) then
              print*,"No vege data: i,j=",i,j,vdata(i,j,1:10)
              vdata(i,j,1:11)=vdata(10,40,1:11)
c**** or these
              snowe(i,j)=0.
              tearth(i,j)=tearth(10,40)
              wearth(i,j)=wearth(10,40)
              aiearth(i,j)=aiearth(10,40)
            end if
            if (top_index_ij(i,j).eq.-1.) then
              print*,"No topo data: i,j=",i,j,top_index_ij(i,j)
              top_index_ij(i,j)=top_index_ij(10,40)
            end if
            if (sum(dz_ij(i,j,1:ngm)).eq.0
     &                .or. wbare(1,i,j) < 1.d-10) then
              print*,"No soil data: i,j=",i,j,dz_ij(i,j,1:ngm),wbare(1,i
     *             ,j)
              dz_ij(i,j,1:ngm)=dz_ij(10,40,1:ngm)
              q_ij(i,j,1:imt,1:ngm)=q_ij(10,40,1:imt,1:ngm)
              qk_ij(i,j,1:imt,1:ngm)=qk_ij(10,40,1:imt,1:ngm)
              sl_ij(i,j)=sl_ij(10,40)
c**** Not really sure what these should be set to.
              wbare(:,i,j) = wbare(:,10,40)
              wvege(:,i,j) = wvege(:,10,40)
              htbare(:,i,j)= htbare(:,10,40)
              htvege(:,i,j)= htvege(:,10,40)
              snowbv(:,i,j)= 0.
c**** or these
              NSN_IJ(:,i,j)=0. ; ISN_IJ(:,i,j)=0 ; DZSN_IJ(:,:,i,j) =0.
              WSN_IJ(:,:,i,j)=0. ; HSN_IJ(:,:,i,j)=0.
              FR_SNOW_IJ(:,i,j)=0
            end if
          end if
          afb(i,j)=vdata(i,j,1)+vdata(i,j,10)
          if(afb(i,j).gt..999) afb(i,j)=1.
          if(pearth.le.0..or.afb(i,j).ge.1.) cycle
c**** calculate lai, cs coefficicents
          sfv=0.
          sla0=0.
          slre=0.
          slim=0.
          scs0=0.
          scsre=0.
          scsim=0.
          svh=0.
          do iv=1,8
            phase=twopi*laday(iv)/365.
            if(j.lt.jeq) phase=phase+twopi/2.
            fv=vdata(i,j,iv+1)
            sfv=sfv+fv
            svh=svh+fv*vhght(iv)
            dif=(alamax(iv) - alamin(iv))
            sla0=sla0+fv*(alamax(iv) + alamin(iv))
            slre=slre+fv*dif*cos(phase)
            slim=slim+fv*dif*sin(phase)
            scs0=scs0+fv*(alamax(iv) + alamin(iv))/rsar(iv)
            scsre=scsre+fv*dif*cos(phase)/rsar(iv)
            scsim=scsim+fv*dif*sin(phase)/rsar(iv)
          end do
          ala(1,i,j)=.5/sfv*sla0
          ala(2,i,j)=.5/sfv*slre
          ala(3,i,j)=.5/sfv*slim
          acs(1,i,j)=.5/sfv*scs0
          acs(2,i,j)=.5/sfv*scsre
          acs(3,i,j)=.5/sfv*scsim
          avh(i,j)=svh/sfv
c**** calculate root fraction afr averaged over vegetation types
          do n=1,ngm
            dz(n)=dz_ij(i,j,n)
            if(dz(n).le.0.) go to 320
          end do
 320      n=n-1
          do iv=1,8
            fv=vdata(i,j,iv+1)
            z=0.
            frup=0.
            do l=1,n
              z=z+dz(l)
              frdn=aroot(iv)*z**broot(iv)
              frdn=min(frdn,one)
              if(l.eq.n)frdn=1.
              afr(l,i,j) = afr(l,i,j) + fv*(frdn-frup)
              frup=frdn
            end do
          end do
          do l=1,n
            afr(l,i,j) = afr(l,i,j)/(1.-afb(i,j))
          end do
        end do
      end do
c****
      print *,' '
      print *,'soils parameters'
      sdstnc=100.
      print *,'interstream distance (m) sdstnc:',sdstnc
      c1=90.
      print *,'canopy conductance related parameter c1:',c1
      prfr=.1d0
      print *,'fraction (by area) of precipitation prfr:',prfr
      print *,' '
      call hl0
c****
c code transplanted from subroutine input
c**** recompute ground hydrology data if necessary (new soils data)
      if (redogh) then
        jday=1+mod(itime/nday,365)
        cosday=cos(twopi/edpery*jday)
        sinday=sin(twopi/edpery*jday)

        do j=1,jm
        do i=1,im
        pearth=fearth(i,j)
        if(pearth.le.0.) then

          wbare(:,i,j)=0.
          wvege(:,i,j)=0.
          htbare(:,i,j)=0.
          htvege(:,i,j)=0.
          snowbv(:,i,j)=0.

        else
ccc??? remove next 5 lines? -check the old version
           w(1:ngm,1) =   wbare(1:ngm,i,j)
           w(0:ngm,2) =   wvege(0:ngm,i,j)
           ht(0:ngm,1) = htbare(0:ngm,i,j)
           ht(0:ngm,2) = htvege(0:ngm,i,j)
           snowd(1:2) =  snowbv(1:2,i,j)

c****     compute soil heat capacity and ground water saturation gws
          call ghinij (i,j,wfc1)
c****     fill in soils common blocks
          snowdp=snowe(i,j)/rhow
          wtr1=wearth(i,j)
          ace1=aiearth(i,j)
          tg1 =tearth(i,j)
          wtr2=wtr1
          ace2=ace1
          tg2 =tg1
c          wtr2=gdata(i,j,9)   ! this cannot be right
c          ace2=gdata(i,j,10)
c          tg2 =gdata(i,j,8)
          call ghinht (snowdp, tg1,tg2, wtr1,wtr2, ace1,ace2)

c****     copy soils prognostic quantities to model variables
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
      do j=1,jm
        do i=1,im
          if (fearth(i,j).gt.0) then
            gtemp(1,4,i,j)=tearth(i,j)
#ifdef TRACERS_WATER
C**** Quick and dirty calculation of water tracer amounts to calculate
C**** gtracer. Should be replaced with proper calculation at some point
C**** Calculate mean tracer ratio
            fb=afb(i,j) ; fv=1.-fb
            wsoil_tot = 0
            wsoil_tot=wsoil_tot+snowbv(1,i,j)*fb+snowbv(2,i,j)*fv+
     *           sum(wbare(1:ngm,i,j))*fb+sum(wvege(0:ngm,i,j))*fv
            trsoil_tot = 0
            do n=1,ntm
              if(itime_tr0(n).le.itime) then
                trsoil_tot=trsoil_tot+trsnowbv(n,1,i,j)*fb
     *               +trsnowbv(n,2,i,j)*fv+sum(trbare(n,1:ngm,i,j))*fb
     *               +sum(trvege(n,0:ngm,i,j))*fv
                gtracer(n,4,i,j)=trsoil_tot/(rhow*wsoil_tot)
              end if
            end do
#endif
          end if
        end do
      end do

ccc   some extra code from snowmodel ghinit
      so_%rosmp = 8.  ! no idea what this number means, but it is used
                      ! in computation of runoff

ccc   init snow here
ccc hope this is the right place to split first layer into soil
ccc and snow  and to set snow arrays
ccc!!! this should be done only when restarting from an old
ccc!!! restart file (without snow model data)

      if (inisnow) then
        do j=1,jm
        do i=1,im
          pearth=fearth(i,j)
          if(pearth.le.0.) then
            nsn_ij(:,i,j)     = 0
            isn_ij(:,i,j)     = 0
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

            call ghinij (i,j,wfc1)
            call set_snow

            nsn_ij    (1:2, i, j)         = nsn(1:2)
            isn_ij    (1:2, i, j)         = isn(1:2)
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

      return
      end subroutine init_gh

      subroutine ghinij (i0,j0, wfcap)
c**** input:
c**** avh(i,j) - array of vegetation heights
c**** spgsn - specific gravity of snow
c**** output:
c**** vh - vegetation height
c**** snowm - snow masking depth
c**** wfcap - water field capacity of top soil layer, m
c****
      use snow_model, only : i_earth,j_earth
      use sle001, only : dz,qk,ngm,imt,ng,zb,zc,fr,q,sl,xklh0 !spgsn,
     *     ,fb,fv,snowm,alai,alaie,rs,prs,ijdebug,n !alaic,vh,
     *     ,thets,thetm,ws,thm,nth,shc,shw,htprs,pr !shcap,shtpr,
     *     ,htpr
     *     ,top_index
      use ghycom, only : dz_ij,sl_ij,q_ij,qk_ij,avh,afr,afb,ala,acs
     *     ,top_index_ij
      implicit none
      integer i0,j0
      real*8 wfcap
      integer l,ibv,k,i
      real*8 aa,one
      real*8 alaic,vh,shtpr
      real*8, parameter :: shcap(imt) = (/2d6,2d6,2d6,2.5d6,2.4d6/)

      one=1.
      ijdebug=i0*1000+j0
      i_earth = i0
      j_earth = j0

ccc passing topmodel parameters
      top_index = top_index_ij(i0, j0)
c**** set up layers
      dz(1:ngm)=dz_ij(i0,j0,1:ngm)
      q(1:imt,1:ngm)=q_ij(i0,j0,1:imt,1:ngm)
      qk(1:imt,1:ngm)=qk_ij(i0,j0,1:imt,1:ngm)
      sl=sl_ij(i0,j0)

      do n=1,ngm
        if(dz(n).le.0.) go to 21
      end do
   21 n=n-1
      if(n.le.0) then
         write (99,*) 'ghinij:  n <= 0:  i,j,n=',i0,j0,n,(dz(k),k=1,43)
         stop
      end if
c**** calculate the boundaries, based on the thicknesses.
      zb(1)=0.
      do l=1,n
        zb(l+1)=zb(l)-dz(l)
      end do
c**** calculate the layer centers, based on the boundaries.
      do l=1,n
        zc(l)=.5*(zb(l)+zb(l+1))
      end do
c**** fr: root fraction in layer l  (1=fr(1)+fr(2)+...+fr(n))
      do l=1,n
        fr(l)=afr(l,i0,j0)
      end do
c**** vh: vegetation height
      vh=avh(i0,j0)
      snowm=vh*spgsn
c**** fb,fv: bare, vegetated fraction (1=fb+fv)
      fb=afb(i0,j0)
      fv=1.-fb
c**** alai: leaf area index
      alai=ala(1,i0,j0)+cosday*ala(2,i0,j0)+sinday*ala(3,i0,j0)
      alai=max(alai,one)
      alaic=5.0
      alaie=alaic*(1.-exp(-alai/alaic))
c**** rs: minimum stomatal resistance
      rs=alai/(acs(1,i0,j0)+cosday*acs(2,i0,j0)+sinday*acs(3,i0,j0))
c???  cnc=alai/rs   redefined before being used (qsbal,cond)
c
cw    write(6,*)'n=',n,'  r=',r
cw    write(6,91)
cw 91 format(1x,5x,'zb',5x,'zc',5x,'dz'/1x,21('-'))
cw    do 95 l=1,n
cw 95 write(6,100)zb(l),zc(l),dz(l)
cw    write(6,100)zb(n+1)
cw100 format(1x,3f7.3)
cw    write(6,*)
c****
      do ibv=1,2
        do l=1,n
          thets(l,ibv)=0.
          thetm(l,ibv)=0.
          do i=1,imt-1
            thets(l,ibv)=thets(l,ibv)+q(i,l)*thm(0,i)
            thetm(l,ibv)=thetm(l,ibv)+q(i,l)*thm(nth,i)
          end do
          ws(l,ibv)=thets(l,ibv)*dz(l)
        end do
      end do
      ws(0,2)=.0001d0*alai
      wfcap=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
c****
      call xklh0
c****
      do ibv=1,2
        do l=1,n
          shc(l,ibv)=0.
          do i=1,imt
            shc(l,ibv)=shc(l,ibv)+q(i,l)*shcap(i)
          end do
          shc(l,ibv)=(1.-thets(l,ibv))*shc(l,ibv)*dz(l)
        end do
      end do
c****
c shc(0,2) is the heat capacity of the canopy
      aa=ala(1,i0,j0)
      shc(0,2)=(.010d0+.002d0*aa+.001d0*aa**2)*shw
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
      integer l, ibv, ll

      wfc1=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
      wfc2=0.
      fbv=fb
      do 30 ibv=1,2
      do 20 l=2,n
      wfc2=wfc2+fbv*ws(l,ibv)
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
        do l=2,n
          fice(l,ibv)=ace2/(wtr2+ace2+1.d-20)
          wmin=thetm(l,ibv)*dz(l)
          wet2=(wtr2+ace2)/(wfc2+1.d-20)
          if(wet2.gt.1.) wet2=1.
          w(l,ibv)=wmin+(ws(l,ibv)-wmin)*wet2
          tp(l,ibv)=tg2
        end do
      end do
c****
      entry ghexht
c****
c**** compute ht (heat w/m+2)
      do ibv=1,2
        ll=2-ibv
        do l=ll,n
          if(tp(l,ibv)) 2,4,6
 2        ht(l,ibv)=tp(l,ibv)*(shc(l,ibv)+w(l,ibv)*shi)-w(l,ibv)*fsn
          cycle
 4        ht(l,ibv)=-fice(l,ibv)*w(l,ibv)*fsn
          cycle
 6        ht(l,ibv)=tp(l,ibv)*(shc(l,ibv)+w(l,ibv)*shw)
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
      integer l, ibv
      tg2av=0.
      wtr2av=0.
      ace2av=0.
      do 3500 ibv=1,2
      wc=0.
      htc=0.
      shcc=0.
      do l=2,n
        wc=wc+w(l,ibv)
        htc=htc+ht(l,ibv)
        shcc=shcc+shc(l,ibv)
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
      integer i,j,k
!@var subr identifies where check was called from
      character*6, intent(in) :: subr

c**** check for nan/inf in earth data
      call check3(wbare ,ngm  ,im,jm,subr,'wb')
      call check3(wvege ,ngm+1,im,jm,subr,'wv')
      call check3(htbare,ngm+1,im,jm,subr,'hb')
      call check3(htvege,ngm+1,im,jm,subr,'hv')
      call check3(snowbv,2    ,im,jm,subr,'sn')

c**** check for reasonable temperatures over earth
      x=1.001
      do j=1,jm
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
      use constant, only : rhow,twopi,edpery,tf
      use model_com, only : nday,nisurf,jday,fearth,wfcs
      use geom, only : imaxj
      use dagcom, only : aij,tdiurn,ij_strngts,ij_dtgdts,ij_tmaxe
     *     ,ij_tdsl,ij_tmnmx,ij_tdcomp
      use ghycom, only : snoage
      implicit none
      real*8 tsavg,wfc1
      integer i,j,itype
      logical, intent(in) :: end_of_day
c****
c**** find leaf-area index & water field capacity for ground layer 1
c****
      cosday=cos(twopi/edpery*jday)
      sinday=sin(twopi/edpery*jday)
      do j=1,jm
        do i=1,im
          wfcs(i,j)=24.
          if (fearth(i,j).gt.0.) then
            call ghinij(i,j,wfc1)
            wfcs(i,j)=rhow*wfc1 ! canopy part changes
          end if
        end do
      end do

      if (end_of_day) then
c****
c**** increase snow age each day (independent of ts)
c****
        do j=1,jm
          do i=1,imaxj(j)
            do itype=1,3
              snoage(itype,i,j)=1.+.98d0*snoage(itype,i,j)
            end do
            tsavg=tdiurn(i,j,5)/(nday*nisurf)
            if(32.+1.8*tsavg.lt.65.)
     *           aij(i,j,ij_strngts)=aij(i,j,ij_strngts)+(33.-1.8*tsavg)
            aij(i,j,ij_dtgdts)=aij(i,j,ij_dtgdts)+18.*((tdiurn(i,j,2)-
     *           tdiurn(i,j,1))/(tdiurn(i,j,4)-tdiurn(i,j,3)+1.d-20)-1.)
            aij(i,j,ij_tdsl)=aij(i,j,ij_tdsl)+
     *           (tdiurn(i,j,4)-tdiurn(i,j,3))
            aij(i,j,ij_tdcomp)=aij(i,j,ij_tdcomp)+
     *           (tdiurn(i,j,6)-tdiurn(i,j,9))
            aij(i,j,ij_tmaxe)=aij(i,j,ij_tmaxe)+(tdiurn(i,j,4)-tf)
            if (tdiurn(i,j,6).lt.aij(i,j,ij_tmnmx))
     *           aij(i,j,ij_tmnmx)=tdiurn(i,j,6)
          end do
        end do
      end if

      return
      end subroutine daily_earth

      subroutine ground_e
!@sum  ground_e driver for applying surface fluxes to land fraction
!@auth original development team
!@ver  1.0
      use model_com, only : fearth,itearth
      use geom, only : imaxj,dxyp
      use ghycom, only : snowe, tearth,wearth,aiearth,wbare,wvege,snowbv
     *     ,fr_snow_ij,afb
      use dagcom, only : aj,areg,aij,jreg,ij_evap,ij_f0e,ij_evape
     *     ,ij_gwtr,ij_tg1,j_wtr1,j_ace1,j_wtr2,j_ace2
     *     ,j_snow,j_evap,j_type,ij_g01,ij_g07,ij_g28
     *     ,ij_g29,j_rsnow,ij_rsnw,ij_rsit,ij_snow
      use fluxes, only : e0,e1,evapor,eprec
      implicit none

      real*8 snow,tg1,tg2,f0dt,f1dt,evap,wtr1,wtr2,ace1,ace2
     *     ,pearth,enrgp,scove
      integer i,j,jr,k

      do j=1,jm
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
        scove = pearth *
     *       ( afb(i,j)*fr_snow_ij(1,i,j)
     *       + (1.-afb(i,j))*fr_snow_ij(2,i,j) )
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
c        aj(j,j_tg1 ,itearth)=aj(j,j_tg1, itearth)+tg1 *pearth
c        aj(j,j_tg2 ,itearth)=aj(j,j_tg2, itearth)+tg2 *pearth
c        aj(j,j_type,itearth)=aj(j,j_type,itearth)+     pearth
        aj(j,j_snow,itearth)=aj(j,j_snow,itearth)+snow*pearth
c        aj(j,j_evap,itearth)=aj(j,j_evap,itearth)+evap*pearth
c        areg(jr,j_tg1) =areg(jr,j_tg1) +tg1 *pearth*dxyp(j)
c        areg(jr,j_tg2) =areg(jr,j_tg2) +tg2 *pearth*dxyp(j)
        areg(jr,j_snow)=areg(jr,j_snow)+snow*pearth*dxyp(j)
        areg(jr,j_wtr1)=areg(jr,j_wtr1)+wtr1*pearth*dxyp(j)
        areg(jr,j_ace1)=areg(jr,j_ace1)+ace1*pearth*dxyp(j)
        areg(jr,j_wtr2)=areg(jr,j_wtr2)+wtr2*pearth*dxyp(j)
        areg(jr,j_ace2)=areg(jr,j_ace2)+ace2*pearth*dxyp(j)

        aij(i,j,ij_f0e)  =aij(i,j,ij_f0e)  +f0dt+enrgp
c        aij(i,j,ij_tg1)  =aij(i,j,ij_tg1)  +tg1 *pearth
        aij(i,j,ij_gwtr) =aij(i,j,ij_gwtr)+(wtr1+ace1+wtr2+ace2)
c        aij(i,j,ij_evap) =aij(i,j,ij_evap) +evap*pearth
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
!@sum  conserv_wtg calculates zonal ground water
!@auth Gavin Schmidt
!@ver  1.0
      use constant, only : rhow
      use model_com, only : fim,fearth
      use geom, only : imaxj
      use sle001, only : ngm
      use ghycom, only : wbare,wvege,afb
      implicit none
!@var waterg zonal ground water (kg/m^2)
      real*8, dimension(jm) :: waterg
      integer i,j,n
      real*8 wij,fb

      do j=1,jm
        waterg(j)=0
        do i=1,imaxj(j)
          if (fearth(i,j).gt.0) then
            fb=afb(i,j)
            wij=(1.-fb)*wvege(0,i,j)
            do n=1,ngm
              wij=wij+fb*wbare(n,i,j)+(1.-fb)*wvege(n,i,j)
            end do
            waterg(j)=waterg(j)+fearth(i,j)*wij*rhow
          end if
        end do
      end do
      waterg(1) =fim*waterg(1)
      waterg(jm)=fim*waterg(jm)
c****
      end subroutine conserv_wtg

      subroutine conserv_htg(heatg)
!@sum  conserv_htg calculates zonal ground energy
!@auth Gavin Schmidt
!@ver  1.0
      use model_com, only : fim,fearth
      use geom, only : imaxj
      use sle001, only : ngm
      use ghycom, only : htbare,htvege,afb
      implicit none
!@var heatg zonal ground heat (J/m^2)
      real*8, dimension(jm) :: heatg
      integer i,j,n
      real*8 hij,fb

      do j=1,jm
        heatg(j)=0
        do i=1,imaxj(j)
          if (fearth(i,j).gt.0) then
            fb=afb(i,j)
            hij=0.
            do n=0,ngm
              hij=hij+fb*htbare(n,i,j)+(1.-fb)*htvege(n,i,j)
            end do
            heatg(j)=heatg(j)+fearth(i,j)*hij
          end if
        end do
      end do
      heatg(1) =fim*heatg(1)
      heatg(jm)=fim*heatg(jm)
c****
      end subroutine conserv_htg

      end module soil_drv
