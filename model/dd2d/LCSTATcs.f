#include "rundeck_opts.h"

c
c initial sketch of cubed-sphere versions of routines DIAGB,DIAG5A,DIAG7A
c which compute eddy and spectral statistics around latitude circles.
c 
c M. Kelley 12/2008
c

      module diag_latcirc
!@sum diag_latcirc a module containing info needed by DIAGB,DIAG5A,DIAG7A
!@auth M. Kelley
      use diag_com, only : imlon=>im,jmlat=>jm,imh
      use cs2ll_utils
      implicit none
      save

      type(cs2ll_type) :: cs2ll

      integer :: jlat_eq,jlat_50n,jlat_45n ! temporarily here

      real*8, dimension(:,:,:), allocatable :: rh,pvort ! temporary

      end module diag_latcirc

      subroutine init_diag_latcirc()
!@sum init_diag_latcirc allocate/initialize contents of module diag_latcirc
!@auth M. Kelley
      use geom, only : lon,lat
      use domain_decomp_1d, only : grid
      use cs2ll_utils
      use diag_latcirc
      call init_cs2ll_type(grid%dd2d,imlon,jmlat,cs2ll,lon,lat)
      return
      end subroutine init_diag_latcirc

      SUBROUTINE DIAGB
!@sum DIAGB calculate constant-pressure latitude-circle diagnostics.
!@+   This version is for a cubed sphere grid.
!@auth M. Kelley
      use constant, only : lhe,sha,tf,teeny
      use model_com, only : mdiag,lm,ls1,sig,sige,dsig
     &     ,u,v,t,p,q,wm
      use geom, only : dxyp
      use diag_latcirc
      use diag_com, only :
     &     ajk=>ajk_loc,speca,nspher,klayer,
     &     jk_dpa,jk_dpb,jk_temp,jk_hght,jk_q,jk_theta,
     &     jk_rh,jk_u,jk_v,jk_zmfke,jk_totke,jk_zmfntsh,
     &     jk_totntsh,jk_zmfntgeo,jk_totntgeo,jk_zmfntlh,
     &     jk_totntlh,jk_zmfntke,jk_totntke,jk_zmfntmom,jk_totntmom,
     &     jk_dpsqr,jk_nptsavg,
     &     jk_vvel,jk_zmfvtdse,jk_totvtdse,jk_zmfvtlh,jk_totvtlh,
     &     jk_vtgeoeddy,jk_potvort,jk_vtpv,
     &     jk_vtpveddy,jk_nptsavg1,jk_totvtke,
     &     jk_vtameddy,jk_totvtam,jk_eddvtpt,jk_cldh2o
      use dynamics, only : phi,plij,sd,pmid,pedn
      use diag_loc, only : w,tx
      use domain_decomp_1d, only : get, grid, am_i_root
      implicit none
      real*8, dimension(imh+1,nspher) :: ke,ke_part
      real*8, dimension(imh+1) :: xu,xv,xke
      real*8, dimension(jmlat,lm) ::
     &     stjk,dpjk,ujk,vjk,wjk,tjk,wtjk,uvjk,wujk

      integer :: i,l,lcp,l1cp,l2cp,ldn,lup,mnow,kq,nq
      integer :: jj,jlat,jcalc,khemi,kspher,kspher_eq,kspher_45n
      integer, dimension(lm) :: klayer_eq,klayer_45n
      integer, parameter :: nqmax=30
      real*8, dimension(lm,nqmax), target :: lonsums
      real*8, dimension(:), pointer ::
     &     dpsum,dpusum,dpvsum,dptxsum,dpqsum,dpthsum,dpphisum,
     &     dpkesum,dppvortsum,nlesum,wlesum,txlesum,ulesum,
     &     philesum,thlesum,qlesum,pvortlesum,kelesum
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     dp,dpu,dpv,dptx,dpq,dprh,dpcw,dpth,dpphi,dpke,dppvort,
     &     wle,txle,ule,vle,phile,thle,qle,pvortle,kele
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     ull,vll,txll,thll,qll,rhll,wmll,phill,pvortll,sdll
      real*8, dimension(cs2ll%isd:cs2ll%ied) :: psll
      real*8, dimension(lm) :: nlesum_loc,bydpsum
      real*8, dimension(imlon,lm) :: usqrtdp,vsqrtdp
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     l1fm,l2fm,lwdn,lwup,lwdn2,lwup2

      real*8 :: dpl,wtdn,wtup,bydp,vi,kedp,uzm,vzm,wzm,nrat
     &     ,bysqrtdp
      real*8 :: dpuvi,dptvi,dpqvi,fimi,dpsqi,dpphivi,dprhi,dpcwi
     &     ,dpkei,dpkevi
     &     ,wui,wti,wthi,wzi,wqi,wkei,wpvorti

c
c set up pointers that facilitate longitudinal sums across PEs
c
      kq=0
c layer midpoint quantities
      kq=kq+1; dpsum => lonsums(:,kq)
      kq=kq+1; dpusum => lonsums(:,kq)
      kq=kq+1; dpvsum => lonsums(:,kq)
      kq=kq+1; dptxsum => lonsums(:,kq)
      kq=kq+1; dpqsum => lonsums(:,kq)
      kq=kq+1; dpthsum => lonsums(:,kq)
      kq=kq+1; dpphisum => lonsums(:,kq)
      kq=kq+1; dpkesum => lonsums(:,kq)
      kq=kq+1; dppvortsum => lonsums(:,kq)
c edge quantities
      kq=kq+1; nlesum => lonsums(:,kq)
      kq=kq+1; wlesum => lonsums(:,kq)
      kq=kq+1; txlesum => lonsums(:,kq)
      kq=kq+1; ulesum => lonsums(:,kq)
      kq=kq+1; philesum => lonsums(:,kq)
      kq=kq+1; thlesum => lonsums(:,kq)
      kq=kq+1; qlesum => lonsums(:,kq)
      kq=kq+1; pvortlesum => lonsums(:,kq)
      kq=kq+1; kelesum => lonsums(:,kq)
      nq = kq

c
c some arrays to make the kspher logic more robust
c
      klayer_eq(:) = klayer(:) + 2
      klayer_45n(:) = klayer(:) + 3

      ke_part(:,:)=0.
c
c loop over latitudes
c
      do jj=1,jmlat

      jlat = cs2ll%jlat_sched(jlat)
      if(cs2ll%ni(jlat).eq.0) cycle

c
c horizontal interpolation to this latitude circle
c
      call interp_xy_3D_jlat(p,psll)
      call interp_xy_3D_jlat(u,ull)
      call interp_xy_3D_jlat(v,vll)
      call interp_xy_3D_jlat(tx,txll)
      call interp_xy_3D_jlat(t,thll)
      call interp_xy_3D_jlat(q,qll)
      call interp_xy_3D_jlat(rh,rhll)
      call interp_xy_3D_jlat(wm,wmll)
      call interp_xy_3D_jlat(phi,phill)
      call interp_xy_3D_jlat(sd,sdll)
      call interp_xy_3D_jlat(pvort,pvortll)

c
c vertical interpolation to the midpoint of cp layers
c
      do lcp=1,lm
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        dp(i,lcp) = 0d0
        dpu(i,lcp) = 0d0
        dpv(i,lcp) = 0d0
        dptx(i,lcp) = 0d0
        dpq(i,lcp) = 0d0
        dpcw(i,lcp) = 0d0
        dprh(i,lcp) = 0d0
        dpth(i,lcp) = 0d0
        dpphi(i,lcp) = 0d0
        dpke(i,lcp) = 0d0
        dppvort(i,lcp) = 0d0
      enddo
      enddo
      call calc_l1l2fm(l1fm,l2fm)
      do l=1,lm
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        if(l.lt.ls1) then
          l1cp = l1fm(i,l)
          l2cp = l2fm(i,l)
          dpl = (l1fm(i,l)-l1cp)*psll(i)*dsig(l)
        else
          l1cp = l
          l2cp = l
          dpl = dpcp(l)
        endif
        do lcp=l1cp,l2cp
          dp(i,lcp) = dp(i,lcp) + dpl
          dpu(i,lcp) = dpu(i,lcp) + dpl*ull(i,l)
          dpv(i,lcp) = dpv(i,lcp) + dpl*vll(i,l)
          dptx(i,lcp) = dptx(i,lcp) + dpl*txll(i,l)
          dpq(i,lcp) = dpq(i,lcp) + dpl*qll(i,l)
          dprh(i,lcp) = dpcw(i,lcp) + dpl*rhll(i,l)
          dpcw(i,lcp) = dpcw(i,lcp) + dpl*wmll(i,l)
          dpth(i,lcp) = dpth(i,lcp) + dpl*thll(i,l)
          dpphi(i,lcp) = dpphi(i,lcp) + dpl*phill(i,l)
          dpke(i,lcp) = dpke(i,lcp) +
     &         dpl*(ull(i,l)**2+vll(i,l)**2)
          dppvort(i,lcp) = dppvort(i,lcp) + dpl*pvortll(i,l)
          if(lcp.eq.l2cp) then
            dpl = (l2fm(i,l)-l2cp)*psll(i)*dsig(l)
          elseif(lcp.lt.l2cp) then
            dpl = dpcp(lcp+1)
          endif
        enddo
      enddo
      enddo

c
c vertical interpolation to cp layer edges
c if a layer edge is underground, wtdn==wtup==0
c
      call calc_lw(lwdn,lwup,lwdn2,lwup2)
      do l=1,lm-1
      nlesum(l) = 0d0
      wlesum(l) = 0d0
      txlesum(l) = 0d0
      ulesum(l) = 0d0
      philesum(l) = 0d0
      thlesum(l) = 0d0
      qlesum(l) = 0d0
      pvortlesum(l) = 0d0
      kelesum(l) = 0d0
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        ldn = lwdn(i,l)
        lup = lwup(i,l)
        wtdn = lwdn(i,l)-ldn
        wtup = lwup(i,l)-lup
        ule(i,l) = wtdn*ull(i,ldn)+wtup*ull(i,lup)
        vle(i,l) = wtdn*vll(i,ldn)+wtup*vll(i,lup)
        txle(i,l) = wtdn*txll(i,ldn)+wtup*txll(i,lup)
        thle(i,l) = wtdn*thll(i,ldn)+wtup*thll(i,lup)
        qle(i,l) = wtdn*qll(i,ldn)+wtup*qll(i,lup)
        phile(i,l) = wtdn*phill(i,ldn)+wtup*phill(i,lup)
        pvortle(i,l) = wtdn*pvortll(i,ldn)+wtup*pvortll(i,lup)
        kele(i,l) = .5*(ule(i,l)**2 + vle(i,l)**2)
c use different weights for wle
        ldn = lwdn2(i,l)
        lup = lwup2(i,l)
        wtdn = lwdn2(i,l)-ldn
        wtup = lwup2(i,l)-lup
        wle(i,l) = wtdn*sdll(i,ldn)+wtup*sdll(i,lup)
c partial longitudinal sums of layer edge quantities
        if(wtdn.ne.0. .and. wtup.ne.0.) nlesum(l) = nlesum(l) + 1d0
        wlesum(l) = wlesum(l) + wle(i,l)
        txlesum(l) = txlesum(l) + txle(i,l)
        ulesum(l) = ulesum(l) + ule(i,l)
        philesum(l) = philesum(l) + phile(i,l)
        thlesum(l) = thlesum(l) + thle(i,l)
        qlesum(l) = qlesum(l) + qle(i,l)
        pvortlesum(l) = pvortlesum(l) + pvortle(i,l)
        kelesum(l) = kelesum(l) + kele(i,l)
      enddo
      enddo

c
c partial longitudinal sums of layer midpoint quantities
c
      do l=1,lm
        dpsum(l) = 0d0
        dpusum(l) = 0d0
        dpvsum(l) = 0d0
        dptxsum(l) = 0d0
        dpqsum(l) = 0d0
        dpthsum(l) = 0d0
        dpphisum(l) = 0d0
        dpkesum(l) = 0d0
        dppvortsum(l) = 0d0
        dpuvi = 0d0
        dptvi = 0d0
        dpqvi = 0d0
        dpphivi = 0d0
        fimi = 0d0
        dpsqi = 0d0
        dprhi = 0d0
        dpcwi = 0d0
        dpkei = 0d0
        dpkevi = 0d0
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          dpsum(l) = dpsum(l) + dp(i,l)
          dpusum(l) = dpusum(l) + dpu(i,l)
          dpvsum(l) = dpvsum(l) + dpv(i,l)
          dptxsum(l) = dptxsum(l) + dptx(i,l)
          dpqsum(l) = dpqsum(l) + dpq(i,l)
          dpthsum(l) = dpthsum(l) + dpth(i,l)
          dpphisum(l) = dpphisum(l) + dpphi(i,l)
          dpkesum(l) = dpkesum(l) + dpke(i,l)
          dppvortsum(l) = dppvortsum(l) + dppvort(i,l)
          if(dp(i,l).gt.0d0) fimi = fimi + 1d0
          bydp = 1d0/(dp(i,l)+teeny)
          vi = dpv(i,l)*bydp
          dpsqi = dpsqi + dp(i,l)*dp(i,l)
          dpuvi = dpuvi + dpu(i,l)*vi
          dptvi = dptvi + dptx(i,l)*vi
          dpqvi = dpqvi + dpq(i,l)*vi
          dpphivi = dpphivi + dpphi(i,l)*vi
          dprhi = dprhi + dprh(i,l)
          dpcwi = dpcwi + dpcw(i,l)
          kedp = (dpu(i,l)**2+dpv(i,l)**2)*bydp
          dpkei = dpkei + kedp
          dpkevi = dpkevi + kedp*vi
        enddo
        ajk(jlat,l,jk_nptsavg)=ajk(jlat,l,jk_nptsavg)+fimi
        ajk(jlat,l,jk_nptsavg1)=ajk(jlat,l,jk_nptsavg1)+fimi
        ajk(jlat,l,jk_dpsqr)=ajk(jlat,l,jk_dpsqr)+dpsqi
        ajk(jlat,l,jk_rh)=ajk(jlat,l,jk_rh)+dprhi
        ajk(jlat,l,jk_cldh2o)=ajk(jlat,l,jk_cldh2o)+dpcwi
        ajk(jlat,l,jk_totntmom) = ajk(jlat,l,jk_totntmom) + dpuvi
        ajk(jlat,l,jk_totntsh) = ajk(jlat,l,jk_totntsh) + dptvi
        ajk(jlat,l,jk_totke)=ajk(jlat,l,jk_totke)+dpkei
        ajk(jlat,l,jk_totntke)=ajk(jlat,l,jk_totntke)+dpkevi
        ajk(jlat,l,jk_totntgeo)=ajk(jlat,l,jk_totntgeo)+dpphivi
        ajk(jlat,l,jk_totntlh)=ajk(jlat,l,jk_totntlh)+dpqvi
c save certain sums for later
        dpjk(jlat,l) = dpsum(l)
        uvjk(jlat,l) = dpuvi
      enddo

c
c save certain local sums before combining sums from different PEs
c
      nlesum_loc(:) = nlesum(:)

c
c add the contributions from different PEs to the sums around this
c latitude circle.  After this call, pointers to partial longitudinal
c sums are pointers to full longitudinal sums.
c
      call sumxpe_zonal(cs2ll,jlat,lonsums,lm*nq)

      bydpsum(:) = 1d0/(dpsum(:)+teeny)

c
c store northward fluxes from the zonal mean flow and
c certain zonal sums/means on the root PE for this jlat
c
      if(cs2ll%am_i_rootj(jlat)) then
      do l=1,lm

        uzm = dpusum(l)*bydpsum(l)
        vzm = dpvsum(l)*bydpsum(l)
        wzm = wlesum(l)/(nlesum(l)+teeny)

        ajk(jlat,l,jk_dpa) = ajk(jlat,l,jk_dpa) + dpsum(l)
        ajk(jlat,l,jk_dpb) = ajk(jlat,l,jk_dpb) + dpsum(l)
        ajk(jlat,l,jk_u)   = ajk(jlat,l,jk_u) + dpusum(l)
        ajk(jlat,l,jk_v)   = ajk(jlat,l,jk_v) + dpvsum(l)
        ajk(jlat,l,jk_temp) = ajk(jlat,l,jk_temp) + dptxsum(l)
        ajk(jlat,l,jk_hght) = ajk(jlat,l,jk_hght) + dpphisum(l)
        ajk(jlat,l,jk_q) = ajk(jlat,l,jk_q) + dpqsum(l)
        ajk(jlat,l,jk_theta) = ajk(jlat,l,jk_theta) + dpthsum(l)
        ajk(jlat,l,jk_potvort) = ajk(jlat,l,jk_potvort) + dppvortsum(l)
        ajk(jlat,l,jk_vvel) = ajk(jlat,l,jk_vvel) + wlesum(l)

        ajk(jlat,l,jk_zmfntmom) = ajk(jlat,l,jk_zmfntmom)
     &       + dpusum(l)*vzm
        ajk(jlat,l,jk_zmfntsh) = ajk(jlat,l,jk_zmfntsh)
     &       + dptxsum(l)*vzm
        ajk(jlat,l,jk_zmfntlh) = ajk(jlat,l,jk_zmfntlh)
     &       + dpqsum(l)*vzm
        ajk(jlat,l,jk_zmfntgeo) = ajk(jlat,l,jk_zmfntgeo)
     &       + dpphisum(l)*vzm

        ajk(jlat,l,jk_zmfke) = ajk(jlat,l,jk_zmfke)
     &       + dpsum(l)*(uzm**2+vzm**2)
        ajk(jlat,l,jk_zmfntke) = ajk(jlat,l,jk_zmfntke)
     &       + dpkesum(l)*vzm

        ajk(jlat,l,jk_zmfvtdse) = ajk(jlat,l,jk_zmfvtdse)
     &       + (sha*txlesum(l)+philesum(l))*wzm
        ajk(jlat,l,jk_zmfvtlh) = ajk(jlat,l,jk_zmfvtlh)
     &       + qlesum(l)*wzm

c save certain sums for later
        ujk(jlat,l) = dpusum(l)*bydpsum(l)
c        vjk(jlat,l) = dpvsum(l)*bydpsum(l)
c        tjk(jlat,l) = dpthsum(l)*bydpsum(l)
c        wjk(jlat,l) = wzm/dxyp(jlat)

      enddo
      endif

c
c partial longitudinal sums of vertical fluxes: total and eddy
c
      do l=1,lm
        wui = 0d0
        wti = 0d0
        wthi = 0d0
        wzi = 0d0
        wqi = 0d0
        wkei = 0d0
        wpvorti = 0d0
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          wui = wui + wle(i,l)*ule(i,l)
          wti = wti + wle(i,l)*txle(i,l)
          wthi = wthi + wle(i,l)*thle(i,l)
          wqi = wqi + wle(i,l)*qle(i,l)
          wzi = wzi + wle(i,l)*phile(i,l)
          wkei = wkei + wle(i,l)*kele(i,l)
          wpvorti = wpvorti + wle(i,l)*pvortle(i,l)
        enddo
        wzm = wlesum(l)/(nlesum(l)+teeny)
        nrat = nlesum_loc(l)/(nlesum(l)+teeny)

        ajk(jlat,l,jk_totvtdse) = ajk(jlat,l,jk_totvtdse)+sha*wti+wzi
        ajk(jlat,l,jk_totvtlh) = ajk(jlat,l,jk_totvtlh)+wqi
        ajk(jlat,l,jk_totvtke) = ajk(jlat,l,jk_totvtke)+wkei

        ajk(jlat,l,jk_vtpv) = ajk(jlat,l,jk_vtpv) + wpvorti
        ajk(jlat,l,jk_vtpveddy) = ajk(jlat,l,jk_vtpveddy)+
     &       wpvorti-pvortlesum(l)*wzm*nrat

        ajk(jlat,l,jk_totvtam) = ajk(jlat,l,jk_totvtam)+wui

        ajk(jlat,l,jk_vtgeoeddy) = ajk(jlat,l,jk_vtgeoeddy)
     &       +wzi-philesum(l)*wzm*nrat

c save certain sums for later
        wujk(jlat,l) = (wui-ulesum(l)*wzm*nrat)/dxyp(jlat)

        ajk(jlat,l,jk_eddvtpt) = ajk(jlat,l,jk_eddvtpt)
     &       +(wthi-thlesum(l)*wzm*nrat)/dxyp(jlat)
        ajk(jlat,l,jk_vtameddy) = ajk(jlat,l,jk_vtameddy)
     &      +wujk(jlat,l)*dxyp(jlat)

      enddo

c
c spectral analysis of KE
c
      do l=1,lm
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          bysqrtdp = 1d0/sqrt(dp(i,l)+teeny)
          dpu(i,l) = dpu(i,l)*bysqrtdp
          dpv(i,l) = dpv(i,l)*bysqrtdp
        enddo
      enddo
      call pack_zonal(cs2ll,jlat,lm,dpu,usqrtdp)
      call pack_zonal(cs2ll,jlat,lm,dpv,vsqrtdp)
      jcalc = cs2ll%jlat_calc(jlat)
      if(jcalc.gt.0) then
        if(jcalc.gt.jmlat/2) then
          khemi = 2
        else
          khemi = 1
        endif
        do l=1,lm
          kspher = klayer(l) + khemi-1
          kspher_eq = klayer_eq(l)
          kspher_45n = klayer_45n(l)
          call ffte(usqrtdp(:,l),xu)
          call ffte(vsqrtdp(:,l),xv)
          xke(:) = (xu(:)+xv(:))*dxyp(jcalc)
          ke_part(:,kspher) = ke_part(:,kspher) + xke(:)
          if(jcalc.eq.jlat_eq) then
            ke_part(:,kspher_eq) = ke_part(:,kspher_eq) + xke(:)
          elseif(jcalc.eq.jlat_45n) then
            ke_part(:,kspher_45n) = ke_part(:,kspher_45n) + xke(:)
          endif
        enddo
      endif

      enddo ! latitude loop

c
c sum up spectral analysis contributions
c
      call sumxpe(ke_part,ke)
      if(am_i_root()) speca(:,18,:)=speca(:,18,:)+ke(:,:)

c
c on the root processor, compute the diagnostics which
c include latitudinal gradients of zonal means.
c for example, the correction to E-P fluxes
c from eddy streamfunction psijk acting on du/dy,du/dz
c
c      call sumxpe(ujk)
c      call sumxpe(psijk)
c      if(am_i_root()) then
c      endif


      call timer (mnow,mdiag)

      RETURN
      END SUBROUTINE DIAGB

      SUBROUTINE DIAG5A (M5,NDT)
!@sum DIAG5A calculate KE/PE spectral diagnostics around latitude circles.
!@+   This version is for a cubed sphere grid.
!@auth M. Kelley
C****
C**** THIS DIAGNOSTICS ROUTINE PRODUCES A SPECTRAL ANALYSIS OF KINETIC
C**** AND AVAILABLE POTENTIAL ENERGIES AND THEIR TRANSFER RATES BY
C**** VARIOUS ATMOSPHERIC PROCESSES.
C****
C**** THE PARAMETER M INDICATES WHAT IS STORED IN SPECA(N,M,KSPHER),
C**** IT ALSO INDICATES WHEN DIAG5A IS BEING CALLED.
C**** M=1  MEAN STANDING KINETIC ENERGY            BEFORE SOURCES
C****   2  MEAN KINETIC ENERGY                     BEFORE DYNAMICS
C****   3  MEAN POTENTIAL ENERGY
C****   4  CONVERSION OF K.E. BY ADVECTION         AFTER ADVECTION
C****   5  CONVERSION OF K.E. BY CORIOLIS FORCE    AFTER CORIOLIS TERM
C****   6  CONVERSION FROM P.E. INTO K.E.          AFTER PRESS GRAD FORC
C****   7  CHANGE OF K.E. BY DYNAMICS              AFTER DYNAMICS
C****   8  CHANGE OF P.E. BY DYNAMICS
C****   9  CHANGE OF K.E. BY CONDENSATION          AFTER CONDENSATION
C****  10  CHANGE OF P.E. BY CONDENSATION
C****  11  CHANGE OF P.E. BY RADIATION             AFTER RADIATION
C****  12  CHANGE OF K.E. BY SURFACE               AFTER SURFACE
C****  13  CHANGE OF P.E. BY SURFACE
C****  14  CHANGE OF K.E. BY FILTER                AFTER FILTER
C****  15  CHANGE OF P.E. BY FILTER
C****  16  CHANGE OF K.E. BY DAILY                 AFTER DAILY
C****  17  CHANGE OF P.E. BY DAILY
C****  18  UNUSED
C****  19  LAST KINETIC ENERGY
C****  20  LAST POTENTIAL ENERGY
C****
      use constant, only : sha
      use model_com, only : lm,dsig,idacc,mdiag,
     &     p,ptop,sig,t,u,v,zatmo
      use geom, only : areag,dxyp,axyp,lat2d
      use diag_com, only : speca,atpe,nspher,kspeca,klayer
c      use diag_com, only : ajl=>ajl_loc,jl_ape
      use diag_loc, only : lupa,ldna
      use dynamics, only : pk,pdsig,sqrtp
      use domain_decomp_1d, only : grid,am_i_root,sumxpe,write_parallel
      use diag_latcirc

      implicit none
      integer :: m5,ndt

      real*8, dimension(imh+1) :: xu,xv,xke,xape
      real*8, dimension(imh+1,nspher) :: ke,ape,ke_part,ape_part
      real*8, dimension(2) :: tpe,tpe_part
      real*8, dimension(lm) :: thgm,gmean,thgm_part,gmean_part
      integer, dimension(kspeca), parameter ::
     &     mtpeof=(/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/)
      integer :: i,j,k,ks,kspher,l,ldn,lup,mape,mke,mnow,mtpe,n,nm
      integer :: jj,jlat,jcalc,khemi,kspher_eq,kspher_45n
      real*8 :: sumt
      integer, dimension(lm) :: klayer_eq,klayer_45n
      logical :: do_ke
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     ull,vll,tll,u_tmp,v_tmp,t_tmp
      real*8, dimension(cs2ll%isd:cs2ll%ied) :: psll,sqrtpa
      real*8, dimension(imlon,lm) :: u_il,v_il,t_il

      integer :: i_0,i_1,j_0,j_1

      j_0=grid%j_strt
      j_1=grid%j_stop
      i_0=grid%i_strt
      i_1=grid%i_stop

      nm=1+imlon/2
      jlat_45n=2.+.75*(jmlat-1.)

C****
C**** Note: KSPHER has been re-arranged from previous models to better
C****       deal with optionally resolved stratosphere. The higher
C****       values are only used if the model top is high enough.
C****
C**** KSPHER=1 SOUTHERN TROPOSPHERE         2 NORTHERN TROPOSPHERE
C****        3 EQUATORIAL TROPOSPHERE       4 45 DEG NORTH TROPOSPHERE
C****
C****        5 SOUTHERN LOW STRATOSPHERE    6 NORTHERN LOW STRATOSPHERE
C****        7 EQUATORIAL LOW STRATOSPHERE  8 45 DEG NORTH LOW STRATOSPH
C****
C****        9 SOUTHERN MID STRATOSPHERE   10 NORTHERN MID STRATOSPHERE
C****       11 EQUATORIAL MID STRATOSPHERE 12 45 DEG NORTH MID STRATOSPH
C****
C****       13 SOUTHERN UPP STRATOSPHERE   14 NORTHERN UPP STRATOSPHERE
C****       15 EQUATORIAL UPP STRATOSPHERE 16 45 DEG NORTH UPP STRATOSPH
C****
      select case (m5)
      case (1,2,7,9,12,14,16) ! both winds and temperature have changed
        do_ke = .true.
        MKE=M5 ! mke,mape are not used for cases 1,2
        MAPE=M5+1
      case (11) ! radiation does not change winds
        do_ke = .false.
        MAPE=M5
      case default
        CALL WRITE_PARALLEL(M5, UNIT=6, format=
     &       "('0INCORRECT VALUE OF M5 IN DIAG5A.  M5=',I5)")
        call stop_model('INCORRECT VALUE OF M5 IN DIAG5A.',255)
      end select

C**** CURRENT TOTAL POTENTIAL ENERGY BY HEMISPHERE
c version w/o zatmo is computed elsewhere, so why the duplication?
      tpe_part(:) = 0d0
      do j=j_0,j_1
      do i=i_0,i_1
        sumt=0d0
        do l=1,lm
          sumt = sumt + t(i,j,l)*pk(l,i,j)*pdsig(l,i,j)
        end do
        sumt=(zatmo(i,j)*(p(i,j)+ptop)+sumt*sha)*axyp(i,j)
        if(lat2d(i,j).lt.0.) then ! southern hemisphere
          tpe_part(1) = tpe_part(1) + sumt
        else                      ! northern hemisphere
          tpe_part(2) = tpe_part(2) + sumt
        endif
      enddo
      enddo
      call sumxpe(tpe_part, tpe)

c
c For APE, first calculate global means for each layer of
c pot. temp (thgm) and static stability (gmean).  For convenience,
c use native-grid thermodynamic variables.  If interpolation
c to the latlon grid is not conservative, there will be
c a small inconsistency that is less important than interpretational
c ambiguities.
c
      do l=1,lm
        ldn=ldna(l)
        lup=lupa(l)
        thgm_part(l)=0d0
        gmean_part(l)=0d0
        do j=j_0,j_1
        do i=i_0,i_1
          thgm_part(l) = thgm_part(l)+t(i,j,l)*sqrtp(i,j)*axyp(i,j)
          gmean_part(l) = gmean_part(l) +axyp(i,j)*
     *         (p(i,j)*sig(l)+ptop)*(t(i,j,lup)-t(i,j,ldn))
     *         /(p(i,j)*pk(l,i,j))
        enddo
        enddo
      enddo
      call sumxpe(thgm_part,thgm)
      call sumxpe(gmean_part,gmean)
      call esmfx_bcast(grid,thgm)
      call esmfx_bcast(grid,gmean)
      do l=1,lm
        thgm(l)=thgm(l)/areag
        ldn=ldna(l)
        lup=lupa(l)
        gmean(l)=areag*(sig(ldn)-sig(lup))/gmean(l)
      enddo

c
c some arrays to make the kspher logic more robust
c
      klayer_eq(:) = klayer(:) + 2
      klayer_45n(:) = klayer(:) + 3

c
c loop over latitudes
c
      ke_part(:,:)=0.
      ape_part(:,:)=0.

      do jj=1,jmlat

      jlat = cs2ll%jlat_sched(jlat)
      if(cs2ll%ni(jlat).eq.0) cycle

c
c interpolate p,u,v,t to this latitude circle
c
      call interp_xy_3D_jlat(p,psll)
      call interp_xy_3D_jlat(u,ull)
      call interp_xy_3D_jlat(v,vll)
      call interp_xy_3D_jlat(t,tll)

      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        sqrtpa(i) = sqrt(psll(i)*dxyp(jlat))
      enddo

c
c send wind data to root for this jlat
c
      if(do_ke) then
        do l=1,lm
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          u_tmp(i,l) = ull(i,l)*sqrtpa(i)
          v_tmp(i,l) = vll(i,l)*sqrtpa(i)
        enddo
        enddo
        call pack_zonal(cs2ll,jlat,lm,u_tmp,u_il)
        call pack_zonal(cs2ll,jlat,lm,v_tmp,v_il)
      endif

c
c send temperature data to root for this jlat
c
      do l=1,lm
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        t_tmp(i,l)=tll(i,l)*sqrtpa(i)-thgm(l)
      enddo
      enddo
      call pack_zonal(cs2ll,jlat,lm,t_tmp,t_il)
c
c if at a jlat for calculation, do the spectral analysis
c
      jcalc = cs2ll%jlat_calc(jlat)
      if(jcalc.gt.0) then
        if(jcalc.gt.jmlat/2) then
          khemi = 2
        else
          khemi = 1
        endif
        do l=1,lm
          kspher = klayer(l) + khemi-1
          kspher_eq = klayer_eq(l)
          kspher_45n = klayer_45n(l)
c kinetic energy
          if(do_ke) then
            call ffte(u_il(:,l),xu)
            call ffte(v_il(:,l),xv)
            xke(:) = (xu(:)+xv(:))*dsig(l)
            ke_part(:,kspher) = ke_part(:,kspher) + xke(:)
            if(jcalc.eq.jlat_eq) then
              ke_part(:,kspher_eq) = ke_part(:,kspher_eq) + xke(:)
            elseif(jcalc.eq.jlat_45n) then
              ke_part(:,kspher_45n) = ke_part(:,kspher_45n) + xke(:)
            endif
          endif
c potential energy
          call ffte(t_il(:,l),xape)
          xape(:) = xape(:)*dsig(l)*gmean(l)
          ape_part(:,kspher) = ape_part(:,kspher) + xape(:)
          if(jcalc.eq.jlat_eq) then
            ape_part(:,kspher_eq) = ape_part(:,kspher_eq) + xape(:)
          elseif(jcalc.eq.jlat_45n) then
            ape_part(:,kspher_45n) = ape_part(:,kspher_45n) + xape(:)
          endif
        enddo
      endif

      enddo ! end latitude loop

c
c combine partial sums onto root
c
      call sumxpe(ke_part,ke)
      call sumxpe(ape_part,ape)

c
c on root, place sums and differences into SPECA
c
      if(am_i_root()) then

        if (ndt /= 0) then
c**** energy transfer rates as differences
          if(do_ke) then
            speca(:,mke,:)=speca(:,mke,:)+(ke(:,:)-speca(:,19,:))/ndt
          endif
          speca(:,mape,:)=speca(:,mape,:)+(ape(:,:)-speca(:,20,:))/ndt
          mtpe=mtpeof(mape)
          atpe(mtpe,:)=atpe(mtpe,:)+(tpe(:)-atpe(8,:))/ndt
        end if
c**** store latest values
        if(do_ke) speca(:,19,:)=ke(:,:)
        speca(:,20,:)=ape(:,:)
        atpe(8,:)=tpe(:)

        if (m5.eq.2) then
c**** accumulate mean kinetic energy and mean potential energy
          speca(:,2,:)=speca(:,2,:)+ke(:,:)
          speca(:,3,:)=speca(:,3,:)+ape(:,:)
          atpe(1,:)=atpe(1,:)+tpe(:)
        end if

      endif ! am_i_root

      call timer (mnow,mdiag)
      return
      end subroutine diag5a

      subroutine diag7a
!@sum DIAG7A calculate wave power around selected latitude circles.
!@+   This version is for a cubed sphere grid.
!@auth M. Kelley
c****
c**** this routine accumulates a time sequence for selected
c**** quantities and from that prints a table of wave frequencies.
c****
      use constant, only : bygrav
      use model_com, only : lm,idacc,mdiag,p,u,v
      use dynamics, only : phi
      use diag_com, only : nwav_dag,wave,max12hr_sequ
     &     ,kwp,re_and_im,ia_12hr
      use diag_loc, only : ldex
      use diag_latcirc
      use domain_decomp_1d, only : sumxpe,am_i_root
      implicit none

      real*8, dimension(0:imh) :: an,bn
      integer, parameter :: km=6,kqmax=12
      integer :: nmax=nwav_dag
      real*8, dimension(imlon,km) :: htrd,htrd_loc
      real*8, dimension(imlon,lm) :: ueq,ueq_loc,veq,veq_loc
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) :: ull,vll
      real*8, dimension(cs2ll%isd:cs2ll%ied,km) :: phill,htrd_tmp
      real*8, dimension(cs2ll%isd:cs2ll%ied) :: psll
      real*8, dimension(km), parameter ::
     &     pmb=(/922.,700.,500.,300.,100.,10./),
     &     ght=(/500.,2600.,5100.,8500.,15400.,30000./)
      real*8, dimension(lm) :: p00,aml,pdsigl,pmidl
      real*8, dimension(lm+1) :: pednl
      real*8 :: slope
      integer i,jlat,idacc9,k,kq,l,mnow,n

      idacc9=idacc(ia_12hr)+1
      idacc(ia_12hr)=idacc9
      if (idacc9.gt.max12hr_sequ) return


c
c interpolate winds to equator and send to root
c
      jlat = jlat_eq
      ueq_loc=0d0
      veq_loc=0d0
      if(cs2ll%ni(jlat).gt.0) then ! this processor has valid lons at jlat
        call interp_xy_3D_jlat(u,ull,jlat)
        call interp_xy_3D_jlat(v,vll,jlat)
        call pack_zonal(cs2ll,jlat,lm,ull,ueq_loc)
        call pack_zonal(cs2ll,jlat,lm,vll,veq_loc)
      endif
      call sumxpe(ueq_loc, ueq)
      call sumxpe(veq_loc, veq)

c
c interpolate geopotential to 50 N and send to root
c
      jlat = jlat_50n
      htrd_loc(:,:)=0d0
      if(cs2ll%ni(jlat).gt.0) then ! this processor has valid lons at jlat
        call interp_xy_3D_jlat(phi,phill,jlat)
        call interp_xy_3D_jlat(p,psll,jlat)
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          call calc_vert_amp(psll(i),lm,p00,aml,pdsigl,pednl,pmidl)
          l=2
          do k=1,km
            do while(pmb(k).lt.pmidl(l) .and. l.lt.lm)
              l=l+1
            enddo
c**** assume that phi is linear in log p
            slope=(phill(i,l-1)-phill(i,l))/log(pmidl(l-1)/pmidl(l))
            htrd_tmp(i,k)=
     &           (phill(i,l)+slope*log(pmb(k)/pmidl(l)))*bygrav-ght(k)
          enddo
        enddo
        call pack_zonal(cs2ll,jlat,km,htrd_tmp,htrd_loc)
      endif
      call sumxpe(htrd_loc, htrd)

c
c do the fourier analyses on root
c
      if(am_i_root()) then
        do kq=1,3
          call fft(ueq(1,ldex(kq)),an,bn)
          do n=1,nmax
            wave(1,idacc9,n,2*kq-1)=an(n)
            wave(2,idacc9,n,2*kq-1)=bn(n)
          enddo
          call fft(veq(1,ldex(kq)),an,bn)
          do n=1,nmax
            wave(1,idacc9,n,2*kq)=an(n)
            wave(2,idacc9,n,2*kq)=bn(n)
          enddo
        enddo
        do kq=7,kqmax
          call fft(htrd(1,kq-6),an,bn)
          do n=1,nmax
            wave(1,idacc9,n,kq)=an(n)
            wave(2,idacc9,n,kq)=bn(n)
          end do
        end do
      endif

      call timer (mnow,mdiag)
      return
      end subroutine diag7a
