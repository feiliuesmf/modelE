c**** SLE001 E001M12 SOMTQ SLB211M9
c**** (same as frank's soils64+2bare_soils+old runoff)
c**** change to evap calculation to prevent negative runoff
c**** soils45 but with snowmelt subroutine snmlt changed
c**** to melt snow before 1st layer ground ice.
ccc   comments from new soils
c**** 8/11/97 - modified to include snow model
c**** 10/29/96 - aroot and broot back to original values; added
c**** call to cpars to change vegetation parameters for pilps.
c**** 9/7/96 - back to heat capacities/10
c**** 5/14/96 - added soils64 surface runoff changes/corrections
c**** 11/13/95 - changed aroot and broot for rainf for 1.5m root depth
c**** 10/11/95 - back to full heat capacity, to avoid cd oscillations.
c**** changes for pilps: (reversed)
c**** use soils100.com instead of soils45.com
c**** set im=36,jm=24
c**** set sdata,vdata and fdata to real*4
c**** divide canopy heat capacities by 10.d0
c**** change aroot of grass to 1.0d0, to reduce root depth.
c**** end of changes for pilps
c**** changes for pilps: (kept)
c**** modify gdtm surface flux timestep limits
c**** define new diagnostics
c**** zero out diagnostics at start of advnc
c**** end of changes for pilps
c**** modified for 2 types of bare soils
c****
c**** soils62 soils45 soils45          cdfxa 04/27/95
c**** same as soils45 but with snowmelt subroutine snmlt changed
c**** to melt snow before 1st layer ground ice.
ccc   end comments from new soils
c**** also corrects evaps calculation.
c**** also includes masking effects in radation fluxes.
c**** modifies timestep for canopy fluxes.
c**** soils45 10/4/93
c**** uses bedrock as a soil texture.  soil depth of 3.5m
c**** everywhere, where layers can have bedrock.
c**** requires sm693.data instead of sm691.data.
c**** sdata needs to be changed in calling program.
c**** soils44b 8/25/93
c**** uses snow conductivity of .088 w m-1 c-1 instead of .3d0
c**** soils44 8/16/93
c**** adds bedrock for heat calculations, to fill out the
c**** number of layers to ngm.
c**** soils43 6/25/93
c**** comments out call to fhlmt heat flux limits.
c**** uses ghinij to return wfc1, eliminates rewfc.
c**** soils42 6/15/93
c**** adds snow insulation
c**** soils41 5/24/93
c**** uses snow masking depth from vegetation height to determine
c**** fraction of snow that is exposed.
c**** reth must be called prior to retp.
c**** soils40 5/10/93
c**** removes snow from canopy and places it on vegetated soil.
c**** soils39 4/19/93
c**** modifications for real*8 or real*4 runs.  common block
c**** ordering changed for efficient alignment.  sdata,fdata,
c**** and vdata are explicitly real*4.  on ibm rs/6000, should
c**** be compiled with -qdpc=e option for real*8 operation.
c**** to run real*4, change implicit statement in include file.
c**** soils38 2/9/93
c**** adds heat flux correction to handle varying coefficients
c**** of drag.
c**** soils37 1/25/93
c**** changes soil crusting parameter ku/d from .05 per hour to .1d0,
c**** to agree with morin et al.
c**** soils36 11/12/92
c**** calculates heat conductivity of soils using devries method.
c**** changes loam material heat capacity and conductivity
c**** to mineral values.
c**** soils35 10/27/92
c**** includes effect of soil crusting for infiltration by
c**** modifying hydraulic conductivity calculation of layer
c**** 1 in hydra.
c**** soils34 8/28/92
c**** uses effective leaf area index alaie for purposes of
c**** canopy conductance calculation.
c**** soils33 8/9/92
c**** changes canopy water storage capacity to .1mm per lai from 1.d0
c**** soils32 7/15/92
c**** 1) for purposes of infiltration only, reduces soil conductivity
c**** by (1-thetr*fice) instead of (1-fice).
c**** 2) betad is reduced by fraction of ice in each layer.
c**** 3) transpired water is removed by betad fraction in each layer,
c**** instead of by fraction of roots. prevents negative runoff.
c**** 4) speeds up hydra by using do loop instead of if check,
c**** by using interpolation point from bisection instead of logs,
c**** and by avoiding unnecessary calls to hydra.  also elimates call
c**** to hydra in ma89ezm9.f.
c**** soils31 7/1/92
c**** 1) fixes fraction of roots when soil depth is less than root
c**** depth, thus fixing non-conservation of water.
c**** soils30 6/4/92
c**** 1) uses actual final snow depth in flux limit calculations,
c**** instead of upper and lower limits.  fixes spurious drying
c**** of first layer.
c****
      module sle001

ccc alai used in printout only

      use constant, only : stbo,tfrz=>tf,sha,lhe,one,zero
      implicit none
      save
      private

c**** public functions:
      public hl0, set_snow, reth, retp, advnc

c**** public parameters:
      public ngm,ng,imt,nlsn  ! used by ghycom

c**** public variables:
      public dz,qk,zb,zc,fr,q,sl,xklh0
     *     ,snowm,alai,alaie,rs,ijdebug,n
     *     ,thets,thetm,ws,thm,nth,shc,shw
     *     ,htpr
     *     ,top_index
      public
     &    pr,htpr,prs,htprs,w,ht,snowd,tp,fice,hour,
     &    fv,fb,atrg,ashg,alhg,
     &    abetad,abetav,abetat,
     &    abetap,abetab,abeta,
     &    acna,acnc,
     &    aevapw,aevapd,aevapb,
     &    aruns,arunu,aeruns,aerunu,
     &    adifs,aedifs,
     &    aepc,aepb,aepp,afhg,af0dt,af1dt,zw,tbcs,
     &    qm1,q1,qs,
     &    pres,rho,ts,vsm,ch,srht,trht,zs,
     &    z1,
     &    tg,t1,vg,eddy,
     &    isn,nsn,dzsn,wsn,hsn,fr_snow
      public dt, fsn, elh, shi, alamw, alami, alama,
     $     alambr, alams, hw,
     $     sdstnc, c1, prfr, so_

      real*8, external :: qsat,dqsatdt
      real*8 beta,betab,betat,betav  ! used only in qsbal (+ maybe accm)
      integer, parameter :: ngm=6, ng=ngm+1, imt=5
      integer, parameter :: igcm=0  !?? do we really need it ?
      real*8 pr,htpr,prs,htprs,w(0:ngm,2),ht(0:ngm,2)
     & ,snowd(2),ws(0:ngm,2),tp(0:ngm,2),fice(0:ngm,2),hour
     & ,dz(ngm),q(imt,ngm),qk(imt,ngm),sl,fv,fb,alai,atrg,ashg,alhg
     & ,abetad,abetav,abetat,abetap,abetab,abeta,acna,acnc,aevapw,aevapd
     & ,aevapb,aruns,arunu,aeruns,aerunu,adifs,aedifs,aepc,aepb,aepp
     & ,afhg,af0dt,af1dt,cnc,zw(2),fd,fw,fm,alaie,tbcs,tcs
     & ,snowm
      integer  ijdebug
      real*8 theta(0:ng,2),thets(0:ng,2),f(0:ng,2)
     &     ,fh(0:ng,2),zb(ng),zc(ng),shw,shi,fsn,elh,sdstnc
     &     ,shc(0:ng,2),alamw,alami,alama,alambr,alams(imt-1)
     &     ,xkh(ng,2),xkhm(ng,2),fr(ng),rnf(2),rnff(ng,2),h(0:ng,2)
     &     ,xk(0:ng,2),xinfc(2),snk(0:ng,2),snkh(0:ng,2),d(0:ng,2)
     &     ,xku(0:ng,2)
      real*8 dt,algdel,dtr
     &     ,dts,betad,hw,dr,drs,rs,prfr,c1
     &     ,pre(2),betadl(ngm)
      integer n
c     common/weight/a(4,imt-1),b(4,imt-1),p(4,imt-1)
      real*8 thetm(0:ng,2),hlm(0:64),
     &     xklm(0:64,imt-1),dlm(0:64,imt-1),thm(0:64,imt-1),alph0
      integer nth
      real*8 qm1,q1,qs,qb,qc,evap(2),evapw,evapd,evaps,epb,epc,
     &     snowf,snowfs
      !dimension fs(2),fs1(2,2)
      real*8 pres,rho,ts,vsm,ch,srht,trht,zs,z1,
     &     tg,t1,vg,eddy
      real*8 thrm(2),snsh(2),xlth(2)
      real*8 top_index, xkus(ng,2), xkusa(2)
      !dimension snowdu(2)

ccc   i hate to do it, but as a quick solution i declare snow variables
ccc   as global ones ...
      integer, parameter :: nlsn=3
      integer isn(2),nsn(2)
      real*8  dzsn(nlsn+1,2),wsn(nlsn,2),hsn(nlsn,2),tsn1(2),fr_snow(2)
      real*8  flmlt(2),fhsng(2),thrmsn(2),hesn(2),snshsn(2)
      real*8  tbs, snshs

ccc   the following looks like diagnistic output
      real*8 etbcs, esnowd, ezw, ewtr1, eice1, etp(0:ngm,2)
      real*8 asnowd,arnff(ngm),aw(0:ngm),atp(0:ngm),af(0:ng),apre,atbcs
      real*8 aevaps

ccc   put some obscure parameters to soil_params structure
      type, public :: soil_params
        real*8 rosmp
      end type soil_params

      type (soil_params) so_

ccc   derivatives of surface fluxes with respect to temperature
      real*8 snsh_dt, epb_dt, evaps_dt

      contains

      subroutine reth
c**** revises values of theta based upon w.
c**** input:
c**** w - water depth, m
c**** ws - saturated water depth, m
c**** dz - layer thickness, m
c**** thets - saturated theta
c**** snowd - snow depth, water equivalent m
c**** snowm - snow masking depth, water equivalent m
c**** output:
c**** theta - water saturation
c**** fw - fraction of wet canopy
c**** fd - fraction of dry canopy
c**** fm - fraction of snow that is exposed, or masking.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      integer lsn,ibv,l
      do ibv=1,2
        do  l=1,n
          theta(l,ibv)=w(l,ibv)/dz(l)
        end do
      end do
c**** do canopy layer
c**** here theta is the fraction of canopy covered by water
      if(ws(0,2).gt.0.d0)then
        theta(0,2)=(w(0,2)/ws(0,2))**(2.d0/3.d0)
      else
        theta(0,2)=0.d0
      endif
      theta(0,2)=min(theta(0,2),one)
c**** set up snowd variables
      do ibv=1,2
        snowd(ibv)=0.d0
        do lsn=1,nsn(ibv)
ccc    we compute snowd as if all snow was distributed uniformly
ccc    over the cell (i.e. snowd = wsn * sn_frac / 1.d0 )
          snowd(ibv) = snowd(ibv) + wsn(lsn,ibv) * fr_snow(ibv)
        enddo
      enddo
c**** fraction of wet canopy fw
      fw=theta(0,2)
c**** determine fm from snowd depth and masking depth
      fm=1.d0-exp(-snowd(2)/(snowm+1d-12))
c**** correct fraction of wet canopy by snow fraction
      fw=fw+fm*(1.d0-fw)
      fd=1.d0-fw
      return
      end subroutine reth

      subroutine hydra
c     routine to return the equlibrium value of h in a mixed soil
c     layer.  the h is such that each soil texture has the same
c     value of h, but differing values of theta.
c     hydra also calculates the conductivity xk and diffussivity d.
c**** input:
c**** theta(l,ibv) - volumetric water concentration
c**** thetm(l,ibv) - minimum theta
c**** thets(l,ibv) - maximum theta
c**** nth - number of h0 intervals in table, a power of two.
c**** hlm(j) - table of h values, from 0 (at j=0) to hmin (at j=nth)
c**** thm(j,i) - value of relative theta at hlm(j) in texture i,
c**** ranging between thets(l,ibv)  at j=0 to thetm(l,ibv) at j=nth.
c**** output:
c**** h - potential, m, including both matric and gravitational
c**** d - diffusivity, dl.
c**** xk - conductivity m/s.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c     solve for h using bisection
c     we assume that if j1.lt.j2 then hlm(j1).gt.hlm(j2)
c     and thm(j1,i).gt.thm(j2,i).
c
c     algdel=alog(1.d0+alph0)
c
      real*8 d1,d2,dl,hl,temp,thr,thr0,thr1,thr2,xk1,xkl,xklu,xku1,xku2
      real*8 xkud
      integer i,j,ibv,l,ith,j1,j2,jcm,jc
      real*8 dz_total
      xkud=2.78d-5
      jcm=nint(log(float(nth))/log(2.d0))
      do ibv=1,2
        xk(n+1,ibv)=0.0d0
        xku(0,ibv)=0.d0
        do l=1,n
          j1=0
          j2=nth
          thr1=thets(l,ibv)
          thr2=thetm(l,ibv)
          thr0=theta(l,ibv)
          thr0=min(thr1,thr0)
          thr0=max(thr2,thr0)
          do jc=1,jcm
            j=(j1+j2)/2
            thr=0.d0
            do i=1,imt-1
              thr=thr+thm(j,i)*q(i,l)
            end do
            if(thr-thr0 .lt. 0.d0) then
c     here thr is too small, bisect on low j end
              j2=j
              thr2=thr
            else if (thr-thr0 .gt. 0.d0) then
c     here thr is too large, bisect on high j end
              j1=j
              thr1=thr
            else                ! i.e. .eq.
c     here thr is equal to thr0
              hl=hlm(j)
              j1=j
              thr1=thr0
c     the strange value for thr2 below is only for calculating temp
              thr2=-10.d0
              go to 500
            end if
          end do                ! jc
c     here theta is between two adjacent thr's. interpolate.
          hl=(hlm(j1)*(thr0-thr2)+hlm(j2)*(thr1-thr0))/(thr1-thr2)
 500      continue
c**** only filling hl array with matric potential (gravitational to be
c**** added later)
          h(l,ibv)=hl
c**** calculate diffusivity
          ith=j1
          temp=(thr1-thr0)/(thr1-thr2)
          d1=0.d0
          d2=0.d0
          xku1=0.d0
          xku2=0.d0
          xkus(l,ibv) = 0.d0
          do i=1,imt-1
            d1=d1+q(i,l)*dlm(ith,i)
            d2=d2+q(i,l)*dlm(ith+1,i)
            xku1=xku1+q(i,l)*xklm(ith,i)
            xku2=xku2+q(i,l)*xklm(ith+1,i)
            xkus(l,ibv) = xkus(l,ibv) + q(i,l)*xklm(0,i)
          end do
          dl=(1.d0-temp)*d1+temp*d2
          dl=(1.d0-fice(l,ibv))*dl
          d(l,ibv)=dl
c**** calculate conductivity
          xklu=(1.d0-temp)*xku1+temp*xku2
          xklu=(1.d0-fice(l,ibv))*xklu
          xku(l,ibv)=xklu
          if(l.eq.1) then
            xk1=0.d0
            do i=1,imt-1
              xk1=xk1+qk(i,1)*xklm(0,i)
            end do
            xkl=xk1
            xkl=xkl/(1.d0+xkl/(-zc(1)*xkud))
            xkl=(1.d0-fice(1,ibv)*theta(1,ibv)/thets(1,ibv))*xkl
            xkl=max(zero,xkl)
            xk(1,ibv)=xkl
          else
            xk(l,ibv)=sqrt(xku(l-1,ibv)*xku(l,ibv))
          end if
        end do                  ! l
      end do                    ! ibv
ccc compute conductivity for topmodel (i.e. mean saturated conductivity)
      do ibv=1,2
        xkusa(ibv) = 0.d0
        dz_total = 0.d0
        do l=1,n
          xkusa(ibv) = xkusa(ibv) + xkus(l,ibv)*dz(l)
          dz_total = dz_total + dz(l)
        enddo
        xkusa(ibv) = xkusa(ibv) / dz_total
      enddo
c     add gravitational potential to hl
      do l=1,n
        do ibv=1,2
          h(l,ibv)=h(l,ibv)+zc(l)
        end do
      end do
      return
      end subroutine hydra

      subroutine hl0
c**** hl0 sets up a table of theta values as a function of matric
c**** potential, h.  h is tabulated in a geometric series from
c**** 0 to hmin, with a first step of delh1.  the theta values
c**** depend not only on the matric potential, but also on the
c**** soil texture.  we solve a cubic equation to determine
c**** theta as a function of h.  hl0 also outputs the conductivity
c**** and diffusivity tables.
c**** input:
c**** a - matric potential function parameters
c**** b - hydraulic conductivity function parameters
c**** p - hydraulic diffusivity function parameters
c**** sat - saturated thetas
c**** output:
c**** thm(j,i) - theta at j'th h point for texture i
c**** xklm(j,i) - conductivity at j'th h point for texture i
c**** dlm(j,i) - diffusivity at j'th h point for texture i
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      integer, parameter :: nexp=6
      real*8, parameter :: c=2.3025851d0
      real*8, dimension(4,imt-1), parameter :: a=reshape(
     &     (/                   ! matric potential coefficients for
     &     .2514d0,  0.0136d0, -2.8319d0,  0.5958d0, ! sand
     &     .1481d0,  1.8726d0,  0.1025d0, -3.6416d0, ! loam
     &     .2484d0,  2.4842d0,  0.4583d0, -3.9470d0, ! clay
     &     .8781d0, -5.1816d0, 13.2385d0,-11.9501d0/), ! peat
     &     (/4,imt-1/))
      real*8, dimension(4,imt-1), parameter :: b=reshape(
     &     (/                   ! conductivity coefficients for
     &     -0.4910d0, -9.8945d0,  9.7976d0, -3.2211d0, ! sand
     &     -0.3238d0,-12.9013d0,  3.4247d0,  4.4929d0, ! loam
     &     -0.5187d0,-13.4246d0,  2.8899d0,  5.0642d0, ! clay
     &     -3.0848d0,  9.5497d0,-26.2868d0, 16.6930d0/), ! peat
     &     (/4,imt-1/))
      real*8, dimension(4,imt-1), parameter :: p=reshape(
     &     (/                   ! diffusivity coefficients for
     &     -0.1800d0, -7.9999d0,  5.5685d0, -1.8868d0, ! sand
     &     -0.1000d0,-10.0085d0,  3.6752d0,  1.2304d0, ! loam
     &     -0.1951d0, -9.7055d0,  2.7418d0,  2.0054d0, ! clay
     &     -2.1220d0,  5.9983d0,-16.9824d0,  8.7615d0/), ! peat
     &     (/4,imt-1/))
      real*8 :: sat(imt-1) = (/.394d0,.537d0,.577d0,.885d0/)
      real*8 a1,a2,a3,alph0o,alpls1,arg(imt-1),delh1,delhn,dfunc
      real*8 diff,func,hmin,hs,s,sxtn,testh,xtol
      integer i,j,l,k,m,mmax

      sxtn=16.d0
      nth=2**nexp
      hlm(0)=0.0d0
      delh1=-0.00625d0
      hmin=-1000.d0
      delhn=delh1
c     solve for alph0 in s=((1+alph0)**n-1)/alph0
      s=hmin/delh1
      alph0=1.d0/8.d0
 10   alph0o=alph0
      alph0=(s*alph0+1.d0)**(1.d0/nth)-1.d0
      if(abs(alph0o-alph0).ge.1d-8) go to 10
      alpls1=1.0d0+alph0
      algdel=log(1.d0+alph0)
      do 100 j=1,nth
        hlm(j)=hlm(j-1)+delhn
        delhn=alpls1*delhn
 100  continue
      mmax=100
      xtol=1d-6
      do 200 i=1,imt-1
        thm(0,i)=1.00d0
        do 150 j=1,nth
          hs=-exp(c*(a(1,i)+a(2,i)+a(3,i)+a(4,i)))
          a1=a(3,i)/a(4,i)
          a2=(a(2,i)-(log(-hlm(j)-hs))/c)/a(4,i)
          a3=a(1,i)/a(4,i)
          testh=thm(j-1,i)
          do 130 m=1,mmax
            func=(testh**3)+(a1*(testh**2))+(a2*(testh))+a3
            dfunc=(3*testh**2)+(2*a1*testh)+a2
            diff=func/dfunc
            testh=testh-diff
            if(abs(diff).lt.xtol) go to 140
 130      continue
          print *,'max # iterations:',mmax
 140      thm(j,i)=testh
 150    continue
 200  continue
      do 280 j=0,nth
        do 245 i=1,imt-1
          xklm(j,i)=0.d0
          arg(i)=0.d0
          do 240 l=-1,2
            arg(i)=arg(i)+b(l+2,i)*thm(j,i)**l
 240      continue
          arg(i)=min(arg(i),sxtn)
          arg(i)=max(arg(i),-sxtn)
          xklm(j,i)=exp(c*arg(i))
 245    continue
        !print *, 'xklm ', j, (xklm(j,i), i=1,imt-1)
        do 265 i=1,imt-1
          dlm(j,i)=0.d0
          arg(i)=0.d0
          do 260 l=-1,2
            arg(i)=arg(i)+p(l+2,i)*thm(j,i)**l
 260      continue
          arg(i)=min(arg(i),sxtn)
          arg(i)=max(arg(i),-sxtn)
          dlm(j,i)=exp(c*arg(i))
 265    continue
 280  continue
      do 350 j=0,nth
        do 310 k=1,imt-1
          thm(j,k)=thm(j,k)*sat(k)
 310    continue
 350  continue
      return
      end subroutine hl0
      subroutine fl
c**** evaluates the flux between layers.
c**** input:
c**** h - soil potential of layers, m
c**** xk - conductivity of layers, m s-1
c**** zc - layer centers, m
c**** output:
c**** f - fluxes between layers, m s-1
c**** xinfc - infiltration capacity, m s-1
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c****
      integer ibv,l
      do ibv=1,2
        f(n+1,ibv)=0.d0
      end do
c****
      do ibv=1,2
        do l=2,n
          f(l,ibv)=-xk(l,ibv)*(h(l-1,ibv)-h(l,ibv))/(zc(l-1)-zc(l))
        end do
      end do
c**** put infiltration maximum into xinfc
      do ibv=1,2
        xinfc(ibv)=xk(1,ibv)*h(1,ibv)/zc(1)
      end do
      return
      end subroutine fl
      subroutine qsbal
c**** finds qs that balances fluxes.
c**** obtains qs by successive approximation.
c**** calculates evaporation.
c**** input:
c**** ch - heat conductivity coefficient from ground to surface
c**** vsm - surface layer wind speed, m s-1
c**** rho - air density, kg m-3
c**** eddy - transfer coefficient from surface to first atmosphere
c**** theta - water saturation of layers and canopy
c**** tp - temperatures of layers and canopy, c
c**** pres - atmospheric pressure at ground
c**** z1 - height of first layer, m
c**** zs - height of surface layer, m
c**** dz - layer thicknesses, m
c**** snowd - snow depths, equivalent water m
c**** pr - precipitation, m s-1
c**** q1 - mixing ratio of first layer
c**** fr - fraction of roots in layer
c**** fb - fraction of bare soil
c**** fv - fraction of vegetated soil
c**** hw - wilting point, m
c**** output:
c**** qs - mixing ratio at surface layer
c**** evap - evaporation from bare and vegetated regions, m s-1
c**** evapw - evaporation from wet canopy, m s-1, including from snow
c**** evapd - evaporation from dry canopy, m s-1
c**** evaps - evaporation from snow from canopy, m s-1
c**** betad - dry canopy beta, based on roots
c****
c**** uses: cond
ccc   include './soils101.com'
ccc   parameter (stbo=5.67032d-8)

c**** soils28   common block     9/25/90
      real*8,parameter :: eps = 2.d-5
      real*8 evap_max(2)
ccc   added declarations for local vars:
      real*8 qm1dt, xkf, tbs1, tcs1, qcv, qcs, epcs
      real*8 cna,dd,ed,qso,rho3,xl
      integer ibv,l,itr

ccc   first compute maximal amount of water available for evaporation
      do ibv=1,2
        evap_max(ibv) = 0.d0
        if ( ibv .eq. 1) evap_max(ibv) = evap_max(ibv) + pr
        do l=1,n
          evap_max(ibv) = evap_max(ibv) +
     &         (w(l,ibv)-dz(l)*thetm(l,ibv))/dt
        enddo
      enddo
c**** qm1 has mass of water vapor in first atmosphere layer, kg m-2
ccc   changing qm1 here is messy - may be fix later
      if(igcm.eq.-1) qm1=1d+7
      qm1dt=.001d0*qm1/dt
c     cna is the conductance of the atmosphere
      cna=ch*vsm
      rho3=.001d0*rho
      if(igcm.ge.0 .and. igcm.le.3) xl=eddy/(z1-zs)
c
c     modify tbs and tcs in the presence of snow to stabilize interface
c     xkf=(149.85d0/sqrt(dts))/(2.d0*sha*rho*cna
c     & +8.d0*stbo*(max(tbs,tcs)+tfrz)**3)
c     if(xkf.gt.1.d0)xkf=1.d0
ccc   debugging!!! - trying to fix evaporation from snow
      xkf=1.d0
      if(isn(1).ne.0.or.snowd(1).ne.0.d0)then
        tbs1=xkf*(tbs-ts+tfrz)+ts-tfrz
      else
        tbs1=tbs
      endif
      if(isn(2).ne.0.or.snowd(2).ne.0.d0)then
        tcs1=xkf*(tcs-ts+tfrz)+ts-tfrz
      else
        tcs1=tcs
      endif
      tbcs=fb*tbs1+fv*tcs1
c     calculate bare soil and canopy mixing ratios
      qb = qsat(tbs1+tfrz,lhe,pres)
c     qc = qsat(tcs1+tfrz,lhe,pres)
      qcv = qsat(tp(0,2)+tfrz,lhe,pres)
      qcs = qsat(tsn1(2)+tfrz,lhe,pres)

      qc = fm*qcs + (1.d0-fm)*qcv
c     on first iteration, assume beta's = 1
      betab=1.d0
      betav=1.d0
      if(igcm.ge.0 .and. igcm.le.3)
     &     qs=(fb*betab*cna*qb+fv*betav*cna*qc+xl*q1)
     &     /(fb*betab*cna+fv*betav*cna+xl+1d-12)
c     bare soils diffusivity dd
      dd=d(1,1)
c     betad is the the root beta for transpiration.
c     hw is the wilting point.
c     fr(l) is the fraction of roots in layer l
      betad=0.d0
      do 30 l=1,n
        betadl(l)=(1.d0-fice(l,2))*fr(l)*max((hw-h(l,2))/hw,zero)
        betad=betad+betadl(l)
 30   continue
      abetad=betad
c     canopy conductivity cnc
      call cond
c     surface layer mixing ratio to balance fluxes
      itr=1
 10   qso=qs
c     potential evaporation for bare soil and canopy
      epb=rho3*cna*(qb-qs)
c     epc=rho3*cna*(qc-qs)
      epc = rho3*cna*(qcv-qs)
      epcs = rho3*cna*(qcs-qs)
c     bare soil correction
c     diffusion limited flux ed
      ed=2.467d0*dd*(theta(1,1)-thetm(1,1))/dz(1)
c     correct if not potential evaporation
c     evap(1) is evaporation from bare soil
      if(snowd(1).gt.0.d0)then
        evap(1)=epb
      else
        evap(1)=min(epb,ed+pr)
      endif
      evap(1) = min( evap(1), evap_max(1) ) !limit to max amount of wate
c     vegetated soil correction
c     evap(2) is evaporation from vegetated land
c     evapd is dry evaporation (transpiration) from canopy
c     evapw is wet evaporation from canopy (from interception)
      if(epc.gt.0) then
        evapw=(1.d0-fm)*epc*fw
c**** limit the wet canopy evaporation to canopy water
        evapw=min(evapw,w(0,2)/dt)
        evapd=(1.d0-fm)*epc*fd
c**** if positive evaporation, limit dry canopy evaporation to trans
        betat=cnc/(cnc+cna+1d-12)
        evapd=min(evapd,evapd*betat)
        evapd = min( evapd, evap_max(2) ) ! limit to max amount of water
      else
        evapw=(1.d0-fm)*max(epc,-qm1dt)
        evapd=0.d0
      end if
c**** evaporation from vegetated snow region is from that part
c**** of the wet canopy that represents snow
c     evaps=epc*fm
      evaps = epcs*fm
ccc   limit it water legt after dry evap
      evaps = min( evaps, evap_max(2)-evapd )
c**** restrict condensation to water available in first atmosphere
      evap(1)=max(evap(1),-qm1dt)
      evap(2)=evapw+evapd+evaps
c**** calculate betas and q of surface layer
      if(epb.le.0.d0) then
        betab=1.0d0
      else
        betab=evap(1)/epb
      end if
      if(epc.le.0.d0) then
        betav=1.0d0
      else
        betav=(evap(2)-evaps)/epc
c     betav=evap(2)/((1.d0-fm)*epc+fm*epcs)
      end if
c**** for overall beta, use weighted average of betab and betav.
c**** don't use total evap over total potential evap.  this avoids
c**** the possibility of negative beta.
      beta=fb*betab+fv*betav
      abetav=betav
      abetat=betat
      abetab=betab
      abeta=beta
      acna=cna
      acnc=cnc
      if(igcm.ge.0 .and. igcm.le.3)
     &     qs=(fb*betab*cna*qb+fv*betav*cna*qcv + fv*fm*cna*qcs +xl*q1)
     &     /(fb*betab*cna+fv*betav*cna + fv*fm*cna  +xl+1.d-12)
c     &  qs=(fb*betab*cna*qb+fv*betav*cna*qc+xl*q1)
c     & /(fb*betab*cna+fv*betav*cna+xl+1.d-12)
 70   continue
c     loop back until qs converged
      if(itr.ge.60)then
        write(99,*)'qsbal:1',ijdebug,itr,qs,qso
        write(99,*)'qsbal:2',fb,betab,cna,qb
        write(99,*)'qsbal:3',fv,betav,qc,xl
        write(99,*)'qsbal:4',evap(1),evap(2),epb,epc
        write(99,*)'qsbal:5',evapw,evapd,ed,pr
        write(99,*)'qsbal:6',w(0,2),qm1dt,dt,cnc
        write(99,*)'qsbal:7',betad,q1,alai,rs
        write(99,*)'qsbal:8',srht,tp(1,1),tcs,ts
      endif
      if(itr.ge.64)then
        call outw(0)
        call abort
        stop 'qsbal'
      endif
      itr=itr+1
      if(abs(qso-qs).gt.eps)go to 10
      do 100 ibv=1,2
        l=2-ibv
c     snsh(ibv)=sha*rho*cna*(tp(l,ibv)-ts+tfrz)
        xlth(ibv)=evap(ibv)*elh
 100  continue
      snsh(1)=sha*rho*cna*(tbs1-ts+tfrz)
c     snsh(2)=sha*rho*cna*(tcs1-ts+tfrz)
      snsh(2)=sha*rho*cna*(tp(0,2)-ts+tfrz)
      snshs = sha*rho*cna*(tsn1(2)-ts+tfrz)
      snsh_dt = sha*rho*cna
      epb_dt = rho3*cna*qsat(tbs1+tfrz,lhe,pres)*dqsatdt(tbs1+tfrz,lhe)
      evaps_dt = rho3*cna*qsat(tsn1(2)+tfrz,lhe,pres)
     *     *dqsatdt(tsn1(2)+tfrz,lhe)
      return
      end subroutine qsbal

      subroutine flg
c**** calculates the ground water fluxes (to the surface)
c**** input:
c**** evap - evaporation from bare and vegetated regions, m s-1
c**** evapw - evaporation from wet canopy, m s-1
c**** pr - precipitation, m s-1
c**** htpr - heat of precipitation
c**** fsn - heat of fusion
c**** pre - extra precipitation, i.e. smowmelt, m s-1
c**** prfr - fraction by area of precipitation
c**** output:
c**** f - water fluxes from ground and canopy
c**** snowf - snow fall, equivalent water m s-1
c**** dr - canopy drip, m s-1
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 ptmp,ptmps,pfac,pm,pmax
c     calculate snow fall.  snowf is snow fall, m s-1 of water depth.
      snowf=0.d0
      if(htpr.lt.0.d0)snowf=min(-htpr/fsn,pr)
c     snowfs is the large scale snow fall.
      snowfs=0.d0
      if(htprs.lt.0.d0)snowfs=min(-htprs/fsn,prs)
c     bare soil
c     upward flux from first layer is evaporation less precipitation
      if(isn(1).ne.0.or.snowd(1).ne.0.d0)then
        pre(1)=flmlt(1)
        f(1,1)=-flmlt(1)
      else
        f(1,1)=-pr+evap(1)
        pre(1)=0.d0
      endif
c     upward flux from wet canopy, including evaporation from snow.
      if(isn(2).ne.0.or.snowd(2).ne.0.d0)then
        f(0,2)=evapw
      else
        f(0,2)=-pr+evapw
      endif
      ptmps=prs-snowfs
      ptmps=ptmps-evapw
      ptmp=pr-prs-(snowf-snowfs)
c     use effects of subgrid scale precipitation to calculate drip
      pm=1d-6
      pmax=fd*pm
      drs=max(ptmps-pmax,zero)
      dr=drs
      if(ptmp.gt.0.d0)then
        pfac=(pmax-ptmps)*prfr/ptmp
        if(pfac.ge.0.d0)then
          if(pfac.lt.30.000d0) dr=ptmp*exp(-pfac)
        else
          dr=ptmp+ptmps-pmax
        endif
      endif
c     vegetated soil
c     upward flux from soil surface is minus drip less snowfall
c     plus the evaporation from snow
      if(isn(2).ne.0.or.snowd(2).ne.0.d0)then
        f(1,2)=-flmlt(2)
        pre(2)=flmlt(2)
        f(0,2)=-flmlt(2)+f(0,2)
      else
        f(1,2)=-dr-snowf+evaps
        pre(2)=0.d0
      endif
      return
      end subroutine flg
      subroutine cond
c**** calculates the canopy conductance
c**** input:
c**** betad - beta due to roots
c**** alaie - effective leaf area index
c**** rs - minimum stomatal resistance, m s-1
c**** xinc - incoming solar radiation, w m-2
c**** tp - temperature of canopy, c
c**** tfrz - freezing point of water, 0 c in k
c**** output:
c**** cnc - canopy conductance, m s-1
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c**** adjust canopy conductance for soil water potential
      real*8 srht0
      cnc=betad*alaie/rs
c**** adjust canopy conductance for incoming solar radiation
      srht0=max(srht,zero)
      cnc=cnc*(srht0/c1)/(1.d0+srht0/c1)
      cnc=cnc/(1.d0+((tp(0,2)+tfrz-296.d0)/15.d0)**4)
      return
      end subroutine cond

      subroutine estimate_zbar
ccc compute average water table for topmodel
      integer l, ibv
      ! will insert this code later
      end subroutine estimate_zbar
      subroutine runoff
c**** calculates surface and underground runoffs.
c**** input:
c**** pre - effective precipitation, m s-1
c**** snowf - snow fall, equivalent water m s-1
c**** evap - evaporation, m s-1
c**** dr - canopy drip, m s-1
c**** xinfc - infiltration capacity, m s-1
c**** prfr - fraction of precipitation
c**** xk - conductivity, m s-1
c**** dz - layer thicknesses, m
c**** sl - slope
c**** sdstnc - interstream distance, m
c**** output:
c**** rnf - surface runoff
c**** rnff - underground runoff, m s-1
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c     use effects of subgrid scale rain
c     use precipitation that includes smow melt
      real*8 ptmp(2),ptmps(2)
      real*8 runfrac,prfac,rnfs
      real*8 f_k0_exp_l ! coefficient for the topmodel
      integer ibv,l
      if(isn(1).eq.0)then
        ptmps(1)=prs+pre(1)
        ptmp(1)=pr-prs
      else
        ptmps(1)=flmlt(1)
        ptmp(1)=0.d0
      endif
      if(isn(2).eq.0)then
        ptmps(2)=drs+pre(2)
        ptmp(2)=dr-drs
      else
        ptmps(2)=flmlt(2)
        ptmp(2)=0.d0
      endif
      do 10 ibv=1,2
        rnfs=max(ptmps(ibv)-xinfc(ibv),zero)
        rnf(ibv)=rnfs
        if(ptmp(ibv).gt.0.d0)then
          prfac=(xinfc(ibv)-ptmps(ibv))*prfr/ptmp(ibv)
          if(prfac.ge.0.d0)then
ccc   !! next line is different in new and old versions
ccc   !! i suppose some of them has a bug - check later
            if(prfac.lt.30.d0) rnf(ibv)=rnf(ibv)+ptmp(ibv)*exp(-prfac)
          else
            rnf(ibv)=ptmp(ibv)+ptmps(ibv)-xinfc(ibv)
          endif
ccc   !! following 3 lines didn't exist in old version
ccc   !! check what it is all about
c**** rosmp is runoff soil moisture parameter. set in ghinit.
          runfrac=(w(1,ibv)/ws(1,ibv))**so_%rosmp
          rnf(ibv)=(1.d0-runfrac)*rnf(ibv)
     $         +runfrac*(ptmp(ibv)+ptmps(ibv))
        endif
ccc   looks like this sometimes creates rnf<0
ccc   don't see anything wrong if i just set it to 0
        rnf(ibv) = max ( rnf(ibv), zero )
 10   continue
c     underground runoff
c     sl is the slope, sdstnc is the interstream distance
      do ibv=1,2
ccc this is some rough estimate for the expression f * k0 exp(zbar/f)
ccc in the topmodel expression for the runoff
!        f_k0_exp = 0.d0
!        do l=1,n
!          f_k0_exp = f_k0_exp + xkus(l,ibv)*w(l,ibv)/ws(l,ibv)*dz(l)
!        enddo
        do l=1,n
          rnff(l,ibv)=xku(l,ibv)*sl*dz(l)/sdstnc
/* #define do_topmodel_runoff */
#ifdef do_topmodel_runoff
          if ( ws(l,ibv) > 1.d-16 ) then
            f_k0_exp_l = (1.d0-fice(l,ibv))
     $           * xkus(l,ibv)*w(l,ibv)/ws(l,ibv)*dz(l)
          else
            f_k0_exp_l = 0.d0
          endif
c         print *,'rnff: ', l, rnff(l,ibv), f_k0_exp_l*exp( -top_index )
c         print *, xkus(l,ibv),w(l,ibv),ws(l,ibv),dz(l),top_index
          rnff(l,ibv)=f_k0_exp_l*exp( -top_index )
#endif
        end do
!        print *,'runoff: ', sl/sdstnc, exp( -top_index )
      end do
      return
      end subroutine runoff

      subroutine sink
c**** calculates water sinks from each soil layer
c**** input:
c**** rnf - surface runoff, m s-1
c**** rnff - underground runoff, m s-1
c**** evapd - evaporation from dry canopy, m s-1
c**** fr - fraction of roots in layers
c**** output:
c**** snk - water sink from layers, m s-1
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c**** underground runoff is a sink
      integer l, ibv
      do ibv=1,2
        do l=1,n
          snk(l,ibv)=rnff(l,ibv)
        end do
      end do
c**** remove transpired water from soil layers
      do l=1,n
        snk(l,2)=snk(l,2)+evapd*betadl(l)/(betad+1d-12)
      end do
c**** include effects of surface runoff in sink from first soil layers
      do ibv=1,2
        snk(1,ibv)=snk(1,ibv)+rnf(ibv)
      end do
      snk(0,2)=0.d0
      return
      end subroutine sink

      subroutine fllmt
c**** places limits on the soil water fluxes
c**** input:
c**** w - water in layers, m
c**** ws - saturated water in layers, m
c**** dts - current time step size, s
c**** f - water fluxes, m s-1
c**** snk - water sink from layers, m s-1
c**** rnf - surface runoff, m s-1
c**** snowd - snow depth, equivalent water m
c**** snowf - snow fall, equivalent water m s-1
c**** output:
c**** f - limited water fluxes, m s-1
c**** snk - limited water sinks, m s-1
c**** rnf - limited surface runoff, m s-1
c**** temp variables:
c**** snowdu - the upper bound on the snow depth at end of time step
c**** snowdl - the lower bound on the snow depth at end of time step
c**** trunc - fix for truncation on ibm mainframes
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 dflux,drnf,flmt,trunc
      integer l, ibv, ll
      trunc=1d-6
      trunc=1d-12
ccc   was 0 in older version - not sure if it is important
ccc   trunc = 0.d0
c     prevent over/undersaturation of layers 2-n
ccc         snowdu(1)=snowd(1)
ccc         snowdu(2)=snowd(2)
ccc         if(ht(1,1).lt.0)snowdu(1)=snowdu(1)+(snowf-evap(1))*dts
ccc         if(ht(1,2).lt.0)snowdu(2)=snowdu(2)+(snowf-evaps  )*dts
      do ibv=1,2
        ll=2-ibv
ccc         snowdu(ibv)=max(zero,snowdu(ibv))
        do l=n,2,-1
          flmt=(w(l,ibv)-ws(l,ibv)+trunc)/dts+f(l+1,ibv)-snk(l,ibv)
          f(l,ibv)=max(f(l,ibv),flmt)
          flmt=(w(l,ibv)-dz(l)*thetm(l,ibv)-trunc)/dts+
     $         f(l+1,ibv)-snk(l,ibv)
          f(l,ibv)=min(f(l,ibv),flmt)
        end do
      end do
c     prevent over/undersaturation of first layer
c     w(1) can include snow layer. - not any more !
c     bare soil
      flmt=(w(1,1)-ws(1,1)+trunc)/dts+f(2,1)-snk(1,1)
      drnf=max(zero,flmt-f(1,1))
      rnf(1)=rnf(1)+drnf
      snk(1,1)=snk(1,1)+drnf
      flmt=(w(1,1)-dz(1)*thetm(1,1)-trunc)/dts+f(2,1)-snk(1,1)
      drnf=min(zero,flmt-f(1,1))
      rnf(1)=rnf(1)+drnf
      snk(1,1)=snk(1,1)+drnf
c     prevent over/undersaturation of canopy layer
      if(isn(2).eq.0)then
        flmt=(ws(0,2)-w(0,2)-trunc)/dts+f(0,2)+snk(0,2)
        f(1,2)=min(flmt,f(1,2))
        flmt=(-w(0,2)+trunc)/dts+f(0,2)+snk(0,2)
        f(1,2)=max(flmt,f(1,2))
        dr=-f(1,2)
        dr=max(zero,dr)
      endif
c     prevent over/undersaturation of first layer
c     vegetated soil
      flmt=(w(1,2)-ws(1,2)+trunc)/dts+f(2,2)-snk(1,2)
      drnf=max(zero,flmt-f(1,2))
      rnf(2)=rnf(2)+drnf
      snk(1,2)=snk(1,2)+drnf
      flmt=(w(1,2)-dz(1)*thetm(1,2)-trunc)/dts+f(2,2)-snk(1,2)
      drnf=min(zero,flmt-f(1,2))
      rnf(2)=rnf(2)+drnf
      snk(1,2)=snk(1,2)+drnf
ccc   now trying to remove negative runoff
      do ibv=1,2
        l = 1
        do while ( rnf(ibv) .lt. 0.d0 .and. l .le. n )
ccc    this is how much water we can take from layer l
          dflux = f(l+1,ibv) + (w(l,ibv)-dz(l)*thetm(l,ibv))/dts
     &         - f(l,ibv) - snk(l,ibv)
          if( l .gt. 1) f(l,ibv) = f(l,ibv) - rnf(ibv)
          rnf(ibv) = rnf(ibv) + min( -rnf(ibv), dflux )
ccc    rnff always >= 0, use it first to compensate rnf<0
          if ( rnff(l,ibv) .lt. 0.d0 ) call abort ! just to be sure
          drnf = min( -rnf(ibv), rnff(l,ibv) )
          rnf(ibv) = rnf(ibv) + drnf
          rnff(l,ibv) = rnff(l,ibv) - drnf
          snk(l,ibv) = snk(l,ibv) - drnf
          l = l + 1
        enddo
ccc    check if rnf==0 up to machine accuracy
        if ( rnf(ibv) .lt. -1d-12 ) then
          print *, 'fllmt: rnf<0, ibv=',ibv,rnf(ibv)
          call abort   ! couldn't redistribute rnf<0 : evap is too big ?
        endif
ccc    if -1d-12 < rnf < 0. put it to 0 to avoid possible problems
ccc    actually for ground hydrology it is not necessary
        if (abs(rnf(ibv)).lt.1d-12) rnf(ibv) = 0.d0
      enddo
      return
      end subroutine fllmt

      subroutine fhlmt
c**** modifies soil heat fluxes to eliminate possible
c**** oscillation in presence of varying coefficient of drag.
c**** input:
c**** fh - heat fluxes
c**** shc - heat capacities
c**** w - water amounts
c**** fice - ice fraction
c**** dt - external time step
c**** output:
c**** fh - corrected heat fluxes
c**** parameter:
c**** dtpl - the max temperature change in a time step
c****
c**** add excess flux to flux of layer below
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 cc, dfh, dtpl, gm
      integer ibv,l,ll
      dtpl=1.0d0
      do ibv=1,2
       ll=2-ibv
       do l=ll,n-1
         cc=(shc(l,ibv)+w(l,ibv)*(fice(l,ibv)*shi
     $        +(1.d0-fice(l,ibv))*shw))
         gm=cc*dtpl/dt
         dfh=fh(l+1,ibv)-fh(l,ibv)
         if(abs(dfh).gt.gm)fh(l+1,ibv)=fh(l,ibv)+sign(gm,dfh)
       end do
      end do
      return
      end subroutine fhlmt

      subroutine xklh
c**** evaluates the heat conductivity between layers
c**** uses the method of devries.
c**** input:
c**** zb - soil layer boundaries, m
c**** zc - soil layer centers, m
c**** theta - soil water saturation
c**** fice - fraction of ice in layers
c**** alami - ice heat conductivity
c**** alamw - water heat conductivity
c**** tp - temperature of layers, c
c**** shw - specific heat of water
c**** shi - specific heat of ice
c**** shc - heat capacity of soil layers
c**** dz - layer thicknesses
c**** l,ibv - soil layer
c**** dts - the current time step
c**** output:
c**** xkh(l,ibv) - heat conductivities in each of the soil layers
c**** xkhm(l,ibv) - average heat conductivity between layer l and l-1
c****
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c     dimension xsha(ng,2),xsh(ng,2),gabc(3),hcwt(imt-1)
c
c     calculate with changing ga for air. ga is the depolarization
c     factor for air, calculated by linear interpolation from .333d0
c     at saturation to .035 at 0 water, following devries.
      real*8 gaa,gabc(3),gca,xa,xb,xden,xi,xnum,xs,xw
      real*8, save :: ba,hcwt(imt-1),hcwta,hcwtb,hcwti,hcwtw
      real*8, save :: xsha(ng,2),xsh(ng,2)
      integer i, j, ibv, l
      do ibv=1,2
        do l=1,n
          gaa=.298d0*theta(l,ibv)/(thets(l,ibv)+1d-6)+.035d0
          gca=1.d0-2.d0*gaa
          hcwta=(2.d0/(1.d0+ba*gaa)+1.d0/(1.d0+ba*gca))/3.d0
c     xw,xi,xa are the volume fractions.  don't count snow in soil lyr 1
          xw=w(l,ibv)*(1.d0-fice(l,ibv))/dz(l)
          xi=w(l,ibv)*fice(l,ibv)/dz(l)
          xa=(thets(l,ibv)-theta(l,ibv))
          xb=q(imt,l)
          xnum=xw*hcwtw*alamw+xi*hcwti*alami+xa*hcwta*alama+xsha(l,ibv)
     &         + xb*hcwtb*alambr
          xden=xw*hcwtw+xi*hcwti+xa*hcwta+xsh(l,ibv)+xb*hcwtb
          xkh(l,ibv)=xnum/xden
          if ( xkh(l,ibv) .lt. 0.d0 ) call abort()
        end do
      end do
c     get the average conductivity between layers
      do ibv=1,2
        do l=2,n
          xkhm(l,ibv)=((zb(l)-zc(l-1))*xkh(l,ibv)
     &         + (zc(l)-zb(l))*xkh(l-1,ibv)
     &         )/(zc(l)-zc(l-1))
        end do
      end do
c****
      return
      entry xklh0
c gabc's are the depolarization factors, or relative spheroidal axes.
      gabc(1)=.125d0
      gabc(2)=gabc(1)
      gabc(3)=1.d0-gabc(1)-gabc(2)
c hcwt's are the heat conductivity weighting factors
      hcwtw=1.d0
      hcwti=0.d0
      hcwtb=1.d0
      do i=1,imt-1
      hcwt(i)=0.d0
      end do
      do j=1,3
        hcwti=hcwti+1.d0/(1.d0+(alami/alamw-1.d0)*gabc(j))
        do i=1,imt-1
          hcwt(i)=hcwt(i)+1.d0/(1.d0+(alams(i)/alamw-1.d0)*gabc(j))
        end do
      end do
      hcwti=hcwti/3.d0
      do i=1,imt-1
        hcwt(i)=hcwt(i)/3.d0
      end do
      do ibv=1,2
        do l=1,n
          xsha(l,ibv)=0.d0
          xsh(l,ibv)=0.d0
          do i=1,imt-1
            xs=(1.d0-thm(0,i))*q(i,l)
            xsha(l,ibv)=xsha(l,ibv)+xs*hcwt(i)*alams(i)
            xsh(l,ibv)=xsh(l,ibv)+xs*hcwt(i)
          end do
        end do
      end do
      ba=alama/alamw-1.d0
      return
      end subroutine xklh

      subroutine flh
c**** evaluates the heat flux between layers
c**** subroutine fl must be called first
c**** input:
c**** zb - soil layer boundaries, m
c**** zc - soil layer centers, m
c**** theta - soil water saturation
c**** fice - fraction of ice in layers
c**** alami - ice heat conductivity
c**** alamw - water heat conductivity
c**** tp - temperature of layers, c
c**** shw - specific heat of water
c**** output:
c**** fh - heat flux between layers
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c****
      integer ibv, l
      do ibv=1,2
        fh(n+1,ibv)=0.d0
c total heat flux is heat carried by water flow plus heat conduction
        do l=2,n
          fh(l,ibv)=-xkhm(l,ibv)*(tp(l-1,ibv)-tp(l,ibv))/(zc(l-1)-zc(l))
          if(f(l,ibv).gt.0)then
            fh(l,ibv)=fh(l,ibv)+f(l,ibv)*tp(l,ibv)*shw
          else
            fh(l,ibv)=fh(l,ibv)+f(l,ibv)*tp(l-1,ibv)*shw
          endif
        end do
      end do
      return
      end subroutine flh

      subroutine flhg
c**** calculates the ground heat fluxes (to the surface)
c**** input:
c**** output:
c**** fh - heat fluxes from bare soil surface, and from canopy,
c****      and between canopy and vegetated soil.
c**** afhg - heat flux from ground to canopy
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c**** bare soil fluxes
      integer ibv, l
      ibv=1
      l=2-ibv
      if(isn(ibv).ne.0.or.snowd(ibv).ne.0.d0)then
        fh(l,ibv)=-fhsng(ibv)
        xlth(ibv)=hesn(ibv)
        thrm(ibv)=thrmsn(ibv)
      else
        fh(l,ibv)=xlth(ibv)+snsh(ibv)
        fh(l,ibv)=fh(l,ibv)-htpr
        fh(l,ibv)=fh(l,ibv)+thrm(ibv)-srht-trht
      endif
c
c**** canopy fluxes, and fluxes from masking snow
      ibv=2
      l=2-ibv
ccc thrm computed elsewhere
ccc      thrm(ibv)=stbo*(tcs+tfrz)**4
      fh(l,ibv)=(xlth(ibv)+snsh(ibv))
      fh(l,ibv)=fh(l,ibv)-htpr
c     fh(l,ibv)=fh(l,ibv)+(1.d0-fm)*(thrm(ibv)-srht-trht)
      fh(l,ibv)=fh(l,ibv)+(thrm(ibv)-srht-trht)
      if(isn(ibv).ne.0.or.snowd(ibv).ne.0.d0)then
        fh(1,2)=-fhsng(2)
        fh(0,2)=-fhsng(2)
     $       +(1.d0-fm)*(2.d0*thrm(ibv)-thrmsn(ibv)-srht-trht)
     &       +elh*(evapw+evapd)+(1.d0-fm)*snsh(ibv)
        xlth(ibv)=elh*(evapw+evapd)+hesn(ibv)
        thrm(ibv)=(1.d0-fm)*thrm(ibv)+fm*thrmsn(ibv)
      else
cccc using old formula as recommended by max
c     fh(l+1,ibv)=fm*fh(l,ibv)
ccc      fh(l+1,ibv)=fm*(thrm(ibv)-srht-trht)
        fh(l+1,ibv)=fm*fh(l,ibv)
c****
        fh(1,2)=fh(1,2)-shw*dr*tp(0,2)
        fh(1,2)=fh(1,2)-stbo*((tp(0,2)+tfrz)**4-(tp(1,2)+tfrz)**4)
      endif
      return
      end subroutine flhg

      subroutine sinkh
c**** calculates the heat removal from each layer
c**** input:
c**** shw - specific heat of water
c**** tp - temperature of layers, c
c**** snk - soil water sink in layers, m s-1
c**** output:
c**** snkh - heat sink from soil layers
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      integer ibv, l
      do ibv=1,2
        do l=1,n
          snkh(l,ibv)=shw*tp(l,ibv)*snk(l,ibv)
        end do
      end do
      return
      end subroutine sinkh

      subroutine retp
c**** evaluates the temperatures in the soil layers based on the
c**** heat values.  also executes snow melt.
c**** input:
c**** w - water in soil layers, m
c**** ht - heat in soil layers
c**** fsn - heat of fusion of water
c**** shc - specific heat capacity of soil
c**** shi - specific heat capacity of ice
c**** shw - specific heat capcity of water
c**** snowd - snow depth, equivalent water m
c**** fb - fraction of bare soil
c**** fv - fraction of vegetation
c**** fm - snow vegetation masking fraction (requires reth called first)
c**** output:
c**** tp - temperature of layers, c
c**** fice - fraction of ice of layers
c**** tbcs - temperature of bare soil, canopy, and snow as seen
c****        by atmosphere, c.  also called ground temperature.
c**** tcs - temperature of canopy and snow, c.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      integer ibv, l, ll
      do ibv=1,2
        ll=2-ibv
        do l=ll,n
          tp(l,ibv)=0.d0
          if(w(l,ibv).ge.1d-12)then
            fice(l,ibv)=-ht(l,ibv)/(fsn*w(l,ibv))
          else
            fice(l,ibv)=0.d0
          endif
          if( fsn*w(l,ibv)+ht(l,ibv) .lt. 0.d0 ) then
            tp(l,ibv)=(ht(l,ibv)+w(l,ibv)*fsn)/(shc(l,ibv)+w(l,ibv)*shi)
            fice(l,ibv)=1.d0
          else if(ht(l,ibv) .gt. 0.d0) then
            tp(l,ibv)=ht(l,ibv)/(shc(l,ibv)+w(l,ibv)*shw)
            fice(l,ibv)=0.d0
          endif
        end do
      end do

ccc this is a fix for undefined tsn1 at the beginning of soil routines
ccc probably should be moved to some other place
      do ibv=1,2
         tsn1(ibv) = 0.d0
         if (  wsn(1,ibv) .gt. 1.d-6 .and.
     &         hsn(1,ibv) + wsn(1,ibv)*fsn .lt. 0.d0  ) then
            tsn1(ibv) = (hsn(1,ibv) + wsn(1,ibv)*fsn)/(wsn(1,ibv)*shi)
         endif
ccc the following is a hack. it is necessary only at the beginning of th
ccc run, when some temperatures are not initialized properly.
ccc should be removed when program is rewritten in a more clean way...
         if ( wsn(1,ibv) .le. 1.d-6 ) then
            tsn1(ibv) = tp(2-ibv,ibv)
         endif
      enddo

      if(isn(2).eq.0)then
       tcs=tp(0,2)
      else
       tcs=(1.d0-fm)*tp(0,2)+fm*tsn1(2)
      endif
      if(isn(1).eq.0)then
       tbs=tp(1,1)
      else
       tbs=tsn1(1)
      endif
      tbcs=fb*tbs+fv*tcs
      etbcs=tbcs
      thrm(1)=stbo*(tp(1,1)+tfrz)**4
      thrm(2)=stbo*(tp(0,2)+tfrz)**4
c****
      if(tp(1,1).gt.100.d0.or.tp(0,2).gt.100.d0)then
      write(99,*)'retp tp bounds error'
      write(99,*)'ijdebug',ijdebug
      call reth
      call hydra
      call outw(1)
      call abort
      stop 'tp'
      endif
      return
      end subroutine retp

      subroutine advnc
c**** advances quantities by one time step.
c**** input:
c**** dt - time step, s
c**** dz - layer thickness, m
c**** tp - layer temperatures, c
c**** tfrz - freezing point of water, k
c**** w - soil water in layers, m
c**** snowd - snow depth, m
c**** f - water flux, m s-1
c**** snk - water sinks, m s-1
c**** ht - heat in soil layers
c**** fh - heat flux in soil layers
c**** snkh - heat sink in layers
c**** snowf - snow fall, m s-1 of equivalent water
c**** evap - evaporation, m s-1
c**** output:
c**** w - updater water in soil layers, m s-1
c**** ht - updated heat in soil layers
c**** snowd - updated snow depth, m s-1 of equivalent water
c**** rus - overall surface runoff, m s-1   replaced by aruns
c**** aruns - overall surface runoff, kg m-2
c**** aeruns - overall surface heat runoff, j m-2
c**** aerunu - underground heat runoff, j m-2
c**** uses:
c**** retp,reth,fl,flg,runoff,sink,sinkh,fllmt,flh,flhg.
c**** also uses surf with its required variables.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 dtm,tb0,tc0
      integer ibv,l,ll,limit,nit
      limit=200
      nit=0
      dtr=dt
      dr=0.d0
      call reth
      call retp
c     call hydra
      tb0=tp(1,1)
      tc0=tp(0,2)
ccc accm0 was not called here in older version - check
      call accm0
      do while ( dtr > 0.d0 )
        nit=nit+1
        if(nit.gt.limit)go to 900
        call hydra
        call wtab
        call qsbal
        call xklh
        call gdtm(dtm)
        dts=min(dtr,dtm)
        dtr=dtr-dts
        call snwlsi( dts )
        call fl
        call flg
        call runoff
        call sink
        call fllmt
        call sinkh
        call flh
        call flhg
c     call fhlmt
        do ibv=1,2
          ll=2-ibv
          do l=ll,n
            w(l,ibv)=w(l,ibv)+(f(l+1,ibv)-f(l,ibv)-snk(l,ibv))*dts
            ht(l,ibv)=ht(l,ibv)+(fh(l+1,ibv)-fh(l,ibv)-snkh(l,ibv))*dts
          end do
        end do
        w(0,2)=max(w(0,2),zero)
        call accm
        call reth
        call retp
c     call hydra
      enddo

      call accmf
c     call wtab
ccc   call outgh
      return
  900 continue
      write(99,*)'limit exceeded'
      write(99,*)'dtr,dtm,dts',dtr,dtm,dts
      write(99,*)'tb0,tc0',tb0,tc0
      call outw(2)
      call abort
      stop 'advnc'
      end subroutine advnc

      subroutine accm
c**** accumulates gcm diagnostics
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
c**** the following lines were originally called before retp,
c**** reth, and hydra.
      real*8 qsats,derun
      real*8 cpfac,dedifs,dqdt,el0,epen,h0
      integer ibv,l

      derun=shw*(fb*tp(1,1)*rnf(1)+fv*tp(1,2)*rnf(2))*dts
c**** need fix for negative energy because rain sometimes runs off
c**** frozen ground without freezing
      if (derun.lt.0) derun=0.
c****
      aeruns=aeruns+derun
      adifs=adifs-dts*(f(2,1)*fb+f(2,2)*fv)
      dedifs=f(2,1)*tp(2,1)
      if(f(2,1).lt.0.d0) dedifs=f(2,1)*tp(1,1)
      aedifs=aedifs-dts*shw*dedifs*fb
      dedifs=f(2,2)*tp(2,2)
      if(f(2,2).lt.0.d0) dedifs=f(2,2)*tp(1,2)
      aedifs=aedifs-dts*shw*dedifs*fv
c     alhg=alhg+dts*(evap(1)*(elh+tp(1,1)*shv)*fb+
c    *               evap(2)*(elh+tp(0,2)*shv)*fv)
      alhg=alhg+(xlth(1)*fb+xlth(2)*fv)*dts
c****
c**** the folowing lines were originally called after retp,
c**** reth, and hydra.
      af0dt=af0dt-dts*(fb*fh(1,1)+fv*fh(0,2)+htpr)
      af1dt=af1dt-dts*(fb*fh(2,1)+fv*fh(2,2))
      aruns=aruns+(fb*rnf(1)+fv*rnf(2))*dts
      do l=1,n
        arunu=arunu+(rnff(l,1)*fb+rnff(l,2)*fv)*dts
        aerunu=aerunu+ shw*( tp(l,1)*rnff(l,1)*fb
     *                     + tp(l,2)*rnff(l,2)*fv )*dts
c        aerunu=aerunu+(snkh(l,1)*fb+snkh(l,2)*fv)*dts ! wrong
ccc some new accumulators were added below this line
ccc check if their results are passed to corresponding programs
        arnff(l)=arnff(l)+(rnff(l,1)*fb+rnff(l,2)*fv)*dts
c**** add new diagnostics
        aw(l)=aw(l)+(w(l,1)*fb+w(l,2)*fv)*dts
        atp(l)=atp(l)+(tp(l,1)*fb+tp(l,2)*fv)*dts
        af(l)=af(l)+(f(l,1)*fb+f(l,2)*fv)*dts
      end do
      aw(0)=aw(0)+w(0,2)*dts
      atp(0)=atp(0)+tp(0,2)*dts
      af(0)=af(0)+f(0,2)*dts
      asnowd=asnowd+(snowd(1)*fb+snowd(2)*fv)*dts
      apre=apre+(pre(1)*fb+pre(2)*fv)*dts
      atbcs=atbcs+tbcs*dts
c**** end of new diagnostics
      if(isn(1).ne.0.or.snowd(1).ne.0.d0)then
        ashg=ashg+snshsn(1)*fb*dts
      else
        ashg=ashg+snsh(1)*fb*dts
      endif
      if(isn(2).ne.0.or.snowd(2).ne.0.d0)then
        ashg=ashg+( snsh(2)*(1.d0-fm) + snshsn(2) )*fv*dts
      else
        ashg=ashg+snsh(2)*fv*dts
      endif
      atrg=atrg+(thrm(1)*fb+thrm(2)*fv)*dts
c****
      aevapw=aevapw+(evapw*fv*dts)
      aevapd=aevapd+(evapd*fv*dts)
      aevaps=aevaps+(evaps*fv*dts)
      aevapb=aevapb+(evap(1)*fb*dts)
      aepc=aepc+(epc*fv*dts)
      aepb=aepb+(epb*fb*dts)
      afhg=afhg+(fh(1,2)*fv*dts)
c
      return
      entry accmf
c provides accumulation units fixups, and calculates
c penman evaporation.  should be called once after
c accumulations are collected.
      aruns=1000.0d0*aruns
      arunu=1000.0d0*arunu
      aevapw=1000.0d0*aevapw
      aevaps=1000.0d0*aevaps
      aevapd=1000.0d0*aevapd
      aevapb=1000.0d0*aevapb
      aepc=1000.0d0*aepc
      aepb=1000.0d0*aepb
      adifs=1000.d0*adifs
      af1dt=af1dt-aedifs
c**** fixup new diagnostics
      do l=1,n
       arnff(l)=1000.0d0*arnff(l)
      enddo
      do l=0,n
       aw(l)=1000.0d0*aw(l)/dt
       atp(l)=atp(l)/dt
       af(l)=1000.0d0*af(l)
      enddo
      apre=1000.0d0*apre
      asnowd=1000.0d0*asnowd/dt
      atbcs=atbcs/dt
c**** end of fuxup new diagnostics
c**** calculation of penman value of potential evaporation, aepp
      h0=fb*(snsh(1)+xlth(1))+fv*(snsh(2)+xlth(2))
c     h0=-atrg/dt+srht+trht
ccc   h0=-thrm(2)+srht+trht
      el0=elh*1d-3
ccc   cna=ch*vsm
      cpfac=sha*rho * ch*vsm
c**** replaced by standard function
c      t0=ts-tfrz
c      edelt=100.d0*pres*(qsat(ts,lhe,pres)-qs)/0.622d0
c      gamma=sha*100.d0*pres/(0.622d0*el0)
c      if(1.8d0*t0+48.0d0 .lt. 0.d0) then
c         delt=33.8639d0*(8.d0*0.00738d0*(0.00738d0*t0+0.8072d0)**7.d0
c     *       +0.000019d0*1.8d0)*100.0d0
c      else
c         delt=33.8639d0*(8.d0*0.00738d0*(0.00738d0*t0+0.8072d0)**7.d0
c     *       -0.000019d0*1.8d0)*100.0d0
c      end if
c      epen=(delt*h0+cpfac*edelt)/(el0*(delt+gamma))
      qsats=qsat(ts,lhe,pres)
      dqdt = dqsatdt(ts,lhe)*qsats
      epen=(dqdt*h0+cpfac*(qsats-qs))/(el0*dqdt+sha)
      aepp=epen*dt
      abetap=1.d0
      if (aepp.gt.0.d0) abetap=(aevapw+aevapd+aevapb)/aepp
      abetap=min(abetap,one)
      abetap=max(abetap,zero)
c     find final values of some derived variables
      esnowd=1000.d0*(fb*snowd(1)+fv*snowd(2))
      ezw=fb*zw(1)+fv*zw(2)
      ewtr1=1000.d0*( fb*w(1,1)*(1.d0-fice(1,1)) +
     +  fv*(w(1,2)*(1.d0-fice(1,2))+w(0,2)*(1.d0-fice(0,2))) )
      eice1=1000.d0*(fb*w(1,1)*fice(1,1) +
     +  fv*(w(1,2)*fice(1,2)+w(0,2)*fice(0,2)) )
      do l=0,n
        do ibv=1,2
          etp(l,ibv)=tp(l,ibv)
        end do
      end do
      return
      entry accm0
c zero out accumulations
c
      atrg=0.d0
      ashg=0.d0
      alhg=0.d0
      abetad=0.d0
      abetav=0.d0
      abetat=0.d0
      abetap=0.d0
      abetab=0.d0
      abeta=0.d0
      acna=0.d0
      acnc=0.d0
      aevapw=0.d0
      aevaps=0.d0
      aevapd=0.d0
      aevapb=0.d0
      aruns=0.d0
      arunu=0.d0
      aeruns=0.d0
      aerunu=0.d0
      adifs=0.d0
      aedifs=0.d0
      aepc=0.d0
      aepb=0.d0
      aepp=0.d0
      afhg=0.d0
      af0dt=0.d0
      af1dt=0.d0
c**** new diagnostics
      asnowd=0.d0
      atbcs=0.d0
      do l=1,n
       arnff(l)=0.d0
      enddo
      do l=0,n
       aw(l)=0.d0
       atp(l)=0.d0
       af(l)=0.d0
      enddo
      apre=0.d0
c**** end of new diagnostics
c
      return
      end subroutine accm

      subroutine gdtm(dtm)
c**** calculates the maximum time step allowed by stability
c**** considerations.
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 ak1,ak2(2),betas(2),cna,xk2(2),dldz2,dqdt,rho3,sgmm
      real*8, intent(out) :: dtm
      real*8 dtm1,dtm2,dtm3,dtm4,t450,xk1
      integer ibv,l
ccc         dimension qg(2),xk2(2),ak2(2),ak3(2)
ccc         dimension betas(2)
      t450=450.d0
c**** replaced with standard function
c      t0=ts-tfrz
c      if(1.8d0*t0+48.0d0 .lt. 0.d0) then
c         delt=33.8639d0*(8.d0*0.00738d0*(0.00738d0*t0+0.8072d0)**7.d0
c     *       +0.000019d0*1.8d0)*100.0d0
c      else
c         delt=33.8639d0*(8.d0*0.00738d0*(0.00738d0*t0+0.8072d0)**7.d0
c     *       -0.000019d0*1.8d0)*100.0d0
c      end if
c      dqdt=.622d0*delt/(100.d0*pres)
      dqdt=dqsatdt(ts,pres)*qsat(ts,lhe,pres)
      !qg(1)=qb
      !qg(2)=qc
c****
c**** first calculate timestep for water movement in soil.
      sgmm=1.0d0
      dldz2=0.d0
      do ibv=1,2
        do l=1,n
          dldz2=max(dldz2,d(l,ibv)/dz(l)**2)
        end do
      end do
      dtm=sgmm/(dldz2+1d-12)
      if(q(4,1).gt.0.d0)dtm=min(dtm,t450)
      dtm1=dtm
      if ( dtm .lt. 0.d0 ) call abort()
c****
c**** next calculate timestep for heat movement in soil.
      do ibv=1,2
        do l=1,n
          xk1=xkh(l,ibv)
          ak1=(shc(l,ibv)+((1.d0-fice(l,ibv))*shw+fice(l,ibv)*shi)
     &         *w(l,ibv))/dz(l)
          dtm=min(dtm,.5d0*ak1*dz(l)**2/(xk1+1d-12))
        end do
      end do
      dtm2=dtm
      if ( dtm .lt. 0.d0 ) call abort()
c****
c**** finally, calculate max time step for top layer bare soil
c**** and canopy interaction with surface layer.
c**** use timestep based on coefficient of drag
      cna=ch*vsm
      rho3=.001d0*rho
      if(epb.le.0.d0)then
       betas(1)=1.0d0
      else
       betas(1)=evap(1)/epb
      endif
      if(epc.le.0.d0)then
       betas(2)=1.0d0
      else
       betas(2)=(evap(2)-evaps)/epc
      endif
      do ibv=1,2
        l=2-ibv
        xk2(ibv)=sha*rho*cna
     &       + betas(ibv)*rho3*cna*elh*dqdt
     &       + 8.d0*stbo*(tp(l,ibv)+tfrz)**3
        ak2(ibv)=shc(l,ibv)+((1.d0-fice(l,ibv))*shw+fice(l,ibv)*shi)
     &       *w(l,ibv)
        dtm=min(dtm,ak2(ibv)/(xk2(ibv)+1d-12))
        if(ibv.eq.1)dtm3=dtm
        if(ibv.eq.2)dtm4=dtm
c
c prevent oscillation of top snow layer
c     if(isn(ibv).ne.0.or.snowd(ibv).ne.0.d0)then
c      ak3(ibv)=.05d0*shi*spgsn
c      dtm=min(dtm,ak3(ibv)/(xk2(ibv)+1.d-12))
c      if(ibv.eq.1)dtm5=dtm
c      if(ibv.eq.2)dtm6=dtm
c     endif
      end do
      if(dtm.lt.1.d0)then
       write(99,*) '*********** gdtm: ijdebug,fb,fv',ijdebug,fb,fv
       write(99,*)'dtm',dtm1,dtm2,dtm3,dtm4
       write(99,*)'xk2',xk2
       write(99,*)'ak2',ak2
       write(99,*)'snsh',snsh
       write(99,*)'xlth',xlth
       write(99,*)'dqdt',dqdt
       write(99,*)'ts,tfrz',ts,tfrz
       write(99,*)'dlt',tp(1,1)-ts+tfrz,tp(0,2)-ts+tfrz
      endif
c****
      return
      end subroutine gdtm


      subroutine outw(i)
c**** prints theta values at time step i
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      use filemanager, only: openunit
      real*8 day,scnds
      integer i,l,iday,ihour
      integer, save :: ichn = 0
      call wtab
      if( ichn == 0 ) call openunit("soil_outw", ichn)
      scnds=i*dt  !??? doesn't make any sense
      day=int(scnds/86400.d0)
      iday=day
      hour=int((scnds-86400.d0*day)/86400.d0)
      ihour=hour
c     print 1001
      write(ichn,1000)
      write(ichn,*)'general quantities (bare soil or vegetation)'
      write(ichn,*)'ij,dts',ijdebug,dts
cc    write(ichn,1021)
      write(ichn,1045)
      write(ichn,1023)'day= ',iday,'pr= ',pr,'ts= ',ts-tfrz,
     *     'q1= ',q1
      write(ichn,1023)'hour= ',ihour,'snowf= ',snowf,'tg= ',tg-tfrz,
     *     'qs= ',qs
      write(ichn,1043)'t1= ',t1-tfrz,'vg= ',vg,'ch= ',ch
      write(ichn,1044)'vsm= ',vsm
      write(ichn,1022)
      write(ichn,1021)
      write(ichn,1014)'bare soil   fb = ',fb
 1014 format(1x,a17,f4.2)
cc    write(ichn,1021)
      write(ichn,1025)
      write(ichn,1026)
      write(ichn,1027) snowd(1),rnf(1),evap(1),xinfc(1),zw(1),qb
      write(ichn,1021)
      write(ichn,1030)
      write(ichn,1031)
      do 100 l=1,n
      write(ichn,1040)l,theta(l,1),tp(l,1),fice(l,1),rnff(l,1),f(l,1),
     & h(l,1),xk(l,1),w(l,1),ws(l,1),shc(l,1),fh(l,1),ht(l,1),
     & q(1,l),q(2,l),q(3,l),q(4,l)
  100 continue
cc    write(ichn,1021)
      write(ichn,1022)
      write(ichn,1021)
      write(ichn,1014)'vegetation  fv = ',fv
cc    write(ichn,1021)
      write(ichn,1035)
      write(ichn,1036)
      write(ichn,1037) snowd(2),rnf(2),evap(2),xinfc(2),zw(2),qc,evapw,
     *     evapd,dr,fw
      write(ichn,1021)
      write(ichn,1030)
      write(ichn,1031)
      l=0
      write(ichn,1049)l,theta(l,2),tp(l,2),fice(l,2),rnff(l,2),f(l,2),
     & w(l,2),ws(l,2),shc(l,2),fh(l,2),ht(l,2)
      do 200 l=1,n
      write(ichn,1040)l,theta(l,2),tp(l,2),fice(l,2),rnff(l,2),f(l,2),
     & h(l,2),xk(l,2),w(l,2),ws(l,2),shc(l,2),fh(l,2),ht(l,2),
     & q(1,l),q(2,l),q(3,l),q(4,l)
  200 continue
cc    write(ichn,1021)
      write(ichn,1022)
      write(ichn,1021)
      write(ichn,1055)
 1055 format(1x,3x,9x,'   kgm-2',2x,9x,'    kgm-2',3x,9x,'   kgm-2',
     *     4x,9x,'1e6jm-2',3x,9x,'1e6jm-2')
      write(ichn,1060) aruns,aevapw,aepc,afhg,af0dt
 1060 format(1x,3x,'aruns = ',f9.4,2x,'aevapw = ',0pf9.4,4x,'aepc = ',
     *     0pf9.4,4x,'afhg = ',-6pf9.4,2x,'af0dt = ',-6pf9.4)
      write(ichn,1065) arunu,aevapd,aepb,atrg,htpr*dts
 1065 format(1x,3x,'arunu = ',f9.4,2x,'aevapd = ',0pf9.4,4x,'aepb = ',
     *     0pf9.4,4x,'atrg = ',-6pf9.4,2x,'aphdt = ',-6pf9.4)
      write(ichn,1070) aevapb,ashg,aeruns
 1070 format(1x,3x,'        ',9x,2x,'aevapb = ',0pf9.4,4x,'       ',
     *     9x,4x,'ashg = ',-6pf9.4,2x,'aerns = ',-6pf9.4)
      write(ichn,1073) alhg,af1dt,aedifs
 1073 format(1x,3x,'        ',9x,2x,'         ',9x,4x,'       ',
     *     9x,4x,'alhg = ',-6pf9.4,2x,'af1dt = ',-6pf9.4/
     *  1x,3x,'        ',9x,2x,'         ',9x,4x,'       ',
     *     9x,4x,'       ',2x,9x,'aedfs = ',-6pf9.4)
c**** more outw outputs
      write(ichn,*)'thrm ',thrm
      write(ichn,*)'xlth ',xlth
      write(ichn,*)'snsh ',snsh
      write(ichn,*)'htpr,srht,trht ',htpr,srht,trht
      return
1000  format(' ',121('='))
1001  format('1')
1010  format(1x,a20,f10.0)
1020  format(1x,a20,1pe12.4)
1021  format('0')
1022  format(' ',60('. '),'.')
1023  format(1x,a10,i10,a10,6pf10.2,2(a10,0pf8.2),a10,0pf8.4)
1043  format(1x,10x,10x,10x,10x,2(a10,0pf8.2),a10,0pf8.4)
1044  format(1x,10x,10x,38x,1(a10,f8.2))
1045  format(1x,20x,10x,'  1e-6ms-1',10x,4x,'t(c)',10x,4x,'ms-1')
1024  format(1x,6(a10,e10.2))
1015  format(1x,4(a8,f8.2))
1019  format(1x,12x,4(a8,f8.2))
1025  format(' ',5x,'snowd',7x,'rnf',6x,'evap',6x,'xinfc',
     *     8x,'zw',8x,'qb')
1026  format(' ','      mh2o',2x,'1e-6ms-1',2x,'1e-6ms-1',3x,
     *     '1e-6ms-1',3x,'      m',5x,'     ')
1027  format(' ',0pf10.4,6pf10.4,6pf10.4,1x,6pf10.1,0pf10.4,
     *     0pf10.4)
1030  format(' ',5x,'theta',3x,'tp',2x,'fice',4x,'runoff'
     & ,8x,'fl',9x,'h',8x,'xk',6x,'w',5x,'ws',8x,
     & 'shc',8x,'fh',8x,'ht',1x,'sand',1x,'loam',1x,'clay',1x,'peat')
1031  format(' ',5x,5x,2x,'(c)',2x,6x,'1e-6ms-1',2x,'1e-6ms-1',
     &     3x,'      m',2x,'1e-6ms-1',1x,'     m',1x,'     m',
     &     1x,'1e6jm-3c-1',4x,'  wm-2',3x,'1e6jm-2',4x,'%',4x,'%',4x,'%'
     &     ,4x,'%'/1x,125('-'))
1035  format(' ',5x,'snowd',7x,'rnf',6x,'evap',6x,'xinfc',
     *     8x,'zw',8x,'qc',5x,'evapw',5x,'evapd',8x,'dr',8x,'fw')
1036  format(' ','      mh2o',2x,'1e-6ms-1',2x,'1e-6ms-1',3x,
     *     '1e-6ms-1',2x,'       m',5x,'     ',2x,'1e-6ms-1',2x,
     *     '1e-6ms-1',2x,'1e-6ms-1','          ')
1037  format(' ',0pf10.4,6pf10.4,6pf10.4,1x,6pf10.1,0pf10.4,
     *     0pf10.4,6pf10.4,6pf10.4,6pf10.4,0pf10.2)
1040  format(1x,i3,f7.3,f5.1,f6.3,1p,6pf10.4,6pf10.4,0pf10.3,6pf10.4,
     *     0pf7.4,0pf7.4,1x,-6pf10.4,0pf10.4,-6pf10.4,4(2pf5.1))
1049  format(1x,i3,f7.3,f5.1,f6.3,1p,6pf10.4,6pf10.4,10x,10x,
     *     0pf7.4,0pf7.4,1x,-6pf10.4,0pf10.4,-6pf10.4,3(2pf5.1))
      end subroutine outw
c****
      subroutine wtab
c**** returns water table zw for ibv=1 and 2.d0
c**** input:
c**** zb - layer boundaries, m
c**** zc - soil centers, m
c**** dz - layer thicknesses, m
c**** h - soil potential of layers, m
c**** f - fluxes between layers, m s-1
c**** xk - conductivities of layers, m s-1
c**** output:
c**** zw(2) - water table for ibv=1 and 2, m
ccc   include 'soils45.com'
c**** soils28   common block     9/25/90
      real*8 denom,hmat,tol
      integer ibv,l
      tol=1d-6
      do 100 ibv=1,2
c**** find non-saturated layer
      do 10 l=n,1,-1
      if(w(l,ibv).lt.ws(l,ibv)*(1.d0-tol))go to 20
   10 continue
      l=1
   20 continue
c**** retrieve matric potential
c     write(6,*)'ij,n,l,hmat,ibv,xkl,ibv',ijdebug,n,l,hmat,ibv,xk(l,ibv)
      hmat=h(l,ibv)-zc(l)
c**** calculate denominator, and keep zw above zb(l+1)
      if(xk(l,ibv).le.1d-20) then
           denom=-2.d0*hmat/dz(l)
           go to 90
           end if
      denom=max(f(l,ibv)/xk(l,ibv)+1.d0,-2.d0*hmat/dz(l))
   90 continue
c**** calculate water table
c     write(6,*) 'denom',denom
      zw(ibv)=zb(l)-sqrt(-2.d0*hmat*dz(l)/(denom+1d-20))
  100 continue
      return
      end subroutine wtab

c***********************************************************************
      subroutine snwlsi( dts )
ccc     $     (epot,snsh,srht,trht,pr,htpr,xkth,cth,
ccc     & tg1,dzg1)
c***********************************************************************
c
c interface routine between soils.f and snow.f
c
c input & output: same as for snow.f, except for ibv variables
c to accomodate soils.f bare soil & vegetation
c
c input:
c dt - time step (s)
c elh - latent heat of evaporation (j m-3)
c fsn - latent heat of fusion (j m-3)
c epot(2) - potential evaporation (m s-1)
c snsh(2) - sensible heat (w m-2)
c pr(2) - precipitation (m s-1)
c htpr(2) - heat of precipitation
c xkth(2) - soil heat conductivity of first ground layer (w m-1 c-1)
c cth(2) - soil heat capacity of first layer (j m-2 c-1)
c tg1(2) - first layer ground temperature (c)
c dzg1 - firlst layer ground thickness (m)
c idgsn(2) - diagnostic print level. 0=no print
c
c output:
c isn(2) - 1 if snow, 0 if no snow
c
c output (common block soilsno):
c flmlt(2) - meltwater from snow to ground (m s-1)
c fhsng(2) - heat flux from snow to ground (w m-2)
c thrmsn(2) - thermal radiation from snow (w m-2)
c tsn(nlsn+1,2) - temperature of snow layers. nlsn+1=ground temp (c)
c rhosn(nlsn,2) - density of snow layers, (kg m-3)
c cvsn(nlsn,2) - specific heat by volume of snow layers (j m-3 k-1)
c xksn(nlsn+1,2) - heat conductivity between snow layers (w m-1 k-1)
c fisn(nlsn,2) - fraction of ice in snow layers (1)
c hesn(2) - heat of evaporating snow (w m-2)
c
c prognostic variables (common block soilsnp):
c dzsn(nlsn,2) - snow layer thicknesses, (m)
c wsn(nlsn,2) - water equivalent depth of snow layers, (m)
c hsn(nlsn,2) - heat in snow layers, (j m-2)
c nsn(2) - number of snow layers
c
ccc         implicit real*8 (a-h,o-z)

      use snow_model, only: snow_adv
      implicit none

      real*8 dts
ccc  snshsn(2) is global !
      real*8 epotsn(2),srhtsn(2),trhtsn(2)
ccc      real*8 xkthsn(2),cthsn(2)

ccc  local vars:
      real*8 tsn_surf, evap_sn_dt(2)
      integer ibv

      elh = 2.50d+9   ! we dont have this common block here

      epotsn(1)=epb
      epotsn(2)=evaps
      snshsn(1)=snsh(1)
      snshsn(2)=fm*snshs

      srhtsn(1)=srht
      srhtsn(2)=fm*srht
      trhtsn(1)=trht
      trhtsn(2)=fm*trht+(1.d0-fm)*thrm(2)

ccc!!! xkthsn cthsn are actually not used at the moment - should include
c!      xkthsn(1)=xkh(1,1)
c!      xkthsn(2)=xkh(1,2)
c!      cthsn(1)=shc(1,1)
c!      cthsn(2)=shc(1,2)

      evap_sn_dt(1) = epb_dt
      evap_sn_dt(2) = evaps_dt
      do ibv=1,2
c!!! should pass ground properties to snow_adv
        call snow_adv(dzsn(1,ibv), wsn(1,ibv), hsn(1,ibv), nsn(ibv),
     &       srhtsn(ibv), trhtsn(ibv), snshsn(ibv), htpr,
     $       epotsn(ibv),
     $       pr, dts,
     &       tp(1,ibv), dz(1), fr_snow(ibv),
     &       tsn_surf, flmlt(ibv), fhsng(ibv),
     &       thrmsn(ibv), snsh_dt, evap_sn_dt(ibv) )

        flmlt(ibv) = flmlt(ibv)/dts
        fhsng(ibv) = fhsng(ibv)/dts
ccc hack to remove flmlt < 0
        if ( flmlt(ibv) < 0.d0 ) then
          epotsn(ibv) = epotsn(ibv) + flmlt(ibv)
          flmlt(ibv) = 0.d0
        endif
c!! fix this later
        hesn(ibv) = epotsn(ibv) * elh
        isn(ibv) = 0
        if ( fr_snow(ibv) .gt. 1.d-12 ) isn(ibv) = 1
ccc!!!        snshsn(ibv) = snsh(ibv)  ! potential confusion
        tsn1(ibv)=tsn_surf
      enddo

      return
      end subroutine snwlsi

      subroutine set_snow
!@sum set_snow extracts snow from the first soil layer and initializes
!@+   snow model prognostic variables
!@+   should be called when model restarts from the old restart file
!@+   ( which doesn't contain new snow model (i.e. 3 layer) data )
c
c input:
c snowd(2) - landsurface snow depth
c w(l,2)   - landsurface water in soil layers
c ht(l,2)  - landsurface heat in soil layers
c fsn      - heat of fusion
c shi      - specific heat of ice
c shc(l,2) - heat capacity of soil layers
c
c output:
c dzsn(lsn,2) - snow layer thicknesses
c wsn(lsn,2)  - snow layer water equivalent depths
c hsn(lsn,2)  - snow layer heat contents
c tsn1(2)     - snow top temperature
c isn(2)      - 0 if no snow, 1 if snow
c nsn(2)      - number of snow layers
c snowd(2)
c w(l,2)
c ht(l,2)
c
c calling sequence:
c
c     assignment of w,ht,snowd
c     call ghinij(i,j,wfc1)
c     call snwin
c note: only to be called when initializing from landsurface
c       prognostic variables without the snow model.
c
ccc         include './soils101.com'
      integer ibv

c outer loop over ibv
      do ibv=1,2

c initalize all cases to nsn=1
       nsn(ibv)=1

ccc since we don't know what kind of data we are dealing with,
ccc better check it

       if( snowd(ibv) .gt. w(1,ibv)-dz(1)*thetm(1,ibv)  ) then
          write(99,*) 'snowd corrected: old=', snowd(ibv)
          snowd(ibv) = w(1,ibv)-dz(1)*thetm(1,ibv) - 1.d-10
          write(99,*) '                 new=', snowd(ibv)
          if ( snowd(ibv) .lt. -0.001d0 ) call abort
          if ( snowd(ibv) .lt. 0.d0 ) snowd(ibv) = 0.d0 ! rounding error
       endif

c if there is no snow, set isn=0.  set snow variables to 0.d0
       if(snowd(ibv).le.0.d0)then
        isn(ibv)=0
        dzsn(1,ibv)=0.d0
        wsn(1,ibv)=0.d0
        hsn(1,ibv)=0.d0
        tsn1(ibv)=0.d0
        fr_snow(ibv) = 0.d0
       else

c given snow, set isn=1.d0
        isn(ibv)=1
c!!!        dzsn(1,ibv)=snowd(ibv)/spgsn
c!!!    replacing prev line considering rho_snow = 200
        dzsn(1,ibv)=snowd(ibv) * 5.d0
        wsn(1,ibv)=snowd(ibv)
c!!! actually have to compute fr_snow and modify dzsn ...
        fr_snow(ibv) = 1.d0

c given snow, temperature of first layer can't be positive.
c the top snow temperature is the temperatre of the first layer.
        if(fsn*w(1,ibv)+ht(1,ibv).lt.0.d0)then
         tsn1(ibv)=(ht(1,ibv)+w(1,ibv)*fsn)/(shc(1,ibv)+w(1,ibv)*shi)
        else
         tsn1(ibv)=0.d0
        endif

c use snow temperature to get the heat of the snow
        hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn

c subtract the snow from the landsurface prognositic variables
        w(1,ibv)=w(1,ibv)-wsn(1,ibv)
        ht(1,ibv)=ht(1,ibv)-hsn(1,ibv)

ccc and now limit all the snow to 5cm water equivalent
        if ( snowd(ibv) .gt. 0.05d0 ) then
          snowd(ibv) = 0.05d0
          dzsn(1,ibv)= snowd(ibv) * 5.d0
          wsn(1,ibv)= snowd(ibv)
          hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
        endif

       endif
      enddo

      return
      end subroutine set_snow


      end module sle001
