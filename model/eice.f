      subroutine eice(m,n,mm,nn,k1m,k1n)
c
c --- estimate the ice formed (odmsi) in order to keep ocean above freezing
c --- and corresponding saltflux (salflx2)
c
      USE SEAICE, only : fsss,tfrez
      USE CONSTANT, only : lhm,shi,shw
c
      implicit none
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      real tmelt,thin,rhoice,kice,fusion,saldif,rate,tmxl,dpth,thkmax,
     .     paybak,borrow,total,totalx,qmax,qmx(jdm),heatfx,dtdflx,
     .     icex(idm,jdm),work(idm,jdm),dm,salflx2(idm,jdm),top,thkinv
      integer imx(jdm),jmx(jdm),k1
      logical dosmoo
c
c --- tmelt  = melting point (deg)
c --- thin   = min. ice thickness (m)
c --- thkmax = max. ice thickness (m)
c --- rhoice = ice density (kg/m^3)
c --- kice   = heat conductivity in ice (W/m/deg)
c --- fusion = latent heat of fusion (J/kg)
c --- saldif = salinity difference water minus ice
c --- rate   = max. ice freezing and melting rate (m/sec)
c --- dtdflx = d (srf.temperature) / d (heat flux)  (deg m^2/W)
c --- sfxice = net total heat flux between atm and ice (W/m^2)
c --- salflx = salt flux (implied by fresh water flux)
c --- heatfx = heat flux through ice sheet (W/m^2)
c --- covice = ice coverage (rel.units)
c --- thkice = grid-box averaged ice thickness (m)
c --- temice = ice surface temperature
c
      data thin/0.2/,rhoice/917./,thkmax/10./
     .    ,kice/2.04/,fusion/334.e3/,rate/5.e-6/,dtdflx/0.05/
c
c --- energy loan: add extra energy to the ocean to keep SST from dropping
c --- below tmelt in winter. return this borrowed energy to the 'energy bank'
c --- in summer as quickly as surflx > 0 allows.
c --- salt loan: analogous to energy loan.
c     
      total=0.
      totalx=0.
      dosmoo=.false.
c$OMP PARALLEL DO PRIVATE(tmxl,tmelt,borrow,paybak,saldif,dpth
c$OMP+  ,top,thkinv,kn)
c$OMP+ REDUCTION(+:total) SHARED(dosmoo,xpivn)
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
cdiag total=total+thkice(i,j)*scp2(i,j)
c
      saldif=max(0.,saln(i,j,k1n)-10.)
      dpth=max(dp(i,j,k1n),thkmin*onem)
c
c --- calculate hypothetical mixed-layer temp due to diab. forcing 
      if (dp(i,j,k1n).le.0.) then
        write (lp,'(i9,2i5,a)') nstep,i,j,'  zero mxlayr thickness'
        stop '(eice error)'
      end if
c
      odmsi(i,j)=0.
      salflx2(i,j)=0.
      tmxl=temp(i,j,k1n)+surflx(i,j)*delt1*g/(spcifh*dpth) 
      tmelt=tfrez(saln(i,j,k1n),0.)
c
      if (tmxl.lt.tmelt) then
        borrow=min((tmelt-tmxl)*spcifh*dpth/(delt1*g),
     .           rate*fusion*rhoice)			! > 0 W/m^2
c
c --- add energy to bring tmxl back to tmelt (only if tmxl < tmelt)
c
        surflx(i,j)=surflx(i,j)+borrow
        odmsi(i,j)=-borrow/                      ! odmsi: kg/m2/sec
     .          (tmelt*shi-lhm*(1.-0.001*fsss*saln(i,j,k1n))) ! > 0
        salflx2(i,j)=                   ! salflx: g/m2/sec fresh WT + salt flux
     .               odmsi(i,j)*saln(i,j,k1n)*(1.-fsss)
      endif
c
c --- build up time integrals of surface fluxes
      eminpav(i,j)=eminpav(i,j)-(salflx(i,j)+salflx2(i,j))
     .                               *thref/saln(i,j,k1n)
      surflav(i,j)=surflav(i,j)+surflx(i,j)
      tauxav(i,j)=tauxav(i,j)+taux(i,j)
      tauyav(i,j)=tauyav(i,j)+tauy(i,j)
c
c --- deposit brine below 500 m or in bottom fourth of column in SH
c
      if (salflx2(i,j).gt.0..and.i.gt.xpivn) then
        top=min(500.*onem,p(i,j,kk+1)*.75)
        thkinv=1./(p(i,j,kk+1)-top)
        do 11 k=1,kk
        kn=k+nn
 11     saln(i,j,kn)=saln(i,j,kn)+salflx2(i,j)*delt1*g*thkinv
     .   *max(0.,min(1.,(p(i,j,k+1)-top+epsil)/(dp(i,j,kn)+epsil)))
      else
        salflx(i,j)=salflx(i,j)+salflx2(i,j)
      end if
 10   continue
c$OMP END PARALLEL DO
c
c     i=itest
c     j=jtest
c     write(*,'(a,/,7e11.2)')
c    .   'chk oflx  emnp  salt  ice  taux  tauy  u*',
c    .  oflxa2o(i,j),oemnp(i,j),osalt(i,j),oice(i,j),
c    .   taux(i,j),tauy(i,j),ustar(i,j)
c
      return
      end
c
c
c> Revision history
c>
c> June 2000 - conversion to SI units
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> Aug. 2000 - revised ice surf.temp. (temice) calculation
c> Mar. 2001 - corrected error in -heatfx- calculation
c> Aug. 2001 - introduced mininum thickness -dpth- in loan calculation
c> Sep. 2001 - corrected -surflx- sign error in -tmxl- calculation
c> Feb. 2002 - re-introduced smoother to spread ice thicker than -thkmax-
c> Oct 2004 - keep tract on ice mass only, add brine subsurface in SH
