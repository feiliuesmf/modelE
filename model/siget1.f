      real function sigocn(t,s)
c
c --- sigma, dsigma/dt, and dsigma/ds as functions of temp and salinity
c
      implicit none
      real t,s
c
#include "state_eqn.h"
c
ccc   sig=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))*1.e-3		! cgs
c7    sigocn=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))		!  SI
      sigocn=(c1+s*(c3+s*(c8+c9*t))+t*(c2+c5*s+t*(c4+c7*s+c6*t)))	!  SI
      return
      end
c
c
      real function dsigdt(t,s)
      implicit none
      real t,s
c
#include "state_eqn.h"
c
ccc   dsigdt=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))*1.e-3		! cgs
      dsigdt=(c2+s*(c5+c9*s)+2.*t*(c4+c7*s+1.5*c6*t))		!  SI
      return
      end
c
c
      real function dsigds(t,s)
      implicit none
      real t,s
c
#include "state_eqn.h"
c
ccc   dsigds=(c3+t*(c5+t*c7))*1.e-3				! cgs
      dsigds=c3+2*s*(c8+c9*t)+t*(c5+t*c7)					!  SI
      return
      end
c
c
      real function tofsig(sigm,salin)
c
c --- temp (deg c) as a function of sigma and salinity
c
      implicit none
      real sigm,salin,sqq,a0,a1,a2,cubr,cubq,cuban,cubrl,cubim,athird
     .    ,sigocn
      parameter (athird=1./3.)
c
#include "state_eqn.h"
c
      a0=(c1+salin*(c3+c8*salin))/c6
      a1=(c2+salin*(c5+c9*salin))/c6
      a2=(c4+c7*salin)/c6
c
      cubq=athird*a1-(athird*a2)**2
ccc   cubr=athird*(.5*a1*a2-1.5*(a0-1.e3*sigm/c6))		! cgs
      cubr=athird*(.5*a1*a2-1.5*(a0-     sigm/c6))		!  SI
     .   -(athird*a2)**3
c
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
c
      cuban=athird*atan2(sqrt(max(0.,-(cubq**3+cubr**2))),cubr)
      sqq=sqrt(-cubq)
      cubrl=sqq*cos(cuban)
      cubim=sqq*sin(cuban)
      tofsig=-cubrl+sqrt(3.)*cubim-athird*a2
      if (abs(sigocn(tofsig,salin)-sigm).gt.1.e-6) write(*,100)
     .   tofsig,salin,sigm,sigocn(tofsig,salin)
 100  format ('tofsig,sal,old/new sig =',4f9.3)
      return
      end
c
c
      real function sofsig(sigma,tem)
      implicit none
      real sigma,tem,bb,aa,cc
c
#include "state_eqn.h"
c
cc    sofsig=(sigma*1.e3				     ! cgs
c7    sofsig=(sigma					     ! SI
c7   .  -c1-tem*(c2+tem*(c4+c6*tem)))/(c3+tem*(c5+c7*tem))
      aa=c8+c9*tem
      bb=c3+c5*tem+c7*tem*tem
      cc=c1+c2*tem+c4*tem*tem+c6*tem*tem*tem-sigma
      sofsig=(-bb+sqrt(bb*bb-4.*aa*cc))/(2.*aa)
      return
      end
c
      real function qsatur(t)
c
c --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)
c
      qsatur=.622e-3*(6.107799961e+00+t*(4.436518521e-01
     .            +t*(1.428945805e-02+t*(2.650648471e-04
     .            +t*(3.031240396e-06+t*(2.034080948e-08
     .            +t* 6.136820929e-11))))))
      return
      end
c
c
      real function kappaf(t,s,prs)
c
c --- compressibility coefficient kappa^(theta) from sun et al. (1999)
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
      real t,s,pres,prs,tdif,sdif
c
#include "state_eqn.h"
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     kappaf=0.				! no thermobaricity
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      sdif=max(30.,min(38.,s))-refsal
      tdif=max(-2.,min(32.,t))-reftem
css   prs=max(pres,200.*onem)			!  experimental
      kappaf=exp(sclkap * (tdif*(qt+tdif*(qtt+tdif*qttt)
     .       +.5*(prs+pref)*(qpt+sdif*qpst+tdif*qptt))
     .       +sdif*(qs+tdif*qst))*(prs-pref)) - 1.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- claes rooth's original formula:
ccc   kappaf=-.075e-12*t*(3.-.06*t+.0004*t*t)*(prs-pref)/thref
c
      return
      end
c
c
      real function tofvir(virt,salin,prs)
c
c --- temp (deg c) as a function of virt.pot.density, salinity, and pressure
c
      implicit none
      real virt,salin,prs,sqq,a0,a1,a2,cubr,cubq,cuban,cubrl,cubim,
     .     b0,b1,b2,fac,sdif,check,athird,sig,kappaf
      external sig,kappaf
      parameter (athird=1./3.)
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
#include "state_eqn.h"
c
ccc      a0=(c1+c3*salin)/c6
ccc      a1=(c2+c5*salin)/c6
ccc      a2=(c4+c7*salin)/c6
c
      fac=sclkap*(prs-pref)/thref
      sdif=max(30.,min(38.,salin))-refsal	! bounds must agree with kappaf
      b2=qtt+.5*(prs+pref)*qptt
      b1=qt+.5*(prs+pref)*(qpt+sdif*qpst)+sdif*qst
      b0=sdif*qs
c
      a0=(c1+c3*salin+fac*
     .  (b0-reftem*(reftem*(reftem*qttt-b2)+b1)))/(c6+fac*qttt)
      a1=(c2+c5*salin+fac*
     .  (b1+reftem*(3.*reftem*qttt-2.*b2)))      /(c6+fac*qttt)
      a2=(c4+c7*salin+fac*
     .  (b2-3.*reftem*qttt))                     /(c6+fac*qttt)
c
      cubq=athird*a1-(athird*a2)**2
      cubr=athird*(.5*a1*a2-1.5*(a0-virt/(c6+fac*qttt)))
     .   -(athird*a2)**3
c
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
c
      cuban=athird*atan2(sqrt(max(0.,-(cubq**3+cubr**2))),cubr)
      sqq=sqrt(-cubq)
      cubrl=sqq*cos(cuban)
      cubim=sqq*sin(cuban)
      tofvir=-cubrl+sqrt(3.)*cubim-athird*a2
      tofvir=max(tofvir,-1.8)
c
ccc      check=sig(tofvir,salin)+kappaf(tofvir,salin,prs)
ccc      if (abs(check-virt).gt.1.e-5)
ccc     .write (lp,100) tofvir,salin,prs/onem,virt,check
 100  format ('tofvir,sal,prs,old/new sig =',2f9.3,f9.2,2f9.3)
c
      return
      end
c
c
c> Revision history:
c>
c> Jan. 2001 - modified kappaf to allow nonzero temp.ref.value
c> Nov. 2002 - corrected error in kappaf
c> June 2003 - added function -tofvir-
