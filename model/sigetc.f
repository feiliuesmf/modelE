      real function sigocn(t,s)
c
c --- sigma, dsigma/dt, and dsigma/ds as functions of temp and salinity
c
      implicit none
      real t,s
c
      include 'state_eqn.h'
c
ccc   sig=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))*1.e-3                ! cgs
c7    sigocn=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))                !  SI
      sigocn=(c1+s*(c3+s*(c8+c9*t))+t*(c2+c5*s+t*(c4+c7*s+c6*t)))       !  SI
      return
      end
c
c
      real function dsigdt(t,s)
      implicit none
      real t,s
c
      include 'state_eqn.h'
c
      dsigdt=(c2+s*(c5+c9*s)+2.*t*(c4+c7*s+1.5*c6*t))           !  SI c9
      return
      end
c
c
      real function dsigds(t,s)
      implicit none
      real t,s
c
      include 'state_eqn.h'
c
c     dsigds=(c3+t*(c5+t*c7))                                        ! SI c7
      dsigds=c3+2*s*(c8+c9*t)+t*(c5+t*c7)                       ! SI c9
      return
      end
c
c
      real function tofsig(sigm,salin)
c
c --- temp (deg c) as a function of sigma and salinity (this routine mimics
c --- the statement function of the same name in state_eqn.h.)
c
      implicit none
      real sigm,salin,s,sqq,a0,a1,a2,cubr,cubq,cuban,cubrl,cubim,athird
      parameter (athird=1./3.)
c
      include 'state_eqn.h'
c
      a0=(c1+salin*(c3+c8*salin))/c6
      a1=(c2+salin*(c5+c9*salin))/c6
      a2=(c4+c7*salin)/c6
c
      cubq=athird*a1-(athird*a2)**2
ccc   cubr=athird*(.5*a1*a2-1.5*(a0-1.e3*sigm/c6))                ! cgs
      cubr=athird*(.5*a1*a2-1.5*(a0-     sigm/c6))                !  SI
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
ccc      if (abs(sig(tofsig,salin)-sigm).gt.1.e-6) write (lp,100)
ccc     .   tofsig,salin,sigm,sig(tofsig,salin)
 100  format ('tofsig,sal,old/new sig =',2f9.3,3p,2f9.3)
      return
      end
c
c
      real function sofsig(sigma,tem)
      implicit none
      real sigma,tem,bb,aa,cc
c
      include 'state_eqn.h'
c
c7    sofsig=(sigma                                          ! SI
c7   .  -c1-tem*(c2+tem*(c4+c6*tem)))/(c3+tem*(c5+c7*tem))
      aa=c8+c9*tem
      bb=c3+c5*tem+c7*tem*tem
      cc=c1+c2*tem+c4*tem*tem+c6*tem*tem*tem-sigma
      sofsig=(-bb+sqrt(bb*bb-4.*aa*cc))/(2.*aa)
      return
      end
c
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
      include 'dimensions.h'
      include 'common_blocks.h'
      real t,s,prs,tdif,sdif
c
      include 'state_eqn.h'
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      kappaf=0.                                ! no thermobaricity
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      sdif=max(30.,min(38.,s))-refsal
ccc      tdif=max(-2.,min(32.,t))-reftem
ccc      kappaf=sclkap * (tdif*(qt+tdif*(qtt+tdif*qttt)
ccc     .       +.5*(prs+pref)*(qpt+sdif*qpst+tdif*qptt))
ccc     .       +sdif*(qs+tdif*qst))*(prs-pref)/thref
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- claes rooth's original formula:
ccc   kappaf=-.075e-12*t*(3.-.06*t+.0004*t*t)*(prs-pref)/thref
c
      return
      end
c
c
      real function sigloc(t,s,prs)
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
c
      include 'state_eqn.h'
c
      real t,s,prs,c1p,c2p,c3p,c4p,c5p,c6p,c7p
      c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
      c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
      c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
      c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
      c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
c
      sigloc=c1p(prs)+c3p(prs)*s+
     &    t*(c2p(prs)+c5p(prs)*s+t*(c4p(prs)+c7p(prs)*s+c6p(prs)*t))
      return
      end
c
c
      real function dsiglocdt(t,s,prs)
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
c
      include 'state_eqn.h'
c
      real t,s,prs,c2p,c4p,c5p,c6p,c7p
ccc   c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
      c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
ccc   c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
      c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
      c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
c
      dsiglocdt=c2p(prs)+c5p(prs)*s+
     &    2.*t*(c4p(prs)+c7p(prs)*s+1.5*c6p(prs)*t)
      return
      end
c
c
      real function dsiglocds(t,s,prs)
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
c
      include 'state_eqn.h'
c
      real t,s,prs,c3p,c5p,c7p
ccc   c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
ccc   c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
      c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
ccc   c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
ccc   c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
c
      dsiglocds=c3p(prs)+t*(c5p(prs)+t*c7p(prs))
      return
      end

c
c
c> Revision history:
c>
c> Jan. 2001 - modified kappaf to allow nonzero temp.ref.value
c> Nov. 2002 - corrected error in kappaf
c> July 2005 - added sig,dsigdt,dsigds functions for in-situ density ('sigloc')
