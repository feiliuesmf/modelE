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
      real function kappaf(t,s,prs,kkf)
c
c --- compressibility coefficient kappa^(theta) from sun et al. (1999)
c
      implicit none
      include 'dimensions.h'
      include 'common_blocks.h'
      real t,s,prs,kappaf1
      integer kkf
      external kappaf1
c
      include 'state_eqn.h'
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- thermobaric compressibility coefficient (integral from prs to pref)
c --- Sun et.al. (1999) JPO 29 pp 2719-2729
c --- kappaf1 used internally to simplify offsetting T and S.
c --- always invoke via kappaf.
c --- t: potential temperature; s: psu; prs: pressure; kkf: ref.state
c ---     example: kappaf(4.5,34.5,1.e7,1) =  0.11411243
c ---     example: kappaf(4.5,34.5,1.e7,2) =  0.03091669
c ---     example: kappaf(4.5,34.5,1.e7,3) = -0.21053029
      kappaf=
     &     kappaf1(max(-2.0,min(32.0,t))-toff(kkf),
     &             max(-4.0,min( 4.0,s-soff(kkf))),
     &             prs,kkf)
c
      return
      end
c
      real function kappaf1(t,s,prs,kkf)
c
c --- compressibility coefficient kappa^(theta) from sun et al. (1999)
c
      implicit none
      include 'dimensions.h'
      include 'common_blocks.h'
      real t,s,prs
      integer kkf
      data pref/2.e7/
c
      include 'state_eqn.h'
c
c --- t: potential temperature; s: psu; prs: pressure; kkf: ref.state
c ---     example: kappaf(4.5,34.5,1.e7,1) =  0.11411243
c ---     example: kappaf(4.5,34.5,1.e7,2) =  0.03091669
c ---     example: kappaf(4.5,34.5,1.e7,3) = -0.21053029
      kappaf1=1.e-11*(prs-pref)*
     . ( s*( qs(kkf)+t* qst(kkf) ) +
     .   t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+
css  .       0.5*(prs+pref)*(qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) )/thref
     .       0.5*(prs+pref)*(qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) )/1.e-3
c
      return
      end
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
      subroutine cpy_p(field)
c
c --- exchange information across bering strait seam
c
      implicit none
      include 'dimensions.h'
      include 'bering.h'
c
      real field(idm,jdm),sign
c
c --- exchange p-point values (half grid size away from seam)
      if (beropn) then
        field(iatls,jatl)=field(ipacs,jpac)
        field(ipacn,jpac)=field(iatln,jatl)
      endif
      return
      end
c
c
c> Revision history:
c>
c> Jan. 2001 - modified kappaf to allow nonzero temp.ref.value
c> Nov. 2002 - corrected error in kappaf
c> July 2005 - added sig,dsigdt,dsigds functions for in-situ density ('sigloc')
c
c
      real function hyc_pechg1(delp,sig,nunit)
c
c --- calculate change in ocean's available potential energy due to some
c --- process 'X'. call pechg1  b e f o r e  , pechg2  a f t e r  process X.
c --- input: 3-d arrays of  h y b r i d  layer thickness and density anomaly.
c --- results from pechg1 are stored in 'nunit' for later use by pechg2.
c --- use different values of 'nunit' for nested APE process diagnostics.
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real slithk
      integer nscli,nunit,lgth
      parameter (slithk=1.,nscli=7000./slithk)
      real delp(idm,jdm,kdm),sig(idm,jdm,kdm),
     .     pbfore(idm,jdm,kdm+1),pafter(idm,jdm,kdm+1),
     .     praw(idm,jdm,kdm+1),unused(idm,jdm,kdm)
      real weight(kdm),sliwgt(nscli),slitop,
     .     slibot,slisum,prefbe(kdm),prefaf(kdm),varian(0:kdm),
     .     weightj(jdm),sliwgtj(jdm),varianj(jdm)
      character flnm*60
      data varian(0),varian(kdm)/0.,0./
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
 11   praw(i,j,1)=0.
      do 10 k=1,kk
      do 10 i=ifp(j,l),ilp(j,l)
 10   praw(i,j,k+1)=praw(i,j,k)+delp(i,j,k)
c$OMP END PARALLEL DO
c
c --- transform pressure to isopycnic interface pressure
c
      call reflux(uflx,vflx,sig,praw,
     .            unused,unused,unused,pbfore,theta,kdm,kdm)
c
c --- in preparation for determining the flat reference state (i.e., the
c --- unavailable pot. energy), find total mass above each isopycnic interface
c
      do 1 k=1,kk
      weight(k)=0.
c
c$OMP PARALLEL DO
      do 2 j=1,jj
      weightj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
 2    weightj(j)=weightj(j)+pbfore(i,j,k+1)*scp2(i,j)
c$OMP END PARALLEL DO
c
      do 1 j=1,jj
 1    weight(k)=weight(k)+weightj(j)
c
c --- divide basin into shallow horizontal slices and find mass of each slice
c
      do 3 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      sliwgt(n)=0.
c$OMP PARALLEL DO
      do 4 j=1,jj
      sliwgtj(j)=0.
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    sliwgtj(j)=sliwgtj(j)+scp2(i,j)*(min(pbfore(i,j,kk+1),slibot)-
     .                                 min(pbfore(i,j,kk+1),slitop))
c$OMP END PARALLEL DO
      do 3 j=1,jj
 3    sliwgt(n)=sliwgt(n)+sliwgtj(j)
c
      do 5 k=1,kk-1
c --- add slices vertically until sum exceeds combined mass of layers 1...k.
c --- this tells us where bottom of layer k will be when flattened
c
      slisum=0.
      do 6 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      slisum=slisum+sliwgt(n)
      if (slisum.ge.weight(k)) go to 7
 6    continue
      write (*,*) 'k =',k,'  error: slisum < weight',slisum,weight(k)
      go to 5
c
c --- interpolate among slices to get precise depth of flattened interface
c
 7    prefbe(k)=(slibot*(weight(k)-slisum+sliwgt(n))-
     .           slitop*(weight(k)-slisum          ))/sliwgt(n)
ccc      write (*,*) 'k =',k,'  flattened interface:',prefbe(k)/onem
 5    continue
c
c --- save results for later use by pechg2
c
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 13
 14   continue
      write (lp,*) 'ape --  cannot find slash in',flnmovt
      stop
 13   write (flnm,'(a,i2.2)') flnmovt(1:lgth)//'ape.',nunit
      open (unit=nunit,file=flnm,form='unformatted',status='unknown')
      write (nunit) pbfore,prefbe
      close (nunit)
c
c --- determine variance of interface pressure relative to flat reference state
c
      do 12 k=1,kk-1
      varian(k)=0.
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      varianj(j)=0.
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    varianj(j)=varianj(j)+scp2(i,j)*
     .  (pbfore(i,j,k+1)-min(pbfore(i,j,kk+1),prefbe(k)))**2
c$OMP END PARALLEL DO
c
      do 12 j=1,jj
 12   varian(k)=varian(k)+varianj(j)
c
      hyc_pechg1=0.
      do 8 k=1,kk
 8    hyc_pechg1=hyc_pechg1+
     .  (varian(k)-varian(k-1))/(1000.+theta(k)+thbase)
c
c --- report result in units of joules (kg m^2/sec^2)
      hyc_pechg1=.5*hyc_pechg1/g
      return
      end
c
c
      real function hyc_pechg2(delp,sig,nunit)
c
c --- calculate change in ocean's available potential energy due to some
c --- process 'X'. call pechg1  b e f o r e  , pechg2  a f t e r  process X.
c --- input: 3-d arrays of  h y b r i d  layer thickness and density anomaly.
c --- results from pechg1 representing 'before' state are read from 'nunit'.
c --- use different values of 'nunit' for nested APE process diagnostics.
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real slithk
      integer nscli,nunit,lgth
      parameter (slithk=1.,nscli=7000./slithk)
      real delp(idm,jdm,kdm),sig(idm,jdm,kdm),
     .     pbfore(idm,jdm,kdm+1),pafter(idm,jdm,kdm+1),
     .     praw(idm,jdm,kdm+1),unused(idm,jdm,kdm)
      real weight(kdm),sliwgt(nscli),slitop,
     .     slibot,slisum,prefbe(kdm),prefaf(kdm),varian(0:kdm),
     .     weightj(jdm),sliwgtj(jdm),varianj(jdm)
      character flnm*60
      data varian(0),varian(kdm)/0.,0./
      common /pechng/ pbfore,pafter,prefbe,prefaf
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
 11   praw(i,j,1)=0.
      do 10 k=1,kk
      do 10 i=ifp(j,l),ilp(j,l)
 10   praw(i,j,k+1)=praw(i,j,k)+delp(i,j,k)
c$OMP END PARALLEL DO
c
c --- transform pressure to isopycnic interface pressure
c
      call reflux(uflx,vflx,sig,praw,
     .            unused,unused,unused,pafter,theta,kdm,kdm)
c
c --- in preparation for determining the flat reference state (i.e., the
c --- unavailable pot. energy), find total mass above each isopycnic interface
c
      do 1 k=1,kk
      weight(k)=0.
c
c$OMP PARALLEL DO
      do 2 j=1,jj
      weightj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
 2    weightj(j)=weightj(j)+pafter(i,j,k+1)*scp2(i,j)
c$OMP END PARALLEL DO
c
      do 1 j=1,jj
 1    weight(k)=weight(k)+weightj(j)
c
c --- divide basin into shallow horizontal slices and find mass of each slice
c
      do 3 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      sliwgt(n)=0.
c$OMP PARALLEL DO
      do 4 j=1,jj
      sliwgtj(j)=0.
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    sliwgtj(j)=sliwgtj(j)+scp2(i,j)*(min(pafter(i,j,kk+1),slibot)-
     .                                 min(pafter(i,j,kk+1),slitop))
c$OMP END PARALLEL DO
      do 3 j=1,jj
 3    sliwgt(n)=sliwgt(n)+sliwgtj(j)
c
      do 5 k=1,kk-1
c --- add slices vertically until sum exceeds combined mass of layers 1...k.
c --- this tells us where bottom of layer k will be when flattened
c
      slisum=0.
      do 6 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      slisum=slisum+sliwgt(n)
      if (slisum.ge.weight(k)) go to 7
 6    continue
      write (*,*) 'k =',k,'  error: slisum < weight',slisum,weight(k)
      go to 5
c
c --- interpolate among slices to get precise depth of flattened interface
c
 7    prefaf(k)=(slibot*(weight(k)-slisum+sliwgt(n))-
     .           slitop*(weight(k)-slisum          ))/sliwgt(n)
ccc      write (*,*) 'k =',k,'  flattened interface:',prefaf(k)/onem
 5    continue
c
c --- read previously created file representing 'before' state
c
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 13
 14   continue
      write (lp,*) 'ape --  cannot find slash in',flnmovt
      stop
 13   write (flnm,'(a,i2.2)') flnmovt(1:lgth)//'ape.',nunit
      open (unit=nunit,file=flnm,form='unformatted',status='old')
      read (nunit) pbfore,prefbe
      rewind (nunit)
c
c --- write out 'after' state (for potential future use as new 'before' state)
c
      write (nunit) pafter,prefaf
      close (nunit)
c
c --- determine 'before/after' difference of interface pressure variances
c
      do 12 k=1,kk-1
      varian(k)=0.
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      varianj(j)=0.
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    varianj(j)=varianj(j)+scp2(i,j)*
     . ((pafter(i,j,k+1)-min(pafter(i,j,kk+1),prefaf(k)))**2
     . -(pbfore(i,j,k+1)-min(pbfore(i,j,kk+1),prefbe(k)))**2)
c$OMP END PARALLEL DO
c
      do 12 j=1,jj
 12   varian(k)=varian(k)+varianj(j)
c
      hyc_pechg2=0.
      do 8 k=1,kk
 8    hyc_pechg2=hyc_pechg2
     .  +(varian(k)-varian(k-1))/(1000.+theta(k)+thbase)
c
c --- report result in units of watts (kg m^2/sec^3)
      hyc_pechg2=.5*hyc_pechg2/(g*delt1)
      return
      end
c
c
c> Revision history:
c>
c> Dec. 2004 - fixed bug in loop 9 (excluded interfaces on shallow bottom)
c
c
      subroutine newbot
c
c --- hycom version 0.9
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      real factor
c
c --- call this routine if bottom topography has been altered in any way.
c --- note: removing grid points entirely (by setting their depth to
c --- zero) should be done elsewhere, i.e., before bigrid is called.
c
ccc      if (idm.eq.181 .and. jdm.eq.180)               !  global 2 deg grid
ccc     .   pbot(166,179)=150.*onem
c
c$OMP PARALLEL DO PRIVATE(factor)
      do 3 j=1,jj
      do 3 l=1,isp(j)
c
      do 1 k=1,kk
      do 1 i=ifp(j,l),ilp(j,l)
 1    p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
c
      do 2 i=ifp(j,l),ilp(j,l)
      pbot(i,j)=max(pbot(i,j),botmin*onem)    ! botmin: minimum water depth
      if (abs(p(i,j,kk+1)-pbot(i,j)).gt.onem)
     .  write (lp,'(2i5,a,2f9.1)') i,j,
     .   '  old/new bottom depth:',p(i,j,kk+1)/onem,pbot(i,j)/onem
 2    continue
c
      do 3 k=1,kk
      do 3 i=ifp(j,l),ilp(j,l)
      factor=pbot(i,j)/p(i,j,kk+1)
      if (factor.lt..99999 .or. factor.gt.1.00001) then
        dp(i,j,k   )=dp(i,j,k   )*factor
        dp(i,j,k+kk)=dp(i,j,k+kk)*factor
c
        if (k.gt.1) then
          psikk(i,j,1)=psikk(i,j,1)-p(i,j,k)*(factor-1.)*
     .     (thstar(i,j,k)-thstar(i,j,k-1))*thref
          psikk(i,j,2)=psikk(i,j,2)-p(i,j,k)*(factor-1.)*
     .     (thstar(i,j,kk+k)-thstar(i,j,kk+k-1))*thref
        end if
      end if
 3    continue
c$OMP END PARALLEL DO
c
      return
      end



