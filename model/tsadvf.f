      subroutine tsadvc(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.3 -- th3d/spiciness advection
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      integer iz,jz,kp
      real flxdiv,spcoef,offset,factor,q,pold,pmid,pnew,snew,tnew,
ccc     .     cortmj(jdm),corsmj(jdm),vlummj(jdm),corrtm,corrsm,vlumm,
     .     cortnj(jdm),corsnj(jdm),vlumnj(jdm),corrtn,corrsn,vlumn,
     .     wlist1(idm,jdm),wlist2(idm,jdm),aux(idm,jdm),tscnsv,fac0
      real uflxn(idm,jdm,kdm),vflxn(idm,jdm,kdm),sign(idm,jdm,kdm),
     .     pn(idm,jdm,kdm+1),spice(idm,jdm,2)
      real sigocn,hfharm,tofsig,tofspi
      external sigocn,hfharm,tofsig,tofspi
      data spcoef/0.13/                !  spiciness coeff.
      logical globcor
      data globcor/.true./                !  transport error correction switch
c
c --- choose whether to conserve T/S (tscnsv=1) or rho/S (tscnsv=0) during
c --- lateral mixing. compromise settings (0 < tscnsv < 1) are allowed
      data tscnsv/.85/
c
c$OMP PARALLEL DO
      do 6 j=1,jj
      do 6 i=1,ii
 6    aux(i,j)=0.                                                ! diapyc.flux
c$OMP END PARALLEL DO
c
      do 1 k=1,kk
      km=k+mm
      kn=k+nn
      kp=min(k+1,kk)
c
c --- smooth mass fluxes in lateral direction
c$OMP PARALLEL DO PRIVATE(ia,ib,ja,jb)
      do 790 j=1,jj
      do 791 l=1,isv(j)
      do 791 i=ifv(j,l),ilv(j,l)
      vflux(i,j)=vflx(i,j,k)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ia=max( 1,i-1)
      ib=min(ii,i+1)
      if (i.eq.1  .or. iv(ia,j).eq.0) ia=i
      if (i.eq.ii .or. iv(ib,j).eq.0) ib=i
      vflux(i,j)=.5*vflx(i,j,k)+.25*(vflx(ia,j,k)+vflx(ib,j,k))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 791  continue
c
      do 790 l=1,isu(j)
      do 790 i=ifu(j,l),ilu(j,l)
      uflux(i,j)=uflx(i,j,k)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      if (iu(i,ja).eq.0) ja=j
      if (iu(i,jb).eq.0) jb=j
      uflux(i,j)=.5*uflx(i,j,k)+.25*(uflx(i,ja,k)+uflx(i,jb,k))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 790  continue
c$OMP END PARALLEL DO
c
c --- --------------------------------------
c --- advection of thermodynamic variable(s)
c --- --------------------------------------
c
c$OMP PARALLEL DO PRIVATE(jb,pold,pmid,flxdiv,offset,iz,jz)
c$OMP. SHARED(spcoef)
      do 2 j=1,jj
      jb=mod(j     ,jj)+1
ccc      corsmj(j)=0.
ccc      cortmj(j)=0.
      corsnj(j)=0.
      cortnj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
c
      if (globcor) then
ccc        cortmj(j)=cortmj(j)+temp(i,j,km)*dp(i,j,km)*scp2(i,j)
ccc        corsmj(j)=corsmj(j)+saln(i,j,km)*dp(i,j,km)*scp2(i,j)
        cortnj(j)=cortnj(j)+temp(i,j,kn)*dpold(i,j,k)*scp2(i,j)
        corsnj(j)=corsnj(j)+saln(i,j,kn)*dpold(i,j,k)*scp2(i,j)
      end if
c
c --- define 2nd thermodyn.variable (to be advected together with density)
      spice(i,j,m)=saln(i,j,km)+spcoef*temp(i,j,km)
      spice(i,j,n)=saln(i,j,kn)+spcoef*temp(i,j,kn)
      util3(i,j)=th3d(i,j,kn)
c
c --- time smoothing of thermodynamic variable(s) (part 1)
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      th3d(i,j,km)=th3d(i,j,km)*(wts1*pmid+epsil)+
     .             th3d(i,j,kn)* wts2*pold
      spice(i,j,m)=spice(i,j,m)*(wts1*pmid+epsil)+
     .             spice(i,j,n)* wts2*pold
c
c --- before calling 'advfct', make sure mass fluxes are consistent
c --- with layer thickness change
      flxdiv=(uflux(i+1,j)-uflux(i,j)
     .       +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
      util2(i,j)=.5*(dpold(i,j,k)+dp(i,j,kn)-flxdiv)
      util1(i,j)=.5*(dpold(i,j,k)+dp(i,j,kn)+flxdiv)
      offset=min(0.,util1(i,j),util2(i,j))
      util2(i,j)=util2(i,j)-offset
      util1(i,j)=util1(i,j)-offset
 2    continue
c$OMP END PARALLEL DO
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag. write (lp,101) nstep,itest,jtest,k,' th,t,s,dp bfore advec ',
cdiag.  th3d(itest,jtest,kn)+thbase,temp(itest,jtest,kn),
cdiag.   saln(itest,jtest,kn),dp(itest,jtest,kn)/onem
 101  format (i9,2i5,i3,a,3f8.3,f9.3)
c
      call advfct(2,th3d(1,1,kn),uflux,vflux,scp2,scp2i,
     .            delt1,util1,util2)
c
      call advfct(2,spice(1,1,n),uflux,vflux,scp2,scp2i,
     .            delt1,util1,util2)
c
c$OMP PARALLEL DO PRIVATE(pold,pmid,pnew,snew,tnew) SHARED(k,spcoef)
      do 3 j=1,jj
      do 3 l=1,isp(j)
      do 3 i=ifp(j,l),ilp(j,l)
c
c --- time smoothing of thickness field
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      pnew=max(0.,dp(i,j,kn))
      dp(i,j,km)=pmid*wts1+(pold+pnew)*wts2
      p(i,j,k+1)=p(i,j,k)+dp(i,j,km)                           ! needed in reflux
      aux(i,j)=aux(i,j)+(dp(i,j,km)-pmid)                   ! diapyc.flux
      diaflx(i,j,k )=diaflx(i,j,k )+aux(i,j)                   ! diapyc.flux
      diaflx(i,j,kp)=diaflx(i,j,kp)-aux(i,j)                   ! diapyc.flux
c
c --- time smoothing of thermodynamic variable(s) (part 2)
      pmid=max(0.,dp(i,j,km))
      th3d(i,j,km)=(th3d(i,j,km)+th3d(i,j,kn)*wts2*pnew)/
     .   (pmid+epsil)
      spice(i,j,m)=(spice(i,j,m)+spice(i,j,n)*wts2*pnew)/
     .   (pmid+epsil)
c
c --- recover t/s from density/spiciness
      if (dp(i,j,km).gt.onemm) then
        temp(i,j,km)=tofspi(th3d(i,j,km)+thbase,spice(i,j,m),spcoef)
        saln(i,j,km)=spice(i,j,m)-spcoef*temp(i,j,km)
      end if
      if (dp(i,j,kn).gt.onemm) then
        tnew=tofspi(th3d(i,j,kn)+thbase,spice(i,j,n),spcoef)
        snew=spice(i,j,n)-spcoef*tnew
        if (tnew.lt.-9. and. dp(i,j,kn).gt.onecm)
     .   write (lp,'(2i5,i3,a,4f8.2/36x,a,4f8.2)') i,j,k,
     .    ' bad th,t,s,dp combo   old:',util3(i,j)+thbase,
     .     temp(i,j,kn),saln(i,j,kn),dpold(i,j,k)/onem,'new:',
     .      th3d(i,j,kn)+thbase,tnew,snew,dp(i,j,kn)/onem
        temp(i,j,kn)=tnew
        saln(i,j,kn)=snew
      end if
c
c --- build up time integral of mass field variables
      pmid=max(0.,dp(i,j,km))
      dpav (i,j,k)=dpav (i,j,k)+pmid
      temav(i,j,k)=temav(i,j,k)+temp(i,j,km)*pmid
      salav(i,j,k)=salav(i,j,k)+saln(i,j,km)*pmid
      th3av(i,j,k)=th3av(i,j,k)+th3d(i,j,km)*pmid
 3    continue
c$OMP END PARALLEL DO 
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag. write (lp,101) nstep,itest,jtest,k,' th,t,s,dp after advec ',
cdiag.  th3d(itest,jtest,kn)+thbase,temp(itest,jtest,kn),
cdiag.   saln(itest,jtest,kn),dp(itest,jtest,kn)/onem
c
c --- --------------------------------------
c --- diffusion of thermodynamic variable(s)
c --- --------------------------------------
c
c$OMP PARALLEL DO PRIVATE(ja,factor)
      do 145 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 144 l=1,isu(j)
      do 144 i=ifu(j,l),ilu(j,l)
      factor=scuy(i,j)*2.*hfharm(max(dp(i-1,j,kn),onemm)
     .                          ,max(dp(i  ,j,kn),onemm))
      uflux (i,j)=factor*(saln(i-1,j,kn)-saln(i,j,kn))
      uflux2(i,j)=factor*(temp(i-1,j,kn)-temp(i,j,kn))
 144  uflux3(i,j)=factor*(th3d(i-1,j,kn)-th3d(i,j,kn))
c
      do 145 l=1,isv(j)
      do 145 i=ifv(j,l),ilv(j,l)
      factor=scvx(i,j)*2.*hfharm(max(dp(i,ja ,kn),onemm)
     .                          ,max(dp(i,j  ,kn),onemm))
      vflux (i,j)=factor*(saln(i,ja ,kn)-saln(i,j,kn))
      vflux2(i,j)=factor*(temp(i,ja ,kn)-temp(i,j,kn))
 145  vflux3(i,j)=factor*(th3d(i,ja ,kn)-th3d(i,j,kn))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jb,factor) SHARED(equatn)
      do 146 j=1,jj
      jb=mod(j     ,jj)+1
      do 146 l=1,isp(j)
      do 146 i=ifp(j,l),ilp(j,l)
      factor=-temdff*delt1/(scp2(i,j)*max(dp(i,j,kn),onemm))
     .     *max(1.,6.-2.*abs(i-equatn))          ! enhance EQ temdff
      saln(i,j,kn)=saln(i,j,kn)+(uflux (i+1,j)-uflux (i,j)
     .                          +vflux (i,jb )-vflux (i,j))*factor
      temp(i,j,kn)=temp(i,j,kn)+(uflux2(i+1,j)-uflux2(i,j)
     .                          +vflux2(i,jb )-vflux2(i,j))*factor
      th3d(i,j,kn)=th3d(i,j,kn)+(uflux3(i+1,j)-uflux3(i,j)
     .                          +vflux3(i,jb )-vflux3(i,j))*factor
c
c --- reconcile -temp- and -th3d-
      if (tscnsv.gt.0.)
     .  th3d(i,j,kn)=tscnsv*(sigocn(temp(i,j,kn),saln(i,j,kn))-thbase)
     .          +(1.-tscnsv)*th3d(i,j,kn)
      if (tscnsv.lt.1.)
     .  temp(i,j,kn)=tofsig(th3d(i,j,kn)+thbase,saln(i,j,kn))
 146  continue
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag. write (lp,101) nstep,itest,jtest,k,' th,t,s,dp after mixing',
cdiag.  th3d(itest,jtest,kn)+thbase,temp(itest,jtest,kn),
cdiag.   saln(itest,jtest,kn),dp(itest,jtest,kn)/onem
c$OMP END PARALLEL DO
c
      if (globcor) then
c --- do global correction on advected fields
c
c$OMP PARALLEL DO
        do 15 j=1,jj
ccc        vlummj(j)=0.
        vlumnj(j)=0.
        do 15 l=1,isp(j)
        do 15 i=ifp(j,l),ilp(j,l)
ccc        cortmj(j)=cortmj(j)-temp(i,j,km)*dp(i,j,km)*scp2(i,j)
ccc        corsmj(j)=corsmj(j)-saln(i,j,km)*dp(i,j,km)*scp2(i,j)
ccc        vlummj(j)=vlummj(j)+dp(i,j,km)*scp2(i,j)
        cortnj(j)=cortnj(j)-temp(i,j,kn)*dp(i,j,kn)*scp2(i,j)
        corsnj(j)=corsnj(j)-saln(i,j,kn)*dp(i,j,kn)*scp2(i,j)
 15     vlumnj(j)=vlumnj(j)+dp(i,j,kn)*scp2(i,j)
c$OMP END PARALLEL DO
c
ccc        corrtm=0.
ccc        corrsm=0.
ccc        vlumm=0.
        corrtn=0.
        corrsn=0.
        vlumn=0.
        do 16 j=1,jj
ccc        corrtm=corrtm+cortmj(j)
ccc        corrsm=corrsm+corsmj(j)
ccc        vlumm=vlumm+vlummj(j)
        corrtn=corrtn+cortnj(j)
        corrsn=corrsn+corsnj(j)
 16     vlumn=vlumn+vlumnj(j)
c
        if (vlumn.gt.area*onecm) then
ccc          corrtm=corrtm/vlumm
ccc          corrsm=corrsm/vlumm
          corrtn=corrtn/vlumn
          corrsn=corrsn/vlumn
        else
ccc          corrtm=0.
ccc          corrsm=0.
          corrtn=0.
          corrsn=0.
        end if
c
c$OMP PARALLEL DO
        do 17 j=1,jj
        do 17 l=1,isp(j)
        do 17 i=ifp(j,l),ilp(j,l)
ccc        temp(i,j,km)=temp(i,j,km)+corrtm
ccc        saln(i,j,km)=saln(i,j,km)+corrsm
ccc        th3d(i,j,km)=sig(temp(i,j,km),saln(i,j,km))-thbase
        temp(i,j,kn)=temp(i,j,kn)+corrtn
        saln(i,j,kn)=saln(i,j,kn)+corrsn
 17     th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))-thbase
c$OMP END PARALLEL DO
c
      end if                                !  globcor = .true.
c
 1    continue
c
c --- convert mass fluxes to density coord. prior to time integration
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call reflux(uflx ,vflx ,th3d(1,1,k1m),p,
     .            uflxn,vflxn,sign,pn,theta,kdm,kdm)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- activate this loop if -reflux- is   n o t   called
ccc      do k=1,kk
ccc      do j=1,jj
ccc      do i=1,ii
ccc      uflxn(i,j,k)=uflx(i,j,k)
ccc      vflxn(i,j,k)=vflx(i,j,k)
ccc      end do
ccc      end do
ccc      end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      fac0=.5*delt1/baclin
c$OMP PARALLEL DO PRIVATE(ja)
      do 155 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 155 k=1,kk
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
      uflxav(i,j,k)=uflxav(i,j,k)+uflxn(i,j,k)                !  uflx time integral
css  .   *.5*min(nstep,2)
     .   *fac0
 154  continue
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
      vflxav(i,j,k)=vflxav(i,j,k)+vflxn(i,j,k)                !  vflx time integral
css  .   *.5*min(nstep,2)
     .   *fac0
 155  continue
c$OMP END PARALLEL DO
c
      return
      end
c
      real function tofspi(sigm,spice,coeff)
c
c --- temp (deg c) as a function of sigma and spiciness
c --- (spiciness defined as  S + coeff * T )
c
      implicit none
      real sigm,salin,sqq,q,a0,a1,a2,cubr,cubq,cuban,cubrl,cubim,athird
     .    ,spice,coeff,sigocn
      external sigocn
      parameter (athird=1./3.)
c
      include 'state_eqn.h'
c
c --- equation of state with 7 coefficients
ccc      q=1./(c6-coeff*c7)
ccc      a0=(c1         +c3*spice)*q
ccc      a1=(c2-coeff*c3+c5*spice)*q
ccc      a2=(c4-coeff*c5+c7*spice)*q
c
c --- equation of state with 9 coefficients
      q=1./(c6-coeff*(c7-                        coeff*c9))
      a0=(c1+spice*(c3+   spice*c8)                       )*q
      a1=(c2-coeff*(c3+2.*spice*c8)+spice*(c5+   spice*c9))*q
      a2=(c4-coeff*(c5-   coeff*c8)+spice*(c7-2.*coeff*c9))*q
c
      cubq=athird*a1-(athird*a2)**2
      cubr=athird*(.5*a1*a2-1.5*(a0-sigm*q))-(athird*a2)**3
c
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
c
      cuban=athird*atan2(sqrt(max(0.,-(cubq**3+cubr**2))),cubr)
      sqq=sqrt(-cubq)
      cubrl=sqq*cos(cuban)
      cubim=sqq*sin(cuban)
      tofspi=-cubrl+sqrt(3.)*cubim-athird*a2
cdiag salin=spice-coeff*tofspi
cdiag if (abs(sigocn(tofspi,salin)-sigm).gt.1.e-9) write (*,100)
cdiag. tofspi,salin,sigm,sigocn(tofspi,salin)
 100  format ('tofspi,sal,old/new sig =',4f11.6)
      return
      end
c
c
c  Revision history:
c
c> June 1995 - eliminated setting of salinity in massless layers (loop 46)
c>             (this is now done in mxlayr.f)
c> Aug. 1995 - omitted t/s/dp time smoothing in case of abrupt mxlayr.thk.change
c> Sep. 1995 - increased temdff if mixed layer occupies >90% of column
c> June 2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Nov. 2000 - added code to make transport globally conservative (globcor)
c> Mar. 2001 - added computation of mass flux time averages
c> Mar. 2001 - converted mass fluxes to density coord. prior to time integr.
c> Apr. 2001 - eliminated stmt_funcs.h
c> Mar. 2002 - added tracer advection
c> Sep. 2002 - added time-averaging of  T/S/dp
c> Sep. 2003 - removed tracer advection (now in trcadv.f)
c> Oct. 2003 - choice between T/S and rho/S conservation during lat. mixing
c> Dec. 2004 - replaced salinity by spiciness as advected variable
c> Dec. 2004 - switched global correction ('globcor') from rho/spice to T/S
c> Jan. 2005 - replaced salinity by spiciness as time-smoothed variable
