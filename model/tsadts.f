      subroutine tsadvc(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 1.0 -- cyclic in j
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real smin,smax,tmin,tmax,sminn,smaxx,tminn,tmaxx,posdef,flxdiv
     .    ,offset,factor,q,pold,pmid,pnew,snew,tnew
      real uflxn(idm,jdm,kdm),vflxn(idm,jdm,kdm),sign(idm,jdm,kdm),
     .     pn(idm,jdm,kdm+1)
      integer kp
      real sigocn,hfharm
      external sigocn,hfharm
c
      do 81 k=1,kk
      km=k+mm
      kn=k+nn
      kp=min(k+1,kk)
c
      smin=999.
      tmin=999.
      smax=-999.
      tmax=-999.
      sminn=999.
      tminn=999.
      smaxx=-999.
      tmaxx=-999.
c
c --- --------------------------------------
c --- advection of thermodynamic variable(s)
c --- --------------------------------------
c
      posdef=100.
c
      call cpy_p(dp(1,1,km))
      call cpy_p(dp(1,1,kn))
c
c$OMP PARALLEL DO PRIVATE(jb,pold,pmid,flxdiv,offset)
      do 49 j=1,jj
      jb=mod(j     ,jj)+1
      do 49 l=1,isp(j)
      do 49 i=ifp(j,l),ilp(j,l)
c
c --- time smoothing of thermodynamic variable(s) (part 1)
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      temp(i,j,km)=temp(i,j,km)*(wts1*pmid+onemm)+
     .             temp(i,j,kn)* wts2*pold
      saln(i,j,km)=saln(i,j,km)*(wts1*pmid+onemm)+
     .             saln(i,j,kn)* wts2*pold
c
c --- before calling 'advem', make sure (a) mass fluxes are consistent
c --- with layer thickness change, and (b) all fields are positive-definite
      flxdiv=(uflx(i+1,j,k)-uflx(i,j,k)
     .       +vflx(i,jb ,k)-vflx(i,j,k))*delt1*scp2i(i,j)
      util2(i,j)=.5*(dpold(i,j,k)+dp(i,j,kn)-flxdiv)
      util1(i,j)=.5*(dpold(i,j,k)+dp(i,j,kn)+flxdiv)
      offset=min(0.,util1(i,j),util2(i,j))
      util2(i,j)=util2(i,j)-offset
      util1(i,j)=util1(i,j)-offset
c
      temp(i,j,kn)=temp(i,j,kn)+posdef
      smin=min(smin,saln(i,j,kn))
      smax=max(smax,saln(i,j,kn))
      tmin=min(tmin,temp(i,j,kn))
      tmax=max(tmax,temp(i,j,kn))
 49   continue
c$OMP END PARALLEL DO
c
      if (tmin.lt.0. .or. smin.lt.0.) then
      do 490 j=1,jj
      do 490 l=1,isp(j)
      do 490 i=ifp(j,l),ilp(j,l)
      if ((tmin.lt.0. and. temp(i,j,kn).eq.tmin) .or.
     .    (smin.lt.0. and. saln(i,j,kn).eq.smin)) then
        write (lp,101) nstep,i,j,k,' neg. temp/saln bfore advem call ',
     .  temp(i,j,kn)-posdef,saln(i,j,kn)
 101  format (i9,' i,j,k =',2i5,i3,a,2f8.2)
        itest=i
        jtest=j
      end if
 490  continue
      call stencl(kn)
      end if
c
      call advem(2,temp(1,1,kn),uflx(1,1,k),vflx(1,1,k),scp2,scp2i,
     .           delt1,util1,util2)
c
      call advem(2,saln(1,1,kn),uflx(1,1,k),vflx(1,1,k),scp2,scp2i,
     .           delt1,util1,util2)
c
c$OMP PARALLEL DO PRIVATE(pold,pmid,pnew)
      do 46 j=1,jj
      do 46 l=1,isp(j)
      do 46 i=ifp(j,l),ilp(j,l)
      temp(i,j,kn)=temp(i,j,kn)-posdef
      sminn=min(sminn,saln(i,j,kn))
      smaxx=max(smaxx,saln(i,j,kn))
      tminn=min(tminn,temp(i,j,kn))
      tmaxx=max(tmaxx,temp(i,j,kn))
c
c --- time smoothing of thickness field
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      pnew=max(0.,dp(i,j,kn))
      dp(i,j,km)=pmid*wts1+(pold+pnew)*wts2
c --- time smoothing of thermodynamic variable(s) (part 2)
      pmid=max(0.,dp(i,j,km))
      temp(i,j,km)=(temp(i,j,km)+temp(i,j,kn)*wts2*pnew)/
     .   (pmid+onemm)
      saln(i,j,km)=(saln(i,j,km)+saln(i,j,kn)*wts2*pnew)/
     .   (pmid+onemm)
      th3d(i,j,km)=sigocn(temp(i,j,km),saln(i,j,km))
c
c --- build up time integral of mass field variables
      pmid=max(0.,dp(i,j,km))
      dpav (i,j,k)=dpav (i,j,k)+pmid
      temav(i,j,k)=temav(i,j,k)+temp(i,j,km)*pmid
      salav(i,j,k)=salav(i,j,k)+saln(i,j,km)*pmid
      th3av(i,j,k)=th3av(i,j,k)+th3d(i,j,km)*pmid
 46   continue
c$OMP END PARALLEL DO
c
      if (tminn+posdef.lt.0. .or. sminn.lt.0.) then
      do 492 j=1,jj
      do 492 l=1,isp(j)
      do 492 i=ifp(j,l),ilp(j,l)
      if (temp(i,j,kn).eq.tminn .or. saln(i,j,kn).eq.sminn)
     .  write (lp,101) nstep,i,j,k,' neg. temp/saln after advem call ',
     .  temp(i,j,kn),saln(i,j,kn)
 492  continue
      end if
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' th,s,dp after advection  '',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,kn),saln(itest,jtest,kn),
cdiag.dp(itest,jtest,kn)/onem
c
      if (diagno)
     . write (lp,'(i9,i3,'' min/max of s after advection:'',4f7.2)')
     . nstep,k,sminn,smaxx
c
c --- --------------------------------------
c --- diffusion of thermodynamic variable(s)
c --- --------------------------------------
c
      call cpy_p(temp(1,1,kn))
      call cpy_p(saln(1,1,kn))
      call cpy_p(th3d(1,1,kn))
c
c$OMP PARALLEL DO PRIVATE(ja,factor)
      do 145 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 144 l=1,isu(j)
      do 144 i=ifu(j,l),ilu(j,l)
      factor=scuy(i,j)*2.*hfharm(max(dp(i-1,j,kn),onemm)
     .                          ,max(dp(i  ,j,kn),onemm))
      uflux (i,j)=factor*(temp(i-1,j,kn)-temp(i,j,kn))
 144  uflux2(i,j)=factor*(saln(i-1,j,kn)-saln(i,j,kn))
c
      do 145 l=1,isv(j)
      do 145 i=ifv(j,l),ilv(j,l)
      factor=scvx(i,j)*2.*hfharm(max(dp(i,ja ,kn),onemm)
     .                          ,max(dp(i,j  ,kn),onemm))
      vflux (i,j)=factor*(temp(i,ja ,kn)-temp(i,j,kn))
 145  vflux2(i,j)=factor*(saln(i,ja ,kn)-saln(i,j,kn))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jb,factor)
      do 146 j=1,jj
      jb=mod(j     ,jj)+1
      do 146 l=1,isp(j)
      do 146 i=ifp(j,l),ilp(j,l)
      factor=-temdff*delt1/(scp2(i,j)*max(dp(i,j,kn),onemm))
      util1(i,j)=(uflux (i+1,j)-uflux (i,j)
     .           +vflux (i,jb )-vflux (i,j))*factor
      util2(i,j)=(uflux2(i+1,j)-uflux2(i,j)
     .           +vflux2(i,jb )-vflux2(i,j))*factor
      temp(i,j,kn)=temp(i,j,kn)+util1(i,j)
      saln(i,j,kn)=saln(i,j,kn)+util2(i,j)
      th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag. write (lp,100) nstep,i,j,k,'t,s,dt,ds,dsigdt,dsigds,cabbl =',
cdiag. temp(i,j,kn),saln(i,j,kn),util1(i,j),util2(i,j),
cdiag. dsigdt(temp(i,j,kn),saln(i,j,kn))*util1(i,j),
cdiag. dsigds(temp(i,j,kn),saln(i,j,kn))*util2(i,j),cabbl(i,j,k)
c
 146  continue
c$OMP END PARALLEL DO
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' t,s,dp after isopyc.mix.'',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,kn),saln(itest,jtest,kn),
cdiag.dp(itest,jtest,kn)/onem
c
      call cpy_p(temp(1,1,km))
      call cpy_p(temp(1,1,kn))
      call cpy_p(saln(1,1,km))
      call cpy_p(saln(1,1,kn))
      call cpy_p(th3d(1,1,km))
      call cpy_p(th3d(1,1,kn))
 81   continue
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
c$OMP PARALLEL DO PRIVATE(ja)
      do 155 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 155 k=1,kk
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
      uflxav(i,j,k)=uflxav(i,j,k)+uflxn(i,j,k)		!  uflx time integral
     .   *.5*min(nstep,2)
 154  continue
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
      vflxav(i,j,k)=vflxav(i,j,k)+vflxn(i,j,k)		!  vflx time integral
     .   *.5*min(nstep,2)
 155  continue
c$OMP END PARALLEL DO
c
      return
      end
c
c
c  Revision history:
c
c> June 1995 - eliminated setting of salinity in massless layers (loop 46)
c>             (this is now done in mxlayr.f)
c> Aug. 1995 - added array -cabbl- to transmit cabbeling info to -diapfl-
c> Aug. 1995 - omitted t/s/dp time smoothing in case of abrupt mxlayr.thk.change
c> Sep. 1995 - increased temdff if mixed layer occupies >90% of column
