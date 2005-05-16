      subroutine inicon
c
c --- hycom version 0.9
      USE FLUXES, only : e0,prec,evapor,flowo,eflowo,dmua,dmva
     .      ,erunosi,runosi,runpsi,dmui,dmvi,dmsi,dhsi,dssi
      USE MODEL_COM, only : focean
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'state_eqn.h'
      include 'a2o.h'
      include 'cpl.h'
c
      integer totlj(jdm,kdm-1),totl(kdm-1),iz,jz,ni
      character text*24,preambl(5)*79
      real tofsig,kappaf,sigocn,cold,temavg,vol,sst,sofsig
      external tofsig,kappaf,sigocn,sofsig
c
c --- set minimum salinity for each isopycnic layer
      cold=-2.5
      do 13 k=2,kk
 13   salmin(k)=sofsig(sigma(k),cold)
      print *,' chk salmin=',salmin
c
c --- subtract constant 'thbase' from sigma to reduce roundoff errors
c
      do 14 k=1,kk
 14   theta(k)=sigma(k)-thbase
c
      if (nstep0.le.0) then                ! start from Levitus
        delt1=baclin
c --- -------------------------
c --- mass field initialization
c --- -------------------------
c
c --- read mixed-layer temperature
c
      write (lp,'(2a)') 'get initial temperature from  ',flnmint
      open(unit=32,file=flnmint,form='formatted',status='old',
     .     action='read')
c
      read (32,'(a79)') (preambl(n),n=1,5)
      write (lp,'(a79)') (preambl(n),n=1,5)
      read (32,'(10f8.4)') ((temp(i,j,1),i=1,idm),j=1,jdm)
      close (32)
c     write (lp,100) 'temp field read, layer',1
c     call zebra(temp,idm,ii1,jj)
 100  format (a,i4)
      close (unit=32)
c
c --- read salinity
c
      write (lp,'(2a)') 'get initial salinity from  ',flnmins
      open(unit=32,file=flnmins,form='formatted',status='old',
     .     action='read')
c
      read (32,'(a79)') (preambl(n),n=1,5)
      write (lp,'(a79)') (preambl(n),n=1,5)
      read(32,'(10f8.4)') (((saln(i,j,k),i=1,idm),j=1,jdm),k=1,kk)
      close (32)
c     write (lp,100) 'saln field read, layers, 1 -',kk
c     call zebra(saln(1,1,1),idm,ii1,jj)
c
c --- read interface pressure
c
      write (lp,'(2a)') 'get initial pressure from  ',flnminp
      open(unit=32,file=flnminp,form='formatted',status='old',
     .     action='read')
c
      read (32,'(a79)') (preambl(n),n=1,5)
      write (lp,'(a79)') (preambl(n),n=1,5)
      read(32,'(10f8.2)') (((p(i,j,k+1),i=1,idm),j=1,jdm),k=1,kk)
      close (32)
c     write (lp,100) 'pres field read, levels 2 -',kk+1
c     call zebra(p(1,1,kk+1),idm,ii1,jj)
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
c
      do 15 i=ifp(j,l),ilp(j,l)
      th3d(i,j,1)=sigocn(temp(i,j,1),saln(i,j,1))-thbase
 15   p(i,j,1)=0.
c
      do 9 k=1,kk
      do 9 i=ifp(j,l),ilp(j,l)
      p(i,j,k+1)=max(p(i,j,k),min(depths(i,j),p(i,j,k+1))*onem)
      dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
      if (k.gt.1) then
        th3d(i,j,k)=theta(k)
        saln(i,j,k)=max(saln(i,j,k),salmin(k))
        temp(i,j,k)=tofsig(theta(k)+thbase,saln(i,j,k))
      end if
 9    continue
c
      do 17 i=ifp(j,l),ilp(j,l)
 17   pbot(i,j)=p(i,j,kk+1)
c
 10   continue
c$OMP END PARALLEL DO
c
      write (lp,'('' sigma(k)     :'',9f7.2/(15x,9f7.2))')
     .   (sigma(k),k=1,kk)
c
cdiag do k=1,kk,3
cdiag write (text,'(''intf.pressure (m), k='',i3)') k+1
cdiag call prtmsk(ip,p(1,1,k+1),util1,idm,ii1,jj,0.,1./onem,text)
cdiag end do
c
      call dpthuv
c     call convec(1,1,0,0,1,1)
c
c$OMP PARALLEL DO
      do 11 j=1,jj
c
      do 12 k=1,kk-1
 12   totlj(j,k)=0
c
      do 11 k=1,kk
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
      tracer(i,j,k)=0.
      dp  (i,j,k+kk)=dp  (i,j,k)
      th3d(i,j,k+kk)=th3d(i,j,k)
      temp(i,j,k+kk)=temp(i,j,k)
      saln(i,j,k+kk)=saln(i,j,k)
      thermb(i,j,k)=kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k))
     .             *(1000.+th3d(i,j,k)+thbase)
      thermb(i,j,k+kk)=thermb(i,j,k)
      thstar(i,j,k)=th3d(i,j,k)+thermb(i,j,k)
c
      if (i.eq.itest.and.j.eq.jtest)
     . write (lp,'(2i5,i3,a,3f7.3,3x,2f7.3,f8.1)')
     .  i,j,k,'  dens,thstar,kappa=',th3d(i,j,k)+thbase,thstar(i,j,k)
     .   +thbase,kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k)),
     .    temp(i,j,k),saln(i,j,k),p(i,j,k)/onem
c
      if (k.gt.1 .and. thstar(i,j,k).lt.thstar(i,j,k-1))
     .   totlj(j,k-1)=totlj(j,k-1)+1
 11   continue
c$OMP END PARALLEL DO
c
      do 18 k=1,kk-1
      totl(k)=0
      do 18 j=1,jj
 18   totl(k)=totl(k)+totlj(j,k)
      write (lp,'(a/(10i7))') 'static instability count by layer:',
     .  totl
c
c$OMP PARALLEL DO
      do 50 j=1,jj
      do 50 l=1,isp(j)
      do 50 i=ifp(j,l),ilp(j,l)
      montg(i,j,1)=0.
c
      do 52 k=1,kk-1
 52   montg(i,j,k+1)=montg(i,j,k)-p(i,j,k+1)*(thstar(i,j,k+1)-
     .                                        thstar(i,j,k  ))*thref**2
c
      thkk(i,j)=thstar(i,j,kk)
 50   psikk(i,j)=montg(i,j,kk)
c$OMP END PARALLEL DO
c
      else				!  nstep0 > 0
c
c --- start from restart file prescribed
c
      write (lp,111) nstep0,time0
 111  format (9x,'chk time step in restart file -',i9,5x,' day ',f9.2)
c
      delt1=baclin+baclin
c
      call dpthuv
c
      do 16 m=1,2
      mm=(m-1)*kk
c
c$OMP PARALLEL DO
      do 19 j=1,jj
      do 19 k=1,kk
      do 19 l=1,isp(j)
      do 19 i=ifp(j,l),ilp(j,l)
 19   p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c$OMP END PARALLEL DO
c
      call dpudpv(mm)
 16   continue
c
cc$OMP PARALLEL DO
c     do 21 j=1,jj
c     do 21 l=1,isp(j)
c     do 21 i=ifp(j,l),ilp(j,l)
css   omlhc(i,j)=spcifh*p(i,j,2)/(onem *thref)           ! J/m*m C
c     thermb(i,j,k)=kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k))
c    .             *(1000.+th3d(i,j,k)+thbase)
c     thermb(i,j,k+kk)=thermb(i,j,k)
c21   continue
cc$OMP END PARALLEL DO
      write(800,'(a,i5)') 'restart',nstep0
      write(800,'(8f13.8)') ((dp(115,170,k)/onem),k=1,2*kk)
      end if				!  nstep0 > 0  or  = 0
c
      if (itest.gt.0.and.jtest.gt.0) write (lp,103) nstep,itest,jtest,
     .  '  init.profile  temp    saln  thstar   thkns    dpth   montg',
     .  (k,temp(itest,jtest,k),saln(itest,jtest,k),
     .  thstar(itest,jtest,k)+thbase,dp(itest,jtest,k)/onem,
     .  p(itest,jtest,k+1)/onem,montg(itest,jtest,k)/g,k=1,kk)
 103  format (i9,2i5,a/(28x,i3,2f8.2,f8.2,2f8.1,f8.3))
c
c
c     print *,' inicon gtemp'
c     call zebra(utila,iia,iia,jja)
c     write(*,'(10f7.2)') ((gtemp(1,1,i,j),i=20,29),j=37,46)
c     print *,' inicon sss'
c     call zebra(sss,iia,iia,jja)
c     write(*,'(10f7.2)') ((sss(i,j),i=20,29),j=37,46)
c     print *,' inicon mlhc'
c     call zebra(mlhc,iia,iia,jja)
c     print *,' sst/sss=',gtemp(1,1,24,4),sss(24,4),' %=',focean(24,4)
c
      return
      end
c
c> Revision history:
c>
c> Mar. 2000 - conversion to SI units
c> Aug. 2000 - added diagnostic count of static instabilities
c> Apr. 2001 - eliminated stmt_funcs.h
