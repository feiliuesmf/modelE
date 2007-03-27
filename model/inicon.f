      subroutine inicon
c
c --- hycom version 0.9
      USE FLUXES, only : e0,prec,evapor,flowo,eflowo,dmua,dmva
     .      ,erunosi,runosi,runpsi,dmui,dmvi,dmsi,dhsi,dssi
     .      ,gtemp,sss,mlhc
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
      integer totlj(jdm,kdm-1),totl(kdm-1),iz,jz,ni,ipo(idm,jdm),kkap
      character text*24,preambl(5)*79
      real tofsig,kappaf,sigocn,cold,temavg,vol,sst,sofsig,spval
      real*4 real4(iold,jdm)
      external tofsig,kappaf,sigocn,sofsig
      data spval/-.03125/
c
c --- set minimum salinity for each isopycnic layer
      cold=-2.5
      do 13 k=2,kk
 13   salmin(k)=sofsig(theta(k),cold)
      write(*,'(a/(10f6.2))') 'chk salmin(2:kk)=',salmin(2:kk)
c
      if (nstep0.eq.0) then                ! start from Levitus
        call geopar
        delt1=baclin
c --- -------------------------
c --- mass field initialization
c --- -------------------------
c
c --- read salinity
c
      write (lp,'(2a)') 'get initial salinity from  ',flnmins
      open(unit=32,file=flnmins,form='formatted',status='old',
     .     action='read')
c
      read (32,'(a79)') (preambl(n),n=1,5)
      write (lp,'(a79)') (preambl(n),n=1,5)
      do k=1,kk
      read(32,'(10f8.4)') ((util1(i,j),i=1,iold),j=1,jdm)
c
      if (k.eq.1) then
        do i=1,iold
        do j=1,jdm
        ipo(i,j)=0.
        if (abs(util1(i,j)-spval).gt.2.) ipo(i,j)=1    ! define old ip
        enddo
        enddo
      endif
c
      call refinp(equato,equatn,ipo,util1,saln(1,1,k))
      if (k.eq.1) 
     .    call prtmsk(ip,saln(1,1,k),util1,idm,139,jj,0.,1.,'s_rfn')
      enddo
      close (32)
      write (lp,100) 'saln field read, layers, 1 -',kk
      call zebra(saln(1,1,1),idm,ii1,jj)
c
c --- read mixed-layer temperature
c
      write (lp,'(2a)') 'get initial temperature from  ',flnmint
      open(unit=32,file=flnmint,form='formatted',status='old',
     .     action='read')
c
      read (32,'(a79)') (preambl(n),n=1,5)
      write (lp,'(a79)') (preambl(n),n=1,5)
      read (32,'(10f8.4)') ((util1(i,j),i=1,iold),j=1,jdm)
      close (32)
 100  format (a,i4)
      close (unit=32)
      write(*,*) 'chk equator =',equato,equatn
      call refinp(equato,equatn,ipo,util1,temp(1,1,1))
      write (lp,100) 'ML temp field '
      call zebra(temp,idm,ii1,jj)
c
c --- read interface pressure
c
      write (lp,'(2a)') 'get initial pressure from  ',flnminp
      open(unit=32,file=flnminp,form='formatted',status='old',
     .     action='read')
c
      read (32,'(a79)') (preambl(n),n=1,5)
      write (lp,'(a79)') (preambl(n),n=1,5)
      do k=1,kk
      read(32,'(10f8.2)') ((util1(i,j),i=1,iold),j=1,jdm)
      call refinp(equato,equatn,ipo,util1,p(1,1,k+1))
      enddo
      close (32)
c     write (lp,100) 'pres field read, levels 2 -',kk+1
c     call zebra(p(1,1,kk+1),idm,ii1,jj)
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
c
      do 15 i=ifp(j,l),ilp(j,l)
      th3d(i,j,1)=sigocn(temp(i,j,1),saln(i,j,1))
 15   p(i,j,1)=0.
c
      do 9 k=1,kk
      do 9 i=ifp(j,l),ilp(j,l)
      p(i,j,k+1)=max(p(i,j,k),min(depths(i,j),p(i,j,k+1))*onem)
      dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
      if (k.gt.1) then
        th3d(i,j,k)=theta(k)
        saln(i,j,k)=max(saln(i,j,k),salmin(k))
        temp(i,j,k)=tofsig(theta(k),saln(i,j,k))
      end if
 9    continue
c
      do 17 i=ifp(j,l),ilp(j,l)
 17   pbot(i,j)=p(i,j,kk+1)
c
 10   continue
c$OMP END PARALLEL DO
c
      write (lp,'('' theta(k)     :'',9f7.2/(15x,9f7.2))')
     .   (theta(k),k=1,kk)
c
cdiag do k=1,kk,3
cdiag write (text,'(''intf.pressure (m), k='',i3)') k+1
cdiag call prtmsk(ip,p(1,1,k+1),util1,idm,ii1,jj,0.,1./onem,text)
cdiag end do
c
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
css      tracer(i,j,k)=0.                 ! moved to hycom.f temperarily
      dp  (i,j,k+kk)=dp  (i,j,k)
      th3d(i,j,k+kk)=th3d(i,j,k)
      temp(i,j,k+kk)=temp(i,j,k)
      saln(i,j,k+kk)=saln(i,j,k)
      thstar(i,j,k)=th3d(i,j,k)+kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k))
c
      if (i.eq.itest.and.j.eq.jtest)
     . write (lp,'(2i5,i3,a,3f7.2,3x,2f7.3,f8.1)')
     .  i,j,k,'  dens,thstar,kappa=',th3d(i,j,k),thstar(i,j,k)
     .   ,kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k),2),
     .    temp(i,j,k),saln(i,j,k),p(i,j,k+1)/onem
c
      if (k.gt.1) then
        if (thstar(i,j,k).lt.thstar(i,j,k-1))
     .    totlj(j,k-1)=totlj(j,k-1)+1
      endif
 11   continue
c$OMP END PARALLEL DO
c
      do k=1,kk
      write (*,'(i5,a,2i5,a/7x,7(i3,3x),3x,7(i3,3x)/
     .  (/(7(i4,7f6.0,3x,7f6.0/))))') 
     .  k,' i,j=',itest,jtest,' input data (t,s,p,depth)'
     . ,        (j,j=jtest-3,jtest+3),(j,j=jtest-3,jtest+3)
     . ,(i,(temp(i,j,k),j=jtest-3,jtest+3)
     . ,   (saln(i,j,k),j=jtest-3,jtest+3),i=itest-3,itest+3)
     . ,(i,(p(i,j,k+1)/onem,j=jtest-3,jtest+3)
     . ,    (depths(i,j),j=jtest-3,jtest+3),i=itest-3,itest+3)
c
c    .write (*,'(2i4,a,i2/(5(5f7.1,3x,5f7.1/)))')
c    . itest,jtest,' input data (t,s,p,depth) at k',k,
c    . ((temp(i,j,k),j=jtest-2,jtest+2)
c    . ,(saln(i,j,k),j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((p(i,j,k+1)/onem,j=jtest-2,jtest+2)
c    . , (depth(i,j),j=jtest-2,jtest+2),i=itest-2,itest+2)
      enddo
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
cc$OMP PARALLEL DO
c     do 21 j=1,jj
c     do 21 l=1,isp(j)
c     do 21 i=ifp(j,l),ilp(j,l)
css   omlhc(i,j)=spcifh*p(i,j,2)/(onem *thref)           ! J/m*m C
c     thermb(i,j,k)=kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k))
c    .             *(1000.+th3d(i,j,k))
c     thermb(i,j,k+kk)=thermb(i,j,k)
c21   continue
cc$OMP END PARALLEL DO
      call ssto2a(temp,asst)
      call ssto2a(saln,sss)
c     call ssto2a(omlhc,mlhc)
c
      call findmx(ip,temp,ii,ii,jj,'ini sst')
      call findmx(ip,saln,ii,ii,jj,'ini sss')
c$OMP PARALLEL DO
      do 22 ja=1,jja
      do 22 ia=1,iia
      if (focean(ia,ja).gt.0.) then
        gtemp(1,1,ia,ja)=asst(ia,ja)
        if (sss(ia,ja).le.10.) then
          write(*,'(a,2i3,3(a,f6.1))')'chk low saln at agcm ',ia,ja
     . ,' sss=',sss(ia,ja),' sst=',asst(ia,ja),' focean=',focean(ia,ja)
          stop 'wrong sss in agcm'
        endif
      endif
 22   continue
c$OMP END PARALLEL DO
c
      else                                !  nstep0 > 0
c
c --- start from restart file prescribed
c
      write (lp,111) nstep0,time0
 111  format (9x,'chk time step in restart file -',i9,5x,' day ',f9.2)
c
      delt1=baclin+baclin
css   call newbot
c
c$OMP PARALLEL DO
      do i=1,ii
      do j=1,jj
      if (i.lt.ii) then
        scp2i(i,j)=1./scp2(i,j)
        scvyi(i,j)=1./scvy(i,j)
      end if
c
      if (i.gt.1) then
        scuxi(i,j)=1./scux(i,j)
        scq2i(i,j)=1./scq2(i,j)
      endif
      enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 21 j=1,jj
      do 21 k=1,kk
      do 21 l=1,isp(j)
      do 21 i=ifp(j,l),ilp(j,l)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
      thstar(i,j,k)=th3d(i,j,k)+kappaf(temp(i,j,k),saln(i,j,k),p(i,j,k))
 21   continue
c$OMP END PARALLEL DO
c
      end if                                !  nstep0 > 0  or  = 0
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
      print *,' focean'
      call zebra(focean,iia,iia,jja)
c
      print *,'chk ini. gtemp at nstep=',nstep0
      call zebra(asst,iia,iia,jja)
c
      print *,'chk ini. sss at nstep=',nstep0
      call zebra(sss,iia,iia,jja)
c
      if (itest.gt.0.and.jtest.gt.0) write (lp,103) nstep,itest,jtest,
     .  '  init.profile  temp    saln  thstar   thkns    dpth   montg',
     .  (k,temp(itest,jtest,k),saln(itest,jtest,k),
     .  thstar(itest,jtest,k),dp(itest,jtest,k)/onem,
     .  p(itest,jtest,k+1)/onem,montg(itest,jtest,k)/g,k=1,kk)
c

 103  format (i9,2i5,a/(28x,i3,2f8.2,f8.2,2f8.1,f8.3))
c
      return
      end
c
      subroutine refinp(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- add grid rows near equator to enhance merid. resolution
c --- eqcrs  -- i index of equator in unexpanded grid
c --- eqfin  -- i index of equator in expanded grid
c --- fieldc -- input field containing ii-2*(eqfin-eqcrs) rows
c --- fieldf -- output field containing ii rows
c
c --- use this entry to operate on -p- points
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
c
      integer numcrs,numfin,newrows,icrs,ifin,ieqcrs,ieqfin,
     .        ipcrs(idm,jdm)
      real fieldc(idm,jdm),fieldf(idm,jdm),
     .     eqcrs,eqfin,wgt,valu0,valu1,x,sumc,sumf
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=1)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0.,1./
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=7)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.25, 0.55, 0.95, 1.50, 2.20, 3.05, 4.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=9)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.400, 0.825, 1.300, 1.850, 2.500, 3.250, 4.100,
ccc     .               5.025, 6.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=11)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.60, 1.20, 1.80, 2.40, 3.02, 3.68, 4.40, 5.20,
ccc     .               6.08, 7.02, 8.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=12)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.375, 0.750, 1.125, 1.500, 1.900, 2.350, 2.875,
ccc     .               3.500, 4.250, 5.100, 6.025, 7.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=13)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.125, 0.250, 0.375, 0.525, 0.725, 1.000, 1.375,
ccc     .               1.875, 2.500, 3.250, 4.100, 5.025, 6.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=13)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.20, 0.40, 0.62, 0.88, 1.20, 1.60, 2.10, 2.70,
ccc     .               3.40, 4.20, 5.08, 6.02, 7.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=13)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.333333, 0.666667, 1.016667, 1.400000, 1.833333,
ccc     .               2.333333, 2.916667, 3.583333, 4.333333, 5.166667,
ccc     .               6.066667, 7.016667, 8.000000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      parameter (numfin=15)
      real coord(0:numfin),xnorm(0:numfin-1)
      data coord/0., 0.30, 0.60, 0.90, 1.20, 1.50, 1.82, 2.18, 2.60,
     .               3.10, 3.70, 4.40, 5.20, 6.08, 7.02, 8.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=16)
ccc      real coord(0:numfin),xnorm(0:numfin-1)
ccc      data coord/0., 0.28, 0.56, 0.84, 1.12, 1.40, 1.68, 1.96, 2.24,
ccc     .               2.52, 2.84, 3.24, 3.76, 4.40, 5.16, 6.04, 7.00/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      parameter (numfin=25)
ccc      real coord(0:numfin)
ccc      data coord/
ccc  . 0.000, 0.300, 0.600, 0.900, 1.200, 1.500, 1.800, 2.100, 2.400
ccc  .,2.700, 3.000, 3.301, 3.608, 3.928, 4.266, 4.630, 5.024, 5.456
ccc  .,5.931, 6.456, 7.037, 7.680, 8.392, 9.178,10.046,11.000/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(refinp)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(refinp)'
      end if
c
c$OMP PARALLEL DO PRIVATE(wgt,icrs,ifin,valu0,valu1,sumc,sumf)
      do 1 j=1,jj
c
c --- make room for 2*(numfin - numcrs) new grid rows
c
c --- rows north of equatorial strip:
      do 2 i=1,ieqcrs-numcrs
 2    fieldf(i,j)=fieldc(i,j)
c
c --- rows south of equatorial strip:
      do 3 i=ieqcrs+numcrs,ii-2*newrows
 3    fieldf(i+2*newrows,j)=fieldc(i,j)
c
c --- equatorial row:
      fieldf(ieqfin,j)=fieldc(ieqcrs,j)
c
c --- fill in remaining rows by interpolation
c
      do 4 i=1,numfin-1
      wgt=coord(i)-int(coord(i))
c
c --- south of equator:
      ifin=ieqfin+i
      icrs=ieqcrs+int(coord(i))
      valu0=fieldc(icrs  ,j)
      valu1=fieldc(icrs+1,j)
c --- make sure no land data are used in interpolation
      if (ipcrs(icrs,j).eq.0 .and. ipcrs(icrs+1,j).eq.1) valu0=valu1
      if (ipcrs(icrs,j).eq.1 .and. ipcrs(icrs+1,j).eq.0) valu1=valu0
      fieldf(ifin,j)=wgt*valu1+(1.-wgt)*valu0
cdiag if (j.eq.6) write (*,'(a,3i4,3f8.2,2i3)')
cdiag. 'i,ifin,icrs,wgt,valu1,valu0:',i,ifin,icrs,wgt,valu1,valu0
cdiag.  ,ipcrs(icrs+1,j),ipcrs(icrs,j)
c
c --- north of equator:
      ifin=ieqfin-i
      icrs=ieqcrs-int(coord(i))
      valu0=fieldc(icrs  ,j)
      valu1=fieldc(icrs-1,j)
c --- make sure no land data are used in interpolation
      if (ipcrs(icrs,j).eq.0 .and. ipcrs(icrs-1,j).eq.1) valu0=valu1
      if (ipcrs(icrs,j).eq.1 .and. ipcrs(icrs-1,j).eq.0) valu1=valu0
      fieldf(ifin,j)=wgt*valu1+(1.-wgt)*valu0
cdiag if (j.eq.6) write (*,'(a,3i4,3f8.2,2i3)')
cdiag. 'i,ifin,icrs,wgt,valu1,valu0:',i,ifin,icrs,wgt,valu1,valu0
cdiag.  ,ipcrs(icrs-1,j),ipcrs(icrs,j)
c
 4    continue
c
 1    continue
c$OMP END PARALLEL DO
      return
c
c
      entry refnap(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- use this entry if 'field' contains values proportional to cell size,
c --- i.e., values that need to be apportioned rather than interpolated
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(refnap)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(refnap)'
      end if
c
c$OMP PARALLEL DO PRIVATE(x)
      do 11 j=1,jj
c
c --- make room for 2*(numfin - numcrs) new grid rows
c
c --- rows north of equatorial strip:
      do 12 i=1,ieqcrs-numcrs
 12   fieldf(i,j)=fieldc(i,j)
c
c --- rows south of equatorial strip:
      do 13 i=ieqcrs+numcrs,ii-2*newrows
 13   fieldf(i+2*newrows,j)=fieldc(i,j)
c
c --- apportion remaining unrefined-grid rows among refined-grid rows
c
      do 15 i=ieqfin-numfin+1,ieqfin+numfin-1
 15   fieldf(i,j)=0.
c
      do 14 icrs=0,numcrs-1
      do 14 ifin=0,numfin-1
      if (icrs+ifin.eq.0) then
        x=.5*coord(1)
      else
        x=min(real(icrs)+.5,max(real(icrs)-.5,
     .     .5*(coord(ifin)+coord(ifin+1))))
     .  - min(real(icrs)+.5,max(real(icrs)-.5,
     .     .5*(coord(ifin)+coord(ifin-1))))
      end if
c
c --- south of equator:
      fieldf(ieqfin+ifin,j)=fieldf(ieqfin+ifin,j)
     .   +fieldc(ieqcrs+icrs,j)*ipcrs(ieqcrs+icrs,j)*x
c
c --- north of equator:
      fieldf(ieqfin-ifin,j)=fieldf(ieqfin-ifin,j)
     .   +fieldc(ieqcrs-icrs,j)*ipcrs(ieqcrs-icrs,j)*x
 14   continue
c
 11   continue
c$OMP END PARALLEL DO
      return
c
c
      entry refinu(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- use this entry to operate on -u- points
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(refinu)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(refinu)'
      end if
c
c$OMP PARALLEL DO PRIVATE(x,wgt,icrs,ifin)
      do 5 j=1,jj
c
c --- make room for 2*(numfin - numcrs) new grid rows
c
c --- rows north of equatorial strip:
      do 6 i=1,ieqcrs-numcrs
 6    fieldf(i,j)=fieldc(i,j)
c
c --- rows south of equatorial strip:
      do 7 i=ieqcrs+numcrs,ii-2*newrows
 7    fieldf(i+2*newrows,j)=fieldc(i,j)
c
c --- fill in remaining rows by interpolation
c
      do 8 i=1,numfin
      x=.5*(coord(i)+coord(i-1)+1.)
      wgt=x-int(x)
c
c --- south of equator:
      ifin=ieqfin+i
      icrs=ieqcrs+int(x)
      fieldf(ifin,j)=wgt*fieldc(icrs+1,j)+(1.-wgt)*fieldc(icrs,j)
cdiag if (j.eq.20) write (*,'(a,3i5,3f6.2)')
cdiag. 'i,ifin,icrs,wgt,fieldc(icrs+1),fieldc(icrs):',
cdiag.  i,ifin,icrs,wgt,fieldc(icrs+1,j),fieldc(icrs,j)
c
c --- north of equator:
      ifin=ieqfin-i+1
      icrs=ieqcrs-int(x)+1
      fieldf(ifin,j)=wgt*fieldc(icrs-1,j)+(1.-wgt)*fieldc(icrs,j)
cdiag if (j.eq.20) write (*,'(a,3i5,3f6.2)')
cdiag. 'i,ifin,icrs,wgt,fieldc(icrs-1),fieldc(icrs):',
cdiag.  i,ifin,icrs,wgt,fieldc(icrs-1,j),fieldc(icrs,j)
 8    continue
c
 5    continue
c$OMP END PARALLEL DO
      return
c
c
      entry unrefp(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- remove extra grid rows introduced by previous calls to 'refinp'
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(unrfin)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(a,2i4,a,i4,f6.2)') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(unrfin)'
      end if
c
c$OMP PARALLEL DO PRIVATE(wgt,icrs,ifin)
      do 21 j=1,jj
c
c --- copy extratropical rows into reduced-size array
c
c --- rows north of equatorial strip:
      do 22 i=1,ieqcrs-numcrs
 22   fieldc(i,j)=fieldf(i,j)
c
c --- rows south of equatorial strip:
      do 23 i=ieqcrs+numcrs,ii-2*newrows
 23   fieldc(i,j)=fieldf(i+2*newrows,j)
c
c --- equatorial row:
      fieldc(ieqcrs,j)=fieldf(ieqfin,j)
c
c --- fill in remaining rows by interpolation
c
      icrs=0
      do 24 i=2,numfin-1
c
      if (int(coord(i)).le.icrs) go to 24
      icrs=coord(i)
      wgt=(coord(i)-icrs)/(coord(i)-coord(i-1))
c
c --- south of equator:
      ifin=ieqfin+i
      fieldc(ieqcrs+icrs,j)=wgt*fieldf(ifin-1,j)+(1.-wgt)*fieldf(ifin,j)
c
c --- north of equator:
      ifin=ieqfin-i
      fieldc(ieqcrs-icrs,j)=wgt*fieldf(ifin+1,j)+(1.-wgt)*fieldf(ifin,j)
 24   continue
      if (icrs.ne.int(coord(numfin-1))) then
        write (lp,'(2(a,i3))') 'icrs=',icrs,'  not',int(coord(numfin-1))
        stop '(unrefp)'
      end if
c
 21   continue
c$OMP END PARALLEL DO
      return
c
      entry unrefu(eqcrs,eqfin,ipcrs,fieldc,fieldf)
c
c --- use this entry to remove extra rows of -u- points
c
      numcrs=coord(numfin)
      if (real(numcrs).ne.coord(numfin)) then
        write (lp,*) '-coord- array must end on whole number'
        stop '(unrefu)'
      end if
      newrows=numfin-numcrs
      ieqcrs=eqcrs
      ieqfin=eqfin
      if (ieqfin-ieqcrs .ne. newrows) then
        write (lp,'(2(a,2i4))') 'ieqfin/old =',ieqfin,ieqcrs,
     .  '  inconsistent with numfin/old =',numfin,numcrs
        stop '(unrefu)'
      end if
c
c$OMP PARALLEL DO PRIVATE(x,wgt,icrs,ifin)
      do 25 j=1,jj
c
c --- copy extratropical rows into reduced-size array
c
c --- rows north of equatorial strip:
      do 26 i=1,ieqcrs-numcrs
 26   fieldc(i,j)=fieldf(i,j)
c
c --- rows south of equatorial strip:
      do 27 i=ieqcrs+numcrs,ii-2*newrows
 27   fieldc(i,j)=fieldf(i+2*newrows,j)
c
c --- fill in remaining rows by interpolation
c
      icrs=0
      do 28 i=2,numfin
      x=.5*(coord(i)+coord(i-1)+1.)
      if (int(x).le.icrs) go to 28
      icrs=x
      wgt=(x-icrs)/(.5*(coord(i)-coord(i-2)))
c
c --- south of equator:
      ifin=ieqfin+i
      fieldc(ieqcrs+icrs  ,j)=
     .   wgt*fieldf(ifin-1,j)+(1.-wgt)*fieldf(ifin,j)
c
c --- north of equator:
      ifin=ieqfin-i+1
      fieldc(ieqcrs-icrs+1,j)=
     .   wgt*fieldf(ifin+1,j)+(1.-wgt)*fieldf(ifin,j)
 28   continue
      if (icrs.ne.int(coord(numfin))) then
        write (lp,'(2(a,i3))') 'icrs=',icrs,'  not',int(coord(numfin))
       stop '(unrefu)'
      end if
c
 25   continue
c$OMP END PARALLEL DO
      return
      end
c> Revision history:
c>
c> Mar. 2000 - conversion to SI units
c> Aug. 2000 - added diagnostic count of static instabilities
c> Apr. 2001 - eliminated stmt_funcs.h
c> Sep. 2005 - added EQ refinement 
