      subroutine getdat(flnm,year1,timav,cnvert)
c
c --- read hybrid model fields (binary) and extracert portion of global fields.
c --- then  t r a n s f o r m   to   i s o p y c n i c   fields
c
      use hycom_arrays, only: ubavg,vbavg,srfhgt,dpmixl,oice,depths,
     .u,v,dp,p,temp,saln,th3d,tracer,uflx,vflx,diaflx,alloc_hycom_arrays
      use hycom_dimen
      use const_proc
      implicit none
c
      character(len=*), intent (in) :: flnm
      logical,   intent (in) :: cnvert,timav

      real*4 real4(idm,jdm),thbase4,time4,theta4(kdm)
      real :: pnew(idm,jdm,kdm+1),theta(kdm),year1
      integer*4 :: lenrec,isiz,jsiz,ksiz,nstep,lgth
      integer :: nrec,lev,imax(kdm),jmax(kdm),nt,iyr,imon,iday
      character what*16,trcrid*12(ntrcr),info*17
      integer, parameter :: ni=11
c
      lenrec=100
      if (diag) write (lp,*) ' open file with record length',lenrec
      open (unit=ni,file=flnm,status='old',access='direct',
     .   recl=100,form='unformatted',err=6)
      read (ni,rec=1,err=6) lgth
      close (unit=ni)
c
c --- f90-compiled scmicom reports file length in bytes
      if (diag) write (lp,*) '  reopen file with record length',lgth
     .  ,' timav=',timav,' cnvert=',cnvert
c
      open (unit=ni,file=flnm,status='old',access='direct',
     .   recl=lgth,form='unformatted',err=6)
c
      nrec=1
      read (ni,rec=nrec) lenrec,isiz,jsiz,ksiz,nstep,time4,
     .     thbase4,theta4,info
      if (diag) 
     . write (lp,*) lenrec,isiz,jsiz,ksiz,nstep,time4,thbase4,theta4
     .   ,info
      if (info(2:11).eq.'          ') then
        year1=time4/365.           !ocean model age, begin at 0 year
      else
        read(info(8:11),'(i4)') iyr
        read(info(5: 6),'(i2)') iday
        read(info(2: 3),'(i2)') imon
        if (imon.le.0 .or. imon.gt.12) then
          print *,' wrong month=',imon
          stop 'wrong month'
        elseif (iday.lt.0 .or. iday.gt.31) then
          print *,' wrong day=',iday
          stop 'wrong day'
        endif
        year1=iyr+jdofm(imon)/365. !coupled model clock, begin at, e.g. yr 1800
      endif    !use model year or model age

      do 35 k=1,kdm
 35   theta(k)=theta4(k)+thbase4
c
      if (diag) write (lp,'(''nstep, year (age/agcm):''
     .  ,i7,f7.2,f10.2,2x,a/'' sigma values:''/
     .   (1x,11f7.2))') nstep,time4/365.,year1,info,(theta(k),k=1,kdm)
c
      nrec=1
      if (timav) nrec=6+(9+ntrcr)*kdm		!  skip instantaneous fields

      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: ubavg        found: ',what,lev
 100  format ('record',i5,5x,2a,i3)
      call extrct(real4,ubavg(1,1,1))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: vbavg        found: ',what,lev
      call extrct(real4,vbavg(1,1,1))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: srfhgt       found: ',what,lev
      call extrct(real4,srfhgt(1,1))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: dpmixl       found: ',what,lev
      call extrct(real4,dpmixl(1,1,1))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: icecover     found: ',what,lev
      call extrct(real4,oice(1,1))
c     call prt9x9(oice,itest,jtest,0.,1.,'ice thknss')
c
 1    do 14 k=1,kdm
c
      if (timav) nrec=1+5*2+(9+ntrcr)*kdm+6*(k-1)   ! reading 3d avg field
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: u            found: ',what,lev
      call extrct(real4,u   (1,1,k))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: v            found: ',what,lev
      call extrct(real4,v   (1,1,k))
    
c     if (k.eq.1)
c    .call prtmsk(iv,v(1,1,1),pnew,idm,idm,jdm,0.,100.,'v(cm/s)')
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: dp           found: ',what,lev
      call extrct(real4,dp  (1,1,k))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: temp         found: ',what,lev
      call extrct(real4,temp(1,1,k))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: saln         found: ',what,lev
      call extrct(real4,saln(1,1,k))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: th3d         found: ',what,lev
      call extrct(real4,th3d(1,1,k))
c     if (k.eq.1) call findmx(im,th3d,idm,idm,jdm,'surf.density')
c
      if (timav) nrec=6+6*k+ntrcr*(k-1)	! reading tracer field
c
      do nt=1,ntrcr			! ntrcr tracers in input archive
        nrec=nrec+1
        read (ni,rec=nrec) what,lev,real4
        if (diag) 
     .  write (lp,100) nrec,'expected: tracer       found: ',what,lev
        call extrct(real4,tracer(1,1,k,nt))
      enddo

 14   continue
c
      if (timav) nrec=6+(6+ntrcr)*kdm
      do 15 k=1,kdm
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: uflxav       found: ',what,lev
      call extrct(real4,uflx(1,1,k))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: vflxav       found: ',what,lev
      call extrct(real4,vflx(1,1,k))
c
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
      if (diag) 
     .  write (lp,100) nrec,'expected: diaflx       found: ',what,lev
      call extrct(real4,diaflx(1,1,k))
c
ccc      write (lp,*) ' shown below: diaflx, layer',k
c
 15   continue
c
      nrec=1+5*2+(6*2+3+ntrcr)*kdm		        !  skip to surface fluxes
      nrec=nrec+1
      read (ni,rec=nrec) what,lev,real4
c     write (lp,100) nrec,'expected: eminp        found: ',what,lev
c     call extrct(real4,eminp )
c
      nrec=nrec+1
c     read (ni,rec=nrec) what,lev,real4
c     write (lp,100) nrec,'expected: htflx        found: ',what,lev
c     call extrct(real4,htflx)
c
      nrec=nrec+1
c     read (ni,rec=nrec) what,lev,real4
c     write (lp,100) nrec,'expected: sflx         found: ',what,lev
c     call extrct(real4,sflx)
c
      nrec=nrec+1
c     read (ni,rec=nrec) what,lev,real4
c     write (lp,100) nrec,'expected: brine        found: ',what,lev
c     call extrct(real4,brine)
c
      close (unit=ni)
c
c --- fix physical dimensions
c
c$OMP PARALLEL DO
      do 20 j=1,jdm
c
      do 21 i=1,idm
      ubavg(i,j,1)=ubavg(i,j,1)*100./dz * iu(i,j)	!  convert to cm/sec
      vbavg(i,j,1)=vbavg(i,j,1)*100./dz * iv(i,j)	!  convert to cm/sec
      if (ip(i,j).gt.0) then
        dpmixl(i,j,1)=dpmixl(i,j,1)/(g*rho*dz)	    ! convert to meters
        srfhgt(i,j)=srfhgt(i,j)*100.        	    ! convert to cm
c       eminp(i,j)=eminp(i,j)*365.*86400.	    ! convert m/s to m/ir
c        sflx(i,j)= sflx(i,j)*1.e-3/35.*365.*86400. ! convert g/m2/s to m/yr
c       brine(i,j)=brine(i,j)*1.e-3/35.*365.*86400. ! convert g/m2/s to m/yr
      else
        dpmixl(i,j,1)=flag
        srfhgt(i,j)=flag
        oice(i,j)=flag
      end if
 21   continue
c
      do 20 i=1,idm
      p(i,j,1)=0.
      do 20 k=1,kdm
      u(i,j,k)=u(i,j,k)*100./dz * iu(i,j)+ubavg(i,j,1)	!  convert to cm/sec
      v(i,j,k)=v(i,j,k)*100./dz * iv(i,j)+vbavg(i,j,1)	!  convert to cm/sec
      uflx(i,j,k)=uflx(i,j,k)*1.e-6/(g*rho*dz**3)       !  convert to Sv*intvl
      vflx(i,j,k)=vflx(i,j,k)*1.e-6/(g*rho*dz**3)       !  convert to Sv*intvl

      if (ip(i,j).gt.0) then
        dp(i,j,k)=dp(i,j,k)/(g*rho*dz)			!  convert to meters
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
        th3d(i,j,k)=(th3d(i,j,k)+thbase4)*1000./rho	!  add thbase
        if (k.ge.2) th3d(i,j,k)=max(th3d(i,j,k-1),th3d(i,j,k))
      else
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
        temp(i,j,k)=flag
        saln(i,j,k)=flag
        th3d(i,j,k)=flag
      end if
 20   continue
c$OMP END PARALLEL DO
c
      if (itest.gt.0.and.jtest.gt.0.and.diag)write (lp,103) itest,jtest,
     .  '  init.profile  temp    saln  thkns    dpth   tracer',
     .  (k,temp(itest,jtest,k),saln(itest,jtest,k),
     .  dp(itest,jtest,k),p(itest,jtest,k+1),tracer(itest,jtest,k,1),
     .  k=1,kdm)
 103  format (2i5,a/(18x,i3,5f8.2))
      if (diag) then
        print *,' itest,jtest=',itest,jtest
        print *,' dpmix=',dpmixl(itest,jtest,1)
        print *,' p=',p(itest,jtest,2)
        print *,' srfhgt=',srfhgt(itest,jtest)
        print *,' depths=',depths(itest,jtest)
c       call prtmsk(ip,dpmixl,real4,idm,idm,jdm,0.,1.,'mld (m)')
c
        print *,' chk input'
        do k=1,kdm
        write(*,'(i2,4(2i4,f7.1,f6.2,2x))')
     .    k,(((i,j,p(i,j,k+1),temp(i,j,k)),
     .    i=itest-1,itest),j=jtest-1,jtest)
        enddo
      endif

      return
 6    stop 'error in reading'
      return
      end
