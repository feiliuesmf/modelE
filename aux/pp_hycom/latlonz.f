      program latlonz3d
c  (1) read in monthly output from hycom
c  (2) convert fields to lat/lon/z grid of 1x1x33.
c
      use hycom_arrays, only: depths,srfhgt,dpmixl,oice,u,v,dp,p
     .   ,temp,saln,th3d,tracer,alloc_hycom_arrays
      use hycom_dimen
      use const_proc
      use hycom_o2a
      implicit none
c
      real :: year1
      integer, allocatable :: im(:,:)
c
      integer mo,ny1,ny2,dcd,mon1,i70,i45,ieq,status,ia,k00,ny,m,mm,nt
      integer :: nrec
      real*8 :: area,avgbot
      character flnmin*80,runid*20,ttl*80,ttl1*80,ttl2*80
     .         ,flnmann*80,flnmout*80
      logical timav,cnvert
      character*26 ttlt(ntrcr)
      character*26,dimension(8),parameter::
     .  ttl0=(/'sea surface height (cm)   '
     .          ,'mixed layer depth (m)     '
     .          ,'sea ice coverage [0:1]    '
     .          ,'eastward  velocity (cm/s) '
     .          ,'northward velocity (cm/s) '
     .          ,'temperature (deg C)       '
     .          ,'salinity                  '
     .          ,'density (sigma2)          '/)
      real :: avg2d(iia,jja),avg3d(iia,jja,k33)
     .  ,sshij(iia,jja),dpmixij(iia,jja),iceij(iia,jja)
     .  ,tz(iia,jja,k33),sz(iia,jja,k33),uz(iia,jja,k33)
     .  ,vz(iia,jja,k33),rz(iia,jja,k33),trz(iia,jja,k33,ntrcr)
     .  ,a2d(iia,jja,12),a3d(iia,jja,k33,12)
     .  ,worka(iia,jja),worko(idm,jdm),depthij(iia,jja)
      real*4 :: real4a(iia,jja)
c
      read(*,'(a/i4/i4)') runid,ny1,ny2
      write(*,'(3a,i4,a,i4)') 
     .   'processing ',runid,' from yr ',ny1,' to ',ny2
      write(*,'(a,i2)') 'number of tracers =',ntrcr
      nrec=8+ntrcr

      do nt=1,ntrcr
      write(ttlt(nt),'(a,i2.2)') 'tracer No.',nt
      enddo

      call alloc_hycom_arrays
      call alloc_hycom_dimen

      allocate (im(idm,jdm),stat=status)
c --- determine do-loop limits for u,v,p,q points
      call gtdpth(depths,im)
      call bigrid(depths)
      if (diag) then
        print *,' depth at i,j=',itest,jtest
        do i=itest-5,itest+5
        write(*,'(11f6.0)') (depths(i,j),j=jtest-5,jtest+5)
        enddo
      endif
c
      dcd=ny1/10
      if (runid(1:1).eq.' ') stop 'empty runid'
      if (dcd.lt.180 .or. dcd.gt.330) then
        print *,' wrong decade=',dcd
        stop 'wrong decade'
      endif
c
      call o2a_wgt
      call o2a_sfc(depths,depthij)
      if (diag) then
        i=iatest
        j=jatest
        print *,'topo after o2a',i,j,depths(itest,jtest),depthij(i,j)
      endif
c
      kij=0
      do 30 i=1,iia
      do 30 j=1,jja
      if (depthij(i,j).le.0.) goto 32
        do 31 k=2,k33
        if (z33(k-1).le.depthij(i,j).and.z33(k).gt.depthij(i,j)) then
          kij(i,j)=k-1
          go to 32
        elseif (z33(k33).le.depthij(i,j)) then
          kij(i,j)=k33
          go to 32
        endif
 31     continue
 32     continue
 30   continue
      if (diag) then
      write(*,'(21f5.0)') ((depthij(i,j),i=iatest-10,iatest+10)
     .                                  ,j=jatest-10,jatest+10)
      print *,' kij main'
      write(*,'(21i5)') ((kij(i,j),i=iatest-10,iatest+10)
     .                            ,j=jatest-10,jatest+10)
      i=iatest
      j=jatest
      write(*,'(a,3i4,f6.0)') 'chk kij=',i,j,kij(i,j),depthij(i,j)
      endif
c
c     real4a=depthij
c     write(ttl,'(a)') 
c    .'Ocean depth(m)  converted from hycom1deg to regular grid of
c    . 360x180'
c     open(41,file='depth_360x180.tbin',form='unformatted',
c    .   status='unknown')
c     write(41) ttl,real4a
c     close(41)

      do 151 ny=ny1,ny2
      do 152 mo=mo1,mo2
c     write(flnmin,'(5a,i3,a,i3,2a,i4.4,2a)')trim(path1),trim(runid)
c    . ,'/out',trim(runid),'_',dcd,'0_',dcd,'9/',amon(mo),ny
c    . ,'.out',trim(runid)
      write(flnmin,'(2a,i4.4,2a)')trim(path1),amon(mo),ny
     .  ,'.out',trim(runid)
      write(flnmout,'(2a,i4.4,2a)')trim(path1),amon(mo),ny
     .  ,'.zout',trim(runid)
      write(ttl1,'(i3,a,i3,3x,2(1x,a),i4)') 
     .            iia,'x',jja,trim(runid),amon(mo),ny
      write(ttl2,'(2(i3,a),i2,2(1x,a),i4)') iia,'x',jja,'x',k33
     .         ,trim(runid),amon(mo),ny
      open(80+mo,file=path2//flnmout,form='unformatted'
     .                              ,status='unknown')
      write (lp,'(2a)') 'reading input from ',trim(flnmin)
      write (lp,'(2a)') 'output file:       ',trim(flnmout)

      mon1=monlg(mod(mo-1,12)+1)
c
c --- read archive data
      timav=.true.
      cnvert=.false.
      call getdat(flnmin,year1,timav,cnvert)
c
      call o2a_sfc(srfhgt,sshij)
      if (diag) then
        i=iatest
        j=jatest
        print *,'ssh after o2a ',i,j,srfhgt(itest,jtest),sshij(i,j)
      endif

      call o2a_sfc(dpmixl,dpmixij)
      if (diag) then
        i=iatest
        j=jatest
        print *,'dpm after o2a ',i,j,dpmixl(itest,jtest,1),dpmixij(i,j)
c       call prtmsk(kij,dpmixij,worka,iia,iia,jja,0.,1.,'dpmxl_ij')
      endif

      call o2a_sfc(oice,iceij)
      call o2a_3dvec(p,u,v,uz,vz)
      call o2a_3d(p,temp,tz)
      call o2a_3d(p,saln,sz)
      call o2a_3d(p,th3d,rz)
      do nt=1,ntrcr
      call o2a_3d(p,tracer(1,1,1,nt),trz(1,1,1,nt))
      enddo

      if (diag) then
c       call prtmsk(kij,temp,worka,iia,iia,jja,0.,1.,'sst')
        i=iatest
        j=jatest
        print *,'chk     t   s   u  v  trc  at ',i,j,depthij(i,j)
        do k=1,k33
        write(*,'(i2,f7.0,5f12.4)')k,z33(k),tz(i,j,k),sz(i,j,k)
     .     ,uz(i,j,k),vz(i,j,k),trz(i,j,k,1)
        enddo

        i=iatest
        j=jatest
        print *,'chk     temp  ',i,j
        do k=1,k33
        write(*,'(i2,5f12.4)') k,tz(i-1,i,k),tz(i,j,k),tz(i+1,j,k)
     .        ,tz(i,j-1,k),tz(i,j+1,k)
        enddo

        do k=1,k33
        print *,'chk  final  temp ', iatest,jatest,k
        do j=jatest+3,jatest-3,-1
        write(*,'(11f7.1)') (tz(i,j,k),i=iatest-3,iatest+3)
        enddo
        enddo
      endif

      write(ttl,'(a,6x,a)')  ttl0(1),trim(ttl1)
      write(80+mo) ttl, sshij
      write(ttl,'(a,6x,a)')  ttl0(2),trim(ttl1)
      write(80+mo) ttl, dpmixij
      write(ttl,'(a,6x,a)')  ttl0(3),trim(ttl1)
      write(80+mo) ttl, iceij
c
      write(ttl,'(a,6x,a)') ttl0(4),trim(ttl2)
      write(80+mo) ttl, uz
      write(ttl,'(a,6x,a)') ttl0(5),trim(ttl2)
      write(80+mo) ttl, vz
      write(ttl,'(a,6x,a)') ttl0(6),trim(ttl2)
      write(80+mo) ttl, tz
      write(ttl,'(a,6x,a)') ttl0(7),trim(ttl2)
      write(80+mo) ttl, sz
      write(ttl,'(a,6x,a)') ttl0(8),trim(ttl2)
      write(80+mo) ttl, rz
      do nt=1,ntrcr
      write(ttl,'(a,6x,a)') ttlt(nt),trim(ttl2)
      write(80+mo) ttl,trz(:,:,:,nt)
      enddo
 152  continue   ! mo=mo1,mo2
c
      if (mo1.ne.1 .or. mo2.ne.12) goto 151     ! no annual mean
      mo=13
      write(ttl1,'(i3,a,i3,3x,2(1x,a),i4)') 
     .            iia,'x',jja,trim(runid),amon(mo),ny
      write(ttl2,'(2(i3,a),i2,2(1x,a),i4)') iia,'x',jja,'x',k33
     .         ,trim(runid),amon(mo),ny
      write(flnmann,'(a,i4.4,2a)') amon(mo),ny,'.zout',trim(runid)
      open(93,file=path2//flnmann,form='unformatted',status='unknown')
      write(*,*) 'output annual =',trim(flnmann)
c
      do 10 n=1,nrec
      do 92 mo=1,12
      if (n==1) rewind(80+mo)
      if (n.le.3) then
        read(80+mo) ttl,((a2d(i,j,mo),i=1,iia),j=1,jja)
      else
        read(80+mo) ttl,(((a3d(i,j,k,mo),i=1,iia),j=1,jja),k=1,k33)
      endif
      if (n.eq.3.and.a3d(iatest,jatest,1,mo).le.flag/9.) then
        print *,' sst=',m,a3d(iatest,jatest,1,m)
        stop 'wrong sst'
      endif
 92   continue
c
      if (n.le.3) then
        avg2d=0.
        do 43 j=1,jja
        do 43 i=1,iia
        if (kij(i,j).ge.1) then
          do m=1,12
            avg2d(i,j)=avg2d(i,j)+a2d(i,j,m)/12.
            if (n.eq.1.and.i.eq.iatest.and.j.eq.jatest.and.diag) 
     .  write(*,'(a,i2,f7.2)') 'chk dpl=',m,avg2d(iatest,jatest)
          enddo
        else
          avg2d(i,j)=flag
        endif
 43     continue
        if (diag) write(*,*) 'got     : ',ttl(1:25)
        write(ttl,'(a,6x,a)')  ttl0(n),trim(ttl1)
        if (diag) write(*,*) 'expected: ',ttl(1:25)
        write(93) ttl,avg2d
      else
        avg3d=0.
        do 44 k=1,k33
        do 44 j=1,jja
        do 44 i=1,iia
        if (kij(i,j).ge.1) then
          do m=1,12
            if (a3d(i,j,k,m).le.flag .or. avg3d(i,j,k).le.flag) then
              avg3d(i,j,k)=flag
            else
              avg3d(i,j,k)=avg3d(i,j,k)+a3d(i,j,k,m)/12.
            endif
          enddo
        else
          avg3d(i,j,k)=flag
        endif
 44     continue
        if (n.le.8) then
          write(ttl,'(a,6x,a)') ttl0(n),trim(ttl2)
        else
          write(ttl,'(a,6x,a)') ttlt(n-8),trim(ttl2)
        endif
        write(93) ttl,avg3d
        if (diag) then
        write(*,*) trim(ttl)
        do k=1,k33
        write(*,'(a,i2,f7.1 )')'chk k=',k,avg3d(iatest,jatest,k)
        enddo
        endif
      endif      ! n<=3

  10  continue   ! n=1,nrec
 151  continue   ! ny=ny1,ny2

      do mo=1,13
      close(80+mo)
      enddo

      stop '(normal finish converting)'
      end

      real function sphdis(x1,y1,x2,y2)
c --- dist.(m) between 2 points on sphere, lat/lon (x1,y1) and lat/lon (x2,y2)
      implicit none
      real x1,y1,x2,y2,ang,radius,radian
      data radius/6376.e3/,radian/57.2957795/
c     
      ang=mod(y2-y1+540.,360.)-180.
      sphdis=radius*acos(min(1.,cosd(90.-x1)*cosd(90.-x2)
     .                         +sind(90.-x1)*sind(90.-x2)*cosd(ang)))
      if (sphdis.eq.0.)        
     .  sphdis=radius*sqrt((x2-x1)**2+(ang*cosd(.5*(x1+x2)))**2)/radian
c     if (sphdis.eq.0.) write (*,'(a,2f8.3,2x,2f8.3)')
c    .  'warning - zero distance between lat/lon points',x1,y1,x2,y2
      sphdis=max(sphdis,1.)
      return
      end
