      program latlonz3d
c  (1) read in monthly output from hycom
c  (2) convert fields to lat/lon/z grid of 1x1x33.
c
      use const_proc
      use hycom_arrays, only: depths,srfhgt,dpmixl,oice,u,v,dp,p
     .   ,temp,saln,th3d,tracer,alloc_hycom_arrays
      use hycom_dimen
      use hycom_o2a

      implicit none
c
      real :: year1
      integer, allocatable :: im(:,:)
c
      integer mo,mon1,i70,i45,ieq,status,ia,k00,ny,m,mm,nt
      integer :: nrec
      real*8 :: area,avgbot
      character flnmin*128,ttl*80,ttl1*80,ttl2*80
     .         ,flnmann*80,flnmout*128
      logical timav,cnvert
      character*26 ttlt(ntrcr)
      character*26,dimension(8),parameter::
     .  ttl0=(/'sea surface height (cm)   '
     .          ,'mixed layer depth (m)     '
     .          ,'sea ice coverage [0:1]    '
     .          ,'eastward  velocity (cm/s) '
     .          ,'northward velocity (cm/s) '
     .          ,'temperature (deg C)       '
     .          ,'salinity (psu)            '
     .          ,'density (sigma2)          '/)
      real :: avg2d(iia,jja),avg3d(iia,jja,k33)
     .  ,sshij(iia,jja),dpmixij(iia,jja),iceij(iia,jja)
     .  ,tz(iia,jja,k33),sz(iia,jja,k33),uz(iia,jja,k33)
     .  ,vz(iia,jja,k33),rz(iia,jja,k33),trz(iia,jja,k33,ntrcr)
     .  ,a2d(iia,jja,12),a3d(iia,jja,k33,12)
     .  ,worka(iia,jja),worko(idm,jdm),depthij(iia,jja)
     .  ,srfhgt_av(idm,jdm),dpmixl_av(idm,jdm)
     .  ,oice_av(idm,jdm),         p_av(idm,jdm,kdm+1)
     .  ,   u_av(idm,jdm,kdm),     v_av(idm,jdm,kdm)
     .  ,temp_av(idm,jdm,kdm),  saln_av(idm,jdm,kdm)
     .  ,th3d_av(idm,jdm,kdm),tracer_av(idm,jdm,kdm,ntrcr)
      real*4 :: real4a(iia,jja)
c
      character(len=1024) :: cmnd_arg   
      character(len=  80) :: cmnd_title 
      integer :: num_args, cmnd_len, cmnd_status, iarg
      real :: factor

      namelist /hdiag_nml/ path0, path1, path2,
     . hycomtopo, latlonij, basinmask, flnmcoso, flnmo2a,
     . runid, ny1, ny2, monave_convert,solo_convert


      call get_command_argument (0, cmnd_arg, cmnd_len, cmnd_status)
      if (cmnd_status /= 0) then
         write (*,*) 'Getting command name failed with cmnd_status = ',
     .                                                 cmnd_status
         stop
      end if
      write (unit=*,fmt='(/a,3x,a)') 'command name =', 
     .      cmnd_arg (1:cmnd_len)

      num_args = command_argument_count ()
      write (unit=*,fmt='(a,i3)') 'number of command arguments = ', 
     .      num_args

       iarg=1
       call get_command_argument (iarg, cmnd_arg, cmnd_len, cmnd_status)
       if (cmnd_status /= 0) then
           write (unit=*,fmt='(a,i2,a)') 
     .          'get_command_argument failed: cmnd_status = ', 
     .          cmnd_status,  ' arg = ', iarg
        stop
       end if
       flnmout=trim(cmnd_arg)
       write (unit=*,fmt='(/a,i2,2a)') 'command arg=', iarg, 
     .       ' (name output file) = ', trim(flnmout)
       open(40,file=trim(flnmout),form='unformatted',status='unknown')
       write (unit=*,fmt='(a)') 'output file is open for writing '

       iarg=2
       call get_command_argument (iarg, cmnd_arg, cmnd_len, cmnd_status)
       if (cmnd_status /= 0) then
          write (*,*) 'get_command_argument failed: cmnd_status = ', 
     .                 cmnd_status,  ' arg = ', iarg
        stop
       end if
       cmnd_title=trim(cmnd_arg)
       write (unit=*,fmt='(/a,i2,3a)') 'command arg=', iarg, 
     .   ' (title output file) ="',trim(cmnd_title),'"'
c
      write(*,'(/a,i2)') 'number of tracers =',ntrcr
      nrec=8+ntrcr

      do nt=1,ntrcr
        write(ttlt(nt),'(a,i2.2)') 'tracer No.',nt
      enddo

      write(ttl1,'(i3,a,i3,3x)')  iia,'x',jja
      write(ttl2,'(2(i3,a),i2)') iia,'x',jja,'x',k33
c  Check the length of output title. Should be lesss 80
c  write(ttl,'(a,3x,a,3x,a)') ttl0(1),trim(ttl1),trim(cmnd_title)
      if( len(ttl0(1))+3+len(trim(ttl2))+3+ 
     .    len(trim(cmnd_title)) > 80 ) then
        write(*,*) "The length of the title at the output file "
        write(*,*) "    exceed 80 characters"
        write(*,*) "You should decrease length title from comand line!"
        write(*,*) "Safe length less than 38 caracters !!!"
C       write(*,*) " STOP at file=",__FILE__," line=", __LINE__  
      end if

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
      end if

      srfhgt_av=0.
      dpmixl_av=0.
        oice_av=0.
           p_av=0.
           u_av=0.
           v_av=0.
        temp_av=0.
        saln_av=0.
        th3d_av=0.
      tracer_av=0.
c     
      runid="Eh100"
c
      factor = 1./ (num_args-2)
      do iarg = 3, num_args
      call get_command_argument (iarg, cmnd_arg, cmnd_len, cmnd_status)
      if (cmnd_status /= 0) then
          write (*,*) 'get_command_argument failed: cmnd_status = ', 
     .                 cmnd_status,  ' arg = ', iarg
        stop
      end if
      flnmin=trim(cmnd_arg)
      write (unit=*,fmt='(/a,i3.3,2a)') 'command arg=', iarg, 
     .       ' (name input(reading)  file) = ', trim(flnmin)

c --- read archive data
      timav=.true.
      cnvert=.false.
      call getdat(flnmin,year1,timav,cnvert)
c
      do  i=1,idm
      do  j=1,jdm
      srfhgt_av(i,j)=srfhgt_av(i,j)+srfhgt(i,j)  *factor
      dpmixl_av(i,j)=dpmixl_av(i,j)+dpmixl(i,j,1)*factor
        oice_av(i,j)=  oice_av(i,j)+  oice(i,j)  *factor
      do k=1,kdm
       p_av(i,j,k+1)= p_av(i,j,k+1)+ p(i,j,k+1)*factor
         u_av(i,j,k)=   u_av(i,j,k)+   u(i,j,k)*factor
         v_av(i,j,k)=   v_av(i,j,k)+   v(i,j,k)*factor
      temp_av(i,j,k)=temp_av(i,j,k)+temp(i,j,k)*factor
      saln_av(i,j,k)=saln_av(i,j,k)+saln(i,j,k)*factor
      th3d_av(i,j,k)=th3d_av(i,j,k)+th3d(i,j,k)*factor

      do nt=1,ntrcr
      tracer_av(i,j,k,nt)=tracer_av(i,j,k,nt)+tracer(i,j,k,nt)*factor
      enddo ! loop nt

      enddo ! loop k

      end do   ! loop j 
      end do   ! loop i
      end do   ! loop arguments 3, .
c
      call o2a_sfc(srfhgt_av,sshij)
      if (diag) then
        i=iatest
        j=jatest
        print *,'ssh after o2a ',i,j,srfhgt_av(itest,jtest),sshij(i,j)
      endif

      call o2a_sfc(dpmixl_av,dpmixij)
      if (diag) then
        i=iatest
        j=jatest
        print *,'ml aft o2a ',i,j,dpmixl_av(itest,jtest),dpmixij(i,j)
c       call prtmsk(kij,dpmixij,worka,iia,iia,jja,0.,1.,'dpmxl_ij')
      endif

      call o2a_sfc(oice_av,iceij)
      call o2a_3dvec(p_av,u_av,v_av,uz,vz)
      call o2a_3d(p_av,temp_av,tz)
      call o2a_3d(p_av,saln_av,sz)
      call o2a_3d(p_av,th3d_av,rz)
      do nt=1,ntrcr
      call o2a_3d(p_av,tracer_av(1,1,1,nt),trz(1,1,1,nt))
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

      write(ttl,'(a,3x,a,3x,a)') ttl0(1),trim(ttl1),trim(cmnd_title)
      write(40) ttl, sshij
      write(ttl,'(a,3x,a,3x,a)') ttl0(2),trim(ttl1),trim(cmnd_title)
      write(40) ttl, dpmixij
      write(ttl,'(a,3x,a,3x,a)') ttl0(3),trim(ttl1),trim(cmnd_title)
      write(40) ttl, iceij
c
      write(ttl,'(a,3x,a,3x,a)') ttl0(4),trim(ttl2),trim(cmnd_title)
      write(40) ttl, uz
      write(ttl,'(a,3x,a,3x,a)') ttl0(5),trim(ttl2),trim(cmnd_title)
      write(40) ttl, vz
      write(ttl,'(a,3x,a,3x,a)') ttl0(6),trim(ttl2),trim(cmnd_title)
      write(40) ttl, tz
      write(ttl,'(a,3x,a,3x,a)') ttl0(7),trim(ttl2),trim(cmnd_title)
      write(40) ttl, sz
      write(ttl,'(a,3x,a,3x,a)') ttl0(8),trim(ttl2),trim(cmnd_title)
      write(40) ttl, rz

      do nt=1,ntrcr
      write(ttl,'(a,3x,a,3x,a)') ttlt(nt),trim(ttl2),trim(cmnd_title)
      write(40) ttl,trz(:,:,:,nt)
      enddo

      close(40)
c
C     write(unit=*,fmt='(/,3a)') " +++ END PROGRAM ", __FILE__," +++ "
      end
