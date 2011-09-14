      program avg
c
c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c --- output key monthly mean variables in mon_[runid]_[decade].txt
c --- output key annual  mean variables in ann_[runid]_[decade].txt
c --- output MOC in lat/rho space in 4 basins averaged over year ny1:ny2  
c     in avg_ov_[runid]_[decade].txt: flux(idm,kdm,4)
c --- output northward heatflux as a function of lat in "heatfl(idm,4)"
c     in 4 basins averaged over year ny1:ny2 in avg_hf_[runid]_[decade].txt
c --- Last index in flux & heatfl: 1: Atl; 2: Indian; 3: Pac; 4: global
c --- Setting rhodot to true will remove model trend during ny1:ny2 period
c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      use const_proc,   only: path0,path1,path2,hycomtopo,latlonij
     .    ,basinmask,flnmcoso,flnmo2a,runid,ny1,ny2,spcifh,julian
     .    ,idrk1,idrk2,jdrk1,jdrk2,indoi,indoj1,indoj2,rhodot,amon
     .    ,monlg,rho,monave_convert,solo_convert,cnvert,timav,mo1,mo2
      use hycom_arrays, only: srfhgt,dpmixl,covice,depths,scp2
     .    ,u,v,dp,p,temp,saln,th3d,tracer,uflxav,vflxav,diaflx
     .    ,uflx,vflx,ubavg,vbavg,alloc_hycom_arrays,latij
     .    ,temav,salav,ubavg,vbavg,ubavav,vbavav
      use hycom_dimen
c
      implicit none
c
      integer :: ntime
      integer :: ny,num,m,k00,ia
      real*8 :: tsum,ssum,sstsum,ssssum,arean3
     . ,annsst,annsss,annssh,anntem,annsal,annht,thin,ssh0
     . ,anndh,anndf,anndb,db,trc
     . ,glbsal,glbtem,glbdep,sum1,sum2,sum3,area,avgbot,vol

      real :: iceextn,iceexts,nino3,day0,day1,flxmax,x1,x2,thrufl,tinvrs
      real, allocatable :: pinit(:,:,:),pfinl(:,:,:),lat(:),flux(:,:,:)
     .       ,sunda(:),heatfl(:,:)
      integer, allocatable :: im(:,:)
c
      real :: dpavav(kdm), heatot
      real, allocatable :: year(:), dpav(:,:)
      character flnm*132,flnmout*30
c
      integer mo,dcd,mon1,i70,i45,ieq,status
      logical :: lexist
c
      character(len=9) :: ayears  ! string n1-n2 example "1905-1955"

      character(len=*), parameter :: FMT1=
     & '(a,         t10,a,     t24,a,     t38,a,      t52, a,      
     &  t66,a,      t80,a,     t94,a,     t108,a,     t122,a)' 
C
      character(len=*), parameter :: FMT2=
     & '(f7.2,sp,t10,es12.5,t24,es12.5,t38,es12.5, t52, es12.5,
     &  t66,es12.5, t80,es12.5,t94,es12.5,t108,es12.5,t122,es12.5)'
c
      character(len=*), parameter :: FMT3=
     & '(a,         t10,a,     t22,a,     t34,a,      t46, a,      
     &  t58,a,      t70,a,     t82,a,     t94,a,      t106,a)' 
C
      character(len=*), parameter :: FMT4=
     & '(f7.2,sp,t10,es10.3,t22,es10.3,t34,es10.3, t46, es10.3,
     &  t58,es10.3, t70,es10.3,t82,es10.3,t94,es10.3,t106,es10.3)'

      character(len=*), parameter :: FMT5=
     & '(1x,a,t6,a,t16,a,t28,a,t40,a,t52,a)'
      character(len=*), parameter :: FMT6=
     & '(1x,i3,sp,t6,f6.1,t16,es10.3,t28,es10.3,t40,es10.3,t52,es10.3)'

      namelist /hdiag_nml/ path0, path1, path2,
     . hycomtopo, latlonij, basinmask, flnmcoso, flnmo2a,
     . runid, ny1, ny2, monave_convert,solo_convert

      open (10,file="hdiag.nml")
      read (10,nml=hdiag_nml)
      write (*,nml=hdiag_nml)
      close(10)
c
      write(*,'(3a,i4,a,i4)')
     .   'processing RunId=',trim(runid),' from yr ',ny1,' to ',ny2
      write(*,'(a,i2)') 'number of tracers =',ntrcr

      ntime=(ny2-ny1+1)*12
      allocate(  year(ntime),dpav(ntime,kdm) )
c
      call alloc_hycom_arrays
      call alloc_hycom_dimen
      allocate (lat(idm),pinit(idm,jdm,kdm+1),pfinl(idm,jdm,kdm+1)
     .    ,flux(idm,kdm,4),im(idm,jdm),heatfl(idm,4),sunda(kdm+1)
     .    ,stat=status)
      if (status/=0) stop 'wrong allocate1'

c --- determine do-loop limits for u,v,p,q points
      call gtdpth(depths,im)
      call bigrid(depths)
c
c --- determine mesh size
      call meshsz
      avgbot=0.
      area=0.
c
      do 10 j=1,jdm
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
      avgbot=avgbot+depths(i,j)*scp2(i,j)
 10   area=area+scp2(i,j)
      avgbot=avgbot/area
      write (*,104) avgbot,area
 104  format(' mean basin depth (m) and area (10^6 km^2):',f9.1,
     .       -12p,f9.1)
c
c --- read archive data
c
      dcd=ny1/10
      if (runid(1:1).eq.' ') stop 'empty runid'
      if (dcd.lt.001 .or. dcd.gt.930) then
        print *,' wrong decade=',dcd
        stop 'wrong decade'
      endif
      write(ayears,'(i4.4,a1,i4.4)') ny1,'-',ny2

c
      do i=1,idm
      lat(i)=latij(i,340,3)
C     write(*,'(a,i3,a,E12.5)') "lat(",i,")=",lat(i)
      enddo
c
      do i=2,idm-1
      if (lat(i+1).lt.70. .and. lat(i).ge.70.) i70=i
      if (lat(i+1).lt.45. .and. lat(i).ge.45.) i45=i
      if (lat(i+1).lt. 0. .and. lat(i).ge. 0.) ieq=i
      enddo
c
      write(flnmout,'(5a)') 'mon_',trim(runid),'_',ayears,'.txt'
      open(301,file=trim(path2)//flnmout, 
     +     form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout) 

      write(flnmout,'(5a)') 'ann_',trim(runid),'_',ayears,'.txt'
      open(302,file=trim(path2)//flnmout,
     +     form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout) 

      write(flnmout,'(5a)') 'avg_ov_',trim(runid),'_',ayears,'.txt'
      open(303,file=trim(path2)//flnmout,
     +     form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout) 

      write(flnmout,'(5a)') 'avg_hf_',trim(runid),'_',ayears,'.txt'
      open(304,file=trim(path2)//flnmout,
     +     form='formatted',status='unknown')
      write(*,'(a,/,a)') 'Open file for writing:',trim(flnmout) 

      write(301,fmt=FMT1) 
     & "Time  ","NINO3","Ice Extent","Ice Extent","Ocean Heat",
     & "Sea Surface","Sea Surface", "Global Ocean","Global Ocean"
      write(301,fmt=FMT1) 
     & "---- ","Index","Arctic","Antarctic","Content",
     & "Temperature","Salinity","Temperature","Salinity"
      write(301,fmt=FMT1) 
     & "Year","    ","Mln.Sq.km","Mln.Sq.km","*1.E6 J/m2",
     & "degC","PSU","degC","PSU"
c
      write(302,fmt=FMT3) 
     & "Time  ","SST","SSS","Tavrg","Savrg","SSH","Ocean Heat", 
     & "Atl (45N)","Indonesian","Drake"
      write(302,fmt=FMT3) 
     & "---- ","Surf_Temp","Surf_Saln","Glob_Temp",
     & "Glob_Saln","Srf_Hght","Content","Max Overt",
     & "Throughfl","Passage"
      write(302,fmt=FMT3) 
     & "Year","degC","PSU","degC","PSU","cm","*1.E6 J/m2","Sv","Sv","Sv"
C
      write(304,'(a)') "North Poleward Ocean Heat Transport" 
      write(304,fmt=FMT5) 
     & "Num","Latitude","Atlantic","Indian","Pacifiq","Global"
      write(304,fmt=FMT5) " ","degr(N)","pWatts","pWatts",
     &  "pWatts","pWatts"

      n=0

      do 151 ny=ny1,ny2

      annsst=0.
      annsss=0.
      anntem=0.
      annsal=0.
      annssh=0.
      annht =0.
      ubavav(:,:)=0.
      vbavav(:,:)=0.
      uflxav(:,:,:)=0.
      vflxav(:,:,:)=0.
c
      timav=.true.
      cnvert=.false.
      do 152 mo=mo1,mo2
      n=n+1
      write(flnm,'(2a,i4.4,2a)') trim(path1),amon(mo),ny
     .  ,'.out',trim(runid)

      inquire(file=flnm,exist=lexist)
      if( .NOT. lexist ) then
        write(*,'(3a)') "!!! ATTENTION !!! file=",trim(flnm),
     +          " does not exist !!!!" 
        write(*,'(3a)') " Calculations continue skip this file "
        cycle
      end if

      write (*,'(2a)') 'reading: ',trim(flnm)
      call get_time(flnm,year(n))
c
      call getdat(flnm,day0,day1,lexist)
      if (.not.lexist) stop '(file open or read error)'
c
      write(*,'(a,f9.2)') 'cpl model yr =',year(n)
c
c --- compute total icea volume and ice extent (for each hemisphere)
c
      iceextn=0.
      iceexts=0.
      ssh0   =0.
      sstsum =0.
      ssssum =0.
      tsum   =0.
      ssum   =0.
      trc    =0.
      nino3  =0.
      arean3 =0.
      vol=0.
      if (idm.ne.387.and.jdm.ne.360) stop 'reset nino3 domain'
      do 5 j=1,jdm
      do 5 k=1,kdm
      do 5 l=1,isp(j)
      do 5 i=ifp(j,l),ilp(j,l)
      if (k.eq.1) then
        if (i.lt.equat) then
          iceextn=iceextn+covice(i,j)*scp2(i,j)
        else
          iceexts=iceexts+covice(i,j)*scp2(i,j)
        end if
c
        ssh0   =   ssh0+srfhgt(i,j)*scp2(i,j)      ! srfhgt in cm
        sstsum = sstsum+temp(i,j,1)*scp2(i,j)
        ssssum = ssssum+saln(i,j,1)*scp2(i,j)
c
c --- nino3 index averaged in 5deg S, 5deg N, 150 W - 90 W
        if (i.ge.227 .and. i.le.259 .and. j.ge.210 .and. j.le.270)then
          nino3=nino3+temp(i,j,1)*scp2(i,j)
          arean3=arean3+scp2(i,j)
        endif
      endif  ! k=1
c 
      tsum=tsum+temp(i,j,k)*dp(i,j,k)*scp2(i,j)
      ssum=ssum+saln(i,j,k)*dp(i,j,k)*scp2(i,j)
      trc=trc+tracer(i,j,k,1)*dp(i,j,k)*scp2(i,j)
      vol=vol+dp(i,j,k)*scp2(i,j)
 5    continue
c
      heatot=spcifh*rho*tsum/area
      write(301,fmt=FMT2)
     . year(n),nino3/arean3,iceextn*1.e-12,iceexts*1.e-12,heatot*1.e-6
     . ,sstsum/area,ssssum/area,tsum/vol,ssum/vol
c
c --- save annual field
      mon1=monlg(mod(n-1,12)+1)
      annsst=annsst+sstsum*mon1/julian
      annsss=annsss+ssssum*mon1/julian
      annssh=annssh+  ssh0*mon1/julian
      anntem=anntem  +tsum*mon1/julian
      annsal=annsal  +ssum*mon1/julian
      annht=annht  +heatot*mon1/julian
c
c --- calculate flux
      timav=.true.
      call getdat(flnm,day0,day1,lexist)
      if (.not.lexist) stop 'file open or read error'
c
      do 13 j=1,jdm
      do 13 i=1,idm
      ubavav(i,j)=ubavav(i,j)+ubavg(i,j,1)*mon1/julian
      vbavav(i,j)=vbavav(i,j)+vbavg(i,j,1)*mon1/julian
      do 13 k=1,kdm
      uflxav(i,j,k)=uflxav(i,j,k)+uflx(i,j,k)		! uflx: Sv*intvl
 13   vflxav(i,j,k)=vflxav(i,j,k)+vflx(i,j,k)		! vflx: Sv*intvl
c
 152  continue
      uflxav(:,:,:)=uflxav(:,:,:)/(365.*86400.)		! => annual in Sv
      vflxav(:,:,:)=vflxav(:,:,:)/(365.*86400.)		! => annual in Sv
c
      flux(:,:,:)=0.
      do 181 j=1,jdm
      do 181 i=1,idm
      if(im(i,j).eq.1.or.im(i,j).eq.2) then          ! Atlantic
        do k=1,kdm
          flux(i,k,1)=flux(i,k,1)-uflxav(i,j,k)
        enddo
      endif
 181   continue
c
      do 184 k=2,kdm
      do 184 i=1,idm
 184  flux(i,k,1)=flux(i,k,1)+flux(i,k-1,1)
c
      i=i45                ! get max overturning rate at 45N in Atlantic
      flxmax=-999.
      do 185 k=1,kdm
      if (flux(i,k,1).gt.flxmax) then
        flxmax=flux(i,k,1)
        k00=k
      endif
 185  continue
c     write(*,*) ' yr n=',n,' flxmax_i45 =',i45,flxmax,' at k=',k00
      x1=thrufl(idrk1,jdrk1,idrk2,jdrk2,'(Drake Passage)')
      x2=thrufl(indoi,indoj1,indoi,indoj2,'(Indonesia)')
c
      write(302,fmt=FMT4) year(n)-0.5 
     .   ,annsst/area,annsss/area,anntem/vol,annsal/vol
     .   ,annssh/area,annht*1.e-6,flxmax,-x2,x1
 151  continue
c
      temav(:,:,:) =0.
      salav(:,:,:) =0.
      uflxav(:,:,:)=0.
      vflxav(:,:,:)=0.
      ubavav(:,:)=0.
      vbavav(:,:)=0.
      heatfl(:,:)=0.
c
      write(*,'(2(a,i6))') 'overturning stream function is averaged over
     . years',ny1,' -',ny2
      timav=.true.
c
      do 153 ny=ny1,ny2
      do 154 mo=mo1,mo2
      write(flnm,'(2a,i4.4,2a)')trim(path1),amon(mo),ny
     .  ,'.out',trim(runid)
      write (*,'(2a)') 'reading: ',trim(flnm)
      mon1=monlg(mod(mo-1,12)+1)
c
c --- read archive data
      call getdat(flnm,day0,day1,lexist)
      if (.not.lexist) stop 'file open or read error'
c
      do 113 i=1,idm
      ia=max(1,i-1)
      do 113 j=1,jdm
      ubavav(i,j)=ubavav(i,j)+ubavg(i,j,1)*mon1
      vbavav(i,j)=vbavav(i,j)+vbavg(i,j,1)*mon1
      do 113 k=1,kdm
      temav(i,j,k)=temav(i,j,k)+temp(i,j,k)*mon1
      salav(i,j,k)=salav(i,j,k)+saln(i,j,k)*mon1
      uflxav(i,j,k)=uflxav(i,j,k)+uflx(i,j,k)		! Sv*intvl
      vflxav(i,j,k)=vflxav(i,j,k)+vflx(i,j,k)		! Sv*intvl
      heatfl(i,4)=heatfl(i,4)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      if (im(i,j).eq.1.or.im(i,j).eq.2) then ! Atlantic
        heatfl(i,1)=heatfl(i,1)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      elseif (im(i,j).eq.3.or.im(i,j).eq.4) then ! Indian
        heatfl(i,2)=heatfl(i,2)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      elseif (im(i,j).eq.5.or.im(i,j).eq.6) then ! Pacific
        heatfl(i,3)=heatfl(i,3)+(temp(i,j,k)+temp(ia,j,k))*uflx(i,j,k)
      endif
 113  continue
c
 154  continue
 153  continue   ! ny=ny1,ny2
      uflx(:,:,:)=uflxav(:,:,:)/((ny2-ny1+1)*julian*86400.)  ! => mean of ny1~ny2 in Sv
      vflx(:,:,:)=vflxav(:,:,:)/((ny2-ny1+1)*julian*86400.)  ! => mean of ny1~ny2 in Sv
      ubavav(:,:)=ubavav(:,:)/((ny2-ny1+1)*julian)       ! => mean of ny1~ny2
      vbavav(:,:)=vbavav(:,:)/((ny2-ny1+1)*julian)       ! => mean of ny1~ny2
      heatfl(:,:)=heatfl(:,:)/((ny2-ny1+1)*julian*86400) ! => mean of ny1~ny2
c
      flux(:,:,:)=0.
c --- global domain
      do j=1,jdm
      do k=1,kdm
      do i=1,idm
        flux(i,k,4)=flux(i,k,4)-uflx(i,j,k)
      enddo
      enddo
      enddo

c --- each basin
      do 81 j=1,jdm
      do 81 i=1,idm
      if (im(i,j).eq.1.or.im(i,j).eq.2) then ! Atlantic
        do k=1,kdm
          flux(i,k,1)=flux(i,k,1)-uflx(i,j,k)
        enddo
      elseif (im(i,j).eq.3.or.im(i,j).eq.4) then ! Indian
        do k=1,kdm
          flux(i,k,2)=flux(i,k,2)-uflx(i,j,k)
        enddo
      elseif (im(i,j).eq.5.or.im(i,j).eq.6) then ! Pacific
        do k=1,kdm
          flux(i,k,3)=flux(i,k,3)-uflx(i,j,k)
        enddo
      endif
 81   continue
c
      do 84 l=1,4
      do 84 i=1,idm
c --- convert to petawatt
      heatfl(i,l)=-.5*heatfl(i,l)*spcifh*rho * 1.e-9          !  N-ward > 0
      do 84 k=2,kdm
 84   flux(i,k,l)=flux(i,k,l)+flux(i,k-1,l)     ! vertical integral in k
c
      i=i45                ! get max overturning rate at 45N in Atlantic
      flxmax=-999.
      do 85 k=1,kdm
      if (flux(i,k,1).gt.flxmax) then
        flxmax=flux(i,k,1)
        k00=k
      endif
 85   continue
      x1=thrufl(idrk1,jdrk1,idrk2,jdrk2,'(Drake Passage)')
      x2=thrufl(indoi,indoj1,indoi,indoj2,'(Indonesia)')
      x2=-x2               ! take only absolute value
c     write(*,'(a,i4,f6.2,a,i2,a,2f6.1)')
c    . 'chk flxmax_i45 =',i45,flxmax,' at k=',k00,'; Drake/Indo=',x1,x2

c --- diagnose indonesian throughflow
c
      sunda=0.
      do 26 k=1,kdm
      i = indoi
      do 35 j=indoj1,indoj2
 35   sunda(k)=sunda(k)+uflx(i,j,k)
 26   sunda(k+1)=sunda(k+1)+sunda(k)
c
c --- subtract out portion due to indonesian throughflow
      do 39 k=1,kdm
      do 39 i=indoi+1,idm
      flux(i,k,2)=flux(i,k,2)+sunda(k)                !  Indian
 39   flux(i,k,3)=flux(i,k,3)-sunda(k)                !  Pacific
c
      if (rhodot) then
        tinvrs=1./((ny2-ny1+1.)*365.*86400.)
c
c --- determine final pressure field (needed for finding rho-dot)
        pfinl(:,:,:)=p(:,:,:)

        write(flnm,'(5a,i3,a,i3,2a,i4.4,2a)')trim(path1),trim(runid)
     .   ,'/out',trim(runid),'_',dcd,'0_',dcd,'9/','DEC',ny1-1
     .   ,'.out',trim(runid)
        write(*,*) ' input 0 file=',trim(flnm)
c --- determine initial pressure field (needed for finding rho-dot)
        call getdat(flnm,day0,day1,lexist)
        if (.not.lexist) stop '(file open or read error)'
c
c --- subtract out portion of flux attributable to interface displacement
        do 34 j=1,jdm
        do 34 k=2,kdm
        do 34 i=1,idm
        flux(i,k,4)=flux(i,k,4)
     .   +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        if (im(i,j).eq.1.or.im(i,j).eq.2) then		! Atlantic
          flux(i,k,1)=flux(i,k,1)
     .     +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        elseif (im(i,j).eq.3.or.im(i,j).eq.4) then	! Indian
          flux(i,k,2)=flux(i,k,2)
     .     +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        elseif (im(i,j).eq.5.or.im(i,j).eq.6) then	! Pacific
          flux(i,k,3)=flux(i,k,3)
     .     +(p(i,j,k)-pfinl(i,j,k))*scp2(i,j)*tinvrs*1.e-6	! => Sv
        endif
 34     continue
      end if
c
      do l=1,4
      write(303,'(387f6.1)') ((flux(i,k,l),i=1,idm),k=1,kdm)
      end do
      close (303)
      do i=1,idm
      write(304,fmt=FMT6) i,lat(i),(heatfl(i,k),k=1,4)
      end do
      close (304)

      stop '(normal finish of avg)'
      end


      subroutine get_time(flnm,year)

! --- extract year and julian day from file name

      character*(*),intent(IN)  :: flnm
      real         ,intent(OUT) :: year

      character*3,parameter :: month(12) =
     .  (/'JAN','FEB','MAR','APR','MAY','JUN',
     .    'JUL','AUG','SEP','OCT','NOV','DEC'/)
      real,parameter :: endmon(0:12) =
     .  (/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.,365./)

      year=0.
      do n=1,12
       i=index(flnm,month(n))
       if (i.gt.0) then
        read(flnm(i+3:i+6),'(f4.0)') year
        year=year+endmon(n)/365.
        exit
       end if
      end do
      print *,'julian day in input file:',year
      return
      end
