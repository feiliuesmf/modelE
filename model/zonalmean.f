c ifort zonalmean.f -o testzone -m64 -I/usr/local/netcdf-64bits-ifort-gcc/include -L/usr/local/netcdf-64bits-ifort-gcc/lib -lnetcdf
c****



c****
      module zonaldata

      integer  :: ncells
      integer,parameter:: im=144,jm=90,ic=48,jc=48,ndomains=6
      integer, pointer, dimension(:,:) :: azonal1,azonal2
      integer, pointer, dimension(:,:) :: azonal3,azonal4
      real*8, pointer, dimension(:,:) :: azonal5
      integer, dimension(6) :: keymax 

      interface init_xgrid
      subroutine init_xgrid(im,jm,ic,jc,ncells,ndomains,
     *     azonal1,azonal2,azonal3,azonal4,azonal5,keymax)
      integer  :: ncells
      integer:: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: azonal1,azonal2
      integer, pointer, dimension(:,:) :: azonal3,azonal4
      real*8, pointer, dimension(:,:) :: azonal5
      integer, dimension(6) :: keymax 
      end subroutine init_xgrid
      end interface 

      end module zonaldata
c****



c****
      program testzone
      use zonaldata
      implicit none
      include 'netcdf.inc'
      real*8 :: tcub(ic,jc,ndomains),zonal_mean_loc(ndomains,jm)
      real*8 :: area_latband_loc(ndomains,jm),area(jc)
      integer :: itile,fid,vid,srt(3),cnt(3),status,ikey,jlat,i
      character*200 :: infi,infile
      character*1 :: istr


      call init_xgrid(im,jm,ic,jc,ncells,ndomains,
     *     azonal1,azonal2,azonal3,azonal4,azonal5,keymax)
      write(*,*) "ncells=",ncells

c
      do i=1,keymax(2)
         write(*,*) "az5=",azonal5(2,i)
      enddo
      
c     read atmosphere data 
c     and calculate zonal mean for each domain
c
      infi='/Users/dgueyffier/fregrid/data/input/atmos_daily.tile'
      do itile=1,ndomains
        write(istr,'(i1)') itile
        infile=trim(infi)//istr//'.nc'
        write(*,*) infile
        status = nf_open(trim(infile),nf_nowrite,fid)
        write(*,*) "status=",status
        status = nf_inq_varid(fid,'t_surf',vid)
        if (status .ne. NF_NOERR) write(*,*) "ERROR"
        srt = (/ 1, 1, 1 /)
        cnt = (/ ic, jc, 1 /)
        status = nf_get_vara_double(fid,vid,srt,cnt,tcub(1,1,itile))
        status = nf_close(fid)
        
        call zonalmean_cs(tcub(:,:,itile),1,ic,1,jc,jm,1,itile,
     *       azonal1(itile,:),azonal2(itile,:),azonal3(itile,:),
     *       azonal4(itile,:),azonal5(itile,:),
     *       zonal_mean_loc(itile,:),
     *       area_latband_loc(itile,:),keymax(itile))

      enddo


      do jlat=1,jm
         area(jlat)=0.d0
         do itile=1,ndomains
c     if (area_latband(jlat) .gt. 1.d-14) then
c     zonal_mean_loc(jlat)=zonal_mean_loc(jlat)
c     *       /area_latband(jlat)
c     endif
c     write(*,*) "zm=",zonal_mean(jlat)
            area(jlat)=area(jlat)+area_latband_loc(itile,jlat)
         enddo
         write(25,*) area(jlat)
      enddo
      

      end program testzone
c****



c****
c     Initialize exchange grid
c****
      subroutine init_xgrid(im,jm,ic,jc,ncells,ndomains,
     *     azonal1,azonal2,azonal3,azonal4,azonal5,keymax)
      implicit none
      include 'netcdf.inc'

      integer  :: ncells
      integer :: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: azonal1,azonal2
      integer, pointer, dimension(:,:) :: azonal3,azonal4
      real*8, pointer, dimension(:,:) :: azonal5
      integer, dimension(6) :: keymax 

      integer :: status,fid,n,vid,ikey,jlat,maxkey
      integer :: itile,j,idomain,icc,jcc
      character*200 :: exchfile
      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:,:) :: ijcub,ijlatlon
      integer, allocatable, dimension(:) :: tile
      real*8 :: checkarea
      character(len=10) :: imch,jmch,icch,jcch 

      write(*,*) "im, jm , ic , jc=",im,jm,ic,jc

      write(imch,'(i10)') im
      write(jmch,'(i10)') jm
      write(icch,'(i10)') ic
      write(jcch,'(i10)') jc

      imch=trim(adjustl(imch))
      jmch=trim(adjustl(jmch))
      icch=trim(adjustl(icch))
      jcch=trim(adjustl(jcch))

      exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *     //"C"//trim(icch)//"-"//trim(jcch)//".nc"
      write(*,*) "filename=",exchfile

c      
c Read weights
c
      status = nf_open(trim(exchfile),nf_nowrite,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN REMAP FILE"
 
      status = nf_inq_dimid(fid,'ncells',vid)
      status = nf_inq_dimlen(fid,vid,ncells)

c
c Allocate arrays with size depending on ncells
c
      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))


      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_get_var_double(fid,vid,xgrid_area)

      status = nf_inq_varid(fid,'tile1',vid)
      status = nf_get_var_int(fid,vid,tile)

      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_get_var_int(fid,vid,ijcub)

      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_get_var_int(fid,vid,ijlatlon)

      status = nf_close(fid)

c****
c*     In order to speed up zonal mean calculation
c*     we create an associative array (equivalent to a C++ multimap)
c*     to store association between latitude_index and (i,j,xcell_index,xarea)
c*     The array is local to each domain - for the moment we 
c*     assume 1 domain per cube face
c*
c*     array azonal1 : ikey->latitude_index
c*     array azonal2 : ikey->xcell_index
c*     array azonal3:  ikey->i
c*     array azonal4:  ikey->j
c      array azonal5:  ikey->xcell_area
c*     keymax:         max(ikey)
c****

      maxkey=0

c     calculate keymax = size of azonal 
      do idomain= 1,ndomains
         ikey=1
         do jlat=1,jm
            do n=1,ncells
               itile=tile(n)
               if (itile .eq. idomain) then
                  j=ijlatlon(2,n)
                  if (j .eq. jlat) then
                     ikey=ikey+1
                  endif
               endif
            enddo
         enddo
         keymax(idomain)=ikey-1
         write(*,*) "keymax=",ikey-1
         maxkey=max(maxkey,ikey-1)
      enddo

      write(*,*) "ndomains maxkey=",ndomains,maxkey

      allocate(azonal1(ndomains,maxkey))
      allocate(azonal2(ndomains,maxkey))
      allocate(azonal3(ndomains,maxkey))
      allocate(azonal4(ndomains,maxkey))
      allocate(azonal5(ndomains,maxkey))

      azonal1(:,:)=0
      azonal2(:,:)=0
      azonal3(:,:)=0
      azonal4(:,:)=0
      azonal5(:,:)=0.d0

      do idomain=1,ndomains
         ikey=1
         do jlat=1,jm
            checkarea=0.d0
            do n=1,ncells
               itile=tile(n)
               if (itile .eq. idomain) then
                  j=ijlatlon(2,n)
                  if (j .eq. jlat) then
                     icc=ijcub(1,n)
                     jcc=ijcub(2,n)
                     azonal1(idomain,ikey)=j
                     azonal2(idomain,ikey)=n
                     azonal3(idomain,ikey)=icc
                     azonal4(idomain,ikey)=jcc
                     azonal5(idomain,ikey)=xgrid_area(n)
                     checkarea=checkarea+xgrid_area(n)
                     ikey=ikey+1
                  endif
               endif
            enddo
            write(*,*) "checkarea=",checkarea
         enddo
      enddo

      end subroutine init_xgrid



c****
c     Calculate Zonal mean for current domain
c****
      subroutine zonalmean_cs(tcub,i0,i1,j0,j1,jm,IT,idomain,
     *     az1,az2,az3,az4,az5,zonal_mean,area_latband,maxkey)
c     INPUT:    tcube          array to be summed
c               i0,i1,j0,j1    boundaries of 2d domain
c               jm             latitude max
c               IT             index of quantity
c               idomain        index of domain
c               azonal1,azonal2,...
c               azonal3,azonal4,...
c               azonal5        contains reordered xarea 
c               maxkey         size of azonal* arrays
      implicit none

      include 'netcdf.inc'
      integer, intent(in) :: i0,i1,j0,j1,jm,IT
      integer, intent(in) :: maxkey,idomain
      integer :: az1(maxkey),az2(maxkey),az3(maxkey),az4(maxkey)
      real*8 :: az5(maxkey)
      real*8 :: tcub(i1-i0+1,j1-j0+1),zonal_mean(j1-j0+1)
      real*8 :: area_latband(j1-j0+1)
      integer :: ikey,jlat,ic,jc

      write(*,*) "*********************idomain=",idomain
      do jlat=1,jm
c         write(*,*) "jlat=",jlat
         zonal_mean(jlat)=0.d0
         area_latband(jlat)=0.d0
         do ikey=1,maxkey
c            write(*,*) "j=",az1(ikey)
            if (az1(ikey) .eq. jlat) then
               ic=az3(ikey)
               jc=az4(ikey)
               area_latband(jlat)=area_latband(jlat)
     *              +az5(ikey)
               zonal_mean(jlat)=zonal_mean(jlat)+az5(ikey)
     *              *tcub(ic,jc)
            endif
         enddo
      enddo

      
      end subroutine zonalmean_cs
