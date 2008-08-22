c ifort zonalmean.f -o testzone -m64 -I/usr/local/netcdf-64bits-ifort-gcc/include -L/usr/local/netcdf-64bits-ifort-gcc/lib -lnetcdf


c****
      program testzone
      use ZONAL_COM
      implicit none
      include 'netcdf.inc'
      real*8 :: tcub(ic,jc,ndomains),zonal_mean_loc(ndomains,jm)
      real*8 :: area_latband_loc(ndomains,jm),area_latband(jm)
      real*8 :: zonal_mean(jm)
      integer :: itile,fid,vid,srt(3),cnt(3),status,ikey,jlat,i,j,itest
      character*200 :: infi,infile
      character*1 :: istr 
      integer, dimension(ndomains) :: keymax_un,keymax_fo

      call init_xgrid_unrolled(im,jm,ic,jc,ncells,ndomains,
     *     az11,az21,az31,az41,az51,keymax_un)
      write(*,*) "ncells=",ncells

      call init_xgrid_loop(im,jm,ic,jc,ncells,ndomains,
     *     az12,az22,az32,keymax_fo)

c      
c     read atmosphere data and calculate zonal mean for each domain
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
        
        call zonalmean_cs_unrolled(tcub(:,:,itile),
     *       1,ic,1,jc,jm,
     *       az11(itile,:),az21(itile,:),az31(itile,:),
     *       az41(itile,:),az51(itile,:),
     *       zonal_mean_loc(itile,:),
     *       area_latband_loc(itile,:),keymax_un(itile) )

      enddo

      
      do jlat=1,jm
         area_latband(jlat)=0.d0
         do itile=1,ndomains
            area_latband(jlat)=area_latband(jlat)
     *           +area_latband_loc(itile,jlat)
         enddo
         write(25,*) area_latband(jlat)
      enddo

      do jlat=1,jm
         zonal_mean(jlat)=0.d0
         do itile=1,ndomains
            if (area_latband(jlat) .gt. 1.d-14) then
               zonal_mean(jlat)=zonal_mean(jlat)
     *              +zonal_mean_loc(itile,jlat)
     *              /area_latband(jlat)
            endif
         enddo
         write(26,*) zonal_mean(jlat)
      enddo


c     zonal mean using folded loop

      do itile=1,ndomains
         itest=0
         zonal_mean_loc(itile,:)=0.d0
         area_latband_loc(itile,:)=0.d0

         do i=1,ic
            do j=1,jc
               call zonalmean_cs_loop(tcub(i,j,itile),i,j,ic,jc,jm,
     *              az12(itile,:),az22(itile,:),az32(itile,:),
     *              zonal_mean_loc(itile,:),area_latband_loc(itile,:),
     *              keymax_fo(itile),itest)
            enddo
         enddo
         write(*,*) "itest=",itest
      enddo
c

     
      do jlat=1,jm
         area_latband(jlat)=0.d0
         do itile=1,ndomains
            area_latband(jlat)=area_latband(jlat)
     *           +area_latband_loc(itile,jlat)
         enddo
         write(27,*) area_latband(jlat)
      enddo

      do jlat=1,jm
         zonal_mean(jlat)=0.d0
         do itile=1,ndomains
            if (area_latband(jlat) .gt. 1.d-14) then
               zonal_mean(jlat)=zonal_mean(jlat)
     *              +zonal_mean_loc(itile,jlat)
     *              /area_latband(jlat)
            endif
         enddo
         write(28,*) zonal_mean(jlat)
      enddo

      end program testzone
c****



c****
c     Initialize exchange grid for unrolled loops
c****
      subroutine init_xgrid_unrolled(im,jm,ic,jc,ncells,
     *     ndomains,az11,az21,az31,az41,az51,keymax)
      implicit none
      include 'netcdf.inc'

      integer  :: ncells
      integer :: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: az11,az21
      integer, pointer, dimension(:,:) :: az31,az41
      real*8, pointer, dimension(:,:) :: az51
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
c*     array az11 : ikey->latitude_index
c*     array az21 : ikey->xcell_index
c*     array az31:  ikey->i
c*     array az41:  ikey->j
c      array az51:  ikey->xcell_area
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

      allocate(az11(ndomains,maxkey))
      allocate(az21(ndomains,maxkey))
      allocate(az31(ndomains,maxkey))
      allocate(az41(ndomains,maxkey))
      allocate(az51(ndomains,maxkey))

      az11(:,:)=0
      az21(:,:)=0
      az31(:,:)=0
      az41(:,:)=0
      az51(:,:)=0.d0

      do idomain=1,ndomains
         checkarea=0.d0
         ikey=1
         do jlat=1,jm
            do n=1,ncells
               itile=tile(n)
               if (itile .eq. idomain) then
                  j=ijlatlon(2,n)
                  if (j .eq. jlat) then
                     icc=ijcub(1,n)
                     jcc=ijcub(2,n)
                     az11(idomain,ikey)=j
                     az21(idomain,ikey)=n
                     az31(idomain,ikey)=icc
                     az41(idomain,ikey)=jcc
                     az51(idomain,ikey)=xgrid_area(n)
                     checkarea=checkarea+xgrid_area(n)
                     ikey=ikey+1
                  endif
               endif
            enddo
         enddo
         write(*,*) "checkarea=",checkarea
      enddo

      end subroutine init_xgrid_unrolled


c****
c     Initialize exchange grid inside do loop
c****
      subroutine init_xgrid_loop(im,jm,ic,jc,ncells,ndomains,
     *     az12,az22,az32,keymax)
      implicit none
      include 'netcdf.inc'

      integer  :: ncells
      integer :: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: az12,az22
      real*8, pointer, dimension(:,:) :: az32
      integer, dimension(6) :: keymax 

      integer :: status,fid,n,vid,ikey,jlat,maxkey
      integer :: itile,j,idomain,icc,jcc,iic,jjc,index,indexc
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
c*     ikey=ikey+1
c*     index=(iic-1)*jcmax+jjc
c*     array az12 : ikey->index
c*     array az22 : ikey->jlat
c*     array az32 : ikey->xgrid_area
c****

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

      allocate(az12(ndomains,maxkey))
      allocate(az22(ndomains,maxkey))
      allocate(az32(ndomains,maxkey))

      az12(:,:)=0
      az22(:,:)=0
      az32(:,:)=0.d0

      do idomain=1,ndomains
         checkarea=0.d0
         ikey=1
         do iic=1,ic
            do jjc=1,jc
               index=(iic-1)*jc+jjc
               do n=1,ncells
                  itile=tile(n)
                  if (itile .eq. idomain) then
                     icc=ijcub(1,n)
                     jcc=ijcub(2,n)
                     indexc=(icc-1)*jc+jcc
                     if ( index .eq. indexc ) then
                        jlat=ijlatlon(2,n)
                        az12(idomain,ikey)=index
                        az22(idomain,ikey)=jlat
                        az32(idomain,ikey)=xgrid_area(n)      
                        ikey=ikey+1
                        checkarea=checkarea+xgrid_area(n)
                     endif
                  endif
               enddo
            enddo
         enddo
         write(*,*) "ikey=",ikey-1
         write(*,*) "->checkarea=",checkarea
      enddo

      end subroutine init_xgrid_loop
c****



c****
c     Calculate Zonal mean for current domain - cubed sphere version - unrolled loop    
c****
      subroutine zonalmean_cs_unrolled(acub,i0,i1,j0,j1,jm,
     *     az1,az2,az3,az4,az5,zonal_mean,area_latband,maxkey)
c     INPUT:    acube          array to be summed
c               i0,i1,j0,j1    boundaries of 2d domain
c               jm             latitude max
c               azonal1,azonal2,...
c               azonal3,azonal4,...
c               azonal5        contains reordered xarea 
c               maxkey         size of azonal* arrays

      implicit none
      include 'netcdf.inc'
      integer, intent(in) :: i0,i1,j0,j1,jm
      integer, intent(in) :: maxkey
      integer :: az1(maxkey),az2(maxkey),az3(maxkey),az4(maxkey)
      real*8 :: az5(maxkey)
      real*8 :: acub(i1-i0+1,j1-j0+1),zonal_mean(jm)
      real*8 :: area_latband(jm)
      integer :: ikey,jlat,ic,jc

      do jlat=1,jm
c         write(*,*) "jlat=",jlat
         zonal_mean(jlat)=0.d0
         area_latband(jlat)=0.d0
         do ikey=1,maxkey
c            write(*,*) "j=",az1(ikey)
            if (az1(ikey) .eq. jlat) then
               ic=az3(ikey)
               jc=az4(ikey)

c     TODO: hardcode calculation of area_latband in init_xgrid or in different subroutine
               area_latband(jlat)=area_latband(jlat)
     *              +az5(ikey)

               zonal_mean(jlat)=zonal_mean(jlat)+az5(ikey)
     *              *acub(ic,jc)
            endif
         enddo
      enddo

c
c     TODO : might be accelerated using matrix-matrix product, kronecker matrix operator and precomputed
c     matrices...Still need to find the right formula 
c
       
      end subroutine zonalmean_cs_unrolled
c****



c****
c     Calculate Zonal mean for current domain - cubed sphere version     
c****
      subroutine zonalmean_cs_loop(value,iinput,jinput,ic,jc,jm,
     *     az1,az2,az3,zonal_mean,area_latband,maxkey,itest)
c     INPUT:    value          value to be added to zonal mean
c               jm             latitude max
c               IT             index of quantity
c               az1,az2,az3    contains reordered xarea 
c               maxkey         size of az* arrays

      implicit none
      include 'netcdf.inc'
      integer, intent(in) :: iinput,jinput,ic,jc,jm
      integer, intent(in) :: maxkey
      integer :: az1(maxkey),az2(maxkey)
      real*8 :: az3(maxkey)
      real*8 :: value,zonal_mean(jm),area_latband(jm),zstore,astore
      integer :: ikey,jlat,current_index,index,itest

      current_index=jc*(iinput-1)+jinput

      do ikey=1,maxkey
         index=az1(ikey)            
         if (index .eq. current_index) then
            jlat=az2(ikey) 
            zstore=zonal_mean(jlat)
            zonal_mean(jlat)=zstore+az3(ikey)*value
            astore=area_latband(jlat)
            area_latband(jlat)=astore+az3(ikey)
            itest=itest+1
         endif
      enddo


      end subroutine zonalmean_cs_loop
c****



c****
c     Calculate Zonal mean for current domain - latitude-longitude version
c****
c      subroutine zonalmean_latlon(aj,X,ftype,zonal_mean)
c
c     add increment to zonal mean : aj(j)=aj(j)+X(i,j)*ftype(i,j)
c
c      implicit none
c      real*8, intend(in) :: ftype,X
c      real*8, intend(out) :: aj
c
c      aj=aj+X*ftype
c
c      end subroutine zonalmean_latlon
      
