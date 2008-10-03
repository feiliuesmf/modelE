
c****
      program testregrid
cc      call test_zonal_loop()    ! tests parallel zonal mean code, using rolled loops
cc      call test_zonal_unrolled() ! tests parallel zonal mean code, using unrolled loops
      call test_regrid()        ! tests parallel regridding code
cc     call offregrid()    ! tests offline regridding for input file
      end program testregrid


c****
      subroutine test_zonal_unrolled()
      use regrid_com

      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop,domain_decomp
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain
      use fv_grid_utils_mod, only : grid_utils_init
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type
      use fv_control_mod, only : uniform_ppm, c2l_ord

      implicit none
      include 'netcdf.inc'
      integer :: npx,npy,npes,ng,ndims,rootpe
      integer :: commID
      integer, dimension(:), allocatable :: pelist

      integer :: l,npz
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline', ofi
      logical :: non_ortho
      integer, parameter :: nedge=4 ! number of edges on each face
      integer :: nij_latlon,nc2,loctile
      real*8 :: gsum

      real*8, allocatable :: tcub_loc(:,:)
      real*8 :: tcub(ic,jc)
      real*8 :: area_latband(jm),zonal_mean(jm)
      integer :: itile,fid,vid,srt(3),cnt(3),status,ikey,jlat,i,j
      character*200 :: infi,infile
      character*1 :: istr 

c
c     Initialize grid and domain decomposition
c
      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
      write(6,*) 'rootpe= ',mpp_root_pe()
c      npx = 90; npy = npx ! 1x1 resolution
      npx = 48; npy = npx 
      
      ng = 3 ! number of ghost zones required
      call domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

c      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

      rootpe=mpp_root_pe()

c
c     Initialize exchange grid and extract information
c
      call init_xgrid_unrolled()

      write(6,*) "HERE"

ccc    remove these lines - this will be initialize somewhere else
      allocate(tcub_loc(isd:ied,jsd:jed))

      infi='atmos_daily.tile'
      loctile=gid/dom_per_tile
      loctile=loctile+1
      write(istr,'(i1)') loctile
      infile=trim(infi)//istr//'.nc'
      write(*,*) infile
      status = nf_open(trim(infile),nf_nowrite,fid)
      write(*,*) "status=",status
      status = nf_inq_varid(fid,'t_surf',vid)
      if (status .ne. NF_NOERR) write(*,*) "ERROR"
      srt = (/ 1, 1, 1 /)
      cnt = (/ ic, jc, 1 /)
      status = nf_get_vara_double(fid,vid,srt,cnt,tcub(1,1))
      status = nf_close(fid)
      tcub_loc(is:ie,js:je)=tcub(is:ie,js:je)

c      tcub_loc(:,:)=1  !gid
ccc    end remove

      call zonalmean_cs_unrolled(tcub_loc(:,:),jm,
     *       az11(:),az21(:),az31(:),
     *       az41(:),az51(:),
     *       zonal_mean(:),
     *       area_latband(:),maxkey )

      write(6,*) "THERE"

c     should be pre-computed
      call mpp_sum(area_latband,jm)

      call mpp_sum(zonal_mean,jm)
         

      if (gid .eq. rootpe) then
      do jlat=1,jm
         if (area_latband(jlat) .gt. 1.d-14) then
            zonal_mean(jlat)=zonal_mean(jlat)
     *           /area_latband(jlat)
         endif
         write(*,*) "zonal_mean unrolled=",zonal_mean(jlat)
      enddo
      endif

      call mpp_exit()
      deallocate(pelist)
      deallocate(tcub_loc)
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine test_zonal_unrolled
c****



      subroutine test_zonal_loop()
      use regrid_com

      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop,domain_decomp
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain
      use fv_grid_utils_mod, only : grid_utils_init
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type
      use fv_control_mod, only : uniform_ppm, c2l_ord

      implicit none
      include 'netcdf.inc'
      integer :: npx,npy,npes,ng,ndims,rootpe
      integer :: commID
      integer, dimension(:), allocatable :: pelist

      integer :: l,npz
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline', ofi
      logical :: non_ortho
      integer, parameter :: nedge=4 ! number of edges on each face
      integer :: nij_latlon,nc2,loctile
      real*8 :: gsum

      real*8, allocatable :: tcub_loc(:,:)
      real*8 :: tcub(ic,jc)
      real*8 :: area_latband(jm),zonal_mean(jm)
      integer :: itile,fid,vid,srt(3),cnt(3),status,ikey,jlat,i,j
      character*200 :: infi,infile
      character*1 :: istr 

c
c     Initialize grid and domain decomposition
c
      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
      write(6,*) 'rootpe= ',mpp_root_pe()
c      npx = 90; npy = npx ! 1x1 resolution
      npx = 48; npy = npx 
      
      ng = 3 ! number of ghost zones required
      call domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

c      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

      rootpe=mpp_root_pe()

c
c     Initialize exchange grid and extract information
c
      call init_xgrid_loop()

      write(6,*) "HERE"

ccc    remove these lines - this will be initialize somewhere else
      allocate(tcub_loc(isd:ied,jsd:jed))

      infi='atmos_daily.tile'
      loctile=gid/dom_per_tile
      loctile=loctile+1
      write(istr,'(i1)') loctile
      infile=trim(infi)//istr//'.nc'
      write(*,*) infile
      status = nf_open(trim(infile),nf_nowrite,fid)
      write(*,*) "status=",status
      status = nf_inq_varid(fid,'t_surf',vid)
      if (status .ne. NF_NOERR) write(*,*) "ERROR"
      srt = (/ 1, 1, 1 /)
      cnt = (/ ic, jc, 1 /)
      status = nf_get_vara_double(fid,vid,srt,cnt,tcub(1,1))
      status = nf_close(fid)
      tcub_loc(is:ie,js:je)=tcub(is:ie,js:je)

c      tcub_loc(:,:)=1  !gid
ccc    end remove


c     zonal mean inside rolled loop

      area_latband(:)=0.0
      zonal_mean(:)=0.0

      do i=is,ie
         do j=js,je
            call zonalmean_cs_loop(tcub_loc(i,j),i,j,jc,jm,
     *           az12(:),az22(:),az32(:),
     *           zonal_mean(:),area_latband(:),maxkey)
            enddo
         enddo

c     sum over all domains (not reducing rank)

      call mpp_sum(area_latband,jm)
      call mpp_sum(zonal_mean,jm)


      if (gid .eq. rootpe) then
      do jlat=1,jm
         if (area_latband(jlat) .gt. 1.d-14) then
            zonal_mean(jlat)=zonal_mean(jlat)
     *           /area_latband(jlat)
         endif
         write(*,*) "zonal_mean rolled=",zonal_mean(jlat)
      enddo
      endif

      call mpp_exit()
      deallocate(pelist)
      deallocate(tcub_loc)
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine test_zonal_loop
c****




      subroutine init_xgrid_unrolled()
c****
c     Initialize exchange grid for unrolled loops
c****
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : gid

      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'

      integer :: status,fid,n,vid,ikey,jlat,nc2,rootpe
      integer :: itile,j,idomain,loctile,ierr,icc,jcc
      character*200 :: exchfile
      real*8 :: checkarea
      character(len=10) :: imch,jmch,icch,jcch 
       
      rootpe=mpp_root_pe()
      
      if (gid .eq. rootpe) then
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

c      exchfile="exexch.nc"


c      
c Read weights
c
      status = nf_open(trim(exchfile),nf_nowrite,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN REMAP FILE"
 
      status = nf_inq_dimid(fid,'ncells',vid)
      status = nf_inq_dimlen(fid,vid,ncells)

      endif
            
c     
c     Broadcast value of ncells & Allocate arrays with size 
c     depending on ncells on each processor
      

      write(6,*) "HERE INIT  UNROL"

      call MPI_BCAST( ncells, 1, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 


      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))


      if (gid .eq. rootpe) then
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)
         
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         
         status = nf_close(fid)
      endif


c
c     Broadcast x-grid area and indices to all procs
c
      nc2=2*ncells
     
      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     rootpe, MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( tile, ncells, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 


c**** 
c*     In order to speed up zonal mean calculation
c*     we create an associative array (equivalent to a C++ multimap)
c*     to store association between latitude index and (i,j,xcell_index,xarea)
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

c     find current face
      loctile=gid/dom_per_tile
      loctile=loctile+1

      write(*,*) "LOCTILE=",loctile

c     calculate maxkey = size of azonal 
      
      ikey=1
      do jlat=1,jm
         do n=1,ncells
            itile=tile(n)
            if (itile .eq. loctile) then
               j=ijlatlon(2,n)
               if (j .eq. jlat) then
                  ikey=ikey+1
               endif
            endif
         enddo
      enddo
      
      maxkey=ikey-1
       
      allocate(az11(maxkey))
      allocate(az21(maxkey))
      allocate(az31(maxkey))
      allocate(az41(maxkey))
      allocate(az51(maxkey))


      az11(:)=0
      az21(:)=0
      az31(:)=0
      az41(:)=0
      az51(:)=0.d0

      checkarea=0.d0
      ikey=1
      do jlat=1,jm
         do n=1,ncells
            itile=tile(n)
            if (itile .eq. loctile) then
               j=ijlatlon(2,n)
               if (j .eq. jlat) then
                  icc=ijcub(1,n)
                  jcc=ijcub(2,n)
                  az11(ikey)=j
                  az21(ikey)=n
                  az31(ikey)=icc
                  az41(ikey)=jcc
                  az51(ikey)=xgrid_area(n)
                  checkarea=checkarea+xgrid_area(n)
                  ikey=ikey+1
               endif
            endif
         enddo
      enddo

      write(*,*) "checkarea=",checkarea

      end subroutine init_xgrid_unrolled



      subroutine init_xgrid_loop()
c****
c     Initialize exchange grid inside do loop
c****
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : gid

      implicit none
      include 'netcdf.inc'
      include 'mpif.h'

      integer :: status,fid,n,vid,ikey,jlat,nc,nc2,rootpe,
     *     loctile,ierr,icc,jcc
      integer :: itile,j,idomain,iic,jjc,index,indexc
      character*200 :: exchfile
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
      write(*,*) "-->filename=",exchfile

c      exchfile="exexch.nc"

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
      
      rootpe=mpp_root_pe()
      
      if (gid .eq. rootpe) then
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)
         
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         
         status = nf_close(fid)
      endif

c
c     Broadcast x-grid area and indices to all procs
c
      nc2=2*ncells
     
      write(6,*) "BEFMPI ROLLED"

      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     rootpe, MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( tile, ncells, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      write(6,*) "AFT BCAST ROLLED"
c****
c*     ikey=ikey+1
c*     index=(iic-1)*jcmax+jjc
c*     array az12 : ikey->index
c*     array az22 : ikey->jlat
c*     array az32 : ikey->xgrid_area
c****

c     find current face
      loctile=gid/dom_per_tile
      loctile=loctile+1

c     calculate maxkey = size of azonal 
 
      ikey=1
      do jlat=1,jm
         do n=1,ncells
            itile=tile(n)
            if (itile .eq. loctile) then
               j=ijlatlon(2,n)
               if (j .eq. jlat) then
                  ikey=ikey+1
               endif
            endif
         enddo
      enddo
      
      maxkey=ikey-1

 
      allocate(az12(maxkey))
      allocate(az22(maxkey))
      allocate(az32(maxkey))

      write(6,*) "AFT ALOC ROLLED, maxkey",maxkey

      az12(:)=0
      az22(:)=0
      az32(:)=0.d0

      checkarea=0.d0
      ikey=1
      do iic=1,ic
         do jjc=1,jc
            index=(iic-1)*jc+jjc
            do n=1,ncells
               itile=tile(n)
c               write(6,*) "itile loctile",itile,loctile
               if (itile .eq. loctile) then
                  icc=ijcub(1,n)
                  jcc=ijcub(2,n)
                  indexc=(icc-1)*jc+jcc
c               write(6,*) "index indexc",index,indexc
                  if ( index .eq. indexc ) then
                     jlat=ijlatlon(2,n)
                     az12(ikey)=index
                     az22(ikey)=jlat
                     az32(ikey)=xgrid_area(n)      
                     write(6,*) "ikey=",ikey
                     ikey=ikey+1
                     checkarea=checkarea+xgrid_area(n)
                  endif
               endif
            enddo
         enddo
      enddo
      write(*,*) "ikey=",ikey-1
      write(*,*) "->checkarea=",checkarea
      

      end subroutine init_xgrid_loop
c****




      subroutine zonalmean_cs_unrolled(acub,jm,
     *     az1,az2,az3,az4,az5,zonal_mean,area_latband,maxkey)
c****
c     Calculate Zonal mean for current domain - cubed sphere version - unrolled loop    
c
c     INPUT:    acube          array to be summed
c               jm             latitude max
c               azonal1,azonal2,...
c               azonal3,azonal4,...
c               azonal5        contains reordered xarea 
c               maxkey         size of azonal* arrays
c****

      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      implicit none
      include 'netcdf.inc'
      integer, intent(in) :: maxkey,jm
      integer :: az1(maxkey),az2(maxkey),az3(maxkey),az4(maxkey)
      real*8 :: az5(maxkey)
      real*8 :: acub(isd:ied,jsd:jed),zonal_mean(jm)
      real*8 :: area_latband(jm)
      integer :: ikey,jlat,icub,jcub

      do jlat=1,jm
c         write(*,*) "jlat=",jlat
         zonal_mean(jlat)=0.d0
         area_latband(jlat)=0.d0
         do ikey=1,maxkey
c            write(*,*) "j=",az1(ikey)
            if (az1(ikey) .eq. jlat) then
               icub=az3(ikey)
               jcub=az4(ikey)
           if ( ( (icub .le. ie) .and. (icub .ge. is)) .and.
     &          ( (jcub .le. je) .and. (jcub .ge. js)) ) then
               area_latband(jlat)=area_latband(jlat)
     *              +az5(ikey)
               zonal_mean(jlat)=zonal_mean(jlat)+az5(ikey)
     *              *acub(icub,jcub)
               endif
            endif
         enddo
      enddo

c
c     TODO : might be accelerated using matrix-matrix product, kronecker matrix operator and precomputed
c     matrices...Still need to find the right formula 
c
       
      end subroutine zonalmean_cs_unrolled
c****




      subroutine zonalmean_cs_loop(value,iinput,jinput,jc,jm,
     *     az1,az2,az3,zonal_mean,area_latband,maxkey)
c****
c     Calculate Zonal mean for current domain - cubed sphere version     
c
c     INPUT:    value          value to be added to zonal mean
c               jm             latitude max
c               IT             index of quantity
c               az1,az2,az3    contains reordered xarea 
c               maxkey         size of az* arrays
c****

      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      implicit none
      include 'netcdf.inc'
      integer, intent(in) :: iinput,jinput
      integer, intent(in) :: jc,jm   ! jc=max j index of cube face
                                     ! jm =max lat for whole globe 
      integer, intent(in) :: maxkey
      integer :: az1(maxkey),az2(maxkey)
      real*8 :: az3(maxkey)
      real*8 :: value,zonal_mean(jm),area_latband(jm),zstore,astore
      integer :: ikey,jlat,current_index,index

      if ( ( (iinput .le. ie) .and. (iinput .ge. is)) .and.
     &     ( (jinput .le. je) .and. (jinput .ge. js)) ) then

      current_index=jc*(iinput-1)+jinput

      do ikey=1,maxkey
         index=az1(ikey)            
         if (index .eq. current_index) then
            jlat=az2(ikey) 
            zstore=zonal_mean(jlat)
            zonal_mean(jlat)=zstore+az3(ikey)*value
            astore=area_latband(jlat)
            area_latband(jlat)=astore+az3(ikey)
         endif
      enddo

      endif

      end subroutine zonalmean_cs_loop
c****



      subroutine test_regrid()
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop,domain_decomp
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain
      use fv_grid_utils_mod, only : grid_utils_init
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type
      use fv_control_mod, only : uniform_ppm, c2l_ord

      implicit none
      include 'netcdf.inc'

      real*8 :: tlatlon(im,jm),alatlon(im,jm)
      integer :: i,j,n,itile,rootpe,ierr
      integer :: nfaces,npx,npy,npes,ng,ndims
      integer :: commID, status, fid, vid
      integer, dimension(:), allocatable :: pelist
      real*8, dimension(:,:), allocatable :: tcub_loc

      integer :: l,npz
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline', ofi
      logical :: non_ortho
      integer, parameter :: nedge=4 ! number of edges on each face
      integer :: nij_latlon,nc2
      real*8 :: gsum


      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
      write(6,*) 'rootpe= ',mpp_root_pe()
c      npx = 90; npy = npx ! 1x1 resolution
      npx = 48; npy = npx 
      
      ng = 3 ! number of ghost zones required
      call domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

c      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

      rootpe=mpp_root_pe()


c read weights and initialize target array to zero
      call init_regrid(pelist)

c
c     local section 
c

ccc    remove these lines - this will be initialize somewhere else
      allocate(tcub_loc(isd:ied,jsd:jed))
      tcub_loc(:,:)=gid
ccc    end remove

c      write(6,*) tcub_loc(:,:)
      call parallel_regrid_cs2ll(tcub_loc,tlatlon,alatlon)   

      ofi='tstout.nc'
      status = nf_open(trim(ofi),nf_write,fid)
      if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
      status = nf_inq_varid(fid,'lwup_sfc',vid)
      write(*,*) NF_STRERROR(status)
      status = nf_put_var_double(fid,vid,tlatlon)
      write(*,*) NF_STRERROR(status)
     
      status = nf_close(fid)

      call mpp_exit()

      deallocate(pelist)
      deallocate(tcub_loc)
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine test_regrid



      subroutine parallel_regrid_cs2ll(tcub_loc,tlatlon,alatlon)
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid
      implicit none
      real*8 :: tlatlon(im,jm),alatlon(im,jm),
     &     tcub_loc(isd:ied,jsd:jed)
      integer :: loctile  !tile index of current domain
      integer :: n,icub,jcub,i,j,itile,nij_latlon,rootpe
      
      loctile=gid/dom_per_tile
      loctile=loctile+1

c      write(6,*) 'GIIIDD=',gid
c      write(6,*) 'LOCTILE=',loctile

      alatlon(:,:) = 0d0
      tlatlon(:,:) = 0d0

      do n=1,ncells
        icub=ijcub(1,n)
        jcub=ijcub(2,n)
        itile=tile(n)
c        write(6,*) 'ie, icub, is=',ie,icub,is
c        write(6,*) 'je, jcub, js=',je,jcub,js

        if (itile .eq. loctile) then
           if ( ( (icub .le. ie) .and. (icub .ge. is)) .and.
     &          ( (jcub .le. je) .and. (jcub .ge. js)) ) then
              write(6,*) 'xg_area',xgrid_area(n)
              i=ijlatlon(1,n)
              j=ijlatlon(2,n)
              alatlon(i,j) = alatlon(i,j) + xgrid_area(n)
              tlatlon(i,j) = tlatlon(i,j) + xgrid_area(n)
     &             *tcub_loc(icub,jcub)
           endif
        endif
      enddo

c
c     sum all contributions
c
      nij_latlon=im*jm
      call mpp_sum(tlatlon,nij_latlon)

c     should be pre-computed in init_regrid
      call mpp_sum(alatlon,nij_latlon)

c
c     root proc section
c
      if (gid .eq. rootpe) then
      tlatlon(:,:) = tlatlon(:,:)/alatlon(:,:)
      endif

      end subroutine parallel_regrid_cs2ll
c****




      subroutine root_regrid_ll2cs(tlatlon,tcubglob)
c
c     root processor regrids data from lat-lon -> cubbed sphere 
c
      use regrid_com
      use mpp_mod
      use fv_mp_mod, only : gid
      implicit none
      real*8 :: tlatlon(im,jm),tcubglob(1:ic,1:jc,6),
     &     acubglob(1:ic,1:jc,6)
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl,rootpe
      

      rootpe=mpp_root_pe()
      
      if (gid .eq. rootpe) then
         
         acubglob(:,:,:) = 0d0
         tcubglob(:,:,:) = 0d0
         
         do n=1,ncells
            itile=tile(n)
            icc=ijcub(1,n)
            jcc=ijcub(2,n)
            il=ijlatlon(1,n)
            jl=ijlatlon(2,n)
            
            acubglob(icc,jcc,itile) = acubglob(icc,jcc,itile) 
     &           + xgrid_area(n)
            tcubglob(icc,jcc,itile) = tcubglob(icc,jcc,itile) 
     &           + xgrid_area(n)*tlatlon(il,jl)
         enddo
         
         do itile=1,6
            do j=1,jc
               do i=1,ic
                  tcubglob(i,j,itile) = tcubglob(i,j,itile)
     &                 /acubglob(i,j,itile)
               enddo
            enddo
         enddo
         
      endif
      
      end subroutine root_regrid_ll2cs
c****
         


      subroutine init_regrid(pelist)
c     
c     Reads regriding file on root proc then
c     broadcasts the data on all procs
c     
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : gid
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'
      
      integer :: status,fid,n,vid,ikey,jlat,rootpe
      integer :: itile,j,idomain,iic,jjc,index,indexc,nc2
      integer :: ierr
      character*200 :: exchfile
      character(len=10) :: imch,jmch,icch,jcch
      integer, dimension(:), allocatable :: pelist

      
      rootpe=mpp_root_pe()
      
      if (gid .eq. rootpe) then
         
         write(6,*) "GID INIT=",gid
         
         write(imch,'(i10)') im
         write(jmch,'(i10)') jm
         write(icch,'(i10)') ic
         write(jcch,'(i10)') jc
         
         imch=trim(adjustl(imch))
         jmch=trim(adjustl(jmch))
         icch=trim(adjustl(icch))
         jcch=trim(adjustl(jcch))
c         exchfile="exexch.nc"
         exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *        //"C"//trim(icch)//"-"//trim(jcch)//".nc"
         write(*,*) "filename=",exchfile
         
c     
c     Read weights
c     
         status = nf_open(trim(exchfile),nf_nowrite,fid)
         if (status .ne. NF_NOERR) write(*,*) 
     *        "UNABLE TO OPEN REMAP FILE"
         
         status = nf_inq_dimid(fid,'ncells',vid)
         status = nf_inq_dimlen(fid,vid,ncells)
         
      endif
            
c     
c     Broadcast value of ncells & Allocate arrays with size 
c     depending on ncells on each processor
      

      call MPI_BCAST( ncells, 1, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      write(6,*) "ncells GID", ncells,gid

      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))
      
      if (gid .eq. rootpe) then
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)
         
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         
         status = nf_close(fid)
      endif
      
c
c     Broadcast x-grid area and indices to all procs
c
      nc2=2*ncells
     

      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     rootpe, MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( tile, ncells, MPI_INTEGER, rootpe, 
     *     MPI_COMM_WORLD, ierr ) 

      end subroutine init_regrid
c**** 




      subroutine init_regrid_rootpe
c     
c     Reads regriding file on root proc 
c     
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : gid
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'
      
      integer :: status,fid,n,vid,ikey,jlat,rootpe
      integer :: itile,j,idomain,iic,jjc,index,indexc,nc2
      integer :: ierr
      character*200 :: exchfile
      character(len=10) :: imch,jmch,icch,jcch
      integer, dimension(:), allocatable :: pelist

      
      rootpe=mpp_root_pe()
      
      if (gid .eq. rootpe) then
         
         write(6,*) "GID INIT=",gid
         
         write(imch,'(i10)') im
         write(jmch,'(i10)') jm
         write(icch,'(i10)') ic
         write(jcch,'(i10)') jc
         
         imch=trim(adjustl(imch))
         jmch=trim(adjustl(jmch))
         icch=trim(adjustl(icch))
         jcch=trim(adjustl(jcch))
c         exchfile="exexch.nc"
         exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *        //"C"//trim(icch)//"-"//trim(jcch)//".nc"
         write(*,*) "filename=",exchfile
         
c     
c     Read weights
c     
         status = nf_open(trim(exchfile),nf_nowrite,fid)
         if (status .ne. NF_NOERR) write(*,*) 
     *        "UNABLE TO OPEN REMAP FILE"
         
         status = nf_inq_dimid(fid,'ncells',vid)
         status = nf_inq_dimlen(fid,vid,ncells)
         

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
      
      endif
      
      end subroutine init_regrid_rootpe
c**** 


ccc
      subroutine offregrid()
      implicit none
      character*200 :: wtfile,ofile,infi,outfile
      character*200 :: outunformat1,outunformat2,outunformat3 
      character*200 :: outunformat4,outunformat5,outunformat6

      character*1 :: ostr
      character*3 :: irstr

      include 'netcdf.inc'
      integer, parameter :: ncells=31344,nmax=500
      integer, parameter :: im=72,jm=46,imc=48,jmc=48
      real*8 :: xgrid_area(ncells),tot_area
      real*8 :: tcub(imc,jmc,6),acub(imc,jmc,6),area_tile(6)
      real*8 :: tll(im,jm,nmax)
      integer :: tile(ncells)
      integer, dimension(2,ncells) :: ijcub,ijlatlon
      integer :: status,vid,fid,itile,srt(2),cnt(2),n,i,j,ic,jc,il,jl
      character*80 TITLE(nmax)
      real*4 data(im,jm,nmax),tout(imc,jmc,6)
      integer iu1,iu2,iu3,iu4,iu5,iu6,irec,recmax,ir,iunit

c read weights
      wtfile='remapC48N4x5.nc'
      status = nf_open(trim(wtfile),nf_nowrite,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN REMAP FILE"
      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_get_var_double(fid,vid,xgrid_area)
      status = nf_inq_varid(fid,'tile1',vid)
      status = nf_get_var_int(fid,vid,tile)
      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_get_var_int(fid,vid,ijcub)
      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_get_var_int(fid,vid,ijlatlon)
      status = nf_close(fid)


c Read data on latlon grid
      infi='/Users/dgueyffier/fregrid/V72X46.1.cor2_no_crops.ext'
      iunit=15

      open( iunit, FILE=infi, 
     &     FORM='unformatted', STATUS='old')

      irec=1

      do
      read(unit=iunit,END=30) TITLE(irec), data(:,:,irec)
 
      write(*,*) TITLE(irec)

      do i=1,im
         do j=1,jm
            tll(i,j,irec)= data(i,j,irec)
         enddo
      enddo

      irec=irec+1

      enddo
 30   continue

      close(iunit)

      write(*,*) "HERE"
      recmax=irec
c
c     prepare unformatted output
c     
      ofile='/Users/dgueyffier/fregrid/output_fregrid/outll4x5-'
      outunformat1=trim(ofile)//'1.dat'            
      outunformat2=trim(ofile)//'2.dat'            
      outunformat3=trim(ofile)//'3.dat'            
      outunformat4=trim(ofile)//'4.dat'            
      outunformat5=trim(ofile)//'5.dat'            
      outunformat6=trim(ofile)//'6.dat'            
      iu1=21
      iu2=22
      iu3=23
      iu4=24
      iu5=25
      iu6=26
      
      open( iu1, FILE=outunformat1, 
     &     FORM='unformatted', STATUS='new')
      open( iu2, FILE=outunformat2, 
     &     FORM='unformatted', STATUS='new')
      open( iu3, FILE=outunformat3, 
     &     FORM='unformatted', STATUS='new')
      open( iu4, FILE=outunformat4, 
     &     FORM='unformatted', STATUS='new')
      open( iu5, FILE=outunformat5, 
     &     FORM='unformatted', STATUS='new')
      open( iu6, FILE=outunformat6, 
     &     FORM='unformatted', STATUS='new')
      
      do ir=1,recmax
         
c     Regrid
         acub(:,:,:) = 0d0
         tcub(:,:,:) = 0d0
         area_tile(:) = 0d0
         do n=1,ncells
            itile=tile(n)
            ic=ijcub(1,n)
            jc=ijcub(2,n)
            il=ijlatlon(1,n)
            jl=ijlatlon(2,n)
c     write(*,*) tll(il,jl)
            area_tile(itile) = area_tile(itile) + xgrid_area(n)
            acub(ic,jc,itile) = acub(ic,jc,itile) + xgrid_area(n)
            tcub(ic,jc,itile) = tcub(ic,jc,itile) 
     &           + xgrid_area(n)*tll(il,jl,ir)
         enddo
         do itile=1,6
            do j=1,jmc
               do i=1,imc
                  tcub(i,j,itile) = tcub(i,j,itile)/acub(i,j,itile)
               enddo
            enddo
         enddo

c     Output area of each tile
         tot_area=0d0
         do itile=1,6
            write(*,*) "area tile",itile,"=",area_tile(itile)
            tot_area=tot_area+area_tile(itile)
         enddo
         
         write(*,*) "tot area=",tot_area
         
c         
c     Write to netcdf output file
c         
         
         do itile=1,6
            write(ostr,'(i1)') itile
            if (ir .lt. 10) then
               write(irstr,'(i1)') ir
            elseif (ir .lt. 100) then
               write(irstr,'(i2)') ir
            else
               write(irstr,'(i3)') ir
            endif
            outfile=trim(ofile)//ostr//'-'//trim(irstr)
     *           //'.nc'
            
            write(*,*) "ni=",im," nj=",jm
            write(*,*) "outfile =",outfile
            status = nf_open(trim(outfile),nf_write,fid)
            if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
            status = nf_inq_varid(fid,'zsurf',vid)
            write(*,*) NF_STRERROR(status)
            status = nf_put_var_double(fid,vid,tcub(:,:,itile))
            write(*,*) NF_STRERROR(status)
            status = nf_close(fid)
            
         enddo

c
c     write to unformated output
c

         tout(:,:,:)=tcub(:,:,:)

         write(unit=iu1) TITLE(ir), tout(:,:,1)
         write(unit=iu2) TITLE(ir), tout(:,:,2)
         write(unit=iu3) TITLE(ir), tout(:,:,3)
         write(unit=iu4) TITLE(ir), tout(:,:,4)
         write(unit=iu5) TITLE(ir), tout(:,:,5)
         write(unit=iu6) TITLE(ir), tout(:,:,6)
         write(*,*) "ir=",ir
                 
      enddo

      close(iu1)
      close(iu2)
      close(iu3)
      close(iu4)
      close(iu5)
      close(iu6)
      

      end subroutine offregrid
c****



      subroutine parallel_read_regrid(iunit,name,nskip,tcubloc,ipos)
c     
c     Read input data on lat-lon grid, regrid to global cubbed sphere grid
c     then scatter to all subdomains
c
      use regrid_com

      use mpp_mod
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use gatscat_mod

      implicit none
      integer, intent(in) :: iunit
      character*16, intent(in) :: name
      integer, intent(in) :: nskip

      real*8, intent(out) :: tcubloc(is:ie,js:ie) !output regridded and scattered data 
      real*4 :: tllr4(IM,JM)  ! latlon real*4 data read from input file
      real*8 :: tdatall(IM,JM)  ! latlon real*8 data 
      real*4 :: X              !dummy arrays
      real*8 :: tcubglob(ic,jc)
      integer, intent(in) :: ipos
      integer :: n,ierr
      character*80 :: TITLE     

c need modele.f to call init_regrid_rootpe

      if (gid .eq. 0) then
         do n=1,ipos-1
            read(iunit,IOSTAT=ierr)
         enddo
         read(iunit,IOSTAT=ierr) TITLE, (X,n=1,nskip), tllr4
c     convert from real*4 to real*8
         tdatall=tllr4

c
c     Regrid from lat-lon to cubbed sphere, form global array
c
         call root_regrid_ll2cs(tdatall,tcubglob)

      endif

c
c     Scatter data to every processor
c
      call unpack_data(tcubglob,tcubloc)

      end subroutine parallel_read_regrid


