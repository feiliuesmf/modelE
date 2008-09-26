      program regrid
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
c      call mp_stop()

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
      call parallel_regrid(tcub_loc,tlatlon,alatlon)   


c
c     sum all contributions
c
      nij_latlon=im*jm
      call mpp_sum(tlatlon,nij_latlon)
      call mpp_sum(alatlon,nij_latlon)
      
c
c     root proc section
c
      if (gid .eq. rootpe) then
      tlatlon(:,:) = tlatlon(:,:)/alatlon(:,:)

c      write(6,*) tlatlon(:,:)

      ofi='tstout.nc'
      status = nf_open(trim(ofi),nf_write,fid)
      if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
      status = nf_inq_varid(fid,'lwup_sfc',vid)
      write(*,*) NF_STRERROR(status)
      status = nf_put_var_double(fid,vid,tlatlon)
      write(*,*) NF_STRERROR(status)
 0    
      status = nf_close(fid)
      endif

      end program regrid



      subroutine parallel_regrid(tcub_loc,tlatlon,alatlon)
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid
      implicit none
      real*8 :: tlatlon(im,jm),alatlon(im,jm),
     &     tcub_loc(isd:ied,jsd:jed)
      integer :: loctile  !tile index of current domain
      integer :: dom_per_tile,n,ic,jc,i,j,itile
      
      dom_per_tile=4            ! change this
      
      loctile=gid/dom_per_tile
      loctile=loctile+1

c      write(6,*) 'GIIIDD=',gid
c      write(6,*) 'LOCTILE=',loctile

      alatlon(:,:) = 0d0
      tlatlon(:,:) = 0d0

      do n=1,ncells
        ic=ijcub(1,n)
        jc=ijcub(2,n)
        i=ijlatlon(1,n)
        j=ijlatlon(2,n)
        itile=tile(n)
c        write(6,*) 'ie, ic, is=',ie,ic,is
c        write(6,*) 'je, jc, js=',je,jc,js

        if (itile .eq. loctile) then
           if ( ( (ic .le. ie) .and. (ic .ge. is)) .and.
     &          ( (jc .le. je) .and. (jc .ge. js)) ) then
              write(6,*) 'xg_area',xgrid_area(n)
              alatlon(i,j) = alatlon(i,j) + xgrid_area(n)
              tlatlon(i,j) = tlatlon(i,j) + xgrid_area(n)
     &             *tcub_loc(ic,jc)
           endif
        endif
      enddo
      end subroutine parallel_regrid



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
      
      integer :: status,fid,n,vid,ikey,jlat,maxkey,rootpe
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
         write(icch,'(i10)') icc
         write(jcch,'(i10)') jcc
         
         imch=trim(adjustl(imch))
         jmch=trim(adjustl(jmch))
         icch=trim(adjustl(icch))
         jcch=trim(adjustl(jcch))
         exchfile="exexch.nc"
c         exchfile="remap"//trim(imch)//"-"//trim(jmch)
c     *        //"C"//trim(icch)//"-"//trim(jcch)//".nc"
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

