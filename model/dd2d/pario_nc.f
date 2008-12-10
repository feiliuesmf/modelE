      module pario
      use dd2d_utils, only : dd2d_grid,pack_row,unpack_row,get_nlnk
      implicit none
      save
      private
      include 'mpif.h'
      include 'netcdf.inc'

c
c i/o interfaces
c
      public :: write_dist_data,read_dist_data
      interface write_dist_data
        module procedure par_write_nc_2D
        module procedure par_write_nc_3D
        module procedure par_write_nc_4D
        module procedure par_write_nc_5D
        module procedure par_write_nc_2D_int
        module procedure par_write_nc_3D_int
        module procedure par_write_nc_4D_int
        module procedure par_write_nc_2D_logical
      end interface write_dist_data
      interface read_dist_data
        module procedure par_read_nc_2D
        module procedure par_read_nc_3D
        module procedure par_read_nc_4D
        module procedure par_read_nc_5D
        module procedure par_read_nc_2D_int
        module procedure par_read_nc_3D_int
        module procedure par_read_nc_4D_int
        module procedure par_read_nc_2D_logical
      end interface read_dist_data

      public :: write_data,read_data
      interface write_data
        module procedure write_nc_0D
        module procedure write_nc_1D
        module procedure write_nc_2D
        module procedure write_nc_3D
        module procedure write_nc_4D
        module procedure write_nc_5D
        module procedure write_nc_0D_int
        module procedure write_nc_1D_int
        module procedure write_nc_2D_int
        module procedure write_nc_3D_int
      end interface write_data
      interface read_data
        module procedure read_nc_0D
        module procedure read_nc_1D
        module procedure read_nc_2D
        module procedure read_nc_3D
        module procedure read_nc_4D
        module procedure read_nc_5D
        module procedure read_nc_0D_int
        module procedure read_nc_1D_int
        module procedure read_nc_2D_int
        module procedure read_nc_3D_int
      end interface read_data

      public :: defvar
      interface defvar
        module procedure defvar_0D
        module procedure defvar_1D
        module procedure defvar_2D
        module procedure defvar_3D
        module procedure defvar_4D
        module procedure defvar_5D
        module procedure defvar_0D_int
        module procedure defvar_1D_int
        module procedure defvar_2D_int
        module procedure defvar_3D_int
        module procedure defvar_4D_int
        module procedure defvar_5D_int
        module procedure defvar_2D_logical
      end interface

      contains

      subroutine par_write_nc_2D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_2D
      subroutine par_write_nc_3D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_3D
      subroutine par_write_nc_4D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_4D
      subroutine par_write_nc_5D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_5D

      subroutine par_read_nc_2D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_2D
      subroutine par_read_nc_3D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_3D
      subroutine par_read_nc_4D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_4D
      subroutine par_read_nc_5D(grid,fid,varname,arr,jdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_5D

      subroutine par_read_nc_2D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      call read_dist_data(grid,fid,varname,arr)
      iarr = arr
      end subroutine par_read_nc_2D_int
      subroutine par_read_nc_3D_int(grid,fid,varname,iarr,jdim)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      if(present(jdim)) then
        call read_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call read_dist_data(grid,fid,varname,arr)
      endif
      iarr = arr
      end subroutine par_read_nc_3D_int
      subroutine par_read_nc_4D_int(grid,fid,varname,iarr,jdim)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3),size(iarr,4))
      if(present(jdim)) then
        call read_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call read_dist_data(grid,fid,varname,arr)
      endif
      iarr = arr
      end subroutine par_read_nc_4D_int

      subroutine par_write_nc_2D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      arr = iarr
      call write_dist_data(grid,fid,varname,arr)
      end subroutine par_write_nc_2D_int
      subroutine par_write_nc_3D_int(grid,fid,varname,iarr,jdim)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      arr = iarr
      if(present(jdim)) then
        call write_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call write_dist_data(grid,fid,varname,arr)
      endif
      end subroutine par_write_nc_3D_int
      subroutine par_write_nc_4D_int(grid,fid,varname,iarr,jdim)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:,:,:)
      integer, intent(in), optional :: jdim
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3),size(iarr,4))
      arr = iarr
      if(present(jdim)) then
        call write_dist_data(grid,fid,varname,arr,jdim=jdim)
      else
        call write_dist_data(grid,fid,varname,arr)
      endif
      end subroutine par_write_nc_4D_int

      subroutine par_read_nc_2D_logical(grid,fid,varname,larr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      call read_dist_data(grid,fid,varname,arr)
      larr = arr.eq.1d0
      end subroutine par_read_nc_2D_logical
      subroutine par_write_nc_2D_logical(grid,fid,varname,larr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      where(larr)
        arr = 1d0
      else where
        arr = 0d0
      end where
      call write_dist_data(grid,fid,varname,arr)
      end subroutine par_write_nc_2D_logical

      subroutine stoprc(rc,rc_ok)
      integer :: rc,rc_ok
      integer :: mpi_err
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      if(rc.ne.rc_ok) call mpi_abort(MPI_COMM_WORLD,1,mpi_err)
      return
      end subroutine stoprc

      subroutine write_nc_0D(grid,fid,varname,arr)
      real*8 :: arr
#include "do_write_nc.inc"
      end subroutine write_nc_0D
      subroutine write_nc_1D(grid,fid,varname,arr)
      real*8 :: arr(:)
#include "do_write_nc.inc"
      end subroutine write_nc_1D
      subroutine write_nc_2D(grid,fid,varname,arr)
      real*8 :: arr(:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_2D
      subroutine write_nc_3D(grid,fid,varname,arr)
      real*8 :: arr(:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_3D
      subroutine write_nc_4D(grid,fid,varname,arr)
      real*8 :: arr(:,:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_4D
      subroutine write_nc_5D(grid,fid,varname,arr)
      real*8 :: arr(:,:,:,:,:)
#include "do_write_nc.inc"
      end subroutine write_nc_5D

      subroutine read_nc_0D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr
#include "do_read_nc.inc"
      end subroutine read_nc_0D
      subroutine read_nc_1D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:)
#include "do_read_nc.inc"
      end subroutine read_nc_1D
      subroutine read_nc_2D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_2D
      subroutine read_nc_3D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_3D
      subroutine read_nc_4D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_4D
      subroutine read_nc_5D(grid,fid,varname,arr,bcast_all)
      real*8 :: arr(:,:,:,:,:)
#include "do_read_nc.inc"
      end subroutine read_nc_5D

      subroutine write_nc_0D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr
      real*8 :: arr
      arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_0D_int
      subroutine write_nc_1D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:)
      real*8 :: arr(size(iarr))
      arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_1D_int
      subroutine write_nc_2D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_2D_int
      subroutine write_nc_3D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_3D_int

      subroutine read_nc_0D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr
      logical, intent(in), optional :: bcast_all
      real*8 :: arr
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(grid%am_i_globalroot .or. bc_all) iarr = arr
      end subroutine read_nc_0D_int
      subroutine read_nc_1D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(iarr))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(grid%am_i_globalroot .or. bc_all) iarr = arr
      end subroutine read_nc_1D_int
      subroutine read_nc_2D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(iarr,1),size(iarr,2))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(grid%am_i_globalroot .or. bc_all) iarr = arr
      end subroutine read_nc_2D_int
      subroutine read_nc_3D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(grid%am_i_globalroot .or. bc_all) iarr = arr
      end subroutine read_nc_3D_int

      subroutine defvar_0D(grid,fid,arr,varinfo)
      real*8 :: arr
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D
      subroutine defvar_1D(grid,fid,arr,varinfo)
      real*8 :: arr(:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D
      subroutine defvar_2D(grid,fid,arr,varinfo)
      real*8 :: arr(:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D
      subroutine defvar_3D(grid,fid,arr,varinfo)
      real*8 :: arr(:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D
      subroutine defvar_4D(grid,fid,arr,varinfo)
      real*8 :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D
      subroutine defvar_5D(grid,fid,arr,varinfo)
      real*8 :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D

      subroutine defvar_0D_int(grid,fid,arr,varinfo)
      integer :: arr
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D_int
      subroutine defvar_1D_int(grid,fid,arr,varinfo)
      integer :: arr(:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D_int
      subroutine defvar_2D_int(grid,fid,arr,varinfo)
      integer :: arr(:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_int
      subroutine defvar_3D_int(grid,fid,arr,varinfo)
      integer :: arr(:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D_int
      subroutine defvar_4D_int(grid,fid,arr,varinfo)
      integer :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D_int
      subroutine defvar_5D_int(grid,fid,arr,varinfo)
      integer :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D_int

      subroutine defvar_2D_logical(grid,fid,arr,varinfo)
      logical :: arr(:,:)
c netcdf file will represent logical as 0/1 int
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_logical


      subroutine define_var(fid,dtype,
     &     varinfo_in,ndims_in,shp_in,im_in,jm_in,
     &     ntiles,rc)
      integer :: fid,dtype,rc
      character(len=*) :: varinfo_in
      integer :: ndims_in,shp_in(ndims_in),im_in,jm_in,ntiles
      character(len=40) :: vname,dname
      character(len=80) :: varinfo
      character*1, dimension(:), allocatable :: char_arr
      integer :: i,ndims,l,l1,l2,xdim,lv,vid,dsize,status
     &     ,dsizx
      integer :: dids(7)
      logical :: is_dist
      rc = 0
      varinfo=trim(varinfo_in)
c remove spaces
      varinfo=''
      l=0
      do i=1,len_trim(varinfo_in)
        if(varinfo_in(i:i).eq.' ') cycle
        l=l+1
        varinfo(l:l)=varinfo_in(i:i)
      enddo
      lv=l
      allocate(char_arr(len_trim(varinfo)))
      do i=1,len_trim(varinfo)
        char_arr(i)=varinfo(i:i)
      enddo
      l=count(char_arr.eq.'('.or.char_arr.eq.')')
      if(l.eq.0) then
        vname=trim(varinfo)
        ndims=0
      elseif(l.eq.2) then
        vname=varinfo(1:index(varinfo,'(')-1)
        ndims=1+count(char_arr.eq.',')
      endif
      deallocate(char_arr)
      if(l.ne.0 .and. l.ne.2) then
        write(6,*) 'parsing error: ',trim(varinfo)
        rc = 1; return
      endif
      if(ndims.ne.ndims_in) then
        write(6,*) 'parsed ndims does not match actual ndims for ',
     &       trim(vname)
        rc = 1; return
      endif
c check whether variable is already defined
      if(nf_inq_varid(fid,trim(vname),vid).eq.nf_noerr) then
        write(6,*) 'error: variable ',trim(vname),
     &       ' is already defined'
        rc = 1; return
      endif
c
c loop through dimensions and define them if necessary
c
      is_dist = .false.
      if(ndims.gt.0) then
        l1=index(varinfo,'(')+1
        do xdim=1,ndims
          if(xdim.lt.ndims) then
            l2=l1+index(varinfo(l1:lv),',')-2
          else
            l2=l1+index(varinfo(l1:lv),')')-2
          endif
          dname=varinfo(l1:l2)
          dsize = shp_in(xdim)
          write(*,*) vname," ",dname," ",dsize
          if(dname(1:6).eq.'dist_i') then
            dname=dname(6:len_trim(dname))
            dsize = im_in
            is_dist = .true.
          elseif(dname(1:6).eq.'dist_j') then
            dname=dname(6:len_trim(dname))
            dsize = jm_in
            is_dist = .true.
          endif
          if(nf_inq_dimid(fid,trim(dname),dids(xdim)).eq.nf_noerr) then
            status = nf_inq_dimlen(fid,dids(xdim),dsizx)
            if(dsizx.ne.dsize) then
              write(6,*) 'illegal operation: changing dimension size ',
     &             trim(dname),' for variable ',trim(vname)
              rc = 1; return
            endif
          else
            status = nf_def_dim(fid,trim(dname),dsize,dids(xdim))
          endif
          l1=l2+2
        enddo
      endif
c
c if more than one tile, add an extra dimension for distributed vars
c
      if(ntiles.gt.1 .and. is_dist) then
        ndims = ndims + 1
        xdim = ndims
        dname='tile'
        dsize = ntiles
        if(nf_inq_dimid(fid,trim(dname),dids(xdim)).eq.nf_noerr) then
          status = nf_inq_dimlen(fid,dids(xdim),dsizx)
          if(dsizx.ne.dsize) then
            write(6,*) 'illegal operation: changing dimension size ',
     &           trim(dname),' for variable ',trim(vname)
            rc = 1; return
          endif
        else
          status = nf_def_dim(fid,trim(dname),dsize,dids(xdim))
        endif
      endif
c
c define the variable
c
      status = nf_def_var(fid,trim(vname),dtype,ndims,dids,vid)
      return
      end subroutine define_var

      end module pario
