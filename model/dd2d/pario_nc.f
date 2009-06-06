      module pario

#ifdef OFFLINE
#else
c see whether model E is running in serial mode
#ifndef USE_ESMF
#define SERIAL_MODE
#endif
#endif

      use dd2d_utils, only : dist_grid
#ifndef SERIAL_MODE
c these routines are only needed when running on multiple CPUs
      use dd2d_utils, only : pack_row,unpack_row,get_nlnk,pack_data
#endif
      implicit none
      save
      private

#ifndef SERIAL_MODE
      include 'mpif.h'
#endif
      include 'netcdf.inc'

c
c i/o interfaces
c
      public :: par_open,par_close

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
        module procedure write_nc_2D_logical
        module procedure write_nc_1D_array_of_strings
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
        module procedure read_nc_2D_logical
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
        module procedure defvar_1D_array_of_strings
      end interface

      public :: write_attr
      interface write_attr
        module procedure write_attr_text
        module procedure write_attr_0D_r8
        module procedure write_attr_1D_r8
        module procedure write_attr_0D_int
        module procedure write_attr_1D_int
      end interface

      public :: read_attr
      interface read_attr
        module procedure read_attr_text
        module procedure read_attr_0D_r8
        module procedure read_attr_1D_r8
        module procedure read_attr_0D_int
        module procedure read_attr_1D_int
      end interface

      interface len_of_obj
        module procedure len_of_text
        module procedure len_of_int0D
        module procedure len_of_int1D
        module procedure len_of_r80D
        module procedure len_of_r81D
      end interface

      interface broadcast
        module procedure broadcast_0D_int
        module procedure broadcast_1D_int
        module procedure broadcast_0D_r8
        module procedure broadcast_1D_r8
      end interface broadcast

#ifndef SERIAL_MODE
c these routines are only needed when running on multiple CPUs
      interface pack_row_no_xdim
        module procedure pack_row_no_xdim_2d
        module procedure pack_row_no_xdim_3d
        module procedure pack_row_no_xdim_4d
        module procedure pack_row_no_xdim_5d
      end interface
      interface unpack_row_no_xdim
        module procedure unpack_row_no_xdim_2d
        module procedure unpack_row_no_xdim_3d
        module procedure unpack_row_no_xdim_4d
        module procedure unpack_row_no_xdim_5d
      end interface

      interface par_write_jdecomp_optimized
        module procedure par_write_ij
        module procedure par_write_ijx
        module procedure par_write_ijxx
        module procedure par_write_ijxxx
      end interface
#endif /* not SERIAL_MODE */

      integer, parameter :: success = 0, fail = -1

      contains

      function par_open(grid,fname,mode)
      type(dist_grid), intent(in) :: grid
      character(len=*) :: fname
      character(len=*) :: mode
      integer :: par_open
      integer :: rc,rc2,fid,vid,wc,idum
      if(grid%am_i_globalroot) then
        if(trim(mode).eq.'create') then
c          rc = nf_create(trim(fname),nf_clobber,fid)
          rc = nf_create(trim(fname),nf_64bit_offset,fid) ! when files get big
          if(rc.ne.nf_noerr) write(6,*)
     &         'error creating ',trim(fname)
        elseif(trim(mode).eq.'write') then
          rc = nf_open(trim(fname),nf_write,fid)
          if(rc.ne.nf_noerr) write(6,*)
     &         'error opening ',trim(fname)
        elseif(trim(mode).eq.'read') then
          rc = nf_open(trim(fname),nf_nowrite,fid)
          if(rc.ne.nf_noerr) then
            write(6,*) 'error opening ',trim(fname)
          else
            rc2 = nf_inq_varid(fid,'write_status',vid)
            if(rc2.eq.nf_noerr) then
              rc2 = nf_get_var_int(fid,vid,wc)
              if(wc.ne.success) then
                write(6,*) 'input file ',trim(fname),
     &            ' does not appear to have been written successfully:'
                write(6,*) 'write_status = ',wc
              endif
            else
              wc = success
            endif
          endif
        else
          write(6,*) 'par_open: invalid mode ',trim(mode)
          write(6,*) 'mode must be one of [create write read]'
          rc = nf_noerr + 1
        endif
      endif
      call stoprc(rc,nf_noerr)
      if(trim(mode).eq.'read') call stoprc(wc,success)
c define/overwrite the success flag for error checking
      if(grid%am_i_globalroot) then
        if(trim(mode).eq.'create') then
          rc = nf_def_var(fid,'write_status',nf_int,0,idum,vid)
          rc = nf_enddef(fid)
          rc = nf_put_var_int(fid,vid,fail)
          rc = nf_redef(fid)
        elseif(trim(mode).eq.'write') then
          rc = nf_inq_varid(fid,'write_status',vid)
          rc = nf_put_var_int(fid,vid,fail)
          rc = nf_sync(fid)
        endif
        par_open = fid
      else
        par_open = -1
      endif
      return
      end function par_open

      subroutine par_close(grid,fid)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      integer :: rc,vid
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,'write_status',vid)
        if(rc.eq.nf_noerr) then
          rc = nf_put_var_int(fid,vid,success)
        endif
        rc = nf_close(fid)
      endif
      return
      end subroutine par_close

      subroutine par_write_nc_2D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_2D
      subroutine par_write_nc_3D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_3D
      subroutine par_write_nc_4D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_4D
      subroutine par_write_nc_5D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_write_nc.inc"
      end subroutine par_write_nc_5D

      subroutine par_read_nc_2D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_2D
      subroutine par_read_nc_3D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_3D
      subroutine par_read_nc_4D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_4D
      subroutine par_read_nc_5D(grid,fid,varname,arr,jdim,no_xdim)
      real*8 :: arr(:,:,:,:,:)
#include "do_par_read_nc.inc"
      end subroutine par_read_nc_5D

      subroutine par_read_nc_2D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      call read_dist_data(grid,fid,varname,arr)
      iarr = arr
      end subroutine par_read_nc_2D_int
      subroutine par_read_nc_3D_int(grid,fid,varname,iarr,jdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      arr = iarr
      call write_dist_data(grid,fid,varname,arr)
      end subroutine par_write_nc_2D_int
      subroutine par_write_nc_3D_int(grid,fid,varname,iarr,jdim)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      call read_dist_data(grid,fid,varname,arr)
      larr = arr.eq.1d0
      end subroutine par_read_nc_2D_logical
      subroutine par_write_nc_2D_logical(grid,fid,varname,larr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
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
#ifndef SERIAL_MODE
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      if(rc.ne.rc_ok) call mpi_abort(MPI_COMM_WORLD,1,mpi_err)
#else
      if(rc.ne.rc_ok) stop
#endif
      return
      end subroutine stoprc

      subroutine broadcast_0D_int(i)
      integer :: i
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_0D_int
      subroutine broadcast_1D_int(i)
      integer :: i(:)
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(i,size(i),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_1D_int
      subroutine broadcast_0D_r8(r8)
      real*8 :: r8
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(r8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_0D_r8
      subroutine broadcast_1D_r8(r8)
      real*8 :: r8(:)
#ifndef SERIAL_MODE
      integer :: ierr
      call mpi_bcast(r8,size(r8),MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_WORLD,ierr)
#endif
      end subroutine broadcast_1D_r8

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
      type(dist_grid), intent(in) :: grid
      integer :: iarr
      real*8 :: arr
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_0D_int
      subroutine write_nc_1D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:)
      real*8 :: arr(size(iarr))
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_1D_int
      subroutine write_nc_2D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2))
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_2D_int
      subroutine write_nc_3D_int(grid,fid,varname,iarr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer :: iarr(:,:,:)
      real*8 :: arr(size(iarr,1),size(iarr,2),size(iarr,3))
      if(grid%am_i_globalroot) arr = iarr
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_3D_int
      subroutine write_nc_2D_logical(grid,fid,varname,larr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      real*8 :: arr(size(larr,1),size(larr,2))
      if(grid%am_i_globalroot) then
        where(larr)
          arr = 1d0
        else where
          arr = 0d0
        end where
      endif
      call write_data(grid,fid,varname,arr)
      end subroutine write_nc_2D_logical
      subroutine write_nc_1D_array_of_strings(grid,fid,varname,arr)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      character(len=*) :: arr(:)
      integer :: rc,vid
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      if(grid%am_i_globalroot) then
        rc = nf_put_var_text(fid,vid,arr)
      endif
      end subroutine write_nc_1D_array_of_strings

      subroutine read_nc_0D_int(grid,fid,varname,iarr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
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
      type(dist_grid), intent(in) :: grid
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
      subroutine read_nc_2D_logical(grid,fid,varname,larr,bcast_all)
      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      logical :: larr(:,:)
      logical, intent(in), optional :: bcast_all
      real*8 :: arr(size(larr,1),size(larr,2))
      logical :: bc_all
      bc_all=.false.
      if(present(bcast_all)) then
        if(bcast_all) bc_all=.true.
      endif
      call read_data(grid,fid,varname,arr,bcast_all=bc_all)
      if(grid%am_i_globalroot .or. bc_all) larr = arr.eq.1d0
      end subroutine read_nc_2D_logical

      subroutine defvar_0D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D
      subroutine defvar_1D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D
      subroutine defvar_2D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D
      subroutine defvar_3D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D
      subroutine defvar_4D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D
      subroutine defvar_5D(grid,fid,arr,varinfo,r4_on_disk,defby)
      real*8 :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_double
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D

      subroutine defvar_0D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_0D_int
      subroutine defvar_1D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_1D_int
      subroutine defvar_2D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_int
      subroutine defvar_3D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_3D_int
      subroutine defvar_4D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_4D_int
      subroutine defvar_5D_int(grid,fid,arr,varinfo,r4_on_disk,defby)
      integer :: arr(:,:,:,:,:)
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_5D_int

      subroutine defvar_2D_logical(grid,fid,arr,varinfo,r4_on_disk,
     &     defby)
      logical :: arr(:,:)
c netcdf file will represent logical as 0/1 int
      integer, parameter :: dtype=nf_int
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_logical

      subroutine defvar_2D_char(grid,fid,arr,varinfo,r4_on_disk,
     &     defby)
      character :: arr(:,:)
      integer, parameter :: dtype=nf_char
      include 'do_defvar_nc.inc'
      return
      end subroutine defvar_2D_char

      subroutine defvar_1D_array_of_strings(grid,fid,arr,varinfo)
      type(dist_grid), intent(in) :: grid
      integer :: fid
      character(len=*) :: arr(:)
      character(len=*) :: varinfo
      character, allocatable :: arr2d(:,:)
      allocate(arr2d(len(arr(1)),size(arr,1)))
      call defvar_2D_char(grid,fid,arr2d,varinfo)
      deallocate(arr2d)
      return
      end subroutine defvar_1D_array_of_strings

      subroutine write_attr_text(grid,fid,varname,attname,attval)
      character(len=*) :: attval
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_text(fid,vid,trim(attname),attlen,attval)
        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_text
      subroutine write_attr_0D_int(grid,fid,varname,attname,attval)
      integer :: attval
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_int(fid,vid,trim(attname),nf_int,attlen,attval)
        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_0D_int
      subroutine write_attr_1D_int(grid,fid,varname,attname,attval)
      integer :: attval(:)
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_int(fid,vid,trim(attname),nf_int,attlen,attval)
        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_1D_int
      subroutine write_attr_0D_r8(grid,fid,varname,attname,attval)
      real*8 :: attval
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_double(fid,vid,trim(attname),nf_double,attlen,
     &       attval)
        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_0D_r8
      subroutine write_attr_1D_r8(grid,fid,varname,attname,attval)
      real*8 :: attval(:)
#include "setup_attput.inc"
      if(grid%am_i_globalroot) then
        rc = nf_put_att_double(fid,vid,trim(attname),nf_double,attlen,
     &       attval)
        if(do_enddef) rc2 = nf_enddef(fid)
      endif
      call stoprc(rc,nf_noerr)
      return
      end subroutine write_attr_1D_r8

      subroutine read_attr_text(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      character(len=*) :: attval
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_text(fid,vid,trim(attname),tmpstr)
      endif
      call stoprc(rc,nf_noerr)
#ifndef SERIAL_MODE
      call mpi_bcast(tmpstr,attlen,MPI_CHARACTER,0,
     &     MPI_COMM_WORLD,ierr)
#endif
      attval=''
      do l=1,attlen
        attval(l:l) = tmpstr(l)
      enddo
      return
      end subroutine read_attr_text
      subroutine read_attr_0D_int(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      integer :: attval
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_int(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_0D_int
      subroutine read_attr_1D_int(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      integer :: attval(:)
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_int(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_1D_int
      subroutine read_attr_0D_r8(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      real*8 :: attval
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_double(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_0D_r8
      subroutine read_attr_1D_r8(grid,fid,varname,attname,attlen,
     &         attval,attnum)
      real*8 :: attval(:)
#include "setup_attget.inc"
      if(grid%am_i_globalroot) then
        rc = nf_get_att_double(fid,vid,trim(attname),attval)
      endif
      call stoprc(rc,nf_noerr)
      call broadcast(attval)
      return
      end subroutine read_attr_1D_r8

      function len_of_text(cstr)
      integer :: len_of_text
      character(len=*) :: cstr
      len_of_text = len_trim(cstr)
      return
      end function len_of_text
      function len_of_int0D(i)
      integer :: len_of_int0D
      integer :: i
      len_of_int0D = 1
      return
      end function len_of_int0D
      function len_of_int1D(i)
      integer :: len_of_int1D
      integer :: i(:)
      len_of_int1D = size(i)
      return
      end function len_of_int1D
      function len_of_r80D(r8)
      integer :: len_of_r80D
      real*8 :: r8
      len_of_r80D = 1
      return
      end function len_of_r80D
      function len_of_r81D(r8)
      integer :: len_of_r81D
      real*8 :: r8(:)
      len_of_r81D = size(r8)
      return
      end function len_of_r81D

      subroutine define_var(fid,dtype,
     &     varinfo_in,ndims_in,shp_in,im_in,jm_in,
     &     ntiles,rc,vid)
      integer :: fid,dtype,rc,vid
      character(len=*) :: varinfo_in
      integer :: ndims_in,shp_in(ndims_in),im_in,jm_in,ntiles
      character(len=40) :: vname,dname
      character(len=80) :: varinfo
      character*1, dimension(:), allocatable :: char_arr
      integer :: i,ndims,l,l1,l2,xdim,lv,dsize,status
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

#ifndef SERIAL_MODE
      subroutine pack_row_no_xdim_2d(grid,arr,arr1d,jdim)
      real*8 arr(:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_2d
      subroutine pack_row_no_xdim_3d(grid,arr,arr1d,jdim)
      real*8 arr(:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_3d
      subroutine pack_row_no_xdim_4d(grid,arr,arr1d,jdim)
      real*8 arr(:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_4d
      subroutine pack_row_no_xdim_5d(grid,arr,arr1d,jdim)
      real*8 arr(:,:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      return
      end subroutine pack_row_no_xdim_5d

      subroutine unpack_row_no_xdim_2d(grid,arr1d,arr,jdim)
      real*8 arr(:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_2d
      subroutine unpack_row_no_xdim_3d(grid,arr1d,arr,jdim)
      real*8 arr(:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_3d
      subroutine unpack_row_no_xdim_4d(grid,arr1d,arr,jdim)
      real*8 arr(:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_4d
      subroutine unpack_row_no_xdim_5d(grid,arr1d,arr,jdim)
      real*8 arr(:,:,:,:,:)
      include 'row_setup_no_xdim.inc'
      call copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      return
      end subroutine unpack_row_no_xdim_5d

      subroutine par_write_ij(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc
      if(grid%am_i_globalroot) allocate(arrgij(grid%npx,grid%npy))
      call pack_data(grid,arr,arrgij)
      if(grid%am_i_globalroot) then
        rc = nf_put_var_double(fid,vid,arrgij)
        deallocate(arrgij)
      endif
      return
      end subroutine par_write_ij
      subroutine par_write_ijx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc,k,srt(3),cnt(3)
      if(grid%am_i_globalroot) allocate(arrgij(grid%npx,grid%npy))
      srt(1:2) = 1; cnt(1:3) = (/ grid%npx, grid%npy, 1 /)
      do k=1,size(arr,3)
        call pack_data(grid,arr(:,:,k),arrgij)
        srt(3) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgij)
      enddo
      if(grid%am_i_globalroot) deallocate(arrgij)
      return
      end subroutine par_write_ijx
      subroutine par_write_ijxx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc,k,l,srt(4),cnt(4)
      if(jdim.eq.3) then
        call par_write_xijx(grid,fid,vid,arr,jdim)
        return
      endif
      if(grid%am_i_globalroot) allocate(arrgij(grid%npx,grid%npy))
      srt(1:2) = 1; cnt(1:4) = (/ grid%npx, grid%npy, 1, 1 /)
      do l=1,size(arr,4)
      srt(4) = l
      do k=1,size(arr,3)
        call pack_data(grid,arr(:,:,k,l),arrgij)
        srt(3) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgij)
      enddo
      enddo
      if(grid%am_i_globalroot) deallocate(arrgij)
      return
      end subroutine par_write_ijxx
      subroutine par_write_ijxxx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgij(:,:)
      integer :: rc,k,l,m,srt(5),cnt(5)
      if(grid%am_i_globalroot) allocate(arrgij(grid%npx,grid%npy))
      srt(1:2) = 1; cnt(1:5) = (/ grid%npx, grid%npy, 1, 1, 1 /)
      do m=1,size(arr,5)
      srt(5) = m
      do l=1,size(arr,4)
      srt(4) = l
      do k=1,size(arr,3)
        call pack_data(grid,arr(:,:,k,l,m),arrgij)
        srt(3) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgij)
      enddo
      enddo
      enddo
      if(grid%am_i_globalroot) deallocate(arrgij)
      return
      end subroutine par_write_ijxxx
      subroutine par_write_xijx(grid,fid,vid,arr,jdim)
      type(dist_grid), intent(in) :: grid
      real*8 :: arr(:,:,:,:)
      integer :: fid,vid,jdim
      real*8, allocatable :: arrgxij(:,:,:)
      integer :: rc,k,srt(4),cnt(4)
      if(grid%am_i_globalroot)
     &     allocate(arrgxij(size(arr,1),grid%npx,grid%npy))
      srt(1:3) = 1; cnt(1:4) = (/ size(arr,1), grid%npx, grid%npy, 1 /)
      do k=1,size(arr,4)
        call pack_data(grid,arr(:,:,:,k),arrgxij,jdim=3)
        srt(4) = k
        if(grid%am_i_globalroot)
     &       rc = nf_put_vara_double(fid,vid,srt,cnt,arrgxij)
      enddo
      if(grid%am_i_globalroot) deallocate(arrgxij)
      return
      end subroutine par_write_xijx
#endif /* not SERIAL_MODE */

      end module pario

      subroutine copy_to_1D_no_xdim(arr,arr1d,nl,nj,nk,j1,j2)
      implicit none
      real*8 arr(nl,nj,nk)
      real*8 arr1d(1)
      integer :: nl,nj,nk,j1,j2
      integer :: j,k,l,n
      n = 0
      do k=1,nk
      do j=j1,j2
      do l=1,nl
        n = n + 1
        arr1d(n) = arr(l,j,k)
      enddo
      enddo
      enddo
      return
      end subroutine copy_to_1D_no_xdim
      subroutine copy_from_1D_no_xdim(arr1d,arr,nl,nj,nk,j1,j2)
      implicit none
      real*8 arr(nl,nj,nk)
      real*8 arr1d(1)
      integer :: nl,nj,nk,j1,j2
      integer :: j,k,l,n
      n = 0
      do k=1,nk
      do j=j1,j2
      do l=1,nl
        n = n + 1
        arr(l,j,k) = arr1d(n)
      enddo
      enddo
      enddo
      return
      end subroutine copy_from_1D_no_xdim
