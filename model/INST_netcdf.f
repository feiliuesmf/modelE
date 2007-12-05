c The following routines in this file can be used to output a model
c variable VAR whose dimensions match those of the variable VARNAME
c in a predefined netcdf file nctemp.nc in the run directory.
c VARNAME is assumed to be a 30-byte character string.
c
c       subroutine inst_ncoutd(VARNAME,VAR) ! arbitrary REAL*8 VAR
c       subroutine inst_ncouti(VARNAME,VAR) ! arbitrary INTEGER VAR
c ESMF-parallelized versions:
c       subroutine inst_ncout_ij(VARNAME,VAR) ! lon-lat VAR
c       subroutine inst_ncout_ijl(VARNAME,VAR) ! lon-lat-level VAR
c       subroutine inst_ncout_lij(VARNAME,VAR) ! level-lon-lat VAR
c
c The module routines inst_ncout_open/inst_ncout_close must be called
c at the beginning/end of model execution (on the root processor
c if in multiprocessor mode).  They take no arguments.
c      use inst_ncout, only : inst_ncout_open,inst_ncout_close
c
c The routines inst_ncoutd and inst_ncouti can output a VAR of arbitrary
c shape and size from a single processor.  To output data distributed
c over multiple processors using the ESMF infrastructure, the ij/ijl/lij
c routines must be used (they first gather the data to the root processor).
c Overloading has not yet been implemented to allow ij/ijl/lij 
c arrays to all be handled by a single routine.

      module inst_ncout
!@sum  inst_ncout facilitates the writing of model variables in netcdf format
!@auth M. Kelley
      use MODEL_COM, only : im,jm,lm
      implicit none

      include 'netcdf.inc'

      private
      public ::
     &     inst_ncout_open,inst_ncout_close,inst_setup_var
     &    ,varid,status,fid,srt,cnt
     &    ,xij_glob,xijl_glob,xlij_glob

      real*8, dimension(im,jm) :: xij_glob
      real*8, dimension(im,jm,lm) :: xijl_glob
      real*8, dimension(lm,im,jm) :: xlij_glob

      character(len=80) :: outfile
!@var status return code from netcdf calls
!@var fid ID-number of output file
!@var ndims number of dimensions of current output variable
!@var varid ID-number of current output variable
      integer :: status,fid,ndims,varid,unlimdimid
c maximum of 100 variables in output file, 7 dimensions for each
      integer, dimension(100) :: varid_list,ndims_list,counts_list
      integer, dimension(7,100) :: dimlens_list,dimids_list
      logical, dimension(100) :: write_all_list
!@var dimids dimension ID-numbers of output variable
!@var dimlens dimension lengths of output variable
!@var srt start indices for output from current output variable
!@var cnt number of elements to output from current output variable
      integer, dimension(7) :: dimids,dimlens,srt,cnt
      integer :: nvar

      contains

      subroutine inst_ncout_open
      outfile='nctemp.nc'
      status = nf_open (trim(outfile), nf_write, fid)
      if(status.ne.nf_noerr) then
         write(6,*) 'cannot open: '//trim(outfile)
         call stop_model('stopped in INST_netcdf.f',255)
      endif
      nvar = 0
      counts_list(:) = 0
      write_all_list(:) = .false.
      status = nf_inq_unlimdim(fid,unlimdimid)
      return
      end subroutine inst_ncout_open

      subroutine inst_ncout_close
      status = nf_close(fid)
      return
      end subroutine inst_ncout_close

      subroutine inst_setup_var(varname)
      implicit none
      include 'netcdf.inc'
      character(len=30) :: varname
      character(len=30) :: yesno
      integer :: n,ivar
c get the id number for this variable
      status = nf_inq_varid(fid,trim(varname),varid)
      if(status.ne.nf_noerr) return
c check whether this variable is already in the database
      ivar=0
      do n=1,nvar
        if(varid.eq.varid_list(n)) then
          ivar=n
          ndims=ndims_list(n)
          dimids(1:ndims)=dimids_list(1:ndims,n)
          dimlens(1:ndims)=dimlens_list(1:ndims,n)
          exit
        endif
      enddo
      if(ivar.eq.0) then ! get info about new variable
        status = nf_inq_varndims(fid,varid,ndims)
        status = nf_inq_vardimid(fid,varid,dimids)
        do n=1,ndims
          status = nf_inq_dimlen(fid,dimids(n),dimlens(n))
        enddo
        nvar = nvar+1
        ivar = nvar
        varid_list(nvar)=varid
        ndims_list(nvar)=ndims
        dimids_list(1:ndims,nvar)=dimids(1:ndims)
        dimlens_list(1:ndims,nvar)=dimlens(1:ndims)
        yesno=''
        status = nf_get_att_text(fid,varid,'write_all',yesno)
        if(status.eq.nf_noerr .and. trim(yesno).eq.'yes') then
          write_all_list(nvar)=.true.
        endif
      endif
      do n=1,ndims-1
         srt(n) = 1
         cnt(n) = dimlens(n)
      enddo
      n = ndims
      if(write_all_list(ivar)) then ! write the entire array
         srt(n) = 1
         cnt(n) = dimlens(n)
      else ! only the next element of the outermost dimension
        cnt(n) = 1
c check whether the last dimension is unlimited
        if(dimids(n).eq.unlimdimid) then ! yes, unlimited.
          counts_list(ivar) = counts_list(ivar)+1
        else                    ! cyclical time dimension
          counts_list(ivar) = counts_list(ivar)+1
          if(counts_list(ivar).gt.dimlens(n)) counts_list(ivar)=1
        endif
        srt(n) = counts_list(ivar)
      endif
      return
      end subroutine inst_setup_var

      end module inst_ncout

      subroutine inst_ncoutd(varname,var)
c write an arbitrary real*8 array predefined in the output file
      use inst_ncout
      implicit none
      include 'netcdf.inc'
      character(len=30) :: varname
      real*8 :: var(1)
      call inst_setup_var(varname)
      status = nf_put_vara_double(fid,varid,srt,cnt,var)
      return
      end subroutine inst_ncoutd

      subroutine inst_ncouti(varname,var)
c write an arbitrary integer array predefined in the output file
      use inst_ncout
      implicit none
      include 'netcdf.inc'
      character(len=30) :: varname
      integer :: var(1)
      call inst_setup_var(varname)
      status = nf_put_vara_int(fid,varid,srt,cnt,var)
      return
      end subroutine inst_ncouti

      subroutine inst_ncout_ij(varname,xij)
!@sum  inst_ncout_ij output lon-lat binary records
!@sum  data is written from the root processor
!@auth M. Kelley
!@ver  1.0
      use inst_ncout
      use MODEL_COM, only : im
      use domain_decomp, only: grid, pack_data, am_i_root
      implicit none
      character(len=30) :: varname
      real*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) :: xij
      call pack_data(grid, xij, xij_glob)
      if(am_i_root()) call inst_ncoutd(varname,xij_glob)
      return
      end subroutine inst_ncout_ij

      subroutine inst_ncout_ijl(varname,xijl)
!@sum  inst_ncout_ijl output lon-lat-height binary records
!@sum  data is written from the root processor
!@auth M. Kelley
!@ver  1.0
      use inst_ncout
      use MODEL_COM, only : im,lm
      use domain_decomp, only: grid, pack_data, am_i_root
      implicit none
      character(len=30) :: varname
      real*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: xijl
      call pack_data(grid, xijl, xijl_glob)
      if(am_i_root()) call inst_ncoutd(varname,xijl_glob)
      return
      end subroutine inst_ncout_ijl

      subroutine inst_ncout_lij(varname,xlij)
!@sum  inst_ncout_lij output height-lon-lat binary records
!@sum  data is written from the root processor
!@auth M. Kelley
!@ver  1.0
      use inst_ncout
      use MODEL_COM, only : im,lm
      use domain_decomp, only: grid, pack_column, am_i_root
      implicit none
      character(len=30) :: varname
      real*8, dimension(lm,im,grid%J_STRT_HALO:grid%J_STOP_HALO) :: xlij
      call pack_column(grid, xlij, xlij_glob)
      if(am_i_root()) call inst_ncoutd(varname,xlij_glob)
      return
      end subroutine inst_ncout_lij
