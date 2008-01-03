c The routines inst_ncoutd and inst_ncouti can be used to output a model
c variable VAR whose dimensions match those of the variable VARNAME
c in a predefined netcdf file nctemp.nc in the run directory.
c VARNAME is assumed to be a 30-byte character string.
c
c       subroutine inst_ncoutd(VARNAME,VAR) ! arbitrary REAL*8 VAR
c       subroutine inst_ncouti(VARNAME,VAR) ! arbitrary INTEGER VAR
c
c These routines output a VAR of arbitrary shape and size from a single
c processor.  To output data distributed over multiple processors using
c the ESMF domain decomposition, two methods are available.
c
c Method 1:  Each processor writes to an appropriately resized local copy
c   of nctemp.nc, called nctempNNN.nc where NNN is the processor number.
c   These local files are automatically generated at runtime and contain
c   halo/domain information allowing them to be stitched together at
c   the completion of a run with a call to inst_ncout_stitch.
c   This method requires that nctemp.nc indicate which of its arbitrarily
c   named dimensions correspond to a longitude or latitude dimension
c   (see the "coordinate variable" section of netcdf_template.cdl).
c
c Method 2:  Distributed data is gathered to the root processor before
c   output.  The following routines must therefore be used:
c       subroutine inst_ncout_ij(VARNAME,VAR) ! lon-lat VAR
c       subroutine inst_ncout_ijl(VARNAME,VAR) ! lon-lat-level VAR
c       subroutine inst_ncout_lij(VARNAME,VAR) ! level-lon-lat VAR
c   Overloading has not yet been implemented to allow ij/ijl/lij arrays
c   to all be handled by a single routine.
c
c The module routines inst_ncout_open/inst_ncout_close must be called
c at the beginning/end of model execution.  They take no arguments.
c If using Method 2 in multiprocessor mode, inst_ncout_open_global
c should be called instead (and the open/close calls should each be
c done from the root processor).
c If using Method 1 in multiprocessor mode, the files nctempNNN.nc can
c be stitched together at the end of model execution with a call to
c inst_ncout_stitch.  This step is optional and in some situations it
c may be preferable to employ an offline version of inst_ncout_stitch.

      module inst_ncout
!@sum  inst_ncout facilitates the writing of model variables in netcdf format
!@auth M. Kelley
      use MODEL_COM, only : im,jm,lm
      implicit none
      save
      include 'netcdf.inc'

      private
      public ::
     &     inst_ncout_open,inst_ncout_close,inst_setup_var
     &    ,inst_ncout_open_global
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
#ifdef USE_ESMF
      USE ESMF_MOD, Only: ESMF_VMGet
      USE ESMF_CUSTOM_MOD, Only: vm => modelE_vm
      use domain_decomp, only: grid
#endif
      integer :: gfid
      character(len=80) :: gfile ! global template file
#ifdef USE_ESMF
      character(len=3) :: pet_str
      integer :: my_pet,jm_local,halo_width,jstart
#endif
      gfile='nctemp.nc'
#ifdef USE_ESMF
      call ESMF_VMGet(vm, localPET = my_pet)
      write(pet_str,'(i3.3)') 1+my_pet
      outfile=gfile(1:len_trim(gfile)-3)//pet_str//'.nc'
      status = nf_open(trim(gfile), nf_nowrite, gfid)
#else
      status = nf_open(trim(gfile), nf_write, fid)
#endif
      if(status.ne.nf_noerr) then
         write(*,*) 'cannot open: '//trim(gfile)
         call stop_model('stopped in INST_netcdf.f',255)
      endif
#ifdef USE_ESMF
c create the local file for this processor
      status = nf_create(trim(outfile), nf_clobber, fid)
      jstart = grid%j_strt
      halo_width = grid%j_strt-grid%j_strt_halo
      jm_local = grid%j_stop_halo-grid%j_strt_halo+1
      call resize_file(gfid,fid,jm_local,jstart,halo_width)
      status = nf_close(gfid)
#endif
      call init_ncout_data
      return
      end subroutine inst_ncout_open

      subroutine inst_ncout_open_global
      outfile='nctemp.nc'
      status = nf_open (trim(outfile), nf_write, fid)
      if(status.ne.nf_noerr) then
         write(*,*) 'cannot open: '//trim(outfile)
         call stop_model('stopped in INST_netcdf.f',255)
      endif
      call init_ncout_data
      return
      end subroutine inst_ncout_open_global

      subroutine init_ncout_data
      nvar = 0
      counts_list(:) = 0
      write_all_list(:) = .false.
      status = nf_inq_unlimdim(fid,unlimdimid)
      end subroutine init_ncout_data

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

      subroutine resize_file(ifid,ofid,jm_local,jstart,halo_width)
c resizes latitudinal dimensions to jm_local
c longitudinal resizing can be easily added when needed
      implicit none
      include 'netcdf.inc'
      integer :: ifid,ofid,jm_local,jstart,halo_width
      integer :: status,nvars,ivarid,dimid,n
      logical :: is_coord,needs_resiz
      integer :: nresiz
      integer, dimension(20) ::
     &     dims2resiz,vars2resiz,newsiz,global_start
      integer, parameter :: nqsj=2
      character(len=2), dimension(nqsj), parameter ::
     &     jgrid_strings=(/ 'cj', 'ej' /)
c
c first find the dimensions that need resizing
c
      status = nf_inq_nvars(ifid,nvars)
      nresiz = 0
      do ivarid=1,nvars
        call resiz_query(ifid,ivarid,is_coord,dimid,needs_resiz
     &       ,jgrid_strings,nqsj)
        if(needs_resiz) then
          nresiz = nresiz + 1
          vars2resiz(nresiz) = ivarid
          dims2resiz(nresiz) = dimid
          newsiz(nresiz) = jm_local
          global_start(nresiz) = jstart
        endif
      enddo
c
c copy dimension and variable definitions from old to new file
c
      call copy_dims(ifid,ofid,nresiz,dims2resiz,newsiz)
      call copy_vars(ifid,ofid)
c
c and add information about the global location and halo width
c
      do n=1,nresiz
        ivarid = vars2resiz(n)
        status = nf_put_att_int(
     &       ofid,ivarid,'global_start',nf_int,1,global_start(n))
        status = nf_put_att_int(
     &       ofid,ivarid,'halo_width',nf_int,1,halo_width)
      enddo
      status = nf_enddef(ofid)
      return
      end subroutine resize_file

      subroutine copy_vars(ifid,ofid)
      implicit none
      include 'netcdf.inc'
      integer :: ifid,ofid
      integer :: status,ivarid,nvars,vtype,ndims,dimids(7),natts,ovid,n
      character(len=20) :: vname,att_name
      status = nf_inq_nvars(ifid,nvars)
      do ivarid=1,nvars
        status = nf_inq_var(ifid,ivarid,vname,vtype,ndims,dimids,natts)
        status = nf_def_var(ofid,vname,vtype,ndims,dimids,ovid)
        do n=1,natts
          status = nf_inq_attname(ifid,ivarid,n,att_name)
          status = nf_copy_att(ifid,ivarid,att_name,ofid,ovid)
        enddo
      enddo
c copy global attributes
      do n=1,natts
        status = nf_inq_attname(ifid,nf_global,n,att_name)
        status = nf_copy_att(ifid,nf_global,att_name,ofid,nf_global)
      enddo
      return
      end subroutine copy_vars

      subroutine copy_dims(ifid,ofid,nresiz,dims2resiz,newsiz)
      implicit none
      include 'netcdf.inc'
      integer :: ifid,ofid,nresiz
      integer, dimension(nresiz) :: dims2resiz,newsiz
      integer :: status,ndims,idimid,dimsiz,n
      character(len=20) :: dimname
      status = nf_inq_ndims(ifid,ndims)
      do idimid=1,ndims
        if(any(idimid.eq.dims2resiz)) then
          n=1
          do while(idimid.ne.dims2resiz(n))
            n=n+1
          enddo
          dimsiz = newsiz(n)
        else
          status = nf_inq_dimlen(ifid,idimid,dimsiz)
        endif
        status = nf_inq_dimname(ifid,idimid,dimname)
        status = nf_def_dim(ofid,dimname,dimsiz,idimid)
      enddo
      return
      end subroutine copy_dims

      subroutine resiz_query(fid,vid,is_coord,dimid,needs_resiz,
     &     query_strings,nqs)
      implicit none
      include 'netcdf.inc'
      integer :: fid,vid,dimid,nqs
      logical :: is_coord,needs_resiz
      character(len=2), dimension(nqs) :: query_strings
      character(len=80) :: varname
      integer :: status,vartype,ndims,dimids(7),natts,att_type,n
      character(len=20) :: dimname,grid_str
      needs_resiz=.false.
      ndims = 0
      status = nf_inq_var(fid,vid,varname,vartype,ndims,dimids,natts)
      if(ndims.eq.0 .or. ndims.gt.1) then
        is_coord=.false.
        return
      endif
      status = nf_inq_dimname(fid,dimids(1),dimname)
      if(trim(dimname).ne.trim(varname)) then
        is_coord=.false.
        return
      endif
c
c see whether this dimension needs resizing
c
      status = nf_inq_atttype(fid,vid,'grid',att_type)
      if(status.eq.nf_noerr .and. att_type.eq.nf_char) then
        grid_str=''
        status = nf_get_att_text(fid,vid,'grid',grid_str)
        do n=1,nqs
          if(grid_str(1:2).eq.query_strings(n)) then
            dimid=dimids(1)
            needs_resiz=.true.
          endif
        enddo
      endif
      return
      end subroutine resiz_query

      subroutine inst_ncout_stitch
#ifdef USE_ESMF
      USE ESMF_MOD, Only: ESMF_VMGet
      USE ESMF_CUSTOM_MOD, Only: vm => modelE_vm
      use domain_decomp, only: grid, am_i_root
      implicit none
      include 'netcdf.inc'
      integer :: status,ifid,gfid,ivarid,ovarid
     &     ,status_read,status_write
      character(len=80) :: varname,gfile,lfile
      character(len=3) :: nstr
      integer :: varsize,v4size,n,nvars,npes
      logical :: is_coord
      integer :: vartype,ndims,srt(7),srtin(7),cnt(7)
      real*4, dimension(:), allocatable :: v4
      if(.not.am_i_root()) return
c open global file
      gfile='nctemp.nc'
      status = nf_open(trim(gfile),nf_write,gfid)
      if(status.ne.nf_noerr) then
         write(*,*) 'cannot open global file '//trim(gfile)
         call stop_model('stopped in INST_netcdf.f',255)
      endif
      v4size=0
c
c Loop over local files
c
      call ESMF_VMGet(vm, petCount = npes)
      file_loop: do n=1,npes
      write(nstr,'(i3.3)') n
      lfile=gfile(1:len_trim(gfile)-3)//nstr(1:3)//'.nc'
      status = nf_open(trim(lfile),nf_nowrite,ifid)
      if(status.ne.nf_noerr) then
         write(*,*) 'cannot open local file '//trim(lfile)
         call stop_model('stopped in INST_netcdf.f',255)
      endif
c
c Loop over the variables in this file
c
      status = nf_inq_nvars(ifid,nvars)
      var_loop: do ivarid=1,nvars
      call get_stitch_info(ifid,ivarid,varname,vartype,is_coord
     &     ,ndims,cnt,srt,srtin)
      if(is_coord) cycle var_loop
      status = nf_inq_varid(gfid,trim(varname),ovarid)
      if(status.ne.nf_noerr) then
        write(*,*) 'nonexistent global file variable '//
     &       trim(varname)
         call stop_model('stopped in INST_netcdf.f',255)
      endif
c allocate space and read/write variable from/to local/global file
      varsize=product(cnt(1:ndims))
      if(vartype.eq.nf_double) varsize=2*varsize
      if(varsize.gt.v4size) then
        if(allocated(v4)) deallocate(v4)
        allocate(v4(varsize),stat=status)
        if(status.ne.0) then
          write(*,*) 'memory allocation error'
          call stop_model('stopped in INST_netcdf.f',255)
        endif
        v4size=varsize
      endif
      if(vartype.eq.nf_int  .or. vartype.eq.nf_float .or.
     &   vartype.eq.nf_int1 .or. vartype.eq.nf_int2) then
        status_read  = nf_get_vara_real(ifid,ivarid,srtin,cnt,v4)
        status_write = nf_put_vara_real(gfid,ovarid,srt,cnt,v4)
      else if(vartype.eq.nf_double) then
        status_read  = nf_get_vara_double(ifid,ivarid,srtin,cnt,v4)
        status_write = nf_put_vara_double(gfid,ovarid,srt,cnt,v4)
      else if(vartype.eq.nf_char) then
        status_read  = nf_get_vara_text(ifid,ivarid,srtin,cnt,v4)
        status_write = nf_put_vara_text(gfid,ovarid,srt,cnt,v4)
      else
        write(*,*) 'unsupported data type for '//trim(varname)
          call stop_model('stopped in INST_netcdf.f',255)
      endif
      if(status_read.ne.nf_noerr) then
        write(*,*) 'unsuccessful read of variable '//trim(varname)//
     &       ' error = ',status
          call stop_model('stopped in INST_netcdf.f',255)
      endif
      if(status_write.ne.nf_noerr) then
        write(*,*) 'unsuccessful write of variable '//trim(varname)//
     &       ' error = ',status
          call stop_model('stopped in INST_netcdf.f',255)
      endif
      enddo var_loop
      status = nf_close(ifid)
      enddo file_loop
      deallocate(v4)
      status = nf_close(gfid)
#endif
      return
      end subroutine inst_ncout_stitch

      subroutine get_stitch_info(fid,vid,varname,vartype,is_coord
     &     ,ndims,dimsizes,srt,srtin)
      implicit none
      include 'netcdf.inc'
      character(len=80) :: varname
      integer :: fid,vid,vartype,ndims,dimsizes(7),srt(7),srtin(7)
      logical :: is_coord
      integer :: status,idim,natts,dimids(7),cvarid,att_type,halo_width
      character(len=20) :: dimname
      status = nf_inq_var(fid,vid,varname,vartype,ndims,dimids,natts)
      do idim=1,ndims
        status = nf_inq_dim(fid,dimids(idim),dimname,dimsizes(idim))
        if(trim(dimname).eq.trim(varname)) then
          is_coord=.true.
          return
        else
          is_coord=.false.
        endif
c
c obtain the offsets for this dimension
c
        srtin(idim) = 1         ! if haloing, will be changed below
        status = nf_inq_varid(fid,trim(dimname),cvarid)
        if(status.eq.nf_noerr) then ! dim may have coordinate info
          status = nf_inq_atttype(fid,cvarid,'global_start',att_type)
          if(status.ne.nf_noerr) then
            srt(idim) = 1       ! default
          else
            if(att_type.ne.nf_int) then
              write(*,*)trim(dimname)//' has invalid offset information'
              call stop_model('stopped in INST_netcdf.f',255)
            endif
            status = nf_get_att_int(fid,cvarid,'global_start',srt(idim))
          endif
          status = nf_inq_atttype(fid,cvarid,'halo_width',att_type)
          if(status.eq.nf_noerr .and. att_type.eq.nf_int) then
            status = nf_get_att_int(fid,cvarid,'halo_width',halo_width)
            dimsizes(idim) = dimsizes(idim) - 2*halo_width
            srtin(idim) = 1+halo_width
          endif
        else
          srt(idim) = 1         ! default
        endif
      enddo
      return
      end subroutine get_stitch_info

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
