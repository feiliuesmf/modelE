      integer :: fid
      character(len=*) :: varname
      type(dist_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      logical, intent(in), optional :: no_xdim
      real*8, dimension(:), allocatable :: arr1d
      integer :: nl,nk,jd,nd,nsend,ierr,src_pe,tag,modpe,vid,rc
      integer :: mpi_status(mpi_status_size)
      integer, dimension(7) :: srt,cnt
      logical :: has_xdim
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in output file - stopping'
      endif
      call stoprc(rc,nf_noerr)
      tag=1024
      if(present(jdim)) then
        jd=jdim
      else
        jd=2
      endif
      if(present(no_xdim)) then
        has_xdim = .not.no_xdim
      elseif(jd.eq.1) then
        has_xdim = .false.
      else
        has_xdim = .true.
      endif
      nd = size(shape(arr))
c if latlon model, call routines better suited for very thin
c zonal strips if they are available
      if(grid%ntiles.eq.1 .and. has_xdim .and. jd.ne.4 .and.
     &     (jd.eq.2 .or. jd.lt.nd)) then
        call par_write_jdecomp_optimized(grid,fid,vid,arr,jd)
        return
      endif
      if(grid%am_i_rowroot) then
        call get_nlnk(arr,jd,nl,nk)
        nsend = nl*nk*grid%maxnj
        if(has_xdim) then
          nsend = nsend*grid%npx
        elseif(jd.gt.1) then
          nsend = nsend*size(arr,jd-1)
        endif
        allocate(arr1d(nsend))
        srt(:) = 1
        cnt(1:nd) = shape(arr)
        if(has_xdim) cnt(jd-1) = grid%npx
        cnt(nd+1) = 1
      endif
      if(has_xdim) then
        call pack_row(grid,arr,arr1d,jdim=jd)
      else
        call pack_row_no_xdim(grid,arr,arr1d,jd)
      endif
      if(grid%am_i_globalroot) then
        do src_pe=0,grid%nproc-grid%nprocx,grid%nprocx
          if(src_pe.gt.0) then
            call mpi_recv(arr1d,nsend,MPI_DOUBLE_PRECISION,src_pe,
     &           tag,MPI_COMM_WORLD,mpi_status,ierr)
          endif
          modpe=mod(src_pe,grid%nproc_tile)+1
          srt(jd) = grid%jsr(modpe)
          cnt(jd) = 1+grid%jer(modpe)-grid%jsr(modpe)
          srt(nd+1) = 1 + src_pe/grid%nproc_tile
          rc = nf_put_vara_double(fid,vid,srt,cnt,arr1d)
        enddo
      elseif(grid%am_i_rowroot) then
        call mpi_send(arr1d,nsend,MPI_DOUBLE_PRECISION,0,
     &       tag,MPI_COMM_WORLD,ierr)
      endif
      if(grid%am_i_rowroot) deallocate(arr1d)
