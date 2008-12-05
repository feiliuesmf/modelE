      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      real*8, dimension(:), allocatable :: arr1d
      integer :: nl,nk,jd,nd,nsend,ierr,src_pe,tag,modpe,vid,rc
      integer :: mpi_status(mpi_status_size)
      integer, dimension(7) :: srt,cnt
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
      if(grid%am_i_rowroot) then
        call get_nlnk(arr,jd,nl,nk)
        nsend = nl*nk*grid%npx*grid%maxnj
        allocate(arr1d(nsend))
        srt(:) = 1
        nd = size(shape(arr))
        cnt(1:nd) = shape(arr)
        cnt(jd-1) = grid%npx
        cnt(nd+1) = 1
      endif
      call pack_row(grid,arr,arr1d,jdim=jd)
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
