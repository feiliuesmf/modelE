      integer :: fid
      character(len=*) :: varname
      type(dd2d_grid), intent(in) :: grid
      integer, intent(in), optional :: jdim
      real*8, dimension(:), allocatable :: arr1d
      integer :: nl,nk,jd,nd,nsend,ierr,dst_pe,tag,modpe,vid,rc
      integer :: mpi_status(mpi_status_size)
      integer, dimension(7) :: srt,cnt
      if(grid%am_i_globalroot) then
        rc = nf_inq_varid(fid,trim(varname),vid)
        if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &       trim(varname),' not found in input file - skipping'
      endif
c      write(*,*) trim(varname)
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(rc.ne.nf_noerr) return
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
      if(grid%am_i_globalroot) then
c loop in reverse order so that arr1d filled with roots row at end of loop
        do dst_pe=grid%nproc-grid%nprocx,0,-grid%nprocx
c read from netcdf file
          modpe=mod(dst_pe,grid%nproc_tile)+1
          srt(jd) = grid%jsr(modpe)
          cnt(jd) = 1+grid%jer(modpe)-grid%jsr(modpe)
          srt(nd+1) = 1 + dst_pe/grid%nproc_tile
          rc = nf_get_vara_double(fid,vid,srt,cnt,arr1d)
          if(dst_pe.gt.0) then
            call mpi_send(arr1d,nsend,MPI_DOUBLE_PRECISION,dst_pe,
     &           tag,MPI_COMM_WORLD,ierr)
c	    write(*,*) "glob root sends packet",grid%gid
          endif
        enddo
      elseif(grid%am_i_rowroot) then
        call mpi_recv(arr1d,nsend,MPI_DOUBLE_PRECISION,0,
     &       tag,MPI_COMM_WORLD,mpi_status,ierr)
c	    write(*,*) "row root received packet",grid%gid
      endif
c      write(*,*) "before unpack_row, gid=",grid%gid
      call unpack_row(grid,arr1d,arr,jdim=jd)
c      write(*,*) "after unpack_row, gid=",grid%gid
      if(grid%am_i_rowroot) deallocate(arr1d)

