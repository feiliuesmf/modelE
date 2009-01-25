      type(dist_grid) :: grid
      integer :: fid
      character(len=*) :: varname,attname
      integer :: attlen
      integer, intent(in), optional :: attnum
      integer :: rc,vid,ierr,lname,l
      character :: tmpstr(1000)
      if(grid%am_i_globalroot) then
        if(trim(varname).eq.'global') then
          rc = nf_noerr
          vid = nf_global
        else
          rc = nf_inq_varid(fid,trim(varname),vid)
          if(rc.ne.nf_noerr) write(6,*) 'variable ',
     &         trim(varname),' not found in input file - skipping'
        endif
      endif
      call mpi_bcast(rc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(rc.ne.nf_noerr) return
c
c Get the name and length of this attribute and broadcast them.
c For some reason, mpi_bcast does not work on strings declared
c as character(len=xxx), so we have to copy attname into/out of
c a conventional array
c
      if(grid%am_i_globalroot) then
        if(present(attnum)) then
          rc = nf_inq_attname(fid,vid,attnum,attname)
          lname = len_trim(attname)
          do l=1,lname
            tmpstr(l) = attname(l:l)
          enddo
        endif
        rc = nf_inq_attlen(fid,vid,trim(attname),attlen)
      endif
      call mpi_bcast(attlen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(present(attnum)) then
        call mpi_bcast(lname,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(tmpstr,lname,MPI_CHARACTER,0,
     &       MPI_COMM_WORLD,ierr)
        attname=''
        do l=1,lname
          attname(l:l) = tmpstr(l)
        enddo
      endif
