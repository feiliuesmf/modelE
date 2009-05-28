      program write_2d_as_giss4d
!@sum For every variable in the input file having the
!@+   two dimensions specified on the command line,
!@+   write a "giss4d" fortran record to
!@+   the output file.
!@+   The title is created from attributes long_name/units.
!@+   If long_name is absent, the netcdf variable name is used.
!@auth M. Kelley
      implicit none
      real*4, dimension(:), allocatable :: coord1,coord2,vmean
      real*4, dimension(:,:), allocatable :: xout,shnhgm
      character(len=80) :: title,infile,outfile,lname
      character(len=40) :: vname,dimname1,dimname2,units
      integer :: dsiz1,dsiz2,nargs,lunit,k,n
      character*16 :: cx,cy
      character*16, parameter :: cblank = '                '
      real*4, parameter :: undef=-1e30,one=1.
      include 'netcdf.inc'
      integer :: fid,status,varid,varid2,nvars,ndims,did1,did2,idim,jdim
      integer, dimension(7) :: dids

      nargs = iargc()
      if(nargs.ne.4) then
        write(6,*)
     &   'usage: write_2d_as_giss4d infile outfile dimname1 dimname2'
        stop
      endif

      call getarg(1,infile)
      call getarg(2,outfile)
      call getarg(3,dimname1)
      call getarg(4,dimname2)

      call handle_err(nf_open(infile,nf_nowrite,fid),
     &     'opening '//trim(infile))

c
c Get info on the two requested dimensions
c
      call get_dimsize(fid,trim(dimname1),dsiz1)
      call get_dimsize(fid,trim(dimname2),dsiz2)
      status = nf_inq_dimid(fid,trim(dimname1),did1)
      status = nf_inq_dimid(fid,trim(dimname2),did2)
      allocate(coord1(dsiz1),coord2(dsiz2))
      call get_var_real(fid,trim(dimname1),coord1)
      call get_var_real(fid,trim(dimname2),coord2)
      status = nf_inq_varid(fid,trim(dimname1),varid)
      cx=''
      cy=''
      if(nf_get_att_text(fid,varid,'giss_name',cx).ne.nf_noerr)
     &     cx = dimname1
      status = nf_inq_varid(fid,trim(dimname2),varid)
      if(nf_get_att_text(fid,varid,'giss_name',cy).ne.nf_noerr)
     &     cy = dimname2

c
c allocate workspace
c
      allocate(xout(dsiz1,dsiz2))
      allocate(vmean(dsiz1+3))
      allocate(shnhgm(3,dsiz2))

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c open fortran output file
c
      lunit = 10
      open(lunit,file=trim(outfile),form='unformatted',
     &     convert='big_endian')

c
c Loop over quantities in the file
c
      do varid=1,nvars
        status = nf_inq_varndims(fid,varid,ndims)
        if(ndims.ne.2) cycle
        dids = -1
        status = nf_inq_vardimid(fid,varid,dids)
        if(count(dids.eq.did1).ne.1) cycle
        if(count(dids.eq.did2).ne.1) cycle
        status = nf_inq_varname(fid,varid,vname)
        if(dids(1).ne.did1) then
          write(6,*) 'dimensions for variable ',trim(vname),
     &         ' match but are reversed: skipping'
          cycle
        endif
        lname = vname
        status = nf_get_att_text(fid,varid,'long_name',lname)
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        if(status.eq.nf_noerr) units = ' ('//trim(units)//') '
        title = trim(lname)//units
        status = nf_get_var_real(fid,varid,xout)
c look for horizontal and vertical means
        shnhgm = undef; vmean = undef
        if(nf_inq_varid(fid,trim(vname)//'_hemis',varid2).eq.nf_noerr)
     &       status = nf_get_var_real(fid,varid2,shnhgm)
        if(nf_inq_varid(fid,trim(vname)//'_vmean',varid2).eq.nf_noerr)
     &       status = nf_get_var_real(fid,varid2,vmean)
c write the record
        write(lunit) title,
     &       dsiz1,dsiz2,1,1,       ! dimension sizes
     &       xout,                  ! the field
     &       coord1,coord2,one,one, ! coordinate axes
     &       cx,cy,cblank,cblank,   ! names of coordinate axes
     &       'NASAGISS',            !
     &       vmean,                 ! vertical mean
     &       shnhgm                 ! horizontal means
      enddo

c
c deallocate workspace
c
      deallocate(xout,coord1,coord2,vmean,shnhgm)

c
c close fortran output file
c
      close(lunit)

      end program write_2d_as_giss4d

      subroutine handle_err(status,errmsg)
      implicit none
      integer :: status
      character(len=*) :: errmsg
      include 'netcdf.inc'
      if(status.ne.nf_noerr) then
        write(6,*) 'error '//trim(errmsg)
        stop
      endif
      end subroutine handle_err
