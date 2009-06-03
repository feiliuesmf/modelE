      program write_giss2d
!@sum For every variable in the input file having the
!@+   two dimensions specified on the command line,
!@+   write a fortran record with title*80,data to
!@+   the output file.  The two dimensions need not
!@+   be the first two of a given variable, nor consecutive.
!@+   Variables with more than these two dimensions are
!@+   written one slab at a time.
!@+   The title is created from attributes long_name/units.
!@+   If long_name is absent, the netcdf variable name is used.
!@auth M. Kelley
      implicit none
      real*4, dimension(:), allocatable :: xout
      character(len=80) :: title,infile,outfile,lname
      character(len=40) :: vname,dimname1,dimname2,units
      integer :: dsiz1,dsiz2,nargs,lunit,k,n,nslab

      include 'netcdf.inc'
      integer :: fid,status,varid,varid2,nvars,ndims,did1,did2,idim,jdim
      integer, dimension(7) :: dids,srt,cnt,dsizes,kmod,p1,p2
      character(len=30), dimension(7) :: dnames
      character(len=30) :: diminfo
      character(len=1) :: str1
      character(len=3) :: str3
      character(len=6) :: ifmt='(ix.x)'
      real*4, parameter :: undef=-1e30
      real*4 :: shnhgm(3)

      nargs = iargc()
      if(nargs.ne.4) then
        write(6,*)
     &       'usage: write_giss2d infile outfile dimname1 dimname2'
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

c
c allocate workspace
c
      allocate(xout(dsiz1*dsiz2))

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
        dids = -1
        status = nf_inq_vardimid(fid,varid,dids)
        if(count(dids.eq.did1).ne.1) cycle
        if(count(dids.eq.did2).ne.1) cycle
        do n=1,ndims
          if(dids(n).eq.did1) idim=n
          if(dids(n).eq.did2) jdim=n
          status = nf_inq_dimlen(fid,dids(n),dsizes(n))
          status = nf_inq_dimname(fid,dids(n),dnames(n))
        enddo
        srt = 1
        cnt = 1
        cnt(idim) = dsizes(idim)
        cnt(jdim) = dsizes(jdim)
        nslab = product(dsizes(1:ndims))/(dsiz1*dsiz2)
        if(ndims.gt.2) then
          k = 1
          diminfo=''
          do n=1,ndims
            if(n.eq.idim .or. n.eq.jdim) cycle
            kmod(n) = k
            k = k*dsizes(n)
            write(str1,'(i1)') int(1.+log10(real(dsizes(n))))
            ifmt(3:3) = str1
            ifmt(5:5) = str1
            write(str3,ifmt) srt(n)
            if(len_trim(diminfo).eq.0) then
              diminfo=trim(dnames(n))//'='//trim(str3)
            else
              diminfo=
     &             trim(diminfo)//' '//trim(dnames(n))//'='//trim(str3)
            endif
            p1(n) = len_trim(diminfo)-len_trim(str3)+1
            p2(n) = len_trim(diminfo)
          enddo
        endif
        status = nf_inq_varname(fid,varid,vname)
        lname = vname
        status = nf_get_att_text(fid,varid,'long_name',lname)
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        if(status.eq.nf_noerr) units = ' ('//trim(units)//') '
        shnhgm = undef
        if(ndims.eq.2) then ! look for global means
          status = nf_inq_varid(fid,trim(vname)//'_hemis',varid2)
          if(status.eq.nf_noerr) then
            status = nf_get_var_real(fid,varid2,shnhgm)
          endif
        endif
        do k=1,nslab
          status = nf_get_vara_real(fid,varid,srt,cnt,xout)
          title = trim(lname)//units
          if(ndims.gt.2) title=trim(title)//' '//trim(diminfo)
          if(shnhgm(3).eq.undef) then ! no global mean available
            write(lunit) title,xout
          else                        ! write with global mean
            write(lunit) title,xout
     &           ,(undef,n=1,dsiz2)   ! have to write means at each lat
     &           ,shnhgm(3)           ! before the global mean
          endif
          if(ndims.gt.2) then
            do n=1,ndims        ! increment the start vector
              if(n.eq.idim .or. n.eq.jdim) cycle
              if(mod(k,kmod(n)).eq.0) then
                srt(n) = srt(n) + 1
                if(srt(n).gt.dsizes(n)) srt(n)=1
                write(str3,ifmt) srt(n)
                diminfo(p1(n):p2(n))=trim(str3)
              endif
            enddo
          endif
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(xout)

c
c close fortran output file
c
      close(lunit)

      end program write_giss2d

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
