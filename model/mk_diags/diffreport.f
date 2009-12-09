!@sum diffreport reports on differences between arrays in 2 netcdf files
!@auth M. Kelley
      program diffreport
      implicit none
      include 'netcdf.inc'
      integer :: status,vtype,varid1,varid2,fid1,fid2
      character(len=80) :: file1,file2
      character(len=40) :: vname,pos_str
      integer :: iargc,nvars,ndims,dsizes(7)
      integer :: arrsize1,arrsize2
      real*8, dimension(:), allocatable :: arr1,arr2,arrdiff
      integer :: nmax_abs(1),nmax_rel(1)
      real*8 :: max_abs,max_rel
c
c get the number of arguments
c
      if(iargc().ne.2) then
        write(6,*) 'usage: diffreport file1 file2'
        stop
      endif

c
c open input files
c
      call getarg(1,file1)
      call getarg(2,file2)
      call handle_err( nf_open(file1,nf_nowrite,fid1),
     &     'nonexistent/non-netcdf input file '//trim(file1) )
      call handle_err( nf_open(file2,nf_nowrite,fid2),
     &     'nonexistent/non-netcdf input file '//trim(file2) )

c
c loop over arrays
c
      status = nf_inq_nvars(fid1,nvars)

      do varid1=1,nvars
c skip character arrays
        status = nf_inq_vartype(fid1,varid1,vtype)
        if(vtype.eq.nf_char) cycle
c get array name
        status = nf_inq_varname(fid1,varid1,vname)
c skip certain modelE arrays
        if(trim(vname).eq.'cputime') cycle
c
        status = nf_inq_varid(fid2,vname,varid2)
c skip arrays not present in both files
        if(status.ne.nf_noerr) then
          write(6,*) 'file1 array '//trim(vname)//
     &         ' is absent in file2: skipping'
          cycle
        endif
c get array sizes and allocate space
        call get_varsize(fid1,vname,arrsize1)
        call get_varsize(fid2,vname,arrsize2)
c if array sizes do not match, skip
        if(arrsize1.ne.arrsize2) then
          write(6,*) 'array '//trim(vname)//
     &         ' has different sizes in file1/file2: skipping'
          cycle
        endif
        allocate(arr1(arrsize1),arr2(arrsize2))
c read arrays
        status = nf_get_var_double(fid1,varid1,arr1)
        status = nf_get_var_double(fid2,varid2,arr2)
c check for differences
        if(any(arr1.ne.arr2)) then
          write(6,*) trim(vname)//' max diffs:'
          call get_vdimsizes(fid1,vname,ndims,dsizes)
          allocate(arrdiff(arrsize1))
c absolute
          arrdiff = abs(arr1-arr2)
          max_abs = maxval(arrdiff)
          nmax_abs = maxloc(arrdiff)
          call get_pos_str(nmax_abs,ndims,dsizes,pos_str)
          write(6,*) '          abs: ',max_abs,trim(pos_str)
c relative
          arrdiff = arrdiff/max(abs(arr1),abs(arr2))
          max_rel = maxval(arrdiff)
          nmax_rel = maxloc(arrdiff)
          call get_pos_str(nmax_rel,ndims,dsizes,pos_str)
          write(6,*) '          rel: ',max_rel,trim(pos_str)
c
          deallocate(arrdiff)
        endif

c deallocate arrays
        deallocate(arr1,arr2)
      enddo

c
c close input files
c
      status = nf_close(fid1)
      status = nf_close(fid2)

      end program diffreport

      subroutine get_pos_str(n,ndims,dsizes,pos_str)
      implicit none
      integer :: n,ndims
      integer :: dsizes(ndims)
      character(len=40) :: pos_str
      integer :: inds(7)
      integer :: idim,nn,denom
      if(ndims.le.0) then
        pos_str=''
      else
        nn = n-1
        denom = product(dsizes)
        do idim=ndims,1,-1
          denom = denom/dsizes(idim)
          inds(idim) = nn/denom
          nn = nn-denom*inds(idim)
        enddo
        write(pos_str,'(7i6)') 1+inds(1:ndims)
        pos_str='at pos. '//trim(pos_str)
      endif
      return
      end subroutine get_pos_str
