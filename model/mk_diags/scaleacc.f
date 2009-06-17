!@sum scaleacc is a generic scaling routine for modelE acc-files.
!@+   See conventions.txt for documentation.
!@auth M. Kelley
      subroutine scaleacc(fid,accname,ofile_base)
      implicit none
      include 'netcdf.inc'
      integer :: fid                  ! input file ID
      character(len=*) :: accname     ! name of acc-array to scale
      character(len=80) :: ofile_base ! basename of output file
c
      character(len=20) :: dcat
      character(len=80) :: cdlfile
      character(len=100), dimension(:), allocatable :: cdl
      real*4, dimension(:), allocatable :: scale_acc
      integer, dimension(:), allocatable :: ia_acc,denom_acc
      character(len=30), dimension(:), allocatable :: sname_acc
      real*4, dimension(:), allocatable :: xout,xden
      real*4, dimension(:), allocatable :: xout_hemis,xden_hemis
      real*4, dimension(:), allocatable :: xout_vmean,xden_vmean
      integer :: idacc(12)
      integer :: k,kd,kacc,kcdl,kgw,kend,arrsize,ndims,ndimsh,sdim,jdim
      integer, dimension(7) :: srt,cnt,accsizes,dimids,hemi_sizes
      integer, dimension(7) :: cnt_hemis,cnt_vmean
      integer :: status,ofid,accid,varid,accid_hemis,jdimid,accid_vmean
      real*4, parameter :: undef=-1.e30
      character(len=132) :: xlabel
      character(len=100) :: fromto
      logical :: do_hemis,do_vmean

      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)

      dcat = trim(accname)

      call handle_err(nf_inq_varid(fid,trim(dcat),accid),
     &     'finding '//trim(dcat)//' in input file')
      call get_vdimsizes(fid,trim(dcat),ndims,accsizes)
      srt(:) = 1
      cnt(1:ndims) = accsizes(1:ndims)

      status = nf_inq_varid(fid,'hemis_'//trim(dcat),accid_hemis)
      do_hemis = status.eq.nf_noerr
      if(do_hemis) then
        call get_vdimsizes(fid,'hemis_'//trim(dcat),ndimsh,hemi_sizes)
        call get_dimsize(fid,'shnhgm',k)
        if(k.ne.3) stop 'bad size for shnhgm dimension'
      endif
      status = nf_inq_varid(fid,'vmean_'//trim(dcat),accid_vmean)
      do_vmean = status.eq.nf_noerr
     &     .and. (index(dcat,'ajl').gt.0 .or. index(dcat,'agc').gt.0)

c
c Find the size of the dimension along which to split the data
c
      status = nf_get_att_int(fid,accid,'split_dim',sdim)
      if(status.eq.nf_noerr) then
        kacc = accsizes(sdim)
        cnt(sdim) = 1
      else
        sdim = 7 ! necessary?
        kacc = 1
      endif

c
c Allocate space for metadata
c
      allocate(ia_acc(kacc),denom_acc(kacc),scale_acc(kacc))
      allocate(sname_acc(kacc))

c
c Calculate the size of one output array, allocate workspace
c
      arrsize = product(cnt(1:ndims))
      allocate(xout(arrsize),xden(arrsize))

c
c Some setup for hemispheric/global means.
c We assume only that shnhgm and the latitude dimension
c have the same index.
c
      if(do_hemis) then
        status = nf_inq_vardimid(fid,accid_hemis,dimids)
        status = nf_inq_dimid(fid,'shnhgm',jdimid)
        do jdim=1,ndimsh
          if(jdimid.eq.dimids(jdim)) exit
        enddo
        cnt_hemis(1:ndimsh) = hemi_sizes(1:ndimsh)
        cnt_hemis(sdim) = 1
        arrsize = product(cnt_hemis(1:ndimsh))
        allocate(xout_hemis(arrsize),xden_hemis(arrsize))
      endif

c
c Setup for vertical means (currently only for jl-type arrays).
c
      do_vmean = do_vmean .and. jdim.eq.1
      if(do_vmean) then
        allocate(xout_vmean(accsizes(1)+3),
     &           xden_vmean(accsizes(1)+3))
        cnt_vmean(:) = cnt(:)
        cnt_vmean(1) = accsizes(1)+3
        cnt_vmean(2) = 1
      endif

c
c Define the output file using the ncgen utility
c
      call get_dimsize(fid,'kcdl_'//trim(dcat),kcdl)      
      allocate(cdl(kcdl))
      cdl = ''
      call get_var_text(fid,'cdl_'//trim(dcat),cdl)
      cdl(1) = 'netcdf '//trim(ofile_base)//' { '
      cdlfile = trim(ofile_base)//'.cdl'
      open(10,file=cdlfile)
      do k=kcdl,1,-1
        if(index(cdl(k),'}').gt.0) then
          kend = k; exit
        endif
      enddo
      kgw = kend
      do k=kend,1,-1
        if(index(cdl(k),'data:').gt.0) then
          kgw = k; exit
        endif
      enddo
      do k=1,kend
        if(k.eq.kgw) then
          write(10,*) '// global attributes:'
          write(10,'(a)') '    :xlabel = "'//trim(xlabel)//'" ;'
          write(10,'(a)') '    :fromto = "'//fromto//'" ;'
        endif
        write(10,'(a)') cdl(k)
      enddo
      close(10)
      call system('ncgen -b '//trim(cdlfile))
      call system('rm -f '//trim(cdlfile))

c
c Read acc metadata needed for scaling
c
      call get_var_int(fid,'idacc',idacc)
      call get_var_real(fid,'scale_'//trim(dcat),scale_acc)
      status = nf_inq_varid(fid,'denom_'//trim(dcat),varid)
      if(status.eq.nf_noerr) then ! this acc array needs denom info
        call get_var_int(fid,'denom_'//trim(dcat),denom_acc)
      else
        denom_acc = 0
      endif
      call get_var_text(fid,'sname_'//trim(dcat),sname_acc)
      status = nf_inq_varid(fid,'ia_'//trim(dcat),varid)
      if(status.eq.nf_noerr) then ! this acc array has idacc-info
        call get_var_int(fid,'ia_'//trim(dcat),ia_acc)
      else ! this acc array has a custom counter.
c Put the counter value into idacc(1)
        call get_var_int(fid,'ntime_'//trim(dcat),idacc(1))
        ia_acc(:) = 1
      endif
c
c open output file
c
      status = nf_open(trim(ofile_base)//'.nc',nf_write,ofid)

c
c copy coordinate and other info from the acc file to the output file
c
      call copy_shared_vars(fid,ofid)

c
c loop over outputs
c
      do k=1,kacc
        status = nf_inq_varid(ofid,trim(sname_acc(k)),varid)
        if(status.ne.nf_noerr) cycle ! this output was not requested

c
c scale this field
c
        srt(sdim) = k
        status = nf_get_vara_real(fid,accid,srt,cnt,xout)
        xout = xout*scale_acc(k)/idacc(ia_acc(k))
        kd = denom_acc(k)
        if(kd.gt.0) then
          srt(sdim) = kd
          status = nf_get_vara_real(fid,accid,srt,cnt,xden)
          where(xden.ne.0.)
            xout = xout*idacc(ia_acc(kd))/xden
          elsewhere
            xout = undef
          end where
        endif

c
c write this field to the output file
c
        call put_var_real(ofid,sname_acc(k),xout)


c
c scale/write the hemispheric/global means of this field if present
c
        if(do_hemis) then
          srt(sdim) = k
          status = nf_get_vara_real(fid,accid_hemis,srt,cnt_hemis,
     &         xout_hemis)
          xout_hemis = xout_hemis*scale_acc(k)/idacc(ia_acc(k))
          if(kd.gt.0) then
            srt(sdim) = kd
            status = nf_get_vara_real(fid,accid_hemis,srt,cnt_hemis,
     &           xden_hemis)
            where(xden_hemis.ne.0.)
              xout_hemis = xout_hemis*idacc(ia_acc(kd))/xden_hemis
            elsewhere
              xout_hemis = undef
            end where
          endif
          call put_var_real(ofid,trim(sname_acc(k))//'_hemis',
     &         xout_hemis)
        endif

c
c scale/write the vertical means of this field if present
c
        if(do_vmean) then
          srt(sdim) = k
          status = nf_get_vara_real(fid,accid_vmean,srt,cnt_vmean,
     &         xout_vmean)
          xout_vmean = xout_vmean*scale_acc(k)/idacc(ia_acc(k))
          if(kd.gt.0) then
            srt(sdim) = kd
            status = nf_get_vara_real(fid,accid_vmean,srt,cnt_vmean,
     &           xden_vmean)
            where(xden_vmean.ne.0.)
              xout_vmean = xout_vmean*idacc(ia_acc(kd))/xden_vmean
            elsewhere
              xout_vmean = undef
            end where
          endif
          status = nf_inq_varid(ofid,trim(sname_acc(k))//'_vmean',varid)
          if(status.eq.nf_noerr) status=nf_put_var_real(ofid,varid,
     &         xout_vmean)
        endif
      enddo

      status = nf_close(ofid)

      return
      end subroutine scaleacc
