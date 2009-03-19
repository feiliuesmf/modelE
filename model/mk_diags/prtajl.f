      subroutine prtajl(fid,progargs)
!@sum prtajl prints the fields in an input file whose
!@+   metadata mark them as having been created from
!@+   a modelE AJL-type array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:), allocatable :: lat_dg,vmean,plm,ple,pm
      real*4, dimension(:,:), allocatable :: xjl,xjl_hemis
      character(len=30) :: units
      character(len=80) :: lname,title
      character(len=40) :: vname,vname_hemis,vname_vmean
      character(len=8) :: tpow
      character(len=4) :: dash='----'
      real*4 :: prtfac,fglob,fnh,fsh
      integer :: j,l,jm,lm,inc,lstr,prtpow,linect
c
c Various string formats
c
      character(len=80), parameter ::
     &     fmtlat = "('  P(MB)   MEAN G      NH      SH  ',24I4)"

      integer :: status,varid,varid_hemis,varid_vmean,nvars,dimids(2),
     &     plm_dimid,ple_dimid
      character(len=132) :: xlabel
      character(len=100) :: fromto

c
c get run ID, time/date info, number of latitudes and levels
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'lat_budg',jm)
      call get_dimsize(fid,'plm',lm)

c
c allocate workspace
c
      allocate(lat_dg(jm),vmean(jm+3),xjl(jm,lm),xjl_hemis(3,lm))
      allocate(plm(lm),ple(lm),pm(lm))

c
c read geometry
c
      call get_var_real(fid,'lat_budg',lat_dg)
      call get_var_real(fid,'plm',plm)
      call get_var_real(fid,'ple',ple)
      status = nf_inq_dimid(fid,'plm',plm_dimid)
      status = nf_inq_dimid(fid,'ple',ple_dimid)

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c Loop over quantities.  An array xyz is printed if there also
c exists an array xyz_hemis containing hemispheric/global averages
c for that quantity.
c
      inc=1+(jm-1)/24
      linect=65

      do varid_hemis=1,nvars
        status = nf_inq_varname(fid,varid_hemis,vname_hemis)
        lstr = len_trim(vname_hemis)
        if(vname_hemis(lstr-5:lstr).ne.'_hemis') cycle
        vname = vname_hemis(1:lstr-6)
        vname_vmean = trim(vname)//'_vmean'
        status = nf_inq_varid(fid,trim(vname),varid)
        status = nf_inq_varid(fid,trim(vname_vmean),varid_vmean)
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        if(trim(units).eq.'unused') cycle
        lname = ''
        status = nf_get_att_text(fid,varid,'long_name',lname)
        prtpow = 0
        status = nf_get_att_int(fid,varid,'prtpow',prtpow)
        status = nf_get_var_real(fid,varid_hemis,xjl_hemis)
        status = nf_get_var_real(fid,varid_vmean,vmean)
        status = nf_get_var_real(fid,varid,xjl)
        where(xjl.eq.-1.e30) xjl=0.
        where(xjl_hemis.eq.-1.e30) xjl_hemis=0.
        where(vmean.eq.-1.e30) vmean=0.

c
c form title string and rescale fields for ASCII output
c
        if(prtpow.ne.0) then
          prtfac = 10.**(-prtpow)
          xjl = xjl*prtfac
          xjl_hemis = xjl_hemis*prtfac
          vmean = vmean*prtfac
          write (tpow, '(i3)') prtpow
          tpow='10**'//trim(adjustl(tpow))
          units = trim(tpow)//' '//trim(units)
        endif
        title = trim(lname)//' ('//trim(units)//')'

c
c retrieve vertical coordinate info
c
        status = nf_inq_vardimid(fid,varid,dimids)
        if(dimids(2).eq.plm_dimid) then
          pm(:) = plm(:)
        else
          pm(:) = ple(:)
        endif

c
c print table
c
        linect = linect + lm + 7
        if(linect.gt.60) then
          write(6,'(a)') xlabel
          write(6,'(a)') fromto
          linect = lm+8
        endif
        write(6,901) title,(dash,j=1,jm,inc)
        write(6,fmtlat) nint(lat_dg(jm:inc:-inc))
        write(6,905)
        do l=lm,1,-1
          fsh  = xjl_hemis(1,l)
          fnh  = xjl_hemis(2,l)
          fglob= xjl_hemis(3,l)
          write(6,902) pm(l),fglob,fnh,fsh,
     &         (nint(xjl(j,l)),j=jm,inc,-inc)
        enddo
        write(6,905)
        fsh  =vmean(jm+1)
        fnh  =vmean(jm+2)
        fglob=vmean(jm+3)
        write(6,903) '    ',fglob,fnh,fsh,(nint(vmean(j)),j=jm,inc,-inc)
      enddo

c
c deallocate workspace
c
      deallocate(lat_dg,vmean,xjl,xjl_hemis)
      deallocate(plm,ple,pm)

      return
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  905 FORMAT (2X,124('-'))
      end subroutine prtajl
