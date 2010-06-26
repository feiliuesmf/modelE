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
      real*4, dimension(:,:), allocatable :: xjl,xjl_hemis,lats_dg
      character(len=30) :: units
      character(len=80) :: lname,title,outfile
      character(len=40) :: vname,vname_hemis,vname_vmean
      character(len=8) :: tpow
      character(len=4) :: dash='----'
      character(len=20) :: latname
      real*4 :: prtfac,fglob,fnh,fsh
      integer :: j,l,jm,lm,inc,lstr,prtpow,linect,nargs,k1,k2,lunit
      integer :: lats_per_zone,j1,j2,zone,nzones,lats_this_zone
      integer :: minj,maxj
      logical :: do_giss,all_lats
      real*4, parameter :: missing=-1.e30
      integer :: vid

c
c Various string formats
c
      character(len=80), parameter :: fmtlat =
     &     "('  P(MB)   MEAN G      NH      SH  ',32I4)"
      integer :: status,varid,varid_hemis,varid_vmean,nvars,dimids(2),
     &     plm_dimid,ple_dimid, lat_dimid,lat2_dimid
      character(len=132) :: xlabel
      character(len=100) :: fromto

      nargs = iargc()

c
c Look at command-line arguments to see whether GISS-convention
c fortran-format binary output was requested
c
      k1 = index(progargs,'gissfmt=')
      do_giss = k1.gt.0
      if(do_giss) then
        do k2=k1+8,len_trim(progargs)
          if(progargs(k2:k2).eq.' ') exit
        enddo
        outfile=progargs(k1+8:k2-1)
        lunit = 10
        open(lunit,file=trim(outfile),form='unformatted',
     &       convert='big_endian')
      endif

c
c Check whether printout of all lats was requested
c
      all_lats = index(progargs,'all_lats').gt.0

c
c get run ID, time/date info, number of latitudes and levels
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)

      latname = 'lat_budg'
      lat2_dimid = -99
      if(nf_inq_dimid(fid,trim(latname),lat_dimid) .ne. nf_noerr) then
        latname = 'lat'
        status = nf_inq_dimid(fid,trim(latname),lat_dimid)
        status = nf_inq_dimid(fid,'lat2',lat2_dimid)
      endif
      call get_dimsize(fid,trim(latname),jm)
      call get_dimsize(fid,'plm',lm)


c
c allocate workspace
c
      allocate(lat_dg(jm),lats_dg(jm,2),
     &     vmean(jm+3),xjl(jm,lm),xjl_hemis(3,lm))
      allocate(plm(lm),ple(lm),pm(lm))

c
c read geometry
c
      call get_var_real(fid,trim(latname),lats_dg(1,1))
      if(lat2_dimid.gt.0) then ! get secondary lats
        call get_var_real(fid,'lat2',lats_dg(1,2))
      endif
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
      nzones = 1
      lats_per_zone = jm
      if(all_lats) then
        do while(lats_per_zone .gt. 32)
          nzones = nzones + 1
          lats_per_zone = jm/nzones
        enddo
        if(lats_per_zone*nzones.ne.jm) stop 'factoring error'
        inc=1
      else
        inc=1+(jm-1)/24
      endif

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
        xjl_hemis = missing
        status = nf_get_var_real(fid,varid_hemis,xjl_hemis)
        vmean = missing
        status = nf_get_var_real(fid,varid_vmean,vmean)
        status = nf_get_var_real(fid,varid,xjl)

c
c form title string and rescale fields for ASCII output
c
        if(prtpow.ne.0) then
          prtfac = 10.**(-prtpow)
          where(xjl.ne.missing) xjl = xjl*prtfac
          where(xjl_hemis.ne.missing) xjl_hemis = xjl_hemis*prtfac
          where(vmean.ne.missing) vmean = vmean*prtfac
          write (tpow, '(i3)') prtpow
          tpow='10**'//trim(adjustl(tpow))
          units = trim(tpow)//' '//trim(units)
        endif
        title = trim(lname)//' ('//trim(units)//')'

c
c retrieve horizontal and vertical coordinate info
c
        status = nf_inq_vardimid(fid,varid,dimids)
        if(dimids(1).eq.lat_dimid) then
          lat_dg(:) = lats_dg(:,1)
          minj = 1; maxj = jm
        else
          lat_dg(:) = lats_dg(:,2)
          if(lat_dg(1).le.-90.) then
            minj = 2; maxj = jm
          else
            minj = 1; maxj = jm-1
          endif
        endif
        if(dimids(2).eq.plm_dimid) then
          pm(:) = plm(:)
        else
          pm(:) = ple(:)
        endif

c
c write binary output
c
        if(do_giss) then
          write(lunit) title,jm,lm,1,1,
     &         xjl,
     &         lat_dg,pm,
     &         real(1.,kind=4),real(1.,kind=4),
     &         'LATITUDE        ',
     &         'PRESSURE (MB)   ',
     &         '                ','                ','NASAGISS',
     &         vmean,
     &         xjl_hemis
        endif

c
c print table
c
        where(xjl.eq.missing) xjl=0.
        where(xjl_hemis.eq.missing) xjl_hemis=0.
        where(vmean.eq.missing) vmean=0.
        linect = linect + lm + 7
        if(linect.gt.60) then
          write(6,'(a)') '1'//xlabel
          write(6,'(a)') ' '//fromto
          linect = lm+8
        endif
        do zone=1,nzones
          j2 = min(maxj, jm -(zone-1)*lats_per_zone)
          j1 = max(minj, j2 -lats_per_zone +1)
          lats_this_zone = 1+(j2-j1)/inc
          write(6,901) title
          call prtdashes(lats_this_zone)
          write(6,fmtlat) int(lat_dg(j2:j1:-inc))
          call prtdashes(lats_this_zone)
          do l=lm,1,-1
            fsh  = xjl_hemis(1,l)
            fnh  = xjl_hemis(2,l)
            fglob= xjl_hemis(3,l)
            write(6,902) pm(l),fglob,fnh,fsh,
     &           (nint(xjl(j,l)),j=j2,j1,-inc)
          enddo
          call prtdashes(lats_this_zone)
          fsh  =vmean(jm+1)
          fnh  =vmean(jm+2)
          fglob=vmean(jm+3)
          write(6,903) '    ',fglob,fnh,fsh,
     &         (nint(vmean(j)),j=j2,j1,-inc)
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(lat_dg,lats_dg,vmean,xjl,xjl_hemis)
      deallocate(plm,ple,pm)

c
c close binary output file
c
      if(do_giss) close(lunit)

      return
  901 FORMAT ('0',30X,A64)
  902 FORMAT (1X,F8.3,3F8.1,1X,32I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,32I4)
      end subroutine prtajl

      subroutine prtdashes(n)
      implicit none
      integer :: n,nn
      write(6,'(a2)',advance='NO') '  '
      do nn=1,8+n-1
        write(6,'(a4)',advance='NO') '----'
      enddo
      write(6,'(a4)') '----'
      return
      end subroutine prtdashes
