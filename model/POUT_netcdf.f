!@sum  POUT_netcdf default output routines for netcdf formats
!@auth M. Kelley
!@ver  1.0

C****
C**** If other formats are desired, please replace these routines
C**** with ones appropriate for your chosen output format, but with the
C**** same interface
C****
C**** Note: it would be nice to amalgamate IL and JL, but that will
C**** have to wait.

      module ncout
      use MODEL_COM, only : xlabel,lrunid
      use DAGCOM, only : acc_period
      implicit none

      include '/usr/local/netcdf-3.4/include/netcdf.inc'

      private
      public ::
     &     outfile
     &    ,open_out,def_dim_out,set_dim_out,close_out
     &    ,status_out,varid_out,out_fid
     &    ,ndims_out,dimids_out,file_dimlens
     &    ,units,long_name,missing,real_att_name,real_att
     &    ,nf_disk_type,nf_arr_type

      character(len=80) :: outfile
      integer :: status_out,out_fid
      integer :: ndims_out
      integer :: varid_out
      integer, dimension(7) :: dimids_out
      character(len=30) :: units=''
      character(len=80) :: long_name=''

      integer, parameter :: ndfmax=100 ! max dims in file
      integer :: ndims_file=0
      integer, dimension(ndfmax) :: file_dimids
      integer, dimension(ndfmax) :: file_dimlens
      character(len=20), dimension(ndfmax) :: file_dimnames

      real, parameter :: missing=-1.e30

      integer, parameter :: real_att_max_size=100
      character(len=30) :: real_att_name=''
      real, dimension(real_att_max_size) :: real_att = missing

c netcdf library will convert prog_dtype to disk_dtype 
c these are the defaults for converting GCM real*8 to real*4 on disk
      integer :: disk_dtype = nf_real
      integer :: prog_dtype = nf_double

      contains
      
      subroutine open_out
      status_out = nf_create (trim(outfile), nf_clobber, out_fid)
      if(status_out.ne.nf_noerr) then
         write(6,*) 'cannot create: '//trim(outfile)
         stop
      endif
c-----------------------------------------------------------------------
c write global header
c-----------------------------------------------------------------------
      status_out=nf_put_att_text(out_fid,nf_global,'run_name',
     &     lrunid,xlabel)
      status_out=nf_put_att_text(out_fid,nf_global,'acc_period',
     &     len_trim(acc_period),acc_period)
      return
      end subroutine open_out

      subroutine def_dim_out(dim_name,dimlen)
      character(len=20) :: dim_name
      integer :: dimlen
      integer :: tmp_id
      if(ndims_file.eq.ndfmax)
     &     stop 'def_dim_out: too many dimensions in file'
      ndims_file = ndims_file + 1
      file_dimlens(ndims_file) = dimlen
      file_dimnames(ndims_file)=dim_name
      status_out = nf_def_dim(out_fid,trim(dim_name),dimlen,tmp_id)
      file_dimids(ndims_file) = tmp_id
      return
      end subroutine def_dim_out

      subroutine set_dim_out(dim_name,dim_num)
      character(len=20) :: dim_name
      integer :: dim_num
      integer :: idim
      if(ndims_file.eq.0) stop 'set_dim_out: no dims defined'
      if(dim_num.gt.ndims_out) stop 'set_dim_out: invalid dim #'
      idim=1
      do while(idim.le.ndims_file)
         if(dim_name.eq.file_dimnames(idim)) then
            dimids_out(dim_num)=file_dimids(idim)
            exit
         endif
         idim = idim + 1
      enddo
      if(idim.gt.ndims_file) stop 'set_dim_out: invalid dim_name'
      return
      end subroutine set_dim_out

      subroutine close_out
      status_out = nf_close(out_fid)
      ndims_file = 0
      return
      end subroutine close_out

      end module ncout

      subroutine wrtarr(var_name,var)
      use ncout
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=30) :: var_name
      double precision :: var(1)
      integer :: var_nelems,n
      status_out = nf_redef(out_fid)
c check to make sure that variable isn't already defined in output file
      if(nf_inq_varid(out_fid,trim(var_name),varid_out)
     &                 .eq.nf_noerr) then
         write(6,*) 'wrtarr: variable already defined: ',trim(var_name)
         stop
      endif
      status_out=nf_def_var(out_fid,trim(var_name),disk_dtype,
     &     ndims_out,dimids_out,varid_out)
      if(len_trim(units).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'units',len_trim(units),units)
      if(len_trim(long_name).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'long_name',len_trim(long_name),long_name)
      if(len_trim(real_att_name).gt.0) status_out =
     &     nf_put_att_real(out_fid,varid_out,trim(real_att_name),
     &     nf_float,count(real_att.ne.missing),real_att)
c define missing value attribute only if there are missing values
      var_nelems = product(file_dimlens(dimids_out(1:ndims_out)))
      do n=1,var_nelems
         if(var(n).eq.missing) then
            status_out = nf_put_att_real(out_fid,varid_out,
     &           'missing_value',nf_float,1,missing)
            exit
         endif
      enddo
      status_out = nf_enddef(out_fid)
      select case (prog_dtype)
      case (nf_double)
         status_out = nf_put_var_double(out_fid,varid_out,var)
      case (nf_real  )
         status_out = nf_put_var_real  (out_fid,varid_out,var)
      case (nf_int   )
         status_out = nf_put_var_int   (out_fid,varid_out,var)
      case (nf_char  )
         status_out = nf_put_var_text  (out_fid,varid_out,var)
      end select
c restore defaults
      units=''
      long_name=''
      real_att_name=''
      real_att=missing
      disk_dtype = nf_real
      prog_dtype = nf_double
      return
      end subroutine wrtarr

      subroutine wrtarrn(var_name,var,outer_pos)
c this subroutine is like wrtarr, but writes the outer_pos-th element
c of an arbitrarily shaped array.  when outer_pos=1, it defines the
c array as well.  the first call to it for a particular array MUST have
c outer_pos=1, otherwise it stops.
c wrtarrn does not currently have the capability to write any real_atts
      use ncout
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=30) :: var_name
      double precision :: var(1)
      integer :: outer_pos
      integer :: var_nelems,n
      integer, dimension(7) :: srt,cnt

      if(ndims_out.lt.2) stop 'wrtarrn: invalid ndims_out'
      if(outer_pos.le.0 .or.
     &     outer_pos.gt.file_dimlens(dimids_out(ndims_out))) then
         write(6,*) 'wrtarrn: invalid outer_pos'
         stop
      endif

      status_out = nf_redef(out_fid)

c define the variable if outer_pos is 1
      if(outer_pos.eq.1) then

c check to make sure that variable isn't already defined in output file
      if(nf_inq_varid(out_fid,trim(var_name),varid_out)
     &                 .eq.nf_noerr) then
         write(6,*) 'wrtarrn: variable already defined: ',trim(var_name)
         stop
      endif
      status_out=nf_def_var(out_fid,trim(var_name),disk_dtype,
     &     ndims_out,dimids_out,varid_out)
      if(len_trim(units).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'units',len_trim(units),units)
      if(len_trim(long_name).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'long_name',len_trim(long_name),long_name)
c define missing value attribute only if there are missing values
      var_nelems = product(file_dimlens(dimids_out(1:ndims_out)))
      do n=1,var_nelems
         if(var(n).eq.missing) then
            status_out = nf_put_att_real(out_fid,varid_out,
     &           'missing_value',nf_float,1,missing)
            exit
         endif
      enddo

      else
c check to make sure that variable is already defined in output file
      if(nf_inq_varid(out_fid,trim(var_name),varid_out)
     &                 .ne.nf_noerr) then
         write(6,*) 'wrtarrn: variable not yet defined: ',trim(var_name)
         stop
      endif

      endif

      status_out = nf_enddef(out_fid)

      srt(1:ndims_out-1) = 1
      srt(ndims_out) = outer_pos
      do n=1,ndims_out-1
         cnt(n) = file_dimlens(dimids_out(n))
      enddo
      cnt(ndims_out) = 1

      select case (prog_dtype)
      case (nf_double)
         status_out = nf_put_vara_double(out_fid,varid_out,srt,cnt,var)
      case (nf_real  )
         status_out = nf_put_vara_real  (out_fid,varid_out,srt,cnt,var)
      case (nf_int   )
         status_out = nf_put_vara_int   (out_fid,varid_out,srt,cnt,var)
      case (nf_char  )
         status_out = nf_put_vara_text  (out_fid,varid_out,srt,cnt,var)
      end select
c restore defaults
      units=''
      long_name=''
      disk_dtype = nf_real
      prog_dtype = nf_double
      return
      end subroutine wrtarrn

      subroutine open_ij(filename)
!@sum  OPEN_IJ opens the lat-lon binary output file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE GEOM, only : lon_dg,lat_dg
      USE DAGCOM, only : iu_ij
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!
      character(len=30) :: var_name,dim_name

      outfile = filename

! define output file
      call open_out
      iu_ij = out_fid

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='lonb'; call def_dim_out(dim_name,im)
      dim_name='latb'; call def_dim_out(dim_name,jm-1)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtarr(var_name,lon_dg(1,1))
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lat_dg(1,1))
      dim_name='lonb'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='lonb'; call wrtarr(var_name,lon_dg(1,2))
      dim_name='latb'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latb'; call wrtarr(var_name,lat_dg(2,2))

      return
      end subroutine open_ij

      subroutine close_ij
!@sum  OPEN_IJ closes the lat-lon binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_ij
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_ij
      call close_out

      return
      end subroutine close_ij

      subroutine POUT_IJ(TITLE,XIJ,XJ,XSUM,IJGRID)
!@sum  POUT_IJ output lat-lon binary records
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : IM,JM
      USE DAGCOM, only : iu_ij
      USE NCOUT
      USE BDIJ, title_ij=>title,units_ij=>units,
     &          lname_ij=>lname,sname_ij=>sname
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var XIJ lat/lon output field 
      REAL*8, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field 
      REAL*8, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field 
      REAL*8, INTENT(IN) :: XSUM
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID

      integer :: nvars, status
      character(len=30) :: var_name,lon_name,lat_name

! (re)set shape of output array
      ndims_out = 2
      if(ijgrid.eq.1) then
         lon_name='longitude'
         lat_name='latitude'
      else if(ijgrid.eq.2) then
         lon_name='lonb'
         lat_name='latb'
      else
         stop 'pout_ij: unrecognized grid'
      endif
      call set_dim_out(lon_name,1)
      call set_dim_out(lat_name,2)

      var_name=sname_ij(nt_ij)
      long_name=lname_ij(nt_ij)
      units=units_ij(nt_ij)
      real_att_name='glb_mean'
      real_att(1)=xsum
      if(ijgrid.eq.1) then
         call wrtarr(var_name,xij(1,1))
      else
         call wrtarr(var_name,xij(1,2)) ! ignore first latitude row
      endif
      return
      end

      subroutine open_jl(filename)
!@sum  OPEN_JL opens the lat-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : jm,lm, sig,sige
      USE DAGCOM, only : iu_jl,lm_req
      use DAGPCOM, only : plm,ple,ple_dn
      USE GEOM, only : lat_dg
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      character(len=30) :: var_name,dim_name

      outfile = filename
      call open_out
      iu_jl = out_fid

      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='latb'; call def_dim_out(dim_name,jm-1)
      dim_name='p'; call def_dim_out(dim_name,lm)
      dim_name='prqt'; call def_dim_out(dim_name,lm+lm_req)
      dim_name='ple_up'; call def_dim_out(dim_name,lm)
      dim_name='ple_dn'; call def_dim_out(dim_name,lm)
      dim_name='ple_int'; call def_dim_out(dim_name,lm-1)

! put lat,ht into output file
      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude'; call wrtarr(var_name,lat_dg(1,1))
      dim_name='latb'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latb'; call wrtarr(var_name,lat_dg(2,2))
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='p'; call wrtarr(var_name,plm)
      dim_name='prqt'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='prqt'; call wrtarr(var_name,plm)
      dim_name='ple_up'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple_up'; call wrtarr(var_name,ple)
      dim_name='ple_dn'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple_dn'; call wrtarr(var_name,ple_dn)
      dim_name='ple_int'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple_int'; call wrtarr(var_name,ple)

      return
      end subroutine open_jl

      subroutine close_jl
!@sum  OPEN_JL closes the lat-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_jl
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_jl
      call close_out

      return
      end subroutine close_jl

      subroutine POUT_JL(TITLE,J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_JL output lat-height binary records
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : JM,LM,LS1
      USE DAGCOM, only : lm_req,iu_jl
      USE DAGPCOM, only : plm,ple,ple_dn
      USE NCOUT
      USE BDJL, title_jl=>title,units_jl=>units,
     &          lname_jl=>lname,sname_jl=>sname
      USE BDJK, title_jk=>title,units_jk=>units,
     &          lname_jk=>lname,sname_jk=>sname
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var J1 minimum j value to output (needed for secondary grid fields)
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field 
!@+       (J1:JM,1:KLMAX) is field
!@+       (JM+1:JM+3,1:KLMAX) are global/NH/SH average over L
!@+       (J1:JM+3,LM+LM_REQ+1) are averages over J
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
      REAL*8, DIMENSION(JM,LM+LM_REQ) :: XJL0
      REAL*8, DIMENSION(JM-1,LM+LM_REQ) :: XJL0B
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM 

      character(len=30) :: var_name,dim_name
      
      CHARACTER*16, INTENT(IN) :: CX,CY
      INTEGER J,L

! (re)set shape of output arrays
      ndims_out = 2

      if(j1.eq.1) then
         dim_name='latitude'
      else if(j1.eq.2) then
         dim_name='latb'
      else
         stop 'pout_jl: unrecognized latitude grid'
      endif
      call set_dim_out(dim_name,1)

      if(klmax.eq.lm+lm_req) then
         dim_name='prqt'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.plm(1:lm))) then
         dim_name='p'
      else if(klmax.eq.ls1-1 .and.all(pm(1:ls1-1).eq.plm(1:ls1-1))) then
! some arrays only defined in troposphere, extend to stratosphere anyway
         dim_name='p'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple(1:lm))) then
         dim_name='ple_up'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple_dn(1:lm))) then
         dim_name='ple_dn'
      else if(klmax.eq.lm-1 .and. all(pm(1:lm-1).eq.ple(1:lm-1))) then
         dim_name='ple_int'
      else
         write(6,*) 'klmax =',klmax,title,pm(1:klmax)
         stop 'pout_jl: unrecognized vertical grid'
      endif
      call set_dim_out(dim_name,2)

      if(in_jlmap) then
         var_name=sname_jl(nt_jl)
         long_name=lname_jl(nt_jl)
         units=units_jl(nt_jl)
      else if(in_jkmap) then
         var_name=sname_jk(nt_jk)
         long_name=lname_jk(nt_jk)
         units=units_jk(nt_jk)
      else
         stop 'pout_jl: unknown output list'
      endif
      real_att_name='g-nh-sh_sums-means'
      real_att(1:3)=XJL(JM+1:JM+3,LM+LM_REQ+1)

      if(j1.eq.1) then
         xjl0(1:jm,1:klmax) = xjl(1:jm,1:klmax)
         call wrtarr(var_name,xjl0)
      else
         xjl0b(1:jm-1,1:klmax) = xjl(2:jm,1:klmax)
         call wrtarr(var_name,xjl0b)
      endif

      return
      end

      subroutine open_il(filename)
!@sum  OPEN_IL opens the lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : im,lm, sig,sige
      USE DAGCOM, only : iu_il,lm_req
      use DAGPCOM, only : plm,ple,ple_dn
      USE GEOM, only : lon_dg
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      character(len=30) :: var_name,dim_name

      outfile = filename
      call open_out
      iu_il = out_fid

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='lonb'; call def_dim_out(dim_name,im)
      dim_name='p'; call def_dim_out(dim_name,lm)
      dim_name='prqt'; call def_dim_out(dim_name,lm+lm_req)
      dim_name='ple_up'; call def_dim_out(dim_name,lm)
      dim_name='ple_dn'; call def_dim_out(dim_name,lm)
      dim_name='ple_int'; call def_dim_out(dim_name,lm-1)

! put lon,ht into output file
      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude'; call wrtarr(var_name,lon_dg(1,1))
      dim_name='lonb'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='lonb'; call wrtarr(var_name,lon_dg(1,2))
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='p'; call wrtarr(var_name,plm)
      dim_name='prqt'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='prqt'; call wrtarr(var_name,plm)
      dim_name='ple_up'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple_up'; call wrtarr(var_name,ple)
      dim_name='ple_dn'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple_dn'; call wrtarr(var_name,ple_dn)
      dim_name='ple_int'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple_int'; call wrtarr(var_name,ple)

      return
      end subroutine open_il

      subroutine close_il
!@sum  OPEN_IL closes the lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_il
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_il
      call close_out

      return
      end subroutine close_il

      subroutine POUT_IL(TITLE,ISHIFT,KLMAX,XIL,PM,CX,CY,
     *     ASUM,GSUM,ZONAL)
!@sum  POUT_IL output lon-height binary records
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : IM,LM,LS1
      USE GEOM, only : lon_dg
      USE DAGCOM, only : lm_req,iu_il
      USE DAGPCOM, only : plm,ple,ple_dn
      USE NCOUT
      USE BDIL, title_il=>title,units_il=>units,
     &          lname_il=>lname,sname_il=>sname
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var ISHIFT flag for secondary grid
      INTEGER, INTENT(IN) :: KLMAX,ISHIFT
!@var XIL output field 
      REAL*8, DIMENSION(IM,LM+LM_REQ+1), INTENT(IN) :: XIL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM 
!@var ASUM vertical mean/sum
      REAL*8, DIMENSION(IM), INTENT(IN) :: ASUM 
!@var GSUM total sum/mean 
      REAL*8, INTENT(IN) :: GSUM
!@var ZONAL zonal sum/mean 
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: ZONAL 
      
      CHARACTER*16, INTENT(IN) :: CX,CY
      INTEGER I,L

      character(len=30) :: var_name,dim_name
      
! (re)set shape of output arrays
      ndims_out = 2

      if(ishift.eq.1) then
         dim_name='longitude'
      else if(ishift.eq.2) then
         dim_name='lonb'
      else
         stop 'pout_il: unrecognized longitude grid'
      endif
      call set_dim_out(dim_name,1)

      if(klmax.eq.lm+lm_req) then
         dim_name='prqt'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.plm(1:lm))) then
         dim_name='p'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple(1:lm))) then
         dim_name='ple_up'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple_dn(1:lm))) then
         dim_name='ple_dn'
      else if(klmax.eq.lm-1 .and. all(pm(1:lm-1).eq.ple(1:lm-1))) then
         dim_name='ple_int'
      else
         write(6,*) 'klmax =',klmax,title,pm(1:klmax)
         stop 'pout_il: unrecognized vertical grid'
      endif
      call set_dim_out(dim_name,2)

      var_name=sname_il(nt_il)
      long_name=lname_il(nt_il)
      units=units_il(nt_il)
      real_att_name='mean'
      real_att(1)=gsum

      call wrtarr(var_name,xil)

      return
      end

      subroutine open_j(filename)
!@sum  OPEN_J opens the latitudinal binary output file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : jm
      USE GEOM, only : lat_dg
      USE DAGCOM, only : iu_j
      USE NCOUT
      USE BDJ, only : nstype_out
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      character(len=30) :: var_name,dim_name

      outfile = filename
      call open_out
      iu_j = out_fid

      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='stype'; call def_dim_out(dim_name,nstype_out)
      dim_name='stype_clen'; call def_dim_out(dim_name,16)

      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lat_dg(1,1))

      return
      end subroutine open_j

      subroutine close_j
!@sum  OPEN_J closes the latitudinal binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_j
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_j
      call close_out

      return
      end subroutine close_j

      subroutine POUT_J(TITLE,BUDG,KMAX,TERRAIN)
!@sum  POUT_J output zonal budget file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : JM
      USE DAGCOM, only : KAJ,iu_j
      USE GEOM, only : lat_dg
      USE NCOUT
      USE BDJ, units_j=>units,
     &         lname_j=>lname,sname_j=>sname
      IMPLICIT NONE
c temporary, to see nf_char
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
c
      CHARACTER*16, DIMENSION(KAJ),INTENT(INOUT) :: TITLE
      CHARACTER*16, INTENT(IN) :: TERRAIN
      REAL*8, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
      INTEGER, INTENT(IN) :: KMAX
      INTEGER K,N,J

      character(len=30) :: var_name,dim_name

! write out surface type name
      ndims_out = 2
      dim_name='stype_clen'; call set_dim_out(dim_name,1)
      dim_name='stype'; call set_dim_out(dim_name,2)
      disk_dtype = nf_char
      prog_dtype = nf_char
      var_name='stype_names';call wrtarrn(var_name,terrain,iotype)


! (re)set shape of output array
      ndims_out = 2
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='stype'; call set_dim_out(dim_name,2)

      DO K=1,KMAX
        var_name=sname_j(nt_j(k))
        long_name=lname_j(nt_j(k))
        units=units_j(nt_j(k))
        call wrtarrn(var_name,budg(1,k),iotype)
      END DO

      return
      end

      subroutine open_ijk(filename)
!@sum  OPEN_IJK opens the lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE GEOM, only : lon_dg,lat_dg
      USE DAGCOM, only : iu_ijk
      use DAGPCOM, only : plm
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!
      character(len=30) :: var_name,dim_name

      outfile = filename

! define output file
      call open_out
      iu_ijk = out_fid

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm-1)
      dim_name='p'; call def_dim_out(dim_name,lm)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtarr(var_name,lon_dg(1,2))
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lat_dg(2,2))
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='p'; call wrtarr(var_name,plm)

      return
      end subroutine open_ijk

      subroutine close_ijk
!@sum  OPEN_IJK closes the lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_ijk
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_ijk
      call close_out

      return
      end subroutine close_ijk

      subroutine POUT_IJK(TITLE,XIJK,XJK,XK)
!@sum  POUT_IJK output lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : IM,JM,LM
      USE DAGCOM, only : iu_ijk
      USE NCOUT
      USE BDIJK, title_ijk=>title,units_ijk=>units,
     &           lname_ijk=>lname,sname_ijk=>sname
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var XIJK lat/lon/height output field 
      REAL*8, DIMENSION(IM,JM-1,LM), INTENT(IN) :: XIJK
!@var XJK lat sum/mean of output field 
      REAL*8, DIMENSION(JM-1,LM), INTENT(IN) :: XJK
!@var XK global sum/mean of output field 
      REAL*8, DIMENSION(LM), INTENT(IN) :: XK

      integer :: nvars, status
      character(len=30) :: var_name,dim_name

! (re)set shape of output array
      ndims_out = 3

      dim_name = 'longitude'
      call set_dim_out(dim_name,1)
      dim_name = 'latitude'
      call set_dim_out(dim_name,2)
      dim_name = 'p'
      call set_dim_out(dim_name,3)

      var_name=sname_ijk(nt_ijk)
      long_name=lname_ijk(nt_ijk)
      units=units_ijk(nt_ijk)

      call wrtarr(var_name,xijk)
      return
      end
