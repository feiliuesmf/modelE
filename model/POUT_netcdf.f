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
      implicit none

      include '/usr/local/netcdf-3.4/include/netcdf.inc'

      private
      public ::
     &     outfile
     &    ,open_out,def_dim_out,set_dim_out,close_out
     &    ,status_out,varid_out,out_fid
     &    ,ndims_out,dimids_out,file_dimlens
     &    ,units,long_name,missing

      character(len=80) :: outfile
      integer :: status_out,out_fid
      integer :: ndims_out
      integer :: varid_out
      integer, dimension(7) :: dimids_out
      character(len=20) :: units=''
      character(len=80) :: long_name=''

      integer, parameter :: ndfmax=100 ! max dims in file
      integer :: ndims_file=0
      integer, dimension(ndfmax) :: file_dimids
      integer, dimension(ndfmax) :: file_dimlens
      character(len=20), dimension(ndfmax) :: file_dimnames

      real :: missing=-1.e30

      contains
      
      subroutine open_out
      status_out = nf_create (trim(outfile), nf_clobber, out_fid)
      if(status_out.ne.nf_noerr) then
         write(6,*) 'cannot create: '//trim(outfile)
         stop
      endif
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
      character(len=20) :: var_name
      double precision :: var(1)
      integer :: var_nelems,n
      status_out = nf_redef(out_fid)
      status_out=nf_def_var(out_fid,trim(var_name),nf_real,
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
      status_out = nf_enddef(out_fid)
      status_out = nf_put_var_double(out_fid,varid_out,var)
      units=''
      long_name=''
      return
      end subroutine wrtarr

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
      character(len=20) :: var_name,dim_name

      outfile = filename

! define output file
      call open_out
      iu_ij = out_fid

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtarr(var_name,lon_dg(1,1))
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lat_dg(1,1))

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

      subroutine POUT_IJ(TITLE,XIJ,XJ,XSUM)
!@sum  POUT_IJ output lat-lon binary records
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : IM,JM
      USE DAGCOM, only : iu_ij
      USE NCOUT
      IMPLICIT NONE
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var XIJ lat/lon output field 
      REAL*8, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field 
      REAL*8, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field 
      REAL*8, INTENT(IN) :: XSUM
      integer :: nvars, status
      character(len=20) :: var_name,dim_name

! set shape of output array
      ndims_out = 2
      dim_name='longitude'; call set_dim_out(dim_name,1)
      dim_name='latitude'; call set_dim_out(dim_name,2)

C numerical naming convention for now
      status = nf_inq_nvars(out_fid, nvars)
      var_name=' '
      write(var_name,'(a3,i3.3)') 'IJ_',nvars
      long_name=title
      call wrtarr(var_name,xij)
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
      character(len=20) :: var_name,dim_name

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
      USE GEOM, only : lat_dg
      USE DAGCOM, only : lm_req,iu_jl
      USE DAGPCOM, only : plm,ple,ple_dn
      USE NCOUT
      IMPLICIT NONE
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
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

      integer :: nvars, status
      character(len=20) :: var_name,dim_name
      
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

C numerical naming convention for now
      status = nf_inq_nvars(out_fid, nvars)
      var_name=' '
      write(var_name,'(a3,i3.3)') 'JL_',nvars
      long_name=title

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
      character(len=20) :: var_name,dim_name

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
      IMPLICIT NONE
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
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

      integer :: nvars, status
      character(len=20) :: var_name,dim_name
      
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

C numerical naming convention for now
      status = nf_inq_nvars(out_fid, nvars)
      var_name=' '
      write(var_name,'(a3,i3.3)') 'IL_',nvars
      long_name=title

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
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      character(len=20) :: var_name,dim_name

      outfile = filename
      call open_out
      iu_j = out_fid

      dim_name='latitude'; call def_dim_out(dim_name,jm)

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
      IMPLICIT NONE
      CHARACTER*16, DIMENSION(KAJ),INTENT(INOUT) :: TITLE
      CHARACTER*16, INTENT(IN) :: TERRAIN
      REAL*8, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
      INTEGER, INTENT(IN) :: KMAX
      INTEGER K,N,J

      character(len=20) :: dim_name
      character(len=32) :: var_name
      character(len=16) :: terr_name

! set shape of output array
      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)

C**** remove blanks, parentheses from TERRAIN
      terr_name=terrain
      if (terr_name(16:16).eq.'(' .or. terr_name(16:16).eq.')')
     &     terr_name(16:16)=' '
      do n=15,1,-1
         if (terr_name(n:n).eq.'(' .or. terr_name(n:n).eq.')')
     &        terr_name(n:15)=terr_name(n+1:16)
      end do
      do n=1,len_trim(terr_name)
         if(terr_name(n:n).eq.' ') cycle
         terr_name(1:16-n+1)=terr_name(n:16)
         exit
      enddo
      do n=1,len_trim(terr_name)
         if (terr_name(n:n).eq.' ') terr_name(n:n)='_'
      end do

      DO K=1,KMAX
C numerical naming convention for now
        var_name=' '
        write(var_name,'(a2,i3.3)') 'J_',k
        var_name=trim(terr_name)//'_'//var_name
        long_name(1:16)=title(k)
        call wrtarr(var_name,budg(1,k))
      END DO

      return
      end

