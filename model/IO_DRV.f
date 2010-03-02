#include "rundeck_opts.h"

      SUBROUTINE io_rsf(fname,it,iaction,ioerr)
!@sum  io_rsf manages the reading/writing of restart and acc files
!@auth M. Kelley
!@ver  beta.

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
      USE DOMAIN_DECOMP_ATM, only : grid
      USE MODEL_COM, only : ioread_single,iowrite_single,irerun,
     *                      ioread,iowrite,iowrite_mon
     &     ,itimei,rsf_file_name
      use pario, only : par_open,par_close
      IMPLICIT NONE
!@var fname name of file to be read or written
      character(len=*) :: fname
!@var iaction flag for reading or writing rsf file
      INTEGER, INTENT(IN) :: iaction
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
      integer :: fid,k,iorw
      logical :: do_io_prog,do_io_acc,do_io_longacc
      character(len=200) :: tmpname

      call set_ioptrs_acc_default

      if(iaction.eq.iowrite_single) then
        call calc_derived_acc
        call set_ioptrs_acc_extended
      endif

      do_io_prog = .true.
      if(iaction.eq.iowrite_single) do_io_prog = .false.
      if(iaction.eq.ioread_single)  do_io_prog = .false.

      do_io_acc = .false.
c this logic would be much simpler if there were fewer
c iaction possibilities related to reruns.  reading
c arrays by name will eliminate a number of the
c rerun cases.
      if(iaction.eq.ioread)         do_io_acc = .true.
      if(iaction.eq.iowrite)        do_io_acc = .true.
      if(iaction.eq.iowrite_single) do_io_acc = .true.
      if(iaction.eq.ioread_single)  do_io_acc = .true.

      do_io_longacc = do_io_acc
      if(iaction.eq.iowrite_mon)    do_io_longacc = .true.
      if(iaction.eq.irerun)         do_io_longacc = .true.

c for routines designed to accept only iowrite or ioread:
      if(iaction.le.iowrite) then
        iorw = iowrite
      else
        iorw = ioread
      endif

      ioerr=-1

c
c create/define file contents if necessary
c
      if(iaction.eq.iowrite_mon .or. iaction.eq.iowrite_single) then
        fid = par_open(grid,trim(fname)//'.nc','create')
        if(iaction.eq.iowrite_mon) then
c end-of-month restart file, no acc arrays
          call def_rsf_label(fid)
          call def_rsf_prog(fid)
          call def_rsf_longacc(fid,.false.)
        else
c end-of-month acc file, no prognostic arrays
          call def_rsf_label(fid)
          call def_acc_all(fid,.true.)     ! real*4 disk storage
          call def_rsf_longacc(fid,.true.) ! real*4 disk storage
          call def_acc_meta(fid)
        endif
        call par_close(grid,fid)!if(am_i_root())status=nf_enddef(fid)
      elseif(iaction.eq.iowrite .and. it.eq.itimei) then
c rsf_file_name(1:2).nc at the beginning of a run
        do k=1,2
          fid=par_open(grid,trim(rsf_file_name(k))//'.nc','create')
          call def_rsf_label(fid)
          call def_rsf_prog(fid)
          call def_acc_all(fid,.false.)
          call def_rsf_longacc(fid,.false.)
          call par_close(grid,fid)
        enddo
      endif

c
c open file using the appropriate mode
c
      if(iaction.le.iowrite) then
        fid = par_open(grid,trim(fname)//'.nc','write')
      elseif(iaction.ge.ioread) then
        tmpname = trim(fname)//'.nc'
        if(fname == 'AIC') tmpname=fname ! symbolic link w/o .nc suffix
        fid = par_open(grid,trim(tmpname),'read')
      endif

      call new_io_label  (fid,iaction,it)

c
c prognostic arrays
c
      if(do_io_prog) then
        call new_io_model  (fid,iorw)
        call new_io_ocean  (fid,iorw)
        call new_io_lakes  (fid,iorw)
        call new_io_seaice (fid,iorw)
        call new_io_earth  (fid,iorw)
        call new_io_soils  (fid,iorw)
        call new_io_vegetation  (fid,iorw)
#ifdef USE_ENT
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
        call new_io_veg_related(fid,iorw)
        !call io_ent    (kunit,iaction,ioerr) ! io_vegetation handles ent
#endif
        call new_io_snow   (fid,iorw)
        call new_io_landice(fid,iorw)
        call new_io_bldat  (fid,iorw)
        call new_io_pbl    (fid,iorw)
        call new_io_clouds (fid,iorw)
        call new_io_somtq  (fid,iorw)
        call new_io_rad    (fid,iorw)
        call new_io_icedyn (fid,iorw)
#ifdef TRACERS_ON
        call new_io_tracer (fid,iorw)
#endif
      end if

c
c acc arrays
c
      if(do_io_acc) then
        if(iaction.eq.iowrite_single) call write_acc_meta(fid)
        call new_io_acc(fid,iaction)
        call new_io_ocdiag(fid,iorw)
        call new_io_icdiag(fid,iorw)
#ifdef TRACERS_ON
        call new_io_trdiag (fid,iaction)
#endif
      endif

c
c long-period acc arrays
c
      if(do_io_longacc) call new_io_longacc(fid,iorw)

c
c close input/output file
c
      call par_close(grid,fid)


      RETURN
      END SUBROUTINE io_rsf

      subroutine read_ground_ic
!@sum   read_ground_ic read initial conditions file for
!@+     sea ice, land surface, land ice.  Transplanted from INPUT.
!@auth  M. Kelley
!@ver   1.0
      use model_com, only : ioread
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close
      implicit none
      integer :: fid
      fid = par_open(grid,'GIC','read')
      call new_io_seaice (fid,ioread)
      call new_io_earth  (fid,ioread)
      call new_io_soils  (fid,ioread)
      call new_io_landice(fid,ioread)
      call par_close(grid,fid)
      return
      end subroutine read_ground_ic

      subroutine find_later_rsf(kdisk)
!@sum set kdisk such that Itime in rsf_file_name(kdisk) is
!@+   the larger of the Itimes in rsf_file_name(1:2).
!@auth  M. Kelley
!@ver   1.0
      use model_com, only : rsf_file_name
      use domain_decomp_atm, only : grid
      use pario, only : read_data,par_open,par_close
      implicit none
      integer, intent(out) :: kdisk
      integer :: itimes(2),fid,k
      itimes(:) = -1
      do k=1,2
        fid = par_open(grid,trim(rsf_file_name(k))//'.nc','read')
        call read_data(grid,fid,'itime',itimes(k),bcast_all=.true.)
        call par_close(grid,fid)
      enddo
      if(maxval(itimes).eq.-1) call stop_model(
     &     'FIND_LATER_RSF: ERRORS ON BOTH RESTART FILES',255)
      kdisk = 1
      if(itimes(2).gt.itimes(1)) kdisk=2
      return
      end subroutine find_later_rsf

      subroutine def_rsf_prog(fid)
!@sum  def_rsf_prog defines prognostic array structure in restart files
!@auth M. Kelley
!@ver  beta
      implicit none
      integer :: fid
      call def_rsf_model  (fid)
      call def_rsf_ocean  (fid)
      call def_rsf_lakes  (fid)
      call def_rsf_seaice (fid)
      call def_rsf_icedyn (fid)
      call def_rsf_earth  (fid)
      call def_rsf_soils  (fid)
      call def_rsf_vegetation(fid)
#ifdef USE_ENT
      call def_rsf_veg_related(fid)
#endif
      call def_rsf_snow   (fid)
      call def_rsf_landice(fid)
      call def_rsf_bldat  (fid)
      call def_rsf_pbl    (fid)
      call def_rsf_clouds (fid)
      call def_rsf_somtq  (fid)
      call def_rsf_rad    (fid)
#ifdef TRACERS_ON
      call def_rsf_tracer (fid)
#endif
      return
      end subroutine def_rsf_prog

      subroutine def_acc_meta(fid)
!@sum  def_acc_meta defines metadata in acc files
!@auth M. Kelley
!@ver  beta
      USE MODEL_COM, only :
     &     jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,itime,itime0,nday
      use domain_decomp_atm, only : grid
      use pario, only : write_attr
      implicit none
      integer :: fid         !@var fid file id
      character(len=100) :: fromto
      real*8 :: days

      days=(itime-itime0)/float(nday)
      write(fromto,902) jyear0,amon0,jdate0,jhour0,
     &     jyear,amon,jdate,jhour,itime,days
      call write_attr(grid,fid,'global','fromto',fromto)

c idacc(5) is not additive
      call write_attr(grid,fid,'idacc','reduction','sum')

      call def_meta_atmacc(fid)
      call def_meta_rvracc(fid)
      call def_meta_ocdiag(fid)
      call def_meta_icdiag(fid)
#ifdef TRACERS_ON
      call def_meta_trdiag(fid)
#endif
      return
  902 FORMAT ('From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      end subroutine def_acc_meta

      subroutine write_acc_meta(fid)
!@sum  def_acc_meta writes metadata to acc files
!@auth M. Kelley
!@ver  beta
      implicit none
      integer :: fid         !@var fid file id
      call write_meta_atmacc(fid)
      call write_meta_rvracc(fid)
      call write_meta_ocdiag(fid)
      call write_meta_icdiag(fid)
#ifdef TRACERS_ON
      call write_meta_trdiag(fid)
#endif
      return
      end subroutine write_acc_meta

      subroutine def_acc_all(fid,r4_on_disk)
!@sum  def_acc_all defines acc array structure in restart files
!@auth M. Kelley
!@ver  beta
      implicit none
      integer :: fid         !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4
      call def_rsf_acc(fid,r4_on_disk)
      call def_rsf_ocdiag(fid,r4_on_disk)
      call def_rsf_icdiag(fid,r4_on_disk)
#ifdef TRACERS_ON
      call def_rsf_trdiag(fid,r4_on_disk)
#endif
      return
      end subroutine def_acc_all

      subroutine def_rsf_label(fid)
!@sum  def_rsf_label defines control info in restart/acc files
!@auth M. Kelley
!@ver  beta
      use model_com, only : xlabel,itime,itimee,itime0,itimei,
     &     iyear1,nday,iowrite
      use timings, only : ntimemax,ntimeacc,timestr,timing
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      implicit none
      integer fid   !@var fid file id
      integer :: n
      integer :: intdum
      call write_attr(grid,fid,'global','xlabel',xlabel)
      call defvar(grid,fid,itime,'itime')
      call write_caldate(fid)
      call defvar(grid,fid,itimee,'itimee')
      call defvar(grid,fid,itime0,'itime0')
      call defvar(grid,fid,itimei,'itimei')
      call defvar(grid,fid,iyear1,'iyear1')
      call defvar(grid,fid,nday,'nday')
      call defvar(grid,fid,timing(1:ntimemax),'cputime(ntimemax)')
      call defvar(grid,fid,ntimeacc,'ntimeacc')
      call io_cputime(fid,iowrite)
c rparam, iparam, cparam are dummy vars which hold the
c parameter database in their attributes.  On disk,
c their actual value will be set to the number of
c parameters that are defined.
      call defvar(grid,fid,intdum,'rparam')
      call defvar(grid,fid,intdum,'iparam')
      call defvar(grid,fid,intdum,'cparam')
      call new_io_param(fid,iowrite,.false.,.false.)
      return
      end subroutine def_rsf_label

      subroutine new_io_label(fid,iaction,ihrX)
!@sum  new_io_label read/write control info from/to restart,acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : xlabel,itime,itimee,itime0,itimei,
     &     iyear1,nday,iowrite,ioread,irsfic,irsficnt,irsficno
      use domain_decomp_atm, only: grid
      use pario, only : write_data,read_data,write_attr,read_attr
      use timings, only : ntimemax,ntimeacc,timestr,timing
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: ihrX !@var ihrX a dummy itime read from IC files
      integer :: idum,nday_dummy
      logical :: is_ic

      select case (iaction)
      case (:iowrite) ! output to restart or acc file
        call write_data(grid, fid, 'itime', itime)
        call write_caldate(fid)
        call write_data(grid, fid, 'itimee', itimee)
        call write_data(grid, fid, 'itime0', itime0)
        call write_data(grid, fid, 'itimei', itimei)
        call write_data(grid, fid, 'iyear1', iyear1)
        call write_data(grid, fid, 'nday', nday)
        call write_data(grid, fid, 'ntimeacc', ntimeacc)
        call write_data(grid, fid, 'cputime', timing(1:ntimemax))
        call io_cputime(fid,iowrite)
        call new_io_param(fid,iowrite,.false.,.true.)
      case (ioread:) ! input from restart or acc file
        call read_attr(grid,fid,'global','xlabel',idum,xlabel)
        call read_data(grid,fid,'nday',nday_dummy,bcast_all=.true.)
        is_ic = .false.
        if(iaction.eq.irsfic)   is_ic = .true.
        if(iaction.eq.irsficnt) is_ic = .true.
        if(iaction.eq.irsficno) is_ic = .true.
        if(.not.is_ic) then
          nday = nday_dummy
          call read_data(grid,fid,'itime', itime, bcast_all=.true.)
          call read_data(grid,fid,'itimee',itimee,bcast_all=.true.)
          call read_data(grid,fid,'itime0',itime0,bcast_all=.true.)
          call read_data(grid,fid,'itimei',itimei,bcast_all=.true.)
          if(iyear1.lt.0) ! this check not present for ioread_single!!!
     &    call read_data(grid,fid,'iyear1',iyear1,bcast_all=.true.)
          call read_data(grid,fid,'ntimeacc',ntimeacc,
     &         bcast_all=.true.)
          call read_data(grid,fid, 'cputime', timing(1:ntimemax),
     &         bcast_all=.true.)
          call io_cputime(fid,ioread)
          call new_io_param(fid,ioread,.false.,.false.)
        else
          call read_data(grid,fid,'itime', IhrX, bcast_all=.true.)
          IhrX=IhrX*24/nday_dummy
        endif
      end select
      end subroutine new_io_label

      subroutine write_caldate(fid)
c write a text version of the date to a restart/acc file
      use model_com
      use pario, only : write_attr
      use domain_decomp_atm, only: grid
      implicit none
      integer :: fid
      character(len=18) :: caldate
      character(len=2) :: cmo,cday
      character(len=4) :: cyr,chr
      write(cmo,'(i2.2)') jmon
      write(cday,'(i2.2)') jdate
      write(cyr,'(i4.4)') jyear
      write(chr,'(f4.1)') real(jhour)
      caldate=cmo//'/'//cday//'/'//cyr//' hr '//chr
      call write_attr(grid,fid,'itime','caldate',caldate)
      return
      end subroutine write_caldate

      subroutine io_cputime(fid,iaction)
c manage the reading/writing of timing information. could be done better
      use model_com, only : iowrite,ioread,nday,itime,itime0
      use timings, only : ntimemax,ntimeacc,timestr,timing
      use domain_decomp_atm, only: grid
      use pario, only : write_attr,read_attr
      implicit none
      integer :: fid,iaction
      character*2 :: c2
      character(len=5) :: cfrac
      character(len=6) :: fracnames(ntimemax)
      character(len=17) :: timestr_pct
      integer :: n,idum
      real*8 :: tsum,min_per_day
      do n=1,ntimemax
        write(c2,'(i2.2)') n
        fracnames(n) = 'frac'//c2
      enddo
      select case (iaction)
      case (:iowrite)
        tsum = sum(timing(1:ntimeacc))+1d-20
        if(itime.gt.itime0) then
          min_per_day = (tsum/60.)*nday/dble(itime-itime0)
        else
          min_per_day = 0.
        endif
        min_per_day = int(min_per_day*100d0)*.01d0
        call write_attr(grid,fid,'cputime','minutes_per_day',
     &       min_per_day)
        do n=1,ntimeacc
          write(cfrac,'(f5.1)') 100.*timing(n)/tsum
          timestr_pct = timestr(n)//cfrac
          call write_attr(grid,fid,'cputime',fracnames(n),
     &         timestr_pct)
        enddo
      case (ioread:)
        do n=1,ntimeacc
          call read_attr(grid,fid,'cputime',fracnames(n),idum,
     &         timestr_pct)
          timestr(n) = timestr_pct(1:12)
        enddo
      end select
      return
      end subroutine io_cputime

      subroutine new_io_param(fid,iaction,ovrwrt,new_possible)
!@sum  new_io_param read/write parameter database from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com
      use domain_decomp_atm, only : grid
      use pario, only : write_attr, read_attr, read_data, write_data,
     &     par_enddef
      use param
      implicit none
      integer fid   !@var fid file id
      integer iaction !@var iaction flag for reading or writing to file
      logical :: ovrwrt,new_possible
      integer :: l,m,n,plen,strlen,niparam,nrparam,ncparam
      character*1 :: ptype
      integer :: pvali(1000)
      real*8 :: pvalr(1000)
      character(len=128) :: pname,pvalc(1000)
      integer, parameter :: lcstr=1000
      character(len=lcstr) :: cstr
      select case (iaction)
      case (iowrite) ! output
        niparam = 0; nrparam = 0; ncparam = 0;
        n = 0
        do
          n = n + 1
          call query_param( n, pname, plen, ptype )
          if(pname(1:5).eq.'EMPTY') exit
          select case(ptype)
          case ('i')
            call get_param(pname,pvali,plen)
            call write_attr(grid,fid,'iparam',trim(pname),
     &           pvali(1:plen))
            niparam = niparam + 1
          case ('r')
            call get_param(pname,pvalr,plen)
            call write_attr(grid,fid,'rparam',trim(pname),
     &           pvalr(1:plen))
            nrparam = nrparam + 1
          case ('c')
            call get_param(pname,pvalc,plen)
            cstr = trim(pvalc(1))
c Arrays of strings are not supported here.  Separate the
c strings with the character "|".
            do m=2,plen
              cstr=trim(cstr)//'|'//trim(pvalc(m))
            enddo
            if(len_trim(cstr).gt.0) then
              call write_attr(grid,fid,'cparam',trim(pname),cstr)
              ncparam = ncparam + 1
            endif
          end select
        enddo
        if(new_possible) call par_enddef(grid,fid)
        call write_data(grid, fid, 'iparam', niparam)
        call write_data(grid, fid, 'rparam', nrparam)
        call write_data(grid, fid, 'cparam', ncparam)
      case (ioread) ! input
        call read_data(grid,fid,'iparam',niparam,bcast_all=.true.)
        call read_data(grid,fid,'rparam',nrparam,bcast_all=.true.)
        call read_data(grid,fid,'cparam',ncparam,bcast_all=.true.)
        do n=1,niparam
          call read_attr(grid,fid,'iparam',pname,plen,pvali,
     &         attnum=n)
          if ( (.not. is_set_param(pname)) .or. ovrwrt )
     &         call set_param( trim(pname), pvali, plen, 'o' )
        enddo
        do n=1,nrparam
          call read_attr(grid,fid,'rparam',pname,plen,pvalr,
     &         attnum=n)
          if ( (.not. is_set_param(pname)) .or. ovrwrt )
     &         call set_param( trim(pname), pvalr, plen, 'o' )
        enddo
        do n=1,ncparam
          cstr(lcstr:lcstr)='x' ! workaround for length checking
          call read_attr(grid,fid,'cparam',pname,strlen,cstr,
     &         attnum=n)
c Arrays of strings are not supported here.  Strings containing
c the character "|" are separated into arrays of strings.
          l = 0
          plen = 1
          pvalc = ''
          do m=1,strlen
            l = l + 1
            if(cstr(m:m).eq.'|') then
              l = 0
              plen = plen + 1
            else
              pvalc(plen)(l:l) = cstr(m:m)
            endif
          enddo
          if ( (.not. is_set_param(pname)) .or. ovrwrt )
     &         call set_param( trim(pname), pvalc, plen, 'o' )
        enddo
      end select
      return
      end subroutine new_io_param

      subroutine set_ioptrs_acc_default
c point i/o pointers for accumlated quantities to the
c instances of the arrays used during normal operation.
      implicit none
      call set_ioptrs_atmacc_default
      call set_ioptrs_ocnacc_default
      call set_ioptrs_iceacc_default
      return
      end subroutine set_ioptrs_acc_default

      subroutine set_ioptrs_acc_extended
c point i/o pointers for accumlated quantities to the
c instances of the arrays holding derived quantities as well
      implicit none
      call set_ioptrs_atmacc_extended
      call set_ioptrs_ocnacc_extended
c      call set_ioptrs_iceacc_extended
      return
      end subroutine set_ioptrs_acc_extended

      subroutine calc_derived_acc
      implicit none
      call calc_derived_acc_atm
      call diag_ocean_prep
#ifdef TRACERS_ON
      call diag_trac_prep
#endif
      return
      end subroutine calc_derived_acc
