#include "rundeck_opts.h"

c
c The routines in this file are demo versions of how to write
c checkpoint files using the pario module.  Eventually they should
c be distributed among the model modules owning the arrays being
c input/output.  They are consolidated here for testing purposes.
c The driver subroutine io_rsf is intended to be plug-compatible
c with the operational version of io_rsf.
c
c M. Kelley 12/2008
c


      SUBROUTINE io_rsf(kunit,it,iaction,ioerr)
      USE DOMAIN_DECOMP, only : am_i_root
      USE pario_fbsa, only : REWIND_PARALLEL
      USE MODEL_COM, only : ioread_single,iowrite_single,Kradia
     *                     ,ioread,ioread_nodiag,iowrite_mon
     &     ,itimei
      IMPLICIT NONE
!@var iaction flag for reading or writing rsf file
!@var kunit Fortran unit number of file i/o
      INTEGER, INTENT(IN) :: iaction,kunit
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var IT1 hour for correct reading check
!@var ITM maximum hour for post-processing
      INTEGER IT1,itm,iact
      logical skip_diag
      include 'netcdf.inc'
      integer :: fid,k,status
      character*1 ck
      character(len=80) :: fname

      iact=iaction   ; skip_diag=.false.
      if(iaction==ioread_nodiag) then
         iact=ioread ; skip_diag=.true.
      end if

      ioerr=-1
      if (kunit.gt.0) CALL REWIND_PARALLEL( kunit )

      if(it.eq.itimei) then
c define restart files fort.[12].nc at the beginning of a run
        do k=1,2
          write(ck,'(i1)') k
          if(am_i_root()) then
            status = nf_create('fort.'//ck//'.nc',nf_clobber,fid)
          endif
          call def_rsf_all(fid,.true.)
          if(am_i_root()) then
            status = nf_enddef(fid)
            status = nf_close(fid)
          endif
        enddo
      endif
c determine which netcdf restart file we are writing and open it
      if(iaction.ne.iowrite_mon) then ! fort.[12]
        if(am_i_root()) then
          inquire(kunit,name=fname)
          status = nf_open(trim(fname)//'.nc',nf_write,fid)
        endif
      else ! end-of-month restart file.  do not define acc arrays.
        if(am_i_root()) then
          status = nf_create(trim(fname)//'.nc',nf_clobber,fid)
        endif
        call def_rsf_all(fid,.false.)
        if(am_i_root()) status = nf_enddef(fid)
      endif

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
C**** Particular values may produce variations in indiv. i/o routines

C**** Calls to individual i/o routines
      call io_label  (kunit,it,itm,iact,ioerr)
      it1=it
      if(iaction.ne.ioread_single.and.iaction.ne.iowrite_single) then
        call par_io_model  (fid,iact) 
        call io_strat  (kunit,iact,ioerr)
        call par_io_ocean  (fid,iact)
        call par_io_lakes  (fid,iact)
        call par_io_seaice (fid,iact)
        call par_io_earth  (fid,iact)
        call par_io_soils  (fid,iact)
        call par_io_vegetation  (fid,iact)
#ifdef USE_ENT
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
        call io_veg_related  (kunit,iact,ioerr)
        !call io_ent    (kunit,iaction,ioerr)
#endif
        call par_io_snow   (fid,iact)
        call par_io_landice(fid,iact)
        call par_io_bldat  (fid,iact)
        call par_io_pbl    (fid,iact)
        call par_io_clouds (fid,iact)
        call par_io_somtq  (fid,iact)
        call par_io_rad    (fid,iact)
        call io_icedyn (kunit,iact,ioerr)
#ifdef TRACERS_ON
        call io_tracer (kunit,iact,ioerr)
#endif
      end if
      if(skip_diag) go to 10
      if(iaction.eq.ioread_single.or.iaction.eq.iowrite_single) then
        call io_diags  (kunit,it,iaction,ioerr)
      else
        call par_io_acc(fid,iaction)
      endif
      call io_ocdiag (kunit,it,iaction,ioerr)
      call io_icdiag (kunit,it,iaction,ioerr)
#ifdef TRACERS_ON
      call io_trdiag (kunit,it,iaction,ioerr)
#endif

      if (it1.ne.it) THEN
        WRITE(6,*) "TIMES DO NOT MATCH READING IN RSF FILE",it,it1
        ioerr=1
      END IF
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN RESTART FILE: KUNIT="
     *     ,kunit

C**** return maximum time
   10 it=itm

c close netcdf restart file
      if(am_i_root()) status = nf_close(fid)

      RETURN
      END SUBROUTINE io_rsf

      subroutine read_ground_ic
!@sum   read_ground_ic read initial conditions file for
!@+     sea ice, land surface, land ice.  Transplanted from INPUT.
!@auth  M. Kelley
!@ver   1.0
!@calls io_seaice,io_earth,io_soils,io_landice
      use model_com, only : ioreadnt,ioread
      use filemanager, only : openunit,closeunit,nameunit
      use domain_decomp, only : am_i_root
      implicit none
      include 'netcdf.inc'
      integer :: iu_GIC,ioerr,status,fid
      character*16 :: name
#ifdef CUBE_GRID
c      if(am_i_root()) then
c     here we use openunit and nameunit only to get the file name
         status = nf_open("GIC",nf_nowrite,fid)
      IF (status .ne. NF_NOERR) write(*,*) "nf_open error"
c      endif
      write(*,*) "Reading ground IC"
      call par_io_seaice (fid,ioread)
      call par_io_earth  (fid,ioread)
      call par_io_soils  (fid,ioread)
      call par_io_landice(fid,ioread)
      status = nf_close(fid)
#else
      call openunit("GIC",iu_GIC,.true.,.true.)
      write(*,*) "Reading ground IC"
      ioerr=-1
      read(iu_GIC) ! ignore first line (ocean ic done in init_OCEAN)
      call io_seaice (iu_GIC,ioreadnt,ioerr)
      call io_earth  (iu_GIC,ioreadnt,ioerr)
      call io_soils  (iu_GIC,ioreadnt,ioerr)
      call io_landice(iu_GIC,ioreadnt,ioerr)
      if (ioerr.eq.1) then
        if (AM_I_ROOT())
     *       WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
        call stop_model("INPUT: GIC READ IN ERROR",255)
      end if
      call closeunit (iu_GIC)
#endif
      return
      end subroutine read_ground_ic

      subroutine def_rsf_all(fid,with_acc)
      integer :: fid
      logical :: with_acc
      call def_rsf_model  (fid) 
      call def_rsf_ocean  (fid)
      call def_rsf_lakes  (fid)
      call def_rsf_seaice (fid)
      call def_rsf_earth  (fid)
      call def_rsf_soils  (fid)
      call def_rsf_vegetation(fid)
      call def_rsf_snow   (fid)
      call def_rsf_landice(fid)
      call def_rsf_bldat  (fid)
      call def_rsf_pbl    (fid)
      call def_rsf_clouds (fid)
      call def_rsf_somtq  (fid)
      call def_rsf_rad    (fid)
      if(with_acc) then
        call def_rsf_acc(fid)
      endif
      return
      end subroutine def_rsf_all

      subroutine def_rsf_model(fid)
      use model_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,u,'u(dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,v,'v(dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,t,'t(dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,q,'q(dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,wm,'wm(dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,p,'p(dist_im,dist_jm)')
#ifdef BLK_2MOM
      call defvar(grid%dd2d,fid,wmice,'wmice(dist_im,dist_jm,lm)')
#endif
      return
      end subroutine def_rsf_model
      subroutine par_io_model(fid,iaction)
      use model_com
      use domain_decomp, only: grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite) ! output to end-of-month restart file
        call write_dist_data(grid%dd2d, fid, 'u', u)
        call write_dist_data(grid%dd2d, fid, 'v', v)
        call write_dist_data(grid%dd2d, fid, 't', t)
        call write_dist_data(grid%dd2d, fid, 'p', p)
        call write_dist_data(grid%dd2d, fid, 'q', q)
        call write_dist_data(grid%dd2d, fid, 'wm', wm)
#ifdef BLK_2MOM
        call write_dist_data(grid%dd2d, fid, 'wmice', wmice)
#endif
      case (ioread)          ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'u', u)
        call read_dist_data(grid%dd2d, fid, 'v', v)
        call read_dist_data(grid%dd2d, fid, 't', t)
        call read_dist_data(grid%dd2d, fid, 'p', p)
        call read_dist_data(grid%dd2d, fid, 'q', q)
        call read_dist_data(grid%dd2d, fid, 'wm', wm)
#ifdef BLK_2MOM
        call read_dist_data(grid%dd2d, fid, 'wmice', wmice)
#endif
      end select
      return
      end subroutine par_io_model

      subroutine def_rsf_somtq(fid)
      use somtq_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,tmom,'tmom(nmom,dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,qmom,'qmom(nmom,dist_im,dist_jm,lm)')
      return
      end subroutine def_rsf_somtq
      subroutine par_io_somtq(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      use somtq_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'tmom', tmom, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'qmom', qmom, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'tmom', tmom, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'qmom', qmom, jdim=3)
      end select
      return
      end subroutine par_io_somtq

      subroutine def_rsf_clouds(fid)
      use clouds_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,ttold,'ttold(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,qtold,'qtold(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,svlhx,'svlhx(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,rhsav,'rhsav(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,cldsav,'cldsav(lm,dist_im,dist_jm)')
#ifdef CLD_AER_CDNC
      call defvar(grid%dd2d,fid,oldnl,'oldnl(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,oldni,'oldni(lm,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_clouds
      subroutine par_io_clouds(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      use clouds_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'ttold', ttold, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'qtold', qtold, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'svlhx', svlhx, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'rhsav', rhsav, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'cldsav', cldsav, jdim=3)
#ifdef CLD_AER_CDNC
        call write_dist_data(grid%dd2d, fid, 'oldnl', oldnl, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'oldni', oldni, jdim=3)
#endif
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'ttold', ttold, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'qtold', qtold, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'svlhx', svlhx, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'rhsav', rhsav, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'cldsav', cldsav, jdim=3)
#ifdef CLD_AER_CDNC
        call read_dist_data(grid%dd2d, fid, 'oldnl', oldnl, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'oldni', oldni, jdim=3)
#endif
      end select
      return
      end subroutine par_io_clouds

      subroutine def_rsf_rad(fid)
      use rad_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,s0,'s0')
      call defvar(grid%dd2d,fid,rqt,'rqt(lm_req,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,kliq,'kliq(lm,d4,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,srhr,'srhr(zero_to_lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,trhr,'trhr(zero_to_lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,trsurf,'trsurf(nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,fsf,'fsf(nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,fsrdir,'fsrdir(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,srvissurf,'srvissurf(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,srdn,'srdn(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,cfrac,'cfrac(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,salb,'salb(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,fsrdif,'fsrdif(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,dirnir,'dirnir(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,difnir,'difnir(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,rcld,'rcld(lm,dist_im,dist_jm)')
#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_Shindell
      call defvar(grid%dd2d,fid,o3_tracer_save,
     &     'o3_tracer_save(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,rad_to_chem,
     &     'rad_to_chem(d5,lm,dist_im,dist_jm)')
#endif
#ifdef TRACERS_DUST
      call defvar(grid%dd2d,fid,srnflb_save,
     &     'srnflb_save(dist_im,dist_jm,lm)')
      call defvar(grid%dd2d,fid,trnflb_save,
     &     'trnflb_save(dist_im,dist_jm,lm)')
#endif
      call defvar(grid%dd2d,fid,ttausv_save,
     &     'ttausv_save(dist_im,dist_jm,ntm,lm)')
      call defvar(grid%dd2d,fid,ttausv_cs_save,
     &     'ttausv_cs_save(dist_im,dist_jm,ntm,lm)')
#endif
#ifdef HTAP_LIKE_DIAGS
      call defvar(grid%dd2d,fid,ttausv_sum,
     &     'ttausv_sum(dist_im,dist_jm,ntm)')
      call defvar(grid%dd2d,fid,ttausv_count,'ttausv_count')
#endif
      return
      end subroutine def_rsf_rad

      subroutine par_io_rad(fid,iaction)
      use model_com, only : ioread,iowrite
#ifdef TRACERS_ON
      USE tracer_com , only : ntm
#endif
      use rad_com
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_data(grid%dd2d, fid,'s0', s0)
        call write_dist_data(grid%dd2d, fid,'rqt',  rqt, jdim=3)
        call write_dist_data(grid%dd2d, fid,'kliq', kliq, jdim=4)
        call write_dist_data(grid%dd2d, fid,'srhr', srhr, jdim=3)
        call write_dist_data(grid%dd2d, fid,'trhr', trhr, jdim=3)
        call write_dist_data(grid%dd2d, fid,'trsurf', trsurf, jdim=3)
        call write_dist_data(grid%dd2d, fid,'fsf',  fsf, jdim=3)
        call write_dist_data(grid%dd2d, fid,'salb', salb)
        call write_dist_data(grid%dd2d, fid,'fsrdir', fsrdir)
        call write_dist_data(grid%dd2d, fid,'srvissurf', srvissurf)
        call write_dist_data(grid%dd2d, fid,'fsrdif', fsrdif)
        call write_dist_data(grid%dd2d, fid,'dirnir', dirnir)
        call write_dist_data(grid%dd2d, fid,'difnir', difnir)
        call write_dist_data(grid%dd2d, fid,'srdn',   srdn)
        call write_dist_data(grid%dd2d, fid,'cfrac',  cfrac)
        call write_dist_data(grid%dd2d, fid,'rcld', rcld, jdim=3)
#ifdef TRACERS_SPECIAL_Shindell
        call write_dist_data(grid%dd2d,fid,
     &       'o3_tracer_save', o3_tracer_save, jdim=3)
        call write_dist_data(grid%dd2d,fid,
     &       'rad_to_chem', rad_to_chem, jdim=4)
#endif
#ifdef TRACERS_DUST
        call write_dist_data(grid%dd2d,fid,'srnflb_save',srnflb_save)
        call write_dist_data(grid%dd2d,fid,'trnflb_save',trnflb_save)
#endif
#ifdef TRACERS_ON
        call write_dist_data(grid%dd2d,fid,'ttausv_save',ttausv_save)
        call write_dist_data(grid%dd2d,fid,'ttausv_cs_save',
     &       ttausv_cs_save)
#endif
#ifdef HTAP_LIKE_DIAGS
        call write_data(grid%dd2d, fid,'ttausv_count',ttausv_count)
        call write_dist_data(grid%dd2d,fid,'ttausv_sum',ttausv_sum)
#endif
      case (ioread)
        call read_data(grid%dd2d, fid,'s0', s0, bcast_all=.true.)
        call read_dist_data(grid%dd2d, fid,'rqt',  rqt, jdim=3)
        call read_dist_data(grid%dd2d, fid,'kliq', kliq, jdim=4)
        call read_dist_data(grid%dd2d, fid,'srhr', srhr, jdim=3)
        call read_dist_data(grid%dd2d, fid,'trhr', trhr, jdim=3)
        call read_dist_data(grid%dd2d, fid,'trsurf', trsurf, jdim=3)
        call read_dist_data(grid%dd2d, fid,'fsf',  fsf, jdim=3)
        call read_dist_data(grid%dd2d, fid,'salb', salb)
        call read_dist_data(grid%dd2d, fid,'fsrdir', fsrdir)
        call read_dist_data(grid%dd2d, fid,'srvissurf', srvissurf)
        call read_dist_data(grid%dd2d, fid,'fsrdif', fsrdif)
        call read_dist_data(grid%dd2d, fid,'dirnir', dirnir)
        call read_dist_data(grid%dd2d, fid,'difnir', difnir)
        call read_dist_data(grid%dd2d, fid,'srdn',   srdn)
        call read_dist_data(grid%dd2d, fid,'cfrac',  cfrac)
        call read_dist_data(grid%dd2d, fid,'rcld', rcld, jdim=3)
#ifdef TRACERS_SPECIAL_Shindell
        call read_dist_data(grid%dd2d,fid,
     &       'o3_tracer_save', o3_tracer_save, jdim=3)
        call read_dist_data(grid%dd2d,fid,
     &       'rad_to_chem', rad_to_chem, jdim=4)
#endif
#ifdef TRACERS_DUST
        call read_dist_data(grid%dd2d,fid,'srnflb_save',srnflb_save)
        call read_dist_data(grid%dd2d,fid,'trnflb_save',trnflb_save)
#endif
#ifdef TRACERS_ON
        call read_dist_data(grid%dd2d,fid,'ttausv_save',ttausv_save)
        call read_dist_data(grid%dd2d,fid,'ttausv_cs_save',
     &       ttausv_cs_save)
#endif
#ifdef HTAP_LIKE_DIAGS
        call read_data(grid%dd2d, fid,'ttausv_count',ttausv_count,
     &       bcast_all=.true.)
        call read_dist_data(grid%dd2d,fid,'ttausv_sum',ttausv_sum)
#endif
      end select
      return
      end subroutine par_io_rad

      subroutine def_rsf_pbl(fid)
      use pblcom
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,uabl,
     &     'uabl(npbl,nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,vabl,
     &     'vabl(npbl,nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tabl,
     &     'tabl(npbl,nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,qabl,
     &     'qabl(npbl,nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,eabl,
     &     'eabl(npbl,nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,cmgs,'cmgs(nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,chgs,'chgs(nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,cqgs,'cqgs(nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,ipbl,'ipbl(nstype,dist_im,dist_jm)')
#ifdef TRACERS_ON
      call defvar(grid%dd2d,fid,trabl,
     &     'trabl(npbl,ntm,nstype,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_pbl

      subroutine par_io_pbl(fid,iaction)
      use model_com, only : ioread,iowrite
      use pblcom
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'uabl', uabl, jdim=4)
        call write_dist_data(grid%dd2d, fid, 'vabl', vabl, jdim=4)
        call write_dist_data(grid%dd2d, fid, 'tabl', tabl, jdim=4)
        call write_dist_data(grid%dd2d, fid, 'qabl', qabl, jdim=4)
        call write_dist_data(grid%dd2d, fid, 'eabl', eabl, jdim=4)
        call write_dist_data(grid%dd2d, fid, 'cmgs', cmgs, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'chgs', chgs, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'cqgs', cqgs, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'ipbl', ipbl, jdim=3)
#ifdef TRACERS_ON
        call write_dist_data(grid%dd2d, fid, 'trabl', trabl, jdim=5)
#endif
      case (ioread)            ! input from restart file or restart
        call read_dist_data(grid%dd2d, fid, 'uabl', uabl, jdim=4)
        call read_dist_data(grid%dd2d, fid, 'vabl', vabl, jdim=4)
        call read_dist_data(grid%dd2d, fid, 'tabl', tabl, jdim=4)
        call read_dist_data(grid%dd2d, fid, 'qabl', qabl, jdim=4)
        call read_dist_data(grid%dd2d, fid, 'eabl', eabl, jdim=4)
        call read_dist_data(grid%dd2d, fid, 'cmgs', cmgs, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'chgs', chgs, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'cqgs', cqgs, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'ipbl', ipbl, jdim=3)
#ifdef TRACERS_ON
        call read_dist_data(grid%dd2d, fid, 'trabl', trabl, jdim=5)
#endif
      end select
      return
      end subroutine par_io_pbl

      subroutine def_rsf_bldat(fid)
      use pblcom
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,wsavg,'wsavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tsavg,'tsavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,qsavg,'qsavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,dclev,'dclev(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,usavg,'usavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,vsavg,'vsavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tauavg,'tauavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tgvavg,'tgvavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,qgavg,'qgavg(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,ustar_pbl,
     &     'ustar_pbl(nstype,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,egcm,'egcm(lm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,w2gcm,'w2gcm(lm,dist_im,dist_jm)')
      return
      end subroutine def_rsf_bldat

      subroutine par_io_bldat(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      use pblcom
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d,fid,'wsavg',wsavg)
        call write_dist_data(grid%dd2d,fid,'tsavg',tsavg)
        call write_dist_data(grid%dd2d,fid,'qsavg',qsavg)
        call write_dist_data(grid%dd2d,fid,'dclev',dclev)
        call write_dist_data(grid%dd2d,fid,'usavg',usavg)
        call write_dist_data(grid%dd2d,fid,'vsavg',vsavg)
        call write_dist_data(grid%dd2d,fid,'tauavg',tauavg)
        call write_dist_data(grid%dd2d,fid,'tgvavg',tgvavg)
        call write_dist_data(grid%dd2d,fid,'qgavg',qgavg)
        call write_dist_data(grid%dd2d,fid,'ustar_pbl',ustar_pbl,jdim=3)
        call write_dist_data(grid%dd2d,fid,'egcm',egcm, jdim=3)
        call write_dist_data(grid%dd2d,fid,'w2gcm',w2gcm, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d,fid,'wsavg',wsavg)
        call read_dist_data(grid%dd2d,fid,'tsavg',tsavg)
        call read_dist_data(grid%dd2d,fid,'qsavg',qsavg)
        call read_dist_data(grid%dd2d,fid,'dclev',dclev)
        call read_dist_data(grid%dd2d,fid,'usavg',usavg)
        call read_dist_data(grid%dd2d,fid,'vsavg',vsavg)
        call read_dist_data(grid%dd2d,fid,'tauavg',tauavg)
        call read_dist_data(grid%dd2d,fid,'tgvavg',tgvavg)
        call read_dist_data(grid%dd2d,fid,'qgavg',qgavg)
        call read_dist_data(grid%dd2d,fid,'ustar_pbl',ustar_pbl,jdim=3)
        call read_dist_data(grid%dd2d,fid,'egcm',egcm, jdim=3)
        call read_dist_data(grid%dd2d,fid,'w2gcm',w2gcm, jdim=3)
      end select
      return
      end subroutine par_io_bldat

      subroutine def_rsf_earth(fid)
      use ghy_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,snowe,'snowe(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tearth,'tearth(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,wearth,'wearth(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,aiearth,'aiearth(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,evap_max_ij,
     &     'evap_max_ij(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,fr_sat_ij,
     &     'fr_sat_ij(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,qg_ij,'qg_ij(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tsns_ij,'tsns_ij(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,snoage,'snoage(d3,dist_im,dist_jm)')
      return
      end subroutine def_rsf_earth

      subroutine par_io_earth(fid,iaction)
      use model_com, only : ioread,iowrite
      use ghy_com
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d,fid,'snowe',snowe)
        call write_dist_data(grid%dd2d,fid,'tearth',tearth)
        call write_dist_data(grid%dd2d,fid,'wearth',wearth)
        call write_dist_data(grid%dd2d,fid,'aiearth',aiearth)
        call write_dist_data(grid%dd2d,fid,'evap_max_ij',evap_max_ij)
        call write_dist_data(grid%dd2d,fid,'fr_sat_ij',fr_sat_ij)
        call write_dist_data(grid%dd2d,fid,'qg_ij',qg_ij)
        call write_dist_data(grid%dd2d,fid,'snoage',snoage,jdim=3)
        call write_dist_data(grid%dd2d,fid,'tsns_ij',tsns_ij)
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d,fid,'snowe',snowe)
        call read_dist_data(grid%dd2d,fid,'tearth',tearth)
        call read_dist_data(grid%dd2d,fid,'wearth',wearth)
        call read_dist_data(grid%dd2d,fid,'aiearth',aiearth)
        call read_dist_data(grid%dd2d,fid,'evap_max_ij',evap_max_ij)
        call read_dist_data(grid%dd2d,fid,'fr_sat_ij',fr_sat_ij)
        call read_dist_data(grid%dd2d,fid,'qg_ij',qg_ij)
        call read_dist_data(grid%dd2d,fid,'snoage',snoage,jdim=3)
        tsns_ij(:,:) = tearth(:,:) ! default if not in input file
        call read_dist_data(grid%dd2d,fid,'tsns_ij',tsns_ij)
      end select
      return
      end subroutine par_io_earth

      subroutine def_rsf_soils(fid)
      use ghy_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,w_ij,
     &     'w_ij(zero_to_ngm,ls_nfrac,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,ht_ij,
     &     'ht_ij(zero_to_ngm,ls_nfrac,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,snowbv,
     &     'snowbv(ls_nfrac,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid%dd2d,fid,w_ij,
     &     'tr_w_ij(ntm,zero_to_ngm,ls_nfrac,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,trsnowbv0,
     &     'trsnowbv0(ntm,bv,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_soils

      subroutine par_io_soils(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only: grid
      use pario, only : write_dist_data,read_dist_data
      use ghy_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d,fid,'w_ij',w_ij,jdim=4)
        call write_dist_data(grid%dd2d,fid,'ht_ij',ht_ij,jdim=4)
        call write_dist_data(grid%dd2d,fid,'snowbv',snowbv,jdim=3)
#ifdef TRACERS_WATER
        call write_dist_data(grid%dd2d,fid,'tr_w_ij',tr_w_ij,jdim=5)
        call write_dist_data(grid%dd2d,fid,'trsnowbv0',trsnowbv0,jdim=4)
#endif
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d,fid,'w_ij',w_ij,jdim=4)
        call read_dist_data(grid%dd2d,fid,'ht_ij',ht_ij,jdim=4)
        call read_dist_data(grid%dd2d,fid,'snowbv',snowbv,jdim=3)
#ifdef TRACERS_WATER
        call read_dist_data(grid%dd2d,fid,'tr_w_ij',tr_w_ij,jdim=5)
        call read_dist_data(grid%dd2d,fid,'trsnowbv0',trsnowbv0,jdim=4)
#endif
      end select
      return
      end subroutine par_io_soils

      subroutine def_rsf_snow(fid)
      use ghy_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,dzsn_ij,
     &     'dzsn_ij(nlsn,bv,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,wsn_ij,
     &     'wsn_ij(nlsn,bv,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,hsn_ij,
     &     'hsn_ij(nlsn,bv,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,nsn_ij,
     &     'nsn_ij(bv,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,fr_snow_ij,
     &     'fr_snow_ij(bv,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid%dd2d,fid,tr_wsn_ij,
     &     'tr_wsn_ij(ntm,nlsn,bv,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_snow

      subroutine par_io_snow(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      use ghy_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'dzsn_ij', dzsn_ij,jdim=4)
        call write_dist_data(grid%dd2d, fid, 'wsn_ij', wsn_ij,jdim=4)
        call write_dist_data(grid%dd2d, fid, 'hsn_ij', hsn_ij,jdim=4)
        call write_dist_data(grid%dd2d, fid, 'nsn_ij', nsn_ij, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'fr_snow_ij',
     &       fr_snow_ij,jdim=3)
#ifdef TRACERS_WATER
        call write_dist_data(grid%dd2d, fid, 'tr_wsn_ij',
     &       tr_wsn_ij,jdim=5)
#endif
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'dzsn_ij', dzsn_ij,jdim=4)
        call read_dist_data(grid%dd2d, fid, 'wsn_ij', wsn_ij,jdim=4)
        call read_dist_data(grid%dd2d, fid, 'hsn_ij', hsn_ij,jdim=4)
        call read_dist_data(grid%dd2d, fid, 'nsn_ij', nsn_ij, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'fr_snow_ij',
     &       fr_snow_ij,jdim=3)
#ifdef TRACERS_WATER
        call read_dist_data(grid%dd2d, fid, 'tr_wsn_ij',
     &       tr_wsn_ij,jdim=5)
#endif
      end select
      return      
      end subroutine par_io_snow

      subroutine def_rsf_vegetation(fid)
      use veg_com, only : Cint, Qfol, cnc_ij
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,cint,'cint(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,qfol,'qfol(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,cnc_ij,'cnc_ij(dist_im,dist_jm)')
      return
      end subroutine def_rsf_vegetation

      subroutine par_io_vegetation(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      use veg_com, only : Cint, Qfol, cnc_ij
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'cint', cint)
        call write_dist_data(grid%dd2d, fid, 'qfol', qfol)
        call write_dist_data(grid%dd2d, fid, 'cnc_ij', cnc_ij)
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'cint', cint)
        call read_dist_data(grid%dd2d, fid, 'qfol', qfol)
        call read_dist_data(grid%dd2d, fid, 'cnc_ij', cnc_ij)
      end select
      return
      end subroutine par_io_vegetation

      subroutine def_rsf_seaice(fid)
      use seaice_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,rsi,'rsi(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,snowi,'snowi(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,msi,'msi(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,pond_melt,'pond_melt(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,flag_dsws,'flag_dsws(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,hsi,'hsi(lmi,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,ssi,'ssi(lmi,dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid%dd2d,fid,trsi,'trsi(ntm,lmi,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_seaice

      subroutine par_io_seaice(fid,iaction)
      use model_com, only : ioread,iowrite
      use seaice_com
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid               !@var fid unit number of read/write
      integer iaction           !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
      call write_dist_data(grid%dd2d, fid, 'rsi', rsi)
      call write_dist_data(grid%dd2d, fid, 'snowi', snowi)
      call write_dist_data(grid%dd2d, fid, 'msi', msi)
      call write_dist_data(grid%dd2d, fid, 'pond_melt', pond_melt)
      call write_dist_data(grid%dd2d, fid, 'flag_dsws', flag_dsws)
      call write_dist_data(grid%dd2d, fid, 'hsi', hsi, jdim=3)
      call write_dist_data(grid%dd2d, fid, 'ssi', ssi, jdim=3)
#ifdef TRACERS_WATER
      call write_dist_data(grid%dd2d, fid, 'trsi', trsi, jdim=4)
#endif
      case (ioread)             ! input from restart file
#ifdef CUBE_GRID
      call read_dist_data(grid%dd2d, fid, 'rsi', rsi,jdim=3)
      write(*,*) "DEBUG after read rsi"
      call read_dist_data(grid%dd2d, fid, 'snowi', snowi, jdim=3)
      call read_dist_data(grid%dd2d, fid, 'msi', msi, jdim=3)
      call read_dist_data(grid%dd2d, fid, 'pond_melt', pond_melt,
     &     jdim=3)
      call read_dist_data(grid%dd2d, fid, 'hsi', hsi, jdim=4)
      call read_dist_data(grid%dd2d, fid, 'ssi', ssi, jdim=4)
c      call read_dist_data(grid%dd2d, fid, 'flag_dsws' 
c     &     ,flag_dsws,jdim=3)

#else 
      call read_dist_data(grid%dd2d, fid, 'rsi', rsi)
     
      call read_dist_data(grid%dd2d, fid, 'snowi', snowi)
      call read_dist_data(grid%dd2d, fid, 'msi', msi)
      call read_dist_data(grid%dd2d, fid, 'pond_melt', pond_melt)
      call read_dist_data(grid%dd2d, fid, 'flag_dsws', flag_dsws)
      call read_dist_data(grid%dd2d, fid, 'hsi', hsi, jdim=3)
      call read_dist_data(grid%dd2d, fid, 'ssi', ssi, jdim=3)
#endif
#ifdef TRACERS_WATER
      call read_dist_data(grid%dd2d, fid, 'trsi', trsi, jdim=4)
#endif
      end select
      return
      end subroutine par_io_seaice

      subroutine def_rsf_landice(fid)
      use landice_com
      use landice
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,snowli,'snowli(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tlandi,'tlandi(d2,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,mdwnimp,'mdwnimp(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,edwnimp,'edwnimp(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,accpda,'accpda')
      call defvar(grid%dd2d,fid,eaccpda,'eaccpda')
      call defvar(grid%dd2d,fid,accpdg,'accpdg')
      call defvar(grid%dd2d,fid,eaccpdg,'eaccpdg')
#ifdef TRACERS_WATER
      call defvar(grid%dd2d,fid,trsnowli,
     &     'trsnowli(ntm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,trlndi,'trlndi(ntm,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,trdwnimp,
     &     'trdwnimp(ntm,dist_im,dist_jm)')
#ifdef TRACERS_OCEAN
      call defvar(grid%dd2d,fid,traccpda,'traccpda(ntm)')
      call defvar(grid%dd2d,fid,traccpdg,'traccpdg(ntm)')
#endif
#endif
      return
      end subroutine def_rsf_landice

      subroutine par_io_landice(fid,iaction)
      use model_com, only : ioread,iowrite
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      use landice_com
      use landice
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d,fid,'snowli',snowli)
        call write_dist_data(grid%dd2d,fid,'tlandi',tlandi,jdim=3)
        call write_dist_data(grid%dd2d,fid,'mdwnimp',mdwnimp)
        call write_dist_data(grid%dd2d,fid,'edwnimp',edwnimp)
        call write_data(grid%dd2d,fid,'accpda',accpda)
        call write_data(grid%dd2d,fid,'eaccpda',eaccpda)
        call write_data(grid%dd2d,fid,'accpdg',accpdg)
        call write_data(grid%dd2d,fid,'eaccpdg',eaccpdg)
#ifdef TRACERS_WATER
        call write_dist_data(grid%dd2d,fid,'trsnowli',trsnowli,jdim=3)
        call write_dist_data(grid%dd2d,fid,'trlndi',trlndi,jdim=3)
        call write_dist_data(grid%dd2d,fid,'trdwnimp',trdwnimp,jdim=3)
#ifdef TRACERS_OCEAN
        call write_data(grid%dd2d,fid,'traccpda',traccpda)
        call write_data(grid%dd2d,fid,'traccpdg',traccpdg)
#endif
#endif
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d,fid,'snowli',snowli)
        call read_dist_data(grid%dd2d,fid,'tlandi',tlandi,jdim=3)
c set some defaults for quantities which may not be in the
c restart file
        mdwnimp(:,:) = 0.; edwnimp(:,:) = 0.
        accpda = 0.; eaccpda = 0.
        accpdg = 0.; eaccpdg = 0.
        call read_dist_data(grid%dd2d,fid,'mdwnimp',mdwnimp)
        call read_dist_data(grid%dd2d,fid,'edwnimp',edwnimp)
        call read_data(grid%dd2d,fid,'accpda',accpda,bcast_all=.true.)
        call read_data(grid%dd2d,fid,'eaccpda',eaccpda,bcast_all=.true.)
        call read_data(grid%dd2d,fid,'accpdg',accpdg,bcast_all=.true.)
        call read_data(grid%dd2d,fid,'eaccpdg',eaccpdg,bcast_all=.true.)
#ifdef TRACERS_WATER
        call read_dist_data(grid%dd2d,fid,'trsnowli',trsnowli,jdim=3)
        call read_dist_data(grid%dd2d,fid,'trlndi',trlndi,jdim=3)
        call read_dist_data(grid%dd2d,fid,'trdwnimp',trdwnimp,jdim=3)
#ifdef TRACERS_OCEAN
        call read_data(grid%dd2d,fid,'traccpda',traccpda,
     &       bcast_all=.true.)
        call read_data(grid%dd2d,fid,'traccpdg',traccpdg,
     &       bcast_all=.true.)
#endif
#endif
      end select
      return
      end subroutine par_io_landice

      subroutine def_rsf_lakes(fid)
      use lakes_com
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,mldlk,'mldlk(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,mwl,'mwl(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,tlake,'tlake(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,gml,'gml(dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,flake,'flake(dist_im,dist_jm)')
#ifdef TRACERS_WATER
      call defvar(grid%dd2d,fid,trlake,'trlake(ntm,d2,dist_im,dist_jm)')
#endif
      return
      end subroutine def_rsf_lakes

      subroutine par_io_lakes(fid,iaction)
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      use lakes_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'mldlk', mldlk)
        call write_dist_data(grid%dd2d, fid, 'mwl',   mwl)
        call write_dist_data(grid%dd2d, fid, 'tlake', tlake)
        call write_dist_data(grid%dd2d, fid, 'gml',   gml)
        call write_dist_data(grid%dd2d, fid, 'flake', flake)
#ifdef TRACERS_WATER
        call write_dist_data(grid%dd2d, fid, 'trlake', trlake, jdim=4)
#endif
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'mldlk', mldlk)
        call read_dist_data(grid%dd2d, fid, 'mwl',   mwl)
        call read_dist_data(grid%dd2d, fid, 'tlake', tlake)
        call read_dist_data(grid%dd2d, fid, 'gml',   gml)
        call read_dist_data(grid%dd2d, fid, 'flake', flake)
#ifdef TRACERS_WATER
        call read_dist_data(grid%dd2d, fid, 'trlake', trlake, jdim=4)
#endif
      end select
      return
      end subroutine par_io_lakes

      subroutine def_rsf_ocean(fid)
      use static_ocean
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,tocean,'tocean(d3,dist_im,dist_jm)')
      call defvar(grid%dd2d,fid,z1o,'z1o(dist_im,dist_jm)')
      return
      end subroutine def_rsf_ocean

      subroutine par_io_ocean(fid,iaction)
      use model_com, only : ioread,iowrite
      use static_ocean
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_dist_data(grid%dd2d, fid, 'tocean', tocean, jdim=3)
        call write_dist_data(grid%dd2d, fid, 'z1o', z1o)
      case (ioread)            ! input from restart file
        call read_dist_data(grid%dd2d, fid, 'tocean', tocean, jdim=3)
        call read_dist_data(grid%dd2d, fid, 'z1o', z1o)
      end select
      return
      end subroutine par_io_ocean

      subroutine def_rsf_acc(fid)
      use model_com, only : idacc
      use diag_com, only : tsfrez=>tsfrez_loc,aij=>aij_loc,
     &     ail=>ail_loc,aijk=>aijk_loc,oa,tdiurn,speca,atpe,
     &     keyct,keynr,atpe,adiurn,hdiurn,energy,wave,aisccp
      use domain_decomp, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid%dd2d,fid,idacc,'idacc(nsampl)')
      call defvar(grid%dd2d,fid,keyct,'keyct')
      call defvar(grid%dd2d,fid,keynr,'keynr(nkeynr,nkeymo)')
      call defvar(grid%dd2d,fid,tsfrez,'tsfrez(dist_im,dist_jm,ktsf)')
      call defvar(grid%dd2d,fid,tdiurn,'tdiurn(dist_im,dist_jm,ktd)')
      call defvar(grid%dd2d,fid,oa,'oa(dist_im,dist_jm,koa)')
      call defvar(grid%dd2d,fid,aij,'aij(dist_im,dist_jm,kaij)')
      call defvar(grid%dd2d,fid,ail,'ail(dist_im,dist_jm,lm,kail)')
      call defvar(grid%dd2d,fid,aijk,'aijk(dist_im,dist_jm,lm,kaijk)')
      call defvar(grid%dd2d,fid,energy,'energy(nehist,hist_days)')
      call defvar(grid%dd2d,fid,speca,
     &     'speca(imh_plus_1,kspeca,nspher)')
      call defvar(grid%dd2d,fid,atpe,'atpe(ktpe,nhemi)')
      call defvar(grid%dd2d,fid,wave,
     &     'wave(re_and_im,max12hr_sequ,nwav_dag,kwp)')
      call defvar(grid%dd2d,fid,aisccp,'aisccp(ntau,npres,nisccp)')
      call defvar(grid%dd2d,fid,adiurn,
     &     'adiurn(ndiuvar,ndiupt,hr_in_day)')
#ifndef NO_HDIURN
      call defvar(grid%dd2d,fid,hdiurn,
     &     'hdiurn(ndiuvar,ndiupt,hr_in_month)')
#endif

c problematic arrays: need to get them from latlon version of diag_com
c     *     AJ,AREG,APJ,AJL,ASJL,CONSRV,AJK

      return
      end subroutine def_rsf_acc

      subroutine par_io_acc(fid,iaction)
      use model_com, only : ioread,iowrite,idacc
      use diag_com, only : tsfrez=>tsfrez_loc,aij=>aij_loc,
     &     ail=>ail_loc,aijk=>aijk_loc,oa,tdiurn,speca,atpe,
     &     keyct,keynr,atpe,adiurn,hdiurn,energy,wave,aisccp
      use domain_decomp, only : grid
      use pario, only : write_dist_data,read_dist_data
     &     ,write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to standard restart file
        call write_data(grid%dd2d,fid,'idacc',idacc)
        call write_data(grid%dd2d,fid,'keyct',keyct)
        call write_data(grid%dd2d,fid,'keynr',keynr)
        call write_data(grid%dd2d,fid,'energy',energy)
        call write_data(grid%dd2d,fid,'speca',speca)
        call write_data(grid%dd2d,fid,'atpe',atpe)
        call write_data(grid%dd2d,fid,'wave',wave)
        call write_data(grid%dd2d,fid,'aisccp',aisccp)
        call write_data(grid%dd2d,fid,'adiurn',adiurn)
#ifndef NO_HDIURN
        call write_data(grid%dd2d,fid,'hdiurn',hdiurn)
#endif
        call write_dist_data(grid%dd2d,fid,'tsfrez',tsfrez)
        call write_dist_data(grid%dd2d,fid,'tdiurn',tdiurn)
        call write_dist_data(grid%dd2d,fid,'oa',oa)
        call write_dist_data(grid%dd2d,fid,'aij',aij)
        call write_dist_data(grid%dd2d,fid,'ail',ail)
        call write_dist_data(grid%dd2d,fid,'aijk',aijk)
      case (ioread)            ! input from restart file
c for which scalars is bcast_all=.true. necessary?
        call read_data(grid%dd2d,fid,'idacc',idacc)
        call read_data(grid%dd2d,fid,'keyct',keyct)
        call read_data(grid%dd2d,fid,'keynr',keynr)
        call read_data(grid%dd2d,fid,'energy',energy)
        call read_data(grid%dd2d,fid,'speca',speca)
        call read_data(grid%dd2d,fid,'atpe',atpe)
        call read_data(grid%dd2d,fid,'wave',wave)
        call read_data(grid%dd2d,fid,'aisccp',aisccp)
        call read_data(grid%dd2d,fid,'adiurn',adiurn)
#ifndef NO_HDIURN
        call read_data(grid%dd2d,fid,'hdiurn',hdiurn)
#endif
        call read_dist_data(grid%dd2d,fid,'tsfrez',tsfrez)
        call read_dist_data(grid%dd2d,fid,'tdiurn',tdiurn)
        call read_dist_data(grid%dd2d,fid,'oa',oa)
        call read_dist_data(grid%dd2d,fid,'aij',aij)
        call read_dist_data(grid%dd2d,fid,'ail',ail)
        call read_dist_data(grid%dd2d,fid,'aijk',aijk)
      end select
      return
      end subroutine par_io_acc
