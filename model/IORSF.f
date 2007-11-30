#include "rundeck_opts.h"

      SUBROUTINE io_rsf(kunit,it,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@ver   1.0
!@calls io_model,io_ocean,io_lakes,io_seaice,io_earth,io_soils,io_snow
!@+     io_landice,io_bldat,io_pbl,io_clouds,io_somtq,io_rad,io_diags
!@+     io_ocdiag,io_icedyn,io_icdiag
      USE DOMAIN_DECOMP, only : REWIND_PARALLEL
      USE MODEL_COM, only : ioread_single,iowrite_single,Kradia
     *                     ,ioread,ioread_nodiag 

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

      iact=iaction   ; skip_diag=.false.
      if(iaction==ioread_nodiag) then
         iact=ioread ; skip_diag=.true.
      end if

      ioerr=-1
      if (kunit.gt.0) CALL REWIND_PARALLEL( kunit )

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
C**** Particular values may produce variations in indiv. i/o routines

C**** Calls to individual i/o routines
      call io_label  (kunit,it,itm,iact,ioerr)
      it1=it
      if (Kradia.gt.0) then
        if (iaction.ne.ioread_single .and.
     *   iaction.ne.iowrite_single) call io_rad (kunit,iact,ioerr)
        call io_diags  (kunit,it,iact,ioerr)
        if(it1.ne.it .or. ioerr.eq.1)
     &       call stop_model('restart problem',255)
        return
      end if
      if(iaction.ne.ioread_single.and.iaction.ne.iowrite_single) then
        call io_model  (kunit,iact,ioerr) 
        call io_strat  (kunit,iact,ioerr)
        call io_ocean  (kunit,iact,ioerr)
        call io_lakes  (kunit,iact,ioerr)
        call io_seaice (kunit,iact,ioerr)
        call io_earth  (kunit,iact,ioerr)
        call io_soils  (kunit,iact,ioerr)
        call io_vegetation  (kunit,iact,ioerr)
#ifdef USE_ENT
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
        call io_veg_related  (kunit,iact,ioerr)
        !call io_ent    (kunit,iaction,ioerr)
#endif
        call io_snow   (kunit,iact,ioerr)
        call io_landice(kunit,iact,ioerr)
        call io_bldat  (kunit,iact,ioerr)
        call io_pbl    (kunit,iact,ioerr)
        call io_clouds (kunit,iact,ioerr)
        call io_somtq  (kunit,iact,ioerr)
        call io_rad    (kunit,iact,ioerr)
        call io_icedyn (kunit,iact,ioerr)
#ifdef TRACERS_ON
        call io_tracer (kunit,iact,ioerr)
#endif
#ifdef NUDGE_ON
        call io_nudge (kunit,iaction,ioerr)
#endif
      end if
      if(skip_diag) go to 10
      call io_diags  (kunit,it,iaction,ioerr)
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

      RETURN
      END SUBROUTINE io_rsf
