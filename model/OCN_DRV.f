#include "rundeck_opts.h"

      subroutine ocean_driver
      use model_com, only : itime,nday
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use seaice_com, only : si_ocn,iceocn
#ifdef CUBED_SPHERE
      use seaice_com, only : si_atm    ! temporary
#endif
      use fluxes, only : atmocn,atmice
      use icedyn_com, only : igice
      implicit none
      call startTimer('OCEANS')
C**** CALCULATE ICE DYNAMICS
#ifdef CUBED_SPHERE
! Temporary fix until seaice moves to ocean grid:
! to match previous results, using ice cover/amounts including
! lakes. This matters along coastlines (gradient calcs), but the
! earlier results are not more accurate.
      call seaice_to_atmgrid(atmice)
      si_ocn%rsisave(:,:) = si_ocn%rsi(:,:)
      CALL DYNSI(atmice,iceocn,si_atm)
#else
      CALL DYNSI(atmice,iceocn,si_ocn)
#endif
C**** CALCULATE BASE ICE-OCEAN/LAKE FLUXES
      CALL UNDERICE(si_ocn,iceocn,atmocn)
C**** APPLY SURFACE/BASE FLUXES TO SEA/LAKE ICE
      CALL GROUND_SI(si_ocn,iceocn,atmice,atmocn)
         CALL CHECKT ('GRNDSI')
C**** set total atmopsheric pressure anomaly in case needed by ocean
      CALL CALC_APRESS(atmice)
C**** APPLY FLUXES TO OCEAN, DO OCEAN DYNAMICS AND CALC. ICE FORMATION
      CALL OCEANS(atmocn,iceocn,igice)
         CALL CHECKT ('OCEANS')
C**** APPLY ICE FORMED IN THE OCEAN/LAKES TO ICE VARIABLES
      CALL FORM_SI(si_ocn,iceocn,atmice)
      call seaice_to_atmgrid(atmice) ! needed only to preserve former result
         CALL CHECKT ('FORMSI')
C**** ADVECT ICE
      CALL ADVSI(atmice)
         CALL CHECKT ('ADVSI ')
      CALL SI_diags(si_ocn,iceocn,atmice)


      call stopTimer('OCEANS')
      end subroutine ocean_driver

      SUBROUTINE INPUT_ocean (istart,istart_fixup,
     &     do_IC_fixups,is_coldstart)
      use fluxes, only : atmocn,atmice
      use icedyn_com, only : igice
      implicit none
!@var istart start(1-8)/restart(>8)  option
      integer :: istart,istart_fixup,do_IC_fixups
      LOGICAL :: is_coldstart

      logical :: iniOcean

      iniOcean = is_coldstart

C**** Initialize sea ice
      if(istart.eq.2) call read_seaice_ic
      CALL init_oceanice(iniOCEAN,do_IC_fixups,atmocn)
      !call seaice_to_atmgrid(atmice)

C**** Initialize ice dynamics code (if required)
      CALL init_icedyn(iniOCEAN,atmice)

C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN(iniOCEAN,istart_fixup,atmocn
     &     ,igice) ! not necessary now that dynsi uosurf,vosurf in rsf?

      return
      end subroutine INPUT_ocean

      subroutine alloc_drv_ocean
c Driver to allocate arrays that become dynamic as a result of
c set-up for MPI implementation
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      USE SEAICE_COM, only : si_atm,si_ocn,sigrid
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      USE TRACER_COM, only : ntm
#endif
#if (defined TRACERS_OCEAN) && !defined(TRACERS_OCEAN_INDEP)
      USE TRACER_COM, only : trname
      USE OCN_TRACER_COM, only :
     &     ntm_ocn    => ntm,
     &     trname_ocn => trname
#endif
#if (defined TRACERS_WATER)
      USE SEAICE, only : ntm_si=>ntm
#endif
      IMPLICIT NONE

#if (defined TRACERS_OCEAN) && !defined(TRACERS_OCEAN_INDEP)
! copy atmosphere-declared tracer info to ocean so that the ocean
! can "inherit" it without referencing atm. code
      ntm_ocn = ntm
! trname is copied here rather than in init_tracer since it is
! needed immediately when reading checkpoint files
      allocate(trname_ocn(ntm))
      trname_ocn(:) = trname(:)
#endif

#if (defined TRACERS_WATER)
! copy atmosphere-declared tracer info to seaice so that the seaice
! can "inherit" it without referencing atm. code
      ntm_si = ntm
      si_ocn % ntm = ntm
#endif

      call alloc_icedyn(grid%im_world,grid%jm_world)
      call alloc_icedyn_com(grid)
      call alloc_seaice_com(grid)
      si_atm%grid => grid
      si_ocn%grid => grid
      sigrid => grid
      call alloc_ocean

      end subroutine alloc_drv_ocean
