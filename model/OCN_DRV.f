#include "rundeck_opts.h"

      subroutine ocean_driver
      use model_com, only : itime,nday
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      implicit none
      call startTimer('OCEANS')
C**** CALCULATE ICE DYNAMICS
      CALL DYNSI
C**** CALCULATE BASE ICE-OCEAN/LAKE FLUXES
      CALL UNDERICE('OCEAN')
C**** APPLY SURFACE/BASE FLUXES TO SEA/LAKE ICE
      CALL GROUND_SI('OCEAN')
         CALL CHECKT ('GRNDSI')
C**** APPLY FLUXES TO OCEAN, DO OCEAN DYNAMICS AND CALC. ICE FORMATION
      CALL OCEANS
         CALL CHECKT ('OCEANS')
C**** APPLY ICE FORMED IN THE OCEAN/LAKES TO ICE VARIABLES
      CALL FORM_SI('OCEAN')
         CALL CHECKT ('FORMSI')
C**** ADVECT ICE
      CALL ADVSI
      CALL ADVSI_DIAG ! needed to update qflux model, dummy otherwise
         CALL CHECKT ('ADVSI ')
C**** SAVE some noon GMT ice quantities
      IF (MOD(Itime+1,NDAY).ne.0 .and. MOD(Itime+1,NDAY/2).eq.0)
     &        call vflx_OCEAN

      call stopTimer('OCEANS')
      end subroutine ocean_driver

      SUBROUTINE INPUT_ocean (istart,istart_fixup,
     &     do_IC_fixups,is_coldstart)
      implicit none
!@var istart start(1-8)/restart(>8)  option
      integer :: istart,istart_fixup,do_IC_fixups
      LOGICAL :: is_coldstart

      logical :: iniOcean

      iniOcean = is_coldstart

C**** Initialize sea ice
      if(istart.eq.2) call read_seaice_ic
      CALL init_ice(iniOCEAN,do_IC_fixups)

C**** Initialize ice dynamics code (if required)
      CALL init_icedyn(iniOCEAN)

C**** Initialize ocean variables
C****  KOCEAN = 1 => ocean heat transports/max. mixed layer depths
C****  KOCEAN = 0 => RSI/MSI factor
      CALL init_OCEAN(iniOCEAN,istart_fixup)

      return
      end subroutine INPUT_ocean

      subroutine alloc_drv_ocean
c Driver to allocate arrays that become dynamic as a result of
c set-up for MPI implementation
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      IMPLICIT NONE

      call alloc_icedyn()
      call alloc_icedyn_com(grid)
      call alloc_seaice_com(grid)
      call alloc_ocean(grid)

      end subroutine alloc_drv_ocean
