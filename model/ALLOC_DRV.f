#include "rundeck_opts.h"
      subroutine alloc_drv()
c Driver to allocate arrays that become dynamic as a result of
c set-up for MPI implementation
      USE DOMAIN_DECOMP, ONLY : grid
      IMPLICIT NONE

      call alloc_dynamics(grid)
      call alloc_model_com(grid)
      call alloc_fluxes(grid)
      call alloc_clouds_com(grid)
      call alloc_ghy_com(grid)
      call alloc_pbl_com(grid)
      call alloc_icedyn(grid)
      call alloc_icedyn_com(grid)
      call alloc_diag_com(grid)
      call alloc_diag_loc(grid)
      call alloc_smomtq(grid)
      call alloc_strat_com(grid)
      call alloc_seaice_com(grid)
      call alloc_rad_com(grid)
      call alloc_lakes(grid)
      call alloc_lakes_com(grid)
      call alloc_landice_com(grid)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call alloc_tracer_com(grid)
#ifdef TRACERS_DRYDEP
      call alloc_trdrydep(grid)
#endif
#ifdef TRACERS_SPECIAL_Lerner
      call alloc_tracer_special_lerner_com(grid)
      call alloc_linoz_chem_com(grid)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call alloc_trchem_shindell_com(grid)
      call alloc_tracer_sources(grid)
      call alloc_lightning(grid)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      call alloc_aerosol_sources(grid)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      CALL alloc_dust(grid)
#endif
#endif
      call alloc_tracer_adv(grid)
#ifdef USE_ENT
!!! should be done in init_module_ent
      call alloc_ent_com(grid)
#else
      call alloc_veg_com(grid)
#endif
      call alloc_ocean(grid)
#ifdef TRACERS_ON
      call alloc_trdiag_com
#endif

      end subroutine alloc_drv
