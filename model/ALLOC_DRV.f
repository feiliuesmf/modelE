#include "rundeck_opts.h"
      subroutine alloc_drv()
c Driver to allocate arrays that become dynamic as a result of
c set-up for MPI implementation
      USE DOMAIN_DECOMP, ONLY : grid
      IMPLICIT NONE

      call alloc_dynamics(grid)
      call alloc_model_com(grid)
c$$$      call alloc_fluxes(grid)
      call alloc_clouds_com(grid)
      call alloc_ghy_com(grid)
      call alloc_pbl_com(grid)
      call alloc_icedyn(grid)
      call alloc_icedyn_com(grid)
c      call alloc_dagcom(grid)
      call alloc_diag_loc(grid)
      call alloc_smomtq(grid)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call alloc_tracer_com(grid)
#endif

      end subroutine alloc_drv
