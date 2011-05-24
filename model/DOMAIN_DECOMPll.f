#include "rundeck_opts.h"

      MODULE DOMAIN_DECOMP_ATM
        use dist_grid_mod
        use Halo_mod
        use SpecialIO_mod
        use GatherScatter_mod
        use GlobalSum_mod
      implicit none
      public
      type(dist_grid), target :: grid
      END MODULE DOMAIN_DECOMP_ATM
