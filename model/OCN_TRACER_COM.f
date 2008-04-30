#include "rundeck_opts.h"

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      MODULE OCN_TRACER_COM
!@sum OCN_TRACER_COM: sets up tracer quantities for ocean tracers
!@+   This information can either be set up directly from the AGCM
!@+   or can be indpendently defined here.
!@+   Thus the ocean can have either more tracers or less than the AGCM 
!@+   depending on application.
!@+   The default behaviour will be to have the same number of 
!@+   TRACERS_WATER, but TRACERS_OCEAN will be independent.
!@+   Use "TRACERS_OCEAN" to allow ocean tracers.
!@+   Use "TRACERS_WATER" for freshwater tracers from ATM (this can be
!@+         used indpendently from TRACERS_OCEAN if a surface boundary 
!@+         condition is all that is required)
!@+   Use "TRACERS_OCEAN_INDEP" for independently defined ocn tracers
!@+        "TRACERS_AGE_OCEAN" is one partciularly case
#ifdef TRACERS_OCEAN_INDEP
C**** this defines tracer parameters that are local to ocean code
#ifdef TRACERS_AGE_OCEAN
      INTEGER, PARAMETER :: ntm=1
      CHARACTER*10 :: trname(ntm) = (/ '       Age'/)
      REAL*8, DIMENSION(ntm) :: trw0=0, trdecay=0
#else
      INTEGER, PARAMETER :: ntm=1
      CHARACTER*10 :: trname(ntm) = (/ '     Water'/)
      REAL*8, DIMENSION(ntm) :: trw0=1d0, trdecay=0
#endif
      LOGICAL, DIMENSION(ntm) :: t_qlimit=.true. 
      INTEGER, DIMENSION(ntm) :: itime_tr0 = 0, ntrocn = 0
      INTEGER :: n_age = 0, n_water = 0
#else
C**** use only agcm data
      USE TRACER_COM, only : ntm, trname, itime_tr0, trdecay,
     *        t_qlimit, trw0, ntrocn, n_age, n_water
#endif
 
      END MODULE OCN_TRACER_COM
#endif
