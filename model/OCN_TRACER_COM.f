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
      CHARACTER*10 :: trname(ntm) = (/ 'Age       '/)
      REAL*8, DIMENSION(ntm) :: trw0=0, trdecay=0, ntrocn=0
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .false.
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
#else
#ifdef TRACERS_OceanBiology
#ifdef TRACERS_Alkalinity
      INTEGER, PARAMETER :: ntm=16
#else
      INTEGER, PARAMETER :: ntm=15
#endif
      CHARACTER*10 :: trname(ntm) = (/ 'Nitr      ', 'Ammo      ', 
     .     'Sili      ', 'Iron      ', 'Diat      ', 'Chlo      ', 
     .     'Cyan      ', 'Cocc      ', 'Herb      ', 'Inert     ',
     .     'N_det     ', 'S_det     ', 'I_det     ', 'DOC       ',
     .     'DIC       '
#ifdef TRACERS_Alkalinity
     .     ,'Alk       '
#endif
     .     /)
      REAL*8, DIMENSION(ntm) :: trw0=0, trdecay=0
      REAL*8  :: obio_tr_mm(ntm)= (/ 14., 14., 28.055, 55.845, 1., 1., 
     .     1., 1., 1., 14., 14., 28.055, 55.845, 12., 12.
#ifdef TRACERS_Alkalinity
     .     , 1.
#endif
     .     /)
!@var ntrocn scaling exponent for tracers
      INTEGER :: ntrocn(ntm) = (/ -4,-6,-4,-8,-8,-8,-8,-8,-8,
     .     -4,-6,-8,-10,-8,-2
#ifdef TRACERS_Alkalinity
     .     ,-6
#endif
     .     /)
!@dbparam to_per_mil For printout of tracer concentration in permil
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
!@param conc_from_fw definition for defining surface ocean conc
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .false.
#else
      INTEGER, PARAMETER :: ntm=1
      CHARACTER*10 :: trname(ntm) = (/ 'Water     '/)
      REAL*8, DIMENSION(ntm) :: trw0=1d0, trdecay=0, ntrocn=0
!@param conc_from_fw definition for defining surface ocean conc
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .true.
#endif  /* TRACERS_OceanBiology */
#endif  /*TRACERS_AGE_OCEAN */

      LOGICAL, DIMENSION(ntm) :: t_qlimit=.true. 
      INTEGER, DIMENSION(ntm) :: itime_tr0 = 0
      INTEGER :: n_water = 0

#else   /* not tracers_ocean_indep */

C**** use only agcm data
      USE TRACER_COM, only : ntm, trname, itime_tr0, trdecay,
     *        t_qlimit, trw0, ntrocn, n_water, conc_from_fw

#endif  /* TRACERS_OCEAN_INDEP */
      INTEGER :: n_age = 0, n_obio = 0

      END MODULE OCN_TRACER_COM
#endif
