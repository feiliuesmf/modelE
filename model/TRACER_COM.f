#include "rundeck_opts.h"

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@+    This version for simple trace gases
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
!@ver  1.0
      USE QUSDEF, only: nmom
      USE MODEL_COM, only: im,jm,lm
      USE CONSTANT, only: mair
      IMPLICIT NONE
      SAVE
!@param NTM number of tracers
      integer, parameter :: ntm=5

C**** Each tracer has a variable name and a unique index
!@var TRNAME: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<
      character*8, parameter :: trname(ntm)=
     * (/ 'Air     ','SF6     ','Rn222   ','CO2     ','14CO2   '/)
!@var N_XXX: variable names of indices for tracers
      integer, parameter :: 
     *     n_air=1,     n_SF6=2, n_Rn222=3,  n_CO2=4,  n_14CO2=5
!@var NTM_POWER: Power of 10 associated with each tracer (for printing)
      integer, parameter, dimension(ntm) :: ntm_power=
     *     (/   -2,         -14,       -21,       -6,  -18 /)
!@var TR_MM: molecular mass of each tracer (g/mole)
      real*8, parameter, dimension(ntm) :: tr_mm=
     *     (/ mair,    146.01d0,    222.d0,    44.d0,   46.d0 /)
!@var T_QLIMIT: if t_qlimit=.true. tracer is maintained as positive
      logical, parameter, dimension (ntm) :: t_qlimit=
     *     (/ .true.,   .true.,    .true.,  .false.,    .true. /)
!@var needtrs: true if surface tracer value from PBL is required
      logical, parameter, dimension (ntm) :: needtrs=
     *     (/.false.,  .false.,   .false.,  .false.,    .false. /)
!@var trdecy radioactive decay constant (1/s) (=0 for stable tracers)
      real*8, parameter, dimension (ntm) :: trdecy=
     *     (/    0d0,      0d0,    2.1d-6,      0d0,     0d0 /)

C**** Note units for these parameters!
C**** Example: clay dust; trpdens=2.5d3, trradius=0.73d-6 
C****          silt dust; trpdens=2.65d3, trradius=6.1d-6 
C****
!@var trpdens tracer particle density (kg/m^3) (=0 for non-particle tracers)
      real*8, parameter, dimension (ntm) :: trpdens=
     *     (/    0d0,      0d0,    0d0,      0d0,     0d0 /)
!@var trradius tracer effective radius (m) (=0 for non particle tracers)
      real*8, parameter, dimension (ntm) :: trradius=
     *     (/    0d0,      0d0,    0d0,      0d0,     0d0 /)

!@var ITIME_TR0: start time for each tracer (hours)
      integer, dimension(ntm) :: itime_tr0=0.
!@var TRM: Tracer array (kg)
      real*8, dimension(im,jm,lm,ntm) :: trm
!@var TRMOM: Second order moments for tracers (kg)
      real*8, dimension(nmom,im,jm,lm,ntm) :: trmom
!@var MTRACE: timing index for tracers
      integer mtrace

#ifdef TRACERS_WATER
!@param nWD_TYPES number of tracer types for wetdep purposes
!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water (original condense method)
      integer, parameter :: nWD_TYPES=3, nGAS=1, nPART=2, nWATER=3
!@var tr_wd_TYPE: tracer wet dep type (gas, particle, water)
      integer, parameter, dimension (ntm) :: tr_wd_TYPE=
     *     (/    nGAS, nGAS, nGAS, nGAS  /)
!@var tr_RKD: Henry's Law coefficient (JOULE/mole !)
      real*8, parameter, dimension (ntm) :: tr_RKD=
     *     (/   0.d0, 0.d0, 0.d0, 0.d0 /)
!@var tr_DHD: coefficient of temperature-dependence term of Henry's
!@+   Law coefficient (joule/mole !)
      real*8, parameter, dimension (ntm) :: tr_DHD=
     *     (/   0.d0 , 0.d0, 0.d0 , 0.d0 /)
!@param tr_evap_fact fraction of re-evaporation by tracer type
C note, tr_evap_fact is not dimensioned as NTM:
      REAL*8, parameter, dimension (nWD_TYPES) :: tr_evap_fact=
     *     (/1.d0, 0.5d0,  1.d0/)
!@var TRWM tracer in cloud liquid water amount (kg)
      real*8, dimension(im,jm,lm,ntm) :: trwm
!@var TRW0 default tracer concentration in water (kg/kg)
      real*8, dimension(ntm) :: trw0 = 1.
#endif

!@var ntsurfsrcmax maximum number of surface sources/sinks
      integer, parameter :: ntsurfsrcmax=10
!@var ntsurfsrc no. of non-interactive surface sources for each tracer
      integer, parameter, dimension(ntm) :: ntsurfsrc = (/
     *     0,       1,       1,       6,      1 /)

!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: nt3Dsrcmax=1
!@var nt3Dsrc no. of 3D sources for each tracer
      integer, parameter, dimension(ntm) :: nt3Dsrc = (/
     *     0,       0,       0,       0,      1 /)
    
    
      END MODULE TRACER_COM

