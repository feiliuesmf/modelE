#include "rundeck_opts.h"

module PlanetParams_mod
#ifdef PLANET_PARAMS
!@type PlanetParams_t a derived type for the collection of
!@+    parameters users may want to change when simulating
!@+    environments other than Earth.  When this module is
!@+    in use, the corresponding parameters in the constants
!@+    module take their values from one of the instances of
!@+    PlanetParams_t declared below.
  type PlanetParams_t
!@param name name of planet being simulated.
!@+   Not intended to be widely used to control physics settings etc.
!@+   May be removed in future.
     character(len=16) :: name
!@param sday  sec per day (s)
     real*8 :: sday
!@param omega rotation rate (1/s)
     real*8 :: omega
! edperd, edpery only seem to be used to set omega.  For now, setting omega directly.
!     real*8 :: edperd,edpery
!@param grav gravitational acceleration (m/s2)
     real*8 :: grav
!@param radius planetary radius (m)
     real*8 :: radius
!@param mair molar mass of "dry" (fixed-composition) air (g/mol)
     real*8 :: mair
!@param srat ratio of specific heats at const. press. and vol. (for fixed-composition air)
     real*8 :: srat
!@param psf (hPa) global mean suface pressure
     real*8 :: psf
  end type PlanetParams_t

!
! Declare the values of parameters for various planets.
!
  type(PlanetParams_t), parameter, private :: &
       notEarth = &
       PlanetParams_t( &
       name = 'notEarth', &
       sday = 86400d0, &
       omega = (2d0*3.1415926535897932d0)*366d0/(365d0*86400d0), &
       grav = 9.80665d0, &
       radius = 6371000d0, &
       mair = 28.9655d0, &
       srat = 1.401d0, &
       psf = 984d0 &
       )
!

!@param PlanetParams the instance of PlanetParams_t to be used. CPP selects
!@+     which instance of PlanetParams_t is copied into PlanetParams.
  type(PlanetParams_t), parameter :: PlanetParams = PLANET_PARAMS
#endif /* PLANET_PARAMS defined */
end module PlanetParams_mod
