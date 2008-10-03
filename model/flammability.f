#include "rundeck_opts.h"

      subroutine calc_flammability(t,p,r,v,flam)
!@sum calculated the flammability of vegetation based on model
!@+ variables for temperature, precipitation, relative humidity,
!@+ and an index of vegetation density read from an input file.
!
!@auth Greg Faluvegi based on information from Olga Pechony
!@ver  1.0 (based on Olga's Word document Flammability.doc)
!
!@var T the surface air temperature passed in Kelvin for current I,J
!@var P the precipitation rate passed in mm/day for current I,J
!@var R the relative humidity passed as fraction for current I,J
!@var V vegetative density (unitless 0 to 1) for current I,J
!@var Z component of the Goff-Gratch saturation vapor pressure equation
!@var tsbyt = reciprocal of the temperature times ts
!@var flam the flammability returned for current I,J
!@param a,b,s,d,f,h,ts,cr coefficients for the parameterization
!@+ e.g. from Goff-Gratch saturation vapor pressure equation
!
      implicit none

      real*8, parameter :: a=-7.90298d0,d=11.344d0,c=-1.3816d-7,
     & b=5.02808,f=8.1328d-3,h=-3.49149d0,ts=373.16d0,cr=-2.d0
      real*8, intent(in) :: t,p,r,v
      real*8, intent(out) :: flam
      real*8 :: z,tsbyt

      tsbyt=ts/t

      z= a*(tsbyt-1.d0) + b*log10(tsbyt) + 
     &   c*(10.d0**(d*(1.d0-tsbyt))-1.d0) +
     &   f*(10.d0**(h*(tsbyt-1.d0))-1.d0)

      flam=(10.d0**z*(1.d0-r)) * exp(cr*p) * v

      return
      end subroutine calc_flammability
