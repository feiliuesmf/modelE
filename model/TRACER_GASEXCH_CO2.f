      MODULE TRACER_GASEXCH_COM

      USE TRACER_COM, only : ntm    !tracers in air-sea gas exch

      implicit none

      SAVE


#include "dimensions.h"
#include "dimension2.h"


      SAVE


      real*8 atracflx(iia,jja,ntm),atrac(iia,jja,ntm)
      common /gasexch1/atracflx,atrac

      real*8   tracflx(idm,jdm,ntm)     !  tracer flux at air-sea intfc
      common /gasexch2/tracflx
  
      real*8 tracflx1d(ntm)

      END MODULE TRACER_GASEXCH_COM

c used with TRACERS_GASEXCH_CO2 to compute transfer velocity for CO2
c
c taken from Watsons' code
c ---------------------------------------------------------------------
      REAL FUNCTION sol_co2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c-------------------------------------------------------------------

      REAL*8    pt,ps,ta,tk100,tk1002


      ta = (pt + 273.16d0)
      tk100 = ta*0.01d0
      tk1002 = tk100*tk100
      sol_co2 = exp(-162.8301d0 + 218.2968d0/tk100  +
     .          90.9241d0*log(tk100) - 1.47696d0*tk1002 +
     .          ps * (.025695d0 - .025225d0*tk100 +
     .          0.0049867d0*tk1002))

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      REAL FUNCTION sc_co2(t)
c---------------------------------------------------
c     CO2 Schmidt number 
c     as a function of temperature. 
c
c     t: temperature (degree Celcius)
c---------------------------------------------------
      IMPLICIT NONE 
      REAL*8  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 t
c
      sc_co2 = 2073.1 - 125.62*t + 3.6276*t**2 - 0.043219*t**3
c
      RETURN 
      END 
