#include "rundeck_opts.h"

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
      common /gasexch3/tracflx1d
!$OMP THREADPRIVATE(/gasexch3/)


      END MODULE TRACER_GASEXCH_COM

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------

      SUBROUTINE TRACERS_GASEXCH_CO2_Natassa_PBL(tg1,ws,
     . alati,psurf,itr,trconstflx,byrho,Kw_gas,alpha_gas,
     . beta_gas,trsf,trcnst,ilong,jlat)
  

      USE CONSTANT, only:    rhows,mair
      USE TRACER_COM, only : ntm,trname,tr_mm
      USE obio_incom, only : awan
      
      implicit none

      integer :: ilong,jlat,itr
      real*8  :: tg1,ws,psurf,trconstflx,byrho,trsf,trcnst
      real*8  :: alati,Kw_gas,alpha_gas,beta_gas
      real*8  :: Sc_gas
      real*8, external :: sc_co2,sol_co2

!@var  psurf surface pressure
!@var  alati SSS at i,j
!@var  Kw_gas  gas exchange transfer velocity at i,j only over ocean
!@var  alpha_gas  solubility of gas
!@var  beta_gas  conversion term  that includes solubility
!@var  trconstflx  constant component of surface tracer flux
!      byrho   1/surface air density

!routine for PBL calculations of CO2 gas exchange

      !---------------------------------------------------------------
      !TRANSFER VELOCITY
      !---------------------------------------------------------------
      !Schmidt number for gas
       Sc_gas=sc_co2(tg1)

      !wind speed ws: magn. of surf. wind modified by buoyancy flux (m/s)
      !compute transfer velocity Kw only over ocean
      if (Sc_gas .le. 0.) then
        write(*,'(a,2i4,a,2f9.3)')
     .          'warning: Sc_gas negtv, at ',ilong,jlat,
     .          ', Sc_gas,temp_c=',Sc_gas,tg1
         Kw_gas=1.e-10
      else
         Kw_gas=(Sc_gas/660.d0)**(-0.5d0) * ws * ws * awan !units of m/s
      endif

      !---------------------------------------------------------------
      !gas SOLUBILITY
      !---------------------------------------------------------------
      !alpha --solubility of CO2 in seawater
      !in mol/m^3/picoatm
       alpha_gas=sol_co2(tg1,alati)
      !convert to mol/m^3/atm
      !alpha_gas=alpha_gas*1.e+12   !mol/m3/pptv --> mol/m3/atm

      !---------------------------------------------------------------
      !psurf is in mb. multiply with 10.197e-4 to get atm
      !include molecular weights for air and CO2
       beta_gas=alpha_gas*(psurf*10.197e-4)*mair*1.e-3
     &                   /(tr_mm(itr)*1.e-3)
       !!atmCO2=368.6D0  !defined in obio_forc
!!!    beta_gas = beta_gas * tr_mm(itr)*1.e-3/rhows * atmCO2

cwatson        xco2 = atmCO2*1013.0/stdslp
cwatson       deltco2 = (xco2-pCO2_ij)*ff*1.0245E-3
cwatson       deltco2=atmCO2*1013.0/stdslp*ff*1.0245E-3  !beta_gas
cwatson               - pCO2_ij *ff*1.0245E-3       !trconstflx
cwatson ff is actually alpha_gas

      !trsf is really sfac = Kw_gas * beta_gas
      !units are such that flux comes out to (m/s)(kg/kg)
      !trsf = Kw_gas * beta_gas
      !trcnst = Kw_gas * beta_gas*trconstflx * byrho   ! convert to (conc * m/s)

       beta_gas = alpha_gas * 1.0245e-3 * psurf/1013.25  !stdslp and psurf in mb, no need to change units
       trsf = Kw_gas * beta_gas

       trcnst = Kw_gas * alpha_gas * 1.024e-3 * trconstflx * byrho   ! convert to (conc * m/s)

cdiag write(*,'(a,2i7,10e12.4)')'PBL, Kw ',
cdiag.  ilong,jlat,tg1,(Sc_gas/660.d0)**(-0.5d0),ws*ws,
cdiag.  Kw_gas,alpha_gas,beta_gas,trsf,trcnst,trconstflx,byrho

      RETURN
      END SUBROUTINE TRACERS_GASEXCH_CO2_Natassa_PBL

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c used with TRACERS_GASEXCH_CO2 to compute transfer velocity for CO2
c
c taken from Watsons' code
c ---------------------------------------------------------------------
      REAL*8 FUNCTION sol_co2(pt,ps)
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
      REAL*8 FUNCTION sc_co2(t)
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
