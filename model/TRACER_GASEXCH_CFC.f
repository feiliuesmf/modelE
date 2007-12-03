      MODULE TRACER_GASEXCH_COM

      USE TRACER_COM, only : ntm    !tracers in air-sea gas exch

      implicit none

      SAVE



#include "dimensions.h"
#include "dimension2.h"


      SAVE


      real*8 atracflx(iia,jja,ntm),atrac(iia,jja,ntm)
      common /gasexch1/atracflx,atrac

      real   tracflx(idm,jdm,ntm)     !  tracer flux at air-sea intfc
      common /gasexch2/tracflx

      END MODULE TRACER_GASEXCH_COM

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------

      SUBROUTINE TRACERS_GASEXCH_CFC_Natassa_PBL(tg1,ws,
     . alati,psurf,itr,trconstflx,byrho,Kw_gas,alpha_gas,
     . beta_gas,trsf,trcnst,ilong,jlat)

      USE CONSTANT, only:    rhows,mair
      USE TRACER_COM, only : ntm,trname,tr_mm
      
      implicit none

      integer :: ilong,jlat,itr
      real*8  :: tg1,ws,psurf,trconstflx,byrho,trsf,trcnst
      real*8  :: alati,Kw_gas,alpha_gas,beta_gas
      real*8  :: Sc_gas
      real*8, external :: sc_cfc,sol_cfc
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !OCMIP implementation www.ipsl.jussieu.fr/OCMIP
      !F=Kw*Csat - Kw*Csurf=
      !  Kw*alpha*trs - Kw*trs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !new treatment for units
      !F=Kw*Csat -Kw*Csurf
      ! =Kw*alpha*mol_weight_air/mol_weight_cfc11*surfp*Cair  !!!TERM_1*Cair
      ! -Kw*rho_water/mol_weight_cfc11*Csurf                  !!!TERM_2
      !
      ! where, Kw                in m/s
      !        alpha                mol/m^3/atm
      !        mol_weight_air       Kg_air
      !        mol_weight_cfc11     Kg_CFC-11
      !        surfp                atm
      !        Cair                 Kg_CFC-11/Kg_air
      !        rho_water            Kg_water/m^3
      !        Csurf                Kg_CFC-11/Kg_water
      !then F is in  (mol/m^3)(m/s)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---------------------------------------------------------------
      !TRANSFER VELOCITY
      !---------------------------------------------------------------
      ! compute Schmidt number for gas
      ! use ground temperature in deg C including skin effects
      Sc_gas=sc_cfc(tg1,11)

      !wind speed ws: magn. of surf. wind modified by buoyancy flux (m/s)
      !compute transfer velocity Kw
      !only over ocean

      if (Sc_gas .le. 0.) then
        write(*,'(a,2i4,a,2f9.3)')
     .          'warning: Sc_gas negtv, at ',ilong,jlat,
     .          ', Sc_gas,temp_c=',Sc_gas,tg1
         Kw_gas=1.e-10
      else
         Kw_gas=
     &       1.d0/3.6e+5*0.337d0*ws*ws*(Sc_gas/660.d0)**(-0.5d0)
      endif


      !---------------------------------------------------------------
      !gas SOLUBILITY
      !---------------------------------------------------------------
      !alpha --solubility of CFC (11 or 12) in seawater
      !in mol/m^3/picoatm
       alpha_gas=sol_cfc(tg1,sss_loc,11)
      !convert to mol/m^3/atm
       alpha_gas=alpha_gas*1.e+12

      !---------------------------------------------------------------
      !psurf is in mb. multiply with 10.197e-4 to get atm
      !include molecular weights for air and CFC-11
       beta_gas=alpha_gas*(psurf*10.197e-4)*mair*1.e-3
     &                   /(tr_mm(itr)*1.e-3)
!!!    beta_gas = beta_gas * tr_mm(itr)*1.e-3/rhows

      !trsf is really sfac = Kw_gas * beta_gas
      !units are such that flux comes out to (m/s)(kg/kg)
       trsf = Kw_gas * beta_gas

       trcnst = Kw_gas * trconstflx(itr)*byrho ! convert to (conc * m/s)

      RETURN
      END SUBROUTINE TRACERS_GASEXCH_CFC_Natassa_PBL

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------

c used with TRACERS_GASEXCH_CFC to compute transfer velocity for CFCs
c
c $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/TRACER_GASEXCH_CFC.f,v $
c $Revision: 2.3 $
c $Date: 2007/12/03 22:40:36 $   ;  $State: Exp $
c $Author: aromanou $ ;  $Locker:  $
c
c ---------------------------------------------------------------------
c $Log: TRACER_GASEXCH_CFC.f,v $
c Revision 2.3  2007/12/03 22:40:36  aromanou
c 1) Ocean_hycom.f : Full gas exchage between CO2 in the ocean and CO2 in the atmosphere. Ocean
c    arays are mapped on the atmospheric grid before the gas exchange computation.
c 2) PBL.f : Separate routines for PBL code for CO2 and CFC exchange: TRACER_GASEXCH_CO2.f,
c TRACER_GASEXCH_CFC.f.
c 3) atmospheric CO2 is kept constant here (very special case: check PBL.f)
c 4) Implementation for using Ron Miller's dust fluxes or GOCART dust fluxes
c 5) Recalcuate phytoplankton sinking rates in obio_update.f
c 6) relevant rundeck is E3arobio7.R
c
c ****my last commit should be the same as this one, but had no commentary, so I redid it!
c
c Revision 2.2  2007/12/03 22:23:20  aromanou
c *** empty log message ***
c
c Revision 2.1  2007/05/01 22:23:22  aromanou
c Ocean Biogeochemistry based on Watson Gregg's model is now implemented.
c (no whales, shrimp or human swimmers, yet though).
c
c Rundeck to use E3arobio.R
c
c Gas exchange between ocean and atmosphere are also implemented.
c Case study1: CO2 exchange
c Case study2: CFC exchange
c
c Notes: 1)ocean hycom routines are NOT "best run" routines.
c        2)modelE is current CVS, better outcome with frozen version.
c
c Revision 1.2  1998/07/17 07:37:02  jomce
c Fixed slight bug in units conversion: converted 1.0*e-12 to 1.0e-12
c following warning from Matthew Hecht at NCAR.
c
c Revision 1.1  1998/07/07 15:22:00  orr
c Initial revision
c
c ---------------------------------------------------------------------
      REAL FUNCTION sol_cfc(pt,ps,kn)
c-------------------------------------------------------------------
c
c     CFC 11 and 12 Solubilities in seawater
c     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
c
c     pt:       temperature (degre Celcius)
c     ps:       salinity    (o/oo)
c     kn:       11 = CFC-11, 12 = CFC-12
c     sol_cfc:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm

c
c     J-C Dutay - LSCE
c-------------------------------------------------------------------

      REAL*8    pt, ps,ta,d
      REAL*8 a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 b1 ( 11: 12), b2 ( 11: 12), b3 ( 11: 12)

      INTEGER kn

cc
cc coefficient for solubility in  mol/l/atm
cc ----------------------------------------
c
c     for CFC 11
c     ----------
      a1 ( 11) = -229.9261d0
      a2 ( 11) =  319.6552d0
      a3 ( 11) =  119.4471d0
      a4 ( 11) =   -1.39165d0
      b1 ( 11) =   -0.142382d0
      b2 ( 11) =    0.091459d0
      b3 ( 11) =   -0.0157274d0
c
c     for CFC/12
c     ----------
      a1 ( 12) = -218.0971d0
      a2 ( 12) =  298.9702d0
      a3 ( 12) =  113.8049d0
      a4 ( 12) =   -1.39165d0
      b1 ( 12) =   -0.143566d0
      b2 ( 12) =    0.091015d0
      b3 ( 12) =   -0.0153924d0
c
      ta       = ( pt + 273.16d0)* 0.01d0
      d    = ( b3 ( kn)* ta + b2 ( kn))* ta + b1 ( kn)

c
c
      sol_cfc
     $    = exp ( a1 ( kn)
     $    +       a2 ( kn)/ ta
     $    +       a3 ( kn)* dlog ( ta )
     $    +       a4 ( kn)* ta * ta  + ps* d )
c
c     conversion from mol/(l * atm) to mol/(m^3 * atm)
c     ------------------------------------------------
      sol_cfc = 1000. * sol_cfc
c
c     conversion from mol/(m^3 * atm) to mol/(m3 * pptv)
c     --------------------------------------------------
      sol_cfc = 1.0e-12 * sol_cfc

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c used with TRACERS_GASEXCH_CFC to compute transfer velocity for CFCs
c
c  $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/TRACER_GASEXCH_CFC.f,v $
c  $Revision: 2.3 $
c  $Date: 2007/12/03 22:40:36 $   ;  $State: Exp $
c  $Author: aromanou $ ;  $Locker:  $
c 
c  ---------------------------------------------------------------------
c  $Log: TRACER_GASEXCH_CFC.f,v $
c  Revision 2.3  2007/12/03 22:40:36  aromanou
c  1) Ocean_hycom.f : Full gas exchage between CO2 in the ocean and CO2 in the atmosphere. Ocean
c     arays are mapped on the atmospheric grid before the gas exchange computation.
c  2) PBL.f : Separate routines for PBL code for CO2 and CFC exchange: TRACER_GASEXCH_CO2.f,
c  TRACER_GASEXCH_CFC.f.
c  3) atmospheric CO2 is kept constant here (very special case: check PBL.f)
c  4) Implementation for using Ron Miller's dust fluxes or GOCART dust fluxes
c  5) Recalcuate phytoplankton sinking rates in obio_update.f
c  6) relevant rundeck is E3arobio7.R
c
c  ****my last commit should be the same as this one, but had no commentary, so I redid it!
c
c  Revision 2.2  2007/12/03 22:23:20  aromanou
c  *** empty log message ***
c
c  Revision 2.1  2007/05/01 22:23:22  aromanou
c  Ocean Biogeochemistry based on Watson Gregg's model is now implemented.
c  (no whales, shrimp or human swimmers, yet though).
c
c  Rundeck to use E3arobio.R
c
c  Gas exchange between ocean and atmosphere are also implemented.
c  Case study1: CO2 exchange
c  Case study2: CFC exchange
c
c  Notes: 1)ocean hycom routines are NOT "best run" routines.
c         2)modelE is current CVS, better outcome with frozen version.
c
c  Revision 1.1  1998/07/07 15:22:00  orr
c  Initial revision
c 
c  ---------------------------------------------------------------------
      REAL FUNCTION sc_cfc(t,kn)
c---------------------------------------------------
c     CFC 11 and 12 Schmidt number 
c     as a function of temperature. 
c
c     ref: Zheng et al (1998), JGR, vol 103,No C1 
c
c     t: temperature (degree Celcius)
c     kn: = 11 for CFC-11,  12 for CFC-12
c
c     J-C Dutay - LSCE
c---------------------------------------------------
      IMPLICIT NONE 
      INTEGER kn
      REAL*8  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 t
c
c   coefficients with t in degre Celcius
c   ------------------------------------
      a1(11) = 3501.8d0
      a2(11) = -210.31d0
      a3(11) =    6.1851d0
      a4(11) =   -0.07513d0
c
      a1(12) = 3845.4d0
      a2(12) = -228.95d0
      a3(12) =    6.1908d0
      a4(12) =   -0.067430d0
c

      sc_cfc = a1(kn) + a2(kn) * t + a3(kn) *t*t  
     &         + a4(kn) *t*t*t
  
      RETURN 
      END 
