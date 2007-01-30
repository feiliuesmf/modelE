c used with TRACERS_GASEXCH_Natassa to compute transfer velocity for CFCs
c
c $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/Attic/TRACER_GASEXCH_Natassa.f,v $
c $Revision: 1.2 $
c $Date: 2007/01/30 16:28:36 $   ;  $State: Exp $
c $Author: ntausnev $ ;  $Locker:  $
c
c ---------------------------------------------------------------------
c $Log: TRACER_GASEXCH_Natassa.f,v $
c Revision 1.2  2007/01/30 16:28:36  ntausnev
c Correction bugs for TRACERS_GASEXCH_Natassa option. Functions must be real*8 type.
c
c Revision 1.1  2007/01/18 23:17:21  ntausnev
c Subroutines used with TRACERS_GASEXCH_Natassa to compute transfer velocity for CFCs
c
c Revision 1.2  1998/07/17 07:37:02  jomce
c Fixed slight bug in units conversion: converted 1.0*e-12 to 1.0e-12
c following warning from Matthew Hecht at NCAR.
c
c Revision 1.1  1998/07/07 15:22:00  orr
c Initial revision
c
c ---------------------------------------------------------------------
      REAL*8 FUNCTION sol_cfc(pt,ps,kn)  
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

c used with TRACERS_GASEXCH_Natassa to compute transfer velocity for CFCs
c
c  $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/Attic/TRACER_GASEXCH_Natassa.f,v $
c  $Revision: 1.2 $
c  $Date: 2007/01/30 16:28:36 $   ;  $State: Exp $
c  $Author: ntausnev $ ;  $Locker:  $
c 
c  ---------------------------------------------------------------------
c  $Log: TRACER_GASEXCH_Natassa.f,v $
c  Revision 1.2  2007/01/30 16:28:36  ntausnev
c  Correction bugs for TRACERS_GASEXCH_Natassa option. Functions must be real*8 type.
c
c  Revision 1.1  2007/01/18 23:17:21  ntausnev
c  Subroutines used with TRACERS_GASEXCH_Natassa to compute transfer velocity for CFCs
c
c  Revision 1.1  1998/07/07 15:22:00  orr
c  Initial revision
c 
c  ---------------------------------------------------------------------
      REAL*8 FUNCTION sc_cfc(t,kn)   
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
