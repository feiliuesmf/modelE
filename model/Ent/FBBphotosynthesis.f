!#define DEBUG 1
      module  photcondmod
      !This version of photcondmod does leaf level
      !photosynthesis (Farquhar and von Caemmerer, 1982) and
      !stomatal conductance (Ball and Berry, 1985, 1987).
      
      use FarquharBBpspar
      use ent_const

      implicit none

      save

      public init_ci, pscondleaf, biophysdrv_setup,calc_Pspar,ciMIN

      !=====CONSTANTS=====!
      real*8,parameter :: ciMIN = 1.d-8  !Small error
!      real*8,parameter :: Kelvin = 273.15d0
      real*8,parameter :: Rgas = gasc !8.314510d0 !gas constant (8.314510 J/mol K)
      real*8,parameter :: Kc = 30.0       !Michaelis-Menten coeff.constant for CO2 (Pa), Collatz (1991)
      real*8,parameter :: Ko = 3.d4       !Michaelis-Menten coeff.constant for O2 (Pa), Collatz (1991)
      real*8,parameter :: KcQ10 = 2.1d0   !Kc Q10 exponent, Collatz (1991)
      real*8,parameter :: KoQ10 = 1.2d0   !Ko Q10 exponent, Collatz (1991)

      !=====DECLARED TYPES======!
      type photosynthpar        !Calculated values from pft-dependent pspartypes
      integer :: pft            !Plant functional type.  1-C3 grassland
      real*8 :: PARabsorb       !Leaf PAR absorptance (fraction)
      real*8 :: Vcmax           !Maximum photosynthetic capacity (umol m-2 s-1)
      real*8 :: Kc              !Michaelis-Menten constant for CO2 (Pa)
      real*8 :: Ko              !Michaelis-Menten constant for O2 (Pa)
      real*8 :: Gammastar       !CO2 compensation point (Pa)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
      real*8 :: Nleaf           !g-N/m^2[leaf] - May want to take this from Tpool instead.
      end type photosynthpar

      private
      !=====GLOBAL VARIABLES (MODULE ONLY)=====!
      type(photosynthpar) :: pspar !Photosynthesis parameters.
!-----------------------------------------------------------------------------

      contains

      subroutine init_ci(ca, ci)
      implicit none
      real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)
      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction  (umol mol-1)

      !ci should be initialized. For rule of thumb, initialize to typically
      !observed ratio:
      ci = 0.7d0*ca

      end subroutine init_ci
      
!-----------------------------------------------------------------------------
      subroutine biophysdrv_setup(cf,ci,Tc,Pa,rh,psdrvpar)
      !* Set up met drivers for photosynthesis.
      implicit none
      real*8,intent(in) :: cf, ci, Tc, Pa, rh
      type(psdrvtype),intent(out) :: psdrvpar

      psdrvpar%cf = cf
      psdrvpar%ci = ci
      psdrvpar%Tc = Tc
      psdrvpar%Pa = Pa
      psdrvpar%rh = rh

      end subroutine biophysdrv_setup
!-----------------------------------------------------------------------------
      subroutine pscondleaf(pft,IPAR,psp,ca,Gb,gsout,Aout,Rdout
     &     ,sunlitshaded)
      implicit none
      integer,intent(in) :: pft
      real*8,intent(in) :: IPAR !umol m-2 s-1
      type(psdrvtype) :: psp
      real*8,intent(in) :: ca !Ambient air CO2 (umol mol-1)
      real*8,intent(in) :: Gb !mol m-2 s-1
      real*8,intent(out) :: gsout, Aout, Rdout !ci in psp
      integer,intent(in) :: sunlitshaded
      !---Local---
      real*8 :: ci, cs
      real*8,parameter :: LOW_LIGHT_LIMIT = 2.5d0 !umol m-2 s-1.  Nobel 1999, lower light limit for green plants is 0.7 W m-2 ~ 3 umol m-2 s-1.
      
!      call Collatz(pft,IPAR,psp%cf,psp%Tc,psp%Pa,psp%rh,psp%ci,
!     &     gsout,Aout)

!      if (IPAR.lt.LOW_LIGHT_LIMIT) then
!        Rdout = Respveg(pftpar(pft)%Nleaf,psp%Tc)  !Should be only leaf respiration!
!        Aout = 0.d0
!        cs = ca - (Aout-Rdout)*1.37d0/Gb
!        gsout = pftpar(pft)%b
!        psp%ci = ca             !Dummy assignment, no need to solve for ci 
!      else
        call Photosynth_analyticsoln(pft,IPAR,ca,ci,
     &     psp%Tc,psp%Pa,psp%rh,Gb,gsout,Aout,Rdout,sunlitshaded)
        psp%ci = ci             !Ball-Berry:  ci is analytically solved.  F-K: ci saved between time steps.
!      endif
        
      !Biological limits for gs - cuticular conductance?
      if(gsout.lt.(0.00006d0*psp%Pa/(gasc*(psp%Tc+KELVIN)))) then
        gsout=0.00006d0*psp%Pa/(gasc*(psp%Tc+KELVIN))
      endif

      end subroutine pscondleaf

!-----------------------------------------------------------------------------
      subroutine Photosynth_analyticsoln(pft,IPAR,ca,ci,Tl,Pa,rh,gb,
     o     gs,Atot,Rd,sunlitshaded)
      !@sum Photosynth_cubic Farquhar photosynthesis and Ball-Berry conductance
      !@sum and autotrophic respiration.  ci is solved for analytically for each
      !@sum of the limiting cases.
      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1) 
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)      
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(out) :: ci   !Leaf internal CO2 concentration (umol mol-1)      
      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (CO2 uptake, micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Dark = above-ground growth + maintenance respiration (umol m-2 s-1)
      integer,intent(in) :: sunlitshaded !For diagnostic outputs only.
        !---Local----
!      type(photosynthpar) :: pspar !Moved to global to module.
      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa) Not exactly .209*101325.
      real*8 :: Cip,Cie,Cic,Cis !Leaf internal CO2 (Pa)
      real*8 :: Je1, Jc1, Js1   !Assimilation of CO2, 3 limiting cases
      real*8 :: Anet            !Net assimilation of CO2 = Atot - aboveground respir (umol m-2 s-1)
      real*8 :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)

!      call calc_Pspar(pft, Pa, Tl, O2pres, pspar) !Moved up a module to reduce computation.
      Rd = Respveg(pspar%Nleaf,Tl)  !Should be only leaf respiration!
      
      call Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, Cie, Je1)
      call Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2pres, Cic, Jc1)
      call Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, Cis, Js1)

      Atot = min(Je1, Jc1, Js1)
      if (Atot.lt.0.d0) then 
#ifdef OFFLINE
        write(997,*) "Error, Atot<0.0:",Je1,Jc1,Js1,ca,gb,rh,IPAR,Pa,
     &       pspar,sunlitshaded
#endif
        Atot = 0.d0
!        return
      endif

      if (Atot.eq.Je1) then
        Cip = Cie
      else if (Atot.eq.Jc1) then
        Cip = Cic
      else if (Atot.eq.Js1) then
        Cip = Cis
      else !Atot was set to 0.d0 due to low light
        Cip = max(0.d0, Cie)
      endif

      Anet = Atot - Rd
      cs = ca - Anet*1.37d0/gb
      gs = BallBerry(Anet, rh, cs, pspar)
      ci = Cip/Pa*1d6

!#ifdef DEBUG
!      if (sunlitshaded.eq.1) then
!        write(991,*) IPAR, Cip/Pa*1d6, Cie/Pa*1d6,Cic/Pa*1d6,Cis/Pa*1d6,
!     &       Je1, Jc1, Js1, Atot, Rd, cs, gs, ca,gb,rh,IPAR,Pa, pspar
!      else
!        write(992,*) IPAR, Cip/Pa*1d6, Cie/Pa*1d6,Cic/Pa*1d6,Cis/Pa*1d6,
!     &       Je1, Jc1, Js1, Atot, Rd, cs, gs, pspar%Gammastar
!      endif
!#endif
      end subroutine Photosynth_analyticsoln
!-----------------------------------------------------------------------------

      subroutine Collatz(pft, IPAR,cs,Tl,Pa,rh,ci,gs,Anet)
!@sum Coupled photosynthesis/stomatal conductance at the leaf level
!@sum after Collatz, G.J., et.al. (1991) AgForMet 54:107-136
      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: IPAR            !Incident PAR (umol m-2 s-1) 
      !If PAR is not directly available, the following conversions may be
      !used:
      !  From total shortwave (W m-2) to PAR (umol m-2 s-1) (Monteith & Unsworth):
      !          PAR(umol m-2 s-1) = 2.3(umol/J)*SW(W m-2)
      !  From PAR (W m-2) to PAR (umol m-2 s-1) (U.Maryland, Dept. of Met., PAR Project),
      !  suggest nominal 485 nm for conversion, which gives:
      !          PAR(umol m-2 s-1) = 4.05(umol/J) * PAR(W m-2)

      !real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)      
      real*8,intent(in) :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: Pa   !Pressure (Pa)
      !real*8,intent(in) :: gb !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: ci              !Leaf internal CO2 mole fraction  (umol mol-1)
      real*8,intent(out) :: gs !Leaf stomatal conductance (mol-H2O m-2 s-1)
      real*8,intent(out) :: Anet !Leaf net photosynthesis (CO2 uptake, micromol m-2 s-1)
      !-----Local-----------------
      type(photosynthpar) :: pspar
!      type(metdatatype) :: mdat
      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa)
      real*8 :: Atot            !Gross assimilation (umol m-2 s-1)
!      real*8 :: Rd              !Dark or maintenance respiration (umol m-2 s-1)
      !ci should be intialized. For rule of thumb, initialize to
      ! initial ci = 0.7*ca

!      call load_metdata(IPAR,ca,Ta,rh,O2conc)

      call calc_Pspar(pft, Pa, Tl, O2pres)

      Atot = Farquhar(IPAR,ci*Pa*1.d-06,O2pres,Tl,pspar) 
      Anet = Atot - Respveg(pspar%Nleaf,Tl)

      gs = BallBerry(Atot, rh, cs, pspar)

!      if (if_ci.eq.1) ci = calc_ci(ca, gb, gs, Anet, IPAR, pspar)

      !Solving:
      !    Stomates adjust to control the gradient of ci to cs, so the coupled
      !photosynthesis/conductance model is generally solved by iterating
      !to obtain consistent Anet, gs, cs, and ci, and also rh at the leaf 
      !surface.  
      !    As this can be computationally intensive, an alternative is to
      !place bounds on how quickly gs can change (biologically measured).
      !Anet would not be affected, since it depends only on light, Vcmax, and
      !ci; however, ci would be affected.

      
      end subroutine Collatz

!-----------------------------------------------------------------------------
!      subroutine load_metdata(IPAR,ca,Ta,rh,O2concmol,mdat)
!      implicit none
!      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
!      real*8 :: ca   !Ambient air CO2 concentration (umol mol-1)
!      real*8 :: Ta   !Ambient air temperature (Celsius)
!      real*8 :: rh   !Relative humidity
!      real*8 :: O2concmol !O2 mole fraction in leaf (mmol mol-1)
!
!      mdat%IPAR = IPAR
!      mdat%ca = ca
!      mdat%Ta = Ta
!      mdat%rh = rh
!      mdat%O2concmol = O2concmol
!
!      end subroutine load_metdata
!-----------------------------------------------------------------------------

      function Respveg(Nleaf,Tl) Result(Rd)
      !@sum Respveg Autotrophic respiration (umol-CO2 m-2[leaf] s-1)
      !Rd = dark respiration = mitochondrial respiration =
      !  growth respiration(activity) + maintenance respiration (biomass)
      !Does not include photorespiration.
      !Need to distinguish aboveground respiration for leaf Ci vs. roots.

      implicit none
!      type(photosynthpar) :: pspar
      real*8,intent(in) :: Nleaf !(g-N/m^2 leaf) leaf nitrogen 
      real*8,intent(in) :: Tl !Leaf temperature (Celsius)
      real*8 :: Rd  !Autotrophic respiration (umol-CO2 m-2[leaf] s-1)
!      integer :: p
      !Collatz, et al. (1991).  No good.  Doesn't rise with temperature.
      !Rd = 0.015 * pspar%Vcmax 

      !* Friend and Kiang (2005) - total autotrophic respiration.
      !Rd  based on temperature and nitrogen content.
      !Rd = 0.2 (umol m-2 s-1) * N (g) (Carswell, et al., 2000)
      !N(g m-2) per Vcmax from Harley, et al. (1992, Fig. 4), cotton.
      !Temperature response from Bernacchi, et al. (2001)
!      Rd = 0.2d0 * (pspar%Vcmax + 9.6d0)/60.d0
!     &     *exp(18.72d0 - 46390.d0/(Rgas*(Tl+Kelvin)))
      !N(g m-2) per LAI from Ponca Ntot/LA, get mean 1st 120 days of season 2.47 g/m-2[leaf]
      !The Harley relation is an order of magnitude too small.
      Rd = Nleaf * exp(18.72d0 - 46390.d0/(Rgas*(Tl+Kelvin)))

!      Rd = exp(pftpar(p)%Rdc - pftpar(p)%RdH/(Rgas*(Tl+Kelvin))) !Harley&Tenhunen, 1991
      end function Respveg
!-----------------------------------------------------------------------------
      function calc_CO2compp(O2,Tl) Result(Gammastar)
!@sum CO2 compensation point in absence of dark respiration (Pa)

      implicit none
      real*8,intent(in) :: O2 !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: Tl !Leaf temperature (Celsius)
      real*8 :: Gammastar  !CO2 compensation point (Pa)
      !----Local-----
      real*8,parameter :: tau=2600.d0  !CO2/O2-specificity ratio

      Gammastar = O2/(2.d0*tau*Q10fn(0.57d0,Tl)) !Collatz (A3) 

      end function calc_CO2compp
!-----------------------------------------------------------------------------

      function Je(IPAR,Cip,pspar) Result(Je_light)
!@sum Photosynthetic rate limited by light electron transport (umol m-2 s-1)
      implicit none
      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: Cip  !Leaf internal CO2 partial pressure (Pa)
      type(photosynthpar) :: pspar
      real*8 :: Je_light !Electron-transport limited rate of photosynth (umol m-2 s-1)
      !----Local---------
      real*8,parameter :: alpha=.08d0 !Intrinsic quantum efficiency for CO2 uptake

      Je_light = (pspar%PARabsorb*IPAR)*alpha*(Cip-pspar%Gammastar)/
     &     (Cip+2*pspar%Gammastar)

      end function Je
!-----------------------------------------------------------------------------
      subroutine Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, Cip, Je_light)
      !@sum Ci_Je Analytical solution for Ci assuming Je is most limiting, 
      !@sum then calculation of Je.
      implicit none
      real*8,intent(in) :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8,intent(in) :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: rh              !Relative humidity
      real*8,intent(in) :: IPAR            !Incident PAR (umol m-2 s-1)
      real*8,intent(in) :: Pa              !
      type(photosynthpar) :: pspar
      real*8,intent(in) :: Rd              !Maintenance or mitochondrial respiration (umol-CO2 m-2[leaf] s-1)
      !---------------
      real*8,intent(out) :: Cip !Leaf internal CO2 partial pressure (Pa)
      real*8,intent(out) :: Je_light !Light-limited assimilation rate (umol m-2 s-1)
      !---Local------
      real*8,parameter :: alpha=.08d0 !Intrinsic quantum efficiency for CO2 uptake
      real*8 :: a1, e1, f1

      !Assimilation is of the form a1*(ci - Gammastar.umol)/(e1*ci + f1)
      a1 = pspar%PARabsorb*IPAR*alpha
      e1 = 1.d0
      f1 = 2*pspar%Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

      Cip = Pa * 1d-06 * ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar)
!      Cip = Pa * 1d-06 * 350.d0 * .7d0  !###Dummy check @350 ppm
      Je_light = Je(IPAR,Cip,pspar)

#ifdef DEBUG_ENT
      write(995,*) ca,rh,gb,IPAR, Pa,Rd,a1,e1,f1,pspar%m,pspar%b,
     &     pspar%Gammastar,Cip, Je_light
#endif
      end subroutine Ci_Je
!-----------------------------------------------------------------------------

      function Jc(Cip,O2,pspar) Result (Jc_RuBP)
!@sum Photosynthetic rate limited by RuBP saturation
      implicit none
      real*8 :: Cip  !Leaf internal CO2 partial pressure (Pa)
      real*8 :: O2  !O2 partial pressure in leaf (Pa)
      type(photosynthpar) :: pspar
      real*8 :: Jc_RuBP

      Jc_RuBP = pspar%Vcmax*(Cip - pspar%Gammastar)/
     &     (Cip + pspar%Kc*(1 + O2/Ko))
      end function Jc
!-----------------------------------------------------------------------------
      subroutine Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2, Cip, Jc_RuBP)
      !@sum Ci_Jc Analytical solution for Ci assuming Jc is most limiting.
      !@sum then calculation of Jc.
      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: rh              !Relative humidity
      real*8 :: IPAR            !Incident PAR (umol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      type(photosynthpar) :: pspar
      real*8 :: Rd              !Maintenance or mitochondrial respiration (umol-CO2 m-2[leaf] s-1)
      real*8 :: O2              !O2 partial pressure in leaf (Pa)
      !---------------
      real*8,intent(out) :: Cip !Leaf internal CO2 partial pressure (Pa)
      real*8,intent(out) :: Jc_RuBP !RuBP-limited assimilation rate (umol m-2 s-1)
      !---Local------
      real*8 :: a1, e1, f1

      !Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
      a1 = pspar%Vcmax
      e1 = 1.d0
      f1 = pspar%Kc*(1.d0 + O2/pspar%Ko) * 1.d06/Pa  !umol/mol

      Cip = Pa * 1d-06 * ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar)
      !Cip = Pa *1.D-06 * 350.d0 *.7d0 !Dummy prescribed ci.
      Jc_RuBP = Jc(Cip,O2,pspar)

#ifdef DEBUG_ENT
      write(993,*) ca,rh,gb,Pa,Rd,a1,e1,f1,pspar%m,pspar%b,
     &     pspar%Gammastar,Cip, Jc_RuBP
#endif      
      end subroutine Ci_Jc
!-----------------------------------------------------------------------------

      function Js(pspar) Result(Js_sucrose)
!@sum Photosynthetic rate limited by "utilization of photosynthetic products."
!@sum (umol m-2 s-1)
      implicit none
      type(photosynthpar) :: pspar
      real*8 :: Js_sucrose

      Js_sucrose = pspar%Vcmax/2.d0
      end function Js
!-----------------------------------------------------------------------------
      
      subroutine Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, Cip, Js_sucrose)
      !@sum Ci_Js Calculates Cip and Js, 
      !@sum Photosynthetic rate limited by "utilization of photosynthetic products."

      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: rh              !Relative humidity
      real*8 :: IPAR            !Incident PAR (umol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      type(photosynthpar) :: pspar
      real*8 :: Rd              !Maintenance or mitochondrial respiration (umol-CO2 m-2[leaf] s-1)
      real*8 :: Cip
      real*8 :: Js_sucrose !umol m-2 s-1
      !---Local----
      real*8 :: X,Y,Z,W,T       !Expressions from solving cubic equ. of ci.
      real*8 :: Anet            !umol m-2 s-1
      real*8 :: m               !Slope of Ball-Berry
      real*8 :: b               !Intercept of Ball-Berry
      real*8 :: rb              !Leaf boundary layer resistance = 1/gb

      m = pspar%m
      b = pspar%b
      rb = 1/gb

      X = (b*1.37d0*rb)**2.d0 - (m*rh - 1.65d0)*1.37d0*rb
      Y = ca*(m*rh - 1.65d0) - b*2.d0*1.37d0*rb*ca
      Z = b*1.37d0*rb - m*rh
      W = -b*ca
      T = b*(ca**2.d0)

      Js_sucrose = Js(pspar)
      Anet = Js_sucrose - Rd
      Cip = Pa * 1.d-06 * (-1.d0*((X*Anet + Y)*Anet + T )/(Z*Anet + W))

      end subroutine Ci_Js
!-----------------------------------------------------------------------------
      
      function arrhenius(Tcelsius,c1,c2) Result(arrh)
      !From David Medvigy's lphys.f90
      implicit none
      real*8 :: Tcelsius,c1,c2
      real*8 :: arrh

      arrh = c1*exp(c2*(1.d0/288.15d0-1.d0/(Tcelsius+Kelvin)))
      return
      end function arrhenius
!=================================================
      function Q10fn(Q10par,Tcelsius) Result(Q10factor)
      !@sum From Collatz, et al. (1991)
      implicit none
      real*8 :: Q10par, Tcelsius,Q10factor

      Q10factor = Q10par**((Tcelsius-25.d0)/10.d0)

      end function Q10fn
!=================================================

      function Tresponse(c,deltaH,Tcelsius) Result(Tfactor)
      !@sum From Bernacchi, et al. (2001).  Also Arrhenius.
      implicit none
      real*8,intent(in) :: c !Scaling factor
      real*8,intent(in) :: deltaH !Activation energy 
      real*8,intent(in) :: Tcelsius !Temperature (Celsius)
      real*8 :: Tfactor

      Tfactor = exp(c - deltaH/(Rgas * (Tcelsius + Kelvin)))

      end function Tresponse
!=================================================
      real*8 function ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar) Result(ci)
      !@sum ci_cubic Analytical solution for Ball-Berry/Farquhar cond/photosynth
      !@sum ci (umol/mol)
      !@sum For the case of assimilation being of the form:
      !@sum         A = a*(Cip - Gammastar)/(e*Cip + f) - Rd
      !@sum Numerator and denominator are converted from (Pa/Pa) to (umol mol-1)/(umol mol-1)
      !@sum         A = a1*(ci - gammamol) /(e1*ci + fmol) - Rd
      !@sum where gammamol = Gammastar*1d06/Pa, fmol = f1 = f*1d06/Pa

      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: rh              !Relative humidity
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      real*8 :: Rd              !Leaf mitochondrial respiration (umol m-2 s-1)
      real*8 :: a1              !Coefficient in linear Farquhar equ.
      real*8 :: e1              !Coefficient in linear Farquhar equ.
      real*8 :: f1              !Coefficient in linear Farquhar equ.
      type(photosynthpar) :: pspar
      !----Local----
      real*8 :: X, Y, Z, W, T   !Expressions from the quadratic of Assim.
      real*8 :: rb
      real*8 :: m
      real*8 :: b
      real*8 :: gammam
      real*8 :: c3, c2, c1, c   !Coefficients of the cubic of ci (c3*ci^3 + c2*ci^2 + c1*ci + c)

      real*8 :: c_1, c_2, c_3, e, fmol, gamol, a, C3_, A_,cs_,A_1,rs_
      real*8 :: res_, C6, Ra,A_2,cixx(3),K, ci_
      integer :: nroots, i
      rb = 1/gb
      m = pspar%m
      b = pspar%b
      gammam = pspar%Gammastar * 1.d06/Pa !Convert Pa to umol/mol

#define USE_IGORS_CUBIC
#ifndef USE_IGORS_CUBIC

      X = b*(1.37d0*rb)**2.d0 - (m*rh - 1.65d0)*1.37d0*rb 
      Y = ca*(m*rh - 1.65d0) - b*2.d0*1.37d0*rb*ca 
      Z = b*1.37d0*rb - m*rh
      W = -b*ca
      T = b*(ca**2.d0)
      
      c3 = e1*(Z*(a1 - Rd*e1) + W*e1)
      c2 = X *(a1-Rd*e1)**2.d0
     &     - Z*e1*( a1*gammam + Rd*f1) + (Y*e1 + Z*f1)* (a1 - Rd*e1)
     &     + W*2.d0*e1*f1 + T*e1**2.d0 
      c1 =  X *(-2.d0*(a1**2.d0)*gammam 
     &     - 2.d0*Rd*a1*(f1 - e1*gammam)
     &     + 2.d0*(Rd**2.d0)*e1*f1)
     &     + Y*f1*(a1 - Rd*e1) - (Y*e1 + Z*f1)*( a1*gammam + Rd) 
     &     + W*f1**2.d0 + T*2.d0*e1*f1 
      c  = X * ( (a1*gammam)**2.d0 + 2.d0*Rd*a1*f1*gammam 
     &     + (Rd*f1)**2.d0)
     &     - Y*f1*( a1*gammam + Rd*f1) + T*f1**2.d0

      !ci = cubicroot(c3, c2, c1, c)
      call cubicroot(c3, c2, c1, c, cixx, nroots)
      if ( nroots<1 ) call stop_model("ci_cubic: no solution",255)
      ci = maxval( cixx(1:3) )

!#DEBUG
      cs_ = ca - A_*C3_*rb
      write(992,*) "A,ca,cs,gs, pspar%pft, cixx,ci",
     &     a1*(ci - gamol)/(e1*ci + f1) - Rd, ca,cs_,
     &     m * (a1*(ci - gamol)/(e1*ci + f1) - Rd)*rh/cs_ + b,
     &     pspar%pft, cixx,ci

#else
      e = e1
      fmol = f1
      gamol = gammam
      a = a1
      C3_ =1.37d0
      C6 = 1.65d0
      Ra = C3_*rb

      c = -((ca*fmol + a*gamol*Ra + fmol*Ra*Rd)*
     -     (b*(ca*fmol + a*gamol*Ra + fmol*Ra*Rd) + 
     -       (a*gamol + fmol*Rd)*(C6 - m*rh)))

      c1 = 2*a**2*gamol*Ra*(C6 + b*Ra - m*rh) + 
     -   a*(2*b*ca*fmol*Ra - 2*b*ca*e*gamol*Ra + b*fmol*gamol*Ra + 
     -      2*b*fmol*Ra**2*Rd - 2*b*e*gamol*Ra**2*Rd + 
     -      C6*(fmol - e*gamol)*(ca + 2*Ra*Rd) - ca*fmol*m*rh + 
     -      ca*e*gamol*m*rh - fmol*gamol*m*rh - 2*fmol*m*Ra*Rd*rh + 
     -      2*e*gamol*m*Ra*Rd*rh) + 
     -   fmol*(-(b*(ca + Ra*Rd)*(2*ca*e - fmol + 2*e*Ra*Rd)) - 
     -      Rd*(2*C6*ca*e + 2*C6*e*Ra*Rd - 2*ca*e*m*rh + fmol*m*rh - 
     -         2*e*m*Ra*Rd*rh))

      c2 = -(a**2*Ra*(C6 + b*Ra - m*rh)) + 
     -   a*(2*b*ca*e*Ra - b*fmol*Ra + b*e*gamol*Ra + 2*b*e*Ra**2*Rd + 
     -      C6*e*(ca + 2*Ra*Rd) - ca*e*m*rh + fmol*m*rh - e*gamol*m*rh - 
     -      2*e*m*Ra*Rd*rh) - e*
     -    (b*(ca + Ra*Rd)*(ca*e - 2*fmol + e*Ra*Rd) + 
     -      Rd*(C6*ca*e + C6*e*Ra*Rd - ca*e*m*rh + 2*fmol*m*rh -
     &     e*m*Ra*Rd*rh)
     -      )

      c3 = e*(b*(ca*e - a*Ra + e*Ra*Rd) + m*(a - e*Rd)*rh)

      !ci = cubicroot(c3, c2, c1, c)
      call cubicroot(c3, c2, c1, c, cixx, nroots)
      if ( nroots<1 ) call stop_model("ci_cubic: no solution",255)
      ci = maxval( cixx(1:nroots) )

#endif

!#define DEBUG_CUBIC
#ifdef DEBUG_CUBIC
      e = e1
      fmol = f1
      gamol = gammam
      a = a1
      C3_ =1.37d0
      C6 = 1.65d0
      Ra = C3_*rb

      !!! uncomment these 3 lines to check all roots
      !c_1 = ci
      !do i=1,nroots
      !ci = cixx(i)
      print *,"-----------------"
      print *,"a,gammamol,fmol,Rd", a,gamol,fmol,Rd
      print *,"ci", ci

      A_ = a*(ci - gamol)/(e*ci + fmol) - Rd
      cs_ = ca - A_*C3_*rb
      rs_ = 1./(m*A_*rh/cs_ + b)
      A_1 = (cs_ - ci)/(C6*rs_)
      A_2 = (ca-cs_)/(C3_*rb)

      print *,"cs", cs_ !, ca-a*Ra, ca - A_*C3_*rb
      print *,"ca", ca !, C6*(ca-a*Ra)/(m*rh*a)
      print *,"rs", rs_ !, C6*(ca-a*Ra)/(m*rh*a)
      print *,"AAAA", A_
!      print *,"RES_A",     A_*C6/(ca-A_*Ra-ci)-m*A_*rh/(ca-A_*Ra)-b
!     &     , -A_*C6/(ca-A_*Ra-ci)-m*A_*rh/(ca-A_*Ra)-b

      ! residual for Igor's cubic
      res_ = -(a**2*(ci - gamol)**2*Ra*(C6 + b*Ra - m*rh)) - 
     -   a*(ci*e + fmol)*(ci - gamol)*
     -    (-2*b*ca*Ra + b*ci*Ra - 2*b*Ra**2*Rd - C6*(ca + 2*Ra*Rd) + 
     -      ca*m*rh - ci*m*rh + 2*m*Ra*Rd*rh) - 
     -   (ci*e + fmol)**2*(b*(ca + Ra*Rd)*(ca - ci + Ra*Rd) + 
     -      Rd*(C6*ca + C6*Ra*Rd - ca*m*rh + ci*m*rh - m*Ra*Rd*rh))
!      print *,"RES= ", res_, c + c1*ci + c2*ci**2 + c3*ci**3
!      print '(a,i3,i3,4e15.4)',"PLOT", i, pspar%pft, ci, A_, cs_, rs_
!      write(994,*) "A_1,A2,cs,rs,res_,i, pspar%pft, ci, A_, cs_, rs_",
!     &     A_1,A_2,cs_,rs_,res_,i, pspar%pft, ci, A_, cs_, rs_
      !!! uncomment these 2 lines to check all roots
      !enddo
      !ci = c_1
#endif

!#define USE_IGORS_CUBIC_2
#ifdef USE_IGORS_CUBIC_2
      b = b/C6
      K = m * rh / C6

      Y=fmol/e
      X= -a/e * (gamol+fmol/e)
      Z= a/e -Rd

      c = -(b*Ca*(X + (Ca + Y)*Z))
      c1 = Ca*Z - K*(X + Ca*Z + Y*Z) + 
     -   b*(Ca**2 + Ca*(Y + 2*Ra*Z) + Ra*(X + Y*Z))
      c2 = Ca*(-1 + K - 2*b*Ra) + K*(Y + Ra*Z) - Ra*(b*Y + Z + b*Ra*Z)
      c3 = Ra*(1 - K + b*Ra)

      call cubicroot(c3, c2, c1, c, cixx, nroots)
      print *,"NNNN ",cixx(1:nroots)
      !print *,"UUUU", A_1-cixx(2)

      call Baldocchi(c3,c2,c1,c)

       !if ( a/e -Rd < 0 ) then
       !  write(502,'(20f15.6)') a1,e1,gamol,f1,Rd,K,b,Ca,Ra
       !endif


#endif

#ifdef DEBUG
      write(994,*) 'rb,m,b,ca,gammam,X,Y,Z,W,T,a1,e1,f1,
     & c3,c2,c1,c,ci,expr(ci)',
     &     rb,m,b,ca,gammam,X,Y,Z,W,T,a1,e1,f1,c3,c2,c1,c,ci,
     &     c3*(ci**3.d0) + c2*(ci**2.d0) + c1*ci + c
#endif
      ! testing new solution for "ci"
      !ci_ = ci_cubic1(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar)
      end function ci_cubic

      real*8 function ci_cubic1(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar)Result(ci)
      !@sum ci_cubic Analytical solution for Ball-Berry/Farquhar cond/photosynth
      !@sum ci (umol/mol)
      !@sum For the case of assimilation being of the form:
      !@sum         A = a*(Cip - Gammastar)/(e*Cip + f) - Rd
      !@sum Numerator and denominator are converted from (Pa/Pa) to (umol mol-1)/(umol mol-1)
      !@sum         A = a1*(ci - gammamol) /(e1*ci + fmol) - Rd
      !@sum where gammamol = Gammastar*1d06/Pa, fmol = f1 = f*1d06/Pa

      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: rh              !Relative humidity
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      real*8 :: Rd              !Leaf mitochondrial respiration (umol m-2 s-1)
      real*8 :: a1              !Coefficient in linear Farquhar equ.
      real*8 :: e1              !Coefficient in linear Farquhar equ.
      real*8 :: f1              !Coefficient in linear Farquhar equ.
      type(photosynthpar) :: pspar
      !----Local----
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 :: Ra, b, K, gamol, A_d_asymp
      real*8 :: X, Y, Z, Y1  ! tmp vars
      real*8 :: c3, c2, c1, c   !Coefficients of the cubic of ci (c3*ci^3 + c2*ci^2 + c1*ci + c)
      real*8 :: cixx(3), A ! solutions of cubic
      real*8 :: cs, Rs ! needed to compute ci
      integer :: nroots, i

      Ra = 1/gb * S_ATM
      b = pspar%b / S_STOM
      K = pspar%m * rh / S_STOM
      gamol = pspar%Gammastar * 1.d06/Pa !Convert Pa to umol/mol
      A_d_asymp = - b*Ca / (K - b*Ra) ! asymptotic val of A from diffusion eq.

      ! first check some special cases
      if ( A_d_asymp >= 0.d0 ) then
        ! this can happen only for very low humidity
        ! probably should never happen in the real world, but if it does,
        ! this case should be considered separately
        print *,"K<b*Ra: m,rh,b,Ra:",pspar%m,rh,b,Ra
        call stop_model("ci_cubic: rh too small ?",255)
      endif

      Y= f1/e1
      X= -a1/e1 * (gamol+f1/e1)
      Z= a1/e1 -Rd

      if ( Z > 0.d0 ) then
        ! Farquhar curve is above zero. May have solution A > 0
        c = -(b*Ca*(X + (Ca + Y)*Z))
        c1 = Ca*Z - K*(X + Ca*Z + Y*Z) + 
     &       b*(Ca**2 + Ca*(Y + 2*Ra*Z) + Ra*(X + Y*Z))
        c2 = Ca*(-1 + K - 2*b*Ra) + K*(Y + Ra*Z) - Ra*(b*Y + Z + b*Ra*Z)
        c3 = Ra*(1 - K + b*Ra)

        call cubicroot(c3, c2, c1, c, cixx, nroots)

        ! find minimal root above the asymptotic value
        A = 1.d30
        do i=1,nroots
          if ( cixx(i) < A .and. cixx(i) > A_d_asymp ) A = cixx(i)
        enddo
        if ( A == 1.d30 ) call stop_model("ci_cubic: no solution",255)

        if ( A >= 0 ) then
          cs = ca - A*Ra
          Rs = 1.d0 / (K*A/cs + b)
          ci = cs - A*Rs
          ! just in case, check consistency
          if ( ci < 0.d0 ) call stop_model("ci_cubic: ci<0",255)
          if ( cs < 0.d0 ) call stop_model("ci_cubic: cs<0",255)
          print *,'QQQQ ',A
          return
        endif

      endif

      ! if we got here then A<0 : have to solve quaratic equation

      Y1 = Y + ca
      c2 = Ra + 1.d0/b
      c1 = - (Y1 + c2*Z)
      c  = X + Y1*Z

      ! just in case,
      if (  c1*c1 - 4.d0*c2*c < 0.d0 )
     &     call stop_model("ci_cubic: no solution to quadratic",255)
      A = ( - c1 - sqrt( c1*c1 - 4.d0*c2*c ) ) / ( 2.d0 * c2 )
      cs = ca - A*Ra
      Rs = 1.d0 / ( b)
      ci = cs - A*Rs
      print *,"q ", ci, cs, A, Rs, Ra
      ! just in case, check consistency
      if ( ci < 0.d0 ) call stop_model("ci_cubic: q: ci<0",255)
      if ( cs < 0.d0 ) call stop_model("ci_cubic: q: cs<0",255)
      print *,'QQQQ ',A

      end function ci_cubic1


      subroutine Baldocchi(c3,c2,c1,c)
      real*8, intent(in) :: c3,c2,c1,c
      !---
      real*8 p,q,r, Q_,R_,thet, x1,x2,x3

      p = c2/c3
      q = c1/c3
      r = c/c3

      Q_ = (p**2 - 3.d0*q)/9.d0
      R_ = ( 2*p**3 - 9*p*q + 27*r )/54.d0

      if ( Q_ < 0.d0 ) then
        print *,"BBBB"," OOPS1"
        return
      endif

      if ( R_**2 > Q_**3 ) then
        print *,"BBBB"," OOPS2"
        return
      endif

      thet = acos(R_/sqrt(Q_**3))

      x1 = -2.d0 * sqrt(Q_) * cos((thet     )/3.d0) - p/3.d0
      x2 = -2.d0 * sqrt(Q_) * cos((thet+2*Pi)/3.d0) - p/3.d0
      x3 = -2.d0 * sqrt(Q_) * cos((thet+4*Pi)/3.d0) - p/3.d0

      print *,"BBBB", x1,x2,x3
      
cddd
cddd      p = (e*bet + fmol*thet - a*alf + e*alf*Rd)/(e*alf)
cddd      q = (e*gam + fmol*gam/Ca - a*bet + a*gamol*thet
cddd     &     + e*Rd*bet + Rd*fmol*thet)
cddd     &     /(e*alf)
cddd      r = (-a*gam + a*gamol*gam/Ca + e*Rd*gam + Rd*fmol*gam/Ca)/(e*alf)
cddd
cddd
cddd
cddd      p = (e*bet + b*thet - a*alf + e*alf*Rd)/(e*alf)
cddd      q = (e*gam + b*gam/Ca - alf*bet + a*d*thet + e*Rd*bet + Rd*b*thet)
cddd     &     /(e*alf)
cddd      r = (-a*gam + a*d*gam
cddd
      end subroutine Baldocchi

!=================================================
      !real*8 function cubicroot(a,b,c,d) Result(x)
      subroutine cubicroot(a,b,c,d,x,n) 
      !* solve cubic equation: a x^3 + b x^2 + c x + d = 0 *!
      !* Written by Igor Aleinov from solution by Cardano in
      !* Korn, Korn, Mathematical Handbook.
      implicit none
      real*8,intent(in) :: a,b,c,d  ! coefficients of cubic
      real*8, intent(out) :: x(:)   ! results ( 0-3 roots )
      integer, intent(out) :: n     ! number of roots
      real*8 :: x0,x1,x2
      real*8 :: a0,a1,a2,Q1,R1,D1
      real*8, parameter :: EPS0 = 1.d-15
      real*8, parameter :: one3rd = 1.d0/3.d0
      real*8 :: arg, S, T
      complex*16 :: ST

      !print *,"cubicroot:",a,b,c,d

      if (abs(a) < (abs(b)+abs(c)+abs(d))*EPS0 ) then
        if (abs(b) < (abs(c)+abs(d))*EPS0) then
#ifdef DEBUG_RZTOOLS
          write(*,*) "Cardano: linear eq: 1 solution"
          write(*,*) "Cardano:   c= %g  d= %g\n", c, d
#endif
          if (abs(c) < abs(d)*EPS0) then
            write(*,*) "Internal Error in Cardano: no solution."
            stop
          endif
          x0 = -d/c
          x(1) = x0
          !write(*,*) "Cardano: returning",x0
          n = 1
        else
          write(*,*) "What's this?"
          D1 = c*c - 4.d0*b*d
          
          if (D1 > 0.d0) then
            Q1 = sqrt(D1)
            x0 = (-c + Q1) / (2.d0 * b)
            x1 = (-c - Q1) / (2.d0 * b)
            !return
            n = 2
          else if (D1.eq.0.) then
            x0 = -c / (2.d0 * b)
            x1 = x0
            n = 1
          else 
            x0 = -c /(2.d0 *b)
            x1 = sqrt(-D1) / (2.d0* b)
            n = 0
          end if
        end if
        !print *,"CX1",x0,x1
        !x = max(x0,x1)
        x(1) = x0
        x(2) = x1
      else
        a2 = b/a
        a1 = c/a
        a0 = d/a
        Q1 = (3.d0 * a1 - a2*a2 ) / 9.d0
        R1 = (9.d0 * a2 * a1 - 27.d0 * a0 - 2.d0 * a2*a2*a2) /54.d0
        D1 = Q1*Q1*Q1 + R1*R1
        !write(*,*) "abcda2a1a0Q1R1D1",a,b,c,d,a2,a1,a0,Q1,R1,D1
        if (D1 > 0.d0) then       !* only one real root *!
          !write(*,*) "One real root."
          arg = R1 + sqrt(D1)
          S = sign(1.d0, arg) * (abs(arg)**one3rd)
          arg = R1 - sqrt(D1)
          T = sign(1.d0, arg) * (abs(arg)**one3rd)
          x0 = -a2/3.d0 + S + T
          x1 = -a2/3.d0 - (S+T)*0.5d0
          x2 = sqrt(3.d0) * (S-T)*0.5d0
          !print *,"CX2",x0,x1,x2
          n = 1
        else if (D1.eq.0.) then !* two roots coincide * *!
          !write(*,*) "Two roots coincide."
          S = sign(1.d0, R1) * (abs(R1)**one3rd)
          x0 = -a2/3.d0 + 2.d0*S
          x1 = -a2/3.d0 - S
          x2 = x1
          !print *,"CX3",x0,x1,x2
          n =2
        else                    !* three different real roots *!
          !call CRtCube( R1, sqrt(-D1), S, T)
          !write(*,*) "Three different real roots. a2R1D1ST",a2,R1,D1,S,T
          ST = ( cmplx(R1, sqrt(-D1),kind(1.d0)) )**one3rd
          S = real (ST)
          T = aimag(ST)
          x0 = -a2/3.d0 + 2.d0*S
          x1 = -a2/3.d0 - S + sqrt(3.d0)*T
          x2 = -a2/3.d0 - S - sqrt(3.d0)*T
          !print *,"CX4",x0,x1,x2
          n = 3
        end if
        !x = max(x0,x1,x2)
        !x = x2
        x(1) = x0
        x(2) = x1
        x(3) = x2
      end if
      end subroutine cubicroot
      
!=================================================

      subroutine calc_Pspar(pft, Pa, Tl, O2pres)
      !@sum calc_Pspar Collatz photosynthesis parameters. GLOBAL TO MODULE.
      !@sum Later need to replace these with von Caemmerer book Arrhenius
      !@sum function sensitivies (her Table 2.3)
      implicit none
      integer,intent(in) :: pft  !Plant functional type, 1=C3 grassland
      real*8 :: Pa  !Atmospheric pressure (Pa)
      real*8 :: Tl  !Leaf temperature (Celsius)
      real*8 :: O2pres !O2 partial pressure in leaf (Pa)
!      type(photosynthpar),intent(inout) :: pspar !Moved to global to module.
      integer :: p
      !Below parameters are declared at top of module, though only used here.
!      real*8,parameter :: Kc              !Michaelis-Menten constant for CO2 (Pa)
!      real*8,parameter :: Ko              !Michaelis-Menten constant for O2 (Pa)
!      real*8,parameter :: KcQ10           !Kc Q10 exponent
!      real*8,parameter :: KoQ10           !Ko Q10 exponent

      p = pft
      pspar%pft = pft
      pspar%PARabsorb = pftpar(p)%PARabsorb !Collatz et al. (1991)
      !pspar%Vcmax = pftpar(p)%Vcmax/(1 + exp((-220.e03+703.*(Tl+Kelvin))
!     &     /(Rgas*(Tl+Kelvin))))
      pspar%Vcmax = pftpar(p)%Vcmax * Q10fn(2.21d0, Tl)

      pspar%Kc = Kc*Q10fn(KcQ10,Tl) !(Collatz, eq. A12)
      pspar%Ko = Ko*Q10fn(KoQ10,Tl) !(Collatz, eq. A12)
      pspar%Gammastar = calc_CO2compp(O2pres,Tl) !(Pa) (Collatz)
      pspar%m = pftpar(p)%m     !Slope of Ball-Berry equation (Collatz)
      pspar%b = pftpar(p)%b     !Intercept of Ball-Berry equation (mol m-2 s-1) (Collatz)
      pspar%Nleaf = pftpar(p)%Nleaf !g-N/m^2[leaf] Needed for foliar respiration.
      !pspar%Nleaf = pftpar(p)%Nleaf * phenology factor !Here can adjust Nleaf according
                                !to foliage N pools or phenology factor.
!      write(*,*) "Tl, pspar:",Tl, pspar

      end subroutine calc_Pspar

!-----------------------------------------------------------------------------
      
      function Farquhar(IPAR,Cip,O2,Tl,pspar) Result(Assim)
!@sum Photosynthesis at the leaf level (umol m-2 s-1)
!@sum after Farquhar, et al.(1980) Planta, 149:78-90.
      implicit none
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
      real*8,intent(in) :: Cip  !Leaf internal CO2 partial pressure (Pa)
      real*8,intent(in) :: O2   !Leaf internal O2 partial pressure (Pa)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      type(photosynthpar) :: pspar
      real*8 :: Assim !Leaf GPP (positive into leaf) (umol m-2 s-1)
      !----Local-----
      real*8 :: Je1, Jc1, Js1
      !Assume pspar already calculated.
      !call calc_Pspar(pft, Tl, pspar)

      Je1 = Je(IPAR,Cip,pspar)
      Jc1 = Jc(Cip,O2,pspar)
      Js1 = Js(pspar)
      Assim = min(Je1, Jc1, Js1)
!#ifdef DEBUG
!      write(997,*) IPAR, Cip, O2,pspar%Gammastar,pspar%Vcmax,
!     &     Tl,Je1,Jc1,Js1,Assim
!#endif
      !Assim = min(Je(IPAR,Cip,pspar), Jc(Cip,O2,pspar), Js(pspar))

      end function Farquhar

!-----------------------------------------------------------------------------


      function Farquhar_standalone(pft,IPAR,ca,Tl,rh,Pa,gb,ci,gs)
     o     Result(Anet)
!@sum Photosynthesis at the leaf level.
!@sum after Farquhar, et al.(1980) Planta, 149:78-90.
!@sum Standalone version.  Can be called directly from external program.

      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1) 
      real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction  (umol mol-1)
      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
      real*8 :: Anet            !Leaf net photosynthesis (CO2 uptake, micromol m-2 s-1)
      !-----Local-----------------
      type(photosynthpar) :: pspar
!      type(metdatatype) :: mdat
      real*8,parameter :: O2conc=209 !O2 mole fraction in leaf (mmol mol-1)
      real*8 :: cs     !CO2 concentration at leaf surface (umol mol-1)

      !ci should be intialized. For rule of thumb, initialize to
      ! initial ci = 0.7*ca

      cs = ca      !Assign CO2 concentration at leaf surface
                   !More rigorously, cs should also be solved for.

      call calc_Pspar(pft, Pa, Tl, O2conc*Pa*1.d-06)

      Anet = Farquhar(IPAR,ci*Pa*1.d-06,O2conc*Pa*1.d-06,Tl,pspar) 
     &     - Respveg(pspar%Nleaf,Tl)

      end function Farquhar_standalone
!-----------------------------------------------------------------------------

      function BallBerry(Anet, rh, cs, pspar) Result (gsw)
!@sum Ball-Berry (1987) model of leaf stomatal conductance of 
!@sum water vapor, gsw (mol m-2 s-1)      
      implicit none
      real*8,intent(in) :: Anet !Net assimilation of CO2 (umol m-2 s-1)
      real*8,intent(in) :: rh   !Relative humidity (fractional ratio)
      real*8,intent(in) :: cs   !Leaf surface CO2 mole fraction (umol mol-1)
      type(photosynthpar) :: pspar
      real*8 :: gsw !Leaf conductance of water vapor (mol m-2 s-1)
      !----Local-----
      
      gsw = pspar%m*Anet*rh/cs + pspar%b
      if (gsw < 0.0) gsw = 0.d0

      end function BallBerry

!-----------------------------------------------------------------------------
      
      function calc_ci(ca,gb, gs,Anet,IPAR,pspar) Result(ci)
!@sum Leaf internal CO2 conc (umol mol-1) assuming diffusive flux of CO2
!@sum is at steady-state with biochemical uptake by photosynthesis and
!@sum that there is zero leaf boundary layer resistance (infinite gb),
!@sum and that there is no leaf cuticular conductance of CO2.
!@sum Full equation:  ci = ca - Anet*(1.6/gb + 1.4/gs)
!@sum 1.6 = ratio of diffusivities of CO2 and water vapor in laminar flow
!@sum       in the leaf boundary layer
!@sum 1.4 = ratio of diffusivities of CO2 and water vapor in still air at
!@sum       the leaf surface
!@sum (Monteith, 1995;  Kiang, 2002)
      implicit none

      real*8 :: ca !Ambient air CO2 mole fraction at surface reference height (umol mol-1)
      real*8 :: gb !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: gs !Stomatal conductance of water vapor(mol m-2 s-1)
      real*8 :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
      type(photosynthpar) :: pspar
      real*8 :: ci !Leaf internal CO2 mole fraction (umol mol-1)
      !----Local------
      real*8,parameter :: MINPARMOL=50.d0  !(umol m-2 s-1)

      if (IPAR.lt.MINPARMOL) then  !Stomates closed
        ci = ca - Anet*1.6d0/pspar%b
      else
        ci = ca - Anet*(1.6d0/gb + 1.4d0/gs)
      endif

      if (ci.lt.ciMIN) ci = ciMIN  !Keep positive definite.

      end function calc_ci

!-----------------------------------------------------------------------------

      end module photcondmod

