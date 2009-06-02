#include "rundeck_opts.h"
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
     &     ,frost_hardiness

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
      real*8 :: Kc              !Michaelis-Menten parameter for CO2 (Pa)
      real*8 :: Ko              !Michaelis-Menten parameter for O2 (Pa)
      real*8 :: Gammastar       !CO2 compensation point (Pa)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
      real*8 :: Nleaf           !g-N/m^2[leaf] - May want to take this from Tpool instead.
      real*8 :: stressH2O       !Water stress factor (fraction, 1=no stress)
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
      subroutine biophysdrv_setup(ca,ci,Tc,Pa,rh,psdrvpar)
      !* Set up met drivers for photosynthesis.
      implicit none
      real*8,intent(in) :: ca, ci, Tc, Pa, rh
      type(psdrvtype),intent(out) :: psdrvpar

      psdrvpar%ca = ca
      psdrvpar%ci = ci
      psdrvpar%Tc = Tc
      psdrvpar%Pa = Pa
      psdrvpar%rh = rh

      end subroutine biophysdrv_setup
!-----------------------------------------------------------------------------
      subroutine pscondleaf(pft,IPAR,psd,Gb,gsout,Aout,Rdout
     &     ,sunlitshaded,ISPout)
      implicit none
      integer,intent(in) :: pft
      real*8,intent(in) :: IPAR !umol m-2 s-1. Absorbed PAR. Should APAR.
      type(psdrvtype) :: psd
      real*8,intent(in) :: Gb !mol m-2 s-1
      real*8,intent(out) :: gsout, Aout, Rdout !ci in psd
      real*8,intent(out) :: ISPout
      integer,intent(in) :: sunlitshaded
      !---Local---
      real*8 :: ci!, cs
      real*8,parameter :: LOW_LIGHT_LIMIT = 2.5d0 !umol m-2 s-1.  Nobel 1999, lower light limit for green plants is 0.7 W m-2 ~ 3 umol m-2 s-1.
      
!      if (IPAR.lt.LOW_LIGHT_LIMIT) then
!        Rdout = Respveg(pftpar(pft)%Nleaf,psd%Tc)  !Should be only leaf respiration!
!        Aout = 0.d0
!        cs = ca - (Aout-Rdout)*1.37d0/Gb
!        gsout = pftpar(pft)%b
!        psd%ci = ca             !Dummy assignment, no need to solve for ci 
!      else
cddd      print *,"called Photosynth_analyticsoln",
cddd     &     pft,IPAR,psd%ca,ci,
cddd     &     psd%Tc,psd%Pa,psd%rh,Gb,gsout,Aout,Rdout,sunlitshaded

        call Photosynth_analyticsoln(pft,IPAR,psd%ca,ci,
     &     psd%Tc,psd%Pa,psd%rh,Gb,gsout,Aout,Rdout,sunlitshaded,
     &  ISPout)
        psd%ci = ci             !Ball-Berry:  ci is analytically solved.  F-K: ci saved between time steps.

!      endif
        
      !Biological limits for gs - cuticular conductance?
      if(gsout.lt.(0.00006d0*psd%Pa/(gasc*(psd%Tc+KELVIN)))) then
        gsout=0.00006d0*psd%Pa/(gasc*(psd%Tc+KELVIN))
      endif

      end subroutine pscondleaf

!-----------------------------------------------------------------------------

      subroutine Photosynth_analyticsoln(pft,IPAR,ca,ci,Tl,Pa,rh,gb,
     o     gs,Atot,Rd,sunlitshaded,isp)
      !@sum Photosynth_cubic Farquhar photosynthesis and Ball-Berry conductance
      !@sum and autotrophic respiration.  ci is solved for analytically for each
      !@sum of the limiting cases.
      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: IPAR !Absorbed PAR.  WRONG OLD COMMENT:Incident PAR (umol m-2 s-1) 
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)      
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(out) :: ci   !Leaf internal CO2 concentration (umol mol-1)      
      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (CO2 uptake, micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Dark = above-ground growth + maintenance respiration (umol m-2 s-1)
      real*8,intent(out) :: isp ! Isoprene emission (umol C m-2 s-1)
      integer,intent(in) :: sunlitshaded !For diagnostic outputs only.
        !---Local----
!      type(photosynthpar) :: pspar !Moved to global to module.
      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa) Not exactly .209*101325.
      real*8 :: cie, cic, cis   !Leaf internal CO2 (umol mol-1)
      real*8 :: Je1, Jc1, Js1   !Assimilation of CO2, 3 limiting cases
      real*8 :: Anet            !Net assimilation of CO2 = Atot - aboveground respir (umol m-2 s-1)
      real*8 :: Aiso            ! Rate of photosynthesis for isoprene emissions (umol m-2 s-1)
      real*8 :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)
      real*8 :: Ae, Ac, As
      real*8 :: a1,f1,e1
      real*8, parameter :: alpha=.08d0 !Intrinsic quantum efficiency for CO2 uptake

      integer, save :: counter = 0
      counter = counter + 1

      !write(888,*) "counter=", counter

!      Rd = Respveg(pspar%Nleaf,Tl)  !Old F&K Respveg is not only leaf respiration.
      Rd = 0.015d0 * pspar%Vcmax    !von Caemmerer book.

!      call Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, cie, Je1)
      ! Photosynthetic rate limited by light electron transport (umol m-2 s-1)
      ! Je_light = (pspar%PARabsorb*IPAR)*alpha*(Cip-pspar%Gammastar)/
      !            (Cip+2*pspar%Gammastar)

      !Assimilation is of the form a1*(ci - Gammastar.umol)/(e1*ci + f1)

!      a1 = pspar%PARabsorb*IPAR*alpha
      a1 = IPAR*alpha  !### HACK:  IPAR from canopyspitters.f is APAR.  When we switch to Wenze's canopyrad, then leaf PARabsorb will be used -NK ###

      e1 = 1.d0
      f1 = 2*pspar%Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

      call ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,Ae)
      !write(888,*) "Ae", ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,Ae

!      call Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2pres, cic, Jc1)
      ! Photosynthetic rate limited by RuBP saturation
      ! Jc_RuBP = pspar%Vcmax*(Cip - pspar%Gammastar)/
      !           (Cip + pspar%Kc*(1 + O2/pspar%Ko))

      !Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
      a1 = pspar%Vcmax
      e1 = 1.d0
      f1 = pspar%Kc*(1.d0 + O2pres/pspar%Ko) * 1.d06/Pa  !umol/mol

      call ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,Ac)
      !write(888,*) "Ac", ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,Ac

!      call Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, cis, Js1)
      !Photosynthetic rate limited by "utilization of photosynthetic products"
      ! (umol m-2 s-1)
      !Js_sucrose = pspar%Vcmax/2.d0
      As = pspar%Vcmax/2.d0 - Rd
      !write(888,*) "As", As

      Anet = min(Ae, Ac, As)
      Atot = Anet + Rd
      Aiso = Ae + Rd

      if (Atot.lt.0.d0) then
        ! can only happen if ca < Gammastar . Does it make sense? -Yes-NK
#ifdef OFFLINE
        write(997,*) "Error, Atot<0.0:",Ae,Ac,As,ca,gb,rh,IPAR,Pa,
     &       pspar,sunlitshaded
#endif
        Atot = 0.d0
        Anet = - Rd
        ci = pspar%Gammastar * 1.d06/Pa  
        return
      endif

      cs = ca - Anet*1.37d0/gb
      gs = BallBerry(Anet, rh, cs, pspar)
      ci = cs - Anet/(gs/1.65d0)

#ifdef PS_BVOC
         call Voccalc(pft,pa,ca,ci,Tl,pspar%Gammastar,
     & isp,Aiso)

#else
       isp=0.0d0
#endif
      !write(888,*) "gs,ci,cs", gs,ci,cs

      end subroutine Photosynth_analyticsoln

!-----------------------------------------------------------------------------
      subroutine Voccalc(pft,pa,ca,ci,Tl,Gammastar,isp,Aiso)
!@sum Isoprene emissions coupled to photosynthesis

      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)      
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: Gammastar   
      real*8,intent(in) :: Aiso   !(umol m-2 s-1)
      real*8,intent(in) :: ci
      real*8,intent(out) :: isp ! isoprene emission (umol C m-2 s-1)
      type(photosynthpar) :: pspar
        !---Local----
      real*8 :: gammamol,fact
      real*8 :: IBASER, Y_alpha, Y_eps, kapco2
      real*8 :: tauiso
      gammamol =  Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

C Y_alpha, Y_eps unitless
  
      Y_alpha=(ci-gammamol)/(6.0*(4.67*ci+9.33*gammamol))

         if(pft.eq.1)then 
         Y_eps = 0.0d0
         endif
         if(pft.eq.2)then 
         Y_eps = 2.10d-02
         endif
         if(pft.eq.3)then 
         Y_eps = 8.24d-02
         endif
         if(pft.eq.4)then 
         Y_eps = 6.48d-02
         endif
         if(pft.eq.5)then 
         Y_eps = 1.08d-01
         endif
         if(pft.eq.6)then 
         Y_eps = 4.44d-02
         endif
         if(pft.eq.7)then 
         Y_eps = 1.38d-01
         endif
         if(pft.eq.8)then 
         Y_eps = 0.0d0
         endif       

       isp = Y_eps*Aiso*Y_alpha

C Include CO2 effects

       kapco2 = (0.7*370.0)/(0.7*ca)

C Include temperature effects

       tauiso = exp(0.1*(Tl-30.0))

C Include seasonal effects? Add later.
C Note can switch on and off kapco2

       isp = isp*kapco2*tauiso

      end subroutine Voccalc

!-----------------------------------------------------------------------------
cddd      subroutine Photosynth_analyticsoln1(pft,IPAR,ca,ci,Tl,Pa,rh,gb,
cddd     o     gs,Atot,Rd,sunlitshaded)
cddd      !@sum Photosynth_cubic Farquhar photosynthesis and Ball-Berry conductance
cddd      !@sum and autotrophic respiration.  ci is solved for analytically for each
cddd      !@sum of the limiting cases.
cddd      implicit none
cddd      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
cddd      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1) 
cddd      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)      
cddd      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
cddd      real*8,intent(in) :: rh   !Relative humidity
cddd      real*8,intent(in) :: Pa   !Pressure (Pa)
cddd      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8,intent(out) :: ci   !Leaf internal CO2 concentration (umol mol-1)      
cddd      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
cddd      real*8,intent(out) :: Atot !Leaf gross photosynthesis (CO2 uptake, micromol m-2 s-1)
cddd      real*8,intent(out) :: Rd  !Dark = above-ground growth + maintenance respiration (umol m-2 s-1)
cddd      integer,intent(in) :: sunlitshaded !For diagnostic outputs only.
cddd        !---Local----
cddd!      type(photosynthpar) :: pspar !Moved to global to module.
cddd      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa) Not exactly .209*101325.
cddd      real*8 :: cie, cic, cis   !Leaf internal CO2 (umol mol-1)
cddd      real*8 :: Je1, Jc1, Js1   !Assimilation of CO2, 3 limiting cases
cddd      real*8 :: Anet            !Net assimilation of CO2 = Atot - aboveground respir (umol m-2 s-1)
cddd      real*8 :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)
cddd
cddd!      call calc_Pspar(pft, Pa, Tl, O2pres, pspar) !Moved up a module to reduce computation.
cddd      Rd = Respveg(pspar%Nleaf,Tl)  !Should be only leaf respiration!
cddd
cddd      call Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, cie, Je1)
cddd      call Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2pres, cic, Jc1)
cddd      call Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, cis, Js1)
cddd
cddd      Atot = min(Je1, Jc1, Js1)
cddd      if (Atot.lt.0.d0) then 
cddd#ifdef OFFLINE
cddd        write(997,*) "Error, Atot<0.0:",Je1,Jc1,Js1,ca,gb,rh,IPAR,Pa,
cddd     &       pspar,sunlitshaded
cddd#endif
cddd        Atot = 0.d0
cddd!        return
cddd      endif
cddd
cddd      if (Atot.eq.Je1) then
cddd        ci = cie
cddd        write(830,*) 1
cddd      else if (Atot.eq.Jc1) then
cddd        ci = cic
cddd        write(830,*) 2
cddd      else if (Atot.eq.Js1) then
cddd        ci = cis
cddd        write(830,*) 3
cddd      else !Atot was set to 0.d0 due to low light
cddd        ci = max(0.d0, cie)
cddd        write(830,*) 4
cddd      endif
cddd
cddd      Anet = Atot - Rd
cddd      cs = ca - Anet*1.37d0/gb
cddd      gs = BallBerry(Anet, rh, cs, pspar)
cddd      !ci = Cip/Pa*1d6
cddd
cddd      ! the following solution seems to be more straightforward
cddd      !cs = ca - Anet/(gb/1.37d0)
cddd      !gs = m*Anet*rh/cs + b
cddd      write(833,*) ci, cs - Anet/(gs/1.65d0),Anet
cddd      ci = cs - Anet/(gs/1.65d0)
cddd
cddd
cddd!#ifdef DEBUG
cddd!      if (sunlitshaded.eq.1) then
cddd!        write(991,*) IPAR, Cip/Pa*1d6, Cie/Pa*1d6,Cic/Pa*1d6,Cis/Pa*1d6,
cddd!     &       Je1, Jc1, Js1, Atot, Rd, cs, gs, ca,gb,rh,IPAR,Pa, pspar
cddd!      else
cddd!        write(992,*) IPAR, Cip/Pa*1d6, Cie/Pa*1d6,Cic/Pa*1d6,Cis/Pa*1d6,
cddd!     &       Je1, Jc1, Js1, Atot, Rd, cs, gs, pspar%Gammastar
cddd!      endif
cddd!#endif
cddd      end subroutine Photosynth_analyticsoln1
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
!       Rd = 0.015 * pspar%Vcmax !Only leaf maintenance respiration.

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
      function calc_CO2compp(O2,Kc,Ko,Tl) Result(Gammastar)
!@sum CO2 compensation point in absence of dark respiration (Pa)

      implicit none
      real*8,intent(in) :: O2 !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: Kc   !Michaelis-Menten parameter for CO2 (Pa)
      real*8,intent(in) :: Ko   !Michaelis-Menten parameter for O2 (Pa)
      real*8,intent(in) :: Tl !Leaf temperature (Celsius)
      real*8 :: Gammastar  !CO2 compensation point (Pa)
      !----Local-----
      real*8,parameter :: tau=2600.d0  !CO2/O2-specificity ratio

!      Gammastar = O2/(2.d0*tau*Q10fn(0.57d0,Tl)) !Collatz (A3)
!      Gammastar = O2*Q10fn(1.75,Tl)/(2.d0*tau) !Collatz (A3) Same as above, !KcQ10/KoQ10 = 2.1/1.2 = 1.75 = 1/.57 
      Gammastar = 0.5d0*(Kc/Ko)*0.21*O2 !CLM Tech Note. Gives smaller Gammastar than Collatz.


      end function calc_CO2compp
!-----------------------------------------------------------------------------

cddd      subroutine Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, ci, Je_light)
cddd      !@sum Ci_Je Analytical solution for Ci assuming Je is most limiting, 
cddd      !@sum then calculation of Je.
cddd      implicit none
cddd      real*8,intent(in) :: ca              !Ambient air CO2 concentration (umol mol-1)
cddd      real*8,intent(in) :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8,intent(in) :: rh              !Relative humidity
cddd      real*8,intent(in) :: IPAR            !Incident PAR (umol m-2 s-1)
cddd      real*8,intent(in) :: Pa              !
cddd      type(photosynthpar) :: pspar
cddd      real*8,intent(in) :: Rd              !Maintenance or mitochondrial respiration (umol-CO2 m-2[leaf] s-1)
cddd      !---------------
cddd      real*8,intent(out) :: ci !Leaf internal CO2 concentration (umol mol-1)
cddd      real*8,intent(out) :: Je_light !Light-limited assimilation rate (umol m-2 s-1)
cddd      !---Local------
cddd      real*8,parameter :: alpha=.08d0 !Intrinsic quantum efficiency for CO2 uptake
cddd      real*8 :: a1, e1, f1, A
cddd
cddd      ! Photosynthetic rate limited by light electron transport (umol m-2 s-1)
cddd      ! Je_light = (pspar%PARabsorb*IPAR)*alpha*(Cip-pspar%Gammastar)/
cddd      !            (Cip+2*pspar%Gammastar)
cddd
cddd      !Assimilation is of the form a1*(ci - Gammastar.umol)/(e1*ci + f1)
cddd      a1 = pspar%PARabsorb*IPAR*alpha
cddd      e1 = 1.d0
cddd      f1 = 2*pspar%Gammastar * 1.d06/Pa !Convert from Pa to umol/mol
cddd
cddd      call ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,ci,A)
cddd!      Cip = Pa * 1d-06 * 350.d0 * .7d0  !###Dummy check @350 ppm
cddd      Je_light = A + Rd
cddd
cddd#ifdef DEBUG_ENT
cddd      write(996,*) ca,rh,gb,IPAR, Pa,Rd,a1,e1,f1,pspar%m,pspar%b,
cddd     &     pspar%Gammastar,Cip, Je_light
cddd#endif
cddd      end subroutine Ci_Je
cddd!-----------------------------------------------------------------------------
cddd
cddd      subroutine Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2, ci, Jc_RuBP)
cddd      !@sum Ci_Jc Analytical solution for Ci assuming Jc is most limiting.
cddd      !@sum then calculation of Jc.
cddd      implicit none
cddd      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
cddd      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8 :: rh              !Relative humidity
cddd      real*8 :: IPAR            !Incident PAR (umol m-2 s-1)
cddd      real*8 :: Pa              !Pressure (Pa)
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: Rd              !Maintenance or mitochondrial respiration (umol-CO2 m-2[leaf] s-1)
cddd      real*8 :: O2              !O2 partial pressure in leaf (Pa)
cddd      !---------------
cddd      real*8,intent(out) :: ci !Leaf internal CO2 concentration (umol mol-1)
cddd      real*8,intent(out) :: Jc_RuBP !RuBP-limited assimilation rate (umol m-2 s-1)
cddd      !---Local------
cddd      real*8 :: a1, e1, f1, A
cddd
cddd      ! Photosynthetic rate limited by RuBP saturation
cddd      ! Jc_RuBP = pspar%Vcmax*(Cip - pspar%Gammastar)/
cddd      !           (Cip + pspar%Kc*(1 + O2/pspar%Ko))
cddd
cddd      !Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
cddd      a1 = pspar%Vcmax
cddd      e1 = 1.d0
cddd      f1 = pspar%Kc*(1.d0 + O2/pspar%Ko) * 1.d06/Pa  !umol/mol
cddd
cddd      call ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,ci, A)
cddd      !Cip = Pa *1.D-06 * 350.d0 *.7d0 !Dummy prescribed ci.
cddd      Jc_RuBP = A + Rd
cddd
cddd#ifdef DEBUG_ENT
cddd      write(993,*) ca,rh,gb,Pa,Rd,a1,e1,f1,pspar%m,pspar%b,
cddd     &     pspar%Gammastar,Cip, Jc_RuBP
cddd#endif      
cddd      end subroutine Ci_Jc
cddd!-----------------------------------------------------------------------------
cddd
cddd      subroutine Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, ci, Js_sucrose)
cddd      !@sum Ci_Js Calculates Cip and Js, 
cddd      !@sum Photosynthetic rate limited by "utilization of photosynthetic products."
cddd
cddd      implicit none
cddd      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
cddd      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8 :: rh              !Relative humidity
cddd      real*8 :: IPAR            !Incident PAR (umol m-2 s-1)
cddd      real*8 :: Pa              !Pressure (Pa)
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: Rd              !Maintenance or mitochondrial respiration (umol-CO2 m-2[leaf] s-1)
cddd      real*8 :: ci              !Leaf internal CO2 concentration (umol mol-1)
cddd      real*8 :: Js_sucrose !umol m-2 s-1
cddd      !---Local----
cddd      real*8 :: X,Y,Z,W,T       !Expressions from solving cubic equ. of ci.
cddd      real*8 :: Anet            !umol m-2 s-1
cddd      real*8 :: m               !Slope of Ball-Berry
cddd      real*8 :: b               !Intercept of Ball-Berry
cddd      real*8 :: rb              !Leaf boundary layer resistance = 1/gb
cddd      real*8 cs,gs
cddd
cddd      m = pspar%m
cddd      b = pspar%b
cddd      rb = 1/gb
cddd
cddd      X = (b*1.37d0*rb)**2 - (m*rh - 1.65d0)*1.37d0*rb
cddd      Y = ca*(m*rh - 1.65d0) - b*2.d0*1.37d0*rb*ca
cddd      Z = b*1.37d0*rb - m*rh
cddd      W = -b*ca
cddd      T = b*(ca**2)
cddd
cddd      !Js_sucrose = Js(pspar)
cddd      !Photosynthetic rate limited by "utilization of photosynthetic products"
cddd      ! (umol m-2 s-1)
cddd      Js_sucrose = pspar%Vcmax/2.d0
cddd      Anet = Js_sucrose - Rd
cddd      ci = -1.d0*((X*Anet + Y)*Anet + T )/(Z*Anet + W)
cddd
cddd      ! the following solution seems to be more straightforward
cddd      cs = ca - Anet/(gb/1.37d0)
cddd      gs = m*Anet*rh/cs + b
cddd      write(834,*) ci, cs - Anet/(gs/1.65d0),cs
cddd      !ci = cs - Anet/(gs/1.65d0)
cddd
cddd      end subroutine Ci_Js
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

      subroutine ci_cubic(ca,rh,gb,Pa,Rd,a1,e1,f1,pspar,A)
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
      real*8, intent(out) :: A
      type(photosynthpar) :: pspar
      !----Local----
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 :: Ra, b, K, gamol, A_d_asymp
      real*8 :: X, Y, Z, Y1  ! tmp vars
      real*8 :: c3, c2, c1, c   !Coefficients of the cubic of ci (c3*ci^3 + c2*ci^2 + c1*ci + c)
      real*8 :: cixx(3) ! solutions of cubic
!      real*8 :: cs, Rs ! needed to compute ci
      integer :: nroots, i
!      real*8 ci

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
        !!print *,"!!! A_d_asymp >= 0.d0 !!!", A_d_asymp
        A_d_asymp = -1.d30 !!! hack
        !!print *,"K<b*Ra: m,rh,b,Ra:",pspar%m,rh,b,Ra
        !!call stop_model("ci_cubic: rh too small ?",255)
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

        !!print *,"roots= ", cixx(1:nroots)

        ! find minimal root above the asymptotic value
        A = 1.d30
        do i=1,nroots
          if ( cixx(i) < A .and. cixx(i) > A_d_asymp ) A = cixx(i)
        enddo
        if ( A == 1.d30 )  then
          print *," m,rh,b,Ra:",pspar%m,rh,b,Ra
          print *,"ca,gb,Pa:",ca,gb,Pa
          print *,"pspar:",pspar
          print *," A_d_asymp,K,gamol,f1,a1,e1,Rd",
     &         A_d_asymp,K,gamol,f1,a1,e1,Rd
          print *,"c3,c2,c1,c", c3,c2,c1,c
          print *,"nroots,cixx",nroots,cixx(1:nroots)
          call stop_model("ci_cubic: no solution",255)
        endif

        if ( A >= 0 ) then
cddd          cs = ca - A*Ra
cddd          Rs = 1.d0 / (K*A/cs + b)
cddd          ci = cs - A*Rs
cddd          ! just in case, check consistency
cddd          if ( ci < 0.d0 ) call stop_model("ci_cubic: ci<0",255)
cddd          if ( cs < 0.d0 ) call stop_model("ci_cubic: cs<0",255)
cddd          !!print *,'QQQQ ',A,ci
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
cddd      cs = ca - A*Ra
cddd      Rs = 1.d0 / ( b)
cddd      ci = cs - A*Rs
cddd      !!print *,"q ", ci, cs, A, Rs, Ra
cddd      ! just in case, check consistency
cddd      if ( ci < 0.d0 ) call stop_model("ci_cubic: q: ci<0",255)
cddd      if ( cs < 0.d0 ) call stop_model("ci_cubic: q: cs<0",255)
cddd      !!print *,'QQQQ ',A,ci

      end subroutine ci_cubic


!=================================================
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
      real*8, parameter :: EPS0 = 1.d-8 ! 1.d-15
      real*8, parameter :: one3rd = 1.d0/3.d0
      real*8 :: arg, S, T
      complex*16 :: ST

      !print *,"cubicroot:",a,b,c,d

      if (abs(a) < (abs(b)+abs(c)+abs(d))*EPS0 ) then
        if (abs(b) < (abs(c)+abs(d))*EPS0) then
          if (abs(c) < abs(d)*EPS0) then
            write(*,*) "Internal Error in Cardano: no solution."
            stop
          endif
          x0 = -d/c
          x(1) = x0
          !write(*,*) "Cardano: returning",x0
          n = 1
        else
          !write(*,*) "What's this?"
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

      subroutine calc_Pspar(dtsec,pft,Pa,Tl,O2pres,stressH2O,
     &                      Sacclim,llspan)
      !@sum calc_Pspar Collatz photosynthesis parameters in type pspar, which
      !@sum is GLOBAL TO MODULE.
      !@sum Later need to replace these with von Caemmerer book Arrhenius
      !@sum function sensitivities (her Table 2.3)
!      use phenology, only : frost_hardiness ! REPEAT: dependency issues
      implicit none
      integer,intent(in) :: pft   !Plant functional type, 1=C3 grassland
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: Pa     !Atmospheric pressure (Pa)
      real*8,intent(in) :: Tl     !Leaf temperature (Celsius)
      real*8,intent(in) :: O2pres !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: stressH2O
      real*8,intent(in) :: Sacclim !state of acclimation/frost hardiness
      real*8,intent(in) :: llspan !mean leaf life span
!      type(photosynthpar),intent(inout) :: pspar !Moved to global to module.
      integer :: p

      !----Local-----
      real*8 :: facclim ! acclimation/forst hardiness factor [-]
      !Below parameters are declared at top of module, though only used here.
!      real*8,parameter :: Kc              !Michaelis-Menten constant for CO2 (Pa)
!      real*8,parameter :: Ko              !Michaelis-Menten constant for O2 (Pa)
!      real*8,parameter :: KcQ10           !Kc Q10 exponent
!      real*8,parameter :: KoQ10           !Ko Q10 exponent
      real*8 :: fparlimit !light(i.e.,PAR) control
      integer, save :: counter = 0
      counter = counter + 1

      facclim = frost_hardiness(Sacclim)
     
!      fparlimit = par_phenology(pft,llspan)
      fparlimit = 1.d0

!!! this var is not reproducible on restart, please figure out why
!      fparlimit = 1.d0 ! seems to be ok now

      !write(877,*) "counter", counter
      !write(877,*) "facclim", facclim
      !write(877,*) "fparlimit", fparlimit

      p = pft
      pspar%pft = pft
      pspar%PARabsorb = pftpar(p)%PARabsorb !Collatz et al. (1991)
!      pspar%Vcmax = pftpar(p)%Vcmax/(1 + exp((-220.e03+703.*(Tl+Kelvin))
!     &     /(Rgas*(Tl+Kelvin))))
      pspar%Vcmax = pftpar(p)%Vcmax * Q10fn(2.21d0, Tl)
     &            * facclim * fparlimit
      pspar%Kc = Kc*Q10fn(KcQ10,Tl) !(Collatz, eq. A12)
      pspar%Ko = Ko*Q10fn(KoQ10,Tl) !(Collatz, eq. A12)
      pspar%Gammastar = calc_CO2compp(O2pres,pspar%Kc,pspar%Ko,Tl) !(Pa) (Collatz)
      pspar%m = stressH2O*pftpar(p)%m     !Slope of Ball-Berry equation (Collatz)
!      pspar%m = pftpar(p)%m     !Slope of Ball-Berry equation (Collatz)
      pspar%b = pftpar(p)%b     !Intercept of Ball-Berry equation (mol m-2 s-1) (Collatz)
      pspar%Nleaf = pftpar(p)%Nleaf !g-N/m^2[leaf] Needed for foliar respiration.
      !pspar%Nleaf = pftpar(p)%Nleaf * phenology factor !Here can adjust Nleaf according
                                !to foliage N pools or phenology factor.
      pspar%stressH2O = stressH2O

      end subroutine calc_Pspar

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
      
      ! just in case check cs (remove after debugging ?)
      if ( cs <= 0.d0 ) call stop_model("BallBerry: cs <= 0", 255)
      gsw = pspar%m*Anet*rh/cs + pspar%b
      if (gsw < pspar%b) gsw = pspar%b

      end function BallBerry

!-----------------------------------------------------------------------------
cddd      function Je(IPAR,Cip,pspar) Result(Je_light)
cddd!@sum Photosynthetic rate limited by light electron transport (umol m-2 s-1)
cddd      implicit none
cddd      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
cddd      real*8 :: Cip  !Leaf internal CO2 partial pressure (Pa)
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: Je_light !Electron-transport limited rate of photosynth (umol m-2 s-1)
cddd      !----Local---------
cddd      real*8,parameter :: alpha=.08d0 !Intrinsic quantum efficiency for CO2 uptake
cddd
cddd      Je_light = (pspar%PARabsorb*IPAR)*alpha*(Cip-pspar%Gammastar)/
cddd     &     (Cip+2*pspar%Gammastar)
cddd
cddd      end function Je
!-----------------------------------------------------------------------------
cddd      function Jc(Cip,O2,pspar) Result (Jc_RuBP)
cddd!@sum Photosynthetic rate limited by RuBP saturation
cddd      implicit none
cddd      real*8 :: Cip  !Leaf internal CO2 partial pressure (Pa)
cddd      real*8 :: O2  !O2 partial pressure in leaf (Pa)
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: Jc_RuBP
cddd
cddd      Jc_RuBP = pspar%Vcmax*(Cip - pspar%Gammastar)/
cddd     &     (Cip + pspar%Kc*(1 + O2/pspar%Ko))
cddd!!!old     &     (Cip + pspar%Kc*(1 + O2/Ko))
cddd      end function Jc
!-----------------------------------------------------------------------------
cddd      function Js(pspar) Result(Js_sucrose)
cddd!@sum Photosynthetic rate limited by "utilization of photosynthetic products."
cddd!@sum (umol m-2 s-1)
cddd      implicit none
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: Js_sucrose
cddd
cddd      Js_sucrose = pspar%Vcmax/2.d0
cddd      end function Js
!-----------------------------------------------------------------------------
cddd      subroutine Collatz(dtsec,pft,IPAR,cs,Tl,Pa,rh,ci,gs,Anet,Sacclim)
cddd!@sum Coupled photosynthesis/stomatal conductance at the leaf level
cddd!@sum after Collatz, G.J., et.al. (1991) AgForMet 54:107-136
cddd      implicit none
cddd      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
cddd      real*8,intent(in) :: IPAR            !Incident PAR (umol m-2 s-1) 
cddd      real*8,intent(in) :: dtsec
cddd      !If PAR is not directly available, the following conversions may be
cddd      !used:
cddd      !  From total shortwave (W m-2) to PAR (umol m-2 s-1) (Monteith & Unsworth):
cddd      !          PAR(umol m-2 s-1) = 2.3(umol/J)*SW(W m-2)
cddd      !  From PAR (W m-2) to PAR (umol m-2 s-1) (U.Maryland, Dept. of Met., PAR Project),
cddd      !  suggest nominal 485 nm for conversion, which gives:
cddd      !          PAR(umol m-2 s-1) = 4.05(umol/J) * PAR(W m-2)
cddd
cddd      !real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)      
cddd      real*8,intent(in) :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)
cddd      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
cddd      real*8,intent(in) :: rh   !Relative humidity
cddd      real*8,intent(in) :: Pa   !Pressure (Pa)
cddd      real*8,intent(in) :: Sacclim ! state of acclimation/frost hardiness [deg C]
cddd      !real*8,intent(in) :: gb !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8 :: ci              !Leaf internal CO2 mole fraction  (umol mol-1)
cddd      real*8,intent(out) :: gs !Leaf stomatal conductance (mol-H2O m-2 s-1)
cddd      real*8,intent(out) :: Anet !Leaf net photosynthesis (CO2 uptake, micromol m-2 s-1)
cddd      !-----Local-----------------
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: stressH2O 
cddd!      type(metdatatype) :: mdat
cddd      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa)
cddd      real*8 :: Atot            !Gross assimilation (umol m-2 s-1)
cddd!      real*8 :: Rd              !Dark or maintenance respiration (umol m-2 s-1)
cddd      !ci should be intialized. For rule of thumb, initialize to
cddd      ! initial ci = 0.7*ca
cddd
cddd!      call load_metdata(IPAR,ca,Ta,rh,O2conc)
cddd
cddd      stressH2O = 1.d0 !Dummy no stress
cddd      call calc_Pspar(dtsec, pft, Pa, Tl, O2pres, stressH2O, Sacclim)
cddd
cddd      Atot = Farquhar(IPAR,ci*Pa*1.d-06,O2pres,Tl,pspar) 
cddd      Anet = Atot - Respveg(pspar%Nleaf,Tl)
cddd
cddd      gs = BallBerry(Anet, rh, cs, pspar)
cddd
cddd!      if (if_ci.eq.1) ci = calc_ci(ca, gb, gs, Anet, IPAR, pspar)
cddd
cddd      !Solving:
cddd      !    Stomates adjust to control the gradient of ci to cs, so the coupled
cddd      !photosynthesis/conductance model is generally solved by iterating
cddd      !to obtain consistent Anet, gs, cs, and ci, and also rh at the leaf 
cddd      !surface.  
cddd      !    As this can be computationally intensive, an alternative is to
cddd      !place bounds on how quickly gs can change (biologically measured).
cddd      !Anet would not be affected, since it depends only on light, Vcmax, and
cddd      !ci; however, ci would be affected.
cddd
cddd      
cddd      end subroutine Collatz

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

cddd      function Farquhar(IPAR,Cip,O2,Tl,pspar) Result(Assim)
cddd!@sum Photosynthesis at the leaf level (umol m-2 s-1)
cddd!@sum after Farquhar, et al.(1980) Planta, 149:78-90.
cddd      implicit none
cddd      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
cddd      real*8,intent(in) :: Cip  !Leaf internal CO2 partial pressure (Pa)
cddd      real*8,intent(in) :: O2   !Leaf internal O2 partial pressure (Pa)
cddd      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: Assim !Leaf GPP (positive into leaf) (umol m-2 s-1)
cddd      !----Local-----
cddd      real*8 :: Je1, Jc1, Js1
cddd      !Assume pspar already calculated.
cddd      !call calc_Pspar(pft, Tl, pspar)
cddd
cddd      Je1 = Je(IPAR,Cip,pspar)
cddd      Jc1 = Jc(Cip,O2,pspar)
cddd      Js1 = Js(pspar)
cddd      Assim = min(Je1, Jc1, Js1)
cddd!#ifdef DEBUG
cddd!      write(997,*) IPAR, Cip, O2,pspar%Gammastar,pspar%Vcmax,
cddd!     &     Tl,Je1,Jc1,Js1,Assim
cddd!#endif
cddd      !Assim = min(Je(IPAR,Cip,pspar), Jc(Cip,O2,pspar), Js(pspar))
cddd
cddd      end function Farquhar
cddd
cddd!-----------------------------------------------------------------------------
cddd
cddd
cddd      function Farquhar_standalone(dtsec,pft,IPAR,ca,Tl,rh,Pa,gb,ci,gs
cddd     i    ,Sacclim )
cddd     o     Result(Anet)
cddd!@sum Photosynthesis at the leaf level.
cddd!@sum after Farquhar, et al.(1980) Planta, 149:78-90.
cddd!@sum Standalone version.  Can be called directly from external program.
cddd
cddd      implicit none
cddd      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
cddd      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1) 
cddd      real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)
cddd      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
cddd      real*8,intent(in) :: rh   !Relative humidity
cddd      real*8,intent(in) :: Pa   !Pressure (Pa)
cddd      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8,intent(in) :: Sacclim ! state of acclimation/frost hardiness
cddd      real*8,intent(in) :: dtsec ! timestep size [sec]
cddd      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction  (umol mol-1)
cddd      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
cddd      real*8 :: Anet            !Leaf net photosynthesis (CO2 uptake, micromol m-2 s-1)
cddd      !-----Local-----------------
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: stressH2O 
cddd!      type(metdatatype) :: mdat
cddd      real*8,parameter :: O2conc=209 !O2 mole fraction in leaf (mmol mol-1)
cddd      real*8 :: cs     !CO2 concentration at leaf surface (umol mol-1)
cddd
cddd      !ci should be intialized. For rule of thumb, initialize to
cddd      ! initial ci = 0.7*ca
cddd
cddd      cs = ca      !Assign CO2 concentration at leaf surface
cddd                   !More rigorously, cs should also be solved for.
cddd      stressH2O = 1.d0 !Dummy no stress
cddd      call calc_Pspar(dtsec,pft,Pa,Tl,O2conc*Pa*1.d-06,stressH2O
cddd     i         ,Sacclim)
cddd      Anet = Farquhar(IPAR,ci*Pa*1.d-06,O2conc*Pa*1.d-06,Tl,pspar) 
cddd     &     - Respveg(pspar%Nleaf,Tl)
cddd
cddd      end function Farquhar_standalone
cddd!-----------------------------------------------------------------------------
cddd

cddd      function calc_ci(ca,gb, gs,Anet,IPAR,pspar) Result(ci)
cddd!@sum Leaf internal CO2 conc (umol mol-1) assuming diffusive flux of CO2
cddd!@sum is at steady-state with biochemical uptake by photosynthesis and
cddd!@sum that there is zero leaf boundary layer resistance (infinite gb),
cddd!@sum and that there is no leaf cuticular conductance of CO2.
cddd!@sum Full equation:  ci = ca - Anet*(1.37/gb + 1.65/gs)
cddd!@sum 1.37 = ratio of diffusivities of CO2 and water vapor in laminar flow
cddd!@sum       in the leaf boundary layer
cddd!@sum 1.65 = ratio of diffusivities of CO2 and water vapor in still air at
cddd!@sum       the leaf surface
cddd!@sum (Monteith, 1995;  Kiang, 2002)
cddd      implicit none
cddd
cddd      real*8 :: ca !Ambient air CO2 mole fraction at surface reference height (umol mol-1)
cddd      real*8 :: gb !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
cddd      real*8 :: gs !Stomatal conductance of water vapor(mol m-2 s-1)
cddd      real*8 :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
cddd      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
cddd      type(photosynthpar) :: pspar
cddd      real*8 :: ci !Leaf internal CO2 mole fraction (umol mol-1)
cddd      !----Local------
cddd      real*8,parameter :: MINPARMOL=50.d0  !(umol m-2 s-1)
cddd
cddd      if (IPAR.lt.MINPARMOL) then  !Stomates closed
cddd        ci = ca - Anet*1.37d0/pspar%b
cddd      else
cddd        ci = ca - Anet*(1.37d0/gb + 1.65d0/gs)
cddd      endif
cddd
cddd      if (ci.lt.ciMIN) ci = ciMIN  !Keep positive definite.
cddd
cddd      end function calc_ci

!-----------------------------------------------------------------------------
!*************************************************************************
      !##### Due to dependency issues, this function is repeated in phenology.f
      !##### Need to put common functions in a different module for both.
      real*8 function frost_hardiness(Sacclim) Result(facclim)
!@sum frost_hardiness.  Calculate factor for adjusting photosynthetic capacity
!@sum  due to frost hardiness phenology.
      real*8,intent(in) :: Sacclim 
      !----Local-----
      real*8,parameter :: Tacclim=-5.93d0 ! threshold temperature for photosynthesis [deg C]
                        ! Site specific thres. temp.: state of photosyn.acclim
                        ! Hyytiala Scots Pine, -5.93 deg C Makela et al (2006)
      real*8,parameter :: a_const=0.0595 ! factor to convert from Sacclim [degC] to facclim [-]
                        ! Site specific; conversion (1/Sacclim_max)=1/16.8115
                        ! estimated by using the max S from Hyytiala 1998
!      real*8 :: facclim ! acclimation/frost hardiness factor [-]

      if (Sacclim > Tacclim) then ! photosynthesis occurs 
         facclim = a_const * (Sacclim-Tacclim) 
         if (facclim > 1.d0) facclim = 1.d0
!      elseif (Sacclim < -1E10)then !UNDEFINED
      elseif (Sacclim.eq.UNDEF)then !UNDEFINED
         facclim = 1.d0   ! no acclimation for this pft and/or simualtion
      else
         facclim = 0.01d0 ! arbitrary min value so that photosyn /= zero
      endif

      end function frost_hardiness
!-----------------------------------------------------------------------------
      real*8 function par_phenology(pft,llspan) Result(fparlimit)  
      integer, intent(in) :: pft
      real*8, intent(in) :: llspan
      real*8, parameter :: vc_tran =12.d0 !9.d0 ! 7.2     !transition
      real*8, parameter :: vc_slop = 15.d0 !10.d0 !16.9    !slope
      real*8, parameter :: vc_amp = 15.d0 !30.d0 !29.8    !amplitude
      real*8, parameter :: vc_min = 10.d0 !25.d0 !7.7     !minimum

      if (llspan > 0.d0) then
         fparlimit = (vc_amp/(1.d0+(llspan/vc_tran)**vc_slop)+vc_min)
     &               /pftpar(pft)%Vcmax
      else
         fparlimit = 1.d0
      endif

      end function par_phenology        
!-----------------------------------------------------------------------------
      end module photcondmod

