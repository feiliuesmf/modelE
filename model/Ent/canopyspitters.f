!#define  DEBUG  1

      module biophysics !canopyspitters
!@sum Spitters (1986) canopy radiative transfer (sunlit/shaded leaves), and
!@sum Simpson's Rule (Price, Numerical Recipes) for canopy layering.
!@sum Call photosynthesis/conductance routines from other module.

      use ent_types
      use ent_const
      use ent_pfts
      use photcondmod, only : pscondleaf, ciMIN
      use FarquharBBpspar

      implicit none
      
      public photosynth_cond !This is interface for Ent.
      !public canopyfluxes

!      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
!      real*8,parameter :: EPS=1.d-8   !Small, to prevent div zero.
      real*8,parameter :: IPARMINMOL=50  !umol m-2 s-1
      real*8,parameter :: O2pres=20900. !Pa, from Collatz? Hanan? Not exactly .209*101325 Pa.

      type :: canraddrv
         !Canopy radiative transfer
         !@var sigma Leaf scattering coefficient (?unitless).
         real*8 :: sigma        !=0.2D0
         real*8  sqrtexpr  !This just calculated from sigma.
         !@var kdf Canopy extinction coeff. for diffuse radiation (unitless).
         real*8 :: kdf          !=0.71D0
         !@var rhor Canopy reflectivity (?unitless).
         real*8 rhor
         !@var kbl Canopy extinction coeff. for direct radiation (unitless).
         real*8 kbl

         !Vegetation geometry, biology
         !UPDATE WITH CANOPY RADIATIVE TRANSFER SCHEME
         integer :: pft
         !real*8 :: leafalbedo
         real*8 :: canalbedo  
         real*8 :: LAI

         !Radiation
!         real*8 :: solarzen !Solar zenith angle (radians)
         real*8 :: CosZen !cos(solarzen) = sin(solarelev)
         real*8 :: I0df, I0dr   !Direct and diffuse incident PAR at top of canopy (umol m-2 s-1)
      end type canraddrv

      contains
!################## MAIN SUBROUTINE #########################################
      subroutine photosynth_cond(dtsec, pp)
      !@sum Farquhar-Ball-Berry version of photosynth_cond.
      !@sum Calculates photosynthesis, autotrophic respiration, conductance,
      !@sum looping through cohorts.
      !@sum Inputs:  met drivers, radiation, Cf, Tcanopy
      !@sum Outputs:  GPP, NPP, respiration components
      use ent_const
      use ent_types
      use FarquharBBpspar !pspartype, psdrvtype
      use photcondmod, only : biophysdrv_setup, calc_Pspar
      use patches, only : patch_print
!      use physutil, only:  QSAT

      implicit none

      real*8, intent(in) :: dtsec
      type(patch),pointer :: pp  
      !----Local----------------!
      type(cohort),pointer :: cop
      type(psdrvtype) :: psdrvpar !Met biophysics drivers, except for radiation.
      real*8 :: ci_umol !umol/mol, Leaf internal CO2 
      real*8 :: ca_umol !umol/mol, Ambient air CO2
      real*8 :: cf_umol !umol/mol, Foliage surface CO2 !NOTE: This is only used by Friend & Kiang biophysics.
      real*8 :: TsurfK, TcanK,Pa !,rh
      real*8 :: CosZen,betad
      real*8 :: IPAR            !Incident PAR 400-700 nm (W m-2)
      real*8 :: fdir            !Fraction of IPAR that is direct
      real*8 :: Gb !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Anet,Atot,Rd    !umol m-2 s-1
      real*8 :: GCANOPY,TRANS_SW ! Ci,NPP !,R_auto
      real*8 :: GCANOPYsum, Ciavg, GPPsum, NPPsum, R_autosum,C_labsum,
     &          R_rootsum  !PK 5/15/07
      real*8 :: molconc_to_umol

#ifdef DEBUG
      print *,"Started photosynth_cond in FBB" ! with patch:"
      !call patch_print(6,pp," ")
#endif

      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
        pp%TRANS_SW = 1.d0
        return
      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *,"photosynth_cond: wrong pft = ", pp%tallest%pft
        call patch_print(6,pp,"ERROR ")
        call stop_model("photosynth_cond: wrong pft",255)
      endif

      !* ZERO SOME OUTPUT VARIABLES AT PATCH LEVEL
      pp%TRANS_SW = 1.d0 !Case of zero LAI.
      !* Time-stepped outputs:  CNC, Ci, Qf.

      !* INITIALIZE SUMMARY OUTPUT VARIABLES *!
      GCANOPYsum = 0.d0
      Ciavg = 0.d0
      GPPsum = 0.d0
      NPPsum = 0.d0
      R_autosum = 0.d0
      R_rootsum = 0.d0
      C_labsum = 0.d0


      !* SET UP DRIVERS *!
      !* Patch-level water stress only needed for Friend&Kiang conductance.
      !* Cohort-level water stress is used for cohort-level photosynthesis.
!      pp%betad = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
!     i     ,pp%fracroot(:)
!     i     ,pp%cellptr%fice(:), pfpar(pp%tallest%pft)%hwilt
!     o     , pp%betadl(:))

      !* Radiation drivers *!
      IPAR = pp%cellptr%IPARdir + pp%cellptr%IPARdif
      if (pp%cellptr%IPARdir.eq.0.d0) then
        fdir = 0.d0
      else
        fdir = pp%cellptr%IPARdir / IPAR
      endif
      CosZen = pp%cellptr%CosZen

      Pa = pp%cellptr%P_mbar * 100.d0

      !* Other photosynthesis drivers *!
      !Set up psdrvpar - pack in met drivers.
      Gb = pp%cellptr%Ch*pp%cellptr%U* Pa/
     &     (gasc*(pp%cellptr%TairC+KELVIN)) !m/s * N/m2 * mol K/J * 1/K = mol/m2/s
!      write(995,*) "canopyspitters Gb:",Gb,pp%cellptr%Ch,pp%cellptr%U,
!     &     Pa,pp%cellptr%TairC,gasc
      molconc_to_umol = gasc * (pp%cellptr%TcanopyC + KELVIN)/Pa * 1d6
      ci_umol = pp%cellptr%Ci * molconc_to_umol
      ca_umol = pp%cellptr%Ca * molconc_to_umol !This isn't used, but may be to calculated next cf.
      cf_umol = pp%cellptr%Cf * molconc_to_umol !##HACK - NEED TO CALC AND SAVE Cf
      TcanK = pp%cellptr%TcanopyC + KELVIN
      TsurfK = pp%cellptr%TairC + KELVIN
      call biophysdrv_setup(cf_umol,ci_umol,
     &     pp%cellptr%TcanopyC,Pa,
     &     min(1.d0,  !RH
     &     pp%cellptr%Qf/QSAT(TcanK*(101325.d0/Pa)**(gasc/cp),
     &     2500800.d0 - 2360.d0*(TsurfK-KELVIN),Pa/100.d0)),
     &     psdrvpar)
!      psdrvpar%cf = cf_umol !##HACK - NEED TO CALC AND SAVE Cf
!      psdrvpar%Tc = pp%cellptr%TcanopyC
!      psdrvpar%Pa = Pa
!      psdrvpar%rh = min(1.d0,
!     &     pp%cellptr%Qf/QSAT(TcanK*(101325.d0/Pa)**(gasc/cp),
!     &     2500800.d0 - 2360.d0*(TsurfK-KELVIN),Pa/100.d0))
!      Ch = pp%cellptr%Ch
!      U = pp%cellptr%U

      
      !* LOOP THROUGH COHORTS *!
      cop => pp%tallest
      do while (ASSOCIATED(cop))
        call calc_Pspar(cop%pft, psdrvpar%Pa, psdrvpar%Tc,O2pres) 

        !* Assign vegpar

        if (cop%LAI.gt.0.d0) then
          cop%stressH2O = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
     i         ,cop%fracroot(:)
     i         ,pp%cellptr%fice(:), pfpar(cop%pft)%hwilt
     o         , cop%stressH2Ol(:))
          betad = cop%stressH2O

          call canopyfluxes(dtsec, cop%pft
     &         ,pp%albedo(1),pp%LAI,cop%LAI,IPAR*4.05 !4.05 see canopyfluxes comments.
     &         ,CosZen,fdir
     &         ,Gb,ca_umol,cf_umol
     &         ,psdrvpar
     &         ,GCANOPY,Anet,Atot,Rd !NOTE: Ci should be cohort level
     &         ,TRANS_SW)       !NOTE:  Should include stressH2O.
!     &       ,if_ci)  

         !* Assign outputs to cohort *!
         !* Multiply in cop%stressH2O
!          cop%GCANOPY = GCANOPY * 18.d-3 / rhow !Convert mol-H2O m-2 s-1 to m/s
          cop%GCANOPY = GCANOPY * (gasc*TsurfK)/Pa  !Convert ?mol-H2O m-2 s-1 to m/s
         ! (mol-H2O m-2 s-1)*(18 g/mol * 1d-3 kg/g) / (rhow kg m-3) = m/s
          cop%Ci = psdrvpar%ci * !ci is in mole fraction
     &         psdrvpar%Pa/(gasc * (psdrvpar%Tc+KELVIN)) !mol m-3
          cop%GPP = Atot * 0.012d-6 !umol m-2 s-1 to kg-C/m2-ground/s
        else
          cop%GCANOPY=0.d0 !May want minimum conductance for stems & cuticle.
          cop%Ci = EPS
          cop%GPP = 0.d0
        endif
        !* Update cohort respiration components, NPP, C_lab
        !## Rd should be removed from pscondleaf, only need total photosynthesis to calculate it.
        call Respiration_autotrophic(dtsec, psdrvpar%Tc, cop)

        !* pp cohort flux summaries
        GCANOPYsum = GCANOPYsum + cop%GCANOPY
        Ciavg = Ciavg + cop%Ci*cop%LAI
        GPPsum = GPPsum + cop%GPP
        NPPsum = NPPsum + cop%NPP
        R_autosum = R_autosum + cop%R_auto
        R_rootsum = R_rootsum + cop%R_root  !PK 5/15/07
        C_labsum = C_labsum + cop%C_lab * cop%n !Sum for cohort.

        cop => cop%shorter
      end do

      !* Patch-level OUTPUTS *!
      pp%GCANOPY = GCANOPYsum
      pp%Ci = Ciavg/pp%LAI
      pp%GPP = GPPsum
      pp%NPP = NPPsum
      pp%R_auto = R_autosum
      pp%R_root = R_rootsum
      pp%TRANS_SW = TRANS_SW 

      !* Accumulate uptake. 
      !* Respiration should be from leaves and not draw down C_lab. ## Need to allocate respiration to leaves.##
!      pp%C_lab = pp%C_lab + max(C_labsum, 0.d0)  !(kg/m2) ###Eventually need to convert to kg/individual.
      pp%C_lab = C_labsum!(kg/m2) ###Eventually need to convert to kg/individual.
      
      end subroutine photosynth_cond


!---------------------------------------------------------------------------
      subroutine canopyfluxes(dt, pft
     i     ,canalbedo,LAIcanopy,LAIcohort,IPAR,CosZen,fdir
     i     ,Gb,ca,cf
     i     ,psdrvpar
     o     ,Gs,Anet,Atot,Rd,TRANS_SW)
!     i     ,if_ci)
!@sum canopyfluxes Calculates photosynthesis and conductance with
!@sum Farqhuar et al. (1980) photosynthesis, Ball-Berry stomatal conductance,
!@um  and Spitters (1986, 1987) canopy radiation (sunlit, shaded leaves).
!@sum Integrates vertically over the canopy with Simpson's Rule.
!@sum Ci is updated at the canopy level using the canopy boundary layer
!@sum conductance as in Friend and Kiang (2005). 

!@sum If PAR is not directly available, the following conversions may be used:
      !  From total shortwave (W m-2) to PAR (umol m-2 s-1) (Monteith & Unsworth):
      !          PAR(umol m-2 s-1) = 2.3(umol/J)*SW(W m-2)
      !  From PAR (W m-2) to PAR (umol m-2 s-1) (U.Maryland, Dept. of Met., PAR Project),
      !  suggest nominal 485 nm for conversion, which gives:
      !          PAR(umol m-2 s-1) = 4.05(umol/J) * PAR(W m-2)
      !  Dye, D.G. (2004) suggests a slightly different conversion:
      !          PAR(umol m-2 s-1) = 4.56(umol/J) * PAR(W m-2)
      
      implicit none
      real*8,intent(in) :: dt   !time step (seconds)
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: canalbedo !Whole-canopy albedo
      real*8,intent(in) :: LAIcanopy  !Whole-canopy LAI
      real*8,intent(in) :: LAIcohort  !cohort LAI
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
!      real*8,intent(in) :: solarzen !Solar zenith angle
      real*8,intent(in) :: CosZen !Cosine ofSolar zenith angle
      real*8,intent(in) :: fdir !Fraction of PAR that is direct
      real*8,intent(in) :: Gb   !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: ca   !CO2 mole fraction at reference height (umol mol-1)
      real*8,intent(in) :: cf   !CO2 mole fraction at canopy surface (umol mol-1)
!      real*8,intent(in) :: Tc   !Canopy (foliage) temperature (Celsius)
!      real*8,intent(in) :: rh   !Relative humidity (fraction)
!      real*8,intent(in) :: Pa   !Atmospheric pressure (Pa)
      type(psdrvtype) :: psdrvpar !Photosynthesis drivers, except for radiation.
!      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction (umol mol-1)
      real*8,intent(inout) :: Gs !Canopy stomatal conductance of water vapor (mol m-2 s-1)
      real*8,intent(out) :: Anet !Leaf net photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Leaf respiration (umol/m2/s)
      real*8,intent(out) :: TRANS_SW !Transmittance of shortwave to ground surface.
!      integer,intent(in) :: if_ci !FLAG 0-don't calculate ci, 1-calculate ci
      
      !-----Local---------------------
      real*8 :: Gsint           !Gs canopy conductance from qsimp integration (mol/m2/s)
      real*8 :: Rdint           !Rd canopy foliage respiration form qsimp integration (umol/m2/s)
      type(canraddrv) :: cradpar
      !real*8 :: Ciconc          !Leaf internal CO2 mole fraction calculate (umol mol-1)


      !Set up cradpar and calculate diffuse and direct PAR.
      call canopy_rad_setup(pft,CosZen,fdir,IPAR,
     &     LAIcanopy,canalbedo,cradpar) 

      !Set up initial conditions for Collatz photosynthesis/conductance.
      !call init_Collatz(cf, ci)

      !Calculate net photosynthesis, integrate vertically with Simpson's Rule.
      !call qsimp(cradpar%LAI,cradpar,cf,ci,Tc,Pa,rh,Anet,Gsint) 
      !### LAIcanopy for radiation and LAIcohort for photosynthesis need to be distinguished.
      call qsimp(cradpar%LAI,cradpar,psdrvpar,ca,Gb,Atot,Gsint,Rdint) 
!      write(993,*) cradpar,psdrvpar,ca,cf,Gb,Atot,Gsint,Rdint
      !Calculate conductance and ci at the canopy level
!      call Gs_bound(dt, LAI,Gsint, Gs) !Limit rate of change of Gs.
!      Ciconc = calc_Ci_canopy(ca,Gb,Gs,Anet,LAI,IPAR)

      !Different options to solve for Gs and Ci:
      !If Ciconc < CiMIN, then Gs was too big, so limit Gs to allow
      !only Ciconc >= CiMIN
!      if (Ciconc < CiMIN) then
!        Ciconc = CiMIN
!        Gs =  Gs_from_Ci(Anet, ca, Gb, Gs, ci,IPAR) 
!      else
!        ci = Ciconc
!      endif

      Gs = Gsint
      Rd = Rdint
      Anet = Atot - Rd
      call canopy_transmittance(TRANS_SW,CosZen,fdir,cradpar)
      !Return:  Ci, Gs, Anet
!#ifdef DEBUG        
!        write(92,*) CosZen,IPAR,cradpar,psdrvpar
!     &       ,Gb,Gsint,Gs,Atot,Anet,Rd,TRANS_SW
!#endif
      end subroutine canopyfluxes

!---------------------------------------------------------------------!
      function water_stress(nlayers, soilmp, fracroot, fice,
     &     hwilt, betadl) Result(betad)
      !1. Rosensweig & Abramopoulos water stress fn.

      implicit none
      integer,intent(in) :: nlayers !Number of soil layers
      real*8,intent(in) ::  soilmp(:) !Soil matric potential (m)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(in) :: hwilt  !Wilting point of pft, matric pot. (m)
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress
      !---Local-----------
      integer :: k
      
      betad = 0.d0
      do k = 1,nlayers
        betadl(k) = (1.d0-fice(k))*fracroot(k)
     &       *max((hwilt-soilmp(k))/hwilt,0.d0)
        betad = betad + betadl(k) 
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress

!################## PHOTOSYNTHESIS #########################################

      subroutine photosynth_sunshd(
      !Spitters parameters
     i     Lcum                 !Cumulative LAI from top of canopy (m2/m2)
     i     ,crp                !Canopy radiation parameters 
     i     ,psp                 !Photosynthesis met drivers
     i     ,ca                  !Ambient air CO2 conc (umol/mol)
     i     ,Gb                  !Leaf boundary layer conductance (mol/m2/s)
     o     ,Alayer              !Leaf Net assimilation of CO2 in layer (umol m-2 s-1)
     o     ,gslayer            !Leaf Conductance of water vapor in layer (mol m-2 s-1)
     o     ,Rdlayer)             !Leaf respiration (umol m-2 s-1)
      implicit none
      real*8,intent(in) :: Lcum
      type(canraddrv) :: crp
      type(psdrvtype) :: psp
      real*8,intent(in) :: ca, Gb
      !real*8,intent(in):: cs,Tl,Pa,rh
      !real*8,intent(inout) :: ci
      real*8,intent(out) :: Alayer !Flux for single leaf
      real*8,intent(out) :: gslayer !Conductance for single leaf
      real*8,intent(out) :: Rdlayer !Respiration for single leaf
      !------Local---------
      real*8 fsl  !Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 Isl !PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 Ish !PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 Asl !Anet from sunlit leaves (umol/m2/s)
      real*8 Ash !Anet from shaded leaves (umol/m2/s)
      real*8 gssl,gssh !Leaf conductance, sunlit,shaded (mol-H2O/m2/s)
      real*8 Rdsl, Rdsh
      integer :: sunlitshaded !1-sunlit, 2-shaded

      call canopy_rad(Lcum,crp,Isl,Ish,fsl)
!      write(997,*) 'Lcum,sigma, sqrtexpr,kdf,rhor,kbl,pft,canalbedo,
!     & LAI,solarzen,I0df,I0dr,Isl,Ish,fsl',Lcum,crp,Isl,Ish,fsl

      !Calculate photosynthesis and stomatal conductance.
!      write(991,*) 'sunlit'
      sunlitshaded = 1
      call pscondleaf(crp%pft,Isl,psp,ca,Gb,gssl,Asl,Rdsl,sunlitshaded)
!      write(992,*) 'shaded'
      sunlitshaded = 2
      call pscondleaf(crp%pft,Ish,psp,ca,Gb,gssh,Ash,Rdsh,sunlitshaded)
      !call Collatz(crp%pft, Isl,cs,Tl,rh, Pa,ci,gssl,Asl)
      !call Collatz(crp%pft, Ish,cs,Tl,rh, Pa,ci,gssh,Ash)
                   
      Alayer = fsl*Asl + (1.0d0 - fsl)*Ash
      gslayer = fsl*gssl + (1.0d0 - fsl)*gssh
      Rdlayer = fsl*Rdsl + (1.0d0 - fsl)*Rdsh
!#ifdef DEBUG
!      write(998,*) Lcum,crp,Isl,Ish,fsl,psp,gssl,gssh,Asl,Ash,Rdsl,Rdsh
!#endif
      end subroutine photosynth_sunshd

!################# CANOPY CONDUCTANCE ######################################
      function calc_Ci_canopy(Ca,Gb, Gs,Anet,LAI,IPAR) Result(ci)
!@sum Foliage internal CO2 conc (mol mol-3) assuming diffusive flux of CO2
!@sum is at steady-state with biochemical uptake by photosynthesis and
!@sum that there is zero leaf boundary layer resistance (infinite gb),
!@sum and that there is no leaf cuticular conductance of CO2.
!@sum Full equation:  ci = ca - Anet*(1.6/gb + 1.4/gs)
!@sum 1.6 = ratio of diffusivities of CO2 and water vapor in laminar flow
!@sum       in the leaf boundary layer
!@sum 1.4 = ratio of diffusivities of CO2 and water vapor in still air at
!@sum       the leaf surface
!@sum (Monteith, 1995;  Kiang, 2003)
!@sum ##### NOTE: Parameter Ball_b should actually be passed in with pspar.
!@sum #####       Have to set up for generic pspar for different photosynthesis
!@sum #####       routines.  (NK)

      implicit none

      real*8,intent(in) :: ca !CO2 mole fraction at surface reference height (umol mol-1)
      real*8,intent(in) :: Gb !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: Gs !Canopy Stomatal conductance of water vapor(mol m-2 s-1)
      real*8,intent(in) :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8,intent(in) :: LAI !Leaf area index 
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: ci              !Leaf internal CO2 concentration (umol mol-1)
      !----Local------
      real*8,parameter :: MINPARMOL=50  !umol m-2 s-1
      real*8,parameter :: Ball_b = 0.01 !mol m-2 s-1

      if (IPAR.lt.MINPARMOL) then  !Stomates closed
        ci = ca - Anet*1.6/Ball_b
      else
        ci = ca - Anet*(1.6/Gb + 1.4/Gs) !LAI cancels in numerator and denominator.
      endif


      end function calc_Ci_canopy

!---------------------------------------------------------------------------
      subroutine Gs_bound(dt, LAI, Gsnew, Gsinout)  
      !@sum Gs_bound Limit rate of change of Gs (umol m-2 s-1)
      ! Required change in canopy conductance to reach equilibrium (m/s).
      real*8,intent(in) :: dt, LAI
      real*8,intent(in) :: Gsnew !Canopy conductance of curr time step (mol m-2 s-1)
      real*8,intent(inout) :: Gsinout !Bounded canopy conductance (mol m-2 s-1)
      !---Local----!
      real*8 :: Gsold           !Canopy conductance at prev time step (mol m-2 s-1)
      real*8, parameter :: rhoH2O = 998.2 !Density of water (1000 kg/m^3)
      real*8, parameter :: MW_H2O = 18.015 !Molecular weight of water (g/mol)
      !real*8, parameter :: ghi = 0.006d0*rhoh2o*1000./MW_H2O !Upper limit of gs leaf (mol m-2 s-1)
      !real*8, parameter :: glo = 0.000001d0*rhoH2O*1000./MW_H2O !Lower limit of gs leaf (mol m-2 s-1), See Ball and Berry paper.
      real*8, parameter :: ghi = 333.0 !Conversion from 6 mm s-1 upper limit.(mol m-2 s-1)
      real*8, parameter :: glo = .015 !Temperature grassland. Korner (1994) (mol m-2 s-1)

      real*8 :: dGs, dGs_max 

      dGs=Gsnew-Gsinout
      Gsold = Gsinout
      Gsinout = Gsnew
      !nu Limit Gs change over timestep because of guard cell mechanics (m/s)
      dGs_max=dt*LAI*(ghi-glo)/1800.0D0
      if( dGs.gt.dGs_max) Gsinout = Gsold + dGs_max
      if(-dGs.gt.dGs_max) Gsinout = Gsold - dGs_max
      ! Biological limits of absolute Gs (m/s).
      if(Gsinout.gt.ghi*LAI) Gsinout=ghi*LAI
      if(Gsinout.lt.glo*LAI) Gsinout=glo*LAI

      end subroutine Gs_bound

!---------------------------------------------------------------------------
      function Gs_from_Ci(Anet,Ca,Gb,Gs,Ci,IPAR) Result(gsout)
      !Inversion of calc_Ci_canopy
      real*8 :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8 :: ca !Ambient air CO2 mole fraction at surface reference height (umol mol-1)
      real*8 :: gb !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: gs !Canopy Stomatal conductance of water vapor(mol m-2 s-1)
      real*8 :: ci !Leaf internal CO2 concentration (umol mol-1)
      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: gsout !Adjusted Gs (mol m-2 s-1)
      !----- Local ------

      if (IPAR.lt.IPARMINMOL) then  !Stomates closed
        gsout = gs
      else
        gsout =  1.4/((Ca - Ci)/(Anet) - 1.6/gb )
      endif

      end function Gs_from_Ci

!################# AUTOTROPHIC RESPIRATION ######################################

      subroutine Respiration_autotrophic(dtsec,TcanopyC,cop)
      !@sum Autotrophic respiration - updates cohort respiration,NPP,C_lab
      !@sum Returns kg-C/m^2/s

      implicit none
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: TcanopyC
      type(cohort),pointer :: cop
      !----Local-----
      real*8 :: Resp_maint, Cdiff

      !NOTE: NEED TO FIX Canopy maintenance respiration for different
      !C:N ratios for the different pools.
      !* Maintenance respiration - root
      cop%R_root = 0.012D-6 * Resp_can_maint(cop%pft,cop%C_froot,
     &     pfpar(cop%pft)%lit_C2N,TcanopyC+Kelvin,cop%n)
      !* Maintenance respiration - leaf + sapwood + storage
      Resp_maint = 0.012D-6 * ( !kg-C/m2/s
!     &       Canopy_resp(vegpar%Ntot, TcanopyC+KELVIN) !Foliage
     &     Resp_can_maint(cop%pft, cop%C_fol,pfpar(cop%pft)%lit_C2N,
     &     TcanopyC+Kelvin,cop%n) !Foliage
     &     + Resp_can_maint(cop%pft,cop%C_sw, !Sapwood - need to get C:N
     &     900.d0,TcanopyC+Kelvin,cop%n) 
     &     + Resp_can_maint(cop%pft,cop%C_lab, !Storage
     &     pfpar(cop%pft)%lit_C2N,TcanopyC+Kelvin,cop%n)  )
      !* Total respiration : maintenance + growth
      cop%R_auto =  
     &     Resp_maint 
     &     + cop%R_root         !PK 5/15/07 
     &     + 0.012D-6 * Resp_can_growth(cop%pft, cop%GPP/0.012D-6,
     &     (Resp_maint+cop%R_root)/0.012D-6 )
!     &       Canopy_resp(vegpar%Ntot, TcanopyC+KELVIN))
      cop%NPP = cop%GPP - cop%R_auto !kg-C/m2-ground/s

      !* Accumulate uptake. - * PUT THIS IN A SUBROUTINE OR FUNCTION C_lab_update
      if ( (cop%NPP*dtsec/cop%n.lt.0.d0).and.
     &     (abs(cop%NPP*dtsec/cop%n).ge.cop%C_lab) ) then 
        !Don't let C_lab go below zero.
        cop%C_lab = EPS
        !* NYK - TEMPORARY HACK - SHOULD CALL STRESS SUBROUTINE FOR SENESCENCE.
        !*       Needs checks for negative NPP loss greater than C pools.
        Cdiff = cop%NPP*dtsec/cop%n + (cop%C_lab - EPS)
        cop%C_fol = cop%C_fol + 0.33d0*Cdiff
        cop%C_froot = cop%C_froot + 0.33d0*Cdiff
        cop%C_sw = cop%C_sw + 0.33d0*Cdiff
        !cop%LAI = !should update
      else
        cop%C_lab = cop%C_lab + cop%NPP*dtsec/cop%n !(kg/individual)
      endif

!      write(96,*) TcanopyC,cop%GPP,cop%R_root,cop%R_auto,cop%C_fol
      end subroutine Respiration_autotrophic

!---------------------------------------------------------------------!
      real*8 function Resp_can_maint(pft,C,CN,T_k,n) 
     &     Result(R_maint)
      !Canopy maintenance respiration (umol/m2-ground/s)
      !Based on biomass amount (total N in pool). From CLM3.0.

      integer :: pft            !Plant functional type.
      real*8 :: C               !g-C/individual  !not per m2 -PK 5/15/07
                                !Can be leaf, stem, or root pools.
      real*8 :: CN              !C:N ratio of the respective pool
!      real*8 :: R_maint !Canopy maintenance respiration rate (umol/m2/s)
      real*8 :: T_k             !Temperature of canopy (Kelvin)
      real*8 :: n               !Density of individuals (no./m2)
      !---Local-------
      real*8,parameter :: k_CLM = 6.34d-07 !(s-1) rate from CLM.

      R_maint = pfpar(pft)%r * k_CLM * (C/CN) * 
     &     exp(308.56d0*(1/56.02d0 - (1/(T_k-227.13d0)))) *
     &     2.d6*n/28.5d0
      !Note:  CLM calculates this per individual*population/area_fraction
      !      to give flux per area of pft cover rather than per ground area.
      end function Resp_can_maint
!---------------------------------------------------------------------!
      real*8 function Resp_root(Tcelsius,froot_kgCm2) Result(Rootresp)
!@sum Frootresp = fine root respiration (kgC/s/m2)
      !From ED model.  Not used.
      real*8 :: Tcelsius, froot_kgCm2
      
      Rootresp = OptCurve(Tcelsius,1.0d0,3000.d0) * froot_kgCm2/SECPY
     &     /((1.d0 + exp(0.4d0*(5.0d0-Tcelsius)))
     &     *(1.d0 + exp(0.4d0*(Tcelsius-45.d0))))

      end function Resp_root
!---------------------------------------------------------------------!
      real*8 function OptCurve(Tcelsius,x,y) Result(OptCurveResult)
!@sum Optimum curve, where OptCurveResult=x at 15 Celsius.
      !From ED model.
      real*8 :: Tcelsius, x, y
      
      OptCurveResult = x * exp(y*(1/288.15d0 - 1/(Tcelsius+KELVIN)))
      end function OptCurve

!---------------------------------------------------------------------!
      real*8 function Resp_can_growth(pft,Acan,Rmaint) Result(R_growth)
      !Canopy growth respiration (umol/m2/s)
      !Based on photosynthetic activity. From CLM3.0.
      !Fixed to min 0.d0 like ED2. - NYK
      integer :: pft
      real*8 :: Acan !Canopy photosynthesis rate (umol/m2/s)
      real*8 :: Rmaint !Canopy maintenance respiration rate (umol/m2/s)
      real*8 :: growth_r !pft-dependent. E.g.CLM3.0-0.25, ED2 conifer-0.53, ED2 hw-0.33

      growth_r = pfpar(pft)%r * 0.5d0  !NYK's hack - chec
      R_growth = max(0.d0, growth_r * (Acan - Rmaint))
      end function Resp_can_growth

!################# RADIATIVE TRANSFER ######################################
      subroutine canopy_rad_setup(
     i     pft, CosZen, fdir,IPAR,LAI,canalbedo,
     o     crp)
!@sum Calculate diffuse and direct PAR, canopy albedo, if unknown, from given
!@sum solar zenith angle, and other canopy radiation parameters.

      implicit none

! NOTES FOR OBTAINING PAR
! If only incident shortwave is available, convert 
!  total shortwave (W/m2) to PAR (umol/m2/s).
!  2.3 umol/J = shortwave to fraction that is PAR (Monteith & Unsworth).
!      PAR=2.3d0*max(srht,zero)/(1.0D0-0.08D0)
! Current: Replaced back-calculation with actual incident PAR at surface.

! GENERAL PAR CONVERSIONS.
! nyk  For 400-700 nm, energy is 3.3-6 umol photons/J.
!      Top-of-atmosphere flux-weighted average of energy is 4.54 umol/J.
!      Univ. Maryland, Dept. of Met., PAR Project, suggests nominal
!      485 nm for conversion, which gives 4.05 umol/J.

      !Input parameters
      integer,intent(in) :: pft
!      real*8,intent(in) :: solarzen !Solar zenith angle
      real*8,intent(in) :: CosZen !Cosine of Solar zenith angle
      real*8,intent(in) :: fdir   !Fraction of shortwave radiation that is direct
      real*8,intent(in) :: IPAR    !PAR incident at top of canopy (umol/m2/s)
      real*8,intent(in) :: LAI    !Leaf area index of whole canopy
      real*8,intent(in) :: canalbedo !Canopy albedo, if known

      !Output parameters
      type(canraddrv) :: crp
      !real*8,intent(out) :: I0df, I0dr !Diffuse and direct PAR incident at top of canopy (umol m-2 s-1)
      !----Local var-------
      real*8 :: sbeta  !Sine of solar zenith angle
      ! CONSTANTS
      real*8,parameter :: sigma=0.2d0  !Leaf scattering coefficient
      real*8,parameter :: kdf=0.71d0  !Canopy extinction coefficient
      real*8,parameter :: EPS=1.d-8   !Small, to prevent div zero.

      sbeta = CosZen ! cos(solarzen) = sin(solarelev)

!      crp%solarzen = acos(CosZen)
      crp%CosZen = CosZen
! Direct beam PAR incident on canopy (umol/m2/s).
      crp%I0dr=fdir*IPAR
! Diffuse PAR incident on canopy (umol/m2/s).
      crp%I0df=(1.0D0-fdir)*IPAR
! Canopy reflectivity depends on sun angle. Or take prescribed albedos.
      crp%sigma = sigma
      crp%kdf = kdf
      crp%sqrtexpr=sqrt(1.0D0-crp%sigma)
! Canopy extinction coefficient for direct beam radiation depends on
! sun angle (unitless).
      crp%kbl=0.5D0*crp%kdf/(0.8D0*crp%sqrtexpr*sbeta+EPS)

      !Hack canopy reflectivity.  This assumes a closed canopy!
      if (canalbedo.eq.0.0) then !Use if radiation gets initialized late.
        crp%rhor=((1.0D0-crp%sqrtexpr)/(1.0D0+crp%sqrtexpr))*
     &       (2.0D0/(1.0D0+1.6D0*sbeta))
        crp%canalbedo = crp%rhor
      else
       !Prescribed veg albedos until have canopy scheme.
        crp%rhor = canalbedo
        crp%canalbedo = crp%rhor
      end if

      crp%pft = pft
      crp%LAI = LAI

      end subroutine canopy_rad_setup

!-----------------------------------------------------------------------------
      subroutine canopy_rad(
!@sum Calculate PAR on sunlit and shaded foliage in layer Lc      
     i     Lc,                  !Cumulative LAI at layer from top of canopy
     i     crp,                !veg parameters
     o     Isla,                !PAR on shaded foliage in layer Lc
     o     Isha,                 !PAR on sunlit foliage in layer Lic
     o     fsl)                 !Fraction of sunlit foliage at Lc (unitless).

      implicit none
! Get incident radiation on layer, diffuse/direct, sunlit/shaded.
      real*8,intent(in) :: Lc
      type(canraddrv) :: crp
!@var Isla PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8,intent(out) :: Isla
!@var Isha PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
      real*8,intent(out) :: Isha
!@var fsl Fraction of sunlit foliage in layer at Lc (unitless).
      real*8,intent(out) :: fsl

      !-------Local vars--------------------------------------------------
      real*8 :: sbeta           !cos(solarzen) = sin(solarelev)
      real*8 :: I0df, I0dr !Diffuse and direct IPAR (umol m-2 s-1)
!@var Idfa Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idfa
!@var Idra Direct PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idra
!@var Idrdra Direct PAR at Lc that remains direct (umol/m[ground]2/s).
      real*8 Idrdra
      real*8,parameter :: zero = 0.d0
      
      sbeta = crp%CosZen !cos(solarzen) = sin(solarelev)
      I0df = crp%I0df
      I0dr = crp%I0dr
! When sun is above the horizon.
      if(sbeta.gt.EPS)then
! Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
        Idfa=(1.0D0-crp%rhor)*I0df*crp%kdf *exp(-crp%kdf*Lc)
! Direct PAR in canopy at Lc (umol/m[ground]2/s).
        Idra=(1.0D0-crp%rhor)*I0dr* crp%sqrtexpr * crp%kbl *
     $        exp(-crp%sqrtexpr*crp%kbl*Lc)
! Direct PAR at Lc that remains direct (umol/m[ground]2/s).
        Idrdra=(1.0D0-crp%sigma)*I0dr * crp%kbl *
     $        exp(-crp%sqrtexpr*crp%kbl*Lc)
! PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
        Isha=Idfa+(Idra-Idrdra)
! PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s)._
        Isla=Isha+(1.0D0-crp%sigma)* crp%kbL * I0dr
! Fraction of sunlit foliage in layer at Lc (unitless).
        fsl=exp(-crp%kbL*Lc)
        if(Isha.lt.0.1D0)Isha=0.1D0
        if(Isla.lt.0.1D0)Isla=0.1D0
        if(fsl.lt.zero)fsl=zero
      else
        Isha=0.1D0
        Isla=0.1D0
        fsl=zero
      endif

      end subroutine canopy_rad

      subroutine canopy_transmittance(TRANS_SW, sbeta,fdir,crp)
! Calculates the transmittance (fraction) of radiation through the canopy
! to the soil surface.      
! Transmission of shortwave radiation through canopy. This has errors.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Parameters -------------------------------------------------------------
!@var trans_sw Total canopy transmittance of shortwave (fraction)
      real*8, intent(out) :: TRANS_SW  
!@var sbeta Cos of solar zenith angle
      real*8, intent(in) :: sbeta
!@var fdir Fraction of surface visible radiation that is direct
      real*8, intent(in) :: fdir
! Canopy radiation parameters
      type(canraddrv) :: crp

!Local variables --------------------------------------------------------
! Get assigned from crp passed in.
      real*8 sigma,sqrtexpr,rhor,kbl,kdf,alai
!@var absdf Absorbance of diffuse SW radiation (fraction)
      real*8 absdf
!@var absdr Absorbance of direct SW radiation (fraction)
      real*8 absdr
!@var absdrdr Absorbance of direct SW that remains direct (fraction)
      real*8 absdrdr
!@var abssh Absorbance of SW through shaded foliage (fraction)
      real*8 abssh
!@var abssl Absorbance of SW through sunlit foliage (fraction)
      real*8 abssl
!@var fracsl Fraction of leaves in canopy that are sunlit (fraction)
      real*8 fracsl
!------------------------------------------------------------------------
      sigma=crp%sigma
      sqrtexpr=crp%sqrtexpr
      rhor=crp%rhor
      kbl=crp%kbl
      kdf=crp%kdf
      alai=crp%LAI

! Diffuse shortwave in canopy (umol/m[ground]2/s).
      if(sbeta.gt.zero)then
        absdf=(1-fdir)*(1.0D0-rhor)*kdf*exp(-kdf*alai)
! Direct shortwave in canopy (umol/m[ground]2/s).
        absdr=fdir*(1.0D0-rhor)*sqrtexpr*kbl*exp(-sqrtexpr*kbl*alai)
! Direct shortwave that remains direct (umol/m[ground]2/s).
        absdrdr=absdr*(1.0D0-sigma)*kbl*exp(-sqrtexpr*kbl*alai)
! Shortwave penetrating shaded foliage (umol/m[foliage]2/s).
        abssh=absdf + (absdr-absdrdr) 
! Shortwave penetrating sunlit foliage (umol/m[foliage]2/s).
        abssl=abssh+(1.0D0-sigma)*kbL*fdir
        if(abssh.lt.0.0D0) abssh=0.0D0
        if(abssl.lt.0.0D0) abssl=0.0D0
        if(abssh.gt.1.0D0) abssh=1.0D0
        if(abssl.gt.0.0D0) abssl=1.0D0
        fracsl=exp(-kbl*alai)
      else
        abssh=0.001D0
        abssl=0.0D0
        fracsl=0
      endif

!      TRANS_SW = 1.d0 - ((1.d0-fracsl)*abssh + fracsl*abssl)
      TRANS_SW = ((1.d0-fracsl)*abssh + fracsl*abssl)

!      write(110,*) TRANS_SW,I0dr,I0df,abssh, abssl, sbeta, sigma, kbl


      end subroutine canopy_transmittance

!#############################################################################
!###################### INTEGRATION ##########################################

      subroutine qsimp(Xlim,crp,psp,ca,Gb,S,Sg,Sr)
!----------------------------------------------------------------------!
! qsimp calculates canopy photosynthesis by increasing the number of
! layers in the integration until the result (S) changes by less than
! 0.1 umol[CO2]/m2/s.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Passed parameters
      real*8,intent(in) :: Xlim
      type(canraddrv) :: crp
      type(psdrvtype) :: psp
      real*8,intent(in) :: ca
      real*8,intent(in) :: Gb
      !real*8,intent(in) :: cf, ci,T,P,RH
      real*8,intent(out) :: S,Sg,Sr
!----------------------------------------------------------------------!
!Local variables
!nu   real*8, parameter :: EPS=1.D-3
      integer, parameter :: MAXIT=6
      real*8,parameter :: ERRLIM=0.1d0
      integer IT, layers
      real*8 :: A, B, OST,OS,ST, ERROR
      real*8 :: Ac,Bc
      real*8 :: OSTg, OSg, STg !NK
      real*8 :: OSTr, OSr, STr !NK
!----------------------------------------------------------------------!
! Calculate canopy radiative transfer and photosynthesis in increasing 
!   number of layers.
! NOTE:  ERROR is only calculated on photosynthesis value (ST).
!        Corresponding conductance value is then returned also (STg).

      A=0.0D0
      B=Xlim      !Canopy LAI.  Ideally, qsimp should not see inside crp.
      Ac=0.0d0
      Bc=crp%LAI  !Cohort LAI, assumes spans height. Later:make different heights.
      OST=-1.D30
      OS= -1.D30
      OSTg=-1.D30 !NK
      OSg= -1.D30 !NK
      OSTr=-1.D30 !NK
      OSr= -1.D30 !NK
      layers=1
      do 11 IT=1,MAXIT
         CALL TRAPZD(A,B,Ac,Bc,ST,STg,STr,IT,layers,crp,psp,ca,Gb)
         S=(4.D0*ST-OST)/3.D0
         Sg=(4.D0*STg-OSTg)/3.D0 !NK
         Sr=(4.D0*STr-OSTr)/3.D0 !NK
         ERROR=ABS(S-OS)
         IF (ERROR.lt.ERRLIM) RETURN
         OS=S
         OST=ST
         OSg=Sg !NK
         OSTg=STg !NK
         OSr=Sr !NK
         OSTr=STr !NK
         if(IT.gt.1) layers=layers*2
   11 enddo
      return
      end subroutine qsimp

!======================================================================
!      subroutine trapzd(L1,L2,S,Sg,N,layers
!     $     ,crp,cf,ci,T,P,RH)
      subroutine trapzd(L1,L2,L1c,L2c,S,Sg,Sr,N,layers,crp,psp,ca,Gb)
!----------------------------------------------------------------------!
! Integrates canopy photosynthesis over canopy layers using Simpson's
! Rule (Press et al., 19??).
!----------------------------------------------------------------------!

      implicit none
!----------------------------------------------------------------------!
      real*8 L1,L2 !LAI of whole canopy
      real*8 L1c,L2c,S,Sg,Sr !LAI and sums of cohort
      integer N,layers
      type(canraddrv) :: crp 
      type(psdrvtype) :: psp
      real*8,intent(in) :: ca,Gb
      !-----Local------------------
      integer L
      real*8 DEL,X,SUM,SUMg,SUMr,RCL
      real*8 DELc
      real*8 Aleaf1 !Mean net photosynthesis at layer bottom (umol[CO2]/m2/s).
      real*8 Aleaf2 !Mean net photosynthesis at layer top (umol[CO2]/m2/s).
      real*8 gleaf1, gleaf2 !Mean conductance at L1, L2 (umol[H2O]/m2/s).
      real*8 Rdleaf1,Rdleaf2 !Mean leaf respiration at L1, L2 (umol[CO2]/m2/s)

      if(N.eq.1)then
         call photosynth_sunshd(L1,crp,psp,ca,Gb,Aleaf1,gleaf1,Rdleaf1)
         call photosynth_sunshd(L2,crp,psp,ca,Gb,Aleaf2,gleaf2,Rdleaf2)
!         S=0.5D0*(L2-L1)*(Aleaf1+Aleaf2)
!         Sg=0.5d0*(L2-L1)*(gleaf1+gleaf2)
!         Sr=0.5d0*(L2-L1)*(Rdleaf1+Rdleaf2)
         S=0.5D0*(L2c-L1c)*(Aleaf1+Aleaf2) !Cohort LAI for leaf processes.
         Sg=0.5d0*(L2c-L1c)*(gleaf1+gleaf2)
         Sr=0.5d0*(L2c-L1c)*(Rdleaf1+Rdleaf2)
      else
         RCL=layers  ! convert to real*8
         DEL=(L2-L1)/RCL !Canopy LAI for radiative transfer.
         DELc=(L2c-L1c)/RCL !Cohort LAI for radiative transfer.
         X=L1+0.5D0*DEL
         SUM=0.D0
         SUMg=0.d0
         SUMr=0.d0
         do 11 L=1,layers
           call photosynth_sunshd(X,crp,psp,ca,Gb,Aleaf1,gleaf1,Rdleaf1)
           SUM = SUM + Aleaf1
           SUMg = SUMg + gleaf1
           SUMr = SUMr + Rdleaf1
           X=X+DEL
   11    continue
!         S=0.5D0*(S+(L2c-L1c)*SUM/RCL)
!         Sg=0.5D0*(Sg+(L2c-L1c)*SUMg/RCL)
!         Sr=0.5D0*(Sr+(L2c-L1c)*SUMr/RCL)
         S=0.5D0*(S+(L2-L1)*SUM/RCL)
         Sg=0.5D0*(Sg+(L2-L1)*SUMg/RCL)
         Sr=0.5D0*(Sr+(L2-L1)*SUMr/RCL)
      endif
      return
      end subroutine trapzd


!============================================================================
      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0 (I think this is at least version 2.0)
!      USE CONSTANT, only : mrat,rvap,tf
      IMPLICIT NONE
!@var Physical constants from GISS GCM CONST.f
      real*8, parameter :: MWAT = 18.015d0 !molecular weight of water vapour (g/mol)
      real*8, parameter :: MAIR = 28.9655d0 !molecular weight of dry air (28.9655 g/mol)
      real*8, parameter :: MRAT = MWAT/MAIR 
      real*8, parameter :: RVAP = 1d3 * GASC/MWAT !gas constant for water vapour (461.5 J/K kg)
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
C**** LH = 0.5*(LH(0)+LH(t)), where LH(0)=
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(LH*(B-C/max(130.d0,TM)))/PR
      RETURN
      END FUNCTION QSAT
!============================================================================
!      FUNCTION QSATold (TM,QL,PR) Result(QSATcalc)
!      implicit none
!!@sum  QSAT calculates saturation vapour mixing ratio (kg/kg)
!!@auth Gary Russell
!!@ver  1.0
!!      USE CONSTANT, only : mrat,rvap,tf
!!      IMPLICIT NONE
!!@var A,B,C   expansion coefficients for QSAT
!      REAL*8, PARAMETER :: A=3.797915d0    !3.797915d0
!      REAL*8, PARAMETER :: B=7.93252d-6    !7.93252d-6
!      REAL*8, PARAMETER :: C=2.166847d-3         !2.166847d-3
!      real*8 :: TM, QL, PR
!      real*8 :: QSATcalc
!!**** Note that if QL is considered to be a function of temperature, the
!!**** correct argument in QSAT is the average QL from t=0 (C) to TM, ie.
!!**** QL = 0.5*(QL(0)+QL(t))
!!      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
!!      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
!!     REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap./sub. (J/kg)
!!      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
!      QSATcalc = A*EXP(QL*(B-C/TM))/PR
!
!    END function QSATold
!============================================================================

      end module biophysics !canopyspitters
