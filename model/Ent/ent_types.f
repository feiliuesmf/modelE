      module ent_types
      use ent_const

      implicit none

! the following two parameters are basically needed for integrity 
! checks (to avoid infinite loops when processing linked lists)
!@var MAX_PATCHES maximal number of patches per cell
!@var MAX_COHORTS maximal number of cohorts per patch
        integer, parameter :: MAX_PATCHES=32, MAX_COHORTS=64


!****************************************************************************
!*       TYPE DECLARATIONS
!****************************************************************************
      type timestruct
         !Structure for holding time variables
         integer :: year
         integer :: month
         integer :: jday  !day of year
         integer :: day   !day of month, not of year
         integer :: hour  !0-23
         integer :: minute 
         integer :: seconds
      end type timestruct
!****************************************************************************
      type pftype
         !Constant parameters specific to vegetation type
         integer :: pst ! Photosynth type 1=C3, 2=C4
         real*8 :: hwilt  !Wilting point matric potential (m)
         real*8 :: sstar  !Rel. soil moist at stress onset (Rodriguez-Iturbe)
         real*8 :: swilt  !Normalized soil water at wilting point (dim'less)
         real*8 :: nf !Canopy nitrogen factor (dimensionless) (Kull and Kruijt)
         !CASA parameters
         real*8 :: sla !Specific leaf area (m^2 leaf area/kg C)
         real*8 :: lrage !CASA Turnover time of leaves and roots (years)
         real*8 :: woodage !CASA Turnover time of stems (years)
         real*8 :: lit_C2N !CASA litcn_casa (C:N ratio) IS THIS FOLIAGE&ROOTS?
         real*8 :: lignin  !CASA lignin (UNITS?  lignin content of ??)
         !real*8 :: !Phenological parameters, other
      end type pftype
!****************************************************************************

!      type :: veg_par_type  !This type is variables specific to biophysics.f
!         !DYNAMICALLY VARYING VEGETATION VARIABLES
!         !from veg_set_cell: VEGETATION SPECIFICATIONS
!         !VEGETATION GEOMETRY BIOLOGY
!         !@var alai  LAI over whole grid cell
!         real*8 :: alai
!         !@var nm   Mean canopy nitrogen (g/m2[leaf]) over whole grid cell
!         real*8 :: nm
!         !@var vh   Mean canopy height (m) over whole grid cell
!         real*8 :: vh
!         !@var Ntot Total canopy nitrogen (g/m[ground]2).
!         real*8 Ntot
         !@var alait  Array of LAI by pft, alait*vfraction=grid LAI contrib.
         !real*8, dimension(8) :: alait 
         !@var vfraction  Array of cover fraction by pft in grid cell
         !real*8, dimension(8) :: vfraction 

         !CANOPY RADIATIVE TRANSFER
         !@var sigma Leaf scattering coefficient (?unitless).
!         real*8 :: sigma        !=0.2D0
!         real*8 :: sqrtexpr
!         !@var kdf Canopy extinction coeff. for diffuse radiation (unitless).
!         real*8 :: kdf          !=0.71D0
!         !@var rhor Canopy reflectivity (?unitless).
!         real*8 :: rhor
!         !@var kbl Canopy extinction coeff. for direct radiation (unitless).
!         real*8 :: kbl
!         !@var leafalbedo Leaf albedo, variable with physiological status 
!         real*8 :: leafalbedo 
!         !@var vegalbedo  Vegetation canopy albedo averaged over the grid cell.
!         !UPDATE WITH CANOPY RADIATIVE TRANSFER SCHEME
!         real*8 :: vegalbedo  
!      end type veg_par_type

!****************************************************************************
      type cohort
         integer :: pft           !* PFT number
         real*8 :: n              ! Density of individuals in cohort (#/m^2)
         type(entcelltype),pointer :: cellptr !Pointer to ent grid cell
         type(patch),pointer :: pptr    !Pointer to patch
         type(cohort),pointer :: taller !Pointer to next tallest cohort
         type(cohort),pointer :: shorter !Pointer to next shortest cohort
         type(cohort),pointer :: csptaller !Pointer to next taller conspecific
         type(cohort),pointer :: cspshorter !Pointer to next shorter conspecfic

         !* PFT PARAMETERS
         ! Only need to index array of pftypes.
         
         !* NITROGEN status */
         !@var nm   Mean cohort nitrogen (g/m2[leaf]) over whole grid cell
         real*8 :: nm
         !@var Ntot Total cohort nitrogen (g/m[ground]2).
         real*8 Ntot
         !@var LAI Total cohort leaf area index (m2[leaf]/m2[ground])
         real*8 LAI               !*
         !@var LMA Leaf mass per leaf area (gC/m2)
         real*8 LMA

         !* ALL QUANTITIES BELOW ARE FOR AN INDIVIDUAL *!

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* Height (m)
         real*8 :: crown_dx       ! Crown horizontal axis length (m)
         real*8 :: crown_dy       ! Crown vertical axis length (m)
         real*8 :: dbh            ! Stem diameter at breast height (m)
         real*8 :: root_d         ! Root half spheroid diameter (m)
         !real*8 :: LA            ! Leaf area (m2[leaf])
         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8,pointer :: froot(:) ! Fraction of roots in soil layer

         !* BIOMASS POOLS (SHOULD CONVERT TO PER INDIVIDUAL LATER)
         real*8 :: C_fol          ! Foliage carbon (=LMA*LAI = kgC/m2-gnd)
         real*8 :: N_fol          ! Foliage nitrogen (gN/m2-gnd)
         real*8 :: C_sw           ! Sapwood carbon (=rho_wood*sw_vol) (units?)
         real*8 :: N_sw           ! Sapwood nitrogen
         real*8 :: C_hw           ! Dead stem (heartwood) carbon
         real*8 :: N_hw           ! Dead stem (heartwood) nitrogen
         real*8 :: C_lab          ! Labile stored carbon 
         real*8 :: N_lab          ! Labile stored nitrogen
         real*8 :: C_froot        ! Fine root carbon
         real*8 :: N_froot        ! Fine root nitrogen
         real*8 :: C_croot        ! Coarse root carbon
         real*8 :: N_croot        ! Coarse root nitrogen

         !* FLUXES
         real*8 :: gcanopy        ! Conductance of water vapor/plant (kg/s)
         real*8 :: GPP            ! GPP flux/plant over time step (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux/plant over time step (kg-C/m2/s)
         real*8 :: R_growth       ! Growth respiration/plant over time step (kg-C/m2/s)
         real*8 :: R_maint        ! Maintenance respir/plant over time step (kg-C/m2/s)
         real*8 :: N_up           ! N uptake from soil over time step (kg-N/m2/s)
!         real*8 :: C_litter       ! C in litterfall
!         real*8 :: N_litter       ! N in litterfall
         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts

         !* REPRODUCTION
         !real*8 :: ------        ! ASK PAUL
         
      end type cohort

!****************************************************************************
      type canradtype
         !Arrays in height levels in the canopy
         real*8,pointer :: heights !height levels (m) betwn layers (#layers+1)
         real*8,pointer :: LAI(:) !LAI within height level
         !Whole-canopy foliage clumping factor
         real*8 :: GORTclump
      end type canradtype

!****************************************************************************
      type patch
         real*8 :: age                !*Patch age (years)
         real*8 :: area               !*Patch area (units?)
         type(entcelltype),pointer:: cellptr !Pointer to grid cell
         type(patch),pointer :: older !Pointer to next older patch
         type(patch),pointer :: younger !Pointer to next younger patch
         type(cohort),pointer :: tallest !Pointer to tallest cohort
         type(cohort),pointer :: shortest !Pointer to shortest cohort
         type(cohort),pointer :: sumcohort !Summary of properties of cohorts.
         !NOTE:  In sumcohort:  
         ! * Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals (may want to do by biomass of cohort)
         ! * Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* Structural variables *!
         real*8 :: nm   !Mean canopy nitrogen (g/m2[leaf]) over patch

         !* Flux variables for GCM/EWB - patch total
         real*8 :: z0              !Roughness length, average over patch
         real*8 :: albedo(N_BANDS) !Spectral albedo, average over patch
         real*8 :: TRANS_SW     !Transmittance of shortwave radiation to the ground (fraction)
         real*8 :: GCANOPY         !Canopy conductance of water vapor (mm/s)
         real*8 :: CO2flux         !Net CO2 flux up (kg-C/m2-gnd/s)
         real*8 :: Ci           !*Internal foliage CO2 (mol/m3) !!Cohort level
         real*8 :: betad  !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) !Water stress in layers.
         real*8 :: Soil_resp       !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8 :: R_can  !Canopy respiration. Temporary hack (kg-C/m2/s) -NYK 7/28/06 ## Will go in sumcohort.
         real*8 :: R_root !Root respiration. Temporary hack (kg-C/m2/s) -NYK 7/28/06 ## Will go in sumcohort.

         !* Variables calculated by GCM/EWB - downscaled from grid cell
!         real*8,pointer :: Soilmoist(:) !Available soil moisture by depth (mm)
         real*8 :: Soilmoist !Soil moisture (volumetric fraction, avg top 30 cm) -PK 6/28/06
         real*8 :: N_deposit          !N deposition (kgN/m2)

         !* Variables for biophysics and biogeochemistry
         type(canradtype) :: crad !Data structure for light profile

         !* Disturbance values
         real*8 :: fuel
         real*8 :: ignition_rate
         real*8 :: lambda1(T_SUB)  !Site-averaged fire dist. rate during year
         real*8 :: disturbance_rate(N_DIST_TYPES)

         !* DIAGNOSTIC SUMMARIES
         !* Biomass pools - patch total
         !* SEE avgcohort and sumcohort
         real*8 :: LAI(N_COVERTYPES)
         real*8 :: C_froot           !Carbon in fine roots.

         !* Soil data (needed for albedo computation)
         integer soil_type      ! 1 - sand (bright) ; 2 - dirt (dark)

         !* CASA carbon pools *!
         !"live" & "dead"
         real*8, DIMENSION(PTRACE,NPOOLS) :: Tpool !(g-C/m^2, CASA Tpools, single cell)

#ifdef NEWDIAG
         !* Soil pools - patch total
         real*8 :: REW               !Relative extractable water (REW)
         real*8 :: soil_labile_C     !Labile soil carbon (kgC/m2)
         real*8 :: soil_slow_C       !Slow soil carbon pool (kgC/m2)
         real*8 :: soil_labile_N     !Labile soil nitrogen (kgN/m2)
         real*8 :: soil_slow_N       !Slow soil nitrogen (kgN/m2)
         real*8 :: mineral_N         !Mineralized N (kgN/m2)
         !real*8 :: other soil pools TBA

         !* Rates of change - patch total
         real*8 :: dadt              !Rate of change of patch age = 1
         real*8 :: dpdt              !Rate of change of patch area
         real*8 :: dwdt              !Rate of change of available soil water
#endif
         !* Activity diagnostics - can be summed by month, year, etc.
         real*8 :: GPP               !Gross primary productivity (kgC/m2/s)
         real*8 :: NPP               !Net primary productivity (kgC/m2/s)
!         real*8 :: Soil_resp         !Soil heterotrophic respiration (kgC/m2/s) !moved to fluxes -PK 6/14/06
      end type patch

!****************************************************************************
!****************************************************************************
      type entcelltype
 !        real*8 :: long, lat      !longitude, latitude
 !        integer :: longi, latj    !grid cell i,j
         real*8 :: area         !Area km^2
         type(patch), pointer:: youngest
         type(patch), pointer:: oldest
         type(patch), pointer:: sumpatch !summary of patches in entcell

         !Cell-level summary values - PHYSICAL
         !EXPORT - from radiative transfer
         real*8 :: albedo(N_BANDS) !Albedo may be in bands or hyperspectral
         real*8 :: TRANS_SW  !Transmittance of shortwave radiation to the ground (fraction)

         !SOIL - IMPORT
         real*8 :: soil_Phi      !Soil porosity (m3/m3)
         real*8 :: soildepth    !Soil depth (m)
         real*8 :: theta_max    !Saturated soil water volume (m/m)
         real*8 :: k_sat        !Saturated hydraulic conductivity
         real*8 :: root_Phi     !Infiltration factor promoted by roots (units?)

         !SOIL - CONSTANTS
         !Soil textures for CASA -PK
         !use siltfrac = 1-(clayfrac+sandfrac)?? (R&A scheme has loam/peat!) -PK 
         real*8 :: soil_texture(N_SOIL_TYPES) ! fractions of soil textures
                                              ! in upper 30 cm of soil
!         real*8 clayfrac  !fractional clay content (passed from GHY.f)
!         real*8 sandfrac  !fractional sand content (also from GHY.f)

         !VEGETATION - EXPORT STATE
         real*8 :: z0           !Roughness length (m)
!         real*8 :: GCANOPY      !Canopy conductance of water vapor (mm s-1)
!         real*8 :: CO2flux      !CO2 flux (umol m-2 s-1)
!         real*8 :: GPP          !GPP
!         !real*8 :: NPP          !NPP
!         !real*8 :: VOCflux     !Other kind of fluxes, aerosols from fire, etc.
         !Cell-level diagnostic values - BIOLOGICAL
         !e.g. LAI, biomass pools, nitrogen pools, PFT fractions, GDD, GPP, etc
!         real*8 :: LAI(N_COVERTYPES)
         real*8 :: LAI  !Total LAI
!         real*8,pointer :: froot(:) !Fraction of roots in soil layer
!         real*8 :: C_froot      !Carbon in fine roots
!         real*8 :: betad  !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
!         real*8,pointer :: betadl(:) !Water stress in layers.
         !-----

         !VEGETATION - PUBLIC
!         real*8 :: Ci           !*Internal foliage CO2 (mol/m3) !!Cohort level

         !METEOROLOGICAL - IMPORT STATE VARIABLES
         !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
         real*8 :: TcanopyC     !Canopy temperatue (Celsius)
         real*8 :: Qf           !*Foliage surface vapor mixing ratio (kg/kg)
         real*8 :: P_mbar       !Atmospheric pressure (mb)
         real*8 :: Ca           !@Atmos CO2 conc at surface height (mol/m3).
          !CASA needs next 2 only for top 30 cm (top 2 = 27 cm; plus 3/dzsoi of lyr 3) -PK
         real*8 :: Soilmoist !Soil moisture (volumetric fraction)
         real*8 :: Soiltemp  !Soil temperature (Celsius)
         real*8,pointer :: Soilmp(:) !Soil matric potential (m)
         real*8,pointer :: fice(:) !Fraction of soil layer that is ice
         real*8 :: Ch           !Ground to surface heat transfer coefficient 
         real*8 :: U            !Surface layer wind speed (m s-1)


         !Radiation - IMPORT STATE VARIABLES
         !may later be broken down into hyperspectral increments.
         ! in an array
!         real*8 :: Ivis          !Incident visible  (W m-2)
!         real*8 :: Idir          !Incident direct visible  (W m-2)
!         real*8 :: IPAR         !Incident PAR 400-700 nm (W m-2)
         real*8 :: IPARdir        !Incident direct PAR (W m-2)
         real*8 :: IPARdif        !Incident diffuse PAR (W m-2)
         real*8 :: Solarzen     !Solar zenith angle

      end type entcelltype


!****************************************************************************
!      type entdatatype
!        real longmin, longmax, latmin, latmax
!        integer longi,latj
!        type(timestruct),pointer :: tt      !Greenwich Mean Time
!        type (entcelltype),pointer :: grid(:,:)
!      end type entdatatype


      end module ent_types
