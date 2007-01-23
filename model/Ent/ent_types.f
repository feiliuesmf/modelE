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
         !CLM respiration parameters
         real*8 :: r      !CLM respiration parameter (gC/gN)
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
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* ALL QUANTITIES BELOW ARE FOR AN INDIVIDUAL *!

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h  !<<<            !* Height (m)
         real*8 :: crown_dx       ! Crown horizontal axis length (m)
         real*8 :: crown_dy       ! Crown vertical axis length (m)
         real*8 :: dbh   !<<<         ! Stem diameter at breast height (m)
         real*8 :: root_d         ! Root half spheroid diameter (m)
         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8,pointer :: fracroot(:) ! Fraction of roots in soil layer

         !* BIOMASS POOLS (g-C/single plant)
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

         !* FLUXES (for whole cohort over area cover)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) !!Cohort level
         real*8 :: gcanopy        ! Conductance of water vapor/cohort (m/s)
         real*8 :: GPP            ! GPP flux/cohort/area cover (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux/cohort/area cover (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration/cohort/area cover (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: N_up           ! N uptake from soil/cohort/area cover (kg-N/m2/s)
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
         real*8 :: area               !*Patch area (fraction of entcell)
         type(entcelltype),pointer:: cellptr !Pointer to grid cell
         type(patch),pointer :: older !Pointer to next older patch
         type(patch),pointer :: younger !Pointer to next younger patch
         type(cohort),pointer :: tallest !Pointer to tallest cohort
         type(cohort),pointer :: shortest !Pointer to shortest cohort

         !*- - - - - - - Cohorts summary variables - - - - - - - - - - -*!
         !  Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals or LAI.
         !  Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* DIAGNOSTICS - NITROGEN and LEAF status */
         !@var nm   Mean leaf nitrogen (g/m2[leaf])
         real*8 :: nm
         !@var Ntot Total cohort nitrogen (g/m[ground]2).
         real*8 Ntot
         !@var LMA Leaf mass per leaf area (gC/m2)
         real*8 LMA
         !@var LAI Total cohort leaf area index (m2[leaf]/m2[ground])
         real*8 LAI               !*
         !@var LAIpft LAI by cover type.
         real*8,pointer :: LAIpft(:) !(N_COVERTYPES)
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* mean Height (m)
         real*8 :: crown_dx       ! Crown horizontal axis length (m)
         real*8 :: crown_dy       ! Crown vertical axis length (m)
!         real*8 :: dbh   !<<<         ! Stem diameter at breast height (m)
!         real*8 :: root_d         ! Root half spheroid diameter (m)
         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8,pointer :: fracroot(:) ! Fraction of roots in soil layer

         !* DIAGNOSTICS - BIOMASS POOLS (g-C/cohort)
         real*8 :: C_fol          ! Foliage carbon (=LMA*LAI = kgC/m2-gnd)
         real*8 :: N_fol          ! Foliage nitrogen (gN/m2-gnd)
         real*8 :: C_w           ! Sapwood+hardwood carbon (=rho_wood*sw_vol) (units?)
         real*8 :: N_w           ! Sapwood+hardwood nitrogen
         real*8 :: C_lab          ! Labile stored carbon 
         real*8 :: N_lab          ! Labile stored nitrogen
         real*8 :: C_froot        ! Fine root carbon
         real*8 :: N_froot        ! Fine root nitrogen
         real*8 :: C_root        ! Fine+coarse root carbon
         real*8 :: N_root        ! Fine+coarse root nitrogen

         !* EXPORT - FLUXES (whole patch)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) 
         real*8 :: GCANOPY        ! Conductance of water vapor (m/s)
         !* DIAGNOSTICS - FLUXES
         real*8 :: GPP            ! GPP flux (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: N_up           ! N uptake from soil(kg-N/m2/s)
!         real*8 :: C_litter       ! C in litterfall
!         real*8 :: N_litter       ! N in litterfall
!         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts
         !- - - - - - - end of cohort summary variables - - - - - - - - - - - - -

         !- - - - - - - Patch total - - - - - - - - - - - - - - - - - - - - - - -

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED
         real*8 :: z0              !Roughness length, average over patch
         !* EXPORT 
         real*8 :: albedo(N_BANDS)   !Spectral albedo, average over patch
         real*8 :: betad             !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) !Water stress in layers.
         real*8 :: TRANS_SW          !Transmittance of shortwave radiation to the ground (fraction)
         real*8 :: CO2flux           !Net CO2 flux up (kg-C/m2-gnd/s)
         !* DIAGNOSTICS - soil
         real*8 :: Soil_resp         !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8, DIMENSION(PTRACE,NPOOLS) :: Tpool !<<< !(g-C/m^2, CASA Tpools, single cell)

         !* IMPORT - Variables calculated by GCM/EWB - downscaled from grid cell
!         real*8,pointer :: Soilmoist(:) !Available soil moisture by depth (mm)
         real*8 :: Soilmoist    !Soil moisture (volumetric fraction, avg top 30 cm) -PK 6/28/06
!         real*8 :: N_deposit    !N deposition (kgN/m2)

         !* Variables for biophysics and biogeochemistry
         type(canradtype) :: crad !Data structure for light profile

         !* Disturbance values
         real*8 :: fuel
         real*8 :: ignition_rate
         real*8 :: lambda1(T_SUB) !Site-averaged fire dist. rate during year
         real*8 :: disturbance_rate(N_DIST_TYPES)


         !* Soil data (needed for albedo computation)
         integer soil_type      ! 1 - sand (bright) ; 2 - dirt (dark)

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
      end type patch


!****************************************************************************
      type entcelltype
 !        real*8 :: long, lat      !longitude, latitude
 !        integer :: longi, latj    !grid cell i,j
         real*8 :: area         !Area km^2
         type(patch), pointer:: youngest
         type(patch), pointer:: oldest

         !*- - - - - - - Cohorts summary variables - - - - - - - - - - -*!
         !* Per vegetated ground area of entcell ** excludes bare soil area.
         !  Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals or leaf area.
         !  Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED - NITROGEN and LEAF status */
         !@var nm   Mean leaf nitrogen (g/m2[leaf]).
         real*8 :: nm
         !@var Ntot Total cohort nitrogen (g/m[ground]2) for vegetated patches only.
         real*8 Ntot
         !@var LMA Leaf mass per leaf area (gC/m2)
         real*8 LMA
         !@var LAI Leaf area index (m2[leaf]/m2[ground]) for vegetated patches only.
         real*8 LAI               
         !@var LAIpft LAI by cover type.
         real*8,pointer :: LAIpft(:) !(N_COVERTYPES)
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED - GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* mean Height (m)
!         real*8 :: crown_dx       ! Crown horizontal axis length (m)
!         real*8 :: crown_dy       ! Crown vertical axis length (m)
!         real*8 :: dbh   !<<<         ! Stem diameter at breast height (m)
!         real*8 :: root_d         ! Root half spheroid diameter (m)
!         real*8 :: clump          ! Leaf clumping parameter (TBA)
         real*8,pointer :: fracroot(:) ! Fraction of roots in soil layer

         !*  IMPORT-PRESCRIBED, EXPORT-SIMULATED - BIOMASS POOLS (g-C/cohort)
         real*8 :: C_fol          ! Foliage carbon (=LMA*LAI = kgC/m2-gnd)
         real*8 :: N_fol          ! Foliage nitrogen (gN/m2-gnd)
         real*8 :: C_w           ! Sapwood+dead wood carbon (=rho_wood*sw_vol) (units?)
         real*8 :: N_w           ! Sapwood+dead wood nitrogen
         real*8 :: C_lab          ! Labile stored carbon 
         real*8 :: N_lab          ! Labile stored nitrogen
         real*8 :: C_froot        ! Fine root carbon
         real*8 :: N_froot        ! Fine root nitrogen
         real*8 :: C_root        ! Fine+coarse root carbon
         real*8 :: N_root        ! Fine+coarse root nitrogen

         !* IMPORT/EXPORT PUBLIC - FLUXES)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) 
         real*8 :: GCANOPY        ! Conductance of water vapor (m/s)
         !* EXPORT - FLUXES 
         real*8 :: GPP            ! GPP flux (kg-C/m2/s)
         real*8 :: NPP            ! NPP flux (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: N_up           ! N uptake from soil(kg-N/m2/s)
!         real*8 :: C_litter       ! C in litterfall
!         real*8 :: N_litter       ! N in litterfall
!         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts
         !- - - - - - - end of cohort summary variables - - - - - - - - - - - - - - - - - - -

         !- - - - - -  Patch-level summary values - PHYSICAL ------------------
         !* EXPORT - from radiative transfer
         real*8 :: z0              !Roughness length, average over patch
         real*8 :: albedo(N_BANDS) !Albedo may be in bands or hyperspectral
         real*8 :: betad             !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) !Water stress in layers.
         real*8 :: TRANS_SW  !Transmittance of shortwave radiation to the ground (fraction)
         real*8 :: CO2flux           !Net CO2 flux up (kg-C/m2-gnd/s)
         !* DIAGNOSTICS - soil
         real*8 :: Soil_resp         !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8, DIMENSION(PTRACE,NPOOLS) :: Tpool !<<< !(g-C/m^2, CASA Tpools, single cell)

         !* Disturbance values
         real*8 :: fuel
         real*8 :: ignition_rate
         real*8 :: lambda1(T_SUB) !Site-averaged fire dist. rate during year
         real*8 :: disturbance_rate(N_DIST_TYPES)

         !- - - - - - Entcell-level variables - - - - - - - - - - - - - - - -
         !* IMPORT/EXPORT - vegetation
         real*8 :: fv  ! vegetation fraction of entcell
         real*8 :: heat_capacity ! total veg. heat capacity
         !* IMPORT - SOIL
         real*8 :: soil_Phi      !Soil porosity (m3/m3)
         real*8 :: soildepth    !Soil depth (m)
         real*8 :: theta_max    !Saturated soil water volume (m/m)
         real*8 :: k_sat        !Saturated hydraulic conductivity
         real*8 :: root_Phi     !Infiltration factor promoted by roots (units?)

         !SOIL - CONSTANTS
         !Soil textures for CASA -PK
         !use siltfrac = 1-(clayfrac+sandfrac)?? (R&A scheme has loam/peat!) -PK 
         real*8 :: soil_texture(N_SOIL_TYPES) ! fractions of soil textures, upper 30 cm of soil
!         real*8 clayfrac  !fractional clay content (passed from GHY.f)
!         real*8 sandfrac  !fractional sand content (also from GHY.f)


         !IMPORT - METEOROLOGICAL STATE VARIABLES
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
         real*8 :: CosZen     !Solar zenith angle

      end type entcelltype


!****************************************************************************
!      type entdatatype
!        real longmin, longmax, latmin, latmax
!        integer longi,latj
!        type(timestruct),pointer :: tt      !Greenwich Mean Time
!        type (entcelltype),pointer :: grid(:,:)
!      end type entdatatype

!****************************************************************************
      type ent_config
      ! this type should contain all parameters that describe the run
      ! i.e. flags, array dimensions etc. They assumed to be constant
      ! during the run but may change from run to run
        logical do_soilresp       ! do soil respiration

      end type ent_config


      end module ent_types
