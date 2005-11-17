      module ent_types
        use ent_const
        implicit none


!        contains
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
         integer :: pft           ! PFT number
         real*8 :: n              ! Number of individuals in cohort
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

         !* ALL QUANTITIES BELOW ARE FOR AN INDIVIDUAL *!

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              ! Height (m)
         real*8 :: crown_dx       ! Crown horizontal axis length (m)
         real*8 :: crown_dy       ! Crown vertical axis length (m)
         real*8 :: dbh            ! Stem diameter at breast height (m)
         real*8 :: root_d         ! Root half spheroid diameter (m)
         real*8 :: LAI
         real*8 :: clump          ! Leaf clumping parameter (TBA)

         !* BIOMASS POOLS
         real*8 :: LMA            ! Leaf mass per leaf area (kgC/m2-leaf)
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
         real*8 :: GPP            ! GPP flux/plant over time step
         real*8 :: NPP            ! NPP flux/plant over time step
         real*8 :: R_growth       ! Growth respiration/plant over time step
         real*8 :: R_maint        ! Maintenance respir/plant over time step
         real*8 :: N_up           ! N uptake from soil over time step
         real*8 :: C_litter       ! C in litterfall
         real*8 :: N_litter       ! N in litterfall
         real*8 :: C_to_Nfix      ! Carbon flux to N fixers symbionts

         !* REPRODUCTION
         !real*8 :: ------        ! ASK PAUL
         
      end type cohort

!****************************************************************************
      type canradtype
         !Arrays in height levels in the canopy
         real*8,pointer :: heights !height levels (m) betwn layers (#layers+1)
         real*8,pointer :: GORTclump(:)
         real*8,pointer :: LAI(:)
         real*8,pointer :: Isw(:) 
         real*8,pointer :: IPAR(:)
         real*8,pointer :: Ibeam(:)
         real*8,pointer :: Idiff(:)
      end type canradtype

!****************************************************************************
      type patch
         real*8 :: age                !Patch age (years)
         real*8 :: area               !Patch area (units?)
         type(entcelltype),pointer:: cellptr !Pointer to grid cell
         type(patch),pointer :: older !Pointer to next older patch
         type(patch),pointer :: younger !Pointer to next younger patch
         type(cohort),pointer :: tallest !Pointer to tallest cohort
         type(cohort),pointer :: shortest !Pointer to shortest cohort
         
         !* Flux variables for GCM/EWB - patch total
         real*8 :: albedo(N_BANDS) !Spectral albedo, average over patch
         real*8 :: z0              !Roughness length, average over patch
         real*8 :: GCANOPY         !Canopy conductance of water vapor (mm/s)
         real*8 :: CO2flux         !Net CO2 flux up (umol-CO2/m2-gnd/s)

         !* Variables calculated by GCM/EWB - downscaled from grid cell
         real*8 :: Soilmoist(N_DEPTH) !Available soil moisture by depth (mm)
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
         real*8 :: LAI_p
         real*8 :: plant_ag_Cp !Total above-gnd biomass carbon (kgC/m2)
         real*8 :: plant_bg_Cp !Total below-gnd biomass carbon (kgC/m2)
         real*8 :: plant_ag_Np !Total above-gnd nitrogen (kgN/m2)
         real*8 :: plant_bg_Np !Total below-gnd nitrogen (kgN/m2)
         !real*8 :: etc. foliage, sapwood, labile, storage, etc.

         !* Biomass pools - by pft
         real*8 :: LAI(N_PFT)
         real*8 :: plant_ag_C(N_PFT) !Total above-gnd biomass carbon (kgC/m2)
         real*8 :: plant_bg_C(N_PFT) !Total below-gnd biomass carbon (kgC/m2)
         real*8 :: plant_ag_N(N_PFT) !Total above-gnd nitrogen (kgN/m2)
         real*8 :: plant_bg_N(N_PFT) !Total below-bnd nitrogen (kgN/m2)
         !real*8 :: etc. foliage, sapwood, labile, storage, etc.

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

         !* Activity diagnostics - can be summed by month, year, etc.
         real*8 :: GPP               !Gross primary productivity (kgC/m2/day)
         real*8 :: NPP               !Net primary productivity (kgC/m2/day)
         real*8 :: Soil_resp         !Soil heterotrophic respiration (kgC/m2/day)
      end type patch

!****************************************************************************
!****************************************************************************
      type entcelltype
         real*8 :: long, lat      !longitude, latitude
         integer :: longi, latj    !grid cell i,j
         real*8 :: area         !Area km^2
         type(patch), pointer:: youngest
         type(patch), pointer:: oldest

         !Cell-level summary values - PHYSICAL
         !EXPORT - from radiative transfer
         real*8 :: albedo(N_BANDS) !Albedo may be in bands or hyperspectral
         real*8 :: TRANS_SW

         !SOIL - IMPORT
         real*8 :: soil_Phi      !Soil porosity (m3/m3)
         real*8 :: soildepth    !Soil depth (m)
         real*8 :: theta_max    !Saturated soil water volume (m/m)
         real*8 :: k_sat        !Saturated hydraulic conductivity
         real*8 :: root_Phi     !Infiltration factor promoted by roots (units?)

         !VEGETATION - EXPORT STATE
         real*8 :: z0           !Roughness length (m)
         real*8 :: GCANOPY      !Canopy conductance of water vapor (mm s-1)
         real*8 :: CO2flux      !CO2 flux (umol m-2 s-1)
         real*8 :: GPP          !GPP
         !real*8 :: NPP          !NPP
         !real*8 :: VOCflux     !Other kind of fluxes, aerosols from fire, etc.
         !Cell-level diagnostic values - BIOLOGICAL
         !e.g. LAI, biomass pools, nitrogen pools, PFT fractions, GDD, GPP, etc
         real*8 :: LAI
         real*8 :: froot        !Fraction of roots in soil layer
         !-----
         

         !VEGETATION - PRIVATE
!         real*8 :: Ci           !Internal foliage CO2 (mol/m3) !!Cohort level
         real*8 :: Qf           !Foliage surface vapor mixing ratio (kg/kg)

         !METEOROLOGICAL - IMPORT STATE VARIABLES
         !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
         real*8 :: TcanopyC     !Canopy temperatue (Celsius)
         real*8 :: Qv           !Canopy air specif humidity (kg vapor/ kg air)
         real*8 :: P_mbar       !Atmospheric pressure (mb)
         real*8 :: Ca           !@Atmos CO2 conc at surface height (mol/m3).
         real*8 :: Soilmoist(N_DEPTH) !May be an array by depth (units TBA)
         real*8 :: fice         !Fraction of soil layer that is ice
         real*8 :: betad  !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8 :: Precip       !Precipitation (mm)
         real*8 :: Ch           !Ground to surface heat transfer coefficient 
         real*8 :: U            !Surface layer wind speed (m s-1)

         !Radiation - IMPORT STATE VARIABLE
         !may later be broken down into hyperspectral increments.
         ! in an array
         real*8 :: Isw          !Incident shortwave 100-2000 nm (W m-2)
         real*8 :: IPAR         !Incident PAR 400-700 nm (W m-2)
         real*8 :: Ibeam        !Incident beam radiation (W m-2)
         real*8 :: Idiff        !Incident diffuse radiation (W m-2),
         real*8 :: Solarzen     !Solar zenith angle
         real*8 :: fdir         !Fraction of surface vis rad that is direct 

      end type entcelltype


!****************************************************************************
      type entdatatype
        real longmin, longmax, latmin, latmax
        integer longi,latj
        type(timestruct) :: tt      !Greenwich Mean Time
        type (entcelltype), pointer :: grid(:,:)
      end type entdatatype


!****************************************************************************
!*    CONSTANTS
!****************************************************************************


      !************************************************************************
       !* ASTRONOMICAL CONSTANTS
      real*8,parameter :: sday = 86400.! sec per day (s)

      !************************************************************************
      !* RUN CONTROL
      integer,parameter :: PATCH_DYNAMICS = 0 ! 0-No, 1=Yes

      !************************************************************************
      !* Ent CONSTANTS
      !***********
      !* BIOLOGY *
      !***********
      integer,parameter :: N_PFT = 13 !Number of plant functional types
      integer,parameter :: N_OTHER = 0
      !* 1 - evergreen broadleaf early successional
      !* 2 - evergreen broadleaf late successional
      !* 3 - evergreen needleleaf early successional
      !* 4 - evergreen needleleaf late successional
      !* 5 - cold deciduous broadleaf early successional
      !* 6 - cold deciduous broadleaf late successional
      !* 7 - drought deciduous broadleaf
      !* 8 - cold adapted shrub
      !* 9 - arid adapted shrub
      !* 10- C3 grass
      !* 11- C4 grass
      !* 12- arctic C3 grass
      !* 13- C4 crops

      integer,parameter :: PST_PFT(N_PFT) !Photosynthetic pathway for pft
     &     = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2/)       !1=C3, 2=C4

      real,parameter :: SSTAR_PFT(N_PFT)!Soil moisture stress point for pft
     &     = (/0.6, 0.6, 0.55, 0.55, 0.5, 0.5, 0.45, 0.5, 0.4, 
     &     0.65, 0.55, 0.6, 0.65/) !Guesses, except for pfts 7 and 10

      real,parameter :: NF_PFT(N_PFT) !Kull&Kruijt ps cap/leaf N param
     &     = (/1.3, 1.2, 0.9, 0.85, 1.5, 1.4, 1.4, 1.4, 1.3, 
     &         0.5, 0.76, 1.4, 0.76/)  ! Guesses from Friend&Kiang (2005)

      integer,parameter :: N_DIST_TYPES = 2 !Number of disturbance types

      integer,parameter :: N_BANDS = 3 !Number of spectral bands (VIS,NIR,MIR)
                                !Expect to adjust to hyperspectral

      !*****************
      !* SOIL / HYDROLOGY
      !*****************
      integer,parameter :: N_DEPTH = 6  !Number of soil layers.  Eventually
                                        !this should be taken from hydrology

      !******************
      !* PATCH DYNAMICS *
      !******************
      real, parameter ::  F_AREA=.01 !* min area of patch as fraction of total
      real, parameter ::  BTOL=.00001!* min cohort biomass for termination (kgC/m2)
      real, parameter ::  NTOL=.001  !* min plant density for whatever
      !real etc -------------



      !***************************************************
      !Temp values for Ent pfts (See ent_const.f for types)
      type(pftype),parameter :: pfpar(N_PFT) =         & !PFT parameters
     !     pst,  hwilt,    sstar, swilt,  nf !
     &           (/                                     &       
     &     pftype(1,   -100.d0,  .60d0, .29d0,  1.3d0), &
     &     pftype(1,   -100.d0,  .60d0, .29d0,  1.2d0), &
     &     pftype(1,   -100.d0,  .55d0, .25d0,  0.9d0), &
     &     pftype(1,   -100.d0,  .55d0, .25d0,  0.85d0),&
     &     pftype(1,   -100.d0,  .50d0, .29d0,  1.5d0), &
     &     pftype(1,   -100.d0,  .50d0, .29d0,  1.4d0), &
     &     pftype(1,   -100.d0,  .45d0, .22d0,  1.4d0), &
     &     pftype(1,   -100.d0,  .50d0, .22d0,  1.4d0), &
     &     pftype(1,   -100.d0,  .40d0, .22d0,  1.3d0), &
     &     pftype(1,   -100.d0,  .65d0, .27d0,  1.5d0), &
     &     pftype(2,   -100.d0,  .55d0, .22d0,  0.76d0),&
     &     pftype(1,   -100.d0,  .60d0, .27d0,  1.4d0), &
     &     pftype(2,   -100.d0,  .65d0, .27d0,  0.76d0) &
     &     /)

      !***************************************************
      !*            GISS VEGETATION TYPES                *
      !***************************************************
      !Types: These are GISS GCM types until we get full data sets for Ent.
      !       In addition, boreal forest from BOREAS SOBS is added.
      !1-tundra, 2-grassland, 3-shrubland, 4-savanna, 5-deciduous forest,
      !6-evergreen needleleaf forest, 7-tropical rainforest, 8-crops
      !9-boreal forest
      integer,parameter :: N_PFT = 8
      integer,parameter :: N_SOILCOV = 2 !light sand, dark dirt (GISS)
      integer,parameter :: N_OTHER = 1
      integer,parameter :: N_COVERTYPES = 11

!      type(pftype),parameter :: pfpar(N_PFT) =         & !PFT parameters
!     !     pst,  hwilt,    sstar, swilt,  nf !
!     &           (/                                     &       
!     &     pftype(1,   -100.d0,  .50d0, .30d0,  1.4d0), &
!     &     pftype(2,   -100.d0,  .45d0, .27d0,  1.5d0),&
!     &     pftype(2,   -100.d0,  .65d0, .22d0,  1.3d0),&
!     &     pftype(2,   -100.d0,  .65d0, .22d0,  1.3d0),&
!     &     pftype(1,   -100.d0,  .55d0, .29d0,  1.5d0),&
!     &     pftype(1,   -100.d0,  .60d0, .25d0,  0.9d0),&
!     &     pftype(1,   -100.d0,  .55d0, .26d0,  1.1d0),&
!     &     pftype(2,   -100.d0,  .45d0, .27d0,  1.3d0),&
!     &     pftype(1,   -100.d0,  .50d0, .30d0,  0.76d0)&
!     &     /)

      !**NOTE:  Above, sstar and swilt are guesses for all except grassland
      ! and savanna.



      !***************************************************
      !*         END - GISS VEGETATION TYPES             *
      !***************************************************

      end module ent_types
