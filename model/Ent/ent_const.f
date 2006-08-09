      module ent_const
!@sum  CONSTANT definitions for physical constants
!@sum  These are constants that would be common to the GCM/EWB, so should
!@sum  have the GCM/EWB constants substituted in for coupled runs.
!@auth N.Kiang
!@ver  1.0

      !* COUPLED RUNS - Replace with values from GCM constants*!
      !use Name_of_GCM_constants_module  

      implicit none
      save

      !************************************************************************
      !* INTEGRATION - TIME STEPS
      integer,parameter :: T_SUB=12  !Number of sub-time steps in a year
      !************************************************************************
      !* NUMERICAL CONSTANTS
      real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: zero = 0d0

      !************************************************************************
      !* PHYSICAL CONSTANTS
      real*8,parameter :: stbo =5.67051d-8 !Stefan-Boltzman (W/m2/K4)
      real*8,parameter :: lhe = 2.5d6 !Latent heat evap (2.5008d6 J/kg)
      real*8,parameter :: rhow = 1d3  !Density of pure water (1000 kg/m^3)
      real*8,parameter :: tfrz = 273.16d0 !freezing pt of H2O at 1 atm (Kelvin)
      real*8,parameter :: KELVIN = tfrz
      real*8,parameter :: gasc = 8.314510d0 !gas constant (8.314510 J/mol K)
      real*8,parameter :: Avogadro=6.023d23 !Avogadro's constant (atmos/mole)
      real*8,parameter :: cp=1012. !Heat capacity of dry air (J kg-1 K-1)
!@param shw heat capacity of water (at 20 C) (4185 J/kg C)
      real*8,parameter :: shw  = 4185.

      !************************************************************************
      !* CONSTANTS FROM GISS CONST.f
      !@param mwat molecular weight of water vapour
      real*8,parameter :: mwat = 18.015d0
      !@param mair molecular weight of dry air (28.9655 g/mol)
      real*8,parameter :: mair = 28.9655d0
      !@param mrat  mass ratio of air to water vapour (0.62197)
      real*8,parameter :: mrat = mwat/mair    ! = 0.62197....
      !@param rvap  gas constant for water vapour (461.5 J/K kg)
      !**** defined as R/M_W = 1000* 8.314510 J/mol K /18.015 g/mol
      real*8,parameter :: rvap = 1d3 * gasc / mwat ! = 461.5...
      !@param tf freezing point of water at 1 atm (273.16 K)
      real*8,parameter :: tf = tfrz

      !************************************************************************
       !* ASTRONOMICAL CONSTANTS

      !************************************************************************
      !* RUN CONTROL
      real*8,parameter :: undef=-1.d30 ! Missing value
      real*8,parameter :: teeny=1.d-30 ! Small positive value to avoid 0/0
      real*8, parameter :: EPS = 1.d-8  !Small error
      real*8, parameter :: EPS2 = 1.d-12 !Smaller error



      
      !************************************************************************
      !* Ent CONSTANTS
      !******************
      !* PATCH DYNAMICS *
      !******************
      real, parameter ::  F_AREA=.01 !* min area of patch as fraction of total
      real, parameter ::  BTOL=.00001!* min cohort biomass for termination (kgC/m2)
      real, parameter ::  NTOL=.001  !* min plant density for whatever
      integer,parameter :: PATCH_DYNAMICS = 0 ! 0-No, 1=Yes

      !************************************************************************
       !* ASTRONOMICAL CONSTANTS
      real*8,parameter :: sday = 86400.! sec per day (s)
      real*8,parameter :: SECPY = 31536000  ! sec per year (s)

      !************************************************************************
      !********************
      !* SOIL / HYDROLOGY *
      !********************
!      integer :: N_DEPTH        !Number of soil layers.  SET IN ENT_INIT
      !!! setting it to constant for the time being to simplify the compilation
      integer, parameter :: N_DEPTH = 6

      !**********************
      !* RADIATIVE TRANSFER *
      !**********************
      integer,parameter :: N_BANDS = 6 !Number of spectral bands (GISS 6)
                                !Expect to adjust to hyperspectral

      !***********************
      !* ECOLOGICAL DYNAMICS *
      !***********************
      integer,parameter :: N_DIST_TYPES = 2 !Number of disturbance types


      !************************************************************************
      !*  GISS VEGETATION CONSTANTS
      !integer :: JEQUATOR !Will get calculated in ent_init.

      integer,parameter :: N_PFT = 8
      integer,parameter :: N_SOILCOV = 2 !light sand, dark dirt (GISS)
      integer,parameter :: N_OTHER = 2 ! algae, c4 grass
      integer,parameter :: N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER

      !************************************************************************
      !* COHORT Biomass pools
      integer,parameter :: N_BPOOLS = 6 !foliage,sapwood,hardwood,labile, fine root,coarse root
      integer,parameter :: FOL = 1   !FOLIAGE Array indices for array(N_BPOOLS)
      integer,parameter :: SW = 2    !SAPWOOD
      integer,parameter :: HW = 3    !HARDWOOD
      integer,parameter :: LABILE = 4 !LABILE
      integer,parameter :: FR = 5    !FINE ROOT
      integer,parameter :: CR = 6    !COARSE ROOT

      !************************************************************************
      !* CASA SOIL CONSTANTS
      !* THESE SHOULD GO IN THE CASA MODULE BUT NEED TO BE USED IN ENT_TYPES.F
      !* Array sizes
      integer,parameter :: N_SOIL_TYPES = 5
      integer,parameter :: NLIVE = 3
      integer,parameter :: NDEAD = 9
      integer,parameter :: NPOOLS = NLIVE + NDEAD
      !* Total pool array indices
      integer,parameter :: PTRACE = 2 !Trace elements in Tpools, C and N
      integer,parameter :: Carbon = 1
      integer,parameter :: Nitrogen = 2
!      integer,parameter :: ptrace = 2  !num. nutrient pools used in CASA resp. routine -PK
      integer,parameter :: nresp_pools = 14  !num. pools used in CASA resp. routine -PK
      real*8,parameter :: Q10 = 2.        !Q10 used in belowground calculations --> value from lit -PK 5/25/06

      !* Live pool array indices
      integer,parameter :: LEAF = 1  !Array index
      integer,parameter :: FROOT = 2  !Array index
      integer,parameter :: WOOD = 3  !Array index

      !* Dead pool array indices
      integer,parameter :: SURFMET = 4 !sfc metabolic
      integer,parameter :: SURFSTR = 5 !sfc structural
      integer,parameter :: SOILMET = 6 !soil metabolic
      integer,parameter :: SOILSTR = 7 !soil structural
      integer,parameter :: CWD = 8     !coarse woody debris
      integer,parameter :: SURFMIC = 9 !sfc microbial
      integer,parameter :: SOILMIC = 10 !soil microbial
      integer,parameter :: SLOW = 11 !slowly decomposing soil o.m. pool (up to a decade)
      integer,parameter :: PASSIVE = 12 !very slowly decomposing soil o.m. pool (decades-centuries)

      !************************************************************************
      real*8 :: CNratio(NPOOLS)
            data CNratio/
     1            30.0,       ! C:N ratio of leaf pool
     2           130.0,       ! C:N ratio of wood pool
     3            55.0,       ! C:N ratio of froot pool
     4            30.0,       ! C:N ratio of surfmet pool
     5            50.0,       ! C:N ratio of surfstr pool
     6            25.0,       ! C:N ratio of soilmet pool
     7            50.0,       ! C:N ratio of soilstr pool
     8           135.0,       ! C:N ratio of cwd pool
     9            12.5,       ! C:N ratio of surfmic pool
     a            12.5,       ! C:N ratio of soilmic pool
     b            12.5,       ! C:N ratio of slow pool
     c             8.5/       ! C:N ratio of passive pool

      real*8,dimension(N_PFT,NPOOLS) :: annK !CASA turnover times (sec-1)
      !real*8,dimension(N_PFT,NPOOLS) :: kdt !CASA turnover times for wood & dead (yr-1)
      real*8,dimension(N_PFT) :: solubfract !Soluble ("metabolic") fraction of litter 
      real*8,dimension(N_PFT) :: structurallignin !fraction of structural C from lignin -PK 7/5/06 
      real*8,dimension(N_PFT) :: lignineffect !effect of lignin on decomp -PK 7/5/06
      real*8,parameter :: woodligninfract = 0.40 !amt lignin in wood C -PK 7/5/06
      end module ent_const
