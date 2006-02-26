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
      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: zero = 0d0

      !************************************************************************
      !* PHYSICAL CONSTANTS
      real*8,parameter :: stbo =5.67051d-8 !Stefan-Boltzman (W/m2/K4)
      real*8,parameter :: lhe = 2.5d6 !Latent heat evap (2.5008d6 J/kg)
      real*8,parameter :: rhow = 1d3  !Density of pure water (1000 kg/m^3)
      real*8,parameter :: tfrz = 273.16d0 !freezing pt of H2O at 1 atm (Kelvin)
      real*8,parameter :: gasc = 8.314510d0 !gas constant (8.314510 J/mol K)
      real*8,parameter :: Avogadro=6.023d23 !Avogadro's constant (atmos/mole)
      real*8,parameter :: cp=1012. !Heat capacity of dry air (J kg-1 K-1)

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
      
      end module ent_const
