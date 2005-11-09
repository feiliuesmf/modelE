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
       !* ASTRONOMICAL CONSTANTS

      !************************************************************************
      !* RUN CONTROL
      real*8,parameter :: undef=-1.d30 ! Missing value
      real*8,parameter :: teeny=1.d-30 ! Small positive value to avoid 0/0
      real*8, parameter :: EPS = 1.d-8  !Small error
      real*8, parameter :: EPS2 = 1.d-12 !Smaller error
      


      end module ent_const
