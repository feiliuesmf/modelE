module TracerConstants_mod
  use CONSTANT, only: molMassWater => mwat
  implicit none
  private

  ! Usage:
  !    mm_CO2 = TRACER_CONSTANTS%CO2%molMass

  public :: TRACER_CONSTANTS
  public :: H2O18

  type Properties
    real*8 :: molMass    ! g
    real*8 :: solubility = 0 ! units?
  end type Properties

  type (Properties), parameter ::  &
       H2O18 = Properties(molMass = 20.d0)

  type TracerMetadata
    type (Properties) :: CO2
    type (Properties) :: O3
    type (Properties) :: water
    type (Properties) :: H2O18
    type (Properties) :: HDO
    type (Properties) :: H2O17
  end type TracerMetadata

  
  type (TracerMetadata), parameter :: TRACER_CONSTANTS = TracerMetadata( &
       & CO2   = Properties(molMass = 44.d0),    &
       & O3    = Properties(molMass = 48.d0),    &
       & water = Properties(molMass = 18.015d0), &
       & H2O18 = H2O18, &
       & HDO   = Properties(molMass = 19d0),     &
       & H2O17 = Properties(molMass = 19d0)      &
       & )

end module TracerConstants_mod

