      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various components
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : IM,JM
      IMPLICIT NONE

!@var RUNOSI run off underneath sea/lake ice (kg/m^2) 
!@var ERUNOSI energy of run off underneath sea/lake ice (J/m^2) 
      REAL*8, DIMENSION(IM,JM) :: RUNOSI, ERUNOSI
!@var RUNOE run off from earth (kg/m^2) 
!@var ERUNOE energy of run off from earth (J/m^2) 
      REAL*8, DIMENSION(IM,JM) :: RUNOE, ERUNOE
!@var DMSI mass change of sea ice under ice and on open water (kg/m^2) 
!@var DHSI energy change of sea ice under ice and on open water (J/m^2) 
      REAL*8, DIMENSION(2,IM,JM) :: DMSI, DHSI
!@var RUNOLI run off from land ice (kg/m^2) (Energy always=0)
      REAL*8, DIMENSION(IM,JM) :: RUNOLI

C**** surface energy fluxes defined over type
!@param NSTYPE number of surface types for radiation purposes
      INTEGER, PARAMETER :: NSTYPE=4
!@var E0 net energy flux at surface for each type (J/m^2)
!@var E1 net energy flux at layer 1 for each type (J/m^2)
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: E0,E1
!@var EVAPOR evaporation over each type (kg/m^2) 
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: EVAPOR
!@var SOLAR solar radiation at surface over water and ice (J/m^2) 
      REAL*8, DIMENSION(2,IM,JM) :: SOLAR

!@var DMUA,DMVA momentum flux from atmosphere over each type (kg/m s) 
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: DMUA,DMVA
!@var DMUI,DMVI momentum flux from sea ice over each type (kg/m s) 
c      REAL*8, DIMENSION(IM,JM) :: DMUI,DMVI

C**** currently saved - should be replaced by fluxed quantities
!@var DU1,DV1 momentum flux from atmosphere summed over type (m/s)
      REAL*8, DIMENSION(IM,JM) :: DU1,DV1
!@var DTH1,DQ1 heat/water flux from atmos. summed over type (C, kg/kg)
      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1

!@var FLOWO,EFLOWO runoff and energy of runoff into ocean (kg, J)
      REAL*8, DIMENSION(IM,JM) :: FLOWO,EFLOWO

!@var PREC precipitation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, DIMENSION(IM,JM) :: EPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PRECSS

!@var GTEMP ground temperature (upper two levels) over surface type (C)
      REAL*8, DIMENSION(2,NSTYPE,IM,JM) :: GTEMP

      END MODULE FLUXES

