      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various components
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : IM,JM
      IMPLICIT NONE

!@var RUNOSI run off underneath sea/lake ice (kg/m^2) 
!@var ERUNOSI energy of run off underneath sea/lake ice (J/m^2) 
      REAL*8, DIMENSION(IM,JM) :: RUNOSI, ERUNOSI
!@var DMSI mass change at base of sea ice (kg/m^2) 
!@var DHSI energy change at base of sea ice (J/m^2) 
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
!@var SOLAR solar radiation at surface over each type (J/m^2) 
c      REAL*8, DIMENSION(IM,JM,NSTYPE) :: SOLAR
!@var SENSHT sensible heat at surface over each type (J/m^2) 
c      REAL*8, DIMENSION(IM,JM,NSTYPE) :: SENSHT

!@var DMUA,DMVA momentum flux from atmosphere over each type (?) 
c      REAL*8, DIMENSION(IM,JM,NSTYPE) :: DMUA,DMVA
!@var DMUI,DMVI momentum flux from sea ice over each type (?) 
c      REAL*8, DIMENSION(IM,JM) :: DMUI,DMVI

C**** currently saved - should be replaced by fluxed quantities
!@var DU1,DV1 momentum flux from atmosphere over each type (?) 
      REAL*8, DIMENSION(IM,JM) :: DU1,DV1
!@var DTH1,DQ1 momentum flux from sea ice over each type (?) 
      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1

!@var FLOWO,EFLOWO runoff and energy of runoff into ocean
      REAL*8, DIMENSION(IM,JM) :: FLOWO,EFLOWO

!@var PREC precipitation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, DIMENSION(IM,JM) :: EPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PRECSS

      END MODULE FLUXES

