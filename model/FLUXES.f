      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various components
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,JM
      IMPLICIT NONE

!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface (J/m^2)
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: RUNOSI, ERUNOSI, SRUNOSI
!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: RUNPSI, SRUNPSI
!@var RUNOE run off from earth (kg/m^2)
!@var ERUNOE energy of run off from earth (J/m^2)
      REAL*8, DIMENSION(IM,JM) :: RUNOE, ERUNOE
!@var DMSI mass flux of sea ice under ice and on open water (kg/m^2)
!@var DHSI energy flux of sea ice under ice and on open water (J/m^2)
!@var DSSI salt flux in sea ice under ice and on open water (kg/m^2)
      REAL*8, DIMENSION(2,IM,JM) :: DMSI, DHSI, DSSI
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
!@var SOLAR absorbed solar radiation (J/m^2)
!@+   SOLAR(1)  absorbed by open water
!@+   SOLAR(2)  absorbed by ice 
!@+   SOLAR(3)  absorbed by water under the ice
      REAL*8, DIMENSION(3,IM,JM) :: SOLAR

C**** Momemtum stresses are calculated as if they were over whole box
!@var DMUA,DMVA momentum flux from atmosphere for each type (kg/m s) 
      REAL*8, DIMENSION(IM,JM,NSTYPE) :: DMUA,DMVA
!@var DMUI,DMVI momentum flux from sea ice to ocean (kg/m s)
      REAL*8, DIMENSION(IM,JM) :: DMUI,DMVI

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

