#include "rundeck_opts.h"

      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various components
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm,ntsurfsrcmax,nt3Dsrcmax
#endif
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
C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice (J/m^2)
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
      REAL*8, DIMENSION(2,IM,JM) :: DMSI, DHSI, DSSI
!@var fmsi_io,fhsi_io,fssi_io basal ice-ocean fluxes (kg or J/m^2)
      REAL*8, DIMENSION(IM,JM) :: fmsi_io,fhsi_io,fssi_io
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
!@var UI2rho Ustar*2*rho ice-ocean friction velocity on atmospheric grid
      REAL*8, DIMENSION(IM,JM) :: UI2rho

C**** currently saved - should be replaced by fluxed quantities
!@var DTH1,DQ1 heat/water flux from atmos. summed over type (C, kg/kg)
      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1

!@var uflux1 surface turbulent u-flux (=-<uw>) 
!@var vflux1 surface turbulent v-flux (=-<vw>)
!@var tflux1 surface turbulent t-flux (=-<tw>)
!@var qflux1 surface turbulent q-flux (=-<qw>)
      real*8, dimension(im,jm) :: uflux1,vflux1,tflux1,qflux1

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
!@var SSS sea surface salinity on atmospheric grid (ppt)
      REAL*8, DIMENSION(IM,JM) :: SSS
!@var MLHC ocean mixed layer heat capacity (J/m^2 C) 
      REAL*8, DIMENSION(IM,JM) :: MLHC

#ifdef TRACERS_ON
!@var TRSOURCE non-interactive surface sources/sinks for tracers (kg/s)
      REAL*8, DIMENSION(IM,JM,ntsurfsrcmax,NTM):: trsource
!@var TRSRFFLX interactive surface sources/sinks for tracers (kg/s)
      REAL*8, DIMENSION(IM,JM,NTM):: trsrfflx
!@var TRFLUX1 total surface flux for each tracer (kg/s)
      REAL*8, DIMENSION(IM,JM,NTM):: trflux1
!@var TRGRDEP gravitationally settled tracers at surface (kg)
      REAL*8, DIMENSION(NTM,IM,JM):: TRGRDEP
!@var GTRACER ground concentration of tracer on atmospheric grid (kg/kg)
      REAL*8, DIMENSION(NTM,NSTYPE,IM,JM):: GTRACER
!@var TR3DSOURCE 3D sources/sinks for tracers (kg/s)
      REAL*8, DIMENSION(IM,JM,LM,nt3Dsrcmax,NTM):: tr3Dsource

#ifdef TRACERS_WATER
!@var TRPREC tracers in precip (kg)
      REAL*8, DIMENSION(NTM,IM,JM):: TRPREC
!@var TREVAPOR tracer evaporation over each type (kg/m^2) 
      REAL*8, DIMENSION(NTM,NSTYPE,IM,JM) :: TREVAPOR
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg)
!@var TRUNOE tracer runoff from earth (kg)
!@var TRUNOLI tracer runoff from land ice (kg)
      REAL*8, DIMENSION(NTM,IM,JM):: TRUNPSI, TRUNOSI, TRUNOE, TRUNOLI
!@var TRFLOWO tracer in river runoff into ocean (kg)
      REAL*8, DIMENSION(NTM,IM,JM) :: TRFLOWO
!@var DTRSI tracer flux in sea ice under ice and on open water (kg)
      REAL*8, DIMENSION(NTM,2,IM,JM) :: DTRSI
!@var ftrsi_io ice-ocean tracer fluxes under ice (kg)
      REAL*8, DIMENSION(NTM,IM,JM) :: ftrsi_io
#endif
#endif

      END MODULE FLUXES

