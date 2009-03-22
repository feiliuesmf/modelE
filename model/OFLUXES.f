#include "rundeck_opts.h"

      MODULE OFLUXES
!@sum  OFLUXES contains the fluxes between various components
!       that are going to be used on ocean grid
!@auth Larissa Nazarenko
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE OCEANRES,  only : imo,jmo
!      USE DOMAIN_DECOMP_1D, ONLY : grid
      USE OCEANR_DIM, only : grid=>ogrid

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      IMPLICIT NONE

!@var SOLAR absorbed solar radiation (J/m^2)
!@+   SOLAR(1)  absorbed by open water
!@+   SOLAR(2)  absorbed by ice
!@+   SOLAR(3)  absorbed by water under the ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oSOLAR
!@param NSTYPE number of surface types for radiation purposes
      INTEGER, PARAMETER :: NSTYPE=4
!@var E0 net energy flux at surface for each type (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oE0
!@var EVAPOR evaporation over each type (kg/m^2) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oEVAPOR
C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice (J/m^2)
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oDMSI, oDHSI, oDSSI
!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface (J/m^2)
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oRUNOSI, oERUNOSI, oSRUNOSI
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg, J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oFLOWO, oEFLOWO
!@var APRESS total atmos + sea ice pressure (at base of sea ice) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oAPRESS
!@var MELTI,EMELTI,SMELTI mass,energy,salt from simelt into ocean (kg,J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oMELTI, oEMELTI, oSMELTI

!@var DMUA,DMVA momentum flux from atmosphere for each type (kg/m s) 
!@+   On OCN A grid (tracer point)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oDMUA, oDMVA
!@var DMUI,DMVI momentum flux from sea ice to ocean (kg/m s)
!@+   On OCN C grid 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oDMUI, oDMVI

!@var GMELT,EGMELT mass,energy from glacial melt into ocean (kg,J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oGMELT, oEGMELT

!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var ERUNPSI energy of run off from sea/lake ice after precip (J/m^2)
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oRUNPSI, oERUNPSI, oSRUNPSI   
!@var PREC precipitation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oPREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oEPREC

!@var RSI fraction of open water area covered in ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oRSI

#ifdef TRACERS_OCEAN

#ifdef TRACERS_GASEXCH_ocean
!@var TRGASEX  tracer gas exchange over each type (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oTRGASEX
#endif

#ifdef TRACERS_WATER
!@var TRFLOWO tracer in river runoff into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oTRFLOWO
!@var TREVAPOR tracer evaporation over each type (kg/m^2) 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oTREVAPOR
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg/m^2)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: oTRUNOSI, oTRUNPSI
!@var TRMELTI tracer from simelt into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oTRMELTI
!@var TRPREC tracers in precip (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: oTRPREC

!@var TRGMELT tracer from glacial melt into ocean (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oTRGMELT
#endif
!@var DTRSI tracer flux in sea ice under ice and on open water (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oDTRSI

#ifdef TRACERS_DRYDEP
!@var TRDRYDEP tracer dry deposition by type (kg/m^2) (positive down)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: oTRDRYDEP 
#endif
#endif

      END MODULE OFLUXES

      SUBROUTINE ALLOC_OFLUXES(grd_dum)
!@sum   Initializes FLUXES''s arrays
!@auth  Larissa Nazarenko
!@ver  1.0
      USE CONSTANT, only : tf
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE OFLUXES
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      !I-J arrays
      ALLOCATE( oRUNOSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          oERUNOSI ( I_0H:I_1H , J_0H:J_1H ), 
     &          oSRUNOSI ( I_0H:I_1H , J_0H:J_1H ),
     &          oRUNPSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          oSRUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          oERUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          oDMUI    ( I_0H:I_1H , J_0H:J_1H ),
     &          oDMVI    ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER )
      ALLOCATE( oFLOWO   ( I_0H:I_1H , J_0H:J_1H ),
     &          oEFLOWO  ( I_0H:I_1H , J_0H:J_1H ),
     &          oMELTI   ( I_0H:I_1H , J_0H:J_1H ),
     &          oEMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          oSMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          oGMELT   ( I_0H:I_1H , J_0H:J_1H ),
     &          oEGMELT  ( I_0H:I_1H , J_0H:J_1H ),
     &          oPREC    ( I_0H:I_1H , J_0H:J_1H ),
     &          oEPREC   ( I_0H:I_1H , J_0H:J_1H ),
     &          oAPRESS  ( I_0H:I_1H , J_0H:J_1H ),
     &          oRSI     ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER)


       !n-I-J arrays
       ALLOCATE( oDMSI    (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           oDHSI    (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           oDSSI    (  2  , I_0H:I_1H , J_0H:J_1H ),
     &           oSOLAR   (  3  , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)


      !I-J-: arrays
      ALLOCATE( oE0      ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &          oEVAPOR  ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &          oDMUA    ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &          oDMVA    ( I_0H:I_1H , J_0H:J_1H , 1 ),
     &   STAT = IER)


!TRACERS **********************
#ifdef TRACERS_OCEAN

      !:-:-I-J arrays
#ifdef TRACERS_GASEXCH_ocean

      ALLOCATE( oTRGASEX( NTM , 1 , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      oTRGASEX=0.     !initialize to zero

#endif

#ifdef TRACERS_WATER
      !(:)-(:)-I-J arrays
      ALLOCATE( oTREVAPOR( NTM , 1, I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)


       !:-I-J arrays
       ALLOCATE( oTRPREC  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRUNPSI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRUNOSI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRFLOWO ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           oTRMELTI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)

       ALLOCATE( oTRGMELT ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif

#ifdef TRACERS_DRYDEP
       ALLOCATE(oTRDRYDEP( NTM , 1 , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
         oTRDRYDEP = 0.   !Initialize to 0.
#endif

      ALLOCATE( oDTRSI( NTM ,    2   , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif


      END SUBROUTINE ALLOC_OFLUXES
