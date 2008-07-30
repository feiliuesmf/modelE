#include "rundeck_opts.h"

      MODULE SEAICE
!@sum  SEAICE contains all the sea ice related subroutines
!@auth Original Development Team
!@ver  1.0
!@cont PREC_SI,SEA_ICE,ADDICE,SIMELT
      USE CONSTANT, only : lhm,rhoi,byrhoi,rhow,shi,shw,byshi,bylhm,sday
     *     ,rhows
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@param LMI number of temperature layers in ice
      INTEGER, PARAMETER :: LMI = 4
!@param XSI fractions of mass layer in each temp. layer
!@param BYXSI recipricol of XSI
      REAL*8, PARAMETER, DIMENSION(LMI) ::
     *     XSI= (/0.5d0, 0.5d0, 0.5d0, 0.5d0/),
     *     BYXSI= (/ 1./XSI(1), 1./XSI(2), 1./XSI(3), 1./XSI(4) /)
!@param Z1I thickness of first layer ice (m)
      REAL*8, PARAMETER :: Z1I = .1d0
!@param ACE1I ice mass first layer (kg/m^2)
      REAL*8, PARAMETER :: ACE1I = Z1I*RHOI
!@param HC1I heat capacity of first layer ice (J/m^2)
      REAL*8, PARAMETER :: HC1I = ACE1I*SHI
!@param Z2OIM min. thickness of 2nd layer ice (m) (if > 0)
      REAL*8, PARAMETER :: Z2OIM = .1d0    ! .4d0
!@param AC2OIM min. ice mass 2nd layer (kg/m^2)   (if > 0)
      REAL*8, PARAMETER :: AC2OIM = Z2OIM*RHOI
!@param ALAMI,ALAMS lambda coefficient for ice/snow J/(m*degC*sec)
      REAL*8, PARAMETER :: ALAMI=2.1762d0, ALAMS=0.35d0
!@param alamdS salinity/temp coefficient for conductivity (J/(m*s)/psu)
      real*8, parameter :: alamdS=0.117d0
!@param RHOS density of snow (kg/m^3)
      REAL*8, PARAMETER :: RHOS = 300.0
!@var FLEADOC lead fraction for ocean ice (%)
      REAL*8, PARAMETER :: FLEADOC = 0.06d0
!@var FLEADLK lead fraction for lakes (%)
      REAL*8, PARAMETER :: FLEADLK = 0.
!@param BYRLI,BYRLS reciprical of density*lambda
      REAL*8, PARAMETER :: BYRLI = 1./(RHOI*ALAMI),
     *     BYRLS = 1./(RHOS*ALAMS)
!@param MU coefficient of seawater freezing point w.r.t. salinity 
      REAL*8, PARAMETER :: MU = 0.054d0   ! C/ppt
!@param SSI0 default value for sea ice salinity (kg/kg) (=3.2ppt)
      REAL*8, PARAMETER :: SSI0 = 0.0032d0
!@param FSSS fraction of ocean salinity found in new-formed ice
      REAL*8, PARAMETER :: FSSS = 13d0/35d0
!@var qsfix flag is true if salinity of sea ice is constant
      LOGICAL :: QSFIX = .false.
!@var alpha implicity for heat diffusion in sea ice (1=fully implicit)
      REAL*8, PARAMETER :: ALPHA = 1.0
!@dbparam oi_ustar0 default ice-ocean friction velocity (m/s)
      REAL*8 :: oi_ustar0 = 1d-3  ! 5d-3 ! not used if ice dynamics is
!@dbparam silmfac factor controlling lateral melt of ocean ice
      REAL*8 :: silmfac = 1.4d-8 ! = pi*(3d-6)/0.66/1000
!@var silmpow exponent for temperature dependence of lateral melt
      REAL*8 :: silmpow = 1.36d0
!@dbparam snow_ice =1 to allow for snow ice formation (default=1)
      INTEGER :: snow_ice = 1
!@var osurf_tilt controls calc. of ocean surface tilt for ice dyn:
!@+       from geostrophy (=0) or from free surface (=1, default)
      INTEGER :: osurf_tilt = 1
!@var DEBUG flag
      LOGICAL DEBUG
!@var FORM formulation of sea ice thermodynamics (BP or SI)
      CHARACTER*2 :: FORM = "SI"  ! "BP"

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SRUN0
     *     ,ERUN0,
#ifdef TRACERS_WATER
     &     TRSIL,TRPRCP,TRUN0,
#endif
     &     WETSNOW)
!@sum  PREC_SI Adds the precipitation to sea/lake ice
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@param SNOMAX maximum allowed snowdepth (1m equivalent) (kg/m^2)
!@param dSNdRN rate of conversion of snow to ice as a function of rain
      REAL*8, PARAMETER :: SNOMAX=1d0*RHOS, dSNdRN=0.
!@var ENRGP total energy of precip (C),(J/m^2)
!@var PRCP amount of precip (kg/m^2)
      REAL*8, INTENT(IN) :: ENRGP, PRCP
!@var SNOW snow mass (kg/m^2)
!@var MSI1 first layer ice mass (= SNOW + ACE1I) (kg/m^2)
!@var MSI2 second layer ice mass (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var TSIL temperature of ice layers (C) (DIAGNOSTIC ONLY)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var WETSNOW true if snow is wet (i.e. has been rained on)
      LOGICAL, INTENT(OUT) :: WETSNOW
!@var RUN0 runoff from ice (kg/m^2)
!@var SRUN0 salt in runoff from ice (kg/m^2)
!@var ERUN0 energy in runoff from ice (J/m^2)
      REAL*8, INTENT(OUT) :: RUN0, SRUN0, ERUN0
      REAL*8 :: FMSI2, FMSI3, FHSI2, FHSI3,
     *     CMPRS, SNWF, RAIN, FREZI, MSI1, FSSI2, SNOW1, 
     *     FSSI3, MELTI, MELTS, SMELTI, SICE, HCMPRS, HSNOW,
     *     HICE, HMELTI, HSNOW2
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
      REAL*8, DIMENSION(NTM) :: FTRSI2,FTRSI3,TRRMF
     *     ,TRMELTI,TRMELTS,TRCMPRS,TRICE,TRSNOW
#endif

      WETSNOW=.FALSE.
      RUN0=0. ; SRUN0=0. ; ERUN0=0.
#ifdef TRACERS_WATER
      TRUN0=0.
#endif
      MSI1=SNOW+ACE1I

      IF (PRCP.gt.0) THEN

C**** Snowfall is calculated from precip energy (0 deg or colder)
      SNWF = MAX(0d0,MIN(PRCP,-ENRGP*BYLHM))
C**** Rain is remaining precip (0 deg or warmer)
      RAIN = PRCP-SNWF
      WETSNOW = RAIN.GT.1d-5*PRCP  ! i.e. a noticeable fraction of prec

C**** New formulation: separate out snow and ice components
C**** Calculate remaining snow and necessary mass flux using
C**** SNOW =SNOW  + SNWF  - CMPRS - MELTS
C**** ACE1I=ACE1I + FREZI + CMPRS - FMSI2
C**** SALTI=SALTI - SMELTI        - FSSI2
      IF (SNOW.gt.0) THEN ! apply fluxes to snow portion first
        SNOW1 = MIN(SNOW,XSI(1)*MSI1)
        HSNOW = SNOW1*(Ti(HSIL(1)/(XSI(1)*MSI1),
     *       1d3*SSIL(1)/(XSI(1)*MSI1))*shi-lhm)
        HICE  = HSIL(1)-HSNOW
        SICE  = SSIL(1)
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSIL(:,1)*MIN(SNOW1/(XSI(1)*MSI1-SSIL(1)),1d0)
        TRICE(:)  = TRSIL(:,1)-TRSNOW(:)
#endif
        MELTS=MAX(0d0,(HSNOW+ENRGP)*BYLHM+SNOW1+SNWF)
        IF (MELTS.gt.SNOW1+SNWF) THEN   ! all snow and some ice melts
          MELTS = SNOW1+SNWF
          MELTI = Mi(HICE+HSNOW+ENRGP,SICE,XSI(1)*MSI1-SNOW1) 
          SMELTI = MELTI*SICE/(XSI(1)*MSI1-SNOW1)
          FREZI = 0.
          HICE = HICE+HSNOW+ENRGP
          HSNOW = 0.
          SNOW1 = 0.
          CMPRS = 0.
#ifdef TRACERS_WATER
          TRRMF(:)  = TRPRCP(:)*(1.-(SNWF+FREZI)/PRCP)
          TRMELTS(:)= TRSNOW(:) + TRPRCP(:)*SNWF/PRCP
          TRMELTI(:)= (MELTI-SMELTI)*TRICE(:)/(XSI(1)*MSI1-SNOW1-SICE)
          TRCMPRS(:)= 0.
          TRSNOW(:) = 0.
          TRICE(:)  = TRICE(:) - TRMELTI(:)
#endif
c         ACE1I = ACE1I - MELTI
        ELSE ! some snow remains
          MELTI =0.
          SMELTI=0.
          FREZI=MIN(RAIN,MAX(-(HSNOW+ENRGP)*BYLHM-SNOW1-SNWF+MELTS,0d0))
          HICE  = HICE - LHM*FREZI
          HSNOW = HSNOW + ENRGP + LHM*FREZI
          IF (SNWF.ge.PRCP .and. SNOW+SNWF .GT. SNOMAX) THEN
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
            CMPRS = SNOW+SNWF-0.9d0*SNOMAX
          ELSE
C**** RAIN and MELT COMPRESSES SNOW INTO ICE
            CMPRS = MIN(dSNdRN*(RAIN+MELTS), SNOW+SNWF-MELTS)
          END IF
          HCMPRS = HSNOW*CMPRS/(SNOW1+SNWF-MELTS)
          HSNOW = HSNOW - HCMPRS
          HICE  = HICE  + HCMPRS
#ifdef TRACERS_WATER
          TRRMF(:)  = TRPRCP(:)*(1.-(SNWF+FREZI)/PRCP)
          TRSNOW(:) = TRSNOW(:)+ TRPRCP(:)*SNWF/PRCP
          TRMELTS(:)= MELTS*TRSNOW(:)/(SNOW1+SNWF)
          TRCMPRS(:)= CMPRS*(TRSNOW(:)-TRMELTS(:))/(SNOW1+SNWF-MELTS)
          TRMELTI(:)= 0.
          TRICE(:)  = TRICE(:) + TRPRCP(:)*FREZI/PRCP + TRCMPRS(:)
          TRSNOW(:) = TRSNOW(:) - TRMELTS(:) - TRCMPRS(:)
#endif
          SNOW1  = SNOW1 + SNWF - MELTS - CMPRS
c         ACE1I = ACE1I + FREZI + CMPRS
        END IF
      ELSE  ! no existing snow
        HSNOW = MIN(ENRGP,0d0)  ! new snow
        HICE  = HSIL(1)+ENRGP-HSNOW
        SICE  = SSIL(1)
        MELTS = 0.
        CMPRS = 0.
        MELTI = Mi(HICE,SICE,XSI(1)*MSI1)
        SMELTI= MELTI*SICE/(XSI(1)*MSI1)
        FREZI = MIN(RAIN,MAX(-HICE*BYLHM-XSI(1)*MSI1+SICE,0d0))
        SNOW1 = SNWF
c       ACE1I = ACE1I - MELTI + FREZI
#ifdef TRACERS_WATER
        TRSNOW(:) = TRPRCP(:)*SNWF/PRCP
        TRICE(:)  = TRSIL(:,1) + TRPRCP(:)*FREZI/PRCP
        TRMELTI(:)= (MELTI-SMELTI)*TRICE(:)/(XSI(1)*MSI1-SICE)
        TRICE(:)  = TRICE(:) - TRMELTI(:)
        TRRMF(:)  = TRPRCP(:)*(1.-(SNWF+FREZI)/PRCP)
        TRCMPRS(:) = 0.
        TRMELTS(:) = 0.
#endif
      END IF
      SICE = SICE - SMELTI

C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = CMPRS+FREZI-MELTI ! either up or down

C***** Calculate consequential heat/salt flux between layers
      IF (FMSI2.gt.0) THEN      ! downward flux to thermal layer 3
        FHSI2 = FMSI2*HSIL(2)/(XSI(2)*MSI1)
        FSSI2 = FMSI2*SSIL(2)/(XSI(2)*MSI1)
#ifdef TRACERS_WATER
        FTRSI2(:)=FMSI2*TRSIL(:,2)/(XSI(2)*MSI1)
#endif
      ELSE                      ! upward flux
        FHSI2  = FMSI2*HSIL(3)/(XSI(3)*MSI2)
        FSSI2  = FMSI2*SSIL(3)/(XSI(3)*MSI2)
#ifdef TRACERS_WATER
        FTRSI2(:) = FMSI2*TRSIL(:,3)/(XSI(3)*MSI2)
#endif
      END IF

C**** Calculate total snow and ice in upper layers
      IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice originally in layer 1
        SNOW = SNOW1
        HICE = HICE + HSIL(2)
        SICE = SICE + SSIL(2)
#ifdef TRACERS_WATER
        TRICE(:) = TRICE(:) + TRSIL(:,2)
#endif
      ELSE ! some snow orignally in second layer
        SNOW = SNOW1 + (XSI(2)*MSI1-ACE1I)
        SICE = SSIL(2)
        HSNOW2 = (XSI(2)*MSI1-ACE1I)*(Ti(HSIL(2)/(XSI(2)*MSI1)
     *       ,1d3*SSIL(2)/(XSI(2)*MSI1))*shi-lhm)
        HSNOW = HSNOW + HSNOW2
        HICE = HICE + HSIL(2)- HSNOW2
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSNOW(:) + TRSIL(:,2)*(XSI(2)*MSI1-ACE1I)
     *       /(XSI(2)*MSI1-SSIL(2))
        TRICE(:) = TRICE(:) + TRSIL(:,2)*(ACE1I-SSIL(2))/(XSI(2)*MSI1
     *       -SSIL(2))
#endif
      END IF

C**** reconstitute upper layers
      MSI1 = ACE1I + SNOW
      IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
        HSIL(1) = (HICE-FHSI2)*(ACE1I-XSI(2)*MSI1)/ACE1I + HSNOW
        SSIL(1) = (SICE-FSSI2)*(ACE1I-XSI(2)*MSI1)/ACE1I ! + SSNOW
#ifdef TRACERS_WATER
        TRSIL(:,1) = (TRICE(:)-FTRSI2(:))*(ACE1I-XSI(2)*MSI1)/ACE1I
     *       +TRSNOW(:)
#endif
      ELSE
        HSIL(1) = HSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
        SSIL(1) = 0.            ! SSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
#ifdef TRACERS_WATER
        TRSIL(:,1) = TRSNOW(:)*XSI(1)*MSI1/(MSI1-ACE1I)
#endif
      END IF
      HSIL(2) = HICE - HSIL(1) - FHSI2 + HSNOW
      SSIL(2) = SICE - SSIL(1) - FSSI2 ! + SSNOW
#ifdef TRACERS_WATER
      TRSIL(:,2) = TRICE(:) - TRSIL(:,1) - FTRSI2(:) + TRSNOW(:)
#endif

C**** Calculate mass/heat flux between layers 3 and 4
      FMSI3 = XSI(4)*FMSI2
      IF (FMSI3.gt.0) THEN      ! downward flux to layer 4
        FHSI3 = FMSI3*HSIL(3)/(XSI(3)*MSI2)
        FSSI3 = FMSI3*SSIL(3)/(XSI(3)*MSI2)
#ifdef TRACERS_WATER
        FTRSI3(:) = FMSI3*TRSIL(:,3)/(XSI(3)*MSI2)
#endif
      ELSE                      ! upward flux
        FHSI3 = FMSI3*HSIL(4)/(XSI(4)*MSI2)
        FSSI3 = FMSI3*SSIL(4)/(XSI(4)*MSI2)
#ifdef TRACERS_WATER
        FTRSI3(:) = FMSI3*TRSIL(:,4)/(XSI(4)*MSI2)
#endif
      END IF

C**** Adjust amounts resulting from fluxes.
      MSI2 = MSI2 + FMSI2

      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3

      SSIL(3)=SSIL(3)+(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)+ FSSI3

#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,3)=TRSIL(:,3)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)+ FTRSI3(:)
#endif

C**** Diagnostics for output
      RUN0 = MELTS+MELTI+(RAIN-FREZI) ! runoff mass to ocean
      IF (RUN0.lt.1d-13) RUN0=0. ! clean up roundoff errors
      SRUN0 = SMELTI             ! salt in runoff
c     ERUN0 = HMELTI             ! energy of runoff

#ifdef TRACERS_WATER
      TRUN0(:)=TRMELTS(:)+TRMELTI(:)+TRRMF(:)  ! tracer flux to ocean
#endif
      END IF

 ! ice temp L=1,..LMI
      call tice(hsil,ssil,msi1,msi2,tsil)

      RETURN
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,EVAP
     &     ,SROX,
#ifdef TRACERS_WATER
     &     TRSIL,TREVAP,FTROC,TRUN,
#endif
     *     FMOC,FHOC,FSOC,RUN,ERUN,SRUN,WETSNOW,MELT12)
!@sum  SEA_ICE applies surface fluxes to ice covered areas
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER :: dSNdML =0.
!@var DTSRCE source time step (s)
      REAL*8, INTENT(IN) :: DTSRCE
!@var F0DT heat flux on the ice top surface (J/m^2)
      REAL*8, INTENT(IN) :: F0DT
!@var F1DT heat flux between the 1st and 2nd ice layers (J/m^2)
      REAL*8, INTENT(IN) :: F1DT
!@var SROX solar radiation at top and bottom ice surface (J/m^2)
      REAL*8, INTENT(INOUT) :: SROX(2)
!@var EVAP evaporation/dew on the top ice surface (kg/m^2)
      REAL*8, INTENT(IN) :: EVAP
!@var FSRI fraction of solar radiation that goes through the bottom
      REAL*8 :: FSRI(LMI)
!@var FMOC,FHOC,FSOC basal fluxes down of mass,heat,salt (J or kg/m^2)
      REAL*8, INTENT(IN) :: FMOC, FHOC, FSOC
!@var RUN,ERUN,SRUN runoff fluxes down of mass,heat,salt (J or kg/m^2)
      REAL*8, INTENT(OUT) :: RUN,ERUN,SRUN
!@var WETSNOW whether snow is wet or not (Used for albedo calc)
      LOGICAL, INTENT(INOUT) :: WETSNOW
!@var MELT12 save surface melt (Used for albedo calc)
      REAL*8, INTENT(OUT) :: MELT12
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TREVAP tracer amount in evap (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TREVAP
!@var FTROC basal tracer flux to ocean (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: FTROC
!@var TRUN tracer runoff flux to ocean (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN
      REAL*8, DIMENSION(NTM) :: FTRSI2,FTRSI3,TRCMPRS,
     *     TRMELTI,TRMELTS,TRMELT3,TRMELT4,TRDEW,TRICE,TRSNOW
      INTEGER N
#endif
!@var ROICE sea ice fraction of open water
      REAL*8, INTENT(IN) :: ROICE
!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var SNOW snow mass (kg/m^2)
!@var MSI1 first layer ice mass (= SNOW + ACE1I) (kg/m^2)
!@var MSI2 second layer ice mass (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
!@var TSIL temperature of ice layers (C)
      REAL*8, DIMENSION(LMI) :: TSIL
      REAL*8 MSI1, MELTS, MELTI, MELT3, MELT4, DSNOW
      REAL*8 DEW, CMPRS, DEWS, DEWI
      REAL*8 SMELTI, SMELT3, SMELT4, FSSI2, FSSI3, SICE
      REAL*8 FMSI2, FMSI3, FHSI2, FHSI3
      REAL*8 HC1, HC2, HC3, HC4, HICE, HSNOW, HCMPRS
      REAL*8 dF2dTI, dF3dTI, F2, F3, FO ! ,dF1dTI, dF4dTI
      REAL*8 ALAM2,ALAM3

      FMSI2=0. ; FHSI2=0. ; FHSI3=0. ; FSSI2=0. ; FSSI3=0.
      MELTS=0. ; MELTI=0. ; MELT3=0. ; MELT4=0.
      SMELTI=0 ; SMELT3=0.; SMELT4=0.
#ifdef TRACERS_WATER
      FTRSI2=0. ; FTRSI3=0.
      TRMELTS=0. ; TRMELTI=0. ; TRMELT3=0. ; TRMELT4=0.
#endif
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C**** Calculate solar fractions
      IF (SROX(1).gt.0) THEN
        call solar_ice_frac(snow,msi2,wetsnow,fsri,4)
      ELSE
        FSRI = 0.
      END IF
C****
C**** OCEAN ICE, CALCULATE TSIL FROM ENTHALPY
C****
 ! ice temp L=1,... Lmi
      call tice(hsil,ssil,msi1,msi2,tsil)

      HC1 = SHI*XSI(1)*MSI1 ! heat capacity of ice L=1 (J/(degC*m^2))
      HC2 = SHI*XSI(2)*MSI1 ! heat capacity of ice L=2 (J/(degC*m^2))
      HC3 = SHI*XSI(3)*MSI2 ! heat capacity of ice L=3 (J/(degC*m^2))
      HC4 = SHI*XSI(4)*MSI2 ! heat capacity of ice L=4 (J/(degC*m^2))
C**** CALCULATE AND APPLY DIFFUSIVE AND SURFACE ENERGY FLUXES
      ALAM2=ALAMI  ; ALAM3=ALAMI  ! ; ALAM4=ALAMI
      IF (FORM.eq."BP") THEN  ! salinity dependence of diffusion
        IF (SSIL(2).gt.0) ALAM2=ALAM2+ALAMDS*1d3*(SSIL(2)*XSI(3)*MSI2
     *       /(XSI(2)*MSI1)+SSIL(3)*XSI(2)*MSI1/(XSI(3)*MSI2))/(TSIL(2)
     *       *XSI(3)*MSI2+TSIL(3)*XSI(2)*MSI1)
        IF (SSIL(3).gt.0) ALAM3=ALAM3+ALAMDS*1d3*(SSIL(3)*XSI(4)/XSI(3)
     *       +SSIL(4)*XSI(3)/XSI(4))/(TSIL(3)*XSI(4)+TSIL(4)*XSI(3))
c        ALAM4=ALAM3+ALAMDS*1d3*SSIL(4)/(XSI(4)*MSI2*TSIL(4))
      END IF
C**** First layer flux calculated already as part of surface calculation
c      dF1dTI = 2.*DTSRCE/(ACE1I*BYRLI+SNOW*BYRLS)
      dF2dTI = ALAM2*RHOI*DTSRCE/(0.5*XSI(2)*MSI1+0.5*XSI(3)*MSI2)
C****          temperature derivative from F2 diffusive flux
      dF3dTI = ALAM3*RHOI*DTSRCE*2./MSI2
C****          temperature derivative from F3 diffusive flux
c      dF4dTI = ALAM4*RHOI*DTSRCE*2.*BYXSI(4)/MSI2  ! now in iceocean
C****          temperature derivative from F4 diffusive flux 
C**** EXPLCIIT DIFFUSIVE FLUXES
CEXP  F2 = dF2dTI*(TSIL(2)-TSIL(3))+SROX(1)*FSRI(2)
CEXP  F3 = dF3dTI*(TSIL(3)-TSIL(4))+SROX(1)*FSRI(3)
CEXP  FO = dF4dTI*(TSIL(4)-TFO)    +SROX(1)*FSRI(4)
C**** DIFFUSIVE FLUXES FROM PARTIALLY IMPLICIT METHOD
c      F1DT=(dF1dTI*(HC1*(TSIL(1)-TSIL(2))+ALPHA*F0DT)
c     *    +HC1*SROX(1)*FSRI(1))/(HC1+ALPHA*dF1dTI) ! already calculated
      F2=(dF2dTI*(HC2*(TSIL(2)-TSIL(3))+ALPHA*F1DT)+HC2*SROX(1)*FSRI(2))
     *     /(HC2+ALPHA*dF2dTI)
      F3=(dF3dTI*(HC3*(TSIL(3)-TSIL(4))+ALPHA*F2)+HC3*SROX(1)*FSRI(3))/
     *     (HC3+ALPHA*dF3dTI)
c      FO=(dF4dTI*(HC4*(TSIL(4)-TFO)+ALPHA*F3)+HC4*SROX(1)*FSRI(4))/
c     *     (HC4+ALPHA*dF4dTI)
c      FHOC=(dF4dTI*ALPHA*F3+HC4*(SROX(1)*FSRI(4)+FHOC))/
c     *     (HC4+ALPHA*dF4dTI)
C**** Add solar flux through bottom to basal heat flux already computed
      SROX(2) = SROX(1)*FSRI(4)

      HSIL(1) = HSIL(1)+(F0DT-F1DT)
      HSIL(2) = HSIL(2)+(F1DT-F2)
      HSIL(3) = HSIL(3)+(F2-F3)
      HSIL(4) = HSIL(4)+(F3-SROX(2)-FHOC)

C**** add basal salt flux
      SSIL(4) = SSIL(4)-FSOC
#ifdef TRACERS_WATER
      TRSIL(:,4) = TRSIL(:,4)-FTROC(:)
      TRDEW(:) = -TREVAP(:)
#endif

      DEW = -EVAP               ! dew/evap to the surface
      SICE=SSIL(1)+SSIL(2)      ! total upper layer salt (only in ice)

C**** Add DEW directly to ice (which can be in either layer)
C**** DEWI adds to ice, DEWS adds to snow (if applicable)
C**** Calculate melting in first two thermal layers
C**** MELTS,MELTI are the amounts of snow and ice melt
      IF (SNOW*XSI(2).gt.XSI(1)*ACE1I) THEN ! first layer is all snow
        DEWI = MAX(0d0,DEW)    ! >0 i.e. dew to second layer ice
        DEWS = DEW - DEWI      ! <0 remaining evap applied to snow
        HSNOW = HSIL(1) + (XSI(2)*MSI1+DEWI-ACE1I)*(Ti(HSIL(2)/(XSI(2)
     *       *MSI1+DEWI),1d3*SSIL(2)/(XSI(1)*MSI1+DEWI))*shi-lhm)
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSIL(:,1) + TRSIL(:,2)*MAX(1.-(ACE1I-SICE)
     *       /(XSI(2)*MSI1-SICE),0d0)
#endif
      ELSE  ! first layer is snow and some ice
        DEWS = -MIN(SNOW,-(DEW-MAX(0d0,DEW))) ! <0 evap upto snow amount
        DEWI = DEW - DEWS       ! either sign, evap/dew applied to ice
        HSNOW = (SNOW+DEWS)*(Ti(HSIL(1)/(XSI(1)*MSI1+DEW),1d3*SSIL(1)
     *       /(XSI(1)*MSI1+DEW))*shi-lhm) 
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSIL(:,1)*MIN(SNOW/(XSI(1)*MSI1-SSIL(1)),1d0)
#endif
      END IF

#ifdef TRACERS_WATER
      DO N=1,NTM
C**** Take into account occasional sign differences when summing
C**** over multiple surface time steps
        TRICE(N) = TRSIL(N,1) + TRSIL(N,2) - TRSNOW(N)
        IF (DEWI.gt.0. .and. TRDEW(N).gt.0) THEN
          TRICE(N) = TRICE(N)+TRDEW(N)
        ELSEIF (TRSNOW(N)+TRDEW(N).gt.0) THEN
          TRSNOW(N) = TRSNOW(N)+TRDEW(N)
        ELSE
          TRICE(N) = TRICE(N)+TRSNOW(N)+TRDEW(N)
          TRSNOW(N) = 0.
        END IF
      END DO
#endif

C**** calculate melt in snow and ice layers
      HICE = HSIL(1)+HSIL(2)-HSNOW
      MELTS = MAX(0d0,HSNOW*BYLHM+SNOW+DEWS)  ! snow melt amount
      MELTI = Mi(HICE,SICE,ACE1I+DEWI) 
      SMELTI= MELTI*SICE/(ACE1I+DEWI)
C**** CMPRS is amount of snow turned to ice during melting
C**** Calculate remaining snow and necessary mass flux using
C**** SNOW =SNOW - MELTS + DEWS - CMPRS
C**** ACE1I=ACE1I- MELTI + DEWI + CMPRS - FMSI2
C**** SICE =SICE- SMELTI                - FSSI2
C****
      IF (SNOW.GT.MELTS-DEWS) THEN ! some snow remains
        CMPRS = MIN(dSNdML*MELTS, SNOW+DEWS-MELTS) ! > 0.
        DSNOW = - (MELTS-DEWS+CMPRS)
        HCMPRS = HSNOW*CMPRS/(SNOW+DEWS-MELTS)
#ifdef TRACERS_WATER
        TRMELTS(:)= MELTS*TRSNOW(:)/(SNOW+DEWS)
        TRCMPRS(:)= CMPRS*(TRSNOW(:)-TRMELTS(:))/(SNOW+DEWS-MELTS)
        TRMELTI(:)= (MELTI-SMELTI)*TRICE(:)/(ACE1I-SICE+DEWI)
#endif
      ELSE ! all snow and some ice melts or evaporates
        CMPRS = 0.
        DSNOW = -SNOW
        HCMPRS=0.
#ifdef TRACERS_WATER
        TRCMPRS(:)= 0.
        TRMELTS(:)= TRSNOW(:)
C**** tracer ice melt is calculated using fresh water concentration
        TRMELTI(:)= (MELTI-SMELTI)*TRICE(:)/(ACE1I-SICE+DEWI)
#endif
      END IF

C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = -MELTI+DEWI+CMPRS  ! either up or down

C**** Check for melting in levels 3 and 4
      MELT3 = Mi(HSIL(3),SSIL(3),XSI(3)*MSI2)
      MELT4 = Mi(HSIL(4),SSIL(4),XSI(4)*MSI2-FMOC)
      SMELT3 = MELT3*SSIL(3)/(XSI(3)*MSI2)
      SMELT4 = MELT4*SSIL(4)/(XSI(4)*MSI2-FMOC)
#ifdef TRACERS_WATER
      TRMELT3(:) = MELT3*TRSIL(:,3)/(XSI(3)*MSI2)
      TRMELT4(:) = MELT4*TRSIL(:,4)/(XSI(4)*MSI2-FMOC)
#endif
C***** Calculate consequential heat flux between layers
      IF (FMSI2.gt.0) THEN ! downward flux to thermal layer 3
        FHSI2 = FMSI2*(HICE+HCMPRS)/(ACE1I+DEWI-MELTI+CMPRS)
        FSSI2 = FMSI2*(SICE-SMELTI)/(ACE1I+DEWI-MELTI+CMPRS)
#ifdef TRACERS_WATER
        FTRSI2(:)=FMSI2*(TRICE(:)-TRMELTI(:)+TRCMPRS(:))/(ACE1I+DEWI
     *       -MELTI+CMPRS)
#endif
      ELSE                 ! upward flux
        FHSI2  = FMSI2*HSIL(3)/(XSI(3)*MSI2-MELT3)
        FSSI2  = FMSI2*(SSIL(3)-SMELT3)/(XSI(3)*MSI2-MELT3)
#ifdef TRACERS_WATER
        FTRSI2(:)=FMSI2*(TRSIL(:,3)-TRMELT3(:))/(XSI(3)*MSI2-MELT3)
#endif
      END IF

C**** reconstitute upper layers
      SNOW = MAX(0d0,SNOW+DSNOW)
      MSI1 = ACE1I + SNOW
      IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
        HSIL(1) = (HICE+HCMPRS-FHSI2)*(ACE1I-XSI(2)*MSI1)/ACE1I + HSNOW
     *       -HCMPRS
        SSIL(1) = (SICE-SMELTI-FSSI2)*(ACE1I-XSI(2)*MSI1)/ACE1I ! +SSNOW
#ifdef TRACERS_WATER
        TRSIL(:,1) = (TRICE(:)-TRMELTI(:)+TRCMPRS(:)-FTRSI2(:))*(ACE1I
     *       -XSI(2)*MSI1)/ACE1I +TRSNOW(:)-TRMELTS(:)-TRCMPRS(:)
#endif
      ELSE
        HSIL(1) = (HSNOW-HCMPRS)*XSI(1)*MSI1/(MSI1-ACE1I)
        SSIL(1) = 0.            ! SSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
#ifdef TRACERS_WATER
        TRSIL(:,1) = (TRSNOW(:)-TRMELTS(:)-TRCMPRS(:))*XSI(1)*MSI1/(MSI1
     *       -ACE1I)
#endif
      END IF
      HSIL(2) = HICE          - HSIL(1) - FHSI2 + HSNOW
      SSIL(2) = SICE - SMELTI - SSIL(1) - FSSI2 ! + SSNOW
#ifdef TRACERS_WATER
      TRSIL(:,2) = MAX(TRICE(:) - TRMELTI(:) - TRSIL(:,1) - FTRSI2(:) +
     *     TRSNOW(:) - TRMELTS(:),0d0)
#endif

C**** Calculate mass/heat flux between layers 3 and 4
      FMSI3 = XSI(3)*(MELT4+FMOC) + XSI(4)*(FMSI2-MELT3)
      IF (FMSI3.gt.0) THEN ! downward flux to layer 4
        FHSI3 = FMSI3*HSIL(3)/(XSI(3)*MSI2-MELT3)
        FSSI3 = FMSI3*(SSIL(3)-SMELT3)/(XSI(3)*MSI2-MELT3)
#ifdef TRACERS_WATER
        FTRSI3(:)=FMSI3*(TRSIL(:,3)-TRMELT3(:))/(XSI(3)*MSI2-MELT3)
#endif
      ELSE                 ! upward flux
        FHSI3 = FMSI3*HSIL(4)/(XSI(4)*MSI2-MELT4-FMOC)
        FSSI3 = FMSI3*(SSIL(4)-SMELT4)/(XSI(4)*MSI2-MELT4-FMOC)
#ifdef TRACERS_WATER
        FTRSI3(:)=FMSI3*(TRSIL(:,4)-TRMELT4(:))/(XSI(4)*MSI2-MELT4-FMOC)
#endif
      END IF

C**** Apply fluxes to lower layers
      MSI2=MSI2-(MELT3+MELT4)+FMSI2-FMOC

      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3

      SSIL(3)=SSIL(3)-SMELT3 +(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)-SMELT4 + FSSI3

#ifdef TRACERS_WATER
      TRSIL(:,3)=TRSIL(:,3)- TRMELT3(:)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)- TRMELT4(:)+ FTRSI3(:)
#endif
C**** Calculate additional runoff output fluxes
      RUN = MELTS + MELTI + MELT3 + MELT4 ! mass flux to ocean
      SRUN=        SMELTI +SMELT3 +SMELT4 ! salt flux to ocean
c HMELT currently assumed to be zero since melting is at 0 deg
      ERUN= SROX(2)        ! + HMELT
#ifdef TRACERS_WATER
      TRUN(:)=TRMELTS(:)+TRMELTI(:)+TRMELT3(:)+TRMELT4(:)
                                ! tracer flux to ocean
#endif

c**** Decide WETSNOW for albedo calculations and save MELT12 for
c**** possible melt ponds
      MELT12=MELTS+MELTI
      WETSNOW = WETSNOW .or. (MELT12.gt.0)
C****
      RETURN
      END SUBROUTINE SEA_ICE

      SUBROUTINE ADDICE (SNOW,ROICE,HSIL,SSIL,MSI2,TSIL,ENRGFO
     &     ,ACEFO,ACEFI,ENRGFI,SALTO,SALTI,
#ifdef TRACERS_WATER
     &     TRSIL,TRO,TRI,DTRIMP,
#endif
     *     DMIMP,DHIMP,DSIMP,FLEAD,QFIXR)
!@sum  ADDICE adds ice formed in the ocean to ice variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER, DIMENSION(LMI) :: YSI =
     *     (/XSI(1)* ACE1I/(ACE1I+AC2OIM),XSI(2)*ACE1I/(ACE1I+AC2OIM),
     *       XSI(3)*AC2OIM/(ACE1I+AC2OIM),XSI(4)*AC2OIM/(ACE1I+AC2OIM)/)
c     REAL*8, PARAMETER :: Z2OIX = 5.d0-Z1I,     ! max ice thickness l=2
c    *                     BYZICX=1./(Z1I+Z2OIX)
!@var QFIXR true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL, INTENT(IN) :: QFIXR
!@var FLEAD minimum lead fraction for ice (%)
      REAL*8, INTENT(IN) :: FLEAD
!@var ACEFO ice mass formed in ocean for open water fraction (kg/m^2)
!@var ACEFI ice mass formed in ocean for ice covered fraction (kg/m^2)
!@var ENRGFO energy of ice formed in ocean for open water frac (J/m^2)
!@var ENRGFI energy of ice formed in ocean for ice covered frac (J/m^2)
!@var SALTO salt in ice formed in ocean for open water frac (kg/m^2)
!@var SALTI salt in ice formed in ocean for ice covered frac (kg/m^2)
      REAL*8, INTENT(IN) ::  ENRGFI, ENRGFO, ACEFO, ACEFI, SALTO, SALTI
!@var ROICE  ice fraction over open water
!@var SNOW snow amount on ice (kg/m^2)
!@var MSI2 second mass layer ice thickness (kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL,SSIL
!@var TSIL temperature of ice layers (C)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var DMIMP,DSIMP,DHIMP implicit mass,salt and heat fluxes required
!@+   to maintain minimum ice thickness if ice fraction is fixed
      REAL*8, INTENT(OUT) :: DMIMP,DSIMP,DHIMP
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: trsil
!@var TRO tracer in ice formed in ocean for open water frac (kg/m^2)
!@var TRI tracer in ice formed in ocean for ice covered frac (kg/m^2)
      REAL*8, INTENT(IN), DIMENSION(NTM) :: tro,tri
!@var DTRIMP implicit tracer flux required to maintain minimum ice
!@+   thickness if ice fraction is fixed
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: DTRIMP
      REAL*8, DIMENSION(NTM) :: FTRSI3,FTRSI4,TRSNOW,TRICE
      INTEGER N
#endif
      REAL*8, DIMENSION(LMI) :: FRI
      REAL*8 FMSI4, FHSI3, FHSI4, FSSI3, FSSI4, HSNOW, HICE, SICE
      REAL*8 ROICEN, OPNOCN, DRSI, MSI1

      DMIMP=0. ; DHIMP=0. ; DSIMP=0.
#ifdef TRACERS_WATER
      DTRIMP=0.
#endif

      MSI1=SNOW+ACE1I
      IF (.not. QFIXR) THEN
      IF (ROICE.LE.0. .and. ACEFO.gt.0) THEN
C**** Create new ice in ice-free ocean
        ROICE=ACEFO/(ACE1I+AC2OIM)
        SNOW=0.
C****   TSIL=Ti(ENRGFO/ACEFO,1d3*SALTO/ACEFO) ! but hsi is primary var.
        HSIL(1:2) =(ENRGFO/ACEFO)*XSI(1:2)*ACE1I
        HSIL(3:4) =(ENRGFO/ACEFO)*XSI(3:4)*AC2OIM
        SSIL(1:2) = (SALTO/ACEFO)*XSI(1:2)*ACE1I
        SSIL(3:4) = (SALTO/ACEFO)*XSI(3:4)*AC2OIM
#ifdef TRACERS_WATER
        DO N=1,NTM
          TRSIL(N,1:2) = (TRO(N)/ACEFO)*XSI(1:2)*ACE1I
          TRSIL(N,3:4) = (TRO(N)/ACEFO)*XSI(3:4)*AC2OIM
        END DO
#endif
        MSI1=ACE1I
        MSI2=AC2OIM
      ELSEIF (ROICE.gt.0) THEN
C**** Create new ice in partially ice-covered ocean
      IF (ACEFI.gt.0) THEN
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
CC      FMSI3 = -XSI(3)*ACEFI ! < 0.
CC      FMSI4 = -ACEFI
        FHSI3 = -HSIL(4)*ACEFI*(XSI(3)/XSI(4))/MSI2
        FSSI3 = -SSIL(4)*ACEFI*(XSI(3)/XSI(4))/MSI2
#ifdef TRACERS_WATER
        FTRSI3(:) = -TRSIL(:,4)*ACEFI*(XSI(3)/XSI(4))/MSI2
#endif
      ELSE
        FHSI3 = 0. ; FSSI3 = 0.
#ifdef TRACERS_WATER
        FTRSI3(:) = 0.
#endif
      END IF

      IF (ACEFO .EQ. 0.) THEN
C**** NEW ICE IS ONLY FORMED BELOW OLD SEA ICE
        HSIL(3) = HSIL(3)-FHSI3
        HSIL(4) = HSIL(4)+(FHSI3+ENRGFI)
        SSIL(3) = SSIL(3)-FSSI3
        SSIL(4) = SSIL(4)+(FSSI3+SALTI)
#ifdef TRACERS_WATER
        TRSIL(:,3) = TRSIL(:,3) - FTRSI3(:)
        TRSIL(:,4) = TRSIL(:,4) +(FTRSI3(:)+TRI(:))
#endif
        MSI2 = MSI2+ACEFI       ! new ice mass of physical layer 2
      ELSE
C**** NEW ICE IS FORMED ON OPEN OCEAN AND POSSIBLY BELOW OLD SEA ICE
        DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
        ROICEN=ROICE+DRSI
C**** New formulation:
C**** Split all upper mass layer fields into ice and snow components
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          HSNOW = SNOW*(Ti(HSIL(1)/(XSI(1)*MSI1),1d3*SSIL(1)/(XSI(1)
     *         *MSI1))*shi-lhm)
        ELSE
          HSNOW = HSIL(1)+(SNOW-XSI(1)*MSI1)*(Ti(HSIL(2)/(XSI(2)*MSI1)
     *         ,1d3*SSIL(2)/(XSI(2)*MSI1))*shi-lhm)
        END IF
        HICE  = HSIL(1)+HSIL(2)-HSNOW
CC      SSNOW = 0.   ! always zero
        SICE  = SSIL(1)+SSIL(2)
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSIL(:,1)*MIN(SNOW/(XSI(1)*MSI1-SSIL(1)),1d0) +
     *       TRSIL(:,2)*MAX(1.-(ACE1I-SICE)/(XSI(2)*MSI1-SSIL(2)),0d0)
        TRICE(:)  = TRSIL(:,1)+TRSIL(:,2)-TRSNOW(:)
#endif
C**** distribute snow variables over new ice extent
        SNOW = SNOW*(ROICE/ROICEN)
        HSNOW = HSNOW*(ROICE/ROICEN)
CC      SSNOW = SSNOW*(ROICE/ROICEN)  ! always zero
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSNOW(:)*(ROICE/ROICEN)
#endif
C**** Add new ice to ice variables
        HICE=((1.-ROICE)*ENRGFO*ACE1I/(ACE1I+AC2OIM)+ROICE*HICE)/ROICEN
        SICE=((1.-ROICE)*SALTO *ACE1I/(ACE1I+AC2OIM)+ROICE*SICE)/ROICEN
#ifdef TRACERS_WATER
        TRICE(:)=((1.-ROICE)*TRO(:)*ACE1I/(ACE1I+AC2OIM)+ROICE*TRICE(:))
     *       /ROICEN
#endif
C**** reconstitute upper thermal layers
        MSI1 = SNOW + ACE1I     ! mass of layer 1
CC      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          HSIL(1) = HICE*(ACE1I-XSI(2)*MSI1)/ACE1I + HSNOW
          SSIL(1) = SICE*(ACE1I-XSI(2)*MSI1)/ACE1I ! + SSNOW
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRICE(:)*(ACE1I-XSI(2)*MSI1)/ACE1I + TRSNOW(:)
#endif
        ELSE
          HSIL(1) = HSNOW*XSI(1)*MSI1/SNOW
          SSIL(1) = 0.          ! SSNOW*XSI(1)*MSI1/SNOW
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRSNOW(:)*XSI(1)*MSI1/SNOW
#endif
        END IF
        HSIL(2) = HICE - HSIL(1) + HSNOW
        SSIL(2) = SICE - SSIL(1) ! + SSNOW
#ifdef TRACERS_WATER
        TRSIL(:,2) = TRICE(:) - TRSIL(:,1) + TRSNOW(:)
#endif
C**** add new ice to old ice (not including snow)
        MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACEFI))/ROICEN ! layer 2
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
        HSIL(3)=((1.-ROICE)*ENRGFO*YSI(3)+ROICE*(HSIL(3)-FHSI3))/ROICEN
        HSIL(4)=((1.-ROICE)*ENRGFO*YSI(4)+ROICE*(HSIL(4)+FHSI3+ENRGFI)
     *       )/ROICEN
        SSIL(3)=((1.-ROICE)*SALTO*YSI(3)+ROICE*(SSIL(3)-FSSI3))/ROICEN
        SSIL(4)=((1.-ROICE)*SALTO*YSI(4)+ROICE*(SSIL(4)+FSSI3+SALTI))
     *       /ROICEN
#ifdef TRACERS_WATER
        DO N=1,NTM
          TRSIL(N,3)=((1.-ROICE)*TRO(N)*YSI(3)+ROICE*(TRSIL(N,3)
     *         -FTRSI3(N)))/ROICEN
          TRSIL(N,4)=((1.-ROICE)*TRO(N)*YSI(4)+ROICE*(TRSIL(N,4)
     *         +FTRSI3(N)+TRI(N)))/ROICEN
        END DO
#endif
        ROICE = ROICEN          ! new ice concentration
      END IF
C**** COMPRESS THE ICE HORIZONTALLY IF TOO THIN OR LEAD FRAC. TOO SMALL
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ROICE*(ACE1I+MSI2)))
      IF ((ROICE*(ACE1I+MSI2)).gt.5.*RHOI) OPNOCN=0.  ! no leads for h>5
      IF (MSI2.LT.AC2OIM .or. ROICE.GT.1.-OPNOCN) THEN
      ROICEN = MIN(ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM),1.-OPNOCN)
      DRSI = ROICEN-ROICE ! < 0. compressed ice concentration
CC    FMSI3 = XSI(3)*FMSI4 ! < 0. upward ice mass flux into layer 3
      FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! upward ice mass into layer 4
      FHSI3 = HSIL(4)*FMSI4*(XSI(3)/XSI(4))/MSI2 ! upward heat flux
      FHSI4 = (HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4))*(DRSI/ROICEN)
      HSIL(3) = HSIL(3)-FHSI3
      HSIL(4) = HSIL(4)+(FHSI3-FHSI4)
      FSSI3 = SSIL(4)*FMSI4*(XSI(3)/XSI(4))/MSI2 ! upward salt flux
      FSSI4 = (SSIL(1)+SSIL(2)+SSIL(3)+SSIL(4))*(DRSI/ROICEN)
      SSIL(3) = SSIL(3)-FSSI3
      SSIL(4) = SSIL(4)+(FSSI3-FSSI4)
#ifdef TRACERS_WATER
      FTRSI3(:) = TRSIL(:,4)*FMSI4*(XSI(3)/XSI(4))/MSI2
      FTRSI4(:) = (TRSIL(:,1)+TRSIL(:,2)+TRSIL(:,3)+TRSIL(:,4))*(DRSI
     *     /ROICEN)
      TRSIL(:,3) = TRSIL(:,3) - FTRSI3(:)
      TRSIL(:,4) = TRSIL(:,4) +(FTRSI3(:)-FTRSI4(:))
#endif
      MSI2 = MSI2-FMSI4 ! new ice mass of second physical layer
CC    SNOW = SNOW   ! snow thickness is conserved
      ROICE = ROICEN
C****
      END IF
C****
      END IF
C**** Clean up ice fraction (if rsi>(1-OPNOCN)-1d-3) => rsi=(1-OPNOCN))
      IF (ROICE.gt.0) THEN
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ROICE*(ACE1I+MSI2)))    ! -BYZICX)
      IF ((ROICE*(ACE1I+MSI2)).gt.5.*RHOI) OPNOCN=0.  ! no leads for h>5
      IF (ROICE.gt.(1.-OPNOCN-1d-3)) THEN
        ROICEN = 1.-OPNOCN
        DRSI = ROICEN-ROICE ! +ve
        FMSI4 = (ACE1I+MSI2)*(DRSI/ROICEN) ! >0 ice mass out of layer 4
        FHSI4 = HSIL(4)*FMSI4/(XSI(4)*MSI2)
        FSSI4 = SSIL(4)*FMSI4/(XSI(4)*MSI2)
        FHSI3 = HSIL(3)*FMSI4/MSI2 ! downward heat flux into layer 4
        FSSI3 = SSIL(3)*FMSI4/MSI2 ! downward salt flux into layer 4
#ifdef TRACERS_WATER
        FTRSI4(:) = TRSIL(:,4)*FMSI4/(XSI(4)*MSI2)
        FTRSI3(:) = TRSIL(:,3)*FMSI4/MSI2
#endif
        MSI2=MSI2-FMSI4         ! new ice mass of second physical layer
        FRI(1)=ACE1I/(ACE1I+MSI2)
        FRI(3:4)=XSI(3:4)*MSI2/(ACE1I+MSI2)
C**** New formulation:
C**** Split all upper mass layer fields into ice and snow components
C**** Assume equal temperatures over snow and ice
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          HSNOW = SNOW*(Ti(HSIL(1)/(XSI(1)*MSI1),1d3*SSIL(1)/(XSI(1)
     *         *MSI1))*shi-lhm) 
        ELSE
          HSNOW = HSIL(1)+(SNOW-XSI(1)*MSI1)*(Ti(HSIL(2)/(XSI(2)*MSI1)
     *         ,1d3*SSIL(2)/(XSI(2)*MSI1))*shi-lhm)
        END IF
        HICE  = HSIL(1)+HSIL(2)-HSNOW
CC      SSNOW = 0.   ! always zero
        SICE  = SSIL(1)+SSIL(2)
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSIL(:,1)*MIN(SNOW/(XSI(1)*MSI1-SSIL(1)),1d0) +
     *       TRSIL(:,2)*MAX(1.-(ACE1I-SICE)/(XSI(2)*MSI1-SSIL(2)),0d0)
        TRICE(:)  = TRSIL(:,1)+TRSIL(:,2)-TRSNOW(:)
#endif
C**** distribute snow variables over new ice extent
        SNOW = SNOW*(ROICE/ROICEN)
        HSNOW = HSNOW*(ROICE/ROICEN)
CC      SSNOW = SSNOW*(ROICE/ROICEN)  ! always zero
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSNOW(:)*(ROICE/ROICEN)
#endif
C**** Add new ice to ice variables
        HICE = (ROICE/ROICEN)*(FHSI4*FRI(1)+HICE)
        SICE = (ROICE/ROICEN)*(FSSI4*FRI(1)+SICE)
#ifdef TRACERS_WATER
        TRICE(:)=(ROICE/ROICEN)*(FTRSI4(:)*FRI(1)+TRICE(:))
#endif
C**** reconstitute upper thermal layers
        MSI1 = SNOW + ACE1I     ! mass of layer 1
CC      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          FRI(1)=(ACE1I-XSI(2)*MSI1)/ACE1I
          HSIL(1) = HICE*FRI(1) + HSNOW
          SSIL(1) = SICE*FRI(1) ! + SSNOW
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRICE(:)*FRI(1) + TRSNOW(:)
#endif
        ELSE
          HSIL(1) = HSNOW*XSI(1)*MSI1/SNOW
          SSIL(1) = 0.          ! SSNOW*XSI(1)*MSI1/SNOW
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRSNOW(:)*XSI(1)*MSI1/SNOW
#endif
        END IF
        HSIL(2) = HICE - HSIL(1) + HSNOW
        SSIL(2) = SICE - SSIL(1) ! + SSNOW
#ifdef TRACERS_WATER
        TRSIL(:,2) = TRICE(:) - TRSIL(:,1) + TRSNOW(:)
#endif
        HSIL(3) = (ROICE/ROICEN)*(HSIL(3)+FHSI4*FRI(3)-FHSI3)
        HSIL(4) = (ROICE/ROICEN)*(HSIL(4)+FHSI4*FRI(4)+FHSI3-FHSI4)
        SSIL(3) = (ROICE/ROICEN)*(SSIL(3)+FSSI4*FRI(3)-FSSI3)
        SSIL(4) = (ROICE/ROICEN)*(SSIL(4)+FSSI4*FRI(4)+FSSI3-FSSI4)
#ifdef TRACERS_WATER
        DO N=1,NTM
          TRSIL(N,3) = (ROICE/ROICEN)*(TRSIL(N,3)+FTRSI4(N)*FRI(3)
     *         -FTRSI3(N))
          TRSIL(N,4) = (ROICE/ROICEN)*(TRSIL(N,4)+FTRSI4(N)*FRI(4)
     *         +FTRSI3(N)-FTRSI4(N))
        END DO
#endif
        ROICE = ROICEN          ! new ice concentration
      END IF
      END IF
      ELSE
C**** Ensure that MSI2 does not get too small for fixed-SST case.
        IF (ROICE.gt.0. and. MSI2.lt.AC2OIM) then
          DMIMP=AC2OIM-MSI2
          DHIMP=SUM(HSIL(3:4)*XSI(3:4)*DMIMP)
          DSIMP=SUM(SSIL(3:4)*XSI(3:4)*DMIMP)
          HSIL(3:4)=HSIL(3:4)*AC2OIM/MSI2
          SSIL(3:4)=SSIL(3:4)*AC2OIM/MSI2
#ifdef TRACERS_WATER
          DTRIMP(:)=(TRSIL(:,3)*XSI(3)+TRSIL(:,4)*XSI(4))*DMIMP
          TRSIL(:,3:4)=TRSIL(:,3:4)*AC2OIM/MSI2
#endif
          MSI2=AC2OIM
        END IF
      END IF
C**** Calculate temperatures for diagnostics and radiation
      call tice(hsil,ssil,msi1,msi2,tsil)

      RETURN
      END SUBROUTINE ADDICE

      SUBROUTINE SIMELT(DT,ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,Tm,TFO,TSIL,
#ifdef TRACERS_WATER
     &     TRSIL,TRUN0,
#endif
     &     ENRGMAX,ENRGUSED,RUN0,SALT)
!@sum  SIMELT melts sea ice laterally and if it is too small
!@+    Note: all amounts are with respect to the ocean/lake fraction
!@auth Original Development Team
!@ver  1.1
      IMPLICIT NONE
!@var POCEAN ocean fraction (zero if lake)
      REAL*8, INTENT(IN) :: POCEAN
!@var TFO freezing temperature of water (C)
      REAL*8, INTENT(IN) :: TFO
!@var ROICE,SNOW,MSI2 ice variables (%,kg/m^2,kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var Tm mixed layer ocean temperature (C)
      REAL*8, INTENT(IN) :: Tm
!@var DT time step (s)
      REAL*8, INTENT(IN) :: DT
!@var HSIL ice enthalpy  (J/m^2)
!@var SSIL ice salt (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var TSIL ice temperature (C)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var ENRGMAX max energy available for melting ice (J/m^2)
      REAL*8, INTENT(IN) :: ENRGMAX
!@var ENRGUSED energy used to melt ice (J/m^2)
      REAL*8, INTENT(OUT) :: ENRGUSED
!@var RUN0 amount of sea ice melt (kg/m^2)
!@var SALT amount of salt in sea ice melt (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0, SALT
!@var DRSI change in RSI fraction
!@var DTEMP temperature diff. used in lateral melt calculation.
      REAL*8 :: DRSI,DTEMP
#ifdef TRACERS_WATER
!@var TRSIL tracer ice amounts (kg/m^2)
!@var TRUN0 tracer amount in sea ice melt (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRUN0
#endif
      INTEGER L

C**** Estimate DRSI
      DRSI=0.
      IF (ROICE.lt.1d-3) THEN
        DRSI=ROICE
      ELSE      ! IF (POCEAN.gt.0) THEN ! now for lakes too
C**** Estimate lateral melt using parameterisation from Makyut/Steele
C**** (via C. Bitz): Rside=dt*pi/(floesize*eta)*(3d-6)*(delT)^(1.36)
        dtemp=MAX(Tm-TFO,0d0)
        DRSI=DT*SILMFAC*dtemp**SILMPOW
        IF (ROICE-DRSI.lt.1d-3) DRSI=ROICE
        IF (ENRGMAX+DRSI*SUM(HSIL).lt.0) DRSI=-ENRGMAX/SUM(HSIL)
      END IF
C**** Remove DRSI amount of ice
      ENRGUSED=-DRSI*SUM(HSIL) !  [J/m^2]
      RUN0=DRSI*(SNOW + ACE1I + MSI2)
      SALT=DRSI*SUM(SSIL)
#ifdef TRACERS_WATER
      TRUN0(:)=DRSI*(TRSIL(:,1)+TRSIL(:,2)+TRSIL(:,3)+TRSIL(:,4))
#endif
      ROICE=ROICE-DRSI
      IF (ROICE.lt.1d-10) THEN
        ROICE=0.                ! deal with possible round off err
C**** set defaults if no ice is left
        SNOW=0.
        MSI2=AC2OIM
        IF (POCEAN.gt.0) THEN
          SSIL(1:2)=SSI0*XSI(1:2)*ACE1I
          SSIL(3:LMI)=SSI0*XSI(3:LMI)*AC2OIM
        ELSE
          SSIL(:) = 0.
        END IF
        DO L=1,2
          HSIL(L)=(XSI(L)*ACE1I)*Ei(TFO,1d3*SSIL(L)/(XSI(L)*ACE1I))
        END DO
        DO L=3,LMI
          HSIL(L)=(XSI(L)*AC2OIM)*Ei(TFO,1d3*SSIL(L)/(XSI(L)*AC2OIM))
        END DO
        TSIL=TFO
#ifdef TRACERS_WATER
        TRSIL(:,:)=0.
#endif
      END IF
C****
      RETURN
      END SUBROUTINE SIMELT

      SUBROUTINE SSIDEC(SNOW,MSI2,HSIL,SSIL,DT,MELT12,
#ifdef TRACERS_WATER
     &     TRSIL,TRFLUX,
#endif
     &     MFLUX,HFLUX,SFLUX)
!@sum  SSIDEC decays salinity in sea ice
!@auth Jiping Liu
!@ver  1.0
      IMPLICIT NONE

!@var HSIL enthalpy of ice layers (J/m^2)
!@var SSIL salt in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: SSIL,HSIL
!@var MSI2 second mass layer ice thickness (kg/m^2)
!@var MSI1 first mass layer ice thickness (kg/m^2)
!@var DT source time step (s)
!@var MELT12 surface melt this time step (kg/m2)
      REAL*8, INTENT(INOUT) :: MSI2
      REAL*8, INTENT(IN) :: DT, SNOW, MELT12
!@var MFLUX,SFLUX,HFLUX mass, salt and heat flux arising from
!@+   sea salinity decay
      REAL*8, INTENT(OUT) :: MFLUX,HFLUX,SFLUX
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
!@var TRFLUX tracer flux arising from sea salinity decay
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRFLUX
      REAL*8, DIMENSION(NTM) :: FTRSI2,FTRSI3,DTR12,TRSNOW,TRICE
      REAL*8, DIMENSION(NTM,3:LMI) :: DTRSI
      INTEGER N
#endif

      REAL*8 DSSI(3:LMI),DMSI(3:LMI),DHSI(3:LMI),TSIL(3:LMI),DS12
     *     ,DM12,DH12
      REAL*8 FMSI2,FMSI3,FSSI2,FSSI3,FHSI2,FHSI3
      REAL*8 HSNOW,HICE,SICE,TICE,MSI1
      INTEGER L
!@var dtssi decay time scale for sea ice salinity (days)
!@var bydtssi decay constant for sea ice salinity (1/s)
      REAL*8, parameter :: dtssi=30.d0, bydtssi=1./(dtssi*sday)
      REAL*8 :: rate,brine_frac

      DSSI = 0. ; DMSI = 0. ; DHSI = 0.
      DS12 = 0. ; DM12 = 0. ; DH12 = 0.
#ifdef TRACERS_WATER
      DTR12(:)  = 0.
      DTRSI(:,:)= 0.
#endif
C**** salinity in sea ice decays if it is above threshold level ssi0

C**** check first layer (default ice and snow)
C**** New formulation:
C**** Split all upper mass layer fields into ice and snow components
C**** Assume equal temperatures over snow and ice
      MSI1=SNOW+ACE1I
      IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
        HSNOW = SNOW*(Ti(HSIL(1)/(XSI(1)*MSI1),1d3*SSIL(1)
     *         /(XSI(1)*MSI1))*shi-lhm)
      ELSE
        HSNOW = HSIL(1)+(XSI(2)*MSI1-ACE1I)*(Ti(HSIL(2)
     *       /(XSI(2)*MSI1),1d3*SSIL(2)/(XSI(2)*MSI1))*shi-lhm)
      END IF
      HICE  = HSIL(1)+HSIL(2)-HSNOW
CC    SSNOW = 0.   ! always zero
      SICE  = SSIL(1)+SSIL(2)
      TICE = Ti(HICE/ACE1I,1d3*SICE/ACE1I)
#ifdef TRACERS_WATER
      TRSNOW(:) = TRSIL(:,1)*MIN(SNOW/(XSI(1)*MSI1-SSIL(1)),1d0)
     *     +TRSIL(:,2)*MAX(1.-(ACE1I-SICE)/(XSI(2)*MSI1-SSIL(2)),0d0)
      TRICE(:)  = TRSIL(:,1)+TRSIL(:,2)-TRSNOW(:)
#endif

C**** brine removal is very sensitive to seaice thermodynamic
C**** formulation
      SELECT CASE (FORM)
      CASE ("SI")               ! salinity only affects mass
C**** calculate removal of excess salinity using simple decay
        IF (SICE.GT.ssi0*ACE1I) THEN
          DS12 = (SICE-ACE1I*ssi0)*DT*BYDTSSI
          DM12 = DS12
          DH12 = -DM12*TICE*SHI
#ifdef TRACERS_WATER
C**** no tracer removed if pure salt is lost
c         DTR12(:) = (DM12-DS12)*TRICE(:)/(ACE1I-SICE)
c         TRICE(:) = TRICE(:)-DTR12(:)
#endif
          SICE=SICE-DS12
          HICE=HICE-DH12
        END IF

C**** check remaining layers
        DO L=3,LMI
          IF (SSIL(L).GT.ssi0*XSI(L)*MSI2) THEN
            DSSI(L)=(SSIL(L)-ssi0*XSI(L)*MSI2)*DT*BYDTSSI
            DMSI(L)=DSSI(L)
            TSIL(L)=Ti(HSIL(L)/(XSI(L)*MSI2),1d3*SSIL(L)/(XSI(L)*MSI2))
            DHSI(L)=-DMSI(L)*TSIL(L)*SHI
#ifdef TRACERS_WATER
C**** no tracer removed if pure salt is lost
c           DTRSI(:,L)=(DMSI(L)-DSSI(L))*TRSIL(:,L)/(XSI(L)*MSI2-SSIL(L))
c           TRSIL(:,L)=TRSIL(:,L)-DTRSI(:,L)
#endif
            SSIL(L) = SSIL(L) - DSSI(L)
            HSIL(L) = HSIL(L) - DHSI(L)
          END IF
        END DO

      CASE ("BP")               ! Brine pocket formulation
C**** Assume basic gravity drainage. 
C**** Flushing is not included yet - but will depend on MELT12
C**** push down a fixed mass (inversely proportional to brine-frac)
C**** first layer ice
        IF (SICE.gt.0) THEN
        brine_frac=-mu*(SICE/ACE1I)/TICE
        if (brine_frac.gt.0.03d0) then ! drain
          rate = DT*BYDTSSI*100.*(brine_frac-0.03d0) ! fractional loss
          DM12 = rate*brine_frac*ACE1I
          DS12 = rate*SICE
          DH12 = rate*brine_frac*ACE1I*shw*Tice
#ifdef TRACERS_WATER
C**** tracers may be fractionated in the brine....(assume not for now)
          DTR12(:) = (DM12-DS12)*TRICE(:)/(ACE1I-SICE)
          TRICE(:) = TRICE(:)-DTR12(:)
#endif
          SICE=SICE-DS12
          HICE=HICE-DH12
        end if
        END IF

C**** check remaining layers
        DO L=3,LMI
          IF (SSIL(L).gt.0.) THEN
          TSIL(L)=Ti(HSIL(L)/(XSI(L)*MSI2),1d3*SSIL(L)/(XSI(L)*MSI2))
          brine_frac=-mu*SSIL(L)/(XSI(L)*MSI2)/TSIL(L)
          if (brine_frac.gt.0.03d0) then ! drain
            rate = DT*BYDTSSI*100.*(brine_frac-0.03d0) ! fractional loss
            DMSI(L) = rate*brine_frac*(XSI(L)*MSI2)
            DSSI(L) = rate*SSIL(L)
            DHSI(L) = rate*brine_frac*(XSI(L)*MSI2)*shw*TSIL(L)
#ifdef TRACERS_WATER
C**** tracers may be fractionated in the brine....(assume not for now)
            DTRSI(:,L)=(DMSI(L)-DSSI(L))*TRSIL(:,L)/
     *          (XSI(L)*MSI2-SSIL(L))
            TRSIL(:,L)=TRSIL(:,L)-DTRSI(:,L)
#endif
            SSIL(L) = SSIL(L) - DSSI(L)
            HSIL(L) = HSIL(L) - DHSI(L)
          end if
          END IF
        END DO
      END SELECT

C**** Calculate fluxes required to rebalance layers
C**** Mass/heat/salt moves from layer 3 to 2
      FMSI2 = -DM12
      FHSI2 = FMSI2*HSIL(3)/(XSI(3)*MSI2-DMSI(3))
      FSSI2 = FMSI2*SSIL(3)/(XSI(3)*MSI2-DMSI(3))
#ifdef TRACERS_WATER
      FTRSI2(:) = FMSI2*TRSIL(:,3)/(XSI(3)*MSI2-DMSI(3))
#endif

C**** reconstitute upper layers
      IF (DM12.gt.0) THEN
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          HSIL(1) = (HICE-FHSI2)*(ACE1I-XSI(2)*MSI1)/ACE1I + HSNOW
          SSIL(1) = (SICE-FSSI2)*(ACE1I-XSI(2)*MSI1)/ACE1I ! + SSNOW
#ifdef TRACERS_WATER
          TRSIL(:,1) = (TRICE(:)-FTRSI2(:))*(ACE1I-XSI(2)*MSI1)/ACE1I +
     *         TRSNOW(:)
#endif
        ELSE
          HSIL(1) = HSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
          SSIL(1) = 0.          ! SSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRSNOW(:)*XSI(1)*MSI1/(MSI1-ACE1I)
#endif
        END IF
        HSIL(2) = HICE - HSIL(1) - FHSI2 + HSNOW
        SSIL(2) = SICE - SSIL(1) - FSSI2 ! + SSNOW
#ifdef TRACERS_WATER
        TRSIL(:,2) = TRICE(:) - TRSIL(:,1) - FTRSI2(:) + TRSNOW(:)
#endif
      END IF

C**** Mass/heat moves between layers 3 to 4
      FMSI3 = XSI(3)*DMSI(4)+XSI(4)*(FMSI2-DMSI(3))
      IF (FMSI3.gt.0) THEN      ! downward flux to layer 4
        FHSI3 = FMSI3*HSIL(3)/(XSI(3)*MSI2-DMSI(3))
        FSSI3 = FMSI3*SSIL(3)/(XSI(3)*MSI2-DMSI(3))
#ifdef TRACERS_WATER
        FTRSI3(:) = FMSI3*TRSIL(:,3)/(XSI(3)*MSI2-DMSI(3))
#endif
      ELSE                      ! upward flux
        FHSI3 = FMSI3*HSIL(4)/(XSI(4)*MSI2-DMSI(4))
        FSSI3 = FMSI3*SSIL(4)/(XSI(4)*MSI2-DMSI(4))
#ifdef TRACERS_WATER
        FTRSI3(:) = FMSI3*TRSIL(:,4)/(XSI(4)*MSI2-DMSI(4))
#endif
      END IF

C**** Apply the second mass layer fluxes
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3
      SSIL(3)=SSIL(3)+(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)+ FSSI3
#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,3)=TRSIL(:,3)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)+ FTRSI3(:)
#endif
c     MSI1 = MSI1 - (DMSI(1)+DMSI(2)) - FMSI2 ! stays fixed
      MSI2 = MSI2 - (DMSI(3)+DMSI(4)) + FMSI2

C**** output fluxes and diagnostics
      MFLUX = DM12+SUM(DMSI(3:LMI))         ! mass flux to ocean
      HFLUX = DH12+SUM(DHSI(3:LMI))         ! energy flux to ocean
      SFLUX = DS12+SUM(DSSI(3:LMI))         ! salt flux to ocean
#ifdef TRACERS_WATER
      DO N=1,NTM
        TRFLUX(N)=DTR12(N)+SUM(DTRSI(N,3:LMI)) ! tracer flux to ocean
      END DO
#endif
C****
      RETURN
      END SUBROUTINE SSIDEC

      subroutine iceocean(Ti,Si,Tm,Sm,dh,ustar,Coriol,dtsrc,mlsh,
#ifdef TRACERS_WATER
     &     Tri,Trm,trflux,tralpha,
#endif
     &     mflux,sflux,hflux)
!@sum  iceocean calculates fluxes at base of sea ice
!@auth Gavin Schmidt
!@ver  1.0
!@usage At the ice-ocean interface the following equations hold:
!@+                                              Tb = -mu Sb         (1)
!@+ -lam_i(Ti,Si)(Ti-Tb)/dh + rho_m shw g_T (Tb-Tm) = -m Lh(Tib,Sib)
!@+                                                  + m shw (Tib-Tb)(2)
!@+                           rho_m     g_S (Sb-Sm) =  m (Sib-Sb )   (3)
!@+  with  Sib=Si, Tib=Ti m>0,  or Sib = fsss Sb, Tib=Tb m<0        (4)
!@+ This routine iteratively solves for m, Tb, Sb and Sib, Tib using
!@+ Newton's method. The form of the latent heat and conductivity
!@+ can optionally include salinity effects.
!@+ The actual fluxes (positive down) are then:
!@+  mflux = m   ; sflux = 1d-3 m Sib (kg/m^2/s)
!@+  hflux = lam_i(Ti,Si) (Ti-Tb)/dh -m Lh(Tib,Sib)+m shw Tib  (J/m^2/s)
!@+
!@+  Formulations of salinty effects PI/SI/BP from Schmidt et al (2004)
!@+
      implicit none
!@var G_mole_T,G_mole_S molecular diffusion terms for heat/salt
C****  G_mole_T = 12.5 * 13.8d0**(2d0/3d0) - 6. = 65.9d0
C****  G_mole_S = 12.5 * 2432d0**(2d0/3d0) - 6. = 2255d0
      real*8, parameter :: G_mole_T = 65.9d0 , G_mole_S = 2255d0
!@var Si,Sm salinity in lowest ice layer and mixed layer (psu)
!@var Ti,Tm temperatures in ice and mixed layer (C)
!@var dh distance from center of bottom ice layer to base of ice (m)
      real*8, intent(in) :: Ti,Si,Tm,Sm,dh
!@var coriol Coriolis parameter (used in turbulent flux calc) (1/s)
!@var ustar friction velocity at ice-ocean interface (m/s)
      real*8, intent(in) :: ustar,Coriol
!@var dtsrc source time step (s)
!@var mlsh mixed layer specific heat capactity (J/m^2 C) (UNUSED)
      real*8, intent(in) :: dtsrc,mlsh
!@var mflux,sflux,hflux mass, salt and heat fluxes at base of ice
      real*8, intent(out) :: mflux,sflux,hflux
#ifdef TRACERS_WATER
!@var Tri,Trm tracer concentration in ice and mixed layer (kg/kg)
!@var tralpha tracer fraction going into ice (1)
      real*8, dimension(ntm), intent(in) :: Trm,Tri,tralpha
!@var trflux tracer mass flux at base
      real*8, dimension(ntm), intent(out) :: trflux
!@var Trib tracer concentration at interface ice (kg/kg)
      real*8, dimension(ntm) :: Trib
#endif

!@var g_T,g_S turbulent exchange velocities (m/s)
!@var G_turb turbulent diffusion term
!@var Sb,Sb0 final and initial estimate for salinity at interface (psu)
      real*8 :: G_turb,g_T,g_S,Sb,Sb0
!@var dSbdSb differential of Sb with respect to initial Sb0
      real*8 :: dSbdSb
!@var m melt rate (negative implies freezing) (kg/m^2/s)
!@var Tb temperature at interface (C)
!@var Tib,Sib actual temperature/salinity in melting/freezing ice (psu)
      real*8 :: m,Tb,Sib  ! ,Tib
      real*8 :: lh,left2,df3dm,dmdTb,dmdSi,alamdh,f0,df,rsg
      integer, parameter :: niter=5 !@param niter number of iterations
      integer :: i

C**** Calculate turbulent exchange velocities g_T, g_S

C**** G_turb = (1/k) ln (u* E n^2 / f h) + 1 / (2 E n) - (1/k)
C****        = 2.5 ln ( 5300*(u*)^2/coriol ) + 9.62 - 2.5 (assuming n=1)
C**** Note: n should theoretically depend on the buoyancy flux and
C****       therfore this should be included in the iteration below.
      G_turb = 2.5d0 * log ( 5300.*ustar*ustar/coriol ) + 7.12d0
C****  g  = u* / ( G_turb + G_mole )
      g_T = ustar / ( G_turb + G_mole_T )
      g_S = ustar / ( G_turb + G_mole_S )

C**** set conductivity term
C**** Diffusive flux is implicit in ice !  + ml temperature (not used)
      rsg=rhows*shw*g_T    ! /(1.+alpha*dtsrc*rhows*shw*g_T/mlsh)

      if (form.eq."SI" .or. Si.eq.0.) then
                                ! no salinity effects on diffusion
        alamdh = alami/(dh+alpha*dtsrc*alami*byshi/(2.*dh*rhoi))
      else                      ! brine pocket impact on diffusion
        alamdh = (alami+alamdS*Si/Ti)/(dh+alpha*dtsrc*
     *           (alami+alamdS*Si/Ti)*byshi/(2.*rhoi*dh))
      end if 

C**** solve for boundary values (uses Newton-Rapheson with bounds)
C**** Estimate initial boundary salinity
      Sb0 = 0.75d0*Sm+0.25d0*Si
      do i=1,niter
C**** we use the full tfrez calc here, but assume for simplicity that
C**** the gradient w.r.t. T is still = -mu below.
        Tb = tfrez(Sb0)   ! -mu*Sb0  ! freezing point at salinity Sb0
C**** calculate left hand side of equation 2
        left2 = -alamdh*(Ti-Tb) + rsg*(Tb-Tm)
C**** depending on whether there is freezing or melting, the sea ice
C**** salinity in the melt/freeze fraction is set to be the original
C**** ice salinity or is a function of the boundary salinity.
C**** Similarly for temperature Tib
        if (left2.gt.0) then    ! freezing
          if (qsfix) then   ! keep salinity in ice constant
            Sib = ssi0*1d3
          else              ! it is a function of boundary value
            Sib = Sb0*fsss
          end if

          select case (form)
c          case ("PI")           ! pure ice
c            lh = lhm + Tb*(shw-shi)
          case ("SI")           ! salinity effects only mass
            lh = lhm*(1.-Sib*1d-3) + Tb*(shw-shi)
          case ("BP")           ! brine pockets 
            if (Sib.gt.0) then
              lh = lhm*(1.+mu*Sib/Tb) + (Tb+mu*Sib)*(shw-shi)
            else
              lh = lhm + Tb*(shw-shi)
            end if
          end select
c
          m = -left2/lh

          select case (form)
c          case ("PI")           ! pure ice
c            dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
c            dmdSi = 0.
          case ("SI")           ! salinity effects only mass
            dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
            dmdSi = -left2*lhm*1d-3/(lh*lh)
          case ("BP")           ! brine pockets 
            if (Sib.gt.0) then
              dmdTb = (left2*(-lhm*mu*Sib/Tb**2+shw-shi)-lh*(alamdh+rsg)
     *             )/(lh*lh)
              dmdSi = (left2*mu*(lhm/Tb+shw-shi))/(lh*lh)
            else
              dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
              dmdSi = 0.
            end if
          end select

          if (qsfix) then   ! keep salinity in ice constant
            Sb = (rhows*g_S*Sm+m*Sib)/(rhows*g_S+m)
            df3dm= rhows*g_S*(Sib-Sm)/(rhows*g_S+m)**2
            dSbdSb= -mu*dmdTb*df3dm
          else              ! it is a function of boundary value
            Sb = rhows*g_S*Sm/(rhows*g_S+m*(1.-fsss))
            df3dm= -rhows*g_S*Sm*(1.-fsss)/(rhows*g_S+m*(1.-fsss))**2
            dSbdSb= (fsss*dmdSi-mu*dmdTb)*df3dm
          end if
        else                    ! melting
          Sib = Si

          select case (form)
c          case ("PI")           ! pure ice
c            lh = lhm + Tb*shw - Ti*shi
          case ("SI")           ! salinity effects only mass
            lh = lhm*(1.-Sib*1d-3) + Tb*shw - Ti*shi
          case ("BP")           ! brine pockets 
            if (Sib.gt.0) then
              lh = lhm*(1.+mu*Sib/Ti)+(Ti+mu*Sib)*(shw-shi)-shw*(Ti-Tb)
            else
              lh = lhm + Tb*shw - Ti*shi
            end if
          end select
c
          m = -left2/lh

          Sb = (m*Sib+rhows*g_S*Sm)/(rhows*g_S+m)
          df3dm= (Sib*(rhows*g_S+m)-(m*Sib+rhows*g_S*Sm))/(rhows*g_S+m)
     *         **2
          dmdTb = -(alamdh + rsg)/lh + left2*shw/lh**2
          dSbdSb= -mu*dmdTb*df3dm
        end if

        f0=Sb-Sb0
        df=dSbdSb-1.+1d-20
        Sb = min(max(Sib,Sb0 - f0/df),40d0)
        Sb0=Sb
      end do
C**** define fluxes (positive down)
C**** Cap mass flux at at 90% of bottom layer
      if (m.gt.0.9d0*2.*dh*rhoi/dtsrc) then
        m=0.9d0*2.*dh*rhoi/dtsrc
      end if

      mflux = m                       ! (kg/m^2 s)
      sflux = 1d-3*m*Sib              ! (kg/m^2 s)
      hflux = alamdh*(Ti-Tb) - m*lh +m*shw*Tb ! (J/m^2 s)
#ifdef TRACERS_WATER
C**** Tracers use salinity turbulent diffusion term
      if (m.gt.0) then
        Trib(:)=Tri(:)
      else
        Trib(:)=tralpha(:)*rhows*g_S*Trm(:)/(rhows*g_S+(m-sflux)*(1.
     *       -tralpha(:)))
      end if
      trflux(:) = (m - sflux) * Trib(:)  ! (kg/m^2 s)
#endif
C****
      return
      end subroutine iceocean

      subroutine icelake(Ti,Tm,dh,dtsrc,mlsh,
#ifdef TRACERS_WATER
     &     Tri,Trm,trflux,tralpha,
#endif
     &     mflux,hflux)
!@sum  icelake calculates fluxes at base of lake ice (no salinity)
!@auth Gavin Schmidt
!@ver  1.0
!@usage At the ice-lake interface the following equations hold:
!@+   interface is at freezing point (Tb=0.)   (1)
!@+   -lam_i Ti/dh - rho_m shw g_T Tm = -m Lh(Tib) + m shw Tib   (2)
!@+     with  Tib=Ti m>0,  or Tib=0. m<0        (4)
      implicit none
!@var mflux,hflux mass and heat fluxes at base of ice
      real*8, intent(out) :: mflux,hflux
!@var dtsrc source time step (s)
!@var Ti,Tm temperatures in ice and upper lake layer (C)
!@var dh distance from center of bottom ice layer to base of ice (m)
!@var mlsh mixed layer specific heat capactity (J/m^2 C) (UNUSED)
      real*8, intent(in) :: Ti,Tm,dh,dtsrc,mlsh
C**** Assume constant g_T = 1.3d-7, g_S = 0.025 * g_T m/s
!@var rsg = rhow * shw * g_T turbulent energy flux (J/m^2 s)
!@var rgS = rhow * g_S turbulent tracer flux (kg/m^2 s)
      real*8, parameter ::  rsg = rhow*shw*1.3d-7, rgS=rhow*3.2d-9
      real*8 left2, lh, m, alamdh
#ifdef TRACERS_WATER
!@var Tri,Trm tracer concentration in ice and upper lake layer (kg/kg)
!@var tralpha tracer fraction going into ice (1)
      real*8, dimension(ntm), intent(in) :: Trm,Tri,tralpha
!@var trflux tracer mass flux at base
      real*8, dimension(ntm), intent(out) :: trflux
!@var Trib tracer concentration at interface ice (kg/kg)
      real*8, dimension(ntm) :: Trib
#endif

C**** Diffusive flux is implicit in ice temperature
      alamdh = alami/(dh+alpha*dtsrc*byshi*alami/(2.*dh*rhoi))
C**** calculate left hand side of equation 2
      left2 = -alamdh*Ti - rsg*Tm    !/(1.+alpha*dtsrc*rsg/mlsh)
      if (left2.gt.0) then      ! freezing
        lh = lhm
      else                      ! melting
        lh = lhm -  Ti*shi
      end if
C**** define fluxes (positive down)
      m = -left2/lh
C**** Cap mass flux at 90% of bottom layer
      if (m.gt.0.9d0*2.*dh*rhoi/dtsrc) then
        m=0.9d0*2.*dh*rhoi/dtsrc
      end if
      mflux = m                 ! (kg/m^2 s)
      hflux = alamdh*Ti - m*lh  ! (J/m^2 s)
#ifdef TRACERS_WATER
C**** Tracers use tracer turbulent diffusion term
      if (m.gt.0) then
        Trib(:) = Tri(:)
      else
        Trib(:) = tralpha(:)*rgS*Trm(:)/(rgS+m*(1.-tralpha(:)))
      end if
      trflux(:) = m*Trib(:)     ! (kg/m^2 s)
#endif
C****
      return
      end subroutine icelake

      subroutine solar_ice_frac(snow,msi2,wetsnow,fsri,lmax)
!@sum solar_ice_frac calculates the fraction of solar radiation
!@+   transmitted through the ice, taking albedo into account, but
!@+   assuming constant proportions of visible and near-ir.
!@+   Based on Ebert et al, 1995.
!@auth Gavin Schmidt
      implicit none
!@var kiextvis,kiextnir1 ice extinction coeff. for vis and nir (1/m)
      real*8, parameter :: kiextvis=1.5d0, kiextnir1=18.d0
!@var snow,msi2 snow and second layer ice amounts (kg/m^2)
      real*8, intent(in) :: snow,msi2
!@var wetsnow =true if snow is wet
      logical, intent(in) :: wetsnow
!@var lmax number of ice layers to do calculation for
      integer, intent(in) :: lmax
!@var fsri fraction of solar energy reaching bottom of layer
      real*8, dimension(lmax), intent(out) :: fsri
      real*8, dimension(lmax) :: fsrvis, fsrnir1, hice
!@var ksextvis,ksextnir1 snow extinction coeff. for vis and nir (1/m)
!@var fracvis,fracnir fraction of absorbed solar in vis and nir
      real*8 fracvis,fracnir1,ksextvis,ksextnir1
      real*8 hsnow,hsnow1,hice12
      integer l

      hsnow = snow/rhos
      hice12= ace1i*byrhoi

C**** old code
c!@param KIEXT extinction coefficient for light in sea ice (1/m)
c      REAL*8, PARAMETER :: KIEXT = 1.5d0
c!@param KSEXT extinction coefficient for light in snow (1/m)
c      REAL*8, PARAMETER :: KSEXT = 15d0
c        IF (ACE1I*XSI(1).gt.SNOW*XSI(2)) THEN
cC**** first thermal layer is snow and ice
c          HICE1 = (ACE1I-XSI(2)*(SNOW+ACE1I))*BYRHOI
c          FSRI(1) = EXP(-KSEXT*HSNOW-KIEXT*HICE1)
c        ELSE ! all snow
c          HSNOW1 = (ACE1I+SNOW)*XSI(1)/RHOS
c          FSRI(1) = EXP(-KSEXT*HSNOW1)
c        END IF
c        FSRI(2) = EXP(-KSEXT*HSNOW-KIEXT*HICE12)
c        FSRI(3) = FSRI(2)*EXP(-KIEXT*MSI2*XSI(3)*BYRHOI)
c        FSRI(4) = FSRI(2)*EXP(-KIEXT*MSI2*BYRHOI)

c**** Set fractions of visible and near-ir bands based on weighting the
c**** solar input by the approximate co-albedo (ignoring melt ponds etc)
c**** VIS:  250 - 690 nm,       49.3% of incoming solar at ground
c**** NIR1: 690 - 1190 nm       34.9%       "       "   "    "
c**** NIR2/3 (> 1190 nm assumed not to be transmitted)
c****
      if (hsnow.gt.0.02) then    ! same cutoff as for albedo
c**** Extinction coefficients for snow depend on wetness
c**** Since albedo does as well, absorbed fractions also vary
        if (wetsnow) then
          fracvis =0.20d0 ;  fracnir1=0.33d0
          ksextvis=10.7d0 ; ksextnir1=118.d0
        else
          fracvis =0.06d0 ;  fracnir1=0.31d0
          ksextvis=19.6d0 ; ksextnir1=196.d0
        end if
      else
        fracvis=0.24d0 ; fracnir1=0.43d0
        ksextvis=10.7d0 ; ksextnir1=118.d0
      end if

C**** calculate sucessive fractions assuming each band is attenuated
C**** independently.
      if (ace1i*xsi(1).gt.snow*xsi(2)) then
c**** First thermal layer is snow and ice
        hice(1)    = (ace1i-xsi(2)*(snow+ace1i))*byrhoi
        fsrvis (1) = exp(-ksextvis *hsnow - kiextvis *hice(1))
        fsrnir1(1) = exp(-ksextnir1*hsnow - kiextnir1*hice(1))
        fsri(1)    = fracvis*fsrvis(1) + fracnir1*fsrnir1(1)
      else                      ! all snow
        hsnow1     = (ace1i+snow)*xsi(1)/rhos
        fsrvis (1) = exp(-ksextvis *hsnow1)
        fsrnir1(1) = exp(-ksextnir1*hsnow1)
        fsri(1)    = fracvis*fsrvis(1) + fracnir1*fsrnir1(1)
      end if
      fsrvis (2) = exp(-ksextvis *hsnow - kiextvis *hice12)
      fsrnir1(2) = exp(-ksextnir1*hsnow - kiextnir1*hice12)
      fsri(2)    = fracvis*fsrvis(2) + fracnir1*fsrnir1(2)

      do l=3,lmax   ! only done if lmax>2
        hice(l)    = xsi(l)*msi2*byrhoi
        fsrvis (l) = exp(-kiextvis *hice(l)) * fsrvis (l-1)
        fsrnir1(l) = exp(-kiextnir1*hice(l)) * fsrnir1(l-1)
        fsri(l)    = fracvis*fsrvis(l) + fracnir1*fsrnir1(l)
      end do
c****
      return
      end subroutine solar_ice_frac

      function tfrez(sss,press)
!@sum tfrez calculates freezing temperature of sea water
!@auth Gavin Schmidt
      implicit none
!@var sss sea surface salinity (psu)
      real*8, intent(in) :: sss
!@var press gauge pressure (default=0) (Pa)
      real*8, intent(in), optional :: press
!@var tfrez approx. freezing point of sea water (C)
      real*8 tfrez,pr
      real*8, parameter :: dtdp = -7.5d-8
      real*8 :: a01 = -.0575d0, a02 = -2.154996d-4, a03 =1.710523d-3

      pr=0.
      if (present(press)) then
        pr=press
      end if
C**** linear approximation
c      tfrez = -mu*sss + dtdp*pr
C**** UNESCO formula (1983)
      tfrez = (a01 + a02*sss)*sss + a03*sss*sqrt(sss)+ dtdp*pr

      return
      end function tfrez

      subroutine snowice(Tm,Sm,SNOW,MSI2,HSIL,SSIL,qsfix,
#ifdef TRACERS_WATER
     *         Trm,Tralpha,TRSIL,TRSNWIC,
#endif
     *         MSNWIC,HSNWIC,SSNWIC)
      IMPLICIT NONE
!@var Tm ocean mixed layer temperature (C)
!@var Sm ocean mixed layer salinity (psu)
      REAL*8, INTENT(IN) :: Tm, Sm
!@var MSI1, MSI2 snow and ice mass in two mass layers (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW,MSI2
!@var HSIL, SSIL energy and salt profiles in ice (J, kg)/m^2
      REAL*8, DIMENSION(LMI), INTENT(INOUT) :: HSIL,SSIL
!@var MSNWIC, HSNWIC, SSNWIC mass, energy and salt flux to ocean (<0) 
      REAL*8, INTENT(OUT) :: MSNWIC, HSNWIC, SSNWIC
!@var qsfix true if sea ice salinity is fixed
      LOGICAL, INTENT(IN) :: qsfix
#ifdef TRACERS_WATER
!@var Trm tracer concentration in mixed layer
!@var Tralpha fractionation for snow ice
C**** Be careful to avoid taking too much tracer from ocean box
      REAL*8, DIMENSION(NTM), INTENT(IN) :: Trm, Tralpha
!@var TRSIL tracer profile in ice (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TRSNWIC tracer flux to ocean (<0) 
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRSNWIC
!@var Tri tracer concentration in snow ice
      REAL*8, DIMENSION(NTM) :: Tri, TRICE, TRSNOW, FTRSI2, FTRSI3
#endif
      REAL*8 DSNOW,Z0,MAXM,MAXME,Eoc,Esnow,Eic,Erat,Si,MSI1,Tf
      REAL*8 SICE,HICE,HSNOW,FMSI2,FMSI3,FHSI2,FSSI2,FHSI3,FSSI3

C**** test for snow ice possibility      
      IF (RHOI*SNOW.gt.(ACE1I+MSI2)*(RHOWS-RHOI)) THEN
        MSI1=SNOW+ACE1I
        Z0=(MSI1+MSI2)/RHOWS-(ACE1I+MSI2)/RHOI
        MAXM=Z0*RHOWS*(RHOI-RHOS)/(RHOWS+RHOS-RHOI)

C**** heat, salt, tracer amounts in snow, first layer ice
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          HSNOW = SNOW*(Ti(HSIL(1)/(XSI(1)*MSI1),1d3*SSIL(1)
     *         /(XSI(1)*MSI1))*shi-lhm)
        ELSE
          HSNOW = HSIL(1)+(XSI(2)*MSI1-ACE1I)*(Ti(HSIL(2)
     *         /(XSI(2)*MSI1),1d3*SSIL(2)/(XSI(2)*MSI1))*shi-lhm)
        END IF
        HICE  = HSIL(1)+HSIL(2)-HSNOW
        SICE  = SSIL(1)+SSIL(2)
#ifdef TRACERS_WATER
        TRSNOW(:) = TRSIL(:,1)*MIN((MSI1-ACE1I)/(XSI(1)*MSI1-SSIL(1))
     *       ,1d0)+TRSIL(:,2)*MAX(1.-(ACE1I-SICE)/(XSI(2)*MSI1-SSIL(2))
     *       ,0d0)
        TRICE(:)  = TRSIL(:,1)+TRSIL(:,2)-TRSNOW(:)
#endif

C**** Check energy constraint on mass flux
C**** Si is the salinity only in the frozen sea water part
        if (qsfix) then
          Si=1d3*ssi0
        else
          Si=fsss*Sm
        end if
        Tf=Tfrez(Sm)
        Eic=Ei(Tf,Si) 
        Esnow=HSNOW/SNOW
        Eoc=Tm*SHW
C**** There is a small error in Schmidt et al (2004) equation (12):
C**** (Esnow-Ei(Tm,Si))/(Ei(Tm,Si)-Eoc) should be 
C**** (Esnow-Ei(Tm, 0))/(Ei(Tm,Si)-Eoc)
        ERAT=(Esnow-Ei(Tf,0d0))/(Eic-Eoc)

        MAXME=Z0*RHOI*ERAT/(1.+ERAT*(RHOWS-RHOI)/RHOWS)
        MAXM=MAX(0d0,MIN(MAXME,MAXM))   ! mass of sea water to be added

        DSNOW=Z0*RHOI-MAXM*(RHOWS-RHOI)/RHOWS  ! loss of snow

        FMSI2 = MAXM+DSNOW      ! mass of ice pushed down to layer 3

C**** Add new ice 
        SICE=SICE+MAXM*0.001d0*Si
        HICE=HICE+MAXM*Eoc+DSNOW*Esnow
#ifdef TRACERS_WATER
C**** need to be careful with fractionation for tracers
        Tri(:)=Tralpha(:)*Trm(:)
        TRICE(:)=TRICE(:)+MAXM*(1.-0.001d0*Si)*Tri(:)+DSNOW*TRSNOW(:)
     *       /SNOW
#endif
        FSSI2 = FMSI2*SICE/(ACE1I+FMSI2)
        FHSI2 = FMSI2*HICE/(ACE1I+FMSI2)
C**** Update ice and snow concentrations
#ifdef TRACERS_WATER
        FTRSI2(:)=FMSI2*TRICE(:)/(ACE1I+FMSI2)
        TRICE(:)=TRICE(:)-FTRSI2(:)
        TRSNOW(:)=TRSNOW(:)*(SNOW-DSNOW)/SNOW
#endif
        SICE=SICE-FSSI2
        HICE=HICE-FHSI2
        HSNOW=HSNOW*(SNOW-DSNOW)/SNOW   ! snow temp remains const
        SNOW=SNOW-DSNOW
        MSI1=ACE1I+SNOW

C**** Reconstitute upper layers
        IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
          HSIL(1) = HICE*(ACE1I-XSI(2)*MSI1)/ACE1I + HSNOW
          SSIL(1) = SICE*(ACE1I-XSI(2)*MSI1)/ACE1I ! + SSNOW
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRICE(:)*(ACE1I-XSI(2)*MSI1)/ACE1I + TRSNOW(:)
#endif
        ELSE
          HSIL(1) = HSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
          SSIL(1) = 0.          ! SSNOW*XSI(1)*MSI1/(MSI1-ACE1I)
#ifdef TRACERS_WATER
          TRSIL(:,1) = TRSNOW(:)*XSI(1)*MSI1/(MSI1-ACE1I)
#endif
        END IF
        HSIL(2) = HICE - HSIL(1) + HSNOW
        SSIL(2) = SICE - SSIL(1) ! + SSNOW
#ifdef TRACERS_WATER
        TRSIL(:,2) = TRICE(:) - TRSIL(:,1) + TRSNOW(:)
#endif
C**** lower level fluxes
        FMSI3 = FMSI2*XSI(4)    ! mass of ice pushed down to layer 4
        FSSI3 = FMSI3*(SSIL(3)+FSSI2)/(XSI(3)*MSI2+FMSI2)
        FHSI3 = FMSI3*(HSIL(3)+FHSI2)/(XSI(3)*MSI2+FMSI2)
#ifdef TRACERS_WATER
        FTRSI3(:)=FMSI3*(TRSIL(:,3)+FTRSI2(:))/(XSI(3)*MSI2+FMSI2)
#endif

C**** Update lower layers
        MSI2=MSI2+FMSI2
        SSIL(3)=SSIL(3)+FSSI2-FSSI3
        SSIL(4)=SSIL(4)      +FSSI3
        HSIL(3)=HSIL(3)+FHSI2-FHSI3
        HSIL(4)=HSIL(4)      +FHSI3
#ifdef TRACERS_WATER
        TRSIL(:,3) = TRSIL(:,3) + FTRSI2(:) - FTRSI3(:)
        TRSIL(:,4) = TRSIL(:,4)             + FTRSI3(:)
#endif

C**** output flux (positive down)
        MSNWIC = -MAXM
        HSNWIC = MSNWIC*Eoc
        SSNWIC = 0.001d0*Si*MSNWIC
#ifdef TRACERS_WATER
        TRSNWIC(:) = Tri(:)*MSNWIC*(1.-0.001d0*Si)
#endif
      ELSE
        MSNWIC = 0. ; HSNWIC = 0. ; SSNWIC = 0.
#ifdef TRACERS_WATER
        TRSNWIC = 0.
#endif
      END IF
      RETURN
      end subroutine snowice

      SUBROUTINE TICE(HSIL,SSIL,MSI1,MSI2,TSIL)
!@sum TICE returns array of ice temperatures from model variables
!@auth Gavin Schmidt
      IMPLICIT NONE
      REAL*8, DIMENSION(LMI), INTENT(IN) :: HSIL, SSIL  ! J/m2, kg/m2
      REAL*8, INTENT(IN) :: MSI1,MSI2  ! kg/m2
      REAL*8, DIMENSION(LMI), INTENT(OUT) :: TSIL  ! deg C
      INTEGER l

      do l=1,2
        TSIL(l)=Ti(HSIL(l)/(XSI(l)*MSI1),1d3*SSIL(l)/(XSI(l)*MSI1))
      end do
      do l=3,lmi
        TSIL(l)=Ti(HSIL(l)/(XSI(l)*MSI2),1d3*SSIL(l)/(XSI(l)*MSI2))
      end do

      RETURN
      END SUBROUTINE TICE

      REAL*8 FUNCTION Ti(Ei,Si)
!@sum Ti calulcates sea ice temperature as a function of internal 
!@+   energy and salinity content depending on ice thermo formulation
!@auth Jiping Liu/Gavin Schmidt
      USE CONSTANT, only : lhm,shw,shi,byshi
      IMPLICIT NONE
!@var Ei internal energy of sea ice (J/kg)
!@var Si salinity  of sea ice (10-3 kg/kg = ppt)
!@var TI temperature (deg C)
      REAL*8, INTENT(IN) :: Ei, Si
      real*8 b,c,det,tm

      select case (form)
c      case ("PI")               ! pure ice 
cc**** solve Ei = shi*Ti-lhm
c        Ti=(Ei+lhm)*byshi
      case ("SI")               ! salinity mass effect
c**** solve Ei = shi*Ti-lhm*(1-1d-3*Si)
        Ti=(Ei+lhm*(1.-1d-3*Si))*byshi
      case ("BP")               ! Brine pocket formulation
c**** solve Ei = shi*(Ti+mu*Si)-lhm*(1+mu*Si/Ti)-mu*Si*shw
        Tm=-mu*Si
        if(Ei .ge. shw*Tm) then ! at or above melting point
          Ti=Tm
        else                    ! solve quadratic, take most negative root
          b=Tm*(shw-shi)-Ei-lhm
          c=lhm*Tm
          det=b*b-4d0*shi*c      ! > 0
          Ti=-5d-1*(b + sqrt(det))*byshi
        endif
      end select

      return
      END FUNCTION Ti

      REAL*8 FUNCTION Ei(Ti,Si)
!@sum Calculation of energy of ice (J/kg) as a function of temp and salinity
!@auth Gavin Schmidt
      USE CONSTANT, only : shi, lhm, shw
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Ti,Si   ! deg C and psu

      select case (form)
c      case ("PI")               ! no salinity effect
c        Ei=Ti*shi-lhm
      case ("SI")               ! salinity affects only mass
        Ei=Ti*shi-lhm*(1.-1d-3*Si)
      case ("BP")               ! brine pocket formulation
        if (Si.gt.0) then  ! is this safe from T=0? (or T>-muS?)
          Ei= (Ti+mu*Si)*shi-lhm*(1.+mu*Si/Ti)-shw*mu*Si
        else
          Ei= Ti*shi-lhm
        end if
      end select

      RETURN
      END FUNCTION Ei

      REAL*8 FUNCTION Mi(HSI,SSI,MSI)
!@sum Calculation of ice melt (kg/m2) as a function of hsi and salinity
!@auth Gavin Schmidt
      USE CONSTANT, only : lhm, shw, bylhm
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HSI,SSI,MSI   ! J/m2, kg/m2, kg/m2

      select case (form)
c      case ("PI")               ! no salinity effect
c        Mi=max(0d0,msi+hsi*bylhm)
      case ("SI")               ! salinity affects only mass
        Mi=max(0d0,msi+hsi*bylhm/(1.-ssi/msi))
      case ("BP")               ! brine pocket formulation
        Mi=0.
        if (hsi+shw*mu*1d3*ssi.gt.0) Mi=msi
      end select

      RETURN
      END FUNCTION Mi

      END MODULE SEAICE

      MODULE SEAICE_COM
!@sum  SEAICE_COM contains the model arrays for seaice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP, only : grid
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      USE SEAICE, only : lmi

      IMPLICIT NONE
!@var RSI fraction of open water area covered in ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RSI
!@var SNOWI snow amount on sea ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SNOWI
!@var MSI mass of ice second layer (layer 1=const) (kg/m^2)
C**** Note that MSI includes the mass of salt in sea ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MSI
!@var HSI enthalpy of each ice layer (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: HSI
!@var SSI sea ice salt content (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SSI
!@var pond_melt amount of melt pond mass (kg/m^2)
C**** Note this is a virtual meltpond and is only used for
C**** albedo calculations
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pond_melt
!@var flag_dsws true if snow on ice is wet (ie. rain or surface melt)
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: flag_dsws

#ifdef TRACERS_WATER
!@var TRSI tracer amount in sea ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSI
!@var TRSI0 default tracer conc. in sea ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: TRSI0
#endif

      END MODULE SEAICE_COM

      SUBROUTINE ALLOC_SEAICE_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at 
!@+    run-time
!@auth Rodger Abel
!@ver  1.0
      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE DOMAIN_DECOMP, ONLY : GET
      USE MODEL_COM, ONLY : IM, JM
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
#endif
      USE SEAICE, only : LMI

      USE SEAICE_COM, ONLY : RSI, SNOWI, MSI, HSI, SSI, pond_melt, 
     *     flag_dsws
#ifdef TRACERS_WATER
      USE SEAICE_COM, ONLY : TRSI, TRSI0
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE( RSI(IM, J_0H:J_1H),
     *     SNOWI(IM, J_0H:J_1H),
     *     MSI(IM, J_0H:J_1H),
     *     pond_melt(IM, J_0H:J_1H),
     *     flag_dsws(IM, J_0H:J_1H),
     *     STAT=IER)

!hack hack hack !!!!!!!! 
      rsi = 0.d0   !!!!  call to io_seaice may be missing during postproc.
      SNOWI = 0.d0
      MSI = 0.d0
      pond_melt = 0.d0
      flag_dsws  = .false.

      ALLOCATE( HSI(LMI, IM, J_0H:J_1H),
     *     SSI(LMI, IM, J_0H:J_1H),
     *     STAT=IER)

      hsi = 0.d0
      ssi = 0.d0

#ifdef TRACERS_WATER
      ALLOCATE( TRSI(NTM, LMI, IM, J_0H:J_1H),
     *     TRSI0(NTM),
     *     STAT=IER)
      TRSI(:, :, :, J_0H:J_1H) = 0.  ! default to prevent unecessary crash
#endif

      RETURN
      END SUBROUTINE ALLOC_SEAICE_COM

      SUBROUTINE io_seaice(kunit,iaction,ioerr)
!@sum  io_seaice reads and writes seaice variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irsfic,irsficno,irerun
      USE SEAICE_COM
      USE DOMAIN_DECOMP, only : GRID, GET, AM_I_ROOT
      USE DOMAIN_DECOMP, only : PACK_COLUMN, PACK_DATA, PACK_BLOCK
      USE DOMAIN_DECOMP, only : UNPACK_COLUMN, UNPACK_DATA
      USE DOMAIN_DECOMP, only :  UNPACK_BLOCK
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SICE02"
      REAL*8, DIMENSION(IM,JM) :: RSI_GLOB, SNOWI_GLOB
      REAL*8, DIMENSION(IM,JM) :: MSI_GLOB, POND_MELT_GLOB
      LOGICAL, DIMENSION(IM,JM) :: FLAG_DSWS_GLOB
      REAL*8 :: HSI_GLOB(LMI,IM,JM), SSI_GLOB(LMI,IM,JM)
      INTEGER :: J_0,J_1
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSICE01"
      REAL*8 ::  TRSI_GLOB(NTM, LMI, IM, JM)

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a8,i3,a1,i3,a)')'R8 TRSI(',ntm,',',lmi,',im,jm)'
#endif

      write (MODULE_HEADER(lhead+1:80),'(a14,i1,a20,i1,a23)')
     * 'R8 F(im,jm),H(',lmi,',im,jm),snw,msi,ssi(',lmi,
     *     '),pond_melt,L flag_dsws'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
         CALL PACK_DATA(grid, RSI, RSI_GLOB)
         CALL PACK_DATA(grid, SNOWI, SNOWI_GLOB)
         CALL PACK_DATA(grid, MSI, MSI_GLOB)
         CALL PACK_DATA(grid, POND_MELT, POND_MELT_GLOB)
         CALL PACK_DATA(grid, FLAG_DSWS, FLAG_DSWS_GLOB)
         CALL PACK_COLUMN(grid, HSI, HSI_GLOB)
         CALL PACK_COLUMN(grid, SSI, SSI_GLOB)
#ifdef TRACERS_WATER
         CALL PACK_BLOCK(grid, TRSI, TRSI_GLOB)
#endif
         IF (AM_I_ROOT()) THEN
           WRITE (kunit,err=10) MODULE_HEADER, RSI_glob, HSI_glob,
     *       SNOWI_glob, MSI_glob, SSI_glob,POND_MELT_glob,
     *       FLAG_DSWS_glob
#ifdef TRACERS_WATER

           WRITE (kunit,err=10) TRMODULE_HEADER,TRSI_glob
#endif
        END IF
      CASE (IOREAD:)            ! input from restart file
         if ( AM_I_ROOT() ) then
           READ (kunit,err=10) HEADER,RSI_GLOB,HSI_GLOB,
     *        SNOWI_GLOB,MSI_GLOB,SSI_GLOB
     *       ,POND_MELT_GLOB,FLAG_DSWS_GLOB
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if

        CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

         CALL UNPACK_DATA( grid, RSI_GLOB,       RSI )
         CALL UNPACK_DATA( grid, SNOWI_GLOB,     SNOWI )
         CALL UNPACK_DATA( grid, MSI_GLOB,       MSI )
         CALL UNPACK_DATA( grid, POND_MELT_GLOB, POND_MELT )
         CALL UNPACK_DATA( grid, FLAG_DSWS_GLOB, FLAG_DSWS )

         CALL UNPACK_COLUMN( grid, HSI_GLOB, HSI )
         CALL UNPACK_COLUMN( grid, SSI_GLOB, SSI )

#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO) ! reruns/restarts
          if  (AM_I_ROOT() ) then
              READ (kunit,err=10) TRHEADER,TRSI_GLOB
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
          end if
         CALL UNPACK_BLOCK(grid, TRSI_GLOB, TRSI)
        END SELECT
#endif
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_seaice

      SUBROUTINE CHECKI(SUBR)
!@sum  CHECKI Checks whether Ice values are reasonable
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM
      USE GEOM, only : imaxj
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE SEAICE, only : lmi,xsi,ace1i,Ti
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
#ifdef TRACERS_WATER
     *     ,trsi
#endif
      USE LAKES_COM, only : flake
      USE FLUXES
      USE DOMAIN_DECOMP, only : GRID
      USE DOMAIN_DECOMP, only : GET
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var QCHECKI true if errors found in seaice
      LOGICAL QCHECKI
      INTEGER I,J,L
      REAL*8 TICE
#ifdef TRACERS_WATER
      integer :: imax,jmax, n
      real*8 relerr,errmax
#endif

      integer :: J_0, J_1, J_0H, J_1H
C**** 
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1, 
     *     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

C**** Check for NaN/INF in ice data
      CALL CHECK3B(RSI,IM,J_0H,J_1H,JM,1,SUBR,'rsi   ')
      CALL CHECK3B(MSI,IM,J_0H,J_1H,JM,1,SUBR,'msi   ')
      CALL CHECK3C(HSI,LMI,IM,J_0H,J_1H,SUBR,'hsi   ')
      CALL CHECK3C(SSI,LMI,IM,J_0H,J_1H,SUBR,'ssi   ')
      CALL CHECK3B(SNOWI,IM,J_0H,J_1H,JM,1,SUBR,'sni   ')

      QCHECKI = .FALSE.
C**** Check for reasonable values for ice variables
      DO J=J_0, J_1
        DO I=1,IMAXJ(J)
          IF (RSI(I,J).lt.0 .or. RSI(I,j).gt.1 .or. MSI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,RSI,MSI=',I,J,RSI(I,J)
     *           ,MSI(I,J)
            QCHECKI = .TRUE.
          END IF
          IF ( (FOCEAN(I,J)+FLAKE(I,J))*RSI(I,J).gt.0) THEN
          DO L=1,LMI
            IF (L.le.2) TICE = Ti(HSI(L,I,J)/(XSI(L)*(ACE1I+SNOWI(I,J)))
     *           ,1d3*SSI(L,I,J)/(XSI(L)*(ACE1I+SNOWI(I,J))))
            IF (L.gt.2) TICE = Ti(HSI(L,I,J)/(XSI(L)*MSI(I,J))
     *           ,1d3*SSI(L,I,J)/(XSI(L)*MSI(I,J)))
            IF (HSI(L,I,J).gt.0.or.TICE.gt.1d-10.or.TICE.lt.-80.) THEN
              WRITE(6,'(3a,3i3,6e12.4/1X,6e12.4)')
     *             'After ',SUBR,': I,J,L,TSI=',I,J,L,TICE,RSI(I,J)
              WRITE(6,*) HSI(:,I,J),MSI(I,J),SNOWI(I,J),SSI(:,I,J)
              IF (TICE.gt.1d-3.or.TICE.lt.-100.) QCHECKI = .TRUE.
            END IF
            IF (SSI(L,I,J).lt.0) THEN
              WRITE(6,*) 'After ',SUBR,': I,J,L,SSI=',I,J,L,SSI(:,I
     *             ,J),MSI(I,J),SNOWI(I,J),RSI(I,J)
              QCHECKI = .TRUE.
            END IF
          END DO
          IF (SNOWI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,SNOWI=',I,J,SNOWI(I,J)
            QCHECKI = .TRUE.
          END IF
          IF (MSI(I,J).gt.10000) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,MSI=',I,J,MSI(I,J),RSI(I,J)
c            QCHECKI = .TRUE.
          END IF
          END IF
        END DO
      END DO

#ifdef TRACERS_WATER
      do n=1,ntm
C**** check negative tracer mass
        if (t_qlimit(n)) then
        do j=J_0, J_1
          do i=1,imaxj(j)
            if ((focean(i,j)+flake(i,j))*rsi(i,j).gt.0) then
              do l=1,lmi
                if (trsi(n,l,i,j).lt.0.) then
                  print*,"Neg Tracer in sea ice after ",subr,i,j,l,
     *                 trname(n),trsi(n,l,i,j),rsi(i,j),msi(i,j),ssi(l,i
     *                 ,j)
                  QCHECKI=.true.
                end if
              end do
            end if
          end do
        end do
        end if
C**** Check conservation of water tracers in sea ice
        if (trname(n).eq.'Water') then
          errmax = 0. ; imax=1 ; jmax=1
          do j=J_0, J_1
          do i=1,imaxj(j)
            if ((focean(i,j)+flake(i,j))*rsi(i,j).gt.0) then
              relerr=max(
     *             abs(trsi(n,1,i,j)-(snowi(i,j)+ace1i)*xsi(1)+ssi(1,i,j
     *             ))/trsi(n,1,i,j),abs(trsi(n,2,i,j)-(snowi(i,j)+ace1i)
     *             *xsi(2)+ssi(2,i,j))/trsi(n,2,i,j),abs(trsi(n,3,i,j)
     *             -msi(i,j)*xsi(3)+ssi(3,i,j))/trsi(n,3,i,j),abs(trsi(n
     *             ,4,i,j)-msi(i,j)*xsi(4)+ssi(4,i,j))/trsi(n,4,i,j))
              if (relerr.gt.errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            end if
          end do
          end do
          print*,"Relative error in sea ice mass after ",trim(subr),":"
     *         ,imax,jmax,errmax,trsi(n,:,imax,jmax),(snowi(imax,jmax)
     *         +ace1i)*xsi(1)-ssi(1,imax,jmax),(snowi(imax,jmax)+ace1i)
     *         *xsi(2)-ssi(2,imax,jmax),msi(imax,jmax)*xsi(3:4)-ssi(3:4
     *         ,imax,jmax),rsi(imax,jmax),msi(imax,jmax)
        end if
      end do
#endif

      IF (QCHECKI)
     &     call stop_model("CHECKI: Ice variables out of bounds",255)

      END SUBROUTINE CHECKI

