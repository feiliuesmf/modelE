#include "rundeck_opts.h"

      MODULE SEAICE
!@sum  SEAICE contains all the sea ice related subroutines
!@auth Original Development Team
!@ver  1.0
!@cont PREC_SI,SEA_ICE,ADDICE,SIMELT
      USE CONSTANT, only : lhm,rhoi,byrhoi,rhow,shi,shw,byshi,bylhm,sday
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
     *     XSI= (/0.5d0, 0.5d0, 0.5d0, 0.5d0/), BYXSI=1./XSI
!@param Z1I thickness of first layer ice (m)
      REAL*8, PARAMETER :: Z1I = .1d0
!@param ACE1I ice mass first layer (kg/m^2)
      REAL*8, PARAMETER :: ACE1I = Z1I*RHOI
!@param HC1I heat capacity of first layer ice (J/m^2)
      REAL*8, PARAMETER :: HC1I = ACE1I*SHI
!@param Z2OIM thickness of 2nd layer ice (m)
      REAL*8, PARAMETER :: Z2OIM = .1d0    ! .4d0
!@param AC2OIM ice mass 2nd layer (kg/m^2)
      REAL*8, PARAMETER :: AC2OIM = Z2OIM*RHOI
!@param ALAMI,ALAMS lambda coefficient for ice/snow J/(m*degC*sec)
      REAL*8, PARAMETER :: ALAMI=2.1762d0, ALAMS=0.35d0
!@param RHOS density of snow (kg/m^3)
      REAL*8, PARAMETER :: RHOS = 300.0
!@param KIEXT extinction coefficient for light in sea ice (1/m)
      REAL*8, PARAMETER :: KIEXT = 1.5d0
!@param KSEXT extinction coefficient for light in snow (1/m)
      REAL*8, PARAMETER :: KSEXT = 15d0
!@var FLEADOC lead fraction for ocean ice (%)
      REAL*8, PARAMETER :: FLEADOC = 0.06d0
!@var FLEADLK lead fraction for lakes (%)
      REAL*8, PARAMETER :: FLEADLK = 0.
!@param BYRLI,BYRLS reciprical of density*lambda
      REAL*8, PARAMETER :: BYRLI = 1./(RHOI*ALAMI),
     *     BYRLS = 1./(RHOS*ALAMS)
!@param MU coefficient of seawater freezing point w.r.t. salinity
      REAL*8, PARAMETER :: MU = 0.054d0
!@param SSI0 default value for sea ice salinity (kg/kg) (=3.2ppt)
      REAL*8, PARAMETER :: SSI0 = 0.0032d0
!@param FSSS fraction of ocean salinity found in new-formed ice
      REAL*8, PARAMETER :: FSSS = 13d0/35d0
!@var qsfix flag is true if salinity of sea ice is constant
      LOGICAL :: QSFIX = .false.
!@var alpha implicitness for heat diffusion in sea ice (1=fully implicit)
      REAL*8, PARAMETER :: ALPHA = 1.0

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SALT,
#ifdef TRACERS_WATER
     *     TRSIL,TRPRCP,TRUN0,
#endif
     *     WETSNOW)
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
!@var TSIL temperature of ice layers (J/m^2) (DIAGNOSTIC ONLY)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var WETSNOW true if snow is wet (i.e. has been rained on)
      LOGICAL, INTENT(OUT) :: WETSNOW
!@var RUN0 runoff from ice (kg/m^2)
!@var SALT salt in runoff from ice (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0, SALT
      REAL*8 :: BYMSI2, FMSI1, FMSI2, FMSI3, FHSI1, FHSI2, FHSI3,
     *     CMPRS, SNWF, RAIN, FREZ1, MSI1, FSSI2,
     *     FSSI3, MELT1, SMELT12, DSNOW, SS12
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
      REAL*8, DIMENSION(NTM) :: FTRSI1,FTRSI2,FTRSI3,FTRSI4,TRRMF
     *     ,TRMELT1
#endif

C**** reciprocal of ice thickness for efficiency
      MSI1=SNOW+ACE1I
c      BYM1=1./(XSI(1)*MSI1)
c      BYM2=1./(XSI(2)*MSI1)
c      BYM3=1./(XSI(3)*MSI2)
c      BYM4=1./(XSI(4)*MSI2)

      HSIL(1) = HSIL(1)+ENRGP  ! add total energy of precipitation

C**** Snowfall is calculated from precip energy (0 deg or colder)
      SNWF = MAX(0d0,MIN(PRCP,-ENRGP*BYLHM))
C**** Rain is remaining precip (0 deg or warmer)
      RAIN = PRCP-SNWF
      WETSNOW = RAIN.GT.0.
C**** Calculate whether rain causes freezing or melting in first layer
      IF (HSIL(1).le.-LHM*(XSI(1)*MSI1+SNWF)) THEN
        FREZ1 = MIN(RAIN,-HSIL(1)*BYLHM-XSI(1)*MSI1-SNWF)
        MELT1 = 0.
      ELSE
        FREZ1 = 0.
        MELT1 = MAX(0d0,HSIL(1)*BYLHM+XSI(1)*MSI1+SNWF)
      END IF

C**** Calculate remaining snow and necessary mass flux using
C**** SNOW =SNOW  + SNWF  - CMPRS - MELT1
C**** ACE1I=ACE1I + FREZ1 + CMPRS - FMSI2
C**** SALTI=SALTI - SMELT12        - FSSI2

C**** Calculate changes to snow
      IF (SNOW+SNWF.GT.MELT1) THEN ! some snow remains
        IF (SNWF.ge.PRCP .and. SNOW+SNWF .GT. SNOMAX) THEN
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
          CMPRS = SNOW+SNWF-0.9d0*SNOMAX
        ELSE
C**** RAIN and MELT COMRESSES SNOW INTO ICE
          CMPRS = MIN(dSNdRN*(RAIN+MELT1), SNOW+SNWF-MELT1) ! cmpression
        END IF
        DSNOW = SNWF - (MELT1+CMPRS)
        SMELT12= 0.
      ELSE
C**** RAIN MELTS ALL SNOW AND SOME ICE
        CMPRS = 0.
        DSNOW = -SNOW
        SMELT12=(MELT1-SNOW-SNWF)*(SSIL(1)+SSIL(2))/ACE1I
      END IF
C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = SNWF+FREZ1-MELT1-DSNOW ! either up or down
      FMSI1 = XSI(1)*FMSI2+XSI(2)*(SNWF+FREZ1-MELT1)
#ifdef TRACERS_WATER
      TRSIL(:,1) = TRSIL(:,1)+TRPRCP(:)*(SNWF+FREZ1)/PRCP
      TRMELT1(:) = MELT1*TRSIL(:,1)/(XSI(1)*MSI1+SNWF+FREZ1)
      TRRMF(:) = TRPRCP(:)-TRPRCP(:)*(SNWF+FREZ1)/PRCP
#endif

C***** Calculate consequential heat/salt flux between layers
      IF (FMSI1.gt.0) THEN      ! downward flux to thermal layer 2
        FHSI1 = FMSI1*HSIL(1)/(XSI(1)*MSI1+SNWF+FREZ1-MELT1)
#ifdef TRACERS_WATER
        FTRSI1(:) = FMSI1*(TRSIL(:,1)-TRMELT1(:))/(XSI(1)*MSI1+SNWF
     *       +FREZ1-MELT1)
#endif
      ELSE                      ! upward flux
        FHSI1 = FMSI1*HSIL(2)/(XSI(2)*MSI1)
#ifdef TRACERS_WATER
        FTRSI1(:) = FMSI1*TRSIL(:,2)/(XSI(2)*MSI1)
#endif
      END IF
      IF (FMSI2.gt.0) THEN      ! downward flux to thermal layer 3
        FHSI2 = FMSI2*HSIL(2)/(XSI(2)*MSI1)
        FSSI2 = FMSI2*(SSIL(1)+SSIL(2))/ACE1I
#ifdef TRACERS_WATER
        FTRSI2(:) = FMSI2*TRSIL(:,2)/(XSI(2)*MSI1)
#endif
      ELSE                      ! upward flux
        FHSI2  = FMSI2*HSIL(3)/(XSI(3)*MSI2)
        FSSI2  = FMSI2*SSIL(3)/(XSI(3)*MSI2)
#ifdef TRACERS_WATER
        FTRSI2(:) = FMSI2*TRSIL(:,3)/(XSI(3)*MSI2)
#endif
      END IF

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
      SNOW = MAX(0d0,SNOW+DSNOW)
      MSI1 = ACE1I + SNOW
      MSI2 = MSI2 + FMSI2

      HSIL(1)=HSIL(1)- FHSI1
      HSIL(2)=HSIL(2)+(FHSI1-FHSI2)
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3

C**** salinity is evenly shared over ice portion
      SS12=(SSIL(1)+SSIL(2))-FSSI2-SMELT12
      IF (SNOW*XSI(2).gt.XSI(1)*ACE1I) THEN ! first layer is all snow
        SSIL(1)=0.
      ELSE
        SSIL(1)=SS12*(XSI(1)*ACE1I-SNOW*XSI(2))/ACE1I
      END IF
      SSIL(2)=SS12-SSIL(1)
      SSIL(3)=SSIL(3)+(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)+ FSSI3

#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,1)=TRSIL(:,1)-TRMELT1(:)-FTRSI1(:)
      TRSIL(:,2)=TRSIL(:,2)+(FTRSI1(:)-FTRSI2(:))
      TRSIL(:,3)=TRSIL(:,3)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)+ FTRSI3(:)
#endif

C**** Diagnostics for output
      RUN0 = MELT1+RAIN-FREZ1   ! runoff mass to ocean
      SALT = SMELT12             ! salt in runoff
c     HFLUX= 0                  ! energy of runoff (currently at 0 deg)      

#ifdef TRACERS_WATER
      TRUN0(:)=TRMELT1(:)+TRRMF(:)  ! tracer flux to ocean
#endif
      TSIL(1:2) = (HSIL(1:2)/(XSI(1:2)*MSI1)+LHM)*BYSHI ! ice temp L=1,2
      TSIL(3:4) = (HSIL(3:4)/(XSI(3:4)*MSI2)+LHM)*BYSHI ! ice temp L=3,4

      RETURN
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,EVAP
     *     ,SROX,
#ifdef TRACERS_WATER
     *     TRSIL,TREVAP,FTROC,
#endif
     *     FMOC,FHOC,FSOC,MELT12)
!@sum  SEA_ICE applies surface fluxes to ice covered areas
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER :: dSNdML =0.
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
      REAL*8, INTENT(INOUT) :: FMOC, FHOC, FSOC
!@var MELT12 amount of surface melting (kg/m^2) (Used for albedo calc)
      REAL*8, INTENT(OUT) :: MELT12
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI), INTENT(INOUT) :: TRSIL
!@var TREVAP tracer amount in evap (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TREVAP
!@var FTROC tracer flux to ocean (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: FTROC
      REAL*8, DIMENSION(NTM) :: FTRSI1,FTRSI2,FTRSI3,
     *     TRMELT1,TRMLET2,TRMELT3,TRMELT4,TRDEW,TRMELT2
#endif

      REAL*8, INTENT(IN) :: ROICE
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
      REAL*8, DIMENSION(LMI) :: TSIL
      REAL*8 MSI1, MELT1, MELT2, MELT3, MELT4, SNMELT, DSNOW
      REAL*8 DEW, CMPRS, DEW1, DEW2, EVAP1
      REAL*8 SMELT12,SMELT3,SMELT4,SALTI, FSSI2, FSSI3, SS12
      REAL*8 FMSI1, FMSI2, FMSI3, FHSI1, FHSI2, FHSI3
      REAL*8 HC1, HC2, HC3, HC4, HICE, HSNOW
      REAL*8 dF1dTI, dF2dTI, dF3dTI, dF4dTI, F1, F2, F3, FO

      IF (ROICE .EQ. 0.) RETURN
      FMSI2=0. ; MELT1=0. ; MELT2=0. ; MELT3=0. ; MELT4=0.
      FHSI1=0. ; FHSI2=0. ; FHSI3=0.
      SMELT12=0; SMELT3=0.; SMELT4=0.
      FSSI2=0. ; FSSI3=0.
#ifdef TRACERS_WATER
      FTRSI1=0. ; FTRSI2=0. ; FTRSI3=0. 
      TRMELT1=0. ; TRMLET2=0. ; TRMELT3=0. ; TRMELT4=0.
#endif
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C**** Calculate solar fractions
      IF (SROX(1).gt.0) THEN
        HICE =ACE1I*BYRHOI
        HSNOW=SNOW/RHOS
        FSRI(2) = EXP(-KSEXT*HSNOW-KIEXT*HICE)
        FSRI(3) = FSRI(2)*EXP(-KIEXT*MSI2*XSI(3)*BYRHOI)
        FSRI(4) = FSRI(2)*EXP(-KIEXT*MSI2*BYRHOI)
      ELSE
        FSRI = 0.
      END IF
C****
C**** OCEAN ICE, CALCULATE TSIL FROM ENTHALPY
C****
      TSIL(1:2) = (HSIL(1:2)/(XSI(1:2)*MSI1) +LHM)*BYSHI ! temp. L=1/2
      TSIL(3:4) = (HSIL(3:4)/(XSI(3:4)*MSI2) +LHM)*BYSHI ! temp. L=3/4

      HC1 = SHI*XSI(1)*MSI1 ! heat capacity of ice L=1 (J/(degC*m^2))
      HC2 = SHI*XSI(2)*MSI1 ! heat capacity of ice L=2 (J/(degC*m^2))
      HC3 = SHI*XSI(3)*MSI2 ! heat capacity of ice L=3 (J/(degC*m^2))
      HC4 = SHI*XSI(4)*MSI2 ! heat capacity of ice L=4 (J/(degC*m^2))
C**** CALCULATE AND APPLY DIFFUSIVE AND SURFACE ENERGY FLUXES
C**** First layer flux calculated already as part of surface calculation
c      dF1dTI = 2.*DTSRCE/(ACE1I*BYRLI+SNOW*BYRLS)
      dF2dTI = ALAMI*RHOI*DTSRCE/(0.5*XSI(2)*MSI1+0.5*XSI(3)*MSI2)
C****          temperature derivative from F2 diffusive flux
      dF3dTI = ALAMI*RHOI*DTSRCE*2./MSI2
C****          temperature derivative from F3 diffusive flux
      dF4dTI = ALAMI*RHOI*DTSRCE*2.*BYXSI(4)/MSI2
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
C**** Add solar flux through bottom to basal heat flux already calculated
      FHOC=SROX(1)*FSRI(4)+FHOC
      SROX(2) = SROX(1)*FSRI(4)

      HSIL(1) = HSIL(1)+(F0DT-F1DT)
      HSIL(2) = HSIL(2)+(F1DT-F2)
      HSIL(3) = HSIL(3)+(F2-F3)
      HSIL(4) = HSIL(4)+(F3-FHOC)

C**** add basal salt flux 
      SSIL(4) = SSIL(4)-FSOC
#ifdef TRACERS_WATER
      TRSIL(:,4) = TRSIL(:,4)-FTROC(:)
      TRDEW(:) = -TREVAP(:)
#endif

      DEW = -EVAP               ! dew/evap to the surface
      EVAP1= MAX(0d0,EVAP)      ! positive evaporation
      SALTI=SSIL(1)+SSIL(2)     ! total upper layer salt (only in ice)

C**** Add DEW directly to ice (which can be in either layer)
C**** DEW1 adds to layer 1, DEW2 adds to layer 2
C**** Calculate melting in first two thermal layers
C**** SNMELT is the melting that is applied to the snow first
      IF (SNOW*XSI(2).gt.XSI(1)*ACE1I) THEN ! first layer is all snow
        DEW1 = -EVAP1           ! <0 i.e. evaporation from snow
        DEW2 =  MAX(0d0,DEW)    ! >0 i.e. dew to second layer ice
        MELT1 = MAX(0d0,HSIL(1)*BYLHM+XSI(1)*MSI1+DEW1)
        MELT2 = MAX(0d0,HSIL(2)*BYLHM+XSI(2)*MSI1+DEW2)
        SNMELT=MELT1+MELT2
      ELSE  ! first layer is snow and some ice
        DEW1 = DEW   ! all fluxes to first layer
        DEW2 = 0.
        MELT1 = MAX(0d0,HSIL(1)*BYLHM+XSI(1)*MSI1+DEW1)
        MELT2 = MAX(0d0,HSIL(2)*BYLHM+XSI(2)*MSI1)
        SNMELT=MELT1
      END IF
#ifdef TRACERS_WATER
      if (DEW1.ne.0.) THEN
        TRSIL(:,1) = TRSIL(:,1)+TRDEW(:)
      ELSE
        TRSIL(:,2) = TRSIL(:,2)+TRDEW(:)
      END IF
      TRMELT1(:) = MELT1*TRSIL(:,1)/(XSI(1)*MSI1+DEW1)
      TRMELT2(:) = MELT2*TRSIL(:,2)/(XSI(2)*MSI1+DEW2)
#endif

C**** CMPRS is amount of snow turned to ice during melting
C**** Calculate remaining snow and necessary mass flux using
C**** SNOW =SNOW -       SNMELT        + DEW1 - CMPRS
C**** ACE1I=ACE1I-(MELT1+MELT2-SNMELT) + DEW2 + CMPRS - FMSI2
C**** SALTI=SALTI- SMELT12                            - FSSI2
C****
      IF (SNOW.GT.SNMELT+EVAP1) THEN ! some snow remains
        CMPRS = MIN(dSNdML*SNMELT, SNOW-EVAP1-SNMELT) ! > 0.
        DSNOW = - (SNMELT+EVAP1+CMPRS)
        SMELT12=(MELT1+MELT2-SNMELT)*SALTI/ACE1I
      ELSE ! all snow and some ice melts
        CMPRS = 0.
        DSNOW = -SNOW
        SMELT12=(MELT1+MELT2-SNOW-EVAP1)*SALTI/ACE1I
      END IF
C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = -DSNOW+DEW-MELT1-MELT2 ! either up or down
      FMSI1 = -XSI(1)*DSNOW+DEW1-MELT1-CMPRS

C**** Check for melting in levels 3 and 4
      MELT3 = MAX(0d0,HSIL(3)*BYLHM+XSI(3)*MSI2)
      MELT4 = MAX(0d0,HSIL(4)*BYLHM+XSI(4)*MSI2-FMOC)
      SMELT3 = MELT3*SSIL(3)/(XSI(3)*MSI2)
      SMELT4 = MELT4*SSIL(4)/(XSI(4)*MSI2-FMOC)
#ifdef TRACERS_WATER
      TRMELT3(:) = MELT3*TRSIL(:,3)/(XSI(3)*MSI2)
      TRMELT4(:) = MELT4*TRSIL(:,4)/(XSI(4)*MSI2-FMOC)
#endif

C***** Calculate consequential heat flux between layers
      IF (FMSI1.gt.0) THEN ! downward flux to thermal layer 2
        FHSI1 = FMSI1*HSIL(1)/(XSI(1)*MSI1+DEW1-MELT1)
#ifdef TRACERS_WATER
        FTRSI1(:)=FMSI1*(TRSIL(:,1)-TRMELT1(:))/(XSI(1)*MSI1+DEW1-MELT1)
#endif
      ELSE                 ! upward flux
        FHSI1 = FMSI1*HSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
#ifdef TRACERS_WATER
        FTRSI1(:)=FMSI1*(TRSIL(:,2)-TRMELT2(:))/(XSI(2)*MSI1+DEW2-MELT2)
#endif
      END IF
      IF (FMSI2.gt.0) THEN ! downward flux to thermal layer 3
        FHSI2 = FMSI2*HSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
        FSSI2 = FMSI2*SALTI/ACE1I
#ifdef TRACERS_WATER
        FTRSI2(:)=FMSI2*(TRSIL(:,2)-TRMELT2(:))/(XSI(2)*MSI1+DEW2-MELT2)
#endif
      ELSE                 ! upward flux
        FHSI2  = FMSI2*HSIL(3)/(XSI(3)*MSI2-MELT3)
        FSSI2  = FMSI2*(SSIL(3)-SMELT3)/(XSI(3)*MSI2-MELT3)
#ifdef TRACERS_WATER
        FTRSI2(:)=FMSI2*(TRSIL(:,3)-TRMELT3(:))/(XSI(3)*MSI2-MELT3)
#endif
      END IF

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
        FTRSI3(:)=FMSI3*(TRSIL(:,4)-TRMELT4(:))/(XSI(4)*MSI2-MELT4)
#endif
      END IF

C**** Apply fluxes
      SNOW = MAX(0d0,SNOW+DSNOW)
      MSI1 = ACE1I+ SNOW
      MSI2=MSI2-(MELT3+MELT4)+FMSI2-FMOC

      HSIL(1)=HSIL(1)- FHSI1
      HSIL(2)=HSIL(2)+(FHSI1-FHSI2)
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3

C**** salinity is evenly shared over ice portion
      SS12=(SSIL(1)+SSIL(2))-SMELT12-FSSI2
      IF (SNOW*XSI(2).gt.XSI(1)*ACE1I) THEN ! first layer is all snow
        SSIL(1)=0.
      ELSE
        SSIL(1)=SS12*(XSI(1)*ACE1I-SNOW*XSI(2))/ACE1I
      END IF
      SSIL(2)=SS12-SSIL(1)
      SSIL(3)=SSIL(3)-SMELT3+(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)-SMELT4+ FSSI3
#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,1)=TRSIL(:,1)- TRMELT1(:)- FTRSI1(:)
      TRSIL(:,2)=TRSIL(:,2)- TRMELT2(:)+(FTRSI1(:)-FTRSI2(:))
      TRSIL(:,3)=TRSIL(:,3)- TRMELT3(:)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)- TRMELT4(:)+ FTRSI3(:)
#endif

C**** Calculate net output fluxes 
      FMOC = FMOC+MELT1+MELT2+MELT3+MELT4 ! mass flux to ocean
      FSOC = FSOC+SMELT12+SMELT3+SMELT4   ! salt flux to ocean
c HMELT currently assumed to be zero since melting is at 0 deg 
c     FHOC = FHOC+HMELT 
#ifdef TRACERS_WATER
      FTROC(:)=FTROC(:)+TRMELT1(:)+TRMELT2(:)+TRMELT3(:)+TRMELT4(:) 
                                ! tracer flux to ocean
#endif
c Save MELT12 seperately for albedo calculations
      MELT12=MELT1+MELT2
C****
      RETURN
      END SUBROUTINE SEA_ICE

      SUBROUTINE ADDICE (SNOW,ROICE,HSIL,SSIL,MSI2,TSIL,ENRGFO
     *     ,ACEFO,ACE2F,ENRGFI,SALTO,SALTI,
#ifdef TRACERS_WATER
     *     TRSIL,TRO,TRI,
#endif
     *     FLEAD,QFIXR)
!@sum  ADDICE adds ice formed in the ocean to ice variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER, DIMENSION(LMI) :: YSI =
     *     (/XSI(1)* ACE1I/(ACE1I+AC2OIM),XSI(2)*ACE1I/(ACE1I+AC2OIM),
     *       XSI(3)*AC2OIM/(ACE1I+AC2OIM),XSI(4)*AC2OIM/(ACE1I+AC2OIM)/)
      REAL*8, PARAMETER :: Z2OIX = 4.9,
     *                     BYZICX=1./(Z1I+Z2OIX)
!@var QFIXR  true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL, INTENT(IN) :: QFIXR
!@var FLEAD minimum lead fraction for ice (%)
      REAL*8, INTENT(IN) :: FLEAD
      REAL*8, INTENT(IN) ::  ENRGFI, ENRGFO, ACEFO, ACE2F, SALTO, SALTI
!@var ROICE,SNOW,MSI1,MSI2
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL,SSIL
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: trsil
      REAL*8, INTENT(IN), DIMENSION(NTM) :: tro,tri
      REAL*8, DIMENSION(NTM) :: FTRSI3,FTRSI4
      INTEGER N
#endif

      REAL*8 FMSI3, FMSI4, FHSI3, FHSI4, FSSI3, FSSI4
      REAL*8 ROICEN, OPNOCN, DRSI, MSI1

      MSI1=SNOW+ACE1I
      IF (.not. QFIXR) THEN
      IF (ROICE.LE.0. .and. ACEFO.gt.0) THEN
        ROICE=ACEFO/(ACE1I+AC2OIM)
        SNOW=0.
C****   TSIL=(ENRGFO/ACEFO + LHM)*BYSHI ! but hsi is primary var.
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

      IF (ACE2F.le.0) GO TO 250 ! go to no freezing case
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
CC    FMSI3 = -XSI(3)*ACE2F ! < 0.
CC    FMSI4 = -ACE2F
      FHSI3 = -HSIL(4)*ACE2F*(XSI(3)/XSI(4))/MSI2
      FSSI3 = -SSIL(4)*ACE2F*(XSI(3)/XSI(4))/MSI2
#ifdef TRACERS_WATER
      FTRSI3(:) = -TRSIL(:,4)*ACE2F*(XSI(3)/XSI(4))/MSI2
#endif
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      IF (ACEFO .GT. 0.) GO TO 240
C**** NEW ICE IS FORMED BELOW OLD SEA ICE
      HSIL(3) = HSIL(3)-FHSI3
      HSIL(4) = HSIL(4)+(FHSI3+ENRGFI)
      SSIL(3) = SSIL(3)-FSSI3
      SSIL(4) = SSIL(4)+(FSSI3+SALTI)
#ifdef TRACERS_WATER
      TRSIL(:,3) = TRSIL(:,3) - FTRSI3(:)
      TRSIL(:,4) = TRSIL(:,4) +(FTRSI3(:)+TRI(:))
#endif
      MSI2 = MSI2+ACE2F ! new ice mass of physical layer 2
      GO TO 270
  240 CONTINUE
C**** NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      MSI1 = SNOW + ACE1I ! mass of layer 1
C**** MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACE2F))/(ROICE+DRSI) ! layer 2
      HSIL(1:2)=((1.-ROICE)*ENRGFO*YSI(1:2)+ROICE*HSIL(1:2))/
     *     (ROICE+DRSI)
      HSIL(3) = ((1.-ROICE)*ENRGFO*YSI(3)+ROICE*(HSIL(3)-FHSI3))/
     A     (ROICE+DRSI)
      HSIL(4) = ((1.-ROICE)*ENRGFO*YSI(4)+ROICE*(HSIL(4)+FHSI3+ENRGFI))/
     A     (ROICE+DRSI)
      SSIL(1:2)=((1.-ROICE)*SALTO*YSI(1:2)+ROICE*SSIL(1:2))/(ROICE+DRSI)
      SSIL(3) = ((1.-ROICE)*SALTO*YSI(3)+ROICE*(SSIL(3)-FSSI3))/
     A     (ROICE+DRSI)
      SSIL(4) = ((1.-ROICE)*SALTO*YSI(4)+ROICE*(SSIL(4)+FSSI3+SALTI))/
     A     (ROICE+DRSI)
#ifdef TRACERS_WATER
      DO N=1,NTM
        TRSIL(N,1:2)=((1.-ROICE)*TRO(N)*YSI(1:2)+ROICE*TRSIL(N,1:2))
     *       /(ROICE+DRSI)
        TRSIL(N,3) = ((1.-ROICE)*TRO(N)*YSI(3)+ROICE*(TRSIL(N,3)
     *       -FTRSI3(N)))/(ROICE+DRSI)
        TRSIL(N,4) = ((1.-ROICE)*TRO(N)*YSI(4)+ROICE*(TRSIL(N,4)
     *       +FTRSI3(N)+TRI(N)))/(ROICE+DRSI)
      END DO
#endif
      ROICE = ROICE+DRSI ! new ice concentration
      GO TO 270
  250 CONTINUE
      IF (ACEFO .GT. 0.) THEN ! new ice on the open ocean
C**** NEW ICE IS FORMED ON THE OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      MSI1 = SNOW + ACE1I ! mass of layer 1
C**** MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*MSI2)/(ROICE+DRSI) ! layer 2
      HSIL(1:4)=((1.-ROICE)*ENRGFO*YSI(1:4)+ROICE*HSIL(1:4))/
     *     (ROICE+DRSI)
      SSIL(1:4)=((1.-ROICE)*SALTO*YSI(1:4)+ROICE*SSIL(1:4))/(ROICE+DRSI)
#ifdef TRACERS_WATER
      DO N=1,NTM
        TRSIL(N,1:4)=((1.-ROICE)*TRO(N)*YSI(1:4)+ROICE*TRSIL(N,1:4))
     *       /(ROICE+DRSI)
      END DO
#endif
      ROICE = ROICE+DRSI ! new ice concentration
      END IF
  270 CONTINUE
C**** COMPRESS THE ICE HORIZONTALLY IF TOO THIN OR LEAD FRAC. TOO SMALL
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ACE1I+MSI2)) ! -BYZICX) sometime -ve!
      IF (MSI2.LT.AC2OIM .or. ROICE.GT.1.-OPNOCN) THEN
      ROICEN = MIN(ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM),1.-OPNOCN)
      DRSI = ROICEN-ROICE ! < 0. compressed ice concentration
C     FMSI3 = XSI(3)*FMSI4 ! < 0. upward ice mass flux into layer 3
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
C     SNOW = SNOW   ! snow thickness is conserved
      ROICE = ROICEN
C****
      END IF
C****
      END IF
C**** Clean up ice fraction (if rsi>(1-OPNOCN)-1d-4) => rsi=(1-OPNOCN))
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ACE1I+MSI2))    ! -BYZICX)
      IF (ROICE.gt.(1.-OPNOCN-1d-4)) THEN
        ROICEN = 1.-OPNOCN
        DRSI = ROICEN-ROICE ! > 0
        FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! upward ice mass into layer 4
        FHSI3 = HSIL(4)*FMSI4*(XSI(3)/XSI(4))/MSI2 ! upward heat flux
        FHSI4 = (HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4))*(DRSI/ROICEN)
        HSIL(3) = HSIL(3)-FHSI3
        HSIL(4) = HSIL(4)+(FHSI3-FHSI4)
        FSSI3 = SSIL(4)*FMSI4*(XSI(3)/XSI(4))/MSI2 ! upward heat flux
        FSSI4 = (SSIL(1)+SSIL(2)+SSIL(3)+SSIL(4))*(DRSI/ROICEN)
        SSIL(3) = SSIL(3)-FSSI3
        SSIL(4) = SSIL(4)+(FSSI3-FSSI4)
#ifdef TRACERS_WATER
        FTRSI3(:) = TRSIL(:,4)*FMSI4*(XSI(3)/XSI(4))/MSI2
        FTRSI4(:) = (TRSIL(:,1)+TRSIL(:,2)+TRSIL(:,3)+TRSIL(:,4))*(DRSI
     *       /ROICEN)
        TRSIL(:,3) = TRSIL(:,3) - FTRSI3(:)
        TRSIL(:,4) = TRSIL(:,4) +(FTRSI3(:)-FTRSI4(:))
#endif
        MSI2 = MSI2-FMSI4       ! new ice mass of second physical layer
C       SNOW = SNOW   ! snow thickness is conserved
        ROICE = ROICEN
      END IF
      END IF
C**** Calculate temperatures for diagnostics and radiation
      TSIL(1:2)=(HSIL(1:2)/(XSI(1:2)*MSI1)+LHM)*BYSHI ! temp. layer 1/2
      TSIL(3:4)=(HSIL(3:4)/(XSI(3:4)*MSI2)+LHM)*BYSHI ! temp. layer 3/4

      RETURN
      END SUBROUTINE ADDICE

      SUBROUTINE SIMELT(ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,TFO,TSIL,
#ifdef TRACERS_WATER
     *     TRSIL,TRUN0,
#endif
     *     ENRGUSED,RUN0,SALT)
!@sum  SIMELT melts sea ice if it is too small
!@auth Original Development Team
!@ver  1.1
      IMPLICIT NONE
!@var POCEAN ocean fraction (zero if lake)
      REAL*8, INTENT(IN) :: POCEAN
!@var TFO freezing temperature of water (C)
      REAL*8, INTENT(IN) :: TFO
!@var ROICE,SNOW,MSI2 ice variables (%,kg/m^2,kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var HSIL ice enthalpy  (J/m^2)
!@var SSIL ice salt (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL, SSIL
!@var TSIL ice temperature (C)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var ENRGUSED energy used to melt ice (J/m^2)
      REAL*8, INTENT(OUT) :: ENRGUSED
!@var RUN0 amount of sea ice melt (kg/m^2)
!@var SALT amount of salt in sea ice melt (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0, SALT
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRUN0
#endif

      IF (ROICE.lt.1d-4) THEN   ! remove too small ice
        ENRGUSED=-ROICE*(HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4)) !  [J/m^2]
        RUN0=ROICE*(SNOW + ACE1I + MSI2)  ! all ice is melted
        SALT=ROICE*(SSIL(1)+SSIL(2)+SSIL(3)+SSIL(4))
        ROICE=0.
C**** set defaults
        SNOW=0.
        MSI2=AC2OIM
        HSIL(1:2)=(SHI*TFO-LHM)*XSI(1:2)*ACE1I
        HSIL(3:4)=(SHI*TFO-LHM)*XSI(3:4)*AC2OIM
        IF (POCEAN.gt.0) THEN
          SSIL(1:2)=SSI0*XSI(1:2)*ACE1I
          SSIL(3:4)=SSI0*XSI(3:4)*AC2OIM
        ELSE
          SSIL(:) = 0.
        END IF
        TSIL=TFO
#ifdef TRACERS_WATER
        TRUN0(:)=TRSIL(:,1)+TRSIL(:,2)+TRSIL(:,3)+TRSIL(:,4)
        TRSIL(:,:)=0.
#endif
      END IF
C****
      RETURN
      END SUBROUTINE SIMELT

      SUBROUTINE SSIDEC(I0,J0,MSI1,MSI2,HSIL,SSIL,DT,
#ifdef TRACERS_WATER
     *             TRSIL,TRFLUX,
#endif
     *     MFLUX,HFLUX,SFLUX)
!@sum  SSIDEC decays salinity in sea ice 
!@auth Jiping Liu
!@ver  1.0
      IMPLICIT NONE

      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: SSIL,HSIL
      REAL*8, INTENT(INOUT) :: MSI2
      REAL*8, INTENT(IN) :: DT, MSI1
      REAL*8, INTENT(OUT) :: MFLUX,HFLUX,SFLUX
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRFLUX
      REAL*8, DIMENSION(NTM) :: FTRSI1,FTRSI2,FTRSI3
      REAL*8, DIMENSION(NTM,LMI) :: DTRSI
      INTEGER N
#endif

      REAL*8 DSSI(LMI),DMSI(LMI),DHSI(LMI),SS12,DS12
      REAL*8 FMSI1,FMSI2,FMSI3,FSSI1,FSSI2,FSSI3,FHSI1,FHSI2,FHSI3
      INTEGER L,I0,J0
!@var dtssi decay time scale for sea ice salinity (days)
!@var bydtssi decay constant for sea ice salinity (1/s)
      REAL*8, parameter :: dtssi=30.d0, bydtssi=1./(dtssi*sday)

      DSSI = 0. ; DMSI = 0. ;  DHSI = 0.          
#ifdef TRACERS_WATER
      DTRSI(:,:)= 0.
#endif
C**** salinity in sea ice decays if it is above threshold level ssi0

C**** check first layer (default ice and snow)
      IF(SSIL(1)+SSIL(2).GT.ssi0*ACE1I) THEN
        DS12 = (SSIL(1)+SSIL(2)-ACE1I*ssi0)*DT*BYDTSSI
        IF (ACE1I.gt.XSI(2)*MSI1) THEN 
          DSSI(1) = (ACE1I-XSI(2)*MSI1)*DS12/ACE1I
          DSSI(2) = DS12-DSSI(1)
        ELSE
          DSSI(1) = 0.
          DSSI(2) = DS12
        END IF
        DMSI(1:2) = DSSI(1:2)
        DHSI(1:2) = DMSI(1:2)*HSIL(1:2)/(XSI(1:2)*MSI1)
        SSIL(1:2) = SSIL(1:2)-DSSI(1:2)
        HSIL(1:2) = HSIL(1:2)-DHSI(1:2)
#ifdef TRACERS_WATER
        DO N=1,NTM
          DTRSI(N,1:2) = DMSI(1:2)*TRSIL(N,1:2)/(XSI(1:2)*MSI1)
          TRSIL(N,1:2) = TRSIL(N,1:2)-DTRSI(N,1:2)
        END DO
#endif
      END IF

C**** check remaining layers
      DO L=3,LMI
        IF (SSIL(L).GT.ssi0*XSI(L)*MSI2) THEN 
          DSSI(L) = (SSIL(L)-ssi0*XSI(L)*MSI2)*DT*BYDTSSI
          DMSI(L) = DSSI(L)
          DHSI(L) = DMSI(L)*HSIL(L)/(XSI(L)*MSI2)
          SSIL(L) = SSIL(L) - DSSI(L)
          HSIL(L) = HSIL(L) - DHSI(L)
#ifdef TRACERS_WATER
          DTRSI(:,L) = DMSI(L)*TRSIL(:,L)/(XSI(L)*MSI2)
          TRSIL(:,L) = TRSIL(:,L)-DTRSI(:,L)
#endif
        END IF
      END DO
      
C**** Calculate fluxes required to rebalance layers
C**** Mass/heat moves from layer 2 to 1 (salt as well, but this is
C**** dealt with later)
      FMSI1 = -DMSI(1)
      FHSI1 = FMSI1*HSIL(2)/(XSI(2)*MSI1-DMSI(2))
#ifdef TRACERS_WATER
      FTRSI1(:) = FMSI1*TRSIL(:,2)/(XSI(2)*MSI1-DMSI(2))
#endif
C**** Mass/heat/salt moves from layer 3 to 2 
      FMSI2 = -(DMSI(1)+DMSI(2))
      FHSI2 = FMSI2*HSIL(3)/(XSI(3)*MSI2-DMSI(3))
      FSSI2 = FMSI2*SSIL(3)/(XSI(3)*MSI2-DMSI(3))
#ifdef TRACERS_WATER
      FTRSI2(:) = FMSI2*TRSIL(:,3)/(XSI(3)*MSI2-DMSI(3))
#endif
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

C**** Apply the fluxes
      HSIL(1)=HSIL(1)- FHSI1
      HSIL(2)=HSIL(2)+(FHSI1-FHSI2)
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3
C**** salinity spread evenly over upper ice layer
      SS12=(SSIL(1)+SSIL(2))-FSSI2
      IF (ACE1I.gt.XSI(2)*MSI1) THEN 
        SSIL(1)=SS12*(ACE1I-XSI(2)*MSI1)/ACE1I
      ELSE
        SSIL(1)= 0.
      END IF
      SSIL(2)=SS12-SSIL(1)
      SSIL(3)=SSIL(3)+(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)+ FSSI3
#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,1)=TRSIL(:,1)- FTRSI1(:)
      TRSIL(:,2)=TRSIL(:,2)+(FTRSI1(:)-FTRSI2(:))
      TRSIL(:,3)=TRSIL(:,3)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)+ FTRSI3(:)
#endif
c     MSI1 = MSI1 - (DMSI(1)+DMSI(2)) - FMSI2 ! stays fixed
      MSI2 = MSI2 - (DMSI(3)+DMSI(4)) + FMSI2

C**** output fluxes and diagnostics 
      MFLUX = SUM(DMSI)         ! mass flux to ocean
      HFLUX = SUM(DHSI)         ! energy flux to ocean
      SFLUX = SUM(DSSI)         ! salt flux to ocean
#ifdef TRACERS_WATER
      DO N=1,NTM
        TRFLUX(N)=SUM(DTRSI(:,N)) ! tracer flux to ocean
      END DO
#endif
C****
      RETURN
      END SUBROUTINE SSIDEC

      subroutine iceocean(Ti,Si,Tm,Sm,dh,ustar,Coriol,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *     Tri,Trm,trflux,tralpha,
#endif
     *     mflux,sflux,hflux)
!@sum  iceocean calculates fluxes at base of sea ice 
!@auth Gavin Schmidt
!@ver  1.0
!@usage At the ice-ocean interface the following equations hold:
!@+                                              Tb = -mu Sb         (1)
!@+ -lam_i(Ti,Si)(Ti-Tb)/dh + rho_m shw g_T (Tb-Tm) = -m Lh(Tib,Sib)
!@+                                                  + m shw (Tib-Tb)(2)
!@+                           rho_m     g_S (Sb-Sm) =  m (Sib-Sb )   (3)
!@+  with  Sib=Si, Tib=Ti m>0,  or Sib = fsss Sm, Tib=Tb m<0        (4)
!@+ This routine iteratively solves for m, Tb, Sb and Sib, Tib using
!@+ Newton's method. The form of the latent heat and conductivity  
!@+ can optionally include salinity effects.
!@+ The actual fluxes (positive down) are then:
!@+     mflux = m   ; sflux = 1d-3 m Sib (kg/m^2/s)
!@+     hflux = lam_i(Ti,Si) (Ti-Tb)/dh -m Lh(Tib,Sib)+m shw Tib  (J/m^2/s)
!@+
      implicit none
!@var alamdS salinity coefficient for conductivity (from common?)
      real*8, parameter :: alamdS=0.117d0

!@var G_mole_T,G_mole_S molecular diffusion terms for heat/salt
C****  G_mole_T = 12.5 * 13.8d0**(2d0/3d0) - 6. = 65.9d0
C****  G_mole_S = 12.5 * 2432d0**(2d0/3d0) - 6. = 2255d0
      real*8, parameter :: G_mole_T = 65.9d0 , G_mole_S = 2255d0
!@var Si salinity in lowest ice layer (psu)
!@var Ti,Tm temperatures in ice and mixed layer (C)
!@var dh distance from center of bottom ice layer to base of ice (m)
      real*8, intent(in) :: Ti,Si,Tm,Sm,dh
!@var coriol Corilos parameter (used in turbulent flux calc) (1/s)
!@var ustar friction velocity at ice-ocean interface (m/s)
      real*8, intent(in) :: ustar,Coriol
!@var dtsrc source time step (s)
!@var mfluxmax maximum melt rate allowed (kg/m^2 s)
!@var mlsh mixed layer specific heat capactity (J/m^2 C) 
      real*8, intent(in) :: dtsrc,mlsh          !,mfluxmax
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
!@var Sb,Sb0 final and initial estimate for salinity at interface (psu)
      real*8 :: G_turb,g_T,g_S,Sb,Sb0
!@var dSbdSb differential of Sb with respect to initial Sb0
      real*8 :: dSbdSb
!@var m melt rate (negative implies freezing) (kg/m^2/s)
!@var Tb temperature at interface (C)
!@var Tib,Sib actual temperature/salinity in melting/freezing ice (psu)
      real*8 :: m,Tb,Sib,Tib
      real*8 :: lh,left2,df3dm,dmdTb,dmdSi,alamdh,f0,df,rsg
      integer, parameter :: niter=5 !@param niter number of iterations 
      integer :: i

C**** Calculate turbulent exchange velocities g_T, g_S

C****  G_turb = (1/k) ln (u* E n^2 / f h) + 1 / (2 E n) - (1/k)
C****         = 2.5 ln ( 5300*(u*)^2/coriol ) + 9.62 - 2.5 (assuming n=1)
C**** Note: n should theoretically depend on the buoyancy flux and
C****       therfore this should be included in the iteration below.
      G_turb = 2.5d0 * log ( 5300.*ustar*ustar/coriol ) + 7.12d0
C****  g  = u* / ( G_turb + G_mole )
      g_T = ustar / ( G_turb + G_mole_T )
      g_S = ustar / ( G_turb + G_mole_S )

C**** set conductivity term
C**** Diffusive flux is implicit in ice + ml temperature
      rsg=rhow*shw*g_T   !/(1.+alpha*dtsrc*rhow*shw*g_T/mlsh)
c no salinity effects
      alamdh = alami/(dh+alpha*dtsrc*alami*byshi/(2.*dh*rhoi))
c S thermo 
c      alamdh = (alami*Ti+alamdS*Si)/(dh*Ti+alpha*dtsrc*(alami*Ti+alamdS*Si)
c     *          *byshi/(2.*rhoi*dh))

C**** solve for boundary values (uses Newton-Rapheson with bounds)
C**** Estimate initial boundary salinity
      Sb0 = 0.75d0*Sm+0.25d0*Si   
      do i=1,niter
        Tb = -mu*Sb0            ! freezing point at salinity Sb0
C**** calculate left hand side of equation 2
        left2 = -alamdh*(Ti-Tb) + rsg*(Tb-Tm)
C**** depending on whether there is freezing or melting, the sea ice
C**** salinity in the melt/freeze fraction is set to be the original
C**** ice salinity or is a function of the boundary salinity. 
C**** Similarly for temperature Tib
        if (left2.gt.0) then    ! freezing   
          if (qsfix) then   ! keep salinity in ice constant
            Sib = ssi0
          else              ! it is a function of boundary value
            Sib = Sb0*fsss
          end if
c no salinity effects
          lh = lhm+ Tb*(shw-shi)
c S thermo
c         lh = lhm*(1.+mu*Sib/Tb) + (Tb+mu*Sib)*(shw-shi) 
          m = -left2/lh

c no salinity effects
          dmdTb = (left2*(shw-shi)-lh*(alamdh + rsg))/(lh*lh)
          dmdSi = 0.
c S thermo
c         dmdTb = (left2*(-lhm*mu*Sib/Tb**2+shw-shi)-lh*(alamdh + rsg
c    *       ))/(lh*lh) 
c         dmdSi = (left2*mu*(lhm/Tb+shw-shi))/(lh*lh)

          if (qsfix) then   ! keep salinity in ice constant
            Sb = (rhow*g_S*Sm+m*Sib)/(rhow*g_S+m)
            df3dm= rhow*g_S*(Sib-Sm)/(rhow*g_S+m)**2
            dSbdSb= -mu*dmdTb*df3dm
          else              ! it is a function of boundary value
            Sb = rhow*g_S*Sm/(rhow*g_S+m*(1.-fsss))
            df3dm= -rhow*g_S*Sm*(1.-fsss)/(rhow*g_S+m*(1.-fsss))**2
            dSbdSb= (fsss*dmdSi-mu*dmdTb)*df3dm
          end if
        else                    ! melting
          Sib = Si
c no salinity effects
          lh = lhm + Tb*shw - Ti*shi
c S thermo
c         lh = lhm*(1.+mu*Sib/Ti) + (Ti+mu*Sib)*(shw-shi) - shw*(Ti-Tb)
          m = -left2/lh

          Sb = (m*Sib+rhow*g_S*Sm)/(rhow*g_S+m)
          df3dm= (Sib*(rhow*g_S+m)-(m*Sib+rhow*g_S*Sm))/(rhow*g_S+m)**2
          dmdTb = -(alamdh + rsg)/lh + left2*shw/lh**2
          dSbdSb= -mu*dmdTb*df3dm
        end if

        f0=Sb-Sb0
        df=dSbdSb-1.+1d-20
        Sb = min(max(Sib,Sb0 - f0/df),40.)
        Sb0=Sb
      end do
#ifdef TRACERS_WATER
C**** Tracers use salinity turbulent diffusion term
      if (m.gt.0) then
        Trib(:)=Tri(:)
      else
        Trib(:)=tralpha(:)*rhow*g_S*Trm(:)/(rhow*g_S+m*(1.-tralpha(:)))
      end if
#endif
C**** define fluxes (positive down)
C**** Cap mass flux at at 90% of bottom layer
      if (m.gt.0.9d0*2.*dh*rhoi/dtsrc) then
        m=0.9d0*2.*dh*rhoi/dtsrc
      end if
      mflux = m                       ! (kg/m^2 s)
      sflux = 1d-3*m*Sib              ! (kg/m^2 s)
      hflux = alamdh*(Ti-Tb) - m*lh +m*shw*Tb ! (J/m^2 s)
#ifdef TRACERS_WATER
      trflux(:) = m * Trib(:)
#endif
C****
      return
      end subroutine iceocean

      subroutine icelake(Ti,Tm,dh,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *     Tri,Trm,trflux,tralpha,
#endif
     *     mflux,hflux)
!@sum  icelake calculates fluxes at base of lake ice (no salinity)
!@auth Gavin Schmidt
!@ver  1.0
!@usage At the ice-lake interface the following equations hold:
!@+   interface is at freezing point (Tb=0.)   (1)
!@+   -lam_i Ti/dh - rho_m shw g_T Tm = -m Lh(Tib) + m shw Tib   (2)
!@+     with  Tib=Ti m>0,  or Tib=0. m<0        (4)
      implicit none
      real*8, intent(out) :: mflux,hflux
      real*8, intent(in) :: Ti,Tm,dh,dtsrc,mlsh     !,mfluxmax
C**** Assume constant g_T = 5d-5, g_S = 0.04 * g_T m/s
!@var rsg = rhow * shw * g_T turbulent energy flux (J/m^2 s)
!@var rgS = rhow * g_S turbulent tracer flux (kg/m^2 s)
      real*8, parameter ::  rsg = rhow*shw*5d-5, rgS=rhow*2d-8
      real*8 left2, lh, m, alamdh
#ifdef TRACERS_WATER
!@var Tri,Trm tracer concentration in ice and mixed layer (kg/kg)
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

      END MODULE SEAICE

      MODULE SEAICE_COM
!@sum  SEAICE_COM contains the model arrays for seaice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE SEAICE, only : lmi
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif

      IMPLICIT NONE
!@var RSI fraction of open water area covered in ice
      REAL*8, DIMENSION(IM,JM) :: RSI
!@var SNOWI snow amount on sea ice (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: SNOWI
!@var MSI mass of ice second layer (layer 1=const) (kg/m^2)
C**** Note that MSI includes the mass of salt in sea ice
      REAL*8, DIMENSION(IM,JM) :: MSI
!@var HSI enthalpy of each ice layer (J/m^2)
      REAL*8, DIMENSION(LMI,IM,JM) :: HSI
!@var SSI sea ice salt content (kg/m^2)
      REAL*8, DIMENSION(LMI,IM,JM) :: SSI
!@var pond_melt amount of melt pond mass (kg/m^2)
C**** Note this is a virtual meltpond and is only used for
C**** albedo calculations
      REAL*8, DIMENSION(IM,JM) :: pond_melt
!@var flag_dsws true if snow on ice is wet (ie. rain or surface melt)
      LOGICAL, DIMENSION(IM,JM) :: flag_dsws

#ifdef TRACERS_WATER
!@var TRSI tracer amount in sea ice (kg/m^2)
      REAL*8, DIMENSION(NTM,LMI,IM,JM) :: TRSI
!@var TRSI0 default tracer conc. in sea ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSI0 = 1.
#endif

      END MODULE SEAICE_COM

      SUBROUTINE io_seaice(kunit,iaction,ioerr)
!@sum  io_seaice reads and writes seaice variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irsfic,irerun
      USE SEAICE_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SICE02"
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSICE01"

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a1,i3,a)')'R8 TRSI(',ntm,',',lmi,',im,jm)'
#endif

      write (MODULE_HEADER(lhead+1:80),'(a14,i1,a20,i1,a23)')
     * 'R8 F(im,jm),H(',lmi,',im,jm),snw,msi,ssi(',lmi,
     *     '),pond_melt,L flag_dsws'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
         WRITE (kunit,err=10) MODULE_HEADER,RSI,HSI,SNOWI,MSI,SSI
     *     ,POND_MELT,FLAG_DSWS
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TRSI
#endif
      CASE (IOREAD:)            ! input from restart file
         READ (kunit,err=10) HEADER,RSI,HSI,SNOWI,MSI,SSI
     *     ,POND_MELT,FLAG_DSWS
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (IRERUN,IOREAD)    ! only need tracers from reruns/restarts
          READ (kunit,err=10) TRHEADER,TRSI
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",TRHEADER
     *           ,TRMODULE_HEADER
            GO TO 10
          END IF
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
      USE CONSTANT, only : lhm,shi,rhow
      USE MODEL_COM
      USE SEAICE, only : lmi,xsi,ace1i
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
      USE FLUXES, only : UI2rho
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var QCHECKI true if errors found in seaice
      LOGICAL QCHECKI
      INTEGER I,J,L
      REAL*8 TICE

C**** Check for NaN/INF in ice data
      CALL CHECK3(RSI,IM,JM,1,SUBR,'rs')
      CALL CHECK3(MSI,IM,JM,1,SUBR,'ms')
      CALL CHECK3(HSI,LMI,IM,JM,SUBR,'hs')
      CALL CHECK3(SSI,LMI,IM,JM,SUBR,'ss')
      CALL CHECK3(SNOWI,IM,JM,1,SUBR,'sn')

      QCHECKI = .FALSE.
C**** Check for reasonable values for ice variables
      DO J=1,JM
        DO I=1,IM
          IF (RSI(I,J).lt.0 .or. RSI(I,j).gt.1 .or. MSI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,RSI,MSI=',I,J,RSI(I,J)
     *           ,MSI(I,J)
            QCHECKI = .TRUE.
          END IF
          DO L=1,LMI
            IF (L.le.2) TICE = (HSI(L,I,J)/(XSI(L)*(ACE1I+SNOWI(I,J)))
     *           +LHM)/SHI
            IF (L.gt.2) TICE = (HSI(L,I,J)/(XSI(L)*MSI(I,J))+LHM)/SHI
            IF (HSI(L,I,J).gt.0.or.TICE.gt.1d-4.or.TICE.lt.-80.) THEN
              WRITE(6,*) 'After ',SUBR,': I,J,L,TSI=',I,J,L,TICE,HSI(:,I
     *             ,J),MSI(I,J),SNOWI(I,J),RSI(I,J),sqrt(UI2RHO(I,J)
     *             /rhow)
c              QCHECKI = .TRUE.
              STOP
            END IF
            IF (SSI(L,I,J).lt.0) THEN
              WRITE(6,*) 'After ',SUBR,': I,J,L,SSI=',I,J,L,SSI(:,I
     *             ,J),MSI(I,J),SNOWI(I,J),RSI(I,J)
              QCHECKI = .TRUE.
              STOP
            END IF
          END DO
          IF (SNOWI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,SNOWI=',I,J,SNOWI(I,J)
            QCHECKI = .TRUE.
            STOP
          END IF
        END DO
      END DO
      IF (QCHECKI) STOP "CHECKI: Ice variables out of bounds"

      END SUBROUTINE CHECKI
