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
!@var alpha implicity for heat diffusion in sea ice (1=fully implicit)
      REAL*8, PARAMETER :: ALPHA = 1.0
!@dbparam oi_ustar0 default ice-ocean friction velocity (m/s) 
      REAL*8 :: oi_ustar0 = 1d-3  ! 5d-3 ! not used if ice dynamics is
!@dbparam silmfac factor controlling lateral melt of ocean ice
      REAL*8 :: silmfac = 1.4d-8 ! = pi*(3d-6)/0.66/1000
!@var silmpow exponent for temperature dependence of lateral melt
      REAL*8 :: silmpow = 1.36d0 

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SRUN0,
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
      REAL*8, INTENT(OUT) :: RUN0, SRUN0
      REAL*8 :: BYMSI2, FMSI1, FMSI2, FMSI3, FHSI1, FHSI2, FHSI3,
     *     CMPRS, SNWF, RAIN, FREZ1, MSI1, FSSI2,
     *     FSSI3, MELT1, SMELT1, DSNOW, SS12, FSSI1
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

      WETSNOW=.FALSE.
      RUN0=0. ; SRUN0=0.
#ifdef TRACERS_WATER
      TRUN0=0.
#endif
      MSI1=SNOW+ACE1I

      IF (PRCP.gt.0) THEN
      HSIL(1) = HSIL(1)+ENRGP  ! add total energy of precipitation

C**** Snowfall is calculated from precip energy (0 deg or colder)
      SNWF = MAX(0d0,MIN(PRCP,-ENRGP*BYLHM))
C**** Rain is remaining precip (0 deg or warmer)
      RAIN = PRCP-SNWF
      WETSNOW = RAIN.GT.1d-5*PRCP  ! i.e. a noticeable fraction of prec
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
C**** SALTI=SALTI - SMELT1        - FSSI2

C**** Calculate changes to snow
      IF (SNOW+SNWF.GT.MELT1) THEN ! some snow remains
        IF (SNWF.ge.PRCP .and. SNOW+SNWF .GT. SNOMAX) THEN
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
          CMPRS = SNOW+SNWF-0.9d0*SNOMAX
        ELSE
C**** RAIN and MELT COMRESSES SNOW INTO ICE
          CMPRS = MIN(dSNdRN*(RAIN+MELT1), SNOW+SNWF-MELT1)
        END IF
        DSNOW = SNWF - (MELT1+CMPRS)
        SMELT1=0.
      ELSE
C**** RAIN MELTS ALL SNOW AND SOME ICE
        CMPRS = 0.
        DSNOW = -SNOW
        SMELT1=(MELT1-SNWF-SNOW)*(SSIL(1)+SSIL(2))/ACE1I
      END IF
C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = SNWF+FREZ1-MELT1-DSNOW ! either up or down
      FMSI1 = XSI(1)*FMSI2+XSI(2)*(SNWF+FREZ1-MELT1)

#ifdef TRACERS_WATER
      TRSIL(:,1) = TRSIL(:,1)+TRPRCP(:)*(SNWF+FREZ1)/PRCP
      TRRMF(:) = TRPRCP(:)*(1.-(SNWF+FREZ1)/PRCP)
      TRMELT1(:) = MELT1*TRSIL(:,1)/(XSI(1)*MSI1+SNWF+FREZ1)
#endif

C***** Calculate consequential heat/salt flux between layers
      IF (FMSI2.gt.0) THEN      ! downward flux to thermal layer 3
        FHSI2 = FMSI2*HSIL(2)/(XSI(2)*MSI1)
        FSSI2 = FMSI2*SSIL(2)/(XSI(2)*MSI1)
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

C**** Fluxes between first two layers are complicated because salinty
C**** is evenly shared over the ice portion, and zero for the snow,
C**** while tracers have distinct concentrations with respect to the
C**** fresh water mass in both layers
      SS12=(SSIL(1)+SSIL(2))-FSSI2-SMELT1
      FSSI1=SSIL(1)-SMELT1-SS12*MAX(0d0,
     *     XSI(1)*ACE1I-(SNOW+DSNOW)*XSI(2))/ACE1I

      IF (FMSI1.gt.0) THEN      ! downward flux to thermal layer 2
        FHSI1 = FMSI1*HSIL(1)/(XSI(1)*MSI1+SNWF+FREZ1-MELT1)
      ELSE                      ! upward flux
        FHSI1 = FMSI1*HSIL(2)/(XSI(2)*MSI1)
      END IF

#ifdef TRACERS_WATER
C**** tracer flux depends on sign of freshwater flux
      IF (FMSI1-FSSI1.gt.0) THEN      ! downward flux to thermal layer 2
        FTRSI1(:) = (FMSI1-FSSI1)*(TRSIL(:,1)-TRMELT1(:))/(XSI(1)*MSI1
     *       +SNWF+FREZ1-MELT1-SSIL(1)+SMELT1)
      ELSE                      ! upward flux
        FTRSI1(:)=(FMSI1-FSSI1)*TRSIL(:,2)/(XSI(2)*MSI1-SSIL(2))
      END IF
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
      SNOW = MAX(0d0,SNOW+DSNOW)
      MSI1 = ACE1I + SNOW
      MSI2 = MSI2 + FMSI2

      HSIL(1)=HSIL(1)- FHSI1
      HSIL(2)=HSIL(2)+(FHSI1-FHSI2)
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3

C**** Note salinity is evenly shared over ice portion
      SSIL(1)=SSIL(1)- FSSI1-SMELT1
      SSIL(2)=SSIL(2)+(FSSI1-FSSI2)
      SSIL(3)=SSIL(3)+(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)+ FSSI3

#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,1)=TRSIL(:,1)- FTRSI1(:)-TRMELT1(:)
      TRSIL(:,2)=TRSIL(:,2)+(FTRSI1(:)-FTRSI2(:))
      TRSIL(:,3)=TRSIL(:,3)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)+ FTRSI3(:)
#endif

C**** Diagnostics for output
      RUN0 = MELT1+(RAIN-FREZ1) ! runoff mass to ocean
      IF (RUN0.lt.1d-13) RUN0=0. ! clean up roundoff errors
      SRUN0 = SMELT1             ! salt in runoff
c     HFLUX= 0                  ! energy of runoff (currently at 0 deg)

#ifdef TRACERS_WATER
      TRUN0(:)=TRMELT1(:)+TRRMF(:)  ! tracer flux to ocean
#endif
      END IF

      TSIL(1:2) = (HSIL(1:2)/(XSI(1:2)*MSI1)+LHM)*BYSHI ! ice temp L=1,2
      TSIL(3:4) = (HSIL(3:4)/(XSI(3:4)*MSI2)+LHM)*BYSHI ! ice temp L=3,4

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
      REAL*8, DIMENSION(NTM) :: FTRSI1,FTRSI2,FTRSI3,
     *     TRMELT1,TRMLET2,TRMELT3,TRMELT4,TRDEW,TRMELT2
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
      REAL*8 MSI1, MELT1, MELT2, MELT3, MELT4, SNMELT, DSNOW
      REAL*8 DEW, CMPRS, DEW1, DEW2, EVAP1
      REAL*8 SMELT1, SMELT2, SMELT3, SMELT4, FSSI1, FSSI2, FSSI3, SS12
      REAL*8 FMSI1, FMSI2, FMSI3, FHSI1, FHSI2, FHSI3
      REAL*8 HC1, HC2, HC3, HC4, HICE, HSNOW
      REAL*8 dF1dTI, dF2dTI, dF3dTI, dF4dTI, F1, F2, F3, FO

      FMSI2=0. ; MELT1=0. ; MELT2=0. ; MELT3=0. ; MELT4=0.
      FHSI1=0. ; FHSI2=0. ; FHSI3=0.
      SMELT1=0 ; SMELT2=0. ; SMELT3=0.; SMELT4=0.
      FSSI2=0. ; FSSI3=0.
#ifdef TRACERS_WATER
      FTRSI1=0. ; FTRSI2=0. ; FTRSI3=0.
      TRMELT1=0. ; TRMLET2=0. ; TRMELT3=0. ; TRMELT4=0.
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
      EVAP1= MAX(0d0,EVAP)      ! positive evaporation
      SS12=SSIL(1)+SSIL(2)      ! total upper layer salt (only in ice)

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
      DO N=1,NTM
C**** Take into account occasional sign differences when summing
C**** over multiple surface time steps
        IF (DEW2.gt.0. .and. TRDEW(N).gt.0) THEN
          TRSIL(N,2) = TRSIL(N,2)+TRDEW(N)
        ELSE
          TRSIL(N,1) = TRSIL(N,1)+TRDEW(N)
        END IF
      END DO
#endif
C**** CMPRS is amount of snow turned to ice during melting
C**** Calculate remaining snow and necessary mass flux using
C**** SNOW =SNOW -       SNMELT        + DEW1 - CMPRS
C**** ACE1I=ACE1I-(MELT1+MELT2-SNMELT) + DEW2 + CMPRS - FMSI2
C**** SS12=SS12- SMELT1 -SMELT2                       - FSSI2
C****
      IF (SNOW.GT.SNMELT+EVAP1) THEN ! some snow remains
        CMPRS = MIN(dSNdML*SNMELT, SNOW-EVAP1-SNMELT) ! > 0.
        DSNOW = - (SNMELT+EVAP1+CMPRS)
        SMELT2=(MELT1+MELT2-SNMELT)*SS12/ACE1I
        SMELT1=0.
      ELSE ! all snow and some ice melts or evaporates
        CMPRS = 0.
        DSNOW = -SNOW
        SMELT1=(MELT1+MIN(0d0,EVAP1-SNOW))*SS12/(ACE1I-MAX(0d0
     *       ,EVAP1-SNOW))
        SMELT2= MELT2*SS12/(ACE1I-MAX(0d0,EVAP1-SNOW))
      END IF
#ifdef TRACERS_WATER
C**** tracer melt is calculated using fresh water concentration
      TRMELT1(:)=(MELT1-SMELT1)*TRSIL(:,1)/(XSI(1)*MSI1+DEW1-SSIL(1))
      TRMELT2(:)=(MELT2-SMELT2)*TRSIL(:,2)/(XSI(2)*MSI1+DEW2-SSIL(2))
#endif
C**** Mass fluxes required to keep first layer ice = ACE1I
      FMSI2 = -DSNOW+DEW-MELT1-MELT2 ! either up or down
      FMSI1 = -XSI(1)*DSNOW+DEW1-MELT1

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
      IF (FMSI2.gt.0) THEN ! downward flux to thermal layer 3
        FHSI2 = FMSI2*HSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
        FSSI2 = FMSI2*SSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
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
C**** Fluxes between first two layers are complicated because salinty
C**** is evenly shared over the ice portion, and zero for the snow,
C**** while tracers have distinct concentrations with respect to the
C**** fresh water mass in both layers
      SS12=(SSIL(1)+SSIL(2))-FSSI2-SMELT1-SMELT2
      FSSI1=SSIL(1)-SMELT1-SS12*MAX(0d0,
     *     XSI(1)*ACE1I-(SNOW+DSNOW)*XSI(2))/ACE1I

      IF (FMSI1.gt.0) THEN ! downward flux to thermal layer 2
        FHSI1 = FMSI1*HSIL(1)/(XSI(1)*MSI1+DEW1-MELT1)
      ELSE                 ! upward flux
        FHSI1 = FMSI1*HSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
      END IF

#ifdef TRACERS_WATER
C**** tracer flux depends on sign of freshwater flux
      IF (FMSI1-FSSI1.gt.0) THEN ! downward flux to thermal layer 2
        FTRSI1(:) = (FMSI1-FSSI1)*(TRSIL(:,1)-TRMELT1(:))/(XSI(1)*MSI1
     *       +DEW1-MELT1-SSIL(1)+SMELT1)
      ELSE                      ! upward flux
        FTRSI1(:)=(FMSI1-FSSI1)*(TRSIL(:,2)-TRMELT2(:))/
     *       (XSI(2)*MSI1+DEW2-MELT2-SSIL(2)+SMELT2)
      END IF
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

C**** Apply fluxes
      SNOW = MAX(0d0,SNOW+DSNOW)
      MSI1 = ACE1I+ SNOW
      MSI2=MSI2-(MELT3+MELT4)+FMSI2-FMOC

      HSIL(1)=HSIL(1)- FHSI1
      HSIL(2)=HSIL(2)+(FHSI1-FHSI2)
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+ FHSI3

C**** Note salinity is evenly shared over ice portion
      SSIL(1)=SSIL(1)-SMELT1 - FSSI1
      SSIL(2)=SSIL(2)-SMELT2 +(FSSI1-FSSI2)
      SSIL(3)=SSIL(3)-SMELT3 +(FSSI2-FSSI3)
      SSIL(4)=SSIL(4)-SMELT4 + FSSI3

#ifdef TRACERS_WATER
C**** Apply the tracer fluxes
      TRSIL(:,1)=TRSIL(:,1)- TRMELT1(:)- FTRSI1(:)
      TRSIL(:,2)=TRSIL(:,2)- TRMELT2(:)+(FTRSI1(:)-FTRSI2(:))
      TRSIL(:,3)=TRSIL(:,3)- TRMELT3(:)+(FTRSI2(:)-FTRSI3(:))
      TRSIL(:,4)=TRSIL(:,4)- TRMELT4(:)+ FTRSI3(:)
#endif

C**** Calculate additional runoff output fluxes
      RUN = MELT1+MELT2+MELT3+MELT4 ! mass flux to ocean
      SRUN= SMELT1+SMELT2+SMELT3+SMELT4   ! salt flux to ocean
c HMELT currently assumed to be zero since melting is at 0 deg
      ERUN= SROX(2)        ! + HMELT
#ifdef TRACERS_WATER
      TRUN(:)=TRMELT1(:)+TRMELT2(:)+TRMELT3(:)+TRMELT4(:)
                                ! tracer flux to ocean
#endif
c**** Decide WETSNOW for albedo calculations and save MELT12 for
c**** possible melt ponds
      MELT12=MELT1+MELT2
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
      REAL*8, PARAMETER :: Z2OIX = 4.9,
     *                     BYZICX=1./(Z1I+Z2OIX)
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
      REAL*8, DIMENSION(NTM) :: FTRSI3,FTRSI4
      INTEGER N
#endif
      REAL*8 FMSI4, FHSI3, FHSI4, FSSI3, FSSI4, DHM, DSM, DTM, SS12
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
C**** Create new ice in partially ice-covered ocean
      IF (ACEFI.le.0) GO TO 250 ! go to no freezing case
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
CC    FMSI3 = -XSI(3)*ACEFI ! < 0.
CC    FMSI4 = -ACEFI
      FHSI3 = -HSIL(4)*ACEFI*(XSI(3)/XSI(4))/MSI2
      FSSI3 = -SSIL(4)*ACEFI*(XSI(3)/XSI(4))/MSI2
#ifdef TRACERS_WATER
      FTRSI3(:) = -TRSIL(:,4)*ACEFI*(XSI(3)/XSI(4))/MSI2
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
      MSI2 = MSI2+ACEFI ! new ice mass of physical layer 2
      GO TO 270
  240 CONTINUE
C**** NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      MSI1 = SNOW + ACE1I ! mass of layer 1
C**** MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACEFI))/(ROICE+DRSI) ! layer 2
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
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ROICE*(ACE1I+MSI2))) 
      IF ((ROICE*(ACE1I+MSI2)).gt.5.*RHOI) OPNOCN=0.  ! no leads for h>5
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
      IF (ROICE.gt.0) THEN
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ROICE*(ACE1I+MSI2)))    ! -BYZICX)
      IF ((ROICE*(ACE1I+MSI2)).gt.5.*RHOI) OPNOCN=0.  ! no leads for h>5
      IF (ROICE.gt.(1.-OPNOCN-1d-4)) THEN
        ROICEN = 1.-OPNOCN
        DRSI = ROICEN-ROICE ! +ve
        FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! >0 ice mass out of layer 4
        FHSI4 = HSIL(4)*FMSI4/(XSI(4)*MSI2)
        FSSI4 = SSIL(4)*FMSI4/(XSI(4)*MSI2)
        FHSI3 = HSIL(3)*FMSI4/MSI2 ! downward heat flux into layer 4 
        FSSI3 = SSIL(3)*FMSI4/MSI2 ! downward salt flux into layer 4
#ifdef TRACERS_WATER
        FTRSI4(:) = TRSIL(:,4)*FMSI4/(XSI(4)*MSI2)
        FTRSI3(:) = TRSIL(:,3)*FMSI4/MSI2
#endif
        MSI2=MSI2-FMSI4         ! new ice mass of second physical layer
C**** redistribute heat/salt/tracer from change of mass
        DHM=FHSI4/(MSI1+MSI2)
        HSIL(1:2)=(ROICE/ROICEN)*(HSIL(1:2)+DHM*XSI(1:2)*MSI1)
        HSIL(3)  =(ROICE/ROICEN)*(HSIL(3)+DHM*XSI(3)*MSI2 -FHSI3)
        HSIL(4)  =(ROICE/ROICEN)*(HSIL(4)+DHM*XSI(4)*MSI2 +FHSI3-FHSI4)
C**** salinity complicated since no salt is put in snow 
        DSM=FSSI4/(MSI1+MSI2)
        SS12=(ROICE/ROICEN)*(SSIL(1)+SSIL(2)+DSM*MSI1)
        IF (ACE1I.gt.XSI(2)*MSI1) THEN
          SSIL(1)=SS12*(ACE1I-XSI(2)*MSI1)/ACE1I
        ELSE
          SSIL(1)= 0.
        END IF
        SSIL(2)=SS12-SSIL(1)
        SSIL(3)  =(ROICE/ROICEN)*(SSIL(3)+DSM*XSI(3)*MSI2 -FSSI3)
        SSIL(4)  =(ROICE/ROICEN)*(SSIL(4)+DSM*XSI(4)*MSI2 +FSSI3-FSSI4)
#ifdef TRACERS_WATER
        DO N=1,NTM
          DTM=FTRSI4(N)/(MSI1+MSI2)
          TRSIL(N,1:2)=(ROICE/ROICEN)*(TRSIL(N,1:2)+DTM*XSI(1:2)*MSI1)
          TRSIL(N,3) = (ROICE/ROICEN)*(TRSIL(N,3)+DTM*XSI(3)*MSI2
     *         -FTRSI3(N))
          TRSIL(N,4) = (ROICE/ROICEN)*(TRSIL(N,4)+DTM*XSI(4)*MSI2
     *         +FTRSI3(N)-FTRSI4(N))
        END DO
#endif
C       SNOW = SNOW   ! snow thickness is conserved
        ROICE = ROICEN
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
      TSIL(1:2)=(HSIL(1:2)/(XSI(1:2)*MSI1)+LHM)*BYSHI ! temp. layer 1/2
      TSIL(3:4)=(HSIL(3:4)/(XSI(3:4)*MSI2)+LHM)*BYSHI ! temp. layer 3/4

      RETURN
      END SUBROUTINE ADDICE

      SUBROUTINE SIMELT(DT,ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,Tm,TFO,TSIL,
#ifdef TRACERS_WATER
     &     TRSIL,TRUN0,
#endif
     &     ENRGUSED,RUN0,SALT)
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

C**** Estimate DRSI
      DRSI=0.
      IF (ROICE.lt.1d-4) THEN
        DRSI=ROICE
      ELSEIF (POCEAN.gt.0) THEN
C**** Estimate lateral melt using parameterisation from Makyut/Steele
C**** (via C. Bitz): Rside=dt*pi/(floesize*eta)*(3d-6)*(delT)^(1.36)
        dtemp=MAX(Tm-TFO,0d0)
        DRSI=DT*SILMFAC*dtemp**SILMPOW
        IF (ROICE-DRSI.lt.1d-4) DRSI=ROICE
      END IF
C**** Remove DRSI amount of ice
      ENRGUSED=-DRSI*(HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4)) !  [J/m^2]
      RUN0=DRSI*(SNOW + ACE1I + MSI2) 
      SALT=DRSI*(SSIL(1)+SSIL(2)+SSIL(3)+SSIL(4))
#ifdef TRACERS_WATER
      TRUN0(:)=DRSI*(TRSIL(:,1)+TRSIL(:,2)+TRSIL(:,3)+TRSIL(:,4))
#endif
      ROICE=ROICE-DRSI
      IF (ROICE.lt.1d-10) THEN
        ROICE=0.                ! deal with possible round off err
C**** set defaults if no ice is left
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
        TRSIL(:,:)=0.
#endif
      END IF
C****
      RETURN
      END SUBROUTINE SIMELT

      SUBROUTINE SSIDEC(I0,J0,MSI1,MSI2,HSIL,SSIL,DT,
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
      REAL*8, INTENT(INOUT) :: MSI2
      REAL*8, INTENT(IN) :: DT, MSI1
!@var MFLUX,SFLUX,HFLUX mass, salt and heat flux arising from 
!@+   sea salinity decay
      REAL*8, INTENT(OUT) :: MFLUX,HFLUX,SFLUX
#ifdef TRACERS_WATER
!@var TRSIL tracer amount in ice layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,LMI) :: TRSIL
!@var TRFLUX tracer flux arising from sea salinity decay
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
C**** no tracer removed if pure salt is lost
c        DO N=1,NTM
c          DTRSI(N,1:2) = (DMSI(1:2)-DSSI(1:2))*TRSIL(N,1:2)/
c     *       (XSI(1:2)*MSI1)
c          TRSIL(N,1:2) = TRSIL(N,1:2)-DTRSI(N,1:2)
c        END DO
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
C**** no tracer removed if pure salt is lost
c          DTRSI(:,L) = (DMSI(L)-DSSI(L))*TRSIL(:,L)/(XSI(L)*MSI2)
c          TRSIL(:,L) = TRSIL(:,L)-DTRSI(:,L)
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
        TRFLUX(N)=SUM(DTRSI(N,:)) ! tracer flux to ocean
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
!@+  with  Sib=Si, Tib=Ti m>0,  or Sib = fsss Sm, Tib=Tb m<0        (4)
!@+ This routine iteratively solves for m, Tb, Sb and Sib, Tib using
!@+ Newton's method. The form of the latent heat and conductivity
!@+ can optionally include salinity effects.
!@+ The actual fluxes (positive down) are then:
!@+  mflux = m   ; sflux = 1d-3 m Sib (kg/m^2/s)
!@+  hflux = lam_i(Ti,Si) (Ti-Tb)/dh -m Lh(Tib,Sib)+m shw Tib  (J/m^2/s)
!@+
      implicit none
!@var alamdS salinity coefficient for conductivity (from common?)
      real*8, parameter :: alamdS=0.117d0

!@var G_mole_T,G_mole_S molecular diffusion terms for heat/salt
C****  G_mole_T = 12.5 * 13.8d0**(2d0/3d0) - 6. = 65.9d0
C****  G_mole_S = 12.5 * 2432d0**(2d0/3d0) - 6. = 2255d0
      real*8, parameter :: G_mole_T = 65.9d0 , G_mole_S = 2255d0
!@var Si,Sm salinity in lowest ice layer and mixed layer (psu)
!@var Ti,Tm temperatures in ice and mixed layer (C)
!@var dh distance from center of bottom ice layer to base of ice (m)
      real*8, intent(in) :: Ti,Si,Tm,Sm,dh
!@var coriol Corilos parameter (used in turbulent flux calc) (1/s)
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
      real*8 :: m,Tb,Sib,Tib
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
      rsg=rhow*shw*g_T    ! /(1.+alpha*dtsrc*rhow*shw*g_T/mlsh)
c no salinity effects
      alamdh = alami/(dh+alpha*dtsrc*alami*byshi/(2.*dh*rhoi))
c S thermo
c      alamdh = (alami*Ti+alamdS*Si)/(dh*Ti+alpha*dtsrc*
c     *          (alami*Ti+alamdS*Si)*byshi/(2.*rhoi*dh))

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
        Trib(:)=tralpha(:)*rhow*g_S*Trm(:)/(rhow*g_S+(m-sflux)*(1.
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
C**** Assume constant g_T = 5d-5, g_S = 0.04 * g_T m/s
!@var rsg = rhow * shw * g_T turbulent energy flux (J/m^2 s)
!@var rgS = rhow * g_S turbulent tracer flux (kg/m^2 s)
      real*8, parameter ::  rsg = rhow*shw*5d-5, rgS=rhow*2d-8
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
c**** solar input by the approximate co-albedo. Could possibly be
c**** retreived from radiation code
c**** VIS: 250 - 690 nm,       NIR1= 690 - 1119 nm
c****         52%,     and          33% of incoming solar
c**** NIR2/3 (> 1119 nm assumed not to be transmitted)
c****
      if (hsnow.gt.0.02) then    ! same cutoff as for albedo
c**** Extinction coefficients for snow depend on wetness
c**** Since albedo does as well, absorbed fractions also vary
        if (wetsnow) then
          fracvis =0.22d0 ;  fracnir1=0.32d0
          ksextvis=10.7d0 ; ksextnir1=118.d0
        else
          fracvis =0.07d0 ;  fracnir1=0.30d0
          ksextvis=19.6d0 ; ksextnir1=196.d0
        end if
      else
        fracvis=0.26d0 ; fracnir1=0.42d0
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
      if (present(press)) pr=press
C**** linear approximation
      tfrez = -mu*sss + dtdp*pr
C**** UNESCO formula (1983)
c     tfrez = (a01 + a02*sss)*sss + a03*sss*sqrt(sss)+ dtdp*pr

      return
      end function tfrez

      END MODULE SEAICE

      MODULE SEAICE_COM
!@sum  SEAICE_COM contains the model arrays for seaice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      USE SEAICE, only : lmi

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
      REAL*8, DIMENSION(NTM) :: TRSI0
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
     *     ,'(a8,i3,a1,i3,a)')'R8 TRSI(',ntm,',',lmi,',im,jm)'
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
      USE GEOM, only : imaxj
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE SEAICE, only : lmi,xsi,ace1i
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
#ifdef TRACERS_WATER
     *     ,trsi
#endif
      USE FLUXES
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

C**** Check for NaN/INF in ice data
      CALL CHECK3(RSI,IM,JM,1,SUBR,'rsi')
      CALL CHECK3(MSI,IM,JM,1,SUBR,'msi')
      CALL CHECK3(HSI,LMI,IM,JM,SUBR,'hsi')
      CALL CHECK3(SSI,LMI,IM,JM,SUBR,'ssi')
      CALL CHECK3(SNOWI,IM,JM,1,SUBR,'sni')

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
     *             ,J),MSI(I,J),SNOWI(I,J),RSI(I,J)
              if (HSI(L,I,J).gt.0.or.TICE.gt.0.01.or.TICE.lt.-80.)
     *             QCHECKI = .TRUE.
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
        END DO
      END DO

#ifdef TRACERS_WATER
      do n=1,ntm
C**** check negative tracer mass        
        if (t_qlimit(n)) then
        do j=1,jm
          do i=1,imaxj(j)
            if (rsi(i,j).gt.0) then
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
          do j=1,jm
          do i=1,imaxj(j)
            if (rsi(i,j).gt.0) then
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

