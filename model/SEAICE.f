      MODULE SEAICE
!@sum  SEAICE contains all the sea ice related subroutines
!@auth Original Development Team
!@ver  1.0
!@cont PREC_SI,SEA_ICE,ADDICE,SIMELT
      USE CONSTANT, only : lhm,rhoi,byrhoi,rhow,shi,shw,byshi
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
!@param KEXT extinction coefficient for light in sea ice (1/m)
      REAL*8, PARAMETER :: KEXT = 1.5d0
!@var FLEADOC lead fraction for ocean ice (%)
      REAL*8, PARAMETER :: FLEADOC = 0.06d0
!@var FLEADLK lead fraction for lakes (%)
      REAL*8, PARAMETER :: FLEADLK = 0.
!@param BYRLI,BYRLS reciprical of density*lambda
      REAL*8, PARAMETER :: BYRLI = 1./(RHOI*ALAMI),
     *     BYRLS = 1./(RHOS*ALAMS)

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,HSIL,TSIL,PRCP,ENRGP,TFO,RUN0,FMOC
     *     ,FHOC,QFIXR)
!@sum  PREC_SI Adds the precipitation to sea/lake ice
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@param SNOMAX maximum allowed snowdepth (1m equivalent) (kg/m^2)
!@param dSNdRN
      REAL*8, PARAMETER :: SNOMAX=1d0*RHOS, dSNdRN=0.
!@var ENRGP total energy of precip (C),(J/m^2)
!@var PRCP amount of precip (kg/m^2)
      REAL*8, INTENT(IN) :: ENRGP, PRCP
!@var QFIXR true if RSI and MSI2 are fixed (ie. for a fixed SST run)
      LOGICAL, INTENT(IN) :: QFIXR
!@var SNOW snow mass (kg/m^2)
!@var MSI1 first layer ice mass (= SNOW + ACE1I) (kg/m^2)
!@var MSI2 second layer ice mass (kg/m^2)
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
!@var HSIL enthalpy of ice layers (J/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL
!@var TSIL temperature of ice layers (J/m^2) (DIAGNOSTIC ONLY)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var TFO temperature of freezing of water below ice (C) 
      REAL*8, INTENT(IN) :: TFO
!@var FHOC implied energy flux at base to keep fixed ice mass (J/m^2)
      REAL*8, INTENT(OUT) :: FHOC
!@var FMOC  implied ice mass flux at base to keep fixed ice (kg/m^2)
      REAL*8, INTENT(OUT) :: FMOC
!@var RUN0 runoff from ice (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0
      REAL*8 :: BYMSI2, FMSI1, FMSI2, FHSI1, FHSI2, FHSI3, CMPRS, SNWF,
     *     MELT1, RAIN, FREZ1, MSI1

C**** initialize fluxes
      RUN0=0. ; FMSI2=0. ; FHSI1=0 ; FHSI2=0. ; FHSI3=0. ;
      FMOC=0. ; FHOC=0.
C**** reciprocal of ice thickness for efficiency
      BYMSI2=1./MSI2
      MSI1=SNOW+ACE1I

C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSIL(1) = HSIL(1)+ENRGP  ! add total energy of precipitation
      IF (ENRGP.LE. -PRCP*LHM) GO TO 180
      IF (ENRGP.LE. 0.) GO TO 160
C**** ALL PRECIPITATION IS RAIN ABOVE 0degC
C**** RAIN COMPRESSES SNOW INTO ICE
      RAIN = PRCP
      IF (HSIL(1)/LHM+XSI(1)*MSI1 .LE. 0.) GO TO 140
C**** WARM RAIN MELTS SOME SNOW OR ICE
      MELT1 = HSIL(1)/LHM+XSI(1)*MSI1 ! melted snow and ice (kg/m^2)
      RUN0=MELT1 + PRCP ! water mass flux to ocean (kg/m^2)
      IF (MSI1-ACE1I .LE. MELT1) GO TO 130
C**** RAIN MELTS SOME SNOW AND COMPRESSES SNOW INTO ICE
      FMSI2 = MIN (dSNdRN*(RAIN+MELT1), MSI1-ACE1I-MELT1) ! > 0.
      FHSI1 = -LHM*(XSI(1)*FMSI2-XSI(2)*MELT1) ! downward heat flux
c     F1 = HSIL(1)*(XSI(1)*FMSI2-XSI(2)*MELT1)/(XSI(1)*MSI1-MELT1)
      FHSI2 = HSIL(2)*FMSI2*BYXSI(2)/MSI1 ! > 0.
      FHSI3 = HSIL(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1-MELT1-FMSI2
      GO TO 210
  130 CONTINUE
C**** RAIN MELTS ALL SNOW AND SOME ICE
      FMSI2 = MSI1-ACE1I-MELT1 ! < 0.(upward ice mass flux)
C     FMSI1 = XSI(1)*(MSI1-ACE1I)-MELT1 ! < 0.(melted snow/ice mass)
      FHSI1 = HSIL(2)*((XSI(1)/XSI(2))*FMSI2-MELT1)/MSI1 !  upw heat flx
      FHSI2 = HSIL(3)*FMSI2*BYXSI(3)*BYMSI2 !  upward heat flux -> lyr 2
      FHSI3 = HSIL(4)*FMSI2*BYMSI2 !  upward heat flux into layer 3
      MSI1 = ACE1I ! Keep the first layer ice mass constant ACE1I
      GO TO 210
  140 CONTINUE
C**** RAIN COMPRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      CMPRS = MIN(dSNdRN*RAIN, MSI1-ACE1I) ! all snow or part of rain
      IF (-HSIL(1)/LHM-XSI(1)*MSI1 .LT. RAIN) GO TO 150
C**** ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
C     FMSI1 = XSI(1)*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSIL(1)*(XSI(1)*CMPRS+RAIN)/(XSI(1)*MSI1+RAIN) ! downward
      FHSI2 = HSIL(2)*FMSI2*BYXSI(2)/MSI1 ! downw heat flux from layer 2
      FHSI3 = HSIL(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1-CMPRS ! first layer ice mass
      GO TO 210
  150 CONTINUE
C**** JUST PART OF RAIN FREEZES IN LAYER 1
      FREZ1 = -HSIL(1)/LHM-XSI(1)*MSI1 ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI(1)*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI(1)*CMPRS+FREZ1) ! downw heat flux from layer 1
      FHSI2 = HSIL(2)*FMSI2*BYXSI(2)/MSI1 ! downw heat flux from layer 2
      FHSI3 = HSIL(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1-CMPRS ! first layer ice mass
      GO TO 210
  160 CONTINUE
C**** PRECIPITATION IS A MIXTURE OF RAIN AND SNOW AT 0 degC
C**** RAIN COMRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      SNWF = -ENRGP/LHM ! snow fall
      RAIN = PRCP-SNWF  ! rain fall
      CMPRS = MIN(dSNdRN*RAIN, MSI1+SNWF-ACE1I) ! compression
      IF (-HSIL(1)/LHM-XSI(1)*MSI1-SNWF .LT. RAIN) GO TO 170
C**** ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
      FMSI1 = XSI(1)*CMPRS+XSI(2)*SNWF+RAIN ! downward ice mass flux
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSIL(1)*FMSI1/(XSI(1)*MSI1+PRCP) ! downward heat flux
      FHSI2 = HSIL(2)*FMSI2*BYXSI(2)/MSI1 ! downw heat flux from layer 2
      FHSI3 = HSIL(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1+SNWF-CMPRS ! first layer ice mass
      GO TO 210
  170 CONTINUE
C**** NOT ALL RAIN FREEZES IN LAYER 1
      FREZ1 = -HSIL(1)/LHM-XSI(1)*MSI1-SNWF ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI(1)*CMPRS+XSI(2)*SNWF+FREZ1 ! downward ice mass flux
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI(1)*CMPRS+XSI(2)*SNWF+FREZ1) ! downward
      FHSI2 = HSIL(2)*FMSI2*BYXSI(2)/MSI1 ! downw heat flux from layer 2
      FHSI3 = HSIL(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1+SNWF-CMPRS ! first layer ice mass
      GO TO 210
  180 CONTINUE
C**** ALL PRECIPITATION IS SNOW, SNOW AMOUNT INCREASES
      IF (MSI1+PRCP .LE. SNOMAX+ACE1I) THEN
C     FMSI1 = XSI(2)*PRCP ! > 0.(snow fall to layer 1)
      FHSI1 = HSIL(1)*XSI(2)*PRCP/(XSI(1)*MSI1+PRCP) ! downwrd heat flux
      MSI1 = MSI1+PRCP ! first layer ice mass
      ELSE
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
      FMSI2 = MSI1+PRCP-(0.9*SNOMAX+ACE1I) ! > 0.(compressed snow)
C     FMSI1 = XSI(2)*PRCP+XSI(1)*FMSI2 ! > 0.(downward ice mass flux)
      FHSI1 = HSIL(1)*(XSI(2)*PRCP+XSI(1)*FMSI2)/(XSI(1)*MSI1+PRCP)
      FHSI2 = HSIL(2)*FMSI2*BYXSI(2)/MSI1 ! downw heat flux from layer 2
      FHSI3 = HSIL(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = 0.9d0*SNOMAX+ACE1I ! first layer ice mass
      END IF
  210 CONTINUE
C**** Adjust amounts resulting from fluxes.
      SNOW = MAX(0d0,MSI1-ACE1I)         ! snow mass
      HSIL(1) = HSIL(1)- FHSI1
      HSIL(2) = HSIL(2)+(FHSI1-FHSI2)
      HSIL(3) = HSIL(3)+(FHSI2-FHSI3)
      IF (QFIXR) THEN           ! Calc. implicit fluxes for fixed ice
        FMOC = FMSI2
        IF (FMOC.gt.0) THEN
          FHOC = FMOC*HSIL(4)/(MSI2*XSI(4))
        ELSE
          FHOC = FMOC*(SHI*TFO-LHM)
        END IF
        HSIL(4) = HSIL(4)+FHSI3-FHOC*XSI(4)
        HSIL(3) = HSIL(3)-FHOC*XSI(3)
      ELSE                      ! ADVECT ICE for predicted ice
        MSI2 = MSI2+FMSI2       ! 2nd layer sea ice mass (kg/m^2)
        HSIL(4) = HSIL(4)+FHSI3
      END IF

      TSIL(1) = (HSIL(1)/(XSI(1)*MSI1)+LHM)*BYSHI ! ice temperature L=1
      TSIL(2) = (HSIL(2)/(XSI(2)*MSI1)+LHM)*BYSHI ! ice temperature L=2
      TSIL(3) = (HSIL(3)/(XSI(3)*MSI2)+LHM)*BYSHI ! ice temperature L=3
      TSIL(4) = (HSIL(4)/(XSI(4)*MSI2)+LHM)*BYSHI ! ice temperature L=4

      RETURN
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,HSIL,MSI2,F0DT,F1DT,EVAP
     *     ,SROX,TFO,RUN0,FMOC,FHOC,FO,QFIXR,QFLUXLIM,FLUXLIM)
!@sum  SEA_ICE applies surface fluxes to ice covered areas
!@auth Gary Russell
!@ver  1.0

      IMPLICIT NONE

      REAL*8, PARAMETER :: ALPHA = 1.0, dSNdML =0.
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
!@var QFIXR  true if RSI and MSI2 are fixed (ie. for fixed SST run)
!@var QFLUXLIM true if the flux at base of ice is limited
      LOGICAL,INTENT(IN) :: QFIXR, QFLUXLIM
!@var FLUXLIM limit of water-ice flux if QFLUXLIM is true 
      REAL*8, INTENT(IN) :: FLUXLIM
!@var TFO freezing temperature of water below ice (C)
      REAL*8, INTENT(IN) :: TFO

      REAL*8, INTENT(IN) :: ROICE
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
      REAL*8, DIMENSION(LMI) :: TSIL
      REAL*8, INTENT(OUT) :: FO, RUN0, FMOC, FHOC
      REAL*8 MSI1, MELT1, MELT2, MELT3, MELT4, SNMELT, DSNOW
      REAl*8 DEW, CMPRS, DEW1, DEW2, EVAP1
      REAL*8 FMSI1, FMSI2, FMSI3, FMSI4, FHSI1, FHSI2, FHSI3, FHSI4
      REAL*8 HC1, HC2, HC3, HC4, HICE2
      REAL*8 dF1dTI, dF2dTI, dF3dTI, dF4dTI, F1, F2, F3
C**** Initiallise output
      FO=0. ; RUN0=0.  ; FMOC=0. ; FHOC=0.

      IF (ROICE .EQ. 0.) RETURN
      FMSI2=0. ; MELT1=0. ; MELT2=0. ; MELT3=0. ; MELT4=0.
      FHSI1=0. ; FHSI2=0. ; FHSI3=0.  
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C**** Calculate solar fractions
      IF (SROX(1).gt.0) THEN
c       FSRI(1) = EXP(-KEXT*HICE1)  ! included in F1DT
        HICE2=ACE1I*BYRHOI+SNOW/RHOS
        FSRI(2) = EXP(-KEXT*HICE2)
        FSRI(3) = EXP(-KEXT*(HICE2+MSI2*XSI(3)*BYRHOI))
        FSRI(4) = EXP(-KEXT*(HICE2+MSI2*BYRHOI))
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
      FO=(dF4dTI*(HC4*(TSIL(4)-TFO)+ALPHA*F3)+HC4*SROX(1)*FSRI(4))/
     *     (HC4+ALPHA*dF4dTI)
      SROX(2) = SROX(1)*FSRI(4)
C**** limit flux if necessary
      IF (QFLUXLIM.and.FO.lt.FLUXLIM) FO = FLUXLIM

      HSIL(1) = HSIL(1)+(F0DT-F1DT)
      HSIL(2) = HSIL(2)+(F1DT-F2)
      HSIL(3) = HSIL(3)+(F2-F3)
      HSIL(4) = HSIL(4)+(F3-FO)

      DEW = -EVAP ! dew/evap to the surface
      EVAP1= MAX(0d0,EVAP)  ! positive evaporation

C**** Add DEW directly to ice (which can be in either layer)
C**** DEW1 adds to layer 1, DEW2 adds to layer 2
C**** Calculate melting in first two thermal layers
C**** SNMELT is the melting that is applied to the snow first
      IF (SNOW*XSI(1).gt.XSI(2)*ACE1I) THEN ! first layer is all snow
        DEW1 = -EVAP1           ! <0 i.e. evaporation from snow
        DEW2 =  MAX(0d0,DEW)    ! >0 i.e. dew to second layer ice
        MELT1 = MAX(0d0,HSIL(1)/LHM+XSI(1)*MSI1+DEW1)
        MELT2 = MAX(0d0,HSIL(2)/LHM+XSI(2)*MSI1+DEW2)
        SNMELT=MELT1+MELT2   
      ELSE  ! first layer is snow and some ice
        DEW1 = DEW   ! all fluxes to first layer
        DEW2 = 0.
        MELT1 = MAX(0d0,HSIL(1)/LHM+XSI(1)*MSI1+DEW1)
        MELT2 = MAX(0d0,HSIL(2)/LHM+XSI(2)*MSI1)
        SNMELT=MELT1
      END IF

C**** Check for melting in levels 3 and 4
      MELT3 = MAX(0d0,HSIL(3)/LHM+XSI(3)*MSI2) 
      MELT4 = MAX(0d0,HSIL(4)/LHM+XSI(4)*MSI2) 

C**** CMPRS is amount of snow turned to ice during melting
C**** Calculate remaining snow and necessary mass flux using
C**** SNOW =SNOW -       SNMELT        + DEW1 - CMPRS
C**** ACE1I=ACE1I-(MELT1+MELT2-SNMELT) + DEW2 + CMPRS - FMSI2
C****
      IF (SNOW.GT.SNMELT+EVAP1) THEN ! some snow remains
        CMPRS = MIN(dSNdML*SNMELT, SNOW-EVAP1-SNMELT) ! > 0.
        DSNOW = - (SNMELT+EVAP1+CMPRS)
      ELSE ! all snow and some ice melts
        CMPRS = 0.
        DSNOW = -SNOW
      END IF
C**** Mass fluxes required to keep first layer ice = ACE1I 
      FMSI2 = -DSNOW+DEW-MELT1-MELT2 ! either up or down
      FMSI1 = -XSI(1)*DSNOW+DEW1-MELT1-CMPRS

C***** Calculate consequential heat flux between layers
      IF (FMSI1.gt.0) THEN ! downward flux to thermal layer 2
        FHSI1 = FMSI1*HSIL(1)/(XSI(1)*MSI1+DEW1-MELT1)
      ELSE                 ! upward flux
        FHSI1 = FMSI1*HSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
      END IF
      IF (FMSI2.gt.0) THEN ! downward flux to thermal layer 3
        FHSI2 = FMSI2*HSIL(2)/(XSI(2)*MSI1+DEW2-MELT2)
      ELSE                 ! upward flux
        FHSI2 = FMSI2*HSIL(3)/(XSI(3)*MSI2-MELT3)
      END IF

C**** Calculate mass/heat flux at base if required
      IF (QFIXR) THEN ! fixed sea ice, calculate implicit base fluxes
        FMSI4 = FMSI2 - MELT3 - MELT4 
        IF (FMSI4.gt.0) THEN ! downward flux to ocean
          FHSI4 = FMSI4*HSIL(4)/(XSI(4)*MSI2-MELT4)
        ELSE ! upward flux, at freezing point of ocean
          FHSI4 = FMSI4*(SHI*TFO-LHM)
        END IF
      ELSE   ! coupled models allow MSI2 to melt
        FMSI4 = 0. 
        FHSI4 = 0.
      END IF

C**** Calculate mass/heat flux between layers 3 and 4
      FMSI3 = XSI(3)*(FMSI4+MELT4) + XSI(4)*(FMSI2-MELT3) 
      IF (FMSI3.gt.0) THEN ! downward flux to layer 4
        FHSI3 = FMSI3*HSIL(3)/(XSI(3)*MSI2-MELT3)
      ELSE                 ! upward flux
        FHSI3 = FMSI3*HSIL(4)/(XSI(4)*MSI2-MELT4)
      END IF

C**** Apply fluxes
      HSIL(1)=HSIL(1)- FHSI1
      HSIL(2)=HSIL(2)+(FHSI1-FHSI2)
      HSIL(3)=HSIL(3)+(FHSI2-FHSI3)
      HSIL(4)=HSIL(4)+(FHSI3-FHSI4)
      SNOW = MAX(0d0,SNOW+DSNOW)
      MSI1 = ACE1I+ SNOW
C**** MSI2 = MSI2 -(MELT3+MELT4)+(FMSI2-FMSI4)
      IF (.not. QFIXR) MSI2=MSI2-(MELT3+MELT4)+FMSI2

C**** Calculate output diagnostics
      RUN0 = MELT1+MELT2+MELT3+MELT4  ! mass flux to ocean
      FMOC = FMSI4                    ! implicit mass flux 
      FHOC = FHSI4                    ! implicit heat flux
C****
      RETURN
      END SUBROUTINE SEA_ICE

      SUBROUTINE ADDICE (SNOW,ROICE,HSIL,MSI2,TSIL,ENRGFO
     *     ,ACEFO,ACE2F,ENRGFI,FLEAD,QFIXR)
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
      REAL*8, INTENT(IN) ::  ENRGFI, ENRGFO, ACEFO, ACE2F
!@var ROICE,SNOW,MSI1,MSI2
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL

      REAl*8 FMSI1, FMSI2, FMSI3, FMSI4, FHSI1, FHSI2, FHSI3, FHSI4
      REAL*8 ROICEN, OPNOCN, DRSI, DRI, MSI1

      MSI1=SNOW+ACE1I
      IF (.not. QFIXR) THEN
      IF (ROICE.LE.0. .and. ACEFO.gt.0) THEN
        ROICE=ACEFO/(ACE1I+AC2OIM)
        SNOW=0.
C****   TSIL=(ENRGFO/ACEFO + LHM)*BYSHI ! but hsi is primary var.
        HSIL(1) = (ENRGFO/ACEFO)*XSI(1)*ACE1I
        HSIL(2) = (ENRGFO/ACEFO)*XSI(2)*ACE1I
        HSIL(3) = (ENRGFO/ACEFO)*XSI(3)*AC2OIM
        HSIL(4) = (ENRGFO/ACEFO)*XSI(4)*AC2OIM
        MSI1=ACE1I
        MSI2=AC2OIM
      ELSEIF (ROICE.gt.0) THEN

      IF (ACE2F.le.0) GO TO 250 ! go to no freezing case
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
C     FMSI3 = -XSI(3)*ACE2F ! < 0.
C     FMSI4 = -ACE2F
      FHSI3 = -HSIL(4)*ACE2F*(XSI(3)/XSI(4))/MSI2
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      IF (ACEFO .GT. 0.) GO TO 240
C**** NEW ICE IS FORMED BELOW OLD SEA ICE
      HSIL(3) = HSIL(3)-FHSI3
      HSIL(4) = HSIL(4)+(FHSI3+ENRGFI)
      MSI2 = MSI2+ACE2F ! new ice mass of physical layer 2
      GO TO 270
  240 CONTINUE
C**** NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      MSI1 = SNOW + ACE1I ! mass of layer 1
C**** MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACE2F))/(ROICE+DRSI) ! layer 2
      HSIL(1) = ((1.-ROICE)*ENRGFO*YSI(1)+ROICE*HSIL(1))/(ROICE+DRSI)
      HSIL(2) = ((1.-ROICE)*ENRGFO*YSI(2)+ROICE*HSIL(2))/(ROICE+DRSI)
      HSIL(3) = ((1.-ROICE)*ENRGFO*YSI(3)+ROICE*(HSIL(3)-FHSI3))/
     A       (ROICE+DRSI)
      HSIL(4) = ((1.-ROICE)*ENRGFO*YSI(4)+ROICE*(HSIL(4)+FHSI3+ENRGFI))/
     A       (ROICE+DRSI)
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
      HSIL(1) = ((1.-ROICE)*ENRGFO*YSI(1)+ROICE*HSIL(1))/(ROICE+DRSI)
      HSIL(2) = ((1.-ROICE)*ENRGFO*YSI(2)+ROICE*HSIL(2))/(ROICE+DRSI)
      HSIL(3) = ((1.-ROICE)*ENRGFO*YSI(3)+ROICE*HSIL(3))/(ROICE+DRSI)
      HSIL(4) = ((1.-ROICE)*ENRGFO*YSI(4)+ROICE*HSIL(4))/(ROICE+DRSI)
      ROICE = ROICE+DRSI ! new ice concentration
      END IF
  270 CONTINUE
C**** COMPRESS THE ICE HORIZONTALLY IF TOO THIN OR LEAD FRAC. TOO SMALL
      OPNOCN=MIN(0.1d0,FLEAD*RHOI/(ACE1I+MSI2))   ! -BYZICX) sometime -ve!
      IF (MSI2.LT.AC2OIM .or. ROICE.GT.1.-OPNOCN) THEN 
      ROICEN = MIN(ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM),1.-OPNOCN)
      DRSI = ROICEN-ROICE ! < 0. compressed ice concentration
      DRI = ROICE-ROICEN ! > 0. for diagnostics
C     FMSI3 = XSI(3)*FMSI4 ! < 0. upward ice mass flux into layer 3
      FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! upward ice mass into layer 4
      FHSI3 = HSIL(4)*FMSI4*(XSI(3)/XSI(4))/MSI2 ! upward heat flux
      FHSI4 = (HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4))*(DRSI/ROICEN)
      HSIL(3) = HSIL(3)-FHSI3
      HSIL(4) = HSIL(4)+(FHSI3-FHSI4)
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

      SUBROUTINE SIMELT(ROICE,SNOW,MSI2,HSIL,TSIL,ENRGW,ENRGUSED,RUN0)
!@sum  SIMELT melts sea ice if surrounding water is warm
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
!@var ENRGW energy available to melt ice (J/m^2)
      REAL*8, INTENT(IN) :: ENRGW
!@var ROICE,SNOW,MSI2 ice variables (%,kg/m^2,kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var HSIL ice enthalpy  (J/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSIL
!@var TSIL ice enthalpy  (J/m^2)
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL
!@var ENRGUSED energy used to melt ice (J/m^2)
      REAL*8, INTENT(OUT) :: ENRGUSED
!@var RUN0 amount of sea ice melt (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0
      REAL*8, DIMENSION(LMI) :: HSI
      REAL*8 MSI1,MELT,DRSI,ROICEN,FHSI4,FHSI3,ENRGI
c      REAL*8 E_BOTTOM,GAMMA,HCRIT,HICE,ACE

      MSI1 = SNOW + ACE1I
      ENRGI = HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4) ! energy in seaice [J/m^2]
      IF (ROICE*ENRGI+ENRGW.GE.0. .or. ROICE.lt.1d-4) THEN
C**** THE WARM OCEAN MELTS ALL THE SNOW AND ICE
      ENRGUSED=-ROICE*ENRGI     ! only some energy is used
      RUN0=ROICE*(MSI1+MSI2)    ! all ice is melted
      ROICE=0.
      SNOW=0.
      MSI2=AC2OIM
      HSIL(1:2)=-LHM*XSI(1:2)*ACE1I
      HSIL(3:4)=-LHM*XSI(3:4)*AC2OIM
      TSIL=0.
      RETURN
      END IF
C**** THE WARM OCEAN COOLS TO 0 DEGREES MELTING SOME SNOW AND ICE
C**** Reduce the ice depth
 230  ENRGUSED=ENRGW            ! all energy is used
C**** This is not currently used, but it could be
C**** I think this prescription is wrong though
c      ACE = SNOW + ACE1I + MSI2 ! total ice mass
c      HCRIT = 1.                ! critical ice thickness [m]
c      HICE = ACE*BYRHOI           ! sea ice thickness [m]
c      GAMMA = MAX(1d0,HICE/HCRIT)
c      E_BOTTOM = ENRGW*(ROICE/(1.+GAMMA)) ! for bottom melting
C**** MELT ICE VERTICALLY, AND THEN HORIZONTALLY
      MELT = -XSI(4)*MSI2*ENRGW/HSIL(4) ! melted ice at the bottom
c      MELT = -XSI(4)*MSI2*E_BOTTOM/HSIL(4) ! melted ice at the bottom
      IF (MSI2-MELT .LT. AC2OIM) MELT = MAX(MSI2-AC2OIM,0d0)
C     FMSI3 = XSI(3)*MELT ! > 0.
      FHSI3 = HSIL(3)*MELT/MSI2
      FHSI4 = HSIL(4)*MELT*BYXSI(4)/MSI2
C**** MELT SOME ICE HORIZONTALLY WITH REMAINING ENERGY
      DRSI = (ENRGW+ROICE*FHSI4)/(HSIL(1)+HSIL(2)+HSIL(3)+HSIL(4)-FHSI4)
      ROICEN = MIN(ROICE+DRSI,1d0)  ! new sea ice concentration (DRSI<0)
c     SNOW = SNOW*(ROICE/ROICEN) ! new snow ice mass
      MSI2 = MSI2-MELT
      HSIL(3) = HSIL(3)-FHSI3
      HSIL(4) = HSIL(4)+FHSI3-FHSI4
C**** CALCULATE SEA ICE TEMPERATURE (FOR OUTPUT ONLY)
      TSIL(1) = (HSIL(1)/(XSI(1)*MSI1) +LHM)*BYSHI ! temperature of L=1
      TSIL(2) = (HSIL(2)/(XSI(2)*MSI1) +LHM)*BYSHI ! temperature of L=2
      TSIL(3) = (HSIL(3)/(XSI(3)*MSI2) +LHM)*BYSHI ! temperature of L=3
      TSIL(4) = (HSIL(4)/(XSI(4)*MSI2) +LHM)*BYSHI ! temperature of L=4
      RUN0=ROICE*MELT
      ROICE=ROICEN
      RETURN
      END SUBROUTINE SIMELT

      END MODULE SEAICE

      MODULE SEAICE_COM
!@sum  SEAICE_COM contains the model arrays for seaice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE SEAICE, only : lmi

      IMPLICIT NONE
!@var RSI fraction of open water area covered in ice
      REAL*8, DIMENSION(IM,JM) :: RSI
!@var SNOWI snow amount on sea ice (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: SNOWI
!@var MSI thickness of ice second layer (layer 1=const) (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: MSI
!@var HSI enthalpy of each ice layer (J/m^2)
      REAL*8, DIMENSION(LMI,IM,JM) :: HSI
!@var SSI sea ice salt content (kg/m^2)
      REAL*8, DIMENSION(LMI,IM,JM) :: SSI
!@var TSI temperature of each ice layer (J/m^2) (DIAGNOSTIC ONLY)
      REAL*8, DIMENSION(LMI,IM,JM) :: TSI

      END MODULE SEAICE_COM

      SUBROUTINE io_seaice(kunit,iaction,ioerr)
!@sum  io_seaice reads and writes seaice variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite
      USE SEAICE_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "SICE01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
         WRITE (kunit,err=10) MODULE_HEADER,RSI,HSI,SNOWI,MSI,SSI
      CASE (IOREAD:)            ! input from restart file
         READ (kunit,err=10) HEADER,RSI,HSI,SNOWI,MSI,SSI
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_seaice

      SUBROUTINE CHECKI(SUBR)
!@sum  CHECKI Checks whether Ice values are reasonable
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : lhm,shi
      USE MODEL_COM
      USE SEAICE, only : lmi,xsi,ace1i
      USE SEAICE_COM, only : rsi,msi,hsi,snowi
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
      CALL CHECK3(HSI,4,IM,JM,SUBR,'hs')
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
            IF (HSI(L,I,J).gt.0.or.TICE.gt.1d-8.or.TICE.lt.-80.) THEN
              WRITE(6,*) 'After ',SUBR,': I,J,L,TSI=',I,J,L,TICE,HSI(:,I
     *             ,J),MSI(I,J),SNOWI(I,J),RSI(I,J)
c              QCHECKI = .TRUE.
            END IF
          END DO
          IF (SNOWI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,SNOWI=',I,J,SNOWI(I,J)
            QCHECKI = .TRUE.
          END IF
        END DO
      END DO
      IF (QCHECKI) STOP "CHECKI: Ice variables out of bounds"
          
      END SUBROUTINE CHECKI
