      MODULE SEAICE
!@sum  SEAICE contains all the sea ice related subroutines
!@auth Original Development Team
!@ver  1.0
!@cont PREC_SI,SEA_ICE.

      USE CONSTANT, only : lhm,rhoi,rhow,shi,shw
      IMPLICIT NONE

!@var XSI1,XSI2,XSI3,XSI4 fractions of mass layer in each temp. layer
      REAL*8, PARAMETER :: XSI1=0.5d0, XSI2=0.5d0,
     *                     XSI3=0.5d0, XSI4=0.5d0
!@var BYXSIn recipricol of XSIn
      REAL*8, PARAMETER :: BYXSI1=1./XSI1, BYXSI2=1./XSI2,
     *                     BYXSI3=1./XSI3, BYXSI4=1./XSI4
!@var Z1I thickness of first layer ice (m)
      REAL*8, PARAMETER :: Z1I = .1d0
!@var ACE1I ice mass first layer (kg/m^2)
      REAL*8, PARAMETER :: ACE1I = Z1I*RHOI
!@var HC1I heat capacity of first layer ice (J/m^2)
      REAL*8, PARAMETER :: HC1I = ACE1I*SHI
!@var Z2OIM thickness of 2nd layer ice (m)
      REAL*8, PARAMETER :: Z2OIM = .4d0
!@var AC2OIM ice mass 2nd layer (kg/m^2)
      REAL*8, PARAMETER :: AC2OIM = Z2OIM*RHOI
!@var ALAMI,ALAMS lambda coefficient for ice/snow J/(m*degC*sec)
      REAL*8, PARAMETER :: ALAMI=2.1762d0, ALAMS=0.35d0
!@var RHOS density of snow (kg/m^3)
      REAL*8, PARAMETER :: RHOS = 300.0
!@var BYRLI,BYRLS reciprical of density*lambda
      REAL*8, PARAMETER :: BYRLI = 1./(RHOI*ALAMI),
     *     BYRLS = 1./(RHOS*ALAMS)
!@var LMI number of temperature layers in ice
      INTEGER, PARAMETER :: LMI = 4

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,MSI1,TG1,TG2,TG3,TG4,PRCP,TPRCP
     *     ,EPRCP,RUN0,DIFS,EDIFS,ERUN2,QFIXR)
!@sum  PREC_SI Adds the precipitation to sea/lake ice
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER :: SNOMAX=100.0, dSNdRN=0.
      REAL*8 SNOW, MSI1, MSI2, MELT1
      REAL*8 TPRCP, EPRE, EPRCP, PRCP, RAIN, FREZ1
      REAL*8 TG1, TG2, TG3, TG4, HSI1, HSI2, HSI3, HSI4,
     *       FMSI1, FMSI2, FHSI1, FHSI2, FHSI3, H2, CMPRS
      REAL*8 DIFS, EDIFS, ERUN2 ! for diagnostics
      REAL*8 RUN0 ! runoff for ocean/lake
      REAL*8 HC_1
!@var QFIXR true if RSI and MSI2 are fixed (ie. for a fixed SST run)
      LOGICAL QFIXR

      RUN0=0. ; DIFS=0. ; EDIFS=0. ; ERUN2=0.
C**** sensible -  latent heat of precipit.
      EPRE = EPRCP - PRCP*LHM ! total energy of precipitation;
C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2
      HSI4 = (SHI*TG4-LHM)*XSI4*MSI2
      HC_1 = XSI1*MSI1*SHI
      HSI1 = HSI1+EPRE
      IF (TPRCP.LT.0.) GO TO 180
      IF (EPRCP.LT.-TG1*HC_1) GO TO 160
C**** ALL PRECIPITATION IS RAIN ABOVE 0degC
C**** RAIN COMPRESSES SNOW INTO ICE
      RAIN = PRCP
      IF (HSI1/LHM+XSI1*MSI1 .LE. 0.) GO TO 140
C**** WARM RAIN MELTS SOME SNOW OR ICE
      MELT1 = HSI1/LHM+XSI1*MSI1 ! melted snow and ice (kg/m^2)
      RUN0=MELT1 + PRCP ! water mass flux to ocean (kg/m^2)
      IF (MSI1-ACE1I .LE. MELT1) GO TO 130
C**** RAIN MELTS SOME SNOW AND COMPRESSES SNOW INTO ICE
      FMSI2 = MIN (dSNdRN*(RAIN+MELT1), MSI1-ACE1I-MELT1) ! > 0.
      FHSI1 = -LHM*(XSI1*FMSI2-XSI2*MELT1) ! downward heat flux
c     F1 = HSI1*(XSI1*FMSI2-XSI2*MELT1)/(XSI1*MSI1-MELT1)
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! > 0.
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = MSI1-MELT1-FMSI2
      SNOW = MSI1-ACE1I
      IF (SNOW .LT. 0.) SNOW = 0.
      HSI1 = HSI1-FHSI1
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      DIFS = FMSI2
      EDIFS = FHSI2
      ERUN2 = HSI3*FMSI2*BYXSI3/MSI2
      IF (QFIXR) RETURN
      GO TO 210
  130 CONTINUE
C**** RAIN MELTS ALL SNOW AND SOME ICE
      FMSI2 = MSI1-ACE1I-MELT1 ! < 0.(upward ice mass flux)
      DIFS = FMSI2 ! < 0.
C     FMSI1 = XSI1*(MSI1-ACE1I)-MELT1 ! < 0.(melted snow/ice mass)
      FHSI1 = HSI2*((XSI1/XSI2)*FMSI2-MELT1)/MSI1 !  upward heat flux
      FHSI2 = HSI3*FMSI2*BYXSI3/MSI2 !  upward heat flux into layer 2
      FHSI3 = HSI4*FMSI2/MSI2 !  upward heat flux into layer 3
      SNOW=0. ! Rain melted all snow
      MSI1 = ACE1I ! Keep the first layer ice mass constant ACE1I
      HSI1 = HSI1-FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      EDIFS=FHSI2 ! for diagnostics
      ERUN2=EDIFS ! for diagnostics
      IF (QFIXR) RETURN
      GO TO 210
  140 CONTINUE
C**** RAIN COMPRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      CMPRS = MIN(dSNdRN*RAIN, MSI1-ACE1I) ! all snow or part of rain
      IF (-HSI1/LHM-XSI1*MSI1 .LT. RAIN) GO TO 150
C**** ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
C     FMSI1 = XSI1*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      RUN0 = 0.
      FHSI1 = HSI1*(XSI1*CMPRS+RAIN)/(XSI1*MSI1+RAIN) ! downward
C      F1 = HSI1*(XSI1*CMPRS+RAIN)/(XSI1*MSI1) ! for prescribed ice
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = MSI1-CMPRS ! first layer ice mass
      SNOW = MSI1-ACE1I ! snow mass
      IF (SNOW .LT. 0.) SNOW=0.
      HSI1 = HSI1-FHSI1
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      DIFS = FMSI2
      EDIFS = FHSI2
      ERUN2 = HSI3*FMSI2*BYXSI3/MSI2
      IF (QFIXR) RETURN
      GO TO 210
  150 CONTINUE
C**** JUST PART OF RAIN FREEZES IN LAYER 1
      FREZ1 = -HSI1/LHM-XSI1*MSI1 ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI1*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI1*CMPRS+FREZ1) ! downward heat flux from layer 1
c      F1 = HSI1*(XSI1*CMPRS+FREZ1)/(XSI1*MSI1) ! for prescribed ice
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = MSI1-CMPRS ! first layer ice mass
      SNOW = MSI1-ACE1I ! snow mass
      IF (SNOW .LT. 0.) SNOW=0.
      HSI1 = HSI1-FHSI1
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      DIFS = FMSI2
      EDIFS = FHSI2
      ERUN2 = HSI3*FMSI2*BYXSI3/MSI2
      IF (QFIXR) RETURN
      GO TO 210
  160 CONTINUE
C**** PRECIPITATION IS A MIXTURE OF RAIN AND SNOW AT 0 degC
C**** RAIN COMRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      SNOW = -EPRCP/LHM ! snow fall
      RAIN = PRCP-SNOW  ! rain fall
      CMPRS = MIN(dSNdRN*RAIN, MSI1+SNOW-ACE1I) ! compression
      IF (-HSI1/LHM-XSI1*MSI1-SNOW .LT. RAIN) GO TO 170
C**** ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
      RUN0 =0.
      FMSI1 = XSI1*CMPRS+XSI2*SNOW+RAIN ! downward ice mass flux
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSI1*FMSI1/(XSI1*MSI1+PRCP) ! downward heat flux
c      F1 = HSI1*FMSI1/(XSI1*MSI1) ! for prescribed ice/ocean
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = MSI1+SNOW-CMPRS ! first layer ice mass
      SNOW = MSI1-ACE1I ! snow mass
      IF (SNOW .LT. 0.) SNOW=0.
      HSI1 = HSI1-FHSI1
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      DIFS = FMSI2
      EDIFS = FHSI2
      ERUN2 = HSI3*FMSI2*BYXSI3/MSI2
      IF (QFIXR) RETURN
      GO TO 210
  170 CONTINUE
C**** NOT ALL RAIN FREEZES IN LAYER 1
      FREZ1 = -HSI1/LHM-XSI1*MSI1-SNOW ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI1*CMPRS+XSI2*SNOW+FREZ1 ! downward ice mass flux
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI1*CMPRS+XSI2*SNOW+FREZ1) ! downward
c      F1 = HSI1*(XSI1*CMPRS+XSI2*SNOW+FREZ1)/(XSI1*MSI1) ! prescribed
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = MSI1+SNOW-CMPRS ! first layer ice mass
      SNOW = MSI1-ACE1I ! snow mass
      IF (SNOW .LT. 0.) SNOW=0.
      HSI1 = HSI1-FHSI1
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      DIFS = FMSI2
      EDIFS = FHSI2
      ERUN2 = HSI3*FMSI2*BYXSI3/MSI2
      IF (QFIXR) RETURN
      GO TO 210
  180 CONTINUE
C**** ALL PRECIPITATION IS SNOW, SNOW AMOUNT INCREASES
      RUN0 = 0.
      IF (MSI1+PRCP .GT. SNOMAX+ACE1I) GO TO 190
C     FMSI1 = XSI2*PRCP ! > 0.(snow fall to layer 1)
      FHSI1 = HSI1*XSI2*PRCP/(XSI1*MSI1+PRCP) ! downward heat flux
      MSI1 = MSI1+PRCP ! first layer ice mass
      SNOW = MSI1-ACE1I ! snow mass
      IF (SNOW .LT. 0.) SNOW=0.
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (HSI2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      RETURN
  190 CONTINUE
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
      FMSI2 = MSI1+PRCP-(0.9*SNOMAX+ACE1I) ! > 0.(compressed snow)
C     FMSI1 = XSI2*PRCP+XSI1*FMSI2 ! > 0.(downward ice mass flux)
      FHSI1 = HSI1*(XSI2*PRCP+XSI1*FMSI2)/(XSI1*MSI1+PRCP) ! downward
c      F1 = HSI1*(XSI2*PRCP+XSI1*FMSI2)/(XSI1*MSI1) ! for prescribed
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = 0.9d0*SNOMAX+ACE1I ! first layer ice mass
      SNOW=.9d0*SNOMAX ! snow mass
      DIFS = FMSI2 ! compressed snow mass
      EDIFS = FHSI2 ! energy of compressed snow mass
      HSI1 = HSI1-FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      IF (QFIXR) THEN
        ERUN2 = HSI3*FMSI2*BYXSI3/MSI2
        H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
        TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
        RETURN
      END IF
  210 CONTINUE
C**** ADVECT ICE (usually downwards)
      HSI2 = HSI2+(FHSI1-FHSI2) ! for predicted ice/ocean
      HSI3 = HSI3+(FHSI2-FHSI3) ! for predicted ice/ocean
      HSI4 = HSI4+FHSI3 ! for predicted ice/ocean
      TG2 = (HSI2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      TG3 = (HSI3/(XSI3*MSI2)+LHM)/SHI ! third layer ice temperature
      TG4 = (HSI4/(XSI4*MSI2)+LHM)/SHI ! fourth layer ice temperature
      MSI2 = MSI2+FMSI2 ! second layer sea ice mass (kg/m^2)

      RETURN
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,TG1,TG2,TG3,TG4,MSI1,MSI2
     *     ,F0DT,F1DT,EVAP,HSI1,HSI2,HSI3,HSI4,TGW,RUN0
     *     ,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,QFIXR)
!@sum  SEA_ICE applies surface fluxes to ice covered areas
!@auth Gary Russell
!@ver  1.0

      IMPLICIT NONE

      REAL*8, PARAMETER :: ALPHA = 1.0, dSNdML =0.
      REAL*8,  INTENT(IN) :: DTSRCE
!@var F0DT heat flux on the ice top surface (W/m^2)
      REAL*8,  INTENT(IN) :: F0DT
!@var F1DT heat flux between the 1st and 2nd ice layers (W/m^2)
      REAL*8,  INTENT(IN) :: F1DT
!@var EVAP evaporation/dew on the top ice surface (kg/m^2)
      REAL*8,  INTENT(IN) :: EVAP
!@var QFIXR  true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL, INTENT(IN) :: QFIXR

      REAL*8 ROICE, SNOW, MSI1, MSI2, ACE2M, MELT1, MELT4, DEW, CMPRS
      REAL*8 TG1, TG2, TG3, TG4, HSI1, HSI2, HSI3, HSI4,
     *       FMSI1, FMSI2, FMSI3, FMSI4, FHSI1, FHSI2, FHSI3, FHSI4
!@var TGW mixed layer temp.(C)
      REAL*8 TGW
      REAL*8 HC1, HC2, HC3, HC4
      REAL*8 dF1dTI, dF2dTI, dF3dTI, dF4dTI, F1, F2, F3
      REAL*8, INTENT(OUT) :: EDIFSI, F2DT, RUN0, DIFSI, DIFS, EDIFS
C**** Initiallise output
      F2DT=0. ; RUN0=0.  ; DIFSI=0. ; EDIFSI=0.
      DIFS=0. ; EDIFS=0. ; ACE2M=0.

      IF (ROICE .EQ. 0.) RETURN
      FMSI2=0 ; FHSI2=0 ; MELT1=0. ; MELT4=0.
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1 ! J/m^2
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1 ! J/m^2
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2 ! J/m^2
      HSI4 = (SHI*TG4-LHM)*XSI4*MSI2 ! J/m^2
C****
C**** OCEAN ICE, CALCULATE TG1 AND
C****
      HC1 = SHI*XSI1*MSI1 ! heat capacity of ice layer 1 (J/(degC*m^2))
      HC2 = SHI*XSI2*MSI1 ! heat capacity of ice layer 2 (J/(degC*m^2))
      HC3 = SHI*XSI3*MSI2 ! heat capacity of ice layer 3 (J/(degC*m^2))
      HC4 = SHI*XSI4*MSI2 ! heat capacity of ice layer 4 (J/(degC*m^2))
C**** CALCULATE AND APPLY DIFFUSIVE AND SURFACE ENERGY FLUXES
c      dF1dTI = 2.*DTSRCE/(ACE1I*BYRLI+SNOW*BYRLS) ! for non-Q-flux
      dF2dTI = ALAMI*RHOI*DTSRCE/(0.5*XSI2*MSI1+0.5*XSI3*MSI2)
C****          temperature derivative from F2 diffusive flux
      dF3dTI = ALAMI*RHOI*DTSRCE*2./MSI2
C****          temperature derivative from F3 diffusive flux
      dF4dTI = ALAMI*RHOI*DTSRCE*2.*BYXSI4/MSI2
C****          temperature derivative from F4 diffusive flux
CEXP  F2 = dF2dTI*(TG2-TG3) ! the diffusive
CEXP  F3 = dF3dTI*(TG3-TG4) ! fluxes from
CEXP  F4 = dF4dTI*(TG4-TGW) ! explicit method
C**** DIFFUSIVE FLUXES FROM IMPLICIT METHOD
c      F1 = dF1dTI*(HC1*(TG1-TG2)+ALPHA*F0DT)/
c     A     (HC1+ALPHA*dF1dTI)   ! for a non-Q-flux ocean model
c      F2 = dF2dTI*(HC2*(TG2-TG3)+ALPHA*F1)/
c     A     (HC2+ALPHA*dF2dTI) ! for a non-Q-flux ocean model

      F2 = dF2dTI*(HC2*(TG2-TG3)+ALPHA*F1DT)/
     A     (HC2+ALPHA*dF2dTI)
      F3 = dF3dTI*(HC3*(TG3-TG4)+ALPHA*F2)/
     A     (HC3+ALPHA*dF3dTI)
      F2DT = dF4dTI*(HC4*(TG4-TGW)+ALPHA*F3)/
     A            (HC4+ALPHA*dF4dTI)
      HSI1 = HSI1+(F0DT-F1DT)
      HSI2 = HSI2+(F1DT-F2)
      HSI3 = HSI3+(F2-F3)
      HSI4 = HSI4+(F3-F2DT)
      DEW = -EVAP ! dew to the surface
      IF (HSI1/LHM+XSI1*MSI1+DEW .LE. 0.) GO TO 160 ! go to freezing
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
      MELT1 = HSI1/LHM+XSI1*MSI1+DEW ! melting + dew to layer 1
C**** EVAPORATION FROM THE SURFACE
      IF (DEW .GT. 0.) GO TO 140 ! go to the dew case
C**** DEW IS EVAPORATION NOW
C**** EVAPORATION REDUCES SNOW OR ICE, MELT1>0.
      IF (MSI1-ACE1I+DEW .GT. MELT1) GO TO 130
C**** ALL SNOW AND SOME ICE EVAPORATE AND MELT
C**** ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 2 into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 3 into layer 2
      FHSI1 = HSI2*FMSI1*BYXSI2/MSI1 ! energy of ice mass FMSI1
      FHSI2 = HSI3*FMSI2*BYXSI3/MSI2 ! energy of ice mass FMSI2
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      SNOW = 0. ! all snow evaporates and melts
      IF (.not.QFIXR) GO TO 170
      GO TO 370
  130 CONTINUE
C**** EVAPORATION AND MELTING REDUCE JUST SNOW AMOUNT
C**** ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
      FMSI2 = MIN(dSNdML*MELT1, MSI1-ACE1I+DEW-MELT1) ! > 0.
      FMSI1 = XSI1*FMSI2+XSI2*(DEW-MELT1)
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*BYXSI2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward from layer 2 -> layer 3
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      MSI1 = MSI1+DEW-MELT1-FMSI2 ! kg/m^2
      SNOW = MSI1-ACE1I
      IF (.not.QFIXR) GO TO 210
      GO TO 370
  140 CONTINUE
C**** DEW INCREASES ICE AMOUNT,  MELT1 > 0.
      IF (MSI1-ACE1I .GT. MELT1) GO TO 150
C**** ALL SNOW AND SOME ICE MELT
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW-MELT1 ! upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1
      SNOW = 0.
      IF (FMSI2 .LE. 0.) THEN ! (if melting is greater than dew)
C****ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
        FHSI1 = HSI2*FMSI1*BYXSI2/MSI1 ! upward -> layer 1 from layer 2
        FHSI2 = HSI3*FMSI2*BYXSI3/MSI2 ! upward -> layer 2 from layer 3
        HSI1 = HSI1-FHSI1
        HSI2 = HSI2+(FHSI1-FHSI2)
        IF (.not.QFIXR) GO TO 170
        GO TO 370
      ENDIF
C**** ICE ADVECTION IS DOWNWARD INTO LAYER 3 FROM LAYER 2
C**** (if dew is greater than melting)
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*BYXSI2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward -> layer 3 from layer 2
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      MSI1 = ACE1I
      IF (.not.QFIXR) GO TO 210
      GO TO 370
  150 CONTINUE
C**** MELTING REDUCES SNOW AMOUNT
C**** ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
      CMPRS = MIN(dSNdML*MELT1, MSI1-ACE1I-MELT1) ! > 0.
      FMSI1 = DEW+XSI1*CMPRS-XSI2*MELT1
      FMSI2 = DEW+CMPRS ! > 0. downward into layer 3 from layer 2
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*BYXSI2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 3
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward energy flux into layer 3
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      MSI1 = MSI1-MELT1-CMPRS ! new mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass after melting and compression
      IF (.not.QFIXR) GO TO 210
      GO TO 370
  160 CONTINUE
C**** THE FIRST LAYER IS BELOW FREEZING
C**** NO SNOW OR ICE MELTS IN LAYER 1
      IF (DEW .GT. 0.) GO TO 200 ! go to the dew case
      IF (MSI1-ACE1I+DEW .GE. 0.) GO TO 190
C**** ALL SNOW AND SOME ICE EVAPORATE
C**** ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW ! < 0. upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW ! < 0. upward into layer 2 from layer 3
      FHSI1 = HSI2*FMSI1*BYXSI2/MSI1 ! upward energy flux into layer 1
      FHSI2 = HSI3*FMSI2*BYXSI3/MSI2 ! upward energy flux into layer 2
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      SNOW = 0. ! all snow evaporated
      IF (QFIXR) GO TO 370
  170 CONTINUE
      MSI1 = ACE1I
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 180 ! go to freezing case
C**** FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
      MELT4 = HSI4/LHM+XSI4*MSI2 ! > 0. melted ice from layer 4
      FMSI3 = XSI4*FMSI2+XSI3*MELT4
      IF (FMSI3 .LE. 0.) FHSI3 = -LHM*FMSI3 ! upward into layer 3
      IF (FMSI3 .GT. 0.) FHSI3 = HSI3*FMSI3*BYXSI3/MSI2 ! downward
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
      MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      GO TO 230
  180 CONTINUE
C**** NO ICE MELTS IN LAYER 4
C     FMSI3 = XSI4*FMSI2 ! upward mass flux into layer 3 from layer 4
      FHSI3 = HSI4*FMSI2/MSI2 ! upward energy flux into layer 3
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
      MSI2 = MSI2+FMSI2
      GO TO 230
  190 CONTINUE
C**** JUST SOME SNOW EVAPORATES
C**** NO ADVECTION BETWEEN LAYER 2 AND 3
C     FMSI1 = XSI2*DEW ! < 0. upward mass flux into layer 1
C     FMSI2 = 0. ! no ice advection between layers 2 and 3
      FHSI1 = HSI2*DEW/MSI1 ! upward energy flux into layer 1
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+FHSI1
      MSI1 = MSI1+DEW ! new ice mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass
      IF (QFIXR) GO TO 370
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 230 ! go to freezing case
C**** FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
      MELT4 = HSI4/LHM+XSI4*MSI2 ! melted ice from layer 4
C     FMSI3 = XSI3*MELT4 ! > 0. downward mass flux into layer 4
      FHSI3 = HSI3*MELT4/MSI2 ! downward heat flux into layer 4
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+FHSI3
      MSI2 = MSI2-MELT4 ! new ice mass of the second physical layer
      GO TO 230
  200 CONTINUE
C**** DEW INCREASES ICE AMOUNT, ADVECT ICE DOWNWARD
C**** ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
C     FMSI1 = DEW ! > 0. downward mass flux into layer 2
      FMSI2 = DEW ! > 0. downward mass flux into layer 3
      FHSI1 = HSI1*FMSI2/(XSI1*MSI1+DEW) ! downward heat flux
      FHSI2 = HSI2*FMSI2*BYXSI2/MSI1 ! downward heat flux into layer 3
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      IF (QFIXR) GO TO 370
  210 CONTINUE
      IF (HSI4/LHM+XSI4*MSI2 .GT. 0.) THEN
C**** FLUXES HEAT UP TG2 TO FREEZING POINT AND MELT SOME ICE
        MELT4 = HSI4/LHM+XSI4*MSI2 ! melted ice from layer 4
C     FMSI3 = XSI4*FMSI2+XSI3*MELT4 ! > 0. downward into layer 4
        FHSI3 = HSI3*((XSI4/XSI3)*FMSI2+MELT4)/MSI2 ! downward
        HSI3 = HSI3+(FHSI2-FHSI3)
        HSI4 = HSI4+FHSI3
        MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      ELSE
C**** NO ICE MELTS IN LAYER 4
C     FMSI3 = XSI4*FMSI2 ! > 0. downward mass flux into layer 4
        FHSI3 = HSI3*(XSI4/XSI3)*FMSI2/MSI2 ! downward heat flux
        HSI3 = HSI3+(FHSI2-FHSI3)
        HSI4 = HSI4+FHSI3
        MSI2 = MSI2+FMSI2
      END IF
 230  CONTINUE
 370  CONTINUE
C**** save output diagnostics
      DIFS  = FMSI2 ! <0 upward ice mass into layer 2
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      DIFSI = ROICE*DIFS
      EDIFSI= ROICE*EDIFS
      ACE2M = MELT4             ! melted ice at the bottom
      RUN0 = MELT1 ! water mass that flows to the ocean (kg/m^2)
      RETURN
      END SUBROUTINE SEA_ICE

      SUBROUTINE ADDICE (SNOW,ROICE,TG1,TG2,TG3,TG4,MSI1,MSI2,HSI1,HSI2
     *     ,HSI3,HSI4,DIFSI,EDIFSI,ENRGFO,ACEFO,ACE2F,ENRGFI,TFW,FLEAD
     *     ,QFIXR,QCMPR)
!@sum  ADDICE adds ice formed in the ocean to ice variables
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER :: YSI1 = XSI1*ACE1I/(ACE1I+AC2OIM),
     *                     YSI2 = XSI2*ACE1I/(ACE1I+AC2OIM),
     *                     YSI3 = XSI3*AC2OIM/(ACE1I+AC2OIM),
     *                     YSI4 = XSI4*AC2OIM/(ACE1I+AC2OIM)
      REAL*8, PARAMETER :: Z2OIX = 4.9,
     *                     BYZICX=1./(Z1I+Z2OIX)
!@var QFIXR  true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL, INTENT(IN) :: QFIXR
!@var QCMPR  true if ice should be compressed due to leads etc.
      LOGICAL, INTENT(IN) :: QCMPR
!@var TFW freezing temperature for water underlying ice (C)
      REAL*8, INTENT(IN) :: TFW
!@var FLEAD minimum lead fraction for ice (%)
      REAL*8, INTENT(IN) :: FLEAD

      REAL*8 ROICE, ROICEN, OPNOCN, SNOW, MSI1, MSI2, DRSI, DRI
      REAL*8, INTENT(OUT) :: TG1, TG2, TG3, TG4
      REAL*8, INTENT(INOUT) :: HSI1, HSI2, HSI3, HSI4
      REAl*8 FMSI1, FMSI2, FMSI3, FMSI4, FHSI1, FHSI2, FHSI3, FHSI4
      REAL*8, INTENT(IN) ::  ENRGFI, ENRGFO, ACEFO, ACE2F
      REAL*8, INTENT(INOUT) :: EDIFSI, DIFSI
      REAL*8 DIFS1

      IF (.not.QFIXR .and. ROICE.LE.0. .and. ACEFO.gt.0) THEN
        ROICE=ACEFO/(ACE1I+AC2OIM)
        SNOW=0.
        TG1=TFW
        TG2=TFW
        TG3=TFW
        TG4=TFW
        MSI2=AC2OIM
        RETURN
      END IF
C****
      IF (ROICE.le.0) RETURN

      IF (QFIXR) GO TO 370
      IF (ACE2F.le.0) GO TO 250 ! go to no freezing case
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
C     FMSI3 = -XSI3*ACE2F ! < 0.
C     FMSI4 = -ACE2F
      FHSI3 = -HSI4*ACE2F*(XSI3/XSI4)/MSI2
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      IF (ACEFO .GT. 0.) GO TO 240
C**** NEW ICE IS FORMED BELOW OLD SEA ICE
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+(FHSI3+ENRGFI)
      MSI2 = MSI2+ACE2F ! new ice mass of physical layer 2
      GO TO 270
  240 CONTINUE
C**** NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACE2F))/(ROICE+DRSI) ! layer 2
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      HSI1 = ((1.-ROICE)*ENRGFO*YSI1+ROICE*HSI1)/(ROICE+DRSI)
      HSI2 = ((1.-ROICE)*ENRGFO*YSI2+ROICE*HSI2)/(ROICE+DRSI)
      HSI3 = ((1.-ROICE)*ENRGFO*YSI3+ROICE*(HSI3-FHSI3))/
     A       (ROICE+DRSI)
      HSI4 = ((1.-ROICE)*ENRGFO*YSI4+ROICE*(HSI4+FHSI3+ENRGFI))/
     A       (ROICE+DRSI)
      ROICE = ROICE+DRSI ! new ice concentration
      GO TO 270
  250 CONTINUE
      IF (ACEFO .GT. 0.) GO TO 260 ! new ice on the open ocean
C**** NO NEW ICE IS FORMED UNDERNEATH THE OLD ONE
      GO TO 270
  260 CONTINUE
C**** NEW ICE IS FORMED ON THE OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*MSI2)/(ROICE+DRSI) ! layer 2
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      HSI1 = ((1.-ROICE)*ENRGFO*YSI1+ROICE*HSI1)/(ROICE+DRSI)
      HSI2 = ((1.-ROICE)*ENRGFO*YSI2+ROICE*HSI2)/(ROICE+DRSI)
      HSI3 = ((1.-ROICE)*ENRGFO*YSI3+ROICE*HSI3)/(ROICE+DRSI)
      HSI4 = ((1.-ROICE)*ENRGFO*YSI4+ROICE*HSI4)/(ROICE+DRSI)
      ROICE = ROICE+DRSI ! new ice concentration
  270 CONTINUE
      IF (QCMPR) THEN           ! COMPRESS THE ICE HORIZONTALLY
      IF (MSI2 .GE. AC2OIM) GO TO 280 ! ice is thick enough
C**** SEA ICE IS TOO THIN
      ROICEN = ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM) ! new ice concentr.
      GO TO 290
  280 CONTINUE
      OPNOCN=FLEAD*(RHOI/(ACE1I+MSI2)-BYZICX)
      IF ((1.-ROICE) .GE. OPNOCN) GO TO 360
      ROICEN = 1.-OPNOCN
 290  CONTINUE
      DRSI = ROICEN-ROICE ! < 0. compressed ice concentration
      DRI = ROICE-ROICEN ! > 0. for diagnostics
C     FMSI3 = XSI3*FMSI4 ! < 0. upward ice mass flux into layer 3
      FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! upward ice mass into layer 4
      FHSI3 = HSI4*FMSI4*(XSI3/XSI4)/MSI2 ! upward heat flux
      FHSI4 = (HSI1+HSI2+HSI3+HSI4)*(DRSI/ROICEN)
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+(FHSI3-FHSI4)
      MSI2 = MSI2-FMSI4 ! new ice mass of second physical layer
C     SNOW = SNOW   ! snow thickness is conserved
      DIFS1 = DRI*ACE1I/ROICE
      EDIFSI = EDIFSI+ROICE*(HSI1+HSI2)*(DIFS1/MSI1)*0.5
      DIFSI = DIFSI+ROICE*DIFS1
      ROICE = ROICEN
      END IF
C**** RESAVE PROGNOSTIC QUANTITIES
  360 CONTINUE
  370 CONTINUE
      TG1 = (HSI1/(XSI1*MSI1) +LHM)/SHI ! temperature of layer 1
      TG2 = (HSI2/(XSI2*MSI1) +LHM)/SHI ! temperature of layer 2
      TG3 = (HSI3/(XSI3*MSI2) +LHM)/SHI ! temperature of layer 3
      TG4 = (HSI4/(XSI4*MSI2) +LHM)/SHI ! temperature of layer 4

      RETURN
      END SUBROUTINE ADDICE

      SUBROUTINE SIMELT(ROICE,SNOW,MSI1,MSI2,ACE,TG1,TG2,TG3,TG4
     *     ,ENRGW,PWATER,ENRGUSED)
!@sum  SIMELT melts sea ice if surrounding water is warm
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
!@var ENRGUSED energy used to melt ice (J/m^2)
      REAL*8, INTENT(OUT) :: ENRGUSED
!@var ENRGW energy available to melt ice (J/m^2)
      REAL*8, INTENT(IN) :: ENRGW
      REAL*8 MSI1, MSI2, MELT,DRSI,ROICEN,FHSI4,FHSI3,ACE
      REAL*8 E_SIDE,E_BOTTOM,GAMMA,H_C,H_ICE,PWATER,POICE
     *     ,ENRGI,HSI4,HSI3,HSI2,HSI1,TG1,TG2,TG3,TG4
     *     ,SNOW,ROICE

C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2
      HSI4 = (SHI*TG4-LHM)*XSI4*MSI2
      ENRGI = HSI1+HSI2+HSI3+HSI4 ! energy in sea ice [J/m^2]
      IF (ROICE*ENRGI+ENRGW.LT.0.) GO TO 230
C**** THE WARM OCEAN MELTS ALL THE SNOW AND ICE
      ENRGUSED=-ROICE*ENRGI     ! only some energy is used
      ROICE=0.
      RETURN
C**** THE WARM OCEAN COOLS TO 0 DEGREES MELTING SOME SNOW AND ICE
C**** Reduce the ice depth
 230  ENRGUSED=ENRGW            ! all energy is used
      POICE = ROICE*PWATER      ! sea ice fraction
      H_C = 1.                  ! critical ice thickness [m]
      H_ICE = ACE/RHOI   ! sea ice thickness [m]
      IF (H_ICE .GT. 1.) THEN
        GAMMA = H_ICE/H_C
      ELSE
        GAMMA = 1.
      END IF
      E_BOTTOM = ENRGW*(POICE/(1.+GAMMA)) ! for bottom melting
      E_SIDE = ENRGW*(PWATER+(POICE*GAMMA)/(1.+GAMMA)) ! for side melt.
      E_SIDE = ENRGW-E_BOTTOM
C**** MELT ICE VERTICALLY, AND THEN HORIZONTALLY
      MELT = -XSI4*MSI2*ENRGW/HSI4 ! melted ice at the bottom
      IF (MSI2-MELT .LT. AC2OIM) MELT = MSI2-AC2OIM
C     FMSI3 = XSI3*MELT ! > 0.
      FHSI3 = HSI3*MELT/MSI2
      FHSI4 = HSI4*MELT*BYXSI4/MSI2
C**** MELT SOME ICE HORIZONTALLY
      DRSI = (ENRGW+ROICE*FHSI4)/(HSI1+HSI2+HSI3+HSI4-FHSI4) ! < 0.
      ROICEN = ROICE+DRSI       ! new sea ice concentration
c     SNOW = SNOW*(ROICE/ROICEN) ! new snow ice mass
      MSI2 = MSI2-MELT
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+FHSI3-FHSI4
C**** CONVERT SEA ICE ENTHALPY MINUS LATENT HEAT INTO TEMPERATURE
      TG1 = (HSI1/(XSI1*MSI1) +LHM)/SHI ! temperatute of layer 1
      TG2 = (HSI2/(XSI2*MSI1) +LHM)/SHI ! temperatute of layer 2
      TG3 = (HSI3/(XSI3*MSI2) +LHM)/SHI ! temperatute of layer 3
      TG4 = (HSI4/(XSI4*MSI2) +LHM)/SHI ! temperatute of layer 4
      ROICE=ROICEN
      RETURN
      END SUBROUTINE SIMELT

      END MODULE SEAICE

      MODULE SEAICE_COM
!@sum  SEAICE_COM contains the model arrays for seaice 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : im,jm
      USE SEAICE, only : lmi

      IMPLICIT NONE
!@var RSI fraction of open water area covered in ice
      REAL*8, DIMENSION(IM,JM) :: RSI
!@var SNOWI snow amount on sea ice (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: SNOWI
!@var MSI thickness of ice second layer (layer 1=const) (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: MSI
!@var TSI temperature of the each ice layer (C)
      REAL*8, DIMENSION(IM,JM,LMI) :: TSI
!@var HSI enthaply of each ice layer (J/m^2) (will replace TSI)
c      REAL*8, DIMENSION(IM,JM,LMI) :: HSI

      END MODULE SEAICE_COM
