      MODULE SEAICE
!@sum  SEAICE contains all the sea ice related subroutines
!@auth Original Development Team
!@ver  1.0
!@cont PREC_SI,SEA_ICE.
      USE CONSTANT, only : lhm,rhoi,rhow,shi,shw,byshi
      IMPLICIT NONE

!@param LMI number of temperature layers in ice
      INTEGER, PARAMETER :: LMI = 4
!@param XSI fractions of mass layer in each temp. layer
!@param BYXSI recipricol of XSI
      REAL*8, PARAMETER, DIMENSION(LMI) ::
     *     XSI= (/0.5d0, 0.5d0, 0.5d0, 0.5d0/),
     *     BYXSI=1./XSI
!@param Z1I thickness of first layer ice (m)
      REAL*8, PARAMETER :: Z1I = .1d0
!@param ACE1I ice mass first layer (kg/m^2)
      REAL*8, PARAMETER :: ACE1I = Z1I*RHOI
!@param HC1I heat capacity of first layer ice (J/m^2)
      REAL*8, PARAMETER :: HC1I = ACE1I*SHI
!@param Z2OIM thickness of 2nd layer ice (m)
      REAL*8, PARAMETER :: Z2OIM = .4d0
!@param AC2OIM ice mass 2nd layer (kg/m^2)
      REAL*8, PARAMETER :: AC2OIM = Z2OIM*RHOI
!@param ALAMI,ALAMS lambda coefficient for ice/snow J/(m*degC*sec)
      REAL*8, PARAMETER :: ALAMI=2.1762d0, ALAMS=0.35d0
!@param RHOS density of snow (kg/m^3)
      REAL*8, PARAMETER :: RHOS = 300.0
!@param BYRLI,BYRLS reciprical of density*lambda
      REAL*8, PARAMETER :: BYRLI = 1./(RHOI*ALAMI),
     *     BYRLS = 1./(RHOS*ALAMS)

      CONTAINS

      SUBROUTINE PREC_SI(SNOW,MSI2,TSIL,PRCP,ENRGP,RUN0,DIFS,EDIFS,ERUN2
     *     ,QFIXR)
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
!@var TSIL temperature of ice layers (C)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: TSIL
!@var DIFS  upward ice mass into layer 2  (kg/m^2)
!@var EDIFS energy of ice mass DIFS (J/m^2)
!@var ERUN2 implied energy flux at base to keep fixed ice mass (J/m^2)
      REAL*8, INTENT(OUT) :: DIFS, EDIFS, ERUN2
!@var RUN0 runoff from ice (kg/m^2)
      REAL*8, INTENT(OUT) :: RUN0
      REAL*8, DIMENSION(LMI) :: HSI
      REAL*8 :: BYMSI2, FMSI1, FMSI2, FHSI1, FHSI2, FHSI3, CMPRS, SNWF,
     *     FHSI2U, MELT1, RAIN, FREZ1, MSI1

C**** initialize fluxes
      RUN0=0. ; FMSI2=0. ; FHSI1=0 ; FHSI2=0. ; FHSI3=0. ; FHSI2U=0.
C**** reciprocal of ice thickness for efficiency
      BYMSI2=1./MSI2
      MSI1=SNOW+ACE1I

C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI(1) = (SHI*TSIL(1)-LHM)*XSI(1)*MSI1
      HSI(2) = (SHI*TSIL(2)-LHM)*XSI(2)*MSI1
      HSI(3) = (SHI*TSIL(3)-LHM)*XSI(3)*MSI2
      HSI(4) = (SHI*TSIL(4)-LHM)*XSI(4)*MSI2
      HSI(1) = HSI(1)+ENRGP  ! add total energy of precipitation
      IF (ENRGP.LE. -PRCP*LHM) GO TO 180
      IF (ENRGP.LE. 0.) GO TO 160
C**** ALL PRECIPITATION IS RAIN ABOVE 0degC
C**** RAIN COMPRESSES SNOW INTO ICE
      RAIN = PRCP
      IF (HSI(1)/LHM+XSI(1)*MSI1 .LE. 0.) GO TO 140
C**** WARM RAIN MELTS SOME SNOW OR ICE
      MELT1 = HSI(1)/LHM+XSI(1)*MSI1 ! melted snow and ice (kg/m^2)
      RUN0=MELT1 + PRCP ! water mass flux to ocean (kg/m^2)
      IF (MSI1-ACE1I .LE. MELT1) GO TO 130
C**** RAIN MELTS SOME SNOW AND COMPRESSES SNOW INTO ICE
      FMSI2 = MIN (dSNdRN*(RAIN+MELT1), MSI1-ACE1I-MELT1) ! > 0.
      FHSI1 = -LHM*(XSI(1)*FMSI2-XSI(2)*MELT1) ! downward heat flux
c     F1 = HSI(1)*(XSI(1)*FMSI2-XSI(2)*MELT1)/(XSI(1)*MSI1-MELT1)
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! > 0.
      FHSI3 = HSI(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1-MELT1-FMSI2
      GO TO 210
  130 CONTINUE
C**** RAIN MELTS ALL SNOW AND SOME ICE
      FMSI2 = MSI1-ACE1I-MELT1 ! < 0.(upward ice mass flux)
C     FMSI1 = XSI(1)*(MSI1-ACE1I)-MELT1 ! < 0.(melted snow/ice mass)
      FHSI1 = HSI(2)*((XSI(1)/XSI(2))*FMSI2-MELT1)/MSI1 !  upward heat flux
      FHSI2U= HSI(3)*FMSI2*BYXSI(3)*BYMSI2 !  upward heat flux into layer 2
      FHSI3 = HSI(4)*FMSI2*BYMSI2 !  upward heat flux into layer 3
      MSI1 = ACE1I ! Keep the first layer ice mass constant ACE1I
      GO TO 210
  140 CONTINUE
C**** RAIN COMPRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      CMPRS = MIN(dSNdRN*RAIN, MSI1-ACE1I) ! all snow or part of rain
      IF (-HSI(1)/LHM-XSI(1)*MSI1 .LT. RAIN) GO TO 150
C**** ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
C     FMSI1 = XSI(1)*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSI(1)*(XSI(1)*CMPRS+RAIN)/(XSI(1)*MSI1+RAIN) ! downward
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1-CMPRS ! first layer ice mass
      GO TO 210
  150 CONTINUE
C**** JUST PART OF RAIN FREEZES IN LAYER 1
      FREZ1 = -HSI(1)/LHM-XSI(1)*MSI1 ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI(1)*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI(1)*CMPRS+FREZ1) ! downward heat flux from layer 1
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1-CMPRS ! first layer ice mass
      GO TO 210
  160 CONTINUE
C**** PRECIPITATION IS A MIXTURE OF RAIN AND SNOW AT 0 degC
C**** RAIN COMRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      SNWF = -ENRGP/LHM ! snow fall
      RAIN = PRCP-SNWF  ! rain fall
      CMPRS = MIN(dSNdRN*RAIN, MSI1+SNWF-ACE1I) ! compression
      IF (-HSI(1)/LHM-XSI(1)*MSI1-SNWF .LT. RAIN) GO TO 170
C**** ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
      FMSI1 = XSI(1)*CMPRS+XSI(2)*SNWF+RAIN ! downward ice mass flux
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSI(1)*FMSI1/(XSI(1)*MSI1+PRCP) ! downward heat flux
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1+SNWF-CMPRS ! first layer ice mass
      GO TO 210
  170 CONTINUE
C**** NOT ALL RAIN FREEZES IN LAYER 1
      FREZ1 = -HSI(1)/LHM-XSI(1)*MSI1-SNWF ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI(1)*CMPRS+XSI(2)*SNWF+FREZ1 ! downward ice mass flux
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI(1)*CMPRS+XSI(2)*SNWF+FREZ1) ! downward
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = MSI1+SNWF-CMPRS ! first layer ice mass
      GO TO 210
  180 CONTINUE
C**** ALL PRECIPITATION IS SNOW, SNOW AMOUNT INCREASES
      IF (MSI1+PRCP .LE. SNOMAX+ACE1I) THEN
C     FMSI1 = XSI(2)*PRCP ! > 0.(snow fall to layer 1)
      FHSI1 = HSI(1)*XSI(2)*PRCP/(XSI(1)*MSI1+PRCP) ! downward heat flux
      MSI1 = MSI1+PRCP ! first layer ice mass
      ELSE
C**** TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
      FMSI2 = MSI1+PRCP-(0.9*SNOMAX+ACE1I) ! > 0.(compressed snow)
C     FMSI1 = XSI(2)*PRCP+XSI(1)*FMSI2 ! > 0.(downward ice mass flux)
      FHSI1 = HSI(1)*(XSI(2)*PRCP+XSI(1)*FMSI2)/(XSI(1)*MSI1+PRCP) ! downward
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI(3)*FMSI2*(XSI(4)/XSI(3))*BYMSI2 ! downward heat flux
      MSI1 = 0.9d0*SNOMAX+ACE1I ! first layer ice mass
      END IF
  210 CONTINUE
C**** Save fluxes for diagnostic output
      DIFS  = FMSI2
      EDIFS = FHSI2+FHSI2U
      ERUN2 = HSI(3)*FMSI2*BYXSI(3)*BYMSI2

C**** Adjust amounts resulting from fluxes.
C**** Note that the flux into layer two is seperated into up and down
C**** components so that the code is valid for fixed sea ice runs also
      SNOW = MSI1-ACE1I         ! snow mass
      IF (SNOW .LT. 0.) SNOW=0. ! to catch roundoff errors
      HSI(1) = HSI(1)- FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)

      IF (.not. QFIXR) THEN     ! ADVECT ICE for predicted ice
        MSI2 = MSI2+FMSI2       ! 2nd layer sea ice mass (kg/m^2)
        HSI(2) = HSI(2)- FHSI2U     ! flux associated with upward advection
        HSI(3) = HSI(3)+(FHSI2U+FHSI2-FHSI3)
        HSI(4) = HSI(4)+ FHSI3
        TSIL(3) = (HSI(3)/(XSI(3)*MSI2)+LHM)*BYSHI ! 3rd layer ice temperature
        TSIL(4) = (HSI(4)/(XSI(4)*MSI2)+LHM)*BYSHI ! 4th layer ice temperature
      END IF
      TSIL(1) = (HSI(1)/(XSI(1)*MSI1)+LHM)*BYSHI ! 1st layer ice temperature
      TSIL(2) = (HSI(2)/(XSI(2)*MSI1)+LHM)*BYSHI ! 2nd layer ice temperature

      RETURN
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,TSIL,MSI2,F0DT,F1DT,EVAP,HSI
     *     ,TGW,RUN0,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,QFIXR)
!@sum  SEA_ICE applies surface fluxes to ice covered areas
!@auth Gary Russell
!@ver  1.0

      IMPLICIT NONE

      REAL*8, PARAMETER :: ALPHA = 1.0, dSNdML =0.
      REAL*8, INTENT(IN) :: DTSRCE
!@var F0DT heat flux on the ice top surface (W/m^2)
      REAL*8, INTENT(IN) :: F0DT
!@var F1DT heat flux between the 1st and 2nd ice layers (W/m^2)
      REAL*8, INTENT(IN) :: F1DT
!@var EVAP evaporation/dew on the top ice surface (kg/m^2)
      REAL*8, INTENT(IN) :: EVAP
!@var QFIXR  true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL,INTENT(IN) :: QFIXR
!@var TGW temperature of water below ice (C)
      REAL*8, INTENT(IN) :: TGW

      REAL*8, INTENT(IN) :: ROICE
      REAL*8, INTENT(IN), DIMENSION(LMI) :: TSIL
      REAL*8, INTENT(INOUT) :: SNOW, MSI2
      REAL*8, INTENT(OUT) :: ACE2M
      REAL*8, INTENT(OUT),DIMENSION(LMI) :: HSI
      REAL*8, INTENT(OUT) :: EDIFSI, F2DT, RUN0, DIFSI, DIFS, EDIFS
      REAL*8 MSI1, MELT1, MELT4, DEW, CMPRS
      REAL*8 FMSI1, FMSI2, FMSI3, FMSI4, FHSI1, FHSI2, FHSI3, FHSI4
      REAL*8 HC1, HC2, HC3, HC4
      REAL*8 dF1dTI, dF2dTI, dF3dTI, dF4dTI, F1, F2, F3
C**** Initiallise output
      F2DT=0. ; RUN0=0.  ; DIFSI=0. ; EDIFSI=0.
      DIFS=0. ; EDIFS=0. ; ACE2M=0.
      HSI(1)=0 ; HSI(2)=0 ; HSI(3)=0 ; HSI(4)=0

      IF (ROICE .EQ. 0.) RETURN
      FMSI2=0 ; FHSI2=0 ; MELT1=0. ; MELT4=0.
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI(1) = (SHI*TSIL(1)-LHM)*XSI(1)*MSI1 ! J/m^2
      HSI(2) = (SHI*TSIL(2)-LHM)*XSI(2)*MSI1 ! J/m^2
      HSI(3) = (SHI*TSIL(3)-LHM)*XSI(3)*MSI2 ! J/m^2
      HSI(4) = (SHI*TSIL(4)-LHM)*XSI(4)*MSI2 ! J/m^2
C****
C**** OCEAN ICE, CALCULATE TSIL(1) AND
C****
      HC1 = SHI*XSI(1)*MSI1 ! heat capacity of ice layer 1 (J/(degC*m^2))
      HC2 = SHI*XSI(2)*MSI1 ! heat capacity of ice layer 2 (J/(degC*m^2))
      HC3 = SHI*XSI(3)*MSI2 ! heat capacity of ice layer 3 (J/(degC*m^2))
      HC4 = SHI*XSI(4)*MSI2 ! heat capacity of ice layer 4 (J/(degC*m^2))
C**** CALCULATE AND APPLY DIFFUSIVE AND SURFACE ENERGY FLUXES
c      dF1dTI = 2.*DTSRCE/(ACE1I*BYRLI+SNOW*BYRLS) ! for non-Q-flux
      dF2dTI = ALAMI*RHOI*DTSRCE/(0.5*XSI(2)*MSI1+0.5*XSI(3)*MSI2)
C****          temperature derivative from F2 diffusive flux
      dF3dTI = ALAMI*RHOI*DTSRCE*2./MSI2
C****          temperature derivative from F3 diffusive flux
      dF4dTI = ALAMI*RHOI*DTSRCE*2.*BYXSI(4)/MSI2
C****          temperature derivative from F4 diffusive flux
CEXP  F2 = dF2dTI*(TSIL(2)-TSIL(3)) ! the diffusive
CEXP  F3 = dF3dTI*(TSIL(3)-TSIL(4)) ! fluxes from
CEXP  F4 = dF4dTI*(TSIL(4)-TGW) ! explicit method
C**** DIFFUSIVE FLUXES FROM IMPLICIT METHOD
c      F1 = dF1dTI*(HC1*(TSIL(1)-TSIL(2))+ALPHA*F0DT)/
c     A     (HC1+ALPHA*dF1dTI)   ! for a non-Q-flux ocean model
c      F2 = dF2dTI*(HC2*(TSIL(2)-TSIL(3))+ALPHA*F1)/
c     A     (HC2+ALPHA*dF2dTI) ! for a non-Q-flux ocean model

      F2 = dF2dTI*(HC2*(TSIL(2)-TSIL(3))+ALPHA*F1DT)/
     A     (HC2+ALPHA*dF2dTI)
      F3 = dF3dTI*(HC3*(TSIL(3)-TSIL(4))+ALPHA*F2)/
     A     (HC3+ALPHA*dF3dTI)
      F2DT = dF4dTI*(HC4*(TSIL(4)-TGW)+ALPHA*F3)/
     A            (HC4+ALPHA*dF4dTI)
      HSI(1) = HSI(1)+(F0DT-F1DT)
      HSI(2) = HSI(2)+(F1DT-F2)
      HSI(3) = HSI(3)+(F2-F3)
      HSI(4) = HSI(4)+(F3-F2DT)
      DEW = -EVAP ! dew to the surface
      IF (HSI(1)/LHM+XSI(1)*MSI1+DEW .LE. 0.) GO TO 160 ! go to freezing
C**** FLUXES HEAT UP TSIL(1) TO FREEZING POINT AND MELT SOME SNOW AND ICE
      MELT1 = HSI(1)/LHM+XSI(1)*MSI1+DEW ! melting + dew to layer 1
C**** EVAPORATION FROM THE SURFACE
      IF (DEW .GT. 0.) GO TO 140 ! go to the dew case
C**** DEW IS EVAPORATION NOW
C**** EVAPORATION REDUCES SNOW OR ICE, MELT1>0.
      IF (MSI1-ACE1I+DEW .GT. MELT1) GO TO 130
C**** ALL SNOW AND SOME ICE EVAPORATE AND MELT
C**** ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
      FMSI1 = XSI(1)*(MSI1-ACE1I)+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 2 into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 3 into layer 2
      FHSI1 = HSI(2)*FMSI1*BYXSI(2)/MSI1 ! energy of ice mass FMSI1
      FHSI2 = HSI(3)*FMSI2*BYXSI(3)/MSI2 ! energy of ice mass FMSI2
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)
      SNOW = 0. ! all snow evaporates and melts
      IF (.not.QFIXR) GO TO 170
      GO TO 370
  130 CONTINUE
C**** EVAPORATION AND MELTING REDUCE JUST SNOW AMOUNT
C**** ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
      FMSI2 = MIN(dSNdML*MELT1, MSI1-ACE1I+DEW-MELT1) ! > 0.
      FMSI1 = XSI(1)*FMSI2+XSI(2)*(DEW-MELT1)
      IF (FMSI1 .LT. 0.) FHSI1 = HSI(2)*FMSI1*BYXSI(2)/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward from layer 2 -> layer 3
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)
      MSI1 = MSI1+DEW-MELT1-FMSI2 ! kg/m^2
      SNOW = MSI1-ACE1I
      IF (.not.QFIXR) GO TO 210
      GO TO 370
  140 CONTINUE
C**** DEW INCREASES ICE AMOUNT,  MELT1 > 0.
      IF (MSI1-ACE1I .GT. MELT1) GO TO 150
C**** ALL SNOW AND SOME ICE MELT
      FMSI1 = XSI(1)*(MSI1-ACE1I)+DEW-MELT1 ! upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1
      SNOW = 0.
      IF (FMSI2 .LE. 0.) THEN ! (if melting is greater than dew)
C****ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
        FHSI1 = HSI(2)*FMSI1*BYXSI(2)/MSI1 ! upward -> layer 1 from layer 2
        FHSI2 = HSI(3)*FMSI2*BYXSI(3)/MSI2 ! upward -> layer 2 from layer 3
        HSI(1) = HSI(1)-FHSI1
        HSI(2) = HSI(2)+(FHSI1-FHSI2)
        IF (.not.QFIXR) GO TO 170
        GO TO 370
      ENDIF
C**** ICE ADVECTION IS DOWNWARD INTO LAYER 3 FROM LAYER 2
C**** (if dew is greater than melting)
      IF (FMSI1 .LT. 0.) FHSI1 = HSI(2)*FMSI1*BYXSI(2)/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward -> layer 3 from layer 2
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)
      MSI1 = ACE1I
      SNOW = 0.    ! necessary?
      IF (.not.QFIXR) GO TO 210
      GO TO 370
  150 CONTINUE
C**** MELTING REDUCES SNOW AMOUNT
C**** ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
      CMPRS = MIN(dSNdML*MELT1, MSI1-ACE1I-MELT1) ! > 0.
      FMSI1 = DEW+XSI(1)*CMPRS-XSI(2)*MELT1
      FMSI2 = DEW+CMPRS ! > 0. downward into layer 3 from layer 2
      IF (FMSI1 .LT. 0.) FHSI1 = HSI(2)*FMSI1*BYXSI(2)/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 3
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward energy flux into layer 3
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)
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
      FMSI1 = XSI(1)*(MSI1-ACE1I)+DEW ! < 0. upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW ! < 0. upward into layer 2 from layer 3
      FHSI1 = HSI(2)*FMSI1*BYXSI(2)/MSI1 ! upward energy flux into layer 1
      FHSI2 = HSI(3)*FMSI2*BYXSI(3)/MSI2 ! upward energy flux into layer 2
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)
      SNOW = 0. ! all snow evaporated
      IF (QFIXR) GO TO 370
  170 CONTINUE
      MSI1 = ACE1I
      IF (HSI(4)/LHM+XSI(4)*MSI2 .LE. 0.) GO TO 180 ! go to freezing case
C**** FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
      MELT4 = HSI(4)/LHM+XSI(4)*MSI2 ! > 0. melted ice from layer 4
      FMSI3 = XSI(4)*FMSI2+XSI(3)*MELT4
      IF (FMSI3 .LE. 0.) FHSI3 = -LHM*FMSI3 ! upward into layer 3
      IF (FMSI3 .GT. 0.) FHSI3 = HSI(3)*FMSI3*BYXSI(3)/MSI2 ! downward
      HSI(3) = HSI(3)+(FHSI2-FHSI3)
      HSI(4) = HSI(4)+FHSI3
      MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      GO TO 230
  180 CONTINUE
C**** NO ICE MELTS IN LAYER 4
C     FMSI3 = XSI(4)*FMSI2 ! upward mass flux into layer 3 from layer 4
      FHSI3 = HSI(4)*FMSI2/MSI2 ! upward energy flux into layer 3
      HSI(3) = HSI(3)+(FHSI2-FHSI3)
      HSI(4) = HSI(4)+FHSI3
      MSI2 = MSI2+FMSI2
      GO TO 230
  190 CONTINUE
C**** JUST SOME SNOW EVAPORATES
C**** NO ADVECTION BETWEEN LAYER 2 AND 3
C     FMSI1 = XSI(2)*DEW ! < 0. upward mass flux into layer 1
C     FMSI2 = 0. ! no ice advection between layers 2 and 3
      FHSI1 = HSI(2)*DEW/MSI1 ! upward energy flux into layer 1
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+FHSI1
      MSI1 = MSI1+DEW ! new ice mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass
      IF (QFIXR) GO TO 370
      IF (HSI(4)/LHM+XSI(4)*MSI2 .LE. 0.) GO TO 230 ! go to freezing case
C**** FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
      MELT4 = HSI(4)/LHM+XSI(4)*MSI2 ! melted ice from layer 4
C     FMSI3 = XSI(3)*MELT4 ! > 0. downward mass flux into layer 4
      FHSI3 = HSI(3)*MELT4/MSI2 ! downward heat flux into layer 4
      HSI(3) = HSI(3)-FHSI3
      HSI(4) = HSI(4)+FHSI3
      MSI2 = MSI2-MELT4 ! new ice mass of the second physical layer
      GO TO 230
  200 CONTINUE
C**** DEW INCREASES ICE AMOUNT, ADVECT ICE DOWNWARD
C**** ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
C     FMSI1 = DEW ! > 0. downward mass flux into layer 2
      FMSI2 = DEW ! > 0. downward mass flux into layer 3
      FHSI1 = HSI(1)*FMSI2/(XSI(1)*MSI1+DEW) ! downward heat flux
      FHSI2 = HSI(2)*FMSI2*BYXSI(2)/MSI1 ! downward heat flux into layer 3
      HSI(1) = HSI(1)-FHSI1
      HSI(2) = HSI(2)+(FHSI1-FHSI2)
      IF (QFIXR) GO TO 370
  210 CONTINUE
      IF (HSI(4)/LHM+XSI(4)*MSI2 .GT. 0.) THEN
C**** FLUXES HEAT UP TSIL(2) TO FREEZING POINT AND MELT SOME ICE
        MELT4 = HSI(4)/LHM+XSI(4)*MSI2 ! melted ice from layer 4
C     FMSI3 = XSI(4)*FMSI2+XSI(3)*MELT4 ! > 0. downward into layer 4
        FHSI3 = HSI(3)*((XSI(4)/XSI(3))*FMSI2+MELT4)/MSI2 ! downward
        HSI(3) = HSI(3)+(FHSI2-FHSI3)
        HSI(4) = HSI(4)+FHSI3
        MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      ELSE
C**** NO ICE MELTS IN LAYER 4
C     FMSI3 = XSI(4)*FMSI2 ! > 0. downward mass flux into layer 4
        FHSI3 = HSI(3)*(XSI(4)/XSI(3))*FMSI2/MSI2 ! downward heat flux
        HSI(3) = HSI(3)+(FHSI2-FHSI3)
        HSI(4) = HSI(4)+FHSI3
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

      SUBROUTINE ADDICE (SNOW,ROICE,TSIL,MSI2,HSI,DIFSI,EDIFSI,ENRGFO
     *     ,ACEFO,ACE2F,ENRGFI,TFW,FLEAD,QFIXR,QCMPR)
!@sum  ADDICE adds ice formed in the ocean to ice variables
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE

      REAL*8, PARAMETER, DIMENSION(LMI) :: YSI =
     *     (/XSI(1)* ACE1I/(ACE1I+AC2OIM),XSI(2)*ACE1I/(ACE1I+AC2OIM),
     *       XSI(3)*AC2OIM/(ACE1I+AC2OIM),XSI(4)*AC2OIM/(ACE1I+AC2OIM)/)
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
      REAL*8, INTENT(IN) ::  ENRGFI, ENRGFO, ACEFO, ACE2F
!@var ROICE,SNOW,MSI1,MSI2
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: HSI
      REAL*8, INTENT(OUT) :: EDIFSI, DIFSI
      REAL*8, INTENT(OUT), DIMENSION(LMI) :: TSIL

      REAl*8 FMSI1, FMSI2, FMSI3, FMSI4, FHSI1, FHSI2, FHSI3, FHSI4
      REAL*8 ROICEN, OPNOCN, DRSI, DRI, MSI1
      REAL*8 DIFS1

      EDIFSI=0. ; DIFSI=0.
      IF (.not.QFIXR .and. ROICE.LE.0. .and. ACEFO.gt.0) THEN
        ROICE=ACEFO/(ACE1I+AC2OIM)
        SNOW=0.
        TSIL=TFW
        MSI2=AC2OIM
        RETURN
      END IF
C****
      TSIL=0
      IF (ROICE.le.0) RETURN

      MSI1=SNOW+ACE1I
      IF (QFIXR) GO TO 370
      IF (ACE2F.le.0) GO TO 250 ! go to no freezing case
C**** CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
C     FMSI3 = -XSI(3)*ACE2F ! < 0.
C     FMSI4 = -ACE2F
      FHSI3 = -HSI(4)*ACE2F*(XSI(3)/XSI(4))/MSI2
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      IF (ACEFO .GT. 0.) GO TO 240
C**** NEW ICE IS FORMED BELOW OLD SEA ICE
      HSI(3) = HSI(3)-FHSI3
      HSI(4) = HSI(4)+(FHSI3+ENRGFI)
      MSI2 = MSI2+ACE2F ! new ice mass of physical layer 2
      GO TO 270
  240 CONTINUE
C**** NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACE2F))/(ROICE+DRSI) ! layer 2
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
      HSI(1) = ((1.-ROICE)*ENRGFO*YSI(1)+ROICE*HSI(1))/(ROICE+DRSI)
      HSI(2) = ((1.-ROICE)*ENRGFO*YSI(2)+ROICE*HSI(2))/(ROICE+DRSI)
      HSI(3) = ((1.-ROICE)*ENRGFO*YSI(3)+ROICE*(HSI(3)-FHSI3))/
     A       (ROICE+DRSI)
      HSI(4) = ((1.-ROICE)*ENRGFO*YSI(4)+ROICE*(HSI(4)+FHSI3+ENRGFI))/
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
      HSI(1) = ((1.-ROICE)*ENRGFO*YSI(1)+ROICE*HSI(1))/(ROICE+DRSI)
      HSI(2) = ((1.-ROICE)*ENRGFO*YSI(2)+ROICE*HSI(2))/(ROICE+DRSI)
      HSI(3) = ((1.-ROICE)*ENRGFO*YSI(3)+ROICE*HSI(3))/(ROICE+DRSI)
      HSI(4) = ((1.-ROICE)*ENRGFO*YSI(4)+ROICE*HSI(4))/(ROICE+DRSI)
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
C     FMSI3 = XSI(3)*FMSI4 ! < 0. upward ice mass flux into layer 3
      FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! upward ice mass into layer 4
      FHSI3 = HSI(4)*FMSI4*(XSI(3)/XSI(4))/MSI2 ! upward heat flux
      FHSI4 = (HSI(1)+HSI(2)+HSI(3)+HSI(4))*(DRSI/ROICEN)
      HSI(3) = HSI(3)-FHSI3
      HSI(4) = HSI(4)+(FHSI3-FHSI4)
      MSI2 = MSI2-FMSI4 ! new ice mass of second physical layer
C     SNOW = SNOW   ! snow thickness is conserved
      DIFS1 = DRI*ACE1I/ROICE
      EDIFSI = ROICE*(HSI(1)+HSI(2))*(DIFS1/MSI1)*0.5
      DIFSI = ROICE*DIFS1
      ROICE = ROICEN
      END IF
C**** RESAVE PROGNOSTIC QUANTITIES
  360 CONTINUE
  370 CONTINUE
      TSIL(1) = (HSI(1)/(XSI(1)*MSI1) +LHM)*BYSHI ! temperature of layer 1
      TSIL(2) = (HSI(2)/(XSI(2)*MSI1) +LHM)*BYSHI ! temperature of layer 2
      TSIL(3) = (HSI(3)/(XSI(3)*MSI2) +LHM)*BYSHI ! temperature of layer 3
      TSIL(4) = (HSI(4)/(XSI(4)*MSI2) +LHM)*BYSHI ! temperature of layer 4

      RETURN
      END SUBROUTINE ADDICE

      SUBROUTINE SIMELT(ROICE,SNOW,MSI2,TSIL,ENRGW,ENRGUSED)
!@sum  SIMELT melts sea ice if surrounding water is warm
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
!@var ENRGW energy available to melt ice (J/m^2)
      REAL*8, INTENT(IN) :: ENRGW
!@var ROICE,SNOW,MSI2 ice variables (%,kg/m^2,kg/m^2)
      REAL*8, INTENT(INOUT) :: ROICE, SNOW, MSI2
!@var TSIL ice temperatures (C)
      REAL*8, INTENT(INOUT), DIMENSION(LMI) :: TSIL
!@var ENRGUSED energy used to melt ice (J/m^2)
      REAL*8, INTENT(OUT) :: ENRGUSED
      REAL*8, DIMENSION(LMI) :: HSI
      REAL*8 MSI1,MSI2,MELT,DRSI,ROICEN,FHSI4,FHSI3,ENRGI
c      REAL*8 E_BOTTOM,GAMMA,HCRIT,HICE,ACE

      MSI1 = SNOW + ACE1I
C**** CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI(1) = (SHI*TSIL(1)-LHM)*XSI(1)*MSI1
      HSI(2) = (SHI*TSIL(2)-LHM)*XSI(2)*MSI1
      HSI(3) = (SHI*TSIL(3)-LHM)*XSI(3)*MSI2
      HSI(4) = (SHI*TSIL(4)-LHM)*XSI(4)*MSI2
      ENRGI = HSI(1)+HSI(2)+HSI(3)+HSI(4) ! energy in sea ice [J/m^2]
      IF (ROICE*ENRGI+ENRGW.LT.0.) GO TO 230
C**** THE WARM OCEAN MELTS ALL THE SNOW AND ICE
      ENRGUSED=-ROICE*ENRGI     ! only some energy is used
      ROICE=0.
      RETURN
C**** THE WARM OCEAN COOLS TO 0 DEGREES MELTING SOME SNOW AND ICE
C**** Reduce the ice depth
 230  ENRGUSED=ENRGW            ! all energy is used
C**** This is not currently used, but it could be
C**** I think this prescription is wrong though
c      ACE = SNOW + ACE1I + MSI2 ! total ice mass
c      HCRIT = 1.                ! critical ice thickness [m]
c      HICE = ACE/RHOI           ! sea ice thickness [m]
c      GAMMA = MAX(1d0,HICE/HCRIT)  
c      E_BOTTOM = ENRGW*(ROICE/(1.+GAMMA)) ! for bottom melting
C**** MELT ICE VERTICALLY, AND THEN HORIZONTALLY
      MELT = -XSI(4)*MSI2*ENRGW/HSI(4) ! melted ice at the bottom
c      MELT = -XSI(4)*MSI2*E_BOTTOM/HSI(4) ! melted ice at the bottom
      IF (MSI2-MELT .LT. AC2OIM) MELT = MSI2-AC2OIM
C     FMSI3 = XSI(3)*MELT ! > 0.
      FHSI3 = HSI(3)*MELT/MSI2
      FHSI4 = HSI(4)*MELT*BYXSI(4)/MSI2
C**** MELT SOME ICE HORIZONTALLY WITH REMAINING ENERGY
      DRSI = (ENRGW+ROICE*FHSI4)/(HSI(1)+HSI(2)+HSI(3)+HSI(4)-FHSI4) ! < 0.
      ROICEN = ROICE+DRSI       ! new sea ice concentration
c     SNOW = SNOW*(ROICE/ROICEN) ! new snow ice mass
      MSI2 = MSI2-MELT
      HSI(3) = HSI(3)-FHSI3
      HSI(4) = HSI(4)+FHSI3-FHSI4
C**** CONVERT SEA ICE ENTHALPY MINUS LATENT HEAT INTO TEMPERATURE
      TSIL(1) = (HSI(1)/(XSI(1)*MSI1) +LHM)*BYSHI ! temperatute of layer 1
      TSIL(2) = (HSI(2)/(XSI(2)*MSI1) +LHM)*BYSHI ! temperatute of layer 2
      TSIL(3) = (HSI(3)/(XSI(3)*MSI2) +LHM)*BYSHI ! temperatute of layer 3
      TSIL(4) = (HSI(4)/(XSI(4)*MSI2) +LHM)*BYSHI ! temperatute of layer 4
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
      REAL*8, DIMENSION(LMI,IM,JM) :: TSI
!@var HSI enthaply of each ice layer (J/m^2) (will replace TSI)
c      REAL*8, DIMENSION(LMI,IM,JM) :: HSI

      END MODULE SEAICE_COM

      SUBROUTINE io_seaice(kunit,iaction,ioerr)
!@sum  io_seaice reads and writes seaice variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
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
         WRITE (kunit,err=10) MODULE_HEADER,RSI,TSI,SNOWI,MSI
      CASE (IOREAD:)            ! input from restart file
         READ (kunit,err=10) HEADER,RSI,TSI,SNOWI,MSI
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
      USE E001M12_COM
      USE SEAICE_COM, only : rsi,msi,tsi,snowi
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in ice data
      CALL CHECK3(RSI,IM,JM,1,SUBR,'rsi')
      CALL CHECK3(MSI,IM,JM,1,SUBR,'msi')
      CALL CHECK3(TSI,4,IM,JM,SUBR,'tsi')
      CALL CHECK3(SNOWI,IM,JM,1,SUBR,'sni')

      END SUBROUTINE CHECKI
