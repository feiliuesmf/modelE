C**** PE001M12 E001M12 SOMTQ PB395M12
C**** OPT(3)
C**** semi-random cloud overlap (+snow age updates+computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C**** Constant pressure at L=LS1 and above (SIGE(LS1)=0., PLE(LS1)=PTOP)
C**** Using 5 harmonics for horizontal ocean heat transport, thinner ice
C**** Routines included:  PRECIP, COSZ0, RADIA,
C****                     GROUND, DRYCNV, SDRAG
*****
C*    Sea ice has four thermal layers
C*    The lead fraction is ice thickness dependent
C*
      SUBROUTINE PRECIP
C****
C**** THIS SUBROUTINE USES THE PRECIPITATION TO CALCULATE THE GROUND
C**** WATER, GROUND ICE, SNOW COVER, AND RUNOFF
C****
C**** RUN1 IS NOT ACUMULATED IN ADAILY FOR DIAG6
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,rhow,rhoi,shw,shi
      USE E001M12_COM
      USE GEOM
      USE CLD01_COM_E001, only : PREC,TPREC
      USE DAGCOM, only : aj,bj,cj,dj,aij,jreg
      USE OCEAN, only : ODATA,OA,XSI1,XSI2,XSI3,XSI4,R2,R3,TTRUNC,Z1I
     *     ,Z2OIM,ACE1I,AC2OIM,TFO

      IMPLICIT REAL*8 (A-H,O-Z)
C*
      PARAMETER (SNOMAX=100.0, dSNdRN=0.)
C*
      REAL*8 MSI1, MSI2, MELT1

      DATA Z1E/.1/,Z2LI/2.9/

      DATA IFIRST/1/
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** ODATA  1  OCEAN TEMPERATURE (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        2  EARTH SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)
C****       15  OCEAN ICE TEMPERATURE OF THIRD LAYER (C)
C****       16  OCEAN ICE TEMPERATURE OF FOURTH LAYER (C)
C****
C*    WTRO   - Liquid ocean mass (kg/m^2)
C*
C*    ROICE  - Horizontal ratio of sea ice to ocean (1)
C*    MSI    - Mass of sea ice (kg/m^2)
C*    HSI    - Enthalpy minus latent heat of sea ice (J/m^2)
C*
C*    POCEAN - Ocean fraction (1)
C*
C*    PREC   - Precipitation from atmosphere (kg/m^2)
C*    EPRCP  - Energy of precipitation (J/m^2)
C*
      IF (IFIRST.NE.1) GO TO 10
      IFIRST=0
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      ACE2LI=Z2LI*RHOI
      HC1I=ACE1I*SHI
      HC1DE=Z1E*1129950.
C*
   10 CONTINUE
C*
      DO 980 J=1,JM
      IMAX=IMAXJ(J)
         AENRGP=0.
         BENRGP=0.
         CENRGP=0.
         BEDIFS=0.
         CEDIFS=0.
         AEFO=0.
         BERUN0=0.
         BERUN2=0.
         CERUN2=0.
         AERUN4=0.
         CERUN4=0.
         BDIFS=0.
         CDIFS=0.
         AIFO=0.
         BRUN0=0.
         CRUN0=0.
         BRUN2=0.
         CRUN2=0.
         ARUN4=0.
         CRUN4=0.
      DO 960 I=1,IMAX
      IF (PREC(I,J).LE.0.) GO TO 960
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      ROICE=ODATA(I,J,2)
      POICE=ROICE*PWATER
      POCEAN=PWATER-POICE
         JR=JREG(I,J)
         DXYPJ=DXYP(J)
         RUN0S=0.
         DIFSS=0.
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
      PRCP=PREC(I,J)
      TPRCP=TPREC(I,J)
      IF (TPRCP.LT.0.) GO TO 30
C     EPRCP=PRCP*TPRCP*SHW
      EPRCP=0.
      ENRGP=EPRCP
      GO TO 50
C     EPRCP=PRCP*TPRCP*SHI
   30 EPRCP=0.
      ENRGP=EPRCP-PRCP*LHM
         AIJ(I,J,70)=AIJ(I,J,70)+PRCP
   50 CONTINUE
      IF (PWATER.LE.0.) GO TO 400
C****
C**** OCEAN
C****
            OA(I,J,4)=OA(I,J,4)+ENRGP
         AENRGP=AENRGP+ENRGP*POCEAN
         AIJ(I,J,65)=AIJ(I,J,65)+ENRGP*POCEAN
      IF (KOCEAN.NE.1) GO TO 100
      TGW=ODATA(I,J,1)
      WTRO=Z1O(I,J)*RHOW
      ENRGO0=WTRO*TGW*SHW
      EOFRZ=WTRO*TFO*SHW
      RUN4=PRCP
      ERUN4=RUN4*TGW*SHW
         AERUN4=AERUN4+ERUN4*POCEAN
         ARUN4=ARUN4+RUN4*POCEAN
      ENRGO=ENRGP-ERUN4
      IF (ENRGO0+ENRGO.LT.EOFRZ) GO TO 80
C**** OCEAN TEMPERATURE IS STILL ABOVE FREEZING, NO ICE IS FORMED
      ACEFO=0.
      ENRGFO=0.
      IF (ROICE.GT.0.) GO TO 110
      ODATA(I,J,1)=TGW+(ENRGO/(WTRO*SHW)+TTRUNC)
      GO TO 400
C**** SNOW COOLS TGO TO FREEZING POINT FOR OCEAN AND FORMS SOME ICE
   80 IF (ROICE.GT.0.) GO TO 110
      ODATA(I,J,1)=TGW+(ENRGO/(WTRO*SHW)+TTRUNC)
      GO TO 400
C****
  100 IF (POICE.LE.0.) GO TO 400
C****
C**** OCEAN ICE
C****
  110 SNOW=GDATA(I,J,1)
      MSI2=ODATA(I,J,3)
      MSI1 = SNOW + ACE1I
      EPRE = EPRCP - PRCP*LHM ! total energy of precipitation;
C***                            EPRCP - sensible heat of precipitation;
C***                            -PRCP*LHM - latent heat of precipit.
C*
      TG1 = GDATA(I,J,3)  ! first layer sea ice temperature
      TG2 = GDATA(I,J,7)  ! second layer sea ice temperature
      TG3 = GDATA(I,J,15) ! third layer sea ice temperature
      TG4 = GDATA(I,J,16) ! fourth layer sea ice temperature
C*
C***  CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
C*
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2
      HSI4 = (SHI*tG4-LHM)*XSI4*MSI2
C*
         CENRGP=CENRGP+ENRGP*POICE
         AIJ(I,J,66)=AIJ(I,J,66)+ENRGP*POICE
      HC1=HC1I+SNOW*SHI
      HC_1 = XSI1*MSI1*SHI
      RUN0=0.
      IF (KOCEAN.NE.1) GO TO 120
         CERUN4=CERUN4+ERUN4*POICE
         CRUN4=CRUN4+RUN4*POICE
      HC2=MSI2*SHI
      WTRW0=WTRO-ROICE*(SNOW+ACE1I+MSI2)
      ENRGW0=WTRW0*TGW*SHW
      DIFS=0.
      EDIFS=0.
  120 HSI1 = HSI1+EPRE
      IF (TPRCP.LT.0.) GO TO 180
      IF (EPRCP.LT.-TG1*HC_1) GO TO 160
C*
C***  ALL PRECIPITATION IS RAIN ABOVE 0degC
C***  RAIN COMPRESSES SNOW INTO ICE
C*
      RAIN = PRCP
      IF (HSI1/LHM+XSI1*MSI1 .LE. 0.) GO TO 140
C*
C***  WARM RAIN MELTS SOME SNOW OR ICE
C*
      MELT1 = HSI1/LHM+XSI1*MSI1
      DWATER = MELT1
      RUN0=DWATER+PRCP
      IF (MSI1-ACE1I .LE. DWATER) GO TO 130
C*
C**** RAIN MELTS SOME SNOW AND COMPRESSES SNOW INTO ICE
C*
      FMSI2 = MIN (dSNdRN*(RAIN+DWATER), MSI1-ACE1I-DWATER) ! > 0.
      FHSI1 = -LHM*(XSI1*FMSI2-XSI2*DWATER) ! downward heat flux
c     F1 = HSI1*(XSI1*FMSI2-XSI2*DWATER)/(XSI1*MSI1-DWATER)
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! > 0.
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = MSI1-DWATER-FMSI2
      SNOW = MSI1-ACE1I
      IF (SNOW .LT. 0.) SNOW = 0.
      HSI1 = HSI1-FHSI1
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      DIFS = FMSI2
      EDIFS = FHSI2
      ERUN2 = HSI3*FMSI2*R3/MSI2
      IF (KOCEAN.NE.1) GO TO 200
      GO TO 210
C*
  130 CONTINUE
C*
C**** RAIN MELTS ALL SNOW AND SOME ICE
C*
      FMSI2 = MSI1-ACE1I-DWATER ! < 0.(upward ice mass flux)
      DIFS = FMSI2 ! < 0.
C     FMSI1 = XSI1*(MSI1-ACE1I)-DWATER ! < 0.(melted snow/ice mass)
      FHSI1 = HSI2*((XSI1/XSI2)*FMSI2-MELT1)/MSI1 !  upward heat flux
      FHSI2 = HSI3*FMSI2*R3/MSI2 !  upward heat flux into layer 2
      FHSI3 = HSI4*FMSI2/MSI2 !  upward heat flux into layer 3
      SNOW=0. ! Rain melted all snow
      MSI1 = ACE1I ! Keep the first layer ice mass constant ACE1I
      HSI1 = HSI1-FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      EDIFS=FHSI2 ! for diagnostics
      ERUN2=EDIFS ! for diagnostics
      IF (KOCEAN.NE.1) GO TO 200
      GO TO 210
C*
  140 CONTINUE
C*
C**** RAIN COMPRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
C*
      CMPRS = MIN(dSNdRN*RAIN, MSI1-ACE1I) ! all snow or part of rain
      IF (-HSI1/LHM-XSI1*MSI1 .LT. RAIN) GO TO 150
C*
C***  ALL RAIN FREEZES IN LAYER 1
C*
C     FREZ1 = RAIN ! frozen rain
C     FMSI1 = XSI1*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      RUN0 = 0.
      FHSI1 = HSI1*(XSI1*CMPRS+RAIN)/(XSI1*MSI1+RAIN) ! downward
      F1 = HSI1*(XSI1*CMPRS+RAIN)/(XSI1*MSI1) ! for prescribed ice
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux from layer 2
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
      ERUN2 = HSI3*FMSI2*R3/MSI2
      IF (KOCEAN.NE.1) GO TO 200
      GO TO 210
C*
  150 CONTINUE
C*
C**** JUST PART OF RAIN FREEZES IN LAYER 1
C*
      FREZ1 = -HSI1/LHM-XSI1*MSI1 ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI1*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI1*CMPRS+FREZ1) ! downward heat flux from layer 1
      F1 = HSI1*(XSI1*CMPRS+FREZ1)/(XSI1*MSI1) ! for prescribed ice
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux from layer 2
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
      ERUN2 = HSI3*FMSI2*R3/MSI2
      IF (KOCEAN.NE.1) GO TO 200
      GO TO 210
C*
  160 CONTINUE
C*
C***  PRECIPITATION IS A MIXTURE OF RAIN AND SNOW AT 0 degC
C***  RAIN COMRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
C*
      SNOW = -EPRCP/LHM ! snow fall
      RAIN = PRCP-SNOW  ! rain fall
      CMPRS = MIN(dSNdRN*RAIN, MSI1+SNOW-ACE1I) ! compression
      IF (-HSI1/LHM-XSI1*MSI1-SNOW .LT. RAIN) GO TO 170
C*
C***  ALL RAIN FREEZES IN LAYER 1
C*
C     FREZ1 = RAIN ! frozen rain
      RUN0 =0.
      FMSI1 = XSI1*CMPRS+XSI2*SNOW+RAIN ! downward ice mass flux
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSI1*FMSI1/(XSI1*MSI1+PRCP) ! downward heat flux
      F1 = HSI1*FMSI1/(XSI1*MSI1) ! for prescribed ice/ocean
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux from layer 2
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
      ERUN2 = HSI3*FMSI2*R3/MSI2
      IF (KOCEAN .NE. 1) GO TO 200
      GO TO 210
C*
  170 CONTINUE
C*
C***  NOT ALL RAIN FREEZES IN LAYER 1
C*
      FREZ1 = -HSI1/LHM-XSI1*MSI1-SNOW ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI1*CMPRS+XSI2*SNOW+FREZ1 ! downward ice mass flux
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI1*CMPRS+XSI2*SNOW+FREZ1) ! downward
      F1 = HSI1*(XSI1*CMPRS+XSI2*SNOW+FREZ1)/(XSI1*MSI1) ! prescribed
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux from layer 2
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
      ERUN2 = HSI3*FMSI2*R3/MSI2
      IF (KOCEAN .NE. 1) GO TO 200
      GO TO 210
C*
  180 CONTINUE
C*
C***  ALL PRECIPITATION IS SNOW, SNOW AMOUNT INCREASES
C*
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
      IF (KOCEAN .NE. 1) GO TO 380
      GO TO 300
C*
  190 CONTINUE
C*
C***  TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
C*
      FMSI2 = MSI1+PRCP-(0.9*SNOMAX+ACE1I) ! > 0.(compressed snow)
C     FMSI1 = XSI2*PRCP+XSI1*FMSI2 ! > 0.(downward ice mass flux)
      FHSI1 = HSI1*(XSI2*PRCP+XSI1*FMSI2)/(XSI1*MSI1+PRCP) ! downward
      F1 = HSI1*(XSI2*PRCP+XSI1*FMSI2)/(XSI1*MSI1) ! for prescribed
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = 0.9*SNOMAX+ACE1I ! first layer ice mass
      SNOW=.9*SNOMAX ! snow mass
C*
      DIFS = FMSI2 ! compressed snow mass
      EDIFS = FHSI2 ! energy of compressed snow mass
C*
      HSI1 = HSI1-FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      IF (KOCEAN .EQ. 1) GO TO 210
C*
      ERUN2 = HSI3*FMSI2*R3/MSI2
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
C*
  200 CONTINUE
C*
C***  ACCUMULATE DIAGNOSTICS
C*
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
         CERUN2=CERUN2+ERUN2*POICE
         CRUN2=CRUN2+DIFS*POICE
      GO TO 380
C*
  210 CONTINUE
C*
C**** ADVECT ICE (usually downwards)
C*
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
C*
  220 CONTINUE
C*
      HSI2 = HSI2+(FHSI1-FHSI2) ! for predicted ice/ocean
      HSI3 = HSI3+(FHSI2-FHSI3) ! for predicted ice/ocean
      HSI4 = HSI4+FHSI3 ! for predicted ice/ocean
      TG2 = (HSI2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      TG3 = (HSI3/(XSI3*MSI2)+LHM)/SHI ! third layer ice temperature
      TG4 = (HSI4/(XSI4*MSI2)+LHM)/SHI ! fourth layer ice temperature
C*
      MSI2 = MSI2+FMSI2 ! second layer sea ice mass (kg/m^2)
C*
C**** CALCULATE THE COMPOSITE WATER MASS AND WATER ENERGY
C*
  300 WTRW=WTRW0+ROICE*(RUN0-RUN4)
      ENRGW=ENRGW0-ROICE*ERUN4+(1.-ROICE)*ENRGP
      TGW=ENRGW/(WTRW*SHW) ! mixed layer temperature
      ODATA(I,J,1)=TGW
      ODATA(I,J,2)=ROICE
      ODATA(I,J,3)=MSI2
      GDATA(I,J,15) = TG3
      GDATA(I,J,16) = TG4
  380 GDATA(I,J,1)=SNOW
      GDATA(I,J,7)=TG2
  390 GDATA(I,J,3)=TG1
         CRUN0=CRUN0+RUN0*POICE
         RUN0S=RUN0S+RUN0*POICE
C****
  400 IF (PLICE.LE.0.) GO TO 600
C****
C**** LAND ICE
C****
      SNOW=GDATA(I,J,12)
      TG1=GDATA(I,J,13)
      TG2=GDATA(I,J,14)
         BENRGP=BENRGP+ENRGP*PLICE
         AIJ(I,J,67)=AIJ(I,J,67)+ENRGP
      HC1=HC1I+SNOW*SHI
      RUN0=0.
      IF (TPRCP.LT.0.) GO TO 480
      IF (EPRCP.LT.-TG1*HC1) GO TO 460
C**** RAIN HEATS UP TG1 TO FREEZING POINT AND MELTS SOME SNOW OR ICE
      DWATER=(TG1*HC1+EPRCP)/LHM
      TG1=0.
      RUN0=DWATER+PRCP
      IF (DWATER.GT.SNOW) GO TO 440
C**** RAIN MELTS SOME SNOW
      SNOW=SNOW-DWATER
      GO TO 580
C**** RAIN MELTS ALL SNOW AND SOME ICE, ICE MOVES UP THROUGH THE LAYERS
  440 DIFS=SNOW-DWATER
      SNOW=0.
      TG1=-TG2*DIFS/ACE1I
      EDIFS=DIFS*(TG2*SHI-LHM)
      ERUN2=EDIFS
      GO TO 560
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1
  460 TG1=TG1+EPRCP/HC1
      RUN0=PRCP
      GO TO 590
C**** SNOW INCREASES SNOW AMOUNT AND SNOW TEMPERATURE RECOMPUTES TG1
  480 TG1=(TG1*HC1+EPRCP)/(HC1+PRCP*SHI)
      SNOW=SNOW+PRCP
      IF (SNOW.LE.ACE1I) GO TO 580
C**** SNOW IS COMPACTED INTO ICE, ICE MOVES DOWN THROUGH THE LAYERS
      DIFS=SNOW-.9*ACE1I
      SNOW=.9*ACE1I
      EDIFS=DIFS*(TG1*SHI-LHM)
      ERUN2=DIFS*(TG2*SHI-LHM)
      GDATA(I,J,14)=TG2+(TG1-TG2)*DIFS/ACE2LI
  560    BEDIFS=BEDIFS+EDIFS*PLICE
         AIJ(I,J,69)=AIJ(I,J,69)+EDIFS
         BDIFS=BDIFS+DIFS*PLICE
         DIFSS=DIFSS+DIFS*PLICE
         BERUN2=BERUN2+ERUN2*PLICE
         AIJ(I,J,72)=AIJ(I,J,72)+ERUN2
         BRUN2=BRUN2+DIFS*PLICE
  580 GDATA(I,J,12)=SNOW
  590 GDATA(I,J,13)=TG1
         BRUN0=BRUN0+RUN0*PLICE
         RUN0S=RUN0S+RUN0*PLICE
         AIJ(I,J,33)=AIJ(I,J,33)+RUN0
C****
  600 IF (PEARTH.LE.0.) GO TO 940
C****
C**** EARTH  (all else is done in subroutine EARTH, called from SURFCE)
C****
         BENRGP=BENRGP+ENRGP*PEARTH
         AIJ(I,J,68)=AIJ(I,J,68)+ENRGP
C****
C**** ACCUMULATE DIAGNOSTICS (ocean, ocean ice, land ice only)
C****
  940    DJ(JR,39)=DJ(JR,39)+ENRGP*DXYPJ
         DJ(JR,45)=DJ(JR,45)+DIFSS*DXYPJ  ! ocn/land ice contribution
         DJ(JR,54)=DJ(JR,54)+RUN0S*DXYPJ
         AIJ(I,J,5)=AIJ(I,J,5)+PREC(I,J)
         AIJ(I,J,23)=AIJ(I,J,23)+ENRGP
  960 CONTINUE
         AJ(J,39)=AJ(J,39)+AENRGP
         BJ(J,39)=BJ(J,39)+BENRGP
         CJ(J,39)=CJ(J,39)+CENRGP
C        BJ(J,40)=BJ(J,40)+BERUN0         ! land ice (Tg=0)
         BJ(J,41)=BJ(J,41)+BEDIFS         ! land ice contribution
         CJ(J,41)=CJ(J,41)+CEDIFS
         AJ(J,43)=AJ(J,43)+AEFO
         BJ(J,43)=BJ(J,43)+BERUN2         ! land ice contribution
         CJ(J,43)=CJ(J,43)+CERUN2
         BJ(J,45)=BJ(J,45)+BDIFS          ! land ice contribution
         CJ(J,45)=CJ(J,45)+CDIFS
         BJ(J,54)=BJ(J,54)+BRUN0          ! land ice contribution
         CJ(J,54)=CJ(J,54)+CRUN0
         AJ(J,46)=AJ(J,46)+AIFO
         BJ(J,46)=BJ(J,46)+BRUN2          ! land ice contribution
         CJ(J,46)=CJ(J,46)+CRUN2
         AJ(J,47)=AJ(J,47)+ARUN4
         CJ(J,47)=CJ(J,47)+CRUN4
         AJ(J,48)=AJ(J,48)+AERUN4
         CJ(J,48)=CJ(J,48)+CERUN4
  980 CONTINUE
      RETURN
      END
      SUBROUTINE COSZ0
C****
C**** calculates the Earth's zenith angle, weighted either
C**** by time or by sun light.
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LT1(IM),LT2(IM)
      DIMENSION COSZ(IM,JM),COSZA(IM,JM)
      DIMENSION SINJ(JM),COSJ(JM),RI(IM),SINI(IM),COSI(IM)
      COMMON/WORK5/LT1,LT2,SLT1(IM),SLT2(IM),S2LT1(IM),S2LT2(IM)
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      DATA ZERO1/1.D-2/
C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      PHIS=-.25*TWOPI
      SPHIS=-1.
      CPHIS=0.
      DO 20 J=1,JM-1
      PHIN=(TWOPI/(JMM1+JMM1))*(J-.5*JM)
      SPHIN=SIN(PHIN)
      CPHIN=COS(PHIN)
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(J)=SIN(PHIM)
      COSJ(J)=COS(PHIM)
      PHIS=PHIN
      SPHIS=SPHIN
   20 CPHIS=CPHIN
      PHIN=.25*TWOPI
      SPHIN=1.
      CPHIN=0.
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(JM)=SIN(PHIM)
      COSJ(JM)=COS(PHIM)
C**** COMPUTE THE SINES AND COSINES OF LONGITUDE
      DO 40 I=1,IM
      RI(I)=(TWOPI/IM)*(I-.5)
      SINI(I)=SIN(RI(I))
   40 COSI(I)=COS(RI(I))
      RETURN
C****
C****
      ENTRY COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
  100 SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO 120 I=1,IM
      LT1(I)=ROT1+RI(I)
      SLT1(I)=SR1*COSI(I)+CR1*SINI(I)
      LT2(I)=ROT2+RI(I)
  120 SLT2(I)=SR2*COSI(I)+CR2*SINI(I)
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO 200 J=1,JM,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 180
      IF (SJSD-CJCD.GE.0.) GO TO 160
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
      COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
      GO TO 200
C**** CONSTANT DAYLIGHT NEAR THE POLES
  160 COSZ(1,J)=SJSD
      GO TO 200
C**** CONSTANT NIGHTIME NEAR THE POLES
  180 COSZ(1,J)=0.
  200 CONTINUE
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 500 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 460
      IF (SJSD-CJCD.GE.0.) GO TO 420
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 400 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 220
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  220 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 240
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      GO TO 400
  240 IF (DAWN.GE.LT1(I)) GO TO 300
      IF (DUSK.LT.LT2(I)) GO TO 260
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 400
  260 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/DROT
      GO TO 400
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  280 COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*
     *  (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/DROT
      GO TO 400
  300 IF (DUSK.LT.LT2(I)) GO TO 320
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/DROT
      GO TO 400
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  320 COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/DROT
  400 CONTINUE
      GO TO 500
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
  420 DO 440 I=1,IM
  440 COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 500
C**** CONSTANT NIGHTIME AT THIS LATITUDE
  460 DO 480 I=1,IM
  480 COSZ(I,J)=0.
  500 CONTINUE
      RETURN
C****
C****
      ENTRY COSZS (ROT1,ROT2,COSZ,COSZA)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.  COSZT MUST HAVE BEEN
C**** CALLED JUST PREVIOUSLY.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO 520 I=1,IM
      LT1(I)=ROT1+RI(I)
      SLT1(I)=SR1*COSI(I)+CR1*SINI(I)
      CLT1=CR1*COSI(I)-SR1*SINI(I)
      S2LT1(I)=2.*SLT1(I)*CLT1
      LT2(I)=ROT2+RI(I)
      SLT2(I)=SR2*COSI(I)+CR2*SINI(I)
      CLT2=CR2*COSI(I)-SR2*SINI(I)
  520 S2LT2(I)=2.*SLT2(I)*CLT2
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO 600 J=1,JM,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 580
      IF (SJSD-CJCD.GE.0.) GO TO 560
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
      ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *  .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
      COSZ(1,J)=ECOSZ/TWOPI
      COSZA(1,J)=ECOSQZ/ECOSZ
      GO TO 600
C**** CONSTANT DAYLIGHT NEAR THE POLES
  560 ECOSZ=SJSD*TWOPI
      ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI
      COSZ(1,J)=ECOSZ/TWOPI
      COSZA(1,J)=ECOSQZ/ECOSZ
      GO TO 600
C**** CONSTANT NIGHTIME NEAR THE POLES
  580 COSZ(1,J)=0.
      COSZA(1,J)=0.
  600 CONTINUE
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 900 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 860
      IF (SJSD-CJCD.GE.0.) GO TO 820
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 800 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 620
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  620 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 640
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      COSZA(I,J)=0.
      GO TO 800
  640 IF (DAWN.GE.LT1(I)) GO TO 700
      IF (DUSK.LT.LT2(I)) GO TO 660
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  660 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 680
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      ECOSZ=SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I))+
     *  .5*CJCD*(DUSK-LT1(I)+.5*(S2DUSK-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  680 ECOSZ=SJSD*(DROT-DAWN-TWOPI+DUSK)+
     *  CJCD*(SLT2(I)-SDAWN+SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I)+SLT2(I)-SDAWN)+
     *  .5*CJCD*(DUSK+DROT-DAWN-TWOPI+
     *  .5*(S2DUSK-S2LT1(I)+S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  700 IF (DUSK.LT.LT2(I)) GO TO 720
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+
     *  .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  720 ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *  .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
  800 CONTINUE
      GO TO 900
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
  820 DO 840 I=1,IM
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
  840 COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 900
C**** CONSTANT NIGHTIME AT THIS LATITUDE
  860 DO 880 I=1,IM
      COSZ(I,J)=0.
  880 COSZA(I,J)=0.
  900 CONTINUE
      RETURN
      END
      SUBROUTINE RADIA
C****
C**** THIS SUBROUTINES ADDS THE RADIATION HEATING TO THE TEMPERATURES
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,tf
      USE E001M12_COM
      USE GEOM
      USE RADNCB
      USE RE001
     &  , only : setnew,rcomp1,writer,rcompx,rcompt ! routines
     &             ,lx
     &             ,FULGAS ,PTLISO ,KTREND ,LMR=>NL ,LMRP=>NLP
C     INPUT DATA
     &             ,PLE=>PLB ,TL=>TLM ,QL=>SHL
     &             ,TAUWC ,TAUIC ,SIZEWC ,SIZEIC
     &             ,POCEAN,PEARTH,POICE,PLICE,AGESN,SNOWE,SNOWOI,SNOWLI
     &             ,TGO,TGE,TGOI,TGLI,TS=>TSL,WS=>WMAG,WEARTH,PTOPTR
     &             ,S00WM2,RATLS0,S0,COSZ,PVT
     &             ,JYEARR=>JYEAR,JDAYR=>JDAY,JLAT,ILON,KCLDEM,KVEGA6
C     OUTPUT DATA
     &             ,TRDFLB ,TRNFLB ,TRFCRL
     &             ,SRDFLB ,SRNFLB ,SRFHRL
     &             ,PLAVIS ,PLANIR ,ALBVIS ,ALBNIR ,FSRNFG
     &             ,SRRVIS ,SRAVIS ,SRRNIR ,SRANIR
     &             ,BTEMPW
C     Some Aerosol Control Parameters not currently referenced
c    &             ,QXAERO ,QSAERO ,QCAERO ,ATAERO
c    &             ,QAER55 ,REAERO ,VEAERO ,ROAERO ,PI0MAX
c    &             ,FSAERO ,FTAERO ,VDGAER ,SSBTAU ,PIAERO
      USE RANDOM
      USE CLD01_COM_E001, only : TAUSS,TAUMC,SVLHX,RHSAV,SVLAT,CLDSAV,
     *     CLDSS,CLDMC,CSIZE
      USE PBLCOM, only : wsavg,tsavg
      USE DAGCOM, only : aj,bj,cj,dj,jreg,aij,ail,ajl,asjl,adaily,
     *     iwrite,jwrite,itwrite
      USE DYNAMICS, only : pk,pedn
      USE OCEAN, only : odata

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/WORK1d/COSZ2(IM,JM),COSZA(IM,JM),
     *  TRINCG(IM,JM),BTMPW(IM,JM),SNFS(IM,JM,4),TNFS(IM,JM,4),
     *  TRHRS(IM,JM,3),SRHRS(IM,JM,3),ALB(IM,JM,9)
      COMMON/WORK2c/ TOTCLD(LM)

      DIMENSION COE(LM+3)
      LOGICAL POLE

      DATA TCIR/258.16/,STBO/.567257D-7/,IFIRST/1/,JDLAST/-9/
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** ODATA  1  OCEAN TEMPERATURE (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        2  EARTH SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)
C****        9  AGE OF SNOW OVER OCEAN ICE (DAYS)
C****       10  AGE OF SNOW OVER LAND ICE (DAYS)
C****       11  AGE OF SNOW OVER EARTH (DAYS)
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)
C****
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****
      QSAT(TM,PR,QLH)=3.797915*DEXP(QLH*(7.93252D-6-2.166847D-3/TM))/PR
         IF (MODRD.EQ.0) IDACC(2)=IDACC(2)+1
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
      CALL COSZ0
      DTCNDS=NCNDS*DT
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
      LMR=LM+3
      LMRP=LMR+1
      COEX=.01*GRAV*KAPA/RGAS
      DO L=1,LM
         COE(L)=DTCNDS*COEX/DSIG(L)
         PLE(L)=SIGE(L)*(PSF-PTOP)+PTOP
      END DO
      PLE(LM+1)=SIGE(LM+1)*(PSF-PTOP)+PTOP
      PLE(LM+2)=.5*PLE(LM+1)
      PLE(LMR)=.2*PLE(LM+1)
      PLE(LMR+1)=1.D-5
      PTOPTR=PTOP ! top of sigma-coord.system
      DO 40 LR=LM+1,LMR
   40 COE(LR)=DT*NRAD*COEX/(PLE(LR)-PLE(LR+1))
      PTLISO=15.
      IF(CO2.LT.0.) KTREND=-NINT(CO2)
C**** Default: time-dependent So/GHG/O3/Trop-Aeros/Dust/Volc-Aeros
C****     For control runs e.g. with Jul 1,1951 atmosphere use
          JYFIX=1951
          JDFIX=182   !  Julian date (if JDFIX=0, annual cycle is used)
        CALL SETNEW(11,JYFIX,JDFIX, 1,0,0.D0) ! fix sol.const. - KSOLAR
        CALL SETNEW( 2,JYFIX,JDFIX, 0,0,0.D0) ! fix GHG (trend KTREND)
        CALL SETNEW(13,0    ,0    , 0,0,0.D0) ! no GHG-resets - MADGAS
        CALL SETNEW( 3,JYFIX,0    , 0,0,0.D0) ! quasi-fix O3 (ann.cycle)
        CALL SETNEW( 4,JYFIX,0    , 0,0,0.D0) ! quasi-fix Trop. Aerosols
        CALL SETNEW( 5,JYFIX,0    , 0,0,0.D0) ! quasi-fix Desert Dust
        CALL SETNEW( 6,JYFIX,JDFIX, 0,0,0.D0) ! fix Volc. Aerosols
C     An annual cycle is used for the data below, to prevent this, use
CnoAC   CALL SETNEW(7, 0,JDFIX, 0,0,0.D0) ! cloud heterogeneity - KCLDEP
CnoAC   CALL SETNEW(8, 0,JDFIX, 0,0,0.D0) !  surface albedo
C**** New options (currently not used)
      KCLDEM=0  ! 0:old 1:new LW cloud scattering scheme  -  KCLDEM
      KVEGA6=0  ! 0:2-band 1:6-band veg.albedo            -  KVEGA6
      MADVEL=123456         ! suppress reading i-th time series by i->0
      CALL RCOMP1 (MADVEL)  ! 1=O3 2=TrAer 3=Dust 4=VAer 5=Clds 6=SoUV
      CO2REF=FULGAS(2)
      IF(CO2.GE.0.) FULGAS(2)=CO2REF*CO2
         CALL WRITER (6,0)
         INCHM=NRAD/NDYN
         JEQ=1+JM/2
         J50N=(50.+90.)*JMM1/180.+1.5
         J70N=(70.+90.)*JMM1/180.+1.5
C**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS
         DO 43 L=1,LM
         LLOW=L
         IF (.5*(PLE(L+1)+PLE(L+2)).LT.750.) GO TO 44 ! was 786. 4/16/97
   43    CONTINUE
   44    LMID1=LLOW+1
         DO 45 L=LMID1,LM
         LMID=L
         IF (.5*(PLE(L+1)+PLE(L+2)).LT.430.) GO TO 46
   45    CONTINUE
   46    LHI1=LMID+1
         LHI=LM
         IF (LHI1.GT.LHI) LHI=LHI1
         WRITE (6,47) LLOW,LMID1,LMID,LHI1,LHI
   47    FORMAT (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN',
     *     ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)
C**** NO RADIATION AVERAGING          IJRA=1   JRA=1   IRA=1
C**** RADIATION AVERAGING IN I             2       1       2
C**** RADIATION AVERAGING IN I AND J       4       2       2
      IF(IJRA.NE.1) STOP 'FSF not ready for radiation averaging'
      JRA=(IJRA+2)/3
      IRA=IJRA/JRA
   50 JALTER=MOD(NSTEP,NRAD*JRA)/NRAD
      IALTER=MOD(NSTEP,NRAD*IJRA)/(NRAD*JRA)
C**** CALCULATE AVERAGE COSINE OF ZENITH ANGLE FOR CURRENT COMP3 STEP
C****   AND RADIATION PERIOD
      ROT1=TWOPI*TOFDAY/24.
      ROT2=ROT1+TWOPI*DTCNDS/SDAY
      CALL COSZT (ROT1,ROT2,COSZ1)
      IF (MODRD.NE.0) GO TO 840
      ROT2=ROT1+TWOPI*NRAD*DT/SDAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
C****
C**** COMPUTE EARTH ALBEDOS AND OTHER PARAMETERS FOR BEGINNING OF DAY
C****
      JDAYR=JDAY
      JYEARR=JYEAR
      IF(JDAY.NE.JDLAST) CALL RCOMPT
      S0=S0X*S00WM2*RATLS0/RSDIST
      JDLAST=JDAY
         IHOUR=1.5+TOFDAY
C****
C**** MAIN J LOOP
C****
      DO 600 J=1,JM
      IF ((J-1)*(JM-J).NE.0) GO TO 140
C**** CONDITIONS AT THE POLES
      POLE=.TRUE.
      MODRJ=0
      IMAX=1
      GO TO 160
C**** CONDITIONS AT NON-POLAR POINTS
  140 POLE=.FALSE.
      MODRJ=MOD(J+JALTER,JRA)
      IMAX=IM
  160 CONTINUE
      JLAT=NINT(1.+(J-1.)*45./(JM-1.))
C****
C**** MAIN I LOOP
C****
      IM1=IM
      DO 500 I=1,IMAX
      MODRIJ=MODRJ+MOD(I+IALTER,IRA)
      IF (POLE) MODRIJ=0
         JR=JREG(I,J)
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      PIJ=P(I,J)
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      RANDSS=RANDU(X)
      RANDMC=RANDU(X)
         CSS=0.
         CMC=0.
         DEPTH=0.
CF       LTOP=0
      DO 240 L=1,LM
      IF(L.EQ.LS1)  PIJ=PSF-PTOP
      QSS=Q(I,J,L)/(RHSAV(I,J,L)+1.D-20)
      QL(L)=QSS
      IF(CLDSAV(I,J,L).LT.1.)
     *  QL(L)=(Q(I,J,L)-QSS*CLDSAV(I,J,L))/(1.-CLDSAV(I,J,L))
      TL(L)=T(I,J,L)*PK(L,I,J)
      IF(CLDSS(I,J,L).EQ.0.) RANDSS=RANDU(X)
      TAUSSL=0.
      TAUMCL=0.
      TAUWC(L)=0.
      TAUIC(L)=0.
      SIZEWC(L)=CSIZE(I,J,L,1)
      SIZEIC(L)=CSIZE(I,J,L,1)
         TOTCLD(L)=0.
      IF (CLDSS(I,J,L).LT.RANDSS.OR.TAUSS(I,J,L).LE.0.) GO TO 220
      TAUSSL=TAUSS(I,J,L)
      QL(L)=QSS
         CSS=1.
         AJL(J,L,28)=AJL(J,L,28)+CSS
         TOTCLD(L)=1.
CF       LTOP=L
  220 IF (CLDMC(I,J,L).LT.RANDMC.OR.TAUMC(I,J,L).LE.0.) GO TO 230
         CMC=1.
         AJL(J,L,29)=AJL(J,L,29)+CMC
         TOTCLD(L)=1.
CF       LTOP=L
         DEPTH=DEPTH+PIJ*DSIG(L)
      IF(TAUMC(I,J,L).LE.TAUSSL) GO TO 230
      TAUMCL=TAUMC(I,J,L)
      ELHX=LHE
      IF(TL(L).LE.TF) ELHX=LHS
      QL(L)=QSAT(TL(L),SIG(L)*PIJ+PTOP,ELHX)
  230    AJL(J,L,19)=AJL(J,L,19)+TOTCLD(L)
      IF(TAUSSL+TAUMCL.GT.0.) THEN
         IF(TAUMCL.GT.TAUSSL) THEN
           IF(SVLAT(I,J,L).EQ.LHE) THEN
             TAUWC(L)=TAUMCL
           ELSE
             TAUIC(L)=TAUMCL
           END IF
         ELSE
           IF(SVLHX(I,J,L).EQ.LHE) THEN
             TAUWC(L)=TAUSSL
             SIZEWC(L)=CSIZE(I,J,L,2)
           ELSE
             TAUIC(L)=TAUSSL
             SIZEIC(L)=CSIZE(I,J,L,2)
           END IF
         END IF
      END IF
  240 CONTINUE
         PIJ=P(I,J)
         AJ(J,57)=AJ(J,57)+CSS*POCEAN
         BJ(J,57)=BJ(J,57)+CSS*PLAND
         CJ(J,57)=CJ(J,57)+CSS*POICE
         DJ(JR,57)=DJ(JR,57)+CSS*DXYP(J)
         AJ(J,58)=AJ(J,58)+CMC*POCEAN
         BJ(J,58)=BJ(J,58)+CMC*PLAND
         CJ(J,58)=CJ(J,58)+CMC*POICE
         DJ(JR,58)=DJ(JR,58)+CMC*DXYP(J)
         AIJ(I,J,17)=AIJ(I,J,17)+CMC
         AJ(J,80)=AJ(J,80)+DEPTH*POCEAN
         BJ(J,80)=BJ(J,80)+DEPTH*PLAND
         CJ(J,80)=CJ(J,80)+DEPTH*POICE
         DJ(JR,80)=DJ(JR,80)+DEPTH*DXYP(J)
         CLDCV=CMC+CSS-CMC*CSS
         AJ(J,59)=AJ(J,59)+CLDCV*POCEAN
         BJ(J,59)=BJ(J,59)+CLDCV*PLAND
         CJ(J,59)=CJ(J,59)+CLDCV*POICE
         DJ(JR,59)=DJ(JR,59)+CLDCV*DXYP(J)
         AIJ(I,J,19)=AIJ(I,J,19)+CLDCV
         DO 250 L=1,LLOW
         IF (TOTCLD(L).NE.1.) GO TO 250
         AIJ(I,J,41)=AIJ(I,J,41)+1.
         GO TO 255
  250    CONTINUE
  255    DO 260 L=LMID1,LMID
         IF (TOTCLD(L).NE.1.) GO TO 260
         AIJ(I,J,42)=AIJ(I,J,42)+1.
         GO TO 265
  260    CONTINUE
  265    DO 270 L=LHI1,LHI
         IF (TOTCLD(L).NE.1.) GO TO 270
         AIJ(I,J,43)=AIJ(I,J,43)+1.
         GO TO 275
  270    CONTINUE
  275    CONTINUE
         PIJ=PSF-PTOP
         DO 280 L=LM,1,-1
         IF(L.EQ.LS1-1) PIJ=P(I,J)
         IF (TOTCLD(L).NE.1.) GO TO 280
         AIJ(I,J,18)=AIJ(I,J,18)+SIGE(L+1)*PIJ+PTOP
         GO TO 285
  280    CONTINUE
  285    DO KR=1,4
            IF (I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               IH=IHOUR
               DO INCH=1,INCHM
                  IF (IH.GT.24) IH=IH-24
                  ADAILY(IH,21,KR)=ADAILY(IH,21,KR)+TOTCLD(6)
                  ADAILY(IH,22,KR)=ADAILY(IH,22,KR)+TOTCLD(5)
                  ADAILY(IH,23,KR)=ADAILY(IH,23,KR)+TOTCLD(4)
                  ADAILY(IH,24,KR)=ADAILY(IH,24,KR)+TOTCLD(3)
                  ADAILY(IH,25,KR)=ADAILY(IH,25,KR)+TOTCLD(2)
                  ADAILY(IH,26,KR)=ADAILY(IH,26,KR)+TOTCLD(1)
                  ADAILY(IH,27,KR)=ADAILY(IH,27,KR)+CLDCV
                  ADAILY(IH,53,KR)=ADAILY(IH,53,KR)+TOTCLD(7)
                  ADAILY(IH,54,KR)=ADAILY(IH,54,KR)+TOTCLD(6)
                  ADAILY(IH,55,KR)=ADAILY(IH,55,KR)+TOTCLD(5)
                  ADAILY(IH,56,KR)=ADAILY(IH,56,KR)+TOTCLD(4)
                  ADAILY(IH,57,KR)=ADAILY(IH,57,KR)+TOTCLD(3)
                  ADAILY(IH,58,KR)=ADAILY(IH,58,KR)+TOTCLD(2)
                  ADAILY(IH,59,KR)=ADAILY(IH,59,KR)+TOTCLD(1)
                  ADAILY(IH,61,KR)=ADAILY(IH,61,KR)+CLDCV
                  IH=IH+1
               END DO
            END IF
         END DO
C****
  300 IF (MODRIJ.NE.0) GO TO 500
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
      PIJ=P(I,J)
      DO 340 L=1,LM
      IF(L.EQ.LS1) PIJ=PSF-PTOP
c      PLE(L)=SIGE(L)*PIJ+PTOP
      PLE(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TL(L)=T(I,J,L)*PK(L,I,J)     ! already defined
      IF(TL(L).LT.130..OR.TL(L).GT.370.) THEN
         WRITE(99,*) 'In Radia: TAU,I,J,L,TL',TAU,I,J,L,TL(L)
         STOP 4255
      END IF
C**** MOISTURE VARIABLES
C---- QL(L)=Q(I,J,L)        ! already defined
  340 CONTINUE
C****
C**** RADIATION, SOLAR AND THERMAL
C****
      DO 420 K=1,3
      IF(RQT(I,J,K).LT.130..OR.RQT(I,J,K).GT.370.) THEN
         WRITE(99,*) 'In Radia: TAU,I,J,L,TL',TAU,I,J,LM+K,RQT(I,J,K)
         STOP 4262
      END IF
  420 TL(LM+K)=RQT(I,J,K)
      COSZ=COSZA(I,J)
      TGO=ODATA(I,J,1)+TF
      TGOI=GDATA(I,J,3)+TF
      TGLI=GDATA(I,J,13)+TF
      TGE=GDATA(I,J,4)+TF
      TS=TSAVG(I,J)
      SNOWOI=GDATA(I,J,1)
      SNOWLI=GDATA(I,J,12)
      SNOWE=GDATA(I,J,2)
      AGESN(1)=GDATA(I,J,11)    ! land
      AGESN(2)=GDATA(I,J,9)     ! ocean ice
      AGESN(3)=GDATA(I,J,10)    ! land ice
      WEARTH=(GDATA(I,J,5)+GDATA(I,J,6))/(WFCS(I,J)+1.D-20)
      DO 430 K=1,11
  430 PVT(K)=VDATA(I,J,K)
      WS=WSAVG(I,J)
C-OLD FGOLDU(2)=XFRADJ*(1.-PEARTH)
C-OLD FGOLDU(3)=XFRADJ*PEARTH
      ILON=NINT(.5+(I-.5)*72./IM)
      CALL RCOMPX
      FSF(I,J,1)=FSRNFG(1)   !  ocean
      FSF(I,J,2)=FSRNFG(3)   !  ocean ice
      FSF(I,J,3)=FSRNFG(4)   !  land ice
      FSF(I,J,4)=FSRNFG(2)   !  soil
      IF(I.EQ.IWRITE.AND.J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      SRHR(I,J,1)=SRNFLB(1)
      TRHR(I,J,1)=STBO*(POCEAN*TGO**4+POICE*TGOI**4+PLICE*TGLI**4
     *  +PEARTH*TGE**4)-TRNFLB(1)
      DO 440 L=1,LM
      SRHR(I,J,L+1)=SRFHRL(L)
  440 TRHR(I,J,L+1)=-TRFCRL(L)
      DO 450 LR=1,3
      SRHRS(I,J,LR)=SRFHRL(LM+LR)
  450 TRHRS(I,J,LR)=-TRFCRL(LM+LR)
      DO 460 K=1,4
      SNFS(I,J,K)=SRNFLB(K+LM)
  460 TNFS(I,J,K)=TRNFLB(K+LM)-TRNFLB(1)
         TRINCG(I,J)=TRDFLB(1)
         BTMPW(I,J)=BTEMPW-TF
         ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
         ALB(I,J,2)=PLAVIS
         ALB(I,J,3)=PLANIR
         ALB(I,J,4)=ALBVIS
         ALB(I,J,5)=ALBNIR
         ALB(I,J,6)=SRRVIS
         ALB(I,J,7)=SRRNIR
         ALB(I,J,8)=SRAVIS
         ALB(I,J,9)=SRANIR
  500 IM1=I
C****
C**** END OF MAIN LOOP FOR I INDEX
C****
      IF (MODRJ.GT.IRA-2) GO TO 600
      IF (POLE) GO TO 600
C**** AVERAGING RADIATION NUMBERS AT ROW J AND EVERY OTHER COLUMN IN I
      IM1=IM-IALTER
      I=IM1+1
      IF (I.GT.IM) I=1
      IP11=2-IALTER
      DO 580 IP1=IP11,IM,2
         JR=JREG(I,J)
      SUMSR=0.
      SUMTR=0.
      DO 520 L=2,LS1-1
      SRHR(I,J,L)=P(I,J)*.5*(SRHR(IM1,J,L)/P(IM1,J)+
     *  SRHR(IP1,J,L)/P(IP1,J))
      SUMSR=SUMSR+SRHR(I,J,L)
      TRHR(I,J,L)=P(I,J)*.5*(TRHR(IM1,J,L)/P(IM1,J)+
     *  TRHR(IP1,J,L)/P(IP1,J))
  520 SUMTR=SUMTR+TRHR(I,J,L)
      DO 525 L=LS1,LM+1
      SRHR(I,J,L)=.5*(SRHR(IM1,J,L)+SRHR(IP1,J,L))
      SUMSR=SUMSR+SRHR(I,J,L)
      TRHR(I,J,L)=.5*(TRHR(IM1,J,L)+TRHR(IP1,J,L))
  525 SUMTR=SUMTR+TRHR(I,J,L)
      DO 530 LR=1,3
C     CORRECTION FOR THE RESTART OF RUN A05 W9 AT YEAR 51
      SRHRS(I,J,LR)=.5*(SRHRS(IM1,J,LR)+SRHRS(IP1,J,LR))
  530 TRHRS(I,J,LR)=.5*(TRHRS(IM1,J,LR)+TRHRS(IP1,J,LR))
         DENOM=1./(COSZ2(IM1,J)+COSZ2(IP1,J)+1.D-20)
         DO 540 K=1,9
  540    ALB(I,J,K)=(ALB(IM1,J,K)*COSZ2(IM1,J)+ALB(IP1,J,K)
     *     *COSZ2(IP1,J))*DENOM
      DTR=SUMTR+.5*(TNFS(IM1,J,1)+TNFS(IP1,J,1))
      DO 560 K=1,4
      SNFS(I,J,K)=.5*(SNFS(IM1,J,K)+SNFS(IP1,J,K))
  560 TNFS(I,J,K)=.5*(TNFS(IM1,J,K)+TNFS(IP1,J,K))-DTR
      SRHR(I,J,1)=SNFS(I,J,1)-SUMSR
      TRHR(I,J,1)=.5*(TRHR(IM1,J,1)+TRHR(IP1,J,1))
         TRINCG(I,J)=.5*(TRINCG(IM1,J)+TRINCG(IP1,J))
         BTMPW(I,J)=.5*(BTMPW(IM1,J)+BTMPW(IP1,J))
      IM1=IP1
  580 I=IM1+1
  600 CONTINUE
C****
C**** END OF MAIN LOOP FOR J INDEX
C****
      IF (JRA.LE.1) GO TO 700
C**** AVERAGING RADIATION NUMBERS AT EVERY OTHER ROW IN J
      LMT2P2=LM*2+2
      DO 620 K=1,LMT2P2
      DO 620 I=2,IM
      SRHR(I,1,K)=SRHR(1,1,K)
  620 SRHR(I,JM,K)=SRHR(1,JM,K)
      DO 640 K=1,14
      DO 640 I=2,IM
      SNFS(I,1,K)=SNFS(1,1,K)
  640 SNFS(I,JM,K)=SNFS(1,JM,K)
      J1=3-JALTER
      DO 690 J=J1,JM-1,2
      JP1=J+1
      JM1=J-1
      DO 690 I=1,IM
         JR=JREG(I,J)
      SUMSR=0.
      SUMTR=0.
      DO 660 L=2,LS1-1
      SRHR(I,J,L)=P(I,J)*.5*(SRHR(I,JP1,L)/P(I,JP1)+
     *  SRHR(I,JM1,L)/P(I,JM1))
      SUMSR=SUMSR+SRHR(I,J,L)
      TRHR(I,J,L)=P(I,J)*.5*(TRHR(I,JP1,L)/P(I,JP1)+
     *  TRHR(I,JM1,L)/P(I,JM1))
  660 SUMTR=SUMTR+TRHR(I,J,L)
      DO 664 L=LS1,LM+1
      SRHR(I,J,L)=.5*(SRHR(I,JP1,L)+SRHR(I,JM1,L))
      SUMSR=SUMSR+SRHR(I,J,L)
      TRHR(I,J,L)=.5*(TRHR(I,JP1,L)+TRHR(I,JM1,L))
  664 SUMTR=SUMTR+TRHR(I,J,L)
      DO 665 LR=1,3
      SRHRS(I,J,LR)=.5*(SRHRS(I,JP1,LR)+SRHRS(I,JM1,LR))
  665 TRHRS(I,J,LR)=.5*(TRHRS(I,JP1,LR)+TRHRS(I,JM1,LR))
         DENOM=1./(COSZ2(I,JP1)+COSZ2(I,JM1)+1.D-20)
         DO 670 K=1,9
  670    ALB(I,J,K)=(ALB(I,JP1,K)*COSZ2(I,JP1)+ALB(I,JM1,K)*
     *     COSZ2(I,JM1))*DENOM
      DTR=SUMTR+.5*(TNFS(I,JP1,1)+TNFS(I,JM1,1))
      DO 680 K=1,4
      SNFS(I,J,K)=.5*(SNFS(I,JP1,K)+SNFS(I,JM1,K))
  680 TNFS(I,J,K)=.5*(TNFS(I,JP1,K)+TNFS(I,JM1,K))-DTR
      SRHR(I,J,1)=SNFS(I,J,1)-SUMSR
      TRHR(I,J,1)=.5*(TRHR(I,JP1,1)+TRHR(I,JM1,1))
         TRINCG(I,J)=.5*(TRINCG(I,JP1)+TRINCG(I,JM1))
         BTMPW(I,J)=.5*(BTMPW(I,JP1)+BTMPW(I,JM1))
  690 CONTINUE
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
  700 CONTINUE
         DO 780 J=1,JM
         DXYPJ=DXYP(J)
         IMAX=IMAXJ(J)
         DO 720 L=1,LM
         ASRHR=0.
         ATRHR=0.
         DO 710 I=1,IMAX
         ASRHR=ASRHR+SRHR(I,J,L+1)*COSZ2(I,J)
  710    ATRHR=ATRHR+TRHR(I,J,L+1)
         AJL(J,L,9)=AJL(J,L,9)+ASRHR
  720    AJL(J,L,10)=AJL(J,L,10)+ATRHR
         ASNFS1=0.
         BSNFS1=0.
         CSNFS1=0.
         ATNFS1=0.
         BTNFS1=0.
         CTNFS1=0.
         DO 770 I=1,IMAX
         COSZ=COSZ2(I,J)
         PLAND=FLAND(I,J)
         POICE=ODATA(I,J,2)*(1.-PLAND)
         POCEAN=(1.-PLAND)-POICE
         JR=JREG(I,J)
         DO 740 LR=1,3
         ASJL(J,LR,3)=ASJL(J,LR,3)+SRHRS(I,J,LR)*COSZ
  740    ASJL(J,LR,4)=ASJL(J,LR,4)+TRHRS(I,J,LR)
         DO KR=1,4
            IF (I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               IH=IHOUR
               DO INCH=1,INCHM
                  IF (IH.GT.24) IH=IH-24
                  ADAILY(IH,2,KR)=ADAILY(IH,2,KR)+(1.-SNFS(I,J,4)/S0)
                  ADAILY(IH,3,KR)=ADAILY(IH,3,KR)+(1.-ALB(I,J,1))
                  ADAILY(IH,4,KR)=ADAILY(IH,4,KR)
     *                 +((SNFS(I,J,4)-SNFS(I,J,1))*COSZ-TNFS(I,J,4)
     *                 +TNFS(I,J,1))
                  IH=IH+1
               END DO
            END IF
         END DO
  750    CONTINUE
         AJ(J,1)=AJ(J,1)+(S0*COSZ)*POCEAN
         BJ(J,1)=BJ(J,1)+(S0*COSZ)*PLAND
         CJ(J,1)=CJ(J,1)+(S0*COSZ)*POICE
         DJ(JR,1)=DJ(JR,1)+(S0*COSZ)*DXYPJ
         AJ(J,2)=AJ(J,2)+(SNFS(I,J,4)*COSZ)*POCEAN
         BJ(J,2)=BJ(J,2)+(SNFS(I,J,4)*COSZ)*PLAND
         CJ(J,2)=CJ(J,2)+(SNFS(I,J,4)*COSZ)*POICE
         DJ(JR,2)=DJ(JR,2)+(SNFS(I,J,4)*COSZ)*DXYPJ
         ASNFS1=ASNFS1+(SNFS(I,J,1)*COSZ)*POCEAN
         BSNFS1=BSNFS1+(SNFS(I,J,1)*COSZ)*PLAND
         CSNFS1=CSNFS1+(SNFS(I,J,1)*COSZ)*POICE
         DJ(JR,3)=DJ(JR,3)+(SNFS(I,J,1)*COSZ)*DXYPJ
         AJ(J,5)=AJ(J,5)+(SRHR(I,J,1)*COSZ/(ALB(I,J,1)+1.D-20))*POCEAN
         BJ(J,5)=BJ(J,5)+(SRHR(I,J,1)*COSZ/(ALB(I,J,1)+1.D-20))*PLAND
         CJ(J,5)=CJ(J,5)+(SRHR(I,J,1)*COSZ/(ALB(I,J,1)+1.D-20))*POICE
         DJ(JR,5)=DJ(JR,5)+(SRHR(I,J,1)*COSZ/(ALB(I,J,1)+1.D-20))*DXYPJ
         AJ(J,6)=AJ(J,6)+(FSF(I,J,1)*COSZ)*POCEAN
         SRNFLG=FSF(I,J,3)*FLICE(I,J)+FSF(I,J,4)*(PLAND-FLICE(I,J))
         BJ(J,6)=BJ(J,6)+(SRNFLG*COSZ)
         CJ(J,6)=CJ(J,6)+(FSF(I,J,2)*COSZ)*POICE
         DJ(JR,6)=DJ(JR,6)+(SRHR(I,J,1)*COSZ)*DXYPJ
         AJ(J,55)=AJ(J,55)+BTMPW(I,J)*POCEAN
         BJ(J,55)=BJ(J,55)+BTMPW(I,J)*PLAND
         CJ(J,55)=CJ(J,55)+BTMPW(I,J)*POICE
         DJ(JR,55)=DJ(JR,55)+BTMPW(I,J)*DXYPJ
         AJ(J,67)=AJ(J,67)+TRINCG(I,J)*POCEAN
         BJ(J,67)=BJ(J,67)+TRINCG(I,J)*PLAND
         CJ(J,67)=CJ(J,67)+TRINCG(I,J)*POICE
         DJ(JR,67)=DJ(JR,67)+TRINCG(I,J)*DXYPJ
         AJ(J,70)=AJ(J,70)-TNFS(I,J,4)*POCEAN
         BJ(J,70)=BJ(J,70)-TNFS(I,J,4)*PLAND
         CJ(J,70)=CJ(J,70)-TNFS(I,J,4)*POICE
         DJ(JR,70)=DJ(JR,70)-TNFS(I,J,4)*DXYPJ
         ATNFS1=ATNFS1-TNFS(I,J,1)*POCEAN
         BTNFS1=BTNFS1-TNFS(I,J,1)*PLAND
         CTNFS1=CTNFS1-TNFS(I,J,1)*POICE
         DJ(JR,71)=DJ(JR,71)-TNFS(I,J,1)*DXYPJ
         DO 760 K=2,9
         AJ(J,K+70)=AJ(J,K+70)+(S0*COSZ)*ALB(I,J,K)*POCEAN
         BJ(J,K+70)=BJ(J,K+70)+(S0*COSZ)*ALB(I,J,K)*PLAND
         CJ(J,K+70)=CJ(J,K+70)+(S0*COSZ)*ALB(I,J,K)*POICE
  760    DJ(JR,K+70)=DJ(JR,K+70)+(S0*COSZ)*ALB(I,J,K)*DXYPJ
         AIJ(I,J,21)=AIJ(I,J,21)-TNFS(I,J,4)
         AIJ(I,J,24)=AIJ(I,J,24)+(SNFS(I,J,4)*COSZ)
         AIJ(I,J,25)=AIJ(I,J,25)+(S0*COSZ)
         AIJ(I,J,26)=AIJ(I,J,26)+(SRHR(I,J,1)*COSZ)
         AIJ(I,J,27)=AIJ(I,J,27)+(SRHR(I,J,1)*COSZ/(ALB(I,J,1)+1.D-20))
         AIJ(I,J,44)=AIJ(I,J,44)+BTMPW(I,J)
         AIJ(I,J,45)=AIJ(I,J,45)+S0*COSZ*ALB(I,J,2)
  770    CONTINUE
         AJ(J,3)=AJ(J,3)+ASNFS1
         BJ(J,3)=BJ(J,3)+BSNFS1
         CJ(J,3)=CJ(J,3)+CSNFS1
         AJ(J,71)=AJ(J,71)+ATNFS1
         BJ(J,71)=BJ(J,71)+BTNFS1
         CJ(J,71)=CJ(J,71)+CTNFS1
  780    CONTINUE
         DO 790 L=1,LM
         DO 790 I=1,IM
         AIL(I,L,7)=AIL(I,L,7)+((SRHR(I,JEQ-2,L+1)*COSZ2(I,JEQ-2)+
     *     TRHR(I,JEQ-2,L+1))*DXYP(JEQ-2)+(SRHR(I,JEQ-1,L+1)*
     *     COSZ2(I,JEQ-1)+TRHR(I,JEQ-1,L+1))*DXYP(JEQ-1)+
     *     (SRHR(I,JEQ,L+1)*COSZ2(I,JEQ)+TRHR(I,JEQ,L+1))*DXYP(JEQ))
         AIL(I,L,11)=AIL(I,L,11)+(SRHR(I,J50N,L+1)*COSZ2(I,J50N)+
     *     TRHR(I,J50N,L+1))*DXYP(J50N)
  790    AIL(I,L,15)=AIL(I,L,15)+(SRHR(I,J70N,L+1)*COSZ2(I,J70N)+
     *     TRHR(I,J70N,L+1))*DXYP(J70N)
C****
C**** UPDATE THE TEMPERATURES BY RADIATION
C****
      DO LR=1,3
         DO J=1,JM
            IMAX=IMAXJ(J)
            DO I=1,IMAX
               RQT(I,J,LR)=RQT(I,J,LR)+(SRHRS(I,J,LR)*COSZ2(I,J)
     *              +TRHRS(I,J,LR))*COE(LR+LM)
            END DO
         END DO
      END DO
  840 DO L=1,LS1-1
         DO J=1,JM
            IMAX=IMAXJ(J)
            DO I=1,IMAX
               T(I,J,L)=T(I,J,L)+(SRHR(I,J,L+1)*COSZ1(I,J)+TRHR(I,J,L+1)
     *              )*COE(L)/(P(I,J)*PK(L,I,J))
            END DO
         END DO
      END DO
      DO L=LS1,LM
         DO J=1,JM
            IMAX=IMAXJ(J)
            DO I=1,IMAX
               T(I,J,L)=T(I,J,L)+(SRHR(I,J,L+1)*COSZ1(I,J)+TRHR(I,J,L+1)
     *              )*COE(L)/((PSF-PTOP)*PK(L,I,J))
            END DO
         END DO
      END DO
C**** daily diagnostics
      DO KR=1,4
         ADAILY(IHOUR,1,KR)=ADAILY(IHOUR,1,KR)+S0*COSZ1(IJD6(1,KR)
     *        ,IJD6(2,KR))
      END DO

      RETURN
      END
      SUBROUTINE GROUND
C****
C**** THIS SUBROUTINE USES THE SURFACE FLUXES TO PREDICT IN TIME THE
C**** GROUND TEMPERATURE, GROUND WATER AND ICE, AND SNOW MELTING.
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,rhow,rhoi,shv,shw,shi
      USE E001M12_COM
      USE GEOM
      USE PBLCOM, only : tsavg
      USE DAGCOM, only : aj,bj,cj,dj,aij,jreg
      USE OCEAN, only : odata,XSI1,XSI2,XSI3,XSI4,R1,R2,R3,R4,TTRUNC,Z1I
     *     ,Z2OIM,ACE1I,AC2OIM,OTA,OTB,OTC,TFO,T50

      IMPLICIT REAL*8 (A-H,O-Z)
C*
      PARAMETER (ALPHA = 1.0, dSNdML = 0.,FLEAD = 0.10d0)
C*
      REAL*8 MELT1, MELT4, MSI1, MSI2
C*
      COMMON/WORK3/E0(IM,JM,4),E1(IM,JM,4),EVAPOR(IM,JM,4)
            COMMON/oldDAG/GDEEP(IM,JM,3)

      DATA ALAMI/2.1762/,Z2LI/2.9/,Z1E/.1/,Z2E/4./,Z2OIX/4.9/

      DATA IFIRST/1/
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** ODATA  1  OCEAN TEMPERATURE (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        2  EARTH SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)
C****        8  currently not used
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)
C****       15  OCEAN ICE TEMPERATURE OF THIRD LAYER (C)
C****       16  OCEAN ICE TEMPERATURE OF FOURTH LAYER (C)
C****
C**** VDATA  9  WATER FIELD CAPACITY OF FIRST LAYER (KG/M**2)
C****       10  WATER FIELD CAPACITY OF SECOND LAYER (KG/M**2)
C****
C*
C***  Make 50 day mean of surface air temperature
C***    for lake ice extent calculation
C*
      FACT_T50 = 1.-1./(24.*50.)
      FACT_TSAVG = 1./(24.*50.)
      DO J = 1,JM
        DO I = 1,IM
          T50(I,J) = T50(I,J)*FACT_T50
     *             + (TSAVG(I,J)-273.16)*FACT_TSAVG
        END DO
      END DO
C*
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
      IF (KOCEAN.NE.1) GO TO 10
      READ (12) OTA,OTB,OTC
C     CALL DREAD (17,OTA,IM*JM*3,OTA)
      REWIND 12
   10 DTSRCE=NDYN*DT
c      ACE1I=Z1I*RHOI
c      AC2OIM=Z2OIM*RHOI
C*
      YSI1 = XSI1*ACE1I/(ACE1I+AC2OIM)
      YSI2 = XSI2*ACE1I/(ACE1I+AC2OIM)
      YSI3 = XSI3*AC2OIM/(ACE1I+AC2OIM)
      YSI4 = XSI4*AC2OIM/(ACE1I+AC2OIM)
C*
      ATRUNC=0.
      BYZICX=1./(Z1I+Z2OIX)
      HC1I=ACE1I*SHI
      HC2LI=Z2LI*RHOI*SHI
      HC1DE=Z1E*1129950.
      HC2DE=Z2E*1129950.+3.5*.125*RHOW*3100.
      DIFFUS=DTSRCE/SDAY
   50 ANGLE=TWOPI*JDAY/365.
      SINANG=SIN(ANGLE)
      SN2ANG=SIN(2*ANGLE)
      SN3ANG=SIN(3*ANGLE)
      SN4ANG=SIN(4*ANGLE)
      COSANG=COS(ANGLE)
      CS2ANG=COS(2*ANGLE)
      CS3ANG=COS(3*ANGLE)
      CS4ANG=COS(4*ANGLE)
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO 980 J=1,JM
      IMAX=IMAXJ(J)
         BF1DT=0.
         CF1DT=0.
         AOTDT=0.
         COTDT=0.
         AEFO=0.
         CEFI=0.
         BEDIFS=0.
         CEDIFS=0.
         BERUN0=0.
         CF2DT=0.
         BERUN2=0.
         CERUN2=0.
         AERUN4=0.
         CERUN4=0.
         ATG1=0.
         BTG1=0.
         CTG1=0.
         ATG2=0.
         BTG2=0.
         CTG2=0.
         ATG3=0.
         AEVAP=0.
         BEVAP=0.
         CEVAP=0.
         BDIFS=0.
         CDIFS=0.
         AIFO=0.
         CIFI=0.
         BRUN0=0.
         CRUN0=0.
         BRUN2=0.
         CRUN2=0.
         ARUN4=0.
         CRUN4=0.
         BWTR1=0.
         BACE1=0.
         BWTR2=0.
         BACE2=0.
         CACE2=0.
         BSNOW=0.
         CSNOW=0.
         CICOV=0.
      DO 960 I=1,IMAX
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      ROICE=ODATA(I,J,2)
      POICE=ROICE*PWATER
      POCEAN=PWATER-POICE
         JR=JREG(I,J)
         DXYPJ=DXYP(J)
         SNOWS=0.
         WTR1S=0.
         ACE1S=0.
         WTR2S=0.
         ACE2S=0.
         TG1S=0.
         TG2S=0.
         EVAPS=0.
         RUN0S=0.
         DIFSS=0.
C****
      IF (PWATER.LE.0.) GO TO 400
C****
C**** OCEAN
C****
      EVAP=EVAPOR(I,J,1)
         ATG1=ATG1+ODATA(I,J,1)*POCEAN
         TG1S=TG1S+ODATA(I,J,1)*POCEAN
         AEVAP=AEVAP+EVAP*POCEAN
         EVAPS=EVAPS+EVAP*POCEAN
         AIJ(I,J,57)=AIJ(I,J,57)+ODATA(I,J,1)
         AIJ(I,J,61)=AIJ(I,J,61)+EVAP*POCEAN
      IF (KOCEAN.EQ.1) GO TO 60
         ATG2=ATG2+ODATA(I,J,1)*POCEAN
         TG2S=TG2S+ODATA(I,J,1)*POCEAN
      IF (POICE.GT.0.) GO TO 110
      GO TO 400
   60 TGW=ODATA(I,J,1)
      WTRO=Z1O(I,J)*RHOW
      ENRGO0=WTRO*TGW*SHW
      EOFRZ=WTRO*TFO*SHW
      F0DT=E0(I,J,1)
         AIJ(I,J,65)=AIJ(I,J,65)+F0DT*POCEAN
      OTDT=DTSRCE*(OTA(I,J,4)*SN4ANG+OTB(I,J,4)*CS4ANG
     *            +OTA(I,J,3)*SN3ANG+OTB(I,J,3)*CS3ANG
     *            +OTA(I,J,2)*SN2ANG+OTB(I,J,2)*CS2ANG
     *            +OTA(I,J,1)*SINANG+OTB(I,J,1)*COSANG+OTC(I,J))
         ATG2=ATG2+ODATA(I,J,4)*POCEAN
         TG2S=TG2S+ODATA(I,J,4)*POCEAN
         ATG3=ATG3+ODATA(I,J,5)*POCEAN
         AOTDT=AOTDT+OTDT*POCEAN
      RUN4=-EVAP
      ERUN4=RUN4*TGW*SHW
         AERUN4=AERUN4+ERUN4*POCEAN
         ARUN4=ARUN4+RUN4*POCEAN
      ENRGO=F0DT+OTDT-ERUN4
      IF (ENRGO0+ENRGO.LT.EOFRZ) GO TO 80
C**** FLUXES RECOMPUTE TGO WHICH IS ABOVE FREEZING POINT FOR OCEAN
      ENRGFO=0.
      ACEFO=0.
      IF (ROICE.GT.0.) GO TO 100
      ODATA(I,J,1)=TGW+(ENRGO/(WTRO*SHW)+TTRUNC)
      GO TO 400
C**** FLUXES COOL TGO TO FREEZING POINT FOR OCEAN AND FORM SOME ICE
c  80 ACEFO=(ENRGO0+ENRGO-EOFRZ)/(TFO*(SHI-SHW)-LHM)
   80 ACEFO=(ENRGO0+ENRGO-EOFRZ)/(-LHM)
      ENRGFO=ACEFO*(TFO*SHI-LHM)
         AEFO=AEFO-ENRGFO*POCEAN
         AIFO=AIFO-ACEFO*POCEAN
      IF (ROICE.GT.0.) GO TO 100
      ROICE=ACEFO/(ACE1I+AC2OIM)
      ODATA(I,J,1)=TFO
      ODATA(I,J,2)=ROICE
      GDATA(I,J,1)=0.
      GDATA(I,J,3)=TFO
      GDATA(I,J,7)=TFO
      GDATA(I,J,15) = TFO
      GDATA(I,J,16) = TFO
      ODATA(I,J,3)=AC2OIM
      GO TO 400
C****
  100 ACE2F=0.
      ACE2M=0.
C****
C**** OCEAN ICE
C****
  110 SNOW = GDATA(I,J,1) ! snow mass (kg/m^2)
      TG1 = GDATA(I,J,3) ! first layer sea ice temperature
      TG2 = GDATA(I,J,7) ! second layer sea ice temperature
      TG3 = GDATA(I,J,15) ! third layer sea ice temperature
      TG4 = GDATA(I,J,16) ! fourth layer sea ice temperature
C*
      MSI2=ODATA(I,J,3) ! second (physical) layer ice mass
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C*
C***  CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
C*
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1 ! J/m^2
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1 ! J/m^2
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2 ! J/m^2
      HSI4 = (SHI*TG4-LHM)*XSI4*MSI2 ! J/m^2
C*
         AIJ(I,J,1)=AIJ(I,J,1)+POICE
         AIJ(I,J,58)=AIJ(I,J,58)+MSI2*POICE
      F0DT=E0(I,J,2)
         AIJ(I,J,66)=AIJ(I,J,66)+F0DT*POICE
      F1DT=E1(I,J,2)
      EVAP=EVAPOR(I,J,2)
         AIJ(I,J,62)=AIJ(I,J,62)+EVAP*POICE
      Z2=MSI2/RHOI
         IF (KOCEAN.NE.1) GO TO 120
      WTRI0=WTRO-(SNOW+ACE1I+MSI2)
      EIW0=WTRI0*TGW*SHW
      WTRW0=WTRO-ROICE*(SNOW+ACE1I+MSI2)
      ENRGW0=WTRW0*TGW*SHW
      RUN0=0.
         DIFSI=0.
         EDIFSI=0.
      RUN4=-EVAP
C*
      WTRW0 = WTRW0-ROICE*RUN4 ! water mass "+-" dew/evaporation
C*
      ERUN4=TGW*RUN4*SHW
         CERUN4=CERUN4+ERUN4*POICE
         CRUN4=CRUN4+RUN4*POICE
C****
C**** OCEAN ICE, CALCULATE TG1
C****
  120 CONTINUE
C*
      HC2 = SHI*XSI2*MSI1 ! heat capacity of ice layer 2 (J/(degC*m^2))
      HC3 = SHI*XSI3*MSI2 ! heat capacity of ice layer 3 (J/(degC*m^2))
      HC4 = SHI*XSI4*MSI2 ! heat capacity of ice layer 4 (J/(degC*m^2))
C*
C***  CALCULATE AND APPLY DIFFUSIVE AND SURFACE ENERGY FLUXES
C*
      dF2dTI = ALAMI*RHOI*DTSRCE/(0.5*XSI2*MSI1+0.5*XSI3*MSI2)
C***           temperature derivative from F2 diffusive flux
      dF3dTI = ALAMI*RHOI*DTSRCE*2./MSI2
C***           temperature derivative from F3 diffusive flux
      dF4dTI = ALAMI*RHOI*DTSRCE*2.*R4/MSI2
C***           temperature derivative from F4 diffusive flux
C*
CEXP  F2 = dF2dTI*(TG2-TG3) ! the diffusive
CEXP  F3 = dF3dTI*(TG3-TG4) ! fluxes from
CEXP  F4 = dF4dTI*(TG4-TGW) ! explicit method
C*
C***  DIFFUSIVE FLUXES FROM IMPLICIT METHOD
C*
      F2 = dF2dTI*(HC2*(TG2-TG3)+ALPHA*E1(I,J,2))/
     A     (HC2+ALPHA*dF2dTI)
      F3 = dF3dTI*(HC3*(TG3-TG4)+ALPHA*F2)/
     A     (HC3+ALPHA*dF3dTI)
      E1(I,J,1) = dF4dTI*(HC4*(TG4-TGW)+ALPHA*F3)/
     A            (HC4+ALPHA*dF4dTI)
      F2DT = E1(I,J,1)
C*
      HSI1 = HSI1+(E0(I,J,2)-E1(I,J,2))
      HSI2 = HSI2+(E1(I,J,2)-F2)
      HSI3 = HSI3+(F2-F3)
      HSI4 = HSI4+(F3-E1(I,J,1))
C*
      DEW = -EVAP ! dew to the surface
C*
      IF (HSI1/LHM+XSI1*MSI1+DEW .LE. 0.) GO TO 160 ! go to freezing
C*
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
C*
      MELT1 = HSI1/LHM+XSI1*MSI1+DEW ! melting + dew to layer 1
      RUN0 = MELT1 ! water mass that flows to the ocean (kg/m^2)
      IF (KOCEAN.EQ.1) WTRW0 = WTRW0+ROICE*RUN0 ! ocean mass (kg/m^2)
         CRUN0=CRUN0+RUN0*POICE
         RUN0S=RUN0S+RUN0*POICE
C*
C***  EVAPORATION FROM THE SURFACE
C*
      IF (DEW .GT. 0.) GO TO 140 ! go to the dew case
C*
C***  DEW IS EVAPORATION NOW
C***  EVAPORATION REDUCES SNOW OR ICE, MELT1>0.
C*
      IF (MSI1-ACE1I+DEW .GT. MELT1) GO TO 130
C*
C***  ALL SNOW AND SOME ICE EVAPORATE AND MELT
C***  ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
C*
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 2 into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 3 into layer 2
C*
      FHSI1 = HSI2*FMSI1*R2/MSI1 ! energy of ice mass FMSI1
      FHSI2 = HSI3*FMSI2*R3/MSI2 ! energy of ice mass FMSI2
C*
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
C*
      SNOW = 0. ! all snow evaporates and melts
      DIFS = FMSI2 ! <0 (upward mass flux into layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN.EQ.1) GO TO 170
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         CERUN2=CERUN2+EDIFS*POICE
         CRUN2=CRUN2+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
      GO TO 370
C*
  130 CONTINUE
C*
C***  EVAPORATION AND MELTING REDUCE JUST SNOW AMOUNT
C***  ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
C*
      FMSI2 = MIN(dSNdML*MELT1, MSI1-ACE1I+DEW-MELT1) ! > 0.
      FMSI1 = XSI1*FMSI2+XSI2*(DEW-MELT1)
C*
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward from layer 2 into layer 3
C*
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
C*
      MSI1 = MSI1+DEW-MELT1-FMSI2 ! kg/m^2
      SNOW = MSI1-ACE1I
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         CERUN2=CERUN2+EDIFS*POICE
         CRUN2=CRUN2+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
      GO TO 370
C*
  140 CONTINUE
C*
C***  DEW INCREASES ICE AMOUNT,  MELT1 > 0.
C*
      IF (MSI1-ACE1I .GT. MELT1) GO TO 150
C*
C***  ALL SNOW AND SOME ICE MELT
C*
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW-MELT1 ! upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1
      SNOW = 0.
C*
      IF (FMSI2 .LE. 0.) THEN ! (if melting is greater than dew)
C*
C*** ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
C*
        FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward into layer 1 from layer 2
        FHSI2 = HSI3*FMSI2*R3/MSI2 ! upward into layer 2 from layer 3
C*
        HSI1 = HSI1-FHSI1
        HSI2 = HSI2+(FHSI1-FHSI2)
C*
        DIFS = FMSI2 ! <0 upward ice mass into layer 2
        EDIFS = FHSI2 ! energy of diffused ice mass DIFS
        IF (KOCEAN .EQ. 1) GO TO 170
          CEDIFS=CEDIFS+EDIFS*POICE
          CDIFS=CDIFS+DIFS*POICE
          CERUN2=CERUN2+EDIFS*POICE
          CRUN2=CRUN2+DIFS*POICE
          DIFSS=DIFSS+DIFS*POICE
        GO TO 370
      ENDIF
C*
C***  ICE ADVECTION IS DOWNWARD INTO LAYER 3 FROM LAYER 2
C***  (if dew is greater than melting)
C*
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward into layer 3 from layer 2
C*
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
C*
      MSI1 = ACE1I
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
          CEDIFS=CEDIFS+EDIFS*POICE
          CDIFS=CDIFS+DIFS*POICE
          CERUN2=CERUN2+EDIFS*POICE
          CRUN2=CRUN2+DIFS*POICE
          DIFSS=DIFSS+DIFS*POICE
      GO TO 370
C*
  150 CONTINUE
C*
C***  MELTING REDUCES SNOW AMOUNT
C***  ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
C*
      CMPRS = MIN(dSNdML*MELT1, MSI1-ACE1I-MELT1) ! > 0.
      FMSI1 = DEW+XSI1*CMPRS-XSI2*MELT1
      FMSI2 = DEW+CMPRS ! > 0. downward into layer 3 from layer 2
C*
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 3
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward energy flux into layer 3
C*
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
C*
      MSI1 = MSI1-MELT1-CMPRS ! new mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass after melting and compression
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         CERUN2=CERUN2+EDIFS*POICE
         CRUN2=CRUN2+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
      GO TO 370
C*
  160 CONTINUE
C*
C***  THE FIRST LAYER IS BELOW FREEZING
C***  NO SNOW OR ICE MELTS IN LAYER 1
C*
      IF (DEW .GT. 0.) GO TO 200 ! go to the dew case
      IF (MSI1-ACE1I+DEW .GE. 0.) GO TO 190
C*
C***  ALL SNOW AND SOME ICE EVAPORATE
C***  ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
C*
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW ! < 0. upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW ! < 0. upward into layer 2 from layer 3
C*
      FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward energy flux into layer 1
      FHSI2 = HSI3*FMSI2*R3/MSI2 ! upward energy flux into layer 2
C*
      HSI1 = HSI1-FMSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
C*
      SNOW = 0. ! all snow evaporated
      DIFS = FMSI2 ! <0 upward ice mass into layer 2
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN.EQ.1) GO TO 170
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         CERUN2=CERUN2+EDIFS*POICE
         CRUN2=CRUN2+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
      GO TO 370
C*
  170 CONTINUE
C*
      MSI1 = ACE1I
         DIFSI = ROICE*DIFS
         EDIFSI = ROICE*EDIFS
C*
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 180 ! go to freezing case
C*
C***  FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
C*
      MELT4 = HSI4/LHM+XSI4*MSI2 ! > 0. melted ice from layer 4
C*
      ACE2M = MELT4
      WTRW0 = WTRW0+ROICE*ACE2M
         AIFI=AIFI+ACE2M*POICE
C*
      FMSI3 = XSI4*FMSI2+XSI3*MELT4
      IF (FMSI3 .LE. 0.) FHSI3 = -LHM*FMSI3 ! upward into layer 3
      IF (FMSI3 .GT. 0.) FHSI3 = HSI3*FMSI3*R3/MSI2 ! downward
C*
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
C*
      MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      GO TO 230
C*
  180 CONTINUE
C*
C***  NO ICE MELTS IN LAYER 4
C*
C     FMSI3 = XSI4*FMSI2 ! upward mass flux into layer 3 from layer 4
C*
      FHSI3 = HSI4*FMSI2/MSI2 ! upward energy flux into layer 3
C*
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
C*
      MSI2 = MSI2+FMSI2
      GO TO 230
C*
  190 CONTINUE
C*
C***  JUST SOME SNOW EVAPORATES
C***  NO ADVECTION BETWEEN LAYER 2 AND 3
C*
C     FMSI1 = XSI2*DEW ! < 0. upward mass flux into layer 1
C     FMSI2 = 0. ! no ice advection between layers 2 and 3
C*
      FHSI1 = HSI2*DEW/MSI1 ! upward energy flux into layer 1
C*
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+FHSI1
C*
      MSI1 = MSI1+DEW ! new ice mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass
C*
      IF (KOCEAN .NE. 1) GO TO 370
C*
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 230 ! go to freezing case
C*
C***  FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
C*
      MELT4 = HSI4/LHM+XSI4*MSI2 ! melted ice from layer 4
C*
      ACE2M = MELT4
      WTRW0 = WTRW0+ROICE*ACE2M
         AIFI=AIFI+ACE2M*POICE
C*
C     FMSI3 = XSI3*MELT4 ! > 0. downward mass flux into layer 4
C*
      FHSI3 = HSI3*MELT4/MSI2 ! downward heat flux into layer 4
C*
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+FHSI3
C*
      MSI2 = MSI2-MELT4 ! new ice mass of the second physical layer
      GO TO 230
C*
  200 CONTINUE
C*
C***  DEW INCREASES ICE AMOUNT, ADVECT ICE DOWNWARD
C***  ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
C*
C     FMSI1 = DEW ! > 0. downward mass flux into layer 2
      FMSI2 = DEW ! > 0. downward mass flux into layer 3
C*
      FHSI1 = HSI1*FMSI2/(XSI1*MSI1+DEW) ! downward heat flux
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux into layer 3
C*
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
C*
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         CERUN2=CERUN2+EDIFS*POICE
         CRUN2=CRUN2+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
      GO TO 370
C*
  210 CONTINUE
C*
         DIFSI = ROICE*DIFS
         EDIFSI = ROICE*EDIFS
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 220 ! go to freezing case
C*
C**** FLUXES HEAT UP TG2 TO FREEZING POINT AND MELT SOME ICE
C*
      MELT4 = HSI4/LHM+XSI4*MSI2 ! melted ice from layer 4
C*
      ACE2M = MELT4
      WTRW0 = WTRW0+ROICE*ACE2M
         AIFI=AIFI+ACE2M*POICE
C*
C     FMSI3 = XSI4*FMSI2+XSI3*MELT4 ! > 0. downward into layer 4
C*
      FHSI3 = HSI3*((XSI4/XSI3)*FMSI2+MELT4)/MSI2 ! downward
C*
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
C*
      MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      GO TO 230
C*
  220 CONTINUE
C*
C***  NO ICE MELTS IN LAYER 4
C*
C     FMSI3 = XSI4*FMSI2 ! > 0. downward mass flux into layer 4
C*
      FHSI3 = HSI3*(XSI4/XSI3)*FMSI2/MSI2 ! downward heat flux
C*
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
C*
      MSI2 = MSI2+FMSI2
C*
  230  CONTINUE
C*
C***  CALCULATE THE ENERGY OF THE WATER BELOW THE ICE
C***          AT FREEZING POINT AND
C***  CHECK WHETHER NEW ICE MUST BE FORMED
C*
      ENRGIW = E1(I,J,1)+OTDT-ERUN4 ! heat flux to the ocean under ice
      ENRGFI = 0.
C*
      WTRI1 = WTRO-(MSI1+MSI2) ! new mass of ocean (kg/m^2)
C*
      EFIW = WTRI1*TFO*SHW ! freezing energy of ocean mass WTRI1
      IF (EIW0+ENRGIW .GT. EFIW) GO TO 250 ! go to no freezing case
C*
C***  FLUXES WOULD COOL TGW TO FREEZING POINT
C***  AND FREEZE SOME MORE ICE
C*
c      ACE2F = (EIW0+ENRGIW-EFIW)/(TFO*(SHI-SHW)-LHM) !
      ACE2F = (EIW0+ENRGIW-EFIW)/(-LHM) !
C*            ocean mass that freezes under the ice
      ENRGFI = ACE2F*(TFO*SHI-LHM) ! energy of frozen ice
         CEFI=CEFI-ENRGFI*POICE
         CIFI=CIFI-ACE2F*POICE
C*
C***  CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
C*
C     FMSI3 = -XSI3*ACE2F ! < 0.
C     FMSI4 = -ACE2F
C*
      FHSI3 = -HSI4*ACE2F*(XSI3/XSI4)/MSI2
C*
C***  COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
C*
      IF (ACEFO .GT. 0.) GO TO 240
C*
C***  NEW ICE IS FORMED BELOW OLD SEA ICE
C*
      WTRW = WTRW0-ROICE*ACE2F ! new ocean mass
      ENRGW = ENRGW0+ROICE*(ENRGIW-ENRGFI)+(1.-ROICE)*ENRGO ! energy
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
c     TGW = TFO
C*
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+(FHSI3+ENRGFI)
C*
      MSI2 = MSI2+ACE2F ! new ice mass of physical layer 2
      GO TO 270
C*
  240 CONTINUE
C*
C***  NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
C*
      WTRW = WTRW0-(1.-ROICE)*ACEFO-ROICE*ACE2F ! new ocean mass
      ENRGW = ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
c     TGW = TFO
C*
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*(MSI2+ACE2F))/(ROICE+DRSI) ! layer 2
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
C*
      HSI1 = ((1.-ROICE)*ENRGFO*YSI1+ROICE*HSI1)/(ROICE+DRSI)
      HSI2 = ((1.-ROICE)*ENRGFO*YSI2+ROICE*HSI2)/(ROICE+DRSI)
      HSI3 = ((1.-ROICE)*ENRGFO*YSI3+ROICE*(HSI3-FHSI3))/
     A       (ROICE+DRSI)
      HSI4 = ((1.-ROICE)*ENRGFO*YSI4+ROICE*(HSI4+FHSI3+ENRGFI))/
     A       (ROICE+DRSI)
C*
      ROICE = ROICE+DRSI ! new ice concentration
      GO TO 270
C*
  250 CONTINUE
C*
      IF (ACEFO .GT. 0.) GO TO 260 ! new ice on the open ocean
C*
C***  NO NEW ICE IS FORMED UNDERNEATH THE OLD ONE
C*
      WTRW = WTRW0 ! new ocean mass
      ENRGW = ENRGW0+ROICE*ENRGIW+(1.-ROICE)*ENRGO ! energy of new oc.
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
c     TGW = TFO
      GO TO 270
C*
  260 CONTINUE
C*
C***  NEW ICE IS FORMED ON THE OPEN OCEAN
C*
      WTRW = WTRW0-(1.-ROICE)*ACEFO ! new ocean mass
      ENRGW = ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*ENRGIW
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
c     TGW = TFO
C*
      DRSI = (1.-ROICE)*ACEFO/(ACE1I+AC2OIM) ! new ice on the open oc.
      MSI1 = (DRSI*ACE1I+ROICE*MSI1)/(ROICE+DRSI) ! mass of layer 1
      MSI2 = (DRSI*AC2OIM+ROICE*MSI2)/(ROICE+DRSI) ! layer 2
      SNOW = SNOW*ROICE/(ROICE+DRSI) ! redistributed over old and new
C*
      HSI1 = ((1.-ROICE)*ENRGFO*YSI1+ROICE*HSI1)/(ROICE+DRSI)
      HSI2 = ((1.-ROICE)*ENRGFO*YSI2+ROICE*HSI2)/(ROICE+DRSI)
      HSI3 = ((1.-ROICE)*ENRGFO*YSI3+ROICE*HSI3)/(ROICE+DRSI)
      HSI4 = ((1.-ROICE)*ENRGFO*YSI4+ROICE*HSI4)/(ROICE+DRSI)
C*
      ROICE = ROICE+DRSI ! new ice concentration
C*
  270 CONTINUE
C*
      IF (FLAKE(I,J) .GT. 0.) GO TO 360 ! no compression for lake ice
C*
C***  COMPRESS THE ICE HORIZONTALLY
C*
      IF (MSI2 .GE. AC2OIM) GO TO 280 ! ice is thick enough
C*
C***  SEA ICE IS TOO THIN
C*
      ROICEN = ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM) ! new ice concentr.
      GO TO 290
C*
  280 CONTINUE
C*
      OPNOCN=.06*(RHOI/(ACE1I+MSI2)-BYZICX)
      IF ((1.-ROICE) .GE. OPNOCN) GO TO 360
      ROICEN = 1.-OPNOCN
C*
 290  CONTINUE
C*
      DRSI = ROICEN-ROICE ! < 0. compressed ice concentration
      DRI = ROICE-ROICEN ! > 0. for diagnostics
C*
C     FMSI3 = XSI3*FMSI4 ! < 0. upward ice mass flux into layer 3
      FMSI4 = (MSI1+MSI2)*(DRSI/ROICEN) ! upward ice mass into layer 4
C*
      FHSI3 = HSI4*FMSI4*(XSI3/XSI4)/MSI2 ! upward heat flux
      FHSI4 = (HSI1+HSI2+HSI3+HSI4)*(DRSI/ROICEN)
C*
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+(FHSI3-FHSI4)
C*
      MSI2 = MSI2-FMSI4 ! new ice mass of second physical layer
      SNOW = SNOW*(ROICE/ROICEN) ! new snow mass
C*
      DIFS = DRI*ACE1I/ROICE
         EDIFSI = EDIFSI+ROICE*(HSI1+HSI2)*(DIFS/MSI1)*0.5
         DIFSI = DIFSI+ROICE*DIFS
C*
      ROICE = ROICEN
C*
C**** RESAVE PROGNOSTIC QUANTITIES
C*
  360 CONTINUE
C*
      ODATA(I,J,1)=TGW
      ODATA(I,J,2)=ROICE
      ODATA(I,J,3)=MSI2
C*
         COTDT=COTDT+OTDT*POICE
         CEDIFS=CEDIFS+EDIFSI*PWATER
         CDIFS=CDIFS+DIFSI*PWATER
         DIFSS=DIFSS+DIFSI*PWATER
C*
  370 CONTINUE
C*
      TG1 = (HSI1/(XSI1*MSI1) +LHM)/SHI ! temperature of layer 1
      TG2 = (HSI2/(XSI2*MSI1) +LHM)/SHI ! temperature of layer 2
      TG3 = (HSI3/(XSI3*MSI2) +LHM)/SHI ! temperature of layer 3
      TG4 = (HSI4/(XSI4*MSI2) +LHM)/SHI ! temperature of layer 4
C*
      GDATA(I,J,1) = SNOW
      GDATA(I,J,3) = TG1
      GDATA(I,J,7) = TG2
      GDATA(I,J,15) = TG3
      GDATA(I,J,16) = TG4
         CSNOW=CSNOW+SNOW*POICE
         CTG1=CTG1+TG1*POICE
         CTG2=CTG2+TG2*POICE
         CACE2=CACE2+MSI2*POICE
         CF1DT=CF1DT+F1DT*POICE
         CF2DT=CF2DT+F2DT*POICE
         CEVAP=CEVAP+EVAP*POICE
         CICOV=CICOV+POICE
         SNOWS=SNOWS+SNOW*POICE
         TG1S=TG1S+TG1*POICE
         ACE1S=ACE1S+ACE1I*POICE
         ACE2S=ACE2S+MSI2*POICE
         TG2S=TG2S+TG2*POICE
         EVAPS=EVAPS+EVAP*POICE
C****
  400 IF (PLICE.LE.0.) GO TO 600
C****
C**** LAND ICE
C****
      SNOW=GDATA(I,J,12)
      TG1=GDATA(I,J,13)
      TG2=GDATA(I,J,14)
      F0DT=E0(I,J,3)
         AIJ(I,J,67)=AIJ(I,J,67)+F0DT
      F1DT=E1(I,J,3)
      EVAP=EVAPOR(I,J,3)
         AIJ(I,J,63)=AIJ(I,J,63)+EVAP
C**** CALCULATE TG1
      SNANDI=SNOW+ACE1I-EVAP
      HC1=SNANDI*SHI
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT
      IF (ENRG1.LE.-TG1*HC1) GO TO 420
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
      RUN0=(ENRG1+TG1*HC1)/LHM
      TG1=0.
      SNANDI=SNANDI-RUN0
         BRUN0=BRUN0+RUN0*PLICE
         RUN0S=RUN0S+RUN0*PLICE
         AIJ(I,J,33)=AIJ(I,J,33)+RUN0
      GO TO 440
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT
  420 TG1=TG1+ENRG1/HC1
  440 IF (SNANDI.GE.ACE1I) GO TO 460
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2
      SNOW=0.
      DIFS=SNANDI-ACE1I
      TG1=(TG1*SNANDI-TG2*DIFS)/ACE1I
      EDIFS=DIFS*(TG2*SHI-LHM)
         BEDIFS=BEDIFS+EDIFS*PLICE
         AIJ(I,J,69)=AIJ(I,J,69)+EDIFS
         BDIFS=BDIFS+DIFS*PLICE
         DIFSS=DIFSS+DIFS*PLICE
         BERUN2=BERUN2+EDIFS*PLICE
         AIJ(I,J,72)=AIJ(I,J,72)+EDIFS
         BRUN2=BRUN2+DIFS*PLICE
      GO TO 500
  460 SNOW=SNANDI-ACE1I
C**** CALCULATE TG2
  500 TG2=TG2+F1DT/HC2LI
C**** RESAVE PROGNOSTIC QUANTITIES
      GDATA(I,J,12)=SNOW
      GDATA(I,J,13)=TG1
      GDATA(I,J,14)=TG2
         BSNOW=BSNOW+SNOW*PLICE
         BTG1=BTG1+TG1*PLICE
         BTG2=BTG2+TG2*PLICE
         BF1DT=BF1DT+F1DT*PLICE
         AIJ(I,J,69)=AIJ(I,J,69)+F1DT
         BEVAP=BEVAP+EVAP*PLICE
         SNOWS=SNOWS+SNOW*PLICE
         TG1S=TG1S+TG1*PLICE
         ACE1S=ACE1S+ACE1I*PLICE
         ACE2S=ACE2S+Z2LI*RHOI*PLICE
         TG2S=TG2S+TG2*PLICE
         EVAPS=EVAPS+EVAP*PLICE
C****
  600 IF (PEARTH.LE.0.) GO TO 940
C****
C**** EARTH (mostly done in subroutine EARTH, called from SURFCE)
C****
         SNOW=GDATA(I,J,2)
         TG1=GDATA(I,J,4)
         WTR1=GDATA(I,J,5)
         ACE1=GDATA(I,J,6)
         TG2=GDEEP(I,J,1)
         WTR2=GDEEP(I,J,2)
         ACE2=GDEEP(I,J,3)
         F0DT=E0(I,J,4)
         AIJ(I,J,68)=AIJ(I,J,68)+F0DT
         F1DT=E1(I,J,4)
         EVAP=EVAPOR(I,J,4)
         EVAPS=EVAPS+EVAP*PEARTH
         AIJ(I,J,64)=AIJ(I,J,64)+EVAP
         BSNOW=BSNOW+SNOW*PEARTH
         BTG1=BTG1+TG1*PEARTH
         BTG2=BTG2+TG2*PEARTH
         BWTR1=BWTR1+WTR1*PEARTH
         BACE1=BACE1+ACE1*PEARTH
         BWTR2=BWTR2+WTR2*PEARTH
         BACE2=BACE2+ACE2*PEARTH
         BF1DT=BF1DT+F1DT*PEARTH
         BEVAP=BEVAP+EVAP*PEARTH
         SNOWS=SNOWS+SNOW*PEARTH
         TG1S=TG1S+TG1*PEARTH
         WTR1S=WTR1S+WTR1*PEARTH
         ACE1S=ACE1S+ACE1*PEARTH
         WTR2S=WTR2S+WTR2*PEARTH
         ACE2S=ACE2S+ACE2*PEARTH
         TG2S=TG2S+TG2*PEARTH
C        AIJ(I,J,7)=AIJ(I,J,7)+(WTR1+ACE1)/WFC1
         AIJ(I,J,50)=AIJ(I,J,50)+(WTR1+ACE1+WTR2+ACE2)
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
  940    IF (JR.EQ.24) GO TO 950
         DJ(JR,17)=DJ(JR,17)+TG2S*DXYPJ
         DJ(JR,18)=DJ(JR,18)+TG1S*DXYPJ
         DJ(JR,30)=DJ(JR,30)+POICE*DXYPJ
         DJ(JR,45)=DJ(JR,45)+DIFSS*DXYPJ ! ocn/land ice contribution
         DJ(JR,49)=DJ(JR,49)+WTR1S*DXYPJ
         DJ(JR,50)=DJ(JR,50)+ACE1S*DXYPJ
         DJ(JR,51)=DJ(JR,51)+WTR2S*DXYPJ
         DJ(JR,52)=DJ(JR,52)+ACE2S*DXYPJ
         DJ(JR,53)=DJ(JR,53)+SNOWS*DXYPJ
         DJ(JR,54)=DJ(JR,54)+RUN0S*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
  950    AIJ(I,J,6)=AIJ(I,J,6)+EVAPS
         AIJ(I,J,28)=AIJ(I,J,28)+TG1S
  960 CONTINUE
C**** LONGITUDINALLY INTEGRATED QUANTITIES FOR DIAGJ
         CJ(J,15)=CJ(J,15)+CF2DT
         AJ(J,17)=AJ(J,17)+ATG2
         BJ(J,17)=BJ(J,17)+BTG2
         CJ(J,17)=CJ(J,17)+CTG2
         AJ(J,18)=AJ(J,18)+ATG1
         BJ(J,18)=BJ(J,18)+BTG1
         CJ(J,18)=CJ(J,18)+CTG1
         AJ(J,19)=AJ(J,19)+AEVAP
         BJ(J,19)=BJ(J,19)+BEVAP
         CJ(J,19)=CJ(J,19)+CEVAP
         CJ(J,30)=CJ(J,30)+CICOV
         AJ(J,33)=AJ(J,33)+AOTDT
         CJ(J,33)=CJ(J,33)+COTDT
         AJ(J,34)=AJ(J,34)+ATG3
C        BJ(J,40)=BJ(J,40)+BERUN0 ! land ice (Tg=0)
         BJ(J,41)=BJ(J,41)+BEDIFS
         CJ(J,41)=CJ(J,41)+CEDIFS
         BJ(J,42)=BJ(J,42)+BF1DT
         CJ(J,42)=CJ(J,42)+CF1DT
         AJ(J,43)=AJ(J,43)+AEFO
         BJ(J,43)=BJ(J,43)+BERUN2
         CJ(J,43)=CJ(J,43)+(CERUN2+CEFI)
         BJ(J,45)=BJ(J,45)+BDIFS
         CJ(J,45)=CJ(J,45)+CDIFS
         AJ(J,46)=AJ(J,46)+AIFO
         BJ(J,46)=BJ(J,46)+BRUN2
         CJ(J,46)=CJ(J,46)+(CRUN2+CIFI)
         AJ(J,47)=AJ(J,47)+ARUN4
         CJ(J,47)=CJ(J,47)+CRUN4
         AJ(J,48)=AJ(J,48)+AERUN4
         CJ(J,48)=CJ(J,48)+CERUN4
         BJ(J,49)=BJ(J,49)+BWTR1
         BJ(J,50)=BJ(J,50)+BACE1
         BJ(J,51)=BJ(J,51)+BWTR2
         BJ(J,52)=BJ(J,52)+BACE2
         CJ(J,52)=CJ(J,52)+CACE2
         BJ(J,53)=BJ(J,53)+BSNOW
         CJ(J,53)=CJ(J,53)+CSNOW
         BJ(J,54)=BJ(J,54)+BRUN0
         CJ(J,54)=CJ(J,54)+CRUN0
  980 CONTINUE
      RETURN
      END
      SUBROUTINE DRYCNV
C****
C**** THIS SUBROUTINE MIXES AIR CAUSED BY DRY CONVECTION.  SINCE DRY
C**** CONVECTION IN THE BOUNDARY LAYER IS DONE IN SUBROUTINE SURFCE,
C**** THIS ROUTINE ONLY CHECKS LAYERS 2 TO LM.
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,sha
      USE E001M12_COM
      USE GEOM
      USE SOMTQ_COM
      USE DAGCOM, only : ajl
      USE DYNAMICS, only : pk
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON/WORK2/UT(IM,JM,LM),VT(IM,JM,LM),DP(LM)
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      LOGICAL POLE

      RVX=0.
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      DTSRCE=DT*NDYN
      DO 50 L=1,LM
      DO 50 J=2,JM
      DO 50 I=1,IM
      UT(I,J,L)=U(I,J,L)
   50 VT(I,J,L)=V(I,J,L)
C**** OUTSIDE LOOPS OVER J AND I
      DO 500 J=1,JM
      POLE=.FALSE.
      IF (J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.

      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)
C****
C**** MAIN LOOP
C****
      IM1=IM
      DO 500 I=1,IMAX
         DO K=1,KMAX
            RA(K)=RAJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
      LMAX=1
  130 LMIN=LMAX+1
      IF (LMIN.GE.LM) GO TO 500
      LMAX=LMIN
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)) GO TO 130
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS
C**** MIX THROUGH TWO LOWER LAYERS
      PIJBOT=P(I,J)
      IF(LMIN.GE.LS1) PIJBOT=PSF-PTOP
      DP(LMIN)=PIJBOT*DSIG(LMIN)
      PIJ=PIJBOT
      IF(LMIN+1.EQ.LS1) PIJ=PSF-PTOP
      DP(LMIN+1)=PIJ*DSIG(LMIN+1)
      PKMS=PK(LMIN,I,J)*DP(LMIN)+PK(LMIN+1,I,J)*DP(LMIN+1)
      THPKMS=T(I,J,LMIN)*(PK(LMIN,I,J)*DP(LMIN))
     *  +T(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMS=Q(I,J,LMIN)*DP(LMIN)+Q(I,J,LMIN+1)*DP(LMIN+1)
C**** sum moments to mix over unstable layers
      TXS=   TX(I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     *       TX(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      TYS=   TY(I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     *       TY(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      TXXS= TXX(I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  )) +
     *      TXX(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      TYYS= TYY(I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  )) +
     *      TYY(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      TXYS= TXY(I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  )) +
     *      TXY(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QXS  =  QX(I,J,LMIN)*DP(LMIN) +  QX(I,J,LMIN+1)*DP(LMIN+1)
      QYS  =  QY(I,J,LMIN)*DP(LMIN) +  QY(I,J,LMIN+1)*DP(LMIN+1)
      QXXS = QXX(I,J,LMIN)*DP(LMIN) + QXX(I,J,LMIN+1)*DP(LMIN+1)
      QYYS = QYY(I,J,LMIN)*DP(LMIN) + QYY(I,J,LMIN+1)*DP(LMIN+1)
      QXYS = QXY(I,J,LMIN)*DP(LMIN) + QXY(I,J,LMIN+1)*DP(LMIN+1)
      IF (LMIN+1.GE.LM) GO TO 150
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO 140 L=LMIN+2,LM
      IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 160
      IF(L.EQ.LS1) PIJ=PSF-PTOP
      DP(L)=PIJ*DSIG(L)
      PKMS=PKMS+(PK(L,I,J)*DP(L))
      THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*DP(L))
      QMS=QMS+Q(I,J,L)*DP(L)
      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*DP(L))
      TXS=   TXS +  TX(I,J,L)*(PK(L,I,J)*DP(L))
      TYS=   TYS +  TY(I,J,L)*(PK(L,I,J)*DP(L))
      TXXS= TXXS + TXX(I,J,L)*(PK(L,I,J)*DP(L))
      TYYS= TYYS + TYY(I,J,L)*(PK(L,I,J)*DP(L))
      TXYS= TXYS + TXY(I,J,L)*(PK(L,I,J)*DP(L))
      QXS  =  QXS +  QX(I,J,L)*DP(L)
      QYS  =  QYS +  QY(I,J,L)*DP(L)
      QXXS = QXXS + QXX(I,J,L)*DP(L)
      QYYS = QYYS + QYY(I,J,L)*DP(L)
      QXYS = QXYS + QXY(I,J,L)*DP(L)
  140 THETA=TVMS/PKMS
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PIJBOT*SIGE(LMIN)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      PIJ=P(I,J)
      DO 180 L=LMIN,LMAX
      IF(L.GE.LS1) PIJ=PSF-PTOP
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(L,I,J)*PIJ
         AJL(J,L,55)=AJL(J,L,55)+(QMS-Q(I,J,L))*PIJ*LHE/SHA
      T(I,J,L)=THM
       TX(I,J,L) = TXS/PKMS
       TY(I,J,L) = TYS/PKMS
       TZ(I,J,L) = 0.
      TXX(I,J,L) = TXXS/PKMS
      TYY(I,J,L) = TYYS/PKMS
      TXY(I,J,L) = TXYS/PKMS
      TZZ(I,J,L) = 0.
      TZX(I,J,L) = 0.
      TYZ(I,J,L) = 0.
      Q(I,J,L)=QMS
       QX(I,J,L) = QXS*RDP
       QY(I,J,L) = QYS*RDP
       QZ(I,J,L) = 0.
      QXX(I,J,L) = QXXS*RDP
      QYY(I,J,L) = QYYS*RDP
      QXY(I,J,L) = QXYS*RDP
      QZZ(I,J,L) = 0.
      QZX(I,J,L) = 0.
      QYZ(I,J,L) = 0.
  180 CONTINUE
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      PIJ=P(I,J)
      DO L=LMIN,LMAX
         IF(L.EQ.LS1) PIJ=PSF-PTOP
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*DP(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*DP(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      PIJ=P(I,J)
      DO L=LMIN,LMAX
         IF(L.EQ.LS1) PIJ=PSF-PTOP
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c the following line gives bytewise different ajl
            AJL(IDJ(K),L,38)=AJL(IDJ(K),L,38)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PIJ*RA(K)
         ENDDO
      ENDDO
      GO TO 130
  500 IM1=I
      RETURN
      END
      SUBROUTINE SDRAG
C****
C**** THIS SUBROUTINE PUTS A DRAG ON THE WINDS ON THE TOP LAYER OF
C**** THE ATMOSPHERE
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
      USE E001M12_COM
      USE GEOM
      USE DAGCOM, only : aij
      USE DYNAMICS, only : pk
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION XCDLM(2)
      NAMELIST/SDRNML/XCDLM
      INTEGER :: IFIRST = 1
      IF(IFIRST.EQ.1) THEN
        READ(5,SDRNML)
        WRITE(6,SDRNML)
        IFIRST=0
      END IF
      DO 100 J=2,JM
      I=IM
      DO 100 IP1=1,IM
      PIJU=PSF-PTOP
      WLM=SQRT(U(I,J,LM)*U(I,J,LM)+V(I,J,LM)*V(I,J,LM))
      RHO=(PIJU*SIGE(LM+1)+PTOP)/(RGAS*T(I,J,LM)*PK(LM,I,J))
      CDN=XCDLM(1)+XCDLM(2)*WLM
         AIJ(I,J,59)=AIJ(I,J,59)+WLM
      X=NDYN*DT*RHO*CDN*WLM*GRAV/(PIJU*DSIG(LM))
      IF(X.GT.1) THEN
        write(99,*) 'SDRAG: TAU,I,J,PIJU,X,RHO,CDN,U(I,J,LM),V(I,J,LM)',
     *   TAU,I,J,PIJU,X,RHO,CDN,U(I,J,LM),V(I,J,LM)
        X=1.
      END IF
      U(I,J,LM)=U(I,J,LM)*(1.-X)
      V(I,J,LM)=V(I,J,LM)*(1.-X)
  100 I=IP1
      RETURN
      END
