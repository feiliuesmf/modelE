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
      USE CONSTANT, only : lhm,rhow,shw,shi
      USE E001M12_COM
      USE GEOM
      USE CLD01_COM_E001, only : prec,tprec
      USE DAGCOM, only : aj,bj,cj,areg,aij,jreg,ij_snwf,ij_f0oc,ij_f0oi
     *     ,ij_f0li,ij_f1li,ij_erun2,ij_runli,ij_f0e,j_eprcp,j_difs
     *     ,j_run1,ij_prec,ij_neth,j_edifs,j_erun2,j_imelt,j_run2
     *     ,j_dwtr2
      USE OCEAN, only : tocean,oa,z1o,prec_oc
      USE SEAICE, only : ace1i,prec_si
      USE SEAICE_COM, only : rsi,msi
      USE LAKES_COM, only : mwl,gml
      USE LANDICE, only : precli
      IMPLICIT NONE
      REAL*8 MSI1, MSI2
!@var QFIXR true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL QFIXR

      INTEGER I,J,L,K,IMAX,JR
      REAL*8 PLAND,PWATER,PLICE,PEARTH
     *     ,ROICE,POICE,POCEAN,DXYPJ,PRCP,TPRCP,EPRCP
     *     ,ENRGP,TGW,WTRO,RUN4,ERUN4,SNOW,WTRW0,ENRGW0,DIFS,EDIFS,WTRW
     *     ,ENRGW,TG1,TG2,TG3,TG4,RUN0,ERUN2,POLAKE,PLKICE
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** TOCEAN(1)  OCEAN TEMPERATURE (C)
C****  RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****  MSI  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
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
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=1,JM
      IMAX=IMAXJ(J)
      DO I=1,IMAX
      IF (PREC(I,J).LE.0.) CYCLE
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      POCEAN=PWATER-POICE
      POLAKE=FLAKE(I,J)*(1.-ROICE)
      PLKICE=FLAKE(I,J)*ROICE
         JR=JREG(I,J)
         DXYPJ=DXYP(J)
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
      PRCP=PREC(I,J)
      TPRCP=TPREC(I,J)
      IF (TPRCP.GE.0.) THEN
C       EPRCP=PRCP*TPRCP*SHW
        EPRCP=0.
        ENRGP=EPRCP
      ELSE
C       EPRCP=PRCP*TPRCP*SHI
        EPRCP=0.
        ENRGP=EPRCP-PRCP*LHM
          AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCP
      END IF
C**** Precipitation diagnostics (put in cloud driver?)
         AREG(JR,J_EPRCP)=AREG(JR,J_EPRCP)+ENRGP*DXYPJ
         AIJ(I,J,IJ_PREC)=AIJ(I,J,IJ_PREC)+PRCP
         AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+ENRGP
      IF (PWATER.LE.0.) GO TO 400
C****
C**** OCEAN AND LAKE
C****
        OA(I,J,4)=OA(I,J,4)+ENRGP
        AJ(J,J_EPRCP)=AJ(J,J_EPRCP)+ENRGP*POCEAN
        AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+ENRGP*POCEAN
      IF (POLAKE.GT.0.) THEN
        MWL(I,J) = MWL(I,J) +  PRCP*POLAKE*DXYP(J)
        GML(I,J) = GML(I,J) + ENRGP*POLAKE*DXYP(J)
      END IF
      IF (KOCEAN .EQ. 1) THEN  ! .and. .not. lake?
        TGW=TOCEAN(1,I,J)
        WTRO=Z1O(I,J)*RHOW
        RUN4=PRCP
        ERUN4=RUN4*TGW*SHW
        AJ(J,J_DWTR2)=AJ(J,J_DWTR2)+ERUN4*POCEAN
        AJ(J,J_RUN2)=AJ(J,J_RUN2)+RUN4*POCEAN
        IF (ROICE.GT.0.) GO TO 110
        CALL PREC_OC(TGW,WTRO,PRCP,ENRGP,ERUN4)
        TOCEAN(1,I,J)=TGW
        GO TO 400
      END IF
C****
  100 IF (POICE.LE.0.) GO TO 400
C****
C**** OCEAN AND LAKE ICE
C****
  110 CONTINUE
      SNOW=GDATA(I,J,1)
      MSI2=MSI(I,J)
      MSI1 = SNOW + ACE1I
         CJ(J,J_EPRCP)=CJ(J,J_EPRCP)+ENRGP*POICE
         AIJ(I,J,IJ_F0OI)=AIJ(I,J,IJ_F0OI)+ENRGP*POICE
      IF (KOCEAN .EQ. 1) THEN
        CJ(J,J_DWTR2)=CJ(J,J_DWTR2)+ERUN4*POICE
        CJ(J,J_RUN2)=CJ(J,J_RUN2)+RUN4*POICE
        WTRW0=WTRO-ROICE*(MSI1+MSI2)
        ENRGW0=WTRW0*TOCEAN(1,I,J)*SHW
      END IF
      TG1 = GDATA(I,J,3)  ! first layer sea ice temperature
      TG2 = GDATA(I,J,7)  ! second layer sea ice temperature
      TG3 = GDATA(I,J,15) ! third layer sea ice temperature
      TG4 = GDATA(I,J,16) ! fourth layer sea ice temperature
      QFIXR = .TRUE.
      IF (KOCEAN.eq.1) QFIXR=.FALSE.
c      IF (FLAKE(I,J).gt.0) QFIXR=.FALSE.   ! will soon be implemented

C***  CALL SUBROUTINE FOR CALCULATION OF PRECIPITATION OVER SEA ICE
      CALL PREC_SI(SNOW,MSI2,MSI1,TG1,TG2,TG3,TG4,PRCP,TPRCP,EPRCP
     *     ,RUN0,DIFS,EDIFS,ERUN2,QFIXR)

      IF (.not. QFIXR) THEN
        MSI(I,J)=MSI2
        GDATA(I,J,15)=TG3
        GDATA(I,J,16)=TG4
      END IF
      GDATA(I,J,1)=SNOW
      GDATA(I,J,7)=TG2
      GDATA(I,J,3)=TG1

C***  ACCUMULATE DIAGNOSTICS
      IF (KOCEAN .NE. 1) THEN
        AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS*POICE*DXYPJ
        CJ(J,J_ERUN2)=CJ(J,J_ERUN2)+ERUN2*POICE
        CJ(J,J_IMELT)=CJ(J,J_IMELT)+DIFS*POICE
      ELSE
C**** CALCULATE THE COMPOSITE WATER MASS AND WATER ENERGY
        WTRW=WTRW0+ROICE*(RUN0-RUN4)
        ENRGW=ENRGW0-ROICE*ERUN4+(1.-ROICE)*ENRGP
        TGW=ENRGW/(WTRW*SHW)    ! mixed layer temperature
        TOCEAN(1,I,J)=TGW
      END IF
C**** If lake ice, add runoff to lake amount (E. of runoff always = 0)
      IF (PLKICE.gt.0) THEN
        MWL(I,J) = MWL(I,J) + RUN0*POICE*DXYP(J)
      END IF
         CJ(J,J_EDIFS)=CJ(J,J_EDIFS)+EDIFS*POICE
         CJ(J,J_DIFS)=CJ(J,J_DIFS)+DIFS*POICE
         CJ(J,J_RUN1)=CJ(J,J_RUN1)+RUN0*POICE
         AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0*POICE*DXYPJ
C****
  400 IF (PLICE.LE.0.) GO TO 600
C****
C**** LAND ICE
C****
      SNOW=GDATA(I,J,12)
      TG1=GDATA(I,J,13)
      TG2=GDATA(I,J,14)
         BJ(J,J_EPRCP)=BJ(J,J_EPRCP)+ENRGP*PLICE
         AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP

      CALL PRECLI(SNOW,TG1,TG2,PRCP,EPRCP,TPRCP,EDIFS,DIFS,ERUN2,RUN0)

      GDATA(I,J,12)=SNOW
      GDATA(I,J,13)=TG1
      GDATA(I,J,14)=TG2
         BJ(J,J_EDIFS)=BJ(J,J_EDIFS)+EDIFS*PLICE
         BJ(J,J_DIFS)=BJ(J,J_DIFS)+DIFS*PLICE
         BJ(J,J_ERUN2)=BJ(J,J_ERUN2)+ERUN2*PLICE
         BJ(J,J_IMELT)=BJ(J,J_IMELT)+DIFS*PLICE
         BJ(J,J_RUN1)=BJ(J,J_RUN1)+RUN0*PLICE
C        BJ(J,J_ERUN1)=BJ(J,J_ERUN1)+ERUN0*PLICE ! land ice (Tg=0)
         AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS*PLICE*DXYPJ
         AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0*PLICE*DXYPJ
         AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+ERUN2
         AIJ(I,J,IJ_F1LI)=AIJ(I,J,IJ_F1LI)+EDIFS
         AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
C**** Add runoff to lake amount (E. of runoff always = 0)
      MWL(I,J) = MWL(I,J) +  RUN0*PLICE*DXYP(J)
C****
  600 IF (PEARTH.GT.0.) THEN
C****
C**** EARTH  (all else is done in subroutine EARTH, called from SURFCE)
C****
        BJ(J,J_EPRCP)=BJ(J,J_EPRCP)+ENRGP*PEARTH
        AIJ(I,J,IJ_F0E)=AIJ(I,J,IJ_F0E)+ENRGP
      END IF
C****
      END DO
      END DO
      RETURN
      END
      SUBROUTINE COSZ0
C****
C**** calculates the Earth's zenith angle, weighted either
C**** by time or by sun light.
C****
      USE CONSTANT, only : twopi
      USE E001M12_COM
      USE GEOM, only : DLAT,DLON
      USE RADNCB, only : COSD,SIND
      IMPLICIT NONE

      REAL*8, DIMENSION(IM) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
      REAL*8, DIMENSION(IM,JM) :: COSZ,COSZA
      REAL*8, DIMENSION(IM) :: RI,SINI,COSI
      REAL*8, DIMENSION(JM) :: SINJ,COSJ
      COMMON/WORK5/LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8 ZERO1
      DATA ZERO1/1.D-2/
      INTEGER I,J,L
      REAL*8 PHIS,PSHIS,CPHIS,SPHIS,PHIN,SPHIN,CPHIN,PHIM,S2DAWN,S2DUSK
     *     ,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK,DAWN,SDUSK,SDAWN
     *     ,CJCD,SJSD,SR1,CR1,SR2,CR2,ROT1,ROT2,DROT

C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      PHIS=-.25*TWOPI
      SPHIS=-1.
      CPHIS=0.
      DO 20 J=1,JM-1
      PHIN=DLAT*(J-.5*JM)
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
      RI(I)=DLON*(I-.5)
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
!@sum  RADIA adds the radiation heating to the temperatures
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,sday,lhe,lhs,twopi,tf,stbo
      USE E001M12_COM
      USE GEOM, Jdummy=>JLAT
      USE RADNCB, only : RQT,SRHR,TRHR,FSF,COSZ1,S0X,CO2,RSDIST
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
      USE CLD01_COM_E001, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss
      USE PBLCOM, only : wsavg,tsavg
      USE DAGCOM, only : aj,bj,cj,areg,jreg,aij,ail,ajl,asjl,adaily,
     *     iwrite,jwrite,itwrite,ndlypt,j_pcldss,j_pcldmc,ij_pmccld,
     *     j_cdldep,j_pcld,ij_cldcv,ij_pcldl,ij_pcldm,ij_pcldh,
     *     ij_cldtppr,lm_req,j_srincp0,j_srnfp0,j_srnfp1,j_srincg,
     *     j_srnfg,j_brtemp,j_trincg,j_hsurf,j_hatm,j_plavis,ij_trnfp0,
     *     ij_srnfp0,ij_srincp0,ij_srnfg,ij_srincg,ij_btmpw,ij_srref
      USE DYNAMICS, only : pk,pedn,plij,pmid,pdsig
      USE OCEAN, only : tocean
      USE SEAICE_COM, only : rsi
      USE FILEMANAGER

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: COSZ2,COSZA,TRINCG,BTMPW
      REAL*8, DIMENSION(4,IM,JM) :: SNFS,TNFS
      REAL*8, DIMENSION(LM_REQ,IM,JM) :: TRHRS,SRHRS
      REAL*8, DIMENSION(IM,JM,9) :: ALB

      REAL*8, DIMENSION(LM) :: TOTCLD

      REAL*8, DIMENSION(LM+3) :: COE

      INTEGER :: IFIRST = 1, JDLAST = -9
      INTEGER I,J,L,K,KR,LR,JYFIX,JDFIX,MADVEL,J50N,J70N,LLOW,LMID
     *     ,LHI,IHOUR,IMAX,IM1,JR,IH,INCH,LMID1,LHI1,JK
      REAL*8 COEX,CO2REF,ROT1,ROT2,PLAND,PIJ,RANDSS,RANDMC,CSS
     *     ,CMC,DEPTH,QSS,TAUSSL,TAUMCL,ELHX,CLDCV,DXYPJ,ASRHR,ATRHR
     *     ,ASNFS1,BSNFS1,CSNFS1,ATNFS1,BTNFS1,CTNFS1,SRNFLG,X
!@var NRFUN indices of unit numbers for radiation routines
      INTEGER NRFUN(14),IU
!@var RUNSTR names of files for radiation routines
      CHARACTER*5 :: RUNSTR(14) = (/"RADN1","RADN2","RADN3",
     *     "RADN4","RADN5","RADN6","RADN7","RADN8",
     *     "RADN9","RADNA","RADNB","RADNC","RADND",
     *     "RADNE"/)
!@var QBIN true if files for radiation input files are binary
      LOGICAL :: QBIN(14)=(/.TRUE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.
     *     ,.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** TOCEAN(1)  OCEAN TEMPERATURE (C)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
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
      REAL*8 QSAT
      IF (MODRD.EQ.0) IDACC(2)=IDACC(2)+1
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
      CALL COSZ0
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
      LMR=LM+3
      LMRP=LMR+1
      COEX=.01*GRAV*KAPA/RGAS
      DO L=1,LM
         COE(L)=DTsrc*COEX/DSIG(L)
         PLE(L)=SIGE(L)*PSFMPT+PTOP
      END DO
      PLE(LM+1)=SIGE(LM+1)*PSFMPT+PTOP
      PLE(LM+2)=.5*PLE(LM+1)
      PLE(LMR)=.2*PLE(LM+1)
      PLE(LMR+1)=1.D-5
      PTOPTR=PTOP ! top of sigma-coord.system
      DO 40 LR=LM+1,LMR
   40 COE(LR)=DTsrc*NRAD*COEX/(PLE(LR)-PLE(LR+1))
      PTLISO=15.
      IF(CO2.LT.0.) KTREND=-NINT(CO2)
C**** Default: time-dependent So/GHG/O3/Trop-Aeros/Dust/Volc-Aeros
C****     For control runs e.g. with Jul 1,1951 atmosphere use
          JYFIX=1951
          JDFIX=182   !  Julian day (if JDFIX=0, annual cycle is used)
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
C**** set up unit numbers for 14 radiation input files
      DO IU=1,14
         IF (IU.EQ.12.OR.IU.EQ.13) CYCLE    ! not used in GCM
         call getunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO
      CALL RCOMP1 (MADVEL,NRFUN) ! MAD 1-6: O3 TrAer Dust VAer Clds SoUV
      CO2REF=FULGAS(2)
      IF(CO2.GE.0.) FULGAS(2)=CO2REF*CO2
         CALL WRITER (6,0)
c         JEQ=1+JM/2
         J50N=(50.+90.)*(JM-1)/180.+1.5
         J70N=(70.+90.)*(JM-1)/180.+1.5
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
C**** Calculate mean cosine of zenith angle for the current physics step
   50 CONTINUE
      ROT1=(TWOPI*MOD(ITIME,NDAY))/NDAY  ! MOD(ITIME,NDAY)*TWOPI/NDAY ??
      ROT2=ROT1+TWOPI*DTsrc/SDAY
      CALL COSZT (ROT1,ROT2,COSZ1)
         IHOUR=1+JHOUR
      IF (MODRD.NE.0) GO TO 900
C****
C**** Interface with radiation routines, done only every NRAD time steps
C****
C**** Calculate mean cosine of zenith angle for the full radiation step
      ROT2=ROT1+TWOPI*NRAD*DTsrc/SDAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
C****
C**** COMPUTE EARTH ALBEDOS AND OTHER PARAMETERS FOR BEGINNING OF DAY
C****
      JDAYR=JDAY
      JYEARR=JYEAR
      IF(JDAY.NE.JDLAST) CALL RCOMPT
      S0=S0X*S00WM2*RATLS0/RSDIST
      JDLAST=JDAY
C****
C**** MAIN J LOOP
C****
      DO 600 J=1,JM
      IMAX=IMAXJ(J)
      JLAT=NINT(1.+(J-1.)*45./(JM-1.))  !  j w.r.to 72x46 grid
C****
C**** MAIN I LOOP
C****
      IM1=IM
      DO 500 I=1,IMAX
         JR=JREG(I,J)
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=RSI(I,J)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      RANDSS=RANDU(X)
      RANDMC=RANDU(X)
         CSS=0.
         CMC=0.
         DEPTH=0.
      DO 240 L=1,LM
      PIJ=PLIJ(L,I,J)
      QSS=Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
      QL(L)=QSS
      IF(CLDSAV(L,I,J).LT.1.)
     *  QL(L)=(Q(I,J,L)-QSS*CLDSAV(L,I,J))/(1.-CLDSAV(L,I,J))
      TL(L)=T(I,J,L)*PK(L,I,J)
      IF(CLDSS(L,I,J).EQ.0.) RANDSS=RANDU(X)
      TAUSSL=0.
      TAUMCL=0.
      TAUWC(L)=0.
      TAUIC(L)=0.
      SIZEWC(L)=CSIZMC(L,I,J)
      SIZEIC(L)=CSIZMC(L,I,J)
         TOTCLD(L)=0.
      IF (CLDSS(L,I,J).LT.RANDSS.OR.TAUSS(L,I,J).LE.0.) GO TO 220
      TAUSSL=TAUSS(L,I,J)
      QL(L)=QSS
         CSS=1.
         AJL(J,L,28)=AJL(J,L,28)+CSS
         TOTCLD(L)=1.
  220 IF (CLDMC(L,I,J).LT.RANDMC.OR.TAUMC(L,I,J).LE.0.) GO TO 230
         CMC=1.
         AJL(J,L,29)=AJL(J,L,29)+CMC
         TOTCLD(L)=1.
         DEPTH=DEPTH+PDSIG(L,I,J)
      IF(TAUMC(L,I,J).LE.TAUSSL) GO TO 230
      TAUMCL=TAUMC(L,I,J)
      ELHX=LHE
      IF(TL(L).LE.TF) ELHX=LHS
      QL(L)=QSAT(TL(L),ELHX,PMID(L,I,J))
  230    AJL(J,L,19)=AJL(J,L,19)+TOTCLD(L)
      IF(TAUSSL+TAUMCL.GT.0.) THEN
        IF(TAUMCL.GT.TAUSSL) THEN
          IF(SVLAT(L,I,J).EQ.LHE) THEN
            TAUWC(L)=TAUMCL
          ELSE
            TAUIC(L)=TAUMCL
          END IF
        ELSE
          IF(SVLHX(L,I,J).EQ.LHE) THEN
            TAUWC(L)=TAUSSL
            SIZEWC(L)=CSIZSS(L,I,J)
          ELSE
            TAUIC(L)=TAUSSL
            SIZEIC(L)=CSIZSS(L,I,J)
          END IF
        END IF
      END IF
  240 CONTINUE
         AJ(J,J_PCLDSS)=AJ(J,J_PCLDSS)+CSS*POCEAN
         BJ(J,J_PCLDSS)=BJ(J,J_PCLDSS)+CSS*PLAND
         CJ(J,J_PCLDSS)=CJ(J,J_PCLDSS)+CSS*POICE
         AREG(JR,J_PCLDSS)=AREG(JR,J_PCLDSS)+CSS*DXYP(J)
         AJ(J,J_PCLDMC)=AJ(J,J_PCLDMC)+CMC*POCEAN
         BJ(J,J_PCLDMC)=BJ(J,J_PCLDMC)+CMC*PLAND
         CJ(J,J_PCLDMC)=CJ(J,J_PCLDMC)+CMC*POICE
         AREG(JR,J_PCLDMC)=AREG(JR,J_PCLDMC)+CMC*DXYP(J)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AJ(J,J_CDLDEP)=AJ(J,J_CDLDEP)+DEPTH*POCEAN
         BJ(J,J_CDLDEP)=BJ(J,J_CDLDEP)+DEPTH*PLAND
         CJ(J,J_CDLDEP)=CJ(J,J_CDLDEP)+DEPTH*POICE
         AREG(JR,J_CDLDEP)=AREG(JR,J_CDLDEP)+DEPTH*DXYP(J)
         CLDCV=CMC+CSS-CMC*CSS
         AJ(J,J_PCLD)=AJ(J,J_PCLD)+CLDCV*POCEAN
         BJ(J,J_PCLD)=BJ(J,J_PCLD)+CLDCV*PLAND
         CJ(J,J_PCLD)=CJ(J,J_PCLD)+CLDCV*POICE
         AREG(JR,J_PCLD)=AREG(JR,J_PCLD)+CLDCV*DXYP(J)
         AIJ(I,J,IJ_CLDCV)=AIJ(I,J,IJ_CLDCV)+CLDCV
         DO 250 L=1,LLOW
         IF (TOTCLD(L).NE.1.) GO TO 250
         AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
         GO TO 255
  250    CONTINUE
  255    DO 260 L=LMID1,LMID
         IF (TOTCLD(L).NE.1.) GO TO 260
         AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
         GO TO 265
  260    CONTINUE
  265    DO 270 L=LHI1,LHI
         IF (TOTCLD(L).NE.1.) GO TO 270
         AIJ(I,J,IJ_PCLDH)=AIJ(I,J,IJ_PCLDH)+1.
         GO TO 275
  270    CONTINUE
  275    CONTINUE
         DO 280 L=LM,1,-1
         PIJ=PLIJ(L,I,J)
         IF (TOTCLD(L).NE.1.) GO TO 280
         AIJ(I,J,IJ_CLDTPPR)=AIJ(I,J,IJ_CLDTPPR)+SIGE(L+1)*PIJ+PTOP
         GO TO 285
  280    CONTINUE
  285    DO KR=1,NDLYPT
           IF (I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
             IH=IHOUR
             DO INCH=1,NRAD
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
  300 CONTINUE
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
      DO 340 L=1,LM
      PLE(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TL(L)=T(I,J,L)*PK(L,I,J)     ! already defined
      IF(TL(L).LT.130..OR.TL(L).GT.370.) THEN
         WRITE(99,*) 'In Radia: Time,I,J,L,TL',ITime,I,J,L,TL(L)
         STOP 'In Radia: Temperature out of range'
      END IF
C**** MOISTURE VARIABLES
C---- QL(L)=Q(I,J,L)        ! already defined
  340 CONTINUE
C****
C**** RADIATION, SOLAR AND THERMAL
C****
      DO K=1,LM_REQ
        IF(RQT(K,I,J).LT.130..OR.RQT(K,I,J).GT.370.) THEN
        WRITE(99,*) 'In RADIA: Time,I,J,L,TL',ITime,I,J,LM+K,RQT(K,I,J)
        STOP 'In Radia: RQT out of range'
        END IF
        TL(LM+K)=RQT(K,I,J)
      END DO
      COSZ=COSZA(I,J)
      TGO=TOCEAN(1,I,J)+TF
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
      DO K=1,11
        PVT(K)=VDATA(I,J,K)
      END DO
      WS=WSAVG(I,J)
C-OLD FGOLDU(2)=XFRADJ*(1.-PEARTH)
C-OLD FGOLDU(3)=XFRADJ*PEARTH
      ILON=NINT(.5+(I-.5)*72./IM)
      CALL RCOMPX
      FSF(1,I,J)=FSRNFG(1)   !  ocean
      FSF(2,I,J)=FSRNFG(3)   !  ocean ice
      FSF(3,I,J)=FSRNFG(4)   !  land ice
      FSF(4,I,J)=FSRNFG(2)   !  soil
      IF(I.EQ.IWRITE.AND.J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      SRHR(1,I,J)=SRNFLB(1)
      TRHR(1,I,J)=STBO*(POCEAN*TGO**4+POICE*TGOI**4+PLICE*TGLI**4
     *  +PEARTH*TGE**4)-TRNFLB(1)
      DO 440 L=1,LM
      SRHR(L+1,I,J)=SRFHRL(L)
  440 TRHR(L+1,I,J)=-TRFCRL(L)
      DO 450 LR=1,LM_REQ
      SRHRS(LR,I,J)= SRFHRL(LM+LR)
  450 TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      DO 460 K=1,4
      SNFS(K,I,J)=SRNFLB(K+LM)
  460 TNFS(K,I,J)=TRNFLB(K+LM)-TRNFLB(1)
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
  600 CONTINUE
C****
C**** END OF MAIN LOOP FOR J INDEX
C****
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
         DO 780 J=1,JM
         DXYPJ=DXYP(J)
         IMAX=IMAXJ(J)
         DO L=1,LM
           ASRHR=0.
           ATRHR=0.
           DO I=1,IMAX
             ASRHR=ASRHR+SRHR(L+1,I,J)*COSZ2(I,J)
             ATRHR=ATRHR+TRHR(L+1,I,J)
           END DO
           AJL(J,L,9)=AJL(J,L,9)+ASRHR
           AJL(J,L,10)=AJL(J,L,10)+ATRHR
         END DO
         ASNFS1=0.
         BSNFS1=0.
         CSNFS1=0.
         ATNFS1=0.
         BTNFS1=0.
         CTNFS1=0.
         DO 770 I=1,IMAX
         COSZ=COSZ2(I,J)
         PLAND=FLAND(I,J)
         POICE=RSI(I,J)*(1.-PLAND)
         POCEAN=(1.-PLAND)-POICE
         JR=JREG(I,J)
         DO LR=1,LM_REQ
           ASJL(J,LR,3)=ASJL(J,LR,3)+SRHRS(LR,I,J)*COSZ
           ASJL(J,LR,4)=ASJL(J,LR,4)+TRHRS(LR,I,J)
         END DO
         DO KR=1,NDLYPT
           IF (I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
             IH=IHOUR
             DO INCH=1,NRAD
               IF (IH.GT.24) IH=IH-24
               ADAILY(IH,2,KR)=ADAILY(IH,2,KR)+(1.-SNFS(4,I,J)/S0)
               ADAILY(IH,3,KR)=ADAILY(IH,3,KR)+(1.-ALB(I,J,1))
               ADAILY(IH,4,KR)=ADAILY(IH,4,KR)
     *              +((SNFS(4,I,J)-SNFS(1,I,J))*COSZ-TNFS(4,I,J)
     *              +TNFS(1,I,J))
               IH=IH+1
             END DO
           END IF
         END DO
  750    CONTINUE
         AJ(J,J_SRINCP0)=AJ(J,J_SRINCP0)+(S0*COSZ)*POCEAN
         BJ(J,J_SRINCP0)=BJ(J,J_SRINCP0)+(S0*COSZ)*PLAND
         CJ(J,J_SRINCP0)=CJ(J,J_SRINCP0)+(S0*COSZ)*POICE
         AREG(JR,J_SRINCP0)=AREG(JR,J_SRINCP0)+(S0*COSZ)*DXYPJ
         AJ(J,J_SRNFP0)=AJ(J,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*POCEAN
         BJ(J,J_SRNFP0)=BJ(J,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*PLAND
         CJ(J,J_SRNFP0)=CJ(J,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*POICE
         AREG(JR,J_SRNFP0)=AREG(JR,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*DXYPJ
         ASNFS1=ASNFS1+(SNFS(1,I,J)*COSZ)*POCEAN
         BSNFS1=BSNFS1+(SNFS(1,I,J)*COSZ)*PLAND
         CSNFS1=CSNFS1+(SNFS(1,I,J)*COSZ)*POICE
         AREG(JR,J_SRNFP1)=AREG(JR,J_SRNFP1)+(SNFS(1,I,J)*COSZ)*DXYPJ
         AJ(J,J_SRINCG)=AJ(J,J_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))*POCEAN
         BJ(J,J_SRINCG)=BJ(J,J_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))*PLAND
         CJ(J,J_SRINCG)=CJ(J,J_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))*POICE
         AREG(JR,J_SRINCG)=AREG(JR,J_SRINCG)+
     *     (SRHR(1,I,J)*COSZ/(ALB(I,J,1)+1.D-20))*DXYPJ
         AJ(J,J_SRNFG)=AJ(J,J_SRNFG)+(FSF(1,I,J)*COSZ)*POCEAN
         SRNFLG=FSF(3,I,J)*FLICE(I,J)+FSF(4,I,J)*(PLAND-FLICE(I,J))
         BJ(J,J_SRNFG)=BJ(J,J_SRNFG)+(SRNFLG*COSZ)
         CJ(J,J_SRNFG)=CJ(J,J_SRNFG)+(FSF(2,I,J)*COSZ)*POICE
         AREG(JR,J_SRNFG)=AREG(JR,J_SRNFG)+(SRHR(1,I,J)*COSZ)*DXYPJ
         AJ(J,J_BRTEMP)=AJ(J,J_BRTEMP)+BTMPW(I,J)*POCEAN
         BJ(J,J_BRTEMP)=BJ(J,J_BRTEMP)+BTMPW(I,J)*PLAND
         CJ(J,J_BRTEMP)=CJ(J,J_BRTEMP)+BTMPW(I,J)*POICE
         AREG(JR,J_BRTEMP)=AREG(JR,J_BRTEMP)+BTMPW(I,J)*DXYPJ
         AJ(J,J_TRINCG)=AJ(J,J_TRINCG)+TRINCG(I,J)*POCEAN
         BJ(J,J_TRINCG)=BJ(J,J_TRINCG)+TRINCG(I,J)*PLAND
         CJ(J,J_TRINCG)=CJ(J,J_TRINCG)+TRINCG(I,J)*POICE
         AREG(JR,J_TRINCG)=AREG(JR,J_TRINCG)+TRINCG(I,J)*DXYPJ
         AJ(J,J_HSURF)=AJ(J,J_HSURF)-TNFS(4,I,J)*POCEAN
         BJ(J,J_HSURF)=BJ(J,J_HSURF)-TNFS(4,I,J)*PLAND
         CJ(J,J_HSURF)=CJ(J,J_HSURF)-TNFS(4,I,J)*POICE
         AREG(JR,J_HSURF)=AREG(JR,J_HSURF)-TNFS(4,I,J)*DXYPJ
         ATNFS1=ATNFS1-TNFS(1,I,J)*POCEAN
         BTNFS1=BTNFS1-TNFS(1,I,J)*PLAND
         CTNFS1=CTNFS1-TNFS(1,I,J)*POICE
         AREG(JR,J_HATM)=AREG(JR,J_HATM)-TNFS(1,I,J)*DXYPJ
         DO 760 K=2,9
           JK=K+J_PLAVIS-2     ! accumulate 8 radiation diags.
         AJ(J,JK)=AJ(J,JK)+(S0*COSZ)*ALB(I,J,K)*POCEAN
         BJ(J,JK)=BJ(J,JK)+(S0*COSZ)*ALB(I,J,K)*PLAND
         CJ(J,JK)=CJ(J,JK)+(S0*COSZ)*ALB(I,J,K)*POICE
  760    AREG(JR,JK)=AREG(JR,JK)+(S0*COSZ)*ALB(I,J,K)*DXYPJ
         AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)-TNFS(4,I,J)
         AIJ(I,J,IJ_SRNFP0)=AIJ(I,J,IJ_SRNFP0)+(SNFS(4,I,J)*COSZ)
         AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*COSZ)
         AIJ(I,J,IJ_SRNFG)=AIJ(I,J,IJ_SRNFG)+(SRHR(1,I,J)*COSZ)
         AIJ(I,J,IJ_SRINCG)=AIJ(I,J,IJ_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_BTMPW)=AIJ(I,J,IJ_BTMPW)+BTMPW(I,J)
         AIJ(I,J,IJ_SRREF)=AIJ(I,J,IJ_SRREF)+S0*COSZ*ALB(I,J,2)
  770    CONTINUE
         AJ(J,J_SRNFP1)=AJ(J,J_SRNFP1)+ASNFS1
         BJ(J,J_SRNFP1)=BJ(J,J_SRNFP1)+BSNFS1
         CJ(J,J_SRNFP1)=CJ(J,J_SRNFP1)+CSNFS1
         AJ(J,J_HATM)=AJ(J,J_HATM)+ATNFS1
         BJ(J,J_HATM)=BJ(J,J_HATM)+BTNFS1
         CJ(J,J_HATM)=CJ(J,J_HATM)+CTNFS1
  780    CONTINUE
         DO 790 L=1,LM
         DO 790 I=1,IM
         AIL(I,L,7)=AIL(I,L,7)+((SRHR(L+1,I,JEQ-2)*COSZ2(I,JEQ-2)+
     *     TRHR(L+1,I,JEQ-2))*DXYP(JEQ-2)+(SRHR(L+1,I,JEQ-1)*
     *     COSZ2(I,JEQ-1)+TRHR(L+1,I,JEQ-1))*DXYP(JEQ-1)+
     *     (SRHR(L+1,I,JEQ)*COSZ2(I,JEQ)+TRHR(L+1,I,JEQ))*DXYP(JEQ))
         AIL(I,L,11)=AIL(I,L,11)+(SRHR(L+1,I,J50N)*COSZ2(I,J50N)+
     *     TRHR(L+1,I,J50N))*DXYP(J50N)
  790    AIL(I,L,15)=AIL(I,L,15)+(SRHR(L+1,I,J70N)*COSZ2(I,J70N)+
     *     TRHR(L+1,I,J70N))*DXYP(J70N)
C****
C**** Update radiative equilibrium temperatures
C****
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          DO LR=1,LM_REQ
            RQT(LR,I,J)=RQT(LR,I,J)+(SRHRS(LR,I,J)*COSZ2(I,J)
     *           +TRHRS(LR,I,J))*COE(LR+LM)
          END DO
        END DO
      END DO
C****
C**** Update other temperatures every physics time step
C****
  900 DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          DO L=1,LM
            T(I,J,L)=T(I,J,L)+(SRHR(L+1,I,J)*COSZ1(I,J)+TRHR(L+1,I,J))*
     *           COE(L)/(PLIJ(L,I,J)*PK(L,I,J))
          END DO
        END DO
      END DO
C**** daily diagnostics
      DO KR=1,NDLYPT
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
      USE CONSTANT, only : sday,twopi,rhow,shw
      USE E001M12_COM
      USE GEOM
      USE PBLCOM, only : tsavg
      USE GHYCOM, only : ghdata
      USE DAGCOM, only : aj,bj,cj,areg,aij,aijg,jreg,ij_evapo,
     *     ij_f0oc,ij_rsoi,ij_msi2,ij_f0oi,ij_evapi,ij_f0li,ij_evapli,
     *     ij_runli,ij_f1li,ij_erun2,ij_f0e,ij_evape,ij_gwtr,ij_evap,
     *     ij_tg1,ij_tgo,j_tg2,j_tg1,j_rsi,j_difs,j_wtr1,j_ace1,j_wtr2,
     *     j_ace2,j_snow,j_run1,j_f2dt,j_evap,j_oht,j_omlt,j_edifs,
     *     j_f1dt,j_erun2,j_imelt,j_run2,j_dwtr2

      USE OCEAN, only : tocean,z1o,ota,otb,otc,tfo,fleadoc,osourc
      USE SEAICE, only : ace1i,sea_ice,addice
      USE SEAICE_COM, only : rsi,msi
      USE LAKES_COM, only : mwl,gml,tfl,fleadlk,t50
      USE LANDICE, only : ace2li,lndice
      IMPLICIT NONE

      REAL*8 MSI2
      REAL*8, DIMENSION(IM,JM,4) :: E0,E1,EVAPOR
      COMMON/WORK3/E0,E1,EVAPOR
      REAL*8, DIMENSION(IM,JM,3) :: GDEEP
      COMMON/oldDAG/GDEEP

      INTEGER I,J,L,KR,JR,IMAX,K
      REAL*8 FACT_T50,FACT_TSAVG
     *     ,DIFFUS,ANGLE,SINANG,SN2ANG,SN3ANG,SN4ANG,COSANG,CS2ANG
     *     ,CS3ANG,CS4ANG,BF1DT,CF1DT,AOTDT,COTDT,AEFO,CEFI,BEDIFS
     *     ,CEDIFS,BERUN0,CF2DT,BERUN2,CERUN2,AERUN4,CERUN4,ATG1,BTG1
     *     ,CTG1,ATG2,BTG2,CTG2,ATG3,AEVAP,BEVAP,CEVAP,BDIFS,CDIFS,AIFO
     *     ,CIFI,BRUN0,CRUN0,BRUN2,CRUN2,ARUN4,CRUN4,BWTR1,BACE1,BWTR2
     *     ,BACE2,CACE2,BSNOW,CSNOW,CICOV,PLAND,PWATER,PLICE,PEARTH
     *     ,ROICE,POICE,POCEAN,DXYPJ,ACE1S,ACE2S,EVAPS,RUN0S,DIFSS,EVAP
     *     ,F0DT,OTDT,RUN4,ERUN4,ENRGO,EDIFSI,DIFSI,F1DT,AIFI,SNOW
     *     ,ACE1,ACE2,SNOWS,WTR1S
     *     ,WTR2S,TG1S,TG2S,TGW,WTRO,TG1,TG2,TG3,TG4,WTR1,WTR2,RUN0,RUN4
     *     ,ERUN4,DIFSI,EDIFSI,DIFS,EDIFS,ENRGFO,ACEFO,ACE2M,ACE2F
     *     ,ENRGFI,F2DT,POLAKE,PLKICE,MSI1,HSI1,HSI2
     *     ,HSI3,HSI4,EIW0,WTRW0,ENRGW0,WTRI0
!@var QFIXR  true if RSI and MSI2 are fixed (ie. for fixed SST run)
      LOGICAL QFIXR
!@var QCMPR  true if ice should be compressed due to leads etc.
      LOGICAL QCMPR
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** TOCEAN(1)  OCEAN TEMPERATURE (C)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****   MSI  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
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

      DIFFUS=DTSRC/SDAY
      ANGLE=TWOPI*JDAY/365.
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
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      POCEAN=PWATER-POICE
      POLAKE=FLAKE(I,J)*(1.-ROICE)
      PLKICE=FLAKE(I,J)*ROICE
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
C**** OCEAN AND LAKE
C****
      EVAP=EVAPOR(I,J,1)
         ATG1=ATG1+TOCEAN(1,I,J)*POCEAN
         TG1S=TG1S+TOCEAN(1,I,J)*POCEAN
         AEVAP=AEVAP+EVAP*POCEAN
         EVAPS=EVAPS+EVAP*POCEAN
         AIJ(I,J,IJ_TGO)=AIJ(I,J,IJ_TGO)+TOCEAN(1,I,J)
         AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+EVAP*POCEAN
      IF (POLAKE.GT.0.) THEN
        MWL(I,J) = MWL(I,J) - EVAP*POLAKE*DXYP(J)
        GML(I,J) = GML(I,J) + E0(I,J,1)*POLAKE*DXYP(J)
      END IF
      TGW=TOCEAN(1,I,J)
      IF (KOCEAN .NE. 1) THEN ! .and. NOT LAKE?
           ATG2=ATG2+TOCEAN(1,I,J)*POCEAN
           TG2S=TG2S+TOCEAN(1,I,J)*POCEAN
      ELSE
        WTRO=Z1O(I,J)*RHOW
        F0DT=E0(I,J,1)
           AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+F0DT*POCEAN
        OTDT=DTSRC*(OTA(I,J,4)*SN4ANG+OTB(I,J,4)*CS4ANG
     *          +OTA(I,J,3)*SN3ANG+OTB(I,J,3)*CS3ANG
     *          +OTA(I,J,2)*SN2ANG+OTB(I,J,2)*CS2ANG
     *          +OTA(I,J,1)*SINANG+OTB(I,J,1)*COSANG+OTC(I,J))
           ATG2=ATG2+TOCEAN(2,I,J)*POCEAN
           TG2S=TG2S+TOCEAN(2,I,J)*POCEAN
           ATG3=ATG3+TOCEAN(3,I,J)*POCEAN
           AOTDT=AOTDT+OTDT*POCEAN
        RUN4=-EVAP
        ERUN4=RUN4*TGW*SHW
           AERUN4=AERUN4+ERUN4*POCEAN
           ARUN4=ARUN4+RUN4*POCEAN
        ENRGO=F0DT+OTDT-ERUN4
      END IF

C***  CALL SEA ICE SUBROUTINE
      F0DT=E0(I,J,2) ! heat flux to the top ice surface (J/m^2)
      F1DT=E1(I,J,2) ! heat flux between 1st and 2nd ice layers (J/m^2)
      EVAP=EVAPOR(I,J,2) ! evaporation/dew at the ice surface (kg/m^2)
      MSI2= MSI(I,J)
      SNOW= GDATA(I,J,1)  ! snow mass (kg/m^2)
      TG1 = GDATA(I,J,3)  ! first layer sea ice temperature
      TG2 = GDATA(I,J,7)  ! second layer sea ice temperature
      TG3 = GDATA(I,J,15) ! third layer sea ice temperature
      TG4 = GDATA(I,J,16) ! fourth layer sea ice temperature

         AIJ(I,J,IJ_RSOI) =AIJ(I,J,IJ_RSOI) +POICE
         AIJ(I,J,IJ_MSI2) =AIJ(I,J,IJ_MSI2) +MSI2*POICE
         AIJ(I,J,IJ_F0OI) =AIJ(I,J,IJ_F0OI) +F0DT*POICE
         AIJ(I,J,IJ_EVAPI)=AIJ(I,J,IJ_EVAPI)+EVAP*POICE

      QFIXR = .TRUE.
      QCMPR = .FALSE.
      IF (KOCEAN.eq.1) QFIXR=.FALSE.
      IF (KOCEAN.eq.1.and.FLAKE(I,J).le.0) QCMPR=.TRUE.
c     IF (FLAKE(I,J).gt.0) QFIXR=.FALSE.   ! will soon be implemented

      IF (KOCEAN .EQ. 1) THEN
        WTRI0=WTRO-(SNOW+ACE1I+MSI2) ! mixd layr mass below ice (kg/m^2)
        EIW0=WTRI0*TGW*SHW     ! energy of mixed layer below ice (J/m^2)
        WTRW0=WTRO-ROICE*(SNOW+ACE1I+MSI2)
        ENRGW0=WTRW0*TGW*SHW
        RUN4=-EVAP
        WTRW0 = WTRW0-ROICE*RUN4 ! water mass "+-" dew/evaporation
        ERUN4=TGW*RUN4*SHW
      END IF

      CALL SEA_ICE(DTSRC,SNOW,ROICE,TG1,TG2,TG3,TG4,MSI1,MSI2
     *     ,F0DT,F1DT,EVAP,HSI1,HSI2,HSI3,HSI4,TGW,RUN0
     *     ,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,QFIXR)

      IF (KOCEAN.EQ.1) WTRW0 = (WTRW0+ROICE*RUN0)+ROICE*ACE2M
      ACEFO=0 ; ACE2F=0. ; ENRGFO=0. ; ENRGFI=0.
      IF (.not.QFIXR) CALL OSOURC (ROICE,MSI1,MSI2,TGW,WTRO,EIW0
     *       ,OTDT,ENRGO,ERUN4,ENRGFO,ACEFO,ACE2F,WTRW0,ENRGW0,ENRGFI
     *       ,F2DT,TFO)

      CALL ADDICE (SNOW,ROICE,TG1,TG2,TG3,TG4,MSI1,MSI2,HSI1,HSI2
     *     ,HSI3,HSI4,DIFSI,EDIFSI,ENRGFO,ACEFO,ACE2F,ENRGFI,TFO,FLEADOC
     *     ,QFIXR,QCMPR)
C****
C**** RESAVE PROGNOSTIC QUANTITIES
C****
      GDATA(I,J,1) =SNOW
      GDATA(I,J,3) =TG1
      GDATA(I,J,7) =TG2
      GDATA(I,J,15)=TG3
      GDATA(I,J,16)=TG4
      IF (.not. QFIXR) THEN
         TOCEAN(1,I,J)=TGW
         RSI(I,J)=ROICE
         MSI(I,J)=MSI2
      END IF

      IF (PLKICE.GT.0.) THEN ! Add mass/energy fluxes to lake variables
        MWL(I,J) = MWL(I,J) +
     *       ((RUN0+ACE2M-ACE2F)*PLKICE - ACEFO*POLAKE)*DXYP(J)
        GML(I,J) = GML(I,J) +
     *       ((F2DT-ENRGFI)*PLKICE - ENRGFO*POLAKE)*DXYP(J)
      END IF
      IF (KOCEAN .EQ. 1) THEN
         AEFO=AEFO-ENRGFO*POCEAN
         AIFO=AIFO-ACEFO*POCEAN
         CERUN4=CERUN4+ERUN4*POICE
         CRUN4=CRUN4+RUN4*POICE
         AIFI=AIFI+ACE2M*POICE
         CEFI=CEFI-ENRGFI*POICE
         CIFI=CIFI-ACE2F*POICE
         COTDT=COTDT+OTDT*POICE
         CEDIFS=CEDIFS+EDIFSI*PWATER
         CDIFS=CDIFS+DIFSI*PWATER
         DIFSS=DIFSS+DIFSI*PWATER
      ELSE
         CEDIFS=CEDIFS+EDIFS*POICE
         CDIFS=CDIFS+DIFS*POICE
         CERUN2=CERUN2+EDIFS*POICE
         CRUN2=CRUN2+DIFS*POICE
         DIFSS=DIFSS+DIFS*POICE
      END IF
         CRUN0=CRUN0+RUN0*POICE
         RUN0S=RUN0S+RUN0*POICE
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
      F1DT=E1(I,J,3)
      EVAP=EVAPOR(I,J,3)
         AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+F0DT
         AIJ(I,J,IJ_EVAPLI)=AIJ(I,J,IJ_EVAPLI)+EVAP

      CALL LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES
      GDATA(I,J,12)=SNOW
      GDATA(I,J,13)=TG1
      GDATA(I,J,14)=TG2
         BRUN0=BRUN0+RUN0*PLICE
         RUN0S=RUN0S+RUN0*PLICE
         AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
         BEDIFS=BEDIFS+EDIFS*PLICE
         AIJ(I,J,IJ_F1LI)=AIJ(I,J,IJ_F1LI)+EDIFS
         BDIFS=BDIFS+DIFS*PLICE
         DIFSS=DIFSS+DIFS*PLICE
         BERUN2=BERUN2+EDIFS*PLICE
         AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+EDIFS
         BRUN2=BRUN2+DIFS*PLICE
         BSNOW=BSNOW+SNOW*PLICE
         BTG1=BTG1+TG1*PLICE
         BTG2=BTG2+TG2*PLICE
         BF1DT=BF1DT+F1DT*PLICE
         AIJ(I,J,IJ_F1LI)=AIJ(I,J,IJ_F1LI)+F1DT
         BEVAP=BEVAP+EVAP*PLICE
         SNOWS=SNOWS+SNOW*PLICE
         TG1S=TG1S+TG1*PLICE
         ACE1S=ACE1S+ACE1I*PLICE
         ACE2S=ACE2S+ACE2LI*PLICE
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
         AIJ(I,J,IJ_F0E)=AIJ(I,J,IJ_F0E)+F0DT
         F1DT=E1(I,J,4)
         EVAP=EVAPOR(I,J,4)
         EVAPS=EVAPS+EVAP*PEARTH
         AIJ(I,J,IJ_EVAPE)=AIJ(I,J,IJ_EVAPE)+EVAP
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
C        AIJ(I,J,IJ_SSAT)=AIJ(I,J,IJ_SSAT)+(WTR1+ACE1)/WFC1
         AIJ(I,J,IJ_GWTR)=AIJ(I,J,IJ_GWTR)+(WTR1+ACE1+WTR2+ACE2)
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
  940    IF (JR.EQ.24) GO TO 950
         AREG(JR,J_TG2)=AREG(JR,J_TG2)+TG2S*DXYPJ
         AREG(JR,J_TG1)=AREG(JR,J_TG1)+TG1S*DXYPJ
         AREG(JR,J_RSI)=AREG(JR,J_RSI)+POICE*DXYPJ
         AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFSS*DXYPJ !ocn/land ice cont.
         AREG(JR,J_WTR1)=AREG(JR,J_WTR1)+WTR1S*DXYPJ
         AREG(JR,J_ACE1)=AREG(JR,J_ACE1)+ACE1S*DXYPJ
         AREG(JR,J_WTR2)=AREG(JR,J_WTR2)+WTR2S*DXYPJ
         AREG(JR,J_ACE2)=AREG(JR,J_ACE2)+ACE2S*DXYPJ
         AREG(JR,J_SNOW)=AREG(JR,J_SNOW)+SNOWS*DXYPJ
         AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0S*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
  950    AIJ(I,J,IJ_EVAP)=AIJ(I,J,IJ_EVAP)+EVAPS
         AIJ(I,J,IJ_TG1) =AIJ(I,J,IJ_TG1)+TG1S
         DO K=1,4
            AIJG(I,J,K)=AIJG(I,J,K)+GHDATA(I,J,K)
         END DO
         DO K=7,10
            AIJG(I,J,K)=AIJG(I,J,K)+GHDATA(I,J,K)
         END DO
         AIJG(I,J,28)=AIJG(I,J,28)+GHDATA(I,J,28)
         AIJG(I,J,29)=AIJG(I,J,29)+GHDATA(I,J,29)
  960 CONTINUE
C**** LONGITUDINALLY INTEGRATED QUANTITIES FOR DIAGJ
         CJ(J,J_F2DT)=CJ(J,J_F2DT)+CF2DT
         AJ(J,J_TG2)=AJ(J,J_TG2)+ATG2
         BJ(J,J_TG2)=BJ(J,J_TG2)+BTG2
         CJ(J,J_TG2)=CJ(J,J_TG2)+CTG2
         AJ(J,J_TG1)=AJ(J,J_TG1)+ATG1
         BJ(J,J_TG1)=BJ(J,J_TG1)+BTG1
         CJ(J,J_TG1)=CJ(J,J_TG1)+CTG1
         AJ(J,J_EVAP)=AJ(J,J_EVAP)+AEVAP
         BJ(J,J_EVAP)=BJ(J,J_EVAP)+BEVAP
         CJ(J,J_EVAP)=CJ(J,J_EVAP)+CEVAP
         CJ(J,J_RSI)=CJ(J,J_RSI)+CICOV
         AJ(J,J_OHT)=AJ(J,J_OHT)+AOTDT
         CJ(J,J_OHT)=CJ(J,J_OHT)+COTDT
         AJ(J,J_OMLT)=AJ(J,J_OMLT)+ATG3
C        BJ(J,J_ERUN1)=BJ(J,J_ERUN1)+BERUN0 ! land ice (Tg=0)
         BJ(J,J_EDIFS)=BJ(J,J_EDIFS)+BEDIFS
         CJ(J,J_EDIFS)=CJ(J,J_EDIFS)+CEDIFS
         BJ(J,J_F1DT)=BJ(J,J_F1DT)+BF1DT
         CJ(J,J_F1DT)=CJ(J,J_F1DT)+CF1DT
         AJ(J,J_ERUN2)=AJ(J,J_ERUN2)+AEFO
         BJ(J,J_ERUN2)=BJ(J,J_ERUN2)+BERUN2
         CJ(J,J_ERUN2)=CJ(J,J_ERUN2)+(CERUN2+CEFI)
         BJ(J,J_DIFS)=BJ(J,J_DIFS)+BDIFS
         CJ(J,J_DIFS)=CJ(J,J_DIFS)+CDIFS
         AJ(J,J_IMELT)=AJ(J,J_IMELT)+AIFO
         BJ(J,J_IMELT)=BJ(J,J_IMELT)+BRUN2
         CJ(J,J_IMELT)=CJ(J,J_IMELT)+(CRUN2+CIFI)
         AJ(J,J_RUN2)=AJ(J,J_RUN2)+ARUN4
         CJ(J,J_RUN2)=CJ(J,J_RUN2)+CRUN4
         AJ(J,J_DWTR2)=AJ(J,J_DWTR2)+AERUN4
         CJ(J,J_DWTR2)=CJ(J,J_DWTR2)+CERUN4
         BJ(J,J_WTR1)=BJ(J,J_WTR1)+BWTR1
         BJ(J,J_ACE1)=BJ(J,J_ACE1)+BACE1
         BJ(J,J_WTR2)=BJ(J,J_WTR2)+BWTR2
         BJ(J,J_ACE2)=BJ(J,J_ACE2)+BACE2
         CJ(J,J_ACE2)=CJ(J,J_ACE2)+CACE2
         BJ(J,J_SNOW)=BJ(J,J_SNOW)+BSNOW
         CJ(J,J_SNOW)=CJ(J,J_SNOW)+CSNOW
         BJ(J,J_RUN1)=BJ(J,J_RUN1)+BRUN0
         CJ(J,J_RUN1)=CJ(J,J_RUN1)+CRUN0
  980 CONTINUE
      RETURN
      END
      SUBROUTINE DRYCNV(LBASE_MIN,LBASE_MAX)
C****
C**** THIS SUBROUTINE MIXES AIR CAUSED BY DRY CONVECTION.
C**** THIS VERSION CHECKS BASE LAYERS LBASE_MIN TO LBASE_MAX.
C****
      USE CONSTANT, only : lhe,sha
      USE E001M12_COM
      USE GEOM
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
      USE DAGCOM, only : ajl
      USE DYNAMICS, only : pk,pdsig,plij
      USE PBLCOM, only : dclev
      IMPLICIT NONE

      integer, intent(in) :: LBASE_MIN,LBASE_MAX
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT
      REAL*8, DIMENSION(LM) :: DP
      COMMON/WORK2/UT,VT,DP
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      LOGICAL POLE
      INTEGER I,J,L,K,IMAX,KMAX,IM1,LMAX,LMIN

      DOUBLE PRECISION, DIMENSION(NMOM) :: TMOMS,QMOMS
      REAL*8 RVX,DOK,PIJBOT,PIJ,PKMS,THPKMS,QMS
     *     ,TVMS,THETA,RDP,THM

      if(LBASE_MAX.GE.LM) stop 'DRYCNV: LBASE_MAX.GE.LM'
      RVX=0.
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      DO 50 L=1,LM
      DO 50 J=2,JM
      DO 50 I=1,IM
      UT(I,J,L)=U(I,J,L)
   50 VT(I,J,L)=V(I,J,L)
C**** OUTSIDE LOOPS OVER J AND I
      JLOOP: DO J=1,JM
      POLE=.FALSE.
      IF (J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.

      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)
C****
C**** MAIN LOOP
C****
      IM1=IM
      ILOOP: DO I=1,IMAX
         DO K=1,KMAX
            RA(K)=RAJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
      LMAX=LBASE_MIN-1
      lbase_loop: do while(lmax.lt.lbase_max)
      LMIN=LMAX+1
      LMAX=LMIN
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)) cycle lbase_loop
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS
C**** MIX THROUGH TWO LOWER LAYERS
      PIJBOT=PLIJ(LMIN,I,J)
      DP(LMIN)=PDSIG(LMIN,I,J)
      PIJ=PLIJ(LMIN+1,I,J)
      DP(LMIN+1)=PDSIG(LMIN+1,I,J)
      PKMS=PK(LMIN,I,J)*DP(LMIN)+PK(LMIN+1,I,J)*DP(LMIN+1)
      THPKMS=T(I,J,LMIN)*(PK(LMIN,I,J)*DP(LMIN))
     *  +T(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMS=Q(I,J,LMIN)*DP(LMIN)+Q(I,J,LMIN+1)*DP(LMIN+1)
C**** sum moments to mix over unstable layers
      TMOMS(XYMOMS) =
     &     TMOM(XYMOMS,I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     &     TMOM(XYMOMS,I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMOMS(XYMOMS) =
     &     QMOM(XYMOMS,I,J,LMIN  )*(DP(LMIN  ))  +
     &     QMOM(XYMOMS,I,J,LMIN+1)*(DP(LMIN+1))
      IF (LMIN+1.GE.LM) GO TO 150
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO 140 L=LMIN+2,LM
      IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 160
      PIJ=PLIJ(L,I,J)
      DP(L)=PDSIG(L,I,J)
      PKMS=PKMS+(PK(L,I,J)*DP(L))
      THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*DP(L))
      QMS=QMS+Q(I,J,L)*DP(L)
      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*DP(L))
      TMOMS(XYMOMS) = TMOMS(XYMOMS) +
     &     TMOM(XYMOMS,I,J,L)*(PK(L,I,J)*DP(L))
      QMOMS(XYMOMS) = QMOMS(XYMOMS) +
     &     QMOM(XYMOMS,I,J,L)*DP(L)
  140 THETA=TVMS/PKMS
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PIJBOT*SIGE(LMIN)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      DO 180 L=LMIN,LMAX
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(L,I,J)*PLIJ(L,I,J)
         AJL(J,L,55)=AJL(J,L,55)+(QMS-Q(I,J,L))*PDSIG(L,I,J)*LHE/SHA
      T(I,J,L)=THM
      TMOM(XYMOMS,I,J,L)=TMOMS(XYMOMS)/PKMS
      TMOM(ZMOMS,I,J,L)=0.
      Q(I,J,L)=QMS
      QMOM(XYMOMS,I,J,L)=QMOMS(XYMOMS)*RDP
      QMOM(ZMOMS,I,J,L)=0.
  180 CONTINUE
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      DO L=LMIN,LMAX
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*DP(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*DP(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      DO L=LMIN,LMAX
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c the following line gives bytewise different ajl
            AJL(IDJ(K),L,38)=AJL(IDJ(K),L,38)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PLIJ(L,I,J)*RA(K)
         ENDDO
      ENDDO
      enddo lbase_loop
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
      if(lbase_min.eq.1) then ! was called from surfce
         DCLEV(I,J)=LMAX
      endif
      IM1=I
      ENDDO ILOOP
      ENDDO JLOOP
      RETURN
      END
      SUBROUTINE SDRAG !(DT1) in stratosphere model
C****
C**** THIS SUBROUTINE PUTS A DRAG ON THE WINDS ON THE TOP LAYER OF
C**** THE ATMOSPHERE
C****
      USE CONSTANT, only : grav,rgas
      USE E001M12_COM
      USE GEOM
      USE DAGCOM, only : aij, ij_wlm,ajl,ij_sdrag
      USE DYNAMICS, only : pk
      IMPLICIT NONE

      INTEGER I,J,IP1,L
      REAL*8 PIJU,WL,RHO,CDN,X,BYPIJU

      PIJU = PSFMPT
      BYPIJU=1./PSFMPT
      DO L=LM,LM
      DO 100 J=2,JM
      I=IM
      DO 100 IP1=1,IM
      WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
      RHO=(PIJU*SIGE(L+1)+PTOP)/(RGAS*T(I,J,L)*PK(L,I,J))
      CDN=XCDLM(1)+XCDLM(2)*WL
         AIJ(I,J,IJ_WLM)=AIJ(I,J,IJ_WLM)+WL
      X=DTsrc*RHO*CDN*WL*GRAV*BYDSIG(L)*BYPIJU
      IF(X.GT.1) THEN
        write(99,*) 'SDRAG: ITime,I,J,PIJU,X,RHO,CDN,U(I,J,L),V(I,J,L)',
     *   ITime,I,J,PIJU,X,RHO,CDN,U(I,J,L),V(I,J,L),
     *   ' If problem persists, winds are too high! ',
     *   'Try setting XCDLM smaller.'
        X=1.
      END IF
         AJL(J,L,52) = AJL(J,L,52)-U(I,J,L)*X
c        IF(L.EQ.LM) AIJ(I,J,IJ_SDRAG)=AIJ(I,J,IJ_SDRAG)-U(I,J,L)*X
         AIJ(I,J,IJ_SDRAG)=AIJ(I,J,IJ_SDRAG)-U(I,J,L)*X
      U(I,J,L)=U(I,J,L)*(1.-X)
      V(I,J,L)=V(I,J,L)*(1.-X)
  100 I=IP1
      END DO
      RETURN
      END
