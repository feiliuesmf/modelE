      MODULE OCEAN
!@sum  OCEAN contains all the ocean (incl. sea ice temporarily) subroutines
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
!@cont OSTRUC,OCLIM,OCLIM0    (soon ,PRECOO,OSOURC,ODIFF,etc.)
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,rhow,rhoi,shw,shi
      USE E001M12_COM, only : im,jm,lm,gdata,focean,flake,fland,fearth
     *     ,flice,Z1O,Z12O,TAU,KOCEAN,JDATE,JDAY,MONTH
      USE PBLCOM
     &     , only : npbl=>n,uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,
     *     cq=>cqgs,ipbl,bldata
      USE GEOM

      IMPLICIT NONE

      integer, parameter :: lmom = 9 ! for deep ocean diffusion
      REAL*8, SAVE,DIMENSION(IM,JM,12) :: OA
c      COMMON/WORKO/OA(IM,JM,12)

      REAL*8, PARAMETER :: XSI1=0.5, XSI2=0.5, XSI3=0.5, XSI4=0.5
      REAL*8, PARAMETER :: R1=1./XSI1, R2=1./XSI2, R3=1./XSI3,
     *     R4=1./XSI4
      REAL*8, PARAMETER :: TTRUNC = 0 !@var TTRUNC why does this exist?
!@var Z1I thickness of first layer ice (m)
      REAL*8, PARAMETER :: Z1I = .1 
!@var ACE1I ice mass first layer (kg/m^2)
      REAL*8, PARAMETER :: ACE1I = Z1I*RHOI 
!@var Z2OIM thickness of 2nd layer ice (m)
      REAL*8, PARAMETER :: Z2OIM = .4  
!@var AC2OIM ice mass 2nd layer (kg/m^2)
      REAL*8, PARAMETER :: AC2OIM = Z2OIM*RHOI
!@var TFO temperature of freezing ocean (C)
      REAL*8, PARAMETER :: TFO = -1.8

c      REAL*8, SAVE,DIMENSION(IM,JM) :: TSO
c      REAL*8, DIMENSION(IM,JM) :: RSI
c      REAL*8, DIMENSION(IM,JM,2) :: MSI
c      REAL*8, DIMENSION(IM,JM,LMSI) :: HSI

      REAL*8, SAVE,DIMENSION(IM,JM,5) :: ODATA

!@var OTA,OTB,OTC ocean heat transport coefficients
      REAL*8, SAVE,DIMENSION(IM,JM,4) :: OTA,OTB
      REAL*8, SAVE,DIMENSION(IM,JM) :: OTC

      
      REAL*8, PRIVATE,DIMENSION(IM,JM) :: DM,AOST,EOST1,EOST0,BOST
     *     ,COST,ARSI,ERSI1,ERSI0,BRSI,CRSI 
      INTEGER, PRIVATE,DIMENSION(IM,JM) ::  KRSI
      COMMON/OOBS/DM,AOST,EOST1,EOST0,BOST,COST,ARSI,ERSI1,ERSI0,BRSI
     *     ,CRSI,KRSI

      CONTAINS

      SUBROUTINE OSTRUC
!@sum  OSTRUC restructures the ocean temperature profile as ML 
!@sum         depths are changed (generally once a day)
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
      IMPLICIT NONE
C*
      INTEGER I,J,IMAX
C*
      REAL*8 MSI1, MSI2, MELT,DRSI,ROICEN,FHSI4,FHSI3
      REAL*8 E_SIDE,E_BOTTOM,GAMMA,H_C,H_ICE,PWATER,POICE,ENRG,ENRGW
     *     ,EW_FREZ,WTRO,WTRW,ENRGI,HSI4,HSI3,HSI2,HSI1,TG1,TG2,TG3,TG4
     *     ,SNOW,ACE,ROICE,TGW,TGAVE,DWTRO,WTR1O,WTR2O
      REAL*8 Z1OOLD
C*
      COMMON/WORK2/Z1OOLD(IM,JM)
C****
C**** FLAND     LAND COVERAGE (1)
C****
C**** ODATA  1  OCEAN TEMPERATURE OF FIRST LAYER (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****        4  MEAN OCEAN TEMPERATURE OF SECOND LAYER (C)
C****        5  OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)
C****
C****
C**** RESTRUCTURE OCEAN LAYERS
C****
      DO 200 J=1,JM
      IMAX=IMAXJ(J)
      DO 200 I=1,IMAX
      IF (FLAND(I,J).GE.1.) GO TO 200
      IF (Z1OOLD(I,J).GE.Z12O(I,J)) GO TO 140
      IF (Z1O(I,J).EQ.Z1OOLD(I,J)) GO TO 200
      WTR1O=RHOW*Z1O(I,J)-ODATA(I,J,2)*(GDATA(I,J,1)+ACE1I+ODATA(I,J,3))
      DWTRO=RHOW*(Z1O(I,J)-Z1OOLD(I,J))
      WTR2O=RHOW*(Z12O(I,J)-Z1O(I,J))
      IF (DWTRO.GT.0.) GO TO 120
C**** MIX LAYER DEPTH IS GETTING SHALLOWER
      ODATA(I,J,4)=ODATA(I,J,4)
     *  +((ODATA(I,J,4)-ODATA(I,J,1))*DWTRO/WTR2O+TTRUNC)
      GO TO 200
C**** MIX LAYER DEPTH IS GETTING DEEPER
  120 TGAVE=(ODATA(I,J,4)*DWTRO+(2.*ODATA(I,J,4)-ODATA(I,J,5))*WTR2O)
     *  /(WTR2O+DWTRO)
      ODATA(I,J,1)=ODATA(I,J,1)+((TGAVE-ODATA(I,J,1))*DWTRO/WTR1O
     *  +TTRUNC)
      IF (Z1O(I,J).GE.Z12O(I,J)) GO TO 140
      ODATA(I,J,4)=ODATA(I,J,4)
     *  +((ODATA(I,J,5)-ODATA(I,J,4))*DWTRO/(WTR2O+DWTRO)+TTRUNC)
      GO TO 200
C**** MIXED LAYER DEPTH IS AT ITS MAXIMUM OR TEMP PROFILE IS UNIFORM
  140 ODATA(I,J,4)=ODATA(I,J,1)
      ODATA(I,J,5)=ODATA(I,J,1)
  200 CONTINUE
C****
C**** REDUCE THE HORIZONTAL EXTENT OF ICE IF OCEAN TEMPERATURE IS WARM
C****
      DO 300 J=1,JM
      IMAX=IMAXJ(J)
      DO 300 I=1,IMAX
      IF (FLAKE(I,J) .GT. 0.) GO TO 300 ! no melting in lakes
      IF (ODATA(I,J,2).LE.0.) GO TO 300
      IF (FLAND(I,J).GE.1.) GO TO 300
C**** REDUCE ICE EXTENT IF OCEAN TEMPERATURE IS GREATER THAN TFO
      IF (ODATA(I,J,1).LE.0.) GO TO 300
      TGW=ODATA(I,J,1)
      ROICE=ODATA(I,J,2)
      MSI2=ODATA(I,J,3)
      ACE=GDATA(I,J,1)+ACE1I+MSI2
C*
        SNOW=GDATA(I,J,1) ! snow mass
        MSI1 = SNOW + ACE1I
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
        HSI4 = (SHI*TG4-LHM)*XSI4*MSI2
C*
      ENRGI = HSI1+HSI2+HSI3+HSI4 ! energy is sea ice [J/m^2]
      WTRO=Z1O(I,J)*RHOW
      WTRW=WTRO-ROICE*ACE
      EW_FREZ = WTRW*TFO*SHW ! energy of water at freezing temperature
      ENRGW=WTRW*TGW*SHW
      ENRG = ENRGW ! energy of water
      IF (ROICE*ENRGI+ENRG.LT.0.) GO TO 230
C**** THE WARM OCEAN MELTS ALL THE SNOW AND ICE
      ODATA(I,J,1)=(ROICE*ENRGI+ENRGW)/(WTRO*SHW)
      GO TO 270
C*
C**** THE WARM OCEAN COOLS TO 0 DEGREES MELTING SOME SNOW AND ICE
C**** Reduce the ice depth
  230 ODATA(I,J,1)=0. ! ODATA(I,J,1)=(ENRGW-ENRG)/(WTRO*SHW)
C*
      PWATER = 1. - FLAND(I,J) ! water fraction
      POICE = ROICE*PWATER       ! sea ice fraction
      H_C = 1.                   ! critical ice thickness [m]
      H_ICE = ACE/RHOI           ! sea ice thickness [m]
C*
      IF (H_ICE .GT. 1.) THEN
         GAMMA = H_ICE/H_C
      ELSE
         GAMMA = 1.
      END IF
C*
      E_BOTTOM = ENRG*(POICE/(1.+GAMMA)) ! for bottom melting
      E_SIDE = ENRG*(PWATER+(POICE*GAMMA)/(1.+GAMMA)) ! for side melt.
      E_SIDE = ENRG-E_BOTTOM
C*
C***  MELT ICE VERTICALLY, AND THEN HORIZONTALLY
C*
        MELT = -XSI4*MSI2*ENRG/HSI4 ! melted ice at the bottom
        IF (MSI2-MELT .LT. AC2OIM) MELT = MSI2-AC2OIM
C       FMSI3 = XSI3*MELT ! > 0.
        FHSI3 = HSI3*MELT/MSI2
        FHSI4 = HSI4*MELT*R4/MSI2
C*
C***  MELT SOME ICE HORIZONTALLY
C*
        DRSI = (ENRG+ROICE*FHSI4)/(HSI1+HSI2+HSI3+HSI4-FHSI4) ! < 0.
C*
        ROICEN = ROICE+DRSI ! new sea ice concentration
c       SNOW = SNOW*(ROICE/ROICEN) ! new snow ice mass
        MSI2 = MSI2-MELT
        HSI3 = HSI3-FHSI3
        HSI4 = HSI4+FHSI3-FHSI4
C*
C***  CONVERT SEA ICE ENTHALPY MINUS LATENT HEAT INTO TEMPERATURE
C*
        TG1 = (HSI1/(XSI1*MSI1) +LHM)/SHI ! temperatute of layer 1
        TG2 = (HSI2/(XSI2*MSI1) +LHM)/SHI ! temperatute of layer 2
        TG3 = (HSI3/(XSI3*MSI2) +LHM)/SHI ! temperatute of layer 3
        TG4 = (HSI4/(XSI4*MSI2) +LHM)/SHI ! temperatute of layer 4
C*
C**** RESAVE PROGNOSTIC QUANTITIES
C*
        ODATA(I,J,2)=ROICEN
        ODATA(I,J,3)=MSI2
C*
        GDATA(I,J,1)  = SNOW
        GDATA(I,J,3)  = TG1
        GDATA(I,J,7)  = TG2
        GDATA(I,J,15) = TG3
        GDATA(I,J,16) = TG4
C*
      GO TO 300
  270 ODATA(I,J,2)=0.
      ODATA(I,J,3)=0.
      GDATA(I,J,1)=0.
      GDATA(I,J,3)=0.
      GDATA(I,J,7)=0.
      GDATA(I,J,15) = 0.
      GDATA(I,J,16) = 0.
C*
  300 CONTINUE
      RETURN
      END SUBROUTINE OSTRUC

      SUBROUTINE OCLIM(DOZ1O)
!@sum  OCLIM calculates daily ocean data from ocean/sea ice climatologies
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean or fixed SST/fixed lakes)
      IMPLICIT NONE

      REAL*8 Z1OOLD,XO,XZO
      COMMON/WORK2/Z1OOLD(IM,JM),XO(IM,JM,3),XZO(IM,JM)

      INTEGER, SAVE :: MONTHO = 0

      INTEGER n,MD,J,I,LSTMON,K,MDMAX,IMAX
      REAL*8 PLICEN,PLICE,POICE,POCEAN,RSICSQ,ZIMIN,ZIMAX,DOZ1O,X1
     *     ,X2,Z1OMIN,RSINEW,TIME 
      INTEGER :: JDOFM(13) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)

C****
C**** OOBS AOST     monthly Average Ocean Surface Temperature
C****      BOST     1st order Moment of Ocean Surface Temperature
C****      COST     2nd order Moment of Ocean Surface Temperature
C****      EOST0/1  Ocean Surface Temperature at beginning/end of month
C****
C****      ARSI     monthly Average of Ratio of Sea Ice to Water
C****      BRSI     1st order Moment of Ratio of Sea Ice to Water
C****        or     .5*((ERSI0-limit)**2+(ERSI1-limit)**2)/(ARSI-limit)
C****      CRSI     2nd order Moment of Ratio of Sea Ice to Water
C****      ERSI0/1  Ratio of Sea Ice to Water at beginning/end of month
C****      KRSI     -1: continuous piecewise linear fit at lower limit
C****                0: quadratic fit
C****                1: continuous piecewise linear fit at upper limit
C****
C****
C**** CALCULATE DAILY OCEAN DATA FROM CLIMATOLOGY
C****
C**** ODATA  1  OCEAN TEMPERATURE (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
C**** READ IN OBSERVED OCEAN DATA
      IF (MONTH.EQ.MONTHO) GO TO 400
      IF (MONTHO.EQ.0) THEN
C****    READ IN LAST MONTH'S END-OF-MONTH DATA
         LSTMON=MONTH-1
         IF(LSTMON.EQ.0) LSTMON=12
         CALL READT (15,IM*JM,EOST0,IM*JM,EOST0,LSTMON)
         CALL READT (17,IM*JM,ERSI0,IM*JM,ERSI0,LSTMON)
      ELSE
C****    COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
         DO 320 I=1,IM*JM
         EOST0(I,1)=EOST1(I,1)
  320    ERSI0(I,1)=ERSI1(I,1)
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      MONTHO=MONTH
      IF (MONTH.EQ.1) THEN
         REWIND 15
         REWIND 17
         READ (17)         ! skip over DM-record
      END IF
      CALL READT (15,0,AOST,2*IM*JM,AOST,1) ! READS AOST,EOST1
      CALL READT (17,0,ARSI,2*IM*JM,ARSI,1) ! READS ARSI,ERSI1
C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO 330 J=1,JM
      IMAX=IMAXJ(J)
      DO 330 I=1,IMAX
      BOST(I,J)=EOST1(I,J)-EOST0(I,J)
      COST(I,J)=3.*(EOST1(I,J)+EOST0(I,J)) - 6.*AOST(I,J)
      BRSI(I,J)=0.
      CRSI(I,J)=0.                                   ! extreme cases
      KRSI(I,J)=0                                    ! ice constant
      IF(ARSI(I,J).LE.0.) GO TO 330
      IF(ARSI(I,J).GT.1.) GO TO 330
      BRSI(I,J)=ERSI1(I,J)-ERSI0(I,J)                ! quadratic fit
      CRSI(I,J)=3.*(ERSI1(I,J)+ERSI0(I,J)) - 6.*ARSI(I,J)
      MD=MDMAX+JDATE-16
      IF(ABS(CRSI(I,J)) .GT. ABS(BRSI(I,J))) THEN    ! linear fits
        RSICSQ=CRSI(I,J)*
     *  (ARSI(I,J)*CRSI(I,J) - .25*BRSI(I,J)**2 - CRSI(I,J)**2/12.)
        IF(RSICSQ.lt.0.)  then
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
          KRSI(I,J) = -1
          BRSI(I,J) = .5*(ERSI0(I,J)**2 + ERSI1(I,J)**2) / ARSI(I,J)
        ELSEIF(RSICSQ.gt.CRSI(I,J)**2)  then
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
          KRSI(I,J) = 1
          BRSI(I,J) = .5*((ERSI0(I,J)-1.)**2 + (ERSI1(I,J)-1.)**2) /
     /                (ARSI(I,J)-1.)
        END IF
      END IF
  330 CONTINUE
C**** Calculate OST, RSI and MSI for current day
  400 TIME=(JDATE-.5)/(JDOFM(MONTH+1)-JDOFM(MONTH))-.5 ! -.5 < TIME < .5
      DO 450 J=1,JM
      ZIMIN=.5
      ZIMAX=2.
      IF(J.GT.JM/2) ZIMAX=3.5
      IMAX=IMAXJ(J)
      DO 450 I=1,IMAX
      IF(FLAND(I,J).GE.1.)  GO TO 450
      IF (KOCEAN.EQ.1 .AND. FLAKE(I,J).LE.0.) GO TO 450
C**** OST always uses quadratic fit
      IF(KOCEAN.EQ.0) ODATA(I,J,1) = AOST(I,J) + BOST(I,J)*TIME +
     +                               COST(I,J)*(TIME*TIME - 1./12.d0)
      IF(KRSI(I,J)) 410,430,420
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
  410 IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .gt. 0.)  then
        RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5)  !  TIME < T0
      ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .gt. 0.)  then
        RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME)  !  T1 < TIME
      ELSE
        RSINEW = 0.                                !  T0 < TIME < T1
          endif
      GO TO 440
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
  420 IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .lt. 1.)  then
        RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5)  !  TIME < T0
      ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .lt. 1.)  then
        RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME)  !  T1 < TIME
      ELSE
        RSINEW = 1.                                !  T0 < TIME < T1
          endif
      GO TO 440
C**** RSI uses quadratic fit
  430 RSINEW = ARSI(I,J) + BRSI(I,J)*TIME +
     +         CRSI(I,J)*(TIME*TIME - 1./12.d0)
  440 ODATA(I,J,2)=RSINEW
      ODATA(I,J,3)=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONABLE VALUE
      IF (KOCEAN.EQ.0.AND.ODATA(I,J,1).LT.-1.8) ODATA(I,J,1)=-1.8
C**** REDUCE THE RATIO OF OCEAN ICE TO WATER BY .1*RHOI/ACEOI
c     IF (ODATA(I,J,2).GT.0.) THEN
c        BYZICE=RHOI/(Z1I*RHOI+ODATA(I,J,3))
c        ODATA(I,J,2)=ODATA(I,J,2)*(1.-.06*(BYZICE-1./5.))
c     END IF
C**** ZERO OUT SNOWOI, TG1OI, TG2OI IF THERE IS NO OCEAN ICE
      IF (ODATA(I,J,2).GT.0.) GO TO 450
      GDATA(I,J,1)=0.
      GDATA(I,J,3)=0.
      GDATA(I,J,7)=0.
      GDATA(I,J,15)=0.
      GDATA(I,J,16)=0.
  450 CONTINUE
C**** REPLICATE VALUES AT POLE
      DO 460 I=2,IM
      GDATA(I,JM,1)=GDATA(1,JM,1)
      GDATA(I,JM,3)=GDATA(1,JM,3)
      GDATA(I,JM,7)=GDATA(1,JM,7)
      GDATA(I,JM,15)=GDATA(1,JM,15)
      GDATA(I,JM,16)=GDATA(1,JM,16)
      DO 460 K=1,3
  460 ODATA(I,JM,K)=ODATA(1,JM,K)
      IF (KOCEAN.EQ.0) RETURN
C****
C**** CALCULATE DAILY OCEAN MIXED LAYER DEPTHS FROM CLIMATOLOGY
C****
C**** SAVE PREVIOUS DAY'S MIXED LAYER DEPTH IN WORK2
  500 DO 510 J=1,JM
      DO 510 I=1,IM
CORR  NEXT LINE NOT NEEDED IF Z1O WERE PART OF THE RESTART FILE
      Z1O(I,J)=BLDATA(I,J,5)
  510 Z1OOLD(I,J)=Z1O(I,J)
C*** COMPUTE Z1O ONLY AT THE BEGINNING OF A DAY (OR AT TAUI)
      IF (DOZ1O.EQ.0.) RETURN
C**** READ IN TWO MONTHS OF OCEAN DATA
      X1=1.
      DO 515 J=1,JM
      DO 515 I=1,IM
  515 XZO(I,J)=0.
      IF (JDAY.GE.16) GO TO 520
      MD=JDATE+15
      GO TO 530
  520 IF (JDAY.LE.350) GO TO 550
      MD=JDATE-16
  530 CALL READT (13,0,XZO,IM*JM,XZO,1)
      MDMAX=31
C     DO 540 MX=1,10
C 540 READ (13) M
      CALL READT (13,0,Z1O,IM*JM,Z1O,11)
      GO TO 600
C 550 DO 560 MX=1,12
  550 CALL READT (13,0,Z1O,IM*JM,Z1O,MONTH-1)
C     IF (M.EQ.MONTH) GO TO 570
      IF (JDATE.LT.16) GO TO 580
      IF (MONTH.GT.1) CALL READT (13,0,Z1O,IM*JM,Z1O,1)
C     STOP 'OCEAN FILE ERROR: MLD NOT FOUND FOR CURRENT MONTH'
  570 IF (JDATE.EQ.16) GO TO 601
      MDMAX=JDOFM(MONTH+1)-JDOFM(MONTH)
      MD=JDATE-16
      GO TO 590
  580 MDMAX=JDOFM(MONTH)-JDOFM(MONTH-1)
      MD=MDMAX+JDATE-16
  590 CALL READT (13,0,XZO,IM*JM,XZO,1)
C**** INTERPOLATE OCEAN DATA TO CURRENT DAY
  600 X1=DFLOAT(MDMAX-MD)/MDMAX
  601 X2=1.-X1
      DO 610 J=1,JM
      DO 610 I=1,IM
      Z1O(I,J)=X1*Z1O(I,J)+X2*XZO(I,J)
      IF (ODATA(I,J,2)*(1.-FLAND(I,J)).LE.0.) GO TO 610
C**** MIXED LAYER DEPTH IS INCREASED TO OCEAN ICE DEPTH + 1 METER
      Z1OMIN=1. +  .09166+.001*(GDATA(I,J,1)+ODATA(I,J,3))
      IF (Z1O(I,J).GE.Z1OMIN) GO TO 605
      WRITE(6,602) TAU,I,J,MONTH,Z1O(I,J),Z1OMIN
  602 FORMAT (' INCREASE OF MIXED LAYER DEPTH ',F9.0,3I4,2F10.3)
      Z1O(I,J)=Z1OMIN
  605 IF (Z1OMIN.LE.Z12O(I,J)) GO TO 610
C**** ICE DEPTH+1>MAX MIXED LAYER DEPTH : CHANGE OCEAN TO LAND ICE
      PLICE=FLICE(I,J)
      PLICEN=1.-FEARTH(I,J)
      POICE=(1.-FLAND(I,J))*ODATA(I,J,2)
      POCEAN=(1.-FLAND(I,J))*(1.-ODATA(I,J,2))
      GDATA(I,J,12)=(GDATA(I,J,12)*PLICE+GDATA(I,J,1)*POICE)/PLICEN
      GDATA(I,J,13)=(GDATA(I,J,13)*PLICE+
     *  (GDATA(I,J,3)*XSI1+GDATA(I,J,7)*XSI2)*POICE+
     *  (LHM+SHW*ODATA(I,J,1))*POCEAN/SHI)/PLICEN
      GDATA(I,J,14)=(GDATA(I,J,14)*PLICE+
     *  (GDATA(I,J,15)*XSI3+GDATA(I,J,16)*XSI4)*POICE+
     *  (LHM+SHW*ODATA(I,J,1))*POCEAN/SHI)/PLICEN
      FLAND(I,J)=1.
      FLICE(I,J)=PLICEN
C**** MARK THE POINT FOR RESTART PURPOSES
      GDATA(I,J,1)=-10000.-GDATA(I,J,1)
C**** Transfer PBL-quantities
      ipbl(i,j,1)=0
      ipbl(i,j,2)=0
      ipbl(i,j,3)=1
      do n=1,npbl
      uabl(n,i,j,3)=(uabl(n,i,j,1)*pocean + uabl(n,i,j,2)*poice +
     +               uabl(n,i,j,3)*plice)/plicen
      vabl(n,i,j,3)=(vabl(n,i,j,1)*pocean + vabl(n,i,j,2)*poice +
     +               vabl(n,i,j,3)*plice)/plicen
      tabl(n,i,j,3)=(tabl(n,i,j,1)*pocean + tabl(n,i,j,2)*poice +
     +               tabl(n,i,j,3)*plice)/plicen
      eabl(n,i,j,3)=(eabl(n,i,j,1)*pocean + eabl(n,i,j,2)*poice +
     +               eabl(n,i,j,3)*plice)/plicen
      end do
      cm(i,j,3)=(cm(i,j,1)*pocean + cm(i,j,2)*poice +
     +           cm(i,j,3)*plice)/plicen
      ch(i,j,3)=(ch(i,j,1)*pocean + ch(i,j,2)*poice +
     +           ch(i,j,3)*plice)/plicen
      cq(i,j,3)=(cq(i,j,1)*pocean + cq(i,j,2)*poice +
     +           cq(i,j,3)*plice)/plicen
      WRITE(6,606) 100.*POICE,100.*POCEAN,TAU,I,J
  606 FORMAT(F6.1,'% OCEAN ICE AND',F6.1,'% OPEN OCEAN WERE',
     *  ' CHANGED TO LAND ICE AT TAU,I,J',F10.1,2I4)
  610 CONTINUE
  620 REWIND 13
C**** PREVENT Z1O, THE MIXED LAYER DEPTH, FROM EXCEEDING Z12O
      DO 630 J=1,JM
      DO 630 I=1,IM
      IF (Z1O(I,J).GT.Z12O(I,J)-.01) Z1O(I,J)=Z12O(I,J)
CORR  NEXT LINE NOT NEEDED IF Z1O WERE PART OF THE RESTART FILE
      BLDATA(I,J,5)=Z1O(I,J)
  630 CONTINUE
C**** SET MARKER INDICATING BLDATA(.,.,5)=Z1O
      BLDATA(IM,1,5)=-9999.
      RETURN
      END SUBROUTINE OCLIM

      SUBROUTINE OCLIM0
!@sum  OCLIM0 reads in mixed layer depths at the start
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean or fixed SST/fixed lakes)
      IMPLICIT NONE

      CALL READT (17,0,DM,IM*JM,DM,1)

      RETURN
      END SUBROUTINE OCLIM0

      END MODULE OCEAN

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean are reasonable
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM
      USE OCEAN
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in ocean data
      CALL CHECK3(ODATA,IM,JM,5,SUBR,'od')

      END SUBROUTINE CHECKO

c      MODULE SEAICE
c!@sum  SEAICE contains all the sea ice related subroutines
c!@auth Original Development Team
c!@ver  1.0 (Q-flux ocean)
c!@cont SIMELT,PRECSI,SEAICE,etc.
c
c      REAL*8, PARAMETER :: XSI1=0.5, XSI2=0.5, XSI3=0.5, XSI4=0.5
c      REAL*8, DIMENSION(IM,JM) :: RSI
c      REAL*8, DIMENSION(IM,JM,2) :: MSI
c      REAL*8, DIMENSION(IM,JM,LMSI) :: HSI
c
c
c      CONTAINS
c
c
c      END MODULE SEAICE
