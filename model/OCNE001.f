      MODULE OCEAN
!@sum  OCEAN contains all the ocean (+ sea ice temporarily) subroutines
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
!@cont OSTRUC,OCLIM,init_OCEAN,uset_OCEAN,daily_OCEAN,PREC_OC,PREC_SI,
!@cont SEA_ICE
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,rhow,rhoi,shw,shi
      USE E001M12_COM, only : im,jm,lm,gdata,focean,flake,fland,fearth
     *     ,flice,tau,kocean,jdate,jday,month
      USE PBLCOM
     &     , only : npbl=>n,uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,
     *     cq=>cqgs,ipbl,Z1O
      USE GEOM

      IMPLICIT NONE
!@param LMOM number of layers for deep ocean diffusion
      INTEGER, PARAMETER :: LMOM = 9 
!@var OA generic array for ocean heat transport calculations
      REAL*8, SAVE,DIMENSION(IM,JM,12) :: OA

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

c      REAL*8, DIMENSION(IM,JM) :: TSO
c      REAL*8, DIMENSION(IM,JM) :: RSI
c      REAL*8, DIMENSION(IM,JM,2) :: MSI
c      REAL*8, DIMENSION(IM,JM,LMSI) :: HSI

      REAL*8, SAVE,DIMENSION(IM,JM,5) :: ODATA

!@var OTA,OTB,OTC ocean heat transport coefficients
      REAL*8, SAVE,DIMENSION(IM,JM,4) :: OTA,OTB
      REAL*8, SAVE,DIMENSION(IM,JM) :: OTC

!@var Z1O ocean mixed layer depth (temporarily in PBLCOM)
!@var Z12O annual maximum ocean mixed layer depth
      REAL*8, SAVE,DIMENSION(IM,JM) :: Z12O

      REAL*8, DIMENSION(IM,JM) :: DM,AOST,EOST1,EOST0,BOST
     *     ,COST,ARSI,ERSI1,ERSI0,BRSI,CRSI
      INTEGER, DIMENSION(IM,JM) ::  KRSI
      COMMON/OOBS/DM,AOST,EOST1,EOST0,BOST,COST,ARSI,ERSI1,ERSI0,BRSI
     *     ,CRSI,KRSI

      REAL*8, SAVE,DIMENSION(IM,JM) :: T50
!@var iu_OSST,iu_SICE,iu_OCNML unit numbers for climatologies
      INTEGER iu_OSST,iu_SICE,iu_OCNML

      CONTAINS

      SUBROUTINE OSTRUC
!@sum  OSTRUC restructures the ocean temperature profile as ML
!@sum         depths are changed (generally once a day)
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
      IMPLICIT NONE
      INTEGER I,J,IMAX
      REAL*8 MSI1, MSI2, MELT,DRSI,ROICEN,FHSI4,FHSI3
      REAL*8 E_SIDE,E_BOTTOM,GAMMA,H_C,H_ICE,PWATER,POICE,ENRG,ENRGW
     *     ,EW_FREZ,WTRO,WTRW,ENRGI,HSI4,HSI3,HSI2,HSI1,TG1,TG2,TG3,TG4
     *     ,SNOW,ACE,ROICE,TGW,TGAVE,DWTRO,WTR1O,WTR2O
      REAL*8 Z1OOLD
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
        SNOW=GDATA(I,J,1) ! snow mass
        MSI1 = SNOW + ACE1I
        TG1 = GDATA(I,J,3)  ! first layer sea ice temperature
        TG2 = GDATA(I,J,7)  ! second layer sea ice temperature
        TG3 = GDATA(I,J,15) ! third layer sea ice temperature
        TG4 = GDATA(I,J,16) ! fourth layer sea ice temperature
C***  CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
        HSI1 = (SHI*TG1-LHM)*XSI1*MSI1
        HSI2 = (SHI*TG2-LHM)*XSI2*MSI1
        HSI3 = (SHI*TG3-LHM)*XSI3*MSI2
        HSI4 = (SHI*TG4-LHM)*XSI4*MSI2
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
C**** THE WARM OCEAN COOLS TO 0 DEGREES MELTING SOME SNOW AND ICE
C**** Reduce the ice depth
  230 ODATA(I,J,1)=0. ! ODATA(I,J,1)=(ENRGW-ENRG)/(WTRO*SHW)
      PWATER = 1. - FLAND(I,J) ! water fraction
      POICE = ROICE*PWATER       ! sea ice fraction
      H_C = 1.                   ! critical ice thickness [m]
      H_ICE = ACE/RHOI           ! sea ice thickness [m]
      IF (H_ICE .GT. 1.) THEN
         GAMMA = H_ICE/H_C
      ELSE
         GAMMA = 1.
      END IF
      E_BOTTOM = ENRG*(POICE/(1.+GAMMA)) ! for bottom melting
      E_SIDE = ENRG*(PWATER+(POICE*GAMMA)/(1.+GAMMA)) ! for side melt.
      E_SIDE = ENRG-E_BOTTOM
C***  MELT ICE VERTICALLY, AND THEN HORIZONTALLY
        MELT = -XSI4*MSI2*ENRG/HSI4 ! melted ice at the bottom
        IF (MSI2-MELT .LT. AC2OIM) MELT = MSI2-AC2OIM
C       FMSI3 = XSI3*MELT ! > 0.
        FHSI3 = HSI3*MELT/MSI2
        FHSI4 = HSI4*MELT*R4/MSI2
C***  MELT SOME ICE HORIZONTALLY
        DRSI = (ENRG+ROICE*FHSI4)/(HSI1+HSI2+HSI3+HSI4-FHSI4) ! < 0.
        ROICEN = ROICE+DRSI ! new sea ice concentration
c       SNOW = SNOW*(ROICE/ROICEN) ! new snow ice mass
        MSI2 = MSI2-MELT
        HSI3 = HSI3-FHSI3
        HSI4 = HSI4+FHSI3-FHSI4
C***  CONVERT SEA ICE ENTHALPY MINUS LATENT HEAT INTO TEMPERATURE
        TG1 = (HSI1/(XSI1*MSI1) +LHM)/SHI ! temperatute of layer 1
        TG2 = (HSI2/(XSI2*MSI1) +LHM)/SHI ! temperatute of layer 2
        TG3 = (HSI3/(XSI3*MSI2) +LHM)/SHI ! temperatute of layer 3
        TG4 = (HSI4/(XSI4*MSI2) +LHM)/SHI ! temperatute of layer 4
C**** RESAVE PROGNOSTIC QUANTITIES
        ODATA(I,J,2)=ROICEN
        ODATA(I,J,3)=MSI2
        GDATA(I,J,1)  = SNOW
        GDATA(I,J,3)  = TG1
        GDATA(I,J,7)  = TG2
        GDATA(I,J,15) = TG3
        GDATA(I,J,16) = TG4
      GO TO 300
  270 ODATA(I,J,2)=0.
      ODATA(I,J,3)=0.
      GDATA(I,J,1)=0.
      GDATA(I,J,3)=0.
      GDATA(I,J,7)=0.
      GDATA(I,J,15) = 0.
      GDATA(I,J,16) = 0.
  300 CONTINUE
      RETURN
      END SUBROUTINE OSTRUC

      SUBROUTINE OCLIM(IDOZ1O)
!@sum OCLIM calculates daily ocean data from ocean/sea ice climatologies
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean or fixed SST/fixed lakes)
      IMPLICIT NONE

      REAL*8 Z1OOLD
      COMMON/WORK2/Z1OOLD(IM,JM)
      REAL*8 XZO(IM,JM),XZN(IM,JM)
      INTEGER IDOFM(0:13)
      DATA IDOFM /-15,16,47,75,106,136,167,197,228,259,289,320,350,381/
      INTEGER, INTENT(IN) :: IDOZ1O

      INTEGER n,MD,J,I,LSTMON,K,MDMAX,IMAX,IMON
      REAL*8 PLICEN,PLICE,POICE,POCEAN,RSICSQ,ZIMIN,ZIMAX,X1
     *     ,X2,Z1OMIN,RSINEW,TIME,FRAC
      INTEGER :: JDOFM(13) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
!@var JDLAST julian day that OCLIM was last called
      INTEGER, SAVE :: JDLAST=0
!@var MONTHO current month for climatology reading
      INTEGER, SAVE :: MONTHO = 0

      IF (KOCEAN.EQ.1) GO TO 500
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
         CALL READT (iu_OSST,IM*JM,EOST0,IM*JM,EOST0,LSTMON)
         CALL READT (iu_SICE,IM*JM,ERSI0,IM*JM,ERSI0,LSTMON)
      ELSE
C****    COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
         DO 320 I=1,IM*JM
         EOST0(I,1)=EOST1(I,1)
  320    ERSI0(I,1)=ERSI1(I,1)
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      MONTHO=MONTH
      IF (MONTH.EQ.1) THEN
         REWIND iu_OSST
         REWIND iu_SICE
         READ (iu_SICE)         ! skip over DM-record
      END IF
      CALL READT (iu_OSST,0,AOST,2*IM*JM,AOST,1) ! READS AOST,EOST1
      CALL READT (iu_SICE,0,ARSI,2*IM*JM,ARSI,1) ! READS ARSI,ERSI1
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
C**** OST always uses quadratic fit
      ODATA(I,J,1) = AOST(I,J) + BOST(I,J)*TIME +
     +               COST(I,J)*(TIME*TIME - 1./12.d0)
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
      IF (ODATA(I,J,1).LT.-1.8) ODATA(I,J,1)=-1.8
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
      RETURN
C****
C**** CALCULATE DAILY OCEAN MIXED LAYER DEPTHS FROM CLIMATOLOGY
C****
C**** SAVE PREVIOUS DAY'S MIXED LAYER DEPTH IN WORK2
 500  Z1OOLD=Z1O
C**** COMPUTE Z1O ONLY AT THE BEGINNING OF A DAY (OR AT TAUI)
C**** Check mark to see if Z1O needs to be set initially
      IF (IDOZ1O.eq.0.and. Z1O(IM,1).eq.-9999.) RETURN
C**** Read in climatological ocean mixed layer depths efficiently
      IF (JDLAST.eq.0) THEN ! need to read in first month climatology
        IMON=1          ! IMON=January
        IF (JDAY.LE.16)  THEN ! JDAY in Jan 1-15, first month is Dec
          CALL READT (iu_OCNML,0,XZO,IM*JM,XZO,12)
          REWIND iu_OCNML
        ELSE            ! JDAY is in Jan 16 to Dec 16, get first month
  520     IMON=IMON+1
          IF (JDAY.GT.IDOFM(IMON) .and. IMON.LE.12) GO TO 520
          CALL READT (iu_OCNML,0,XZO,IM*JM,XZO,IMON-1)
          IF (IMON.EQ.13)  REWIND iu_OCNML
        END IF
      ELSE                      ! Do we need to read in second month?
        IF (JDAY.NE.JDLAST+1) THEN ! Check that data is read in daily
          IF (JDAY.NE.1 .or. JDLAST.NE.365) THEN
            WRITE (6,*) 'Incorrect values in OCLIM: JDAY,JDLAST=',JDAY
     *           ,JDLAST
            STOP 'ERROR READING IN SETTING OCEAN CLIMATOLOGY'
          END IF
          IMON=IMON-12          ! New year
          GO TO 530
        END IF
        IF (JDAY.LE.IDOFM(IMON)) GO TO 530
        IMON=IMON+1          ! read in new month of climatological data
        XZO = XZN
        IF (IMON.EQ.13) REWIND iu_OCNML
      END IF
      CALL READT (iu_OCNML,0,XZN,IM*JM,XZN,1)
 530  JDLAST=JDAY
      REWIND iu_OCNML
C**** INTERPOLATE OCEAN DATA TO CURRENT DAY
      FRAC = DBLE(IDOFM(IMON)-JDAY)/(IDOFM(IMON)-IDOFM(IMON-1))
      DO 610 J=1,JM
      DO 610 I=1,IM
      Z1O(I,J)=FRAC*XZO(I,J)+(1.-FRAC)*XZN(I,J)

      IF (ODATA(I,J,2)*FOCEAN(I,J).le.0) GO TO 610
C**** MIXED LAYER DEPTH IS INCREASED TO OCEAN ICE DEPTH + 1 METER
      Z1OMIN=1.+(RHOI*Z1I+GDATA(I,J,1)+ODATA(I,J,3))/RHOW
      IF (Z1O(I,J).GE.Z1OMIN) GO TO 605
      WRITE(6,602) TAU,I,J,MONTH,Z1O(I,J),Z1OMIN
  602 FORMAT (' INCREASE OF MIXED LAYER DEPTH ',F9.0,3I4,2F10.3)
      Z1O(I,J)=Z1OMIN
  605 IF (Z1OMIN.LE.Z12O(I,J)) GO TO 610
C**** ICE DEPTH+1>MAX MIXED LAYER DEPTH : CHANGE OCEAN TO LAND ICE
      PLICE=FLICE(I,J)
      PLICEN=1.-FEARTH(I,J)
      POICE=FOCEAN(I,J)*ODATA(I,J,2)
      POCEAN=FOCEAN(I,J)*(1.-ODATA(I,J,2))
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
C**** PREVENT Z1O, THE MIXED LAYER DEPTH, FROM EXCEEDING Z12O
      DO 630 J=1,JM
      DO 630 I=1,IM
      IF (Z1O(I,J).GT.Z12O(I,J)-.01) Z1O(I,J)=Z12O(I,J)
  630 CONTINUE
C**** SET MARKER INDICATING THAT Z1O HAS BEEN SET
      Z1O(IM,1)=-9999.
      RETURN
      END SUBROUTINE OCLIM

      SUBROUTINE PREC_OC(TGW,WTRO,PRCP,ENRGP,ERUN4)
C****                                                                      
C**** THIS SUBROUTINE USES THE PRECIPITATION TO CALCULATE                 
C**** THE OCEAN SURFACE TEMPERATURE                                
C****                                                                
      IMPLICIT NONE
      REAL*8 TGW, PRCP, WTRO, ENRGP, ERUN4, ENRGO 

      ENRGO=ENRGP-ERUN4
      TGW=TGW+(ENRGO/(WTRO*SHW)+TTRUNC)

      RETURN 
      END SUBROUTINE PREC_OC

      SUBROUTINE PREC_SI(SNOW,MSI2,MSI1,TG1,TG2,TG3,TG4,PRCP,TPRCP
     *     ,EPRCP,RUN0,DIFS,EDIFS,ERUN2)
C****                                                                      
C**** THIS SUBROUTINE USES THE PRECIPITATION TO CALCULATE                 
C**** THE OCEAN ICE SNOW COVER AND RUNOFF                               
C****                                                                
      IMPLICIT NONE

      REAL*8, PARAMETER :: SNOMAX=100.0, dSNdRN=0.
      REAL*8 SNOW, MSI1, MSI2, MELT1 
      REAL*8 TPRCP, EPRE, EPRCP, PRCP, RAIN, FREZ1
      REAL*8 TG1, TG2, TG3, TG4, HSI1, HSI2, HSI3, HSI4, 
     *       FMSI1, FMSI2, FHSI1, FHSI2, FHSI3, H2, CMPRS
      REAL*8 DIFS, EDIFS, ERUN2 ! for diagnostics 
      REAL*8 RUN0 ! for SST calculation 
      REAL*8 HC_1 

      RUN0=0. ; DIFS=0. ; EDIFS=0. ; ERUN2=0.
C*** sensible -  latent heat of precipit.
      EPRE = EPRCP - PRCP*LHM ! total energy of precipitation;
C***  CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2
      HSI4 = (SHI*tG4-LHM)*XSI4*MSI2
      HC_1 = XSI1*MSI1*SHI
      HSI1 = HSI1+EPRE
      IF (TPRCP.LT.0.) GO TO 180
      IF (EPRCP.LT.-TG1*HC_1) GO TO 160
C***  ALL PRECIPITATION IS RAIN ABOVE 0degC
C***  RAIN COMPRESSES SNOW INTO ICE
      RAIN = PRCP
      IF (HSI1/LHM+XSI1*MSI1 .LE. 0.) GO TO 140
C***  WARM RAIN MELTS SOME SNOW OR ICE
      MELT1 = HSI1/LHM+XSI1*MSI1 ! melted snow and ice (kg/m^2)
      RUN0=MELT1 + PRCP ! water mass flux to ocean (kg/m^2)
      IF (MSI1-ACE1I .LE. MELT1) GO TO 130
C**** RAIN MELTS SOME SNOW AND COMPRESSES SNOW INTO ICE
      FMSI2 = MIN (dSNdRN*(RAIN+MELT1), MSI1-ACE1I-MELT1) ! > 0.
      FHSI1 = -LHM*(XSI1*FMSI2-XSI2*MELT1) ! downward heat flux
c     F1 = HSI1*(XSI1*FMSI2-XSI2*MELT1)/(XSI1*MSI1-MELT1)
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! > 0.
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
      ERUN2 = HSI3*FMSI2*R3/MSI2
      IF (KOCEAN.NE.1) GO TO 380
      GO TO 210
  130 CONTINUE
C**** RAIN MELTS ALL SNOW AND SOME ICE
      FMSI2 = MSI1-ACE1I-MELT1 ! < 0.(upward ice mass flux)
      DIFS = FMSI2 ! < 0.
C     FMSI1 = XSI1*(MSI1-ACE1I)-MELT1 ! < 0.(melted snow/ice mass)
      FHSI1 = HSI2*((XSI1/XSI2)*FMSI2-MELT1)/MSI1 !  upward heat flux
      FHSI2 = HSI3*FMSI2*R3/MSI2 !  upward heat flux into layer 2
      FHSI3 = HSI4*FMSI2/MSI2 !  upward heat flux into layer 3
      SNOW=0. ! Rain melted all snow
      MSI1 = ACE1I ! Keep the first layer ice mass constant ACE1I
      HSI1 = HSI1-FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      EDIFS=FHSI2 ! for diagnostics
      ERUN2=EDIFS ! for diagnostics
      IF (KOCEAN.NE.1) GO TO 380
      GO TO 210
  140 CONTINUE
C**** RAIN COMPRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      CMPRS = MIN(dSNdRN*RAIN, MSI1-ACE1I) ! all snow or part of rain
      IF (-HSI1/LHM-XSI1*MSI1 .LT. RAIN) GO TO 150
C***  ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
C     FMSI1 = XSI1*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      RUN0 = 0.
      FHSI1 = HSI1*(XSI1*CMPRS+RAIN)/(XSI1*MSI1+RAIN) ! downward
C      F1 = HSI1*(XSI1*CMPRS+RAIN)/(XSI1*MSI1) ! for prescribed ice
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
      IF (KOCEAN.NE.1) GO TO 380
      GO TO 210
  150 CONTINUE
C**** JUST PART OF RAIN FREEZES IN LAYER 1
      FREZ1 = -HSI1/LHM-XSI1*MSI1 ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI1*CMPRS+FREZ1 ! downward ice mass flux from layer 1
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI1*CMPRS+FREZ1) ! downward heat flux from layer 1
c      F1 = HSI1*(XSI1*CMPRS+FREZ1)/(XSI1*MSI1) ! for prescribed ice
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
      IF (KOCEAN.NE.1) GO TO 380
      GO TO 210
  160 CONTINUE
C***  PRECIPITATION IS A MIXTURE OF RAIN AND SNOW AT 0 degC
C***  RAIN COMRESSES SNOW INTO ICE, SOME RAIN WILL FREEZE
      SNOW = -EPRCP/LHM ! snow fall
      RAIN = PRCP-SNOW  ! rain fall
      CMPRS = MIN(dSNdRN*RAIN, MSI1+SNOW-ACE1I) ! compression
      IF (-HSI1/LHM-XSI1*MSI1-SNOW .LT. RAIN) GO TO 170
C***  ALL RAIN FREEZES IN LAYER 1
C     FREZ1 = RAIN ! frozen rain
      RUN0 =0.
      FMSI1 = XSI1*CMPRS+XSI2*SNOW+RAIN ! downward ice mass flux
      FMSI2 = CMPRS+RAIN ! downward ice mass flux from layer 2
      FHSI1 = HSI1*FMSI1/(XSI1*MSI1+PRCP) ! downward heat flux
c      F1 = HSI1*FMSI1/(XSI1*MSI1) ! for prescribed ice/ocean
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
      IF (KOCEAN .NE. 1) GO TO 380
      GO TO 210
  170 CONTINUE
C***  NOT ALL RAIN FREEZES IN LAYER 1
      FREZ1 = -HSI1/LHM-XSI1*MSI1-SNOW ! part of rain that freezes
      RUN0 = RAIN-FREZ1 ! water mass flux into the ocean
C     FMSI1 = XSI1*CMPRS+XSI2*SNOW+FREZ1 ! downward ice mass flux
      FMSI2 = CMPRS+FREZ1 ! downward ice mass flux from layer 2
      FHSI1 = -LHM*(XSI1*CMPRS+XSI2*SNOW+FREZ1) ! downward
c      F1 = HSI1*(XSI1*CMPRS+XSI2*SNOW+FREZ1)/(XSI1*MSI1) ! prescribed
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
      IF (KOCEAN .NE. 1) GO TO 380
      GO TO 210
  180 CONTINUE
C***  ALL PRECIPITATION IS SNOW, SNOW AMOUNT INCREASES
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
      GO TO 380
  190 CONTINUE
C***  TOO MUCH SNOW HAS ACCUMULATED, SOME SNOW IS COMPACTED INTO ICE
      FMSI2 = MSI1+PRCP-(0.9*SNOMAX+ACE1I) ! > 0.(compressed snow)
C     FMSI1 = XSI2*PRCP+XSI1*FMSI2 ! > 0.(downward ice mass flux)
      FHSI1 = HSI1*(XSI2*PRCP+XSI1*FMSI2)/(XSI1*MSI1+PRCP) ! downward
c      F1 = HSI1*(XSI2*PRCP+XSI1*FMSI2)/(XSI1*MSI1) ! for prescribed
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux from layer 2
      FHSI3 = HSI3*FMSI2*(XSI4/XSI3)/MSI2 ! downward heat flux
      MSI1 = 0.9*SNOMAX+ACE1I ! first layer ice mass
      SNOW=.9*SNOMAX ! snow mass
      DIFS = FMSI2 ! compressed snow mass
      EDIFS = FHSI2 ! energy of compressed snow mass
      HSI1 = HSI1-FHSI1
      TG1 = (HSI1/(XSI1*MSI1)+LHM)/SHI ! first layer ice temperature
      IF (KOCEAN .EQ. 1) GO TO 210
      ERUN2 = HSI3*FMSI2*R3/MSI2
      H2 = HSI2+(FHSI1-FHSI2) ! for prescribed ice/ocean
      TG2 = (H2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      GO TO 380
  210 CONTINUE
C**** ADVECT ICE (usually downwards)
      HSI2 = HSI2+(FHSI1-FHSI2) ! for predicted ice/ocean
      HSI3 = HSI3+(FHSI2-FHSI3) ! for predicted ice/ocean
      HSI4 = HSI4+FHSI3 ! for predicted ice/ocean
      TG2 = (HSI2/(XSI2*MSI1)+LHM)/SHI ! second layer ice temperature
      TG3 = (HSI3/(XSI3*MSI2)+LHM)/SHI ! third layer ice temperature
      TG4 = (HSI4/(XSI4*MSI2)+LHM)/SHI ! fourth layer ice temperature
      MSI2 = MSI2+FMSI2 ! second layer sea ice mass (kg/m^2)
 380  CONTINUE

      RETURN 
      END SUBROUTINE PREC_SI

      SUBROUTINE SEA_ICE(DTSRCE,SNOW,ROICE,PLAKE,TG1,TG2,TG3,TG4,
     *     MSI2,F0DT,F1DT,EVAP,TGW,WTRO,OTDT,ENRGO,RUN0,RUN4,ERUN4,
     *     DIFSI,EDIFSI,DIFS,EDIFS,ENRGFO,ACEFO,ACE2M,ACE2F,ENRGFI,F2DT)
C****                                                                      
C**** THIS SUBROUTINE USES THE SURFACE FLUXES TO PREDICT IN TIME THE       
C**** SEA ICE TEMPERATURE, SEA ICE AND SNOW GROWTH/MELTING,          
C**** AND FREEZING ON THE OPEN OCEAN AND BELOW SEA ICE 
C****                                                                     
      IMPLICIT NONE

      REAL*8, PARAMETER :: ALAMI=2.1762, ALAMS=0.35, ! J/(m*degC*sec)
     *        RHOS = 300.0, ! kg/m^3 (snow density) 
     *        BYRLI = 1./(RHOI*ALAMI), 
     *        BYRLS = 1./(RHOS*ALAMS), ! (m^4*degC*sec)/(J*kg)
     *        ALPHA = 1.0, dSNdML =0. 
      REAL*8, PARAMETER :: YSI1 = XSI1*ACE1I/(ACE1I+AC2OIM),
     *                     YSI2 = XSI2*ACE1I/(ACE1I+AC2OIM),
     *                     YSI3 = XSI3*AC2OIM/(ACE1I+AC2OIM), 
     *                     YSI4 = XSI4*AC2OIM/(ACE1I+AC2OIM)
      REAL*8, PARAMETER :: Z2OIX = 4.9, 
     *                     BYZICX=1./(Z1I+Z2OIX) 
      REAL*8,  INTENT(IN) :: DTSRCE
!@var F0DT heat flux on the ice top surface (W/m^2)
      REAL*8,  INTENT(IN) :: F0DT
!@var F1DT heat flux between the 1st and 2nd ice layers (W/m^2)
      REAL*8,  INTENT(IN) :: F1DT
!@var EVAP evaporation/dew on the top ice surface (kg/m^2)
      REAL*8,  INTENT(IN) :: EVAP 
      REAL*8 PLAKE, ROICE, ROICEN, OPNOCN,  
     *       SNOW, MSI1, MSI2, ACE2F, ACE2M, MELT1, MELT4,  
     *       DEW, CMPRS, DRSI, DRI
      REAL*8 TG1, TG2, TG3, TG4, HSI1, HSI2, HSI3, HSI4, 
     *       FMSI1, FMSI2, FMSI3, FMSI4, 
     *       FHSI1, FHSI2, FHSI3, FHSI4  
!@var TGW and WTRO mixed layer temp.(C) and mass (kg/m^2)  
      REAL*8 TGW, WTRO, OTDT, ENRGO
      REAL*8 WTRI0, EIW0, ENRGIW, ENRGFI, WTRI1, EFIW, 
     *       ENRGO0, EOFRZ, ENRGFO, ACEFO
      REAL*8 WTRW0, ENRGW0, WTRW, ENRGW ! only for ocean 
      REAL*8 HC1, HC2, HC3, HC4 
      REAL*8 dF1dTI, dF2dTI, dF3dTI, dF4dTI, F1, F2, F3
      REAL*8, INTENT(OUT) :: EDIFSI, ERUN4, F2DT ! for diagnostics 
      REAL*8, INTENT(OUT) :: RUN0, RUN4, DIFSI, DIFS, EDIFS 
C**** Initiallise output 
      F2DT=0. ; RUN0=0.  ; DIFSI=0. ; EDIFSI=0.
      DIFS=0. ; EDIFS=0. ; ERUN4=0. ; RUN4=0.

      IF (KOCEAN .EQ. 1 ) THEN 
        ENRGO0=WTRO*TGW*SHW
        EOFRZ=WTRO*TFO*SHW
        IF (ENRGO0+ENRGO.GE.EOFRZ) THEN
C**** FLUXES RECOMPUTE TGW WHICH IS ABOVE FREEZING POINT FOR OCEAN
          ENRGFO=0.
          ACEFO=0.
          IF (ROICE.LE.0.) THEN
            TGW=TGW+(ENRGO/(WTRO*SHW)+TTRUNC)
            RETURN
          END IF
        ELSE
C**** FLUXES COOL TGO TO FREEZING POINT FOR OCEAN AND FORM SOME ICE
c 80      ACEFO=(ENRGO0+ENRGO-EOFRZ)/(TFO*(SHI-SHW)-LHM)
 80       ACEFO=(ENRGO0+ENRGO-EOFRZ)/(-LHM)
          ENRGFO=ACEFO*(TFO*SHI-LHM)
          IF (ROICE.LE.0.) THEN
            ROICE=ACEFO/(ACE1I+AC2OIM)
            TGW=TFO
            SNOW=0.
            TG1=TFO
            TG2=TFO
            TG3=TFO
            TG4=TFO
            MSI2=AC2OIM
            RETURN
          END IF
        END IF
      END IF 
C****
      IF (ROICE .EQ. 0.) RETURN 
      ACE2F=0. ! frozen ice at the bottom
      ACE2M=0. ! melted ice at the bottom 
C****
      MSI1 = SNOW+ACE1I ! snow and first (physical) layer ice mass
C***  CONVERT SEA ICE TEMPERATURE INTO ENTHALPY MINUS LATENT HEAT
      HSI1 = (SHI*TG1-LHM)*XSI1*MSI1 ! J/m^2
      HSI2 = (SHI*TG2-LHM)*XSI2*MSI1 ! J/m^2
      HSI3 = (SHI*TG3-LHM)*XSI3*MSI2 ! J/m^2
      HSI4 = (SHI*TG4-LHM)*XSI4*MSI2 ! J/m^2

      IF (KOCEAN .EQ. 1) THEN
        WTRI0=WTRO-(SNOW+ACE1I+MSI2) ! mixed layer mass below ice (kg/m^2)
        EIW0=WTRI0*TGW*SHW      ! energy of mixed layer below ice (J/m^2) 
        WTRW0=WTRO-ROICE*(SNOW+ACE1I+MSI2)
        ENRGW0=WTRW0*TGW*SHW
        RUN4=-EVAP
        WTRW0 = WTRW0-ROICE*RUN4 ! water mass "+-" dew/evaporation
        ERUN4=TGW*RUN4*SHW
      END IF 
C****
C**** OCEAN ICE, CALCULATE TG1 AND 
C****
      HC1 = SHI*XSI1*MSI1 ! heat capacity of ice layer 1 (J/(degC*m^2))
      HC2 = SHI*XSI2*MSI1 ! heat capacity of ice layer 2 (J/(degC*m^2))
      HC3 = SHI*XSI3*MSI2 ! heat capacity of ice layer 3 (J/(degC*m^2))
      HC4 = SHI*XSI4*MSI2 ! heat capacity of ice layer 4 (J/(degC*m^2))
C***  CALCULATE AND APPLY DIFFUSIVE AND SURFACE ENERGY FLUXES
c      dF1dTI = 2.*DTSRCE/(ACE1I*BYRLI+SNOW*BYRLS) ! for non-Q-flux      
      dF2dTI = ALAMI*RHOI*DTSRCE/(0.5*XSI2*MSI1+0.5*XSI3*MSI2)
C***           temperature derivative from F2 diffusive flux
      dF3dTI = ALAMI*RHOI*DTSRCE*2./MSI2
C***           temperature derivative from F3 diffusive flux
      dF4dTI = ALAMI*RHOI*DTSRCE*2.*R4/MSI2
C***           temperature derivative from F4 diffusive flux
CEXP  F2 = dF2dTI*(TG2-TG3) ! the diffusive
CEXP  F3 = dF3dTI*(TG3-TG4) ! fluxes from
CEXP  F4 = dF4dTI*(TG4-TGW) ! explicit method
C***  DIFFUSIVE FLUXES FROM IMPLICIT METHOD
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
      RUN0 = MELT1 ! water mass that flows to the ocean (kg/m^2)
      IF (KOCEAN.EQ.1) WTRW0 = WTRW0+ROICE*RUN0 ! ocean mass (kg/m^2)
C***  EVAPORATION FROM THE SURFACE
      IF (DEW .GT. 0.) GO TO 140 ! go to the dew case
C***  DEW IS EVAPORATION NOW
C***  EVAPORATION REDUCES SNOW OR ICE, MELT1>0.
      IF (MSI1-ACE1I+DEW .GT. MELT1) GO TO 130
C***  ALL SNOW AND SOME ICE EVAPORATE AND MELT
C***  ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 2 into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1 ! < 0.
C*            upward mass flux from layer 3 into layer 2
      FHSI1 = HSI2*FMSI1*R2/MSI1 ! energy of ice mass FMSI1
      FHSI2 = HSI3*FMSI2*R3/MSI2 ! energy of ice mass FMSI2
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      SNOW = 0. ! all snow evaporates and melts
      DIFS = FMSI2 ! <0 (upward mass flux into layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN.EQ.1) GO TO 170
      GO TO 370
  130 CONTINUE
C***  EVAPORATION AND MELTING REDUCE JUST SNOW AMOUNT
C***  ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
      FMSI2 = MIN(dSNdML*MELT1, MSI1-ACE1I+DEW-MELT1) ! > 0.
      FMSI1 = XSI1*FMSI2+XSI2*(DEW-MELT1)
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward from layer 2 into layer 3
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      MSI1 = MSI1+DEW-MELT1-FMSI2 ! kg/m^2
      SNOW = MSI1-ACE1I
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
      GO TO 370
  140 CONTINUE
C***  DEW INCREASES ICE AMOUNT,  MELT1 > 0.
      IF (MSI1-ACE1I .GT. MELT1) GO TO 150
C***  ALL SNOW AND SOME ICE MELT
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW-MELT1 ! upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW-MELT1
      SNOW = 0.
      IF (FMSI2 .LE. 0.) THEN ! (if melting is greater than dew)
C*** ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
        FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward into layer 1 from layer 2
        FHSI2 = HSI3*FMSI2*R3/MSI2 ! upward into layer 2 from layer 3
        HSI1 = HSI1-FHSI1
        HSI2 = HSI2+(FHSI1-FHSI2)
        DIFS = FMSI2 ! <0 upward ice mass into layer 2
        EDIFS = FHSI2 ! energy of diffused ice mass DIFS
        IF (KOCEAN .EQ. 1) GO TO 170
        GO TO 370
      ENDIF
C***  ICE ADVECTION IS DOWNWARD INTO LAYER 3 FROM LAYER 2
C***  (if dew is greater than melting)
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 2
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward into layer 3 from layer 2
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      MSI1 = ACE1I
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
      GO TO 370
  150 CONTINUE
C***  MELTING REDUCES SNOW AMOUNT
C***  ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
      CMPRS = MIN(dSNdML*MELT1, MSI1-ACE1I-MELT1) ! > 0.
      FMSI1 = DEW+XSI1*CMPRS-XSI2*MELT1
      FMSI2 = DEW+CMPRS ! > 0. downward into layer 3 from layer 2
      IF (FMSI1 .LT. 0.) FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward
      IF (FMSI1 .GE. 0.) FHSI1 = -LHM*FMSI1 ! downward into layer 3
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward energy flux into layer 3
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      MSI1 = MSI1-MELT1-CMPRS ! new mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass after melting and compression
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 1) GO TO 210
      GO TO 370
  160 CONTINUE
C***  THE FIRST LAYER IS BELOW FREEZING
C***  NO SNOW OR ICE MELTS IN LAYER 1
      IF (DEW .GT. 0.) GO TO 200 ! go to the dew case
      IF (MSI1-ACE1I+DEW .GE. 0.) GO TO 190
C***  ALL SNOW AND SOME ICE EVAPORATE
C***  ICE ADVECTION IS UPWARD INTO LAYER 2 FROM LAYER 3
      FMSI1 = XSI1*(MSI1-ACE1I)+DEW ! < 0. upward into layer 1
      FMSI2 = MSI1-ACE1I+DEW ! < 0. upward into layer 2 from layer 3
      FHSI1 = HSI2*FMSI1*R2/MSI1 ! upward energy flux into layer 1
      FHSI2 = HSI3*FMSI2*R3/MSI2 ! upward energy flux into layer 2
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      SNOW = 0. ! all snow evaporated
      DIFS = FMSI2 ! <0 upward ice mass into layer 2
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN.EQ.1) GO TO 170
      GO TO 370
  170 CONTINUE
      MSI1 = ACE1I
         DIFSI = ROICE*DIFS
         EDIFSI = ROICE*EDIFS
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 180 ! go to freezing case
C***  FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
      MELT4 = HSI4/LHM+XSI4*MSI2 ! > 0. melted ice from layer 4
      ACE2M = MELT4
      WTRW0 = WTRW0+ROICE*ACE2M
      FMSI3 = XSI4*FMSI2+XSI3*MELT4
      IF (FMSI3 .LE. 0.) FHSI3 = -LHM*FMSI3 ! upward into layer 3
      IF (FMSI3 .GT. 0.) FHSI3 = HSI3*FMSI3*R3/MSI2 ! downward
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
      MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      GO TO 230
  180 CONTINUE
C***  NO ICE MELTS IN LAYER 4
C     FMSI3 = XSI4*FMSI2 ! upward mass flux into layer 3 from layer 4
      FHSI3 = HSI4*FMSI2/MSI2 ! upward energy flux into layer 3
      HSI3 = HSI3+(FHSI2-FHSI3)
      HSI4 = HSI4+FHSI3
      MSI2 = MSI2+FMSI2
      GO TO 230
  190 CONTINUE
C***  JUST SOME SNOW EVAPORATES
C***  NO ADVECTION BETWEEN LAYER 2 AND 3
C     FMSI1 = XSI2*DEW ! < 0. upward mass flux into layer 1
C     FMSI2 = 0. ! no ice advection between layers 2 and 3
      FHSI1 = HSI2*DEW/MSI1 ! upward energy flux into layer 1
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+FHSI1
      MSI1 = MSI1+DEW ! new ice mass of the first physical layer
      SNOW = MSI1-ACE1I ! new snow mass

      IF (KOCEAN .NE. 1) GO TO 370
      IF (HSI4/LHM+XSI4*MSI2 .LE. 0.) GO TO 230 ! go to freezing case
C***  FLUXES HEAT LAYER 4 TO FREEZING POINT AND MELT SOME ICE
      MELT4 = HSI4/LHM+XSI4*MSI2 ! melted ice from layer 4
      ACE2M = MELT4
      WTRW0 = WTRW0+ROICE*ACE2M
C     FMSI3 = XSI3*MELT4 ! > 0. downward mass flux into layer 4
      FHSI3 = HSI3*MELT4/MSI2 ! downward heat flux into layer 4
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+FHSI3
      MSI2 = MSI2-MELT4 ! new ice mass of the second physical layer
      GO TO 230
  200 CONTINUE
C***  DEW INCREASES ICE AMOUNT, ADVECT ICE DOWNWARD
C***  ICE ADVECTION IS DOWNWARD FROM LAYER 2 INTO LAYER 3
C     FMSI1 = DEW ! > 0. downward mass flux into layer 2
      FMSI2 = DEW ! > 0. downward mass flux into layer 3
      FHSI1 = HSI1*FMSI2/(XSI1*MSI1+DEW) ! downward heat flux
      FHSI2 = HSI2*FMSI2*R2/MSI1 ! downward heat flux into layer 3
      HSI1 = HSI1-FHSI1
      HSI2 = HSI2+(FHSI1-FHSI2)
      DIFS = FMSI2 ! >0 (downward mass flux from layer 2)
      EDIFS = FHSI2 ! energy of diffused ice mass DIFS
      IF (KOCEAN .EQ. 0) GO TO 370
  210 CONTINUE
         DIFSI = ROICE*DIFS
         EDIFSI = ROICE*EDIFS
      IF (HSI4/LHM+XSI4*MSI2 .GT. 0.) THEN
C**** FLUXES HEAT UP TG2 TO FREEZING POINT AND MELT SOME ICE
        MELT4 = HSI4/LHM+XSI4*MSI2 ! melted ice from layer 4
        ACE2M = MELT4
        WTRW0 = WTRW0+ROICE*ACE2M
C     FMSI3 = XSI4*FMSI2+XSI3*MELT4 ! > 0. downward into layer 4
        FHSI3 = HSI3*((XSI4/XSI3)*FMSI2+MELT4)/MSI2 ! downward
        HSI3 = HSI3+(FHSI2-FHSI3)
        HSI4 = HSI4+FHSI3
        MSI2 = MSI2+FMSI2-MELT4 ! new ice mass of physical layer 2
      ELSE
C***  NO ICE MELTS IN LAYER 4
C     FMSI3 = XSI4*FMSI2 ! > 0. downward mass flux into layer 4
        FHSI3 = HSI3*(XSI4/XSI3)*FMSI2/MSI2 ! downward heat flux
        HSI3 = HSI3+(FHSI2-FHSI3)
        HSI4 = HSI4+FHSI3
        MSI2 = MSI2+FMSI2
      END IF
  230  CONTINUE
C***  CALCULATE THE ENERGY OF THE WATER BELOW THE ICE AT FREEZING POINT
C***  AND CHECK WHETHER NEW ICE MUST BE FORMED
      ENRGIW = F2DT+OTDT-ERUN4 ! heat flux to the ocean under ice
      ENRGFI = 0.
      WTRI1 = WTRO-(MSI1+MSI2) ! new mass of ocean (kg/m^2)
      EFIW = WTRI1*TFO*SHW ! freezing energy of ocean mass WTRI1
      IF (EIW0+ENRGIW .GT. EFIW) GO TO 250 ! go to no freezing case
C***  FLUXES WOULD COOL TGW TO FREEZING POINT
C***  AND FREEZE SOME MORE ICE
c      ACE2F = (EIW0+ENRGIW-EFIW)/(TFO*(SHI-SHW)-LHM) !
      ACE2F = (EIW0+ENRGIW-EFIW)/(-LHM) !
C*            ocean mass that freezes under the ice
      ENRGFI = ACE2F*(TFO*SHI-LHM) ! energy of frozen ice
C***  CALCULATE ADVECTIVE HEAT FLUX FROM LAYER 3 TO LAYER 4 OF ICE
C     FMSI3 = -XSI3*ACE2F ! < 0.
C     FMSI4 = -ACE2F
      FHSI3 = -HSI4*ACE2F*(XSI3/XSI4)/MSI2
C***  COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      IF (ACEFO .GT. 0.) GO TO 240
C***  NEW ICE IS FORMED BELOW OLD SEA ICE
      WTRW = WTRW0-ROICE*ACE2F ! new ocean mass
      ENRGW = ENRGW0+ROICE*(ENRGIW-ENRGFI)+(1.-ROICE)*ENRGO ! energy
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
      HSI3 = HSI3-FHSI3
      HSI4 = HSI4+(FHSI3+ENRGFI)
      MSI2 = MSI2+ACE2F ! new ice mass of physical layer 2
      GO TO 270
  240 CONTINUE
C***  NEW ICE IS FORMED BELOW OLD SEA ICE AND ON OPEN OCEAN
      WTRW = WTRW0-(1.-ROICE)*ACEFO-ROICE*ACE2F ! new ocean mass
      ENRGW = ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
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
C***  NO NEW ICE IS FORMED UNDERNEATH THE OLD ONE
      WTRW = WTRW0 ! new ocean mass
      ENRGW = ENRGW0+ROICE*ENRGIW+(1.-ROICE)*ENRGO ! energy of new oc.
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
      GO TO 270
  260 CONTINUE
C***  NEW ICE IS FORMED ON THE OPEN OCEAN
      WTRW = WTRW0-(1.-ROICE)*ACEFO ! new ocean mass
      ENRGW = ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*ENRGIW
      TGW = ENRGW/(WTRW*SHW)+TTRUNC ! new ocean temperature
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
      IF (PLAKE .GT. 0.) GO TO 360 ! no compression for lake ice
C***  COMPRESS THE ICE HORIZONTALLY
      IF (MSI2 .GE. AC2OIM) GO TO 280 ! ice is thick enough
C***  SEA ICE IS TOO THIN
      ROICEN = ROICE*(ACE1I+MSI2)/(ACE1I+AC2OIM) ! new ice concentr.
      GO TO 290
  280 CONTINUE
      OPNOCN=.06*(RHOI/(ACE1I+MSI2)-BYZICX)
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
      DIFS = DRI*ACE1I/ROICE
         EDIFSI = EDIFSI+ROICE*(HSI1+HSI2)*(DIFS/MSI1)*0.5
         DIFSI = DIFSI+ROICE*DIFS
      ROICE = ROICEN
C**** RESAVE PROGNOSTIC QUANTITIES
  360 CONTINUE
  370 CONTINUE
      TG1 = (HSI1/(XSI1*MSI1) +LHM)/SHI ! temperature of layer 1
      TG2 = (HSI2/(XSI2*MSI1) +LHM)/SHI ! temperature of layer 2
      TG3 = (HSI3/(XSI3*MSI2) +LHM)/SHI ! temperature of layer 3
      TG4 = (HSI4/(XSI4*MSI2) +LHM)/SHI ! temperature of layer 4

      RETURN 
      END SUBROUTINE SEA_ICE 

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

      SUBROUTINE init_OCEAN
!@sum init_OCEAN initiallises ocean variables 
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : im,jm,fland,flice,gdata,kocean
      USE OCEAN, only : ota,otb,otc,z12o,dm,iu_osst,iu_sice,iu_ocnml
      USE FILEMANAGER
      IMPLICIT NONE
!@var iu_OHT,iu_MLMAX unit numbers for reading in input files
      INTEGER :: iu_OHT,iu_MLMAX
      INTEGER :: I,J

C**** set up unit numbers for ocean climatologies
      call getunit("OSST",iu_OSST,.TRUE.)
      call getunit("SICE",iu_SICE,.TRUE.)

C**** Read in constant factor relating RSI to MSI from sea ice clim.
      CALL READT (iu_SICE,0,DM,IM*JM,DM,1)

      IF (KOCEAN.eq.1) THEN
C**** Set up unit number of mixed layer depth climatogies
      call getunit("OCNML",iu_OCNML,.TRUE.)

C**** DATA FOR QFLUX MIXED LAYER OCEAN RUNS
C**** read in ocean heat transport coefficients
      call getunit("OHT",iu_OHT,.TRUE.)
      READ (iu_OHT) OTA,OTB,OTC
      CLOSE (iu_OHT)

C**** read in ocean max mix layer depth
      call getunit("MLMAX",iu_MLMAX,.TRUE.)
      CALL READT (iu_MLMAX,0,Z12O,IM*JM,Z12O,1)
      CLOSE (iu_MLMAX)

C**** IF GDATA(I,J,1)<0, THE OCEAN PART WAS CHANGED TO LAND ICE
C**** BECAUSE THE OCEAN ICE REACHED THE MAX MIXED LAYER DEPTH
      DO J=1,JM
        DO I=1,IM
          IF(GDATA(I,J,1).GE.-1.) CYCLE
          FLICE(I,J)=1-FLAND(I,J)+FLICE(I,J)
          FLAND(I,J)=1.
          WRITE(6,'(2I3,'' OCEAN WAS CHANGED TO LAND ICE'')') I,J
        END DO
      END DO
      END IF

      END SUBROUTINE init_OCEAN

      SUBROUTINE daily_OCEAN(IEND)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : IM,JM,KOCEAN
      USE OCEAN, only : ODATA,OSTRUC,OCLIM
      USE DAGCOM, only : AIJ, IJ_ODATA4,IJ_TGO2
      IMPLICIT NONE
      INTEGER I,J,IEND

C**** update ocean related climatologies      
      CALL OCLIM(IEND)

C**** Only do this at end of the day
      IF (KOCEAN.EQ.1.and.IEND.eq.1) THEN
         DO J=1,JM
            DO I=1,IM
               AIJ(I,J,IJ_ODATA4)=AIJ(I,J,IJ_ODATA4)+ODATA(I,J,4)
               AIJ(I,J,IJ_TGO2)=AIJ(I,J,IJ_TGO2)+ODATA(I,J,5)
            END DO
         END DO
C**** RESTRUCTURE THE OCEAN LAYERS AND ELIMINATE SMALL ICE BERGS
         CALL OSTRUC
      END IF

      RETURN
      END

      SUBROUTINE uset_OCEAN
!@sum  uset_OCEAN saves quantities for OHT calculations
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : RHOI
      USE E001M12_COM, only : IM,JM,GDATA
      USE OCEAN, only : OA,Z1I,XSI1,XSI2,XSI3,XSI4
      IMPLICIT NONE
      INTEGER I,J
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  TG1OI  (INSTANTANEOUS AT NOON GMT)
C****       3  TG2OI  (INSTANTANEOUS AT NOON GMT)
C****       4  ENRGP  (INTEGRATED OVER THE DAY)
C****       5  SRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       6  TRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       7  SHDT   (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       8  EVHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       9  TRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      10  SHDT   (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      11  EVHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      12  SRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****
      DO J=1,JM
         DO I=1,IM
            OA(I,J,1)=Z1I*RHOI+GDATA(I,J,1)
            OA(I,J,2)=GDATA(I,J,3)*XSI1+GDATA(I,J,7)*XSI2
            OA(I,J,3)=GDATA(I,J,15)*XSI3+GDATA(I,J,16)*XSI4
         END DO
      END DO

      RETURN
      END

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
