      MODULE OCEAN
!@sum  OCEAN contains all the ocean subroutines
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
!@cont OSTRUC,OCLIM,init_OCEAN,daily_OCEAN
      USE CONSTANT, only : lhm,rhow,rhoi,shw,shi,by12,byshi
      USE E001M12_COM, only : im,jm,lm,focean,flake,fland,fearth
     *     ,flice,kocean,Itime,jmon,jdate,jday,JDendOfM,JDmidOfM,ftype
     *     ,itocean,itlake,itoice,itlkice,itlandi,itearth
      USE PBLCOM
     &     , only : npbl=>n,uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,
     *     cq=>cqgs,ipbl
      USE GEOM
      USE SEAICE_COM, only : rsi,msi,hsi,snowi
      USE SEAICE, only : xsi,ace1i,z1i,ac2oim
      USE LANDICE_COM, only : snowli,tlandi
      USE LAKES_COM, only : tlake,tfl
      USE FLUXES, only : gtemp

      IMPLICIT NONE
!@param LMOM number of layers for deep ocean diffusion
      INTEGER, PARAMETER :: LMOM = 9
!@var OA generic array for ocean heat transport calculations
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
      REAL*8, SAVE,DIMENSION(IM,JM,12) :: OA

!@var TFO temperature of freezing ocean (C)
      REAL*8, PARAMETER :: TFO = -1.8d0
!@var FLEADOC lead fraction for ocean ice (%)
      REAL*8, PARAMETER :: FLEADOC = 0.06d0
!@var TOCEAN temperature of the ocean (C)
      REAL*8, DIMENSION(3,IM,JM) :: TOCEAN

!@var OTA,OTB,OTC ocean heat transport coefficients
      REAL*8, SAVE,DIMENSION(IM,JM,4) :: OTA,OTB
      REAL*8, SAVE,DIMENSION(IM,JM) :: OTC
!@var SINANG,COSANG Fourier coefficients for heat transports
      REAL*8, SAVE :: SINANG,SN2ANG,SN3ANG,SN4ANG,
     *                COSANG,CS2ANG,CS3ANG,CS4ANG

!@var Z1O ocean mixed layer depth
      REAL*8, SAVE,DIMENSION(IM,JM) :: Z1O
!@var Z1OOLD previous ocean mixed layer depth
      REAL*8, SAVE,DIMENSION(IM,JM) :: Z1OOLD
!@var Z12O annual maximum ocean mixed layer depth
      REAL*8, SAVE,DIMENSION(IM,JM) :: Z12O

      REAL*8, DIMENSION(IM,JM) :: DM,AOST,EOST1,EOST0,BOST
     *     ,COST,ARSI,ERSI1,ERSI0,BRSI,CRSI
      INTEGER, DIMENSION(IM,JM) ::  KRSI
      COMMON/OOBS/DM,AOST,EOST1,EOST0,BOST,COST,ARSI,ERSI1,ERSI0,BRSI
     *     ,CRSI,KRSI

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
      REAL*8 TGAVE,DWTRO,WTR1O,WTR2O
C****
C**** FLAND     LAND COVERAGE (1)
C****
C**** TOCEAN 1  OCEAN TEMPERATURE OF FIRST LAYER (C)
C****        2  MEAN OCEAN TEMPERATURE OF SECOND LAYER (C)
C****        3  OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)
C****     RSI   RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****     MSI   OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
C**** SNOWI     OCEAN ICE SNOW AMOUNT (KG/M**2)
C**** HSI(1:2)  OCEAN ICE ENTHALPY OF FIRST/SECOND LAYER (C)
C****
C**** RESTRUCTURE OCEAN LAYERS
C****
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF (FLAND(I,J).GE.1.) CYCLE
          IF (Z1OOLD(I,J).GE.Z12O(I,J)) GO TO 140
          IF (Z1O(I,J).EQ.Z1OOLD(I,J)) CYCLE
          WTR1O=RHOW*Z1O(I,J)-RSI(I,J)*(SNOWI(I,J)
     *         +ACE1I+MSI(I,J))
          DWTRO=RHOW*(Z1O(I,J)-Z1OOLD(I,J))
          WTR2O=RHOW*(Z12O(I,J)-Z1O(I,J))
          IF (DWTRO.GT.0.) GO TO 120
C**** MIX LAYER DEPTH IS GETTING SHALLOWER
          TOCEAN(2,I,J)=TOCEAN(2,I,J)
     *         +((TOCEAN(2,I,J)-TOCEAN(1,I,J))*DWTRO/WTR2O)
          CYCLE
C**** MIX LAYER DEPTH IS GETTING DEEPER
 120      TGAVE=(TOCEAN(2,I,J)*DWTRO+(2.*TOCEAN(2,I,J)-TOCEAN(3,I,J))
     *         *WTR2O)/(WTR2O+DWTRO)
          TOCEAN(1,I,J)=TOCEAN(1,I,J)+((TGAVE-TOCEAN(1,I,J))
     *         *DWTRO/WTR1O)
          IF (Z1O(I,J).GE.Z12O(I,J)) GO TO 140
          TOCEAN(2,I,J)=TOCEAN(2,I,J)
     *         +((TOCEAN(3,I,J)-TOCEAN(2,I,J))*DWTRO/(WTR2O+DWTRO))
          CYCLE
C**** MIXED LAYER DEPTH IS AT ITS MAXIMUM OR TEMP PROFILE IS UNIFORM
 140      TOCEAN(2,I,J)=TOCEAN(1,I,J)
          TOCEAN(3,I,J)=TOCEAN(1,I,J)
        END DO
      END DO
      RETURN
      END SUBROUTINE OSTRUC

      SUBROUTINE OCLIM(IDOZ1O)
!@sum OCLIM calculates daily ocean data from ocean/sea ice climatologies
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean or fixed SST/fixed lakes)
      IMPLICIT NONE

      REAL*8 XZO(IM,JM),XZN(IM,JM)
      INTEGER, INTENT(IN) :: IDOZ1O

      INTEGER n,MD,J,I,LSTMON,K,MDMAX,IMAX,IMON
      REAL*8 PLICEN,PLICE,POICE,POCEAN,RSICSQ,ZIMIN,ZIMAX,X1
     *     ,X2,Z1OMIN,RSINEW,TIME,FRAC,MSINEW
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
C**** CALCULATE DAILY OCEAN DATA FROM CLIMATOLOGY
C****
C**** TOCEAN(1)  OCEAN TEMPERATURE (C)
C****  RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****  MSI  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
C**** READ IN OBSERVED OCEAN DATA
      IF (JMON.EQ.MONTHO) GO TO 400
      IF (MONTHO.EQ.0) THEN
C****    READ IN LAST MONTH'S END-OF-MONTH DATA
         LSTMON=JMON-1
         IF(LSTMON.EQ.0) LSTMON=12
         CALL READT (iu_OSST,IM*JM,EOST0,IM*JM,EOST0,LSTMON)
         CALL READT (iu_SICE,IM*JM,ERSI0,IM*JM,ERSI0,LSTMON)
      ELSE
C****    COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
        EOST0=EOST1
        ERSI0=ERSI1
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      MONTHO=JMON
      IF (JMON.EQ.1) THEN
         REWIND iu_OSST
         REWIND iu_SICE
         READ (iu_SICE)         ! skip over DM-record
      END IF
      CALL READT (iu_OSST,0,AOST,2*IM*JM,AOST,1) ! READS AOST,EOST1
      CALL READT (iu_SICE,0,ARSI,2*IM*JM,ARSI,1) ! READS ARSI,ERSI1
C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          BOST(I,J)=EOST1(I,J)-EOST0(I,J)
          COST(I,J)=3.*(EOST1(I,J)+EOST0(I,J)) - 6.*AOST(I,J)
          BRSI(I,J)=0.
          CRSI(I,J)=0.          ! extreme cases
          KRSI(I,J)=0           ! ice constant
          IF(ARSI(I,J).LE.0. .or. ARSI(I,J).GT.1.) CYCLE
          BRSI(I,J)=ERSI1(I,J)-ERSI0(I,J) ! quadratic fit
          CRSI(I,J)=3.*(ERSI1(I,J)+ERSI0(I,J)) - 6.*ARSI(I,J)
          MD=MDMAX+JDATE-16
          IF(ABS(CRSI(I,J)) .GT. ABS(BRSI(I,J))) THEN ! linear fits
            RSICSQ=CRSI(I,J)*(ARSI(I,J)*CRSI(I,J) - .25*BRSI(I,J)**2 -
     *           CRSI(I,J)**2*BY12)
            IF(RSICSQ.lt.0.) THEN
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
              KRSI(I,J) = -1
              BRSI(I,J) = .5*(ERSI0(I,J)**2 + ERSI1(I,J)**2) / ARSI(I,J)
            ELSEIF(RSICSQ.gt.CRSI(I,J)**2)  then
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
              KRSI(I,J) = 1
              BRSI(I,J) = .5*((ERSI0(I,J)-1.)**2 + (ERSI1(I,J)-1.)**2) /
     /             (ARSI(I,J)-1.)
            END IF
          END IF
        END DO
      END DO
C**** Calculate OST, RSI and MSI for current day
  400 TIME=(JDATE-.5)/(JDendOFM(JMON)-JDendOFM(JMON-1))-.5 ! -.5<TIME<.5
      DO J=1,JM
        ZIMIN=.5d0
        ZIMAX=2d0
        IF(J.GT.JM/2) ZIMAX=3.5d0
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF(FLAND(I,J).GE.1.) CYCLE
C**** OST always uses quadratic fit
          IF (FOCEAN(I,J).gt.0) THEN
            TOCEAN(1,I,J)=AOST(I,J)+BOST(I,J)*TIME
     *                   +COST(I,J)*(TIME**2-BY12)
          ELSE
            TLAKE(I,J)=AOST(I,J)+BOST(I,J)*TIME
     *                +COST(I,J)*(TIME**2-BY12)
          END IF
          SELECT CASE (KRSI(I,J))
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
          CASE (-1)
            IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .gt. 0.)  then
              RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5) !  TIME < T0
            ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .gt. 0.)  then
              RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME) !  T1 < TIME
            ELSE
              RSINEW = 0.       !  T0 < TIME < T1
            END IF
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
          CASE (1)
            IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .lt. 1.)  then
              RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5) !  TIME < T0
            ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .lt. 1.)  then
              RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME) !  T1 < TIME
            ELSE
              RSINEW = 1.       !  T0 < TIME < T1
            END IF
C**** RSI uses quadratic fit
          CASE (0)
            RSINEW=ARSI(I,J)+BRSI(I,J)*TIME+CRSI(I,J)*(TIME**2-BY12)
          END SELECT
          RSI(I,J)=RSINEW
          MSINEW=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** adjust enthalpy so that temperature remains constant
          HSI(3:4,I,J)=HSI(3:4,I,J)*MSINEW/MSI(I,J)
          MSI(I,J)=MSINEW
C**** set ftype arrays
          IF (FOCEAN(I,J).gt.0) THEN
            FTYPE(ITOICE ,I,J)=FOCEAN(I,J)*    RSI(I,J)
            FTYPE(ITOCEAN,I,J)=FOCEAN(I,J)*(1.-RSI(I,J))
          ELSE
            FTYPE(ITLKICE,I,J)= FLAKE(I,J)*    RSI(I,J)
            FTYPE(ITLAKE ,I,J)= FLAKE(I,J)*(1.-RSI(I,J))
          END IF
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONABLE VALUE
          IF (TOCEAN(1,I,J).LT.TFO) TOCEAN(1,I,J)=TFO
          IF (TLAKE(I,J).LT.TFL) TLAKE(I,J)=TFL
C**** REDUCE THE RATIO OF OCEAN ICE TO WATER BY .1*RHOI/ACEOI
c     IF (RSI(I,J).GT.0.) THEN
c        BYZICE=RHOI/(Z1I*RHOI+MSI(I,J))
c        RSI(I,J)=RSI(I,J)*(1.-.06d0*(BYZICE-0.2d0))
c     END IF
C**** ZERO OUT SNOWOI, TG1OI, TG2OI IF THERE IS NO OCEAN ICE
          IF (RSI(I,J).LE.0.) THEN
            HSI(1:2,I,J)=-LHM*XSI(1:2)*ACE1I
            HSI(3:4,I,J)=-LHM*XSI(3:4)*AC2OIM
            SNOWI(I,J)=0.
            GTEMP(1:2,2,I,J)=0.
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE
      DO I=2,IM
        SNOWI(I,JM)=SNOWI(1,JM)
        TOCEAN(1,I,JM)=TOCEAN(1,1,JM)
        TLAKE(I,JM)=TLAKE(1,JM)
        RSI(I,JM)=RSI(1,JM)
        MSI(I,JM)=MSI(1,JM)
        HSI(:,I,JM)=HSI(:,1,JM)
        GTEMP(1:2,2,I,JM)=GTEMP(1:2,2,1,JM)
C**** set ftype arrays
        IF (FOCEAN(1,JM).gt.0) THEN
          FTYPE(ITOICE ,I,JM)=FOCEAN(1,JM)*    RSI(1,JM)
          FTYPE(ITOCEAN,I,JM)=FOCEAN(1,JM)*(1.-RSI(1,JM))
        ELSE
          FTYPE(ITLKICE,I,JM)= FLAKE(1,JM)*    RSI(1,JM)
          FTYPE(ITLAKE ,I,JM)= FLAKE(1,JM)*(1.-RSI(1,JM))
        END IF
      END DO
      RETURN
C****
C**** CALCULATE DAILY OCEAN MIXED LAYER DEPTHS FROM CLIMATOLOGY
C****
C**** SAVE PREVIOUS DAY'S MIXED LAYER DEPTH IN WORK2
 500  Z1OOLD=Z1O
C**** COMPUTE Z1O ONLY AT THE BEGINNING OF A DAY (OR AT ItimeI)
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
          IF (JDAY.GT.JDmidOFM(IMON) .and. IMON.LE.12) GO TO 520
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
        IF (JDAY.LE.JDmidOFM(IMON)) GO TO 530
        IMON=IMON+1          ! read in new month of climatological data
        XZO = XZN
        IF (IMON.EQ.13) REWIND iu_OCNML
      END IF
      CALL READT (iu_OCNML,0,XZN,IM*JM,XZN,1)
 530  JDLAST=JDAY
      REWIND iu_OCNML
C**** INTERPOLATE OCEAN DATA TO CURRENT DAY
      FRAC = DBLE(JDmidOFM(IMON)-JDAY)/(JDmidOFM(IMON)-JDmidOFM(IMON-1))
      DO J=1,JM
      DO I=1,IM
      Z1O(I,J)=FRAC*XZO(I,J)+(1.-FRAC)*XZN(I,J)

      IF (RSI(I,J)*FOCEAN(I,J).le.0) CYCLE
C**** MIXED LAYER DEPTH IS INCREASED TO OCEAN ICE DEPTH + 1 METER
      Z1OMIN=1.+(RHOI*Z1I+SNOWI(I,J)+MSI(I,J))/RHOW
      IF (Z1O(I,J).GE.Z1OMIN) GO TO 605
      WRITE(6,602) ITime,I,J,JMON,Z1O(I,J),Z1OMIN
  602 FORMAT (' INCREASE OF MIXED LAYER DEPTH ',I10,3I4,2F10.3)
      Z1O(I,J)=Z1OMIN
  605 IF (Z1OMIN.LE.Z12O(I,J)) CYCLE
C**** ICE DEPTH+1>MAX MIXED LAYER DEPTH : CHANGE OCEAN TO LAND ICE
      PLICE=FLICE(I,J)
      PLICEN=1.-FEARTH(I,J)
      POICE=FOCEAN(I,J)*RSI(I,J)
      POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
      SNOWLI(I,J)=(SNOWLI(I,J)*PLICE+SNOWI(I,J)*POICE)/PLICEN
      TLANDI(1,I,J)=(TLANDI(1,I,J)*PLICE+
     *     ((HSI(1,I,J)+HSI(2,I,J))/(SNOWI(I,J)+ACE1I)+LHM)*BYSHI*POICE+
     *     (LHM+SHW*TOCEAN(1,I,J))*POCEAN/SHI)/PLICEN
      TLANDI(2,I,J)=(TLANDI(2,I,J)*PLICE+
     *     ((HSI(3,I,J)+HSI(4,I,J))/MSI(I,J)+LHM)*BYSHI*POICE+
     *     (LHM+SHW*TOCEAN(1,I,J))*POCEAN/SHI)/PLICEN
      FLAND(I,J)=1.
      FLICE(I,J)=PLICEN
      FOCEAN(I,J)=0.
C**** set ftype/gtemp array. Summation is necessary for cases where
C**** Earth and Land Ice are lumped together
      FTYPE(ITOICE ,I,J)=0.
      FTYPE(ITOCEAN,I,J)=0.
      FTYPE(ITLANDI,I,J)=0.
      FTYPE(ITEARTH,I,J)=FEARTH(I,J)
      FTYPE(ITLANDI,I,J)=FTYPE(ITLANDI,I,J)+FLICE(I,J)
      GTEMP(1:2,3,I,J)  =TLANDI(1:2,I,J)
C**** MARK THE POINT FOR RESTART PURPOSES
      SNOWI(I,J)=-10000.-SNOWI(I,J)
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
      WRITE(6,606) 100.*POICE,100.*POCEAN,ITime,I,J
  606 FORMAT(F6.1,'% OCEAN ICE AND',F6.1,'% OPEN OCEAN WERE',
     *  ' CHANGED TO LAND ICE AT (I)Time,I,J',I10,2I4)
      END DO
      END DO
C**** PREVENT Z1O, THE MIXED LAYER DEPTH, FROM EXCEEDING Z12O
      DO J=1,JM
        DO I=1,IM
          IF (Z1O(I,J).GT.Z12O(I,J)-.01) Z1O(I,J)=Z12O(I,J)
        END DO
      END DO
C**** SET MARKER INDICATING THAT Z1O HAS BEEN SET
      Z1O(IM,1)=-9999.
      RETURN
      END SUBROUTINE OCLIM

      SUBROUTINE OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,F0DT,F2DT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFW,RUN4O,ERUN4O,RUN4I,ERUN4I
     *           ,ENRGFO,ACEFO,ACE2F,ENRGFI)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@var TFW freezing temperature for water underlying ice (C)
      REAL*8, INTENT(IN) :: TFW

!@var TGW mixed layer temp.(C)
      REAL*8, INTENT(INOUT) :: TGW
      REAL*8, INTENT(IN) :: ROICE, SMSI, WTRO, EVAPO, EVAPI, RUN0
      REAL*8, INTENT(IN) :: F0DT, F2DT, OTDT, RVRRUN, RVRERUN
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACE2F, RUN4O, RUN4I
     *     , ERUN4O, ERUN4I
      REAL*8 EIW0, ENRGIW, WTRI1, EFIW, ENRGO0, EOFRZ, ENRGO
      REAL*8 WTRW0, ENRGW0, WTRW, ENRGW, WTRI0, SMSI0
C**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACE2F=0. ; ENRGFI=0.

C**** Calculate extra mass flux to ocean, balanced by deep removal
      RUN4O=-EVAPO+RVRRUN  ! open ocean
      RUN4I=-EVAPI+RVRRUN  ! under ice
      ERUN4O=RUN4O*TGW*SHW ! corresponding heat flux at bottom (open)
      ERUN4I=RUN4I*TGW*SHW !                              (under ice)

C**** Calculate heat fluxes to ocean
      ENRGO  = F0DT+OTDT+RVRERUN-ERUN4O ! in open water
      ENRGIW = F2DT+OTDT+RVRERUN-ERUN4I ! under ice

C**** Calculate energy in mixed layer (open ocean)
      ENRGO0=WTRO*TGW*SHW
      EOFRZ=WTRO*TFW*SHW
      IF (ENRGO0+ENRGO.GE.EOFRZ) THEN
C**** FLUXES RECOMPUTE TGW WHICH IS ABOVE FREEZING POINT FOR OCEAN
        IF (ROICE.LE.0.) THEN
          TGW=TGW+(ENRGO/(WTRO*SHW))
          RETURN
        END IF
      ELSE
C**** FLUXES COOL TGO TO FREEZING POINT FOR OCEAN AND FORM SOME ICE
        ACEFO=(ENRGO0+ENRGO-EOFRZ)/(TFW*(SHI-SHW)-LHM)
        ENRGFO=ACEFO*(TFW*SHI-LHM)
        IF (ROICE.LE.0.) THEN
          TGW=TFW
          RETURN
        END IF
      END IF
C****
      IF (ROICE.le.0) RETURN

C**** Calculate initial energy of mixed layer (before fluxes applied)
      SMSI0=SMSI+RUN0+EVAPI  ! previous ice mass
      WTRI0=WTRO-SMSI0       ! mixed layer mass below ice (kg/m^2)
      EIW0=WTRI0*TGW*SHW     ! energy of mixed layer below ice (J/m^2)

      WTRW0=WTRO-ROICE*SMSI0 ! initial total mixed layer mass
      ENRGW0=WTRW0*TGW*SHW   ! initial energy in mixed layer

C**** CALCULATE THE ENERGY OF THE WATER BELOW THE ICE AT FREEZING POINT
C**** AND CHECK WHETHER NEW ICE MUST BE FORMED
      WTRI1 = WTRO-SMSI         ! new mass of ocean (kg/m^2)
      EFIW = WTRI1*TFW*SHW ! freezing energy of ocean mass WTRI1
      IF (EIW0+ENRGIW .LE. EFIW) THEN ! freezing case
C**** FLUXES WOULD COOL TGW TO FREEZING POINT AND FREEZE SOME MORE ICE
        ACE2F = (EIW0+ENRGIW-EFIW)/(TFW*(SHI-SHW)-LHM)
        ENRGFI = ACE2F*(TFW*SHI-LHM) ! energy of frozen ice
      END IF
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      WTRW = WTRO  -(1.-ROICE)* ACEFO        -ROICE*(SMSI+ACE2F)
      ENRGW= ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)
      TGW  = ENRGW/(WTRW*SHW) ! new ocean temperature

      RETURN
      END SUBROUTINE OSOURC

      END MODULE OCEAN

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean are reasonable
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : im,jm
      USE OCEAN, only : tocean
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in ocean data
      CALL CHECK3(TOCEAN,3,IM,JM,SUBR,'toc')

      END SUBROUTINE CHECKO

      SUBROUTINE init_OCEAN
!@sum init_OCEAN initiallises ocean variables
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : im,jm,fland,flice,kocean,ftype,focean
     *     ,itocean,itoice,itearth,itlandi,fearth
      USE OCEAN, only : ota,otb,otc,z12o,dm,iu_osst,iu_sice,iu_ocnml
     *     ,tocean
      USE SEAICE_COM, only : snowi,rsi
      USE FLUXES, only : gtemp
      USE FILEMANAGER
      IMPLICIT NONE
!@var iu_OHT,iu_MLMAX unit numbers for reading in input files
      INTEGER :: iu_OHT,iu_MLMAX
      INTEGER :: I,J

C**** set up unit numbers for ocean climatologies
      call getunit("OSST",iu_OSST,.true.,.true.)
      call getunit("SICE",iu_SICE,.true.,.true.)

C**** Read in constant factor relating RSI to MSI from sea ice clim.
      CALL READT (iu_SICE,0,DM,IM*JM,DM,1)

      IF (KOCEAN.eq.1) THEN
C**** Set up unit number of mixed layer depth climatogies
      call getunit("OCNML",iu_OCNML,.true.,.true.)

C**** DATA FOR QFLUX MIXED LAYER OCEAN RUNS
C**** read in ocean heat transport coefficients
      call getunit("OHT",iu_OHT,.true.,.true.)
      READ (iu_OHT) OTA,OTB,OTC
      CLOSE (iu_OHT)

C**** read in ocean max mix layer depth
      call getunit("MLMAX",iu_MLMAX,.true.,.true.)
      CALL READT (iu_MLMAX,0,Z12O,IM*JM,Z12O,1)
      CLOSE (iu_MLMAX)

C**** IF SNOWI(I,J)<0, THE OCEAN PART WAS CHANGED TO LAND ICE
C**** BECAUSE THE OCEAN ICE REACHED THE MAX MIXED LAYER DEPTH
      DO J=1,JM
        DO I=1,IM
          IF(SNOWI(I,J).GE.-1.) CYCLE
          FLICE(I,J)=1-FLAND(I,J)+FLICE(I,J)
          FLAND(I,J)=1.
          FEARTH(I,J)=1.-FLICE(I,J)
          FOCEAN(I,J)=0.
          WRITE(6,'(2I3,'' OCEAN WAS CHANGED TO LAND ICE'')') I,J
C**** Reset ftype array. Summation is necessary for cases where Earth
C**** and Land Ice are lumped together
          FTYPE(ITLANDI,I,J)=0.
          FTYPE(ITEARTH,I,J)=FEARTH(I,J)
          FTYPE(ITLANDI,I,J)=FTYPE(ITLANDI,I,J)+FLICE(I,J)
        END DO
      END DO
      END IF
C**** Set ftype array for oceans
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          FTYPE(ITOICE ,I,J)=FOCEAN(I,J)*RSI(I,J)
          FTYPE(ITOCEAN,I,J)=FOCEAN(I,J)*(1.-RSI(I,J))
          GTEMP(1:2,1,I,J)=TOCEAN(1:2,I,J)
        END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE init_OCEAN

      SUBROUTINE daily_OCEAN(IEND)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : rhow,shw,twopi,edpery
      USE E001M12_COM, only : im,jm,kocean,focean,jday,ftype,itocean
     *     ,itoice,fland,flake
      USE OCEAN, only : tocean,ostruc,oclim,z1O,
     *     sinang,sn2ang,sn3ang,sn4ang,cosang,cs2ang,cs3ang,cs4ang
      USE DAGCOM, only : aij,ij_toc2,ij_tgo2
      USE SEAICE_COM, only : rsi,msi,hsi,snowi
      USE SEAICE, only : simelt,ace1i,lmi
      USE LAKES_COM, only : tlake
      USE GEOM, only : imaxj
      USE FLUXES, only : gtemp
      IMPLICIT NONE
      INTEGER I,J,IEND,IMAX
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL
      REAL*8 MSI2,ROICE,SNOW,TGW,WTRO,WTRW,ENRGW,ENRGUSED,ANGLE

C**** update ocean related climatologies
      CALL OCLIM(IEND)

C**** Set fourier coefficients for heat transport calculations
      IF (KOCEAN.eq.1) THEN
        ANGLE=TWOPI*JDAY/EDPERY
        SINANG=SIN(ANGLE)
        SN2ANG=SIN(2*ANGLE)
        SN3ANG=SIN(3*ANGLE)
        SN4ANG=SIN(4*ANGLE)
        COSANG=COS(ANGLE)
        CS2ANG=COS(2*ANGLE)
        CS3ANG=COS(3*ANGLE)
        CS4ANG=COS(4*ANGLE)
      END IF

C**** Only do this at end of the day
      IF (KOCEAN.EQ.1.and.IEND.eq.1) THEN
        DO J=1,JM
          DO I=1,IM
            AIJ(I,J,IJ_TOC2)=AIJ(I,J,IJ_TOC2)+TOCEAN(2,I,J)
            AIJ(I,J,IJ_TGO2)=AIJ(I,J,IJ_TGO2)+TOCEAN(3,I,J)
          END DO
        END DO
C**** DO DEEP DIFFUSION
c        IF (QDEEP) CALL ODIFS
C**** RESTRUCTURE THE OCEAN LAYERS
        CALL OSTRUC
C**** AND ELIMINATE SMALL AMOUNTS OF SEA ICE
        DO J=1,JM
          IMAX=IMAXJ(J)
          DO I=1,IMAX
C**** Only melting of ocean ice (not lakes)
            IF (FOCEAN(I,J)*RSI(I,J) .GT. 0.) THEN
C**** REDUCE ICE EXTENT IF OCEAN TEMPERATURE IS GREATER THAN ZERO
C**** (MELTING POINT OF ICE)
            TGW=TOCEAN(1,I,J)
            IF (TGW.LE.0.) CYCLE
            ROICE=RSI(I,J)
            MSI2=MSI(I,J)
            SNOW=SNOWI(I,J)   ! snow mass
            HSIL(:)= HSI(:,I,J) ! sea ice enthalpy
            WTRO=Z1O(I,J)*RHOW
            WTRW=WTRO-ROICE*(SNOW + ACE1I + MSI2)
            ENRGW=WTRW*TGW*SHW  ! energy of water available for melting
            CALL SIMELT(ROICE,SNOW,MSI2,HSIL,TSIL,ENRGW,ENRGUSED)
C**** RESAVE PROGNOSTIC QUANTITIES
            TGW=(ENRGW-ENRGUSED)/(WTRO*SHW)
            TOCEAN(1,I,J)=TGW
            RSI(I,J)=ROICE
            MSI(I,J)=MSI2
            SNOWI(I,J)=SNOW
            HSI(:,I,J)=HSIL(:)
C**** set ftype/gtemp arrays
            FTYPE(ITOCEAN,I,J)=FOCEAN(I,J)*(1.-RSI(I,J))
            FTYPE(ITOICE ,I,J)=FOCEAN(I,J)*    RSI(I,J)
            GTEMP(1:2,2,I,J) = TSIL(1:2)
            END IF
          END DO
        END DO
      END IF

C**** set gtemp array for ocean temperature
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          GTEMP(1:2,1,I,J) = TOCEAN(1:2,I,J)
        ELSEIF (FLAKE(I,J).gt.0) THEN
          GTEMP(1  ,1,I,J) = TLAKE(I,J)
        END IF
      END DO
      END DO
C****
      RETURN
      END SUBROUTINE daily_OCEAN

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean reads and writes ocean arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
      USE OCEAN
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "OCN01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,TOCEAN,OA,Z1O
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,TOCEAN,OA,Z1O
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocean

      SUBROUTINE PRECIP_OC
!@sum  PRECIP_OC driver for applying precipitation to ocean fraction
!@auth Original Development Team
!@ver  1.0
!@calls
      USE CONSTANT, only : rhow,shw
      USE E001M12_COM, only : im,jm,focean,kocean,itocean,itoice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi,prec,eprec,gtemp
      USE OCEAN, only : oa,tocean,z1o
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : ace1i
      USE DAGCOM, only : aj,aij,ij_f0oc,j_run2,j_dwtr2
      IMPLICIT NONE
      REAL*8 TGW,PRCP,WTRO,ENRGP,ERUN4,ENRGO,POCEAN,POICE,SNOW
     *     ,SMSI,ENRGW,WTRW0,WTRW,RUN0,RUN4,DXYPJ,ROICE
      INTEGER I,J,IMAX

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
        ROICE=RSI(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        IF (FOCEAN(I,J).gt.0) THEN

          PRCP=PREC(I,J)
          ENRGP=EPREC(I,J)
          OA(I,J,4)=OA(I,J,4)+ENRGP
          AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+ENRGP*POCEAN

          IF (KOCEAN .EQ. 1) THEN
            TGW=TOCEAN(1,I,J)
            WTRO=Z1O(I,J)*RHOW
            RUN0=RUNOSI(I,J)
            SNOW=SNOWI(I,J)
            SMSI=MSI(I,J)+ACE1I+SNOW
            RUN4=PRCP
            ERUN4=RUN4*TGW*SHW

            IF (POICE.LE.0.) THEN
              ENRGW=TGW*WTRO*SHW + ENRGP - ERUN4
              WTRW =WTRO
            ELSE
              WTRW0=WTRO -ROICE*SMSI
              ENRGW=WTRW0*TGW*SHW + (1.-ROICE)*ENRGP - ERUN4
              WTRW =WTRW0+ROICE*(RUN0-RUN4)
            END IF
            TGW=ENRGW/(WTRW*SHW)
            TOCEAN(1,I,J)=TGW
            AJ(J,J_RUN2 ,ITOCEAN)=AJ(J,J_RUN2 ,ITOCEAN)+RUN4 *POCEAN
            AJ(J,J_DWTR2,ITOCEAN)=AJ(J,J_DWTR2,ITOCEAN)+ERUN4*POCEAN
            AJ(J,J_RUN2 ,ITOICE) =AJ(J,J_RUN2 ,ITOICE) +RUN4 *POICE
            AJ(J,J_DWTR2,ITOICE) =AJ(J,J_DWTR2,ITOICE) +ERUN4*POICE
          END IF
          GTEMP(1,1,I,J)=TOCEAN(1,I,J)
        END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE PRECIP_OC

      SUBROUTINE GROUND_OC
!@sum  GROUND_OC driver for applying surface fluxes to ocean fraction
!@auth Original Development Team
!@ver  1.0
!@calls
      USE CONSTANT, only : rhow,shw
      USE E001M12_COM, only : im,jm,focean,kocean,jday,dtsrc,itocean
     *     ,itoice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi,erunosi,e0,e1,evapor,dmsi,dhsi,
     *     flowo,eflowo,gtemp
      USE OCEAN, only : tocean,z1o,oa,ota,otb,otc,tfo,osourc,
     *     sinang,sn2ang,sn3ang,sn4ang,cosang,cs2ang,cs3ang,cs4ang
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : ace1i
      USE DAGCOM, only : aj,aij,areg,jreg,ij_f0oc,j_run2
     *     ,j_dwtr2,j_tg1,j_tg2,j_evap,j_oht,j_omlt,j_erun2,j_imelt
     *     ,ij_tgo,ij_tg1,ij_evap,ij_evapo,j_type
      IMPLICIT NONE
C**** grid box variables
      REAL*8 POCEAN, POICE, DXYPJ
C**** prognostic variables
      REAL*8 TGW, WTRO, SMSI, ROICE
C**** fluxes
      REAL*8 EVAPO, EVAPI, F2DT, F0DT, OTDT, RVRRUN, RVRERUN, RUN0
C**** output from OSOURC
      REAL*8 ERUN4I, ERUN4O, RUN4I, RUN4O, ENRGFO, ACEFO, ACE2F, ENRGFI

      INTEGER I,J,IMAX,JR

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
        JR=JREG(I,J)
        ROICE=RSI(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        IF (FOCEAN(I,J).gt.0) THEN

          TGW  =TOCEAN(1,I,J)
          EVAPO=EVAPOR(I,J,1)
          EVAPI=EVAPOR(I,J,2)   ! evapor/dew at the ice surface (kg/m^2)
          F0DT =E0(I,J,1)
          SMSI =MSI(I,J)+ACE1I+SNOWI(I,J)
C**** get ice-ocean fluxes from sea ice routine
          RUN0=RUNOSI(I,J) ! includes ACE2M
          F2DT=ERUNOSI(I,J)
C**** get river runoff
          RVRRUN = FLOWO(I,J)/(FOCEAN(I,J)*DXYPJ)
          RVRERUN=EFLOWO(I,J)/(FOCEAN(I,J)*DXYPJ)
          OA(I,J,4)=OA(I,J,4)+RVRERUN    ! add to surface energy budget

          AJ(J,J_TG1, ITOCEAN)=AJ(J,J_TG1, ITOCEAN)+TGW  *POCEAN
          AJ(J,J_EVAP,ITOCEAN)=AJ(J,J_EVAP,ITOCEAN)+EVAPO*POCEAN
          AJ(J,J_TYPE,ITOCEAN)=AJ(J,J_TYPE,ITOCEAN)+      POCEAN
          AJ(J,J_TYPE,ITOICE) =AJ(J,J_TYPE,ITOICE) +      POICE
          IF (JR.ne.24) AREG(JR,J_TG1)=AREG(JR,J_TG1)+TGW*POCEAN*DXYPJ
          AIJ(I,J,IJ_TGO)  =AIJ(I,J,IJ_TGO)  +TGW
          AIJ(I,J,IJ_TG1)  =AIJ(I,J,IJ_TG1)  +TGW  *POCEAN
          AIJ(I,J,IJ_EVAP) =AIJ(I,J,IJ_EVAP) +EVAPO*POCEAN
          AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+EVAPO*POCEAN

          IF (KOCEAN .EQ. 1) THEN
            WTRO=Z1O(I,J)*RHOW
            OTDT=DTSRC*(OTA(I,J,4)*SN4ANG+OTB(I,J,4)*CS4ANG
     *           +OTA(I,J,3)*SN3ANG+OTB(I,J,3)*CS3ANG
     *           +OTA(I,J,2)*SN2ANG+OTB(I,J,2)*CS2ANG
     *           +OTA(I,J,1)*SINANG+OTB(I,J,1)*COSANG+OTC(I,J))

C**** Calculate the amount of ice formation
            CALL OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,F0DT,F2DT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFO,RUN4O,ERUN4O,RUN4I,ERUN4I
     *           ,ENRGFO,ACEFO,ACE2F,ENRGFI)

C**** Resave prognostic variables
            TOCEAN(1,I,J)=TGW
C**** Open Ocean diagnostics
            AJ(J,J_TG2  ,ITOCEAN)=AJ(J,J_TG2  ,ITOCEAN)+TOCEAN(2,I,J)
     *           *POCEAN
            AJ(J,J_OHT  ,ITOCEAN)=AJ(J,J_OHT  ,ITOCEAN)+OTDT  *POCEAN
            AJ(J,J_RUN2 ,ITOCEAN)=AJ(J,J_RUN2 ,ITOCEAN)+RUN4O *POCEAN
            AJ(J,J_OMLT ,ITOCEAN)=AJ(J,J_OMLT ,ITOCEAN)+TOCEAN(3,I,J)
     *           *POCEAN
            AJ(J,J_DWTR2,ITOCEAN)=AJ(J,J_DWTR2,ITOCEAN)+ERUN4O*POCEAN
            AJ(J,J_ERUN2,ITOCEAN)=AJ(J,J_ERUN2,ITOCEAN)-ENRGFO*POCEAN
            AJ(J,J_IMELT,ITOCEAN)=AJ(J,J_IMELT,ITOCEAN)-ACEFO *POCEAN
            IF (JR.ne.24) AREG(JR,J_TG2)=AREG(JR,J_TG2)+TOCEAN(2,I,J)
     *           *POCEAN*DXYPJ
            AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+F0DT*POCEAN
C**** Ice-covered ocean diagnostics
            AJ(J,J_OHT  ,ITOICE)=AJ(J,J_OHT  ,ITOICE)+OTDT  *POICE
            AJ(J,J_RUN2 ,ITOICE)=AJ(J,J_RUN2 ,ITOICE)+RUN4I *POICE
            AJ(J,J_DWTR2,ITOICE)=AJ(J,J_DWTR2,ITOICE)+ERUN4I*POICE
            AJ(J,J_ERUN2,ITOICE)=AJ(J,J_ERUN2,ITOICE)-ENRGFI*POICE
            AJ(J,J_IMELT,ITOICE)=AJ(J,J_IMELT,ITOICE)-ACE2F *POICE
          ELSE
            ACEFO=0 ; ACE2F=0. ; ENRGFO=0. ; ENRGFI=0.
            AJ(J,J_TG2,ITOCEAN)  =AJ(J,J_TG2,ITOCEAN)  +TGW   *POCEAN
            IF (JR.ne.24) AREG(JR,J_TG2)=AREG(JR,J_TG2)+TGW*POCEAN*DXYPJ
          END IF

C**** Store mass and energy fluxes for formation of sea ice
          DMSI(1,I,J)=ACEFO
          DMSI(2,I,J)=ACE2F
          DHSI(1,I,J)=ENRGFO
          DHSI(2,I,J)=ENRGFI
C**** store surface temperatures
          GTEMP(1:2,1,I,J)=TOCEAN(1:2,I,J)

        END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE GROUND_OC

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates zonal ocean energy
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shw,rhow
      USE E001M12_COM, only : im,jm,fim,focean,kocean
      USE OCEAN, only : tocean,z1o,z12o
c      USE OCEAN_COM, only : dz,rtgo
      USE GEOM, only : imaxj
      IMPLICIT NONE
!@var OCEANE zonal ocean energy (J/M^2)
      REAL*8, DIMENSION(JM) :: OCEANE
      INTEGER I,J

      OCEANE=0
      IF (KOCEAN.ne.0) THEN
      DO J=1,JM
        DO I=I,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
            OCEANE(J)=OCEANE(J)+(TOCEAN(1,I,J)*Z1O(I,J)
     *           +TOCEAN(2,I,J)*(Z12O(I,J)-Z1O(I,J)))*SHW*RHOW
c     IF (QDEEPO) THEN
c     DO L=1,LMOM
c     OCEANE(J)=OCEANE(J)+(RTGO(L,I,J)*DZ(L)*SHW*RHOW)
c     END DO
c     END IF
          END IF
        END DO
      END DO
      OCEANE(1) =FIM*OCEANE(1)
      OCEANE(JM)=FIM*OCEANE(JM)
      END IF
C****
      END SUBROUTINE conserv_OCE

C**** Things to do:
C**** Add QDEEPO option
C**** CALL ODIFS ! before OSTRUC
C**** Add RTGO to acc file
c    *    (((SNGL(RTGO(L,I,J)),L=2,lmom),I=1,IM),J=1,JM)
C**** INITIALIZE DEEP OCEAN ARRAYS
c      STG3=0. ; DTG3=0 ; RTGO=0
C**** TG3M must be set from a previous ML run?


      MODULE OCEAN_COM
!@sum  OCEAN_COM defines the model variables relating to the ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      USE E001M12_COM, only : im,jm
      USE OCEAN, only : lmom
      IMPLICIT NONE

!@var TG3M Monthly accumulation of temperatures at base of mixed layer
      REAL*8, DIMENSION(IM,JM,12) :: TG3M
!@var RTGO Temperature anomaly in thermocline
      REAL*8, DIMENSION(LMOM,IM,JM) :: RTGO
!@var STG3 accumulated temperature at base of mixed layer
      REAL*8, DIMENSION(IM,JM) :: STG3
!@var DTG3 accumulated temperature diff. from initial monthly values
      REAL*8, DIMENSION(IM,JM) :: DTG3
!@var EDO ocean vertical diffusion (m^2/s)
      REAL*8, DIMENSION(IM,JM) :: EDO
!@var DZ thermocline layer thickness (m)
      REAL*8, DIMENSION(LMOM) :: DZ
!@var DZO,BYDZO distance between centres in thermocline layer (m)
      REAL*8, DIMENSION(LMOM-1) :: DZO,BYDZO

      END MODULE OCEAN_COM

      SUBROUTINE ODIFS
!@sum  ODFIS calculates heat diffusion at the base of the mixed layer
!@auth Gary Russell
!@ver  1.0
!@calls ODFFUS
C****
C**** THIS SUBROUTINE CALCULATES THE ANNUAL OCEAN TEMPERATURE AT THE
C**** MAXIMUM MIXED LAYER, COMPARES THAT TO THE CONTROL RUN'S
C**** TEMPERATURE, CALLS SUBROUTINE ODFFUS, AND REDUCES THE UPPER
C**** OCEAN TEMPERATURES BY THE AMOUNT OF HEAT THAT IS DIFFUSED INTO
C**** THE THERMOCLINE
C****
      USE CONSTANT, only : sday
      USE E001M12_COM, only : im,jm,focean,jmon,jday,jdate,itocean
     *     ,itoice
      USE GEOM, only : imaxj
      USE OCEAN_COM, only : tg3m,rtgo,stg3,dtg3,edo,dz,dzo,bydzo
      USE OCEAN, only : z12o,lmom,tocean
      USE SEAICE_COM, only : rsi
      USE DAGCOM, only : aj,j_ftherm
      USE FLUXES, only : gtemp
      USE FILEMANAGER
      IMPLICIT NONE

      REAL*8 :: ADTG3
      INTEGER I,J,L,IMAX,IFIRST,iu_EDDY
      REAL*8, PARAMETER :: PERDAY=1./365d0
!@param ALPHA degree of implicitness (1 fully implicit,0 fully explicit)
      REAL*8, PARAMETER :: ALPHA=.5d0
!@param FAC ratio of adjacent deep ocean layers so total depth is 1000m
C**** NOTE: This assumes that LMOM is 9. For any different number of
C**** layers, the equation 1000=(1-x^(LMOM-1))/(1-x) should be solved.
      REAL*8, PARAMETER :: FAC=1.705357255658901d0
      DATA IFIRST/1/
C****
C**** READ IN EDDY DIFFUSIVITY AT BASE OF MIXED LAYER
C****
      IF (IFIRST.EQ.1) THEN
        CALL getunit("EDDY",iu_EDDY,.TRUE.,.TRUE.)
        CALL READT (iu_EDDY,0,EDO,IM*JM,EDO,1)
C**** DEFINE THE VERTICAL LAYERING EVERYWHERE EXCEPT LAYER 1 THICKNESS
        DZ(2)=10.
        DZO(1)=0.5*DZ(2)   ! 10./SQRT(FAC)
        BYDZO(1)=1./DZO(1)
        DO L=2,LMOM-1
          DZ(L+1)=DZ(L)*FAC
          DZO(L)=0.5*(DZ(L+1)+DZ(L))     !DZO(L-1)*FAC
          BYDZO(L)=1./DZO(L)
        END DO
        IFIRST=0
      END IF
C****
C**** ACCUMULATE OCEAN TEMPERATURE AT MAXIMUM MIXED LAYER
C****
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          STG3(I,J)=STG3(I,J)+TOCEAN(3,I,J)
        END DO
      END DO
C****
C**** AT THE END OF EACH MONTH, UPDATE THE OCEAN TEMPERATURE
C**** DIFFERENCE AND REPLACE THE MONTHLY SUMMED TEMPERATURE
C****
      IF(JDATE.EQ.1) THEN
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          DTG3(I,J)=DTG3(I,J)+(STG3(I,J)-TG3M(I,J,JMON))
          TG3M(I,J,JMON)=STG3(I,J)
          STG3(I,J)=0.
        END DO
      END DO
      END IF
C****
C**** DIFFUSE THE OCEAN TEMPERATURE DIFFERENCE OF THE UPPER LAYERS
C**** INTO THE THERMOCLINE AND REDUCE THE UPPER TEMPERATURES BY THE
C**** HEAT THAT IS DIFFUSED DOWNWARD
C****
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF(FOCEAN(I,J).GT.0.) THEN

            ADTG3=DTG3(I,J)*PERDAY
            RTGO(1,I,J)=ADTG3
C**** Set first layer thickness
            DZ(1)=Z12O(I,J)

            CALL ODFFUS (SDAY,ALPHA,EDO(I,J),DZ,BYDZO,RTGO(1,I,J),LMOM)

            DO L=1,3
              TOCEAN(L,I,J)=TOCEAN(L,I,J)+(RTGO(1,I,J)-ADTG3)
            END DO
            AJ(J,J_FTHERM,ITOCEAN)=AJ(J,J_FTHERM,ITOCEAN)-(RTGO(1,I,J)
     *           -ADTG3)*Z12O(I,J)*FOCEAN(I,J)*(1.-RSI(I,J))
            AJ(J,J_FTHERM,ITOICE )=AJ(J,J_FTHERM,ITOICE )-(RTGO(1,I,J)
     *           -ADTG3)*Z12O(I,J)*FOCEAN(I,J)*RSI(I,J)
            GTEMP(1:2,1,I,J) = TOCEAN(1:2,I,J)
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE ODIFS

      SUBROUTINE ODFFUS (DT,ALPHA,ED,DZ,BYDZO,R,LMIJ)
!@sum  ODFFUS calculates the vertical mixing of a tracer
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
!@calls TRIDIAG
      IMPLICIT NONE
!@var LMIJ IS THE NUMBER OF VERTICAL LAYERS
      INTEGER, INTENT(IN) :: LMIJ
!@var ED diffusion coefficient between adjacent layers (m**2/s)
!@var ALPHA determines the time scheme (0 explicit,1 fully implicit)
!@var DT time step (s)
      REAL*8, INTENT(IN) :: ED,ALPHA,DT
!@var DZ the depth of the layers (m)
!@var BYDZO is the inverse of depth between layer centers (1/m)
      REAL*8, INTENT(IN) :: DZ(LMIJ),BYDZO(LMIJ-1)
!@var R tracer concentration
      REAL*8, INTENT(INOUT) :: R(LMIJ)

      REAL*8 AM(LMIJ),BM(LMIJ),CM(LMIJ),DM(LMIJ)
      INTEGER L
C**** SET UP TRIDIAGONAL MATRIX ENTRIES AND RIGHT HAND SIDE
      AM(1)=0
      BM(1)=DZ(1)+ALPHA*DT*ED*BYDZO(1)
      CM(1)=     -ALPHA*DT*ED*BYDZO(1)
      DM(1)=DZ(1)*R(1)-(1.-ALPHA)*DT*ED*(R(1)-R(2))*BYDZO(1)

      DO L=2,LMIJ-1
        AM(L)=     -ALPHA*DT* ED*BYDZO(L-1)
        BM(L)=DZ(L)+ALPHA*DT*(ED*BYDZO(L-1)+ED*BYDZO(L))
        CM(L)=     -ALPHA*DT*               ED*BYDZO(L)
        DM(L)=DZ(L)*R(L)+(1.-ALPHA)*DT*(ED*(R(L-1)-R(L))*BYDZO(L-1)
     *                                 -ED*(R(L)-R(L+1))*BYDZO(L))
      END DO

      AM(LMIJ)=        -ALPHA*DT*ED*BYDZO(LMIJ-1)
      BM(LMIJ)=DZ(LMIJ)+ALPHA*DT*ED*BYDZO(LMIJ-1)
      CM(LMIJ)=0.
      DM(LMIJ)=DZ(LMIJ)*R(LMIJ)+(1.-ALPHA)*DT*ED*
     *         (R(LMIJ-1)-R(LMIJ))*BYDZO(LMIJ-1)

      CALL TRIDIAG(AM,BM,CM,DM,R,LMIJ)

      RETURN
      END SUBROUTINE ODFFUS
