      MODULE OCEAN
!@sum  OCEAN contains all the ocean subroutines
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
!@cont OSTRUC,OCLIM,init_OCEAN,vflx_OCEAN,daily_OCEAN,PREC_OC
      USE CONSTANT, only : lhm,rhow,rhoi,shw,shi
      USE E001M12_COM, only : im,jm,lm,gdata,focean,flake,fland,fearth
     *     ,flice,kocean,Itime,jmon,jdate,jday,JDendOfM,JDmidOfM
      USE PBLCOM
     &     , only : npbl=>n,uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,
     *     cq=>cqgs,ipbl
      USE GEOM
      USE SEAICE_COM, only : rsi,msi
      USE SEAICE, only : xsi1,xsi2,xsi3,xsi4,ace1i,z1i

      IMPLICIT NONE
!@param LMOM number of layers for deep ocean diffusion
      INTEGER, PARAMETER :: LMOM = 9
!@var OA generic array for ocean heat transport calculations
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

!@var Z1O ocean mixed layer depth
      REAL*8, SAVE,DIMENSION(IM,JM) :: Z1O
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
      REAL*8 Z1OOLD
      COMMON/WORK2/Z1OOLD(IM,JM)
C****
C**** FLAND     LAND COVERAGE (1)
C****
C**** TOCEAN 1  OCEAN TEMPERATURE OF FIRST LAYER (C)
C****        2  MEAN OCEAN TEMPERATURE OF SECOND LAYER (C)
C****        3  OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)
C****     RSI   RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****     MSI   OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)
C****
C****
C**** RESTRUCTURE OCEAN LAYERS
C****
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF (FLAND(I,J).GE.1.) CYCLE
          IF (Z1OOLD(I,J).GE.Z12O(I,J)) GO TO 140
          IF (Z1O(I,J).EQ.Z1OOLD(I,J)) CYCLE
          WTR1O=RHOW*Z1O(I,J)-RSI(I,J)*(GDATA(I,J,1)
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

      REAL*8 Z1OOLD
      COMMON/WORK2/Z1OOLD(IM,JM)
      REAL*8 XZO(IM,JM),XZN(IM,JM)
      INTEGER, INTENT(IN) :: IDOZ1O

      INTEGER n,MD,J,I,LSTMON,K,MDMAX,IMAX,IMON
      REAL*8 PLICEN,PLICE,POICE,POCEAN,RSICSQ,ZIMIN,ZIMAX,X1
     *     ,X2,Z1OMIN,RSINEW,TIME,FRAC
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
     *           CRSI(I,J)**2/12.)
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
          TOCEAN(1,I,J) = AOST(I,J) + BOST(I,J)*TIME +
     +         COST(I,J)*(TIME*TIME - 1./12.d0)            !-BY12
          IF(KRSI(I,J)) 410,430,420
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
 410      IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .gt. 0.)  then
            RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5) !  TIME < T0
          ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .gt. 0.)  then
            RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME) !  T1 < TIME
          ELSE
            RSINEW = 0.         !  T0 < TIME < T1
          END IF
          GO TO 440
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
 420      IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .lt. 1.)  then
            RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5) !  TIME < T0
          ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .lt. 1.)  then
            RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME) !  T1 < TIME
          ELSE
            RSINEW = 1.         !  T0 < TIME < T1
          END IF
          GO TO 440
C**** RSI uses quadratic fit
 430      RSINEW = ARSI(I,J) + BRSI(I,J)*TIME +
     +         CRSI(I,J)*(TIME*TIME - 1./12.d0)    !-BY12
 440      RSI(I,J)=RSINEW
          MSI(I,J)=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONABLE VALUE
          IF (TOCEAN(1,I,J).LT.-1.8) TOCEAN(1,I,J)=TFO
C**** REDUCE THE RATIO OF OCEAN ICE TO WATER BY .1*RHOI/ACEOI
c     IF (RSI(I,J).GT.0.) THEN
c        BYZICE=RHOI/(Z1I*RHOI+MSI(I,J))
c        RSI(I,J)=RSI(I,J)*(1.-.06*(BYZICE-1./5.))
c     END IF
C**** ZERO OUT SNOWOI, TG1OI, TG2OI IF THERE IS NO OCEAN ICE
          IF (RSI(I,J).LE.0.) THEN
            GDATA(I,J,1)=0.
            GDATA(I,J,3)=0.
            GDATA(I,J,7)=0.
            GDATA(I,J,15)=0.
            GDATA(I,J,16)=0.
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE
      DO I=2,IM
        GDATA(I,JM,1)=GDATA(1,JM,1)
        GDATA(I,JM,3)=GDATA(1,JM,3)
        GDATA(I,JM,7)=GDATA(1,JM,7)
        GDATA(I,JM,15)=GDATA(1,JM,15)
        GDATA(I,JM,16)=GDATA(1,JM,16)
        TOCEAN(1,I,JM)=TOCEAN(1,1,JM)
        RSI(I,JM)=RSI(1,JM)
        MSI(I,JM)=MSI(1,JM)
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
      Z1OMIN=1.+(RHOI*Z1I+GDATA(I,J,1)+MSI(I,J))/RHOW
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
      GDATA(I,J,12)=(GDATA(I,J,12)*PLICE+GDATA(I,J,1)*POICE)/PLICEN
      GDATA(I,J,13)=(GDATA(I,J,13)*PLICE+
     *  (GDATA(I,J,3)*XSI1+GDATA(I,J,7)*XSI2)*POICE+
     *  (LHM+SHW*TOCEAN(1,I,J))*POCEAN/SHI)/PLICEN
      GDATA(I,J,14)=(GDATA(I,J,14)*PLICE+
     *  (GDATA(I,J,15)*XSI3+GDATA(I,J,16)*XSI4)*POICE+
     *  (LHM+SHW*TOCEAN(1,I,J))*POCEAN/SHI)/PLICEN
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

      SUBROUTINE PREC_OC(TGW,WTRO,PRCP,ENRGP,ERUN4)
!@sum  PREC_OC Adds the precipitation to ocean arrays
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
      REAL*8 TGW, PRCP, WTRO, ENRGP, ERUN4, ENRGO

      ENRGO=ENRGP-ERUN4
      TGW=TGW+(ENRGO/(WTRO*SHW))

      RETURN
      END SUBROUTINE PREC_OC


      SUBROUTINE OSOURC (ROICE,MSI1,MSI2,TGW,WTRO,EIW0,OTDT,ENRGO,ERUN4
     *     ,ENRGFO,ACEFO,ACE2F,WTRW0,ENRGW0,ENRGFI,F2DT,TFW)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@var TFW freezing temperature for water underlying ice (C)
      REAL*8, INTENT(IN) :: TFW

      REAL*8, INTENT(IN) :: ROICE, MSI1, MSI2
!@var TGW and WTRO mixed layer temp.(C) and mass (kg/m^2)
      REAL*8 TGW, WTRO, ENRGO
      REAL*8 EIW0, ENRGIW, WTRI1, EFIW, ENRGO0, EOFRZ
      REAL*8 WTRW0, ENRGW0, WTRW, ENRGW
      REAL*8, INTENT(IN) :: F2DT, ERUN4, OTDT
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACE2F
C**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACE2F=0. ; ENRGFI=0.

C**** Calculate energy in mixed layer
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
c 80      ACEFO=(ENRGO0+ENRGO-EOFRZ)/(TFW*(SHI-SHW)-LHM)
 80     ACEFO=(ENRGO0+ENRGO-EOFRZ)/(-LHM)
        ENRGFO=ACEFO*(TFW*SHI-LHM)
        IF (ROICE.LE.0.) THEN
          TGW=TFW
          RETURN
        END IF
      END IF
C****
      IF (ROICE.le.0) RETURN

C**** CALCULATE THE ENERGY OF THE WATER BELOW THE ICE AT FREEZING POINT
C**** AND CHECK WHETHER NEW ICE MUST BE FORMED
      ENRGIW = F2DT+OTDT-ERUN4 ! heat flux to the ocean under ice
      WTRI1 = WTRO-(MSI1+MSI2) ! new mass of ocean (kg/m^2)
      EFIW = WTRI1*TFW*SHW ! freezing energy of ocean mass WTRI1
      IF (EIW0+ENRGIW .LE. EFIW) THEN ! freezing case
C**** FLUXES WOULD COOL TGW TO FREEZING POINT AND FREEZE SOME MORE ICE
c       ACE2F = (EIW0+ENRGIW-EFIW)/(TFW*(SHI-SHW)-LHM) !
        ACE2F = (EIW0+ENRGIW-EFIW)/(-LHM) ! mass freezes under the ice
        ENRGFI = ACE2F*(TFW*SHI-LHM) ! energy of frozen ice
      END IF
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      WTRW = WTRW0-(1.-ROICE)*ACEFO-ROICE*ACE2F ! new ocean mass
      ENRGW = ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)
      TGW = ENRGW/(WTRW*SHW) ! new ocean temperature

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
      USE E001M12_COM, only : im,jm,fland,flice,gdata,kocean
      USE OCEAN, only : ota,otb,otc,z12o,dm,iu_osst,iu_sice,iu_ocnml
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
      USE CONSTANT, only : rhow,shw
      USE E001M12_COM, only : im,jm,kocean,focean,gdata
      USE OCEAN, only : tocean,ostruc,oclim,tfo,z1O
      USE DAGCOM, only : aij,ij_toc2,ij_tgo2
      USE SEAICE_COM, only : rsi,msi
      USE SEAICE, only : simelt,ace1i
      USE GEOM, only : imaxj
      IMPLICIT NONE
      INTEGER I,J,IEND,IMAX
      REAL*8 MSI1, MSI2, ROICE, SNOW, TG1, TG2, TG3, TG4, TGW, WTRO,
     *     PWATER, ACE, WTRW, ENRGW, ENRGUSED

C**** update ocean related climatologies
      CALL OCLIM(IEND)

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
C**** Only melting of ocean ice (no lakes) (?)
            IF (FOCEAN(I,J)*RSI(I,J) .LE. 0.) CYCLE
C**** REDUCE ICE EXTENT IF OCEAN TEMPERATURE IS GREATER THAN ZERO
C**** (MELTING POINT OF ICE)
            TGW=TOCEAN(1,I,J)
            IF (TGW.LE.0.) CYCLE
            ROICE=RSI(I,J)
            MSI2=MSI(I,J)
            SNOW=GDATA(I,J,1)   ! snow mass
            MSI1 = SNOW + ACE1I
            TG1 = GDATA(I,J,3)  ! first layer sea ice temperature
            TG2 = GDATA(I,J,7)  ! second layer sea ice temperature
            TG3 = GDATA(I,J,15) ! third layer sea ice temperature
            TG4 = GDATA(I,J,16) ! fourth layer sea ice temperature
            PWATER = FOCEAN(I,J) ! ocean fraction
            ACE = SNOW + ACE1I + MSI2
            WTRO=Z1O(I,J)*RHOW
            WTRW=WTRO-ROICE*ACE
            ENRGW=WTRW*TGW*SHW  ! energy of water available for melting
            CALL SIMELT(ROICE,SNOW,MSI1,MSI2,ACE,TG1,TG2,TG3,TG4
     *           ,ENRGW,PWATER,ENRGUSED)
C**** RESAVE PROGNOSTIC QUANTITIES
            TGW=(ENRGW-ENRGUSED)/(WTRO*SHW)
            TOCEAN(1,I,J)=TGW
            IF (ROICE.le.0.) THEN
              RSI(I,J)=0.
              MSI(I,J)=0.
              GDATA(I,J,1)=0.
              GDATA(I,J,3)=0.
              GDATA(I,J,7)=0.
              GDATA(I,J,15)=0.
              GDATA(I,J,16)=0.
            ELSE
              RSI(I,J)=ROICE
              MSI(I,J)=MSI2
              GDATA(I,J,1)=SNOW
              GDATA(I,J,3)=TG1
              GDATA(I,J,7)=TG2
              GDATA(I,J,15)=TG3
              GDATA(I,J,16)=TG4
            END IF
          END DO
        END DO
      END IF

      RETURN
      END

      SUBROUTINE vflx_OCEAN
!@sum  vflx_OCEAN saves quantities for OHT calculations
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : RHOI
      USE E001M12_COM, only : IM,JM,GDATA
      USE OCEAN, only : OA
      USE SEAICE, only : Z1I,XSI1,XSI2,XSI3,XSI4
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
      END SUBROUTINE vflx_OCEAN

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
      END SUBROUTINE io_ocean

