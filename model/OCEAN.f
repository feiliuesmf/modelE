#include "rundeck_opts.h"

      MODULE STATIC_OCEAN
!@sum  STATIC_OCEAN contains the ocean subroutines common to all Q-flux
!@+    and fixed SST runs
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
!@cont OSTRUC,OCLIM,init_OCEAN,daily_OCEAN,DIAGCO
!@+    PRECIP_OC,OCEANS
      USE CONSTANT, only : lhm,rhow,rhoi,shw,shi,by12,byshi
      USE MODEL_COM, only : im,jm,lm,focean,fland,fearth,flice
     *     ,Iyear1,Itime,jmon,jdate,jday,jyear,jmpery,JDendOfM,JDmidOfM
     *     ,ItimeI,kocean,itocean,itoice
      USE GEOM
      USE PBLCOM, only : npbl,uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs
     *     ,cq=>cqgs,ipbl,roughl
      USE SEAICE, only : xsi,ace1i,z1i,ac2oim,z2oim,ssi0,tfrez,fleadoc
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
#ifdef TRACERS_WATER
     *     ,trsi,trsi0,ntm
#endif
      USE LANDICE_COM, only : snowli,tlandi
      USE FLUXES, only : gtemp,sss
      USE DAGCOM, only : aij, ij_smfx, aj, j_implh, j_implm
      IMPLICIT NONE
      SAVE

!@dbparam ocn_cycl determines whether prescribed ocean data repeat
!@+       after 1 year - 0:no 1:yes 2:seaice repeats,sst does not
      integer :: ocn_cycl = 1
!@var TOCEAN temperature of the ocean (C)
      REAL*8, DIMENSION(3,IM,JM) :: TOCEAN
!@var SSS0 default sea surface salinity (psu)
      REAL*8 :: SSS0=34.7d0

!@var OTA,OTB,OTC ocean heat transport coefficients
      REAL*8, DIMENSION(IM,JM,4) :: OTA,OTB
      REAL*8, DIMENSION(IM,JM) :: OTC
!@var SINANG,COSANG Fourier coefficients for heat transports
      REAL*8 :: SINANG,SN2ANG,SN3ANG,SN4ANG,
     *          COSANG,CS2ANG,CS3ANG,CS4ANG

!@var Z1O ocean mixed layer depth
      REAL*8, DIMENSION(IM,JM) :: Z1O
!@var Z1OOLD previous ocean mixed layer depth
      REAL*8, DIMENSION(IM,JM) :: Z1OOLD
!@var Z12O annual maximum ocean mixed layer depth
      REAL*8, DIMENSION(IM,JM) :: Z12O

      REAL*8, DIMENSION(IM,JM) :: DM,AOST,EOST1,EOST0,BOST
     *     ,COST,ARSI,ERSI1,ERSI0,BRSI,CRSI
      INTEGER, DIMENSION(IM,JM) ::  KRSI
      COMMON/OOBS/DM,AOST,EOST1,EOST0,BOST,COST,ARSI,ERSI1,ERSI0,BRSI
     *     ,CRSI,KRSI

!@var iu_OSST,iu_SICE,iu_OCNML unit numbers for climatologies
      INTEGER iu_OSST,iu_SICE,iu_OCNML

      CONTAINS

      SUBROUTINE OSTRUC(QTCHNG)
!@sum  OSTRUC restructures the ocean temperature profile as ML
!@sum         depths are changed (generally once a day)
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 TGAVE,DWTRO,WTR1O,WTR2O
!@var QTCHNG true if TOCEAN(1) is changed (needed for qflux calculation)
      LOGICAL, INTENT(IN) :: QTCHNG
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
        DO I=1,IMAXJ(J)
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
 120      IF (QTCHNG) THEN  ! change TOCEAN(1)
            TGAVE=(TOCEAN(2,I,J)*DWTRO+(2.*TOCEAN(2,I,J)-TOCEAN(3,I,J))
     *           *WTR2O)/(WTR2O+DWTRO)
            TOCEAN(1,I,J)=TOCEAN(1,I,J)+((TGAVE-TOCEAN(1,I,J))
     *           *DWTRO/WTR1O)
          END IF
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

      SUBROUTINE OCLIM(end_of_day)
!@sum OCLIM calculates daily ocean data from ocean/sea ice climatologies
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean or fixed SST)
      IMPLICIT NONE

      REAL*8, SAVE :: XZO(IM,JM),XZN(IM,JM)
      LOGICAL, INTENT(IN) :: end_of_day

      INTEGER n,J,I,LSTMON,K,m,m1
      REAL*8 PLICEN,PLICE,POICE,POCEAN,RSICSQ,ZIMIN,ZIMAX,X1
     *     ,X2,Z1OMIN,RSINEW,TIME,FRAC,MSINEW,OPNOCN,TFO
!@var JDLAST julian day that OCLIM was last called
      INTEGER, SAVE :: JDLAST=0
!@var IMON0 current month for SST climatology reading
      INTEGER, SAVE :: IMON0 = 0
!@var IMON current month for ocean mixed layer climatology reading
      INTEGER, SAVE :: IMON = 0

      IF (KOCEAN.EQ.1) GO TO 500
      if(.not.(end_of_day.or.itime.eq.itimei)) return
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
      IF (JMON.EQ.IMON0) GO TO 400
      IF (IMON0.EQ.0) THEN
C****   READ IN LAST MONTH'S END-OF-MONTH DATA
        if (ocn_cycl.ge.1) then
          LSTMON=JMON-1
          if (lstmon.eq.0) lstmon = 12
          call readt (iu_SICE,IM*JM,ERSI0,IM*JM,ERSI0,LSTMON)
          if (ocn_cycl.eq.1) then
            call readt (iu_OSST,IM*JM,EOST0,IM*JM,EOST0,LSTMON)
          else ! if (ocn_cycl.eq.2) then
            LSTMON=JMON-1+(JYEAR-IYEAR1)*JMperY
  290       read (iu_OSST) M
            if (m.lt.lstmon) go to 290
            backspace iu_OSST
            CALL MREAD (iu_OSST, m,IM*JM,EOST0,IM*JM,EOST0)
            WRITE(6,*) 'Read End-of-month ocean data from ',JMON-1,M
            IF(M.NE.LSTMON)
     &         call stop_model('Read error: ocean data',255)
          end if
        else   !  if (ocn_cycl.eq.0) then
          LSTMON=JMON-1+(JYEAR-IYEAR1)*JMperY
  300     read (iu_OSST) M
          if (m.lt.lstmon) go to 300
          backspace iu_OSST
  310     read (iu_SICE) m
          if (m.lt.lstmon) go to 310
          backspace iu_SICE
          CALL MREAD (iu_OSST, m,IM*JM,EOST0,IM*JM,EOST0)
          CALL MREAD (iu_SICE,m1,IM*JM,ERSI0,IM*JM,ERSI0)
          WRITE(6,*) 'Read End-of-month ocean data from ',JMON-1,M,M1
          IF(M.NE.M1.OR.M.NE.LSTMON)
     &         call stop_model('Read error: ocean data',255)
        end if
      ELSE
C****   COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
        EOST0=EOST1
        ERSI0=ERSI1
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      IMON0=JMON
      if (ocn_cycl.ge.1) then
        if (jmon.eq.1) then
          if (ocn_cycl.eq.1) rewind iu_OSST
          rewind iu_SICE
          read (iu_SICE)  ! skip over DM-record
        end if
        call readt (iu_SICE,0,ARSI,2*im*jm,ARSI,1) ! reads ARSI,ERSI1
        if (ocn_cycl.eq.1) then
          call readt (iu_OSST,0,AOST,2*im*jm,AOST,1) ! reads AOST,EOST1
        else  ! if (ocn_cycl.eq.2) then
          CALL MREAD (iu_OSST, M,0,AOST,2*IM*JM,AOST) ! READS AOST,EOST1
          WRITE(6,*) 'Read in ocean data for month',JMON,M
          IF(JMON.NE.MOD(M-1,12)+1)
     &       call stop_model('Error: Ocean data',255)
        end if
      else   !  if (ocn_cycl.eq.0) then
        CALL MREAD (iu_OSST, M,0,AOST,2*IM*JM,AOST) ! READS AOST,EOST1
        CALL MREAD (iu_SICE,M1,0,ARSI,2*IM*JM,ARSI) ! READS ARSI,ERSI1
        WRITE(6,*) 'Read in ocean data for month',JMON,M,M1
        IF(M.NE.M1.OR.JMON.NE.MOD(M-1,12)+1)
     &       call stop_model('Error: Ocean data',255)
      end if
C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO J=1,JM
        DO I=1,IMAXJ(J)
          BOST(I,J)=EOST1(I,J)-EOST0(I,J)
          COST(I,J)=3.*(EOST1(I,J)+EOST0(I,J)) - 6.*AOST(I,J)
          BRSI(I,J)=0.
          CRSI(I,J)=0.          ! extreme cases
          KRSI(I,J)=0           ! ice constant
          IF(ARSI(I,J).LE.0. .or. ARSI(I,J).GE.1.) CYCLE
          BRSI(I,J)=ERSI1(I,J)-ERSI0(I,J) ! quadratic fit
          CRSI(I,J)=3.*(ERSI1(I,J)+ERSI0(I,J)) - 6.*ARSI(I,J)
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
        ZIMIN=Z1I+Z2OIM
        ZIMAX=2d0
        IF(J.GT.JM/2) ZIMAX=3.5d0
        DO I=1,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
C**** OST always uses quadratic fit
            TOCEAN(1,I,J)=AOST(I,J)+BOST(I,J)*TIME
     *                   +COST(I,J)*(TIME**2-BY12)
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
                RSINEW = 1.     !  T0 < TIME < T1
            END IF
C**** RSI uses quadratic fit
            CASE (0)
              RSINEW=ARSI(I,J)+BRSI(I,J)*TIME+CRSI(I,J)*(TIME**2-BY12)
            END SELECT
C**** Set new mass
            MSINEW=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** Ensure that lead fraction is consistent with kocean=1 case
            IF (RSINEW.gt.0) THEN
              OPNOCN=MIN(0.1d0,FLEADOC*RHOI/(RSINEW*(ACE1I+MSINEW)))
              IF (RSINEW.GT.1.-OPNOCN) THEN
                RSINEW = 1.-OPNOCN
                MSINEW=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
              END IF
            END IF
C**** accumulate diagnostics
            IF(MSI(I,J).eq.0) MSI(I,J)=MSINEW
            IF (end_of_day) THEN
              AIJ(I,J,IJ_SMFX)=AIJ(I,J,IJ_SMFX)+
     *           (SNOWI(I,J)+ACE1I)*(RSINEW-RSI(I,J))+
     *           RSINEW*MSINEW-RSI(I,J)*MSI(I,J)
              AJ(J,J_IMPLM,ITOICE)=AJ(J,J_IMPLM,ITOICE)-
     *             (MSINEW-MSI(I,J))*RSI(I,J)*FOCEAN(I,J)
              AJ(J,J_IMPLH,ITOICE)=AJ(J,J_IMPLH,ITOICE)-SUM(HSI(3:4,I,J)
     *             )*(MSINEW/MSI(I,J)-1.)*RSI(I,J)*FOCEAN(I,J)
              AJ(J,J_IMPLM,ITOCEAN)=AJ(J,J_IMPLM,ITOCEAN)-(1.-RSI(I,J))
     *             *FOCEAN(I,J)*(RSINEW-RSI(I,J))*(MSINEW+ACE1I+SNOWI(I
     *             ,J))
              AJ(J,J_IMPLH,ITOCEAN)=AJ(J,J_IMPLH,ITOCEAN)-(1.-RSI(I,J))
     *             *FOCEAN(I,J)*(RSINEW-RSI(I,J))*SUM(HSI(1:4,I,J))
            END IF
C**** adjust enthalpy and salt so temperature/salinity remain constant
            HSI(3:4,I,J)=HSI(3:4,I,J)*(MSINEW/MSI(I,J))
            SSI(3:4,I,J)=SSI(3:4,I,J)*(MSINEW/MSI(I,J))
#ifdef TRACERS_WATER
            TRSI(:,3:4,I,J)=TRSI(:,3:4,I,J)*(MSINEW/MSI(I,J))
#endif
            RSI(I,J)=RSINEW
            MSI(I,J)=MSINEW
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONABLE VALUE
            TFO=tfrez(sss(i,j))
            IF (TOCEAN(1,I,J).LT.TFO) TOCEAN(1,I,J)=TFO
C**** SET DEFAULTS IF NO OCEAN ICE
            IF (RSI(I,J).LE.0.) THEN
              HSI(1:2,I,J)=(SHI*TFO-LHM)*XSI(1:2)*ACE1I
              HSI(3:4,I,J)=(SHI*TFO-LHM)*XSI(3:4)*AC2OIM
              SSI(1:2,I,J)=SSI0*XSI(1:2)*ACE1I
              SSI(3:4,I,J)=SSI0*XSI(3:4)*AC2OIM
#ifdef TRACERS_WATER
              DO N=1,NTM
                TRSI(N,1:2,I,J)=TRSI0(N)*(1.-SSI0)*XSI(1:2)*ACE1I
                TRSI(N,3:4,I,J)=TRSI0(N)*(1.-SSI0)*XSI(3:4)*AC2OIM
              END DO
#endif
              SNOWI(I,J)=0.
              GTEMP(1:2,2,I,J)=TFO
            END IF
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE (for prescribed data only)
      IF (FOCEAN(1,JM).gt.0) THEN
        DO I=2,IM
          SNOWI(I,JM)=SNOWI(1,JM)
          TOCEAN(1,I,JM)=TOCEAN(1,1,JM)
          RSI(I,JM)=RSI(1,JM)
          MSI(I,JM)=MSI(1,JM)
          HSI(:,I,JM)=HSI(:,1,JM)
          SSI(:,I,JM)=SSI(:,1,JM)
#ifdef TRACERS_WATER
          TRSI(:,:,I,JM)=TRSI(:,:,1,JM)
#endif
          GTEMP(1:2,2,I,JM)=GTEMP(1:2,2,1,JM)
        END DO
      END IF
      RETURN
C****
C**** CALCULATE DAILY OCEAN MIXED LAYER DEPTHS FROM CLIMATOLOGY
C****
C**** SAVE PREVIOUS DAY'S MIXED LAYER DEPTH
 500  Z1OOLD=Z1O
C**** COMPUTE Z1O ONLY AT THE END OF A DAY OR AT ItimeI
      IF (.not.(end_of_day.or.itime.eq.itimei)) RETURN
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
            call stop_model(
     &           'ERROR READING IN SETTING OCEAN CLIMATOLOGY',255)
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

C**** INTERPOLATE OCEAN DATA TO CURRENT DAY
      FRAC = REAL(JDmidOFM(IMON)-JDAY,KIND=8)/
     a           (JDmidOFM(IMON)-JDmidOFM(IMON-1))
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
C**** Always define a new roughness length (to prevent restart problems)
      ROUGHL(I,J)=1.84d0 ! typical for Antarctica
C****
      FLAND(I,J)=1.
      FLICE(I,J)=PLICEN
      FOCEAN(I,J)=0.
C**** set gtemp array
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
      RETURN
      END SUBROUTINE OCLIM

      SUBROUTINE OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,FODT,FIDT,RVRRUN
     *     ,RVRERUN,SIMELT,ESIMELT,EVAPO,EVAPI,TFW,RUN4O,ERUN4O,RUN4I
     *     ,ERUN4I,ENRGFO,ACEFO,ACEFI,ENRGFI)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@var TFW freezing temperature for water underlying ice (C)
      REAL*8, INTENT(IN) :: TFW

!@var TGW mixed layer temp.(C)
      REAL*8, INTENT(INOUT) :: TGW
      REAL*8, INTENT(IN) :: ROICE, SMSI, WTRO, EVAPO, EVAPI, RUN0
      REAL*8, INTENT(IN) :: FODT, FIDT, OTDT, RVRRUN, RVRERUN, SIMELT,
     *     ESIMELT
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACEFI, RUN4O, RUN4I
     *     , ERUN4O, ERUN4I
      REAL*8 EIW0, ENRGIW, WTRI1, EFIW, ENRGO0, EOFRZ, ENRGO
      REAL*8 WTRW0, ENRGW0, WTRW, ENRGW, WTRI0, SMSI0
C**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACEFI=0. ; ENRGFI=0.

C**** Calculate previous ice mass (before fluxes applied)
      SMSI0=SMSI+RUN0+EVAPI

C**** Calculate effect of lateral melt of sea ice
      IF (SIMELT.gt.0) THEN
        TGW=TGW + (ESIMELT-SIMELT*TGW*SHW)/
     *       (SHW*(1-ROICE)*WTRO+ROICE*(WTRO-SMSI0))
      END IF

C**** Calculate extra mass flux to ocean, balanced by deep removal
      RUN4O=-EVAPO+RVRRUN  ! open ocean
      RUN4I=-EVAPI+RVRRUN  ! under ice
      ERUN4O=RUN4O*TGW*SHW ! corresponding heat flux at bottom (open)
      ERUN4I=RUN4I*TGW*SHW !                              (under ice)

C**** Calculate heat fluxes to ocean
      ENRGO  = FODT+OTDT+RVRERUN-ERUN4O ! in open water
      ENRGIW = FIDT+OTDT+RVRERUN-ERUN4I ! under ice

C**** Calculate energy in mixed layer (open ocean)
      ENRGO0=WTRO*TGW*SHW
      EOFRZ=WTRO*TFW*SHW
      IF (ENRGO0+ENRGO.GE.EOFRZ) THEN
C**** FLUXES RECOMPUTE TGW WHICH IS ABOVE FREEZING POINT FOR OCEAN
        IF (ROICE.LE.0.) THEN
          TGW=TGW+ENRGO/(WTRO*SHW)
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
        ACEFI = (EIW0+ENRGIW-EFIW)/(TFW*(SHI-SHW)-LHM)
        ENRGFI = ACEFI*(TFW*SHI-LHM) ! energy of frozen ice
      END IF
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      WTRW = WTRO  -(1.-ROICE)* ACEFO        -ROICE*(SMSI+ACEFI)
      ENRGW= ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)
      TGW  = ENRGW/(WTRW*SHW) ! new ocean temperature

      RETURN
      END SUBROUTINE OSOURC

      END MODULE STATIC_OCEAN


      SUBROUTINE init_OCEAN(iniOCEAN)
!@sum init_OCEAN initiallises ocean variables
!@auth Original Development Team
!@ver  1.0
      USE FILEMANAGER
      USE PARAM
      USE MODEL_COM, only : im,jm,fland,flice,kocean,focean
     *     ,fearth,iyear1
#ifdef TRACERS_WATER
      USE TRACER_COM, only : trw0
      USE FLUXES, only : gtracer
#endif
      USE FLUXES, only : gtemp,sss
      USE SEAICE, only : qsfix
      USE SEAICE_COM, only : snowi
      USE STATIC_OCEAN, only : ota,otb,otc,z12o,dm,iu_osst,iu_sice
     *     ,iu_ocnml,tocean,ocn_cycl,sss0
      USE DAGCOM, only : npts,icon_OCE,conpt0
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      LOGICAL, INTENT(IN) :: iniOCEAN  ! true if starting from ic.
      CHARACTER CONPT(NPTS)*10
!@var iu_OHT,iu_MLMAX unit numbers for reading in input files
      INTEGER :: iu_OHT,iu_MLMAX
      INTEGER :: I,J

      if (kocean.eq.0) then
        call sync_param( "ocn_cycl", ocn_cycl ) ! 1:cycle data 0:dont
C****   set up unit numbers for ocean climatologies
        call openunit("OSST",iu_OSST,.true.,.true.)
        call openunit("SICE",iu_SICE,.true.,.true.)
        if (ocn_cycl.ne.1) then
          write(6,*) '********************************************'
          write(6,*) '* Make sure that IYEAR1 is consistent with *'
          write(6,*) '*    the ocean data files OSST and SICE    *'
          write(6,*) '********************************************'
          write(6,*) 'IYEAR1=',IYEAR1
        end if

C****   Read in constant factor relating RSI to MSI from sea ice clim.
        CALL READT (iu_SICE,0,DM,IM*JM,DM,1)

      else !  IF (KOCEAN.eq.1) THEN
C****   Set up unit number of observed mixed layer depth data
        call openunit("OCNML",iu_OCNML,.true.,.true.)

C****   DATA FOR QFLUX MIXED LAYER OCEAN RUNS
C****   read in ocean heat transport coefficients
        call openunit("OHT",iu_OHT,.true.,.true.)
        READ (iu_OHT) OTA,OTB,OTC
        WRITE(6,*) "Read ocean heat transports from OHT"
        call closeunit (iu_OHT)

C****   read in ocean max mix layer depth
        call openunit("MLMAX",iu_MLMAX,.true.,.true.)
        CALL READT (iu_MLMAX,0,Z12O,IM*JM,Z12O,1)
        call closeunit (iu_MLMAX)

C****   initialise deep ocean arrays if required
        call init_ODEEP(iniOCEAN)

C****   IF SNOWI(I,J)<0, THE OCEAN PART WAS CHANGED TO LAND ICE
C****   BECAUSE THE OCEAN ICE REACHED THE MAX MIXED LAYER DEPTH
        DO J=1,JM
        DO I=1,IM
          IF(SNOWI(I,J).GE.-1.) CYCLE
          FLICE(I,J)=1-FLAND(I,J)+FLICE(I,J)
          FLAND(I,J)=1.
          FEARTH(I,J)=1.-FLICE(I,J)
          FOCEAN(I,J)=0.
          WRITE(6,'(2I3,'' OCEAN WAS CHANGED TO LAND ICE'')') I,J
        END DO
        END DO

C*****  set conservation diagnostic for ocean heat
        CONPT=CONPT0
        QCON=(/ F, F, F, T, F, T, F, T, T, F, F/)
        CALL SET_CON(QCON,CONPT,"OCN HEAT","(10^6 J/M**2)   ",
     *       "(W/M^2)         ",1d-6,1d0,icon_OCE)

      END IF
C**** Set gtemp array for oceans
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          GTEMP(1:2,1,I,J)=TOCEAN(1:2,I,J)
          SSS(I,J) = SSS0
#ifdef TRACERS_WATER
          gtracer(:,1,i,j)=trw0(:)
#endif
        ELSE
          SSS(I,J) = 0.
        END IF
      END DO
      END DO
C**** keep salinity in sea ice constant for fixed-SST and qflux models
      qsfix = .true.
C****
      RETURN
      END SUBROUTINE init_OCEAN

      SUBROUTINE daily_OCEAN(end_of_day)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi,edpery
      USE MODEL_COM, only : im,jm,kocean,focean,jday
      USE DAGCOM, only : aij,ij_toc2,ij_tgo2
      USE FLUXES, only : gtemp
      USE STATIC_OCEAN, only : tocean,ostruc,oclim,z1O,
     *     sinang,sn2ang,sn3ang,sn4ang,cosang,cs2ang,cs3ang,cs4ang
      IMPLICIT NONE
      INTEGER I,J
      LOGICAL, INTENT(IN) :: end_of_day
      REAL*8 ANGLE

C**** update ocean related climatologies
      CALL OCLIM(end_of_day)

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
      IF (KOCEAN.EQ.1.and.end_of_day) THEN
        DO J=1,JM
          DO I=1,IM
            AIJ(I,J,IJ_TOC2)=AIJ(I,J,IJ_TOC2)+TOCEAN(2,I,J)
            AIJ(I,J,IJ_TGO2)=AIJ(I,J,IJ_TGO2)+TOCEAN(3,I,J)
          END DO
        END DO
C**** DO DEEP DIFFUSION IF REQUIRED
        CALL ODIFS
C**** RESTRUCTURE THE OCEAN LAYERS
        CALL OSTRUC(.TRUE.)
      END IF

C**** set gtemp array for ocean temperature
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) GTEMP(1:2,1,I,J) = TOCEAN(1:2,I,J)
      END DO
      END DO
C****
      RETURN
      END SUBROUTINE daily_OCEAN

      SUBROUTINE PRECIP_OC
!@sum  PRECIP_OC driver for applying precipitation to ocean fraction
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : rhow,shw
      USE MODEL_COM, only : im,jm,focean,kocean,itocean,itoice
      USE GEOM, only : imaxj,dxyp
      USE DAGCOM, only : aj,j_implm,j_implh,oa,areg,jreg
      USE FLUXES, only : runpsi,prec,eprec,gtemp,mlhc
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      USE STATIC_OCEAN, only : tocean,z1o
      IMPLICIT NONE
      REAL*8 TGW,PRCP,WTRO,ENRGP,ERUN4,ENRGO,POCEAN,POICE,SNOW
     *     ,SMSI,ENRGW,WTRW0,WTRW,RUN0,RUN4,ROICE
      INTEGER I,J,JR

      DO J=1,JM
      DO I=1,IMAXJ(J)
        ROICE=RSI(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        JR=JREG(I,J)
        IF (FOCEAN(I,J).gt.0) THEN

          PRCP=PREC(I,J)
          ENRGP=EPREC(I,J)
          RUN0=RUNPSI(I,J)
          OA(I,J,4)=OA(I,J,4)+ENRGP

          IF (KOCEAN .EQ. 1) THEN
            TGW=TOCEAN(1,I,J)
            WTRO=Z1O(I,J)*RHOW
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
            AJ(J,J_IMPLM,ITOCEAN)=AJ(J,J_IMPLM,ITOCEAN)+RUN4 *POCEAN
            AJ(J,J_IMPLH,ITOCEAN)=AJ(J,J_IMPLH,ITOCEAN)+ERUN4*POCEAN
            AJ(J,J_IMPLM,ITOICE) =AJ(J,J_IMPLM,ITOICE) +RUN4 *POICE
            AJ(J,J_IMPLH,ITOICE) =AJ(J,J_IMPLH,ITOICE) +ERUN4*POICE
            AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+RUN4 *FOCEAN(I,J)*DXYP(J)
            AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+ERUN4*FOCEAN(I,J)*DXYP(J)
            MLHC(I,J)=WTRW*SHW  ! needed for underice fluxes
          END IF
          GTEMP(1,1,I,J)=TOCEAN(1,I,J)
        END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE PRECIP_OC

      SUBROUTINE OCEANS
!@sum  OCEANS driver for applying surface fluxes to ocean fraction
!@auth Original Development Team
!@ver  1.0
!@calls OCEAN:OSOURC
      USE CONSTANT, only : rhow,shw
      USE MODEL_COM, only : im,jm,focean,kocean,jday,dtsrc,itocean
     *     ,itoice
      USE GEOM, only : imaxj,dxyp
      USE DAGCOM, only : aj,areg,jreg,j_implm,j_implh,
     *     j_oht,j_imelt,j_hmelt,j_smelt,oa
      USE FLUXES, only : runosi,erunosi,srunosi,e0,e1,evapor,dmsi,dhsi
     *     ,dssi,flowo,eflowo,gtemp,sss,melti,emelti
#ifdef TRACERS_WATER
     *     ,dtrsi
#endif
      USE SEAICE, only : ace1i,ssi0,tfrez
      USE SEAICE_COM, only : rsi,msi,snowi
#ifdef TRACERS_WATER
     *     ,trsi0
#endif
      USE STATIC_OCEAN, only : tocean,z1o,ota,otb,otc,osourc,
     *     sinang,sn2ang,sn3ang,sn4ang,cosang,cs2ang,cs3ang,cs4ang
      IMPLICIT NONE
C**** grid box variables
      REAL*8 POCEAN, POICE, DXYPJ, TFO
C**** prognostic variables
      REAL*8 TGW, WTRO, SMSI, ROICE
C**** fluxes
      REAL*8 EVAPO, EVAPI, FIDT, FODT, OTDT, RVRRUN, RVRERUN, RUN0, 
     *     SALT, SIMELT, ESIMELT
C**** output from OSOURC
      REAL*8 ERUN4I, ERUN4O, RUN4I, RUN4O, ENRGFO, ACEFO, ACEFI, ENRGFI

      INTEGER I,J,JR

      DO J=1,JM
      DXYPJ=DXYP(J)
      DO I=1,IMAXJ(J)
        JR=JREG(I,J)
        ROICE=RSI(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        IF (FOCEAN(I,J).gt.0) THEN
          TGW  =TOCEAN(1,I,J)
          EVAPO=EVAPOR(I,J,1)
          EVAPI=EVAPOR(I,J,2)   ! evapor/dew at the ice surface (kg/m^2)
          FODT =E0(I,J,1)
          SMSI =MSI(I,J)+ACE1I+SNOWI(I,J)
          TFO  =tfrez(sss(i,j))
C**** get ice-ocean fluxes from sea ice routine
          RUN0=RUNOSI(I,J)  ! includes ACE2M + basal term
          FIDT=ERUNOSI(I,J)
          SALT=SRUNOSI(I,J)
C**** get river runoff/simelt flux
          RVRRUN = FLOWO(I,J)/(FOCEAN(I,J)*DXYPJ)
          RVRERUN=EFLOWO(I,J)/(FOCEAN(I,J)*DXYPJ)
          SIMELT = MELTI(I,J)/(FOCEAN(I,J)*DXYPJ)
          ESIMELT=EMELTI(I,J)/(FOCEAN(I,J)*DXYPJ)
c         SSIMELT=SMELTI(I,J)/(FOCEAN(I,J)*DXYPJ)
          OA(I,J,4)=OA(I,J,4)+RVRERUN ! add rvr E to surf. energy budget

          AJ(J,J_IMELT,ITOICE)=AJ(J,J_IMELT,ITOICE)+RUN0   *POICE
          AJ(J,J_SMELT,ITOICE)=AJ(J,J_SMELT,ITOICE)+SALT   *POICE
          AREG(JR,J_IMELT)=AREG(JR,J_IMELT)+RUN0*POICE*DXYP(J)
          AREG(JR,J_SMELT)=AREG(JR,J_SMELT)+SALT*POICE*DXYP(J)

          IF (KOCEAN .EQ. 1) THEN
            WTRO=Z1O(I,J)*RHOW
            OTDT=DTSRC*(OTA(I,J,4)*SN4ANG+OTB(I,J,4)*CS4ANG
     *           +OTA(I,J,3)*SN3ANG+OTB(I,J,3)*CS3ANG
     *           +OTA(I,J,2)*SN2ANG+OTB(I,J,2)*CS2ANG
     *           +OTA(I,J,1)*SINANG+OTB(I,J,1)*COSANG+OTC(I,J))

C**** Calculate the amount of ice formation
            CALL OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,FODT,FIDT,RVRRUN
     *           ,RVRERUN,SIMELT,ESIMELT,EVAPO,EVAPI,TFO,RUN4O,ERUN4O
     *           ,RUN4I,ERUN4I,ENRGFO,ACEFO,ACEFI,ENRGFI)

C**** Resave prognostic variables
            TOCEAN(1,I,J)=TGW
C**** Open Ocean diagnostics
            AJ(J,J_OHT  ,ITOCEAN)=AJ(J,J_OHT  ,ITOCEAN)+OTDT  *POCEAN
            AJ(J,J_IMPLM,ITOCEAN)=AJ(J,J_IMPLM,ITOCEAN)+RUN4O *POCEAN
            AJ(J,J_IMPLH,ITOCEAN)=AJ(J,J_IMPLH,ITOCEAN)+ERUN4O*POCEAN
C**** Ice-covered ocean diagnostics
            AJ(J,J_OHT  ,ITOICE)=AJ(J,J_OHT  ,ITOICE)+OTDT  *POICE
            AJ(J,J_IMPLM,ITOICE)=AJ(J,J_IMPLM,ITOICE)+RUN4I *POICE
            AJ(J,J_IMPLH,ITOICE)=AJ(J,J_IMPLH,ITOICE)+ERUN4I*POICE
C**** regional diagnostics
            AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+
     *             (RUN4O *POCEAN+RUN4I *POICE)*DXYPJ
            AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+
     *             (ERUN4O*POCEAN+ERUN4I*POICE)*DXYPJ
          ELSE
            ACEFO=0 ; ACEFI=0. ; ENRGFO=0. ; ENRGFI=0.
          END IF

C**** Store mass and energy fluxes for formation of sea ice
          DMSI(1,I,J)=ACEFO
          DMSI(2,I,J)=ACEFI
          DHSI(1,I,J)=ENRGFO
          DHSI(2,I,J)=ENRGFI
          DSSI(1,I,J)=SSI0*ACEFO   ! assume constant mean salinity
          DSSI(2,I,J)=SSI0*ACEFI
#ifdef TRACERS_WATER
C**** assume const mean tracer conc over freshwater amount
          DTRSI(:,1,I,J)=TRSI0(:)*(1.-SSI0)*ACEFO
          DTRSI(:,2,I,J)=TRSI0(:)*(1.-SSI0)*ACEFI
#endif
C**** store surface temperatures
          GTEMP(1:2,1,I,J)=TOCEAN(1:2,I,J)

        END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE OCEANS

      SUBROUTINE DIAGCO (M)
!@sum  DIAGCO Keeps track of the ocean conservation properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : kocean
      USE DAGCOM, only : icon_OCE
      IMPLICIT NONE
!@var M index denoting from where DIAGCO is called
      INTEGER, INTENT(IN) :: M
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCO IS BEING CALLED
C****     (see DIAGCA)
      REAL*8, EXTERNAL :: conserv_OCE

C**** OCEAN POTENTIAL ENTHALPY
      IF (KOCEAN.gt.0) CALL conserv_DIAG(M,conserv_OCE,icon_OCE)
C****
      RETURN
      END SUBROUTINE DIAGCO


      SUBROUTINE io_oda(kunit,it,iaction,ioerr)
!@sum  io_oda reads/writes ocean/ice data for initializing deep ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,Itime,im,jm
      USE DAGCOM, only : ij_tgo2,aij
      USE SEAICE_COM, only : rsi,msi,hsi,ssi
      USE STATIC_OCEAN, only : tocean
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
      INTEGER I,J

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output
        WRITE (kunit,err=10) it,TOCEAN,RSI,MSI,HSI,SSI,
     *     ((AIJ(I,J,IJ_TGO2),I=1,IM),J=1,JM)
      CASE (IOREAD:)            ! input
        READ (kunit,err=10) it,TOCEAN,RSI,MSI,HSI,SSI,
     *     ((AIJ(I,J,IJ_TGO2),I=1,IM),J=1,JM)
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_oda

#ifdef TRACERS_WATER
      subroutine tracer_ic_ocean
!@sum tracer_ic_ocean initialise ocean tracer concentration
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,focean,itime
      USE GEOM, only : imaxj
      USE TRACER_COM, only : trw0,ntm,itime_tr0
      USE FLUXES, only : gtracer
      IMPLICIT NONE
      INTEGER i,j,n

      do n=1,ntm
        if (itime.eq.itime_tr0(n)) then
          do j=1,jm
          do i=1,imaxj(j)
            if (focean(i,j).gt.0) gtracer(n,1,i,j)=trw0(n)
          end do
          end do
        end if
      end do

      return
      end subroutine tracer_ic_ocean
#endif
