#include "rundeck_opts.h"

      MODULE STATIC_OCEAN
!@sum  STATIC_OCEAN contains the ocean subroutines common to all Q-flux
!@+    and fixed SST runs
!@auth Original Development Team
!@ver  1.0 (Q-flux ocean)
!@cont OSTRUC,OCLIM,init_OCEAN,daily_OCEAN,DIAGCO
!@+    PRECIP_OC,OCEANS
      USE CONSTANT, only : lhm,rhows,rhoi,shw,shi,by12,byshi
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP, only : GRID, GET,
     *                          DREAD_PARALLEL,
     *                          MREAD_PARALLEL,
     *                          READT_PARALLEL,
     *                          REWIND_PARALLEL,
     *                          BACKSPACE_PARALLEL
      USE MODEL_COM, only : im,jm,lm,focean,fland,fearth,flice
     *     ,Iyear1,Itime,jmon,jdate,jday,jyear,jmpery,JDendOfM,JDmidOfM
     *     ,ItimeI,kocean,itocean,itoice
      USE GEOM
      USE PBLCOM, only : npbl,uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs
     *     ,cq=>cqgs,ipbl,roughl
      USE SEAICE, only : xsi,ace1i,z1i,ac2oim,z2oim,ssi0,tfrez,fleadoc
     *     ,lmi
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
#ifdef TRACERS_WATER
     *     ,trsi,trsi0,ntm
#endif
      USE LANDICE_COM, only : snowli,tlandi
      USE FLUXES, only : gtemp,sss,fwsim,mlhc
      USE DIAG_COM, only : aij=>aij_loc, areg, jreg,ij_smfx,aj=>aj_loc
     *     ,j_implh, j_implm, j_imelt, j_hmelt, j_smelt, NREG, KAJ
      IMPLICIT NONE
      SAVE
      logical :: off_line=.false.
!@dbparam ocn_cycl determines whether prescribed ocean data repeat
!@+       after 1 year - 0:no 1:yes 2:seaice repeats,sst does not
      integer :: ocn_cycl = 1
!@var TOCEAN temperature of the ocean (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TOCEAN
!@var SSS0 default sea surface salinity (psu)
      REAL*8 :: SSS0=34.7d0

!@var OTA,OTB,OTC ocean heat transport coefficients
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: OTA,OTB
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OTC
!@var SINANG,COSANG Fourier coefficients for heat transports
      REAL*8 :: SINANG,SN2ANG,SN3ANG,SN4ANG,
     *          COSANG,CS2ANG,CS3ANG,CS4ANG

!@var Z1O ocean mixed layer depth
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Z1O
!@var Z1OOLD previous ocean mixed layer depth
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Z1OOLD
!@var Z12O annual maximum ocean mixed layer depth
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: Z12O

      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::DM,AOST,EOST1,EOST0,BOST, 
     *     COST,ARSI,ERSI1,ERSI0,BRSI,CRSI
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  KRSI
C      COMMON/OOBS/DM,AOST,EOST1,EOST0,BOST,COST,ARSI,ERSI1,ERSI0,BRSI
C     *     ,CRSI,KRSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: XZO,XZN

!@var iu_OSST,iu_SICE,iu_OCNML unit numbers for climatologies
      INTEGER iu_OSST,iu_SICE,iu_OCNML

!@dbparam qflux_fix an energy leak (default=0; no fix for continuity)
      INTEGER :: qflux_fix=0

!@dbparam qfluxX multiplying factor for qfluxes 
      REAL*8 :: qfluxX=1.

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
      INTEGER :: J_0, J_1
C****
C**** FLAND     LAND COVERAGE (1)
C****
C**** TOCEAN 1  OCEAN TEMPERATURE OF FIRST LAYER (C)
C****        2  MEAN OCEAN TEMPERATURE OF SECOND LAYER (C)
C****        3  OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)
C****     RSI   RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** FWSIM     Fresh water sea ice amount (KG/M**2)
C****
C**** RESTRUCTURE OCEAN LAYERS
C****
      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF (FLAND(I,J).GE.1.) CYCLE
          IF (Z1OOLD(I,J).GE.Z12O(I,J)) GO TO 140
          IF (Z1O(I,J).EQ.Z1OOLD(I,J)) CYCLE
          WTR1O=RHOWS*Z1O(I,J)-FWSIM(I,J)
          DWTRO=RHOWS*(Z1O(I,J)-Z1OOLD(I,J))
          WTR2O=RHOWS*(Z12O(I,J)-Z1O(I,J))
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
      USE DOMAIN_DECOMP, only : HERE
      USE DOMAIN_DECOMP, ONLY : GLOBALSUM,AM_I_ROOT
      IMPLICIT NONE

C now allocated from ALLOC_STATIC OCEAN      REAL*8, SAVE :: XZO(IM,JM),XZN(IM,JM)
      LOGICAL, INTENT(IN) :: end_of_day

      INTEGER n,J,I,LSTMON,K,m,m1,JR
      REAL*8 RSICSQ,ZIMIN,ZIMAX,Z1OMIN,RSINEW,TIME,FRAC,MSINEW,OPNOCN
     *     ,TFO
!@var JDLAST julian day that OCLIM was last called
      INTEGER, SAVE :: JDLAST=0
!@var IMON0 current month for SST climatology reading
      INTEGER, SAVE :: IMON0 = 0
!@var IMON current month for ocean mixed layer climatology reading
      INTEGER, SAVE :: IMON = 0
!@var TEMP_LOCAL stores AOST+EOST1 or ARSI+ERST1 to avoid the use of common block OOBS in MODULE STATIC_OCEAN
      REAL*8 :: TEMP_LOCAL(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,2)
      REAL*8 :: AREG_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO,NREG,KAJ)
      REAL*8 :: gsum

      INTEGER :: J_0,J_1
      LOGICAL :: HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      IF (KOCEAN.EQ.1) GO TO 500
      if(.not.(end_of_day.or.itime.eq.itimei.or.off_line)) return
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
          CALL READT_PARALLEL
     *           (grid,iu_SICE,NAMEUNIT(iu_SICE),IM*JM,ERSI0,LSTMON)
          if (ocn_cycl.eq.1) then
            CALL READT_PARALLEL
     *           (grid,iu_OSST,NAMEUNIT(iu_OSST),IM*JM,EOST0,LSTMON)
          else ! if (ocn_cycl.eq.2) then
            LSTMON=JMON-1+(JYEAR-IYEAR1)*JMperY
  290       read (iu_OSST) M
            if (m.lt.lstmon) go to 290
            CALL BACKSPACE_PARALLEL( iu_OSST )
            CALL MREAD_PARALLEL
     *           (GRID,iu_OSST,NAMEUNIT(iu_OSST),m,IM*JM,EOST0)
            WRITE(6,*) 'Read End-of-month ocean data from ',JMON-1,M
            IF(M.NE.LSTMON)
     &         call stop_model('Read error: ocean data',255)
          end if
        else   !  if (ocn_cycl.eq.0) then
          LSTMON=JMON-1+(JYEAR-IYEAR1)*JMperY
  300     read (iu_OSST) M
          if (m.lt.lstmon) go to 300
          CALL BACKSPACE_PARALLEL( iu_OSST )
  310     read (iu_SICE) m
          if (m.lt.lstmon) go to 310
          CALL BACKSPACE_PARALLEL( iu_SICE )
          CALL MREAD_PARALLEL
     *           (GRID,iu_OSST,NAMEUNIT(iu_OSST), m,IM*JM,EOST0)
          CALL MREAD_PARALLEL
     *           (GRID,iu_SICE,NAMEUNIT(iu_SICE),m1,IM*JM,ERSI0)
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
          if (ocn_cycl.eq.1) CALL REWIND_PARALLEL( iu_OSST )
          CALL REWIND_PARALLEL( iu_SICE )
          read (iu_SICE)  
        end if
        CALL READT_PARALLEL
     *           (grid,iu_SICE,NAMEUNIT(iu_SICE),0,TEMP_LOCAL,1)
        ARSI  = TEMP_LOCAL(:,:,1)
        ERSI1 = TEMP_LOCAL(:,:,2)
        if (ocn_cycl.eq.1) then
          CALL READT_PARALLEL
     *           (grid,iu_OSST,NAMEUNIT(iu_OSST),0,TEMP_LOCAL,1)
          AOST  = TEMP_LOCAL(:,:,1)
          EOST1 = TEMP_LOCAL(:,:,2)
        else  ! if (ocn_cycl.eq.2) then
          CALL MREAD_PARALLEL
     *           (GRID,iu_OSST,NAMEUNIT(iu_OSST),M,0,TEMP_LOCAL)
          AOST  = TEMP_LOCAL(:,:,1)
          EOST1 = TEMP_LOCAL(:,:,2)
          WRITE(6,*) 'Read in ocean data for month',JMON,M
          IF(JMON.NE.MOD(M-1,12)+1)
     &       call stop_model('Error: Ocean data',255)
        end if
      else   !  if (ocn_cycl.eq.0) then
        CALL MREAD_PARALLEL
     *           (GRID,iu_OSST,NAMEUNIT(iu_OSST),M,0,TEMP_LOCAL)
        AOST  = TEMP_LOCAL(:,:,1)
        EOST1 = TEMP_LOCAL(:,:,2)
        CALL MREAD_PARALLEL
     *           (GRID,iu_SICE,NAMEUNIT(iu_SICE),M1,0,TEMP_LOCAL)
        ARSI  = TEMP_LOCAL(:,:,1)
        ERSI1 = TEMP_LOCAL(:,:,2)
        WRITE(6,*) 'Read in ocean data for month',JMON,M,M1
        IF(M.NE.M1.OR.JMON.NE.MOD(M-1,12)+1)
     &       call stop_model('Error: Ocean data',255)
      end if
C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO J=J_0,J_1
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
      DO J=J_0,J_1
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
            IF(MSI(I,J).eq.0) MSI(I,J)=MSINEW  ! does this happen?
            IF (end_of_day.and..not.off_line) THEN
              AIJ(I,J,IJ_SMFX)=AIJ(I,J,IJ_SMFX)+
     *           (SNOWI(I,J)+ACE1I)*(RSINEW-RSI(I,J))+
     *           RSINEW*MSINEW-RSI(I,J)*MSI(I,J)
              AJ(J,J_IMPLM,ITOICE)=AJ(J,J_IMPLM,ITOICE)-FOCEAN(I,J)
     *             *RSI(I,J)*(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))
     *               /MSI(I,J))
              AJ(J,J_IMPLH,ITOICE)=AJ(J,J_IMPLH,ITOICE)-FOCEAN(I,J)
     *             *SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)*RSI(I,J)
              AJ(J,J_IMPLM,ITOCEAN)=AJ(J,J_IMPLM,ITOCEAN)-FOCEAN(I,J)
     *             *(RSINEW-RSI(I,J))*(MSINEW+ACE1I+SNOWI(I,J)
     *             -SUM(SSI(1:2,I,J))
     *             -SUM(SSI(3:4,I,J))*(MSINEW/MSI(I,J)))
              AJ(J,J_IMPLH,ITOCEAN)=AJ(J,J_IMPLH,ITOCEAN)-FOCEAN(I,J)
     *             *(RSINEW-RSI(I,J))*(SUM(HSI(1:2,I,J))+
     *             SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)))
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
              HSI(1:2,I,J)=(SHI*TFO-LHM*(1.-SSI0))*XSI(1:2)*ACE1I
              HSI(3:4,I,J)=(SHI*TFO-LHM*(1.-SSI0))*XSI(3:4)*AC2OIM
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
            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-SUM(SSI(1:LMI
     *           ,I,J)))
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE (for prescribed data only)
      IF(HAVE_NORTH_POLE) THEN
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
            FWSIM(I,JM)=FWSIM(1,JM)
          END DO
        END IF
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
          CALL READT_PARALLEL(grid,iu_OCNML,NAMEUNIT(iu_OCNML),0,XZO,12)
          CALL REWIND_PARALLEL( iu_OCNML )
        ELSE            ! JDAY is in Jan 16 to Dec 16, get first month
  520     IMON=IMON+1
          IF (JDAY.GT.JDmidOFM(IMON) .and. IMON.LE.12) GO TO 520
          CALL READT_PARALLEL
     *           (grid,iu_OCNML,NAMEUNIT(iu_OCNML),0,XZO,IMON-1)
          IF (IMON.EQ.13)  CALL REWIND_PARALLEL( iu_OCNML )
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
        IF (IMON.EQ.13) CALL REWIND_PARALLEL( iu_OCNML )
      END IF
      CALL READT_PARALLEL(grid,iu_OCNML,NAMEUNIT(iu_OCNML),0,XZN,1)
 530  JDLAST=JDAY

C**** Interpolate the mixed layer depth z1o to the current day and
C**** limit it to the annual maxmimal mixed layer depth z12o
      FRAC = REAL(JDmidOFM(IMON)-JDAY,KIND=8)/
     a           (JDmidOFM(IMON)-JDmidOFM(IMON-1))
      areg_part = 0

      DO J=J_0,J_1
      DO I=1,IM
      Z1O(I,J)=min( z12o(i,j) , FRAC*XZO(I,J)+(1.-FRAC)*XZN(I,J) )

      IF (RSI(I,J)*FOCEAN(I,J).GT.0.and.off_line) THEN
        Z1OMIN=1.+FWSIM(I,J)/(RHOWS*RSI(I,J))
        IF (Z1OMIN.GT.Z1O(I,J)) THEN
C**** MIXED LAYER DEPTH IS INCREASED TO OCEAN ICE DEPTH + 1 METER
          WRITE(6,602) ITime,I,J,JMON,Z1O(I,J),Z1OMIN,z12o(i,j)
 602      FORMAT (' INCREASE OF MIXED LAYER DEPTH ',I10,3I4,3F10.3)
          Z1O(I,J)=MIN(Z1OMIN, z12o(i,j))
          IF (Z1OMIN.GT.Z12O(I,J)) THEN
C****       ICE DEPTH+1>MAX MIXED LAYER DEPTH :
C****       lose the excess mass to the deep ocean
C**** Calculate freshwater mass to be removed, and then any energy/salt
            MSINEW=MSI(I,J)*(1.-RHOWS*(Z1OMIN-Z12O(I,J))/(FWSIM(I,J)
     *           -RSI(I,J)*(ACE1I+SNOWI(I,J)-SUM(SSI(1:2,I,J)))))
            HSI(3:4,I,J) = HSI(3:4,I,J)*(MSINEW/MSI(I,J))
            SSI(3:4,I,J) = SSI(3:4,I,J)*(MSINEW/MSI(I,J))
#ifdef TRACERS_WATER
            TRSI(:,3:4,I,J) = TRSI(:,3:4,I,J)*(MSINEW/MSI(I,J))
#endif
            AJ(J,J_IMELT,ITOICE)=AJ(J,J_IMELT,ITOICE)-FOCEAN(I,J)
     *           *RSI(I,J)*(MSINEW-MSI(I,J))
            AJ(J,J_HMELT,ITOICE)=AJ(J,J_HMELT,ITOICE)-FOCEAN(I,J)
     *           *RSI(I,J)*SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)
            AJ(J,J_SMELT,ITOICE)=AJ(J,J_SMELT,ITOICE)-FOCEAN(I,J)
     *           *RSI(I,J)*SUM(SSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)
            AJ(J,J_IMPLM,ITOICE)=AJ(J,J_IMPLM,ITOICE)-FOCEAN(I,J)
     *           *RSI(I,J)*(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))
     *           /MSI(I,J))
            AJ(J,J_IMPLH,ITOICE)=AJ(J,J_IMPLH,ITOICE)-FOCEAN(I,J)
     *           *RSI(I,J)*SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)
            JR=JREG(I,J)
            AREG_part(J,JR,J_IMPLM)=AREG_part(J,JR,J_IMPLM)-
     *           FOCEAN(I,J)*RSI(I,J)
     *           *(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))
     *               /MSI(I,J))*DXYP(J)
            AREG_part(J,JR,J_IMPLH)=AREG_part(J,JR,J_IMPLH)-
     *           FOCEAN(I,J)*RSI(I,J)
     *           *SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)*DXYP(J)
            MSI(I,J)=MSINEW
            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-SUM(SSI(1:LMI
     *           ,I,J)))
          END IF
        END IF
      END IF
      END DO
      END DO
      DO JR = 1, NREG
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLM), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+gsum
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLH), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+gsum
      END DO

      RETURN
      END SUBROUTINE OCLIM

      SUBROUTINE OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,FODT,FIDT,RVRRUN
     *     ,RVRERUN,EVAPO,EVAPI,TFW,WTRW,RUN4O,ERUN4O
     *     ,RUN4I,ERUN4I,ENRGFO,ACEFO,ACEFI,ENRGFI)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@+    ACEFO/I are the freshwater ice amounts,
!@+    ENRGFO/I is total energy (including a salt component)
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@var TFW freezing temperature for water underlying ice (C)
      REAL*8, INTENT(IN) :: TFW

!@var TGW mixed layer temp.(C)
      REAL*8, INTENT(INOUT) :: TGW
      REAL*8, INTENT(IN) :: ROICE, SMSI, WTRO, EVAPO, EVAPI, RUN0
      REAL*8, INTENT(IN) :: FODT, FIDT, OTDT, RVRRUN, RVRERUN
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACEFI, RUN4O, RUN4I
     *     , ERUN4O, ERUN4I, WTRW
      REAL*8 EIW0, ENRGIW, WTRI1, EFIW, ENRGO0, EOFRZ, ENRGO
      REAL*8 WTRW0, ENRGW0, ENRGW, WTRI0, SMSI0
!@var emin min energy deficit required before ice forms (J/m^2)
      REAL*8, PARAMETER :: emin=-1d-10
C**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACEFI=0. ; ENRGFI=0. ; WTRW=WTRO

C**** Calculate previous ice mass (before fluxes applied)
      SMSI0=SMSI+RUN0+EVAPI

C**** Calculate extra mass flux to ocean, balanced by deep removal
      RUN4O=-EVAPO+RVRRUN  ! open ocean
      RUN4I=-EVAPI+RVRRUN  ! under ice
! force energy consveration
      ERUN4O=0. ! RUN4O*TGW*SHW ! corresponding heat flux at bottom (open)
      ERUN4I=0. ! RUN4I*TGW*SHW !                              (under ice)

C**** Calculate heat fluxes to ocean
      ENRGO  = FODT+OTDT+RVRERUN-ERUN4O ! in open water
      ENRGIW = FIDT+OTDT+RVRERUN-ERUN4I ! under ice

C**** Calculate energy in mixed layer (open ocean)
      ENRGO0=WTRO*TGW*SHW
      EOFRZ=WTRO*TFW*SHW
      IF (ENRGO0+ENRGO.GE.EOFRZ+emin) THEN
C**** FLUXES RECOMPUTE TGW WHICH IS ABOVE FREEZING POINT FOR OCEAN
        IF (ROICE.LE.0.) THEN
          TGW=TGW+ENRGO/(WTRO*SHW)
          RETURN
        END IF
      ELSE
C**** FLUXES COOL TGO TO FREEZING POINT FOR OCEAN AND FORM SOME ICE
        ACEFO=(ENRGO0+ENRGO-EOFRZ)/(TFW*(SHI-SHW)-LHM)
        ENRGFO=ACEFO*(TFW*SHI/(1.-SSI0)-LHM)
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
      IF (EIW0+ENRGIW .LE. EFIW+emin) THEN ! freezing case
C**** FLUXES WOULD COOL TGW TO FREEZING POINT AND FREEZE SOME MORE ICE
        ACEFI = (EIW0+ENRGIW-EFIW)/(TFW*(SHI-SHW)-LHM)
        ENRGFI = ACEFI*(TFW*SHI/(1.-SSI0)-LHM) ! energy of frozen ice
      END IF
C**** COMBINE OPEN OCEAN AND SEA ICE FRACTIONS TO FORM NEW VARIABLES
      WTRW = WTRO  -(1.-ROICE)* ACEFO        -ROICE*(SMSI+ACEFI)
      ENRGW= ENRGW0+(1.-ROICE)*(ENRGO-ENRGFO)+ROICE*(ENRGIW-ENRGFI)
      TGW  = ENRGW/(WTRW*SHW) ! new ocean temperature

      RETURN
      END SUBROUTINE OSOURC

      END MODULE STATIC_OCEAN

      SUBROUTINE ALLOC_STATIC_OCEAN(grid)
      USE MODEL_COM, only : im
      USE STATIC_OCEAN, only  : TOCEAN,OTA,OTB,OTC,Z1O,Z1OOLD,Z12O,
     &                          DM,AOST,EOST1,EOST0,BOST,
     &                          COST,ARSI,ERSI1,ERSI0,BRSI,CRSI,XZO,XZN
      USE STATIC_OCEAN, only  : KRSI
      USE DOMAIN_DECOMP, only : DIST_GRID,GET
      IMPLICIT NONE
      INTEGER :: J_0H,J_1H,IER
      TYPE (DIST_GRID), INTENT(IN) :: grid

      CALL GET(GRID,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      ALLOCATE(TOCEAN(3,IM,J_0H:J_1H),
     &    STAT=IER)
      ALLOCATE(OTA(IM,J_0H:J_1H,4),
     &         OTB(IM,J_0H:J_1H,4),
     &    STAT=IER)
      ALLOCATE(KRSI(IM,J_0H:J_1H),
     &    STAT=IER)
      ALLOCATE(OTC(IM,J_0H:J_1H),
     &         Z1O(IM,J_0H:J_1H),
     &         Z1OOLD(IM,J_0H:J_1H),
     &         Z12O(IM,J_0H:J_1H),
     &         DM(IM,J_0H:J_1H),
     &         AOST(IM,J_0H:J_1H),
     &         EOST1(IM,J_0H:J_1H),
     &         EOST0(IM,J_0H:J_1H),
     &         BOST(IM,J_0H:J_1H),
     &         COST(IM,J_0H:J_1H),
     &         ARSI(IM,J_0H:J_1H),
     &         ERSI1(IM,J_0H:J_1H),
     &         ERSI0(IM,J_0H:J_1H),
     &         BRSI(IM,J_0H:J_1H),
     &         CRSI(IM,J_0H:J_1H),
     &         XZO(IM,J_0H:J_1H),
     &         XZN(IM,J_0H:J_1H),
     &    STAT=IER)
      END SUBROUTINE ALLOC_STATIC_OCEAN


      SUBROUTINE init_OCEAN(iniOCEAN,istart)
!@sum init_OCEAN initiallises ocean variables
!@auth Original Development Team
!@ver  1.0
      USE FILEMANAGER
      USE PARAM
      USE DOMAIN_DECOMP, only : GRID, GET,
     *                          DREAD_PARALLEL,
     *                          MREAD_PARALLEL,
     *                          READT_PARALLEL,
     *                          REWIND_PARALLEL,
     *                          BACKSPACE_PARALLEL
      USE MODEL_COM, only : im,jm,fland,flice,kocean,focean
     *     ,fearth,iyear1,ioreadnt,jmpery
#ifdef TRACERS_WATER
      USE TRACER_COM, only : trw0
      USE FLUXES, only : gtracer
#endif
      USE FLUXES, only : gtemp,sss,uosurf,vosurf,uisurf,visurf,ogeoza
      USE SEAICE, only : qsfix, osurf_tilt
      USE SEAICE_COM, only : snowi
      USE STATIC_OCEAN, only : ota,otb,otc,z12o,dm,iu_osst,iu_sice
     *     ,iu_ocnml,tocean,ocn_cycl,sss0,qflux_fix,qfluxX
      USE DIAG_COM, only : npts,icon_OCE,conpt0
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      LOGICAL, INTENT(IN) :: iniOCEAN  ! true if starting from ic.
      CHARACTER CONPT(NPTS)*10
!@var iu_OHT unit number for reading in ocean heat transports & z12o_max
      INTEGER :: iu_OHT,iu_GIC
      INTEGER :: I,J,m,ISTART,ioerr
!@var z12o_max maximal mixed layer depth (m) for qflux model
      real*8 :: z12o_max
      real*8 z1ox(im,grid%j_strt_halo:grid%j_stop_halo)
      integer :: j_0,j_1

      call get(grid,j_strt=j_0,j_stop=j_1)

      if (kocean.ge.1) then
C****   set conservation diagnostic for ocean heat
        CONPT=CONPT0
        QCON=(/ F, F, F, T, F, T, F, T, T, F, F/)
        CALL SET_CON(QCON,CONPT,"OCN HEAT","(10^6 J/M**2)   ",
     *       "(W/M^2)         ",1d-6,1d0,icon_OCE)
      end if

      if (istart.le.0) return

      call sync_param( "qflux_fix",qflux_fix)

      call sync_param( "qfluxX"   ,qfluxX)

C**** if starting from AIC/GIC files need additional read for ocean
      if (istart.le.2) then
        call openunit("GIC",iu_GIC,.true.,.true.)
        ioerr=-1
        call io_ocean (iu_GIC,ioreadnt,ioerr)
        if (ioerr.eq.1) then
          WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
          call stop_model("INPUT: GIC READ IN ERROR",255)
        end if
        call closeunit (iu_GIC)
      end if

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
        CALL READT_PARALLEL(grid,iu_SICE,NAMEUNIT(iu_SICE),0,DM,1)

      else !  IF (KOCEAN.eq.1) THEN
C****   DATA FOR QFLUX MIXED LAYER OCEAN RUNS
C****   read in ocean heat transport coefficients
        call openunit("OHT",iu_OHT,.true.,.true.)
        READ (iu_OHT) OTA,OTB,OTC,z12o_max
        WRITE(6,*) "Read ocean heat transports from OHT"
        call closeunit (iu_OHT)

C****   Set up unit number of observed mixed layer depth data
        call openunit("OCNML",iu_OCNML,.true.,.true.)

C****   find and limit ocean ann max mix layer depths
        z12o = 0.
        do m=1,jmpery
          CALL READT_PARALLEL(grid,iu_OCNML,NAMEUNIT(iu_OCNML),0,z1ox,1)
          do j=j_0,j_1
          do i=1,im
ccc         z12o(i,j)=min( z12o_max , max(z12o(i,j),z1ox(i,j)) )
ccc     the above line could substitute for next 3 lines w/o any change
            if (focean(i,j).gt.0. .and. z1ox(i,j).gt.z12o_max)
     *          z1ox(i,j)=z12o_max
            if (z1ox(i,j).gt.z12o(i,j)) z12o(i,j)=z1ox(i,j)
          end do
          end do
        end do
        CALL REWIND_PARALLEL( iu_OCNML )
        write(6,*) 'Mixed Layer Depths limited to',z12o_max

C****   initialise deep ocean arrays if required
        call init_ODEEP(iniOCEAN)

      END IF
C**** Set fluxed arrays for oceans
      DO J=J_0,J_1
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
C**** For the time being assume zero surface velocities for drag calc
        uosurf(i,j)=0. ; vosurf(i,j)=0.
        uisurf(i,j)=0. ; visurf(i,j)=0.
C**** Also zero out surface height variations
        ogeoza(i,j)=0.
      END DO
      END DO
C**** keep salinity in sea ice constant for fixed-SST and qflux models
      qsfix = .true.
C**** Make sure to use geostrophy for ocean tilt term in ice dynamics
C**** (if required). Since ocean currents are zero, this implies no sea
C**** surface tilt term.
      osurf_tilt = 0
C****
      RETURN
      END SUBROUTINE init_OCEAN

      SUBROUTINE daily_OCEAN(end_of_day)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi,edpery,shw,rhows
      USE MODEL_COM, only : im,jm,kocean,focean,jday
      USE DIAG_COM, only : aij=>aij_loc,ij_toc2,ij_tgo2
      USE FLUXES, only : gtemp,mlhc,fwsim
      USE STATIC_OCEAN, only : tocean,ostruc,oclim,z1o,
     *     sinang,sn2ang,sn3ang,sn4ang,cosang,cs2ang,cs3ang,cs4ang
      USE DOMAIN_DECOMP, only : GRID,GET
      IMPLICIT NONE
      INTEGER I,J
      LOGICAL, INTENT(IN) :: end_of_day
      REAL*8 ANGLE
      INTEGER :: J_0,J_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

C**** update ocean related climatologies
      CALL OCLIM(end_of_day)

C**** Set fourier coefficients for heat transport calculations
      IF (KOCEAN.ge.1) THEN
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
      IF (KOCEAN.ge.1.and.end_of_day) THEN
        DO J=J_0,J_1
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
      DO J=J_0,J_1
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          GTEMP(1:2,1,I,J) = TOCEAN(1:2,I,J)
          MLHC(I,J) = SHW*(Z1O(I,J)*RHOWS-FWSIM(I,J))
        END IF
      END DO
      END DO
C****
      RETURN
      END SUBROUTINE daily_OCEAN

      SUBROUTINE PRECIP_OC
!@sum  PRECIP_OC driver for applying precipitation to ocean fraction
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : rhows,shw
      USE MODEL_COM, only : im,jm,focean,kocean,itocean,itoice
      USE GEOM, only : imaxj,dxyp
      USE DIAG_COM, only : aj=>aj_loc,j_implm,j_implh,oa,areg,jreg
      USE FLUXES, only : runpsi,srunpsi,prec,eprec,gtemp,mlhc,melti
     *     ,emelti,smelti,fwsim
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      USE STATIC_OCEAN, only : tocean,z1o
      USE DOMAIN_DECOMP, only : GRID,GET,AM_I_ROOT,GLOBALSUM
      USE DIAG_COM, only : NREG, KAJ
      IMPLICIT NONE
      REAL*8 TGW,PRCP,WTRO,ENRGP,ERUN4,POCEAN,POICE,SNOW
     *     ,SMSI0,ENRGW,WTRW0,WTRW,RUN0,RUN4,ROICE,SIMELT,ESIMELT
      INTEGER I,J,JR
      INTEGER :: J_0,J_1
      REAL*8 :: AREG_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO,NREG,KAJ)
      REAL*8 :: gsum

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

      areg_part = 0

      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        ROICE=RSI(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        JR=JREG(I,J)
        IF (FOCEAN(I,J).gt.0) THEN
          PRCP=PREC(I,J)
          ENRGP=EPREC(I,J)
          RUN0=RUNPSI(I,J)-SRUNPSI(I,J) ! fresh water runoff
          SIMELT =(MELTI(I,J)-SMELTI(I,J))/(FOCEAN(I,J)*DXYP(J))
          ESIMELT=EMELTI(I,J)/(FOCEAN(I,J)*DXYP(J))
          OA(I,J,4)=OA(I,J,4)+ENRGP

          IF (KOCEAN .ge. 1) THEN
            TGW=TOCEAN(1,I,J)
            WTRO=Z1O(I,J)*RHOWS
            SNOW=SNOWI(I,J)
            SMSI0=FWSIM(I,J)+ROICE*(RUN0-PRCP) ! initial ice

C**** Calculate effect of lateral melt of sea ice
            IF (SIMELT.gt.0) THEN
              TGW=TGW + (ESIMELT-SIMELT*TGW*SHW)/
     *             (SHW*(WTRO-SMSI0))
            END IF

C**** Additional mass (precip) is balanced by deep removal
            RUN4=PRCP
            ERUN4=0. ! RUN4*TGW*SHW ! force energy consveration

            IF (POICE.LE.0.) THEN
              ENRGW=TGW*WTRO*SHW + ENRGP - ERUN4
              WTRW =WTRO
            ELSE
              WTRW0=WTRO -SMSI0
              ENRGW=WTRW0*TGW*SHW + (1.-ROICE)*ENRGP - ERUN4
              WTRW =WTRW0+ROICE*(RUN0-RUN4)
            END IF
            TGW=ENRGW/(WTRW*SHW)
            TOCEAN(1,I,J)=TGW
            AJ(J,J_IMPLM,ITOCEAN)=AJ(J,J_IMPLM,ITOCEAN)+RUN4 *POCEAN
            AJ(J,J_IMPLH,ITOCEAN)=AJ(J,J_IMPLH,ITOCEAN)+ERUN4*POCEAN
            AJ(J,J_IMPLM,ITOICE) =AJ(J,J_IMPLM,ITOICE) +RUN4 *POICE
            AJ(J,J_IMPLH,ITOICE) =AJ(J,J_IMPLH,ITOICE) +ERUN4*POICE
            AREG_part(J,JR,J_IMPLM)=AREG_part(J,JR,J_IMPLM)+
     &           RUN4 *FOCEAN(I,J)*DXYP(J)
            AREG_part(J,JR,J_IMPLH)=AREG_part(J,JR,J_IMPLH)+
     &           ERUN4*FOCEAN(I,J)*DXYP(J)
            MLHC(I,J)=WTRW*SHW  ! needed for underice fluxes
          END IF
          GTEMP(1,1,I,J)=TOCEAN(1,I,J)
        END IF
      END DO
      END DO
      DO JR = 1, NREG
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLM), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+gsum
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLH), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+gsum
      END DO

      RETURN
C****
      END SUBROUTINE PRECIP_OC

      SUBROUTINE OCEANS
!@sum  OCEANS driver for applying surface fluxes to ocean fraction
!@auth Original Development Team
!@ver  1.0
!@calls OCEAN:OSOURC
      USE CONSTANT, only : rhows,shw
      USE MODEL_COM, only : im,jm,focean,kocean,jday,dtsrc,itocean
     *     ,itoice
      USE GEOM, only : imaxj,dxyp
      USE DIAG_COM, only : aj=>aj_loc,areg,jreg,j_implm,j_implh,j_oht,oa
      USE FLUXES, only : runosi,erunosi,srunosi,e0,e1,evapor,dmsi,dhsi
     *     ,dssi,flowo,eflowo,gtemp,sss,fwsim,mlhc,gmelt,egmelt
#ifdef TRACERS_WATER
     *     ,dtrsi
#endif
      USE SEAICE, only : ace1i,ssi0,tfrez
      USE SEAICE_COM, only : rsi,msi,snowi
#ifdef TRACERS_WATER
     *     ,trsi0
#endif
      USE STATIC_OCEAN, only : tocean,z1o,ota,otb,otc,osourc,qfluxX,
     *     sinang,sn2ang,sn3ang,sn4ang,cosang,cs2ang,cs3ang,cs4ang
      USE DOMAIN_DECOMP, only : GRID,GET,AM_I_ROOT,GLOBALSUM
      USE DIAG_COM, only : NREG, KAJ
      IMPLICIT NONE
C**** grid box variables
      REAL*8 POCEAN, POICE, DXYPJ, TFO
C**** prognostic variables
      REAL*8 TGW, WTRO, SMSI, ROICE
C**** fluxes
      REAL*8 EVAPO, EVAPI, FIDT, FODT, OTDT, RVRRUN, RVRERUN, RUN0
C**** output from OSOURC
      REAL*8 ERUN4I, ERUN4O, RUN4I, RUN4O, ENRGFO, ACEFO, ACEFI, ENRGFI,
     *     WTRW

      INTEGER I,J,JR
      INTEGER :: J_0,J_1
      REAL*8 :: AREG_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO,NREG,KAJ)
      REAL*8 :: gsum

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

      areg_part = 0
      DO J=J_0,J_1
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
          SMSI = 0.
          IF (ROICE.gt.0) SMSI =FWSIM(I,J)/ROICE
          TFO  =tfrez(sss(i,j))
C**** get ice-ocean fluxes from sea ice routine
          RUN0=RUNOSI(I,J)-SRUNOSI(I,J) ! fw, includes runoff+basal
          FIDT=ERUNOSI(I,J)
c          SALT=SRUNOSI(I,J)
C**** get river runoff/iceberg melt flux
          RVRRUN =( FLOWO(I,J)+ GMELT(I,J))/(FOCEAN(I,J)*DXYPJ)
          RVRERUN=(EFLOWO(I,J)+EGMELT(I,J))/(FOCEAN(I,J)*DXYPJ)
          OA(I,J,4)=OA(I,J,4)+RVRERUN ! add rvr E to surf. energy budget

          IF (KOCEAN .ge. 1) THEN
            WTRO=Z1O(I,J)*RHOWS
            OTDT=qfluxX*DTSRC*(OTA(I,J,4)*SN4ANG+OTB(I,J,4)*CS4ANG
     *           +OTA(I,J,3)*SN3ANG+OTB(I,J,3)*CS3ANG
     *           +OTA(I,J,2)*SN2ANG+OTB(I,J,2)*CS2ANG
     *           +OTA(I,J,1)*SINANG+OTB(I,J,1)*COSANG+OTC(I,J))

C**** Calculate the amount of ice formation
            CALL OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,FODT,FIDT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFO,WTRW,RUN4O,ERUN4O
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
            AREG_part(J,JR,J_IMPLM)=AREG_part(J,JR,J_IMPLM)+
     *             (RUN4O *POCEAN+RUN4I *POICE)*DXYPJ
            AREG_part(J,JR,J_IMPLH)=AREG_part(J,JR,J_IMPLH)+
     *             (ERUN4O*POCEAN+ERUN4I*POICE)*DXYPJ
            MLHC(I,J)=SHW*WTRW
          ELSE
            ACEFO=0 ; ACEFI=0. ; ENRGFO=0. ; ENRGFI=0.
          END IF

C**** Store mass and energy fluxes for formation of sea ice
C**** Note ACEFO/I is the freshwater amount of ice.
C**** Add salt for total mass
          DMSI(1,I,J)=ACEFO/(1.-SSI0)
          DMSI(2,I,J)=ACEFI/(1.-SSI0)
          DHSI(1,I,J)=ENRGFO
          DHSI(2,I,J)=ENRGFI
          DSSI(1,I,J)=SSI0*DMSI(1,I,J)
          DSSI(2,I,J)=SSI0*DMSI(2,I,J)
#ifdef TRACERS_WATER
C**** assume const mean tracer conc over freshwater amount
          DTRSI(:,1,I,J)=TRSI0(:)*ACEFO
          DTRSI(:,2,I,J)=TRSI0(:)*ACEFI
#endif
C**** store surface temperatures
          GTEMP(1:2,1,I,J)=TOCEAN(1:2,I,J)

        END IF
      END DO
      END DO
      DO JR = 1, NREG
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLM), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+gsum
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLH), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+gsum
      END DO
      RETURN
C****
      END SUBROUTINE OCEANS

      SUBROUTINE ADVSI_DIAG
!@sum  ADVSI_DIAG adjust diagnostics + mlhc for qflux
!@auth Gavin Schmidt
      USE CONSTANT, only : shw,rhows
      USE MODEL_COM, only : focean,im,jm,itocean,itoice,kocean,itime
     *     ,jmon
      USE GEOM, only : dxyp,imaxj
      USE STATIC_OCEAN, only : tocean,z1o,z12o
      USE SEAICE, only : ace1i,lmi
      USE SEAICE_COM, only : rsi,msi,hsi,ssi,snowi
#ifdef TRACERS_WATER
      USE SEAICE_COM, only : trsi
#endif

      USE DIAG_COM, only : NREG, KAJ
      USE FLUXES, only : fwsim,msicnv,mlhc
      USE DIAG_COM, only : aj=>aj_loc,areg,J_IMPLM,J_IMPLH,jreg
     *     ,aij=>aij_loc,j_imelt,j_hmelt,j_smelt
      USE DOMAIN_DECOMP, only : GRID,GET,AM_I_ROOT,GLOBALSUM
      IMPLICIT NONE
      INTEGER I,J,JR
      REAL*8 DXYPJ,RUN4,ERUN4,TGW,POICE,POCEAN,Z1OMIN,MSINEW
      INTEGER :: J_0,J_1
      REAL*8 :: AREG_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO,NREG,KAJ)
      REAL*8 :: gsum

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

      IF (KOCEAN.ge.1) THEN     ! qflux model
        areg_part = 0
      DO J=J_0,J_1
      DXYPJ=DXYP(J)
      DO I=1,IMAXJ(J)
        JR=JREG(I,J)
        POICE=FOCEAN(I,J)*RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-RSI(I,J))
        IF (FOCEAN(I,J).gt.0) THEN
          TGW  = TOCEAN(1,I,J)
          RUN4  = MSICNV(I,J)
          ERUN4 = 0.  ! TGW*SHW*RUN4 ! force energy consveration
C**** Ensure that we don't run out of ocean if ice gets too thick
          IF (POICE.GT.0) THEN
            Z1OMIN=1.+FWSIM(I,J)/(RHOWS*RSI(I,J))
            IF (Z1OMIN.GT.Z1O(I,J)) THEN
C**** MIXED LAYER DEPTH IS INCREASED TO OCEAN ICE DEPTH + 1 METER
              WRITE(6,602) ITime,I,J,JMON,Z1O(I,J),Z1OMIN,z12o(i,j)
 602          FORMAT (' INCREASE OF MIXED LAYER DEPTH ',I10,3I4,3F10.3)
              Z1O(I,J)=MIN(Z1OMIN, z12o(i,j))
              IF (Z1OMIN.GT.Z12O(I,J)) THEN
C****       ICE DEPTH+1>MAX MIXED LAYER DEPTH :
C****       lose the excess mass to the deep ocean
C**** Calculate freshwater mass to be removed, and then any energy/salt
                MSINEW=MSI(I,J)*(1.-RHOWS*(Z1OMIN-Z12O(I,J))/(FWSIM(I,J)
     *               -RSI(I,J)*(ACE1I+SNOWI(I,J)-SUM(SSI(1:2,I,J)))))
                HSI(3:4,I,J) = HSI(3:4,I,J)*(MSINEW/MSI(I,J))
                SSI(3:4,I,J) = SSI(3:4,I,J)*(MSINEW/MSI(I,J))
#ifdef TRACERS_WATER
                TRSI(:,3:4,I,J) = TRSI(:,3:4,I,J)*(MSINEW/MSI(I,J))
#endif
                AJ(J,J_IMELT,ITOICE)=AJ(J,J_IMELT,ITOICE)-FOCEAN(I,J)
     *               *RSI(I,J)*(MSINEW-MSI(I,J))
                AJ(J,J_HMELT,ITOICE)=AJ(J,J_HMELT,ITOICE)-FOCEAN(I,J)
     *               *RSI(I,J)*SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)
                AJ(J,J_SMELT,ITOICE)=AJ(J,J_SMELT,ITOICE)-FOCEAN(I,J)
     *               *RSI(I,J)*SUM(SSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)
                AJ(J,J_IMPLM,ITOICE)=AJ(J,J_IMPLM,ITOICE)-FOCEAN(I,J)
     *               *RSI(I,J)*(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))
     *               /MSI(I,J))
                AJ(J,J_IMPLH,ITOICE)=AJ(J,J_IMPLH,ITOICE)-FOCEAN(I,J)
     *               *RSI(I,J)*SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)
                AREG_part(J,JR,J_IMPLM)=AREG_part(J,JR,J_IMPLM)-
     *               FOCEAN(I,J)*RSI(I,J)
     *               *(MSINEW-MSI(I,J))*(1.-SUM(SSI(3:4,I,J))
     *               /MSI(I,J))*DXYPJ
                AREG_part(J,JR,J_IMPLH)=AREG_part(J,JR,J_IMPLH)-
     *               FOCEAN(I,J)*RSI(I,J)
     *               *SUM(HSI(3:4,I,J))*(MSINEW/MSI(I,J)-1.)*DXYPJ
                MSI(I,J)=MSINEW
                FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)
     *               -SUM(SSI(1:LMI,I,J)))
              END IF
            END IF
          END IF
          MLHC(I,J) = SHW*(Z1O(I,J)*RHOWS-FWSIM(I,J))
C**** Open Ocean diagnostics
          AJ(J,J_IMPLM,ITOCEAN)=AJ(J,J_IMPLM,ITOCEAN)+RUN4 *POCEAN
          AJ(J,J_IMPLH,ITOCEAN)=AJ(J,J_IMPLH,ITOCEAN)+ERUN4*POCEAN
C**** Ice-covered ocean diagnostics
          AJ(J,J_IMPLM,ITOICE)=AJ(J,J_IMPLM,ITOICE)+RUN4 *POICE
          AJ(J,J_IMPLH,ITOICE)=AJ(J,J_IMPLH,ITOICE)+ERUN4*POICE
C**** regional diagnostics
          AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+RUN4 *FOCEAN(I,J)*DXYPJ
          AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+ERUN4*FOCEAN(I,J)*DXYPJ
        END IF
      END DO
      END DO
      DO JR = 1, NREG
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLM), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+gsum
        CALL GLOBALSUM(grid, AREG_part(:,JR,J_IMPLH), gsum)
        IF (AM_I_ROOT()) AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+gsum
      END DO
      END IF

      RETURN
      END

      SUBROUTINE DIAGCO (M)
!@sum  DIAGCO Keeps track of the ocean conservation properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : kocean
      USE DIAG_COM, only : icon_OCE
      IMPLICIT NONE
!@var M index denoting from where DIAGCO is called
      INTEGER, INTENT(IN) :: M
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCO IS BEING CALLED
C****     (see DIAGCA)
      REAL*8, EXTERNAL :: conserv_OCE

C**** OCEAN POTENTIAL ENTHALPY
      IF (KOCEAN.ge.1) CALL conserv_DIAG(M,conserv_OCE,icon_OCE)
C****
      RETURN
      END SUBROUTINE DIAGCO


      SUBROUTINE io_oda(kunit,it,iaction,ioerr)
!@sum  io_oda reads/writes ocean/ice data for initializing deep ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,Itime,im,jm
      USE DIAG_COM, only : ij_tgo2,aij=>aij_loc
      USE SEAICE_COM, only : rsi,msi,hsi,ssi
      USE STATIC_OCEAN, only : tocean
      USE DOMAIN_DECOMP, ONLY : GRID, PACK_DATA, UNPACK_DATA
      USE DOMAIN_DECOMP, ONLY : AM_I_ROOT
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
      INTEGER I,J
      REAL*8 :: AIJ_tmp_glob(IM,JM)

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output
        CALL PACK_DATA(grid, AIJ(:,:,IJ_TGO2), AIJ_tmp_glob)
        IF (AM_I_ROOT()) WRITE (kunit,err=10) it,AIJ_tmp_glob
      CASE (IOREAD:)            ! input
        READ (kunit,err=10) it,((AIJ(I,J,IJ_TGO2),I=1,IM),J=1,JM)
        CALL UNPACK_DATA(grid, AIJ_tmp_glob, AIJ(:,:,IJ_TGO2), 
     *       local=.true.)
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
      USE DOMAIN_DECOMP, only : grid,get
      IMPLICIT NONE
      INTEGER i,j,n
      INTEGER :: j_0,j_1

      call get(grid,j_strt=j_0,j_stop=j_1)

      do n=1,ntm
        if (itime.eq.itime_tr0(n)) then
          do j=j_0,j_1
          do i=1,imaxj(j)
            if (focean(i,j).gt.0) gtracer(n,1,i,j)=trw0(n)
          end do
          end do
        end if
      end do

      return
      end subroutine tracer_ic_ocean
#endif
