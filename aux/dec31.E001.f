C**** For other runs change:  OM70Y to your climatol.ocn.file (15)
C****                         Z12OM65.Y  to your max.MLD file   (25)
C****                         ZM70  to your topography file   (26)
C****                         SEP30 to last day of your run (modify)
C****                         all PARAMETER statements (3)
C**** Changes for B08 M9:  OM70Z,Z12OM250.Y,ZM70,IM=36
C**** Changes for B23 M9:  OM70Z,Z12OM250.Y,ZM70,DEC31,IM=36,JDAY0=5
C**** SEP30.OCN writes an output disk file containing ALL ocean data
C**** on SEP30.  The values are obtained
C**** by integrating in time from Day 1 and applying subroutine
C**** OSTRUC.
C**** Output: TGO = mixed layer temperature
C****       ROICE = ratio ice/water in ocean
C****        ICE2 = ice amount in ocean layer 2
C****   27   TG2O = mean temperature from mixed layer to annual maximum
C****       TG12O = ocean temperature at annual maximum mixed layer
C****         Z1O = current mixed layer depth (on Sep 30)
C****
C**** Input: 15 = OM70Y      climatological ocean data
C****        17 = sea ice data
C****        18 = mixed layer depth
C****        25 = Z12OM65.Y annual maximal mixed layer depths
C****        26 = ZM70       topography
C****
      PARAMETER (IM=72,JM=46,JDAY0=1)
      COMMON /FIXED/  PWATER(IM,JM),DM(IM,JM)
      COMMON /OCLIMD/ TGO(IM,JM),ROICE(IM,JM),ACE2(IM,JM),Z1O(IM,JM),
     *  Z12O(IM,JM),TG2O(IM,JM),TG12O(IM,JM),Z1OOLD(IM,JM)
      REAL ODATA(IM,JM,7),PLAKE(IM,JM),POCN(IM,JM)
      REAL*8 OUTDBL(IM,JM,6)
      EQUIVALENCE (ODATA(1,1,1),TGO(1,1))
      CHARACTER*80 TITLE(7),TITLE0
C****
C**** Contents of OCLIMD
C****
      DATA TITLE/
     1  '  TGO = Ocean Temperature of Mixed Layer (C)',
     2  'ROICE = Ratio of Ocean Ice to Water (1)',
     3  ' ACE2 = Ocean Ice Amount in Second Layer (kg/m**2)',
     4  '  Z1O = MIXED LAYER DEPTH (M)',
     5  ' Z12O = DEPTH OF BOTTOM OF SECOND LAYER (M)',
     6  ' TG2O = OCEAN TEMPERATURE OF SECOND LAYER (C)',
     7  'TG12O = OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)'/
C****
C**** Read in PLAND, calculate PWATER
C****
      READ (26) TITLE0,PLAKE
      write(6,*) TITLE0
      READ (26) TITLE0,POCN
      write(6,*) TITLE0
      REWIND 26
      WRITE (6,*) 'PLAND read in on unit 26, ZM70'
      DO 110 J=1,JM
      DO 110 I=1,IM
  110 PWATER(I,J) = PLAKE(I,J)+POCN(I,J)
C**** Read in aux. sea-ice file
      READ(17) TITLE0,DM  ! don't rewind this file
      WRITE(6,*) TITLE0
C**** Read in Z12O, the annual maximum mixed layer depth
      READ (25) TITLE0,Z12O
      REWIND 25
      WRITE (6,*) TITLE0
C****
C**** Loop over days of the year
C****
      DO 340 JDAY1=JDAY0,JDAY0+364
      JDAY=1+MOD(JDAY1-1,365)
C**** Save yesterday's Z1O in Z1OOLD
      DO 310 J=1,JM
      DO 310 I=1,IM
  310 Z1OOLD(I,J) = Z1O(I,J)
C**** Interpolate daily ocean data from the monthly climatology
      CALL OCLIM (JDAY)
      IF(JDAY1.GT.JDAY0) GO TO 330
C**** Initialize TG2O and TG12O on Day 1
      DO 320 J=1,JM
      DO 320 I=1,IM
      TG2O(I,J)  = TGO(I,J)
  320 TG12O(I,J) = TGO(I,J)
      GO TO 340
C**** Restructure the ocean temperature profile on subsequent days
  330 CALL OSTRUC
  340 CONTINUE
C****
C**** Write Z1O, TG2O and TG12O to a disk file in REAL*8
C****
      DO 350 J=1,JM
      DO 350 I=1,IM
      OUTDBL(I,J,6)=Z1O(I,J)
      OUTDBL(I,J,1)=ODATA(I,J,1)
      OUTDBL(I,J,2)=ODATA(I,J,2)
      OUTDBL(I,J,3)=ODATA(I,J,3)
      OUTDBL(I,J,4)=TG2O(I,J)
  350 OUTDBL(I,J,5)=TG12O(I,J)
      WRITE (27) OUTDBL
      REWIND 27
      WRITE (6,940)
C****
C**** PRODUCE MAPS OF OCEAN DATA ON NOVEMBER 30
C****
      IH=24*(JDAY0-1)
      DO 700 K=1,7
      IF((K.EQ.2).OR.(K.EQ.3))  GO TO 700
      CALL MAP1 (IM,JM,IH,TITLE(K),ODATA(1,1,K),PWATER,1.,0.,26)
  700 CONTINUE
      STOP
C****
  940 FORMAT ('0Z1O, TG2O and TG12O written on unit 2,',
     *  ' SEP30.M65.Z')
      END

      SUBROUTINE OCLIM (JDAY)
C****
C**** OCLIM calculates the daily ocean data by interpolating from
C**** the monthly ocean climatology.
C**** Input: JDAY = Julian day of the year
C****      OCLIMM = monthly ocean climatology data
C**** Output: interpolated ocean climatology
C****         TGO = surface ocean temperature (C)
C****       ROICE = ratio of ocean ice to water (1)
C****        ACE2 = second layer ocean ice amount (kg/m**2)
C****         Z1O = mixed layer depth (m)
C**** Variables: JDAY1 = day associated with earlier month
C****            JDAY2 = day associated with later month
C****
      PARAMETER (IM=72,JM=46)
      INTEGER*4 JDOFM(0:12)
      COMMON /OCLIMM/AOST(IM,JM),EOST1(IM,JM),EOST0(IM,JM),
     *  BOST(IM,JM),COST(IM,JM),ARSI(IM,JM),ERSI1(IM,JM),ERSI0(IM,JM),
     *  BRSI(IM,JM),CRSI(IM,JM),KRSI(IM,JM),Z1O1(IM,JM),Z1O2(IM,JM)
      COMMON /OCLIMD/ TGO(IM,JM),ROICE(IM,JM),ACE2(IM,JM),
     *  Z1O(IM,JM),Z12O(IM,JM)
      COMMON /FIXED/  PWATER(IM,JM),DM(IM,JM)
      CHARACTER*80 TITLE,TITLE1
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA MONTH1/13/,MONTH2/13/,MONTHO/0/
      DATA Z1I/.1/,RHOI/916.6/
C****
C**** Calculate calendar
C****
      DO 110 MONTH=1,12
      IF(JDAY.LE.JDOFM(MONTH))  GO TO 120
  110 CONTINUE
      WRITE (6,*)  'JDAY is not between 1 and 365.'
      STOP 110
  120 JDATE = JDAY-JDOFM(MONTH-1)
C****
C**** Interpolate Mixed layer depths to current day (lin. interp.)
C****
      IF(JDATE.GT.16)  GO TO 130
C**** JDATE is 1 to 16
      IF(MONTH2.EQ.MONTH)  GO TO  400
      IF(MONTH2.EQ.MONTH-1)  GO TO  300
      MONTH1 = MONTH-1
      IF(MONTH.EQ.1)  MONTH1 = 12
      GO TO 200
C**** JDATE is 17 to 31
  130 IF(MONTH1.EQ.MONTH)  GO TO  400
      IF(MONTH2.EQ.MONTH)  GO TO  300
      MONTH1 = MONTH
C****
C**** Ocean climatology file is improperly aligned, read in data
C**** for month MONTH1 into OCLIM1
C****
  200 IF(MONTH1.GT.MONTH2)  GO TO 210
      REWIND 18
      MONTH2 = 0
  210 DO 220 M=MONTH2+1,MONTH1-1
  220 READ (18)
      READ (18) TITLE,Z1O1
      JDAY1 = JDOFM(MONTH1-1)+16
      WRITE (6,*) MONTH1,TITLE
      GO TO 320
C****
C**** Copy OCLIM2 into OCLIM1 and then read the next record
C****
  300 DO 310 J=1,JM
      DO 310 I=1,IM
  310 Z1O1(I,J) = Z1O2(I,J)
      MONTH1 = MONTH2
      JDAY1  = JDAY2
  320 IF(MONTH1.EQ.12)  REWIND 18
      READ (18) TITLE,Z1O2
      MONTH2=MONTH1+1
      IF(MONTH2.GT.12) MONTH2=1
      JDAY2 = JDOFM(MONTH2-1)+16
      WRITE (6,*) MONTH2,TITLE
      IF(JDAY1.GT.JDAY)  JDAY1 = JDAY1-365
      IF(JDAY2.LT.JDAY)  JDAY2 = JDAY2+365
C****
C**** Interpolate data to current day
C****
  400 IF(JDAY1.GT.JDAY)  JDAY1 = JDAY1-365
      IF(JDAY2-365.GT.JDAY)  JDAY2 = JDAY2-365
      X1 = FLOAT(JDAY2-JDAY)/(JDAY2-JDAY1)
      X2 = 1.-X1
      DO 410 J=1,JM
      DO 410 I=1,IM
  410 Z1O(I,J) = X1*Z1O1(I,J)+X2*Z1O2(I,J)
C****
C**** Read in and interpolate other Observed Ocean Data - TGO, ROICE
C****
      IF (MONTH.EQ.MONTHO) GO TO 800
      IF (MONTHO.EQ.0) THEN
C****    READ IN LAST MONTH'S END-OF-MONTH DATA
         LSTMON=MONTH-1
         IF(LSTMON.EQ.0) LSTMON=12
         DO 710 M=1,LSTMON
         READ(15) TITLE,EOST0,EOST0
  710    READ(17) TITLE1,ERSI0,ERSI0
         WRITE(*,*) LSTMON,TITLE
         WRITE(*,*) LSTMON,TITLE1
      ELSE
C****    COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
         DO 720 I=1,IM*JM
         EOST0(I,1)=EOST1(I,1)
  720    ERSI0(I,1)=ERSI1(I,1)
      END IF
C**** Read in current month's data: MEAN and END-OF-MONTH
      MONTHO=MONTH
      IF (MONTH.EQ.1) THEN
         REWIND 15
         REWIND 17
         READ (17)         ! skip over DM-record
      END IF
      READ (15) TITLE,AOST,EOST1
      WRITE(*,*) TITLE
      READ (17) TITLE,ARSI,ERSI1
      WRITE(*,*) TITLE
C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO 730 J=1,JM
      IMAX=IM
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 730 I=1,IMAX
      BOST(I,J)=EOST1(I,J)-EOST0(I,J)
      COST(I,J)=3.*(EOST1(I,J)+EOST0(I,J)) - 6.*AOST(I,J)
      BRSI(I,J)=ERSI1(I,J)-ERSI0(I,J)                ! quadratic fit
      CRSI(I,J)=3.*(ERSI1(I,J)+ERSI0(I,J)) - 6.*ARSI(I,J)
      KRSI(I,J)=0
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
  730 CONTINUE
C**** Calculate TGO, ROICE and ACE2 for current day
  800 TIME=(JDATE-.5)/(JDOFM(MONTH)-JDOFM(MONTH-1))-.5 ! -.5 < TIME < .5
      DO 850 J=1,JM
      ZIMIN=.5
      ZIMAX=2.
      IF(J.GT.JM/2) ZIMAX=3.5
      IMAX=IM
      IF(J.eq.1 .or. J.eq.JM)  IMAX=1
      DO 850 I=1,IMAX
      IF(PWATER(I,J).LE.0.)  GO TO 850
C**** OST always uses quadratic fit
      TGO(I,J) = AOST(I,J) + BOST(I,J)*TIME +
     +           COST(I,J)*(TIME*TIME - 1./12.d0)
      IF(KRSI(I,J)) 810,830,820
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
  810 IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .gt. 0.)  then
        RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5)  !  TIME < T0
      ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .gt. 0.)  then
        RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME)  !  T1 < TIME
      ELSE
        RSINEW = 0.                                !  T0 < TIME < T1
          endif
      GO TO 840
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
  820 IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .lt. 1.)  then
        RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5)  !  TIME < T0
      ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .lt. 1.)  then
        RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME)  !  T1 < TIME
      ELSE
        RSINEW = 1.                                !  T0 < TIME < T1
          endif
      GO TO 840
C**** RSI uses quadratic fit
  830 RSINEW = ARSI(I,J) + BRSI(I,J)*TIME +
     +         CRSI(I,J)*(TIME*TIME - 1./12.d0)
  840 ROICE(I,J)=RSINEW
      ACE2(I,J)=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** WHEN TGO IS NOT DEFINED, MAKE IT A REASONABLE VALUE
      IF (TGO(I,J).LT.-1.8) TGO(I,J)=-1.8
C**** Reduce the ratio of ocean ice to water by .06*(1/ZICE-1/5)
      IF (ROICE(I,J).GT.0.) THEN
         BYZICE=RHOI/(Z1I*RHOI+ACE2(I,J))
         ROICE(I,J)=ROICE(I,J)*(1.-.06*(BYZICE-1./5.))
      END IF
C**** Reduce mixed layer depth to max mixed layer depth
      IF(Z1O(I,J).GT.Z12O(I,J)-.01) Z1O(I,J)=Z12O(I,J)
  850 CONTINUE
      RETURN
      END

      SUBROUTINE OSTRUC
C****
C**** OSTRUC restructures the ocean temperature profile when the
C**** mixed layers depths are changed (normally once a day).
C****
      PARAMETER (IM=72,JM=46)
      COMMON /FIXED/  PWATER(IM,JM),DM(IM,JM)
      COMMON /OCLIMD/ TGO(IM,JM),ROICE(IM,JM),ACE2(IM,JM),
     *  Z1O(IM,JM),Z12O(IM,JM),TG2O(IM,JM),TG12O(IM,JM),Z1OOLD(IM,JM)
      COMMON /OVF/ OA(IM,JM,12)
      DATA RHOW/1000./,RHOI/916.6/,Z1I/.1/
C****
C**** OCLIM  1  Ocean temperature of mixed layer (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****        4  MIXED LAYER DEPTH (M)
C****        5  OCEAN TEMPERATURE OF SECOND LAYER (C)
C****        6  OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)
C****        7  DEPTH OF BOTTOM OF SECOND LAYER (M)
C****
C**** Ocean vertical flux data saved from run B210M9
C****  4  ACE1I+SNOWOI  (instantaneous at noon GMT)
C****
      ACE1I = RHOI*Z1I
C****
C**** RESTRUCTURE OCEAN LAYERS
C****
      DO 200 J=1,JM
      IMAX=IM
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 200 I=1,IMAX
      IF(PWATER(I,J).LE.0.)  GO TO 200
      IF(Z1OOLD(I,J).GE.Z12O(I,J))  GO TO 140
      IF(Z1O(I,J).EQ.Z1OOLD(I,J))  GO TO 200
      WTR1O = RHOW*Z1O(I,J)-ROICE(I,J)*(OA(I,J,1)+ACE2(I,J))
      DWTRO = RHOW*(Z1O(I,J)-Z1OOLD(I,J))
      WTR2O = RHOW*(Z12O(I,J)-Z1O(I,J))
      IF(DWTRO.GT.0.)  GO TO 120
C**** Mixed layer depth is getting shallower
      TG2O(I,J) = TG2O(I,J) - (TGO(I,J)-TG2O(I,J))*DWTRO/WTR2O
      GO TO 200
C**** Mixed layer depth is getting deeper, but TGO is not changed
C     TGAVE = (TG2O(I,J)*DWTRO+(2.*TG2O(I,J)-TG12O(I,J))*WTR2O)
C    *  /(WTR2O+DWTRO)
C     TGO(I,J) = TGO(I,J) + (TGAVE-TGO(I,J))*DWTRO/WTR1O
  120 IF(Z1O(I,J).GE.Z12O(I,J))  GO TO 140
      TG2O(I,J) = TG2O(I,J) + (TG12O(I,J)-TG2O(I,J))*DWTRO/(WTR2O+DWTRO)
      GO TO 200
C**** Mixed layer depth is at its maximum or temp profile is uniform
  140 TG2O(I,J)  = TGO(I,J)
      TG12O(I,J) = TGO(I,J)
  200 CONTINUE
      RETURN
      END
