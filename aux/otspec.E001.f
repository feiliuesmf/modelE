C**** For other runs change:  O4X5... to your climatol.ocn.file (15)
C****  (all but one affect    Z12OW250.Y  to your max.MLD file   (25)
C****   printout comments     ZW40  to your topography file   (26)
C****   only)                 NOV30 to last day of your run
C****  ONLY ESSENTIAL CHANGE  all PARAMETER statements (3)
C****                         B05  to run number
C**** Changes for B08 M9:  OM70Z,Z12OM250.Y,ZM70,IM=36,JYEARS=5
C**** Changes for B23 M9:  OM70Z,Z12OM250.Y,ZM70,DEC31,IM=36,JYEARS=5
C**** OTSPEC.OCN reads the ocean vertical fluxes saved from a GCM
C**** run, calculates the ocean energy from the climatological
C**** ocean file, and accumulates the zeroeth and first harmonics
C**** of each of those quantities.  The spectral coefficients are
C**** used to calculate the zeroeth and first harmonics of the
C**** ocean transports.
C****
C**** Input:
C****         4 - XCORR produced by program   vertflux.B05.S
C****        20 - flux data from run B05 W9
C****        15 - O4X5.1979-93means    climatological ocean data
C****        25 - Z12O.B4X5.AMIP.250M  ann max mixed layer depths
C****        26 - Z4X500      topography
C****        27 - dec31.O1-5.MLD.250M.DBL      O(1-5),Z1O at JAN 1
C**** Output:
C****        30 - OTSPEC.RB150.M250D
C****
      PARAMETER (IM=72,JM=46)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COT(IM,JM),AOT(IM,JM,4),BOT(IM,JM,4)
      REAL*4  AMPOT(IM,JM),PHAOT(IM,JM),COTS(IM,JM),DM,TAU4
      REAL*8 TWOPI,ARG,COSDAY(4),SINDAY(4),VFX,
     *  CV(IM,JM),AV(IM,JM,4),BV(IM,JM,4),AE(IM,JM,4),BE(IM,JM,4)
      CHARACTER*80 TITLE(5),TITLE0
      REAL*4 PWATER,PLAKE(IM,JM),POCN(IM,JM),Z12O
      COMMON /FIXED/  PWATER(IM,JM),DM(IM,JM)
      COMMON /OCLIMD/ TGO(IM,JM),ROICE(IM,JM),ACE2(IM,JM),
     *  Z1O(IM,JM),Z12O(IM,JM),TG2O(IM,JM),TG12O(IM,JM),Z1OOLD(IM,JM)
      COMMON /OVF/ OA(IM,JM,12)
C****
      DATA TWOPI/6.283185307179586477D0/
      DATA SHW/4185./, SHI/2060./, ELHM/334000./, RHOW/1000./
      DATA TITLE/
     1  'Amplitude of Ocean Transport (W/m**2) - 1st Harmonic',
     2  'Phase of Ocean Transport (days) - 1st Harmonic',
     3  'Mean Ocean Transport (W/m**2)',
     4  'Cosine coefficient of Ocean Transport (W/m**2)',
     5  'Sine coefficient of Ocean Transport (W/m**2)'/
C****
C**** Ocean vertical flux data saved from run B210M9
C****    K= 1  ACE1I+SNOWOI  (instantaneous at noon GMT)
C****       2  TG1OI  (instantaneous at noon GMT)
C****       3  TG2OI  (instantaneous at noon GMT)
C****       4  ENRGP  (integrated over the day)
C****       5  SRHDT  (integrated over the day)
C****       6  TRHDT  (for ocean, integrated over the day)
C****       7  SHDT   (for ocean, integrated over the day)
C****       8  EVHDT  (for ocean, integrated over the day)
C****       9  TRHDT  (for ocean ice, integrated over the day)
C****      10  SHDT   (for ocean ice, integrated over the day)
C****      11  EVHDT  (for ocean ice, integrated over the day)
C****
      call getarg(1,title0)
      read(title0,*) nskip
      call getarg(2,title0)
      read(title0,*) JYEARS
      READ(4) XCORR
      WRITE(6,*) ' XCORR=',XCORR
C****
C**** Read in input files
C****
      READ(17) TITLE0,DM
C**** Read in PLAND, calculate PWATER
      READ (26) TITLE0,POCN
      write(6,*) TITLE0
      READ (26) TITLE0,PLAKE
      write(6,*) TITLE0
      REWIND 26
      WRITE (6,*) 'PLAND read on unit 26, ZW40'
      DO 110 J=1,JM
      DO 110 I=1,IM
  110 PWATER(I,J) = POCN(I,J)    ! not PLAKE(I,J)+POCN(I,J)
C**** Read in annual maximum mixed layer depth
      READ (25) TITLE0,Z12O
      REWIND 25
      WRITE (6,*) TITLE0
C**** Read in ocean data below mixed layer on December 31
      READ (27) Z1O,Z1O,Z1O,TG2O,TG12O,Z1O
      WRITE(6,*)'Z1O, TG2O, TG12O read on unit 27',Z1O(71,23),Z1O(71,24)
C**** Zero out spectral coefficients
      DO 120 J=1,JM
      DO 120 I=1,IM
      CV(I,J) = 0.
      DO 120 K=1,4
      AV(I,J,K) = 0.
      BV(I,J,K) = 0.
      AE(I,J,K) = 0.
  120 BE(I,J,K) = 0.
C****
C**** Loop over years of data and days in the year
C****
      do 230 n=1,nskip+1
  230 read(20) TAU
      write(*,*) 'first tau used:',TAU
      BACKSPACE (20)
      KDAY1=NINT(TAU/24.)
      IF(KDAY1.EQ.0) KDAY1=365
      DO 320 JYEAR=1,JYEARS
      DO 320 KDAY=KDAY1,KDAY1+364
      JDAY=1+MOD(KDAY-1,365)
      ARG  = JDAY*TWOPI/365.
      DO 205 K=1,4
      COSDAY(K) = DCOS(K*ARG)
  205 SINDAY(K) = DSIN(K*ARG)
C**** Read input tape
      READ (20) TAU,OA
      WRITE (6,920) JDAY,TAU
C**** Interpolate ocean data and restructure temperature profile
      DO 210 J=1,JM
      DO 210 I=1,IM
  210 Z1OOLD(I,J) = Z1O(I,J)
      CALL OCLIM (JDAY)
      CALL OSTRUC
C****
C**** Calculate the vertical flux (J/m**2) and ocean energy (J/m**2)
C****
      DO 310 J=1,JM
      IMAX=IM
      IF((J.EQ.1).OR.(J.EQ.JM))  IMAX=1
      DO 310 I=1,IMAX
      IF(PWATER(I,J).LE.0.)  GO TO 310
      VFX = OA(I,J,4)
     * + (1.-ROICE(I,J))*(OA(I,J,6)+OA(I,J,7)+OA(I,J,8)+XCORR*OA(I,J,5))
     * + ROICE(I,J)*(OA(I,J,9)+OA(I,J,10)+OA(I,J,11)+XCORR*OA(I,J,12))
C     OE = ROICE*(ACE1S*(TG1OI*SHI-ELHM) + ACE2*(TG2OI*SHI-ELHM))
C    *   + ((Z1O*RHOW-(ACE1S+ACE2)*ROICE)*TGO
C    *   +  (Z12O-Z1O)*RHOW              *TG2O)*SHW
      OE = ROICE(I,J)*(OA(I,J,1)*(OA(I,J,2)*SHI-ELHM)
     *               + ACE2(I,J)*(OA(I,J,3)*SHI-ELHM))
     *   + ((Z1O(I,J)*RHOW-(OA(I,J,1)+ACE2(I,J))*ROICE(I,J))*TGO(I,J)
     *   +  (Z12O(I,J)-Z1O(I,J))*RHOW*TG2O(I,J))*SHW
C**** Accumulate the spectral coefficients
      CV(I,J) = CV(I,J) + VFX
      DO 305 K=1,4
      AV(I,J,K) = AV(I,J,K) + VFX*COSDAY(K)
      BV(I,J,K) = BV(I,J,K) + VFX*SINDAY(K)
      AE(I,J,K) = AE(I,J,K) + OE *COSDAY(K)
  305 BE(I,J,K) = BE(I,J,K) + OE *SINDAY(K)
  310 CONTINUE
c     I=71
c     J1=23
c     J2=24
c     write(*,*) I,J1,PWATER(I,J1),ROICE(I,J1),(OA(I,J1,K),K=1,3)
c     write(*,*) I,J2,PWATER(I,J2),ROICE(I,J2),(OA(I,J2,K),K=1,3)
c     write(*,*) I,J1,Z1O(I,J1),TGO(I,J1),Z12O(I,J1),TG2O(I,J1)
c     write(*,*) I,J2,Z1O(I,J2),TGO(I,J2),Z12O(I,J2),TG2O(I,J2)
  320 CONTINUE
      SYEARS = 86400.*365.*JYEARS
C****
C**** SCALE AV TO W/M**2 , AE TO J/M**2 TO CALCULATE SPECTRAL COEFF
C****
      DO 330 J=1,JM
      DO 330 I=1,IM
      CV(I,J) = CV(I,J)/SYEARS
      DO 330 K=1,4
      AV(I,J,K) = AV(I,J,K)*2./SYEARS
      BV(I,J,K) = BV(I,J,K)*2./SYEARS
      AE(I,J,K) = AE(I,J,K)*2./(365.*JYEARS)
  330 BE(I,J,K) = BE(I,J,K)*2./(365.*JYEARS)
C****
C**** Calculate the ocean transports spectral coefficients
C****
      OMEG = TWOPI/(86400.*365.)
      DO 410 J=1,JM
      IMAX=IM
      IF((J.EQ.1).OR.(J.EQ.JM))  IMAX=1
      DO 410 I=1,IMAX
      IF(PWATER(I,J).LE.0.)  GO TO 410
      COT(I,J) =               - CV(I,J)
      DO 405 K=1,4
      AOT(I,J,K) =  BE(I,J,K)*K*OMEG - AV(I,J,K)
  405 BOT(I,J,K) = -AE(I,J,K)*K*OMEG - BV(I,J,K)
  410 CONTINUE
C**** Compute phase and amplitude of ocean transports
      DO 420 J=1,JM
      IMAX=IM
      IF((J.EQ.1).OR.(J.EQ.JM))  IMAX=1
      DO 420 I=1,IMAX
      IF(PWATER(I,J).LE.0.)  GO TO 420
      AMPOT(I,J) = SQRT(AOT(I,J,1)*AOT(I,J,1)+BOT(I,J,1)*BOT(I,J,1))
      PHAOT(I,J) = ATAN2(BOT(I,J,1),AOT(I,J,1))*365./TWOPI
  420 COTS(I,J)  = COT(I,J)
      TAU4 = TAU
      CALL MAP1 (IM,JM,TAU4,TITLE(1),AMPOT,PWATER,1.,0.,26)
      CALL MAP1 (IM,JM,TAU4,TITLE(2),PHAOT,PWATER,1.,0.,26)
      CALL MAP1 (IM,JM,TAU4,TITLE(3),COTS,PWATER,1.,0.,26)
C****
C**** Write ocean transports spectral coefficients to disk
C****
      WRITE (30) BOT,AOT,COT
      STOP
  920 FORMAT (' GCM ocean vertical fluxes read from tape:',I5,F9.1)
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
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 Z12O,DM,PWATER,AOST,EOST1,EOST0,ARSI,ERSI1,ERSI0,Z1O1,Z1O2
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
      BRSI(I,J)=0.                              ! extreme cases
      CRSI(I,J)=0.                              ! ice constant
      KRSI(I,J)=0
      IF(ARSI(I,J).LE.0.) GO TO 730
      IF(ARSI(I,J).GT.1.) GO TO 730
      BRSI(I,J)=ERSI1(I,J)-ERSI0(I,J)                ! quadratic fit
      CRSI(I,J)=3.*(ERSI1(I,J)+ERSI0(I,J)) - 6.*ARSI(I,J)
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
      IF(PWATER(I,J).LE.0.)   GO TO 850
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
C**** DON'T Reduce the ratio of ocean ice to water by .06*(1/ZICE-1/5)
cc    IF (ROICE(I,J).GT.0.) THEN
cc       BYZICE=RHOI/(Z1I*RHOI+ACE2(I,J))
cc       ROICE(I,J)=ROICE(I,J)*(1.-.06*(BYZICE-1./5.))
cc    END IF
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
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 PWATER,Z12O
      COMMON /FIXED/  PWATER(IM,JM)
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
