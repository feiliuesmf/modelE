C**** For other runs change:  both PARAMETER statements
C**** Changes for B08M9/B23M9:  OM70Z,ZM70,IM=36,JYEARS=5
C**** VERTFLUX.OCN reads the ocean data saved from a GCM run,
C**** calculates the ocean vertical flux (with solar correction),
C**** VFX(IM,JM,365), and writes the data onto a disk file.  ROICE
C**** and ACE2 are interpolated in time from the climatological
C**** ocean file for the run (OW40Z)
C****
C**** Input files : 20 = OCEAN.FLUX.DATA from model run
C****               17 = SICE4X5.B.1979-90avg seaice file
C****               26 = Z4X500          topography
C**** Output files:  4 = XCORR
C****               99 = SRCOR-line for rundeck
C****            [ 2/3 = VFX/AVFX  vert flux/ann vflux (commented out) ]
C****
      PARAMETER (IM=72,JM=46)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COSP(JM),DXV(JM),DXYP(JM)
      REAL*8 PWATER,OA(IM,JM,12),VFX(IM,JM,365),AVFX(IM,JM)
      REAL*8 ASR(IM,JM),AVFXSR(IM,JM),GSR,GVFXSR,SYEAR
      REAL*4 PLAKE(IM,JM),POCN(IM,JM),OAS(IM,JM,5),DM
      EQUIVALENCE (OA(1,1,1),OAS(1,1,1),PLAKE(1,1))
      CHARACTER*80 TITLE(4),TITLE0
      COMMON /OCLIMD/ PWATER(IM,JM),DM(IM,JM),ROICE(IM,JM),ACE2(IM,JM)
C****
      DATA SDAY/86400./,SYEAR/31536000./
      DATA TITLE/'Uncorrected annual solar radiation (W/m**2)',
     2           'Vertical flux excluding solar radiation (W/m**2)',
     3           'Uncorrected annual vertical flux (W/m**2)',
     4           'Corrected annual vertical flux (W/m**2), run B113I'/
C****
C**** Contents of OA(I,J,K) saved from the GCM
C****
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
      SYEARS = SYEAR*JYEARS
      READ(17) TITLE0,DM
C****
C**** Calculate spherical geometry
C****
      ONE=1.
      TWOPI  = 8.*ATAN(ONE)
      WRITE(6,*) 'TWOPI',TWOPI
      RADIUS = 6375000.
      DLON = TWOPI/IM
      DLAT = .5*TWOPI/(JM-1)
      FJEQ = .5*(1+JM)
      DO 10 J=2,JM-1
   10 COSP(J) = COS((J-FJEQ)*DLAT)
      COSP(1) = 0.
      COSP(JM)= 0.
      DO 20 J=2,JM
   20 DXV(J) = RADIUS*DLON*.5*(COSP(J-1)+COSP(J))
      DYP = RADIUS*DLAT
      DO 30 J=2,JM-1
   30 DXYP(J) = .5*(DXV(J)+DXV(J+1))*DYP
      DXYP(1) = .25*DXV(2)*DYP
      DXYP(JM)= .25*DXV(JM)*DYP
C****
C**** Read in PLAND and calculate PWATER
C****
      READ (26) TITLE0,POCN
      WRITE(6,*) TITLE0
      READ (26) TITLE0,PLAKE
      WRITE(6,*) TITLE0
      REWIND 26
C     WRITE (6,*) 'PLAND read on unit 26, ZW40'
      DO 110 J=1,JM
      DO 110 I=1,IM
  110 PWATER(I,J) = POCN(I,J)   ! not PLAKE(I,J)+POCN(I,J)
C****
C**** Zero out the vertical flux and its components
C****
      DO 220 J=1,JM
      DO 220 I=1,IM
      ASR(I,J)    = 0.
      AVFXSR(I,J) = 0.
C     DO 210 JDAY=1,365
C 210 VFX(I,J,JDAY) = 0.
  220 CONTINUE
C****
C**** Loop over days of the year and read the tape for each day
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
      READ (20) TAU,OA
      WRITE (6,932) JDAY,TAU
      CALL OCLIM (JDAY)
C****
C**** Accumulate ASR, the daily flux of solar radiation,
C****     and AVFXSR, the vertical flux excluding solar radiation
C****
      DO 310 J=1,JM
      DO 310 I=1,IM
      IF(PWATER(I,J).LE.0.) GO TO 310
      VFSR = (1.-ROICE(I,J))*OA(I,J,5) + ROICE(I,J)*OA(I,J,12)
      VF = OA(I,J,4)
     *  + (1.-ROICE(I,J))*(OA(I,J,6)+OA(I,J,7)+OA(I,J,8))
     *  + ROICE(I,J)*(OA(I,J,9)+OA(I,J,10)+OA(I,J,11))
      ASR(I,J)      = ASR(I,J)      + VFSR
      AVFXSR(I,J)   = AVFXSR(I,J)   + VF
C     VFX(I,J,JDAY) = VFX(I,J,JDAY) + VF
  310 CONTINUE
  320 CONTINUE
C****
C**** DETERMINE THE SOLAR RADIATION CORRECTION FACTOR, XCORR, SO THAT
C**** THE CORRECTED GLOBAL VERTICAL FLUX IS EQUAL TO ZERO
C****
      GSR    = ASR(1,JM)   *PWATER(1,JM)*IM*DXYP(JM)
      GVFXSR = AVFXSR(1,JM)*PWATER(1,JM)*IM*DXYP(JM)
      DO 510 J=2,JM-1
      DO 510 I=1,IM
      GSR    = GSR    + ASR(I,J)   *PWATER(I,J)*DXYP(J)
  510 GVFXSR = GVFXSR + AVFXSR(I,J)*PWATER(I,J)*DXYP(J)
      GSR    = GSR/JYEARS
      GVFXSR = GVFXSR/JYEARS
      XCORR = -GVFXSR/GSR
      WRITE(6,*)  'GSR,GVFXSR,XCORR=',GSR,GVFXSR,XCORR
      WRITE(99,*) '  KOCEAN=1,SRCOR=',XCORR,','
      WRITE(4) XCORR
C****
C**** REREAD THE TAPE TO GET THE SOLAR RADIATION AND CALCULATE THE
C**** CORRECTED VERTICAL FLUX FOR EACH DAY
C****
C     REWIND 20
C     DO 610 JDAY=1,365                   change this section
C     READ (20) TAU,OA                    to not reread the tape
C     DO 610 J=1,JM
C     DO 610 I=1,IM
C     IF(PWATER(I,J).LE.0.) GO TO 610
C     VFSR = (1.-ROICE(I,J))*OA(I,J,5) + ROICE(I,J)*OA(I,J,12)
C     VFX(I,J,JDAY) = (XCORR*VFSR+VFX(I,J,JDAY))/SDAY
C 610 CONTINUE
      DO 620 J=1,JM
      DO 620 I=1,IM
  620 AVFX(I,J) = (XCORR*ASR(I,J)+AVFXSR(I,J))/SYEARS
C****
C**** Write VFX and AVFX onto disk after fixing the poles
C****
C     DO 710 JDAY=1,365
C     DO 710 I=2,IM
C 710 VFX(I,JM,JDAY) = VFX(1,JM,JDAY)
C     WRITE (2) VFX
C     WRITE (6,*) 'VFX written to disk file '
      DO 720 I=2,IM
  720 AVFX(I,JM) = AVFX(1,JM)
C     WRITE (3) TITLE(4),AVFX
C     WRITE (6,*) 'AVFX written to disk file '
C****
C**** Produce maps of ASR, AVFXSR, ASR+AVFXSR, and XCORR*ASR+AVFXSR
C****
      DO 810 J=1,JM
      DO 810 I=1,IM
      OAS(I,J,5) = PWATER(I,J)
      OAS(I,J,4) = AVFX(I,J)
      OAS(I,J,1) = ASR(I,J)/SYEARS
      OAS(I,J,2) = AVFXSR(I,J)/SYEARS
  810 OAS(I,J,3) = (ASR(I,J)+AVFXSR(I,J))/SYEARS
      CALL MAP1 (IM,JM,0.,TITLE(1),OAS(1,1,1),OAS(1,1,5),1.,0.,26)
      CALL MAP1 (IM,JM,0.,TITLE(2),OAS(1,1,2),OAS(1,1,5),1.,0.,26)
      CALL MAP1 (IM,JM,0.,TITLE(3),OAS(1,1,3),OAS(1,1,5),1.,0.,26)
      CALL MAP1 (IM,JM,0.,TITLE(4),OAS(1,1,4),OAS(1,1,5),1.,0.,26)
C****
      STOP
  932 FORMAT (' Data read from GCM output tape.  JDAY,TAU=',I6,F8.1)
      END

      SUBROUTINE OCLIM (JDAY)
C****
C**** OCLIM calculates the daily ocean data by interpolating from
C**** the monthly ocean climatology.
C**** Input:    JDAY = Julian day of the year
C****         OCLIMM = monthly ocean climatology data
C**** Output:  ROICE = ratio of ocean ice to water (1)
C****           ACE2 = second layer ocean ice amount (kg/m**2)
C****
      PARAMETER (IM=72,JM=46)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 JDOFM(0:12),KRSI(IM,JM)
      REAL*4 DM,ARSI(IM,JM),ERSI1(IM,JM),ERSI0(IM,JM),BRSI,CRSI
      CHARACTER*80 TITLE
      COMMON /OCLIMM/ARSI,ERSI1,ERSI0,BRSI(IM,JM),CRSI(IM,JM)
      COMMON /OCLIMD/ PWATER(IM,JM),DM(IM,JM),ROICE(IM,JM),ACE2(IM,JM)
      DATA JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA MONTHO/0/
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
C**** READ IN OBSERVED OCEAN DATA
      IF(MONTH.EQ.MONTHO) GO TO 400
      IF(MONTHO.EQ.0) THEN
C****    READ IN LAST MONTH'S END-OF-MONTH DATA
         MONTH1 = MONTH-1
         IF(MONTH.EQ.1)  MONTH1 = 12
         DO 130 MM=1,MONTH1
  130    READ(17) TITLE,ERSI0,ERSI0
         WRITE(*,*) 'read:',TITLE,' EofMonth'
      ELSE
C****    COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
         DO 140 I=1,IM*JM
  140    ERSI0(I,1)=ERSI1(I,1)
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      MONTHO=MONTH
      IF (MONTH.EQ.1) THEN
         REWIND 17
         READ (17) ! skip over DM-record
      END IF
      READ (17) TITLE,ARSI,ERSI1
      WRITE(*,*) 'read:',TITLE,' Av&EofMonth'
C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO 330 J=1,JM
      IMAX=IM
      IF(J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 330 I=1,IMAX
      BRSI(I,J)=0.                              ! extreme cases
      CRSI(I,J)=0.                              ! ice constant
      KRSI(I,J)=0
      IF(ARSI(I,J).LE.0.) GO TO 330
      IF(ARSI(I,J).GT.1.) GO TO 330
      BRSI(I,J)=ERSI1(I,J)-ERSI0(I,J)           ! quadratic fit
      CRSI(I,J)=3.*(ERSI1(I,J)+ERSI0(I,J)) - 6.*ARSI(I,J)
      IF(ABS(CRSI(I,J)) .GT. ABS(BRSI(I,J))) THEN    ! linear fits
        RSICSQ=CRSI(I,J)*
     *  (ARSI(I,J)*CRSI(I,J) - .25*BRSI(I,J)**2 - CRSI(I,J)**2/12.)
        IF(RSICSQ.lt.0.)  then
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
          KRSI(I,J) = -1
          BRSI(I,J) = .5*(ERSI0(I,J)**2 + ERSI1(I,J)**2) / ARSI(I,J)
        ELSE IF(RSICSQ.gt.CRSI(I,J)**2)  then
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
          KRSI(I,J) = 1
          BRSI(I,J) = .5*((ERSI0(I,J)-1.)**2 + (ERSI1(I,J)-1.)**2) /
     /                (ARSI(I,J)-1.)
        END IF
      END IF
  330 CONTINUE
C****
C**** Calculate RSI and MSI for current day
C****
  400 TIME=(JDATE-.5)/(JDOFM(MONTH)-JDOFM(MONTH-1))-.5 ! -.5 < TIME < .5
      DO 450 J=1,JM
      ZIMIN=.5
      ZIMAX=2.
      IF(J.GT.JM/2) ZIMAX=3.5
      IMAX=IM
      IF(J.eq.1 .or. J.eq.JM)  IMAX=1
      DO 450 I=1,IMAX
      IF(PWATER(I,J).LE.0.) GO TO 450
      IF(KRSI(I,J)) 410,430,420
C**** RSI uses piecewise linear fit because quadratic fit at apex < 0
  410 IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .gt. 0.)  then
        RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5)  !  TIME < T0
      ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .gt. 0.)  then
        RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME)  !  T1 < TIME
      ELSE
        RSINEW = 0.
      END IF
      GO TO 440
C**** RSI uses piecewise linear fit because quadratic fit at apex > 1
  420 IF(ERSI0(I,J)-BRSI(I,J)*(TIME+.5) .lt. 1.)  then
        RSINEW = ERSI0(I,J) - BRSI(I,J)*(TIME+.5)  !  TIME < T0
      ELSEIF(ERSI1(I,J)-BRSI(I,J)*(.5-TIME) .lt. 1.)  then
        RSINEW = ERSI1(I,J) - BRSI(I,J)*(.5-TIME)  !  T1 < TIME
      ELSE
        RSINEW = 1.
      END IF
      GO TO 440
C**** RSI uses quadratic fit
  430 RSINEW = ARSI(I,J) + BRSI(I,J)*TIME +
     +         CRSI(I,J)*(TIME*TIME - 1./12.d0)
  440 ROICE(I,J)=RSINEW
      ACE2(I,J)=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
C**** DON'T Reduce the ratio of ocean ice to water by .06*(1/ZICE-1/5)
cc    IF(ROICE(I,J).LE.0.)  GO TO 450
cc    BYZICE=RHOI/(Z1I*RHOI+ACE2(I,J))
cc    ROICE(I,J) = ROICE(I,J)*(1.-.06*(BYZICE-1./5.))
  450 CONTINUE
      RETURN
      END
