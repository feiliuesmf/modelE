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
      USE OCEAN 
      USE DAGCOM, only : OA
      USE SEAICE_COM, only : rsi,msi, snowi   
      USE SEAICE, only : ace1i 
      USE GEOM  
      implicit none
      integer first_month, first_year, last_month, 
     *        last_year, jyears, jyear, i, j, k, imax
      integer months, monthe, itime1, month, kday, last_day
      REAL*8 COT(IM,JM),AOT(IM,JM,4),BOT(IM,JM,4)
      REAL*4  AMPOT(IM,JM),PHAOT(IM,JM),COTS(IM,JM),TAU4
      REAL*8 ARG,COSDAY(4),SINDAY(4),VFX, XCORR, SYEARS, OMEG, OE, 
     *  CV(IM,JM),AV(IM,JM,4),BV(IM,JM,4),AE(IM,JM,4),BE(IM,JM,4)
      CHARACTER*80 TITLE(5),TITLE0, RunID, file_name
      REAL*4 PWATER(IM,JM), month_day(12)
      REAL*8 TGO(IM,JM),ROICE(IM,JM),ACE2(IM,JM),
     *       TG2O(IM,JM),TG12O(IM,JM)   
C****
      character*4 month_name(12), tmonth, tyear  
      data month_name /'JAN ','FEB ','MAR ','APR ','MAY ','JUN ',
     *                 'JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/ 
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/ 

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
      call getarg(1,RunID )
      call getarg(2,title0)
      read(title0,*) first_month
      call getarg(3,title0)
      read(title0,*) first_year
      call getarg(4,title0)
      read(title0,*) last_month
      call getarg(5,title0)
      read(title0,*) last_year
      jyears = last_year-first_year+1 
C* 
      READ(4) XCORR
      WRITE(6,*) ' XCORR=',XCORR
C****
C**** Read in input files
C****
      CALL READT (17,0,DM,IM*JM,DM,1)
C* 
C**** Read in FOCEAN - ocean fraction 
C* 
      CALL READT (26,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      REWIND 26
C* 
C**** Read in annual maximum mixed layer depth
C* 
      CALL READT (25,0,Z12O,IM*JM,Z12O,1)
      REWIND 25
C* 
C**** Read in ocean data below mixed layer on December 31
C*
      READ (27) Z1O,Z1O,Z1O,TG2O,TG12O,Z1O
      WRITE(6,*)'Z1O, TG2O, TG12O read on unit 27',Z1O(71,23),Z1O(71,24)
C* 
C**** Zero out spectral coefficients
C* 
      DO 120 J=1,JM
      DO 120 I=1,IM
      PWATER(I,J) = FOCEAN(I,J) 
      FLAND (I,J) = 1. - FOCEAN(I,J) 
      CV(I,J) = 0.
      DO 120 K=1,4
      AV(I,J,K) = 0.
      BV(I,J,K) = 0.
      AE(I,J,K) = 0.
  120 BE(I,J,K) = 0.
C****
C**** Calculate spherical geometry
C****
      call GEOM_B
C****
C**** Loop over years of data and days in the year
C****
      DO 320 JYEAR=first_year,last_year
         months = 1 
         monthe = 12 
         if (jyear .eq. first_year) months = first_month 
         if (jyear .eq. last_year)  monthe = last_month 
         jday = JDendOfM(first_month-1)
         do month = months, monthe 
            tmonth = month_name(month) 
            write (tyear, '(i4.4)') jyear
            file_name= 
     *  '/u/cmrun/'//trim(RunID)//'/VFLXO'//trim(tmonth)//tyear
            open (20,file=file_name,form='unformatted')
            last_day = month_day(month) 
            do kday = 1,last_day
               jday = jday + 1 
               ARG  = JDAY*TWOPI/365.
               DO K=1,4
                  COSDAY(K) = DCOS(K*ARG)
                  SINDAY(K) = DSIN(K*ARG) 
               END DO
C*
               READ (20, END=555) itime,OA,itime1 
C* 
C**** Interpolate daily ocean data from the monthly climatology
C*
               kocean = 0 
               jmon = month 
               jdate = kday
               iu_OSST =  15 
               iu_SICE =  17 
               iu_OCNML = 18
               CALL OCLIM (1)
C*             
               do j = 1,jm 
                  do i = 1,im 
                     roice(i,j) = rsi(i,j) 
                     ace2(i,j)  = msi(i,j) 
                     tgo(i,j) = tocean(1,i,j) 
                  end do 
               end do 
C* 
C***  Read in the ocean mixed layer depth data 
C***  and interpolate them for the current day 
C* 
               kocean = 1 
               CALL OCLIM (1)          
C* 
               DO J = 1,JM 
                  DO I = 1,IM 
                     SNOWI(I,J) = OA(I,J,1) - ACE1I  
                  END DO
               END DO
           CALL OSTRUC(.false.)
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
C* 
      OE = ROICE(I,J)*(OA(I,J,1)*(OA(I,J,2)*SHI-LHM)
     *               + ACE2(I,J)*(OA(I,J,3)*SHI-LHM))
     *   + ((Z1O(I,J)*RHOW-(OA(I,J,1)+ACE2(I,J))*ROICE(I,J))*TGO(I,J)
     *   +  (Z12O(I,J)-Z1O(I,J))*RHOW*TG2O(I,J))*SHW
C* 
C**** Accumulate the spectral coefficients
C* 
      CV(I,J) = CV(I,J) + VFX
      DO 305 K=1,4
      AV(I,J,K) = AV(I,J,K) + VFX*COSDAY(K)
      BV(I,J,K) = BV(I,J,K) + VFX*SINDAY(K)
      AE(I,J,K) = AE(I,J,K) + OE *COSDAY(K)
  305 BE(I,J,K) = BE(I,J,K) + OE *SINDAY(K)
  310 CONTINUE
      end do
      end do
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
      TAU4 = itime
      CALL MAP1 (IM,JM,TAU4,TITLE(1),AMPOT,PWATER,1.,0.,26)
      CALL MAP1 (IM,JM,TAU4,TITLE(2),PHAOT,PWATER,1.,0.,26)
      CALL MAP1 (IM,JM,TAU4,TITLE(3),COTS,PWATER,1.,0.,26)
C****
C**** Write ocean transports spectral coefficients to disk
C****
      WRITE (30) BOT,AOT,COT
      STOP
 555  write (*,*) ' Reached end of file ',file_name  
      END

      SUBROUTINE CHECK3(A,IN,JN,LN,SUBR,FIELD)
!@sum  CHECK3 Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: IN,JN,LN
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*2, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(IN,JN,LN),INTENT(IN) :: A

      INTEGER I,J,L !@var I,J,L loop variables

      DO L=1,LN
         DO J=1,JN
            DO I=1,IN
               IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.)) THEN
                  WRITE (6,*) FIELD,': ',I,J,L,A(I,J,L),'after '
     *                 ,SUBR
                  IF (J.LT.JN.AND.J.GT.1) STOP 'CHECK3'
               END IF
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE CHECK3


