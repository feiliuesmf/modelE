      program vertflux
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
C****               26 = Z72X46N.cor4          topography
C**** Output files:  4 = XCORR
C****               99 = SRCOR-line for rundeck
C****            [ 2/3 = VFX/AVFX  vert flux/ann vflux (commented out) ]
C****
      USE OCEAN 
      USE DAGCOM, only : OA
      USE SEAICE_COM, only : rsi,msi, snowi 
      USE SEAICE, only : ace1i 
      USE GEOM  
      implicit none 
      integer first_month, first_year, last_month, last_year, jyears
      integer months, monthe, jyear, month, itime1, 
     *        last_day, kday, i, j 
      REAL*8 PWATER(im,jm),ROICE(IM,JM),AVFX(IM,JM)  
      REAL*8 GSR,GVFXSR,SYEAR,SYEARS
      real*8 VFSR, VF, XCORR 
      REAL*8 ASR(im,jm),AVFXSR(im,jm)
      REAL*4 OAS(IM,JM,5), month_day(12)
      CHARACTER*80 TITLE(4),TITLE0, RunID, file_name
      character*4 month_name(12), tmonth, tyear  
      character*2 tday 
      data month_name /'JAN ','FEB ','MAR ','APR ','MAY ','JUN ',
     *                 'JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/ 
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/ 
C****
      DATA SYEAR/31536000./
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
      SYEARS = SYEAR*JYEARS
      CALL READT (17,0,DM,IM*JM,DM,1)
C****
C**** Calculate spherical geometry
C****
      call GEOM_B
C****
C**** Read in FOCEAN - ocean fraction
C****
      CALL READT (26,0,PWATER,IM*JM,PWATER,1) ! Ocean fraction
      REWIND 26
C****
C**** Zero out the vertical flux and its components
C****
      DO 220 J=1,JM
      DO 220 I=1,IM
      ASR(I,J)    = 0.
      AVFXSR(I,J) = 0.
C*
      focean(i,j) = pwater(i,j) 
      fland (i,j) = 1.-focean(i,j) 
  220 CONTINUE
C****
C**** Loop over days of the year and read the data for each day
C****
      DO 320 JYEAR=first_year,last_year
         months = 1 
         monthe = 12 
         if (jyear .eq. first_year) months = first_month 
         if (jyear .eq. last_year)  monthe = last_month 
         do month = months, monthe 
            tmonth = month_name(month) 
            write (tyear, '(i4.4)') jyear
            file_name= 
     *  '/u/cmrun/'//trim(RunID)//'/VFLXO'//trim(tmonth)//tyear
            open (20,file=file_name,form='unformatted')
            last_day = month_day(month) 
            do kday = 1,last_day
               READ (20, END=555) itime,OA,itime1 
C*
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
                  end do 
               end do 
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
  310 CONTINUE
C* 
C**** WRITE OUT OA(I,J,1) FOR THE LAST YEAR  
C* 
      IF (JYEAR .EQ. LAST_YEAR) THEN 
         write (tday, '(i2.2)') kday
         TITLE0='SNOW AMOUNT ON SEA ICE, '//tyear//', '//tmonth//
     *   ', '//tday//', (kg/m^2)'
C*
         DO J = 1,JM 
         DO I = 1,IM 
            SNOWI(I,J) = OA(I,J,1) - ACE1I 
         END DO
         END DO 
C* 
         write(*,*) title0
         WRITE(77) TITLE0,SNOWI 
      end if 
      end do
      end do 
  320 CONTINUE
C* 
      CLOSE(77) 
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
      DO 620 J=1,JM
      DO 620 I=1,IM
  620 AVFX(I,J) = (XCORR*ASR(I,J)+AVFXSR(I,J))/SYEARS
C****
      DO 720 I=2,IM
  720 AVFX(I,JM) = AVFX(1,JM)
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
 555  write (*,*) ' Reached end of file ',file_name  
c      write (*,*) ' jyear= ',jyear,' month= ',month,' day= ',kday  
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









