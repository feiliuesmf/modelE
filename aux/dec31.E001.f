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
      USE OCEAN 
      USE SEAICE_COM, only : rsi,msi,snowi 
      USE FILEMANAGER
      implicit none 
      integer i, j, k, last_day, kday, jday0, IH,   
     *     months, monthe, month,iu_TOPO,iu_SICE,iu_MLMAX,iu_OSST
     *     ,iu_OCNML

      REAL*4 TGO(IM,JM),TG2O(IM,JM),TG12O(IM,JM),
     *       ROICE(IM,JM),ACE2(IM,JM)
      REAL*4 ODATA(IM,JM,7),PWATER(im,jm)
      REAL*8 OUTDBL(IM,JM,6)
      REAL*4 month_day(12)
      CHARACTER*80 TITLE(7),TITLE0
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/ 
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
C**** Read in FOCEAN - ocean fraction
C****
      call getunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      REWIND iu_TOPO
C* 
      DO J = 1,JM
         DO I = 1,IM
            PWATER(I,J) = FOCEAN(I,J)
            fland (i,j) = 1.- focean(i,j) 
         end do
      end do
C* 
C**** Read in aux. sea-ice file
C* 
      call getunit("SICE",iu_SICE,.true.,.true.)
      CALL READT (iu_SICE,0,DM,IM*JM,DM,1)
C*
C**** Read in Z12O, the annual maximum mixed layer depth
C* 
      call getunit("MLMAX",iu_MLMAX,.true.,.true.)
      CALL READT (iu_MLMAX,0,Z12O,IM*JM,Z12O,1)
      REWIND iu_MLMAX
C****
C**** Calculate spherical geometry
C****
      call GEOM_B
C**** set up unit numbers for ocean climatologies
      call getunit("OSST",iu_OSST,.true.,.true.)
C**** Set up unit number of mixed layer depth climatogies
      call getunit("OCNML",iu_OCNML,.true.,.true.)
C****
C**** Loop over days of the year
C****
      jday = 0 
      months = 1 
      monthe = 12 
      do month = months, monthe 
         last_day = month_day(month) 
         do kday = 1,last_day
C*
C**** Interpolate daily ocean data from the monthly climatology
C*
         kocean = 0
         jmon = month 
         jdate = kday
         CALL OCLIM (1)
c* 
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
         jday = jday + 1  
         CALL OCLIM (1)          
C* 
C***  Read in ocean ice snow(+ 10 cm ice) data
C* 
         READ(77) TITLE0,SNOWI 
         WRITE (6,*) TITLE0
C* 
         IF(month.eq.1 .and. kday.eq.1) THEN
C*
C**** Initialize TG2O and TG12O on Day 1
C* 
           DO J = 1,JM
           DO I = 1,IM
              TG2O(I,J) = TGO(I,J) 
              TG12O(I,J) = TGO(I,J) 
           END DO
           END DO
         ELSE 
C*
C**** Restructure the ocean temperature profile on subsequent days
C* 
           CALL OSTRUC(.false.)
C* 
           DO J = 1,JM
           DO I = 1,IM
              TOCEAN(1,I,J) = TGO(I,J) 
              TOCEAN(2,I,J) = TG2O(I,J)
              TOCEAN(3,I,J) = TG12O(I,J) 
           END DO
           END DO
         END IF 
         end do
      end do
      CLOSE(77) 
C****
C**** Write Z1O, TG2O and TG12O to a disk file in REAL*8
C****
      DO 350 J=1,JM
      DO 350 I=1,IM
      ODATA(I,J,1) =TOCEAN(1,I,J) 
      ODATA(I,J,4) =Z1O(I,J) 
      ODATA(I,J,5) =Z12O(I,J) 
      ODATA(I,J,6) =TOCEAN(2,I,J) 
      ODATA(I,J,7) =TOCEAN(3,I,J) 
      OUTDBL(I,J,6)=Z1O(I,J)
      OUTDBL(I,J,1)=TOCEAN(1,I,J)
      OUTDBL(I,J,2)=ODATA(I,J,2)
      OUTDBL(I,J,3)=ODATA(I,J,3)
      OUTDBL(I,J,4)=TOCEAN(2,I,J)
  350 OUTDBL(I,J,5)=TOCEAN(3,I,J)
      WRITE (27) OUTDBL
      REWIND 27
      WRITE (6,940)
C****
C**** PRODUCE MAPS OF OCEAN DATA ON NOVEMBER 30
C****
      jday0=1
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


