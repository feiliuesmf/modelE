      MODULE LANDICE
!@sum  LANDICE contains variables and routines dealing with land ice
!@auth Original Development Team
!@ver  1.0
!@cont PRECLI,LNDICE
      USE CONSTANT, only : lhm,rhoi,shi
      IMPLICIT NONE

!@var Z1E,Z2LI first and second layer thicknesses for land ice
      REAL*8, PARAMETER :: Z1E=.1d0, Z2LI=2.9d0
!@var ACE1LI first layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE1LI = Z1E*RHOI
!@var HC1LI heat capacity of first layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC1LI = ACE1LI*SHI
!@var ACE2LI second layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE2LI = Z2LI*RHOI
!@var HC2LI heat capacity of second layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC2LI = ACE2LI*SHI

      CONTAINS

      SUBROUTINE PRECLI(SNOW,TG1,TG2,PRCP,EPRCP,TPRCP,EDIFS,DIFS,ERUN2
     *     ,RUN0)
!@sum  PRECLI apply precipitation to land ice fraction
!@auth Original Development Team
!@ver  1.0
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: PRCP,EPRCP,TPRCP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,ERUN2,RUN0

      REAL*8 HC1,DWATER
C**** initiallize output
      EDIFS=0. ; DIFS=0. ; ERUN2=0. ; RUN0=0.

      HC1=HC1LI+SNOW*SHI
      IF (TPRCP.GE.0.) THEN
        IF (EPRCP.GE.-TG1*HC1) THEN
C**** RAIN HEATS UP TG1 TO FREEZING POINT AND MELTS SOME SNOW OR ICE
          DWATER=(TG1*HC1+EPRCP)/LHM
          TG1=0.
          RUN0=DWATER+PRCP
          IF (DWATER.LE.SNOW) THEN
C**** RAIN MELTS SOME SNOW
            SNOW=SNOW-DWATER
          ELSE
C**** RAIN MELTS ALL SNOW AND SOME ICE, ICE MOVES UP THROUGH THE LAYERS
            DIFS=SNOW-DWATER
            SNOW=0.
            TG1=-TG2*DIFS/ACE1LI
            EDIFS=DIFS*(TG2*SHI-LHM)
            ERUN2=EDIFS
          END IF
        ELSE
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1
          TG1=TG1+EPRCP/HC1
          RUN0=PRCP
        END IF
      ELSE
C**** SNOW INCREASES SNOW AMOUNT AND SNOW TEMPERATURE RECOMPUTES TG1
        TG1=(TG1*HC1+EPRCP)/(HC1+PRCP*SHI)
        SNOW=SNOW+PRCP
        IF (SNOW.gt.ACE1LI) THEN
C**** SNOW IS COMPACTED INTO ICE, ICE MOVES DOWN THROUGH THE LAYERS
          DIFS=SNOW-.9d0*ACE1LI
          SNOW=.9d0*ACE1LI
          EDIFS=DIFS*(TG1*SHI-LHM)
          ERUN2=DIFS*(TG2*SHI-LHM)
          TG2=TG2+(TG1-TG2)*DIFS/ACE2LI
        END IF
      END IF
      RETURN
      END SUBROUTINE PRECLI

      SUBROUTINE LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0)
!@sum  LNDICE apply surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: F0DT,F1DT,EVAP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,RUN0

      REAL*8 :: SNANDI,HC1,ENRG1
C**** initiallize output
      EDIFS=0. ; DIFS=0. ; RUN0=0.

C**** CALCULATE TG1
      SNANDI=SNOW+ACE1LI-EVAP
      HC1=SNANDI*SHI
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT
      IF (ENRG1.GT.-TG1*HC1) THEN
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
        RUN0=(ENRG1+TG1*HC1)/LHM
        TG1=0.
        SNANDI=SNANDI-RUN0
      ELSE
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT
        TG1=TG1+ENRG1/HC1
      END IF
      IF (SNANDI.LT.ACE1LI) THEN
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2
        SNOW=0.
        DIFS=SNANDI-ACE1LI
        TG1=(TG1*SNANDI-TG2*DIFS)/ACE1LI
        EDIFS=DIFS*(TG2*SHI-LHM)
      ELSE
        SNOW=SNANDI-ACE1LI
      END IF
C**** CALCULATE TG2
      TG2=TG2+F1DT/HC2LI

      RETURN
      END SUBROUTINE LNDICE

      END MODULE LANDICE

      MODULE LANDICE_COM
!@sum  LANDICE_COM contains the model arrays for land ice
!@auth Gavin Schmidt
!@ver  1.0
!@cont io_landice
      USE E001M12_COM, only : im,jm,ioread,iowrite

      IMPLICIT NONE
!@var SNOWLI snow amount on land ice (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: SNOWLI
!@var TLANDI temperature of each land ice layer (C)
      REAL*8, DIMENSION(2,IM,JM) :: TLANDI

      END MODULE LANDICE_COM

      SUBROUTINE io_landice(kunit,iaction,ioerr)
!@sum  io_diag reads and writes diagnostics to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
      USE LANDICE_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "GLAIC01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,SNOWLI,TLANDI
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,SNOWLI,TLANDI
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_landice


