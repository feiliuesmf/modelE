      MODULE RADNCB
!@sum  RADNCB Model radiation arrays and parameters
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : im,jm,lm
      USE RE001, only : S0
!@var S0 solar constant needs to be saved between calls to radiation 
      IMPLICIT NONE
      SAVE

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3 
!@var RQT Radiative equilibrium temperatures above model top
      DOUBLE PRECISION, DIMENSION(LM_REQ,IM,JM) :: RQT
!@var SRHR Solar radition absorbed every hour (W/m^2)
!@var TRHR Thermal radition absorbed every hour (W/m^2)
      DOUBLE PRECISION, DIMENSION(LM+1,IM,JM) :: SRHR,TRHR
!@var FSF Solar Forcing over each type (W/m^2)
      DOUBLE PRECISION, DIMENSION(4,IM,JM) :: FSF
!@var COSZ1 Integral of Solar Zenith angle over time (????)
      DOUBLE PRECISION, DIMENSION(IM,JM) :: COSZ1
!@var S0X solar constant multiplication factor (set-able in NAMELIST)
      REAL*8 :: S0X = 1.
!@var CO2 carbon dioxide multiplication factor (set-able in NAMELIST)
      REAL*8 :: CO2 = 1.
!@var RSDIST,SIND,COSD orbit related variables
      REAl*8 :: RSDIST,SIND,COSD 

      END MODULE RADNCB

      SUBROUTINE io_rad(kunit,iaction,ioerr)
!@sum  io_rad reads and writes radiation arrays to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite,irestart,irsfic,irerun
      USE RADNCB
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "RAD01"
!@var RADDUM dummy variable
      REAL*8, DIMENSION(6) :: RADDUM

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,S0,S0X,CO2,RSDIST,SIND,COSD,
     *     RQT,SRHR,TRHR,FSF
      CASE (IOREAD:)
        SELECT CASE  (IACTION)
        CASE (IRESTART,IRERUN)  ! input for restart, rerun or extension
          READ (kunit,err=10) HEADER,S0,S0X,CO2,RSDIST,SIND,COSD,
     *       RQT,SRHR,TRHR,FSF
        CASE (IRSFIC)           ! start from old restart file
          READ (kunit,err=10) HEADER,RADDUM,RQT
        END SELECT
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_rad

