      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : IM,JM,ioread,iowrite

      IMPLICIT NONE

!@var MWL mass of lake water (kg)
      REAL*8, DIMENSION(IM,JM) :: MWL 
!@var GML total enthalpy of lake (J)
      REAL*8, DIMENSION(IM,JM) :: GML 
!@var TLAKE temperature of lake (C)
      REAL*8, DIMENSION(IM,JM) :: TLAKE 
!@var TFL freezing temperature for lakes (=0 C)
      REAL*8, PARAMETER :: TFL = 0
!@var FLEADLK lead fraction for lakes
      REAL*8, PARAMETER :: FLEADLK=.05d0   ! = 0?
!@var T50 50 day mean temperature (used for estimating lake ice cover)
      REAL*8, SAVE,DIMENSION(IM,JM) :: T50

      END MODULE LAKES_COM

      SUBROUTINE io_lakes(kunit,iaction,ioerr)
!@sum  io_lakes reads and writes lake arrays to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
      USE LAKES_COM
      IMPLICIT NONE
      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "LAKE01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,T50,MWL,TLAKE,GML
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,T50,MWL,TLAKE,GML
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_lakes


