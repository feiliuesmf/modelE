      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,JM,ioread,iowrite,lhead

      IMPLICIT NONE
      SAVE
!@var MWL mass of lake water (kg)
      REAL*8, DIMENSION(IM,JM) :: MWL
!@var GML total enthalpy of lake (J)
      REAL*8, DIMENSION(IM,JM) :: GML
!@var TLAKE temperature of lake (C)
      REAL*8, DIMENSION(IM,JM) :: TLAKE
!@var MLDLK mixed layer depth in lake (m)
      REAL*8, DIMENSION(IM,JM) :: MLDLK
!@var FLAKE variable lake fraction (1)
      REAL*8, DIMENSION(IM,JM) :: FLAKE
!@var TANLK tan(alpha) = slope for conical lake (1)
      REAL*8, DIMENSION(IM,JM) :: TANLK

      END MODULE LAKES_COM

      SUBROUTINE io_lakes(kunit,iaction,ioerr)
!@sum  io_lakes reads and writes lake arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE LAKES_COM
      IMPLICIT NONE
      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "LAKE01"

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm):MixLD,MWtr,Tlk,Enth'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,MLDLK,MWL,TLAKE,GML !,FLAKE
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,MLDLK,MWL,TLAKE,GML   !,FLAKE
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_lakes


