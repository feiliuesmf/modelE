      MODULE GHYCOM
!@sum  GHYCOM contains the areas used by the Ground Hydrology routines
!@auth Frank Abramopolus/Igor Aleinov
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE SLE001, only : ngm,imt,nlsn
      IMPLICIT NONE
      SAVE
C bare/veg not in merged array because WBARE does not contain
C 0 index for legacy reasons
      DOUBLE PRECISION, DIMENSION(  NGM,IM,JM) :: WBARE
      DOUBLE PRECISION, DIMENSION(0:NGM,IM,JM) :: WVEGE
      DOUBLE PRECISION, DIMENSION(0:NGM,IM,JM) :: HTBARE
      DOUBLE PRECISION, DIMENSION(0:NGM,IM,JM) :: HTVEGE
      DOUBLE PRECISION, DIMENSION(2,IM,JM) :: SNOWBV

      DOUBLE PRECISION, DIMENSION(IM,JM,NGM) :: DZ_IJ
      DOUBLE PRECISION, DIMENSION(IM,JM,IMT,NGM) :: Q_IJ,QK_IJ
      DOUBLE PRECISION, DIMENSION(IM,JM) :: SL_IJ
C this common only for purpose of reading its contents from a
C file opened in fortran unformatted sequential access mode
C containing its contents in a contiguous real*4 block
      COMMON/SDATA/ DZ_IJ,Q_IJ,QK_IJ,SL_IJ

      DOUBLE PRECISION, DIMENSION(NGM,IM,JM) :: AFR
      DOUBLE PRECISION, DIMENSION(3,IM,JM) :: ALA,ACS
      DOUBLE PRECISION, DIMENSION(IM,JM) :: AFB,AVH

ccc the following arrays contain prognostic variables for the snow model
ccc ( ISN can be eliminated later, since FR_SNOW contains similar info )
      INTEGER, DIMENSION(2,IM,JM)     :: NSN_IJ
      INTEGER, DIMENSION(2,IM,JM)     :: ISN_IJ   
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: DZSN_IJ
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: WSN_IJ
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: HSN_IJ
      REAL*8, DIMENSION(2,IM,JM)      :: FR_SNOW_IJ
C**** replacements for GDATA
      REAL*8, DIMENSION(IM,JM) :: SNOWE
      REAL*8, DIMENSION(IM,JM) :: TEARTH
      REAL*8, DIMENSION(IM,JM) :: WEARTH
      REAL*8, DIMENSION(IM,JM) :: AIEARTH
      REAL*8, DIMENSION(3,IM,JM) :: SNOAGE

      END MODULE GHYCOM

      SUBROUTINE io_earth(kunit,iaction,ioerr)
!@sum  io_earth reads and writes ground data to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite
      USE GHYCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "EARTH01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,SNOWE,TEARTH,WEARTH,AIEARTH
     *       ,SNOAGE
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,SNOWE,TEARTH,WEARTH,AIEARTH,SNOAGE
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_earth

      SUBROUTINE io_soils(kunit,iaction,ioerr)
!@sum  io_soils reads and writes soil arrays to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite
      USE GHYCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "SOILS01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,wbare,wvege,htbare,htvege
     *       ,snowbv
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,wbare,wvege,htbare,htvege,snowbv
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_soils

      SUBROUTINE io_snow(kunit,iaction,ioerr)
!@sum  io_snow reads and writes snow model arrays to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite
      USE GHYCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "SNOW01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ
     *       ,HSN_IJ,FR_SNOW_IJ 
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ
     *       ,HSN_IJ,FR_SNOW_IJ 
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_snow

