      MODULE CLD01_COM_E001
!@sum  CLD01_COM_E001 model variables for moist convction and
!@sum          large-scale condensation
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE E001M12_COM, only : IM,JM,LM
C**** Note that we USE the NAMELIST set-able parameters from CLD01 only
C**** to be able to pass them to INPUT and io_clouds
      USE CLD01, only : U00wtr,U00ice,LMCM

      IMPLICIT NONE
      SAVE
!@var PREC precipitation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PREC
!@var TPREC temperature of preciptiation (C)
      REAL*8, DIMENSION(IM,JM) :: TPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PRECSS
!@var TTOLD,QTOLD previous potential temperature, humidity
      REAL*8, DIMENSION(LM,IM,JM) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
      REAL*8, DIMENSION(LM,IM,JM) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
      REAL*8, DIMENSION(LM,IM,JM) :: RHSAV
!@var CLDSAV previous cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDSAV
!@var PBLTOP height of PBL (m) (NOT USED)
      REAL*8, DIMENSION(IM,JM) :: PBLTOP

C**** variables saved for radiation calculations
!@var TAUSS optical depth from super-saturated clouds
      REAL*8, DIMENSION(LM,IM,JM) :: TAUSS
!@var TAUMC optical depth from moist-convective clouds
      REAL*8, DIMENSION(LM,IM,JM) :: TAUMC
!@var CLDSS super-saturated cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDMC
!@var CSIZMC,CSIZSS mc,ss effective cloud droplet radius (microns)
      REAL*8, DIMENSION(LM,IM,JM) :: CSIZMC,CSIZSS

      END MODULE CLD01_COM_E001

      SUBROUTINE io_clouds(kunit,iaction,ioerr)
!@sum  io_clouds reads and writes cloud arrays to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite,irestart,irsfic,irerun
      USE CLD01_COM_E001
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "CLD01"
!@var UDUM1,UDUM2,LDUM dummy variables
      REAL*8 UDUM1,UDUM2
      INTEGER LDUM1

      SELECT CASE (IACTION)
      CASE (IOWRITE)           ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,U00wtr,U00ice,LMCM,TTOLD,
     *     QTOLD,SVLHX,RHSAV,CLDSAV
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE  (IACTION)
        CASE (IRESTART,IRERUN) ! input for restart, rerun or extension
          READ (kunit,err=10) HEADER,U00wtr,U00ice,LMCM,TTOLD,QTOLD
     *         ,SVLHX,RHSAV,CLDSAV
        CASE (IRSFIC)           ! start from old restart file
          READ (kunit,err=10) HEADER,Udum1,Udum2,Ldum1,TTOLD,QTOLD
     *         ,SVLHX,RHSAV,CLDSAV
        END SELECT
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_clouds


