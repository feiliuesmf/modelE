!@sum OCNML contains routines used for Qflux mixed layer,no deep diff.
!@auth G. Schmidt
!@ver  1.0

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean are reasonable
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE OCEAN, only : tocean
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in ocean data
      CALL CHECK3(TOCEAN,3,IM,JM,SUBR,'toc')

      END SUBROUTINE CHECKO

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean reads and writes ocean arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE OCEAN
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCN01"

      MODULE_HEADER(lhead+1:80) = 'R8 Tocn(3,im,jm),MixLD(im,jm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,TOCEAN,Z1O
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,TOCEAN,Z1O
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocean

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates zonal ocean energy for Qflux ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shw,rhow
      USE MODEL_COM, only : im,jm,fim,focean
      USE OCEAN, only : tocean,z1o,z12o
      USE GEOM, only : imaxj
      IMPLICIT NONE
!@var OCEANE zonal ocean energy (J/M^2)
      REAL*8, DIMENSION(JM) :: OCEANE
      INTEGER I,J

      OCEANE=0
      DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
            OCEANE(J)=OCEANE(J)+(TOCEAN(1,I,J)*Z1O(I,J)
     *           +TOCEAN(2,I,J)*(Z12O(I,J)-Z1O(I,J)))*SHW*RHOW
          END IF
        END DO
      END DO
      OCEANE(1) =FIM*OCEANE(1)
      OCEANE(JM)=FIM*OCEANE(JM)
C****
      END SUBROUTINE conserv_OCE

      SUBROUTINE DUMMY_OCN
!@sum  DUMMY necessary entry points for non-dynamic/non-deep oceans
!@auth Gavin Schmidt
!@ver  1.0
      ENTRY ODYNAM
      ENTRY DYNSI
      ENTRY ADVSI
      ENTRY ODIFS
      ENTRY io_ocdiag
      ENTRY init_ODEEP
      ENTRY reset_ODIAG

      RETURN
      END SUBROUTINE DUMMY_OCN
