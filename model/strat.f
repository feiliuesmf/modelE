C**** Temporary file until STRATDYN is working

      SUBROUTINE DUMMY_STRAT
!@sum DUMMY dummy routines for non-stratospheric models
      ENTRY GWDRAG
      ENTRY VDIFF
      ENTRY EPFLUX
      ENTRY EPFLXI
      ENTRY diaga0
      RETURN
      END SUBROUTINE DUMMY_STRAT

      SUBROUTINE io_strat(kunit,iaction,ioerr)
!@sum  io_strat reads and writes strat. model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "STRAT01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,AIRX,LMC
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) HEADER,AIRX,LMC
        IF (HEADER.ne.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_strat
      
