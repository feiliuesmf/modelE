      MODULE SOMTQ_COM
!@sum  SOMTQ_COM contains the arrays containing second order moments
!@auth Gary Russell
!@ver  1.0
      USE E001M12_COM, only : im,jm,lm

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) ::
     *  TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ
      COMMON/TADV/
     *  TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) ::
     *  QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
      COMMON/QADV/
     *  QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,9) :: TMOM,QMOM
      EQUIVALENCE (TX,TMOM),(QX,QMOM)

      END MODULE SOMTQ_COM

      SUBROUTINE io_somtq(kunit,iaction,ioerr)
!@sum  io_somtq reads and writes second order moments to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
      USE SOMTQ_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "QUS01"

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        WRITE (KUNIT,ERR=10) MODULE_HEADER,TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX
     *     ,TYZ,QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
      CASE (IOREAD:)            ! input from restart file
        READ (KUNIT,ERR=10) HEADER,TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ
     *     ,QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_somtq

