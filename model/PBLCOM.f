      MODULE PBLCOM
!@sum  PBLCOM contains the arrays used by the Boundary Layer code
!@auth Greg Hartke/Ye Cheng
!@ver  1.0
      USE E001M12_COM, only : im,jm
      USE SOCPBL, only : n
      IMPLICIT NONE

!@var uabl boundary layer profile for zonal wind
!@var vabl boundary layer profile for meridional wind
!@var tabl boundary layer profile for temperature
!@var qabl boundary layer profile for humidity
      real*8, dimension(n,im,jm,4) :: uabl,vabl,tabl,qabl
!@var eabl boundary layer profile for turbulent KE (calc. on sec. grid)
      real*8, dimension(n,im,jm,4) :: eabl

!@var cmgs drag coefficient (dimensionless surface momentum flux)
!@var chgs Stanton number   (dimensionless surface heat flux)
!@var cqgs Dalton number    (dimensionless surface moisture flux)
      real*8, dimension(im,jm,4) :: cmgs,chgs,cqgs

!@var ipbl flag for whether pbl properties were found at last timestep
      integer, dimension(im,jm,4) :: ipbl

!@var ROUGHL log10(zgs/roughness length), prescribed with zgs=30 m.
      double precision, dimension(im,jm) :: roughl

!@var WSAVG 1  COMPOSITE SURFACE WIND MAGNITUDE (M/S)
!@var TSAVG 2  COMPOSITE SURFACE AIR TEMPERATURE (K)
!@var QSAVG 3  COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
!@var DCLEV 4  LAYER TO WHICH DRY CONVECTION MIXES (1)
!@var USAVG 6  COMPOSITE SURFACE U WIND
!@var VSAVG 7  COMPOSITE SURFACE V WIND
!@var TAUAVG 8  COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
!@var USTAR friction velocity for each ITYPE (sqrt of srfc momentum flux) (m/s)
      double precision, dimension(im,jm) ::
     &     wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
      double precision, dimension(im,jm,4) :: ustar

c      common /bleq/
c     &     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg,ustar
c!@var BLDATA handle for referring to all boundary layer data
c      DOUBLE PRECISION, DIMENSION(IM,JM,12) :: BLDATA
c      EQUIVALENCE(BLDATA(1,1,1),WSAVG(1,1))

C**** model related constants (should really be taken from E001M12_COM)
      integer, parameter ::  iq1=im/4+1,iq2=im/2+1,iq3=3*im/4+1

C**** pressure gradient arrays
      double precision, dimension(im,jm) ::
     &     dpdxr,dpdyr,phi,dpdxr0,dpdyr0

      END MODULE PBLCOM

      SUBROUTINE io_pbl(kunit,iaction,ioerr)
!@sum  io_pbl reads and writes model variables to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "PBL01"

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        WRITE (KUNIT,ERR=10) MODULE_HEADER,UABL,VABL,TABL,QABL,EABL,
     *     CMGS,CHGS,CQGS,IPBL
      CASE (IOREAD:)            ! input from restart file
        READ (KUNIT,ERR=10) HEADER,UABL,VABL,TABL,QABL,EABL,CMGS,CHGS,
     *     CQGS,IPBL
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_pbl

      SUBROUTINE io_bldat(kunit,iaction,ioerr)
!@sum  io_bldat reads and writes boundary layer data to file 
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : ioread,iowrite
      USE PBLCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "BLD01"

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,wsavg,tsavg,qsavg,dclev
     *       ,usavg,vsavg,tauavg,ustar
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,wsavg,tsavg,qsavg,dclev,usavg
     *       ,vsavg,tauavg,ustar
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_bldat

