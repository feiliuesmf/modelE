#include "rundeck_opts.h"

      MODULE GHYCOM
!@sum  GHYCOM contains the areas used by the Ground Hydrology routines
!@auth Frank Abramopolus/Igor Aleinov
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE SLE001, only : ngm,imt,nlsn
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
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

ccc topmodel input data
      REAL*8, DIMENSION(IM,JM) :: TOP_INDEX_IJ

ccc evaporation limits from previous time step
      real*8, dimension(IM,JM) :: evap_max_ij=1., fr_sat_ij=1.,
     &     qg_ij = 0.

#ifdef TRACERS_WATER
!@var TRBARE,TRVEGE tracers in bare and veg. soil fraction (kg/m^2)
      REAL*8, DIMENSION(NTM,  NGM,IM,JM) :: TRBARE
      REAL*8, DIMENSION(NTM,0:NGM,IM,JM) :: TRVEGE
C**** What is the prognostic variable for snow here?
!@var TRSNOWBV tracer amount in snow over bare and veg. soil (kg/m^2)
      REAL*8, DIMENSION(NTM,2,IM,JM) :: TRSNOWBV
!@var TRSN_IJ tracer amount in snow on earth (kg)  ???? Is this correct?
      REAL*8, DIMENSION(NTM,NLSN,2,IM,JM) :: TRSN_IJ
#endif

      END MODULE GHYCOM

      SUBROUTINE io_earth(kunit,iaction,ioerr)
!@sum  io_earth reads and writes ground data to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE GHYCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "EARTH01"

      MODULE_HEADER(lhead+1:80) =
     *   'R8 dim(ijm) : SNOWe,Te,WTRe,ICEe, SNOage(3,.),evmax,fsat,gq'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,SNOWE,TEARTH,WEARTH,AIEARTH
     *       ,SNOAGE,evap_max_ij,fr_sat_ij,qg_ij
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,SNOWE,TEARTH,WEARTH,AIEARTH
     &       ,SNOAGE,evap_max_ij,fr_sat_ij,qg_ij
        IF (HEADER(1:lhead).NE.MODULE_HEADER(1:lhead)) THEN
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
      USE MODEL_COM, only : ioread,iowrite,lhead,irerun,irsfic
      USE GHYCOM
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SOILS01"
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSOILS01"

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a21,i3,a1,i2,a9,i3,a1,i2,a11,i3,a2)')
     *     'R8 dim(im,jm) TRBARE(',NTM,',',NGM,'),TRVEGE(',NTM,',',NGM+1
     *     ,'),TRSNOWBV(',ntm,'2)'
#endif

      write(MODULE_HEADER(lhead+1:80),'(a6,i1,a11,i1,a)') 'R8 Wb(',
     *   ngm,',ijm), dim(',ngm+1,',ijm):Wv,HTb,HTv, SNWbv(2,ijm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,wbare,wvege,htbare,htvege
     *       ,snowbv
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TRBARE,TRVEGE,TRSNOWBV
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,wbare,wvege,htbare,htvege,snowbv
        IF (HEADER(1:lhead).NE.MODULE_HEADER(1:lhead)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (IRERUN,IOREAD)    ! only need tracers from reruns/restarts
          READ (kunit,err=10) TRHEADER,TRBARE,TRVEGE,TRSNOWBV
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",TRHEADER
     *           ,TRMODULE_HEADER
            GO TO 10
          END IF
        END SELECT
#endif
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_soils

      SUBROUTINE io_snow(kunit,iaction,ioerr)
!@sum  io_snow reads and writes snow model arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irerun,irsfic
      USE GHYCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SNOW01"
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSNOW01"

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a1,i3,a)')'R8 dim(',NTM,',',NLSN,',2,IM,JM):TRSNW'
#endif

      write (MODULE_HEADER(lhead+1:80),'(a29,I1,a)') 'I dim(2,ijm):'//
     *  'Nsn,Isn, R8 dim(',NLSN,',2,ijm):dz,w,ht, Fsn(2,ijm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ
     *       ,HSN_IJ,FR_SNOW_IJ
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TRSN_IJ
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,NSN_IJ,ISN_IJ,DZSN_IJ,WSN_IJ
     *       ,HSN_IJ,FR_SNOW_IJ
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (IRERUN,IOREAD)    ! only need tracers from reruns/restarts
          READ (kunit,err=10) TRHEADER,TRSN_IJ
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",TRHEADER
     *           ,TRMODULE_HEADER
            GO TO 10
          END IF
        END SELECT
#endif
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_snow

