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
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WBARE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WVEGE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: HTBARE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: HTVEGE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOWBV

C GISS-ESMF EXCEPTIONAL CASE
C-BMP Should be allocatable, but not sure about common block
      REAL*8, DIMENSION(IM,JM,NGM) :: DZ_IJ
      REAL*8, DIMENSION(IM,JM,IMT,NGM) :: Q_IJ
      REAL*8, DIMENSION(IM,JM,IMT,NGM) :: QK_IJ
      REAL*8, DIMENSION(IM,JM) :: SL_IJ
C this common only for purpose of reading its contents from a
C file opened in fortran unformatted sequential access mode
C containing its contents in a contiguous real*4 block
      COMMON/SDATA/ DZ_IJ,Q_IJ,QK_IJ,SL_IJ

ccc the following arrays contain prognostic variables for the snow model
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:)     :: NSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: DZSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: WSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: HSN_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)      :: FR_SNOW_IJ
ccc FR_SNOW_RAD_IJ is snow fraction for albedo computations
ccc actually it should be the same as FR_SNOW_IJ but currently the snow
ccc model can't handle fractional cover for thick snow (will fix later)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)      :: FR_SNOW_RAD_IJ
C**** Canopy temperature (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)      :: CANOPY_TEMP_IJ
C**** replacements for GDATA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SNOWE
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TEARTH
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: WEARTH
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: AIEARTH
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOAGE

!@var GDEEP keeps average (2:n) values of temperature, water and ice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GDEEP

!@dbparam snoage_def determines how snowage is calculated:
!@+       = 0     independent of temperature
!@+       = 1     only when max daily local temp. over type > 0
      integer :: snoage_def = 0

ccc topmodel input data and standard deviation of the elevation
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TOP_INDEX_IJ, top_dev_ij

ccc evaporation limits from previous time step
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: evap_max_ij, fr_sat_ij,
     &     qg_ij

#ifdef TRACERS_WATER_OLD
!@var TRBARE,TRVEGE tracers in bare and veg. soil fraction (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRBARE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRVEGE
C**** What is the prognostic variable for snow here?
!@var TRSNOWBV tracer amount in snow over bare and veg. soil (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSNOWBV
#endif
#ifdef TRACERS_WATER
ccc new tracers
      !integer, parameter :: NTM = 3
!@var TR_WBARE tracers in bare soil fraction (kg/m^2)
!@var TR_WVEGE tracers in vegetated soil fraction (kg/m^2)
!@var TR_WSN_IJ tracer amount in snow (multiplied by fr_snow) (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TR_WBARE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TR_WVEGE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TR_WSN_IJ
ccc TRSNOWBV is not used
!@var TRSNOWBV tracer amount in snow over bare and veg. soil (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSNOWBV0
#endif

      END MODULE GHYCOM

      SUBROUTINE ALLOC_GHY_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE GHYCOM
      USE DOMAIN_DECOMP, ONLY : DYN_GRID, GET
      IMPLICIT NONE
      TYPE (DYN_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(     WBARE(  NGM,IM,J_0H:J_1H),
     *               WVEGE(0:NGM,IM,J_0H:J_1H),
     *              HTBARE(0:NGM,IM,J_0H:J_1H),
     *              HTVEGE(0:NGM,IM,J_0H:J_1H),
     *              SNOWBV(    2,IM,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(      NSN_IJ(     2,IM,J_0H:J_1H),
     *              DZSN_IJ(NLSN,2,IM,J_0H:J_1H),
     *               WSN_IJ(NLSN,2,IM,J_0H:J_1H),
     *               HSN_IJ(NLSN,2,IM,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(         FR_SNOW_IJ(2,IM,J_0H:J_1H),
     *              FR_SNOW_RAD_IJ(2,IM,J_0H:J_1H),
     *              CANOPY_TEMP_IJ(  IM,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(      SNOWE(  IM,J_0H:J_1H),
     *              TEARTH(  IM,J_0H:J_1H),
     *              WEARTH(  IM,J_0H:J_1H),
     *             AIEARTH(  IM,J_0H:J_1H),
     *              SNOAGE(3,IM,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(     GDEEP(IM,J_0H:J_1H,3),
     *         STAT=IER)

        ALLOCATE(    TOP_INDEX_IJ(IM,J_0H:J_1H),
     *                 top_dev_ij(IM,J_0H:J_1H),
     *                evap_max_ij(IM,J_0H:J_1H),
     *                  fr_sat_ij(IM,J_0H:J_1H),
     *                      qg_ij(IM,J_0H:J_1H),
     *           STAT=IER)

C**** Initialize evaporation limits
      evap_max_ij(:,J_0H:J_1H)=1.
      fr_sat_ij(:,J_0H:J_1H)=1.
      qg_ij(:,J_0H:J_1H)=0.

#ifdef TRACERS_WATER_OLD
      ALLOCATE(     TRBARE(NTM,  NGM,IM,J_0H:J_1H),
     *              TRVEGE(NTM,0:NGM,IM,J_0H:J_1H),
     *            TRSNOWBV(NTM,    2,IM,J_0H:J_1H),
     *         STAT=IER)
#endif
#ifdef TRACERS_WATER
      ALLOCATE(     TR_WBARE(NTM,  NGM,IM,J_0H:J_1H),
     *              TR_WVEGE(NTM,0:NGM,IM,J_0H:J_1H),
     *             TR_WSN_IJ(NTM,NLSN,2,IM,J_0H:J_1H),
     *             TRSNOWBV0(NTM,2,IM,J_0H:J_1H),
     *         STAT=IER)
C**** Initialize to zero
      TR_WBARE(:,:,:,J_0H:J_1H)=0.d0
      TR_WVEGE(:,:,:,J_0H:J_1H)=0.d0
      TR_WSN_IJ(:,:,:,:,J_0H:J_1H)=0.d0
#endif

      END SUBROUTINE ALLOC_GHY_COM

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
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
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
      USE MODEL_COM, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      USE GHYCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SOILS02"
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSOILS01"

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a21,i3,a1,i2,a9,i3,a1,i2,a11,i3,a2)')
     *     'R8 dim(im,jm) TRBARE(',NTM,',',NGM,'),TRVEGE(',NTM,',',NGM+1
     *     ,'),TRSNOWBV(',ntm,'2)'
#endif

      write(MODULE_HEADER(lhead+1:80),'(a6,i1,a11,i1,a,a)') 'R8 Wb(',
     *   ngm,',ijm), dim(',ngm+1,',ijm):Wv,HTb,HTv, SNWbv(2,ijm),'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,wbare,wvege,htbare,htvege
     *       ,snowbv
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TR_WBARE,TR_WVEGE,TRSNOWBV0
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,wbare,wvege,htbare,htvege,snowbv
        IF (HEADER(1:lhead).NE.MODULE_HEADER(1:lhead)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)  ! reruns/restarts
          READ (kunit,err=10) TRHEADER,TR_WBARE,TR_WVEGE,TRSNOWBV0
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TRHEADER
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
      USE MODEL_COM, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
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
     *  'Nsn, R8 dim(',NLSN,',2,ijm):dz,w,ht, Fsn(2,ijm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,NSN_IJ,DZSN_IJ,WSN_IJ
     *       ,HSN_IJ,FR_SNOW_IJ
#ifdef TRACERS_WATER
        WRITE (kunit,err=10) TRMODULE_HEADER,TR_WSN_IJ
#endif
      CASE (IOREAD:)            ! input from restart file
        READ (kunit,err=10) HEADER,NSN_IJ,DZSN_IJ,WSN_IJ
     *       ,HSN_IJ,FR_SNOW_IJ
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO) ! reruns/restarts
          READ (kunit,err=10) TRHEADER,TR_WSN_IJ
          IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TRHEADER
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

