#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

      MODULE GHY_COM
!@sum  GHY_COM contains the areas used by the Ground Hydrology routines
!@auth Frank Abramopolus/Igor Aleinov
!@ver  1.0
      USE MODEL_COM, only : im,jm
!!!      USE SLE001, only : ngm,imt,nlsn
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif


      IMPLICIT NONE
      SAVE

ccc CONSTANTS

ccc dimensions of the GHY arrays
!@var ngm number of soil layers
!@var imt number of soil textures
!@var nlsn max number of snow layers
      integer, parameter, public :: ngm=6, imt=5, nlsn=3
!@dbparam WSN_MAX snow amount limit (m, water equivalent)
      real*8, public :: WSN_MAX = 2.d0

!@var LS_NFRAC number of land surface fractions
      integer, parameter, public :: LS_NFRAC=3

!@var shc_soil_texture specific heat capacity of soil texture (J/K/M^3)
      real*8, parameter,public :: shc_soil_texture(imt)
     &     = (/2d6,2d6,2d6,2.5d6,2.4d6/)


ccc variable earth fraction
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FEARTH

ccc bare/veg not in merged array because WBARE does not contain
ccc 0 index for legacy reasons
      !REAL*8, POINTER, DIMENSION(:,:,:) :: WBARE
      !REAL*8, POINTER, DIMENSION(:,:,:) :: WVEGE
      !REAL*8, POINTER, DIMENSION(:,:,:) :: HTBARE
      !REAL*8, POINTER, DIMENSION(:,:,:) :: HTVEGE
      REAL*8, ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: W_IJ
      REAL*8, ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: HT_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOWBV

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: DZ_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: QK_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:)     :: SL_IJ

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

!@var tsns_ij surface temperature corresponding to sensible heat flux (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tsns_ij


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
!!!@var TR_WBARE tracers in bare soil fraction (kg/m^2)
!!!@var TR_WVEGE tracers in vegetated soil fraction (kg/m^2)
!@var TR_W_IJ water tracers in soil (kg/m^2)
!@var TR_WSN_IJ tracer amount in snow (multiplied by fr_snow) (kg/m^2)
      !REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TR_WBARE
      !REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TR_WVEGE
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TR_W_IJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TR_WSN_IJ
ccc TRSNOWBV is not used
!@var TRSNOWBV tracer amount in snow over bare and veg. soil (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSNOWBV0
!@var ntixw index array for tracers (shared by OMP threads)
      integer ntixw(ntm)
#endif
#ifdef USE_ENT
ccc stuff that got back from VEG_COM, maybe should be relocated to Ent
!@var Cint Internal foliage CO2 concentration (mol/m3)
      real*8, ALLOCATABLE, dimension(:,:) :: Ci_ij
!@var Qfol Foliage surface mixing ratio (kg/kg)
      real*8, ALLOCATABLE, dimension(:,:) :: Qf_ij
!@var cnc_ij canopy conductance
      real*8, ALLOCATABLE, dimension(:,:) :: cnc_ij

!@var aalbveg vegetation albedo, eventually should be moved to a
!@+   better place DO WE NEED THIS ???
      real*8, ALLOCATABLE, dimension(:,:) :: aalbveg
#endif
      END MODULE GHY_COM

      SUBROUTINE ALLOC_GHY_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE GHY_COM
      USE DOMAIN_DECOMP, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(     FEARTH(IM,J_0H:J_1H),
     *         STAT=IER)

cddd      ALLOCATE(      WBARE(  NGM,IM,J_0H:J_1H),
cddd     *               WVEGE(0:NGM,IM,J_0H:J_1H),
cddd     *              HTBARE(0:NGM,IM,J_0H:J_1H),
cddd     *              HTVEGE(0:NGM,IM,J_0H:J_1H),
cddd     *              SNOWBV(    2,IM,J_0H:J_1H),
cddd     *         STAT=IER)

      ALLOCATE(      W_IJ(0:NGM,LS_NFRAC,IM,J_0H:J_1H),
     *              HT_IJ(0:NGM,LS_NFRAC,IM,J_0H:J_1H),
     *              SNOWBV(     LS_NFRAC,IM,J_0H:J_1H),
     *         STAT=IER)

      !WBARE => W_IJ(1:,1,1:,1:)
      !WVEGE => W_IJ(1:,1,1:,1:)
      !HTBARE => HT_IJ(1:,1,1:,1:)
      !HTVEGE => HT_IJ(1:,1,1:,1:)

      ALLOCATE(     DZ_IJ(IM,J_0H:J_1H,NGM),
     *               Q_IJ(IM,J_0H:J_1H,IMT,NGM),
     *              QK_IJ(IM,J_0H:J_1H,IMT,NGM),
     *              SL_IJ(IM,J_0H:J_1H),
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

      ALLOCATE(      TOP_INDEX_IJ(IM,J_0H:J_1H),
     *                 top_dev_ij(IM,J_0H:J_1H),
     *                evap_max_ij(IM,J_0H:J_1H),
     *                  fr_sat_ij(IM,J_0H:J_1H),
     *                      qg_ij(IM,J_0H:J_1H),
     *           STAT=IER)
      ALLOCATE(     tsns_ij(IM,J_0H:J_1H),
     *         STAT=IER)

#ifdef USE_ENT
      ALLOCATE(     aalbveg(im,J_0H:J_1H),
     *              Ci_ij(im,J_0H:J_1H),
     *              Qf_ij(im,J_0H:J_1H),
     *              cnc_ij(im,J_0H:J_1H),
     *           STAT=IER)
#endif

C**** Initialize evaporation limits
      evap_max_ij(:,J_0H:J_1H)=1.
      fr_sat_ij(:,J_0H:J_1H)=1.
      qg_ij(:,J_0H:J_1H)=0.

ccc init snow arrays to prevent addressing uninitialized vars
      nsn_ij    (:,:,J_0H:J_1H)=0
      fr_snow_ij(:,:,J_0H:J_1H)=0.
      dzsn_ij (:,:,:,J_0H:J_1H)=0.
      hsn_ij  (:,:,:,J_0H:J_1H)=0.
      wsn_ij  (:,:,:,J_0H:J_1H)=0.

#ifdef TRACERS_WATER_OLD
      ALLOCATE(     TRBARE(NTM,  NGM,IM,J_0H:J_1H),
     *              TRVEGE(NTM,0:NGM,IM,J_0H:J_1H),
     *            TRSNOWBV(NTM,    2,IM,J_0H:J_1H),
     *         STAT=IER)
#endif
#ifdef TRACERS_WATER
cddd      ALLOCATE(     TR_WBARE(NTM,  NGM,IM,J_0H:J_1H),
cddd     *              TR_WVEGE(NTM,0:NGM,IM,J_0H:J_1H),
      ALLOCATE(      TR_W_IJ(NTM,0:NGM,LS_NFRAC,IM,J_0H:J_1H),
     *             TR_WSN_IJ(NTM,NLSN,        2,IM,J_0H:J_1H),
     *             TRSNOWBV0(NTM,             2,IM,J_0H:J_1H),
     *         STAT=IER)
C**** Initialize to zero
      !TR_WBARE(:,:,:,J_0H:J_1H)=0.d0
      !TR_WVEGE(:,:,:,J_0H:J_1H)=0.d0
      TR_W_IJ(:,:,:,:,J_0H:J_1H)=0.d0
      TR_WSN_IJ(:,:,:,:,J_0H:J_1H)=0.d0
#endif

      END SUBROUTINE ALLOC_GHY_COM

      SUBROUTINE io_earth(kunit,iaction,ioerr)
!@sum  io_earth reads and writes ground data to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE GHY_COM
      USE DOMAIN_DECOMP, only : GRID, GET, AM_I_ROOT
      USE DOMAIN_DECOMP, only : PACK_DATA, PACK_COLUMN
      USE DOMAIN_DECOMP, only : UNPACK_DATA, UNPACK_COLUMN
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "EARTH02"

      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SNOWE_glob, TEARTH_glob,
     &     WEARTH_glob, AIEARTH_GLOB,evap_max_ij_glob, fr_sat_ij_glob,
     &     qg_ij_glob
      REAL*8, ALLOCATABLE :: SNOAGE_glob(:,:,:)
      REAL*8, ALLOCATABLE :: tsns_ij_glob(:,:)

      call allocate_me

      MODULE_HEADER(lhead+1:80) =
     *   'R8 dim(ijm) : SNOWe,Te,WTRe,ICEe, SNOage(3,.),evmax,fsat,gq'
!     *     //',fe'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_DATA(grid, SNOWE       , SNOWE_glob)
        CALL PACK_DATA(grid, TEARTH      , TEARTH_glob)
        CALL PACK_DATA(grid, WEARTH      , WEARTH_glob)
        CALL PACK_DATA(grid, AIEARTH     , AIEARTH_GLOB)
        CALL PACK_DATA(grid, evap_max_ij , evap_max_ij_glob)
        CALL PACK_DATA(grid, fr_sat_ij   , fr_sat_ij_glob)
        CALL PACK_DATA(grid, qg_ij       , qg_ij_glob)
        CALL PACK_COLUMN(grid, SNOAGE    , SNOAGE_glob)
        CALL PACK_DATA(grid, tsns_ij     , tsns_ij_glob)
        IF (AM_I_ROOT())
     *     WRITE (kunit,err=10) MODULE_HEADER,SNOWE_glob,TEARTH_glob
     *       ,WEARTH_glob,AIEARTH_glob
     *       ,SNOAGE_glob,evap_max_ij_glob,fr_sat_ij_glob,qg_ij_glob
     &       ,tsns_ij_glob
      CASE (IOREAD:)            ! input from restart file
cgsfc        READ (kunit,err=10) HEADER,SNOWE,TEARTH,WEARTH,AIEARTH
cgsfc     &       ,SNOAGE,evap_max_ij,fr_sat_ij,qg_ij
        if (AM_I_ROOT()) then
          READ(kunit,err=10) HEADER
          BACKSPACE kunit
          if (HEADER(1:lhead) == "EARTH01" ) then ! hack to read old format
              READ (kunit,err=10) HEADER,SNOWE_glob,TEARTH_glob
     &         ,WEARTH_glob ,AIEARTH_glob,SNOAGE_glob
     &         ,evap_max_ij_glob,fr_sat_ij_glob ,qg_ij_glob
            tsns_ij_glob(:,:) = TEARTH_glob
          else if (HEADER(1:lhead) == MODULE_HEADER(1:lhead)) then
            READ(kunit,err=10) HEADER,SNOWE_glob,TEARTH_glob
     &           ,WEARTH_glob ,AIEARTH_glob,SNOAGE_glob
     &           ,evap_max_ij_glob,fr_sat_ij_glob ,qg_ij_glob
     &           ,tsns_ij_glob
          else
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          end if
        endif ! I_AM_ROOT

        CALL UNPACK_DATA(grid, SNOWE_glob       , SNOWE      )
        CALL UNPACK_DATA(grid, TEARTH_glob      , TEARTH     )
        CALL UNPACK_DATA(grid, WEARTH_glob      , WEARTH     )
        CALL UNPACK_DATA(grid, AIEARTH_glob     , AIEARTH    )
        CALL UNPACK_DATA(grid, evap_max_ij_glob , evap_max_ij)
        CALL UNPACK_DATA(grid, fr_sat_ij_glob   , fr_sat_ij  )
        CALL UNPACK_DATA(grid, qg_ij_glob       , qg_ij      )
        CALL UNPACK_COLUMN(grid, SNOAGE_glob    , SNOAGE     )
        CALL UNPACK_DATA(grid, tsns_ij_glob     , tsns_ij      )

      END SELECT

      call deallocate_me
      RETURN
 10   IOERR=1
      call deallocate_me
      RETURN

      contains
      subroutine allocate_me
      if (AM_I_ROOT()) then
        ALLOCATE( SNOAGE_glob(3,IM,JM),
     &       tsns_ij_glob(IM,JM),
     &       SNOWE_glob(IM,JM),
     &       TEARTH_glob(IM,JM),
     &       WEARTH_glob(IM,JM),
     &       AIEARTH_GLOB(IM,JM),
     &       evap_max_ij_glob(IM,JM),
     &       fr_sat_ij_glob(IM,JM),
     &       qg_ij_glob(IM,JM) )
      endif
      end subroutine allocate_me
      subroutine deallocate_me
      if (AM_I_ROOT()) then
        DEALLOCATE( SNOAGE_glob,
     &       tsns_ij_glob,
     &       SNOWE_glob,
     &       TEARTH_glob,
     &       WEARTH_glob,
     &       AIEARTH_GLOB,
     &       evap_max_ij_glob,
     &       fr_sat_ij_glob,
     &       qg_ij_glob )
      endif
      end subroutine deallocate_me

      END SUBROUTINE io_earth


      SUBROUTINE io_soils(kunit,iaction,ioerr)
!@sum  io_soils reads and writes soil arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
      USE DOMAIN_DECOMP, ONLY: GRID, GET, CHECKSUM_COLUMN
      USE DOMAIN_DECOMP, ONLY:  PACK_COLUMN, AM_I_ROOT
      USE DOMAIN_DECOMP, ONLY: PACK_BLOCK, UNPACK_BLOCK
      USE DOMAIN_DECOMP, ONLY: UNPACK_COLUMN
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      USE GHY_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
      REAL*8, ALLOCATABLE ::  SNOWBV_GLOB(:,:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: W_GLOB,HT_GLOB
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SOILS03"
      INTEGER :: J_0H, J_1H

#ifdef TRACERS_WATER
      integer m
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSOILS03"
      REAL*8, ALLOCATABLE :: TRSNOWBV0_GLOB(:,:,:,:)
      REAL*8, ALLOCATABLE, TARGET :: TR_W_GLOB  (:,:,:,:,:)
      REAL*8, POINTER :: ptr4(:,:,:,:)
      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a21,i3,a1,i2,a1,i2,a11,i3,a2)')
     *     'R8 dim(im,jm) TR_W(',NTM,',',NGM,',',LS_NFRAC
     *     ,'),TRSNOWBV(',ntm,'2)'
#endif

      call allocate_me

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      write(MODULE_HEADER(lhead+1:80),'(a7,i1,a1,i1,a23)')
     *   'R8 dim(',ngm+1,',',LS_NFRAC,',ijm):W,HT, SNWbv(2,ijm)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_BLOCK(grid, W_IJ(0:NGM,:,1:IM,J_0H:J_1H) , W_GLOB)
        CALL PACK_BLOCK(grid, HT_IJ(0:NGM,:,1:IM,J_0H:J_1H),HT_GLOB)
        CALL PACK_COLUMN(grid, SNOWBV(1:2,1:IM,J_0H:J_1H), SNOWBV_GLOB)
#ifdef TRACERS_WATER
        do m=1,LS_NFRAC
          if (AM_I_ROOT()) ptr4 => TR_W_GLOB(1:NTM,0:NGM,m,1:IM,1:JM)
          CALL PACK_BLOCK(grid, TR_W_IJ(1:NTM,0:NGM,m,1:IM,J_0H:J_1H),
     &         ptr4 )
        enddo
        CALL PACK_BLOCK(grid, TRSNOWBV0, TRSNOWBV0_GLOB)
#endif
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,w_glob,
     *       ht_glob,snowbv_glob
#ifdef TRACERS_WATER
          WRITE (kunit,err=10) TRMODULE_HEADER,TR_W_GLOB
     &       ,TRSNOWBV0_GLOB
#endif
        END IF
      CASE (IOREAD:)            ! input from restart file
        if (AM_I_ROOT()) then
          READ(kunit,err=10) HEADER
          BACKSPACE kunit
          if (HEADER(1:lhead) == "SOILS02" ) then ! hack to read old format
            w_glob = 0.d0; ht_glob = 0.d0
            READ(kunit,err=10) HEADER, w_glob(1:NGM,1,:,:),
     &           w_glob(0:NGM,2,:,:), ht_glob(0:NGM,1,:,:),
     &           ht_glob(0:NGM,2,:,:), snowbv_glob
          else if (HEADER(1:lhead) == MODULE_HEADER(1:lhead)) then
            READ(kunit,err=10) HEADER,w_glob,ht_glob,snowbv_glob
          else
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          end if
        end if    !...am_i_root

        call unpack_block(grid, w_glob,
     &       w_ij(0:NGM,:,1:IM,J_0H:J_1H))
         call unpack_block(grid, ht_glob,
     &       ht_ij(0:NGM,:,1:IM,J_0H:J_1H))
        call unpack_column(grid, snowbv_glob,
     &       snowbv(1:2,1:IM,J_0H:J_1H))

#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)  ! reruns/restarts
          if (AM_I_ROOT()) then
            READ (kunit,err=10) TRHEADER, TR_W_GLOB
     &         ,TRSNOWBV0_GLOB
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
          end if
          do m=1,LS_NFRAC
            if (AM_I_ROOT()) ptr4 => TR_W_GLOB(1:NTM,0:NGM,m,1:IM,1:JM)
            CALL UNPACK_BLOCK(grid,ptr4,
     &           TR_W_IJ(1:NTM,0:NGM,m,1:IM,J_0H:J_1H) )
          enddo
          CALL UNPACK_BLOCK(grid,TRSNOWBV0_GLOB,TRSNOWBV0)

        END SELECT
#endif
      END SELECT

      call deallocate_me
      RETURN
 10   IOERR=1
      call deallocate_me
      RETURN

      contains
      subroutine allocate_me
      if (AM_I_ROOT()) then
        ALLOCATE( SNOWBV_GLOB(2,IM,JM),
     &       W_GLOB(0:NGM,LS_NFRAC,IM,JM),
     &       HT_GLOB(0:NGM,LS_NFRAC,IM,JM) )
#ifdef TRACERS_WATER
        ALLOCATE( TRSNOWBV0_GLOB(NTM,2,IM,JM),
     &       TR_W_GLOB  (NTM,0:NGM,LS_NFRAC,IM,JM) )
#endif
      endif
      end subroutine allocate_me
      subroutine deallocate_me
      if (AM_I_ROOT()) then
        DEALLOCATE(SNOWBV_GLOB,
     &       W_GLOB,
     &       HT_GLOB )
#ifdef TRACERS_WATER
        DEALLOCATE( TRSNOWBV0_GLOB,
     &       TR_W_GLOB )
#endif
      endif
      end subroutine deallocate_me

      END SUBROUTINE io_soils


      SUBROUTINE io_snow(kunit,iaction,ioerr)
!@sum  io_snow reads and writes snow model arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
      USE DOMAIN_DECOMP, only : grid, AM_I_ROOT
      USE DOMAIN_DECOMP, only : PACK_BLOCK  , PACK_COLUMN
      USE DOMAIN_DECOMP, only : UNPACK_BLOCK, UNPACK_COLUMN
      USE GHY_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "SNOW01"
      INTEGER, ALLOCATABLE ::  NSN_IJ_GLOB(:,:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: DZSN_IJ_GLOB,
     &     WSN_IJ_GLOB, HSN_IJ_GLOB
      REAL*8, ALLOCATABLE :: FR_SNOW_IJ_GLOB(:,:,:)
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRSNOW01"
      REAL*8, ALLOCATABLE, TARGET ::  TR_WSN_IJ_GLOB(:,:,:,:,:)
      REAL*8, POINTER :: ptr4(:,:,:,:)

      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a1,i3,a)')'R8 dim(',NTM,',',NLSN,',2,IM,JM):TRSNW'
#endif

      write (MODULE_HEADER(lhead+1:80),'(a29,I1,a)') 'I dim(2,ijm):'//
     *  'Nsn, R8 dim(',NLSN,',2,ijm):dz,w,ht, Fsn(2,ijm)'

      call allocate_me

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_BLOCK(grid, DZSN_IJ, DZSN_IJ_GLOB)
        CALL PACK_BLOCK(grid,  WSN_IJ,  WSN_IJ_GLOB)
        CALL PACK_BLOCK(grid,  HSN_IJ,  HSN_IJ_GLOB)
        CALL PACK_COLUMN(grid, NSN_IJ,  NSN_IJ_GLOB)
        CALL PACK_COLUMN(grid, FR_SNOW_IJ,  FR_SNOW_IJ_GLOB)
#ifdef TRACERS_WATER
        if (AM_I_ROOT()) ptr4 => TR_WSN_IJ_GLOB(:,:,1,:,:)
        CALL PACK_BLOCK(grid, TR_WSN_IJ(     :,:,1,:,:), ptr4 )
        if (AM_I_ROOT()) ptr4 => TR_WSN_IJ_GLOB(:,:,2,:,:)
        CALL PACK_BLOCK(grid, TR_WSN_IJ(     :,:,2,:,:), ptr4 )
#endif
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER, NSN_IJ_glob, DZSN_IJ_glob
     *         ,WSN_IJ_glob, HSN_IJ_glob, FR_SNOW_IJ_glob
#ifdef TRACERS_WATER
          WRITE (kunit,err=10) TRMODULE_HEADER,TR_WSN_IJ_glob
#endif
        END IF
      CASE (IOREAD:)            ! input from restart file
        if (AM_I_ROOT()) then
          READ (kunit,err=10) HEADER,NSN_IJ_glob, DZSN_IJ_glob
     *           ,WSN_IJ_glob, HSN_IJ_glob, FR_SNOW_IJ_glob
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if
        CALL UNPACK_BLOCK(grid, DZSN_IJ_GLOB, DZSN_IJ)
        CALL UNPACK_BLOCK(grid,  WSN_IJ_GLOB,  WSN_IJ)
        CALL UNPACK_BLOCK(grid,  HSN_IJ_GLOB,  HSN_IJ)
        CALL UNPACK_COLUMN(grid, NSN_IJ_GLOB,  NSN_IJ)
        CALL UNPACK_COLUMN(grid, FR_SNOW_IJ_GLOB,  FR_SNOW_IJ)
#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO) ! reruns/restarts
          if (AM_I_ROOT()) then
            READ (kunit,err=10) TRHEADER,TR_WSN_IJ_GLOB
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
          end if

          if (AM_I_ROOT()) ptr4 => TR_WSN_IJ_GLOB(:,:,1,:,:)
          CALL UNPACK_BLOCK(grid, ptr4, TR_WSN_IJ(:,:,1,:,:))
          if (AM_I_ROOT()) ptr4 => TR_WSN_IJ_GLOB(:,:,2,:,:)
          CALL UNPACK_BLOCK(grid, ptr4, TR_WSN_IJ(:,:,2,:,:))
        END SELECT
#endif
      END SELECT

      call deallocate_me
      RETURN
 10   IOERR=1
      call deallocate_me
      RETURN      

      contains
      subroutine allocate_me
      if (AM_I_ROOT()) then
        ALLOCATE( NSN_IJ_GLOB(2,IM,JM),
     &       DZSN_IJ_GLOB(NLSN,2,IM,JM),
     &       WSN_IJ_GLOB(NLSN,2,IM,JM),
     &       HSN_IJ_GLOB(NLSN,2,IM,JM),
     &       FR_SNOW_IJ_GLOB(2,IM,JM) )
#ifdef TRACERS_WATER
        ALLOCATE( TR_WSN_IJ_GLOB(NTM,NLSN,2,IM,JM) )
#endif
      endif
      end subroutine allocate_me
      subroutine deallocate_me
      if (AM_I_ROOT()) then
        DEALLOCATE( NSN_IJ_GLOB,
     &       DZSN_IJ_GLOB,
     &       WSN_IJ_GLOB,
     &       HSN_IJ_GLOB,
     &       FR_SNOW_IJ_GLOB )
#ifdef TRACERS_WATER
        DEALLOCATE( TR_WSN_IJ_GLOB )
#endif
      endif
      end subroutine deallocate_me

      END SUBROUTINE io_snow

#ifdef USE_ENT
      subroutine io_veg_related(kunit,iaction,ioerr)
!@sum  reads and writes data needed to drive vegetation module
!@auth I. Aleinov
!@ver  1.0
      use model_com, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
      use model_com, only : im,jm
      use domain_decomp, only : grid, am_i_root
      use domain_decomp, only : pack_data, unpack_data
      use ghy_com, only : Ci_ij, Qf_ij, cnc_ij
      use param
      implicit none

      integer kunit   !@var kunit unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr
!@var header character string label for individual records
      character*80 :: header, module_header = "vegetation01"
!@var Ci_ij_glob work array for parallel_io
!@var Qf_ij_glob work array for parallel_io
!@var cnc_ij_glob work array for parallel_io
      real*8, dimension(im,jm) :: Ci_ij_glob, Qf_ij_glob, cnc_ij_glob
      integer :: force_init_ent=0

      call allocate_me
!!! hack
      call sync_param( "init_ent", force_init_ent)

      write(module_header(lhead+1:80),'(a)') 'Ci_ij,Qf_ij,cnc_ij'

      select case (iaction)
      case (:iowrite)            ! output to standard restart file
        call pack_data(grid, Ci_ij, Ci_ij_glob)
        call pack_data(grid, Qf_ij, Qf_ij_glob)
        call pack_data(grid, cnc_ij, cnc_ij_glob)
        if (am_i_root())
     &    write (kunit,err=10) module_header,Ci_ij_glob,Qf_ij_glob,
     &                         cnc_ij_glob
      case (ioread:)            ! input from restart file
       if ( force_init_ent .ne. 1 ) then
        if ( AM_I_ROOT() ) then
          read(kunit,err=10) header, Ci_ij_glob, Qf_ij_glob, cnc_ij_glob
          if (header(1:lhead).ne.module_header(1:lhead)) then
            print*,"discrepancy in module version ",header,module_header
            go to 10
          end if
        end if
        call unpack_data(grid, Ci_ij_glob,   Ci_ij)
        call unpack_data(grid, Qf_ij_glob,   Qf_ij)
        call unpack_data(grid, cnc_ij_glob, cnc_ij)
       endif
      end select

      call deallocate_me
      return
 10   ioerr=1
      call deallocate_me
      return

      contains
      subroutine allocate_me
      if (AM_I_ROOT()) then
        ALLOCATE( Ci_ij_glob(im,jm),
     &       Qf_ij_glob(im,jm),
     &       cnc_ij_glob(im,jm) )
      endif
      end subroutine allocate_me
      subroutine deallocate_me
      if (AM_I_ROOT()) then
        DEALLOCATE( Ci_ij_glob,
     &       Qf_ij_glob,
     &       cnc_ij_glob )
      endif
      end subroutine deallocate_me

      end subroutine io_veg_related
#endif
