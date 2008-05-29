#include "rundeck_opts.h"
      MODULE RAD_COM
!@sum  RAD_COM Model radiation arrays and parameters
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,lm_req
      USE DOMAIN_DECOMP, only : grid
      USE RADPAR, only : S0,ITRMAX
!@var S0 solar 'constant' needs to be saved between calls to radiation
      IMPLICIT NONE
      SAVE

C**** DEFAULT ORBITAL PARAMETERS FOR EARTH
C**** Note PMIP runs had specified values that do not necesarily
C**** coincide with those used as the default, or the output of ORBPAR.
C****                    OMEGT          OBLIQ        ECCEN
C**** DEFAULT (2000 AD): 282.9          23.44        0.0167
C**** PMIP CONTROL:      282.04         23.446       0.016724
C**** PMIP 6kyr BP:      180.87         24.105       0.018682
C**** PMIP LGM (21k):    294.42         22.949       0.018994
!@param OMEGT_def precession angle (degrees from vernal equinox)
      real*8, parameter :: omegt_def = 282.9d0
!@param OBLIQ_def obliquity angle  (degrees)
      real*8, parameter :: obliq_def = 23.44d0
!@param ECCN_def eccentricity
      real*8, parameter :: eccn_def  = .0167d0
!@var OMEGT,OBLIQ,ECCN actual orbital parameters used
      real*8 OMEGT,OBLIQ,ECCN

C**** Database parameters to control orbital parameter calculation
C**** Note: setting calc_orb_par with paleo_orb_yr=2000 does not produce
C**** exactly the same as the default values.
!@dbparam calc_orb_par = 1 to calc orbital parameters
      integer :: calc_orb_par = 0
!@dbparam paleo_orb_yr is paleo year (BP) for orbital calc
      real*8 :: paleo_orb_yr = -50.  ! (i.e. 2000AD)
!@dbparam calc_orb_par_sp = 1 to directly specify orbital parameters
      integer :: calc_orb_par_sp = 0
!@dbparam paleo_orb_par :: directly specifies orbital parameters
      real*8, dimension(3) :: paleo_orb_par (/ eccn_def, obliq_def,
     *     omegt_def /)

!@var dimrad_sv dimension sum of input fields saved for radia_only runs
      INTEGER, PARAMETER :: dimrad_sv=IM*JM*(7*LM+3*LM_REQ+24)
!@var RQT Radiative equilibrium temperatures above model top
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RQT
!@var Tchg Total temperature change in adjusted forcing runs
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: Tchg
!@var SRHR(0) Solar   raditive net flux into the ground          (W/m^2)
!@var TRHR(0) Thermal raditive net flux into ground(W/O -StB*T^4)(W/m^2)
!@*   Note: -StB*T^4 is MISSING, since T may vary a lot betw. rad. calls
!@var SRHR(1->LM) Solar   raditive heating rate (W/m^2)  (short wave)
!@var TRHR(1->LM) Thermal raditive heating rate (W/m^2)  (long wave)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SRHR,TRHR
!@var TRSURF upward thermal radiation at the surface from rad step W/m2
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRSURF
!@var FSF Solar Forcing over each type (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FSF
!@var FSRDIR Direct beam solar incident at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FSRDIR
!@var SRVISSURF Incident solar direct+diffuse visible at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRVISSURF
!@var SRDN Total incident solar at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRDN  ! saved in rsf
!@var FSRDIF diffuse visible incident solar at surface
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FSRDIF
!@var DIRNIR direct  nir     incident solar at surface
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DIRNIR
!@var DIFNIR diffuse nir     incident solar at surface
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DIFNIR
!@var srnflb_save  Net solar radiation (W/m^2)
!@var trnflb_save  Net thermal radiation (W/m^2)
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: srnflb_save,trnflb_save
!@var ttausv_save  Tracer optical thickness
!@var ttausv_cs_save  Tracer optical thickness clear sky
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: ttausv_save,
     &     ttausv_cs_save

!@var CFRAC Total cloud fraction as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: CFRAC ! saved in rsf
!@var RCLD Total cloud optical depth as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RCLD ! saved in rsf
!@var O3_tracer_save 3D ozone saved elsewhere for use in radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: O3_tracer_save!saved rsf
!@var rad_to_chem save 3D quantities from radiation code for use in
!@+   chemistry (or rest of model). 1=Ozone, 2=aerosol ext, 3=N2O, 4=CH4,
!@+   5=CFC11+CFC12
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: rad_to_chem !saved in rsf
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: rad_to_file
!@var KLIQ Flag indicating dry(0)/wet(1) atmosphere (memory feature)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: KLIQ ! saved in rsf
!@dbparam Ikliq 0,1,-1 initialize kliq as dry,equil,current model state
      INTEGER :: Ikliq = -1  !  get kliq-array from restart file
!@dbparam RHfix const.rel.humidity passed to radiation for aeros. tests
      REAL*8 :: RHfix = -1.  !  pass the current model rel.humidity
!@dbparam dalbsnX global coeff for snow alb change by black carbon depos
      REAL*8 ::  dalbsnX = 0.
!@dbparam albsn_yr year of blk carb depos used for snow alb. reduction
      INTEGER ::  albsn_yr = 1951

!     variables related to aerosol indirect effects:
!     (CDNC=cloud droplet number concentration)
!@dbparam CC_CDNCx scaling factor relating cld cvr change and CDNC change
      REAL*8 :: CC_CDNCX = .0000d0  ! .0036d0
!@dbparam OC_CDNCx scaling factor relating cld opt depth and CDNC change
      REAL*8 :: OD_CDNCX = .0000d0  ! .007d0
!@var pcdnc,vcdnc pressure,vertical profile for cld.cvr change
      real*8, parameter, dimension(7) ::
     * pcdnc=(/984.d0, 964.d0, 934.d0, 884.d0, 810.d0, 710.d0, 550.d0/)
     *,vcdnc=(/ .35d0,  .20d0,  .10d0,  .17d0,  .10d0,  .08d0,   0.d0/)
!@var cdncl = vcdnc interpolated to current vertical resolution
      real*8 cdncl(LM)

!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ1
!@dbparam S0X solar constant multiplication factor
      REAL*8 :: S0X = 1.
!@dbparam S0_yr,S0_day obs.date of solar constant (if 0: time var)
      INTEGER :: S0_yr = 1951 , S0_day = 182
!@dbparam CO2X,... scaling factors for CO2 N2O CH4 CFC11 CFC12 XGHG
      REAL*8 :: CO2X=1.,N2OX=1.,CH4X=1., CFC11X=1.,CFC12X=1.,XGHGX=1.
!@dbparm ref_mult factor to control REFDRY from rundeck
      REAL*8 :: ref_mult = 1.
!@dbparam GHG_yr,GHG_day obs.date of well-mixed GHgases (if 0: time var)
      INTEGER :: GHG_yr = 1951 , GHG_day = 182
!@dbparam Volc_yr,Volc_day obs.date of Volc.Aerosols (if 0: time var)
      INTEGER :: Volc_yr = 1951 , Volc_day = 182
!@dbparam Aero_yr obs.year of troposph.Aerosols (if 0: time var)
      INTEGER :: Aero_yr = 1951    ! always use annual cycle
!@dbparam O3_yr obs.year of Ozone (if 0: time var, -yyyy: use decyyyy,
!@+   janyyyy to find jan1-15,dec16-31,not dec of yyyy-1, jan of yyyy+1)
      INTEGER :: O3_yr = 1951      ! always use annual cycle
!@dbparam crops_yr obs.year of crops (if 0: time var, -1: default)
      INTEGER :: crops_yr = -1
!@dbparam H2OstratX strat_water_vapor, cloud, Ozone scaling factor
      REAL*8 :: H2OstratX = 1. , cldX = 1. , O3X = 1.
!@dbparam H2ObyCH4 if not 0: add CH4 produced H2O into layers 1->LM
      REAL*8 :: H2ObyCH4 = 1.
!@var dH2O  zonal H2O-prod.rate in kg/m^2/ppm_CH4/second in layer L
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: dH2O
!@var RSDIST,SIND,COSD orbit related variables computed once a day
      REAL*8 :: RSDIST,SIND,COSD
!@var depoBC,depoBC_1990 observed black carbon deposition (curr,1990)
      REAL*8, DIMENSION(72,46) :: depoBC,depoBC_1990  ! rad.grid
!@var ALB is SRNFLB(1)/(SRDFLB(1)+1.D-20),PLAVIS,PLANIR,ALBVIS,ALBNIR,
!@+       SRRVIS,SRRNIR,SRAVIS,SRANIR (see RADIATION)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: ALB

!@var SALB (1.-broadband surface albedo) - saved in rsf
      REAL*8, POINTER, DIMENSION(:,:) :: SALB   ! = ALB(:,:,1)
!      EQUIVALENCE (SALB,ALB)

#ifdef ALTER_RADF_BY_LAT
!@var FULGAS_lat multiplicative factors for altering FULGAS by latitude
!@+ for non-transient runs. (greenhouse gas regional forcing)
      real*8, dimension(13,46):: FULGAS_lat !rad not model grid, 13 gasses
!@var FS8OPX_lat multiplicative factors for altering FS8OPX by latitude
!@+ for non-transient runs. (aerosol regional forcing) SOLAR
!@var FT8OPX_lat multiplicative factors for altering FT8OPX by latitude
!@+ for non-transient runs. (aerosol regional forcing) THERMAL
      real*8, dimension(8,46):: FS8OPX_lat,FT8OPX_lat !rad not model grid
                                                !8 groups of aerosols
#endif
!@dbparam rad_interact_tr =1 for radiatively active tracers (default=0)
      INTEGER :: rad_interact_tr = 0
C**** the radiative forcing level for instantaneous forcing calcs is set
C**** using the rad_forc_lev parameter.
!@dbparam rad_forc_lev = 0 for TOA, 1 for LTROPO (default=0)
      INTEGER :: rad_forc_lev = 0

!@dbparam cloud_rad_forc = 1 for calculation of cloud radiative forcing
      INTEGER :: cloud_rad_forc = 0

!@var co2ppm Current CO2 level as seen by radiation
      REAL*8 :: co2ppm = 280.    ! set a resaonable default value

#ifdef CHL_from_SeaWIFs
!@var ACHL,ECHL1,ECHL0,BCHL,CCHL arrays for the reading in chlorophyll
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ACHL,ECHL1,ECHL0,BCHL,
     *     CCHL
!@var iu_CHL unit for chlorophyll file
      INTEGER iu_CHL
#endif
#ifdef OBIO_RAD_coupling
      REAL*8, ALLOCATABLE, DIMENSION(:)   ::  wfac
#endif

C**** Local variables initialised in init_RAD
!@var PLB0,QL0 global parts of local arrays (to avoid OMP-copyin)
      REAL*8, DIMENSION(LM_REQ)       :: PLB0,SHL0
!@var SINJ,COSJ sines and cosines for zenith angle calculation
      REAL*8, ALLOCATABLE, DIMENSION(:) :: SINJ,COSJ

!@var NTRIX Indexing array for optional aerosol interaction
      INTEGER, DIMENSION(ITRMAX) :: NTRIX = 0
!@var WTTR weighting array for optional aerosol interaction
      REAL*8, DIMENSION(ITRMAX) :: WTTR = 1.
!@var nrad_clay index of clay in arrays for optional aerosol interaction
      INTEGER :: nrad_clay

      END MODULE RAD_COM

      SUBROUTINE ALLOC_RAD_COM(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel
!@ver  1.0

      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE DOMAIN_DECOMP, ONLY : GET
      USE MODEL_COM, ONLY : IM, JM, LM, LM_REQ
#ifdef TRACERS_ON
      USE tracer_com,ONLY : Ntm
#endif
      USE RAD_COM, ONLY : RQT,Tchg,SRHR,TRHR,FSF,FSRDIR,SRVISSURF,TRSURF
     *     ,SRDN, CFRAC, RCLD, O3_tracer_save,rad_to_chem,rad_to_file
     *     ,KLIQ, COSZ1, dH2O, ALB, SALB, SINJ, COSJ,srnflb_save,
     *     trnflb_save, ttausv_save, ttausv_cs_save
     *     ,FSRDIF,DIRNIR,DIFNIR
#ifdef CHL_from_SeaWIFs
     *     ,achl,echl1,echl0,bchl,cchl
#endif
#ifdef OBIO_RAD_coupling
     *     ,wfac
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_0H, J_1H
      INTEGER :: IER

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE( RQT(LM_REQ, IM, J_0H:J_1H),
     *     Tchg(LM+LM_REQ, IM, J_0H:J_1H),
     *     SRHR(0:LM, IM, J_0H:J_1H),
     *     TRHR(0:LM, IM, J_0H:J_1H),
     *     TRSURF(1:4, IM, J_0H:J_1H),
     *     FSF(4,IM, J_0H:J_1H),
     *     FSRDIR(IM, J_0H:J_1H),
     *     SRVISSURF(IM, J_0H:J_1H),
     *     FSRDIF(IM, J_0H:J_1H),
     *     DIRNIR(IM, J_0H:J_1H),
     *     DIFNIR(IM, J_0H:J_1H),
     *     SRDN(IM, J_0H:J_1H),
     *     CFRAC(IM, J_0H:J_1H),
     *     RCLD(LM, IM, J_0H:J_1H),
     *     O3_tracer_save(LM, IM, J_0H:J_1H),
     *     rad_to_chem(LM, IM, J_0H:J_1H,5),
     *     rad_to_file(LM, IM, J_0H:J_1H,5),
     *     KLIQ(LM,4, IM, J_0H:J_1H),
     *     COSZ1(IM, J_0H:J_1H),
     *     dH2O(J_0H:J_1H, LM, 12),
     *     ALB(IM, J_0H:J_1H, 9),
     *     SINJ(J_0H:J_1H),
     *     COSJ(J_0H:J_1H),
     &     srnflb_save(Im,J_0H:J_1H,Lm),
     &     trnflb_save(Im,J_0H:J_1H,Lm),
#ifdef TRACERS_ON
     &     ttausv_save(Im,J_0H:J_1H,Ntm,Lm),
     &     ttausv_cs_save(Im,J_0H:J_1H,Ntm,Lm),
#endif
#ifdef CHL_from_SeaWIFs
     &         ACHL(IM,J_0H:J_1H),   
     &         ECHL1(IM,J_0H:J_1H),
     &         ECHL0(IM,J_0H:J_1H),
     &         BCHL(IM,J_0H:J_1H),
     &         CCHL(IM,J_0H:J_1H),
#endif
#ifdef OBIO_RAD_coupling
     &         wfac(1:33),
#endif
     *     STAT=IER)

      KLIQ = 1
      dH2O = 0.
      SALB => ALB(:,:,1)
      SRVISSURF = 0
      RETURN
      END SUBROUTINE ALLOC_RAD_COM

      SUBROUTINE io_rad(kunit,iaction,ioerr)
!@sum  io_rad reads and writes radiation arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun,ioread_single
     *         ,lhead,Kradia,irsficnt,irsficno
#ifdef TRACERS_ON
      USE tracer_com,ONLY : Ntm
#endif
      USE RAD_COM
      USE PARAM
      USE DOMAIN_DECOMP, ONLY : GRID, GET, AM_I_ROOT
      USE DOMAIN_DECOMP, ONLY : UNPACK_COLUMN, PACK_COLUMN
      USE DOMAIN_DECOMP, ONLY : UNPACK_BLOCK , PACK_BLOCK
      USE DOMAIN_DECOMP, ONLY : UNPACK_DATA  , PACK_DATA
      USE DOMAIN_DECOMP, ONLY : ESMF_BCAST
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "RAD05"
!@var HEADER_F Character string label for records (forcing runs)
      CHARACTER*80 :: HEADER_F, MODULE_HEADER_F = "RADF"

      REAL*8 :: Tchg_glob(LM+LM_REQ,IM,JM)
      REAL*8 :: RQT_glob(LM_REQ,IM,JM)
      INTEGER :: kliq_glob(LM,4,IM,JM)
      REAL*8, DIMENSION(0:LM,IM,JM) :: SRHR_GLOB, TRHR_GLOB
      REAL*8, DIMENSION(1:4,IM,JM) :: TRSURF_GLOB
      REAL*8  :: FSF_GLOB(4,IM,JM)
      REAL*8, DIMENSION(IM, JM) :: FSRDIR_GLOB, SRVISSURF_GLOB,
     &           SRDN_GLOB, CFRAC_GLOB, SALB_GLOB
      REAL*8, DIMENSION(LM, IM, JM) :: RCLD_GLOB,O3_tracer_save_GLOB
      REAL*8, DIMENSION(LM, IM, JM, 5) :: rad_to_chem_GLOB
      REAL*8,DIMENSION(Im,Jm,Lm) :: srnflb_save_glob,trnflb_save_glob
#ifdef TRACERS_ON
      REAL*8,DIMENSION(Im,Jm,Ntm,Lm) :: ttausv_save_glob,
     &     ttausv_cs_save_glob
#endif
      INTEGER :: J_0,J_1
      INTEGER :: k

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if (kradia.gt.0) then
        write (MODULE_HEADER_F(lhead+1:80),'(a8,i2,a15,i2,a7)')
     *    'R8 Tchg(',lm+LM_REQ,',ijm), I4 KLIQ(',lm,',4,ijm)'
        SELECT CASE (IACTION)
        CASE (:IOWRITE)            ! output to standard restart file
          CALL PACK_COLUMN(grid, Tchg, Tchg_glob)
          Call PACK_BLOCK( grid, kliq, kliq_glob)
          IF (AM_I_ROOT())
     &       WRITE (kunit,err=10) MODULE_HEADER_F,Tchg_glob,kliq_glob
        CASE (IOREAD:)
          SELECT CASE  (IACTION)
          CASE (ioread,irerun,ioread_single)  ! input for restart
            if (AM_I_ROOT() )
     &        READ (kunit,err=10) HEADER_F,Tchg_glob,kliq_glob
            CALL UNPACK_COLUMN(grid, Tchg_glob, Tchg)
            CALL UNPACK_BLOCK( grid, kliq_glob, kliq)
          CASE (IRSFIC)  ! only way to start frc. runs
            if (AM_I_ROOT() )
     &        READ (kunit,err=10) HEADER,RQT_glob,kliq_glob ; Tchg = 0.
            CALL UNPACK_COLUMN(grid, RQT_glob, RQT)
            CALL UNPACK_BLOCK( grid, kliq_glob, kliq)
            call sync_param( "Ikliq", Ikliq )
            if(Ikliq.eq.1) kliq=1  ! hysteresis init: equilibrium
            if(Ikliq.eq.0) kliq=0  ! hysteresis init: dry
          END SELECT
        END SELECT
      else

      MODULE_HEADER(lhead+1:80) = 'R8:Tq(3,ijm), I:LIQ(lm,4,ijm),'//
     *  ' S0,s+tHr,fs(4),fd,sv,al,sd,cf,o3'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_COLUMN(grid,  RQT,  RQT_glob)
        Call PACK_BLOCK( grid, kliq, kliq_glob)
        CALL PACK_COLUMN(grid, SRHR, SRHR_GLOB)
        CALL PACK_COLUMN(grid, TRHR, TRHR_GLOB)
        CALL PACK_COLUMN(grid, TRSURF, TRSURF_GLOB)
        CALL PACK_COLUMN(grid,  FSF,  FSF_GLOB)

        CALL PACK_DATA( grid, SALB     , SALB_GLOB)
        CALL PACK_DATA( grid, FSRDIR   , FSRDIR_GLOB)
        CALL PACK_DATA( grid, SRVISSURF, SRVISSURF_GLOB)
        CALL PACK_DATA( grid, SRDN     ,   SRDN_GLOB)
        CALL PACK_DATA( grid, CFRAC    ,  CFRAC_GLOB)

        CALL PACK_COLUMN(grid, RCLD          , RCLD_GLOB)
#ifdef TRACERS_SPECIAL_Shindell
        CALL PACK_COLUMN(grid, O3_tracer_save, O3_tracer_save_GLOB)
        CALL PACK_COLUMN(grid, rad_to_chem, rad_to_chem_GLOB)
#endif
        CALL PACK_DATA(grid,srnflb_save,srnflb_save_glob)
        CALL PACK_DATA(grid,trnflb_save,trnflb_save_glob)
#ifdef TRACERS_ON
        CALL PACK_DATA(grid,ttausv_save,ttausv_save_glob)
        CALL PACK_DATA(grid,ttausv_cs_save,ttausv_cs_save_glob)
#endif

        IF (AM_I_ROOT())
     *     WRITE (kunit,err=10) MODULE_HEADER,RQT_GLOB,KLIQ_GLOB
  ! rest needed only if MODRAD>0 at restart
     *      ,S0,  SRHR_GLOB, TRHR_GLOB, FSF_GLOB , FSRDIR_GLOB
     *       ,SRVISSURF_GLOB, SALB_GLOB, SRDN_GLOB, CFRAC_GLOB,
     *       RCLD_GLOB,  TRSURF_GLOB
#ifdef TRACERS_SPECIAL_Shindell
     *      ,O3_tracer_save_GLOB,rad_to_chem_GLOB
#endif
#ifdef TRACERS_DUST
     &      ,srnflb_save_glob,trnflb_save_glob,ttausv_save_glob
     &      ,ttausv_cs_save_glob
#endif
      CASE (IOREAD:)
        SELECT CASE  (IACTION)
        CASE (ioread,IRERUN)  ! input for restart, rerun or extension
          if (AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,RQT_glob,KLIQ_glob
  ! rest needed only if MODRAD>0 at restart
  ! rest needed only if MODRAD>0 at restart
     *       ,S0,  SRHR_GLOB, TRHR_GLOB, FSF_GLOB , FSRDIR_GLOB
     *       ,SRVISSURF_GLOB, SALB_GLOB, SRDN_GLOB, CFRAC_GLOB,
     *        RCLD_GLOB ,TRSURF_GLOB
#ifdef TRACERS_SPECIAL_Shindell
     *       ,O3_tracer_save_GLOB,rad_to_chem_GLOB
#endif
#ifdef TRACERS_DUST
     &       ,srnflb_save_glob,trnflb_save_glob,ttausv_save_glob
     &       ,ttausv_cs_save_glob
#endif
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER,
     *               MODULE_HEADER
              GO TO 10
            END IF
          end if

          CALL ESMF_BCAST(grid, S0)
          CALL UNPACK_COLUMN(grid,  RQT_glob,  RQT)
          Call UNPACK_BLOCK( grid, kliq_glob, kliq)
          CALL UNPACK_COLUMN(grid, SRHR_glob, SRHR)
          CALL UNPACK_COLUMN(grid, TRHR_glob, TRHR)
          CALL UNPACK_COLUMN(grid, TRSURF_glob, TRSURF)
          CALL UNPACK_COLUMN(grid,  FSF_glob,  FSF)

          CALL UNPACK_DATA( grid, SALB_glob     , SALB)
          CALL UNPACK_DATA( grid, FSRDIR_glob   , FSRDIR)
          CALL UNPACK_DATA( grid, SRVISSURF_glob, SRVISSURF)
          CALL UNPACK_DATA( grid, SRDN_glob     ,   SRDN)
          CALL UNPACK_DATA( grid, CFRAC_glob    ,  CFRAC)

          CALL UNPACK_COLUMN(grid, RCLD_glob          ,RCLD)
#ifdef TRACERS_SPECIAL_Shindell
          CALL UNPACK_COLUMN(grid, O3_tracer_save_glob,
     &         O3_tracer_save)
          CALL UNPACK_COLUMN(grid, rad_to_chem_glob, rad_to_chem)
#endif
          CALL UNPACK_DATA(grid,srnflb_save_glob,srnflb_save)
          CALL UNPACK_DATA(grid,trnflb_save_glob,trnflb_save)
#ifdef TRACERS_ON
          CALL UNPACK_DATA(grid,ttausv_save_glob,ttausv_save)
          CALL UNPACK_DATA(grid,ttausv_cs_save_glob,ttausv_cs_save)
#endif


        CASE (IRSFIC,irsficnt,IRSFICNO)  ! restart file of prev. run
          if (AM_I_ROOT() )
     &      READ (kunit,err=10) HEADER,RQT_glob,KLIQ_glob
          CALL UNPACK_COLUMN(grid,  RQT_glob,  RQT)
          Call UNPACK_BLOCK( grid, kliq_glob, kliq)
          call sync_param( "Ikliq", Ikliq )
          if(Ikliq.eq.1) kliq=1  ! hysteresis init: equilibrium
          if(Ikliq.eq.0) kliq=0  ! hysteresis init: dry
        END SELECT
      END SELECT
      end if

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_rad
