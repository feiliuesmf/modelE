#include "rundeck_opts.h"
#ifdef SKIP_TRACERS_RAD
#undef TRACERS_ON
#endif
      MODULE RAD_COM
!@sum  RAD_COM Model radiation arrays and parameters
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,lm_req
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
C**** Note: setting variable_orb_par=0, orb_par_year_bp=-50 (=year 2000)
C**** does not produce exactly the same as the default values.
!@dbparam variable_orb_par 1 if orbital parameters are time dependent
!@+       1 : use orb par from year "JYEAR - orb_par_year_bp"
!@+       0 : use orb par from year orb_par_year_bp (BP=before 1950)
!@+      -1 : set eccn/obliq/omegt to orb_par(1:3)
!@+    else : set eccn/obliq/omegt to defaults of orb_par
      integer :: variable_orb_par = -2
!@dbparam orb_par_year_bp = offset from model_year or 1950 (fixed case)
      integer :: orb_par_year_bp = 0
!@dbparam orb_par :: directly specifies orbital parameters
      real*8, dimension(3) :: orb_par = (/ eccn_def, obliq_def,
     *     omegt_def /)

!@var dimrad_sv dimension sum of input fields saved for radia_only runs
      INTEGER, PARAMETER :: dimrad_sv=IM*JM*(7*LM+3*LM_REQ+24)
!@var RQT Radiative equilibrium temperatures above model top
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RQT
!@var Tchg Total temperature change in adjusted forcing runs
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: Tchg
!@var SRHR(0) Solar   raditive net flux into the ground          (W/m^2)
!@var TRHR(0) Thermal raditive downward flux into ground(W/O -StB*T^4)(W/m^2)
!@*   Note: -StB*T^4 is added in SURFACE, since T varies betw. rad. calls
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
!@var TAUSUMW,TAUSUMI column-sum water,ice cloud opt. depths (for diags)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: TAUSUMW,TAUSUMI
!@var ttausv_save  Tracer optical thickness
!@var ttausv_cs_save  Tracer optical thickness clear sky
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: ttausv_save,
     &     ttausv_cs_save
#ifdef mjo_subdd
!@var OLR_acc, OLR_cnt --  Net thermal radiation at TOA (W/m^2) for SUBDD
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: OLR_acc
      REAL*8 :: OLR_cnt = 0.d0
!@var SWHR,LWHR,SWHR_cnt,LWHR_cnt -- shortwave/longwave heating rates for SUBDD (C/d)
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: SWHR,LWHR
      REAL*8 :: SWHR_cnt = 0.d0
      REAL*8 :: LWHR_cnt = 0.d0
!@var swu_avg,swu_cnt -- upward shortwave fluxes at srf for SUBDD (C/d)
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: swu_avg
      REAL*8 :: swu_cnt = 0.d0
#endif
#ifdef TRACERS_ON
!@var ttausv_sum(_cs) daily sum opt depth by tracer, all (clear) sky
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: ttausv_sum,ttausv_sum_cs
      real*8 :: ttausv_count = 0.d0
!@var nTracerRadiaActive number of radiatively active tracers (always .le. Ntm)
      integer :: nTracerRadiaActive
!@var tracerRadiaActiveFlag array of flags of dimension Ntm, whose elements
!@+                         are set to .true. for radiatively active tracers
      logical,allocatable,dimension(:) :: tracerRadiaActiveFlag
!@var aerAbs6SaveInst Band 6 sum over aerosols  of extinction-scattering,
!@+   saved for instantaneous SUBDDiag output
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: aerAbs6SaveInst
#ifdef TRACERS_SPECIAL_Shindell
!@var maxNtraceFastj max expected rad code tracers passed to photolysis
!@var ttausv_ntrace Tracer optical thickness saved 1:NTRACE not 1:NTM
!@+   This is so clays are separate. Only needed for chemistry on.
      integer, parameter :: maxNtraceFastj=17 
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: ttausv_ntrace
#endif
#endif
!@var CFRAC Total cloud fraction as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: CFRAC ! saved in rsf
!@var RCLD Total cloud optical depth as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RCLD ! saved in rsf
!@var chem_tracer_save 3D O3, CH4 saved elsewhere for use in radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: chem_tracer_save!saved rsf
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
!@var stratO3_tracer_save 3D stratOx saved elsewhere for use in rad code
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::stratO3_tracer_save!saved rsf
#endif
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
!@var COSZ_day Mean Solar Zenith angle for current day
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ_day
!@var SUNSET Time of sunset for current day (radians from local noon)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SUNSET
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
!@+   special cases: Volc_yr=-1    : 150-yr mean 1850-1999
!@+                  Volc_yr=-2010 : current year up to 2010 then
!@+                                  repeat volcanos from 100 yrs ago
!@+                  Volc_yr=-2000 : older way of creating future volc
      INTEGER :: Volc_yr = 1951 , Volc_day = 182
!@dbparam Aero_yr obs.year of troposph.Aerosols (if 0: use current yr)
      INTEGER :: Aero_yr = 1951    ! always use annual cycle
!@dbparam O3_yr obs.year of Ozone (if 0: use current year)
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
!@dbparam rad_interact_aer =1 for radiatively active non-chem tracers
      INTEGER :: rad_interact_aer = 0  ! defaults to 0
!@dbparam clim_interact_chem=1 for radiatively active chem tracers
!@+ also affects chemisty to humidity feedback
      INTEGER :: clim_interact_chem= 0  ! defaults to 0

!@dbparam nradfrc sets frequency of inst. rad. forcing calculations
      INTEGER :: nradfrc=1 ! do them every nrad*nradfrc physics steps
!                nradfrc=0: skip all, no repeated radiation calculations
C**** the radiative forcing level for instantaneous forcing calcs is set
C**** using the rad_forc_lev parameter.
!@dbparam rad_forc_lev = 0 for TOA, 1 for LTROPO (default=0)
      INTEGER :: rad_forc_lev = 0

!@dbparam cloud_rad_forc = 1 for calculation of cloud radiative forcing
      INTEGER :: cloud_rad_forc = 0

!@dbparam aer_rad_forc = 1 for calculation of aerosol radiative forcing
!@+   note this only works for background aerosols, not tracers
      INTEGER :: aer_rad_forc = 0

!@var co2ppm Current CO2 level as seen by radiation
      REAL*8 :: co2ppm = 280.    ! set a reasonable default value

#ifdef CHL_from_SeaWIFs
!@var ACHL,ECHL1,ECHL0,BCHL,CCHL arrays for the reading in chlorophyll
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ACHL,ECHL1,ECHL0,BCHL,
     *     CCHL
!@var iu_CHL unit for chlorophyll file
      INTEGER iu_CHL
#endif
#if (defined CHL_from_SeaWIFs) || (defined OBIO_RAD_coupling)
      REAL*8,  DIMENSION(33)   ::  wfac
#endif

C**** Local variables initialised in init_RAD
!@var PLB0,QL0 global parts of local arrays (to avoid OMP-copyin)
      REAL*8, DIMENSION(LM_REQ)       :: PLB0,SHL0
!@var NTRIX Indexing array for optional aerosol interaction
      INTEGER, DIMENSION(ITRMAX) :: NTRIX = 0
!@var WTTR weighting array for optional aerosol interaction
      REAL*8, DIMENSION(ITRMAX) :: WTTR = 1.
!@var nrad_clay index of clay in arrays for optional aerosol interaction
      INTEGER :: nrad_clay

#ifdef CUBED_SPHERE
!@var JM_DH2O number of latitudes in CH4->H2O input file
!@var LAT_DH2O latitudes in CH4->H2O input file (converted to radians)
      integer, parameter :: jm_dh2o=18
      real*8 :: lat_dh2o(jm_dh2o)
#endif

      END MODULE RAD_COM

      SUBROUTINE ALLOC_RAD_COM(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel
!@ver  1.0

      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE DOMAIN_DECOMP_ATM, ONLY : GET
      USE MODEL_COM, ONLY : IM, JM, LM, LM_REQ
#ifdef TRACERS_ON
      USE tracer_com,ONLY : Ntm
#endif
      USE RAD_COM, ONLY : RQT,Tchg,SRHR,TRHR,FSF,FSRDIR,SRVISSURF,TRSURF
     *     ,SRDN, CFRAC, RCLD, chem_tracer_save,rad_to_chem,rad_to_file
     *     ,KLIQ, COSZ1, COSZ_day, SUNSET, dH2O, ALB, SALB
     *     ,srnflb_save, trnflb_save, ttausv_save, ttausv_cs_save
     *     ,FSRDIF,DIRNIR,DIFNIR,TAUSUMW,TAUSUMI
#ifdef mjo_subdd
     *     ,SWHR_cnt,LWHR_cnt,SWHR,LWHR,OLR_acc,OLR_cnt
     *     ,swu_avg,swu_cnt
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *     ,ttausv_ntrace,maxNtraceFastj
#endif
#ifdef CHL_from_SeaWIFs
     *     ,achl,echl1,echl0,bchl,cchl
#endif
#if (defined CHL_from_SeaWIFs) || (defined OBIO_RAD_coupling)
     *     ,wfac
#endif
#ifdef TRACERS_ON
     *     ,ttausv_sum,ttausv_sum_cs,ttausv_count,nTracerRadiaActive
     &     ,tracerRadiaActiveFlag,aerAbs6SaveInst
#endif
#ifdef CUBED_SPHERE
     &     ,JM_DH2O
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &     ,stratO3_tracer_save
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_0H, J_1H
      INTEGER :: IER

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE( RQT(LM_REQ, I_0H:I_1H, J_0H:J_1H),
     *     Tchg(LM+LM_REQ, I_0H:I_1H, J_0H:J_1H),
     *     SRHR(0:LM, I_0H:I_1H, J_0H:J_1H),
     *     TRHR(0:LM, I_0H:I_1H, J_0H:J_1H),
     *     TRSURF(4, I_0H:I_1H, J_0H:J_1H),
     *     FSF(4,I_0H:I_1H, J_0H:J_1H),
     *     FSRDIR(I_0H:I_1H, J_0H:J_1H),
     *     SRVISSURF(I_0H:I_1H, J_0H:J_1H),
     *     FSRDIF(I_0H:I_1H, J_0H:J_1H),
     *     DIRNIR(I_0H:I_1H, J_0H:J_1H),
     *     DIFNIR(I_0H:I_1H, J_0H:J_1H),
     *     TAUSUMW(I_0H:I_1H, J_0H:J_1H),
     *     TAUSUMI(I_0H:I_1H, J_0H:J_1H),
     *     SRDN(I_0H:I_1H, J_0H:J_1H),
     *     CFRAC(I_0H:I_1H, J_0H:J_1H),
     *     RCLD(LM, I_0H:I_1H, J_0H:J_1H),
     *     chem_tracer_save(2,LM, I_0H:I_1H, J_0H:J_1H),
     *     rad_to_chem(5, LM, I_0H:I_1H, J_0H:J_1H),
     *     rad_to_file(5, LM, I_0H:I_1H, J_0H:J_1H),
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     *     stratO3_tracer_save(LM, I_0H:I_1H, J_0H:J_1H),
#endif
     *     KLIQ(LM,4, I_0H:I_1H, J_0H:J_1H),
     *     COSZ1   (I_0H:I_1H, J_0H:J_1H),
     *     COSZ_day(I_0H:I_1H, J_0H:J_1H),
     *     SUNSET  (I_0H:I_1H, J_0H:J_1H),
#ifdef CUBED_SPHERE
     *     dH2O(JM_DH2O, LM, 12),
#else
     *     dH2O(J_0H:J_1H, LM, 12),
#endif
     *     ALB(I_0H:I_1H, J_0H:J_1H, 9),
     &     srnflb_save(I_0H:I_1H,J_0H:J_1H,Lm),
     &     trnflb_save(I_0H:I_1H,J_0H:J_1H,Lm),
#ifdef mjo_subdd
     *     OLR_acc(I_0H:I_1H,J_0H:J_1H),
     *     SWHR(I_0H:I_1H,J_0H:J_1H,Lm),
     *     LWHR(I_0H:I_1H,J_0H:J_1H,Lm),
     *     swu_avg(I_0H:I_1H,J_0H:J_1H),
#endif
#ifdef TRACERS_ON
     &     ttausv_save(I_0H:I_1H,J_0H:J_1H,Ntm,Lm),
     &     ttausv_cs_save(I_0H:I_1H,J_0H:J_1H,Ntm,Lm),
     &     ttausv_sum(I_0H:I_1H,J_0H:J_1H,Ntm),
     &     ttausv_sum_cs(I_0H:I_1H,J_0H:J_1H,Ntm),
     &     tracerRadiaActiveFlag(Ntm),
     &     aerAbs6SaveInst(I_0H:I_1H,J_0H:J_1H,Lm),
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &     ttausv_ntrace(I_0H:I_1H,J_0H:J_1H,maxNtraceFastj,Lm),
#endif
#ifdef CHL_from_SeaWIFs
     &         ACHL(I_0H:I_1H,J_0H:J_1H),
     &         ECHL1(I_0H:I_1H,J_0H:J_1H),
     &         ECHL0(I_0H:I_1H,J_0H:J_1H),
     &         BCHL(I_0H:I_1H,J_0H:J_1H),
     &         CCHL(I_0H:I_1H,J_0H:J_1H),
#endif
     *     STAT=IER)

#ifdef mjo_subdd
      OLR_acc=0.
      OLR_cnt=0.
      SWHR=0.
      LWHR=0.
      SWHR_cnt=0.
      LWHR_cnt=0.
      swu_avg=0.
      swu_cnt=0.
#endif
      KLIQ = 1
      dH2O = 0.
      SALB => ALB(:,:,1)
      SRVISSURF = 0
      FSF=0
      TRSURF=0
#ifdef TRACERS_ON
      nTracerRadiaActive=0
      tracerRadiaActiveFlag=.false.
#endif

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
      USE Dictionary_mod
      USE DOMAIN_DECOMP_1D, ONLY : GRID, GET, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, ONLY : UNPACK_COLUMN, PACK_COLUMN
      USE DOMAIN_DECOMP_1D, ONLY : UNPACK_BLOCK , PACK_BLOCK
      USE DOMAIN_DECOMP_1D, ONLY : UNPACK_DATA  , PACK_DATA
      USE DOMAIN_DECOMP_1D, ONLY : ESMF_BCAST
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "RAD05"
!@var HEADER_F Character string label for records (forcing runs)
      CHARACTER*80 :: HEADER_F, MODULE_HEADER_F = "RADF"

      real*8, dimension(:,:,:), allocatable ::
     &     Tchg_glob !(LM+LM_REQ,IM,JM)
     &    ,RQT_glob  !(LM_REQ,IM,JM)
     &    ,RCLD_GLOB !(LM, IM, JM)
     &    ,SRHR_GLOB !(0:LM,IM,JM)
     &    ,TRHR_GLOB !(0:LM,IM,JM)
     &    ,TRSURF_GLOB !(1:4,IM,JM)
     &    ,FSF_GLOB    !(4,IM,JM)
      integer, dimension(:,:,:,:), allocatable :: kliq_glob!(LM,4,IM,JM)
      real*8, dimension(:,:), allocatable :: ! IM,JM arrays
     &     FSRDIR_GLOB, SRVISSURF_GLOB,
     &     SRDN_GLOB, CFRAC_GLOB, SALB_GLOB,
     &     FSRDIF_GLOB, DIRNIR_GLOB, DIFNIR_GLOB
#ifdef TRACERS_ON
      REAL*8,DIMENSION(:,:,:), allocatable :: !(Im,Jm,Lm)
     &     aerAbs6SaveInst_glob
#ifdef TRACERS_SPECIAL_Shindell
      REAL*8,DIMENSION(:,:,:,:),allocatable::chem_tracer_save_GLOB !(2,LM,IM,JM)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      REAL*8,DIMENSION(:,:,:),allocatable::stratO3_tracer_save_GLOB !(LM,IM,JM)
#endif
      REAL*8,DIMENSION(:,:,:,:),allocatable::rad_to_chem_GLOB ! !(5,LM,IM,JM)
#endif
#ifdef TRACERS_DUST
      REAL*8,DIMENSION(:,:,:), allocatable :: !(Im,Jm,Lm)
     &     srnflb_save_glob,trnflb_save_glob
#endif
      REAL*8,DIMENSION(:,:,:,:), allocatable :: ! (Im,Jm,Ntm,Lm)
     &     ttausv_save_glob,ttausv_cs_save_glob
#ifdef TRACERS_SPECIAL_Shindell
      REAL*8,DIMENSION(:,:,:,:), allocatable :: ttausv_ntrace_glob
#endif
      REAL*8,DIMENSION(:,:,:), allocatable ::   !(Im,Jm,Ntm)
     &     ttausv_sum_glob,ttausv_sum_cs_glob
#endif
      INTEGER :: J_0,J_1
      integer :: img, jmg, lmg

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if(am_i_root()) then
         img = IM
         jmg = JM
         lmg = LM
      else ! small unused global array (but must be allocated for MPI)
         img = 1
         jmg = 1
         lmg = 1
      end if
      allocate(Tchg_glob(lmg+LM_REQ,img,jmg))
      allocate(RQT_glob(LM_REQ,img,jmg))
      allocate(kliq_glob(lmg,4,img,jmg))
      allocate(SRHR_GLOB(0:lmg,img,jmg))
      allocate(TRHR_GLOB(0:lmg,img,jmg))
      allocate(TRSURF_GLOB(1:4,img,jmg))
      allocate(FSF_GLOB(4,img,jmg))
      allocate(FSRDIR_GLOB(img,jmg),
     &     SRVISSURF_GLOB(img,jmg),
     &     SRDN_GLOB(img,jmg),
     &     CFRAC_GLOB(img,jmg),
     &     SALB_GLOB(img,jmg),
     &     FSRDIF_GLOB(img,jmg),
     &     DIRNIR_GLOB(img,jmg),
     &     DIFNIR_GLOB(img,jmg))
      allocate(RCLD_GLOB(lmg,img,jmg))
#ifdef TRACERS_ON
      allocate(aerAbs6SaveInst_glob(img,jmg,lmg))
#ifdef TRACERS_SPECIAL_Shindell
      allocate(chem_tracer_save_GLOB(2,lmg, img,jmg))
      allocate(rad_to_chem_GLOB(5,lmg,img,jmg))
      allocate(ttausv_ntrace_glob(img,jmg,maxNtraceFastj,lmg))
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      allocate(stratO3_tracer_save_GLOB(lmg, img, jmg))
#endif
#endif
#ifdef TRACERS_DUST
      allocate(srnflb_save_glob(img,jmg,lmg))
      allocate(trnflb_save_glob(img,jmg,lmg))
#endif
      allocate(ttausv_save_glob(img,jmg,Ntm,lmg))
      allocate(ttausv_cs_save_glob(Im,Jm,Ntm,lmg))
      allocate(ttausv_sum_glob(img,jmg,Ntm))
      allocate(ttausv_sum_cs_glob(img,jmg,Ntm))
#endif

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
!!!
        CALL PACK_DATA( grid, FSRDIF, FSRDIF_GLOB)
        CALL PACK_DATA( grid, DIRNIR, DIRNIR_GLOB)
        CALL PACK_DATA( grid, DIFNIR, DIFNIR_GLOB)
!!!
        CALL PACK_DATA( grid, SRDN     ,   SRDN_GLOB)
        CALL PACK_DATA( grid, CFRAC    ,  CFRAC_GLOB)

        CALL PACK_COLUMN(grid, RCLD          , RCLD_GLOB)
#ifdef TRACERS_SPECIAL_Shindell
        CALL PACK_BLOCK(grid, chem_tracer_save, chem_tracer_save_GLOB)
        CALL PACK_BLOCK(grid, rad_to_chem, rad_to_chem_GLOB)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        CALL PACK_COLUMN
     &  (grid, stratO3_tracer_save, stratO3_tracer_save_GLOB)
#endif
#endif
#ifdef TRACERS_DUST
        CALL PACK_DATA(grid,srnflb_save,srnflb_save_glob)
        CALL PACK_DATA(grid,trnflb_save,trnflb_save_glob)
#endif
#ifdef TRACERS_ON
        CALL PACK_DATA(grid,aerAbs6SaveInst,aerAbs6SaveInst_glob)
        CALL PACK_DATA(grid,ttausv_save,ttausv_save_glob)
        CALL PACK_DATA(grid,ttausv_cs_save,ttausv_cs_save_glob)
        CALL PACK_DATA(grid,ttausv_sum,ttausv_sum_glob)
        CALL PACK_DATA(grid,ttausv_sum_cs,ttausv_sum_cs_glob)
#endif
#ifdef TRACERS_SPECIAL_Shindell
        CALL PACK_DATA(grid,ttausv_ntrace,ttausv_ntrace_glob)
#endif

        IF (AM_I_ROOT())
     *     WRITE (kunit,err=10) MODULE_HEADER,RQT_GLOB,KLIQ_GLOB
  ! rest needed only if MODRAD>0 at restart
     *      ,S0,  SRHR_GLOB, TRHR_GLOB, FSF_GLOB , FSRDIR_GLOB
     *       ,SRVISSURF_GLOB, SALB_GLOB, SRDN_GLOB, CFRAC_GLOB,
     *       RCLD_GLOB,  TRSURF_GLOB
     .      ,FSRDIF_GLOB, DIRNIR_GLOB, DIFNIR_GLOB
#ifdef TRACERS_SPECIAL_Shindell
     *      ,chem_tracer_save_GLOB,rad_to_chem_GLOB
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     *      ,stratO3_tracer_save_GLOB
#endif
#ifdef TRACERS_DUST
     &      ,srnflb_save_glob,trnflb_save_glob
#endif
#ifdef TRACERS_ON
     &      ,ttausv_sum_glob,ttausv_sum_cs_glob,ttausv_count
     &      ,ttausv_save_glob,ttausv_cs_save_glob,aerAbs6SaveInst_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &      ,ttausv_ntrace_glob
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
     .      ,FSRDIF_GLOB, DIRNIR_GLOB, DIFNIR_GLOB
#ifdef TRACERS_SPECIAL_Shindell
     *       ,chem_tracer_save_GLOB,rad_to_chem_GLOB
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     *       ,stratO3_tracer_save_GLOB
#endif
#ifdef TRACERS_DUST
     &       ,srnflb_save_glob,trnflb_save_glob
#endif
#ifdef TRACERS_ON
     &       ,ttausv_sum_glob,ttausv_sum_cs_glob,ttausv_count
     &       ,ttausv_save_glob,ttausv_cs_save_glob,aerAbs6SaveInst_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &       ,ttausv_ntrace_glob
#endif
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER,
     *               MODULE_HEADER
              GO TO 10
            END IF
          end if

#ifdef TRACERS_ON
          CALL ESMF_BCAST(grid, ttausv_count)
#endif
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
          CALL UNPACK_DATA( grid, FSRDIF_glob, FSRDIF)
          CALL UNPACK_DATA( grid, DIRNIR_glob, DIRNIR)
          CALL UNPACK_DATA( grid, DIFNIR_glob, DIFNIR)
          CALL UNPACK_DATA( grid, SRDN_glob     ,   SRDN)
          CALL UNPACK_DATA( grid, CFRAC_glob    ,  CFRAC)

          CALL UNPACK_COLUMN(grid, RCLD_glob          ,RCLD)
#ifdef TRACERS_SPECIAL_Shindell
          CALL UNPACK_BLOCK(grid, chem_tracer_save_glob,
     &         chem_tracer_save)
          CALL UNPACK_BLOCK(grid, rad_to_chem_glob, rad_to_chem)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
          CALL UNPACK_COLUMN(grid, stratO3_tracer_save_glob,
     &         stratO3_tracer_save)
#endif
#endif
#ifdef TRACERS_DUST
          CALL UNPACK_DATA(grid,srnflb_save_glob,srnflb_save)
          CALL UNPACK_DATA(grid,trnflb_save_glob,trnflb_save)
#endif
#ifdef TRACERS_ON
          CALL UNPACK_DATA(grid,ttausv_save_glob,ttausv_save)
          CALL UNPACK_DATA(grid,ttausv_cs_save_glob,ttausv_cs_save)
          CALL UNPACK_DATA(grid,ttausv_sum_glob,ttausv_sum)
          CALL UNPACK_DATA(grid,ttausv_sum_cs_glob,ttausv_sum_cs)
          CALL UNPACK_DATA(grid,aerAbs6SaveInst_glob,aerAbs6SaveInst)
#endif
#ifdef TRACERS_SPECIAL_Shindell
          CALL UNPACK_DATA(grid,ttausv_ntrace_glob,ttausv_ntrace)
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
      call freemem
      RETURN
 10   IOERR=1
      call freemem
      RETURN
      contains
      subroutine freemem
      deallocate(Tchg_glob)
      deallocate(RQT_glob)
      deallocate(kliq_glob)
      deallocate(SRHR_GLOB)
      deallocate(TRHR_GLOB)
      deallocate(TRSURF_GLOB)
      deallocate(FSF_GLOB)
      deallocate(FSRDIR_GLOB,
     &     SRVISSURF_GLOB,
     &     SRDN_GLOB,
     &     CFRAC_GLOB,
     &     SALB_GLOB,
     &     FSRDIF_GLOB,
     &     DIRNIR_GLOB,
     &     DIFNIR_GLOB)
      deallocate(RCLD_GLOB)
#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_Shindell
      deallocate(chem_tracer_save_GLOB)
      deallocate(rad_to_chem_GLOB)
      deallocate(ttausv_ntrace_glob)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      deallocate(stratO3_tracer_save_GLOB)
#endif
#endif
#ifdef TRACERS_DUST
      deallocate(srnflb_save_glob)
      deallocate(trnflb_save_glob)
#endif
      deallocate(ttausv_save_glob)
      deallocate(ttausv_cs_save_glob)
      deallocate(ttausv_sum_glob)
      deallocate(ttausv_sum_cs_glob)
      deallocate(aerAbs6SaveInst_glob)
#endif
      end subroutine freemem
      END SUBROUTINE io_rad

#ifdef NEW_IO
      subroutine def_rsf_rad(fid)
!@sum  def_rsf_rad defines radiation array structure in restart files
!@auth M. Kelley
!@ver  beta
      use rad_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id

      call defvar(grid,fid,s0,'s0')
      call defvar(grid,fid,rqt,'rqt(lm_req,dist_im,dist_jm)')
      call defvar(grid,fid,kliq,'kliq(lm,four,dist_im,dist_jm)')
      call defvar(grid,fid,srhr,'srhr(zero_to_lm,dist_im,dist_jm)')
      call defvar(grid,fid,trhr,'trhr(zero_to_lm,dist_im,dist_jm)')
      call defvar(grid,fid,trsurf,'trsurf(nstype,dist_im,dist_jm)')
      call defvar(grid,fid,fsf,'fsf(nstype,dist_im,dist_jm)')
      call defvar(grid,fid,fsrdir,'fsrdir(dist_im,dist_jm)')
      call defvar(grid,fid,srvissurf,'srvissurf(dist_im,dist_jm)')
      call defvar(grid,fid,srdn,'srdn(dist_im,dist_jm)')
      call defvar(grid,fid,cfrac,'cfrac(dist_im,dist_jm)')
      call defvar(grid,fid,salb,'salb(dist_im,dist_jm)')
      call defvar(grid,fid,fsrdif,'fsrdif(dist_im,dist_jm)')
      call defvar(grid,fid,dirnir,'dirnir(dist_im,dist_jm)')
      call defvar(grid,fid,difnir,'difnir(dist_im,dist_jm)')
      call defvar(grid,fid,rcld,'rcld(lm,dist_im,dist_jm)')
#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_Shindell
      call defvar(grid,fid,chem_tracer_save,
     &     'chem_tracer_save(two,lm,dist_im,dist_jm)')
      call defvar(grid,fid,rad_to_chem,
     &     'rad_to_chem(five,lm,dist_im,dist_jm)')
      call defvar(grid,fid,ttausv_ntrace,
     &     'ttausv_ntrace(dist_im,dist_jm,maxNtraceFastj,lm)')
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      call defvar(grid,fid,strato3_tracer_save,
     &     'strato3_tracer_save(lm,dist_im,dist_jm)')
#endif
#endif
#ifdef TRACERS_DUST
      call defvar(grid,fid,srnflb_save,
     &     'srnflb_save(dist_im,dist_jm,lm)')
      call defvar(grid,fid,trnflb_save,
     &     'trnflb_save(dist_im,dist_jm,lm)')
#endif
      call defvar(grid,fid,ttausv_save,
     &     'ttausv_save(dist_im,dist_jm,ntm,lm)')
      call defvar(grid,fid,ttausv_cs_save,
     &     'ttausv_cs_save(dist_im,dist_jm,ntm,lm)')
      call defvar(grid,fid,ttausv_sum,
     &     'ttausv_sum(dist_im,dist_jm,ntm)')
      call defvar(grid,fid,ttausv_sum_cs,
     &     'ttausv_sum_cs(dist_im,dist_jm,ntm)')
      call defvar(grid,fid,ttausv_count,'ttausv_count')
      call defvar(grid,fid,aerAbs6SaveInst,
     &     'aerAbs6SaveInst(dist_im,dist_jm,lm)')
#endif
      return
      end subroutine def_rsf_rad

      subroutine new_io_rad(fid,iaction)
!@sum  new_io_rad read/write radiation arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
#ifdef TRACERS_ON
      USE tracer_com , only : ntm
#endif
      use rad_com
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid, fid,'s0', s0)
        call write_dist_data(grid, fid,'rqt',  rqt, jdim=3)
        call write_dist_data(grid, fid,'kliq', kliq, jdim=4)
        call write_dist_data(grid, fid,'srhr', srhr, jdim=3)
        call write_dist_data(grid, fid,'trhr', trhr, jdim=3)
        call write_dist_data(grid, fid,'trsurf', trsurf, jdim=3)
        call write_dist_data(grid, fid,'fsf',  fsf, jdim=3)
        call write_dist_data(grid, fid,'salb', salb)
        call write_dist_data(grid, fid,'fsrdir', fsrdir)
        call write_dist_data(grid, fid,'srvissurf', srvissurf)
        call write_dist_data(grid, fid,'fsrdif', fsrdif)
        call write_dist_data(grid, fid,'dirnir', dirnir)
        call write_dist_data(grid, fid,'difnir', difnir)
        call write_dist_data(grid, fid,'srdn',   srdn)
        call write_dist_data(grid, fid,'cfrac',  cfrac)
        call write_dist_data(grid, fid,'rcld', rcld, jdim=3)
#ifdef TRACERS_SPECIAL_Shindell
        call write_dist_data(grid,fid,
     &       'chem_tracer_save', chem_tracer_save, jdim=4)
        call write_dist_data(grid,fid,
     &       'rad_to_chem', rad_to_chem, jdim=4)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        call write_dist_data(grid,fid,
     &       'strato3_tracer_save', strato3_tracer_save, jdim=3)
#endif
#endif
#ifdef TRACERS_DUST
        call write_dist_data(grid,fid,'srnflb_save',srnflb_save)
        call write_dist_data(grid,fid,'trnflb_save',trnflb_save)
#endif
#ifdef TRACERS_ON
        call write_dist_data(grid,fid,'ttausv_save',ttausv_save)
        call write_dist_data(grid,fid,'ttausv_cs_save',
     &       ttausv_cs_save)
        call write_data(grid, fid,'ttausv_count',ttausv_count)
        call write_dist_data(grid,fid,'ttausv_sum',ttausv_sum)
        call write_dist_data(grid,fid,'ttausv_sum_cs',ttausv_sum_cs)
        call write_dist_data(grid,fid,'aerAbs6SaveInst',aerAbs6SaveInst)
#endif
#ifdef TRACERS_SPECIAL_Shindell
        call write_dist_data(grid,fid,'ttausv_ntrace',ttausv_ntrace)
#endif
      case (ioread)
        call read_data(grid, fid,'s0', s0, bcast_all=.true.)
        call read_dist_data(grid, fid,'rqt',  rqt, jdim=3)
        call read_dist_data(grid, fid,'kliq', kliq, jdim=4)
        call read_dist_data(grid, fid,'srhr', srhr, jdim=3)
        call read_dist_data(grid, fid,'trhr', trhr, jdim=3)
        call read_dist_data(grid, fid,'trsurf', trsurf, jdim=3)
        call read_dist_data(grid, fid,'fsf',  fsf, jdim=3)
        call read_dist_data(grid, fid,'salb', salb)
        call read_dist_data(grid, fid,'fsrdir', fsrdir)
        call read_dist_data(grid, fid,'srvissurf', srvissurf)
        call read_dist_data(grid, fid,'fsrdif', fsrdif)
        call read_dist_data(grid, fid,'dirnir', dirnir)
        call read_dist_data(grid, fid,'difnir', difnir)
        call read_dist_data(grid, fid,'srdn',   srdn)
        call read_dist_data(grid, fid,'cfrac',  cfrac)
        call read_dist_data(grid, fid,'rcld', rcld, jdim=3)
#ifdef TRACERS_SPECIAL_Shindell
        call read_dist_data(grid,fid,
     &       'chem_tracer_save', chem_tracer_save, jdim=4)
        call read_dist_data(grid,fid,
     &       'rad_to_chem', rad_to_chem, jdim=4)
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        call read_dist_data(grid,fid,
     &       'strato3_tracer_save', strato3_tracer_save, jdim=3)
#endif
#endif
#ifdef TRACERS_DUST
        call read_dist_data(grid,fid,'srnflb_save',srnflb_save)
        call read_dist_data(grid,fid,'trnflb_save',trnflb_save)
#endif
#ifdef TRACERS_ON
        call read_dist_data(grid,fid,'ttausv_save',ttausv_save)
        call read_dist_data(grid,fid,'ttausv_cs_save',
     &       ttausv_cs_save)
        call read_data(grid, fid,'ttausv_count',ttausv_count,
     &       bcast_all=.true.)
        call read_dist_data(grid,fid,'ttausv_sum',ttausv_sum)
        call read_dist_data(grid,fid,'ttausv_sum_cs',ttausv_sum_cs)
        call read_dist_data(grid,fid,'aerAbs6SaveInst',aerAbs6SaveInst)
#endif
#ifdef TRACERS_SPECIAL_Shindell
        call read_dist_data(grid,fid,'ttausv_ntrace',ttausv_ntrace)
#endif
      end select
      return
      end subroutine new_io_rad
#endif /* NEW_IO */
