E4CubeTest.R GISS Model E coupled version         kelley 02/25/2010

E4CubeTest: C32 cubed-sphere version of E4M20 coupled to 4x5 13-layer ocean.
            Intended to be used for regression testing and demonstration only.

E4M20.R GISS Model E  2004 modelE                     rar     07/15/2009
E4M20: modelE as frozen in July 2009 without gravity wave drag
atmospheric composition from year 1850
time steps: physics 30 min.; radiation 2.5 hrs

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define CHECK_OCEAN                 ! needed to compile aux/file CMPE002
#define AG2OG_PRECIP_BUNDLE
#define OG2AG_TOC2SST_BUNDLE
#define AG2OG_OCEANS_BUNDLE
#define OG2AG_OCEANS_BUNDLE
#define BUNDLE_INTERP
#define NEW_IO
!#define CALC_GWDRAG                 ! need to make ZVAR for C32 res first
#define USE_ENT                  ! include dynamic vegetation model
#define SET_SOILCARBON_GLOBAL_TO_ZERO
!#define ROUGHL_HACK ! no longer needed?

End Preprocessor Options

Object modules: (in order of decreasing priority)
!
! Atmosphere/land/ice models:
!
RES_C32M20AT                        ! C32L20 resolution
MODEL_COM GNOM_CS IO_DRV            ! model variables, geometry, I/O driver
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! allocate distributed arrays
DOMAIN_DECOMPcs                     ! cubed sphere domain decomposition for atm. routines
ATMDYN_COM ATM_DUM                  ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
FV_UTILS FV_CS_Mod FV_INTERFACE     ! FV dynamical core wrapper
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
QUScubed                            ! cubed-sphere adaptation of QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV ! + component giss_LSM: land surface and soils
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + component Ent: new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB_E1                            ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
!SNOW_DRV SNOW                      ! snow model
RAD_COM RAD_DRV RADIATION COSZ_2D   ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC QUICKPRT       ! diagnostics
DIAG_RES_M                          ! diagnostics (resolution dependent)
GCDIAGcs DIAG_ZONALcs               ! grid-dependent code for lat-circle diags
POUT                                ! not used for cubed sphere
!
! Sea ice dynamics
!
ICEDYN ICEDYN_DRV ADVSIcs           ! note ADVSIcs runs on the atm grid
!
! Ocean model
!
RES_5x4_L13                         ! ocean horiz res 4x5deg, 13 layers
ODIAG_COM OCEAN_COM OSTRAITS_COM OGEOM   ! dynamic ocean modules
OCNDYN OCNDYN2                           ! dynamic ocean routines
OCN_Interp OSTRAITS OCNGM OCNKPP         ! dynamic ocean routines
OCEANR_DIM AFLUXES OFLUXES               ! dynamic ocean routines
ODIAG_PRT                              ! ocean diagnostic print out
OCNFUNTAB                           ! ocean function look up table
SparseCommunicator_mod              ! sparse gather/scatter module used for straits
!
! Utilities
!
cs2ll_utils                         ! CS utilities for diags/flux coupler/regrids
regrid_com                          ! more regrid routines
FFTW_COM OFFT72E                    ! FFTs
pario_fbsa                          ! fortran-style I/O

Components:
Ent ESMF_Interface shared solvers giss_LSM dd2d

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
  ! atmosphere, land surface, and sea ice
AIC=AIC_CS32                     ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC_CS32.nc                  ! I.C. for land surface and sea ice
TOPO=Z_CS32_4X5                  ! topography
SOIL=SOIL_CS32                   ! soil types
VEG=V_CS32_144X90_1percent       ! veg. fractions
CROPS=CROPS_CS32_4X5             ! crops history
CDN=CD_CS32                      ! surf.drag coefficient
REG=REG.txt                      ! special regions-diag
RVR=RDdistocean_CS32_4X5.bin     ! river directions
TOP_INDEX=top_index_CS32         ! only used if #define do_topmodel_runoff

  ! ocean
OIC=OIC4X5LD.Z12.gas1.CLEV94.DEC01   ! ocean initial conditions
OFTAB=OFTABLE_NEW                    ! ocean function table
AVR=AVR72X46.L13.gas1.modelE         ! ocean filter
KBASIN=KB4X513.OCN.gas1              ! ocean basin designations
TOPO_OC=Z72X46N_gas.1_nocasp         ! ocean bdy.cond
GLMELT=GLMELT_CS32                   ! glacial melt distribution (on atm grid)
REMAP=remap72-46C32-32.nc            ! weights for atm-ocean coupling

  ! resolution independent files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k.1a              ! rad.tables and history files
RADN4=LWTables33k.1b              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_Ma_2008
RADN3=miescatpar.abcdv2
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr      ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files + 1 trend file)
O3file_01=mar2004_o3_shindelltrop_72x46x49x12_1850
O3file_02=mar2004_o3_shindelltrop_72x46x49x12_1890
O3file_03=mar2004_o3_shindelltrop_72x46x49x12_1910
O3file_04=mar2004_o3_shindelltrop_72x46x49x12_1930
O3file_05=mar2004_o3_shindelltrop_72x46x49x12_1950
O3file_06=mar2004_o3_shindelltrop_72x46x49x12_1960
O3file_07=mar2004_o3_shindelltrop_72x46x49x12_1970
O3file_08=mar2005_o3_shindelltrop_72x46x49x12_1980
O3file_09=mar2005_o3_shindelltrop_72x46x49x12_1990
O3trend=mar2005_o3timetrend_46x49x2412_1850_2050
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
MSU_wts=MSU.RSS.weights.data

Label and Namelist:
E4CubeTest (cubed-sphere version of E4M20 coupled to 4x5 13-layer ocean)

&&PARAMETERS
! parameters set for coupled ocean runs:
KOCEAN=1        ! ocn is prognostic

! parameters usually not changed when switching to coupled ocean:

variable_lk=1 ! let lakes grow or shrink in horizontal extent
wsn_max=2.   ! restrict snow depth to 2 m-h2o (if 0. snow depth is NOT restricted)

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_SDRAG=1         ! conserve ang. mom.

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

UOdrag=1 ! for flux-coupler verifications, atmosphere sees nonzero ocean/seaice velocity

U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.59      ! U00ice up => nethtz0 down (alb down); goals: nethtz0=0,plan.alb=30%
U00wtrX=1.39    ! U00wtrX+.01=>nethtz0+.7                      for global annual mean
!?1979 U00wtrX=1.38
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

madaer=3        ! updated aerosols
! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1850 ! if -1, crops in VEG-file is used
s0_yr=1850
s0_day=182
ghg_yr=1850
ghg_day=182
volc_yr=-1  ! 1850-1999 mean strat.aeros
volc_day=182
aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (was .0036)
albsn_yr=1850
dalbsnX=.024
o3_yr=1850

calc_orb_par=1
paleo_orb_yr=100.  !  BP i.e. 1950-paleo_orb_yr AD = 1850 AD

DTsrc=1800.      ! physics timestep
DT=1800.         ! for FV dynamics, set same as DTsrc

! parameters that may have to be changed in emergencies:
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=480       ! use =48 except on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time if isccp-diags are not essential
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc

&&END_PARAMETERS

 &INPUTZ
   YEARI=1901,MONTHI=1,DATEI=1,HOURI=0, !  from default: IYEAR1=YEARI
   YEARE=1901,MONTHE=1,DATEE=2,HOURE=0, KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1
 &END
