E001TdsM23.R GISS Model E                                gas 06/00

modelE1 (3.0) with Drew Shindell tropospheric and stratospheric
chemistry. (based upon E012TdsSjM23.R)

************************************************************
Please note this rundeck is provided to show how to run with
chemistry tracers on. You should check on the GCM tuning and
boundary conditions before using in production simulations.
************************************************************

You might need to increase your stacksize to get it to run in
parallel. Try 'unlimit stacksize' or increase to the proper amount.

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! tracers can interact with water
#define TRACERS_DRYDEP              ! include tracer dry deposition
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
#define SHINDELL_STRAT_CHEM         ! turns on stratospheric chemistry
#define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INITIAL_GHG_SETUP        ! only for setup hour to get ghg IC file
!  OFF #define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
!  OFF #define EDGAR_HYDE_SOURCES       ! use EDGAR-HYDE tracers sources instead
!  OFF #define regional_Ox_tracers      ! turns on regional Ox tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
!  OFF #define GFED_3D_BIOMASS          ! turns on IIASA AR4 GFED biomass burning
!  OFF #define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions

End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M23                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
DOMAIN_DECOMP ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! common and driver for tracers
TRDRYDEP                            ! tracer dry deposition from Harvard CTM
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
! TRCHEM_fastj ! troposphere-only version of tracer chem photlysis code/rad transf
TRCHEM_fastj2                       ! used for trop+strat chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
! TRACERS_AEROSOLS_Koch_e4
! COSMO_SOURCES
! BIOGENIC_EMISSIONS                  ! old N.Unger interactive isoprene emissions
! ----------------------------------
CLOUDS2_E1 CLOUDS2_DRV CLOUDS_COM   ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL_E1              ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB_E1                            ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN ICEDYN_DRV                   ! dynamic ice modules
SparseCommunicator_mod              ! sparse gather/scatter module
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV_E1 RADIATION_E1        ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_RES_M                          ! diagnostics (resolution dependent)
CONST FFT72 UTILDBL SYSTEM          ! utilities
! NUDGE                             ! S. Bauer's code to nudge meteorology
POUT !_netcdf                       ! post-processing output

Data input files:
AIC=AIC.RES_M23.D771201
GIC=GIC.E046D3M20A.1DEC1955
OCNML=Z1O.B4X5.cor ! needed for post-processing only
OSST=OST4X5.B.1975-84avg.Hadl1.1
SICE=SICE4X5.B.1975-84avg.Hadl1.1
CDN=CD4X500S
VEG=V72X46.1.cor2_no_crops
CROPS=CROPS_72X46N.cor4
SOIL=S4X50093 TOPO=Z72X46N.cor4_nocasp ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD_modelE_M.RVR    ! river direction file
ZVAR=ZVAR4X5         ! topographic variation for gwdrag
RADN1=sgpgxg.table8  ! rad.tables
RADN2=radfil33k      ! 8/2003 version
RADN3=miescatpar.abcdv2
MSU_wts=MSU.RSS.weights.data
TAero_PRE=sep2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux    ! need KSOLAR=2
! RADN9=Solar_spectrum.1500-2004_gsf ! need KSOLAR=9
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files)
! warning: I know there are newer ones of these:
O3file_01=jan2004_o3_shindelltrop_72x46x49x12_1850
O3file_02=jan2004_o3_shindelltrop_72x46x49x12_1890
O3file_03=jan2004_o3_shindelltrop_72x46x49x12_1910
O3file_04=jan2004_o3_shindelltrop_72x46x49x12_1930
O3file_05=jan2004_o3_shindelltrop_72x46x49x12_1950
O3file_06=jan2004_o3_shindelltrop_72x46x49x12_1960
O3file_07=jan2004_o3_shindelltrop_72x46x49x12_1970
O3file_08=jan2004_o3_shindelltrop_72x46x49x12_1980
O3file_09=jan2004_o3_shindelltrop_72x46x49x12_1990
O3trend=jan2004_o3timetrend_46x49x2412_1850_2050
delta_O3=do3_shindell_72x46x49x12_2037-2107_E009TdsSxHp-SnHpM23
GHG=GHG_A1B.June2004.txt
! e.g. GHGic=E009TdsG1701FM23/GHG_IC_1701
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46.ij.ext
BC_dep=BC.Dry+Wet.depositions.ann
!-----------------------------------------------
! choose these for trop-only chem model
!-----------------------------------------------
!  MOLEC=chem_files/ds3ch4_moleculesE
!  JPLRX=chem_files/gs_jpl00_trop_15_fix
!  JPLPH=chem_files/ds_photlist_trop_15
!  RATJ=chem_files/ratj.giss_15
!  SPECFJ=chem_files/jv_spec00_15.dat
!-----------------------------------------------
! choose these for strat+trop chem model
!-----------------------------------------------
MOLEC=chem_files/ds4_moleculesE
JPLRX=chem_files/ds4_jpl00_T25_SSS
JPLPH=chem_files/ds4_photlist_T25
RATJ=chem_files/ratj.giss_25
SPECFJ=chem_files/jv_spec00_25.dat
N2O_IC=gsin/N2O_IC_M23_4x5_6.17
CFC_IC=gsin/CFC_IC_M23_4x5_6.17
CH4_IC=gsin/CH4_IC_M23_4x5_6.17
!-----------------------------------------------
ATMFJ=chem_files/jv_atms.dat
DRYCOEFF=chem_files/drydep.coef
VEGTYPE=chem_files/vegtype.global
OLSON=chem_files/drydep.table
LAI01=chem_files/lai01.global
LAI02=chem_files/lai02.global
LAI03=chem_files/lai03.global
LAI04=chem_files/lai04.global
LAI05=chem_files/lai05.global
LAI06=chem_files/lai06.global
LAI07=chem_files/lai07.global
LAI08=chem_files/lai08.global
LAI09=chem_files/lai09.global
LAI10=chem_files/lai10.global
LAI11=chem_files/lai11.global
LAI12=chem_files/lai12.global
Ox_IC=gsin/Ox_init_cond_M23_4x5 !see README in /usr/people/cmrun/gsin
! fltran file used if rad_FL.ne.0:
! FLTRAN=chem_files/Solar_spectrum.1500-2004_fastj2 ! KSOLAR=9
! FLTRAN=chem_files/solar.lean02.ann.uvflux_fastj2  ! KSOLAR=2
! next one needed only if correct_strat_Ox=.true.
Ox_corr=gsin/corrOx_modelE_v4
!----------Default emissions case (giss/geia)------------------------
CO_01=CO_sources/CO_GEIA_industrial_head
CO_02=CO_sources/CO_GEIA_biomass_burning_head
Alkenes_01=gsin/Alkenes_GEIA_industrial_head
Alkenes_02=gsin/Alkenes_GEIA_biomass_burning_head
Alkenes_03=gsin/Alkenes_GEIA_vegetation_head
Paraffin_01=gsin/Paraffin_GEIA_industrial_head
Paraffin_02=gsin/Paraffin_GEIA_biomass_burning_head
Paraffin_03=gsin/Paraffin_GEIA_vegetation_head
NOx_01=NOy_sources/NOx_GEIA_fossil_fuels_head
NOx_02=NOy_sources/NOx_GEIA_biomass_burning_head
NOx_03=NOy_sources/NOx_GEIA_soil_head
NOx_AIRCRAFT=NOy_sources/aircraft_4x5
CH4_01=methane/gcm_data/CH4_GEIA_Animals_header
CH4_02=methane/gcm_data/CH4_GEIA_Coal_Mining_header
CH4_03=methane/gcm_data/CH4_GEIA_Gas_Leak_header
CH4_04=methane/gcm_data/CH4_GEIA_Gas_Vent_header
CH4_05=methane/gcm_data/CH4_GEIA_Landfill_header
CH4_06=methane/gcm_data/CH4_GEIA_Soil_Absorption_header
CH4_07=methane/gcm_data/CH4_GEIA_Termites_header
CH4_08=methane/gcm_data/CH4_GEIA_Coal_Burning_header
CH4_09=methane/gcm_data/CH4_GEIA_Biomass_Burning_header
CH4_10=methane/gcm_data/CH4_GEIA_Rice_header
CH4_11=methane/gcm_data/CH4_GEIA_Wetlands_and_Tundra_header
Isoprene_01=gsin/Isoprene_GEIA_vegetation_head
SULFATE_SA=NOy_sinks/sulfate_fakeM23_M_SA
DMS_FIELD=dms_conc
SO2_FIELD=so2_conc
! ----- for interactive wetlands -----
PREC_NCEP=gsin/mean_prec_ncep_4x5
TEMP_NCEP=gsin/mean_temp_ncep_4x5
BETA_NCEP=gsin/beta_p_ch4_4x5
ALPHA_NCEP=gsin/alpha_t_ch4_4x5
!-----------------------------------------------

Label and Namelist:
E001TdsM23 (sample rundeck with Shindell chemistry tracers)
R=00BG/B
DTFIX=300
&&PARAMETERS

!--- define emission sectors above files belong to ---
CO_01_sect='CO FFUEL'
CO_02_sect='CO BBURN'
Alkenes_01_sect='ALK FFUEL'
Alkenes_02_sect='ALK BBURN'
Alkenes_03_sect='ALK VEG'
Paraffin_01_sect='PAR FFUEL'
Paraffin_02_sect='PAR BBURN'
Paraffin_03_sect='PAR VEG'
NOx_01_sect='NOX FFUEL'
NOx_02_sect='NOX BBURN'
NOx_03_sect='NOX'
NOx_AIRCRAFT_sect='NOX FFUEL' ! special 3D source case
CH4_01_sect='CH4'
CH4_02_sect='CH4 FFUEL'
CH4_03_sect='CH4 FFUEL'
CH4_04_sect='CH4 FFUEL'
CH4_05_sect='CH4'
CH4_06_sect='CH4'
CH4_07_sect='CH4'
CH4_08_sect='CH4 FFUEL'
CH4_09_sect='CH4 BBURN'
CH4_10_sect='CH4'
CH4_11_sect='CH4 WETL'
Isoprene_01_sect='ISO VEG'
!      (careful; they're allowed to overlap):
!       ---------define-REGIONS------------
!        global S.Asia E.Asia Europe N.Amer
REG_S=    -90.,    5.,   15.,   25.,   15.
REG_N=     90.,   35.,   50.,   65.,   55.
REG_W=   -180.,   50.,   95.,  -10., -125.
REG_E=    180.,   95.,  160.,   50.,  -60.
!       ---define-regions-names/order------
REGIONS_ARE='global S_Asia E_Asia Europe N_America'
!-fit-here--|                                                              |---
!       ---define-factors-by-sector--------
!        global S.Asia E.Asia Europe N.Amer
SECT_01= 1.000, 1.000, 1.000, 1.000, 1.000 ! CO
SECT_02= 1.000, 1.000, 1.000, 1.000, 1.000 ! ALK
SECT_03= 1.000, 1.000, 1.000, 1.000, 1.000 ! PAR
SECT_04= 1.000, 1.000, 1.000, 1.000, 1.000 ! CH4
SECT_05= 1.000, 1.000, 1.000, 1.000, 1.000 ! NOX
SECT_06= 1.000, 1.000, 1.000, 1.000, 1.000 ! ISO
SECT_07= 1.000, 1.000, 1.000, 1.000, 1.000 ! FFUEL
SECT_08= 1.000, 1.000, 1.000, 1.000, 1.000 ! BBURN
SECT_09= 1.000, 1.000, 1.000, 1.000, 1.000 ! WETL
SECT_10= 1.000, 1.000, 1.000, 1.000, 1.000 ! VEG
!       ---define-sectors-names/order------
SECTORS_ARE='CO ALK PAR CH4 NOX ISO FFUEL BBURN WETL VEG'
!-fit-here--|                                                              |---
!-----
aircraft_Tyr1=1990 ! for non-transient emissions,
aircraft_Tyr2=1990 ! set these two equal or omit them.
biomass_Tyr1= 1990 ! for non-transient emissions,
biomass_Tyr2= 1990 ! set these two equal or omit them.

! factor to tune the base isoprene emissions globally,
! only when #defined BIOGENIC_EMISSIONS, otherwise use
! a sector/region method above...
base_isopreneX=1.d0

! ---- for interactive wetlands -----
nn_or_zon=1     ! int dist method 1=zonal avg, 0=nearest neighbor
int_wet_dist=1  ! turn on(1)/off(0) interacive SPATIAL wetlands
ice_age=0.      ! if not 0 no wetl emis for lats poleward of +/- this in deg
ns_wet=11       ! index of CH4 source that is the wetlands (dumb, I know)
exclude_us_eu=0 ! to exclude (=1) the U.S. and E.U. from inter wetl dist
topo_lim=205.d0 ! upper limit of topographic variation for new wetlands 
sat_lim=-9.d0   ! lower limit on surf air temp for new wetlants
gw_ulim=100.d0  ! upper limit on ground wetness for new wetlands
gw_llim=18.d0   ! lower limit on ground wetness for new wetlands
SW_lim=27.d0    ! lower limit on SW downward flux for new wetlands
! -----------------------------------
X_SDRAG=.0005,.00005  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.01     ! lin. sdrag above p_sdrag mb (top layer for M23) except near poles
PP_SDRAG=4.6   ! lin. sdrag above  pp_sdrag mb near poles (top 5 layers for M23)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
WMAX=1000.     ! maximum wind velocity in sdrag; default=200 when GW drag not used
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000030 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=1.      ! Shear drag coefficient
CMTN=0.25      ! default is 0.5
CDEF=1.5       ! deformation drag coefficient
LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters

KOCEAN=0
ocn_cycl=1      ! =0 if ocean varies from year to year
U00ice  = .60   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
u00wtrx = 1.43  ! 1980 conditions (NEEDS TO BE VERIFIED!)
cond_scheme=2   ! more elaborate conduction scheme (GHY, Nancy Kiang)
H2ObyCH4=0.     ! (=1) activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2 ! 9
! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction
DTsrc = 1800.        ! half-hour physics time step (default: DTsrc=3600.)
DT=450.,             ! dynamics time step
NIsurf=1,            ! number of surface time steps
Kvflxo=0             ! saving VFLXO (daily)
Ndisk=480            ! e.g. every 10 days
NSUBDD=0             ! saving sub-daily diags :
! pick sub-daily frequency diags, next line ::
SUBDD=' '
SUBDD1=' '
LmaxSUBDD=18         ! only need to save up to level LmaxSUBDD
KCOPY=2              ! saving acc + rsf
isccp_diags=0        ! use =0 to save cpu time
! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=-1 ! if -1, crops in VEG-file is used
s0_yr   = 1979
s0_day  = 182
ghg_yr  = 1979
ghg_day = 182
volc_yr = 1984 !1979
volc_day=182
aero_yr =1979
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1979
dalbsnX=.015
o3_yr   =1979
aer_int_yr=1995
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
COUPLED_CHEM=0     ! to couple chemistry and aerosols
imAER=0
use_sol_Ox_cycle=0 ! (=1) apply ozone changes in radiation, based on solar cycle
rad_interact_tr=1  ! 1=use calculated Ox in radiation, 0=use climatology
                   ! (either case does the rad-forcing calculation)
rad_forc_lev=1     ! use LTROPO(I,J) level for rad forcing diags.
use_rad_n2o=0      ! use the radiation code's N2O 
use_rad_cfc=0      ! use rad code cfc11+cfc12, adjusted 
use_rad_ch4=0      ! use rad code CH4, shut off sfc sources 
rad_FL=0           ! use rad code insolation getting fastj2 photon flux
prather_limits=1   ! to avoid some negative tracers in sub-gridscale
which_trop=0       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
fix_CH4_chemistry=-1   ! for setting fixed methane value for chemistry:
pfix_CH4_S=1.750d-6    ! Southern Hemisphere (fix_CH4_chemistry=1)
pfix_CH4_N=1.855d-6    ! Northern Hemisphere (fix_CH4_chemistry=1)
ch4_init_sh=1.750      ! init cond, S.Hemi (fix_CH4_chemistry=0 ?)
ch4_init_nh=1.855      ! init cond, N.Hemi (fix_CH4_chemistry=0 ?)
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)

! To run a preindustrial case, alter the sulfate surface area, SST, & seaice
! files and the radiation years in this rundeck. Also, fix the ch4 chemistry
! to the appropriate value (0.73ppmv in both hemispheres?) Then, set PI_run=1
! and choose the various ratios for altering initial conditions, stratospheric
! overwriting. Probably best to alter sources via the sectors set in the
! rundeck instead, use these for initial conditions/overwriting:
PI_run        =    0    ! 0 =no, 1=yes for running pre-industrial cases
PIratio_N     = 0.667d0 ! {NOx, HNO3, N2O5, HO2NO2}
PIratio_CO_T  = 0.333d0 ! CO in troposphere
PIratio_CO_S  = 0.500d0 ! CO in stratosphere
PIratio_other = 0.500d0 ! {PAN,Isoprene,AlkyNit,Alkenes,Paraffin}
PIratio_N2O   = 0.896d0 ! {N2O IC's and L=1 overwriting}
PIratio_CFC   = 0.000d0 ! {CFC IC's and L=1 overwriting}
!----------------------------------------------------------------------
!----------------------------------------------------------------------
&&END_PARAMETERS

 &INPUTZ
   QCHECK=.false.
   kdiag = 0,0,0,0,0,0,0,0,0,0,0,0,0,
   YEARI=2000,MONTHI=1,DATEI=1,HOURI=0,
   YEARE=2001,MONTHE=1,DATEE=1,HOURE=1,
   ISTART=2,IRANDI=0, YEARE=2000, MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END

