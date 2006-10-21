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
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INITIAL_GHG_SETUP        ! only for setup hour to get ghg IC file 
!  OFF #define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
!  OFF #define EDGAR_HYDE_SOURCES       ! use EDGAR-HYDE tracers sources instead
!  OFF #define regional_Ox_tracers      ! turns on regional Ox tracers   
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
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
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
! TRCHEM_fastj (not ready for MPI yet)! tracer chem photlysis code/rad transf
TRCHEM_fastj2                       ! used for stratosphere chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
! TRACERS_AEROSOLS_Koch_e4
! COSMO_SOURCES   
! ----------------------------------
TRDRYDEP                            ! tracer dry deposition from Harvard CTM
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
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
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION_E1        ! radiation modules
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
GHG=GHG.1850-2050.Mar2002
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
CO_INDUSTRIAL=CO_sources/co_industrial
CO_BIOMASS=CO_sources/co_biomass_new
Alkenes_INDUSTRIAL=gsin/alkenes_industrial_4x5   !monthly  KG/hr/grid
Alkenes_BIOMASS=gsin/alkenes_biomass_4x5         !monthly  KG/hr/grid
Alkenes_VEGETATION=gsin/alkenes_vegetation_4x5   !monthly  KG/hr/grid
Paraffin_INDUSTRIAL=gsin/paraffin_industrial_4x5 !monthly  KG/hr/grid
Paraffin_BIOMASS=gsin/paraffin_biomass_4x5       !monthly  KG/hr/grid
Paraffin_VEGETATION=gsin/paraffin_vegetation_4x5 !monthly  KG/hr/grid
NOx_FOSSIL_FUEL=NOy_sources/fossil_fuel_4x5
NOx_BIOMASS=NOy_sources/bio_burning_4x5
NOx_SOIL=NOy_sources/soil_nox_4x5
NOx_AIRCRAFT=NOy_sources/aircraft_4x5
CH4_ANIMALS=methane/gcm_data/CH4ANIMLS_4X5    ! Annual 1.3847
CH4_COALMINE=methane/gcm_data/CH4COAL_4X5      ! Annual 1.0285
CH4_GASLEAK=methane/gcm_data/CH4GASLEAK_4X5   ! Annual 3.904
CH4_GASVENT=methane/gcm_data/CH4GASVENT_4X5   ! Annual 1.659
CH4_CITYDUMP=methane/gcm_data/CH4MSW_4X5       ! Annual 1.233
CH4_SOIL_ABS=methane/gcm_data/CH4SOILABS_4X5   ! Annual 1.194
CH4_TERMITES=methane/gcm_data/CH4TRMITE_4X5    ! Annual 0.999
CH4_COALBURN=methane/gcm_data/COAL_BURN_BY_POP84_4X5  ! Annual 7.2154
CH4_BURN=methane/gcm_data/CH4BURN_4X5      ! Monthly 0.4369
CH4_RICE=methane/gcm_data/CH4RICEC_4X5     ! Monthly 0.7533
CH4_WETL=methane/gcm_data/CH4WETL+TUNDRA_4X5  ! Monthly 0.9818; zonal also
Isoprene_VEGETATION=gsin/isoprene_vegetation_4x5 !monthly  KG/hr/grid
SULFATE_SA=NOy_sinks/sulfate_fakeM23_M_SA
DMS_FIELD=dms_conc
SO2_FIELD=so2_conc
! --------Dorothy's------------------------------
! DMS_SEA=DMS.dat
! AER_CHEM=trace_gas_3D_fields_E
! AER_NO3=NO3_12
! SO2_IND=SO2_1990.EDGAR3.2
! SO2_BIOMASS=bioburn
! SO2_VOLCANO=SO2_volc_conti2000.AEROCOM_FEB12
! AIRCRAFT=MM_fuel_2015_subsonic_M31
! BC_BIOFUEL=BC_bf2000.AEROCOM_DEC03
! BC_FOSSIL_FUEL=BC_ff2000.AEROCOM_DEC03
! BC_BIOMASS=carb_biom
! OC_BIOFUEL=OC_bf2000.AEROCOM_DEC03
! OC_FOSSIL_FUEL=OC_ff2000.AEROCOM_DEC03
! OC_BIOMASS=OC_fire2000.AEROCOM_DEC03
! TERPENE=terp
! BE7_COSMO=Be23m1phi700.dat 

Label and Namelist:
E001TdsM23 (sample rundeck with Shindell chemistry tracers)
R=00BG/B
DTFIX=300
&&PARAMETERS

ocn_cycl=1      ! =0 if ocean varies from year to year
X_SDRAG=.00025,.000025  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.      ! no constant sdrag
P_SDRAG=.1      ! lin. sdrag above .1mb (top 2 layers) except near poles
PP_SDRAG=.1     ! lin. sdrag above 1.mb near poles (top 4 layers)
ANG_SDRAG=1     ! if =1: sdrag conserves ang mom.
PBREAK = 200.   ! The level for GW breaking above.
DEFTHRESH=0.000035 !the default is 15d-6
PCONPEN=500.    ! penetrating convection defn for GWDRAG
CMC = 0.0000003
KOCEAN=0
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
LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DTsrc = 1800.        ! half-hour physics time step (default: DTsrc=3600.)
DT=450.,             ! dynamics time step
NIsurf=1,            ! number of surface time steps
Kvflxo=0             ! saving VFLXO (daily)
Ndisk=480            ! i.e. 10 days on halem
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
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
COUPLED_CHEM=0     ! to couple chemistry and aerosols
imAER=0
rad_interact_tr=1  ! 1=use calculated Ox in radiation, 0=use climatology
                   ! (either case does the rad-forcing calculation)
rad_forc_lev=1     ! use LTROPO(I,J) level for rad forcing diags.
use_rad_n2o=0      ! use the radiation code's N2O (set ghg_yr=0)
use_rad_cfc=0      ! use rad code cfc11+cfc12, adjusted  (set ghg_yr=0)
use_rad_ch4=0      ! use rad code CH4, shut off sfc sources (set ghg_yr=0)
rad_FL=0           ! use rad code insolation getting fastj2 photon flux (s0_yr=0)
prather_limits=1   ! to avoid some negative tracers in sub-gridscale
which_trop=0       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
fix_CH4_chemistry=-1   ! for setting fixed methane value for chemistry:
pfix_CH4_S=1.750d-6    ! Southern Hemisphere
pfix_CH4_N=1.855d-6    ! Northern Hemisphere
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)

! To run a preindustrial case, alter the sulfate surface area, SST, & seaice
! files and the radiation years in this rundeck. Also, fix the ch4 chemistry
! to the appropriate value (0.73ppmv in both hemispheres?) Then, set PI_run=1
! and choose the various ratios for altering initial conditions, stratospheric
! overwriting, and sources of the tracers below:
PI_run        =    0    ! 0 =no, 1=yes for running pre-industrial cases
PIratio_N     = 0.667d0 ! {NOx, HNO3, N2O5, HO2NO2}
PIratio_CO_T  = 0.333d0 ! CO in troposphere
PIratio_CO_S  = 0.500d0 ! CO in stratosphere
PIratio_other = 0.500d0 ! {PAN,Isoprene,AlkyNit,Alkenes,Paraffin}
PIratio_indus = 0.000d0 ! factor for industrial sources
PIratio_bburn = 0.100d0 ! factor for biomass burning sources
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

