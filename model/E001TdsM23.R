E001TdsM23.R GISS Model E                                 gas 06/00

WARNING: The boundary conditions used here may not be what you want
         and no tuning has yet been done.
  Please check and see before running
E001TdsM23: modelE1 (3.0) strat. version with Shindell tropospheric chemistry
Balanced for 1950 SST + forcings
 T1=Ox
 T2=NOx
 T3=N2O5
 T4=HNO3
 T5=H2O2
 T6=CH3OOH
 T7=HCHO
 T8=HO2NO2
 T9=CO
 T10=CH4
 T11=PAN
 T12=Isoprene
 T13=AlkylNit
 T14=Alkenes
 T15=Paraffin

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! tracers can interact with water
#define TRACERS_DRYDEP              ! include tracer dry deposition
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
!  OFF #define EDGAR_HYDE_SOURCES          ! use EDGAR-HYDE tracers sources instead
!  OFF #define SHINDELL_STRAT_CHEM         ! to turn on stratospheric chemistry
!  OFF #define regional_Ox_tracers         ! to turn on regional Ox tracers
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M23                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
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
TRCHEM_fastj                        ! tracer chem photlysis code/rad transf
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
! ----------------------------------
TRDRYDEP                            ! tracer dry deposition from Harvard CTM
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
CLOUDS2 CLOUDS2_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB                               ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DUM ! or: ICEDYN ICEDYN_DRV  ! dynamic ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
AIC=AIC.RES_M23.D771201
GIC=GIC.E046D3M20A.1DEC1955
! AIC=DEC1981.rsfB396M23_AIC_2.3.3
! OHT=OTSPEC.RB399AM12.M250D   !  not needed if KOCEAN=0
OCNML=Z1O.B4X5.cor ! needed for post-processing only
OSST=OST4X5.B.1990-99avg.Hadl1.1
SICE=SICE4X5.B.1990-99avg.Hadl1.1
CDN=CD4X500S VEG=V72X46.1.cor2
SOIL=S4X50093 TOPO=Z72X46N.cor4_nocasp ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR.2      ! river direction file
ZVAR=ZVAR4X5         ! topographic variation for gwdrag
RADN1=sgpgxg.table8    ! rad.tables
RADN2=radfil33k                   !     8/2003 version
RADN3=miescatpar.abcdv2
! RADN4,RADN5,RADNA,RADNB are no longer used
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux    ! need KSOLAR=2
RADNE=topcld.trscat8
! new ozone files (minimum 1, maximum 9 files)
O3file_01=aug2003_o3_shindelltrop_72x46x49x12_1850
O3file_02=aug2003_o3_shindelltrop_72x46x49x12_1890
O3file_03=aug2003_o3_shindelltrop_72x46x49x12_1910
O3file_04=aug2003_o3_shindelltrop_72x46x49x12_1930
O3file_05=aug2003_o3_shindelltrop_72x46x49x12_1950
O3file_06=aug2003_o3_shindelltrop_72x46x49x12_1960
O3file_07=aug2003_o3_shindelltrop_72x46x49x12_1970
O3file_08=aug2003_o3_shindelltrop_72x46x49x12_1980
O3file_09=aug2003_o3_shindelltrop_72x46x49x12_1990
O3trend=aug2003_o3timetrend_46x49x2412_1850_2050
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4_monthly
TOP_INDEX=top_index_72x46.ij.ext
!-----------------------------------------------
MOLEC=chem_files/ds3ch4_moleculesE
JPLRX=chem_files/gs_jpl00_trop_15_fix
JPLPH=chem_files/ds_photlist_trop_15
RATJ=chem_files/ratj.giss_15
SPECFJ=chem_files/jv_spec00_15.dat
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
! next one needed only if correct_strat_Ox=.true.
Ox_corr=gsin/corrOx_modelE_v4
CO_INDUSTRIAL=CO_sources/co_industrial
CO_BIOMASS=CO_sources/co_biomass_new
Alkenes_INDUSTRIAL=gsin/alkenes_industrial_4x5   !monthly  KG/hr/grid
Alkenes_BIOMASS=gsin/alkenes_biomass_4x5         !monthly  KG/hr/grid
Alkenes_VEGETATION=gsin/alkenes_vegetation_4x5   !monthly  KG/hr/grid
Paraffin_INDUSTRIAL=gsin/paraffin_industrial_4x5 !monthly  KG/hr/grid
Paraffin_BIOMASS=gsin/paraffin_biomass_4x5       !monthly  KG/hr/grid
Paraffin_VEGETATION=gsin/paraffin_vegetation_4x5 !monthly  KG/hr/grid
Isoprene_VEGETATION=gsin/isoprene_vegetation_4x5 !monthly  KG/hr/grid
NOx_FOSSIL_FUEL=NOy_sources/fossil_fuel_4x5
NOx_BIOMASS=NOy_sources/bio_burning_4x5
NOx_SOIL=NOy_sources/soil_nox_4x5
NOx_AIRCRAFT=NOy_sources/aircraft_4x5
SULFATE_SA=NOy_sinks/sulfate_fakeM23_M_SA
!-----------------------------------------------
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

Label and Namelist:
E001TdsM23 (modelE1 (3.0) strat version with tropospheric chemsitry - 1950 conditions)
R=00BG/B
DTFIX=180
&&PARAMETERS
X_SDRAG=.00025,.000025  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.1     ! lin. sdrag above .1mb (top 2 layers) except near poles
PP_SDRAG=1.    ! lin. sdrag above 1.mb near poles (top 4 layers)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000035 !the default is 15d-6
PCONPEN=500.   ! penetrating convection defn for GWDRAG
CMC = 0.0000003

KOCEAN=0
U00ice  = .61   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
u00wtrx = 1.259
RWCLDOX= 1.5  !  wtr cld particle size *3/2 over ocean
RICLDX = 0.3333!  ice cld particle size * 1(at 0mb)->1/3(at 1000mb)

H2ObyCH4=1.  ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the Shapiro filter
DT_XUfilter=300. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=300. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DT=300.,             ! from default: DTsrc=3600.,
NIsurf=3,            ! number of surface time steps

NSUBDD=0        ! saving sub-daily diags 12hrly
Kvflxo=0        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=0

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=-1 ! if -1, crops in VEG-file is used
s0_yr=1950
ghg_yr=1950
ghg_day=182
s0_day=182
volc_yr=1950
volc_day=182
aero_yr=1950
o3_yr=1950
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0,  !  from default: IYEAR1=YEARI
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,  KDIAG=0,2,2,9*0,9,
   YEARE=1950,MONTHE=2,
   ISTART=2,IRANDI=0, YEARE=1950,MONTHE=1,DATEE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
