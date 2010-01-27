E4TcadF40.R GISS Model E  1850 ocn/atm   jan perlwitz  07/2009
 Template for setting up simulations using the E4F40 model with dust,chemistry,aerosol tracers

E4TcadF40: Chemistry, aerosols and dust
E4TdusF40: E4F40 + dust tracers

E4F40: modelE as frozen (or not yet) in July 2009
modelE4 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 2000
ocean data: prescribed, 1996-2005 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define TRAC_ADV_CPU
#define USE_ENT                  ! include dynamic vegetation model
#define ROUGHL_HACK
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define NO_HDIURN                ! exclude hdiurn diagnostics
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
#define SHINDELL_STRAT_CHEM         ! turns on stratospheric chemistry
#define RAD_O3_GCM_HRES     ! Use GCM horiz resl to input rad code clim Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
#define INITIAL_GHG_SETUP        ! only for setup hour to get ghg IC file
#define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
!#define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
#define TRACERS_NITRATE
#define TRACERS_HETCHEM
#define BC_ALB                      !optional tracer BC affects snow albedo
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define CALCULATE_LIGHTNING ! turn on Colin Price lightning when TRACERS_SPECIAL_Shindell off
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
!  OFF #define GFED_3D_BIOMASS          ! turns on IIASA AR4 GFED biomass burning
!  OFF #define HTAP_LIKE_DIAGS    ! adds many diags, changes OH diag, adds Air tracer
!  OFF #define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_stratF40                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
QUS3D                               ! advection of Q and tracer gases
TRACER_COM TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRDRYDEP                            ! dry deposition of tracers
TRDUST_COM TRDUST                   ! dust tracer specific code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
TRCHEM_fastj2                       ! used for trop+strat chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
BIOGENIC_EMISSIONS                  ! old N.Unger interactive isoprene
! ----------------------------------
TRACERS_AEROSOLS_Koch_e4            ! BC/OC/sulfate/seasalt
!TRACERS_AEROSOLS_SOA                ! Secondary Organic Aerosols
TRACER_NITRATE                      ! Nitrate aerosol
TRAMP_eqsam_v03d                    ! EQSAM module for inorganic aerosol thermodynamic equilibrium
TRACER_HETCHEM
! ----------------------------------
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV ! + component giss_LSM: land surface and soils
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + component Ent: new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
lightning                           ! Colin Price lightning model
! flammability_drv flammability       ! Olga's fire model
ATURB_E1                            ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
RAD_native_O3                       ! for reading ozone to rad code at native GCM horiz res.
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_F                          ! diagnostics (resolution dependent)
FFT144                              ! utilities
POUT_netcdf                         ! post-processing output

Components:
Ent shared ESMF_Interface solvers giss_LSM

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
    ! resolution dependent files
    ! start up from the restart file of an earlier run ...
! AIC=1....rsfE... ! initial conditions, no GIC needed, use ISTART=8
    ! ... or from observed conditions AIC and model ground data GIC
AIC=AIC.RES_F40.D771201  ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1   ! initial ground conditions      ISTART=2
OSST=OST_144x90.1876-1885avg.HadISST1.1    ! prescr. climatological ocean (1 yr data)
SICE=SICE_144x90.1876-1885avg.HadISST1.1   ! prescr. climatological sea ice
! For q-flux ocean, replace all the above by the next 2 lines, set KOCEAN=1, ISTART=8
!! AIC=1JAN1961.rsfE4F40.MXL65m        ! end of preliminary run with KOCEAN=0
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960   ! ocean horizontal heat transports
OCNML=Z1O.B144x90                    ! mixed layer depth (not used if KOCEAN=0)
CDN=CD144X90.ext
VEG=V144X90_no_crops.ext
CROPS=CROPS2007_144X90N_nocasp
SOIL=S144X900098M.ext
TOPO=Z144X90N_nocasp              ! bdy.cond
REG=REG2X2.5                      ! special regions-diag
RVR=RD_modelE_Fa.RVR.bin          ! river direction file
TOP_INDEX=top_index_144x90_a.ij.ext
ZVAR=ZVAR2X25A             ! topographic variation for gwdrag
GLMELT=GLMELT_144X90_gas.OCN   ! glacial melt distribution
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_nmaps_2x2.5bin.dat
VEG_DENSE=gsin/veg_dense_2x2.5 ! vegetation density for flammability calculation

    ! resolution independent files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k.1a              ! rad.tables and history files
RADN4=LWTables33k.1b              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux       ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
GHG=GHG.Mar2009.txt ! use GHG.Jul2009.txt for runs that start before 1850
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU.RSS.weights.data

! ozone files (minimum 1, maximum 12 files )
! this set for #defined RAD_O3_GCM_HRES
O3file_01=jan2010_o3_shindell_144x90x49x12_1850
O3file_02=jan2010_o3_shindell_144x90x49x12_1870
O3file_03=jan2010_o3_shindell_144x90x49x12_1890
O3file_04=jan2010_o3_shindell_144x90x49x12_1910
O3file_05=jan2010_o3_shindell_144x90x49x12_1930
O3file_06=jan2010_o3_shindell_144x90x49x12_1940
O3file_07=jan2010_o3_shindell_144x90x49x12_1950
O3file_08=jan2010_o3_shindell_144x90x49x12_1960
O3file_09=jan2010_o3_shindell_144x90x49x12_1970
O3file_10=jan2010_o3_shindell_144x90x49x12_1980
O3file_11=jan2010_o3_shindell_144x90x49x12_1990
O3file_12=jan2010_o3_shindell_144x90x49x12_2000
Ox_ref=jan2010_o3_shindell_144x90x49x12_April1850 ! for radiative forcing reference

!! ozone files (minimum 1, maximum 9 files + 1 trend file)
!! leaving commented files for NOT #defined RAD_O3_GCM_HRES:
!! O3file_01=mar2004_o3_shindelltrop_72x46x49x12_1850
!! O3file_02=mar2004_o3_shindelltrop_72x46x49x12_1890
!! O3file_03=mar2004_o3_shindelltrop_72x46x49x12_1910
!! O3file_04=mar2004_o3_shindelltrop_72x46x49x12_1930
!! O3file_05=mar2004_o3_shindelltrop_72x46x49x12_1950
!! O3file_06=mar2004_o3_shindelltrop_72x46x49x12_1960
!! O3file_07=mar2004_o3_shindelltrop_72x46x49x12_1970
!! O3file_08=mar2005_o3_shindelltrop_72x46x49x12_1980
!! O3file_09=mar2005_o3_shindelltrop_72x46x49x12_1990
!! O3trend=mar2005_o3timetrend_46x49x2412_1850_2050
!! Ox_ref=gsin/O3ref_O3JDAY_1850_182.dat

! sulfate+black carbon files:
!  MADAER=1 (default) needs:
!TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
!TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
!TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
!TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
BC_dep=BC.Dry+Wet.depositions.ann

!-----------------------------------------------
!  Start tracer code files:
!-----------------------------------------------
!  full-atmosphere chemistry model files:
!-----------------------------------------------
MOLEC=chem_files/ds4_moleculesE_terp!_soa
JPLRX=chem_files/jpl00_T25_SEP08_fastterp
JPLPH=chem_files/ds4_photlist_T25
RATJ=chem_files/ratj.giss_25
SPECFJ=chem_files/jv_spec00_25.dat
ATMFJ=chem_files/jv_atms.dat
N2O_IC=gsin/N2O_IC_M23_4x5_6.17_conc_2x2.5_conc
CFC_IC=gsin/CFC_IC_M23_4x5_6.17_conc_2x2.5_conc
CH4_IC=gsin/CH4_IC_M23_4x5_6.17_conc_2x2.5_conc
Ox_IC=gsin/Ox_init_cond_M23_4x5_conc_2x2.5_conc
CO_IC=gsin/CO_init_cond_M23_conc_2x2.5_conc
! fltran file used if rad_FL.ne.0:
! FLTRAN=chem_files/Solar_spectrum.1500-2004_fastj2 ! KSOLAR=9
! FLTRAN=chem_files/solar.lean02.ann.uvflux_fastj2  ! KSOLAR=2
!-----------------------------------------------
SULFATE_SA=temp_2x2.5/sulfate_pi_fakeM23_M_SA_2x2.5gf ! really 4x5 and 9-layer
DMS_FIELD=temp_2x2.5/dms_conc_2x2.5gf ! really 4x5
SO2_FIELD=temp_2x2.5/so2_conc_2x2.5gf ! really 4x5

! files for dust tracers
ERS=ERS1_1993_MONTHLY.144x90.threshold-13 ! ERS data
GIN=Ginoux2001_source_VegMask_144x90      ! preferred sources
LKTAB=log_dust_emission_60ms-1 ! look up table for emission calculations
LKTAB1=table_wspdf             ! look up table for wind speed probabilities

!------- Needed for dry deposition ---------
VEGTYPE=chem_files/vegtype.global_2x2.5gf ! really 4x5
OLSON=chem_files/drydep.table
DRYCOEFF=chem_files/drydep.coef
LAI01=chem_files/lai01.global_2x2.5gf ! really 4x5
LAI02=chem_files/lai02.global_2x2.5gf ! really 4x5
LAI03=chem_files/lai03.global_2x2.5gf ! really 4x5
LAI04=chem_files/lai04.global_2x2.5gf ! really 4x5
LAI05=chem_files/lai05.global_2x2.5gf ! really 4x5
LAI06=chem_files/lai06.global_2x2.5gf ! really 4x5
LAI07=chem_files/lai07.global_2x2.5gf ! really 4x5
LAI08=chem_files/lai08.global_2x2.5gf ! really 4x5
LAI09=chem_files/lai09.global_2x2.5gf ! really 4x5
LAI10=chem_files/lai10.global_2x2.5gf ! really 4x5
LAI11=chem_files/lai11.global_2x2.5gf ! really 4x5
LAI12=chem_files/lai12.global_2x2.5gf ! really 4x5

!---------- mostly transient, mostly AR5 gas tracer emissions ------------------
CO_01=AR5_emis/F/T/CO_ind_AR5_1850-2000_2x2.5_h
CO_02=AR5_emis/F/T/CO_tra_AR5_1850-2000_2x2.5_h
CO_03=AR5_emis/F/T/CO_wst_AR5_1850-2000_2x2.5_h
CO_04=AR5_emis/F/T/CO_awb_AR5_1850-2000_2x2.5_h
CO_05=AR5_emis/F/T/CO_dom_AR5_1850-2000_2x2.5_h
CO_06=AR5_emis/F/T/CO_forestfire_AR5_1900-2000_2x2.5_h
CO_07=AR5_emis/F/T/CO_grassfire_AR5_1900-2000_2x2.5_h
CO_08=AR5_emis/F/T/m_CO_shp_AR5_1850-2000_2x2.5_h
CO_09=AR5_emis/F/T/CO_slv_AR5_1990-2000_2x2.5_h
CO_10=AR5_emis/F/T/CO_ene_AR5_1850-2000_2x2.5_h
CO_11=AR5_emis/F/T/CO_agr_AR5_1990-2000_2x2.5_h
NOx_AIRC=AR5_emis/F/T/NOx_air_AR5_1910-2000_2x2.5
NOx_01=AR5_emis/F/NAT/NOx_Soil_GEIA_2x2.5_HALF_h ! half because we have ag source
NOx_02=AR5_emis/F/T/NOx_awb_AR5_1850-2000_2x2.5_h
NOx_03=AR5_emis/F/T/NOx_dom_AR5_1850-2000_2x2.5_h
NOx_04=AR5_emis/F/T/NOx_ene_AR5_1850-2000_2x2.5_h
NOx_05=AR5_emis/F/T/NOx_forestfire_AR5_1900-2000_2x2.5_h
NOx_06=AR5_emis/F/T/NOx_grassfire_AR5_1900-2000_2x2.5_h
NOx_07=AR5_emis/F/T/NOx_ind_AR5_1850-2000_2x2.5_h
NOx_08=AR5_emis/F/T/m_NOx_shp_AR5_1850-2000_2x2.5_h
NOx_09=AR5_emis/F/T/NOx_tra_AR5_1850-2000_2x2.5_h
NOx_10=AR5_emis/F/T/NOx_wst_AR5_1850-2000_2x2.5_h
NOx_11=AR5_emis/F/T/NOx_agr_AR5_1850-2000_2x2.5_h
! Note that the Isoprene emis file is ignored when BIOGENIC_EMISSIONS
! directive is on. But I am commenting it anyway.
! if BIOGENIC_EMISSIONS or PS_BVOC are defined, and only one
! Terpenes file is available, the Isoprene file is needed too
! Isoprene_01=ORCHIDEE_Isoprene_1990_2x2.5_h
Terpenes_01=ORCHIDEE_Terpenes_1990_2x2.5_h
Terpenes_02=ORCHIDEE_ORVOC_1990_2x2.5_h
! ========= please remember that Alkenes =================
! ========= and Paraffin emissions files =================
! ========= must now be in Kmole units,  =================
! ========= not Kg units ...             =================
Alkenes_01=AR5_emis/F/NAT/Alkenes_vegetation_GEIA_2x2.5_h_1
Alkenes_02=AR5_emis/F/T/m_Alkenes_shp_AR5_1850-2000_2x2.5_h
Alkenes_03=AR5_emis/F/T/Alkenes_wst_AR5_1850-2000_2x2.5_h
Alkenes_04=AR5_emis/F/T/Alkenes_dom_AR5_1850-2000_2x2.5_h
Alkenes_05=AR5_emis/F/T/Alkenes_forestfire_AR5_1900-2000_2x2.5_h
Alkenes_06=AR5_emis/F/T/Alkenes_grassfire_AR5_1900-2000_2x2.5_h
Alkenes_07=AR5_emis/F/T/Alkenes_ind_AR5_1850-2000_2x2.5_h
Alkenes_08=AR5_emis/F/T/Alkenes_tra_AR5_1850-2000_2x2.5_h
Alkenes_09=AR5_emis/F/T/Alkenes_ene_AR5_1850-2000_2x2.5_h
Alkenes_10=AR5_emis/F/T/Alkenes_awb_AR5_1890-2000_2x2.5_h
Alkenes_11=AR5_emis/F/T/Alkenes_agr_AR5_1990-2000_2x2.5_h
Paraffin_01=AR5_emis/F/NAT/Paraffin_vegetation_GEIA_2x2.5_h_1
Paraffin_02=AR5_emis/F/T/m_Paraffin_shp_AR5_1850-2000_2x2.5_h
Paraffin_03=AR5_emis/F/T/Paraffin_wst_AR5_1850-2000_2x2.5_h
Paraffin_04=AR5_emis/F/T/Paraffin_dom_AR5_1850-2000_2x2.5_h
Paraffin_05=AR5_emis/F/T/Paraffin_forestfire_AR5_1900-2000_2x2.5_h
Paraffin_06=AR5_emis/F/T/Paraffin_grassfire_AR5_1900-2000_2x2.5_h
Paraffin_07=AR5_emis/F/T/Paraffin_ind_AR5_1850-2000_2x2.5_h
Paraffin_08=AR5_emis/F/T/Paraffin_tra_AR5_1850-2000_2x2.5_h
Paraffin_09=AR5_emis/F/T/Paraffin_slv_AR5_1860-2000_2x2.5_h
Paraffin_10=AR5_emis/F/T/Paraffin_ene_AR5_1850-2000_2x2.5_h
Paraffin_11=AR5_emis/F/T/Paraffin_awb_AR5_1890-2000_2x2.5_h
Paraffin_12=AR5_emis/F/T/Paraffin_agr_AR5_1990-2000_2x2.5_h
codirect_01=AR5_emis/F/HTAP_codirect_emissions_2x2.5_h
!------------ end of chem emissions files ---------------

!-------Aerosol inputs, for imAER=5 --------------
!----oxidants needed if not coupled to chem ----
!---these are place-holders, right now must run coupled ----
!AER_CHEM=OXID_E__1TgfF40_2x2.5
!AER_OH_STRAT=Strat_OH_drewE_20
!O3_FIELD=Ox_3D_field_bell
!OFFLINE_HNO3.nc=HNO3_E70_GISS4x5.nc
!OFFLINE_SEAS.nc=SEASALT_EK1su_GISS4x5.nc
!------nitrate inputs--------------
NH3_CON_01=NH3hCON_OCEANflux_Jan10_2x2.5_h
NH3_CON_02=NH3_forestfire_AR5_1900-2000_2x2.5_h
NH3_CON_03=NH3_grassfire_AR5_1900-2000_2x2.5_h
NH3_CON_04=NH3_agr_AR5_1850-2000_2x2.5_h
NH3_CON_05=NH3_awb_AR5_1850-2000_2x2.5_h
NH3_CON_06=NH3_dom_AR5_1850-2000_2x2.5_h
NH3_CON_07=NH3_ind_AR5_1850-2000_2x2.5_h
NH3_CON_08=NH3_ene_AR5_1850-2000_2x2.5_h
NH3_CON_09=NH3_tra_AR5_1850-2000_2x2.5_h
! ------- aerosol, needed for imAER= 1 or 3 or 5-----------
SO2_VOLCANO=SO2_volc_conti2000_HR2x2.5.AEROCOM
! -------- aerosol, needed for imAER = 3 or 5 --------------
DMS_SEA=DMS_Kettle_Andreae_2x2.5
! The following line is only needed for imAER=1
! Note that if gas-phase chemistry is enabled, both
! TERPENE and Terpenes_01 are used, which differ significantly
!TERPENE=terp_Guenther_2x2.5
!SO2_AIRCRAFT=NOy_sources/aircraft_4x5_1940-2000 ! zero in 1940 and before.
!--------AR5 inputs for imAER=5 --------------------
BC_EM_1=BCII_awb_AR5_1850-2000_2x2.5_h
BC_EM_2=BCII_dom_AR5_1850-2000_2x2.5_h
BC_EM_3=BCII_ene_AR5_1850-2000_2x2.5_h
BC_EM_4=BCII_ind_AR5_1850-2000_2x2.5_h
BC_EM_5=BCII_tra_AR5_1850-2000_2x2.5_h
BC_EM_6=BCII_wst_AR5_1850-2000_2x2.5_h
BC_EM_7=BCII_shp_AR5_1850-2000_2x2.5_h
OC_EM_1=OCII_awb_AR5_1850-2000_2x2.5_h
OC_EM_2=OCII_dom_AR5_1850-2000_2x2.5_h
OC_EM_3=OCII_ene_AR5_1850-2000_2x2.5_h
OC_EM_4=OCII_ind_AR5_1850-2000_2x2.5_h
OC_EM_5=OCII_tra_AR5_1850-2000_2x2.5_h
OC_EM_6=OCII_wst_AR5_1850-2000_2x2.5_h
OC_EM_7=OCII_shp_AR5_1850-2000_2x2.5_h
SO2_EM_1=SO2_shp_AR5_1850-2000_2x2.5_h
SO2_EM_2=SO2_awb_AR5_1850-2000_2x2.5_h
SO2_EM_3=SO2_dom_AR5_1850-2000_2x2.5_h
SO2_EM_4=SO2_ind_AR5_1850-2000_2x2.5_h
SO2_EM_5=SO2_tra_AR5_1850-2000_2x2.5_h
SO2_EM_6=SO2_wst_AR5_1850-2000_2x2.5_h
SO2_EM_7=SO2_ene_AR5_1850-2000_2x2.5_h
BCB_EM_1=BCB_forestfire_AR5_1900-2000_2x2.5_h
BCB_EM_2=BCB_grassfire_AR5_1900-2000_2x2.5_h
OCB_EM_1=OCB_forestfire_AR5_1900-2000_2x2.5_h
OCB_EM_2=OCB_grassfire_AR5_1900-2000_2x2.5_h
SO2B_EM_1=SO2_forestfire_AR5_1900-2000_2x2.5_h
SO2B_EM_2=SO2_grassfire_AR5_1900-2000_2x2.5_h


Label and Namelist:  (next 2 lines)
E4TcadF40 (modelE4 2x2.5 hor., 40 lyrs, 1850 atm., 1996-2005 clim ocn, aer,chem,dust tracers)

&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! 0: ocean is prescribed 1: ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)
ocn_cycl=1      ! 0: ocean varies from year to year 1: ocean doesn't vary

variable_lk=1      ! 0: fixed lakes 1: variable lakes
roughl_from_file=0 ! 1: read roughness length from file

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000045 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.2       ! default is 0.5
CDEF=1.95      ! deformation drag coefficient
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.74 ! affects clouds above 850mb w/o MC region;  tune this first to get about 30% high cloud
U00b=1.65 ! affects clouds below 850mb and MC regions; tune this last  to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.57      ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.46    ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

H2ObyCH4=0.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
cloud_rad_forc=1 ! 1: calculate cloud radiative forcing

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1850  ! if -1, crops in VEG-file is used
s0_yr=1850
s0_day=182
ghg_yr=1850
ghg_day=182
volc_yr=-1 ! -1 means 150-yr mean 1850-1999
volc_day=182
aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
albsn_yr=1850
dalbsnX=.024
o3_yr=-1850
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
CO2X=1.
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

calc_orb_par=1
paleo_orb_yr=100.  !  BP i.e. 1950-paleo_orb_yr AD = 1850 AD

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=225.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=960
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags 0hrly
Nssw=2   ! 'til diurnal diags fixed,has to be even ; best to set=NSUBDD if not 0
KCOPY=2         ! saving acc + rsf
isccp_diags=0   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc

!-----------------------------------------------
!  Start tracer code parameters:
!-----------------------------------------------
!--- define emission sectors above files belong to ---
! example: CH4_13_sect='WET'

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
SECT_01= 1.000, 1.000, 1.000, 1.000, 1.000 ! WET (for example)
!       ---define-sectors-names/order------
SECTORS_ARE='WET'
!-fit-here--|                                                              |---
!-----
aircraft_Tyr1=1910 ! for non-transient emissions,
aircraft_Tyr2=2000 ! set these two equal or omit them.
! biomass_Tyr1= 0 ! for non-transient emissions,
! biomass_Tyr2= 0 ! set these two equal or omit them.

! Colin Price lightning model needs resolution-dependant tuning:
tune_lt_land=1.4322d0 ! =2.2d0*2.17d0 then *0.5*1.2*0.5 for 2x2.5 model
tune_lt_sea= 3.1736d0 ! =3.9d0*2.17d0 then *0.25*1.5 for 2x2.5 model

! -----------------------------------
! Pressure above which Ox, NOx, BrOx, and ClOx will be
! overwritten with climatology. Default of 0.1 is (From 23-
! layer model) won't choose any 40-layer model levels):
PltOx=0.2

Tpsc_offset_N=0.d0 ! pol.strat.cloud temperature offset NH
Tpsc_offset_S=0.d0 ! pol.strat.cloud temperature offset SH

COUPLED_CHEM=1     ! to couple chemistry and aerosols
use_sol_Ox_cycle=0 ! (=1) apply ozone changes in radiation, based on solar cycle
rad_interact_chem=1! 1=use calculated Ox/CH4 in radiation, 0=use climatology
                   ! (either case does the rad-forcing calculation)
! Lmax_rad_O3=0    ! Ox levels used in rad code default is LM
! Lmax_rad_CH4=0   ! CH4 levels used in rad code default is LM
use_rad_n2o=1      ! use the radiation code's N2O
use_rad_cfc=1      ! use rad code cfc11+cfc12, adjusted
use_rad_ch4=1      ! use rad code CH4, shut off sfc sources
rad_FL=0           ! use rad code insolation getting fastj2 photon flux
which_trop=0       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
fix_CH4_chemistry=0    ! for setting fixed methane value for chemistry:
ch4_init_sh=0.791      ! init cond/fixed conditions SH CH4 ppmv
ch4_init_nh=0.791      ! init cond/fixed conditions NH CH4 ppmv
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)

! For altering tracer initial conditions and overwriting by a factor:
! set PI_run=1 and change the corresponding factors below. [For altering
! emissions, use the sectors above in the rundeck.
PI_run        = 1       ! =1 turns on below factors:
PIratio_N     = 0.667d0 ! {NOx, HNO3, N2O5, HO2NO2}
PIratio_CO_T  = 0.333d0 ! CO in troposphere
PIratio_CO_S  = 0.500d0 ! CO in stratosphere
PIratio_other = 0.500d0 ! {PAN,Isoprene,AlkyNit,Alkenes,Paraffin}
PIratio_N2O   = 1.000d0 ! {N2O ICs, L=1 overwrit}, set to 1 for use_rad_n2o=1
PIratio_CFC   = 1.000d0 ! {CFC ICs, L=1 overwrit}, set to 1 for use_rad_cfc=1
!--------- Aerosol parameters----------------
imAER=5         !3 historic; 1 AEROCOM ; 0,2 for standard or sector inputs (not working)
imPI=0          !for pre-industrial aerosols (natural-only) use imPI=1, imAER=5, aer_int_yr=1850
aer_int_yr=1850    !used for imAER=3, select desired year (1890 to 2000) or 0 to use JYEAR
rad_interact_aer=1  ! 1=couples aerosols to radiation, 0=use climatology

!--------- general aerosol parameters-------------
aer_rad_forc=0     ! 1: calculate aerosol radiative forcing
rad_forc_lev=1     ! 0: for TOA, 1: for tropopause for rad forcing diags.
prather_limits=1   ! 1: to avoid some negative tracers in sub-gridscale
diag_rad=1         ! 1: additional radiation diagnostics
diag_wetdep=1      ! 1: additional wet deposition diagnostics

!--------- sulfate and carbon aerosol parameters -----
madaer=3           ! 1: default sulfate and carbon aerosol 3: updated aerosols

!--------- dust aerosol parameters----------------
imDust=0           ! 0: PDF emission scheme, 1: AEROCOM
adiurn_dust=0      ! 1: daily dust diagnostics at certain grid points
!-------------------------------------------------
!----------------------------------------------------------------------
!  End tracer code parameters.
!----------------------------------------------------------------------

&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
 YEARE=1951,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 &END

