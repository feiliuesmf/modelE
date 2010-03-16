E4TdusF40.R GISS Model E  2000 ocn/atm   jan perlwitz  07/2009
 Template for setting up simulations using the E4F40 model with dust tracers

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
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define USE_ENT                  ! include dynamic vegetation model
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer 
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define NO_HDIURN                ! exclude hdiurn diagnostics
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
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
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
OSST=OST_144x90.1996-2005avg.HadISST1.1    ! prescr. climatological ocean
SICE=SICE_144x90.1996-2005avg.HadISST1.1   ! prescr. climatological sea ice
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
RADN9=solar.lean02.ann.uvflux_hdr       ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU.RSS.weights.data

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

! files for dust tracers
ERS=ERS1_1993_MONTHLY.144x90.threshold-13 ! ERS data
DSRC=Ginoux_source_v2009_VegMask_144x90   ! preferred dust sources
                                          ! optimized use: prefDustSources=1
! alternative preferred dust source files:
!  Ginoux2001_source_VegMask_144x90     (optimized use: prefDustSources=0)
!  Ginoux_source_v2009_NoVegMask_144x90 (optimized use: prefDustSources=2)
!  GriniZender_DustSources_144x90       (optimized use: prefDustSources=3)
!  Tegen_DustSources_144x90             (optimized use: prefDustSources=4)
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

Label and Namelist:  (next 2 lines)
E4TdusF40 (modelE4 2x2.5 hor., 40 lyrs, 2000 atm., 1996-2005 clim ocn, dust tracers)
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

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2
cloud_rad_forc=1 ! 1: calculate cloud radiative forcing

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=2000  ! if -1, crops in VEG-file is used
s0_yr=2000
s0_day=182
ghg_yr=2000
ghg_day=182
volc_yr=-1
volc_day=182
aero_yr=2000
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
albsn_yr=2000
dalbsnX=.024
o3_yr=2000
CO2X=1.
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

calc_orb_par=1
paleo_orb_yr=-50.  !  BP i.e. 1950-paleo_orb_yr AD = 1850 AD

!--------- general aerosol parameters-------------
aer_rad_forc=0     ! 1: calculate aerosol radiative forcing
rad_forc_lev=1     ! 0: for TOA, 1: for tropopause for rad forcing diags.
rad_interact_aer=1 ! 1: couples aerosols to radiation, 0: use climatology
prather_limits=1   ! 1: to avoid some negative tracers in sub-gridscale
diag_rad=1         ! 1: additional radiation diagnostics
diag_wetdep=1      ! 1: additional wet deposition diagnostics
to_conc=0,1,1,1,1,1 ! 1: taijln diags as concentration; 0: as mixing ratio

!--------- sulfate and carbon aerosol parameters -----
madaer=3           ! 1: default sulfate and carbon aerosol 3: updated aerosols

!--------- dust aerosol parameters----------------
imDust=0           ! 0: PDF emission scheme, 1: AEROCOM
adiurn_dust=0      ! 1: daily dust diagnostics at certain grid points
prefDustSources=1     ! 0: Ginoux 2001 w/ vegetation mask
                      ! 1: Ginoux 2009 w/ vegetation mask (current default)
                      ! 2: Ginoux 2009 w/o vegetation, 3: Grini/Zender sources
                      ! 4: Tegen sources, >4: Free choice of emis. parameters
!fracClayPDFscheme=1. ! Frac. clay emis, only effective for prefDustSources > 4
!fracSiltPDFscheme=1. ! Frac. silt emis, only effective for prefDustSources > 4
                      ! set internally for 0-4

!-------------------------------------------------

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
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
Nssw=2   ! until diurnal diags are fixed, Nssw has to be even

&&END_PARAMETERS

 &INPUTZ
 YEARI=1999,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
 YEARE=2005,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1999,MONTHE=12,DATEE=1,HOURE=1,
 &END
