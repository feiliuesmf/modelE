E001M20A.R GISS Model E  2002 modelE                 rar 12/01/03

modelE 2.4 with 20 lyrs, top at .1 mb - 1880 atmosphere/ocean
!M12 modelE 2.4 with 12 lyrs, top at 10 mb - 1880 atmosphere/ocean
no gravity wave drag;     uses turbulence (not dry convection)
Sdrag: weak linear strat. drag in top layer, near poles down to 20 mb
!M12 Sdrag: weak linear strat. drag in top layer only
       ang.mom loss is added in below 150 mb
sea level pressure filter applied every hour, UV-filter used
6-band oice albedo; Hogstrom(1984) pbl drag
Note: Some of these choices may be changed using the PARAMETERs below.

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M20AT                           ! horiz/vert resolution, 4x5deg, 20 layers -> .1mb
!M12 RES_M12                        ! horiz/vert resolution, 4x5deg, 20 layers -> .1mb
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
DOMAIN_DECOMP ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
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
ICEDYN_DRV ICEDYN  ! or: ICEDYN_DUM ! dynamic sea ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
    ! the first group of files is specific to prescribed ocean runs
! AIC=1DEC1951.rsfE000   ! or:    ! initial conditions (atm./ground), no GIC, ISTART=8
AIC=AIC.RES_M20A.D771201          ! initial conditions (atm.)      needs GIC, ISTART=2
!M12 AIC=AIC.RES_M12.D771201      ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955       ! initial conditions (ground)
OSST=OST4X5.B.1876-85avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1876-85avg.Hadl1.1 ! prescr. climatological sea ice
!1979 OSST=OST4X5.B.1975-84avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
!1979 SICE=SICE4X5.B.1975-84avg.Hadl1.1 ! prescr. climatological sea ice
    ! if the prescr. ocean varies from year to year use instead:
! OSST=OST4X5.B.1871.M02.Hadl1.1  ! ocean data   Feb 1871 - 2002
! SICE=SICE4X5.B.1871.M02.Hadl1.1 ! ocean data   Feb 1871 - 2002
    ! the next group files are specific to q-flux ocean runs
! AIC=E001M20A/1JAN1960.rsfE001M20A.MXL65m   ! AIC/OHT made by aux/mkOTSPEC
! OHT=E001M20A/OTSPEC.E001M20A.MXL65m.1951-1960 ! horizontal ocean heat transport
OCNML=Z1O.B4X5.cor                ! mixed layer depth (use for post processing)
    ! files needed for all models
CDN=CD4X500S                      ! surf.drag coefficient
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops CROPS=CROPS_72X46N.cor4  ! veg. fractions, crops history
SOIL=S4X50093 TOPO=Z72X46N.cor4   ! soil/topography bdy.conds
REG=REG4X5                        ! special regions-diag
RVR=RD4X525.RVR                   ! river direction file
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=radfil33k                   !     8/2003 version
RADN3=miescatpar.abcdv
! RADN4,RADN5,RADNA,RADNB are no longer used
TAero_PRE=sep2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux      ! need KSOLAR=2
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
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E001M20A (ModelE 1880 atm/ocn)

DTFIX=180

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0        ! ocn is prescribed
Kvflxo=0        ! do NOT save VFLXO (daily) (use 1 to prepare for q-flux run)
ocn_cycl=1      ! =0 if ocean varies from year to year

! parameters usually not changed when switching to q-flux ocean:

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.001,.0001  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0001       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=20.        ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)

! drag params if top is at 10mb
!M12 X_SDRAG=.00025,.000025  ! very weak linear drag
!M12 C_SDRAG=0.      ! no constant drag
!M12 P_sdrag=0.      ! linear SDRAG only in top layer
!M12 PP_sdrag=0.     ! linear SDRAG only in top layer even at poles

ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

U00ice=.60      ! U00ice up => nethtz0 down (alb down); goals: nethtz0=0,plan.alb=30%
U00wtrX=1.25    ! U00wtrX+.01=>nethtz0+.5   (alb down);        for global annual mean
!        U00wtrX=1.25    ! use with 1880 atmosphere/ocean
!1979    U00wtrX=1.22    ! use with 1979 atmosphere/ocean
!M12     U00wtrX=1.18    ! use with 1880 atmosphere/ocean AND RES_M12
!M121979 U00wtrX=1.15    ! use with 1979 atmosphere/ocean AND RES_M12
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

RWCLDOX=1.5  !  wtr cld particle size *3/2 over ocean
RICLDX=.3333 !  ice cld particle size * 1(at 0mb)->1/3(at 1000mb)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1880 ! if -1, crops in VEG-file is used
s0_yr=1880
s0_day=182
ghg_yr=1880
ghg_day=182
volc_yr=1880
volc_day=182
aero_yr=1880
o3_yr=1880

! parameters that control the Shapiro filter
DT_XUfilter=300. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=300. ! Shapiro filter on V in E-W direction; usually same as DT (below)
!M12 DT_XUfilter=450.
!M12 DT_XVfilter=450.
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DT=300.         ! from default: DTsrc=3600.,
!M12 DT=450.    ! from default: DTsrc=3600., (also use DTFIX=300 , not 180)
NIsurf=2        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=24        ! use =240 on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags 0hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=1         ! use =7 to save cpu time, but energy cons. diag may not be accurate
nda5s=1         ! use =7 to save cpu time, but energy cons. diag may not be accurate
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default)
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=0,2,2,9*0,9,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 &END
