E1M20.R GISS Model E  2004 modelE                 rar 12/01/03

E1M20: replace this section by a description of what distinguishes this run   ?
       Use as many lines as you need. Look carefully at all the possible      ?
       choices, particularly the lines containing '?'.
       The final rundeck should contain no '?'
       Check and modify the rest of the description below:                    ?
modelE1 (3.0) 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)         ?
atmospheric composition from year 1880 (or 1979)                   (see _yr)  ?
ocean data: prescribed, 1876-1885 (or 1975-1984) climatology  (see OSST/SICE) ?
uses turbulence scheme (not dry conv), simple strat.drag (not grav.wave drag) ?
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs    ?
filters: U,V in E-W direction (after every dynamics time step)                ?
         sea level pressure (after every physics time step)                   ?

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M20AT                           ! horiz/vert resolution, 4x5deg, 20 layers -> .1mb
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
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
    ! start up (from restart file of earlier run or) from observed conditions
! AIC=1DEC????.rsfE???   ! or:    ! initial conditions (atm./ground), no GIC, ISTART=8
AIC=AIC.RES_M20A.D771201          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955       ! initial conditions (ground)
    ! ocean data for "prescribed ocean" runs : climatological ocean
OSST=OST4X5.B.1876-85avg.Hadl1.1  ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1876-85avg.Hadl1.1 ! prescr. climatological sea ice
!1979 OSST=OST4X5.B.1975-84avg.Hadl1.1
!1979 SICE=SICE4X5.B.1975-84avg.Hadl1.1
    ! or:            annually varying ocean  (pick IYEAR1 appropriately, here 1871)
! OSST=OST4X5.B.1871.M02.Hadl1.1  ! ocean data   Feb 1871 - 2002
! SICE=SICE4X5.B.1871.M02.Hadl1.1 ! ocean data   Feb 1871 - 2002
    ! the next files are specific to q-flux ocean runs; replace files above
! AIC=E1M20/1JAN1960.rsfE1M20.MXL65m   ! AIC/OHT made by aux/mkOTSPEC
! OHT=E1M20/OTSPEC.E1M20.MXL65m.1951-1960 ! horizontal ocean heat transport
OCNML=Z1O.B4X5.cor                ! mixed layer depth (needed for post processing)
    ! files needed for all models
CDN=CD4X500S                      ! surf.drag coefficient
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops CROPS=CROPS_72X46N.cor4  ! veg. fractions, crops history
SOIL=S4X50093 TOPO=Z72X46N.cor4_nocasp   ! soil/topography bdy.conds
REG=REG4X5                        ! special regions-diag
RVR=RD4X525.gas2.RVR              ! river direction file
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=radfil33k                   !     8/2003 version
RADN3=miescatpar.abcdv2
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux      ! need KSOLAR=2
RADNE=topcld.trscat8
! ozone files (minimum 1, maximum 9 files + 1 trend file)
O3file_01=mar2004_o3_shindelltrop_72x46x49x12_1850
O3file_02=mar2004_o3_shindelltrop_72x46x49x12_1890
O3file_03=mar2004_o3_shindelltrop_72x46x49x12_1910
O3file_04=mar2004_o3_shindelltrop_72x46x49x12_1930
O3file_05=mar2004_o3_shindelltrop_72x46x49x12_1950
O3file_06=mar2004_o3_shindelltrop_72x46x49x12_1960
O3file_07=mar2004_o3_shindelltrop_72x46x49x12_1970
O3file_08=mar2004_o3_shindelltrop_72x46x49x12_1980
O3file_09=mar2004_o3_shindelltrop_72x46x49x12_1990
O3trend=mar2004_o3timetrend_46x49x2412_1850_2050
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
TOP_INDEX=top_index_72x46.ij
MSU_wts=MSU.RSS.weights.data

Label and Namelist:
E1M20 (ModelE1 4x5, 20 lyrs, 1880 atm/ocn; use up to 72 (or 80) columns and ??
up to 60 (or 52) columns here to describe your run)?<- col 53  to  72 ->   80 ->
DTFIX=300

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0 ! 0 or 1 , use =0 if ocn is prescribed, use =1 if ocn is predicted
Kvflxo=0 ! use 1 ONLY to save VFLXO daily to prepare for q-flux run ?
ocn_cycl=1  ! use =0 if ocean varies from year to year; irrelevant for pred. ocn

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

U00ice=.59      ! increase U00ice to decrease albedo    goals: NetHtz0=0,plan.alb=30%
U00wtrX=1.40    ! U00wtrX+.01=>nethtz0+.7                      for global annual mean
!use  U00wtrX=1.39    for 1979 atmosphere/ocean
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

CO2X=1.
H2OstratX=1.

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1880 ! if -1, crops in VEG-file is used   ! =1979 , also change OSST,SICE,U00wtrX
s0_yr=1880                                         ! =1979 , also change OSST,SICE,U00wtrX
s0_day=182
ghg_yr=1880                                        ! =1979 , also change OSST,SICE,U00wtrX
ghg_day=182
volc_yr=1880                                       ! =1979 , also change OSST,SICE,U00wtrX
volc_day=182
aero_yr=1880                                       ! =1979 , also change OSST,SICE,U00wtrX
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1880
dalbsnX=.015
o3_yr=1880
o3_yr=1880                                         ! =1979 , also change OSST,SICE,U00wtrX

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=450.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=48        ! use =480 on halem
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
cloud_rad_forc=0 ! use =1 to activate this diagnostic (doubles radiation calls !)
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default) or earlier
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 &END
!for q-flux run, also replace the above namelist by
! &INPUTZ
!   YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!   YEARE=1931,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
!   ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,
! &END
