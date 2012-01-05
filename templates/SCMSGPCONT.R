SCMSGPCONT.R GISS Model E      awolf 08/2009   

SCM: RUN GISS Model E as SCM  using one latitude band                   
scm run using SGP Continuous Forcing Data from Jan 1999 - Dec 2001  
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)     ?
atmospheric composition from year 1979
ocean data: prescribed, 1975-1984 climatology
uses turbulence scheme (no dry conv), simple strat.drag (no grav.wave drag) ?
time steps: physics 30 min.; radiation 30 min.
filters: U,V in E-W direction (after every dynamics time step)              ?
         sea level pressure (after every physics time step)                 ?

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define SCM                          ! run as Single Column Model
#define NEW_IO
End Preprocessor Options

Object modules: (in order of decreasing priority)
#include "latlon_source_files"

ATM_COM
RES_F40                             ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
MODEL_COM                           ! model variables and geometry
IO_DRV                              ! new i/o
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
ALLOC_DRV                           ! domain decomposition, allocate 
                                    ! global distributed arrays
ATM_DRV                             ! driver for atmosphere-grid components
OCN_DRV                             ! driver for ocean-grid components
ATMDYN_COM ATMDYN_SCM MOMEN2ND  ! replace atmospheric dynamics with SCM routines
ATMDYN_SCM_EXT ATM_UTILS
SCM_COM SCMDATA_SGPCONT             ! routines for reading and processing SCM forcings and IC's
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE  FLUXES                     ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils and use ARM surface forcings     
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! ATURB - turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DUM                          ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO                    ! radiation and albedo
DIAG_COM DEFACC DIAG                ! diagnostics (diag, diag_prt dummies in scm_diag) 
DIAG_RES_F                          ! diagnostics (resolution dependent)
SCM_DIAG_COM SCM_DIAG               ! SCM diagnostics
      FFT144                        ! utilities

Components:
MPI_Support shared dd2d

Data input files:
AIC=AIC.RES_F40.D771201  ! observed init cond (atm. only) ISTART=2
GIC=GIC.144X90.DEC01.1.ext.nc   ! initial ground conditions      ISTART=2
OSST=OST_144x90.B.1975-1984avg.Hadl1 ! prescr. climatological ocean (1 yr data)
SICE=SICE_144x90.B.1975-1984avg.Hadl1 ! prescr. climatological sea ice
CDN=CD144X90 VEG=V144X90_no_crops CROPS=CROPS2007_144X90N_nocasp
SOIL=S144X900098M TOPO=Z144X90N_nocasp ! bdy.cond
REG=REG2X2.5          ! special regions-diag
RVR=RD_modelE_F.RVR.bin      ! river direction file

#include "rad_input_files"

#include "TAero2003_input_files"

#include "O3_2005_input_files"


GHG=GHG.Mar2004.txt
TOP_INDEX=top_index_144x90_a.ij.ext
MSU_wts=MSU.RSS.weights.data
GLMELT=GLMELT_144X90_gas.OCN   ! glacial melt distribution
SCMSRF=scm_sgpcont_0001_surface.dat
SCMLAY=scm_sgpcont_0001_layer.dat

Label and Namelist:
SCMSGPCONT (ModelE 2x2.5, 40 lyrs, 1979 atm/ocn;
up to 60 (or 52) columns here to describe your run)?<--col 53  to  72-->to 80-->
DTFIX=180.
&&PARAMETERS
! Target Coordinates for SCM
I_TARG=34        !Southern Great Plains 
J_TARG=64

! parameters set for prescribed ocean runs:
KOCEAN=0        ! ocn is prescribed
Kvflxo=0        ! use =1 to save VFLXO daily ONLY to prepare for q-flux runs
ocn_cycl=1      ! =0 if ocean varies from year to year

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

 
U00a=.74    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=2.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.57      ! tune this first to get: glob. ann. mean plan.alb=30%   (U00ice up=>albedo down)
U00wtrX=1.46    ! this to get: glob. ann. mean net heat at surf. = 0   (U00wtrX+.01=>NetHtSrf+.7)

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1979  ! if -1, crops in VEG-file is used
s0_yr=1979
s0_day=182
ghg_yr=1979
ghg_day=182
volc_yr=1979
volc_day=182
aero_yr=1979
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1979
dalbsnX=.024
o3_yr=-1979

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=225.
NIsurf=1        ! increase as layer 1 gets thinner
NRAD=1          ! Radiation called every dynamics time step
! parameters that affect at most diagn. output:
Ndisk=48 
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
   YEARI=2000,MONTHI=1,DATEI=3,HOURI=0, ! IYEAR1=YEARI (default) or earlier
   YEARE=2000,MONTHE=1,DATEE=31,HOURE=23,     KDIAG=12*0,9,
   ISTART=2,IRANDI=0, YEARE=2000,MONTHE=1,DATEE=3,HOURE=3
 &END
