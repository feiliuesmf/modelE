E001tr.R GISS Model E with sample tracers               jal 10/01

E001tr: new modelE (based on B402A with 4 sample tracers)
 T1=Air mass
 T2=SF6
 T3=RN222
 T4=CO2

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
!#define TRACERS_WATER               ! include water tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution
MODEL_COM GEOM_B                    ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q and tracer gases
! pick one of the next two lines
TRACER_COM TRACERS_Air              ! air mass tracers
! TRACER_Water_COM TRACERS_Water    ! water mass tracers
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
CLOUDS CLOUDS_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices ATURB or DRYCNV
! ATURB                               ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
OCEAN OCNML                         ! ocean modules
SNOW                                ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT                                ! post-processing output

Data input files:
AIC=DEC1958.rsfB394M12.modelE.13 ! initial conditions (atm. and ground)
! OHT=OTSPEC.RunIDM12.M250D  ! hor.heat transp.  not needed if ocn prescribed
OCNML=Z1O.B4X5.cor         ! mixed layer depth,needed for post-processing only
MLMAX=Z1OMAX.B4X5.250M.cor ! ann max mix.l.dp.,needed for post-processing only
OSST=OST4X5.B.1946-55avg.Hadl1.1 ! prescr. climatological ocean (1 yr of data)
SICE=SICE4X5.B.1946-55avg.Hadl1.1 ! prescr. climatological sea ice
CDN=CD4X500S VEG=V72X46.1.cor
SOIL=S4X50093 TOPO=Z72X46N.cor4 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR      ! river direction file
RADN1=sgpgxg.table8    ! rad.tables
RADN2=kdist33.tautabs4
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=trop8aer.tau5090.minimum
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1950-2000.Jul99
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean99.uvflux
RADNA=o3trend.1951-2050.2
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Oct2000
TOP_INDEX=top_index_72x46.ij
CO2_IC=CO2ijl_IC_Jan1_scale334_M23  !wofsy+B140TQaM9
CO2_FOS_FUEL=CO2_sources/gcm_data/CO2FOS_MRL_4X5
CO2_FERT=CO2_sources/gcm_data/CO2fert01_4X5
CO2_REGROWTH=CO2_sources/gcm_data/CO2_Nforest_4X5
CO2_LAND_USE=CO2_sources/gcm_data/CO2DEF_HOU_4X5
CO2_VEG=CO2_sources/gcm_data/CO2VEG_MON_4X5        ! Monthly source
CO2_OCEAN=CO2_sources/gcm_data/CO2_4X5_Ocean_flux02  ! Monthly source  x 1.749

Label and Namelist:
E001tr (new modelE based on B402A, uses dry adiab. adjustment; tracers)
R=00BG/B

&&PARAMETERS
CO2=-6.
XCDLM=.0005,.00005
KOCEAN=0
U00wtr=.49
U00ice=.50

DT=450.,        ! from default: DTsrc=3600.,
to_volume_MixRat=1   ! for tracer printout
NSLP=12         ! saving SLP 12hrly
Kvflxo=1        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default)
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,
   YEARE=1950,MONTHE=2,
   ISTART=7,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,
 &END
