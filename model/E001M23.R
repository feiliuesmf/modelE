E001M23.R GISS Model E                                 gas 06/00

E001M23: new modelE (based on B402A - strat. version)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M23                             ! horiz/vert resolution
MODEL_COM GEOM_B                    ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database 
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
CLOUDS CLOUDS_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
OCEAN                               ! ocean modules
SNOW                                ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL RAND~           ! utilities
POUT                                ! post-processing output

Data input files:
AIC=DEC1986.rsfB436TM23.modelE
OHT=OTSPEC.RB399AM12.M250D OCNML=Z1O.B4X5.cor
MLMAX=Z1OMAX.B4X5.250M.cor ! ocn data
OSST=OST4X5.B.1946-55avg.Hadl1.1 SICE=SICE4X5.B.1946-55avg.Hadl1.1 ! ocn
CDN=CD4X500S VEG=V72X46.1.cor
SOIL=S4X50093 TOPO=Z72X46N.cor4 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR      ! river direction file
ZVAR=ZVAR4X5         ! topogrphic variation for gwdrag
RADN1=sgpgxg.table8    ! rad.tables
RADN2=kdist33.tautabs4
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=trop8aer.tau5090
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1950-2000.Jul99
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean99.uvflux
RADNA=o3trend.1951-2050
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E001M23 (new modelE based on B402A - strat. version)
R=00BG/B

&&PARAMETERS
CO2=-6.
XCDLM=.0005,.00005
KOCEAN=0
U00wtr=.50
U00ice=.50

LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DT=180.,             ! from default: DTsrc=3600.,
NIsurf=4,            ! number of surface time steps

NSLP=0          ! saving SLP 12hrly
Kvflxo=0        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0,
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,
                       !  from default: IYEAR1=YEARI
   YEARE=1950,MONTHE=2,
   ISTART=7,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,
 &END
