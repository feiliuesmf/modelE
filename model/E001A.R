E001A.R GISS Model E                                 gas 06/00

E001A: new modelE (Qflux version)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution
MODEL_COM GEOM_B                    ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
CLOUDS CLOUDS_DRV CLOUDS_COM        ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
! ATURB                             ! turbulence in whole atmosphere
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
!
! The next 2 lines depend on the current control with prescr. ocean
AIC=IMJAN.E001.O250D               ! initial conditions
OHT=OTSPEC.E001.M250D.1951-1955    ! horizontal ocean heat transports
! base the above 2 lines on the preliminary run with prescribed ocean
!
OCNML=Z1O.B4X5.cor                 ! ocn data: mixed layer depth
MLMAX=Z1OMAX.B4X5.250M.cor         ! ocn data: ann. max. mixed layer depth
! OSST=OST4X5.B.1946-55avg.Hadl1.1 ! not needed since ocn is predicted
! SICE=SICE4X5.B.1946-55avg.Hadl1.1 !not needed since ocn is predicted
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
dH2O=dH2O_by_CH4
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E001A (new modelE based on B402A - Qflux version)
R=00BG/B

&&PARAMETERS
XCDLM=.0005,.00005
KOCEAN=1
U00wtr=.49
U00ice=.50

DT=450.         ! from default: DTsrc=3600.,
NSLP=12         ! saving SLP 12hrly
Kvflxo=0        ! don't save VFLXO (daily) if ocn is predicted
KCOPY=2         ! saving acc + rsf
isccp_diags=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default)
   YEARE=1980,MONTHE=1,DATEE=1,HOURE=0,
   ISTART=7,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,
 &END
