E001B.R GISS Model E                                 gas 06/00

E001B: new modelE (based on B402A - coupled version)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution
MODEL_COM GEOM_B FLUXES             ! model modules
ME001M12                            ! Main and model overhead
MOMEN2ND DYNE001 DYNCOM             ! dynamics
SOMTQ_COM QUSDEF QUSEM12            ! advection of tracers
CLD01 CLD01_DRV_E001 CLD01_COM_E001 ! clouds modules
SE001M12                            ! surface calculation
GHYCOM EE001M12 SLE001              ! land surface and soils
PBLCOM PBLDRV PBLE001               ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ODIAG_COM OCEAN_COM OSTRAITS_COM OGEOM ! dynamic ocean modules
OCNDYN OSTRAITS OCNGM OCNKPP ICEDYN    ! dynamic ocean routines
OCNFUNTAB                           ! ocean function look up table
snowmodel                           ! snow model
RADNCB PE001M12 RE001               ! radiation modules
DAGCOM DE001M12 DEFACC DAGPRT       ! diagnostics
const FFT72 UTILDBL RAND~           ! utilities
POUT                                ! for post-processing
PARAM PARSER

Data input files:
AIC=DEC1958.rsfB394M12.modelE.12
OIC=OIC4X5LD.Z12.CLEV94.DEC01S  ! ocean initial conditions
OFTAB=OFTABLE_NEW               ! ocean function table
AVR=AVR4X5LD.Z12                ! ocean filter
TOPO_OC=Z72X46N.cor4 ! ocean bdy.cond
CDN=CD4X500S VEG=V72X46.1.cor
SOIL=S4X50093 TOPO=Z72X46N.cor4 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR      ! river direction file
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
E001B (new modelE based on B402A - coupled version)
R=00BG/B

&&PARAMETERS
CO2=-6.
XCDLM=.0005,.00005
KOCEAN=1
U00wtr=.50
U00ice=.50

DT=450.,        ! from default: DTsrc=3600.,
NSLP=0          ! saving SLP 0hrly
Kvflxo=0        ! not saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0,
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,
                       !  from default: IYEAR1=YEARI
   YEARE=1950,MONTHE=2,
   ISTART=7,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,
 &END
