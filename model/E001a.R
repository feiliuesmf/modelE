E001a.R GISS Model E                                 gas 06/00

E001a: new modelE (based on B402A) prescribed ann.varying ocean data

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M12                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
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
AIC=DEC1958.rsfB394M12.modelE.16 ! init cond (atm. and ground)     ISTART=7
! AIC=AIC.RES_M12.D771201 ! initial conditions (atm. only)         ISTART=2
! GIC=GIC.E005gasA.1DEC1956 ! initial ground conditions            ISTART=2
! OSST/SICE: data rec. start with an integer M, where
!     M=1 means Jan of the year mentioned in the file names (here 1950)
!     In order that data and date agree, set IYEAR1=that year (1950)
OSST=OST4X5.B.1950.M02.Hadl1.1  ! ocean data   Feb 1950 - 1999
SICE=SICE4X5.B.1950.M02.Hadl1.1 ! sea ice data Feb 1950 - 1999
! OHT=OTSPEC.RunIDM12.M250D  ! hor.heat transp.  not needed if ocn prescribed
! OCNML=Z1O.B4X5.cor         ! mixed layer depth not needed if ocn prescribed
! EDDY=ED4X5 ! Eddy diffusivity for deep ocean mixing - not needed
CDN=CD4X500S VEG=V72X46.1.cor
SOIL=S4X50093 TOPO=Z72X46N.cor4 ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR      ! river direction file
RADN1=sgpgxg.table8    ! rad.tables
RADN2=kdist33.tautab8
RADN3=miescatpar.abcdv
RADN4=o3Prather1979-80.London1957-70
RADN5=TROAER.1875-1990.Jun2002
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean99.uvflux           ! need KSOLAR<2
! RADN9=solar.lean02.ann.uvflux     ! need KSOLAR=2
RADNA=O3.1850-2050.depl.rec       ! with recovery of O3 after 2000
!  RADNA=O3.1850-2050.depl.con    ! O3 'constant' after 2000
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Mar2002
dH2O=dH2O_by_CH4
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E001a (E001 - prescribed ann. varying ocean='ocean A')

DTFIX=300
&&PARAMETERS
X_SDRAG=.00025,.000025
C_SDRAG=0.
KOCEAN=0
ocn_cycl=0
U00ice=.55   ! use same values as in corr. run with climatological ocean
HRMAX=1000.  ! use same values as in corr. run with climatological ocean
KSOLAR=1
DT=450.,        ! from default: DTsrc=3600.,
dt_UVfilter=450.
SUBDD='SLP'     ! save SLP at sub-daily frequency
NSUBDD=12       ! saving sub-daily diags 12hrly
Kvflxo=1        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1951,MONTHI=1,DATEI=1,HOURI=0,
   YEARE=1957,MONTHE=1,DATEE=1,HOURE=0,    KDIAG=0,2,2,9*0,
   IYEAR1=1950  ! has to be consistent with OSST/SICE files !!!
   YEARE=1951,MONTHE=2,
   ISTART=7,IRANDI=0, YEARE=1951,MONTHE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
