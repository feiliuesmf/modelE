E001.R GISS Model E                                 gas 06/00

E001: new modelE (based on B402A)

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
AIC=DEC1958.rsfB394M12.modelE.15 ! initial conditions (atm. and ground)
! GIC=GIC.E005gasA.1DEC1956 ! initial ground conditions (needed if ISTART=2)
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
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean99.uvflux         ! need KSOLAR<2
!  RADN9=solar.lean02.ann.uvflux  ! need KSOLAR=2
RADNA=o3trend.1850-2050
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Sep2001
dH2O=dH2O_by_CH4
TOP_INDEX=top_index_72x46.ij

Label and Namelist:
E001 (new modelE based on B402A, uses dry adiab. adjustment)
R=00BG/B

&&PARAMETERS
XCDLM=.0005,.00005
KOCEAN=0
U00ice=.50   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
HRMAX=1000.  ! tune this to get rad.equilibrium (range: 100.-1500. meters)
KSOLAR=1
DT=450.,        ! from default: DTsrc=3600.,
SUBDD='SLP'     ! save SLP at sub-daily frequency
NSUBDD=12       ! saving sub-daily diags 12hrly
Kvflxo=1        ! saving VFLXO (daily)
KCOPY=2         ! saving acc + rsf
isccp_diags=1
&&END_PARAMETERS

 &INPUTZ
   YEARI=1950,MONTHI=1,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default)
   YEARE=1956,MONTHE=1,DATEE=1,HOURE=0,
   YEARE=1950,MONTHE=2,
   ISTART=7,IRANDI=0, YEARE=1950,MONTHE=1,HOURE=1,IWRITE=1,JWRITE=1,
 &END
