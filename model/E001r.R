E001r.R GISS Model E  strat.H2O added inst.frc          rar 2/01/02

E001r: radiation only, instantaneous forcing

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
AIC=1JAN1951.rsfE001 ! just label records of control run are needed
! GIC=GIC.E005gasA.1DEC1956 ! initial conditions (ground)
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
RADN9=solar.lean99.uvflux             ! need KSOLAR<2
! RADN9=solar.lean02.ann.uvflux       ! need KSOLAR=2
RADNA=o3trend.1850-2050 
RADNB=o3WangJacob.1890.1979
RADNE=topcld.trscat8
GHG=GHG.1850-2050.Sep2001
dH2O=dH2O_by_CH4
TOP_INDEX=top_index_72x46.ij
RADJAN=../E001/RADJAN1950    ! replace E001 by the control run
RADFEB=../E001/RADFEB1950    ! assuming that the saved data are
RADMAR=../E001/RADMAR1950    ! still in the original directory
RADAPR=../E001/RADAPR1950    ! and that the forcing run is on
RADMAY=../E001/RADMAY1950    ! the same file system
RADJUN=../E001/RADJUN1950
RADJUL=../E001/RADJUL1950
RADAUG=../E001/RADAUG1950
RADSEP=../E001/RADSEP1950
RADOCT=../E001/RADOCT1950
RADNOV=../E001/RADNOV1950
RADDEC=../E001/RADDEC1950

Label and Namelist:
E001r (inst.forcing run - control)
R=00BG/B

&&PARAMETERS
KSOLAR=1 ! all params affecting radia are relevant
         ! most other parameters are irrelevant
NSUBDD=0            ! don't touch this line
Kvflxo=0            ! don't touch this line
KCOPY=2             ! saving acc + rsf
Kradia=1            ! use Kradia=2 for adj. forcing run
&&END_PARAMETERS

 &INPUTZ
   YEARE=1952,MONTHE=1,DATEE=1,HOURE=0,  ! assumed start: 1/1/1951
   ISTART=8, YEARE=1951,MONTHE=1,HOURE=1,IWRITE=1,JWRITE=1,ITWRITE=23,
 &END

It's important that DTsrc, NIrad, ItimeI are left as in the control
run so the input data are available when they are needed ! That's
why in this case ISTART=8 (start) looks more like ISTART=9 (restart).
The only parameters/input files that matter are the ones that affect
the radiation. It is their change whose instantaneous/adjusted
forcing is computed.
