E001A.R GISS Model E                                 gas 06/00

E001A: new modelE (Qflux version)

Object modules: (in order of decreasing priority)
E001M12_COM SOMTQ_COM GEOM_B        ! model modules
GHYCOM const DAGCOM RADNCB          ! more model modules
PBLCOM LAKES_COM                    ! more model modules
ME001M12 QUSEM12 DYNE001            ! daily dyn,Filt adv/avrx (2d ord mom.)
CLD01 CLD01_DRV_E001 CLD01_COM_E001 ! clouds modules
SUBSIDEM12 PE001M12                 ! phys(no surfce)
SE001M12 EE001M12   PBLE001  SLE001 ! surfce and its subr
PBLDRV                              ! surfce and its subr
OCNE001 SEAICE LAKES LANDICE        ! ocean, lake and sea ice modules
RE001                               ! setsur  rad.subr_incl._forcings 
DE001M12 FFT72 UTILDBL RANDSGI      ! diag,utilities
POUT                                ! for post-processing
snowmodel sweep3d                   ! snow model

Data input files:
AIC=DEC1958.rsfB394M12
OHT=OTSPEC.RB399AM12.M250D OCNML=Z1O.B4X5.cor 
MLMAX=Z1OMAX.B4X5.250M.cor ! ocn data
OSST=OST4X5.B.1946-55avg.Hadl1.1 SICE=SICE4X5.B.1946-55avg.Hadl1.1 ! ocn
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

Label and Namelist:
E001A (new modelE based on B402A)
R=00BG/B
 &INPUTZ
   TAUI=0.,IYEAR=1950, CO2=-6.,
   KOCEAN=1,  S0X=1., U00wtr=.60, U00ice=.60,
   SIGE=1.0000000,.9400480,.8441247,.6834533,.4796163,.2877698,
        .1618705, .0719424, .0000000,-.0599520,-.1079137,-.1438849,
        -.1678657,
   TAUE=52560.,
   TAUE=744.,
   DT=450.,NDYN=8,NCNDS=8,NRAD=40,NFILTR=16,NDAA=58,NDA5D=56,
   NDA5S=56,NDA4=192,NDA5K=58,XINT=120.,
   PTOP=150.,  USET=24., TAUO=0.,
   XCDLM=.0005,.00005,
   ISTART=3,KCOPY=2,NDPRNT(1)=-1,  TAUE=1.,USESLP=-12.,  USET=0.,USESLP=0.,
 &END
