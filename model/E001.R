E001.R GISS Model E                                 gas 06/00

E001: new modelE (based on B402A)

Object modules: (in order of decreasing priority)
E001M12_COM SOMTQ_COM GEOM_B        ! model modules
GHYCOM const DAGCOM RADNCB          ! more model modules
PBLCOM                              ! more model modules
ME001M12 QUSEM12 DYNE001            ! daily dyn,Filt adv/avrx (2d ord mom.)
CE001M12 SUBSIDEM12 PE001M12        ! mstcnv,condse phys(no surfce)
SE001M12 EE001M12   PBLE001  SLE001 ! surfce and its subr
PBLDRV                              ! surfce and its subr
RE001                               ! setsur  rad.subr_incl._forcings 
DE001M12  FFT72  UTILDBL            ! diag,utilities

Data input files:
9=DEC1958.rsfB394M12
15=OST4X5.B.1946-55avg.Hadl1.1 17=SICE4X5.B.1946-55avg.Hadl1.1 ! ocn
19=CD4X500S 23=V72X46.1.cor 25=S4X50093 26=Z72X46N.cor4 ! bdy.cond
29=REG4X5           ! special regions-diag
71=sgpgxg.table8    ! rad.tables
72=kdist33.tautabs4 
73=miescatpar.abcdv 
74=o3Prather1979-80.London1957-70
75=trop8aer.tau5090 
76=dust8.tau9x8x13  
77=STRATAER.VOL.1950-2000.Jul99
78=cloud.epsilon4.72x46
79=solar.lean99.uvflux
80=o3trend.1951-2050  
81=o3WangJacob.1890.1979
84=topcld.trscat8

Label and Namelist:
E001 (new modelE based on B402A)
R=00BG/B
 &INPUTZ
   TAUI=0.,IYEAR=1950, CO2=-6.,
   KOCEAN=0,  CCMCX=1.,  S0X=1.,  IJRA=1, U00=.60,
   SIGE=1.0000000,.9400480,.8441247,.6834533,.4796163,.2877698,
        .1618705, .0719424, .0000000,-.0599520,-.1079137,-.1438849,
        -.1678657,
   TAUE=52560.,
   TAUE=744.,
   DT=450.,NDYN=8,NCNDS=8,NRAD=40,NFILTR=16,NDAA=58,NDA5D=56,
   NDA5S=56,NDA4=192,NDA5K=58,XINT=120.,
   PTOP=150.,  USET=24.,TAUO=0.,
   ISTART=3,KCOPY=2,NDPRNT=0,  TAUE=1.,USESLP=-12.,  USET=0.,USESLP=0.,
 &END
 &SDRNML
   XCDLM=.0005,.00005
 &END
