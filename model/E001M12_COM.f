C**** E001M12_COM
C**** COMMON BLOCK    4x5 deg Model  - Real*8 - VDATA:10->11+1

      MODULE E001M12_COM

C**** TO CHANGE THE GRID, MODIFY THE NEXT LINE ONLY
      INTEGER, PARAMETER :: IM=72,JM=46,LM=12, KTD=8,KAIJ=100,KAJK=51
      integer, parameter :: lmom = 9 ! for deep ocean diffusion  
C**** IM,JM,LM LIMITED TO 72,46,36 RESPECTIVELY BY RADCOM & SIGmas

C**** THERE ARE 100 INTEGER PARAMETERS IN COMMON (JC-ARRAY)
      INTEGER ::
     *                IM0,JM0,LM0,JMM1,LMM1,    LS1,LTM,LBLM,LMCM,LSSM,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRA,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0
      INTEGER, DIMENSION(3) :: IDUM
      INTEGER, DIMENSION(13) :: NDZERO,NDPRNT
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC
      COMMON /IPARMB/ IM0,JM0,LM0,JMM1,LMM1,    LS1,LTM,LBLM,LMCM,LSSM,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRA,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0,IDUM   ,NDZERO    ,NDPRNT    ,  IJD6     ,IDACC    

! handle for referring to integer parameters
      INTEGER, DIMENSION(100) :: JC
      EQUIVALENCE (JC,IM0)

C**** THERE ARE 161 REAL NUMBERS IN COMMON (RC-ARRAY)
      DOUBLE PRECISION ::
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  TWOPI,SDAY,LHE,LHM,LHS,         RADIUS,GRAV,RGAS,KAPA,OMEGA,
     *  CCMCX,U00,S0X,CO2,SRCOR,        PTOP,PSF,PSDRAG,PTRUNC,AREAG,
     *            XINT,DLAT,DLON,       SKIPSE,USESLP,USEP,USET,FIM,
     *  RSDIST,SIND,COSD,DOPK
      DOUBLE PRECISION, DIMENSION(2) :: XCDNST
      DOUBLE PRECISION, DIMENSION(36) :: SIG
      DOUBLE PRECISION, DIMENSION(37) :: SIGE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(40) :: RDM2
      COMMON /RPARMB/
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  TWOPI,SDAY,LHE,LHM,LHS,         RADIUS,GRAV,RGAS,KAPA,OMEGA,
     *  CCMCX,U00,S0X,CO2,SRCOR,        PTOP,PSF,PSDRAG,PTRUNC,AREAG,
     *  XCDNST   ,XINT,DLAT,DLON,       SKIPSE,USESLP,USEP,USET,FIM,
     *  RSDIST,SIND,COSD,DOPK,     SIG    ,SIGE    ,TAUTR0   ,RDM2    

! handle for referring to real parameters
      DOUBLE PRECISION, DIMENSION(161) :: RC
      EQUIVALENCE (RC,TAU)

      CHARACTER*4 XLABEL,NAMD6,JMONTH,JMNTH0
      COMMON /TEXT/ XLABEL(33),NAMD6(4),JMONTH,JMNTH0

! handles for referring to text parameters
      CHARACTER*8 LABEL1*16
      CHARACTER CLABEL*4
      DIMENSION CLABEL(39)
      EQUIVALENCE (CLABEL,XLABEL,LABEL1)

      DOUBLE PRECISION, DIMENSION(JM) ::
     *  RAPVS,RAPVN,RAVPS,RAVPN,FCOR,
     *  DXYP,DXP,DYP,DXYS,SINP,LAT,
     *  DXYV,DXV,DYV,DXYN,COSP,COSV
      COMMON /GEOMCB/
     *  RAPVS,RAPVN,RAVPS,RAVPN,FCOR,
     *  DXYP,DXP,DYP,DXYS,SINP,LAT,
     *  DXYV,DXV,DYV,DXYN,COSP,COSV

      DOUBLE PRECISION, DIMENSION(IM,JM,3) :: FDATA
      DOUBLE PRECISION, DIMENSION(IM,JM,5) :: ODATA
      DOUBLE PRECISION, DIMENSION(IM,JM,16) :: GDATA
      DOUBLE PRECISION, DIMENSION(IM,JM,12) :: BLDATA
      DOUBLE PRECISION, DIMENSION(IM,JM,11) :: VDATA
      DOUBLE PRECISION, DIMENSION(IM,JM) :: WFCS,Z1O,Z12O
      COMMON /BNDYCB/ FDATA,ODATA,GDATA,BLDATA,VDATA,WFCS,Z1O,Z12O

      DOUBLE PRECISION, DIMENSION(IM,JM,3) :: RQT
      DOUBLE PRECISION, DIMENSION(IM,JM,LM+1) :: SRHR,TRHR
      COMMON /RADNCB/ RQT,SRHR,TRHR

      DOUBLE PRECISION, DIMENSION(36) :: DSIG
      DOUBLE PRECISION, DIMENSION(35) :: DSIGO
      COMMON /LAYACB/ DSIG,DSIGO

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P
      COMMON /UVTPQ/ U,V,T,P,Q

C**** DIAGNOSTIC ARRAYS
      INTEGER, PARAMETER :: IMH=IM/2,     KACC=JM*94*3 + 24*94 + JM*3 +
     *   JM*LM*57 + JM*3*4 + IM*JM*KAIJ + IM*LM*16 +
     *   IM*JM*29 + 20*100 + JM*36 + (IMH+1)*20*8 +
     *   8*2 + 24*63*4 + 2*62*10*12 + JM*LM*KAJK +
     *   IM*JM*LM*6 + IM*JM*LM*5 + JM*LM*10*3

      DOUBLE PRECISION, DIMENSION(JM,94) :: AJ,BJ,CJ
      DOUBLE PRECISION, DIMENSION(24,94) :: DJ
      DOUBLE PRECISION, DIMENSION(JM,3) :: APJ
      DOUBLE PRECISION, DIMENSION(JM,LM,57) :: AJL
      DOUBLE PRECISION, DIMENSION(JM,3,4) :: ASJL
      DOUBLE PRECISION, DIMENSION(IM,JM,KAIJ) :: AIJ
      DOUBLE PRECISION, DIMENSION(IM,LM,16) :: AIL
      DOUBLE PRECISION, DIMENSION(IM,JM,29) :: AIJG
      DOUBLE PRECISION, DIMENSION(20,100) :: ENERGY
      DOUBLE PRECISION, DIMENSION(JM,36) :: CONSRV
      DOUBLE PRECISION, DIMENSION((IMH+1),20,8) :: SPECA
      DOUBLE PRECISION, DIMENSION(8,2) :: ATPE
      DOUBLE PRECISION, DIMENSION(24,63,4) :: ADAILY
      DOUBLE PRECISION, DIMENSION(2,62,10,12) :: WAVE
      DOUBLE PRECISION, DIMENSION(JM,LM,KAJK) :: AJK
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,6) :: AIJK
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,5) :: AIJL
      DOUBLE PRECISION, DIMENSION(JM,LM,10,3) :: AJLSP
      DOUBLE PRECISION, DIMENSION(IM,JM,2) :: TSFREZ
      DOUBLE PRECISION, DIMENSION(IM,JM,KTD) :: TDIURN
      COMMON /ACCUM/ AJ,BJ,CJ,DJ,APJ,AJL,ASJL,AIJ,AIL,
     *  AIJG,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,
     *  AJK,AIJK,AIJL,AJLSP,TSFREZ,TDIURN

      INTEGER, DIMENSION(IM,JM) :: JREG
      COMMON /REGION/ JREG
      INTEGER, DIMENSION(42,50) :: KEYNR
      INTEGER, DIMENSION(12) :: KDIAG
      COMMON /KEYS/ KEYNR,KDIAG

C**** Some helpful arrays 

!@var  AM  Air mass of each box (kg)
!      REAL*8, DIMENSION(IM,JM,LM) :: AM
!@var  BYAM  1/Air mass (kg**-1)
!      REAL*8, DIMENSION(IM,JM,LM) :: BYAM
!@var  PK   P**KAPA 
      REAL*8, SAVE,DIMENSION(IM,JM,LM) :: PK
!@var  PEUP  Pressure at upper edge of box (mb)
!      REAL*8, DIMENSION(IM,JM,LM) :: PEUP
!@var  PMID  Pressure at mid point of box (mb)
!      REAL*8, DIMENSION(IM,JM,LM) :: PMID

      END MODULE E001M12_COM
