C**** E001M12_COM
C**** COMMON BLOCK    4x5 deg Model  - Real*8 - VDATA:10->11+1

      MODULE E001M12_COM

C**** TO CHANGE THE GRID, MODIFY THE NEXT LINE ONLY
      INTEGER, PARAMETER :: IM=72,JM=46,LM=12, IMH=IM/2

C**** IM,JM,LM LIMITED TO 72,46,36 RESPECTIVELY BY RADCOM & SIGmas

C**** THERE ARE 100 INTEGER PARAMETERS IN COMMON (JC-ARRAY)
      INTEGER ::
     *                IM0,JM0,LM0,JMM1,LMM1,    LS1,LTM,LBLM,LMCM,LSSM,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRA,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0,MONTH
      INTEGER, DIMENSION(2) :: IDUM
      INTEGER, DIMENSION(13) :: NDZERO,NDPRNT
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC
      COMMON /IPARMB/ IM0,JM0,LM0,JMM1,LMM1,    LS1,LTM,LBLM,LMCM,LSSM,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRA,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0,MONTH,IDUM   ,NDZERO    ,NDPRNT    ,  IJD6     ,IDACC

! handle for referring to integer parameters
      INTEGER, DIMENSION(100) :: JC
      EQUIVALENCE (JC,IM0)

C**** THERE ARE 161 REAL NUMBERS IN COMMON (RC-ARRAY)
      DOUBLE PRECISION ::
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  CCMCX,U00,S0X,CO2,SRCOR,        PTOP,PSF,PSDRAG,PTRUNC,XINT,
     *  SKIPSE,USESLP,USEP,USET,FIM,    RSDIST,SIND,COSD
      DOUBLE PRECISION, DIMENSION(2) :: XCDNST
      DOUBLE PRECISION, DIMENSION(36) :: SIG
      DOUBLE PRECISION, DIMENSION(37) :: SIGE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(54) :: RDM2
      COMMON /RPARMB/
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  CCMCX,U00,S0X,CO2,SRCOR,        PTOP,PSF,PSDRAG,PTRUNC,
     *  XCDNST,XINT,                    SKIPSE,USESLP,USEP,USET,FIM,
     *  RSDIST,SIND,COSD,       SIG    ,SIGE    ,TAUTR0   ,RDM2

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

      DOUBLE PRECISION, DIMENSION(IM,JM) :: FLAND,FOCEAN,FLICE,FLAKE
     *     ,FEARTH,ZATMO
c      DOUBLE PRECISION, DIMENSION(IM,JM,5) :: ODATA
      DOUBLE PRECISION, DIMENSION(IM,JM,16) :: GDATA
      DOUBLE PRECISION, DIMENSION(IM,JM,11) :: VDATA
      DOUBLE PRECISION, DIMENSION(IM,JM) :: WFCS,Z1O,Z12O
      COMMON /BNDYCB/ FLAND,FOCEAN,FLICE,FLAKE,FEARTH,ZATMO,GDATA,
     *     VDATA,WFCS,Z1O,Z12O

      DOUBLE PRECISION, DIMENSION(36) :: DSIG,BYDSIG
      DOUBLE PRECISION, DIMENSION(35) :: DSIGO
c      COMMON /LAYACB/ DSIG,DSIGO

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q,WM
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P
c      COMMON /UVTPQ/ U,V,T,P,Q,WM

c      CONTAINS
c
c      SUBROUTINE io_model(iunit,irw)
c
c      INTEGER iunit  !@var iunit unit number of read/write
c      INTEGER irw    !@var irw   read or write flag
c      INTEGER, PARAMETER :: iwrite=1,iread=2
c      CHARACTER*20 HEADER,MODULE_HEADER
c      DATA MODULE_HEADER/"E001M12"
c
c      select case (irw)
c      case (iwrite)
c         write(iunit) TAU,MODULE_HEADER,JC,CLABEL,RC,KEYNR,U,V,T,P,Q
c      case (iread)
c         read(iunit) TAU,HEADER,JC,CLABEL,RC,KEYNR,U,V,T,P,Q
c         if (HEADER.ne.MODULE_HEADER)
c     &      print*,"Discrepancy in module version",HEADER,MODULE_HEADER
c      end select
c
c      RETURN
c      END SUBROUTINE io_model

      END MODULE E001M12_COM
