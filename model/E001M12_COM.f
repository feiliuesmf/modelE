C**** E001M12_COM
C**** COMMON BLOCK    4x5 deg Model  - Real*8 - VDATA:10->11+1

      MODULE E001M12_COM

C**** TO CHANGE THE GRID, MODIFY THE NEXT LINE ONLY
C**** IM,JM,LM LIMITED TO 72,46,36 RESPECTIVELY BY RADCOM & SIGmas

!@var IM,JM longituidinal and latitudinal number of grid boxes
!@var LM number of vertical levels
      INTEGER, PARAMETER :: IM=72,JM=46,LM=12

!@var IMH half the number of latitudinal boxes
      INTEGER, PARAMETER :: IMH=IM/2

!@var FIM,BYIM real parameter values for the number of long. grid boxes
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM

!@var JEQ grid box immediately north of the equator
      INTEGER, PARAMETER :: JEQ=1+JM/2

C**** THERE ARE 100 INTEGER PARAMETERS IN COMMON (JC-ARRAY)
      INTEGER ::
     *       IM0,JM0,LM0,JMM1x,LMM1,    LS1,LTMx,LBLMx,LMCMx,LSSMx,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRAx,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  KTACC0,MONTH
      INTEGER, DIMENSION(2) :: IDUM
      INTEGER, DIMENSION(13) :: NDZERO,NDPRNT
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC
      COMMON /IPARMB/IM0,JM0,LM0,JMM1x,LMM1,LS1,LTMx,LBLMx,LMCMx,LSSMx,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRAx,MFILTR,NDYN,NCNDS,
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
     *  PTOP,PSF,PSDRAG,PTRUNC,XINT,
     *  SKIPSE,USESLP,USEP,USET,    RSDIST,SIND,COSD
      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(60+2*(36-LM)) :: RDM2
!@var PSFMPT,PSTRAT derived pressure constants
      REAL*8 :: PSFMPT,PSTRAT
      COMMON /RPARMB/
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  PTOP,PSF,PSDRAG,PTRUNC,
     *         XINT,                    SKIPSE,USESLP,USEP,USET,
     *  RSDIST,SIND,COSD,SIG,SIGE,      TAUTR0,PSFMPT,PSTRAT,RDM2

!@var RC handle for referring to real parameters
      DOUBLE PRECISION, DIMENSION(161) :: RC
      EQUIVALENCE (RC,TAU)

      CHARACTER*4 NAMD6,JMONTH,JMNTH0
      CHARACTER*132 XLABEL
      COMMON /TEXT/ XLABEL,NAMD6(4),JMONTH,JMNTH0

!@var LABLE1,CLABEL,XLABEL handles for referring to text parameters
      CHARACTER LABEL1*16
      CHARACTER CLABEL*156
      EQUIVALENCE (CLABEL,XLABEL,LABEL1)

      DOUBLE PRECISION, DIMENSION(IM,JM) :: FLAND,FOCEAN,FLICE,FLAKE
     *     ,FEARTH,ZATMO,HLAKE

      DOUBLE PRECISION, DIMENSION(IM,JM,16) :: GDATA
      DOUBLE PRECISION, DIMENSION(IM,JM,11) :: VDATA
      DOUBLE PRECISION, DIMENSION(IM,JM) :: WFCS
      COMMON /BNDYCB/ FLAND,FOCEAN,FLICE,FLAKE,FEARTH,ZATMO,GDATA,
     *     VDATA,WFCS

!@var DSIG,BYDSIG,DSIGO sigma level cooridinates
      DOUBLE PRECISION, DIMENSION(LM) :: DSIG,BYDSIG
      DOUBLE PRECISION, DIMENSION(LM-1) :: DSIGO

!@var Q_GISS,Q_HDF,Q_PRT,Q_NETCDF are switches for post-processing
      LOGICAL Q_GISS,Q_HDF,Q_PRT,Q_NETCDF
      COMMON /Q_PP/Q_GISS,Q_HDF,Q_PRT,Q_NETCDF

C**** Main model prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
!@var P surface pressure (hecto-Pascals - PTOP) 

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q,WM
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P

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
