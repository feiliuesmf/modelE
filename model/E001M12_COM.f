      MODULE E001M12_COM
!@sum  E001M12_COM Main model variables 4x5 deg Model, 12 layers
!@auth Original Development Team
!@ver  1.0

!@var IM,JM longitudinal and latitudinal number of grid boxes
!@var LM number of vertical levels (limited to 40 by BR00B.COM)
      INTEGER, PARAMETER :: IM=72,JM=46,LM=12

!@var IMH half the number of latitudinal boxes
      INTEGER, PARAMETER :: IMH=IM/2
!@var FIM,BYIM real parameter values for the number of long. grid boxes
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM
!@var JEQ grid box immediately north of the equator
      INTEGER, PARAMETER :: JEQ=1+JM/2

C**** THERE ARE 100 INTEGER PARAMETERS IN COMMON (JC-ARRAY)
      INTEGER ::
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Kvflxo,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH
      INTEGER, DIMENSION(25) :: IDUM
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(12) :: IDACC
      COMMON /IPARMB/
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Kvflxo,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH,NIPRNT,NMONAV,  IDUM    ,  IJD6     ,IDACC

! handle for referring to integer parameters
      INTEGER, DIMENSION(100) :: JC
      EQUIVALENCE (JC,IM0)

C**** THERE ARE 161 REAL NUMBERS IN COMMON (RC-ARRAY)
      DOUBLE PRECISION ::
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG, SKIPSE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(161-13-2*LM) :: RDM2
!@var PSFMPT,PSTRAT derived pressure constants
      COMMON /RPARMB/
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG,  SKIPSE,
     *  TAUTR0,SIG,SIGE,  RDM2

!@var RC handle for referring to real parameters
      DOUBLE PRECISION, DIMENSION(161) :: RC
      EQUIVALENCE (RC,DTsrc)

      CHARACTER*4 NAMD6,AMON,AMON0
      CHARACTER*132 XLABEL
      COMMON /TEXT/ XLABEL,NAMD6(4),AMON,AMON0

!@var LABEL1,CLABEL,XLABEL handles for referring to text parameters
      CHARACTER LABEL1*16
      CHARACTER CLABEL*156
      EQUIVALENCE (CLABEL,XLABEL,LABEL1)
      DOUBLE PRECISION, DIMENSION(IM,JM) :: FLAND,FOCEAN,FLICE,FLAKE
     *     ,FEARTH,ZATMO,HLAKE

      DOUBLE PRECISION, DIMENSION(IM,JM,11) :: VDATA
      DOUBLE PRECISION, DIMENSION(IM,JM) :: WFCS
c      COMMON /BNDYCB/ FLAND,FOCEAN,FLICE,FLAKE,FEARTH,ZATMO,VDATA,WFCS

!@var DSIG,BYDSIG,DSIGO sigma level cooridinates
      DOUBLE PRECISION, DIMENSION(LM) :: DSIG,BYDSIG
      DOUBLE PRECISION, DIMENSION(LM-1) :: DSIGO

!@var XCDLM.  SDRAG ~XCDLM(1)+XCDLM(2)*wind_magnitude
      DOUBLE PRECISION, DIMENSION(2) :: XCDLM

C**** (Simplified) Calendar Related Terms
!@VAR JDperY,JMperY    number of days,months per year
!@VAR JDendOfM(0:12)   last Julian day in month
!@VAR JDmidOfM(0:13)   middle Julian day in month
      INTEGER, PARAMETER :: JDPERY = 365, JMPERY = 12
      INTEGER :: JDendOfM(0:JMPERY) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      INTEGER :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,47,75,106,136,167,197,228,259,289,320,350,381/)
!@VAR AMONTH   (3-4 letter) names for months
      CHARACTER*4 :: AMONTH(0:12) = (/'IC  ',
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)

!@var Q_GISS,Q_HDF,Q_PRT,Q_NETCDF are switches for post-processing
      LOGICAL :: Q_GISS=.FALSE.,Q_HDF=.FALSE.,
     *           Q_PRT =.FALSE.,Q_NETCDF=.FALSE.

C**** IO read/write flags used by the io_xyz routines
!@param IOWRITE,IOREAD Flags set for writing or reading files
!@param IOWRITE_SINGLE Flag used for writing out in single precision
!@param IOWRITE_MON Flag used for writing out monthly rsf file
!@param IRESTART Flag used for reading in normal restart files
!@param IRSFIC Flag used for reading in a restart file I.C.
!@param IRERUN Flag used for reading in a restart file I.C. for rerun
      INTEGER, PARAMETER :: iowrite=-1,ioread=1,iowrite_single=-2
     *     ,iowrite_mon=-3,irestart=1,irsfic=2,irerun=3

C**** Main model prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q,WM
!@var P surface pressure (hecto-Pascals - PTOP)
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P

      END MODULE E001M12_COM

      SUBROUTINE io_model(kunit,it,iaction,ioerr)
!@sum  io_model reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "E001M12"
!@var itime input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var JC1 dummy array
      INTEGER JC1(100)
!@var RC1 dummy array
      REAL*8 RC1(161)
!@var CLABEL1 dummy label
      CHARACTER*156 CLABEL1

C**** Possible additions to this file: FTYPE, (remove rsi from seaice?)
C****  size of common block arrays (and/or should we be explicit?)
C****  timing info as named array?
      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON) ! output to end-of-month restart file
        WRITE (kunit,err=10) it,JC,CLABEL,RC
C**** need a blank line to fool 'qrsfnt' etc.
        WRITE (kunit,err=10)
        WRITE (kunit,err=10) MODULE_HEADER,U,V,T,P,Q,WM
      CASE (IOWRITE_SINGLE)   ! output to single precision ACC file
        WRITE (kunit,err=10) it,JC,CLABEL,SNGL(RC)
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) it,JC1,CLABEL1,RC1
        READ (kunit,err=10)
        READ (kunit,err=10) HEADER,U,V,T,P,Q,WM
        IF (HEADER.ne.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
        SELECT CASE (IACTION)   ! set model common according to iaction
        CASE (IRESTART)         ! normal restart
          JC=JC1 ; CLABEL=CLABEL1 ; RC=RC1
        CASE (IRSFIC)        ! start from old restart file => do nothing
        CASE (IRERUN) ! rerun or extension
          JC=JC1 ; RC=RC1
        END SELECT
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_model

