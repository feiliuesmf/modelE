      MODULE MODEL_COM
!@sum  MODEL_COM Main model variables, independent of resolution
!@auth Original Development Team
!@ver  1.0
      USE RESOLUTION, only : im,jm,lm
      IMPLICIT NONE
      SAVE

!@var IMH half the number of latitudinal boxes
      INTEGER, PARAMETER :: IMH=IM/2
!@var FIM,BYIM real parameter values for the number of long. grid boxes
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM
!@var JEQ grid box immediately north of the equator
      INTEGER, PARAMETER :: JEQ=1+JM/2

C**** THERE ARE 100 INTEGER PARAMETERS IN COMMON (JC-ARRAY)
c      pointer KOCEAN, LS1, NRAD, NIsurf, NFILTR, MFILTR, NDAA
c     *     ,NDA5D,NDA5K,NDA5S,NDA4,NDASF,IRAND,Nslp,Kvflxo,KCOPY,Ndisk
c     *     ,NSSW,KEYCT,NIPRNT,NMONAV
      INTEGER ::
     *  IM0,JM0,LM0,LS1,KACC0,        KTACC0,Itime,ItimeI,ItimeE,Itime0,
     *  KOCEAN,KDISK,KEYCT,KCOPY,IRAND,  MFILTR,Ndisk,Kvflxo,Nslp,NIdyn,
     *  NRAD,NIsurf,NFILTR,NDAY,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR0,JYEAR,JYEAR0,JMON,JMON0, JDATE,JDATE0,JHOUR,JHOUR0,JDAY,
     *  NSSW,NSTEP,MRCH,NIPRNT,NMONAV
      INTEGER, DIMENSION(32) :: IDUM
      INTEGER, DIMENSION(2,4) :: IJD6
      INTEGER, DIMENSION(:), pointer :: IDACC ! dim = 12

! handle for referring to integer parameters
ccc JC doesh't contain any data and will be removed soon
      INTEGER, DIMENSION(100) :: JC
      !EQUIVALENCE (JC,IM0)

C**** THERE ARE 161 REAL NUMBERS IN COMMON (RC-ARRAY)
c      pointer PTOP, PSF, DTsrc, DT, SKIPSE
      DOUBLE PRECISION ::
     *  DTsrc,DT,  PTOP,PSF,PSFMPT,PSTRAT,PSDRAG, SKIPSE
      DOUBLE PRECISION, DIMENSION(4) :: TAUTR0
      DOUBLE PRECISION, DIMENSION(LM) :: SIG
      DOUBLE PRECISION, DIMENSION(LM+1) :: SIGE
      DOUBLE PRECISION, DIMENSION(161-13-2*LM) :: RDM2
!@var PSFMPT,PSTRAT derived pressure constants

!@var RC handle for referring to real parameters
ccc RC doesh't contain any data and will be removed soon
      DOUBLE PRECISION, DIMENSION(161) :: RC
      !EQUIVALENCE (RC,DTsrc)

c      pointer NAMD6
      CHARACTER*4 NAMD6(4),AMON,AMON0
      CHARACTER*132 XLABEL

!@var LABEL1,CLABEL,XLABEL handles for referring to text parameters
ccc CLABEL doesh't contain any data and will be removed soon
      !CHARACTER LABEL1*16
      CHARACTER CLABEL*156
      !EQUIVALENCE (CLABEL,XLABEL,LABEL1)

      DOUBLE PRECISION, DIMENSION(IM,JM) :: FLAND,FOCEAN,FLICE,FLAKE0
     *     ,FEARTH,ZATMO,HLAKE

      DOUBLE PRECISION, DIMENSION(IM,JM,11) :: VDATA
      DOUBLE PRECISION, DIMENSION(IM,JM) :: WFCS

!@var DSIG,BYDSIG,DSIGO sigma level cooridinates
      DOUBLE PRECISION, DIMENSION(LM) :: DSIG,BYDSIG
      DOUBLE PRECISION, DIMENSION(LM-1) :: DSIGO

!@var XCDLM.  SDRAG ~XCDLM(1)+XCDLM(2)*wind_magnitude
      DOUBLE PRECISION, DIMENSION(2) :: XCDLM = (/5.D-4,5.D-5/)

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

!@var aturb_on switches diffus on/off, drycnv off/on
      LOGICAL :: aturb_on=.true.

!@var vt_on switches on/off surface and turbulence temperature virtual.
c**** not working yet for EARTH
      LOGICAL :: vt_on=.true.

C**** IO read/write flags used by the io_xyz routines
!@param IOWRITE Flag used for writing normal restart files
!@param IOWRITE_SINGLE Flag used for saving diags in single precision
!@param IOWRITE_MON Flag used for saving restart part only (no diags)
!@param IOREAD Flag used for reading in (composite) restart files
!@param IRSFIC Flag used for reading in restart part to start NEW run
!@param IRERUN Flag used for reading in restart part to extend OLD run
      INTEGER, PARAMETER :: ioread=1,ioread_single=2,irerun=3,irsfic=4
     *    ,iowrite=-1,iowrite_single=-2,iowrite_mon=-3

C**** Main model prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q,WM
!@var P surface pressure (hecto-Pascals - PTOP)
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P

C**** Define surface types (mostly used for weighting diagnostics)
!@param NTYPE number of different surface types
      INTEGER, PARAMETER :: NTYPE=6   ! orig = 3
!@var FTYPE fractions of each surface type
      REAL*8, DIMENSION(NTYPE,IM,JM) :: FTYPE
!@var ITxx indices of various types (used only when it matters)
      INTEGER, PARAMETER :: ITOCEAN=1, ITOICE=2, ITEARTH=3,
     *                      ITLANDI=4, ITLAKE=5, ITLKICE=6

!@var AIRX, AIRMX*DXYP(J) Used in stratosphere model. (kg?)
      REAL*8 AIRX(IM,JM)
!@var LMC max layer of mc convective mass flux. (Strat model)
      INTEGER, DIMENSION(2,IM,JM) :: LMC

      DATA IM0,JM0,LM0, KACC0/         ! KTACC0 should be here too ???
     *     IM ,JM ,LM , -131313 /,   ! KACC0 - should be removed !!!
     *  KOCEAN,KDISK,KEYCT,KCOPY,     IRAND,MFILTR,Ndisk,Kvflxo,Nslp/
     *       1,    1,    1,    2, 123456789,     1,   24,   0,   0/,
     *  Nrad, Nfiltr, NIsurf, Nssw, NIPRNT, NMONAV,         IYEAR0/
     *     5,      2,      2,    1,      1,      1,           1976/,
     *  NDAa,   NDA5d, NDA5k, NDA5s, NDA4, NDAsf/
     *     7,       7,     7,     7,   24,     1/,
     *  MODRD,MODD5K,MODD5S/
     *      0,     0,     0/
      DATA  DT, DTsrc/
     *    450., 3600./,
C****
C**** Note:           DT = DTdyn and NIdyn = DTsrc/DTdyn (set in INPUT)
C**** In general      DTxxx = Nxxx*DTsrc  and  DTxxx = DTsrc/NIxxx
C**** except that the time steps related to NDAa, NDA5k, NDAsf are
C**** slightly larger:     NDAa:   NDAa*DTsrc + 2*DT(dyn),
C****                      NDA5k: NDA5k*DTsrc + 2*DT(dyn),
C****                      NDAsf: NDAsf*DTsrc + DTsrc/NIsurf
C****
     *  PTOP, PSF, PSDRAG,SKIPSE/
     *  150.,984.,  500.,     0./
      DATA SIGE /1.0000000,LM*0./                    ! Define in rundeck
      DATA NAMD6 /'AUSD','MWST','SAHL','EPAC'/,
     *  IJD6/63,17, 17,34, 37,27, 13,23/

      END MODULE MODEL_COM

      MODULE TIMINGS
!@sum  TIMINGS contains variables for keeping track of computing time
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@param NTIMEMAX maximum number of possible time accumulators
      INTEGER, PARAMETER :: NTIMEMAX=10
!@var NTIMEACC actual number of time accumulators
      INTEGER :: NTIMEACC = 0
!@var TIMING array that holds timing info
      INTEGER, DIMENSION(0:NTIMEMAX) :: TIMING
!@var TIMESTR array that holds timing info description
      CHARACTER*12, DIMENSION(NTIMEMAX) :: TIMESTR

      END MODULE TIMINGS

      SUBROUTINE SET_TIMER(STR,MINDEX)
!@sum  SET_TIMER sets an index of TIMING for a particular description
!@auth Gavin Schmidt
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
!@var STR string that describes timing accumulator
      CHARACTER*12, INTENT(IN) :: STR
!@var MINDEX index for that accumulator
      INTEGER, INTENT(OUT) :: MINDEX
      INTEGER N

C**** Check whether index has been set
      DO N=1,NTIMEACC
        IF (STR.EQ.TIMESTR(N)) THEN
          MINDEX=N
          RETURN
        END IF
      END DO
C**** Otherwise increase number of indexes
      NTIMEACC = NTIMEACC + 1
      IF (NTIMEACC.gt.NTIMEMAX)
     *     STOP "Too many timing indices: increase NTIMEMAX"
      MINDEX = NTIMEACC
      TIMESTR(MINDEX) = STR
C****
      RETURN
      END SUBROUTINE SET_TIMER

      SUBROUTINE TIMER (MNOW,MSUM)
!@sum  TIMER keeps track of elapsed CPU time in hundredths of seconds
!@auth Gary Russell
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW   !@var MNOW current CPU time (.01 s)
      INTEGER, INTENT(INOUT) :: MSUM !@var MSUM index for running total
      INTEGER :: MINC                !@var MINC time since last call
      INTEGER, SAVE :: MLAST = 0     !@var MLAST  last CPU time

      CALL GETTIME(MNOW)
      MINC  = MNOW - MLAST
      TIMING(MSUM)  = TIMING(MSUM) + MINC
      MLAST = MNOW
      RETURN
      END SUBROUTINE TIMER

      SUBROUTINE TIMEOUT (MBEGIN,MIN,MOUT)
!@sum  TIMEOUT redistributes timing info between counters
!@auth Gary Russell
!@ver  1.0
      USE TIMINGS
      IMPLICIT NONE
!@var MBEGIN CPU time start of section (.01 s)
      INTEGER, INTENT(IN) :: MBEGIN
      INTEGER, INTENT(INOUT) :: MIN  !@var MIN index to be added to
      INTEGER, INTENT(INOUT) :: MOUT !@var MOUT index to be taken from
      INTEGER :: MINC                !@var MINC time since MBEGIN
      INTEGER :: MNOW                !@var MNOW current CPU time (.01 s)

      CALL GETTIME(MNOW)
      MINC  = MNOW - MBEGIN
      TIMING(MIN)  = TIMING(MIN)  + MINC
      TIMING(MOUT) = TIMING(MOUT) - MINC
      RETURN
      END SUBROUTINE TIMEOUT

      SUBROUTINE io_label(kunit,it,iaction,ioerr)
!@sum  io_label reads and writes label/parameters to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      USE TIMINGS, only : ntimemax,ntimeacc,timestr,timing
      USE PARAM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var itime input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var XLABEL1 dummy label
      CHARACTER*132 XLABEL1
!@var NTIM1,TSTR1,TIM1 timing related dummy arrays
      INTEGER NTIM1,TIM1(NTIMEMAX)
      CHARACTER*12 TSTR1(NTIMEMAX)
!@var ITmin,ITmax minimal/maximal Itime of acc files to be summed up
      INTEGER, SAVE :: ITmin=999999, ITmax=-1, IT0min
      INTEGER it1,it0,idacc0(12) !@var it1,it0,idacc1 dummy variables
      INTEGER, DIMENSION(11) :: idind = (/1,2,3,4,6,7,8,9,10,11,12/)

C**** Possible additions to this file: FTYPE, (remove rsi from seaice?)
C****  size of common block arrays (and/or should we be explicit?)
C****  timing info as named array?
      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to end-of-month restart file
        !WRITE (kunit,err=10) it,JC,CLABEL,RC,
        WRITE (kunit,err=10) it,XLABEL,
     *       NTIMEACC,TIMING(1:NTIMEACC),TIMESTR(1:NTIMEACC)
C**** need a blank line to fool 'qrsfnt' etc. (to be dropped soon)
        WRITE (kunit,err=10)
C**** write parameters database here
        call write_param(kunit)
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) it,XLABEL1,
     *        NTIM1,TIM1(1:NTIM1),TSTR1(1:NTIM1)
C**** need a blank line to fool 'qrsfnt' etc. (to be dropped soon)
        READ (kunit,err=10)
        SELECT CASE (IACTION)   ! set model common according to iaction
        CASE (ioread) ! parameters/label from restart file
          call read_param(kunit,.true.)
          XLABEL=XLABEL1
          NTIMEACC=NTIM1
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIM1(1:NTIM1)
        CASE (IRSFIC)       ! use defaults, rundeck label
          read(kunit,err=10) ! skip parameters
         ! switch 'it' to 'ihour' using 'nday' of restart file ?????
        CASE (IRERUN)  ! params: rsfile, label: rundeck
          call read_param(kunit,.false.)
        CASE (IOREAD_SINGLE) ! parameters/label from restart file
                             ! accumulate idacc and keep track of Itmax/min
          IDACC0=IDACC
          XLABEL=XLABEL1
          NTIMEACC=NTIM1
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIMING(1:NTIM1)+TIM1(1:NTIM1)
          call read_param(kunit,.true.)

C**** keep track of min/max time and earliest diagnostic period
          call get_param("Itime",it1)
          call get_param("Itime0",it0)
          if (it1.gt.ITmax) ITmax=it1
          if (it1.lt.ITmin) THEN
            ITmin=it1
            IT0min=it0
          end if
          call set_param("Itime",ITmax,'o')
          call set_param("Itime0",IT0min,'o')

C**** This should probably be moved to io_diags
          IDACC(idind) = IDACC(idind) + IDACC0(idind)
          IF (IDACC0(5).gt.0) IDACC(5) = MIN(IDACC(5),IDACC0(5))

        END SELECT ! namelist parameters may still be changed in rundeck
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_label

       SUBROUTINE io_model(kunit,iaction,ioerr)
!@sum  io_model reads and writes model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "MODEL01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,U,V,T,P,Q,WM
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) HEADER,U,V,T,P,Q,WM
        IF (HEADER.ne.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_model

      subroutine getdte(It,Nday,Iyr0,Jyr,Jmn,Jd,Jdate,Jhour,amn)
!@sum  getdte gets julian calander info from internal timing info
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : hrday
      USE MODEL_COM, only : JDperY,JDendOfM,amonth
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: It,Nday,Iyr0
      INTEGER, INTENT(OUT) :: Jyr,Jmn,Jd,Jdate,Jhour
      CHARACTER*4, INTENT(OUT) :: amn

      Jyr=Iyr0+It/(Nday*JDperY)
      Jd=1+It/Nday-(Jyr-Iyr0)*JDperY
      Jmn=1
      do while (Jd.GT.JDendOfM(Jmn))
        Jmn=Jmn+1
      end do
      Jdate=Jd-JDendOfM(Jmn-1)
      Jhour=nint(mod(It*hrday/Nday,hrday))
      amn=amonth(Jmn)

      return
      end subroutine getdte
