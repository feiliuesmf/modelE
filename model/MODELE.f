#include "rundeck_opts.h"
      subroutine modelE_mainDriver()
!@sum Acquire configuration options from the command line and pass to
!@+ the model.
!@auth T. Clune
C**** Command line options
      logical :: qcRestart=.false.
      integer, parameter :: MAX_LEN_IFILE = 32
      character(len=MAX_LEN_IFILE) :: iFile

      call read_options(qcRestart, iFile )
      call GISS_modelE(qcRestart, iFile)

      contains

      subroutine read_options(qcRestart, iFile )
!@sum Reads options from the command line
!@auth I. Aleinov
      implicit none
!@var qcRestart true if "-r" is present
!@var iFile is name of the file containing run configuration data
      logical, intent(inout) :: qcRestart
      character(*),intent(out)  :: ifile
      integer, parameter :: MAX_LEN_ARG = 80
      character(len=MAX_LEN_ARG) :: arg, value

      iFile = "";
      do
        call nextarg( arg, 1 )
        if ( arg == "" ) exit          ! end of args
        select case (arg)
        case ("-r")
          qcRestart = .true.
        case ("-i")
          call nextarg( value, 0 )
          iFile=value
        ! new options can be included here
        case default
          print *,'Unknown option specified: ', arg
          print *,'Aborting...'
          call stop_model("Unknown option on a command line",255)
        end select
      enddo

      if (iFile == "") then
        print*, 'No configuration file specified on command line: '
        print*, 'Aborting ...'
        call stop_model("No configuration file on command line.",255)
      end if

      return
      end subroutine read_options

      end subroutine modelE_mainDriver

      subroutine GISS_modelE(qcRestart, iFile)
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  2009/05/11 (Based originally on B399)
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE Dictionary_mod
      Use Parser_mod
      USE MODEL_COM
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,ESMF_BCAST,sumxpe
      USE RANDOM
      USE GETTIME_MOD
      USE DIAG_COM, only : MODD5S
      USE MDIAG_COM, only : monacc,acc_period
      USE DIAG_SERIAL, only : print_diags
#ifdef USE_MPP
      USE fms_mod,         only : fms_init, fms_end
#endif
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: fvstate
      USE FV_INTERFACE_MOD, only: Checkpoint,Compute_Tendencies
#endif
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use SystemTimers_mod
      implicit none
C**** Command line options
      logical, intent(in) :: qcRestart
      character(len=*), intent(in) :: iFile

      INTEGER K,M,MSTART,MNOW,months,ioerr,Ldate,istart
      INTEGER :: MDUM = 0

      character(len=80) :: filenm

      REAL*8, DIMENSION(NTIMEMAX) :: PERCENT
      REAL*8, DIMENSION(0:NTIMEMAX) ::TIMING_glob
      REAL*8 start,now, DTIME,TOTALT

      CHARACTER aDATE*14
      CHARACTER*8 :: flg_go='___GO___'      ! green light
      integer :: iflag=1
      external sig_stop_model
      logical :: start9

      integer :: iu_IFILE
      real*8 :: tloopbegin, tloopend

#ifdef USE_SYSUSAGE
      do i_su=0,max_su
        call sysusage(i_su,0)
      enddo
#endif

C****
C**** Reading rundeck (I-file) options
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

      call initializeModelE

C****
C**** INITIALIZATIONS
C****
         CALL TIMER (NOW,MDUM)

C**** Read input/ic files
      CALL INPUT (istart,ifile)

C**** Set run_status to "run in progress"
      call write_run_status("Run in progress...",1)

      IF (AM_I_ROOT()) Then
         open(3,file='flagGoStop',form='FORMATTED',status='REPLACE')
         write (3,'(A8)') flg_go
         close (3)
      END IF
      call sys_signal( 15, sig_stop_model )  ! works only on single CPU
      START=NOW
      DO M=1,NTIMEACC
        START= START-TIMING(M)
      END DO

      if (AM_I_ROOT())
     *   WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)')
     *   '0NASA/GISS Climate Model (re)started',
     *   'Year',JYEAR,aMON,JDATE,', Hr',JHOUR,
     *   'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME

         CALL TIMER (NOW,MELSE)

      call sys_flush(6)

C****
C**** MAIN LOOP
C****
      call gettime(tloopbegin)
      start9 = (istart == 9)

      main_loop: DO WHILE (Itime.lt.ItimeE)
        call startTimer('Main Loop')

      if (Ndisk > 0) then
        if (mod(Itime-ItimeI,Ndisk).eq.0 .or. start9) then
         start9 = .false.
         call checkpointModelE()
         call timer(NOW,MELSE)
        END IF
      end if
      
      if (isBeginningOfDay(modelEclock)) then
        call startNewDay()
      end if

      call atm_phase1

C****
C**** SURFACE INTERACTION AND GROUND CALCULATION
C****
C**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
C**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
C****
C**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      call startTimer('Surface')
      CALL PRECIP_SI('OCEAN')  ! move to ocean_driver
      CALL PRECIP_OC           ! move to ocean_driver

C**** CALCULATE SURFACE FLUXES (and, for now, this procedure
C**** also drives "surface" components that are on the atm grid)
      CALL SURFACE
      call stopTimer('Surface')
         CALL CHECKT ('SURFACE')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (5)

      call ocean_driver

! phase 2 changes surf pressure which affects the ocean
      call atm_phase2

C****
C**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
C****
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1
      Jhour=MOD(Itime*24/NDAY,24)         ! Hour (0-23)

      if (isBeginningOfDay(modelEclock)) THEN ! NEW DAY
        months=(Jyear-Jyear0)*JMperY + JMON-JMON0
        call startTimer('Daily')
        call dailyUpdates
        call TIMER (NOW,MELSE)
        call stopTimer('Daily')
      end if                                  !  NEW DAY
       
#ifdef USE_FVCORE
! Since dailyUpdates currently adjusts surf pressure,
! moving this call to the atm driver will change results.
! todo 2: fold this into the fv run procedure.
       Call Compute_Tendencies(fvstate)
#endif

C****
C**** CALL DIAGNOSTIC ROUTINES
C****
        call startTimer('Diagnostics')

C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
      IF (NIPRNT.GT.0) THEN
        acc_period='PARTIAL      '
#ifdef NEW_IO
        filenm='PARTIAL.acc'//XLABEL(1:LRUNID)
        call io_rsf (filenm,Itime,iowrite_single,ioerr)
#endif
        call print_diags(1)
        NIPRNT=NIPRNT-1
        call set_param( "NIPRNT", NIPRNT, 'o' )
      END IF

C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
C**** (after the end of a diagn. accumulation period)
      if (isBeginningAccumPeriod(modelEclock)) then

C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
        call aPERIOD (JMON0,JYEAR0,months,1,0, aDATE(1:12),Ldate)
        acc_period=aDATE(1:12)
        WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),JYEAR
        call print_diags(0)
C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.12) m = m-12
            monacc(m) = 1
          end do
          filenm=aDATE(1:7)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
C**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
          IF (KCOPY.GT.1) THEN
            CALL RFINAL (IRAND)
            call set_param( "IRAND", IRAND, 'o' )
            filenm='1'//aDATE(8:14)//'.rsf'//XLABEL(1:LRUNID)
            call io_rsf(filenm,Itime,iowrite_mon,ioerr)
#if defined( USE_FVCORE )
            call Checkpoint(fvstate, filenm)
#endif
          END IF
        END IF

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (NOW,MDIAG)
        CALL SUMXPE(TIMING, TIMING_glob, increment=.true.)
        if (am_i_root()) then
          TOTALT=SUM(TIMING_glob(1:NTIMEACC)) ! over all processors
          DO M=1,NTIMEACC
            PERCENT(M) = 100d0*TIMING_glob(M)/(TOTALT+.00001)
          END DO
          TOTALT=SUM(TIMING(1:NTIMEACC)) ! on the root processor
          TOTALT=TOTALT/60.     ! seconds -> minutes
          DTIME = NDAY*TOTALT/(Itime-Itime0) ! minutes/day
          WRITE (6,'(/A,F7.2,A,/(8(A13,F5.1/))//)')
     *         '0TIME',DTIME,'(MINUTES) ',
     *         (TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        end if
        TIMING = 0
        START= NOW
        
      END IF  ! beginning of accumulation period

C**** CPU TIME FOR CALLING DIAGNOSTICS
      call stopTimer('Diagnostics')
      CALL TIMER (NOW,MDIAG)

C**** TEST FOR TERMINATION OF RUN
      IF (MOD(Itime,Nssw).eq.0) then
       IF (AM_I_ROOT()) then
        flg_go = '__STOP__'     ! stop if flagGoStop if missing
        iflag=0
        open(3,file='flagGoStop',form='FORMATTED',status='OLD',err=210)
        read (3,'(A8)',end=210) flg_go
        close (3)
 210    continue
        IF (flg_go .eq. '___GO___') iflag=1
        call ESMF_BCAST( iflag)
       else
        call ESMF_BCAST( iflag)
        if (iflag .eq. 1) flg_go = '___GO___'
        if (iflag .eq. 0) flg_go = '__STOP__'
       end if
      endif
      IF (flg_go.ne.'___GO___' .or. stop_on) THEN
C**** Flag to continue run has been turned off
         WRITE (6,'("0Flag to continue run has been turned off.")')
         EXIT main_loop
      END IF

      call stopTimer('Main Loop')
      END DO main_loop
C****
C**** END OF MAIN LOOP
C****

      call gettime(tloopend)
      if (AM_I_ROOT())
     *     write(*,*) "Time spent in the main loop in seconds:",
     *     tloopend-tloopbegin

C**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)

      call finalize_atm

      if (AM_I_ROOT()) then
      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *  '0Restart file written on fort.',KDISK,'Year',JYEAR,
     *     aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
      end if

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)

      call printSysTimers()

      IF (AM_I_ROOT())
     *   WRITE (6,'(/////4(1X,33("****")/)//,A,I8
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME

      IF (Itime.ge.ItimeE) then
         call reportProfile((itimee-itimei)*dtSRC)
         if (AM_I_ROOT()) call print_unused_param(6)
         CALL stop_model (
     &     'Terminated normally (reached maximum time)',13)
      END IF

      CALL stop_model ('Run stopped with sswE',12)  ! voluntary stop
#ifdef USE_MPP
      call fms_end( )
#endif

      contains

      subroutine initializeModelE
      USE DOMAIN_DECOMP_1D, ONLY : init_app

      call initializeSysTimers()

#ifdef USE_MPP
      call fms_init( )
#endif
      call init_app()
      call initializeDefaultTimers()

      call alloc_drv_atm()
      call alloc_drv_ocean()

      end subroutine initializeModelE

      subroutine startNewDay()
C**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      logical :: newmonth
      newmonth = (JDAY == 1+JDendOfM(Jmon-1))
      call daily_DIAG(newmonth) ! atmosphere
      if(newmonth) then         ! ocean
        call reset_ODIAG(0)
      endif
      end subroutine startNewDay

!TODO fv, fv_fname, and fv_dfname are  not yet passed as arguments
!TODO exist except when building an FV version
      subroutine checkpointModelE
!@sum Every Ndisk Time Steps (DTsrc), starting with the first one,
!@+ write restart information alternately onto 2 disk files
      use MODEL_COM, only: rsf_file_name,kdisk,irand
      use MODEL_COM, only: Jyear, aMon, Jdate, Jhour, itime
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: Checkpoint,fvstate
#endif
      
      CALL rfinal(IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)
#if defined( USE_FVCORE )
      call checkpoint(fvstate, rsf_file_name(KDISK))
#endif
      if (AM_I_ROOT())
     *     WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *     '0Restart file written on fort.',KDISK,'Year',
     *     JYEAR,aMON,JDATE,', Hr',JHOUR,'  Internal clock time:',ITIME
      kdisk=3-kdisk

      end subroutine checkpointModelE

      subroutine initializeDefaultTimers()
      use TimerPackage_mod, only: initialize
      use TimerPackage_mod, only: addTimer

      call initialize()
      call addTimer('Main Loop')
      call addTimer(' Atm. Dynamics')
      call addTimer(' MELT_SI()')
      call addTimer(' CONDSE()')
      call addTimer(' RADIA()')
      call addTimer(' Surface')
      call addTimer('  Precip')
      call addTimer('  SURFACE()')
      call addTimer(' DYNSI()')
      call addTimer(' UNDERICE()')
      call addTimer(' GROUND_SI()')
      call addTimer(' GROUND_LI()')
      call addTimer(' GROUND_LK()')
      call addTimer(' RIVERF()')
      call addTimer(' OCEANS')
      call addTimer(' Daily')
      call addTimer(' Diagnostics')

      end subroutine initializeDefaultTimers

      subroutine reportProfile(elapsedTimeInSeconds)
      use TimerPackage_mod
      use DOMAIN_DECOMP_1D, only: AM_I_ROOT

#if  defined(USE_ESMF) || defined(USE_MPP)
      include 'mpif.h'
#endif
      real*8, intent(in) :: elapsedTimeInSeconds
      character(len=MAX_RECORD_LENGTH), pointer :: lines(:)
      integer :: i

      integer, parameter :: SECONDS_PER_MINUTE = 60
      integer, parameter :: MINUTES_PER_DAY = 24*60
      type (ProfileReport_type) :: report
      type (ReportColumn_type) :: column
      real*8 :: totalTime

      call finalize()
      report = newReport()
      call addColumn(report, newColumn(NAME_COLUMN, fieldWidth=20))

      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=7)
      totalTime = getInclusiveTime(getTimer('main'))
      call setPrecision(column, 2)
      call setScale(column, 100/totalTime, ' %')
      call addColumn(report, column)

      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, MINUTES_PER_DAY/elapsedTimeInSeconds,
     &     'min/day')
      call addColumn(report, column)

      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, MINUTES_PER_DAY/elapsedTimeInSeconds,
     &     'min/day')
      call addColumn(report, column)

      call addColumn(report, newColumn(TRIP_COUNTS_COLUMN,fieldWidth=6))

      column = newColumn(MAXIMUM_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)

      call addColumn(report, newColumn(MAX_PROCESS_COLUMN,fieldWidth=4))

      column = newColumn(MINIMUM_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)
      call addColumn(report, newColumn(MIN_PROCESS_COLUMN,fieldWidth=4))

      column = newColumn(AVERAGE_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)

#if  defined(USE_ESMF) || defined(USE_MPP)
      lines => generateParallelReport(report, MPI_COMM_WORLD)
#else
      lines => generateReport(report)
#endif

      if (AM_I_ROOT()) then
         do i = 1, size(lines)
            write(*,'(a)') trim(lines(i))
         end do
      end if
      deallocate(lines)
      call delete(report)

      end subroutine reportProfile
      
      end subroutine GISS_modelE

      subroutine dailyUpdates
      implicit none
      
      call daily_CAL(.true.)    ! end_of_day
      call daily_OCEAN(.true.)  ! end_of_day
      call daily_ATM(.true.)

      return
      end subroutine dailyUpdates

      subroutine sig_stop_model
      USE MODEL_COM, only : stop_on
      implicit none
      stop_on = .true.
      end subroutine sig_stop_model


      subroutine init_Model
!@sum This program reads most of parameters from the database (DB)
!@+   get_param( "A", X ) reads parameter A into variable X
!@+   if "A" is not in the database, it will generate an error
!@+   message and stop
!@+   sync_param( "B", Y ) reads parameter B into variable Y
!@+   if "B" is not in the database, then Y is unchanged and its
!@+   value is saved in the database as "B" (here sync = synchronize)
      USE ATM_COM, only : ij_debug
      USE MODEL_COM, only : NIPRNT
     *     ,NMONAV,Ndisk,Nssw,KCOPY,KOCEAN,IRAND,ItimeI
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE Dictionary_mod
      use FLUXES, only : UOdrag,NIsurf
      implicit none

C**** Rundeck parameters:
      call sync_param( "NMONAV", NMONAV )
      call sync_param( "NIPRNT", NIPRNT )
      call sync_param( "Ndisk", Ndisk )
      call sync_param( "Nssw", Nssw )
      call sync_param( "KCOPY", KCOPY )
      call sync_param( "KOCEAN", KOCEAN )
      call sync_param( "NIsurf", NIsurf )
      call sync_param( "UOdrag", UOdrag )
      call sync_param( "IRAND", IRAND )
      call sync_param( "ij_debug",ij_debug , 2)

      RETURN
C****
      end subroutine init_Model

      SUBROUTINE INPUT (istart,ifile)

C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : timing,ntimeacc
      USE Dictionary_mod
      USE CONSTANT, only : sday,hrday
      USE MODEL_COM, only :
     *      xlabel,lrunid,nmonav,qcheck,irand
     *     ,nday,dtsrc,kdisk,jmon0,jyear0
     *     ,iyear1,itime,itimei,itimee
     *     ,idacc,jyear,jmon,jday,jdate,jhour
     *     ,aMONTH,jdendofm,jdpery,aMON,aMON0
     *     ,ioread,irerun,irsfic
     *     ,melse,Itime0,Jdate0
     *     ,Jhour0,rsf_file_name
     *     ,HOURI,DATEI,MONTHI,YEARI ,HOURE,DATEE,MONTHE,YEARE
     &     ,iwrite_sv,jwrite_sv,itwrite_sv,kdiag_sv
      USE RANDOM
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE RESOLUTION, only : LM ! atm reference for init_tracer hack
#endif
      IMPLICIT NONE
!@var istart  postprocessing(-1)/start(1-8)/restart(>8)  option
      integer, intent(out) :: istart
!@dbparam init_topog_related : set = 1 if IC and topography are incompatible
      integer :: init_topog_related = 0
!@dbparam do_IC_fixups : set = 1 if surface IC are to be checked/corrected
      integer :: do_IC_fixups = 0
      character(*), intent(in) :: ifile
!@var iu_AIC,iu_IFILE unit numbers for input files
      INTEGER iu_AIC,iu_IFILE
      INTEGER I,J,L,K,LID1,LID2,NOFF,ioerr
!@param INTHRDAY number of hours in a day in integer form
      INTEGER, PARAMETER :: INTHRDAY=INT(HRDAY)
!@nlparam IHRI,TIMEE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
!@nlparam IRANDI  random number seed to perturb init.state (if>0)
      INTEGER :: IHRI=-1,TIMEE=-1,IHOURE=-1,IRANDI=0
      INTEGER IhrX, KDISK_restart
      LOGICAL :: is_coldstart
      CHARACTER NLREC*80,RLABEL*132

      INTEGER :: IWRITE=0,JWRITE=0,ITWRITE=23
      INTEGER, DIMENSION(13) :: KDIAG

      NAMELIST/INPUTZ/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,KDIAG
     *     ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI
      integer istart_fixup
      character*132 :: bufs
      integer, parameter :: MAXLEN_RUNID = 32

C****
C**** Default setting for ISTART : restart from latest save-file (10)
C****
      ISTART=10

C**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
      KDIAG(13)=9

C**** Set global default timing descriptions
C**** Other speciality descriptions can be added/used locally
      NTIMEACC = 0

C****
C**** Print Header and Label (2 lines) from rundeck
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      if (AM_I_ROOT()) WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(iu_IFILE,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      if (AM_I_ROOT()) WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label

C****
C**** Print preprocessing options (if any are defined)
C****
      IF(AM_I_ROOT()) call print_and_check_PPopts

C****
C**** Read parameters from the rundeck to database and namelist
C****
      !call parse_params(iu_IFILE)
      ! skip "&&PARAMETERS" section
      do
        read( iu_IFILE, *, err=910, end=910 ) bufs
        if ( bufs == '&&END_PARAMETERS' ) exit
      enddo

      READ (iu_IFILE,NML=INPUTZ,ERR=900)
      call closeunit(iu_IFILE)

      IWRITE_sv = IWRITE
      JWRITE_sv = JWRITE
      ITWRITE_sv = ITWRITE
      KDIAG_sv = KDIAG

      if (istart.le.0) then
        call stop_model('pdE not supported',255)
      end if

C**** Get those parameters which are needed in this subroutine
      call get_param( "DTsrc", DTsrc )
      if(is_set_param("IRAND"))  call get_param( "IRAND", IRAND )

      if (istart.lt.9) then
C***********************************************************************
C****                                                               ****
C****                  INITIAL STARTS - ISTART: 2, 8                ****
C****                                                               ****
C****   Current settings: 2 - from observed data                    ****
C****                     8 - from current model M-file - no resets ****
C****                                                               ****
C***********************************************************************
C****
C**** Set quantities that are derived from the namelist parameters
C****
!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc later if necessary
      NDAY = 2*NINT(.5*SDAY/DTsrc)

C**** Get Start Time; at least YearI HAS to be specified in the rundeck
      IF (YearI.lt.0) then
        IF (AM_I_ROOT())
     *   WRITE(6,*) 'Please choose a proper start year yearI, not',yearI
        call stop_model('INPUT: yearI not provided',255)
      END IF
      IF (Iyear1.lt.0) Iyear1 = yearI
      IhrI = HourI +
     +  INTHRDAY*(dateI-1 + JDendofM(monthI-1) + JDperY*(yearI-Iyear1))
      ITimeI = IhrI*NDAY/INTHRDAY ! internal clock counts DTsrc-steps
      Itime=ItimeI
      IF (IhrI.lt.0) then
        IF (AM_I_ROOT())
     *  WRITE(6,*) 'Improper start time OR Iyear1=',Iyear1,' > yearI;',
     *  ' yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
        call stop_model(
     &       'INPUT: Improper start date or base year Iyear1',255)
      END IF

      IF (ISTART.EQ.2) THEN
C****
C**** Cold Start: ISTART=2
C****
        XLABEL(1:80)='Observed atmospheric data from NMC tape'

C**** Set flag to initialise topography-related variables
        init_topog_related = 1

      ELSE IF (ISTART==8) THEN
C****
C****   Data from current type of RESTART FILE
C****
      ! no need to read SRHR,TRHR,FSF,TSFREZ,diag.arrays
        call io_rsf("AIC",IhrX,irsfic,ioerr)

C**** Check consistency of starting time
        IF( (MOD(IHRI-IHRX,8760).ne.0) ) THEN
         WRITE (6,*) ' Difference in hours between ',
     *       'Starting date and Data date:',MOD(IHRI-IHRX,8760)
         WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
         call stop_model('INPUT: start date inconsistent with data',255)
        ENDIF
      END IF

C**** Set flags to initialise some variables related to topography
      call sync_param( "init_topog_related", init_topog_related )

      IF (init_topog_related == 1) then
        do_IC_fixups = 1        ! new default, not necessarily final
      ENDIF

      IF (AM_I_ROOT())
     *     WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
     *  '0Model started on',datei,aMONTH(monthi),yeari,' Hr',houri,
     *  'ISTART =',ISTART,XLABEL(1:80)    ! report input file label
      XLABEL = RLABEL                     ! switch to rundeck label

      else ! initial versus restart
C***********************************************************************
C****                                                               ****
C****                  RESTARTS: ISTART > 8                         ****
C****                                                               ****
C****   Current settings: 9 - from own model M-file                 ****
C****                    10 - from later of fort.1 or fort.2        ****
C****                    11 - from fort.1                           ****
C****                    12 - from fort.2                           ****
C****               13 & up - from earlier of fort.1 or fort.2      ****
C****                                                               ****
C***********************************************************************
C****
C****   DATA FROM end-of-month RESTART FILE     ISTART=9
C****        mainly used for REPEATS and delayed EXTENSIONS
      IF(ISTART==9) THEN                     !  diag.arrays are not read in
        call io_rsf("AIC",Itime,irerun,ioerr)
        WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *    ISTART,', TIME=',Itime,' ',XLABEL(1:80) ! sho input file label
        XLABEL = RLABEL                        ! switch to rundeck label
        TIMING = 0
      ELSE
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
C****
C**** CHOOSE DATA SET TO RESTART ON
        IF(ISTART==11 .OR. ISTART==12) THEN
          KDISK=ISTART-10
        ELSEIF(ISTART==10 .OR. ISTART==13) THEN
          call find_later_rsf(kdisk)
          IF (ISTART.GE.13)     KDISK=3-KDISK
        ENDIF
        call io_rsf(rsf_file_name(KDISK),Itime,ioread,ioerr)
        KDISK_restart = KDISK
        if (AM_I_ROOT())
     *      WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
     *      KDISK,', Time=',Itime,' ',XLABEL(1:80)

C**** Switch KDISK if the other file is (or may be) bad (istart>10)
C****     so both files will be fine after the next write execution
        IF (istart.gt.10) KDISK=3-KDISK
C**** Keep KDISK after reading from the later restart file, so that
C****     the same file is overwritten first; in case of trouble,
C****     the earlier restart file will still be available

      ENDIF

      endif ! initial versus restart

C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************

C**** initialize Lrunid (length of the identifying part of XLABEL)
C****
      lid1 = INDEX(XLABEL,'(') -1
      if (lid1.lt.1) lid1=MAXLEN_RUNID+1
      lid2 = INDEX(XLABEL,' ') -1
      if (lid2.lt.1) lid2=MAXLEN_RUNID+1
      LRUNID = min(lid1,lid2)
      IF (LRUNID.gt.MAXLEN_RUNID) call stop_model
     *     ('INPUT: Rundeck name too long. Shorten to 32 char or less'
     *     ,255)

C**** Update ItimeE only if YearE or IhourE is specified in the rundeck
C****
      if(timee.lt.0) timee=houre*nday/INTHRDAY
      IF(yearE.ge.0) ItimeE = (( (yearE-iyear1)*JDperY +
     *  JDendofM(monthE-1)+dateE-1 )*INTHRDAY )*NDAY/INTHRDAY + TIMEE
C**** Alternate (old) way of specifying end time
      if(IHOURE.gt.0) ItimeE=IHOURE*NDAY/INTHRDAY

C**** Check consistency of DTsrc with NDAY
      if (is_set_param("DTsrc") .and. nint(sday/DTsrc).ne.NDAY) then
        if (AM_I_ROOT())
     *        write(6,*) 'DTsrc=',DTsrc,' has to stay at/be set to',SDAY/NDAY
        call stop_model('INPUT: DTsrc inappropriately set',255)
      end if
      DTsrc = SDAY/NDAY
      call set_param( "DTsrc", DTsrc, 'o' )   ! copy DTsrc into DB

C**** NMONAV has to be 1(default),2,3,4,6,12, i.e. a factor of 12
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      if (NMONAV.lt.1 .or. MOD(12,NMONAV).ne.0) then
        write (6,*) 'NMONAV has to be 1,2,3,4,6 or 12, not',NMONAV
        call stop_model('INPUT: nmonav inappropriately set',255)
      end if
      if (AM_I_ROOT())
     *     write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'

C**** Updating Parameters: If any of them changed beyond this line
C**** use set_param(.., .., 'o') to update them in the database (DB)

C**** Get the rest of parameters from DB or put defaults to DB
      call init_Model

C**** Set julian date information
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      call getdte(Itime0,Nday,iyear1,Jyear0,Jmon0,J,Jdate0,Jhour0,amon0)

      CALL DAILY_cal(.false.)                  ! not end_of_day

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** Initialise tracer parameters and diagnostics
C**** MUST be before other init routines
C**** TODO: split init_tracer into general definitions and
C**** component-specific ops, folding the latter into component inits
      if(istart.eq.2) call read_nmc()  ! hack, see TODO
      CALL CALC_AMPK(LM)               ! hack
      call init_tracer
#endif

!!! hack: may be prevented if post-processing option is eliminated
      istart_fixup = istart
      if (istart==8 .and. do_IC_fixups==1) istart_fixup = 9

      is_coldstart = (istart<9 .and. init_topog_related == 1)
! long version: 
!      is_coldstart = istart==2 .or. (istart==8 .and. init_topog_related == 1)

      call INPUT_ocean (istart,istart_fixup,
     &     do_IC_fixups,is_coldstart)

      call INPUT_atm(istart,istart_fixup,is_coldstart,
     &     KDISK_restart,IRANDI)

      if (istart.le.9) then
        call reset_adiag(0)
        call reset_odiag(0)
      endif

      CALL SET_TIMER("       OTHER",MELSE)

      if (AM_I_ROOT()) then
         WRITE (6,INPUTZ)
         call print_param( 6 )
         WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
      end if

C****
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
  900 write (6,*) 'Error in NAMELIST parameters'
      call stop_model('Error in NAMELIST parameters',255)
  910 write (6,*) 'Error readin I-file'
      call stop_model('Error reading I-file',255)
      END SUBROUTINE INPUT

      subroutine print_and_check_PPopts
!@sum prints preprocessor options in english and checks some
!@+  interdependencies. (moved from subroutine INPUT).
!@+  Called by root thread only.

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      write(6,*) 'This program includes tracer code'
#endif
#ifdef TRACERS_WATER
      write(6,*) '...and water tracer code'
#ifndef TRACERS_ON
      call stop_model(
     &' Water tracers need TRACERS_ON as well as TRACERS_WATER',255)
#endif
#endif
#ifdef TRACERS_OCEAN
      write(6,*) '...and ocean tracer code'
#endif
#ifdef TRACERS_SPECIAL_O18
      write(6,*) '...and water isotope code'
#ifndef TRACERS_WATER
      call stop_model('Water isotope tracers need TRACERS_WATER '//
     *'as well as TRACERS_SPECIAL_O18',255)
#endif
#endif
#ifdef TRACERS_SPECIAL_Lerner
      write(6,*) '...and Jean/David tracers and chemistry'
#endif
#ifdef TRACERS_GASEXCH_ocean
      write(6,*) '          '
      write(6,*) '...and Natassa Romanou air-sea GAS EXCHANGE'
#ifdef TRACERS_OceanBiology
      write(6,*) '          '
      write(6,*) '...and Natassa Romanou/Watson Gregg ocean biology '
#endif
#ifdef TRACERS_GASEXCH_ocean_CFC
      write(6,*) '****CFC flux across air/sea interface****'
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      write(6,*) '****CO2 flux across air/sea interface****'
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Drew Shindell tracers and chemistry'
#endif
#ifdef TRACERS_TERP
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Terpenes tracer'
#else
      call stop_model('Terpenes tracer needs tropo chemistry',255)
#endif
#endif  /* TRACERS_TERP */
#ifdef CALCULATE_FLAMMABILITY
      write(6,*) '...and calculating sfc veg flammability'
#endif
#ifdef CALCULATE_LIGHTNING
      write(6,*) '...and calculating lightning flash rate'
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
      write(6,*) '...and dynamic biomass burning srcs by flammability'
#endif
#ifdef TRACERS_AEROSOLS_Koch
      write(6,*) '...and Dorothy Koch aerosols'
#endif
#ifdef TRACERS_AMP
      write(6,*) '...and aerosol microphysics'
#endif
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and secondary organic aerosols'
#else
      call stop_model('SOA version needs tropo chemistry',255)
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef SOA_DIAGS
#ifdef TRACERS_AEROSOLS_SOA
      write(6,*) '...and additional SOA diagnostics'
#else
      call stop_model('SOA_DIAGS needs TRACERS_AEROSOLS_SOA',255)
#endif  /* TRACERS_AEROSOLS_SOA */
#endif  /* SOA_DIAGS */
#ifdef TRACERS_AEROSOLS_OCEAN
      write(6,*) '...and oceanic organic aerosol sources'
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_DRYDEP
      write(6,*) '...and tracer dry deposition'
#endif
#ifdef EDGAR_HYDE_SOURCES
      write(6,*) '...and EDGAR HYDE sources instead of GISS'
#endif
#ifdef SHINDELL_STRAT_EXTRA
      write(6,*) '...and Drew Shindell extra strat tracers'
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      write(6,*) '...and interactive CH4 wetlands emissions'
#endif
#ifdef NUDGE_ON
      write(6,*) '...and nudging of meteorology'
#endif
#ifdef HTAP_LIKE_DIAGS
      write(6,*) '...and HTAP set of diagnostics'
#endif
#ifdef ACCMIP_LIKE_DIAGS
      write(6,*) '...and ACCMIP set of diagnostics'
#ifndef SHINDELL_STRAT_EXTRA
      call stop_model
     & ('SHINDELL_STRAT_EXTRA should be on for ACCMIP_LIKE_DIAGS',255)
#endif
#ifndef HTAP_LIKE_DIAGS
      call stop_model
     & ('HTAP_LIKE_DIAGS should be on for ACCMIP_LIKE_DIAGS',255)
#endif
#endif /* ACCMIP_LIKE_DIAGS */

      return
      end subroutine print_and_check_PPopts

