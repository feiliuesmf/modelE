module MODELE

!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  2009/05/11 (Based originally on B399)
!
! ----------------------------------------------------------
! Restructured modele component
! Common use statements are moved to the begining of the code.
! These modules are used by both NUOPC and non-NUOPC version of the
! modelE code
!
! Fei Liu, 8/26/2013
! ----------------------------------------------------------
!
  USE FILEMANAGER, only : openunit,closeunit
  USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
  USE Dictionary_mod
  Use Parser_mod
  USE MODEL_COM, only: modelEclock, ItimeI, Itime, Ndisk, &
  Jyear0, JMON0, Iyear1, ItimeE, Itime0, &
  NIPRNT, XLABEL, LRUNID, MELSE, Nssw, stop_on, &
  iowrite_single, isBeginningAccumPeriod, &
  KCOPY, NMONAV, IRAND, iowrite_mon, MDIAG, NDAY, &
  rsf_file_name, iowrite, KDISK, dtSRC, MSURF, &
  JDendOfM
  USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,broadcast,sumxpe
  USE RANDOM
  USE GETTIME_MOD
  USE MDIAG_COM, only : monacc,acc_period
#ifdef USE_MPP
  USE fms_mod,         only : fms_init, fms_end
#endif
#ifdef USE_FVCORE
  USE FV_INTERFACE_MOD, only: fvstate
  USE FV_INTERFACE_MOD, only: Checkpoint,Compute_Tendencies
#endif
  use TimeConstants_mod, only: SECONDS_PER_MINUTE, &
                               INT_MONTHS_PER_YEAR
  use TimerPackage_mod, only: startTimer => start
  use TimerPackage_mod, only: stopTimer => stop
  use SystemTimers_mod
  use seaice_com, only : si_ocn,iceocn ! temporary until precip_si,
  use fluxes, only : atmocn,atmice     ! precip_oc calls are moved
  use Month_mod, only: LEN_MONTH_ABBREVIATION

#ifdef USE_ESMF_LIB
  !-----------------------------------------------------------------------------
  ! MODEL Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS    => routine_SetServices, &
    model_label_Advance => label_Advance
  
  implicit none
  private
  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call model_routine_SS(gcomp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    call ESMF_MethodAdd(gcomp, label=model_label_Advance, &
      userRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -> optional Finalize overwrite
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
      userRoutine=NFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
      implicit none
 
      ! Input Arguments
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc
      
      ! Local Arguments
      logical :: qcRestart=.false.
      logical :: coldRestart=.false.
      integer, parameter :: MAX_LEN_IFILE = 32
      character(len=MAX_LEN_IFILE) :: iFile
      !logical, intent(in) :: qcRestart
      !logical, intent(in) :: coldRestart
      !character(len=*), intent(in) :: iFile

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
      integer :: hour, month, day, date, year
      character(len=LEN_MONTH_ABBREVIATION) :: amon

      rc = ESMF_SUCCESS

      call read_options(qcRestart, coldRestart, iFile )

#ifdef USE_SYSUSAGE
      do i_su=0,max_su
        call sysusage(i_su,0)
      enddo
#endif

!****
!**** Reading rundeck (I-file) options
!****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

      call initializeModelE

!****
!**** INITIALIZATIONS
!****
         CALL TIMER (NOW,MDUM)

!**** Read input/ic files
      CALL INPUT (istart,ifile,coldRestart)

!**** Set run_status to "run in progress"
      call write_run_status("Run in progress...",1)

      IF (AM_I_ROOT()) Then
         open(3,file='flagGoStop',form='FORMATTED',status='REPLACE')
         write (3,'(A8)') flg_go
         close (3)
      END IF
      ! FEI: Not sure why a system signal is raised, comment out for now
      ! because the external requirement on sig_stop_model causes compilation
      ! problems.
      !call sys_signal( 15, sig_stop_model )  ! works only on single CPU
      START=NOW
      DO M=1,NTIMEACC
        START= START-TIMING(M)
      END DO

      call modelEclock%getDate(hour=hour, date=date, year=year,amn=amon)

      if (AM_I_ROOT()) &
        WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)') &
        '0NASA/GISS Climate Model (re)started', &
        'Year', year, aMon, date, ', Hr', hour, &
        'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME

         CALL TIMER (NOW,MELSE)

      call sys_flush(6)
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState

    logical :: qcRestart=.false.
    logical :: coldRestart=.false.
    integer, parameter :: MAX_LEN_IFILE = 32
    character(len=MAX_LEN_IFILE) :: iFile
    !logical, intent(in) :: qcRestart
    !logical, intent(in) :: coldRestart
    !character(len=*), intent(in) :: iFile

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
    integer :: hour, month, day, date, year
    character(len=LEN_MONTH_ABBREVIATION) :: amon

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.
    
    call NUOPC_ClockPrintCurrTime(clock, &
      "------>Advancing MODEL from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call NUOPC_ClockPrintStopTime(clock, &
      "--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!****
!**** MAIN LOOP
!****
      call gettime(tloopbegin)
      start9 = (istart == 9)

! Start Main Loop:
      call startTimer('Main Loop')

      if (Ndisk > 0) then
        if (mod(Itime-ItimeI,Ndisk).eq.0 .or. start9) then
         start9 = .false.
         call checkpointModelE()
         call timer(NOW,MELSE)
        END IF
      end if

      if (modelEclock%isBeginningOfDay()) then
        call startNewDay()
      end if

      call atm_phase1

!****
!**** SURFACE INTERACTION AND GROUND CALCULATION
!****
!**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
!**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
!****
!**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      call startTimer('Surface')
      CALL PRECIP_SI(si_ocn,iceocn,atmice)  ! move to ocean_driver
      CALL PRECIP_OC(atmocn,iceocn)         ! move to ocean_driver

!**** CALCULATE SURFACE FLUXES (and, for now, this procedure
!**** also drives "surface" components that are on the atm grid)
      CALL SURFACE
      call stopTimer('Surface')
         CALL TIMER (NOW,MSURF)

      call ocean_driver

! phase 2 changes surf pressure which affects the ocean
      call atm_phase2
!****
!**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
!****
      call modelEclock%nextTick()
      call modelEclock%getDate(year, month, day, date, hour, amon)
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1

      if (modelEclock%isBeginningOfDay()) THEN ! NEW DAY
        months=(year-Jyear0)*INT_MONTHS_PER_YEAR + month - JMON0
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

!****
!**** CALL DIAGNOSTIC ROUTINES
!****
        call startTimer('Diagnostics')

!**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
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

!**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
!**** (after the end of a diagn. accumulation period)
      if (isBeginningAccumPeriod(modelEClock)) then

!**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
        call aPERIOD (JMON0,JYEAR0,months,1,0, aDATE(1:12),Ldate)
        acc_period=aDATE(1:12)
        WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),year
        call print_diags(0)
!**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
!**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.INT_MONTHS_PER_YEAR) m = m - INT_MONTHS_PER_YEAR
            monacc(m) = 1
          end do
          filenm=aDATE(1:7)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
!**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
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

!**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (NOW,MDIAG)
        CALL SUMXPE(TIMING, TIMING_glob, increment=.true.)
        if (am_i_root()) then
          TOTALT=SUM(TIMING_glob(1:NTIMEACC)) ! over all processors
          DO M=1,NTIMEACC
            PERCENT(M) = 100d0*TIMING_glob(M)/(TOTALT+.00001)
          END DO
          TOTALT=SUM(TIMING(1:NTIMEACC)) ! on the root processor
          TOTALT=TOTALT/SECONDS_PER_MINUTE     ! seconds -> minutes
          DTIME = NDAY*TOTALT/(Itime-Itime0) ! minutes/day
          WRITE (6,'(/A,F7.2,A,/(8(A13,F5.1/))//)'), &
              '0TIME',DTIME,'(MINUTES) ', &
              (TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        end if
        TIMING = 0
        START= NOW
        
      END IF  ! beginning of accumulation period

!**** CPU TIME FOR CALLING DIAGNOSTICS
      call stopTimer('Diagnostics')
      CALL TIMER (NOW,MDIAG)

!**** TEST FOR TERMINATION OF RUN
      IF (MOD(Itime,Nssw).eq.0) then
       IF (AM_I_ROOT()) then
        flg_go = '__STOP__'     ! stop if flagGoStop if missing
        iflag=0
        open(3,file='flagGoStop',form='FORMATTED',status='OLD',err=210)
        read (3,'(A8)',end=210) flg_go
        close (3)
 210    continue
        IF (flg_go .eq. '___GO___') iflag=1
        call broadcast(iflag)
       else
        call broadcast(iflag)
        if (iflag .eq. 1) flg_go = '___GO___'
        if (iflag .eq. 0) flg_go = '__STOP__'
       end if
      endif
      IF (flg_go.ne.'___GO___' .or. stop_on) THEN
!**** Flag to continue run has been turned off
         WRITE (6,'("0Flag to continue run has been turned off.")')
         ! FEI:
         ! Somehow needs to stop model advancing if this condition is true
      END IF

      call stopTimer('Main Loop')
!****
!**** END OF MAIN LOOP
!****

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine NFinalize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer :: iu_IFILE
    real*8 :: tloopbegin, tloopend
    integer :: hour, month, day, date, year, ioerr
    character(len=LEN_MONTH_ABBREVIATION) :: amon
    
    rc = ESMF_SUCCESS

    call gettime(tloopend)
    if (AM_I_ROOT()) &
         write(*,*) "Time spent in the main loop in seconds:", &
         tloopend-tloopbegin

!**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
    CALL RFINAL (IRAND)
    call set_param( "IRAND", IRAND, 'o' )
    call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)

    call finalize_atm

    if (AM_I_ROOT()) then
    WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)') &
      '0Restart file written on fort.',KDISK,'Year',year, &
         aMON, date,', Hr',hour,'  Internal clock time:',ITIME
    end if

!**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)

    call printSysTimers()

    IF (AM_I_ROOT()) &
       WRITE (6,'(/////4(1X,33("****")/)//,A,I8 &
                 ///4(1X,33("****")/))') &
      ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME

    IF (Itime.ge.ItimeE) then
       call reportProfile((itimee-itimei)*dtSRC)
       if (AM_I_ROOT()) call print_unused_param(6)
       CALL stop_model ('Terminated normally (reached maximum time)',13)
    END IF

    CALL stop_model ('Run stopped with sswE',12)  ! voluntary stop
#ifdef USE_MPP
    call fms_end( )
#endif
    
  end subroutine

!-----------------------------------------------------------------------------
! Here comes the original version of MODELE code
!-----------------------------------------------------------------------------
#else

  implicit none
  private
  public modelE_mainDriver

  contains

#include "rundeck_opts.h"
      subroutine modelE_mainDriver()
!@sum Acquire configuration options from the command line and pass to
!@+ the model.
!@auth T. Clune
!**** Command line options
      logical :: qcRestart=.false.
      logical :: coldRestart=.false.
      integer, parameter :: MAX_LEN_IFILE = 32
      character(len=MAX_LEN_IFILE) :: iFile

      call read_options(qcRestart, coldRestart, iFile )
      call GISS_modelE(qcRestart, coldRestart, iFile)

      end subroutine modelE_mainDriver

      subroutine GISS_modelE(qcRestart, coldRestart, iFile)
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  2009/05/11 (Based originally on B399)
      implicit none
!**** Command line options
      logical, intent(in) :: qcRestart
      logical, intent(in) :: coldRestart
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
      integer :: hour, month, day, date, year, nstep
      character(len=LEN_MONTH_ABBREVIATION) :: amon

#ifdef USE_SYSUSAGE
      do i_su=0,max_su
        call sysusage(i_su,0)
      enddo
#endif

!****
!**** Reading rundeck (I-file) options
!****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

      call initializeModelE

!****
!**** INITIALIZATIONS
!****
         CALL TIMER (NOW,MDUM)

!**** Read input/ic files
      CALL INPUT (istart,ifile,coldRestart)

!**** Set run_status to "run in progress"
      call write_run_status("Run in progress...",1)

      IF (AM_I_ROOT()) Then
         open(3,file='flagGoStop',form='FORMATTED',status='REPLACE')
         write (3,'(A8)') flg_go
         close (3)
      END IF
      ! FEI: Not sure why a system signal is raised, comment out for now
      ! because the external requirement on sig_stop_model causes compilation
      ! problems.
      !call sys_signal( 15, sig_stop_model )  ! works only on single CPU
      START=NOW
      DO M=1,NTIMEACC
        START= START-TIMING(M)
      END DO

      call modelEclock%getDate(hour=hour, date=date, year=year,amn=amon)

      if (AM_I_ROOT()) &
        WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)') &
        '0NASA/GISS Climate Model (re)started', &
        'Year', year, aMon, date, ', Hr', hour, &
        'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME

         CALL TIMER (NOW,MELSE)

      call sys_flush(6)

!****
!**** MAIN LOOP
!****
      call gettime(tloopbegin)
      start9 = (istart == 9)
      nstep = 0

      main_loop: DO WHILE (Itime.lt.ItimeE)
        call startTimer('Main Loop')
        nstep = nstep + 1

      if (Ndisk > 0) then
        if (mod(Itime-ItimeI,Ndisk).eq.0 .or. start9) then
         start9 = .false.
         call checkpointModelE()
         call timer(NOW,MELSE)
        END IF
      end if

      if (modelEclock%isBeginningOfDay()) then
        call startNewDay()
      end if

      call atm_phase1

!****
!**** SURFACE INTERACTION AND GROUND CALCULATION
!****
!**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
!**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
!****
!**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      call startTimer('Surface')
      CALL PRECIP_SI(si_ocn,iceocn,atmice)  ! move to ocean_driver
      CALL PRECIP_OC(atmocn,iceocn)         ! move to ocean_driver

!**** CALCULATE SURFACE FLUXES (and, for now, this procedure
!**** also drives "surface" components that are on the atm grid)
      CALL SURFACE
      call stopTimer('Surface')
         CALL TIMER (NOW,MSURF)

      call ocean_driver

! phase 2 changes surf pressure which affects the ocean
      call atm_phase2
!****
!**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
!****
      call modelEclock%nextTick()
      call modelEclock%getDate(year, month, day, date, hour, amon)
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1

      if (modelEclock%isBeginningOfDay()) THEN ! NEW DAY
        months=(year-Jyear0)*INT_MONTHS_PER_YEAR + month - JMON0
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

!****
!**** CALL DIAGNOSTIC ROUTINES
!****
        call startTimer('Diagnostics')

!**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
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

!**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
!**** (after the end of a diagn. accumulation period)
      if (isBeginningAccumPeriod(modelEClock)) then

!**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
        call aPERIOD (JMON0,JYEAR0,months,1,0, aDATE(1:12),Ldate)
        acc_period=aDATE(1:12)
        WRITE (aDATE(8:14),'(A3,I4.4)') aMON(1:3),year
        call print_diags(0)
!**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
!**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.INT_MONTHS_PER_YEAR) m = m - INT_MONTHS_PER_YEAR
            monacc(m) = 1
          end do
          filenm=aDATE(1:7)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
!**** KCOPY > 1 : ALSO SAVE THE RESTART INFORMATION
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

!**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (NOW,MDIAG)
        CALL SUMXPE(TIMING, TIMING_glob, increment=.true.)
        if (am_i_root()) then
          TOTALT=SUM(TIMING_glob(1:NTIMEACC)) ! over all processors
          DO M=1,NTIMEACC
            PERCENT(M) = 100d0*TIMING_glob(M)/(TOTALT+.00001)
          END DO
          TOTALT=SUM(TIMING(1:NTIMEACC)) ! on the root processor
          TOTALT=TOTALT/SECONDS_PER_MINUTE     ! seconds -> minutes
          DTIME = NDAY*TOTALT/(Itime-Itime0) ! minutes/day
          WRITE (6,'(/A,F7.2,A,/(8(A13,F5.1/))//)'), &
              '0TIME',DTIME,'(MINUTES) ', &
              (TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        end if
        TIMING = 0
        START= NOW
        
      END IF  ! beginning of accumulation period

!**** CPU TIME FOR CALLING DIAGNOSTICS
      call stopTimer('Diagnostics')
      CALL TIMER (NOW,MDIAG)

!**** TEST FOR TERMINATION OF RUN
      IF (MOD(Itime,Nssw).eq.0) then
       IF (AM_I_ROOT()) then
        flg_go = '__STOP__'     ! stop if flagGoStop if missing
        iflag=0
        open(3,file='flagGoStop',form='FORMATTED',status='OLD',err=210)
        read (3,'(A8)',end=210) flg_go
        close (3)
 210    continue
        IF (flg_go .eq. '___GO___') iflag=1
        call broadcast(iflag)
       else
        call broadcast(iflag)
        if (iflag .eq. 1) flg_go = '___GO___'
        if (iflag .eq. 0) flg_go = '__STOP__'
       end if
      endif
      IF (flg_go.ne.'___GO___' .or. stop_on) THEN
!**** Flag to continue run has been turned off
         WRITE (6,'("0Flag to continue run has been turned off.")')
         EXIT main_loop
      END IF

      call stopTimer('Main Loop')
      END DO main_loop
!****
!**** END OF MAIN LOOP
!****

      call gettime(tloopend)
      if (AM_I_ROOT()) &
           write(*,*) "Time spent in the main loop in seconds:", &
           tloopend-tloopbegin
      if (AM_I_ROOT()) &
           write(*,*) "Total number of time steps:", &
           nstep

!**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)

      call finalize_atm

      if (AM_I_ROOT()) then
      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)') &
        '0Restart file written on fort.',KDISK,'Year',year, &
           aMON, date,', Hr',hour,'  Internal clock time:',ITIME
      end if

!**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)

      call printSysTimers()

      IF (AM_I_ROOT()) &
         WRITE (6,'(/////4(1X,33("****")/)//,A,I8 &
                   ///4(1X,33("****")/))') &
        ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME

      IF (Itime.ge.ItimeE) then
         call reportProfile((itimee-itimei)*dtSRC)
         if (AM_I_ROOT()) call print_unused_param(6)
         CALL stop_model ('Terminated normally (reached maximum time)',13)
      END IF

      CALL stop_model ('Run stopped with sswE',12)  ! voluntary stop
#ifdef USE_MPP
      call fms_end( )
#endif

      end subroutine GISS_modelE

! USE_ESMF_LIB macro endif
#endif
!-----------------------------------------------------------------------------
! End of NUOPC and non-NUOPC specific code
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
! COMMON subroutines used by both NUOPC and non-NUOPC version of MODELE
!
!-----------------------------------------------------------------------------

      subroutine read_options(qcRestart, coldRestart, iFile )
!@sum Reads options from the command line
!@auth I. Aleinov
      implicit none
!@var qcRestart true if "-r" is present
!@var iFile is name of the file containing run configuration data
      logical, intent(inout) :: qcRestart
      logical, intent(inout) :: coldRestart
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
        case ("-cold-restart")
          coldRestart = .true.
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

      print *, 'Configuration file name: ', iFile

      return
      end subroutine read_options

      subroutine initializeModelE
      USE DOMAIN_DECOMP_1D, ONLY : init_app, setCommunicator
#ifdef USE_ESMF_LIB
      type(ESMF_VM)              :: vm
      integer                    :: mpicom, rc
#endif

      call initializeSysTimers()

#ifdef USE_MPP
      call fms_init( )
#endif
#ifdef USE_ESMF_LIB
#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort(__LINE__,rc)
      call ESMF_VMGetCurrent(vm, rc=rc)
      VERIFY_(rc)
      call ESMF_VMGet(vm, mpicommunicator=mpicom, rc=rc)
      VERIFY_(rc)
      call setCommunicator(mpicom)
#else
      call init_app()
#endif
      call initializeDefaultTimers()

      call alloc_drv_atm()
      call alloc_drv_ocean()

      end subroutine initializeModelE

      subroutine startNewDay()
      use model_com, only: modelEclock
!**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      logical :: newmonth
      integer :: month, day

      month = modelEclock%month()
      day = modelEclock%dayOfYear()
      newmonth = (day == 1+JDendOfM(month-1))
      call daily_DIAG(newmonth) ! atmosphere
      if(newmonth) then         ! ocean
        call reset_ODIAG(0)
#ifndef STANDALONE_OCEAN
        call reset_glaacc
#endif
      endif
      end subroutine startNewDay

!TODO fv, fv_fname, and fv_dfname are  not yet passed as arguments
!TODO exist except when building an FV version
      subroutine checkpointModelE
!@sum Every Ndisk Time Steps (DTsrc), starting with the first one,
!@+ write restart information alternately onto 2 disk files
      use MODEL_COM, only: rsf_file_name,kdisk,irand
      use MODEL_COM, only: itime
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: Checkpoint,fvstate
#endif
      
      integer :: hour, date, year, ioerr
      character(len=LEN_MONTH_ABBREVIATION) :: amon

      call modelEclock%getDate(hour=hour, date=date, amn=amon)

      CALL rfinal(IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)
#if defined( USE_FVCORE )
      call checkpoint(fvstate, rsf_file_name(KDISK))
#endif
      if (AM_I_ROOT()) &
          WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)') &
          '0Restart file written on fort.',KDISK,'Year', &
          year,aMon,date,', Hr',hour,'  Internal clock time:',ITIME
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
      use TimeConstants_mod, only: INT_MINUTES_PER_DAY

#ifdef USE_MPI
      include 'mpif.h'
#endif
      real*8, intent(in) :: elapsedTimeInSeconds
      character(len=MAX_RECORD_LENGTH), pointer :: lines(:)
      integer :: i

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
      call setScale(column, INT_MINUTES_PER_DAY/elapsedTimeInSeconds, &
          'min/day')
      call addColumn(report, column)

      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, INT_MINUTES_PER_DAY/elapsedTimeInSeconds, &
           'min/day')
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

#ifdef USE_MPI
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

      subroutine dailyUpdates
      use fluxes, only : atmocn
      implicit none
      
      call daily_CAL(.true.)    ! end_of_day
      call daily_OCEAN(.true.,atmocn)  ! end_of_day
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
      USE MODEL_COM, only : NIPRNT,master_yr &
          ,NMONAV,Ndisk,Nssw,KCOPY,KOCEAN,IRAND,ItimeI
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE Dictionary_mod
      implicit none

!**** Rundeck parameters:
      call sync_param( "NMONAV", NMONAV )
      call sync_param( "NIPRNT", NIPRNT )
      call sync_param( "Ndisk", Ndisk )
      call sync_param( "Nssw", Nssw )
      call sync_param( "KCOPY", KCOPY )
      call sync_param( "KOCEAN", KOCEAN )
      call sync_param( "IRAND", IRAND )
      if (is_set_param("master_yr")) then
        call get_param( "master_yr", master_yr )
      else
        call stop_model('Please define master_yr in the rundeck.',255)
      endif
      RETURN
!****
      end subroutine init_Model

      SUBROUTINE INPUT (istart,ifile,coldRestart)

!****
!**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
!**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
!****
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : timing,ntimeacc
      USE Dictionary_mod
      USE MODEL_COM, only : &
            xlabel,lrunid,nmonav,qcheck,irand &
           ,nday,dtsrc,kdisk,jmon0,jyear0 &
           ,iyear1,itime,itimei,itimee &
           ,idacc,modelEclock &
           ,aMONTH,jdendofm,aMON0 &
           ,ioread,irerun,irsfic &
           ,melse,Itime0,Jdate0 &
           ,Jhour0,rsf_file_name &
           ,HOURI,DATEI,MONTHI,YEARI ,HOURE,DATEE,MONTHE,YEARE &
           ,iwrite_sv,jwrite_sv,itwrite_sv,kdiag_sv
      USE RANDOM
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
#ifndef STANDALONE_OCEAN
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE RESOLUTION, only : LM ! atm reference for init_tracer hack
#endif
#endif

      use TimeConstants_mod, only : SECONDS_PER_DAY, INT_HOURS_PER_DAY, &
                                    INT_DAYS_PER_YEAR
      use ModelClock_mod, only: ModelClock, newModelClock
      use Time_mod, only: Time, newTime
      use Calendar_mod, only: Calendar
      use JulianCalendar_mod, only: makeJulianCalendar
      use Month_mod, only: LEN_MONTH_ABBREVIATION

      IMPLICIT NONE
!@var istart  postprocessing(-1)/start(1-8)/restart(>8)  option
      integer, intent(out) :: istart
      character(*), intent(in) :: ifile
      logical, intent(in) :: coldRestart
!@dbparam init_topog_related : set = 1 if IC and topography are incompatible
      integer :: init_topog_related = 0
!@dbparam do_IC_fixups : set = 1 if surface IC are to be checked/corrected
      integer :: do_IC_fixups = 0
!@var iu_AIC,iu_IFILE unit numbers for input files
      INTEGER iu_AIC,iu_IFILE
      INTEGER I,J,L,K,LID1,LID2,NOFF,ioerr
!@nlparam IHRI,TIMEE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
!@nlparam IRANDI  random number seed to perturb init.state (if>0)
      INTEGER :: IHRI=-1,TIMEE=-1,IHOURE=-1,IRANDI=0
      INTEGER IhrX, KDISK_restart
      LOGICAL :: is_coldstart
      CHARACTER NLREC*80,RLABEL*132

      INTEGER :: IWRITE=0,JWRITE=0,ITWRITE=23
      INTEGER, DIMENSION(13) :: KDIAG

      NAMELIST/INPUTZ/ ISTART,IRANDI &
           ,IWRITE,JWRITE,ITWRITE,QCHECK,KDIAG &
           ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1 &
!****    List of parameters that are disregarded at restarts
           ,        HOURI,DATEI,MONTHI,YEARI
      NAMELIST/INPUTZ_cold/ ISTART,IRANDI &
           ,IWRITE,JWRITE,ITWRITE,QCHECK,KDIAG &
           ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1 &
!****    List of parameters that are disregarded at restarts
           ,        HOURI,DATEI,MONTHI,YEARI
      integer istart_fixup
      character*132 :: bufs
      integer, parameter :: MAXLEN_RUNID = 32

      type (Time) :: modelETimeI
      class (Calendar), pointer :: pCalendar
      integer :: hour, month, day, date, year
      character(len=LEN_MONTH_ABBREVIATION) :: amon

!****
!**** Default setting for ISTART : restart from latest save-file (10)
!****
      ISTART=10

!**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
      KDIAG(13)=9

!**** Set global default timing descriptions
!**** Other speciality descriptions can be added/used locally
      NTIMEACC = 0

!****
!**** Print Header and Label (2 lines) from rundeck
!****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      if (AM_I_ROOT()) WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(iu_IFILE,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      if (AM_I_ROOT()) WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label

!****
!**** Print preprocessing options (if any are defined)
!****
      IF(AM_I_ROOT()) call print_and_check_PPopts

!****
!**** Read parameters from the rundeck to database and namelist
!****
      !call parse_params(iu_IFILE)
      ! skip "&&PARAMETERS" section
      do
        read( iu_IFILE, *, err=910, end=910 ) bufs
        if ( bufs == '&&END_PARAMETERS' ) exit
      enddo

      READ (iu_IFILE,NML=INPUTZ,ERR=900)
      if (coldRestart) READ (iu_IFILE,NML=INPUTZ_cold,ERR=900)
      call closeunit(iu_IFILE)

      IWRITE_sv = IWRITE
      JWRITE_sv = JWRITE
      ITWRITE_sv = ITWRITE
      KDIAG_sv = KDIAG

      if (istart.le.0) then
        call stop_model('pdE not supported',255)
      end if

!**** Get those parameters which are needed in this subroutine
      call get_param( "DTsrc", DTsrc )
      if(is_set_param("IRAND"))  call get_param( "IRAND", IRAND )

      if (istart.lt.9) then
!***********************************************************************
!****                                                               ****
!****                  INITIAL STARTS - ISTART: 2, 8                ****
!****                                                               ****
!****   Current settings: 2 - from observed data                    ****
!****                     8 - from current model M-file - no resets ****
!****                                                               ****
!***********************************************************************
!****
!**** Set quantities that are derived from the namelist parameters
!****
!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc later if necessary
      NDAY = 2*NINT(.5*SECONDS_PER_DAY/DTsrc)

!**** Get Start Time; at least YearI HAS to be specified in the rundeck
      IF (YearI.lt.0) then
        IF (AM_I_ROOT()) &
         WRITE(6,*) 'Please choose a proper start year yearI, not',yearI
        call stop_model('INPUT: yearI not provided',255)
      END IF
      IF (Iyear1.lt.0) Iyear1 = yearI
      IhrI = HourI + INT_HOURS_PER_DAY*(dateI-1 + JDendofM(monthI-1) + &
             INT_DAYS_PER_YEAR*(yearI-Iyear1))
      ITimeI = IhrI*NDAY/INT_HOURS_PER_DAY ! internal clock counts DTsrc-steps
      Itime=ItimeI
      IF (IhrI.lt.0) then
        IF (AM_I_ROOT()) &
        WRITE(6,*) 'Improper start time OR Iyear1=',Iyear1,' > yearI;', &
        ' yearI,monthI,dateI,hourI=',yearI,monthI,dateI,hourI
        call stop_model('INPUT: Improper start date or base year Iyear1',255)
      END IF

      IF (ISTART.EQ.2) THEN
!****
!**** Cold Start: ISTART=2
!****
        XLABEL(1:80)='Observed atmospheric data from NMC tape'

!**** Set flag to initialise topography-related variables
        init_topog_related = 1

      ELSE IF (ISTART==8) THEN
!****
!****   Data from current type of RESTART FILE
!****
      ! no need to read SRHR,TRHR,FSF,TSFREZ,diag.arrays
        call io_rsf("AIC",IhrX,irsfic,ioerr)

!**** Check consistency of starting time
        IF( (MOD(IHRI-IHRX,INT_HOURS_PER_DAY*INT_DAYS_PER_YEAR).ne.0) )  &
        THEN
         WRITE (6,*) ' Difference in hours between ', &
             'Starting date and Data date:', &
             MOD(IHRI-IHRX,INT_HOURS_PER_DAY*INT_DAYS_PER_YEAR)
         WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
         call stop_model('INPUT: start date inconsistent with data',255)
        ENDIF
      END IF

!**** Set flags to initialise some variables related to topography
      call sync_param( "init_topog_related", init_topog_related )

      IF (init_topog_related == 1) then
        do_IC_fixups = 1        ! new default, not necessarily final
      ENDIF

      IF (AM_I_ROOT()) &
           WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)') &
        '0Model started on',datei,aMONTH(monthi),yeari,' Hr',houri, &
        'ISTART =',ISTART,XLABEL(1:80)    ! report input file label
      XLABEL = RLABEL                     ! switch to rundeck label

      else ! initial versus restart
!***********************************************************************
!****                                                               ****
!****                  RESTARTS: ISTART > 8                         ****
!****                                                               ****
!****   Current settings: 9 - from own model M-file                 ****
!****                    10 - from later of fort.1 or fort.2        ****
!****                    11 - from fort.1                           ****
!****                    12 - from fort.2                           ****
!****               13 & up - from earlier of fort.1 or fort.2      ****
!****                                                               ****
!***********************************************************************
!****
!****   DATA FROM end-of-month RESTART FILE     ISTART=9
!****        mainly used for REPEATS and delayed EXTENSIONS
      IF(ISTART==9) THEN                     !  diag.arrays are not read in
        call io_rsf("AIC",Itime,irerun,ioerr)
        WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=', &
          ISTART,', TIME=',Itime,' ',XLABEL(1:80) ! sho input file label
        XLABEL = RLABEL                        ! switch to rundeck label
        TIMING = 0
      ELSE
!****
!**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
!****
!**** CHOOSE DATA SET TO RESTART ON
        IF(ISTART==11 .OR. ISTART==12) THEN
          KDISK=ISTART-10
        ELSEIF(ISTART==10 .OR. ISTART==13) THEN
          call find_later_rsf(kdisk)
          IF (ISTART.GE.13)     KDISK=3-KDISK
        ENDIF
        call io_rsf(rsf_file_name(KDISK),Itime,ioread,ioerr)
        KDISK_restart = KDISK
        if (AM_I_ROOT()) &
            WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT', &
            KDISK,', Time=',Itime,' ',XLABEL(1:80)

!**** Switch KDISK if the other file is (or may be) bad (istart>10)
!****     so both files will be fine after the next write execution
        IF (istart.gt.10) KDISK=3-KDISK
!**** Keep KDISK after reading from the later restart file, so that
!****     the same file is overwritten first; in case of trouble,
!****     the earlier restart file will still be available

      ENDIF

      endif ! initial versus restart

!***********************************************************************
!****                                                              *****
!****       INITIAL- AND RESTARTS: Final Initialization steps      *****
!****                                                              *****
!***********************************************************************

!**** initialize Lrunid (length of the identifying part of XLABEL)
!****
      lid1 = INDEX(XLABEL,'(') -1
      if (lid1.lt.1) lid1=MAXLEN_RUNID+1
      lid2 = INDEX(XLABEL,' ') -1
      if (lid2.lt.1) lid2=MAXLEN_RUNID+1
      LRUNID = min(lid1,lid2)
      IF (LRUNID.gt.MAXLEN_RUNID) call stop_model &
           ('INPUT: Rundeck name too long. Shorten to 32 char or less' &
           ,255)

!**** Update ItimeE only if YearE or IhourE is specified in the rundeck
!****
      if(timee.lt.0) timee=houre*nday/INT_HOURS_PER_DAY
      IF(yearE.ge.0) ItimeE = (( (yearE-iyear1)*INT_DAYS_PER_YEAR +  &
                    JDendofM(monthE-1) + dateE-1) * INT_HOURS_PER_DAY) * &
                    NDAY/INT_HOURS_PER_DAY + TIMEE
!**** Alternate (old) way of specifying end time
      if(IHOURE.gt.0) ItimeE=IHOURE*NDAY/INT_HOURS_PER_DAY

!**** Check consistency of DTsrc with NDAY
      if (is_set_param("DTsrc") .and. nint(SECONDS_PER_DAY/DTsrc) &
          .ne. NDAY) then
        if (AM_I_ROOT()) then
          write(6,*) 'DTsrc=',DTsrc,' has to stay at/be set to',  &
                     SECONDS_PER_DAY/NDAY
        end if
        call stop_model('INPUT: DTsrc inappropriately set',255)
      end if
      DTsrc = SECONDS_PER_DAY/NDAY
      call set_param( "DTsrc", DTsrc, 'o' )   ! copy DTsrc into DB

!**** NMONAV has to be 1(default),2,3,4,6,12, i.e. a factor of 12
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      if (NMONAV.lt.1 .or. MOD(12,NMONAV).ne.0) then
        write (6,*) 'NMONAV has to be 1,2,3,4,6 or 12, not',NMONAV
        call stop_model('INPUT: nmonav inappropriately set',255)
      end if
      if (AM_I_ROOT()) &
           write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'

!**** Updating Parameters: If any of them changed beyond this line
!**** use set_param(.., .., 'o') to update them in the database (DB)

!**** Get the rest of parameters from DB or put defaults to DB

      call init_Model

!**** Set julian date information
      call getdte(Itime,Nday,Iyear1,year,month,day,date,hour,amon)
      call getdte(Itime0,Nday,iyear1,Jyear0,Jmon0,J,Jdate0,Jhour0,amon0)

      pCalendar => makeJulianCalendar()
      modelETimeI = newTime(pCalendar)
      call getdte(Itime,Nday,Iyear1,year,month,day,date,hour,amon)
      call modelEtimeI%setByDate(year, month, date, hour)
      modelEclock = newModelClock(modelEtimeI,itime,Nday)

      CALL DAILY_cal(.false.)                  ! not end_of_day

#ifndef STANDALONE_OCEAN
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
!**** Initialise tracer parameters and diagnostics
!**** MUST be before other init routines
!**** TODO: split init_tracer into general definitions and
!**** component-specific ops, folding the latter into component inits
      if(istart.eq.2) call read_nmc()  ! hack, see TODO
      CALL CALC_AMPK(LM)               ! hack
      call init_tracer
#endif
#endif

!!! hack: may be prevented if post-processing option is eliminated
      istart_fixup = istart
      if (istart==8 .and. do_IC_fixups==1) istart_fixup = 9

      is_coldstart = (istart<9 .and. init_topog_related == 1)
! long version: 
!      is_coldstart = istart==2 .or. (istart==8 .and. init_topog_related == 1)

      call INPUT_ocean (istart,istart_fixup, & 
           do_IC_fixups,is_coldstart)

      call INPUT_atm(istart,istart_fixup, &
           do_IC_fixups,is_coldstart, &
           KDISK_restart,IRANDI)

      if (istart.le.9) then
        call reset_adiag(0)
        call reset_odiag(0)
#ifndef STANDALONE_OCEAN
        call reset_glaacc
#endif
      endif

      CALL SET_TIMER("       OTHER",MELSE)

      if (AM_I_ROOT()) then
         WRITE (6,INPUTZ)
         call print_param( 6 )
         WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
      end if

!****
      RETURN
!****
!**** TERMINATE BECAUSE OF IMPROPER PICK-UP
!****
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
      call stop_model( &
      ' Water tracers need TRACERS_ON as well as TRACERS_WATER',255)
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
#ifdef TRACERS_AEROSOLS_VBS
      write(6,*) '   with VBS organics'
#endif
#endif
#ifdef TRACERS_AMP
      write(6,*) '...and aerosol microphysics'
#endif
#ifdef TRACERS_TOMAS
      write(6,*) '...and TOMAS aerosol microphysics'
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
#ifdef MERRA_NUDGING
      write(6,*) '   WITH MERRA WINDS.'
#endif
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

end module
