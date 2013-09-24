#ifdef USE_ESMF_LIB
module modele_utils

! Provides basic utilities for the NUOPC modele driver to initialize
! its clock.

  use ESMF

  implicit none
  private

  public initialize_modele_clock

contains

      ! Called by modele_driver to initialize nuopc modele driver clock
      subroutine initialize_modele_clock(gcomp, startTime, stopTime, timeStep, rc)

        ! Arguments
        type(ESMF_GridComp),     intent(inout)     :: gcomp
        type(ESMF_Time),         intent(inout)     :: startTime
        type(ESMF_Time),         intent(inout)     :: stopTime
        type(ESMF_TimeInterval), intent(inout)     :: timeStep
        integer, optional,       intent(out)       :: rc

        ! Local Variables
        logical                                    :: qcRestart=.false.
        logical                                    :: coldRestart=.false.
        integer                                    :: istart
        character(len=32)                          :: ifile

        if(present(rc)) rc = ESMF_SUCCESS

        call read_options(qcRestart, coldRestart, iFile )

        call ESMF_AttributeSet(gcomp, name="iFile", value=iFile, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_AttributeSet(gcomp, name="coldRestart", value=coldRestart, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_AttributeSet(gcomp, name="qcRestart", value=qcRestart, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call initialize_clocks(istart, iFile, coldRestart, &
          startTime, stopTime, timeStep)

      end subroutine initialize_modele_clock

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

      SUBROUTINE initialize_clocks(istart,ifile,coldRestart, &
        startTime, stopTime, timeStep)

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
      integer,                 intent(out)       :: istart
      character(*),            intent(in)        :: ifile
      logical,                 intent(in)        :: coldRestart
      type(ESMF_Time),         intent(inout)     :: startTime
      type(ESMF_Time),         intent(inout)     :: stopTime
      type(ESMF_TimeInterval), intent(inout)     :: timeStep
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
      integer :: hour, month, day, date, year, rc
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
      if(is_set_param("DTsrc"))  call get_param( "DTsrc", DTsrc )
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

!**** Set julian date information
      call getdte(Itime,Nday,Iyear1,year,month,day,date,hour,amon)
      call getdte(Itime0,Nday,iyear1,Jyear0,Jmon0,J,Jdate0,Jhour0,amon0)

      pCalendar => makeJulianCalendar()
      modelETimeI = newTime(pCalendar)
      call getdte(Itime,Nday,Iyear1,year,month,day,date,hour,amon)
      call modelEtimeI%setByDate(year, month, date, hour)
      modelEclock = newModelClock(modelEtimeI,itime,Nday)

!**** Initialize ESMF clocks needed by the driver
      if (AM_I_ROOT()) then
        write(6,*) 'itimei=',itimei, 'itimee=', itimee, 'dtsrc=', dtsrc
      end if

      call getdte(Itimei,Nday,Iyear1,year,month,day,date,hour,amon)
      call ESMF_TimeSet(startTime, yy=year, mm=month, d=date, h=hour, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_TimePrint(startTime)

      call getdte(Itimee,Nday,Iyear1,year,month,day,date,hour,amon)
      call ESMF_TimeSet(stopTime, yy=year, mm=month, d=date, h=hour, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_TimePrint(stopTime)

      call ESMF_TimeIntervalSet(timeStep, s_r8=dtsrc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_TimeIntervalPrint(timeStep)

!****
      RETURN
!****
!**** TERMINATE BECAUSE OF IMPROPER PICK-UP
!****
  900 write (6,*) 'Error in NAMELIST parameters'
      call stop_model('Error in NAMELIST parameters',255)
  910 write (6,*) 'Error readin I-file'
      call stop_model('Error reading I-file',255)
      END SUBROUTINE initialize_clocks

end module
#endif
