! This module implements a singleton intended to manage the set
! of timers associated with an entire application.  A special
! timer "main" is created at initialization. 
module TimerList_mod
   use Timer_mod, only: Timer_type
   implicit none
   private

   public :: TimerList_type
   public :: initialize
   public :: finalize

   public :: getNumTimers
   public :: isTimer
   public :: addTimer
   public :: getTimer
   public :: getName

   public :: start
   public :: stop
   public :: printSummary
#ifdef USE_MPI
   public :: printParallelSummary
#endif

   public :: reset
   public :: getDefaultList

   integer, parameter :: MAX_NAME_LENGTH = 20

   type NamedTimer_type
      private
      type (Timer_type) :: timer
      character(len=MAX_NAME_LENGTH) :: name
   end type NamedTimer_type

   type TimerList_type
      private
      type (NamedTimer_type), pointer :: list(:) => null()
   end type TimerList_type

   type (TimerList_type), target, save :: defaultList ! singleton

   integer, parameter :: MAX_TIMERS = 100
   integer, parameter :: NOT_FOUND =-1

   interface initialize
      module procedure initialize_
      module procedure initializeDefault
   end interface

   interface finalize
      module procedure finalize_
      module procedure finalizeDefault
   end interface

   interface start
      module procedure startByName
      module procedure startByNameAtTime
      module procedure startByNameDefault
      module procedure startByNameAtTimeDefault 
   end interface

   interface stop
      module procedure stopByName
      module procedure stopByNameAtTime
      module procedure stopByNameDefault
      module procedure stopByNameAtTimeDefault
   end interface

   interface getName
      module procedure getNameByIndex
   end interface

   interface getTimer
      module procedure getTimer_
      module procedure getTimerDefault
      module procedure getTimerByIndex
   end interface

   interface isTimer
      module procedure isTimer_
      module procedure isTimerDefault
   end interface

   interface getNumTimers
      module procedure getNumTimers_
      module procedure getNumTimersDefault
   end interface

   interface addTimer
      module procedure addTimer_
      module procedure addTimerDefault
   end interface

   interface reset
      module procedure reset_
      module procedure resetDefault
      module procedure resetTimer
   end interface

   interface printSummary
      module procedure printSummaryDefault
      module procedure printSummaryList
   end interface

   integer, parameter :: r64 = selected_real_kind(14)
   
   logical, save :: test_ = .false.

contains

   function getDefaultList()
      type(TimerList_type), pointer :: getDefaultList
      getDefaultList => defaultList
   end function getDefaultList

   subroutine initialize_(this, mockTime)
      type (TimerList_type), intent(inout) :: this
      real(kind=r64), optional, intent(in) :: mockTime
      real(kind=r64) :: mockTime_

      if (present(mockTime)) mockTime_ = mockTime

      if (associated(this%list)) deallocate(this%list)
      allocate(this%list(0))

      call addTimer(this, 'main')
      if (present(mockTime)) then
         call start(this, 'main', time = mockTime_)
      else
         call start(this, 'main')
      end if
   end subroutine initialize_
      
   ! test enables
   subroutine initializeDefault(mockTime)
      real(kind=r64), optional, intent(in) :: mockTime
      call initialize(defaultList, mockTime)
   end subroutine initializeDefault
   
   subroutine finalize_(this, mockTime)
      type (TimerList_type), target, intent(inout) :: this
      real(kind=r64), optional, intent(in) :: mockTime

      real(kind=r64) :: mockTime_

      if (present(mockTime)) mockTime_ = mockTime

      if (present(mockTime)) then
         call stop(this, 'main', time = mockTime_)
      else
         call stop(this, 'main')
      end if

      call checkTimerConsistenncy()

   contains

      subroutine checkTimerConsistenncy()
#ifdef USE_PFUNIT
         use pFUnit
#endif
         use Timer_mod, only: isActive
         type (NamedTimer_type), pointer :: namedTimer
         integer :: i
         character(len=50) :: message
         do i = 1, getNumTimers(this)
            namedTimer => this%list(i)
            if (isActive(namedTimer%timer)) then
               message = 'Unbalanced start/stop for timer <'//trim(namedTimer%name)//'>.'
#ifdef USE_PFUNIT
               call throw(Exception(message))
#else
               write(*,*) message
#endif
               exit
            end if
         end do
      end subroutine checkTimerConsistenncy
       
   end subroutine finalize_

   subroutine finalizeDefault(mockTime)
      real(kind=r64), optional, intent(in) :: mockTime
      call finalize(defaultList, mockTime)
   end subroutine finalizeDefault

   subroutine reset_(this)
      type (TimerList_type), intent(inOut) :: this
      if (associated(this%list)) deallocate(this%list)
      allocate(this%list(0))
   end subroutine reset_

   subroutine resetDefault()
      call reset(defaultList)
   end subroutine resetDefault

   subroutine resetTimer(name)
      use Timer_mod, only: reset
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: timer
      timer => getTimer(name)
      if (associated(timer)) call reset(timer)
   end subroutine resetTimer

   logical function isTimer_(this, name)
      type (TimerList_type), intent(in) :: this
      character(len=*), intent(in) :: name

      isTimer_ = getIndex(this, name) /= NOT_FOUND

   end function isTimer_
   
   logical function isTimerDefault(name)
      character(len=*), intent(in) :: name
      isTimerDefault = isTimer(defaultList, name)
   end function isTimerDefault
   
   integer function getNumTimers_(this)
      type (TimerList_type), intent(in) :: this
      getNumTimers_ = size(this%list)
   end function getNumTimers_

   integer function getNumTimersDefault()
      getNumTimersDefault = getNumTimers(defaultList)
   end function getNumTimersDefault
   
   subroutine addTimer_(this, name)
      use Timer_mod, only: reset
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      type (NamedTimer_type), allocatable :: tmpList(:)

      integer :: n

      n = size(this%list)

      if (n > 0) then
         allocate(tmpList(n))
         tmpList = this%list
      end if
      if (associated(this%list)) deallocate(this%list)
      allocate(this%list(n+1))
      if (n > 0) then
         this%list(:n) = tmpList
         deallocate(tmpList)
      end if
      this%list(n+1)%name = name
      call reset(this%list(n+1)%timer)

   end subroutine addTimer_

   subroutine addTimerDefault(name)
      character(len=*), intent(in) :: name
      call addTimer(defaultList, name)
   end subroutine addTimerDefault

   subroutine startByName(this, name)
      use Timer_mod, only: start
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      integer :: index

      index = getIndex(this, name)
      if (index == NOT_FOUND) then
         call addTimer(this, name)
         index = getNumTimers(this)
      end if

      call start(this%list(index)%timer)

   end subroutine startByName

   subroutine startByNameDefault(name)
      character(len=*), intent(in) :: name
      call startByName(defaultList, name)
   end subroutine startByNameDefault

   subroutine startByNameAtTime(this, name, time)
      use Timer_mod, only: start
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time
      call start(this%list(getIndex(this, name))%timer, time)
   end subroutine startByNameAtTime

   subroutine startByNameAtTimeDefault(name, time)
      use Timer_mod, only: start
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time
      call startByNameAtTime(defaultList, name, time)
   end subroutine startByNameAtTimeDefault

   subroutine stopByName(this, name)
#ifdef USE_PFUNIT
      use pFUnit
#endif
      use Timer_mod, only: stop
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name

      character(len=70) :: message
      integer :: index

      index = getIndex(this, name)
      if (index /= NOT_FOUND) then
         call stop(this%list(index)%timer)
      else
         message = 'Timer <'//trim(name)//'> has not been declared prior to use.'
#ifdef USE_PFUNIT
         call throw(Exception(trim(message)))
#else
         write(*,*) trim(message)
#endif
      end if

   end subroutine stopByName

   subroutine stopByNameDefault(name)
      character(len=*), intent(in) :: name
      call stopByName(defaultList, name)
   end subroutine stopByNameDefault

   subroutine stopByNameAtTime(this, name, time)
      use Timer_mod, only: stop
      type (TimerList_type), intent(inOut) :: this
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time

      integer :: index
      index = getIndex(this, name)
      call stop(this%list(index)%timer, time)

   end subroutine stopByNameAtTime

   subroutine stopByNameAtTimeDefault(name, time)
      use Timer_mod, only: stop
      character(len=*), intent(in) :: name
      real(kind=r64), intent(in) :: time
      call stopByNameAtTime(defaultList, name, time)
   end subroutine stopByNameAtTimeDefault

   integer function getIndex(this, name)
      type (TimerList_type), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: i

      do i = 1, getNumTimers(this)
         if (trim(name) == trim(this%list(i)%name)) then
            getIndex = i
            return
         end if
      end do

      getIndex = NOT_FOUND

   end function getIndex

   function getTimer_(this, name) result(timer)
      type (TimerList_type), target, intent(in) :: this
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: timer

      timer => this%list(getIndex(this, name))%timer

   end function getTimer_

   function getTimerDefault(name) result(timer)
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: timer
      timer => getTimer(defaultList, name)
   end function getTimerDefault

   function getTimerByIndex(this, index) result(timer)
      type (TimerList_type), target, intent(in) :: this
      integer, intent(in) :: index
      type (Timer_type), pointer :: timer

      timer => this%list(index)%timer

   end function getTimerByIndex

   function getNameByIndex(this, index) result(name)
      type (TimerList_type), target, intent(in) :: this
      integer, intent(in) :: index
      character(len=MAX_NAME_LENGTH) :: name

      name = this%list(index)%name

   end function getNameByIndex

   subroutine printSummaryDefault(report, scale, units)
      use Timer_mod, only: summary
      character(len=*), pointer :: report(:)
      real(kind=r64), optional, intent(in) :: scale
      character(len=*), optional, intent(in) :: units

      call printSummary(report, defaultList, scale, units)

   end subroutine printSummaryDefault

   subroutine printSummaryList(report, timerList, scale, units)
      use Timer_mod, only: summary
      use Timer_mod, only: getInclusiveTime, getExclusiveTime
      character(len=*), pointer :: report(:)
      type (TimerList_type), target, intent(in) :: timerList
      real(kind=r64), optional, intent(in) :: scale
      character(len=*), optional, intent(in) :: units

      character(len=7) :: units_
      integer :: i,idx
      real (kind=r64) :: timeMain
      real (kind=r64) :: fraction
      integer :: numTimers
      type (NamedTimer_type), pointer :: namedTimer

      if (present(units)) then
         units_ = units(1:min(7,len(units)))
      else
         units_ = 'seconds'
      end if

      numTimers = getNumTimers(timerList)
      allocate(report(4 + numTimers))
      write(report(1),"(T46,a,T88,a)")'Total Time (hhh:mm:ss.hh)','Trip Time (seconds)'
      write(report(2),"(T8,a,T27,a,T34,a,T47,a,T57,a,T75,a,T85,a,T99,a,T113,a)") &
           & 'Timer', '%', units_, 'Inclusive','(  Exclusive )',  &
           & 'Cycles','Average', 'Maximum', 'Minimum'
      write(report(3),"(122('-'))")

      timeMain = getInclusiveTime(timerList%list(1)%timer)
      do i = 1, numTimers
         namedTimer => timerList%list(i)
         idx = 3+i
         report(idx)(1:22) = namedTimer%name
         fraction = getExclusiveTime(namedTimer%timer) / timeMain
         call writePercentage(report(idx)(23:29), fraction)
         report(idx)(30:)  = trim(summary(namedTimer%timer, scale))
      end do

      write(report(4+numTimers),"(122('-'))")

   contains

      subroutine writePercentage(string, fraction)
         character(len=*), intent(inOut) :: string
         real (kind=r64), intent(in) :: fraction
         
         integer :: ticks, percentage
         
         ticks = nint(fraction * 100**2)
         percentage = ticks/100
         write(string,'(1x,i3.1,".",i2.2)') percentage, ticks - percentage*100
         
   end subroutine writePercentage

   end subroutine printSummaryList

#ifdef USE_MPI
   subroutine printParallelSummary(comm, report)
      use Timer_mod, only: summary, gather
      integer, intent(in) :: comm
      character(len=*), pointer :: report(:)
      character(len=220) :: line
      integer :: i
      integer :: rank, ier
      integer, parameter :: root = 0
      integer :: idx

      type (TimerList_type) :: globalTimerList
      integer :: numTimers

      numTimers = getNumTimers(defaultList)

      allocate(globalTimerList%list(numTimers))
      do i = 1, numTimers
         globalTimerList%list(i)%name = defaultList%list(i)%name
         globalTimerList%list(i)%timer = gather(defaultList%list(i)%timer, comm)
      end do

      call mpi_comm_rank(comm, rank, ier)
      if (rank == root) then
         call printSummary(report, globalTimerList)
      else
         allocate(report(1))
         report(1) = ' '
      end if

      deallocate(globalTimerList%list)
   end subroutine printParallelSummary
#endif

end module TimerList_mod
