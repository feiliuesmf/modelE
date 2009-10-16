! This module implements a singleton intended to manage the set
! of timers associated with an entire application.  A special
! timer "main" is created at initialization. 
module TimerList_mod
   use Timer_mod, only: Timer_type
   implicit none
   private

   public :: initialize
   public :: finalize

   public :: getNumTimers
   public :: isTimer
   public :: addTimer
   public :: getTimer

   public :: start
   public :: stop
   public :: printSummary
   public :: printParallelSummary

   public :: reset

   integer :: numTimers = 0

   integer, parameter :: MAX_TIMERS = 100
   integer, parameter :: MAX_NAME_LENGTH = 20
   integer, parameter :: NOT_FOUND =-1

   character(len=MAX_NAME_LENGTH), save :: names(MAX_TIMERS)
   type (Timer_type), target,      save :: timers(MAX_TIMERS)

   interface start
      module procedure startByName
   end interface

   interface stop
      module procedure stopByName
   end interface

   interface reset
      module procedure reset_
   end interface

   interface printSummary
      module procedure printSummaryDefault
      module procedure printSummaryList
   end interface

   integer, parameter :: r64 = selected_real_kind(14)
   
   logical, save :: test_ = .false.

contains

   ! test enables
   subroutine initialize()
      numTimers = 0
      call addTimer('main')
      call start('main')
   end subroutine initialize
   
   subroutine finalize()
      call stop('main')
   end subroutine finalize

   subroutine reset_()
      numTimers = 0
   end subroutine reset_

   logical function isTimer(name)
      character(len=*), intent(in) :: name

      integer :: i

      isTimer = getIndex(name) /= NOT_FOUND

   end function isTimer
   
   integer function getNumTimers()

      getNumTimers = numTimers

   end function getNumTimers
   
   subroutine addTimer(name)
      use Timer_mod, only: reset
      character(len=*), intent(in) :: name

      numTimers = numTimers + 1
      names(numTimers) = trim(name)
      call reset(timers(numTimers))

   end subroutine addTimer

   subroutine startByName(name)
      use Timer_mod, only: start
      character(len=*), intent(in) :: name

      call start(timers(getIndex(name)))

   end subroutine startByName

   subroutine stopByName(name)
      use Timer_mod, only: stop
      character(len=*), intent(in) :: name

      call stop(timers(getIndex(name)))

   end subroutine stopByName

   integer function getIndex(name)
      character(len=*), intent(in) :: name
      integer :: i

      do i = 1, numTimers
         if (trim(name) == trim(names(i))) then
            getIndex = i
            return
         end if
      end do

      getIndex = NOT_FOUND

   end function getIndex

   function getTimer(name)
      character(len=*), intent(in) :: name
      type (Timer_type), pointer :: getTimer

      getTimer => timers(getIndex(name))

   end function getTimer

   subroutine printSummaryDefault(report)
      use Timer_mod, only: summary
      character(len=*), pointer :: report(:)
      integer :: i,idx

      call printSummary(report, timers)

   end subroutine printSummaryDefault

   subroutine printSummaryList(report, timers)
      use Timer_mod, only: summary
      use Timer_mod, only: getInclusiveTime, getExclusiveTime
      character(len=*), pointer :: report(:)
      type (Timer_type), intent(in) :: timers(:)
      integer :: i,idx
      real (kind=r64) :: timeMain
      real (kind=r64) :: fraction

      allocate(report(4 + numTimers))
      write(report(1),"(T39,a,T77,a)")'Total Time','Trip Time'
      write(report(2),"(T8,a,T27,a,T34,a,T45,a,T56,a,T73,a,T81,a,T91,a,T101,a)") &
           & 'Timer', '%', 'seconds', 'Inclusive','(  Exclusive )',  &
           & 'Cycles','Average', 'Maximum', 'Minimum'
      write(report(3),"(120('-'))")

      timeMain = getInclusiveTime(timers(1))
      do i = 1, numTimers
         idx = 3+i
         report(idx)(1:22) = trim(names(i))
         fraction = getExclusiveTime(timers(i)) / timeMain
         call writePercentage(report(idx)(23:29), fraction)
         report(idx)(30:)  = trim(summary(timers(i)))
      end do

      write(report(4+numTimers),"(120('-'))")

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

   subroutine printParallelSummary(comm, report)
      use Timer_mod, only: summary, gather
      integer, intent(in) :: comm
      character(len=*), pointer :: report(:)
      character(len=220) :: line
      integer :: i
      integer :: rank, ier
      integer, parameter :: root = 0
      integer :: idx

      type (Timer_type) :: globalTimers(numTimers)

      do i = 1, numTimers
         globalTimers(i) = gather(timers(i), comm)
      end do

      call mpi_comm_rank(comm, rank, ier)
      if (rank == root) then
         call printSummary(report, globalTimers)
      else
         allocate(report(1))
         report(1) = ' '
      end if
   end subroutine printParallelSummary

end module TimerList_mod
