! This module implements a minimal capability for profiling fortran
! applications.  It is meant to be directly included in other source
! trees and thereby be very easily used.  Products which are far
! superior in their capabilities are available, but generally require
! a larger investment to integrate with a user's application.
! Perhaps eventually, this module could be used as an interface to
! other libraries so that the user need not modify their interfaces.

module Timer_mod
   implicit none
   private

   public :: Timer_type ! derived type

   public :: start
   public :: stop
   public :: reset
   public :: getAverageTripTime

   ! Accessor methods
   public :: getNumTrips   
   public :: isActive
   public :: getInclusiveTime
   public :: getExclusiveTime
   public :: getMaximumTripTime
   public :: getMinimumTripTime

   ! public paramteers
   public :: TIMER_SUMMARY_LENGTH

   ! private
   public :: addTrip ! still public for testing purposes
   private :: addTime ! for testing purposes
   private :: setActive
   private :: setInactive

   public :: r64

   integer, parameter :: TIMER_SUMMARY_LENGTH =  200
   integer, parameter :: r64 = selected_real_kind(14)

   public :: summary
#ifdef USE_MPI
   public :: gather
   public :: parallelSummary
#endif

   type Timer_type
      private
      real (kind=r64) :: inclusiveTime  = 0
      real (kind=r64) :: exclusiveTime  = 0
      real (kind=r64) :: maximumTripTime  = 0
      real (kind=r64) :: minimumTripTime  = huge(1._r64)
      integer         :: numTrips     = 0
      real (kind=r64) :: startTime    = 0.
      real (kind=r64) :: startExclusiveTime = 0
      logical         :: isActive     = .false.
   end type Timer_type

   interface start
      module procedure start_
      module procedure startAtTime
   end interface

   interface stop
      module procedure stop_
      module procedure stopAtTime
   end interface

   interface reset
      module procedure reset_
   end interface

   interface summary
      module procedure summary_
   end interface

   ! private shared variable for computing exclusive
   ! time
   real (kind=r64), save :: globalExclusiveTime = 0

contains

   integer function getNumTrips(this)
      type (Timer_type), intent(in) :: this
      getNumTrips = this%numTrips
   end function getNumTrips

   real(kind=r64) function getWTime() result(time)
#ifdef USE_MPI
      include 'mpif.h'
      time = mpi_Wtime()
#else      
      integer :: counter, rate
      call system_clock(counter, rate)
      time = real(counter,kind=r64) / rate
#endif
   end function getWTime

   subroutine start_(this)
      type (Timer_type), intent(inout) :: this
      call startAtTime(this, getWTime())
   end subroutine start_

   subroutine startAtTime(this, time)
      type (Timer_type), intent(inout) :: this
      real(kind=r64), intent(in) :: time

      this%isActive = .true.
      this%startTime = time
      this%startExclusiveTime = globalExclusiveTime
      call addTrip(this)

   end subroutine startAtTime

   subroutine stop_(this)
      type (Timer_type), intent(inout) :: this

      call stopAtTime(this, getWTime())

   end subroutine stop_

   subroutine stopAtTime(this, time)
      type (Timer_type), intent(inout) :: this
      real(kind=r64), intent(in) :: time

      real(kind=r64) :: dtInclusive, dtExclusive

      dtInclusive = (time - this%startTime)
      dtExclusive = dtInclusive - (globalExclusiveTime - this%startExclusiveTime)
      call addTime_(this, dtInclusive, dtExclusive)
      this%isActive = .false.

   end subroutine stopAtTime

   function clockTick()
      real (kind=r64) :: clockTick
#ifdef USE_MPI
      include 'mpif.h'
      clockTick = mpi_Wtick()
      write(*,*)'MPI: tick = ', clockTick
#else      
      integer :: clockRate
      call system_clock(count_rate = clockRate)
      clockTick = 1 / real(clockRate, kind=r64)
      write(*,*)'serial: tick = ', clockTick
#endif
   end function clockTick

   subroutine setActive(this)
      type (Timer_type), intent(inout) :: this
      this%isActive = .true.
   end subroutine setActive

   subroutine setInactive(this)
      type (Timer_type), intent(inout) :: this
      this%isActive = .false.
   end subroutine setInactive

   subroutine reset_(this)
      type (Timer_type), intent(inout) :: this
      call setInactive(this)
      this%numTrips = 0
      this%inclusiveTime = 0
      this%exclusiveTime = 0
      this%maximumTripTime = 0
      this%minimumTripTime = huge(1._r64)
   end subroutine reset_

   subroutine addTime(this, dtInclusive, dtExclusive)
      type (Timer_type), intent(inout) :: this
      real (kind=r64),   intent(in)    :: dtInclusive
      real (kind=r64), optional, intent(in) :: dtExclusive

      real (kind=r64) :: dtExclusive_

      dtExclusive_ = dtInclusive
      if (present(dtExclusive)) dtExclusive_ = dtExclusive
      call addTime_(this, dtInclusive, dtExclusive_)
      
   end subroutine addTime

   subroutine addTime_(this, dtInclusive, dtExclusive)
      type (Timer_type), intent(inout) :: this
      real (kind=r64),   intent(in)    :: dtInclusive
      real (kind=r64),   intent(in)    :: dtExclusive

      this%inclusiveTime = this%inclusiveTime + dtInclusive
      this%exclusiveTime = this%exclusiveTime + dtExclusive

      this%maximumTripTime = max(this%maximumTripTime, dtExclusive)
      this%minimumTripTime = min(this%minimumTripTime, dtExclusive)

      globalExclusiveTime = globalExclusiveTime + dtExclusive

   end subroutine addTime_

   function getInclusiveTime(this) result (inclusiveTime)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: inclusiveTime

      inclusiveTime = this%inclusiveTime

   end function getInclusiveTime

   function getExclusiveTime(this) result (exclusiveTime)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: exclusiveTime

      exclusiveTime = this%exclusiveTime

   end function getExclusiveTime

   logical function isActive(this)
      type (Timer_type), intent(in) :: this

      isActive = this%isActive

   end function isActive

   function getMaximumTripTime(this) result(time)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: time

      time = this%maximumTripTime

   end function getMaximumTripTime

   function getMinimumTripTime(this) result(time)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: time

      if (this%numTrips > 0) then
         time = this%minimumTripTime
      else
         time = 0
      end if

   end function getMinimumTripTime

   function getAverageTripTime(this) result(time)
      type (Timer_type), intent(in) :: this
      real (kind=r64) :: time

      if (this%numTrips == 0) then
         time = 0
      else
         time = this%exclusiveTime / this%numTrips
      end if

   end function getAverageTripTime

   subroutine addTrip(this)
      type (Timer_type), intent(inout) :: this
      this%numTrips = this%numTrips + 1
   end subroutine addTrip

   function summary_(this, scale) result(report)
      use TimeFormatUtilities_mod, only: formatHMS, formatSeconds
      type (Timer_type), intent(in) :: this
      real(kind=r64), optional, intent(in) :: scale
      character(len=TIMER_SUMMARY_LENGTH) :: report
      character(len=15) :: hmsInclusive
      character(len=15) :: hmsExclusive
      character(len=15) :: seconds
      character(len=60) :: tripStats

      real(kind=r64) :: scale_
      if (present(scale)) then
         scale_ = scale
      else
         scale_ = 1 ! seconds
      end if

      seconds = formatSeconds(scale_ * getInclusiveTime(this), &
           & decimalsBeforePoint = 8)
      hmsInclusive = formatHMS(getInclusiveTime(this))
      hmsExclusive = formatHMS(getExclusiveTime(this))
      tripStats = countMinMaxAvg(getNumTrips(this), getMinimumTripTime(this), getMaximumTripTime(this), getAverageTripTime(this))
      report = trim(seconds) // ' | ' // trim(hmsInclusive) // ' (' // trim(hmsExclusive) // ')' // '    '  // trim(tripStats)

   end function summary_

   function countMinMaxAvg(count, minTime, maxTime, avgTime) result (string)
      use TimeFormatUtilities_mod, only: formatSeconds
      integer, intent(in) :: count
      real (kind=r64), intent(in) :: minTime, maxTime, avgTime

      character(len=100) :: string
      character(len=20) :: avgStr, maxStr, minStr

      avgStr = formatSeconds(avgTime, decimalsAfterPoint=6)
      maxStr = formatSeconds(maxTime, decimalsAfterPoint=6)
      minStr = formatSeconds(minTime, decimalsAfterPoint=6)
      
      write(string,'(i4,3(2x,a12))') count, &
           & avgStr, maxStr, minStr

   end function countMinMaxAvg

#ifdef USE_MPI
   function parallelSummary(this, comm) result(report)
      use TimeFormatUtilities_mod, only: formatHMS, formatSeconds
      type (Timer_type), intent(in) :: this
      integer, intent(in) :: comm

      character(len=TIMER_SUMMARY_LENGTH) :: report
      character(len=15) :: hmsInclusive, hmsExclusive
      character(len=60) :: tripStats
      character(len=150) :: seconds
      integer, parameter :: root = 0

      real (kind=r64) :: inclusiveTime, exclusiveTime, maxTime, minTime, avgTime

      integer :: globalCount
      integer :: rank, npes, ier

      include 'mpif.h'

      call MPI_Comm_size(comm, npes, ier)
      call MPI_Comm_rank(comm, rank, ier)

      call MPI_Reduce(getInclusiveTime(this), inclusiveTime, 1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, root, comm, ier)
      call MPI_Reduce(getExclusiveTime(this), exclusiveTime, 1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, root, comm, ier)
      call MPI_Reduce(getMaximumTripTime(this), maxTime, 1, MPI_DOUBLE_PRECISION, &
           & MPI_MAX, root, comm, ier)
      call MPI_Reduce(getMinimumTripTime(this), minTime, 1, MPI_DOUBLE_PRECISION, &
           & MPI_MIN, root, comm, ier)
      call MPI_Reduce(getNumTrips(this), globalCount, 1, MPI_INTEGER, &
           & MPI_SUM, root, comm, ier)

      if (rank == root) then
         seconds = formatSeconds(inclusiveTime)
         hmsInclusive = formatHMS(inclusiveTime)
         hmsExclusive = formatHMS(exclusiveTime)
         avgTime = exclusiveTime/npes
         tripStats = countMinMaxAvg(globalCount, minTime, maxTime, avgTime)

         report = trim(seconds) // ' | ' // trim(hmsInclusive) // ' (' // trim(hmsExclusive) // ')' // '    '  // trim(tripStats)
      else
         report = ' '
      end if

   end function parallelSummary

   function gather(this, comm) result(globalTimer)
      type (Timer_type), intent(in) :: this
      integer, intent(in) :: comm
      type (Timer_type) :: globalTimer

      integer :: npes, rank, ier
      integer, parameter :: root = 0

      include 'mpif.h'

      call MPI_Comm_size(comm, npes, ier)
      call MPI_Comm_rank(comm, rank, ier)

      call MPI_Allreduce(getNumTrips(this),        globalTimer%numTrips,        1, MPI_INTEGER, &
           & MPI_SUM, comm, ier)
      call MPI_Allreduce(getInclusiveTime(this),   globalTimer%inclusiveTime,   1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, comm, ier)
      call MPI_Allreduce(getexclusiveTime(this),   globalTimer%exclusiveTime,   1, MPI_DOUBLE_PRECISION, &
           & MPI_SUM, comm, ier)
      call MPI_Allreduce(getMaximumTripTime(this), globalTimer%maximumTripTime, 1, MPI_DOUBLE_PRECISION, &
           & MPI_MAX, comm, ier)
      call MPI_Allreduce(getMinimumTripTime(this), globalTimer%minimumTripTime, 1, MPI_DOUBLE_PRECISION, &
           & MPI_MIN, comm, ier)

   end function gather
#endif

end module Timer_mod
