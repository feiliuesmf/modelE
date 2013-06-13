module ModelClock_mod
  use Time_mod
  implicit none
  private

  public :: ModelClock
  public :: newModelClock

  type :: ModelClock
    type (Time) :: currentTime
    type (Time) :: startTime

    ! modelE legacy representation
    integer :: timeStep
    integer :: stepsPerDay

  contains
    procedure :: getTimeInSecondsFromDate
    procedure :: getAbsoluteTimeInSeconds
    procedure :: getTimeStep
    procedure :: isBeginningOfDay
    procedure :: nextTick
    procedure :: getDate
    procedure :: year
    procedure :: month
    procedure :: date
    procedure :: dayOfYear
    procedure :: hour
    procedure :: abbrev ! month abbreviation
  end type ModelClock

contains

  ! constructor
  function newModelClock(startTime, timeStep, stepsPerDay) result(clock)
    type (ModelClock) :: clock
    type (Time), intent(in) :: startTime
    integer, intent(in) :: timeStep
    integer, intent(in) :: stepsPerDay

    clock%timeStep = timeStep
    clock%stepsPerDay = stepsPerDay

    clock%startTime= startTime
    clock%currentTime = startTime
  end function newModelClock

  subroutine nextTick(this)
    use TimeConstants_mod, only: INT_SECONDS_PER_DAY
    class (ModelClock), intent(inout) :: this
    integer*8 :: stepSize

    this%timeStep = this%timeStep + 1
    stepSize = INT_SECONDS_PER_DAY / this%stepsPerDay
    call this%currentTime%add(stepSize)
  end subroutine nextTick

  integer function getTimeStep(this)
    class (ModelClock), intent(in) :: this
    getTimeStep = this%timeStep
  end function getTimeStep

  logical function isBeginningOfDay(this)
    class (ModelClock), intent(in) :: this
    
    isBeginningOfDay = mod(this%timeStep, this%stepsPerDay) == 0
  end function isBeginningOfDay

  function getAbsoluteTimeInSeconds(this) result (secs)
    integer*8 :: secs
    class (ModelClock), intent(inout) :: this
    secs = this%currentTime%get()
  end function getAbsoluteTimeInSeconds

  function getTimeInSecondsFromDate(this, year, month, date, hour) result (secs)
    use Calendar_mod, only: Calendar
    use BaseTime_mod
    integer*8 :: secs
    class (ModelClock), intent(inout) :: this
    integer, intent(in) :: year, month, date, hour
    type (Time) :: aTime
    class (Calendar), pointer :: pCalendar

    pCalendar => this%currentTime%calendar
    aTime = newTime(pCalendar)
    call aTime%setByDate(year, month, date, hour)
    secs = this%currentTime%get() - aTime%get()

  end function getTimeInSecondsFromDate


  subroutine getDate(this, year, month, dayOfYear, date, hour, amn)
!@sum  getDate gets Calendar info from internal timing info
!@auth Gavin Schmidt (updated by Tom CLune)
    use TimeConstants_mod, only: INT_SECONDS_PER_HOUR
    use JulianCalendar_mod, only: LAST_JULIAN_DAY_IN_MONTH
    use JulianCalendar_mod, only: JULIAN_MONTHS
    use Month_mod, only: LEN_MONTH_ABBREVIATION, Month_type

    class (ModelClock), intent(in) :: this
    integer, optional, intent(out) :: year
    integer, optional, intent(out) :: month
    integer, optional, intent(out) :: dayOfYear
    integer, optional, intent(out) :: date
    integer, optional, intent(out) :: hour
    character(len=LEN_MONTH_ABBREVIATION), optional, intent(out) :: amn
    integer :: mnth

    if (present(year)) year = this%currentTime%getYear()
    if (present(dayOfYear)) dayOfYear = this%currentTime%getDayOfYear()
    if (present(month)) month = this%currentTime%getMonth()
    if (present(amn)) then
      mnth = this%currentTime%getMonth()
      amn = JULIAN_MONTHS(mnth)%abbreviation
    end if
    if (present(date)) date = this%currentTime%getDate()
    if (present(hour)) hour = this%currentTime%getHour()

    return
  end subroutine getDate

  integer function year(this)
    class (ModelClock), intent(in) :: this
    year = this%currentTime%getYear()
  end function year

  integer function month(this)
    class (ModelClock), intent(in) :: this
    month = this%currentTime%getMonth()
  end function month

  integer function date(this)
    class (ModelClock), intent(in) :: this
    date = this%currentTime%getDate()
  end function date

  integer function dayOfYear(this)
    class (ModelClock), intent(in) :: this
    dayOfYear = this%currentTime%getDayOfYear()
  end function dayOfYear

  integer function hour(this)
    class (ModelClock), intent(in) :: this
    hour = this%currentTime%getHour()
  end function hour

  function abbrev(this)
    use Month_mod, only: LEN_MONTH_ABBREVIATION
    character(len=LEN_MONTH_ABBREVIATION) abbrev
    class (ModelClock), intent(in) :: this

    abbrev = this%currentTime%getAbbreviation()
  end function abbrev

end module ModelClock_mod
