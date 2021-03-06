! This module is an intentional kludge.  Time and Calendar classes
! collaborate with each other, but the rules of F2003 do not allow
! circular module dependencies.  Thus, the classes must be declared
! in the same module.  We pretend that this is not true by treating this
! as a private module and re-exposing the classes through the modules
! Time_mod and Calendar_mod

module TimeAndCalendar_mod
  use BaseTime_mod
  implicit none
  private

  public :: Time
  public :: newTime

  public :: Calendar

  type, abstract :: Calendar
    integer :: placeholder
  contains
    procedure(Iget), deferred :: getYear
    procedure(Iget), deferred :: getDayOfYear
    procedure(Iget), deferred :: getMonth
    procedure(IgetAbbrev), deferred :: getAbbreviation
    procedure(Iget), deferred :: getDate
    procedure(Iconvert), deferred :: convertToTime
    procedure :: getTimeOfDay ! seconds since 0h:0m:0.0s
  end type Calendar

  type, extends(BaseTime) :: Time
    class (Calendar), pointer :: calendar => null()
  contains
    procedure :: setByDate
    procedure :: getYear
    procedure :: getMonth
    procedure :: getAbbreviation
    procedure :: getDate
    procedure :: getDayOfyear
    procedure :: getHour
  end type Time

  abstract interface

    integer function Iget(this, t)
      import Time
      import Calendar
      class (Calendar), intent(in) :: this
      class (Time), intent(in) :: t
    end function Iget

    function IgetAbbrev(this, t) result(abbrev)
      use Month_mod, only: LEN_MONTH_ABBREVIATION
      import Time
      import Calendar
      character(len=LEN_MONTH_ABBREVIATION) :: abbrev
      class (Calendar), intent(in) :: this
      class (Time), intent(in) :: t
    end function IgetAbbrev

    function Iconvert(this, year, month, date, hour) result(t)
      use BaseTime_mod
      import Calendar
      class (Calendar), intent(in) :: this
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: date
      integer, intent(in) :: hour
      type (BaseTime) :: t
    end function Iconvert

    function IgetTime(this, t) result(timeOfDay)
      use BaseTime_mod
      import Time
      import Calendar
      class (Calendar), intent(in) :: this
      class (Time), intent(in) :: t
      type (BaseTime) :: timeOfDay
    end function IgetTime

  end interface

contains

  function newTime(aCalendar) result(t)
    class (Calendar), pointer :: aCalendar
    type (Time) :: t
    t%calendar => aCalendar
  end function newTime

  subroutine setByDate(this, year, month, date, hour)
    class (Time), intent(inout) :: this
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: date
    integer, intent(in) :: hour
    
    this%BaseTime = this%calendar%convertToTime(year, month, date, hour)
    
  end subroutine setByDate

  integer function getYear(this) result(year)
    class (Time), intent(in) :: this
    year = this%calendar%getYear(this)
  end function getYear

  integer function getMonth(this) result(month)
    class (Time), intent(in) :: this
    month = this%calendar%getMonth(this)
  end function getMonth

  function getAbbreviation(this) result(abbrev)
    use Month_mod, only: LEN_MONTH_ABBREVIATION
    character(len=LEN_MONTH_ABBREVIATION) :: abbrev
    class (Time), intent(in) :: this
    abbrev = this%calendar%getAbbreviation(this)
  end function getAbbreviation

  integer function getDate(this) result(date)
    class (Time), intent(in) :: this
    date = this%calendar%getDate(this)
  end function getDate

  integer function getDayOfYear(this) result(dayOfYear)
    class (Time), intent(in) :: this
    dayOfYear = this%calendar%getDayOfYear(this)
  end function getDayOfYear

  integer function getHour(this) result(hour)
    use TimeConstants_mod, only: INT_SECONDS_PER_HOUR
    class (Time), intent(in) :: this
    type (BaseTime) :: t

    t = this%calendar%getTimeOfDay(this)
    hour = t%get() / INT_SECONDS_PER_HOUR
  end function getHour

  function getTimeOfDay(this, t) result(timeOfDay)
    use BaseTime_mod
    class (Calendar), intent(in) :: this
    class (Time), intent(in) :: t
    type (BaseTime) :: timeOfDay

    type (Time) :: timeAtBeginningOfDAy
    integer :: year, month, date

    year = t%getYear()
    month = t%getMonth()
    date = t%getDate()

    timeAtBeginningOfDay%calendar => t%calendar
    call timeAtBeginningOfDay%setByDate(year, month, date, hour=0)
    call timeOfDay%set(t%get() - timeAtBeginningOfDay%get())

  end function GetTimeOfDay

end module TimeAndCalendar_mod
