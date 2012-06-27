module Time_mod
  use TimeAndCalendar_mod, only: Time
  use TimeAndCalendar_mod, only: newTime
  implicit none
  private
  public :: Time
  public :: newTime
end module Time_mod
