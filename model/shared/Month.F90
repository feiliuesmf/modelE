module Month_mod
  implicit none
  private

  public :: Month_type
  integer, parameter, public :: LEN_MONTH_ABBREVIATION = 4
  integer, parameter, public :: LEN_MONTH_NAME = 20

  type Month_type
    character(len=LEN_MONTH_ABBREVIATION) :: abbreviation
    character(len=LEN_MONTH_NAME) :: longName
    integer :: numDays
    integer :: lastJulianDay
    integer :: middleJulianDay
  end type Month_type

end module Month_mod
