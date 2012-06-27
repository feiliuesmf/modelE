module BaseTime_mod
  implicit none
  private

  public :: BaseTime

  type BaseTime
    private
    integer*8 :: wholeSeconds
  contains
    procedure :: set
    procedure :: get
    procedure :: add
  end type BaseTime

contains

  subroutine set(this, wholeSeconds)
    class (BaseTime), intent(out) :: this
    integer*8, intent(in) :: wholeSeconds
    this%wholeSeconds = wholeSeconds
  end subroutine set

  integer*8 function get(this) result(wholeSeconds)
    class (BaseTime), intent(in) :: this
    wholeSeconds = this%wholeSeconds
  end function get

  subroutine add(this, increment)
    class (BaseTime), intent(inout) :: this
    integer*8, intent(in) :: increment
    this%wholeSeconds = this%wholeSeconds + increment
  end subroutine add

end module BaseTime_mod
