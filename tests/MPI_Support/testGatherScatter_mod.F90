module testGatherScatter_mod
  use pFUnit
  use GatherScatter_mod
  implicit none
  private

  public :: testIntegerGather
  public :: testIntegerScatter
  public :: testReal8Gather
  public :: testReal8Scatter

contains

  subroutine testIntegerGather()
  integer :: expected, found
  expected = 1
  found = 0
  call assertEqual(expected, found)
  end subroutine testIntegerGather

  subroutine testIntegerScatter()
  integer :: expected, found
  expected = 1
  found = 0
  call assertEqual(expected, found)
  end subroutine testIntegerScatter

  subroutine testReal8Gather()
  real(8) :: expected, found
  real(8), parameter :: eps=1.0e-8
  expected = 1.0
  found = 0.0 
  call assertEqual(expected, found, tolerance=eps)
  end subroutine testReal8Gather

  subroutine testReal8Scatter()
  real(8) :: expected, found
  real(8), parameter :: eps=1.0e-8
  expected = 1.0
  found = 0.0 
  call assertEqual(expected, found, tolerance=eps)
  end subroutine testReal8Scatter

end module testGatherScatter_mod
