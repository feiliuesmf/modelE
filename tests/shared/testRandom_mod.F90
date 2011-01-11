module testRandom_mod
  use pFUnit
  use RANDOM
  implicit none
  private

  public :: testSameAsIfort
  public :: testInRange
  public :: testBurnRandom

contains

  subroutine testSameAsIfort()
#ifdef COMPILER_Intel8
    real*4 :: ran ! intel internal
    integer :: ix

    ix = 1;
    call RINIT(ix)

    call assertEqual(ran(IX), RANDU(1.d+0))
    call assertEqual(ran(IX), RANDU(1.d+0))
    call assertEqual(ran(IX), RANDU(1.d+0))
    call assertEqual(ran(IX), RANDU(1.d+0))

#endif

  end subroutine testSameAsIfort

  subroutine testInRange()
    real*8 :: x

    call RINIT(1)

    x = RANDU(1.d+0)
    call assertTrue(0 <= x .and. x < 1)
    x = RANDU(1.d+0)
    call assertTrue(0 <= x .and. x < 1)
    x = RANDU(1.d+0)
    call assertTrue(0 <= x .and. x < 1)

  end subroutine testInRange

  subroutine testBurnRandom()

    call checkSkip(1)
    call checkSkip(5)

  contains
    
    subroutine checkSkip(n)
      integer, intent(in) :: n

      real*8 :: x, xSkip
      integer :: i

      call RINIT(1)

      do i = 1, n+1
        x = RANDU(1.d+0)
      end do

      call RINIT(1)
      call burn_random(n)
      xSkip = RANDU(1.d+0)

      call assertEqual(x, xSkip)
    end subroutine checkSkip

  end subroutine testBurnRandom

end module testRandom_mod

