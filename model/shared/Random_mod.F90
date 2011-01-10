module RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI,IBM,Linux,DEC)
!@cont  RANDU, RINIT, RFINAL
  implicit none
  integer, save :: IX            !@var IX     random number seed

  ! Parameters used for "burning" sequences of random numbers
  integer, parameter :: A_linear = 69069
!!$$  integer, parameter :: A_linear = 65539  ! an alternate seed

  integer, parameter :: MAX_BITS = 31
  integer, parameter :: B_Half = 2**(MAX_BITS-1)
  integer, parameter :: B_linear = B_Half + (B_Half-1) ! 2147483647 = 2^31-1
  real*8, parameter :: DENOMINATOR = (1.d+0/B_HALF)/2

contains

  real function RANDU (x)
!@sum   RANDU calculates a random number based on the seed IX
    real*8 :: x ! unused

    ix = ix * A_LINEAR + 1
    if (ix < 0) ix = ix + B_LINEAR
    RANDU = ix / DENOMINATOR
  end function RANDU

  subroutine RINIT (INIT)
!@sum   RINIT sets the initial seed IX
    integer, intent(IN)  :: INIT   !@var INIT   first random no. seed
    IX=INIT
    return
  end subroutine RINIT

  subroutine RFINAL (IFINAL)
!@sum   RFINAL retrieves final seed value
    integer, intent(OUT) :: IFINAL !@var IFINAL last random no. seed
    IFINAL=IX
    return
  end subroutine RFINAL

  subroutine BURN_RANDOM(n)
!@sum  BURN_RANDOM burns a set number of random numbers. It is used to
    !                  maintain bit-wise correspondence on parallel runs.
    implicit none
    integer, intent(in) :: n
    integer :: i
    real*8 x, randss
    integer :: a, b ! linear coefficient
    integer :: nn
    if (n.eq.0) return
    a = A_linear
    nn = n
    b = 1
    do i=1,MAX_BITS
      if (mod(nn,2) == 1) ix = ix * a + b
      b=(a+1)*b
      a=a*a
      nn=nn/2
      if (nn == 0) exit
    end do
    return
  end subroutine burn_random

end module RANDOM
