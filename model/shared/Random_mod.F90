module RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI,IBM,Linux,DEC)
!@cont  RANDU, RINIT, RFINAL
  implicit none
  integer, save :: IX            !@var IX     random number seed

  ! Parameters used for "burning" sequences of random numbers
#if defined( MACHINE_DEC ) || ( MACHINE_Linux)
  integer, parameter :: A_linear = 69069
#elif defined(MACHINE_SGI) \
  || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) ) \
       || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
  ! don't know how RAN() is implemented, so burn
#else
  integer, parameter :: A_linear = 65539
#endif
  integer, parameter :: MAX_BITS = 31
  integer, parameter :: B_Half = 2**(MAX_BITS-1)
  integer, parameter :: B_linear = B_Half + (B_Half-1) ! 2147483647 = 2^31-1

contains

#if defined(MACHINE_SGI) \
  || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) && ! defined(COMPILER_NAG) ) \
       || defined(MACHINE_DEC) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_Intel8) ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
  function RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
!@calls RAN
    real*8 X                       !@var X      dummy variable
    real*4 RAN                     !@fun RAN    SGI intrinsic func.
    real*8 :: RANDU                !@var RANDU  random number
    RANDU=RAN(IX)
    return
  end function RANDU
#elif defined( MACHINE_IBM ) \
  || ( defined(COMPILER_NAG) ) \
  || ( defined(MACHINE_Linux) && defined(COMPILER_G95) ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_G95) ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
  function RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
    real*8 X                       !@var X      dummy variable
    real*8 :: RANDU                !@var RANDU  random number
    integer :: IY                  !@var IY     dummy integer
10  IY=IX*A_linear
    select case (IY)
    case (:-1)
      IY=(IY+B_linear)+1
    case (0)
      IX=1
      GO TO 10
    end select
    IX=IY
    RANDU=dble(IY)*.465661287308D-9
    return
  end function RANDU
#else
  none of supported architectures was specified.
  This will crash the compiling process.
#endif

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

#if defined( MACHINE_DEC ) \
  || ( defined(MACHINE_Linux) && defined(COMPILER_Intel8) )
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
#else
  subroutine burn_random(n)
    integer :: n
    integer :: i
    real*8  :: x
    integer, save :: first_call=1

#ifdef USE_ESMF
    if ( first_call .ne. 0 ) then
      write(6,*) ' ***********************************************'
      write(6,*) ' Warning: slow implementation of burn_random()  '
      write(6,*) ' on this platform.  Better performance can be   '
      write(6,*) ' achieved by using a recursion relation for most'
      write(6,*) ' random number generators. (contact Tom Clune)  '
      write(6,*) ' ***********************************************'
      first_call = 0
    endif
#endif

    do i = 1, n
      x = RANDU(x)
    end do

  end subroutine burn_random
#endif

end module RANDOM
