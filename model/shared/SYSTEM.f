!@sum This file contains architecture specific code for SGI, IBM, Linux, DEC
      MODULE RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI,IBM,Linux,DEC)
!@cont  RANDU, RINIT, RFINAL
      IMPLICIT NONE
      INTEGER, SAVE :: IX            !@var IX     random number seed

! Parameters used for "burning" sequences of random numbers
#if defined( MACHINE_DEC ) || ( MACHINE_Linux)
      INTEGER, PARAMETER :: A_linear = 69069
#elif defined(MACHINE_SGI) \
 || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
! don't know how RAN() is implemented, so burn
#else
      INTEGER, PARAMETER :: A_linear = 65539
#endif
      INTEGER, PARAMETER :: MAX_BITS = 31
      INTEGER, PARAMETER :: B_Half = 2**(MAX_BITS-1)
      INTEGER, PARAMETER :: B_linear = B_Half + (B_Half-1) ! 2147483647 = 2^31-1

      CONTAINS

#if defined(MACHINE_SGI) \
 || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) ) \
 || defined(MACHINE_DEC) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
!@calls RAN
      REAL*8 X                       !@var X      dummy variable
      REAL*4 RAN                     !@fun RAN    SGI intrinsic func.
      REAL*8 :: RANDU                !@var RANDU  random number
      RANDU=RAN(IX)
      RETURN
      END FUNCTION RANDU
#elif defined( MACHINE_IBM ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_NAG) ) \
 || ( defined(MACHINE_Linux) && defined(COMPILER_G95) ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_G95) ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
      REAL*8 X                       !@var X      dummy variable
      REAL*8 :: RANDU                !@var RANDU  random number
      INTEGER :: IY                  !@var IY     dummy integer
   10 IY=IX*A_linear
      SELECT CASE (IY)
      CASE (:-1)
         IY=(IY+B_linear)+1
      CASE (0)
         IX=1
         GO TO 10
      END SELECT
      IX=IY
      RANDU=DFLOAT(IY)*.465661287308D-9
      RETURN
      END FUNCTION RANDU
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif

      SUBROUTINE RINIT (INIT)
!@sum   RINIT sets the initial seed IX
      INTEGER, INTENT(IN)  :: INIT   !@var INIT   first random no. seed
      IX=INIT
      RETURN
      END SUBROUTINE RINIT

      SUBROUTINE RFINAL (IFINAL)
!@sum   RFINAL retrieves final seed value
      INTEGER, INTENT(OUT) :: IFINAL !@var IFINAL last random no. seed
      IFINAL=IX
      RETURN
      END SUBROUTINE RFINAL

#if defined( MACHINE_DEC ) \
 || ( defined(MACHINE_Linux) && defined(COMPILER_Intel8) )
      SUBROUTINE BURN_RANDOM(n)
!@sum  BURN_RANDOM burns a set number of random numbers. It is used to
!                  maintain bit-wise correspondence on parallel runs.
      implicit none
      integer, intent(in) :: n
      integer :: i
      real*8 x, randss
      integer :: a, b ! linear coefficient
      integer :: nn
      a = A_linear
      nn = n
      b = 1
      do i=1,MAX_BITS
        If (mod(nn,2) == 1) ix = ix * a + b
        b=(a+1)*b
        a=a*a
        nn=nn/2
        If (nn == 0) Exit
      end do
      return
      end subroutine burn_random
#else
      Subroutine burn_random(n)
      Integer :: n
      Integer :: i
      Real*8  :: x

      Write(6,*) ' ***********************************************'
      Write(6,*) ' Warning: slow implementation of burn_random()  '
      Write(6,*) ' on this platform.  Better performance can be   '
      Write(6,*) ' achieved by using a recursion relation for most'
      Write(6,*) ' random number generators. (contact Tom Clune)  '
      Write(6,*) ' ***********************************************'

      Do i = 1, n
        x = RANDU(x)
      End Do

      End Subroutine burn_random
#endif

      END MODULE RANDOM

      ! Use F90 system_clock for portable accuracy
      module GETTIME_MOD
      contains
      subroutine GETTIME(counter, count_rate_out)
      implicit none
      integer, intent(out) :: counter
      integer, intent(out), optional :: count_rate_out
      integer :: count_rate
      call system_clock(counter,count_rate)
      if( present(count_rate_out) ) count_rate_out = count_rate
      counter=100*(counter/count_rate)  ! force 100ths of seconds
      end subroutine GETTIME
      end module GETTIME_MOD

      SUBROUTINE sys_flush (unit)
!@sum system call to flush corresponding I/O unit
!@auth I. Aleinov
!@ver  1.0 (SGI,IBM,Linux,DEC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: unit !@var unit
      INTEGER status
#if defined(MACHINE_SGI)
      call flush(unit,status)
#elif defined(MACHINE_Linux) || defined(MACHINE_DEC) \
 || ( defined(MACHINE_MAC) && ! defined(COMPILER_XLF) )
      call flush(unit) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
      call flush_(unit)
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE sys_flush

      SUBROUTINE sys_signal (sig, prog)
!@sum system call to "signal"
!@auth I. Aleinov
!@ver  1.0 (SGI,IBM,Linux,Dec) !! should check if works with DEC !!
      IMPLICIT NONE
!@var unit signal number to catch
      INTEGER, INTENT(IN) :: sig
!@var prog handler subroutine for given signal
      EXTERNAL prog
#if defined(MACHINE_SGI) \
 || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) ) \
 || defined(MACHINE_DEC) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
      call signal( sig, prog, -1 )
#elif defined( MACHINE_IBM ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
      call signal( sig, prog )
#elif ( defined(MACHINE_MAC) && defined(COMPILER_NAG) ) \
 || ( defined(MACHINE_Linux) && defined(COMPILER_G95) ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_G95) )
      ! do nothing if "signal" is not supported by NAG
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE sys_signal


      SUBROUTINE sys_abort
!@sum system call to "abort" (to dump core)
#if ( defined(MACHINE_MAC) && defined(COMPILER_NAG) )
      use f90_unix_proc
#endif
      call abort
      END SUBROUTINE sys_abort


      subroutine nextarg( arg, opt )
!@sum returns next argument on the command line
!@+  arg - returned argument, or returns "" if no more arguments
!@+  if opt==1 return arg only if it is an option (starts with -)
#if ( defined(MACHINE_MAC) && defined(COMPILER_NAG) )
      use f90_unix_env
#endif
      implicit none
      character(*), intent(out) :: arg
      integer, external :: iargc
      integer, intent(in) :: opt
      integer, save :: count = 1
      if ( count > iargc() ) then
        arg=""
        return
      endif
      call getarg( count, arg )
      !if ( present(opt) ) then
        if ( opt == 1 .and. arg(1:1) .ne. '-' ) then
          arg=""
          return
        endif
      !endif
      count = count + 1
      return
      end subroutine nextarg



