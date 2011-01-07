!@sum This file contains architecture specific code for SGI, IBM, Linux, DEC

! hack for PGI : use the same settings as G95 (except for iargc)
#ifdef COMPILER_PGI
#define COMPILER_G95
#endif

      SUBROUTINE exit_rc (code)
!@sum  exit_rc stops the run and sets a return code
!@auth Reto A Ruedy
!@ver  1.0 (SGI,IBM,Linux,DEC)
#if ( defined(COMPILER_NAG) )
      use f90_unix_proc
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: code !@var code return code set by user
#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC) \
 || ( defined(MACHINE_MAC) && ! defined(COMPILER_XLF) )
      call exit(code) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
      call exit_(code)
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE exit_rc

      SUBROUTINE sys_flush (unit)
!@sum system call to flush corresponding I/O unit
!@auth I. Aleinov
!@ver  1.0 (SGI,IBM,Linux,DEC)
#if ( defined(COMPILER_NAG) )
      use f90_unix_io
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: unit !@var unit
#if defined(MACHINE_SGI)
#if defined(COMPILER_G95)
      call flush(unit)
#else
      INTEGER status
      call flush(unit,status)
#endif
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
 || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) && ! defined(COMPILER_NAG) ) \
 || defined(MACHINE_DEC) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_Intel8) ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
      call signal( sig, prog, -1 )
#elif defined( MACHINE_IBM ) \
 || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
      call signal( sig, prog )
#elif ( defined(COMPILER_NAG) ) \
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
#if ( defined(COMPILER_NAG) )
      use f90_unix_proc
#endif
      call abort
      END SUBROUTINE sys_abort


      subroutine nextarg( arg, opt )
!@sum returns next argument on the command line
!@+  arg - returned argument, or returns "" if no more arguments
!@+  if opt==1 return arg only if it is an option (starts with -)
#if ( defined(COMPILER_NAG) )
      use f90_unix_env
#endif
      implicit none
      character(*), intent(out) :: arg
#if ((! defined(COMPILER_NAG) ) && (! defined(COMPILER_G95) )) || (defined COMPILER_PGI)
      integer, external :: iargc
#endif
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



