!@sum This file contains architecture specific code for SGI, IBM and Linux
      MODULE RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI version)
!@cont  RANDU, RINIT, RFINAL
      IMPLICIT NONE
      INTEGER, SAVE :: IX            !@var IX     random number seed

      CONTAINS

#if defined( MACHINE_SGI ) || defined( MACHINE_Linux )
      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
      REAL*8 X                       !@var X      dummy variable
      REAL*4 RAN                     !@fun RAN    SGI intrinsic func.
      REAL*8 :: RANDU                !@var RANDU  random number
      RANDU=RAN(IX)
      RETURN
      END FUNCTION RANDU
#elif defined( MACHINE_IBM )
      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
      REAL*8 X                       !@var X      dummy variable
      REAL*8 :: RANDU                !@var RANDU  random number
      INTEGER :: IY                  !@var IY     dummy integer
   10 IY=IX*65539
      SELECT CASE (IY)
      CASE (:-1)
         IY=(IY+2147483647)+1
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

      END MODULE RANDOM

#if defined( MACHINE_SGI ) || defined( MACHINE_IBM )
      SUBROUTINE GETTIME (MNOW)
!@sum  GETTIME returns current CPU time
!@auth Gary Russell
!@ver  1.0 (SGI version)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW !@var MNOW current CPU time (.01 s)
      INTEGER :: MCLOCK            !@var MCLOCK intrinsic function
C**** Note this routine is only here so that all MCLOCK related
C**** functions are in the same place, for ease of change on other
C**** platforms
      MNOW = MCLOCK()
      RETURN
      END SUBROUTINE GETTIME
#elif defined( MACHINE_Linux )
      SUBROUTINE GETTIME (MNOW)
!@sum  GETTIME returns current CPU time
!@auth Gary Russell
!@ver  1.0 (Absoft version)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW !@var MNOW current CPU time (.01 s)
      REAL*4 :: ETIME, TARR        !@var ETIME intrinsic function
      MNOW = NINT(ETIME(TARR)*100.)
      RETURN
      END SUBROUTINE GETTIME
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif

      SUBROUTINE exit_rc (code)
!@sum  exit_rc stops the run and sets a return code
!@auth Reto A Ruedy
!@ver  1.0 (SGI version)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: code !@var code return code set by user
#if defined( MACHINE_SGI )
      call exit(code)
#elif defined( MACHINE_IBM )
      call exit_(code)
#elif defined( MACHINE_Linux )
      call exit(code) !!! should check if it works for Absoft
#else
      None of supported architectures was specified.
      This will crash the compiling process.
#endif
      RETURN
      END SUBROUTINE exit_rc
