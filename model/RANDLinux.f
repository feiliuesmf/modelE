C**** This file contains all SGI specific code
      MODULE RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Reto Ruedy
!@ver   1.0 (SGI version)
!@cont  RANDU, RINIT, RFINAL
      IMPLICIT NONE
      INTEGER, SAVE :: IX            !@var IX     random number seed

      CONTAINS

      FUNCTION RANDU (X)
!@sum   RANDU calculates a random number based on the seed IX
!@calls RAN
      REAL*8 X                       !@var X      dummy variable
      REAL*4 RAN                     !@fun RAN    SGI intrinsic func.
      REAL*8 :: RANDU                !@var RANDU  random number
      RANDU=RAN(IX)
      RETURN
      END FUNCTION RANDU

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

      SUBROUTINE GETTIME (MNOW)
!@sum  GETTIME returns current CPU time
!@auth Gary Russell
!@ver  1.0 (Absoft version)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW !@var MNOW current CPU time (.01 s)
      REAL*4 :: ETIME, TARR        !@var ETIME intrinsic function
C**** Note this routine is only here so that all MCLOCK related
C**** functions are in the same place, for ease of change on other
C**** platforms
      MNOW = NINT(ETIME(TARR)*100.)
      RETURN
      END SUBROUTINE GETTIME

      SUBROUTINE exit_rc (code)
!@sum  exit_rc stops the run and sets a return code
!@auth Reto A Ruedy
!@ver  1.0 (SGI version)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: code !@var return code set by user
C**** Note this routine is only here so that all MCLOCK related
C**** functions are in the same place, for ease of change on other
C**** platforms
      SELECT CASE (code)
      CASE (:10)
        stop     ! not used yet: rc=0
      CASE (11)
        call exit (11)
      CASE (12)
        call exit (12)
      CASE (13)
        call exit (13)
      CASE (14:250,252:)
        stop ! not used yet: rc=0
      CASE (251)
        call exit (251)
      END SELECT
      RETURN
      END SUBROUTINE exit_rc
