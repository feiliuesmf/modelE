      MODULE RANDOM
!@sum   RANDOM generates random numbers: 0<RANDom_nUmber<1
!@auth  Gary L. Russell
!@ver   1.0 (IBM version, needs 32-bit architecture for integers!)
!@cont  RANDU, RINIT, RFINAL
      IMPLICIT NONE
      INTEGER, SAVE :: IX            !@var IX     random number seed

      CONTAINS

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
