   
      FUNCTION THBAR (X,Y)
C****
C**** TH-mean used for vertical differencing (Arakawa)
C**** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
C****              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
C****      g(x) is replaced by a rational function
C****           (a+bx+cxx+dxxx+cx**4)/(e+fx+gxx)
C****      approx.error <1.E-6 for x between .9 and 1.7
C****
!@sum   THBAR calculates mean temperature in the vertical
!@auth  Gary Russell, Jean Lerner
!@ver   1.0
      IMPLICIT NONE
!@var A,B,C,D,E,F,G   expansion coefficients for THBAR
      REAL*8, PARAMETER :: A=113.4977618974100d0
      REAL*8, PARAMETER :: B=438.5012518098521d0
      REAL*8, PARAMETER :: C=88.49964112645850d0
      REAL*8, PARAMETER :: D=-11.50111432385882d0
      REAL*8, PARAMETER :: E=30.00033943846368d0
      REAL*8, PARAMETER :: F=299.9975118132485d0
      REAL*8, PARAMETER :: G=299.9994728900967d0
      REAL*8 :: Q,AL                 !@var Q,AL   working variables
      REAL*8, INTENT(IN) :: X,Y      !@var X,Y    input temperatures
      REAl*8 :: THBAR                !@var THBAR  averaged temperature
      Q=X/Y
      AL=(A+Q*(B+Q*(C+Q*(D+Q))))/(E+Q*(F+G*Q))
      THBAR=X*AL
      RETURN
      END

      FUNCTION EXPBYK (X)
!@sum   EXPBYK exponentiates pressure by KAPA
!@auth  Compiler development made original version obsolete
!@ver   1.0
      USE CONSTANT, only : kapa
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: X      !@var X          input pressure
      REAL*8 :: EXPBYK             !@var EXPBYK     output P^KAPA
      EXPBYK=X**KAPA
      RETURN
      END

      FUNCTION QSAT_NEW (TM,QL,PR)
!@sum   QSAT calculates saturation vapour mixing ratio
!@auth  Gary Russell
!@ver   1.0
      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=3.797915d0
      REAL*8, PARAMETER :: B=7.93252d-6
      REAL*8, PARAMETER :: C=2.166847d-3
      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
      REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap. (J/kg)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8 :: QSAT_NEW        !@var QSAT sat. vapour mixing ratio
      QSAT_NEW = A*EXP(QL*(B-C/TM))/PR
      RETURN
      END

      SUBROUTINE DREAD (IUNIT,AIN,LENGTH,AOUT)
!@sum   DREAD   read in real*4 array and convert to real*8
!@auth  Reto Ruedy (should become obsolete)
!@ver   1.0
      IMPLICIT NONE
      INTEGER :: IUNIT                    !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      INTEGER :: N                        !@var  N      loop variable

      READ (IUNIT) AIN
C**** do transfer backwards in case AOUT and AIN are same workspace
      DO N=LENGTH,1,-1
         AOUT(N)=AIN(N)
      END DO
      RETURN
      END

      SUBROUTINE MREAD (IUNIT,M,NSKIP,AIN,LENGTH,AOUT)
!@sum   MREAD   read in integer and real*4 array and convert to real*8
!@auth  Reto Ruedy (should become obsolete)
!@ver   1.0
      IMPLICIT NONE
      INTEGER :: IUNIT                    !@var  IUNIT  file unit number
      INTEGER, INTENT(OUT) :: M           !@var  M      initial integer
      INTEGER, INTENT(IN) :: NSKIP        !@var  NSKIP  words to skip
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      REAL*4 :: X                         !@var  X      dummy variable
      INTEGER :: N                        !@var  N      loop variable

      READ (IUNIT) M,(X,N=1,NSKIP),AIN
C**** do transfer backwards in case AOUT and AIN are same workspace
      DO N=LENGTH,1,-1
         AOUT(N)=AIN(N)
      END DO
      RETURN
      END

      SUBROUTINE READT (IUNIT,NSKIP,AIN,LENGTH,AOUT,IPOS)
!@sum   READT  read in title and real*4 array and convert to real*8
!@auth  Reto Ruedy (should become obsolete)
!@ver   1.0
      IMPLICIT NONE
      INTEGER :: IUNIT             !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: NSKIP !@var  NSKIP  no. of chars. to skip
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      INTEGER, INTENT(IN) :: IPOS  !@var  IPOS   no. of recs. to advance
      REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      REAL*4 :: X                  !@var  X      dummy variable
      INTEGER :: N                 !@var  N      loop variable
      CHARACTER*80 TITLE           !@var  TITLE  title of file record

      DO N=1,IPOS-1
         READ (IUNIT,END=920)
      END DO
      READ (IUNIT,ERR=910,END=920) TITLE,(X,N=1,NSKIP),AIN
C**** do transfer backwards in case AOUT and AIN are same workspace
      DO N=LENGTH,1,-1
         AOUT(N)=AIN(N)
      END DO
      WRITE(6,'('' Read from Unit '',I3,'':'',A80)') IUNIT,TITLE
      RETURN
  910 WRITE(6,*) 'READ ERROR ON UNIT',IUNIT
      STOP 'READ ERROR'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON UNIT',IUNIT
      STOP 'NO DATA TO READ'
      END

      SUBROUTINE TIMER (MNOW,MINC,MSUM)
!@sum  TIMER keeps track of elapsed CPU time in hundredths of seconds
!@auth Gary Russell
!@ver  1.0 (SGI version)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: MNOW   !@var MNOW current CPU time (.01 s)
      INTEGER, INTENT(OUT) :: MINC   !@var MINC time since last call
      INTEGER, INTENT(INOUT) :: MSUM !@var MSUM running total
      INTEGER :: MCLOCK              !@var MCLOCK intrinsic function
      INTEGER, SAVE :: MLAST = 0     !@var MLAST  last CPU time

      MNOW  = MCLOCK ()
      MINC  = MNOW - MLAST
      MSUM  = MSUM + MINC
      MLAST = MNOW
      RETURN
      END

      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
!@auth Numerical Recipes
!@ver  1.0
      IMPLICIT NONE
      INTEGER :: N                 !@var N    dimension of arrays
      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
      REAL*8, INTENT(OUT) :: U(N)  !@var U    solution vector
      REAL*8 :: BET                !@var BET  work variable
      REAL*8 :: GAM(N)             !@var GAM  work array
      INTEGER :: J                 !@var J    loop variable

      BET=B(1)
      U(1)=R(1)/BET
      DO J=2,N
         GAM(J)=C(J-1)/BET
         BET=B(J)-A(J)*GAM(J)
         U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J=N-1,1,-1
         U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END

