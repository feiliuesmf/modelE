C CMS SYSTEM ROUTINES EMULATION FOR IBM RS/6000
C
      SUBROUTINE CLOCKS(IHSC)
C THIS VERSION OF CLOCKS RETURNS PROCESS TIME OF USER AND
C SYSTEM TIME OF CHILD PROCESSES
C NOTE: MCLOCK IS REALLY IN HUNDREDTHS OF A SECOND, NOT SIXTIETHS.
      IHSC=-MCLOCK()
      RETURN
      END

      FUNCTION THBAR (X,Y)
C****
C**** TH-mean used for vertical differencing (Arakawa)
C**** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
C****              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
C****      g(x) is replaced by a rational function
C****           (a+bx+cxx+dxxx+cx**4)/(e+fx+gxx)
C****      approx.error <1.E-6 for x between .9 and 1.7
C****
CV    REAL*8 A,B,C,D,E,F,G,Q,AL
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (
     *  A=113.4977618974100d0,B=438.5012518098521d0,
     *  C=88.49964112645850d0,D=-11.50111432385882d0,
     *  E=30.00033943846368d0,F=299.9975118132485d0,
     *  G=299.9994728900967d0)
C     DATA A,B,C,D,E,F,G/113.4977618974100,438.5012518098521,
C    *  88.49964112645850,-11.50111432385882,
C    *  30.00033943846368,299.9975118132485,299.9994728900967/
      Q=X/Y
      AL=(A+Q*(B+Q*(C+Q*(D+Q))))/(E+Q*(F+G*Q))
      THBAR=X*AL
      RETURN
      END

      FUNCTION EXPBYK (X)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXPBYK=X**.286
      RETURN
      END

      SUBROUTINE DREAD (IUNIT,AIN,LENGTH,AOUT)
C****
C**** READ IN REAL*4 ARRAY AND CONVERT TO REAL*8
C****
      REAL*4 AIN(LENGTH)
      REAL*8 AOUT(LENGTH)
      READ (IUNIT) AIN
      DO 10 N=LENGTH,1,-1
   10 AOUT(N)=AIN(N)
      RETURN
      END

      SUBROUTINE MREAD (IUNIT,M,NSKIP,AIN,LENGTH,AOUT)
C****
C**** READ IN INTEGER & REAL*4 ARRAY AND CONVERT TO REAL*8
C****
      REAL*4 AIN(LENGTH),X
      REAL*8 AOUT(LENGTH)
      READ (IUNIT) M,(X,N=1,NSKIP),AIN
      DO 10 N=LENGTH,1,-1
   10 AOUT(N)=AIN(N)
      RETURN
      END

      SUBROUTINE READT (IUNIT,NSKIP,AIN,LENGTH,AOUT,IPOS)
C****
C**** READ IN TITLE & REAL*4 ARRAY AND CONVERT TO REAL*8
C****
      REAL*4 AIN(LENGTH),X
      REAL*8 AOUT(LENGTH)
      CHARACTER*80 TITLE
      DO 10 N=1,IPOS-1
   10 READ (IUNIT,END=920)
      READ (IUNIT,ERR=910,END=920) TITLE,(X,N=1,NSKIP),AIN
CCC   IF(LEN.LT.4*(20+NSKIP+LENGTH)) GO TO 930
      DO 100 N=LENGTH,1,-1
  100 AOUT(N)=AIN(N)
      WRITE(6,'('' Read from Unit '',I3,'':'',A80)') IUNIT,TITLE
      RETURN
  910 WRITE(6,*) 'READ ERROR ON UNIT',IUNIT
      STOP 'READ ERROR'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON UNIT',IUNIT
      STOP 'NO DATA TO READ'
C 930 WRITE(6,*) LEN/4,' RATHER THAN',20+NSKIP+LENGTH,' WORDS ON UNIT',
CC   *  IUNIT
CCC   STOP 'NOT ENOUGH DATA FOUND'
      END
      SUBROUTINE TIMER (MNOW,MINC,MSUM)
C****
C**** OUTPUT: MNOW (.01 S) = CURRENT CPU TIME
C****         MINC (.01 S) = TIME SINCE LAST CALL TO SUBROUTINE TIMER
C****         MSUM (.01 S) = MSUM + MINC
C****
      MNOW  = MCLOCK ()
      MINC  = MNOW - MLAST
      MSUM  = MSUM + MINC
      MLAST = MNOW
      RETURN
      END
