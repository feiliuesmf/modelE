!@sum  UTILDBL Model Independent Utilities
!@auth Original Development Team
!@ver  1.0
!@cont THBAR,QSAT,DQSATDT,TRIDAIG,TIMER,FILEMANAGER,READT,DREAD,MREAD

      FUNCTION THBAR (X,Y)
!@sum  THBAR calculates mean temperature used in vertical differencing
!@auth Gary Russell, Jean Lerner, Arakawa
!@ver  1.0
C****
C**** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
C****              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
C****      g(x) is replaced by a rational function
C****           (a+bx+cxx+dxxx+cx**4)/(e+fx+gxx)
C****      approx.error <1.E-6 for x between .9 and 1.7
C****
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

      FUNCTION QSAT (TM,QL,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : mrat,rvap,tf
      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
C**** Note that if QL is considered to be a function of temperature, the
C**** correct argument in QSAT is the average QL from t=0 (C) to TM, ie.
C**** QL = 0.5*(QL(0)+QL(t))
      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(QL*(B-C/TM))/PR
      RETURN
      END

      FUNCTION DQSATDT (TM,QL)
!@sum  DQSATDT calculates change of sat. vapour mixing ratio with temp.
!@auth Gary Russell
!@ver  1.0
C**** Note that d(qsat)/dt = qsat * ql * c / T*T
C**** Only the factor of qsat is given here
      USE CONSTANT, only : rvap
      IMPLICIT NONE
!@var C coefficient for QSAT
      REAL*8, PARAMETER :: C = 1./RVAP        !2.166847d-3
C**** Note that if QL is considered to be a function of temperature, the
C**** correct argument in DQSATDT is the actual QL at TM i.e. QL=QL(TM)
      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
      REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap./sub. (J/kg)
      REAL*8 :: DQSATDT         !@var DQSATDT d(qsat)/dT factor only.
      DQSATDT = QL*C/(TM*TM)    ! * QSAT(TM,QL,PR)
      RETURN
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

      SUBROUTINE TIMEOUT (MBEGIN,MIN,MOUT)
!@sum  TIMEOUT redistributes timing info between counters
!@auth Gary Russell
!@ver  1.0 (SGI version)
      IMPLICIT NONE
!@var MBEGIN CPU time start of section (.01 s)
      INTEGER, INTENT(IN) :: MBEGIN
      INTEGER, INTENT(INOUT) :: MIN  !@var MIN counter to be added to
      INTEGER, INTENT(INOUT) :: MOUT !@var MOUT counter to be taken from
      INTEGER :: MINC                !@var MINC time since MBEGIN
      INTEGER :: MNOW                !@var MNOW current CPU time (.01 s)
      INTEGER :: MCLOCK              !@var MCLOCK intrinsic function

      MNOW  = MCLOCK()
      MINC  = MNOW - MBEGIN
      MIN   = MIN  + MINC
      MOUT  = MOUT - MINC
      RETURN
      END

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

      MODULE FILEMANAGER
!@sum  FILEMANAGER keeps data concerning the files and unit numbers
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

!@param IUNIT0 starting unit for GCM files
      INTEGER, PARAMETER :: IUNIT0 = 50
!@param IUNITMX maxiumum unit number allowed
      INTEGER, PARAMETER :: IUNITMX = 98
!@var NUNIT number of already opened files (initially zero)
      INTEGER, SAVE :: NUNIT = 0
!@var NAME name of unit number (or at least first ten characters)
      CHARACTER*10 :: NAME(IUNIT0:IUNITMX)

      CONTAINS

      SUBROUTINE GETUNIT(FILENM,IUNIT,QBIN,QOLD)
!@sum  GETUNIT sets a unit number for a requested file and opens it
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

!@var IUNIT unit number for file in current request
      INTEGER, INTENT(OUT) :: IUNIT
!@var FILENAME name of file to open
      CHARACTER*(*), INTENT(IN) :: FILENM
!@var QOLD,QBIN true if (old) binary file is to be opened (UNFORMATTED)
      LOGICAL, INTENT(IN) :: QBIN,QOLD

      IF (IUNIT0+NUNIT.gt.IUNITMX)
     *     STOP "Maximum file number reached"
C**** Set unit number
      IUNIT = IUNIT0 + NUNIT
C**** Open file
      IF (QOLD) THEN
        IF (QBIN) THEN
          OPEN(IUNIT,FILE=FILENM,FORM="UNFORMATTED",
     *       STATUS="OLD",ERR=10)
        ELSE
          OPEN(IUNIT,FILE=FILENM,FORM="FORMATTED",
     *       STATUS="OLD",ERR=10)
        END IF
      ELSE
        IF (QBIN) THEN
          OPEN(IUNIT,FILE=FILENM,FORM="UNFORMATTED",ERR=10)
        ELSE
          OPEN(IUNIT,FILE=FILENM,FORM="FORMATTED",ERR=10)
        END IF
      END IF
C**** set NAME for error tracking purposes
      NAME (IUNIT) = FILENM
C**** increment no of files
      NUNIT = NUNIT + 1
      RETURN
 10   WRITE(6,*) "Error opening file ",TRIM(FILENM)
      STOP 'FILE OPENING ERROR IN GETUNIT'
      END SUBROUTINE GETUNIT

      SUBROUTINE GETUNITS(FILENM,IUNIT,QBIN,NREQ)
!@sum  GETUNITS sets unit number for requested files and opens them
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

!@var NREQ number of file unit numbers requested
      INTEGER, INTENT(IN) :: NREQ
!@var IUNIT unit numbers for files in current request
      INTEGER, INTENT(OUT), DIMENSION(*) :: IUNIT
!@var FILENAME name of file to open
      CHARACTER*(*), INTENT(IN), DIMENSION(*) :: FILENM
!@var QBIN true if binary file is to be opened (UNFORMATTED)
      LOGICAL, INTENT(IN), DIMENSION(*) :: QBIN
      INTEGER I !@var I loop variable

      IF (IUNIT0+NUNIT+NREQ-1.gt.IUNITMX)
     *     STOP "Maximum file number reached"

      DO I=1,NREQ
C**** Set unit number
        IUNIT(I) = IUNIT0 + NUNIT
C**** Open file
        IF (QBIN(I)) THEN
          OPEN(IUNIT(I),FILE=FILENM(I),FORM="UNFORMATTED",
     *         STATUS="OLD",ERR=10)
        ELSE
          OPEN(IUNIT(I),FILE=FILENM(I),FORM="FORMATTED",
     *         STATUS="OLD",ERR=10)
        END IF
C**** set NAME for error tracking purposes
        WRITE(NAME(IUNIT(I)),'(A8,I2)') FILENM(1:8),I
C**** increment no of files
        NUNIT = NUNIT + 1
      END DO
      RETURN
 10   WRITE(6,*) "Error opening file ",TRIM(FILENM(I))
      STOP 'FILE OPENING ERROR IN GETUNITS'
      END SUBROUTINE GETUNITS

      SUBROUTINE CLOSEUNITS
!@sum  CLOSEUNIT closes all previously opened files
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE
      INTEGER IUNIT

      DO IUNIT=IUNIT0,IUNIT0+NUNIT-1
        CLOSE(IUNIT)
      END DO

      RETURN
      END SUBROUTINE CLOSEUNITS

      END MODULE FILEMANAGER

      SUBROUTINE DREAD (IUNIT,AIN,LENGTH,AOUT)
!@sum   DREAD   read in real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
      USE FILEMANAGER, only : NAME !@var NAME name of record being read
      IMPLICIT NONE
      INTEGER :: IUNIT                    !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      INTEGER :: N                        !@var  N      loop variable

      READ (IUNIT,ERR=910,END=920) AIN
C**** do transfer backwards in case AOUT and AIN are same workspace
      DO N=LENGTH,1,-1
        AOUT(N)=AIN(N)
      END DO
      WRITE(6,*) "Sucessful read from file ",NAME(IUNIT)
      RETURN
  910 WRITE(6,*) 'READ ERROR ON FILE ',TRIM(NAME(IUNIT))
      STOP 'READ ERROR IN DREAD'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',TRIM(NAME(IUNIT))
      STOP 'NO DATA TO READ'
      RETURN
      END

      SUBROUTINE MREAD (IUNIT,M,NSKIP,AIN,LENGTH,AOUT)
!@sum   MREAD   read in integer and real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
      USE FILEMANAGER, only : NAME !@var NAME name of record being read
      IMPLICIT NONE
      INTEGER :: IUNIT                    !@var  IUNIT  file unit number
      INTEGER, INTENT(OUT) :: M           !@var  M      initial integer
      INTEGER, INTENT(IN) :: NSKIP        !@var  NSKIP  words to skip
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      REAL*4 :: X                         !@var  X      dummy variable
      INTEGER :: N                        !@var  N      loop variable

      READ (IUNIT,ERR=910,END=920) M,(X,N=1,NSKIP),AIN
C**** do transfer backwards in case AOUT and AIN are same workspace
      DO N=LENGTH,1,-1
        AOUT(N)=AIN(N)
      END DO
      WRITE(6,*) "Sucessful read from file ",NAME(IUNIT)
      RETURN
  910 WRITE(6,*) 'READ ERROR ON FILE ',NAME(IUNIT)
      STOP 'READ ERROR IN MREAD'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      STOP 'NO DATA TO READ'
      RETURN
      END

      SUBROUTINE READT (IUNIT,NSKIP,AIN,LENGTH,AOUT,IPOS)
!@sum   READT  read in title and real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
      USE FILEMANAGER, only : NAME !@var NAME name of record being read
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT        !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: NSKIP    !@var  NSKIP  no. of R*4's to skip
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      INTEGER, INTENT(IN) :: IPOS  !@var  IPOS   no. of recs. to advance
      REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      REAL*4 :: X               !@var  X      dummy variable
      INTEGER :: N              !@var  N      loop variable
      CHARACTER*80 TITLE        !@var  TITLE  title of file record

      DO N=1,IPOS-1
        READ (IUNIT,END=920)
      END DO
      READ (IUNIT,ERR=910,END=920) TITLE,(X,N=1,NSKIP),AIN
C**** do transfer backwards in case AOUT and AIN are same workspace
      DO N=LENGTH,1,-1
        AOUT(N)=AIN(N)
      END DO
      WRITE(6,*) "Read from file ",TRIM(NAME(IUNIT)),": ",TRIM(TITLE)
      RETURN
  910 WRITE(6,*) 'READ ERROR ON FILE ',NAME(IUNIT)
      STOP 'READ ERROR IN READT'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      STOP 'NO DATA TO READ'
      END

      subroutine WRITEI (iunit,it,aout,len4)
!@sum   WRITEI  writes array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
      use FILEMANAGER, only : NAME !@var NAME name of record being read
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT        !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT         !@var  time_tag, 1st & last word
      INTEGER, INTENT(IN) :: LEN4         !@var  LENGTH size of array
      REAL*4,  INTENT(IN) :: AOUT(LEN4)   !@var  AOUT   real*4 array

      write (iunit) it,aout,it
      endfile iunit
      backspace iunit
      write (6,*) "Wrote to file ",TRIM(NAME(IUNIT)),", time=",it
      return
      END subroutine WRITEI

      subroutine io_POS (iunit,it,len4,itdif)
!@sum   io_POS  positions a seq. output file for the next write operat'n
!@auth  Original Development Team
!@ver   1.0
      use FILEMANAGER, only : NAME !@var NAME name of record being read
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT        !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT,ITdif   !@var  current time,time step
      INTEGER, INTENT(IN) :: LEN4       !@var  LENGTH of array in words
      INTEGER :: IT1,IT2   !@var  time_tags at start,end of each record
      INTEGER :: N                      !@var  N      loop variable

      read (iunit,end=10,err=50) it1,(it2,n=1,len4+1)
      if(it.gt.it1) go to 30
   10 write(6,*) "Starting a new file ",TRIM(NAME(IUNIT)),", time=",it
      rewind iunit
      return

   20 read (iunit,end=40,err=50) it1,(it2,n=1,len4+1)
   30 if (it2 .ne. it1) then
        write(6,*) 'file ',TRIM(NAME(IUNIT)),' damaged: it/it1/it2=',
     *    it,it1,it2
        stop 1
      end if
      if (it .ge. it1+itdif) go to 20
      write (6,*) "positioned ",TRIM(NAME(IUNIT)),", it1/itime=",it1,it
      return
   40 write (6,*) "file ",TRIM(NAME(IUNIT))," too short, it1/it=",it1,it
      stop 1
   50 write (6,*) "Read error on: ",TRIM(NAME(IUNIT)),", it1/it=",it1,it
      stop 1
      END subroutine io_POS
