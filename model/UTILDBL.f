!@sum  UTILDBL Model Independent Utilities
!@auth Original Development Team
!@ver  1.0
!@cont THBAR,QSAT,DQSATDT,TRIDIAG,FILEMANAGER,READT,DREAD,MREAD

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

      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
!@auth Numerical Recipes
!@ver  1.0
      IMPLICIT NONE
      INTEGER, PARAMETER :: NMAX = 5000  !@var NMAX workspace
      INTEGER, INTENT(IN):: N         !@var N    dimension of arrays
      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(NMAX)             !@var GAM  work array
      INTEGER :: J                    !@var J    loop variable

      IF ( N > NMAX ) STOP "TRIDIAG: N > NMAX, increase NMAX"
      BET=B(1)
      IF (BET.eq.0) STOP "TRIDIAG: DENOMINATOR = ZERO"
      U(1)=R(1)/BET
      DO J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF (BET.eq.0) STOP "TRIDIAG: DENOMINATOR = ZERO"
        U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END


      MODULE FILEMANAGER
!@sum  FILEMANAGER keeps data concerning the files and unit numbers
!@ver  2.0
      implicit none
      save
      private

      public openunit, closeunit, print_open_units
      public openunits, closeunits, nameunit

!@param MINUNIT, MAXUNIT - min and max unit number allowed
      integer, parameter :: MINUNIT = 50, MAXUNIT = 98
      integer, parameter :: TOTALUNITS = MAXUNIT-MINUNIT+1

      type UnitStr
        logical in_use                     ! is the unit in use or not
        character*16 filename              ! the name on the file
      end type UnitStr

      type (UnitStr) :: Units( MINUNIT : MAXUNIT )
      data  Units  / TOTALUNITS*UnitStr(.false.," ") /

      contains

      function nameunit( unit )
!@sum nameunit returns the name of the file corresponding to unit <unit>
      implicit none
      character*16 nameunit
      integer, intent(in) :: unit

      if ( unit>MAXUNIT .or. unit<MINUNIT
     $     .or. .not. Units(unit)%in_use ) then
        write(6,*) "FILEMANAGER: asked name of a wrong unit: ",unit
        stop "FILEMANAGER: asked name of a wrong unit"
      endif
      nameunit = Units(unit)%filename
      end function nameunit


      subroutine findunit( unit )
!@sum findunit finds available unit
      implicit none
      integer, intent(out) :: unit

      do unit=MINUNIT,MAXUNIT
        if ( .not. Units(unit)%in_use ) return
      enddo
      write(6,*) "FILEMANAGER: Maximum file number reached"
      call print_open_units
      stop "FILEMANAGER: Maximum file number reached"
      end subroutine findunit


      subroutine openunit( filename, iunit, qbin, qold )
!@sum openunit opens the file <filename> and returns its unit in <unit>
      implicit none
!@var unit - unit of opened file
      integer, intent(out) :: iunit
!@var filename - name of the file to open
      character*(*), intent(in) :: filename
!@var qbin, qold - .true. if file is binary, old (.false. by default)
      logical, optional, intent(in) :: qbin, qold
      character*11 form, status
      integer name_len

C**** set default options
      form = "FORMATTED"
      status = "UNKNOWN"
C**** parse options
      if( present(qbin) .and. qbin ) form = "UNFORMATTED"
      if( present(qold) .and. qold ) status = "OLD"

      call findunit( iunit )

#ifdef CONVERT_BIGENDIAN
      open( iunit, FILE=filename, FORM=form, STATUS=status,
     *     CONVERT='BIG_ENDIAN', ERR=10 )
#else
      open( iunit, FILE=filename, FORM=form, STATUS=status,
     *     ERR=10 )
#endif

      Units(iunit)%in_use = .true.
      name_len = len_trim(filename)
      Units(iunit)%filename = filename( max(1,name_len-15) : name_len )
      return

 10   write(6,*) "FILEMANAGER: Error opening file ",trim(filename)
      stop 'FILEMANAGER: FILE OPENING ERROR'
      end subroutine openunit

      subroutine closeunit( unit )
!@sum closeunit closes the file the file corresponding to unit <unit>
      implicit none
!@var unit - unit of the file to close
      integer, intent(in) :: unit

      if ( unit>MAXUNIT .or. unit<MINUNIT
     $     .or. .not. Units(unit)%in_use ) then
        write(6,*) "FILEMANAGER: attempt to close wrong unit: ",unit
        stop "FILEMANAGER: attempt to close wrong unit"
      endif

      close( unit )
      Units(unit)%in_use = .false.
      Units(unit)%filename = " "
      end subroutine closeunit


      subroutine print_open_units
!@sum print_open_units prints info on open units (for debugging)
      implicit none
      integer unit

      write(6,*) "FILEMANAGER: Open Units:"
      do unit=MINUNIT,MAXUNIT
        if ( Units(unit)%in_use ) then
          write(6,*) "unit = ", unit, "file = ", Units(unit)%filename
        endif
      enddo
      end subroutine print_open_units


      SUBROUTINE OPENUNITS(FILENM,IUNIT,QBIN,NREQ)
!@sum  OPENUNITS sets unit number for requested files and opens them
!@auth Gavin Schmidt, simplified by I. Aleinov
!@ver  1.0
      IMPLICIT NONE
!@var IUNIT unit numbers for files in current request
      INTEGER, INTENT(OUT), DIMENSION(*) :: IUNIT
!@var FILENM name of file to open
      CHARACTER*(*), INTENT(IN), DIMENSION(*) :: FILENM
!@var QBIN true if binary file is to be opened (UNFORMATTED)
      LOGICAL, INTENT(IN), DIMENSION(*) :: QBIN
!@var NREQ number of file unit numbers requested
      INTEGER, INTENT(IN) :: NREQ
      INTEGER I !@var I loop variable

      DO I=1,NREQ
        call openunit( FILENM(I), IUNIT(I), QBIN(I) )
      END DO
      END SUBROUTINE OPENUNITS


      SUBROUTINE CLOSEUNITS( IUNIT, NREQ)
!@sum  CLOSEUNITS closes files corresponding to units IUNIT
!@auth Gavin Schmidt, I. Aleinov
!@ver  1.0
      IMPLICIT NONE
!@var IUNIT unit numbers for files in current request
      INTEGER, INTENT(IN), DIMENSION(*) :: IUNIT
!@var NREQ number of file unit numbers
      INTEGER, INTENT(IN) :: NREQ
      INTEGER I

      DO I=1,NREQ
        call closeunit( IUNIT(I) )
      END DO
      END SUBROUTINE CLOSEUNITS

      END MODULE FILEMANAGER

      SUBROUTINE DREAD (IUNIT,AIN,LENGTH,AOUT)
!@sum   DREAD   read in real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
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
      STOP 'dREAD: READ ERROR'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',TRIM(NAME(IUNIT))
      STOP 'dREAD: No data found'
      RETURN
      END

      SUBROUTINE MREAD (IUNIT,M,NSKIP,AIN,LENGTH,AOUT)
!@sum   MREAD   read in integer and real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
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
      STOP 'mREAD: READ ERROR'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      STOP 'mREAD: No data found'
      RETURN
      END

      SUBROUTINE READT (IUNIT,NSKIP,AIN,LENGTH,AOUT,IPOS)
!@sum   READT  read in title and real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
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
      STOP 'tREAD: READ ERROR'
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      STOP 'tREAD: No data found'
      END

      subroutine WRITEI (iunit,it,aout,len4)
!@sum   WRITEI  writes array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT       !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT          !@var  IT time, 1st & last word
      INTEGER, INTENT(IN) :: LEN4        !@var  LENGTH size of array
      REAL*4,  INTENT(IN) :: AOUT(LEN4)  !@var  AOUT   real*4 array

      write (iunit) it,aout,it
      call sys_flush(iunit)
      write (6,*) "Wrote to file ",TRIM(NAME(IUNIT)),", time=",it
      return
      END subroutine WRITEI

      subroutine io_POS (iunit,it,len4,itdif)
!@sum   io_POS  positions a seq. output file for the next write operat'n
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT !@var IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT,ITdif !@var IT,ITdif current time,dt
      INTEGER, INTENT(IN) :: LEN4 !@var LENGTH of array in words
      INTEGER :: IT1,IT2   !@var  time_tags at start,end of each record
      INTEGER :: N              !@var  N      loop variable

      read (iunit,end=10,err=50) it1,(it2,n=1,len4+1)
      if(it.ge.it1) go to 30
   10 write(6,*) "Starting a new file ",TRIM(NAME(IUNIT)),", time=",it
      rewind iunit
      return

   20 read (iunit,end=40,err=50) it1,(it2,n=1,len4+1)
   30 if (it2 .ne. it1) then
        write(6,*) 'file ',TRIM(NAME(IUNIT)),' damaged: it/it1/it2=',
     *    it,it1,it2
        stop 'io_POS: damaged file'
      end if
      if (it .ge. it1+itdif) go to 20
      write (6,*) "positioned ",TRIM(NAME(IUNIT)),", it1/itime=",it1,it
      return
   40 write (6,*) "file ",TRIM(NAME(IUNIT))," too short, it1/it=",it1,it
      stop 'io_POS: file too short'
   50 write (6,*) "Read error on: ",TRIM(NAME(IUNIT)),", it1/it=",it1,it
      stop 'io_POS: read error'
      END subroutine io_POS

      SUBROUTINE CHECK3(A,IN,JN,LN,SUBR,FIELD)
!@sum  CHECK3 Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: IN,JN,LN
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(IN,JN,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK3 = .FALSE.
      INTEGER I,J,L !@var I,J,L loop variables

      DO L=1,LN
      DO J=1,JN
      DO I=1,IN
        IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.)) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,L,A(I,J,L),'after ',SUBR
          IF (J.LT.JN.AND.J.GT.1) QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) STOP 'CHECK3'
      RETURN
      END SUBROUTINE CHECK3

#ifdef TEST_FM
      program test_FM
      use filemanager
      implicit none
      integer :: iu_AIC, iunit = 10

      call openunit("AIC",iu_AIC,.true.,.true.)
      print *, 'opened ', iu_AIC
      call closeunit(iu_AIC)
      call openunit("AIC",iu_AIC,.true.,.true.)
      print *, 'opened ', iu_AIC
      call closeunit(iu_AIC)

      stop

      open( iunit, FILE="AIC",
     $     FORM="UNFORMATTED", STATUS="OLD", ERR=10 )
      print *, 'opened '
      stop

          OPEN(IUNIT,FILE='AIC',FORM="UNFORMATTED",
     *       STATUS="OLD",ERR=10)

          print *, 'opened '


 10   write(6,*) "test_FM: Error opening file ", "AIC"
      stop 'FM: FILE OPENING ERROR'


      print *, 'closed ', iu_AIC
      end
#endif

      function unit_string (pow10,ending)
!@sum Construct a units string with nice properties (no embedded blanks)
!@auth G. Schmidt, J. Lerner
C**** If a trailing ')' is supplied, it is assumed that a leading
C****      '(' is required, so it is inserted
      character*(*) ending,unit_string
      character*10 tpow
      integer pow10

      tpow = ' '
      if(pow10.ne.0) then
        write(tpow,'(i3)') pow10
        if (index(ending,')') .ne.0) then
          tpow='(10^'//trim(adjustl(tpow))
        else
          tpow= '10^'//trim(adjustl(tpow))
        end if
      endif
      unit_string = adjustl(trim(tpow)//" "//trim(ending))
      return
      end function unit_string

