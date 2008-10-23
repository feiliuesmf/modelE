#include "rundeck_opts.h"
!@sum  UTILDBL Model Independent Utilities
!@auth Original Development Team
!@ver  1.0
!@cont THBAR,QSAT,DQSATDT,FILEMANAGER,READT,DREAD,MREAD

#if ( defined USE_ESMF )  || ( defined USE_MPP )
#define USE_MPI
#endif

      FUNCTION THBAR (X,Y)
!@sum  THBAR calculates mean temperature used in vertical differencing
!@auth Gary Russell, Jean Lerner, Arakawa
!@ver  1.0
C****
C**** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
C****              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
C****      g(x) is replaced by a rational function
C****           (a+bx+cxx+dxxx+xxxx)/(e+fx+gxx)
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
      REAL*8 :: THBAR                !@var THBAR  averaged temperature
      Q=X/Y
      AL=(A+Q*(B+Q*(C+Q*(D+Q))))/(E+Q*(F+G*Q))
      THBAR=X*AL
      RETURN
      END

      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : mrat,rvap,tf
      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
C**** LH = 0.5*(LH(0)+LH(t))
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(LH*(B-C/max(130.d0,TM)))/PR
      RETURN
      END

      FUNCTION DQSATDT (TM,LH)
!@sum  DQSATDT calculates change of sat. vapour mixing ratio with temp.
!@auth Gary Russell
!@ver  1.0
C**** Note that d(qsat)/dt = qsat * lh * c / T*T
C**** Only the factor of qsat is given here
      USE CONSTANT, only : rvap
      IMPLICIT NONE
!@var C coefficient for QSAT
      REAL*8, PARAMETER :: C = 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in DQSATDT is the actual LH at TM i.e. LH=LH(TM)
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: DQSATDT         !@var DQSATDT d(qsat)/dT factor only.
      DQSATDT = LH*C/(TM*TM)    ! * QSAT(TM,LH,PR)
      RETURN
      END

      FUNCTION SLP(PS,TAS,ZS)
!@sum SLP estimates sea level pressure in the presence of topography
!@+   for better match to reanalyses.
      USE CONSTANT, only: bmoist, grav, rgas, by3
      IMPLICIT NONE
!@var PS surface pressure (mb)
!@var TAS surface temperature (K)
!@var ZS surface elevation (m)
      REAL*8, INTENT(IN) :: PS, TAS, ZS
      REAL*8 :: SLP, TSL, BETA, BZBYT, GBYRB, TASn

      IF (ZS.ne.0.) THEN
        TSL= TAS+BMOIST*ZS
        TASn=TAS
        BETA=BMOIST
        IF (TAS < 290.5 .and. TSL > 290.5) BETA= (290.5d0 - TAS)/ZS
        IF (TAS > 290.5 .and. TSL > 290.5) TASn = 0.5*(290.5d0 + TAS)
        IF (TAS < 255) TASn = 0.5*(255d0 + TAS)
        BZBYT=BETA*ZS/TASn
        GBYRB=GRAV/(RGAS*BETA)
        IF (BETA > 1d-6 ) THEN
          SLP=PS*(1.+BZBYT)**GBYRB
        ELSE
          SLP=PS*EXP((1.-0.5*BZBYT+BZBYT**(2.*by3))*GBYRB*BZBYT)
        END IF
      ELSE
        SLP=PS
      END IF
      RETURN
      END FUNCTION SLP

      MODULE FILEMANAGER
!@sum  FILEMANAGER keeps data concerning the files and unit numbers
!@ver  2.0
      implicit none
      save
      private

      public openunit, closeunit, print_open_units,findunit
      public nameunit

      interface openunit
        module procedure openunit_0d
        module procedure openunit_1d
      end interface openunit

      interface closeunit
        module procedure closeunit_0d
        module procedure closeunit_1d
      end interface closeunit

!@param MINUNIT, MAXUNIT - min and max unit number allowed
      integer, parameter :: MINUNIT = 50, MAXUNIT = 128
      integer, parameter :: TOTALUNITS = MAXUNIT-MINUNIT+1

      type UnitStr
        logical in_use                     ! is the unit in use or not
        character*16 filename              ! the name on the file
      end type UnitStr

      type (UnitStr) :: Units( MINUNIT:MAXUNIT ) = UnitStr(.false.," ")

      contains

      function nameunit( unit )
!@sum nameunit returns the name of the file corresponding to unit <unit>
      implicit none
      character*16 nameunit
      integer, intent(in) :: unit

      if ( unit>MAXUNIT .or. unit<MINUNIT
     $     .or. .not. Units(unit)%in_use ) then
        write(6,*) "FILEMANAGER: asked name of a wrong unit: ",unit
        call stop_model("FILEMANAGER: asked name of a wrong unit",255)
      endif
      nameunit = Units(unit)%filename
      end function nameunit


      subroutine findunit( unit )
!@sum findunit finds available unit
      implicit none
      integer, intent(out) :: unit
      logical :: opened

!ia next line is needed for Absoft compiler on Linux (compiler bug?)
      unit = 0
      do unit=MINUNIT,MAXUNIT
        if ( unit > 98 .and. unit < 103 ) cycle
        !if ( .not. Units(unit)%in_use ) return
        inquire( unit=unit, opened=opened )
        if ( .not. opened ) return
      enddo
      write(6,*) "FILEMANAGER: Maximum file number reached"
      call print_open_units
      call stop_model("FILEMANAGER: Maximum file number reached",255)
      end subroutine findunit


      subroutine openunit_0d( filename, iunit, qbin, qold )
!@sum openunit opens the file <filename> and returns its unit in <unit>
      implicit none
!@var unit - unit of opened file
      integer, intent(out) :: iunit
!@var filename - name of the file to open
      character*(*), intent(in) :: filename
!@var qbin - .true. if file is binary (.false. by default)
!@var qold - .true. if file is old (.false. by default)
      logical, optional, intent(in) :: qbin, qold
      character*11 form, status
      integer name_len

C**** set default options
      form = "FORMATTED"
      status = "UNKNOWN"
C**** parse options
      if( present(qbin) ) then
         if( qbin ) form = "UNFORMATTED"
      end if
      if( present(qold) ) then
         if( qold ) status = "OLD"
      end if

      call findunit( iunit )

!dbg  write(6,*) "FILEMANAGER: Before Opening file ",trim(filename) !RKF debug
      if ( form == "FORMATTED" ) then
        open( iunit, FILE=filename, FORM=form, STATUS=status,
#ifdef CONVERT_BIGENDIAN
     *       CONVERT='BIG_ENDIAN',
#endif
     *       RECL=65534, ERR=10 )
      else
        open( iunit, FILE=filename, FORM=form, STATUS=status,
#ifdef CONVERT_BIGENDIAN
     *       CONVERT='BIG_ENDIAN',
#endif
     *       ERR=10 )
      endif

!dbg  write(6,*) "FILEMANAGER: Opened file ",trim(filename) !RKF debug

      Units(iunit)%in_use = .true.
      name_len = len_trim(filename)
      Units(iunit)%filename = filename( max(1,name_len-15) : name_len )
      return

 10   write(6,*) "FILEMANAGER: Error opening file ",trim(filename)
      call stop_model('FILEMANAGER: FILE OPENING ERROR',255)
      end subroutine openunit_0d

      subroutine closeunit_0d( unit )
!@sum closeunit closes the file the file corresponding to unit <unit>
      implicit none
!@var unit - unit of the file to close
      integer, intent(in) :: unit

      if ( unit>MAXUNIT .or. unit<MINUNIT
     $     .or. .not. Units(unit)%in_use ) then
        write(6,*) "FILEMANAGER: attempt to close wrong unit: ",unit
        call stop_model("FILEMANAGER: attempt to close wrong unit",255)
      endif

      close( unit )
      Units(unit)%in_use = .false.
      Units(unit)%filename = " "
      end subroutine closeunit_0d


      subroutine openunit_1d( filename, iunit, qbin, qold )
!@sum  openunit_1d sets unit number for requested files and opens them
!@auth Gavin Schmidt, simplified by I. Aleinov
!@ver  1.0
!@var unit - unit of opened file
      integer, intent(out), dimension(:) :: iunit
!@var filename - name of the file to open
      character*(*), intent(in), dimension(:) :: filename
!@var qbin - .true. if file is binary (.false. by default)
!@var qold - .true. if file is old (.false. by default)
      logical, optional, intent(in) :: qbin, qold
      integer i !@var i loop variable

      do i=1,size(filename)
        call openunit( filename(i), iunit(i), qbin, qold )
      end do
      end subroutine openunit_1d


      subroutine closeunit_1d( iunit )
!@sum  closeunit_1d closes files corresponding to units iunit
!@auth Gavin Schmidt, I. Aleinov
!@ver  1.0
      implicit none
!@var iunit unit numbers for files in current request
      integer, intent(in), dimension(:) :: iunit
!@var nreq number of file unit numbers
      integer i

      do i=1,size(iunit)
        call closeunit( iunit(i) )
      end do
      end subroutine closeunit_1d


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
      call stop_model('dREAD: READ ERROR',255)
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',TRIM(NAME(IUNIT))
      call stop_model('dREAD: No data found',255)
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
      call stop_model('mREAD: READ ERROR',255)
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      call stop_model('mREAD: No data found',255)
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
      call stop_model('tREAD: READ ERROR',255)
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      call stop_model('tREAD: No data found',255)
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

      subroutine READI (iunit,it,ain,it1,len4,iok)
!@sum  READI reads array surrounded by IT (for post processing)
!@auth  Original Development Team
!@ver   1.0
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT,LEN4
      INTEGER, INTENT(out) :: IT,it1,iok ! iok: 0='ok',1='not ok'
      real*4  AIN(LEN4)
      iok = 0
      read(iunit,end=555) it,ain,it1
      return
  555 iok=1
      return
      end subroutine readi

      subroutine WRITEI8 (iunit,it,aout,len8)
!@sum   WRITEI8 writes real*8 array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT       !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT          !@var  IT time, 1st & last word
      INTEGER, INTENT(IN) :: LEN8        !@var  LENGTH size of array
      REAL*8,  INTENT(IN) :: AOUT(LEN8)  !@var  AOUT   real*8 array

      write (iunit) it,aout,it
      call sys_flush(iunit)
      write (6,*) "Wrote to file ",TRIM(NAME(IUNIT)),", time=",it
      return
      END subroutine WRITEI8

      subroutine READI8 (iunit,it,ain,it1,len8,iok)
!@sum  READI reads array surrounded by IT (for post processing)
!@auth  Original Development Team
!@ver   1.0
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT,LEN8
      INTEGER, INTENT(out) :: IT,it1,iok ! iok: 0='ok',1='not ok'
      real*8  AIN(LEN8)
      iok = 0
      read(iunit,end=555) it,ain,it1
      return
  555 iok=1
      return
      end subroutine readi8

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
      if(it1 .le. it) go to 30
   10 write(6,*) "Starting a new file ",TRIM(NAME(IUNIT)),", time=",it
      rewind iunit
      return

   20 read (iunit,end=35,err=50) it1,(it2,n=1,len4+1)
   30 if (it2 .ne. it1) then
        write(6,*) 'file ',TRIM(NAME(IUNIT)),' damaged: it/it1/it2=',
     *    it,it1,it2
        call stop_model('io_POS: damaged file',255)
      end if
      if (it1 .le. it) go to 20
      it1=it1-itdif
   35 backspace iunit
      if (it1+itdif .le. it) go to 40
      write (6,*) "positioned ",TRIM(NAME(IUNIT)),", it1/itime=",it1,it
      return
   40 write (6,*) "file ",TRIM(NAME(IUNIT))," too short, it1/it=",it1,it
      call stop_model('io_POS: file too short',255)
   50 write (6,*) "Read error on: ",TRIM(NAME(IUNIT)),", it1/it=",it1,it
      call stop_model('io_POS: read error',255)
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

!$OMP PARALLEL DO PRIVATE (L,J,I) SHARED (QCHECK3)
      DO L=1,LN
      DO J=1,JN
      DO I=1,IN
        IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.) .or.
     *       ABS(A(I,J,L)) .gt.HUGE(A(I,J,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,L,A(I,J,L),'after ',SUBR
          IF (J.LT.JN.AND.J.GT.1) QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) call stop_model('CHECK3',255)
      RETURN
      END SUBROUTINE CHECK3

      SUBROUTINE CHECK3B(A,I1,I2,J1,J2,NJPOL,LN,SUBR,FIELD)
!@sum  CHECK3B Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: I1,I2,J1,J2,NJPOL,LN
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(I1:I2,J1:J2,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK3 = .FALSE.
      INTEGER I,J,L !@var I,J,L loop variables

!$OMP PARALLEL DO PRIVATE (L,J,I) SHARED (QCHECK3)
      DO L=1,LN
      DO J=J1+NJPOL,J2-NJPOL
      DO I=I1,I2
        IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.) .or.
     *       ABS(A(I,J,L)) .gt.HUGE(A(I,J,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,L,A(I,J,L),'after ',SUBR
          QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) call stop_model('CHECK3',255)
      RETURN
      END SUBROUTINE CHECK3B

      SUBROUTINE CHECK3C(A,LN,I1,I2,J1,J2,NJPOL,SUBR,FIELD)
!@sum  CHECK3B Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: LN,I1,I2,J1,J2,NJPOL
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(LN,I1:I2,J1:J2),INTENT(IN) :: A
      LOGICAL :: QCHECK3 = .FALSE.
      INTEGER I,J,L !@var I,J,L loop variables

!$OMP PARALLEL DO PRIVATE (L,J,I) SHARED (QCHECK3)
      DO J=J1+NJPOL,J2-NJPOL
      DO I=I1,I2
      DO L=1,LN
        IF (.NOT.(A(L,I,J).GT.0..OR.A(L,I,J).LE.0.) .or.
     *       ABS(A(L,I,J)) .gt.HUGE(A(L,I,J)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',L,I,J,A(L,I,J),'after ',SUBR
          QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) call stop_model('CHECK3',255)
      RETURN
      END SUBROUTINE CHECK3C

      SUBROUTINE CHECK4(A,IN,JN,KN,LN,SUBR,FIELD)
!@sum  CHECK4 Checks for NaN/INF in real 4-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,KN,LN size of 4-D array
      INTEGER, INTENT(IN) :: IN,JN,KN,LN
!@var SUBR identifies where CHECK4 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(IN,JN,KN,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK4 = .FALSE.
      INTEGER I,J,K,L !@var I,J,K,L loop variables

!$OMP PARALLEL DO PRIVATE (L,K,J,I) SHARED (QCHECK4)
      DO L=1,LN
      DO K=1,KN
      DO J=1,JN
      DO I=1,IN
        IF (.NOT.(A(I,J,K,L).GT.0..OR.A(I,J,K,L).LE.0.) .or.
     *       ABS(A(I,J,K,L)) .gt.HUGE(A(I,J,K,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,K,L,A(I,J,K,L),'after ',SUBR
          IF (J.LT.JN.AND.J.GT.1) QCHECK4 = .TRUE.
        END IF
      END DO
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK4) call stop_model('CHECK4',255)
      RETURN
      END SUBROUTINE CHECK4

      SUBROUTINE CHECK4B(A,I1,I2,J1,J2,NJPOL,KN,LN,SUBR,FIELD)
!@sum  CHECK4 Checks for NaN/INF in real 4-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,KN,LN size of 4-D array
      INTEGER, INTENT(IN) :: I1,I2,J1,J2,NJPOL,KN,LN
!@var SUBR identifies where CHECK4 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(I1:I2,J1:J2,KN,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK4 = .FALSE.
      INTEGER I,J,K,L !@var I,J,K,L loop variables

!$OMP PARALLEL DO PRIVATE (L,K,J,I) SHARED (QCHECK4)
      DO L=1,LN
      DO K=1,KN
      DO J=J1+NJPOL,J2-NJPOL
      DO I=I1,I2
        IF (.NOT.(A(I,J,K,L).GT.0..OR.A(I,J,K,L).LE.0.) .or.
     *       ABS(A(I,J,K,L)) .gt.HUGE(A(I,J,K,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,K,L,A(I,J,K,L),'after ',SUBR
          QCHECK4 = .TRUE.
        END IF
      END DO
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK4) call stop_model('CHECK4',255)
      RETURN
      END SUBROUTINE CHECK4B

      function unit_string (pow10,ending)
!@sum Construct a units string with nice properties (no embedded blanks)
!@auth G. Schmidt, J. Lerner
C**** If a trailing ')' is supplied, it is assumed that a leading
C****      '(' is required, so it is inserted
      implicit none
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


      subroutine stop_model( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
      USE PARAM
      implicit none
!@var message an error message (reason to stop)
      character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
      integer, intent(in) :: retcode
      integer, parameter :: iu_err = 9
      integer :: mpi_err
      integer :: rank
#ifdef USE_MPI
#include "mpi_defs.h"
#include "mpif.h"
#endif

#ifdef USE_MPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
#else
      rank =0
#endif
      if (rank == 0) then
        if ( retcode >= 0 ) ! skip writing status file for retcode<0
     &       call write_run_status( message, retcode )
        write (6,'(//2(" ",132("*")/))')
        write (6,*) ' Program terminated due to the following reason:'
        write (6,*) ' >>  ', message, '  <<'
        write (6,'(/2(" ",132("*")/))')
      endif

      call sys_flush(6)

      if ( retcode > 13 ) then
#ifdef USE_MPI
        call mpi_abort(MPI_COMM_WORLD, retcode,iu_err)
#else
        call sys_abort
#endif
      else
#ifdef USE_MPI
        call mpi_finalize(mpi_err)
#endif
        call exit_rc (0)
      endif

      end subroutine stop_model


      subroutine write_run_status( message, retcode )
      implicit none
      character*(*), intent (in) :: message
      integer, intent(in) :: retcode
      integer, parameter :: iu_status = 9
      character*10 :: form_str
      integer num_digits

      ! construct format string in such a way that retcode is printed
      ! at the beginning of the line with no extra spaces
      if ( retcode .ne. 0 ) then
        num_digits = log10( real( abs(retcode), kind(1.d0) ) ) + 1
      else
        num_digits = 1
      endif
      if ( retcode < 0 ) num_digits = num_digits + 1

      write(form_str,"('(I',I1,')')") num_digits

      open( iu_status, file='run_status', form='FORMATTED',
     &     status='UNKNOWN', ERR=10 )
      write( iu_status, form_str, ERR=10 ) retcode
      write( iu_status, '(A)', ERR=10 ) message
      close( iu_status )

      return
 10   continue
      write( 0, * ) "ERROR: Can't write to the run_status file"
      write( 0, * ) "STATUS:", message
      end subroutine write_run_status

#ifdef NEW_PGRAD_PBL
      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@sum  This version works for a nonorthogonal grid
!@auth M. Kelley
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas,radius
      USE MODEL_COM, only : im,jm,t,p,zatmo,sig,byim
      USE GEOM, only : ddx_ci,ddx_cj,ddy_ci,ddy_cj
     &     ,cosip,sinip
      USE DYNAMICS, only : phi,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      USE DOMAIN_DECOMP, only : grid, GET, HALO_UPDATE
      IMPLICIT NONE
      REAL*8 by_rho1,dpx1,dpy1,dpx0,dpy0,hemi
      real*8 :: dpsi,dpsj,dg1i,dg1j,dgsi,dgsj
      real*8 :: dpsdx,dpsdy,dg1dx,dg1dy,dgsdx,dgsdy
      INTEGER I,J,K,IP1,IM1,J1,IM1S,IP1S,IP1E

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, I_0, I_1, I_0H, I_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      I_0H = GRID%I_STRT_HALO
      I_1H = GRID%I_STOP_HALO


C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the primary grids

      ! for dPdy/rho at non-pole grids
      CALL HALO_UPDATE(grid, P)
      CALL HALO_UPDATE(grid, PHI)
      CALL HALO_UPDATE(grid, ZATMO)

      IF(I_0H.LT.I_0) THEN      ! halo cells exist in I direction
        IM1S=I_0-1; IP1S=I_0+1; IP1E=I_1+1
      ELSE                      ! periodic
        IM1S=IM-1; IP1S=1; IP1E=IM
      ENDIF

      DO J=J_0S,J_1S
      IM1=IM1S
      I=IM1S+1
      DO IP1=IP1S,IP1E
        dpsi = p(ip1,j)-p(im1,j)
        dg1i = phi(ip1,j,1)-phi(im1,j,1)
        dgsi = zatmo(ip1,j)-zatmo(im1,j)
        dpsj = p(i,j+1)-p(i,j-1)
        dg1j = phi(i,j+1,1)-phi(i,j-1,1)
        dgsj = zatmo(i,j+1)-zatmo(i,j-1)
        dpsdx = dpsi*ddx_ci(i,j) + dpsj*ddx_cj(i,j)
        dpsdy = dpsi*ddy_ci(i,j) + dpsj*ddy_cj(i,j)
        dg1dx = dg1i*ddx_ci(i,j) + dg1j*ddx_cj(i,j)
        dg1dy = dg1i*ddy_ci(i,j) + dg1j*ddy_cj(i,j)
        dgsdx = dgsi*ddx_ci(i,j) + dgsj*ddx_cj(i,j)
        dgsdy = dgsi*ddy_ci(i,j) + dgsj*ddy_cj(i,j)
        by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(pmid(1,I,J))
        DPDX_BY_RHO(I,J)=dpsdx*sig(1)*by_rho1+dg1dx
        DPDX_BY_RHO_0(I,J)=dpsdx*by_rho1+dgsdx
        DPDY_BY_RHO(I,J)=dpsdy*sig(1)*by_rho1+dg1dy
        DPDY_BY_RHO_0(I,J)=dpsdy*by_rho1+dgsdy
        IM1=I
        I=IP1
      ENDDO
      ENDDO

c
c the following are for a latlon grid only
c
      IF (have_south_pole) THEN
        hemi = -1.; J1 = 2
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,1)  =dpx1*BYIM
        DPDY_BY_RHO(1,1)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,1)=dpx0*BYIM
        DPDY_BY_RHO_0(1,1)=dpy0*BYIM
      END IF

      If (have_north_pole) THEN
          hemi= 1.; J1=JM-1
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,JM)  =dpx1*BYIM
        DPDY_BY_RHO(1,JM)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,JM)=dpx0*BYIM
        DPDY_BY_RHO_0(1,JM)=dpy0*BYIM
      END IF

      END SUBROUTINE PGRAD_PBL

#else
      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@auth Ye Cheng
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE MODEL_COM, only : im,jm,t,p,zatmo,sig,byim
      USE GEOM, only : bydyp,bydxp,cosip,sinip
      USE DYNAMICS, only : phi,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE
      REAL*8 by_rho1,dpx1,dpy1,dpx0,dpy0,hemi
      INTEGER I,J,K,IP1,IM1,J1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)


C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the primary grids

      ! for dPdy/rho at non-pole grids
      CALL HALO_UPDATE(grid, P,   FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, ZATMO, FROM=SOUTH+NORTH)

      DO I=1,IM
        DO J=J_0S,J_1S
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDY_BY_RHO(I,J)=(100.*(P(I,J+1)-P(I,J-1))*SIG(1)*by_rho1
     2         +PHI(I,J+1,1)-PHI(I,J-1,1))*BYDYP(J)*.5d0
          DPDY_BY_RHO_0(I,J)=(100.*(P(I,J+1)-P(I,J-1))*by_rho1
     2         +ZATMO(I,J+1)-ZATMO(I,J-1))*BYDYP(J)*.5d0
        END DO
      END DO

      ! for dPdx/rho at non-pole grids

      DO J=J_0S,J_1S
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDX_BY_RHO(I,J)=(100.*(P(IP1,J)-P(IM1,J))*SIG(1)*by_rho1
     2         +PHI(IP1,J,1)-PHI(IM1,J,1))*BYDXP(J)*.5d0
          DPDX_BY_RHO_0(I,J)=(100.*(P(IP1,J)-P(IM1,J))*by_rho1
     2         +ZATMO(IP1,J)-ZATMO(IM1,J))*BYDXP(J)*.5d0
          IM1=I
          I=IP1
        END DO
      END DO

      ! at poles

      IF (haveLatitude(grid, J=1)) THEN
        hemi = -1.; J1 = 2
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,1)  =dpx1*BYIM
        DPDY_BY_RHO(1,1)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,1)=dpx0*BYIM
        DPDY_BY_RHO_0(1,1)=dpy0*BYIM
      END IF

      If (haveLatitude(grid, J=JM)) THEN
          hemi= 1.; J1=JM-1
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,JM)  =dpx1*BYIM
        DPDY_BY_RHO(1,JM)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,JM)=dpx0*BYIM
        DPDY_BY_RHO_0(1,JM)=dpy0*BYIM
      END IF

      END SUBROUTINE PGRAD_PBL
#endif

      SUBROUTINE CALC_PIJL(lmax,p,pijl)
!@sum  CALC_PIJL Fills in P as 3-D
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only : lm,ls1,psfmpt
C****
      USE DOMAIN_DECOMP, Only : grid, GET
      implicit none
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: pijl
      integer :: l,lmax

      do l=1,ls1-1
        pijl(:,:,l) = p(:,:)
      enddo
      do l=ls1,lmax
        pijl(:,:,l) = PSFMPT
      enddo
      return
      end subroutine calc_pijl

      SUBROUTINE CALC_AMPK(LMAX)
!@sum  CALC_AMPK calculate air mass and pressure arrays
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,p
      USE DYNAMICS, only : plij,pdsig,pmid,pk,pedn,pek,sqrtp,am,byam
      USE DOMAIN_DECOMP, Only : grid, GET, HALO_UPDATE, SOUTH
      IMPLICIT NONE

      INTEGER :: I,J,L  !@var I,J,L  loop variables
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max. level for update
      REAL*8, DIMENSION(LMAX) :: PL,AML,PDSIGL,PMIDL
      REAL*8, DIMENSION(LMAX+1) :: PEDNL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, I_0H, I_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_HALO= J_0H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

C**** Calculate air mass, layer pressures, P**K, and sqrt(P)
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure. Routine should be called with LMAX=LM at start, and
C**** subsequentaly with LMAX=LS1-1
C**** Note Air mass is calculated in (kg/m^2)

C**** Fill in polar boxes
      IF (have_south_pole) P(2:IM,1) = P(1,1)
      IF (have_north_pole) P(2:IM,JM)= P(1,JM)
      Call HALO_UPDATE(grid, P, FROM=SOUTH)

!$OMP  PARALLEL DO PRIVATE (I,J,L,PL,AML,PDSIGL,PEDNL,PMIDL)
      DO J=J_0H,J_1 ! filling halo for P is faster than PDSIG
        DO I=I_0H,I_1H

          CALL CALC_VERT_AMP(P(I,J),LMAX,PL,AML,PDSIGL,PEDNL,PMIDL)

          DO L=1,MIN(LMAX,LM)
            PLIJ (L,I,J) = PL    (L)
            PDSIG(L,I,J) = PDSIGL(L)
            PMID (L,I,J) = PMIDL (L)
            PEDN (L,I,J) = PEDNL (L)
            AM   (L,I,J) = AML   (L)
            PK   (L,I,J) = PMIDL (L)**KAPA
            PEK  (L,I,J) = PEDNL (L)**KAPA
            BYAM (L,I,J) = 1./AM(L,I,J)
          END DO

          IF (LMAX.ge.LM) THEN
            PEDN(LM+1:LMAX+1,I,J) = PEDNL(LM+1:LMAX+1)
            PEK (LM+1:LMAX+1,I,J) = PEDN(LM+1:LMAX+1,I,J)**KAPA
          END IF
          SQRTP(I,J) = SQRT(P(I,J))
        END DO
      END DO
!$OMP  END PARALLEL DO

      RETURN
      END SUBROUTINE CALC_AMPK

      SUBROUTINE CALC_AMP(p,amp)
!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
!@auth Jean Lerner/Max Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,dsig,psf,ptop
      USE GEOM, only : axyp
C****
      USE DOMAIN_DECOMP, Only : grid, GET
      implicit none
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: amp
      integer :: j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

C
!$OMP  PARALLEL DO PRIVATE(J,L)
      DO L=1,LM
        IF(L.LT.LS1) THEN
ccc   do l=1,ls1-1
          do j=J_0,J_1
            amp(:,j,l) = p(:,j)*axyp(:,j)*dsig(l)
          enddo
ccc   enddo
        ELSE
ccc   do l=ls1,lm
          do j=J_0,J_1
            amp(:,j,l) = (psf-ptop)*axyp(:,j)*dsig(l)
          enddo
        END IF
ccc   enddo
      enddo
!$OMP  END PARALLEL DO
C
      return
C****
      end subroutine calc_amp

      SUBROUTINE CALC_TROP
!@sum  CALC_TROP (to calculate tropopause height and layer)
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,t
      USE GEOM, only : imaxj
      USE DIAG_COM, only : aij => aij_loc, ij_ptrop, ij_ttrop
      USE DYNAMICS, only : pk, pmid, PTROPO, LTROPO
      USE DOMAIN_DECOMP, Only : grid, GET
      IMPLICIT NONE
      INTEGER I,J,L,IERR
      REAL*8, DIMENSION(LM) :: TL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Find WMO Definition of Tropopause to Nearest L
!$OMP  PARALLEL DO PRIVATE (I,J,L,TL,IERR)
      do j=J_0,J_1
      do i=I_0,imaxj(j)
        do l=1,lm
          TL(L)=T(I,J,L)*PK(L,I,J)
        end do
        CALL TROPWMO(TL,PMID(1,I,J),PK(1,I,J),PTROPO(I,J),LTROPO(I,J)
     *       ,IERR)
        IF (IERR.gt.0) print*,"TROPWMO error: ",i,j
        AIJ(I,J,IJ_PTROP)=AIJ(I,J,IJ_PTROP)+PTROPO(I,J)
        AIJ(I,J,IJ_TTROP)=AIJ(I,J,IJ_TTROP)+TL(LTROPO(I,J))
      end do
      end do
!$OMP  END PARALLEL DO
      IF (have_south_pole) THEN
        PTROPO(2:IM,1) = PTROPO(1,1)
        LTROPO(2:IM,1) = LTROPO(1,1)
      END IF
      IF (have_north_pole) THEN
        PTROPO(2:IM,JM)= PTROPO(1,JM)
        LTROPO(2:IM,JM)= LTROPO(1,JM)
      END IF

      END SUBROUTINE CALC_TROP


      subroutine tropwmo(ptm1, papm1, pk, ptropo, ltropp,ierr)
!@sum  tropwmo calculates tropopasue height according to WMO formula
!@auth D. Nodorp/T. Reichler/C. Land
!@+    GISS Modifications by Jean Lerner/Gavin Schmidt
!@ver  1.0
!@alg  WMO Tropopause Definition
!@+
!@+ From A Temperature Lapse Rate Definition of the Tropopause Based on
!@+ Ozone, J. M. Roe and W. H. Jasperson, 1981
!@+
!@+ In the following discussion the lapse rate is defined as -dT/dz.
!@+
!@+ The main features of the WMO tropopause definition are as follows:
!@+ * The first tropopause (i.e., the conventional tropopause) is
!@+   defined as the lowest level at which the lapse rate decreases to 2
!@+   K/km or less, and the average lapse rate from this level to any
!@+   level within the next higher 2 km does not exceed 2 K/km.
!@+ * If above the first tropopause the average lapse rate between any
!@+   level and all higher levels within 1 km exceed 3 K/km, then a
!@+   second tropopause is defined by the same criterion as under the
!@+   statement above. This tropopause may be either within or above the
!@+   1 km layer.
!@+ * A level otherwise satisfying the definition of tropopause, but
!@+   occuring at an altitude below that of the 500 mb level will not be
!@+   designated a tropopause unless it is the only level satisfying the
!@+   definition and the average lapse rate fails to exceed 3 K/km over
!@+   at least 1 km in any higher layer.
!@+ * (GISS failsafe) Some cases occur when the lapse rate never falls
!@+   below 2 K/km. In such cases the failsafe level is that where the
!@+   lapse rate first falls below 3 K/km. If this still doesn't work
!@+   (ever?), the level is set to the pressure level below 30mb.
!@+
      USE MODEL_COM, only : klev=>lm
      USE CONSTANT, only : zkappa=>kapa,zzkap=>bykapa,grav,rgas
      USE DOMAIN_DECOMP, Only : grid, GET
      implicit none

      real*8, intent(in), dimension(klev) :: ptm1, papm1, pk
      real*8, intent(out) :: ptropo
      integer, intent(out) :: ltropp,ierr
      real*8, dimension(klev) :: zpmk, zpm, za, zb, ztm, zdtdz
!@param zgwmo min lapse rate (* -1) needed for trop. defn. (-K/km)
!@param zgwmo2 GISS failsafe minimum lapse rate (* -1) (-K/km)
!@param zdeltaz distance to check for lapse rate changes (km)
!@param zfaktor factor for caluclating height from pressure (-rgas/grav)
!@param zplimb min pressure at which to define tropopause (mb)
      real*8, parameter :: zgwmo  = -2d-3, zgwmo2=-3d-3,
     *     zdeltaz = 2000.0, zfaktor = -GRAV/RGAS, zplimb=500.
      real*8 zptph, zp2km, zag, zbg, zasum, zaquer, zptf
      integer iplimb,iplimt, jk, jj, kcount, ltset,l
      logical ldtdz
c****
c****  2. Calculate the height of the tropopause
c****  -----------------------------------------
      ltset = -999
      ierr=0
      iplimb=1
c**** set limits based on pressure
      do jk=2,klev-1
        if (papm1(jk-1).gt.600d0) then
          iplimb=jk
        else
          if (papm1(jk).lt.30d0) exit
        end if
      end do
      iplimt=jk
c****
c****  2.1 compute dt/dz
c****  -----------------
c****       ztm  lineare Interpolation in p**kappa
c****     gamma  dt/dp = a * kappa + papm1(jx,jk)**(kappa-1.)

      do jk=iplimb+1,iplimt       ! -1 ?????
        zpmk(jk)=0.5*(pk(jk-1)+pk(jk))

        zpm(jk)=zpmk(jk)**zzkap ! p mitte

        za(jk)=(ptm1(jk-1)-ptm1(jk))/(pk(jk-1)-pk(jk))
        zb(jk) = ptm1(jk)-(za(jk)*pk(jk))

        ztm(jk)=za(jk)*zpmk(jk)+zb(jk) ! T mitte
        zdtdz(jk)=zfaktor*zkappa*za(jk)*zpmk(jk)/ztm(jk)
      end do
c****
c****  2.2 First test: valid dt/dz ?
c****  -----------------------------
c****
      do 1000 jk=iplimb+1,iplimt-1

c**** GISS failsafe test
        if (zdtdz(jk).gt.zgwmo2.and.ltset.ne.1) then
          ltropp=jk
          ltset =1
        end if
c****
        if (zdtdz(jk).gt.zgwmo .and. ! dt/dz > -2K/km
     &       zpm(jk).le.zplimb) then ! zpm not too low
          ltropp = jk
          ltset = 1
c****
c****  2.3 dtdz is valid > something in German
c****  ----------------------------------------
c****    1.lineare in p^kappa (= Dieters neue Methode)

          zag = (zdtdz(jk)-zdtdz(jk+1))/
     &         (zpmk(jk)-zpmk(jk+1)) ! a-gamma
          zbg = zdtdz(jk+1) - zag*zpmk(jk+1) ! b-gamma
          if(((zgwmo-zbg)/zag).lt.0.) then
            zptf=0.
          else
            zptf=1.
          end if
          zptph = zptf*abs((zgwmo-zbg)/zag)**zzkap
          ldtdz=zdtdz(jk+1).lt.zgwmo
          if(.not.ldtdz) zptph=zpm(jk)
c****
c****  2.4 2nd test: dt/dz above 2km must not be lower than -2K/km
c****  -----------------------------------------------------------
c****
          zp2km = zptph + zdeltaz*zpm(jk)
     &         / ztm(jk)*zfaktor ! p at ptph + 2km
          zasum = 0.0           ! zdtdz above
          kcount = 0            ! number of levels above
c****
c****  2.5 Test until pm < p2km
c****  --------------------------
c****
          do jj=jk,iplimt-1
            if(zpm(jj).gt.zptph) cycle ! doesn't happen
            if(zpm(jj).lt.zp2km) goto 2000 ! ptropo valid
            zasum = zasum+zdtdz(jj)
            kcount = kcount+1
            zaquer = zasum/float(kcount) ! dt/dz mean
            if(zaquer.le.zgwmo) goto 1000 ! dt/dz above < 2K/1000
                                          ! discard it
          end do                ! test next level
          goto 2000
        endif
 1000 continue                  ! next level
 2000 continue

      if (ltset.eq.-999) then
        ltropp=iplimt-1  ! default = last level below 30mb
        print*,"In tropwmo ltropp not set, using default: ltropp ="
     *       ,ltropp
        write(6,'(12(I4,5F10.5,/))') (l,ptm1(l),papm1(l),pk(l),zdtdz(l)
     *       ,zpm(l),l=iplimb+1,iplimt-1)
        ierr=1
      end if
      ptropo = papm1(ltropp)
c****
      return
      end subroutine tropwmo

      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy returns the sum of kinetic and potential energy.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use GEOM, only: AXYP, AREAG
      use DOMAIN_DECOMP, only: grid, GLOBALSUM, get
      REAL*8 :: totalEnergy
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     KEIJ,PEIJ,TEIJ
      INTEGER :: I,J
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

      call conserv_PE(PEIJ)
      call conserv_KE(KEIJ)
      DO J=J_0,J_1
      DO I=I_0,I_1
        TEIJ(I,J)= (KEIJ(I,J) + PEIJ(I,J)*AXYP(I,J))/AREAG
      ENDDO
      ENDDO
      CALL GLOBALSUM(grid, TEIJ, totalEnergy, ALL=.true.)

      end function getTotalEnergy

      subroutine addEnergyAsDiffuseHeat(deltaEnergy)
!@sum  addEnergyAsDiffuseHeat adds in energy increase as diffuse heat.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use CONSTANT, only: sha, mb2kg
      use MODEL_COM, only: T, PSF, PMTOP, LM
      use DYNAMICS, only: PK
      use DOMAIN_DECOMP, only: grid, get
      implicit none
      real*8, intent(in) :: deltaEnergy

      real*8 :: ediff
      integer :: l
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      ediff = deltaEnergy / ((PSF-PMTOP)*SHA*mb2kg)
!$OMP  PARALLEL DO PRIVATE (L)
      do l=1,lm
        T(I_0:I_1,J_0:J_1,L)=T(I_0:I_1,J_0:J_1,L)
     &       -ediff/PK(L,I_0:I_1,J_0:J_1)
      end do
!$OMP  END PARALLEL DO

      end subroutine addEnergyAsDiffuseHeat

      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE (m^2/s^2) as heat locally
!@auth Gavin Schmidt
      USE MODEL_COM, only : t
      USE DYNAMICS, only : dke,kea,pk
      IMPLICIT NONE
C**** temporarily store latest KE in DKE array
      call calc_kea_3d(dke)
      dke(:,:,:) = dke(:,:,:) - kea(:,:,:)
      call addEnergyAsLocalHeat(DKE, T, PK)

      END SUBROUTINE DISSIP

C***** Add in dissipiated KE as heat locally
      subroutine addEnergyAsLocalHeat(deltaKE, T, PK)!, diagIndex)
!@sum  addEnergyAsLocalHeat adds in dissipated kinetic energy as heat locally.
!@sum  deltaKE is now on the A grid!!!
!@auth Tom Clune (SIVO)
!@ver  1.0
      use CONSTANT, only: SHA
      use GEOM, only: IMAXJ
      use MODEL_COM, only: LM
      use DOMAIN_DECOMP, only: grid, get, HALO_UPDATE, NORTH
     &     ,am_i_root
      use DIAG_COM, only: ajl => ajl_loc
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     deltaKE,T
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo) :: PK
c      integer, optional, intent(in) :: diagIndex

      integer :: i, j, k, l
      real*8 :: ediff
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff,K)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ediff = deltaKE(I,J,L) / (SHA*PK(L,I,J))
        T(I,J,L)=T(I,J,L)-ediff
c        if (present(diagIndex)) then
c          AJL(J,L,diagIndex) = AJL(J,L,diagIndex) - ediff
c        end if
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
      end subroutine addEnergyAsLocalHeat

C**** Calculate 3D vertical velocity (take SDA which has units
C**** mb*m2/s (but needs averaging over no. of leap frog timesteps)
C**** and convert to WSAVE, units of m/s):

      subroutine COMPUTE_WSAVE(wsave, sda, T, PK, PEDN, NIdyn)
      use CONSTANT, only: rgas, bygrav
      !use MODEL_COM, only: NIdyn
      use DOMAIN_DECOMP, only: grid, GET
      use GEOM, only: byaxyp
      use MODEL_COM, only: IM,JM,LM
      implicit none

      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO, lm),
     &     intent(out) :: WSAVE
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO, lm),
     &     intent(in)  :: SDA, T
      real*8, dimension(lm, grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     intent(in)  :: PK,PEDN
      integer, intent(in) :: NIdyn

      integer :: i, j, l
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!$OMP PARALLEL DO PRIVATE (l,i)
      do l=1,lm-1
        do j=J_0,J_1
        do i=I_0,I_1
         wsave(i,j,l)=2.*sda(i,j,l)*byaxyp(i,j)*
     &   rgas*0.5*(T(i,j,l)*pk(l,i,j)+T(i,j,l+1)*
     &   pk(l+1,i,j))*bygrav/(NIdyn*pedn(l+1,i,j))
        end do
        end do
      end do
!$OMP END PARALLEL DO

      end subroutine COMPUTE_WSAVE

      function nij_before_j0(j0)
#ifdef CUBED_SPHERE
      use resolution, only : im,jm
      use domain_decomp, only : tile
#else
      use geom, only : imaxj
#endif
      implicit none
      integer :: nij_before_j0,j0
#ifdef CUBED_SPHERE
      nij_before_j0 = im*((tile-1)*jm + (j0-1))
#else
      nij_before_j0 = SUM(IMAXJ(1:J0-1))
#endif
      return
      end function nij_before_j0

      function nij_after_j1(j1)
      use resolution, only : im,jm
#ifdef CUBED_SPHERE
      use domain_decomp, only : tile
#else
      use geom, only : imaxj
#endif
      implicit none
      integer :: nij_after_j1,j1
#ifdef CUBED_SPHERE
      nij_after_j1 = im*((6-tile)*jm + (jm-j1))
#else
      nij_after_j1 = SUM(IMAXJ(J1+1:JM))
#endif
      return
      end function nij_after_j1

      function nij_after_i1(i1)
      use resolution, only : im,jm
      implicit none
      integer :: nij_after_i1,i1
#ifdef CUBED_SPHERE
      nij_after_i1 = im-i1
#else
      nij_after_i1 = 0
#endif
      return
      end function nij_after_i1
