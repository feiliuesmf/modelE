!@sum  Model Independent Utilities from GCM UTILDBL file
!@auth Original Development Team
!@ver  1.0
!@cont QSAT,DQSATDT,FILEMANAGER,READT,STOP_MODEL


      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
!      USE CONSTANT, only : mrat,rvap,tf
      USE ent_const, only : mrat,rvap,tf
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
      !USE CONSTANT, only : rvap
      USE ent_const, only : rvap
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


      MODULE FILEMANAGER
!@sum  FILEMANAGER keeps data concerning the files and unit numbers
!@ver  2.0
      implicit none
      save
      private

      public openunit, closeunit, print_open_units
      public openunits, closeunits, nameunit

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

!ia next line is needed for Absoft compiler on Linux (compiler bug?)
      unit = 0
      do unit=MINUNIT,MAXUNIT
        if ( unit > 98 .and. unit < 103 ) cycle
        if ( .not. Units(unit)%in_use ) return
      enddo
      write(6,*) "FILEMANAGER: Maximum file number reached"
      call print_open_units
      call stop_model("FILEMANAGER: Maximum file number reached",255)
      end subroutine findunit


      subroutine openunit( filename, iunit, qbin, qold )
!@sum openunit opens the file <filename> and returns its unit in <unit>
      implicit none
!@var unit - unit of opened file
      integer, intent(out) :: iunit
!@var filename - name of the file to open
      character*(*), intent(in) :: filename
!@var qbin - .true. if file is binary (.false. by default)
!@var qold - .true. if file is old (.false. by default)
      logical, optional, intent(in) :: qbin, qold
      character*11 form, status, action
      integer name_len

C**** set default options
      form = "FORMATTED"
      status = "UNKNOWN"
      action = "READWRITE"
C**** parse options
      if( present(qbin) ) then
         if( qbin ) form = "UNFORMATTED"
      end if
      if( present(qold) ) then
         if( qold ) then
           status = "OLD"
           action = "READ"
         endif
      end if

      call findunit( iunit )

!      write(6,*) "FILEMANAGER: Before Opening file ",trim(filename) !RKF debug
      if ( form == "FORMATTED" ) then
        open( iunit, FILE=filename, FORM=form, STATUS=status,
     &       ACTION=action,
#ifdef CONVERT_BIGENDIAN
     *       CONVERT='BIG_ENDIAN',
#endif
     *       RECL=65534, ERR=10 )
      else
        open( iunit, FILE=filename, FORM=form, STATUS=status,
     &       ACTION=action,
#ifdef CONVERT_BIGENDIAN
     *       CONVERT='BIG_ENDIAN',
#endif
     *       ERR=10 )
      endif

      write(6,*) "FILEMANAGER: Opened file ",trim(filename) !RKF debug

      Units(iunit)%in_use = .true.
      name_len = len_trim(filename)
      Units(iunit)%filename = filename( max(1,name_len-15) : name_len )
      return

 10   write(6,*) "FILEMANAGER: Error opening file ",trim(filename)
      call stop_model('FILEMANAGER: FILE OPENING ERROR',255)
      end subroutine openunit

      subroutine closeunit( unit )
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

      subroutine stop_model( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
!@var message an error message (reason to stop)
      character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
      integer, intent(in) :: retcode
!@dbparam dump_core if set to 1 dump core for debugging
      integer :: dump_core = 0
      integer, parameter :: iu_err = 9
      !##Not used integer :: mpi_err
#ifdef USE_ESMF
      include 'mpif.h'
      integer :: rank
#endif

c**** don't call sync_param if the error is in 'PARAM' to avoid loops
      if ( message(1:6) .ne. 'PARAM:' ) then
        ! don't worry about dumping core in standalone mode
        !call sync_param( "dump_core", dump_core )
      else
        write (6,*) " Error in PARAM: Can't use sync_param."
        write (6,*) " Will use default dump_core = ", dump_core
      endif
#ifdef USE_ESMF      
      call mpi_comm_rank(MPI_COMM_WORLD, rank, mpi_err)
#else
      rank =0
#endif
      If (rank == 0) Then
        write (6,'(//2(" ",132("*")/))')
        write (6,*) ' Program terminated due to the following reason:'
        write (6,*) ' >>  ', message, '  <<'
        write (6,'(/2(" ",132("*")/))')
      End If

      if ( retcode .ne. 12 .and. retcode .ne. 13 ) then
        open( iu_err, file='error_message',
     &       form='FORMATTED', status='REPLACE', ERR=10 )
        write ( iu_err, *, ERR=10 ) message
        close( iu_err )
      endif

      call sys_flush(6)

      if ( retcode > 13 .and. dump_core > 0 ) then
#ifdef USE_ESMF
        call mpi_abort(MPI_COMM_WORLD, retcode,iu_err)
#else
        call sys_abort
#endif
      else
#ifdef USE_ESMF
        call mpi_finalize(mpi_err)
#endif
        call exit_rc ( retcode )
      endif

 10   continue
c**** something is terribly wrong, writing to stderr instead ...
      write( 0, * ) "Can't write to the error_message file"
      write( 0, * ) message
      call exit_rc ( retcode )
      end subroutine stop_model
