      program qc
!@sum  qc query c-array and other parameters
!@auth G. Schmidt/R. Ruedy
!@ver  1.0
      USE MODEL_COM
      USE TIMINGS
      USE PARAM
      IMPLICIT NONE
      CHARACTER*80 FILEIN
      INTEGER N,NARGS,K,iargc,KFILE,I,days_togo
      INTEGER :: ioerr=0, KSTART=1
      REAL*8 TOT,yrs_togo,FAC,FACT
      LOGICAL :: QCALL = .FALSE., QCMIN=.FALSE.

      NARGS = IARGC()
      IF(NARGS.LE.0)  GO TO 800
C**** check for arguments
 10   CALL GETARG(KSTART,FILEIN)
      IF (FILEIN(1:1).eq."-") THEN
        SELECT CASE (FILEIN(2:2))
        CASE ("a","A")
          QCALL=.TRUE.
        CASE ("t","T")
          QCMIN=.TRUE.
        END SELECT
        KSTART=KSTART+1
        GOTO 10
      END IF

      DO K=KSTART,NARGS
      if (qcall .and. k.gt.kstart) write (6,*)
      CALL GETARG (K,FILEIN)
      OPEN (10,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',err=850)
      ioerr=0
      call io_label(10,Itime,ioread,ioerr)
      if (ioerr.eq.1) go to 860
      CLOSE (10)

      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

      WRITE (6,900) ITIME,JMON,JDATE,JYEAR,JHOUR,XLABEL(1:50)
      TOT=0
      DO N=1,NTIMEACC
        TOT = TOT + TIMING(N)
      END DO
      IF (Itime-Itime0.gt.0) THEN
        FACT = 24./(6000*(Itime-Itime0))
      ELSE
        FACT = 0.
      END IF
      IF (QCMIN) THEN ! output in minutes
        FAC = FACT
      ELSE            ! output in percentages
        FAC = 100./TOT
      END IF
      IF (TOT.gt.0) THEN
        WRITE (6,906) FACT*TOT,(TIMESTR(N),FAC*TIMING(N),N=1,3)
        WRITE (6,907) (TIMESTR(N),FAC*TIMING(N),N=4,NTIMEACC)
      END IF
      IF (QCALL) THEN
        call print_param(6)
c       write (6,*) "IDACC = ",(IDACC(I),I=1,12)
        call getdte(ItimeI,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
        WRITE (6,900) ITIMEI,JMON,JDATE,JYEAR,JHOUR,' = start of run'
        call getdte(ItimeE,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
        WRITE (6,900) ITIMEE,JMON,JDATE,JYEAR,JHOUR,' =   end of run'
        if(itimee.gt.itime) then
          days_togo = (Itimee-itime+nday-1)/nday
          yrs_togo  = (Itimee-itime)/(365.*nday)
          write(XLABEL(29:50),'(I10,a12)') days_togo,'  days to go'
          if (days_togo.eq.1) XLABEL(44:44) = ' '
         call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
          if (yrs_togo.ge.2.)
     *    write(XLABEL(29:50),'(f10.1,a12)') yrs_togo,' years to go'
        WRITE (6,900) ITIME,JMON,JDATE,JYEAR,JHOUR,XLABEL(1:50)
      end if
      END IF
      END DO
      Stop
  800 continue
      write(6,*) " Usage: qc [-a] [-t] FILE1 [FILE2...]"
      write(6,*) " Where FILE1,2.. are modelE files ",
     *     "(rsf/acc/fort.[12])"
      write(6,*) "  -a  output all header information, otherwise"
      write(6,*) "      only timing info is printed"
      write(6,*) "  -t  output absolute minutes for each sub-section,"
      write(6,*) "      otherwise percentage time is printed"
      STOP
  850 continue
      write(6,*) "Cannot open file ",FILEIN
      STOP
  860 continue
      write(6,*) "Error reading file ", FILEIN
      STOP
 900  FORMAT (I10,1X,I2,'/',I2.2,'/',I4.4,' hr',I2,2X,A)
 906  FORMAT (' TIME',F7.2,' (MINUTES)',3(X,A12,F5.1))
 907  FORMAT (10(22X,3(X,A12,F5.1) / ))
      end

