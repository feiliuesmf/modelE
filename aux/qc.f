      program qc
!@sum  qc query c-array and other parameters
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM
      USE CLD01_COM_E001
      USE RADNCB
      USE DAGCOM
      USE TIMINGS
      USE PARAM
      IMPLICIT NONE
      CHARACTER*80 FILEIN,DIR*60
      INTEGER N,NARGS,K,iargc,KSTART,KFILE,I,days_togo
      INTEGER :: ioerr=0
      REAL*8 TOT,yrs_togo
      LOGICAL QCALL

      QCALL=.FALSE.
      KSTART=1

      NARGS = IARGC()
      IF(NARGS.LE.0)  GO TO 800
C**** check for argument
      CALL GETARG(1,FILEIN)
      IF (FILEIN(1:1).eq."-") THEN
        QCALL=.TRUE.
        KSTART=2
      END IF
C**** loop over remaining arguments
c     DO KFILE=KSTART,NARGS
c       CALL GETARG (KFILE,FILEIN(KFILE))
c     END DO
C**** Check first file
c     K=KSTART
c     OPEN (10,FILE=FILEIN(K),FORM='UNFORMATTED',STATUS='OLD',err=100)
c     call io_label(10,Itime,ioread,ioerr)
c     CLOSE (10)
c     if (ioerr.eq.0) goto 150
C**** problem reading file => assume a run directory in /u/cmrun
c100  NARGS=KSTART+1
c     DIR = "/u/cmrun/"//TRIM(FILEIN(KSTART))
c     FILEIN(KSTART) = TRIM(DIR)//"/fort.1"
c     FILEIN(KSTART+1) = TRIM(DIR)//"/fort.2"

 150  DO K=KSTART,NARGS
      if (qcall .and. k.gt.kstart) write (6,*)
      CALL GETARG (K,FILEIN)
      OPEN (10,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',err=850)
      ioerr=0
      call io_label(10,Itime,ioread,ioerr)
      if (ioerr.eq.1) go to 860
      CLOSE (10)

c     call get_param( "Itime", Itime)
c     call get_param( "IYEAR0", IYEAR0)
c     call get_param( "NDAY", NDAY )
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

      WRITE (6,900) ITIME,JMON,JDATE,JYEAR,JHOUR,XLABEL(1:50)
      TOT=0
      DO N=1,NTIMEACC
        TOT = TOT + TIMING(N)
      END DO
      TOT = TOT/100.
      IF (TOT.gt.0) THEN
      WRITE (6,906) 24.*TOT/(60.*(Itime-Itime0)),(TIMESTR(N),TIMING(N)
     *     /TOT,N=1,3)
      WRITE (6,907) (TIMESTR(N),TIMING(N)/TOT,N=4,NTIMEACC)
      END IF
      IF (QCALL) THEN
        write (6,*) "IM,JM,LM = ",IM,JM,LM," LS1 = ",LS1," P_cp = ",PTOP
        write (6,*) "PLbot(mean values) = ",PTOP+PSFMPT*SIGE
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
          if (yrs_togo.ge.2.)
     *    write(XLABEL(29:50),'(f10.1,a12)') yrs_togo,' years to go'
        WRITE (6,900) ITIME,JMON,JDATE,JYEAR,JHOUR,XLABEL(1:50)
      end if
      END IF
      END DO
      Stop
  800 continue
      write(6,*) " Usage: qc [-a] FILE1 [FILE2...]"
      write(6,*) " Where FILE1,2.. are modelE files ",
     *     "(rsf/acc/fort.[12])"
      write(6,*) "  -a  output all header information, otherwise"
      write(6,*) "      only timing info is printed"
c     write(6,*) " If first FILE argument is a directory, then ",
c    *     "DIR/fort.[12] are opened"
      STOP
  850 continue
      write(6,*) "Cannot open file ",FILEIN
      STOP
  860 continue
      write(6,*) "Error reading file ", FILEIN
      STOP
 900  FORMAT (I10,1X,I2,'/',I2,'/',I4,' hr',I2,2X,A)
 906  FORMAT (' TIME',F7.2,' (MINUTES)',3(X,A12,F5.1))
 907  FORMAT (10(22X,3(X,A12,F5.1) / ))
      end

