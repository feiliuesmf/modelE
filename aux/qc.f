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
      NAMELIST/INPUTZ/ KOCEAN,PTOP,PSF,LS1,DTsrc
C****    List of parameters that COULD be changed during a run:
     *     ,DT,  NRAD,NIsurf,NFILTR, MFILTR
     *     ,SKIPSE,  NDAA,NDA5D,NDA5K,NDA5S,NDA4,NDASF
     *     ,U00wtr,U00ice,LMCM, S0X,CO2,XCDLM,IRAND
     *     ,Nslp,Kvflxo,KCOPY,Ndisk, Nssw        ,keyct ! keyct??
     *     ,KDIAG,NIPRNT,NMONAV,IJD6,NAMD6
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK
     *     ,IYEAR0,KACC0,KTACC0

      CHARACTER*80 FILEIN(5),DIR*60
      INTEGER N,Itime1,NARGS,K,iargc,KSTART,KFILE,I
      INTEGER :: ioerr=0 
      REAL*8 TOT
      LOGICAL QCALL

      QCALL=.FALSE.
      KSTART=1
      call alloc_param( "IDACC", IDACC, (/12*0/), 12) 
     
      NARGS = IARGC()
      IF(NARGS.LE.0)  GO TO 800
C**** check for argument
      CALL GETARG(1,FILEIN(1))
      IF (FILEIN(1).eq."-a") THEN
        QCALL=.TRUE.
        KSTART=2
      END IF
C**** loop over remaining arguments
      DO KFILE=KSTART,NARGS
        CALL GETARG (KFILE,FILEIN(KFILE))
      END DO
C**** Check first file
      K=KSTART
      OPEN (10,FILE=FILEIN(K),FORM='UNFORMATTED',STATUS='OLD',err=100)
      call io_label(10,Itime1,ioread,ioerr)
      CLOSE (10)
      if (ioerr.eq.0) goto 150
C**** problem reading file => assume a run directory in /u/cmrun
 100  NARGS=KSTART+1
      DIR = "/u/cmrun/"//TRIM(FILEIN(KSTART))
      FILEIN(KSTART) = TRIM(DIR)//"/fort.1"
      FILEIN(KSTART+1) = TRIM(DIR)//"/fort.2"
        
 150  DO K=KSTART,NARGS
      OPEN (10,FILE=FILEIN(K),FORM='UNFORMATTED',STATUS='OLD',err=850)
      ioerr=0
      call io_label(10,Itime1,ioread,ioerr)
      if (ioerr.eq.1) go to 860
      CLOSE (10)
      
      call get_param( "Itime", Itime) 
      call get_param( "IYEAR0", IYEAR0) 
      call get_param( "NDAY", NDAY )
      call getdte(Itime,Nday,Iyear0,Jyear,Jmon,Jday,Jdate,Jhour,amon)

      WRITE (6,900) ITIME,JMON,JDATE,JYEAR,JHOUR,XLABEL(1:50)
      TOT=0
      DO N=1,NTIMEACC
        TOT = TOT + TIMING(N) 
      END DO
      TOT = TOT/100.
      IF (TOT.gt.0) THEN
      WRITE (6,906) 24.*TOT/(60.*(Itime1-Itime0)),(TIMESTR(N),TIMING(N)
     *     /TOT,N=1,3) 
      WRITE (6,907) (TIMESTR(N),TIMING(N)/TOT,N=4,NTIMEACC)
      END IF
      IF (QCALL) THEN
        write (6,*) "IM,JM,LM = ",IM,JM,LM
        write (6,*) "IDACC = ",(IDACC(I),I=1,12)
        write (6,inputz)
        write (6,*) 
      END IF
      END DO
      Stop
  800 continue
      write(6,*) " Usage: qc [-a] FILE1 [FILE2...]"
      write(6,*) " Where FILE1,2..(up to 5) are modelE files ",
     *     "(rsf/acc/fort.[12])"
      write(6,*) "  -a  output all header information, otherwise"
      write(6,*) "      only timing info is printed"
      write(6,*) " If first FILE argument is a directory, then ",
     *     "DIR/fort.[12] are opened"
      STOP
  850 continue
      write(6,*) "Cannot open file ",FILEIN(K)
      STOP
  860 continue
      write(6,*) "Error reading file ", FILEIN(K)
      STOP
 900  FORMAT (I10,1X,I2,'/',I2,'/',I4,' hr',I2,2X,A50)
 906  FORMAT (' TIME',F7.2,' (MINUTES)',3(X,A12,F5.1))
 907  FORMAT (10(22X,3(X,A12,F5.1) / ))
      end

