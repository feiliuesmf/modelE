      program qc
!@sum  qc query c-array and other parameters
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM
      USE CLD01_COM_E001
      USE RADNCB
      USE DAGCOM
      USE TIMINGS
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

      CHARACTER*80 FILEIN
      INTEGER N,Itime1,NARGS,K,ioerr,iargc,KSTART,KFILE
      REAL*8 TOT
      LOGICAL QCALL

      QCALL=.FALSE.
      KSTART=1
      NARGS = IARGC()
      IF(NARGS.LE.0)  GO TO 800
C**** check for argument
      CALL GETARG(1,FILEIN)
      IF (FILEIN.eq."-a") THEN
        QCALL=.TRUE.
        KSTART=2
      END IF
C**** loop over remaining arguments (assumed to be files)
      DO KFILE=KSTART,NARGS
      CALL GETARG (KFILE,FILEIN)
      OPEN (10,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',err=850)
      call io_label(10,Itime1,ioread,ioerr)
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
        write (6,inputz)
        write (6,*) 
      END IF
      END DO
      Stop
  800 continue
      write(6,*) " Usage: qc [-a] FILE1 [FILE2...]"
      write(6,*) " Where FILE1,2 are modelE files (rsf/acc/fort.[12])"
      write(6,*) "  -a  output all header information, otherwise"
      write(6,*) "      only timing info is printed"
      Stop
  850 continue
      write(6,'(" cannot open file",a)') FILEIN
      Stop
 900  FORMAT (I10,1X,I2,'/',I2,'/',I4,' hr',I2,2X,A50)
 906  FORMAT (' TIME',F7.2,' (MINUTES)',3(X,A12,F5.1))
 907  FORMAT (10(22X,3(X,A12,F5.1) / ))
      end

