      program qc
!@sum  qc query c-array and other parameters
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM
      USE CLD01_COM_E001
      USE RADNCB
      USE DAGCOM
      IMPLICIT NONE
      NAMELIST/INPUTZ/ KOCEAN,PTOP,PSF,LS1,DTsrc
C****    List of parameters that COULD be changed during a run:
     *     ,DT,  NRAD,NIsurf,NFILTR, MFILTR
     *     ,SKIPSE,  NDAA,NDA5D,NDA5K,NDA5S,NDA4,NDASF
     *     ,U00wtr,U00ice,LMCM, S0X,CO2,XCDLM,IRAND
     *     ,Nslp,Kvflxo,KCOPY,Ndisk, Nssw        ,keyct ! keyct??
     *     ,KDIAG,NIPRNT,NMONAV,IJD6,NAMD6
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK
     *     ,IYEAR0,mdyn,mcnds,mrad,msurf,mdiag,melse
     *     ,KACC0,KTACC0

      CHARACTER*80 FILEIN
      INTEGER Itime1,NARGS,K,ioerr,iargc
      REAL*8 TOT

      NARGS = IARGC()
      IF(NARGS.LE.0)  GO TO 800
      CALL GETARG (1,FILEIN)
      OPEN (10,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',err=850)
      call io_rsf(10,Itime1,ioread,ioerr)
      WRITE (6,900) ITIME,JMON,JDATE,JYEAR,JHOUR,XLABEL(1:50)
      TOT=(MDYN+MCNDS+MRAD+MSURF+MDIAG+MELSE)/100.
      WRITE (6,906) 24.*TOT/(60.*(Itime1-Itime0)),MDYN/TOT,MCNDS/TOT
     *     ,MRAD/TOT,MSURF/TOT,MDIAG/TOT,MELSE/TOT
      write (6,*) "IM,JM,LM = ",IM,JM,LM
      write (6,inputz)
      Stop
  800 continue
      write(6,'(" Usage: qc FILENAME(rsf-file modelE)")')
      Stop
  850 continue
      write(6,'(" cannot open file",a)') FILEIN
      Stop
 900  FORMAT (I10,1X,I2,'/',I2,'/',I4,' hr',I2,2X,A50)
 906  FORMAT (' TIME',F7.2,'(MINUTES)    DYNAMICS',F5.1,
     *  '    CONDENSATION',F5.1,'    RADIATION',F5.1/25X,'SURFACE ',
     *  F5.1,'    DIAGNOSTICS ',F5.1,'    OTHER',F9.1)
      end
