      PROGRAM mkdeep
!@sum  mkdeep creates the mixed layer temperature climatology needed of
!@+    deep ocean diffusion runs
!@auth G. Russell/L. Nazarenko/G. Schmidt
C**** Note this must be compiled using the deep ocean gcmlib
      USE MODEL_COM, only : ioread,Itime,im,jm,amonth
      USE OCEAN, only : tocean
      USE ODEEP_COM, only : lmom
      USE SEAICE_COM, only : rsi,msi,hsi,ssi
      USE SEAICE, only : lmi
      USE DAGCOM, only : ij_tgo2,aij
      USE FILEMANAGER, only : closeunit,openunit
      IMPLICIT NONE
      CHARACTER TITLE*80,FNAME*60,RUNID*20
      REAL*8 TG3M(IM,JM,12),TOCEANS(3,IM,JM),RSIS(IM,JM),MSIS(IM,JM)
     *     ,HSIS(LMI,IM,JM),SSIS(LMI,IM,JM)
      INTEGER M,MONACC(12),K,MTG3,ioerr,IYR,IYRE,IYRI,NFILES,IARGC
     *     ,iu_oda,lf

C**** ZERO OUT ACCUMULATING ARRAYS
      TOCEANS=0. ; RSIS=0. ; MSIS=0. ; HSIS=0. ; SSIS=0. ; TG3M=0.
C**** Initiallise AIJ array (needed due to the way io_oda works)
      IJ_TGO2 = 1
      AIJ(:,:,1) = 0.

C**** LOOP OVER INPUT FILES 
      NFILES=IARGC()
      IF (NFILES.le.0) THEN
        PRINT*,"mkdeep: make climatology to initialise deep ocean runs"
        PRINT*,"Usage: mkdeep *.oda.*"
        STOP
      END IF

      IYRE=0
      IYRI=9999
      monacc=0

C**** loop over input files
      DO K=1,NFILES
        CALL GETARG(K,FNAME)
C**** FNAME = monYYYY.oda.RUNID
C**** get relevant month 
        DO M=1,12
          IF (FNAME(1:3).eq.AMONTH(M)(1:3)) THEN
            monacc(M)=monacc(M)+1
            MTG3=M
          END IF
        END DO
C**** years
        READ(FNAME(4:7),'(I4)') IYR
        IYRE=MAX(IYRE,IYR)
        IYRI=MIN(IYRI,IYR)

        call openunit(FNAME,iu_oda,.true.,.true.)
        call io_oda(iu_oda,itime,ioread,ioerr)
        if (ioerr.eq.1) goto 910
        call closeunit(iu_oda)
          
C**** ACCUMULATE TG3M SEPARATELY FOR EACH MONTH
        TG3M(:,:,MTG3)=TG3M(:,:,MTG3)+AIJ(:,:,IJ_TGO2)

C**** ACCUMULATE OCEAN DATA AT DEC 1 (from NOV file)
        IF (MTG3.eq.11) THEN
          TOCEANS=TOCEANS+TOCEAN
          RSIS=RSIS+RSI
          MSIS=MSIS+MSI
          HSIS=HSIS+HSI
          SSIS=SSIS+SSI
        END IF
      END DO
      RUNID=FNAME(13:len_trim(FNAME))

C**** DIVIDE EACH MONTH BY NUMBER OF ACCUMULATIONS
      DO M=1,12
        TG3M(:,:,M)=TG3M(:,:,M)/monacc(M)
      END DO
      TOCEANS=TOCEANS/monacc(11)
      RSIS=RSIS/monacc(11)
      MSIS=MSIS/monacc(11)
      HSIS=HSIS/monacc(11)
      SSIS=SSIS/monacc(11)
      
      TITLE="TOCEAN,RSI,MSI,HSI,SSI mean initial conditions Dec.1, TG3M"
      FNAME="TG3M."//RUNID//"."
      lf=len_trim(FNAME)
      WRITE(FNAME(lf:lf+7),'(2I4)') IYRI,IYRE

      call openunit(FNAME,iu_oda,.true.,.false.)
      WRITE(iu_oda) TITLE,TOCEANS,RSIS,MSIS,HSIS,SSIS,TG3M
      call closeunit(iu_oda)
      
      STOP
  910 WRITE(6,*) 'READ ERROR ON FILE ',FNAME
      STOP
      END




