      program convert_rsf_files
C**** convert rsf files from model II' (B399) format to modelE format 
C**** compile with: gmake conv_rsf.o
C**** f90 -o conv_rsf conv_rsf.o *.o -O2 -64 -mips4 -static \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE E001M12_COM, only : im,jm,lm,wm,u,v,t,p,q,gdata,jc,rc,clabel
     *     ,iowrite_mon 
      USE SOMTQ_COM
      USE GHYCOM, only : ghdata
      USE RADNCB, only : rqt,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav
      USE DAGCOM, only : keynr
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
     *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
      USE OCEAN, only : tocean,oa,z1o
      USE SEAICE_COM, only : rsi,msi
      USE LAKES_COM, only : t50
      IMPLICIT NONE
      CHARACTER infile*60
      INTEGER IARGC,iu_AIC,I,J,L,ioerr
      REAL*8 TAUX   ! ? temporary for compatibility only
      INTEGER ItimeX

      IF (IARGC().eq.0) THEN
        PRINT*,"Convert rsf files from old format to new"
        PRINT*,"conv_rsf filename"
        STOP
      END IF

      CALL GETARG(1,infile)
      iu_AIC=9
      OPEN(iu_AIC,FILE=trim(infile),FORM="UNFORMATTED",STATUS="OLD")

      READ (iu_AIC,ERR=800,END=810) TAUX,JC,CLABEL,RC,KEYNR,
     *     U,V,T,P,Q,
     2     ((TOCEAN(1,I,J),I=1,IM),J=1,JM),RSI,MSI,
     *     (((TOCEAN(L,I,J),I=1,IM),J=1,JM),L=2,3),
     *     GDATA,GHDATA,
     *     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg,ustar,
     *     uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     A     (((TTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     B     (((QTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     C     (((SVLHX(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     D     (((RHSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),WM,
     E     (((CLDSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     4     TX,TY,TZ,TXX,TYY,TZZ,TXY,TZX,TYZ,
     5     QX,QY,QZ,QXX,QYY,QZZ,QXY,QZX,QYZ,
     6     (((RQT(L,I,J),I=1,IM),J=1,JM),L=1,LM_REQ),T50
      CLOSE (iu_AIC)

      ItimeX=NINT(TAUX)
      print*,ItimeX

      OPEN(iu_AIC,FILE=trim(infile)//".modelE",
     *     FORM="UNFORMATTED",STATUS="UNKNOWN")

      call io_rsf(iu_AIC,ItimeX,iowrite_mon,ioerr)
      close (iu_AIC)
      print*,ioerr
      print*,"New rsf file written out to ",trim(infile)//".modelE"
      stop
 800  print*,"Error reading in file"
 810  print*,"Error reading in file"
      stop
      end 







