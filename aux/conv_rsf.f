      program convert_rsf_files
C**** convert rsf files from model II' (B399) format to modelE format 
C**** compile with: gmake conv_rsf.o
C**** f90 -o conv_rsf conv_rsf.o *.o -O2 -64 -mips4 -static \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE E001M12_COM, only : im,jm,lm,wm,u,v,t,p,q,jc,rc,clabel
     *     ,iowrite_mon 
      USE SOMTQ_COM
      USE GHYCOM, only : ghdata,snowe,tearth,wearth,aiearth,snoage
      USE RADNCB, only : rqt,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav
      USE DAGCOM, only : keynr
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
     *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
      USE OCEAN, only : tocean,oa,z1o
      USE SEAICE_COM, only : rsi,msi,tsi,snowi
      USE LANDICE_COM, only : tlandi,snowli
      USE LAKES_COM, only : t50
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60
      INTEGER IARGC,iu_AIC,I,J,L,N,ioerr
      REAL*8 TAUX,X   ! ? temporary for compatibility only
      INTEGER ItimeX

      IF (IARGC().lt.2) THEN
        PRINT*,"Convert rsf files from old format to new"
        PRINT*,"conv_rsf filename output_file"
        STOP
      END IF

      CALL GETARG(1,infile)
      CALL GETARG(2,outfile)

      iu_AIC=9
      OPEN(iu_AIC,FILE=trim(infile),FORM="UNFORMATTED",STATUS="OLD")

      READ (iu_AIC,ERR=800,END=810) TAUX,JC,CLABEL,RC,KEYNR,
     *     U,V,T,P,Q,
     2     ((TOCEAN(1,I,J),I=1,IM),J=1,JM),RSI,MSI,
     *     (((TOCEAN(L,I,J),I=1,IM),J=1,JM),L=2,3),
     *     SNOWI,SNOWE,
     *     ((TSI(1,I,J),I=1,IM),J=1,JM),TEARTH,WEARTH,AIEARTH,
     *     ((TSI(2,I,J),I=1,IM),J=1,JM),((X,I=1,IM),J=1,JM),
     *     (((SNOAGE(L,I,J),I=1,IM),J=1,JM),L=1,3),SNOWLI,
     *     (((TLANDI(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     (((TSI(L,I,J),I=1,IM),J=1,JM),L=3,4),GHDATA,
     *     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg,ustar,
     *     uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     A     (((TTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     B     (((QTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     C     (((SVLHX(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     D     (((RHSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),WM,
     E     (((CLDSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     4     ((((TMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9),
     4     ((((QMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9),
     6     (((RQT(L,I,J),I=1,IM),J=1,JM),L=1,LM_REQ),T50
      CLOSE (iu_AIC)

      ItimeX=NINT(TAUX)
      print*,ItimeX

      OPEN(iu_AIC,FILE=trim(outfile),
     *     FORM="UNFORMATTED",STATUS="UNKNOWN")

      call io_rsf(iu_AIC,ItimeX,iowrite_mon,ioerr)
      close (iu_AIC)
      print*,ioerr
      print*,"New rsf file written out to ",trim(outfile)
      stop
 800  print*,"Error reading in file"
 810  print*,"Error reading in file"
      stop
      end 
