C**** CMPE001.F    CoMPare restartfiles for modelE          6/00
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IM=72,JM=46,LM=12,NGM=6,KTD=8,KAIJ=100,KAJK=51)
      parameter (npbl=8)
      PARAMETER (IMH=IM/2,  KACC=JM*94*3 + 24*94 + JM*3 +
     *   JM*LM*57 + JM*3*4 + IM*JM*KAIJ + IM*LM*16 +
     *   IM*JM*29 + 20*100 + JM*36 + (IMH+1)*20*8 +
     *   8*2 + 24*63*4 + 2*62*10*12 + JM*LM*KAJK +
     *   IM*JM*LM*6 + IM*JM*LM*5 + JM*LM*10*3)
      COMMON/TADV/TMOM1(IM,JM,9*LM),TMOM2(IM,JM,9*LM)
      COMMON/QADV/QMOM1(IM,JM,9*LM),QMOM2(IM,JM,9*LM)
      COMMON/BNDYCB1/ ODATA1(IM,JM,5),GDATA1(IM,JM,16),BLDATA1(IM,JM,12)
      COMMON/BNDYCB2/ ODATA2(IM,JM,5),GDATA2(IM,JM,16),BLDATA2(IM,JM,12)
      COMMON/RADNCB/ RADN1(IM,JM,9+2*LM),RADN2(IM,JM,9+2*LM)
      COMMON/WORKO/  OA(IM,JM,12)
      DIMENSION U1(IM,JM,LM),V1(IM,JM,LM),T1(IM,JM,LM),P1(IM,JM),Q1(IM
     *     ,JM,LM)
      DIMENSION U2(IM,JM,LM),V2(IM,JM,LM),T2(IM,JM,LM),P2(IM,JM),Q2(IM
     *     ,JM,LM)
      common /socabl1/pbl1(npbl,im,jm,5*4),pblb1(im,jm,3*4),ipbl1(im,jm
     *     ,4)
      common /socabl2/pbl2(npbl,im,jm,5*4),pblb2(im,jm,3*4),ipbl2(im,jm
     *     ,4)
      COMMON/CLDCOM/CLOUD1(IM,JM,6*LM),CLOUD2(IM,JM,6*LM)
      COMMON/SOILS3/GHDATA1(IM,JM,4*NGM+5),GHDATA2(IM,JM,4*NGM+5)
      CHARACTER C*4,XLABEL*4,LABEL*16,FILEIN*60
      DIMENSION C(39),JC(100),RC(161)
      EQUIVALENCE (C,XLABEL,LABEL)
      REAL*8 DIAG(KACC),TSFREZ(IM,JM,2),TDIURN(IM,JM,KTD)
      COMMON /KEYS/ KEYNR(42,50)
C****
      IF(IARGC().NE.2)  GO TO 800
C****
C**** Read ReStartFiles
C****
      CALL GETARG (1,FILEIN)
      OPEN (1,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
      READ (1) TAU1,JC,C,RC,KEYNR,U1,V1,T1,P1,Q1,ODATA1,GDATA1,
     *  GHDATA1,BLDATA1,PBL1,pblb1,ipbl1,CLOUD1,TMOM1,QMOM1,RADN1,
     *  TSFREZ,DIAG,TDIURN,OA,TAU2
      IF (TAU1.ne.TAU2) then
         WRITE (6,*) 'FILE 1 NOT READ CORRECTLY. IHOUR,IHOURB =',tau1
     *        ,tau2
         STOP
      END IF
      CLOSE (1)

      WRITE (6,*) 'Unit 1 read.  IHOUR1 =',tau1,'   ',LABEL
C****
      CALL GETARG (2,FILEIN)
      OPEN (2,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
      READ (2) TAU1,JC,C,RC,KEYNR,U2,V2,T2,P2,Q2,ODATA2,GDATA2,
     *  GHDATA2,BLDATA2,PBL2,pblb2,ipbl2,CLOUD2,TMOM2,QMOM2,RADN2,
     *  TSFREZ,DIAG,TDIURN,OA,TAU2
      IF (tau1.ne.tau2) then
         WRITE (6,*) 'FILE 2 NOT READ CORRECTLY. IHOUR,IHOURB =',tau1
     *        ,tau2
         STOP
      END IF
      CLOSE (2)

      WRITE (6,*) 'Unit 2 read.  IHOUR2 =',tau1,'   ',LABEL
C****
C**** Compare arrays
C****
      WRITE (6,*) 'FIELD IM JM LM     VAL1        VAL2     RELERR'
      CALL COMP8 ('U     ',IM,JM,LM     ,U1 ,U2 )
      CALL COMP8 ('V     ',IM,JM,LM     ,V1 ,V2 )
      CALL COMP8 ('T     ',IM,JM,LM     ,T1 ,T2 )
      CALL COMP8 ('Q     ',IM,JM,LM     ,Q1 ,Q2 )
      CALL COMP8 ('P     ',IM,JM,1      ,P1 ,P2 )

      CALL COMP8 ('ODATA ',IM,JM,5      ,ODATA1 ,ODATA2 )
      CALL COMP8 ('GDATA ',IM,JM,16     ,GDATA1 ,GDATA2 )
      CALL COMP8 ('GHDATA',IM,JM,4*NGM+5,GHDATA1,GHDATA2)
      CALL COMP8 ('BLDATA',IM,JM,12     ,BLDATA1,BLDATA2)

      CALL COMP8 ('PBL   ',npbl,IM*JM,5*4 ,PBL1  ,PBL2  )
      CALL COMP8 ('PBLB  ',IM,JM,3*4      ,PBLB1 ,PBLB2 )
      CALL COMP8 ('CLOUD ',IM,JM,6*LM   ,CLOUD1,CLOUD2)
      CALL COMP8 ('TMOM  ',IM,JM,9*LM   ,TMOM1 ,TMOM2 )
      CALL COMP8 ('QMOM  ',IM,JM,9*LM   ,QMOM1 ,QMOM2 )
      CALL COMP8 ('RADN  ',IM,JM,9+2*LM ,RADN1 ,RADN2 )
      STOP
C****
  800 WRITE (0,*) 'Example: CMPE001 E001.rsf_1 E001.rsf_2   ',
     *            'CoMPare restartfiles     7/18/94'
      STOP
  810 WRITE (0,*) 'Error accessing ',FILEIN
      STOP
      END

      SUBROUTINE COMP4 (LABEL,IM,JM,LM,X1,X2)
C****
C**** COMP compares two arrays and prints any discrepancies to the
C**** line printer.
C****
      REAL*4 X1(IM,JM,LM),X2(IM,JM,LM),X1MAX,X2MAX
      CHARACTER*6 LABEL
      NP = 0
      DIFMAX = 0.
      IMAX=0
      JMAX=0
      LMAX=0
      X1MAX=0.
      X2MAX=0.
      WRITE (6,900) LABEL
      DO 10 L=1,LM
      DO 10 J=1,JM
      DO 10 I=1,IM
         IF (X2(I,J,L)*X1(I,J,L).EQ.0 .AND. X2(I,J,L)+X1(I,J,L).NE.0
     *        .and. X2(I,J,L).NE.X1(I,J,L)) THEN
         WRITE (6,901) I,J,L,X1(I,J,L),X2(I,J,L)
         NP = NP+1
         IF(NP.GE.10)  RETURN
      ELSE
         DIF = ABS(X2(I,J,L)-X1(I,J,L)) / (ABS(X1(I,J,L)) + 1.E-30)
         IF(DIF.LE.DIFMAX)  GO TO 10
         DIFMAX = DIF
         IMAX=I
         JMAX=J
         LMAX=L
         X1MAX=X1(I,J,L)
         X2MAX=X2(I,J,L)
c         WRITE (6,901) I,J,L,X1(I,J,L),X2(I,J,L),DIF
c         NP = NP+1
c         IF(NP.GE.10)  RETURN
      END IF
   10 CONTINUE
      IF (IMAX.ne.0) WRITE (6,901) IMAX,JMAX,LMAX,X1MAX,X2MAX,DIFMAX
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I6,3E18.8)
      END

      SUBROUTINE COMP8 (LABEL,IM,JM,LM,X1,X2)
C****
C**** COMP8 compares two arrays and prints any discrepancies to the
C**** line printer.
C****
      REAL*8 X1(IM,JM,LM),X2(IM,JM,LM), DIF,DIFMAX,X2MAX,X1MAX
      CHARACTER*6 LABEL
      NP = 0
      WRITE (6,900) LABEL
      DO 20 L=1,LM
      DIFMAX = 0.
      IMAX=0
      JMAX=0
      LMAX=0
      X1MAX=0.
      X2MAX=0.
      DO 10 J=1,JM
      DO 10 I=1,IM
      IF (X2(I,J,L)*X1(I,J,L).EQ.0 .AND. X2(I,J,L)+X1(I,J,L).NE.0 
     *        .and. X2(I,J,L).NE.X1(I,J,L)) THEN
         WRITE (6,901) I,J,L,X1(I,J,L),X2(I,J,L)
         NP = NP+1
         IF(NP.GE.10)  RETURN
      ELSE
         DIF = ABS(X2(I,J,L)-X1(I,J,L)) / (ABS(X1(I,J,L)) + 1.D-30)
c         IF (DIF.NE.0.) PRINT*,I,J,L,DIF
         IF(DIF.LE.DIFMAX)  GO TO 10
         DIFMAX = DIF
         IMAX=I
         JMAX=J
         LMAX=L
         X1MAX=X1(I,J,L)
         X2MAX=X2(I,J,L)
c      WRITE (6,901) I,J,L,X1(I,J,L),X2(I,J,L),DIF
c      NP = NP+1
c      IF(NP.GE.10)  RETURN
      END IF
 10   CONTINUE
      IF (IMAX.ne.0) WRITE (6,901) IMAX,JMAX,LMAX,X1MAX,X2MAX,DIFMAX
 20   CONTINUE
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I4,E30.20,2E30.20)
      END

