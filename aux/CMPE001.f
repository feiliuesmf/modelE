C**** CMPE001.F    CoMPare restartfiles for modelE          6/00
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IM=72,JM=46,LM=12,NGM=6,KTD=8,KAIJ=150,KAJK=51,KCON=125
     *     )
      parameter (npbl=8)
      PARAMETER (IMH=IM/2, NWD=MIN(IMH,9),
     *     KACC=JM*94*6 + 24*94 +
     *     JM*3 +JM*LM*57 + JM*3*4 + IM*JM*KAIJ + IM*LM*16 +
     *     20*100 + JM*KCON + (IMH+1)*20*8 +8*2 + 24*63*4 + 2*62*NWD*12
     *     + JM*LM*KAJK +IM*JM*LM*6 + IM*JM*LM*5 + JM*LM*(NWD+1)*3)
      COMMON/TADV/TMOM1(IM,JM,9*LM),TMOM2(IM,JM,9*LM)
      COMMON/QADV/QMOM1(IM,JM,9*LM),QMOM2(IM,JM,9*LM)
      COMMON/BNDYCB1/TOCEAN1(3,IM,JM),GDATA1(IM,JM,7),
     *     BLDATA1(IM,JM,11+LM),Z1(IM,JM)
      COMMON/BNDYCB2/TOCEAN2(3,IM,JM),GDATA2(IM,JM,7),
     *     BLDATA2(IM,JM,11+LM),Z2(IM,JM)
ccc   ice data:
      INTEGER, PARAMETER  :: LMI=4
      REAL*8 MSI1,MSI2
      COMMON/SICECB/ RSI1(IM,JM),TSI1(LMI,IM,JM),MSI1(IM,JM),SNOWI1(IM
     *     ,JM),RSI2(IM,JM),TSI2(LMI,IM,JM),MSI2(IM,JM),SNOWI2(IM,JM)
      COMMON/RADNCB/ RQT1( 3,IM,JM),RQT2( 3,IM,JM),
     *               SRHR1(1+LM,IM,JM),SRHR2(1+LM,IM,JM),
     *               TRHR1(1+LM,IM,JM),TRHR2(1+LM,IM,JM),
     *               FSF1(   4,IM,JM),FSF2(   4,IM,JM)
ccc   snow data:
      INTEGER, PARAMETER  :: NLSN=3
      INTEGER, DIMENSION(2,IM,JM)     :: NSN1, NSN2
      INTEGER, DIMENSION(2,IM,JM)     :: ISN1, ISN2
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: DZSN1, DZSN2
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: WSN1, WSN2
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: HSN1, HSN2
      REAL*8, DIMENSION(2,IM,JM)      :: FR_SNOW1, FR_SNOW2
ccc   land ice data
      REAL*8 LANDI1,LANDI2
      COMMON/LANDICB/ LANDI1(IM,JM,3),LANDI2(IM,JM,3)
      COMMON/WORKO/  OA1(IM,JM,12),OA2(IM,JM,12)
      DIMENSION U1(IM,JM,LM),V1(IM,JM,LM),T1(IM,JM,LM),P1(IM,JM),Q1(IM
     *     ,JM,LM),WM1(IM,JM,LM)
      DIMENSION U2(IM,JM,LM),V2(IM,JM,LM),T2(IM,JM,LM),P2(IM,JM),Q2(IM
     *     ,JM,LM),WM2(IM,JM,LM)
      common /socabl1/pbl1(npbl,im,jm,5*4),pblb1(im,jm,3*4),ipbl1(im,jm
     *     ,4)
      common /socabl2/pbl2(npbl,im,jm,5*4),pblb2(im,jm,3*4),ipbl2(im,jm
     *     ,4)
      COMMON/CLDCOM/CLOUD2(IM,JM,5*LM),CLOUD1(IM,JM,5*LM)
      COMMON/SOILS3/GHDATA1(IM,JM,4*NGM+5),GHDATA2(IM,JM,4*NGM+5)
      CHARACTER C*4,XLABEL*132,LABEL*16,FILEIN*60,HEADER*8
      DIMENSION C(39),JC(100),RC(161)
      EQUIVALENCE (C,XLABEL,LABEL)
      COMMON/DAG1/DIAG1(KACC),TSFREZ1(IM,JM,4),TDIURN1(IM,JM,KTD)
      COMMON/DAG2/DIAG2(KACC),TSFREZ2(IM,JM,4),TDIURN2(IM,JM,KTD)
      COMMON /KEYS/ KEYNR(42,50)
      REAL*8 LAKE1,LAKE2
      COMMON /LAKE/LAKE1(IM,JM,4),LAKE2(IM,JM,4)
      INTEGER DAGPOS,DAGPOS1,DAGPOS2
      LOGICAL ERRQ,COMP8,COMP8p,COMPI,COMP8LIJp,COMPILIJ
      INTEGER itau1,itau2
C****
      IF(IARGC().NE.2)  GO TO 800
C****
C**** Read ReStartFiles
C****
      CALL GETARG (1,FILEIN)
      OPEN (1,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
         READ (1) ITAU1,JC,C,RC
         READ (1)
         READ (1) HEADER,U1,V1,T1,P1,Q1,WM1
         READ (1) HEADER,TOCEAN1,Z1
         READ (1) HEADER,LAKE1
         READ (1) HEADER,RSI1,TSI1,SNOWI1,MSI1
         READ (1) HEADER,GDATA1
         READ (1) HEADER,GHDATA1
         READ (1) HEADER,NSN1,ISN1,DZSN1,WSN1,HSN1,FR_SNOW1
         READ (1) HEADER,LANDI1
         READ (1) HEADER,BLDATA1
         READ (1) HEADER,PBL1,pblb1,ipbl1
         READ (1) HEADER,U00wtr,U00ice,LMCM,CLOUD1
         READ (1) HEADER,TMOM1,QMOM1
         READ (1) HEADER,S0,S0X,CO2,RSD,SIND,COSD,RQT1,SRHR1,TRHR1,FSF1
         READ (1) HEADER,KEYNR,TSFREZ1,DIAG1,TDIURN1,OA1,ITAU2

      IF (ITAU1.ne.ITAU2) then
         WRITE (6,*) 'FILE 1 NOT READ CORRECTLY. IHOUR,IHOURB =',itau1
     *        ,itau2
         STOP
      END IF
      CLOSE (1)

      WRITE (6,*) 'Unit 1 read.  IHOUR1 =',itau1,'   ',LABEL
C****
      CALL GETARG (2,FILEIN)
      OPEN (2,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
         READ (2) ITAU1,JC,C,RC
         READ (2)
         READ (2) HEADER,U2,V2,T2,P2,Q2,WM2
         READ (2) HEADER,TOCEAN2,Z2
         READ (2) HEADER,LAKE2
         READ (2) HEADER,RSI2,TSI2,SNOWI2,MSI2
         READ (2) HEADER,GDATA2
         READ (2) HEADER,GHDATA2
         READ (2) HEADER,NSN2,ISN2,DZSN2,WSN2,HSN2,FR_SNOW2
         READ (2) HEADER,LANDI2
         READ (2) HEADER,BLDATA2
         READ (2) HEADER,PBL2,pblb2,ipbl2
         READ (2) HEADER,U00wtr,U00ice,LMCM,CLOUD2
         READ (2) HEADER,TMOM2,QMOM2
         READ (2) HEADER,S0,S0X,CO2,RSD,SIND,COSD,RQT2,SRHR2,TRHR2,FSF2
         READ (2) HEADER,KEYNR,TSFREZ2,DIAG2,TDIURN2,OA2,ITAU2

      IF (itau1.ne.itau2) then
         WRITE (6,*) 'FILE 2 NOT READ CORRECTLY. IHOUR,IHOURB =',itau1
     *        ,itau2
         STOP
      END IF
      CLOSE (2)

      WRITE (6,*) 'Unit 2 read.  IHOUR2 =',itau1,'   ',LABEL
C****
C**** Compare arrays
C**** ERRQ flags whether any discrepancies have occurred
C****
      ERRQ = .FALSE.
      WRITE (6,*) 'FIELD IM JM LM     VAL1        VAL2     RELERR'
      ERRQ=COMP8 ('U     ',IM,JM,LM     ,U1 ,U2 ) .or. ERRQ
      ERRQ=COMP8 ('V     ',IM,JM,LM     ,V1 ,V2 ) .or. ERRQ
      ERRQ=COMP8 ('T     ',IM,JM,LM     ,T1 ,T2 ) .or. ERRQ
      ERRQ=COMP8 ('Q     ',IM,JM,LM     ,Q1 ,Q2 ) .or. ERRQ
      ERRQ=COMP8 ('P     ',IM,JM,1      ,P1 ,P2 ) .or. ERRQ

      ERRQ=COMP8LIJp('TOCEAN',3,IM,JM      ,TOCEAN1,TOCEAN2) .or. ERRQ
      ERRQ=COMP8 ('Z10   ',IM,JM,1      ,     Z1,     Z2) .or. ERRQ
      ERRQ=COMP8 ('RSI   ',IM,JM,1      ,   RSI1,   RSI2) .or. ERRQ
      ERRQ=COMP8LIJp('TEMPSI',LMI,IM,JM    ,   TSI1,   TSI2) .or. ERRQ
      ERRQ=COMP8p('SNOWI ',IM,JM,1      , SNOWI1, SNOWI2) .or. ERRQ
      ERRQ=COMP8 ('MSI2  ',IM,JM,1      ,   MSI1,   MSI2) .or. ERRQ
      ERRQ=COMP8 ('GDATA ',IM,JM,7      ,GDATA1 ,GDATA2 ) .or. ERRQ
      ERRQ=COMP8 ('GHDATA',IM,JM,4*NGM+5,GHDATA1,GHDATA2) .or. ERRQ
      ERRQ=COMP8 ('BLDATA',IM,JM,11+LM  ,BLDATA1,BLDATA2) .or. ERRQ
      ERRQ=COMPILIJ ('NSN   ',2,IM,JM       ,NSN1  ,NSN2   ) .or. ERRQ
      ERRQ=COMPILIJ ('ISN   ',2,IM,JM       ,ISN1  ,ISN2   ) .or. ERRQ
      ERRQ=COMP8LIJp('DZSN  ',2*NLSN,IM,JM  ,DZSN1 ,DZSN2  ) .or. ERRQ
      ERRQ=COMP8LIJp('WSN   ',2*NLSN,IM,JM  ,WSN1  ,WSN2   ) .or. ERRQ
      ERRQ=COMP8LIJp('HSN   ',2*NLSN,IM,JM  ,HSN1  ,HSN2   ) .or. ERRQ
      ERRQ=COMP8LIJp('FR_SNO',2  ,IM,JM  ,FR_SNOW1,FR_SNOW2) .or. ERRQ
      ERRQ=COMP8 ('LANDI ',IM,JM,3      ,LANDI1 ,LANDI2 ) .or. ERRQ

      ERRQ=COMP8 ('PBL   ',npbl,IM*JM,5*4,PBL1  ,PBL2  ) .or. ERRQ
      ERRQ=COMP8 ('PBLB  ',IM,JM,3*4     ,PBLB1 ,PBLB2 ) .or. ERRQ
      ERRQ=COMP8LIJp('CLOUD ',5*LM,IM,JM    ,CLOUD1,CLOUD2) .or. ERRQ
      ERRQ=COMP8 ('WM    ',IM,JM,LM      ,WM1   ,WM2   ) .or. ERRQ
      ERRQ=COMP8 ('TMOM  ',IM,JM,9*LM    ,TMOM1 ,TMOM2 ) .or. ERRQ
      ERRQ=COMP8 ('QMOM  ',IM,JM,9*LM    ,QMOM1 ,QMOM2 ) .or. ERRQ
      ERRQ=COMP8LIJp('RQT   ',3     ,IM,JM  ,RQT1  ,RQT2  ) .or. ERRQ
      ERRQ=COMP8LIJp('SRHR  ',1+  LM,IM,JM  ,SRHR1 ,SRHR2 ) .or. ERRQ
      ERRQ=COMP8LIJp('TRHR  ',1+  LM,IM,JM  ,TRHR1 ,TRHR2 ) .or. ERRQ
      ERRQ=COMP8LIJp('FSF   ',4     ,IM,JM  ,FSF1  ,FSF2  ) .or. ERRQ
      ERRQ=COMP8 ('LAKE  ',IM,JM,4       ,LAKE1 ,LAKE2 ) .or. ERRQ

c      if(errq) then
c      write(6,*) 'errors in prognostic vars: not checking diagnostics'
c only check diagnostics if no prognostic errors
c      else
      DAGPOS=1
      ERRQ=COMP8 ('AJ    ',JM,94,6  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*94*6
      ERRQ=COMP8 ('AREG  ',24,94,1  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+24*94
      ERRQ=COMP8 ('APJ   ',JM,3,1   ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*3
      ERRQ=COMP8 ('AJL   ',JM,LM,57 ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*LM*57
      ERRQ=COMP8 ('ASJL  ',JM,3,4   ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*3*4
      ERRQ=COMP8 ('AIJ   ',IM,JM,KAIJ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*KAIJ
      ERRQ=COMP8 ('AIL   ',IM,LM,16 ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*LM*16
      ERRQ=COMP8 ('ENERGY',20,100,1 ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+20*100
      ERRQ=COMP8 ('CONSRV',JM,KCON,1  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*KCON
      ERRQ=COMP8 ('SPECA ',IMH+1,20,8,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+(IMH+1)*20*8
      ERRQ=COMP8 ('ATPE  ',8,2,1    ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+8*2
      ERRQ=COMP8 ('ADAILY',24,63,4  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+24*63*4
      ERRQ=COMP8 ('WAVE  ',124,NWD,12,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+2*62*NWD*12
      ERRQ=COMP8 ('AJK   ',JM,LM,KAJK,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*LM*KAJK
      ERRQ=COMP8 ('AIJK  ',IM,JM,LM*6,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*LM*6
      ERRQ=COMP8 ('AIJL  ',IM,JM,LM*5,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*LM*5
      ERRQ=COMP8 ('AJLSP ',JM,LM,(NWD+1)*3 ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*LM*(NWD+1)*3
      ERRQ=COMP8 ('TSFREZ',IM,JM,4      ,TSFREZ1,TSFREZ2)
      ERRQ=COMP8 ('TDIURN',IM,JM,KTD    ,TDIURN1,TDIURN2)
      ERRQ=COMP8p('OA    ',IM,JM,12     ,OA1,OA2)
c      endif
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

      FUNCTION COMP8 (LABEL,IM,JM,LM,X1,X2)
C****
C**** COMP8 compares two arrays and prints any discrepancies to the
C**** line printer.
C**** RETURNS .FALSE. IF NO DISCREPANCIES
C****
      LOGICAL COMP8
      REAL*8 X1(IM,JM,LM),X2(IM,JM,LM), DIF,DIFMAX,X2MAX,X1MAX
      CHARACTER*6 LABEL
      COMP8 = .FALSE.
      NP = 0
      ICNT = 0
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
         IF(DIF.NE.0.) ICNT = ICNT + 1
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
      IF (IMAX.ne.0) then
         WRITE (6,901) IMAX,JMAX,LMAX,X1MAX,X2MAX,DIFMAX
         COMP8 = .TRUE.
      endif
 20   CONTINUE
c     IF(ICNT.GT.0) WRITE(6,*) '#pts = ',ICNT
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I4,E30.20,2E30.20)
      END FUNCTION COMP8

      FUNCTION COMPI (LABEL,IM,JM,LM,I1,I2)
C****
C**** COMPI compares two integer arrays and prints any discrepancies to
C**** the line printer.
C**** RETURNS .FALSE. IF NO DISCREPANCIES
C****
      LOGICAL COMPI
      INTEGER I1(IM,JM,LM),I2(IM,JM,LM),I2MAX,I1MAX
      REAL*8 DIF,DIFMAX
      CHARACTER*6 LABEL
      COMPI = .FALSE.
      NP = 0
      ICNT = 0
      WRITE (6,900) LABEL
      DO 20 L=1,LM
      DIFMAX = 0.
      IMAX=0
      JMAX=0
      LMAX=0
      I1MAX=0
      I2MAX=0
      DO 10 J=1,JM
      DO 10 I=1,IM
c     write(0,*) i,j,l,I2(I,J,L),I1(I,J,L)
      IF (I2(I,J,L)*I1(I,J,L).EQ.0 .AND. I2(I,J,L)+I1(I,J,L).NE.0
     *        .and. I2(I,J,L).NE.I1(I,J,L)) THEN
         WRITE (6,901) I,J,L,I1(I,J,L),I2(I,J,L)
         NP = NP+1
         IF(NP.GE.10)  RETURN
      ELSE
         DIF = ABS(I2(I,J,L)-I1(I,J,L)) / (ABS(I1(I,J,L)) + 1.D-30)
c        IF (DIF.NE.0.) PRINT*,I,J,L,DIF
         IF(DIF.NE.0.) ICNT = ICNT + 1
         IF(DIF.LE.DIFMAX)  GO TO 10
         DIFMAX = DIF
         IMAX=I
         JMAX=J
         LMAX=L
         I1MAX=I1(I,J,L)
         I2MAX=I2(I,J,L)
c      WRITE (6,901) I,J,L,I1(I,J,L),I2(I,J,L),DIF
c      NP = NP+1
c      IF(NP.GE.10)  RETURN
      END IF
 10   CONTINUE
      IF (IMAX.ne.0) then
         WRITE (6,901) IMAX,JMAX,LMAX,I1MAX,I2MAX,DIFMAX
         COMPI = .TRUE.
      endif
 20   CONTINUE
c     IF(ICNT.GT.0) WRITE(6,*) '#pts = ',ICNT
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I4,2I30,E30.20)
      END FUNCTION COMPI

      FUNCTION COMP8LIJp(LABEL,LM,IM,JM,X1,X2)
C****
C**** COMP8LIJp compares two arrays and prints any discrepancies to the
C**** line printer.
C**** RETURNS .FALSE. IF NO DISCREPANCIES
C****
      LOGICAL COMP8LIJp
      REAL*8 X1(LM,IM,JM),X2(LM,IM,JM), DIF,DIFMAX,X2MAX,X1MAX
      CHARACTER*6 LABEL
      COMP8LIJp = .FALSE.
      NP = 0
      ICNT = 0
      WRITE (6,900) LABEL
      DO 20 L=1,LM
      DIFMAX = 0.
      IMAX=0
      JMAX=0
      LMAX=0
      X1MAX=0.
      X2MAX=0.
      DO 10 J=1,JM
      IMofJ=IM
      if(j.eq.1.or.j.eq.jm) IMofJ=1
      DO 10 I=1,IMofJ
      IF (X2(L,I,J)*X1(L,I,J).EQ.0 .AND. X2(L,I,J)+X1(L,I,J).NE.0
     *        .and. X2(L,I,J).NE.X1(L,I,J)) THEN
         WRITE (6,901) L,I,J,X1(L,I,J),X2(L,I,J)
         NP = NP+1
         IF(NP.GE.10)  RETURN
      ELSE
         DIF = ABS(X2(L,I,J)-X1(L,I,J)) / (ABS(X1(L,I,J)) + 1.D-30)
c         IF (DIF.NE.0.) PRINT*,L,I,J,DIF
         IF(DIF.NE.0.) ICNT = ICNT + 1
         IF(DIF.LE.DIFMAX)  GO TO 10
         DIFMAX = DIF
         IMAX=I
         JMAX=J
         LMAX=L
         X1MAX=X1(L,I,J)
         X2MAX=X2(L,I,J)
c      WRITE (6,901) L,I,J,X1(L,I,J),X2(L,I,J),DIF
c      NP = NP+1
c      IF(NP.GE.10)  RETURN
      END IF
 10   CONTINUE
      IF (IMAX.ne.0) then
         WRITE (6,901) IMAX,JMAX,LMAX,X1MAX,X2MAX,DIFMAX
         COMP8LIJp = .TRUE.
      endif
 20   CONTINUE
c     IF(ICNT.GT.0) WRITE(6,*) '#pts = ',ICNT
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I4,E30.20,2E30.20)
      END FUNCTION COMP8LIJp

      FUNCTION COMPILIJ (LABEL,LM,IM,JM,I1,I2)
C****
C**** COMPILIJ compares two integer arrays and prints any discrepancies
C**** to the line printer.
C**** RETURNS .FALSE. IF NO DISCREPANCIES
C****
      LOGICAL COMPILIJ
      INTEGER I1(LM,IM,JM),I2(LM,IM,JM),I2MAX,I1MAX
      REAL*8 DIF,DIFMAX
      CHARACTER*6 LABEL
      COMPILIJ = .FALSE.
      NP = 0
      ICNT = 0
      WRITE (6,900) LABEL
      DO 20 L=1,LM
      DIFMAX = 0.
      IMAX=0
      JMAX=0
      LMAX=0
      I1MAX=0
      I2MAX=0
      DO 10 J=1,JM
      DO 10 I=1,IM
c     write(0,*) L,I,J,I2(L,I,J),I1(L,I,J)
      IF (I2(L,I,J)*I1(L,I,J).EQ.0 .AND. I2(L,I,J)+I1(L,I,J).NE.0
     *        .and. I2(L,I,J).NE.I1(L,I,J)) THEN
         WRITE (6,901) L,I,J,I1(L,I,J),I2(L,I,J)
         NP = NP+1
         IF(NP.GE.10)  RETURN
      ELSE
         DIF = ABS(I2(L,I,J)-I1(L,I,J)) / (ABS(I1(L,I,J)) + 1.D-30)
c        IF (DIF.NE.0.) PRINT*,L,I,J,DIF
         IF(DIF.NE.0.) ICNT = ICNT + 1
         IF(DIF.LE.DIFMAX)  GO TO 10
         DIFMAX = DIF
         IMAX=I
         JMAX=J
         LMAX=L
         I1MAX=I1(L,I,J)
         I2MAX=I2(L,I,J)
c      WRITE (6,901) L,I,J,I1(L,I,J),I2(L,I,J),DIF
c      NP = NP+1
c      IF(NP.GE.10)  RETURN
      END IF
 10   CONTINUE
      IF (IMAX.ne.0) then
         WRITE (6,901) IMAX,JMAX,LMAX,I1MAX,I2MAX,DIFMAX
         COMPILIJ = .TRUE.
      endif
 20   CONTINUE
c     IF(ICNT.GT.0) WRITE(6,*) '#pts = ',ICNT
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I4,2I30,E30.20)
      END FUNCTION COMPILIJ

      FUNCTION COMP8p (LABEL,IM,JM,LM,X1,X2)
C****
C**** COMP8p compares two arrays and prints any discrepancies to the
C**** line printer.
C**** RETURNS .FALSE. IF NO DISCREPANCIES
C****
      LOGICAL COMP8p
      REAL*8 X1(IM,JM,LM),X2(IM,JM,LM), DIF,DIFMAX,X2MAX,X1MAX
      CHARACTER*6 LABEL
      COMP8p = .FALSE.
      NP = 0
      ICNT = 0
      WRITE (6,900) LABEL
      DO 20 L=1,LM
      DIFMAX = 0.
      IMAX=0
      JMAX=0
      LMAX=0
      X1MAX=0.
      X2MAX=0.
      DO 10 J=1,JM
      IMofJ=IM
      IF(J.EQ.1.or.J.eq.JM) IMofJ=1
      DO 10 I=1,IMofJ
      IF (X2(I,J,L)*X1(I,J,L).EQ.0 .AND. X2(I,J,L)+X1(I,J,L).NE.0
     *        .and. X2(I,J,L).NE.X1(I,J,L)) THEN
         WRITE (6,901) I,J,L,X1(I,J,L),X2(I,J,L)
         NP = NP+1
         IF(NP.GE.10)  RETURN
      ELSE
         DIF = ABS(X2(I,J,L)-X1(I,J,L)) / (ABS(X1(I,J,L)) + 1.D-30)
c         IF (DIF.NE.0.) PRINT*,I,J,L,DIF
         IF(DIF.NE.0.) ICNT = ICNT + 1
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
      IF (IMAX.ne.0) then
         WRITE (6,901) IMAX,JMAX,LMAX,X1MAX,X2MAX,DIFMAX
         COMP8p = .TRUE.
      endif
 20   CONTINUE
c     IF(ICNT.GT.0) WRITE(6,*) '#pts = ',ICNT
      RETURN
  900 FORMAT (' ',A6)
  901 FORMAT (3I4,E30.20,2E30.20)
      END FUNCTION COMP8p
