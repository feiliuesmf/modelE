C**** CMPE001.F    CoMPare restartfiles for modelE          6/00
C****
      USE MODEL_COM, only : im,jm,lm,ntype,imh
      USE DAGCOM, ONLY: kacc,ktsf,KTD,kaj,nreg,kapj
     *  ,kajl,lm_req,KASJL,KAIJ,KAIL,NEHIST,HIST_DAYS,KCON
     *  ,KSPECA,NSPHER,KTPE,NHEMI,HR_IN_DAY,NDIUVAR,NDIUPT
     *  ,Max12HR_sequ,NWAV_DAG,KWP,KAJK,KAIJK
      USE SLE001, ONLY: NGM,nlsn
!!!   USE RADNCB, only : LM_REQ
      USE SEAICE_COM, only : lmi
      USE PBLCOM, only : npbl
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE

      COMMON/TADV/TMOM1(IM,JM,9*LM),TMOM2(IM,JM,9*LM)
      COMMON/QADV/QMOM1(IM,JM,9*LM),QMOM2(IM,JM,9*LM)
      COMMON/BNDYCB1/TOCEAN1(3,IM,JM),GDATA1(IM,JM,7),
     *     BLDATA1(IM,JM,11+LM),Z1(IM,JM)
      COMMON/BNDYCB2/TOCEAN2(3,IM,JM),GDATA2(IM,JM,7),
     *     BLDATA2(IM,JM,11+LM),Z2(IM,JM)
ccc   ice data:
      REAL*8 MSI1,MSI2
      LOGICAL IFLAG1(IM,JM),IFLAG2(IM,JM)
      COMMON/SICECB/ RSI1(IM,JM),HSI1(LMI,IM,JM),MSI1(IM,JM),SNOWI1(IM
     *     ,JM),SSI1(LMI,IM,JM),PM1(IM,JM),RSI2(IM,JM),HSI2(LMI,IM,JM)
     *     ,MSI2(IM,JM),SNOWI2(IM,JM),SSI2(LMI,IM,JM),PM2(IM,JM)
      COMMON/RADNCB1/ RQT1( 3,IM,JM),RQT2( 3,IM,JM),
     *               SRHR1(1+LM,IM,JM),SRHR2(1+LM,IM,JM),
     *               TRHR1(1+LM,IM,JM),TRHR2(1+LM,IM,JM),
     *               FSF1(   4,IM,JM),FSF2(   4,IM,JM)
ccc   snow data:
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
      CHARACTER XLABEL*132,LABEL*16,FILEIN*60,HEADER*80
      EQUIVALENCE (XLABEL,LABEL)
      COMMON/DAG1/DIAG1(KACC),TSFREZ1(IM,JM,ktsf),TDIURN1(IM,JM,KTD)
      COMMON/DAG2/DIAG2(KACC),TSFREZ2(IM,JM,ktsf),TDIURN2(IM,JM,KTD)
      COMMON /KEYS/ KEYNR(1+42*50)  !  also incl. keyct
      REAL*8 LAKE1,LAKE2
      COMMON /LAKE/LAKE1(IM,JM,4),LAKE2(IM,JM,4)
C**** coupled model ocean data (cannot be read in from module for
C**** compatibility across model configurations)
      INTEGER, PARAMETER :: NMST=12, LMO=13
      INTEGER, PARAMETER :: KOIJ=5,KOIJL=22,KOL=6,KOLNST=8,
     *     KACCO=IM*JM*KOIJ + IM*JM*LMO*KOIJL + LMO*KOL + LMO*NMST
     *     *KOLNST
      REAL*8 OCEAN1(IM,JM,LMO*11+1),OCEAN2(IM,JM,LMO*11+1)
      REAL*8 STRAITS1(LMO,NMST,7),STRAITS2(LMO,NMST,7)
      REAL*8 STRAITI1(NMST,4+LMI),STRAITI2(NMST,4+LMI)
      REAL*8 ODIAG1(KACCO),ODIAG2(KACCO)
      COMMON /ODAG1/OIJ1(IM,JM,KOIJ),OIJL1(IM,JM,LMO,KOIJL),
     *     OL1(LMO,KOL),OLNST1(LMO,NMST,KOLNST)
      COMMON /ODAG2/OIJ2(IM,JM,KOIJ),OIJL2(IM,JM,LMO,KOIJL),
     *     OL2(LMO,KOL),OLNST2(LMO,NMST,KOLNST)
      EQUIVALENCE (ODIAG1,OIJ1),(ODIAG2,OIJ2)
C**** ice dynamic data (cannot be read in from module for
C**** compatibility across model configurations)
      INTEGER, PARAMETER :: KICIJ=12, IMIC=IM,JMIC=JM
      REAL*8 ICDIAG1(IMIC,JMIC,KICIJ),ICDIAG2(IMIC,JMIC,KICIJ)
      REAL*8 ICEDYN1(IMIC,JMIC,4),ICEDYN2(IMIC,JMIC,4)
C****
      INTEGER DAGPOS,DAGPOS1,DAGPOS2
      LOGICAL ERRQ,COMP8,COMP8p,COMPI,COMP8LIJp,COMPILIJ
      INTEGER itau1,itau2,idacc1(12),idacc2(12)

      real*8 strat1(im,jm),strat2(im,jm)
      integer istrat1(2,im,jm),istrat2(2,im,jm)
C****
      IF(IARGC().NE.2)  GO TO 800
C****
C**** Read ReStartFiles
C****
      CALL GETARG (1,FILEIN)
      OPEN (1,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
c        write(0,*) 'trying to read label'
         READ (1) ITAU1,XLABEL
c        write(0,*) 'trying to skip label'
         READ (1)
c        write(0,*) 'trying to skip param'
         READ (1) ! - skip parameters
c        write(0,*) 'trying to read model'
         READ (1) HEADER,U1,V1,T1,P1,Q1,WM1
C**** check whether stratosphere
c        write(0,*) 'trying to read stratosphere'
         READ (1) HEADER
         BACKSPACE(1)
         IF (HEADER(1:8).eq."STRAT01") THEN
c        write(0,*) 'trying to read stratosphere'
           READ(1) HEADER,STRAT1,istrat1
         END IF
C**** check which ocean
c        write(0,*) 'trying to read ocean'
         READ (1) HEADER
         BACKSPACE(1)
         IF (HEADER(1:8).eq."OCN01") THEN ! Qflux or fixed SST
           KOCEAN1 = 1
c        write(0,*) 'trying to read ocea1'
           READ(1) HEADER,TOCEAN1,Z1
         ELSE
           KOCEAN1 = 2
c        write(0,*) 'trying to read ocea2'
           READ(1) HEADER,OCEAN1
           READ(1) HEADER,STRAITS1,STRAITI1
         END IF
c        write(0,*) 'trying to read lake'
         READ (1) HEADER,LAKE1
c        write(0,*) 'trying to read sice'
         READ (1) HEADER,RSI1,HSI1,SNOWI1,MSI1,SSI1,PM1,IFLAG1
c        write(0,*) 'trying to read gdata'
         READ (1) HEADER,GDATA1
c        write(0,*) 'trying to read soils'
         READ (1) HEADER,GHDATA1
c        write(0,*) 'trying to read snow'
         READ (1) HEADER,NSN1,ISN1,DZSN1,WSN1,HSN1,FR_SNOW1
c        write(0,*) 'trying to read landi'
         READ (1) HEADER,LANDI1
c        write(0,*) 'trying to read bldat'
         READ (1) HEADER,BLDATA1
c        write(0,*) 'trying to read pbl'
         READ (1) HEADER,PBL1,pblb1,ipbl1
c        write(0,*) 'trying to read clds'
         READ (1) HEADER,CLOUD1
c        write(0,*) 'trying to read mom'
         READ (1) HEADER,TMOM1,QMOM1
c        write(0,*) 'trying to read radia'
         READ (1) HEADER,RQT1, S0,SRHR1,TRHR1,FSF1
c        write(0,*) 'trying to read diag'
         IF (KOCEAN1.eq.2) READ(1) HEADER,ICEDYN1
         READ (1,ERR=100) HEADER,KEYNR,TSFREZ1,idacc1,DIAG1,TDIURN1,OA1
     *        ,ITAU2
         GOTO 200
 100     BACKSPACE(1)
         READ (1) HEADER,KEYNR,TSFREZ1,ITAU2
 200     IF (KOCEAN1.eq.2) THEN
c        write(0,*) 'trying to read ocn3'
           READ(1) HEADER,ODIAG1,itau2    !OIJ,OIJL,OL,OLNST,it
           READ(1) HEADER,ICDIAG1
         END IF

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
c        write(0,*) 'trying to read label'
         READ (2) ITAU1,XLABEL
c        write(0,*) 'trying to skip label'
         READ (2)
c        write(0,*) 'trying to skip param'
         READ (2) ! - skip parameters
c        write(0,*) 'trying to read model'
         READ (2) HEADER,U2,V2,T2,P2,Q2,WM2
C**** check whether stratosphere
c        write(0,*) 'trying to read stratosphere'
         READ (2) HEADER
         BACKSPACE(2)
         IF (HEADER(1:8).eq."STRAT01") THEN
c        write(0,*) 'trying to read stratosphere'
           READ(2) HEADER,STRAT2,istrat2
         END IF
C**** check which ocean
c        write(0,*) 'trying to read ocean'
         READ (2) HEADER
         BACKSPACE(2)
         IF (HEADER(1:8).eq."OCN01") THEN ! Qflux or fixed SST
           KOCEAN2 = 1
c        write(0,*) 'trying to read ocea1'
           READ(2) HEADER,TOCEAN2,Z2
         ELSE
           KOCEAN2 = 2
c        write(0,*) 'trying to read ocea2'
           READ(2) HEADER,OCEAN2
           READ(2) HEADER,STRAITS2,STRAITI2
         END IF
c        write(0,*) 'trying to read lake'
         READ (2) HEADER,LAKE2
c        write(0,*) 'trying to read sice'
         READ (2) HEADER,RSI2,HSI2,SNOWI2,MSI2,SSI2,PM2,IFLAG2
c        write(0,*) 'trying to read gdata'
         READ (2) HEADER,GDATA2
c        write(0,*) 'trying to read soils'
         READ (2) HEADER,GHDATA2
c        write(0,*) 'trying to read snow'
         READ (2) HEADER,NSN2,ISN2,DZSN2,WSN2,HSN2,FR_SNOW2
c        write(0,*) 'trying to read landi'
         READ (2) HEADER,LANDI2
c        write(0,*) 'trying to read bldat'
         READ (2) HEADER,BLDATA2
c        write(0,*) 'trying to read pbl'
         READ (2) HEADER,PBL2,pblb2,ipbl2
c        write(0,*) 'trying to read clds'
         READ (2) HEADER,CLOUD2
c        write(0,*) 'trying to read mom'
         READ (2) HEADER,TMOM2,QMOM2
c        write(0,*) 'trying to read radia'
         READ (2) HEADER,RQT2, S0,SRHR2,TRHR2,FSF2
c        write(0,*) 'trying to read diag'
         IF (KOCEAN2.eq.2) READ(2) HEADER,ICEDYN2
         READ (2,ERR=300) HEADER,KEYNR,TSFREZ2,idacc2,DIAG2,TDIURN2,OA2
     *        ,ITAU2
         GOTO 400
 300     BACKSPACE(2)
         READ (2) HEADER,KEYNR,TSFREZ2,ITAU2
 400     IF (KOCEAN2.eq.2) THEN
c        write(0,*) 'trying to read ocn3'
           READ(2) HEADER,ODIAG2,itau2    !OIJ,OIJL,OL,OLNST,it
           READ(2) HEADER,ICDIAG2
         END IF

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
      ERRQ=COMP8 ('STRAT ',IM,JM,1      ,strat1 ,strat2 ) .or. ERRQ
      ERRQ=COMPi ('ISTRAT',2,IM,JM      ,istrat1,istrat2 ) .or. ERRQ
      IF (KOCEAN1.eq.KOCEAN2) THEN ! compare oceans else don't
        IF (KOCEAN1.eq.1) THEN  ! Qflux/fixed
          ERRQ=COMP8LIJp('TOCEAN',3,IM,JM ,TOCEAN1,TOCEAN2).or.ERRQ
          ERRQ=COMP8 ('Z10   ',IM,JM,1    ,     Z1,     Z2).or.ERRQ
        ELSE                    ! coupled
          ERRQ=COMP8LIJp('OCEAN ',IM,JM,11*LMO+1,OCEAN1,OCEAN2).or.ERRQ
          ERRQ=COMP8LIJp('STRATO',LMO,NMST,7,STRAITS1,STRAITS2).or.ERRQ
          ERRQ=COMP8LIJp('STRATI',NMST,4+LMI,1,STRAITI1,STRAITI2).or
     *         .ERRQ
          ERRQ=COMP8LIJp('ICEDYN',IMIC,JMIC,4,ICEDYN1,ICEDYN2).or.ERRQ
          ERRQ=COMP8LIJp('ICEDAG',IMIC,JMIC,KICIJ,ICDIAG1,ICDIAG2).or
     *         .ERRQ
        END IF
      END IF
      ERRQ=COMP8 ('RSI   ',IM,JM,1      ,   RSI1,   RSI2) .or. ERRQ
      ERRQ=COMP8LIJp('HSI  ',LMI,IM,JM    ,   HSI1,   HSI2) .or. ERRQ
      ERRQ=COMP8LIJp('SSI  ',LMI,IM,JM    ,   SSI1,   SSI2) .or. ERRQ
      ERRQ=COMP8p('SNOWI ',IM,JM,1      , SNOWI1, SNOWI2) .or. ERRQ
      ERRQ=COMP8 ('MSI2  ',IM,JM,1      ,   MSI1,   MSI2) .or. ERRQ
      ERRQ=COMP8 ('MPOND ',IM,JM,1      ,    PM1,    PM2) .or. ERRQ
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
      ERRQ=COMP8 ('AJ    ',JM,kaj,ntype,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*kaj*ntype
      ERRQ=COMP8 ('AREG  ',nreg,kaj,1  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+nreg*kaj
      ERRQ=COMP8 ('APJ   ',JM,kapj,1   ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*kapj
      ERRQ=COMP8 ('AJL   ',JM,LM,kajl ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*LM*kajl
      ERRQ=COMP8 ('ASJL  ',JM,lm_req,KASJL,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*lm_req*KASJL
      ERRQ=COMP8 ('AIJ   ',IM,JM,KAIJ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*KAIJ
      ERRQ=COMP8 ('AIL   ',IM,LM,KAIL,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*LM*KAIL
      ERRQ=COMP8 ('ENERGY',NEHIST,NEHIST,HIST_DAYS
     *                                  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+NEHIST*HIST_DAYS
      ERRQ=COMP8 ('CONSRV',JM,KCON,1  ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*KCON
      ERRQ=COMP8 ('SPECA ',IMH+1,KSPECA,NSPHER
     *                               ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+(IMH+1)*KSPECA*NSPHER
      ERRQ=COMP8 ('ATPE  ',KTPE,NHEMI,1,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+KTPE*NHEMI
      ERRQ=COMP8 ('ADAILY',HR_IN_DAY,NDIUVAR,NDIUPT
     *                              ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+HR_IN_DAY*NDIUVAR*NDIUPT
      ERRQ=COMP8 ('WAVE  ',2*Max12HR_sequ,NWAV_DAG,KWP
     *                               ,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+2*Max12HR_sequ*NWAV_DAG*KWP
      ERRQ=COMP8 ('AJK   ',JM,LM,KAJK,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+JM*LM*KAJK
      ERRQ=COMP8 ('AIJK  ',IM,JM,LM*6,DIAG1(DAGPOS),DIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*LM*6
      ERRQ=COMP8 ('TSFREZ',IM,JM,ktsf   ,TSFREZ1,TSFREZ2)
      ERRQ=COMP8 ('TDIURN',IM,JM,KTD    ,TDIURN1,TDIURN2)
      ERRQ=COMP8p('OA    ',IM,JM,12     ,OA1,OA2)

      IF (KOCEAN1.eq.KOCEAN2.and.KOCEAN1.eq.2) THEN ! compare ocn diags
      DAGPOS=1
      ERRQ=COMP8('OIJ   ',IM,JM,KOIJ,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*KOIJ
      ERRQ=COMP8('OIJL  ',IM,JM,LMO*KOIJL,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*LMO*KOIJL
      ERRQ=COMP8('OL    ',LMO,KOL,1   ,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      DAGPOS=DAGPOS+LMO*KOL
      ERRQ=COMP8('OLNST ',LMO,NMST,KOLNST,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      END IF

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
