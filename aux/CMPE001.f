#include "../model/rundeck_opts.h"

C**** CMPE001.F    CoMPare restartfiles for modelE          6/00
C****
      USE MODEL_COM, only : im,jm,lm,ntype,imh
      USE DAGCOM, ONLY: kacc,ktsf,KTD,KOA,kaj,nreg,kapj
     *  ,kajl,lm_req,KASJL,KAIJ,KAIL,NEHIST,HIST_DAYS,KCON
     *  ,KSPECA,NSPHER,KTPE,NHEMI,HR_IN_DAY,NDIUVAR,NDIUPT
     *  ,Max12HR_sequ,NWAV_DAG,KWP,KAJK,KAIJK
      USE SLE001, ONLY: NGM,nlsn
!!!   USE RADNCB, only : LM_REQ
      USE SEAICE_COM, only : lmi
      USE PBLCOM, only : npbl
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm
      USE TRACER_DIAG_COM, only: ktacc,ktaij,ktaijs,ktajlx,ktajls,ktcon
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: JPPJ, n_rx
#endif
#endif
      IMPLICIT NONE
      SAVE

      REAL*8 TMOM1,TMOM2,QMOM1,QMOM2,TOCN1,z1,sne1,te1,wtre1,ace1,snag1
     *     ,evm1,fsat1,qge1,bld1,eg1,we1,tg1,qg1,TOCN2,z2,sne2,te2,wtre2
     *     ,ace2,snag2,evm2,fsat2,qge2,bld2,eg2,we2,tg2,qg2
      COMMON/TADV/TMOM1(IM,JM,9*LM),TMOM2(IM,JM,9*LM)
      COMMON/QADV/QMOM1(IM,JM,9*LM),QMOM2(IM,JM,9*LM)
      COMMON/BNDYCB1/TOCN1(3,IM,JM),z1(IM,JM)
     *     ,sne1(IM,JM),te1(im,jm),wtre1(im,jm),ace1(im,jm),snag1(3,im
     *     ,jm),evm1(im,jm),fsat1(im,jm),qge1(im,jm)       ,BLD1(IM,JM
     *     ,11),eg1(lm,im,jm),we1(lm,im,jm),tg1(im,jm),qg1(im,jm)
      COMMON/BNDYCB2/TOCN2(3,IM,JM),z2(IM,JM)
     *     ,sne2(IM,JM),te2(im,jm),wtre2(im,jm),ace2(im,jm),snag2(3,im
     *     ,jm),evm2(im,jm),fsat2(im,jm),qge2(im,jm)       ,BLD2(IM,JM
     *     ,11),eg2(lm,im,jm),we2(lm,im,jm),tg2(im,jm),qg2(im,jm)
ccc   ice data:
      REAL*8 MSI1,MSI2,RSI1,RSI2,HSI1,HSI2,SSI1,PM1,SSI2,PM2,SNOWI1
     *     ,SNOWI2
      LOGICAL IFLAG1(IM,JM),IFLAG2(IM,JM)
      COMMON/SICECB/ RSI1(IM,JM),HSI1(LMI,IM,JM),MSI1(IM,JM),SNOWI1(IM
     *     ,JM),SSI1(LMI,IM,JM),PM1(IM,JM),RSI2(IM,JM),HSI2(LMI,IM,JM)
     *     ,MSI2(IM,JM),SNOWI2(IM,JM),SSI2(LMI,IM,JM),PM2(IM,JM)
      REAL*8 RQT1,RQT2,SRHR1,SRHR2,TRHR1,TRHR2,FSF1,FSF2,FSD1,FSD2,S01
     *     ,S02,RCLD1,RCLD2
      COMMON/RADNCB1/ RQT1( 3,IM,JM),RQT2( 3,IM,JM),
     *               SRHR1(1+LM,IM,JM),SRHR2(1+LM,IM,JM),
     *               TRHR1(1+LM,IM,JM),TRHR2(1+LM,IM,JM),
     *               FSF1(   4,IM,JM),FSF2(   4,IM,JM),
     *               FSD1(IM,JM,5),FSD2(IM,JM,5),
     *               RCLD1(LM,IM,JM),RCLD2(LM,IM,JM)
ccc   snow data:
      INTEGER, DIMENSION(2,IM,JM)     :: NSN1, NSN2
      INTEGER, DIMENSION(2,IM,JM)     :: ISN1, ISN2
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: DZSN1, DZSN2
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: WSN1, WSN2
      REAL*8, DIMENSION(NLSN,2,IM,JM) :: HSN1, HSN2
      REAL*8, DIMENSION(2,IM,JM)      :: FR_SNOW1, FR_SNOW2
ccc   land ice data
      REAL*8 SNLI1,TLI1,SNLI2,TLI2,OA1,OA2
      COMMON/LNDICB1/SNLI1(IM,JM),TLI1(2,IM,JM)
      COMMON/LNDICB2/SNLI2(IM,JM),TLI2(2,IM,JM)
      COMMON/WORKO/  OA1(IM,JM,KOA),OA2(IM,JM,KOA)
      REAL*8 :: U1(IM,JM,LM),V1(IM,JM,LM),T1(IM,JM,LM),P1(IM,JM),Q1(IM
     *     ,JM,LM),WM1(IM,JM,LM)
      REAL*8 :: U2(IM,JM,LM),V2(IM,JM,LM),T2(IM,JM,LM),P2(IM,JM),Q2(IM
     *     ,JM,LM),WM2(IM,JM,LM)
      REAL*8 pbl1,pblb1,pbl2,pblb2,CLOUD1,CLOUD2
      INTEGER ipbl1,ipbl2
      common /socabl1/pbl1(npbl,im,jm,5*4),pblb1(im,jm,3*4),ipbl1(im,jm
     *     ,4)
      common /socabl2/pbl2(npbl,im,jm,5*4),pblb2(im,jm,3*4),ipbl2(im,jm
     *     ,4)
      COMMON/CLDCOM/CLOUD2(IM,JM,5*LM),CLOUD1(IM,JM,5*LM)
      REAL*8 wb1,wv1,htb1,htv1,snbv1,wb2,wv2,htb2,htv2,snbv2
      real*8 ci1(im,jm),qfol1(im,jm),cnc1(im,jm)
     $     ,ci2(im,jm),qfol2(im,jm),cnc2(im,jm)
      COMMON/SOILS1/wb1(ngm,IM,JM),wv1(ngm+1,im,jm),htb1(ngm+1,im,jm),
     *  htv1(ngm+1,im,jm),snbv1(2,im,jm)
      COMMON/SOILS2/wb2(ngm,IM,JM),wv2(ngm+1,im,jm),htb2(ngm+1,im,jm),
     *  htv2(ngm+1,im,jm),snbv2(2,im,jm)
      CHARACTER XLABEL*132,LABEL*16,FILEIN*60,HEADER*80
      EQUIVALENCE (XLABEL,LABEL)
      REAL*8 DIAG1,DIAG2,TSFREZ1,TSFREZ2,TDIURN1,TDIURN2
      INTEGER KEYNR
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
      REAL*8 ODIAG1(KACCO),ODIAG2(KACCO),OIJ1,OIJ2,OIJL1,OIJL2,OL1,OL2
     *     ,OLNST1,OLNST2
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
#ifdef TRACERS_ON
      REAL*8 TR1(IM,JM,LM*NTM),TR2(IM,JM,LM*NTM),TRMOM1(9,IM,JM,LM*NTM)
     *     ,TRMOM2(9,IM,JM,LM*NTM),
     *     TRABL1(npbl*ntm,im,jm,4),TRABL2(npbl*ntm,im,jm,4)
      REAL*8 TRACC1(KTACC),TRACC2(KTACC),TAIJLN1,TAIJN1,TAIJS1,TAJLN1
     *     ,TAJLS1,TCON1,TAIJLN2,TAIJN2,TAIJS2,TAJLN2,TAJLS2,TCON2 
      COMMON /TACCUM1/ TAIJLN1(IM,JM,LM,NTM),TAIJN1(IM,JM,KTAIJ,NTM)
     *     ,TAIJS1(IM,JM,ktaijs),TAJLN1(JM,LM,ktajlx,NTM),TAJLS1(JM,LM
     *     ,ktajls),TCON1(JM,ktcon,ntm)
      COMMON /TACCUM2/ TAIJLN2(IM,JM,LM,NTM),TAIJN2(IM,JM,KTAIJ,NTM)
     *     ,TAIJS2(IM,JM,ktaijs),TAJLN2(JM,LM,ktajlx,NTM),TAJLS2(JM,LM
     *     ,ktajls),TCON2(JM,ktcon,ntm)
      EQUIVALENCE (TRACC1,TAIJLN1),(TRACC2,TAIJLN2)
#ifdef TRACERS_WATER
      INTEGER, PARAMETER :: KTICIJ=2
      REAL*8 TRW1(IM,JM,LM*NTM),TRW2(IM,JM,LM*NTM),
     *     TRLK1(2*NTM,IM,JM),TRLK2(2*NTM,IM,JM),
     *     TRSI1(LMI*NTM,IM,JM),TRSI2(LMI*NTM,IM,JM),
     *     TRLI1(2*NTM,IM,JM),TRLI2(2*NTM,IM,JM),
     *     TRSOIL1(NTM*(2*NGM+3),IM,JM),TRSOIL2(NTM*(2*NGM+3),IM,JM),
     *     TRSN1(2*NTM*NLSN,IM,JM),TRSN2(2*NTM*NLSN,IM,JM)
      REAL*8 TRICDG1(IMIC,JMIC,KTICIJ*NTM),TRICDG2(IMIC,JMIC,KTICIJ*NTM)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      INTEGER K1
      REAL*8, DIMENSION(IM,JM,LM):: yNO31,pHOx1,pNOx1,pOx1,yCH3O21,
     &     yC2O31,yROR1,yXO21,yAldehyde1,yXO2N1,yRXPAR1,
     &     yNO32,pHOx2,pNOx2,pOx2,yCH3O22,yC2O32,yROR2,yXO22,
     &     yAldehyde2,yXO2N2,yRXPAR2,temp1,temp2
      REAL*8, DIMENSION(JM,4,12) :: corrOx1,corrOx2
      REAL*8, DIMENSION(LM,JM,IM):: O3DLJI1,O3DLJI_clim1,
     &     O3DLJI2,O3DLJI_clim2
      REAL*8, DIMENSION(jppj,IM,JM,LM) :: ss1,ss2
#endif
#endif
C****
      INTEGER DAGPOS,DAGPOS1,DAGPOS2,KOCEAN1,KOCEAN2,IARGC
      LOGICAL ERRQ,COMP8,COMP8p,COMPI,COMP8LIJp,COMPILIJ
      INTEGER itau1,itau2,idacc1(12),idacc2(12)

      real*8 strat1(im,jm),strat2(im,jm)
      integer istrat1(2,im,jm),istrat2(2,im,jm)
      logical :: debug = .false.
C****
      IF(IARGC().NE.2)  GO TO 800
C****
C**** Read ReStartFiles
C****
      CALL GETARG (1,FILEIN)
      OPEN (1,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
         if (debug) write(0,*) 'trying to read label'
         READ (1) ITAU1,XLABEL
         if (debug) write(0,*) 'trying to skip label'
         READ (1)
         if (debug) write(0,*) 'trying to skip param'
         READ (1) ! - skip parameters
         if (debug) write(0,*) 'trying to read model'
         READ (1) HEADER,U1,V1,T1,P1,Q1,WM1
C**** check whether stratosphere
         if (debug) write(0,*) 'trying to read stratosphere'
         READ (1) HEADER
         BACKSPACE(1)
         IF (HEADER(1:8).eq."STRAT01") THEN
           if (debug) write(0,*) 'trying to read stratosphere'
           READ(1) HEADER,STRAT1,istrat1
         END IF
C**** check which ocean
         if (debug) write(0,*) 'trying to read ocean'
         READ (1) HEADER
         BACKSPACE(1)
         IF (HEADER(1:8).eq."OCN01") THEN ! Qflux or fixed SST
           KOCEAN1 = 0
           if (debug) write(0,*) 'trying to read ocea1'
           READ(1) HEADER,TOCN1,Z1
         ELSE
           KOCEAN1 = 2
           if (debug) write(0,*) 'trying to read ocea2'
           READ(1) HEADER,OCEAN1
           READ(1) HEADER,STRAITS1,STRAITI1
         END IF
         if (debug) write(0,*) 'trying to read lake'
         READ (1) HEADER,LAKE1
#ifdef TRACERS_WATER
         READ (1) HEADER,TRLK1
#endif
         if (debug) write(0,*) 'trying to read sice'
         READ (1) HEADER,RSI1,HSI1,SNOWI1,MSI1,SSI1,PM1,IFLAG1
         if (debug) write(0,*) 'trying to read gdata'
#ifdef TRACERS_WATER
         READ (1) HEADER,TRSI1 
#endif
         READ (1) HEADER,sne1,te1,wtre1,ace1,snag1,fsat1,qge1
         if (debug) write(0,*) 'trying to read soils'
         READ (1) HEADER,wb1,wv1,htb1,htv1,snbv1
         if (debug) write(0,*) 'trying to read vegetation'
         READ (1) HEADER,ci1,qfol1,cnc1
#ifdef TRACERS_WATER
         READ (1) HEADER,TRSOIL1
#endif
         if (debug) write(0,*) 'trying to read snow'
         READ (1) HEADER,NSN1,ISN1,DZSN1,WSN1,HSN1,FR_SNOW1
#ifdef TRACERS_WATER
         READ (1) HEADER,TRSN1
#endif
         if (debug) write(0,*) 'trying to read landi'
         READ (1) HEADER,SNLI1,TLI1
#ifdef TRACERS_WATER
         READ (1) HEADER,TRLI1
#endif
         if (debug) write(0,*) 'trying to read bldat'
         READ (1) HEADER,BLD1,eg1,we1,tg1,qg1
         if (debug) write(0,*) 'trying to read pbl'
         READ (1) HEADER,PBL1,pblb1,ipbl1
#ifdef TRACERS_ON
         READ (1) HEADER,TRABL1
#endif
         if (debug) write(0,*) 'trying to read clds'
         READ (1) HEADER,CLOUD1
         if (debug) write(0,*) 'trying to read mom'
         READ (1) HEADER,TMOM1,QMOM1
         if (debug) write(0,*) 'trying to read radia'
         READ (1) HEADER,RQT1, S01,SRHR1,TRHR1,FSF1,fsd1,RCLD1
         if (debug) write(0,*) 'trying to read icedyn'
         READ (1) HEADER
         BACKSPACE(1)
         IF (HEADER(1:6).eq."ICEDYN".and.KOCEAN1.eq.0) KOCEAN1=1
         IF (HEADER(1:6).eq."ICEDYN") READ(1) HEADER,ICEDYN1
#ifdef TRACERS_ON
         if (debug) write(0,*) 'trying to read tracers'
         READ (1) HEADER,TR1,TRMOM1
#ifdef TRACERS_WATER
     *        ,TRW1
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *        ,yNO31,pHOx1,pNOx1,pOx1,yCH3O21,yC2O31,yROR1,yXO21
     *        ,yAldehyde1,yXO2N1,yRXPAR1,corrOx1,O3DLJI1,O3DLJI_clim1
     *        ,ss1
#endif
#endif
         if (debug) write(0,*) 'trying to read diag'
         READ (1,ERR=100) HEADER,KEYNR,TSFREZ1,idacc1,DIAG1,TDIURN1,OA1
     *        ,ITAU2
         GOTO 200
 100     BACKSPACE(1)
         READ (1) HEADER,KEYNR,TSFREZ1,ITAU2
 200     IF (KOCEAN1.gt.0) THEN
           if (debug) write(0,*) 'trying to read ocn3'
           IF (KOCEAN1.eq.2) READ(1) HEADER,ODIAG1,itau2
           READ(1) HEADER,ICDIAG1
#ifdef TRACERS_WATER
           READ (1) HEADER,TRICDG1
#endif
         END IF
#ifdef TRACERS_ON
         READ (1) HEADER,TRACC1
#endif
         IF (ITAU1.ne.ITAU2) then
           WRITE (6,*) 'FILE 1 NOT READ CORRECTLY. IHOUR,IHOURB =',itau1
     *          ,itau2
         STOP
      END IF
      CLOSE (1)

      WRITE (6,*) 'Unit 1 read.  IHOUR1 =',itau1,'   ',LABEL
C****
      CALL GETARG (2,FILEIN)
      OPEN (2,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=810)
c        if (debug) write(0,*) 'trying to read label'
         READ (2) ITAU1,XLABEL
c        if (debug) write(0,*) 'trying to skip label'
         READ (2)
c        if (debug) write(0,*) 'trying to skip param'
         READ (2) ! - skip parameters
c        if (debug) write(0,*) 'trying to read model'
         READ (2) HEADER,U2,V2,T2,P2,Q2,WM2
C**** check whether stratosphere
c        if (debug) write(0,*) 'trying to read stratosphere'
         READ (2) HEADER
         BACKSPACE(2)
         IF (HEADER(1:8).eq."STRAT01") THEN
c        if (debug) write(0,*) 'trying to read stratosphere'
           READ(2) HEADER,STRAT2,istrat2
         END IF
C**** check which ocean
c        if (debug) write(0,*) 'trying to read ocean'
         READ (2) HEADER
         BACKSPACE(2)
         IF (HEADER(1:8).eq."OCN01") THEN ! Qflux or fixed SST
           KOCEAN2 = 0
c        if (debug) write(0,*) 'trying to read ocea1'
           READ(2) HEADER,TOCN2,Z2
         ELSE
           KOCEAN2 = 2
c        if (debug) write(0,*) 'trying to read ocea2'
           READ(2) HEADER,OCEAN2
           READ(2) HEADER,STRAITS2,STRAITI2
         END IF
c        if (debug) write(0,*) 'trying to read lake'
         READ (2) HEADER,LAKE2
#ifdef TRACERS_WATER
         READ (2) HEADER,TRLK2
#endif
c        if (debug) write(0,*) 'trying to read sice'
         READ (2) HEADER,RSI2,HSI2,SNOWI2,MSI2,SSI2,PM2,IFLAG2
#ifdef TRACERS_WATER
         READ (2) HEADER,TRSI2
#endif
c        if (debug) write(0,*) 'trying to read gdata'
         READ (2) HEADER,sne2,te2,wtre2,ace2,snag2,fsat2,qge2
c        if (debug) write(0,*) 'trying to read soils'
         READ (2) HEADER,wb2,wv2,htb2,htv2,snbv2
c        if (debug) write(0,*) 'trying to read vegetation'
         READ (2) HEADER,ci2,qfol2,cnc2
#ifdef TRACERS_WATER
         READ (2) HEADER,TRSOIL2
#endif
c        if (debug) write(0,*) 'trying to read snow'
         READ (2) HEADER,NSN2,ISN2,DZSN2,WSN2,HSN2,FR_SNOW2
#ifdef TRACERS_WATER
         READ (2) HEADER,TRSN2
#endif
c        if (debug) write(0,*) 'trying to read landi'
         READ (2) HEADER,SNLI2,TLI2
#ifdef TRACERS_WATER
         READ (2) HEADER,TRLI2
#endif
c        if (debug) write(0,*) 'trying to read bldat'
         READ (2) HEADER,BLD2,eg2,we2,tg2,qg2
c        if (debug) write(0,*) 'trying to read pbl'
         READ (2) HEADER,PBL2,pblb2,ipbl2
#ifdef TRACERS_ON
         READ (2)  HEADER,TRABL2
#endif
c        if (debug) write(0,*) 'trying to read clds'
         READ (2) HEADER,CLOUD2
c        if (debug) write(0,*) 'trying to read mom'
         READ (2) HEADER,TMOM2,QMOM2
c        if (debug) write(0,*) 'trying to read radia'
         READ (2) HEADER,RQT2, S02,SRHR2,TRHR2,FSF2,fsd2,rcld2
c        if (debug) write(0,*) 'trying to read icedyn'
         READ (2) HEADER
         BACKSPACE(2)
         IF (HEADER(1:6).eq.'ICEDYN'.and.KOCEAN2.eq.0) KOCEAN2=1
         IF (HEADER(1:6).eq.'ICEDYN') READ(2) HEADER,ICEDYN2
#ifdef TRACERS_ON
         READ (2) HEADER,TR2,TRMOM2
#ifdef TRACERS_WATER
     *        ,TRW2
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *        ,yNO32,pHOx2,pNOx2,pOx2,yCH3O22,yC2O32,yROR2,yXO22
     *        ,yAldehyde2,yXO2N2,yRXPAR2,corrOx2,O3DLJI2,O3DLJI_clim2
     *        ,ss2
#endif
#endif
c        if (debug) write(0,*) 'trying to read diag'
         READ (2,ERR=300) HEADER,KEYNR,TSFREZ2,idacc2,DIAG2,TDIURN2,OA2
     *        ,ITAU2
         GOTO 400
 300     BACKSPACE(2)
         READ (2) HEADER,KEYNR,TSFREZ2,ITAU2
 400     IF (KOCEAN2.gt.0) THEN
c        if (debug) write(0,*) 'trying to read ocn3'
           IF (KOCEAN2.eq.2) READ(2) HEADER,ODIAG2,itau2
           READ(2) HEADER,ICDIAG2
#ifdef TRACERS_WATER
           READ (2) HEADER,TRICDG2
#endif
         END IF
#ifdef TRACERS_ON
         READ (2) HEADER,TRACC2
#endif

      IF (itau1.ne.itau2) then
         WRITE (6,*) 'FILE 2 NOT READ CORRECTLY. IHOUR,IHOURB =',itau1
     *        ,itau2
         STOP
      END IF
      CLOSE (2)

      WRITE (6,*) 'Unit 2 read.  IHOUR2 =',itau1,'   ',LABEL
      WRITE (6,*) 'solar constants ',s01,s02
C****
C**** Compare arrays
C**** ERRQ flags whether any discrepancies have occurred
C****
      print*," Prognostic variables:"
      ERRQ = .FALSE.
      WRITE (6,*) 'FIELD IM JM LM     VAL1        VAL2     RELERR'
      ERRQ=COMP8 ('U     ',IM,JM,LM     ,U1 ,U2 ) .or. ERRQ
      ERRQ=COMP8 ('V     ',IM,JM,LM     ,V1 ,V2 ) .or. ERRQ
      ERRQ=COMP8 ('T     ',IM,JM,LM     ,T1 ,T2 ) .or. ERRQ
      ERRQ=COMP8 ('Q     ',IM,JM,LM     ,Q1 ,Q2 ) .or. ERRQ
      ERRQ=COMP8 ('P     ',IM,JM,1      ,P1 ,P2 ) .or. ERRQ
      ERRQ=COMP8 ('STRAT ',IM,JM,1      ,strat1 ,strat2 ) .or. ERRQ
      ERRQ=COMPiLIJ ('ISTRAT',2,IM,JM      ,istrat1,istrat2 ) .or. ERRQ
      IF (KOCEAN1.eq.KOCEAN2) THEN ! compare oceans else don't
        IF (KOCEAN1.le.1) THEN  ! Qflux/fixed
          ERRQ=COMP8LIJp('TOCN  ',3,IM,JM ,TOCN1,TOCN2).or.ERRQ
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
        IF (KOCEAN1.eq.1) THEN
          ERRQ=COMP8LIJp('ICEDAG',IMIC,JMIC,KICIJ,ICDIAG1,ICDIAG2).or
     *         .ERRQ
        END IF
      END IF
      ERRQ=COMP8 ('LAKE  ',IM,JM,4       ,LAKE1 ,LAKE2 ) .or. ERRQ
      ERRQ=COMP8 ('RSI   ',IM,JM,1      ,   RSI1,   RSI2) .or. ERRQ
      ERRQ=COMP8LIJp('HSI  ',LMI,IM,JM    ,   HSI1,   HSI2) .or. ERRQ
      ERRQ=COMP8LIJp('SSI  ',LMI,IM,JM    ,   SSI1,   SSI2) .or. ERRQ
      ERRQ=COMP8p('SNOWI ',IM,JM,1      , SNOWI1, SNOWI2) .or. ERRQ
      ERRQ=COMP8 ('MSI2  ',IM,JM,1      ,   MSI1,   MSI2) .or. ERRQ
      ERRQ=COMP8 ('MPOND ',IM,JM,1      ,    PM1,    PM2) .or. ERRQ
      ERRQ=COMP8 ('snwe  ',IM,JM,1      ,   sne1,   sne2).or.ERRQ
      ERRQ=COMP8 ('Tearth',IM,JM,1      ,    te1,    te2).or.ERRQ
      ERRQ=COMP8 ('WTRe  ',IM,JM,1      ,  wtre1,  wtre2).or.ERRQ
      ERRQ=COMP8 ('ICEe  ',IM,JM,1      ,   ace1,   ace2).or.ERRQ
      ERRQ=COMP8LIJp('SNWage',3,IM,JM,     snag1,  snag2).or.ERRQ
      ERRQ=COMP8 ('evmax ',IM,JM,1      ,   evm1,   evm2).or.ERRQ
      ERRQ=COMP8 ('fsat  ',IM,JM,1      ,  fsat1,  fsat2).or.ERRQ
      ERRQ=COMP8 ('qge   ',IM,JM,1      ,   qge1,   qge2).or.ERRQ
      ERRQ=COMP8LIJp ('WTRb  ',ngm,IM,JM,       wb1,    wb2).or.ERRQ
      ERRQ=COMP8LIJp ('WTRv  ',1+ngm,IM,JM,     wv1,    wv2).or.ERRQ
      ERRQ=COMP8LIJp ('HEATb ',1+ngm,IM,JM,    htb1,   htb2).or.ERRQ
      ERRQ=COMP8LIJp ('HEATv ',1+ngm,IM,JM,    htv1,   htv2).or.ERRQ
      ERRQ=COMP8LIJp ('SNOWbv',2    ,IM,JM,   snbv1,  snbv2).or.ERRQ
      ERRQ=COMP8 ('Cint  ',IM,JM,1      ,    ci1,    ci2) .or. ERRQ
      ERRQ=COMP8 ('Qfol  ',IM,JM,1      ,  qfol1,  qfol2) .or. ERRQ
      ERRQ=COMP8 ('CNC  ',IM,JM,1      ,  cnc1,  cnc2) .or. ERRQ
      ERRQ=COMPILIJ ('NSN   ',2,IM,JM       ,NSN1  ,NSN2   ) .or. ERRQ
      ERRQ=COMPILIJ ('ISN   ',2,IM,JM       ,ISN1  ,ISN2   ) .or. ERRQ
      ERRQ=COMP8LIJp('DZSN  ',2*NLSN,IM,JM  ,DZSN1 ,DZSN2  ) .or. ERRQ
      ERRQ=COMP8LIJp('WSN   ',2*NLSN,IM,JM  ,WSN1  ,WSN2   ) .or. ERRQ
      ERRQ=COMP8LIJp('HSN   ',2*NLSN,IM,JM  ,HSN1  ,HSN2   ) .or. ERRQ
      ERRQ=COMP8LIJp('FR_SNO',2  ,IM,JM  ,FR_SNOW1,FR_SNOW2) .or. ERRQ
      ERRQ=COMP8 ('snowLI',IM,JM,3      ,SNLI1  ,SNLI2  ) .or. ERRQ
      ERRQ=COMP8LIJp('TLI   ',2     ,IM,JM  ,TLI1  ,TLI2  ) .or. ERRQ

      ERRQ=COMP8 ('BLDATA',IM,JM,11,BLD1,BLD2) .or. ERRQ
      ERRQ=COMP8LIJp('KEtrb',LMI,IM,JM    ,    eg1,    eg2) .or. ERRQ
      ERRQ=COMP8LIJp('KEvrt',LMI,IM,JM    ,    we1,    we2) .or. ERRQ
      ERRQ=COMP8 ('TG1av ',IM,JM,1      ,    tg1,    tg2) .or. ERRQ
      ERRQ=COMP8 ('QG1av ',IM,JM,1      ,    qg1,    qg2) .or. ERRQ
      ERRQ=COMP8 ('PBL   ',npbl,IM*JM,5*4,PBL1  ,PBL2  ) .or. ERRQ
      ERRQ=COMP8 ('PBLB  ',IM,JM,3*4     ,PBLB1 ,PBLB2 ) .or. ERRQ
      ERRQ=COMPI ('ipbl  ',IM,JM,4       ,ipbl1 ,ipbl2 ) .or. ERRQ
      ERRQ=COMP8LIJp('CLOUD ',5*LM,IM,JM    ,CLOUD1,CLOUD2) .or. ERRQ
      ERRQ=COMP8 ('WM    ',IM,JM,LM      ,WM1   ,WM2   ) .or. ERRQ
      ERRQ=COMP8 ('TMOM  ',9,IM*JM,LM    ,TMOM1 ,TMOM2 ) .or. ERRQ
      ERRQ=COMP8 ('QMOM  ',9,IM*JM,LM    ,QMOM1 ,QMOM2 ) .or. ERRQ
      ERRQ=COMP8LIJp('RQT   ',3     ,IM,JM  ,RQT1  ,RQT2  ) .or. ERRQ
      ERRQ=COMP8LIJp('SRHR  ',1+  LM,IM,JM  ,SRHR1 ,SRHR2 ) .or. ERRQ
      ERRQ=COMP8LIJp('TRHR  ',1+  LM,IM,JM  ,TRHR1 ,TRHR2 ) .or. ERRQ
      ERRQ=COMP8LIJp('FSF   ',4     ,IM,JM  ,FSF1  ,FSF2  ) .or. ERRQ
      ERRQ=COMP8 ('FSdir ',IM,JM,5,        fsd1 , fsd2 ) .or. ERRQ
      ERRQ=COMP8Lijp('RCLD  ',LM,IM,JM,    RCLD1, RCLD2 ) .or. ERRQ

#ifdef TRACERS_ON
      ERRQ=COMP8('TR    ',IM,JM,LM*NTM    , TR1  , TR2  ) .or. ERRQ
      ERRQ=COMP8('TRMOM ',9*IM,JM,LM*NTM  ,TRMOM1,TRMOM2) .or. ERRQ
      ERRQ=COMP8('TRpbl ',npbl*NTM,IM*JM,4,TRABL1,TRABL2) .or. ERRQ
#ifdef TRACERS_WATER
      ERRQ=COMP8('TRW   ',IM,JM,LM*NTM    , TRW1 , TRW2 ) .or. ERRQ
      ERRQ=COMP8Lijp('TRLK  ',2*NTM,IM,JM  , TRLK1, TRLK2) .or. ERRQ
      ERRQ=COMP8Lijp('TRSI  ',LMI*NTM,IM,JM, TRSI1, TRSI2) .or. ERRQ
      ERRQ=COMP8Lijp('TRLI  ',2*NTM,IM,JM  , TRLI1, TRLI2) .or. ERRQ
      ERRQ=COMP8Lijp('TRSOIL',NTM*(2*NGM+3),IM,JM  , TRSOIL1,
     *     TRSOIL2) .or. ERRQ
      ERRQ=COMP8Lijp('TRSNOW',2*NTM*NLSN,IM,JM, TRSN1, TRSN2) .or. ERRQ
#endif
#ifdef TRACERS_SPECIAL_Shindell
      ERRQ=COMP8('yNO3  ',IM,JM,LM    , yNO31   , yNO32   )  .or. ERRQ
      ERRQ=COMP8('pHOx  ',IM,JM,LM    , pHOx1   , pHOx2   )  .or. ERRQ 
      ERRQ=COMP8('pNOx  ',IM,JM,LM    , pNOx1   , pNOx2   )  .or. ERRQ
      ERRQ=COMP8('pOx   ',IM,JM,LM    , pOx1    , pOx2    )  .or. ERRQ
      ERRQ=COMP8('yCH3O2',IM,JM,LM    , yCH3O21 , yCH3O22 )  .or. ERRQ 
      ERRQ=COMP8('yC2O3 ',IM,JM,LM    , yC2O31  , yC2O32  )  .or. ERRQ   
      ERRQ=COMP8('yROR  ',IM,JM,LM    , yROR1   , yROR2   )  .or. ERRQ
      ERRQ=COMP8('yXO2  ',IM,JM,LM    , yXO21   , yXO22   )  .or. ERRQ
      ERRQ=COMP8('yAldeh',IM,JM,LM    ,yAldehyde1,yAldehyde2).or. ERRQ
      ERRQ=COMP8('yXO2N ',IM,JM,LM    , yXO2N1  , yXO2N2  )  .or. ERRQ
      ERRQ=COMP8('yRXPAR',IM,JM,LM    , yRXPAR1 , yRXPAR2 )  .or. ERRQ
      ERRQ=COMP8('corrOx',JM, 4,12    , corrOx1 , corrOx2 )  .or. ERRQ
      ERRQ=COMP8('O3DLJI',40,JM,IM    , O3DLJI1 , O3DLJI2 )  .or. ERRQ
      ERRQ=COMP8('O3clim',40,JM,IM,O3DLJI_clim1,O3DLJI_clim2).or. ERRQ
      DO K1=1,JPPJ
        WRITE(6,*)'K1=',K1
        temp1(:,:,:) = ss1(K1,:,:,:)
        temp2(:,:,:) = ss2(K1,:,:,:)
        ERRQ=COMP8('ss    ',IM,JM,LM,temp1,temp2) .or. ERRQ
      END DO
#endif
#endif

      print*," Diagnostic variables:"
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
      ERRQ=COMP8p('OA    ',IM,JM,KOA    ,OA1,OA2)

#ifdef TRACERS_ON
      ERRQ=COMP8('TAIJLN',IM,JM,LM*NTM    ,TAIJLN1,TAIJLN2) .or. ERRQ
      ERRQ=COMP8('TAIJN ',IM,JM,KTAIJ*NTM ,TAIJN1 ,TAIJN2) .or. ERRQ
      ERRQ=COMP8('TAIJS ',IM,JM,KTAIJS    ,TAIJS1 ,TAIJS2) .or. ERRQ
      ERRQ=COMP8('TAJLN ',JM,LM,KTAJLX*NTM,TAJLN1 ,TAJLN2) .or. ERRQ
      ERRQ=COMP8('TAJLS ',JM,LM,KTAJLS    ,TAJLS1 ,TAJLS2) .or. ERRQ
      ERRQ=COMP8('TCONS ',JM,KTCON,NTM    ,TCON1  ,TCON2 ) .or. ERRQ
#endif
      IF (KOCEAN1.eq.KOCEAN2.and.KOCEAN1.eq.2) THEN ! compare ocn diags
      DAGPOS=1
      ERRQ=COMP8('OIJ   ',IM,JM,KOIJ,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*KOIJ
      ERRQ=COMP8('OIJL  ',IM,JM,LMO*KOIJL,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      DAGPOS=DAGPOS+IM*JM*LMO*KOIJL
      ERRQ=COMP8('OL    ',LMO,KOL,1   ,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
      DAGPOS=DAGPOS+LMO*KOL
      ERRQ=COMP8('OLNST ',LMO,NMST,KOLNST,ODIAG1(DAGPOS),ODIAG2(DAGPOS))
#ifdef TRACERS_WATER
      ERRQ=COMP8('TRICDG',IMIC,JMIC,KTICIJ*NTM,TRICDG1,TRICDG2) .or.
     *     ERRQ
#endif
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
