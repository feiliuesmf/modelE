C**** EE001M12 E001M12 SOMTQ EB357M12
C****
C**** Subroutine EARTH used by new land surface model.
C**** Coded to use new SOC PBL routines.
      SUBROUTINE EARTH (NS,MODDSF,MODD6)
C****
C**** This subroutine calculates surface fluxes of sensible heat,
C**** evaporation, thermal radiation, and momentum drag.
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,edpery
      USE E001M12_COM
      USE GEOM
      USE RADNCB, only : trhr,fsf,cosz1
      USE GHYCOM, only : wbare,wvege,htbare,htvege,snowbv
      USE SLE001
     &    , only : reth,retp,advnc,
     &    NGM,
     &    PR,HTPR,PRS,HTPRS,GW=>W,HT,SNOWD,TP,FICE,GHOUR=>HOUR,
     &    FV,FB,ATRG,ASHG,ALHG,
     &    BETAD=>ABETAD,BETAV=>ABETAV,BETAT=>ABETAT,
     &    BETAP=>ABETAP,BETAB=>ABETAB,BETA=>ABETA,
     &    ACNA,ACNC,
     &    EVAPW=>AEVAPW,EVAPD=>AEVAPD,EVAPB=>AEVAPB,
     &    ARUNS,ARUNU,AERUNS,AERUNU,
     &    DIFS=>ADIFS,EDIFS=>AEDIFS,
     &    AEPC,AEPB,AEPP,AFHG,AF0DT,AF1DT,ZW,TBCS,
     &    QM1,Q1,QS,
     &    PRES,RHO,TSPASS=>TS,VSM,CH,CD,SNHT,SRHT,TRHT,ZS,
     &    ZMIXE=>Z1,CN=>CDN,P1,PBLP=>PPBL,
     &    TGPASS=>TG,TKPASS=>T1,VGM=>VG,EDDY
      USE CLD01_COM_E001, only : PREC,TPREC
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs,tsavg,qsavg
      USE SOCPBL, only : zgs
      USE DAGCOM  !, only : aijg,aij,tsfrez,tdiurn,bj,areg,adaily,jreg
      USE DYNAMICS, only : pmid,pk,pek,pedn

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NS,MODDSF,MODD6
      INTEGER I,J,L,KR,JR,ITYPE,IM1,IMAX,IHOUR
      REAL*8 SHDT,QSATS,EVAP,EVHDT,TG2AV,ACE2AV,TRHDT,CDN,TK,RCDMWS
     *     ,RCDHWS,DHGS,CDQ,CDM,CDH,ELHX,TG,SRHEAT,TG1,PTYPE,QSATSS
     *     ,EVAPS,PPBLS,EDS1S,DBLS,DGSS,CDHS,CDMS,QGS,SRHDTS,SHDTS
     *     ,EVHDTS,TG1S,DXYPJ,TRHEAT,WTR2AV,WFC1,WGS,VGS,UGS,TRHDTS
     *     ,RTAUVS,RTAUUS,RTAUS,TAUS,WSS,TSS,QSS,USS,VSS,RHOSRF,RMBYA
     *     ,TFS,TH1,THV1,P1K,PSK,TS,PS,PIJ,PSOIL,PEARTH,WARMER,BRUN0
     *     ,BERUN0,BDIFS,BEDIFS,BTS,BEVHDT,BRUNU,BERUNU,BSHDT,BTRHDT
     *     ,TIMEZ,SPRING,ZS1CO,RVX,DTSRCE,DTCNDS,QLH,PM,TM,QSAT,STBO

      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1,DU1,DV1
      COMMON /WORK1d/DTH1,DQ1
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT
      COMMON/WORK2/UT,VT,DU1,DV1
      REAL*8, DIMENSION(IM,JM,4) :: E0,E1,EVAPOR,TGRND
      COMMON/WORK3/E0,E1,EVAPOR,TGRND
      REAL*8, DIMENSION(IM,JM,3) :: GDEEP
      COMMON/oldDAG/GDEEP
C****
      REAL*8, DIMENSION(IM,JM) :: PRCSS
      COMMON /WORKLS/PRCSS
C**** Interface to PBL
      REAL*8 ZS1,TGV,TKV,QG,HEMI,DTSURF,US,VS,WS,TSV,QSRF,PSI,DBL,EDVISC
     *     ,EDS1,PPBL,UG,VG,WG,ZMIX
      LOGICAL POLE
      COMMON /PBLPAR/ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE

      COMMON /PBLOUT/US,VS,WS,TSV,QSRF,PSI,DBL,EDVISC,EDS1,
     2               PPBL,UG,VG,WG,ZMIX

      DATA STBO/.5672573E-7/

      QSAT(TM,PM,QLH)=3.797915*EXP(QLH*(7.93252D-6-2.166847D-3/TM))/PM
C****
C**** FEARTH    SOIL COVERED LAND FRACTION (1)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        2  EARTH SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)
C****        7  OCEAN ICE TEMPERATURE OF SECOND LAYER (C)
C****        8  currently unused
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)
C****
C**** GHDATA
C****      1-6  WATER OF BARE SOIL LAYER 1-6 (M)
C****        7  WATER OF VEGETATION CANOPY (M)
C****     8-13  WATER OF VEGETATED SOIL LAYER 1-6 (M)
C****       14  BARE SOIL LAYER 0 IS UNUSED
C****    15-20  HEAT CONTENT OF BARE SOIL LAYER 1-6 (J M-2)
C****       21  HEAT CONTENT OF VEGETATION CANOPY (J M-2)
C****    22-27  HEAT CONTENT OF VEGETATED SOIL LAYER 1-6 (J M-2)
C****       28  SNOW DEPTH OVER BARE SOIL (M)
C****       29  SNOW DEPTH OVER VEGETATED SOIL (M)
C****

      DTSURF=NDYN*DT/NSURF
      DTCNDS=NDYN*DT
      DTSRCE=DT*NDYN
      RVX=0.
      ZS1CO=.5*DSIG(1)*RGAS/GRAV

         SPRING=-1.
         IF((JDAY.GE.32).AND.(JDAY.LE.212)) SPRING=1.
         IHOUR=1.5+TOFDAY
C        COSDAY=COS(TWOPI*JDAY/EDPERY)
C        SINDAY=SIN(TWOPI*JDAY/EDPERY)
C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NSURF TIMES EVERY HOUR
C****
         TIMEZ=JDAY+(TOFDAY+(NS-1.)/NSURF)/24.
         IF(JDAY.LE.31) TIMEZ=TIMEZ+365.
         GHOUR=TAU+(NS-1.)/NSURF
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO 7000 J=1,JM
      HEMI=1.
      IF(J.LE.JM/2) HEMI=-1.
      IMAX=IMAXJ(J)
      POLE=.FALSE.
C**** CONDITIONS AT THE SOUTH POLE
      IF(J.EQ.1) THEN
         POLE = .TRUE.
      ENDIF
C**** CONDITIONS AT THE NORTH POLE
      IF(J.EQ.JM) THEN
         POLE=.TRUE.
      ENDIF
C**** ZERO OUT SURFACE DIAGNOSTICS WHICH WILL BE SUMMED OVER LONGITUDE
  100    BTRHDT=0.
         BSHDT=0.
         BEVHDT=0.
         BTS=0.
         BRUN0=0.
         BERUN0=0.
         BDIFS=0.
         BEDIFS=0.
         BERUNU=0.
         BRUNU=0.
c         JEQ=1+JM/2
         IF(J.LT.JEQ) WARMER=-SPRING
         IF(J.GE.JEQ) WARMER=SPRING
      IM1=IM
      DO 6000 I=1,IMAX
C****
C**** DETERMINE SURFACE CONDITIONS
C****
c     PLAND=FLAND(I,J)
c     PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      PSOIL=PEARTH
      PIJ=P(I,J)
      PS=PEDN(1,I,J)
      PSK=PEK(1,I,J)
      P1=PMID(1,I,J)
      P1K=PK(1,I,J)
      TH1=T(I,J,1)
      Q1=Q(I,J,1)
      THV1=TH1*(1.+Q1*RVX)
      TKV=THV1*PSK
         TFS=TF*PSOIL
      RMBYA=100.*PIJ*DSIG(1)/GRAV
      QM1=Q1*RMBYA
      RHOSRF=100.*PS/(RGAS*TKV)
C**** ZERO OUT QUANTITIES TO BE SUMMED OVER SURFACE TYPES
      USS=0.
      VSS=0.
      WSS=0.
      TSS=0.
      QSS=0.
      TAUS=0.
         RTAUS=0.
         RTAUUS=0.
         RTAUVS=0.
C        SINAPS=0.
C        COSAPS=0.
         JR=JREG(I,J)
         DXYPJ=DXYP(J)
         TG1S=0.
         QGS=0.
         SRHDTS=0.
         TRHDTS=0.
         SHDTS=0.
         EVHDTS=0.
         UGS=0.
         VGS=0.
         WGS=0.
C        USRS=0.
         CDMS=0.
         CDHS=0.
         DGSS=0.
         EDS1S=0.
         PPBLS=0.
         DBLS=0.
         EVAPS=0.
         QSATSS=0.
C**** New quantities to be zeroed out over ground timesteps
         ARUNS=0.
         ARUNU=0.
         AERUNS=0.
         AERUNU=0.
         DIFS=0.
         EDIFS=0.
      EVAPW=0.
      EVAPD=0.
      EVAPB=0.
         ALHG=0.
         AEPC=0.
         AEPB=0.
         AFHG=0.
         ATRG=0.
      ASHG=0.
         AF0DT=0.
         AF1DT=0.
C****
C**** EARTH
C****
 2400 IF (PEARTH.LE.0.) then
        ipbl(i,j,4)=0
        go to 5000
      endif
C     NGRNDZ=NGRND HAS TO BE 1
c      ZGS=10.
      ITYPE=4
      PTYPE=PEARTH
      PR=PREC(I,J)/(DTCNDS*RHOW)
      PRS=PRCSS(I,J)/(DTCNDS*RHOW)
      HTPR=0.
      IF(TPREC(I,J).LT.0.) HTPR=-LHM*PREC(I,J)/DTCNDS
c      DO 2410 L=1,4*NGM+5
c 2410 GW(L,1)=GHDATA(I,J,L)
      GW(1:NGM,1) = WBARE(I,J,1:NGM)
      GW(0:NGM,2) = WVEGE(I,J,0:NGM)
      HT(0:NGM,1) = HTBARE(I,J,0:NGM)
      HT(0:NGM,2) = HTVEGE(I,J,0:NGM)
      SNOWD(1:2) = SNOWBV(I,J,1:2)
      CALL GHINIJ (I,J,WFC1)
      CALL RETH
      CALL RETP
C     CALL HYDRA
C??   SNOW = SNOWD(1)*FB + SNOWD(2)*FV
      TG1=TBCS
      SRHEAT=FSF(I,J,ITYPE)*COSZ1(I,J)
         SRHDTS=SRHDTS+SRHEAT*DTSURF*PTYPE
C****
C**** BOUNDARY LAYER INTERACTION
C****
 3000 ZS1=ZS1CO*TKV*PIJ/PS
      P1=PMID(1,I,J)    ! SIG(1)*PIJ+PTOP
C**** LOOP OVER GROUND TIME STEPS
      TG=TG1+TF
      ELHX=LHE
      IF(TG1.LT.0.)  ELHX=LHS
      QG=QSAT(TG,PS,ELHX)
      TGV=TG*(1.+QG*RVX)
C***********************************************************************
C***
      CALL PBL(I,J,ITYPE,PTYPE)
C***
      CDM = cmgs(i,j,itype)
      CDH = chgs(i,j,itype)
      CDQ = cqgs(i,j,itype)
C***********************************************************************
C**** CALCULATE QS
         dhgs=(zmix-zgs)*cdh*ws
C****    dqgs=(zmix-zgs)*cdq*ws
      QS=QSRF
      TS=TSV/(1.+QS*RVX)
      TK=TKV/(1.+QS*RVX)
C**** CALCULATE RHOSRF*CDM*WS
      RCDMWS=CDM*WS*RHOSRF
      RCDHWS=CDH*WS*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)
      SNHT=-SHA*RCDHWS*(TS-TG)
      TRHEAT=TRHR(I,J,1)
C **********************************************************************
C *****
C  Define extra variables to be passed in SURFC:
      PRES  =PS
      RHO   =RHOSRF
      VSM   =WS
      CH    =CDH
      CD    =CDM
      SRHT  =SRHEAT
      TRHT  =TRHEAT
      ZS    =ZGS
      ZMIXE =ZMIX
      PBLP  =PPBL
      VGM   =WG
      EDDY  =EDS1
      CN    =CDN
      TGPASS=TGV
      TSPASS=TSV
      TKPASS=TKV
C **********************************************************************
C *****
C**** Calculate ground fluxes
C     CALL QSBAL
      CALL ADVNC
      TG1=TBCS
C      DO 3410 L=1,4*NGM+5
C 3410 GHDATA(I,J,L)=GW(L,1)
      WBARE(I,J,1:NGM) = GW(1:NGM,1)
      WVEGE(I,J,0:NGM) = GW(0:NGM,2)
      HTBARE(I,J,0:NGM) = HT(0:NGM,1)
      HTVEGE(I,J,0:NGM) = HT(0:NGM,2)
      SNOWBV(I,J,1:2) = SNOWD(1:2)
      AIJG(I,J, 5)=AIJG(I,J, 5)+BETAB/NSURF
      AIJG(I,J, 6)=AIJG(I,J, 6)+BETAP/NSURF
      AIJG(I,J,11)=AIJG(I,J,11)+BETA/NSURF
      AIJG(I,J,12)=AIJG(I,J,12)+ACNA/NSURF
      AIJG(I,J,13)=AIJG(I,J,13)+ACNC/NSURF
      AIJG(I,J,14)=AIJG(I,J,14)+AEPP
      AIJG(I,J,15)=AIJG(I,J,15)+TP(1,1)
      AIJG(I,J,16)=AIJG(I,J,16)+TP(2,1)
      AIJG(I,J,17)=AIJG(I,J,17)+TP(3,1)
      AIJG(I,J,18)=AIJG(I,J,18)+EVAPB
      AIJG(I,J,19)=AIJG(I,J,19)+EVAPD
      AIJG(I,J,20)=AIJG(I,J,20)+EVAPW
      AIJG(I,J,21)=AIJG(I,J,21)+TP(0,2)
      AIJG(I,J,22)=AIJG(I,J,22)+TP(1,2)
      AIJG(I,J,23)=AIJG(I,J,23)+TP(2,2)
      AIJG(I,J,24)=AIJG(I,J,24)+TP(3,2)
      AIJG(I,J,25)=AIJG(I,J,25)+FB*ZW(1)+FV*ZW(2)
      AIJG(I,J,26)=AIJG(I,J,26)+BETAV/NSURF
      AIJG(I,J,27)=AIJG(I,J,27)+BETAT/NSURF
          TRHDT=TRHEAT*DTSURF-ATRG
C     FOR RADIATION FIND COMPOSITE VALUES  GDATA 2 4 5 6
C           FOR DIAGNOSTIC PURPOSES ALSO COMPUTE GDEEP 1 2 3
      GDATA(I,J,2)=1000.*(SNOWD(1)*FB+SNOWD(2)*FV)
      GDATA(I,J,4)=TG1
      GDATA(I,J,5)=1000.*( FB*GW(1,1)*(1.-FICE(1,1)) +
     +  FV*(GW(1,2)*(1.-FICE(1,2))+GW(0,2)*(1.-FICE(0,2))) )
      GDATA(I,J,6)=1000.*(FB*(GW(1,1)*FICE(1,1)-SNOWD(1)) +
     +  FV*(GW(1,2)*FICE(1,2)+GW(0,2)*FICE(0,2)-SNOWD(2)) )
            CALL RETP2 (TG2AV,WTR2AV,ACE2AV)
            GDEEP(I,J,1)=TG2AV
            GDEEP(I,J,2)=WTR2AV
            GDEEP(I,J,3)=ACE2AV
C**** CALCULATE FLUXES USING IMPLICIT TIME STEP FOR NON-OCEAN POINTS
      DU1(I,J)=DU1(I,J)+PTYPE*DTSURF*RCDMWS*US/RMBYA
 3600 DV1(I,J)=DV1(I,J)+PTYPE*DTSURF*RCDMWS*VS/RMBYA
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      EVAP=EVAPW+EVAPD+EVAPB
         EVAPOR(I,J,4)=EVAPOR(I,J,4)+EVAP
         EVHDT=-ALHG
      SHDT=-ASHG
      DTH1(I,J)=DTH1(I,J)-SHDT*PTYPE/(SHA*RMBYA*P1K)
      DQ1(I,J) =DQ1(I,J)+EVAP*PTYPE/RMBYA
      USS=USS+US*PTYPE
      VSS=VSS+VS*PTYPE
      WSS=WSS+WS*PTYPE
      TSS=TSS+TS*PTYPE
      QSS=QSS+QS*PTYPE
      QSAVG(I,J)=QSAVG(I,J)+QSS
      TAUS=TAUS+CDM*WS*WS*PTYPE
         QSATS=QSAT(TS,PS,ELHX)
         RTAUS=RTAUS+RCDMWS*WS*PTYPE
         RTAUUS=RTAUUS+RCDMWS*US*PTYPE
         RTAUVS=RTAUVS+RCDMWS*VS*PTYPE
C        SINAPS=SINAPS+SINAP*PTYPE
C        COSAPS=COSAPS+COSAP*PTYPE
         TG1S=TG1S+TG1*PTYPE
         QGS=QGS+QG*PTYPE
         TRHDTS=TRHDTS+TRHDT*PTYPE
         SHDTS=SHDTS+SHDT*PTYPE
         EVHDTS=EVHDTS+EVHDT*PTYPE
         UGS=UGS+UG*PTYPE
         VGS=VGS+VG*PTYPE
         WGS=WGS+WG*PTYPE
C        USRS=USRS+USR*PTYPE
         CDMS=CDMS+CDM*PTYPE
         CDHS=CDHS+CDH*PTYPE
         DGSS=DGSS+DHGS*PTYPE
         EDS1S=EDS1S+EDS1*PTYPE
         PPBLS=PPBLS+PPBL*PTYPE
         DBLS=DBLS+DBL*PTYPE
         EVAPS=EVAPS+EVAP*PTYPE
         QSATSS=QSATSS+QSATS*PTYPE
C****
C**** EARTH
C****
 4600    BSHDT=BSHDT+SHDT*PEARTH
         BEVHDT=BEVHDT+EVHDT*PEARTH
         BTRHDT=BTRHDT+TRHDT*PEARTH
         BTS=BTS+(TS-TF)*PEARTH
         BRUN0=BRUN0+ARUNS*PEARTH
         BERUN0=BERUN0+AERUNS*PEARTH
         BRUNU=BRUNU+ARUNU*PEARTH
         BERUNU=BERUNU+AERUNU*PEARTH
         AIJ(I,J,IJ_RUNE)=AIJ(I,J,IJ_RUNE)+ARUNS
         AIJ(I,J,IJ_ARUNU)=AIJ(I,J,IJ_ARUNU)+ARUNU
         AIJ(I,J,IJ_PEVAP)=AIJ(I,J,IJ_PEVAP)+(AEPC+AEPB)
         BDIFS=BDIFS+DIFS*PEARTH
         BEDIFS=BEDIFS+EDIFS*PEARTH
         E0(I,J,4)=E0(I,J,4)+AF0DT
         E1(I,J,4)=E1(I,J,4)+AF1DT
         IF(WARMER.LT.0.) GO TO 4610
         IF(TS.LT.TF) TSFREZ(I,J,1)=TIMEZ
         TSFREZ(I,J,2)=TIMEZ
         GO TO 4620
 4610    IF(TSFREZ(I,J,2)+.03.LT.TIMEZ) GO TO 4620
         IF(TS.GE.TF) TSFREZ(I,J,2)=TIMEZ
 4620    IF(TG1.LT.TDIURN(I,J,1)) TDIURN(I,J,1)=TG1
         IF(TG1.GT.TDIURN(I,J,2)) TDIURN(I,J,2)=TG1
         IF(TS.LT.TDIURN(I,J,3)) TDIURN(I,J,3)=TS
         IF(TS.GT.TDIURN(I,J,4)) TDIURN(I,J,4)=TS
C**** NON-OCEAN POINTS WHICH ARE NOT MELTING OR FREEZING WATER USE
C****   IMPLICIT TIME STEPS
C****
C**** UPDATE SURFACE AND FIRST LAYER QUANTITIES
C****
 5000 CONTINUE
         TDIURN(I,J,5)=TDIURN(I,J,5)+(TSAVG(I,J)-TF)
         IF(TSAVG(I,J).GT.TDIURN(I,J,6)) TDIURN(I,J,6)=TSAVG(I,J)
C        IF(TSAVG(I,J).GT.AIJ(I,J,IJ_TMAX)) AIJ(I,J,IJ_TMAX)=TSAVG(I,J)
C        IF(TSAVG(I,J).LT.AIJ(I,J,IJ_TMIN)) AIJ(I,J,IJ_TMIN)=TSAVG(I,J)
      IF(PEARTH.LE.0.)  GO TO 6000
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
         IF(JR.EQ.24) GO TO 5700
         AREG(JR,9)=AREG(JR,9)+TRHDTS*DXYPJ
         AREG(JR,13)=AREG(JR,13)+SHDTS*DXYPJ
         AREG(JR,14)=AREG(JR,14)+EVHDTS*DXYPJ
         AREG(JR,19)=AREG(JR,19)+EVAPS*DXYPJ
         AREG(JR,40)=AREG(JR,40)+AERUNS*PEARTH*DXYPJ
         AREG(JR,45)=AREG(JR,45)+DIFS*PEARTH*DXYPJ
         AREG(JR,47)=AREG(JR,47)+ARUNU*PEARTH*DXYPJ
         AREG(JR,48)=AREG(JR,48)+AERUNU*PEARTH*DXYPJ
         AREG(JR,54)=AREG(JR,54)+ARUNS*PEARTH*DXYPJ
         IF(MODDSF.NE.0) GO TO 5700
         AREG(JR,23)=AREG(JR,23)+(TSS-TFS)*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
 5700    AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDTS
         AIJ(I,J,IJ_BETA)=AIJ(I,J,IJ_BETA)+BETAD/NSURF
         IF(MODRD.EQ.0) AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)+TRHDTS
     *        /DTSRCE
         AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHDTS+TRHDTS)
         AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHDTS+TRHDTS+SHDTS+EVHDTS)
         IF(MODDSF.NE.0) GO TO 5800
         AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WSS                ! added 3/3/95 -rar-
         AIJ(I,J,IJ_TS)=AIJ(I,J,IJ_TS)+(TSS-TFS)
         AIJ(I,J,IJ_US)=AIJ(I,J,IJ_US)+USS
         AIJ(I,J,IJ_VS)=AIJ(I,J,IJ_VS)+VSS
         AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RTAUS
         AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RTAUUS
         AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RTAUVS
         AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QSS
CHYD     AIJ(I,J,IJ_ARUNU)=AIJ(I,J,IJ_ARUNU)
CHYD *  +   (40.6*PSOIL+.72*(2.*(TSS-TFS)-(QSATSS-QSS)*LHE/SHA))
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAG6
 5800    IF(MODD6.NE.0) GO TO 6000
         DO KR=1,4
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(IHOUR,12,KR)=ADAILY(IHOUR,12,KR)+TSS
               ADAILY(IHOUR,13,KR)=ADAILY(IHOUR,13,KR)+(TG1S+TFS)
               ADAILY(IHOUR,19,KR)=ADAILY(IHOUR,19,KR)+QSS
               ADAILY(IHOUR,20,KR)=ADAILY(IHOUR,20,KR)+QGS
               ADAILY(IHOUR,28,KR)=ADAILY(IHOUR,28,KR)+SRHDTS
               ADAILY(IHOUR,29,KR)=ADAILY(IHOUR,29,KR)+TRHDTS
               ADAILY(IHOUR,30,KR)=ADAILY(IHOUR,30,KR)+SHDTS
               ADAILY(IHOUR,31,KR)=ADAILY(IHOUR,31,KR)+EVHDTS
               ADAILY(IHOUR,32,KR)=ADAILY(IHOUR,32,KR)
     *              +SRHDTS+TRHDTS+SHDTS+EVHDTS
               ADAILY(IHOUR,33,KR)=ADAILY(IHOUR,33,KR)+UGS
               ADAILY(IHOUR,34,KR)=ADAILY(IHOUR,34,KR)+VGS
               ADAILY(IHOUR,35,KR)=ADAILY(IHOUR,35,KR)+WGS
               ADAILY(IHOUR,36,KR)=ADAILY(IHOUR,36,KR)+USS
               ADAILY(IHOUR,37,KR)=ADAILY(IHOUR,37,KR)+VSS
               ADAILY(IHOUR,38,KR)=ADAILY(IHOUR,38,KR)+WSS
               ADAILY(IHOUR,42,KR)=ADAILY(IHOUR,42,KR)+CDMS
               ADAILY(IHOUR,43,KR)=ADAILY(IHOUR,43,KR)+CDHS
               ADAILY(IHOUR,44,KR)=ADAILY(IHOUR,44,KR)+DGSS
               ADAILY(IHOUR,45,KR)=ADAILY(IHOUR,45,KR)+EDS1S
               ADAILY(IHOUR,46,KR)=ADAILY(IHOUR,46,KR)+DBLS
               ADAILY(IHOUR,50,KR)=ADAILY(IHOUR,50,KR)+EVAPS
            END IF
         END DO
 6000 IM1=I
C**** QUANTITIES ACCUMULATED FOR SURFACE TYPE TABLES IN DIAGJ
         BJ(J,9)=BJ(J,9)+BTRHDT
         BJ(J,13)=BJ(J,13)+BSHDT
         BJ(J,14)=BJ(J,14)+BEVHDT
         BJ(J,40)=BJ(J,40)+BERUN0
         BJ(J,41)=BJ(J,41)+BEDIFS
         BJ(J,45)=BJ(J,45)+BDIFS
         BJ(J,47)=BJ(J,47)+BRUNU
         BJ(J,48)=BJ(J,48)+BERUNU
         BJ(J,54)=BJ(J,54)+BRUN0
         IF(MODDSF.NE.0) GO TO 7000
         BJ(J,23)=BJ(J,23)+BTS
 7000 CONTINUE
      RETURN
      END

      SUBROUTINE GHINIT (DTSURF,redoGH)
C**** Modifications needed for split of bare soils into 2 types
      USE CONSTANT, only : twopix=>twopi,rhow,edpery
      USE E001M12_COM, only : im,jm,fearth,vdata,gdata,tau,jeq
      USE GHYCOM
      USE SLE001, sinday=>sint,cosday=>cost
      USE FILEMANAGER
      IMPLICIT NONE

      REAL*8 DTSURF
      INTEGER iu_SOIL
      INTEGER JDAY
      REAL*8 SNOWDP,WTR1,WTR2,ACE1,ACE2,TG1,TG2
      LOGICAL redoGH
C****             TUNDR GRASS SHRUB TREES DECID EVRGR RAINF CROPS
C****
C**** LADAY(veg type, lat belt) = day of peak LAI
C OLD PEAK LAI:  2ND LINE IS FOR LATITUDES < 23.5 DEG
C****    1  temperate latitudes
C****    2  non-temperate latitudes
C     DATA  LADAY/ 196,  196,  196,  196,  196,  196,  105,  196/
C     DATA  LADAY/ 196,  288,  288,  288,  288,  196,  105,  288/
C****
C**** CONTENTS OF ALA(K,I,J),  LAI coefficients
C****   1  AVERAGE LEAF AREA INDEX
C****   2  REAL AMPLITUDE OF LEAF AREA INDEX
C****   3  IMAGINARY AMPLITUDE OF LEAF AREA INDEX
C****
C**** CONTENTS OF ACS(K,I,J),  CS coefficients
C****   1  AVERAGE STOMATAL CONDUCTANCE
C****   2  REAL AMPLITUDE OF STOMATAL CONDUCTANCE
C****   3  IMAGINARY AMPLITUDE OF STOMATAL CONDUCTANCE
C****
C**** CONTENTS OF SDATA(I,J,K):
C****       1 -   NGM   DZ(NGM)
C****   NGM+1 - 6*NGM   Q(IS,NGM)
C**** 6*NGM+1 - 11*NGM   QK(IS,NGM)
C**** 11*NGM+1           SL
C
C READ SOILS PARAMETERS
      call getunit("SOIL",iu_SOIL,.TRUE.)
      CALL DREAD (iu_SOIL,DZ_IJ,IM*JM*(11*NGM+1),DZ_IJ)
      CLOSE (iu_SOIL)
C
      ONE=1.
   10 CONTINUE
C****
C**** INITIALIZE CONSTANTS
C****
C**** Time step for ground hydrology
      DT=DTSURF
C**** UNITS ARE MKS
C**** WATER QUANTITIES ARE DENSITY TIMES USUAL VALUES IN MKS
C**** TO GET VOLUMETRIC UNITS
C**** 1M WATER = 1000 KG M-2; 1M3 WATER = 1000 KG
C FSN IS THE HEAT OF FUSION
      FSN=3.34 E+8
C ELH IS THE HEAT OF VAPORIZATION
      ELH=2.50 E+9
C THE SH'S ARE THE SPECIFIC HEAT CAPACATIES
      SHW=4.185 E+6
      SHI=2.060 E+6
      SHA=1003.4965
      SHV=1911.
C THE ALAM'S ARE THE HEAT CONDUCTIVITIES
      ALAMW=.573345
      ALAMI=2.1762
      ALAMA=.025
      ALAMSN=0.088
      ALAMBR=2.9
      ALAMS(1)=8.8
      ALAMS(2)=2.9
      ALAMS(3)=2.9
      ALAMS(4)=.25
C HW IS THE WILTING POINT IN METERS
      HW=-100
C TFRZ IS 0 C IN K
c      TFRZ=273.16
C ZHTB IS DEPTH FOR COMBINING HEAT LAYERS FOR STABILITY
      IF(Q(4,1).LT..01)THEN
      ZHTB=6.
      ELSE
      ZHTB=6.
      ENDIF
C SPGSN IS THE SPECIFIG GRAVITY OF SNOW
      SPGSN=.1
C
C****
C**** Initialize global arrays  ALA, ACS, AFB, AFR
C****
      TWOPI=6.283185

      ALA(:,:,:)=0.
      ACS(:,:,:)=0.
      AFB(:,:)=0.
      AFR(:,:,:)=0.
      DO 220 J=1,JM
      DO 220 I=1,IM
  220 ACS(1,I,J)=.01
      DO 400 J=1,JM
      DO 400 I=1,IM
      PEARTH=FEARTH(I,J)
      AFB(I,J)=VDATA(I,J,1)+VDATA(I,J,10)
      IF(AFB(I,J).GT..999) AFB(I,J)=1.
      IF(PEARTH.LE.0..OR.AFB(I,J).GE.1.) GO TO 400
C**** CALCULATE LAI, CS COEFFICICENTS
      SFV=0.
      SLA0=0.
      SLRE=0.
      SLIM=0.
      SCS0=0.
      SCSRE=0.
      SCSIM=0.
      SVH=0.
      DO 250 IV=1,8
      PHASE=TWOPI*LADAY(IV)/365.
      IF(J.LT.JEQ) PHASE=PHASE+TWOPI/2.
      FV=VDATA(I,J,IV+1)
      SFV=SFV+FV
      SVH=SVH+FV*VHGHT(IV)
      DIF=(ALAMAX(IV) - ALAMIN(IV))
      SLA0=SLA0+FV*(ALAMAX(IV) + ALAMIN(IV))
      SLRE=SLRE+FV*DIF*COS(PHASE)
      SLIM=SLIM+FV*DIF*SIN(PHASE)
      SCS0=SCS0+FV*(ALAMAX(IV) + ALAMIN(IV))/RSAR(IV)
      SCSRE=SCSRE+FV*DIF*COS(PHASE)/RSAR(IV)
  250 SCSIM=SCSIM+FV*DIF*SIN(PHASE)/RSAR(IV)
      ALA(1,I,J)=.5/SFV*SLA0
      ALA(2,I,J)=.5/SFV*SLRE
      ALA(3,I,J)=.5/SFV*SLIM
      ACS(1,I,J)=.5/SFV*SCS0
      ACS(2,I,J)=.5/SFV*SCSRE
      ACS(3,I,J)=.5/SFV*SCSIM
      AVH(I,J)=SVH/SFV
C**** CALCULATE ROOT FRACTION AFR AVERAGED OVER VEGETATION TYPES
      DO 310 N=1,NGM
      DZ(N)=DZ_IJ(I,J,N)
      IF(DZ(N).LE.0.) GO TO 320
  310 CONTINUE
  320 N=N-1
      DO 350 IV=1,8
      FV=VDATA(I,J,IV+1)
      Z=0.
      FRUP=0.
      DO 350 L=1,N
      Z=Z+DZ(L)
      FRDN=AROOT(IV)*Z**BROOT(IV)
      FRDN=MIN(FRDN,ONE)
      IF(L.EQ.N)FRDN=1.
      AFR(I,J,L) = AFR(I,J,L) + FV*(FRDN-FRUP)
  350 FRUP=FRDN
      DO 370 L=1,N
  370 AFR(I,J,L) = AFR(I,J,L)/(1.-AFB(I,J))
  400 CONTINUE
C****
      SDSTNC=100.
      PRINT *,'SDSTNC:',SDSTNC
      C1=90.
      PRINT *,'C1:',C1
      PRFR=.1
      PRINT *,'PRFR:',PRFR
      CALL HL0
C****
C code transplanted from subroutine INPUT
C**** Recompute GHDATA if necessary (new soils data)
      IF (redoGH) THEN
        JDAY=1+MOD(NINT(TAU/24.),365)
        COSDAY=COS(TWOPIx/EDPERY*JDAY)
        SINDAY=SIN(TWOPIx/EDPERY*JDAY)

        DO 930 J=1,JM
        DO 930 I=1,IM
        PEARTH=FEARTH(I,J)
        IF(PEARTH.LE.0.) THEN

          WBARE(I,J,:)=0.
          WVEGE(I,J,:)=0.
          HTBARE(I,J,:)=0.
          HTVEGE(I,J,:)=0.
          SNOWBV(I,J,:)=0.
        ELSE
C****     COMPUTE SOIL HEAT CAPACITY AND GROUND WATER SATURATION GWS
          CALL GHINIJ (I,J,WFC1)
C****     FILL IN SOILS COMMON BLOCKS
          SNOWDP=GDATA(I,J,2)/RHOW
          WTR1=GDATA(I,J,5)
          ACE1=GDATA(I,J,6)
          TG1 =GDATA(I,J,4)
          WTR2=GDATA(I,J,9)
          ACE2=GDATA(I,J,10)
          TG2 =GDATA(I,J,8)
          CALL GHINHT (SNOWDP, TG1,TG2, WTR1,WTR2, ACE1,ACE2)
C****     COPY SOILS PROGNOSTIC QUANTITIES TO EXTENDED GHDATA
          WBARE(I,J,1:NGM) = W(1:NGM,1)
          WVEGE(I,J,0:NGM) = W(0:NGM,2)
          HTBARE(I,J,0:NGM) = HT(0:NGM,1)
          HTVEGE(I,J,0:NGM) = HT(0:NGM,2)
          SNOWBV(I,J,1:2) = SNOWD(1:2)
        END IF
  930   CONTINUE
      END IF

      RETURN
      END SUBROUTINE GHINIT

      SUBROUTINE GHINIJ (I0,J0, WFCAP)
C**** INPUT:
C**** AVH(I,J) - ARRAY OF VEGETATION HEIGHTS
C**** SPGSN - SPECIFIC GRAVITY OF SNOW
C**** OUTPUT:
C**** VH - VEGETATION HEIGHT
C**** SNOWM - SNOW MASKING DEPTH
C**** WFCAP - WATER FIELD CAPACITY OF TOP SOIL LAYER, M
C****
      USE GHYCOM
      USE SLE001
      IMPLICIT NONE
      INTEGER I0,J0
      REAL*8 WFCAP

      ONE=1.
      ID=I0*100+J0
C**** SET UP LAYERS
      DZ(1:NGM)=DZ_IJ(I0,J0,1:NGM)
      Q(1:IMT,1:NGM)=Q_IJ(I0,J0,1:IMT,1:NGM)
      QK(1:IMT,1:NGM)=QK_IJ(I0,J0,1:IMT,1:NGM)
      SL=SL_IJ(I0,J0)
      DO 20 N=1,NGM
      IF(DZ(N).LE.0.) GO TO 21
   20 CONTINUE
   21 N=N-1
      IF(N.LE.0) THEN
         WRITE (99,*) 'GHINIJ:  N <= 0:  I,J,N=',I0,J0,N,(DZ(K),K=1,43)
         STOP
      END IF
C**** CALCULATE THE BOUNDARIES, BASED ON THE THICKNESSES.
      ZB(1)=0.
      DO 30 L=1,N
   30 ZB(L+1)=ZB(L)-DZ(L)
C**** CALCULATE THE LAYER CENTERS, BASED ON THE BOUNDARIES.
      DO 40 L=1,N
   40 ZC(L)=.5*(ZB(L)+ZB(L+1))
C**** FR: ROOT FRACTION IN LAYER L  (1=FR(1)+FR(2)+...+FR(N))
      DO 45 L=1,N
      FR(L)=AFR(I0,J0,L)
   45 CONTINUE
C**** VH: VEGETATION HEIGHT
      VH=AVH(I0,J0)
      SNOWM=VH*SPGSN
C**** FB,FV: BARE, VEGETATED FRACTION (1=FB+FV)
      FB=AFB(I0,J0)
      FV=1.-FB
C**** ALAI: LEAF AREA INDEX
      ALAI=ALA(1,I0,J0)+COST*ALA(2,I0,J0)+SINT*ALA(3,I0,J0)
      ALAI=MAX(ALAI,ONE)
      ALAIC=5.0
      ALAIE=ALAIC*(1.-EXP(-ALAI/ALAIC))
C**** RS: MINIMUM STOMATAL RESISTANCE
      RS=ALAI/(ACS(1,I0,J0)+COST*ACS(2,I0,J0)+SINT*ACS(3,I0,J0))
C???  CNC=ALAI/RS   REDEFINED BEFORE BEING USED (QSBAL,COND)
C
CW    WRITE(6,*)'N=',N,'  R=',R
CW    WRITE(6,91)
CW 91 FORMAT(1X,5X,'ZB',5X,'ZC',5X,'DZ'/1X,21('-'))
CW    DO 95 L=1,N
CW 95 WRITE(6,100)ZB(L),ZC(L),DZ(L)
CW    WRITE(6,100)ZB(N+1)
CW100 FORMAT(1X,3F7.3)
CW    WRITE(6,*)
C****
      DO 60 IBV=1,2
      DO 60 L=1,N
      THETS(L,IBV)=0.
      THETM(L,IBV)=0.
      DO 50 I=1,IMT-1
      THETS(L,IBV)=THETS(L,IBV)+Q(I,L)*THM(0,I)
      THETM(L,IBV)=THETM(L,IBV)+Q(I,L)*THM(NTH,I)
   50 CONTINUE
      WS(L,IBV)=THETS(L,IBV)*DZ(L)
   60 CONTINUE
      WS(0,2)=.0001*ALAI
      WFCAP=FB*WS(1,1)+FV*(WS(0,2)+WS(1,2))
C****
      CALL XKLH0
C****
      DO 90 IBV=1,2
      DO 90 L=1,N
      SHC(L,IBV)=0.
      DO 80 I=1,IMT
      SHC(L,IBV)=SHC(L,IBV)+Q(I,L)*SHCAP(I)
   80 CONTINUE
      SHC(L,IBV)=(1.-THETS(L,IBV))*SHC(L,IBV)*DZ(L)
   90 CONTINUE
C****
C SHC(0,2) IS THE HEAT CAPACITY OF THE CANOPY
      AA=ALA(1,I0,J0)
      SHC(0,2)=(.010+.002*AA+.001*AA**2)*SHW
C****
C HTPR IS THE HEAT OF PRECIPITATION.
C SHTPR IS THE SPECIFIC HEAT OF PRECIPITATION.
      SHTPR=0.
      IF(PR.GT.0.)SHTPR=HTPR/PR
C HTPRS IS THE HEAT OF LARGE SCALE PRECIPITATION
      HTPRS=SHTPR*PRS
C****
      RETURN
      END SUBROUTINE GHINIJ

      SUBROUTINE GHINHT (SNOWDP,TG1,TG2,WTR1,WTR2,ACE1,ACE2)
C**** INITIALIZES NEW GROUND (W,HT,SNW) FROM OLD (T,W,ICE,SNW)
C**** EVALUATES THE HEAT IN THE SOIL LAYERS BASED ON THE
C**** TEMPERATURES.
C**** INPUT:
C**** W - WATER IN SOIL LAYERS, M
C**** TP - TEMPERATURE OF LAYERS, C
C**** FICE - FRACTION OF ICE OF LAYERS
C**** FSN - HEAT OF FUSION OF WATER
C**** SHC - SPECIFIC HEAT CAPACITY OF SOIL
C**** SHI - SPECIFIC HEAT CAPACITY OF ICE
C**** SHW - SPECIFIC HEAT CAPCITY OF WATER
C**** SNOWD - SNOW DEPTH, EQUIVALENT WATER M
C**** OUTPUT:
C**** HT - HEAT IN SOIL LAYERS
C**** ADD CALCULATION OF WFC2
C**** BASED ON COMBINATION OF LAYERS 2-N, AS IN RETP2
      USE SLE001
      IMPLICIT NONE

      REAL*8 SNOWDP,TG1,TG2,WTR1,WTR2,ACE1,ACE2

      WFC1=FB*WS(1,1)+FV*(WS(0,2)+WS(1,2))
      WFC2=0.
      FBV=FB
      DO 30 IBV=1,2
      DO 20 L=2,N
      WFC2=WFC2+FBV*WS(L,IBV)
   20 CONTINUE
      FBV=FV
   30 CONTINUE
      WFC1=1000.*WFC1
      WFC2=1000.*WFC2
      FICE(0,2)=1.
      FICE(1,1)=(ACE1+SNOWDP*1000.)/(WTR1+ACE1+SNOWDP*1000.+1.D-20)
      FICE(1,2)=FICE(1,1)
      TP(0,2)=TG1
C**** W = SNOW(IF TOP LAYER) + WMIN + (WMAX-WMIN)*(WTR+ICE)/WFC
      W(0,2)=0.
      DO 1 IBV=1,2
      W(1,IBV)=SNOWDP
      WMIN=THETM(1,IBV)*DZ(1)
      WET1=(WTR1+ACE1)/(WFC1+1.D-20)
      IF(WET1.GT.1.) WET1=1.
      W(1,IBV)=W(1,IBV)+WMIN+(WS(1,IBV)-WMIN)*WET1
      SNOWD(IBV)=SNOWDP
      TP(1,IBV)=TG1
      DO 1 L=2,N
      FICE(L,IBV)=ACE2/(WTR2+ACE2+1.D-20)
      WMIN=THETM(L,IBV)*DZ(L)
      WET2=(WTR2+ACE2)/(WFC2+1.D-20)
      IF(WET2.GT.1.) WET2=1.
      W(L,IBV)=WMIN+(WS(L,IBV)-WMIN)*WET2
      TP(L,IBV)=TG2
    1 CONTINUE
C****
      ENTRY GHEXHT
C****
C**** COMPUTE HT (HEAT W/M+2)
      DO 10 IBV=1,2
      LL=2-IBV
      DO 10 L=LL,N
      IF(TP(L,IBV))2,4,6
    2 HT(L,IBV)=TP(L,IBV)*(SHC(L,IBV)+W(L,IBV)*SHI)-W(L,IBV)*FSN
      GO TO 10
    4 HT(L,IBV)=-FICE(L,IBV)*W(L,IBV)*FSN
      GO TO 10
    6 HT(L,IBV)=TP(L,IBV)*(SHC(L,IBV)+W(L,IBV)*SHW)
   10 CONTINUE
      IF(ID.EQ.0)THEN
       WRITE(99,*)'GHINHT ID CHECK',ID
       WRITE(99,*)'TG1,TG2',TG1,TG2
       WRITE(99,*)'TP',TP
       WRITE(99,*)'HT',HT
       WRITE(99,*)'W',W
       WRITE(99,*)'WTR1,WTR2',WTR1,WTR2
       WRITE(99,*)'ACE1,ACE2',ACE1,ACE2
       WRITE(99,*)'WFC1,WFC2',WFC1,WFC2
       WRITE(99,*)'SHC',SHC
       WRITE(99,*)'FICE',FICE
      ENDIF
      RETURN
      END SUBROUTINE GHINHT

      SUBROUTINE RETP2 (TG2AV,WTR2AV,ACE2AV)
C**** EVALUATES THE MEAN TEMPERATURE IN THE SOIL LAYERS 2-NGM
C**** AS WELL AS THE WATER AND ICE CONTENT.
C**** INPUT:
C**** W - WATER IN SOIL LAYERS, M
C**** HT - HEAT IN SOIL LAYERS
C**** FSN - HEAT OF FUSION OF WATER
C**** SHC - SPECIFIC HEAT CAPACITY OF SOIL
C**** SHI - SPECIFIC HEAT CAPACITY OF ICE
C**** SHW - SPECIFIC HEAT CAPCITY OF WATER
C**** OUTPUT:
C**** TG2AV - TEMPERATURE OF LAYERS 2 TO NGM, C
C**** ICE2AV - ICE AMOUNT IN LAYERS 2 TO NGM, KG/M+2
C**** WTR2AV - WATER IN LAYERS 2 TO NGM, KG/M+2
      USE SLE001
      IMPLICIT NONE
      REAL*8 TG2AV,WTR2AV,ACE2AV

      TG2AV=0.
      WTR2AV=0.
      ACE2AV=0.
      DO 3500 IBV=1,2
      WC=0.
      HTC=0.
      SHCC=0.
      DO 3420 L=2,N
      WC=WC+W(L,IBV)
      HTC=HTC+HT(L,IBV)
 3420 SHCC=SHCC+SHC(L,IBV)
      TPC=0.
      FICEC=0.
      IF(WC.NE.0.)  FICEC=-HTC/(FSN*WC)
      IF(FSN*WC+HTC.GE.0.)GO TO 3430
      TPC=(HTC+WC*FSN)/(SHCC+WC*SHI)
      FICEC=1.
      GO TO 3440
 3430 IF(HTC.LE.0.) GO TO 3440
      TPC=HTC/(SHCC+WC*SHW)
      FICEC=0.
 3440 CONTINUE
      FTP=FB
      IF(IBV.EQ.2) FTP=FV
      TG2AV=TG2AV+TPC*FTP
      WTR2AV=WTR2AV+WC*FTP*1000.*(1.-FICEC)
      ACE2AV=ACE2AV+WC*FTP*1000.*FICEC
 3500 CONTINUE
      RETURN
      END SUBROUTINE RETP2

      SUBROUTINE CHECKE(SUBR)
!@sum  CHECKE Checks whether gdata are reasonable over earth
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM
      USE GEOM
      USE GHYCOM
      IMPLICIT NONE

      REAL*8 X,TGL,WTRL,ACEL,PEARTH
      INTEGER I,J,K,IMAX
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in earth data
      CALL CHECK3(GHDATA,IM,JM,29,SUBR,'gh')

C**** Check for reasonable temperatures over earth
      X=1.001
      DO J=1,JM
         IMAX=IMAXJ(J)
         DO I=1,IMAX
            PEARTH = FEARTH(I,J)
            IF (PEARTH.GT.0.) THEN
               TGL=GDATA(I,J,4)
               WTRL=GDATA(I,J,5)
               ACEL=GDATA(I,J,6)
               IF ((TGL+60.)*(60.-TGL).LE.0.) WRITE (6,901) SUBR,I,J,TAU
     *              ,PEARTH,'TG1 ',(GDATA(I,J,K),K=2,6)
               IF (WTRL.LT.0..OR.ACEL.LT.0..OR.(WTRL+ACEL).GT.X*WFCS(I
     *              ,J)) WRITE(6,901) SUBR,I,J,TAU,PEARTH,'WTR '
     *              ,(GDATA(I,J,K),K=2,6),WFCS(I,J)
            END IF
         END DO
      END DO

      RETURN
 901  FORMAT ('0GDATA OFF, SUBR,I,J,TAU,PEARTH,',A7,2I4,F14.1,F5.2,1X
     *     ,A4/' SNW,x,TG1,WTR1,ICE1, WFC1 ',6F12.4)

      END SUBROUTINE CHECKE

      SUBROUTINE daily_SNOW
!@sum  daily_SNOW updates the snow ages every day
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : IM,JM,GDATA,NSURF
      USE GEOM, only : IMAXJ
      USE DAGCOM, only : AIJ,TDIURN,IJ_STRNGTS,IJ_DTGDTS,IJ_TMAXE
     *     ,IJ_TDSL,IJ_TMNMX
      IMPLICIT NONE
      REAL*8 TSAVG
      INTEGER I,J,IMAX

C**** INCREASE SNOW AGE EACH DAY (independent of Ts)
      DO J=1,JM
         IMAX=IMAXJ(J)
         DO I=1,IMAX
            GDATA(I,J,9)=1.+.98*GDATA(I,J,9)
            GDATA(I,J,10)=1.+.98*GDATA(I,J,10)
            GDATA(I,J,11)=1.+.98*GDATA(I,J,11)
            TSAVG=TDIURN(I,J,5)/(24.*NSURF)
            IF(32.+1.8*TSAVG.LT.65.)
     *           AIJ(I,J,IJ_STRNGTS)=AIJ(I,J,IJ_STRNGTS)+(33.-1.8*TSAVG)
            AIJ(I,J,IJ_DTGDTS)=AIJ(I,J,IJ_DTGDTS)+18.*((TDIURN(I,J,2)-
     *           TDIURN(I,J,1))/(TDIURN(I,J,4)-TDIURN(I,J,3)+1.D-20)-1.)
         AIJ(I,J,IJ_TDSL)=AIJ(I,J,IJ_TDSL)+(TDIURN(I,J,4)-TDIURN(I,J,3))
            AIJ(I,J,IJ_TMAXE)=AIJ(I,J,IJ_TMAXE)+(TDIURN(I,J,4)-273.16)
            IF (TDIURN(I,J,6).LT.AIJ(I,J,IJ_TMNMX))
     *           AIJ(I,J,IJ_TMNMX)=TDIURN(I,J,6)
         END DO
      END DO

      RETURN
      END SUBROUTINE daily_SNOW
