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
     *     ,sha,tf,rhow,rhoi,shw,shi,edpery,deltx
      USE MODEL_COM, only : im,jm,lm,t,p,q,DTsrc,NIsurf,dsig
     *     ,jday,JHOUR,NDAY,ITime,jeq,fearth,modrd,ijd6,itearth
     *     ,vt_on
      USE GEOM, only : imaxj,dxyp
      USE RADNCB, only : trhr,fsf,cosz1
      USE GHYCOM, only : wbare,wvege,htbare,htvege,snowbv,
     &     nsn_ij,isn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij,
     *     snowe,tearth,wearth,aiearth
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
     &    TGPASS=>TG,TKPASS=>T1,VGM=>VG,EDDY,
     &    NLSN,ISN,NSN,DZSN,WSN,HSN,FR_SNOW
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs,tsavg,qsavg
      USE SOCPBL, only : zgs
      USE DAGCOM , only : aij,tsfrez,tdiurn,aj,areg,adaily,jreg,
     *     ij_rune, ij_arunu, ij_pevap, ij_shdt, ij_beta, ij_trnfp0,
     *     ij_srtr, ij_neth, ij_ws, ij_ts, ij_us, ij_vs, ij_taus,
     *     ij_tauus, ij_tauvs, ij_qs, j_edifs, j_trhdt, j_shdt, j_evhdt,
     *     j_evap,j_erun1,j_difs,j_run2,j_dwtr2,j_run1,j_tsrf,j_f1dt,
     *     ij_g05,ij_g06,ij_g11,ij_g12,ij_g13,ij_g14,ij_g15,
     *     ij_g16,ij_g17,ij_g18,ij_g19,ij_g20,ij_g21,ij_g22,ij_g23,
     *     ij_g24,ij_g25,ij_g26,ij_g27
      USE DYNAMICS, only : pmid,pk,pek,pedn,pdsig
      USE FLUXES, only : dth1,dq1,du1,dv1,e0,e1,evapor,prec,eprec,runoe
     *     ,erunoe,gtemp

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NS,MODDSF,MODD6
      INTEGER I,J,L,KR,JR,ITYPE,IM1,IMAX,IHOUR
      REAL*8 SHDT,QSATS,EVAP,EVHDT,TG2AV,ACE2AV,TRHDT,CDN,RCDMWS
     *     ,RCDHWS,DHGS,CDQ,CDM,CDH,ELHX,TG,SRHEAT,TG1,PTYPE,QSATSS
     *     ,EVAPS,PPBLS,EDS1S,DBLS,DGSS,CDHS,CDMS,QGS,SRHDTS,SHDTS
     *     ,EVHDTS,TG1S,DXYPJ,TRHEAT,WTR2AV,WFC1,WGS,VGS,UGS,TRHDTS
     *     ,RTAUVS,RTAUUS,RTAUS,TAUS,WSS,TSS,QSS,USS,VSS,RHOSRF,RMBYA
     *     ,TFS,TH1,THV1,P1K,PSK,TS,PS,PIJ,PSOIL,PEARTH,WARMER,BRUN0
     *     ,BERUN0,BDIFS,BEDIFS,BTS,BEVHDT,BRUNU,BERUNU,BSHDT,BTRHDT
     *     ,TIMEZ,SPRING,ZS1CO,F1DTS
     *     ,rvx

c      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1,DU1,DV1
c      COMMON /WORK1d/DTH1,DQ1
c      COMMON /WORK2/DU1,DV1
c      REAL*8, DIMENSION(IM,JM,4) :: E0,E1,EVAPOR,TGRND
c      COMMON/WORK3/E0,E1,EVAPOR,TGRND
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

      REAL*8 QSAT
C****
C**** FEARTH    SOIL COVERED LAND FRACTION (1)
C****
C**** SNOWI     OCEAN ICE SNOW AMOUNT (KG/M**2)
C**** SNOWE     EARTH SNOW AMOUNT (KG/M**2)
C**** TSI(1:2)  OCEAN ICE TEMPERATURE OF FIRST/SECOND LAYER (C)
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)
C**** SNOWLI      LAND ICE SNOW AMOUNT (KG/M**2)
C**** TLANDI(1:2) LAND ICE TEMPERATURE OF FIRST/SECOND LAYER (C)
C****
C**** WBARE  1-6 WATER OF BARE SOIL LAYER 1-6 (M)
C**** WVEGE   0  WATER OF VEGETATION CANOPY (M)
C****        1-6 WATER OF VEGETATED SOIL LAYER 1-6 (M)
C**** HTBARE  0  BARE SOIL LAYER 0 IS UNUSED
C****        1-6 HEAT CONTENT OF BARE SOIL LAYER 1-6 (J M-2)
C**** HTVEGE  0  HEAT CONTENT OF VEGETATION CANOPY (J M-2)
C****        1-6 HEAT CONTENT OF VEGETATED SOIL LAYER 1-6 (J M-2)
C**** SNOWBV  1  SNOW DEPTH OVER BARE SOIL (M)
C****         2  SNOW DEPTH OVER VEGETATED SOIL (M)
C****
c     if(.not. vt_on) then
          rvx=0.
c     else
c         rvx=deltx !not working yet
c     endif

      DTSURF=DTsrc/NISURF
      ZS1CO=.5*DSIG(1)*RGAS/GRAV

         SPRING=-1.
         IF((JDAY.GE.32).AND.(JDAY.LE.212)) SPRING=1.
         IHOUR=1+JHOUR
C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NISURF TIMES EVERY HOUR
C****
         TIMEZ=JDAY+(MOD(ITime,NDAY)+(NS-1.)/NISURF)/NDAY  ! -1 ??
         IF(JDAY.LE.31) TIMEZ=TIMEZ+365.
         GHOUR=(ITime+(NS-1.)/NISURF) ! *(24./NDAY)
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=1,JM
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
         IF(J.LT.JEQ) WARMER=-SPRING
         IF(J.GE.JEQ) WARMER=SPRING
      IM1=IM
      DO I=1,IMAX
C****
C**** DETERMINE SURFACE CONDITIONS
C****
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
      RMBYA=100.*PDSIG(1,I,J)/GRAV
      QM1=Q1*RMBYA
c     RHOSRF=100.*PS/(RGAS*TKV)
      RHOSRF=100.*PS/(RGAS*TSV)
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
         F1DTS=0.
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
      PR=PREC(I,J)/(DTsrc*RHOW)
      PRS=PRCSS(I,J)/(DTsrc*RHOW)
      HTPR=EPREC(I,J)/DTsrc
      GW(1:NGM,1) =  WBARE(1:NGM,I,J)
      GW(0:NGM,2) =  WVEGE(0:NGM,I,J)
      HT(0:NGM,1) = HTBARE(0:NGM,I,J)
      HT(0:NGM,2) = HTVEGE(0:NGM,I,J)
      SNOWD(1:2)  = SNOWBV(1:2,I,J)  
ccc extracting snow variables
      NSN(1:2)          = NSN_IJ    (1:2, I, J)
      ISN(1:2)          = ISN_IJ    (1:2, I, J)
      DZSN(1:NLSN, 1:2) = DZSN_IJ   (1:NLSN, 1:2, I, J)
      WSN(1:NLSN, 1:2)  = WSN_IJ    (1:NLSN, 1:2, I, J)
      HSN(1:NLSN, 1:2)  = HSN_IJ    (1:NLSN, 1:2, I, J)
      FR_SNOW(1:2)      = FR_SNOW_IJ(1:2, I, J)
      CALL GHINIJ (I,J,WFC1)
      CALL RETH
      CALL RETP
C     CALL HYDRA
C??   SNOW = SNOWD(1)*FB + SNOWD(2)*FV
      TG1=TBCS
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
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
      QG=QSAT(TG,ELHX,PS)
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
C**** CALCULATE RHOSRF*CDM*WS
      RCDMWS=CDM*WS*RHOSRF
      RCDHWS=CDH*WS*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)
      SNHT=-SHA*RCDHWS*(TS-TG)
      TRHEAT=TRHR(1,I,J)
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

        WBARE(1:NGM,I,J) = GW(1:NGM,1)
        WVEGE(0:NGM,I,J) = GW(0:NGM,2)
       HTBARE(0:NGM,I,J) = HT(0:NGM,1)
       HTVEGE(0:NGM,I,J) = HT(0:NGM,2)
       SNOWBV(1:2,I,J)   = SNOWD(1:2)
ccc copy snow variables back to storage
      NSN_IJ    (1:2, I, J)         = NSN(1:2)
      ISN_IJ    (1:2, I, J)         = ISN(1:2)
      DZSN_IJ   (1:NLSN, 1:2, I, J) = DZSN(1:NLSN,1:2)
      WSN_IJ    (1:NLSN, 1:2, I, J) = WSN(1:NLSN,1:2)
      HSN_IJ    (1:NLSN, 1:2, I, J) = HSN(1:NLSN,1:2)
      FR_SNOW_IJ(1:2, I, J)         = FR_SNOW(1:2)

      AIJ(I,J,IJ_G05)=AIJ(I,J,IJ_G05)+BETAB/NIsurf
      AIJ(I,J,IJ_G06)=AIJ(I,J,IJ_G06)+BETAP/NIsurf
      AIJ(I,J,IJ_G11)=AIJ(I,J,IJ_G11)+BETA/NIsurf
      AIJ(I,J,IJ_G12)=AIJ(I,J,IJ_G12)+ACNA/NIsurf
      AIJ(I,J,IJ_G13)=AIJ(I,J,IJ_G13)+ACNC/NIsurf
      AIJ(I,J,IJ_G14)=AIJ(I,J,IJ_G14)+AEPP
      AIJ(I,J,IJ_G15)=AIJ(I,J,IJ_G15)+TP(1,1)
      AIJ(I,J,IJ_G16)=AIJ(I,J,IJ_G16)+TP(2,1)
      AIJ(I,J,IJ_G17)=AIJ(I,J,IJ_G17)+TP(3,1)
      AIJ(I,J,IJ_G18)=AIJ(I,J,IJ_G18)+EVAPB
      AIJ(I,J,IJ_G19)=AIJ(I,J,IJ_G19)+EVAPD
      AIJ(I,J,IJ_G20)=AIJ(I,J,IJ_G20)+EVAPW
      AIJ(I,J,IJ_G21)=AIJ(I,J,IJ_G21)+TP(0,2)
      AIJ(I,J,IJ_G22)=AIJ(I,J,IJ_G22)+TP(1,2)
      AIJ(I,J,IJ_G23)=AIJ(I,J,IJ_G23)+TP(2,2)
      AIJ(I,J,IJ_G24)=AIJ(I,J,IJ_G24)+TP(3,2)
      AIJ(I,J,IJ_G25)=AIJ(I,J,IJ_G25)+FB*ZW(1)+FV*ZW(2)
      AIJ(I,J,IJ_G26)=AIJ(I,J,IJ_G26)+BETAV/NIsurf
      AIJ(I,J,IJ_G27)=AIJ(I,J,IJ_G27)+BETAT/NIsurf
          TRHDT=TRHEAT*DTSURF-ATRG
C     FOR RADIATION FIND COMPOSITE VALUES OVER EARTH
C           FOR DIAGNOSTIC PURPOSES ALSO COMPUTE GDEEP 1 2 3
      SNOWE(I,J)=1000.*(SNOWD(1)*FB+SNOWD(2)*FV)
      TEARTH(I,J)=TG1
      WEARTH(I,J)=1000.*( FB*GW(1,1)*(1.-FICE(1,1)) +
     +  FV*(GW(1,2)*(1.-FICE(1,2))+GW(0,2)*(1.-FICE(0,2))) )
      AIEARTH(I,J)=1000.*( FB*GW(1,1)*FICE(1,1) +
     +  FV*(GW(1,2)*FICE(1,2)+GW(0,2)*FICE(0,2)) )
            CALL RETP2 (TG2AV,WTR2AV,ACE2AV)
            GDEEP(I,J,1)=TG2AV
            GDEEP(I,J,2)=WTR2AV
            GDEEP(I,J,3)=ACE2AV
            GTEMP(1,4,I,J)=TEARTH(I,J)
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
         QSATS=QSAT(TS,ELHX,PS)
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
         F1DTS=F1DTS+AF1DT*PTYPE
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
C**** Save runoff for addition to lake mass/energy resevoirs
         RUNOE (I,J)=RUNOE (I,J)+ ARUNS+ ARUNU
         ERUNOE(I,J)=ERUNOE(I,J)+AERUNS+AERUNU
C****
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
         AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+TRHDTS*DXYPJ
         AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+SHDTS*DXYPJ
         AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+EVHDTS*DXYPJ
         AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+EVAPS*DXYPJ
         AREG(JR,J_ERUN1)=AREG(JR,J_ERUN1)+AERUNS*PEARTH*DXYPJ
         AREG(JR,J_DIFS )=AREG(JR,J_DIFS )+DIFS*PEARTH*DXYPJ
         AREG(JR,J_RUN2 )=AREG(JR,J_RUN2 )+ARUNU*PEARTH*DXYPJ
         AREG(JR,J_DWTR2)=AREG(JR,J_DWTR2)+AERUNU*PEARTH*DXYPJ
         AREG(JR,J_RUN1 )=AREG(JR,J_RUN1 )+ARUNS*PEARTH*DXYPJ
         AREG(JR,J_F1DT )=AREG(JR,J_F1DT )+F1DTS*DXYPJ
         IF(MODDSF.NE.0) GO TO 5700
         AREG(JR,J_TSRF )=AREG(JR,J_TSRF )+(TSS-TFS)*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
 5700    AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDTS
         AIJ(I,J,IJ_BETA)=AIJ(I,J,IJ_BETA)+BETAD/NIsurf
         IF(MODRD.EQ.0) AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)+TRHDTS
     *        /DTSRC
         AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHDTS+TRHDTS)
         AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHDTS+TRHDTS+SHDTS+EVHDTS)
         IF(MODDSF.NE.0) GO TO 5800
         AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WSS          ! added 3/3/95 -rar-
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
 6000    IM1=I
      END DO
C**** QUANTITIES ACCUMULATED FOR SURFACE TYPE TABLES IN DIAGJ
         AJ(J,J_TRHDT,ITEARTH)=AJ(J,J_TRHDT,ITEARTH)+BTRHDT
         AJ(J,J_SHDT ,ITEARTH)=AJ(J,J_SHDT ,ITEARTH)+BSHDT
         AJ(J,J_EVHDT,ITEARTH)=AJ(J,J_EVHDT,ITEARTH)+BEVHDT
         AJ(J,J_ERUN1,ITEARTH)=AJ(J,J_ERUN1,ITEARTH)+BERUN0
         AJ(J,J_EDIFS,ITEARTH)=AJ(J,J_EDIFS,ITEARTH)+BEDIFS
         AJ(J,J_DIFS ,ITEARTH)=AJ(J,J_DIFS ,ITEARTH)+BDIFS
         AJ(J,J_RUN2 ,ITEARTH)=AJ(J,J_RUN2 ,ITEARTH)+BRUNU
         AJ(J,J_DWTR2,ITEARTH)=AJ(J,J_DWTR2,ITEARTH)+BERUNU
         AJ(J,J_RUN1 ,ITEARTH)=AJ(J,J_RUN1 ,ITEARTH)+BRUN0
         IF(MODDSF.EQ.0) AJ(J,J_TSRF,ITEARTH)=AJ(J,J_TSRF,ITEARTH)+BTS
      END DO
      RETURN
      END SUBROUTINE EARTH

      SUBROUTINE init_GH(DTSURF,redoGH,iniSNOW)
C**** Modifications needed for split of bare soils into 2 types
      USE CONSTANT, only : twopi,rhow,edpery,sha,shw_const=>shw,
     *     shi_const=>shi,lhe,lhm
      USE MODEL_COM, only : im,jm,fearth,vdata,Itime,Nday,jeq
      USE GHYCOM
      USE SLE001, sinday=>sint,cosday=>cost,fglob=>f
      USE FLUXES, only : gtemp
      USE DAGCOM, only : npts,icon_WTG,icon_HTG
      USE FILEMANAGER
      IMPLICIT NONE

      REAL*8 DTSURF
      INTEGER iu_SOIL
      INTEGER JDAY
      REAL*8 SNOWDP,WTR1,WTR2,ACE1,ACE2,TG1,TG2
      LOGICAL redoGH, iniSNOW
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER I, J
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
      call getunit("SOIL",iu_SOIL,.true.,.true.)
      CALL DREAD (iu_SOIL,DZ_IJ,IM*JM*(11*NGM+1),DZ_IJ)
      CLOSE (iu_SOIL)
C
      ONE=1.
      igcm = 0
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
      FSN= lhm * rhow
C ELH IS THE HEAT OF VAPORIZATION
      ELH= lhe * rhow
C THE SH'S ARE THE SPECIFIC HEAT CAPACATIES
      SHW= shw_const * rhow
      SHI= shi_const * rhow
c      SHA= sha_const
c      SHV=1911.
C THE ALAM'S ARE THE HEAT CONDUCTIVITIES
      ALAMW=.573345d0
      ALAMI=2.1762d0
      ALAMA=.025d0
      ALAMSN=0.088d0
      ALAMBR=2.9d0
      ALAMS(1)=8.8d0
      ALAMS(2)=2.9d0
      ALAMS(3)=2.9d0
      ALAMS(4)=.25d0
C HW IS THE WILTING POINT IN METERS
      HW=-100
C TFRZ IS 0 C IN K
c      TFRZ=273.16d0
C ZHTB IS DEPTH FOR COMBINING HEAT LAYERS FOR STABILITY
      IF(Q(4,1).LT..01)THEN
      ZHTB=6.
      ELSE
      ZHTB=6.
      ENDIF
C SPGSN IS THE SPECIFIG GRAVITY OF SNOW
      SPGSN=.1d0
C
C****
C**** Initialize global arrays  ALA, ACS, AFB, AFR
C****
c      TWOPI=6.283185    ! should be taken from CONSTANT

      ALA(:,:,:)=0.
      ACS(:,:,:)=0.
      AFB(:,:)=0.
      AFR(:,:,:)=0.
      ACS(1,:,:)=.01d0
      DO J=1,JM
        DO I=1,IM
          PEARTH=FEARTH(I,J)
          AFB(I,J)=VDATA(I,J,1)+VDATA(I,J,10)
          IF(AFB(I,J).GT..999) AFB(I,J)=1.
          IF(PEARTH.LE.0..OR.AFB(I,J).GE.1.) CYCLE
C**** CALCULATE LAI, CS COEFFICICENTS
          SFV=0.
          SLA0=0.
          SLRE=0.
          SLIM=0.
          SCS0=0.
          SCSRE=0.
          SCSIM=0.
          SVH=0.
          DO IV=1,8
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
            SCSIM=SCSIM+FV*DIF*SIN(PHASE)/RSAR(IV)
          END DO
          ALA(1,I,J)=.5/SFV*SLA0
          ALA(2,I,J)=.5/SFV*SLRE
          ALA(3,I,J)=.5/SFV*SLIM
          ACS(1,I,J)=.5/SFV*SCS0
          ACS(2,I,J)=.5/SFV*SCSRE
          ACS(3,I,J)=.5/SFV*SCSIM
          AVH(I,J)=SVH/SFV
C**** CALCULATE ROOT FRACTION AFR AVERAGED OVER VEGETATION TYPES
          DO N=1,NGM
            DZ(N)=DZ_IJ(I,J,N)
            IF(DZ(N).LE.0.) GO TO 320
          END DO
 320      N=N-1
          DO IV=1,8
            FV=VDATA(I,J,IV+1)
            Z=0.
            FRUP=0.
            DO L=1,N
              Z=Z+DZ(L)
              FRDN=AROOT(IV)*Z**BROOT(IV)
              FRDN=MIN(FRDN,ONE)
              IF(L.EQ.N)FRDN=1.
              AFR(L,I,J) = AFR(L,I,J) + FV*(FRDN-FRUP)
              FRUP=FRDN
            END DO
          END DO
          DO L=1,N
            AFR(L,I,J) = AFR(L,I,J)/(1.-AFB(I,J))
          END DO
        END DO
      END DO
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
C**** Recompute ground hydrology data if necessary (new soils data)
      IF (redoGH) THEN
        JDAY=1+MOD(ITime/NDAY,365)
        COSDAY=COS(TWOPI/EDPERY*JDAY)
        SINDAY=SIN(TWOPI/EDPERY*JDAY)

        DO J=1,JM
        DO I=1,IM
        PEARTH=FEARTH(I,J)
        IF(PEARTH.LE.0.) THEN

          WBARE(:,I,J)=0.
          WVEGE(:,I,J)=0.
          HTBARE(:,I,J)=0.
          HTVEGE(:,I,J)=0.
          SNOWBV(:,I,J)=0.

        ELSE
ccc??? remove next 5 lines? -check the old version
           W(1:NGM,1) =   WBARE(1:NGM,I,J)
           W(0:NGM,2) =   WVEGE(0:NGM,I,J)
           HT(0:NGM,1) = HTBARE(0:NGM,I,J)
           HT(0:NGM,2) = HTVEGE(0:NGM,I,J)
           SNOWD(1:2) =  SNOWBV(1:2,I,J)  

C****     COMPUTE SOIL HEAT CAPACITY AND GROUND WATER SATURATION GWS
          CALL GHINIJ (I,J,WFC1)
C****     FILL IN SOILS COMMON BLOCKS
          SNOWDP=SNOWE(I,J)/RHOW
          WTR1=WEARTH(I,J)
          ACE1=AIEARTH(I,J)
          TG1 =TEARTH(I,J)
          WTR2=WTR1
          ACE2=ACE1
          TG2 =TG1
c          WTR2=GDATA(I,J,9)   ! this cannot be right
c          ACE2=GDATA(I,J,10)
c          TG2 =GDATA(I,J,8)
          CALL GHINHT (SNOWDP, TG1,TG2, WTR1,WTR2, ACE1,ACE2)

C****     COPY SOILS PROGNOSTIC QUANTITIES TO MODEL VARIABLES
            WBARE(1:NGM,I,J) = W(1:NGM,1)
            WVEGE(0:NGM,I,J) = W(0:NGM,2)
           HTBARE(0:NGM,I,J) = HT(0:NGM,1)
           HTVEGE(0:NGM,I,J) = HT(0:NGM,2)
           SNOWBV(1:2,I,J)   = SNOWD(1:2)
        END IF
      END DO
      END DO
      WRITE (*,*) 'GROUND HYDROLOGY DATA WAS MADE FROM GROUND DATA'
      END IF
C**** set gtemp array
      DO J=1,JM
        DO I=1,IM
          IF (FEARTH(I,J).gt.0) THEN
            GTEMP(1,4,I,J)=TEARTH(I,J)
          END IF
        END DO
      END DO

ccc   some extra code from snowmodel GHINIT
      SO_%ROSMP = 8.  ! no idea what this number means, but it is used
                      ! in computation of RUNOFF

ccc   init snow here
ccc hope this is the right place to split first layer into soil
ccc and snow  and to set snow arrays
ccc!!! this should be done only when restarting from an old
ccc!!! restart file (without snow model data)

      IF (iniSNOW) THEN
      DO I=1,IM
        DO J=1,JM
          PEARTH=FEARTH(I,J)
          IF(PEARTH.LE.0.) THEN
            NSN_IJ(:,I,J)     = 0
            ISN_IJ(:,I,J)     = 0
            DZSN_IJ(:,:,I,J)  = 0.
            WSN_IJ(:,:,I,J)   = 0.
            HSN_IJ(:,:,I,J)   = 0.
            FR_SNOW_IJ(:,I,J) = 0.
          ELSE
            JDAY=1+MOD(ITime/NDAY,365)
            COSDAY=COS(TWOPI/EDPERY*JDAY)
            SINDAY=SIN(TWOPI/EDPERY*JDAY)

            W(1:NGM,1) =   WBARE(1:NGM,I,J)
            W(0:NGM,2) =   WVEGE(0:NGM,I,J)
            HT(0:NGM,1) = HTBARE(0:NGM,I,J)
            HT(0:NGM,2) = HTVEGE(0:NGM,I,J)
            SNOWD(1:2) =  SNOWBV(1:2,I,J)  

            CALL GHINIJ (I,J,WFC1)
            CALL SET_SNOW

            NSN_IJ    (1:2, I, J)         = NSN(1:2)
            ISN_IJ    (1:2, I, J)         = ISN(1:2)
            DZSN_IJ   (1:NLSN, 1:2, I, J) = DZSN(1:NLSN,1:2)
            WSN_IJ    (1:NLSN, 1:2, I, J) = WSN(1:NLSN,1:2)
            HSN_IJ    (1:NLSN, 1:2, I, J) = HSN(1:NLSN,1:2)
            FR_SNOW_IJ(1:2, I, J)         = FR_SNOW(1:2)

C****     COPY SOILS PROGNOSTIC QUANTITIES TO MODEL VARIABLES
              WBARE(1:NGM,I,J) = W(1:NGM,1)
              WVEGE(0:NGM,I,J) = W(0:NGM,2)
             HTBARE(0:NGM,I,J) = HT(0:NGM,1)
             HTVEGE(0:NGM,I,J) = HT(0:NGM,2)
             SNOWBV(1:2,I,J)   = SNOWD(1:2)

          END IF
        END DO
      END DO
      END IF

C**** Set conservation diagnostics for ground water mass and energy
      QCON=(/ F, F, F, F, F, T, F, F, T, F, F/)
      CALL SET_CON(QCON,"GRND WTR","(KG/M^2)        ",
     *     "(10^-9 KG/S/M^2)",1d0,1d9,icon_WTG)
      QCON=(/ F, F, F, F, F, T, F, F, T, F, F/)
      CALL SET_CON(QCON,"GRND ENG","(10**6 J/M^2)   ",
     *     "(10^-3 J/S/M^2) ",1d-6,1d3,icon_HTG)

      RETURN
      END SUBROUTINE init_GH

      SUBROUTINE GHINIJ (I0,J0, WFCAP)
C**** INPUT:
C**** AVH(I,J) - ARRAY OF VEGETATION HEIGHTS
C**** SPGSN - SPECIFIC GRAVITY OF SNOW
C**** OUTPUT:
C**** VH - VEGETATION HEIGHT
C**** SNOWM - SNOW MASKING DEPTH
C**** WFCAP - WATER FIELD CAPACITY OF TOP SOIL LAYER, M
C****
      USE GHYCOM, only : dz_ij,sl_ij,q_ij,qk_ij,avh,afr,afb,ala,acs
      USE SLE001, only : dz,qk,ngm,imt,ng,zb,zc,fr,spgsn,q,sl,xklh0
     *     ,fb,fv,snowm,vh,alai,alaic,alaie,cost,sint,rs,prs,ijdebug,n
     *     ,thets,thetm,ws,thm,nth,shc,shcap,shw,shtpr,htprs,pr
     *     ,htpr
      USE snow_model, only : i_earth,j_earth
      IMPLICIT NONE
      INTEGER I0,J0
      REAL*8 WFCAP
      INTEGER L,IBV,K,I
      REAL*8 AA,ONE

      ONE=1.
      IJdebug=I0*100+J0
      i_earth = I0
      j_earth = J0
C**** SET UP LAYERS
      DZ(1:NGM)=DZ_IJ(I0,J0,1:NGM)
      Q(1:IMT,1:NGM)=Q_IJ(I0,J0,1:IMT,1:NGM)
      QK(1:IMT,1:NGM)=QK_IJ(I0,J0,1:IMT,1:NGM)
      SL=SL_IJ(I0,J0)
      DO N=1,NGM
        IF(DZ(N).LE.0.) GO TO 21
      END DO
   21 N=N-1
      IF(N.LE.0) THEN
         WRITE (99,*) 'GHINIJ:  N <= 0:  I,J,N=',I0,J0,N,(DZ(K),K=1,43)
         STOP
      END IF
C**** CALCULATE THE BOUNDARIES, BASED ON THE THICKNESSES.
      ZB(1)=0.
      DO L=1,N
        ZB(L+1)=ZB(L)-DZ(L)
      END DO
C**** CALCULATE THE LAYER CENTERS, BASED ON THE BOUNDARIES.
      DO L=1,N
        ZC(L)=.5*(ZB(L)+ZB(L+1))
      END DO
C**** FR: ROOT FRACTION IN LAYER L  (1=FR(1)+FR(2)+...+FR(N))
      DO L=1,N
        FR(L)=AFR(L,I0,J0)
      END DO
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
      DO IBV=1,2
        DO L=1,N
          THETS(L,IBV)=0.
          THETM(L,IBV)=0.
          DO I=1,IMT-1
            THETS(L,IBV)=THETS(L,IBV)+Q(I,L)*THM(0,I)
            THETM(L,IBV)=THETM(L,IBV)+Q(I,L)*THM(NTH,I)
          END DO
          WS(L,IBV)=THETS(L,IBV)*DZ(L)
        END DO
      END DO
      WS(0,2)=.0001*ALAI
      WFCAP=FB*WS(1,1)+FV*(WS(0,2)+WS(1,2))
C****
      CALL XKLH0
C****
      DO IBV=1,2
        DO L=1,N
          SHC(L,IBV)=0.
          DO I=1,IMT
            SHC(L,IBV)=SHC(L,IBV)+Q(I,L)*SHCAP(I)
          END DO
          SHC(L,IBV)=(1.-THETS(L,IBV))*SHC(L,IBV)*DZ(L)
        END DO
      END DO
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
      DO IBV=1,2
        W(1,IBV)=SNOWDP
        WMIN=THETM(1,IBV)*DZ(1)
        WET1=(WTR1+ACE1)/(WFC1+1.D-20)
        IF(WET1.GT.1.) WET1=1.
        W(1,IBV)=W(1,IBV)+WMIN+(WS(1,IBV)-WMIN)*WET1
        SNOWD(IBV)=SNOWDP
        TP(1,IBV)=TG1
        DO L=2,N
          FICE(L,IBV)=ACE2/(WTR2+ACE2+1.D-20)
          WMIN=THETM(L,IBV)*DZ(L)
          WET2=(WTR2+ACE2)/(WFC2+1.D-20)
          IF(WET2.GT.1.) WET2=1.
          W(L,IBV)=WMIN+(WS(L,IBV)-WMIN)*WET2
          TP(L,IBV)=TG2
        END DO
      END DO
C****
      ENTRY GHEXHT
C****
C**** COMPUTE HT (HEAT W/M+2)
      DO IBV=1,2
        LL=2-IBV
        DO L=LL,N
          IF(TP(L,IBV)) 2,4,6
 2        HT(L,IBV)=TP(L,IBV)*(SHC(L,IBV)+W(L,IBV)*SHI)-W(L,IBV)*FSN
          CYCLE
 4        HT(L,IBV)=-FICE(L,IBV)*W(L,IBV)*FSN
          CYCLE
 6        HT(L,IBV)=TP(L,IBV)*(SHC(L,IBV)+W(L,IBV)*SHW)
        END DO
      END DO
      IF(IJdebug.EQ.0)THEN
       WRITE(99,*)'GHINHT ID CHECK',IJdebug
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
      DO L=2,N
        WC=WC+W(L,IBV)
        HTC=HTC+HT(L,IBV)
        SHCC=SHCC+SHC(L,IBV)
      END DO
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
!@sum  CHECKE Checks whether arrays are reasonable over earth
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,fearth,ITime,wfcs
      USE GEOM, only : imaxj
      USE GHYCOM, only : tearth,wearth,aiearth,snowe,wbare,wvege,htbare
     *     ,htvege,snowbv,ngm
      IMPLICIT NONE

      REAL*8 X,TGL,WTRL,ACEL
      INTEGER I,J,K,IMAX
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in earth data
      CALL CHECK3(WBARE ,NGM  ,IM,JM,SUBR,'wb')
      CALL CHECK3(WVEGE ,NGM+1,IM,JM,SUBR,'wv')
      CALL CHECK3(HTBARE,NGM+1,IM,JM,SUBR,'hb')
      CALL CHECK3(HTVEGE,NGM+1,IM,JM,SUBR,'hv')
      CALL CHECK3(SNOWBV,2    ,IM,JM,SUBR,'sn')

C**** Check for reasonable temperatures over earth
      X=1.001
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF (FEARTH(I,J).GT.0.) THEN
            TGL=TEARTH(I,J)
            WTRL=WEARTH(I,J)
            ACEL=AIEARTH(I,J)
            IF ((TGL+60.)*(60.-TGL).LE.0.) WRITE (6,901) SUBR,I,J,ITime
     *           ,FEARTH(I,J),'TG1 ',SNOWE(I,J),TGL,WTRL,ACEL
            IF (WTRL.LT.0..OR.ACEL.LT.0..OR.(WTRL+ACEL).GT.X*WFCS(I
     *           ,J)) WRITE(6,901) SUBR,I,J,ITime,FEARTH(I,J),'WTR '
     *           ,SNOWE(I,J),TGL,WTRL,ACEL,WFCS(I,J)
          END IF
        END DO
      END DO

      RETURN
 901  FORMAT ('0GDATA OFF, SUBR,I,J,I-Time,PEARTH,',A7,2I4,I10,F5.2,1X
     *     ,A4/' SNW,x,TG1,WTR1,ICE1, WFC1 ',6F12.4)

      END SUBROUTINE CHECKE

      SUBROUTINE daily_EARTH(IEND)
!@sum  daily_EARTH performs daily tasks for EARTH related functions
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : rhow,twopi,edpery,tf
      USE MODEL_COM, only : im,jm,NDAY,NIsurf,jday,fearth,wfcs
      USE GEOM, only : imaxj
      USE DAGCOM, only : aij,tdiurn,ij_strngts,ij_dtgdts,ij_tmaxe
     *     ,ij_tdsl,ij_tmnmx
      USE SLE001
     &  , only : cosday=>cost, sinday=>sint
      USE GHYCOM, only : snoage
      IMPLICIT NONE
      REAL*8 TSAVG,WFC1
      INTEGER I,J,IMAX,ITYPE
      INTEGER, INTENT(IN) :: IEND  !@var IEND 1 if at end of day
C****
C**** FIND LEAF-AREA INDEX & WATER FIELD CAPACITY FOR GROUND LAYER 1
C****
      COSDAY=COS(TWOPI/EDPERY*JDAY)
      SINDAY=SIN(TWOPI/EDPERY*JDAY)
      DO J=1,JM
        DO I=1,IM
          WFCS(I,J)=24.
          IF (FEARTH(I,J).GT.0.) THEN
            CALL GHINIJ(I,J,WFC1)
            WFCS(I,J)=RHOW*WFC1 ! canopy part changes
          END IF
        END DO
      END DO

      IF (IEND.eq.1) THEN
C****
C**** INCREASE SNOW AGE EACH DAY (independent of Ts)
C****
        DO J=1,JM
          IMAX=IMAXJ(J)
          DO I=1,IMAX
            DO ITYPE=1,3
              SNOAGE(ITYPE,I,J)=1.+.98*SNOAGE(ITYPE,I,J)
            END DO
            TSAVG=TDIURN(I,J,5)/(NDAY*NIsurf)
            IF(32.+1.8*TSAVG.LT.65.)
     *           AIJ(I,J,IJ_STRNGTS)=AIJ(I,J,IJ_STRNGTS)+(33.-1.8*TSAVG)
            AIJ(I,J,IJ_DTGDTS)=AIJ(I,J,IJ_DTGDTS)+18.*((TDIURN(I,J,2)-
     *           TDIURN(I,J,1))/(TDIURN(I,J,4)-TDIURN(I,J,3)+1.D-20)-1.)
            AIJ(I,J,IJ_TDSL)=AIJ(I,J,IJ_TDSL)+
     *           (TDIURN(I,J,4)-TDIURN(I,J,3))
            AIJ(I,J,IJ_TMAXE)=AIJ(I,J,IJ_TMAXE)+(TDIURN(I,J,4)-TF)
            IF (TDIURN(I,J,6).LT.AIJ(I,J,IJ_TMNMX))
     *           AIJ(I,J,IJ_TMNMX)=TDIURN(I,J,6)
          END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE daily_EARTH

      SUBROUTINE GROUND_E
!@sum  GROUND_E driver for applying surface fluxes to land fraction
!@auth Original Development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,fearth,itearth
      USE GEOM, only : imaxj,dxyp
      USE GHYCOM, only : snowe, tearth,wearth,aiearth,wbare,wvege,snowbv
      USE DAGCOM, only : aj,areg,aij,jreg,ij_evap,ij_f0e,ij_evape
     *     ,ij_gwtr,ij_tg1,j_tg2,j_tg1,j_wtr1,j_ace1,j_wtr2,j_ace2
     *     ,j_snow,j_f2dt,j_f1dt,j_evap,j_type,ij_g01,ij_g07,ij_g28
     *     ,ij_g29,j_rsnow,ij_rsnw,ij_rsit,ij_snow
      USE FLUXES, only : e0,e1,evapor,eprec
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,DXYPJ,WTR1,WTR2,ACE1,ACE2
     *     ,PEARTH,ENRGP,SCOVE
      INTEGER I,J,IMAX,JR,K
      REAL*8, DIMENSION(IM,JM,3) :: GDEEP
      COMMON/oldDAG/GDEEP

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      PEARTH=FEARTH(I,J)
      JR=JREG(I,J)
      IF (PEARTH.gt.0) THEN

        SNOW=SNOWE(I,J)
        TG1 = TEARTH(I,J)
        WTR1= WEARTH(I,J)
        ACE1=AIEARTH(I,J)
        TG2=GDEEP(I,J,1)
        WTR2=GDEEP(I,J,2)
        ACE2=GDEEP(I,J,3)
        F0DT=E0(I,J,4)
        F1DT=E1(I,J,4)
        EVAP=EVAPOR(I,J,4)
        ENRGP=EPREC(I,J)      ! including latent heat

C**** ACCUMULATE DIAGNOSTICS
        SCOVE=0.
        IF (SNOWE(I,J).GT.0.) SCOVE=PEARTH
        AJ(J,J_RSNOW,ITEARTH)=AJ(J,J_RSNOW,ITEARTH)+SCOVE
        AREG(JR,J_RSNOW)=AREG(JR,J_RSNOW)+SCOVE*DXYPJ
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVE
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*PEARTH
        AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+SCOVE

        AJ(J,J_WTR1,ITEARTH)=AJ(J,J_WTR1,ITEARTH)+WTR1*PEARTH
        AJ(J,J_ACE1,ITEARTH)=AJ(J,J_ACE1,ITEARTH)+ACE1*PEARTH
        AJ(J,J_WTR2,ITEARTH)=AJ(J,J_WTR2,ITEARTH)+WTR2*PEARTH
        AJ(J,J_ACE2,ITEARTH)=AJ(J,J_ACE2,ITEARTH)+ACE2*PEARTH
        AJ(J,J_TG1 ,ITEARTH)=AJ(J,J_TG1, ITEARTH)+TG1 *PEARTH
        AJ(J,J_TG2 ,ITEARTH)=AJ(J,J_TG2, ITEARTH)+TG2 *PEARTH
        AJ(J,J_TYPE,ITEARTH)=AJ(J,J_TYPE,ITEARTH)+     PEARTH
        AJ(J,J_SNOW,ITEARTH)=AJ(J,J_SNOW,ITEARTH)+SNOW*PEARTH
        AJ(J,J_F1DT,ITEARTH)=AJ(J,J_F1DT,ITEARTH)+F1DT*PEARTH
        AJ(J,J_EVAP,ITEARTH)=AJ(J,J_EVAP,ITEARTH)+EVAP*PEARTH
        IF (JR.ne.24) THEN
        AREG(JR,J_TG1) =AREG(JR,J_TG1) +TG1 *PEARTH*DXYPJ
        AREG(JR,J_TG2) =AREG(JR,J_TG2) +TG2 *PEARTH*DXYPJ
        AREG(JR,J_SNOW)=AREG(JR,J_SNOW)+SNOW*PEARTH*DXYPJ
        AREG(JR,J_WTR1)=AREG(JR,J_WTR1)+WTR1*PEARTH*DXYPJ
        AREG(JR,J_ACE1)=AREG(JR,J_ACE1)+ACE1*PEARTH*DXYPJ
        AREG(JR,J_WTR2)=AREG(JR,J_WTR2)+WTR2*PEARTH*DXYPJ
        AREG(JR,J_ACE2)=AREG(JR,J_ACE2)+ACE2*PEARTH*DXYPJ
        END IF
        AIJ(I,J,IJ_F0E)  =AIJ(I,J,IJ_F0E)  +F0DT+ENRGP
        AIJ(I,J,IJ_TG1)  =AIJ(I,J,IJ_TG1)  +TG1 *PEARTH
        AIJ(I,J,IJ_GWTR) =AIJ(I,J,IJ_GWTR)+(WTR1+ACE1+WTR2+ACE2)
        AIJ(I,J,IJ_EVAP) =AIJ(I,J,IJ_EVAP) +EVAP*PEARTH
        AIJ(I,J,IJ_EVAPE)=AIJ(I,J,IJ_EVAPE)+EVAP
        DO K=1,4
          AIJ(I,J,IJ_G01+K-1)=AIJ(I,J,IJ_G01+K-1)+WBARE(K,I,J)
          AIJ(I,J,IJ_G07+K-1)=AIJ(I,J,IJ_G07+K-1)+WVEGE(K-1,I,J)
        END DO
        AIJ(I,J,IJ_G28)=AIJ(I,J,IJ_G28)+SNOWBV(1,I,J)
        AIJ(I,J,IJ_G29)=AIJ(I,J,IJ_G29)+SNOWBV(2,I,J)
      END IF
C****
      END DO
      END DO
      END SUBROUTINE GROUND_E

      SUBROUTINE conserv_WTG(WATERG)
!@sum  conserv_WTG calculates zonal ground water
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow
      USE MODEL_COM, only : im,jm,fim,fearth
      USE GEOM, only : imaxj
      USE GHYCOM, only : wbare,wvege,afb
      USE SLE001, only : ngm
      IMPLICIT NONE
!@var WATERG zonal ground water (kg/m^2)
      REAL*8, DIMENSION(JM) :: WATERG
      INTEGER I,J,N
      REAL*8 WIJ,FB

      DO J=1,JM
        WATERG(J)=0
        DO I=1,IMAXJ(J)
          IF (FEARTH(I,J).gt.0) THEN
            FB=AFB(I,J)
            WIJ=(1.-FB)*WVEGE(0,I,J)
            DO N=1,NGM
              WIJ=WIJ+FB*WBARE(N,I,J)+(1.-FB)*WVEGE(N,I,J)
            END DO
            WATERG(J)=WATERG(J)+FEARTH(I,J)*WIJ*RHOW
          END IF
        END DO
      END DO
      WATERG(1) =FIM*WATERG(1)
      WATERG(JM)=FIM*WATERG(JM)
C****
      END SUBROUTINE conserv_WTG

      SUBROUTINE conserv_HTG(HEATG)
!@sum  conserv_HTG calculates zonal ground energy
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,fearth
      USE GEOM, only : imaxj
      USE GHYCOM, only : htbare,htvege,afb
      USE SLE001, only : ngm
      IMPLICIT NONE
!@var HEATG zonal ground heat (J/m^2)
      REAL*8, DIMENSION(JM) :: HEATG
      INTEGER I,J,N
      REAL*8 HIJ,FB

      DO J=1,JM
        HEATG(J)=0
        DO I=1,IMAXJ(J)
          IF (FEARTH(I,J).gt.0) THEN
            FB=AFB(I,J)
            HIJ=0.
            DO N=0,NGM
              HIJ=HIJ+FB*HTBARE(N,I,J)+(1.-FB)*HTVEGE(N,I,J)
            END DO
            HEATG(J)=HEATG(J)+FEARTH(I,J)*HIJ
          END IF
        END DO
      END DO
      HEATG(1) =FIM*HEATG(1)
      HEATG(JM)=FIM*HEATG(JM)
C****
      END SUBROUTINE conserv_HTG
