C**** EE001M12 E001M12 SOMTQ EB357M12
C****
C**** Subroutine EARTH used by new land surface model.
C**** Coded to use new SOC PBL routines.
      SUBROUTINE EARTH (NS,MODDSF,MODD6)
C****
C**** This subroutine calculates surface fluxes of sensible heat,
C**** evaporation, thermal radiation, and momentum drag.
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'E001M12.COM'
      COMMON U,V,T,P,Q
      COMMON/WORK1/CONV(IM,JM,LM),PK(IM,JM,LM),PREC(IM,JM),
     2             TPREC(IM,JM),TAUSS(IM,JM,LM),TAUMC(IM,JM,LM),
     4             COSZ1(IM,JM),DTH1(IM,JM),DQ1(IM,JM)
      COMMON/WORK2/UT(IM,JM,LM),VT(IM,JM,LM),DU1(IM,JM),
     *  DV1(IM,JM),RA(8),ID(8),UMS(8)
      COMMON/WORK3/E0(IM,JM,4),E1(IM,JM,4),EVAPOR(IM,JM,4),
     *  TGRND(IM,JM,4),BLTEMP(IM,JM,8)
         COMMON/oldDAG/GDEEP(IM,JM,3)
      COMMON/RDATA/ROUGHL(IM,JM),FSF(IM,JM,4)
C**** Interface to SOILS subroutines
      PARAMETER (NGM=6)
      COMMON/SOILS2/SDATA(IM,JM,9*NGM+1)
      COMMON/SOILS3/GHDATA(IM,JM,4*NGM+5)
      COMMON/SOILS/PR,HTPR,PRS,HTPRS,GW(0:NGM,2),HT(0:NGM,2)
     & ,SNOWD(2),GWS(0:NGM,2),TP(0:NGM,2),FICE(0:NGM,2),GHOUR,COSDAY
     * ,SINDAY,GDZ(NGM),QSOIL(5,NGM),QKSOIL(5,NGM),SL,FV,FB,ALAI,ATRG
     & ,ASHG,ALHG,BETAD,BETAV,BETAT,BETAP,BETAB,BETA,ACNA,ACNC,EVAPW
     & ,EVAPD,EVAPB,ARUNS,ARUNU,AERUNS,AERUNU,DIFS,EDIFS,AEPC,AEPB
     & ,AEPP,AFHG,AF0DT,AF1DT,CNC,ZW(2),FD,FW,FM,VH,ALAIE,TBCS
     & ,TCS,SNOWM,IDPASS,IGCM
      COMMON/EVPTR/QM1,Q1,QS
      COMMON /WORKLS/PRCSS(IM,JM)
C**** Interface to SURFCE and land surface subroutines:

      COMMON/SURFC/PRES,RHO,TSPASS,VSM,CH,CD,SNHT,SRHT,TRHT,ZS,
     *             ZMIXE,CN,P1,PBLP,FCR,FMAG,TGPASS,TKPASS,VGM,SINA,
     *             EDDY,xxxRI1,xxxRIS,FPAD,FS(2),FS1(2,2),T1G,TSG,IZ,IZX

      COMMON /PBLPAR/ZGS,ZS1,PIJ,PSK,TGV,TKV,THV1,QG,HEMI,SHA,
     2               OMEGA2,DTSURF,JVPO,IQ1,IQ2,IQ3,IM1,POLE
      COMMON /PBLOUT/US,VS,WS,TSV,QSRF,PSI,DBL,EDVISC,EDS1,
     2               USTAR,PPBL,CDM,CDH,CDQ,UG,VG,WG,ZMIX

      parameter (npbl=8)
      common /socabl/uabl(npbl,im,jm,4),vabl(npbl,im,jm,4),
     2               tabl(npbl,im,jm,4),qabl(npbl,im,jm,4),
     3               eabl(npbl,im,jm,4),
     4               cmgs(im,jm,4),chgs(im,jm,4),cqgs(im,jm,4),
     5               ipbl(im,jm,4)

      LOGICAL POLE
      DATA SHV/0./,SHW/4185./,SHI/2060./,RHOW/1000./,RHOI/916.6/,
     *  STBO/.5672573E-7/,TF/273.16/,TFO/-1.80/, EDPERY/365./
      QSAT(TM,PM,QLH)=3.797915*EXP(QLH*(7.93252D-6-2.166847D-3/TM))/PM
      DATA IFIRST/1/
C****
C**** FDATA  2  LAND COVERAGE (1)
C****        3  RATIO OF LAND ICE COVERAGE TO LAND COVERAGE (1)
C****
C**** ODATA  1  OCEAN TEMPERATURE (C)
C****        2  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****        3  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
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
C**** BLDATA 1  COMPOSITE SURFACE WIND MAGNITUDE (M/S)
C****        2  COMPOSITE SURFACE AIR TEMPERATURE (K)
C****        3  COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
C****        4  LAYER TO WHICH DRY CONVECTION MIXES (1)
C****        5  MIXED LAYER DEPTH (Z1O NOT YET PART OF RESTART FILE)
C****        6  COMPOSITE SURFACE U WIND
C****        7  COMPOSITE SURFACE V WIND
C****        8  COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
C****
C**** ROUGHL    LOG(30./ROUGHNESS LENGTH) (LOGARITHM TO BASE 10)
C****
      IF(IFIRST.EQ.0) GO TO 30
      IFIRST=0
      IQ1=IM/4+1
      IQ2=IM/2+1
      IQ3=3*IM/4+1
      DTSURF=NDYN*DT/NSURF
      DTCNDS=NDYN*DT
         DTSRCE=DT*NDYN
      SHA=RGAS/KAPA
      RVX=0.
      ZS1CO=.5*DSIG(1)*RGAS/GRAV
   30 S0=S0X*1367./RSDIST
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
      POLE=.FALSE.
      IF(J.EQ.1) THEN
C**** CONDITIONS AT THE SOUTH POLE
        POLE=.TRUE.
        IMAX=1
        JVPO=2
        RAPO=2.*RAPVN(1)
        GO TO 100
      ENDIF
      IF(J.EQ.JM) THEN
C**** CONDITIONS AT THE NORTH POLE
        POLE=.TRUE.
        IMAX=1
        JVPO=JM
        RAPO=2.*RAPVS(JM)
        GO TO 100
      ENDIF
      POLE=.FALSE.
      IMAX=IM
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
         JEQ=1+JM/2
         IF(J.LT.JEQ) WARMER=-SPRING
         IF(J.GE.JEQ) WARMER=SPRING
      IM1=IM
      DO 6000 I=1,IMAX
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FDATA(I,J,2)
      PLICE=FDATA(I,J,3)*PLAND
      PEARTH=PLAND-PLICE
      PSOIL=PEARTH
      PIJ=P(I,J)
      PS=PIJ+PTOP
      PSK=EXPBYK(PS)
      P1=SIG(1)*PIJ+PTOP
      P1K=PK(I,J,1)
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
      ZGS=10.
      ITYPE=4
      PTYPE=PEARTH
      PR=PREC(I,J)/(DTCNDS*RHOW)
      PRS=PRCSS(I,J)/(DTCNDS*RHOW)
      HTPR=0.
      IF(TPREC(I,J).LT.0.) HTPR=-LHM*PREC(I,J)/DTCNDS
      DO 2410 L=1,4*NGM+5
 2410 GW(L,1)=GHDATA(I,J,L)
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
      P1=SIG(1)*PIJ+PTOP
C**** LOOP OVER GROUND TIME STEPS
      TG=TG1+TF
      ELHX=LHE
      IF(TG1.LT.0.)  ELHX=LHS
      QG=QSAT(TG,PS,ELHX)
      TGV=TG*(1.+QG*RVX)
C***********************************************************************
C***
      CALL PBL(I,J,ITYPE)
C***
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
      DO 3410 L=1,4*NGM+5
 3410 GHDATA(I,J,L)=GW(L,1)
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
         AIJ(I,J,32)=AIJ(I,J,32)+ARUNS
         AIJ(I,J,53)=AIJ(I,J,53)+ARUNU
         AIJ(I,J,79)=AIJ(I,J,79)+(AEPC+AEPB)
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
 5000 BLDATA(I,J,1)=WSS+BLTEMP(I,J,1)
      BLDATA(I,J,2)=TSS+BLTEMP(I,J,2)
      BLDATA(I,J,3)=QSS+BLTEMP(I,J,3)
      BLDATA(I,J,6)=USS+BLTEMP(I,J,6)
      BLDATA(I,J,7)=VSS+BLTEMP(I,J,7)
      BLDATA(I,J,8)=TAUS+BLTEMP(I,J,8)
         TDIURN(I,J,5)=TDIURN(I,J,5)+(BLDATA(I,J,2)-TF)
         IF(BLDATA(I,J,2).GT.TDIURN(I,J,6)) TDIURN(I,J,6)=BLDATA(I,J,2)
C        IF(BLDATA(I,J,2).GT.AIJ(I,J,76)) AIJ(I,J,76)=BLDATA(I,J,2)
C        IF(BLDATA(I,J,2).LT.AIJ(I,J,77)) AIJ(I,J,77)=BLDATA(I,J,2)
      IF(PEARTH.LE.0.)  GO TO 6000
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
         IF(JR.EQ.24) GO TO 5700
         DJ(JR,9)=DJ(JR,9)+TRHDTS*DXYPJ
         DJ(JR,13)=DJ(JR,13)+SHDTS*DXYPJ
         DJ(JR,14)=DJ(JR,14)+EVHDTS*DXYPJ
         DJ(JR,19)=DJ(JR,19)+EVAPS*DXYPJ
         DJ(JR,40)=DJ(JR,40)+AERUNS*PEARTH*DXYPJ
         DJ(JR,45)=DJ(JR,45)+DIFS*PEARTH*DXYPJ
         DJ(JR,47)=DJ(JR,47)+ARUNU*PEARTH*DXYPJ
         DJ(JR,48)=DJ(JR,48)+AERUNU*PEARTH*DXYPJ
         DJ(JR,54)=DJ(JR,54)+ARUNS*PEARTH*DXYPJ
         IF(MODDSF.NE.0) GO TO 5700
         DJ(JR,23)=DJ(JR,23)+(TSS-TFS)*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
 5700    AIJ(I,J,4)=AIJ(I,J,4)+SHDTS
         AIJ(I,J,7)=AIJ(I,J,7)+BETAD/NSURF
         IF(MODRD.EQ.0) AIJ(I,J,21)=AIJ(I,J,21)+TRHDTS/DTSRCE
         AIJ(I,J,22)=AIJ(I,J,22)+(SRHDTS+TRHDTS)
         AIJ(I,J,23)=AIJ(I,J,23)+(SRHDTS+TRHDTS+SHDTS+EVHDTS)
         IF(MODDSF.NE.0) GO TO 5800
         AIJ(I,J,34)=AIJ(I,J,34)+WSS                ! added 3/3/95 -rar-
         AIJ(I,J,35)=AIJ(I,J,35)+(TSS-TFS)
         AIJ(I,J,36)=AIJ(I,J,36)+USS
         AIJ(I,J,37)=AIJ(I,J,37)+VSS
         AIJ(I,J,47)=AIJ(I,J,47)+RTAUS
         AIJ(I,J,48)=AIJ(I,J,48)+RTAUUS
         AIJ(I,J,49)=AIJ(I,J,49)+RTAUVS
         AIJ(I,J,51)=AIJ(I,J,51)+QSS
CHYD     AIJ(I,J,53)=AIJ(I,J,53)
CHYD *  +   (40.6*PSOIL+.72*(2.*(TSS-TFS)-(QSATSS-QSS)*LHE/SHA))
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAG6
 5800    IF(MODD6.NE.0) GO TO 6000
         DO 5820 KR=1,4
         IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) GO TO 5840
 5820    CONTINUE
         GO TO 6000
 5840    ADAILY(IHOUR,12,KR)=ADAILY(IHOUR,12,KR)+TSS
         ADAILY(IHOUR,13,KR)=ADAILY(IHOUR,13,KR)+(TG1S+TFS)
         ADAILY(IHOUR,19,KR)=ADAILY(IHOUR,19,KR)+QSS
         ADAILY(IHOUR,20,KR)=ADAILY(IHOUR,20,KR)+QGS
         ADAILY(IHOUR,28,KR)=ADAILY(IHOUR,28,KR)+SRHDTS
         ADAILY(IHOUR,29,KR)=ADAILY(IHOUR,29,KR)+TRHDTS
         ADAILY(IHOUR,30,KR)=ADAILY(IHOUR,30,KR)+SHDTS
         ADAILY(IHOUR,31,KR)=ADAILY(IHOUR,31,KR)+EVHDTS
         ADAILY(IHOUR,32,KR)=ADAILY(IHOUR,32,KR)
     *       +SRHDTS+TRHDTS+SHDTS+EVHDTS
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
