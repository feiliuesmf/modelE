C**** SE001M12 E001M12 SOMTQ SB394M12
C****
C**** NEW CORRECTED PBL (Feb.19, 1998)
C****
C**** SURFCE for second order closure BL - Greg Hartke
C**** Changes for constant pressure above LS1 in Dry Conv.
C****   Boundary layer properties computed in subroutine PBL.
C**** SNOW AGING  - 11 vegetation types
C**** C version changes the way in which melting and freezing
C****  ice is handled in the boundary layer.
C**** D version corrects computation of TKV
C**** J version restores fluxes into first model to Model 2 method.
      SUBROUTINE SURFCE
C****
C**** THIS SUBROUTINE CALCULATES THE SURFACE FLUXES WHICH INCLUDE
C**** SENSIBLE HEAT, EVAPORATION, THERMAL RADIATION, AND MOMENTUM
C**** DRAG.  IT ALSO CALCULATES INSTANTANEOUS SURFACE TEMPERATURE,
C**** SURFACE SPECIFIC HUMIDITY, AND SURFACE WIND COMPONENTS.
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,rvap,stbo,bygrav,by6,byshi
     *     ,byrhoi,deltx,byrt3
      USE MODEL_COM, only : im,jm,lm,fim,DTsrc,NIsurf,u,v,t,p,q
     *     ,idacc,dsig,jday,ndasf,jeq,fland,flice,focean
     *     ,fearth,nday,modrd,ijd6,ITime,JHOUR,sige,byim,itocean
     *     ,itoice,itlake,itlkice,itlandi,aturb_on,vt_on
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : dxyp,imaxj,kmaxj,ravj,idij,idjj,siniv,cosiv
      USE RADNCB, only : trhr,fsf,cosz1
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs
     &     ,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     &     ,uflux,vflux,tflux,qflux
      USE SOCPBL, only : zgs
      USE DAGCOM, only : oa,aij,tdiurn,aj,areg,adaily,ndlypt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,j_tsrf,j_evap,j_evhdt,j_shdt,j_trhdt,j_f1dt
      USE DYNAMICS, only : pmid,pk,pedn,pek,pdsig,plij
      USE LANDICE, only : hc2li,z1e,z2li,hc1li
      USE LANDICE_COM, only : snowli
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : xsi,z1i,ace1i,hc1i,alami,byrli,byrls,rhos,kiext
     *     ,ksext 
      USE LAKES_COM, only : mwl,mldlk,gml,flake
      USE LAKES, only : minmld
      USE FLUXES, only : dth1,dq1,du1,dv1,e0,e1,evapor,runoe,erunoe
     *     ,solar,dmua,dmva,gtemp,nstype
      USE SOIL_DRV, only: earth
      IMPLICIT NONE

      INTEGER I,J,K,IM1,IP1,KR,JR,NS,NSTEPS,MODDSF,MODD6
     *     ,KMAX,IMAX,ITYPE,NGRNDZ,NG
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1,P1K,H0M1,HZM1
     *     ,PGK,HS,PKDN,DXYPJ,BETAS,EVHDTS,CDMS,CDHS,CDQS,EDS1S,PPBLS
     *     ,EVAPS,DBLS,BETA,ELHX,ACE2,CDTERM,CDENOM,HC1,dF1dTG,HCG1,HCG2
     *     ,DTGRND,EVHDT,F1DT,CM,CH,CQ,DHGS,DQGS,DGS,BETAUP,EVHEAT,F0,F1
     *     ,DSHDTG,DQGDTG,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG,dEVdQS
     *     ,HSDEN,HSCON,HSMUL,dHS,dQS,dT2,dTS,DQ1X,EVHDT0,EVAP,F0DT
     *     ,FTEVAP,VAP,TIMEZ,PWATER,PXSOIL,PSK,TH1,Q1,THV1,TFS
     *     ,RMBYA,HZM1,Q0M1,QZM1,TSS,QSS,TAUS,RTAUS,RTAUUS,RTAUVS,TG1S
     *     ,QGS,SRHDTS,TRHDTS,SHDTS,UGS,PTYPE,TG1,SRHEAT,SNOW,TG2,SHDT
     *     ,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,SHEAT,TRHEAT,QSDEN
     *     ,QSCON,QSMUL,T2DEN,T2CON,T2MUL,TGDEN,FQEVAP,ZS1CO,USS
     *     ,VSS,WSS,VGS,WGS,USRS,VSRS,Z2,Z2BY4L,Z1BY6L,THZ1,QZ1,POC,POI
     *     ,PLK,PLKI,EVAPLIM,F1DTS,HICE,HSNOW,HICE1,HSNOW1,F2,FSRI(2)
     *     ,rvx,PSIS

      REAL*8 MSUM, MA1, MSI1, MSI2
      REAL*8, DIMENSION(NSTYPE,IM,JM) :: TGRND,TGRN2

C**** Interface to PBL
      REAL*8 ZS1,TGV,TKV,QG,HEMI,DTSURF,US,VS,WS,TSV,QS,PSI,DBL,KM,KH,
     *     PPBL,UG,VG,WG,ZMIX
      LOGICAL POLE
      COMMON /PBLPAR/ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE

      COMMON /PBLOUT/US,VS,WS,TSV,QS,PSI,DBL,KM,KH,PPBL,UG,VG,WG,ZMIX

      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: S1BYG1 = BYRT3, 
     *     Z1IBYL=Z1I/ALAMI, Z2LI3L=Z2LI/(3.*ALAMI), Z1LIBYL=Z1E/ALAMI
      REAL*8 QSAT,DQSATDT,TOFREZ
C****
C**** ZATMO     GEOPOTENTIAL (G*Z)
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** GTEMP(1)  GROUND TEMPERATURE ARRAY OVER ALL SURFACE TYPES (C)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****   MSI  OCEAN ICE AMOUNT OF SECOND LAYER (KG/M**2)
C****
      if(.not. vt_on) then
          rvx=0.
      else
          rvx=deltx
      endif

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf

      ZS1CO=.5*DSIG(1)*RGAS*BYGRAV

C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND
      DO J=1,JM
      DO I=1,IM
        TGRND(2,I,J)=GTEMP(1,2,I,J)
        TGRND(3,I,J)=GTEMP(1,3,I,J)
C       TGRND(4,I,J)=GTEMP(1,4,I,J)
        TGRN2(2,I,J)=GTEMP(2,2,I,J)
        TGRN2(3,I,J)=GTEMP(2,3,I,J)
      END DO
      END DO

C**** Zero out fluxes summed over type
      E0=0. ; E1=0. ; EVAPOR=0. ; RUNOE=0. ; ERUNOE=0.
      DMUA=0. ; DMVA=0. ; SOLAR=0.

C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NIsurf TIMES EVERY HOUR
C****
      DO NS=1,NIsurf
         MODDSF=MOD(NSTEPS+NS-1,NDASF*NIsurf+1)
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1
         MODD6=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??
         TIMEZ=JDAY+(MOD(Itime,nday)+(NS-1.)/NIsurf)/NDAY ! -1 ??
         IF(JDAY.LE.31) TIMEZ=TIMEZ+365.
C**** ZERO OUT LAYER 1 WIND INCREMENTS
      DTH1=0. ;  DQ1 =0. ;  DU1=0. ; DV1=0.

      call loadbl
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO 7000 J=1,JM
         KMAX=KMAXJ(J)
         IMAX=IMAXJ(J)
      HEMI=1.
      IF(J.LE.JM/2) HEMI=-1.
      POLE=.FALSE.
      IF(J.EQ.1 .or. J.EQ.JM) POLE = .TRUE.
      IM1=IM
      DO I=1,IMAX

      ! until pbl loops over i,j,itype
      WSAVG(I,J)=0.
      TSAVG(I,J)=0.
      QSAVG(I,J)=0.
      USAVG(I,J)=0.
      VSAVG(I,J)=0.
      TAUAVG(I,J)=0.

      ! initialize fluxes before calling PBL subroutine
      uflux(I,J)=0.
      vflux(I,J)=0.
      tflux(I,J)=0.
      qflux(I,J)=0.
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      POICE=RSI(I,J)*PWATER
      POCEAN=PWATER-POICE
      POC =FOCEAN(I,J)*(1.-RSI(I,J))
      POI =FOCEAN(I,J)*    RSI(I,J)
      PLK = FLAKE(I,J)*(1.-RSI(I,J))
      PLKI= FLAKE(I,J)*    RSI(I,J)
      PXSOIL=POCEAN+POICE+PLICE
      PIJ=P(I,J)
      PS=PEDN(1,I,J)    ! PIJ+PTOP
      PSK=PEK(1,I,J)    ! EXPBYK(PS)
      P1=PMID(1,I,J)    ! SIG(1)*PIJ+PTOP
      P1K=PK(1,I,J)     ! EXPBYK(P1)
      TH1=T(I,J,1)
      Q1=Q(I,J,1)
      THV1=TH1*(1.+Q1*RVX)
         TFS=TF*PXSOIL
      RMBYA=100.*PDSIG(1,I,J)/GRAV
C      THZ1 = TZ(I,J,1) ! vertical gradient of potential temperature
C      QZ1 = QZ(I,J,1) ! vertical gradient of specific humidity
      MSUM = (PS*100.)/GRAV ! total column mass of atmosphere (kg/m^2)
      MA1 = RMBYA ! mass of lowest atmospheric layer (kg/m^2)
      H0M1 = TH1*SHA*MA1*DXYP(J) ! mean pot.enthalpy of lowest atm. (J)
      HZM1 = THZ1*SHA*MA1*DXYP(J) ! vert. grad. of lowest pot. enth.(J)
      Q0M1 = Q1*MA1*DXYP(J) ! mean water vapor of lowest atmosphere (kg)
      QZM1 = QZ1*MA1*DXYP(J) ! vert. grad. of lowest layer  vapor (kg)
      PGK = (PS*100.)**KAPA
      HS = (H0M1-HZM1*S1BYG1)*PGK/(DXYP(J)*MA1) ! pot. spec. enth.(J/kg)
      PKDN = (GRAV*(MSUM-MA1*0.25))**KAPA
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
         JR=JREG(I,J)
         DXYPJ=DXYP(J)
         TG1S=0.
         QGS=0.
         BETAS=0.
         SRHDTS=0.
         TRHDTS=0.
         SHDTS=0.
         EVHDTS=0.
         UGS=0.
         VGS=0.
         WGS=0.
         USRS=0.
         VSRS=0.
         CDMS=0.
         CDHS=0.
         CDQS=0.
         EDS1S=0.
         PPBLS=0.
         EVAPS=0.
         F1DTS=0.
         DBLS=0.
         PSIS=0.
C****
      IF (POCEAN.LE.0.) then
        ipbl(i,j,1)=0.
        GO TO 2200
      endif
C****
C**** OCEAN
C****
      ITYPE=1
      PTYPE=POCEAN
      NGRNDZ=1
      TG1=GTEMP(1,1,I,J)
      IF (FLAKE(I,J).gt.0) THEN
C**** limit evap if between MINMLD and 40cm, no evap below 40cm
        IF (MWL(I,J).lt.MINMLD*RHOW*FLAKE(I,J)*DXYP(J)) THEN
          EVAPLIM=MAX(0.5*(MWL(I,J)/(FLAKE(I,J)*DXYP(J))-0.4d0*RHOW),
     *         0d0)
        ELSE
          EVAPLIM=MWL(I,J)/(FLAKE(I,J)*DXYP(J))-(0.5*MINMLD+0.2d0)*RHOW
        END IF
      END IF
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(1,I,J)=SOLAR(1,I,J)+DTSURF*SRHEAT
            OA(I,J,5)=OA(I,J,5)+SRHEAT*DTSURF
      BETA=1.
      ELHX=LHE
      GO TO 3000
C****
 2200 IF (POICE.LE.0.) then
        ipbl(i,j,2)=0.
        GO TO 2400
      endif
C****
C**** OCEAN ICE
C****
      ITYPE=2
      PTYPE=POICE
      NGRNDZ=1    ! NIgrnd>1 currently not an option
      SNOW=SNOWI(I,J)
      TG1=TGRND(2,I,J)
      TG2=TGRN2(2,I,J)
      ACE2=MSI(I,J)
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(2,I,J) = 0.
      SOLAR(2,I,J)=SOLAR(2,I,J)+DTSURF*SRHEAT
C**** fraction of solar radiation leaving layer 1 and 2
      IF (SRHEAT.gt.0) THEN ! don't bother if there is no sun
        HSNOW = SNOW/RHOS
        HICE  = ACE1I*BYRHOI
        IF (ACE1I*XSI(1).gt.SNOW*XSI(2)) THEN 
C**** first thermal layer is snow and ice
          HICE1 = (ACE1I-XSI(2)*(SNOW+ACE1I))*BYRHOI
	  FSRI(1) = EXP(-KSEXT*HSNOW-KIEXT*HICE1)
        ELSE ! all snow			   
          HSNOW1 = (ACE1I+SNOW)*XSI(1)/RHOS
          FSRI(1) = EXP(-KSEXT*HSNOW1)
        END IF
        FSRI(2) = EXP(-KSEXT*HSNOW-KIEXT*HICE)
      ELSE
        FSRI(1:2) = 0
      END IF
            OA(I,J,12)=OA(I,J,12)+SRHEAT*DTSURF
      Z2=ACE2/RHOI
      Z2BY4L=Z2/(4.*ALAMI)
      Z1BY6L=(Z1IBYL+SNOW*BYRLS)*BY6
      CDTERM=1.5*TG2-.5*TOFREZ(I,J)
      CDENOM=1./(2.*Z1BY6L+Z2BY4L)
      HC1=HC1I+SNOW*SHI
      BETA=1.
      ELHX=LHS
C*
      MSI1 = SNOW+ACE1I ! snow and first layer ice mass (kg/m^2)
      MSI2 = ACE2 ! second (physical) layer ice mass (kg/m^2)
      dF1dTG = 2./(ACE1I*BYRLI+SNOW*BYRLS)
      HCG1 = SHI*XSI(1)*MSI1 ! heat capacity of top ice layer (J/C*m^2)
      HCG2 = SHI*XSI(2)*MSI1 ! heat capacity of second layer ice
      GO TO 3000
C****
 2400 IF (PLICE.LE.0.) then
        ipbl(i,j,3)=0.
        GO TO 5000
      endif
      NGRNDZ=1    ! NIgrnd>1 currently not an option
C****
C**** LAND ICE
C****
      ITYPE=3
      PTYPE=PLICE
      SNOW=SNOWLI(I,J)
      TG1=TGRND(3,I,J)
      TG2=GTEMP(2,3,I,J)
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
      CDTERM=TG2
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)
      HCG1=HC1LI+SNOW*SHI
      HC1 = HCG1
      BETA=1.
      ELHX=LHS
C****
C**** BOUNDARY LAYER INTERACTION
C****
 3000 TKV=THV1*PSK
      ZS1=ZS1CO*TKV*PIJ/PS
      DTGRND=DTSURF/NGRNDZ
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.
C**********************************************************************
C**** LOOP OVER GROUND TIME STEPS *************************************
C**********************************************************************
      DO 3600 NG=1,NGRNDZ
      TG=TG1+TF
      QG=QSAT(TG,ELHX,PS)
      TGV=TG*(1.+QG*RVX)
C =====================================================================
      CALL PBL(I,J,ITYPE,PTYPE)
      CM = cmgs(i,j,itype)
      CH = chgs(i,j,itype)
      CQ = cqgs(i,j,itype)
      DHGS=(ZMIX-ZGS)*CH*WS
      DQGS=(ZMIX-ZGS)*CQ*WS
      DGS =DQGS
C =====================================================================
      TS=TSV/(1.+QS*RVX)
C**** CALCULATE RHOS*CM*WS AND RHOS*CH*WS
3500  CONTINUE
      RHOSRF=100.*PS/(RGAS*TSV)
      RCDMWS=CM*WS*RHOSRF
      RCDHWS=CH*WS*RHOSRF
      RCDQWS=CQ*WS*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)
      SHEAT=SHA*RCDHWS*(TS-TG)
      BETAUP = BETA
      IF (QS .GT. QG) BETAUP = 1.
      EVHEAT=(LHE+TG1*SHV)*BETAUP*RCDQWS*(QS-QG)
      TRHEAT=TRHR(1,I,J)-STBO*(TG*TG)*(TG*TG)
      IF(ITYPE.EQ.1) GO TO 3620
C**** CALCULATE FLUXES USING IMPLICIT TIME STEP FOR NON-OCEAN POINTS
      IF (ITYPE .EQ. 2) GO TO 3550
C*
      F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
      F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
      DSHDTG=-RCDHWS*KH*SHA/(DHGS+KH)
      DQGDTG=QG*DQSATDT(TG,ELHX)
      DEVDTG=-RCDQWS*KH*LHE*BETAUP*DQGDTG/(DQGS+KH)
      DTRDTG=-4.*STBO*TG*TG*TG
      DF0DTG=DSHDTG+DEVDTG+DTRDTG
      DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
      DTG=(F0-F1)*DTGRND/(HC1-DTGRND*DFDTG)
      SHDT=SHDT+DTGRND*(SHEAT+DTG*DSHDTG)
      EVHDT=EVHDT+DTGRND*(EVHEAT+DTG*DEVDTG)
      TRHDT=TRHDT+DTGRND*(TRHEAT+DTG*DTRDTG)
      F1DT=F1DT+DTGRND*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
      DU1(I,J)=DU1(I,J)+PTYPE*DTGRND*RCDMWS*US/RMBYA
      DV1(I,J)=DV1(I,J)+PTYPE*DTGRND*RCDMWS*VS/RMBYA
      TG1=TG1+DTG
      GO TO 3600
C*
 3550 CONTINUE
! heat flux on first/second/third layers (W/m^2)
      F1 = (TG1-TG2)*dF1dTG + SRHEAT*FSRI(1)
      F2 = SRHEAT*FSRI(2)
      EVHEAT = LHE*RCDQWS*(QS-QG) ! latent heat flux (W/m^2)
      F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
      dSNdTG=-RCDHWS*KH*SHA/(DHGS+KH)
      dQGdTG=QG*DQSATDT(TG,ELHX) ! d(QG)/dTG
      dEVdTG = -dQGdTG*LHE*RCDQWS*KH/(DGS+KH) ! d(EVHEAH)/dTG
      dTRdTG = -4*STBO*TG*TG*TG ! d(TRHEAT)/dTG
      dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG
C     dSNdHS = RCDHWS ! d(SHEAT)/dHS - kg/(sec*m^2)
      dEVdQS = LHE*RCDQWS ! d(EVHEAH)/dQS
      HSDEN = -(1.+2.*S1BYG1)*DTGRND*PGK*BETA*dSNdTG+MA1*PKDN*SHA
      HSCON = -(1.+2.*S1BYG1)*DTGRND*PGK*SHEAT/HSDEN ! (J*sec)/kg
      HSMUL = -(1.+2.*S1BYG1)*DTGRND*PGK*BETA*dSNdTG/HSDEN ! J/(kg*degC)
      QSDEN = (1.+2.*S1BYG1)*BETA*DTGRND*dEVdQS+MA1*LHE
      QSCON = -(1.+2.*S1BYG1)*DTGRND*EVHEAT/QSDEN
      QSMUL = -(1.+2.*S1BYG1)*DTGRND*BETA*dEVdTG/QSDEN
      T2DEN = HCG2+BETA*DTGRND*dF1dTG
      T2CON = DTGRND*(F1-F2)/T2DEN
      T2MUL = BETA*DTGRND*dF1dTG/T2DEN
      TGDEN = HCG1-BETA*DTGRND*(dF0dTG-dF1dTG-
     A        HSMUL*dSNdTG+QSMUL*dEVdQS+T2MUL*dF1dTG) ! W/(m^2*degC)
      dTG = DTGRND*(F0-F1+BETA*
     A      (QSCON*dEVdQS-HSCON*dSNdTG+T2CON*dF1dTG))/TGDEN ! degC
      IF (TG1+dTG .GT. 0.) dTG = -TG1
      dHS = HSCON+HSMUL*dTG ! (J*sec)/kg
      dQS = QSCON+QSMUL*dTG
      dT2 = T2CON+T2MUL*dTG
      SHDT = DTGRND*(SHEAT+BETA*((dTG-dHS)*dSNdTG)) ! sensible
      EVHDT = DTGRND*(EVHEAT+BETA*(dTG*dEVdTG+dQS*dEVdQS)) ! latent
      TRHDT = DTGRND*(TRHEAT+BETA*dTG*dTRdTG) ! thermal flux (J/m^2)
      F1DT = DTGRND*(F1+BETA*(dTG*dF1dTG-dT2*dF1dTG))
      DU1(I,J)=DU1(I,J)+PTYPE*DTGRND*RCDMWS*US/RMBYA
      DV1(I,J)=DV1(I,J)+PTYPE*DTGRND*RCDMWS*VS/RMBYA
      TG1 = TG1+dTG ! first layer sea ice temperature (degC)
      TG2 = TG2+dT2 ! second layer sea ice temperature (degC)
      TGRN2(ITYPE,I,J) = TG2
 3600 CONTINUE
      GO TO 3700
C**** CALCULATE FLUXES USING IMPLICIT TIME STEP ALSO FOR OCEAN POINTS
 3620 CONTINUE
      DSHDTG=-RCDHWS*SHA
      dEVdQS = LHE*RCDQWS
      dHS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A      ((1.+2.*S1BYG1)*DTSURF*PGK*RCDHWS+MA1*PKDN)
      dTS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A      (MA1*PKDN*SHA-(1.+2.*S1BYG1)*DTSURF*PGK*DSHDTG)
      dQS = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/
     A      ((1.+2.*S1BYG1)*DTSURF*dEVdQS+MA1*LHE)
      SHDT = DTSURF*(SHEAT-dTS*DSHDTG)
      EVHDT=DTSURF*(EVHEAT+dQS*dEVdQS) ! latent heat flux
      TRHDT=DTSURF*TRHEAT
      DU1(I,J)=DU1(I,J)+PTYPE*DTGRND*RCDMWS*US/RMBYA
      DV1(I,J)=DV1(I,J)+PTYPE*DTGRND*RCDMWS*VS/RMBYA
C**** CALCULATE EVAPORATION
 3700 CONTINUE
      DQ1X =EVHDT/((LHE+TG1*SHV)*RMBYA)
      EVHDT0=EVHDT
C**** Limit evaporation if lake mass is at minimum 
      IF (ITYPE.EQ.1 .and. PLK.GT.0 .and.
     *     (EVAPOR(I,J,1)-DQ1X*RMBYA).gt.EVAPLIM) THEN
        WRITE(99,*) "Lake EVAP limited: I,J,EVAP,MWL",I,J,EVAPOR(I,J,1)
     *       -DQ1X*RMBYA, MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J))
        DQ1X=(EVAPOR(I,J,1)-EVAPLIM)/RMBYA
      ELSEIF (DQ1X.GT.Q1+DQ1(I,J)) THEN
        DQ1X=(Q1+DQ1(I,J))
      ELSE
        GO TO 3720
      END IF
      EVHDT=DQ1X*(LHE+TG1*SHV)*RMBYA
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*RMBYA
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP
      TGRND(ITYPE,I,J)=TG1
      DTH1(I,J)=DTH1(I,J)-SHDT*PTYPE/(SHA*RMBYA*P1K)
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
      DMUA(I,J,ITYPE)=DMUA(I,J,ITYPE)+PTYPE*DTGRND*RCDMWS*US
      DMVA(I,J,ITYPE)=DMVA(I,J,ITYPE)+PTYPE*DTGRND*RCDMWS*VS
      USS=USS+US*PTYPE
      VSS=VSS+VS*PTYPE
      WSS=WSS+WS*PTYPE
      TSS=TSS+TS*PTYPE
      QSS=QSS+QS*PTYPE
      TAUS=TAUS+CM*WS*WS*PTYPE
         RTAUS=RTAUS+RCDMWS*WS*PTYPE
         RTAUUS=RTAUUS+RCDMWS*US*PTYPE
         RTAUVS=RTAUVS+RCDMWS*VS*PTYPE
         TG1S=TG1S+TG1*PTYPE
         QGS=QGS+QG*PTYPE
         BETAS=BETAS + BETA*PTYPE
         SRHDTS=SRHDTS+SRHEAT*DTSURF*PTYPE
         TRHDTS=TRHDTS+TRHDT*PTYPE
         SHDTS=SHDTS+SHDT*PTYPE
         EVHDTS=EVHDTS+EVHDT*PTYPE
         UGS=UGS+UG*PTYPE
         VGS=VGS+VG*PTYPE
         WGS=WGS+WG*PTYPE
         CDMS=CDMS+CM*PTYPE
         CDHS=CDHS+CH*PTYPE
         CDQS=CDQS+CQ*PTYPE
         EDS1S=EDS1S+KH*PTYPE
         PPBLS=PPBLS+PPBL*PTYPE
         EVAPS=EVAPS+EVAP*PTYPE
         F1DTS=F1DTS+F1DT*PTYPE
         DBLS=DBLS+DBL*PTYPE
         PSIS=PSIS+PSI*PTYPE
5666  GO TO (4000,4100,4400),ITYPE
C****
C**** OCEAN
C****
 4000 CONTINUE
         AJ(J,J_EVHDT,ITOCEAN)=AJ(J,J_EVHDT,ITOCEAN)+EVHDT  *POC
         AJ(J,J_SHDT ,ITOCEAN)=AJ(J,J_SHDT ,ITOCEAN)+SHDT   *POC
         AJ(J,J_TRHDT,ITOCEAN)=AJ(J,J_TRHDT,ITOCEAN)+TRHDT  *POC
         AJ(J,J_TRHDT,ITLAKE) =AJ(J,J_TRHDT,ITLAKE) +TRHDT  *PLK
         AJ(J,J_SHDT ,ITLAKE) =AJ(J,J_SHDT ,ITLAKE) +SHDT   *PLK
         AJ(J,J_EVHDT,ITLAKE) =AJ(J,J_EVHDT,ITLAKE) +EVHDT  *PLK
         IF(MODDSF.EQ.0) THEN
           AJ(J,J_TSRF,ITOCEAN)=AJ(J,J_TSRF,ITOCEAN)+(TS-TF)*POC
           AJ(J,J_TSRF,ITLAKE) =AJ(J,J_TSRF,ITLAKE) +(TS-TF)*PLK
         END IF
            OA(I,J,6)=OA(I,J,6)+TRHDT
            OA(I,J,7)=OA(I,J,7)+SHDT
            OA(I,J,8)=OA(I,J,8)+EVHDT
      GO TO 2200
C****
C**** OCEAN ICE
C****
 4100    CONTINUE
         AJ(J,J_EVHDT,ITOICE) =AJ(J,J_EVHDT,ITOICE) +EVHDT  *POI
         AJ(J,J_SHDT ,ITOICE) =AJ(J,J_SHDT ,ITOICE) +SHDT   *POI
         AJ(J,J_TRHDT,ITOICE) =AJ(J,J_TRHDT,ITOICE) +TRHDT  *POI
         AJ(J,J_TRHDT,ITLKICE)=AJ(J,J_TRHDT,ITLKICE)+TRHDT  *PLKI
         AJ(J,J_SHDT ,ITLKICE)=AJ(J,J_SHDT ,ITLKICE)+SHDT   *PLKI
         AJ(J,J_EVHDT,ITLKICE)=AJ(J,J_EVHDT,ITLKICE)+EVHDT  *PLKI
         IF(MODDSF.EQ.0) THEN
           AJ(J,J_TSRF,ITOICE) =AJ(J,J_TSRF,ITOICE) +(TS-TF)*POI
           AJ(J,J_TSRF,ITLKICE)=AJ(J,J_TSRF,ITLKICE)+(TS-TF)*PLKI
         END IF
         IF (TG1.GT.TDIURN(I,J,7)) TDIURN(I,J,7) = TG1
            OA(I,J,9)=OA(I,J,9)+TRHDT
            OA(I,J,10)=OA(I,J,10)+SHDT
            OA(I,J,11)=OA(I,J,11)+EVHDT
      GO TO 2400
C****
C**** LAND ICE
C****
 4400    CONTINUE
         AJ(J,J_EVHDT,ITLANDI) =AJ(J,J_EVHDT,ITLANDI) +EVHDT  *PLICE
         AJ(J,J_SHDT ,ITLANDI) =AJ(J,J_SHDT ,ITLANDI) +SHDT   *PLICE
         AJ(J,J_TRHDT,ITLANDI) =AJ(J,J_TRHDT,ITLANDI) +TRHDT  *PLICE
         IF(MODDSF.EQ.0) AJ(J,J_TSRF,ITLANDI)=AJ(J,J_TSRF,ITLANDI)
     *        +(TS-TF)*PLICE
         IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
         AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(TS-TF)
         AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)+SHDT
         AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)+EVHDT
         AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)+TRHDT
C**** NON-OCEAN POINTS WHICH ARE NOT MELTING OR FREEZING WATER USE
C****   IMPLICIT TIME STEPS
C****
C**** UPDATE SURFACE AND FIRST LAYER QUANTITIES
C****
 5000    CONTINUE
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
         IF(JR.EQ.24) GO TO 5700
         AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+TRHDTS*DXYPJ
         AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+SHDTS*DXYPJ
         AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+EVHDTS*DXYPJ
         AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+EVAPS*DXYPJ
         AREG(JR,J_F1DT )=AREG(JR,J_F1DT )+F1DTS*DXYPJ
         IF(MODDSF.NE.0) GO TO 5700
         AREG(JR,J_TSRF)=AREG(JR,J_TSRF)+(TSS-TFS)*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
 5700    AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDTS
         IF(MODRD.EQ.0) AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)+TRHDTS
     *        /DTSRC
         AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHDTS+TRHDTS)
         AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHDTS+TRHDTS+SHDTS+EVHDTS)
         IF(MODDSF.NE.0) GO TO 5800
         AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WSS        ! added 12/29/96 -rar-
         AIJ(I,J,IJ_TS)=AIJ(I,J,IJ_TS)+(TSS-TFS)
         AIJ(I,J,IJ_US)=AIJ(I,J,IJ_US)+USS
         AIJ(I,J,IJ_VS)=AIJ(I,J,IJ_VS)+VSS
         AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RTAUS
         AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RTAUUS
         AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RTAUVS
         AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QSS
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAG6
 5800    IF(MODD6.EQ.0) THEN
         DO KR=1,NDLYPT
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(JHOUR+1,6,KR)=ADAILY(JHOUR+1,6,KR)+PS  ! ???
               ADAILY(JHOUR+1,7,KR)=ADAILY(JHOUR+1,7,KR)+PSK*T(I,J,5)
               ADAILY(JHOUR+1,8,KR)=ADAILY(JHOUR+1,8,KR)+PSK*T(I,J,4)
               ADAILY(JHOUR+1,9,KR)=ADAILY(JHOUR+1,9,KR)+PSK*T(I,J,3)
               ADAILY(JHOUR+1,10,KR)=ADAILY(JHOUR+1,10,KR)+PSK*T(I,J,2)
               ADAILY(JHOUR+1,11,KR)=ADAILY(JHOUR+1,11,KR)+PSK*T(I,J,1)
               ADAILY(JHOUR+1,12,KR)=ADAILY(JHOUR+1,12,KR)+TSS
               ADAILY(JHOUR+1,13,KR)=ADAILY(JHOUR+1,13,KR)+(TG1S+TFS)
               ADAILY(JHOUR+1,14,KR)=ADAILY(JHOUR+1,14,KR)+Q(I,J,5)
               ADAILY(JHOUR+1,15,KR)=ADAILY(JHOUR+1,15,KR)+Q(I,J,4)
               ADAILY(JHOUR+1,16,KR)=ADAILY(JHOUR+1,16,KR)+Q(I,J,3)
               ADAILY(JHOUR+1,17,KR)=ADAILY(JHOUR+1,17,KR)+Q(I,J,2)
               ADAILY(JHOUR+1,18,KR)=ADAILY(JHOUR+1,18,KR)+Q1
               ADAILY(JHOUR+1,19,KR)=ADAILY(JHOUR+1,19,KR)+QSS
               ADAILY(JHOUR+1,20,KR)=ADAILY(JHOUR+1,20,KR)+QGS
               ADAILY(JHOUR+1,28,KR)=ADAILY(JHOUR+1,28,KR)+SRHDTS
               ADAILY(JHOUR+1,29,KR)=ADAILY(JHOUR+1,29,KR)+TRHDTS
               ADAILY(JHOUR+1,30,KR)=ADAILY(JHOUR+1,30,KR)+SHDTS
               ADAILY(JHOUR+1,31,KR)=ADAILY(JHOUR+1,31,KR)+EVHDTS
               ADAILY(JHOUR+1,32,KR)=ADAILY(JHOUR+1,32,KR)
     *              +(SRHDTS+TRHDTS+SHDTS+EVHDTS)
               ADAILY(JHOUR+1,33,KR)=ADAILY(JHOUR+1,33,KR)+UGS
               ADAILY(JHOUR+1,34,KR)=ADAILY(JHOUR+1,34,KR)+VGS
               ADAILY(JHOUR+1,35,KR)=ADAILY(JHOUR+1,35,KR)+WGS
               ADAILY(JHOUR+1,36,KR)=ADAILY(JHOUR+1,36,KR)+USS
               ADAILY(JHOUR+1,37,KR)=ADAILY(JHOUR+1,37,KR)+VSS
               ADAILY(JHOUR+1,38,KR)=ADAILY(JHOUR+1,38,KR)+WSS
               ADAILY(JHOUR+1,39,KR)=ADAILY(JHOUR+1,39,KR)+PSIS
               ADAILY(JHOUR+1,42,KR)=ADAILY(JHOUR+1,42,KR)+CDMS
               ADAILY(JHOUR+1,43,KR)=ADAILY(JHOUR+1,43,KR)+CDHS
               ADAILY(JHOUR+1,44,KR)=ADAILY(JHOUR+1,44,KR)+CDQS
               ADAILY(JHOUR+1,45,KR)=ADAILY(JHOUR+1,45,KR)+EDS1S
               ADAILY(JHOUR+1,46,KR)=ADAILY(JHOUR+1,46,KR)+DBLS
               ADAILY(JHOUR+1,50,KR)=ADAILY(JHOUR+1,50,KR)+EVAPS
            END IF
         END DO
       END IF
       IM1=I
      END DO

 7000 CONTINUE
C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODD6)
C****
C**** UPDATE FIRST LAYER QUANTITIES
C****
      DO J=1,JM
      IMAX=IMAXJ(J)
      DO I=1,IMAX
        FTEVAP=0
        IF (DTH1(I,J)*T(I,J,1).lt.0) FTEVAP=-DTH1(I,J)/T(I,J,1)
        FQEVAP=0
        IF (DQ1(I,J).lt.0.and.Q(I,J,1).gt.0) FQEVAP=-DQ1(I,J)/Q(I,J,1)
        T(I,J,1)=  T(I,J,1)+DTH1(I,J)
        Q(I,J,1)=  Q(I,J,1)+DQ1(I,J)
! Z-moments should be set from PBL
        TMOM(:,I,J,1) = TMOM(:,I,J,1)*(1.-FTEVAP)
        QMOM(:,I,J,1) = QMOM(:,I,J,1)*(1.-FQEVAP)
        IF (Q(I,J,1).LT.qmin) THEN
          WRITE(99,*) ITime,'I,J:',I,J,' Q1:',Q(I,J,1),'->',qmin
          Q(I,J,1)=qmin
          QMOM(:,I,J,1)=0.
        ENDIF
      END DO
      END DO
C****
C**** ADD IN SURFACE FRICTION TO FIRST LAYER WIND
C****
C**** Polar boxes
      DO J=1,JM,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        HEMI=1.
        IF(J.LE.JM/2) HEMI=-1.
        DO I=1,IMAX
        DO K=1,KMAX
          U(IDIJ(K,I,J),IDJJ(K,J),1)=U(IDIJ(K,I,J),IDJJ(K,J),1) -
     *           RAVJ(K,J)*(DU1(I,J)*COSIV(K)+DV1(I,J)*SINIV(K)*HEMI)
          V(IDIJ(K,I,J),IDJJ(K,J),1)=V(IDIJ(K,I,J),IDJJ(K,J),1) -
     *           RAVJ(K,J)*(DV1(I,J)*COSIV(K)-DU1(I,J)*SINIV(K)*HEMI)
        END DO
        END DO
      END DO
C**** non polar boxes
      DO J=2,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
        DO K=1,KMAX
          U(IDIJ(K,I,J),IDJJ(K,J),1)=U(IDIJ(K,I,J),IDJJ(K,J),1) -
     *           RAVJ(K,J)*DU1(I,J)
          V(IDIJ(K,I,J),IDJJ(K,J),1)=V(IDIJ(K,I,J),IDJJ(K,J),1) -
     *           RAVJ(K,J)*DV1(I,J)
        END DO
        END DO
      END DO
C****
      if(.not. aturb_on) then
C****     DRY CONVECTION ORIGINATING FROM THE FIRST LAYER
          CALL DRYCNV(1,1)
      else
C****     turbulence throughout whole atmosphere
          call diffus(dtsurf)
      endif
C****
C**** ACCUMULATE SOME ADDITIONAL BOUNDARY LAYER DIAGNOSTICS
C****
      IF(MODD6.EQ.0) THEN
        DO KR=1,4
C**** CHECK IF DRY CONV HAS HAPPENED FOR THIS DIAGNOSTIC
          IF(DCLEV(IJD6(1,KR),IJD6(2,KR)).GT.1.) THEN
            ADAILY(JHOUR+1,47,KR)=ADAILY(JHOUR+1,47,KR)+1.
            ADAILY(JHOUR+1,48,KR)=ADAILY(JHOUR+1,48,KR)
     *           +DCLEV(IJD6(1,KR),IJD6(2,KR))
          END IF
        END DO
      END IF
C****
      END DO
      RETURN
 9991 FORMAT ('0SURFACE ',4I4,5F10.4,3F11.7)
 9992 FORMAT ('0',I2,10F10.4/23X,4F10.4,10X,2F10.4/
     *  33X,3F10.4,10X,2F10.4)
 9993 FORMAT ('0',I2,10F10.4/23X,7F10.4/33X,7F10.4)
 9994 FORMAT ('0',I2,11F10.4)
      END SUBROUTINE SURFCE

      SUBROUTINE DRYCNV(LBASE_MIN,LBASE_MAX)
!@sum  DRYCNV mixes air caused by dry convection.
!@+    this version checks base layers lbase_min to lbase_max.
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : lhe,sha,deltx
      USE MODEL_COM
      USE GEOM
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
      USE DAGCOM, only : ajl
      USE DYNAMICS, only : pk,pdsig,plij
      USE PBLCOM, only : dclev
      IMPLICIT NONE

      integer, intent(in) :: LBASE_MIN,LBASE_MAX
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT
      REAL*8, DIMENSION(LM) :: DP
      COMMON/WORK2/UT,VT,DP
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      LOGICAL POLE
      INTEGER I,J,L,K,IMAX,KMAX,IM1,LMAX,LMIN

      DOUBLE PRECISION, DIMENSION(NMOM) :: TMOMS,QMOMS
      REAL*8 DOK,PIJBOT,PIJ,PKMS,THPKMS,QMS
     *     ,TVMS,THETA,RDP,THM
     *     ,rvx

      if(.not. vt_on) then
          rvx=0.
      else
          rvx=deltx
      endif

      if(LBASE_MAX.GE.LM) stop 'DRYCNV: LBASE_MAX.GE.LM'
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      UT=U ; VT=V
C**** OUTSIDE LOOPS OVER J AND I
      JLOOP: DO J=1,JM
      POLE=.FALSE.
      IF (J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.

      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)
C****
C**** MAIN LOOP
C****
      IM1=IM
      ILOOP: DO I=1,IMAX
         DO K=1,KMAX
            RA(K)=RAVJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
      LMAX=LBASE_MIN-1
      lbase_loop: do while(lmax.lt.lbase_max)
      LMIN=LMAX+1
      LMAX=LMIN
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)) cycle lbase_loop
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS
C**** MIX THROUGH TWO LOWER LAYERS
      PIJBOT=PLIJ(LMIN,I,J)
      DP(LMIN)=PDSIG(LMIN,I,J)
      PIJ=PLIJ(LMIN+1,I,J)
      DP(LMIN+1)=PDSIG(LMIN+1,I,J)
      PKMS=PK(LMIN,I,J)*DP(LMIN)+PK(LMIN+1,I,J)*DP(LMIN+1)
      THPKMS=T(I,J,LMIN)*(PK(LMIN,I,J)*DP(LMIN))
     *  +T(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMS=Q(I,J,LMIN)*DP(LMIN)+Q(I,J,LMIN+1)*DP(LMIN+1)
C**** sum moments to mix over unstable layers
      TMOMS(XYMOMS) =
     &     TMOM(XYMOMS,I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     &     TMOM(XYMOMS,I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMOMS(XYMOMS) =
     &     QMOM(XYMOMS,I,J,LMIN  )*(DP(LMIN  ))  +
     &     QMOM(XYMOMS,I,J,LMIN+1)*(DP(LMIN+1))
      IF (LMIN+1.GE.LM) GO TO 150
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO L=LMIN+2,LM
        IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 160
        PIJ=PLIJ(L,I,J)
        DP(L)=PDSIG(L,I,J)
        PKMS=PKMS+(PK(L,I,J)*DP(L))
        THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*DP(L))
        QMS=QMS+Q(I,J,L)*DP(L)
        TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*DP(L))
        TMOMS(XYMOMS) = TMOMS(XYMOMS) +
     &       TMOM(XYMOMS,I,J,L)*(PK(L,I,J)*DP(L))
        QMOMS(XYMOMS) = QMOMS(XYMOMS) +
     &       QMOM(XYMOMS,I,J,L)*DP(L)
        THETA=TVMS/PKMS
      END DO
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PIJBOT*SIGE(LMIN)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      DO L=LMIN,LMAX
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(L,I,J)*PLIJ(L,I,J)
         AJL(J,L,55)=AJL(J,L,55)+(QMS-Q(I,J,L))*PDSIG(L,I,J)*LHE/SHA
      T(I,J,L)=THM
      TMOM(XYMOMS,I,J,L)=TMOMS(XYMOMS)/PKMS
      TMOM(ZMOMS,I,J,L)=0.
      Q(I,J,L)=QMS
      QMOM(XYMOMS,I,J,L)=QMOMS(XYMOMS)*RDP
      QMOM(ZMOMS,I,J,L)=0.
      END DO
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      DO L=LMIN,LMAX
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*DP(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*DP(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      DO L=LMIN,LMAX
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c the following line gives bytewise different ajl
            AJL(IDJ(K),L,38)=AJL(IDJ(K),L,38)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PLIJ(L,I,J)*RA(K)
         ENDDO
      ENDDO
      enddo lbase_loop
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
      if(lbase_min.eq.1) then ! was called from surfce
         DCLEV(I,J)=LMAX
      endif
      IM1=I
      ENDDO ILOOP
      ENDDO JLOOP
      RETURN
      END SUBROUTINE DRYCNV

