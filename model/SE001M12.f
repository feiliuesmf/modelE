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
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,rvap
      USE E001M12_COM
      USE SOMTQ_COM
      USE GEOM
      USE SOCPBL, only : ipbl,pbl,omega2,zgs
      IMPLICIT REAL*8 (A-H,O-Z)
C*
      REAL*8 KM, KH, MSUM, MA1, MSI1, MSI2
      COMMON/WORK1d/COSZ1(IM,JM),DTH1(IM,JM),DQ1(IM,JM)
      COMMON/WORK2/UT(IM,JM,LM),VT(IM,JM,LM),DU1(IM,JM),
     *  DV1(IM,JM)
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      COMMON/WORK3/E0(IM,JM,4),E1(IM,JM,4),EVAPOR(IM,JM,4),
     *  TGRND(IM,JM,4),BLTEMP(IM,JM,8)
            COMMON/WORKO/OA(IM,JM,12)
      COMMON /CIRCLE/PI,RADIAN,DEGREE
      COMMON /RDATA/ROUGHL(IM,JM)
      DIMENSION SINI(IM),COSI(IM),               TGRN2(IM,JM,4)
      LOGICAL POLE

      COMMON /PBLPAR/ZS1,PIJ,PSK,TGV,TKV,THV1,QG,HEMI,
     2               DTSURF,JVPO,IM1,POLE

      COMMON /PBLOUT/US,VS,WS,TSV,QS,PSI,DBL,KM,KH,USTAR,PPBL,
     2               CM,CH,CQ,UG,VG,WG,ZMIX

      parameter (qmin=1.e-12)

C      common /dflux1/fluxu1,fluxv1,fluxt1,fluxq1

c      DATA RVAP/461.5/
      DATA ALAMI/2.1762/,STBO/.5672573E-7/
c ,TF/273.16/,,RHOW/1000./,RHOI/916.6/,SHV/0./,SHW/4185./,SHI/2060./,
      DATA TFO/-1.80/
      DATA Z1I/.1/,Z2LI/2.9/,Z1E/.1/,Z2E/4./,RHOS/300.0/,ALAMS/.35/,
     A     S1BYG1 /0.57735/, XSI1 /0.5/, XSI2 /0.5/
      QSAT(TM,PR,QLH)=3.797915*DEXP(QLH*(7.93252D-6-2.166847D-3/TM))/PR
      DATA IFIRST/1/
C****
C**** ZATMO     GEOPOTENTIAL (G*Z)
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
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
C****        8  EARTH TEMPERATURE OF SECOND LAYER (C)
C****        9  EARTH WATER OF SECOND LAYER (KG/M**2)
C****       10  EARTH ICE OF SECOND LAYER (KG/M**2)
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)
C****       14  LAND ICE TEMPERATURE OF SECOND LAYER (C)
C****       15  OCEAN ICE TEMPERATURE OF THIRD LAYER (C)
C****       16  OCEAN ICE TEMPERATURE OF FOURTH LAYER (C)
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
C**** ROUGHL    LOG(ZGS/ROUGHNESS LENGTH) (LOGARITHM TO BASE 10)
C****
      NSTEPS=NSURF*NSTEP/NDYN
      IF(IFIRST.NE.1) GO TO 30
      IFIRST=0
      PI=ACOS(-1.D0)
      DEGREE=180./PI
      RADIAN=PI/180.
c      CALL DREAD (19,ROUGHL,IM*JM,ROUGHL)
      CALL READT (19,0,ROUGHL,IM*JM,ROUGHL,1)
      REWIND 19
c      ZGS=10.
      OMEGA2=2.*OMEGA

      DTSURF=NDYN*DT/NSURF
      DTSRCE=DT*NDYN
c      SHA=RGAS/KAPA
C*
      SHCD = SHA ! specific heat capacity for dry air - J/(kg*degC)
      BYRLI = 1./(RHOI*ALAMI) ! (m^4*degC*sec)/(J*kg)
      BYRLS = 1./(RHOS*ALAMS) ! (m^4*degC*sec)/(J*kg)
C*
      RVX=0.
      ACE1I=Z1I*RHOI
      HC1I=ACE1I*SHI
      HC2LI=Z2LI*RHOI*SHI
      HC1DE=Z1E*1129950.
      HC2DE=Z2E*1129950.+3.5*.125*RHOW*3100.
      Z1IBYL=Z1I/ALAMI
      Z2LI3L=Z2LI/(3.*ALAMI)
      ZS1CO=.5*DSIG(1)*RGAS/GRAV
      DO 20 I=1,IM
      SINI(I)=SIN((I-1)*TWOPI/FIM)
   20 COSI(I)=COS((I-1)*TWOPI/FIM)
   30 S0=S0X*1367./RSDIST
         SPRING=-1.
         IF((JDAY.GE.32).AND.(JDAY.LE.212)) SPRING=1.
C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND
      DO 40 J=1,JM
      DO 40 I=1,IM
      TGRND(I,J,2)=GDATA(I,J,3)
      TGRND(I,J,3)=GDATA(I,J,13)
C     TGRND(I,J,4)=GDATA(I,J,4)
      TGRN2(I,J,2) = GDATA(I,J,7)
      TGRN2(I,J,3) = GDATA(I,J,14)
C*
C*
      DO 40 K=1,12
   40 E0(I,J,K)=0.
      DO 45 I=1,IM*JM*8
        BLTEMP(I,1,1)=0.
45    CONTINUE
         IHOUR=1.5+TOFDAY
      call pgrads1
C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NSURF TIMES EVERY HOUR
C****
      DO 9000 NS=1,NSURF
         MODDSF=MOD(NSTEPS+NS-1,NDASF)
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1
         MODD6=MOD(IDAY+NS,NSURF)
         TIMEZ=JDAY+(TOFDAY+(NS-1.)/NSURF)/24.
         IF(JDAY.LE.31) TIMEZ=TIMEZ+365.
C**** ZERO OUT LAYER 1 WIND INCREMENTS
      DO 60 J=1,JM
      DO 60 I=1,IM
      DTH1(I,J)=0.
      DQ1(I,J) =0.
      DU1(I,J)=0.
   60 DV1(I,J)=0.
      call loadbl
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO 7000 J=1,JM
         KMAX=KMAXJ(J)
         IMAX=IMAXJ(J)
      HEMI=1.
      IF(J.LE.JM/2) HEMI=-1.
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
  100    ATRHDT=0.
         BTRHDT=0.
         CTRHDT=0.
         ASHDT=0.
         BSHDT=0.
         CSHDT=0.
         AEVHDT=0.
         BEVHDT=0.
         CEVHDT=0.
         ATS=0.
         BTS=0.
         CTS=0.
         JEQ=1+JM/2
         IF(J.LT.JEQ) WARMER=-SPRING
         IF(J.GE.JEQ) WARMER=SPRING
      IM1=IM
      DO 6000 I=1,IMAX
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      POICE=ODATA(I,J,2)*PWATER
      POCEAN=PWATER-POICE
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
      RMBYA=100.*PIJ*DSIG(1)/GRAV
C*
C      THZ1 = TZ(I,J,1) ! vertical gradient of potential temperature
C      QZ1 = QZ(I,J,1) ! vertical gradient of specific humidity
      MSUM = (PS*100.)/GRAV ! total column mass of atmosphere (kg/m^2)
      MA1 = RMBYA ! mass of lowest atmospheric layer (kg/m^2)
      H0M1 = TH1*SHCD*MA1*DXYP(J) ! mean pot.enthalpy of lowest atm. (J)
      HZM1 = THZ1*SHCD*MA1*DXYP(J) ! vert. grad. of lowest pot. enth.(J)
      Q0M1 = Q1*MA1*DXYP(J) ! mean water vapor of lowest atmosphere (kg)
      QZM1 = QZ1*MA1*DXYP(J) ! vert. grad. of lowest layer  vapor (kg)
      PGK = (PS*100.)**KAPA
      HS = (H0M1-HZM1*S1BYG1)*PGK/(DXYP(J)*MA1) ! pot. spec. enth.(J/kg)
C      QS = (Q0M1-QZM1*S1BYG1)/(DXYP(J)*MA1)
      PKDN = (GRAV*(MSUM-MA1*0.25))**KAPA
C*
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
         DGSS=0.
         EDS1S=0.
         PPBLS=0.
         EVAPS=0.
C         EKMS=0.
C         RNLDS=0.
C         FRCTVS=0.
C         PSIS=0.
         DBLS=0.
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
      TG1=ODATA(I,J,1)
      SRHEAT=SRCOR*FSF(I,J,ITYPE)*COSZ1(I,J)
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
      NGRNDZ=NGRND
      SNOW=GDATA(I,J,1)
      TG1=TGRND(I,J,2)
      TG2=TGRN2(I,J,2)
      ACE2=ODATA(I,J,3)
      SRHEAT=SRCOR*FSF(I,J,ITYPE)*COSZ1(I,J)
            OA(I,J,12)=OA(I,J,12)+SRHEAT*DTSURF
      Z2=ACE2/RHOI
      Z2BY4L=Z2/(4.*ALAMI)
      Z1BY6L=(Z1IBYL+SNOW*BYRLS)*.1666667
      CDTERM=1.5*TG2-.5*TFO
      CDENOM=1./(2.*Z1BY6L+Z2BY4L)
      HC1=HC1I+SNOW*SHI
      BETA=1.
      ELHX=LHS
C*
      MSI1 = SNOW+ACE1I ! snow and first layer ice mass (kg/m^2)
      MSI2 = ACE2 ! second (physical) layer ice mass (kg/m^2)
      dF1dTG = 2./(ACE1I*BYRLI+SNOW*BYRLS)
      HCG1 = SHI*XSI1*MSI1 ! heat capacity of first layer ice (J/C*m^2)
      HCG2 = SHI*XSI2*MSI1 ! heat capacity of second layer ice
      GO TO 3000
C****
 2400 IF (PLICE.LE.0.) then
        ipbl(i,j,3)=0.
        GO TO 5000
      endif
      NGRNDZ=NGRND
C****
C**** LAND ICE
C****
      ITYPE=3
      PTYPE=PLICE
      SNOW=GDATA(I,J,12)
      TG1=TGRND(I,J,3)
      TG2=GDATA(I,J,14)
      SRHEAT=FSF(I,J,ITYPE)*COSZ1(I,J)
      Z1BY6L=(Z1IBYL+SNOW*BYRLS)*.1666667
      CDTERM=TG2
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)
      HCG1=HC1I+SNOW*SHI
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
      QG=QSAT(TG,PS,ELHX)
      TGV=TG*(1.+QG*RVX)
C =====================================================================
      CALL PBL(I,J,ITYPE)
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
C*
      BETAUP = BETA
      IF (QS .GT. QG) BETAUP = 1.
      EVHEAT=(LHE+TG1*SHV)*BETAUP*RCDQWS*(QS-QG)
      TRHEAT=TRHR(I,J,1)-STBO*(TG*TG)*(TG*TG)
      IF(ITYPE.EQ.1) GO TO 3620
C**** CALCULATE FLUXES USING IMPLICIT TIME STEP FOR NON-OCEAN POINTS
C*
      IF (ITYPE .EQ. 2) GO TO 3550
C*
      F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
      F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
      DSHDTG=-RCDHWS*KH*SHA/(DHGS+KH)
      DQGDTG=QG*ELHX/(RVAP*TG*TG)
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
C*
      F1 = (TG1-TG2)*dF1dTG ! heat flux on first/second layers (W/m^2)
C*
      EVHEAT = LHE*RCDQWS*(QS-QG) ! latent heat flux (W/m^2)
C*
      F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
C*
      dSNdTG=-RCDHWS*KH*SHA/(DHGS+KH)
C*
      dQGdTG = QG*ELHX/(RVAP*TG*TG) ! d(QG)/dTG
      dEVdTG = -dQGdTG*LHE*RCDQWS*KH/(DGS+KH) ! d(EVHEAH)/dTG
C*
      dTRdTG = -4*STBO*TG*TG*TG ! d(TRHEAT)/dTG
      dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG
C     dSNdHS = RCDHWS ! d(SHEAT)/dHS - kg/(sec*m^2)
      dEVdQS = LHE*RCDQWS ! d(EVHEAH)/dQS
C*
      HSDEN = -(1.+2.*S1BYG1)*DTGRND*PGK*BETA*dSNdTG+MA1*PKDN*SHCD
      HSCON = -(1.+2.*S1BYG1)*DTGRND*PGK*SHEAT/HSDEN ! (J*sec)/kg
      HSMUL = -(1.+2.*S1BYG1)*DTGRND*PGK*BETA*dSNdTG/HSDEN ! J/(kg*degC)
C*
      QSDEN = (1.+2.*S1BYG1)*BETA*DTGRND*dEVdQS+MA1*LHE
      QSCON = -(1.+2.*S1BYG1)*DTGRND*EVHEAT/QSDEN
      QSMUL = -(1.+2.*S1BYG1)*DTGRND*BETA*dEVdTG/QSDEN
C*
      T2DEN = HCG2+BETA*DTGRND*dF1dTG
      T2CON = DTGRND*F1/T2DEN
      T2MUL = BETA*DTGRND*dF1dTG/T2DEN
C*
      TGDEN = HCG1-BETA*DTGRND*(dF0dTG-dF1dTG-
     A        HSMUL*dSNdTG+QSMUL*dEVdQS+T2MUL*dF1dTG) ! W/(m^2*degC)
      dTG = DTGRND*(F0-F1+BETA*
     A      (QSCON*dEVdQS-HSCON*dSNdTG+T2CON*dF1dTG))/TGDEN ! degC
C*
      IF (TG1+dTG .GT. 0.) dTG = -TG1
C*
      dHS = HSCON+HSMUL*dTG ! (J*sec)/kg
      dQS = QSCON+QSMUL*dTG
      dT2 = T2CON+T2MUL*dTG
C*
      SHDT = DTGRND*(SHEAT+BETA*((dTG-dHS)*dSNdTG)) ! sensible
C*
      EVHDT = DTGRND*(EVHEAT+BETA*(dTG*dEVdTG+dQS*dEVdQS)) ! latent
      TRHDT = DTGRND*(TRHEAT+BETA*dTG*dTRdTG) ! thermal flux (J/m^2)
      F1DT = DTGRND*(F1+BETA*(dTG*dF1dTG-dT2*dF1dTG))
      DU1(I,J)=DU1(I,J)+PTYPE*DTGRND*RCDMWS*US/RMBYA
      DV1(I,J)=DV1(I,J)+PTYPE*DTGRND*RCDMWS*VS/RMBYA
C*
      TG1 = TG1+dTG ! first layer sea ice temperature (degC)
      TG2 = TG2+dT2 ! second layer sea ice temperature (degC)
      TGRN2(I,J,ITYPE) = TG2
 3600 CONTINUE
      GO TO 3700
C**** CALCULATE FLUXES USING IMPLICIT TIME STEP ALSO FOR OCEAN POINTS
 3620 CONTINUE
C*
      DSHDTG=-RCDHWS*SHA
      dEVdQS = LHE*RCDQWS
      dHS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A      ((1.+2.*S1BYG1)*DTSURF*PGK*RCDHWS+MA1*PKDN)
C*
      dTS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A      (MA1*PKDN*SHCD-(1.+2.*S1BYG1)*DTSURF*PGK*DSHDTG)
C*
      dQS = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/
     A      ((1.+2.*S1BYG1)*DTSURF*dEVdQS+MA1*LHE)
C*
      SHDT = DTSURF*(SHEAT-dTS*DSHDTG)
      EVHDT=DTSURF*(EVHEAT+dQS*dEVdQS) ! latent heat flux
      TRHDT=DTSURF*TRHEAT
C*
      DU1(I,J)=DU1(I,J)+PTYPE*DTGRND*RCDMWS*US/RMBYA
      DV1(I,J)=DV1(I,J)+PTYPE*DTGRND*RCDMWS*VS/RMBYA
C
C**** CALCULATE EVAPORATION
 3700 CONTINUE
      DQ1X =EVHDT/((LHE+TG1*SHV)*RMBYA)
      EVHDT0=EVHDT
C*
      IF (DQ1X .LE. Q1+DQ1(I,J)) GO TO 3720
      DQ1X = Q1+DQ1(I,J)
      EVHDT=DQ1X*(LHE+TG1*SHV)*RMBYA
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*RMBYA
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP
      TGRND(I,J,ITYPE)=TG1
      DTH1(I,J)=DTH1(I,J)-SHDT*PTYPE/(SHA*RMBYA*P1K)
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
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
         DGSS=DGSS+DGS*PTYPE
         EDS1S=EDS1S+KH*PTYPE
         PPBLS=PPBLS+PPBL*PTYPE
         EVAPS=EVAPS+EVAP*PTYPE
C****
C         EKMS  =EKMS  +EKMAN*PTYPE
C         RNLDS =RNLDS +REYNLD*PTYPE
C         FRCTVS=FRCTVS+USTAR*PTYPE
C         PSIS=PSIS+PSI*DEGREE*PTYPE
         DBLS  =DBLS  +DBL*PTYPE
5666  GO TO (4000,4100,4400),ITYPE
C****
C**** OCEAN
C****
 4000    ASHDT=ASHDT+SHDT*POCEAN
         AEVHDT=AEVHDT+EVHDT*POCEAN
         ATRHDT=ATRHDT+TRHDT*POCEAN
         ATS=ATS+(TS-TF)*POCEAN
            OA(I,J,6)=OA(I,J,6)+TRHDT
            OA(I,J,7)=OA(I,J,7)+SHDT
            OA(I,J,8)=OA(I,J,8)+EVHDT
      GO TO 2200
C****
C**** OCEAN ICE
C****
 4100    CSHDT=CSHDT+SHDT*POICE
         CEVHDT=CEVHDT+EVHDT*POICE
         CTRHDT=CTRHDT+TRHDT*POICE
         CTS=CTS+(TS-TF)*POICE
         IF (TG1.GT.TDIURN(I,J,7)) TDIURN(I,J,7) = TG1
            OA(I,J,9)=OA(I,J,9)+TRHDT
            OA(I,J,10)=OA(I,J,10)+SHDT
            OA(I,J,11)=OA(I,J,11)+EVHDT
      GO TO 2400
C****
C**** LAND ICE
C****
 4400    BSHDT=BSHDT+SHDT*PLICE
         BEVHDT=BEVHDT+EVHDT*PLICE
         BTRHDT=BTRHDT+TRHDT*PLICE
         BTS=BTS+(TS-TF)*PLICE
         IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
         AIJ(I,J,71)=AIJ(I,J,71)+(TS-TF)
         AIJ(I,J,73)=AIJ(I,J,73)+SHDT
         AIJ(I,J,74)=AIJ(I,J,74)+EVHDT
         AIJ(I,J,75)=AIJ(I,J,75)+TRHDT
C**** NON-OCEAN POINTS WHICH ARE NOT MELTING OR FREEZING WATER USE
C****   IMPLICIT TIME STEPS
C****
C**** UPDATE SURFACE AND FIRST LAYER QUANTITIES
C****
5000  BLTEMP(I,J,1)=WSS
      BLTEMP(I,J,2)=TSS
      BLTEMP(I,J,3)=QSS
      BLTEMP(I,J,6)=USS
      BLTEMP(I,J,7)=VSS
      BLTEMP(I,J,8)=TAUS
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
         IF(JR.EQ.24) GO TO 5700
         DJ(JR,9)=DJ(JR,9)+TRHDTS*DXYPJ
         DJ(JR,13)=DJ(JR,13)+SHDTS*DXYPJ
         DJ(JR,14)=DJ(JR,14)+EVHDTS*DXYPJ
         DJ(JR,19)=DJ(JR,19)+EVAPS*DXYPJ
         IF(MODDSF.NE.0) GO TO 5700
         DJ(JR,23)=DJ(JR,23)+(TSS-TFS)*DXYPJ
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
 5700    AIJ(I,J,4)=AIJ(I,J,4)+SHDTS
         IF(MODRD.EQ.0) AIJ(I,J,21)=AIJ(I,J,21)+TRHDTS/DTSRCE
         AIJ(I,J,22)=AIJ(I,J,22)+(SRHDTS+TRHDTS)
         AIJ(I,J,23)=AIJ(I,J,23)+(SRHDTS+TRHDTS+SHDTS+EVHDTS)
         IF(MODDSF.NE.0) GO TO 5800
         AIJ(I,J,34)=AIJ(I,J,34)+WSS              ! added 12/29/96 -rar-
         AIJ(I,J,35)=AIJ(I,J,35)+(TSS-TFS)
         AIJ(I,J,36)=AIJ(I,J,36)+USS
         AIJ(I,J,37)=AIJ(I,J,37)+VSS
         AIJ(I,J,47)=AIJ(I,J,47)+RTAUS
         AIJ(I,J,48)=AIJ(I,J,48)+RTAUUS
         AIJ(I,J,49)=AIJ(I,J,49)+RTAUVS
         AIJ(I,J,51)=AIJ(I,J,51)+QSS
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAG6
 5800    IF(MODD6.NE.0) GO TO 6000
         DO 5820 KR=1,4
         IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) GO TO 5840
 5820    CONTINUE
         GO TO 6000
 5840    ADAILY(IHOUR,1,KR)=ADAILY(IHOUR,1,KR)+S0*COSZ1(I,J)
         ADAILY(IHOUR,6,KR)=ADAILY(IHOUR,6,KR)+PS
         ADAILY(IHOUR,7,KR)=ADAILY(IHOUR,7,KR)+PSK*T(I,J,5)
         ADAILY(IHOUR,8,KR)=ADAILY(IHOUR,8,KR)+PSK*T(I,J,4)
         ADAILY(IHOUR,9,KR)=ADAILY(IHOUR,9,KR)+PSK*T(I,J,3)
         ADAILY(IHOUR,10,KR)=ADAILY(IHOUR,10,KR)+PSK*T(I,J,2)
         ADAILY(IHOUR,11,KR)=ADAILY(IHOUR,11,KR)+PSK*T(I,J,1)
         ADAILY(IHOUR,12,KR)=ADAILY(IHOUR,12,KR)+TSS
         ADAILY(IHOUR,13,KR)=ADAILY(IHOUR,13,KR)+(TG1S+TFS)
         ADAILY(IHOUR,14,KR)=ADAILY(IHOUR,14,KR)+Q(I,J,5)
         ADAILY(IHOUR,15,KR)=ADAILY(IHOUR,15,KR)+Q(I,J,4)
         ADAILY(IHOUR,16,KR)=ADAILY(IHOUR,16,KR)+Q(I,J,3)
         ADAILY(IHOUR,17,KR)=ADAILY(IHOUR,17,KR)+Q(I,J,2)
         ADAILY(IHOUR,18,KR)=ADAILY(IHOUR,18,KR)+Q1
         ADAILY(IHOUR,19,KR)=ADAILY(IHOUR,19,KR)+QSS
         ADAILY(IHOUR,20,KR)=ADAILY(IHOUR,20,KR)+QGS
         ADAILY(IHOUR,28,KR)=ADAILY(IHOUR,28,KR)+SRHDTS
         ADAILY(IHOUR,29,KR)=ADAILY(IHOUR,29,KR)+TRHDTS
         ADAILY(IHOUR,30,KR)=ADAILY(IHOUR,30,KR)+SHDTS
         ADAILY(IHOUR,31,KR)=ADAILY(IHOUR,31,KR)+EVHDTS
         ADAILY(IHOUR,32,KR)=ADAILY(IHOUR,32,KR)
     *                        +(SRHDTS+TRHDTS+SHDTS+EVHDTS)
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
         AJ(J,9)=AJ(J,9)+ATRHDT
         BJ(J,9)=BJ(J,9)+BTRHDT
         CJ(J,9)=CJ(J,9)+CTRHDT
         AJ(J,13)=AJ(J,13)+ASHDT
         BJ(J,13)=BJ(J,13)+BSHDT
         CJ(J,13)=CJ(J,13)+CSHDT
         AJ(J,14)=AJ(J,14)+AEVHDT
         BJ(J,14)=BJ(J,14)+BEVHDT
         CJ(J,14)=CJ(J,14)+CEVHDT
         IF(MODDSF.NE.0) GO TO 7000
         AJ(J,23)=AJ(J,23)+ATS
         BJ(J,23)=BJ(J,23)+BTS
         CJ(J,23)=CJ(J,23)+CTS
 7000 CONTINUE
C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODD6)
C****
C**** UPDATE FIRST LAYER QUANTITIES
C****
      DO 7205 J=1,JM
        IMAX=IMAXJ(J)
        DO 7200 I=1,IMAX
          FTEVAP=0
          IF (DTH1(I,J)*T(I,J,1).lt.0) FTEVAP=-DTH1(I,J)/T(I,J,1)
          FQEVAP=0
          IF (DQ1(I,J).lt.0.and.Q(I,J,1).gt.0) FQEVAP=-DQ1(I,J)/Q(I,J,1)
            T(I,J,1)=  T(I,J,1)+DTH1(I,J)
           TX(I,J,1)= TX(I,J,1)*(1.-FTEVAP)
           TY(I,J,1)= TY(I,J,1)*(1.-FTEVAP)
          TXX(I,J,1)=TXX(I,J,1)*(1.-FTEVAP)
          TXY(I,J,1)=TXY(I,J,1)*(1.-FTEVAP)
          TYY(I,J,1)=TYY(I,J,1)*(1.-FTEVAP)
          TYZ(I,J,1)=TYZ(I,J,1)*(1.-FTEVAP)
          TZX(I,J,1)=TZX(I,J,1)*(1.-FTEVAP)
           TZ(I,J,1)= TZ(I,J,1)*(1.-FTEVAP) ! should be set from PBL
          TZZ(I,J,1)=TZZ(I,J,1)*(1.-FTEVAP) ! should be set from PBL
            Q(I,J,1)=  Q(I,J,1)+DQ1(I,J)
           QX(I,J,1)= QX(I,J,1)*(1.-FQEVAP)
           QY(I,J,1)= QY(I,J,1)*(1.-FQEVAP)
          QXX(I,J,1)=QXX(I,J,1)*(1.-FQEVAP)
          QXY(I,J,1)=QXY(I,J,1)*(1.-FQEVAP)
          QYY(I,J,1)=QYY(I,J,1)*(1.-FQEVAP)
          QYZ(I,J,1)=QYZ(I,J,1)*(1.-FQEVAP)
          QZX(I,J,1)=QZX(I,J,1)*(1.-FQEVAP)
           QZ(I,J,1)= QZ(I,J,1)*(1.-FQEVAP) ! should be set from PBL
          QZZ(I,J,1)=QZZ(I,J,1)*(1.-FQEVAP) ! should be set from PBL
          IF (Q(I,J,1).LT.qmin) THEN
             WRITE (99,*) TAU,I,J,' Q1:',Q(I,J,1),'-> 0',DQ1(I,J)
             Q(I,J,1)=qmin
             QX(I,J,1)=0.
             QY(I,J,1)=0.
             QZ(I,J,1)=0.
             QXX(I,J,1)=0.
             QYY(I,J,1)=0.
             QZZ(I,J,1)=0.
             QXY(I,J,1)=0.
             QYZ(I,J,1)=0.
             QZX(I,J,1)=0.
          ENDIF
7200    CONTINUE
7205  CONTINUE
C****
C**** ADD IN SURFACE FRICTION TO FIRST LAYER WIND
C****
      DO 7600 I=1,IM
      U(I,2,1)=U(I,2,1)-2.*(DU1(1,1)*COSI(I)-DV1(1,1)*SINI(I))*RAPVN(1)
      V(I,2,1)=V(I,2,1)-2.*(DV1(1,1)*COSI(I)+DU1(1,1)*SINI(I))*RAPVN(1)
      U(I,JM,1)=U(I,JM,1)
     *  -2.*(DU1(1,JM)*COSI(I)+DV1(1,JM)*SINI(I))*RAPVS(JM)
 7600 V(I,JM,1)=V(I,JM,1)
     *  -2.*(DV1(1,JM)*COSI(I)-DU1(1,JM)*SINI(I))*RAPVS(JM)
      DO 7700 J=2,JMM1
      I=IM
      DO 7700 IP1=1,IM
      U(I,J,1)=U(I,J,1)-(DU1(I,J)+DU1(IP1,J))*RAPVS(J)
      V(I,J,1)=V(I,J,1)-(DV1(I,J)+DV1(IP1,J))*RAPVS(J)
      U(I,J+1,1)=U(I,J+1,1)-(DU1(I,J)+DU1(IP1,J))*RAPVN(J)
      V(I,J+1,1)=V(I,J+1,1)-(DV1(I,J)+DV1(IP1,J))*RAPVN(J)
 7700 I=IP1
C****
C**** DRY CONVECTION ORIGINATING FROM THE FIRST LAYER
C****
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.
      DO 8050 L=1,LM
      DO 8050 J=2,JM
      DO 8050 I=1,IM
      UT(I,J,L)=U(I,J,L)
 8050 VT(I,J,L)=V(I,J,L)
C**** OUTSIDE LOOPS OVER J AND I
      DO 8500 J=1,JM
      POLE=.FALSE.
      IF(J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.
      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)

      IM1=IM
      DO 8500 I=1,IMAX
         DO K=1,KMAX
            RA(K)=RAJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
      BLDATA(I,J,4)=1.
      IF(T(I,J,1)*(1.+Q(I,J,1)*RVX).LE.
     *   T(I,J,2)*(1.+Q(I,J,2)*RVX)) GO TO 8500
C**** MIX HEAT AND MOISTURE THROUGHOUT THE BOUNDARY LAYER
      PIJ=P(I,J)
      PKMS=(PK(1,I,J)*DSIG(1)+PK(2,I,J)*DSIG(2))*PIJ
      THPKMS=(T(I,J,1)*(PK(1,I,J)*DSIG(1))+T(I,J,2)*(PK(2,I,J)*DSIG(2)))
     *   *PIJ
      TXS= (TX(I,J,1)*(PK(1,I,J)*DSIG(1)) + TX(I,J,2)*(PK(2,I,J)*
     *     DSIG(2)))*PIJ
      TYS= (TY(I,J,1)*(PK(1,I,J)*DSIG(1)) + TY(I,J,2)*(PK(2,I,J)*
     *     DSIG(2)))*PIJ
      TXXS=(TXX(I,J,1)*(PK(1,I,J)*DSIG(1))+TXX(I,J,2)*(PK(2,I,J)*
     *     DSIG(2)))*PIJ
      TYYS=(TYY(I,J,1)*(PK(1,I,J)*DSIG(1))+TYY(I,J,2)*(PK(2,I,J)*
     *     DSIG(2)))*PIJ
      TXYS=(TXY(I,J,1)*(PK(1,I,J)*DSIG(1))+TXY(I,J,2)*(PK(2,I,J)*
     *     DSIG(2)))*PIJ
      QMS=(Q(I,J,1)*DSIG(1)+Q(I,J,2)*DSIG(2))*PIJ
      QXS  = (QX(I,J,1)*DSIG(1) +  QX(I,J,2)*DSIG(2))*PIJ
      QYS  = (QY(I,J,1)*DSIG(1) +  QY(I,J,2)*DSIG(2))*PIJ
      QXXS =(QXX(I,J,1)*DSIG(1) + QXX(I,J,2)*DSIG(2))*PIJ
      QYYS =(QYY(I,J,1)*DSIG(1) + QYY(I,J,2)*DSIG(2))*PIJ
      QXYS =(QXY(I,J,1)*DSIG(1) + QXY(I,J,2)*DSIG(2))*PIJ
      TVMS=(T(I,J,1)*(1.+Q(I,J,1)*RVX)*(PK(1,I,J)*DSIG(1))
     *    +T(I,J,2)*(1.+Q(I,J,2)*RVX)*(PK(2,I,J)*DSIG(2)))*PIJ
      THETA=TVMS/PKMS
C**** MIX THROUGH SUMSEQUENT LAYERS
      DO 8140 L=3,LM
      IF(THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 8160
      IF(L.EQ.LS1) PIJ=PSF-PTOP
      PKMS=PKMS+(PK(L,I,J)*(DSIG(L)*PIJ))
      THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*(DSIG(L)*PIJ))
       TXS =  TXS +  TX(I,J,L)*(PK(L,I,J)*DSIG(L))*PIJ
       TYS =  TYS +  TY(I,J,L)*(PK(L,I,J)*DSIG(L))*PIJ
      TXXS = TXXS + TXX(I,J,L)*(PK(L,I,J)*DSIG(L))*PIJ
      TYYS = TYYS + TYY(I,J,L)*(PK(L,I,J)*DSIG(L))*PIJ
      TXYS = TXYS + TXY(I,J,L)*(PK(L,I,J)*DSIG(L))*PIJ
      QMS=QMS+Q(I,J,L)*(DSIG(L)*PIJ)
       QXS =  QXS +  QX(I,J,L)*(DSIG(L)*PIJ)
       QYS =  QYS +  QY(I,J,L)*(DSIG(L)*PIJ)
      QXXS = QXXS + QXX(I,J,L)*(DSIG(L)*PIJ)
      QYYS = QYYS + QYY(I,J,L)*(DSIG(L)*PIJ)
      QXYS = QXYS + QXY(I,J,L)*(DSIG(L)*PIJ)
      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*(DSIG(L)*PIJ))
 8140 THETA=TVMS/PKMS
      L=LM+1
 8160 LMAX=L-1
      RDP=1./(P(I,J)*SIGE(1)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      BLDATA(I,J,4)=LMAX
      PIJ=P(I,J)
      DO 8180 L=1,LMAX
      IF(L.EQ.LS1) PIJ=(PSF-PTOP)
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(L,I,J)*PIJ
      T(I,J,L)=THM
       TX(I,J,L) = TXS/PKMS
       TY(I,J,L) = TYS/PKMS
       TZ(I,J,L) = 0.
      TXX(I,J,L) = TXXS/PKMS
      TYY(I,J,L) = TYYS/PKMS
      TXY(I,J,L) = TXYS/PKMS
      TZZ(I,J,L) = 0.
      TZX(I,J,L) = 0.
      TYZ(I,J,L) = 0.
      Q(I,J,L)=QMS
       QX(I,J,L) = QXS*RDP
       QY(I,J,L) = QYS*RDP
       QZ(I,J,L) = 0.
      QXX(I,J,L) = QXXS*RDP
      QYY(I,J,L) = QYYS*RDP
      QXY(I,J,L) = QXYS*RDP
      QZZ(I,J,L) = 0.
      QZX(I,J,L) = 0.
      QYZ(I,J,L) = 0.
 8180 CONTINUE
C**** MIX MOMENTUM THROUGHOUT THE BOUNDARY LAYER
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      PIJ=P(I,J)
      DO L=1,LMAX
         IF(L.EQ.LS1) PIJ=PSF-PTOP
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*PIJ*DSIG(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*PIJ*DSIG(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      PIJ=P(I,J)
      DO L=1,LMAX
         IF(L.EQ.LS1) PIJ=PSF-PTOP
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c the following line gives bytewise different ajl
            AJL(IDJ(K),L,38)=AJL(IDJ(K),L,38)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PIJ*RA(K)
         ENDDO
      ENDDO
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
 8400    IF(MODD6.NE.0) GO TO 8500
         DO KR=1,4
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(IHOUR,47,KR)=ADAILY(IHOUR,47,KR)+1.
               ADAILY(IHOUR,48,KR)=ADAILY(IHOUR,48,KR)+LMAX
            END IF
         END DO
 8500 IM1=I
C****
 9000 CONTINUE
      RETURN
 9991 FORMAT ('0SURFACE ',4I4,5F10.4,3F11.7)
 9992 FORMAT ('0',I2,10F10.4/23X,4F10.4,10X,2F10.4/
     *  33X,3F10.4,10X,2F10.4)
 9993 FORMAT ('0',I2,10F10.4/23X,7F10.4/33X,7F10.4)
 9994 FORMAT ('0',I2,11F10.4)
      END
