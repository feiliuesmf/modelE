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
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,rvap,stbo,bygrav,by6
      USE E001M12_COM, only : im,jm,lm,fim,DTsrc,NIsurf,u,v,t,p,q
     *     ,idacc,dsig,jday,gdata,ndasf,jeq,fland,flice
     *     ,fearth,nday,modrd,ijd6,ITime,JHOUR,sige,byim
      USE SOMTQ_COM, only : tx,ty,tz,txx,tyy,tzz,txy,tzx,tyz,qx,qy,qz
     *     ,qxx,qyy,qzz,qxy,qzx,qyz
      USE GEOM, only : dxyp,imaxj,kmaxj,raj,idij,idjj,rapvn,rapvs,sini
     *     ,cosi
      USE RADNCB, only : trhr,fsf,cosz1
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs
     &     ,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
      USE SOCPBL, only : zgs
      USE DAGCOM, only : aij,tdiurn,aj,bj,cj,areg,ajl,adaily,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,j_tsrf,j_evap,j_evhdt,j_shdt,j_trhdt
      USE DYNAMICS, only : pmid,pk,pedn,pek,pdsig,plij
      USE LANDICE, only : hc2li,z1e,z2li,hc1li
      USE OCEAN, only : odata,oa,tfo
      USE SEAICE, only : xsi1,xsi2,z1i,ace1i,hc1i,alami,byrli,byrls,rhos

      IMPLICIT NONE

      INTEGER I,J,L,K,IM1,IP1,LMAX,KR,JR,NS,NSTEPS,MODDSF,MODD6
     *     ,KMAX,IMAX,ITYPE,NGRNDZ,NG
      REAL*8 ATRHDT,BTRHDT
     *     ,CTRHDT,ASHDT,BSHDT,CSHDT,AEVHDT,BEVHDT,CEVHDT,ATS,BTS,CTS
     *     ,PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1,P1K,H0M1,HZM1,PGK,HS,PKDN
     *     ,DXYPJ,BETAS,EVHDTS,CDMS,CDHS,DGSS,EDS1S,PPBLS,EVAPS,DBLS
     *     ,BETA,ELHX,ACE2,CDTERM,CDENOM,HC1,dF1dTG,HCG1,HCG2,DTGRND
     *     ,EVHDT,F1DT,CM,CH,CQ,DHGS,DQGS,DGS,BETAUP,EVHEAT,F0
     *     ,F1,DSHDTG,DQGDTG,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG
     *     ,dEVdQS,HSDEN,HSCON,HSMUL,dHS,dQS,dT2,dTS,DQ1X,EVHDT0,EVAP
     *     ,F0DT,FTEVAP,VAP,PKMS,SPRING,TIMEZ,PWATER
     *     ,PXSOIL,PSK,TH1,Q1,THV1,TFS,RMBYA,HZM1,Q0M1,QZM1,TSS,QSS,TAUS
     *     ,RTAUS,RTAUUS,RTAUVS,TG1S,QGS,SRHDTS,TRHDTS,SHDTS,UGS,PTYPE
     *     ,TG1,SRHEAT,SNOW,TG2,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS
     *     ,RCDHWS,RCDQWS,SHEAT,TRHEAT,QSDEN,QSCON,QSMUL,T2DEN,T2CON
     *     ,T2MUL,TGDEN,FQEVAP,THPKMS,TXS,TYS,TXXS,TYYS,TXYS,QMS,QXS,QYS
     *     ,QXXS,QYYS,QXYS,TVMS,THETA,RDP,THM,ZS1CO,WARMER
     *     ,USS,VSS,WSS,VGS,WGS,USRS,VSRS,Z2,Z2BY4L,Z1BY6L
     *     ,THZ1,QZ1

      REAL*8 MSUM, MA1, MSI1, MSI2
      REAL*8, DIMENSION(IM,JM) :: DTH1,DQ1,DU1,DV1
      COMMON /WORK1d/DTH1,DQ1
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT
      COMMON/WORK2/UT,VT,DU1,DV1

      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      REAL*8, DIMENSION(IM,JM,4) :: E0,E1,EVAPOR,TGRND,TGRN2
      COMMON/WORK3/E0,E1,EVAPOR,TGRND

C**** Interface to PBL
      REAL*8 ZS1,TGV,TKV,QG,HEMI,DTSURF,US,VS,WS,TSV,QS,PSI,DBL,KM,KH,
     *     PPBL,UG,VG,WG,ZMIX
      LOGICAL POLE
      COMMON /PBLPAR/ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE

      COMMON /PBLOUT/US,VS,WS,TSV,QS,PSI,DBL,KM,KH,PPBL,
     2               UG,VG,WG,ZMIX

      REAL*8, PARAMETER :: qmin=1.e-12

      REAL*8, PARAMETER :: S1BYG1 = 0.57735, RVX=0.,
     *     Z1IBYL=Z1I/ALAMI, Z2LI3L=Z2LI/(3.*ALAMI), Z1LIBYL=Z1E/ALAMI

      REAL*8 QSAT,DQSATDT
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
      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf

      ZS1CO=.5*DSIG(1)*RGAS*BYGRAV

      SPRING=-1.
      IF((JDAY.GE.32).AND.(JDAY.LE.212)) SPRING=1.
C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND
      DO J=1,JM
        DO I=1,IM
          TGRND(I,J,2)=GDATA(I,J,3)
          TGRND(I,J,3)=GDATA(I,J,13)
C         TGRND(I,J,4)=GDATA(I,J,4)
          TGRN2(I,J,2) = GDATA(I,J,7)
          TGRN2(I,J,3) = GDATA(I,J,14)
        END DO
      END DO
C*
C**** Zero out fluxes summed over type
      E0=0. ; E1=0. ; EVAPOR=0.

      call pgrads1
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
C**** ZERO OUT SURFACE DIAGNOSTICS WHICH WILL BE SUMMED OVER LONGITUDE
         ATRHDT=0.
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
         IF(J.LT.JEQ) WARMER=-SPRING
         IF(J.GE.JEQ) WARMER=SPRING
      IM1=IM
      DO I=1,IMAX

      ! until pbl loops over i,j,itype
      WSAVG(I,J)=0.
      TSAVG(I,J)=0.
      QSAVG(I,J)=0.
      USAVG(I,J)=0.
      VSAVG(I,J)=0.
      TAUAVG(I,J)=0.

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
      SRHEAT=FSF(I,J,ITYPE)*COSZ1(I,J)
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
      SNOW=GDATA(I,J,1)
      TG1=TGRND(I,J,2)
      TG2=TGRN2(I,J,2)
      ACE2=ODATA(I,J,3)
      SRHEAT=FSF(I,J,ITYPE)*COSZ1(I,J)
            OA(I,J,12)=OA(I,J,12)+SRHEAT*DTSURF
      Z2=ACE2/RHOI
      Z2BY4L=Z2/(4.*ALAMI)
      Z1BY6L=(Z1IBYL+SNOW*BYRLS)*BY6
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
      NGRNDZ=1    ! NIgrnd>1 currently not an option
C****
C**** LAND ICE
C****
      ITYPE=3
      PTYPE=PLICE
      SNOW=GDATA(I,J,12)
      TG1=TGRND(I,J,3)
      TG2=GDATA(I,J,14)
      SRHEAT=FSF(I,J,ITYPE)*COSZ1(I,J)
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
      TRHEAT=TRHR(I,J,1)-STBO*(TG*TG)*(TG*TG)
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
      F1 = (TG1-TG2)*dF1dTG ! heat flux on first/second layers (W/m^2)
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
      T2CON = DTGRND*F1/T2DEN
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
      TGRN2(I,J,ITYPE) = TG2
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
         AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(TS-TF)
         AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)+SHDT
         AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)+EVHDT
         AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)+TRHDT
C**** NON-OCEAN POINTS WHICH ARE NOT MELTING OR FREEZING WATER USE
C****   IMPLICIT TIME STEPS
C****
C**** UPDATE SURFACE AND FIRST LAYER QUANTITIES
C****
5000  CONTINUE
C****
C**** ACCUMULATE DIAGNOSTICS
C****
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
         IF(JR.EQ.24) GO TO 5700
         AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+TRHDTS*DXYPJ
         AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+SHDTS*DXYPJ
         AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+EVHDTS*DXYPJ
         AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+EVAPS*DXYPJ
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
         DO KR=1,4
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
               ADAILY(JHOUR+1,42,KR)=ADAILY(JHOUR+1,42,KR)+CDMS
               ADAILY(JHOUR+1,43,KR)=ADAILY(JHOUR+1,43,KR)+CDHS
               ADAILY(JHOUR+1,44,KR)=ADAILY(JHOUR+1,44,KR)+DGSS
               ADAILY(JHOUR+1,45,KR)=ADAILY(JHOUR+1,45,KR)+EDS1S
               ADAILY(JHOUR+1,46,KR)=ADAILY(JHOUR+1,46,KR)+DBLS
               ADAILY(JHOUR+1,50,KR)=ADAILY(JHOUR+1,50,KR)+EVAPS
            END IF
         END DO
       END IF
       IM1=I
      END DO
C**** QUANTITIES ACCUMULATED FOR SURFACE TYPE TABLES IN DIAGJ
         AJ(J,J_TRHDT)=AJ(J,J_TRHDT)+ATRHDT
         BJ(J,J_TRHDT)=BJ(J,J_TRHDT)+BTRHDT
         CJ(J,J_TRHDT)=CJ(J,J_TRHDT)+CTRHDT
         AJ(J,J_SHDT )=AJ(J,J_SHDT )+ASHDT
         BJ(J,J_SHDT )=BJ(J,J_SHDT )+BSHDT
         CJ(J,J_SHDT )=CJ(J,J_SHDT )+CSHDT
         AJ(J,J_EVHDT)=AJ(J,J_EVHDT)+AEVHDT
         BJ(J,J_EVHDT)=BJ(J,J_EVHDT)+BEVHDT
         CJ(J,J_EVHDT)=CJ(J,J_EVHDT)+CEVHDT
         IF(MODDSF.NE.0) GO TO 7000
         AJ(J,J_TSRF)=AJ(J,J_TSRF)+ATS
         BJ(J,J_TSRF)=BJ(J,J_TSRF)+BTS
         CJ(J,J_TSRF)=CJ(J,J_TSRF)+CTS
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
             WRITE(99,*) ITime,'I,J:',I,J,' Q1:',Q(I,J,1),'->0',DQ1(I,J)
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
            U(IDIJ(K,I,J),IDJJ(K,J),1)=U(IDIJ(K,I,J),IDJJ(K,J),1)-RAJ(K
     *           ,J)*(DU1(I,J)*COSI(K)+DV1(I,J)*SINI(K)*HEMI)
            V(IDIJ(K,I,J),IDJJ(K,J),1)=V(IDIJ(K,I,J),IDJJ(K,J),1)-RAJ(K
     *           ,J)*(DV1(I,J)*COSI(K)-DU1(I,J)*SINI(K)*HEMI)
          END DO
        END DO
      END DO
C**** non polar boxes
      DO J=2,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
          DO K=1,KMAX
            U(IDIJ(K,I,J),IDJJ(K,J),1)=U(IDIJ(K,I,J),IDJJ(K,J),1)
     *           - RAJ(K,J)*DU1(I,J)
            V(IDIJ(K,I,J),IDJJ(K,J),1)=V(IDIJ(K,I,J),IDJJ(K,J),1)
     *           - RAJ(K,J)*DV1(I,J)
          END DO
        END DO
      END DO
C****
C**** DRY CONVECTION ORIGINATING FROM THE FIRST LAYER
C****
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.
      UT=U ; VT=V
C**** OUTSIDE LOOPS OVER J AND I
      DO J=1,JM
      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)

      DO I=1,IMAX
        DO K=1,KMAX
          RA(K)=RAJ(K,J)
          IDI(K)=IDIJ(K,I,J)
          IDJ(K)=IDJJ(K,J)
        END DO
      DCLEV(I,J)=1.
      IF(T(I,J,1)*(1.+Q(I,J,1)*RVX).GT.
     *   T(I,J,2)*(1.+Q(I,J,2)*RVX)) THEN
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
     *     +T(I,J,2)*(1.+Q(I,J,2)*RVX)*(PK(2,I,J)*DSIG(2)))*PIJ
      THETA=TVMS/PKMS
C**** MIX THROUGH SUMSEQUENT LAYERS
      DO L=3,LM
      IF(THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 8160
      PIJ=PLIJ(L,I,J)
      PKMS=PKMS+(PK(L,I,J)*PDSIG(L,I,J))
      THPKMS=THPKMS+  T(I,J,L)*(PK(L,I,J)*PDSIG(L,I,J))
       TXS =  TXS +  TX(I,J,L)*(PK(L,I,J)*PDSIG(L,I,J))
       TYS =  TYS +  TY(I,J,L)*(PK(L,I,J)*PDSIG(L,I,J))
      TXXS = TXXS + TXX(I,J,L)*(PK(L,I,J)*PDSIG(L,I,J))
      TYYS = TYYS + TYY(I,J,L)*(PK(L,I,J)*PDSIG(L,I,J))
      TXYS = TXYS + TXY(I,J,L)*(PK(L,I,J)*PDSIG(L,I,J))
      QMS=    QMS +   Q(I,J,L)*PDSIG(L,I,J)
       QXS =  QXS +  QX(I,J,L)*PDSIG(L,I,J)
       QYS =  QYS +  QY(I,J,L)*PDSIG(L,I,J)
      QXXS = QXXS + QXX(I,J,L)*PDSIG(L,I,J)
      QYYS = QYYS + QYY(I,J,L)*PDSIG(L,I,J)
      QXYS = QXYS + QXY(I,J,L)*PDSIG(L,I,J)
      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*PDSIG(L,I,J))
      THETA=TVMS/PKMS
      END DO
      L=LM+1
 8160 LMAX=L-1
      RDP=1./(P(I,J)*SIGE(1)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      DCLEV(I,J)=LMAX
      DO 8180 L=1,LMAX
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(L,I,J)*PLIJ(L,I,J)
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
      DO L=1,LMAX
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*PDSIG(L,I,J)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*PDSIG(L,I,J)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      DO L=1,LMAX
        PIJ=PLIJ(L,I,J)
        DO K=1,KMAX
          U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &         +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
          V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &         +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
          AJL(IDJ(K),L,38)=AJL(IDJ(K),L,38)
     &         +(UMS(K)-UT(IDI(K),IDJ(K),L))*PIJ*RA(K)
        END DO
      END DO
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
 8400   IF(MODD6.EQ.0) THEN
        DO KR=1,4
          IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
            ADAILY(JHOUR+1,47,KR)=ADAILY(JHOUR+1,47,KR)+1.
            ADAILY(JHOUR+1,48,KR)=ADAILY(JHOUR+1,48,KR)+LMAX
          END IF
        END DO
      END IF
      END IF
      END DO
      END DO
C****
      END DO
      RETURN
 9991 FORMAT ('0SURFACE ',4I4,5F10.4,3F11.7)
 9992 FORMAT ('0',I2,10F10.4/23X,4F10.4,10X,2F10.4/
     *  33X,3F10.4,10X,2F10.4)
 9993 FORMAT ('0',I2,10F10.4/23X,7F10.4/33X,7F10.4)
 9994 FORMAT ('0',I2,11F10.4)
      END
