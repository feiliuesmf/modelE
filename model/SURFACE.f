#include "rundeck_opts.h"

      SUBROUTINE SURFCE
!@sum SURFCE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,rvap,stbo,bygrav,by6,byshi
     *     ,byrhoi,deltx,byrt3
      USE MODEL_COM, only : im,jm,lm,fim,dtsrc,nisurf,u,v,t,p,q
     *     ,idacc,dsig,jday,ndasf,jeq,fland,flice,focean
     *     ,fearth,nday,modrd,itime,jhour,sige,byim,itocean
     *     ,itoice,itlake,itlkice,itlandi,qcheck
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : dxyp,imaxj,kmaxj,ravj,idij,idjj,siniv,cosiv
      USE RADNCB, only : trhr,fsf,cosz1
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs
     &     ,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     &     ,uflux,vflux,tflux,qflux
     &     ,uflux1,vflux1,tflux1,qflux1
C**** Interface to PBL
      USE SOCPBL, only : zgs
     &     ,ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE
     &     ,US,VS,WS,WSH,WSQ,TSV,QS,PSI,DBL,KMS,KHS,KQS,PPBL
     &     ,UG,VG,WG,ZMIX
      USE DAGCOM, only : oa,aij,tdiurn,aj,areg,adiurn,ndiupt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,j_tsrf,j_evap,j_evhdt,j_shdt,j_trhdt,j_f1dt
     *     ,ijdd,idd_spr,idd_pt5,idd_pt4,idd_pt3,idd_pt2,idd_pt1,idd_ts
     *     ,idd_tg1,idd_q5,idd_q4,idd_q3,idd_q2,idd_q1,idd_qs,idd_qg
     *     ,idd_swg,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg
     *     ,idd_wg,idd_us,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq
     *     ,idd_eds,idd_dbl,idd_ev,idd_ldc,idd_dcf
      USE DYNAMICS, only : pmid,pk,pedn,pek,pdsig,plij,am
      USE LANDICE, only : hc2li,z1e,z2li,hc1li
      USE LANDICE_COM, only : snowli
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : xsi,z1i,ace1i,hc1i,alami,byrli,byrls,rhos,kiext
     *     ,ksext
      USE LAKES_COM, only : mwl,mldlk,gml,flake
      USE LAKES, only : minmld
      USE FLUXES, only : dth1,dq1,du1,dv1,e0,e1,evapor,runoe,erunoe
     *     ,solar,dmua,dmva,gtemp,nstype
#ifdef TRACERS_ON
     *     ,tot_trsource
      USE TRACER_COM, only : ntm,itime_tr0,needtrs,trm
#endif
      USE SOIL_DRV, only: earth
      IMPLICIT NONE

      INTEGER I,J,K,IM1,IP1,KR,JR,NS,NSTEPS,MODDSF,MODDD
     *     ,KMAX,IMAX,ITYPE,NGRNDZ,NG,IH
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1,P1K,H0M1
     *     ,PGK,PKDN,DXYPJ,BETAS,EVHDTS,CDMS,CDHS,CDQS,EDS1S,PPBLS
     *     ,EVAPS,DBLS,BETA,ELHX,ACE2,CDTERM,CDENOM,HC1,dF1dTG,HCG1,HCG2
     *     ,DTGRND,EVHDT,F1DT,CM,CH,CQ,DHGS,DQGS,DGS,BETAUP,EVHEAT,F0,F1
     *     ,DSHDTG,DQGDTG,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG,dEVdQS
     *     ,HSDEN,HSCON,HSMUL,dHS,dQS,dT2,dTS,DQ1X,EVHDT0,EVAP,F0DT
     *     ,FTEVAP,VAP,TIMEZ,PWATER,PXSOIL,PSK,TH1,Q1,THV1,TFS
     *     ,RMBYA,Q0M1,TSS,QSS,TAUS,RTAUS,RTAUUS,RTAUVS,TG1S
     *     ,QGS,SRHDTS,TRHDTS,SHDTS,UGS,PTYPE,TG1,SRHEAT,SNOW,TG2,SHDT
     *     ,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,SHEAT,TRHEAT,QSDEN
     *     ,QSCON,QSMUL,T2DEN,T2CON,T2MUL,TGDEN,FQEVAP,ZS1CO,USS
     *     ,VSS,WSS,VGS,WGS,USRS,VSRS,Z2,Z2BY4L,Z1BY6L,QZ1,POC,POI
     *     ,PLK,PLKI,EVAPLIM,F1DTS,HICE,HSNOW,HICE1,HSNOW1,F2,FSRI(2)
     *     ,PSIS,HTLIM   ! THZ1,HZM1,HS,QZM1,

      REAL*8 MSUM, MA1, MSI1, MSI2,tmp
      REAL*8, DIMENSION(NSTYPE,IM,JM) :: TGRND,TGRN2
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: S1BYG1 = BYRT3,
     *     Z1IBYL=Z1I/ALAMI, Z2LI3L=Z2LI/(3.*ALAMI), Z1LIBYL=Z1E/ALAMI
      REAL*8 QSAT,DQSATDT,TOFREZ

#ifdef TRACERS_ON
C**** Tracer input/output common block for PBL
!@var trsfac, trconstflx factors in surface flux boundary cond.
!@var ntx number of tracers that need pbl calculation
      real*8, dimension(ntm) :: trtop,trs,trsfac,trconstflx
      real*8 rhosrf0
      integer itr,n,ntx
      common /trspec/trtop,trs,trsfac,trconstflx,ntx
#endif

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      IH=JHOUR+1

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
         MODDD=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??
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
      EVAPLIM = 0. ; HTLIM=0.  ! need initialisation
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
      THV1=TH1*(1.+Q1*deltx)
         TFS=TF*PXSOIL
      RMBYA=100.*PDSIG(1,I,J)/GRAV
C      THZ1 = TZ(I,J,1) ! vertical gradient of potential temperature
C      QZ1 = QZ(I,J,1) ! vertical gradient of specific humidity
      MSUM = (PS*100.)/GRAV ! total column mass of atmosphere (kg/m^2)
      MA1 = RMBYA ! mass of lowest atmospheric layer (kg/m^2)
      H0M1 = TH1*SHA*MA1*DXYP(J) ! mean pot.enthalpy of lowest atm. (J)
c      HZM1 = THZ1*SHA*MA1*DXYP(J) ! vert. grad. of lowest pot. enth.(J)
      Q0M1 = Q1*MA1*DXYP(J) ! mean water vapor of lowest atmosphere (kg)
c      QZM1 = QZ1*MA1*DXYP(J) ! vert. grad. of lowest layer  vapor (kg)
      PGK = (PS*100.)**KAPA
c      HS = (H0M1-HZM1*S1BYG1)*PGK/(DXYP(J)*MA1) ! pot. spec. enth.(J/kg)
      PKDN = (GRAV*(MSUM-MA1*0.25))**KAPA
#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      n=0
      do itr=1,ntm
        if (itime_tr0(itr).le.itime .and. needtrs(itr)) then
          n=n+1
C**** Calculate first layer tracer concentration
          trtop(n)=trm(i,j,1,itr)/(AM(I,J,1)*DXYP(J))
        end if
      end do
      ntx = n
#endif
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
C**** limit evap/cooling if between MINMLD and 40cm, no evap below 40cm
        IF (MWL(I,J).lt.MINMLD*RHOW*FLAKE(I,J)*DXYP(J)) THEN
          EVAPLIM=MAX(0.5*(MWL(I,J)/(FLAKE(I,J)*DXYP(J))-0.4d0*RHOW),
     *         0d0)
        ELSE
          EVAPLIM=MWL(I,J)/(FLAKE(I,J)*DXYP(J))-(0.5*MINMLD+0.2d0)*RHOW
        END IF
        HTLIM = GML(I,J)/(FLAKE(I,J)*DXYP(J)) + 0.5*LHM*EVAPLIM
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
      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG=0.98d0*QG
      TGV=TG*(1.+QG*deltx)
#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
      n=0
      do itr=1,ntm
        if (itime_tr0(itr).le.itime .and. needtrs(itr)) then
          n=n+1
C**** Calculate trsfac (set to zero for const flux)
          trsfac(n)=0.
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
          rhosrf0=100.*PS/(RGAS*TGV) ! estimated surface density
          trconstflx(n)=tot_trsource(i,j,itr)/(DXYP(J)*rhosrf0)
        end if
      end do
#endif
C =====================================================================
      CALL PBL(I,J,ITYPE,PTYPE)
      CM = cmgs(i,j,itype)
      CH = chgs(i,j,itype)
      CQ = cqgs(i,j,itype)
      DHGS=(ZMIX-ZGS)*CH*WSH
      DQGS=(ZMIX-ZGS)*CQ*WSQ
      DGS =DQGS
C =====================================================================
      TS=TSV/(1.+QS*deltx)
C**** CALCULATE RHOS*CM*WS AND RHOS*CH*WS
3500  CONTINUE
      RHOSRF=100.*PS/(RGAS*TSV)
      RCDMWS=CM*WS*RHOSRF
      RCDHWS=CH*WSH*RHOSRF
      RCDQWS=CQ*WSQ*RHOSRF
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
      DSHDTG=-RCDHWS*KHS*SHA/(DHGS+KHS)
      DQGDTG=QG*DQSATDT(TG,ELHX)
      DEVDTG=-RCDQWS*KHS*LHE*BETAUP*DQGDTG/(DQGS+KHS)
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
      dSNdTG=-RCDHWS*KHS*SHA/(DHGS+KHS)
      dQGdTG=QG*DQSATDT(TG,ELHX) ! d(QG)/dTG
      dEVdTG = -dQGdTG*LHE*RCDQWS*KHS/(DGS+KHS) ! d(EVHEAH)/dTG
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
        if (QCHECK) WRITE(99,*) "Lake EVAP limited: I,J,EVAP,MWL",I,J
     *       ,EVAPOR(I,J,1)-DQ1X*RMBYA, MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J
     *       ))
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
C**** Limit heat fluxes out of lakes if near minimum depth
      IF (ITYPE.eq.1 .and. PLK.gt.0 .and. E0(I,J,1)+F0DT+HTLIM.lt.0)
     *     THEN
        if (QCHECK) write(6,*) "Limiting heat flux from lake",i,j,SHDT
     *       ,F0DT,E0(I,J,1),DTSURF*SRHEAT,TRHDT,EVHDT,HTLIM
        SHDT = -(HTLIM+E0(I,J,1)+DTSURF*SRHEAT+TRHDT+EVHDT)
        F0DT = -E0(I,J,1)-HTLIM
        if (QCHECK) write(6,*) "New SHDT,F0DT",i,j,SHDT,F0DT
      END IF
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
         EDS1S=EDS1S+KHS*PTYPE
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
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
 5800    IF(MODDD.EQ.0) THEN
         DO KR=1,NDIUPT
            IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
              ADIURN(IH,IDD_SPR,KR)=ADIURN(IH,IDD_SPR,KR)+PS
              ADIURN(IH,IDD_PT5,KR)=ADIURN(IH,IDD_PT5,KR)+PSK*T(I,J,5)
              ADIURN(IH,IDD_PT4,KR)=ADIURN(IH,IDD_PT4,KR)+PSK*T(I,J,4)
              ADIURN(IH,IDD_PT3,KR)=ADIURN(IH,IDD_PT3,KR)+PSK*T(I,J,3)
              ADIURN(IH,IDD_PT2,KR)=ADIURN(IH,IDD_PT2,KR)+PSK*T(I,J,2)
              ADIURN(IH,IDD_PT1,KR)=ADIURN(IH,IDD_PT1,KR)+PSK*T(I,J,1)
              ADIURN(IH,IDD_TS,KR)=ADIURN(IH,IDD_TS,KR)+TSS
              ADIURN(IH,IDD_TG1,KR)=ADIURN(IH,IDD_TG1,KR)+(TG1S+TFS)
              ADIURN(IH,IDD_Q5,KR)=ADIURN(IH,IDD_Q5,KR)+Q(I,J,5)
              ADIURN(IH,IDD_Q4,KR)=ADIURN(IH,IDD_Q4,KR)+Q(I,J,4)
              ADIURN(IH,IDD_Q3,KR)=ADIURN(IH,IDD_Q3,KR)+Q(I,J,3)
              ADIURN(IH,IDD_Q2,KR)=ADIURN(IH,IDD_Q2,KR)+Q(I,J,2)
              ADIURN(IH,IDD_Q1,KR)=ADIURN(IH,IDD_Q1,KR)+Q1
              ADIURN(IH,IDD_QS,KR)=ADIURN(IH,IDD_QS,KR)+QSS
              ADIURN(IH,IDD_QG,KR)=ADIURN(IH,IDD_QG,KR)+QGS
              ADIURN(IH,IDD_SWG,KR)=ADIURN(IH,IDD_SWG,KR)+SRHDTS
              ADIURN(IH,IDD_LWG,KR)=ADIURN(IH,IDD_LWG,KR)+TRHDTS
              ADIURN(IH,IDD_SH,KR)=ADIURN(IH,IDD_SH,KR)+SHDTS
              ADIURN(IH,IDD_LH,KR)=ADIURN(IH,IDD_LH,KR)+EVHDTS
              ADIURN(IH,IDD_HZ0,KR)=ADIURN(IH,IDD_HZ0,KR)
     *             +(SRHDTS+TRHDTS+SHDTS+EVHDTS)
              ADIURN(IH,IDD_UG,KR)=ADIURN(IH,IDD_UG,KR)+UGS
              ADIURN(IH,IDD_VG,KR)=ADIURN(IH,IDD_VG,KR)+VGS
              ADIURN(IH,IDD_WG,KR)=ADIURN(IH,IDD_WG,KR)+WGS
              ADIURN(IH,IDD_US,KR)=ADIURN(IH,IDD_US,KR)+USS
              ADIURN(IH,IDD_VS,KR)=ADIURN(IH,IDD_VS,KR)+VSS
              ADIURN(IH,IDD_WS,KR)=ADIURN(IH,IDD_WS,KR)+WSS
              ADIURN(IH,IDD_CIA,KR)=ADIURN(IH,IDD_CIA,KR)+PSIS
              ADIURN(IH,IDD_CM,KR)=ADIURN(IH,IDD_CM,KR)+CDMS
              ADIURN(IH,IDD_CH,KR)=ADIURN(IH,IDD_CH,KR)+CDHS
              ADIURN(IH,IDD_CQ,KR)=ADIURN(IH,IDD_CQ,KR)+CDQS
              ADIURN(IH,IDD_EDS,KR)=ADIURN(IH,IDD_EDS,KR)+EDS1S
              ADIURN(IH,IDD_DBL,KR)=ADIURN(IH,IDD_DBL,KR)+DBLS
              ADIURN(IH,IDD_EV,KR)=ADIURN(IH,IDD_EV,KR)+EVAPS
            END IF
         END DO
       END IF
       IM1=I
      END DO

 7000 CONTINUE
C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODDD)
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
          tmp=q(i,j,1)
          Q(I,J,1)=qmin
          dq1(i,j)=qmin-tmp
          QMOM(:,I,J,1)=0.
        ENDIF
c****   retrieve fluxes
        RMBYA=100.*PDSIG(1,I,j)/GRAv
        P1K=PK(1,I,J)     ! EXPBYK(P1)
        tflux1(i,j)=dth1(i,j)*(-RMBYA*P1K)/(dtsurf)
        qflux1(i,j)=dq1(i,j)*(-RMBYA)/(dtsurf)
C**** Diurnal cycle of temperature diagnostics
        tdiurn(i,j,5)=tdiurn(i,j,5)+(tsavg(i,j)-tf)
        if(tsavg(i,j).gt.tdiurn(i,j,6)) tdiurn(i,j,6)=tsavg(i,j)
        if(tsavg(i,j).lt.tdiurn(i,j,9)) tdiurn(i,j,9)=tsavg(i,j)
      END DO
      END DO
#ifdef TRACERS_ON
C****
C**** Apply tracer surface sources and sinks
C****
      call apply_tracer_source(dtsurf)
#endif
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
c****
      DO J=1,JM
      IMAX=IMAXJ(J)
      DO I=1,IMAX
        RMBYA=100.*PDSIG(1,I,j)/GRAV
        uflux1(i,j)=du1(i,j)*RMBYA/(dtsurf)
        vflux1(i,j)=dv1(i,j)*RMBYA/(dtsurf)
      end do
      end do
C**** Call dry convection or aturb depending on rundeck
      CALL ATM_DIFFUS(1,1,dtsurf)
C****
C**** ACCUMULATE SOME ADDITIONAL BOUNDARY LAYER DIAGNOSTICS
C****
      IF(MODDD.EQ.0) THEN
        DO KR=1,4
C**** CHECK IF DRY CONV HAS HAPPENED FOR THIS DIAGNOSTIC
          IF(DCLEV(IJDD(1,KR),IJDD(2,KR)).GT.1.) THEN
            ADIURN(IH,IDD_DCF,KR)=ADIURN(IH,IDD_DCF,KR)+1.
            ADIURN(IH,IDD_LDC,KR)=ADIURN(IH,IDD_LDC,KR)
     *           +DCLEV(IJDD(1,KR),IJDD(2,KR))
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

