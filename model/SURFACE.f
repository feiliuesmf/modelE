#include "rundeck_opts.h"

      SUBROUTINE SURFCE
!@sum SURFCE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi
     *     ,sha,tf,rhow,rhoi,shv,shw,shi,rvap,stbo,bygrav,by6,byshi
     *     ,byrhoi,deltx,byrt3,teeny
      USE MODEL_COM, only : im,jm,lm,fim,dtsrc,nisurf,u,v,t,p,q
     *     ,idacc,dsig,jday,ndasf,jeq,fland,flice,focean
     *     ,fearth,nday,modrd,itime,jhour,sige,byim,itocean
     *     ,itoice,itlake,itlkice,itlandi,qcheck
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE GEOM, only : dxyp,imaxj,bydxyp
      USE RADNCB, only : trhr,fsf,cosz1
      USE PBLCOM, only : ipbl,cmgs,chgs,cqgs
     &     ,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     &     ,uflux,vflux,tflux,qflux,tgvavg,qgavg
C**** Interface to PBL
      USE SOCPBL, only : zgs,ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE
     &     ,US,VS,WS,WSH,TSV,QS,PSI,DBL,KMS,KHS,KQS,PPBL
     &     ,UG,VG,WG,ZMIX
      USE PBL_DRV, only : pbl,evap_max,fr_sat
      USE DAGCOM, only : oa,aij,tdiurn,aj,areg,adiurn,ndiupt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,j_tsrf,j_evap,j_evhdt,j_shdt,j_trhdt,j_f1dt
     *     ,ijdd,idd_spr,idd_pt5,idd_pt4,idd_pt3,idd_pt2,idd_pt1,idd_ts
     *     ,idd_tg1,idd_q5,idd_q4,idd_q3,idd_q2,idd_q1,idd_qs,idd_qg
     *     ,idd_swg,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg
     *     ,idd_wg,idd_us,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq
     *     ,idd_eds,idd_dbl,idd_ev,idd_ldc,idd_dcf
      USE DYNAMICS, only : pmid,pk,pedn,pek,pdsig,plij,am,byam
      USE LANDICE, only : hc2li,z1e,z2li,hc1li
      USE LANDICE_COM, only : snowli
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : xsi,z1i,ace1i,hc1i,alami,byrli,byrls,rhos,kiext
     *     ,ksext
      USE LAKES_COM, only : mwl,mldlk,gml,flake
      USE LAKES, only : minmld
      USE FLUXES, only : dth1,dq1,e0,e1,evapor,runoe,erunoe
     *     ,solar,dmua,dmva,gtemp,nstype,uflux1,vflux1,tflux1,qflux1
#ifdef TRACERS_ON
     *     ,trsrfflx,trsource
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer
#endif
      USE TRACER_COM, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
      USE TRACER_DIAG_COM, only : taijn,tij_surf
#endif
      USE SOIL_DRV, only: earth
      IMPLICIT NONE

      INTEGER I,J,K,KR,JR,NS,NSTEPS,MODDSF,MODDD,ITYPE,IH,IDTYPE 
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1,P1K,H0M1,PGK,PKDN
     *     ,BETA,ELHX,ACE2,CDTERM,CDENOM,HC1,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,DHGS,DQGS,DGS,BETAUP,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG,dEVdQS,HSDEN,HSCON
     *     ,HSMUL,dHS,dQS,dT2,dTS,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER
     *     ,PXSOIL,PSK,TH1,Q1,THV1,TFS,Q0M1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,SHEAT,TRHEAT
     *     ,QSDEN,QSCON,QSMUL,T2DEN,T2CON,T2MUL,TGDEN,FQEVAP,Z2BY4L
     *     ,Z1BY6L,QZ1,EVAPLIM,HICE,HSNOW,HICE1,HSNOW1,F2,FSRI(2),HTLIM

      REAL*8 MSUM, MA1, MSI1, MSI2
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
      real*8 rhosrf0, totflux
      integer n,nx,ntx,nsrc
      integer, dimension(ntm) :: ntix
      common /trspec/trtop,trs,trsfac,trconstflx,ntx
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: tevaplim,tevap,trgrnd
      real*8  TEV, dTEVdTQS, dTQS, TDP, TDT1
#endif
#endif

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      IH=JHOUR+1

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

C**** Zero out fluxes summed over type and surface time step
      E0=0. ; E1=0. ; EVAPOR=0. ; RUNOE=0. ; ERUNOE=0.
      DMUA=0. ; DMVA=0. ; SOLAR=0.
#ifdef TRACERS_WATER
      TREVAPOR = 0. ; TRUNOE = 0.
#endif

C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NIsurf TIMES EVERY HOUR
C****
      DO NS=1,NIsurf
         MODDSF=MOD(NSTEPS+NS-1,NDASF*NIsurf+1)
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1
         MODDD=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??
C**** ZERO OUT FLUXES ACCUMULATED OVER SURFACE TYPES
         DTH1=0. ;  DQ1 =0. ;  uflux1=0. ; vflux1=0.
#ifdef TRACERS_WATER
         trsrfflx = 0.
#endif

      call loadbl
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=1,JM
      HEMI=1.
      IF(J.LE.JM/2) HEMI=-1.
      POLE= (J.EQ.1 .or. J.EQ.JM)

      DO I=1,IMAXJ(J)

      ! until pbl loops over i,j,itype
      WSAVG(I,J)=0.
      TSAVG(I,J)=0.
      QSAVG(I,J)=0.
      USAVG(I,J)=0.
      VSAVG(I,J)=0.
      TAUAVG(I,J)=0.
      TGVAVG(I,J)=0.
      QGAVG(I,J)=0.

      ! initialize fluxes before calling PBL subroutine
      uflux(I,J)=0.
      vflux(I,J)=0.
      tflux(I,J)=0.
      qflux(I,J)=0.
      EVAPLIM = 0. ; HTLIM=0.  ! need initialisation
#ifdef TRACERS_WATER
      tevaplim = 0.
#endif
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      POICE=RSI(I,J)*PWATER
      POCEAN=PWATER-POICE
      PXSOIL=POCEAN+POICE+PLICE
      PIJ=P(I,J)
      PS=PEDN(1,I,J)
      PSK=PEK(1,I,J)
      P1=PMID(1,I,J)
      P1K=PK(1,I,J) 
      TH1=T(I,J,1)
      Q1=Q(I,J,1)
      THV1=TH1*(1.+Q1*deltx)
      TFS=TF*PXSOIL
      JR=JREG(I,J)
      MA1=AM(1,I,J) !@var MA1 mass of lowest atmospheric layer (kg/m^2)
      MSUM = (PS*100.)/GRAV !@var MSUM total column mass of atmosphere (kg/m^2)
      H0M1 = TH1*SHA*MA1*DXYP(J) !@var H0M1 mean pot.enthalpy of lowest atm. (J)
      Q0M1 = Q1*MA1*DXYP(J) !@var Q0M1 mean water vapor of lowest atmosphere (kg)
      PGK = (PS*100.)**KAPA
      PKDN = (GRAV*(MSUM-MA1*0.25))**KAPA
#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          nx=nx+1
          ntix(nx) = n
C**** Calculate first layer tracer concentration
          trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
        end if
      end do
      ntx = nx
#endif
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
         IF(MODDD.EQ.0) THEN
         DO KR=1,NDIUPT
           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             ADIURN(IH,IDD_SPR,KR)=ADIURN(IH,IDD_SPR,KR)+PS
             ADIURN(IH,IDD_PT5,KR)=ADIURN(IH,IDD_PT5,KR)+PSK*T(I,J,5)
             ADIURN(IH,IDD_PT4,KR)=ADIURN(IH,IDD_PT4,KR)+PSK*T(I,J,4)
             ADIURN(IH,IDD_PT3,KR)=ADIURN(IH,IDD_PT3,KR)+PSK*T(I,J,3)
             ADIURN(IH,IDD_PT2,KR)=ADIURN(IH,IDD_PT2,KR)+PSK*T(I,J,2)
             ADIURN(IH,IDD_PT1,KR)=ADIURN(IH,IDD_PT1,KR)+PSK*T(I,J,1)
             ADIURN(IH,IDD_Q5,KR)=ADIURN(IH,IDD_Q5,KR)+Q(I,J,5)
             ADIURN(IH,IDD_Q4,KR)=ADIURN(IH,IDD_Q4,KR)+Q(I,J,4)
             ADIURN(IH,IDD_Q3,KR)=ADIURN(IH,IDD_Q3,KR)+Q(I,J,3)
             ADIURN(IH,IDD_Q2,KR)=ADIURN(IH,IDD_Q2,KR)+Q(I,J,2)
             ADIURN(IH,IDD_Q1,KR)=ADIURN(IH,IDD_Q1,KR)+Q1
           END IF
         END DO
         END IF
C****
      DO ITYPE=1,3       ! no earth type
      ipbl(i,j,itype)=0
      SELECT CASE (ITYPE)
C****
C**** OCEAN
C****
      CASE (1)

      PTYPE=POCEAN
      IF (PTYPE.gt.0) THEN
      TG1=GTEMP(1,1,I,J)
      IF (FLAKE(I,J).gt.0) THEN
C**** limit evap/cooling if between MINMLD and 40cm, no evap below 40cm
        IF (MWL(I,J).lt.MINMLD*RHOW*FLAKE(I,J)*DXYP(J)) THEN
          EVAPLIM=MAX(0.5*(MWL(I,J)/(FLAKE(I,J)*DXYP(J))-0.4d0*RHOW),
     *         0d0)
        ELSE
          EVAPLIM=MWL(I,J)/(FLAKE(I,J)*DXYP(J))-(0.5*MINMLD+0.2d0)*RHOW
        END IF
#ifdef TRACERS_WATER
C**** limit on tracer evporation from lake
        TEVAPLIM(1:NTX)=EVAPLIM*GTRACER(NTIX(1:NTX),1,I,J)
#endif
        HTLIM = GML(I,J)/(FLAKE(I,J)*DXYP(J)) + 0.5*LHM*EVAPLIM
        IDTYPE=ITLAKE
      ELSE
        IDTYPE=ITOCEAN
      END IF
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(1,I,J)=SOLAR(1,I,J)+DTSURF*SRHEAT
            OA(I,J,5)=OA(I,J,5)+SRHEAT*DTSURF
      BETA=1.
      ELHX=LHE
      END IF
C****
C**** OCEAN ICE
C****
      CASE (2) 

      PTYPE=POICE
      IF (PTYPE.gt.0) THEN
      IF (FLAKE(I,J).gt.0) THEN
        IDTYPE=ITLKICE
      ELSE
        IDTYPE=ITOICE
      END IF
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
      Z2BY4L=ACE2/(RHOI*(4.*ALAMI))
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
      END IF
C****
C**** LAND ICE
C****
      CASE (3)

      PTYPE=PLICE
      IF (PTYPE.gt.0) THEN
      IDTYPE=ITLANDI
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
      END IF
      END SELECT
C****
      IF (PTYPE.gt.0) THEN
C****
C**** BOUNDARY LAYER INTERACTION
C****
      TKV=THV1*PSK
      ZS1=.5*DSIG(1)*RGAS*BYGRAV*TKV*PIJ/PS
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.

      TG=TG1+TF
      QG=QSAT(TG,ELHX,PS)
      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG=0.98d0*QG
      TGV=TG*(1.+QG*deltx)
#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
      do nx=1,ntx
#ifndef TRACERS_WATER
C**** Calculate trsfac (set to zero for const flux)
        trsfac(nx)=0.
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
        rhosrf0=100.*ps/(rgas*tgv) ! estimated surface density
        totflux=0.
        do nsrc=1,ntsurfsrc(ntix(nx))
          totflux = totflux+trsource(i,j,nsrc,ntix(nx))
        end do
        trconstflx(nx)=totflux/(dxyp(j)*rhosrf0)
#else
C**** Set surface boundary conditions for water tracers
C**** trsfac and trconstflx are multiplied by cq*wsh in PBL
        trsfac(nx)=1.
        trconstflx(nx)=gtracer(ntix(nx),itype,i,j)*QG
        trgrnd(nx)=gtracer(ntix(nx),itype,i,j)*QG
#endif
      end do
#endif
C =====================================================================
      fr_sat = 1. ! entire surface is saturated
      evap_max = 1.
      CALL PBL(I,J,ITYPE,PTYPE)
      CM = cmgs(i,j,itype)
      CH = chgs(i,j,itype)
      CQ = cqgs(i,j,itype)
      DHGS=(ZMIX-ZGS)*CH*WSH
      DQGS=(ZMIX-ZGS)*CQ*WSH
      DGS =DQGS
C =====================================================================
      TS=TSV/(1.+QS*deltx)
C**** CALCULATE RHOS*CM*WS AND RHOS*CH*WS
      RHOSRF=100.*PS/(RGAS*TSV)
      RCDMWS=CM*WS*RHOSRF
      RCDHWS=CH*WSH*RHOSRF
      RCDQWS=CQ*WSH*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)
      SHEAT=SHA*RCDHWS*(TS-TG)
      BETAUP = BETA
      IF (QS .GT. QG) BETAUP = 1.
      EVHEAT=(LHE+TG1*SHV)*BETAUP*RCDQWS*(QS-QG)
      TRHEAT=TRHR(1,I,J)-STBO*(TG*TG)*(TG*TG)
C****
      SELECT CASE (ITYPE)

      CASE (1) ! FLUXES USING IMPLICIT TIME STEP FOR OCEAN POINTS
        DSHDTG=-RCDHWS*SHA
        dEVdQS = LHE*RCDQWS
        dHS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A       ((1.+2.*S1BYG1)*DTSURF*PGK*RCDHWS+MA1*PKDN)
        dTS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
     A       (MA1*PKDN*SHA-(1.+2.*S1BYG1)*DTSURF*PGK*DSHDTG)
        dQS = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/
     A       ((1.+2.*S1BYG1)*DTSURF*dEVdQS+MA1*LHE)
        SHDT = DTSURF*(SHEAT-dTS*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+dQS*dEVdQS) ! latent heat flux
        TRHDT=DTSURF*TRHEAT
C**** 
      CASE (2) ! FLUXES USING IMPLICIT TIME STEP FOR ICE POINTS

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
C       dSNdHS = RCDHWS ! d(SHEAT)/dHS - kg/(sec*m^2)
        dEVdQS = LHE*RCDQWS     ! d(EVHEAH)/dQS
        HSDEN = -(1.+2.*S1BYG1)*DTSURF*PGK*BETA*dSNdTG+MA1*PKDN*SHA
        HSCON = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/HSDEN ! (J*sec)/kg
        HSMUL = -(1.+2.*S1BYG1)*DTSURF*PGK*BETA*dSNdTG/HSDEN ! J/(kg*degC)
        QSDEN = (1.+2.*S1BYG1)*BETA*DTSURF*dEVdQS+MA1*LHE
        QSCON = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/QSDEN
        QSMUL = -(1.+2.*S1BYG1)*DTSURF*BETA*dEVdTG/QSDEN
        T2DEN = HCG2+BETA*DTSURF*dF1dTG
        T2CON = DTSURF*(F1-F2)/T2DEN
        T2MUL = BETA*DTSURF*dF1dTG/T2DEN
        TGDEN = HCG1-BETA*DTSURF*(dF0dTG-dF1dTG-
     A       HSMUL*dSNdTG+QSMUL*dEVdQS+T2MUL*dF1dTG) ! W/(m^2*degC)
        dTG = DTSURF*(F0-F1+BETA*
     A       (QSCON*dEVdQS-HSCON*dSNdTG+T2CON*dF1dTG))/TGDEN ! degC
        IF (TG1+dTG .GT. 0.) dTG = -TG1
        dHS = HSCON+HSMUL*dTG   ! (J*sec)/kg
        dQS = QSCON+QSMUL*dTG
        dT2 = T2CON+T2MUL*dTG
        SHDT = DTSURF*(SHEAT+BETA*((dTG-dHS)*dSNdTG)) ! sensible
        EVHDT = DTSURF*(EVHEAT+BETA*(dTG*dEVdTG+dQS*dEVdQS)) ! latent
        TRHDT = DTSURF*(TRHEAT+BETA*dTG*dTRdTG) ! thermal flux (J/m^2)
        F1DT = DTSURF*(F1+BETA*(dTG*dF1dTG-dT2*dF1dTG))
        TG1 = TG1+dTG           ! first layer sea ice temperature (degC)
        TG2 = TG2+dT2           ! second layer sea ice temperature (degC)
        TGRN2(ITYPE,I,J) = TG2
C****
      CASE (3) ! IMPLICIT TIME STEP OVER LANDICE

        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
        DSHDTG=-RCDHWS*KHS*SHA/(DHGS+KHS)
        DQGDTG=QG*DQSATDT(TG,ELHX)
        DEVDTG=-RCDQWS*KHS*LHE*BETAUP*DQGDTG/(DQGS+KHS)
        DTRDTG=-4.*STBO*TG*TG*TG
        DF0DTG=DSHDTG+DEVDTG+DTRDTG
        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HC1-DTSURF*DFDTG)
        SHDT=SHDT+DTSURF*(SHEAT+DTG*DSHDTG)
        EVHDT=EVHDT+DTSURF*(EVHEAT+DTG*DEVDTG)
        TRHDT=TRHDT+DTSURF*(TRHEAT+DTG*DTRDTG)
        F1DT=F1DT+DTSURF*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
        TG1=TG1+DTG

      END SELECT

C**** CALCULATE EVAPORATION
      DQ1X =EVHDT/((LHE+TG1*SHV)*MA1)
      EVHDT0=EVHDT
C**** Limit evaporation if lake mass is at minimum
      IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and.
     *     (EVAPOR(I,J,1)-DQ1X*MA1).gt.EVAPLIM) THEN
        if (QCHECK) WRITE(99,*) "Lake EVAP limited: I,J,EVAP,MWL",I,J
     *     ,EVAPOR(I,J,1)-DQ1X*MA1,MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J))
        DQ1X=(EVAPOR(I,J,1)-EVAPLIM)*BYAM(1,I,J)
      ELSEIF (DQ1X.GT.Q1+DQ1(I,J)) THEN
        DQ1X=(Q1+DQ1(I,J))
      ELSE
        GO TO 3720
      END IF
      EVHDT=DQ1X*(LHE+TG1*SHV)*MA1
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*MA1
#ifdef TRACERS_WATER
C****
C**** Calculate Water Tracer Evaporation
C****                                                            
      DO NX=1,NTX                                                 
        N=NTIX(NX)
        IF (ITYPE.EQ.1) THEN ! OCEAN                             
C**** do calculation implicitly for TQS                          
          TEV=-RCDQWS*(trs(nx)-trgrnd(nx)) !*FRACLK(WS,ITRSPC)
          dTEVdTQS =-RCDQWS                !*FRACLK(WS,ITRSPC) 
          dTQS = -(1.+2.*S1BYG1)*DTSURF*TEV/                     
     *           ((1.+2.*S1BYG1)*DTSURF*dTEVdTQS-MA1)            
          TEVAP(NX)=DTSURF*(TEV+dTQS*dTEVdTQS)                    
        ELSE ! ICE AND LAND ICE
C**** tracer flux is set by source tracer concentration
          IF (EVAP.GE.0) THEN   ! EVAPORATION                        
            TEVAP(NX)=EVAP*trgrnd(nx)/QG
          ELSE                  ! DEW                                
c           TEVAP(NX)=EVAP*trs(nx)/(QS*FRACVL(TG1,ITRSPC))
            TEVAP(NX)=EVAP*trs(nx)/(QS+teeny)
          END IF            
        END IF
C**** Limit evaporation if lake mass is at minimum
        IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and. 
     *       (TREVAPOR(n,1,I,J)+TEVAP(NX).gt.TEVAPLIM(NX))) THEN
          WRITE(99,*) "Lake TEVAP limited: I,J,TEVAP,TMWL",N
     *         ,TREVAPOR(n,1,I,J)+TEVAP(NX),TEVAPLIM(NX)
          TEVAP(NX)= TEVAPLIM(NX)-TREVAPOR(n,1,I,J)
        END IF                                                          
        TDP = TEVAP(NX)*DXYP(J)*ptype
        TDT1 = trsrfflx(I,J,n)*DTSURF
        IF (TRM(I,J,1,n)+TDT1+TDP.lt.1d-5) THEN
          WRITE(99,*) "LIMITING TEVAP",I,J,N,TDP,TRM(I,J,1,n)
          TEVAP(NX) = - (TRM(I,J,1,n)+TDT1-1d-5)/(DXYP(J)*ptype)
          trsrfflx(I,J,n)= - TRM(I,J,1,n)+1d-5
        ELSE                                                            
          trsrfflx(I,J,n)=trsrfflx(I,J,n)+TDP/DTSURF
        END IF                                                          
        TREVAPOR(n,ITYPE,I,J)=TREVAPOR(n,ITYPE,I,J)+TEVAP(NX)
      END DO
#endif
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
C**** Limit heat fluxes out of lakes if near minimum depth
      IF (ITYPE.eq.1 .and. FLAKE(I,J).gt.0 .and.
     *     E0(I,J,1)+F0DT+HTLIM.lt.0) THEN
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
      DTH1(I,J)=DTH1(I,J)-SHDT*PTYPE/(SHA*MA1*P1K)
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
      DMUA(I,J,ITYPE)=DMUA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*US
      DMVA(I,J,ITYPE)=DMVA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*VS
      uflux1(i,j)=uflux1(i,j)+PTYPE*RCDMWS*US
      vflux1(i,j)=vflux1(i,j)+PTYPE*RCDMWS*VS
C****
C**** ACCUMULATE DIAGNOSTICS FOR EACH SURFACE TIME STEP AND ITYPE
C****
        AJ(J,J_EVHDT,IDTYPE)=AJ(J,J_EVHDT,IDTYPE)+EVHDT*PTYPE
        AJ(J,J_SHDT ,IDTYPE)=AJ(J,J_SHDT ,IDTYPE)+SHDT *PTYPE
        AJ(J,J_TRHDT,IDTYPE)=AJ(J,J_TRHDT,IDTYPE)+TRHDT*PTYPE
        IF(MODDSF.EQ.0)
     *       AJ(J,J_TSRF,IDTYPE)=AJ(J,J_TSRF,IDTYPE)+(TS-TF)*PTYPE
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
        IF(JR.LT.24) THEN
          AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+TRHDT*PTYPE*DXYP(J)
          AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+SHDT *PTYPE*DXYP(J)
          AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+EVHDT*PTYPE*DXYP(J)
          AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+EVAP *PTYPE*DXYP(J)
          AREG(JR,J_F1DT )=AREG(JR,J_F1DT )+F1DT *PTYPE*DXYP(J)
          IF(MODDSF.EQ.0)
     *         AREG(JR,J_TSRF)=AREG(JR,J_TSRF)+(TS-TF)*PTYPE*DXYP(J)
        END IF
        IF (PLICE.gt.0) THEN
          AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(TS-TF)
          AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)+SHDT
          AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)+EVHDT
          AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)+TRHDT
        END IF
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
        AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDT*PTYPE
        IF(MODRD.EQ.0) AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)+TRHDT
     *       *PTYPE/DTSRC
        AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHEAT*DTSURF+TRHDT)*PTYPE
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHEAT*DTSURF+TRHDT+SHDT
     *       +EVHDT)*PTYPE
        IF(MODDSF.EQ.0) THEN
          AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WS*PTYPE ! added 12/29/96 -rar-
          AIJ(I,J,IJ_TS)=AIJ(I,J,IJ_TS)+(TS-TF)*PTYPE
          AIJ(I,J,IJ_US)=AIJ(I,J,IJ_US)+US*PTYPE
          AIJ(I,J,IJ_VS)=AIJ(I,J,IJ_VS)+VS*PTYPE
          AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RCDMWS*WS*PTYPE
          AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RCDMWS*US*PTYPE
          AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RCDMWS*VS*PTYPE
          AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QS*PTYPE
        END IF
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
        IF(MODDD.EQ.0) THEN
          DO KR=1,NDIUPT
            IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
              ADIURN(IH,IDD_TS,KR)=ADIURN(IH,IDD_TS,KR)+TS*PTYPE
              ADIURN(IH,IDD_TG1,KR)=ADIURN(IH,IDD_TG1,KR)+(TG1+TF)*PTYPE
              ADIURN(IH,IDD_QS,KR)=ADIURN(IH,IDD_QS,KR)+QS*PTYPE
              ADIURN(IH,IDD_QG,KR)=ADIURN(IH,IDD_QG,KR)+QG*PTYPE
              ADIURN(IH,IDD_SWG,KR)=ADIURN(IH,IDD_SWG,KR)+SRHEAT*DTSURF
     *             *PTYPE
              ADIURN(IH,IDD_LWG,KR)=ADIURN(IH,IDD_LWG,KR)+TRHDT*PTYPE
              ADIURN(IH,IDD_SH,KR)=ADIURN(IH,IDD_SH,KR)+SHDT*PTYPE
              ADIURN(IH,IDD_LH,KR)=ADIURN(IH,IDD_LH,KR)+EVHDT*PTYPE
              ADIURN(IH,IDD_HZ0,KR)=ADIURN(IH,IDD_HZ0,KR)
     *             +(SRHEAT*DTSURF+TRHDT+SHDT+EVHDT)*PTYPE
              ADIURN(IH,IDD_UG,KR)=ADIURN(IH,IDD_UG,KR)+UG*PTYPE
              ADIURN(IH,IDD_VG,KR)=ADIURN(IH,IDD_VG,KR)+VG*PTYPE
              ADIURN(IH,IDD_WG,KR)=ADIURN(IH,IDD_WG,KR)+WG*PTYPE
              ADIURN(IH,IDD_US,KR)=ADIURN(IH,IDD_US,KR)+US*PTYPE
              ADIURN(IH,IDD_VS,KR)=ADIURN(IH,IDD_VS,KR)+VS*PTYPE
              ADIURN(IH,IDD_WS,KR)=ADIURN(IH,IDD_WS,KR)+WS*PTYPE
              ADIURN(IH,IDD_CIA,KR)=ADIURN(IH,IDD_CIA,KR)+PSI*PTYPE
              ADIURN(IH,IDD_CM,KR)=ADIURN(IH,IDD_CM,KR)+CM*PTYPE
              ADIURN(IH,IDD_CH,KR)=ADIURN(IH,IDD_CH,KR)+CH*PTYPE
              ADIURN(IH,IDD_CQ,KR)=ADIURN(IH,IDD_CQ,KR)+CQ*PTYPE
              ADIURN(IH,IDD_EDS,KR)=ADIURN(IH,IDD_EDS,KR)+KHS*PTYPE
              ADIURN(IH,IDD_DBL,KR)=ADIURN(IH,IDD_DBL,KR)+DBL*PTYPE
              ADIURN(IH,IDD_EV,KR)=ADIURN(IH,IDD_EV,KR)+EVAP*PTYPE
            END IF
          END DO
        END IF
C****
#ifdef TRACERS_ON
C**** Save surface tracer concentration whether calculated or not
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          nx=nx+1
          taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)+trs(nx)*ptype
        else
          taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)
     *         +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *         *bydxyp(j),0d0)*ptype
        end if
      end do
#endif
C****
C**** SAVE SOME TYPE DEPENDENT FLUXES
C****
      SELECT CASE (ITYPE)
      CASE (1)  ! ocean
        OA(I,J,6)=OA(I,J,6)+TRHDT
        OA(I,J,7)=OA(I,J,7)+SHDT
        OA(I,J,8)=OA(I,J,8)+EVHDT
C****
      CASE (2)  ! seaice
        IF (TG1.GT.TDIURN(I,J,7)) TDIURN(I,J,7) = TG1
        OA(I,J,9)=OA(I,J,9)+TRHDT
        OA(I,J,10)=OA(I,J,10)+SHDT
        OA(I,J,11)=OA(I,J,11)+EVHDT
C**** 
      CASE (3) ! land ice
        IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
C****
      END SELECT
C****
      END IF
      END DO   ! end of itype loop
 
      END DO   ! end of I loop

      END DO   ! end of J loop
C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODDD)
C****
C**** UPDATE FIRST LAYER QUANTITIES
C****
      DO J=1,JM
      DO I=1,IMAXJ(J)
        FTEVAP=0
        IF (DTH1(I,J)*T(I,J,1).lt.0) FTEVAP=-DTH1(I,J)/T(I,J,1)
        FQEVAP=0
        IF (DQ1(I,J).lt.0.and.Q(I,J,1).gt.0) FQEVAP=-DQ1(I,J)/Q(I,J,1)
! Z-moments should be set from PBL
        TMOM(:,I,J,1) = TMOM(:,I,J,1)*(1.-FTEVAP)
        QMOM(:,I,J,1) = QMOM(:,I,J,1)*(1.-FQEVAP)
        IF ( Q(I,J,1)+DQ1(I,J) .LT. qmin ) THEN
          WRITE(99,*)
     &         ITime,'I,J:',I,J,' Q1:',Q(I,J,1)+DQ1(I,J),'->',qmin
          dq1(i,j)=qmin-q(i,j,1)
          QMOM(:,I,J,1)=0.
        ENDIF
c****   retrieve fluxes
        P1K=PK(1,I,J)
        tflux1(i,j)=dth1(i,j)*(-AM(1,I,J)*P1K)/(dtsurf)
        qflux1(i,j)=dq1(i,j)*(-AM(1,I,J))/(dtsurf)
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
c****
c**** apply earth fluxes to the first layer of the atmosphere
c****  (replaced with dummy subroutine when ATURB is used)
c****
      call apply_fluxes_to_atm(dtsurf)
c****
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
      END DO   ! end of surface time step
      RETURN
C****
      END SUBROUTINE SURFCE

