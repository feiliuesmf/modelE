C**** VDIFF: Vertical Diffusion                                         
C****
C**** Vertical diffusion coefficient depends on wind shear.
C**** Uses TRIDIAG for implicit scheme (MU=1) as in diffuse53.
C**** This version only does diffusion for lowest LDIFM layers.
C****
      SUBROUTINE VDIFF (P,U,V,T,Q,DT1)
      INCLUDE 'BB396M23.COM'
      PARAMETER (LDIFM=LM)
      COMMON /DRGCOM/ AIRX(IM,JM),LMC(IM,JM,2),DEFRM(IM,JM)
      COMMON /WORK3/ PK(IM,JM,LM),UBAR(IM,JM,LM), ! PK is from GWDRAG
     *   RHO(IM,JM,LM)
      COMMON /WORK03/ VKEDDY(IM,JM,LM+1)  ! WORK03 is used in AADVT
      DIMENSION USURF(IM,JM),VSURF(IM,JM)
      DIMENSION TSURF(IM,JM),QSURF(IM,JM)
      DIMENSION PL(0:LDIFM+1),RHOL(0:LDIFM+1),AIRM(LDIFM)
      DIMENSION UL(0:LDIFM+1),VL(0:LDIFM+1),TL(0:LDIFM+1),QL(0:LDIFM+1)
      DIMENSION TE(LDIFM+1)
      DIMENSION PLE(LDIFM+1),RHOE(LDIFM+1),DPE(LDIFM+1),DFLX(LDIFM+1)
      DIMENSION AM(LDIFM),AL(LDIFM),AU(LDIFM),B(LDIFM)
      DIMENSION DU(LDIFM),DV(LDIFM)
      DIMENSION DTEMP(LDIFM),DQ(LDIFM)
      DIMENSION SINI(IM),COSI(IM)
      REAL*8 KMEDGE(LDIFM+1),KHEDGE(LDIFM+1)
      REAL*8 LMEDGE(LDIFM+1),LHEDGE(LDIFM+1),MU
      PARAMETER (MU=1.)
      EQUIVALENCE (BLDATA(1,1,2),TSURF(1,1))
      EQUIVALENCE (BLDATA(1,1,3),QSURF(1,1))
      EQUIVALENCE (BLDATA(1,1,6),USURF(1,1))
      EQUIVALENCE (BLDATA(1,1,7),VSURF(1,1))
      DATA IFIRST/1/
      IF (IFIRST.EQ.1) THEN
        PSFMPT=PSF-PTOP
        IFIRST=0
        DO 20 I=1,IM
        SINI(I)=SIN((I-1)*TWOPI/FIM)
   20   COSI(I)=COS((I-1)*TWOPI/FIM)
        BYRGAS = 1./RGAS
      ENDIF
C        WRITE (*,*) ' VDIFF: MRCH=',MRCH
      G2DT=GRAV*GRAV*DT1
C**** Fill in USURF,VSURF at poles (really BLDATA(6),BLDATA(7))
      DO 60 I=2,IM
      USURF(I,1)=USURF(1,1)*COSI(I)-VSURF(1,1)*SINI(I)
      VSURF(I,1)=VSURF(1,1)*COSI(I)+USURF(1,1)*SINI(I)
      USURF(I,JM)=USURF(1,JM)*COSI(I)+VSURF(1,JM)*SINI(I)
   60 VSURF(I,JM)=VSURF(1,JM)*COSI(I)-USURF(1,JM)*SINI(I)
C**** Calculate RHO(I,J,L)
      DO 115 L=1,LTM
      DO 115 J=1,JM
      IMAX=IM
      IF (J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 110 I=1,IMAX
  110 RHO(I,J,L)=   BYRGAS*(P(I,J)*SIG(L)+PTOP)/(T(I,J,L)*PK(I,J,L))
  115 CONTINUE
      DO 125 L=LS1,LM
      DO 125 J=1,JM
      IMAX=IM
      IF (J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 120 I=1,IMAX
  120 RHO(I,J,L)=   BYRGAS*(PSFMPT*SIG(L)+PTOP)/(T(I,J,L)*PK(I,J,L))
  125 CONTINUE
C**** Fill in T,Q,RHO at poles
      DO 130 L=1,LM
      DO 130 I=1,IM
      T(I,1,L)=T(1,1,L)
      T(I,JM,L)=T(1,JM,L)
      Q(I,1,L)=Q(1,1,L)
      Q(I,JM,L)=Q(1,JM,L)
      RHO(I,1,L)=RHO(1,1,L)
      RHO(I,JM,L)=RHO(1,JM,L)
  130 CONTINUE
C**** Get Vertical Diffusion Coefficient for this timestep
      CALL GETVK (U,V,VKEDDY,LDIFM)
C****
C**** U,V Diffusion
C****
      DO 300 J=2,JM
      I=IM
      DO 300 IP1=1,IM
C     IF(DEFRM(I,J).LE.25.E-6) GO TO 300
C**** Surface values are used for F(0)
      PIJ      =.25*(P(I,J-1) + P(IP1,J-1) + P(I,J) + P(IP1,J))
      PL(0)    =(PIJ+PTOP)
      UL(0)    =.25*(USURF(I  ,J-1) + USURF(IP1,J-1) +
     *               USURF(I  ,J)   + USURF(IP1,J))
      VL(0)    =.25*(VSURF(I  ,J-1) + VSURF(IP1,J-1) +
     *               VSURF(I  ,J  ) + VSURF(IP1,J  ))
      QL(0)    =QSURF(I,J)
      TPHYS    =.25*(TSURF(I  ,J-1) + TSURF(IP1,J-1) +
     *               TSURF(I  ,J)   + TSURF(IP1,J))
      TL(0)    =TPHYS*EXPBYK(PL(0))
      RHOL(0)  =  PL(0)/(RGAS*TPHYS)
      DO 200 L=1,MIN(LDIFM+1,LM)
      IF (L.GE.LS1) PIJ=PSF-PTOP
      PL(L)   =PIJ*SIG(L)+PTOP
      UL(L)   =U(I,J,L)
      VL(L)   =V(I,J,L)
      RHOL(L) =.25*(RHO(I,J-1,L)+RHO(IP1,J-1,L)+RHO(I,J,L)+RHO(IP1,J,L))
      PLE(L)  =PIJ*SIGE(L)+PTOP
  200 CONTINUE
C**** Edge values at LM+1 don't matter since diffusiv flx=F(L)-F(L-1)=0
      IF (L.EQ.LM+1) THEN
         PL(L)   =PIJ*SIGE(L)+PTOP
         UL(L)   =UL(L-1)
         VL(L)   =VL(L-1)
         RHOL(L) =RHOL(L-1)
         PLE(L)  =PIJ*SIGE(L)+PTOP
       ENDIF
      DO 220 L=1,LDIFM
  220 AIRM(L)    =PLE(L)-PLE(L+1)
      DO 230 L=1,LDIFM+1
      KMEDGE(L)=VKEDDY(I,J,L)
      KHEDGE(L)=VKEDDY(I,J,L)
  230 DPE(L)     =PL(L-1)-PL(L)
C**** RHOE is obtained from average of T (=P/RHO)
      DO 240 L=1,LDIFM+1
  240 RHOE(L)=2.*PLE(L)*RHOL(L-1)*RHOL(L)/
     *           (PL(L-1)*RHOL(L)+PL(L)*RHOL(L-1))
C**** Calculate diffusive flux and number of timesteps NDT
      CALL DFUSEF (DPE,RHOE,KMEDGE,G2DT,DFLX,NDT,LDIFM)
      DO 290 N=1,NDT
C**** dq/dt by diffusion as tridiagonal matrix
      CALL DFUSEQ(AIRM,DFLX,UL(0),MU,
     *            AM,AL,AU,B,LDIFM)
      CALL TRIDIAG(Al,Am,AU,B,DU,LDIFM)
      CALL DFUSEQ(AIRM,DFLX,VL(0),MU,
     *            AM,AL,AU,B,LDIFM)
      CALL TRIDIAG(Al,Am,AU,B,DV,LDIFM)
C**** Update model winds
      IF (MRCH.GT.0) THEN
         DO 275 L=1,LM
  275    AJL(J,L,32) = AJL(J,L,32) + DU(L)
      ENDIF
      DO 280 L=1,LM
      U(I,J,L) = U(I,J,L) + DU(L)
  280 V(I,J,L) = V(I,J,L) + DV(L)
  290 CONTINUE
  300 I=IP1
      RETURN
      END SUBROUTINE VDIFF

      SUBROUTINE DFUSEF(DPE,RHOE,DEDGE,G2DT,DFLX,NDT,LM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DPE(LM+1),AIRM(LM),RHOE(LM+1),DEDGE(LM+1),DFLX(LM+1)
      DATA EPS/.3/
      DO 190 L=1,LM+1
  190 DFLX(L)=G2DT*RHOE(L)**2*DEDGE(L)/DPE(L)
      NDT=1
      RETURN
      END SUBROUTINE DFUSEF

      SUBROUTINE DFUSEQ(AIRM,DFLX,F,MU,AM,AL,AU,B,LM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MU
      DIMENSION AIRM(LM),F(0:LM+1),DFLX(LM+1)
      DIMENSION AM(LM),AL(LM),AU(LM),B(LM)
      DO 200 L=1,LM
      BYAIRM = 1./AIRM(L)
      B(L)=     BYAIRM*(DFLX(L+1)*(F(L+1)-F(L))
     *                - DFLX(L)*(F(L)-F(L-1)))
      AM(L)=1.+ (MU*BYAIRM)*(DFLX(L+1)+DFLX(L))
      if (l.lt.lm) AL(L+1)=-(MU*BYAIRM)*DFLX(L)
  200 AU(L)=-(MU*BYAIRM)*DFLX(L+1)
      RETURN
      END SUBROUTINE DFUSEQ

      SUBROUTINE GETVK (U,V,VKEDDY,LDIFM)
C****
C**** Generate vertical diffusion coefficient - use high wind shear
C****
      INCLUDE 'BB396M23.COM'
      DIMENSION VKEDDY(IM,JM,LM+1) ! Vert. Diffusion coefficent
      DIMENSION USURF(IM,JM),VSURF(IM,JM)
      EQUIVALENCE (USURF,BLDATA(1,1,6)),(VSURF,BLDATA(1,1,7))
      DATA IFIRST/1/, XEDDY/10./
      DO 150 J=2,JM
      I=IM
      DO 150 IP1=1,IM
      US=.25*(USURF(I,J-1)+USURF(IP1,J-1)+USURF(I,J)+USURF(IP1,J))
      VS=.25*(VSURF(I,J-1)+VSURF(IP1,J-1)+VSURF(I,J)+VSURF(IP1,J))
      DELV2   =(U(I,J,1)-US)**2 + (V(I,J,1)-VS)**2
      VKEDDY(I,J,1)=0.
      IF (DELV2.GT.25.**2)  VKEDDY(I,J,1)=XEDDY
  150 I=IP1
      DO 200 L=2,MIN(LDIFM+1,LM)
      DO 200 J=2,JM
      DO 200 I=1,IM
      DELV2   =(U(I,J,L)-U(I,J,L-1))**2 + (V(I,J,L)-V(I,J,L-1))**2
      VKEDDY(I,J,L)=0.
      IF (DELV2.GT.25.**2) THEN
         VKEDDY(I,J,L)=XEDDY
C        WRITE (99,901) I,J,L,DELV2
  901 FORMAT (1X,'GETVK:  I,J,L =',3I3,4X,'DELV2=',F12.2)
      ENDIF
  200 CONTINUE
      IF (LDIFM.EQ.LM) THEN
        DO 250 J=2,JM
        DO 250 I=1,IM
  250   VKEDDY(I,J,LM+1)=0.
      ENDIF
      DO 300 J=2,JM
      DO 300 I=1,IM
      VKEDDY(I,J,1)=0.
      VKEDDY(I,J,2)=0.
  300 CONTINUE
      RETURN
      END SUBROUTINE GETVK


      SUBROUTINE GWDRAG (P,U,V,T,SZ,DT1)
C****
C**** THIS SUBROUTINE PUTS A MOMENTUM DRAG IN THE STRATOSPHERE AND
C**** THE MESOSPHERE.
C**** GWDRAG is called from DYNAM with arguments:
C****      P = Pressure (mb) at end of timestep
C****      U,V,T = wind and Potential Temp. at beginning of timestep
C****      DT1 = timestep (s)
C****
C**** LD: incident layer - lowest layer at which wave is allowed to
C****   break - set to LM+1 if wave doesn't exist
C**** MU: momentum flux of wave - conserved if the wave doesn't break
C**** MUB: breaking momentum flux - proportional to -(W-CN)**3
C****    MU is the Greek letter mu, for momentum.
C****    MU and MUB must be negative.
C**** CN: wave phase speed (m/s)
C**** W-CN: (local wind in direction of wave) - (phase speed of wave)
C****    calculated as WMC = (U*UR + V*VR - CN)  must be positive!!!
C**** For propagation, we need MUB < MU < 0. When MUB > 0, the wave
C****    cannot propagate, instead, it breaks and deposits all its
C****    momentum.  The level where MUB changes sign is a critical
C****    level.  W-CN also changes sign.
C**** (UR,VR): unit vector in direction of the wave
C**** There is some trickiness in the way (W-CN), MU, and UR,VR
C**** are gotten for MC waves. (UR,VR) is in the direction of the mean
C**** wind in the source region. But adding +/- 10, 20, and 40 m/s means
C**** sometimes W-CN is negative, unless the direction of (UR,VR) is
C**** changed.  7159-7159.9 fixes this so that (UR,VR) agrees with CN.
C**** NM:
C**** 1   Mountain waves
C**** 2   Shear waves
C**** 3-8 Convective waves
C**** 9   Deformation wave
C****
      INCLUDE 'BB396M23.COM'
      PARAMETER (NM=9)
C     LOGICAL QWRITE
      REAL*8 MDN,MUP,MU,MUR,MUB
      DIMENSION SZ(IM,JM,LM)
      COMMON/WORK3/PK(IM,JM,LM),UBAR(IM,JM,LM),
     *  DKE(IM,JM,LM),
     *  PLE(LM+1),PL(LM),DP(LM),TL(LM),THL(LM),RHO(LM),
     *  BVF(LM),WL(LM),UL(LM),VL(LM),DL(LM),DUT(LM),DVT(LM),
     *  MUB(LM+1,NM),MU(NM),UR(NM),VR(NM),
     *  DQT(LM),DTT(LM),
     *  RDI(LM),DFTL(LM),DFM(LM),DFR(LM),WMC(LM),
     *  RA(4),VARXS(IM),VARXN(IM),VARYS(IM),VARYN(IM)
     *  ,LD(NM),IO(4),JO(4),WT(NM)
     *  ,UEDGE(LM),VEDGE(LM),BYFACS(LM)
      COMMON/DRGCOM/AIRX(IM,JM),LMC(IM,JM,2),DEFRM(IM,JM),
     *   ZVART(IM,JM),ZVARX(IM,JM),ZVARY(IM,JM),ZWT(IM,JM),
     *   PDEF,LDEF,LDEFM
      DIMENSION EK(NM,JM),CN(NM),PKS(LM)
      DIMENSION BYDXYV(JM)
      DATA CN(1)/0./
      DATA IFIRST/1/
C****
      IF (IFIRST.NE.1) GO TO 50
      ERR=1.E-20
      GRAVS=GRAV*GRAV
      G2DT=GRAVS*DT1
      DTHR=2./NDYN
      SHA=RGAS/KAPA
      BYDT1 = 1./DT1
      BYIM=1./IM
      H0=8000.
      FMC=2.E-7
      VARMIN=XCDNST(1)*XCDNST(1)
      XFROUD=1.
C     USEDEF=1.        ! Use Deformation drag if USEDEF = 1
      XLIMIT=.1      ! per timestep limit on mixing and drag
      ROTK=1.5
      RKBY3=ROTK*ROTK*ROTK
      EKS=0.
      DO 5 J=2,JM
      BYDXYV(J)=1./DXYV(J)
      EK1=TWOPI/SQRT(DXYV(J))
      EK2=EK1*SQRT((360./IM)*(180./(JM-1)))
      EKX=(EK2-EK1)/LOG(EK2/EK1)
      IF (EKX.LT.0.) EKX=0.
      EKS=EKS+EKX*DXYV(J)
      EK(1,J)=EKX
    5 EK(2,J)=EKX
      EKS=EKS*IM/AREAG
      DO 6 J=1,JM
      DO 6 N=3,NM
    6 EK(N,J)=EKS
      WRITE (6,970) (J,EK(1,J),J=2,JM)
      WRITE (6,971) EKS
  970 FORMAT ('0  J,EK:',9X,1P,7(I4,E12.2)/,9(1X,8(I4,E12.2)/))
  971 FORMAT ('   AVG EK: ',4X,E12.2)
C**** Coefficients for Radiative Damping (not valid below 500 mb)
      DO 12 L500=1,LM-2
      IF((PSF-PTOP)*SIG(L500)+PTOP.LT.500.) GO TO 13
   12 CONTINUE
   13 L500=L500+1
      WRITE (*,*) ' L500=',L500
C**** LEVEL FOR WIND SHEAR WAVE GENERATION
      LSHR=9    ! ~300 MB
      LD2=19    ! ~0.2 MB
C****
C**** TOPOGRAPHY VARIANCE FOR MOUNTAIN WAVES
C****
      CALL READT (27,0,ZVART,IM*JM*4,ZVART,1)
      REWIND 27
      DO 20 J=JM-1,1,-1
      DO 20 I=1,IM
      ZVART(I,J+1)=ZVART(I,J)
      ZVARX(I,J+1)=ZVARX(I,J)
      ZVARY(I,J+1)=ZVARY(I,J)
   20 ZWT(I,J+1)=ZWT(I,J)
      DO 40 L=LS1,LM
   40 PKS(L)=((PSF-PTOP)*SIG(L)+PTOP)**KAPA
   50 IFIRST=0
C**** P**KAPA IN THE STRATOSPHERE
      DO 55 L=LS1,LM
      DO 52 J=1,JM
      IMAX=IM
      IF (J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 52 I=1,IMAX
      PK(I,J,L)=PKS(L)
   52 CONTINUE
   55 CONTINUE
C**** P**KAPA IN THE TROPOSPHERE
C     DOPK??
      DO 65 L=1,LS1-1
      DO 60 J=1,JM
      IMAX=IM
      IF (J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 60 I=1,IMAX
      PK(I,J,L)=EXPBYK(P(I,J)*SIG(L)+PTOP)
   60 CONTINUE
   65 CONTINUE
C****
C**** FILL IN QUANTITIES AT POLES
C****
      DO 70 I=2,IM
      AIRX(I,1)=AIRX(1,1)
      LMC(I,1,1)=LMC(1,1,1)
      LMC(I,1,2)=LMC(1,1,2)
      AIRX(I,JM)=AIRX(1,JM)
      LMC(I,JM,1)=LMC(1,JM,1)
   70 LMC(I,JM,2)=LMC(1,JM,2)
      DO 80 L=1,LM
      DO 80 I=2,IM
      T(I,1,L)=T(1,1,L)
      SZ(I,1,L)=SZ(1,1,L)
      PK(I,1,L)=PK(1,1,L)
      T(I,JM,L)=T(1,JM,L)
      SZ(I,JM,L)=SZ(1,JM,L)
   80 PK(I,JM,L)=PK(1,JM,L)
C****
C**** DEFORMATION
C****
      IF(MRCH.EQ.0)  CALL DEFORM (P,U,V)
      DO 85 L=1,LM
      DO 85 J=1,JM
      DO 85 I=1,IM
      DKE(I,J,L)=0.
   85 CONTINUE
C****
C**** BEGINNING OF OUTER LOOP OVER I,J
C****
      DO 500 J=2,JM
      CN(1)=0.
      FCORU=(ABS(FCOR(J-1))+ABS(FCOR(J)))*BYDXYV(J)
      DO 90 K=1,2
      RA(K)=RAVPN(J-1)
      RA(K+2)=RAVPS(J)
      JO(K)=J-1
   90 JO(K+2)=J
      IF (J.NE.2) GO TO 92
      RA(1)=RA(1)*BYIM
      RA(2)=RA(2)*BYIM
   92 IF (J.NE.JM) GO TO 94
      RA(3)=RA(3)*BYIM
      RA(4)=RA(4)*BYIM
   94 I=IM
      DO 500 IP1=1,IM
C     QWRITE=.FALSE.
C     IF (M.EQ.0.AND.J.EQ.23) QWRITE=.TRUE.
C     IF (QWRITE) WRITE (6,980) I,J,(U(I,J,L),L=1,LM)
C 980 FORMAT (1H0,' U BEFORE DRAG AT ',2I6/(1X,12F10.3))
C****
C**** CALCULATE VERTICAL ARRAYS
C****
      PIJ=(P(I,J-1)+P(IP1,J-1))*RAPVN(J-1)+(P(I,J)+P(IP1,J))*RAPVS(J)
      SP=PIJ
      DO 110 L=1,LM
      IF (L.GE.LS1) PIJ=PSF-PTOP
      TL(L)=.25*(PK(I,J-1,L)*T(I,J-1,L)+PK(IP1,J-1,L)*T(IP1,J-1,L)+
     *  PK(I,J,L)*T(I,J,L)+PK(IP1,J,L)*T(IP1,J,L))
      THL(L)=.25*(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))
      PLE(L)=PIJ*SIGE(L)+PTOP
      PL(L)=PIJ*SIG(L)+PTOP
      DP(L)=PIJ*DSIG(L)
      RHO(L)=PL(L)/(RGAS*TL(L))
      BVFSQ=.5*(SZ(I,J-1,L)+SZ(IP1,J-1,L)+SZ(I,J,L)+SZ(IP1,J,L))/
     *  (DP(L)*THL(L))*GRAVS*RHO(L)
      IF (PL(L).GE..4d0) THEN
         BVF(L)=SQRT(MAX(BVFSQ,1.d-10))
        ELSE
         BVF(L)=SQRT(MAX(BVFSQ,1.d-4))
      END IF
CRAD  RDI(L)=960.*960./(TL(L)*TL(L))*EXP(-960./TL(L))
      DL(L)=0.
      DUT(L)=0.
      DVT(L)=0.
      UL(L)=U(I,J,L)
      VL(L)=V(I,J,L)
      WL(L)=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
  110 CONTINUE
      PLE(LM+1)=PIJ*SIGE(LM+1)+PTOP
      PIJ=SP
C****
C**** INITIALIZE THE MOMENTUM FLUX FOR VARIOUS WAVES
C****
      DO 120 N=1,NM
      MU(N)=0.
      LD(N)=LM+1     ! set drag level = LM+1 (turns off drag)
  120 WT(N)=1.       ! initialize area weight to 1
      DO 122 N=1,NM
      DO 122 L=1,LM+1
  122 MUB(L,N)=0.
C**** MOUNTAIN WAVES generate at 1 s.d. above topography ...
      LD(1)=L500
      U0=(UL(1)*DSIG(1)+UL(2)*DSIG(2))/(SIGE(1)-SIGE(3))
      V0=(VL(1)*DSIG(1)+VL(2)*DSIG(2))/(SIGE(1)-SIGE(3))
      W0=SQRT(U0*U0+V0*V0)
      UR(1)=U0/(W0+ ERR)
      VR(1)=V0/(W0+ ERR)
      BV0=(BVF(1)*DSIG(1)+BVF(2)*DSIG(2))/(SIGE(1)-SIGE(3))
      ZVAR=ABS(UR(1))*ZVARX(I,J)+ABS(VR(1))*ZVARY(I,J)
         IF(ZVAR.LT.VARMIN) ZVAR=0.
C.... if Froude number (U0/BV0*ZSD) > 1
C.... limit ZSD to be consistent with Froude no. (U0/BV0*ZSD) > 1
      IF (ZVAR.GT.(XFROUD*W0/BV0)**2) ZVAR=(XFROUD*W0/BV0)**2
      P0=(SIG(1)*DSIG(1)+SIG(2)*DSIG(2))/(SIGE(1)-SIGE(3))*PIJ+PTOP
      WT(1)=ZWT(I,J)
      MU(1)=-.5*EK(1,J)/(H0*ROTK)*P0*BV0*W0*ZVAR
      IF(MU(1)*(UL(L500)*UR(1)+VL(L500)*VR(1)).GE.0.) MU(1)=0.
C     IF(QWRITE) WRITE (6,996) I,J,U0,MU(1),UL(L500),P0
C**** DEFORMATION WAVE (X    d cm-2)
      IF (9.GT.NM)  GO TO 155
      IF (FDATA(I,J,1)/GRAV.GT.1000.) GO TO 155
      IF (DEFRM(I,J).LT.15.E-6) GO TO 155         !  Threshold= 15e-6
      FDEFRM=- 3.*GRAV/(1000.*15.E-6)*DEFRM(I,J)  !  3 d cm-2 @ 15e-6
      DU=UL(LDEF +1)-UL(LDEF )
      DV=VL(LDEF +1)-VL(LDEF )
      DW=SQRT(DU**2        + DV**2       )
      UR(9)=DU       /(DW+ ERR)
      VR(9)=DV       /(DW+ ERR)
      CN(9)=UL(LDEF )*UR(9)+VL(LDEF )*VR(9)
      MU(9)=FDEFRM
      LD(9)=L500
CW    WRITE (99,*) 'DEFORM: I,J,C,MU,DEFRM=',I,J,CN(1),MU(1),DEFRM(I,J)
C**** WIND SHEAR: USE SHEAR BETWEEN 7 AND 8 UNLESS CRIT. LEVEL ABOVE..
  155 IF (2.GT.NM)  GO TO 200
      LN=LSHR
  160 L=LN
      CU=.5*(UL(L)+UL(L+1))
  170 LN=LN+1
      IF (LN.GE.LD2) GO TO 172
      IF ((UL(LN)-CU)*(UL(LN+1)-CU).LT.0.) GO TO 160
      GO TO 170
  172 LD(2)=L+2
      DU=UL(L+1)-UL(L)
      DV=VL(L+1)-VL(L)
      DW=SQRT(DU*DU+DV*DV)
      UR(2)=DU/(DW+ ERR)
      VR(2)=DV/(DW+ ERR)
      CN(2)=.5*((UL(L+1)+UL(L))*UR(2)+(VL(L+1)+VL(L))*VR(2))
      MU(2)=-FCORU*PLE(L+1)*DW*DW/(240.*H0*(BVF(L+1)))
C     IF (QWRITE) WRITE (6,997) I,J,L,LN,DU2K,MU(2),BVF(L+1)
C**** MOIST CONVECTIVE MASS FLUX BEGINS TWO LEVELS ABOVE CLOUD...
C**** AMPLITUDE DEPENDS ON |U(SOURCE)-C|.
C**** NOTE:  NM LE 2, NM EQ 4, NM GE 8  ARE ALLOWED FOR MC DRAG
      IF (4.GT.NM)  GO TO 200
      USRC=0.
      VSRC=0.
      AIRX4=AIRX(I,J-1)+AIRX(IP1,J-1)+AIRX(I,J)+AIRX(IP1,J)
      AIRXS= ((AIRX(I,J-1)+AIRX(IP1,J-1))*RAVPN(J-1)
     *     +  (AIRX(I,J)+AIRX(IP1,J))*RAVPS(J))
      IF (AIRX4.LE.0.) GO TO 200
      LMC0=.5+(LMC(I,J-1,1)*AIRX(I,J-1)+LMC(IP1,J-1,1)*AIRX(IP1,J-1)+
     *  LMC(I,J,1)*AIRX(I,J)+LMC(IP1,J,1)*AIRX(IP1,J))/(AIRX4+ ERR)
C**** Note: LMC1 was defined in CB245M31 as LMAX+1
      LMC1=.5+(LMC(I,J-1,2)*AIRX(I,J-1)+LMC(IP1,J-1,2)*AIRX(IP1,J-1)+
     *  LMC(I,J,2)*AIRX(I,J)+LMC(IP1,J,2)*AIRX(IP1,J))/(AIRX4+ ERR)
      IF (LMC1.LE.4) GO TO 200
      NMX=4
Cwmc     IF(M.EQ.0.AND.J.EQ.23) QWRITE=.TRUE.
      CLDDEP=PIJ*(SIGE(LMC0)-SIGE(LMC1))
      FPLUME=AIRXS/(DXYV(J)*CLDDEP)
      CLDHT=H0*LOG((PIJ*SIGE(LMC0)+PTOP)/(PIJ*SIGE(LMC1)+PTOP))
      WTX=FPLUME
      IF (WTX.GT.1. OR. WTX.LT.0.) THEN
        PRINT *, 'WARNING IN GWDRAG, WTX INCORRECT',WTX,FPLUME,CLDDEP
        IF (WTX.GT.1.) WTX=1.
        IF (WTX.LT.0.) STOP ' WTX <0 IN GWDRAG'
      END IF
      DO 177 L=LMC0,LMC1-1
      USRC=USRC+UL(L)
  177 VSRC=VSRC+VL(L)
      USRC=USRC/(LMC1-LMC0)
      VSRC=VSRC/(LMC1-LMC0)
      WSRC=SQRT(USRC*USRC+VSRC*VSRC)
      UR(3)=USRC/(WSRC+ ERR)
      VR(3)=VSRC/(WSRC+ ERR)
      MU(3)=-EK(3,J)*FMC*BVF(LMC1-1)*PL(LMC1-1)*CLDHT**2
      MU(4)=MU(3)
      CN(3)=WSRC-10.
      CN(4)=WSRC+10.
      UR(4)=UR(3)
      VR(4)=VR(3)
      LD(3)=10
      LD(4)=10
      WT(3)=WTX
      WT(4)=WTX
C     IF (QWRITE) WRITE (6,998) I,J,LMC0,LMC1,CLDHT,WT(3),
C    *  MU(3),CN(3),USRC,VSRC
C     IF (QWRITE) WRITE (*,'(''AIRX,FPL,WTX,CLDDP,DXYV='',1P,5E12.2)')
C    *      AIRXS,FPLUME,WTX,CLDDEP,DXYV(J)
      IF (LMC1.GT.9.AND.NM.GE.8) THEN
      NMX=8
      DO 182 N=3,NMX
      WT(N)=WTX
  182 LD(N)=LMC1+1
      CN(5)=WSRC-40.
      CN(6)=WSRC+40.
      CN(7)=WSRC-20.
      CN(8)=WSRC+20.
      DO 184 N=5,NMX
      MU(N)=MU(3)
      UR(N)=UR(3)
  184 VR(N)=VR(3)
C     ELSE
Cwmc  QWRITE=.FALSE.
      ENDIF
  190 WCHECK=UL(LD(3))*UR(3)+VL(LD(3))*VR(3)
      DO 195 N=3,NMX
      IF (WCHECK.GT.CN(N)) GO TO 195
      UR(N)=-UR(N)
      VR(N)=-VR(N)
      CN(N)=-CN(N)
  195 CONTINUE
C****
C**** BREAKING MOMENTUM FLUX AT LAYER EDGES
C****
  200 CONTINUE
      DO L=2,LM
        UEDGE(L)=.5*(UL(L-1)+UL(L))
        VEDGE(L)=.5*(VL(L-1)+VL(L))
        TEDGE=.5*(TL(L-1)+TL(L))
        BVEDGE=.5*(BVF(L-1)+BVF(L))
        BYFACS(L)=-.5*GRAV*PLE(L)/(RGAS*RKBY3*BVEDGE*TEDGE)
      END DO
      DO 210 N=1,NM
        DO 205 L=LD(N),LM
        WMCE=UEDGE(L)*UR(N)+VEDGE(L)*VR(N) - CN(N)
        IF (WMCE.GE.0.) THEN
          MUB(L,N)=BYFACS(L)*EK(N,J)*WMCE**3
          IF (MUB(L,N).GE.-ERR) MUB(L,N)=0.
        END IF
  205   CONTINUE
  210 CONTINUE
C**** DEPOSIT REMAINING MOM. FLUX IN TOP LAYER
      DO 215 N=1,NM
  215 MUB(LM+1,N)=0.
C**** DISTRIBUTE CRIT LEVEL NEAR TOP
      DO 220 N=1,NM
      IF (MUB(LM,N).EQ.0.) MUB(LM,N)=.3*MUB(LM-1,N)
  220 CONTINUE
Cw*** APPLY AREA WEIGHTING OF MTN WAVE TO BREAKING MOM. FLUX
Cw    DO 225 L=LD(1),LM
Cw225 MUB(L,1)=MUB(L,1)*ZWT(I,J)
Cw*** APPLY AREA WEIGHTING OF MC WAVE TO BREAKING MOM. FLUX
Cw    DO 230 N=3,NMX
Cw    DO 230 L=LD(N),LM
Cw    MUB(L,N)=MUB(L,N)*WT(N)
Cw230 CONTINUE
C****
C**** DETERMINE INCIDENT MOMENTUM FLUX
C****
C**** INCIDENT FLUX FOR MTN WAVES
      LN=LD(1)
      IF (MU(1).LT.MUB(LN,1)) MU(1)=MUB(LN,1)
      IF (LN.LE.LM.AND.BVF(LN).LT.1.E-5) MU(1)=0.
C**** INCIDENT FLUX FOR SHR AND MC WAVES
      DO 250 N=2,NM
      LN=LD(N)
      IF (LN.GT.LM) GO TO 250
      IF (MU(N).LT.MUB(LN,N)) MU(N)=MUB(LN,N)
  250 CONTINUE
C**** LOCATE LDRAG - LOWEST LAYER TO APPLY DRAG
      LDRAG=LM+1
      DO 252 N=1,NM
      IF (MU(N).LE.-ERR) LDRAG=MIN(LD(N),LDRAG)
  252 CONTINUE
C**** DIAGNOSTICS
C     IF (QWRITE) WRITE (6,983) (N,N=1,6),(EK(N,J),N=1,6),
C    *  (UR(N),N=1,6),(VR(N),N=1,6),
C    *  (CN(N),N=1,6),(MU(N),N=1,6),(LD(N),N=1,6)
C 983 FORMAT (1H0,10X,6I11/4X,'EK',4X,3P,6F11.4/
C    * 4X,'UR',5X,0P,6F11.3/4X,'VR',5X,6F11.3/4X,'CN',5X,6F11.1/
C    * 4X,'MU',8X,1P,6E11.2/4X,'LD',10X,6I11)
C     IF (QWRITE) WRITE (6,994) I,J,(L,UL(L),VL(L),(MUB(L,N),N=1,6),
C    *  PLE(L),BVF(L),L=1,LM+1)
         IF (MRCH.NE.2) GO TO 251
         AIJ(I,J,46)=AIJ(I,J,46)+MU(9)*UR(9)*DTHR
         AIJ(I,J,92)=AIJ(I,J,92)+MU(1)*UR(1)*DTHR  *WT(1)
         AIJ(I,J,93)=AIJ(I,J,93)+MU(2)*UR(2)*DTHR
         AIJ(I,J,94)=AIJ(I,J,94)+MU(3)*UR(3)*DTHR  *WT(3)
         AIJ(I,J,95)=AIJ(I,J,95)+MU(7)*UR(7)*DTHR  *WT(7)
         AIJ(I,J,96)=AIJ(I,J,96)+MU(5)*UR(5)*DTHR  *WT(5)
         AIJ(I,J,97)=AIJ(I,J,97)+CN(2)*UR(2)*DTHR
         AIJ(I,J,98)=AIJ(I,J,98)+USRC*DTHR
C****
C**** CALCULATE THE DRAG
C****
  251 IF (LDRAG.GT.LM) GO TO 500
      DO 400 N=1,NM
      IF (LD(N).GT.LM) GO TO 400
      LD1=LD(N)
      WMC(LD1-1)=UL(LD1-1)*UR(N)+VL(LD1-1)*VR(N)-CN(N)
      DO 260 L=LD1,LM-1
      DFM(L)=0.
      DFR(L)=0.
      WMC(L)=UL(L)*UR(N)+VL(L)*VR(N)-CN(N)
      IF (WMC(L).GT..01) GO TO 255
      WMC(L)=.5*(WMC(L-1)+WMC(L))
      IF (WMC(L).LT..01) WMC(L)=.01
  255 LTOP=L
      IF (MUB(L+1,N).EQ.0.) GO TO 270
  260 CONTINUE
      LTOP=L
      DFM(L)=0.
      DFR(L)=0.
      WMC(L)=UL(L)*UR(N)+VL(L)*VR(N)-CN(N)
      IF (WMC(L).GT..01) GO TO 270
      WMC(L)=.5*(WMC(L-1)+WMC(L))
      IF (WMC(L).LT..01) WMC(L)=.01
  270 CONTINUE
Cd    IF (QWRITE) WRITE (6,988) LD1,LTOP,RAREA
C 988 FORMAT (1X,'LD1=',I6,4X,'LTOP=',I6,4X,'RAREA=',1P,E12.3)
      DO 300 L=LD1,LTOP
      IF (L.EQ.LM.AND.MRCH.EQ.2)
     *   AIJ(I,J,99)=AIJ(I,J,99)+MU(N)*UR(N)*DTHR  *WT(N)
C**** RADIATIVE DAMPING
      MUR=MU(N)
CRAD  DTW=DP(L)*BVF(L)/(EK(N,J)*WMC(L)*WMC(L)*RHO(L)*GRAV+ ERR)
CRAD  WAVEM=1.E3*BVF(L)/(ABS(WMC(L))+ ERR)
CRAD  IF (WAVEM.GT.1.6) WAVEM=1.6
CRAD  WM3BY2=WAVEM*SQRT(WAVEM)
CRAD  DTR=86400./(RDI(L)*(AZ(L)+BZ(L)*WAVEM*WAVEM/(CZ(L)+WM3BY2))+ERR)
CRAD  MUR=MU(N)*EXP(-2.*DTW/(DTR+ ERR))
      DFR(L)=MU(N)-MUR
C**** MECHANICAL (TURBULENT) DRAG
      MU(N)=MUR
      IF (N.NE.1.AND.N.NE.9) THEN    ! not mtn or deformation
      IF (MUR.LT.MUB(L+1,N)) MU(N)=MUB(L+1,N)  ! saturation drag
      ELSE
         IF (MUR.LT.MUB(L+1,N)) MU(N)=0.     ! sensitivity test
      ENDIF
      DFM(L)=MUR-MU(N)
C     IF (QWRITE) WRITE (6,986) L,N,WMC(L),DFR(L),DFM(L),MU(N),DP(L)
C 986 FORMAT (1X,'L,N=',2I4,' WMC=',F9.1,' DFR,DFM,MU=',1P,3E12.3,
C    *  '  DP=',E12.3)
  300 CONTINUE
C**** LIMIT THE CONVERGENCE TO XLIMIT*(U-C)
  320 DO 350 L=LTOP,LD1,-1
      DFT=DFR(L)+DFM(L)
      DFMAX=-XLIMIT*WMC(L)*DP(L)*BYDT1
      IF (DFT.GE.DFMAX) GO TO 350
      EXCESS=DFT-DFMAX
      ALFA=DFM(L)/DFT
      DFR(L)=DFR(L)-EXCESS*(1.-ALFA)
      DFM(L)=DFM(L)-EXCESS*ALFA
      IF (L.GT.LD1) THEN
        DFR(L-1)=DFR(L-1)+EXCESS*(1.-ALFA)
        DFM(L-1)=DFM(L-1)+EXCESS*ALFA
      ENDIF
  350 CONTINUE
C**** COMBINE MECHANICAL AND RADIATIVE DRAG
      DO 360 L=LD1,LTOP
  360 DFTL(L)=DFM(L)+DFR(L)
C**** ACCUMULATE SATURATION FLUX FOR MTN WAVES DURING BRKING
C     IF (N.NE.1.OR.MRCH.NE.2) GO TO 375
C     DO 370 L=LD1,LTOP
C     IF (DFM(L).EQ.0.) GO TO 370
C        AJL(J,L,30)=AJL(J,L,30)+ZWT(I,J)
C        AJL(J,L,18)=AJL(J,L,18)+MUB(L+1,1)*UR(1)*ZWT(I,J)
C 370 CONTINUE
C**** CALCULATE DIFFUSION COEFFICIENT AND ADD DRAG TO THE WINDS
C**** (DEFORMATION DIFFUSION IS * XDIFF)
      XDIFF=1.
      DO 390 L=LD1,LTOP
      DL(L)=DL(L)+ABS(DFM(L)*WMC(L))*XDIFF*WT(N)
      DFT=DFTL(L)
      DWT=DFT*DT1/DP(L)
      IF (DFT.GT.-ERR) GO TO 390
C**** SHEAR AND MOUNTAIN ORIENTATION
      DUTN=DWT*UR(N)*WT(N)
      DUT(L)=DUT(L)+DUTN
      DVT(L)=DVT(L)+DWT*VR(N)
C     GO TO 380
C 380 IF (QWRITE) WRITE (6,985) L,N,DWT,UL(L),VL(L),WL(L),DFT
C 985 FORMAT (1X,'L,N=',2I4,'  DWT=',1P,E12.3,'  UL,VL,WL=',
C    *  3E12.3,'  DFT=',E12.3)
      IF (MRCH.NE.2) GO TO 390
        IF (N.LT.9) THEN
         AJL(J,L,N+19)=AJL(J,L,N+19)+DUTN
        ELSE
         AJL(J,L,18)=AJL(J,L,18)+DUTN
        ENDIF
  390 CONTINUE
  400 CONTINUE
C     IF (QWRITE) WRITE (6,987) (L,DUT(L),DVT(L),DL(L),L=LDRAG,LM)
C 987 FORMAT (1H0,'  L      DUT         DVT         DL'/
C    *  (1X,I4,1P,3E12.3))
C****
C**** MOMENTUM DIFFUSION   (DOUBLED)
C****    (limited to XLIMIT per timestep)
      PDN=PL(LDRAG-1)
      MDN=DP(LDRAG-1)
      YDN=0.
      FLUXUD=0.
      FLUXVD=0.
      DO 430 L=LDRAG,LM
      PUP=PL(L)
      MUP=DP(L)
      YUP=G2DT*RHO(L)*RHO(L)*DL(L)/(DP(L)*BVF(L)*BVF(L))
      DX=   (YDN+YUP)/(PDN-PUP)     ! double diffusion coefficient
      DLIMIT=XLIMIT*MIN(MUP,MDN)
      IF (DX.GT.DLIMIT) THEN
C       IF (DX.GT.DLIMIT)
C    *   WRITE (99,9430) MRCH,I,J,L,DX,MUP,MDN,UL(L),UL(L-1)
C9430 FORMAT (' 7267 DIFFX > LIMIT:',4I3,1P,3E10.2,0P,2F7.1)
         DX=DLIMIT
      ENDIF
      FLUXU=DX/(1.+DX*(MUP+MDN)/(MUP*MDN))*(UL(L)-UL(L-1))
      FLUXV=DX/(1.+DX*(MUP+MDN)/(MUP*MDN))*(VL(L)-VL(L-1))
      DUT(L-1)=DUT(L-1)-(FLUXUD-FLUXU)/MDN
      DVT(L-1)=DVT(L-1)-(FLUXVD-FLUXV)/MDN
C        IF (MRCH.EQ.2) AJL(J,L-1,32)=AJL(J,L-1,32)-(FLUXUD-FLUXU)/MDN
      FLUXUD=FLUXU
      FLUXVD=FLUXV
      YDN=YUP
      MDN=MUP
      PDN=PUP
  430 CONTINUE
C**** DIFFUSION IN THE TOP LAYER COMES ONLY FROM BELOW.
      DUT(LM)=DUT(LM)-FLUXUD/MDN
      DVT(LM)=DVT(LM)-FLUXVD/MDN
C****    ACCUMULATE DIAGNOSTICS  (DU, DIFFUSION COEFFICIENT)
         IF (MRCH.EQ.2) THEN
C        AJL(J,LM,32)=AJL(J,LM,32)-FLUXUD/MDN
         DO 440 L=LDRAG,LM
  440    AJL(J,L,31)=AJL(J,L,31)+DL(L)/(BVF(L)*BVF(L))*DTHR
         ENDIF
C****
C**** Save KE change and diffusion coefficient on A-grid
C****
      IF (MRCH.NE.2) GO TO 470
      DO 450 K=1,2
      IO(2*K-1)=I
  450 IO(2*K)=IP1
      IF (J.NE.2) GO TO 452
      IO(1)=1
      IO(2)=1
  452 IF (J.NE.JM) GO TO 454
      IO(3)=1
      IO(4)=1
  454 CONTINUE
      DO 460 L=LDRAG-1,LM
      DKEX=.5*(((DUT(L)+U(I,J,L))**2+(DVT(L)+V(I,J,L))**2)-
     *     (U(I,J,L)**2+V(I,J,L)**2))
      DO 460 K=1,4
      IA=IO(K)
      JA=JO(K)
      DKE(IA,JA,L)=DKE(IA,JA,L) + RA(K)*DKEX
  460 CONTINUE
C****
C**** UPDATE THE U AND V WINDS
C****
  470 DO 480 L=LDRAG-1,LM
      U(I,J,L)=U(I,J,L)+DUT(L)
  480 V(I,J,L)=V(I,J,L)+DVT(L)
C     IF (QWRITE) WRITE (6,951) (U(I,J,L),L=1,LM)
C 951 FORMAT (1H0,'U AFTER DRAG'/(1X,12F10.3))
C**** END OF LOOP OVER I,J
  500 I=IP1
      IF (MRCH.NE.2) RETURN
C****
C**** PUT THE KINETIC ENERGY BACK IN AS HEAT
C****
      DO 600 L=LDRAG-1,LM
      DO 560 J=1,JM
      IMAX=IM
      IF (J.EQ.1.OR.J.EQ.JM) IMAX=1
      DO 550 I=1,IMAX
      T(I,J,L)=T(I,J,L)-DKE(I,J,L)/(SHA*PK(I,J,L))
  550    AJL(J,L,33)=AJL(J,L,33)-DKE(I,J,L)/(SHA*PK(I,J,L))
  560 CONTINUE
  600 CONTINUE
      RETURN
  911 FORMAT ('0TOPOGRAPHY VARIANCE    J=',I6,' :')
  912 FORMAT (1X,1P,12E11.1)
  993 FORMAT (1H0,'N=993','I,J,N,LN,MU,MUB=',4I4,2P,2F10.4)
  994 FORMAT (1H0,'N=994  I,J=',2I4,'   L    U  V   MUB(1)  MUB(2)...',
     *  14X,'PLE',6X,'BVF'/(21X,I4,2F8.1,2P,6F11.4,0P,F11.3,F11.4))
  995 FORMAT (' DUT TOO BIG: I,J,L,N,UL,CN,DWTT,MUR,MU,DP=',4I4,3F8.1,
     *  2P,2F8.4,0P,F10.4)
  996 FORMAT (1H0,'MTN WAVE GEN.: I,J,L1,U0,M,U(L1),P0=',
     *  3I4,F7.1,2P,F8.4,0P,F7.1,F11.1)
  997 FORMAT (1X,'SHR WAVE GEN.: I,J,LUP,LDN,DU2K,M,BVF=',
     *  4I4,F7.1,2P,F8.4,1P,2E11.1)
  998 FORMAT (1X,'MC GEN.: I,J,L0,1,HT,W3,M3,C3,US,VS=',
     *  4I4,2X,-3P,F7.2,0P,F7.3,1P,E12.2,1P,F7.2,2X,2F7.2)
  920 FORMAT (1X)
      END SUBROUTINE GWDRAG

      SUBROUTINE DEFORM (P,U,V)
C****
C**** Deformation terms  DEFRM1=du/dx-dv/dy   DEFRM2=du/dy+dv/dx
C**** For spherical coordinates, we are treating DEFRM1 like DIV
C**** and DEFRM2 like CURL (i.e., like FLUX and CIRCULATION),
C**** except the "V" signs are switched.  DEFRM is RMS on u,v grid
C**** MUST be called while PIT,SD,PU,PV are still in WORK1 !!!
C****
      INCLUDE 'BB396M23.COM'
      CHARACTER*80 TITLE
      COMMON/WORK1/PIT(IM,JM),SD(IM,JM,LM-1),PU(IM,JM,LM),PV(IM,JM,LM)
      COMMON/WORK3/PK(IM,JM,LM),UBAR(IM,JM,LM),WSQ(IM,JM,LM),
     *                UDXS(IM),DUMS1(IM),DUMS2(IM),DUMN1(IM),DUMN2(IM)
      COMMON/DRGCOM/AIRX(IM,JM),LMC(IM,JM,2),DEFRM(IM,JM),
     *   ZVART(IM,JM),ZVARX(IM,JM),ZVARY(IM,JM),ZWT(IM,JM),
     *   PDEF,LDEF,LDEFM
      DIMENSION DEFRM1(IM,JM),DEFRM2(IM,JM),DEF1A(IM,JM),DEF2A(IM,JM)
      DIMENSION BYDXYV(JM)
      EQUIVALENCE (WSQ(1,1,1),DEFRM1), (WSQ(1,1,2),DEFRM2)
      EQUIVALENCE (WSQ(1,1,3),DEF1A), (WSQ(1,1,4),DEF2A)
      DATA IFIRST/1/
      SAVE IFIRST,BYDXYV
      IF (IFIRST.EQ.1) THEN
         DO 10 L=1,LM
         PL=(PSF-PTOP)*SIG(L)+PTOP
         IF (PL.GE.700) THEN
           LDEF=L
           PDEF=PL
         ENDIF
         IF (PL.GE.200.) LDEFM=L
   10    CONTINUE
         WRITE (*,*) ' LEVEL FOR DEFORMATION IS: LDEF,PDEF= ',LDEF,PDEF
     *             ,' LDEFM=',LDEFM
         DO J=2,JM
           BYDXYV(J)=1./DXYV(J)
         END DO
         IFIRST=0
      ENDIF
C****
C**** Deformation terms  DEFRM1=du/dx-dv/dy   DEFRM2=du/dy+dv/dx
C**** For spherical coordinates, we are treating DEFRM1 like DIV
C**** and DEFRM2 like CURL (i.e., like FLUX and CIRCULATION),
C**** except the "V" signs are switched
C****
      L=LDEF
C**** U-terms
      DO 91 J=1,JM
      IM1=IM
      DO 91 I=1,IM
      DEFRM1(I,J) = PU(I,J,L)-PU(IM1,J,L)
   91 IM1=I
      DO 92 I=1,IM
   92 UDXS(I)=0.
      IM1=IM
      DO 93 J=1,JM-1
      DO 93 I=1,IM
      UDXN=.5*(U(IM1,J+1,L)+U(I,J+1,L))*DXV(J+1)
      DEFRM2(I,J ) = UDXN-UDXS(I)
      UDXS(I)=UDXN
   93 IM1=I
      DO 94 I=1,IM
   94 DEFRM2(I,JM) =     -UDXS(I)
C**** V-terms
      IM1=IM
      DO 98 I=1,IM
      DO 96 J=1,JM-1
   96 DEFRM1(I,J ) = DEFRM1(I,J ) + (PV(I,J+1,L)-PV(I,J ,L))
      DEFRM1(I,1 ) = DEFRM1(I,1 ) + (PV(I,2  ,L)           )
      DEFRM1(I,JM) = DEFRM1(I,JM) + (           -PV(I,JM,L))
      DO 97 J=2,JM-1
   97 DEFRM2(I,J ) = DEFRM2(I,J ) +
     *  .5*((V(I,J,L)+V(I,J+1,L))-(V(IM1,J,L)+V(IM1,J+1,L)))*DYP(J)
      DEFRM2(I,1 ) = DEFRM2(I,1 ) + (V(I,2 ,L)-V(IM1,2 ,L))*DYP(1)
      DEFRM2(I,JM) = DEFRM2(I,JM) + (V(I,JM,L)-V(IM1,JM,L))*DYP(JM)
   98 IM1=I
C**** Convert to UV-grid
      I=IM
      DO 99 IP1=1,IM
      DUMS1(I)=DEFRM1(I,1)+DEFRM1(IP1,1)
      DUMS2(I)=DEFRM2(I,1)+DEFRM2(IP1,1)
   99 I=IP1
      DO 110 J=2,JM
      I=IM
      DO 100 IP1=1,IM
      DUMN1(I)=DEFRM1(I,J)+DEFRM1(IP1,J)
      DUMN2(I)=DEFRM2(I,J)+DEFRM2(IP1,J)
  100 I=IP1
      DO 105 I=1,IM
      DEFRM1(I,J)=DUMS1(I)+DUMN1(I)
  105 DEFRM2(I,J)=DUMS2(I)+DUMN2(I)
      DO 110 I=1,IM
      DUMS1(I)=DUMN1(I)
      DUMS2(I)=DUMN2(I)
  110 CONTINUE
      DO 120 J=2,JM
      I=IM
      DO 120 IP1=1,IM
      DEFRM1(I,J)=DEFRM1(I,J)/
     *   (DXYV(J)*DSIG(L)*(P(I,J)+P(IP1,J)+P(I,J-1)+P(IP1,J-1)))
      DEFRM2(I,J)=.25*DEFRM2(I,J)*BYDXYV(J)
      DEFRM(I,J)=SQRT(DEFRM1(I,J)**2+DEFRM2(I,J)**2)
  120 I=IP1
C**** Set deformation to zero near the poles
      DO 130 I=1,IM
      DEFRM(I,2)=0.
      DEFRM(I,3)=0.
      DEFRM(I,JM-1)=0.
  130 DEFRM(I,JM  )=0.
C****
C**** Print Deformation terms
C****
CWF   IF (MOD(TAU,8.D0).EQ.0.)  RETURN
CW    D1MAX=0.
CW    D2MAX=0.
CWF   DFMAX=0.
CWF   DO 220 J=4,JM-2
CWF   DO 220 I=1,IM
CW    DEF1A(I,J)=ABS(DEFRM1(I,J))
CW    DEF2A(I,J)=ABS(DEFRM2(I,J))
CW    IF(D1MAX.LT.DEF1A(I,J)) D1MAX=DEF1A(I,J)
CW    IF(D2MAX.LT.DEF2A(I,J)) D2MAX=DEF2A(I,J)
CWF   IF(DFMAX.LT.DEFRM(I,J)) DFMAX=DEFRM(I,J)
  220 CONTINUE
CWF   IHR=NINT(TAU)
CW    SCALE=100./10.**(INT(LOG10(D1MAX)))
CW    TITLE='ABS DEFORMATION TERM 1 = DU/DX - DV/DY (XXXXXXXX S-1)'
CW    WRITE (TITLE(41:48),'(1PE8.0)') SCALE
CCW   CALL MAP0(IM,JM,IHR,TITLE,DEF1A,SCALE,0., 0)
CW    SCALE=100./10.**(INT(LOG10(D2MAX)))
CW    TITLE='ABS DEFORMATION TERM 2 = DU/DY - DV/DX (XXXXXXXX S-1)'
CW    WRITE (TITLE(41:48),'(1PE8.0)') SCALE
CCW   CALL MAP0(IM,JM,IHR,TITLE,DEF2A,SCALE,0., 0)
CWF   SCALE=100./10.**(INT(LOG10(DFMAX)))
CWF   TITLE='RMS OF DEFORMATION TERMS 1 AND 2     (XXXXXXXX S-1)'
CWF   WRITE (TITLE(38:45),'(1PE8.0)') SCALE
CWF   CALL MAP0(IM,JM,IHR,TITLE,DEFRM,SCALE,0., 0)
      RETURN
      END SUBROUTINE DEFORM

      SUBROUTINE io_strat(kunit,iaction,ioerr)
!@sum  io_strat reads and writes strat. model variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "STRAT01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,AIRX,LMC
      CASE (IOREAD:)          ! input from restart file
        READ (kunit,err=10) HEADER,AIRX,LMC
        IF (HEADER.ne.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_strat
      
