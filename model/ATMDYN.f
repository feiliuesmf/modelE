      SUBROUTINE AFLUX (U,V,PIJL)
!@sum  AFLUX Calculates horizontal/vertical air mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,psfmpt,dsig,bydsig,byim
      USE GEOM
      USE DYNAMICS, only : pit,sd,conv,pu,pv,sd_clouds,spa
      IMPLICIT NONE
C**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
      REAL*8 U(IM,JM,LM),V(IM,JM,LM),P(IM,JM) ! p is just workspace
      REAL*8, DIMENSION(IM) :: DUMMYS,DUMMYN
      REAL*8 PIJL(IM,JM,LM)
      INTEGER I,J,L,IP1,IM1
      REAL*8 PUS,PUN,PVS,PVN,PBS,PBN,SDNP,SDSP
C****
C**** BEGINNING OF LAYER LOOP
C****
      DO 2000 L=1,LM
C****
C**** COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
C**** ARAKAWA'S SCHEME B             PV   U,V    SECONDARY GRID ROW
C****
C**** COMPUTE PU, THE WEST-EAST MASS FLUX, AT NON-POLAR POINTS
      DO 2154 J=2,JM-1
      DO 2154 I=1,IM
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
      CALL AVRX (SPA(1,1,L))
      I=IM
      DO 2166 J=2,JM-1
      DO 2165 IP1=1,IM
      PU(I,J,L)=.25*DYP(J)*SPA(I,J,L)*(PIJL(I,J,L)+PIJL(IP1,J,L))*
     *  DSIG(L)
 2165 I=IP1
 2166 CONTINUE
C**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX
      IM1=IM
      DO 2172 J=2,JM
      DO 2170 I=1,IM
      PV(I,J,L)=.25*DXV(J)*(V(I,J,L)+V(IM1,J,L))*
     *   (PIJL(I,J,L)+PIJL(I,J-1,L))*DSIG(L)
 2170 IM1=I
 2172 CONTINUE
C**** COMPUTE PU*3 AT THE POLES
      PUS=0.
      PUN=0.
      PVS=0.
      PVN=0.
      DO 1110 I=1,IM
      PUS=PUS+U(I,2,L)
      PUN=PUN+U(I,JM,L)
      PVS=PVS+PV(I,2,L)
 1110 PVN=PVN+PV(I,JM,L)
      PUS=.25*DYP(2)*PUS*PIJL(1,1,L)*BYIM
      PUN=.25*DYP(JM-1)*PUN*PIJL(1,JM,L)*BYIM
      PVS=PVS*BYIM
      PVN=PVN*BYIM
      DUMMYS(1)=0.
      DUMMYN(1)=0.
      DO 1120 I=2,IM
      DUMMYS(I)=DUMMYS(I-1)+(PV(I,2,L) -PVS)*BYDSIG(L)
 1120 DUMMYN(I)=DUMMYN(I-1)+(PV(I,JM,L)-PVN)*BYDSIG(L)
      PBS=0.
      PBN=0.
      DO 1130 I=1,IM
      PBS=PBS+DUMMYS(I)
 1130 PBN=PBN+DUMMYN(I)
      PBS=PBS*BYIM
      PBN=PBN*BYIM
      DO 1140 I=1,IM
      SPA(I,1,L)=4.*(PBS-DUMMYS(I)+PUS)/(DYP(2)*PIJL(1,1,L))
      SPA(I,JM,L)=4.*(DUMMYN(I)-PBN+PUN)/(DYP(JM-1)*PIJL(1,JM,L))
      PU(I,1,L)=3.*(PBS-DUMMYS(I)+PUS)*DSIG(L)
 1140 PU(I,JM,L)=3.*(DUMMYN(I)-PBN+PUN)*DSIG(L)
C****
C**** CONTINUITY EQUATION
C****
C**** COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
      DO 1510 J=2,JM-1
      IM1=IM
      DO 1510 I=1,IM
      CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
 1510 IM1=I
      CONV(1,1,L)=-PVS
      CONV(1,JM,L)=PVN
 2000 CONTINUE
C****
C**** END OF HORIZONTAL ADVECTION LAYER LOOP
C****
C**** COMPUTE PIT, THE PRESSURE TENDENCY
C     PIT(I,J)=CONV(I,J,1)
      DO 2420 L=LM,2,-1
      PIT(1,1)=PIT(1,1)+CONV(1,1,L)
      PIT(1,JM)=PIT(1,JM)+CONV(1,JM,L)
      DO 2420 J=2,JM-1
      DO 2420 I=1,IM
 2420 PIT(I,J)=PIT(I,J)+CONV(I,J,L)
C**** COMPUTE SD, SIGMA DOT
      SD(1, 1,LM-1)=CONV(1, 1,LM)
      SD(1,JM,LM-1)=CONV(1,JM,LM)
      DO 2430 J=2,JM-1
      DO 2430 I=1,IM
 2430 SD(I,J,LM-1)=CONV(I,J,LM)
      DO 2435 L=LM-2,LS1-1,-1
      SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)
      SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)
      DO 2435 J=2,JM-1
      DO 2435 I=1,IM
      SD(I, J,L)=SD(I, J,L+1)+CONV(I, J,L+1)
 2435 CONTINUE
      DO 2440 L=LS1-2,1,-1
      SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)-DSIG(L+1)*PIT(1, 1)
      SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)-DSIG(L+1)*PIT(1,JM)
      DO 2440 J=2,JM-1
      DO 2440 I=1,IM
      SD(I, J,L)=SD(I, J,L+1)+CONV(I, J,L+1)-DSIG(L+1)*PIT(I, J)
 2440 CONTINUE
      DO 2450 L=1,LM-1
      DO 2450 I=2,IM
      SD(I,1,L)=SD(1,1,L)
 2450 SD(I,JM,L)=SD(1,JM,L)
C**** temporary fix for CLOUDS module
      SD_CLOUDS(:,:,1)    = PIT
      SD_CLOUDS(:,:,2:LM) = SD(:,:,1:LM-1)
C****
      RETURN
      END SUBROUTINE AFLUX

      SUBROUTINE ADVECM (P,PA,DT1)
!@sum  ADVECM Calculates updated column pressures using mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,mrch,zatmo,u,v,t,q
      USE GEOM, only : bydxyp,imaxj
      USE DYNAMICS, only : pit
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: P(IM,JM)
      REAL*8, INTENT(OUT) :: PA(IM,JM)
      REAL*8, INTENT(IN) :: DT1
      INTEGER I,J,L,K,IMAX  !@var I,J,L,K  loop variables

C**** COMPUTE PA, THE NEW SURFACE PRESSURE
      DO J=1,JM
        DO I=1,IMAXJ(J)
          PA(I,J)=P(I,J)+(DT1*PIT(I,J)*BYDXYP(J))
          IF (PA(I,J).GT.1150.) WRITE (6,990) I,J,MRCH,P(I,J),PA(I,J),
     *         ZATMO(I,J),(U(I-1,J,L),U(I,J,L),U(I-1,J+1,L),U(I,J+1,L),
     *         V(I-1,J,L),V(I,J,L),V(I-1,J+1,L),V(I,J+1,L),
     *         T(I,J,L),Q(I,J,L),L=1,LM)
        END DO
      END DO
      PA(2:IM, 1)=PA(1,1)
      PA(2:IM,JM)=PA(1,JM)
C****
      RETURN
  990 FORMAT (/'0PRESSURE DIAGNOSTIC     I,J,MRCH,P,PA=',3I4,2F10.2/
     *  '     ZATMO=',F10.3/
     *  '0    U(I-1,J)     U(I,J)   U(I-1,J+1)    U(I,J+1)    V(I-1,J)',
     *   '     V(I,J)   V(I-1,J+1)    V(I,J+1)     T(I,J)     Q(I,J)'/
     *  (1X,9F12.3,F12.6))
      END SUBROUTINE ADVECM

      SUBROUTINE PGF (UT,VT,PB,U,V,T,SZ,P,DT1)
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,sige,ptop
     *     ,zatmo,sig,modd5k,bydsig
      USE GEOM, only : imaxj,dxyv,dxv,dyv,dxyp,dyp,dxp
      USE DYNAMICS, only : gz,pu,pit,phi,spa,dut,dvt
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM,LM) :: U,V,T
      REAL*8, DIMENSION(IM,JM) :: FD,RFDUX
      REAL*8 UT(IM,JM,LM),VT(IM,JM,LM),TT(IM,JM,LM),
     *  PA(IM,JM),PB(IM,JM),QT(IM,JM,LM),P(IM,JM,LM)
      REAL*8, DIMENSION(IM,JM,1) :: PU0

      REAL*8 PKE(LM+1)
      REAL*8 SZ(IM,JM,LM),DT4,DT1
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDU,PHIDN,FACTOR
      INTEGER I,J,L,IM1,IP1,IMAX  !@var I,J,IP1,IM1,L,IMAX loop variab.
C****
      DT4=DT1/4.
      DO 10 L=1,LM+1
   10 PKE(L)=(SIGE(L)*PSFMPT+PTOP)**KAPA
C****
C**** VERTICAL DIFFERENCING
C****
      DO 100 L=LS1,LM
      DO 100 J=1,JM
      DO 100 I=1,IM
  100 SPA(I,J,L)=0.
      DO 330 J=1,JM
      DO 330 I=1,IMAXJ(J)
      PIJ=P(I,J,1)
      PDN=PIJ+PTOP
      PKDN=PDN**KAPA
      PHIDN=ZATMO(I,J)
C**** LOOP OVER THE LAYERS
      DO 310 L=1,LM
      IF(L.LE.LS1-1) GO TO 290
      PKPDN=PKDN*PDN
      PKPPDN=PKPDN*PDN
      DP=DSIG(L)*PSFMPT
      BYDP=1./DP
      P0=SIG(L)*PSFMPT+PTOP
      TZBYDP=2.*SZ(I,J,L)*BYDP
      X=T(I,J,L)+TZBYDP*P0
      PUP=SIGE(L+1)*PSFMPT+PTOP
      PKUP=PKE(L+1)
      PKPUP=PKUP*PUP
      PKPPUP=PKPUP*PUP
      GO TO 300
  290 PKPDN=PKDN*PDN
      PKPPDN=PKPDN*PDN
      DP=DSIG(L)*PIJ
      BYDP=1./DP
      P0=SIG(L)*PIJ+PTOP
      TZBYDP=2.*SZ(I,J,L)*BYDP
      X=T(I,J,L)+TZBYDP*P0
      PUP=SIGE(L+1)*PIJ+PTOP
      PKUP=PUP**KAPA
      PKPUP=PKUP*PUP
      PKPPUP=PKPUP*PUP
C**** CALCULATE SPA, MASS WEIGHTED THROUGHOUT THE LAYER
      SPA(I,J,L)=RGAS*((X+TZBYDP*PTOP)*(PKPDN-PKPUP)*BYKAPAP1
     *     -X*PTOP*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *     *BYDP
C**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
  300 PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1
     *  -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *  *BYDP*BYKAPAP1)
C**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
      PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP)
     *     *BYKAPAP1)
      PDN=PUP
  310 PKDN=PKUP
  330 CONTINUE
C**** SET POLAR VALUES FROM THOSE AT I=1
      DO 340 L=1,LM
      DO 340 I=2,IM
      SPA(I,1,L)=SPA(1,1,L)
      SPA(I,JM,L)=SPA(1,JM,L)
      PHI(I,1,L)=PHI(1,1,L)
  340 PHI(I,JM,L)=PHI(1,JM,L)
         DO 3081 L=1,LM
         DO 3081 J=1,JM
         DO 3081 I=1,IM
 3081    GZ(I,J,L)=PHI(I,J,L)
C****
C**** PRESSURE GRADIENT FORCE
C****
C**** NORTH-SOUTH DERIVATIVE AFFECTS THE V-COMPONENT OF MOMENTUM
      DO 3236 L=1,LM
      DO 3236 J=2,JM
      FACTOR = DT4*DXV(J)*DSIG(L)
      IM1=IM
      DO 3234 I=1,IM
      FLUX=    ((P(I,J,L)+P(I,J-1,L))*(PHI(I,J,L)-PHI(I,J-1,L))+
     *  (SPA(I,J,L)+SPA(I,J-1,L))*(P(I,J,L)-P(I,J-1,L)))*FACTOR
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3234 IM1=I
 3236 CONTINUE
C**** SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
      DO 3300 L=1,LM
      DO 3275 I=1,IM
      PU(I,1,L)=0.
 3275 PU(I,JM,L)=0.
      I=IM
      DO 3290 J=2,JM-1
      DO 3280 IP1=1,IM
      PU(I,J,L)=(P(IP1,J,L)+P(I,J,L))*(PHI(IP1,J,L)-PHI(I,J,L))+
     *  (SPA(IP1,J,L)+SPA(I,J,L))*(P(IP1,J,L)-P(I,J,L))
 3280 I=IP1
 3290 CONTINUE
      CALL AVRX (PU(1,1,L))
      DO 3294 J=2,JM
      FACTOR = -DT4*DYV(J)*DSIG(L)
      DO 3294 I=1,IM
 3294 DUT(I,J,L)=DUT(I,J,L)+FACTOR*(PU(I,J,L)+PU(I,J-1,L))
 3300 CONTINUE

C**** CALL DIAGNOSTICS
      IF(MRCH.GT.0) THEN
         IF(MODD5K.LT.MRCH) CALL DIAG5D (6,MRCH,DUT,DVT)
         IF(MODD5K.LT.MRCH) CALL DIAG9D (3,DT1,U,V,DUT,DVT,PIT)
      ENDIF
C****
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
C****
      DO 3410 J=2,JM-1
      DO 3410 I=1,IM
 3410 FD(I,J)=PB(I,J)*DXYP(J)
      FDSP=PB(1, 1)*DXYP( 1)
      FDNP=PB(1,JM)*DXYP(JM)
      FDSP=FDSP+FDSP
      FDNP=FDNP+FDNP
      DO 3520 I=1,IM
      FD(I, 1)=FDSP
 3520 FD(I,JM)=FDNP

      DO 3530 J=2,JM
      I=IM
      DO 3525 IP1=1,IM
      RFDUX(I,J)=4./(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
 3525 I = IP1
 3530 CONTINUE
      DO 3550 L=1,LM
      DO 3550 J=2,JM
      RFDU=1./(PSFMPT*DXYV(J)*DSIG(L))
      DO 3540 I=1,IM
      IF(L.LT.LS1) RFDU=RFDUX(I,J)*BYDSIG(L)
      VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)*RFDU
      UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)*RFDU
 3540 CONTINUE
 3550 CONTINUE
      RETURN
      END SUBROUTINE PGF

      SUBROUTINE AVRX (X)
!@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,imh
      USE GEOM, only : dlon,dxp,dyp
      USE DYNAMICS, only : xAVRX
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: X(IM,JM)
      REAL*8, SAVE :: SM(IMH,JM),DRAT(JM)
      REAL*8, SAVE, DIMENSION(IMH) :: BYSN
      REAL*8, DIMENSION(0:IMH) :: AN,BN
      INTEGER, SAVE :: NMIN(JM)
      INTEGER, SAVE :: IFIRST = 1
      INTEGER J,N

      IF (IFIRST.EQ.1) THEN
      IFIRST=0
C     CALL FFT0(IM)
      DO N=1,IMH
        BYSN(N)=xAVRX/SIN(.5*DLON*N)
      END DO
      DO 50 J=2,JM-1
        DRAT(J) = DXP(J)/DYP(3)
        DO 40 N=IMH,1,-1
          SM(N,J) = BYSN(N)*DRAT(J)
          IF(SM(N,J).GT.1.) THEN
            NMIN(J) = N+1
            GO TO 50
          ENDIF
 40     CONTINUE
 50   CONTINUE
      END IF
C****
      DO 140 J=2,JM-1
      IF (DRAT(J).GT.1) GO TO 140
      CALL FFT (X(1,J),AN,BN)
      DO N=NMIN(J),IMH-1
        AN(N)=SM(N,J)*AN(N)
        BN(N)=SM(N,J)*BN(N)
      END DO
      AN(IMH) = SM(IMH,J)*AN(IMH)
      CALL FFTI(AN,BN,X(1,J))
  140 CONTINUE
      RETURN
      END SUBROUTINE AVRX

      SUBROUTINE FILTER
!@sum  FILTER Performs 8-th order shapiro filter in zonal direction
!@auth Original development team
!@ver  1.0
C****
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER
C****        3  SMOOTH P AND T
C****
      USE CONSTANT, only : bbyg,gbyrb,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,t,p,q,wm,mfiltr,zatmo,ptop
     *     ,byim,sig
      USE SOMTQ_COM, only : tmom,qmom
      USE PBLCOM, only : tsavg
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: X,Y
      REAL*8 PSUMO(JM)

      REAL*8 POLD(IM,JM),PRAT(IM,JM)
      REAL*8 PSUMN,PDIF,AKAP
      INTEGER I,J,L  !@var I,J,L  loop variables

      IF (MOD(MFILTR,2).NE.1) GO TO 200
C****
C**** SEA LEVEL PRESSURE FILTER ON P
C****
      DO 120 J=2,JM-1
         PSUMO(J)=0.
      DO 120 I=1,IM
         PSUMO(J)=PSUMO(J)+P(I,J)
         POLD(I,J)=P(I,J)      ! Save old pressure
      Y(I,J)=(1.+BBYG*ZATMO(I,J)/TSAVG(I,J))**GBYRB
  120 X(I,J)=(P(I,J)+PTOP)*Y(I,J)
      CALL SHAP1D (8,X)
      DO 150 J=2,JM-1
         PSUMN=0.
      DO 140 I=1,IM
      P(I,J)=X(I,J)/Y(I,J)-PTOP
  140 PSUMN=PSUMN+P(I,J)
         PDIF=(PSUMN-PSUMO(J))*BYIM
      DO 145 I=1,IM
  145    P(I,J)=P(I,J)-PDIF
  150 CONTINUE
C**** Scale mixing ratios (incl moments) to conserve mass/heat
      DO J=2,JM-1
        DO I=1,IM
          PRAT(I,J)=POLD(I,J)/P(I,J)
        END DO
      END DO
      DO L=1,LS1-1
        DO J=2,JM-1
          DO I=1,IM
             Q(I,J,L)=  Q(I,J,L)*PRAT(I,J)
             QMOM(:,I,J,L)=QMOM(:,I,J,L)*PRAT(I,J)
             T(I,J,L)=  T(I,J,L)*PRAT(I,J)
             TMOM(:,I,J,L)=TMOM(:,I,J,L)*PRAT(I,J)
             WM(I,J,L)= WM(I,J,L)*PRAT(I,J)
          END DO
        END DO
      END DO

      CALL CALC_AMPK(LS1-1)

  200 IF (MFILTR.LT.2) RETURN
C****
C**** TEMPERATURE STRATIFICATION FILTER ON T
C****
      AKAP=KAPA-.205
      DO 260 L=1,LS1-1
      DO 220 J=2,JM-1
      DO 220 I=1,IM
      Y(I,J)=(SIG(L)*P(I,J)+PTOP)**AKAP
  220 X(I,J)=T(I,J,L)*Y(I,J)
      CALL SHAP1D (8,X)
      DO 240 J=2,JM-1
      DO 240 I=1,IM
  240 T(I,J,L)=X(I,J)/Y(I,J)
  260 CONTINUE
      DO 280 L=LS1,LM
      DO 270 J=2,JM-1
      DO 270 I=1,IM
  270 X(I,J)=T(I,J,L)
      CALL SHAP1D (8,X)
      DO 280 J=2,JM-1
      DO 280 I=1,IM
  280 T(I,J,L)=X(I,J)
      RETURN
      END SUBROUTINE FILTER

      SUBROUTINE FLTRUV(U,V)
!@sum  FLTRUV Filters 2 gridpoint noise from the velocity fields
!@auth Original development team
!@ver  1.0
!@calls SHAP1D
      USE MODEL_COM, only : im,jm,lm,byim
C**********************************************************************
C**** FILTERING IS DONE IN BOTH DIMENSIONS WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES.
C**********************************************************************
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,LM), INTENT(INOUT) :: U,V
      REAL*8 X(IM),Y(0:JM+1),F2D(IM,JM)
      LOGICAL*4, SAVE :: QFILY,QFILX,QFIL2D
      INTEGER, SAVE :: NSHAP,ISIGN
      REAL*8, SAVE :: by4toN
      INTEGER I,J,L,N  !@var I,J,L,N  loop variables
      REAL*8 Y4TO8,YJ,YJM1,X1,XI,XIM1, BY16
      INTEGER,SAVE :: IFIRST = 1
      PARAMETER (BY16=1./16.)
C****
      IF(IFIRST.EQ.1) THEN
         IFIRST = 0
C        CALL FFT0(IM)
         QFILY  = .FALSE.
         QFILX  = .TRUE.
         QFIL2D = .FALSE.
         NSHAP  = 8
         ISIGN  = (-1.0)**(NSHAP-1)
         by4toN  = 1./(4.**NSHAP)
      ENDIF
C****
C**** Filtering in east-west direction
C****
      DO 350 L=1,LM
C**** Filter U component of momentum
      IF(QFIL2D) GOTO 250
      DO 240 J=2,JM
      DO 210 I=1,IM
  210 X(I) = U(I,J,L)
      DO 230 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 220 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  220 XIM1 = XI
  230 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 240 I=1,IM
  240 U(I,J,L) = U(I,J,L) + ISIGN*X(I)*by4toN
      GOTO 270
  250 DO 255,J=2,JM-1
      F2D(1,J)=.125*(U(IM,J,L)+U(1,J+1,L)+U(2,J,L)+
     *                  U(1,J-1,L)-4.*U(1,J,L))+
     *        BY16*(U(IM,J-1,L)+U(IM,J+1,L)+U(2,J+1,L)+
     *              U(2,J-1,L)-4.*U(1,J,L))
      DO 260,I=2,IM-1
      F2D(I,J)=.125*(U(I-1,J,L)+U(I,J+1,L)+U(I+1,J,L)+
     *                  U(I,J-1,L)-4.*U(I,J,L))+
     *        BY16*(U(I-1,J-1,L)+U(I-1,J+1,L)+U(I+1,J+1,L)+
     *              U(I+1,J-1,L)-4.*U(I,J,L))
  260 CONTINUE
      F2D(IM,J)=.125*(U(IM-1,J,L)+U(IM,J+1,L)+U(1,J,L)+
     *                  U(IM,J-1,L)-4.*U(IM,J,L))+
     *        BY16*(U(IM-1,J-1,L)+U(IM-1,J+1,L)+U(1,J+1,L)+
     *              U(1,J-1,L)-4.*U(IM,J,L))
  255 CONTINUE
      DO 265,J=2,JM-1
      DO 265,I=1,IM
        U(I,J,L)=U(I,J,L)+F2D(I,J)
  265 CONTINUE
  270 CONTINUE
C**** Filter V component of momentum
      DO 340 J=2,JM
      DO 310 I=1,IM
  310 X(I) = V(I,J,L)
      DO 330 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 320 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  320 XIM1 = XI
  330 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 340 I=1,IM
  340 V(I,J,L) = V(I,J,L) + ISIGN*X(I)*by4toN
  350 CONTINUE
C****
C**** Filtering in north-south direction
C****
      IF(.NOT.QFILY) GOTO 651
  400 Y4TO8 = 1./(4.**NSHAP)
C**** Filter U component of momentum
      DO 650 L=1,LM
      DO 540 I=1,IM
      DO 510 J=1,JM
  510 Y(J) = U(I,J,L)
      DO 530 N=1,NSHAP
         Y(0)    = Y(2)
         Y(JM+1) = Y(JM-1)
         YJM1 = Y(0)
         DO 520 J=1,JM
            YJ   = Y(J)
            Y(J) = YJM1-YJ-YJ+Y(J+1)
  520    YJM1 = YJ
  530 CONTINUE
      DO 540 J=1,JM
  540 U(I,J,L) = U(I,J,L) + ISIGN*Y(J)*Y4TO8
C**** Make U winds at poles to be uniform
      DO 560 J=1,JM,JM-1
      YJ = 0.
      DO 550 I=1,IM
  550 YJ = YJ + U(I,J,L)
      DO 560 I=1,IM
  560 U(I,J,L) = YJ*BYIM
C**** Filter V component of momentum
      DO 640 I=1,IM
      DO 610 J=1,JM-1
  610 Y(J) = V(I,J,L)
      DO 630 N=1,NSHAP
      YJM1 = Y(1)
      Y(1) = Y(2)-Y(1)
      DO 620 J=2,JM-2
      YJ   = Y(J)
      Y(J) = YJM1-YJ-YJ+Y(J+1)
  620 YJM1 = YJ
  630 Y(JM-1)= YJM1-Y(JM-1)
      DO 640 J=1,JM-1
  640 V(I,J,L) = V(I,J,L) + ISIGN*Y(J)*Y4TO8
  650 CONTINUE
C****
  651 RETURN
      END SUBROUTINE FLTRUV

      SUBROUTINE SHAP1D (NORDER,X)
!@sum  SHAP1D Smoothes in zonal direction use n-th order shapiro filter
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only :im,jm
      IMPLICIT NONE
!@var NORDER order of shapiro filter (must be even)
      INTEGER, INTENT(IN) :: NORDER
      REAL*8, INTENT(INOUT), DIMENSION(IM,JM) :: X
      REAL*8, DIMENSION(IM)::XS
      REAL*8 by4toN,XS1,XSIM1,XSI
      INTEGER I,J,N   !@var I,J,N  loop variables

      by4toN=1./4.**NORDER
      DO J=2,JM-1
        XS(:) = X(:,J)
        DO N=1,NORDER
          XS1=XS(1)
          XSIM1=XS(IM)
          DO I=1,IM-1
            XSI=XS(I)
            XS(I)=XSIM1-XSI-XSI+XS(I+1)
            XSIM1=XSI
          END DO
          XS(IM)=XSIM1-XS(IM)-XS(IM)+XS1
        END DO
        X(:,J)=X(:,J)-XS(:)*by4toN
      END DO
      RETURN
      END SUBROUTINE SHAP1D

      SUBROUTINE CALC_PIJL(lmax,p,pijl)
!@sum  CALC_PIJL Fills in P as 3-D
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,psfmpt
C****
      implicit none
      double precision, dimension(im,jm) :: p
      double precision, dimension(im,jm,lm) :: pijl
      integer :: l,lmax

      do l=1,ls1-1
        pijl(:,:,l) = p(:,:)
      enddo
      do l=ls1,lmax
        pijl(:,:,l) = PSFMPT
      enddo
      return
      end subroutine calc_pijl


      SUBROUTINE CALC_AMPK(LMAX)
!@sum  CALC_AMPK calculate air mass and pressure functions
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,p,dsig,sig,sige,ptop,psfmpt
      USE DYNAMICS, only : plij,pdsig,pmid,pk,pedn,pek,sqrtp
      IMPLICIT NONE

      INTEGER :: I,J,IMAX,L  !@var I,J,IMAX,L  loop variables
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max. level for update

C**** Calculate air mass, layer pressures, P**K, and sqrt(P)
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure. Routine should be called with LMAX=LM at start, and
C**** subsequentaly with LMAX=LS1-1

C**** Fill in polar boxes
      P(2:IM,1) = P(1,1)
      P(2:IM,JM)= P(1,JM)

      DO J=1,JM
         DO I=1,IM
            DO L=1,LS1-1
               PLIJ(L,I,J) = P(I,J)
               PDSIG(L,I,J) = P(I,J)*DSIG(L)
               PMID(L,I,J) = SIG(L)*P(I,J)+PTOP
               PK  (L,I,J) = PMID(L,I,J)**KAPA
               PEDN(L,I,J) = SIGE(L)*P(I,J)+PTOP
               PEK (L,I,J) = PEDN(L,I,J)**KAPA
c               AM  (L,I,J) = P(I,J)*DSIG(L)*1d2*BYGRAV
c               BYAM(L,I,J) = 1./AM(L,I,J)
            END DO
            DO L=LS1,LMAX
               PLIJ(L,I,J) = PSFMPT
               PDSIG(L,I,J) = PSFMPT*DSIG(L)
               PMID(L,I,J) = SIG(L)*PSFMPT+PTOP
               PK  (L,I,J) = PMID(L,I,J)**KAPA
               PEDN(L,I,J) = SIGE(L)*PSFMPT+PTOP
               PEK (L,I,J) = PEDN(L,I,J)**KAPA
c               AM  (L,I,J) = PSFMPT*DSIG(L)*1d2*BYGRAV
c               BYAM(L,I,J) = 1./AM(L,I,J)
            END DO
            IF (LMAX.eq.LM) THEN
               PEDN(LM+1,I,J) = SIGE(LM+1)*PSFMPT+PTOP
               PEK (LM+1,I,J) = PEDN(LM+1,I,J)**KAPA
            END IF
            SQRTP(I,J) = SQRT(P(I,J))
         END DO
      END DO

      RETURN
      END SUBROUTINE CALC_AMPK

      SUBROUTINE CALC_AMP(p,amp)
!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
!@auth Jean Lerner/Max Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,dsig,psf,ptop
      USE GEOM, only : dxyp
C****
      implicit none
      double precision, dimension(im,jm) :: p
      double precision, dimension(im,jm,lm) :: amp
      integer :: j,l
      do l=1,ls1-1
      do j=1,jm
        amp(:,j,l) = p(:,j)*dxyp(j)*dsig(l)
      enddo
      enddo
      do l=ls1,lm
      do j=1,jm
        amp(:,j,l) = (psf-ptop)*dxyp(j)*dsig(l)
      enddo
      enddo
      return
C****
      end subroutine calc_amp

      SUBROUTINE DYNAM
!@sum  DYNAM Integrate dynamic terms
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : by3
      USE MODEL_COM, only : im,jm,lm,u,v,t,p,q,wm,dsig,NIdyn,dt,MODD5K
     *     ,NSTEP,NDA5K,ndaa,mrch,psfmpt,ls1,lsdrag
      USE GEOM, only : dyv,dxv
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE DYNAMICS, only : ptold,pu,pv,pit,sd,phi,dut,dvt
      USE DAGCOM, only : aij,ij_fpeu,ij_fpev,ij_fqu,ij_fqv,ij_fmv,ij_fmu
     *     ,ij_fgzu,ij_fgzv
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: PRAT
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT,TT,TZ,TZT,MA
      REAL*8, DIMENSION(IM,JM,LM) :: UX,VX,PIJL
      REAL*8 PA(IM,JM),PB(IM,JM),PC(IM,JM),FPEU(IM,JM),FPEV(IM,JM),
     *          FWVU(IM,JM),FWVV(IM,JM)

      REAL*8 DTFS,DTLF,PP,UU,VV
      INTEGER I,J,L,IP1,IM1   !@var I,J,L,IP1,IM1  loop variables
      INTEGER NS, NSOLD,MODDA    !? ,NIdynO

!?    NIdynO=MOD(NIdyn,2)   ! NIdyn odd is currently not an option
      DTFS=DT*2./3.
      DTLF=2.*DT
      NS=0
      NSOLD=0                            ! strat
      PTOLD(:,:) = P(:,:)
  300 CONTINUE
      UX(:,:,:)  = U(:,:,:)
      UT(:,:,:)  = U(:,:,:)
      VX(:,:,:)  = V(:,:,:)
      VT(:,:,:)  = V(:,:,:)
C*** copy z-moment of temperature into contiguous memory
      tz(:,:,:) = tmom(mz,:,:,:)
      PA(:,:) = P(:,:)
      PB(:,:) = P(:,:)
      PC(:,:) = P(:,:)
C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      MRCH=0
C     CALL DYNAM (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS)
      CALL CALC_PIJL(LM,P,PIJL)
      CALL AFLUX (U,V,PIJL)
      CALL ADVECM (P,PB,DTFS)
      CALL GWDRAG (PB,UX,VX,T,TZ,DTFS)   ! strat
      CALL VDIFF (PB,UX,VX,T,Q,DTFS)     ! strat
      CALL ADVECV (P,UX,VX,PB,U,V,Pijl,DTFS)  !P->pijl
      CALL PGF (UX,VX,PB,U,V,T,TZ,Pijl,DTFS)
      CALL FLTRUV(UX,VX)
C**** INITIAL BACKWARD STEP IS ODD, QT = Q + DT*F(QX)
      MRCH=-1
C     CALL DYNAM (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT)
      CALL CALC_PIJL(LS1-1,PB,PIJL)
      CALL AFLUX (UX,VX,PIJL)
      CALL ADVECM (P,PA,DT)
      CALL GWDRAG (PA,UT,VT,T,TZ,DT)   ! strat
      CALL VDIFF (PA,UT,VT,T,Q,DT)     ! strat
      CALL ADVECV (P,UT,VT,PA,UX,VX,Pijl,DT)   !PB->pijl
      CALL PGF (UT,VT,PA,UX,VX,T,TZ,Pijl,DT)
      CALL FLTRUV(UT,VT)
      GO TO 360
C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)
  340 MRCH=-2
C     CALL DYNAM (UT,VT,TT,PT,QT,U,V,T,P,Q,DTLF)
      CALL CALC_PIJL(LS1-1,P,PIJL)
      CALL AFLUX (U,V,PIJL)
      CALL ADVECM (PA,PB,DTLF)
      CALL GWDRAG (PB,UT,VT,T,TZ,DTLF)   ! strat
      CALL VDIFF (PB,UT,VT,T,Q,DTLF)     ! strat
      CALL ADVECV (PA,UT,VT,PB,U,V,Pijl,DTLF)   !P->pijl
      CALL PGF (UT,VT,PB,U,V,T,TZ,Pijl,DTLF)
      CALL FLTRUV(UT,VT)
C**** LOAD PB TO PA
      PA(:,:) = PB(:,:)
C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)
  360 NS=NS+2
         MODD5K=MOD(NSTEP+NS-NIdyn+NDA5K*NIdyn+2,NDA5K*NIdyn+2)
      MRCH=2
C     CALL DYNAM (U,V,T,P,Q,UT,VT,TT,PT,QT,DTLF)
      CALL CALC_PIJL(LS1-1,PA,PIJL)
      CALL AFLUX (UT,VT,PIJL)
      CALL ADVECM (PC,P,DTLF)
      CALL GWDRAG (P,U,V,T,TZ,DTLF)   ! strat
      CALL ADVECV (PC,U,V,P,UT,VT,Pijl,DTLF)     !PA->pijl
         MODDA=MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2)  ! strat
         IF(MODDA.LT.MRCH) CALL DIAGA0 (T)   ! strat
         FPEU(:,:) = 0.
         FPEV(:,:) = 0.
         FWVU(:,:) = 0.
         FWVV(:,:) = 0.
      TT(:,:,:) = T(:,:,:)
      TZT(:,:,:)= TZ(:,:,:)
      call calc_amp(pc,ma)
      CALL AADVT (MA,T,TMOM, SD,PU,PV, DTLF,.FALSE.,FPEU,FPEV)
C*** copy z-moment of temperature into contiguous memory
      tz(:,:,:) = tmom(mz,:,:,:)
      call calc_amp(pc,ma)
      CALL AADVT (MA,Q,QMOM, SD,PU,PV, DTLF,.TRUE. ,FWVU,FWVV)
      CALL VDIFF (P,U,V,T,Q,DTLF)        ! strat
      PC(:,:)    = .5*(P(:,:)+PC(:,:))
      TT(:,:,:)  = .5*(T(:,:,:)+TT(:,:,:))
      TZT(:,:,:) = .5*(TZ(:,:,:)+TZT(:,:,:))
         AIJ(:,:,IJ_FPEU) = AIJ(:,:,IJ_FPEU)+FPEU(:,:)
         AIJ(:,:,IJ_FPEV) = AIJ(:,:,IJ_FPEV)+FPEV(:,:)
         AIJ(:,:,IJ_FQU)  = AIJ(:,:,IJ_FQU )+FWVU(:,:)
         AIJ(:,:,IJ_FQV)  = AIJ(:,:,IJ_FQV )+FWVV(:,:)
         DO L=1,LM
           AIJ(:,2:JM,IJ_FMV)  = AIJ(:,2:JM,IJ_FMV )+PV(:,2:JM,L)*DTLF
           AIJ(:,1,IJ_FMU)  = AIJ(:, 1,IJ_FMU )+PU(:, 1,L)*DTLF*BY3
           AIJ(:,JM,IJ_FMU) = AIJ(:,JM,IJ_FMU )+PU(:,JM,L)*DTLF*BY3
           AIJ(:,2:JM-1,IJ_FMU)=AIJ(:,2:JM-1,IJ_FMU)+PU(:,2:JM-1,L)*DTLF
         END DO
      CALL CALC_PIJL(LS1-1,PC,PIJL)
      CALL PGF (U,V,P,UT,VT,TT,TZT,Pijl,DTLF)    !PC->pijl
         DO L=1,LM
         DO J=2,JM
         IM1=IM
         DO I=1,IM
           PP=.5*(PHI(I,J-1,L)+PHI(I,J,L))
           UU=.5*(U(I,J,L)+U(IM1,J,L))
           AIJ(I,J,IJ_FGZU)=AIJ(I,J,IJ_FGZU)+PP*UU*DYV(J)*DTLF
           VV=.5*(V(I,J,L)+V(IM1,J,L))
           AIJ(I,J,IJ_FGZV)=AIJ(I,J,IJ_FGZV)+PP*VV*DXV(J)*DTLF
         IM1=I
         END DO
         END DO
         END DO
      CALL FLTRUV(U,V)
C**** LOAD P TO PC
      PC(:,:)=P(:,:)
         IF (MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2).LT.MRCH) THEN
           CALL DIAGA (U,V,T,P,Q,PIT,SD)
           CALL DIAGB (U,V,T,P,Q,WM,DUT,DVT)
           CALL EPFLUX (U,V,T,P)
         ENDIF
C**** Restart after 8 steps due to divergence of solutions
C**** STRATOSPHERIC MOMENTUM DRAG must be called at least once
      IF (NS-NSOLD.LT.8 .AND. NS.LT.NIdyn) GO TO 340
      CALL CALC_AMPK(LS1-1)
      CALL SDRAG (LSDRAG,DT*(NS-NSOLD))
      NSOLD=NS
      IF (NS.LT.NIdyn) GO TO 300
C**** Scale WM mixing ratios to conserve liquid water
      PRAT(:,:)=PTOLD(:,:)/P(:,:)
      DO L=1,LS1-1
        WM(:,:,L)=WM(:,:,L)*PRAT(:,:)
      END DO
      RETURN
      END SUBROUTINE DYNAM

      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@auth Ye Cheng
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE MODEL_COM, only : im,jm,t,p,zatmo,sig
      USE GEOM, only : dyp,dxp
      USE DYNAMICS, only : phi,pu,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,1) :: PU0
      REAL*8 rho1,FLUX
      INTEGER I,J,IP1,IM1

C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL
      DPDY_BY_RHO=0.
      DPDY_BY_RHO_0=0.
      IM1=IM
      DO I=1,IM
        DO J=2,JM
          rho1=100.*pmid(1,i,j)/(rgas*t(i,j,1)*pk(1,i,j))
          FLUX=(PHI(I,J,1)-PHI(I,J-1,1))/DYP(J)
     2         +100.*(P(I,J)-P(I,J-1))*SIG(1)/(rho1*DYP(J))
          DPDY_BY_RHO(I,J)=DPDY_BY_RHO(I,J)+FLUX
          DPDY_BY_RHO(IM1,J)=DPDY_BY_RHO(IM1,J)+FLUX
          FLUX=(ZATMO(I,J)-ZATMO(I,J-1))/DYP(J)
     2         +100.*(P(I,J)-P(I,J-1))/(rho1*DYP(J))
          DPDY_BY_RHO_0(I,J)=DPDY_BY_RHO_0(I,J)+FLUX
          DPDY_BY_RHO_0(IM1,J)=DPDY_BY_RHO_0(IM1,J)+FLUX
        END DO
        IM1=I
      END DO
c
      DPDX_BY_RHO=0.
      DPDX_BY_RHO_0=0.
      I=IM
      DO IP1=1,IM
        PU(I,1,1)=0.
        PU(I,JM,1)=0.
        PU0(I,1,1)=0.
        PU0(I,JM,1)=0.
        DO J=2,JM-1
          rho1=100.*pmid(1,i,j)/(rgas*t(i,j,1)*pk(1,i,j))
          PU(I,J,1)=(PHI(IP1,J,1)-PHI(I,J,1))/DXP(J)
     2              +100.*(P(IP1,J)-P(I,J))*SIG(1)/(rho1*DXP(J))
          PU0(I,J,1)=(ZATMO(IP1,J)-ZATMO(I,J))/DXP(J)
     2               +100.*(P(IP1,J)-P(I,J))/(rho1*DXP(J))
        END DO
        I=IP1
      END DO
      CALL AVRX (PU(1,1,1))
      CALL AVRX (PU0(1,1,1))
      DO J=2,JM
        DO I=1,IM
          DPDX_BY_RHO(I,J)=DPDX_BY_RHO(I,J)+0.5*(PU(I,J,1)+PU(I,J-1,1))
          DPDX_BY_RHO_0(I,J)=DPDX_BY_RHO_0(I,J)+
     2                       0.5*(PU0(I,J,1)+PU0(I,J-1,1))
        END DO
      END DO
C****
      END SUBROUTINE PGRAD_PBL

      SUBROUTINE SDRAG(LMIN,DT1)
!@sum  SDRAG puts a drag on the winds on the top layers of atmosphere
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas
      USE MODEL_COM, only : im,jm,lm,psfmpt,u,v,sige,ptop,t,xcdlm
     *     ,bydsig,itime
      USE GEOM, only : cosv
      USE DAGCOM, only : aij, ij_wlm,ajl,ij_sdrag
      USE DYNAMICS, only : pk
      IMPLICIT NONE

!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var LMIN lowest level at which SDRAG is applied
C**** Normally 1 mb and up (inclusive)
C**** Note that this low level is only applied at the pole, elsewhere
C**** only top two layers are done (unless LMIN=LM i.e. only top layer)
      INTEGER, INTENT(IN) :: LMIN
      REAL*8 WL,RHO,CDN,X,BYPIJU
      INTEGER I,J,IP1,L

      BYPIJU=1./PSFMPT
      DO L=LMIN,LM 
      DO J=2,JM
      I=IM
      IF (COSV(J).LE..15.OR.L.GE.LM-1) THEN
      DO IP1=1,IM
        WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
        RHO=(PSFMPT*SIGE(L+1)+PTOP)/(RGAS*T(I,J,L)*PK(L,I,J))
        CDN=XCDLM(1)+XCDLM(2)*WL
        AIJ(I,J,IJ_WLM)=AIJ(I,J,IJ_WLM)+WL
        X=DT1*RHO*CDN*WL*GRAV*BYDSIG(L)*BYPIJU
        IF(X.GT.1) THEN
          write(99,*)'SDRAG: ITime,I,J,PSFMPT,X,RHO,CDN,U,V'
     *         ,ITime,I,J,PSFMPT,X,RHO,CDN,U(I,J,L),V(I,J,L)
     *         ,' If problem persists, winds are too high! '
     *         ,'Try setting XCDLM smaller.'
          X=1.
        END IF
        AJL(J,L,52) = AJL(J,L,52)-U(I,J,L)*X
        AIJ(I,J,IJ_SDRAG)=AIJ(I,J,IJ_SDRAG)-U(I,J,L)*X
        U(I,J,L)=U(I,J,L)*(1.-X)
        V(I,J,L)=V(I,J,L)*(1.-X)
        I=IP1
      END DO
      END IF
      END DO
      END DO
      RETURN
      END
