      SUBROUTINE AFLUX (U,V,PA)
C****
C**** CALCULATE THE HORIZONTAL AIR MASS FLUXES
C**** AND VERTICAL AIR MASS FLUXES AS DETERMINED BY U, V AND P.
C**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
C****
      USE DYNAMICS, UGLOB=>U, VGLOB=>V, PGLOB=>P
      IMPLICIT NONE
      REAL*8 U(IM,JM,LM),V(IM,JM,LM),P(IM,JM) ! p is just workspace
      REAL*8 PHI,SPA
      COMMON/WORK3/PHI(IM,JM,LM),SPA(IM,JM,LM)
      REAL*8 FD,FLUXQ,DUMMYS,DUMMYN
      COMMON/WORK4/FD(IM,JM),FLUXQ(IM),DUMMYS(IM),DUMMYN(IM)
      REAL*8 PA(IM,JM)
      INTEGER I,J,L,IP1,IM1
      REAL*8 PUS,PUN,PVS,PVN,PBS,PBN,SDNP,SDSP
C****
C**** BEGINNING OF LAYER LOOP
C****
      L=LM
      P(:,:)=PSFMPT
C****
 2150 IF(L.EQ.LS1-1) P(:,:)=PA(:,:)
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
      PU(I,J,L)=.25*DYP(J)*SPA(I,J,L)*(P(I,J)+P(IP1,J))*DSIG(L)
 2165 I=IP1
 2166 CONTINUE
C**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX
      IM1=IM
      DO 2172 J=2,JM
      DO 2170 I=1,IM
      PV(I,J,L)=.25*DXV(J)*(V(I,J,L)+V(IM1,J,L))*
     *   (P(I,J)+P(I,J-1))*DSIG(L)
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
      PUS=.25*DYP(2)*PUS*P(1,1)*BYIM
      PUN=.25*DYP(JM-1)*PUN*P(1,JM)*BYIM
      PVS=PVS/FIM
      PVN=PVN/FIM
      DUMMYS(1)=0.
      DUMMYN(1)=0.
      DO 1120 I=2,IM
      DUMMYS(I)=DUMMYS(I-1)+(PV(I,2,L)-PVS)*BYDSIG(L)
 1120 DUMMYN(I)=DUMMYN(I-1)+(PV(I,JM,L)-PVN)*BYDSIG(L)
      PBS=0.
      PBN=0.
      DO 1130 I=1,IM
      PBS=PBS+DUMMYS(I)
 1130 PBN=PBN+DUMMYN(I)
      PBS=PBS*BYIM
      PBN=PBN*BYIM
      DO 1140 I=1,IM
      SPA(I,1,L)=4.*(PBS-DUMMYS(I)+PUS)/(DYP(2)*P(1,1))
      SPA(I,JM,L)=4.*(DUMMYN(I)-PBN+PUN)/(DYP(JM-1)*P(1,JM))
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
      L=L-1
      IF (L.GE.1) GO TO 2150
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
C****
C**** CALCULATE UPDATED COLUMN PRESSURES AS
C**** DETERMINED BY DT1 AND THE CURRENT AIR MASS FLUXES.
C****
      USE DYNAMICS, PGLOB=>P
      IMPLICIT NONE
      REAL*8 P(IM,JM)
      REAL*8 FD
      COMMON/WORK4/FD(IM,JM)
      REAL*8 PA(IM,JM)
      INTEGER I,J,L,K,IMAX  !@var I,J,L,K  loop variables
      REAL*8 DT1

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

      SUBROUTINE ADVECV (PA,UT,VT,PB,U,V,P,DT1)
C****
C**** ADVECT MOMENTUM (INCLUDING THE CORIOLIS FORCE)
C**** AS DETERMINED BY DT1 AND THE CURRENT AIR MASS FLUXES
C****
      USE DYNAMICS, UGLOB=>U, VGLOB=>V, PGLOB=>P
      IMPLICIT NONE
      REAL*8 U(IM,JM,LM),V(IM,JM,LM),P(IM,JM)
      REAL*8 PHI,SPA
      COMMON/WORK3/PHI(IM,JM,LM),SPA(IM,JM,LM)
      REAL*8 FD,DUT,DVT
      COMMON/WORK4/FD(IM,JM)
      COMMON/WORK5/DUT(IM,JM,LM),DVT(IM,JM,LM)
      REAL*8 UT(IM,JM,LM),VT(IM,JM,LM),PA(IM,JM),PB(IM,JM)
      REAL*8, SAVE :: SMASS(JM)

      INTEGER,SAVE :: IFIRST = 1
      INTEGER I,J,IP1,IM1,L,K  !@var I,J,IP1,IM1,L,K  loop variables
      REAL*8 VMASS,RVMASS,ALPH,PDT4,SDU,DT1,DT2,DT4,DT8,DT12,DT24
     *     ,FLUX,FLUXU,FLUXV
C****
         IF(MODD5K.LT.MRCH) CALL DIAG5F (U,V)
      IF(IFIRST.NE.1) GO TO 50
      IFIRST=0
c      JMM2=JM-2
c      IJL2=IM*JM*LM*2
      DO 10 J=2,JM
   10 SMASS(J)=PSFMPT*DXYV(J)

   50 CONTINUE
      DT2=DT1/2.
      DT4=DT1/4.
      DT8=DT1/8.
      DT12=DT1/12.
      DT24=DT1/24.
C****
C**** SCALE UT AND VT WHICH MAY THEN BE PHYSICALLY INTERPRETED AS
C**** MOMENTUM COMPONENTS
C****
      I=IM
      DO 120 J=2,JM
      DO 120 IP1=1,IM
      VMASS=.5*((PA(I,J-1)+PA(IP1,J-1))*DXYN(J-1)
     *  +(PA(I,J)+PA(IP1,J))*DXYS(J))
      DO 110 L=1,LS1-1
      UT(I,J,L)=UT(I,J,L)*VMASS*DSIG(L)
  110 VT(I,J,L)=VT(I,J,L)*VMASS*DSIG(L)
  120 I=IP1
      DO 150 L=LS1,LM
      DO 150 J=2,JM
      VMASS=SMASS(J)*DSIG(L)
      DO 150 I=1,IM
      UT(I,J,L)=UT(I,J,L)*VMASS
  150 VT(I,J,L)=VT(I,J,L)*VMASS
      DUT=0.
      DVT=0.
C****
C**** BEGINNING OF LAYER LOOP
C****
      DO 300 L=1,LM
C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****
      I=IM
      DO 230 IP1=1,IM
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX
      DO 210 J=2,JM
      FLUX=DT12*(PU(IP1,J,L)+PU(IP1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
      FLUXU=FLUX*(U(IP1,J,L)+U(I,J,L))
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU
      DUT(I,J,L)  =DUT(I,J,L)  -FLUXU
      FLUXV=FLUX*(V(IP1,J,L)+V(I,J,L))
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV
  210 DVT(I,J,L)  =DVT(I,J,L)  -FLUXV
      DO 220 J=2,JM-1
C**** CONTRIBUTION FROM THE SOUTH-NORTH MASS FLUX
      FLUX=DT12*(PV(I,J,L)+PV(IP1,J,L)+PV(I,J+1,L)+PV(IP1,J+1,L))
      FLUXU=FLUX*(U(I,J,L)+U(I,J+1,L))
      DUT(I,J+1,L)=DUT(I,J+1,L)+FLUXU
      DUT(I,J,L)  =DUT(I,J,L)  -FLUXU
      FLUXV=FLUX*(V(I,J,L)+V(I,J+1,L))
      DVT(I,J+1,L)=DVT(I,J+1,L)+FLUXV
      DVT(I,J,L)=  DVT(I,J,L)  -FLUXV
C**** CONTRIBUTION FROM THE SOUTHWEST-NORTHEAST MASS FLUX
      FLUX=DT24*(PU(IP1,J,L)+PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
      FLUXU=FLUX*(U(IP1,J+1,L)+U(I,J,L))
      DUT(IP1,J+1,L)=DUT(IP1,J+1,L)+FLUXU
      DUT(I,J,L)=    DUT(I,J,L)    -FLUXU
      FLUXV=FLUX*(V(IP1,J+1,L)+V(I,J,L))
      DVT(IP1,J+1,L)=DVT(IP1,J+1,L)+FLUXV
      DVT(I,J,L)=    DVT(I,J,L)    -FLUXV
C**** CONTRIBUTION FROM THE SOUTHEAST-NORTHWEST MASS FLUX
      FLUX=DT24*(-PU(IP1,J,L)-PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
      FLUXU=FLUX*(U(I,J+1,L)+U(IP1,J,L))
      DUT(I,J+1,L)=DUT(I,J+1,L)+FLUXU
      DUT(IP1,J,L)=DUT(IP1,J,L)-FLUXU
      FLUXV=FLUX*(V(I,J+1,L)+V(IP1,J,L))
      DVT(I,J+1,L)=DVT(I,J+1,L)+FLUXV
  220 DVT(IP1,J,L)=DVT(IP1,J,L)-FLUXV
  230 I=IP1
  300 CONTINUE
C****
C**** VERTICAL ADVECTION OF MOMENTUM
C****
      DO 310 L=1,LM-1
      DO 310 J=2,JM
      I=IM
      DO 310 IP1=1,IM
      SDU=DT2*((SD(I,J-1,L)+SD(IP1,J-1,L))*RAVPN(J-1)+
     *  (SD(I,J,L)+SD(IP1,J,L))*RAVPS(J))
      DUT(I,J,L)  =DUT(I,J,L)  +SDU*(U(I,J,L)+U(I,J,L+1))
      DUT(I,J,L+1)=DUT(I,J,L+1)-SDU*(U(I,J,L)+U(I,J,L+1))
      DVT(I,J,L)  =DVT(I,J,L)  +SDU*(V(I,J,L)+V(I,J,L+1))
      DVT(I,J,L+1)=DVT(I,J,L+1)-SDU*(V(I,J,L)+V(I,J,L+1))
  310 I=IP1
C**** CALL DIAGNOSTICS
         IF(MODD5K.LT.MRCH) CALL DIAG5D (4,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAG9D (1,DT1,U,V,DUT,DVT,PIT)
      DO L=1,LM
         DO J=2,JM
            DO I=1,IM
               UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
               VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
               DUT(I,J,L)=0.
               DVT(I,J,L)=0.
            ENDDO
         ENDDO
      ENDDO
C****
C**** CORIOLIS FORCE
C****
      DO 430 L=1,LM
      IM1=IM
      DO 405 I=1,IM
C     FD(I,1)=FCOR(1)*2.-.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
C     FD(I,JM)=FCOR(JM)*2.+.5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
C**** Set the Coriolis term to zero at the Poles:
      FD(I,1)=  -.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
      FD(I,JM)=  .5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
  405 IM1=I
      DO 410 J=2,JM-1
      DO 410 I=1,IM
      FD(I,J)=FCOR(J)+.25*(SPA(IM1,J,L)+SPA(I,J,L))*(DXV(J)-DXV(J+1))
  410 IM1=I
  415 CONTINUE
      DO 425 J=2,JM
      DO 420 I=1,IM
      PDT4=DT8*(P(I,J-1)+P(I,J))
      IF(L.GE.LS1) PDT4=DT4*PSFMPT
      ALPH=PDT4*(FD(I,J)+FD(I,J-1))*DSIG(L)
      DUT(I,J,L)=DUT(I,J,L)+ALPH*V(I,J,L)
      DUT(IM1,J,L)=DUT(IM1,J,L)+ALPH*V(IM1,J,L)
      DVT(I,J,L)=DVT(I,J,L)-ALPH*U(I,J,L)
      DVT(IM1,J,L)=DVT(IM1,J,L)-ALPH*U(IM1,J,L)
  420 IM1=I
  425 CONTINUE
  430 CONTINUE
C**** CALL DIAGNOSTICS, ADD CORIOLIS FORCE INCREMENTS TO UT AND VT
         IF(MODD5K.LT.MRCH) CALL DIAG5D (5,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAG9D (2,DT1,U,V,DUT,DVT,PIT)
      DO L=1,LM
         DO J=2,JM
            DO I=1,IM
               UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
               VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
               DUT(I,J,L)=0.
               DVT(I,J,L)=0.
            ENDDO
         ENDDO
      ENDDO
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF ADVECV
C****
      I=IM
      DO 520 J=2,JM
      DO 520 IP1=1,IM
      VMASS=.5*((PB(I,J-1)+PB(IP1,J-1))*DXYN(J-1)
     *  +(PB(I,J)+PB(IP1,J))*DXYS(J))
      DO 510 L=1,LS1-1
      VT(I,J,L)=VT(I,J,L)/(VMASS*DSIG(L))
  510 UT(I,J,L)=UT(I,J,L)/(VMASS*DSIG(L))
  520 I=IP1
      DO 550 L=LS1,LM
      DO 550 J=2,JM
      RVMASS=1./(SMASS(J)*DSIG(L))
      DO 550 I=1,IM
      UT(I,J,L)=UT(I,J,L)*RVMASS
  550 VT(I,J,L)=VT(I,J,L)*RVMASS
      RETURN
      END SUBROUTINE ADVECV

      SUBROUTINE PGF (UT,VT,PB,U,V,T,SZ,P,DT1)
C****
C**** ADD TO MOMENTUM THE TENDENCIES DETERMINED BY
C**** THE PRESSURE GRADIENT FORCE
C****
      USE DYNAMICS, UGLOB=>U, VGLOB=>V, PGLOB=>P, TGLOB=>T
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,LM) :: U,V,T
      REAL*8, DIMENSION(IM,JM) :: P,RFDUX

      REAL*8, DIMENSION(IM,JM,LM) :: PHI,SPA
      COMMON/WORK3/PHI,SPA

      REAL*8 FD(IM,JM),FLUXQ(IM),DUMMYS(IM),DUMMYN(IM)
      COMMON/WORK4/FD,FLUXQ,DUMMYS,DUMMYN

      REAL*8 DUT(IM,JM,LM),DVT(IM,JM,LM)
      COMMON/WORK5/DUT,DVT

      REAL*8 UT(IM,JM,LM),VT(IM,JM,LM),TT(IM,JM,LM),
     *  PA(IM,JM),PB(IM,JM),QT(IM,JM,LM)

      REAL*8 PKE(LM+1)
      REAL*8 SZ(IM,JM,LM),DT4,DT1
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
!     REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDUX,RFDU,PHIDN,FACTORS
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,      RFDU,PHIDN,FACTORS
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
      IMAX=IMAXJ(J)
      DO 330 I=1,IMAX
      PIJ=P(I,J)
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
      DO 3236 L=1,LS1-1
      DO 3236 J=2,JM
      IM1=IM
      DO 3234 I=1,IM
      FLUX=DT4*((P(I,J)+P(I,J-1))*(PHI(I,J,L)-PHI(I,J-1,L))+
     *  (SPA(I,J,L)+SPA(I,J-1,L))*(P(I,J)-P(I,J-1)))*DXV(J)*DSIG(L)
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3234 IM1=I
 3236 CONTINUE
      DO 3246 L=LS1,LM
      DO 3246 J=2,JM
!     FACTORS = 2.*DT4*PSFMPT*DXV(J)*DSIG(L)
!     FACTORS =               DXV(J)*DSIG(L)
      IM1=IM
      DO 3244 I=1,IM
      FLUX=2.*DT4*PSFMPT*(PHI(I,J,L)-PHI(I,J-1,L))*DXV(J)*DSIG(L)
!     FLUX=2.*DT4*PSFMPT*(PHI(I,J,L)-PHI(I,J-1,L))*FACTORS
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3244 IM1=I
 3246 CONTINUE
C**** SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
      DO 3295 L=1,LS1-1
      DO 3275 I=1,IM
      PU(I,1,L)=0.
 3275 PU(I,JM,L)=0.
      I=IM
      DO 3290 J=2,JM-1
      DO 3280 IP1=1,IM
      PU(I,J,L)=(P(IP1,J)+P(I,J))*(PHI(IP1,J,L)-PHI(I,J,L))+
     *  (SPA(IP1,J,L)+SPA(I,J,L))*(P(IP1,J)-P(I,J))
 3280 I=IP1
 3290 CONTINUE
      CALL AVRX (PU(1,1,L))
      DO 3294 J=2,JM
!     FACTORS = -DT4*DYV(J)*DSIG(L)
      DO 3294 I=1,IM
 3294 DUT(I,J,L)=DUT(I,J,L)-DT4*DYV(J)*(PU(I,J,L)+PU(I,J-1,L))*DSIG(L)
!3294 DUT(I,J,L)=DUT(I,J,L)-FACTORS   *(PU(I,J,L)+PU(I,J-1,L))
 3295 CONTINUE
      DO 3340 L=LS1,LM
      DO 3315 I=1,IM
      PU(I,1,L)=0.
 3315 PU(I,JM,L)=0.
      I=IM
      DO 3325 J=2,JM-1
      DO 3320 IP1=1,IM
      PU(I,J,L)=2.*PSFMPT*(PHI(IP1,J,L)-PHI(I,J,L))
 3320 I=IP1
 3325 CONTINUE
      CALL AVRX (PU(1,1,L))
      DO 3330 J=2,JM
      DO 3330 I=1,IM
 3330 DUT(I,J,L)=DUT(I,J,L)-DT4*DYV(J)*(PU(I,J,L)+PU(I,J-1,L))*DSIG(L)
 3340 CONTINUE
C**** CALL DIAGNOSTICS
      IF(MRCH.LE.0) GO TO 500
         IF(MODD5K.LT.MRCH) CALL DIAG5D (6,MRCH,DUT,DVT)
         IF(MODD5K.LT.MRCH) CALL DIAG9D (3,DT1,U,V,DUT,DVT,PIT)
C****
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
C****
  500 CONTINUE
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
      I=IM
      DO 3540 I=1,IM
      IF(L.LT.LS1) RFDU=RFDUX(I,J)*BYDSIG(L)
      VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)*RFDU
      UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)*RFDU
 3540 CONTINUE
 3550 CONTINUE
      RETURN
      END SUBROUTINE PGF

      SUBROUTINE AVRX (X)
C****
C**** THIS SUBROUTINE SMOOTHES THE ZONAL MASS FLUX AND GEOPOTENTIAL
C**** GRADIENTS NEAR THE POLES TO HELP AVOID COMPUTATIONAL INSTABILITY.
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
C****
      USE DYNAMICS
      IMPLICIT NONE
      REAL*8 X(IM,JM)
      REAL*8, SAVE :: SM(IMH,JM),DRAT(JM)
      INTEGER, SAVE :: NMIN(JM)
      REAL*8, SAVE, DIMENSION(IMH) :: BYSN
      REAL*8, DIMENSION(0:IMH) :: AN,BN
c      COMMON/AVRXX/BYSN(IMH),AN(0:IMH),BN(0:IMH)

      INTEGER,SAVE :: IFIRST = 1
      INTEGER J,N

      IF (IFIRST.NE.1) GO TO 100
      IFIRST=0
C     CALL FFT0(IM)
      DO 30 N=1,IMH
   30 BYSN(N)=1./SIN(.5*DLON*N)
      DO 50 J=2,JM-1
      DRAT(J) = DXP(J)/DYP(3)
      DO 40 N=IMH,1,-1
      SM(N,J) = BYSN(N)*DRAT(J)
      IF(SM(N,J).GT.1.) THEN
         NMIN(J) = N+1
         GO TO 50
      ENDIF
   40 CONTINUE
   50 CONTINUE
C****
  100 DO 140 J=2,JM-1
      IF (DRAT(J).GT.1) GO TO 140
      CALL FFT (X(1,J),AN,BN)
      DO 130 N=NMIN(J),IMH-1
      AN(N)=SM(N,J)*AN(N)
  130 BN(N)=SM(N,J)*BN(N)
      AN(IMH) = SM(IMH,J)*AN(IMH)
      CALL FFTI(AN,BN,X(1,J))
  140 CONTINUE
      RETURN
      END SUBROUTINE AVRX

      SUBROUTINE FILTER
C****
C**** THIS SUBROUTINE PERFORMS AN 8-TH ORDER SHAPIRO FILTER ON
C**** SELECTED PROGNOSTIC QUANTITIES IN THE ZONAL DIRECTION
C****
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER
C****        3  SMOOTH P AND T
C****
      USE DYNAMICS
      IMPLICIT NONE
      REAL*8 X,XS,Y
      COMMON/WORK2/X(IM,JM),XS(IM),Y(IM,JM)
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
      CALL SHAP1D (8)
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
      CALL SHAP1D (8)
      DO 240 J=2,JM-1
      DO 240 I=1,IM
  240 T(I,J,L)=X(I,J)/Y(I,J)
  260 CONTINUE
      DO 280 L=LS1,LM
      DO 270 J=2,JM-1
      DO 270 I=1,IM
  270 X(I,J)=T(I,J,L)
      CALL SHAP1D (8)
      DO 280 J=2,JM-1
      DO 280 I=1,IM
  280 T(I,J,L)=X(I,J)
      RETURN
      END SUBROUTINE FILTER

      SUBROUTINE FLTRUV(U,V)
C**********************************************************************
C**** THIS ROUTINE FILTERS THE TWO GRIDLENGTH NOISE FROM THE
C**** VELOCITY FIELDS (U,V) IN BOTH DIMENSIONS WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES.
C**********************************************************************
      USE DYNAMICS, UGLOB=>U, VGLOB=>V
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM,LM) :: U,V
      REAL*8 X(IM),Y(0:JM+1),F2D(IM,JM)
      LOGICAL*4, SAVE :: QFILY,QFILX,QFIL2D
      INTEGER, SAVE :: NSHAP,ISIGN
      INTEGER I,J,L,N  !@var I,J,L,N  loop variables
      REAL*8, SAVE :: X4TON
      REAL*8 Y4TO8,YJ,YJM1,X1,XI,XIM1
      INTEGER,SAVE :: IFIRST = 1

      IF(IFIRST.EQ.1) THEN
         IFIRST = 0
         CALL FFT0(IM)
         QFILY  = .FALSE.
         QFILX  = .TRUE.
         QFIL2D = .FALSE.
         NSHAP  = 8
         ISIGN  = (-1.0)**(NSHAP-1)
         X4TON  = 1./(4.**NSHAP)
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
  240 U(I,J,L) = U(I,J,L) + ISIGN*X(I)*X4TON
      GOTO 270
  250 DO 255,J=2,JM-1
      F2D(1,J)=.125*(U(IM,J,L)+U(1,J+1,L)+U(2,J,L)+
     *                  U(1,J-1,L)-4.*U(1,J,L))+
     *      1./16.*(U(IM,J-1,L)+U(IM,J+1,L)+U(2,J+1,L)+
     *              U(2,J-1,L)-4.*U(1,J,L))
      DO 260,I=2,IM-1
      F2D(I,J)=.125*(U(I-1,J,L)+U(I,J+1,L)+U(I+1,J,L)+
     *                  U(I,J-1,L)-4.*U(I,J,L))+
     *      1./16.*(U(I-1,J-1,L)+U(I-1,J+1,L)+U(I+1,J+1,L)+
     *              U(I+1,J-1,L)-4.*U(I,J,L))
  260 CONTINUE
      F2D(IM,J)=.125*(U(IM-1,J,L)+U(IM,J+1,L)+U(1,J,L)+
     *                  U(IM,J-1,L)-4.*U(IM,J,L))+
     *      1./16.*(U(IM-1,J-1,L)+U(IM-1,J+1,L)+U(1,J+1,L)+
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
  340 V(I,J,L) = V(I,J,L) + ISIGN*X(I)*X4TON
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

      SUBROUTINE SHAP1D (NORDER)
C****
C**** THIS SUBROUTINE SMOOTHES THE ARRAY X IN THE ZONAL DIRECTION
C**** USING AN N-TH ORDER SHAPIRO FILTER.  N MUST BE EVEN.
C****
      USE E001M12_COM, only :im,jm
      IMPLICIT NONE
      REAL*8 X,XS,RE4TON,XS1,XSIM1,XSI
      INTEGER I,J,N,NORDER    !@var I,J,N  loop variables

      COMMON/WORK2/X(IM,JM),XS(IM)
      RE4TON=1./4.**NORDER
      DO 180 J=2,JM-1
      DO 120 I=1,IM
  120 XS(I)=X(I,J)
      DO 160 N=1,NORDER
      XS1=XS(1)
      XSIM1=XS(IM)
      DO 140 I=1,IM-1
      XSI=XS(I)
      XS(I)=XSIM1-XSI-XSI+XS(I+1)
  140 XSIM1=XSI
  160 XS(IM)=XSIM1-XS(IM)-XS(IM)+XS1
      DO 180 I=1,IM
  180 X(I,J)=X(I,J)-XS(I)*RE4TON
      RETURN
      END SUBROUTINE SHAP1D

      SUBROUTINE CALC_AMPK(LMAX)
!@sum  CALC_AMPK calculate air mass and pressure functions
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE DYNAMICS
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
C**** Compute AMP: kg air * grav/100, including constant pressure strat
      USE DYNAMICS, PGLOB=>P
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
      end subroutine calc_amp


      SUBROUTINE DYNAM
C****
C**** INTEGRATE DYNAMIC TERMS
C****
      USE CONSTANT, only : by3
      USE E001M12_COM, only : im,jm,lm,u,v,t,p,q,wm,dsig,NIdyn,dt,MODD5K
     *     ,NSTEP,NDA5K,ndaa,mrch,psfmpt,ls1
      USE GEOM, only : dyv,dxv
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE DYNAMICS, only : ptold,pu,pv,pit,sd

      USE DAGCOM, only : aij,ij_fpeu,ij_fpev,ij_fqu,ij_fqv,ij_fmv,ij_fmu
     *     ,ij_fgzu,ij_fgzv

      IMPLICIT NONE
      REAL*8 PRAT
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: UT,VT,TT,TZ,TZT,WMT,MA
      COMMON/WORK6/UT,VT,TT,TZ,TZT,WMT,MA,PRAT(IM,JM)
      REAL*8 UX,VX
      COMMON/WORK2/UX(IM,JM,LM),VX(IM,JM,LM)
      REAL*8 PA(IM,JM),PB(IM,JM),PC(IM,JM),FPEU(IM,JM),FPEV(IM,JM),
     *          FWVU(IM,JM),FWVV(IM,JM)

      REAL*8 PHI
      COMMON/WORK3/PHI(IM,JM,LM)
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: DUT,DVT
      COMMON/WORK5/DUT,DVT

      REAL*8 DTFS,DTLF,PP,UU,VV
      INTEGER I,J,L,IP1,IM1   !@var I,J,L,IP1,IM1  loop variables
      INTEGER NS    !? ,NIdynO

!?    NIdynO=MOD(NIdyn,2)   ! NIdyn odd is currently not an option
      DTFS=DT*2./3.
      DTLF=2.*DT
      PTOLD(:,:) = P(:,:)
      UX(:,:,:)  = U(:,:,:)
      UT(:,:,:)  = U(:,:,:)
      VX(:,:,:)  = V(:,:,:)
      VT(:,:,:)  = V(:,:,:)
      WMT(:,:,:) = WM(:,:,:)
C*** copy z-moment of temperature into contiguous memory
      tz(:,:,:) = tmom(mz,:,:,:)
      PA(:,:) = P(:,:)
      PB(:,:) = P(:,:)
      PC(:,:) = P(:,:)
C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      NS=0
      MRCH=0
C     CALL DYNAM (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS)
      CALL AFLUX (U,V,P)
      CALL ADVECM (P,PB,DTFS)
      CALL ADVECV (P,UX,VX,PB,U,V,P,DTFS)
      CALL PGF (UX,VX,PB,U,V,T,TZ,P,DTFS)
      CALL FLTRUV(UX,VX)
!?    IF (NIdynO.EQ.1) GO TO 320
C**** INITIAL BACKWARD STEP IS ODD, QT = Q + DT*F(QX)
      MRCH=-1
C     CALL DYNAM (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT)
      CALL AFLUX (UX,VX,PB)
      CALL ADVECM (P,PA,DT)
      CALL ADVECV (P,UT,VT,PA,UX,VX,PB,DT)
      CALL PGF (UT,VT,PA,UX,VX,T,TZ,PB,DT)
      CALL FLTRUV(UT,VT)
      GO TO 360
C**** INITIAL BACKWARD STEP IS EVEN, Q = Q + DT*F(QX)
!?320 NS=1
!?       MODD5K=MOD(NSTEP+NS-NIdyn+NDA5K*NIdyn+2,NDA5K*NIdyn+2)
!?    MRCH=1
C     CALL DYNAM (U,V,T,P,Q,UX,VX,TX,PX,QT,DT)
CD       DIAGA SHOULD BE CALLED HERE BUT THEN ARRAYS MUST BE CHANGED
C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)
  340 MRCH=-2
C     CALL DYNAM (UT,VT,TT,PT,QT,U,V,T,P,Q,DTLF)
      CALL AFLUX (U,V,P)
      CALL ADVECM (PA,PB,DTLF)
      CALL ADVECV (PA,UT,VT,PB,U,V,P,DTLF)
      CALL PGF (UT,VT,PB,U,V,T,TZ,P,DTLF)
      CALL FLTRUV(UT,VT)
C**** LOAD PB TO PA
      PA(:,:) = PB(:,:)
C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)
  360 NS=NS+2
         MODD5K=MOD(NSTEP+NS-NIdyn+NDA5K*NIdyn+2,NDA5K*NIdyn+2)
      MRCH=2
C     CALL DYNAM (U,V,T,P,Q,UT,VT,TT,PT,QT,DTLF)
      CALL AFLUX (UT,VT,PA)
      CALL ADVECM (PC,P,DTLF)
      CALL ADVECV (PC,U,V,P,UT,VT,PA,DTLF)
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
      CALL PGF (U,V,P,UT,VT,TT,TZT,PC,DTLF)
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
           CALL DIAGA (UT,VT,TT,PB,Q,PIT,SD)
           CALL DIAGB (UT,VT,TT,PB,Q,WMT,DUT,DVT)
         ENDIF
      IF (NS.LT.NIdyn) GO TO 340
C**** Scale WM mixing ratios to conserve liquid water
      PRAT(:,:)=PTOLD(:,:)/P(:,:)
      DO L=1,LS1-1
        WM(:,:,L)=WM(:,:,L)*PRAT(:,:)
      END DO
      RETURN
      END SUBROUTINE DYNAM

