      subroutine init_MOM
!@sum  init_MOM sets an order dependent coefficient for AVRX
      USE DYNAMICS, only : XAVRX
      XAVRX = 1. ! for second order scheme, byrt2 for 4th order scheme
      CALL AVRX0
      RETURN
      END subroutine init_MOM

      SUBROUTINE ADVECV (PA,UT,VT,PB,U,V,P,DT1)
!@sum  ADVECV Advects momentum (incl. coriolis) using mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,modd5k
      USE DOMAIN_DECOMP, only : HALO_UPDATE, GRID,NORTH,SOUTH
      USE GEOM, only : fcor,dxyv,dxyn,dxys,dxv,ravpn,ravps
      USE DYNAMICS, only : pu,pv,pit,sd,spa,dut,dvt
      USE DIAG, only : diagcd
      IMPLICIT NONE
      REAL*8 U(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * V(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * P(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
      REAL*8 UT(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM), 
     * VT(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     * PA(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     * PB(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     * FD(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      REAL*8, SAVE :: SMASS(JM)

      INTEGER :: FROM,J_0,J_1,J_0H,J_1H,J_0S,J_1S
      INTEGER,SAVE :: IFIRST = 1
      INTEGER I,J,IP1,IM1,L,K  !@var I,J,IP1,IM1,L,K  loop variables
      REAL*8 VMASS,RVMASS,ALPH,PDT4,SDU,DT1,DT2,DT4,DT8,DT12,DT24
     *     ,FLUX,FLUXU,FLUXV

      REAL*8   VMASS2(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *         ASDU(IM,GRID%J_STRT_SKP:GRID%J_STOP_HALO,LM-1)
C****
      J_0 =GRID%J_STRT
      J_1 =GRID%J_STOP
      J_0H=GRID%J_STRT_HALO
      J_1H=GRID%J_STOP_HALO
      J_0S=GRID%J_STRT_SKP
      J_1S=GRID%J_STOP_SKP

         IF(MODD5K.LT.MRCH) CALL DIAG5F (U,V)
      IF(IFIRST.EQ.1) THEN
        IFIRST=0
        DO 10 J=J_0S,J_1
 10     SMASS(J)=PSFMPT*DXYV(J)
      END IF

      DT2=DT1/2.
      DT4=DT1/4.
      DT8=DT1/8.
      DT12=DT1/12.
      DT24=DT1/24.
C****
C**** SCALE UT AND VT WHICH MAY THEN BE PHYSICALLY INTERPRETED AS
C**** MOMENTUM COMPONENTS
C****
C     I=IM
C     DO 120 J=2,JM
C     DO 120 IP1=1,IM
C     VMASS=.5*((PA(I,J-1)+PA(IP1,J-1))*DXYN(J-1)
C    *  +(PA(I,J)+PA(IP1,J))*DXYS(J))
C     DO 110 L=1,LS1-1
C     UT(I,J,L)=UT(I,J,L)*VMASS*DSIG(L)
C 110 VT(I,J,L)=VT(I,J,L)*VMASS*DSIG(L)
C 120 I=IP1
C     DO 150 L=LS1,LM
C     DO 150 J=2,JM
C     VMASS=SMASS(J)*DSIG(L)
C     DO 150 I=1,IM
C     UT(I,J,L)=UT(I,J,L)*VMASS
C 150 VT(I,J,L)=VT(I,J,L)*VMASS
C     DUT=0.
C     DVT=0.
C
      FROM = SOUTH
      CALL HALO_UPDATE(GRID,PA,FROM)
CAOO Leave 1-D halo update for now   CALL HALO_UPDATE(GRID,DXYN,FROM)
      DO 110 J=J_0S,J_1
      I=IM
      DO 110 IP1=1,IM
         VMASS2(I,J)=.5*((PA(I,J-1)+PA(IP1,J-1))*DXYN(J-1)
     *                 +(PA(I,J)+PA(IP1,J))*DXYS(J))
         I=IP1
  110 CONTINUE
!$OMP  PARALLEL DO PRIVATE(J,L,VMASS)
      DO L=1,LM
        IF(L.LT.LS1) THEN  !  DO L=1,LS1-1
          DO J=J_0S,J_1
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
            UT(:,J,L)=UT(:,J,L)*VMASS2(:,J)*DSIG(L)
            VT(:,J,L)=VT(:,J,L)*VMASS2(:,J)*DSIG(L)
          END DO
        ELSE               !  DO L=LS1,LM
          DO J=J_0S,J_1
            VMASS=SMASS(J)*DSIG(L)
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
            UT(:,J,L)=UT(:,J,L)*VMASS
            VT(:,J,L)=VT(:,J,L)*VMASS
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C****
C**** BEGINNING OF LAYER LOOP
C****
!$OMP  PARALLEL DO PRIVATE(I,IP1,J,L,FLUX,FLUXU,FLUXV)
      DO 300 L=1,LM
C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****
      FROM = SOUTH
      CALL HALO_UPDATE(GRID,PU,FROM)
      I=IM
      DO 230 IP1=1,IM
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX
      DO 210 J=J_0S,J_1
      FLUX=DT12*(PU(IP1,J,L)+PU(IP1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
      FLUXU=FLUX*(U(IP1,J,L)+U(I,J,L))
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU
      DUT(I,J,L)  =DUT(I,J,L)  -FLUXU
      FLUXV=FLUX*(V(IP1,J,L)+V(I,J,L))
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV
  210 DVT(I,J,L)  =DVT(I,J,L)  -FLUXV
      FROM = NORTH
      CALL HALO_UPDATE(GRID,PV ,FROM)
      CALL HALO_UPDATE(GRID,U  ,FROM)
      CALL HALO_UPDATE(GRID,V  ,FROM)
CAOO no need to communicate, local compute      CALL HALO_UPDATE(GRID,DUT,FROM)
CAOO no need to communicate, local compute      CALL HALO_UPDATE(GRID,DVT,FROM)
      DO 220 J=J_0S,J_1S
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
!$OMP  END PARALLEL DO
C****
C**** VERTICAL ADVECTION OF MOMENTUM
C****
C     DO 310 L=1,LM-1
C     DO 310 J=2,JM
C     I=IM
C     DO 310 IP1=1,IM
C     SDU=DT2*((SD(I,J-1,L)+SD(IP1,J-1,L))*RAVPN(J-1)+
C    *  (SD(I,J,L)+SD(IP1,J,L))*RAVPS(J))
C     DUT(I,J,L)  =DUT(I,J,L)  +SDU*(U(I,J,L)+U(I,J,L+1))
C     DUT(I,J,L+1)=DUT(I,J,L+1)-SDU*(U(I,J,L)+U(I,J,L+1))
C     DVT(I,J,L)  =DVT(I,J,L)  +SDU*(V(I,J,L)+V(I,J,L+1))
C     DVT(I,J,L+1)=DVT(I,J,L+1)-SDU*(V(I,J,L)+V(I,J,L+1))
C 310 I=IP1
      FROM = SOUTH
      CALL HALO_UPDATE(GRID,SD,FROM)
CAOO Leave 1-D halo for now   CALL HALO_UPDATE(GRID,RAVPN,FROM)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM-1
      DO J=J_0S,J_1
         DO I=1,IM-1
           ASDU(I,J,L)=DT2*((SD(I,J-1,L)+SD(I+1,J-1,L))*RAVPN(J-1)+
     *                (SD(I,J,L)+SD(I+1,J,L))*RAVPS(J))
         END DO
           ASDU(IM,J,L)=DT2*((SD(IM,J-1,L)+SD(1,J-1,L))*RAVPN(J-1)+
     *                (SD(IM,J,L)+SD(1,J,L))*RAVPS(J))
      END DO
      END DO
!$OMP  END PARALLEL DO
      L=1
      DO J=J_0S,J_1
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
      END DO
!$OMP  PARALLEL DO PRIVATE(J,L)
      DO L=2,LM-1
        DO J=J_0S,J_1
         DUT(:,J,L)  =DUT(:,J,L)  -ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  -ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
        END DO
      END DO
!$OMP  END PARALLEL DO
      L=LM
      DO J=J_0S,J_1
         DUT(:,J,L)=DUT(:,J,L)-ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DVT(:,J,L)=DVT(:,J,L)-ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
      END DO
C**** CALL DIAGNOSTICS
         IF(MODD5K.LT.MRCH) CALL DIAG5D (4,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAGCD (1,U,V,DUT,DVT,DT1,PIT)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM
      DO J=J_0S,J_1
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
        VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
C****
C**** CORIOLIS FORCE
C****
!$OMP  PARALLEL DO PRIVATE(I,IM1,J,L,FD,PDT4,ALPH)
      DO L=1,LM
        IM1=IM
        DO I=1,IM
C         FD(I,1)=FCOR(1)*2.  -.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
C         FD(I,JM)=FCOR(JM)*2.+.5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
C****     Set the Coriolis term to zero at the Poles:
          IF(J_0 .EQ. 1)
     *    FD(I,J_0)= -.5*(SPA(IM1,J_0,L)+
     *                 SPA(I,J_0,L))*
     *                 DXV(J_0+1)
          IF(J_1 .EQ. JM)
     *    FD(I,J_1)=  .5*(SPA(IM1,J_1,L)+
     *                 SPA(I,J_1,L))*
     *                 DXV(J_1)
          IM1=I
        END DO
CAOO      FROM = NORTH
CAOO Leave 1-D hale update for now  CALL HALO_UPDATE(GRID,DXV ,FROM)
        DO J=J_0S,J_1S
          IM1=IM
          DO I=1,IM
            FD(I,J)=FCOR(J)+.25*(SPA(IM1,J,L)+SPA(I,J,L))
     *             *(DXV(J)-DXV(J+1))
            IM1=I
          END DO
        END DO
        FROM = SOUTH
        CALL HALO_UPDATE(GRID,P ,FROM)
        CALL HALO_UPDATE(GRID,FD,FROM)
        DO J=J_0S,J_1
          IM1=IM
          DO I=1,IM
            PDT4=DT8*(P(I,J-1,L)+P(I,J,L))
            IF(L.GE.LS1) PDT4=DT4*PSFMPT
            ALPH=PDT4*(FD(I,J)+FD(I,J-1))*DSIG(L)
            DUT(I,J,L)=DUT(I,J,L)+ALPH*V(I,J,L)
            DUT(IM1,J,L)=DUT(IM1,J,L)+ALPH*V(IM1,J,L)
            DVT(I,J,L)=DVT(I,J,L)-ALPH*U(I,J,L)
            DVT(IM1,J,L)=DVT(IM1,J,L)-ALPH*U(IM1,J,L)
            IM1=I
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO
C**** CALL DIAGNOSTICS
         IF(MODD5K.LT.MRCH) CALL DIAG5D (5,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAGCD (2,U,V,DUT,DVT,DT1)
C****
C**** ADD CORIOLIS FORCE INCREMENTS TO UT AND VT
C**** AND UNDO SCALING PERFORMED AT BEGINNING OF ADVECV
C****
      FROM = SOUTH
      CALL HALO_UPDATE(GRID,PB,FROM)
CAOO Leave 1-D halo update for now      CALL HALO_UPDATE(GRID,DXYN,FROM)
      DO J=J_0S,J_1
        I=IM
        DO IP1=1,IM
          VMASS2(I,J)=.5*((PB(I,J-1)+PB(IP1,J-1))*DXYN(J-1)
     *                 +(PB(I,J)+PB(IP1,J))*DXYS(J))
          I=IP1
        END DO
      END DO
!$OMP  PARALLEL DO PRIVATE(J,L,RVMASS)
      DO L=1,LM
        IF(L.LT.LS1) THEN  !  DO L=1,LS1-1
          DO J=J_0S,J_1
            VT(:,J,L)=(VT(:,J,L)+DVT(:,J,L))/(VMASS2(:,J)*DSIG(L))
            UT(:,J,L)=(UT(:,J,L)+DUT(:,J,L))/(VMASS2(:,J)*DSIG(L))
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
          END DO
        ELSE               !  DO L=LS1,LM
          DO J=J_0S,J_1
            RVMASS=1./(SMASS(J)*DSIG(L))
            UT(:,J,L)=(UT(:,J,L)+DUT(:,J,L))*RVMASS
            VT(:,J,L)=(VT(:,J,L)+DVT(:,J,L))*RVMASS
            DUT(:,J,L)=0.0
            DVT(:,J,L)=0.0
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE ADVECV
