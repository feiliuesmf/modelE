      subroutine init_MOM
!@sum  init_MOM sets an order dependent coefficient for AVRX
      USE DYNAMICS, only : XAVRX
      USE constant, only : byrt2
      XAVRX = byrt2 ! for 4th order scheme, 1 for 2nd order scheme
      RETURN
      END subroutine init_MOM

      SUBROUTINE ADVECV (PA,UT,VT,PB,U,V,P,DT1)
!@sum  ADVECV Advects momentum (incl. coriolis) using 4-th order scheme
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,modd5k
      USE GEOM, only : fcor,dxyp,dxv,ravpn,ravps
      USE DYNAMICS, only : pu,pv,pit,sd,spa,dut,dvt
      USE constant, only : by3,by6,by12
      IMPLICIT NONE

      REAL*8 U(IM,JM,LM),V(IM,JM,LM),P(IM,JM,LM)
      REAL*8 UT(IM,JM,LM),VT(IM,JM,LM),PA(IM,JM),PB(IM,JM),FD(IM,0:JM+1)
      REAL*8 FLUXQ(IM),DUMMYS(IM),DUMMYN(IM)
      REAL*8 FDU1(IM,JM+1),FDU2(IM,JM)
      REAL*8 FX(IM,0:JM+1),GY(IM,JM),FX1(IM,JM),GY1(IM,JM)
      REAL*8,SAVE :: FDU(IM,JM+1,LM)

      INTEGER,SAVE :: IFIRST = 1
      INTEGER I,J,IP1,IM1,L,K  !@var I,J,IP1,IM1,L,K  loop variables
      REAL*8 VMASS,RVMASS,ALPH,PDT4,SDU,DT1,DT2,DT3,DT6,DT8,DT24
     *     ,FLUX,FLUXU,FLUXV
C****
         IF(MODD5K.LT.MRCH) CALL DIAG5F (U,V)
C****
C**** Find scaling factors FDU for const. pressure layers (time indep)
C****
      IF (IFIRST .eq. 1) then
        DO J=2,JM-1
        DO I=1,IM
          FD(I,J)=PSFMPT*DXYP(J)
        END DO
        END DO
        FD(:,0)    = 0.
        FD(:,JM+1) = 0.
        FD(:,1)    = 2.*PSFMPT*DXYP(1)
        FD(:,JM)   = 2.*PSFMPT*DXYP(JM)

        I=IM
        DO J=1,JM+1
        DO IP1=1,IM
          FDU1(I,J) = .25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
          I=IP1
        END DO
        END DO

        I=IM
        IM1=I-1
        DO J=2,JM
        DO IP1=1,IM
         FDU2(I,J)=.25*(FDU1(IP1,J)+FDU1(I,J+1)+FDU1(I,J-1)+FDU1(IM1,J))
         IM1=I
         I=IP1
        END DO
        END DO

        DO J=2,JM
        DO I=1,IM
          FDU2(I,J)=(5.*FDU1(I,J)-2.*FDU2(I,J))*by3
        END DO
        END DO

        DO L=LS1,LM
        DO J=2,JM
        DO I=1,IM
          FDU(I,J,L)=FDU2(I,J)*DSIG(L)
        END DO
        END DO
        END DO
      END IF

      DT2=DT1/2.
      DT3=DT1*by3
      DT6=DT1*by6
      DT8=DT1/8.
      DT24=DT1*by12/2.
C****
C**** Scale UT and VT which may then be physically interpreted as
C**** the momentum components
C****
C**** USE FOURTH ORDER ON THE INTERIOR.
C**** CONSIDER FLUXES PU AND PV TO BE ZERO BEYOND POLES.
C**** CONSIDER P TO BE ZERO BEYOND THE POLES
      DO J=2,JM-1
      DO I=1,IM
        FD(I,J)   = PA(I,J)*DXYP(J)
      END DO
      END DO
      FD(:,0)     = 0.
      FD(:,JM+1)  = 0.
      FD(:,1)     = 2.*PA(1,1)*DXYP(1)
      FD(:,JM)    = 2.*PA(1,JM)*DXYP(JM)

      I=IM
      DO J=1,JM+1
      DO IP1=1,IM
        FDU1(I,J) = .25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
        I=IP1
      END DO
      END DO

      I=IM
      IM1=I-1
      DO J=2,JM
      DO IP1=1,IM
        FDU2(I,J)=.25*(FDU1(IP1,J)+FDU1(I,J+1)+FDU1(I,J-1)+FDU1(IM1,J))
        IM1=I
        I=IP1
      END DO
      END DO

      DO J=2,JM
      DO I=1,IM
        FDU2(I,J)=(5.*FDU1(I,J)-2.*FDU2(I,J))*by3
      END DO
      END DO

      DO L=1,LS1-1
      DO J=2,JM
      DO I=1,IM
        FDU(I,J,L)=FDU2(I,J)*DSIG(L)
      END DO
      END DO
      END DO

      DUT=0.
      DVT=0.
      DO L=1,LM
      DO J=2,JM
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)*FDU(I,J,L)
        VT(I,J,L)=VT(I,J,L)*FDU(I,J,L)
      END DO
      END DO
      END DO
C****
C**** BEGINNING OF LAYER LOOP : FIND HORIZONTAL FLUXES
C****
      DO 300 L=1,LM
      I=IM
      DO 175 J=1,JM
      DO 170 IP1=1,IM
      FX(I,J)=.5*(PU(IP1,J,L)+PU(I,J,L))
  170 I=IP1
  175 CONTINUE
C****
      DO 180 J=2,JM-1
      DO 180 I=1,IM
      GY(I,J)=.5*(PV(I,J,L)+PV(I,J+1,L))
  180 CONTINUE
C**** ASSUME FLUXES PU AND PV ARE ZERO BEYOND POLE
      DO 185 I=1,IM
      GY(I,1)=.5*PV(I,2,L)
      GY(I,JM)=.5*PV(I,JM,L)
      FX(I,0)=0.
      FX(I,JM+1)=0.
  185 CONTINUE
C****
      I=IM
      IM1=I-1
      DO 192 J=1,JM
      DO 190 IP1=1,IM
      FX1(I,J)=.25*(FX(IP1,J)+FX(IM1,J)+FX(I,J+1)+FX(I,J-1))
      IM1=I
  190 I=IP1
  192 CONTINUE
      DO 195 J=2,JM-1
      DO 193 IP1=1,IM
      GY1(I,J)=.25*(GY(IP1,J)+GY(IM1,J)+GY(I,J+1)+GY(I,J-1))
      IM1=I
  193 I=IP1
  195 CONTINUE
C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****
      I=IM
C**** Contribution from the West-East mass flux
      DO 215 J=2,JM
      DO 210 IP1=1,IM
      FLUX = DT3*(FX(I,J)+FX(I,J-1))
      FLUXU = FLUX*(U(IP1,J,L)+U(I,J,L))
      DUT(IP1,J,L) = DUT(IP1,J,L) + FLUXU
      DUT(I,J,L)   = DUT(I,J,L)   - FLUXU
      FLUXV = FLUX*(V(IP1,J,L)+V(I,J,L))
      DVT(IP1,J,L) = DVT(IP1,J,L) + FLUXV
      DVT(I,J,L)   = DVT(I,J,L)   - FLUXV
  210 I = IP1
  215 CONTINUE
C**** Contribution from the South-North mass flux
      DO 230 J=2,JM-1
      DO 220 IP1=1,IM
      FLUX = DT3*(GY(I,J)+GY(IP1,J))
      FLUXU = FLUX*(U(I,J,L)+U(I,J+1,L))
      DUT(I,J+1,L) = DUT(I,J+1,L) + FLUXU
      DUT(I,J,L)   = DUT(I,J,L)   - FLUXU
      FLUXV = FLUX*(V(I,J,L)+V(I,J+1,L))
      DVT(I,J+1,L) = DVT(I,J+1,L) + FLUXV
      DVT(I,J,L)   = DVT(I,J,L)   - FLUXV
C**** Contribution from the Southwest-Northeast mass flux
      FLUX = DT6*(FX(I,J)+GY(IP1,J)-FX1(I,J)-GY1(IP1,J))
      FLUXU = FLUX*(U(IP1,J+1,L)+U(I,J,L))
      DUT(IP1,J+1,L) = DUT(IP1,J+1,L) + FLUXU
      DUT(I,J,L)     = DUT(I,J,L)     - FLUXU
      FLUXV = FLUX*(V(IP1,J+1,L)+V(I,J,L))
      DVT(IP1,J+1,L) = DVT(IP1,J+1,L) + FLUXV
      DVT(I,J,L)     = DVT(I,J,L)     - FLUXV
C**** Contribution from the Southeast-Northwest mass flux
      FLUX = DT6*(-FX(I,J)+GY(IP1,J)+FX1(I,J)-GY1(IP1,J))
      FLUXU = FLUX*(U(I,J+1,L)+U(IP1,J,L))
      DUT(I,J+1,L) = DUT(I,J+1,L) + FLUXU
      DUT(IP1,J,L) = DUT(IP1,J,L) - FLUXU
      FLUXV = FLUX*(V(I,J+1,L)+V(IP1,J,L))
      DVT(I,J+1,L) = DVT(I,J+1,L) + FLUXV
      DVT(IP1,J,L) = DVT(IP1,J,L) - FLUXV
  220 I=IP1
  230 CONTINUE
      I=IM
      IM1=I-1
C**** CONTRIBUTION FROM THE BIG STEP WEST-EAST MASS FLUX
      DO 245 J=2,JM
      DO 240 IP1=1,IM
      FLUX=-DT24*(PU(I,J,L)+PU(I,J-1,L))
      FLUXU=FLUX*(U(IP1,J,L)+U(IM1,J,L))
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU
      DUT(IM1,J,L)=DUT(IM1,J,L)-FLUXU
      FLUXV=FLUX*(V(IP1,J,L)+V(IM1,J,L))
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUXV
      IM1=I
  240 I=IP1
  245 CONTINUE
C**** CONTRIBUTION FROM THE BIG STEP SOUTH-NORTH MASS FLUX
      DO 260 J=3,JM-1
      DO 250 IP1=1,IM
      FLUX=-DT24*(PV(I,J,L)+PV(IP1,J,L))
      FLUXU=FLUX*(U(I,J-1,L)+U(I,J+1,L))
      DUT(I,J+1,L)=DUT(I,J+1,L)+FLUXU
      DUT(I,J-1,L)=DUT(I,J-1,L)-FLUXU
      FLUXV=FLUX*(V(I,J-1,L)+V(I,J+1,L))
      DVT(I,J+1,L)=DVT(I,J+1,L)+FLUXV
      DVT(I,J-1,L)=DVT(I,J-1,L)-FLUXV
  250 I=IP1
  260 CONTINUE
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
         IF(MRCH.GT.0) CALL DIAGCD (1,U,V,DUT,DVT,DT1,PIT)
      DO L=1,LM
      DO J=2,JM
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
        VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO
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
      PDT4=DT8*(P(I,J-1,L)+P(I,J,L))
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
         IF(MRCH.GT.0) CALL DIAGCD (2,U,V,DUT,DVT,DT1)
      DO L=1,LM
      DO J=2,JM
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
        VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF ADVECV
C****
      DO J=2,JM-1
      DO I=1,IM
        FD(I,J)   = PB(I,J)*DXYP(J)
      END DO
      END DO
      FD(:,0)     = 0.
      FD(:,JM+1)  = 0.
      FD(:,1)     = 2.*PB(1,1)*DXYP(1)
      FD(:,JM)    = 2.*PB(1,JM)*DXYP(JM)

      I=IM
      DO J=1,JM+1
      DO IP1=1,IM
        FDU1(I,J) = .25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
        I=IP1
      END DO
      END DO

      I=IM
      IM1=I-1
      DO J=2,JM
      DO IP1=1,IM
        FDU2(I,J)=.25*(FDU1(IP1,J)+FDU1(I,J+1)+FDU1(I,J-1)+FDU1(IM1,J))
        IM1=I
        I=IP1
      END DO
      END DO

      DO J=2,JM
      DO I=1,IM
        FDU2(I,J)=(5.*FDU1(I,J)-2.*FDU2(I,J))*by3
      END DO
      END DO

      DO L=1,LS1-1
      DO J=2,JM
      DO I=1,IM
        FDU(I,J,L)=FDU2(I,J)*DSIG(L)
      END DO
      END DO
      END DO

      DUT=0.
      DVT=0.
      DO L=1,LM
      DO J=2,JM
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)/FDU(I,J,L)
        VT(I,J,L)=VT(I,J,L)/FDU(I,J,L)
      END DO
      END DO
      END DO
      RETURN
      END SUBROUTINE ADVECV
