C**** QUSEM12 E001M12 SOMTQ QUSB261AM12
C**** cannot use OPT(3)
C**** QUSBM9=QUSB140M9 with correction to second order moment calc.
C**** Changes for constant pressure above LS1 + double precision
C**** QUS   is Russell quadratic upstream scheme for temperature
C**** and water vapor advection, with limits applied to water vapor.
C**** Changes for constant pressure above LS1
C**** FQU,FQV for additional diagnostics
C**** Routines included: AADVT, AADVTX, AADVTY, AADVTZ
      SUBROUTINE AADVT (PC,P,RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RZXM,
     *                  RYZM,  DXYP,DSIG,PSMPTP,LS1,DT,QLIMIT,FQU,FQV)
C****
C**** AADVT advects tracers using the Quadradic Upstream Scheme.
C****
C**** Input: MA (kg) = fluid mass before advection
C****         DT (s) = time step
C****      MU (kg/s) = west to east mass flux
C****      MV (kg/s) = south to north mass flux
C****      MW (kg/s) = vertical mass flux
C****         QLIMIT = whether moment limitations should be used
C****
C**** Output:    RM (kg) = mean tracer mass
C****   RXM,RYM,RZM (kg) = first moments of tracer mass
C****   RXXM,RYYM,. (kg) = second moments of tracer mass
C****
      USE E001M12_COM
     &   , ONLY : IM,JM,LM
      IMPLICIT REAL*8 (A-H,M,O-Z)
      LOGICAL*4 QLIMIT
      REAL*8 RXM(IM,JM,LM), RYM(IM,JM,LM), RZM(IM,JM,LM),
     *      RXXM(IM,JM,LM),RYYM(IM,JM,LM),RZZM(IM,JM,LM),
     *      RXYM(IM,JM,LM),RYZM(IM,JM,LM),RZXM(IM,JM,LM)
      REAL*8 RM(IM,JM,LM),PC(IM,JM),DXYP(JM),DSIG(LM)
      COMMON /WORK03/ MA(IM,JM,LM)
      COMMON /FLUXCB/ MU(IM,JM,LM),MV(IM,JM,LM),MW(IM,JM,LM-1)
      COMMON /WORK1/PIT(IM,JM),SD(IM,JM,LM-1),PU(IM,JM,LM),
     *        PV(IM,JM,LM)
      DIMENSION P(IM,JM),FQU(IM,JM),FQV(IM,JM)

C**** Fill in values at the poles
      DO 120 L=1,LM
      DO 120 I=2,IM
        RM(I, 1,L) =   RM(1,1,L)
       RZM(I, 1,L) =  RZM(1,1,L)
      RZZM(I, 1,L) = RZZM(1,1,L)
        RM(I,JM,L) =   RM(1,JM,L)
       RZM(I,JM,L) =  RZM(1,JM,L)
  120 RZZM(I,JM,L) = RZZM(1,JM,L)
C****
C**** Load mass after advection from mass before advection
C****
      DO 200 L=1,LS1-1
      DO 200 J=1,JM
      DO 200 I=1,IM
  200 MA(I,J,L)=PC(I,J)*DSIG(L)*DXYP(J)
      DO 220 L=LS1,LM
      DO 220 J=1,JM
      DO 220 I=1,IM
  220 MA(I,J,L)=PSMPTP*DXYP(J)*DSIG(L)
      DO 300 L=1,LM
      DO 300 J=1,JM
      DO 300 I=1,IM
      RM(I,J,L)=RM(I,J,L)*MA(I,J,L)
      RXM(I,J,L)=RXM(I,J,L)*MA(I,J,L)
      RYM(I,J,L)=RYM(I,J,L)*MA(I,J,L)
      RZM(I,J,L)=RZM(I,J,L)*MA(I,J,L)
      RXXM(I,J,L)=RXXM(I,J,L)*MA(I,J,L)
      RYYM(I,J,L)=RYYM(I,J,L)*MA(I,J,L)
      RZZM(I,J,L)=RZZM(I,J,L)*MA(I,J,L)
      RXYM(I,J,L)=RXYM(I,J,L)*MA(I,J,L)
      RYZM(I,J,L)=RYZM(I,J,L)*MA(I,J,L)
      RZXM(I,J,L)=RZXM(I,J,L)*MA(I,J,L)
  300 MU(I,J,L)=PU(I,J,L)
      DO 310 L=1,LM
      DO 310 J=1,JM-1
      DO 310 I=1,IM
  310 MV(I,J,L)=PV(I,J+1,L)
      DO 320 L=1,LM-1
      DO 320 J=1,JM
      DO 320 I=1,IM
  320 MW(I,J,L)=-SD(I,J,L)
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
      CALL AADVTX (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,MA,
     *  .5*DT,QLIMIT,FQU)
      CALL AADVTY (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,MA,
     *     DT,QLIMIT,FQV)
      CALL AADVTZ (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,MA,
     *     DT,QLIMIT)
      CALL AADVTX (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,MA,
     *  .5*DT,QLIMIT,FQU)
C**** Fill in values at the poles
      DO 210 L=1,LM
      DO 210 I=1,IM
        MA(I, 1,L) =   MA(IM,1,L)
        RM(I, 1,L) =   RM(IM,1,L)
       RZM(I, 1,L) =  RZM(IM,1,L)
      RZZM(I, 1,L) = RZZM(IM,1,L)
        MA(I,JM,L) =   MA(1,JM,L)
        RM(I,JM,L) =   RM(1,JM,L)
       RZM(I,JM,L) =  RZM(1,JM,L)
  210 RZZM(I,JM,L) = RZZM(1,JM,L)
      DO 330 L=1,LM
      DO 330 J=1,JM
      DO 330 I=1,IM
      BYMA = 1.D0/MA(I,J,L)
      RM(I,J,L)=RM(I,J,L)    *BYMA
      RXM(I,J,L)=RXM(I,J,L)  *BYMA
      RYM(I,J,L)=RYM(I,J,L)  *BYMA
      RZM(I,J,L)=RZM(I,J,L)  *BYMA
      RXXM(I,J,L)=RXXM(I,J,L)*BYMA
      RYYM(I,J,L)=RYYM(I,J,L)*BYMA
      RZZM(I,J,L)=RZZM(I,J,L)*BYMA
      RXYM(I,J,L)=RXYM(I,J,L)*BYMA
      RYZM(I,J,L)=RYZM(I,J,L)*BYMA
  330 RZXM(I,J,L)=RZXM(I,J,L)*BYMA
C**** CHECK MA(I,J,L)   (debugging only)
c     DO 340 L=1,LS1-1
c     DO 340 J=1,JM
c     DO 340 I=1,IM
c     CHECK=P(I,J)-MA(I,J,L)/DSIG(L)/DXYP(J)
c 340 IF(ABS(CHECK).GT..01)WRITE(6,9000)I,J,L,MA(I,J,L),PC(I,J),DSIG(L)
c    *     ,DXYP(J),CHECK
c     DO 345 L=LS1,LM
c     DO 345 J=1,JM
c     DO 345 I=1,IM
c     CHECK=PSMPTP-MA(I,J,L)/DSIG(L)/DXYP(J)
c 345 IF(ABS(CHECK).GT..01)WRITE(6,9000)I,J,L,MA(I,J,L),PSMPTP,DSIG(L)
c    *     ,DXYP(J),CHECK
 9000 FORMAT(' MA DIFFERS FROM P*DSIG*DXYP BY TOO MUCH',3I3,5E11.3)
      RETURN
      END
      SUBROUTINE AADVTX
     *  (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,M,DT,QLIMIT,FQU)
C****
C**** AADVTX advects tracers in the west to east direction using the
C**** Quadratic Upstream Scheme.  If QLIMIT is true, the moments are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MU (kg/s) = west to east mass flux
C****        QLIMIT = whether moment limitations should be used
C****
C**** Output:       RM (kg) = tracer mass
C****      RXM,RYM,RZM (kg) = first moments of tracer mass
C****   RXXM,RYYM,RZZM (kg) = second moments of tracer mass
C****                M (kg) = fluid mass
C****
      USE E001M12_COM
     &   , ONLY : IM,JM,LM
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 RXM(IM,JM,LM), RYM(IM,JM,LM), RZM(IM,JM,LM),
     *      RXXM(IM,JM,LM),RYYM(IM,JM,LM),RZZM(IM,JM,LM),
     *      RXYM(IM,JM,LM),RYZM(IM,JM,LM),RZXM(IM,JM,LM)
      REAL*8  RM(IM,JM,LM),   M(IM,JM,LM)
      LOGICAL*4 QLIMIT
      COMMON /FLUXCB/ MU(IM,JM,LM)
      COMMON /WORK04/ A(IM),AM(IM), F(IM),FX(IM),FY(IM),
     *  FZ(IM),FXX(IM),FYY(IM),FZZ(IM),FXY(IM),FZX(IM),FYZ(IM)
      DIMENSION FQU(IM,JM)
C**** Loop over layers and latitudes
      DO 420 L=1,LM
      DO 420 J=2,JM-1
C****
C**** Calculate tracer mass flux F (kg)
C****
      I=IM
      DO 120 IP1=1,IM
      AM(I) = DT*MU(I,J,L)
      IF(AM(I).LT.0.)  GO TO 110
      A(I) = AM(I)/M(I,J,L)
      IF(A(I).GT.1.)  WRITE (6,*) 'A>1:',I,J,L,A(I),M(I,J,L)
      F(I) = A(I)*(RM(I,J,L)+(1.-A(I))*(RXM(I,J,L)
     *                   +(1.-2.*A(I))*RXXM(I,J,L)))
      GO TO 120
  110 A(I) = AM(I)/M(IP1,J,L)
      IF(A(I).LT.-1.)  WRITE (6,*) 'A<-1:',I,J,L,A(I),M(IP1,J,L)
      F(I) = A(I)*(RM(IP1,J,L)-(1.+A(I))*(RXM(IP1,J,L)
     *                     -(1.+2.*A(I))*RXXM(IP1,J,L)))
  120 I=IP1
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      IM1=IM
      DO 295 I=1,IM
      IF(A(IM1).GE.0.)  GO TO 280
C**** Air is leaving through the left edge: 2 or 3 divisions
      IF(A(I).LE.0.)  GO TO 260
C**** Air is leaving through the right edge: 3 divisions
      IF(F(IM1).GT.0.)  GO TO 230
C**** Tracer leaving to the left is non-negative: case 1, 3, 4 or 6
      IF(F(I).LT.0.)  GO TO 210
C**** Tracer leaving to the right is non-negative: case 1 or 4
      RMCEN = RM(I,J,L) + (F(IM1)-F(I))
      IF(RMCEN.GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 4
      GAMMA = 1.+(A(IM1)-A(I))
      G13AB = GAMMA*GAMMA - 1. + 3.*(A(IM1)+A(I))**2
      RXXM(I,J,L) = RXXM(I,J,L) - RMCEN*10.*G13AB /
     /  (GAMMA*(12.*(A(IM1)+A(I))**2 + 5.*G13AB*G13AB))
      RXM(I,J,L) = RXM(I,J,L) + RMCEN*12.*(A(IM1)+A(I)) /
     /  (GAMMA*(12.*(A(IM1)+A(I))**2 + 5.*G13AB*G13AB))
      F(IM1) = A(IM1)*(RM(I,J,L)-(1.+A(IM1))*(RXM(I,J,L)
     -                       -(1.+2.*A(IM1))*RXXM(I,J,L)))
      IF(F(IM1).GT.0.)  GO TO 240
      IF(RM(I,J,L)+F(IM1).LT.0.)  GO TO 220
      F(I) = RM(I,J,L)+F(IM1)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3 or 6
  210 IF(RM(I,J,L).LT.F(I)-F(IM1))  GO TO 220
C**** Only the tracer leaving to the right is negative: case 3
      RXXM(I,J,L) = RXXM(I,J,L) - F(I)*(1.-2.*A(I)) /
     /         (A(I)*(1.-A(I))*(.6d0+(1.-2.*A(I))**2))
      RXM(I,J,L) = RM(I,J,L)/(A(I)-1.) - (1.-2.*A(I))*RXXM(I,J,L)
      F(IM1) = A(IM1)*(RM(I,J,L)-(1.+A(IM1))*(RXM(I,J,L)
     -                       -(1.+2.*A(IM1))*RXXM(I,J,L)))
      IF(F(IM1).GT.0.)  GO TO 250
      F(I) = 0.
      IF(RM(I,J,L)+F(IM1).GE.0.)  GO TO 295
C**** The two right divisions are negative: case 6
  220 RXXM(I,J,L) = RM(I,J,L)/(2.*A(IM1)*(A(I)-1.))
      RXM (I,J,L) = RM(I,J,L)*(.5-A(IM1)-A(I))/(A(IM1)*(1.-A(I)))
      F(IM1) = -RM(I,J,L)
      F(I)   = 0.
      GO TO 295
C**** Tracer leaving to the left is negative: case 2, 5 or 7
  230 IF(F(I).LT.0.)  GO TO 250
C**** Tracer leaving to the right is non-negative: case 2 or 5
      IF(RM(I,J,L).LT.F(I)-F(IM1))  GO TO 240
C**** Only the tracer leaving to the left is negative: case 2
      RXXM(I,J,L) = RXXM(I,J,L) - F(IM1)*(1.+2.*A(IM1)) /
     /  (A(IM1)*(1.+A(IM1))*(.6d0+(1.+2.*A(IM1))**2))
      RXM(I,J,L) = RM(I,J,L)/(1.+A(IM1)) + (1.+2.*A(IM1))*RXXM(I,J,L)
      F(I) = A(I)*(RM(I,J,L)+(1.-A(I))*(RXM(I,J,L)
     +                   +(1.-2.*A(I))*RXXM(I,J,L)))
      IF(F(I).LT.0.)  GO TO 250
      F(IM1) = 0.
      IF(RM(I,J,L)-F(I).GE.0.)  GO TO 295
C**** The two left divisions are negative: case 5
  240 RXXM(I,J,L) = RM(I,J,L)/(2.*A(I)*(1.+A(IM1)))
      RXM (I,J,L) = RM(I,J,L)*(.5+A(IM1)+A(I))/(A(I)*(1.+A(IM1)))
      F(IM1) = 0.
      F(I)   = RM(I,J,L)
      GO TO 295
C**** Tracer leaving to both the left and right is negative: case 7
  250 GAMMA = 1.+(A(IM1)-A(I))
      RXXM(I,J,L) = -RM(I,J,L)*(1.+GAMMA)/
     /  (2.*GAMMA*(1.+A(IM1))*(1.-A(I)))
      RXM (I,J,L) = -RM(I,J,L)*(A(IM1)+A(I))*(1.+2.*GAMMA)/
     /  (2.*GAMMA*(1.+A(IM1))*(1.-A(I)))
      F(IM1) = 0.
      F(I  ) = 0.
      GO TO 295
C**** No air is leaving through the right edge: 2 divisions
  260 IF(F(IM1).GT.0.)  GO TO 270
C**** Tracer leaving to the left is non-negative: case 1 or 3
      IF(RM(I,J,L)+F(IM1).GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 3
      IF(A(IM1).LT.-1.)  STOP ' AADVTX 260'
      RXXM(I,J,L) = RXXM(I,J,L)-(RM(I,J,L)+F(IM1))*(1.+2.*A(IM1))/
     /        (A(IM1)*(1.+A(IM1))*(.6d0+(1.+2.*A(IM1))**2))
      RXM(I,J,L) = RM(I,J,L)/A(IM1) + (1.+2.*A(IM1))*RXXM(I,J,L)
      F(IM1) = -RM(I,J,L)
      GO TO 295
C**** Tracer leaving to the left is negative: case 2
  270 IF(A(IM1).LT.-1.)  STOP ' AADVTX 270'
      RXXM(I,J,L) = RXXM(I,J,L) - F(IM1)*(1.+2.*A(IM1)) /
     /           (A(IM1)*(1.+A(IM1))*(.6d0+(1.+2.*A(IM1))**2))
      RXM(I,J,L) = RM(I,J,L)/(1.+A(IM1))+(1.+2.*A(IM1))*RXXM(I,J,L)
      F(IM1) = 0.
      GO TO 295
C**** No air is leaving through the left edge: 1 or 2 divisions
  280 IF(A(I).LE.0.)  GO TO 295
C**** Air is leaving through the right edge: 2 divisions
      IF(F(I).LT.0.)  GO TO 290
C**** Tracer leaving to the right is non-negative: case 1 or 2
      IF(RM(I,J,L)-F(I).GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 2
      IF(A(I).GT.1.)  STOP ' AADVTX 280'
      RXXM(I,J,L) = RXXM(I,J,L)+(RM(I,J,L)-F(I))*(1.-2.*A(I))/
     /        (A(I)*(1.-A(I))*(.6d0+(1.-2.*A(I))**2))
      RXM(I,J,L) = RM(I,J,L)/A(I)-(1.-2.*A(I))*RXXM(I,J,L)
      F(I) = RM(I,J,L)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3
  290 IF(A(I).GT.1.)  STOP ' AADVTX 290'
      RXXM(I,J,L) = RXXM(I,J,L) - F(I)*(1.-2.*A(I)) /
     *             (A(I)*(1.-A(I))*(.6d0+(1.-2.*A(I))**2))
      RXM(I,J,L) = RM(I,J,L)/(A(I)-1.)-(1.-2.*A(I))*RXXM(I,J,L)
      F(I) = 0.
  295 IM1=I
C****
C**** Calculate FX (kg**2), FY (kg), and FZ (kg)
C****
  300 I=IM
      DO 320 IP1=1,IM
      IF(AM(I).LT.0.)  GO TO 310
C**** Air mass flux is positive
      FYY(I) = A(I)*RYYM(I,J,L)
      FZZ(I) = A(I)*RZZM(I,J,L)
      FYZ(I) = A(I)*RYZM(I,J,L)
      FY(I)  = A(I)*(RYM(I,J,L)+(1.-A(I))*RXYM(I,J,L))
      FZ(I)  = A(I)*(RZM(I,J,L)+(1.-A(I))*RZXM(I,J,L))
      FXY(I) = AM(I)*(A(I)*A(I)*RXYM(I,J,L)-3.*FY(I))
      FZX(I) = AM(I)*(A(I)*A(I)*RZXM(I,J,L)-3.*FZ(I))
      FX(I)  = AM(I)*(A(I)*A(I)*(RXM(I,J,L)
     *            +3.*(1.-A(I))*RXXM(I,J,L))-3.*F(I))
      FXX(I) = AM(I)*(AM(I)*A(I)**3*RXXM(I,J,L)-5.*(AM(I)*F(I)+FX(I)))
      GO TO 320
C**** Air mass flux is negative
  310 FYY(I) = A(I)*RYYM(IP1,J,L)
      FZZ(I) = A(I)*RZZM(IP1,J,L)
      FYZ(I) = A(I)*RYZM(IP1,J,L)
      FY(I)  = A(I)*(RYM(IP1,J,L)-(1.+A(I))*RXYM(IP1,J,L))
      FZ(I)  = A(I)*(RZM(IP1,J,L)-(1.+A(I))*RZXM(IP1,J,L))
      FXY(I) = AM(I)*(A(I)*A(I)*RXYM(IP1,J,L)-3.*FY(I))
      FZX(I) = AM(I)*(A(I)*A(I)*RZXM(IP1,J,L)-3.*FZ(I))
      FX(I)  = AM(I)*(A(I)*A(I)*(RXM(IP1,J,L)
     *            -3.*(1.+A(I))*RXXM(IP1,J,L))-3.*F(I))
      FXX(I) = AM(I)*(AM(I)*A(I)**3*RXXM(IP1,J,L)-5.*(AM(I)*F(I)+FX(I)))
  320 I=IP1
C****
C**** Calculate new tracer mass and moments of tracer mass
C**** Calculate new air mass distribution
C****
      IM1 = IM
      DO 410 I=1,IM
      MNEW = M(I,J,L) + AM(IM1)-AM(I)
      BYMNEW = 1./MNEW
      MOT2 = -(AM(IM1)+AM(I))
      RYYM(I,J,L) = RYYM(I,J,L) + FYY(IM1)-FYY(I)
      RZZM(I,J,L) = RZZM(I,J,L) + FZZ(IM1)-FZZ(I)
      RYZM(I,J,L) = RYZM(I,J,L) + FYZ(IM1)-FYZ(I)
      RYM (I,J,L) = RYM (I,J,L) +  FY(IM1)- FY(I)
      RZM (I,J,L) = RZM (I,J,L) +  FZ(IM1)- FZ(I)
      RMFIM1 = RM(I,J,L)+F(IM1)      ! 7/94 These lines fix
      RM  (I,J,L) = RMFIM1-F(I)      !  truncation error if opt2
C     RM  (I,J,L) = RM  (I,J,L) +   F(IM1)-  F(I)
      RXYM(I,J,L) = (RXYM(I,J,L)*M(I,J,L)-3.*(MOT2*RYM(I,J,L) +
     +   M(I,J,L)*(FY(IM1)+FY(I)))+(FXY(IM1)-FXY(I)))*BYMNEW
      RZXM(I,J,L) = (RZXM(I,J,L)*M(I,J,L)-3.*(MOT2*RZM(I,J,L) +
     +   M(I,J,L)*(FZ(IM1)+FZ(I)))+(FZX(IM1)-FZX(I)))*BYMNEW
      RXM (I,J,L) = (RXM (I,J,L)*M(I,J,L)-3.*(MOT2* RM(I,J,L) +
     +   M(I,J,L)*( F(IM1)+ F(I)))+( FX(IM1)- FX(I)))*BYMNEW
      RXXM(I,J,L) = (RXXM(I,J,L)*M(I,J,L)*M(I,J,L) +
     +  2.5*RM(I,J,L)*(M(I,J,L)*M(I,J,L)-MNEW*MNEW -
     -  3.*MOT2*MOT2)+5.*(M(I,J,L)*(M(I,J,L)*(F(IM1)-F(I)) -
     -  FX(IM1)-FX(I)) - MNEW*MOT2*RXM(I,J,L)) + (FXX(IM1)-FXX(I)))*
     *  (BYMNEW*BYMNEW)
C**** Store new air mass in M array
      M(I,J,L) = MNEW
         IF(M(I,J,L).LT.0.)  GO TO 800
         IF(QLIMIT .AND. RM(I,J,L).LT.0.) GO TO 810
C**** COLLECT EAST-WEST TRACER FLUX
      FQU(I,J)=FQU(I,J)+F(I)
  410 IM1=I
  420 CONTINUE
      RETURN
C****
  800 WRITE (6,*) 'M<0 IN AADVTX:',I,J,L,M(I,J,L)
  810 WRITE (6,*) 'RM IN AADVTX:',I,J,L,RM(I,J,L),RMFIM1
      WRITE (6,*) 'A=',(I,A(I),I=1,IM)
      STOP
      END
      SUBROUTINE AADVTY
     *  (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,M,DT,QLIMIT,FQV)
C****
C**** AADVTY advects tracers in the south to north direction using the
C**** Quadratic Upstream Scheme.  If QLIMIT is true, the moments are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MV (kg/s) = south to north mass flux
C****        QLIMIT = whether slope limitations should be used
C****
C**** Output:       RM (kg) = tracer mass
C****      RXM,RYM,RZM (kg) = first moments of tracer mass
C****   RXXM,RYYM,RZZM (kg) = second moments of tracer mass
C****                M (kg) = fluid mass
C****
      USE E001M12_COM
     &   , ONLY : IM,JM,LM
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 RXM(IM,JM,LM), RYM(IM,JM,LM), RZM(IM,JM,LM),
     *      RXXM(IM,JM,LM),RYYM(IM,JM,LM),RZZM(IM,JM,LM),
     *      RXYM(IM,JM,LM),RYZM(IM,JM,LM),RZXM(IM,JM,LM)
      REAL*8  RM(IM,JM,LM),   M(IM,JM,LM)
      LOGICAL*4 QLIMIT
      COMMON /FLUXCB/ MU(IM,JM,LM),MV(IM,JM,LM)
      COMMON /WORK04/ B(JM),BM(JM), F(JM),FX(JM),FY(JM),
     *  FZ(JM),FXX(JM),FYY(JM),FZZ(JM),FXY(JM),FZX(JM),FYZ(JM)
      DIMENSION FQV(IM,JM)
      BYIM = 1./IM
C**** Loop over layers and longitudes
      DO 440 L=1,LM
      SBMS  = 0.
      SFS   = 0.
      SFZS  = 0.
      SFZZS = 0.
      SBMN  = 0.
      SFN   = 0.
      SFZN  = 0.
      SFZZN = 0.
      DO 430 I=1,IM
C****
C**** Calculate tracer mass flux F (kg)
C****
C**** Near the South Pole
      BM(1) = DT*MV(I,1,L)
      IF(BM(1).LT.0.)  GO TO 110
      B(1) = BM(1)/M(IM,1,L)
C     IF(B(1).GT.1.)  WRITE (6,*) 'B>1:',I,1,L,B(1),M(IM,1,L)
      F(1) = B(1)*RM(IM,1,L)
      GO TO 120
  110 B(1) = BM(1)/M(I,2,L)
      IF(B(1).LT.-1.)  WRITE (6,*) 'B<-1:',I,1,L,B(1),M(1,2,L)
      F(1) = B(1)*(RM(I,2,L)-(1.+B(1))*(RYM(I,2,L)
     *                   -(1.+2.*B(1))*RYYM(I,2,L)))
C**** In the interior
  120 DO 140 J=2,JM-2
      BM(J) = DT*MV(I,J,L)
      IF(BM(J).LT.0.)  GO TO 130
      B(J) = BM(J)/M(I,J,L)
      IF(B(J).GT.1.)  WRITE (6,*) 'B>1:',I,J,L,B(J),M(I,J,L)
      F(J) =  B(J)*(RM(I,J,L)+(1.-B(J))*(RYM(I,J,L)
     *                    +(1.-2.*B(J))*RYYM(I,J,L)))
      GO TO 140
  130 B(J) = BM(J)/M(I,J+1,L)
      IF(B(J).LT.-1.)  WRITE (6,*) 'B<-1:',I,J,L,B(J),M(I,J+1,L)
      F(J) = B(J)*(RM(I,J+1,L)-(1.+B(J))*(RYM(I,J+1,L)
     *                     -(1.+2.*B(J))*RYYM(I,J+1,L)))
  140 CONTINUE
C**** Near the North Pole
      J=JM-1
      BM(J) = DT*MV(I,J,L)
      IF(BM(J).LT.0.)  GO TO 150
      B(J) = BM(J)/M(I,J,L)
      IF(B(J).GT.1.)  WRITE (6,*) 'B>1:',I,J,L,B(J),M(I,J,L)
      F(J) = B(J)*(RM(I,J,L)+(1.-B(J))*(RYM(I,J,L)
     *                   +(1.-2.*B(J))*RYYM(I,J,L)))
      GO TO 200
  150 B(J) = BM(J)/M(1,JM,L)
C     IF(B(J).LT.-1.)  WRITE (6,*) 'B<-1:',I,J,L,B(J),M(1,J+1,L)
      F(J) = B(J)*RM(1,JM,L)
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
  200 IF(.NOT.QLIMIT)  GO TO 300
      DO 295 J=2,JM-1
      IF(B(J-1).GE.0.)  GO TO 280
C**** Air is leaving through the left edge: 2 or 3 divisions
      IF(B(J).LE.0.)  GO TO 260
C**** Air is leaving through the right edge: 3 divisions
      IF(F(J-1).GT.0.)  GO TO 230
C**** Tracer leaving to the left is non-negative: case 1, 3, 4 or 6
      IF(F(J).LT.0.)  GO TO 210
C**** Tracer leaving to the right is non-negative: case 1 or 4
      RMCEN = RM(I,J,L) + (F(J-1)-F(J))
      IF(RMCEN.GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 4
      GAMMA = 1.+(B(J-1)-B(J))
      G13AB = GAMMA*GAMMA - 1. + 3.*(B(J-1)+B(J))**2
      RYYM(I,J,L) = RYYM(I,J,L) - RMCEN*10.*G13AB /
     /  (GAMMA*(12.*(B(J-1)+B(J))**2 + 5.*G13AB*G13AB))
      RYM(I,J,L) = RYM(I,J,L) + RMCEN*12.*(B(J-1)+B(J)) /
     /  (GAMMA*(12.*(B(J-1)+B(J))**2 + 5.*G13AB*G13AB))
      F(J-1) = B(J-1)*(RM(I,J,L)-(1.+B(J-1))*(RYM(I,J,L)
     *                       -(1.+2.*B(J-1))*RYYM(I,J,L)))
      IF(F(J-1).GT.0.)  GO TO 240
      IF(RM(I,J,L)+F(J-1).LT.0.)  GO TO 220
      F(J) = RM(I,J,L)+F(J-1)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3 or 6
  210 IF(RM(I,J,L).LT.F(J)-F(J-1))  GO TO 220
C**** Only the tracer leaving to the right is negative: case 3
      RYYM(I,J,L) = RYYM(I,J,L) - F(J)*(1.-2.*B(J)) /
     /             (B(J)*(1.-B(J))*(.6d0+(1.-2.*B(J))**2))
      RYM(I,J,L) = RM(I,J,L)/(B(J)-1.)-(1.-2.*B(J))*RYYM(I,J,L)
      F(J-1) = B(J-1)*(RM(I,J,L)-(1.+B(J-1))*(RYM(I,J,L)
     -                       -(1.+2.*B(J-1))*RYYM(I,J,L)))
      IF(F(J-1).GT.0.)  GO TO 250
      F(J)  = 0.
      IF(RM(I,J,L)+F(J-1).GE.0.)  GO TO 295
C**** The two right divisions are negative: case 6
  220 RYYM(I,J,L) = RM(I,J,L)/(2.*B(J-1)*(B(J)-1.))
      RYM(I,J,L)  = RM(I,J,L)*(.5-B(J-1)-B(J))/(B(J-1)*(1.-B(J)))
      F(J-1) = -RM(I,J,L)
      F(J)   = 0.
      GO TO 295
C**** Tracer leaving to the left is negative: case 2, 5 or 7
  230 IF(F(J).LT.0.)  GO TO 250
C**** Tracer leaving to the right is non-negative: case 2 or 5
      IF(RM(I,J,L).LT.F(J)-F(J-1))  GO TO 240
C**** Only the tracer leaving to the left is negative: case 2
      RYYM(I,J,L) = RYYM(I,J,L) - F(J-1)*(1.+2.*B(J-1)) /
     /           (B(J-1)*(1.+B(J-1))*(.6d0+(1.+2.*B(J-1))**2))
      RYM(I,J,L) = RM(I,J,L)/(1.+B(J-1))+(1.+2.*B(J-1))*RYYM(I,J,L)
      F(J) = B(J)*(RM(I,J,L)+(1.-B(J))*(RYM(I,J,L)
     +                   +(1.-2.*B(J))*RYYM(I,J,L)))
      IF(F(J).LT.0.)  GO TO 250
      F(J-1) = 0.
      IF(RM(I,J,L)-F(J).GE.0.)  GO TO 295
C**** The two left divisions are negative: case 5
  240 RYYM(I,J,L) = RM(I,J,L)/(2.*B(J)*(1.+B(J-1)))
      RYM(I,J,L)  = RM(I,J,L)*(.5+B(J-1)+B(J))/(B(J)*(1.+B(J-1)))
      F(J-1) = 0.
      F(J)   = RM(I,J,L)
      GO TO 295
C**** Tracer leaving to both the left and right is negative: case 7
  250 GAMMA = 1.+(B(J-1)-B(J))
      RYYM(I,J,L) = -RM(I,J,L)*(1.+GAMMA)/
     /  (2.*GAMMA*(1.+B(J-1))*(1.-B(J)))
      RYM (I,J,L) = -RM(I,J,L)*(B(J-1)+B(J))*(1.+2.*GAMMA)/
     /  (2.*GAMMA*(1.+B(J-1))*(1.-B(J)))
      F(J-1) = 0.
      F(J)   = 0.
      GO TO 295
C**** No air is leaving through the right edge: 2 divisions
  260 IF(F(J-1).GT.0.)  GO TO 270
C**** Tracer leaving to the left is non-negative: case 1 or 3
      IF(RM(I,J,L)+F(J-1).GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 3
      IF(B(J-1).LT.-1.)  STOP ' AADVTY 260'
      RYYM(I,J,L) = RYYM(I,J,L) - (RM(I,J,L)+F(J-1)) *
     *  (1.+2.*B(J-1))/ (B(J-1)*(1.+B(J-1))*(.6d0+(1.+2.*B(J-1))**2))
      RYM(I,J,L) = RM(I,J,L)/B(J-1) + (1.+2.*B(J-1))*RYYM(I,J,L)
      F(J-1) = -RM(I,J,L)
      GO TO 295
C**** Tracer leaving to the left is negative: case 2
  270 CONTINUE
      IF(B(J-1).LT.-1.)  STOP ' AADVTY 270'
      RYYM(I,J,L) = RYYM(I,J,L) - F(J-1)*(1.+2.*B(J-1)) /
     /  (B(J-1)*(1.+B(J-1))*(.6d0+(1.+2.*B(J-1))**2))
      RYM(I,J,L) = RM(I,J,L)/(1.+B(J-1)) + (1.+2.*B(J-1))*RYYM(I,J,L)
      F(J-1) = 0.
      GO TO 295
C**** No air is leaving through the left edge: 1 or 2 divisions
  280 IF(B(J).LE.0.)  GO TO 295
C**** Air is leaving through the right edge: 2 divisions
      IF(F(J).LT.0.)  GO TO 290
C**** Tracer leaving to the right is non-negative: case 1 or 2
      IF(RM(I,J,L)-F(J).GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 2
      IF(B(J).GT.1.)  STOP ' AADVTY 280'
      RYYM(I,J,L) = RYYM(I,J,L)+(RM(I,J,L)-F(J))*(1.-2.*B(J))/
     /        (B(J)*(1.-B(J))*(.6d0+(1.-2.*B(J))**2))
      RYM(I,J,L) = RM(I,J,L)/B(J)-(1.-2.*B(J))*RYYM(I,J,L)
      F(J) = RM(I,J,L)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3
  290 IF(B(J).GT.1.)  STOP ' AADVTY 290'
      RYYM(I,J,L) = RYYM(I,J,L) - F(J)*(1.-2.*B(J)) /
     *             (B(J)*(1.-B(J))*(.6d0+(1.-2.*B(J))**2))
      RYM(I,J,L) = RM(I,J,L)/(B(J)-1.)-(1.-2.*B(J))*RYYM(I,J,L)
      F(J) = 0.
  295 CONTINUE
C****
C**** Calculate F? (kg), FY? (kg**2), and FYY (kg**3)
C****
C**** Near the South Pole
  300 IF(BM(1).LT.0.)  GO TO 310
      FXX(1) = 0.
      FZZ(1) = B(1)*RZZM(IM,1,L)
      FZX(1) = 0.
      FX(1)  = 0.
      FZ(1)  = B(1)*RZM(IM,1,L)
      FXY(1) = 0.
      FYZ(1) = -3.*BM(1)*FZ(1)
      FY(1)  = -3.*BM(1)*F (1)
      FYY(1) = -5.*BM(1)*(BM(1)*F(1)+FY(1))
      GO TO 320
  310 FXX(1) = B(1)*RXXM(I,2,L)
      FZZ(1) = B(1)*RZZM(I,2,L)
      FZX(1) = B(1)*RZXM(I,2,L)
      FX(1)  = B(1)*(RXM(I,2,L)-(1.+B(1))*RXYM(I,2,L))
      FZ(1)  = B(1)*(RZM(I,2,L)-(1.+B(1))*RYZM(I,2,L))
      FXY(1) = BM(1)*(B(1)*B(1)*RXYM(I,2,L)-3.*FX(1))
      FYZ(1) = BM(1)*(B(1)*B(1)*RYZM(I,2,L)-3.*FZ(1))
      FY(1)  = BM(1)*(B(1)*B(1)*(RYM(I,2,L)
     *            -3.*(1.+B(1))*RYYM(I,2,L))-3.*F(1))
      FYY(1) = BM(1)*(BM(1)*B(1)**3*RYYM(I,2,L)-5.*(BM(1)*F(1)+FY(1)))
C**** In the interior
  320 DO 340 J=2,JM-2
      IF(BM(J).LT.0.)  GO TO 330
      FXX(J) = B(J)*RXXM(I,J,L)
      FZZ(J) = B(J)*RZZM(I,J,L)
      FZX(J) = B(J)*RZXM(I,J,L)
      FX(J)  = B(J)*(RXM(I,J,L)+(1.-B(J))*RXYM(I,J,L))
      FZ(J)  = B(J)*(RZM(I,J,L)+(1.-B(J))*RYZM(I,J,L))
      FXY(J) = BM(J)*(B(J)*B(J)*RXYM(I,J,L)-3.*FX(J))
      FYZ(J) = BM(J)*(B(J)*B(J)*RYZM(I,J,L)-3.*FZ(J))
      FY(J)  = BM(J)*(B(J)*B(J)*(RYM(I,J,L)
     *            +3.*(1.-B(J))*RYYM(I,J,L))-3.*F(J))
      FYY(J) = BM(J)*(BM(J)*B(J)**3*RYYM(I,J,L)  -5.*(BM(J)*F(J)+FY(J)))
      GO TO 340
  330 FXX(J) = B(J)*RXXM(I,J+1,L)
      FZZ(J) = B(J)*RZZM(I,J+1,L)
      FZX(J) = B(J)*RZXM(I,J+1,L)
      FX(J)  = B(J)*(RXM(I,J+1,L)-(1.+B(J))*RXYM(I,J+1,L))
      FZ(J)  = B(J)*(RZM(I,J+1,L)-(1.+B(J))*RYZM(I,J+1,L))
      FXY(J) = BM(J)*(B(J)*B(J)*RXYM(I,J+1,L)-3.*FX(J))
      FYZ(J) = BM(J)*(B(J)*B(J)*RYZM(I,J+1,L)-3.*FZ(J))
      FY(J)  = BM(J)*(B(J)*B(J)*(RYM(I,J+1,L)
     *            -3.*(1.+B(J))*RYYM(I,J+1,L))-3.*F(J))
      FYY(J) = BM(J)*(BM(J)*B(J)**3*RYYM(I,J+1,L)-5.*(BM(J)*F(J)+FY(J)))
  340 CONTINUE
C**** Near the North Pole
      J=JM-1
      IF(BM(J).LT.0.) GO TO 350
      FXX(J) = B(J)*RXXM(I,J,L)
      FZZ(J) = B(J)*RZZM(I,J,L)
      FZX(J) = B(J)*RZXM(I,J,L)
      FX(J)  = B(J)*(RXM(I,J,L)+(1.-B(J))*RXYM(I,J,L))
      FZ(J)  = B(J)*(RZM(I,J,L)+(1.-B(J))*RYZM(I,J,L))
      FXY(J) = BM(J)*(B(J)*B(J)*RXYM(I,J,L)-3.*FX(J))
      FYZ(J) = BM(J)*(B(J)*B(J)*RYZM(I,J,L)-3.*FZ(J))
      FY(J)  = BM(J)*(B(J)*B(J)*(RYM(I,J,L)
     *            +3.*(1.-B(J))*RYYM(I,J,L))-3.*F(J))
      FYY(J) = BM(J)*(BM(J)*B(J)**3*RYYM(I,J,L)  -5.*(BM(J)*F(J)+FY(J)))
      GO TO 360
  350 FXX(J) = 0.
      FZZ(J) = B(J)*RZZM(1,JM,L)
      FZX(J) = 0.
      FX(J)  = 0.
      FZ(J)  = B(J)*RZM(1,JM,L)
      FXY(J) = 0.
      FYZ(J) = -3.*BM(J)*FZ(J)
      FY(J)  = -3.*BM(J)*F(J)
      FYY(J) = -5.*BM(J)*(BM(J)*F(J)+FY(J))
C****
  360 SBMS  = SBMS  + BM(1)
      SFS   = SFS   + F(1)
      SFZS  = SFZS  + FZ(1)
      SFZZS = SFZZS + FZZ(1)
      SBMN  = SBMN  + BM(JM-1)
      SFN   = SFN   + F(JM-1)
      SFZN  = SFZN  + FZ(JM-1)
      SFZZN = SFZZN + FZZ(JM-1)
C****
C**** Calculate new tracer mass and moments of tracer mass
C****
C**** Non-polar regions
      DO 410 J=2,JM-1
      MNEW = M(I,J,L) + (BM(J-1)-BM(J))
      BYMNEW = 1./MNEW
      MOT2 = -(BM(J-1)+BM(J))
      RXXM(I,J,L) = RXXM(I,J,L)+(FXX(J-1)-FXX(J))
      RZZM(I,J,L) = RZZM(I,J,L)+(FZZ(J-1)-FZZ(J))
      RZXM(I,J,L) = RZXM(I,J,L)+(FZX(J-1)-FZX(J))
      RXM (I,J,L) = RXM (I,J,L)+ (FX(J-1)-FX(J))
      RZM (I,J,L) = RZM (I,J,L)+ (FZ(J-1)-FZ(J))
      RMFJM1 = RM(I,J,L)+F(J-1)      ! 7/94 These lines fix
      RM  (I,J,L) = RMFJM1-F(J)      !  truncation error if opt2
C     RM  (I,J,L) = RM  (I,J,L)+  (F(J-1)-F(J))
      RXYM(I,J,L) = (RXYM(I,J,L)*M(I,J,L)-3.*(MOT2*RXM(I,J,L) +
     +   M(I,J,L)*(FX(J-1)+FX(J)))+(FXY(J-1)-FXY(J)))*BYMNEW
      RYZM(I,J,L) = (RYZM(I,J,L)*M(I,J,L)-3.*(MOT2*RZM(I,J,L) +
     +   M(I,J,L)*(FZ(J-1)+FZ(J)))+(FYZ(J-1)-FYZ(J)))*BYMNEW
      RYM (I,J,L)  = (RY M(I,J,L)*M(I,J,L)-3.*(MOT2*R M(I,J,L) +
     +   M(I,J,L)*( F(J-1)+ F(J)))+( FY(J-1)- FY(J)))*BYMNEW
      RYYM(I,J,L) = (RYYM(I,J,L)*M(I,J,L)*M(I,J,L) +
     +  2.5*RM(I,J,L)*(M(I,J,L)*M(I,J,L)-MNEW*MNEW -
     -  3.*MOT2*MOT2)+5.*(M(I,J,L)*(M(I,J,L)*(F(J-1)-F(J)) -
     -  FY(J-1)-FY(J)) - MNEW*MOT2*RYM(I,J,L)) + (FYY(J-1)-FYY(J)))*
     *  (BYMNEW*BYMNEW)
      M(I,J,L) = MNEW
         IF(M(I,J,L).LT.0.)  GO TO 800
         IF(QLIMIT .AND. RM(I,J,L).LT.0.)  GO TO 810
C**** COLLECT NORTH-SOUTH TRACER FLUX
      FQV(I,J)=FQV(I,J)+F(J)
  410 CONTINUE
      FQV(I,1)=FQV(I,1)+F(1)
      FQV(I,JM)=FQV(I,JM)+F(JM)
  430 CONTINUE
C**** At the South Pole
      RZZM(IM,1,L) = RZZM(IM,1,L) - SFZZS*BYIM
      RZM (IM,1,L) =  RZM(IM,1,L) - SFZS*BYIM
      RM  (IM,1,L) =   RM(IM,1,L) - SFS*BYIM
      M   (IM,1,L) =    M(IM,1,L) - SBMS*BYIM
         IF(M(IM,1,L).LT.0.) GO TO 800
         IF(QLIMIT .AND. RM(IM,1,L).LT.0.) GO TO 810
C**** At the North Pole
      RZZM(1,JM,L) = RZZM(1,JM,L) + SFZZN*BYIM
      RZM (1,JM,L) =  RZM(1,JM,L) + SFZN*BYIM
      RM  (1,JM,L) =   RM(1,JM,L) + SFN*BYIM
      M   (1,JM,L) =    M(1,JM,L) + SBMN*BYIM
         IF(M(1,JM,L).LT.0.) GO TO 800
         IF(QLIMIT .AND. RM(1,JM,L).LT.0.) GO TO 810
  440 CONTINUE
      RETURN
C****
  800 WRITE (6,*) 'M<0 IN AADVTY:',I,J,L,M(I,J,L)
  810 WRITE (6,*) 'RM IN AADVTY:',I,J,L,RM(I,J,L),RMFJM1
      WRITE (6,*) 'B=',(J,B(J),J=1,JM-1)
      STOP
      END
      SUBROUTINE AADVTZ
     *  (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,M,DT,QLIMIT)
C****
C**** AADVTZ advects tracers in the upward vertical direction using the
C**** Quadratic Upstream Scheme.  If QLIMIT is true, the moments are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MW (kg/s) = vertical mass flux
C****        QLIMIT = whether slope limitations should be used
C****
C**** Output:       RM (kg) = tracer mass
C****      RXM,RYM,RZM (kg) = first moments of tracer mass
C****   RXXM,RYYM,RZZM (kg) = second moments of tracer mass
C****                M (kg) = air mass
C****
      USE E001M12_COM
     &   , ONLY : IM,JM,LM
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 RXM(IM,JM,LM), RYM(IM,JM,LM), RZM(IM,JM,LM),
     *      RXXM(IM,JM,LM),RYYM(IM,JM,LM),RZZM(IM,JM,LM),
     *      RXYM(IM,JM,LM),RYZM(IM,JM,LM),RZXM(IM,JM,LM)
      REAL*8  RM(IM,JM,LM),   M(IM,JM,LM)
      LOGICAL*4 QLIMIT
      COMMON /FLUXCB/ MU(IM,JM,LM),MV(IM,JM,LM),MW(IM,JM,LM-1)
      COMMON /WORK04/ C(0:LM),CM(0:LM), F(0:LM),FX(LM),FY(LM),
     *  FZ(LM),FXX(LM),FYY(LM),FZZ(LM),FXY(LM),FZX(LM),FYZ(LM)
      C(0)   = 0.
      C(LM)  = 0.
      CM(0)  = 0.
      CM(LM) = 0.
      F(0)   = 0.
      F(LM)  = 0.
C**** Loop over latitudes and longitudes
      DO 430 I=IM,IM*(JM-1)+1
C****
C**** Calculate the tracer mass flux F (kg)
C****
      DO 120 L=1,LM-1
      CM(L) = DT*MW(I,1,L)
      IF(CM(L).LT.0.)  GO TO 110
      C(L) = CM(L)/M(I,1,L)
      IF(C(L).GT.1.)  WRITE (6,*) 'C>1:',I,1,L,C(L),M(I,1,L)
      F(L) = C(L)*(RM(I,1,L)+(1.-C(L))*(RZM(I,1,L)
     *                   +(1.-2.*C(L))*RZZM(I,1,L)))
      GO TO 120
  110 C(L) = CM(L)/M(I,1,L+1)
      IF(C(L).LT.-1.)  WRITE (6,*) 'C<-1:',I,1,L,C(L),M(I,1,L+1)
      F(L) = C(L)*(RM(I,1,L+1)-(1.+C(L))*(RZM(I,1,L+1)
     *                     -(1.+2.*C(L))*RZZM(I,1,L+1)))
  120 CONTINUE
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      DO 295 L=1,LM
      IF(C(L-1).GE.0.)  GO TO 280
C**** Air is leaving through the left edge: 2 or 3 divisions
      IF(C(L).LE.0.)  GO TO 260
C**** Air is leaving through the right edge: 3 divisions
      IF(F(L-1).GT.0.)  GO TO 230
C**** Tracer leaving to the left is non-negative: case 1, 3, 4 or 6
      IF(F(L).LT.0.)  GO TO 210
C**** Tracer leaving to the right is non-negative: case 1 or 4
      RMCEN = RM(I,1,L) + (F(L-1)-F(L))
      IF(RMCEN.GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 4
      GAMMA = 1.+(C(L-1)-C(L))
      G13AB = GAMMA*GAMMA - 1. + 3.*(C(L-1)+C(L))**2
      RZZM(I,1,L) = RZZM(I,1,L) - RMCEN*10.*G13AB /
     /  (GAMMA*(12.*(C(L-1)+C(L))**2 + 5.*G13AB*G13AB))
      RZM(I,1,L) = RZM(I,1,L) + RMCEN*12.*(C(L-1)+C(L)) /
     /  (GAMMA*(12.*(C(L-1)+C(L))**2 + 5.*G13AB*G13AB))
      F(L-1) = C(L-1)*(RM(I,1,L)-(1.+C(L-1))*(RZM(I,1,L)
     -                       -(1.+2.*C(L-1))*RZZM(I,1,L)))
      IF(F(L-1).GT.0.)  GO TO 240
      IF(RM(I,1,L)+F(L-1).LT.0.)  GO TO 220
      F(L)  = RM(I,1,L)+F(L-1)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3 or 6
  210 IF(RM(I,1,L).LT.F(L)-F(L-1))  GO TO 220
C**** Only the tracer leaving to the right is negative: case 3
      RZZM(I,1,L) = RZZM(I,1,L) - F(L)*(1.-2.*C(L)) /
     /                 (C(L)*(1.-C(L))*(.6d0+(1.-2.*C(L))**2))
      RZM(I,1,L) = RM(I,1,L)/(C(L)-1.)-(1.-2.*C(L))*RZZM(I,1,L)
      F(L-1) = C(L-1)*(RM(I,1,L)-(1.+C(L-1))*(RZM(I,1,L)
     -                       -(1.+2.*C(L-1))*RZZM(I,1,L)))
      IF(F(L-1).GT.0.)  GO TO 250
      F(L)  = 0.
      IF(RM(I,1,L)+F(L-1).GE.0.)  GO TO 295
C**** The two right divisions are negative: case 6
  220 RZZM(I,1,L) = RM(I,1,L)/(2.*C(L-1)*(C(L)-1.))
      RZM(I,1,L)  = RM(I,1,L)*(.5-C(L-1)-C(L)) / (C(L-1)*(1.-C(L)))
      F(L-1) = -RM(I,1,L)
      F(L)  = 0.
      GO TO 295
C**** Tracer leaving to the left is negative: case 2, 5 or 7
  230 IF(F(L).LT.0.)  GO TO 250
C**** Tracer leaving to the right is non-negative: case 2 or 5
      IF(RM(I,1,L).LT.F(L)-F(L-1))  GO TO 240
C**** Only the tracer leaving to the left is negative: case 2
      RZZM(I,1,L) = RZZM(I,1,L) - F(L-1)*(1.+2.*C(L-1)) /
     /  (C(L-1)*(1.+C(L-1))*(.6d0+(1.+2.*C(L-1))**2))
      RZM(I,1,L) = RM(I,1,L)/(1.+C(L-1)) + (1.+2.*C(L-1))*RZZM(I,1,L)
      F(L)  = C(L)*(RM(I,1,L)+(1.-C(L))*(RZM(I,1,L) +
     +                     (1.-2.*C(L))*RZZM(I,1,L)))
      IF(F(L).LT.0.)  GO TO 250
      F(L-1) = 0.
      IF(RM(I,1,L)-F(L).GE.0.)  GO TO 295
C**** The two left divisions are negative: case 5
  240 RZZM(I,1,L) = RM(I,1,L)/(2.*C(L)*(1.+C(L-1)))
      RZM(I,1,L)  = RM(I,1,L)*(.5+C(L-1)+C(L)) / (C(L)*(1.+C(L-1)))
      F(L-1) = 0.
      F(L)   = RM(I,1,L)
      GO TO 295
C**** Tracer leaving to both the left and right is negative: case 7
  250 GAMMA = 1.+(C(L-1)-C(L))
      RZZM(I,1,L) = -RM(I,1,L)*(1.+GAMMA)/
     /   (2.*GAMMA*(1.+C(L-1))*(1.-C(L)))
      RZ M(I,1,L) = -RM(I,1,L)*(C(L-1)+C(L))*(1.+2.*GAMMA)/
     /   (2.*GAMMA*(1.+C(L-1))*(1.-C(L)))
      F(L-1) = 0.
      F(L)   = 0.
      GO TO 295
C**** No air is leaving through the right edge: 2 divisions
  260 IF(F(L-1).GT.0.)  GO TO 270
C**** Tracer leaving to the left is non-negative: case 1 or 3
      IF(RM(I,1,L)+F(L-1).GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 3
      IF(C(L-1).LT.-1.)  STOP ' AADVTZ 260'
      RZZM(I,1,L) = RZZM(I,1,L) -(RM(I,1,L)+F(L-1))*
     *  (1.+2.*C(L-1)) / (C(L-1)*(1.+C(L-1))*(.6d0+(1.+2.*C(L-1))**2))
      RZM(I,1,L) = RM(I,1,L)/C(L-1) + (1.+2.*C(L-1))*RZZM(I,1,L)
      F(L-1) = -RM(I,1,L)
      GO TO 295
C**** Tracer leaving to the left is negative: case 2
  270 IF(C(L-1).LT.-1.)  STOP ' AADVTZ 270'
      RZZM(I,1,L) = RZZM(I,1,L) - F(L-1)*(1.+2.*C(L-1)) /
     /             (C(L-1)*(1.+C(L-1))*(.6d0+(1.+2.*C(L-1))**2))
      RZM(I,1,L) = RM(I,1,L)/(1.+C(L-1)) + (1.+2.*C(L-1))*RZZM(I,1,L)
      F(L-1) = 0.
      GO TO 295
C**** No air is leaving through the left edge: 1 or 2 divisions
  280 IF(C(L).LE.0.)  GO TO 295
C**** Air is leaving through the right edge: 2 divisions
      IF(F(L).LT.0.)  GO TO 290
C**** Tracer leaving to the right is non-negative: case 1 or 2
      IF(RM(I,1,L)-F(L).GE.0.)  GO TO 295
C**** Tracer remaining in the box is negative: case 2
      IF(C(L).GT.1.)  STOP ' AADVTZ 280'
      RZZM(I,1,L) = RZZM(I,1,L)+(RM(I,1,L)-F(L))*(1.-2.*C(L))/
     /        (C(L)*(1.-C(L))*(.6d0+(1.-2.*C(L))**2))
      RZM(I,1,L) = RM(I,1,L)/C(L)-(1.-2.*C(L))*RZZM(I,1,L)
      F(L) = RM(I,1,L)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3
  290 IF(C(L).GT.1.)  STOP ' AADVTZ 290'
      RZZM(I,1,L) = RZZM(I,1,L) - F(L)*(1.-2.*C(L)) /
     /             (C(L)*(1.-C(L))*(.6d0+(1.-2.*C(L))**2))
      RZM(I,1,L) = RM(I,1,L)/(C(L)-1.)-(1.-2.*C(L))*RZZM(I,1,L)
      F(L) = 0.
  295 CONTINUE
C****
C**** Calculate F? (kg), FZ? (kg**2), and FZZ (kg**3)
C****
  300 DO 320 L=1,LM-1
      IF(CM(L).LT.0.)  GO TO 310
C**** Air mass flux is positive
      FYY(L) = C(L)*RYYM(I,1,L)
      FXX(L) = C(L)*RXXM(I,1,L)
      FXY(L) = C(L)*RXYM(I,1,L)
      FY (L) = C(L)*(RYM(I,1,L)+(1.-C(L))*RYZM(I,1,L))
      FX (L) = C(L)*(RXM(I,1,L)+(1.-C(L))*RZXM(I,1,L))
      FYZ(L) = CM(L)*(C(L)*C(L)*RYZM(I,1,L)-3.*FY(L))
      FZX(L) = CM(L)*(C(L)*C(L)*RZXM(I,1,L)-3.*FX(L))
      FZ (L) = CM(L)*(C(L)*C(L)*(RZM(I,1,L)
     *            +3.*(1.-C(L))*RZZM(I,1,L))-3.*F(L))
      FZZ(L) = CM(L)*(CM(L)*C(L)**3*RZZM(I,1,L)-5.*(CM(L)*F(L)+FZ(L)))
      GO TO 320
C**** Air mass flux is negative
  310 FYY(L) = C(L)*RYYM(I,1,L+1)
      FXX(L) = C(L)*RXXM(I,1,L+1)
      FXY(L) = C(L)*RXYM(I,1,L+1)
      FY (L) = C(L)*(RYM(I,1,L+1)-(1.+C(L))*RYZM(I,1,L+1))
      FX (L) = C(L)*(RXM(I,1,L+1)-(1.+C(L))*RZXM(I,1,L+1))
      FYZ(L) = CM(L)*(C(L)*C(L)*RYZM(I,1,L+1)-3.*FY(L))
      FZX(L) = CM(L)*(C(L)*C(L)*RZXM(I,1,L+1)-3.*FX(L))
      FZ (L) = CM(L)*(C(L)*C(L)*(RZM(I,1,L+1)
     *            -3.*(1.+C(L))*RZZM(I,1,L+1))-3.*F(L))
      FZZ(L) = CM(L)*(CM(L)*C(L)**3*RZZM(I,1,L+1)-5.*(CM(L)*F(L)+FZ(L)))
  320 CONTINUE
C****
C**** Calculate new tracer mass and moments of tracer mass
C**** Calculate new air mass distribution
C****
C**** Calculation in the first layer
      L=1
      MNEW = M(I,1,L)-CM(L)
      BYMNEW = 1./MNEW
      MOT2 = -CM(L)
      RXXM(I,1,L) = RXXM(I,1,L)-FXX(L)
      RYYM(I,1,L) = RYYM(I,1,L)-FYY(L)
      RXYM(I,1,L) = RXYM(I,1,L)-FXY(L)
      RXM (I,1,L) = RXM (I,1,L)-FX(L)
      RYM (I,1,L) = RYM (I,1,L)-FY(L)
      RM  (I,1,L) = RM  (I,1,L)-F(L)
      RZXM(I,1,L) = (RZXM(I,1,L)*M(I,1,L)-3.*(MOT2*RXM(I,1,L) +
     +   M(I,1,L)*FX(L))-FZX(L))*BYMNEW
      RYZM(I,1,L) = (RYZM(I,1,L)*M(I,1,L)-3.*(MOT2*RYM(I,1,L) +
     +   M(I,1,L)*FY(L))-FYZ(L))*BYMNEW
      RZM (I,1,L) = (RZM (I,1,L)*M(I,1,L)-3.*(MOT2*RM(I,1,L) +
     +   M(I,1,L)* F(L))- FZ(L))*BYMNEW
      RZZM(I,1,L) = (RZZM(I,1,L)*M(I,1,L)*M(I,1,L) +
     +  2.5*RM(I,1,L)*(M(I,1,L)*M(I,1,L)-MNEW*MNEW-3.*MOT2*MOT2) -
     -  5.*(M(I,1,L)*(M(I,1,L)*F(L) + FZ(L)) + MNEW*MOT2*RZM(I,1,L)) -
     -  FZZ(L)) * (BYMNEW*BYMNEW)
      M(I,1,L) = MNEW
         IF(M(I,1,L).LE.0.) GO TO 800
         IF(QLIMIT .AND. RM(I,1,L).LT.0.) GO TO 810
C**** Calculation in the interior layers
      DO 410 L=2,LM-1
      MNEW = M(I,1,L)+(CM(L-1)-CM(L))
      BYMNEW = 1./MNEW
      MOT2 = -(CM(L-1)+CM(L))
      RXXM(I,1,L) = RXXM(I,1,L)+(FXX(L-1)-FXX(L))
      RYYM(I,1,L) = RYYM(I,1,L)+(FYY(L-1)-FYY(L))
      RXYM(I,1,L) = RXYM(I,1,L)+(FXY(L-1)-FXY(L))
      RXM (I,1,L) = RXM (I,1,L)+ (FX(L-1)-FX(L))
      RYM (I,1,L) = RYM (I,1,L)+ (FY(L-1)-FY(L))
      RMFLM1 = RM(I,1,L)+F(L-1)      ! 7/94 These lines fix
      RM  (I,1,L) = RMFLM1-F(L)      !  truncation error if opt2
C     RM  (I,1,L) = RM  (I,1,L)+  (F(L-1)-F(L))
      RZXM(I,1,L) = (RZXM(I,1,L)*M(I,1,L)-3.*(MOT2*RXM(I,1,L) +
     +   M(I,1,L)*(FX(L-1)+FX(L)))+(FZX(L-1)-FZX(L)))*BYMNEW
      RYZM(I,1,L) = (RYZM(I,1,L)*M(I,1,L)-3.*(MOT2*RYM(I,1,L) +
     +   M(I,1,L)*(FY(L-1)+FY(L)))+(FYZ(L-1)-FYZ(L)))*BYMNEW
      RZM (I,1,L) = (RZM (I,1,L)*M(I,1,L)-3.*(MOT2*RM (I,1,L) +
     +   M(I,1,L)*( F(L-1)+ F(L)))+( FZ(L-1)- FZ(L)))*BYMNEW
      RZZM(I,1,L) = (RZZM(I,1,L)*M(I,1,L)*M(I,1,L) +
     +  2.5*RM(I,1,L)*(M(I,1,L)*M(I,1,L)-MNEW*MNEW-3.*MOT2*MOT2) +
     +  5.*(M(I,1,L)*(M(I,1,L)*(F(L-1)-F(L)) - FZ(L-1)-FZ(L)) -
     -  MNEW*MOT2*RZM(I,1,L)) + (FZZ(L-1)-FZZ(L))) * (BYMNEW*BYMNEW)
      M(I,1,L) = MNEW
         IF(M(I,1,L).LE.0.)  GO TO 800
         IF(QLIMIT .AND. RM(I,1,L).LT.0.)  GO TO 810
  410 CONTINUE
C**** Calculation in the top layer
      L=LM
      MNEW = M(I,1,L)+CM(L-1)
      BYMNEW = 1./MNEW
      MOT2 = -CM(L-1)
      RXXM(I,1,L) = RXXM(I,1,L)+FXX(L-1)
      RYYM(I,1,L) = RYYM(I,1,L)+FYY(L-1)
      RXYM(I,1,L) = RXYM(I,1,L)+FXY(L-1)
      RXM (I,1,L) = RXM (I,1,L)+FX(L-1)
      RYM (I,1,L) = RYM (I,1,L)+FY(L-1)
      RM  (I,1,L) = RM  (I,1,L)+F (L-1)
      RZXM(I,1,L) = (RZXM(I,1,L)*M(I,1,L)-3.*(MOT2*RXM(I,1,L)
     *  + M(I,1,L)*FX(L-1))+FZX(L-1))*BYMNEW
      RYZM(I,1,L) = (RYZM(I,1,L)*M(I,1,L)-3.*(MOT2*RYM(I,1,L)
     *  + M(I,1,L)*FY(L-1))+FYZ(L-1))*BYMNEW
      RZM (I,1,L) = (RZM (I,1,L)*M(I,1,L)-3.*(MOT2*R M(I,1,L)
     *  + M(I,1,L)* F(L-1))+ FZ(L-1))*BYMNEW
      RZZM(I,1,L) = (RZZM(I,1,L)*M(I,1,L)*M(I,1,L)
     *  + 2.5*RM(I,1,L)*(M(I,1,L)*M(I,1,L)-MNEW*MNEW
     *  - 3.*MOT2*MOT2)+5.*(M(I,1,L)*(M(I,1,L)*F(L-1)
     *  - FZ(L-1)) - MNEW*MOT2*RZM(I,1,L))+FZZ(L-1))*(BYMNEW*BYMNEW)
      M(I,1,L) = MNEW
         IF(M(I,1,L).LE.0.) GO TO 800
         IF(QLIMIT .AND. RM(I,1,L).LT.0.) GO TO 810
  430 CONTINUE
      RETURN
C****
  800 WRITE (6,*) 'M<0 IN AADVTZ:',I,1,L,M(I,1,L),QLIMIT
  810 WRITE (6,*) 'RM IN AADVTZ:',I,1,L,RM(I,1,L),RMFLM1
      WRITE( 6,*) 'C=',(L,C(L),L=0,LM)
      STOP
      END
