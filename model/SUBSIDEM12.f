C**** SUBSIDEM12 E001M12 SOMTQ SUBSIDgM12
C**** copied from SUBSIDgM53
C****      debugging statements added, along with an arg
C****   NQ>0 TRACER NUMBER
C****   NQ=0 Q
C****   NQ=-1 T
C**** Russell quadratic upstream scheme for TRACER advection
C**** This version stops if Q becomes < 0 (should not happen with opt 2)
C**** Routines included: SUBSID (Called by MSTCNV to do subsidence)
      SUBROUTINE SUBSID
     *  (RM,RXM,RYM,RZM,RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,M,QLIMIT,
     *   CMN,LMIN,LMAX,NQ)
C**** SUBSID is a modification of AADVQZ
C**** AADVQZ advects tracers in the upward vertical direction using the
C**** Quadratic Upstream Scheme.  If QLIMIT is true, the moments are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: QLIMIT = whether slope limitations should be used
C****
C**** Output:       RM (kg) = tracer mass
C****      RXM,RYM,RZM (kg) = first moments of tracer mass
C****   RXXM,RYYM,RZZM (kg) = second moments of tracer mass
C****                M (kg) = air mass
C****
      USE E001M12_COM
     &   , ONLY : LM
      IMPLICIT NONE
      INTEGER, PARAMETER :: IM=1,JM=1

      INTEGER, INTENT(IN) :: LMIN,LMAX,NQ
      REAL*8, INTENT(INOUT),DIMENSION(IM,JM,LM) :: RM,RXM,RYM,RZM,
     *      RXXM,RYYM,RZZM,RXYM,RYZM,RZXM,RM,M
      LOGICAL*4 QLIMIT
      REAL*8, DIMENSION(LM) :: FX,FY,FZ,FXX,FYY,FZZ,FXY,FZX,FYZ
      REAL*8, DIMENSION(0:LM) :: C,CM,F
      INTEGER I,L,LL
      REAL*8 MNEW,BYMNEW,MOT2,RMFLM1,G13AB,GAMMA,RMCEN
      REAL*8  CMN(0:LM)

      C(LMIN-1) = 0.
      C(LMAX)   = 0.
      CM(LMIN-1)= 0.
      CM(LMAX)  = 0.
      F(LMIN-1) = 0.
      F(LMAX)   = 0.
      DO L=LMIN,LMAX-1
      CM(L) = -CMN(L)
      END DO
C**** Loop over latitudes and longitudes
      IF(NQ.GE.1) THEN
         DO L=1,LM
          IF(RM(1,1,L).LT.0.) WRITE(*,*) 'RM NEG',(RM(1,1,LL),LL=1,LM)
         END DO
       END IF
C     DO 430 I=IM,IM*(JM-1)+1
      I = 1
C****
C**** Calculate the tracer mass flux F (kg)
C****
      DO 120 L=LMIN,LMAX-1
      IF(CM(L).LT.0.)  GO TO 110
      C(L) = CM(L)/M(I,1,L)
      IF(C(L).GT.1.)  WRITE (6,*) 'C>1:',I,1,L,C(L),M(I,1,L)
      F(L) = C(L)*(RM(I,1,L)+(1.-C(L))*(RZM(I,1,L)
     *                   +(1.-2.*C(L))*RZZM(I,1,L)))
      GO TO 120
  110 C(L) = CM(L)/M(I,1,L+1)
      IF(C(L).LT.-1.) then
         WRITE(6,*) 'C<-1 in SUBSID:',I,1,L,C(L),M(I,1,L+1)
         stop ' C<-1 in SUBSID'
      endif
      F(L) = C(L)*(RM(I,1,L+1)-(1.+C(L))*(RZM(I,1,L+1)
     *                     -(1.+2.*C(L))*RZZM(I,1,L+1)))
  120 CONTINUE
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      DO 295 L=LMIN,LMAX
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
      IF(C(L-1).LT.-1.)  STOP ' SUBSID 260'
      RZZM(I,1,L) = RZZM(I,1,L) -(RM(I,1,L)+F(L-1))*
     *  (1.+2.*C(L-1)) / (C(L-1)*(1.+C(L-1))*(.6d0+(1.+2.*C(L-1))**2))
      RZM(I,1,L) = RM(I,1,L)/C(L-1) + (1.+2.*C(L-1))*RZZM(I,1,L)
      F(L-1) = -RM(I,1,L)
      GO TO 295
C**** Tracer leaving to the left is negative: case 2
  270 IF(C(L-1).LT.-1.)  STOP ' SUBSID 270'
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
      IF(C(L).GT.1.)  STOP ' SUBSID 280'
      RZZM(I,1,L) = RZZM(I,1,L)+(RM(I,1,L)-F(L))*(1.-2.*C(L))/
     /        (C(L)*(1.-C(L))*(.6d0+(1.-2.*C(L))**2))
      RZM(I,1,L) = RM(I,1,L)/C(L)-(1.-2.*C(L))*RZZM(I,1,L)
      F(L) = RM(I,1,L)
      GO TO 295
C**** Tracer leaving to the right is negative: case 3
  290 IF(C(L).GT.1.)  STOP ' SUBSID 290'
      RZZM(I,1,L) = RZZM(I,1,L) - F(L)*(1.-2.*C(L)) /
     /             (C(L)*(1.-C(L))*(.6d0+(1.-2.*C(L))**2))
      RZM(I,1,L) = RM(I,1,L)/(C(L)-1.)-(1.-2.*C(L))*RZZM(I,1,L)
      F(L) = 0.
  295 CONTINUE
C****
C**** Calculate F? (kg), FZ? (kg**2), and FZZ (kg**3)
C****
  300 DO 320 L=LMIN,LMAX-1
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
      L = LMIN
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
      DO 410 L=LMIN+1,LMAX-1
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
      L=LMAX
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
C 430 CONTINUE
      RETURN
C****
  800 WRITE (6,*) 'M<0 IN SUBSID at L:',L,M(I,1,L)
  810 WRITE (6,*) 'RM IN SUBSID at L:',L,RM(I,1,L),RMFLM1
      WRITE (6,*) 'C=',(LL,C(LL),LL=0,LM)
      WRITE (6,*) 'F=',(LL,F(LL),LL=0,LM)
      WRITE (6,*) 'NQ (-1=T;0=Q;N=TRACER)=',NQ
      WRITE (6,*) 'RM IN SUBSID at L:',L,(RM(I,1,LL),LL=1,LM)
      STOP ' SUBSID END'
      END
