!@sum  STRAT_DIAG file for special E-P flux diagnostics from strat model
!@auth B. Suozzo/J/ Lerner
!@ver  1.0

      SUBROUTINE EPFLUX (U,V,T,P)
!@sum  EPFLUX calculates finite difference EP Flux on B-grid
!@auth B. Suozzo/J. Lerner
!@ver  1.0
C**** B-grid
C**** INPUT:   U - Zonal wind (corners) (m s-1)
C****          V - Meridional wind (corners) (m s-1)
C****          W - dp/dt (downward) (P-T grid level edges) (mb m2 s-1)
C****                (L=1 is ground level)
C****          T - Potential Temperature (K)  real T = TH*(P(mb))**K
C****          P - Pressure (mb) of column from 100mb to surface
C**** OUTPUT:  FMY,FEY - North. flux of mean,eddy ang. mom.   m3 s-2
C****          FMZ,FEZ - Vert. flux of mean,eddy ang. mom.    m3 mb s-2
C****         COR,CORR - Coriolis term and transf. coriolis term  m3 s-2
C****          FER1 - North. flux of error term 1   m2 s-2
C****          ER21,ER22 - error term 2 parts 1 & 2  m s-2
C****       FMYR,FEYR - Nor. flux of transf. mean,eddy ang. mom. m3 s-2
C****     FMZR,FEZR - Vert. flx of transf. mean,eddy ang. mom. m3 mb s-2
C****          RX - <v'th'>/<dth/dp> (U-V grid level edges)  (mb m s-1)
C****                (L=1 is ground level)
C****
C**** Note: W(,,1) is really PIT, the pressure tendency, but is not used
C****
      USE MODEL_COM, only : im,jm,lm,sig,dsig,sige,psfmpt,byim
      USE DOMAIN_DECOMP, only : GRID, GET, HALO_UPDATE, CHECKSUM, 
     *                          SOUTH, NORTH
      USE GEOM, only : dxv,rapvn,rapvs,fcor,dxyv,cosv,cosp
      USE  DAGCOM, only : ajl,kajl,kep,pl=>plm
      USE DYNAMICS, only : w=>conv     ! I think this is right....?

      IMPLICIT NONE
      REAL*8, INTENT(INOUT), 
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        U,V,T
      REAL*8, INTENT(IN), 
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: P

C**** NOTE: AEP was a separate array but is now saved in AJL (pointer?)
c      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,KEP) ::
c     *                                                            AEP

C**** ARRAYS CALCULATED HERE:
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        UV,UW
      REAL*8, POINTER, DIMENSION(:,:) :: 
C**** The following quantities are added into AEP(,,1..19)
     *        FMY,FEY,FMZ,FEZ,FMYR,FEYR,FMZR,FEZR,COR,CORR,
     *        FER1,ER21,ER22,VR,WR,RX,UI,VI,WI
C**** End of quantities added into AEP
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *        STB,TI,DUT,AX
      REAL*8 FD(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
      REAL*8  RXCL(LM),UCL(LM), WXXS(IM),WXXN(IM)
C**** The following quantities are added into AEP(,,1..19)
C**** (equivalenced to XEP)
cBMP  COMMON /EPCOM1/ FMY, FEY,  FMZ, FEZ, FMYR, FEYR, FMZR, FEZR, COR,
cBMP *     CORR, FER1, ER21, ER22, VR, WR, RX, UI, VI,  WI
      REAL*8, TARGET :: XEP(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,KEP)
cBMP EQUIVALENCE replaced by pointers      EQUIVALENCE (XEP,FMY)

      INTEGER I,J,L,N,IM1,IP1
      REAL*8 UCB,UCT,UCM,ALPH,RXC

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0H, J_1H, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      ALLOCATE ( FMY  (J_0H:J_1H,LM),
     &           FEY  (J_0H:J_1H,LM),
     &           FMZ  (J_0H:J_1H,LM),
     &           FEZ  (J_0H:J_1H,LM),
     &           FMYR (J_0H:J_1H,LM),
     &           FEYR (J_0H:J_1H,LM),
     &           FMZR (J_0H:J_1H,LM),
     &           FEZR (J_0H:J_1H,LM),
     &           COR  (J_0H:J_1H,LM),
     &           CORR (J_0H:J_1H,LM),
     &           FER1 (J_0H:J_1H,LM),
     &           ER21 (J_0H:J_1H,LM),
     &           ER22 (J_0H:J_1H,LM),
     &           VR   (J_0H:J_1H,LM),
     &           WR   (J_0H:J_1H,LM),
     &           RX   (J_0H:J_1H,LM),
     &           UI   (J_0H:J_1H,LM),
     &           VI   (J_0H:J_1H,LM),
     &           WI   (J_0H:J_1H,LM) )

      FMY  => XEP(:,:,1 )
      FEY  => XEP(:,:,2 )
      FMZ  => XEP(:,:,3 )
      FEZ  => XEP(:,:,4 )
      FMYR => XEP(:,:,5 )
      FEYR => XEP(:,:,6 )
      FMZR => XEP(:,:,7 )
      FEZR => XEP(:,:,8 )
      COR  => XEP(:,:,9 )
      CORR => XEP(:,:,10)
      FER1 => XEP(:,:,11)
      ER21 => XEP(:,:,12)
      ER22 => XEP(:,:,13)
      VR   => XEP(:,:,14)
      WR   => XEP(:,:,15)
      RX   => XEP(:,:,16)
      UI   => XEP(:,:,17)
      VI   => XEP(:,:,18)
      WI   => XEP(:,:,19)

C     Use IDACC(4) for calling frequency

C**** Initialise for this call
      XEP = 0

C****
C**** ZONAL AVERAGE QUANTITIES
C****
C**** UI(J,L) ... ZONAL AVERAGE U-WIND (m s-1)
C**** VI(J,L) ... ZONAL AVERAGE V-WIND (m s-1)
C**** TI(J,L) ... ZONAL AVERAGE POTENTIAL TEMPERATURE (K)
C**** WI(J,L) ... ZONAL AVERAGE VERTICAL WIND (mb m2 s-1)
      CALL AVGVI (U,UI)
      CALL AVGVI (V,VI)
      CALL AVGI (W,WI)
      CALL AVGI (T,TI)
C****
C**** STB(J,L) ... DELTA THETA             (PRESSURE EDGES)
C****
      DO L=2,LM
        DO J=J_0,J_1
          STB(J,L) = TI(J,L)-TI(J,L-1)
          IF(STB(J,L).LT.1d-2) THEN
CW         IF(STB(J,L).LT.1d-3)
CW   *     WRITE (*,*) 'STB < .001 AT J,L,STB:',J,L,STB(J,L)
            STB(J,L)=1d-2
          ENDIF
        END DO
      END DO
      STB(:,1)  = 1d-2

      Call CHECKSUM(grid, DXV  , __LINE__, __FILE__)
      Call CHECKSUM(grid, VI   , __LINE__, __FILE__)
      Call CHECKSUM(grid, UI   , __LINE__, __FILE__)
      Call CHECKSUM(grid, WI   , __LINE__, __FILE__)
      Call CHECKSUM(grid, U    , __LINE__, __FILE__)
      Call CHECKSUM(grid, V    , __LINE__, __FILE__)
      Call CHECKSUM(grid, RAPVN, __LINE__, __FILE__)
      Call CHECKSUM(grid, FCOR , __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, DXV  , from = NORTH)
      CALL HALO_UPDATE(grid, VI   , from = NORTH)
      CALL HALO_UPDATE(grid, UI   , from = NORTH)
      CALL HALO_UPDATE(grid, WI   , from = SOUTH)
      CALL HALO_UPDATE(grid, U    , from = NORTH)
      CALL HALO_UPDATE(grid, V    , from = NORTH)
      CALL HALO_UPDATE(grid, RAPVN, from = SOUTH)
      CALL HALO_UPDATE(grid, FCOR , from = SOUTH)

C****
C**** CORRELATIONS OF VARIOUS QUANTITIES
C****
C****
C**** FMY .............. (VUdx) NORTHWARD MEAN MOMENTUM TRANSPORT
C****                         (m3 s-1)
      DO L=1,LM
        IF (HAVE_SOUTH_POLE) FMY(1,L)=0.
        IF (HAVE_NORTH_POLE) FMY(JM,L)=0.
        DO J=J_0S,J_1S
          FMY(J,L)=.25*(DXV(J)*VI(J,L)+DXV(J+1)*VI(J+1,L))
     *         * (UI(J,L)+UI(J+1,L))
        END DO
C****
C**** FEY, UV(I,J,L) ... (V'U'dx) NORTHWARD EDDY MOMENTUM TRANSPORT
C****                           (m3 s-2)
        DO I=1,IM
          IF (HAVE_SOUTH_POLE) UV(I,1,L)=0.
          IF (HAVE_NORTH_POLE) UV(I,JM,L)=0.
        END DO
        DO J=J_0S,J_1S
          IM1=IM-1
          I=IM
          DO IP1=1,IM
            UV(I,J,L) = (
     *           DXV(J)*(V(IM1,J,L)+2*V(I,J,L)+V(IP1,J,L)-4*VI(J,L)) +
     *           DXV(J+1)*(V(IM1,J+1,L)+2*V(I,J+1,L)+V(IP1,J+1,L)-4*VI(J
     *           +1,L)))* (U(I,J,L)+U(I,J+1,L)-(UI(J,L)+UI(J+1,L)))
            IM1=I
            I=IP1
          END DO
        END DO
      END DO
      CALL AVGI (UV,FEY)
      DO L=1,LM
      DO J=J_0,J_1
        FEY(J,L)=.0625d0*FEY(J,L)
      END DO
      END DO
C****
C**** FMZ .............. (WUdA) VERTICAL MEAN MOMENTUM TRANSPORT
C**** FEZ, UW(I,J,L) ... (U'W'dA) VERTICAL EDDY MOMENTUM TRANSPORT
C****                           (m3 mb s-2)
      UW(:,:,1)=0.
      DO L=2,LM
        WXXS(1:IM)=0.
        UW(1:IM,1,L)=0.
        DO J=J_0STG,J_1STG
          FMZ(J,L) = (WI(J-1,L)*RAPVN(J-1)+WI(J,L)*RAPVS(J))
     *         * (UI(J,L-1)+UI(J,L))
          I=IM
          DO IP1=1,IM
            WXXN(I)=(W(I,J,L)+W(IP1,J,L))-2*WI(J,L)
            UW(I,J,L) = (WXXS(I)*RAPVN(J-1)+WXXN(I)*RAPVS(J))
     *           * ((U(I,J,L-1)+U(I,J,L))-(UI(J,L-1)+UI(J,L)))
            WXXS(I)=WXXN(I)
            I=IP1
          END DO
        END DO
      END DO
      CALL AVGI (UW,FEZ)
      FMZ(:,1)=0.
      FEZ(:,1)=0.
      FMZ(1,2:LM)=0.
      FEZ(1,2:LM)=0.
      DO L=1,LM
      DO J=J_0,J_1
        FEZ(J,L)=.5*FEZ(J,L)
      END DO
      END DO
C****
C**** COR(J,L) .......... CORIOLIS FORCE
C****                        (m3 s-2)
      DO L=1,LM
        IF (HAVE_SOUTH_POLE) 
     *      COR( 2,L)=.5*(2.*FCOR( 1)+FCOR(   2))*VI(2,L)
        IF (HAVE_NORTH_POLE) 
     *      COR(JM,L)=.5*(2.*FCOR(JM)+FCOR(JM-1))*VI(JM,L)
        DO J=3,JM-1
          COR(J,L)=.5*(FCOR(J-1)+FCOR(J))*VI(J,L)
        END DO
      END DO

      Call CHECKSUM(grid, UV  , __LINE__, __FILE__)
      Call CHECKSUM(grid, UW  , __LINE__, __FILE__)
      Call CHECKSUM(grid, COR , __LINE__, __FILE__)

cBMP FD calc moved out of next loop for ghosting purposes
      DO L=1,LM
        IM1=IM
        DO I=1,IM
          IF (HAVE_SOUTH_POLE) FD(I,1,L)= 0.
          IF (HAVE_NORTH_POLE) FD(I,JM,L)=0.
          DO J=J_0S,J_1S
            FD(I,J,L) = (DXV(J)-DXV(J+1)) *
     *            .25*(U(IM1,J,L)+U(I,J,L)+U(IM1,J+1,L)+U(I,J+1,L))
          END DO
          IM1=I
        END DO
      END DO
      Call CHECKSUM(grid, FD , __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, FD , from = SOUTH)
cBMP

C****
C**** ERROR TERMS
C****
      DO L=1,LM
        FER1(:,L)=0.
        ER21(:,L)=0.
        ER22(:,L)=0.
C**** ERROR TERM 1 IS DIVERGENCE OF FER1 (m2 s-2)
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          IF (HAVE_SOUTH_POLE) 
     *        FER1(1,L) =FER1(1,L) +(U(I,2,L)**2   - U(I,2,L)**2)
          IF (HAVE_NORTH_POLE)
     *        FER1(JM,L)=FER1(JM,L)+(U(IM1,JM,L)**2- U(IP1,JM,L)**2)
          DO J=J_0S,J_1S
C? FER1 does not include COSP, yet ER1 = 1/COSV*(FER1(J)-FER1(J-1))
            FER1(J,L)= FER1(J,L) + ((U(IM1,J,L)+U(I,J+1,L))**2
     *           - (U(IP1,J,L)+U(I,J+1,L))**2)
          END DO
          IM1=I
          I=IP1
        END DO
C**** ERROR TERM 2 - THE METRIC TERM AS IN RUN 999
        IM1=IM
        DO I=1,IM
c         IF (HAVE_SOUTH_POLE) FD(I,1,L)= 0.
c         IF (HAVE_NORTH_POLE) FD(I,JM,L)=0.
c         DO J=J_0S,J_1S
c           FD(I,J,L) = (DXV(J)-DXV(J+1)) * 
c    *           .25*(U(IM1,J,L)+U(I,J,L)+U(IM1,J+1,L)+U(I,J+1,L))
c         END DO
          DO J=J_0STG,J_1STG
            ALPH=.25*(FD(I,J-1,L)+FD(I,J,L))
            ER22(J,L)=ER22(J,L)+ALPH*(V(IM1,J,L)+V(I,J,L))
          END DO
          IM1=I
        END DO
CE*** ERROR TERM 2, PART 2  -  ER22  (m s-2)
CE    IM1=IM-1
CE    I=IM
CE    DO IP1=1,IM
CE    IF (HAVE_SOUTH_POLE) UXXNS=(U(IM1,2,L)+2*U(I,2,L)+U(IP1,2,L))
CE    IF (HAVE_SOUTH_POLE) XUXXYS=2*UXXNS*DXP(2)
CE    DO J=J_0S,J_1S
CE    UXXNN=(U(IM1,J+1,L)+2*U(I,J+1,L)+U(IP1,J+1,L))
CE    XUXXYN=(UXXNS+UXXNN)*(DXP(J+1)-DXP(J-1))
CE    ER22(J,L) = ER22(J,L)-V(I,J,L)*(XUXXYS+XUXXYN)
CE    XUXXYS=XUXXYN
CE    UXXNS=UXXNN
CE    END DO
CE    IF (HAVE_NORTH_POLE) XUXXYN=-2*UXXNS*DXP(JM-1)
CE    IF (HAVE_NORTH_POLE) ER22(JM,L) = ER22(JM,L)-V(I,JM,L)*(XUXXYS+XUXXYN)
CE    IM1=I
CE    I=IP1
CE    END DO
C**** NORMALIZE ERROR TERMS AND CONVERT TO OUTPUT UNITS
      DO J=J_0,J_1
        FER1(J,L)=1./(48*IM)*FER1(J,L)
      END DO
CE      DO J=J_0STG,J_1STG
CE        ER21(J,L) = -1./(2*DXYV(J)*COSV(J))*
CE   *  (FMY(J-1,L)+FEY(J-1,L)+FMY(J,L)+FEY(J,L))*(COSP(J)-COSP(J-1))
CE        ER22(J,L) = 1./(32*IM*DXYV(J))*ER22(J,L)
CE      END DO
      DO J=J_0STG,J_1STG
        ER22(J,L) = 1./(IM*DXYV(J))*ER22(J,L)
      END DO
      END DO
C****
C****  TRANSFORMED CIRCULATION
C****

      Call CHECKSUM(grid, FER1  , __LINE__, __FILE__)
      Call CHECKSUM(grid, ER22  , __LINE__, __FILE__)

C****
C**** RX(J,L) ..... ([V'TH']/[DTH/DP]) TRANSFORMATION GAUGE
C****               (U-wind grid, level edges)  (m mb s-1)

      CALL CHECKSUM(grid, T   , __LINE__, __FILE__)
      CALL CHECKSUM(grid, TI  , __LINE__, __FILE__)
      CALL CHECKSUM(grid, STB , __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, T   , from = SOUTH)
      CALL HALO_UPDATE(grid, TI  , from = SOUTH)
      CALL HALO_UPDATE(grid, STB , from = SOUTH)

      DO L=2,LM
      DO J=J_0STG,J_1STG
        RX(J,L) = 0.
        I=IM
        DO IP1=1,IM
          RX(J,L) = RX(J,L) + (V(I,J,L)+V(IP1,J,L)-2*VI(J,L))
     *         * (T(I,J-1,L)+T(I,J,L) - (TI(J-1,L)+TI(J,L)))
     *         +              (V(I,J,L-1)+V(IP1,J,L-1)-2*VI(J,L-1))
     *         * (T(I,J-1,L-1)+T(I,J,L-1) - (TI(J-1,L-1)+TI(J,L-1)))
          I=IP1
        END DO
        RX(J,L) = .25*BYIM * RX(J,L) *
     *       (PL(L)-PL(L-1))/(STB(J-1,L)+STB(J,L))
      END DO
      END DO
      RX(2:JM,1)=0.

      CALL CHECKSUM(grid, RX, __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, RX, from = NORTH)

C****
C**** VR(J,L) ........  TRANSFORMED WIND
C**** WR(J,L) ........
C****
      DO L=1,LM-1
      DO J=J_0STG,J_1STG
        VR(J,L)=VI(J,L)+(RX(J,L+1)-RX(J,L))/(PSFMPT*DSIG(L))
      END DO
      END DO
      DO J=J_0STG,J_1STG
        VR(J,LM)=VI(J,LM)-RX(J,LM)/(PSFMPT*DSIG(LM))
      END DO
      DO L=1,LM
      DO J=J_0S,J_1S
        WR(J,L)=WI(J,L)+(RX(J+1,L)*DXV(J+1)-RX(J,L)*DXV(J))
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) THEN
         DO L=1,LM
            WR(1,L)=WI(1,L) + RX(2,L)*DXV(2)
         END DO
      ENDIF
      IF (HAVE_NORTH_POLE) THEN
         DO L=1,LM
            WR(JM,L)=WI(JM,L) - RX(JM,L)*DXV(JM)
         END DO
      ENDIF

      CALL CHECKSUM(grid, VR  , __LINE__, __FILE__)
      CALL CHECKSUM(grid, WR  , __LINE__, __FILE__)
      CALL CHECKSUM(grid, COSV, __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, VR  , from = NORTH)
      CALL HALO_UPDATE(grid, WR  , from = SOUTH)
      CALL HALO_UPDATE(grid, COSV, from = NORTH)

C****
C**** TRANSFORMED MEAN FLUXES
C****
C**** FMYR(J,L)     (m3 s-2)
C**** FMZR(J,L)     (m3 mb s-2)
C**** CORR(J,L)     (m3 s-2)
C****
      DO L=1,LM
      DO J=J_0S,J_1S
        FMYR(J,L)=.25*(VR(J,L)*DXV(J)+VR(J+1,L)*DXV(J+1))
     *       *  (UI(J,L)+UI(J+1,L))
      END DO
      END DO
      DO L=2,LM
      DO J=J_0STG,J_1STG
        FMZR(J,L) = (WR(J-1,L)*RAPVN(J-1)+WR(J,L)*RAPVS(J))
     *       *  (UI(J,L-1)+UI(J,L))
      END DO
      END DO
      DO L=1,LM
        IF (HAVE_SOUTH_POLE) 
     *      CORR( 2,L)=.5*(2.*FCOR( 1)+FCOR(   2))*VR(2,L)
        IF (HAVE_NORTH_POLE) 
     *      CORR(JM,L)=.5*(2.*FCOR(JM)+FCOR(JM-1))*VR(JM,L)
        DO J=3,JM-1
          CORR(J,L)=.5*(FCOR(J-1)+FCOR(J))*VR(J,L)
        END DO
      END DO
C****
C**** TRANSFORMED EDDY FLUXES
C****
C**** FEYR(J,L)
C****
      DO L=1,LM
        RXCL(L)=RX(2,L)*COSV(2)
        UCL(L)=UI(2,L)*DXV(2)
      END DO
      DO J=J_0S,J_1S
        UCB=(UI(J+1,1)*DXV(J+1)+UCL(1))
        UCM=(UI(J+1,2)*DXV(J+1)+UCL(2))
        FEYR(J,1)=0.
        FEYR(J,LM)=0.
        DO L=2,LM-1
          UCT=(UI(J+1,L+1)*DXV(J+1)+UCL(L+1))
          FEYR(J,L) = .125d0/(COSP(J)*PSFMPT*DSIG(L))
     *         *    ((RX(J+1,L)*COSV(J+1)+RXCL(L))*(UCM-UCB)
     *         +    (RX(J+1,L+1)*COSV(J+1)+RXCL(L+1))*(UCT-UCM))
          UCB=UCM
          UCM=UCT
        END DO
        DO L=1,LM
          RXCL(L)=RX(J+1,L)*COSV(J+1)
          UCL(L)=UI(J+1,L)*DXV(J+1)
        END DO
      END DO
C****
C**** FEZR
C****
      RXCL(:)=0.
      UCL(:)=0.
      DO J=J_0S,J_1S
      DO L=2,LM
        RXC = RX(J,L)*COSV(J)
        UCB = (UI(J,L)+UI(J,L-1))*DXV(J)
        FEZR(J,L) = -.5*(FCOR(J-1)+FCOR(J))*RX(J,L)
     *       + .125d0/COSV(J)* ((RXC+RXCL(L))
     *       *     ((UI(J,L)+UI(J,L-1))*DXV(J) - (UCL(L)+UCL(L-1)))
     *       +     (RX(J+1,L)*COSV(J+1)+RXC)
     *       *     ((UI(J+1,L)+UI(J+1,L-1))*DXV(J+1) - UCB))
      END DO
      DO L=1,LM
        RXCL(L)=RX(J,L)*COSV(J)
        UCL(L)=UI(J,L)*DXV(J)
      END DO
      END DO
      IF (HAVE_NORTH_POLE) THEN
         DO L=2,LM
           RXC=RX(JM,L)*COSV(JM)
           UCB=(UI(JM,L)+UI(JM,L-1))*DXV(JM)
           FEZR(JM,L) = -FCOR(JM)*RX(JM,L)
     *          + .125d0/COSV(JM)* (RXC+RXCL(L))
     *          *     ((UI(JM,L)+UI(JM,L-1))*DXV(JM) 
     *          -      (UCL(L)+UCL(L-1)))
          END DO
      ENDIF
C**** Add Eulerian circulation to transformed eddy components
      DO L=1,LM
      DO J=J_0,J_1
        FEYR(J,L)=FEY(J,L)+FEYR(J,L)
        FEZR(J,L)=FEZ(J,L)+FEZR(J,L)
      END DO
      END DO

      CALL CHECKSUM(grid, FMYR, __LINE__, __FILE__)
      CALL CHECKSUM(grid, FMYR, __LINE__, __FILE__)
      CALL CHECKSUM(grid, FEYR, __LINE__, __FILE__)
      CALL CHECKSUM(grid, FEZR, __LINE__, __FILE__)

C****
C**** ACCUMULATE EP FLUXES
C****
      DO N=1,KEP-2
      DO L=1,LM
      DO J=J_0,J_1
        AJL(J,L,KAJL-KEP+N)=AJL(J,L,KAJL-KEP+N)+XEP(J,L,N)
c        AEP(J,L,N)=AEP(J,L,N)+XEP(J,L,N)
      END DO
      END DO
      END DO
      RETURN
C****
      ENTRY EPFLXI (U)
      CALL AVGVI (U,AJL(1,1,KAJL))
c      CALL AVGVI (U,AEP(1,1,KEP))
CW          WRITE (36,'('' TAU='',F12.0)') TAU
CW          CALL WRITJL ('U - INITIAL     ',AEP(1,1,KEP),1.)
      RETURN
C****
      ENTRY EPFLXF (U)
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)
      CALL AVGVI (U,AX)
      DO L=1,LM
      DO J=J_0STG,J_1STG
        AJL(J,L,KAJL-1) = AX(J,L)-AJL(J,L,KAJL)
c       AEP(J,L,KEP-1) = AX(J,L)-AEP(J,L,KEP)
      END DO
      END DO
CW          CALL WRITJL ('U - FINAL       ',AX,1.)
CW          CALL WRITJL ('DU - TOTAL      ',AEP(1,1,KEP-1,1.)
CW         WRITE (36,'('' TAU='',F12.0,''  IDUM(1)='',I6)') TAU,IDUM(1)
      RETURN
      END

      SUBROUTINE EPFLXP
!@sum  EPFLXP prints out diagnostics of E-P Fluxes
!@auth B. Suozzo/J. Lerner
!@ver  1.0
C****
C**** B-grid
C**** INPUT: AEP contains the following quantities summed over time:
C****        FMY,FEY - North. flux of mean,eddy ang. mom.   m3 s-2
C****        FMZ,FEZ - Vert. flux of mean,eddy ang. mom.    m3 mb s-2
C****        COR,CORR - coriolis term and transf. coriolis term  m3 s-2
C****        FER1 - North. flux of error term 1   m2 s-2
C****        ER21,ER22 - error term 2 parts 1 & 2  m s-2
C****     FMYR,FEYR - Nor. flux of transf. mean,eddy ang. mom. m3 s-2
C****   FMZR,FEZR - Vert. flx of transf. mean,eddy ang. mom. m3 mb s-2
C**** OUTPUT:
C****        DMF,DEF - Divergence of mean,eddy ang. mom.  m3 s-2
C****        DMFR,DEFR - Div. of transf. mean,eddy ang. mom.  m3 s-2
C****        COR,CORR - same as above   m3 s-2
C****        ER1,ER2 - error terms 1 and 2   m s-2
C****   DUD,DUR - Delta U by Eulerian and transf. circulation  m s-2
C****
      USE MODEL_COM, only : im,jm,lm,sig,dsig,dtsrce=>dtsrc,psfmpt,fim
     *     ,idacc,ndaa,ls1,ptop
      USE DOMAIN_DECOMP, only : GRID, GET, HALO_UPDATE, CHECKSUM,
     *                          SOUTH, NORTH
      USE GEOM, only : dxyv,bydxyv,cosv,cosp,dxv,dyv
      USE DAGCOM, only : ajl,kajl,kep,apj
     &     ,jl_dudfmdrg,jl_dumtndrg,jl_dushrdrg
     &     ,jl_dumcdrgm10,jl_dumcdrgp10
     &     ,jl_dumcdrgm40,jl_dumcdrgp40
     &     ,jl_dumcdrgm20,jl_dumcdrgp20
     &     ,jl_dudtsdif,jl_damdc,jl_dammc,jl_dudtvdif
      IMPLICIT NONE
C**** NOTE: AEP was a separate array but is now saved in AJL (pointer?)
c      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,KEP) ::
c     *                                                             AEP

c      COMMON /PROGCB/ U,V,T,SX,SY,SZ,P,Q   !not used?
C**** diagnostic information for print out
C**** this should be in an init_ep routine or something
      integer, parameter :: njl_out=7
      character(len=50) :: lname(njl_out) = (/
     *     'DU/DT BY EULER CIRC. + CONVEC + DRAG+DIF+ER2',
     *     'DU/DT BY MEAN ADVECTION                     ',
     *     'DU/DT BY EDDY CONVERGENCE                   ',
     *     'DU/DT BY TRANSFORMED ADVECTION              ',
     *     'DU/DT BY ELIASSEN-PALM DIVERGENCE           ',
     *     'DU/DT BY F.D. ERROR TERM 1                  ',
     *     'DU/DT BY F.D. ERROR TERM 2                  '/)
      character(len=30) :: sname(njl_out) = (/
     *     'dudt_sum1    ','dudt_meanadv ','dudt_eddycnv '
     *     ,'dudt_trnsadv ','dudt_epflxdiv','dudt_fderr1  '
     *     ,'dudt_fderr2  '/)
      character(len=50) :: units(njl_out) = 'm/s^2'
      integer, dimension(njl_out) :: pow = (/ -6,-6,-6,-6,-6,-6,-6 /)

C**** ARRAYS CALCULATED HERE:
      REAL*8 ONES(JM+LM),PMO(LM),DP(LM)
      REAL*8 DXCOSV(GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      REAL*8, POINTER, DIMENSION(:,:) :: 
     *   FMY,FEY,FMZ,FEZ,FMYR,FEYR,FMZR,FEZR,COR,CORR,FER1,ER21
     *  ,ER22,VR,WR,RX,UI,VI,WI,DUT,DUD
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *   DUDS,DUR,DMF,DEF,DMFR,DEFR
CW   *  ,DMY,DMZ,DEY,DEZ,DMYR,DMZR,DEYR,DEZR
     *  ,ER1,ER2,BYDPJL
C**** common needed to pass from XEP to individual arrays (better way??)
C**** Not clear how many of these are actually used.
cBMP  COMMON /EPCOM/ FMY, FEY, FMZ, FEZ, FMYR, FEYR, FMZR, FEZR, COR,
cBMP *     CORR, FER1, ER21, ER22, VR, WR, RX, UI, VI, WI, DUT, DUD
      REAL*8, TARGET, 
     *        DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,KEP) :: XEP
cBMP  EQUIVALENCE (XEP,FMY)
      REAL*8 DTAEP,BYDT,SCALEP,SCALE1,BYIAEP
      INTEGER I,J,L,N,JL

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0H, J_1H, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      ALLOCATE ( FMY  (J_0H:J_1H,LM),
     &           FEY  (J_0H:J_1H,LM),
     &           FMZ  (J_0H:J_1H,LM),
     &           FEZ  (J_0H:J_1H,LM),
     &           FMYR (J_0H:J_1H,LM),
     &           FEYR (J_0H:J_1H,LM),
     &           FMZR (J_0H:J_1H,LM),
     &           FEZR (J_0H:J_1H,LM),
     &           COR  (J_0H:J_1H,LM),
     &           CORR (J_0H:J_1H,LM),
     &           FER1 (J_0H:J_1H,LM),
     &           ER21 (J_0H:J_1H,LM),
     &           ER22 (J_0H:J_1H,LM),
     &           VR   (J_0H:J_1H,LM),
     &           WR   (J_0H:J_1H,LM),
     &           RX   (J_0H:J_1H,LM),
     &           UI   (J_0H:J_1H,LM),
     &           VI   (J_0H:J_1H,LM),
     &           WI   (J_0H:J_1H,LM),
     &           DUT  (J_0H:J_1H,LM),
     &           DUD  (J_0H:J_1H,LM) )

      FMY  => XEP(:,:,1 )
      FEY  => XEP(:,:,2 )
      FMZ  => XEP(:,:,3 )
      FEZ  => XEP(:,:,4 )
      FMYR => XEP(:,:,5 )
      FEYR => XEP(:,:,6 )
      FMZR => XEP(:,:,7 )
      FEZR => XEP(:,:,8 )
      COR  => XEP(:,:,9 )
      CORR => XEP(:,:,10)
      FER1 => XEP(:,:,11)
      ER21 => XEP(:,:,12)
      ER22 => XEP(:,:,13)
      VR   => XEP(:,:,14)
      WR   => XEP(:,:,15)
      RX   => XEP(:,:,16)
      UI   => XEP(:,:,17)
      VI   => XEP(:,:,18)
      WI   => XEP(:,:,19)
      DUT  => XEP(:,:,20)
      DUD  => XEP(:,:,21)

C**** Initialize constants
      DTAEP = DTsrce*NDAA   ! change of definition of NDAA
      BYIAEP=1./(IDACC(4)+1.D-20)
      BYDT=1./(DTSRCE*IDACC(1)+1.D-20)
      DO J=J_0STG,J_1STG
        DXCOSV(J) = DXV(J)*COSV(J)
      END DO

c GISS-ESMF EXCEPTIONAL CASE
      DO JL=1,JM+LM
        ONES(JL)=1.
      END DO
      DO L=1,LM
        DP(L)  = PSFMPT*DSIG(L)
        PMO(L) = PSFMPT*SIG(L)+PTOP
      END DO
      DO L=1,LS1-1
      DO J=J_0STG,J_1STG
        BYDPJL(J,L)=(FIM*IDACC(4))/(0.25*DSIG(L)*APJ(J,2)+1.D-20)
      END DO
      END DO
      DO L=LS1,LM
      DO J=J_0STG,J_1STG
        BYDPJL(J,L)=1./DP(L)
      END DO
      END DO
C**** Normalize AEP
C     CALL EPFLXF (U)
      DO N=1,KEP-2
      DO L=1,LM
      DO J=J_0,J_1
        XEP(J,L,N)=AJL(J,L,KAJL-KEP+N)*BYIAEP
c        XEP(J,L,N)=AEP(J,L,N)*BYIAEP
      END DO
      END DO
      END DO
      DO N=KEP-1,KEP
      DO L=1,LM
      DO J=J_0,J_1
        XEP(J,L,N)=AJL(J,L,KAJL-KEP+N)
c        XEP(J,L,N)=AEP(J,L,N)
      END DO
      END DO
      END DO
C****
C**** DMF,DEF by horizontal convergence
C****

      Call CHECKSUM(grid, FMY  , __LINE__, __FILE__)
      Call CHECKSUM(grid, FEY  , __LINE__, __FILE__)
      Call CHECKSUM(grid, FMYR , __LINE__, __FILE__)
      Call CHECKSUM(grid, FEYR , __LINE__, __FILE__)
      Call CHECKSUM(grid, COSP , __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, FMY  , from = SOUTH)
      CALL HALO_UPDATE(grid, FEY  , from = SOUTH)
      CALL HALO_UPDATE(grid, FMYR , from = SOUTH)
      CALL HALO_UPDATE(grid, FEYR , from = SOUTH)
      CALL HALO_UPDATE(grid, COSP , from = SOUTH)

      DO L=1,LM
      DO J=J_0STG,J_1STG
        DMF(J,L)  = (FMY(J-1,L)*COSP(J-1)-FMY(J,L)*COSP(J))/COSV(J)
        DEF(J,L)  = (FEY(J-1,L)*COSP(J-1)-FEY(J,L)*COSP(J))/COSV(J)
C       DMY(J,L)  = (FMY(J-1,L)*COSP(J-1)-FMY(J,L)*COSP(J))/COSV(J)
C       DEY(J,L)  = (FEY(J-1,L)*COSP(J-1)-FEY(J,L)*COSP(J))/COSV(J)
C       DMYR(J,L) = (FMYR(J-1,L)*COSP(J-1)-FMYR(J,L)*COSP(J))/COSV(J)
C       DEYR(J,L) = (FEYR(J-1,L)*COSP(J-1)-FEYR(J,L)*COSP(J))/COSV(J)
        DMFR(J,L) = (FMYR(J-1,L)*COSP(J-1)-FMYR(J,L)*COSP(J))/COSV(J)
        DEFR(J,L) = (FEYR(J-1,L)*COSP(J-1)-FEYR(J,L)*COSP(J))/COSV(J)
      END DO
      END DO
C****
C**** Add DMF,DEF by vertical convergence
C****
      DO J=J_0STG,J_1STG
        DMF(J,LM)  = DMF(J,LM)  - FMZ(J,LM)/BYDPJL(J,LM)
        DEF(J,LM)  = DEF(J,LM)  - FEZ(J,LM)/BYDPJL(J,LM)
CW      DMZ(J,LM)  = - FMZ(J,LM)/BYDPJL(J,LM)
CW      DEZ(J,LM)  = - FEZ(J,LM)/BYDPJL(J,LM)
CW      DMZR(J,LM) = - FMZR(J,LM)/BYDPJL(J,LM)
CW      DEZR(J,LM) = - FEZR(J,LM)/BYDPJL(J,LM)
        DMFR(J,LM) = DMFR(J,LM) - FMZR(J,LM)/BYDPJL(J,LM)
        DEFR(J,LM) = DEFR(J,LM) - FEZR(J,LM)/BYDPJL(J,LM)
        DO L=1,LM-1
          DMF(J,L)  = DMF(J,L)  + (FMZ(J,L+1)-FMZ(J,L))*BYDPJL(J,L)
          DEF(J,L)  = DEF(J,L)  + (FEZ(J,L+1)-FEZ(J,L))*BYDPJL(J,L)
CW        DMZ(J,L)  = (FMZ(J,L+1)-FMZ(J,L))*BYDPJL(J,L)
CW        DEZ(J,L)  = (FEZ(J,L+1)-FEZ(J,L))*BYDPJL(J,L)
CW        DMZR(J,L) = (FMZR(J,L+1)-FMZR(J,L))*BYDPJL(J,L)
CW        DEZR(J,L) = (FEZR(J,L+1)-FEZR(J,L))*BYDPJL(J,L)
          DMFR(J,L) = DMFR(J,L) + (FMZR(J,L+1)-FMZR(J,L))*BYDPJL(J,L)
          DEFR(J,L) = DEFR(J,L) + (FEZR(J,L+1)-FEZR(J,L))*BYDPJL(J,L)
        END DO
      END DO
C****
C**** ADD COR,CORR
C****
      DO L=1,LM
        DO J=J_0STG,J_1STG
          DMF(J,L)  = DMF(J,L)  + COR(J,L)
C         DMY(J,L)  = DMY(J,L)  + COR(J,L)
C         DMYR(J,L) = DMYR(J,L) + CORR(J,L)
          DMFR(J,L) = DMFR(J,L) + CORR(J,L)
        END DO
      END DO
C****
C**** ER1 by horizontal convergence (m s-2)
C****

      Call CHECKSUM(grid, FER1 , __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, FER1 , from = SOUTH)

      DO L=1,LM
      DO J=J_0STG,J_1STG
        ER1(J,L) = 1./(DYV(J)*COSV(J))*(FER1(J-1,L)-FER1(J,L))
        ER2(J,L) = ER21(J,L)+ER22(J,L)
      END DO
      END DO
C****
C**** DUD,DUR total change in Eulerian, Transformed wind  (m s-2)
C****
      DO L=1,LM
      DO J=J_0STG,J_1STG
        DUR(J,L)  = BYDXYV(J)*(DMFR(J,L) +DEFR(J,L))
        DUD(J,L)  = BYDXYV(J)*(DMF(J,L)+DEF(J,L))
      END DO
      END DO
C****
C**** Let Transformed == Transformed - Eulerian
C****
CW      DO L=1,LM
CW      DO J=J_0,J_1
CW        DUR(J,L) = DUR(J,L) - DUD(J,L)
CW        DMFR(J,L) = DMFR(J,L) - DMF(J,L)
CW        DEFR(J,L) = DEFR(J,L) - DEF(J,L)
CW      END DO
CW      END DO
C****
C**** Print maps of EP fluxes
C**** note: JLMAP (lname,sname,units,power,Pres,Array,ScalP,ScalJ,ScalL)
C****          prints maps of  (Array * SCALEP * ScaleJ * ScaleL)
C****
      DO L=1,LM
      DO J=J_0STG,J_1STG
        DUDS(J,L)=((AJL(J,L,JL_DUDFMDRG)+AJL(J,L,JL_DUMTNDRG))+
     &             (AJL(J,L,JL_DUSHRDRG)+AJL(J,L,JL_DUMCDRGM10))+
     *       ((AJL(J,L,JL_DUMCDRGP10)+AJL(J,L,JL_DUMCDRGM40))+
     &        (AJL(J,L,JL_DUMCDRGP40)+AJL(J,L,JL_DUMCDRGM20)))+
     *        (AJL(J,L,JL_DUMCDRGP20)+
     &         AJL(J,L,JL_DUDTSDIF)+AJL(J,L,JL_DUDTVDIF)))
      END DO
      END DO
      SCALE1=1./(FIM*DTSRCE*IDACC(1)+1.D-20)
      DO L=1,LM
      DO J=J_0STG,J_1STG
        DUDS(J,L) = DUDS(J,L)+
     &        (AJL(J,L,JL_DAMDC)+AJL(J,L,JL_DAMMC))*BYDPJL(J,L)
      END DO
      END DO
      DO L=1,LM
      DO J=J_0STG,J_1STG
        DUDS(J,L) = (DUD(J,L)-ER2(J,L)) + DUDS(J,L)*SCALE1
      END DO
      END DO
      SCALEP=1
CW      CALL WRITJL ('DUDT: EUL+SOURCE',DUDS,SCALEP)
CW      CALL WRITJL ('DUDT: ENTIRE GCM',DUT,SCALEP) ! AJK-47 DIAGJK
      CALL JLMAP (LNAME(1),SNAME(1),UNITS(1),POW(1),PMO,DUDS,SCALEP,ONES
     *     ,ONES,LM,2,2)
C****
      CALL JLMAP (LNAME(2),SNAME(2),UNITS(2),POW(2),PMO,DMF,SCALEP
     *     ,BYDXYV,ONES,LM,2,2)
      CALL JLMAP (LNAME(3),SNAME(3),UNITS(3),POW(3),PMO,DEF,SCALEP
     *     ,BYDXYV,ONES,LM,2,2)
      CALL JLMAP (LNAME(4),SNAME(4),UNITS(4),POW(4),PMO,DMFR,SCALEP
     *     ,BYDXYV,ONES,LM,2,2)
      CALL JLMAP (LNAME(5),SNAME(5),UNITS(5),POW(5),PMO,DEFR,SCALEP
     *     ,BYDXYV,ONES,LM,2,2)
      CALL JLMAP (LNAME(6),SNAME(6),UNITS(6),POW(6),PMO,ER1,SCALEP,ONES
     *     ,ONES,LM,2,2)
      CALL JLMAP (LNAME(7),SNAME(7),UNITS(7),POW(7),PMO,ER2,SCALEP,ONES
     *     ,ONES,LM,2,2)
C****
CW      DO L=1,LM
CW      DO J=J_0STG,J_1STG
CW      DMF(J,L)=DMF(J,L)*BYDXYV(J)
CW      DEF(J,L)=DEF(J,L)*BYDXYV(J)
CW      DMFR(J,L)=DMFR(J,L)*BYDXYV(J)
CW      DEFR(J,L)=DEFR(J,L)*BYDXYV(J)
CW      END DO
CW      END DO
CW      WRITE (36,'(''DU/DT = 10**-6 M S-2'')')
CW      WRITE (36,'(''TR: TRANSFORMED - EULERIAN'')')
CW      SCALEP = 1.E6
CW      CALL WRITJL ('DUDT: MEAN EULER',DMF,SCALEP)
CW      CALL WRITJL ('DUDT: EDDY EULER',DEF,SCALEP)
CW      CALL WRITJL ('DUDT: EULER CIRC',DUD,SCALEP)
CW      CALL WRITJL ('DUDT: MEAN TRANS',DMFR,SCALEP)
CW      CALL WRITJL ('DUDT: EDDY TRANS',DEFR,SCALEP)
CW      CALL WRITJL ('DUDT: TRANS-EULE',DUR,SCALEP)
      RETURN
      END

      SUBROUTINE AVGI (X,XI)
!@sum  AVGI average a 3-dimensional array in the x-direction
!@auth B. Suozzo
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,BYIM
      USE DOMAIN_DECOMP, only : GRID, GET
      USE GEOM, only : imaxj
      IMPLICIT NONE

!@var X input 3-D array
      REAL*8, INTENT(IN), 
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: X
!@var XI output zonally averaged 2-D array
      REAL*8, INTENT(OUT), 
     &        DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: XI
      INTEGER I,J,L
      REAL*8 XXI

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      DO L=1,LM
        DO J=J_0,J_1
          XXI=0.
          DO I=1,IMAXJ(J)
            XXI = XXI + X(I,J,L)
          END DO
          XI(J,L) = XXI/IMAXJ(J)
        END DO
      END DO
      RETURN
C****
      ENTRY AVGVI (X,XI)
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)
!@sum  AVGVI average a 3-dimensional array in the x-direction (no pole)
!@auth B. Suozzo
!@ver  1.0
C****
      DO L=1,LM
        DO J=J_0STG,J_1STG
          XXI=0.
          DO I=1,IM
            XXI = XXI + X(I,J,L)
          END DO
          XI(J,L) = XXI*BYIM
        END DO
      END DO
C****
      RETURN
      END

