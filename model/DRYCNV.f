      SUBROUTINE DIFFUS(LBASE_MIN,LBASE_MAX,DTIME)
!@sum  DIFFUS(DRYCNV) mixes air caused by dry convection.
!@+    this version checks base layers lbase_min to lbase_max.
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : lhe,sha,deltx
      USE MODEL_COM
      USE GEOM
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
      USE DAGCOM, only : ajl,jl_trbhr,jl_damdc,jl_trbdlht
      USE DYNAMICS, only : pk,pdsig,plij
      USE PBLCOM, only : dclev
      IMPLICIT NONE

      integer, intent(in) :: LBASE_MIN,LBASE_MAX
      real*8, intent(in) :: dtime  ! dummy variable
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT
      REAL*8, DIMENSION(LM) :: DP
      COMMON/WORK2/UT,VT,DP
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      LOGICAL POLE
      INTEGER I,J,L,K,IMAX,KMAX,IM1,LMAX,LMIN

      DOUBLE PRECISION, DIMENSION(NMOM) :: TMOMS,QMOMS
      REAL*8 DOK,PIJBOT,PIJ,PKMS,THPKMS,QMS
     *     ,TVMS,THETA,RDP,THM
     *     ,rvx

      if(.not. vt_on) then
          rvx=0.
      else
          rvx=deltx
      endif

      if(LBASE_MAX.GE.LM) stop 'DRYCNV: LBASE_MAX.GE.LM'
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      UT=U ; VT=V
C**** OUTSIDE LOOPS OVER J AND I
      JLOOP: DO J=1,JM
      POLE=.FALSE.
      IF (J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.

      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)
C****
C**** MAIN LOOP
C****
      IM1=IM
      ILOOP: DO I=1,IMAX
         DO K=1,KMAX
            RA(K)=RAVJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
      LMAX=LBASE_MIN-1
      lbase_loop: do while(lmax.lt.lbase_max)
      LMIN=LMAX+1
      LMAX=LMIN
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)) cycle lbase_loop
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS
C**** MIX THROUGH TWO LOWER LAYERS
      PIJBOT=PLIJ(LMIN,I,J)
      DP(LMIN)=PDSIG(LMIN,I,J)
      PIJ=PLIJ(LMIN+1,I,J)
      DP(LMIN+1)=PDSIG(LMIN+1,I,J)
      PKMS=PK(LMIN,I,J)*DP(LMIN)+PK(LMIN+1,I,J)*DP(LMIN+1)
      THPKMS=T(I,J,LMIN)*(PK(LMIN,I,J)*DP(LMIN))
     *  +T(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMS=Q(I,J,LMIN)*DP(LMIN)+Q(I,J,LMIN+1)*DP(LMIN+1)
C**** sum moments to mix over unstable layers
      TMOMS(XYMOMS) =
     &     TMOM(XYMOMS,I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     &     TMOM(XYMOMS,I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMOMS(XYMOMS) =
     &     QMOM(XYMOMS,I,J,LMIN  )*(DP(LMIN  ))  +
     &     QMOM(XYMOMS,I,J,LMIN+1)*(DP(LMIN+1))
      IF (LMIN+1.GE.LM) GO TO 150
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO L=LMIN+2,LM
        IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 160
        PIJ=PLIJ(L,I,J)
        DP(L)=PDSIG(L,I,J)
        PKMS=PKMS+(PK(L,I,J)*DP(L))
        THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*DP(L))
        QMS=QMS+Q(I,J,L)*DP(L)
        TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*DP(L))
        TMOMS(XYMOMS) = TMOMS(XYMOMS) +
     &       TMOM(XYMOMS,I,J,L)*(PK(L,I,J)*DP(L))
        QMOMS(XYMOMS) = QMOMS(XYMOMS) +
     &       QMOM(XYMOMS,I,J,L)*DP(L)
        THETA=TVMS/PKMS
      END DO
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PIJBOT*SIGE(LMIN)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      DO L=LMIN,LMAX
         AJL(J,L,JL_TRBHR)=AJL(J,L,JL_TRBHR)+
     &        (THM-T(I,J,L))*PK(L,I,J)*PLIJ(L,I,J)
         AJL(J,L,JL_TRBDLHT)=AJL(J,L,JL_TRBDLHT)+
     &        (QMS-Q(I,J,L))*PDSIG(L,I,J)*LHE/SHA
      T(I,J,L)=THM
      TMOM(XYMOMS,I,J,L)=TMOMS(XYMOMS)/PKMS
      TMOM(ZMOMS,I,J,L)=0.
      Q(I,J,L)=QMS
      QMOM(XYMOMS,I,J,L)=QMOMS(XYMOMS)*RDP
      QMOM(ZMOMS,I,J,L)=0.
      END DO
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      DO L=LMIN,LMAX
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*DP(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*DP(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      DO L=LMIN,LMAX
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c the following line gives bytewise different ajl
            AJL(IDJ(K),L,JL_DAMDC)=AJL(IDJ(K),L,JL_DAMDC)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PLIJ(L,I,J)*RA(K)
         ENDDO
      ENDDO
      enddo lbase_loop
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
      if(lbase_min.eq.1) then ! was called from surfce
         DCLEV(I,J)=LMAX
      endif
      IM1=I
      ENDDO ILOOP
      ENDDO JLOOP
      RETURN
      END SUBROUTINE DIFFUS

