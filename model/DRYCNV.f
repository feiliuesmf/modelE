#include "rundeck_opts.h"

      SUBROUTINE ATM_DIFFUS(LBASE_MIN,LBASE_MAX,DTIME)
!@sum  ATM_DIFFUS(DRYCNV) mixes air caused by dry convection.
!@+    this version checks base layers lbase_min to lbase_max.
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : lhe,sha,deltx
      USE MODEL_COM
      USE GEOM
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
#ifdef TRACERS_ON
      USE TRACER_COM, only: TRM,TRMOM,NTM
      USE TRACER_DIAG_COM, only: TAJLN,JLNT_TURB
#endif
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

#ifdef TRACERS_ON
      DOUBLE PRECISION, DIMENSION(NMOM,NTM) :: TRMOMS
      DOUBLE PRECISION, DIMENSION(     NTM) :: TRMS
      REAL*8 SDPL,BYSDPL
#endif

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
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*deltx).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*deltx)) cycle lbase_loop
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
#ifdef TRACERS_ON
      TRMS(:) = TRM(I,J,LMIN,:)+TRM(I,J,LMIN+1,:)
      TRMOMS(XYMOMS,:) = 
     &     TRMOM(XYMOMS,I,J,LMIN,:)+TRMOM(XYMOMS,I,J,LMIN+1,:)
#endif
      IF (LMIN+1.GE.LM) GO TO 150
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*deltx)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*deltx)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO L=LMIN+2,LM
        IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*deltx)) GO TO 160
        PIJ=PLIJ(L,I,J)
        DP(L)=PDSIG(L,I,J)
        PKMS=PKMS+(PK(L,I,J)*DP(L))
        THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*DP(L))
        QMS=QMS+Q(I,J,L)*DP(L)
        TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*deltx)*(PK(L,I,J)*DP(L))
        TMOMS(XYMOMS) = TMOMS(XYMOMS) +
     &       TMOM(XYMOMS,I,J,L)*(PK(L,I,J)*DP(L))
        QMOMS(XYMOMS) = QMOMS(XYMOMS) +
     &       QMOM(XYMOMS,I,J,L)*DP(L)
        THETA=TVMS/PKMS
#ifdef TRACERS_ON
      TRMS(:) = TRMS(:) + TRM(I,J,L,:)
      TRMOMS(XYMOMS,:) = TRMOMS(XYMOMS,:) + TRMOM(XYMOMS,I,J,L,:)
#endif
      END DO
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PIJBOT*SIGE(LMIN)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
#ifdef TRACERS_ON
        SDPL = 0.d0
        DO L=LMIN,LMAX
          SDPL = SDPL+DP(L)
        ENDDO
        BYSDPL = 1.D0/SDPL
#endif
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
#ifdef TRACERS_ON
        TAJLN(J,L,JLNT_TURB,:)=TAJLN(J,L,JLNT_TURB,:) +
     &     (TRMS(:)*(DP(L)*BYSDPL)-TRM(I,J,L,:))
      TRM(I,J,L,:) = TRMS(:)*(DP(L)*BYSDPL)
      TRMOM(XYMOMS,I,J,L,:) = TRMOMS(XYMOMS,:)*(DP(L)*BYSDPL)
      TRMOM(ZMOMS,I,J,L,:) = 0.
#endif
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
      END SUBROUTINE ATM_DIFFUS


      subroutine apply_fluxes_to_atm
!@sum applies earth fluxes to the first layer of the atmosphere
!@auth Original Development Team
!@ver  1.0
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj,siniv,cosiv
      USE FLUXES, only : dth1,dq1,du1,dv1
      USE MODEL_COM, only : im,jm,u,v,t,q
      implicit none
      integer i,j,k,imax,kmax
      real*8 hemi

      do j=1,jm
        imax=imaxj(j)
        do i=1,imax
          t(i,j,1)=  t(i,j,1)+dth1(i,j)
          q(i,j,1)=  q(i,j,1)+dq1(i,j)
        end do
      end do
c****
c**** add in surface friction to first layer wind
c****
c**** polar boxes
      do j=1,jm,jm-1
        imax=imaxj(j)
        kmax=kmaxj(j)
        hemi=1.
        if(j.le.jm/2) hemi=-1.
        do i=1,imax
        do k=1,kmax
          u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*(du1(i,j)*cosiv(k)+dv1(i,j)*siniv(k)*hemi)
          v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*(dv1(i,j)*cosiv(k)-du1(i,j)*siniv(k)*hemi)
        end do
        end do
      end do
c**** non polar boxes
      do j=2,jm-1
        imax=imaxj(j)
        kmax=kmaxj(j)
        do i=1,imax
        do k=1,kmax
          u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*du1(i,j)
          v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*dv1(i,j)
        end do
        end do
      end do
c****
      return
      end subroutine apply_fluxes_to_atm
