#include "rundeck_opts.h"

      SUBROUTINE ATM_DIFFUS(LBASE_MIN,LBASE_MAX,DTIME)
!@sum  ATM_DIFFUS(DRYCNV) mixes air caused by dry convection.
!@+    this version checks base layers lbase_min to lbase_max.
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : lhe,sha,deltx
      USE MODEL_COM
      USE DOMAIN_DECOMP, only : grid
      USE DOMAIN_DECOMP, only : halo_update,NORTH,checksum
      USE GEOM
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
      USE DAGCOM, only : ajl,jl_trbhr,jl_damdc,jl_trbdlht
#ifdef TRACERS_ON
      USE TRACER_COM, only: TRM,TRMOM,NTM
      USE TRACER_DIAG_COM, only: TAJLN,JLNT_TURB
#endif
      USE DYNAMICS, only : pk,pdsig,plij,dke
      USE PBLCOM, only : dclev,w2gcm,w2_l1
      IMPLICIT NONE

      integer, intent(in) :: LBASE_MIN,LBASE_MAX
      real*8, intent(in) :: dtime  ! dummy variable
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) :: 
     &        UT,VT
      REAL*8, DIMENSION(LM) :: DP
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      INTEGER I,J,L,K,IMAX,KMAX,IM1,LMAX,LMIN
C
      REAL*8  UKP1(IM,LM), VKP1(IM,LM), UKPJM(IM,LM),VKPJM(IM,LM)

!RKF: No halo needed by UKM or VMK. Not used at the poles.
      REAL*8  UKM(4,IM,grid%J_STRT_SKP:grid%J_STOP_SKP,LM)
      REAL*8  VKM(4,IM,grid%J_STRT_SKP:grid%J_STOP_SKP,LM)
      INTEGER  LRANG(2,IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
C
      REAL*8, DIMENSION(NMOM) :: TMOMS,QMOMS
      REAL*8 DOK,PIJBOT,PIJ,PKMS,QMS
     *     ,TVMS,THETA,RDP,THM

#ifdef TRACERS_ON
      REAL*8, DIMENSION(NMOM,NTM) :: TRMOMS
      REAL*8, DIMENSION(     NTM) :: TRMS
      REAL*8 SDPL,BYSDPL
#endif

      REAL*8, DIMENSION(IM,LM) :: DEL_U, DEL_V
      REAL*8  DEL_AJL(LM)

      INTEGER ::  J_1, J_0
      INTEGER ::  J_1H, J_0H
      INTEGER ::  J_1S, J_0S
      INTEGER ::  J_1STG, J_0STG

      J_0   = grid%J_STRT
      J_1   = grid%J_STOP

      J_0H  = grid%J_STRT_HALO
      J_1H  = grid%J_STOP_HALO

      J_0S  = grid%J_STRT_SKP 
      J_1S  = grid%J_STOP_SKP 

      J_0STG= grid%J_STRT_STGR
      J_1STG= grid%J_STOP_STGR

      if(LBASE_MAX.GE.LM) call stop_model('DRYCNV: LBASE_MAX.GE.LM',255)

      ! update w2gcm at 1st GCM layer
      w2gcm=0.d0
      do j=j_0,j_1
      do i=1,im
         w2gcm(1,i,j)=w2_l1(i,j)
      end do
      end do

C****
C**** Update north halos for arrays U and V
C****
      CALL CHECKSUM   (grid, U, __LINE__,__FILE__)
      CALL HALO_UPDATE(grid, U, from=NORTH)
      CALL CHECKSUM   (grid, V, __LINE__,__FILE__)
      CALL HALO_UPDATE(grid, V, FROM=NORTH)

C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      UT=U ; VT=V
C**** OUTSIDE LOOPS OVER J AND I
!$OMP  PARALLEL DO PRIVATE (I,IM1,IMAX,J,K,KMAX,L,LMIN,LMAX,IDI,IDJ,
!$OMP*   DP,PIJ,PIJBOT,PKMS, QMS,QMOMS, RA,RDP,
#ifdef TRACERS_ON
!$OMP*   TRMS,TRMOMS,SDPL,BYSDPL,
#endif
!$OMP*   THM,TVMS,THETA,TMOMS,UMS,VMS)
!$OMP*   SCHEDULE(DYNAMIC,2)
      JLOOP: DO J=J_0,J_1

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
C
         LRANG(1,I,J)=-1
         LRANG(2,I,J)=-2
C
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
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*deltx)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*deltx)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
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
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO L=LMIN+2,LM
        IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*deltx)) GO TO 160
        PIJ=PLIJ(L,I,J)
        DP(L)=PDSIG(L,I,J)
        PKMS=PKMS+(PK(L,I,J)*DP(L))
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
      QMS=QMS*RDP
      THM=TVMS/(PKMS*(1.+QMS*deltx))
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
      LRANG(1,I,J)=LMIN
      LRANG(2,I,J)=LMAX
c     DO L=LMIN,LMAX
c        DO K=1,KMAX
c           U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
c    &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
c           V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
c    &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c           AJL(IDJ(K),L,JL_DAMDC)=AJL(IDJ(K),L,JL_DAMDC)
c    &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PLIJ(L,I,J)*RA(K)
c        ENDDO
c     ENDDO
      DO L=LMIN,LMAX
         IF(J.EQ.1)  THEN
            DO K=1,KMAX
               UKP1(K,L)=(UMS(K)-UT(IDI(K),IDJ(K),L))
               VKP1(K,L)=(VMS(K)-VT(IDI(K),IDJ(K),L))
            END DO
         ELSE IF(J.EQ.JM)  THEN
            DO K=1,KMAX
               UKPJM(K,L)=(UMS(K)-UT(IDI(K),IDJ(K),L))
               VKPJM(K,L)=(VMS(K)-VT(IDI(K),IDJ(K),L))
            END DO
         ELSE
            DO K=1,KMAX
               UKM(K,I,J,L)=(UMS(K)-UT(IDI(K),IDJ(K),L))
               VKM(K,I,J,L)=(VMS(K)-VT(IDI(K),IDJ(K),L))
            END DO
         END IF
      ENDDO
C
      enddo lbase_loop
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
      if(lbase_min.eq.1) then ! was called from surfce
         DCLEV(I,J)=LMAX
      endif
      IM1=I
      ENDDO ILOOP
      ENDDO JLOOP
!$OMP  END PARALLEL DO
C
C     NOW REALLY UPDATE THE MODEL WINDS
C
      IF (grid%HAVE_SOUTH_POLE) then
       J=1
       DO K=1,KMAXJ(J)
         IDI(K)=IDIJ(K,1,J)
         IDJ(K)=IDJJ(K,J)
         RA(K) =RAVJ(K,J)
       END DO
       LMIN=LRANG(1,1,J)
       LMAX=LRANG(2,1,J)
       DO L=LMIN,LMAX
        DO K=1,KMAXJ(J)
          U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKP1(K,L)*RA(K)
          V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKP1(K,L)*RA(K)
          AJL(IDJ(K),L,JL_DAMDC)=AJL(IDJ(K),L,JL_DAMDC)+
     *       UKP1(K,L)*PLIJ(L,1,J)*RA(K)
       END DO ; END DO
      END IF   !END SOUTH POLE
C
      DO J=J_0S, J_1-1       !J_1S computed below
        KMAX=KMAXJ(J)
        DO K=1,KMAX
           IDJ(K)=IDJJ(K,J)
           RA(K) =RAVJ(K,J)
        END DO
        DO I=1,IM
          LMIN=LRANG(1,I,J)
          LMAX=LRANG(2,I,J)
          DO L=LMIN,LMAX
          DO K=1,KMAX
            IDI(K)=IDIJ(K,I,J)
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)*RA(K)
            AJL(IDJ(K),L,JL_DAMDC)=AJL(IDJ(K),L,JL_DAMDC)+
     *            UKM(K,I,J,L)*PLIJ(L,I,J)*RA(K)
          END DO ; END DO
        END DO
      END DO
C
      IF (grid%HAVE_NORTH_POLE) THEN
       J=JM
       KMAX=KMAXJ(J)
       DO K=1,KMAX
         IDI(K)=IDIJ(K,1,J)
         IDJ(K)=IDJJ(K,J)
         RA(K) =RAVJ(K,J)
       END DO
       LMIN=LRANG(1,1,J)
       LMAX=LRANG(2,1,J)
       DO L=LMIN,LMAX
        DO K=1,KMAX
          U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKPJM(K,L)*RA(K)
          V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKPJM(K,L)*RA(K)
          AJL(IDJ(K),L,JL_DAMDC)=AJL(IDJ(K),L,JL_DAMDC)+
     *        UKPJM(K,L)*PLIJ(L,1,J)*RA(K)
       END DO ; END DO

      ELSE
C**** Loop cycle for j=j_1 for internal blocks
        J=J_1
C***  Initialize dummy work arrays
        DEL_U(1:IM,1:LM)=0.
        DEL_V(1:IM,1:LM)=0.
        DEL_AJL(1:LM)=0.

        KMAX=KMAXJ(J)
        DO K=1,KMAX
           IDJ(K)=IDJJ(K,J)
           RA(K) =RAVJ(K,J)
        END DO
        DO I=1,IM
          LMIN=LRANG(1,I,J)
          LMAX=LRANG(2,I,J)
          DO L=LMIN,LMAX
            DO K=1,2
              IDI(K)=IDIJ(K,I,J)
              U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)*RA(K)
              V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)*RA(K)
              AJL(IDJ(K),L,JL_DAMDC)=AJL(IDJ(K),L,JL_DAMDC)+
     *              UKM(K,I,J,L)*PLIJ(L,I,J)*RA(K)
            END DO 
            DO K=3,4
              IDI(K)=IDIJ(K,I,J)
              DEL_U(IDI(K),L) = DEL_U(IDI(K),L) + UKM(K,I,J,L)*RA(K)  
              DEL_V(IDI(K),L) = DEL_V(IDI(K),L) + VKM(K,I,J,L)*RA(K)  
              DEL_AJL(L)   = DEL_AJL(L) + UKM(K,I,J,L)*PLIJ(L,I,J)*RA(K)
            END DO
          END DO
        END DO
      ENDIF   !END NORTH POLE

C****EXCEPTIONS!!
C****Transfer the sums above: DEL_U DEL_V DEL_AJL to the block to the north
C*** and add them to U(I,J_0,L), V(I,J_0,L) and AJL(J_0,L,JL_DAMDC)
C*** respectively.
!PSEUDOCODE:
!     IF (.not. HAVE_NORTH_POLE ) 
!      ===>  send del_u del_v del_ajl to neighbor to the north
!     if (.not. HAVE_SOUTH_POLE) 
!      ===>  receive del_u_r del_v_r del_ajl_r from neighbor to the south 
!      ===>    U(1:IM,J_0,1:LM) = U(1:IM,J_0,1:LM) + DEL_U_R(1:IM,1:LM)
!              V(1:IM,J_0,1:LM) = V(1:IM,J_0,1:LM) + DEL_V_R(1:IM,1:LM)
!              AJL(J_0,1:LM,JL_DAMDC) = AJL(J_0,1:LM,JL_DAMDC) +
!                                                    DEL_AJL_R(L)
!END_PSEUDOCODE

C***

C**** Save additional changes in KE for addition as heat later
!$OMP  PARALLEL DO PRIVATE (L,I,J)
      DO L=1,LM
      DO J=J_0STG, J_1STG
      DO I=1,IM
        DKE(I,J,L)=DKE(I,J,L)+0.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L)
     *       -UT(I,J,L)*UT(I,J,L)-VT(I,J,L)*VT(I,J,L))
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

      RETURN
      END SUBROUTINE ATM_DIFFUS


      subroutine apply_fluxes_to_atm(dt)
!@sum applies earth fluxes to the first layer of the atmosphere
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,u,v,t,q,qcheck
      USE DOMAIN_DECOMP, only : grid, get
      USE DOMAIN_DECOMP, only : halo_update,NORTH,checksum
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj,siniv,cosiv,dxyp
      USE DYNAMICS, only : byam,am,dke
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,trm,trmom,trname,t_qlimit
#ifdef TRACERS_WATER
     *     ,trw0,tr_wd_TYPE,nWATER
#endif
      USE FLUXES, only : trflux1
#endif
      USE FLUXES, only : dth1,dq1,uflux1,vflux1,qflux1
      implicit none
      REAL*8, PARAMETER :: qmin=1.d-12
      integer i,j,k,n
      real*8, intent(in) :: dt
      real*8 hemi,trmin
      real*8, dimension(im,jm) :: usave,vsave
      INTEGER :: J_0,J_1,J_0S,J_1S,J_0STG,J_1STG
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

C****For distributed parallelization
            REAL*8, DIMENSION(IM) :: DEL_U, DEL_V

      del_u(1:im)=0.
      del_v(1:im)=0.
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE       )

      do j=j_0,j_1
        do i=1,imaxj(j)
          t(i,j,1) = t(i,j,1) + dth1(i,j)
          q(i,j,1) = q(i,j,1) + dq1(i,j)
        end do
      end do

#ifdef TRACERS_ON
      do n=1,ntm
        do j=j_0,j_1
          do i=1,imaxj(j)
            trm(i,j,1,n) = trm(i,j,1,n) + trflux1(i,j,n)*dt
            trmin=0.d0
#ifdef TRACERS_WATER
            IF(tr_wd_TYPE(n).eq.nWATER) trmin = 
     &      qmin*trw0(n)*am(1,i,j)*dxyp(j)
#endif
            if (t_qlimit(n).and.trm(i,j,1,n).lt.trmin) then
              if (qcheck) write(99,*) trname(n),I,J,' TR1:',
     *        trm(i,j,1,n),'->',trmin
              trm(i,j,1,n) = trmin
              trmom(:,i,j,1,n)=0.
            end if
          end do
        end do
      end do
#endif
c****
c**** add in surface friction to first layer wind
c****
      usave=u(:,:,1) ; vsave=v(:,:,1)
c**** SOUTH POLE BOX
      if (HAVE_SOUTH_POLE) then
      j=1
        hemi=-1.
        do i=1,imaxj(j)
        do k=1,kmaxj(j)
          u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *     ravj(k,j)*(uflux1(i,j)*cosiv(k)+vflux1(i,j)*siniv(k)*hemi)
     *     *dt*byam(1,I,J)
          v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *     ravj(k,j)*(vflux1(i,j)*cosiv(k)-uflux1(i,j)*siniv(k)*hemi)
     *     *dt*byam(1,I,J)
        end do
        end do
       end if  !SOUTH POLE

      IF (HAVE_NORTH_POLE) then
        j=jm
        hemi=1.
        do i=1,imaxj(j)
          do k=1,kmaxj(j)
            u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *       ravj(k,j)*(uflux1(i,j)*cosiv(k)+vflux1(i,j)*siniv(k)*hemi)
     *       *dt*byam(1,I,J)
            v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *       ravj(k,j)*(vflux1(i,j)*cosiv(k)-uflux1(i,j)*siniv(k)*hemi)
     *       *dt*byam(1,I,J)
          end do
        end do
      END IF   !NORTH POLE

c**** non polar boxes
      do j=J_0S,J_1-1
        do i=1,imaxj(j)
        do k=1,kmaxj(j)
          u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*uflux1(i,j)*dt*byam(1,I,J)
          v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*vflux1(i,j)*dt*byam(1,I,J)
        end do
        end do
      end do

C****For distr. parallelization: North-most lattitude of internal blocks.
        IF (.not. HAVE_NORTH_POLE) then
          j=j_1
          do i=1,imaxj(j)
            do k=1,2
              u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *               ravj(k,j)*uflux1(i,j)*dt*byam(1,I,J)
              v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *               ravj(k,j)*vflux1(i,j)*dt*byam(1,I,J)
            end do
            do k=3,4
              del_u(idij(k,i,j)) = del_u(idij(k,i,j)) -
     *               ravj(k,j)*uflux1(i,j)*dt*byam(1,I,J)
              del_v(idij(k,i,j)) = del_v(idij(k,i,j)) -
     *               ravj(k,j)*vflux1(i,j)*dt*byam(1,I,J)
            end do
          end do
        ENDIF   !.not. NORTH POLE


C****EXCEPTIONS!!
!****Transfer the sums above: DEL_U DEL_V to the block to the north
!*** and add them to U(I,J_0,L), V(I,J_0,L) respectively
!PSEUDOCODE:
!     IF (.not. HAVE_NORTH_POLE )
!      ===>  send del_u del_v to neighbor to the north
!     if (.not. HAVE_SOUTH_POLE)
!      ===>  receive del_u_r del_v_r from neighbor to the south
!      ===>    U(1:IM,J_0,1:LM) = U(1:IM,J_0,1:LM) + DEL_U_R(1:IM)
!              V(1:IM,J_0,1:LM) = V(1:IM,J_0,1:LM) + DEL_V_R(1:IM)
!END_PSEUDOCODE

C**** save change of KE for addition as heat later
      do j=J_0STG, J_1STG
        do i=1,im
          dke(i,j,1)=dke(i,j,1)+0.5*(u(i,j,1)*u(i,j,1)+v(i,j,1)*v(i,j,1)
     *         -usave(i,j)*usave(i,j)-vsave(i,j)*vsave(i,j))
        end do
      end do

c****
      return
      end subroutine apply_fluxes_to_atm
