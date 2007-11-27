#include "rundeck_opts.h"

      module ATMDYN
      implicit none
      private

      public init_ATMDYN, DYNAM,CALC_TROP,PGRAD_PBL
     &     ,DISSIP,FILTER,CALC_AMPK,CALC_AMP
     &     ,COMPUTE_DYNAM_AIJ_DIAGNOSTICS, COMPUTE_WSAVE
     &     ,AFLUX, CALC_PIJL, COMPUTE_MASS_FLUX_DIAGS
     &     ,getTotalEnergy,SDRAG
     &     ,addEnergyAsDiffuseHeat, addEnergyAsLocalHeat
#ifdef TRACERS_ON
     &     ,trdynam
#endif

      contains

      SUBROUTINE init_ATMDYN
      CALL AVRX
      end SUBROUTINE init_ATMDYN


      SUBROUTINE DYNAM
!@sum  DYNAM Integrate dynamic terms
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : by3,sha,mb2kg,rgas,bygrav
      USE MODEL_COM, only : im,jm,lm,u,v,t,p,q,wm,NIdyn,dt,MODD5K
     *     ,NSTEP,NDA5K,ndaa,mrch,ls1,byim,QUVfilter
      USE GEOM, only : dyv,dxv,dxyp,areag,bydxyp
      USE SOMTQ_COM, only : tmom,mz
      USE DYNAMICS, only : ptold,pu,pv,sd,phi,dut,dvt
     &    ,pua,pva,sda,ps,mb,pk,pmid,sd_clouds,pedn
      USE DIAG_COM, only : aij => aij_loc,ij_fmv,ij_fgzv
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE, GLOBALSUM
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      USE MOMENTS, only : advecv
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     PA, PB, PC, FPEU, FPEV
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     UT,VT,TT,TZ,TZT,MA,
     &     UX,VX,PIJL

      REAL*8 DTFS,DTLF
      INTEGER I,J,L,IP1,IM1   !@var I,J,L,IP1,IM1  loop variables
      INTEGER NS, NSOLD,MODDA    !? ,NIdynO

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

!?    NIdynO=MOD(NIdyn,2)   ! NIdyn odd is currently not an option
      DTFS=DT*2./3.
      DTLF=2.*DT
      NS=0
      NSOLD=0                            ! strat

!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         PUA(:,:,L) = 0.
         PVA(:,:,L) = 0.
         SDA(:,:,L) = 0.
      ENDDO
!$OMP  END PARALLEL DO
C**** Leap-frog re-initialization: IF (NS.LT.NIdyn)
  300 CONTINUE
c     UX(:,:,:)  = U(:,:,:)
c     UT(:,:,:)  = U(:,:,:)
c     VX(:,:,:)  = V(:,:,:)
c     VT(:,:,:)  = V(:,:,:)
!     copy z-moment of temperature into contiguous memory
c     tz(:,:,:) = tmom(mz,:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         UX(:,:,L)  = U(:,:,L)
         UT(:,:,L)  = U(:,:,L)
         VX(:,:,L)  = V(:,:,L)
         VT(:,:,L)  = V(:,:,L)
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
C
#ifdef NUDGE_ON
      CALL NUDGE_PREP
#endif
      PA(:,:) = P(:,:)
      PB(:,:) = P(:,:)
      PC(:,:) = P(:,:)

C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      MRCH=0
C     CALL DYNAM (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS)
      CALL CALC_PIJL(LM,P,PIJL)
      CALL AFLUX (U,V,PIJL)
      CALL ADVECM (P,PB,DTFS)
      CALL GWDRAG (PB,UX,VX,U,V,T,TZ,DTFS)   ! strat
#ifdef NUDGE_ON
      CALL NUDGE (UX,VX,DTFS)
#endif
      CALL VDIFF (PB,UX,VX,U,V,T,DTFS)       ! strat
      CALL ADVECV (P,UX,VX,PB,U,V,Pijl,DTFS)  !P->pijl
      CALL PGF (UX,VX,PB,U,V,T,TZ,Pijl,DTFS)
      if (QUVfilter) CALL FLTRUV(UX,VX,U,V)
C**** INITIAL BACKWARD STEP IS ODD, QT = Q + DT*F(QX)
      MRCH=-1
C     CALL DYNAM (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT)
      CALL CALC_PIJL(LS1-1,PB,PIJL)
      CALL AFLUX (UX,VX,PIJL)
      CALL ADVECM (P,PA,DT)
#ifdef NUDGE_ON
      CALL NUDGE  (UT,VT,DT)
#endif
      CALL GWDRAG (PA,UT,VT,UX,VX,T,TZ,DT)   ! strat
      CALL VDIFF (PA,UT,VT,UX,VX,T,DT)       ! strat
      CALL ADVECV (P,UT,VT,PA,UX,VX,Pijl,DT)   !PB->pijl
      CALL PGF (UT,VT,PA,UX,VX,T,TZ,Pijl,DT)

      if (QUVfilter) CALL FLTRUV(UT,VT,UX,VX)
      GO TO 360
C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)
  340 MRCH=-2
C     CALL DYNAM (UT,VT,TT,PT,QT,U,V,T,P,Q,DTLF)
      CALL CALC_PIJL(LS1-1,P,PIJL)
      CALL AFLUX (U,V,PIJL)
      CALL ADVECM (PA,PB,DTLF)
#ifdef NUDGE_ON
      CALL NUDGE  (UT,VT,DTLF)
#endif
      CALL GWDRAG (PB,UT,VT,U,V,T,TZ,DTLF)   ! strat
      CALL VDIFF (PB,UT,VT,U,V,T,DTLF)       ! strat
      CALL ADVECV (PA,UT,VT,PB,U,V,Pijl,DTLF)   !P->pijl
      CALL PGF (UT,VT,PB,U,V,T,TZ,Pijl,DTLF)
      if (QUVfilter) CALL FLTRUV(UT,VT,U,V)
      PA(:,:) = PB(:,:)     ! LOAD PB TO PA
C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)
  360 NS=NS+2
         MODD5K=MOD(NSTEP+NS-NIdyn+NDA5K*NIdyn+2,NDA5K*NIdyn+2)
      MRCH=2
C     CALL DYNAM (U,V,T,P,Q,UT,VT,TT,PT,QT,DTLF)
      CALL CALC_PIJL(LS1-1,PA,PIJL)
      CALL AFLUX (UT,VT,PIJL)
      CALL ADVECM (PC,P,DTLF)
#ifdef NUDGE_ON
      CALL NUDGE  (U,V,DTLF)
#endif
      CALL GWDRAG (P,U,V,UT,VT,T,TZ,DTLF)   ! strat
      CALL ADVECV (PC,U,V,P,UT,VT,Pijl,DTLF)     !PA->pijl
         MODDA=MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2)  ! strat
         IF(MODDA.LT.MRCH) CALL DIAGA0   ! strat
C**** ACCUMULATE MASS FLUXES FOR TRACERS and Q
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
        PUA(:,:,L)=PUA(:,:,L)+PU(:,:,L)
        PVA(:,:,L)=PVA(:,:,L)+PV(:,:,L)
        IF (L.LE.LM-1) SDA(:,:,L)=SDA(:,:,L)+SD(:,:,L)
      END DO
!$OMP  END PARALLEL DO
C**** ADVECT Q AND T
CCC   TT(:,:,:) = T(:,:,:)
CCC   TZT(:,:,:)= TZ(:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TT(:,:,L)  = T(:,:,L)
         TZT(:,:,L) = TZ(:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      call calc_amp(pc,ma)
C**** Update halo of PV from the north (needed in AADVT).
      CALL HALO_UPDATE(grid, PV, from=NORTH)
      CALL AADVT (MA,T,TMOM, SD,PU,PV, DTLF,.FALSE.,FPEU,FPEV)
!     save z-moment of temperature in contiguous memory for later
CCC   tz(:,:,:) = tmom(mz,:,:,:)
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TZ(:,:,L) = TMOM(MZ,:,:,L)
      ENDDO
!$OMP  END PARALLEL DO
      CALL VDIFF (P,U,V,UT,VT,T,DTLF)          ! strat
      PC(:,:)    = .5*(P(:,:)+PC(:,:))
CCC   TT(:,:,:)  = .5*(T(:,:,:)+TT(:,:,:))
CCC   TZT(:,:,:) = .5*(TZ(:,:,:)+TZT(:,:,:))
!$OMP PARALLEL DO PRIVATE (L)
      DO L=1,LM
         TT(:,:,L)  = .5*(T(:,:,L)+TT(:,:,L))
         TZT(:,:,L) = .5*(TZ(:,:,L)+TZT(:,:,L))
      ENDDO
!$OMP END PARALLEL DO

      CALL CALC_PIJL(LS1-1,PC,PIJL)
      CALL PGF (U,V,P,UT,VT,TT,TZT,Pijl,DTLF)    !PC->pijl

      call compute_mass_flux_diags(PHI, PU, PV, dt)

      CALL CALC_AMPK(LS1-1)
      if (QUVfilter) CALL FLTRUV(U,V,UT,VT)
      PC(:,:) = P(:,:)      ! LOAD P TO PC
      CALL SDRAG (DTLF)
         IF (MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2).LT.MRCH) THEN
           CALL DIAGA
           CALL DIAGB
           CALL EPFLUX (U,V,T,P)
         ENDIF
C**** Restart after 8 steps due to divergence of solutions
      IF (NS-NSOLD.LT.8 .AND. NS.LT.NIdyn) GO TO 340
      NSOLD=NS
      IF (NS.LT.NIdyn) GO TO 300

      PUA = PUA * DTLF
      PVA = PVA * DTLF
      SDA(:,:,1:LM-1) = SDA(:,:,1:LM-1) * DTLF

      RETURN
      END SUBROUTINE DYNAM

C**** Calculate 3D vertical velocity (take SDA which has units
C**** mb*m2/s (but needs averaging over no. of leap frog timesteps)
C**** and convert to WSAVE, units of m/s):

      subroutine COMPUTE_WSAVE(wsave, sda, T, PK, PEDN, NIdyn)
      use CONSTANT, only: rgas, bygrav
      !use MODEL_COM, only: NIdyn
      use DOMAIN_DECOMP, only: grid, GET
      use GEOM, only: bydxyp
      use MODEL_COM, only: IM,JM,LM

      real*8, dimension(:, grid%J_STRT_HALO:, :), intent(out) :: WSAVE
      real*8, dimension(:, grid%J_STRT_HALO:, :), intent(in)  :: SDA, T
      real*8, dimension(:, :, grid%J_STRT_HALO:), intent(in)  :: PK,PEDN
      integer, intent(in) :: NIdyn

      integer :: i, l
      integer :: J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)

!$OMP PARALLEL DO PRIVATE (l,i)
      do l=1,lm-1
        do i=1,im
         wsave(i,J_0:J_1,l)=2.*sda(i,J_0:J_1,l)*bydxyp(J_0:J_1)*
     &   rgas*0.5*(T(i,J_0:J_1,l)*pk(l,i,J_0:J_1)+T(i,J_0:J_1,l+1)*
     &   pk(l+1,i,J_0:J_1))*bygrav/(NIdyn*pedn(l+1,i,J_0:J_1))
        end do
      end do
!$OMP END PARALLEL DO

      end subroutine COMPUTE_WSAVE

      Subroutine compute_mass_flux_diags(PHI, PU, PV, dt)
      use MODEL_COM, only: IM, LM
      use DOMAIN_DECOMP, only: grid, halo_update, SOUTH, get
      use DIAG_COM, only: AIJ => AIJ_loc, IJ_FGZU, IJ_FGZV

      real*8, intent(inout) :: PHI(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: PU(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: PV(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt


      integer :: J_0S, J_1S
      integer :: J_0STG, J_1STG
      integer :: I, IP1, L, J

      call get(grid, J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG,
     &               J_STRT_SKP =J_0S,  J_STOP_SKP =J_1S)

      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE (J,L,I,IP1)
      DO J=J_0S,J_1S ! eastward transports
      DO L=1,LM
         I=IM
         DO IP1=1,IM
            AIJ(I,J,IJ_FGZU)=AIJ(I,J,IJ_FGZU)+
     &           (PHI(I,J,L)+PHI(IP1,J,L))*PU(I,J,L)*DT ! use DT=DTLF/2
            I=IP1
         END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE (J,L,I)
      DO J=J_0STG,J_1STG ! northward transports
      DO L=1,LM
         DO I=1,IM
            AIJ(I,J,IJ_FGZV)=AIJ(I,J,IJ_FGZV)+
     &           (PHI(I,J-1,L)+PHI(I,J,L))*PV(I,J,L)*DT ! use DT=DTLF/2
         END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

      end subroutine compute_mass_flux_diags

      subroutine COMPUTE_DYNAM_AIJ_DIAGNOSTICS(
     &    PUA, PVA, dt)
      use CONSTANT,      only: BY3
      use DOMAIN_DECOMP, only: grid, get, halo_update, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      use DIAG_COM, only: AIJ => AIJ_loc,
     &     IJ_FGZU, IJ_FGZV, IJ_FMV, IJ_FMU
      use MODEL_COM, only: IM,JM,LM

      real*8, intent(in) :: PUA(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: PVA(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      integer :: I, IP1, J, L
      integer :: J_0STG, J_1STG, J_0S, J_1S
      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
      real*8 :: dtlf

      dtlf = 2.*dt

      call get(grid, J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG,
     &               J_STRT_SKP =J_0S,  J_STOP_SKP =J_1S,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

!$OMP  PARALLEL DO PRIVATE (J,L)
      do j=J_0STG,J_1STG
      do L=1,LM
         AIJ(:,J,IJ_FMV)  = AIJ(:,J,IJ_FMV )+PVA(:,J,L)*DTLF
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO PRIVATE (J,L)
      do j=J_0S,J_1S
      do l=1,lm
         AIJ(:,J,IJ_FMU) = AIJ(:,J,IJ_FMU)+PUA(:,J,L)*DTLF
      enddo
      enddo
!$OMP  END PARALLEL DO

      if (haveLatitude(grid, J=1)) then
         do l=1,lm
            AIJ(:,1,IJ_FMU)  = AIJ(:, 1,IJ_FMU )+PUA(:, 1,L)*DTLF*BY3
         enddo
      endif
      if(haveLatitude(grid, J=JM)) then
         do l=1,lm
            AIJ(:,JM,IJ_FMU) = AIJ(:,JM,IJ_FMU )+PUA(:,JM,L)*DTLF*BY3
         enddo
      endif
      end subroutine COMPUTE_DYNAM_AIJ_DIAGNOSTICS


#ifdef TRACERS_ON
      SUBROUTINE TrDYNAM
!@sum  TrDYNAM is the driver to integrate tracer dynamic terms
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,itime,dt,byim
      USE TRACER_COM, only: itime_tr0,trm,trmom,trname,t_qlimit,ntm
      USE TRACER_ADV
      USE TRDIAG_COM, only: TAJLN=>TAJLN_loc, TAIJN=>TAIJN_LOC,
     *     jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,jlnt_vt_mm,
     *     tij_uflx,tij_vflx 
      IMPLICIT NONE
      REAL*8 DTLF,byncyc
      INTEGER N

C**** uses the fluxes pua,pva,sda from DYNAM and QDYNAM
      DO N=1,NTM
        IF (itime.LT.itime_tr0(N)) cycle
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
        safv = 0.; sbfv = 0. 

        CALL AADVQ (TRM(:,:,:,n),TrMOM(:,:,:,:,n),t_qlimit(n),trname(n))

C**** Flux diagnostics
        byncyc = 1./ncyc
        TAJLN(:,:,jlnt_nt_tot,n) = TAJLN(:,:,jlnt_nt_tot,n) + sbf(:,:)
        TAJLN(:,:,jlnt_nt_mm, n) = TAJLN(:,:,jlnt_nt_mm, n)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        TAJLN(:,:,jlnt_vt_tot,n) = TAJLN(:,:,jlnt_vt_tot,n) + scf(:,:)
        TAJLN(:,:,jlnt_vt_mm, n) = TAJLN(:,:,jlnt_vt_mm, n)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc

#ifdef TRACERS_WATER
C**** vertically integrated atmospheric fluxes
        TAIJN(:,:,tij_uflx,n) = TAIJN(:,:,tij_uflx,n) + safv(:,:)
        TAIJN(:,:,tij_vflx,n) = TAIJN(:,:,tij_vflx,n) + sbfv(:,:)
#endif

      ENDDO
      RETURN
      END SUBROUTINE TrDYNAM
#endif

      SUBROUTINE AFLUX (U,V,PIJL)
!@sum  AFLUX Calculates horizontal/vertical air mass fluxes
!@+    Input: U,V velocities, PIJL pressure
!@+    Output: PIT  pressure tendency (mb m^2/s)
!@+            SD   sigma dot (mb m^2/s)
!@+            PU,PV horizontal mass fluxes (mb m^2/s)
!@+            CONV  horizontal mass convergence (mb m^2/s)
!@+            SPA
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,imh,jm,lm,ls1,dsig,bydsig,byim
     &     ,zatmo,sige,do_polefix
      USE GEOM
      USE DYNAMICS, only : pit,sd,conv,pu,pv,sd_clouds,spa
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE
C**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
!@var U,V input velocities (m/s)
!@var PIJL input 3-D pressure field (mb) (no DSIG)
      REAL*8, INTENT(INOUT),
     &  DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) :: U,V
      REAL*8, INTENT(INOUT),
     &  DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) :: PIJL
      REAL*8, DIMENSION(IM) :: DUMMYS,DUMMYN
      INTEGER I,J,L,IP1,IM1,IPOLE
      REAL*8 PUS,PUN,PVS,PVN,PBS,PBN
      REAL*8 WT,DXDSIG,DYDSIG,PVSA(LM),PVNA(LM),xx,twoby3
      real*8, dimension(im,2,LM) :: usv0,vsv0
      integer :: jvs,jvn,jv
      real*8 :: wts
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****
C**** BEGINNING OF LAYER LOOP
C****
      CALL HALO_UPDATE(grid, U,FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid, V,FROM=NORTH+SOUTH)

c
c interpolate polar velocities to the appropriate latitude
c
      do ipole=1,2
        if (haveLatitude(grid, J=2) .and. ipole == 1) then
          jv  = 2
          jvs = 2          ! jvs is the southernmost velocity row
          jvn = jvs+1      ! jvs is the northernmost velocity row
          wts = polwt
        else if(haveLatitude(grid,JM) .and. ipole == 2) then
          jv = JM
          jvs = jv - 1
          jvn = jvs + 1
          wts = 1.-polwt
        else
          cycle
        endif

        do L = 1, LM
          usv0(:,ipole,l) = u(:,jv,l)
          vsv0(:,ipole,l) = v(:,jv,l)
          u(:,jv,l) = wts*u(:,jvs,l) + (1.-wts)*u(:,jvn,l)
          v(:,jv,l) = wts*v(:,jvs,l) + (1.-wts)*v(:,jvn,l)
        end do
      enddo

C****
C**** COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
C**** ARAKAWA'S SCHEME B             PV   U,V    SECONDARY GRID ROW
C****
C**** COMPUTE PU, THE WEST-EAST MASS FLUX, AT NON-POLAR POINTS
      do L = 1, LM
      DO 2154 J=J_0S,J_1S
      DO 2154 I=1,IM
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
      CALL AVRX (SPA(1,J_0H,L))
      I=IM
      DO 2166 J=J_0S,J_1S
      DYDSIG = 0.25D0*DYP(J)*DSIG(L)
      DO 2165 IP1=1,IM
      PU(I,J,L)=DYDSIG*SPA(I,J,L)*(PIJL(I,J,L)+PIJL(IP1,J,L))
 2165 I=IP1
 2166 CONTINUE
      end do
C**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX
      CALL HALO_UPDATE(grid, PIJL, FROM=SOUTH)
      do L = 1, LM
      IM1=IM
      DO 2172 J=J_0STG, J_1STG
      DXDSIG = 0.25D0*DXV(J)*DSIG(L)
      DO 2170 I=1,IM
      PV(I,J,L)=DXDSIG*(V(I,J,L)+V(IM1,J,L))*(PIJL(I,J,L)+PIJL(I,J-1,L))

 2170 IM1=I
 2172 CONTINUE
      end do

c restore uninterpolated values of u,v at the pole
      do ipole=1,2
        if (haveLatitude(grid, J=2) .and. ipole == 1) then
          jv = 2           ! why not staggered grid
        else if(haveLatitude(grid, J=JM) .and. ipole.eq.2) then
          jv = JM            ! why not staggered grid
        else
          cycle
        endif
        do L = 1, LM
          u(:,jv,l) = usv0(:,ipole,l)
          v(:,jv,l) = vsv0(:,ipole,l)
        end do
      enddo

C**** COMPUTE PU*3 AT THE POLES
      IF (haveLatitude(grid, J=1)) Then
        do L = 1, LM
        PUS=0.
        PVS=0.
        DO I=1,IM
          PUS=PUS+U(I,2,L)
          PVS=PVS+PV(I,2,L)
        END DO
        PUS=.25*DYP(2)*PUS*PIJL(1,1,L)*BYIM
        PVS=PVS*BYIM
        PVSA(L)=PVS
        DUMMYS(1)=0.
        DO I=2,IM
          DUMMYS(I)=DUMMYS(I-1)+(PV(I,2,L) -PVS)*BYDSIG(L)
        END DO
        PBS=0.
        PBN=0.
        DO I=1,IM
          PBS=PBS+DUMMYS(I)
        END DO
        PBS=PBS*BYIM
        DO I=1,IM
          SPA(I,1,L)=4.*(PBS-DUMMYS(I)+PUS)/(DYP(2)*PIJL(1,1,L))
          PU(I,1,L)=3.*(PBS-DUMMYS(I)+PUS)*DSIG(L)
        END DO
        end do
      END IF

      IF (haveLatitude(grid, J=JM)) THEN
        do L = 1,LM
        PUN=0.
        PVN=0.
        DO I=1,IM
          PUN=PUN+U(I,JM,L)
          PVN=PVN+PV(I,JM,L)
        END DO
        PUN=.25*DYP(JM-1)*PUN*PIJL(1,JM,L)*BYIM
        PVN=PVN*BYIM
        PVNA(L)=PVN
        DUMMYN(1)=0.
        DO I=2,IM
          DUMMYN(I)=DUMMYN(I-1)+(PV(I,JM,L)-PVN)*BYDSIG(L)
        END DO
        PBN=0.
        DO I=1,IM
          PBN=PBN+DUMMYN(I)
        END DO
        PBN=PBN*BYIM
        DO  I=1,IM
          SPA(I,JM,L)=4.*(DUMMYN(I)-PBN+PUN)/(DYP(JM-1)*PIJL(1,JM,L))
          PU(I,JM,L)=3.*(DUMMYN(I)-PBN+PUN)*DSIG(L)
        END DO
        END DO
      END IF
C****
C**** CONTINUITY EQUATION
C****
C**** COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
c     CALL HALO_UPDATE(grid, PV, FROM=NORTH)
c     DO 1510 J=J_0S,J_1S
c     IM1=IM
c     DO 1510 I=1,IM
c     CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
c1510 IM1=I
c     IF (HAVE_SOUTH_POLE) CONV(1,1,L)=-PVS
c     IF (HAVE_NORTH_POLE) CONV(1,JM,L)=PVN

      if(do_polefix.eq.1) then
c To maintain consistency with subroutine ADVECV,
c adjust pu at the pole if no corner fluxes are used there
c in ADVECV.
      do ipole=1,2
        if(haveLatitude(grid,J=1) .and. ipole.eq.1) then
          j = 1
        else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
          j = JM
        else
          cycle
        endif
        twoby3 = 2d0/3d0
        do l=1,lm
           pu(:,j,l) = pu(:,j,l)*twoby3
        enddo
      enddo
      endif
C****
C**** END OF HORIZONTAL ADVECTION LAYER LOOP
C****
c
c modify uphill air mass fluxes around steep topography
      do 2015 j=J_0S,J_1S
      i = im
      do 2010 ip1=1,im
         xx = zatmo(ip1,j)-zatmo(i,j)
         if(xx.eq.0.0)  go to 2007
         DO 2005 L=1,LS1-1
ccc         if((zatmo(ip1,j)-zatmo(i,j))*pu(i,j,l).gt.0.) then
            if(xx*pu(i,j,l).gt.0.) then
               if(pu(i,j,l).gt.0.) then
                  wt = (pijl(ip1,j,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j,l)/pijl(ip1,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.le.0.) then
                  pu(i,j,l+1) = pu(i,j,l+1) + pu(i,j,l)
                  pu(i,j,l) = 0.
               else
                  go to 2007
               endif
            endif
 2005    CONTINUE
 2007    CONTINUE
         i = ip1
 2010 CONTINUE
 2015 CONTINUE

ccc   do j=2,jm
c**** Exceptional J loop boundaries
      do 2035 j=max(J_0S,3), J_1S
      do 2035 i=1,im
         xx = zatmo(i,j)-zatmo(i,j-1)
         if(xx.eq.0.0)  go to 2035
         DO 2020 L=1,LS1-1
ccc         if((zatmo(i,j)-zatmo(i,j-1))*pv(i,j,l).gt.0.) then
            if(xx*pv(i,j,l).gt.0.) then
               if(pv(i,j,l).gt.0.) then
                  wt = (pijl(i,j,l)/pijl(i,j-1,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j-1,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.le.0.) then
                  pv(i,j,l+1) = pv(i,j,l+1) + pv(i,j,l)
                  pv(i,j,l) = 0.
               else
                  go to 2035
               endif
            endif
 2020    CONTINUE
 2035 CONTINUE
C
C     Now Really Do  CONTINUITY EQUATION
C
C     COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
C
      CALL HALO_UPDATE(grid, PV, FROM=NORTH)
!$OMP  PARALLEL DO PRIVATE (I,J,L,IM1)
      DO 2400 L=1,LM
      DO 1510 J=J_0S,J_1S
      IM1=IM
      DO 1510 I=1,IM
      CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
 1510 IM1=I
      IF (haveLatitude(grid,J=1)) CONV(1,1,L)=-PVSA(L)
      If (haveLatitude(grid,J=JM)) CONV(1,JM,L)=PVNA(L)
 2400 CONTINUE
!$OMP  END PARALLEL DO
C
C**** COMPUTE PIT, THE PRESSURE TENDENCY
      PIT(:,J_0:J_1) = CONV(:,J_0:J_1,1)
      SD(:,J_0:J_1,1:LM-1) = CONV(:,J_0:J_1,2:LM)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2420 J=J_0,J_1
         DO 2410 I=1,IMAXJ(J)
         DO 2410 L=LM-1,1,-1
           PIT(I,J) = PIT(I,J) + SD(I,J,L)
 2410    CONTINUE
 2420 CONTINUE
!$OMP  END PARALLEL DO
C**** COMPUTE SD, SIGMA DOT                                             -------
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2435 J=J_0,J_1
         DO 2430 I=1,IMAXJ(J)
         DO 2430 L=LM-2,LS1-1,-1
            SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)
 2430    CONTINUE
 2435 CONTINUE
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2440 J=J_0,J_1
         DO 2438 I=1,IMAXJ(J)
         DO 2438 L=LS1-2,1,-1
           SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)-
     &           DSIG(L+1)*PIT(I,J)
 2438    CONTINUE
 2440 CONTINUE
!$OMP  END PARALLEL DO
      DO 2450 L=1,LM-1
      DO 2450 I=2,IM
        IF (haveLatitude(grid,J=1)) 
     *     SD(I,1,L)=SD(1,1,L)
 2450   IF (haveLatitude(grid,J=JM))
     *       SD(I,JM,L)=SD(1,JM,L)

        ! Recopy into CONV to support prior usage
      CONV(:,J_0:J_1,1) = PIT(:,J_0:J_1)
      CONV(:,J_0:J_1,2:LM) = SD(:,J_0:J_1,1:LM-1)

      RETURN
      END SUBROUTINE AFLUX


      SUBROUTINE ADVECM (P,PA,DT1)
!@sum  ADVECM Calculates updated column pressures using mass fluxes
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,mrch,zatmo,u,v,t,q,ptop
      USE GEOM, only : bydxyp,imaxj
      USE DYNAMICS, only : pit
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE, GLOBALSUM
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: P(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: PA(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8, INTENT(IN) :: DT1
      INTEGER I,J,L  !@var I,J,L  loop variables
      INTEGER IM1 ! @var IM1 = I - 1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      INTEGER :: n_exception, n_exception_all

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1, J_STRT_HALO = J_0H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE )

C**** COMPUTE PA, THE NEW SURFACE PRESSURE
      ! 1st pass count warning/termination events
      ! This avoides the need for 2 halo fills during normal
      ! execution.
      n_exception = 0
      outer_loop:  DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          PA(I,J)=P(I,J)+(DT1*PIT(I,J)*BYDXYP(J))
          IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
            n_exception = n_exception + 1
            IF (PA(I,J)+PTOP.lt.250. .or. PA(I,J)+PTOP.GT.1200.)
     *           Exit outer_loop
          End If
        END DO
      END DO outer_loop

      Call GLOBALSUM(grid, n_exception, n_exception_all, all=.true.)
      IF (n_exception_all > 0) Then ! need halos
        CALL HALO_UPDATE(grid, U, FROM=NORTH)
        CALL HALO_UPDATE(grid, V, FROM=NORTH)
      END IF

      IF (n_exception > 0)  Then ! 2nd pass report problems
        Do J = J_0, J_1
          DO I = 1, IMAXJ(J)
            IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
              IM1 = 1 + MOD(I-1,IM)
              WRITE (6,990) I,J,MRCH,P(I,J),PA(I,J),ZATMO(I,J),DT1,
     *             (U(IM1,J,L),U(I,J,L),U(IM1,J+1,L),U(I,J+1,L),
     *             V(IM1,J,L),V(I,J,L),V(IM1,J+1,L),V(I,J+1,L),
     *             T(I,J,L),Q(I,J,L),L=1,LM)
              write(6,*) "Pressure diagnostic error"
              IF (PA(I,J)+PTOP.lt.250. .or. PA(I,J)+PTOP.GT.1200.)
     &          call stop_model('ADVECM: Pressure diagnostic error',11)
            END IF
          END DO
        END DO
      END IF

      IF (haveLatitude(grid, J=1)) PA(2:IM, 1)=PA(1,1)
      IF (haveLatitude(grid, J=JM)) PA(2:IM,JM)=PA(1,JM)

C****
      RETURN
  990 FORMAT (/'0PRESSURE DIAGNOSTIC     I,J,MRCH,P,PA=',3I4,2F10.2/
     *  '     ZATMO=',F10.3,' DT=',F6.1/
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
     &     ,do_polefix
      USE GEOM, only : imaxj,dxyv,dxv,dyv,dxyp,dyp,dxp,acor,acor2
      USE DYNAMICS, only : gz,pu,pit,phi,spa,dut,dvt
      USE DIAG, only : diagcd
      USE DOMAIN_DECOMP, Only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM):: U,V,T
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO):: FD,RFDUX
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  UT, VT, QT, P, SZ
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *  PB

      REAL*8 PKE(LS1:LM+1)
      REAL*8 DT4,DT1
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDU,PHIDN,FACTOR
      INTEGER I,J,L,IM1,IP1,IPOLE  !@var I,J,IP1,IM1,L,IPOLE loop variab.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      DT4=DT1/4.
      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO
C****
C**** VERTICAL DIFFERENCING
C****
!$OMP  PARALLEL DO PRIVATE (L)
      DO L=LS1,LM
      SPA(:,:,L)=0.
      END DO
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO PRIVATE(I,J,L,DP,P0,PIJ,PHIDN,TZBYDP,X,
!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        PIJ=P(I,J,1)
        PDN=PIJ+PTOP
        PKDN=PDN**KAPA
        PHIDN=ZATMO(I,J)
C**** LOOP OVER THE LAYERS
        DO L=1,LM
          PKPDN=PKDN*PDN
          PKPPDN=PKPDN*PDN
          IF(L.GE.LS1) THEN
            DP=DSIG(L)*PSFMPT
            BYDP=1./DP
            P0=SIG(L)*PSFMPT+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PSFMPT+PTOP
            PKUP=PKE(L+1)
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
          ELSE
            DP=DSIG(L)*PIJ
            BYDP=1./DP
            P0=SIG(L)*PIJ+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PIJ+PTOP
            PKUP=PUP**KAPA
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
C****   CALCULATE SPA, MASS WEIGHTED THROUGHOUT THE LAYER
            SPA(I,J,L)=RGAS*((X+TZBYDP*PTOP)*(PKPDN-PKPUP)*BYKAPAP1
     *      -X*PTOP*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP
          END IF
C**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
          PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1
     *      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP*BYKAPAP1)
C**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
          PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP)
     *     *BYKAPAP1)
          PDN=PUP
          PKDN=PKUP
        END DO
      END DO
      END DO
!$OMP END PARALLEL DO
C**** SET POLAR VALUES FROM THOSE AT I=1
      IF (haveLatitude(grid, J=1)) THEN
        DO L=1,LM
          SPA(2:IM,1,L)=SPA(1,1,L)
          PHI(2:IM,1,L)=PHI(1,1,L)
        END DO
      END IF
      IF (haveLatitude(grid, J=JM)) THEN
        DO L=1,LM
          SPA(2:IM,JM,L)=SPA(1,JM,L)
          PHI(2:IM,JM,L)=PHI(1,JM,L)
        END DO
      END IF

!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
        GZ(:,:,L)=PHI(:,:,L)
      END DO
!$OMP END PARALLEL DO
C****
C**** PRESSURE GRADIENT FORCE
C****
C**** NORTH-SOUTH DERIVATIVE AFFECTS THE V-COMPONENT OF MOMENTUM
C
      CALL HALO_UPDATE(grid, P,   FROM=SOUTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, SPA, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE(I,IM1,J,L,FACTOR,FLUX)
      DO 3236 L=1,LM
      DO 3236 J=J_0STG,J_1STG
      FACTOR = DT4*DXV(J)*DSIG(L)
      IM1=IM
      DO 3234 I=1,IM
      FLUX=    ((P(I,J,L)+P(I,J-1,L))*(PHI(I,J,L)-PHI(I,J-1,L))+
     *  (SPA(I,J,L)+SPA(I,J-1,L))*(P(I,J,L)-P(I,J-1,L)))*FACTOR
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3234 IM1=I
 3236 CONTINUE
!$OMP  END PARALLEL DO
C
C**** SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
C
C Although PU appears to require a halo update, the halos
C of PHI, SPA, and P enable implementation without the additional halo.
C
!$OMP  PARALLEL DO PRIVATE(I,IP1,J,L,FACTOR)
      DO L=1,LM
        IF (haveLatitude(grid, J=1)) PU(:,1,L)=0.
        IF (haveLatitude(grid, J=JM)) PU(:,JM,L)=0.
        I=IM

        DO J=Max(2,J_0STG-1),J_1STG
          DO IP1=1,IM
            PU(I,J,L)=(P(IP1,J,L)+P(I,J,L))*(PHI(IP1,J,L)-PHI(I,J,L))+
     *           (SPA(IP1,J,L)+SPA(I,J,L))*(P(IP1,J,L)-P(I,J,L))
            I=IP1
          END DO
        END DO

        CALL AVRX (PU(1,J_0H,L),jrange=(/MAX(2,J_0H),MIN(JM-1,J_1H)/))

        DO J=J_0STG,J_1STG
          FACTOR = -DT4*DYV(J)*DSIG(L)
          DO I=1,IM
            DUT(I,J,L)=DUT(I,J,L)+FACTOR*(PU(I,J,L)+PU(I,J-1,L))
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

c correct for erroneous dxyv at the poles
      if(do_polefix.eq.1) then
         do ipole=1,2
            if(haveLatitude(grid,J=2) .and. ipole.eq.1) then
               j = 2
            else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
               j = JM
            else
               cycle
            endif
            dut(:,j,:) = dut(:,j,:)*acor
            dvt(:,j,:) = dvt(:,j,:)*acor2
         enddo
      endif
C
C**** CALL DIAGNOSTICS
      IF(MRCH.GT.0) THEN
         IF(MODD5K.LT.MRCH) CALL DIAG5D (6,MRCH,DUT,DVT)
         CALL DIAGCD (grid,3,U,V,DUT,DVT,DT1)
      ENDIF
C****
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
C****
      DO 3410 J=J_0STG,J_1STG
      DO 3410 I=1,IM
 3410 FD(I,J)=PB(I,J)*DXYP(J)
      IF (haveLatitude(grid, J=1)) THEN
        FDSP=PB(1, 1)*DXYP( 1)
        FDSP=FDSP+FDSP
        DO I=1,IM
          FD(I, 1)=FDSP
        END DO
      END IF
      IF (haveLatitude(grid, J=JM)) THEN
        FDNP=PB(1,JM)*DXYP(JM)
        FDNP=FDNP+FDNP
        DO I=1,IM
          FD(I,JM)=FDNP
        END DO
      END IF
C
      CALL HALO_UPDATE(grid, FD, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE(I,IP1,J)
      DO 3530 J=J_0STG,J_1STG
      I=IM
      DO 3525 IP1=1,IM
      RFDUX(I,J)=4./(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
 3525 I = IP1
 3530 CONTINUE
!$OMP  END PARALLEL DO
C
!$OMP  PARALLEL DO PRIVATE(I,J,L,RFDU)
      DO 3550 L=1,LM
      DO 3550 J=J_0STG,J_1STG
      RFDU=1./(PSFMPT*DXYV(J)*DSIG(L))
      DO 3540 I=1,IM
      IF(L.LT.LS1) RFDU=RFDUX(I,J)*BYDSIG(L)
      VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)*RFDU
      UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)*RFDU
 3540 CONTINUE
 3550 CONTINUE
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE PGF

      SUBROUTINE AVRX(X,jrange)
!@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only : im,jm,imh
      USE GEOM, only : dlon,dxp,dyp,bydyp
      !USE DYNAMICS, only : xAVRX
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
      USE DOMAIN_DECOMP, Only : grid, GET
      USE MOMENTS, only : moment_enq_order
      USE constant, only : byrt2
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), optional ::
     &     X(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      Integer, Intent(In), optional :: jrange(2)
      REAL*8, ALLOCATABLE, SAVE  :: DRAT(:)
      REAL*8, SAVE ::  BYSN(IMH)
      REAL*8, DIMENSION(0:IMH) :: AN,BN
CCC   INTEGER, SAVE :: NMIN(grid%J_STRT_HALO:grid%J_STOP_HALO)
CCC   INTEGER, SAVE :: IFIRST = 1
      INTEGER, ALLOCATABLE, SAVE :: NMIN(:)
      INTEGER J,N
      LOGICAL, SAVE :: init = .false.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, J0, J1
      REAL*8, SAVE :: xAVRX
      INTEGER order

      if ( present(X) ) goto 1000

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO = J_0H, J_STOP_HALO = J_1H,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
C
      call moment_enq_order(order)
      if ( order==4 ) then
        xAVRX = byrt2
      else if ( order==2 ) then
        xAVRX = 1.d0
      else
        call stop_model("unsupported scheme order in AVRX0",255)
      endif
C
      IF (.NOT. init) THEN
        init = .true.
C       CALL FFT0(IM)
        j0 = MAX(1,J_0H)
        j1 = MIN(JM,J_1H)
        ALLOCATE(DRAT(j0:j1), NMIN(j0:j1))
        DO N=1,IMH
          BYSN(N)=xAVRX/SIN(.5*DLON*N)
        END DO
        DO J=j0,j1
          DRAT(J) = DXP(J)*BYDYP(3)
          DO N=IMH,1,-1
            IF(BYSN(N)*DRAT(J) .GT.1.) THEN
              NMIN(J) = N+1
              EXIT
            ENDIF
          END DO
        END DO
      END IF
      RETURN
C****
!!!      ENTRY AVRX (X)
 1000 continue
C****

      If (Present(jrange)) Then
        j0 = jrange(1)
        j1 = jrange(2)
      Else
        CALL GET(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
        j0=J_0S
        j1=J_1S
      End If

      DO J=j0,j1
        IF (DRAT(J).GT.1) CYCLE
        CALL FFT (X(1,J),AN,BN)
        DO N=NMIN(J),IMH-1
          AN(N)=BYSN(N)*DRAT(J) * AN(N)
          BN(N)=BYSN(N)*DRAT(J) * BN(N)
        END DO
        AN(IMH) = BYSN(IMH)*DRAT(J) * AN(IMH)
        CALL FFTI(AN,BN,X(1,J))
      END DO

      RETURN
      END SUBROUTINE AVRX

      SUBROUTINE FILTER
!@sum  FILTER Performs 8-th order shapiro filter in zonal direction
!@auth Original development team
!@ver  1.0
!@calls SHAP1D
C****
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER
C****        3  SMOOTH P AND T
C****
      USE CONSTANT, only : bygrav,kapa,sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,ls1,t,p,q,wm,mfiltr,zatmo,ptop
     *     ,byim,sig,itime,psf,pmtop
      USE GEOM, only : areag,dxyp
      USE SOMTQ_COM, only : tmom,qmom
      USE DYNAMICS, only : pk, COS_LIMIT
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm,trname,ITIME_TR0,trm,trmom
#endif
      USE PBLCOM, only : tsavg
      USE DOMAIN_DECOMP, Only : grid, GET, GLOBALSUM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: X,Y
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *        POLD, PRAT
      REAL*8 PSUMO,PSUMN,PDIF,AKAP,PS,ZS
      REAL*8, EXTERNAL :: SLP
      INTEGER I,J,L,N  !@var I,J,L  loop variables
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: KEJ,PEJ
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      REAL*8 initialTotalEnergy, finalTotalEnergy

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      IF (MOD(MFILTR,2).NE.1) GO TO 200
C**** Initialise total energy (J/m^2)
      initialTotalEnergy = getTotalEnergy()
C****
C**** SEA LEVEL PRESSURE FILTER ON P
C****
!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,PS,ZS)
      DO J=J_0S,J_1S
        DO I=1,IM
          POLD(I,J)=P(I,J)      ! Save old pressure
          PS=P(I,J)+PTOP
          ZS=ZATMO(I,J)*BYGRAV
          X(I,J)=SLP(PS,TSAVG(I,J),ZS)
          Y(I,J)=X(I,J)/PS
        END DO
      END DO
!$OMP  END PARALLEL DO
      CALL SHAP1D (8,X)
      call isotropslp(x,COS_LIMIT)
!$OMP  PARALLEL DO PRIVATE(I,J,PSUMO,PSUMN,PDIF)
      DO J=J_0S,J_1S
        PSUMO=0.
        PSUMN=0.
        DO I=1,IM
          PSUMO=PSUMO+P(I,J)
          P(I,J)=X(I,J)/Y(I,J)-PTOP
C**** reduce large variations (mainly due to topography)
          P(I,J)=MIN(MAX(P(I,J),0.99d0*POLD(I,J)),1.01d0*POLD(I,J))
          PSUMN=PSUMN+P(I,J)
        END DO
        PDIF=(PSUMN-PSUMO)*BYIM
        DO I=1,IM
          P(I,J)=P(I,J)-PDIF
        END DO
      END DO
!$OMP  END PARALLEL DO
C**** Scale mixing ratios (incl moments) to conserve mass/heat
!$OMP  PARALLEL DO PRIVATE(I,J)
      DO J=J_0S,J_1S
        DO I=1,IM
          PRAT(I,J)=POLD(I,J)/P(I,J)
        END DO
      END DO
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO L=1,LS1-1
      DO J=J_0S,J_1S
      DO I=1,IM
         Q(I,J,L)= Q(I,J,L)*PRAT(I,J)
         T(I,J,L)= T(I,J,L)*PRAT(I,J)
        WM(I,J,L)=WM(I,J,L)*PRAT(I,J)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*PRAT(I,J)
        IF (PRAT(I,J).lt.1.) THEN
          TMOM(:,I,J,L)=TMOM(:,I,J,L)*PRAT(I,J)
        END IF
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO
#ifdef TRACERS_ON
C**** In general, only an air tracer is affected by the filter
C**** This fix conserves tracer concentration, BUT NOT MASS!
C****   It is not wanted for most tracers.
C**** Thus, this code slows model and should be removed if not wanted
C**** Instead of if(trname...) could use n=n_air and avoid the IF-test
C**** But if n_air=0 this will cause problems...
      do n=1,ntm
      if (trname(n).ne.'Air') cycle
!     if (itime.lt.itime_tr0(n)) cycle   !probably not needed
!$OMP  PARALLEL DO PRIVATE (I,J,L)
      DO L=1,LS1-1
        DO J=J_0S,J_1S
          DO I=1,IM
             trm(I,J,L,n)=  trm(I,J,L,n)/PRAT(I,J)
             trmom(:,I,J,L,n)=trmom(:,I,J,L,n)/PRAT(I,J)
      end do; end do; end do
!$OMP  END PARALLEL DO
      end do
#endif
      CALL CALC_AMPK(LS1-1)

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)

  200 IF (MFILTR.LT.2) RETURN
C****
C**** TEMPERATURE STRATIFICATION FILTER ON T
C****
      AKAP=KAPA-.205d0    ! what is this number?
!$OMP  PARALLEL DO PRIVATE (J,L,X,Y)
      DO L=1,LM
        IF(L.LT.LS1) THEN
          DO J=J_0S,J_1S
            Y(:,J)=(SIG(L)*P(:,J)+PTOP)**AKAP
            X(:,J)=T(:,J,L)*Y(:,J)
          END DO
          CALL SHAP1D (8,X)
          DO J=J_0S,J_1S
            T(:,J,L)=X(:,J)/Y(:,J)
          END DO
        ELSE
          DO J=J_0S,J_1S
            X(:,J)=T(:,J,L)
          END DO
          CALL SHAP1D (8,X)
          DO J=J_0S,J_1S
            T(:,J,L)=X(:,J)
          END DO
        END IF
      END DO
!$OMP  END PARALLEL DO
C
      RETURN
      END SUBROUTINE FILTER

      SUBROUTINE FLTRUV(U,V,UT,VT)
!@sum  FLTRUV Filters 2 gridpoint noise from the velocity fields
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : sha
      USE MODEL_COM, only : im,jm,lm,byim,mrch,dt,t,ang_uv
     *  ,DT_XUfilter,DT_XVfilter,DT_YVfilter,DT_YUfilter
     &  ,do_polefix
      USE GEOM, only : dxyn,dxys,idij,idjj,rapj,imaxj,kmaxj
      USE DYNAMICS, only : pdsig,pk, COS_LIMIT
      USE DIAG, only : diagcd
C**********************************************************************
C**** FILTERING IS DONE IN X-DIRECTION WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES.
C**********************************************************************
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     *     INTENT(INOUT) :: U,V
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     *     INTENT(IN) :: UT,VT
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *     DUT,DVT,USAVE,VSAVE,DKE
      REAL*8 X(IM),YV(max(2*JM,IM)),DP(IM)
      REAL*8 XUby4toN,XVby4toN,YVby4toN,YUby4toN
      REAL*8 :: DT1=0.
      INTEGER I,J,K,L,N,IP1  !@var I,J,L,N  loop variables
      REAL*8 YV2,YVJ,YVJM1,X1,XI,XIM1
      INTEGER, PARAMETER :: NSHAP=8  ! NSHAP MUST BE EVEN
      REAL*8, PARAMETER :: BY16=1./16., by4toN=1./(4.**NSHAP)
      REAL*8 angm,dpt,D2V,D2U
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      INTEGER :: II
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      USAVE=U ; VSAVE=V
      if (DT_XUfilter.gt.0.) then
        XUby4toN = (DT/DT_XUfilter)*by4toN
      else
        XUby4toN = 0.
      end if
      if (DT_XVfilter.gt.0.) then
        XVby4toN = (DT/DT_XVfilter)*by4toN
      else
        XVby4toN = 0.
      end if
C****
C**** Filtering in east-west direction
C****
!$OMP  PARALLEL DO PRIVATE (I,J,L,N,X,X1,XI,XIM1)
      DO 350 L=1,LM
C**** Filter U component of velocity
      DO 240 J=J_0STG,J_1STG
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
  240 U(I,J,L) = U(I,J,L) - X(I)*XUby4toN
C**** Filter V component of velocity
      DO 340 J=J_0STG,J_1STG
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
  340 V(I,J,L) = V(I,J,L) - X(I)*XVby4toN
  350 CONTINUE
!$OMP  END PARALLEL DO
C****
C**** Filtering in north-south direction
C****

      IF (DT_YVfilter.gt.0.) THEN
      YVby4toN = (DT/DT_YVfilter)*by4toN
C**** Filter V component of velocity
      call halo_update(grid, V, from=NORTH)
      IF (haveLatitude(grid,J=2)) THEN
        J = 2
        DO L = 1, LM
          DO I = 1, IM/2
            II = I + IM/2
            D2V = V(I,J+1,L) - V(I,J,L) - V(I,J,L) - V(II,J,L)
            V(I, J,L) = V(I, J,L) - YVby4toN * D2V
            V(II,J,L) = V(II,J,L) - YVby4toN * D2V
          END DO
        END DO
      END IF
      IF (haveLatitude(grid,J=JM)) THEN
        J = JM
        DO L = 1, LM
          DO I = 1, IM/2
            II = I + IM/2
            D2V = -V(II,J-1,L) - V(I,J,L) - V(I,J,L) + V(I,J,L)
            V(I, J,L) = V(I, J,L) - YVby4toN * D2V
            V(II,J,L) = V(II,J,L) - YVby4toN * D2V
          END DO
        END DO
      END IF

      CALL HALO_UPDATE(grid, V, FROM=SOUTH+NORTH)
!$OMP  PARALLEL DO PRIVATE (I,J,L,D2V)
      DO L=1,LM
        DO J=MAX(3,J_0S),J_1S
          DO I=1,IM
            D2V= V(I,J-1,L)-V(I,J,L)-V(I,J,L)+V(I,J+1,L)
            V(I,J,L)=V(I,J,L) - YVby4toN * D2V
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

C**** Filtering longitudes on opposite sides of the globe simultaneously
C**** This is designed for resolutions with 1/2 boxes at poles and must
C****   be re-thought for others.
!$OMP  PARALLEL DO PRIVATE (I,J,L,N,YV,YV2,YVJ,YVJm1)
      DO 650 L=1,LM
      DO 650 I=1,IM/2
      DO 610 J=J_0STG,J_1STG
      YV(J) = V(I,J,L)
      YV(2*JM+1-J) = -V(I+IM/2,J,L)
  610 CONTINUE

      DO 630 N=1,NSHAP
      YV2   = YV(2)
      YVJm1 = YV(2*JM-1)
      DO 620 J=2,2*JM-2
      YVJ   = YV(J)
      YV(J) = YVJm1-YVJ-YVJ+YV(J+1)
      YVJm1 = YVJ
  620 CONTINUE
      J=2*JM-1
      YV(J)= YVJm1-YV(J)-YV(J)+YV2
  630 CONTINUE

      DO 640 J=J_0STG,J_1STG
      V(I,J,L) = V(I,J,L) - YV(J)*YVby4toN
      V(I+IM/2,J,L) = V(I+IM/2,J,L) + YV(2*JM+1-J)*YVby4toN
  640 CONTINUE
  650 CONTINUE
!$OMP  END PARALLEL DO
      END IF

      IF (DT_YUfilter.gt.0.) THEN
      YUby4toN = (DT/DT_YUfilter)*by4toN
C**** Filter U component of velocity

      CALL HALO_UPDATE(grid, U, FROM=SOUTH+NORTH)
      IF (haveLatitude(grid, J=2)) THEN
        J = 2
        DO L = 1, LM
          DO I = 1, IM/2
            II = I + IM/2
            D2U = U(I,J+1,L) - U(I,J,L) - U(I,J,L) - U(II,J,L)
            U(I, J,L) = U(I, J,L) - YUby4toN * D2U
            U(II,J,L) = U(II,J,L) - YUby4toN * D2U
          END DO
        END DO
      END IF
      IF (haveLatitude(grid,J=JM)) THEN
        J = JM
        DO L = 1, LM
          DO I = 1, IM/2
            II = I + IM/2
            D2U = -U(II,J-1,L) - U(I,J,L) - U(I,J,L) + U(I,J,L)
            U(I, J,L) = U(I, J,L) - YUby4toN * D2U
            U(II,J,L) = U(II,J,L) - YUby4toN * D2U
          END DO
        END DO
      END IF

      CALL HALO_UPDATE(grid, U, FROM=SOUTH+NORTH)
!$OMP  PARALLEL DO PRIVATE (I,J,L,D2U)
      DO L=1,LM
        DO J=MAX(3,J_0S),J_1S
          DO I=1,IM
            D2U= U(I,J-1,L)-U(I,J,L)-U(I,J,L)+U(I,J+1,L)
            U(I,J,L)=U(I,J,L) - YUby4toN * D2U
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

C**** Filtering longitudes on opposite sides of the globe simultaneously
C**** This is designed for resolutions with 1/2 boxes at poles and must
C****   be re-thought for others.
!$OMP  PARALLEL DO PRIVATE (I,J,L,N,YV,YV2,YVJ,YVJm1)
      DO 750 L=1,LM
      DO 750 I=1,IM/2
      DO 710 J=J_0STG,J_1STG
      YV(J) = U(I,J,L)
      YV(2*JM+1-J) = -U(I+IM/2,J,L)
  710 CONTINUE

      DO 730 N=1,NSHAP
      YV2   = YV(2)
      YVJm1 = YV(2*JM-1)
      DO 720 J=2,2*JM-2
      YVJ   = YV(J)
      YV(J) = YVJm1-YVJ-YVJ+YV(J+1)
      YVJm1 = YVJ
  720 CONTINUE
      J=2*JM-1
      YV(J)= YVJm1-YV(J)-YV(J)+YV2
  730 CONTINUE

      DO 740 J=J_0STG,J_1STG
      U(I,J,L) = U(I,J,L) - YV(J)*YUby4toN
      U(I+IM/2,J,L) = U(I+IM/2,J,L) + YV(2*JM+1-J)*YUby4toN
  740 CONTINUE
  750 CONTINUE
!$OMP  END PARALLEL DO
      END IF

      if(do_polefix.eq.1) then
         call isotropuv(u,v,COS_LIMIT)
c         if(have_south_pole) call isotropuv(u,v,-1)
c         if(have_north_pole) call isotropuv(u,v,+1)
      endif

C**** Conserve angular momentum along latitudes
c***  The following halo is not needed because PDSIG halo is up to date
c***      CALL HALO_UPDATE_COLUMN(grid, PDSIG, FROM=SOUTH)
!$OMP  PARALLEL DO PRIVATE (I,IP1,J,L,DP,ANGM,DPT)
      DO L=1,LM
        DO J=J_0STG,J_1STG
          ANGM=0.
          DPT=0.
          I=IM
          DO IP1=1,IM
            DP(I)=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *           +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
            ANGM=ANGM-DP(I)*(U(I,J,L)-USAVE(I,J,L))
            DPT=DPT+DP(I)
            I=IP1
          END DO
          DO I=1,IM
            if (ang_uv.eq.1) U(I,J,L)=U(I,J,L)+ANGM/DPT
            DKE(I,J,L)=0.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L)
     *           -USAVE(I,J,L)*USAVE(I,J,L)-VSAVE(I,J,L)*VSAVE(I,J,L))
            DUT(I,J,L)=(U(I,J,L)-USAVE(I,J,L))*DP(I)
            DVT(I,J,L)=(V(I,J,L)-VSAVE(I,J,L))*DP(I)
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO

C**** Call diagnostics and KE dissipation only for even time step
      IF (MRCH.eq.2) THEN
        CALL DIAGCD(grid,5,UT,VT,DUT,DVT,DT1)
        call addEnergyAsLocalHeat(DKE, T, PK)
      END IF

      RETURN
      END SUBROUTINE FLTRUV

      subroutine isotropslp(slp,coscut)
      use MODEL_COM, only : im,jm,dt
      USE DOMAIN_DECOMP, Only : GET,grid
      use GEOM, only : cosp,dxp
      implicit none
      real*8, parameter :: k=1d3
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: slp
      real*8 :: coscut,fac
      integer :: ipole,j,jcut,jinc,jp
      integer :: j_0s, j_1s, hemi

      CALL GET(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      Do j = j_0s, j_1s
        If (far_from_pole(j, cosp, coscut)) Cycle
        fac = k*dt/(dxp(j)*dxp(j))
        Call shap1(slp(1,j),im,fac)
      enddo

      return
      end subroutine isotropslp

      subroutine isotropuv(u,v,coscut)
!@sum  isotropuv isotropizes the velocity field in the near-polar row(s)
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : im,imh,jm,lm,dt
      USE DOMAIN_DECOMP, Only : GET, grid
      USE GEOM, only : cosv,dxv,cosi=>cosiv,sini=>siniv
      implicit none
      real*8, parameter :: klo=1d3,khi=1d7
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  U, V
      real*8 :: coscut,fac,k
      real*8, dimension(im) :: ua,va
      real*8, dimension(0:imh) :: an,bn
      integer :: i,j,l,hemi,jp,jinc,jcut,ipole
      Integer :: J_0STG, J_1STG

      CALL GET(grid, J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)

      do l=1,lm
        do j = J_0STG, J_1STG
          hemi = Hemisphere(j)
          If (far_from_pole(j, cosv, coscut)) Cycle

c compute xy velocities
          do i=1,im
            ua(i) = cosi(i)*u(i,j,l)-hemi*sini(i)*v(i,j,l)
            va(i) = cosi(i)*v(i,j,l)+hemi*sini(i)*u(i,j,l)
          enddo
c filter the xy velocities

          k = maxval(abs(u(:,j,l)))*2.*dt/dxv(j)
          if(k.lt.0.5) then
            k = klo
          else if(k.gt.1d0) then
            k = khi
          else
            k = klo + 2d0*(k-0.5)*(khi-klo)
          endif
          fac = k*dt/(dxv(j)*dxv(j))
          call shap1(ua,im,fac)
          call shap1(va,im,fac)
          if(at_pole(j)) then   ! really strong filtering right at the pole
            call fft(ua,an,bn)
            an(2:imh) = 0.
            bn(2:imh) = 0.
            call ffti(an,bn,ua)
            call fft(va,an,bn)
            an(2:imh) = 0.
            bn(2:imh) = 0.
            call ffti(an,bn,va)
          endif
c convert xy velocities back to polar coordinates
          do i=1,im
            u(i,j,l) = cosi(i)*ua(i)+hemi*sini(i)*va(i)
            v(i,j,l) = cosi(i)*va(i)-hemi*sini(i)*ua(i)
          enddo
        enddo                   ! j
      enddo                     ! l

      return
      end subroutine isotropuv

      Integer Function Hemisphere(j)
      Use GEOM, only: FJEQ
      Integer :: j

      If (J < FJEQ) Then
        hemisphere = -1
      Else
        hemisphere = +1
      End If
      End Function Hemisphere

      ! Detect whether at the pole on staggered grid
      Logical Function at_pole(j)
      Use model_com, only : jm
      Integer :: j
      If (j == jm .or. j == 2) Then
        at_pole = .true.
      else
        at_pole = .false.
      end if
      End Function at_pole

      Logical Function far_from_pole(j, cosj, coscut)
      Use MODEL_COM, only: JM
        Integer :: j
        Real*8 :: cosj(JM), coscut

        far_from_pole = (cosj(j) >= coscut)

      End Function far_from_pole

      subroutine shap1(x,im,fac)
      implicit none
      integer :: im
      real*8, dimension(im) :: x
      real*8 :: fac,facby4,x1,xim1,xi
      integer :: i,n,nn
      n = int(fac) + 1
      facby4 = fac*.25d0/n
      do nn=1,n
      x1 = x(1)
      xim1 = x(im)
      do i=1,im-1
         xi = x(i)
         x(i) = x(i) + facby4*(xim1-xi-xi+x(i+1))
         xim1 = xi
      enddo
      i = im
      x(i) = x(i) + facby4*(xim1-x(i)-x(i)+x1)
      enddo
      return
      end subroutine shap1

      SUBROUTINE SHAP1D (NORDER,X)
!@sum  SHAP1D Smoothes in zonal direction use n-th order shapiro filter
!@auth Original development team
!@ver  1.0
      USE MODEL_COM, only :im,jm
      USE DOMAIN_DECOMP, Only : grid, GET
      IMPLICIT NONE
!@var NORDER order of shapiro filter (must be even)
      INTEGER, INTENT(IN) :: NORDER
      REAL*8, INTENT(INOUT),
     *        DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: X
      REAL*8, DIMENSION(IM)::XS
      REAL*8 by4toN,XS1,XSIM1,XSI
      INTEGER I,J,N   !@var I,J,N  loop variables
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      by4toN=1./4.**NORDER
      DO J=J_0S,J_1S
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
      USE DOMAIN_DECOMP, Only : grid, GET
      implicit none
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: pijl
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
!@sum  CALC_AMPK calculate air mass and pressure arrays
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav,kapa
      USE MODEL_COM, only : im,jm,lm,ls1,p
      USE DYNAMICS, only : plij,pdsig,pmid,pk,pedn,pek,sqrtp,am,byam
      USE DOMAIN_DECOMP, Only : grid, GET, HALO_UPDATE, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE

      INTEGER :: I,J,L  !@var I,J,L  loop variables
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max. level for update
      REAL*8, DIMENSION(LMAX) :: PL,AML,PDSIGL,PMIDL
      REAL*8, DIMENSION(LMAX+1) :: PEDNL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_HALO= J_0H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Calculate air mass, layer pressures, P**K, and sqrt(P)
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure. Routine should be called with LMAX=LM at start, and
C**** subsequentaly with LMAX=LS1-1
C**** Note Air mass is calculated in (kg/m^2)

C**** Fill in polar boxes
      IF (haveLatitude(grid, J=1)) P(2:IM,1) = P(1,1)
      IF (haveLatitude(grid, J=JM)) P(2:IM,JM)= P(1,JM)
      Call HALO_UPDATE(grid, P, FROM=SOUTH)

!$OMP  PARALLEL DO PRIVATE (I,J,L,PL,AML,PDSIGL,PEDNL,PMIDL)
      DO J=J_0H,J_1 ! filling halo for P is faster than PDSIG
        DO I=1,IM

          CALL CALC_VERT_AMP(P(I,J),LMAX,PL,AML,PDSIGL,PEDNL,PMIDL)

          DO L=1,MIN(LMAX,LM)
            PLIJ (L,I,J) = PL    (L)
            PDSIG(L,I,J) = PDSIGL(L)
            PMID (L,I,J) = PMIDL (L)
            PEDN (L,I,J) = PEDNL (L)
            AM   (L,I,J) = AML   (L)
            PK   (L,I,J) = PMIDL (L)**KAPA
            PEK  (L,I,J) = PEDNL (L)**KAPA
            BYAM (L,I,J) = 1./AM(L,I,J)
          END DO

          IF (LMAX.ge.LM) THEN
            PEDN(LM+1:LMAX+1,I,J) = PEDNL(LM+1:LMAX+1)
            PEK (LM+1:LMAX+1,I,J) = PEDN(LM+1:LMAX+1,I,J)**KAPA
          END IF
          SQRTP(I,J) = SQRT(P(I,J))
        END DO
      END DO
!$OMP  END PARALLEL DO

      RETURN
      END SUBROUTINE CALC_AMPK

      SUBROUTINE CALC_AMP(p,amp)
!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
!@auth Jean Lerner/Max Kelley
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,ls1,dsig,psf,ptop
      USE GEOM, only : dxyp
C****
      USE DOMAIN_DECOMP, Only : grid, GET
      implicit none
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: amp
      integer :: j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

C
!$OMP  PARALLEL DO PRIVATE(J,L)
      DO L=1,LM
        IF(L.LT.LS1) THEN
ccc   do l=1,ls1-1
          do j=J_0,J_1
            amp(:,j,l) = p(:,j)*dxyp(j)*dsig(l)
          enddo
ccc   enddo
        ELSE
ccc   do l=ls1,lm
          do j=J_0,J_1
            amp(:,j,l) = (psf-ptop)*dxyp(j)*dsig(l)
          enddo
        END IF
ccc   enddo
      enddo
!$OMP  END PARALLEL DO
C
      return
C****
      end subroutine calc_amp

      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@auth Ye Cheng
!@ver  1.0
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE MODEL_COM, only : im,jm,t,p,zatmo,sig,byim
      USE GEOM, only : bydyp,bydxp,cosip,sinip
      USE DYNAMICS, only : phi,dpdy_by_rho,dpdy_by_rho_0,dpdx_by_rho
     *     ,dpdx_by_rho_0,pmid,pk
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE
      REAL*8 by_rho1,dpx1,dpy1,dpx0,dpy0,hemi
      INTEGER I,J,K,IP1,IM1,J1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)


C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the promary grids

      ! for dPdy/rho at non-pole grids
      CALL HALO_UPDATE(grid, P,   FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, ZATMO, FROM=SOUTH+NORTH)

      DO I=1,IM
        DO J=J_0S,J_1S
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDY_BY_RHO(I,J)=(100.*(P(I,J+1)-P(I,J-1))*SIG(1)*by_rho1
     2         +PHI(I,J+1,1)-PHI(I,J-1,1))*BYDYP(J)*.5d0
          DPDY_BY_RHO_0(I,J)=(100.*(P(I,J+1)-P(I,J-1))*by_rho1
     2         +ZATMO(I,J+1)-ZATMO(I,J-1))*BYDYP(J)*.5d0
        END DO
      END DO

      ! for dPdx/rho at non-pole grids

      DO J=J_0S,J_1S
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDX_BY_RHO(I,J)=(100.*(P(IP1,J)-P(IM1,J))*SIG(1)*by_rho1
     2         +PHI(IP1,J,1)-PHI(IM1,J,1))*BYDXP(J)*.5d0
          DPDX_BY_RHO_0(I,J)=(100.*(P(IP1,J)-P(IM1,J))*by_rho1
     2         +ZATMO(IP1,J)-ZATMO(IM1,J))*BYDXP(J)*.5d0
          IM1=I
          I=IP1
        END DO
      END DO

      ! at poles

      IF (haveLatitude(grid, J=1)) THEN
        hemi = -1.; J1 = 2
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,1)  =dpx1*BYIM
        DPDY_BY_RHO(1,1)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,1)=dpx0*BYIM
        DPDY_BY_RHO_0(1,1)=dpy0*BYIM
      END IF

      If (haveLatitude(grid, J=JM)) THEN
          hemi= 1.; J1=JM-1
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,JM)  =dpx1*BYIM
        DPDY_BY_RHO(1,JM)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,JM)=dpx0*BYIM
        DPDY_BY_RHO_0(1,JM)=dpy0*BYIM
      END IF

      END SUBROUTINE PGRAD_PBL

      SUBROUTINE SDRAG(DT1)
!@sum  SDRAG puts a drag on the winds in the top layers of atmosphere
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,sha
      USE MODEL_COM, only : im,jm,lm,ls1,u,v,t,q,x_sdrag,csdragl,lsdrag
     *     ,lpsdrag,ang_sdrag,itime,Wc_Jdrag
      USE GEOM, only : cosv,imaxj,kmaxj,idij,idjj,rapj,dxyv,dxyn,dxys
     *     ,rapvs,rapvn
      USE DIAG_COM, only : ajl=>ajl_loc,jl_dudtsdrg
      USE DYNAMICS, only : pk,pdsig,pedn
      USE DIAG, only : diagcd
      USE DOMAIN_DECOMP, only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      IMPLICIT NONE

!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var L(P)SDRAG lowest level at which SDRAG_lin is applied (near poles)
C**** SDRAG_const is applied above PTOP (150 mb) and below the SDRAG_lin
C**** regime (but not above P_CSDRAG)
      REAL*8 WL,TL,RHO,CDN,X,DP,DPL(LM),du,dps
!@var DUT,DVT change in momentum (mb m^3/s)
!@var DKE change in kinetic energy (m^2/s^2)
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *        DUT,DVT,DKE
      INTEGER I,J,L,IP1,K,Lmax
      logical cd_lin
!@var ang_mom is the sum of angular momentun at layers LS1 to LM
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)    ::
     *        ang_mom, sum_airm
!@param wmax imposed limit for stratospheric winds (m/s)
      real*8, parameter :: wmax = 200.d0 , wmaxp = wmax*3.d0/4.d0
      real*8 wmaxj,xjud
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      ang_mom=0. ;  sum_airm=0. ; dke=0. ; dut=0.
C*
      DUT=0. ; DVT=0.
c***  The following halo is not needed because PDSIG halo is up to date
c***      CALL HALO_UPDATE_COLUMN(grid, PDSIG, FROM=SOUTH)
      DO L=LS1,LM
      DO J=J_0STG, J_1STG
      cd_lin=.false.
      IF( L.ge.LSDRAG .or.
     *   (L.ge.LPSDRAG.and.COSV(J).LE..15) ) cd_lin=.true.
      wmaxj=wmax
      if(COSV(J).LE..15) wmaxj=wmaxp
      I=IM
      DO IP1=1,IM
        TL=T(I,J,L)*PK(L,I,J)   ! not quite correct - should be on UV grid
C**** check T to make sure it stayed within physical bounds
        if (TL.lt.100..or.TL.gt.373.) then
          write(99,*) 'SDRAG:',itime,i,j,l,'T,U,V=',TL,U(I,J,L),V(I,J,L)
          call stop_model('Stopped in ATMDYN::SDRAG',11)
        end if
        RHO=PEDN(L+1,I,J)/(RGAS*TL)   ! not quite correct - should be on UV grid
        WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
        xjud=1.
        if(Wc_JDRAG.gt.0.) xjud=(Wc_JDRAG/(Wc_JDRAG+min(WL,wmaxj)))**2
C**** WL is restricted to Wmax by adjusting X, if necessary;
C**** the following is equivalent to first reducing (U,V), if necessary,
C**** then finding the drag and applying it to the reduced winds
                    CDN=CSDRAGl(l)*xjud
        IF (cd_lin) CDN=(X_SDRAG(1)+X_SDRAG(2)*min(WL,wmaxj))*xjud
        DPS= (PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*RAPVN(J-1)+
     *       (PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*RAPVS(J)
        X=DT1*RHO*CDN*min(WL,wmaxj)*GRAV/DPS
        if (wl.gt.wmaxj) X = 1. - (1.-X)*wmaxj/wl
C**** adjust diags for possible difference between DT1 and DTSRC
        AJL(J,L,JL_DUDTSDRG) = AJL(J,L,JL_DUDTSDRG)-U(I,J,L)*X
        DP=DPS*DXYV(J)
        ang_mom(i,j) = ang_mom(i,j)+U(I,J,L)*X*DP
        DUT(I,J,L)=-X*U(I,J,L)*DP
        DVT(I,J,L)=-X*V(I,J,L)*DP
        DKE(I,J,L)=0.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*(X*X-2.*X)
        U(I,J,L)=U(I,J,L)*(1.-X)
        V(I,J,L)=V(I,J,L)*(1.-X)
        I=IP1
      END DO
      END DO
      END DO

C*
C***  Add the lost angular momentum uniformly back in if ang_sdrag>0
C***  only below 150mb if ang_sdrag=1, into whole column if ang_sdrag>1
C*
      if (ang_sdrag.gt.0) then
        lmax=ls1-1
        if (ang_sdrag.gt.1) lmax=lm
        do j = J_0STG,J_1STG
        I=IM
        do ip1 = 1,im
          do l = 1,lmax
            DPL(L)=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *        +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
            sum_airm(i,j) = sum_airm(i,j)+DPL(L)
          end do
C*
          do l = 1,lmax
            du = ang_mom(i,j)/sum_airm(i,j)
            DUT(I,J,L) = DUT(I,J,L) + du*dpl(l)
            AJL(J,L,JL_DUDTSDRG) = AJL(J,L,JL_DUDTSDRG) + du
            dke(i,j,l) = dke(i,j,l) + du*(u(i,j,l)+0.5*du)
            U(I,J,L)=U(I,J,L) + du
          end do
          I=IP1
        end do
        end do
      end if

      call addEnergyAsLocalHeat(DKE, T, PK)

C**** conservation diagnostic
C**** (technically we should use U,V from before but this is ok)
      CALL DIAGCD (grid,4,U,V,DUT,DVT,DT1)
      RETURN
      END SUBROUTINE SDRAG

      SUBROUTINE CALC_TROP
!@sum  CALC_TROP (to calculate tropopause height and layer)
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,t
      USE GEOM, only : imaxj
      USE DIAG_COM, only : aij => aij_loc, ij_ptrop, ij_ttrop
      USE DYNAMICS, only : pk, pmid, PTROPO, LTROPO
      USE DOMAIN_DECOMP, Only : grid, GET
      USE DOMAIN_DECOMP, only : haveLatitude
      IMPLICIT NONE
      INTEGER I,J,L,IERR
      REAL*8, DIMENSION(LM) :: TL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)


C**** Find WMO Definition of Tropopause to Nearest L
!$OMP  PARALLEL DO PRIVATE (I,J,L,TL,IERR)
      do j=J_0,J_1
      do i=1,imaxj(j)
        do l=1,lm
          TL(L)=T(I,J,L)*PK(L,I,J)
        end do
        CALL TROPWMO(TL,PMID(1,I,J),PK(1,I,J),PTROPO(I,J),LTROPO(I,J)
     *       ,IERR)
        IF (IERR.gt.0) print*,"TROPWMO error: ",i,j
        AIJ(I,J,IJ_PTROP)=AIJ(I,J,IJ_PTROP)+PTROPO(I,J)
        AIJ(I,J,IJ_TTROP)=AIJ(I,J,IJ_TTROP)+TL(LTROPO(I,J))
      end do
      end do
!$OMP  END PARALLEL DO
      IF (haveLatitude(grid, J=1)) THEN
        PTROPO(2:IM,1) = PTROPO(1,1)
        LTROPO(2:IM,1) = LTROPO(1,1)
      END IF
      IF (haveLatitude(grid,J=JM)) THEN
        PTROPO(2:IM,JM)= PTROPO(1,JM)
        LTROPO(2:IM,JM)= LTROPO(1,JM)
      END IF

      END SUBROUTINE CALC_TROP

      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE as heat locally
!@auth Gavin Schmidt
      USE MODEL_COM, only : t
      USE DYNAMICS, only : dke,pk
      IMPLICIT NONE
C**** DKE (m^2/s^2) is saved from surf,dry conv,aturb and m.c

      call addEnergyAsLocalHeat(DKE, T, PK)

      END SUBROUTINE DISSIP

      subroutine tropwmo(ptm1, papm1, pk, ptropo, ltropp,ierr)
!@sum  tropwmo calculates tropopasue height according to WMO formula
!@auth D. Nodorp/T. Reichler/C. Land
!@+    GISS Modifications by Jean Lerner/Gavin Schmidt
!@ver  1.0
!@alg  WMO Tropopause Definition
!@+
!@+ From A Temperature Lapse Rate Definition of the Tropopause Based on
!@+ Ozone, J. M. Roe and W. H. Jasperson, 1981
!@+
!@+ In the following discussion the lapse rate is defined as -dT/dz.
!@+
!@+ The main features of the WMO tropopause definition are as follows:
!@+ * The first tropopause (i.e., the conventional tropopause) is
!@+   defined as the lowest level at which the lapse rate decreases to 2
!@+   K/km or less, and the average lapse rate from this level to any
!@+   level within the next higher 2 km does not exceed 2 K/km.
!@+ * If above the first tropopause the average lapse rate between any
!@+   level and all higher levels within 1 km exceed 3 K/km, then a
!@+   second tropopause is defined by the same criterion as under the
!@+   statement above. This tropopause may be either within or above the
!@+   1 km layer.
!@+ * A level otherwise satisfying the definition of tropopause, but
!@+   occuring at an altitude below that of the 500 mb level will not be
!@+   designated a tropopause unless it is the only level satisfying the
!@+   definition and the average lapse rate fails to exceed 3 K/km over
!@+   at least 1 km in any higher layer.
!@+ * (GISS failsafe) Some cases occur when the lapse rate never falls
!@+   below 2 K/km. In such cases the failsafe level is that where the
!@+   lapse rate first falls below 3 K/km. If this still doesn't work
!@+   (ever?), the level is set to the pressure level below 30mb.
!@+
      USE MODEL_COM, only : klev=>lm
      USE CONSTANT, only : zkappa=>kapa,zzkap=>bykapa,grav,rgas
      USE DOMAIN_DECOMP, Only : grid, GET
      implicit none

      real*8, intent(in), dimension(klev) :: ptm1, papm1, pk
      real*8, intent(out) :: ptropo
      integer, intent(out) :: ltropp,ierr
      real*8, dimension(klev) :: zpmk, zpm, za, zb, ztm, zdtdz
!@param zgwmo min lapse rate (* -1) needed for trop. defn. (-K/km)
!@param zgwmo2 GISS failsafe minimum lapse rate (* -1) (-K/km)
!@param zdeltaz distance to check for lapse rate changes (km)
!@param zfaktor factor for caluclating height from pressure (-rgas/grav)
!@param zplimb min pressure at which to define tropopause (mb)
      real*8, parameter :: zgwmo  = -2d-3, zgwmo2=-3d-3,
     *     zdeltaz = 2000.0, zfaktor = -GRAV/RGAS, zplimb=500.
      real*8 zptph, zp2km, zag, zbg, zasum, zaquer, zptf
      integer iplimb,iplimt, jk, jj, kcount, ltset,l
      logical ldtdz
c****
c****  2. Calculate the height of the tropopause
c****  -----------------------------------------
      ltset = -999
      ierr=0
      iplimb=1
c**** set limits based on pressure
      do jk=2,klev-1
        if (papm1(jk-1).gt.600d0) then
          iplimb=jk
        else
          if (papm1(jk).lt.30d0) exit
        end if
      end do
      iplimt=jk
c****
c****  2.1 compute dt/dz
c****  -----------------
c****       ztm  lineare Interpolation in p**kappa
c****     gamma  dt/dp = a * kappa + papm1(jx,jk)**(kappa-1.)

      do jk=iplimb+1,iplimt       ! -1 ?????
        zpmk(jk)=0.5*(pk(jk-1)+pk(jk))

        zpm(jk)=zpmk(jk)**zzkap ! p mitte

        za(jk)=(ptm1(jk-1)-ptm1(jk))/(pk(jk-1)-pk(jk))
        zb(jk) = ptm1(jk)-(za(jk)*pk(jk))

        ztm(jk)=za(jk)*zpmk(jk)+zb(jk) ! T mitte
        zdtdz(jk)=zfaktor*zkappa*za(jk)*zpmk(jk)/ztm(jk)
      end do
c****
c****  2.2 First test: valid dt/dz ?
c****  -----------------------------
c****
      do 1000 jk=iplimb+1,iplimt-1

c**** GISS failsafe test
        if (zdtdz(jk).gt.zgwmo2.and.ltset.ne.1) then
          ltropp=jk
          ltset =1
        end if
c****
        if (zdtdz(jk).gt.zgwmo .and. ! dt/dz > -2K/km
     &       zpm(jk).le.zplimb) then ! zpm not too low
          ltropp = jk
          ltset = 1
c****
c****  2.3 dtdz is valid > something in German
c****  ----------------------------------------
c****    1.lineare in p^kappa (= Dieters neue Methode)

          zag = (zdtdz(jk)-zdtdz(jk+1))/
     &         (zpmk(jk)-zpmk(jk+1)) ! a-gamma
          zbg = zdtdz(jk+1) - zag*zpmk(jk+1) ! b-gamma
          if(((zgwmo-zbg)/zag).lt.0.) then
            zptf=0.
          else
            zptf=1.
          end if
          zptph = zptf*abs((zgwmo-zbg)/zag)**zzkap
          ldtdz=zdtdz(jk+1).lt.zgwmo
          if(.not.ldtdz) zptph=zpm(jk)
c****
c****  2.4 2nd test: dt/dz above 2km must not be lower than -2K/km
c****  -----------------------------------------------------------
c****
          zp2km = zptph + zdeltaz*zpm(jk)
     &         / ztm(jk)*zfaktor ! p at ptph + 2km
          zasum = 0.0           ! zdtdz above
          kcount = 0            ! number of levels above
c****
c****  2.5 Test until pm < p2km
c****  --------------------------
c****
          do jj=jk,iplimt-1
            if(zpm(jj).gt.zptph) cycle ! doesn't happen
            if(zpm(jj).lt.zp2km) goto 2000 ! ptropo valid
            zasum = zasum+zdtdz(jj)
            kcount = kcount+1
            zaquer = zasum/float(kcount) ! dt/dz mean
            if(zaquer.le.zgwmo) goto 1000 ! dt/dz above < 2K/1000
                                          ! discard it
          end do                ! test next level
          goto 2000
        endif
 1000 continue                  ! next level
 2000 continue

      if (ltset.eq.-999) then
        ltropp=iplimt-1  ! default = last level below 30mb
        print*,"In tropwmo ltropp not set, using default: ltropp ="
     *       ,ltropp
        write(6,'(12(I4,5F10.5,/))') (l,ptm1(l),papm1(l),pk(l),zdtdz(l)
     *       ,zpm(l),l=iplimb+1,iplimt-1)
        ierr=1
      end if
      ptropo = papm1(ltropp)
c****
      return
      end subroutine tropwmo

      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy returns the sum of kinetic and potential energy.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use GEOM, only: DXYP, AREAG
      use DOMAIN_DECOMP, only: grid, GLOBALSUM, get
      USE DOMAIN_DECOMP, only : haveLatitude
      REAL*8 :: totalEnergy

      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) ::KEJ,PEJ,TEJ
      REAL*8 :: totalPotentialEnergy
      REAL*8 :: totalKineticEnergy
      integer :: J_0, J_1
      logical :: HAVE_SOUTH_POLE

      call get(grid, J_STRT=J_0, J_STOP=J_1,
     *     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

      call conserv_PE(PEJ)
      call conserv_KE(KEJ)
      if (haveLatitude(grid, J=1)) KEJ(1) = 0

      TEJ(J_0:J_1)= (KEJ(J_0:J_1) + PEJ(J_0:J_1)*DXYP(J_0:J_1))/AREAG

      CALL GLOBALSUM(grid, TEJ, totalEnergy, ALL=.true.)

      end function getTotalEnergy

      subroutine addEnergyAsDiffuseHeat(deltaEnergy)
!@sum  addEnergyAsDiffuseHeat adds in energy increase as diffuse heat.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use CONSTANT, only: sha, mb2kg
      use MODEL_COM, only: T, PSF, PMTOP, LM
      use DYNAMICS, only: PK
      use DOMAIN_DECOMP, only: grid, get
      real*8, intent(in) :: deltaEnergy

      real*8 :: ediff
      integer :: l
      integer :: J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)

      ediff = deltaEnergy / ((PSF-PMTOP)*SHA*mb2kg)
!$OMP  PARALLEL DO PRIVATE (L)
      do l=1,lm
        T(:,J_0:J_1,L)=T(:,J_0:J_1,L)-ediff/PK(L,:,J_0:J_1)
      end do
!$OMP  END PARALLEL DO

      end subroutine addEnergyAsDiffuseHeat

C***** Add in dissipiated KE as heat locally
      subroutine addEnergyAsLocalHeat(deltaKE, T, PK, diagIndex)
!@sum  addEnergyAsLocalHeat adds in dissipated kinetic energy as heat locally.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use CONSTANT, only: SHA
      use GEOM, only: IDIJ, IDJJ, RAPJ, IMAXJ, KMAXJ
      use MODEL_COM, only: LM
      use DOMAIN_DECOMP, only: grid, get, HALO_UPDATE, NORTH
      use DIAG_COM, only: ajl => ajl_loc
      implicit none
      real*8 :: deltaKE(:,grid%j_strt_halo:,:)
      real*8 :: T(:,grid%j_strt_halo:,:)
      real*8 :: PK(:,:,grid%j_strt_halo:)
      integer, optional, intent(in) :: diagIndex

      integer :: i, j, k, l
      real*8 :: ediff
      integer :: J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      CALL HALO_UPDATE(grid, deltaKE, FROM=NORTH)

!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff,K)
      DO L=1,LM
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            ediff=0.
            DO K=1,KMAXJ(J)     ! loop over surrounding vel points
              ediff=ediff+deltaKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
            END DO
            ediff = ediff / (SHA*PK(L,I,J))
            T(I,J,L)=T(I,J,L)-ediff
            if (present(diagIndex)) then
              AJL(J,L,diagIndex) = AJL(J,L,diagIndex) - ediff
            end if
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO
      end subroutine addEnergyAsLocalHeat

      end module ATMDYN

      module ATMDYN_QDYNAM
      USE ATMDYN
      implicit none
      private

      public QDYNAM

      contains

      SUBROUTINE QDYNAM
!@sum  QDYNAM is the driver to integrate dynamic terms by the method
!@+          of pre-computing Courant limits using mean fluxes
!@+    It replaces CALL AADVT (MA,Q,QMOM, SD,PU,PV, DTLF,.TRUE.,
!@auth J. Lerner
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm,q,dt,byim
      USE SOMTQ_COM, only : qmom
      USE DIAG_COM, only: ajl=>ajl_loc,jl_totntlh,jl_zmfntlh,jl_totvtlh
     *     ,jl_zmfvtlh
      USE DYNAMICS, only: ps,mb,ma
      USE TRACER_ADV, only:
     *    AADVQ,AADVQ0,sbf,sbm,sfbm,scf,scm,sfcm,ncyc
      USE DOMAIN_DECOMP, only : grid, GET
      IMPLICIT NONE
      REAL*8 DTLF,byncyc,byma
      INTEGER I,J,L   !@var I,J,L loop variables

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)


      DTLF=2.*DT
      CALL CALC_AMP(PS,MB)
      CALL AADVQ0 (1._8)  ! uses the fluxes pua,pva,sda from DYNAM
C****
C**** convert from concentration to mass units
C****
!$OMP PARALLEL DO PRIVATE (L,J,I)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
        Q(I,J,L)=Q(I,J,L)*MB(I,J,L)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*MB(I,J,L)
      enddo; enddo; enddo
!$OMP END PARALLEL DO
C**** ADVECT
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
      CALL AADVQ (Q,QMOM, .TRUE. ,'q       ')
        byncyc = 1./ncyc
        AJL(:,:,jl_totntlh) = AJL(:,:,jl_totntlh) + sbf(:,:)
        AJL(:,:,jl_zmfntlh) = AJL(:,:,jl_zmfntlh)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        AJL(:,:,jl_totvtlh) = AJL(:,:,jl_totvtlh) + scf(:,:)
        AJL(:,:,jl_zmfvtlh)  = AJL(:,:,jl_zmfvtlh)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc
C****
C**** convert from mass to concentration units (using updated MA)
C****
!$OMP PARALLEL DO PRIVATE (L,I,J,BYMA)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
        BYMA = 1.D0/MA(I,J,L)
        Q(I,J,L)=Q(I,J,L)*BYMA
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*BYMA
      enddo; enddo; enddo
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE QDYNAM


      end module ATMDYN_QDYNAM
