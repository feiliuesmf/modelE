#include "rundeck_opts.h"
#define JJ(J) (J)-J_0H+1

      module ATMDYN
      implicit none
      private

      public init_ATMDYN, DYNAM,CALC_TROP,PGRAD_PBL
     &     ,DISSIP,FILTER,CALC_AMPK
     &     ,COMPUTE_DYNAM_AIJ_DIAGNOSTICS, COMPUTE_WSAVE
     &     ,getTotalEnergy,SDRAG   
     &     ,addEnergyAsDiffuseHeat,addEnergyAsLocalHeat     
#ifdef TRACERS_ON
     &     ,trdynam
#endif
 

      contains

      SUBROUTINE init_ATMDYN
      return
      end SUBROUTINE init_ATMDYN





      SUBROUTINE DYNAM
      USE MODEL_COM, only : im,lm,t,p,q,ls1,NSTEPSCM     
      USE SOMTQ_COM, only : tmom,mz
      USE DOMAIN_DECOMP, only : grid
      USE DYNAMICS, only : PMID,PEDN     


      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     TZ,PIJL

      INTEGER L


      call pass_SCMDATA

      CALL CALC_PIJL(LM,P,PIJL)
      CALL CALC_AMPK(LM)

      call SCM_FORCN

      CALL tq_zmom_init(T,Q,PMID,PEDN)

      DO L=1,LM
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO

      CALL PGF_SCM(T,TZ,PIJL)

      call FCONV

      return
      END SUBROUTINE DYNAM

      SUBROUTINE SCM_FORCN
c     apply advective forcings from ARM Variational analysis to T and Q

      USE RESOLUTION , only : LM
      USE MODEL_COM , only : P,T,Q,PTOP,SIG,NSTEPSCM,DTSRC        
     &                ,I_TARG,J_TARG
      USE DYNAMICS, only : PK
      USE CONSTANT , only : KAPA 

      USE SCMCOM , only : SG_HOR_TMP_ADV, SG_VER_S_ADV, SG_HOR_Q_ADV,
     &              SG_VER_Q_ADV,iu_scm_prt      
      USE CLOUDS , only : SCM_DEL_T, SCM_DEL_Q

      IMPLICIT NONE



      INTEGER L

cccccc is there some other variable they keep or function for doing this
c     write(0,*) 'enter SCM_FORCN     dtsrc ',dtsrc 
      do L = 1,LM
         T(I_TARG,J_TARG,L) = T(I_TARG,J_TARG,L)*PK(L,I_TARG,J_TARG) 
c        write(iu_scm_prt,*) 'FORCN -old tq  ',L,T(I_TARG,J_TARG,L),
c    &               Q(I_TARG,J_TARG,L)*1000.0 
      enddo
     

      do L = 1,LM
c        write(iu_scm_prt,*) 'tadvs ',L,SG_HOR_TMP_ADV(L),
c    *                       SG_VER_S_ADV(L)
         SCM_DEL_T(L) = SG_HOR_TMP_ADV(L)*DTSRC + SG_VER_S_ADV(L)*DTSRC   
         T(I_TARG,J_TARG,L) = T(I_TARG,J_TARG,L) + SCM_DEL_T(L)
c        write(iu_scm_prt,*) 'add tadv delT T ',L,SCM_DEL_T(L),
c    &               T(I_TARG,J_TARG,L)    
      enddo   
      do L = 1,LM
c        write(iu_scm_prt,*) 'qadvs ',L,SG_HOR_Q_ADV(L),SG_VER_Q_ADV(L) 
         SCM_DEL_Q(L) = SG_HOR_Q_ADV(L)*DTSRC + SG_VER_Q_ADV(L)*DTSRC    
         Q(I_TARG,J_TARG,L) = Q(I_TARG,J_TARG,L) + SCM_DEL_Q(L)
         if (Q(I_TARG,J_TARG,L).lt.0.0) then
            write(99,51) NSTEPSCM,I_TARG,J_TARG,L,Q(I_TARG,J_TARG,L)
  51        format(1x,'SCM_FORCN NSTEP  I_TARG J_TARG L Q ',
     &               4(i5),f10.7) 
            SCM_DEL_Q(L) = -Q(I_TARG,J_TARG,L)
            Q(I_TARG,J_TARG,L) = 0.0
         endif
      enddo
    
      do L = 1,LM
c        write(iu_scm_prt,*) 'FORCN - new tq  ',L,T(I_TARG,J_TARG,L),
c    &               Q(I_TARG,J_TARG,L)*1000.0 
         T(I_TARG,J_TARG,L) = T(I_TARG,J_TARG,L)/PK(L,I_TARG,J_TARG)    
      enddo




      RETURN

      END SUBROUTINE SCM_FORCN 

  
      SUBROUTINE FCONV
C*****
C     for single column model
C     compute CONV=Horizontal Mass Convergence
C     as filled in subroutine AFLUX in the GCM for use in
C     the CONDSE and MSTCNV Subroutines
C     Use the Wind Divergence from the ARM data
C     CONV = Wind Divergence*dSigma*P*DelArea
C
C     NOTE:    Wind Divergence is calculated for the area of the
C              ARM site. Therefore we need to take into account the
C              difference between the GCM grid box area and the ARM
C              Site.   Oklahoma site (SGP)  300 x 365 KM = 109500KM**2
C                      GCM 2 x 2.5 degrees (smaller for SGP)
C                          ~ 222.63 * 2223.42 = 49739.01 KM**2
C                     SGP/GCM = 2.2
C                     
c              Note: for NSA  domain for the variational analysis
c                    is  230KM (longitudinal) x 100KM (latitudinal)
c                          230x100 = 23000
c                    GCM 2x2.5 degree grid box ~ 21266
c               area/box = (sin(q1)-sin(q2))*2(pi)R**2/144
c                        72 degrees-71degrees
c                    ARMFAC = NSA/GCM = 23000/21266 ~ 1.08
c   
c              What about for TWP site ? ? ?
c
c
c
  
      USE RESOLUTION , only : LM

      USE MODEL_COM , only : P,DSIG, I_TARG, J_TARG   
      USe GEOM , only : DXYP
   
      USE SCMCOM , only : SG_WINDIV, SG_CONV    
   
      IMPLICIT NONE


      real*4 ARMFAC 
      integer L

      DATA ARMFAC/1.0/
c     DATA ARMFAC/2.2/
c     DATA ARMFAC/1.08/
      

c     want to fill SD (IDUM,JDUM)  check out 

      DO L=1,LM
         SG_CONV(L) = SG_WINDIV(L)*DSIG(L)*P(I_TARG,J_TARG)
     &                 *DXYP(J_TARG)*ARMFAC
      ENDDO

      return

      end SUBROUTINE FCONV  

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


      RETURN

      END SUBROUTINE PGRAD_PBL
      
      SUBROUTINE SDRAG(DT1)
      REAL*8, INTENT(IN) :: DT1 
      return
      END SUBROUTINE SDRAG 



      SUBROUTINE PGF_SCM (T,SZ,P)
!@SCM-version    For SCM need to calculate geopotential height. 
!                Remove other calculations.
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE MODEL_COM, only : im,jm,lm,ls1,mrch,dsig,psfmpt,sige,ptop
     *     ,zatmo,sig,modd5k,bydsig
     &     ,do_polefix,I_TARG,J_TARG   
      USE GEOM, only : imaxj,dxyv,dxv,dyv,dxyp,dyp,dxp,acor,acor2
      USE DYNAMICS, only : gz,pu,pit,phi,spa,dut,dvt
      USE DOMAIN_DECOMP, Only : grid, GET
      USE DOMAIN_DECOMP, only : HALO_UPDATE
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      USE SCMCOM, only : iu_scm_prt
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM):: T
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO):: FD,RFDUX
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  P, SZ

      REAL*8 PKE(LS1:LM+1)
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
c     do L=1,LM
c        write(iu_scm_prt,*) 'PGF_SCM  L GZ ',L,GZ(I_TARG,J_TARG,L)
c     enddo
!$OMP END PARALLEL DO
C****
C
      RETURN
      END SUBROUTINE PGF_SCM


c     SUBROUTINE AFLUX (U,V,PIJL)
c     END SUBROUTINE AFLUX


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
c     IF (haveLatitude(grid, J=1)) THEN
c       PTROPO(2:IM,1) = PTROPO(1,1)
c       LTROPO(2:IM,1) = LTROPO(1,1)
c     END IF
c     IF (haveLatitude(grid,J=JM)) THEN
c       PTROPO(2:IM,JM)= PTROPO(1,JM)
c       LTROPO(2:IM,JM)= LTROPO(1,JM)
c     END IF

      END SUBROUTINE CALC_TROP


      SUBROUTINE DISSIP
      return
      END SUBROUTINE DISSIP

      SUBROUTINE FILTER
      return
      END SUBROUTINE FILTER


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


      SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS( PUA,PVA,dt)
      USE CONSTANT,      only: BY3
      use DOMAIN_DECOMP, only: grid, get, halo_update, SOUTH
      USE DOMAIN_DECOMP, only : haveLatitude
      use DIAG_COM, only: AIJ => AIJ_loc,
     &     IJ_FGZU, IJ_FGZV, IJ_FMV, IJ_FMU
      use MODEL_COM, only: IM,JM,LM

      real*8, intent(in) :: PUA(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: PVA(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      return
      END SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS
     


      SUBROUTINE COMPUTE_WSAVE(wsave, sda, T, PK, PEDN, NIdyn)
      USE CONSTANT, only: rgas, bygrav
      !use MODEL_COM, only: NIdyn
      use DOMAIN_DECOMP, only: grid, GET
      use GEOM, only: bydxyp
      use MODEL_COM, only: IM,JM,LM

      real*8, dimension(:, grid%J_STRT_HALO:, :), intent(out) :: WSAVE
      real*8, dimension(:, grid%J_STRT_HALO:, :), intent(in)  :: SDA, T
      real*8, dimension(:, :, grid%J_STRT_HALO:), intent(in)  :: PK,PEDN
      integer, intent(in) :: NIdyn

      return
      END SUBROUTINE COMPUTE_WSAVE

      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy returns the sum of kinetic and potential energy.
!@auth Tom Clune (SIVO)
!@ver  1.0
      use GEOM, only: DXYP, AREAG
      use DOMAIN_DECOMP, only: grid, GLOBALSUM, get
      USE DOMAIN_DECOMP, only : haveLatitude
      REAL*8 :: totalEnergy

c     REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) ::KEJ,PEJ,TEJ
c     REAL*8 :: totalPotentialEnergy
c     REAL*8 :: totalKineticEnergy
c     integer :: J_0, J_1
c     logical :: HAVE_SOUTH_POLE

      totalEnergy = 0.

c     call get(grid, J_STRT=J_0, J_STOP=J_1,
c    *     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

c     call conserv_PE(PEJ)
c     call conserv_KE(KEJ)
c     if (haveLatitude(grid, J=1)) KEJ(1) = 0

c     TEJ(J_0:J_1)= (KEJ(J_0:J_1) + PEJ(J_0:J_1)*DXYP(J_0:J_1))/AREAG

c     CALL GLOBALSUM(grid, TEJ, totalEnergy, ALL=.true.)

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

c     real*8 :: ediff
c     integer :: l
c     integer :: J_0, J_1

c     call get(grid, J_STRT=J_0, J_STOP=J_1)

c     ediff = deltaEnergy / ((PSF-PMTOP)*SHA*mb2kg)
!$OMP  PARALLEL DO PRIVATE (L)
c     do l=1,lm
c       T(:,J_0:J_1,L)=T(:,J_0:J_1,L)-ediff/PK(L,:,J_0:J_1)
c     end do
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

c     integer :: i, j, k, l
c     real*8 :: ediff
c     integer :: J_0, J_1

c     call get(grid, J_STRT=J_0, J_STOP=J_1)
c     CALL HALO_UPDATE(grid, deltaKE, FROM=NORTH)
!$OMP  PARALLEL DO PRIVATE(I,J,L,ediff,K)
c     DO L=1,LM
c       DO J=J_0,J_1
c         DO I=1,IMAXJ(J)
c           ediff=0.
c           DO K=1,KMAXJ(J)     ! loop over surrounding vel points
c             ediff=ediff+deltaKE(IDIJ(K,I,J),IDJJ(K,J),L)*RAPJ(K,J)
c           END DO
c           ediff = ediff / (SHA*PK(L,I,J))
c           T(I,J,L)=T(I,J,L)-ediff
c           if (present(diagIndex)) then
c             AJL(J,L,diagIndex) = AJL(J,L,diagIndex) - ediff
c           end if
c         END DO
c       END DO
c     END DO
!$OMP  END PARALLEL DO
      end subroutine addEnergyAsLocalHeat

      end module ATMDYN


      module ATMDYN_QDYNAM
      private

      public QDYNAM

      contains

      SUBROUTINE QDYNAM
      return
      END SUBROUTINE QDYNAM


      end module ATMDYN_QDYNAM

