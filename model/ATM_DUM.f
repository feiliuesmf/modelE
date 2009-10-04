#include "rundeck_opts.h"

      SUBROUTINE conserv_AM(AM)
      USE DOMAIN_DECOMP_ATM, only : GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
c empty - not needed yet.
      am = 0d0
      return
      end subroutine


      SUBROUTINE conserv_KE(KE)
c calculate column sum of kinetic energy on the A grid, using A-grid winds.
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : lm
      USE DYNAMICS, only : ua=>ualij,va=>valij,pdsig
      USE DOMAIN_DECOMP_ATM, only : GRID,get
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: KE
      integer :: i,j,l,I_0,I_1,J_0,J_1
      call get(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      call recalc_agrid_uv
      do j=j_0,j_1
      do i=i_0,i_1
        ke(i,j) = 0d0
        do l=1,lm
          ke(i,j) = ke(i,j)+(ua(l,i,j)**2+va(l,i,j)**2)*pdsig(l,i,j)
        enddo
        ke(i,j) = .5*ke(i,j)*mb2kg
      enddo
      enddo
      return
      end subroutine conserv_KE

      SUBROUTINE calc_kea_3d(kea)
c calculates .5*(u**2+v**2) on the A grid, using A-grid winds.
      USE MODEL_COM, only : lm
      USE DOMAIN_DECOMP_ATM, only :  GRID,get
      USE DYNAMICS, only : ua=>ualij,va=>valij
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: KEA
      integer :: i,j,l,I_0,I_1,J_0,J_1
      call get(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      call recalc_agrid_uv
      do j=j_0,j_1
      do i=i_0,i_1
        do l=1,lm
          kea(i,j,l)=.5*(ua(l,i,j)**2+va(l,i,j)**2)
        enddo
      enddo
      enddo
      return
      end subroutine calc_kea_3d

      subroutine recalc_agrid_uv
!@sum recalc_agrid_uv Computes a-grid u,v from native-grid u,v
!@var u   x-component of wind, D grid
!@var v   y-component of wind, D grid
!@var ua  x-component of wind, A grid, latlon orientation, l,i,j order
!@var va  y-component of wind, A grid, latlon orientation, l,i,j order
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : lm,u,v
      USE DYNAMICS, only : ua=>ualij,va=>valij
      USE DOMAIN_DECOMP_ATM, only : grid,get
      use FV_StateMod, only : INTERP_DGRID_TO_AGRID
      implicit none
      real*8, dimension(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,
     &     lm) :: ua_tmp,va_tmp,ud_tmp,vd_tmp
      integer :: i,j,l,I_0,I_1,J_0,J_1
      call get(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          ud_tmp(i,j,l) = u(i,j,l) ! strip halo cells from u
          vd_tmp(i,j,l) = v(i,j,l) ! strip halo cells from v
        enddo
        enddo
      enddo
      call interp_dgrid_to_agrid(ud_tmp, vd_tmp, ua_tmp, va_tmp,
     &       rotate=.true.)
      do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          ua(l,i,j) = ua_tmp(i,j,l)
          va(l,i,j) = va_tmp(i,j,l)
        enddo
        enddo
      enddo
      return
      end subroutine recalc_agrid_uv

      subroutine replicate_uv_to_agrid(u_r,v_r,k,ursp,vrsp,urnp,vrnp)
c replicate the d-grid u's and v's surrounding a-grid cells into
c a-grid arrays u_r, v_r containing 2 u's and v's at each cell.
c
c           ---U_d---            --------- 
c          |    |    |          |         |
c          |  U_r(2) |          |/V_r(1)  |
c          |         |         V_d       V_d 
c          |  U_r(1) |          |  V_r(2)/|
c          |    |    |          |         |
c           ---U_d---            --------- 
c
c
      USE MODEL_COM, only : im,jm,lm,u,v ! u,v are u_d and v_d
      USE DOMAIN_DECOMP_ATM, only : get,grid,halo_update
      implicit none
      integer :: k
      real*8, dimension(k,lm,
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo) :: u_r,v_r
c ursp,vrsp,urnp,vrnp are only for compatibility with the latlon configuration
c and are not set here
      real*8, dimension(im,lm) :: ursp,vrsp,urnp,vrnp
      integer :: i,j,l,i_0,i_1,j_0,j_1
      if(k.ne.2) call stop_model('replicate_uv_to_agrid: bad k',255)
c
      call get(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
c for now, using symmetric scalar halo_update plus trivial edge adjustments.
c in future, this could be done better using a vector version, but see
c the discussion below in avg_replicated_duv_to_vgrid()
      call halo_update(grid, u)
      call halo_update(grid, v)
c at the top edge of an odd face, -u is the halo v
      if(mod(grid%tile,2).eq.1 .and. j_1.eq.jm) then
        u(i_0:i_1,jm+1,:) = -v(i_0:i_1,jm+1,:)
      endif
c at the right edge of an even face, -v is the halo u
      if(mod(grid%tile,2).eq.0 .and. i_1.eq.im) then
        v(im+1,j_0:j_1,:) = -u(im+1,j_0:j_1,:)
      endif
      do j=j_0,j_1
      do i=i_0,i_1
      do l=1,lm
        u_r(1,l,i,j) = u(i,j  ,l)
        u_r(2,l,i,j) = u(i,j+1,l)
        v_r(1,l,i,j) = v(i,j  ,l)
        v_r(2,l,i,j) = v(i+1,j,l)
      enddo
      enddo
      enddo
      return
      end subroutine replicate_uv_to_agrid


      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
c average the a-grid tendencies du,dv of replicated d-grid winds to
c d-grid points.  the .5 averaging factor is currently applied by
c the column physics, but this will soon change
c
c     ---------
c    |         |         --------- --------- 
c    |  du(1)  |        |         |         |
c    |    |    |        |         |/dv(1)   |
c     ---U_d---         |        V_d        |  
c    |    |    |        |   dv(2)/|         |
c    |  du(2)  |        |         |         |
c    |         |         --------- --------- 
c     ---------
c
      USE MODEL_COM, only : im,jm,lm,u,v  ! u,v are u_d,v_d
      USE DOMAIN_DECOMP_ATM, only : get,grid,halo_update
      implicit none
      integer :: k
      real*8, dimension(k,lm,grid%i_strt_halo:grid%i_stop_halo,
     &                       grid%j_strt_halo:grid%j_stop_halo) :: du,dv
c dusp,dvsp,dunp,dvnp are only for compatibility with the latlon configuration
c and are not used here.
      real*8, dimension(im,lm) :: dusp,dvsp,dunp,dvnp
      integer :: i,j,l,i_0,i_1,j_0,j_1
      call get(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
c for now, using symmetric scalar halo_update plus trivial edge adjustments.
c the storage pattern in du,dv prevents use of a vector halo_update here.
      call halo_update(grid, du, jdim=4)
      call halo_update(grid, dv, jdim=4)
c in the left halo of an odd face, dv is -du
      if(mod(grid%tile,2).eq.1 .and. i_0.eq.1) then
        dv(2,:,0,j_0:j_1) = -du(2,:,0,j_0:j_1)
      endif
c in the bottom halo of an even face, du is -dv
      if(mod(grid%tile,2).eq.0 .and. j_0.eq.1) then
        du(2,:,i_0:i_1,0) = -dv(2,:,i_0:i_1,0)
      endif
      do j=j_0,j_1
      do i=i_0,i_1
      do l=1,lm
        u(i,j,l) = u(i,j,l) + !.5*
     &       (du(2,l,i,j-1)+du(1,l,i,j))
        v(i,j,l) = v(i,j,l) + !.5*
     &       (dv(2,l,i-1,j)+dv(1,l,i,j))
      enddo
      enddo
      enddo
      return
      end subroutine avg_replicated_duv_to_vgrid

      subroutine regrid_atov_1d(u_a,v_a,uv1d)
c regrid an a-grid vector u_a,v_a into its d-grid components u,d
c and pack those components into a 1-dimensional array uv1d
c
c           ----------- ----------- 
c          |           |           |      U_d is stored in odd n
c          |           |           |      V_d              even
c          |           |           |
c      uv1d(n+1)   uv1d(n+3)       |
c          |           |           |  
c          |           |           |
c          |           |           |
c           --uv1d(n)-- -uv1d(n+2)- 
c
      USE DOMAIN_DECOMP_ATM, only : grid,get
      use FV_StateMod, only : INTERP_AGRID_TO_DGRID
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo)  :: u_a,v_a
      real*8, dimension((1+grid%i_stop-grid%i_strt)*
     &                  (1+grid%j_stop-grid%j_strt)*2) :: uv1d
      real*8, dimension(:,:), allocatable :: u_d,v_d
      integer :: i,j,n,I_0,I_1,J_0,J_1
      call get(grid, I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
c allocate u_d, v_d with the extra row/column for interp routine
      allocate(u_d(i_0:i_1,j_0:j_1+1),v_d(i_0:i_1+1,j_0:j_1))
      call interp_agrid_to_dgrid(
     &     u_a(i_0:i_1,j_0:j_1), v_a(i_0:i_1,j_0:j_1), u_d, v_d)
      n = 0
      do j=j_0,j_1
      do i=i_0,i_1
        n = n + 1
        uv1d(n) = u_d(i,j)
        n = n + 1
        uv1d(n) = v_d(i,j)
      enddo
      enddo
      deallocate(u_d,v_d)
      return
      end subroutine regrid_atov_1d

      subroutine get_nuv(nuv)
c reports the number of u's and v's in the domain
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: nuv
      nuv = 2*(1+grid%i_stop-grid%i_strt)*(1+grid%j_stop-grid%j_strt)
      return
      end subroutine get_nuv

      subroutine get_uv_of_n(n,uv)
c u is associated with odd n, v with even n
      USE MODEL_COM, only : u,v,lm
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: i,j
      call get_ij_of_n(n,i,j)
      if(mod(n,2).eq.1) then
        uv(1:lm) = u(i,j,1:lm)
      else
        uv(1:lm) = v(i,j,1:lm)
      endif
      return
      end subroutine get_uv_of_n

      subroutine store_uv_of_n(n,uv)
c u is associated with odd n, v with even n
      USE MODEL_COM, only : u,v,lm
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: i,j
      call get_ij_of_n(n,i,j)
      if(mod(n,2).eq.1) then
        u(i,j,1:lm) = uv(1:lm)
      else
        v(i,j,1:lm) = uv(1:lm)
      endif
      return
      end subroutine store_uv_of_n

      subroutine get_vpkey_of_n(n,vpkey)
c return an identifier of the spatial location of velocity index n.
c for the d-grid, every n is a different spatial location, as its
c velocities are not collocated.
      implicit none
      integer :: n,vpkey
      vpkey = n
      return
      end subroutine get_vpkey_of_n

      subroutine get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
c u is associated with odd n, v with even n
      implicit none
      integer :: n
      integer, dimension(2) :: ilist,jlist
      real*8, dimension(2) :: wts
      integer :: i,j
      integer :: nnbr
      nnbr = 2
      wts = .5
      call get_ij_of_n(n,i,j)
      if(mod(n,2).eq.1) then ! u_d
        ilist = i; jlist = (/ j-1, j /)
      else                   ! v_d
        ilist = (/ i-1, i /); jlist = j
      endif
      return
      end subroutine get_regrid_info_for_n

      subroutine get_ij_of_n(n,i,j)
      use model_com, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: n
      integer :: i,j
      integer :: nv,njm1,ni
      ni = (1+grid%i_stop-grid%i_strt)
      nv = 1+(n-1)/2
      njm1 = (nv-1)/ni
      j = grid%j_strt + njm1
      i = grid%i_strt + nv - njm1*ni - 1
      return
      end subroutine get_ij_of_n

      SUBROUTINE SDRAG(DT1)
!@sum  SDRAG puts a drag on the winds in the top layers of the atmosphere
!@auth Coefficients from the Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,sha
      USE MODEL_COM, only : lm,ls1,u,v,p,t,x_sdrag,csdragl,lsdrag
     &     ,ang_sdrag,Wc_Jdrag,wmax,vsdragl,psfmpt,dsig
      USE DYNAMICS, only : pk,pdsig,pedn,ualij,valij
      USE DOMAIN_DECOMP_ATM, only : grid, GET, HALO_UPDATE
      IMPLICIT NONE
!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var LSDRAG lowest level at which SDRAG_lin is applied
C**** SDRAG_const is applied above PTOP (150 mb) and below the SDRAG_lin
C**** regime (but not above P_CSDRAG)
      real*8 wl,tl,rho,cdn
      integer i,j,l
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     dpdu,dpdv,du,dv
      real*8, dimension(ls1:lm,grid%i_strt_halo:grid%i_stop_halo,
     &                         grid%j_strt_halo:grid%j_stop_halo) :: x
      real*8 xjud,xmid

      integer :: i_0, i_1, j_0, j_1
      call get(grid, i_strt = i_0, i_stop = i_1,
     &               j_strt = j_0, j_stop = j_1)

c
c calculate the slowdown factor X on the A grid
c
c      call recalc_agrid_uv ! not needed

      DO J=J_0,J_1
      DO I=I_0,I_1
        DO L=LS1,LM
          TL=T(I,J,L)*PK(L,I,J)
          RHO=PEDN(L+1,I,J)/(RGAS*TL)
          WL=SQRT(UALIJ(L,I,J)**2 +VALIJ(L,I,J)**2)
          xjud=1.
          if(Wc_JDRAG.gt.0.) xjud=(Wc_JDRAG/(Wc_JDRAG+min(WL,wmax)))**2
C**** WL is restricted to Wmax by adjusting X, if necessary;
C**** the following is equivalent to first reducing (U,V), if necessary,
C**** then finding the drag and applying it to the reduced winds
          IF(L.ge.LSDRAG) then
            CDN=(X_SDRAG(1)+X_SDRAG(2)*min(WL,wmax))*xjud
          else
            CDN=CSDRAGl(l)*xjud
          endif
          X(L,I,J)=DT1*RHO*CDN*min(WL,wmax)*GRAV*VSDRAGL(L)/PDSIG(L,I,J)
          if (wl.gt.wmax) X(L,I,J) = 1. - (1.-X(L,I,J))*wmax/wl
        END DO
      END DO
      END DO

c
c Apply the slowdown factor on the D grid
c
      call halo_update(grid,x,jdim=3)
      dpdu(:,:)=0.
      dpdv(:,:)=0.
      do l=ls1,lm
      do j=j_0,j_1
      do i=i_0,i_1
        xmid = min(.5*(x(l,i,j-1)+x(l,i,j)),1d0)
        dpdu(i,j) = dpdu(i,j) + u(i,j,l)*xmid*dsig(l)
        u(i,j,l) = u(i,j,l)*(1.-xmid)
        xmid = min(.5*(x(l,i-1,j)+x(l,i,j)),1d0)
        dpdv(i,j) = dpdv(i,j) + v(i,j,l)*xmid*dsig(l)
        v(i,j,l) = v(i,j,l)*(1.-xmid)
      enddo
      enddo
      enddo

c
c Conserve the column integrals of U and V by uniformly adding
c the lost stratospheric momentum to the tropospheric layers
c
      if(ang_sdrag.eq.1) then
        do j=j_0,j_1
        do i=i_0,i_1
          du(i,j) = dpdu(i,j)*psfmpt*2./(p(i,j)+p(i,j-1))
          dv(i,j) = dpdv(i,j)*psfmpt*2./(p(i-1,j)+p(i,j))
        enddo
        enddo
        do l=1,ls1-1
        do j=j_0,j_1
        do i=i_0,i_1
          u(i,j,l) = u(i,j,l) + du(i,j)
          v(i,j,l) = v(i,j,l) + dv(i,j)
        enddo
        enddo
        enddo
      endif

      return
      end subroutine sdrag

#ifdef CALC_GWDRAG

      SUBROUTINE GWDRAG
      USE CONSTANT, only : omega,grav,sha,kapa,rgas,radius
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, GET
      USE MODEL_COM, only : lm,zatmo,dtsrc,p,t,u,v
      USE DYNAMICS, only : ualij,valij,pk
      USE SOMTQ_COM, only : tmom,mz
      USE CLOUDS_COM,       ONLY : AIRX,LMC   
      USE STRAT,            ONLY : GWDCOL, NM, ZVARX,ZVARY,ZWT,DEFRM   
     *     ,QGWCNV,EK_globavg,dfm_type,ldef
     *     ,pbreaktop,defthresh,pconpen,ang_gwd
      USE STRAT, only : PL,PLE,TL,THL,RHO,BVF,DL,DUT,DVT,UL,VL,
     &     DP, DTIME, BYDTIME, CORIOL, AIRXS,
     &     UR, WT, CN, USRC, DUSDIF, MU_TOP, DUGWD, FLXUDM,    
     &     LDRAG, LMC0, LMC1,
     &     MU_INC, EK,
     &     RANMTN_CELL=>RANMTN,
     &     DEFRM_CELL,ZVARX_CELL,ZVARY_CELL,ZWT_CELL,ZATMO_CELL,
     &     XLIMIT
      USE GEOM, only : axyp,sinlat2d
      USE DIAG_COM, only : aij=>aij_loc
     *     ,jl_gwfirst,jl_dudtsdif,jl_sdifcoef
     *     ,ij_gw1,ij_gw2,ij_gw3,ij_gw4,ij_gw5
     *     ,ij_gw6,ij_gw7,ij_gw8,ij_gw9
      USE RANDOM
      use cs2ll_utils, only : uv_derivs_cs_agrid
      use FV_StateMod, only : INTERP_AGRID_TO_DGRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,LM) :: DUA,DVA ! no halo!
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP,   ! extra j for interp
     &                  GRID%J_STRT:GRID%J_STOP+1,LM) :: DUD
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP+1, ! extra i for interp
     &                  GRID%J_STRT:GRID%J_STOP,LM) :: DVD
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     RANMTN,ULDEF,VLDEF
      REAL*8, DIMENSION(LM) :: P00,AML
      INTEGER I,J,L,N,LTOP
      REAL*8 BVFSQ,ANGM,DPT,DUANG,XY
      integer lmax_angm(nm),lp10,lp2040,lpshr
C
      INTEGER :: I_0,I_1, J_0,J_1
      integer :: nij_before_j0,nij_after_j1,nij_after_i1 ! funcs
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, I_STRT=I_0,I_STOP=I_1, J_STRT=J_0,J_STOP=J_1)

      BYDTIME = 1./DTsrc
      DTIME = DTsrc
      EK(:) = EK_globavg ! for now, no horizontal variation
      xlimit = .1d0

c
c recalculate A-grid winds
c
      call recalc_agrid_uv

C****
C**** DEFORMATION
C****
      do j=j_0,j_1
        do i=i_0,i_1
          uldef(i,j) = ualij(ldef,i,j)
          vldef(i,j) = valij(ldef,i,j)
        enddo
      enddo
      call uv_derivs_cs_agrid(grid,dfm_type,uldef,vldef,deform=defrm)
      defrm = defrm/radius ! move division by radius to cs2ll_utils

c
c Set Random numbers for Mountain Drag intermittency
c
      CALL BURN_RANDOM(nij_before_j0(J_0))
      DO J=J_0,J_1
        CALL BURN_RANDOM((I_0-1))
        DO I=I_0,I_1
          RANMTN(I,J) = RANDU(XY)
        ENDDO
        CALL BURN_RANDOM(nij_after_i1(I_1))
      ENDDO
      CALL BURN_RANDOM((nij_after_j1(J_1)))

      DO L=1,LM
        DUA(:,:,L)=0.
        DVA(:,:,L)=0.
      ENDDO

C****
C**** LOOP OVER I,J
C****
      DO J=J_0,J_1
      DO I=I_0,I_1

      CN(:)=0.
      CORIOL=2.*omega*abs(sinlat2d(i,j))

C****
C**** CALCULATE 1D ARRAYS
C****
      CALL CALC_VERT_AMP(P(I,J),LM,P00,AML,DP,PLE,PL)
      DO L=1,LM
        TL(L)=PK(L,I,J)*T(I,J,L)
        THL(L)=T(I,J,L)
        RHO(L)=PL(L)/(RGAS*TL(L))
        BVFSQ=2.*TMOM(MZ,I,J,L)/(DP(L)*THL(L))*GRAV*GRAV*RHO(L)
        IF (PL(L).GE..4d0) THEN
          BVF(L)=SQRT(MAX(BVFSQ,1.d-10))
        ELSE
          BVF(L)=SQRT(MAX(BVFSQ,1.d-4))
        ENDIF
        DL(L)=0.
        DUT(L)=0.
        DVT(L)=0.
        UL(L)=UALIJ(L,I,J)
        VL(L)=VALIJ(L,I,J)
      ENDDO
C**** Levels for angular momentum restoration
      do l=lm,1,-1
        if (pl(l).lt.500.) lp10 = l
        if (pl(l).lt.400.) lp2040 = l
        if (pl(l).lt.200.) lpshr = l
      enddo
      if (pl(lpshr)-ple(lpshr-1).lt.pl(lpshr)-200.) lpshr=lpshr-1
      if (pl(lp2040)-ple(lp2040-1).lt.pl(lp2040)-400.) lp2040=lp2040-1
      if (pl(lp10)-ple(lp10-1).lt.pl(lp10)-500.) lp10=lp10-1
      lmax_angm(1) = lpshr  !Mountain
      lmax_angm(9) = lpshr  !Deformation
      lmax_angm(2) = lpshr  !Shear
      lmax_angm(3) = lp10
      lmax_angm(4) = lp10
      lmax_angm(5) = lp2040
      lmax_angm(6) = lp2040
      lmax_angm(7) = lp2040
      lmax_angm(8) = lp2040
C 
      RANMTN_CELL = RANMTN(I,J) 
      ZVARX_CELL  = ZVARX(I,J) 
      ZVARY_CELL  = ZVARY(I,J) 
      ZWT_CELL    = ZWT(I,J) 
      ZATMO_CELL  = ZATMO(I,J) 
      DEFRM_CELL  = DEFRM(I,J) 
C 
      AIRXS = 0.0 
      LMC0  = 0 
      LMC1  = 0 
      IF(QGWCNV.EQ.1) THEN
        AIRXS=AIRX(I,J)/AXYP(I,J)
        IF (AIRXS.GT.0.)  THEN
          LMC0=LMC(1,I,J)
          LMC1=LMC(2,I,J)
        ENDIF 
      ENDIF 

c
c     Call the column GWD routine
c 
      CALL GWDCOL

      IF (LDRAG.LE.LM) THEN

C**** AM conservation
      if (ang_gwd.gt.0) then    ! add in ang mom
        DO N=1,NM
          ANGM = 0.
          DO L=LDRAG-1,LM
            ANGM = ANGM - DUGWD(L,N)*DP(L)
          ENDDO
          DPT=0
          DO L=1,LMAX_ANGM(N)
            DPT=DPT+DP(L)
          ENDDO
          DUANG = ANGM/DPT
          DO L=1,LMAX_ANGM(N)
             DUT(L) = DUT(L) + DUANG
          ENDDO
        ENDDO
      endif

C****
C**** Store the U and V increments
C****
      DO L=1,LM ! LDRAG-1,LM?
        DUA(I,J,L)=DUT(L)
        DVA(I,J,L)=DVT(L)
      ENDDO

c
c     UPDATE DIAGNOSTICS
c
      AIJ(I,J,IJ_GW1)=AIJ(I,J,IJ_GW1)+MU_INC(9)*UR(9)
      AIJ(I,J,IJ_GW2)=AIJ(I,J,IJ_GW2)+MU_INC(1)*UR(1)  *WT(1)
      AIJ(I,J,IJ_GW3)=AIJ(I,J,IJ_GW3)+MU_INC(2)*UR(2)
      AIJ(I,J,IJ_GW4)=AIJ(I,J,IJ_GW4)+MU_INC(3)*UR(3)  *WT(3)
      AIJ(I,J,IJ_GW5)=AIJ(I,J,IJ_GW5)+MU_INC(7)*UR(7)  *WT(7)
      AIJ(I,J,IJ_GW6)=AIJ(I,J,IJ_GW6)+MU_INC(5)*UR(5)  *WT(5)
      AIJ(I,J,IJ_GW7)=AIJ(I,J,IJ_GW7)+CN(2)*UR(2)
      AIJ(I,J,IJ_GW8)=AIJ(I,J,IJ_GW8)+USRC
      DO N=1,NM 
        AIJ(I,J,IJ_GW9)=AIJ(I,J,IJ_GW9)+MU_TOP(N) 
      ENDDO 

      DO N=1,NM 
        DO L=1,LM 
          call inc_ajl(i,j,l,N+JL_gwFirst-1,DUGWD(L,N))
        ENDDO 
      ENDDO 

      DO L=LDRAG-1,LM
        call inc_ajl(i,j,l-1,JL_DUDTSDIF,DUSDIF(L))
        call inc_ajl(i,j,l-1,JL_SDIFCOEF,DL(L)/(BVF(L)*BVF(L)))
      ENDDO

      ENDIF ! LDRAG.LE.LM
      ENDDO ! END OF LOOP OVER I
      ENDDO ! END OF LOOP OVER J

c
c Interpolate wind increments to the D grid
c
      call interp_agrid_to_dgrid(dua,dva, dud,dvd)

c
c Update D-grid winds
c
      do l=1,lm
        do j=j_0,j_1
          do i=i_0,i_1
            u(i,j,l) = u(i,j,l) + dud(i,j,l)
            v(i,j,l) = v(i,j,l) + dvd(i,j,l)
          enddo
        enddo
      enddo

C**** conservation diagnostic
c      CALL DIAGCD (GRID,6,UT,VT,DUT3,DVT3,DT1)

      RETURN
      END SUBROUTINE GWDRAG

#else

c      SUBROUTINE DUMMY_STRAT
c!@sum DUMMY dummy routines for non-stratospheric models
cC**** Dummy routines in place of STRATDYN
c      ENTRY init_GWDRAG
c      ENTRY GWDRAG
c      ENTRY VDIFF
c      ENTRY alloc_strat_com
c      RETURN
c      END SUBROUTINE DUMMY_STRAT

#endif
