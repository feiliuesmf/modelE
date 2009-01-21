      module ATMDYN
      implicit none
      end module atmdyn


      SUBROUTINE conserv_AM(AM)
      USE MODEL_COM, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
      end subroutine


      SUBROUTINE conserv_KE(RKE)
      USE MODEL_COM, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RKE
      end subroutine

      SUBROUTINE calc_kea_3d(kea)
      USE MODEL_COM, only : im,lm
      USE DOMAIN_DECOMP_ATM, only :  GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: KEA
      end subroutine


      subroutine recalc_agrid_uv
      end subroutine


      subroutine replicate_uv_to_agrid(ur,vr,k,ursp,vrsp,urnp,vrnp)
      USE MODEL_COM, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     UR,VR
      real*8, dimension(im,lm) :: ursp,vrsp,urnp,vrnp
      end subroutine


      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
      USE MODEL_COM, only : im,jm,lm,u,v
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     DU,DV
      real*8, dimension(im,lm) :: dusp,dvsp,dunp,dvnp
      end subroutine


      subroutine regrid_atov_1d(u_a,v_a,uv1d)
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP_ATM, only : grid
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo)  ::
     &          u_a,v_a
      real*8, dimension(2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)),
     &        intent(out) :: uv1d
      end subroutine

      subroutine get_nuv(nuv)
      integer :: nuv
      end subroutine


      subroutine get_uv_of_n(n,uv)
      USE MODEL_COM, only : lm
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      end subroutine

      
      subroutine get_vpkey_of_n(n,vpkey)
      implicit none
      integer :: n,vpkey
      end subroutine


      subroutine get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
      implicit none
      integer :: n
      integer, dimension(4) :: ilist,jlist
      real*8, dimension(4) :: wts
      integer :: nnbr
      end subroutine


      subroutine store_uv_of_n(n,uv)
      USE MODEL_COM, only : lm
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      end subroutine

