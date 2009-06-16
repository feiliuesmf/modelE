#include "rundeck_opts.h"
#define JJ(J) (J)-J_0H+1
!
!dummy routines for running SCM   --   A. Wolf
!
      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy Dummy
!@auth Tom Clune (SIVO)
!@ver  1.0
!    
      use GEOM, only: DXYP, AREAG
      use DOMAIN_DECOMP_1D, only: grid, GLOBALSUM, get
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      REAL*8 :: totalEnergy

      totalEnergy = 0.

      return
      end function getTotalEnergy


      subroutine addEnergyAsDiffuseHeat(deltaEnergy)
!@sum  addEnergyAsDiffuseHeat Dummy
!@auth Tom Clune (SIVO)
!@ver  1.0
      real*8, intent(in) :: deltaEnergy

      return
      end subroutine addEnergyAsDiffuseHeat
     
      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE (m^2/s^2) as heat locally
!@auth Gavin Schmidt
      USE MODEL_COM, only : t
      USE DYNAMICS, only : dke,kea,pk
      IMPLICIT NONE
     
      return

      END SUBROUTINE DISSIP

C***** Add in dissipiated KE as heat locally
      subroutine addEnergyAsLocalHeat(deltaKE, T, PK, diagIndex)
!@sum  addEnergyAsLocalHeat Dummy
!@auth Tom Clune (SIVO)
!@ver  1.0
      use DOMAIN_DECOMP_1D, only: grid, get
      implicit none
      real*8 :: deltaKE(:,grid%j_strt_halo:,:)
      real*8 :: T(:,grid%j_strt_halo:,:)
      real*8 :: PK(:,:,grid%j_strt_halo:)
      integer, optional, intent(in) :: diagIndex

      return
      end subroutine addEnergyAsLocalHeat

C**** Calculate 3D vertical velocity (take SDA which has units
C**** mb*m2/s (but needs averaging over no. of leap frog timesteps)
C**** and convert to WSAVE, units of m/s):

      subroutine COMPUTE_WSAVE(wsave, sda, T, PK, PEDN)
!@sum COMPUTE_WSAVE Dummy
      use DOMAIN_DECOMP_1D, only: grid, GET
      use MODEL_COM, only: IM,JM,LM
      implicit none

      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO, lm),
     &     intent(out) :: WSAVE
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO, lm),
     &     intent(in)  :: SDA, T
      real*8, dimension(lm, grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     intent(in)  :: PK,PEDN
      integer, intent(in) :: NIdyn

      integer :: i, j, l
      integer :: I_0, I_1, J_0, J_1

      call get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      end subroutine COMPUTE_WSAVE
