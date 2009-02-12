#include "rundeck_opts.h"

      subroutine obio_alkalinity(kmax)

!@sum  online computation of alkalinity
!@auth Natassa Romanou

! compute source/sink term for alkalinity
! zc     compensation depth: photosynthesis = respiration
! J_PO4  source/sinks of phosphate
! J_Ca   source/sinks of calcium carbonate
! Jprod  total productivity. Really, production of organic phosphorus 
!        above the compensation depth
! F_Ca   downward flux of CaCO3
! Fc     instantaneous downward flux of particulate organic phosphorus at 
!        compensation depth
! alkalinity units should be umol/kg

      USE obio_dim
      USE obio_incom, only: rain_ratio,cpratio,sigma_Ca,d_Ca,zc,npratio
      USE obio_com, only: P_tend,p1d,pp2_1d,dp1d,A_tend

      implicit none

      integer nt,k,kmax
      real*8 J_PO4(kmax),pp,Jprod,Jprod_sum,Fc,zz,F_Ca(kmax),J_Ca
!--------------------------------------------------------------------------


!integrate net primary production down to zc
      do nt=1,nchl
      do k=1,kmax
         if (p1d(kmax+1).gt.zc)     !total depth must be greater than 200m
     .        pp=pp+pp2_1d(k,nt)
      enddo
      enddo
      Jprod = pp

      J_PO4 = 0.1 * P_tend(:,1)    !approximate by nitrate conc
                                   !NO3/PO4 ratio from Conkright et al, 1994

      
      Jprod_sum = 0.
      do k=1,kmax
      if (p1d(k+1) .le. zc)  then
       Jprod_sum = Jprod_sum + Jprod*dp1d(k)
      endif
      enddo

      Fc = (1-sigma_Ca)* Jprod_sum
      do k=1,kmax
         if (p1d(k+1) .le. zc)  then
           F_Ca(k) = rain_ratio*cpratio*Fc*exp((p1d(k+1)-zc)/d_Ca)
         endif
      enddo

      do k=1,kmax

      zz = (p1d(k+1)-p1d(k))/2.

      if (zz .le. zc)  then

        J_Ca = rain_ratio*cpratio*(1.-sigma_Ca)*Jprod

      else

        J_Ca = -1* (F_Ca(k-1)-F_Ca(k))/ (p1d(k)-zc)

      endif

      A_tend(k) = -1.*npratio * J_PO4(k) + 2* J_Ca

      enddo

!!!!!!!!!! NEED TO ADD BOTTOM BOUNDARY CONDITIONS AND SURFACE FORCING FOR ALKALINITY

      end subroutine obio_alkalinity
