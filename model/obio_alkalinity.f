#include "rundeck_opts.h"

      subroutine obio_alkalinity(vrbos,kmax,i,j)

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
! alk tendency is also computed under ice.  this is because 
! alk changes due to sal and temp changes.

      USE obio_dim
      USE obio_incom, only: rain_ratio,cpratio,sigma_Ca,d_Ca,
     .      npratio,uMtomgm3,cnratio
      USE obio_com, only: P_tend,p1d,pp2_1d,dp1d,A_tend,
     .      rhs,zc,alk1d,bn,pnoice

#ifdef OBIO_ON_GARYocean
      USE MODEL_COM, only: nstep=> itime
#else
      USE hycom_scalars, only: nstep
#endif

      implicit none

      integer nt,k,kmax,nchl1,nchl2,i,j
      real*8 J_PO4(kmax),pp,Jprod,Jprod_sum,Fc,zz,F_Ca(kmax+1),
     .       J_Ca(kmax),term
      logical vrbos
!--------------------------------------------------------------------------
!pp and therefore Jprod are defined at dp points (mid level)
!F_Ca is defined at pressure interfaces
!J_Ca is therefore defined mid-layer (dp points)

!compute sources/sinks of phosphate
!J_PO4 units uM/hr
      do k=1,kmax
      J_PO4(k) = 0.1 * P_tend(k,1)    !approximate by nitrate conc
                                  !NO3/PO4 ratio from Conkright et al, 1994
      term = -1.*npratio * J_PO4(k) 
      rhs(k,15,1) = term
      A_tend(k)= term
      enddo

!distinguish two cases: OCMIP uses total pp for Jprod, whereas here we also
!consider the case where Jprod only includes coccolithophores
#ifdef Jprod_based_on_pp
      nchl1=1
      nchl2=nchl
#endif
#ifdef Jprod_based_on_cocc
      nchl1=4
      nchl2=4   
#endif

      if (nchl1.le.0 .or. nchl2.le.0) then
          print*, nchl1, nchl2, 
     .   'MUST SET Jprod_based_on_pp or Jprod_based_on_cocc in rundeck'
          stop  
      endif

!compute total primary production
      pp=0.
      do nt=nchl1,nchl2
      do k=1,kmax
         if (p1d(kmax+1).gt.zc)     !total depth must be greater than 200m
     .        pp=pp+pp2_1d(k,nt)/max(p1d(k+1),1.e-3)/bn(k)/cnratio   
              
              !negative values are not accepted
!             if (pp.lt.0.) then
!                 pp=0.
!                 write(*,'(a,5i5,4e12.4)')'obio_alkalinity, pp neg:',
!    .          nstep,i,j,k,nt,pp2_1d(k,nt),max(p1d(k+1),1.e-3),
!    .          bn(k),cnratio
!             endif
      enddo
      enddo
      Jprod = pp/uMtomgm3    !convert to uM C
                             !Jprod units uM/hr


!integrate net primary production down to zc
      Jprod_sum = 0.
      do k=1,kmax
      if (p1d(k+1) .le. zc)  then
       Jprod_sum = Jprod_sum + Jprod*dp1d(k)
      endif
      enddo

!compute downward flux of CaCO3
      Fc = (1-sigma_Ca)* Jprod_sum
      do k=1,kmax+1
         if (p1d(k) .le. zc)  then
           F_Ca(k) = rain_ratio*cpratio*Fc*exp((p1d(k)-zc)/d_Ca)
         endif
      enddo

!compute sources/sinks of CaCO3
      do k=1,kmax
         if (p1d(k) .le. zc)  then
             J_Ca(k) = rain_ratio*cpratio*(1.-sigma_Ca)*Jprod
         else
             J_Ca(k) = -1* (F_Ca(k+1)-F_Ca(k)) / dp1d(k)
         endif
       term = 2.* J_Ca(k) 
       rhs(k,15,4) = term
       A_tend(k) = A_tend(k) + term      !A_term = -rN:P JPO4 + 2* JCa
      enddo

! surface boundary condition 
! we probably do not need one here since the effects of buoyancy 
! changes at the surface i.e. the E-P term are included in the 
! total mass changes of the surface model boxes. YES/NO?

!!!!!!!!!! NEED TO ADD BOTTOM BOUNDARY CONDITIONS 

      if (vrbos) then
      do k=1,kmax
      write(*,'(a,4i5,10e12.4)')'obio_alkalinity; ',
     .    nstep,i,j,k
     .   ,p1d(k),zc,J_PO4(k),pp,Jprod_sum,Fc
     .   ,F_Ca(k),J_Ca(k),alk1d(k),A_tend(k)
      enddo
      endif

      end subroutine obio_alkalinity
