#include "rundeck_opts.h"

#ifdef OBIO_ON_GARYocean
      subroutine obio_trint
#else
      subroutine obio_trint(nn)
#endif

#ifdef OBIO_ON_GARYocean
      Use OCEANRES,  only: idm=>imo,jdm=>jmo,kdm=>lmo,dZO
      USE OCEAN, only : focean
      USE OCEANR_DIM, only : ogrid
      USE obio_com, only: tracer, tracer_loc
      USE OCN_TRACER_COM, only : ntrcr=>ntm
      USE MODEL_COM, only : nstep=>itime
      USE GEOM, only : dxyp
#else
      USE hycom_dim_glob, only : jj,isp,ifp,ilp,kk,ntrcr,idm,jdm,kdm
      USE hycom_arrays_glob, only : tracer,dp,scp2
      USE hycom_scalars, only : nstep,huge
#endif

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,GLOBALSUM 


      implicit none

      integer i,j,k,l,kn,nn
      integer ntr
      real sumo,sum_area,summ
      real :: arraySum
      real :: sumtrac(idm,jdm,kdm)

      !!! need to gather tracer etc. before calling this routine
      if ( AM_I_ROOT() ) then

      do ntr=1,ntrcr
      summ=0.
      sumo=0.
        do k=1,kdm
        sum_area=0.
        do j=1,jdm
        do i=1,idm
#ifdef OBIO_ON_GARYocean
         kn=k
         if (FOCEAN(i,j).gt.0) then
           sumo=sumo+dzo(kn)*dxyp(j)     !sum volume
           sum_area=sum_area+dxyp(j)     !sum area
           summ=summ+dzo(kn)*tracer(i,j,k,ntr)*dxyp(j)   !sum over volume
         endif !focean
#else
         kn=k+nn
         if (dp(i,j,kn)<huge) then
           sumo=sumo+dp(i,j,kn)*scp2(i,j)
           sum_area=sum_area+scp2(i,j)
           summ=summ+dp(i,j,kn)*tracer(i,j,k,ntr)*scp2(i,j)
         endif
#endif
          
        enddo
        enddo
        enddo
      write(*,*)
     .  'total intgrl tracer:',ntr,nstep,summ,summ/sumo
      enddo
      print*,'   '

      endif   !if am_i_root

!     do ntr=1,ntrcr
!       write(*,*),'ntr =',ntr
!       sumtrac=tracer_loc(:,:,:,ntr)
!       write(*,*),'doing sumtrac'
!       call GLOBALSUM(ogrid,sum(sumtrac,dim=3),arraySum)
!       if(AM_I_ROOT()) write(*,*) __FILE__,__LINE__, ntr, arraySum
!     enddo

      return
      end
