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
      USE obio_com, only: tracer       !use global array
      USE OCN_TRACER_COM, only : ntrcr=>ntm
      USE MODEL_COM, only : nstep=>itime
      USE GEOM, only : dxyp
#else
      USE hycom_dim_glob, only : jj,isp,ifp,ilp,kk,ntrcr,idm,jdm,kdm
      USE hycom_arrays_glob, only : tracer,dpinit,scp2
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
         if (dpinit(i,j,k)<huge) then
           sumo=sumo+dpinit(i,j,k)*scp2(i,j)
           sum_area=sum_area+scp2(i,j)
           summ=summ+dpinit(i,j,k)*tracer(i,j,k,ntr)*scp2(i,j)
          if (nstep.eq.97.and.ntr.eq.1)
     .     write(*,*)'obio_trint: ',
     .         nstep,i,j,k,ntr,
     .         dpinit(i,j,k),scp2(i,j),tracer(i,j,k,ntr)
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

      return
      end
