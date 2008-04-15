      subroutine obio_trint

      USE hycom_dim_glob
      USE hycom_arrays_glob
      USE hycom_scalars, only : lp, nstep
      implicit none
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "common_blocks.h"

      integer ntr
      real sum,sum_area,summ

      print*,'   '
      do ntr=1,ntrcr
      summ=0.
        do k=1,kk
        sum=0.
        sum_area=0.
        do j=1,jj
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
           kn=k+nn
           sum=sum+dp(i,j,kn)*tracer(i,j,k,ntr)*scp2(i,j)
           sum_area=sum_area+dp(i,j,kn)*scp2(i,j)

           summ=summ+dp(i,j,kn)*tracer(i,j,k,ntr)*scp2(i,j)
        enddo
        enddo
        enddo
cdiag write(lp,'(a,i3,i5,a,i3,2e12.4)')
cdiag.  'total intgrl and mean for tracer ',ntr,nstep,', k= ',
cdiag.   k,sum,sum/sum_area
        enddo
!     write(lp,'(a,i3,i5,e12.4)')
      write(lp,*)
     .  'total intgrl for tracer ',ntr,nstep,summ
      enddo

      return
      end
