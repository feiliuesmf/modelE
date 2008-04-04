      subroutine obio_limits(name)

! compute min max for each tracer and check if in bounds

!!#ifdef TRACERS_OceanBiology
      USE obio_dim
!!#endif


      USE hycom_dim
      USE hycom_arrays
      USE hycom_scalars, only : lp
      implicit none
 
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "common_blocks.h"

      data trminima/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.2,1913./
      data trmaxima/40.,0.6,108.,3.5,3.,0.6,0.4,0.3,0.8,16.3,2400./

      integer nt,ineg1,jneg1,ipos1,jpos1
      integer kneg1,kpos1
      real tracer_min,tracer_max
      real trminima(11)
      real trmaxima(11)
      character name*(*)


      do nt=1,ntrcr
        tracer_min= 1.e30
         tracer_max=-1.e30
          ineg1=-1
          jneg1=-1
          kneg1=-1
          ipos1=-1
          ipos1=-1
          kpos1=-1
          do k=1,kdm
           do 1000 j=1,jj
           do 1000 l=1,isp(j)
           do 1000 i=ifp(j,l),ilp(j,l)
            if (ip(i,j).gt.0 .and. tracer(i,j,k,nt).le.tracer_min) then
                tracer_min = tracer(i,j,k,nt)
                 ineg1=i
                 jneg1=j
                 kneg1=k
            endif
            if (ip(i,j).gt.0 .and. tracer(i,j,k,nt).ge.tracer_max) then
                tracer_max = tracer(i,j,k,nt)
                 ipos1=i
                 jpos1=j
                 kpos1=k
            endif
 1000      continue
          enddo
       write(lp,'(a,a,i3,a,2(es9.2,1x,3i4))')
     . name,', tracer=',nt,', min,max=',
     . tracer_min,ineg1,jneg1,kneg1,
     . tracer_max,ipos1,jpos1,kpos1
      enddo

      return
      end
