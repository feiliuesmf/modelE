      subroutine newbot
c
c --- hycom version 0.9
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      real factor
c
c --- call this routine if bottom topography has been altered in any way.
c --- note: removing grid points entirely (by setting their depth to
c --- zero) should be done elsewhere, i.e., before bigrid is called.
c
ccc      if (idm.eq.181 .and. jdm.eq.180)               !  global 2 deg grid
ccc     .   pbot(166,179)=150.*onem
c
c$OMP PARALLEL DO PRIVATE(factor)
      do 3 j=1,jj
      do 3 l=1,isp(j)
c
      do 1 k=1,kk
      do 1 i=ifp(j,l),ilp(j,l)
 1    p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
c
      do 2 i=ifp(j,l),ilp(j,l)
      pbot(i,j)=max(pbot(i,j),30.*onem)    ! 30m minimum
      if (abs(p(i,j,kk+1)-pbot(i,j)).gt.onem)
     .  write (lp,'(2i5,a,2f9.1)') i,j,
     .   '  old/new bottom depth:',p(i,j,kk+1)/onem,pbot(i,j)/onem
 2    continue
c
      do 3 k=1,kk
      do 3 i=ifp(j,l),ilp(j,l)
      factor=pbot(i,j)/p(i,j,kk+1)
      if (factor.lt..99999 .or. factor.gt.1.00001) then
        dp(i,j,k   )=dp(i,j,k   )*factor
        dp(i,j,k+kk)=dp(i,j,k+kk)*factor
c
        if (k.gt.1) then
          psikk(i,j)=psikk(i,j)-p(i,j,k)*(factor-1.)*
     .     (thstar(i,j,k)-thstar(i,j,k-1))*thref
        end if
      end if
 3    continue
c$OMP END PARALLEL DO
c
      return
      end
