      subroutine stencl(iz,jz,k,mn)
c
c --- write 5 x 5 point cluster of grid point values centered on (iz,jz)
c --- input parameters: k = layer index; mn = time slot index, i.e., mm or nn
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      integer mn,ks,iz,jz
c
 99   format(13x,a10,30x,a10/i9,4i7,i12,4i7)
 100  format(13x,a7,i3,30x,a7,i3/i9,4i7,i12,4i7)
 101  format(i3,1p,5e7.0,5x,5e7.0)
 102  format(i3,5f7.1,5x,5f7.1)
 103  format(i3,2p,5f7.2,5x,5f7.2)			!  SI units
c103  format(i3,   5f7.2,5x,5f7.2)			!  cgs units
 104  format(i3,5f7.2,5x,5f7.2)
 105  format(i3,   5f7.1,5x,   5f7.2)			!  SI units
c105  format(i3,0p,5f7.1,5x,3p,5f7.2)			!  cgs units
 106  format(i3,5f7.1,5x,2p,5f7.1)			!  SI units
c106  format(i3,-2p,5f7.1,5x,2p,5f7.1)			!  cgs units
c
      write (lp,99) 'ice cover',' ice cover',(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,106)
     .  (i,(oice(i,j),j=jz-2,jz+2),
     .     (oice(i,j),j=jz-2,jz+2),i=iz-2,iz+2)
c
      write (lp,99) '  ubavg   ','   vbavg  ',(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,103)
     .  (i,(ubavg(i,j,1),j=jz-2,jz+2),
     .     (vbavg(i,j,1),j=jz-2,jz+2),i=iz-2,iz+2)
c
      do 1 ks=1,kk				!  print out all layers
ccc   do 1 ks=max(1,k-1),min(kk,k+1)		!  print out 3 adjacent layers
c
      write (lp,100) 'u at k=',ks,'v at k=',ks,(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,103)
     .  (i,(u(i,j,ks+mn),j=jz-2,jz+2),
     .     (v(i,j,ks+mn),j=jz-2,jz+2),i=iz-2,iz+2)
c
      write (lp,100) 'uflx k=',ks,'vflx k=',ks,(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,101)
     .  (i,(uflx(i,j,ks),j=jz-2,jz+2),
     .     (vflx(i,j,ks),j=jz-2,jz+2),i=iz-2,iz+2)
c
      write (lp,100) 'temp k=',ks,'saln k=',ks,(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,104)
     .  (i,(temp(i,j,ks+mn),j=jz-2,jz+2),
     .     (saln(i,j,ks+mn),j=jz-2,jz+2),i=iz-2,iz+2)
c
      write (lp,100) 'dp_o k=',ks,'dp_n k=',ks,(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,102)
     .  (i,(dpold(i,j,ks)/onem,j=jz-2,jz+2),
     .     (dp(i,j,ks+mn)/onem,j=jz-2,jz+2),i=iz-2,iz+2)
c
      write (lp,100) 'pres k=',ks+1,'th3d k=',ks,(j,j=jz-2,jz+2)
     .   ,(j,j=jz-2,jz+2)
      write (lp,105)
     .  (i,(p(i,j,ks+1)/onem,j=jz-2,jz+2),
     .  (th3d(i,j,ks)+thbase,j=jz-2,jz+2),i=iz-2,iz+2)
c
 1    continue
ccc      if (1.gt.0) stop '(stencl)'		!  optional
      return
      end

