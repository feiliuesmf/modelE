      subroutine psmoo(alist,blist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- alist is smoothed array, blist is work array
c --- this routine is set up to smooth data carried at -p- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
c
      real alist(idm,jdm),blist(idm,jdm),wgt
      parameter (wgt=.25)
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 1 j=1,jj
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      ja=mod(j-2+jj,jj)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (ip(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib)
      do 2 j=1,jj
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
