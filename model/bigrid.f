      subroutine bigrid(depth)
c
c --- set loop bounds for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
c --- this version works for both cyclic and noncyclic domains.
c --- land barrier at i=ii and/or j=jj signals closed-basin conditions
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
c
      real depth(idm,jdm)
      integer nfill,nzero,jsec,jfrst,jlast
      character fmt*12,char2*2
      data fmt/'(i4,1x,75i1)'/
c
c$OMP PARALLEL DO
      do 17 j=1,jj
      do 17 i=1,ii
      ip(i,j)=0
      iq(i,j)=0
      iu(i,j)=0
 17   iv(i,j)=0
c$OMP END PARALLEL DO
c
c --- fill single-width inlets
 16   nfill=0
c$OMP PARALLEL DO PRIVATE(ja,jb,ia,ib,nzero) REDUCTION(+:nfill)
      do 15 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 15 i=1,ii
      ia=mod(i-2+ii,ii)+1
      ib=mod(i     ,ii)+1
      nzero=0
      if (depth(i,j).gt.0.) then
        if (depth(ia,j).le.0.) nzero=nzero+1
        if (depth(ib,j).le.0.) nzero=nzero+1
        if (depth(i,ja).le.0.) nzero=nzero+1
        if (depth(i,jb).le.0.) nzero=nzero+1
        if (nzero.eq.3) then
          write (lp,'(a,i4,a,i4,a)') ' depth(',i,',',j,') set to zero'
          depth(i,j)=0.
          nfill=nfill+1
        end if
      end if
 15   continue
c$OMP END PARALLEL DO
      if (nfill.gt.0) go to 16
c
c --- mass points are defined where water depth is greater than zero
c$OMP PARALLEL DO
      do 2 j=1,jj
      do 2 i=1,ii
      if (depth(i,j).gt.0.) ip(i,j)=1
  2   continue
c$OMP END PARALLEL DO
c
c --- u,v points are located halfway between any 2 adjoining mass points
c$OMP PARALLEL DO PRIVATE(ia)
      do 3 j=1,jj
      do 3 i=1,ii
      ia=mod(i-2+ii,ii)+1
      if (ip(ia,j).gt.0.and.ip(i,j).gt.0) iu(i,j)=1
  3   continue
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(ja)
      do 4 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 4 i=1,ii
      if (ip(i,ja).gt.0.and.ip(i,j).gt.0) iv(i,j)=1
  4   continue
c$OMP END PARALLEL DO
c
c --- 'interior' q points require water on all 4 sides.
c$OMP PARALLEL DO PRIVATE(ja,ia)
      do 5 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 5 i=1,ii
      ia=mod(i-2+ii,ii)+1
      if (min(ip(i,j),ip(ia,j),ip(i,ja),ip(ia,ja)).gt.0) iq(i,j)=1
  5   continue
c$OMP END PARALLEL DO
c
c --- 'promontory' q points require water on 3 (or at least 2 diametrically 
c --- opposed) sides
c$OMP PARALLEL DO PRIVATE(ja,ia)
      do 10 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 10 i=1,ii
      ia=mod(i-2+ii,ii)+1
      if ((ip(i ,j).gt.0.and.ip(ia,ja).gt.0).or.
     .    (ip(ia,j).gt.0.and.ip(i ,ja).gt.0)) iq(i,j)=1
 10   continue
c$OMP END PARALLEL DO
c
c --- determine loop bounds for vorticity points, including interior and
c --- promontory points
      call indxi(iq,ifq,ilq,isq)
      call indxj(iq,jfq,jlq,jsq)
c
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
      call indxj(ip,jfp,jlp,jsp)
      call indxi(iu,ifu,ilu,isu)
      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
      call indxj(iv,jfv,jlv,jsv)
c
c --- write out  -ip-  array
c --- data are written in strips 75 points wide
      jsec=(jj-1)/75
      do 9 jfrst=0,75*jsec,75
      jlast=min(jj,jfrst+75)
      write (char2,'(i2)') jlast-jfrst
      fmt(8:9)=char2
      write (lp,'(''ip array, cols'',i5,'' --'',i5)') jfrst+1,jlast
 9    write (lp,fmt) (i,(10*ip(i,j),j=jfrst+1,jlast),i=1,ii)
c
      return
      end
c
c
      subroutine indxi(ipt,if,il,is)
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
c
      integer ipt(idm,jdm),if(jdm,ms),il(jdm,ms),is(jdm)
      do 1 j=1,jj
      is(j)=0
      do 4 k=1,ms
      if(j,k)=0
 4    il(j,k)=0
      i=1
      k=1
 3    if (ipt(i,j).ne.0) go to 2
      i=i+1
      if (i.le.ii) go to 3
      go to 1
 2    if (k.gt.ms) then
      write (lp,'('' error in indxi - ms too small at i,j ='',2i5)') i,j
      write (lp,'('' j-th line of ipt array:'',/(7(1x,10i1)))')
     .   (ipt(l,j),l=1,ii)
      stop '(indxi)'
      end if
      if(j,k)=i
 6    i=i+1
      if (i.le.ii) go to 5
      il(j,k)=ii
      is(j)=k
      go to 1
 5    if (ipt(i,j).ne.0) go to 6
      il(j,k)=i-1
      is(j)=k
      k=k+1
      go to 3
 1    continue
      return
      end
c
c
      subroutine indxj(jpt,jf,jl,js)
c
c --- input array jpt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays jf, jl, js  where
c --- jf(i,k) gives column index of first point in row i for k-th section
c --- jl(i,k) gives column index of last point
c --- js(i) gives number of sections in row i (maximum: ms)
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
c
      integer jpt(idm,jdm),jf(idm,ms),jl(idm,ms),js(idm)
      do 1 i=1,ii
      js(i)=0
      do 4 k=1,ms
      jf(i,k)=0
 4    jl(i,k)=0
      j=1
      k=1
 3    if (jpt(i,j).ne.0) go to 2
      j=j+1
      if (j.le.jj) go to 3
      go to 1
 2    if (k.gt.ms) then
      write (lp,'('' error in indxj - ms too small at i,j ='',2i5)') i,j
      write (lp,'('' i-th line of jpt array:'',/(7(1x,10i1)))')
     .   (jpt(i,l),l=1,jj)
      stop '(indxj)'
      end if
      jf(i,k)=j
 6    j=j+1
      if (j.le.jj) go to 5
      jl(i,k)=jj
      js(i)=k
      go to 1
 5    if (jpt(i,j).ne.0) go to 6
      jl(i,k)=j-1
      js(i)=k
      k=k+1
      go to 3
 1    continue
      return
      end
c
c
c> Revision history:
c>
c> May  2000 - routine generalized to accomodate cyclic & noncyclic b.c.
