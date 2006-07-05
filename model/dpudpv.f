      subroutine dpudpv(mmnn)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
      implicit none
c
c --- ----------------------------------
c --- define layer depth at  u,v  points
c --- ----------------------------------
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      integer mmnn,kmn
c
      do k=1,kk
        call cpy_p(p(1,1,k+1))
      end do
c
c$OMP PARALLEL DO PRIVATE(ja,kmn)
      do j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 155 k=1,kk
      kmn=k+mmnn
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
 154  dpu(i,j,kmn)=max(0.,
     .           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     .           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
 155  dpv(i,j,kmn)=max(0.,
     .           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,ja ,k+1)))-
     .           min(depthv(i,j),.5*(p(i,j,k  )+p(i,ja ,k  ))))
c
      end do
c$OMP END PARALLEL DO
      return
      end
c
c
c> Revision history:
c>
c> Sep. 2000 - (dpu,dpv,depthu,depthv,p) no longer passed as arguments
c> Mar. 2006 - added bering strait exchange logic
