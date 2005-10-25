      subroutine prtmsk(mask,array,work,idm,ii,jj,offset,scale,title)
c
c --- Delete 'array' elements outside 'mask'. Then
c --- break 'array' into sections, each 'nchar' characters wide, for printing.
c
      common/linepr/lp
c
      character title*(*)
      real array(idm,1),work(idm,1)
      integer mask(idm,1)
      data nchar/76/
ccc   data nchar/80/
ccc   data nchar/132/
c
      cvmgp(a,b,c)=a*(.5+sign(.5,c))+b*(.5-sign(.5,c))
      cvmgz(a,b,ic)=cvmgp(a,b,-1.*iabs(ic))
c
      ncols=nchar/4
      do 1 n=1,jj/ncols+1
      j1=ncols*(n-1)+1
      j2=min0(ncols*n,jj)
      if (j1.gt.j2) go to 1
      write (lp,'(/'' Sec.'',i2,'' (cols'',i4,'' -'',i4,'') -- '',a)')
     .   n,j1,j2,title
ccc      if (j2.lt.j1+5) then
ccc      write (lp,'('' (Not printed. Too few columns. Save paper.)'')')
ccc      go to 1
ccc      end if
      do 2 i=1,ii
      do 3 j=j1,j2
 3    work(i,j)=cvmgz(0.,array(i,j),mask(i,j))
      do 4 j=j1,j2
 4    work(i,j)=cvmgz(0.,(work(i,j)-offset)*scale,mask(i,j))
ccc 2    write (lp,'(32i4)') (int(work(i,j)),j=j1,j2)
 2    write (lp,'(32i4)') i,(int(work(i,j)),j=j1,j2)
 1    continue
      return
      end
c
      subroutine linout(value,char,index)
c
      common/linepr/lp
c
      parameter (length=77)
      character*1 char,line(length)
      common/lnout/line
c
c --- replace n-th element of array 'line' by character 'char', where
c --- n = 'value' modulo 'length'
c --- index < 0  -- initialize 'line' by blanks before adding 'char'
c --- index > 0  -- output 'line' after adding 'char'
c
      if (index.lt.0) then
        do 1 l=1,length
 1      line(l)=' '
      end if
      if (value.gt.0.) then
        n=int(mod(value,float(length)))
        line(n)=char
      end if
      if (index.gt.0) write (lp,'(''=+='',80a1)') (line(l),l=1,length)
      return
      end
c
      subroutine pipe_init(master)
c
c --- this set of routines facilitates output comparison from two micom
c --- versions running side by side. one model, the 'slave', writes its
c --- output into a named pipe. the other model, the 'master', reads from
c --- the pipe and compares. differences are recorded in 'base.out'.
c
c --- call 'pipe_init' initially from both code versions undergoing testing.
c --- one version must set master=.true., the other must set master=.false.
c
      logical master,slave
c
      common/cmp_pipe/iunit,lpunit,slave
      iunit=39
      lpunit=38
      slave=.not.master
c
c --- open the pipe and some output files
c
      open (unit=iunit,file='cmp_pipe',status='unknown',
     .   form='unformatted')
      if (master) then
        open (unit=lpunit,file='base.out',status='unknown')
      else
        open (unit=lpunit,file='test.out',status='unknown')
      end if
c
      return
      end
c
c
      subroutine compare(field,mask,what)
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
c --- call this routine from anywhere in the code (from both versions, of
c --- course) to check whether data stored in 'field' are identical
c
      real field(idm,jdm),field1(idm,jdm)
      integer mask(idm,jdm)
      character*20 what,which
      logical slave,fail
      common/cmp_pipe/iunit,lpunit,slave
c
      if (nstep.le.450) return                !  don't start right away
c
      if (slave) then
      write (lpunit,'(2a)') 'writing for comparison: ',what
      write (lp    ,'(2a)') 'writing for comparison: ',what
      write (iunit) what,field
c
      else                                !  slave = .false.
c
      read (iunit) which,field1
      write (lpunit,'(2a)') 'reading for comparison: ',which
      write (lp    ,'(2a)') 'reading for comparison: ',which
      if (what.ne.which) then
        write (lpunit,'(4a)') 'out of sync -- trying to compare ',what,
     .     '  to  ',which
        stop
      end if
c
      fail=.false.
      do 1 j=1,jdm
      do 1 i=1,idm
      if (mask(i,j).gt.0 .and. field(i,j).ne.field1(i,j)) then
        write (lpunit,'(a,2i5,1p,3(a,e14.7))') 'i,j=',i,j,
     .    '  master:',field(i,j),'  slave:',field1(i,j),
     .    '  diff:',field(i,j)-field1(i,j)
        fail=.true.
        if (fail) return                !  optional
      end if
 1    continue
      if (fail) stop                        !  optional
c
      end if
      return
      end
c
c
      subroutine comparall(m,n,mm,nn,info)
c
c --- write out a standard menu of arrays for testing
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      common/cmp_pipe/iunit,lpunit,slave
      character text*20,info*(*)
c
      if (nstep.le.450) return                !  don't start right away
c
      write (lpunit,'(2a)') 'comparall called: ',info
      do 1 k=1,kk
      km=k+mm
      kn=k+nn
 100  format (a9,i3,a8)
      write (text,100) 'u(kn)  k=',k,info(13:20)
      call compare(u(1,1,kn),iu,text)
      write (text,100) 'v(kn)  k=',k,info(13:20)
      call compare(v(1,1,kn),iv,text)
      write (text,100) 'dp(kn) k=',k,info(13:20)
      call compare(dp(1,1,kn),ip,text)
      write (text,100) 'temp(kn) ',k,info(13:20)
      call compare(temp(1,1,kn),ip,text)
      write (text,100) 'saln(kn) ',k,info(13:20)
      call compare(saln(1,1,kn),ip,text)
      write (text,100) 'th3d(kn) ',k,info(13:20)
      call compare(th3d(1,1,kn),ip,text)
 1    continue
c
      return
      end
c
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
c
      subroutine findmx(mask,array,idm,ii,jj,name)
c
c --- find maximum and minimum in 'array'. only check points where mask > 0
c
      implicit none
      integer jmax
      parameter (jmax=2000)
      integer lp,i,j,idm,ii,jj,mask(idm,jj),ipos,jpos,ineg,jneg,
     .        jpoj(jmax),ipoj(jmax),jnej(jmax),inej(jmax)
      real array(idm,jj),difpos,difneg,difpoj(jmax),difnej(jmax),huge
      character name*(*)
      data huge/1.e33/
      common/linepr/lp
c
      if (jj.gt.jmax) stop '(error subr.findmx -- jmax < jj)'
c
c$OMP PARALLEL DO PRIVATE(difpos,difneg,ipos,jpos,ineg,jneg)
      do j=1,jj
         difpos=-huge
         difneg= huge
         do i=1,ii
           if (mask(i,j).gt.0) then
             if (array(i,j).gt.difpos) then
               difpos=array(i,j)
               ipos=i
               jpos=j
             end if
             if (array(i,j).lt.difneg) then
               difneg=array(i,j)
               ineg=i
               jneg=j
             end if
           end if
         end do
c
         difpoj(j)=difpos
         difnej(j)=difneg
         ipoj(j)=ipos
         jpoj(j)=jpos
         inej(j)=ineg
         jnej(j)=jneg
      end do
c$OMP END PARALLEL DO
c
      difpos=-huge
      difneg= huge
      ipos=-1
      jpos=-1
      ineg=-1
      jneg=-1
c
      do j=1,jj
        if (difpoj(j).gt.difpos) then
          difpos=difpoj(j)
          ipos=ipoj(j)
          jpos=jpoj(j)
        end if
        if (difnej(j).lt.difneg) then
          difneg=difnej(j)
          ineg=inej(j)
          jneg=jnej(j)
        end if
      end do
c
      write (lp,'(2a,1p,2(e11.2,2i5))') name,'  min,max =',
     .    difneg,ineg,jneg,difpos,ipos,jpos
c
      return
      end
c
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
 103  format(i3,2p,5f7.2,5x,5f7.2)                        !  SI units
c103  format(i3,   5f7.2,5x,5f7.2)                        !  cgs units
 104  format(i3,5f7.2,5x,5f7.2)
 105  format(i3,   5f7.1,5x,   5f7.2)                        !  SI units
c105  format(i3,0p,5f7.1,5x,3p,5f7.2)                        !  cgs units
 106  format(i3,5f7.1,5x,2p,5f7.1)                        !  SI units
c106  format(i3,-2p,5f7.1,5x,2p,5f7.1)                        !  cgs units
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
      do 1 ks=1,kk                                !  print out all layers
ccc   do 1 ks=max(1,k-1),min(kk,k+1)                !  print out 3 adjacent layers
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
ccc      if (1.gt.0) stop '(stencl)'                !  optional
      return
      end

