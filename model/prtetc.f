#include "hycom_mpi_hacks.h"
      subroutine prtmsk(mask,array,work,idm,ii,jj,offset,scale,title)
c
c --- Delete 'array' elements outside 'mask'. Then
c --- break 'array' into sections, each 'nchar' characters wide, for printing.
c
      implicit none
      integer lp
      common/linepr/lp
c
      integer idm,ii,jj
      character title*(*)
      real array(idm,1),work(idm,1),offset,scale
      integer mask(idm,1)
      !--- local vars
      real cvmgp,cvmgz,a,b,c
      integer nchar,ncols,n,j1,j2,i,j,ic
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
      implicit none
      real value
      integer index
      integer lp
      common/linepr/lp
c
      integer, parameter :: length=77
      character*1 char,line(length)
      integer l,n
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
      implicit none
      logical master,slave
c
      integer iunit,lpunit
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
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS
      USE HYCOM_ARRAYS_GLOB
      implicit none
!!#include "dimensions.h"
!!#include "dimension2.h"    ! TNL
      integer i,j
!!#include "common_blocks.h"
c
c --- call this routine from anywhere in the code (from both versions, of
c --- course) to check whether data stored in 'field' are identical
c
      real field(idm,jdm),field1(idm,jdm)
      integer mask(idm,jdm)
      character*20 what,which
      logical slave,fail
      integer iunit,lpunit
      common/cmp_pipe/iunit,lpunit,slave
c
      if (nstep.le.24) return                ! don't start right away
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
      subroutine compare_r4(field,mask,what)
      USE HYCOM_DIM_GLOB, only : idm,jdm
c
!!#include "dimensions.h"
!!#include "dimension2.h"    ! TNL
      integer i,j
!!#include "common_blocks.h"
c
c --- call this routine from anywhere in the code (from both versions, of
c --- course) to check whether data stored in 'field' are identical
c
      real*4 field(idm,jdm),field1(idm,jdm)
      integer mask(idm,jdm)
      character*20 what,which
      logical slave,fail
      common/cmp_pipe/iunit,lpunit,slave
c
      if (nstep.le.24) return                ! don't start right away
c
      if (slave) then
      write (lpunit,'(2a)') 'writing r4 for comparison: ',what
      write (lp    ,'(2a)') 'writing r4 for comparison: ',what
      write (iunit) what,field
c
      else                                !  slave = .false.
c
      read (iunit) which,field1
      write (lpunit,'(2a)') 'reading r4 for comparison: ',which
      write (lp    ,'(2a)') 'reading r4 for comparison: ',which
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
      subroutine comparall(m,n,mm,nn,info)
c
c --- write out a standard menu of arrays for testing
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : nstep
      USE HYCOM_ARRAYS_GLOB
      implicit none
!!#include "dimensions.h"
!!#include "dimension2.h"    ! TNL
      integer k,m,n,mm,nn,km,kn
!!#include "common_blocks.h"
c
      integer iunit,lpunit
      logical slave
      common/cmp_pipe/iunit,lpunit,slave
      character text*20,info*(*)
c
      if (nstep.le.24) return                !  don't start right away
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
      subroutine findmx(mask,array,idm,ii,jj,name)
c
c --- find maximum and minimum in 'array'. only check points where mask > 0
c
      USE HYCOM_DIM_GLOB, only : jchunk
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
!!! hack : not defined ?
      ipos = 0 ; jpos = 0
      ineg = 0 ;jneg = 0

      if (jj.gt.jmax) stop '(error subr.findmx -- jmax < jj)'
c
c$OMP PARALLEL DO PRIVATE(difpos,difneg,ipos,jpos,ineg,jneg)
c$OMP+ SCHEDULE(STATIC,jchunk)
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
      USE HYCOM_DIM_GLOB, only : kk
      USE HYCOM_SCALARS, only : lp,onem
      USE HYCOM_ARRAYS_GLOB
      implicit none
!!#include "dimensions.h"
!!#include "dimension2.h"    ! TNL
      integer i,j,k
!!#include "common_blocks.h"
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
     .  (th3d(i,j,ks),j=jz-2,jz+2),i=iz-2,iz+2)
c
 1    continue
ccc      if (1.gt.0) stop '(stencl)'                !  optional
      return
      end
c
      subroutine prt9x9(array,iz,jz,offset,scale,what)
c
c --- write 9 x 9 point cluster of 'array' values centered on (iz,jz).
c --- the printed numbers actually represent (array(i,j) + offset) * scale
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : lp
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j
c
      real array(idm,jdm),scale,offset
      character what*12
      integer iz,jz,jwrap
      jwrap(j)=mod(j-1+jj,jj)+1		!  for use in cyclic domain
c
 100  format(a10,2x,9i7)
 101  format(i10,3x,9f7.1)
c
      write (lp,100) what,(jwrap(j),j=jz-4,jz+4)
      write (lp,101) (i,(scale*(array(i,jwrap(j))+offset),
     .   j=jz-4,jz+4),i=iz-4,iz+4)
c
      return
      end
c
c
      subroutine pr_9x9(array,idm,jdm,iz,jz,offset,scale,what)
c
c --- (sames as prt9x9 except that array dimensions are passed as arguments)
c
c --- write 9 x 9 point cluster of 'array' values centered on (iz,jz).
c --- the printed numbers actually represent (array(i,j) + offset) * scale
c
      USE HYCOM_SCALARS, only : lp
      implicit none
c
      real array(idm,jdm),scale,offset
      character what*12
      integer lp,idm,jdm,i,j,iz,jz,jwrap
      common /linepr/ lp
      jwrap(j)=mod(j-1+jdm,jdm)+1               !  for use in cyclic domain
c
 100  format(a12,9i7)
 101  format(i10,3x,9f7.1)
c
      write (lp,100) what,(jwrap(j),j=jz-4,jz+4)
      write (lp,101) (i,(scale*(array(i,jwrap(j))+offset),
     .   j=jz-4,jz+4),i=iz-4,iz+4)
c
      return
      end
c
c
      subroutine prt7x7(array,iz,jz,what)
c
c --- write 7 x 7 point cluster of 'array' values centered on (iz,jz).
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : lp
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j
c
      real array(idm,jdm),scale,offset
      character what*12
      integer iz,jz,jwrap
      jwrap(j)=mod(j-1+jj,jj)+1		!  for use in cyclic domain
c
 100  format(a10,2x,7i9)
 101  format(i10,3x,1p,7e9.1)
c
      write (lp,100) what,(jwrap(j),j=jz-3,jz+3)
      write (lp,101) (i,(array(i,jwrap(j)),
     .   j=jz-3,jz+3),i=iz-3,iz+3)
c
      return
      end
c
      subroutine psmo1(alist,pbot)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c
c --- psmo1 is specially set up for interface smoothing.
c --- it only alters -alist- values that don't coincide with -pbot-.
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real,intent(INOUT) :: alist(idm,jdm)
      real,intent(IN)    :: pbot(idm,jdm)
      real blist(idm,jdm),flxlo,flxhi
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ia,flxhi,flxlo) SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      do 5 i=1,ii
5     blist(i,j)=0.
      do 1 l=1,isu(j)
      do 1 i=ifu(j,l),ilu(j,l)
      ia=mod(i-2+ii,ii)+1
      flxhi= .25*(pbot(i ,j)-alist(i ,j))
      flxlo=-.25*(pbot(ia,j)-alist(ia,j))
 1    blist(i,j)=min(flxhi,max(flxlo,wgt*(alist(ia,j)-alist(i,j))))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ib) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ib=mod(i,ii)+1
 2    alist(i,j)=alist(i,j)-(blist(ib,j)-blist(i,j))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja,flxhi,flxlo) SCHEDULE(STATIC,jchunk)
      do 3 j=1,jj
      do 6 i=1,ii
 6    blist(i,j)=0.
      ja=mod(j-2+jj,jj)+1
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
      flxhi= .25*(pbot(i,j )-alist(i,j ))
      flxlo=-.25*(pbot(i,ja)-alist(i,ja))
 3    blist(i,j)=min(flxhi,max(flxlo,wgt*(alist(i,ja)-alist(i,j))))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jb) SCHEDULE(STATIC,jchunk)
      do 4 j=1,jj
      jb=mod(j,jj)+1
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    alist(i,j)=alist(i,j)-(blist(i,jb)-blist(i,j))
c$OMP END PARALLEL DO
c
      return
      end
c
c
      subroutine psmoo(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c --- this routine is set up to smooth data carried at -p- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real,intent(INOUT) :: alist(idm,jdm)
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
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
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
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
c
      subroutine usmoo(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where iu > 0
c --- this routine is set up to smooth data carried at -u- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real,intent(INOUT) :: alist(idm,jdm)
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      do 1 l=1,isu(j)
      do 1 i=ifu(j,l),ilu(j,l)
      ja=mod(j-2+jj,jj)+1
      if (iu(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (iu(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      ia=max( 1,i-1)
      if (iu(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (iu(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine vsmoo(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at poins where iv > 0
c --- this routine is set up to smooth data carried at -v- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real,intent(INOUT) :: alist(idm,jdm)
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      do 1 l=1,isv(j)
      do 1 i=ifv(j,l),ilv(j,l)
      ja=mod(j-2+jj,jj)+1
      if (iv(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (iv(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      do 2 l=1,isv(j)
      do 2 i=ifv(j,l),ilv(j,l)
      ia=max( 1,i-1)
      if (iv(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (iv(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
      subroutine psmooo(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- alist is smoothed array, blist is work array
c --- this routine is set up to smooth data carried at -p- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real,intent(INOUT) :: alist(idm,jdm)
      real blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      ja=mod(j-2+jj,jj)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (ip(i,jb).eq.0) jb=j
      blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
 1    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 3 j=1,jj
      do 3 l=1,isp(j)
      do 3 i=ifp(j,l),ilp(j,l)
      ja=mod(j-2+jj,jj)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (ip(i,jb).eq.0) jb=j
      blist(i,j)=2.*wgt*(alist(i,ja)+alist(i,jb))
 3    continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 4 j=1,jj
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      alist(i,j)=2.*wgt*(blist(ia,j)+blist(ib,j))
 4    continue
c$OMP END PARALLEL DO
c
      return
      end
c
c
      subroutine psmoo4(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c --- this routine is set up to smooth data carried at -p- points
c 
c --- this version works for both cyclic-in-j and noncyclic domains
c 
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real*4,intent(INOUT) :: alist(idm,jdm)
      real*4 blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
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
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
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
c
      subroutine usmoo4(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where iu > 0
c --- this routine is set up to smooth data carried at -u- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real*4,intent(INOUT) :: alist(idm,jdm)
      real*4 blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      do 1 l=1,isu(j)
      do 1 i=ifu(j,l),ilu(j,l)
      ja=mod(j-2+jj,jj)+1
      if (iu(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (iu(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      ia=max( 1,i-1)
      if (iu(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (iu(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine vsmoo4(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at poins where iv > 0
c --- this routine is set up to smooth data carried at -v- points
c 
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM_GLOB
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real*4,intent(INOUT) :: alist(idm,jdm)
      real*4 blist(idm,jdm)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 1 j=1,jj
      do 1 l=1,isv(j)
      do 1 i=ifv(j,l),ilv(j,l)
      ja=mod(j-2+jj,jj)+1
      if (iv(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (iv(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 2 j=1,jj
      do 2 l=1,isv(j)
      do 2 i=ifv(j,l),ilv(j,l)
      ia=max( 1,i-1)
      if (iv(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (iv(ib,j).eq.0) ib=i
 2    alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
c$OMP END PARALLEL DO
      return
      end
c
c
      subroutine zebra(array,idim,ii,jj)
c
c --- find nice contour interval resulting in 7 to 10 contour lines and
c --- draw contours on line printer through the following set of grid points:
c
c         array( 1, 1) . . . . . . . . .  array( 1,jj)
c              .                               .
c              .                               .          (plot will appear
c              .                               .          on paper as shown,
c              .                               .          i down, j across)
c              .                               .
c         array(ii,jj) . . . . . . . . .  array(ii,jj)
c
c --- ii  may be smaller than  idim, the first (row) dimension of 'array'
c --- in the calling program. thus, plotting of partial arrays is possible.
c
      USE HYCOM_SCALARS, only : lp
      implicit none
      integer lgth,idim,ii,jj,i,j
      parameter (lgth=1600)
      integer lp,imn,imx,jmn,jmx,
     .        imnj(lgth),imxj(lgth),jmnj(lgth),jmxj(lgth)
      real sqrt2,contur,q,ratio,amn,amx,amnj(lgth),amxj(lgth)
      common/linepr/lp
      real array(idim,jj)
      data sqrt2/1.414/
c
      write (lp,'(a,3i6)') 'ZEBRA call with arguments',idim,ii,jj
      if (jj.gt.lgth) stop '(insuff. workspace in zebra: increase lgth)'
c
c$OMP PARALLEL DO 
      do 1 j=1,jj
      amxj(j)=-1.e33
      amnj(j)= 1.e33
      imnj(j)=-1
      imxj(j)=-1
      do 1 i=1,ii
      if (amxj(j).lt.array(i,j)) then
        amxj(j)=array(i,j)
        imxj(j)=i
        jmxj(j)=j
      end if
      if (amnj(j).gt.array(i,j)) then
        amnj(j)=array(i,j)
        imnj(j)=i
        jmnj(j)=j
      end if
 1    continue
c$OMP END PARALLEL DO
c
      amx=-1.e33
      amn= 1.e33
      imn=-1
      imx=-1
      jmn=-1
      jmx=-1
      do 2 j=1,jj
      if (amx.lt.amxj(j)) then
        amx=amxj(j)
        imx=imxj(j)
        jmx=jmxj(j)
      end if
      if (amn.gt.amnj(j)) then
        amn=amnj(j)
        imn=imnj(j)
        jmn=jmnj(j)
      end if
 2    continue
c
      if (amx.gt.amn) go to 3
      write (lp,100) array(1,1)
 100  format (//' field to be contoured is constant ...',1pe15.5/)
      return
c
 3    contur=(amx-amn)/6.
      q=10.**int(log10(contur))
      if (contur.lt.1.) q=q/10.
      ratio=contur/q
      if (ratio.gt.sqrt2*5.)  contur=q*10.
      if (ratio.le.sqrt2*5.)  contur=q*5.
      if (ratio.le.sqrt2*2.)  contur=q*2.
      if (ratio.le.sqrt2)     contur=q
      write (lp,101) contur,amn,imn,jmn,amx,imx,jmx
 101  format ('contour interval in plot below is',1pe9.1,
     . 6x,'min =',e11.3,'  at',2i5/48x,'max =',e11.3,'  at',2i5)
      call digplt(array,idim,ii,jj,contur)
c
      return
      end
c
c
      subroutine digplt(array,idim,ii,jj,dec)
c
c --- simulate a contour line plot on the printer
c
      USE HYCOM_SCALARS, only : lp
      implicit none
      integer lp
      common/linepr/lp
c
      integer idim,ii,jj
      real array(idim,jj),dec
      integer i,j,n,k,nchar,ia,ja
      real ratio,xinc,yinc,x,y,dx,dy,dxdy,value
      character*1 digit(130),dig(20)
      data dig/'0',' ','1',' ','2',' ','3',' ','4',' ',
     .         '5',' ','6',' ','7',' ','8',' ','9',' '/
c
c     nchar = number of character increments in 'j' direction
c     ratio = character width / line spacing
c
      data nchar/74/,ratio/.58/
      xinc=float(jj-1)/(float(nchar)*ratio)
      yinc=float(jj-1)/ float(nchar)
      k=float(nchar)*ratio*float(ii-1)/float(jj-1)+1.00001
      do 1 i=1,k
      x=1.+float(i-1)*xinc
      ia=min(ii-1,int(x))
      dx=x-float(ia)
      do 2 j=1,nchar+1
      y=1.+float(j-1)*yinc
      ja=min(jj-1,int(y))
      dy=y-float(ja)
      dxdy=dx*dy
      value=array(ia,ja)*(1.-dx-dy+dxdy)
     .     +array(ia+1,ja)*(dx-dxdy)
     .     +array(ia,ja+1)*(dy-dxdy)
     .     +array(ia+1,ja+1)*dxdy
      n=mod(mod(int(2.*value/dec+sign(.5,value)),20)+20,20)+1
 2    digit(j)=dig(n)
 1    write (lp,100) 'i',' ',(digit(j),j=1,nchar+1),' ','i'
 100  format(1x,130a1)
      return
      end
c
      subroutine ParPsmoo(alist)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist- at points where ip > 0
c --- this routine is set up to smooth data carried at -p- points
c
c --- this version works for both cyclic-in-j and noncyclic domains
c
      USE HYCOM_DIM
      implicit none
!!      include 'dimensions.h'
!!    include "dimension2.h"    ! TNL
      integer i,j,l,ia,ib,ja,jb
c
      real,intent(INOUT) :: alist(idm,J_0H:J_1H)
      real blist(idm,J_0H:J_1H)
      real,parameter :: wgt=.25
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 1 j=J_0,J_1
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
      ja = PERIODIC_INDEX(j-1, jj)
      if (ip(i,ja).eq.0) ja=j
      jb = PERIODIC_INDEX(j+1, jj)
      if (ip(i,jb).eq.0) jb=j
 1    blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ia,ib) SCHEDULE(STATIC,jchunk)
      do 2 j=J_0,J_1
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

