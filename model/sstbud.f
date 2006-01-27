      subroutine sstbud(nproc,text,sstnew)
      implicit none
c
c --- record time-averaged SST tendencies caused by various physical processes
c
c --- nproc    - process counter (normally 1,,...,maxpro, except:)
c                nproc = 99 initializes time integrals
c                nproc = 0 sets sstold = sstnew (initialization)
c --- text     - character string identifying process associated with 'nproc'
c --- sstnew   - sea surface temp. after completion of process 'nproc'
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      integer maxpro,nproc,lgth,no
      parameter (maxpro=6)
      common/sst/sstold,tndcy,protxt
      data no/29/
      real sstold(idm,jdm),sstnew(idm,jdm),q
      real*4 tndcy(idm,jdm,maxpro)
      character text*16,protxt(maxpro)*16,flnm*64,string*20
c
      if (nproc.eq.99) then
        do n=1,maxpro
c$OMP PARALLEL DO SHARED(n)
        do 5 j=1,jdm
        do 5 i=1,idm
 5      tndcy(i,j,n)=0.
c$OMP END PARALLEL DO
        end do
        return
c
      else if (nproc.eq.0) then
c$OMP PARALLEL DO
        do 1 j=1,jdm
        do 1 i=1,idm
 1      sstold(i,j)=sstnew(i,j)
c$OMP END PARALLEL DO
        return
c
      else                                !  nproc =/= 0
        if (nproc.gt.maxpro) stop 'sstbud: maxpro too small'
        protxt(nproc)=text
c
cdiag   write (lp,'(i9,2i4,a,f11.5)') nstep,itest,jtest,protxt(nproc),
cdiag.   sstold(itest,jtest)
c
c$OMP PARALLEL DO
        do 2 j=1,jdm
        do 2 l=1,isp(j)
        do 2 i=ifp(j,l),ilp(j,l)
c
c --- build up time average of sst change attributed to process 'nproc'
        tndcy(i,j,nproc)=tndcy(i,j,nproc)+sstnew(i,j)-sstold(i,j)
c
c --- 'new' sst becomes 'old' sst for next process
 2      sstold(i,j)=sstnew(i,j)
c$OMP END PARALLEL DO
c
cdiag   write (lp,'(33x,f11.5)') sstold(itest,jtest)
      end if
c
      if (mod(time-tavini+.001,diagfq).gt..002 .or. nproc.lt.maxpro)
     .  return
c
c --- output time-averaged sst tencendies attributed to processes 1,...,maxpro
c
      do 3 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 4
 3    continue
      write (lp,*) 'sstbud  --  cannot find slash in ',flnmovt
      stop
 4    write (string,'(a7,i6.6,''-'',i6.6)') 'sstbud.',int(tavini+.001),
     .   int(time+.001)
      flnm=flnmovt(1:lgth)//string
      write (lp,'(a/9x,a)') 'save sst budget in ',flnm
      open (unit=no,file=flnm,status='unknown',form='unformatted')
      write (no) idm,jdm,maxpro,tndcy,protxt
      close (no)
c
      do 7 n=1,maxpro
c
c$OMP PARALLEL DO SHARED(n)
      do 6 j=1,jdm
      do 6 i=1,idm
      util1(i,j)=tndcy(i,j,n)
 6    tndcy(i,j,n)=0.
c$OMP END PARALLEL DO
c
      write (lp,'(2a)') 'shown below: sst change due to ',
     .   protxt(n)
      call zebra(util1,idm,ii1,jj)
 7    continue
c
      return
      end
c
c
c> Revision history:
c>
c> Dec. 2000 - eliminated nproc < 0 option related to surflx prescription
