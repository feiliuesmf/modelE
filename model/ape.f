      real function hyc_pechg1(delp,sig,nunit)
c
c --- calculate change in ocean's available potential energy due to some
c --- process 'X'. call pechg1  b e f o r e  , pechg2  a f t e r  process X.
c --- input: 3-d arrays of  h y b r i d  layer thickness and density anomaly.
c --- results from pechg1 are stored in 'nunit' for later use by pechg2.
c --- use different values of 'nunit' for nested APE process diagnostics.
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real slithk
      integer nscli,nunit,lgth
      parameter (slithk=1.,nscli=7000./slithk)
      real delp(idm,jdm,kdm),sig(idm,jdm,kdm),
     .     pbfore(idm,jdm,kdm+1),pafter(idm,jdm,kdm+1),
     .     praw(idm,jdm,kdm+1),unused(idm,jdm,kdm)
      real weight(kdm),sliwgt(nscli),slitop,
     .     slibot,slisum,prefbe(kdm),prefaf(kdm),varian(0:kdm),
     .     weightj(jdm),sliwgtj(jdm),varianj(jdm)
      character flnm*60
      data varian(0),varian(kdm)/0.,0./
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
 11   praw(i,j,1)=0.
      do 10 k=1,kk
      do 10 i=ifp(j,l),ilp(j,l)
 10   praw(i,j,k+1)=praw(i,j,k)+delp(i,j,k)
c$OMP END PARALLEL DO
c
c --- transform pressure to isopycnic interface pressure
c
      call reflux(uflx,vflx,sig,praw,
     .            unused,unused,unused,pbfore,theta,kdm,kdm)
c
c --- in preparation for determining the flat reference state (i.e., the
c --- unavailable pot. energy), find total mass above each isopycnic interface
c
      do 1 k=1,kk
      weight(k)=0.
c
c$OMP PARALLEL DO
      do 2 j=1,jj
      weightj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
 2    weightj(j)=weightj(j)+pbfore(i,j,k+1)*scp2(i,j)
c$OMP END PARALLEL DO
c
      do 1 j=1,jj
 1    weight(k)=weight(k)+weightj(j)
c
c --- divide basin into shallow horizontal slices and find mass of each slice
c
      do 3 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      sliwgt(n)=0.
c$OMP PARALLEL DO
      do 4 j=1,jj
      sliwgtj(j)=0.
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    sliwgtj(j)=sliwgtj(j)+scp2(i,j)*(min(pbfore(i,j,kk+1),slibot)-
     .                                 min(pbfore(i,j,kk+1),slitop))
c$OMP END PARALLEL DO
      do 3 j=1,jj
 3    sliwgt(n)=sliwgt(n)+sliwgtj(j)
c
      do 5 k=1,kk-1
c --- add slices vertically until sum exceeds combined mass of layers 1...k.
c --- this tells us where bottom of layer k will be when flattened
c
      slisum=0.
      do 6 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      slisum=slisum+sliwgt(n)
      if (slisum.ge.weight(k)) go to 7
 6    continue
      write (*,*) 'k =',k,'  error: slisum < weight',slisum,weight(k)
      go to 5
c
c --- interpolate among slices to get precise depth of flattened interface
c
 7    prefbe(k)=(slibot*(weight(k)-slisum+sliwgt(n))-
     .           slitop*(weight(k)-slisum          ))/sliwgt(n)
ccc      write (*,*) 'k =',k,'  flattened interface:',prefbe(k)/onem
 5    continue
c
c --- save results for later use by pechg2
c
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 13
 14   continue
      write (lp,*) 'ape --  cannot find slash in',flnmovt
      stop
 13   write (flnm,'(a,i2.2)') flnmovt(1:lgth)//'ape.',nunit
      open (unit=nunit,file=flnm,form='unformatted',status='unknown')
      write (nunit) pbfore,prefbe
      close (nunit)
c
c --- determine variance of interface pressure relative to flat reference state
c
      do 12 k=1,kk-1
      varian(k)=0.
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      varianj(j)=0.
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    varianj(j)=varianj(j)+scp2(i,j)*
     .  (pbfore(i,j,k+1)-min(pbfore(i,j,kk+1),prefbe(k)))**2
c$OMP END PARALLEL DO
c
      do 12 j=1,jj
 12   varian(k)=varian(k)+varianj(j)
c
      hyc_pechg1=0.
      do 8 k=1,kk
 8    hyc_pechg1=hyc_pechg1+
     .  (varian(k)-varian(k-1))/(1000.+theta(k)+thbase)
c
c --- report result in units of joules (kg m^2/sec^2)
      hyc_pechg1=.5*hyc_pechg1/g
      return
      end
c
c
      real function hyc_pechg2(delp,sig,nunit)
c
c --- calculate change in ocean's available potential energy due to some
c --- process 'X'. call pechg1  b e f o r e  , pechg2  a f t e r  process X.
c --- input: 3-d arrays of  h y b r i d  layer thickness and density anomaly.
c --- results from pechg1 representing 'before' state are read from 'nunit'.
c --- use different values of 'nunit' for nested APE process diagnostics.
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real slithk
      integer nscli,nunit,lgth
      parameter (slithk=1.,nscli=7000./slithk)
      real delp(idm,jdm,kdm),sig(idm,jdm,kdm),
     .     pbfore(idm,jdm,kdm+1),pafter(idm,jdm,kdm+1),
     .     praw(idm,jdm,kdm+1),unused(idm,jdm,kdm)
      real weight(kdm),sliwgt(nscli),slitop,
     .     slibot,slisum,prefbe(kdm),prefaf(kdm),varian(0:kdm),
     .     weightj(jdm),sliwgtj(jdm),varianj(jdm)
      character flnm*60
      data varian(0),varian(kdm)/0.,0./
      common /pechng/ pbfore,pafter,prefbe,prefaf
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
 11   praw(i,j,1)=0.
      do 10 k=1,kk
      do 10 i=ifp(j,l),ilp(j,l)
 10   praw(i,j,k+1)=praw(i,j,k)+delp(i,j,k)
c$OMP END PARALLEL DO
c
c --- transform pressure to isopycnic interface pressure
c
      call reflux(uflx,vflx,sig,praw,
     .            unused,unused,unused,pafter,theta,kdm,kdm)
c
c --- in preparation for determining the flat reference state (i.e., the
c --- unavailable pot. energy), find total mass above each isopycnic interface
c
      do 1 k=1,kk
      weight(k)=0.
c
c$OMP PARALLEL DO
      do 2 j=1,jj
      weightj(j)=0.
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
 2    weightj(j)=weightj(j)+pafter(i,j,k+1)*scp2(i,j)
c$OMP END PARALLEL DO
c
      do 1 j=1,jj
 1    weight(k)=weight(k)+weightj(j)
c
c --- divide basin into shallow horizontal slices and find mass of each slice
c
      do 3 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      sliwgt(n)=0.
c$OMP PARALLEL DO
      do 4 j=1,jj
      sliwgtj(j)=0.
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    sliwgtj(j)=sliwgtj(j)+scp2(i,j)*(min(pafter(i,j,kk+1),slibot)-
     .                                 min(pafter(i,j,kk+1),slitop))
c$OMP END PARALLEL DO
      do 3 j=1,jj
 3    sliwgt(n)=sliwgt(n)+sliwgtj(j)
c
      do 5 k=1,kk-1
c --- add slices vertically until sum exceeds combined mass of layers 1...k.
c --- this tells us where bottom of layer k will be when flattened
c
      slisum=0.
      do 6 n=1,nscli
      slitop=float(n-1)*slithk*onem
      slibot=float(n  )*slithk*onem
      slisum=slisum+sliwgt(n)
      if (slisum.ge.weight(k)) go to 7
 6    continue
      write (*,*) 'k =',k,'  error: slisum < weight',slisum,weight(k)
      go to 5
c
c --- interpolate among slices to get precise depth of flattened interface
c
 7    prefaf(k)=(slibot*(weight(k)-slisum+sliwgt(n))-
     .           slitop*(weight(k)-slisum          ))/sliwgt(n)
ccc      write (*,*) 'k =',k,'  flattened interface:',prefaf(k)/onem
 5    continue
c
c --- read previously created file representing 'before' state
c
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 13
 14   continue
      write (lp,*) 'ape --  cannot find slash in',flnmovt
      stop
 13   write (flnm,'(a,i2.2)') flnmovt(1:lgth)//'ape.',nunit
      open (unit=nunit,file=flnm,form='unformatted',status='old')
      read (nunit) pbfore,prefbe
      rewind (nunit)
c
c --- write out 'after' state (for potential future use as new 'before' state)
c
      write (nunit) pafter,prefaf
      close (nunit)
c
c --- determine 'before/after' difference of interface pressure variances
c
      do 12 k=1,kk-1
      varian(k)=0.
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      varianj(j)=0.
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    varianj(j)=varianj(j)+scp2(i,j)*
     . ((pafter(i,j,k+1)-min(pafter(i,j,kk+1),prefaf(k)))**2
     . -(pbfore(i,j,k+1)-min(pbfore(i,j,kk+1),prefbe(k)))**2)
c$OMP END PARALLEL DO
c
      do 12 j=1,jj
 12   varian(k)=varian(k)+varianj(j)
c
      hyc_pechg2=0.
      do 8 k=1,kk
 8    hyc_pechg2=hyc_pechg2
     .  +(varian(k)-varian(k-1))/(1000.+theta(k)+thbase)
c
c --- report result in units of watts (kg m^2/sec^3)
      hyc_pechg2=.5*hyc_pechg2/(g*delt1)
      return
      end
c
c
c> Revision history:
c>
c> Dec. 2004 - fixed bug in loop 9 (excluded interfaces on shallow bottom)

