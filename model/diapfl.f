      subroutine diapfl(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.1
      implicit none
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      real pres(kdm+1),delp(kdm),olddp(kdm),ssal(kdm),ttem(kdm),
     .     dens(kdm),virt(kdm),trac(kdm),flxu(kdm),flxl(kdm),pdot(kdm),
     .     flngth(kdm),ennsq,alfa,beta,q,qmin,qmax,amount,salt,
     .     froglp,small,thk,trflxu(0:kdm+1),trflxl(0:kdm+1),
     .     sflxu(0:kdm+1),sflxl(0:kdm+1),clip,sold(2),trold(2),
     .     tosal,tndcys			!  col.integrals (diag.use only)
      real kappaf,tofsig,dsigdt,dsigds
      external kappaf,tofsig,dsigdt,dsigds
      integer kmax,ka,kan
      character text*20
      parameter (small=1.e-22)
c
      if (diapyc.eq.0. .or. (mod(nstep  ,mixfrq).ne.0 .and.
     .                       mod(nstep+1,mixfrq).ne.0)) return
c
c --- ----------------
c --- diapycnal mixing
c --- ----------------
c
c --- if mixfrq > 1, apply mixing algorithm to both time levels
      froglp=max(2,mixfrq)
c
ccc   salt=0.
c
c$OMP PARALLEL DO PRIVATE(flxu,flxl,pdot,flngth,ennsq,alfa,beta,q,
c$OMP. qmin,qmax,amount,salt,trflxu,trflxl,sflxu,sflxl,clip,thk,
c$OMP. pres,delp,olddp,ssal,ttem,dens,virt,trac,sold,trold,tosal,
c$OMP. tndcys,kmax,ka,kan,kn)
      do 3 j=1,jj
      do 3 l=1,isp(j)
c
      do 4 k=1,kk
      do 4 i=ifp(j,l),ilp(j,l)
 4    p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c
      do 3 i=ifp(j,l),ilp(j,l)
c
c --- extract 1-d column from 3-d fields
c
      tosal=0.
      pres(1)=0.
      do 1 k=1,kk
      kn=k+nn
      dpold(i,j,k)=dp(i,j,kn)
      tosal=tosal+saln(i,j,kn)*dp(i,j,kn)
      if (p(i,j,k).lt.p(i,j,kk+1)) then
        kmax=k
        delp(k)=dp(i,j,kn)
        olddp(k)=delp(k)
        ssal(k)=saln(i,j,kn)
        ttem(k)=temp(i,j,kn)
        dens(k)=th3d(i,j,kn)
        virt(k)=dens(k)+thermb(i,j,kn)
        trac(k)=tracer(i,j,k)
        pres(k+1)=pres(k)+delp(k)
      end if
 1    continue
c
      sold(1)=ssal(kmax)
      trold(1)=trac(kmax)
      sflxl(     0)=0.
      sflxu(kmax+1)=0.
      trflxl(     0)=0.
      trflxu(kmax+1)=0.
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag. write (lp,103) nstep,i,j,
cdiag.  '  entering diapfl:  temp    saln    dens    thkns',
cdiag.  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn)+thbase,
cdiag.       dp(i,j,k+nn)/onem,k=1,kk)
 103  format(i9,2i5,a/(33x,i3,3f8.3,f8.2))
c
c --- find buoyancy frequency for each layer
c
      do 43 k=2,kmax-1
c --- ennsq = buoy.freq.^2 / g^2
        ennsq=max(0.,min(virt(k+1)-virt(k  ),
     .                   virt(k  )-virt(k-1)))
     .    /max(pres(k+1)-pres(k),onecm)
c --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
c --- (dimensions of flngth: length in pressure units)
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc / buoyancy frequency
        flngth(k)=diapyc*sqrt(ennsq) * baclin*froglp * onem
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc
ccc        flngth(k)=diapyc*ennsq*g * baclin*froglp * onem
c -----------------------------------------------------------------------
 43   continue
c
c --- find t/s fluxes at the upper and lower interface of each layer
c --- (compute only the flux part common to t and s fluxes)
c
      flxu(   1)=0.
      flxl(   1)=0.
      flxu(kmax)=0.
      flxl(kmax)=0.
c
      do 37 k=2,kmax-1
      alfa=-thref*dsigdt(ttem(k),ssal(k))
      beta= thref*dsigds(ttem(k),ssal(k))
c
      flxu(k)=flngth(k)/
     .  max(beta*(ssal(k)-ssal(k-1))
     .     -alfa*(ttem(k)-ttem(k-1)),small)
      flxl(k)=flngth(k)/
     .  max(beta*(ssal(k+1)-ssal(k))
     .     -alfa*(ttem(k+1)-ttem(k)),small)
c
      q=min(1.,.5*min(pres(k)-pres(k-1),pres(k+2)-pres(k+1))/
     .  max(flxu(k),flxl(k),epsil))
c
ccc        if (q.ne.1.) write (lp,'(i9,2i5,i3,a,1p,2e10.2,0p,2f7.2,f5.2)') 
ccc     .    nstep,k,' flxu/l,dpu/l,q=',flxu(k),flxl(k),
ccc     .    (pres(k)-pres(k-1))/onem,(pres(k+2)-pres(k+1))/onem,q
c
      flxu(k)=flxu(k)*q
      flxl(k)=flxl(k)*q
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag.   write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)') nstep,
cdiag.   i,j,k,'thknss   temp   saln    flngth      flxu      flxl',
cdiag.   delp(k)/onem,ttem(k),ssal(k),flngth(k),
cdiag.   flxu(k)/onem,flxl(k)/onem
c
 37   continue
c
c --- determine mass flux -pdot- implied by t/s fluxes.
c
      pdot(1)=0.
      do 38 k=2,kmax
 38   pdot(k)=flxu(k)-flxl(k-1)
c
c --- convert flxu,flxl into actual salt (and tracer) fluxes
c --- cut tracer fluxes in half to compensate for execution at 2 consecutive
c --- time steps
c
      sflxu(   1)=0.
      sflxl(   1)=0.
      sflxu(kmax)=0.
      sflxl(kmax)=0.
      trflxu(   1)=0.
      trflxl(   1)=0.
      trflxu(kmax)=0.
      trflxl(kmax)=0.
c
      do 35 k=2,kmax-1
      sflxu(k)=flxu(k)*ssal(k-1)
      trflxu(k)=flxu(k)*trac(k-1) * .5
c
      sflxl(k)=flxl(k)*ssal(k+1)
      trflxl(k)=flxl(k)*trac(k+1) * .5
 35   continue
c
      clip=0.
c
c --- update interface pressure
      do 39 k=kmax,1,-1
      ka=max(1,k-1)
c
      sold(2)=sold(1)
      sold(1)=ssal(k)
      trold(2)=trold(1)
      trold(1)=trac(k)
c
      pres(k)=pres(k)-pdot(k)
      delp(k)=pres(k+1)-pres(k)
c
      thk=delp(k)
      if (thk.gt.0.) then
        amount=ssal(k)*olddp(k)
     .    -(sflxu(k+1)-sflxu(k)+sflxl(k-1)-sflxl(k))
        q=amount
        qmax=max(ssal(ka),sold(1),sold(2))
        qmin=min(ssal(ka),sold(1),sold(2))
        amount=max(qmin*thk,min(amount,qmax*thk))
        clip=clip+(q-amount)
        ssal(k)=amount/thk
c
        amount=trac(k)*olddp(k)
     .    -(trflxu(k+1)-trflxu(k)+trflxl(k-1)-trflxl(k))
        q=amount
        qmax=max(trac(ka),trold(1),trold(2))
        qmin=min(trac(ka),trold(1),trold(2))
        amount=max(qmin*thk,min(amount,qmax*thk))
        trac(k)=amount/thk
      end if
 39   continue
c
c --- put 1-d column back into 3-d field
c
      do 2 k=1,kmax
      kn=k+nn
      dp(i,j,kn)=delp(k)
      saln(i,j,kn)=ssal(k)
 2    tracer(i,j,k)=trac(k)
c
      clip=clip/pbot(i,j)
c
      do 41 k=1,kk
      kn=k+nn
c
c --- return 'clipped' amount to column
      saln(i,j,kn)=saln(i,j,kn)+clip
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-dpold(i,j,k))	! diapyc.flux
      p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- update temperature
      temp(i,j,kn)=tofsig(th3d(i,j,kn)+thbase,saln(i,j,kn))
ccc      thermb(i,j,kn)=kappaf(temp(i,j,kn),saln(i,j,kn),p(i,j,k))
 41   continue
c
c --- column conservation diagnostics (optional):
      tndcys=-tosal
      tosal=35.*pbot(i,j)
      do k=1,kk
        kn=k+nn
        tndcys=tndcys+saln(i,j,kn)*dp(i,j,kn)
      end do
      if (abs(tndcys).gt.acurcy*tosal)
     . write (lp,100) i,j,'  diapfl - bad saln.intgl.:',
     .  tosal,tndcys,clip
 100    format(2i4,a,1p,e16.8,2e13.5)
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,103) nstep,i,j,
cdiag.  '  exiting  diapfl:  temp    saln    dens    thkns',
cdiag.  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn)+thbase,
cdiag.       dp(i,j,k+nn)/onem,k=1,kk)
c
 3    continue
c$OMP END PARALLEL DO
c
ccc   write (lp,'(i9,7x,1p,e9.2,a)') nstep,salt*1.e-6/g,
ccc  .  ' kg salt added in diapfl'
c
      call dpudpv(nn)
c
      return
      end
c
c
c> Revision history:
c>
c> May  2000 - converted T/S advection equations to flux form
c> Mar. 2001 - eliminated stmt_funcs.h
c> Jan. 2002 - conversion to SI units
c> Ape. 2002 - fixed formula for diagnosing -diaflx- (loop 41)
c> Mar. 2003 - n-square now computed from virt.pot.density 
c> Mar. 2003 - conversion to 'column physics' format (i,j loop outside)
