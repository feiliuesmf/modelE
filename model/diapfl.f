      subroutine diapfl(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real flxu(idm,kdm),flxl(idm,kdm),pdot(idm,kdm),flngth(idm,kdm),
     .     ennsq,alfa,beta,q,qmin,qmax,amount,salt,froglp,small,delp,
     .     trflxu(idm,0:kdm+1,ntrcr),trflxl(idm,0:kdm+1,ntrcr),
     .     cliptr(idm,ntrcr),
     .     tflxu(idm,0:kdm+1),tflxl(idm,0:kdm+1),clipt(idm),
     .     sflxu(idm,0:kdm+1),sflxl(idm,0:kdm+1),clips(idm),
     .     told(idm,2),sold(idm,2),trold(idm,2,ntrcr),scale(idm),
     .     totem(idm),tosal(idm),totra(idm),tndcyt,tndcys,tndtra
      integer kmin(idm),kmax(idm),ka,kan,nt
      character text*20
      data small/1.e-22/
      real sigocn,dsigdt,dsigds
      external sigocn,dsigdt,dsigds
c
      if (diapyc.eq.0. .or. mod(nstep,mixfrq).gt.1) return
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
c$OMP PARALLEL DO PRIVATE(kn,q,ennsq,alfa,beta,totem,tosal,totra,sold,
c$OMP. told,trold,tflxl,tflxu,sflxl,sflxu,trflxl,trflxu,kmin,kmax,
c$OMP. flngth,flxu,flxl,pdot,clipt,clips,cliptr,ka,kan,delp,amount,
c$OMP. qmax,qmin,tndcyt,tndcys,tndtra,scale)
      do 31 j=1,jj
      do 31 l=1,isp(j)
c
c --- t/s conservation diagnostics (optional):
      do i=ifp(j,l),ilp(j,l)
        totem(i)=0.
        tosal(i)=0.
        totra(i)=0.
        scale(i)=1.e-99
        do k=1,kk
          kn=k+nn
          totem(i)=totem(i)+temp(i,j,kn)*dp(i,j,kn)
          tosal(i)=tosal(i)+saln(i,j,kn)*dp(i,j,kn)
          if (dotrcr) then
            totra(i)=totra(i)+tracer(i,j,k,1)*dp(i,j,kn)
            scale(i)=scale(i)+abs(tracer(i,j,k,1))
          end if
        end do
      end do
c
      do 33 k=1,kk
      kn=k+nn
      do 33 i=ifp(j,l),ilp(j,l)
 33   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
      do 32 i=ifp(j,l),ilp(j,l)
      sold(i,1)=saln(i,j,kk+nn)
      told(i,1)=temp(i,j,kk+nn)
      tflxl(i,   0)=0.
      tflxu(i,kk+1)=0.
      sflxl(i,   0)=0.
      sflxu(i,kk+1)=0.
      if (dotrcr) then
        trold(i,1,:)=tracer(i,j,kk,:)
        trflxl(i,   0,:)=0.
        trflxu(i,kk+1,:)=0.
      end if
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,103) nstep,i,j,
cdiag. '  entering diapfl:  temp    saln    dens    thkns   tracer',
cdiag.  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
cdiag.   dp(i,j,k+nn)/onem,tracer(i,j,k,1),k=1,kk)
 103  format(i9,2i5,a/(33x,i3,3f8.3,f8.2,f9.4))
c
      kmin(i)=kk+1
 32   kmax(i)=1
c
      do 36 k=2,kk
      kn=k+nn
c
      do 36 i=ifp(j,l),ilp(j,l)
c
c --- locate lowest mass-containing layer and upper edge of stratified region
      if (p(i,j,k).lt.p(i,j,kk+1)-onemm)  then
        kmax(i)=k
        if (kmin(i).eq.kk+1 .and. th3d(i,j,kn).gt.th3d(i,j,kn-1)+sigjmp)
     .      kmin(i)=k-1
      end if
 36   continue
c
cdiag if (j.eq.jtest.and.itest.ge.ifp(j,l).and.itest.le.ilp(j,l))
cdiag.  write (lp,'(i9,2i5,a,2i5)') nstep,itest,j,' kmin,kmax =',
cdiag.  kmin(itest),kmax(itest)
c
c --- find buoyancy frequency for each layer
c
      do 43 k=2,kk-1
      kn=k+nn
c
      do 43 i=ifp(j,l),ilp(j,l)
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
c --- ennsq = buoy.freq.^2 / g^2
        ennsq=max(0.,min(th3d(i,j,kn+1)-th3d(i,j,kn  ),
     .                   th3d(i,j,kn  )-th3d(i,j,kn-1)))
     .    /max(p(i,j,k+1)-p(i,j,k),onemm)
c --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
c --- (dimensions of flngth: length in pressure units)
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc / buoyancy frequency
        flngth(i,k)=diapyc*sqrt(ennsq) * baclin*froglp * onem
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc
ccc        flngth(i,k)=diapyc*ennsq*g * baclin*froglp * onem
c -----------------------------------------------------------------------
c
      end if
 43   continue
c
c --- find t/s fluxes at the upper and lower interface of each layer
c --- (compute only the part common to t and s fluxes)
c
      do 37 k=1,kk
      kn=k+nn
c
      do 37 i=ifp(j,l),ilp(j,l)
      flxu(i,k)=0.
      flxl(i,k)=0.
c
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
c
        alfa=-thref*dsigdt(temp(i,j,kn),saln(i,j,kn))
        beta= thref*dsigds(temp(i,j,kn),saln(i,j,kn))
c
        flxu(i,k)=flngth(i,k)/
     .    max(beta*(saln(i,j,kn)-saln(i,j,kn-1))
     .       -alfa*(temp(i,j,kn)-temp(i,j,kn-1)),small)
        flxl(i,k)=flngth(i,k)/
     .    max(beta*(saln(i,j,kn+1)-saln(i,j,kn))
     .       -alfa*(temp(i,j,kn+1)-temp(i,j,kn)),small)
c
        q=min(1.,.5*min(p(i,j,k)-p(i,j,k-1),p(i,j,k+2)-p(i,j,k+1))/
     .    max(flxu(i,k),flxl(i,k),epsil))
c
cdiag   if (q.ne.1.) write (lp,'(i9,2i5,i3,a,1p,2e10.2,0p,2f7.2,f5.2)') 
cdiag.    nstep,i,j,k,' flxu/l,dpu/l,q=',flxu(i,k),flxl(i,k),
cdiag.    (p(i,j,k)-p(i,j,k-1))/onem,(p(i,j,k+2)-p(i,j,k+1))/onem,q
c
        flxu(i,k)=flxu(i,k)*q
        flxl(i,k)=flxl(i,k)*q
c
      end if				!  kmin < k < kmax
c
cdiag if (i.eq.itest.and.j.eq.jtest.and.k.ge.kmin(i).and.k.le.kmax(i))
cdiag.   write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)')nstep,
cdiag.   i,j,k,'thknss   temp   saln    flngth      flxu      flxl',
cdiag.   dp(i,j,kn)/onem,temp(i,j,kn),saln(i,j,kn),flngth(i,k),
cdiag.   flxu(i,k)/onem,flxl(i,k)/onem
c
 37   continue
c
c --- determine mass flux -pdot- implied by t/s fluxes.
c
      do 38 k=1,kk
      kn=k+nn
c
      do 38 i=ifp(j,l),ilp(j,l)
      pdot(i,k)=0.
      if (k.gt.kmin(i) .and. k.le.kmax(i))
     .    pdot(i,k)=flxu(i,k)-flxl(i,k-1)
 38   continue
c
c --- convert flxu,flxl into actual t/s (and tracer) fluxes
c --- cut tracer fluxes in half to compensate for execution at 2 consecutive
c --- time steps
c
      do 35 k=1,kk
      kn=k+nn
c
      do 35 i=ifp(j,l),ilp(j,l)
      tflxu(i,k)=0.
      tflxl(i,k)=0.
      sflxu(i,k)=0.
      sflxl(i,k)=0.
      if (dotrcr) then
        trflxu(i,k,:)=0.
        trflxl(i,k,:)=0.
      end if
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
        tflxu(i,k)=flxu(i,k)*temp(i,j,kn-1)
        tflxl(i,k)=flxl(i,k)*temp(i,j,kn+1)
        sflxu(i,k)=flxu(i,k)*saln(i,j,kn-1)
        sflxl(i,k)=flxl(i,k)*saln(i,j,kn+1)
        if (dotrcr) then
          trflxu(i,k,:)=flxu(i,k)*tracer(i,j,k-1,:)
          trflxl(i,k,:)=flxl(i,k)*tracer(i,j,k+1,:)
        end if
      end if
 35   continue
c
      do 34 i=ifp(j,l),ilp(j,l)
      if (dotrcr) cliptr(i,:)=0.
      clipt(i)=0.
 34   clips(i)=0.
c
c --- update interface pressure and layer temperature/salinity
      do 39 k=kk,1,-1
      kn=k+nn
      ka=max(1,k-1)
      kan=ka+nn
c
      do 39 i=ifp(j,l),ilp(j,l)
      sold(i,2)=sold(i,1)
      sold(i,1)=saln(i,j,kn)
      told(i,2)=told(i,1)
      told(i,1)=temp(i,j,kn)
      if (dotrcr) then
        trold(i,2,:)=trold(i,1,:)
        trold(i,1,:)=tracer(i,j,k,:)
      end if
c
      dpold(i,j,k)=dp(i,j,kn)
      p(i,j,k)=p(i,j,k)-pdot(i,k)
      dp(i,j,kn)=p(i,j,k+1)-p(i,j,k)
c
      if (k.ge.kmin(i) .and. k.le.kmax(i)) then
        delp=dp(i,j,kn)
        if (delp.gt.0.) then
          amount=temp(i,j,kn)*dpold(i,j,k)
     .      -(tflxu(i,k+1)-tflxu(i,k)+tflxl(i,k-1)-tflxl(i,k))
          q=amount
          qmax=max(temp(i,j,kan),told(i,1),told(i,2))
          qmin=min(temp(i,j,kan),told(i,1),told(i,2))
          amount=max(qmin*delp,min(amount,qmax*delp))
          clipt(i)=clipt(i)+(q-amount)
          temp(i,j,kn)=amount/delp
c
          amount=saln(i,j,kn)*dpold(i,j,k)
     .      -(sflxu(i,k+1)-sflxu(i,k)+sflxl(i,k-1)-sflxl(i,k))
          q=amount
          qmax=max(saln(i,j,kan),sold(i,1),sold(i,2))
          qmin=min(saln(i,j,kan),sold(i,1),sold(i,2))
          amount=max(qmin*delp,min(amount,qmax*delp))
          clips(i)=clips(i)+(q-amount)
          saln(i,j,kn)=amount/delp
c
          if (dotrcr) then
            do nt=1,ntrcr
              amount=tracer(i,j,k,nt)*dpold(i,j,k)
     .            -(trflxu(i,k+1,nt)-trflxu(i,k,nt)
     .             +trflxl(i,k-1,nt)-trflxl(i,k,nt))
              q=amount
              qmax=max(tracer(i,j,ka,nt),trold(i,1,nt),trold(i,2,nt))
              qmin=min(tracer(i,j,ka,nt),trold(i,1,nt),trold(i,2,nt))
              amount=max(qmin*delp,min(amount,qmax*delp))
              cliptr(i,nt)=cliptr(i,nt)+(q-amount)
              tracer(i,j,k,nt)=amount/delp
            end do
          end if			!  dotrcr
        end if				!  delp > 0
      end if
 39   continue
c
      do 30 i=ifp(j,l),ilp(j,l)
      if (dotrcr) cliptr(i,:)=cliptr(i,:)/pbot(i,j)
      clipt(i)=clipt(i)/pbot(i,j)
 30   clips(i)=clips(i)/pbot(i,j)
c
      do 41 k=1,kk
      kn=k+nn
      do 41 i=ifp(j,l),ilp(j,l)
c
c --- restore 'clipped' t/s amount to column
      temp(i,j,kn)=temp(i,j,kn)+clipt(i)
      saln(i,j,kn)=saln(i,j,kn)+clips(i)
      th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))
      if (dotrcr) tracer(i,j,k,:)=tracer(i,j,k,:)+cliptr(i,:)
c
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-dpold(i,j,k))	! diapyc.flux
c --- make sure p is computed from dp, not the other way around (roundoff!)
 41   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- t/s conservation diagnostics (optional):
      do i=ifp(j,l),ilp(j,l)
        tndcys=-tosal(i)
        tndcyt=-totem(i)
        tndtra=-totra(i)
        do k=1,kk
          kn=k+nn
          tndcys=tndcys+saln(i,j,kn)*dp(i,j,kn)
          tndcyt=tndcyt+temp(i,j,kn)*dp(i,j,kn)
          if (dotrcr) tndtra=tndtra+tracer(i,j,k,1)*dp(i,j,kn)
        end do
        if (abs(tndcyt).gt.acurcy*10.*pbot(i,j))
     .  write (lp,100) i,j,
     .   '  diapfl - bad temp.intgl.',totem(i),tndcyt,clipt(i)
        if (abs(tndcys).gt.acurcy*35.*pbot(i,j))
     .  write (lp,100) i,j,
     .   '  diapfl - bad saln.intgl.',tosal(i),tndcys,clips(i)
        if (dotrcr) then
          if (abs(tndtra)*kk.gt.acurcy*scale(i)*pbot(i,j))
     .    write (lp,100) i,j,
     .     '  diapfl - bad trcr.intgl.',totra(i),tndtra,cliptr(i,1)
        end if
 100    format(2i5,a,1p,e16.8,2e13.5)
      end do
c
cdiag do 31 i=ifp(j,l),ilp(j,l)
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,103) nstep,i,j,
cdiag. '  exiting  diapfl:  temp    saln    dens    thkns   tracer',
cdiag.  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
cdiag.   dp(i,j,k+nn)/onem,tracer(i,j,k,1),k=1,kk)
c
 31   continue
c$OMP END PARALLEL DO
c
ccc   write (lp,'(i9,7x,1p,e9.2,a)') nstep,salt*1.e-6/g,
ccc  .  ' kg salt added in diapfl'
c
      call dpudpv(nn)
c
      if (dotrcr) write (lp,'(a)') 'tracer diapycnal mixing done'
      return
      end
c
c
c> Revision history:
c>
c> Mar. 2000 - conversion to SI units
c> May  2000 - converted T/S advection equations to flux form
c> Apr. 2001 - eliminated stmt_funcs.h
c> Sep. 2003 - added code to return 'clipped' tracer amount to column (loop 41)
c> Sep. 2003 - added logical switch to enable/disable tracer diffusion
c> Feb. 2005 - added multiple tracer capability
