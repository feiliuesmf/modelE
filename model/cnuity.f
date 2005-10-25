      subroutine cnuity(m,n,mm,nn,k1m,k1n)
c
c --- micom version 2.9
c --- hycom version 0.9
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      integer mask(idm,jdm),jcyc,iz,jz
      real pold(idm,jdm),q,dpmin,dpmn(jdm),clip,flxhi,flxlo,dtinv,old,
     .     hymar,pbot1,pbot2,p1,p2,hyc_pechg1,hyc_pechg2
ccc  .    ,boost
      character text*20
      external hyc_pechg1,hyc_pechg2
c
ccc   boost(pbot1,pbot2,p1,p2)=max(1.,1.5-min(pbot1-p1,pbot2-p2)
ccc  .  /min(3.*tenm,.125*(pbot1+pbot2)))
c
c --- ------------------------------------------------------
c --- continuity equation (flux-corrected transport version)
c --- ------------------------------------------------------
c
c$OMP PARALLEL DO
      do 41 j=1,jj
c
      do 74 l=1,isp(j)
      do 74 i=ifp(j,l),ilp(j,l)
 74   pold(i,j)=0.
c
      do 40 l=1,isu(j)
      do 40 i=ifu(j,l),ilu(j,l)
 40   utotn(i,j)=0.
c
      do 41 l=1,isv(j)
      do 41 i=ifv(j,l),ilv(j,l)
 41   vtotn(i,j)=0.
c$OMP END PARALLEL DO
c
      do 76 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- uflux/vflux = low-order (diffusive) mass fluxes at old time level.
c --- uflux2/vflux2 = 'antidiffusive' fluxes, defined as high-order minus low-
c --- order fluxes. high-order fluxes are second-order in space, time-centered.
c
c$OMP PARALLEL DO PRIVATE(ja,q)
      do 12 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 11 l=1,isu(j)
      do 11 i=ifu(j,l),ilu(j,l)
      utotm(i,j)=(u(i,j,km)+ubavg(i,j,m))*scuy(i,j)
      if (utotm(i,j).ge.0.) then
        q=min(dp(i-1,j,kn),max(0.,depthu(i,j)-pold(i-1,j)))
      else
        q=min(dp(i  ,j,kn),max(0.,depthu(i,j)-pold(i  ,j)))
      end if
      uflux(i,j)=utotm(i,j)*q
      uflux2(i,j)=utotm(i,j)*dpu(i,j,km)-uflux(i,j)
 11   uflx(i,j,k)=uflux(i,j)
c
      do 12 l=1,isv(j)
      do 12 i=ifv(j,l),ilv(j,l)
      vtotm(i,j)=(v(i,j,km)+vbavg(i,j,m))*scvx(i,j)
      if (vtotm(i,j).ge.0.) then
        q=min(dp(i,ja ,kn),max(0.,depthv(i,j)-pold(i,ja )))
      else
        q=min(dp(i,j  ,kn),max(0.,depthv(i,j)-pold(i,j  )))
      end if
      vflux(i,j)=vtotm(i,j)*q
      vflux2(i,j)=vtotm(i,j)*dpv(i,j,km)-vflux(i,j)
 12   vflx(i,j,k)=vflux(i,j)
c$OMP END PARALLEL DO
c
c --- advance -dp- field using low-order (diffusive) flux values
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 19 j=1,jj
      jb=mod(j     ,jj)+1
      dpmn(j)=999.
      do 19 l=1,isp(j)
      do 19 i=ifp(j,l),ilp(j,l)
      dpold(i,j,k)=dp(i,j,kn)
      pold(i,j)=pold(i,j)+dp(i,j,kn)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
 19   dpmn(j)=min(dpmn(j),dp(i,j,kn))
c$OMP END PARALLEL DO
c
ccc      do j=1,jdm
ccc      do i=1,idm
ccc      mask(i,j)=iu(i,j)
ccc      if (i.gt. 1 ) mask(i,j)=mask(i,j)+iu(i-1,j)
ccc      if (i.lt.idm) mask(i,j)=mask(i,j)+iu(i+1,j)
ccc      end do
ccc      end do
ccc      write (text,'(a9,i3,i8)') 'uflux  k=',k,nstep
ccc      call compare(uflux,mask,text)
ccc      do j=1,jdm
ccc      do i=1,idm
ccc      mask(i,j)=iv(i,j)
ccc      if (j.gt. 1 ) mask(i,j)=mask(i,j)+iv(i,ja )
ccc      if (j.lt.jdm) mask(i,j)=mask(i,j)+iv(i,jb )
ccc      end do
ccc      end do
ccc      write (text,'(a9,i3,i8)') 'vflux  k=',k,nstep
ccc      call compare(vflux,mask,text)
c
      dpmin=999.
c$OMP PARALLEL DO REDUCTION(min:dpmin)
      do 191 j=1,jj
 191  dpmin=min(dpmin,dpmn(j))
c$OMP END PARALLEL DO
c
      if (dpmin.lt.-onem) then
        do 190 j=1,jj
        do 190 l=1,isp(j)
        do 190 i=ifp(j,l),ilp(j,l)
        if (dp(i,j,kn).eq.dpmin) then
          write (lp,100) nstep,i,j,k,19,dpmin/onem
 100      format (i9,' i,j,k=',2i5,i3,' neg. dp (m) in loop ',i2,f9.2)
          iz=i
          jz=j
        end if
 190    continue
        call stencl(iz,jz,k,nn)
      end if
c
cdiag write (lp,*) 'time step',nstep,'    layer',k
cdiag do jcyc=jtest-1,jtest+1
cdiag j =mod(jcyc-1+jj,jj)+1
cdiag ja=mod(jcyc-2+jj,jj)+1
cdiag jb=mod(jcyc     ,jj)+1
cdiag do i=itest-1,itest+1
cdiag write (lp,101) i,j,k,'old thknss','mid thknss,vel.',
cdiag.'new thknss,fluxes',
cdiag.dpold(i-1,j,k)/onem,u(i,j,km)+ubavg(i,j,m),uflux(i,j),
cdiag.dpold(i,ja,k)/onem,dpold(i,j,k)/onem,dpold(i,jb,k)/onem,
cdiag.v(i,j,km)+vbavg(i,j,m),dp(i,j,km)/onem,v(i,jb,km)+vbavg(i,jb,m),
cdiag.vflux(i,j),dp(i,j,kn)/onem,vflux(i,jb),
cdiag.dpold(i+1,j,k)/onem,u(i,jb,km)+ubavg(i,jb,m),uflux(i+1,j)
cdiag end do
cdiag end do
 101  format (2i5,i3,2x,a,7x,a,10x,a/0p,f14.1,f26.1,1p,e30.2/
     .   0p,3f7.1,f12.1,f7.1,f6.1,1p,e15.2,0p,f8.1,1p,e10.2/
     .   0p,f14.1,f26.1,1p,e30.2)
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -dp- to the sum of all incoming (outgoing) fluxes
c
c$OMP PARALLEL DO PRIVATE(ia,ib,ja,jb)
      do 26 j=1,jj
      do 26 l=1,isp(j)
      do 26 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      ja=mod(j-2+jj,jj)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (ip(i,jb).eq.0) jb=j
      util1(i,j)=max(dp(i,j,kn),dp(ia,j,kn),dp(ib,j,kn),
     .                          dp(i,ja,kn),dp(i,jb,kn))
      util2(i,j)=max(0.,
     .           min(dp(i,j,kn),dp(ia,j,kn),dp(ib,j,kn),
     .                          dp(i,ja,kn),dp(i,jb,kn)))
c
      jb=mod(j     ,jj)+1
      util1(i,j)=(util1(i,j)-dp(i,j,kn))
     ./((max(0.,uflux2(i,j))-min(0.,uflux2(i+1,j))
     .  +max(0.,vflux2(i,j))-min(0.,vflux2(i,jb ))+epsil)
     .*delt1*scp2i(i,j))
c
 26   util2(i,j)=(util2(i,j)-dp(i,j,kn))
     ./((min(0.,uflux2(i,j))-max(0.,uflux2(i+1,j))
     .  +min(0.,vflux2(i,j))-max(0.,vflux2(i,jb ))-epsil)
     .*delt1*scp2i(i,j))
c$OMP END PARALLEL DO
c
c --- limit antidiffusive fluxes
c --- (keep track in -utotn,vtotn- of discrepancy between high-order
c --- fluxes and the sum of low-order and clipped antidiffusive fluxes.
c --- this will be used later to restore nondivergence of barotropic flow)
c
c$OMP PARALLEL DO PRIVATE(ja,clip)
      do 29 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 28 l=1,isu(j)
      do 28 i=ifu(j,l),ilu(j,l)
      if (uflux2(i,j).ge.0.) then
        clip=min(1.,util1(i,j),util2(i-1,j))
      else
        clip=min(1.,util2(i,j),util1(i-1,j))
      end if
      utotn(i,j)=utotn(i,j)+uflux2(i,j)*(1.-clip)
      uflux(i,j)=uflux2(i,j)*clip
 28   uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
c
      do 29 l=1,isv(j)
      do 29 i=ifv(j,l),ilv(j,l)
      if (vflux2(i,j).ge.0.) then
        clip=min(1.,util1(i,j),util2(i,ja ))
      else
        clip=min(1.,util2(i,j),util1(i,ja ))
      end if
      vtotn(i,j)=vtotn(i,j)+vflux2(i,j)*(1.-clip)
      vflux(i,j)=vflux2(i,j)*clip
 29   vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
c$OMP END PARALLEL DO
c
c --- evaluate effect of antidiffusive fluxes on -dp- field
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 15 j=1,jj
      jb=mod(j     ,jj)+1
      dpmn(j)=999.
      do 15  l=1,isp(j)
      do 15  i=ifp(j,l),ilp(j,l)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
 15   dpmn(j)=min(dpmn(j),dp(i,j,kn))
c$OMP END PARALLEL DO
c
      dpmin=999.
c$OMP PARALLEL DO REDUCTION(min:dpmin)
      do 149 j=1,jj
 149  dpmin=min(dpmin,dpmn(j))
c$OMP END PARALLEL DO
c
      if (dpmin.lt.-onem) then
      do 150 j=1,jj
      do 150 l=1,isp(j)
      do 150 i=ifp(j,l),ilp(j,l)
      if (dp(i,j,kn).eq.dpmin) write (lp,100) nstep,i,j,k,15,dpmin/onem
 150  continue
      end if
c
 76   continue
c
c --- restore nondivergence of vertically integrated mass flow by
c --- recovering fluxes lost in the flux limiting process.
c --- treat these fluxes as an 'upstream' barotropic correction to
c --- the sum of diffusive and antidiffusive fluxes obtained so far.
c
      do 77 k=1,kk
      km=k+mm
      kn=k+nn
c
c$OMP PARALLEL DO PRIVATE(ja,q)
      do 45 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 44 l=1,isu(j)
      do 44 i=ifu(j,l),ilu(j,l)
      if (utotn(i,j).ge.0.) then
        q=dp(i-1,j,kn)/p(i-1,j,kk+1)
      else
        q=dp(i  ,j,kn)/p(i  ,j,kk+1)
      end if
      uflux(i,j)=utotn(i,j)*q
 44   uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
c
      do 45 l=1,isv(j)
      do 45 i=ifv(j,l),ilv(j,l)
      if (vtotn(i,j).ge.0.) then
        q=dp(i,ja ,kn)/p(i,ja ,kk+1)
      else
        q=dp(i,j  ,kn)/p(i,j  ,kk+1)
      end if
      vflux(i,j)=vtotn(i,j)*q
 45   vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 14 j=1,jj
      jb=mod(j     ,jj)+1
      dpmn(j)=999.
      do 14 l=1,isp(j)
      do 14 i=ifp(j,l),ilp(j,l)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
 14   dpmn(j)=min(dpmn(j),dp(i,j,kn))
c$OMP END PARALLEL DO
c
      dpmin=999.
c$OMP PARALLEL DO REDUCTION(min:dpmin)
      do 139 j=1,jj
 139  dpmin=min(dpmin,dpmn(j))
c$OMP END PARALLEL DO
c
      if (dpmin.lt.-onem) then
      do 140 j=1,jj
      do 140 l=1,isp(j)
      do 140 i=ifp(j,l),ilp(j,l)
      if (dp(i,j,kn).eq.dpmin) write (lp,100) nstep,i,j,k,14,dpmin/onem
 140  continue
      end if
c
 77   continue
c
c --- add bottom-pressure restoring term arising from split-explicit treatment
c --- of continuity equation (step 4 in appendix B of 1992 BRHS paper)
c
c$OMP PARALLEL DO
      do 36 j=1,jj
      dpmn(j)=0.
      do 36 l=1,isp(j)
      do 36 i=ifp(j,l),ilp(j,l)
      util3(i,j)=abs(p(i,j,kk+1)-pbot(i,j))
 36   dpmn(j)=max(dpmn(j),util3(i,j))
c$OMP END PARALLEL DO
c
      dpmin=0.
c$OMP PARALLEL DO REDUCTION(max:dpmin)
      do 37 j=1,jj
 37   dpmin=max(dpmin,dpmn(j))
c$OMP END PARALLEL DO
c
      if (dpmin.gt.2.*onem) then
      do 38 j=1,jj
      do 38 l=1,isp(j)
      do 38 i=ifp(j,l),ilp(j,l)
      if (util3(i,j).eq.dpmin) write (lp,'(i9,2i5,a,f9.2,a)')
     .   nstep,i,j,'  warning - large pbot correction:',dpmin/onem,' m'
 38     continue
      end if
c
c$OMP PARALLEL DO PRIVATE(kn,old)
      do 39 j=1,jj
      do 39 l=1,isp(j)
      do 39 k=1,kk
      kn=k+nn
      do 39 i=ifp(j,l),ilp(j,l)
      old=dp(i,j,kn)
      dp(i,j,kn)=dp(i,j,kn)*pbot(i,j)/p(i,j,kk+1)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-old)                ! diapyc.flux
 39   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c$OMP END PARALLEL DO
c
c --- ---------------------------------------------
c --- laplacian interface smoothing (=> bolus flux)
c --- ---------------------------------------------
c
      if (thkdff.eq.0.) return
c
      if (nstep.eq.1 .or. diagno)
     .  q=hyc_pechg1(dp(1,1,k1n),th3d(1,1,k1n),31)
c
c --- diagnostic use only
c     if (nstep.eq.1 .or. diagno) then
c       call pechng1(dp(1,1,k1n),th3d(1,1,k1n))
c       call ape    (dp(1,1,k1n),th3d(1,1,k1n))
c     end if
c
      dtinv=1./delt1
      hymar=1./(10.*sigjmp)
c
      do 13 k=kk,2,-1
      km=k+mm
      kn=k+nn
c
c$OMP PARALLEL DO PRIVATE(flxhi,flxlo,ia,ib,ja,jb)
      do 151 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 141 l=1,isu(j)
      do 141 i=ifu(j,l),ilu(j,l)
      uflux(i,j)=delt1*thkdff*(p(i-1,j,k)-p(i,j,k))*scuy(i,j)
ccc     .   *boost(pbot(i,j),pbot(i-1,j),p(i,j,k),p(i-1,j,k))
c --- confine interface smoothing to isopycnic coord. subdomain
      uflux(i,j)=uflux(i,j)*min(1.,max(.1,
     .   2.-(max(th3d(i,j,km-1),th3d(i-1,j,km-1))-theta(k-1))*hymar))
ccc  .   2.-(max(th3d(i,j,km  ),th3d(i-1,j,km  ))-theta(k  ))*hymar))
c --- prevent intertwining of interfaces
      flxhi= .25*min((p(i  ,j,k+1)-p(i  ,j,k  ))*scp2(i  ,j),
     .               (p(i-1,j,k  )-p(i-1,j,k-1))*scp2(i-1,j))
      flxlo=-.25*min((p(i-1,j,k+1)-p(i-1,j,k  ))*scp2(i-1,j),
     .               (p(i  ,j,k  )-p(i  ,j,k-1))*scp2(i  ,j))
      uflux(i,j)=min(flxhi,max(flxlo,uflux(i,j)))
c --- add bolus flux to mass flux
      uflx(i,j,k-1)=uflx(i,j,k-1)+uflux(i,j)*dtinv
 141  uflx(i,j,k  )=uflx(i,j,k  )-uflux(i,j)*dtinv
c
      do 151 l=1,isv(j)
      do 151 i=ifv(j,l),ilv(j,l)
      vflux(i,j)=delt1*thkdff*(p(i,ja ,k)-p(i,j,k))*scvx(i,j)
ccc     .   *boost(pbot(i,j),pbot(i,ja ),p(i,j,k),p(i,ja ,k))
c --- confine interface smoothing to isopycnic coord. subdomain
      vflux(i,j)=vflux(i,j)*min(1.,max(.1,
     .   2.-(max(th3d(i,j,km-1),th3d(i,ja ,km-1))-theta(k-1))*hymar))
ccc  .   2.-(max(th3d(i,j,km  ),th3d(i,ja ,km  ))-theta(k  ))*hymar))
c --- prevent intertwining of interfaces
      flxhi= .25*min((p(i,j  ,k+1)-p(i,j  ,k  ))*scp2(i,j  ),
     .               (p(i,ja ,k  )-p(i,ja ,k-1))*scp2(i,ja ))
      flxlo=-.25*min((p(i,ja ,k+1)-p(i,ja ,k  ))*scp2(i,ja ),
     .               (p(i,j  ,k  )-p(i,j  ,k-1))*scp2(i,j  ))
      vflux(i,j)=min(flxhi,max(flxlo,vflux(i,j)))
c --- add bolus flux to mass flux
      vflx(i,j,k-1)=vflx(i,j,k-1)+vflux(i,j)*dtinv
 151  vflx(i,j,k  )=vflx(i,j,k  )-vflux(i,j)*dtinv
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 18 j=1,jj
      jb=mod(j     ,jj)+1
      do 18 l=1,isp(j)
      do 18 i=ifp(j,l),ilp(j,l)
      pold(i,j)=p(i,j,k)
 18   p(i,j,k)=p(i,j,k)-(uflux(i+1,j)-uflux(i,j)
     .                  +vflux(i,jb )-vflux(i,j))*scp2i(i,j)
c$OMP END PARALLEL DO
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,3i4,'' intfc.depth diffusion -- p_old,p_new ='',
cdiag.2f9.3)') nstep,itest,jtest,k,pold(itest,jtest)/onem,p(itest,
cdiag.jtest,k)/onem
c
 13   continue
c
c$OMP PARALLEL DO PRIVATE(kn)
      do 17 j=1,jj
      do 17 l=1,isp(j)
      do 17 k=1,kk
      kn=k+nn
      do 17 i=ifp(j,l),ilp(j,l)
ccc      p(i,j,k+1)=max(p(i,j,k),min(p(i,j,k+1),p(i,j,kk+1)))
      dp(i,j,kn)=p(i,j,k+1)-p(i,j,k)
 17   continue
c$OMP END PARALLEL DO
      if (nstep.eq.1 .or. diagno)
     . write (501,103) time,'  APE change due to intfc smoothing:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
 103  format (f9.1,a,-12p,f9.3,' TW')
c
c --- diagnostic use only
c     if (nstep.eq.1 .or. diagno) then
c       call ape    (dp(1,1,k1n),th3d(1,1,k1n))
c       call pechng2(dp(1,1,k1n),th3d(1,1,k1n))
c     end if
c
c --- to conserve tracer, project it from mid (mm) onto new (nn) thknss field
c
c     if (trcout) call recast(dp(1,1,1+mm),dp(1,1,1+nn))
c
      return
      end
c
c
c> Revision history:
c>
c> July 1997 - combined diff. and antidiff.flux calc. (eliminated loops 20,21)
c> Aug. 1997 - set u/vflux=0 before entering thickness smoothing k-loop 13
c> Jul. 1998 - reversed i/j loop nesting in loop 26
c> Apr. 2000 - changed i/j loop nesting to j/i
c> May  2000 - added code to eliminate neg. dp resulting from intfc.smoothing
c> May  2000 - added provisions to enhance -thkdff- in regions of strong flow
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Oct. 2000 - added limiter to intfc smoothing to maintain finite mxlyr thknss
c> Feb. 2001 - further enhanced -thkdff- enhancer by adding -glue-
c> June 2001 - added -dp- change in loop 39 to diapycnal flux diagnostics
c> Mar. 2002 - changed thickness smoothing from biharmonic to laplacian
c> Sep. 2003 - added diagnostics to track magnitude of -pbot- restoration
c> Sep. 2003 - added call to -recast- to improve tracer conservation
c> Dec. 2003 - confined interface smoothing to isopycnic coord. subdomain
c> Jan. 2004 - boosted thickness diffusion near sea floor (loops 141,151)
c> June 2004 - amended def'n of flxlo,flxhi in intfc.smoothing (loops 141,151)
