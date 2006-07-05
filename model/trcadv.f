      subroutine tradv0(n,nn)
c
c --- online tracer advection routine designed for intermittent
c --- (i.e., long time step) execution. there are 3 entries:
c
c --- tradv0 - initializes mass flux arrays and saves initial -dp- field
c ---          (should be called immediately  a f t e r  diapfl)
c --- tradv1 - builds up time integral of horizontal mass fluxes
c --- tradv2 - performs the actual transport operation
c ---          (should be called immediately  b e f o r e  diapfl)
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
c --- initialize arrays
c
c$OMP PARALLEL DO
      do 1 j=1,jj
      do 1 k=1,kk
c
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
 2    ufxcum(i,j,k)=0.
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
 3    vfxcum(i,j,k)=0.
c
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    dpinit(i,j,k)=dp(i,j,k+nn)
 1    continue
c$OMP END PARALLEL DO
      oddev=n
      write (lp,'(a)') 'tracer transport arrays initialized'
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tradv1(n,nn)
c
c --- build up time integrals of horiz. mass fluxes
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      if (n.ne.oddev) then
        write (lp,'(2(a,i2/))')
     .   'tracer advection time interval begins/ends at n =',oddev,
     .    '=> mass fluxes must be accumulated when n =',oddev
        stop '(n=oddev required in tradv1)'
      end if
c
c$OMP PARALLEL DO
      do 5 j=1,jj
      do 5 k=1,kk
c
      do 6 l=1,isu(j)
      do 6 i=ifu(j,l),ilu(j,l)
 6    ufxcum(i,j,k)=ufxcum(i,j,k)+uflx(i,j,k)*delt1

      do 7 l=1,isv(j)
      do 7 i=ifv(j,l),ilv(j,l)
 7    vfxcum(i,j,k)=vfxcum(i,j,k)+vflx(i,j,k)*delt1
 5    continue
c$OMP END PARALLEL DO
      write (lp,'(a)') 'mass fluxes saved for tracer transport'
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tradv2(n,nn)
c
c --- advect tracer over 'mixfrq' time steps
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real vertfx(idm,jdm,kdm),hordiv(idm,jdm,kdm),
     .     coldiv(idm),verdiv,q,fluxdv,thkchg,
     .     trcold(jdm,kdm),prold(jdm,kdm+1),
     .     trcnew(jdm,kdm),prnew(jdm,kdm+1)
      integer ka,nt
      character string*18
c
      if (n.ne.oddev) then
        write (lp,'(2(a,i2))')
     .   'tracer advection interval began at n =',oddev,
     .    '  and must end at n=',oddev
        stop '(n=oddev required in tradv2)'
      end if
c
c --- compute horizontal flux divergence (units: per transport time step)
c
      do k=1,kk
        ufxcum(iatls,jatl,k)=-ufxcum(ipacn,jpac,k)
      end do
c
c$OMP PARALLEL DO PRIVATE(ib,jb,q,coldiv)
      do 9 j=1,jj
      jb=mod(j,jj)+1
      do 9 l=1,isp(j)
c
      do 10 i=ifp(j,l),ilp(j,l)
 10   coldiv(i)=0.
c
      do 11 k=1,kk
      do 11 i=ifp(j,l),ilp(j,l)
      ib=mod(i,ii)+1
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
      hordiv(i,j,k)=(ufxcum(ib,j,k)-ufxcum(i,j,k)
     .              +vfxcum(i,jb,k)-vfxcum(i,j,k))*scp2i(i,j)
      coldiv(i)=coldiv(i)+hordiv(i,j,k)
c
cdiag if (i.eq.itest .and. j.eq.jtest) write (lp,103) i,j,k,
cdiag. 'mass flux',ufxcum(i,j,k)*scp2i(i,j),vfxcum(i,j,k)*scp2i(i,j),
cdiag.  coldiv(i),vfxcum(i,jb,k)*scp2i(i,j),ufxcum(ib,j,k)*scp2i(i,j)
 11   continue
 103  format (2i5,i3,2x,a/f21.2/f14.2,2f7.2/f21.2)
c
c --- adjust initial layer thickness to cancel effect of column divergence
c
      do 12 i=ifp(j,l),ilp(j,l)
      util1(i,j)=(p(i,j,kk+1)+coldiv(i))/p(i,j,kk+1)
 12   util2(i,j)=1./util1(i,j)
c
      do 13 k=1,kk
      do 13 i=ifp(j,l),ilp(j,l)
      tracer(i,j,k,:)=tracer(i,j,k,:)*util2(i,j)
 13   dpinit(i,j,k)=dpinit(i,j,k)*util1(i,j)
c
c --- compute the various terms in the continuity equation integrated
c --- over time interval since last call to -tradv0-
c --- the continuity eqn is split into horiz. and vert. terms as follows:
c ---        (dpfinl-dpinit) + hordiv + verdiv = 0
c
      do 15 i=ifp(j,l),ilp(j,l)
 15   vertfx(i,j,1)=0.
c
      do 16 k=1,kk
      ka=max(1,k-1)
      do 16 i=ifp(j,l),ilp(j,l)
      verdiv=(dpinit(i,j,k)-dp(i,j,k+nn))-hordiv(i,j,k)
 16   vertfx(i,j,k)=vertfx(i,j,ka)+verdiv	!  flx thru botm of lyr k
 9    continue
c$OMP END PARALLEL DO
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: check balance of terms in continuity eqn.
cdiag i=itest
cdiag j=jtest
cdiag ib=mod(i,ii)+1
cdiag jb=mod(j,jj)+1
cdiag write (lp,'(2i5,a/a)') i,j,
cdiag. '  trcadv -- time-integrated continuity eqn diagnostics:',
cdiag.  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
cdiag do k=1,kk
cdiag   thkchg=dp(i,j,k+nn)-dpinit(i,j,k)
cdiag   fluxdv=(ufxcum(ib,j,k)-ufxcum(i,j,k)
cdiag.         +vfxcum(i,jb,k)-vfxcum(i,j,k))*scp2i(i,j)
cdiag   if (k.eq.1) then
cdiag     write (lp,104) k,thkchg,fluxdv,vertfx(i,j,k),
cdiag.    thkchg+fluxdv+vertfx(i,j,k)
cdiag   else
cdiag     write (lp,104) k,thkchg,fluxdv,vertfx(i,j,k)-vertfx(i,j,k-1),
cdiag.    thkchg+fluxdv+vertfx(i,j,k)-vertfx(i,j,k-1)
cdiag   end if
cdiag end do
 104  format (i3,4f14.1)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
cdiag call totals(dpinit,tracer(1,1,1,1),
cdiag.            dpinit,tracer(1,1,1,2),'bfore fct3d')
c
      do nt=1,ntrcr
cdiag   do k=1,kk
cdiag     write (string,'(a,i2,a,i3)') 'trcr',nt,' bfor adv',k
cdiag     call findmx(ip,tracer(1,1,k,nt),idm,ii1,jj,string)
cdiag   end do
c
        call fct3d(2,tracer(1,1,1,nt),ufxcum,vfxcum,vertfx,
     .             scp2,scp2i,dpinit,dp(1,1,1+nn))
c
cdiag   do k=1,kk
cdiag     write (string,'(a,i2,a,i3)') 'trcr',nt,' aftr adv',k
cdiag     call findmx(ip,tracer(1,1,k,nt),idm,ii1,jj,string)
cdiag   end do
c
        do k=1,kk
          call cpy_p(tracer(1,1,k,nt))
        end do
      end do
c
cdiag call totals(dp(1,1,1+nn),tracer(1,1,1,1),
cdiag.            dp(1,1,1+nn),tracer(1,1,1,2),'after fct3d')
c
      write (lp,'(a)') 'tracer transport done'
c
      return
      end
c
c
c> Revision history:
c>
c> Dec. 2004 - made code cyclic in both -i- and -j- direction
c> Feb. 2005 - added multiple tracer capability
c> Mar. 2006 - added bering strait exchange logic
c
      subroutine fct3d(iord,fld,u,v,w,scal,scali,fco,fc)
c
c --- fully 3-d version of advfct.f
c
c  fld    - transported mixing ratio, e.g., salinity or temperature
c  u,v,w  - mass fluxes (x time step) satisfying continuity equation
c           (w(k) = mass flux per unit area across  b o t t o m  of layer k)
c  scal   - grid cell size
c  scali  - inverse of scal
c  fco,fc - depth of the layer at previous and new time step
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'bering.h'
c
      real fld(idm,jdm,kdm),u(idm,jdm,kdm),v(idm,jdm,kdm),
     .     w(idm,jdm,kdm),fco(idm,jdm,kdm),fc(idm,jdm,kdm),
     .     vertfx(idm,jdm,kdm),vertdv(idm,jdm,kdm),
     .     scal(idm,jdm),scali(idm,jdm)
      real fmx(idm,jdm),fmn(idm,jdm),flp(idm,jdm),fln(idm,jdm),
     .     flx(idm,jdm),fly(idm,jdm),uan(idm,jdm),van(idm,jdm),
     .     flxdiv(idm,jdm),clipj(jdm),vlumj(jdm)
      real a(kdm),b(kdm),c(kdm),athird,dx,fcdx,yl,yr
      real onemu,q,clip,vlume,amount,bfore,after,slab,dslab,thkchg,
     .     fluxdv,epsil,bforej(jdm),afterj(jdm)
      integer iord,ip1,im1,jp1,jm1,kp,jaa,itest,jtest
      character string*16
      logical wrap,recovr
      data recovr/.false./
      parameter (athird=1./3.)
      common/testpt/itest,jtest

c
c --- if iord=1, scheme reduces to simple donor cell scheme.
      parameter (epsil=1.e-11,onemu=1.e-6)
c
c --- optional: check mass conservation
cdiag i=itest
cdiag j=jtest
cdiag jb=mod(j,jj)+1
cdiag write (lp,'(2i5,a/a)') i,j,
cdiag. '  fct3d -- time-integrated continuity eqn diagnostics:',
cdiag.  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
cdiag do k=1,kk
cdiag thkchg=fc(i,j,k)-fco(i,j,k)
cdiag fluxdv=(u(i+1,j,k)-u(i,j,k)
cdiag.       +v(i,jb ,k)-v(i,j,k))*scali(i,j)
cdiag if (k.eq.1) then
cdiag   write (lp,103) k,thkchg,fluxdv,w(i,j,k),
cdiag.  thkchg+fluxdv+w(i,j,k)
cdiag else
cdiag   write (lp,103) k,thkchg,fluxdv,w(i,j,k)-w(i,j,k-1),
cdiag.  thkchg+fluxdv+w(i,j,k)-w(i,j,k-1)
cdiag end if
cdiag end do
 103  format (i3,4f14.1)
cc c
cc c$OMP PARALLEL DO
cc       do 6 j=1,jdm
cc       do 6 i=1,idm
cc  6    flxdiv(i,j)=0.
cc c$OMP END PARALLEL DO
cc       do 5 k=1,kk
cc c$OMP PARALLEL DO PRIVATE(jb)
cc       do 9 j=1,jj
cc       jb=mod(j,jj)+1
cc       do 9 l=1,isp(j)
cc       if (k.eq.1) then
cc         do 10 i=ifp(j,l),ilp(j,l)
cc  10     flxdiv(i,j)=(u(i+1,j,k)-u(i,j,k)+v(i,jb,k)-v(i,j,k))*scali(i,j)
cc      .    +w(i,j,k)           +fc(i,j,k)-fco(i,j,k)
cc       else
cc         do 12 i=ifp(j,l),ilp(j,l)
cc  12     flxdiv(i,j)=(u(i+1,j,k)-u(i,j,k)+v(i,jb,k)-v(i,j,k))*scali(i,j)
cc      .    +w(i,j,k)-w(i,j,k-1)+fc(i,j,k)-fco(i,j,k)
cc       end if
cc  9    continue
cc c$OMP END PARALLEL DO
cc       call findmx(ip,flxdiv,idm,ii1,jj,'mass consv')
cc       write (lp,*) 'shown below: mass consv. residual in layer',k
cc       call zebra(flxdiv,idm,idm-1,jdm)
cc  5    continue
c
c --- get vertical flux by summing -fld- over upstream slab of thickness -w-
c
c$OMP PARALLEL DO PRIVATE(jb,amount,slab,dslab,kp,a,b,c,dx,fcdx,yl,yr)
      do 26 j=1,jj
      jb=mod(j     ,jj)+1
      do 26 l=1,isp(j)
c
c --- fill massless cells with data from layer above or below
      do 17 k=kk-1,1,-1
      do 17 i=ifp(j,l),ilp(j,l)
 17   fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k+1)*onemu)
     .          /(           fco(i,j,k)+             onemu)
      do 18 k=2,kk
      do 18 i=ifp(j,l),ilp(j,l)
 18   fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k-1)*onemu)
     .          /(           fco(i,j,k)+             onemu)
c
      do 26 i=ifp(j,l),ilp(j,l)
c
c --- fit 0th, 1st, or 2nd deg. polynomial to tracer in each cell
      a(1 )=fld(i,j,1 )
      b(1 )=0.
      c(1 )=0.
      a(kk)=fld(i,j,kk)
      b(kk)=0.
      c(kk)=0.
      do 16 k=2,kk-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise constant method:
ccc      a(k)=fld(i,j,k)
ccc      b(k)=0.
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise linear method:
c --- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
ccc      a(k)=fld(i,j,k)
ccc      b(k)=0.
ccc      if (fld(i,j,k).le.min(fld(i,j,k-1),fld(i,j,k+1)) .or.
ccc     .    fld(i,j,k).ge.max(fld(i,j,k-1),fld(i,j,k+1))) then
ccc        b(k)=0.
ccc      else if ((fld(i,j,k+1)-fld(i,j,k-1))*(fld(i,j,k-1)+fld(i,j,k+1)
ccc     .  -2.*fld(i,j,k)).gt.0.) then
ccc        b(k)=fld(i,j,k)-fld(i,j,k-1)
ccc      else
ccc        b(k)=fld(i,j,k+1)-fld(i,j,k)
ccc      end if
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise parabolic method:
c --- fit parabola a+bx+cx^2 to tracer in each cell (-.5 < x < +.5)
      yl=.5*(fld(i,j,k-1)+fld(i,j,k))
      yr=.5*(fld(i,j,k+1)+fld(i,j,k))
      a(k)=1.5*fld(i,j,k)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-fld(i,j,k))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fld(i,j,k))) then
c --- apex of parabola occurs inside interval [-.5,+.5], implying an over-
c --- or undershoot situation. change curve to prevent over/undershoots.
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fld(i,j,k))) then
c --- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-fld(i,j,k)) .gt. 0.) then
c --- apex at x=-.5
            a(k)=.25*(3.*fld(i,j,k)+yl)
            c(k)=3.*(fld(i,j,k)-yl)
            b(k)=c(k)
          else
c --- apex at x=+.5
            a(k)=.25*(3.*fld(i,j,k)+yr)
            c(k)=3.*(fld(i,j,k)-yr)
            b(k)=-c(k)
          end if
        else			!  -1/6 < x < +1/6
c --- moving apex won't help. replace parabola by constant.
          a(k)=fld(i,j,k)
          b(k)=0.
          c(k)=0.
        end if
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 16   continue
c
      do 23 k=1,kk-1
      slab=onemu
      if (w(i,j,k).lt.0.) then			! interface moves down
        amount=slab*fld(i,j,k+1)
        kp=k
 24     kp=kp+1
        if (slab.ge.-w(i,j,k)) goto 23
        if (fco(i,j,kp).gt.0.) then
          dslab=min(slab+fco(i,j,kp),-w(i,j,k))
     .         -min(slab            ,-w(i,j,k))
          dx=dslab/fco(i,j,kp)
          fcdx=a(kp)
     .        +b(kp)*.5*(dx-1.)			!  not needed in pcm
     .        +c(kp)*(.25-dx*(.5-dx*athird))	!  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (kp.lt.kk) go to 24
      else if (w(i,j,k).gt.0.) then		! interface moves up
        amount=slab*fld(i,j,k)
        kp=k+1
 25     kp=kp-1
        if (slab.ge.w(i,j,k)) goto 23
        if (fco(i,j,kp).gt.0.) then
          dslab=min(slab+fco(i,j,kp), w(i,j,k))
     .         -min(slab            , w(i,j,k))
          dx=dslab/fco(i,j,kp)
          fcdx=a(kp)
     .        +b(kp)*.5*(1.-dx)			!  not needed in pcm
     .        +c(kp)*(.25-dx*(.5-dx*athird))	!  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (kp.gt.2) go to 25
      end if
 23   vertfx(i,j,k)=w(i,j,k)*amount/slab
c
      vertfx(i,j,kk)=0.			!  don't allow flux through bottom
      vertdv(i,j,1)=vertfx(i,j,1)
      do 26 k=2,kk
      vertdv(i,j,k)=vertfx(i,j,k)-vertfx(i,j,k-1)
c
 26   continue
c$OMP END PARALLEL DO
c
      bfore=0.
      after=0.
c
      do 4 k=1,kk
c
      if (beropn) then
        fld(ipacn,jpac,k)=0.
        fld(iatls,jatl,k)=0.
      end if
c
c$OMP PARALLEL DO SHARED(k)
      do 14 j=1,jj
      bforej(j)=0.
      do 14 l=1,isp(j)
      do 14 i=ifp(j,l),ilp(j,l)
 14   bforej(j)=bforej(j)+fld(i,j,k)*fco(i,j,k)*scal(i,j)
c$OMP END PARALLEL DO
c
c --- compute antidiffusive (high- minus low-order) fluxes
c
      call cpy_p(fld(1,1,k))
c
c$OMP PARALLEL DO PRIVATE(ja,jaa,jb,q,ia,ib) SHARED(k)
      do 11 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      jaa=mod(j-3+jj,jj)+1
c
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      if (u(i,j,k).ge.0.) then
        q=fld(i-1,j,k)
      else
        q=fld(i  ,j,k)
      end if
      flx(i,j)=u(i,j,k)*q
      q=fld(i,j,k)+fld(i-1,j,k)				!  2nd order
      if (ip(i+1,j)+iu(i-1,j).eq.2)
     .  q=1.125*q-.125*(fld(i+1,j,k)+fld(i-2,j,k))	!  4th order
 2    uan(i,j)=.5*q*u(i,j,k)-flx(i,j)
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
      if (v(i,j,k).ge.0.) then
        q=fld(i,ja ,k)
      else
        q=fld(i,j  ,k)
      end if
      fly(i,j)=v(i,j,k)*q
      q=fld(i,ja ,k)+fld(i,j,k)				!  2nd order
      if (ip(i,jb )+iv(i,ja).eq.2)
     .  q=1.125*q-.125*(fld(i,jb ,k)+fld(i,jaa,k))	!  4th order
 3    van(i,j)=.5*q*v(i,j,k)-fly(i,j)
c
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      ja=mod(j-2+jj,jj)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (ip(i,jb).eq.0) jb=j
      fmx(i,j)=max(fld(i,j,k),
     .   fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
      fmn(i,j)=min(fld(i,j,k),
     .   fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
      if (k.lt.kk) then
        if (w(i,j,k  ).lt.0.) then
          fmx(i,j)=max(fmx(i,j),vertfx(i,j,k  )/w(i,j,k  ))
          fmn(i,j)=min(fmn(i,j),vertfx(i,j,k  )/w(i,j,k  ))
        end if
      end if
      if (k.gt.1) then
        if (w(i,j,k-1).gt.0.) then
          fmx(i,j)=max(fmx(i,j),vertfx(i,j,k-1)/w(i,j,k-1))
          fmn(i,j)=min(fmn(i,j),vertfx(i,j,k-1)/w(i,j,k-1))
        end if
      end if
 11   continue
c$OMP END PARALLEL DO
c
cdiag write (string,'(a6,i2)') 'fmx k=',k
cdiag call findmx(ip,fmx,idm,ii1,jj,string(1:8))
cdiag write (string,'(a6,i2)') 'fmn k=',k
cdiag call findmx(ip,fmx,idm,ii1,jj,string(1:8))

c$OMP PARALLEL DO
      do 22 j=1,jj
      do 22 l=1,isp(j)
      flx(ifp(j,l)  ,j)=0.
      flx(ilp(j,l)+1,j)=0.
      uan(ifp(j,l)  ,j)=0.
      uan(ilp(j,l)+1,j)=0.
  22  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(j,wrap)
      do 33 i=1,ii1
      wrap=jfv(i,1).eq.1	! true if j=1 and j=jj are both water points
      do 33 l=1,jsp(i)
      j=jfp(i,l)
      if (j.gt.1 .or. .not.wrap) then
        fly(i,j)=0.
        van(i,j)=0.
      end if
      j=mod(jlp(i,l),jj)+1
      if (j.gt.1 .or. .not.wrap) then
        fly(i,j)=0.
        van(i,j)=0.
      end if
   33 continue
c$OMP END PARALLEL DO
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,101) 'advem(1)',i,j,k,fld(i-1,j,k),u(i,j,k),
cdiag. fld(i,j-1,k),v(i,j,k),fld(i,j,k),v(i,j+1,k),fld(i,j+1,k),
cdiag.  u(i+1,j,k),fld(i+1,j,k)
 101  format(a,2i5,i3,f18.3/1pe39.2/0pf19.3,1pe11.2,0pf9.3,
     .1pe11.2,0pf9.3/1pe39.2/0pf39.3)
c
c$OMP PARALLEL DO PRIVATE(jb,q,amount) SHARED(k)
      do 61 j=1,jj
      jb=mod(j     ,jj)+1
      if (recovr) vlumj(j)=0.
      if (recovr) clipj(j)=0.
      do 61 l=1,isp(j)
      do 61 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*scali(i,j)
c
cdiag if (i.eq.itest .and. j.eq.jtest)
cdiag. write (lp,'(2i5,i3,a,4f10.5,1pe9.2)') i,j,k,'  fc,fco,divs:',
cdiag.  fc(i,j,k),fco(i,j,k),flxdiv(i,j),vertdv(i,j,k),
cdiag.  fc(i,j,k)-fco(i,j,k)+flxdiv(i,j)+vertdv(i,j,k)
c
      q=fld(i,j,k)*fco(i,j,k)-flxdiv(i,j)-vertdv(i,j,k)
      amount=max(0.,fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
      if (recovr) then
        vlumj(j)=vlumj(j)+scal(i,j)*fc(i,j,k)
        clipj(j)=clipj(j)+(q-amount)*scal(i,j)
      end if
 61   fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
c$OMP END PARALLEL DO
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 61')
c
      if (iord.le.1) go to 100
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes
c
c$OMP PARALLEL DO PRIVATE(jb) SHARED(k)
      do 12 j=1,jj
      jb=mod(j     ,jj)+1
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      flp(i,j)=(fmx(i,j)-fld(i,j,k))*fc(i,j,k)
     ./((max(0.,uan(i,j))-min(0.,uan(i+1,j))
     .  +max(0.,van(i,j))-min(0.,van(i,jb ))+epsil)*scali(i,j))
c
 12   fln(i,j)=(fmn(i,j)-fld(i,j,k))*fc(i,j,k)
     ./((min(0.,uan(i,j))-max(0.,uan(i+1,j))
     .  +min(0.,van(i,j))-max(0.,van(i,jb ))-epsil)*scali(i,j))
c$OMP END PARALLEL DO
c
c---- limit antidiffusive fluxes
c
      call cpy_p(flp)
      call cpy_p(fln)
c
c$OMP PARALLEL DO PRIVATE(ja,clip)
      do 8 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 7 l=1,isu(j)
      do 7 i=ifu(j,l),ilu(j,l)
      if (uan(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i-1,j))
      else
        clip=min(1.,fln(i,j),flp(i-1,j))
      end if
 7    flx(i,j)=uan(i,j)*clip
c
      do 8 l=1,isv(j)
      do 8 i=ifv(j,l),ilv(j,l)
      if (van(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i,ja ))
      else
        clip=min(1.,fln(i,j),flp(i,ja ))
      end if
 8    fly(i,j)=van(i,j)*clip
c$OMP END PARALLEL DO
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,101) 'advem(2)',i,j,k,fld(i-1,j,k),u(i,j,k),
cdiag. fld(i,j-1,k),v(i,j,k),fld(i,j,k),v(i,j+1,k),fld(i,j+1,k),
cdiag.  u(i+1,j,k),fld(i+1,j,k)
c
c$OMP PARALLEL DO PRIVATE(jb,amount,q) SHARED(k)
      do 62 j=1,jj
      jb=mod(j     ,jj)+1
      do 62 l=1,isp(j)
      do 62 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*scali(i,j)
      q=fld(i,j,k)*fc(i,j,k)-flxdiv(i,j)
      amount=max(0.,fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
      if (recovr) clipj(j)=clipj(j)+(q-amount)*scal(i,j)
 62   fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
c$OMP END PARALLEL DO
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 62')
c
 100  continue
c
c --- recover 'clipped' amount and return to field layer by layer
c
      if (recovr) then
        vlume=0.
        clip=0.
c
        do 19 j=1,jj
        vlume=vlume+vlumj(j)
 19     clip=clip+clipj(j)
c
        if (vlume.ne.0.) then
          clip=clip/vlume
          write (lp,'(i2,a,1pe11.3)') k,'  tracer drift in fct3d',-clip
c$OMP PARALLEL DO SHARED(k)
          do 13 j=1,jj
          do 13 l=1,isp(j)
          do 13 i=ifp(j,l),ilp(j,l)
 13       fld(i,j,k)=fld(i,j,k)+clip
c$OMP END PARALLEL DO
        end if
      end if
c
      if (beropn) then
        fld(ipacn,jpac,k)=0.
        fld(iatls,jatl,k)=0.
      end if
c
c$OMP PARALLEL DO SHARED(k)
      do 15 j=1,jj
      afterj(j)=0.
      do 15 l=1,isp(j)
      do 15 i=ifp(j,l),ilp(j,l)
 15   afterj(j)=afterj(j)+fld(i,j,k)*fc(i,j,k)*scal(i,j)
c$OMP END PARALLEL DO
      do j=1,jj
        bfore=bfore+bforej(j)
        after=after+afterj(j)
      end do
c
 4    continue
c
      if (bfore.ne.0.)
     . write (lp,'(a,1p,3e14.6,e11.1)') 'fct3d conservation:',
     .  bfore,after,after-bfore,(after-bfore)/bfore
      q=1.
      if (after.ne.0.) q=bfore/after
      write (lp,'(a,f11.6)') 'fct3d: multiply tracer field by',q
      if (q.gt.1.1 .or. q.lt..9) stop '(excessive nonconservation)'
c
c$OMP PARALLEL DO
      do 20 j=1,jj
      do 20 k=1,kk
      do 20 l=1,isp(j)
      do 20 i=ifp(j,l),ilp(j,l)
 20   fld(i,j,k)=fld(i,j,k)*q
c$OMP END PARALLEL DO
c
      return
      end
c
c
c> Revision history:
c>
c> Feb. 2005 - added plm,ppm options to vertical flux calc'n (so far: pcm)
c> Mar. 2006 - added bering strait exchange logic
