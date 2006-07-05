      subroutine barotp(m,n,mm,nn,k1m,k1n)
c
c --- version 2.8.2
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real q,utndcy,vtndcy,damp,uglue,vglue,
     .     u_ja,u_jb,v_ia,v_ib,dpu_ja,dpu_jb,dpv_ia,dpv_ib,hfharm
      external hfharm
      integer lll,ml,nl,mn,ll,kan,jcyc
      logical vthenu,mxing
      character text*20
      data damp/1.e-6/			!  newtonian damping
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations from baroclinic time level -m- to level -n-
c --- ------------------------------------------------------------------------
c
      ml=n
      nl=3
c
c --- explicit time integration of barotropic flow (forward-backward scheme)
c --- in order to combine forward-backward scheme with leapfrog treatment of
c --- coriolis term, v-eqn must be solved before u-eqn every other time step
      vthenu=.false.
c
      do 840 lll=1,lstep
      mxing=.false.
ccc      if (mod(4*lll,lstep).eq.0) mxing=.true.
      if (mod(2*lll,lstep).eq.0) mxing=.true.
ccc      mxing=.true.				!  mix every time step
c
c --- continuity equation
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 843 j=1,jj
      jb=mod(j     ,jj)+1
      do 843 l=1,isp(j)
      do 843 i=ifp(j,l),ilp(j,l)
 843  pbavg(i,j,nl)=(1.-wbaro)*pbavg(i,j,ml)+wbaro*pbavg(i,j,nl)
     . -(1.+wbaro)*dlt*(ubavg(i+1,j,ml)*depthu(i+1,j)*scuy(i+1,j)
     .                 -ubavg(i  ,j,ml)*depthu(i  ,j)*scuy(i  ,j)
     .                 +vbavg(i,jb ,ml)*depthv(i,jb )*scvx(i,jb )
     .                 -vbavg(i,j  ,ml)*depthv(i,j  )*scvx(i,j  ))
     .  *scp2i(i,j)
c$OMP END PARALLEL DO
c
      mn=ml
      if (vthenu) go to 901
c
c --- u momentum equation
c
 900  continue
c
      if (mxing) then
c$OMP PARALLEL DO PRIVATE(ja,jb,u_ja,u_jb,dpu_ja,dpu_jb)
      do 822 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
c
      do 824 l=1,isp(j)
      do 824 i=ifp(j,l),ilp(j,l)
c --- longitudinal turb. momentum flux (at mass points)
 824  if (iu(i,j)+iu(i+1,j).gt.0)
     . uflux1(i,j)=(ubavg(i,j,ml)-ubavg(i+1,j,ml))
     .             *hfharm(depthu(i,j),depthu(i+1,j))
     .             *scp2(i,j)*2./(scux(i,j)+scux(i+1,j))
c
      do 822 l=1,isu(j)
      do 822 i=ifu(j,l),ilu(j,l)
      if (depthu(i,ja ).gt.0.) then
        u_ja=ubavg(i,ja ,ml)
        dpu_ja=depthu(i,ja )
      else
        u_ja=slip*ubavg(i,j,ml)
        dpu_ja=depthu(i,j)
      end if
      if (depthu(i,jb ).gt.0.) then
        u_jb=ubavg(i,jb ,ml)
        dpu_jb=depthu(i,jb )
      else
        u_jb=slip*ubavg(i,j,ml)
        dpu_jb=depthu(i,j)
      end if
c --- lateral turb. momentum flux (at vorticity points)
      uflux2(i,j)=(u_ja-ubavg(i,j,ml))*hfharm(depthu(i,j),dpu_ja)
     .            *scq2(i,j )*2./(scuy(i,j)+scuy(i,ja))
      uflux3(i,j)=(ubavg(i,j,ml)-u_jb)*hfharm(depthu(i,j),dpu_jb)
     .            *scq2(i,jb)*2./(scuy(i,j)+scuy(i,jb))
 822  continue
c$OMP END PARALLEL DO
c
      call cpy_p(uflux1)
      end if				!  mxing = .true.
c
      call cpy_p(pbavg(1,1,nl))
c
c$OMP PARALLEL DO PRIVATE(jb,utndcy,uglue)
      do 841 j=1,jj
      jb=mod(j     ,jj)+1
      do 841 l=1,isu(j)
      do 841 i=ifu(j,l),ilu(j,l)
      utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)
     .+(vbavg(i  ,j,mn)*depthv(i  ,j)+vbavg(i  ,jb ,mn)*depthv(i  ,jb )
     . +vbavg(i-1,j,mn)*depthv(i-1,j)+vbavg(i-1,jb ,mn)*depthv(i-1,jb ))
     . *(pvtrop(i,j)+pvtrop(i,jb ))*.125
      util1(i,j)=(damp*(glue(i,j)+glue(i-1,j)-2.0)+dampu(i,j))
     . *ubavg(i,j,ml)
      if (mxing) util1(i,j)=util1(i,j)
     . +veldff*(uflux1(i,j)-uflux1(i-1,j)+uflux3(i,j)-uflux2(i,j))/
     .  (scu2(i,j)*depthu(i,j)) * (1.+10.*(glue(i,j)+glue(i-1,j)-2.0))
 841  ubavg(i,j,nl)=(1.-wbaro)*ubavg(i,j,ml)+wbaro*ubavg(i,j,nl)
     . +(1.+wbaro)*dlt*(utndcy+utotn(i,j)-util1(i,j))
c$OMP END PARALLEL DO
c
      if (abs(ubavg(ipacs,jpac,nl)+ubavg(iatln,jatl,nl)).gt.
     .    abs(ubavg(ipacs,jpac,nl)-ubavg(iatln,jatl,nl))*1.e-12)
     .  write(*,'(2i4,a,1p,2e15.7)') nstep,lll,' barotp WRONG ubavg_nl '
     .  ,ubavg(ipacs,jpac,nl),ubavg(iatln,jatl,nl)
c
cdiag write (lp,100) nstep
cdiag do jcyc=jtest,jtest+1
cdiag j =mod(jcyc-1+jj,jj)+1
cdiag jb=mod(jcyc     ,jj)+1
cdiag do i=itest,itest+1
cdiag if (iu(i,j).gt.0) then
cdiag write (lp,'(i3,2i5,2p,6f9.3)') lll,i,j,ubavg(i,j,ml),ubavg(i,j,nl)
cdiag. ,-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
cdiag. (vbavg(i  ,j,mn)*depthv(i  ,j)+vbavg(i  ,jb ,mn)*depthv(i  ,jb )
cdiag. +vbavg(i-1,j,mn)*depthv(i-1,j)+vbavg(i-1,jb ,mn)*depthv(i-1,jb ))
cdiag. *(pvtrop(i,j)+pvtrop(i,jb ))*.125*dlt,utotn(i,j)*dlt,
cdiag. -util1(i,j)*dlt
cdiag end if
cdiag end do
cdiag end do
 100  format(i9,8x,'ubold    ubnew    gradp    corio    ustar     fric')
c
      mn=nl
      if (vthenu) go to 902
c
c --- v momentum equation
c
 901  continue
c
      if (mxing) then
c$OMP PARALLEL DO PRIVATE(ia,ib,ja,jb,v_ia,v_ib,dpv_ia,dpv_ib)
      do 823 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
c
      do 825 l=1,isp(j)
      do 825 i=ifp(j,l),ilp(j,l)
c --- longitudinal turb. momentum flux (at mass points)
 825  if (iv(i,j)+iv(i,jb ).gt.0)
     . vflux1(i,j)=(vbavg(i,j,ml)-vbavg(i,jb ,ml))
     .             *hfharm(depthv(i,j),depthv(i,jb ))
     .             *scp2(i,j)*2./(scvy(i,j)+scvy(i,jb ))
c
      do 823 l=1,isv(j)
      do 823 i=ifv(j,l),ilv(j,l)
      ia=mod(i-2+ii,ii)+1
      ib=i+1
      if (depthv(ia ,j).gt.0.) then
        v_ia=vbavg(ia ,j,ml)
        dpv_ia=depthv(ia ,j)
      else
        v_ia=slip*vbavg(i,j,ml)
        dpv_ia=depthv(i,j)
      end if
      if (depthv(ib ,j).gt.0.) then
        v_ib=vbavg(ib ,j,ml)
        dpv_ib=depthv(ib ,j)
      else
        v_ib=slip*vbavg(i,j,ml)
        dpv_ib=depthv(i,j)
      end if
c --- lateral turb. momentum flux (at vorticity points)
      vflux2(i,j)=(v_ia-vbavg(i,j,ml))*hfharm(depthv(i,j),dpv_ia)
     .            *scq2(i ,j)*2./(scvx(i,j)+scvx(ia,j))
      vflux3(i,j)=(vbavg(i,j,ml)-v_ib)*hfharm(depthv(i,j),dpv_ib)
     .            *scq2(ib,j)*2./(scvx(i,j)+scvx(ib,j))
 823  continue
c$OMP END PARALLEL DO
      end if				!  mxing = .true.
c
c$OMP PARALLEL DO PRIVATE(ja,vtndcy,vglue)
      do 842 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 842 l=1,isv(j)
      do 842 i=ifv(j,l),ilv(j,l)
      vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,ja ,nl))*scvyi(i,j)
     .-(ubavg(i,j  ,mn)*depthu(i,j  )+ubavg(i+1,j  ,mn)*depthu(i+1,j  )
     . +ubavg(i,ja ,mn)*depthu(i,ja )+ubavg(i+1,ja ,mn)*depthu(i+1,ja ))
     . *(pvtrop(i,j)+pvtrop(i+1,j))*.125
      util2(i,j)=(damp*(glue(i,j)+glue(i,ja )-2.0)+dampv(i,j))
     . *vbavg(i,j,ml)
      if (mxing) util2(i,j)=util2(i,j)
     . +veldff*(vflux1(i,j)-vflux1(i,ja )+vflux3(i,j)-vflux2(i,j))/
     .  (scv2(i,j)*depthv(i,j)) * (1.+10.*(glue(i,j)+glue(i,ja )-2.0))
c
 842   vbavg(i,j,nl)=(1.-wbaro)*vbavg(i,j,ml)+wbaro*vbavg(i,j,nl)
     . +(1.+wbaro)*dlt*(vtndcy+vtotn(i,j)-util2(i,j))
c$OMP END PARALLEL DO
c
cdiag write (lp,101) nstep
cdiag do jcyc=jtest,jtest+1
cdiag ja=mod(jcyc-2+jj,jj)+1
cdiag j =mod(jcyc-1+jj,jj)+1
cdiag jb=mod(jcyc     ,jj)+1
cdiag do i=itest,itest+1
cdiag if (iv(i,j).gt.0) then
cdiag write (lp,'(i3,2i5,2p,6f9.3)') lll,i,j,vbavg(i,j,ml),vbavg(i,j,nl)
cdiag. ,-thref*(pbavg(i,j,nl)-pbavg(i,ja ,nl))*scvyi(i,j)*dlt,
cdiag.-(ubavg(i,j  ,mn)*depthu(i,j  )+ubavg(i+1,j  ,mn)*depthu(i+1,j  )
cdiag. +ubavg(i,ja ,mn)*depthu(i,ja )+ubavg(i+1,ja ,mn)*depthu(i+1,ja ))
cdiag. *(pvtrop(i,j)+pvtrop(i+1,j))*.125*dlt,vtotn(i,j)*dlt,
cdiag. -util2(i,j)*dlt
cdiag end if
cdiag end do
cdiag end do
 101  format(i9,8x,'vbold    vbnew    gradp    corio    ustar     fric')
c
      mn=nl
      if (vthenu) go to 900
c
c --- switch order in which -u,v- equations are solved
 902  vthenu=.not.vthenu
c
      ll=ml
      ml=nl
      nl=ll
c
 840  continue
c
      return
      end
c
c
c> Revision history:
c>
c> Mar. 1995 - changed vertical velocity averaging interval from 10 cm to 1 m
c>             (loops 33,35)
c> Mar. 1995 - changed order of loop nesting in loop 842
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - transferred loops preceding loop 840 to momeq2.f
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Oct. 2000 - added regional viscosity enhancement ('glue(i,j)')
c> Oct. 2002 - added lateral mixing terms to the momentum equations
