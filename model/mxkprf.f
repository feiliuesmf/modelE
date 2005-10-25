      subroutine mxkprf(m,n,mm,nn,k1m,k1n)
ccc   use mod_xc    ! HYCOM communication interface
ccc   use mod_pipe  ! HYCOM debugging interface
c
c --- hycom version 2.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
      include 'kprf_arrays.h'
c
c ------------------------------------
c --- k-profile vertical mixing model
c --- (large, mc williams, doney)
c ------------------------------------
c
      integer ktr
      real    delp,sigmlj,thtop,thbot,thjmp(kdm),sigocn,dsigdt,dsigds
      external sigocn,dsigdt,dsigds
      character text*12
c
      include 'state_eqn.h'
c
 108  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f8.2,f8.1))
 109  format (i9,2i5,a/(33x,i3,f8.3,f8.1,3x,f8.3,f8.1))
 110  format (i9,2i5,a/(33x,i3,0pf8.3,1p,2e11.3))
      i=itest
      j=jtest
cdiag write (lp,108) nstep,itest,jtest,
cdiag.'  entering mxkprf:  temp    saln    dens    thkns    dpth',(k,
cdiag. temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn)+thbase,
cdiag.  dp(i,j,k+nn)/onem,p(i,j,k+1)/onem,k=1,kk)
cdiag if (dotrcr) write (lp,110) nstep,itest,jtest,
cdiag.'  entering mxkprf:  thkns   tracer1    tracer2    tracer3',
cdiag. (k,dp(i,j,k+nn)/onem,(tracer(i,j,k,ktr),ktr=1,2),k=1,kk)
cdiag write (lp,109) nstep,itest,jtest,
cdiag.'  entering mxkprf:    u     dp_u         v     dp_v',(k,
cdiag. u(i,j,k+nn),dpu(i,j,k+nn)/onem,v(i,j,k+nn),dpv(i,j,k+nn)/onem,
cdiag.  k=1,kk)
c
c --- diffusivity/viscosity calculation
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n,mm,nn,k1m,k1n)
cc !$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxkprfaj(m,n,mm,nn,k1m,k1n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- optional spatial smoothing of viscosity and diffusivities on interior
c --- interfaces.
c
ccc      if(difsmo) then
ccc        do k=2,kk
ccc          call psmooth_max(vcty(1,1,k), 1)
ccc          call psmooth_max(dift(1,1,k), 0)
ccc          call psmooth_max(difs(1,1,k), 0)
ccc        enddo
ccc      endif
c
c ---   final mixing of variables at p points
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n,mm,nn,k1m,k1n)
cc !$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxkprfbj(m,n,mm,nn,k1m,k1n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- final velocity mixing at u,v points
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n,mm,nn,k1m,k1n)
cc !$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxkprfcj(m,n,mm,nn,k1m,k1n, j)
      enddo
!$OMP END PARALLEL DO
c
      if (dotrcr) write (lp,'(a)') 'tracer kpp mixing done'
c
c --- mixed layer diagnostics
c
c --- diagnose new mixed layer depth based on density jump criterion
!$OMP   PARALLEL DO PRIVATE(kn,sigmlj,thtop,thbot,thjmp)
cc !$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
c
c --- depth of mixed layer base set to interpolated depth where
c --- the density jump is equivalent to a tmljmp temperature jump.
c --- this may not vectorize, but is used infrequently.
            do i=ifp(j,l),ilp(j,l)
              sigmlj = -tmljmp*dsigdt(temp(i,j,k1n),saln(i,j,k1n))
              dpmixl(i,j) = -zgrid(i,j,klist(i,j)+1)*onem  !bottom
              thjmp(1) = 0.0
              do k=2,klist(i,j)
                kn=k+nn
                thjmp(k) = max(th3d(i,j,kn)-th3d(i,j,k1n),
     &                         thjmp(k-1)) !stable profile simplifies the code
c
                if (thjmp(k).ge.sigmlj) then
c
c ---             find the two densities on the interface between layers
c ---             k-1 and k, using PLM assuming layers k-2 and k+1 are PCM.
c
                  if     (k.eq.2) then
                    thtop = thjmp(1)
                  else
                    thtop = thjmp(k-1) + 
     &                        min(thjmp(k-1)-thjmp(k-2),
     &                            thjmp(k)  -thjmp(k-1) )
                  endif !k.eq.2:else
                  if     (k.eq.klist(i,j)) then
                    thbot = thjmp(k)
                  else
                    thjmp(k+1) = max(th3d(i,j,kn+1)-th3d(i,j,k1n),
     &                               thjmp(k))
                    thbot = thjmp(k) - 
     &                        min(thjmp(k+1)-thjmp(k),
     &                            thjmp(k)  -thjmp(k-1) )
                  endif !k.eq.klist:else
                  if     (thtop.gt.thbot) then
                    thtop = 0.5*(thtop+thbot)  !make stable at interface
                    thbot = thtop
                  endif
c
                  if      (thtop.gt.sigmlj) then
c
c ---               in bottom half of layer k-1
c
                    dpmixl(i,j) =
     &                -zgrid(i,j,k-1)*onem +
     &                         0.5*dp(i,j,kn-1)*
     &                         (sigmlj+epsil-thjmp(k-1))/
     &                         (thtop +epsil-thjmp(k-1))
                  elseif (thbot.ge.sigmlj) then
c
c ---               at layer interface
c
                    dpmixl(i,j) =
     &                -zgrid(i,j,k-1)*onem + 0.5*dp(i,j,kn-1)
                  else
c
c ---               in top half of layer k
c
                    dpmixl(i,j) =
     &                -zgrid(i,j,k)*onem -
     &                         0.5*dp(i,j,kn)*
     &                         (1.0-(sigmlj  +epsil-thbot)/
     &                              (thjmp(k)+epsil-thbot) )
                  endif !part of layer
                  exit  !calculated dpmixl
                endif  !found dpmixl layer
              enddo
              dpmxav(i,j)=dpmxav(i,j)+dpmixl(i,j)
              oiceav(i,j)=oiceav(i,j)+oice(i,j)
            enddo
          enddo
        enddo
c
!$OMP   END PARALLEL DO
c
      if (diagno) then
c
c --- calculate bulk mixed layer t, s, theta
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,delp)
cc !$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
c
            do i=ifp(j,l),ilp(j,l)
              tmix(i,j)=temp(i,j,k1n)*dp(i,j,k1n)
              smix(i,j)=saln(i,j,k1n)*dp(i,j,k1n)
            enddo
c
            do k=2,kk
              do i=ifp(j,l),ilp(j,l)
                delp=min(p(i,j,k+1),dpmixl(i,j))
     &              -min(p(i,j,k  ),dpmixl(i,j))
                tmix(i,j)=tmix(i,j)+delp*temp(i,j,k+nn)
                smix(i,j)=smix(i,j)+delp*saln(i,j,k+nn)
              enddo
            enddo
c
            do i=ifp(j,l),ilp(j,l)
              tmix(i,j)=tmix(i,j)/dpmixl(i,j)
              smix(i,j)=smix(i,j)/dpmixl(i,j)
              thmix(i,j)=sigocn(tmix(i,j),smix(i,j))-thbase
           enddo
c
          enddo
        enddo
c
!$OMP   END PARALLEL DO
c
c --- calculate bulk mixed layer u
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,delp)
cc !$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isu(j)
c
            do i=ifu(j,l),ilu(j,l)
              umix(i,j)=u(i,j,k1n)*2.*dpu(i,j,k1n)
            enddo
c
            do k=2,kk
              do i=ifu(j,l),ilu(j,l)
                delp=
     &             (min(p(i,j,k+1)+p(i-1,j,k+1),
     &                  dpmixl(i,j)+dpmixl(i-1,j))
     &             -min(p(i,j,k  )+p(i-1,j,k  ),
     &                  dpmixl(i,j)+dpmixl(i-1,j)))
                umix(i,j)=umix(i,j)+delp*u(i,j,k+nn)
              enddo
            enddo
c
            do i=ifu(j,l),ilu(j,l)
              umix(i,j)=umix(i,j)/(dpmixl(i,j)+dpmixl(i-1,j))
            enddo
c
          enddo
        enddo
!$OMP   END PARALLEL DO
c
c --- calculate bulk mixed layer v
c
!$OMP   PARALLEL DO PRIVATE(j,ja,l,i,k,delp)
cc !$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          ja=mod(j-2+jj,jj)+1
          do l=1,isv(j)
c
            do i=ifv(j,l),ilv(j,l)
              vmix(i,j)=v(i,j,k1n)*2.*dpv(i,j,k1n)
            enddo
c
            do k=2,kk
              do i=ifv(j,l),ilv(j,l)
                delp=
     &             (min(p(i,j,k+1)+p(i,ja ,k+1),
     &                  dpmixl(i,j)+dpmixl(i,ja ))
     &             -min(p(i,j,k  )+p(i,ja ,k  ),
     &                  dpmixl(i,j)+dpmixl(i,ja )))
                vmix(i,j)=vmix(i,j)+delp*v(i,j,k+nn)
              enddo
            enddo
c
            do i=ifv(j,l),ilv(j,l)
              vmix(i,j)=vmix(i,j)/(dpmixl(i,j)+dpmixl(i,ja ))
            enddo
c
          enddo
        enddo
!$OMP   END PARALLEL DO
      endif                                           ! diagno
c
      i=itest
      j=jtest
cdiag write (lp,108) nstep,itest,jtest,
cdiag.'   exiting mxkprf:  temp    saln    dens    thkns    dpth',(k,
cdiag. temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn)+thbase,
cdiag.  dp(i,j,k+nn)/onem,p(i,j,k+1)/onem,k=1,kk)
cdiag if (dotrcr) write (lp,110) nstep,itest,jtest,
cdiag.'   exiting mxkprf:  thkns   tracer1    tracer2    tracer3',
cdiag. (k,dp(i,j,k+nn)/onem,(tracer(i,j,k,ktr),ktr=1,2),k=1,kk)
cdiag write (lp,109) nstep,itest,jtest,
cdiag.'   exiting mxkprf:    u     dp_u         v     dp_v',(k,
cdiag. u(i,j,k+nn),dpu(i,j,k+nn)/onem,v(i,j,k+nn),dpv(i,j,k+nn)/onem,
cdiag.  k=1,kk)
c
      return
      end
c
c***************************************************************************
c
      subroutine mxkprfaj(m,n,mm,nn,k1m,k1n, j)
ccc   use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
c
c --- calculate viscosity and diffusivity
c
      do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          if (mxlkpp) then
            call mxkppaij(m,n,mm,nn,k1m,k1n, i,j)
          else if (mxlgis) then
            call mxgisaij(m,n,mm,nn,k1m,k1n, i,j)
          else
            print *,'error - either mxlkpp or mxlgis must be true'
            stop '(mxlkpp,mxlgis)'
          end if
        enddo
      enddo
c
      return
      end
c
c***************************************************************************
      subroutine mxgisaij(m,n,mm,nn,k1m,k1n, i,j)
ccc   use mod_xc  ! HYCOM communication interface
c
c --- hycom version 2.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
      include 'kprf_arrays.h'
c
c--------------------------------------------------------------------
c --- nasa giss vertical mixing model, single j-row (part A)
c
c --- V.M. Canuto, A. Howard, Y. Cheng, and M.S. Dubovikov, 2001:
c --- Ocean turbulence, part I: One-point closure model -- momentum
c --- and heat vertical diffusivities.  JPO, 31, 1413-1426.
c --- V.M. Canuto, A. Howard, Y. Cheng, and M.S. Dubovikov, 2002:
c --- Ocean turbulence, part II: Vertical diffusivities of momentum,
c --- heat, salt, and passive tracers.  JPO, 32, 240-264.
c
c --- Modified by Armando Howard to implement the latitude dependent
c --- background mixing due to waves formula from:
c --- Gregg et. al. (2003): Reduced mixing from the breaking of
c --- internal waves in equatorial waters, Nature 422 pp 513-515.
c--------------------------------------------------------------------
c
c     1D turbulence calculation routine adapted by A.Romanou from A.Howard
c     from the 2000 original model.
c
c     For a discussion of the turbulence model used here see "Ocean
c     Turbulence. Part II" referenced above.
c
c     In general ria(k),rid(k) are calculated from the difference between
c     level k and k+1 and ak{m,h,s}(k) should be used to mix levels k and k+1.
c     n is (nlayers-1) because there are nlayers-1 ocean interfaces to
c     be mixed.
c
c     In the mixed layer the model diffusivity for each field is a product of 
c     a dimensionless function of the two variables ria and rid 
c     and the Shear and the square of a lengthscale.
c     The lengthscale is proportional to depth near the surface but asymptotes
c     towards a fixed fraction of the Mixed Layer Depth deeper in the mixed
c     layer. The MLD is thus a necessary ingredient for calculating model
c     diffusivities.
c
c     latitude dependent background mixing ("latdiw") option:
c       background mixing depends on |f| and N,
c       where f is the coriolis parameter and N brunt vaisala frequency.
c       note Gregg et al. use "f" when they mean the absolute value of
c       the coriolis parameter.
c     latitude dependent deep background mixing ("botdiw") option:
c       additional factor multiplying the `epsilon/N^2' for deep mixing
c       based on the formula cited as from Henyey et. al,
c       JGR vol.91 8487-8495,1986) in Gregg et al. where it is shown
c       confirmed by observations for lower latitudes except for being
c       low very near the equator. 
c       I place a minimum, "eplatidepmin", on the Gregg et al. factor, "L".
c       Note that Gregg et. al.'s formula:
c         L(\theta,N) = 
c         (|f| cosh^{-1} (N/|f|))/(f_30^o cosh^{-1} (N_0/f_30^o)
c       is only defined as a real number when N > |f|, since arccosh
c       can only be defined as a real for arguments of at least 1.
c       At 1 arccosh is zero.
c       I decide to set "L(\theta,N)" to "eplatidepmin"FOR (N/f < 1).
c       This corresponds to setting a floor of 1 on (N/f).
c       for foreground mixing at depth detached from the mixed-layer
c       I revert to the "deep" lengthscale, which uses density gradients,
c       in case (N/f)<1 to try not to make deep arctic&subarctic mixing 
c       too small.
c
c-----------------------------------------------------------------------------
c --- this is a level 2 turbulence model
c --- sm and sh depends only on the richardson number
c --- In salinity model case level 2 means S_M,S_H,S_C depend only on Ria,Ri_d
c-----------------------------------------------------------------------------
c                                                                             
      real, parameter :: difmax = 9999.0e-4  !maximum diffusion/viscosity     
      real, parameter :: acormin= 2.5453e-6  !minimum abs(corio), i.e. 1 degN     
c
c --- local variables for giss mixing
c
      real z1d(kdm),th1d(kdm),u1d(kdm),v1d(kdm)
      real ria(kdm),rid(kdm),s2(kdm),v_back(kdm),
     &     t_back(kdm),s_back(kdm),dtemp,dsaln,sflux1,
     &     alfadt,betads,al0,ri1,rid1,slq2,sm,sh,ss,akz,al,al2,anlq2,
     &     back_ri1,back_rid1,back_ra_r1,back_rit1,back_ric1,rit,ric,
     &     ra_r,theta_r,theta_r0,theta_r1,theta_r_deg,deltheta_r1,
     &     delback_ra_r,dback_ra_r_o_dtheta,slq2_back,sm_back,sh_back,
     &     ss_back,delsm_back,dsm_back_o_dtheta,delsh_back,
     &     dsh_back_o_dtheta,delss_back,dss_back_o_dtheta,delslq2_back,
     &     dslq2_back_o_dtheta,s2_back,al0_back,al_back,
     &     al2_back,anlq2_back,tmp_back,tmp,delz,delth,del2th,
     &     dzth,d2zth,rdzlndzth,al0deep,thsum,dens,
     &     beta_b,beta_r,frac_b,frac_r,qspcifh,hbl,
     &     epson2,epson2_bot,eplatidep
      real akm(kdm),akh(kdm),aks(kdm),aldeep(kdm),tmpk(kdm)
      real an2(kdm),an,acorio,zbot
c
      real swfrac(kdm+1)       ! fractional surface shortwave radiation flux
      real swfrml              ! fractional surface sw rad flux at ml base
      real hwide(kdm)          ! layer thicknesses in m (minimum 1mm)
      real dpmm(kdm)           !     max(onemm,dp(i,j,:,n))
      real qdpmm(kdm)          ! 1.0/max(onemm,dp(i,j,:,n))
      real pij(kdm+1)          ! local copy of p(i,j,:)
c
      integer ifbelow,ifrafglt,jtheta_r
      integer ifnofsmall
c
      integer            jtheta_r0,jtheta_r1,itheta_r0,itheta_r1
c     common/mxgissij_b/ jtheta_r0,jtheta_r1,itheta_r0,itheta_r1
c     save  /mxgissij_b/ !bugfix, reduces optimization of *theta_r*
c
      integer k1,jrlv
c
      integer iglobal,jglobal
c
      real      sigocn,dsigdt,dsigds,sigloc,dsiglocdt,dsiglocds
      external  sigocn,dsigdt,dsigds,sigloc,dsiglocdt,dsiglocds
      include 'state_eqn.h'
      real    acosh1,xx
c
      acosh1(xx) = log(xx+sqrt((xx**2)-1.0))
c
      jb=mod(j,jj)+1
c
c --- set mid-time pressure array
c --- locate lowest substantial mass-containing layer.
      dpmm( 1)=max(onemm,dp(i,j,k1n))
      qdpmm(1)=1.0/dpmm(1)
      pij(  1)=p(i,j,1)
      pij(  2)=pij(1)+dp(i,j,k1n)
      p(i,j,2)=pij(2)
      do k=2,kk
        kn=k+nn
        dpmm( k)  =max(onemm,dp(i,j,kn))
        qdpmm(k)  =1.0/dpmm(k)
        pij(  k+1)=pij(k)+dp(i,j,kn)
        p(i,j,k+1)=pij(k+1)
      enddo
      do k=kk,1,-1
        if (dpmm(k).gt.tencm) then
          exit
        endif
      enddo
      klist(i,j)=max(k,2)  !always consider at least 2 layers
c
c --- nominal surface bld is the mld, note that this depends on tmljmp
css   hbl=0.5*(dpmixl(i,j,n)+dpmixl(i,j,m))  !reduce time splitting
      hbl=dpmixl(i,j)
c
c --- forcing of t,s by surface fluxes. flux positive into ocean.
c --- shortwave flux penetration depends on kpar or jerlov water type.
c
      if     (jerlv0.eq.0) then
        beta_r = 2.0/onem
        beta_b = ( akpar(i,j,ls0)*ws0+akpar(i,j,ls1)*ws1
     &            +akpar(i,j,ls2)*ws2+akpar(i,j,ls3)*ws3)/onem
        frac_b = max( 0.27, 0.695 - 5.7*onem*beta_b )
        frac_r = 1.0 - frac_b
      else
        jrlv   = jerlov(i,j)
        beta_r = betard(jrlv)
        beta_b = betabl(jrlv)
        frac_r = redfac(jrlv)
        frac_b = 1.0 - frac_r
      endif
      qspcifh=1.0/spcifh
c
      do k=1,kk
        kn=k+nn
c
        if (thermo) then
          if     (-pij(k+1)*beta_r.gt.-10.0) then
            swfrac(k+1)=frac_r*exp(-pij(k+1)*beta_r)+
     &                  frac_b*exp(-pij(k+1)*beta_b)
          elseif (-pij(k+1)*beta_b.gt.-10.0) then
            swfrac(k+1)=frac_b*exp(-pij(k+1)*beta_b)
          else
            swfrac(k+1)=0.0
          endif
          if (k.eq.1) then
            sflux1=surflx(i,j)-sswflx(i,j)
            dtemp=(sflux1+(1.-swfrac(k+1))*sswflx(i,j))*
     &            delt1*g*qspcifh*qdpmm(k)
            dsaln=salflx(i,j)*
     &            delt1*g*        qdpmm(k)
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,101) nstep,i+i0,j+j0,k,
cdiag&          0.,1.-swfrac(k+1),dtemp,dsaln
cdiag         call flush(lp)
cdiag       endif
          elseif (k.le.klist(i,j)) then
            if(k.ne.klist(i,j)) then
              dtemp=(swfrac(k)-swfrac(k+1))*sswflx(i,j)*
     &              delt1*g*qspcifh*qdpmm(k)
            else
              dtemp= swfrac(k)             *sswflx(i,j)*
     &               delt1*g*qspcifh*qdpmm(k)
            endif
            dsaln=0.0
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,101) nstep,i+i0,j+j0,k,
cdiag&          1.-swfrac(k),1.-swfrac(k+1),dtemp
cdiag         call flush(lp)
cdiag       endif
          else !k.gt.klist(i,j)
            dtemp=0.0
            dsaln=0.0
          endif
        else !.not.thermo
          dtemp=0.0
          dsaln=0.0
        endif
c
c --- modify t and s
        temp(i,j,kn)=temp(i,j,kn)+dtemp
        saln(i,j,kn)=saln(i,j,kn)+dsaln
        th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))-thbase
c
      enddo
c
c --- calculate swfrml, the fraction of solar radiation left at depth hbl
      if     (-hbl*beta_r.gt.-10.0) then
        swfrml=frac_r*exp(-hbl*beta_r)+
     &         frac_b*exp(-hbl*beta_b)
      elseif (-hbl*beta_b.gt.-10.0) then
        swfrml=frac_b*exp(-hbl*beta_b)
      else
        swfrml=0.0
      endif
      mixflx(i,j)=surflx(i,j)-swfrml*sswflx(i,j)
c
c --- save array dpbl=hbl for ice, output and diagnosis
      dpbl(i,j)=hbl
c
c --- calculate z at vertical grid levels - this array is the z values in m
c --- at the mid-depth of each model layer except for index klist+1, where it
c --- is the z value of the bottom
c
c --- calculate layer thicknesses
      do k=1,kk
        kn=k+nn
        if (k.eq.1) then
          hwide(k)=dpmm(k)/onem
          zgrid(i,j,k)=-.5*hwide(k)
        else if (k.lt.klist(i,j)) then
          hwide(k)=dpmm(k)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
        else if (k.eq.klist(i,j)) then
          hwide(k)=dpmm(k)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
          zgrid(i,j,k+1)=zgrid(i,j,k)-.5*hwide(k)
        else
          hwide(k)=0.
        endif
c
c --- set 1-d array values; use cgs units
        if (k.le.klist(i,j)) then
          z1d (k)=-100.0*zgrid(i,j,k)
          u1d (k)=50.0*(u(i,j,kn)+u(i+1,j  ,kn))
          v1d (k)=50.0*(v(i,j,kn)+v(i  ,jb ,kn))
          if (k.eq.1) then
            thsum=th3d(i,j,kn)
            th1d(k)=1.0+0.001*(thsum+thbase)
          else
            deltaz=z1d(k)-z1d(k-1)
            s2(k-1)=((u1d(k-1)-u1d(k))**2+(v1d(k-1)-v1d(k))**2)/
     &              (deltaz*deltaz)
            if (locsig) then
              alfadt=0.0005*
     &           (dsiglocdt(temp(i,j,kn-1),saln(i,j,kn-1),p(i,j,k))+
     &            dsiglocdt(temp(i,j,kn  ),saln(i,j,kn  ),p(i,j,k)))*
     &           (temp(i,j,kn-1)-temp(i,j,kn))
              betads=0.0005*
     &           (dsiglocds(temp(i,j,kn-1),saln(i,j,kn-1),p(i,j,k))+
     &            dsiglocds(temp(i,j,kn  ),saln(i,j,kn  ),p(i,j,k)))*
     &           (saln(i,j,kn-1)-saln(i,j,kn))
              thsum=thsum-alfadt-betads
              th1d(k)=1.0+0.001*(thsum+thbase)
            else
              alfadt=0.0005*
     &           (dsigdt(temp(i,j,kn-1),saln(i,j,kn-1))+
     &            dsigdt(temp(i,j,kn  ),saln(i,j,kn  )))*
     &           (temp(i,j,kn-1)-temp(i,j,kn))
              betads=0.0005*
     &           (dsigds(temp(i,j,kn-1),saln(i,j,kn-1))+
     &            dsigds(temp(i,j,kn  ),saln(i,j,kn  )))*
     &           (saln(i,j,kn-1)-saln(i,j,kn))
              th1d(k)=1.0+0.001*(th3d(i,j,kn)+thbase)
            endif
            dens=1.0+0.001*(sigloc(temp(i,j,kn),saln(i,j,kn),
     &           p(i,j,k)+0.5*dp(i,j,kn))+thbase)
            an2(k-1)=-980.0*(alfadt+betads)/
     &             (deltaz*dens)
            ria(k-1)=-980.0*(alfadt+betads)/
     &             (deltaz*dens*max(epsil,s2(k-1)))
            rid(k-1)=-980.0*(alfadt-betads)/
     &             (deltaz*dens*max(epsil,s2(k-1)))
          endif
        endif
      enddo
      k=klist(i,j)
      s2 (k)=0.0
      ria(k)=0.0
      rid(k)=0.0
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   do k=1,klist(i,j)
cdiag   write(6,'(a,i9,i3,2f10.3,f9.1,f8.5,2f8.2)') 'giss1din1',
cdiag&        nstep,k,zgrid(i,j,k),hwide(k),z1d(k),
cdiag&        th1d(k),u1d(k),v1d(k)
cdiag   enddo
cdiag   write(6,'(a,a9,a3,3a13)') 
cdiag&    'giss1din2','    nstep','  k',
cdiag&    '           s2','          ria','          rid'
cdiag   do k=1,klist(i,j)
cdiag   write(6,'(a,i9,i3,1p,3e13.5)') 'giss1din2',
cdiag&        nstep,k,s2(k),ria(k),rid(k)
cdiag   enddo
cdiag endif
c
      al0=0.17*hbl/onecm
c
c --- Write internal turbulence quantities to fort.91 when writing enabled.
c ---  Headers for each outputstep.
c ---  Add S_M,H,S to outputs.
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write(lp,'(a,i9)') 'nstep = ',nstep
cdiag   write(lp,*) 'hbl,al0 = ',hbl/onecm,al0
cdiag   write(lp,*) 'b1 = ',b1
cdiag   write(lp,'(a)')    "  z          al         slq2       "//
cdiag&                     "ri1        rid1       "//
cdiag&                     "sm         sh         ss         "//
cdiag&                     "v_back     t_back     s_back     "
cdiag endif
c
      if     (latdiw) then
        acorio = max(acormin, !f(1degN)
     &               0.25*(abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     &                     abs(corio(i,jb ))+abs(corio(i+1,jb )) ))
      endif !latdiw
      if     (botdiw) then
        zbot = -100.0*zgrid(i,j,klist(i,j)+1)  !sea depth
      endif !botdiw
c
c --- START OF FIRST LOOP THROUGH LEVELS
c
      if (ifepson2.eq.2) then
c ---   Initialize switch for sub(background-only) depth. 
        ifbelow=0      
      endif
c
c --- depth-grid dooloop starts here
      do 22 k=1,klist(i,j)-1
        ri1=ria(k)
c
c --- Use Ri_d = Ri_C - Ri_T in salinity-temperature turbulence model.
        rid1=rid(k)
c
c --- Interpolate 2D table for salinity-temperature model case.
        call interp2d_expabs(ri1,rid1,slq2,sm,sh,ss,mt,mt0,dri,rri)
c
c --- Check that "slq2" has been set to 0 where it might have been negative.
      if (slq2.lt.0.) then
        write(lp,*) "************************************************"
        write(lp,*) "Error detected in turbulence module." 
        write(lp,*) "'slq2' negative in turb_2 subroutine"
     &               //" after interpolation."
        write(lp,*) "k=",k,"     slq2=",slq2
        write(lp,*) "sm=",sm,"   sh=",sh,"   ss=",ss
        write(lp,*) "ri1=",ri1,"    rid1=",rid1
        write(lp,*) "dri=",dri
        write(lp,*) "Program will stop."
        call flush(lp)
css     call xchalt('(mxgissaij)')
               stop '(mxgissaij)'
      endif
c
c
c --- Assume region contiguous with surface where foreground model is
c --- realizable has ended when get 0 "slq2".
      if (slq2.eq.0.) then
        ifbelow = 1
      endif
c
      akz=0.4*z1d(k)
      al=akz*al0/(al0+akz)
      al2=al*al
c
c --- Do not use Deardorff limitation when use (\epsilon/N^2) dimensionalization.
      if (.NOT.((ifepson2.EQ.2).AND.(ifbelow.EQ.1))) then
c --- length scale reduction by buoyancy
          if(ri1.gt.0.) then
            anlq2=slq2*ri1
            if(anlq2.gt.0.281) then  !0.281=0.53**2
              al2=0.281/anlq2*al2
              slq2=0.281/(ri1+1.E-20)
            endif
          endif
      endif !length scale reduction by buoyancy
c
      if     (.not.latdiw) then
c ---   use constant epson2 from inigiss.
        epson2 =  epson2_ref
      else
c
c ---   latitude dependent internal wave diffusion/viscosity
c ---   Gregg et. al. (2003): Reduced mixing from the breaking of
c ---   internal waves in equatorial waters, Nature 422 pp 513-515.
c
        if     (an2(k).le.0.0) then
          eplatidep  = eplatidepmin
        else
          an = SQRT(an2(k)) !N from N^2
          if    (an/acorio.lt.1.0) then   !arccosh(N/|f|) undefined
            eplatidep = eplatidepmin
          else
            eplatidep = (acorio*acosh1(an/acorio))/wave_30
            eplatidep = max(eplatidep,eplatidepmin)
          endif
        endif
        epson2 = epson2_ref*eplatidep
c
        if     (botdiw .and. an2(k).gt.epsil) then
c ---     enhanced bottom mixing
          epson2_bot = eps_bot0/an2(k) * exp((z1d(k) - zbot)/scale_bot)
          epson2     = max(epson2,epson2_bot)
        endif !botdiw
      endif !latdiw
c
c------------------------------------------------------------------------
c
c --- BEGIN SECTION .or.SALINITY MODEL BACKGROUND DIFFUSIVITY CALCULATION.
      if (ifsali.gt.0) then
c
      if (ifsalback.ge.4) then
c --- Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
c --- diffusivities calculated using the turbulence model
c --- with Ri and l_0 replaced by constants 'ri_internal' and 'back_l_0' 
c --- and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
c --- to represent a modified Dubovikov internal wave generated turbulence
c --- with constant Richardson number for ifsalback=4 case.
c
c --- Use a constant background Ri estimate. 
      if (ifsalback.EQ.4) then
        back_rit1 = 0.
        back_ric1 = 0.
        back_ri1  = ri_internal
        back_rid1 = (rid1/ri1)*ri_internal
      else
c
c --- Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
c --- diffusivities calculated using the turbulence model
c --- with l_0 replaced by a constant 'back_l_0' and 
c --- Ri by a function of Ri_d
c --- and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
c --- to represent a modified Dubovikov internal wave generated turbulence
c --- with stability-ratio dependent Richardson number for ifsalback>4 case.
c
c --- When Ri_T = 0 and Ri_C \ne 0,
c --- correctly set the angle 'theta_r' in the (Ri_T,Ri_C) plane to 'pi'/2 . 
c
c --- Skip background ra_r calculation in unstable .or.NEUTRAL* case.
c --- Set background ra_r arbitrarily to zero in these cases.
          if (ria(k).le.0.) then
            back_ra_r1 = 0.
            back_rit1  = 0.
            back_ric1  = 0.
            back_ri1   = 0.
            back_rid1  = 0.
            go to 19
          endif
c
c --- Linearly interpolate back_ra_r array to this angle in (Ri_C,Ri_T) space.
c --- Ri \equiv Ri_T + Ri_C 	; Ri_d \equiv Ri_T - Ri_C .
          rit = (ria(k) + rid(k))/2.
          ric = (ria(k) - rid(k))/2.
          ra_r = sqrt((rit**2) + (ric**2))
c
c ---  use newer better treatment of zero thermal gradient case.
c ---  find \theta_r for the Ri_T = 0 case. Treat "0/0 = 1".
          if(rit.eq.0.0) then
            if(ric.eq.0.0) then
              theta_r = atan(1.0)
            else
              theta_r = pi/2.0       ! Arctangent of infinity.
            endif
          else
            theta_r = atan(ric/rit)
          endif
c
c --- Make sure the right choice of arctan(Ri_C/Ri_T) [\theta_r] is made.
c --- Arctan only covers the range (-pi/2,pi/2) which theta_r may be outside.
c --- Want to consider statically stable case only: Ri > 0.
          if (abs(theta_r).gt.(pi/2.)) then
            write(lp,*) 
     &       "************************************************"
            write(lp,*) "Error detected in turbulence module." 
            write(lp,*) "theta_r (=",abs(theta_r),") too large"
            call flush(lp)
css         call xchalt('(mxgissaij)')
                   stop '(mxgissaij)'
          endif
          if (theta_r.lt.(-pi)/4.) then
            theta_r = theta_r + pi
          endif
c
c --- MAKE 'jtheta' A NON-NEGATIVE INDEX - ZERO AT THETA = -PI/4 .
c --- The fortran function "INT" rounds to the integer *NEAREST TO ZERO*
c --- **I.E. ROUNDS **UP** .or.NEGATIVE NUMBERS**, DOWN ONLY .or.POSITIVES.
          jtheta_r0 = INT((theta_r + (pi/4.))/deltheta_r)
          jtheta_r1 = jtheta_r0+1
c
c --- INTRODUCE 'itheta' HERE .or.THE INDEX THAT IS ZERO AT THETA=0.
          itheta_r0 = jtheta_r0 - n_theta_r_oct
          itheta_r1 = itheta_r0+1
c
c --- ***WHEN THE ANGLE IS BETWEEN THE ANGLE .or.REALIZABILITY AT INFINITY***
c --- ***AND THE LAST TABLE ANGLE BE.or. THAT CRITICAL ANGLE, *** 
c --- ***SET IT TO THE LAST TABLE ANGLE BE.or. THE CRITICAL ANGLE.****
          theta_r0 = itheta_r0*deltheta_r
          theta_r1 = itheta_r1*deltheta_r
c
          if     ((theta_r0.le.theta_rcrp).AND.
     &            (theta_r .gt.theta_rcrp)     ) then
            theta_r = theta_r1
            theta_r0 = theta_r1
            itheta_r0 = itheta_r1 
            itheta_r1 = itheta_r1+1
            theta_r1 = theta_r1 + deltheta_r
          elseif ((theta_r1.ge.theta_rcrn).AND.
     &            (theta_r .lt.theta_rcrn)     ) then
            theta_r = theta_r0
            theta_r1 = theta_r0
            itheta_r1 = itheta_r0 
            itheta_r0 = itheta_r0-1
            theta_r0 = theta_r0 - deltheta_r
          endif
c
c --- Angle in degrees.
          theta_r_deg = theta_r*180./pi
c
c --- Sound the alarm if have unrealizability outside expected range in angle.
          if ((itheta_r1.gt.3*n_theta_r_oct).or.
     &        (itheta_r0.lt. -n_theta_r_oct)    ) then
               write(lp,*) 
     &         "************************************************"
            write(lp,*) "Problem in turbulence module!"
            write(lp,*) "Unrealizability outside Ri>0 region. "
            write(lp,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
            write(lp,*) "k=",k,"  ria(k)=",ria(k),"  rid(k)=",rid(k)
            write(lp,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
            write(lp,*) "theta_r_deg =",theta_r_deg
            write(lp,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
            write(lp,*) "n_theta_r_oct=",n_theta_r_oct 
            write(lp,*) " "
            write(lp,*) "Program will stop."
            call flush(lp)
css         call xchalt('(mxgissaij)')
                   stop '(mxgissaij)'
          endif
c
          deltheta_r1 = theta_r - theta_r0
          delback_ra_r = back_ra_r(itheta_r1) - back_ra_r(itheta_r0)
          dback_ra_r_o_dtheta = delback_ra_r/deltheta_r
          back_ra_r1 = back_ra_r(itheta_r0) + 
     &                   deltheta_r1*dback_ra_r_o_dtheta
c
c --- In case choose ifrafgmax=1, ra_r is at maximum the ForeGround ra_r
c --- at the "strong" double diffusive \theta_r's 
c --- where have turbulence as Ri+> infinity. 
         ifrafglt=0
         if (ifrafgmax.EQ.1) then
           if ((theta_r.le.theta_rcrp).or.(theta_r.ge.theta_rcrn)) then
             if (back_ra_r1.gt.ra_r) then
               ifrafglt=1
               back_ra_r1=ra_r
             endif
           endif
         endif
c
        if (back_ra_r1.lt.0.) then
          write(lp,*) 
     &       "************************************************"
          write(lp,*) "Problem in turbulence module!"
          write(lp,*) "Negative bg ra_r \\equiv (Ri_T^2+Ri_C^2)^(1/2)"
          write(lp,*) "back_ra_r1 =", back_ra_r1
          write(lp,*) "theta_r =", theta_r
          write(lp,*) " "
          write(lp,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
          write(lp,*) "k=",k,"  ria(k)=",ria(k),"  rid(k)=",rid(k)
          write(lp,*) "rit=",rit,"ric=",ric
          write(lp,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
          write(lp,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
          write(lp,*) "theta_r_deg =",theta_r_deg
          write(lp,*) "n_theta_r_oct=",n_theta_r_oct 
          write(lp,*) " "
          write(lp,*) "Program will stop."
          call flush(lp)
css       call xchalt('(mxgissaij)')
                 stop '(mxgissaij)'
        endif 
c
c --- Calculate the background Ri and Ri_d .
        back_rit1 = cos(theta_r)*back_ra_r1
        back_ric1 = sin(theta_r)*back_ra_r1
        back_ri1  = back_rit1 + back_ric1
        back_rid1 = back_rit1 - back_ric1
c
      endif !ifsalback.EQ.4:else
c
c --- CALCULATE THE BACKGROUND DIMENSIONLESS TURBULENCE FUNCTIONS
c --- USING TABLE OF VALUES .or.BACKGROUND "\theta_r"'S 
c --- .or."ifbg_theta_interp"=1.
c --- Can only use theta_r table when do *not* reduce ra_r_BackGround
c --- to a smaller ra_r_ForeGround.
      if ((ifbg_theta_interp.EQ.0).or.(ifrafglt.EQ.1)) then
c
c --- Use the calculated background Ri and Ri_d in the turbulence model.
c --- Interpolate 2D table for salinity-temperature model case.
c 
        call interp2d_expabs(back_ri1,back_rid1,
     &               slq2_back,sm_back,sh_back,ss_back,mt,mt0,dri,rri)
c 
      elseif(ifbg_theta_interp.EQ.1) then
c --- Interpolate 1D table of background vs. theta_r instead.
*       if     (mnproc.eq.-99) then !always .false.
*         ! bugfix, potential I/O reduces the level of optimization
*       if     ((iglobal.eq.344.and.jglobal.eq.  1) .or.
*    &          (iglobal.eq.378.and.jglobal.eq. 32)     ) then
*         write(lp,*) 'i,j,k   = ',iglobal,jglobal,k
*         write(lp,*) '  ithet = ',itheta_r0,itheta_r1
*         write(lp,*) '  theta = ',theta_r,deltheta_r
*         write(lp,*) '  sm_r  = ',sm_r1(itheta_r0),sm_r1(itheta_r1)
*         write(lp,*) '  sh_r  = ',sh_r1(itheta_r0),sh_r1(itheta_r1)
*         write(lp,*) '  ss_r  = ',ss_r1(itheta_r0),ss_r1(itheta_r1)
*         write(lp,*) 'slq2_r  = ',slq2_r1(itheta_r0),slq2_r1(itheta_r1)
*         call flush(lp)
*       endif
        deltheta_r1 = theta_r - itheta_r0*deltheta_r
        delsm_back = sm_r1(itheta_r1) - sm_r1(itheta_r0)
        dsm_back_o_dtheta = delsm_back/deltheta_r
        sm_back = sm_r1(itheta_r0) + 
     &                   deltheta_r1*dsm_back_o_dtheta
        delsh_back = sh_r1(itheta_r1) - sh_r1(itheta_r0)
        dsh_back_o_dtheta = delsh_back/deltheta_r
        sh_back = sh_r1(itheta_r0) + 
     &                   deltheta_r1*dsh_back_o_dtheta
        delss_back = ss_r1(itheta_r1) - ss_r1(itheta_r0)
        dss_back_o_dtheta = delss_back/deltheta_r
        ss_back = ss_r1(itheta_r0) + 
     &                   deltheta_r1*dss_back_o_dtheta
        delslq2_back = slq2_r1(itheta_r1) - slq2_r1(itheta_r0)
        dslq2_back_o_dtheta = delslq2_back/deltheta_r
        slq2_back = slq2_r1(itheta_r0) + 
     &                   deltheta_r1*dslq2_back_o_dtheta
      else
        write(lp,*) "Problem with choice of background interpolation."
        write(lp,*) "ifbg_theta_interp=",ifbg_theta_interp
        write(lp,*) "ifrafglt=",ifrafglt
        write(lp,*) "Program is stopping."
        call flush(lp)
css     call xchalt('(mxgissaij)')
               stop '(mxgissaij)'
      endif
c
c --- Calculate the square of the shear from the background Richardson number.
c --- s2_back   = N^2 / ri_internal = (N^2 / S_ext^2) (S_ext^2 /ri_internal) 
c ---           = (Ri_ext / ri_internal) S_ext^2
        s2_back = (ri1/back_ri1)*s2(k)
c
c --- Set square of shear to zero for unstable density stratification.
   19   continue
        if (ri1.le.0.) then
          s2_back = 0.
        endif
c
c --- Set ill-defined S_M,H,S for unstable density stratification to zero.
        if (ri1.lt.0.) then
c
          sm_back = 0
          sh_back = 0
          ss_back = 0
        endif
c
        if     (  sm_back.lt.0.0 .or.
     &            sh_back.lt.0.0 .or.
     &            ss_back.lt.0.0 .or.
     &          slq2_back.lt.0.0     ) then
           v_back=2.0e-1
           t_back=5.0e-2
           s_back=5.0e-2
*          write(lp,'(i9,a,2i5,i3,a)')
*    &       nstep,' i,j,k=',i+i0,j+j0,k,' GISS neg. sX_back'
c
c --- Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
        elseif (ifepson2.eq.0) then
c
c --- Use the constant background l_0 lengthscale in the turbulence model.
          al0_back = back_l_0
          akz=0.4*z1d(k)
          al_back=akz*al0_back/(al0_back+akz)
          al2_back=al_back*al_back
c --- length scale reduction by buoyancy
          if(back_ri1.gt.0.) then
            anlq2_back=slq2_back*back_ri1
            if(anlq2_back.gt.0.281) then  !0.281=0.53**2
              al2_back=0.281/anlq2_back*al2_back
              slq2_back=0.281/(back_ri1+1.E-20)
            endif
          endif
c
c --- Calculate the background diffusivities.
          tmp_back=0.5*b1*al2_back*sqrt(s2_back/(slq2_back+1.E-40))
          v_back(k)=tmp_back*sm_back
          t_back(k)=tmp_back*sh_back
          s_back(k)=tmp_back*ss_back
c
c --- Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
c --- From NBp.000215-5, Volume IX : 
c --- K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
c --- K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
        else !if(ifepson2.gt.0) then
c 09/22/05 inserted if condition to fix error in calculation of tmp_back
ccc       tmp_back=0.5*b1**2*back_ri1*slq2_back*epson2
          if (back_ri1.gt.0.) then
            tmp_back=0.5*b1**2*back_ri1*slq2_back*epson2
          else
            tmp_back=0.
          endif
          v_back(k)=tmp_back*sm_back
          t_back(k)=tmp_back*sh_back
          s_back(k)=tmp_back*ss_back
        endif
c
c --- Stop if background diffusivities are negative.
        if ((v_back(k).lt.0.).or.
     &      (t_back(k).lt.0.).or.
     &      (s_back(k).lt.0.)    ) then
               write(lp,*) 
     &         "************************************************"
            write(lp,*) "Problem in turbulence module!"
            write(lp,*) "Negative Background Diffusivity."
            write(lp,*) "v_back=",v_back(k)
            write(lp,*) "t_back=",t_back(k)
            write(lp,*) "s_back=",s_back(k)
            write(lp,*) " "
            write(lp,*) "slq2_back=",slq2_back
            write(lp,*) "sm_back=",sm_back
            write(lp,*) "sh_back=",sh_back
            write(lp,*) "ss_back=",ss_back
            write(lp,*) " "
            write(lp,*) "back_ra_r1 =", back_ra_r1
            write(lp,*) "theta_r =", theta_r,
     &                   "   theta_r_deg=",theta_r_deg
            write(lp,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
            write(lp,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
            write(lp,*) " "
            write(lp,*) "k=",k,"  ria(k)=",ria(k),"  rid(k)=",rid(k)
            write(lp,*) "rit=",rit,"ric=",ric
            write(lp,*) " "
            write(lp,*) "Program will stop."
            call flush(lp)
css         call xchalt('(mxgissaij)')
                   stop '(mxgissaij)'
          endif
c
c --- Stop if background diffusivities are zero at positive Ri.
          if ((ria(k).gt.0.).and.
     &        ((v_back(k).eq.0.).or.
     &         (t_back(k).EQ.0.).or.
     &         (s_back(k).EQ.0.)    )) then
c
               write(lp,*) 
     &         "************************************************"
            write(lp,*) "Problem in turbulence module!"
            write(lp,*) "Zero Background Diffusivity in stable case."
            write(lp,*) "v_back=",v_back(k),
     &                   " t_back=",t_back(k),
     &                   " s_back=",s_back(k)
c Natassa
              write(lp,*) "tmp_back=",tmp_back
              write(lp,*) "b1=",b1
              write(lp,*) "back_ri1=",back_ri1
              write(lp,*) "epson2=",epson2
              write(lp,*) "ri1=", ri1
c
c
            write(lp,*) " "
            write(lp,*) "slq2_back=",slq2_back
            write(lp,*) "sm_back=",sm_back,
     &                   " sh_back=",sh_back,
     &                   " ss_back=",ss_back
            write(lp,*) " "
            write(lp,*) "slq2_r1(itheta_r0)=",slq2_r1(itheta_r0),
     &                   " slq2_r1(itheta_r1)=",slq2_r1(itheta_r1)
            write(lp,*) "sm_r1(itheta_r0)=",sm_r1(itheta_r0),
     &                   " sm_r1(itheta_r1)=",sm_r1(itheta_r1)
            write(lp,*) "sh_r1(itheta_r0)=",sh_r1(itheta_r0),
     &                   " sh_r1(itheta_r1)=",sh_r1(itheta_r1)
            write(lp,*) "ss_r1(itheta_r0)=",ss_r1(itheta_r0),
     &                   " ss_r1(itheta_r1)=",ss_r1(itheta_r1)
            write(lp,*) " "
            write(lp,*) "back_ra_r1 =", back_ra_r1
            write(lp,*) "theta_r =", theta_r
            write(lp,*) "theta_r_deg =", theta_r_deg
            write(lp,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
            write(lp,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
            write(lp,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
            write(lp,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
            write(lp,*) "n_theta_r_oct=",n_theta_r_oct 
            write(lp,*) "deltheta_r=",deltheta_r
            write(lp,*) " "
            write(lp,*) "k=",k,"  ria(k)=",ria(k),"  rid(k)=",rid(k)
            write(lp,*) "rit=",rit,"ric=",ric
            write(lp,*) " "
            write(lp,*) "Program will stop."
            call flush(lp)
css         call xchalt('(mxgissaij)')
                   stop '(mxgissaij)'
          endif
        endif
 20     continue
c
      endif !ifsali.gt.0
c
c --- end SECTION .or.SALINITY MODEL BACKGROUND DIFFUSIVITY CALCULATION.
c
c
c --- Write internal turbulence quantities
c --- Add S_M,H,S to outputs
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write(lp,9000) z1d(k),al,slq2,ri1,rid1,sm,sh,ss,
cdiag&                 v_back(k),t_back(k),s_back(k)
cdiag endif
c
c --- introduce foreground minimum shear squared due to internal waves
c --- to allow mixing in the unstable zero shear case
c --- Reversion to ifsalback=3 model for this purpose,
c --- based on Gargett et. al. JPO Vol.11 p.1258-71 "deep record".
         s2(k) = max(s2(k),back_s2)
c
c --- In the case where the model is realizable at 
c --- the Ri obtained from the external Shear,
c --- *but* there is a level above where it is NOT thus realizable, 
c --- USE THE "epsilon/(N^2)" DIMENSIONALIZATION .or."ifepson=2".
c --- EXCEPT if  "Ri<0" do *NOT* USE "epsilon/(N^2)" DIMENSIONALIZATION
c --- BECAUSE IT PRODUCES NEGATIVE DIFFUSIVITIES IN THIS CASE.
c --- *INSTEAD USE "l_deep^2 S", WHERE "l_deep" IS DERIVED FROM "rho" PROFILE.
c --- "|{{d \rho / dz} \over {d2 \rho / dz^2}}| takes place of MLD in l_deep".
c --- BUT *REVERT* TO "MLD" IN CASES OF FIRST TWO LEVELS.
cdiag   aldeep(k)=0.0
        if ((ifepson2.EQ.2).AND.(ifbelow.EQ.1)) then
          if (ri1.ge.0) then
              tmp=0.5*b1**2*ri1*slq2*epson2
          elseif(k.gt.2) then
                delz = z1d(k+1) - z1d(k-1)
                delth = th1d(k+1) - th1d(k-1)
                del2th = th1d(k+1) - 2.*th1d(k) + th1d(k-1)
              dzth = delth/delz
              d2zth = 4.*del2th/(delz**2)
c
c --- rdzlndzth = *Reciprocal* of Dz_{ln(Dz_{th})} = Dz_{th}/Dz2_{th} .
              if (d2zth.eq.0.0) d2zth=epsil
              rdzlndzth = dzth/d2zth
c
c --- introduce deep foreground minimum length scale due to internal waves
c --- to prevent zero lengthscale in deep zero density gradient case
c --- reversion to ifsalback=3 model for this purpose
              al0deep=max(0.17*abs(rdzlndzth),back_l_0)
              akz=0.4*z1d(k)
              aldeep(k)=akz*al0deep/(al0deep+akz)
              al2=aldeep(k)*aldeep(k)
              tmp=0.5*b1*al2*sqrt(s2(k)/(slq2+1.E-40))
          else
              tmp=0.5*b1*al2*sqrt(s2(k)/(slq2+1.E-40))
          endif
        else
            tmp=0.5*b1*al2*sqrt(s2(k)/(slq2+1.E-40))
        endif
          akm(k)=tmp*sm+v_back(k)
          akh(k)=tmp*sh+t_back(k)
          aks(k)=tmp*ss+s_back(k)
c
cdiag     tmpk(k)=tmp
c
 22   continue
c
c --- stop if DIFFUSIVITY IS NEGATIVE.
      do k =1,klist(i,j)-1
      if ((akm(k).lt.0.).or.(akh(k).lt.0.).or.(aks(k).lt.0.)) then
          write(lp,*) "Diffusivity is negative."
        write(lp,*) "k=",k
        write(lp,*) "z[cm]      tem[C]     sal[ppt]   rho[g/cm3] "//
     &                "Ri         Ri_d	   S^2[/s2]   "//
     &                "K_M[cm2/s] K_H[cm2/s] K_S[cm2/s] "
          write(*,9000) z1d(k),th1d(k),ria(k),rid(k),s2(k),
     &                  akm(k),akh(k),aks(k)
        write(lp,*) " "
        write(lp,*) "Program will stop."
        call flush(lp)
css     call xchalt('(mxgissaij)')
               stop '(mxgissaij)'
      endif
      enddo
c
c --- store new k values in the 3-d arrays
      do k=1,kk
        k1=k+1
        if(k.lt.klist(i,j)) then
          vcty(i,j,k1)=min(akm(k)*1.0e-4,difmax)
          dift(i,j,k1)=min(akh(k)*1.0e-4,difmax)
          difs(i,j,k1)=min(aks(k)*1.0e-4,difmax)
        else
          vcty(i,j,k1)=vcty(i,j,klist(i,j))
          dift(i,j,k1)=dift(i,j,klist(i,j))
          difs(i,j,k1)=difs(i,j,klist(i,j))
        endif
      enddo
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write(6,'(a,a9,a3,5a13)') 
cdiag&    'giss1dout','    nstep','  k',
cdiag&    '          tmp','       aldeep',
cdiag&    '          akm','          akh','          aks'
cdiag   do k=1,klist(i,j)
cdiag     write(6,'(a,i9,i3,1p,6e13.5)') 'giss1dout',
cdiag&          nstep,k,tmpk(k),aldeep(k),akm(k),akh(k),aks(k)
cdiag   enddo
cdiag endif
c
 9000 format(12(1pe11.3))
c
 101  format(i9,3i4,'absorbup,dn,dtemp,dsaln ',2f6.3,2f10.6)
c
      return   
      end
c***************************************************************************
c
      subroutine mxkprfbj(m,n,mm,nn,k1m,k1n, j)
ccc   use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
c --- final mixing at p points
c
      do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          call mxkprfbij(m,n,mm,nn,k1m,k1n, i,j)
        enddo
      enddo
c
      return
      end
c
c***************************************************************************
c
      subroutine mxkprfcj(m,n,mm,nn,k1m,k1n, j)
ccc   use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
c --- final velocity mixing at u,v points
c
      do l=1,isu(j)
        do i=ifu(j,l),ilu(j,l)
          call mxkprfciju(m,n,mm,nn,k1m,k1n, i,j)
        enddo
      enddo
c
      do l=1,isv(j)
        do i=ifv(j,l),ilv(j,l)
          call mxkprfcijv(m,n,mm,nn,k1m,k1n, i,j)
        enddo
      enddo
c
      return
      end
c
c***************************************************************************
c
      subroutine mxkppaij(m,n,mm,nn,k1m,k1n, i,j)
ccc   use mod_xc  ! HYCOM communication interface
c
c --- hycom version 2.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
      include 'kprf_arrays.h'
c
c -----------------------------------------------------------------
c --- kpp vertical diffusion, single j-row (part A)
c --- vertical coordinate is z negative below the ocean surface
c
c --- Large, W.C., J.C. McWilliams, and S.C. Doney, 1994: Oceanic
c --- vertical mixing: a review and a model with a nonlocal
c --- boundary layer paramterization. Rev. Geophys., 32, 363-403.
c
c --- quadratic interpolation and variable Cv from a presentation
c --- at the March 2003 CCSM Ocean Model Working Group Meeting
c --- on KPP Vertical Mixing by Gokhan Danabasoglu and Bill Large
c --- http://www.ccsm.ucar.edu/working_groups/Ocean/agendas/030320.html
c --- quadratic interpolation implemented here by 3-pt collocation, 
c --- which is slightly different to the Danabasoglu/Large approach.
c -----------------------------------------------------------------
c
      real, parameter :: difmax = 9999.0e-4  !maximum diffusion/viscosity
      real, parameter :: cv_max =    2.1     !maximum cv
      real, parameter :: cv_min =    1.7     !minimum cv
      real, parameter :: cv_bfq =  200.0     !cv scale factor
c
c local variables for kpp mixing
      real delta               ! fraction hbl lies beteen zgrid neighbors
      real zrefmn              ! nearsurface reference z, minimum
      real zref                ! nearsurface reference z
      real wref,qwref          ! nearsurface reference width,inverse
      real uref                ! nearsurface reference u
      real vref                ! nearsurface reference v
      real bref                ! nearsurface reference buoyancy
      real swfrac(kdm+1)       ! fractional surface shortwave radiation flux
      real shsq(kdm+1)         ! velocity shear squared
      real alfadt(kdm+1)       ! t contribution to density jump
      real betads(kdm+1)       ! s contribution to density jump
      real swfrml              ! fractional surface sw rad flux at ml base
      real ritop(kdm)          ! numerator of bulk richardson number
      real dbloc(kdm+1)        ! buoyancy jump across interface 
      real dvsq(kdm)           ! squared current shear for bulk richardson no.
      real zgridb(kdm+1)       ! zgrid for bottom boundary layer
      real hwide(kdm)          ! layer thicknesses in m (minimum 1mm)
      real dpmm(kdm)           !     max(onemm,dp(i,j,:,n))
      real qdpmm(kdm)          ! 1.0/max(onemm,dp(i,j,:,n))
      real pij(kdm+1)          ! local copy of p(i,j,:)
      real case                ! 1 in case A; =0 in case B
      real hbl                 ! boundary layer depth
      real hbbl                ! bottom boundary layer depth
      real rib(3)              ! bulk richardson number
      real rrho                ! double diffusion parameter
      real diffdd              ! double diffusion diffusivity scale
      real prandtl             ! prandtl number
      real fri                 ! function of Rig for KPP shear instability
      real stable              ! = 1 in stable forcing; =0 in unstable
      real dkm1(3)             ! boundary layer diffusions at nbl-1 level
      real gat1(3)             ! shape functions at dnorm=1
      real dat1(3)             ! derivative of shape functions at dnorm=1
      real blmc(kdm+1,3)       ! boundary layer mixing coefficients
      real wm                  ! momentum velocity scale
      real ws                  ! scalar velocity scale
      real dnorm               ! normalized depth
      real tmn                 ! time averaged SST
      real smn                 ! time averaged SSS
      real dsgdt               ! dsigdt(tmn,smn)
      real buoyfs              ! salinity  surface buoyancy (into atmos.)
      real buoyfl              ! total     surface buoyancy (into atmos.)
      real buoysw              ! shortwave surface buoyancy (into atmos.)
      real bfsfc               ! surface buoyancy forcing   (into atmos.)
      real bfbot               ! bottom buoyancy forcing
      real hekmanb             ! bottom ekman layer thickness
      real cormn4              ! = 4 x min. coriolis magnitude (at 4N, 4S)
      real dflsiw              ! lat.dep. internal wave diffusivity
      real dflmiw              ! lat.dep. internal wave viscosity
      real bfq                 ! buoyancy frequency
      real cvk                 ! ratio of buoyancy frequencies
      real ahbl,bhbl,chbl,dhbl ! coefficients for quadratic hbl calculation
c
      logical lhbl             ! safe to use quadratic hbl calculation
c
      integer nbl              ! layer containing boundary layer base
      integer nbbl             ! layer containing bottom boundary layer base
      integer kup2,kup,kdn     ! bulk richardson number indices
c
c --- local 1-d arrays for matrix solution
      real u1do(kdm+1),u1dn(kdm+1),v1do(kdm+1),v1dn(kdm+1),t1do(kdm+1),
     &     t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),
     &     diffm(kdm+1),difft(kdm+1),diffs(kdm+1),
     &     ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- local 1-d arrays for iteration loops
      real uold(kdm+1),vold (kdm+1),told (kdm+1),
     &     sold(kdm+1),thold(kdm+1)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
c
      real dtemp,dsaln,wq,wt,ratio,q,ghatflux,
     &     dvdzup,dvdzdn,viscp,difsp,diftp,f1,sigg,aa1,aa2,aa3,gm,gs,gt,
     &     dkmp2,dstar,hblmin,hblmax,sflux1,vtsq,
     &     vctyh,difsh,difth,zrefo,qspcifh,hbblmin,hbblmax,
     &     beta_b,beta_r,frac_b,frac_r,
     &     x0,x1,x2,y0,y1,y2
      real     sigocn,dsigdt,dsigds,dsiglocdt,dsiglocds
      external sigocn,dsigdt,dsigds,dsiglocdt,dsiglocds
c
      integer ka,kb,nlayer,ksave,iter,jrlv
c
      include 'state_eqn.h'
c
      jb=mod(j,jj)+1
      cormn4 = 4.0e-5  !4 x min. coriolis magnitude (at 4N, 4S)
c
      if     (latdiw) then
c
c ---   latitude dependent internal wave diffusion/viscosity
c ---   Gregg et. al. (2003): Reduced mixing from the breaking of
c ---   internal waves in equatorial waters, Nature 422 pp 513-515.
c ---   q is a quadratic fit to eqn 2 of Gregg, assuming N=N0 (lat=1-45).
c
        q      = max(1.0,min(45.0,abs(latij(i,j,3))))
        q      = 0.0424 + 0.04*q - 26.e-05*q**2  !q=1 at 30degN
        dflsiw = q*difsiw  !difsiw is ref.value at 30degN
        dflmiw = q*difmiw  !difmiw is ref.value at 30degN
      else
c ---   constant internal wave diffusion/viscosity
        dflsiw =   difsiw
        dflmiw =   difmiw
      endif
c
c --- locate lowest substantial mass-containing layer.
      pij(1)=p(i,j,1)
      do k=1,kk
        kn=k+nn
        dpmm( k)  =max(onemm,dp(i,j,kn))
        qdpmm(k)  =1.0/dpmm(k)
        pij(  k+1)=pij(k)+dp(i,j,kn)
        p(i,j,k+1)=pij(k+1)
      enddo
      do k=kk,1,-1
        if (dpmm(k).gt.tencm) then
          exit
        endif
      enddo
      klist(i,j)=max(k,2)  !always consider at least 2 layers
c
      if     (bblkpp) then
c ---   save klist in another array so that momtum.f can properly
c ---   distribute bottom stress in the vertical
        klstsv(i,j)=klist(i,j)
      endif
c
c --- forcing of t,s by surface fluxes. flux positive into ocean.
c --- shortwave flux penetration depends on kpar or jerlov water type.
c
      if     (jerlv0.eq.0) then
        beta_r = 2.0/onem
        beta_b = ( akpar(i,j,ls0)*ws0+akpar(i,j,ls1)*ws1
     &            +akpar(i,j,ls2)*ws2+akpar(i,j,ls3)*ws3)/onem
        frac_b = max( 0.27, 0.695 - 5.7*onem*beta_b )
        frac_r = 1.0 - frac_b
      else
        jrlv   = jerlov(i,j)
        beta_r = betard(jrlv)
        beta_b = betabl(jrlv)
        frac_r = redfac(jrlv)
        frac_b = 1.0 - frac_r
      endif
      qspcifh=1.0/spcifh
c
cdiag if (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,'(a,4f10.4)')
cdiag&   'frac[rb],beta[rb] =',
cdiag&   frac_r,frac_b,onem*beta_r,onem*beta_b
cdiag   call flush(lp)
cdiag endif
c
      do k=1,kk
        kn=k+nn
        if (thermo) then
          if     (-pij(k+1)*beta_r.gt.-10.0) then
            swfrac(k+1)=frac_r*exp(-pij(k+1)*beta_r)+
     &                  frac_b*exp(-pij(k+1)*beta_b)
          elseif (-pij(k+1)*beta_b.gt.-10.0) then
            swfrac(k+1)=frac_b*exp(-pij(k+1)*beta_b)
          else
            swfrac(k+1)=0.0
          endif
          if (k.eq.1) then
            sflux1=surflx(i,j)-sswflx(i,j)
            dtemp=(sflux1+(1.-swfrac(k+1))*sswflx(i,j))*
     &            delt1*g*qspcifh*qdpmm(k)
            dsaln=salflx(i,j)*
     &            delt1*g*        qdpmm(k) 
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,101) nstep,i+i0,j+j0,k, 
cdiag&          0.,1.-swfrac(k+1),dtemp,dsaln
cdiag         call flush(lp)
cdiag       endif
          elseif (k.le.klist(i,j)) then
            if(k.ne.klist(i,j)) then
              dtemp=(swfrac(k)-swfrac(k+1))*sswflx(i,j)*
     &              delt1*g*qspcifh*qdpmm(k)
            else
              dtemp= swfrac(k)             *sswflx(i,j)*
     &               delt1*g*qspcifh*qdpmm(k)
            endif
            dsaln=0.0
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,101) nstep,i+i0,j+j0,k,
cdiag&          1.-swfrac(k),1.-swfrac(k+1),dtemp
cdiag         call flush(lp)
cdiag       endif
          else !k.gt.klist(i,j)
            dtemp=0.0
            dsaln=0.0
          endif 
        else !.not.thermo
          dtemp=0.0
          dsaln=0.0
        endif
c
c --- modify t and s; set old value arrays at p points for initial iteration
        if (k.le.klist(i,j)) then
          temp(i,j,kn)=temp(i,j,kn)+dtemp
          saln(i,j,kn)=saln(i,j,kn)+dsaln
          th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))-thbase
          told (k)=temp(i,j,kn)
          sold (k)=saln(i,j,kn)
          if (locsig) then
            if (k.eq.1) then
              thold(k)=th3d(i,j,kn)
            else
              alfadt(k)=0.5*
     &                 (dsiglocdt(told(k-1),sold(k-1),p(i,j,k))+
     &                  dsiglocdt(told(k  ),sold(k  ),p(i,j,k)))*
     &                 (told(k-1)-told(k))
              betads(k)=0.5*
     &                 (dsiglocds(told(k-1),sold(k-1),p(i,j,k))+
     &                  dsiglocds(told(k  ),sold(k  ),p(i,j,k)))*
     &                 (sold(k-1)-sold(k))
              thold(k)=thold(k-1)-alfadt(k)-betads(k)
            endif
          else
            thold(k)=th3d(i,j,kn)
          endif
          uold(k)=.5*(u(i,j,kn)+u(i+1,j  ,kn))
          vold(k)=.5*(v(i,j,kn)+v(i  ,jb ,kn))
        endif
      enddo
c
      k=klist(i,j)
      kn=k+nn
      kb=min(k+1,kk)
      told (k+1)=temp(i,j,kb+nn)
      sold (k+1)=saln(i,j,kb+nn)
      if (locsig) then
        alfadt(k+1)=0.5*
     &            (dsiglocdt(told(k  ),sold(k  ),p(i,j,k+1))+
     &             dsiglocdt(told(k+1),sold(k+1),p(i,j,k+1)))*
     &            (told(k)-told(k+1))
        betads(k+1)=0.5*
     &            (dsiglocds(told(k  ),sold(k  ),p(i,j,k+1))+
     &             dsiglocds(told(k+1),sold(k+1),p(i,j,k+1)))*
     &            (sold(k)-sold(k+1))
        thold(k+1)=thold(k)-alfadt(k+1)-betads(k+1)
      else
        thold(k+1)=th3d(i,j,kb+nn)
      endif
      uold(k+1)=.5*(u(i,j,kn)+u(i+1,j  ,kn))
      vold(k+1)=.5*(v(i,j,kn)+v(i  ,jb ,kn))
c
c --- calculate z at vertical grid levels - this array is the z values in m
c --- at the mid-depth of each micom layer except for index klist+1, where it
c --- is the z value of the bottom
c
c --- calculate layer thicknesses in m
      do k=1,kk
        if (k.eq.1) then
          hwide(k)=dpmm(k)/onem
          zgrid(i,j,k)=-.5*hwide(k)
        else if (k.lt.klist(i,j)) then
          hwide(k)=dpmm(k)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
        else if (k.eq.klist(i,j)) then
          hwide(k)=dpmm(k)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
          zgrid(i,j,k+1)=zgrid(i,j,k)-.5*hwide(k)
        else
          hwide(k)=0.
        endif
      enddo
c
c --- perform niter iterations to execute the semi-implicit solution of the
c --- diffusion equation. at least two iterations are recommended
c
      do iter=1,niter
c
c --- calculate layer variables required to estimate bulk richardson number
c
c --- calculate nearsurface reference variables,
c --- averaged over -2*epsilon*zgrid, but no more than 8m.
        zrefmn = -4.0
        zrefo  =  1.0  ! impossible value
        do k=1,klist(i,j)
          zref=max(epsilon*zgrid(i,j,k),zrefmn)  ! nearest to zero
          if     (zref.ne.zrefo) then  ! new zref
            wref =-2.0*zref
            qwref=1.0/wref
            wq=min(hwide(1),wref)*qwref
            uref=uold(1)*wq
            vref=vold(1)*wq
            bref=-g*thref*(thold(1)+thbase)*wq
            wt=0.0
            do ka=2,k
              wt=wt+wq
              if (wt.ge.1.0) then
                exit
              endif
              wq=min(1.0-wt,hwide(ka)*qwref)
              uref=uref+uold(ka)*wq
              vref=vref+vold(ka)*wq
              bref=bref-g*thref*(thold(ka)+thbase)*wq
            enddo
          endif
          zrefo=zref
c
          ritop(k)=(zref-zgrid(i,j,k))*
     &             (bref+g*thref*(thold(k)+thbase))
          dvsq(k)=(uref-uold(k))**2+(vref-vold(k))**2
*
*         if (i.eq.itest.and.j.eq.jtest) then
*           if     (k.eq.1) then
*             write(lp,'(3a)')
*    &          ' k        z  zref',
*    &          '      u   uref      v   vref',
*    &          '      b   bref    ritop   dvsq'
*           endif
*           write(lp,'(i2,f9.2,f6.2,4f7.3,2f7.3,f9.4,f7.4)')
*    &         k,zgrid(i,j,k),zref,
*    &         uold(k),uref,vold(k),vref,
*    &         -g*thref*(thold(k)+thbase),bref,
*    &         ritop(k),dvsq(k)
*           call flush(lp)
*         endif
c
          if     (zgrid(i,j,k)*onem*beta_r.gt.-10.0) then
            swfrac(k)=frac_r*exp(zgrid(i,j,k)*onem*beta_r)+
     &                frac_b*exp(zgrid(i,j,k)*onem*beta_b)
          elseif (zgrid(i,j,k)*onem*beta_b.gt.-10.0) then
            swfrac(k)=frac_b*exp(zgrid(i,j,k)*onem*beta_b)
          else
            swfrac(k)=0.0
          endif
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,'(i9,2i5,i3,a,f8.2,f6.3)')
cdiag&          nstep,i+i0,j+j0,k,
cdiag&          '  z,swfrac =',zgrid(i,j,k),swfrac(k)
cdiag       call flush(lp)
cdiag     endif
        enddo  !k=1,klist
c
c --- calculate interface variables required to estimate interior diffusivities
        do k=1,klist(i,j)
          kb=min(k+1,kk)
          shsq  (k+1)=(uold(k)-uold(k+1))**2+(vold(k)-vold(k+1))**2
          if (.not.locsig) then
            alfadt(k+1)=.5*(dsigdt(told(k  ),sold(k  ))+
     &                      dsigdt(told(k+1),sold(k+1)))*
     &                (told(k)-told(k+1))
            betads(k+1)=.5*(dsigds(told(k  ),sold(k  ))+
     &                      dsigds(told(k+1),sold(k+1)))*
     &                (sold(k)-sold(k+1))
            dbloc(k+1)=-g*thref*(th3d(i,j,k+nn)-th3d(i,j,kb+nn))
          else
            dbloc(k+1)=-g*thref*(alfadt(k+1)+betads(k+1))
          endif
        enddo
c
c --- determine interior diffusivity profiles throughout the water column
c
        do k=1,kk+1
          vcty (i,j,k)  =0.0
          blmc(     k,1)=0.0
          blmc(     k,2)=0.0
          blmc(     k,3)=0.0
          ghats(i,j,k)  =0.0
          dift (i,j,k)  =0.0
          difs (i,j,k)  =0.0
        enddo
c
c --- shear instability plus background internal wave contributions
        do k=2,klist(i,j)+1
          if     (shinst) then
c ---       standard KPP shear instability profile
            rigr=max(0.0,dbloc(k)*(zgrid(i,j,k-1)-zgrid(i,j,k))/
     &                            (shsq(k)+epsil))
            ratio=min(rigr/rinfty,1.0)
            fri=(1.0-ratio*ratio)
            fri=fri*fri*fri
            vcty(i,j,k)=min(difm0*fri+dflmiw,difmax)
            difs(i,j,k)=min(difs0*fri+dflsiw,difmax)
          else
            vcty(i,j,k)=dflmiw
            difs(i,j,k)=dflsiw
          endif
          dift(i,j,k)=difs(i,j,k)
        enddo 
c
c --- double-diffusion (salt fingering and diffusive convection)
        if (dbdiff) then
          do k=2,klist(i,j)+1
c
c --- salt fingering case
            if (-alfadt(k).gt.betads(k) .and. betads(k).gt.0.) then
              rrho= min(-alfadt(k)/betads(k),rrho0)
              diffdd=1.-((rrho-1.)/(rrho0-1.))**2
              diffdd=dsfmax*diffdd*diffdd*diffdd
              dift(i,j,k)=dift(i,j,k)+0.7*diffdd
              difs(i,j,k)=difs(i,j,k)+diffdd
c
c --- diffusive convection case
            else if ( alfadt(k).gt.0.0 .and. betads(k).lt.0.0
     &         .and. -alfadt(k).gt.betads(k)) then
              rrho=-alfadt(k)/betads(k)
              diffdd=1.5e-6*9.*.101*exp(4.6*exp(-.54*(1./rrho-1.)))
              if (rrho.gt.0.5) then
                prandtl=(1.85-.85/rrho)*rrho
              else
                prandtl=.15*rrho
              endif
              dift(i,j,k)=dift(i,j,k)+diffdd
              difs(i,j,k)=difs(i,j,k)+prandtl*diffdd
            endif
          enddo
        endif
c
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag      write (lp,102) (nstep,iter,i+i0,j+j0,k,
cdiag&     hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag&       k=1,kk+1)
cdiag      call flush(lp)
cdiag   endif
c
c --- calculate boundary layer diffusivity profiles and match these to the
c --- previously-calculated interior diffusivity profiles
c
c --- diffusivities within the surface boundary layer are parameterized
c --- as a function of boundary layer thickness times a depth-dependent
c --- turbulent velocity scale (proportional to ustar) times a third-order
c --- polynomial shape function of depth. boundary layer diffusivities depend
c --- on surface forcing (the magnitude of this forcing and whether it is
c --- stabilizing or de-stabilizing) and the magnitude and gradient of interior
c --- mixing at the boundary layer base. boundary layer diffusivity profiles
c --- are smoothly matched to interior diffusivity profiles at the boundary
c --- layer base (the profiles and their first derivatives are continuous
c --- at z=-hbl). the turbulent boundary layer depth is diagnosed first, the
c --- boundary layer diffusivity profiles are calculated, then the boundary
c --- and interior diffusivity profiles are combined.
c
c --- minimum hbl is    top mid-layer + 1 cm or bldmin,
c --- maximum hbl is bottom mid-layer - 1 cm or bldmax.
c
        hblmin=max(              hwide(1)+0.01,bldmin)
        hblmax=min(-zgrid(i,j,klist(i,j))-0.01,bldmax)
c
c --- buoyfl = buoyancy flux (m**2/sec**3) into atmos.
c --- note: surface density increases (column is destabilized) if buoyfl > 0
c --- buoysw = shortwave radiation buoyancy flux (m**2/sec**3) into atmos.
c --- salflx, sswflx and surflx are positive into the ocean
      tmn=.5*(temp(i,j,k1m)+temp(i,j,k1n))
      smn=.5*(saln(i,j,k1m)+saln(i,j,k1n))
      dsgdt=          dsigdt(tmn,smn)
      buoyfs=g*thref*(dsigds(tmn,smn)*salflx(i,j)*thref)
      buoyfl=buoyfs+
     &       g*thref*(dsgdt          *surflx(i,j)*thref/spcifh)
      buoysw=g*thref*(dsgdt          *sswflx(i,j)*thref/spcifh)
c 
c --- diagnose the new boundary layer depth as the depth where a bulk
c --- richardson number exceeds ric
c
c --- initialize hbl and nbl to bottomed out values
        kup2=1
        kup =2
        kdn =3
        rib(kup2)=0.0
        rib(kup) =0.0
        nbl=klist(i,j)
        hbl=hblmax
c
c --- diagnose hbl and nbl
        do k=2,nbl
          case=-zgrid(i,j,k)
          bfsfc=buoyfl-swfrac(k)*buoysw
          if     (bfsfc.le.0.0) then
            stable=1.0
            dnorm =1.0
          else
            stable=0.0
            dnorm =epsilon
          endif
c
c --- compute turbulent velocity scales at dnorm, for
c --- hbl = case = -zgrid(i,j,k)
          call wscale(i,j,case,dnorm,bfsfc,wm,ws,1)
c
c --- compute the turbulent shear contribution to rib
          if     (dbloc(k).gt.0.0) then
            bfq =sqrt((dbloc(k)+epsil)/(zgrid(i,j,k-1)-zgrid(i,j,k)))
          else
            bfq =0.0
          endif
          if     (cv.ne.0.0) then
            cvk=cv
          else !frequency dependent version
            cvk=max(cv_max-cv_bfq*bfq,cv_min) !between cv_min and cv_max
          endif
          vtsq=-zgrid(i,j,k)*ws*bfq*vtc*cvk
c
c --- compute bulk richardson number at new level
          rib(kdn)=ritop(k)/(dvsq(k)+vtsq+epsil)
          if (nbl.eq.klist(i,j).and.rib(kdn).ge.ricr) then
c ---       interpolate to find hbl as the depth where rib = ricr
            if     (k.eq.2 .or. hblflg.eq.0) then  !nearest interface
              hbl = -zgrid(i,j,k-1)+0.5*hwide(k-1)
            elseif (k.lt.4 .or. hblflg.eq.1) then  !linear
              hbl = -zgrid(i,j,k-1)+
     &                 (zgrid(i,j,k-1)-zgrid(i,j,k))*
     &                 (ricr-rib(kup))/(rib(kdn)-rib(kup)+epsil)
            else !quadratic
c
c ---         Determine the coefficients A,B,C of the polynomial
c ---           Y(X) = A * (X-X2)**2 + B * (X-X2) + C
c ---         which goes through the data: (X[012],Y[012])
c
              x0 = zgrid(i,j,k-2)
              x1 = zgrid(i,j,k-1)
              x2 = zgrid(i,j,k)
              y0 = rib(kup2)
              y1 = rib(kup)
              y2 = rib(kdn)
              ahbl = ( (y0-y2)*(x1-x2) -
     &                 (y1-y2)*(x0-x2)  )/
     &               ( (x0-x2)*(x1-x2)*(x0-x1) )
              bhbl = ( (y1-y2)*(x0-x2)**2 -
     &                 (y0-y2)*(x1-x2)**2  ) /
     &               ( (x0-x2)*(x1-x2)*(x0-x1) )
              if     (abs(bhbl).gt.epsil) then
                lhbl = abs(ahbl)/abs(bhbl).gt.epsil
              else
                lhbl = .true.
              endif
              if     (lhbl) then !quadratic
c ---           find root of Y(X)-RICR nearest to X2
                chbl = y2 - ricr
                dhbl = bhbl**2 - 4.0*ahbl*chbl
                if     (dhbl.lt.0.0) then !linear
                  hbl = -(x2 + (x1-x2)*(y2-ricr)/(y2-y1+epsil))
                else
                  dhbl = sqrt(dhbl)
                  if     (abs(bhbl+dhbl).ge.
     &                    abs(bhbl-dhbl)    ) then
                    hbl = -(x2 - 2.0*chbl/(bhbl+dhbl))
                  else
                    hbl = -(x2 - 2.0*chbl/(bhbl-dhbl))
                  endif !nearest root
                endif !bhbl**2-4.0*ahbl*chbl.lt.0.0:else
              else !linear
                hbl = -(x2 + (x1-x2)*(y2-ricr)/(y2-y1+epsil))
              endif !quadratic:linear
            endif !linear:quadratic
            nbl=k
            if (hbl.lt.hblmin) then
              hbl=hblmin
              nbl=2
            endif
            if (hbl.gt.hblmax) then
              hbl=hblmax
              nbl=klist(i,j)
            endif
            exit !k-loop
          endif
c
          ksave=kup2
          kup2=kup
          kup =kdn
          kdn =ksave
        enddo  !k=1,nbl
c
c --- calculate swfrml, the fraction of solar radiation left at depth hbl
        if     (-hbl*onem*beta_r.gt.-10.0) then
          swfrml=frac_r*exp(-hbl*onem*beta_r)+
     &           frac_b*exp(-hbl*onem*beta_b)
        elseif (-hbl*onem*beta_b.gt.-10.0) then
          swfrml=frac_b*exp(-hbl*onem*beta_b)
        else
          swfrml=0.0
        endif
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag     write (lp,'(i9,2i5,i3,a,f8.2,f6.3)')
cdiag&        nstep,i+i0,j+j0,nbl,
cdiag&        '  hbl,swfrml =',hbl,swfrml
cdiag     call flush(lp)
cdiag   endif
c
c --- limit check on hbl for negative (stablizing) surface buoyancy forcing
        bfsfc=buoyfl-swfrml*buoysw
        if (bfsfc.le.0.0) then
          bfsfc=bfsfc-epsil  !insures bfsfc never=0
          hmonob(i,j)=min(-cmonob*ustar(i,j)**3/(vonk*bfsfc), hblmax)
          hbl=max(hblmin,
     &            min(hbl,
     &                hekman(i,j),
     &                hmonob(i,j)))
        else
          hmonob(i,j)=hblmax
        endif
c
c --- find new nbl and re-calculate swfrml
        nbl=klist(i,j)
        do k=2,klist(i,j)
          if (-zgrid(i,j,k).gt.hbl) then
            nbl=k
            exit
          endif
        enddo
        if     (-hbl*onem*beta_r.gt.-10.0) then
          swfrml=frac_r*exp(-hbl*onem*beta_r)+
     &           frac_b*exp(-hbl*onem*beta_b)
        elseif (-hbl*onem*beta_b.gt.-10.0) then
          swfrml=frac_b*exp(-hbl*onem*beta_b)
        else
          swfrml=0.0
        endif
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag     write (lp,'(i9,2i5,i3,a,f8.2,f6.3)')
cdiag&        nstep,i+i0,j+j0,nbl,
cdiag&        '  hbl,swfrml =',hbl,swfrml
cdiag     call flush(lp)
cdiag   endif
c
c --- find forcing stability and buoyancy forcing for final hbl values
c --- determine case (for case=0., hbl lies between -zgrid(i,j,nbl)
c --- and the interface above. for case=1., hbl lies between 
c --- -zgrid(i,j,nbl-1) and the interface below)
c
c --- velocity scales at hbl
        bfsfc=buoyfl-swfrml*buoysw
        if     (bfsfc.le.0.0) then
          bfsfc=bfsfc-epsil  !insures bfsfc never=0
          stable=1.0
          dnorm =1.0
        else
          stable=0.0
          dnorm =epsilon
        endif
        case=.5+sign(.5,-zgrid(i,j,nbl)-.5*hwide(nbl)-hbl)
c
        buoflx(i,j)=bfsfc
        bhtflx(i,j)=bfsfc-buoyfs
        mixflx(i,j)=surflx(i,j)-swfrml*sswflx(i,j)
c
        call wscale(i,j,hbl,dnorm,bfsfc,wm,ws,1)
c
c --- compute the boundary layer diffusivity profiles. first, find interior
c --- viscosities and their vertical derivatives at hbl
        ka=ifix(case+epsil)*(nbl-1)+(1-ifix(case+epsil))*nbl
        q=(hbl*onem-p(i,j,ka))*qdpmm(ka)
        vctyh=vcty(i,j,ka)+q*(vcty(i,j,ka+1)-vcty(i,j,ka))
        difsh=difs(i,j,ka)+q*(difs(i,j,ka+1)-difs(i,j,ka))
        difth=dift(i,j,ka)+q*(dift(i,j,ka+1)-dift(i,j,ka))
c
        q=(hbl+zgrid(i,j,nbl-1))/(zgrid(i,j,nbl-1)-zgrid(i,j,nbl))
        dvdzup=(vcty(i,j,nbl-1)-vcty(i,j,nbl  ))/hwide(nbl-1)
        dvdzdn=(vcty(i,j,nbl  )-vcty(i,j,nbl+1))/hwide(nbl  )
        viscp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(difs(i,j,nbl-1)-difs(i,j,nbl  ))/hwide(nbl-1)
        dvdzdn=(difs(i,j,nbl  )-difs(i,j,nbl+1))/hwide(nbl  )
        difsp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(dift(i,j,nbl-1)-dift(i,j,nbl  ))/hwide(nbl-1) 
        dvdzdn=(dift(i,j,nbl  )-dift(i,j,nbl+1))/hwide(nbl  )
        diftp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
c
        f1=-stable*c11*bfsfc/(ustar(i,j)**4+epsil) 
c
        gat1(1)=vctyh/hbl/(wm+epsil)
        dat1(1)=min(0.,-viscp/(wm+epsil)+f1*vctyh)
c
        gat1(2)=difsh/hbl/(ws+epsil)
        dat1(2)=min(0.,-difsp/(ws+epsil)+f1*difsh) 
c
        gat1(3)=difth/hbl/(ws+epsil)
        dat1(3)=min(0.,-diftp/(ws+epsil)+f1*difth)
c
c --- compute turbulent velocity scales on the interfaces
        do k=2,kk+1
          if (k.le.min(nbl+1,klist(i,j))) then
            sigg=p(i,j,k)/(hbl*onem)
            dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)
c
            call wscale(i,j,hbl,dnorm,bfsfc,wm,ws,1)
c
c --- compute the dimensionless shape functions at the interfaces
            aa1=sigg-2.
            aa2=3.-2.*sigg
            aa3=sigg-1.
c
            gm=aa1+aa2*gat1(1)+aa3*dat1(1) 
            gs=aa1+aa2*gat1(2)+aa3*dat1(2)
            gt=aa1+aa2*gat1(3)+aa3*dat1(3)
c
c --- compute boundary layer diffusivities at the interfaces
            blmc(k,1)=hbl*wm*sigg*(1.+sigg*gm)
            blmc(k,2)=hbl*ws*sigg*(1.+sigg*gs)
            blmc(k,3)=hbl*ws*sigg*(1.+sigg*gt)
c
c --- compute nonlocal transport forcing term = ghats * <ws>o
            if (nonloc) ghats(i,j,k)=(1.-stable)*cg/(ws*hbl+epsil)
          endif
        enddo
c
c
c --- enhance diffusivities on the interface closest to hbl
c
c --- first compute diffusivities at nbl-1 grid level 
        sigg=-zgrid(i,j,nbl-1)/hbl
        dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)
c
        call wscale(i,j,hbl,dnorm,bfsfc,wm,ws,1)
c
        sigg=-zgrid(i,j,nbl-1)/hbl
        aa1=sigg-2.
        aa2=3.-2.*sigg
        aa3=sigg-1.
        gm=aa1+aa2*gat1(1)+aa3*dat1(1)
        gs=aa1+aa2*gat1(2)+aa3*dat1(2)
        gt=aa1+aa2*gat1(3)+aa3*dat1(3)
        dkm1(1)=hbl*wm*sigg*(1.+sigg*gm)
        dkm1(2)=hbl*ws*sigg*(1.+sigg*gs)
        dkm1(3)=hbl*ws*sigg*(1 +sigg*gt)
c
c --- now enhance diffusivity at interface nbl
c
c --- this procedure was altered for hycom to reduce diffusivity enhancement
c --- if the interface in question is located more than dp0enh below hbl.
c --- this prevents enhanced boundary layer mixing from penetrating too far
c --- below hbl when the m.l. base is located in a very thick layer
        k=nbl-1
        delta=(hbl+zgrid(i,j,k))/(zgrid(i,j,k)-zgrid(i,j,k+1))
        q=1.-case*max(0.,min(1.,(p(i,j,k+1)-hbl*onem-dp0enh)/dp0enh))
c
        dkmp2=case*vcty(i,j,k+1)+(1.-case)*blmc(k+1,1)
        dstar=(1.-delta)**2*dkm1(1)+delta**2*dkmp2      
        blmc(k+1,1)=(1.-delta)*vcty(i,j,k+1)+delta*dstar
c
        dkmp2=case*difs(i,j,k+1)+(1.-case)*blmc(k+1,2)
        dstar=(1.-delta)**2*dkm1(2)+delta**2*dkmp2    
        blmc(k+1,2)=(1.-delta)*difs(i,j,k+1)+delta*dstar
c
        dkmp2=case*dift(i,j,k+1)+(1.-case)*blmc(k+1,3)
        dstar=(1.-delta)**2*dkm1(3)+delta**2*dkmp2     
        blmc(k+1,3)=(1.-delta)*dift(i,j,k+1)+delta*dstar
c
        if (case.eq.1.) then
          blmc(k+1,1)=max(vcty(i,j,k+1),q*blmc(k+1,1))
          blmc(k+1,2)=max(difs(i,j,k+1),q*blmc(k+1,2))
          blmc(k+1,3)=max(dift(i,j,k+1),q*blmc(k+1,3))
        endif
c
        if (nonloc) ghats(i,j,k+1)=(1.-case)*ghats(i,j,k+1)
c
c --- combine interior and boundary layer coefficients and nonlocal term
c --- zero scalar bottom diffusivities to assure no bottom flux
        do k=2,nbl
          vcty(i,j,k)=min(blmc(k,1),difmax)
          difs(i,j,k)=min(blmc(k,2),difmax)
          dift(i,j,k)=min(blmc(k,3),difmax)
        enddo
        do k=nbl+1,klist(i,j)
          ghats(i,j,k)=0.0
        enddo
        do k=klist(i,j)+1,kk+1
          vcty(i,j,k)=dflmiw
          difs(i,j,k)=dflsiw
          dift(i,j,k)=dflsiw
          ghats(i,j,k)=0.0
        enddo
c
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag     write (lp,103) (nstep,iter,i+i0,j+j0,k,
cdiag&    hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag&      ghats(i,j,k),k=1,kk+1)
cdiag     call flush(lp)
cdiag   endif
c
c --- save array dpbl=onem*hbl for ice, output and diagnosis
        dpbl(i,j)=onem*hbl
c
        if (bblkpp) then
c
c ------------------------------------------------
c
c --- begin bottom boundary layer parameterization
c
c ------------------------------------------------
c
c --- this bottom boundary algorithm follows the kpp algorithm included
c --- in the rutgers roms model. it is essentially an adaptation of the
c --- algorithm used for the surface boundary layer, involving diagnosis
c --- of the bottom boundary layer thickness hbbl using a bulk
c --- richardson number
c
c --- limit execution of this algorithm to shallow water and shelf
c --- regions and to grid points where reasonable vertical resolution
c --- exists near the bottom.
c
c --- the vertical resolution test performed here is that at least two
c --- model layers must exist within 10*thkbot of the bottom. a more
c --- stringent vertical resolution test will be performed after hbbl
c --- and nbbl are estimated
c
        if (pij(kk+1)-pij(klist(i,j)-1) .gt. thkbot*tenm .or.
     &      pij(kk+1)-pij(klist(i,j)  ) .lt. onem .or.
     &      pij(kk+1).gt.300.0*onem) then
          klstsv(i,j)=0
          go to 201
        endif
c
c --- calculate zgridb
        do k=klist(i,j)+1,1,-1
          zgridb(k)=zgrid(i,j,k)-zgrid(i,j,klist(i,j)+1)
        enddo
c
c --- calculate layer variables required to estimate bulk richardson number
c
c --- calculate nearbottom reference variables,
c --- averaged over -2*epsilon*zgrid, but no more than 8m
c
c --- assume velocity decreases linearly with depth to zero across the
c --- bottom layer
c
        zrefmn=4.0
        zrefo=-1.0  ! impossible value
        do k=klist(i,j),1,-1
          zref=min(epsilon*zgridb(k),zrefmn)  ! nearest to zero
          if     (zref.ne.zrefo) then  ! new zref
            wref =2.0*zref
            qwref=1.0/wref
            wq=min(hwide(klist(i,j)),wref)*qwref
            uref=uold(klist(i,j))*wq*wq
            vref=vold(klist(i,j))*wq*wq
            bref=-g*thref*(thold(klist(i,j))+thbase)*wq
            wt=0.0
            q=wq
            do ka=k,2,-1
              wt=wt+wq
              if (wt.ge.1.0) then
                exit
              endif
              wq=min(1.0-wt,hwide(ka)*qwref)
              if (ka.eq.klist(i,j)) then
                q=q+0.5*wq
                uref=uref+uold(ka)*wq*q
                vref=vref+vold(ka)*wq*q
              else
                uref=uref+uold(ka)*wq
                vref=vref+vold(ka)*wq
              endif
              bref=bref-g*thref*(thold(ka)+thbase)*wq
            enddo !ka
          endif
          zrefo=zref
c
          ritop(k)=max(epsil,-(zgridb(k)-zref)*
     &            (bref+g*thref*(thold(k)+thbase)))
          dvsq(k)=(uref-uold(k))**2+(vref-vold(k))**2
          dbloc(k+1)=-g*thref*(thold(k)-thold(k+1))
c
cdiag     if (i.eq.itest .and. j.eq.jtest) then
cdiag       if     (k.eq.klist(i,j)) then
cdiag         write(lp,'(3a)')
cdiag&          ' k        z  zref',
cdiag&          '      u   uref      v   vref',
cdiag&          '      b   bref    ritop   dvsq'
cdiag       endif
cdiag       write(lp,'(i2,f9.2,f6.2,4f7.3,2f7.3,f9.4,f7.4)')
cdiag&         k,zgridb(k),zref,
cdiag&         uold(k),uref,vold(k),vref,
cdiag&         -g*thref*(thold(k)+thbase),bref,
cdiag&         ritop(k),dvsq(k)
cdiag       call flush(lp)
cdiag     endif
c
        enddo  !k=klist,1,-1
c
c --- calculate bottom boundary layer diffusivity profiles and match these
c --- to the existing profiles
c
c --- minimum hbbl is 1 m, maximum is distance between bottom and one meter
c --- below the base of model layer 1
c
        hbblmin=1.0
        hbblmax=(pij(kk+1)-p(i,j,2)-onem)/onem
c
c --- buoyfl = buoyancy flux (m**2/sec**3) at the bottom due to heating by
c ---          the penetrating shortwave radiation
c --- note: bottom density increases (column is destabilized) if buoyfl > 0
c --- buoysw = shortwave radiation buoyancy flux (m**2/sec**3) at the surface
c
c --- NOTE: the convention for the bottom bl in the roms model was to use the
c --- surface net (turbulent plus radiative) heat flux to represent buoyfl
c --- this convention is not used here - instead, bottom turbulent heat
c --- flux arises entirely due to heating of the bottom by the penetrating
c --- shortwave radiation. as a result, net heat flux (buoyfl) at the bottom
c --- is zero since upward turbulent heat flux due to bottom heating is
c --- opposed by the downward shortwave radiative heat flux (it is presently
c --- assumed that no heat is absorbed by the bottom). Moving upward from
c --- the bottom, penetrating shortwave radiation acts to stabilize the
c --- water column.
c
        tmn=.5*(temp(i,j,k1m)+temp(i,j,k1n))
        smn=.5*(saln(i,j,k1m)+saln(i,j,k1n))
        dsgdt=dsigdt(tmn,smn)
        buoysw=-g*thref*(dsgdt          *sswflx(i,j)*thref/spcifh)
        buoyfl=buoysw*swfrac(klist(i,j)+1)
c
c --- diagnose the new boundary layer depth as the depth where a bulk
c --- richardson number exceeds ric
c
c --- initialize hbbl and nbbl to extreme values
        kup2=1
        kup =2
        kdn =3
        rib(kup2)=0.0
        rib(kup) =0.0
        nbbl=2
        hbbl=hbblmax
c
c --- diagnose hbbl and nbbl
        do k=klist(i,j),nbbl,-1
          case=zgridb(k)
          bfbot=buoyfl-swfrml*buoysw
          stable=.5+sign(.5,-bfbot)
          dnorm=stable+(1.-stable)*epsilon
c
c --- compute turbulent velocity scales at dnorm, for
c --- hbbl = case = zgridb(k)
          call wscale(i,j,case,dnorm,bfbot,wm,ws,2)
c
c --- compute the turbulent shear contribution to rib
          if     (dbloc(k).gt.0.0) then
            bfq =sqrt((dbloc(k)+epsil)/(zgridb(k-1)-zgridb(k)))
          else
            bfq =0.0
          endif
          if     (cv.ne.0.0) then
            cvk=cv
          else
            cvk=max(cv_max-cv_bfq*bfq,cv_min) !between cv_min and cv_max
          endif
          vtsq=zgridb(k)*ws*bfq*vtc*cvk
c
c --- compute bulk richardson number at new level
c --- interpolate to find hbbl as the depth where rib = ricr
c --- in stable or neutral conditions, hbbl can be no thicker than the
c --- minimum of the bottom ekman layer and monin-obukhov layer thicknesses
c --- ustarb is estimated in momtum.f
c
          rib(kdn)=ritop(k)/(dvsq(k)+vtsq+epsil)
          if (rib(kdn).ge.ricr) then
            if (stable.eq.0.0) then
              hekmanb=1.0e10
              hmonob(i,j) =1.0e10
            else
              hekmanb=ustarb(i,j)*(cekman*4.0)/max( cormn4,
     &                abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     &                abs(corio(i,jb ))+abs(corio(i+1,jb )))
              hmonob(i,j)=-cmonob*ustarb(i,j)**3/(vonk*bfbot)
            endif
            if     (hblflg.eq.0) then                         !nearest intf.
              hbbl = zgridb(k+1)-0.5*hwide(k+1)
            elseif (k.gt.klist(i,j)-2 .or. hblflg.eq.1) then  !linear
              hbbl = zgridb(k+1)-
     &                 (zgridb(k+1)-zgridb(k))*
     &                 (ricr-rib(kup))/(rib(kdn)-rib(kup)+epsil)
            else !quadratic
c
c ---         Determine the coefficients A,B,C of the polynomial
c ---           Y(X) = A * (X-X2)**2 + B * (X-X2) + C
c ---         which goes through the data: (X[012],Y[012])
c
              x0 = zgridb(k+2)
              x1 = zgridb(k+1)
              x2 = zgridb(k)
              y0 = rib(kup2)
              y1 = rib(kup)
              y2 = rib(kdn)
              ahbl = ( (y0-y2)*(x1-x2) -
     &                 (y1-y2)*(x0-x2)  )/
     &               ( (x0-x2)*(x1-x2)*(x0-x1) )
              bhbl = ( (y1-y2)*(x0-x2)**2 -
     &                 (y0-y2)*(x1-x2)**2  ) /
     &               ( (x0-x2)*(x1-x2)*(x0-x1) )
              if     (abs(bhbl).gt.epsil) then
                lhbl = abs(ahbl)/abs(bhbl).gt.epsil
              else
                lhbl = .true.
              endif
              if     (lhbl) then !quadratic
c ---           find root of Y(X)-RICR nearest to X2
                chbl = y2 - ricr
                dhbl = bhbl**2 - 4.0*ahbl*chbl
                if     (dhbl.lt.0.0) then !linear
                  hbbl = x2 + (x1-x2)*(y2-ricr)/(y2-y1+epsil)
                else
                  dhbl = sqrt(dhbl)
                  if     (abs(bhbl+dhbl).ge.
     &                    abs(bhbl-dhbl)    ) then
                    hbbl = x2 - 2.0*chbl/(bhbl+dhbl)
                  else
                    hbbl = x2 - 2.0*chbl/(bhbl-dhbl)
                  endif !nearest root
                endif !bhbl**2-4.0*ahbl*chbl.lt.0.0:else
              else !linear
                hbbl = x2 + (x1-x2)*(y2-ricr)/(y2-y1+epsil)
              endif !quadratic:linear
            endif
            hbbl=max(hbblmin,min(hekmanb,hmonob(i,j),hbblmax,hbbl))
            exit !k-loop
          endif
c
          ksave=kup2
          kup2=kup
          kup =kdn
          kdn =ksave
        enddo  !k=klist(i,j),nbbl,-1
c
c--- find new nbbl and re-calculate swfrml
        nbbl=2
        do k=klist(i,j),2,-1
          if (zgridb(k-1).gt.hbbl) then
            nbbl=k
            q=(zgridb(k-1)-hbbl)/(zgridb(k-1)-zgridb(k))
            swfrml=swfrac(k-1)+q*(swfrac(k)-swfrac(k-1))
            exit
          endif
        enddo  !k=klist(i,j),2,-1
c
c --- do not execute this bottom boundary layer algorithm if vertical
c --- resolution is insufficient or if very thin model layers exist
c --- in the boundary layer (the latter avoids numerical problems
c
        if (nbbl.eq.klist(i,j)+1 .or.
     &      pij(kk+1)-pij(klist(i,j)-1)+onecm.gt.hbbl*onem) then
          klstsv(i,j)=0
          go to 201
        endif
        do ka=nbbl,klist(i,j)
          if (hwide(ka).lt.0.5) then
            klstsv(i,j)=0
            go to 201
          endif
        enddo
c
c --- calculate swfrml, the fraction of solar radiation absorbed by depth hbbl
        dpbbl(i,j)=hbbl*onem
        k=nbbl
        q=(zgridb(k-1)-hbbl)/(zgridb(k-1)-zgridb(k))
        swfrml=swfrac(k-1)+q*(swfrac(k)-swfrac(k-1))
c
c --- find forcing stability and buoyancy forcing for final hbbl values
c --- determine case (for case=0., hbbl lies between -zgridb(nbbl)
c --- and the interface below. for case=1., hbbl lies between 
c --- -zgrid(nbbl+1) and the interface above)
c
c --- velocity scales at hbbl
        bfbot=buoyfl-swfrml*buoysw
        stable=.5+sign(.5,-bfbot)
        bfbot=bfbot-stable*epsil                  !insures bfbot never=0
        dnorm=stable+(1.-stable)*epsilon
        case=.5+sign(.5,zgridb(nbbl)+.5*hwide(nbbl)-hbbl)
c
        call wscale(i,j,hbbl,dnorm,bfbot,wm,ws,2)
c
c --- compute the boundary layer diffusivity profiles. first, find interior
c --- viscosities and their vertical derivatives at hbbl
        ka=ifix(case+epsil)*nbbl+(1-ifix(case+epsil))*(nbbl-1)
        q=(pij(ka)-hbbl*onem)*qdpmm(ka)
        vctyh=vcty(i,j,ka)+q*(vcty(i,j,ka+1)-vcty(i,j,ka))
        difsh=difs(i,j,ka)+q*(difs(i,j,ka+1)-difs(i,j,ka))
        difth=dift(i,j,ka)+q*(dift(i,j,ka+1)-dift(i,j,ka))
c
        q=(zgridb(nbbl-1)-hbbl)/(zgridb(nbbl-1)-zgridb(nbbl))
        dvdzup=(vcty(i,j,nbbl-1)-vcty(i,j,nbbl  ))/hwide(nbbl-1)
        dvdzdn=(vcty(i,j,nbbl  )-vcty(i,j,nbbl+1))/hwide(nbbl  )
        viscp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(difs(i,j,nbbl-1)-difs(i,j,nbbl  ))/hwide(nbbl-1)
        dvdzdn=(difs(i,j,nbbl  )-difs(i,j,nbbl+1))/hwide(nbbl  )
        difsp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(dift(i,j,nbbl-1)-dift(i,j,nbbl  ))/hwide(nbbl-1) 
        dvdzdn=(dift(i,j,nbbl  )-dift(i,j,nbbl+1))/hwide(nbbl  )
        diftp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
c
        f1=-stable*c11*bfbot/(ustarb(i,j)**4+epsil) 
c
        gat1(1)=vctyh/hbbl/(wm+epsil)
        dat1(1)=min(0.,-viscp/(wm+epsil)+f1*vctyh)
c
        gat1(2)=difsh/hbbl/(ws+epsil)
        dat1(2)=min(0.,-difsp/(ws+epsil)+f1*difsh) 
c
        gat1(3)=difth/hbbl/(ws+epsil)
        dat1(3)=min(0.,-diftp/(ws+epsil)+f1*difth)
c
c --- compute turbulent velocity scales on the interfaces
        do k=klist(i,j),2,-1
          if (k.ge.nbbl) then
            sigg=(pij(klist(i,j)+1)-pij(k-1))/(hbbl*onem)
            dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)
c
            call wscale(i,j,hbbl,dnorm,bfbot,wm,ws,2)
c
c --- compute the dimensionless shape functions at the interfaces
            aa1=sigg-2.
            aa2=3.-2.*sigg
            aa3=sigg-1.
c
            gm=aa1+aa2*gat1(1)+aa3*dat1(1) 
            gs=aa1+aa2*gat1(2)+aa3*dat1(2)
            gt=aa1+aa2*gat1(3)+aa3*dat1(3)
c
c --- compute boundary layer diffusivities at the interfaces
            blmc(k,1)=hbbl*wm*sigg*(1.+sigg*gm)
            blmc(k,2)=hbbl*ws*sigg*(1.+sigg*gs)
            blmc(k,3)=hbbl*ws*sigg*(1.+sigg*gt)
c
          else
            blmc(k,1)=0.0
            blmc(k,2)=0.0
            blmc(k,3)=0.0
          endif
        enddo
c
c --- combine interior and boundary layer coefficients and nonlocal term
c --- select maximum viscosity/diffusivity over depth range where the
c --- surface and bottom boundary layers overlap
        do k=klist(i,j),nbbl,-1
          vcty(i,j,k)=max(vcty(i,j,k),min(0.05,difmax,blmc(k,1)))
          difs(i,j,k)=max(difs(i,j,k),min(0.05,difmax,blmc(k,2)))
          dift(i,j,k)=max(dift(i,j,k),min(0.05,difmax,blmc(k,3)))
        enddo
c
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag     write (lp,103) (nstep,iter,i+i0,j+j0,k,
cdiag&    hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag&      ghats(i,j,k),k=kk,1,-1)
cdiag     call flush(lp)
cdiag      print *,'hbl,hbbl = ',hbl,hbbl
cdiag   endif
c
        endif    ! bblkpp
c
 201    continue ! exited from bblkpp
c
        if (iter.lt.niter) then
c
c --- perform the vertical mixing at p points
c
          do k=1,klist(i,j)
            kn=k+nn
            difft(k+1)=dift(i,j,k+1)
            diffs(k+1)=difs(i,j,k+1)
            diffm(k+1)=vcty(i,j,k+1)
            ghat(k+1)=ghats(i,j,k+1)
            t1do(k)=temp(i,j,kn)
            s1do(k)=saln(i,j,kn)
            u1do(k)=.5*(u(i,j,kn)+u(i+1,j  ,kn))
            v1do(k)=.5*(v(i,j,kn)+v(i  ,jb ,kn))
            hm(k)=hwide(k)
            zm(k)=zgrid(i,j,k)
          enddo
c
          nlayer=klist(i,j)
          k=nlayer+1
          ka=min(k,kk)
          difft(k)=0.0
          diffs(k)=0.0
          diffm(k)=0.0
          ghat(k)=0.0
          t1do(k)=temp(i,j,ka+nn)
          s1do(k)=saln(i,j,ka+nn)
          u1do(k)=u1do(k-1)
          v1do(k)=v1do(k-1)
          zm(k)=zgrid(i,j,k)
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c         tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
c         tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
c
          do k=1,nlayer
            dzb(k)=zm(k)-zm(k+1)
          enddo
c
          tri(1,1)=delt1/(hm(1)*dzb(1))
          tri(1,0)=0.
          do k=2,nlayer
            tri(k,1)=delt1/(hm(k)*dzb(k))
            tri(k,0)=delt1/(hm(k)*dzb(k-1))
          enddo
c
c --- solve the diffusion equation
c --- salflx, sswflx and surflx are positive into the ocean
c
c --- t solution
          ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
          call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
          call tridrhs(hm,t1do,difft,ghat,ghatflux,tri,nlayer,rhs,delt1)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)
c
c --- s solution
          ghatflux=-salflx(i,j)*thref
          call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
          call tridrhs(hm,s1do,diffs,ghat,ghatflux,tri,nlayer,rhs,delt1)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)
c
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,104) (nstep,iter,i+i0,j+j0,k,
cdiag&        hm(k),t1do(k),t1dn(k),s1do(k),s1dn(k),
cdiag&        0.0,0.0,
cdiag&        k=1,nlayer)
cdiag       call flush(lp)
cdiag     endif
c
c --- u solution
          call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
          do k=1,nlayer
            rhs(k)=u1do(k)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm)
c
c --- v solution
          do k=1,nlayer
            rhs(k)=v1do(k)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm)
c
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,105) (nstep,iter,i+i0,j+j0,k,
cdiag&        hm(k),u1do(k),u1dn(k),v1do(k),v1dn(k),k=1,nlayer)
cdiag       call flush(lp)
cdiag     endif
c
c --- reset old variables in preparation for next iteration
          do k=1,nlayer+1
            told(k)=t1dn(k)
            sold(k)=s1dn(k)
            if (locsig) then
              if (k.eq.1) then
                thold(k)=sigocn(told(k),sold(k))-thbase
              else
                alfadt(k)=0.5*
     &                   (dsiglocdt(told(k-1),sold(k-1),p(i,j,k))+
     &                    dsiglocdt(told(k  ),sold(k  ),p(i,j,k)))*
     &                   (told(k-1)-told(k))
                betads(k)=0.5*
     &                   (dsiglocds(told(k-1),sold(k-1),p(i,j,k))+
     &                    dsiglocds(told(k  ),sold(k  ),p(i,j,k)))*
     &                   (sold(k-1)-sold(k))
                thold(k)=thold(k-1)-alfadt(k)-betads(k)
              endif
            else
              thold(k)=sigocn(told(k),sold(k))-thbase
            endif
            if (iter.lt.niter) then
              uold(k)=u1dn(k)
              vold(k)=v1dn(k)
            endif
          enddo
        endif                         ! iter < niter
c
      enddo                           ! iteration loop
c
 101  format(i9,2i5,i3,'absorbup,dn,dtemp,dsaln ',2f6.3,2f10.6)
 102  format(25x,'   thick      viscty    t diff    s diff  '
     &     /(i9,i2,2i5,i3,2x,4f10.2))
 103  format(25x,'   thick      viscty    t diff    s diff   nonlocal'
     &     /(i9,i2,2i5,i3,2x,4f10.2,f11.6))
 104  format(25x,
     &     '  thick   t old   t new   s old   s new trc old trc new'
     &     /(i9,i2,2i5,i3,1x,f9.2,4f8.3,2f7.4))
 105  format(25x,'   thick   u old   u new   v old   v new'
     &     /(i9,i2,2i5,i3,1x,f10.2,4f8.3))
c
      return
      end
c
c***************************************************************************
c
      subroutine mxkprfbij(m,n,mm,nn,k1m,k1n, i,j)
ccc   use mod_xc  ! HYCOM communication interface
c
c --- hycom version 2.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
      include 'kprf_arrays.h'
c
c -------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row (part B)
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------
c
c --- perform the final vertical mixing at p points
c
c --- local 1-d arrays for matrix solution
      real t1do(kdm+1),t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),
     &     tr1do(kdm+1,ntrcr),tr1dn(kdm+1,ntrcr),
     &     difft(kdm+1),diffs(kdm+1),
     &     ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm),
     &     totemo,totemn,tosalo,tosaln,tndcyt,tndcys,
     &     totrco(ntrcr),totrcn(ntrcr),trscal(ntrcr)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
c
      real    ghatflux,sigocn
      external sigocn
      integer ka,ktr,nlayer
c
      include 'state_eqn.h'
c
      nlayer=klist(i,j)
c
      do k=1,nlayer
        kn=k+nn
        difft(k+1)=dift(i,j,k+1)
        diffs(k+1)=difs(i,j,k+1)
        ghat(k+1)=ghats(i,j,k+1)
        t1do(k)=temp(i,j,kn)
        s1do(k)=saln(i,j,kn)
        t1dn(k)=t1do(k)
        s1dn(k)=s1do(k)
        if (dotrcr) then
          do ktr=1,ntrcr
            tr1do(k,ktr)=tracer(i,j,k,ktr)
            tr1dn(k,ktr)=tr1do(k,ktr)
          enddo
        end if
        hm(k)=max(onemm,dp(i,j,kn))/onem
        zm(k)=zgrid(i,j,k)
      enddo
c
      do k=nlayer+1,kk+1
        ka=min(k,kk)
        difft(k)=0.0
        diffs(k)=0.0
        ghat(k)=0.0
        t1do(k)=temp(i,j,ka+nn)
        s1do(k)=saln(i,j,ka+nn)
        t1dn(k)=t1do(k)
        s1dn(k)=s1do(k)
        if (dotrcr) then
          do ktr=1,ntrcr
            tr1do(k,ktr)=tracer(i,j,ka,ktr)
            tr1dn(k,ktr)=tr1do(k,ktr)
          enddo
        end if
        zm(k)=zm(k-1)-0.001
      end do
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c     tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
c     tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
c
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
c --- salflx, sswflx and surflx are positive into the ocean
c
c --- t solution
      ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
      call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,t1do,difft,ghat,ghatflux,tri,nlayer,rhs,delt1)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)
c
c --- s solution
      ghatflux=-salflx(i,j)*thref
      call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,s1do,diffs,ghat,ghatflux,tri,nlayer,rhs,delt1)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)
c
      if (dotrcr) then
c --- tracer solution (using time step trcfrq*baclin)
c
      tri(1,1)=trcfrq*baclin/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=trcfrq*baclin/(hm(k)*dzb(k))
        tri(k,0)=trcfrq*baclin/(hm(k)*dzb(k-1))
      enddo
c
        do ktr=1,ntrcr
cdiag     if (i.eq.itest .and. j.eq.jtest) write (*,'(a,i2/1p(8e10.1))')
cdiag.      'old tracer',ktr,(tr1do(k,ktr),k=1,nlayer)
          ghatflux=0.
          call tridrhs(hm,
     &        tr1do(1,ktr),diffs,ghat,ghatflux,tri,nlayer,rhs,
     &              trcfrq*baclin)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,
     &                 tr1do(1,ktr),tr1dn(1,ktr),diffs)
cdiag     if (i.eq.itest .and. j.eq.jtest) write (*,'(a,i2/1p(8e10.1))')
cdiag.      'new tracer',ktr,(tr1dn(k,ktr),k=1,nlayer)
        enddo
      end if
c
c --- check conservation of column integrals
      totemo=t1do(1)*dp(i,j,k1n)
      totemn=t1dn(1)*dp(i,j,k1n)
      tosalo=s1do(1)*dp(i,j,k1n)
      tosaln=s1dn(1)*dp(i,j,k1n)
      if (dotrcr) then
        do ktr=1,ntrcr
          totrco(ktr)=tr1do(1,ktr)*dp(i,j,k1n)
          totrcn(ktr)=tr1dn(1,ktr)*dp(i,j,k1n)
          trscal(ktr)=max(abs(tr1do(1,ktr)),abs(tr1dn(1,ktr)))
        end do
      end if
      do k=2,kk
        kn=k+nn
        totemo=totemo+t1do(k)*dp(i,j,kn)
        totemn=totemn+t1dn(k)*dp(i,j,kn)
        tosalo=tosalo+s1do(k)*dp(i,j,kn)
        tosaln=tosaln+s1dn(k)*dp(i,j,kn)
        if (dotrcr) then
          do ktr=1,ntrcr
            totrco(ktr)=totrco(ktr)+tr1do(k,ktr)*dp(i,j,kn)
            totrcn(ktr)=totrcn(ktr)+tr1dn(k,ktr)*dp(i,j,kn)
            trscal(ktr)=max(trscal(ktr),abs(tr1do(k,ktr)),
     &                                  abs(tr1dn(k,ktr)))
          end do
        end if                !  dotrcr
      end do
      tndcyt=totemn-totemo
      tndcys=tosaln-tosalo
      totemn=10.*pbot(i,j)
      tosaln=35.*pbot(i,j)
      if (abs(tndcyt).gt.acurcy*totemn) write (lp,101) i,j,
     .  '  mxkprf - bad temp.intgl.',totemo,tndcyt,tndcyt/totemn
      if (abs(tndcys).gt.acurcy*tosaln) write (lp,101) i,j,
     .  '  mxkprf - bad saln.intgl.',tosalo,tndcys,tndcys/tosaln
      if (dotrcr) then
        do ktr=1,ntrcr
          tndcyt=totrcn(ktr)-totrco(ktr)
          if (abs(tndcyt).lt.1.e-199) tndcyt=0.
          totemn=trscal(ktr)*pbot(i,j)
          if (abs(tndcyt).gt.acurcy*totemn) write (lp,101) i,j,
     .    '  mxkprf - bad trcr.intgl.',totrco(ktr),tndcyt,tndcyt/totemn
        end do
      end if
 101  format (2i5,a,1p,2e16.8,e9.1)
c
c --- adjust t, s, th, arrays
      do k=1,klist(i,j)
        kn=k+nn
        temp(i,j,kn)=t1dn(k)
        saln(i,j,kn)=s1dn(k)
        th3d(i,j,kn)=sigocn(t1dn(k),s1dn(k))-thbase
        if (dotrcr) then
          do ktr=1,ntrcr
            tracer(i,j,k,ktr)=tr1dn(k,ktr)
          enddo
        end if
      enddo
c
      return
      end
c
c***************************************************************************
c
      subroutine mxkprfciju(m,n,mm,nn,k1m,k1n, i,j)
ccc   use mod_xc  ! HYCOM communication interface
c
c --- hycom version 2.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
      include 'kprf_arrays.h'
c
c -------------------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row, momentum at u grid points
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------------------
c
c local variables for kpp mixing
c
c --- local 1-d arrays for matrix solution
      real u1do(kdm+1),u1dn(kdm+1),
     &     diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm),
     &     toto,totn,tndcy
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
c
      real presu
      integer nlayer
c
      nlayer=1
      presu=0.
      do k=1,kk
        u1do(k)=0.
        if (presu.lt.depthu(i,j)-tencm) then
          diffm(k+1)=.5*(vcty(i,j,k+1)+vcty(i-1,j,k+1))
          u1do(k)=u(i,j,k+nn)
          hm(k)=max(onemm,dpu(i,j,k+nn))/onem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presu=presu+dpu(i,j,k+nn)
          nlayer=k
ccc        else
ccc          exit
        endif
        u1dn(k)=u1do(k)
      enddo
c
      k=nlayer+1
      diffm(k)=0.
      u1do(k)=u1do(k-1)
      zm(k)=zm(k-1)-.5*hm(k-1)
c
c --- compute factors for coefficients of tridiagonal matrix elements.
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
      call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)= u1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm)
      do k=1,nlayer
        u(i,j,k+nn)=u1dn(k)
      enddo
c
c --- check conservation of column integrals
      toto=u1do(1)*dpu(i,j,k1n)
      totn=u1dn(1)*dpu(i,j,k1n)
      do k=2,kk
        kn=k+nn
        toto=toto+u1do(k)*dpu(i,j,kn)
        totn=totn+u1dn(k)*dpu(i,j,kn)
      end do
      tndcy=totn-toto
      totn=0.1*pbot(i,j)
      if (abs(tndcy).gt.acurcy*totn) write (lp,101) i,j,
     .  '  mxkprf - bad u intgl.',toto,tndcy,tndcy/totn
 101  format (2i5,a,1p,2e16.8,e9.1)
c
cdiag if (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,106) (nstep,i+i0,j+j0,k,
cdiag&    hm(k),u1do(k),u1dn(k),k=1,nlayer)
cdiag   call flush(lp)
cdiag endif
      return
 106  format(23x,'   thick   u old   u new'/(i9,2i5,i3,1x,f10.3,2f8.3))
      end
c
c***************************************************************************
c
      subroutine mxkprfcijv(m,n,mm,nn,k1m,k1n, i,j)
ccc   use mod_xc  ! HYCOM communication interface
c
c --- hycom version 2.1
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
      include 'kprf_arrays.h'
c
c --------------------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row , momentum at v grid points
c --- vertical coordinate is z negative below the ocean surface
c --------------------------------------------------------------------------
c
c local variables for kpp mixing
c
c --- local 1-d arrays for matrix solution
      real v1do(kdm+1),v1dn(kdm+1),
     &     diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm),
     &     toto,totn,tndcy
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
c
      real presv
      integer nlayer
c
      ja=mod(j-2+jj,jj)+1
      nlayer=1
      presv=0.
      do k=1,kk
        v1do(k)=0.
        if (presv.lt.depthv(i,j)-tencm) then
          diffm(k+1)=.5*(vcty(i,j,k+1)+vcty(i,ja ,k+1))
          v1do(k)=v(i,j,k+nn)
          hm(k)=max(onemm,dpv(i,j,k+nn))/onem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presv=presv+dpv(i,j,k+nn)
          nlayer=k
ccc        else
ccc          exit
        endif
        v1dn(k)=v1do(k)
      enddo
c
      k=nlayer+1
      diffm(k)=0.
      v1do(k)=v1do(k-1)
      zm(k)=zm(k-1)-.5*hm(k-1)
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
      call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=v1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm)
      do k=1,nlayer
        v(i,j,k+nn)=v1dn(k)
      enddo
c
c --- check conservation of column integrals
      toto=v1do(1)*dpv(i,j,k1n)
      totn=v1dn(1)*dpv(i,j,k1n)
      do k=2,kk
        kn=k+nn
        toto=toto+v1do(k)*dpv(i,j,kn)
        totn=totn+v1dn(k)*dpv(i,j,kn)
      end do
      tndcy=totn-toto
      totn=0.1*pbot(i,j)
      if (abs(tndcy).gt.acurcy*totn) write (lp,101) i,j,
     .  '  mxkprf - bad v intgl.',toto,tndcy,tndcy/totn
 101  format (2i5,a,1p,2e16.8,e9.1)
c
cdiag if (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,107) (nstep,i+i0,j+j0,k,
cdiag&    hm(k),v1do(k),v1dn(k),k=1,nlayer)
cdiag   call flush(lp)
cdiag endif
      return
 107  format(23x,'   thick   v old   v new'/(i9,2i5,i3,1x,f10.3,2f8.3))
      end
c
c***************************************************************************
c
      subroutine wscale(i,j,zlevel,dnorm,bfsfc,wm,ws,isb)
ccc   use mod_xc  ! HYCOM communication interface
c
      implicit none
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'kprf_scalars.h'
c
      integer isb
      real    zlevel,dnorm,bfsfc,wm,ws
c
c -------------------------------------------------------------------------
c --- subroutine to compute turbulent velocity scales for kpp mixing scheme
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------------------
c
c --- see inikpp for initialization of /kppltr/ and other constants.
c
      integer    nzehat,nustar
      parameter (nzehat=890,nustar=192)
c
      real, dimension (0:nzehat+1,0:nustar+1) ::
     & wmt            ! momentum velocity scale table
     &,wst            ! scalar   velocity scale table
      common/kppltr/ wmt,wst
      save  /kppltr/
c
      real    zdiff,udiff,zfrac,ufrac,
     &        wam,wbm,was,wbs,ucube,zehat,ust,zlv
      integer iz,izp1,ju,jup1
c
c --- isb=1 for surface, isb=2 for bottom
c
      if(isb.eq.1) then
        ust=ustar(i,j)
        zlv=zlevel
      else
        ust=ustarb(i,j)
        zlv=-zlevel
      endif
c
c --- use lookup table for zehat < zmax  only;  otherwise use stable formulae
c
      zehat=-vonk*dnorm*zlv*bfsfc
      if (zehat.le.zmax) then
        zdiff=zehat-zmin
        iz=int(zdiff/deltaz)
        iz=max(min(iz,nzehat),0)
        izp1=iz+1
c
        udiff=ust-umin
        ju=int(udiff/deltau)
        ju=max(min(ju,nustar),0)
        jup1=ju+1
c
        zfrac=zdiff/deltaz-iz
        ufrac=udiff/deltau-ju
c
        wam=(1.-zfrac)*wmt(iz,jup1)+zfrac*wmt(izp1,jup1)
        wbm=(1.-zfrac)*wmt(iz,ju  )+zfrac*wmt(izp1,ju  )
        wm =(1.-ufrac)*wbm         +ufrac*wam
c
        was=(1.-zfrac)*wst(iz,jup1)+zfrac*wst(izp1,jup1)
        wbs=(1.-zfrac)*wst(iz,ju  )+zfrac*wst(izp1,ju  )
        ws =(1.-ufrac)*wbs         +ufrac*was
c
      else
c
        ucube=ust**3
        wm=vonk*ust*ucube/(ucube+c11*zehat)
        ws=wm
c
      endif
c
      return
      end
c
c
c> Revision history:
c>
c> Jun  2000 - conversion to SI units.
c> Jul  2000 - included wscale in this file to facilitate in-lining
c> May  2002 - buoyfl (into the atmos.), calculated here
c> Nov  2002 - added kPAR based turbidity
c> Nov  2002 - hmonob,mixflx,buoflx,bhtflx saved for diagnostics
c> Mar  2003 - added GISS mixed layer
c> May  2003 - added bldmin and bldmax to KPP
c> Jan  2004 - added latdiw to KPP and MY
c> Jan  2004 - added bblkpp to KPP (bottom boundary layer)
c> Jan  2004 - cv can now depend on bouyancy freqency
c> Jan  2004 - added hblflg to KPP
