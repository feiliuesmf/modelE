      subroutine diffus(dtime)
!@sum  diffus updates u,v,t,q due to 
!@sum  turbulent transport throughout all GCM layers
!@sum  using a second order closure (SOC)
!@sum  turbulence model developed at GISS, 2000.
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from diffB347D6M20)
!@cont diffus,getdz,dout,diff_uv,diff_tq,diff_e
!@cont lgcm,kgcm,ave_uv_to_tcell,ave_to_ucell
!@var u 3d west-east wind component
!@var v 3d south-north wind component 
!@var t 3d potential temperature
!@var q 3d relative humidity
!@var p 2-d pressure
!@var dtime time step
!@var mout 0:don't call dout; 1:call dout
!@var itest longitude at which to call dout
!@var jtest latitude at which to call dout
      USE DYNAMICS, only : pk,plij
      USE MODEL_COM, only :
     *      im,jm,lm,sig,sige,psf,ptop,ls1,u,v,t,q,p
     *     ,vt_on      
      USE CONSTANT, only : kapa,deltx
      USE PBLCOM, only : tsavg,qsavg,dclev,uflux,vflux,tflux,qflux,egcm
      USE GEOM, only : imaxj

      IMPLICIT NONE

      real*8, intent(in) :: dtime

      real*8, dimension(lm) :: uij,vij,tij,pij,qij,eij
      real*8, dimension(lm) :: u0ij,v0ij,t0ij,q0ij,e0ij
      real*8, dimension(lm) :: km,kh,ke,lscale,gm,gh,rhoij,rhoeij
      real*8, dimension(lm-1) :: dzij,dzeij
      real*8, dimension(im,jm,lm) :: rho,rhoe
      real*8, dimension(im,jm,lm-1) :: dz,dzedge
      real*8, dimension(lm) :: peij,ttest
      real*8, dimension(im,jm,lm) :: u_tcell,v_tcell,tv_ucell,t_virtual
      real*8, dimension(im,jm,lm) :: km_gcm,km_gcm_ucell
     2        ,dz_ucell,dzedge_ucell,rho_ucell,rhoe_ucell
      real*8, dimension(im,jm) :: p_ucell,tsavg_ucell,qsavg_ucell
      real*8, dimension(im,jm) :: uflux_ucell,vflux_ucell

      integer, parameter :: mout=0,itest=52,jtest=33,itmax=2
      real*8, parameter :: tol=1.d-4,qmin=1.d-12,p00=1000.d0
      integer, save :: ifirst=1
      real*8, save :: p1000k=1.d0
      real*8 :: uflx,vflx,tflx,qflx,pijgcm,pl,rvx
      real*8 :: temp0,ps,ustar2,dbll,reserv,test,check
      integer :: loc,icount
      integer :: i,j,l,iter  !@var i,j,l,iter loop variable

      if (ifirst.eq.1) then
        p1000k=p00**kapa
        ifirst=0
      endif

      if(.not. vt_on) then
          rvx=0.d0
      else
          rvx=deltx
      endif

      !  convert input T to virtual T
      do l=1,lm
        do j=1,jm
          do i=1,im
            t_virtual(i,j,l)=t(i,j,l)*(1.d0+RVX*Q(i,j,l))
          end do
        end do
      end do

c     integrate T,Q equations at tcells

      ! get u_tcell and v_tcell at t-cells
      call ave_uv_to_tcell(u,v,u_tcell,v_tcell,im,jm,lm)

      call getdz(t_virtual,p,sig,sige,ptop,psf,
     2           dz,dzedge,rho,rhoe,tsavg,qsavg,rvx,ls1,im,jm,lm)

      loop_j_tq: do j=1,jm
        loop_i_tq: do i=1,imaxj(j)
          do l=1,lm
            uij(l)=u_tcell(i,j,l)
            vij(l)=v_tcell(i,j,l)
            tij(l)=t_virtual(i,j,l)*p1000k  !virtual,potential temp.
            pij(l)=100.d0*(plij(l,i,j)*sig(l)+ptop)
            qij(l)=q(i,j,l)
            eij(l)=egcm(l,i,j)
            rhoeij(l)=rhoe(i,j,l)
            rhoij(l)=rho(i,j,l)
            u0ij(l)=uij(l)
            v0ij(l)=vij(l)
            t0ij(l)=tij(l)
            q0ij(l)=qij(l)
            e0ij(l)=eij(l)
          end do
          do l=1,lm-1
            dzij(l)=dz(i,j,l)
            dzeij(l)=dzedge(i,j,l)
          end do

          uflx  =uflux(i,j)
          vflx  =vflux(i,j)
          ustar2=sqrt(uflx*uflx+vflx*vflx)
          tflx  =tflux(i,j)
          qflx  =qflux(i,j)

          call lgcm(lscale,uij,vij,tij,eij,dzij,dzeij,rhoij,lm)
          call kgcm(km,kh,ke,gm,gh,uij,vij,tij,eij,lscale,dzij
     2             ,dbll,lm)
          call diff_e(e0ij,eij,km,kh,ke,lscale,uij,vij,tij,
     2           dzij,dzeij,rhoij,rhoeij,dtime,ustar2,lm)
          call diff_tq(t0ij,tij,kh,dzij,dzeij,
     2                 rhoij,rhoeij,tflx,dtime,lm)
          call diff_tq(q0ij,qij,kh,dzij,dzeij,
     2                 rhoij,rhoeij,qflx,dtime,lm)

          do 300 iter=1,itmax

            do l=2,lm-1
              ttest(l)=tij(l)
            end do

            call lgcm(lscale,uij,vij,tij,eij,dzij,dzeij,rhoij,lm)
            call kgcm(km,kh,ke,gm,gh,uij,vij,tij,eij,lscale,dzij
     2               ,dbll,lm)
            call diff_e(e0ij,eij,km,kh,ke,lscale,uij,vij,tij,
     2           dzij,dzeij,rhoij,rhoeij,dtime,ustar2,lm)
            call diff_tq(t0ij,tij,kh,dzij,dzeij,
     2                   rhoij,rhoeij,tflx,dtime,lm)
            call diff_tq(q0ij,qij,kh,dzij,dzeij,
     2                   rhoij,rhoeij,qflx,dtime,lm)

            test=0.d0
            icount=0
            do l=2,lm-1
              check=abs(ttest(l)-tij(l))
              if (check.gt.0.d0) then
                test=test+check/abs(ttest(l)+tij(l))
                icount=icount+1
              endif
            end do
c           test=test/float(lm-2)
            test=test/(float(icount)+1.d-40)
            if (test.lt.tol) exit

300       continue

          do l=1,lm
            t_virtual(i,j,l)=tij(l)/p1000k
            q(i,j,l)=max(qij(l),qmin)
            t(i,j,l)=t_virtual(i,j,l)/(1.d0+RVX*Q(i,j,l))
            egcm(l,i,j)=eij(l)
            km_gcm(i,j,l)=km(l)
          end do
          dclev(i,j)=dbll
c
c         Write out diagnostics if at the correct grid point:
c
          if (mout.eq.1.and.(i.eq.itest).and.(j.eq.jtest)) then
            ps=100.d0*(p(i,j)+ptop)
            do l=1,lm
              if (l.ge.ls1) then
                peij(l)=100.d0*((psf-ptop)*sige(l)+ptop)
              else
                peij(l)=100.d0*(p(i,j)*sige(l)+ptop)
              endif
            end do
            call dout(uij,vij,tij,pij,peij,qij,eij,dzij,dzeij,
     2                rhoij,rhoeij,u0ij,v0ij,t0ij,q0ij,
     3                km,kh,ke,gm,gh,lscale,ps,reserv,tsavg(i,j),
     4                uflx,vflx,tflx,qflx,itest,jtest,lm)
          endif

        end do loop_i_tq
      end do loop_j_tq

c     integrate U,V equations at ucells

      ! average some quantities at u-cells
c     call ave_to_ucell(p,p_ucell,im,jm,1)
      call ave_to_ucell(uflux,uflux_ucell,im,jm,1)
      call ave_to_ucell(vflux,vflux_ucell,im,jm,1)
c     call ave_to_ucell(tsavg,tsavg_ucell,im,jm,1)
c     call ave_to_ucell(qsavg,qsavg_ucell,im,jm,1)

c     call ave_to_ucell(t_virtual,tv_ucell,im,jm,lm)
      call ave_to_ucell(km_gcm,km_gcm_ucell,im,jm,lm)
 
      call ave_to_ucell(dz,dz_ucell,im,jm,lm)
      call ave_to_ucell(dzedge,dzedge_ucell,im,jm,lm)
      call ave_to_ucell(rho,rho_ucell,im,jm,lm)
      call ave_to_ucell(rhoe,rhoe_ucell,im,jm,lm)

      loop_j_uv: do j=2,jm
        loop_i_uv: do i=1,im
          do l=1,lm
            uij(l)=u(i,j,l)
            vij(l)=v(i,j,l)
            rhoeij(l)=rhoe_ucell(i,j,l)
            rhoij(l)=rho_ucell(i,j,l)
            u0ij(l)=uij(l)
            v0ij(l)=vij(l)
          end do
          do l=1,lm-1
            dzij(l)=dz_ucell(i,j,l)
            dzeij(l)=dzedge_ucell(i,j,l)
            km(l)=km_gcm_ucell(i,j,l)
          end do

          uflx  =uflux_ucell(i,j)
          vflx  =vflux_ucell(i,j)

          do iter=1,itmax
            call diff_uv(u0ij,v0ij,uij,vij,km,dzij,dzeij,rhoij,rhoeij,
     2                   uflx,vflx,dtime,lm)
          end do !loop iter

          do l=1,lm
            u(i,j,l)= uij(l)
            v(i,j,l)= vij(l)
          end do
 
        end do loop_i_uv 
      end do loop_j_uv

      return
      end subroutine diffus

      subroutine getdz(tv,p,sig,sige,ptop,psf,
     2   dz,dzedge,rho,rhoe,tsavg,qsavg,rvx,ls1,im,jm,lm)
!@sum  getdz computes the 3d finite difference dz and dzedge
!@sum  as well as the 3d density rho and rhoe
!@sum  called at the primary cells (t-cells)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var  dz(i,j,l) z(i,j,l+1) - z(i,j,l)
!@var  dzedge(i,j,l) zedge(i,j,l+1) - zedge(i,j,l)
!@var  z vertical coordinate associated with SIG(l)
!@var  zedge vertical coordinate associated with SIGE(l)
!@var  temp0  actual temperature (K) at (i,j) and SIG(l)
!@var  temp1  actual temperature (K) at (i,j) and SIG(l+1)
!@var  temp1e actual temperature (K) at (i,j) and SIGE(l+1)
!@var  tsavg(i,j) COMPOSITE SURFACE AIR TEMPERATURE (K)

c
c     Grids:
c
c                   -------------------------  lm+1
c                lm  - - - - - - - - - - - - -
c                   -------------------------  lm
c                l+1 - - - - - - - - - - - - -
c                    -------------------------  l+1
C     (main)     l   - - - - - - - - - - - - -
c                    -------------------------  l     (edge)
c                l-1 - - - - - - - - - - - - -
c                    -------------------------  l-1
c                2   - - - - - - - - - - - - -
c                    -------------------------    2
c                1   - - - - - - - - - - - - -
c                    -------------------------    1
c           rhoe(i,j,l+1)=100.d0*(pl-pl1)/(grav*dz(i,j,l))
c           rho(i,j,l)=100.d0*(ple-pl1e)/(grav*dzedge(i,j,l))
c
c     at main: u,v,tv,q,ke
c     at edge: e,lscale,km,kh,gm,gh
c
      USE CONSTANT, only : grav,rgas,kapa
      USE GEOM, only : imaxj
      USE DYNAMICS, only : gz
      implicit none

      integer, intent(in) :: ls1,im,jm,lm
      real*8, intent(in) :: ptop,psf,rvx
      real*8, dimension(im,jm,lm), intent(in) :: tv
      real*8, dimension(im,jm), intent(in) :: p,tsavg,qsavg
      real*8, dimension(lm), intent(in) :: sig
      real*8, dimension(lm+1), intent(in) :: sige
      real*8, dimension(im,jm,lm), intent(out) :: rho,rhoe
      real*8, dimension(im,jm,lm-1), intent(out) :: dz,dzedge

      real*8 :: temp0,temp1,temp1e,pl1,pl,pl1e,ple,pijgcm
      integer :: i,j,l  !@var i,j,l loop variable
      integer :: imax
      real*8 :: plm1e

      do l=1,lm-1
        do j=1,jm
          do i=1,imaxj(j)

            pijgcm=p(i,j)
            if (l.ge.ls1) then
              pl1 =(psf-ptop)*sig(l+1)+ptop
              pl  =(psf-ptop)*sig(l)  +ptop
              pl1e=(psf-ptop)*sige(l+1)+ptop
              ple =(psf-ptop)*sige(l)  +ptop
            else
              pl1 =pijgcm*sig(l+1)+ptop
              pl  =pijgcm*sig(l)  +ptop
              pl1e=pijgcm*sige(l+1)+ptop
              ple =pijgcm*sige(l)  +ptop
            endif
            temp0 =tv(i,j,l)*pl**kapa
            temp1 =tv(i,j,l+1)*pl1**kapa
            temp1e=((sig (l+1)-sige(l+1))*temp0+
     2              (sige(l+1)-sig (l))  *temp1)/
     3              (sig(l+1)-sig(l))
            dz(i,j,l)    =-(rgas/grav)*temp1e*log(pl1/pl)
            dzedge(i,j,l)=-(rgas/grav)*temp0 *log(pl1e/ple)
            rhoe(i,j,l+1)=100.d0*(pl-pl1)/(grav*dz(i,j,l))
            rho(i,j,l)=100.d0*(ple-pl1e)/(grav*dzedge(i,j,l))

            if(l.eq.1) then
              rhoe(i,j,1)=100.d0*ple/(tsavg(i,j)*
     2                    (1.d0+RVX*qsavg(i,j))*rgas)
c             rhoe(i,j,1)=2.d0*rho(i,j,1)-rhoe(i,j,2)       
            endif
            if(l.eq.lm-1) then
c             rho(i,j,lm)=100.d0*pl1/(temp1*rgas)
              if (l.ge.ls1) then
                plm1e=(psf-ptop)*sige(lm+1)+ptop
              else
                plm1e=pijgcm*sige(lm+1)+ptop
              endif
              rho(i,j,lm)=100.d0*(pl1e-plm1e)/(grav*dzedge(i,j,lm))
            endif
          end do
        end do
      end do
      return
      end subroutine getdz

      subroutine dout(u,v,t,pres,prese,q,e,dz,dzedge,
     2                rho,rhoe,u0,v0,t0,q0,
     3                km,kh,ke,gm,gh,lscale,ps,reserv,tsurf,
     4                uflx,vflx,tflx,qflx,itest,jtest,n)
!@sum writes out diagnostics at i=itest, j=jtest
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of potential temperature T
!@var q z-profle of potential temperature Q
!@var e z-profle of turbulent kinetic energy
!@var u0 z-profle of u at previous time step
!@var v0 z-profle of v at previous time step
!@var t0 z-profle of t at previous time step
!@var q0 z-profle of q at previous time step
!@var pres 1-d pressure at z
!@var prese 1-d pressure at zedge
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var lscale z-profile of turbulent length scale
!@var dz(i) z(i+1)-z(i)
!@var dzedge(i) zedge(i+1)-zedge(i)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var ps surface pressure
!@var reserv reserved for furture use
!@var tsurf surface temperature
!@var uflx momentun flux -uw at surface, zedge(1) 
!@var vflx momentun flux -vw at surface, zedge(1) 
!@var tflx heat flux -wt at surface, zedge(1) 
!@var qflx moisture flux -wq at surface, zedge(1) 
!@var itest,jtest (i,j) location at which the output is wriiten
!@var dtime time step
!@var n number of vertical main layers

      USE CONSTANT, only : grav,rgas
      implicit none

      integer, intent(in) :: n,itest,jtest
      real*8, dimension(n), intent(in) :: u,v,t,pres,prese,q,e
      real*8, dimension(n), intent(in) :: rho,rhoe,u0,v0,t0,q0
      real*8, dimension(n), intent(in) :: km,kh,ke,gm,gh,lscale
      real*8, dimension(n-1), intent(in) :: dz,dzedge
      real*8, intent(in) :: ps,reserv,tsurf
      real*8, intent(in) :: uflx,vflx,tflx,qflx

      real*8 :: z,utotal,du,dv,dt,dq,zedge,qturb,dmdz,dtdz,dqdz
      real*8 :: galpha,an2,dudz,dvdz,as2,uf,hf,qf
      real*8 :: ri,rif,sigmat,reserv2
      integer :: i,j,l  !@var i,j,l loop variable

c     Fields on main vertical grid:

      z=(rgas/grav)*0.5d0*(tsurf+t(1))*log(ps/pres(1))+10.d0
      Write (67,1000) itest,jtest,reserv,ps,uflx,vflx,tflx,qflx
      do l=1,n-1
        utotal=sqrt(u(l)*u(l)+v(l)*v(l))
c       du=u(l)-u0(l)
c       dv=v(l)-v0(l)
        dt=t(l)-t0(l)
c       dq=q(l)-q0(l)
        write (67,2000) l,z,pres(l),dz(l),u(l),v(l),t(l),q(l),utotal,
     2                    rho(l),rhoe(l),dt
        z=z+dz(l)
      end do
      utotal=sqrt(u(n)*u(n)+v(n)*v(n))
c       du=u(n)-u0(n)
c       dv=v(n)-v0(n)
        dt=t(n)-t0(n)
c       dq=q(n)-q0(n)
      write (67,2500) l,z,pres(n),u(n),v(n),t(n),q(n),utotal,
     2                  rho(n),rhoe(n),dt
      write (67,9000)
c
c Fields on secondary vertical grid:
c
      write (67,3000)
      l=1
      zedge=10.d0
      qturb=sqrt(2.d0*e(l))
      dq=q(1)-q0(1)
      write (67,2100) l,zedge,prese(l),dzedge(l),dq,
     2                lscale(l),e(l),qturb
      zedge=10.d0+dzedge(1)
      do l=2,n-1
        qturb=sqrt(2.d0*e(l))
        dq=q(l)-q0(l)
        write (67,2000) l,zedge,prese(l),dzedge(l),km(l),kh(l),dq,
     2                  gm(l),gh(l),lscale(l),e(l),qturb
        zedge=zedge+dzedge(l)
      end do
      qturb=sqrt(2.d0*e(n))
      dq=q(n)-q0(n)
      write (67,2500) n,zedge,prese(n),km(n),kh(n),dq,
     2                gm(n),gh(n),lscale(n),e(n),qturb
      write (67,9000)
c
c Fluxes on the secondary grid:
c
      zedge=10.d0
      write (67,4000)
      do l=2,n
        zedge=zedge+dzedge(l-1)
        dtdz=(t(l)-t(l-1))/dz(l-1)
        galpha=grav*2.d0/(t(l)+t(l-1))
        dqdz =(q(l)-q(l-1))/dz(l-1)
        an2=galpha*dtdz
        dudz=(u(l)-u(l-1))/dz(l-1)
        dvdz=(v(l)-v(l-1))/dz(l-1)
        as2=dudz*dudz+dvdz*dvdz
        dmdz=sqrt(as2)
        uf=km(l)*dmdz
        hf=kh(l)*dtdz
        qf=kh(l)*dqdz
        ri=an2/as2
        sigmat=km(l)/kh(l)
        rif=ri/sigmat
        reserv2=0.d0
        write (67,2000) l,zedge,uf,hf,qf,dmdz,dtdz,dqdz,
     2                  ri,rif,sigmat,reserv2
      end do
c
      write (67,9000)
      write (67,9000)
      return
1000  format (1h ,'i    = ',9x,i2,/,1x,
     2            'j    = ',9x,i2,/,1x,
     3            'reserv = ',1pe11.4,/,1x,
     3            'ps   = ',1pe11.4,/,1x,
     3            'uflx = ',1pe11.4,/,1x,
     3            'vflx = ',1pe11.4,/,1x,
     3            'tflx = ',1pe11.4,/,1x,
     3            'qflx = ',1pe11.4,//,1x,' l',1x,
     4            '     z     ',1x,
     5            '     p     ',1x,'     dz    ',1x,'     u     ',1x,
     6            '     v     ',1x,'     t     ',1x,'     q     ',1x,
     7            '   utotal  ',1x,'    rho   ',1x,'   rhoe   ',1x,
     8            '     dt    ',/)
2000  format (1h ,i2,1x,1pe11.4,9(1x,1pe11.4),1x,1pe10.3)
2100  format (1h ,i2,3(1x,1pe11.4),24x,1x,1pe11.4,24x,2(1x,1pe11.4),
     2                 1x,1pe10.3)
2500  format (1h ,i2,2(1x,1pe11.4),12x,7(1x,1pe11.4),1x,1pe10.3)
3000  format (1h ,' l',1x,
     2            '  z (edge) ',1x,
     2            '  p (edge) ',1x,' dz (edge) ',1x,'     km    ',1x,
     3            '     kh    ',1x,'     dq    ',1x,'     gm    ',1x,
     4            '     gh    ',1x,'   lscale  ',1x,'     e     ',1x,
     5            '    vturb  ',/)
4000  format (1h ,' l',1x,
     2            '  z (edge) ',1x,
     3            '    uflux  ',1x,'    hflux  ',1x,'    qflux  ',1x,
     4            '     dmdz  ',1x,'    dtdz   ',1x,'    dqdz   ',1x,
     5            '     Ri    ',1x,'     Rif   ',1x,'   Sigmat  ',1x,
     6            '  reserve  ',/)
9000  format (1h )
      end subroutine dout

      subroutine diff_uv(u0,v0,u,v,km,dz,dzedge,rho,rhoe,
     2                   uflx,vflx,dtime,n)
!@sum integrates differential eqns for u and v (tridiagonal method)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var u0 z-profle of u at previous time step
!@var v0 z-profle of v at previous time step
!@var km z-profile of turbulent viscosity 
!@var kh z-profile of turbulent conductivity
!@var dz(i) z(i+1)-z(i)
!@var dzedge(i) zedge(i+1)-zedge(i)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var uflx momentun flux -uw at surface, zedge(1) 
!@var vflx momentun flux -vw at surface, zedge(1) 
!@var dtime time step
!@var n number of vertical main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: u0,v0,km,rho,rhoe
      real*8, dimension(n), intent(inout) :: u,v
      real*8, dimension(n-1), intent(in) :: dz,dzedge
      real*8, intent(in) :: uflx,vflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1
      real*8 :: alpha
      integer :: j  !@var j loop variable
c
c     sub(j)*u_jm1_kp1+dia(j)*u_j_kp1+sup(j)*u_jp1_kp1 = rhs(j)
c     sub(j)*v_jm1_kp1+dia(j)*v_j_kp1+sup(j)*v_jp1_kp1 = rhs1(j)
c     note: the j on the leftmost refers to the primary grid
c       dxi/dz(j)=dxi/(zedge(j+1)-zedge(j))==dxi/dzedge(j)
c       dxi/dz(j-1/2)=dxi/(z(j)-z(j-1))==dxi/dz(j-1)
c       dxi/dz(j+1/2)=dxi/(z(j+1)-z(j))==dxi/dz(j)
c       p1(j-1/2)=km(j)
c       p1(j+1/2)=km(j+1)
c       similarly in subroutine diff_tq
c
      do j=2,n-1
          sub(j)=-dtime*km(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
          sup(j)=-dtime*km(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=u0(j)
          rhs1(j)=v0(j)
      end do
c
c     Lower boundary conditions:
c
c     d/dt U = -d/dz uw where
c     d/dt U = (u(1)-u0(1))/dtime
c     d/dz uw = (uw(2)-uw(1))/dze(1), dze(1)=ze(2)-ze(1)
c     uw(2)=-km(2)*(u(2)-u(1))/dz(1), dz(1)=z(2)-z(1)
c     uw(1)=-uflx
c     if U at first gcm layer has been updated, then uflx=0
c     this is for U, similarly for V
c     in addition, rho(1) and rhoe(2) are in place to balance the
c     mass
c     from the above, the following follow
c
      alpha=dtime*km(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=u0(1)
      rhs1(1)=v0(1)
c     rhs(1)=u0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*uflx
c     rhs1(1)=v0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*vflx
c
c     Upper boundary conditions:
c
c     d/dt U = -d/dz uw where
c     d/dt U = (u(n)-u0(n))/dtime
c     d/dz uw = (uw(n+1)-uw(n))/dze(n), dze(n)=ze(n+1)-ze(n)
c     uw(n)=-km(n)*(u(n)-u(n-1))/dz(n-1), dz(n-1)=z(n)-z(n-1)
c     uw(n+1)=0 
c
      alpha=dtime*km(n)/(dzedge(n)*dz(n-1)*rho(n))*rhoe(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=u0(n)
      rhs1(n)=v0(n)
c
      call tridiag(sub,dia,sup,rhs,u,n)
      call tridiag(sub,dia,sup,rhs1,v,n)
c
      return
      end subroutine diff_uv

      subroutine diff_tq(tq0,tq,khq,dz,dzedge,rho,rhoe,sflx,dtime,n)
!@sum integrates differential eqns for t and q (tridiagonal method)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var tq z-profle of potential temperature T or relative humidity Q
!@var tq0 z-profle of T or Q at previous time step
!@var khq z-profile of turbulent diffusivity Kh or Kq
!@var dz(i) z(i+1)-z(i)
!@var dzedge(i) zedge(i+1)-zedge(i)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var sflx heat flux -wt or humidity flux -wq at surface, zedge(1) 
!@var dtime time step
!@var n number of vertical main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: tq0,khq,rho,rhoe
      real*8, dimension(n), intent(inout) :: tq
      real*8, dimension(n-1), intent(in) :: dz,dzedge
      real*8, intent(in) :: sflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable
c
c     tq = t or q; khq = kh or kq
c     sub(j)*tq_jm1_kp1+dia(j)*tq_j_kp1+sup(j)*tq_jp1_kp1 = rhs(j)
c
      do j=2,n-1
          sub(j)=-dtime*khq(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
          sup(j)=-dtime*khq(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=tq0(j)
      end do
c
c     Lower boundary conditions:
c     d/dt T = -d/dz wt where
c     d/dt T = (T(1)-T0(1))/dtime
c     d/dz wt = (wt(2)-wt(1))/dze(1), dze(1)=ze(2)-ze(1)
c     wt(2)=-kh(2)*(T(2)-T(1))/dz(1), dz(1)=z(2)-z(1)
c     wt(1)=-tflx
c     if T at first gcm layer has been updated, then tflx=0
c     this is for T, similarly for Q
c     in addition, rho(1) and rhoe(2) are in place to balance the
c     mass
c     from the above, the following follow
c
      alpha=dtime*khq(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=tq0(1)
c     rhs(1)=tq0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*sflx
c
c     Upper boundary conditions:
c
c     d/dt T = -d/dz wt where
c     d/dt T = (T(n)-T0(n))/dtime
c     d/dz wt = (wt(n+1)-wt(n))/dze(n), dze(n)=ze(n+1)-ze(n)
c     wt(n)=-kh(n)*(T(n)-T(n-1))/dz(n-1), dz(n-1)=z(n)-z(n-1)
c     wt(n+1)=0 
c
      alpha=dtime*khq(n)/(dzedge(n)*dz(n-1)*rho(n))*rhoe(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=tq0(n)
c
      call tridiag(sub,dia,sup,rhs,tq,n)
c
      return
      end subroutine diff_tq

      subroutine diff_e(e0,e,km,kh,ke,lscale,u,v,t,
     2           dz,dzedge,rho,rhoe,dtime,ustar2,n)
!@sum integrates differential eqns for e (tridiagonal method)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var e z-profle of turbulent kinetic energy
!@var e0 z-profle of e at previous time step
!@var km z-profile of turbulent diffusivity
!@var kh z-profile of turbulent conductivity
!@var ke z-profile of turbulent diffusivity for e-equation
!@var lscale z-profile of turbulent length scale
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of potential temperature
!@var lscale z-profile of turbulent length scale
!@var dz(i) z(i+1)-z(i)
!@var dzedge(i) zedge(i+1)-zedge(i)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var grav gravitational acceleration
!@var ustar2 sqrt(uw^2+vw^2)
!@var dtime time step
!@var n number of vertical main layers
      USE CONSTANT, only : grav
      USE SOCPBL, only : b1,b123
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: e0,km,kh,ke,lscale,
     &        u,v,t,rho,rhoe
      real*8, dimension(n), intent(inout) :: e
      real*8, dimension(n-1), intent(in) :: dz,dzedge
      real*8, intent(in) :: dtime,ustar2

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: qturb,tmp,an2,dudz,dvdz,as2
      integer :: j  !@var j loop variable

      real*8, parameter :: emin=5.d-5,emax=2.d0
c
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     note: the j on the leftmost refers to the secondary grid
c       dxi/dz(j)=dxi/(z(j)-z(j-1))==dxi/dz(j-1)
c       dxi/dz(j-1/2)=dxi/(zedge(j)-zedge(j-1))==dxi/dzedge(j-1)
c       dxi/dz(j+1/2)=dxi/(zedge(j+1)-zedge(j))==dxi/dzedge(j)
c       p1(j-1/2)=0.5d0*(ke(j)+ke(j-1))
c       p1(j+1/2)=0.5d0*(ke(j)+ke(j+1))
c       ke(j)=sq*lscale(j)*qturb, qturb=sqrt(2.d0*e(j))
c
      do j=2,n-1
          qturb=sqrt(2.d0*e(j))
          tmp=rho(j-1)/rhoe(j)
          sub(j)=-dtime*0.5d0*(ke(j)+ke(j-1))/(dz(j-1)*dzedge(j-1))*tmp
          tmp=rho(j)/rhoe(j)
          sup(j)=-dtime*0.5d0*(ke(j)+ke(j+1))/(dz(j-1)*dzedge(j))*tmp
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*2*qturb/(b1*lscale(j))
          an2=grav*2.d0/(t(j)+t(j-1))*(t(j)-t(j-1))/dz(j-1)
          dudz=(u(j)-u(j-1))/dz(j-1)
          dvdz=(v(j)-v(j-1))/dz(j-1)
          as2=dudz*dudz+dvdz*dvdz
          rhs(j)=e0(j)+dtime*(km(j)*as2-kh(j)*an2)
      end do
c
c     Boundary conditions:
c
      dia(1)=1.d0
      sup(1)=0.d0
      rhs(1)=0.5d0*b123*ustar2
c
      sub(n)=-1.d0
      dia(n)=1.d0
      rhs(n)=0.d0
c
      call tridiag(sub,dia,sup,rhs,e,n)
c
      do j=1,n
         if(e(j).lt.emin) e(j)=emin
          if(e(j).gt.emax) e(j)=emax
      end do

      return
      end subroutine diff_e

      subroutine lgcm(lscale,u,v,t,e,dz,dzedge,rho,n)
!@sum lgcm calculates the z-profle of the turbulent length scale
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var lscale z-profle of turbulent length scale
!@var u z-profle of due-east wind component
!@var v z-profle of due-north wind component
!@var t z-profle of potential temperature
!@var e z-profle of pturbulent kinetic energy
!@var dz(i) z(i+1)-z(i)
!@var dzedge(i) zedge(i+1)-zedge(i)
!@var rho the z-profile of the density 
!@var n number of GCM layers
!@var zgs height of surface layer (m), imported from SOCPBL
      USE CONSTANT, only : grav
      USE SOCPBL, only : kappa,zgs

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: u,v,t,e,rho
      real*8, dimension(n), intent(out) :: lscale
      real*8, dimension(n-1), intent(in) :: dz,dzedge

      real*8, dimension(n) :: zedge
      real*8, parameter :: alpha0=0.2d0,emin=5.d-5,emax=2.d0
      real*8 :: dudz,dvdz,as2,lmax2
      real*8 :: sum1,sum2,qi,qim1,l0,l1,kz,an2,lmax
      integer :: i  !@var i loop variable
 
      zedge(1)=zgs
      do i=2,n
          zedge(i)=zedge(i-1)+dzedge(i-1)
      end do
 
c     integration of monotonically tabulated function by
c     trapezoidal rule
 
      sum1=0.d0
      sum2=0.d0
      do i=2,n
        qi=sqrt(2.d0*e(i))
        qim1=sqrt(2.d0*e(i-1))
        sum1=sum1+.5d0*dzedge(i-1)*(qi+qim1)*rho(i-1)
        sum2=sum2+.5d0*dzedge(i-1)*(qi*zedge(i)+qim1*zedge(i-1))
     &           *rho(i-1)
      end do
      l0=alpha0*sum2/sum1

      kz=kappa*zedge(1)
      lscale(1)=l0*kz/(l0+kz)
      if (lscale(1).gt.dzedge(1)) lscale(1)=dzedge(1)
      do i=2,n-1
        kz=kappa*zedge(i)
        l1=l0*kz/(l0+kz)
        if (t(i).gt.t(i-1)) then
          an2=grav*2.d0/(t(i)+t(i-1))*(t(i)-t(i-1))/dz(i-1)
          dudz=(u(i)-u(i-1))/dz(i-1)
          dvdz=(v(i)-v(i-1))/dz(i-1)
          as2=dudz*dudz+dvdz*dvdz
          lmax  =0.53d0*sqrt(2.d0*e(i)/(an2+1.d-40))
          lmax2 =1.95d0*sqrt(2.d0*e(i)/(as2+1.d-40))
          lmax=min(lmax,lmax2)
          if (l1.gt.lmax) l1=lmax
        endif
        if (l1.gt.dzedge(i)) l1=dzedge(i)
        lscale(i)=l1
      end do
      lscale(n)=lscale(n-1)

      return
      end subroutine lgcm


      subroutine kgcm(km,kh,ke,gm,gh,u,v,t,e,lscale,dz,dbll,n)
c
c     Grids:
c
c                j+1 - - - - - - - - - - - - -
c                    -------------------------  j+1
c     (main)     j   - - - - - - - - - - - - -
c                    -------------------------  j     (edge)
c                j-1 - - - - - - - - - - - - -
c                    -------------------------  j-1
c                2   - - - - - - - - - - - - -
c                    -------------------------    2
c                1   - - - - - - - - - - - - -
c                    -------------------------    1
c
c     dz(j)==z(j+1)-z(j), dzedge(j)==ze(j+1)-ze(j)
c     at main: u,v,t,q,ke
c     at edge: e,lscale,km,kh,gm,gh
!@sum kgcm computes the turbulent stability functions Km, Kh and Ke
!@sum using the GISS second order closure model (2000)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u,v,t,e,lscale z-profiles
!@var dz(i) z(i+1)-z(i)
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var dbll number of layers where the PBL height corresponds to (real*8)
!@var n number of GCM layers
!@var tau B1*lscale/sqrt(2*e) 
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@car an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var sq stability constant for e, adjustable
      USE CONSTANT, only : grav
      USE SOCPBL, only : ghmin,ghmax,gmmax0,d0,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s3,s4,s5,s6,b1

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,e,lscale
      real*8, dimension(n-1), intent(in) :: dz
      real*8, dimension(n), intent(out) :: km,kh,ke,gm,gh
      real*8, intent(out) :: dbll

      ! 0.2*lscale*qturb \approx 0.02*e*tau, because b1/2 \approx 10
      real*8, parameter ::  sq=0.02d0
      real*8 :: an2,dudz,dvdz,as2,ell,den,qturb,tau,ghi,gmi,gmmax
      real*8 :: sm,sh,taue,e_lpbl,e_main_i
      integer :: i  !@var i loop variable

      do i=2,n
        an2=grav*2.d0/(t(i)+t(i-1))*(t(i)-t(i-1))/dz(i-1)
        dudz=(u(i)-u(i-1))/dz(i-1)
        dvdz=(v(i)-v(i-1))/dz(i-1)
        as2=dudz*dudz+dvdz*dvdz
        ell=lscale(i)
        qturb=sqrt(2.d0*e(i))
        tau=B1*ell/qturb
        ghi=tau*tau*an2
        gmi=tau*tau*as2
        if(ghi.lt.ghmin) ghi=ghmin
        if(ghi.gt.ghmax) ghi=ghmax
        gmmax=(d0+d1*ghi+d3*ghi*ghi)/(d2+d4*ghi)
        gmmax=min(gmmax,gmmax0)
        if(gmi.gt.gmmax) gmi=gmmax
        den=d0+d1*ghi+d2*gmi+d3*ghi*ghi+d4*ghi*gmi+d5*gmi*gmi
        sm=(s0+s1*ghi+s2*gmi)/den
        sh=(s4+s5*ghi+s6*gmi)/den
        taue=tau*e(i)
        km(i)=min(max(taue*sm,1.5d-5),100.d0)
        kh(i)=min(max(taue*sh,2.5d-5),100.d0)
        ke(i)=min(max(taue*sq,1.5d-5),100.d0)
        gm(i)=gmi
        gh(i)=ghi
      end do
      ke(1)=b1/sqrt(2.d0)*lscale(1)*sqrt(e(1))*sq
      do i=1,n-1
        ke(i)=0.5d0*(ke(i)+ke(i+1))
      end do

c     find the pbl top (at main level lpbl)

      e_lpbl=0.1d0*e(1) ! if e(main i) < e_lpbl, i is the pbl top
      do i=1,n-1
        e_main_i = 0.5d0*(e(i)+e(i+1))
        dbll=i        ! dbll is real*8
        if (e_main_i.lt.e_lpbl) exit
      end do

      return
      end subroutine kgcm

      subroutine ave_uv_to_tcell(u,v,u_tcell,v_tcell,im,jm,lm)
!@sum Computes u_tcell,v_tcell from u,v
!@sum u_tcell is the average of 4 nearest u
!@sum v_tcell is the average of 4 nearest v
!@auth  Ye Cheng
!@ver   1.0
!@var u west-east   velocity component at secondary grids (ucell)
!@var v south-north velocity component at secondary grids (ucell)
!@var u_tcell u at primary grids (tcell)
!@var v_tcell v at primary grids (tcell)
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm,lm), intent(in) :: u,v
      real*8, dimension(im,jm,lm), intent(out) :: u_tcell,v_tcell

      integer :: im1,iq1,iq2,iq3
      integer :: i,j,l

      do j=2,jm-1
        do i=1,im
          if(i.eq.1) then
             im1=im
          else
             im1=i-1
          endif

          do l=1,lm
            u_tcell(i,j,l)=0.25d0*(u(im1,j+1,l)+u(i,j+1,l)
     2                          +u(im1,j,l)+u(i,j,l))
            v_tcell(i,j,l)=0.25d0*(v(im1,j+1,l)+v(i,j+1,l)
     2                          +v(im1,j,l)+v(i,j,l))
          end do
        end do
      end do

      ! for j=1 (south pole) and j=jm (north pole)
      iq1=nint(0.25d0*im)+1
      iq2=nint(0.50d0*im)+1
      iq3=nint(0.75d0*im)+1
      do l=1,lm
            u_tcell(1,1,l)=0.25d0*(u(1,2,l)-u(iq2,2,l)
     2                            +v(iq1,2,l)-v(iq3,2,l))
            u_tcell(1,jm,l)=0.25d0*(u(1,jm,l)-u(iq2,jm,l)
     2                             -v(iq1,jm,l)+v(iq3,jm,l))
            v_tcell(1,1,l)=0.25d0*(v(1,2,l)-v(iq2,2,l)
     2                            -u(iq1,2,l)+u(iq3,2,l))
            v_tcell(1,jm,l)=0.25d0*(v(1,jm,l)-v(iq2,jm,l)
     2                             +u(iq1,jm,l)-u(iq3,jm,l))
      end do

      return
      end subroutine ave_uv_to_tcell

      subroutine ave_to_ucell(t,t_ucell,im,jm,lm)
!@sum Computes t_ucell from t
!@sum t_ucell is the average of 4 nearest t
!@var t scalar at primary grids (tcell)
!@var t_ucell t at secondary grids (ucell)
!@aut  Ye Cheng/G. Hartke
!@ver  1.0
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm,lm), intent(in) :: t
      real*8, dimension(im,jm,lm), intent(out) :: t_ucell

      integer :: ip1
      integer :: i,j,l
      ! ucell is from j=2 to j=jm, no ucell on the poles
      do j=3,jm-1
        do i=1,im
          if(i.eq.im) then
             ip1=1
          else
             ip1=i+1
          endif
          do l=1,lm
            t_ucell(i,j,l)=0.25d0*(t(i,j,l)+t(ip1,j,l)
     2                         + t(i,j-1,l)+t(ip1,j-1,l))
          end do
        end do
      end do
      ! for j = 2 and j = jm
        do i=1,im
          if(i.eq.im) then
             ip1=1
          else
             ip1=i+1
          endif
          do l=1,lm
            t_ucell(i,2,l)=0.333333333333333d0*(t(1,1,l)
     2               +t(i,2,l)+t(ip1,2,l))
            t_ucell(i,jm,l)=0.333333333333333d0*(t(1,jm,l)
     2               +t(i,jm-1,l)+t(ip1,jm-1,l))
          end do
        end do

      return
      end subroutine ave_to_ucell
