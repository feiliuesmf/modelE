      subroutine diffus(dtime)
!@sum  diffus updates u,v,t,q due to 
!@+  turbulent transport throughout all GCM layers
!@+  using a second order closure (SOC)
!@+  turbulence model developed at GISS, 2000.
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
      USE DYNAMICS, only : pmid,pk,pedn
      USE MODEL_COM, only :
     *      im,jm,lm,sig,sige,u,v,t,q,p
     *     ,vt_on      
      USE CONSTANT, only : kapa,deltx
      USE PBLCOM, only : tsavg,qsavg,dclev,uflux,vflux,tflux,qflux,egcm
      USE GEOM, only : imaxj

      IMPLICIT NONE

      real*8, intent(in) :: dtime

      real*8, dimension(lm) :: uij,vij,tij,pij,qij,eij
      real*8, dimension(lm) :: u0ij,v0ij,t0ij,q0ij,e0ij
      real*8, dimension(lm) :: rhoebydz,bydzerho
      real*8, dimension(lm) :: km,kh,ke,lscale,gm,gh,rhoij,rhoeij
      real*8, dimension(lm) :: dzij,dzeij
      real*8, dimension(lm,im,jm) :: rho,rhoe
      real*8, dimension(lm,im,jm) :: dz,dzedge
      real*8, dimension(lm) :: peij,ttest
      real*8, dimension(lm,im,jm) :: u_tcell,v_tcell,tv_ucell,t_virtual
      real*8, dimension(lm,im,jm) :: km_gcm,km_gcm_ucell
     2        ,dz_ucell,dzedge_ucell,rho_ucell,rhoe_ucell
      real*8, dimension(im,jm) :: p_ucell,tsavg_ucell,qsavg_ucell
      real*8, dimension(im,jm) :: uflux_ucell,vflux_ucell

      integer, parameter :: mout=0,itest=52,jtest=33,itmax=2
      real*8, parameter :: tol=1.d-4,qmin=1.d-20,p00=1000.d0
      integer, save :: ifirst=1
      real*8, save :: p1000k=1.d0
      real*8 :: uflx,vflx,tflx,qflx,pijgcm,pl,rvx
      real*8 :: temp0,ustar2,dbll,reserv,test,check
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
            t_virtual(l,i,j)=t(i,j,l)*(1.d0+RVX*Q(i,j,l))
          end do
        end do
      end do

c     integrate T,Q equations at tcells

      ! get u_tcell and v_tcell at t-cells
      call ave_uv_to_tcell1(u,v,u_tcell,v_tcell,im,jm,lm)

      call getdz(t_virtual,p,dz,dzedge,rho,rhoe
     2           ,tsavg,qsavg,rvx,im,jm,lm)

      loop_j_tq: do j=1,jm
        loop_i_tq: do i=1,imaxj(j)
          do l=1,lm
            uij(l)=u_tcell(l,i,j)
            vij(l)=v_tcell(l,i,j)
            tij(l)=t_virtual(l,i,j)*p1000k  !virtual,potential temp.
            pij(l)=100.d0*pmid(l,i,j)
            qij(l)=q(i,j,l)
            eij(l)=egcm(l,i,j)
            rhoeij(l)=rhoe(l,i,j)
            rhoij(l)=rho(l,i,j)
            u0ij(l)=uij(l)
            v0ij(l)=vij(l)
            t0ij(l)=tij(l)
            q0ij(l)=qij(l)
            e0ij(l)=eij(l)
          end do
          do l=1,lm-1
            dzij(l)=dz(l,i,j)
            rhoebydz(l+1)=rhoeij(l+1)/dzij(l)
          end do
          do l=1,lm
            dzeij(l)=dzedge(l,i,j)
            bydzerho(l)=1.d0/(dzeij(l)*rhoij(l))
          end do

          uflx  =uflux(i,j)
          vflx  =vflux(i,j)
          ustar2=sqrt(uflx*uflx+vflx*vflx)
          tflx  =tflux(i,j)
          qflx  =qflux(i,j)

          call lgcm(lscale,uij,vij,tij,eij,dzij,dzeij,rhoij,lm)
          call kgcm(km,kh,ke,gm,gh,uij,vij,tij,eij,lscale,dzij,lm)
          call diff_e(e0ij,eij,km,kh,ke,lscale,uij,vij,tij,dzij,dzeij
     2         ,rhoij,rhoeij,dtime,ustar2,lm)
          call diff_tq(t0ij,tij,kh,dzij,dzeij,
     2                 rhoij,rhoeij,rhoebydz,bydzerho,tflx,dtime,lm)
          call diff_tq(q0ij,qij,kh,dzij,dzeij,
     2                 rhoij,rhoeij,rhoebydz,bydzerho,qflx,dtime,lm)

          do 300 iter=1,itmax

            do l=2,lm-1
              ttest(l)=tij(l)
            end do

            call lgcm(lscale,uij,vij,tij,eij,dzij,dzeij,rhoij,lm)
            call kgcm(km,kh,ke,gm,gh,uij,vij,tij,eij,lscale,dzij,lm)
            call diff_e(e0ij,eij,km,kh,ke,lscale,uij,vij,tij,dzij,dzeij
     2           ,rhoij,rhoeij,dtime,ustar2,lm)
            call diff_tq(t0ij,tij,kh,dzij,dzeij
     2           ,rhoij,rhoeij,rhoebydz,bydzerho,tflx,dtime,lm)
            call diff_tq(q0ij,qij,kh,dzij,dzeij
     2           ,rhoij,rhoeij,rhoebydz,bydzerho,qflx,dtime,lm)
            call find_pbl_top(eij,dbll,lm)

            test=0.d0
            icount=0
            do l=2,lm-1
              check=abs(ttest(l)-tij(l))
              if (check.gt.0.d0) then
                test=test+check/abs(ttest(l)+tij(l))
                icount=icount+1
              endif
            end do
            test=test/(float(icount)+1.d-40)
            if (test.lt.tol) exit

300       continue

          do l=1,lm
            t_virtual(l,i,j)=tij(l)/p1000k
            q(i,j,l)=max(qij(l),qmin)
            t(i,j,l)=t_virtual(l,i,j)/(1.d0+RVX*Q(i,j,l))
            egcm(l,i,j)=eij(l)
            km_gcm(l,i,j)=km(l)
          end do
          dclev(i,j)=dbll
c
c         Write out diagnostics if at the correct grid point:
c
          if (mout.eq.1.and.(i.eq.itest).and.(j.eq.jtest)) then
            do l=1,lm
                peij(l)=100.d0*pedn(l,i,j)
            end do
            call dout(uij,vij,tij,pij,peij,qij,eij,dzij,dzeij,
     2                rhoij,rhoeij,u0ij,v0ij,t0ij,q0ij,
     3                km,kh,ke,gm,gh,lscale,reserv,tsavg(i,j),
     4                uflx,vflx,tflx,qflx,itest,jtest,lm)
          endif

        end do loop_i_tq
      end do loop_j_tq

c     integrate U,V equations at ucells

      ! average some quantities at u-cells
c     call ave_to_ucell(p,p_ucell,im,jm,1)
c     call ave_ufvf_to_ucell(uflux,vflux,uflux_ucell,vflux_ucell,im,jm)
c     call ave_to_ucell(uflux,uflux_ucell,im,jm,1)
c     call ave_to_ucell(vflux,vflux_ucell,im,jm,1)
c     call ave_to_ucell(tsavg,tsavg_ucell,im,jm,1)
c     call ave_to_ucell(qsavg,qsavg_ucell,im,jm,1)

c     call ave_to_ucell(t_virtual,tv_ucell,im,jm,lm)

      call ave_to_ucell1(km_gcm,km_gcm_ucell,im,jm,lm)
      call ave_to_ucell1(dz,dz_ucell,im,jm,lm)
      call ave_to_ucell1(dzedge,dzedge_ucell,im,jm,lm)
      call ave_to_ucell1(rho,rho_ucell,im,jm,lm)
      call ave_to_ucell1(rhoe,rhoe_ucell,im,jm,lm)

      loop_j_uv: do j=2,jm
        loop_i_uv: do i=1,im
          do l=1,lm
            uij(l)=u(i,j,l)
            vij(l)=v(i,j,l)
            rhoeij(l)=rhoe_ucell(l,i,j)
            rhoij(l)=rho_ucell(l,i,j)
            u0ij(l)=uij(l)
            v0ij(l)=vij(l)
          end do
          do l=1,lm-1
            dzij(l)=dz_ucell(l,i,j)
            rhoebydz(l+1)=rhoeij(l+1)/dzij(l)
            km(l)=km_gcm_ucell(l,i,j)
          end do
          do l=1,lm
            dzeij(l)=dzedge_ucell(l,i,j)
            bydzerho(l)=1.d0/(dzeij(l)*rhoij(l))
          end do

c         uflx  =uflux_ucell(i,j)
c         vflx  =vflux_ucell(i,j)

          call diff_uv(u0ij,v0ij,uij,vij,km,dzij,dzeij,rhoij,rhoeij
     2                   ,rhoebydz,bydzerho,uflx,vflx,dtime,lm)

          do l=1,lm
            u(i,j,l)= uij(l)
            v(i,j,l)= vij(l)
          end do
 
        end do loop_i_uv 
      end do loop_j_uv

      return
      end subroutine diffus

      subroutine getdz(tv,p,dz,dzedge,rho,rhoe
     2          ,tsavg,qsavg,rvx,im,jm,lm)
!@sum  getdz computes the 3d finite difference dz and dzedge
!@+    as well as the 3d density rho and rhoe
!@+    called at the primary cells (t-cells)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var  dz(l,i,j) z(l+1,i,j) - z(l,i,j)
!@var  dzedge(l,i,j) zedge(l+1,i,j) - zedge(l,i,j)
!@var  z vertical coordinate associated with SIG(l)
!@var  zedge vertical coordinate associated with SIGE(l)
!@var  temp0 virtual temperature (K) at (i,j) and SIG(l)
!@var  temp1 virtual temperature (K) at (i,j) and SIG(l+1)
!@var  temp1e average of temp0 and temp1
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
c           rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
c           rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dzedge(l,i,j))
c
c     at main: u,v,tv,q,ke
c     at edge: e,lscale,km,kh,gm,gh
c
      USE CONSTANT, only : grav,rgas,kapa
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pmid,pk,pedn
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, intent(in) :: rvx
      real*8, dimension(lm,im,jm), intent(in) :: tv
      real*8, dimension(im,jm), intent(in) :: p,tsavg,qsavg
      real*8, dimension(lm,im,jm), intent(out) :: rho,rhoe
      real*8, dimension(lm,im,jm), intent(out) :: dz,dzedge

      real*8 :: temp0,temp1,temp1e,pl1,pl,pl1e,ple,pijgcm
      integer :: i,j,l  !@var i,j,l loop variable
      integer :: imax
      real*8 :: plm1e

      do j=1,jm
        do i=1,imaxj(j)
          do l=1,lm-1

            pl1 =pmid(l+1,i,j)
            pl  =pmid(l,i,j)
            pl1e=pedn(l+1,i,j)
            ple =pedn(l,i,j)
            temp0 =tv(l,i,j)*pk(l,i,j)
            temp1 =tv(l+1,i,j)*pk(l+1,i,j)
            temp1e=0.5d0*(temp0+temp1)
            dz(l,i,j)    =-(rgas/grav)*temp1e*log(pl1/pl)
            dzedge(l,i,j)=-(rgas/grav)*temp0 *log(pl1e/ple)
            rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
            rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dzedge(l,i,j))
            if(l.eq.1) then
              rhoe(1,i,j)=100.d0*ple/(tsavg(i,j)*
     2                    (1.d0+RVX*qsavg(i,j))*rgas)
c             rhoe(1,i,j)=2.d0*rho(1,i,j)-rhoe(2,i,j)       
            endif
            if(l.eq.lm-1) then
c             rho(lm,i,j)=100.d0*pl1/(temp1*rgas)
              plm1e=pedn(lm+1,i,j)
              dzedge(lm,i,j)=-(rgas/grav)*temp1 *log(plm1e/pl1e)
              rho(lm,i,j)=100.d0*(pl1e-plm1e)/(grav*dzedge(lm,i,j))
            endif
          end do
        end do
      end do
      return
      end subroutine getdz

      subroutine dout(u,v,t,pres,prese,q,e,dz,dzedge,
     2                rho,rhoe,u0,v0,t0,q0,
     3                km,kh,ke,gm,gh,lscale,reserv,tsurf,
     4                uflx,vflx,tflx,qflx,itest,jtest,n)
!@sum dout writes out diagnostics at i=itest, j=jtest
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
      real*8, dimension(n), intent(in) :: dz,dzedge
      real*8, intent(in) :: reserv,tsurf
      real*8, intent(in) :: uflx,vflx,tflx,qflx

      real*8 :: z,utotal,du,dv,dt,dq,zedge,qturb,dmdz,dtdz,dqdz
      real*8 :: galpha,an2,dudz,dvdz,as2,uf,hf,qf
      real*8 :: ri,rif,sigmat,reserv2,ps
      integer :: i,j,l  !@var i,j,l loop variable

c     Fields on main vertical grid:
      ps=prese(1)
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

      subroutine diff_uv(u0,v0,u,v,km,dz,dzedge,rho,rhoe
     2                   ,rhoebydz,bydzerho,uflx,vflx,dtime,n)
!@sum diff_uv integrates differential eqns for u and v (tridiagonal method)
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
     2        ,rhoebydz,bydzerho      
      real*8, dimension(n), intent(inout) :: u,v
      real*8, dimension(n), intent(in) :: dz,dzedge
      real*8, intent(in) :: uflx,vflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1
      real*8 :: alpha
      integer :: j  !@var j loop variable
c
c     sub(j)*u_jm1_kp1+dia(j)*u_j_kp1+sup(j)*u_jp1_kp1 = rhs(j)
c     sub(j)*v_jm1_kp1+dia(j)*v_j_kp1+sup(j)*v_jp1_kp1 = rhs1(j)
c     note: j refers to the layer middle
c     except for km(j), which is defined on the layer edge
c     similarly in subroutine diff_tq
c
      do j=2,n-1
c         sub(j)=-dtime*km(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
c         sup(j)=-dtime*km(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          sub(j)=-dtime*km(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*km(j+1)*rhoebydz(j+1)*bydzerho(j)
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
c     alpha=dtime*km(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
      alpha=dtime*km(2)*rhoebydz(2)*bydzerho(1)
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
c     alpha=dtime*km(n)/(dzedge(n)*dz(n-1)*rho(n))*rhoe(n)
      alpha=dtime*km(n)*rhoebydz(n)*bydzerho(n)
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

      subroutine diff_tq(tq0,tq,khq,dz,dzedge,rho,rhoe
     2                   ,rhoebydz,bydzerho,sflx,dtime,n)
!@sum diff_tq integrates differential eqns for t and q (tridiagonal method)
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
     2        ,rhoebydz,bydzerho      
      real*8, dimension(n), intent(inout) :: tq
      real*8, dimension(n), intent(in) :: dz,dzedge
      real*8, intent(in) :: sflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable
c
c     tq = t or q; khq = kh or kq
c     sub(j)*tq_jm1_kp1+dia(j)*tq_j_kp1+sup(j)*tq_jp1_kp1 = rhs(j)
c
      do j=2,n-1
c         sub(j)=-dtime*khq(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
c         sup(j)=-dtime*khq(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          sub(j)=-dtime*khq(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*khq(j+1)*rhoebydz(j+1)*bydzerho(j)
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
c     alpha=dtime*khq(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
      alpha=dtime*khq(2)*rhoebydz(2)*bydzerho(1)
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
c     alpha=dtime*khq(n)/(dzedge(n)*dz(n-1)*rho(n))*rhoe(n)
      alpha=dtime*khq(n)*rhoebydz(n)*bydzerho(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=tq0(n)
c
      call tridiag(sub,dia,sup,rhs,tq,n)
c
      return
      end subroutine diff_tq

      subroutine diff_e(e0,e,km,kh,ke,lscale,u,v,t,dz,dzedge
     2          ,rho,rhoe,dtime,ustar2,n)
!@sum diff_e integrates differential eqns for e (tridiagonal method)
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
      real*8, dimension(n), intent(in) :: dz,dzedge
      real*8, intent(in) :: dtime,ustar2

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: qturb,tmp,an2,dudz,dvdz,as2
      integer :: j  !@var j loop variable

      real*8, parameter :: emin=1.d-20,emax=10.d0
c
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     j refers to the layer edge
c     except ke(j) which is defined on the layer middle
      do j=2,n-1
          qturb=sqrt(2.d0*e(j))
          tmp=rho(j-1)/rhoe(j)
          sub(j)=-dtime*ke(j-1)/(dz(j-1)*dzedge(j-1))*tmp
          tmp=rho(j)/rhoe(j)
          sup(j)=-dtime*ke(j)/(dz(j-1)*dzedge(j))*tmp
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
      real*8, dimension(n), intent(in) :: dz,dzedge

      real*8, dimension(n) :: zedge
      real*8, parameter :: alpha0=0.1d0
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
      do i=2,n
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

      return
      end subroutine lgcm


      subroutine kgcm(km,kh,ke,gm,gh,u,v,t,e,lscale,dz,n)
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
!@+   using the GISS second order closure model (2000)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u,v,t,e,lscale z-profiles
!@var dz(i) z(i+1)-z(i)
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var n number of GCM layers
!@var tau B1*lscale/sqrt(2*e) 
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@var an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var sq stability constant for e, adjustable
      USE CONSTANT, only : grav
      USE SOCPBL, only : ghmin,ghmax,gmmax0,d0,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s3,s4,s5,s6,b1

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,e,lscale
      real*8, dimension(n), intent(in) :: dz
      real*8, dimension(n), intent(out) :: km,kh,ke,gm,gh

      ! note e *tau = b1/2 *lscale * qturb
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
      ke(1)=b1*lscale(1)*sqrt(0.5d0*e(1))*sq
      do i=1,n-1
        ke(i)=0.5d0*(ke(i)+ke(i+1)) !defined on main levels
      end do
      return
      end subroutine kgcm

      subroutine find_pbl_top(e,dbll,n)
!@sum find_pbl_top Find the pbl top (at main level lpbl)
!@+   if e(i) le certain fraction of e(1), real(i) is pbl top
!@auth  Ye Cheng
!@ver   1.0
!@var dbll the (real*8) layer number corresponds to the top of the pbl
      real*8, dimension(n), intent(in) :: e
      real*8, intent(out) :: dbll

      real*8, parameter :: fraction = 0.1d0
      real*8 :: e1p    ! certain percent of e_1
      integer i

      e1p=fraction*e(1)
      do i=2,n
        if (e(i).lt.e1p) exit
      end do
      dbll=i   ! dbll is real*8
      return
      end subroutine find_pbl_top

      subroutine ave_uv_to_tcell(u,v,u_tcell,v_tcell,im,jm,lm)
!@sum ave_uv_to_tcell Computes u_tcell,v_tcell from u,v, where u and v
!@+   may be the x and y components of a vector defined at secondary grids 
!@+   Note that u_tcell,v_tcell are of dimension (lm,im,jm)
!@var u an x-component at secondary grids (ucell)
!@var v a  y-component at secondary grids (ucell)
!@var u_tcell an x-component at primary grids (tcell)
!@var v_tcell a  y-component at primary grids (tcell)
!@auth Ye Cheng
!@ver  1.0
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm,lm), intent(in) ::  u,v
      real*8, dimension(lm,im,jm), intent(out) :: u_tcell,v_tcell

      real*8 :: HEMI,u_t,v_t
      integer :: i,j,l,k,IMAX,KMAX

c     polar boxes
      DO J=1,JM,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        HEMI=1.
        IF(J.LE.JM/2) HEMI=-1.
        DO I=1,IMAX
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              u_t=u_t+rapj(k,j)*(u(idij(k,i,j),idjj(k,j),L)*cosiv(k)-
     2                      hemi*v(idij(k,i,j),idjj(k,j),L)*siniv(k))
              v_t=v_t+rapj(k,j)*(v(idij(k,i,j),idjj(k,j),L)*cosiv(k)+
     2                      hemi*u(idij(k,i,j),idjj(k,j),L)*siniv(k))
            END DO
            u_tcell(l,i,j)=u_t
            v_tcell(l,i,j)=v_t
          END DO
        END DO
      END DO
c     non polar boxes
      DO J=2,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              u_t=u_t+u(idij(k,i,j),idjj(k,j),L)*rapj(k,j)
              v_t=v_t+v(idij(k,i,j),idjj(k,j),L)*rapj(k,j)
            END DO
            u_tcell(l,i,j)=u_t
            v_tcell(l,i,j)=v_t
          END DO
        END DO
      END DO
C****
      return
      end subroutine ave_uv_to_tcell

      subroutine ave_ufvf_to_ucell(uf,vf,uf_ucell,vf_ucell,im,jm)
!@sum ave_ufvf_to_ucell Computes uf_ucell,vf_ucell from uf,vf where uf and vf
!@+   may be the x and y components of a vector defined at primary grids 
!@var uf an x-component at primary grids (tcell)
!@var vf a  y-component at primary grids (tcell)
!@var uf_ucell an x-component at secondary grids (ucell)
!@var vf_ucell an y-component at secondary grids (ucell)
!@auth Ye Cheng
!@ver  1.0
      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj,cosiv,siniv
      implicit none

      integer, intent(in) :: im,jm
      real*8, dimension(im,jm), intent(in) ::  uf,vf
      real*8, dimension(im,jm), intent(out) :: uf_ucell,vf_ucell

      real*8 :: HEMI
      integer :: i,j,l,k,IMAX,KMAX

      uf_ucell=0.d0; vf_ucell=0.d0
c     polar boxes
      DO J=1,JM,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        HEMI=1.
        IF(J.LE.JM/2) HEMI=-1.
        DO I=1,IMAX
        DO K=1,KMAX
          uf_ucell(IDIJ(K,I,J),IDJJ(K,J))=
     *    uf_ucell(IDIJ(K,I,J),IDJJ(K,J)) -
     *      RAVJ(K,J)*(uf(I,J)*COSIV(K)+vf(I,J)*SINIV(K)*HEMI)
          vf_ucell(IDIJ(K,I,J),IDJJ(K,J))=
     *    vf_ucell(IDIJ(K,I,J),IDJJ(K,J)) -
     *      RAVJ(K,J)*(vf(I,J)*COSIV(K)-uf(I,J)*SINIV(K)*HEMI)
        END DO
        END DO
      END DO
c     non polar boxes
      DO J=2,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
        DO K=1,KMAX
          uf_ucell(IDIJ(K,I,J),IDJJ(K,J))=
     *    uf_ucell(IDIJ(K,I,J),IDJJ(K,J)) -
     *           RAVJ(K,J)*uf(I,J)
          vf_ucell(IDIJ(K,I,J),IDJJ(K,J))=
     *    vf_ucell(IDIJ(K,I,J),IDJJ(K,J)) -
     *           RAVJ(K,J)*vf(I,J)
        END DO
        END DO
      END DO
C****

      return
      end subroutine ave_ufvf_to_ucell

      subroutine ave_to_ucell(s,s_ucell,im,jm,lm)
!@sum ave_to_ucell Computes s_ucell from s
!@+   Note that all arrays here are of dimension (lm,im,jm)
!@var s a scalar at primary grids (tcell)
!@var s_ucell a scalar at secondary grids (ucell)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) :: s
      real*8, dimension(lm,im,jm), intent(out) :: s_ucell

      real*8 :: s_u
      integer :: i,j,l,k,ip1

      ! ucell is from j=2 to j=jm, no ucell on the poles
      DO J=3,JM-1
        DO I=1,IMAXJ(J)
          DO L=1,LM
            s_u=0.d0
            DO K=1,KMAXJ(J)
              s_u=s_u+s(l,idij(k,i,j),idjj(k,j))*ravj(k,j)
            END DO
            s_ucell(l,i,j)=s_u
          END DO
        END DO
      END DO
      ! for j = 2 and j = jm
        do i=1,im
          if(i.eq.im) then
             ip1=1
          else
             ip1=i+1
          endif
          do l=1,lm
            s_ucell(l,i,2)=0.333333333333333d0*(s(l,1,1)
     2               +s(l,i,2)+s(l,ip1,2))
            s_ucell(l,i,jm)=0.333333333333333d0*(s(l,1,jm)
     2               +s(l,i,jm-1)+s(l,ip1,jm-1))
          end do
        end do

      return
      end subroutine ave_to_ucell

ccccccccccccccccccccccccccccc

      subroutine ave_uv_to_tcell1(u,v,u_tcell,v_tcell,im,jm,lm)
!@sum ave_uv_to_tcell1 Computes u_tcell,v_tcell from u,v
!@+   u_tcell is the average of 4 nearest u
!@+   v_tcell is the average of 4 nearest v
!@auth  Ye Cheng
!@ver   1.0
!@var u west-east   velocity component at secondary grids (ucell)
!@var v south-north velocity component at secondary grids (ucell)
!@var u_tcell u at primary grids (tcell)
!@var v_tcell v at primary grids (tcell)
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm,lm), intent(in) :: u,v
      real*8, dimension(lm,im,jm), intent(out) :: u_tcell,v_tcell

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
            u_tcell(l,i,j)=0.25d0*(u(im1,j+1,l)+u(i,j+1,l)
     2                          +u(im1,j,l)+u(i,j,l))
            v_tcell(l,i,j)=0.25d0*(v(im1,j+1,l)+v(i,j+1,l)
     2                          +v(im1,j,l)+v(i,j,l))
          end do
        end do
      end do

      ! for j=1 (south pole) and j=jm (north pole)
      iq1=nint(0.25d0*im)+1
      iq2=nint(0.50d0*im)+1
      iq3=nint(0.75d0*im)+1
      do l=1,lm
            u_tcell(l,1,1)=0.25d0*(u(1,2,l)-u(iq2,2,l)
     2                            +v(iq1,2,l)-v(iq3,2,l))
            u_tcell(l,1,jm)=0.25d0*(u(1,jm,l)-u(iq2,jm,l)
     2                             -v(iq1,jm,l)+v(iq3,jm,l))
            v_tcell(l,1,1)=0.25d0*(v(1,2,l)-v(iq2,2,l)
     2                            -u(iq1,2,l)+u(iq3,2,l))
            v_tcell(l,1,jm)=0.25d0*(v(1,jm,l)-v(iq2,jm,l)
     2                             +u(iq1,jm,l)-u(iq3,jm,l))
      end do

      return
      end subroutine ave_uv_to_tcell1

      subroutine ave_to_ucell1(s,s_ucell,im,jm,lm)
!@sum ave_to_ucell1 Computes s_ucell from s
!@+   s_ucell is the average of 4 nearest s
!@var s scalar at primary grids (tcell)
!@var s_ucell s at secondary grids (ucell)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) :: s
      real*8, dimension(lm,im,jm), intent(out) :: s_ucell

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
            s_ucell(l,i,j)=0.25d0*(s(l,i,j)+s(l,ip1,j)
     2                         + s(l,i,j-1)+s(l,ip1,j-1))
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
            s_ucell(l,i,2)=0.333333333333333d0*(s(l,1,1)
     2               +s(l,i,2)+s(l,ip1,2))
            s_ucell(l,i,jm)=0.333333333333333d0*(s(l,1,jm)
     2               +s(l,i,jm-1)+s(l,ip1,jm-1))
          end do
        end do

      return
      end subroutine ave_to_ucell1
