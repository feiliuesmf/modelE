      subroutine diffus(lbase_min,lbase_max,dtime)
!@sum  diffus updates u,v,t,q due to
!@+  turbulent transport throughout all GCM layers
!@+  using a second order closure (SOC)
!@+  turbulence model developed at GISS, 2000.
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from diffB347D6M20)
!@cont diffus,getdz,dout,diff_uv,diff_t,diff_q,diff_e
!@cont lgcm,kgcm,ave_uv_to_tcell,ave_to_ucell
!@var u 3d west-east wind component
!@var v 3d south-north wind component
!@var t 3d potential temperature
!@var q 3d relative humidity
!@var p 2-d pressure
!@var dtime time step
!@var lbase_min/max levels through which to apply turbulence (dummy)
!@var mout 0:don't call dout; 1:call dout
!@var itest longitude at which to call dout
!@var jtest latitude at which to call dout
      USE DYNAMICS, only : pmid,pk,pedn,pdsig,plij
      USE MODEL_COM, only :
     *      im,jm,lm,sig,sige,u,v,t,q,p,vt_on
      USE CONSTANT, only : grav,kapa,deltx,lhe,sha
      USE PBLCOM, only : tsavg,qsavg,dclev,uflux,vflux,tflux,qflux,egcm
     *     ,t2gcm
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj
      USE DAGCOM, only : ajl,
     &     jl_trbhr,jl_damdc,jl_trbke,jl_trbdlht
      USE DAGPCOM, only : p1000k
      USE SOCPBL, only : nlevel

      IMPLICIT NONE

      integer, intent(in) :: lbase_min,lbase_max
      real*8, intent(in) :: dtime

      real*8, dimension(lm) :: uij,vij,tij,pij,qij,eij,t2ij
      real*8, dimension(lm) :: u0ij,v0ij,t0ij,q0ij,e0ij,t20ij
      real*8, dimension(lm) :: dudz,dvdz,dtdz,dqdz,g_alpha,as2,an2
      real*8, dimension(lm) :: rhoebydz,bydzerho,gc_t2
      real*8, dimension(lm) :: km,kh,kq,ke,kt2,gc,lscale,gm,gh
      real*8, dimension(lm) :: rhoij,rhoeij,dzij,dzeij
      real*8, dimension(im,jm,lm) :: uold
      real*8, dimension(lm,im,jm) :: rho,rhoe
      real*8, dimension(lm,im,jm) :: dz,dzedge
      real*8, dimension(lm) :: peij
      real*8, dimension(lm,im,jm) :: u_tcell,v_tcell,tv_ucell,t_virtual
      real*8, dimension(lm,im,jm) :: km_gcm,km_gcm_ucell
     2        ,dz_ucell,dzedge_ucell,rho_ucell,rhoe_ucell
      real*8, dimension(im,jm) :: p_ucell,tsavg_ucell,qsavg_ucell
      real*8, dimension(im,jm) :: uflux_ucell,vflux_ucell

      integer, parameter :: mout=0,itest=52,jtest=33
      real*8, parameter :: tol=1.d-4,qmin=1.d-20
      real*8 :: uflx,vflx,tflx,qflx,pl,rvx
      real*8 :: temp0,ustar2,dbll,reserv,test,check,t0ijl,rak
      integer :: loc,icount,imax,kmax,idik,idjk
      integer :: i,j,l,k,iter,ifirst !@i,j,l,k loop variable
      data ifirst/0/

C**** Note that lbase_min/max are here for backwards compatibility with
C**** original drycnv. They are only used to determine where the
C**** routine has been called from.
      if (lbase_min.eq.2) return       ! quit if called from main

      if(.not. vt_on) then
          rvx=0.d0
      else
          rvx=deltx
      endif

      !  convert input T to virtual T
      do j=1,jm
        do i=1,imaxj(j)
          do l=1,lm
            t_virtual(l,i,j)=t(i,j,l)*(1.d0+RVX*Q(i,j,l))
          end do
        end do
      end do

      ! integrate T,Q equations at tcells

      ! get u_tcell and v_tcell at t-cells
      call ave_uv_to_tcell1(u,v,u_tcell,v_tcell,im,jm,lm)

      call getdz(t_virtual,p,dz,dzedge,rho,rhoe,
     2           tsavg,qsavg,rvx,im,jm,lm)
 
      if(ifirst.eq.0) then
        ifirst=1
        do j=1,jm
          do i=1,imaxj(j)
            do l=1,lm
               t2gcm(l,i,j)=egcm(l,i,j)*1.d-3
            end do
          end do
        end do
      endif

      loop_j_tq: do j=1,jm
        loop_i_tq: do i=1,imaxj(j)
          do l=1,lm
            uij(l)=u_tcell(l,i,j)
            vij(l)=v_tcell(l,i,j)
            tij(l)=t_virtual(l,i,j)*p1000k  !virtual,potential temp.
            pij(l)=100.d0*pmid(l,i,j)
            qij(l)=q(i,j,l)
            eij(l)=egcm(l,i,j)
            t2ij(l)=t2gcm(l,i,j)
            rhoeij(l)=rhoe(l,i,j)
            rhoij(l)=rho(l,i,j)
            t0ij(l)=tij(l)
            q0ij(l)=qij(l)
            e0ij(l)=eij(l)
            t20ij(l)=t2ij(l)
          end do
          do l=1,lm-1
            dzij(l)=dz(l,i,j)
            rhoebydz(l+1)=rhoeij(l+1)/dzij(l)
          end do
          do l=1,lm
            dzeij(l)=dzedge(l,i,j)
            bydzerho(l)=1.d0/(dzeij(l)*rhoij(l))
          end do

          ! calculate z-derivatives on the edges of the layers
          do l=2,lm
            dudz(l)=(uij(l)-uij(l-1))/dzij(l-1)
            dvdz(l)=(vij(l)-vij(l-1))/dzij(l-1)
            dtdz(l)=(tij(l)-tij(l-1))/dzij(l-1)
            dqdz(l)=(qij(l)-qij(l-1))/dzij(l-1)
            g_alpha(l)=grav*2.d0/(tij(l)+tij(l-1))
            an2(l)=g_alpha(l)*dtdz(l)
            as2(l)=dudz(l)*dudz(l)+dvdz(l)*dvdz(l)
          end do

          uflx=uflux(i,j)
          vflx=vflux(i,j)
          tflx=tflux(i,j)
          qflx=qflux(i,j)
          ustar2=sqrt(uflx*uflx+vflx*vflx)

          call lgcm(lscale,uij,vij,tij,eij,dudz,dvdz,as2,an2,
     2              dzij,dzeij,rhoij,lm)
c          call init_t2(uij,vij,tij,eij,t2ij,t20ij,
c    2          dudz,dvdz,as2,dtdz,g_alpha,an2,lscale,dzij,lm)
          call kgcm(km,kh,kq,ke,kt2,gc,gc_t2,gm,gh,uij,vij,tij,
     2             eij,t20ij,dudz,dvdz,as2,dtdz,g_alpha,an2,
     3             lscale,dzij,dzeij,0,lm)
          call diff_e(e0ij,eij,km,kh,ke,gc,lscale,uij,vij,tij,
     2         dzij,dzeij,dudz,dvdz,as2,dtdz,g_alpha,an2,
     3         rhoij,rhoeij,ustar2,dtime,lm)
          if(nlevel.eq.3) call diff_t2(t20ij,t2ij,eij,kh,kt2,gc,gc_t2,
     2      lscale,tij,dtdz,dzij,dzeij,rhoij,rhoeij,
     3      ustar2,tflx,dtime,lm)
          call diff_t(t0ij,tij,kh,gc,dzij,dzeij,
     2                 rhoij,rhoeij,rhoebydz,bydzerho,tflx,dtime,lm)
          call diff_q(q0ij,qij,kq,dzij,dzeij,
     2                 rhoij,rhoeij,rhoebydz,bydzerho,qflx,dtime,lm)

          call find_pbl_top(eij,dbll,lm)

          do l=1,lm
            ! update 3-d q,t,egcm,t2gcm and km_gcm
            q(i,j,l)=max(qij(l),qmin)
            t0ijl=t(i,j,l)
            t(i,j,l)=tij(l)/(p1000k*(1.d0+RVX*Q(i,j,l)))
            egcm(l,i,j)=eij(l)
            t2gcm(l,i,j)=t2ij(l)
            km_gcm(l,i,j)=km(l)
            ! ACCUMULATE DIAGNOSTICS for t and q
            AJL(J,L,JL_TRBHR)=AJL(J,L,JL_TRBHR)
     2                 +(T(I,J,L)-t0ijl)*PK(L,I,J)*PLIJ(L,I,J)
            AJL(J,L,JL_TRBDLHT)=AJL(J,L,JL_TRBDLHT)
     2                 +(Q(I,J,L)-q0ij(l))*PDSIG(L,I,J)*LHE/SHA
            AJL(J,L,JL_TRBKE)=AJL(J,L,JL_TRBKE)+eij(l)
          end do
          ! update pbl height (layer number counted from below, real*8)
          dclev(i,j)=dbll

          ! Write out diagnostics if at selected grid point:

          if (mout.eq.1.and.(i.eq.itest).and.(j.eq.jtest)) then
            do l=1,lm
                peij(l)=100.d0*pedn(l,i,j)
            end do
            call dout(uij,vij,tij,pij,peij,qij,eij,t2ij,dzij,dzeij,
     1           dudz,dvdz,as2,dtdz,g_alpha,an2,dqdz,
     2           rhoij,rhoeij,t0ij,q0ij,
     3           km,kh,kq,ke,kt2,gc,gm,gh,lscale,reserv,tsavg(i,j),
     4           uflx,vflx,tflx,qflx,itest,jtest,iter,lm)
          endif

        end do loop_i_tq
      end do loop_j_tq

      ! integrate U,V equations at ucells

      ! average some quantities at u-cells
c     call ave_to_ucell(p,p_ucell,im,jm,1)
c     call ave_to_ucell(tsavg,tsavg_ucell,im,jm,1)
c     call ave_to_ucell(qsavg,qsavg_ucell,im,jm,1)

      call ave_ufvf_to_ucell(uflux,vflux,uflux_ucell,vflux_ucell,im,jm)
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
            uold(i,j,l)=uij(l)
            dzeij(l)=dzedge_ucell(l,i,j)
            bydzerho(l)=1.d0/(dzeij(l)*rhoij(l))
            km(l)=km_gcm_ucell(l,i,j)
          end do
          do l=1,lm-1
            dzij(l)=dz_ucell(l,i,j)
            rhoebydz(l+1)=rhoeij(l+1)/dzij(l)
          end do

          uflx  =uflux_ucell(i,j)
          vflx  =vflux_ucell(i,j)

          call diff_uv(u0ij,v0ij,uij,vij,km,dzij,dzeij,rhoij,rhoeij,
     2                 rhoebydz,bydzerho,uflx,vflx,dtime,lm)

          ! update u and v
          do l=1,lm
            u(i,j,l)= uij(l)
            v(i,j,l)= vij(l)
          end do

        end do loop_i_uv
      end do loop_j_uv

      ! ACCUMULATE DIAGNOSTICS for u and v
      DO J=1,JM
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
           DO L=1,LM
             DO K=1,KMAX
               RAK=RAVJ(K,J)
               IDIK=IDIJ(K,I,J)
               IDJK=IDJJ(K,J)
               AJL(IDJK,L,JL_DAMDC)=AJL(IDJK,L,JL_DAMDC)
     &           +(U(IDIK,IDJK,L)-uold(IDIK,IDJK,L))*PLIJ(L,I,J)*RAK
             ENDDO
          ENDDO
        ENDDO
      ENDDO

      return
      end subroutine diffus

      subroutine getdz(tv,p,dz,dzedge,rho,rhoe,
     2           tsavg,qsavg,rvx,im,jm,lm)
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

      real*8 :: temp0,temp1,temp1e,pl1,pl,pl1e,ple
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

      subroutine dout(u,v,t,pres,prese,q,e,t2,dz,dzedge,
     1                dudz,dvdz,as2,dtdz,g_alpha,an2,dqdz,
     2                rho,rhoe,t0,q0,
     3                km,kh,kq,ke,kt2,gc,gm,gh,lscale,reserv,tsurf,
     4                uflx,vflx,tflx,qflx,itest,jtest,iter,n)
!@sum dout writes out diagnostics at i=itest, j=jtest
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var t z-profle of potential temperature T
!@var q z-profle of potential temperature Q
!@var e z-profle of turbulent kinetic energy
!@var t0 z-profle of t at previous time step
!@var q0 z-profle of q at previous time step
!@var pres 1-d pressure at z
!@var prese 1-d pressure at zedge
!@var km turbulent viscosity for u and v equations
!@var kh turbulent diffusivity for t
!@var kq turbulent diffusivity for q
!@var ke turbulent diffusivity for e equation
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var lscale z-profile of turbulent length scale
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
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

      integer, intent(in) :: n,iter,itest,jtest
      real*8, dimension(n), intent(in) :: u,v,t,pres,prese,q,e,t2
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2,dqdz
      real*8, dimension(n), intent(in) :: rho,rhoe,t0,q0
      real*8, dimension(n), intent(in) :: km,kh,kq,ke,kt2,gc,gm,gh
      real*8, dimension(n), intent(in) :: lscale,dz,dzedge
      real*8, intent(in) :: reserv,tsurf
      real*8, intent(in) :: uflx,vflx,tflx,qflx

      real*8 :: z,utotal,dt,dq,zedge,qturb,dmdz
      real*8 :: uf,hf,qf
      real*8 :: ri,rif,sigmat,reserv2,ps,wt1
      integer :: i,j,l  !@var i,j,l loop variable

c     Fields on main vertical grid:
      ps=prese(1)
      z=(rgas/grav)*0.5d0*(tsurf+t(1))*log(ps/pres(1))+10.d0
      Write (67,1001) "iter=",iter
      Write (67,1000) itest,jtest,reserv,ps,uflx,vflx,tflx,qflx
      do l=1,n-1
        utotal=sqrt(u(l)*u(l)+v(l)*v(l))
        dt=t(l)-t0(l)
c       dq=q(l)-q0(l)
        write (67,2000) l,z,pres(l),dz(l),u(l),v(l),t(l),q(l),utotal,
     2                    rho(l),rhoe(l),dt
        z=z+dz(l)
      end do
      utotal=sqrt(u(n)*u(n)+v(n)*v(n))
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
      write (67,2100) l,zedge,prese(l),
     2                lscale(l),e(l),t2(l)
      zedge=10.d0+dzedge(1)
      do l=2,n
        wt1=-kh(l)*(t(l)-t(l-1))/dz(l-1)
        write (67,2000) l,zedge,prese(l),gc(l),wt1,kh(l),kt2(l),
     2                  gm(l),gh(l),lscale(l),e(l),t2(l)
        zedge=zedge+dzedge(l)
      end do
      write (67,9000)
c
c Fluxes on the secondary grid:
c
      zedge=10.d0
      write (67,4000)
      do l=2,n
        zedge=zedge+dzedge(l-1)
        dmdz=sqrt(as2(l))
        uf=km(l)*dmdz
        hf=kh(l)*dtdz(l)
        qf=kq(l)*dqdz(l)
        ri=an2(l)/as2(l)
        sigmat=km(l)/kh(l)
        rif=ri/sigmat
        reserv2=0.d0
        write (67,2000) l,zedge,uf,hf,qf,dmdz,dtdz(l),dqdz(l),
     2                  ri,rif,sigmat,reserv2
      end do
c
      write (67,9000)
      write (67,9000)
      return
1001  format(a,i4)
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
2100  format (1h ,i2,2(1x,1pe11.4),24x,1x,11x,36x,2(1x,1pe11.4),
     2           ,1x,1pe10.3)
2500  format (1h ,i2,2(1x,1pe11.4),12x,7(1x,1pe11.4),1x,1pe10.3)
3000  format (1h ,' l',1x,
     2            '  z (edge) ',1x,
     2            '  p (edge) ',1x,'     gc    ',1x,' -kh*dtdz  ',1x,
     3            '     kh    ',1x,'    kt2    ',1x,'     gm    ',1x,
     4            '     gh    ',1x,'   lscale  ',1x,'     e     ',1x,
     5            '     t2    ',/)
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
!@sum diff_uv integrates differential eqns for u and v (tridiag. method)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var u z-profle of west-east   velocity component
!@var v z-profle of south-north velocity component
!@var u0 z-profle of u at previous time step
!@var v0 z-profle of v at previous time step
!@var km z-profile of turbulent viscosity
!@var kh z-profile of turbulent conductivity
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var uflx momentun flux -uw at surface, zedge(1)
!@var vflx momentun flux -vw at surface, zedge(1)
!@var dtime time step
!@var n number of vertical main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: u0,v0,km,dz,dzedge
     2        ,rho,rhoe,rhoebydz,bydzerho
      real*8, dimension(n), intent(out) :: u,v
      real*8, intent(in) :: uflx,vflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs,rhs1
      real*8 :: alpha
      integer :: j  !@var j loop variable
c
c     sub(j)*u_jm1_kp1+dia(j)*u_j_kp1+sup(j)*u_jp1_kp1 = rhs(j)
c     sub(j)*v_jm1_kp1+dia(j)*v_j_kp1+sup(j)*v_jp1_kp1 = rhs1(j)
c     note: j refers to the layer middle
c     except for km(j), which is defined on the layer edge
c     similarly in subroutine diff_t, diff_q
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
      alpha=dtime*km(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
c     alpha=dtime*km(2)*rhoebydz(2)*bydzerho(1)
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

      subroutine diff_t(t0,t,kh,gc,dz,dzedge,rho,rhoe
     2                   ,rhoebydz,bydzerho,tflx,dtime,n)
!@sum diff_t integrates differential eqns for t (tridiag. method)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var t z-profle of potential temperature T
!@var t0 z-profle of T at previous time step
!@var kh z-profile of heat diffusivity
!@var gc z-profile of countergradient term in heat flux (level 3)
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var tflx heat flux -wt at surface, zedge(1)
!@var dtime time step
!@var n number of vertical main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: t0,kh,gc,dz,dzedge
     2        ,rho,rhoe,rhoebydz,bydzerho
      real*8, dimension(n), intent(out) :: t
      real*8, intent(in) :: tflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha,p4j,dgcdz
      integer :: j  !@var j loop variable

c     sub(j)*t_jm1_kp1+dia(j)*t_j_kp1+sup(j)*t_jp1_kp1 = rhs(j)

      do j=2,n-1
c         sub(j)=-dtime*kh(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
c         sup(j)=-dtime*kh(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          sub(j)=-dtime*kh(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*kh(j+1)*rhoebydz(j+1)*bydzerho(j)
          dgcdz=(gc(j+1)-gc(j))/dzedge(j)
          p4j=-dgcdz
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=t0(j)+dtime*p4j
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
      alpha=dtime*kh(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
c     alpha=dtime*kh(2)*rhoebydz(2)*bydzerho(1)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
c     rhs(1)=t0(1)
      rhs(1)=t0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*gc(2)
c     rhs(1)=t0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*(gc(2)+tflx)
c
c     Upper boundary conditions:
c
c     d/dt T = -d/dz wt where
c     d/dt T = (T(n)-T0(n))/dtime
c     d/dz wt = (wt(n+1)-wt(n))/dze(n), dze(n)=ze(n+1)-ze(n)
c     wt(n)=-kh(n)*(T(n)-T(n-1))/dz(n-1), dz(n-1)=z(n)-z(n-1)
c     wt(n+1)=0
c
c     alpha=dtime*kh(n)/(dzedge(n)*dz(n-1)*rho(n))*rhoe(n)
      alpha=dtime*kh(n)*rhoebydz(n)*bydzerho(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=t0(n)+dtime*gc(n)/(dzedge(n)*rho(n))*rhoe(n)
c
      call tridiag(sub,dia,sup,rhs,t,n)
c
      return
      end subroutine diff_t

      subroutine diff_q(q0,q,kq,dz,dzedge,rho,rhoe
     2                   ,rhoebydz,bydzerho,qflx,dtime,n)
!@sum diff_q integrates differential eqns for q (tridiag. method)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var q z-profle of relative humidity Q
!@var q0 z-profle of Q at previous time step
!@var kq z-profile of moisture diffusivity
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var qflx heat flux -wt or humidity flux -wq at surface, zedge(1)
!@var dtime time step
!@var n number of vertical main layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: q0,kq,dz,dzedge
     2        ,rho,rhoe,rhoebydz,bydzerho
      real*8, dimension(n), intent(out) :: q
      real*8, intent(in) :: qflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable

c     sub(j)*q_jm1_kp1+dia(j)*q_j_kp1+sup(j)*q_jp1_kp1 = rhs(j)

      do j=2,n-1
c         sub(j)=-dtime*kq(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
c         sup(j)=-dtime*kq(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          sub(j)=-dtime*kq(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*kq(j+1)*rhoebydz(j+1)*bydzerho(j)
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=q0(j)
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
      alpha=dtime*kq(2)/(dzedge(1)*dz(1)*rho(1))*rhoe(2)
c     alpha=dtime*kq(2)*rhoebydz(2)*bydzerho(1)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=q0(1)
c     rhs(1)=q0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*qflx
c
c     Upper boundary conditions:
c
c     d/dt T = -d/dz wt where
c     d/dt T = (T(n)-T0(n))/dtime
c     d/dz wt = (wt(n+1)-wt(n))/dze(n), dze(n)=ze(n+1)-ze(n)
c     wt(n)=-kh(n)*(T(n)-T(n-1))/dz(n-1), dz(n-1)=z(n)-z(n-1)
c     wt(n+1)=0
c
c     alpha=dtime*kq(n)/(dzedge(n)*dz(n-1)*rho(n))*rhoe(n)
      alpha=dtime*kq(n)*rhoebydz(n)*bydzerho(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=q0(n)
c
      call tridiag(sub,dia,sup,rhs,q,n)
c
      return
      end subroutine diff_q

      subroutine diff_e(e0,e,km,kh,ke,gc,lscale,u,v,t,dz,dzedge,
     1           dudz,dvdz,as2,dtdz,g_alpha,an2,
     2           rho,rhoe,ustar2,dtime,n)
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
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
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
      real*8, dimension(n), intent(in) :: e0,km,kh,ke,gc,lscale
      real*8, dimension(n), intent(in) :: u,v,t,dz,dzedge,rho,rhoe
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2
      real*8, dimension(n), intent(out) :: e
      real*8, intent(in) :: ustar2,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: qturb,tmp
      integer :: j  !@var j loop variable

      real*8, parameter :: emin=1.d-20,emax=100.d0
c
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     j refers to the layer edge
c     except ke(j) which is defined on the layer middle
      do j=2,n-1
          qturb=sqrt(2.d0*e0(j))
          tmp=rho(j-1)/rhoe(j)
          sub(j)=-dtime*ke(j-1)/(dz(j-1)*dzedge(j-1))*tmp
          tmp=rho(j)/rhoe(j)
          sup(j)=-dtime*ke(j)/(dz(j-1)*dzedge(j))*tmp
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*2*qturb/(b1*lscale(j))
          rhs(j)=e0(j)+dtime*(km(j)*as2(j)-kh(j)*an2(j)
     2           +g_alpha(j)*gc(j))
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

      subroutine diff_t2(t20,t2,e,kh,kt2,gc,gc_t2,lscale,t,dtdz,
     2                   dz,dzedge,rho,rhoe,ustar2,tflx,dtime,n)
!@sum diff_t2 integrates differential eqns for t2 (tridiagonal method)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var t2 z-profle of turbulent temperature varince
!@var t20 z-profle of t2 at previous time step
!@var kh z-profile of turbulent conductivity
!@var kt2 z-profile of turbulent diffusivity for t2-equation
!@var lscale z-profile of turbulent length scale
!@var t z-profle of potential temperature
!@var lscale z-profile of turbulent length scale
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var grav gravitational acceleration
!@var ustar2 sqrt(uw^2+vw^2)
!@var tflx turbulent heat flux
!@var dtime time step
!@var n number of vertical main layers
      USE CONSTANT, only : grav
      USE SOCPBL, only : b1,b123,b2,prt
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: t20,e,kh,kt2,gc,lscale
      real*8, dimension(n), intent(in) :: gc_t2,t,rho,dz,dzedge,rhoe
      real*8, dimension(n), intent(in) :: dtdz
      real*8, dimension(n), intent(out) :: t2
      real*8, intent(in) :: ustar2,tflx,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: qturb,tmp,p3j,p4j
      integer :: j  !@var j loop variable

      real*8, parameter :: t2min=1.d-20,t2max=10.d0
c
c     sub(j)*t2_jm1_kp1+dia(j)*t2_j_kp1+sup(j)*t2_jp1_kp1 = rhs(j)
c     j refers to the layer edge
c     except ke(j) which is defined on the layer middle
      do j=2,n-1
          qturb=sqrt(2.d0*e(j))
          tmp=rho(j-1)/rhoe(j)
          sub(j)=-dtime*kt2(j-1)/(dz(j-1)*dzedge(j-1))*tmp
          tmp=rho(j)/rhoe(j)
          sup(j)=-dtime*kt2(j)/(dz(j-1)*dzedge(j))*tmp
c         p3j=2.d0*(qturb/(b2*lscale(j))+gc(j)/t20(j)*dtdz(j))
          p3j=2.d0*(qturb/(b2*lscale(j))+gc_t2(j)*dtdz(j))
c         p3j=max(p3j,1.d-20)
c         p3j=2.d0*(qturb/(b2*lscale(j)))
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*p3j
c         p4j=2.d0*(kh(j)*dtdz(j)-gc(j))*dtdz(j)
          p4j=2.d0*(kh(j)*dtdz(j)*dtdz(j))
          rhs(j)=t20(j)+dtime*p4j
      end do
c
c     Boundary conditions:
c
      dia(1)=1.d0
      sup(1)=0.d0
      rhs(1)=tflx*tflx*b2*prt/(ustar2*b1**(1./3.))
C
      sub(n)=-1.d0
      dia(n)=1.d0
      rhs(n)=0.d0
c
      call tridiag(sub,dia,sup,rhs,t2,n)
c
      do j=1,n
         if(t2(j).lt.t2min) t2(j)=t2min
         if(t2(j).gt.t2max) t2(j)=t2max
      end do

      return
      end subroutine diff_t2

      subroutine lgcm(lscale,u,v,t,e,dudz,dvdz,as2,an2,dz,dzedge,rho,n)
!@sum lgcm calculates the z-profle of the turbulent length scale
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var lscale z-profle of turbulent length scale
!@var u z-profle of due-east wind component
!@var v z-profle of due-north wind component
!@var t z-profle of potential temperature
!@var e z-profle of pturbulent kinetic energy
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
!@var rho the z-profile of the density
!@var n number of GCM layers
!@var zgs height of surface layer (m), imported from SOCPBL
      USE CONSTANT, only : grav
      USE SOCPBL, only : kappa,zgs

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: u,v,t,e,rho
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2,an2
      real*8, dimension(n), intent(out) :: lscale
      real*8, dimension(n), intent(in) :: dz,dzedge

      real*8, dimension(n) :: zedge
      real*8, parameter :: alpha0=0.1d0
      real*8 :: lmax2
      real*8 :: sum1,sum2,qj,qjm1,l0,l1,kz,lmax
      integer :: j  !@var j loop variable

      zedge(1)=zgs
      do j=2,n
          zedge(j)=zedge(j-1)+dzedge(j-1)
      end do

c     integration of monotonically tabulated function by
c     trapezoidal rule

      sum1=0.d0
      sum2=0.d0
      do j=2,n
        qj=sqrt(2.d0*e(j))
        qjm1=sqrt(2.d0*e(j-1))
        sum1=sum1+.5d0*dzedge(j-1)*(qj+qjm1)*rho(j-1)
        sum2=sum2+.5d0*dzedge(j-1)*(qj*zedge(j)+qjm1*zedge(j-1))
     &           *rho(j-1)
      end do
      l0=alpha0*sum2/sum1

      kz=kappa*zedge(1)
      lscale(1)=l0*kz/(l0+kz)
      if (lscale(1).gt.dzedge(1)) lscale(1)=dzedge(1)
      do j=2,n
        kz=kappa*zedge(j)
        l1=l0*kz/(l0+kz)
        if (t(j).gt.t(j-1)) then
          lmax  =0.53d0*sqrt(2.d0*e(j)/(an2(j)+1.d-40))
          lmax2 =1.95d0*sqrt(2.d0*e(j)/(as2(j)+1.d-40))
          lmax=min(lmax,lmax2)
          if (l1.gt.lmax) l1=lmax
        endif
        if (l1.gt.dzedge(j)) l1=dzedge(j)
        lscale(j)=l1
      end do

      return
      end subroutine lgcm

      subroutine level2(u,v,t,e,t2,t20,dudz,dvdz,as2,
     2                  dtdz,g_alpha,an2,lscale,dz,n)
      USE CONSTANT, only : grav
      USE SOCPBL, only : rimax,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,b1,b2
     *     ,c1,c2,c3,c4,c5
      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,lscale,dz
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2
      real*8, dimension(n), intent(out) :: e,t2,t20

      real*8 :: ri,ell,qturb,ghj,gmj
      real*8 :: den,sm,sh
      real*8 :: aa,bb,cc,tmp
      integer :: j  !@var i loop variable

      do j=2,n
        ri=an2(j)/max(as2(j),1.d-20)
        if(ri.gt.rimax) ri=rimax
        aa=c1*ri*ri-c2*ri+c3
        bb=c4*ri+c5
        cc=2.
        if(abs(aa).lt.1.e-8) then
            gmj= -cc/bb
        else
            tmp=bb*bb-4.*aa*cc
            gmj=(-bb-sqrt(tmp))/(2.*aa)
        endif
        ghj=ri*gmj
        ell=lscale(j)
        e(j)=0.5d0*(B1*ell)**2*as2(j)/(gmj+1.e-20)
        if(e(j).lt.1.d-20) e(j)=1.d-20
        qturb=sqrt(2.d0*e(j))
        den=1+d1*ghj+d2*gmj+d3*ghj*ghj+d4*ghj*gmj+d5*gmj*gmj
        sm=(s0+s1*ghj+s2*gmj)/den
        sh=(s4+s5*ghj+s6*gmj)/den
        t2(j)=0.5d0*b1*b2*sh*(ell*dtdz(j))**2
        t20(j)=t2(j)
      enddo

      return
      end subroutine level2

      subroutine init_t2(u,v,t,e,t2,t20,dudz,dvdz,as2,
     2                   dtdz,g_alpha,an2,lscale,dz,n)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u,v,t,e,lscale z-profiles
!@var dz(j) z(j+1)-z(j)
!@var n number of GCM layers
!@var tau B1*lscale/sqrt(2*e)
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@var an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var se stability constant for e, adjustable
!@var st2 stability constant for t2, adjustable
      USE CONSTANT, only : grav
      USE SOCPBL, only : ghmin,ghmax,gmmax0,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,b1,b2

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,e,lscale,dz
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2
      real*8, dimension(n), intent(out) :: t2,t20

      real*8 :: ell,den,qturb,tau,ghj,gmj,gmmax
      real*8 :: sm,sh,taue,e_lpbl,t21,tmp
      integer :: j  !@var j loop variable

      do j=2,n
        ell=lscale(j)
        qturb=sqrt(2.d0*e(j))
        tau=B1*ell/qturb
        ghj=tau*tau*an2(j)
        gmj=tau*tau*as2(j)
        if(ghj.lt.ghmin) ghj=ghmin
        if(ghj.gt.ghmax) ghj=ghmax
        gmmax=(1+d1*ghj+d3*ghj*ghj)/(d2+d4*ghj)
        gmmax=min(gmmax,gmmax0)
        if(gmj.gt.gmmax) gmj=gmmax
        den=1+d1*ghj+d2*gmj+d3*ghj*ghj+d4*ghj*gmj+d5*gmj*gmj
        sh=(s4+s5*ghj+s6*gmj)/den
        t2(j)=0.5d0*b1*b2*sh*(ell*dtdz(j))**2
        t20(j)=t2(j)
      end do

      return
      end subroutine init_t2 

      subroutine kgcm(km,kh,kq,ke,kt2,gc,gc_t2,
     2                gm,gh,u,v,t,e,t2,dudz,dvdz,as2,
     3                dtdz,g_alpha,an2,lscale,dz,dze,iter,n)
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
!@var dz(j) z(j+1)-z(j)
!@var km turbulent viscosity for u and v equations
!@var kh turbulent conductivity for t and q equations
!@var ke turbulent diffusivity for e equation
!@var gc countergradient term in heat flux (level 3)
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var n number of GCM layers
!@var tau B1*lscale/sqrt(2*e)
!@var as2 shear squared, (dudz)**2+(dvdz)**2
!@var an2 Brunt-Vaisala frequency, grav/T*dTdz
!@var se stability constant for e, adjustable
!@var st2 stability constant for t2, adjustable
      USE CONSTANT, only : grav
      USE SOCPBL, only : ghmin,ghmax,gmmax0,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,b1
     *     ,g0,d1_3,d2_3,d3_3,d4_3,d5_3
     *     ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *     ,g2,g3,g4,g5,nlevel

      implicit none

      integer, intent(in) :: n,iter    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,e,t2,lscale
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2
      real*8, dimension(n), intent(in) :: dz,dze
      real*8, dimension(n), intent(out) :: km,kh,kq,ke,kt2,gc,gm,gh
     2                                    ,gc_t2

      ! note e *tau = b1/2 *lscale * qturb
      real*8, parameter ::  se=0.1d0,st2=0.06d0
      real*8 :: ell,den,qturb,tau,ghj,gmj,gmmax
      real*8 :: sm,sh,taue,e_lpbl,t21,tmp
      real*8 kh_canuto,c8,sig,sw,tpj,tpjm1,tppj,w3pj,taupj,m
      real*8, dimension(n) :: taua,delta,w2,w3
      integer :: j  !@var j loop variable

      do j=2,n
        ell=lscale(j)
        qturb=sqrt(2.d0*e(j))
        tau=B1*ell/qturb
        ghj=tau*tau*an2(j)
        gmj=tau*tau*as2(j)
        if(ghj.lt.ghmin) ghj=ghmin
        if(ghj.gt.ghmax) ghj=ghmax
        gmmax=(1+d1*ghj+d3*ghj*ghj)/(d2+d4*ghj)
        gmmax=min(gmmax,gmmax0)
        if(gmj.gt.gmmax) gmj=gmmax
        if(nlevel.eq.25) then
          den=1+d1*ghj+d2*gmj+d3*ghj*ghj+d4*ghj*gmj+d5*gmj*gmj
          sm=(s0+s1*ghj+s2*gmj)/den
          sh=(s4+s5*ghj+s6*gmj)/den
          gc(j)=0.d0
        elseif(nlevel.eq.3) then
          den=1+d1_3*ghj+d2_3*gmj+d3_3*ghj*ghj+d4_3*ghj*gmj
     &         +d5_3*gmj*gmj
          t21=(g_alpha(j)*tau)**2*t2(j)/e(j)
          sm=(s0_3+s1_3*ghj+s2_3*gmj+s3_3*t21)/den
          sh=(s4_3+s5_3*ghj+s6_3*gmj)/den
          tmp=(1+g4/g5*ghj+(g3**2-1./3.*g2**2)*gmj)/(g5*den)
          gc_t2(j)=tmp*g0*g_alpha(j)*tau
          gc(j)=tmp*g0*g_alpha(j)*tau*t2(j)
        endif
        taue=tau*e(j)
        km(j)=max(taue*sm,1.d-30)
        kh(j)=max(taue*sh,1.d-30)
        kq(j)=kh(j)
        ke(j)=max(taue*se,1.d-30)
        kt2(j)=max(taue*st2,1.d-30)
        gm(j)=gmj
        gh(j)=ghj

c!! Canuto modifications (6-01):
c        w2(j)=2./3.*e(j)+tau/3.*(-(3*g3-g2)*km(j)*as2(j)
c     2        +4*g4*g_alpha(j)*(-kh(j)*dtdz(j)+gc(j)))
c        taua(j)=tau
c        delta(j)=-2./3.*g_alpha(j)*dtdz(j)*0.82/11.04*tau*tau
c!! end of Canuto modifications (6-01):
       end do
       
c!! Canuto modifications (6-01):
c      sig=0.2d0
c      sw=(1-2*sig)/sqrt(sig*(1-sig))
c      do j=2,n
c        w3(j)=sw*abs(w2(j))**1.5d0
c      end do
c      c8=8.
c      do j=3,n-1
c        tpj=(t(j+1)-t(j-1))/(dz(j-1)+dz(j))
c        tpjm1=(t(j)-t(j-2))/(dz(j-2)+dz(j-1))
c        tppj=(tpj-tpjm1)/dz(j-1)
c        w3pj=(w3(j+1)-w3(j-1))/(dze(j-1)+dze(j))
c        taupj=(taua(j+1)-taua(j-1))/(dze(j-1)+dze(j))
c        m=1./(2*c8*(1.-delta(j)))*tau*tau/(11.04*kh(j))
c     2    *(tppj/dtdz(j)*w3(j)+w3pj+w3(j)/taua(j)*taupj)
c        kh_canuto=taua(j)/11.04*w2(j)/(1-delta(j))
c        write(68,1001) iter,j,m,kh(j),kh_canuto
c      end do
c        write(68,*)
c 1001 format(i4,1x,i4,1x,9(1pe14.4))
c!! end of Canuto modifications (6-01):

      ke(1)=b1*lscale(1)*sqrt(0.5d0*e(1))*se
      kt2(1)=ke(1)*st2/se
      do j=1,n-1
        ke(j)=0.5d0*(ke(j)+ke(j+1)) !defined on main levels
        kt2(j)=0.5d0*(kt2(j)+kt2(j+1)) !defined on main levels
      end do

      return
      end subroutine kgcm

      subroutine find_pbl_top(e,dbll,n)
!@sum find_pbl_top Find the pbl top (at main level lpbl)
!@+   if e(j) le certain fraction of e(1), real(j) is pbl top
!@auth  Ye Cheng
!@ver   1.0
!@var dbll the (real*8) layer number corresponds to the top of the pbl

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: e
      real*8, intent(out) :: dbll

      real*8, parameter :: fraction = 0.1d0
      real*8 :: e1p    ! certain percent of e_1
      integer j

      e1p=fraction*e(1)
      do j=2,n
        if (e(j).lt.e1p) exit
      end do
      dbll=j   ! dbll is real*8
      return
      end subroutine find_pbl_top

      subroutine ave_uv_to_tcell(u,v,u_tcell,v_tcell,im,jm,lm)
!@sum ave_uv_to_tcell Computes u_tcell,v_tcell from u and v,
!@+   the x and y components of a vector defined on the secondary grid
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
!@sum ave_ufvf_to_ucell Computes uf_ucell and vf_ucell from uf and vf,
!@+   the x and y components of a vector defined on the primary grid
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
      USE CONSTANT, only : by3
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
            s_ucell(l,i,2)=(s(l,1,1)+s(l,i,2)+s(l,ip1,2))*by3
            s_ucell(l,i,jm)=(s(l,1,jm)+s(l,i,jm-1)+s(l,ip1,jm-1))*by3
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
      use constant, only : by3
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
            s_ucell(l,i,2)=(s(l,1,1)+s(l,i,2)+s(l,ip1,2))*by3
            s_ucell(l,i,jm)=(s(l,1,jm)+s(l,i,jm-1)+s(l,ip1,jm-1))*by3
          end do
        end do

      return
      end subroutine ave_to_ucell1
