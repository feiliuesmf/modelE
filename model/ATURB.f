#include "rundeck_opts.h"

      Subroutine atm_diffus(lbase_min,lbase_max,dtime)
!@sum  atm_diffus updates u,v,t,q due to
!@+  turbulent transport throughout all GCM layers
!@+  using a second order closure (SOC)
!@+  turbulence model developed at GISS, 2000.
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from diffB347D6M20)
!@cont atm_diffus,getdz,dout,diff_uv,diff_t,diff_q,diff_,diff_t2
!@cont lgcm,kgcm,find_pbl_top,find_ew,ave_uv_to_agrid,ave_s_to_bgrid
!@var lbase_min/max levels through which to apply turbulence (dummy)
!@var dtime time step
!@var qmin minimum value of specific humidity 
!@var itest longitude at which to call dout
!@var jtest latitude at which to call dout
!@var call_diag logical variable whether dout is called

      USE DYNAMICS, only : pk,pdsig,plij,pek
      USE MODEL_COM, only :
     *      im,jm,lm,sig,sige,u,v,t,q,p,itime
      USE CONSTANT, only : grav,deltx,lhe,sha
      USE PBLCOM, only : tsavg,qsavg,dclev,uflux,vflux,tflux,qflux,egcm
     *     ,t2gcm,uflux1,vflux1,tflux1,qflux1
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj
      USE DAGCOM, only : ajl,
     &     jl_trbhr,jl_damdc,jl_trbke,jl_trbdlht
      USE SOCPBL, only : b2,prt,kappa,zgs,nlevel
cc      USE QUSDEF, only : nmom,zmoms,xymoms
cc      USE SOMTQ_COM, only : tmom,qmom
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,trm  !,trmom
      USE FLUXES, only : trflux1=>tot_trsource
      USE TRACER_DIAG_COM, only: tajln,jlnt_turb
#endif

      IMPLICIT NONE

      integer, intent(in) :: lbase_min,lbase_max
      real*8, intent(in) :: dtime

      real*8, parameter :: qmin=1.d-20
      integer, parameter :: itest= 1,jtest=11
      logical, parameter :: call_diag=.false.
      integer, SAVE :: ifirst=0

      real*8, dimension(lm) :: uij,vij,tij,qij,eij,t2ij
      real*8, dimension(lm) :: u0ij,v0ij,t0ij,q0ij,e0ij,t20ij
      real*8, dimension(lm) :: dudz,dvdz,dtdz,dqdz,g_alpha,as2,an2
      real*8, dimension(lm) :: rhoebydz,bydzerho
      real*8, dimension(lm) :: km,kh,kq,ke,kt2,gc,gc_t2,ew_rest,lscale
      real*8, dimension(lm) :: rhoij,rhoeij,dzij,dzeij,gm,gh
      real*8, dimension(im,jm,lm) :: uold
      real*8, dimension(lm,im,jm) :: rho,rhoe,dz,dzedge
      real*8, dimension(lm,im,jm) :: u_agrid,v_agrid,t_virtual
      real*8, dimension(lm,im,jm) :: km_gcm,km_gcm_bgrid
     2        ,dz_bgrid,dzedge_bgrid,rho_bgrid,rhoe_bgrid
      real*8, dimension(im,jm) :: tvsurf,uflux_bgrid,vflux_bgrid
cc      real*8, dimension(nmom,lm) :: tmomij,qmomij

      real*8 :: uflx,vflx,tflx,qflx,tvs
      real*8 :: ustar2,dbll,reserv,t0ijl,rak,alpha1,thes,ustar,qs
      integer :: imax,kmax,idik,idjk
      integer :: i,j,l,k,n,iter !@i,j,l,k loop variable
#ifdef TRACERS_ON
!@var tr0ij initial vertical tracer concentration profile (kg/kg)
!@var trij vertical tracer concentration profile (kg/kg)
!@var trmomij vertical tracer concentration moment profile (kg/kg)
      real*8, dimension(lm,ntm) :: tr0ij,trij
cc      real*8, dimension(nmom,lm,ntm) :: trmomij
!@var trflx surface tracer flux (-w tr) (kg/kg m/s)
      real*8, dimension(ntm) :: trflx
#endif
cccc
c     real*8, save :: p1000k
c     REAL*8,save, DIMENSION(IM,JM,4) :: sum=0.
c     integer, save :: n_save=0
cccccc
      ! Note that lbase_min/max are here for backwards compatibility with
      ! original drycnv. They are only used to determine where the
      ! routine has been called from.

      if (lbase_min.eq.2) return       ! quit if called from main

      !  convert input T to virtual T
      do j=1,jm
        do i=1,imaxj(j)
            !@var tvsurf(i,j) surface virtual temperature 
            !@var tsavg(i,j) COMPOSITE SURFACE AIR TEMPERATURE (K)
          tvsurf(i,j)=tsavg(i,j)*(1.d0+deltx*qsavg(i,j))
          do l=1,lm
            ! t_virtual is virtual potential temp. referenced at 1 mb
            t_virtual(l,i,j)=t(i,j,l)*(1.d0+deltx*Q(i,j,l))
          end do
        end do
      end do

      ! integrate T,Q equations at agrids

      ! get u_agrid and v_agrid
      call ave_uv_to_agrid(u,v,u_agrid,v_agrid,im,jm,lm)

      call getdz(t_virtual,p,dz,dzedge,rho,rhoe,tvsurf,im,jm,lm)
 
c     if(ifirst.eq.0) then
c       ifirst=1
c       do j=1,jm
c         do i=1,imaxj(j)
c           do l=1,lm
c              t2gcm(l,i,j)=egcm(l,i,j)*1.d-3
c            write(97,*) j,i,l,egcm(l,i,j)
c           end do
c         end do
c       end do
c     endif

      loop_j_tq: do j=1,jm
        loop_i_tq: do i=1,imaxj(j)

          do l=1,lm
            uij(l)=u_agrid(l,i,j)
            vij(l)=v_agrid(l,i,j)
            ! virtual potential temp. referenced at 1 mb
            tij(l)=t_virtual(l,i,j)
            if(q(i,j,l).lt.0.d0) q(i,j,l)=qmin
            qij(l)=q(i,j,l)
cc            qmomij(:,l)=qmom(:,i,j,l)
cc            tmomij(:,l)=tmom(:,i,j,l) ! vertical gradients should virtual ?
            eij(l)=egcm(l,i,j)
c           t2ij(l)=t2gcm(l,i,j)
            rhoij(l)=rho(l,i,j)
            rhoeij(l)=rhoe(l,i,j)
            t0ij(l)=tij(l)
            q0ij(l)=qij(l)
            e0ij(l)=eij(l)
c           t20ij(l)=t2ij(l)
#ifdef TRACERS_ON
            do n=1,ntm
              if (itime_tr0(n).le.itime) then
                trij(l,n)=trm(i,j,l,n)
cc                trmomij(:,l,n)=trmom(:,i,j,l,n)
                tr0ij(l,n)=trij(l,n)
              end if
            end do
#endif
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

          ! tvs is surface virtual potential temp. referenced at 1 mb
          tvs=tvsurf(i,j)/pek(1,i,j)
          uflx=uflux1(i,j)/rhoe(1,i,j)
          vflx=vflux1(i,j)/rhoe(1,i,j)
          tflx=tflux1(i,j)/(rhoe(1,i,j)*pek(1,i,j)) !referenced at 1 mb
          qflx=qflux1(i,j)/rhoe(1,i,j)
          ! redefine uflux1,vflux1 for later use
          uflux1(i,j)=uflx
          vflux1(i,j)=vflx
#ifdef TRACERS_ON
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              trflx(n)=trflux1(n,i,j)/rhoe(1,i,j)
            end if
          end do
#endif

c         if ((i.eq.itest).and.(j.eq.jtest)) then
c           write(99,*) i,j
c           write(99,1001) "uflux: ", uflux(i,j),uflx,uflux(i,j)/uflx
c           write(99,1001) "vflux: ", vflux(i,j),vflx,vflux(i,j)/vflx
c           write(99,1001) "tflux: ", tflux(i,j)/pek(1,i,j),tflx
c    &                            , (tflux(i,j)/pek(1,i,j))/tflx
c           write(99,1001) "qflux: ", qflux(i,j),qflx,qflux(i,j)/qflx
c           write(99,*) 
c         endif

          ustar=(uflx*uflx+vflx*vflx)**(0.25d0)
          ustar2=ustar*ustar
          alpha1=atan2(vflx,uflx)

          ! calculate z-derivatives at the surface 
          ! @var zgs height of surface layer (m), imported from SOCPBL
          thes=tflx*prt/ustar
          qs=qflx*prt/ustar
          dudz(1)=ustar/(kappa*zgs)*cos(alpha1)
          dvdz(1)=ustar/(kappa*zgs)*sin(alpha1)
          dtdz(1)=thes/(kappa*zgs)
          dqdz(1)=qs/(kappa*zgs)
          g_alpha(1)=grav/tvs
          !@var an2 brunt-vassala frequency
          !@var as2 shear number squared
          an2(1)=g_alpha(1)*dtdz(1)
          as2(1)=dudz(1)*dudz(1)+dvdz(1)*dvdz(1)

          ! calculate turbulence length scale lscale
          call lgcm(lscale,eij,as2,an2,dzij,dzeij,rhoeij,lm)

          ! calculate turbulent diffusivities km,kh,kq,ke and kt2
          call kgcm(km,kh,kq,ke,kt2,gc,gc_t2,ew_rest,gm,gh,uij,vij,tij,
     2             eij,t20ij,dudz,dvdz,as2,dtdz,g_alpha,an2,
     3             lscale,dzij,dzeij,tvs,lm)

          ! integrate eqn for T.K.E, e.
          call diff_e(e0ij,eij,km,kh,ke,gc,ew_rest,lscale,uij,vij,tij,
     2         dzij,dzeij,dudz,dvdz,as2,dtdz,g_alpha,an2,
     3         rhoij,rhoeij,ustar2,dtime,lm)

          ! integrate eqn for turb. virt. pot. temp. variance t2
c         if(nlevel.eq.3) call diff_t2(t20ij,t2ij,eij,kh,kt2,gc,gc_t2,
c    2      lscale,tij,dtdz,dzij,dzeij,rhoij,rhoeij,
c    3      ustar2,tflx,dtime,lm)

          ! integrate eqn for virtual potential temperature tij
          call diff_t(t0ij,tij,kh,gc,dzij,dzeij,
     2                 rhoij,rhoeij,rhoebydz,bydzerho,tflx,dtime,lm)

C**** also diffuse moments
cc        call diff_mom(tmomij)

          ! integrate eqn for specific humidity qij
          call diff_q(q0ij,qij,kq,dzij,dzeij,
     2         rhoij,rhoeij,rhoebydz,bydzerho,qflx,dtime,lm)

C**** also diffuse moments
cc        call diff_mom(qmomij)

#ifdef TRACERS_ON
C**** Use q diffusion coefficient for tracers
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              call diff_q(tr0ij(1,n),trij(1,n),kq,dzij,dzeij,
     2             rhoij,rhoeij,rhoebydz,bydzerho,trflx(1,n),dtime,lm)
cc        call diff_mom(trmomij)
            end if
          end do
#endif

          !@var dbll PBL top layer number counted from below, real*8
          call find_pbl_top(eij,dbll,lm)
          dclev(i,j)=dbll

          do l=1,lm
            ! update 3-d q,t,egcm,t2gcm and km_gcm
            q(i,j,l)=max(qij(l),qmin)
            t0ijl=t(i,j,l)
            t(i,j,l)=tij(l)/(1.d0+deltx*Q(i,j,l))
C**** moment variation to be added
cc            qmom(:,i,j,l)=qmomij(:,l)
cc            tmom(:,i,j,l)=tmomij(:,l)
            egcm(l,i,j)=eij(l)
c           t2gcm(l,i,j)=t2ij(l)
            km_gcm(l,i,j)=km(l)
            ! ACCUMULATE DIAGNOSTICS for t and q
            AJL(J,L,JL_TRBHR)=AJL(J,L,JL_TRBHR)
     2                 +(T(I,J,L)-t0ijl)*PK(L,I,J)*PLIJ(L,I,J)
            AJL(J,L,JL_TRBDLHT)=AJL(J,L,JL_TRBDLHT)
     2                 +(Q(I,J,L)-q0ij(l))*PDSIG(L,I,J)*LHE/SHA
            AJL(J,L,JL_TRBKE)=AJL(J,L,JL_TRBKE)+eij(l)
#ifdef TRACERS_ON
            do n=1,ntm
              if (itime_tr0(n).le.itime) then
                tajln(j,l,jlnt_turb,n)=tajln(j,l,jlnt_turb,n) +
     &               (trij(l,n)-tr0ij(l,n))
                trm(i,j,l,n)=trij(l,n)
cc                trmom(:,i,j,l,n)=trmomij(:,l,n)
              end if
            end do
#endif
          end do

          ! Write out diagnostics if at selected grid point:

          if (call_diag.and.(i.eq.itest).and.(j.eq.jtest)) then
            call dout(uij,vij,tij,qij,eij,ew_rest,
     1           dzij,dzeij,dudz,dvdz,as2,dtdz,g_alpha,an2,dqdz,
     2           rhoij,rhoeij,
     3           km,kh,kq,ke,kt2,gc,gm,gh,lscale,reserv,tvs,
     4           uflx,vflx,tflx,qflx,i,j,iter,lm)
          endif

        end do loop_i_tq
      end do loop_j_tq

      ! integrate U,V equations at bgrids

cccc testing:
cccc testing:
c     do j=1,jm
c       do i=1,imaxj(j)
c         uflux(i,j)=1.d0
c         vflux(i,j)=1.d0
c         do l=1,lm
c           rho(l,i,j)=1.d0
c         end do
c       end do
c     end do
c     call ave_uv_to_bgrid(uflux,vflux,uflux_bgrid,vflux_bgrid,im,jm,1)
c     call ave_S_to_bgrid(rho,rho_bgrid,im,jm,lm)
c     do j=2,jm
c       do i=1,im
c         if(abs(uflux_bgrid(i,j)-1.d0).gt.1.d-4.or.
c    &       abs(vflux_bgrid(i,j)-1.d0).gt.1.d-4) then
c            write(70,1003) i,j,uflux_bgrid(i,j),vflux_bgrid(i,j)
c         endif
c         if(abs(rho_bgrid(1,i,j)-1.d0).gt.1.d-4.or.
c    &       abs(rho_bgrid(2,i,j)-1.d0).gt.1.d-4.or.
c    &       abs(rho_bgrid(3,i,j)-1.d0).gt.1.d-4.or.
c    &       abs(rho_bgrid(lm-1,i,j)-1.d0).gt.1.d-4.or.
c    &       abs(rho_bgrid(lm,i,j)-1.d0).gt.1.d-4) then
c           write(71,1003) i,j,rho_bgrid(1,i,j),rho_bgrid(2,i,j),
c    &       rho_bgrid(3,i,j), rho_bgrid(lm-1,i,j),rho_bgrid(lm,i,j)
c         endif
c       end do
c     end do
c1003 format(2x,i4,2x,i4,2x,9(1pe14.4))
c
c     stop
cccc end testing:
cccc end testing:

      call ave_uv_to_bgrid(uflux1,vflux1,uflux_bgrid,vflux_bgrid,
     &                     im,jm,1)
      call ave_s_to_bgrid(km_gcm,km_gcm_bgrid,im,jm,lm)
      call ave_s_to_bgrid(dz,dz_bgrid,im,jm,lm)
      call ave_s_to_bgrid(dzedge,dzedge_bgrid,im,jm,lm)
      call ave_s_to_bgrid(rho,rho_bgrid,im,jm,lm)
      call ave_s_to_bgrid(rhoe,rhoe_bgrid,im,jm,lm)

      loop_j_uv: do j=2,jm
        loop_i_uv: do i=1,im

          do l=1,lm
            uij(l)=u(i,j,l)
            vij(l)=v(i,j,l)
            rhoij(l)=rho_bgrid(l,i,j)
            rhoeij(l)=rhoe_bgrid(l,i,j)
            u0ij(l)=uij(l)
            v0ij(l)=vij(l)
            uold(i,j,l)=uij(l)
            dzeij(l)=dzedge_bgrid(l,i,j)
            bydzerho(l)=1.d0/(dzeij(l)*rhoij(l))
            km(l)=km_gcm_bgrid(l,i,j)
          end do
          do l=1,lm-1
            dzij(l)=dz_bgrid(l,i,j)
            rhoebydz(l+1)=rhoeij(l+1)/dzij(l)
          end do

          uflx  =uflux_bgrid(i,j)
          vflx  =vflux_bgrid(i,j)

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

 1001 format(a,3(1pe14.4))

cccccc testing 9-6-01
c     p1000k=1000.**0.2862d0
c     n_save=n_save+1
c     rewind(71)
c     do i=1,im
c     do j=1,jm
c       do l=1,4
c         sum(i,j,l)=sum(i,j,l)+t(i,j,l)*p1000k*
c    &               (1.+0.6078d0*q(i,j,l))
c       end do
c       write(71,1011) n_save,i,j,
c    &                 sum(i,j,1)/n_save,sum(i,j,2)/n_save,
c    &                 sum(i,j,3)/n_save,sum(i,j,4)/n_save
c     end do
c     end do
c1011 format(3(2x,i4),2x,,9(1pe14.6))
c
cccccc end testing

      return
      end subroutine atm_diffus

      subroutine getdz(tv,p,dz,dzedge,rho,rhoe,tvsurf,im,jm,lm)
!@sum  getdz computes the 3d finite difference dz and dzedge
!@+    as well as the 3d density rho and rhoe
!@+    called at the primary grid (A-grid)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var  dz(l,i,j) z(l+1,i,j) - z(l,i,j)
!@var  dzedge(l,i,j) zedge(l+1,i,j) - zedge(l,i,j)
!@var  z vertical coordinate associated with SIG(l)
!@var  zedge vertical coordinate associated with SIGE(l)
!@var  temp0 virtual temperature (K) at (i,j) and SIG(l)
!@var  temp1 virtual temperature (K) at (i,j) and SIG(l+1)
!@var  temp1e average of temp0 and temp1

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
      USE CONSTANT, only : grav,rgas
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pmid,pk,pedn
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) :: tv
      real*8, dimension(im,jm), intent(in) :: p,tvsurf
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
              rhoe(1,i,j)=100.d0*ple/(tvsurf(i,j)*rgas)
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

      subroutine dout(u,v,t,q,e,ew_rest,dz,dze,
     1                dudz,dvdz,as2,dtdz,g_alpha,an2,dqdz,rho,rhoe,
     2                km,kh,kq,ke,kt2,gc,gm,gh,lscale,reserv,tvs,
     3                uflx,vflx,tflx,qflx,i,j,iter,n)
!@sum dout writes out diagnostics at i,j
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var u z-profile of west-east   velocity component
!@var v z-profile of south-north velocity component
!@var t z-profile of virt. pot. temperature (referenced to 1mb)
!@var p z-profile of pressure at main grid z
!@var pe z-profile of pressure at secondary grid zedge
!@var q z-profile of specific humidity
!@var e z-profile of turbulent kinetic energy
!@var t2 z-profile of turbulent virt. pot. temperature variance
!@var dz(j) z(j+1)-z(j)
!@var dze(j) zedge(j+1)-zedge(j)
!@var dudz z-derivative of u at secondary grid zedge
!@var dvdz z-derivative of v at secondary grid zedge
!@var as2 dudz^2+dvdz^2
!@var dtdz z-derivative of t at secondary grid zedge
!@var g_alpha grav*(thermal expansion coefficient)
!@var an2 g_alpha*dtdz, brunt-vasala frequency
!@var dqdz z-derivative of q at secondary grid zedge
!@var rho z-profile of density at z
!@var rhoe z-profile of density at zedge
!@var km turbulent viscosity for u and v equations
!@var kh turbulent diffusivity for t
!@var kq turbulent diffusivity for q
!@var ke turbulent diffusivity for e
!@var kt2 turbulent diffusivity for t2
!@var gc countergradient term is vertical heat flux (at model level 3)
!@var gm normalized velocity gradient, tau**2*as2
!@var gh normalized temperature gradient, tau**2*an2
!@var lscale z-profile of turbulent length scale
!@var tvs surface virtual temperature
!@var uflx momentun flux -uw at surface, zedge(1)
!@var vflx momentun flux -vw at surface, zedge(1)
!@var tflx heat flux -wt at surface, zedge(1)
!@var qflx moisture flux -wq at surface, zedge(1)
!@var i/j location at which the output is wriiten
!@var n number of vertical main layers

      USE CONSTANT, only : grav,rgas
      USE SOCPBL, only : b1,b2
      USE MODEL_COM, only : sig,sige
      USE RESOLUTION, only : PTOP
      USE DYNAMICS, only : pmid,pedn,plij,pk,pek

      implicit none

      integer, intent(in) :: n,iter,i,j
      real*8, dimension(n), intent(in) :: u,v,t,q,e,ew_rest
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2,dqdz
      real*8, dimension(n), intent(in) :: rho,rhoe
      real*8, dimension(n), intent(in) :: km,kh,kq,ke,kt2,gc,gm,gh
      real*8, dimension(n), intent(in) :: lscale,dz,dze
      real*8, intent(in) :: reserv,tvs
      real*8, intent(in) :: uflx,vflx,tflx,qflx

      real*8, dimension(n) :: p,pe,t2
      real*8 :: z,utotal,zedge,dmdz
      real*8 :: uf,hf,qf
      real*8 :: ri,rif,sigmat,reserv2,wt1
      integer :: l  !@var l loop variable

      do l=1,n
          p(l)=100.d0*pmid(l,i,j)
          pe(l)=100.d0*pedn(l,i,j)
          t2(l)=b2*lscale(l)/sqrt(2.d0*e(l))*kh(l)*dtdz(l)**2
      end do
      Write (67,1000) "iter=",iter, "i=",i,"j=",j
      Write (67,1100) "pe(1)=",pe(1)
      Write (67,1100) "uflx=",uflx,"vflx=",vflx
      Write (67,1100) "tflx=",tflx*pek(1,i,j),"qflx=",qflx

      ! Fields on main vertical grid:
      z=(rgas/grav)*0.5d0*(tvs*pek(1,i,j)+t(1)*pk(1,i,j))
     2  *log(pe(1)/p(1))+10.d0
      write (67,1500)
      do l=1,n-1
        write (67,2000) l,z,p(l),dz(l),u(l),v(l),t(l)*pk(l,i,j),
     2                  q(l),ke(l),rho(l),rhoe(l),ew_rest(l)
        z=z+dz(l)
      end do
      utotal=sqrt(u(n)*u(n)+v(n)*v(n))
      write (67,2500) l,z,p(n),u(n),v(n),t(n)*pk(n,i,j),q(n),0.,
     2                  rho(n),rhoe(n),0.
      write (67,*)
 
      do l=1,n
        write(68,1001) plij(l,i,j)*sig(l)+ptop,
     &       t(l)*pk(l,i,j),t(l)*pk(l,i,j)/(1.+0.61*q(l)),q(l)
      end do

      ! Fields on secondary vertical grid:
 
      write (67,3000)
      l=1
      zedge=10.d0
      write (67,2100) l,zedge,pe(l),lscale(l),e(l),t2(l)*pek(l,i,j)**2
      zedge=10.d0+dze(1)
      do l=2,n
        wt1=-kh(l)*dtdz(l)*pek(l,i,j)
        write (67,2000) l,zedge,pe(l),gc(l),wt1,kh(l),kt2(l),
     2       gm(l),gh(l),lscale(l),e(l),t2(l)*pek(l,i,j)**2
        zedge=zedge+dze(l)
      end do
      write (67,*)
 
      ! Fluxes on the secondary grid:
 
      zedge=10.d0
      write (67,4000)
      do l=2,n
        zedge=zedge+dze(l-1)
        dmdz=sqrt(as2(l))
        uf=km(l)*dmdz
        hf=kh(l)*dtdz(l)*pek(l,i,j)
        qf=kq(l)*dqdz(l)
        ri=an2(l)/as2(l)
        sigmat=km(l)/kh(l)
        rif=ri/sigmat
        reserv2=0.d0
        write (67,2000) l,zedge,uf,hf,qf,dmdz,dtdz(l)*pek(l,i,j),
     2                  dqdz(l),ri,rif,sigmat,reserv2
      end do
 
      write (67,*) "------------"
      write (67,*)

      return
  
1000  format(3(4x,a,i4))
1001  format(9(1pe14.4))
1100  format(2(4x,a,1pe14.4))
1500  format (1h ,' l',1x,
     1            '     z     ',1x,
     2            '     p     ',1x,'     dz    ',1x,'     u     ',1x,
     3            '     v     ',1x,'     t     ',1x,'     q     ',1x,
     4            '    ke     ',1x,'    rho   ',1x,'   rhoe   ',1x,
     5            '   ew_rest ',/)
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
c     alpha=0.d0
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
      Real*8 :: alpha,p4j,dgcdz
      integer :: j  !@var j loop variable

c     sub(j)*t_jm1_kp1+dia(j)*t_j_kp1+sup(j)*t_jp1_kp1 = rhs(j)

      do j=2,n-1
c         sub(j)=-dtime*kh(j)/(dz(j-1)*dzedge(j)*rho(j))*rhoe(j)
c         sup(j)=-dtime*kh(j+1)/(dz(j)*dzedge(j)*rho(j))*rhoe(j+1)
          sub(j)=-dtime*kh(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*kh(j+1)*rhoebydz(j+1)*bydzerho(j)
c         dgcdz=(gc(j+1)-gc(j))/dzedge(j)
c         p4j=-dgcdz
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=t0(j)
c         rhs(j)=t0(j)+dtime*p4j
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
c     alpha=0.d0
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=t0(1)
c     rhs(1)=t0(1)-dtime/(dzedge(1)*rho(1))*rhoe(1)*gc(2)
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
!@var q z-profle of specific humidity Q
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
c     alpha=0.d0
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

      subroutine diff_e(e0,e,km,kh,ke,gc,ew_rest,lscale,u,v,t,
     1           dz,dzedge,dudz,dvdz,as2,dtdz,g_alpha,an2,
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
      real*8, dimension(n), intent(in) :: e0,km,kh,ke,gc,ew_rest,lscale
      real*8, dimension(n), intent(in) :: u,v,t,dz,dzedge,rho,rhoe
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2
      real*8, dimension(n), intent(out) :: e
      real*8, intent(in) :: ustar2,dtime

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: qturb,tmp,d_ew_rest_dz
      integer :: j  !@var j loop variable

      real*8, parameter :: emin=1.d-20,emax=100.d0
c
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     j refers to the layer edge
c     except ke(j) which is defined on the layer middle

      do j=2,n-1
          d_ew_rest_dz=(ew_rest(j)-ew_rest(j-1))/dz(j-1)
          qturb=sqrt(2.d0*e0(j))
          tmp=rho(j-1)/rhoe(j)
          sub(j)=-dtime*ke(j-1)/(dz(j-1)*dzedge(j-1))*tmp
          tmp=rho(j)/rhoe(j)
          sup(j)=-dtime*ke(j)/(dz(j-1)*dzedge(j))*tmp
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*2*qturb/(b1*lscale(j))
          rhs(j)=e0(j)+dtime*(km(j)*as2(j)-kh(j)*an2(j)-d_ew_rest_dz)
c         rhs(j)=e0(j)+dtime*(km(j)*as2(j)-kh(j)*an2(j)
c    2           +g_alpha(j)*gc(j))
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

      subroutine lgcm(lscale,e,as2,an2,dz,dzedge,rhoe,n)
!@sum lgcm calculates the z-profle of the turbulent length scale
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var lscale z-profle of turbulent length scale
!@var e z-profle of pturbulent kinetic energy
!@var as2 z-profle of dudz^2+dvdz^2
!@var an2 z-profle of g_alpha*dtdz
!@var dz(j) z(j+1)-z(j)
!@var dzedge(j) zedge(j+1)-zedge(j)
!@var rho the z-profile of the density
!@var n number of GCM layers
!@var zgs height of surface layer (m), imported from SOCPBL

      USE CONSTANT, only : grav
      USE SOCPBL, only : kappa,zgs

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: e,rhoe
      real*8, dimension(n), intent(in) :: as2,an2
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
        sum1=sum1+dzedge(j-1)*(qj*rhoe(j)+qjm1*rhoe(j-1))
        sum2=sum2+dzedge(j-1)*(qj*rhoe(j)*zedge(j)+
     &                         qjm1*rhoe(j-1)*zedge(j-1))
      end do
      l0=alpha0*sum2/sum1

      kz=kappa*zedge(1)
      lscale(1)=l0*kz/(l0+kz)
       ! if (lscale(1).gt.dzedge(1)) lscale(1)=dzedge(1)
      do j=2,n
        kz=kappa*zedge(j)
        l1=l0*kz/(l0+kz)
        if (an2(j).gt.0.d0) then
          lmax  =0.53d0*sqrt(2.d0*e(j)/(an2(j)+1.d-40))
c         lmax2 =1.95d0*sqrt(2.d0*e(j)/(as2(j)+1.d-40))
c         lmax=min(lmax,lmax2)
          if (l1.gt.lmax) l1=lmax
        endif
        ! if (l1.gt.dzedge(j)) l1=dzedge(j)
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

      subroutine kgcm(km,kh,kq,ke,kt2,gc,gc_t2,ew_rest,
     2                gm,gh,u,v,t,e,t2,dudz,dvdz,as2,
     3                dtdz,g_alpha,an2,lscale,dz,dze,tvs,n)
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
!@var tvs surface virtual temperature
      USE CONSTANT, only : grav
      USE SOCPBL, only : ghmin,ghmax,gmmax0,d1,d2,d3,d4,d5
     *     ,s0,s1,s2,s4,s5,s6,b1,b2
     *     ,g0,d1_3,d2_3,d3_3,d4_3,d5_3
     *     ,s0_3,s1_3,s2_3,s3_3,s4_3,s5_3,s6_3
     *     ,g2,g3,g4,g5,g6,g7,nlevel

      implicit none

      integer, intent(in) :: n    !@var n  array dimension
      real*8, dimension(n), intent(in) :: u,v,t,e,lscale
      real*8, dimension(n), intent(in) :: dudz,dvdz,as2
      real*8, dimension(n), intent(in) :: dtdz,g_alpha,an2
      real*8, dimension(n), intent(in) :: dz,dze
      real*8, dimension(n), intent(inout) :: t2
      real*8, dimension(n), intent(out) :: km,kh,kq,ke,kt2,gc,gm,gh
     2                                    ,gc_t2,ew_rest
      real*8, intent(in) :: tvs

      ! note e *tau = b1/2 *lscale * qturb
      real*8, parameter ::  se=0.1d0,st2=0.06d0
      real*8 :: ell,den,qturb,tau,ghj,gmj,gmmax
      real*8 :: sm,sh,taue,e_lpbl,t21,tmp
      real*8 :: kh_canuto,c8,sig,sw,tpj,tpjm1,tppj,w3pj,taupj,m
      real*8, dimension(n) :: taua,uw,vw,w2,wt
      real*8, dimension(n) :: u2,v2,uv,ut,vt
      real*8 :: g_alphaj,tauj,dudzj,dvdzj,as2j,an2j,dtdzj
      real*8 :: uwj,vwj,w2j,wtj,du2dz,dv2dz,dw2dz 
      real*8 :: duvdz,duwdz,dvwdz,dutdz,dvtdz,dwtdz,dt2dz,ke0
      integer :: j  !@var j loop variable

      logical, parameter :: non_local=.false.

      do j=1,n
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
        km(j)=min(max(taue*sm,1.d-30),300.d0)
        kh(j)=min(max(taue*sh,1.d-30),300.d0)
        kq(j)=kh(j)
        ke(j)=min(max(taue*se,1.d-30),300.d0)
        kt2(j)=min(max(taue*st2,1.d-30),300.d0)
        taua(j)=tau
        gm(j)=gmj
        gh(j)=ghj
        ew_rest(j)=0.d0
      end do
      do j=1,n-1
        ke(j)=0.5d0*(ke(j)+ke(j+1)) !defined on main levels
        kt2(j)=0.5d0*(kt2(j)+kt2(j+1)) !defined on main levels
      end do

      if(non_local) then

        do j=1,n
          tau=taua(j)
          uw(j)=-km(j)*dudz(j)
          vw(j)=-km(j)*dvdz(j)
          wt(j)=-kh(j)*dtdz(j)
          u2(j)=2.d0/3*e(j)-tau/3*((g2+3*g3)*dudz(j)*uw(j)
     2        -2*g2*dvdz(j)*vw(j)+2*g4*g_alpha(j)*wt(j))
          v2(j)=2.d0/3*e(j)-tau/3*((g2+3*g3)*dvdz(j)*vw(j)
     2        -2*g2*dudz(j)*uw(j)+2*g4*g_alpha(j)*wt(j))
          w2(j)=2.d0/3*e(j)+tau/3*((3*g3-g2)*(
     2     dudz(j)*uw(j)+dvdz(j)*vw(j))+4*g4*g_alpha(j)*wt(j))
          uv(j)=-(g2+g3)/2*tau*(dvdz(j)*uw(j)+dudz(j)*vw(j))
          ut(j)=-tau/g5*(dtdz(j)*uw(j)+(g6+g7)/2*dudz(j)*wt(j))
          vt(j)=-tau/g5*(dtdz(j)*vw(j)+(g6+g7)/2*dvdz(j)*wt(j))
          t2(j)=b2*lscale(j)/sqrt(2.d0*e(j))*kh(j)*dtdz(j)**2
        end do
      
        do j=1,n-1  ! on main grids
          g_alphaj=0.5d0*(g_alpha(j)+g_alpha(j+1))
          tauj=0.5d0*(taua(j)+taua(j+1))
          dudzj=0.5d0*(dudz(j)+dudz(j+1))
          dvdzj=0.5d0*(dvdz(j)+dvdz(j+1))
          dtdzj=0.5d0*(dtdz(j)+dtdz(j+1))
          as2j=dudzj*dudzj+dvdzj*dvdzj
          an2j=g_alphaj*dtdzj
c         if(an2j.lt.0.d0) tauj=tauj/(1.d0-0.04d0*an2j*tauj*tauj)
          uwj=0.5d0*(uw(j)+uw(j+1))
          vwj=0.5d0*(vw(j)+vw(j+1))
          w2j=0.5d0*(w2(j)+w2(j+1))
          wtj=0.5d0*(wt(j)+wt(j+1))
          du2dz=(u2(j+1)-u2(j))/dze(j)
          dv2dz=(v2(j+1)-v2(j))/dze(j)
          dw2dz=(w2(j+1)-w2(j))/dze(j)
          duvdz=(uv(j+1)-uv(j))/dze(j)
          duwdz=(uw(j+1)-uw(j))/dze(j)
          dvwdz=(vw(j+1)-vw(j))/dze(j)
          dutdz=(ut(j+1)-ut(j))/dze(j)
          dvtdz=(vt(j+1)-vt(j))/dze(j)
          dwtdz=(wt(j+1)-wt(j))/dze(j)
          dt2dz=(t2(j+1)-t2(j))/dze(j)
          call find_ew(g_alphaj,tauj,dudzj,dvdzj,as2j,an2j,
     &      uwj,vwj,w2j,wtj,du2dz,dv2dz,dw2dz,
     &      duvdz,duwdz,dvwdz,dutdz,dvtdz,dwtdz,dt2dz,
     &      Ke(j),ew_rest(j))
          ke(j)=min(max(ke(j),1.d-30),300.d0)
        end do
      
      endif

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
      dbll=j-1   ! dbll is real*8
      return
      end subroutine find_pbl_top

      subroutine find_ew(ga,tau,dudz,dvdz,as2,an2,uw,vw,w2,wt,
     &   du2dz,dv2dz,dw2dz,duvdz,duwdz,dvwdz,dutdz,dvtdz,dwtdz,dt2dz,
     &   Ke,ew_rest)
c
c      q2w  = - K*dq2dz + q2w_rest    
c      ew  = - K*dedz + ew_rest, ew_rest = q2w_rest/2    
c output of 3m_eqns,3m_solve_sb0,3m_solve_sb0_more,more2,more3,more31,
c more32, more33 on kirk:/u/acyxc/papers/third/3m_publication
c tau=q2/epsilon=B1*ell/q, ell->0.4*z for small z (height)
c each t obsorbs a lamda=g*alpha*tau, ga==g*alpha       --- 04/6/00
c
      implicit none

      real*8, intent(in) :: ga,tau,dudz,dvdz,as2,an2,uw,vw,w2
      real*8, intent(in) :: du2dz,dv2dz,dw2dz,duvdz,duwdz,dvwdz 
      real*8, intent(inout) :: wt,dutdz,dvtdz,dwtdz,dt2dz 
      real*8, intent(out) :: Ke,ew_rest
     
      real*8, parameter :: c=1.d0/8.d0
      real*8, SAVE :: c0,c1,c2,c3,c4,c5,c6,c7,c8,c9
      real*8, SAVE :: c10,c11,c12,c13,c14,c15
      real*8, SAVE :: c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27
      real*8, SAVE :: d0,d1,d2,d3,d4,d5,d6,d7,d8
      real*8 :: U,V,S2,N2
      real*8 :: n0,n1,s1,b1,b2,b3,b4,b5,d
      real*8 :: m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,num_rest,q2w_rest
      integer, SAVE :: ifirst=0

      U=dudz*tau
      V=dvdz*tau
      S2=as2*tau*tau
      N2=an2*tau*tau
      wt    = ga*tau * wt
      dutdz = ga*tau * dutdz
      dvtdz = ga*tau * dvtdz
      dwtdz = ga*tau * dwtdz
      dt2dz = (ga*tau)**2 * dt2dz

      if (ifirst.eq.0) then
          ifirst=1
          c0 = -2*c**2*(13*c+9)
          c1 = -900
          c2 = -30*c*(44*c + 3)
          c3 = -54*c**3
          c4 = 25*c
          c5 = 5*c**3
          c6 = 25*c + 15
          c7 =  5*c**3
          c8 =  6*c**2*(3*c+2)
          c9   =  1500*(5*c+3)
          c10   = 750*c**2*(7*c+3)
          c11   = 50*c*(268*c**2+183*c+9)
          c12   = 25*c**3*(146*c**2+63*c+9)
          c13   = 30*c**3*(242*c**2+181*c+15)
          c14   = 30*c**5*(57*c**2+28*c-3)
          c15   = 36*c**5*(3*c+2)*(10*c+1)
          c16   = 4500*(5*c+3)
          c17   = 1500*c**2*(17*c+9)
          c18   = 50*c*(394*c**2+303*c+27)
          c19   = 50*c**3*(-17*c**2+24*c+27)
          c20   = 30*c**3*( 84*c**2+125*c+45 )
          c21   = 45*c**5*(7*c+3)
          c22   = 108*c**5*(3*c+2)
          c23   = 2050*c+1230
          c24   = 90*c**2*(29*c+16)
          c25   = 6*c**2*(321*c+209)
          c26   = 9*c**4*(57*c+31)
          c27   = 108*c**4*(3*c+2)
    
          d0 = 7500*(5*c+3)**2
          d1 = 3750*c**3*(23*c+18)
          d2 = 250*c*(5*c+3)*(292*c**2+237*c+9)
          d3 = 125*c**4*(526*c**2+645*c+54)
          d4 = 150*c**3*(1254*c**3+2081*c**2+933*c+72)
          d5 = 1350*c**6*(19*c**2+17*c+1)
          d6 = 90*c**5*(180*c**3+438*c**2+299*c+57)
          D7 = 405*c**8
          d8 = 324*c**7*(3*c+2)
      endif

      n0 = c1 + c2*N2 + c3*N2**2
      n1 = c4 + c5*N2
      s1 = -2*c5*U*V
      b1 = c6 +  c7*V**2 + c8*N2
      b2 = c6 +  c7*U**2 + c8*N2
      b3 = c9+c10*S2+(c11+c12*S2)*N2+(c13+c14*S2)*N2**2+c15*N2**3
      b4 = c16+c17*S2+(c18+c19*S2)*N2+(c20+c21*S2)*N2**2+c22*N2**3
      b5 = c23 + c24*S2 + (c25 + c26*S2)*N2 + c27*N2**2
      D = (d0 + d1*S2) + (d2 + d3*S2)*N2 + (d4 + d5*S2)*N2**2 + 
     &    (d6 + d7*S2)*N2**3 + d8*N2**4

      m1 = n0*(n1*c*(3*U*uw+V*vw)+b1*(5*w2+2*c*wt))
      m2 = n0*(n1*c*(3*V*vw+U*uw)+b2*(5*w2+2*c*wt))

      m3 = -n0*n1*(4*c+3)*( U*uw+V*vw )-45*b3*w2-6*b4*c*wt
      m4 = n0*(2*n1*c*(V*uw+U*vw)+s1*(5*w2+2*c*wt))

      m5 = n0*(10*b1*uw+5*s1*vw-2*n1*(4*c+3)*U*w2+5*c0*U*wt)
      m6 = n0*(10*b2*vw+5*s1*uw-2*n1*(4*c+3)*V*w2+5*c0*V*wt)

      m7 = n0*(2*c*(2*b1*uw+s1*vw)+c0*U*(5*w2+4*c*wt))
      m8 = n0*(2*c*(2*b2*vw+s1*uw)+c0*V*(5*w2+4*c*wt))

      m9 = 5*n0*c0*(U*uw+V*vw)+12*c*(-b4*w2-10*b5*c*wt)
      m10 = 2*c*(n0*c0*(U*uw+V*vw)-6*c*b5*(5*w2+wt))

c     Num = m1*du2dz+m2*dv2dz+m3*dw2dz+m4*duvdz+m5*duwdz+
c    &      m6*dvwdz+m7*dutdz+m8*dvtdz+m9*dwtdz+m10*dt2dz 
c     q2w = tau*c/(2*D)*Num
c         = tau*c/(2*D)*(m3*dq2dz + Num_rest)
c         = - Ke*dq2dz + q2w_rest    
      Num_rest = (m1-m3)*du2dz+(m2-m3)*dv2dz+m4*duvdz+m5*duwdz+
     &      m6*dvwdz+m7*dutdz+m8*dvtdz+m9*dwtdz+m10*dt2dz 
      Ke   = -tau*c/(2*D)*m3
      q2w_rest = tau*c/(2*D)*Num_rest
      ew_rest  = q2w_rest/2.d0
c     note: if you don't want q2w itself, but just want Ke and 
c           q2w_rest, then dw2dz is not actually used.

      return
      end subroutine find_ew

      subroutine ave_uv_to_agrid(u,v,u_a,v_a,im,jm,lm)
!@sum Computes u_a,v_a from u and v
!@var u x-component at secondary grids (B_grid)
!@var v y-component at secondary grids (B_grid)
!@var u_a x-component at primary grids (A_grid)
!@var v_a y-component at primary grids (A_grid)
!@auth Ye Cheng
!@ver  1.0
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(im,jm,lm), intent(in) ::  u,v
      real*8, dimension(lm,im,jm), intent(out) :: u_a,v_a

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
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
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
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
          END DO
        END DO
      END DO
C****
      return
      end subroutine ave_uv_to_agrid

      subroutine ave_uv_to_bgrid(u_a,v_a,u,v,im,jm,lm)
!@sum Computes u and v from u_a and v_a
!@var u_a x-component of wind at primary grids (A_grid)
!@var v_a y-component of wind at primary grids (A_grid)
!@var u x-component of wind at secondary grids (B_grid)
!@var v y-component of wind at secondary grids (B_grid)
!@auth Ye Cheng
!@ver  1.0

      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj,cosiv,siniv

      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) ::  u_a,v_a
      real*8, dimension(im,jm,lm), intent(out) :: u,v

      real*8 :: HEMI
      integer :: i,j,l,k,IMAX,KMAX

      u=0.d0; v=0.d0
c     polar boxes
      DO J=1,JM,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        HEMI=1.
        IF(J.LE.JM/2) HEMI=-1.
        DO I=1,IMAX
        DO L=1,LM
        DO K=1,KMAX
          U(IDIJ(K,I,J),IDJJ(K,J),L)=U(IDIJ(K,I,J),IDJJ(K,J),L)
     *      +RAVJ(K,J)*(U_A(L,I,J)*COSIV(K)+V_A(L,I,J)*SINIV(K)*HEMI)
          V(IDIJ(K,I,J),IDJJ(K,J),L)=V(IDIJ(K,I,J),IDJJ(K,J),L)
     *      +RAVJ(K,J)*(V_A(L,I,J)*COSIV(K)-U_A(L,I,J)*SINIV(K)*HEMI)
        END DO
        END DO
        END DO
      END DO
c     non polar boxes
      DO J=2,JM-1
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
        DO L=1,LM
        DO K=1,KMAX
          U(IDIJ(K,I,J),IDJJ(K,J),L)=U(IDIJ(K,I,J),IDJJ(K,J),L)
     *          +RAVJ(K,J)*U_A(L,I,J)
          V(IDIJ(K,I,J),IDJJ(K,J),L)=V(IDIJ(K,I,J),IDJJ(K,J),L)
     *          +RAVJ(K,J)*V_A(L,I,J)
        END DO
        END DO
        END DO
      END DO

      return
      end subroutine ave_uv_to_bgrid

      subroutine ave_s_to_bgrid(s,s_b,im,jm,lm)
!@sum Computes s_b from s
!@var s scalar at primary grid (A_grid)
!@var s_b scalar at secondary grid (B_grid)
!@auth Ye Cheng
!@ver  1.0

      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj

      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(lm,im,jm), intent(in) ::  s
      real*8, dimension(lm,im,jm), intent(out) :: s_b

      integer :: i,j,l,k,IMAX,KMAX

      S_B=0.d0
      DO J=1,JM
        IMAX=IMAXJ(J)
        KMAX=KMAXJ(J)
        DO I=1,IMAX
        DO L=1,LM
        DO K=1,KMAX
          S_B(L,IDIJ(K,I,J),IDJJ(K,J))=S_B(L,IDIJ(K,I,J),IDJJ(K,J))
     *          +RAVJ(K,J)*S(L,I,J)
        END DO
        END DO
        END DO
      END DO
1003  format(4(i4,1x),3(1pe14.4))
      return
      end subroutine ave_s_to_bgrid
