#include "rundeck_opts.h"

      subroutine atm_diffus(lbase_min,lbase_max,dtime)
!@sum  atm_diffus updates u,v,t,q due to turbulent transport throughout
!@+    all GCM layers using a non-local turbulence model
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@ver  1.0 (from diffB347D6M20)
!@cont atm_diffus,getdz,dout,de_solver_main,de_solver_edge,l_gcm,k_gcm,
!@+    e_gcm,find_pbl_top,zze,ave_uv_to_agrid,
!@+    ave_uv_to_bgrid,ave_s_to_bgrid,apply_fluxes_to_atm
!@var lbase_min/max levels through which to apply turbulence (dummy)
!@var dtime time step

!@var qmin minimum value of specific humidity
!@var itest/jtest longitude/latitude at which dout may be called
!@var call_diag logical variable whether dout is called

      USE CONSTANT, only : grav,deltx,lhe,sha,by3,teeny,mb2kg
      USE MODEL_COM, only :
     *     im,jm,lm,sig,sige,u_3d=>u,v_3d=>v,t_3d=>t,q_3d=>q,itime,psf
     *     ,pmtop
cc      USE QUSDEF, only : nmom,zmoms,xymoms
cc      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj,bydxyp,dxyp
      USE DYNAMICS, only : pk,pdsig,plij,pek,byam,am,dke
      USE DOMAIN_DECOMP, ONLY : grid, get, SOUTH, NORTH
      USE DOMAIN_DECOMP, ONLY : halo_update_column, checksum_column
      USE DOMAIN_DECOMP, ONLY : halo_update       , checksum
      USE DIAG_COM, only : ajl=>ajl_loc,jl_trbhr,jl_damdc
     *     ,jl_trbke,jl_trbdlht
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,trm,t_qlimit  !,trmom
      USE TRDIAG_COM, only: tajln=>tajln_loc,jlnt_turb
      USE FLUXES, only : trflux1
#endif
      USE SOCPBL, only : b1,b123,prt,kappa,zgs
      USE PBLCOM, only : tsavg,qsavg,dclev,uflux,vflux,tflux,qflux
     *     ,e_3d=>egcm,w2_3d=>w2gcm !,t2_3d=>t2gcm
      USE FLUXES, only : uflux1,vflux1,tflux1,qflux1


      IMPLICIT NONE

      integer, intent(in) :: lbase_min,lbase_max
      real*8, intent(in) :: dtime

      real*8, parameter :: qmin=1.d-20
      integer, parameter :: itest= 1,jtest=11
      logical, parameter :: call_diag=.false.

      real*8, dimension(lm) :: u,v,t,q,e,u0,v0,t0,q0,e0
     &    ,dudz,dvdz,dtdz,dqdz,g_alpha,as2,an2
     &    ,rhoebydz,bydzerho,rho,rhoe,dz,dze
     &    ,km,kh,ke,wt_nl,wq_nl
     &    ,lscale,qturb,p3,p4,rhobydze,bydzrhoe,uw,vw,wt,wq
      real*8, dimension(lm+1) :: ze

      real*8, dimension(lm,im,grid%j_strt_halo:grid%j_stop_halo) :: 
     &     u_3d_old,rho_3d,rhoe_3d,dz_3d
     &    ,dze_3d,u_3d_agrid,v_3d_agrid,t_3d_virtual,km_3d,km_3d_bgrid
     &    ,dz_3d_bgrid,dze_3d_bgrid,rho_3d_bgrid,rhoe_3d_bgrid
     &    ,wt_3d,v_3d_old
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: 
     &     tvsurf,uflux_bgrid,vflux_bgrid
cc      real*8, dimension(nmom,lm) :: tmomij,qmomij

      real*8 :: uflx,vflx,tvflx,qflx,tvs
     &   ,ustar2,t0ijl,tijl,rak,alpha1,ustar
     &   ,flux_bot,flux_top,x_surf
     &   ,wstar,dbl,lmonin,tpe0,tpe1,ediff
      integer :: idik,idjk,ldbl,kmax,
     &    i,j,l,k,n,iter !@i,j,l,k,n,iter loop variable
#ifdef TRACERS_ON
!@var tr0ij initial vertical tracer concentration profile (kg/kg)
!@var trij vertical tracer concentration profile (kg/kg)
!@var trmomij vertical tracer concentration moment profile (kg/kg)
!@var wc_nl non-local fluxes of tracers
      real*8, dimension(lm,ntm) :: tr0ij,trij,wc_nl
cc      real*8, dimension(nmom,lm,ntm) :: trmomij
!@var trflx surface tracer flux (-w tr) (kg/kg m/s)
      real*8, dimension(ntm) :: trflx
      integer nta,nx,ntix(ntm)
#endif

      INTEGER :: I_0, I_1, J_1, J_0, J_0H, J_1H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      ! Note that lbase_min/max are here for backwards compatibility
      ! with original drycnv. They are only used to determine where the
      ! routine has been called from.

      if (lbase_min.eq.2) return       ! quit if called from main

      !  convert input T to virtual T

!$OMP  PARALLEL DO PRIVATE (L,I,J)
      do j=J_0, J_1
        do i=1,imaxj(j)
          !@var tvsurf(i,j) surface virtual temperature
          !@var tsavg(i,j) composite surface air temperature (k)
          tvsurf(i,j)=tsavg(i,j)*(1.d0+deltx*qsavg(i,j))
          do l=1,lm
            ! t_3d_virtual is virtual potential temp. referenced at 1 mb
            t_3d_virtual(l,i,j)=t_3d(i,j,l)*(1.d0+deltx*q_3d(i,j,l))
          end do
        end do
      end do
!$OMP  END PARALLEL DO

#ifdef TRACERS_ON
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          nx=nx+1
          ntix(nx)=n
        end if
      end do
      nta=nx
#endif

      ! integrate equations other than u,v at agrids

      ! get u_3d_agrid and v_3d_agrid
      call ave_uv_to_agrid(u_3d,v_3d,u_3d_agrid,v_3d_agrid,lm)

      call getdz(t_3d_virtual,dz_3d,dze_3d,rho_3d,rhoe_3d,tvsurf
     &     ,lm)

!$OMP  PARALLEL DO PRIVATE (L,I,J,u,v,t,q,e,rho,rhoe,t0,q0,e0,qturb,
!$OMP*   dze,dz,bydzerho,rhobydze,bydzrhoe,rhoebydz,tvs,uflx,vflx,
!$OMP*   qflx,tvflx,ustar,ustar2,alpha1,dudz,dvdz,dtdz,dqdz,g_alpha,
!$OMP*   an2,as2,ze,lscale,dbl,ldbl,wstar,kh,km,ke,wt,wq,uw,vw,
!$OMP*   wt_nl,wq_nl,lmonin,p3,p4,x_surf,flux_bot,flux_top,t0ijl,tijl,
!$OMP*   tpe0,tpe1,ediff
#ifdef TRACERS_ON
!$OMP*   ,n,nx,trij,tr0ij,trflx,wc_nl
#endif
!$OMP*    ) SHARED(dtime
#ifdef TRACERS_ON
!$OMP*    ,nta
#endif
!$OMP*    ) SCHEDULE(DYNAMIC,2)

      loop_j_tq: do j=J_0, J_1
        loop_i_tq: do i=1,imaxj(j)

          do l=1,lm
            u(l)=u_3d_agrid(l,i,j)
            v(l)=v_3d_agrid(l,i,j)
            ! virtual potential temp. referenced at 1 mb
            t(l)=t_3d_virtual(l,i,j)
            q(l)=q_3d(i,j,l)
cc            qmomij(:,l)=qmom(:,i,j,l)
cc            tmomij(:,l)=tmom(:,i,j,l) ! vert. grad. should virtual ?
            if(q(l).lt.qmin) q(l)=qmin
            e(l)=e_3d(l,i,j) !e_3d was called egcM
            ! t2(l)=t2_3d(l,i,j)  ! not in use
            rho(l)=rho_3d(l,i,j)
            rhoe(l)=rhoe_3d(l,i,j)
            t0(l)=t(l)
            q0(l)=q(l)
            e0(l)=e(l)
            qturb(l)=sqrt(2.d0*e(l))
            dze(l)=dze_3d(l,i,j)
            dz(l)=dz_3d(l,i,j)
            bydzerho(l)=1.d0/(dze(l)*rho(l))
            rhobydze(l)=rho(l)/dze(l)
#ifdef TRACERS_ON
            do nx=1,nta
              n=ntix(nx)
              trij(l,nx)=trm(i,j,l,n)*byam(l,i,j)*bydxyp(j)
cc                trmomij(:,l,nx)=trmom(:,i,j,l,n)
              tr0ij(l,nx)=trij(l,nx)
            end do
#endif
          end do
          bydzrhoe(1)=0.
          rhoebydz(1)=0.
          do l=1,lm-1
            bydzrhoe(l+1)=1.d0/(dz(l)*rhoe(l+1))
            rhoebydz(l+1)=rhoe(l+1)/dz(l)
          end do

          ! tvs is surface virtual potential temp. referenced at 1 mb
          tvs=tvsurf(i,j)/pek(1,i,j)
          uflx=uflux1(i,j)/rhoe(1)
          vflx=vflux1(i,j)/rhoe(1)
          qflx=qflux1(i,j)/rhoe(1)
          ! tvflx is virtual, potential temp. flux referenced at 1 mb
          tvflx=tflux1(i,j)*(1.d0+deltx*qsavg(i,j))/(rhoe(1)*pek(1,i,j))
     &         +deltx*tsavg(i,j)/pek(1,i,j)*qflx
          ! redefine uflux1,vflux1 for later use
          uflux1(i,j)=uflx
          vflux1(i,j)=vflx
#ifdef TRACERS_ON
          do nx=1,nta
            n=ntix(nx)
C**** minus sign needed for ATURB conventions
            trflx(nx)=-trflux1(i,j,n)*bydxyp(j)/rhoe(1)
          end do
#endif

          if(abs(uflx).lt.teeny) uflx=sign(teeny,uflx)
          if(abs(vflx).lt.teeny) vflx=sign(teeny,vflx)
          ustar=(uflx*uflx+vflx*vflx)**(0.25d0)
          ustar2=ustar*ustar
          alpha1=atan2(vflx,uflx)

          ! calculate z-derivatives at the surface

          ! @var zgs height of surface layer (m), imported from SOCPBL
          dudz(1)=ustar/(kappa*zgs)*cos(alpha1)
          dvdz(1)=ustar/(kappa*zgs)*sin(alpha1)
          dtdz(1)=(tvflx*prt/ustar)/(kappa*zgs)
          dqdz(1)=(qflx*prt/ustar)/(kappa*zgs)

          g_alpha(1)=grav/tvs

          ! calculate z-derivatives on the edges of the layers

          do l=2,lm
            dudz(l)=(u(l)-u(l-1))/dz(l-1)
            dvdz(l)=(v(l)-v(l-1))/dz(l-1)
            dtdz(l)=(t(l)-t(l-1))/dz(l-1)
            dqdz(l)=(q(l)-q(l-1))/dz(l-1)
            g_alpha(l)=grav*2.d0/(t(l)+t(l-1))
          end do

          !@var an2 brunt-vassala frequency
          !@var as2 shear number squared

          do l=1,lm
            an2(l)=g_alpha(l)*dtdz(l)
            as2(l)=dudz(l)*dudz(l)+dvdz(l)*dvdz(l)
          end do

          ! calculate turbulence length scale lscale
          call zze(dze,ze,lm)
          call l_gcm(ze,lscale,lm)

          call find_pbl_top(e,ze,dbl,ldbl,lm)
          if(tvflx.lt.0.) then ! convective case
            wstar=(-g_alpha(1)*tvflx*dbl)**by3
          else
            wstar=0.
          endif
          ! calculate turbulent diffusivities km,kh and ke
          call k_gcm(tvflx,qflx,ustar,wstar,dbl
     &        ,ze,lscale,e,qturb,an2,as2,dtdz,dqdz,dudz,dvdz
     &        ,kh,km,ke,wt,wq,uw,vw,wt_nl,wq_nl
#ifdef TRACERS_ON
     &        ,trflx,wc_nl,nta
#endif
     &        ,lm)

          call e_gcm(tvflx,wstar,ustar,dbl,lmonin,ze,g_alpha
     &              ,an2,as2,lscale,e,lm)

          ! integrate differential eqn for e
          p3(1)=0. ; p3(lm)=0.
          p4(1)=0. ; p4(lm)=0.
c         do l=2,lm-1
c             p3(l)=2.d0*(qturb(l)/(b1*lscale(l)))
c             p4(l)=-uw(l)*dudz(l)-vw(l)*dvdz(l)+g_alpha(l)*wt(l)
c             ! p4(l)=km(l)*as2(l)-kh(l)*an2(l)
c         end do
c         x_surf=0.5d0*b123*ustar2
c         call de_solver_edge(e,e0,ke,p3,p4,
c    &        rhobydze,bydzrhoe,x_surf,dtime,lm,.true.)

          do l=1,lm
              qturb(l)=sqrt(2.d0*e(l))
          end do

          ! integrate differential eqn for T
          do l=2,lm-1
              p4(l)=-(rhoe(l+1)*wt_nl(l+1)-rhoe(l)*wt_nl(l))
     &              *bydzerho(l)
          end do
          flux_bot=rhoe(1)*tvflx+rhoe(2)*wt_nl(2)
          flux_top=0.
          call de_solver_main(t,t0,kh,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)

C**** also diffuse moments
cc        call diff_mom(tmomij)

          ! integrate differential eqn for Q
          do l=2,lm-1
              p4(l)=-(rhoe(l+1)*wq_nl(l+1)-rhoe(l)*wq_nl(l))
     &              *bydzerho(l)
C**** check on physicality of non-local fluxes....
              if ( p4(l)*dtime+q0(l).lt.0 ) then
                p4(l)=-q(l)/dtime
                wq_nl(l+1)=(q0(l)/(dtime*bydzerho(l))+rhoe(l)
     *               *wq_nl(l))/rhoe(l+1)
              end if
          end do
          flux_bot=rhoe(1)*qflx+rhoe(2)*wq_nl(2)
C**** fix first layer for rare tracer problems
C**** Does this ever happen for q? (put this in just in case)
            if ( q0(1)-dtime*bydzerho(1)*flux_bot.lt.0 ) then
              flux_bot=-q0(1)/(dtime*bydzerho(1))
              wq_nl(2)=(flux_bot-rhoe(1)*qflx)/rhoe(2)
            end if
          flux_top=0.

          call de_solver_main(q,q0,kh,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.true.)
          do l=1,lm
              if(q(l).lt.qmin) q(l)=qmin
          end do

C**** also diffuse moments
cc        call diff_mom(qmomij)

#ifdef TRACERS_ON
C**** Use q diffusion coefficient for tracers
C**** Note that non-local effects for tracers can be included
C**** parallel to the case of Q
          do n=1,nta
            do l=2,lm-1
              p4(l)=-(rhoe(l+1)*wc_nl(l+1,n)-rhoe(l)*wc_nl(l,n))
     &             *bydzerho(l)
C**** check on physicality of non-local fluxes....
              if ( p4(l)*dtime+tr0ij(l,n).lt.0 ) then
                p4(l)=-tr0ij(l,n)/dtime
                wc_nl(l+1,n)=(tr0ij(l,n)/(dtime*bydzerho(l))+rhoe(l)
     *               *wc_nl(l,n))/rhoe(l+1)
              end if
            end do
            flux_bot=rhoe(1)*trflx(n)+rhoe(2)*wc_nl(2,n) !tr0ij(1,n)
C**** fix first layer for rare tracer problems
            if ( tr0ij(1,n)-dtime*bydzerho(1)*flux_bot.lt.0 ) then
              flux_bot=-tr0ij(1,n)/(dtime*bydzerho(1))
              wc_nl(2,n)=(flux_bot-rhoe(1)*trflx(n))/rhoe(2)
            end if
            flux_top=0.

            call de_solver_main(trij(1,n),tr0ij(1,n),kh,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,t_qlimit(n))
cc          call diff_mom(trmomij)
          end do
#endif
          call find_pbl_top(e,ze,dbl,ldbl,lm)
          dclev(i,j)=real(ldbl)

C**** calculate possible energy loss
          tpe0=-tflux1(i,j)*dtime*sha
          tpe1=0.
          do l=1,lm
            tpe0=tpe0+t_3d(i,j,l)*pk(l,i,j)*pdsig(l,i,j)*sha*mb2kg    
            tpe1=tpe1+t(l)*pk(l,i,j)*pdsig(l,i,j)*sha*mb2kg/(1.d0+deltx
     *           *q(l))
          end do
          ediff=(tpe1-tpe0)/((psf-pmtop)*sha*mb2kg)        ! C

          do l=1,lm
C**** correct virt pot t for energy error
            t(l) = t(l) - ediff*(1.d0+deltx*q(l))/pk(l,i,j)
            ! update 3-d t,q,e and km
            t0ijl=t_3d(i,j,l)
            tijl=t(l)/(1.d0+deltx*q(l))
            t_3d(i,j,l)=tijl
C**** moment variation to be added
cc            qmom(:,i,j,l)=qmomij(:,l)
cc            tmom(:,i,j,l)=tmomij(:,l)

            q_3d(i,j,l)=q(l)
            e_3d(l,i,j)=e(l)
            w2_3d(l,i,j)=2.*by3*e(l)
            ! t2_3d(l,i,j)=t2(l)  ! not in use
            km_3d(l,i,j)=km(l)
            ! ACCUMULATE DIAGNOSTICS for t and q
            AJL(J,L,JL_TRBHR)=AJL(J,L,JL_TRBHR)
     2                 +(tijl-t0ijl)*PK(L,I,J)*PLIJ(L,I,J)
            AJL(J,L,JL_TRBDLHT)=AJL(J,L,JL_TRBDLHT)
     2                 +(q(l)-q0(l))*PDSIG(L,I,J)*LHE/SHA
            AJL(J,L,JL_TRBKE)=AJL(J,L,JL_TRBKE)+e(l)
#ifdef TRACERS_ON
            do nx=1,nta
              n=ntix(nx)
              tajln(j,l,jlnt_turb,n)=tajln(j,l,jlnt_turb,n) +
     &             (trij(l,nx)-tr0ij(l,nx))*am(l,i,j)*dxyp(j)
              trm(i,j,l,n)=trij(l,nx)*am(l,i,j)*dxyp(j)
cc            trmom(:,i,j,l,n)=trmomij(:,l,nx)
            end do
#endif
          end do

          ! Write out diagnostics if at selected grid point:

          if (call_diag.and.(i.eq.itest).and.(j.eq.jtest)) then
            call dout(ze,dz,u,v,t,q,ke,dtdz,dqdz,an2,as2
     &         ,wt_nl,wq_nl,kh,km,e,lscale
     &         ,uflx,vflx,tvflx,qflx,dbl,ldbl,i,j,lm)
          endif

        end do loop_i_tq
      end do loop_j_tq
!$OMP  END PARALLEL DO

      ! integrate U,V equations at bgrids
      call ave_uv_to_bgrid(uflux1,vflux1,uflux_bgrid,vflux_bgrid,1)
      call ave_s_to_bgrid(km_3d,dz_3d,dze_3d,rho_3d,rhoe_3d,
     &  km_3d_bgrid,dz_3d_bgrid,dze_3d_bgrid,rho_3d_bgrid,
     &  rhoe_3d_bgrid,lm)

!$OMP  PARALLEL DO PRIVATE (L,I,J,u,v,rho,rhoe,u0,v0,km,dze,dz,
!$OMP*   bydzerho,rhoebydz,uflx,vflx,p4,flux_bot,flux_top)
!$OMP*    SHARED(dtime)
!$OMP*    SCHEDULE(DYNAMIC,2)
      loop_j_uv: do j=J_0STG, J_1STG
        loop_i_uv: do i=1,im

          do l=1,lm
            u(l)=u_3d(i,j,l)
            v(l)=v_3d(i,j,l)
            rho(l)=rho_3d_bgrid(l,i,j)
            rhoe(l)=rhoe_3d_bgrid(l,i,j)
            u0(l)=u(l)
            v0(l)=v(l)
            u_3d_old(l,i,j)=u(l)
            v_3d_old(l,i,J)=v(l)
            km(l)=km_3d_bgrid(l,i,j)
            dze(l)=dze_3d_bgrid(l,i,j)
            dz(l)=dz_3d_bgrid(l,i,j)
            bydzerho(l)=1.d0/(dze(l)*rho(l))
          end do
          rhoebydz(1)=0.
          do l=1,lm-1
            rhoebydz(l+1)=rhoe(l+1)/dz(l)
          end do

          uflx=uflux_bgrid(i,j)
          vflx=vflux_bgrid(i,j)

          ! integrate differential eqns for U and V
          do l=2,lm-1
              p4(l)=0.d0
          end do
          flux_bot=rhoe(1)*uflx
          flux_top=0.d0
          call de_solver_main(u,u0,km,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)
          flux_bot=rhoe(1)*vflx
          call de_solver_main(v,v0,km,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)

          ! update u and v
          do l=1,lm
            u_3d(i,j,l)= u(l)
            v_3d(i,j,l)= v(l)
          end do

        end do loop_i_uv
      end do loop_j_uv
!$OMP  END PARALLEL DO
      CALL CHECKSUM(GRID, u_3d,__LINE__,__FILE__,STGR=.true.)
      CALL CHECKSUM(GRID, v_3d,__LINE__,__FILE__,STGR=.true.)

      CALL HALO_UPDATE(grid, u_3d, from=NORTH)

C**** Use halo values of plij from southern neighbor to complete that
C     neighbor's update of AJL in this process (i.e. using PLIJ 
C     contribution at J_0-1 to finish accumulation into AJL(J_0,:,:) )
        CALL HALO_UPDATE_COLUMN(grid, PLIJ, from=SOUTH)
 
! ACCUMULATE DIAGNOSTICS for u and v
C*** Adapt  use of idjj stencil to do correct calculation in distributed
C    parallel version.

C Contribution at southern border        
      IF (HAVE_SOUTH_POLE) THEN
        J=1
      ELSE
        J=J_0H
      END IF
      KMAX=KMAXJ(J)
      DO I=1,IMAXJ(J)
        DO L=1,LM
          DO K=1,KMAX
            RAK=RAVJ(K,J)
            IDIK=IDIJ(K,I,J)
            IDJK=IDJJ(K,J)
            IF (IDJK >= J_0)
     &        AJL(IDJK,L,JL_DAMDC)=AJL(IDJK,L,JL_DAMDC)
     &        +(U_3d(IDIK,IDJK,L)-u_3d_old(L,IDIK,IDJK))*PLIJ(L,I,J)*RAK
          ENDDO
        ENDDO
      ENDDO

      do j=J_0S,J_1-1
        KMAX=KMAXJ(J)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,KMAX
              RAK=RAVJ(K,J)
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              AJL(IDJK,L,JL_DAMDC)=AJL(IDJK,L,JL_DAMDC)
     &        +(U_3d(IDIK,IDJK,L)-u_3d_old(L,IDIK,IDJK))*PLIJ(L,I,J)*RAK
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C**** J_1 computation (add only k=1,2 contributions to ajl when j_1 is
C     NOT the north pole -- k=3,4 contr. to be added by southern neighbor).
      if (HAVE_NORTH_POLE) THEN
        J=JM 
        KMAX=KMAXJ(J)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,KMAX
              RAK=RAVJ(K,J)
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              AJL(IDJK,L,JL_DAMDC)=AJL(IDJK,L,JL_DAMDC)
     &        +(U_3d(IDIK,IDJK,L)-u_3d_old(L,IDIK,IDJK))*PLIJ(L,I,J)*RAK
            ENDDO
          ENDDO
        ENDDO
      else
        J=J_1
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,2
              RAK=RAVJ(K,J)
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              AJL(IDJK,L,JL_DAMDC)=AJL(IDJK,L,JL_DAMDC)
     &        +(U_3d(IDIK,IDJK,L)-u_3d_old(L,IDIK,IDJK))*PLIJ(L,I,J)*RAK
            ENDDO
          ENDDO
        ENDDO
      end if

C**** Save additional changes in KE for addition as heat later
!$OMP  PARALLEL DO PRIVATE (L,I,J)
      DO L=1,LM
      DO J=J_0STG, J_1STG
      DO I=1,IM
        DKE(I,J,L)=DKE(I,J,L)+0.5*(U_3d(I,J,L)*U_3d(I,J,L)
     *       +V_3d(I,J,L)*V_3d(I,J,L)-U_3d_old(L,I,J)*U_3d_old(L,I,J)
     *       -V_3d_old(L,I,J)*V_3d_old(L,I,J))
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

      return
      end subroutine atm_diffus

      subroutine getdz(tv,dz,dze,rho,rhoe,tvsurf,lm)
!@sum  getdz computes the 3-d finite difference dz and dze
!@+    as well as the 3-d density rho and rhoe
!@+    called at the primary grid (A-grid)
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var  tv virtual potential temp. referenced at 1 mb
!@var  dz main grid spacing
!@var  dze edge grid spacing
!@var  rho,rhoe air density at the main/edge grids
!@var  tvsurf surface virtual temperature
!@var  im,jm,lm 3-d grids

      !
      !     Grids:
      !
      !                   -------------------------  lm+1
      !                lm  - - - - - - - - - - - - -
      !                   -------------------------  lm
      !                l+1 - - - - - - - - - - - - -
      !                    -------------------------  l+1
      !     (main)     l   - - - - - - - - - - - - -
      !                    -------------------------  l     (edge)
      !                l-1 - - - - - - - - - - - - -
      !                    -------------------------  l-1
      !                2   - - - - - - - - - - - - -
      !                    -------------------------    2
      !                1   - - - - - - - - - - - - -
      !                    -------------------------    1
      !           dz(l,i,j) z(l+1,i,j) - z(l,i,j)
      !           dze(l,i,j) ze(l+1,i,j) - ze(l,i,j)
      !           rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
      !           rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dze(l,i,j))
      !
      !     at main: u,v,tv,q,ke
      !     at edge: e,lscale,km,kh,gm,gh
      !
      USE CONSTANT, only : grav,rgas
      USE GEOM, only : imaxj
      USE MODEL_COM, only : im,jm
      USE DYNAMICS, only : pmid,pk,pedn
      USE DOMAIN_DECOMP, ONLY : grid

      implicit none

      integer, intent(in) :: lm
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo), 
     &        intent(in) :: tvsurf
      real*8, dimension(lm,im,grid%j_strt_halo:grid%j_stop_halo), 
     &        intent(in) :: tv
      real*8, dimension(lm,im,grid%j_strt_halo:grid%j_stop_halo), 
     &        intent(out) :: rho,rhoe,dz,dze

      real*8 :: temp0,temp1,temp1e,pl1,pl,pl1e,ple,plm1e
      integer :: i,j,l  !@var i,j,l loop variable

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG

C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP
      J_0STG = grid%J_STRT_STGR
      J_1STG = grid%J_STOP_STGR

      !@ temp0 virtual temperature (K) at (i,j) and SIG(l)
      !@ temp1 virtual temperature (K) at (i,j) and SIG(l+1)
      !@ temp1e average of temp0 and temp1
!$OMP  PARALLEL DO PRIVATE (J,I,L,pl1,pl,pl1e,ple,temp0,temp1,temp1e,
!$OMP*  plm1e)
!$OMP*    SCHEDULE(DYNAMIC,2)
      do j=J_0, J_1
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
            dze(l,i,j)=-(rgas/grav)*temp0 *log(pl1e/ple)
            rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
            rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dze(l,i,j))
            if(l.eq.1) then
              rhoe(1,i,j)=100.d0*ple/(tvsurf(i,j)*rgas)
            endif
            if(l.eq.lm-1) then
              plm1e=pedn(lm+1,i,j)
              dze(lm,i,j)=-(rgas/grav)*temp1 *log(plm1e/pl1e)
              dz(lm,i,j)=0.
              rho(lm,i,j)=100.d0*(pl1e-plm1e)/(grav*dze(lm,i,j))
            endif

          end do
        end do
      end do
!$OMP  END PARALLEL DO

      return
      end subroutine getdz

      subroutine dout(ze,dz,u,v,t,q,ke,dtdz,dqdz,an2,as2
     &      ,wt_nl,wq_nl,kh,km,e,lscale
     &      ,uflx,vflx,tvflx,qflx,dbl,ldbl,i,j,n)
!@sum dout writes out diagnostics at (i,j)
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var p  pressure at main grid z
!@var pe  pressure at secondary grid ze
!@var u  west-east   velocity component
!@var v  south-north velocity component
!@var t  virt. pot. temperature (referenced to 1mb)
!@var q  specific humidity
!@var e  turbulent kinetic energy
!@var ze hight at edge lyer
!@var dz(j) z(j+1)-z(j)
!@var dze(j) ze(j+1)-ze(j)
!@var as2 dudz^2+dvdz^2
!@var dtdz z-derivative of t at edge grid ze
!@var dqdz z-derivative of q at edge grid ze
!@var an2 g_alpha*dtdz, brunt-vasala frequency
!@var km turbulent viscosity for u and v equations
!@var kh turbulent diffusivity for t,q
!@var ke turbulent diffusivity for e
!@var lscale  turbulent length scale
!@var uflx momentun flux -uw at surface, ze(1)
!@var vflx momentun flux -vw at surface, ze(1)
!@var tvflx heat flux -wt at surface, ze(1)
!@var qflx moisture flux -wq at surface, ze(1)
!@var i/j horizontal location at which the output is written
!@var n number of vertical main layers

      USE DYNAMICS, only : pmid,pedn,pk,pek

      implicit none

      integer, intent(in) :: ldbl,i,j,n 
      real*8, dimension(n), intent(in) ::
     &   dz,u,v,t,q,ke,dtdz,dqdz,an2,as2
     &   ,wt_nl,wq_nl,kh,km,e,lscale
      real*8, dimension(n+1), intent(in) :: ze
      real*8, intent(in) :: uflx,vflx,tvflx,qflx,dbl

      real*8, dimension(n) :: p,pe
      real*8 :: z,ri,wt_lcl,wq_lcl
      integer :: l  !@var l loop variable

      Write (67,1100) "i=",i,"j=",j
      Write (67,1200) "uflx=",uflx,"vflx=",vflx
      Write (67,1200) "tvflx=",tvflx*pek(1,i,j),"qflx=",qflx
      Write (67,1250) "ldbl=",ldbl,"dbl=",dbl

      ! pressure at main and edge layers

      do l=1,n
          p(l)=100.*pmid(l,i,j)
          pe(l)=100.*pedn(l,i,j)
      end do

      ! fields at layer center

      write (67,1300)
      do l=1,n
        z=.5d0*(ze(l)+ze(l+1))
        write (67,2000) l,p(l),z,dz(l),u(l),v(l),t(l)*pk(l,i,j),q(l)
     &                  ,ke(l)
      end do
      write (67,*)

      ! fields at layer edge

      write (67,1400)
      do l=1,n
        wt_lcl=-kh(l)*dtdz(l)*pek(l,i,j)
        wq_lcl=-kh(l)*dqdz(l)
        ri=an2(l)/as2(l)
        write (67,2000) l,pe(l),ze(l),wt_lcl,wt_nl(l)*pek(l,i,j)
     &    ,wq_lcl,wq_nl(l),kh(l),km(l),lscale(l),ri,e(l)
      end do
      write (67,1500)

      return
1100  format(3(2x,a,i6))
1200  format(2(2x,a,1pe14.4))
1250  format(2x,a,i6,2x,a,1pe14.4)
1300  format (' ',' l',1x,
     1            '     p     ',1x,
     2            '     z     ',1x,'     dz    ',1x,'     u     ',1x,
     3            '     v     ',1x,'     t     ',1x,'     q     ',1x,
     4            '     ke    ',1x,'           ',1x,'           ',1x,
     5            '           ',/)
        write (67,2000) l,pe(l),ze(l),wt_lcl,wt_nl(l),wq_lcl,wq_nl(l)
     2    ,kh(l),km(l),lscale(l),ri,e(l)
1400  format (' ',' l',1x,
     2            '  p (edge) ',1x,
     2            '  z (edge) ',1x,'  wt_lcl   ',1x,'   wt_nl   ',1x,
     3            '   wq_lcl  ',1x,'  wq_nl    ',1x,'     kh    ',1x,
     4            '     km    ',1x,'  lscale   ',1x,'     ri    ',1x,
     5            '     e     ',/)
1500  format(132('-'))
2000  format (1x,i2,1x,1pe11.4,9(1x,1pe11.4),1x,1pe10.3)
      end subroutine dout

      subroutine de_solver_main(x,x0,p1,p4,
     &    rhoebydz,bydzerho,flux_bot,flux_top,dtime,n,qlimit)
!@sum differential eqn solver for x using tridiagonal method
!@+   d/dt x = d/dz (P1 d/dz x) + P4
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var x the unknown to be solved (at main drid)
!@var x0 x at previous time step
!@var p1,p4 coeff. of the d.e.
!@var rhoebydz(j) rhoe(j+1)/dz(j)
!@var bydzerho(j) 1/(dze(j)*rho(j))
!@var flux_bot flux from the bottom
!@var flux_top flux at the top
!@var dtime time step
!@var n number of main grid
!@var qlimit true if tracer must be positive definite
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) ::
     &        x0,p1,p4,rhoebydz,bydzerho
      real*8, intent(in) :: flux_bot,flux_top,dtime
      real*8, dimension(n), intent(out) :: x

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable
      logical, intent(in) :: qlimit 

      !sub(j)*x_jm1_kp1+dia(j)*x_j_kp1+sup(j)*x_jp1_kp1 = rhs(j)
      !k refers to time step, j refers to main grid
      !except p1(j) which is defined on the edge grid

c GISS-ESMF EXCEPTIONAL CASE
c This looks like a loop which is always over levels for all
c calls to this routine.
      do j=2,n-1
          sub(j)=-dtime*p1(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*p1(j+1)*rhoebydz(j+1)*bydzerho(j)
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=x0(j)+dtime*p4(j)
          if (qlimit.and. rhs(j).lt.0) rhs(j)=0.  ! prevent roundoff error
      end do

      ! Lower boundary conditions(x=T) :
      ! d/dt T = -(1/rho)*d/dz(rho*wt)
      ! d/dt T = (T(1)-T0(1))/dtime
      ! -d/dz(rho*wt)=-(rhoe(2)*wt(2)-rhoe(1)*wt(1))/dze(1)
      ! wt(2)=-p1(2)*(T(2)-T(1))/dz(1)+wt_nl(2)
      ! wt(1)=-tvflx, therefore,flux_bot(the quantity used below)
      !       flux_bot=rhoe(1)*tvflx+rhoe(2)*wt_nl(2)
      
      alpha=dtime*p1(2)*rhoebydz(2)*bydzerho(1)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=x0(1)-dtime*bydzerho(1)*flux_bot

      ! Upper boundary conditions:

      ! Upper boundary conditions(x=T) :
      ! d/dt T = -(1/rho)*d/dz(rho*wt)
      ! d/dt T = (T(n)-T0(n))/dtime
      ! -d/dz(rho*wt)=-(rhoe(n+1)*wt(n+1)-rhoe(n)*wt(n))/dze(n)
      ! wt(n)=-p1(n)*(T(n)-T(n-1))/dz(n-1)+wt_nl(n)
      ! wt(n+1)=0, therefore,flux_top
      ! flux_top=0.

      alpha=dtime*p1(n)*rhoebydz(n)*bydzerho(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=x0(n)+dtime*bydzerho(n)*flux_top

      call tridiag(sub,dia,sup,rhs,x,n)

      return
      end subroutine de_solver_main

      subroutine de_solver_edge(x,x0,p1,p3,p4,
     &    rhobydze,bydzrhoe,x_surf,dtime,n)
!@sum differential eqn solver for x using tridiagonal method
!@+   d/dt x = d/dz (P1 d/dz x) - P3 x + P4
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var x the unknown to be solved (at edge drid)
!@var x0 x at previous time step
!@var p1,p3,p4 coeff. of the d.e.
!@var rhobydze(j) rho(j)/dze(j)
!@var bydzrhoe(j) 1/(dz(j-1)*rhoe(j))
!@var x_surf surface value of x
!@var dtime time step
!@var n number of vertical edge grid

      use CONSTANT, only : teeny
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) ::
     &        x0,p1,p3,p4,rhobydze,bydzrhoe
      real*8, intent(in) :: x_surf,dtime
      real*8, dimension(n), intent(out) :: x

      real*8, dimension(n) :: sub,dia,sup,rhs
      integer :: j  !@var j loop variable

      ! sub(j)*x_jm1_kp1+dia(j)*x_j_kp1+sup(j)*x_jp1_kp1 = rhs(j)
      ! k refers to time step, j refers to edge grid
      ! except p1(j) which is defined on the main grid

      do j=2,n-1
          sub(j)=-dtime*p1(j-1)*rhobydze(j-1)*bydzrhoe(j)
          sup(j)=-dtime*p1(j)*rhobydze(j)*bydzrhoe(j)
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*p3(j)
          rhs(j)=x0(j)+dtime*p4(j)
      end do

      ! Boundary conditions:

      dia(1)=1.d0
      sup(1)=0.d0
      rhs(1)=x_surf

      sub(n)=-1.d0
      dia(n)=1.d0
      rhs(n)=0.d0

      call tridiag(sub,dia,sup,rhs,x,n)
      ! assuming x is a non-negative scalar:
      do j=1,n
         if(x(j).lt.teeny) x(j)=teeny
      end do

      return
      end subroutine de_solver_edge

      subroutine ave_uv_to_agrid(u,v,u_a,v_a,lm)
!@sum Computes u_a,v_a from u and v
!@var u x-component at secondary grids (B_grid)
!@var v y-component at secondary grids (B_grid)
!@var u_a x-component at primary grids (A_grid)
!@var v_a y-component at primary grids (A_grid)
!@auth Ye Cheng
!@ver  1.0

      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP, only : grid,get,NORTH, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP, only : halo_update,checksum, checksum_column
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none

      integer, intent(in) :: lm
      real*8, dimension(im, grid%j_strt_halo:grid%j_stop_halo,lm),
     &                  intent(inout) ::  u,v
      real*8, dimension(lm,im, grid%j_strt_halo:grid%j_stop_halo), 
     &                  intent(out) :: u_a,v_a

      real*8, dimension(im) :: ra
      integer, dimension(im) :: idj
      real*8 :: HEMI,u_t,v_t,rak,ck,sk,uk,vk
      integer :: i,j,l,k,idik,idjk,kmax

      integer :: J_0S,J_1S
      logical :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      call get(grid, J_STRT_SKP=J_0S,   J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
!     polar boxes

C**** Update halos of U and V 
      CALL CHECKSUM(grid, u, __LINE__,__FILE__,STGR=.true.)
      CALL CHECKSUM(grid, v, __LINE__,__FILE__,STGR=.true.)

      CALL HALO_UPDATE(grid,u, from=NORTH)
      CALL HALO_UPDATE(grid,v, from=NORTH)


      if (HAVE_SOUTH_POLE) then
        J=1
        KMAX=KMAXJ(J)
        HEMI=-1.
!$OMP  PARALLEL DO PRIVATE (I,L,u_t,v_t,K,IDIK,IDJK,RAK,ck,sk,uk,vk)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
          END DO
        END DO
!$OMP  END PARALLEL DO
      end if              !south pole
!
      if (HAVE_NORTH_POLE) then
        J=JM
        KMAX=KMAXJ(J)
        HEMI=1.
!$OMP  PARALLEL DO PRIVATE (I,L,u_t,v_t,K,IDIK,IDJK,RAK,ck,sk,uk,vk)
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
          END DO
        END DO
!$OMP  END PARALLEL DO
      end if                !north pole

!     non polar boxes
C**** Update halos of u and v. (Needed bcs. IDJJ(3:4,J_1S)=J_1S+1)
C     ---> done by calling routine...

      CALL CHECKSUM(grid, u, __LINE__,__FILE__,STGR=.true.)
      CALL CHECKSUM(grid, v, __LINE__,__FILE__,STGR=.true.)

      CALL HALO_UPDATE(grid, u, FROM=NORTH)
      CALL HALO_UPDATE(grid, v, FROM=NORTH)
!$OMP  PARALLEL DO PRIVATE (J,I,L,u_t,v_t,K,KMAX,IDJ,RA,IDIK,IDJK,RAK)
!$OMP*    SCHEDULE(DYNAMIC,2)
      DO J=J_0S,J_1S
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAPJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              u_t=u_t+u(IDIK,IDJK,L)*RAK
              v_t=v_t+v(IDIK,IDJK,L)*RAK
            END DO
            u_a(l,i,j)=u_t
            v_a(l,i,j)=v_t
          END DO
        END DO
      END DO
!$OMP  END PARALLEL DO
C****
      return
      end subroutine ave_uv_to_agrid

      subroutine ave_uv_to_bgrid(u_a,v_a,u,v,lm)
!@sum Computes u and v from u_a and v_a
!@var u_a x-component of wind at primary grids (A_grid)
!@var v_a y-component of wind at primary grids (A_grid)
!@var u x-component of wind at secondary grids (B_grid)
!@var v y-component of wind at secondary grids (B_grid)
!@auth Ye Cheng
!@ver  1.0

      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP, only : grid, get
      USE DOMAIN_DECOMP, only : halo_update_column, checksum_column
      USE DOMAIN_DECOMP, only : halo_update       , checksum
      USE DOMAIN_DECOMP, only : NORTH, SOUTH
      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj,cosiv,siniv

      implicit none

      integer, intent(in) :: lm
      real*8, dimension(lm,im,grid%j_strt_halo:grid%j_stop_halo)  ::
     &          u_a,v_a
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm), 
     &        intent(out) :: u,v

      real*8, dimension(im) :: ra
      integer, dimension(im) :: idj
      real*8 :: HEMI,rak,ck,sk,uk,vk
      integer :: i,j,l,k,idik,idjk,kmax
      integer :: j_0, j_1, j_0s, j_1S  
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call get(grid, J_STRT=J_0,   J_STOP=J_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )


      u=0.d0; v=0.d0
      CALL HALO_UPDATE_COLUMN(grid,U_A, from=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid,V_A, from=SOUTH)

      CALL CHECKSUM(grid, U,__LINE__,__FILE__,STGR=.true.)
      CALL CHECKSUM(grid, V,__LINE__,__FILE__,STGR=.true.)

      if (HAVE_SOUTH_POLE) then
        J=1
        KMAX=KMAXJ(J)
        HEMI=-1.
        DO I=1,IMAXJ(J)
        DO L=1,LM
        DO K=1,KMAX
          IDIK=IDIJ(K,I,J)
          IDJK=IDJJ(K,J)
          RAK=RAVJ(K,J)
          ck=cosiv(k)
          sk=siniv(k)
          uk=u_a(L,I,J)
          vk=v_a(L,I,J)
          U(IDIK,IDJK,L)=U(IDIK,IDJK,L)+RAK*(UK*CK+VK*SK*HEMI)
          V(IDIK,IDJK,L)=V(IDIK,IDJK,L)+RAK*(VK*CK-UK*SK*HEMI)
        END DO
        END DO
        END DO
      Else
! carry over from southern neighbor
        J=J_0-1
        DO K=3,4
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=3,4
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              U(IDIK,IDJK,L)=U(IDIK,IDJK,L)+RAK*U_A(L,I,J)
              V(IDIK,IDJK,L)=V(IDIK,IDJK,L)+RAK*V_A(L,I,J)
            END DO
          END DO
        END DO
      end if ! south pole

      if (HAVE_NORTH_POLE) then
        J=JM
        KMAX=KMAXJ(J)
        HEMI=1.
        DO I=1,IMAXJ(J)
        DO L=1,LM
        DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAVJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u_a(L,I,J)
              vk=v_a(L,I,J)
          U(IDIK,IDJK,L)=U(IDIK,IDJK,L)+RAK*(UK*CK+VK*SK*HEMI)
          V(IDIK,IDJK,L)=V(IDIK,IDJK,L)+RAK*(VK*CK-UK*SK*HEMI)
        END DO
        END DO
        END DO
      end if     !north pole     

!     non polar boxes (skip J_1 box --> computed below)
      DO J=J_0S,J_1-1
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
        DO L=1,LM
        DO K=1,KMAX
          IDIK=IDIJ(K,I,J)
          IDJK=IDJ(K)
          RAK=RA(K)
          U(IDIK,IDJK,L)=U(IDIK,IDJK,L)+RAK*U_A(L,I,J)
          V(IDIK,IDJK,L)=V(IDIK,IDJK,L)+RAK*V_A(L,I,J)
        END DO
        END DO
        END DO
      END DO

c**** J_1 box
      if (.not. HAVE_NORTH_POLE) then
        J=J_1
        DO K=1,2
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,2
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              U(IDIK,IDJK,L)=U(IDIK,IDJK,L)+RAK*U_A(L,I,J)
              V(IDIK,IDJK,L)=V(IDIK,IDJK,L)+RAK*V_A(L,I,J)
            END DO
          END DO
        END DO
      end if

C**** Update halos of U_A and V_A
      CALL CHECKSUM_COLUMN(grid, U_A, __LINE__,__FILE__,STGR=.true.)
      CALL CHECKSUM_COLUMN(grid, V_A, __LINE__,__FILE__,STGR=.true.)
                                                               
      CALL CHECKSUM(grid, U,__LINE__,__FILE__,STGR=.true.)
      CALL CHECKSUM(grid, V,__LINE__,__FILE__,STGR=.true.)


      return
      end subroutine ave_uv_to_bgrid

      subroutine ave_s_to_bgrid(s1,s2,s3,s4,s5,
     &  sb1,sb2,sb3,sb4,sb5,lm)
!@sum Computes s_b from s
!@var s  scalar at primary grid (A_grid)
!@var sb scalar at secondary grid (B_grid)
!@auth Ye Cheng
!@ver  1.0

      USE MODEL_COM, only : im,jm
      USE GEOM, only : imaxj,idij,idjj,kmaxj,ravj
      USE DOMAIN_DECOMP, only : grid,get,NORTH,SOUTH
      USE DOMAIN_DECOMP, only : HALO_UPDATE_COLUMN, CHECKSUM_COLUMN

      implicit none

      integer, intent(in) :: lm
      real*8, dimension(lm,im,grid%j_strt_halo:grid%j_stop_halo) 
     &         ::  s1,s2,s3,s4,s5
      real*8, dimension(lm,im,grid%j_strt_halo:grid%j_stop_halo), 
     &        intent(out) :: sb1,sb2,sb3,sb4,sb5

      real*8, dimension(im) :: ra
      integer, dimension(im) :: idj
      integer :: idik,idjk,kmax
      integer :: i,j,l,k
      real*8 :: rak

C*** GET USEFUL DOMAIN DECOMPOSITION PARAMETERS
      integer :: J_0,J_1
      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      SB1=0.;SB2=0.;SB3=0.;SB4=0.;SB5=0.
C*** Adapt  use of idjj stencil to do correct calculation in distributed
C    parallel version.
C   ......first skip j_1 (do later below)

      CALL HALO_UPDATE_COLUMN(grid, S1, from=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, S2, from=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, S3, from=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, S4, from=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, S5, from=SOUTH)

      If (.not. HAVE_SOUTH_POLE) THEN
        j = J_0-1
        KMAX=KMAXJ(J)
        DO K=3,4
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=3,4
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              SB1(L,IDIK,IDJK)=SB1(L,IDIK,IDJK)+RAK*S1(L,I,J)
              SB2(L,IDIK,IDJK)=SB2(L,IDIK,IDJK)+RAK*S2(L,I,J)
              SB3(L,IDIK,IDJK)=SB3(L,IDIK,IDJK)+RAK*S3(L,I,J)
              SB4(L,IDIK,IDJK)=SB4(L,IDIK,IDJK)+RAK*S4(L,I,J)
              SB5(L,IDIK,IDJK)=SB5(L,IDIK,IDJK)+RAK*S5(L,I,J)
            END DO
          END DO
        END DO
      End If        

      DO j=J_0,J_1-1
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              SB1(L,IDIK,IDJK)=SB1(L,IDIK,IDJK)+RAK*S1(L,I,J)
              SB2(L,IDIK,IDJK)=SB2(L,IDIK,IDJK)+RAK*S2(L,I,J)
              SB3(L,IDIK,IDJK)=SB3(L,IDIK,IDJK)+RAK*S3(L,I,J)
              SB4(L,IDIK,IDJK)=SB4(L,IDIK,IDJK)+RAK*S4(L,I,J)
              SB5(L,IDIK,IDJK)=SB5(L,IDIK,IDJK)+RAK*S5(L,I,J)
            END DO
          END DO
        END DO
      END DO
C**** J_1 computation (add only k=1,2 contributions to ajl when j_1 is
C     NOT the north pole -- k=3,4 contr. to be added by southern neighbor).
      if (HAVE_NORTH_POLE) THEN
        J=JM
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              SB1(L,IDIK,IDJK)=SB1(L,IDIK,IDJK)+RAK*S1(L,I,J)
              SB2(L,IDIK,IDJK)=SB2(L,IDIK,IDJK)+RAK*S2(L,I,J)
              SB3(L,IDIK,IDJK)=SB3(L,IDIK,IDJK)+RAK*S3(L,I,J)
              SB4(L,IDIK,IDJK)=SB4(L,IDIK,IDJK)+RAK*S4(L,I,J)
              SB5(L,IDIK,IDJK)=SB5(L,IDIK,IDJK)+RAK*S5(L,I,J)
            END DO
          END DO
        END DO
      else
        J=J_1
        DO K=1,2
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAVJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            DO K=1,2
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              SB1(L,IDIK,IDJK)=SB1(L,IDIK,IDJK)+RAK*S1(L,I,J)
              SB2(L,IDIK,IDJK)=SB2(L,IDIK,IDJK)+RAK*S2(L,I,J)
              SB3(L,IDIK,IDJK)=SB3(L,IDIK,IDJK)+RAK*S3(L,I,J)
              SB4(L,IDIK,IDJK)=SB4(L,IDIK,IDJK)+RAK*S4(L,I,J)
              SB5(L,IDIK,IDJK)=SB5(L,IDIK,IDJK)+RAK*S5(L,I,J)
            END DO
          END DO
        END DO
      end if

C**** Use halo values of S1 -- S5 from southern neighbor to complete that
C     neighbor's updates of SB1 -- SB5  in this process (i.e. using S[n]
C     contribution at J_0 -1 to finish accumulation into SB[n](:,:,J_0) )
C**** Halo updates from the south
        CALL CHECKSUM_COLUMN(grid, SB1,  __LINE__, __FILE__)
        CALL CHECKSUM_COLUMN(grid, SB2,  __LINE__, __FILE__)
        CALL CHECKSUM_COLUMN(grid, SB3,  __LINE__, __FILE__)
        CALL CHECKSUM_COLUMN(grid, SB4,  __LINE__, __FILE__)
        CALL CHECKSUM_COLUMN(grid, SB5,  __LINE__, __FILE__)


      return
      end subroutine ave_s_to_bgrid

      subroutine apply_fluxes_to_atm
!@sum dummy subroutine - replaces the real one needed by DRYCNV
!@auth I. Aleinov
!@ver  1.0
      return
      end subroutine apply_fluxes_to_atm

      subroutine zze(dze,ze,n)
!@sum finds the layer edge height ze
!@var  ze vertical coordinate associated with SIGE(l)
!@var  dze(l) ze(l+1) - ze(l)
!@var  n number of layers

      USE SOCPBL, only : zgs

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: dze
      real*8, dimension(n+1), intent(out) :: ze
      integer :: j  !@var j loop variable

c GISS-ESMF EXCEPTIONAL CASE
c This looks like a loop which is always over levels for all
c calls to this routine.
      ze(1)=zgs
      do j=2,n+1
          ze(j)=ze(j-1)+dze(j-1)
      end do
      return
      end subroutine zze


      subroutine l_gcm(ze,lscale,n)
!@sum l_gcm calculates the turbulent length scale
!@auth Ye Cheng/G. Hartke
!@ver  1.0
!@var ze height (meters) of layer edge
!@var lscale turbulent length scale
!@var n number of layers

      USE CONSTANT, only : teeny
      USE SOCPBL, only : kappa

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n+1), intent(in) :: ze
      real*8, dimension(n), intent(out) :: lscale

      real*8 :: l0,kz
      integer :: j  !@var j loop variable

c GISS-ESMF EXCEPTIONAL CASE
c This looks like a loop which is always over levels for all
c calls to this routine.
      !@ Ref: Holtslag and Boville 1993, J. Climate, 6, 1825-1842.
      do j=1,n
        l0=30.+270.*exp(1.-1.d-3*ze(j))
        kz=kappa*ze(j)
        lscale(j)=l0*kz/(l0+kz)
      end do

      return
      end subroutine l_gcm

      subroutine k_gcm(tvflx,qflx,ustar,wstar,dbl
     &  ,ze,lscale,e,qturb,an2,as2,dtdz,dqdz,dudz,dvdz
     &  ,kh,km,ke,wt,wq,uw,vw,wt_nl,wq_nl
#ifdef TRACERS_ON
     &  ,trflx,wc_nl,nta
#endif
     &  ,n)

!@sum k_gcm computes the turbulent stability functions Km, Kc
!@+   and the non-local part of the fluxes
!@auth  Ye Cheng/G. Hartke
!@ver   1.0
!@var tvflx virtual potential temperature flux at surface
!@var qflx moisture flux at surface
!@var ustar friction velocity
!@var wstar Deardorff convective velocity scale
!@var dbl the height of the PBL (real*8, in meters)
!@var kh turbulent diffusivity for scalars (heat, moisture,...)
!@var ze height (in meters) of layer edges
!@var lscale turbulent length scale
!@var e turbulent kinetic energy
!@var qturb sqrt(2*e)
!@var an2 g_alpha*dTdz
!@var as2 (dUdz)**2+(dVdz)**2
!@var dtdz dt/dz
!@var dqdz dq/dz
!@var dudz du/dz
!@var dvdz dv/dz
!@var km turbulent diffusivity for momentun
!@var ke turbulent diffusivity for e
!@var wt,wq,uw,vw turbulent fluxes
!@var wt_nl non-local part of heat flux wt
!@var wq_nl non-local part of moisture flux wq
!@var wc_nl non-local part of tracer flux
!@var n number of layers

      USE CONSTANT, only : teeny,by3
      USE SOCPBL, only : kappa,prt,ghmin,ghmax,d1,d2,d3,d4,d5
     &                  ,s0,s1,s2,s4,s5,s6,b1,g5

      implicit none

      integer, intent(in) :: n
#ifdef TRACERS_ON
      integer, intent(in) :: nta
      real*8, dimension(n,nta),intent(out) :: wc_nl
      real*8, dimension(nta),intent(in) :: trflx
      integer nt
#endif

      real*8, intent(in) :: tvflx,qflx,ustar,wstar,dbl
      real*8, dimension(n), intent(in) :: lscale,e,qturb,an2,as2
     &        ,dtdz,dqdz,dudz,dvdz
      real*8, dimension(n+1), intent(in) :: ze
      real*8, dimension(n), intent(out) ::
     &  kh,km,ke,wt,wq,uw,vw,wt_nl,wq_nl

      real*8, parameter :: kmmin=1.5d-5,khmin=2.5d-5,kmax=600.d0
      real*8 :: tmp0,tmp,tau,gm,gh,gmmax,byden,sm,sh
     &    ,ustar2,wstar3,zzi,tau_pt,w2j
      integer :: j  !@var j loop variable
      
      !@ Ref: Holtslag and Moeng 1991, JAS, 48, 1690-1698.
      ustar2=ustar*ustar
      wstar3=wstar*wstar*wstar

      tmp0=-wstar/(g5*dbl)
      do j=1,n
          tau=b1*lscale(j)/(qturb(j)+teeny)
          gh=tau*tau*an2(j)
          gm=tau*tau*as2(j)
          if(gh.lt.ghmin) gh=ghmin
          if(gh.gt.ghmax) gh=ghmax
          gmmax=(1.+d1*gh+d3*gh*gh)/(d2+d4*gh)
          if(gm.gt.gmmax) gm=gmmax
          byden=1./(1.+d1*gh+d2*gm+d3*gh*gh+d4*gh*gm+d5*gm*gm)
          sm=(s0+s1*gh+s2*gm)*byden
          km(j)=min(max(tau*e(j)*sm,kmmin),kmax)
          if(ze(j).le.dbl) then
              tau_pt=tau/g5
              zzi=ze(j)/dbl
              tmp=(1.6d0*ustar2*(1.-zzi)+teeny)**1.5d0
     &           +1.2d0*wstar3*zzi*(1.-.9d0*zzi)**1.5d0
              w2j=tmp**(2.*by3)
              kh(j)=min(max(.5d0*w2j*tau_pt,khmin),kmax)
              tmp=tmp0*tau
              wt_nl(j)=tmp*tvflx
              wq_nl(j)=tmp*qflx
#ifdef TRACERS_ON
              do nt=1,nta
                wc_nl(j,nt)=tmp*trflx(nt)
              end do
#endif
          else
              sh=(s4+s5*gh+s6*gm)*byden
              kh(j)=min(max(tau*e(j)*sh,khmin),kmax)
              wt_nl(j)=0.
              wq_nl(j)=0.
#ifdef TRACERS_ON
              do nt=1,nta
                wc_nl(j,nt)=0.
              end do
#endif
          endif
          ke(j)=5.*km(j)
          wt(j) = -kh(j)*dtdz(j)+wt_nl(j)
          wq(j) = -kh(j)*dqdz(j)+wq_nl(j)
          uw(j) = -km(j)*dudz(j)
          vw(j) = -km(j)*dvdz(j)
      end do

      return
      end subroutine k_gcm

      subroutine e_gcm(tvflx,wstar,ustar,dbl,lmonin,ze,g_alpha
     &    ,an2,as2,lscale,e,n)
!@sum e_gcm finds e according to the parameterization of les data
!@+   (Moeng and Sullivan 1994) for ze<=dbl and using the giss
!@+   soc model (level 2) for ze>dbl
!@auth  Ye Cheng
!@ver   1.0
!@var (see subroutine k_gcm)
!@var lmonin Monin-Obukov length
!@var g_alpha grav*alpha
      USE CONSTANT, only : teeny,by3
      USE SOCPBL, only : kappa,emax,rimax,b1,c1,c2,c3,c4,c5

      implicit none

      integer, intent(in) :: n   !@var n  array dimension
      real*8, intent(in) :: tvflx,wstar,ustar,dbl
      real*8, dimension(n), intent(in)   :: g_alpha,an2,as2,lscale
      real*8, dimension(n+1), intent(in) :: ze
      real*8, dimension(n), intent(out) :: e
      real*8, intent(out) :: lmonin
      integer :: j !@var j loop variable
      real*8 :: ri,gm,aa,bb,cc,phi_m,tmp,wstar3,ustar3,zj,zeta

      ustar3=ustar*ustar*ustar
      wstar3=wstar*wstar*wstar
      lmonin=ustar3/(kappa*g_alpha(1)*tvflx)  ! tvflx=-(wtv)_0
      do j=1,n   ! Dyer 1974
        zj=ze(j)
        if(zj.le.dbl) then
          zeta=zj/lmonin
          if(zeta.ge.0.) then ! stable or neutral
            if(zeta.le.1.) then
              phi_m=1.+5.*zeta
            else
             phi_m=5.+zeta
            endif
          else                ! unstable
            phi_m=(1.-15.*zeta)**(-.25d0)
          endif
          tmp=.4d0*wstar3+ustar3*(dbl-zj)*phi_m/(kappa*zj)
          e(j)=min(max(tmp**(2.*by3),teeny),emax)
        else
          ri=an2(j)/max(as2(j),teeny)
          if(ri.lt.rimax) then
            aa=c1*ri*ri-c2*ri+c3
            bb=c4*ri+c5
            cc=2.d0
            if(abs(aa).lt.1d-8) then
              gm= -cc/bb
            else
              tmp=bb*bb-4.*aa*cc
              gm=(-bb-sqrt(tmp))/(2.*aa)
            endif
            tmp=0.5d0*(b1*lscale(j))**2*as2(j)/max(gm,teeny)
            e(j)=min(max(tmp,teeny),emax)
          else
            e(j)=teeny
          endif
        endif
      end do
      return
      end subroutine e_gcm

      subroutine find_pbl_top(e,ze,dbl,ldbl,n)
!@sum find_pbl_top finds the pbl top (at main level)
!@auth  Ye Cheng
!@ver   1.0
!@var e turbulent kinetic energy
!@var ze height at the edge level (meters)
!@var ldbl the (main) layer corresponding to top of pbl
!@var dbl the height (in meters) of the pbl, at main layer
!@+   this dbl is different from the dbl in module socpbl,
!@+   the latter is itype dependent
!@var n number of layers

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: e
      real*8, dimension(n+1), intent(in) :: ze
      real*8, intent(out) :: dbl
      integer, intent(out) :: ldbl

      real*8, parameter :: dbl_max=3000.,fraction = 0.1d0
      real*8 :: e1p    ! a fraction of e(1)
      integer :: l

      e1p=fraction*e(1)
      do l=2,n
        if (e(l).lt.e1p) exit
      end do
      ldbl=l-1
      dbl=.5d0*(ze(l-1)+ze(l))  ! dbl is at main layer (l-1)
      dbl=min(dbl,dbl_max)
    
      return
      end subroutine find_pbl_top
