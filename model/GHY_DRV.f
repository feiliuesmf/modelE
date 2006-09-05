#include "rundeck_opts.h"

c******************   TRACERS             ******************************
#ifdef TRACERS_ON
      module ghy_tracers

      use sle001, only :
#ifdef TRACERS_WATER
     &     tr_w,tr_wsn,trpr,tr_surf,ntg,ntgm,atr_evap,atr_rnff,atr_g
#endif
#ifdef TRACERS_SPECIAL_O18
     &     ,tr_name
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      use sle001,ONLY : aevap
#endif
      use tracer_com, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,Ntm_dust
#endif
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,trname
#endif
#ifdef TRACERS_DUST
     &     ,n_clay,imDust
#else
#ifdef TRACERS_MINERALS
     &     ,n_clayilli
#else
#ifdef TRACERS_QUARZHEM
     &     ,n_sil1quhe
#endif
#endif
#endif
      use trdiag_com, only : taijn=>taijn_loc,tij_surf
     *  ,taijs=>taijs_loc,ijts_isrc,jls_isrc,tajls=>tajls_loc
#ifdef TRACERS_WATER
     *     ,tij_evap,tij_grnd,tij_soil
#endif
#ifdef TRACERS_DRYDEP
     *     ,tij_drydep,tij_gsdep,itcon_dd,dtr_dd
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     ,jls_source
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,ijts_source,nDustEmij,nDustEmjl
     &     ,ijts_spec,nDustEv1ij,nDustEv2ij,nDustWthij
     &     ,jls_spec,nDustEv1jl,nDustEv2jl,nDustWthjl
#endif
#ifdef TRACERS_DUST
     &     ,nDustEm2ij,nDustEm2jl
#endif
      use fluxes, only : trsource,trsrfflx
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer,trprec
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,prec,pprec,pevap,dust_flux_glob,trs_glob
#ifdef TRACERS_DRYDEP
     &     ,depo_turb_glob,depo_grav_glob
#endif
#endif
#ifdef TRACERS_DUST
     &     ,dust_flux2_glob
#endif
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
#endif

#ifdef TRACERS_WATER
      use ghy_com, only : tr_wbare,tr_wvege,tr_wsn_ij
#endif

      use pbl_drv, only : trtop,trs,trsfac,trconstflx,ntx,ntix
#ifdef TRACERS_WATER
     *     ,tr_evap_max
#endif
#ifdef TRACERS_DRYDEP
     *     ,dep_vel,gs_vel
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     *     ,DMS_flux, ss1_flux, ss2_flux
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,dust_flux,dust_flux2,z,km,gh,gm,zhat,lmonin,wsubtke,wsubwd
     &     ,wsubwm,dust_event1,dust_event2,wtrsh
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only : EMIS_SOURCE
#endif
ccc extra stuff which was present in "earth" by default
#ifdef TRACERS_WATER
      use ghy_com, only : ngm,nlsn
      use constant, only : rhow
#endif
      use dynamics, only : byam
      use geom, only : bydxyp

      implicit none
      private

      public ghy_tracers_set_step
      public ghy_tracers_set_cell
      public ghy_tracers_save_cell

      real*8 totflux(ntm)
#ifdef TRACERS_WATER
!@var ntixw index array for tracers (shared by OMP threads)
      real*8 ntixw(ntm)
#endif
      integer, parameter :: itype=4

      common /ghy_tracers_tp/ totflux
!$OMP  THREADPRIVATE (/ghy_tracers_tp/)

      contains


      subroutine ghy_tracers_set_step
!@sum tracers stuff to be called at the beginning of ghy time step
!@+   i.e. before the i,j loop
ccc set i,j - independent stuff for tracers
      use model_com, only : itime
      implicit none
      integer n

      ntx=0
#ifdef TRACERS_WATER
      ntg=0
#endif
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          ntx=ntx+1
          ntix(ntx) = n
#ifdef TRACERS_WATER
          if ( tr_wd_TYPE(n) == nWATER ) then
            ntg = ntg + 1
            if ( ntg > ntgm ) call stop_model("ghy_drv: ntg > ntgm",255)
            ntixw(ntg) = n
#ifdef TRACERS_SPECIAL_O18
            tr_name(ntg) = trname(n)
#endif
          endif
#endif
        end if
      end do

      end subroutine ghy_tracers_set_step


      subroutine ghy_tracers_set_cell(i,j,qg,evap_max)
!@sum tracers code to be called before the i,j cell is processed
      use model_com, only : dtsrc
#ifdef TRACERS_WATER
      use sle001, only : qm1
#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: qg,evap_max
      integer n,nx,nsrc

C**** Set up tracers for PBL calculation if required
      do nx=1,ntx
        n = ntix(nx)
C**** Calculate first layer tracer concentration
        trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
      end do

ccc tracers variables
#ifdef TRACERS_WATER
      do nx=1,ntg
        n = ntixw(nx)
        ! prognostic vars
        tr_w(nx,0,1) = 0.d0
        tr_w(nx,1:ngm,1) = tr_wbare(n,1:ngm,i,j)
        tr_w(nx,0:ngm,2) = tr_wvege(n,0:ngm,i,j)
        tr_wsn(nx,1:nlsn,1:2) = tr_wsn_ij(n,1:nlsn, 1:2, i, j)
        ! flux in
        trpr(nx)=(trprec(n,i,j)*bydxyp(j))/dtsrc ! kg/m^2 s
        ! concentration of tracers in atm. water at the surface
        if (qm1.gt.0) then  ! avoid very rare error
          tr_surf(nx) = trm(i,j,1,n)*bydxyp(j)*rhow/qm1 ! kg/m^3
        else
          tr_surf(nx) = 0.
        end if
      enddo
#endif

      do nx=1,ntx
        n=ntix(nx)
C**** set defaults
        trsfac(nx)=0.
        totflux(nx)=0.
        trconstflx(nx)=0.
#ifdef TRACERS_WATER
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
C**** The select is used to distinguish water from gases or particle
! select removed because of OMP compiler bug
!        select case (tr_wd_TYPE(n))
!        case (nWATER)
        if (tr_wd_TYPE(n) .eq. nWATER) then
C**** no fractionation from ground (yet)
C**** trsfac and trconstflx are multiplied by cq*wsh in PBL
          trsfac(nx)=1.
          trconstflx(nx)=gtracer(n,itype,i,j)*QG
!        case (nGAS, nPART)
        elseif (tr_wd_TYPE(n).eq.nGAS .or. tr_wd_TYPE(n).eq.nPART) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) trsfac(nx) = 1.
          !then multiplied by deposition velocity in PBL
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
          do nsrc=1,ntsurfsrc(n)
            totflux(nx) = totflux(nx)+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux(nx)*bydxyp(j)   ! kg/m^2/s
#ifdef TRACERS_WATER
!        end select
        end if
#endif
      end do

#ifdef TRACERS_WATER
c**** water tracers are also flux limited
      do nx=1,ntx
        n=ntix(nx)
C       tr_evap_max(nx) = evap_max * trsoil_rat(nx)
        tr_evap_max(nx) = evap_max * gtracer(n,itype,i,j)
#ifdef TRACERS_DRYDEP
        if(dodrydep(n)) tr_evap_max(nx) = 1.d30
#endif
      end do
#endif

      end subroutine ghy_tracers_set_cell


      subroutine ghy_tracers_save_cell(i,j,ptype,dtsurf,rhosrf
#ifdef INTERACTIVE_WETLANDS_CH4
     & ,ra_temp,ra_sat,ra_gwet
#endif
     & )
!@sum tracers code to be called after the i,j cell is processed
      use model_com, only : itime,qcheck,nisurf
      use ghy_com, only : tearth
      use somtq_com, only : mz
 !     use socpbl, only : dtsurf
      use geom, only : dxyp
      use sle001, only : nsn,fb,fv
#if (defined TRACERS_DUST) && (defined TRACERS_DRYDEP)
      USE trdiag_com,ONLY : rts_save
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      use constant, only : tf
      use tracer_sources, only : n__temp,n__sat,n__gwet
#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: ptype,dtsurf,rhosrf
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      integer :: n1
#endif
      integer n,nx
#ifdef TRACERS_DRYDEP
      real*8 tdryd,tdd,td1,rtsdt,rts
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      real*8, intent(IN) :: ra_temp,ra_sat,ra_gwet
#endif

ccc tracers
#ifdef TRACERS_WATER
      do nx=1,ntg
        n = ntixw(nx)
        tr_wbare(n,1:ngm,i,j) = tr_w(nx,1:ngm,1)
        tr_wvege(n,0:ngm,i,j) = tr_w(nx,0:ngm,2)
        tr_wsn_ij(n,1:nlsn, 1:2, i, j) = tr_wsn(nx,1:nlsn,1:2)
      enddo
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
c     saves precipitation for dust emission calculation at next time step
      pprec(i,j)=prec(i,j)
c     saves evaporation for dust emission calculation at next time step
      pevap(i,j,itype)=aevap
#endif
#ifdef TRACERS_WATER
ccc accumulate tracer evaporation and runoff
      do nx=1,ntg
        n=ntixw(nx)
        trevapor(n,itype,i,j) = trevapor(n,itype,i,j) + atr_evap(nx)
        !trevapor(n,itype,i,j) = trevapor(n,itype,i,j) + aevap  !*rhow
        trunoe(n,i,j) = trunoe(n,i,j) + atr_rnff(nx)
        !trunoe(n,i,j) = trunoe(n,i,j) + (aruns+arunu)  !*rhow
        gtracer(n,itype,i,j) = atr_g(nx)   ! /dtsurf
        trsrfflx(i,j,n)=trsrfflx(i,j,n)+
     &       atr_evap(nx)/dtsurf *dxyp(j)*ptype
      enddo
#endif

      DO nx=1,ntx
        n=ntix(nx)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
C**** technicallly these are ocean emissions, but if fixed datasets
C**** are used, it can happen over land as well.
        select case (trname(n))
        case ('DMS')
          trsrfflx(i,j,n)=trsrfflx(i,j,n)+DMS_flux*dxyp(j)*ptype
          taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &         DMS_flux*dxyp(j)*ptype*dtsurf
          tajls(j,1,jls_isrc(1,n)) = tajls(j,1,jls_isrc(1,n))+
     *         DMS_flux*dxyp(j)*ptype*dtsurf
        case ('seasalt1')
          trsrfflx(i,j,n)=trsrfflx(i,j,n)+ss1_flux*dxyp(j)*ptype
          taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &         ss1_flux*dxyp(j)*ptype*dtsurf
          tajls(j,1,jls_isrc(1,n)) = tajls(j,1,jls_isrc(1,n))+
     *         ss1_flux*dxyp(j)*ptype*dtsurf
        case ('seasalt2')
          trsrfflx(i,j,n)=trsrfflx(i,j,n)+ss2_flux*dxyp(j)*ptype
          taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &         ss2_flux*dxyp(j)*ptype*dtsurf
          tajls(j,1,jls_isrc(1,n)) = tajls(j,1,jls_isrc(1,n))+
     *         ss2_flux*dxyp(j)*ptype*dtsurf
#ifdef TRACERS_AMP
        case ('M_SSA_SS')
       EMIS_SOURCE(i,j,1,5)=EMIS_SOURCE(i,j,1,5)+ss1_flux*dxyp(j)*ptype
       EMIS_SOURCE(i,j,1,6)=EMIS_SOURCE(i,j,1,6)+ss2_flux*dxyp(j)*ptype
        case ('M_DD1_DU')
       EMIS_SOURCE(i,j,1,7)=EMIS_SOURCE(i,j,1,7)
     * +(dust_flux(1)+dust_flux(2))*dxyp(j)*ptype
       EMIS_SOURCE(i,j,1,8)=EMIS_SOURCE(i,j,1,8)
     * +(dust_flux(3)+dust_flux(4))*dxyp(j)*ptype
#endif
        end select
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
ccc dust emission from earth
        SELECT CASE (trname(n))
#ifdef TRACERS_DUST
        CASE ('Clay','Silt1','Silt2','Silt3','Silt4')
          n1=n-n_clay+1
#else
#ifdef TRACERS_MINERALS
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &        'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &        'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &        'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &        'Sil1QuHe','Sil2QuHe','Sil3QuHe')
          n1=n-n_clayilli+1
#else
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          n1=n-n_sil1quhe+1
#endif
#endif
#endif
          trsrfflx(i,j,n)=trsrfflx(i,j,n)+dust_flux(n1)*dxyp(j)*ptype
          taijs(i,j,ijts_source(nDustEmij,n))
     &         =taijs(i,j,ijts_source(nDustEmij,n))+dust_flux(n1)
     &         *dxyp(j)*ptype*dtsurf
          tajls(j,1,jls_source(nDustEmjl,n))
     &         =tajls(j,1,jls_source(nDustEmjl,n))+dust_flux(n1)*dxyp(j)
     &         *ptype*dtsurf
#ifdef TRACERS_DUST
          IF (imDust == 0) THEN
            taijs(i,j,ijts_source(nDustEm2ij,n))
     &           =taijs(i,j,ijts_source(nDustEm2ij,n))+dust_flux2(n1)
     &           *dxyp(j)*ptype*dtsurf
            tajls(j,1,jls_source(nDustEm2jl,n))
     &           =tajls(j,1,jls_source(nDustEm2jl,n))
     &           +dust_flux2(n1)*dxyp(j)*ptype*dtsurf
          END IF
#endif

        END SELECT
#endif

#ifdef TRACERS_DRYDEP

ccc accumulate tracer dry deposition
        if(dodrydep(n)) then
          rts=rhosrf*trs(nx)
          rtsdt=rts*dtsurf                             ! kg*s/m^3
          tdryd=-rtsdt*(dep_vel(n)+gs_vel(n))          ! kg/m2
          tdd = tdryd*dxyp(j)*ptype                    ! kg
          td1 = (trsrfflx(i,j,n)+totflux(nx))*dtsurf   ! kg
          if (trm(i,j,1,n)+td1+tdd.le.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd earth",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,trs(nx),trtop(nx)
            tdd= -max(trm(i,j,1,n)+td1,0d0)
            tdryd= tdd/(dxyp(j)*ptype)
            trsrfflx(i,j,n)= - trm(i,j,1,n)/dtsurf
          else
            trsrfflx(i,j,n)=trsrfflx(i,j,n)+tdd/dtsurf
          end if
#ifdef TRACERS_DUST
          rts_save(i,j)=rts
#endif
! trdrydep downward flux by surface type (kg/m^2)
          trdrydep(n,itype,i,j)=trdrydep(n,itype,i,j) - tdryd
! diagnose turbulent and settling fluxes separately
          taijn(i,j,tij_drydep,n)=taijn(i,j,tij_drydep,n) +
     &         ptype*rtsdt*dep_vel(n)
          taijn(i,j,tij_gsdep ,n)=taijn(i,j,tij_gsdep ,n) +
     &         ptype*rtsdt* gs_vel(n)
          dtr_dd(j,n,1)=dtr_dd(j,n,1)-ptype*rtsdt*dxyp(j)*dep_vel(n)
          dtr_dd(j,n,2)=dtr_dd(j,n,2)-ptype*rtsdt*dxyp(j)* gs_vel(n)
        end if
#endif
      end do

C**** Save surface tracer concentration whether calculated or not
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if (needtrs(n)) then
            nx=nx+1
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)+trs(nx)*ptype
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
            trs_glob(i,j,itype,n)=trs(nx)*ptype
#endif
          else
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *bydxyp(j),0d0)*ptype
          end if
        end if
      end do

ccc not sure about the code below. hopefully that''s what is meant above
#ifdef TRACERS_WATER
      do nx=1,ntg
        n=ntixw(nx)
        taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *       trevapor(n,itype,i,j)*ptype
        taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     *         gtracer(n,itype,i,j)*ptype
        taijn(i,j,tij_soil,n)=taijn(i,j,tij_soil,n) + (
     &       fb*(sum( tr_w(nx,1:ngm,1) ) + sum( tr_wsn(nx,1:nsn(1),1)))+
     &       fv*(sum( tr_w(nx,0:ngm,2) ) + sum( tr_wsn(nx,1:nsn(2),2) ))
     *       )
        tajls(j,1,jls_source(1,n))=tajls(j,1,jls_source(1,n))
     *       +trevapor(n,itype,i,j)*ptype
      enddo
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c ..........
c Accumulates dust events. One diagnostic field for all dust tracers
c ..........
      taijs(i,j,ijts_spec(nDustEv1ij))=taijs(i,j,ijts_spec(nDustEv1ij))
     &     +dust_event1
      tajls(j,1,jls_spec(nDustEv1jl))=tajls(j,1,jls_spec(nDustEv1jl))
     &     +dust_event1
      taijs(i,j,ijts_spec(nDustEv2ij))=taijs(i,j,ijts_spec(nDustEv2ij))
     &     +dust_event2
      tajls(j,1,jls_spec(nDustEv2jl))=tajls(j,1,jls_spec(nDustEv2jl))
     &     +dust_event2
      taijs(i,j,ijts_spec(nDustWthij))=taijs(i,j,ijts_spec(nDustWthij))
     &     +wtrsh*ptype
      tajls(j,1,jls_spec(nDustWthjl))=tajls(j,1,jls_spec(nDustWthjl))
     &     +wtrsh*ptype
#endif

c     ..........
c     save global variables for subdaily diagnostics
c     ..........

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      DO n=1,Ntm_dust
        dust_flux_glob(i,j,n)=dust_flux(n)*ptype
#ifdef TRACERS_DUST
        dust_flux2_glob(i,j,n)=dust_flux2(n)*ptype
#endif
      END DO
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
#ifdef TRACERS_DRYDEP
      DO n=1,Ntm
        IF (dodrydep(n)) THEN
          depo_turb_glob(i,j,itype,n)=ptype*rts*dep_vel(n)
          depo_grav_glob(i,j,itype,n)=ptype*rts*gs_vel(n)
        END IF
      END DO
#endif
#endif

#ifdef INTERACTIVE_WETLANDS_CH4
C**** update running-average of ground temperature:
      call running_average(ra_temp,I,J,dble(nisurf),n__temp)
      call running_average(ra_sat,I,J,dble(nisurf),n__sat)
      call running_average(ra_gwet,I,J,dble(nisurf),n__gwet)
#endif
      end subroutine ghy_tracers_save_cell

      end module ghy_tracers
#endif
c******************   END   TRACERS             ************************


c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************


      module soil_drv
!@sum soil_drv contains variables and routines for the ground
!@+   hydrology driver
!@auth I. Alienov/F. Abramopolous
      use model_com, only : im,jm
      use veg_drv, only : cosday,sinday
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE pbl_drv,ONLY : dust_flux,dust_flux2,z,km,gh,gm,zhat,lmonin,
     &     wsubtke,wsubwd,wsubwm,dust_event1,dust_event2,wtrsh
#endif
      use diag_com, only : idd_ts,idd_tg1,idd_qs
     *     ,idd_qg,idd_swg,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg
     *     ,idd_wg,idd_us,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq
     *     ,idd_eds,idd_dbl,idd_ev
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     *     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf,idd_wtrsh
#endif
#ifdef TRACERS_DUST
     *     ,idd_ws2,idd_ustar,idd_us3,idd_stress,idd_lmon
     *     ,idd_rifl,idd_zpbl1,idd_zpbl2,idd_zpbl3,idd_zpbl4
     *     ,idd_zpbl5,idd_zpbl6,idd_zpbl7,idd_zpbl8
     *     ,idd_uabl1,idd_uabl2,idd_uabl3,idd_uabl4,idd_uabl5
     *     ,idd_uabl6,idd_uabl7,idd_uabl8,idd_vabl1,idd_vabl2
     *     ,idd_vabl3,idd_vabl4,idd_vabl5,idd_vabl6,idd_vabl7
     *     ,idd_vabl8,idd_uvabl1,idd_uvabl2,idd_uvabl3
     *     ,idd_uvabl4,idd_uvabl5,idd_uvabl6,idd_uvabl7
     *     ,idd_uvabl8,idd_tabl1,idd_tabl2,idd_tabl3,idd_tabl4
     *     ,idd_tabl5,idd_tabl6,idd_tabl7,idd_tabl8,idd_qabl1
     *     ,idd_qabl2,idd_qabl3,idd_qabl4,idd_qabl5,idd_qabl6
     *     ,idd_qabl7,idd_qabl8,idd_zhat1,idd_zhat2,idd_zhat3
     *     ,idd_zhat4,idd_zhat5,idd_zhat6,idd_zhat7,idd_e1,idd_e2
     *     ,idd_e3,idd_e4,idd_e5,idd_e6,idd_e7,idd_km1,idd_km2
     *     ,idd_km3,idd_km4,idd_km5,idd_km6,idd_km7,idd_ri1,idd_ri2
     *     ,idd_ri3,idd_ri4,idd_ri5,idd_ri6,idd_ri7
     &     ,idd_emis,idd_emis2,idd_grav,idd_turb
#endif
      implicit none
      private
      save

      public daily_earth
      public ground_e
      public init_gh
      public earth
      public conserv_wtg
      public conserv_htg

      !real*8 cosday,sinday
      !real*8 cosdaym1, sindaym1               !nyk TEMPORARY for jday-1
      real*8 adlmass          ! accumulator for dleafmass in daily_earth

      real*8 spgsn !@var spgsn specific gravity of snow
!@dbparam snow_cover_coef coefficient for topography variance in
!@+       snow cover parameterisation for albedo
      real*8 :: snow_cover_coef = .15d0

      ! Indexes used for adiurn and hdiurn
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER,PARAMETER :: n_idx=28
#else
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idx=114
#else
      INTEGER, PARAMETER :: n_idx = 22
#endif
#endif

      contains

      subroutine earth (ns,moddsf,moddd)
!@sum EARTH calculates surface fluxes of sensible heat,
!@+   evaporation, thermal radiation, and momentum drag over earth.
!@auth I. Alienov/F. Abramopolous
c****
      use constant, only : grav,rgas,lhe,lhs
     *     ,sha,tf,rhow,deltx
      use model_com, only : t,p,q,dtsrc,nisurf,dsig,jdate
     *     ,jday,jhour,nday,itime,jeq,u,v
      use DOMAIN_DECOMP, only : GRID, GET
      use DOMAIN_DECOMP, only : HALO_UPDATE, CHECKSUM, NORTH
      use DOMAIN_DECOMP, only : GLOBALSUM, AM_I_ROOT
      use geom, only : imaxj
      use dynamics, only : pmid,pk,pek,pedn,am
      use rad_com, only : trhr,fsf, cosz1

      !use surf_albedo, only: albvnh   ! added 5/23/03 from RADIATION.f
      !albvnh(9,6,2)=albvnh(sand+8veg,6bands,2hemi) - only need 1st band
      use sle001, only : advnc,evap_limits,
     &    pr,htpr,prs,htprs,w,snowd,tp,fice,
     &    fv,fb,ashg,alhg,
     &    aevap,abetad,
     &    aruns,arunu,aeruns,aerunu,
     &    tbcs,af0dt,af1dt,
     &    qm1,qs,
     &    pres,rho,ts,vsm,ch,srht,trht

      use veg_drv, only: veg_save_cell,veg_set_cell

      use fluxes, only : dth1,dq1,uflux1,vflux1,e0,e1,evapor,prec,eprec
     *     ,runoe,erunoe,gtemp,precss
      use ghy_com, only : snowbv, fearth,
     &     fr_snow_ij,
     *     canopy_temp_ij,snowe,tearth,wearth,aiearth,
     &     evap_max_ij, fr_sat_ij, qg_ij, fr_snow_rad_ij,top_dev_ij

      use vegetation, only :
     &    veg_srht=>srht,veg_pres=>pres,veg_ch=>ch,veg_vsm=>vsm !ia

  !    USE SOCPBL, only :
  !   &     dtsurf               ! zgs,     ! global
  !   &     ,zs1,tgv,tkv,qg_sat,qg_aver,hemi,pole,     ! rest local
  !   &     us,vs,ws,wsm,wsh,tsv,qsrf
  !    use pblcom, only : ipbl,cmgs,chgs,cqgs,qsavg
      use pbl_drv, only : pbl, t_pbl_args
  !   &     , evap_max,fr_sat,uocean,vocean,psurf,trhr0


      use snow_drvm, only : snow_cover_same_as_rad

      use diag_com , only : j_trhdt,j_shdt,j_evhdt,j_evap,j_erun,j_run
     &     ,j_tsrf,j_tg1,j_tg2,areg,jreg,HR_IN_DAY,HR_IN_MONTH,NDIUVAR
     &     ,adiurn,NDIUPT
#ifndef NO_HDIURN
     &     ,hdiurn
#endif
#ifdef TRACERS_ON
      use ghy_tracers, only : ghy_tracers_set_step,ghy_tracers_set_cell,
     &     ghy_tracers_save_cell
#endif
      implicit none

      integer, intent(in) :: ns,moddsf,moddd
      integer i,j,itype,ibv
      real*8 shdt,evhdt,rcdmws,rcdhws
     *     ,cdq,cdm,cdh,elhx,tg,srheat,tg1,ptype,trheat    !,dhgs
     *     ,rhosrf,ma1,tfs,th1,thv1,p1k,psk,ps,pij
     *     ,spring,zs1co,q1

!@var rhosrf0 estimated surface air density
      real*8 rhosrf0


      real*8 qsat
ccc hack for openMP: need temporary array
      real*8 aregij(9,im,jm)
      integer jr
!@var qg rel. humidity at the ground, defined: total_evap = Cq V (qg-qs)
!@var qg_nsat rel. humidity at non-saturated fraction of soil
      real*8 qg, qg_nsat
c****
c**** fearth    soil covered land fraction (1)
c****
c**** snowe     earth snow amount (kg/m**2)
c**** tearth    earth temperature of first layer (c)
c**** wearth    earth water of first layer (kg/m**2)
c**** aiearth   earth ice of first layer (kg/m**2)
c****
c**** wbare  1-6 water of bare soil layer 1-6 (m)
c**** wvege   0  water of vegetation canopy (m)
c****        1-6 water of vegetated soil layer 1-6 (m)
c**** htbare  0  bare soil layer 0 is unused
c****        1-6 heat content of bare soil layer 1-6 (j m-2)
c**** htvege  0  heat content of vegetation canopy (j m-2)
c****        1-6 heat content of vegetated soil layer 1-6 (j m-2)
c**** snowbv  1  snow depth over bare soil (m)
c****         2  snow depth over vegetated soil (m)
c****

c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 qg_sat

C**** Work array for regional diagnostic accumulation
      real*8, dimension(size(AREG,1),9) :: areg_sum
      real*8, DIMENSION(
     &        size(AREG,1),grid%J_STRT_HALO:grid%J_STOP_HALO,9 )
     &        :: AREG_PART

      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUVAR, NDIUPT) :: adiurn_part
#ifndef NO_HDIURN
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUVAR, NDIUPT) :: hdiurn_part
#endif
      REAL*8, DIMENSION(N_IDX,grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUPT) :: adiurn_temp
#ifndef NO_HDIURN
      REAL*8, DIMENSION(N_IDX,grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUPT) :: hdiurn_temp
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      REAL*8, DIMENSION(N_IDX, NDIUPT) :: ADIURNSUM
#else
      REAL*8, DIMENSION(N_IDX, NDIUPT) :: ADIURNSUM, HDIURNSUM
#endif
      INTEGER :: ih, ihm, ii, ivar, kr
      INTEGER :: idx(n_idx)



C****   define local grid
      integer J_0, J_1, J_0H, J_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      ! dtsurf=dtsrc/nisurf

      zs1co=.5*dsig(1)*rgas/grav

      spring=-1.
      if((jday.ge.32).and.(jday.le.212)) spring=1.
c****
c**** outside loop over time steps, executed nisurf times every hour
c****
ccc set i,j - independent stuff for tracers
#ifdef TRACERS_ON
      call ghy_tracers_set_step
#endif

      aregij = 0.
c****
c**** outside loop over j and i, executed once for each grid point
c****
C**** halo update u and v for distributed parallelization
       call halo_update(grid, U, from=NORTH)
       call halo_update(grid, V, from=NORTH)

       adiurn_part = 0
#ifndef NO_HDIURN
       hdiurn_part = 0
#endif

       idx =
     &     (/ idd_ts,  idd_tg1, idd_qs,  idd_qg,  idd_swg,
     &        idd_lwg, idd_sh,  idd_lh,  idd_hz0, idd_ug,
     &        idd_vg,  idd_wg,  idd_us,  idd_vs,  idd_ws,
     &        idd_cia, idd_cm,  idd_ch,  idd_cq,  idd_eds,
     &        idd_dbl, idd_ev
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &      ,idd_wsgcm,idd_wspdf,idd_wtke,idd_wd,idd_wm,idd_wtrsh
#endif
#ifdef TRACERS_DUST
     &      ,idd_emis,idd_emis2,idd_turb,idd_grav
     &      ,idd_ws2,idd_ustar,idd_us3,idd_stress,idd_lmon,idd_rifl
     &      ,idd_zpbl1,idd_zpbl2,idd_zpbl3,idd_zpbl4,idd_zpbl5,idd_zpbl6
     &      ,idd_zpbl7,idd_zpbl8
     &      ,idd_uabl1,idd_uabl2,idd_uabl3,idd_uabl4,idd_uabl5,idd_uabl6
     &      ,idd_uabl7,idd_uabl8
     &      ,idd_vabl1,idd_vabl2,idd_vabl3,idd_vabl4,idd_vabl5,idd_vabl6
     &      ,idd_vabl7,idd_vabl8
     &      ,idd_uvabl1,idd_uvabl2,idd_uvabl3,idd_uvabl4,idd_uvabl5
     &      ,idd_uvabl6,idd_uvabl7,idd_uvabl8
     &      ,idd_tabl1,idd_tabl2,idd_tabl3,idd_tabl4,idd_tabl5,idd_tabl6
     &      ,idd_tabl7,idd_tabl8
     &      ,idd_qabl1,idd_qabl2,idd_qabl3,idd_qabl4,idd_qabl5,idd_qabl6
     &      ,idd_qabl7,idd_qabl8
     &      ,idd_zhat1,idd_zhat2,idd_zhat3,idd_zhat4,idd_zhat5,idd_zhat6
     &      ,idd_zhat7
     &      ,idd_e1,idd_e2,idd_e3,idd_e4,idd_e5,idd_e6,idd_e7
     &      ,idd_km1,idd_km2,idd_km3,idd_km4,idd_km5,idd_km6,idd_km7
     &      ,idd_ri1,idd_ri2,idd_ri3,idd_ri4,idd_ri5,idd_ri6,idd_ri7
#endif
     &      /)

!$OMP  PARALLEL DO PRIVATE
!$OMP*  (ELHX,EVHDT, CDM,CDH,CDQ,
!$OMP*   I,ITYPE,ibv, J, MA1,PIJ,PSK,PS,P1K,PTYPE, QG,
!$OMP*   QG_NSAT, RHOSRF,RHOSRF0,RCDMWS,RCDHWS,SRHEAT,SHDT,
!$OMP*   TRHEAT, TH1,TFS,THV1,TG1,TG,q1,pbl_args,qg_sat
!$OMP*   )
!$OMP*   SHARED(ns,moddsf,moddd)
!$OMP*   SCHEDULE(DYNAMIC,2)

      loop_j: do j=J_0,J_1
  !    hemi=1.
  !    if(j.le.jm/2) hemi=-1.

c**** conditions at the south/north pole
  !    pole= ( j.eq.1 .or. j.eq.jm )

      loop_i: do i=1,imaxj(j)

      pbl_args%dtsurf=dtsrc/nisurf
      pbl_args%hemi=1.
      if(j.le.jm/2) pbl_args%hemi=-1.
      pbl_args%pole= ( j.eq.1 .or. j.eq.jm )

c****
c**** determine surface conditions
c****
      ptype=fearth(i,j)
      pij=p(i,j)
      ps=pedn(1,i,j)
      psk=pek(1,i,j)
      p1k=pk(1,i,j)
      th1=t(i,j,1)
      q1=q(i,j,1)
      thv1=th1*(1.+q1*deltx)
      pbl_args%tkv=thv1*psk
  !    tkv=thv1*psk

      tfs=tf*ptype
      ma1=am(1,i,j)
      qm1=q1*ma1
c     rhosrf=100.*ps/(rgas*tsv)
c     rhosrf=100.*ps/(rgas*tkv)
c****
c**** earth
c****
      if (ptype.le.0.) then
  !      ipbl(i,j,4)=0
        cycle loop_i
      endif
      itype=4

      pr=prec(i,j)/(dtsrc*rhow)
C**** This variable was originally assoicated with super-saturated
C**** large-scale precip, but has not worked for many moons.
C**** If you want to reinstate it, uncomment this calculation.
c      prs=precss(i,j)/(dtsrc*rhow)
      prs=0.
      htpr=eprec(i,j)/dtsrc

ccc tracers variables

      tg1 = tearth(i,j)
      srheat=fsf(itype,i,j)*cosz1(i,j)
c****
c**** boundary layer interaction
c****
      pbl_args%zs1=zs1co*pbl_args%tkv*pij/pmid(1,i,j)
  !    zs1=zs1co*tkv*pij/pmid(1,i,j)

c**** loop over ground time steps
      tg=tg1+tf
      elhx=lhe
      if(tg1.lt.0.)  elhx=lhs
      pbl_args%qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step
  !    qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step

      qg = qg_ij(i,j)
      ! if ( qg > 999.d0 ) qg = qg_sat
      pbl_args%qg_aver = qg
  !    qg_aver = qg

      pbl_args%tgv=tg*(1.+qg*deltx)
  !    tgv=tg*(1.+qg*deltx)

      pbl_args%psurf=ps
  !    psurf=ps

      pbl_args%trhr0 = TRHR(0,I,J)
  !    trhr0 = TRHR(0,I,J)

      rhosrf0=100.*ps/(rgas*pbl_args%tgv) ! estimated surface density
C**** Obviously there are no ocean currents for earth points, but
C**** variables set for consistency with surfce
      pbl_args%uocean=0 ; pbl_args%vocean=0
  !    uocean=0 ; vocean=0

c***********************************************************************
c****
ccc actually PBL needs evap (kg/m^2*s) / rho_air
      pbl_args%evap_max = evap_max_ij(i,j) * 1000.d0 / rhosrf0
      pbl_args%fr_sat = fr_sat_ij(i,j)
  !    evap_max = evap_max_ij(i,j) * 1000.d0 / rhosrf0
  !    fr_sat = fr_sat_ij(i,j)

c**** call tracers stuff
#ifdef TRACERS_ON
      call ghy_tracers_set_cell(i,j,qg,pbl_args%evap_max)
#endif
      call pbl(i,j,itype,ptype,pbl_args)
c****
      cdm = pbl_args%cm ! cmgs(i,j,itype)
      cdh = pbl_args%ch ! chgs(i,j,itype)
      cdq = pbl_args%cq ! cqgs(i,j,itype)
c***********************************************************************
c**** calculate qs
      qs=pbl_args%qsrf
      ts=pbl_args%tsv/(1.+qs*deltx)
c**** calculate rhosrf*cdm*ws
      rhosrf=100.*ps/(rgas*pbl_args%tsv)
      rcdmws=cdm*pbl_args%wsm*rhosrf
      rcdhws=cdh*pbl_args%wsh*rhosrf
      trheat=trhr(0,i,j)
c***********************************************************************
c****
c  define extra variables to be passed in surfc:
      pres  =ps
      veg_pres = ps
      rho   =rhosrf
      vsm   =pbl_args%ws
      veg_vsm = pbl_args%ws
      ch    =cdh
      veg_ch = cdh
      srht  =srheat
      veg_srht = srheat
      trht  =trheat
c***********************************************************************
c****
c**** calculate ground fluxes
c     call qsbal
!!! insert htprs here ???


      call ghinij (i,j)
      call veg_set_cell(i,j)
      call advnc
      call evap_limits( .false., evap_max_ij(i,j), fr_sat_ij(i,j) )

      call veg_save_cell(i,j)
      call ghy_save_cell(i,j)

      tg1=tbcs

c**** computing ground humidity to be used on next time step
      !qg_ij(i,j) = qs  !!! - this seemed to work ok
      !! trying more precise value for qg :  qsat(tg1+tf,elhx,ps)
      qg_sat = qsat(tg1+tf,elhx,ps) ! saturated soil
      qg_nsat = qs              ! non-sat soil, no evap
      if ( rcdhws > 1.d-30 ) then   ! correction to non-sat, due to evap
        qg_nsat = qg_nsat + evap_max_ij(i,j)/(0.001*rcdhws)
      endif
      qg_nsat = min( qg_nsat, qg_sat )
      qg_ij(i,j) = fr_sat_ij(i,j) * qg_sat
     &     + (1.d0 -fr_sat_ij(i,j)) * qg_nsat

ccc Save canopy temperature.
ccc canopy_temp_ij is not used so far... do we need it?
      canopy_temp_ij(i,j) = tp(0,2)  !nyk

c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      if ( snow_cover_same_as_rad == 0 ) then
        ! recompute snow fraction using different formula
        do ibv=1,2
          call snow_cover(fr_snow_rad_ij(ibv,i,j),
     &         snowbv(ibv,i,j), top_dev_ij(i,j) )
          fr_snow_rad_ij(ibv,i,j) = min (
     &         fr_snow_rad_ij(ibv,i,j), fr_snow_ij(ibv, i, j) )
        enddo
      else
        ! snow fraction same as in snow model
        fr_snow_rad_ij(:,i,j) = fr_snow_ij(:, i, j)
      endif

c**** snowe used in RADIATION
      snowe(i,j)=1000.*(snowd(1)*fb+snowd(2)*fv)
c**** tearth used only internaly in GHY_DRV
      tearth(i,j)=tg1
c**** wearth+aiearth are used in radiation only
      wearth(i,j)=1000.*( fb*w(1,1)*(1.-fice(1,1)) +
     &     fv*(w(1,2)*(1.-fice(1,2))+w(0,2)*(1.-fice(0,2))) )
      aiearth(i,j)=1000.*( fb*w(1,1)*fice(1,1) +
     &     fv*(w(1,2)*fice(1,2)+w(0,2)*fice(0,2)) )
      gtemp(1,4,i,j)=tearth(i,j)
c**** calculate fluxes using implicit time step for non-ocean points
      uflux1(i,j)=uflux1(i,j)+ptype*rcdmws*(pbl_args%us) !-uocean)
      vflux1(i,j)=vflux1(i,j)+ptype*rcdmws*(pbl_args%vs) !-vocean)
c**** accumulate surface fluxes and prognostic and diagnostic quantities
      evapor(i,j,4)=evapor(i,j,4)+aevap
      shdt=-ashg
      dth1(i,j)=dth1(i,j)-shdt*ptype/(sha*ma1*p1k)
      dq1(i,j) =dq1(i,j)+aevap*ptype/ma1
  !    qsavg(i,j)=qsavg(i,j)+qs*ptype
c**** save runoff for addition to lake mass/energy resevoirs
      runoe (i,j)=runoe (i,j)+ aruns+ arunu
      erunoe(i,j)=erunoe(i,j)+aeruns+aerunu
c****
      e0(i,j,4)=e0(i,j,4)+af0dt
      e1(i,j,4)=e1(i,j,4)+af1dt

      call ghy_diag( i,j,ns,moddsf,moddd
     &     ,rcdmws,cdm,cdh,cdq,qg
     &     ,aregij, pbl_args, pbl_args%dtsurf
     &     ,adiurn_part
#ifndef NO_HDIURN
     &     ,hdiurn_part
#endif
     &     ,idx)

c**** update tracers
#ifdef TRACERS_ON
      call ghy_tracers_save_cell(i,j,ptype,pbl_args%dtsurf,rhosrf
#ifdef INTERACTIVE_WETLANDS_CH4
     & ,tg1,ts-tf,1.d2*abetad/real(nisurf)
#endif
     & )
#endif

      end do loop_i
      end do loop_j
!$OMP  END PARALLEL DO

! Accumulate contributions to ADIURN and HDIURN
      ih=1+jhour
      ihm = ih+(jdate-1)*24


      DO kr = 1, ndiupt
        DO ii = 1, N_IDX
          ivar = idx(ii)
          ADIURN_temp(ii,:,kr)=ADIURN_part(:,ivar,kr)
#ifndef NO_HDIURN
          HDIURN_temp(ii,:,kr)=HDIURN_part(:,ivar,kr)
#endif
        END DO
      END DO
      CALL GLOBALSUM(grid, ADIURN_temp(1:N_IDX,:,1:ndiupt),
     &    ADIURNSUM(1:N_IDX,1:ndiupt))
#ifndef NO_HDIURN
      CALL GLOBALSUM(grid, HDIURN_temp(1:N_IDX,:,1:ndiupt),
     &    HDIURNSUM(1:N_IDX,1:ndiupt))
#endif
      DO kr = 1, ndiupt
        DO ii = 1, N_IDX
          ivar = idx(ii)
          IF (AM_I_ROOT()) THEN
            ADIURN(ih,ivar,kr)=ADIURN(ih,ivar,kr) + ADIURNSUM(ii,kr)
#ifndef NO_HDIURN
            HDIURN(ihm,ivar,kr)=HDIURN(ihm,ivar,kr) + HDIURNSUM(ii,kr)
#endif
          END IF
        END DO
      END DO

C***Initialize work array
      areg_part(:,:,1:9) = 0.

      DO 825 J=J_0,J_1
      DO 825 I=1,IMAXJ(J)
         IF(FEARTH(I,J).LE.0.0)  GO TO 825
         JR=JREG(I,J)
         areg_part(jr,j,1 )=areg_part(jr,j,1 )+AREGIJ(1,I,J)
         areg_part(jr,j,2 )=areg_part(jr,j,2 )+AREGIJ(2,I,J)
         areg_part(jr,j,3 )=areg_part(jr,j,3 )+AREGIJ(3,I,J)
         areg_part(jr,j,4 )=areg_part(jr,j,4 )+AREGIJ(4,I,J)
         areg_part(jr,j,5 )=areg_part(jr,j,5 )+AREGIJ(5,I,J)
         areg_part(jr,j,6 )=areg_part(jr,j,6 )+AREGIJ(6,I,J)
         if( moddsf == 0 ) then
           areg_part(jr,j,7)=areg_part(jr,j,7 )+AREGIJ(7,I,J)
           areg_part(jr,j,8)=areg_part(jr,j,8 )+AREGIJ(8,I,J)
           areg_part(jr,j,9)=areg_part(jr,j,9 )+AREGIJ(9,I,J)
         end if
  825 CONTINUE


      call globalsum(grid,areg_part(1:size(areg,1),:,1:9),
     &    areg_sum(1:size(areg,1),1:9), all=.true.)
      areg(1:size(areg,1),j_trhdt)=areg(1:size(areg,1),j_trhdt)
     &    + areg_sum(1:size(areg,1),1)

      areg(1:size(areg,1),j_shdt)=areg(1:size(areg,1),j_shdt)
     &    + areg_sum(1:size(areg,1),2)

      areg(1:size(areg,1),j_evhdt)=areg(1:size(areg,1),j_evhdt)
     &    + areg_sum(1:size(areg,1),3)

      areg(1:size(areg,1),j_evap)=areg(1:size(areg,1),j_evap)
     &    + areg_sum(1:size(areg,1),4)

      areg(1:size(areg,1),j_erun)=areg(1:size(areg,1),j_erun)
     &    + areg_sum(1:size(areg,1),5)

      areg(1:size(areg,1),j_run)=areg(1:size(areg,1),j_run)
     &    + areg_sum(1:size(areg,1),6)

      areg(1:size(areg,1),j_tsrf)=areg(1:size(areg,1),j_tsrf)
     &    + areg_sum(1:size(areg,1),7)

      areg(1:size(areg,1),j_tg1)=areg(1:size(areg,1),j_tg1)
     &    + areg_sum(1:size(areg,1),8)

      areg(1:size(areg,1),j_tg2)=areg(1:size(areg,1),j_tg2)
     &    + areg_sum(1:size(areg,1),9)


C
      return
      end subroutine earth

c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine ghy_diag( i,j,ns,moddsf,moddd
     &     ,rcdmws,cdm,cdh,cdq,qg
     &     ,aregij, pbl_args, dtsurf
     &     ,adiurn_part
#ifndef NO_HDIURN
     &     ,hdiurn_part
#endif
     &     ,idx)

      use diag_com , only : aij=>aij_loc
     *     ,tsfrez=>tsfrez_loc,tdiurn,aj=>aj_loc,areg,jreg
     *     ,ij_rune, ij_arunu, ij_pevap, ij_shdt, ij_beta, ij_trnfp0
     *     ,ij_srtr, ij_neth, ij_ws, ij_ts, ij_us, ij_vs, ij_taus
     *     ,ij_tauus, ij_tauvs, ij_qs, ij_tg1, ij_evap, j_trhdt, j_shdt
     *     ,j_evhdt,j_evap,j_erun,j_run,j_tsrf,j_type,j_tg1,j_tg2,ij_g05
     *     ,ij_g06,ij_g11,ij_g12,ij_g13,ij_g14,ij_g15,ij_g16,ij_g17
     *     ,ij_gpp,ij_pblht,ij_g18,ij_g19,ij_g20,ij_g21,ij_g22,ij_g23
     *     ,ij_g24,ij_g25,ij_g26,ij_g27,ijdd,idd_ts,idd_tg1,idd_qs
     *     ,idd_qg,idd_swg,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg
     *     ,idd_wg,idd_us,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq
     *     ,idd_eds,idd_dbl,idd_ev,tf_day1,tf_last,ndiupt
     *     ,HR_IN_DAY,HR_IN_MONTH,NDIUVAR
     &     ,ij_aflmlt,ij_aeruns,ij_aerunu
     &     ,ij_htsoil,ij_htsnow,ij_aintrcp,ij_trsdn,ij_trsup,adiurn_dust
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
#endif
#ifdef TRACERS_DUST
      USE PBLCOM,ONLY : eabl,uabl,vabl,tabl,qabl
#endif
      use constant, only : tf
#ifdef TRACERS_DUST
     &     ,grav
#endif
      use model_com, only : dtsrc,nisurf,jdate
     *     ,jday,jhour,nday,itime,jeq,modrd,itearth
      use DOMAIN_DECOMP, only : grid
      use geom, only : dxyp
      use rad_com, only : trhr,fsf, cosz1

      use sle001, only :
     &     tp
     &    ,fv,fb,atrg,ashg,alhg
     &    ,abetad,abetav,abetat
     &    ,abetap,abetab,abeta
     &    ,acna,acnc,agpp
     &    ,aevap,aevapw,aevapd,aevapb
     &    ,aruns,arunu,aeruns,aerunu,aflmlt,aintercep
     &    ,aepc,aepb,aepp,zw,tbcs
     &    ,qs,ts,ngr=>n,ht,hsn,fr_snow,nsn

      use ghy_com, only : gdeep, fearth

 !     USE SOCPBL, only : dtsurf         ! zgs,     ! global
 !    &     ,us,vs,ws,wsm,psi,dbl    ! ,edvisc=>kms
 !    &     ,khs,ug,vg,wg
 !     use pbl_drv, only : uocean,vocean
      use pbl_drv, only : t_pbl_args
#if (defined TRACERS_DUST) && (defined TRACERS_DRYDEP)
     &     ,dep_vel,gs_vel
#endif

#ifdef TRACERS_DUST
      USE tracer_com,ONLY : Ntm_dust,n_clay
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#endif
#if (defined TRACERS_DUST) && (defined TRACERS_DRYDEP)
      USE trdiag_com,ONLY : rts_save
#endif

      implicit none
      integer, intent(in) :: i,j,ns,moddsf,moddd
      real*8, intent(in) :: rcdmws,cdm,cdh,cdq,qg
      real*8, intent(out) :: aregij(:,:,:)
      type (t_pbl_args) :: pbl_args
      real*8, intent(in) :: dtsurf
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUVAR, NDIUPT) :: adiurn_part
#ifndef NO_HDIURN
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUVAR, NDIUPT) :: hdiurn_part
#endif
      INTEGER, INTENT(IN) :: idx(:)

      real*8 timez
      real*8 trhdt,tg2av,wtr2av,ace2av,tg1,shdt,ptype,srheat,srhdt
      real*8 warmer,spring,trheat,evhdt
      integer, parameter :: itype=4
      integer kr,ih,ihm,jr
#ifdef TRACERS_DUST
      INTEGER :: n,n1
#endif

      REAL*8 :: tmp(NDIUVAR)

ccc the following values are returned by PBL
      real*8 us,vs,ws,wsm,psi,dbl,khs,ug,vg,wg
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,wsgcm,wspdf
#endif
      us = pbl_args%us
      vs = pbl_args%vs
      ws = pbl_args%ws
      wsm = pbl_args%wsm
      psi = pbl_args%psi
      dbl = pbl_args%dbl
      khs = pbl_args%khs
      ug = pbl_args%ug
      vg = pbl_args%vg
      wg = pbl_args%wg
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      wsgcm=pbl_args%wsgcm
      wspdf=pbl_args%wspdf
#endif

      timez=jday+(mod(itime,nday)+(ns-1.)/nisurf)/nday ! -1 ??
      if(jday.le.31) timez=timez+365.


      spring=-1.
      if((jday.ge.32).and.(jday.le.212)) spring=1.

      if(j.lt.jeq)  then
         warmer=-spring
       else
         warmer=spring
      end if
      ih=1+jhour
      ihm = ih+(jdate-1)*24

      ptype=fearth(i,j)
      jr=jreg(i,j)

      srheat=fsf(itype,i,j)*cosz1(i,j)
      srhdt=srheat*dtsurf
      trheat=trhr(0,i,j)

      tg1=tbcs
      shdt=-ashg
      evhdt=-alhg


      aij(i,j,ij_g18)=aij(i,j,ij_g18)+aevapb
      aij(i,j,ij_g19)=aij(i,j,ij_g19)+aevapd
      aij(i,j,ij_g20)=aij(i,j,ij_g20)+aevapw
      aij(i,j,ij_g05)=aij(i,j,ij_g05)+abetab/nisurf
      aij(i,j,ij_g06)=aij(i,j,ij_g06)+abetap/nisurf
      aij(i,j,ij_g11)=aij(i,j,ij_g11)+abeta/nisurf
      aij(i,j,ij_g12)=aij(i,j,ij_g12)+acna/nisurf
      aij(i,j,ij_g13)=aij(i,j,ij_g13)+acnc/nisurf
      aij(i,j,ij_gpp)=aij(i,j,ij_gpp)+agpp
      aij(i,j,ij_g26)=aij(i,j,ij_g26)+abetav/nisurf
      aij(i,j,ij_g27)=aij(i,j,ij_g27)+abetat/nisurf
      aij(i,j,ij_g14)=aij(i,j,ij_g14)+aepp
      if (moddsf.eq.0) then
        aij(i,j,ij_g15)=aij(i,j,ij_g15)+tp(1,1)
        aij(i,j,ij_g16)=aij(i,j,ij_g16)+tp(2,1)
        aij(i,j,ij_g17)=aij(i,j,ij_g17)+tp(6,1)
        aij(i,j,ij_g21)=aij(i,j,ij_g21)+tp(0,2)
        aij(i,j,ij_g22)=aij(i,j,ij_g22)+tp(1,2)
        aij(i,j,ij_g23)=aij(i,j,ij_g23)+tp(2,2)
        aij(i,j,ij_g24)=aij(i,j,ij_g24)+tp(6,2)
        aij(i,j,ij_g25)=aij(i,j,ij_g25)+fb*zw(1)+fv*zw(2)
      end if
ccc accumulate total heat storage
      if (moddsf.eq.0) then
        aij(i,j,ij_htsoil)=aij(i,j,ij_htsoil) +
     &       fb*sum(ht(1:ngr,1)) + fv*sum(ht(0:ngr,2))
        aij(i,j,ij_htsnow)=aij(i,j,ij_htsnow)
     &       + fb*fr_snow(1)*sum(hsn(1:nsn(1),1))
     &       + fv*fr_snow(2)*sum(hsn(1:nsn(2),2))
      endif
      trhdt=trheat*dtsurf-atrg
c           for radiation find composite values over earth
c           for diagnostic purposes also compute gdeep 1 2 3
      call retp2 (tg2av,wtr2av,ace2av)
      gdeep(i,j,1)=tg2av
      gdeep(i,j,2)=wtr2av
      gdeep(i,j,3)=ace2av

      aij(i,j,ij_rune)=aij(i,j,ij_rune)+aruns
      aij(i,j,ij_arunu)=aij(i,j,ij_arunu)+arunu
      aij(i,j,ij_aeruns)=aij(i,j,ij_aeruns)+aeruns
      aij(i,j,ij_aerunu)=aij(i,j,ij_aerunu)+aerunu
      aij(i,j,ij_pevap)=aij(i,j,ij_pevap)+(aepc+aepb)
      aij(i,j,ij_aflmlt)=aij(i,j,ij_aflmlt)+aflmlt
      aij(i,j,ij_aintrcp)= aij(i,j,ij_aintrcp)+aintercep

      if ( warmer >= 0 ) then
        if(ts.lt.tf) tsfrez(i,j,tf_day1)=timez
        tsfrez(i,j,tf_last)=timez
      else
        if ( tsfrez(i,j,tf_last)+.03 >= timez .and. ts >= tf )
     $       tsfrez(i,j,tf_last)=timez
      endif

      if(tg1.lt.tdiurn(i,j,1)) tdiurn(i,j,1)=tg1
      if(tg1.gt.tdiurn(i,j,2)) tdiurn(i,j,2)=tg1
      if(ts.lt.tdiurn(i,j,3)) tdiurn(i,j,3)=ts
      if(ts.gt.tdiurn(i,j,4)) tdiurn(i,j,4)=ts

c**** quantities accumulated for regions in diagj
cddd      areg(jr,j_trhdt)=areg(jr,j_trhdt)+trhdt*ptype*dxyp(j)
cddd      areg(jr,j_shdt )=areg(jr,j_shdt )+shdt*ptype*dxyp(j)
cddd      areg(jr,j_evhdt)=areg(jr,j_evhdt)+evhdt*ptype*dxyp(j)
cddd      areg(jr,j_evap )=areg(jr,j_evap )+aevap*ptype*dxyp(j)
cddd      areg(jr,j_erun)=areg(jr,j_erun)+(aeruns+aerunu)*ptype*dxyp(j)
cddd      areg(jr,j_run )=areg(jr,j_run )+(aruns+arunu)*ptype*dxyp(j)
cddd      if ( moddsf == 0 ) then
cddd        areg(jr,j_tsrf )=areg(jr,j_tsrf )+(ts-tf)*ptype*dxyp(j)
cddd        areg(jr,j_tg1 ) =areg(jr,j_tg1 ) + tg1   *ptype*dxyp(j)
cddd        areg(jr,j_tg2 ) =areg(jr,j_tg2 ) + tg2av *ptype*dxyp(j)
cddd      end if

c!!! do something with regional diag !!!
      AREGIJ(1,I,J)  = trhdt*ptype*dxyp(j)
      AREGIJ(2,I,J)  = shdt*ptype*dxyp(j)
      AREGIJ(3,I,J)  = evhdt*ptype*dxyp(j)
      AREGIJ(4,I,J)  = aevap*ptype*dxyp(j)
      AREGIJ(5,I,J)  = (aeruns+aerunu)*ptype*dxyp(j)
      AREGIJ(6,I,J)  = (aruns+arunu)*ptype*dxyp(j)
      if ( moddsf == 0 ) THEN
        AREGIJ(7,I,J)  = (ts-tf)*ptype*dxyp(j)
        AREGIJ(8,I,J)  = tg1    *ptype*dxyp(j)
        AREGIJ(9,I,J)  = tg2av  *ptype*dxyp(j)
      end if
c**** quantities accumulated for latitude-longitude maps in diagij
      aij(i,j,ij_shdt)=aij(i,j,ij_shdt)+shdt*ptype
      aij(i,j,ij_beta)=aij(i,j,ij_beta)+abetad/nisurf
      IF (MODDSF.EQ.0) THEN
        AIJ(I,J,IJ_TRSDN)=AIJ(I,J,IJ_TRSDN)+TRHR(0,I,J)*PTYPE
        AIJ(I,J,IJ_TRSUP)=AIJ(I,J,IJ_TRSUP)+(TRHR(0,I,J)-TRHDT/DTSURF)
     *       *PTYPE
      END IF
      if(modrd.eq.0)aij(i,j,ij_trnfp0)=aij(i,j,ij_trnfp0)+trhdt*ptype
     *     /dtsrc
      aij(i,j,ij_srtr)=aij(i,j,ij_srtr)+(srhdt+trhdt)*ptype
      aij(i,j,ij_neth)=aij(i,j,ij_neth)+(srhdt+trhdt+shdt+evhdt)*ptype
      aij(i,j,ij_evap)=aij(i,j,ij_evap)+aevap*ptype
      if ( moddsf == 0 ) then
        aij(i,j,ij_ws)=aij(i,j,ij_ws)+ws*ptype
        aij(i,j,ij_ts)=aij(i,j,ij_ts)+(ts-tf)*ptype
        aij(i,j,ij_us)=aij(i,j,ij_us)+us*ptype
        aij(i,j,ij_vs)=aij(i,j,ij_vs)+vs*ptype
        aij(i,j,ij_taus)=aij(i,j,ij_taus)+rcdmws*wsm*ptype
        aij(i,j,ij_tauus)=aij(i,j,ij_tauus)+rcdmws*(us)*ptype !-uocean
        aij(i,j,ij_tauvs)=aij(i,j,ij_tauvs)+rcdmws*(vs)*ptype !-vocean
        aij(i,j,ij_qs)=aij(i,j,ij_qs)+qs*ptype
        aij(i,j,ij_tg1)=aij(i,j,ij_tg1)+tg1*ptype
        aij(i,j,ij_pblht)=aij(i,j,ij_pblht)+dbl*ptype
chyd       aij(i,j,ij_arunu)=aij(i,j,ij_arunu)
chyd      *  +   (40.6*psoil+.72*(2.*(tss-tfs)-(qsatss-qss)*lhe/sha))

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        aij(i,j,ij_wsgcm)=aij(i,j,ij_wsgcm)+wsgcm*ptype
        aij(i,j,ij_wspdf)=aij(i,j,ij_wspdf)+wspdf*ptype
        aij(i,j,ij_wdry)=aij(i,j,ij_wdry)+wsubwd*ptype
        aij(i,j,ij_wtke)=aij(i,j,ij_wtke)+wsubtke*ptype
        aij(i,j,ij_wmoist)=aij(i,j,ij_wmoist)+wsubwm*ptype
#endif

      endif

c**** quantities accumulated hourly for diagDD
      if ( moddd == 0 ) then
        do kr=1,ndiupt
          if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then

            tmp(idd_ts)=+ts*ptype
            tmp(idd_tg1)=+(tg1+tf)*ptype
            tmp(idd_qs)=+qs*ptype
            tmp(idd_qg)=+qg*ptype
            tmp(idd_swg)=+srhdt*ptype
            tmp(idd_lwg)=+trhdt*ptype
            tmp(idd_sh)=+shdt*ptype
            tmp(idd_lh)=+evhdt*ptype
            tmp(idd_hz0)=
     *           +(srhdt+trhdt+shdt+evhdt)*ptype
            tmp(idd_ug)=+ug*ptype
            tmp(idd_vg)=+vg*ptype
            tmp(idd_wg)=+wg*ptype
            tmp(idd_us)=+us*ptype
            tmp(idd_vs)=+vs*ptype
            tmp(idd_ws)=+ws*ptype
            tmp(idd_cia)=+psi*ptype
            tmp(idd_cm)=+cdm*ptype
            tmp(idd_ch)=+cdh*ptype
            tmp(idd_cq)=+cdq*ptype
            tmp(idd_eds)=+khs*ptype
            tmp(idd_dbl)=+dbl*ptype
            tmp(idd_ev)=+aevap*ptype
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
            IF (adiurn_dust == 1) THEN
              tmp(idd_wsgcm)=+wsgcm*ptype
              tmp(idd_wspdf)=+wspdf*ptype
              tmp(idd_wtke)=+wsubtke*ptype
              tmp(idd_wd)=+wsubwd*ptype
              tmp(idd_wm)=+wsubwm*ptype
              tmp(idd_wtrsh)=+wtrsh*ptype
            END IF
#endif
#ifdef TRACERS_DUST
            IF (adiurn_dust == 1) THEN

              tmp(idd_emis)=0.D0
              tmp(idd_emis2)=0.D0
              tmp(idd_turb)=0.D0
              tmp(idd_grav)=0.D0
              DO n=1,Ntm_dust
                n1=n_clay+n-1
                tmp(idd_emis)=tmp(idd_emis)+dust_flux(n)*ptype
                tmp(idd_emis2)=tmp(idd_emis2)+dust_flux2(n)*ptype
#ifdef TRACERS_DRYDEP
                IF (dodrydep(n1)) THEN
                  tmp(idd_turb)=tmp(idd_turb)+ptype*rts_save(i,j)
     &                 *dep_vel(n1)
                  tmp(idd_grav)=tmp(idd_grav)+ptype*rts_save(i,j)
     &                 *gs_vel(n1)
                END IF
#endif
              END DO

              tmp(idd_ws2)=+ws*ws*ptype
              tmp(idd_ustar)=+pbl_args%ustar*ptype
              tmp(idd_us3)=+ptype*pbl_args%ustar*pbl_args%ustar
     &             *pbl_args%ustar
              tmp(idd_stress)=+rcdmws*wsm*ptype
              tmp(idd_lmon)=+lmonin*ptype
              tmp(idd_rifl)=
     &             +ptype*grav*(ts-(tg1+tf))*pbl_args%zgs
     &             /(ws*ws*(tg1+tf))

              tmp(idd_zpbl1)=+ptype*z(1)
              tmp(idd_zpbl2)=+ptype*z(2)
              tmp(idd_zpbl3)=+ptype*z(3)
              tmp(idd_zpbl4)=+ptype*z(4)
              tmp(idd_zpbl5)=+ptype*z(5)
              tmp(idd_zpbl6)=+ptype*z(6)
              tmp(idd_zpbl7)=+ptype*z(7)
              tmp(idd_zpbl8)=+ptype*z(8)

              tmp(idd_uabl1)=+ptype*uabl(1,i,j,itype)
              tmp(idd_uabl2)=+ptype*uabl(2,i,j,itype)
              tmp(idd_uabl3)=+ptype*uabl(3,i,j,itype)
              tmp(idd_uabl4)=+ptype*uabl(4,i,j,itype)
              tmp(idd_uabl5)=+ptype*uabl(5,i,j,itype)
              tmp(idd_uabl6)=+ptype*uabl(6,i,j,itype)
              tmp(idd_uabl7)=+ptype*uabl(7,i,j,itype)
              tmp(idd_uabl8)=+ptype*uabl(8,i,j,itype)

              tmp(idd_vabl1)=+ptype*vabl(1,i,j,itype)
              tmp(idd_vabl2)=+ptype*vabl(2,i,j,itype)
              tmp(idd_vabl3)=+ptype*vabl(3,i,j,itype)
              tmp(idd_vabl4)=+ptype*vabl(4,i,j,itype)
              tmp(idd_vabl5)=+ptype*vabl(5,i,j,itype)
              tmp(idd_vabl6)=+ptype*vabl(6,i,j,itype)
              tmp(idd_vabl7)=+ptype*vabl(7,i,j,itype)
              tmp(idd_vabl8)=+ptype*vabl(8,i,j,itype)

              tmp(idd_uvabl1)=
     *             +ptype*sqrt( uabl(1,i,j,itype)*uabl(1,i,j,itype)
     *             +            vabl(1,i,j,itype)*vabl(1,i,j,itype))
              tmp(idd_uvabl2)=
     *             +ptype*sqrt( uabl(2,i,j,itype)*uabl(2,i,j,itype)
     *             +            vabl(2,i,j,itype)*vabl(2,i,j,itype))
              tmp(idd_uvabl3)=
     *             +ptype*sqrt( uabl(3,i,j,itype)*uabl(3,i,j,itype)
     *             +            vabl(3,i,j,itype)*vabl(3,i,j,itype))
              tmp(idd_uvabl4)=
     *             +ptype*sqrt( uabl(4,i,j,itype)*uabl(4,i,j,itype)
     *             +            vabl(4,i,j,itype)*vabl(4,i,j,itype))
              tmp(idd_uvabl5)=
     *             +ptype*sqrt( uabl(5,i,j,itype)*uabl(5,i,j,itype)
     *             +            vabl(5,i,j,itype)*vabl(5,i,j,itype))
              tmp(idd_uvabl6)=
     *             +ptype*sqrt( uabl(6,i,j,itype)*uabl(6,i,j,itype)
     *             +            vabl(6,i,j,itype)*vabl(6,i,j,itype))
              tmp(idd_uvabl7)=
     *             +ptype*sqrt( uabl(7,i,j,itype)*uabl(7,i,j,itype)
     *             +            vabl(7,i,j,itype)*vabl(7,i,j,itype))
              tmp(idd_uvabl8)=
     *             +ptype*sqrt( uabl(8,i,j,itype)*uabl(8,i,j,itype)
     *             +            vabl(8,i,j,itype)*vabl(8,i,j,itype))

              tmp(idd_tabl1)=+ptype*tabl(1,i,j,itype)
              tmp(idd_tabl2)=+ptype*tabl(2,i,j,itype)
              tmp(idd_tabl3)=+ptype*tabl(3,i,j,itype)
              tmp(idd_tabl4)=+ptype*tabl(4,i,j,itype)
              tmp(idd_tabl5)=+ptype*tabl(5,i,j,itype)
              tmp(idd_tabl6)=+ptype*tabl(6,i,j,itype)
              tmp(idd_tabl7)=+ptype*tabl(7,i,j,itype)
              tmp(idd_tabl8)=+ptype*tabl(8,i,j,itype)

              tmp(idd_qabl1)=+ptype*qabl(1,i,j,itype)
              tmp(idd_qabl2)=+ptype*qabl(2,i,j,itype)
              tmp(idd_qabl3)=+ptype*qabl(3,i,j,itype)
              tmp(idd_qabl4)=+ptype*qabl(4,i,j,itype)
              tmp(idd_qabl5)=+ptype*qabl(5,i,j,itype)
              tmp(idd_qabl6)=+ptype*qabl(6,i,j,itype)
              tmp(idd_qabl7)=+ptype*qabl(7,i,j,itype)
              tmp(idd_qabl8)=+ptype*qabl(8,i,j,itype)

              tmp(idd_zhat1)=+ptype*zhat(1)
              tmp(idd_zhat2)=+ptype*zhat(2)
              tmp(idd_zhat3)=+ptype*zhat(3)
              tmp(idd_zhat4)=+ptype*zhat(4)
              tmp(idd_zhat5)=+ptype*zhat(5)
              tmp(idd_zhat6)=+ptype*zhat(6)
              tmp(idd_zhat7)=+ptype*zhat(7)

              tmp(idd_e1)=+eabl(1,i,j,itype)*ptype
              tmp(idd_e2)=+eabl(2,i,j,itype)*ptype
              tmp(idd_e3)=+eabl(3,i,j,itype)*ptype
              tmp(idd_e4)=+eabl(4,i,j,itype)*ptype
              tmp(idd_e5)=+eabl(5,i,j,itype)*ptype
              tmp(idd_e6)=+eabl(6,i,j,itype)*ptype
              tmp(idd_e7)=+eabl(7,i,j,itype)*ptype

              tmp(idd_km1)=+ptype*km(1)
              tmp(idd_km2)=+ptype*km(2)
              tmp(idd_km3)=+ptype*km(3)
              tmp(idd_km4)=+ptype*km(4)
              tmp(idd_km5)=+ptype*km(5)
              tmp(idd_km6)=+ptype*km(6)
              tmp(idd_km7)=+ptype*km(7)

              tmp(idd_ri1)=+ptype*gh(1)/(gm(1)+1.d-20)
              tmp(idd_ri2)=+ptype*gh(2)/(gm(2)+1.d-20)
              tmp(idd_ri3)=+ptype*gh(3)/(gm(3)+1.d-20)
              tmp(idd_ri4)=+ptype*gh(4)/(gm(4)+1.d-20)
              tmp(idd_ri5)=+ptype*gh(5)/(gm(5)+1.d-20)
              tmp(idd_ri6)=+ptype*gh(6)/(gm(6)+1.d-20)
              tmp(idd_ri7)=+ptype*gh(7)/(gm(7)+1.d-20)

            END IF
#endif

            ADIURN_part(J,idx(:),kr)=ADIURN_part(J,idx(:),kr) +
     *           tmp(idx(:))
#ifndef NO_HDIURN
            HDIURN_part(J,idx(:),kr)=HDIURN_part(J,idx(:),kr) +
     *           tmp(idx(:))
#endif

          end if
        end do
      endif
c**** quantities accumulated for surface type tables in diagj
      aj(j,j_evap ,itearth)=aj(j,j_evap ,itearth)+ aevap*ptype
      aj(j,j_trhdt,itearth)=aj(j,j_trhdt,itearth)+trhdt*ptype
      aj(j,j_evhdt,itearth)=aj(j,j_evhdt,itearth)+evhdt*ptype
      aj(j,j_shdt ,itearth)=aj(j,j_shdt ,itearth)+ shdt*ptype
      aj(j,j_erun ,itearth)=aj(j,j_erun ,itearth)+(aeruns+aerunu)*ptype
      aj(j,j_run  ,itearth)=aj(j,j_run  ,itearth)+(aruns+arunu)*ptype
      if(moddsf.eq.0) then
        aj(j,j_tsrf,itearth)=aj(j,j_tsrf,itearth)+(ts-tf)*ptype
        aj(j,j_tg1 ,itearth)=aj(j,j_tg1 ,itearth)+    tg1*ptype
        aj(j,j_tg2 ,itearth)=aj(j,j_tg2 ,itearth)+  tg2av*ptype
        aj(j,j_type,itearth)=aj(j,j_type,itearth)+        ptype
      end if

      end subroutine ghy_diag

c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine snow_cover( fract_snow, snow_water, top_dev )
!@sum computes snow cover from snow water eq. and topography
!@var fract_snow snow cover fraction (0-1)
!@var snow_water snow water equivalent (m)
!@var top_dev standard deviation of the surface elevation
      use DOMAIN_DECOMP, only : GRID, GET
      use constant, only : teeny
      real*8, intent(out) :: fract_snow
      real*8, intent(in) :: snow_water, top_dev

      ! using formula from the paper by A. Roesch et al
      ! (Climate Dynamics (2001), 17: 933-946)
      fract_snow =
ccc     $     .95d0 * tanh( 100.d0 * snow_water ) *
ccc                               currently using only topography part
     $     sqrt ( 1000.d0 * snow_water /
     $     (1000.d0 * snow_water + teeny + snow_cover_coef * top_dev) )

      end subroutine snow_cover


      subroutine init_gh(dtsurf,redogh,inisnow,istart)
c**** modifications needed for split of bare soils into 2 types
      use filemanager
      use param
      use constant, only : twopi,rhow,edpery,sha,lhe,tf
      use DOMAIN_DECOMP, only : GRID, DIST_GRID
      use DOMAIN_DECOMP, only : GET,READT_PARALLEL, DREAD_PARALLEL
      use DOMAIN_DECOMP, only : CHECKSUM, HERE, CHECKSUM_COLUMN
      use DOMAIN_DECOMP, only : GLOBALSUM
      use model_com, only : fearth0,itime,nday,jeq,jyear,fland,flice
      use lakes_com, only : flake
      use diag_com, only : npts,icon_wtg,icon_htg,conpt0
      use sle001
#ifdef TRACERS_WATER
      use tracer_com, only : ntm,tr_wd_TYPE,nwater,itime_tr0,needtrs
      use fluxes, only : gtracer
      use veg_com, only:  afb, avh
#endif
      use fluxes, only : gtemp
      use ghy_com
      use dynamics, only : pedn
      use snow_drvm, only : snow_cover_coef2=>snow_cover_coef
     &     ,snow_cover_same_as_rad
      use veg_drv, only : init_vegetation
      use veg_com, only : vdata

      implicit none

      real*8, intent(in) :: dtsurf
      integer, intent(in) :: istart
      logical, intent(in) :: redogh, inisnow
      integer iu_soil,iu_top_index
      integer jday
      real*8 snowdp,wtr1,wtr2,ace1,ace2,tg1,tg2
      logical :: qcon(npts)
      integer i, j, ibv
      real*8 pearth
      logical ghy_data_missing
      character conpt(npts)*10
#ifdef TRACERS_WATER
      real*8 trsoil_tot,wsoil_tot,fm
#endif
c****
cgsfc      REAL*8::TEMP_LOCAL(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,11*NGM+1)
      REAL*8, Allocatable, DIMENSION(:,:,:) :: TEMP_LOCAL
c**** contents of TEMP_LOCAL used for reading the following in a block:
c****       1 -   ngm   dz(ngm)
c****   ngm+1 - 6*ngm   q(is,ngm)
c**** 6*ngm+1 - 11*ngm   qk(is,ngm)
c**** 11*ngm+1           sl
      real*8, external :: qsat
!@dbparam ghy_default_data if == 1 reset all GHY data to defaults
!@+ (do not read it from files)
      integer :: ghy_default_data = 0

      real*8 :: evap_max_ij_sum
C**** define local grid
      integer J_0, J_1
      integer J_0H, J_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

c**** set conservation diagnostics for ground water mass and energy
      conpt=conpt0
      conpt(4)="EARTH"
      qcon=(/ .false., .false., .false., .true., .false., .false.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,conpt,"GRND WTR","(kg/m^2)        ",
     *     "(10^-9 kg/s/m^2)",1d0,1d9,icon_wtg)
      qcon=(/ .false., .false., .false., .true., .false., .false.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,conpt,"GRND ENG","(10**6 J/m^2)   ",
     *     "(10^-3 W/m^2)   ",1d-6,1d3,icon_htg)

c**** read rundeck parameters
      call sync_param( "snow_cover_coef", snow_cover_coef )
! hack. snow_cover_coef should be moved to snow_drvm
      snow_cover_coef2 = snow_cover_coef
      call sync_param( "snow_cover_same_as_rad", snow_cover_same_as_rad)
      call sync_param( "snoage_def", snoage_def )
      call sync_param( "ghy_default_data", ghy_default_data )

c**** read land surface parameters or use defaults
      if ( ghy_default_data == 0 ) then ! read from files

        !!!if (istart.le.0) return ! avoid reading unneeded files
c**** read soils parameters
        call openunit("SOIL",iu_SOIL,.true.,.true.)
        ALLOCATE(TEMP_LOCAL(IM,J_0H:J_1H,11*NGM+1))
        call DREAD_PARALLEL(grid,iu_SOIL,NAMEUNIT(iu_SOIL),TEMP_LOCAL)
        DZ_IJ(:,:,:)   = TEMP_LOCAL(:,:,1:NGM)
         Q_IJ(:,J_0:J_1,:,:) = RESHAPE( TEMP_LOCAL(:,J_0:J_1,1+NGM:) ,
     *                   (/im,J_1-J_0+1,imt,ngm/) )
        QK_IJ(:,J_0:J_1,:,:) =
     *                 RESHAPE( TEMP_LOCAL(:,J_0:J_1,1+NGM+NGM*IMT:) ,
     *                   (/im,J_1-J_0+1,imt,ngm/) )
        SL_IJ(:,J_0:J_1)  = TEMP_LOCAL(:,J_0:J_1,1+NGM+NGM*IMT+NGM*IMT)
        DEALLOCATE(TEMP_LOCAL)
        call closeunit (iu_SOIL)
c**** read topmodel parameters
        call openunit("TOP_INDEX",iu_TOP_INDEX,.true.,.true.)
        call READT_PARALLEL
     *    (grid,iu_TOP_INDEX,NAMEUNIT(iu_TOP_INDEX),0,top_index_ij,1)
        call READT_PARALLEL
     *    (grid,iu_TOP_INDEX,NAMEUNIT(iu_TOP_INDEX),0,top_dev_ij  ,1)
        call closeunit (iu_TOP_INDEX)
      else  ! reset to default data
        if ( istart>0 .and. istart<10 ) then ! reset all
          call reset_gh_to_defaults( .true. )
        else   ! do not reset ghy prognostic variables
          call reset_gh_to_defaults( .false. )
        endif
        !!!if (istart.le.0) return
      endif

c****
c**** initialize constants
c****
c**** time step for ground hydrology
      dt=dtsurf
c spgsn is the specific gravity of snow
      spgsn=.1d0

ccc read and initialize vegetation here
      call init_vegetation(redogh,istart)

      ! no need to continue computations for postprocessing
      if (istart.le.0) return

c**** check whether ground hydrology data exist at this point.
      ghy_data_missing = .false.
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).gt.0) then
            if ( top_index_ij(i,j).eq.-1. ) then
              print *,"No top_index data: i,j=",i,j,top_index_ij(i,j)
              ghy_data_missing = .true.
            end if
            if ( sum(dz_ij(i,j,1:ngm)).eq.0 ) then
              print *, "No soil data: i,j=",i,j,dz_ij(i,j,1:ngm)
              ghy_data_missing = .true.
            endif
            if (wbare(1,i,j) < 1.d-10 .and. wvege(1,i,j) < 1.d-10) then
              print*,"No gh data in restart file: i,j=",i,j,
     &             wbare(:,i,j),wvege(:,i,j)
              ghy_data_missing = .true.
            endif
          end if
        enddo
      enddo
      if ( ghy_data_missing ) then
        write(6,*) 'Ground Hydrology data is missing at some pts'
        write(6,*) 'If you have a non-standard land mask, please'
        write(6,*) 'consider using extended GH data and rfs file.'
        call stop_model(
     &       'Ground Hydrology data is missing at some cells',255)
      endif

      call hl0

c****
c      print *,' '
c      print *,'soils parameters'
c      sdstnc=100.
c      print *,'interstream distance (m) sdstnc:',sdstnc
c      c1=90.
c      print *,'canopy conductance related parameter c1:',c1
c      prfr=.1d0
c      print *,'fraction (by area) of precipitation prfr:',prfr
c      print *,' '
c****
c code transplanted from subroutine input
c**** recompute ground hydrology data if necessary (new soils data)
      if (redogh) then
        jday=1+mod(itime/nday,365)
        cosday=cos(twopi/edpery*jday)
        sinday=sin(twopi/edpery*jday)

        do j=J_0,J_1
        do i=1,im
          pearth=fearth(i,j)
          if(pearth.le.0.) then

            wbare(:,i,j)=0.
            wvege(:,i,j)=0.
            htbare(:,i,j)=0.
            htvege(:,i,j)=0.
            snowbv(:,i,j)=0.

          else
ccc   ??? remove next 5 lines? -check the old version
            w(1:ngm,1) =   wbare(1:ngm,i,j)
            w(0:ngm,2) =   wvege(0:ngm,i,j)
            ht(0:ngm,1) = htbare(0:ngm,i,j)
            ht(0:ngm,2) = htvege(0:ngm,i,j)
            snowd(1:2) =  snowbv(1:2,i,j)

c**** compute soil heat capacity and ground water saturation gws
            call ghinij (i,j)
c**** fill in soils common blocks
            snowdp=snowe(i,j)/rhow
            wtr1=wearth(i,j)
            ace1=aiearth(i,j)
            tg1 =tearth(i,j)
            wtr2=wtr1
            ace2=ace1
            tg2 =tg1
c           wtr2=gdata(i,j,9)   ! this cannot be right
c           ace2=gdata(i,j,10)
c           tg2 =gdata(i,j,8)
            call ghinht (snowdp, tg1,tg2, wtr1,wtr2, ace1,ace2)

c**** copy soils prognostic quantities to model variables
            wbare(1:ngm,i,j) = w(1:ngm,1)
            wvege(0:ngm,i,j) = w(0:ngm,2)
            htbare(0:ngm,i,j) = ht(0:ngm,1)
            htvege(0:ngm,i,j) = ht(0:ngm,2)
            snowbv(1:2,i,j)   = snowd(1:2)
          end if
        end do
        end do
        write (*,*) 'ground hydrology data was made from ground data'
      end if
c**** set gtemp array
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).gt.0) then
            gtemp(1,4,i,j)=tearth(i,j)
          end if
        end do
      end do

C GISS-ESMF EXCEPTIONAL CASE
C-BMP Global sum on evap_max_ij

ccc if not initialized yet, set evap_max_ij, fr_sat_ij, qg_ij
ccc to something more appropriate

       call globalsum(grid, evap_max_ij,evap_max_ij_sum,
     &                all=.true.)
      if ( evap_max_ij_sum > im*jm-1.d0 ) then ! old default
        do j=J_0,J_1
          do i=1,im
            if ( fearth(i,j) .le. 0.d0 ) cycle
            qg_ij(i,j) = qsat(tearth(i,j)+tf,lhe,pedn(1,i,j))
          enddo
        enddo
        fr_sat_ij(:,:) = 0.d0
        evap_max_ij(:,:) = 0.d0
      endif

ccc   init snow here
ccc hope this is the right place to split first layer into soil
ccc and snow  and to set snow arrays
ccc!!! this should be done only when restarting from an old
ccc!!! restart file (without snow model data)

      if (inisnow) then
        do j=J_0,J_1
        do i=1,im
          pearth=fearth(i,j)
          if(pearth.le.0.) then
            nsn_ij(:,i,j)     = 0
            !isn_ij(:,i,j)     = 0
            dzsn_ij(:,:,i,j)  = 0.
            wsn_ij(:,:,i,j)   = 0.
            hsn_ij(:,:,i,j)   = 0.
            fr_snow_ij(:,i,j) = 0.
          else
            jday=1+mod(itime/nday,365)
            cosday=cos(twopi/edpery*jday)
            sinday=sin(twopi/edpery*jday)

            w(1:ngm,1) =   wbare(1:ngm,i,j)
            w(0:ngm,2) =   wvege(0:ngm,i,j)
            ht(0:ngm,1) = htbare(0:ngm,i,j)
            ht(0:ngm,2) = htvege(0:ngm,i,j)
            snowd(1:2) =  snowbv(1:2,i,j)

            call ghinij (i,j)
            call set_snow

            nsn_ij    (1:2, i, j)         = nsn(1:2)
            !isn_ij    (1:2, i, j)         = isn(1:2)
            dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
            wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
            hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
            fr_snow_ij(1:2, i, j)         = fr_snow(1:2)

c****     copy soils prognostic quantities to model variables
             wbare(1:ngm,i,j) = w(1:ngm,1)
             wvege(0:ngm,i,j) = w(0:ngm,2)
            htbare(0:ngm,i,j) = ht(0:ngm,1)
            htvege(0:ngm,i,j) = ht(0:ngm,2)
            snowbv(1:2,i,j)   = snowd(1:2)

          end if
        end do
        end do
      end if

c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      fr_snow_rad_ij(:,:,:) = 0.d0
      do j=J_0,J_1
        do i=1,im
          if ( fearth(i,j) > 0.d0 ) then
            do ibv=1,2
              call snow_cover(fr_snow_rad_ij(ibv,i,j),
     &             snowbv(ibv,i,j), top_dev_ij(i,j) )
              fr_snow_rad_ij(ibv,i,j) = min (
     &             fr_snow_rad_ij(ibv,i,j), fr_snow_ij(ibv, i, j) )
            enddo
          endif
        enddo
      enddo

      ! land water deficit for changing lake fractions
      call compute_water_deficit(jday)

#ifdef TRACERS_WATER
ccc still not quite correct (assumes fw=1)
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).le.0.d0) cycle
          fb=afb(i,j) ; fv=1.-fb
          fm=1.d0-exp(-snowbv(2,i,j)/((avh(i,j)*spgsn) + 1d-12))
          if ( fm < 1.d-3 ) fm=0.d0
          wsoil_tot=fb*( wbare(1,i,j)*(1.d0-fr_snow_ij(1,i,j))
     &     + wsn_ij(1,1,i,j)*fr_snow_ij(1,i,j) )
     &     + fv*( wvege(0,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))   !*1.d0
     &     + wsn_ij(1,2,i,j)*fm*fr_snow_ij(2,i,j) )
          do n=1,ntm
            if (itime_tr0(n).gt.itime) cycle
            if ( .not. needtrs(n) ) cycle
            ! should also restrict to TYPE=nWATER ?
            if ( wsoil_tot > 1.d-30 ) then
            gtracer(n,4,i,j) = (
     &           fb*( tr_wbare(n,1,i,j)*(1.d0-fr_snow_ij(1,i,j))
     &           + tr_wsn_ij(n,1,1,i,j) )         !*fr_snow_ij(1,i,j)
     &           + fv*( tr_wvege(n,0,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))
     &           + tr_wsn_ij(n,1,2,i,j)*fm ) )    !*fr_snow_ij(2,i,j)
     &           /(rhow*wsoil_tot)
            else
              gtracer(n,4,i,j) = 0.
            end if
          enddo
        end do
      end do
#endif

      return
      end subroutine init_gh


      subroutine reset_gh_to_defaults( reset_prognostic )
      !use model_com, only: vdata
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use ghy_com
      use veg_drv, only : reset_veg_to_defaults
      logical, intent(in) :: reset_prognostic
      integer i,j

C**** define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

ccc ugly, should fix later
      call reset_veg_to_defaults( reset_prognostic )

      do j=J_0,J_1
      do i=1,im

      dz_ij(i,j,1:ngm)= (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      q_ij(i,j,1:imt,1:ngm)=
     &     reshape( (/  0.33491623d+00,  0.52958947d+00,
     &     0.13549370d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.32995611d+00,  0.52192056d+00,  0.14812243d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.32145596d+00,
     &     0.48299056d+00,  0.19555295d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.47638881d+00,  0.40400982d+00,
     &     0.11959970d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.99985123d-01,  0.95771909d-01,  0.41175738d-01,
     &     0.00000000d+00,  0.76306665d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.10000000d+01 /), (/imt,ngm/) )
      qk_ij(i,j,1:imt,1:ngm)=
     &     reshape( (/  0.34238762d+00,  0.52882469d+00,
     &     0.12878728d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.32943058d+00,  0.52857041d+00,  0.14199871d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.30698991d+00,
     &     0.52528000d+00,  0.16772974d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.39890009d+00,  0.43742162d+00,
     &     0.16367787d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.46536058d+00,  0.39922065d+00,  0.13541836d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.10000000d+01 /), (/imt,ngm/) )
      sl_ij(i,j)= 0.22695422d+00
      top_index_ij(i,j)= 0.10832934d+02
      top_dev_ij(i,j)= 0.21665636d+03

      if ( .not. reset_prognostic ) cycle

      snowe(i,j)= 0.65458111d-01
      tearth(i,j)= -0.12476520d+00
      wearth(i,j)=  0.29203081d+02
      aiearth(i,j)=  0.93720329d-01
      wbare(:,i,j) = (/  0.17837750d-01,  0.40924843d-01,
     &     0.77932012d-01,  0.11919649d+00,  0.57237469d-01,
     &     0.10000000d-11 /)
      wvege(:,i,j) = (/  0.10000000d-11,  0.29362259d-01,
     &     0.50065177d-01,  0.82533140d-01,  0.10383620d+00,
     &     0.31552459d-01,  0.10000000d-11 /)
      htbare(:,i,j)= (/  0.00000000d+00, -0.15487181d+07,
     &     -0.50720067d+07,  0.18917623d+07,  0.77174974d+07,
     &     0.21716040d+08,  0.44723067d+08 /)
      htvege(:,i,j)= (/ -0.13991376d+05, -0.53165599d+05,
     &     0.65443775d+06,  0.29276050d+07,  0.81455096d+07,
     &     0.21575081d+08,  0.45952255d+08 /)
      snowbv(:,i,j)= (/  0.00000000d+00,  0.65458111d-04 /)

      enddo
      enddo

      end subroutine reset_gh_to_defaults


      subroutine ghinij (i0,j0)
c**** input:
c**** avh(i,j) - array of vegetation heights
c**** spgsn - specific gravity of snow
c**** output:
c**** vh - vegetation height
c**** snowm - snow masking depth
c**** wfcap - water field capacity of top soil layer, m
c****
      use snow_model, only : i_earth,j_earth
      use sle001, only : dz,qk,ng,zb,zc,q,sl,xklh !spgsn,
     *     ,fb,fv,prs,ijdebug,n
     *     ,thets,thetm,ws,thm,nth,shc,shw,htprs,pr !shcap,shtpr,
     *     ,htpr
     *     ,top_index,top_stdev
     &     ,w,ht,snowd,nsn,dzsn,wsn,hsn,fr_snow
      use ghy_com, only : ngm,imt,nlsn,dz_ij,sl_ij,q_ij,qk_ij
     *     ,top_index_ij,top_dev_ij
     &     ,wbare,wvege,htbare,htvege,snowbv,nsn_ij,dzsn_ij,wsn_ij
     &     ,hsn_ij,fr_snow_ij
      use veg_com, only: afb
      USE DOMAIN_DECOMP, ONLY : GRID, GET
!      use veg_drv, only : veg_set_cell

      implicit none
      integer i0,j0
!      real*8 wfcap
      integer k,ibv,i
      real*8 shtpr
!----------------------------------------------------------------------!
      real*8, parameter :: shcap(imt) = (/2d6,2d6,2d6,2.5d6,2.4d6/)


      ijdebug=i0*1000+j0
      i_earth = i0
      j_earth = j0

ccc extracting ghy prognostic vars
      w(1:ngm,1) =  wbare(1:ngm,i0,j0)
      w(0:ngm,2) =  wvege(0:ngm,i0,j0)
      ht(0:ngm,1) = htbare(0:ngm,i0,j0)
      ht(0:ngm,2) = htvege(0:ngm,i0,j0)
      snowd(1:2)  = snowbv(1:2,i0,j0)
ccc extracting snow variables
      nsn(1:2)          = nsn_ij    (1:2, i0, j0)
      !isn(1:2)          = isn_ij    (1:2, i0, j)
      dzsn(1:nlsn, 1:2) = dzsn_ij   (1:nlsn, 1:2, i0, j0)
      wsn(1:nlsn, 1:2)  = wsn_ij    (1:nlsn, 1:2, i0, j0)
      hsn(1:nlsn, 1:2)  = hsn_ij    (1:nlsn, 1:2, i0, j0)
      fr_snow(1:2)      = fr_snow_ij(1:2, i0, j0)

ccc setting vegetation
 !     call veg_set_cell(i0,j0)

ccc passing topmodel parameters
      top_index = top_index_ij(i0, j0)
      top_stdev = top_dev_ij(i0, j0)
c**** set up layers
      dz(1:ngm)=dz_ij(i0,j0,1:ngm)
      q(1:imt,1:ngm)=q_ij(i0,j0,1:imt,1:ngm)
      qk(1:imt,1:ngm)=qk_ij(i0,j0,1:imt,1:ngm)
      sl=sl_ij(i0,j0)

      n=0
      do k=1,ngm
        if(dz(k).le.0.) exit
        n=k
      end do

      if(n.le.0) then
         write (99,*) 'ghinij:  n <= 0:  i,j,n=',i0,j0,n,(dz(k),k=1,43)
         call stop_model('stopped in GHY_DRV.f',255)
      end if

c**** calculate the boundaries, based on the thicknesses.
      zb(1)=0.
      do k=1,n
        zb(k+1)=zb(k)-dz(k)
      end do
c**** calculate the layer centers, based on the boundaries.
      do k=1,n
        zc(k)=.5*(zb(k)+zb(k+1))
      end do
c**** fb,fv: bare, vegetated fraction (1=fb+fv)
      fb=afb(i0,j0)
      fv=1.-fb

c****
      do ibv=1,2
        do k=1,n
          thets(k,ibv)=0.
          thetm(k,ibv)=0.
          do i=1,imt-1
            thets(k,ibv)=thets(k,ibv)+q(i,k)*thm(0,i)
            thetm(k,ibv)=thetm(k,ibv)+q(i,k)*thm(nth,i)
          end do
          ws(k,ibv)=thets(k,ibv)*dz(k)
        end do
      end do
!veg      ws(0,2)=.0001d0*alai  ! from call veg_set_cell above
  !    wfcap=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
c****
      call xklh(1)
c****

      do ibv=1,2
        do k=1,n
          shc(k,ibv)=0.
          do i=1,imt
            shc(k,ibv)=shc(k,ibv)+q(i,k)*shcap(i)
          end do
          shc(k,ibv)=(1.-thets(k,ibv))*shc(k,ibv)*dz(k)
        end do
      end do
c****
c shc(0,2) is the heat capacity of the canopy
!veg      aa=ala(1,i0,j0)
!veg      shc(0,2)=(.010d0+.002d0*aa+.001d0*aa**2)*shw
c****
c htpr is the heat of precipitation.
c shtpr is the specific heat of precipitation.
      shtpr=0.
      if(pr.gt.0.)shtpr=htpr/pr
c htprs is the heat of large scale precipitation
      htprs=shtpr*prs
c****
      return
      end subroutine ghinij


      subroutine ghy_save_cell(i,j)
      use sle001, only : w,ht,snowd,nsn,dzsn,wsn,hsn,fr_snow
      use ghy_com, only : ngm,nlsn
     &     ,dz_ij,wbare,wvege,htbare,htvege,snowbv
     &     ,nsn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij
      implicit none
      integer, intent(in) :: i,j

      wbare(1:ngm,i,j) = w(1:ngm,1)
      wvege(0:ngm,i,j) = w(0:ngm,2)
      htbare(0:ngm,i,j) = ht(0:ngm,1)
      htvege(0:ngm,i,j) = ht(0:ngm,2)
      snowbv(1:2,i,j)   = snowd(1:2)
ccc copy snow variables back to storage
      nsn_ij    (1:2, i, j)         = nsn(1:2)
      !isn_ij    (1:2, i, j)         = isn(1:2)
      dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
      wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
      hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
      fr_snow_ij(1:2, i, j)         = fr_snow(1:2)

      end subroutine ghy_save_cell


      subroutine ghinht (snowdp,tg1,tg2,wtr1,wtr2,ace1,ace2)
c**** initializes new ground (w,ht,snw) from old (t,w,ice,snw)
c**** evaluates the heat in the soil layers based on the
c**** temperatures.
c**** input:
c**** w - water in soil layers, m
c**** tp - temperature of layers, c
c**** fice - fraction of ice of layers
c**** fsn - heat of fusion of water
c**** shc - specific heat capacity of soil
c**** shi - specific heat capacity of ice
c**** shw - specific heat capcity of water
c**** snowd - snow depth, equivalent water m
c**** output:
c**** ht - heat in soil layers
c**** add calculation of wfc2
c**** based on combination of layers 2-n, as in retp2
      use sle001, only : tp, ht, w, shc, fice, snowd, ws, fb, fv,
     &    n, dz, fsn, thetm, shi, shw, ijdebug
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      implicit none

      real*8 snowdp,tg1,tg2,wtr1,wtr2,ace1,ace2
      real*8 wfc1, wfc2, wet1, wet2, wmin, fbv
      integer k, ibv, ll

      wfc1=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
      wfc2=0.
      fbv=fb
      do 30 ibv=1,2
      do 20 k=2,n
      wfc2=wfc2+fbv*ws(k,ibv)
   20 continue
      fbv=fv
   30 continue
      wfc1=1000.*wfc1
      wfc2=1000.*wfc2
      fice(0,2)=1.
      fice(1,1)=(ace1+snowdp*1000.)/(wtr1+ace1+snowdp*1000.+1.d-20)
      fice(1,2)=fice(1,1)
      tp(0,2)=tg1
c**** w = snow(if top layer) + wmin + (wmax-wmin)*(wtr+ice)/wfc
      w(0,2)=0.
      do ibv=1,2
        w(1,ibv)=snowdp
        wmin=thetm(1,ibv)*dz(1)
        wet1=(wtr1+ace1)/(wfc1+1.d-20)
        if(wet1.gt.1.) wet1=1.
        w(1,ibv)=w(1,ibv)+wmin+(ws(1,ibv)-wmin)*wet1
        snowd(ibv)=snowdp
        tp(1,ibv)=tg1
        do k=2,n
          fice(k,ibv)=ace2/(wtr2+ace2+1.d-20)
          wmin=thetm(k,ibv)*dz(k)
          wet2=(wtr2+ace2)/(wfc2+1.d-20)
          if(wet2.gt.1.) wet2=1.
          w(k,ibv)=wmin+(ws(k,ibv)-wmin)*wet2
          tp(k,ibv)=tg2
        end do
      end do
c****
!!!      entry ghexht
c****
c**** compute ht (heat w/m+2)
      do ibv=1,2
        ll=2-ibv
        do k=ll,n
          if(tp(k,ibv)) 2,4,6
 2        ht(k,ibv)=tp(k,ibv)*(shc(k,ibv)+w(k,ibv)*shi)-w(k,ibv)*fsn
          cycle
 4        ht(k,ibv)=-fice(k,ibv)*w(k,ibv)*fsn
          cycle
 6        ht(k,ibv)=tp(k,ibv)*(shc(k,ibv)+w(k,ibv)*shw)
        end do
      end do
      if(ijdebug.eq.0)then
       write(99,*)'ghinht id check',ijdebug
       write(99,*)'tg1,tg2',tg1,tg2
       write(99,*)'tp',tp
       write(99,*)'ht',ht
       write(99,*)'w',w
       write(99,*)'wtr1,wtr2',wtr1,wtr2
       write(99,*)'ace1,ace2',ace1,ace2
       write(99,*)'wfc1,wfc2',wfc1,wfc2
       write(99,*)'shc',shc
       write(99,*)'fice',fice
      endif
      return
      end subroutine ghinht

      subroutine retp2 (tg2av,wtr2av,ace2av)
c**** evaluates the mean temperature in the soil layers 2-ngm
c**** as well as the water and ice content.
c**** input:
c**** w - water in soil layers, m
c**** ht - heat in soil layers
c**** fsn - heat of fusion of water
c**** shc - specific heat capacity of soil
c**** shi - specific heat capacity of ice
c**** shw - specific heat capcity of water
c**** output:
c**** tg2av - temperature of layers 2 to ngm, c
c**** ice2av - ice amount in layers 2 to ngm, kg/m+2
c**** wtr2av - water in layers 2 to ngm, kg/m+2
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use sle001
      implicit none
      real*8 tg2av,wtr2av,ace2av, wc,htc,shcc,tpc,ficec,ftp
      integer k, ibv
      tg2av=0.
      wtr2av=0.
      ace2av=0.
      do 3500 ibv=1,2
      wc=0.
      htc=0.
      shcc=0.
      do k=2,n
        wc=wc+w(k,ibv)
        htc=htc+ht(k,ibv)
        shcc=shcc+shc(k,ibv)
      end do
      tpc=0.
      ficec=0.
      if(wc.ne.0.)  ficec=-htc/(fsn*wc)
      if(fsn*wc+htc.ge.0.)go to 3430
      tpc=(htc+wc*fsn)/(shcc+wc*shi)
      ficec=1.
      go to 3440
 3430 if(htc.le.0.) go to 3440
      tpc=htc/(shcc+wc*shw)
      ficec=0.
 3440 continue
      ftp=fb
      if(ibv.eq.2) ftp=fv
      tg2av=tg2av+tpc*ftp
      wtr2av=wtr2av+wc*ftp*1000.*(1.-ficec)
      ace2av=ace2av+wc*ftp*1000.*ficec
 3500 continue
      return
      end subroutine retp2

      subroutine checke(subr)
!@sum  checke checks whether arrays are reasonable over earth
!@auth original development team
!@ver  1.0
      use model_com, only : itime,wfcs
      use geom, only : imaxj
      use ghy_com, only : tearth,wearth,aiearth,snowe,wbare,wvege,htbare
     *     ,htvege,snowbv,ngm,fearth
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      implicit none

      real*8 x,tgl,wtrl,acel
      integer i,j
!@var subr identifies where check was called from
      character*6, intent(in) :: subr

C**** define local grid
      integer I_0, I_1
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

c**** check for nan/inf in earth data
      call check3(wbare(1:ngm,I_0:I_1,J_0:J_1) ,ngm  ,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'wb')
      call check3(wvege(0:ngm,I_0:I_1,J_0:J_1) ,ngm+1,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'wv')
      call check3(htbare(0:ngm,I_0:I_1,J_0:J_1),ngm+1,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'hb')
      call check3(htvege(0:ngm,I_0:I_1,J_0:J_1),ngm+1,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'hv')
      call check3(snowbv(1:ngm,I_0:I_1,J_0:J_1),2    ,
     *                 (I_1-I_0+1),(J_1-J_0+1),subr,'sn')

c**** check for reasonable temperatures over earth
      x=1.001
      do j=J_0,J_1
        do i=1,imaxj(j)
          if (fearth(i,j).gt.0.) then
            tgl=tearth(i,j)
            wtrl=wearth(i,j)
            acel=aiearth(i,j)
            if ((tgl+60.)*(60.-tgl).le.0.) write (6,901) subr,i,j,itime
     *           ,fearth(i,j),'tg1 ',snowe(i,j),tgl,wtrl,acel
            if (wtrl.lt.0..or.acel.lt.0..or.(wtrl+acel).gt.x*wfcs(i
     *           ,j)) write(6,901) subr,i,j,itime,fearth(i,j),'wtr '
     *           ,snowe(i,j),tgl,wtrl,acel,wfcs(i,j)
          end if
        end do
      end do

      return
 901  format ('0gdata off, subr,i,j,i-time,pearth,',a7,2i4,i10,f5.2,1x
     *     ,a4/' snw,x,tg1,wtr1,ice1, wfc1 ',6f12.4)

      end subroutine checke

      subroutine daily_earth(end_of_day)
!@sum  daily_earth performs daily tasks for earth related functions
!@auth original development team
!@ver  1.0
!@calls RDLAI
      use constant, only : rhow,twopi,edpery,tf
      use model_com, only : nday,nisurf,jday,jyear,wfcs
      use veg_com, only : vdata                 !nyk
      use geom, only : imaxj
      use diag_com, only : aij=>aij_loc
     *     ,tdiurn,ij_strngts,ij_dtgdts,ij_tmaxe
     *     ,ij_tdsl,ij_tmnmx,ij_tdcomp, ij_dleaf
      use ghy_com, only : snoage, snoage_def,fearth
      use veg_com, only : almass,aalbveg       !nyk
      use vegetation, only: crops_yr,cond_scheme,vegCO2X_off !nyk
      use surf_albedo, only: albvnh, updsur  !nyk
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use sle001, only : fb,fv,ws
      use veg_drv, only : veg_set_cell

      implicit none
      real*8 tsavg,wfc1
      real*8 aleafmass, aalbveg0, fvp, sfv  !nyk veg ! , aleafmasslast
      integer i,j,itype
      integer northsouth,iv  !nyk
      logical, intent(in) :: end_of_day

C**** define local grid
      integer J_0, J_1

C**** Extract useful local domain parameters from "grid"
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C**** Update vegetation file if necessary  (i.e. if crops_yr=0)
      if(crops_yr.eq.0) call updveg(jyear,.true.)

c**** find leaf-area index & water field capacity for ground layer 1
      if(cond_scheme.eq.2) call updsur (0,jday) ! Update vegn albedos
            !albvnh(9,6,2)=albvnh(1+8veg,6bands,2hemi), band 1 is VIS.
      cosday=cos(twopi/edpery*jday)
      sinday=sin(twopi/edpery*jday)
      do j=J_0,J_1
        if(j.le.jm/2) then      !nyk added northsouth
          northsouth=1          !southern hemisphere
        else
          northsouth=2          !northern hemisphere
        end if
        do i=1,im
          wfcs(i,j)=24.
          if (fearth(i,j).gt.0.) then
            if (cond_scheme.eq.2) then
              aalbveg0 = 0.d0
              sfv=0.d0
              do iv=1,11
                if ( iv==9 .or. iv==10 ) cycle
                fvp=vdata(i,j,iv+1)
                sfv=sfv+fvp
                aalbveg0 = aalbveg0 + fvp*(ALBVNH(iv+1,1,northsouth))
                !write (99,*) 'fvp',fvp
                !write (99,*) 'ALBVNH',ALBVNH(iv+1,1,northsouth)
              end do
              aalbveg(i,j) = 0.08D0
              if(sfv.gt.0.) aalbveg(i,j) = aalbveg0/sfv !nyk
             !write (99,*) 'daily aalbveg', aalbveg(i,j)
            end if

            call ghinij(i,j)
            call veg_set_cell(i,j,.true.)
            wfc1=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
            wfcs(i,j)=rhow*wfc1 ! canopy part changes

            !-----------------------------------------------------------
            !nyk - TEMPORARY calculate change in leaf mass per day
            !get aleafmass(i,j) at jday
            aleafmass=
     $           almass(1,i,j)+cosday*almass(2,i,j)+sinday*almass(3,i,j)

            !Calculate dlmass(i,j) increment from last jday
            !cosdaym1=cos(twopi/edpery*(jday-1))
            !sindaym1=sin(twopi/edpery*(jday-1))
            !aleafmasslast=almass(1,i,j)+cosdaym1*almass(2,i,j)+
!     $      !     sindaym1*almass(3,i,j)
            !accumulate dlmass
            !adlmass = aleafmass - aleafmasslast
            adlmass = aleafmass
            !aij(i,j,ij_dleaf)=aij(i,j,ij_dleaf)+adlmass
            aij(i,j,ij_dleaf)=adlmass  !accumulate just instant. value
            !PRINT '(F4.4)',adlmass                            !DEBUG
            !call stop_model('Just did adlmass',255)           !DEBUG
          end if
        end do
      end do

      if (end_of_day) then
        do j=J_0,J_1
        do i=1,imaxj(j)
c****
c**** increase snow age depending on snoage_def
c****
          if (snoage_def.eq.0) then ! update indep. of ts
            do itype=1,3
              snoage(itype,i,j)=1.+.98d0*snoage(itype,i,j)
            end do
          elseif (snoage_def.eq.1) then ! update if max T>0
            if (tdiurn(i,j,7).gt.0) snoage(1,i,j)=1.+.98d0
     *           *snoage(1,i,j) ! ocean ice (not currently used)
            if (tdiurn(i,j,8).gt.0) snoage(2,i,j)=1.+.98d0
     *           *snoage(2,i,j) ! land ice
            if (tdiurn(i,j,2).gt.0) snoage(3,i,j)=1.+.98d0
     *           *snoage(3,i,j) ! land
          else
            write(6,*) "This snoage_def is not defined: ",snoage_def
            write(6,*) "Please use: 0 (update indep of T)"
            write(6,*) "            1 (update if T>0)"
            call stop_model('stopped in GHY_DRV.f',255)
          end if
          tsavg=tdiurn(i,j,5)/(nday*nisurf)
          if(32.+1.8*tsavg.lt.65.)
     *         aij(i,j,ij_strngts)=aij(i,j,ij_strngts)+(33.-1.8*tsavg)
          aij(i,j,ij_dtgdts)=aij(i,j,ij_dtgdts)+18.*((tdiurn(i,j,2)-
     *         tdiurn(i,j,1))/(tdiurn(i,j,4)-tdiurn(i,j,3)+1.d-20)-1.)
          aij(i,j,ij_tdsl)=aij(i,j,ij_tdsl)+
     *         (tdiurn(i,j,4)-tdiurn(i,j,3))
          aij(i,j,ij_tdcomp)=aij(i,j,ij_tdcomp)+
     *         (tdiurn(i,j,6)-tdiurn(i,j,9))
          aij(i,j,ij_tmaxe)=aij(i,j,ij_tmaxe)+(tdiurn(i,j,4)-tf)
          if (tdiurn(i,j,6).lt.aij(i,j,ij_tmnmx))
     *         aij(i,j,ij_tmnmx)=tdiurn(i,j,6)
        end do
        end do
      end if

      ! land water deficit for changing lake fractions
      call compute_water_deficit(jday)

#ifdef TRACERS_DRYDEP
      CALL RDLAI ! read leaf area indices for tracer dry deposition
#endif

      return
      end subroutine daily_earth

      subroutine ground_e
!@sum  ground_e driver for applying surface fluxes to land fraction
!@auth original development team
!@ver  1.0
      use model_com, only : itearth
      use geom, only : imaxj,dxyp
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use DOMAIN_DECOMP, only : GLOBALSUM
      use ghy_com, only : snowe, tearth,wearth,aiearth,wbare,wvege
     *     ,snowbv,fr_snow_ij,fr_snow_rad_ij, gdeep, dzsn_ij, nsn_ij,
     *     fearth 
      use veg_com, only : afb
      use diag_com, only : aj=>aj_loc,areg,aij=>aij_loc
     *     ,jreg,ij_evap,ij_f0e,ij_evape
     *     ,ij_gwtr,ij_tg1,j_wtr1,j_ace1,j_wtr2,j_ace2
     *     ,j_snow,j_evap,j_type,ij_g01,ij_g07,ij_g04,ij_g10,ij_g28
     *     ,ij_g29,j_rsnow,ij_rsnw,ij_rsit,ij_snow,ij_gice, ij_gwtr1
     &     ,ij_zsnow
      use fluxes, only : e0,e1,evapor,eprec
      implicit none

      real*8 snow,tg1,tg2,f0dt,f1dt,evap,wtr1,wtr2,ace1,ace2
     *     ,pearth,enrgp,scove
      integer i,j,jr,k

C**** Work array for regional diagnostic accumulation
      real*8, DIMENSION(size(AREG,1),6) :: areg_sum
      real*8, DIMENSION(
     &        size(AREG,1),grid%j_strt_halo:grid%j_stop_halo,6 )
     &        :: AREG_PART

C**** define local grid
      integer J_0, J_1, J_0H, J_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0      , J_STOP=J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )

C***Initialize work array
      areg_part(:,J_0H:j_1H,1:6) = 0.
      do j=J_0,J_1
      do i=1,imaxj(j)
      pearth=fearth(i,j)
      jr=jreg(i,j)
      if (pearth.gt.0) then

        snow=snowe(i,j)
        tg1 = tearth(i,j)
        wtr1= wearth(i,j)
        ace1=aiearth(i,j)
        tg2=gdeep(i,j,1)
        wtr2=gdeep(i,j,2)
        ace2=gdeep(i,j,3)
        f0dt=e0(i,j,4)
        f1dt=e1(i,j,4)
        evap=evapor(i,j,4)
        enrgp=eprec(i,j)      ! including latent heat

c**** accumulate diagnostics
c**** the following is the actual snow cover of the snow model
c        scove = pearth *
c     *       ( afb(i,j)*fr_snow_ij(1,i,j)
c     *       + (1.-afb(i,j))*fr_snow_ij(2,i,j) )
c**** the following computes the snow cover as it is used in RAD_DRV.f
        scove = pearth *
     *       ( afb(i,j)*fr_snow_rad_ij(1,i,j)
     *       + (1.-afb(i,j))*fr_snow_rad_ij(2,i,j) )

        !if (snowe(i,j).gt.0.) scove=pearth
        aj(j,j_rsnow,itearth)=aj(j,j_rsnow,itearth)+scove
        areg_part(jr,j,1)=areg_part(jr,j,1)+scove*dxyp(j)
        aij(i,j,ij_rsnw)=aij(i,j,ij_rsnw)+scove
        aij(i,j,ij_snow)=aij(i,j,ij_snow)+snow*pearth
        aij(i,j,ij_rsit)=aij(i,j,ij_rsit)+scove

        aj(j,j_wtr1,itearth)=aj(j,j_wtr1,itearth)+wtr1*pearth
        aj(j,j_ace1,itearth)=aj(j,j_ace1,itearth)+ace1*pearth
        aj(j,j_wtr2,itearth)=aj(j,j_wtr2,itearth)+wtr2*pearth
        aj(j,j_ace2,itearth)=aj(j,j_ace2,itearth)+ace2*pearth
        aj(j,j_snow,itearth)=aj(j,j_snow,itearth)+snow*pearth
        areg_part(jr,j,2)=areg_part(jr,j,2)+snow*pearth*dxyp(j)
        areg_part(jr,j,3)=areg_part(jr,j,3)+wtr1*pearth*dxyp(j)
        areg_part(jr,j,4)=areg_part(jr,j,4)+ace1*pearth*dxyp(j)
        areg_part(jr,j,5)=areg_part(jr,j,5)+wtr2*pearth*dxyp(j)
        areg_part(jr,j,6)=areg_part(jr,j,6)+ace2*pearth*dxyp(j)

        aij(i,j,ij_f0e)  =aij(i,j,ij_f0e)  +f0dt+enrgp
        aij(i,j,ij_gwtr) =aij(i,j,ij_gwtr)+(wtr1+ace1+wtr2+ace2)
        aij(i,j,ij_gwtr1) =aij(i,j,ij_gwtr1)+(wtr1+ace1)
        aij(i,j,ij_gice) =aij(i,j,ij_gice)+(ace1+ace2)
        aij(i,j,ij_evape)=aij(i,j,ij_evape)+evap
        do k=1,3
          aij(i,j,ij_g01+k-1)=aij(i,j,ij_g01+k-1)+wbare(k,i,j)
          aij(i,j,ij_g07+k-1)=aij(i,j,ij_g07+k-1)+wvege(k-1,i,j)
        end do
        aij(i,j,ij_g04)=aij(i,j,ij_g04)+wbare(6,i,j)
        aij(i,j,ij_g10)=aij(i,j,ij_g10)+wvege(6,i,j)
        aij(i,j,ij_g28)=aij(i,j,ij_g28)+snowbv(1,i,j)
        aij(i,j,ij_g29)=aij(i,j,ij_g29)+snowbv(2,i,j)
        aij(i,j,ij_zsnow)=aij(i,j,ij_zsnow) + pearth *
     &       ( afb(i,j)*fr_snow_ij(1,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(1,i,j),1,i,j) )
     &       + (1.-afb(i,j))*fr_snow_ij(2,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(2,i,j),2,i,j) ) )
      end if
c****
      end do
      end do

      call globalsum(grid,areg_part(1:size(areg,1),:,1:6),
     &    areg_sum(1:size(areg,1),1:6), all=.true.)
      areg(1:size(areg,1),j_rsnow)=areg(1:size(areg,1),j_rsnow)
     &    + areg_sum(1:size(areg,1),1)
      areg(1:size(areg,1),j_snow)=areg(1:size(areg,1),j_snow)
     &    + areg_sum(1:size(areg,1),2)
      areg(1:size(areg,1),j_wtr1)=areg(1:size(areg,1),j_wtr1)
     &    + areg_sum(1:size(areg,1),3)
      areg(1:size(areg,1),j_ace1)=areg(1:size(areg,1),j_ace1)
     &    + areg_sum(1:size(areg,1),4)
      areg(1:size(areg,1),j_wtr2)=areg(1:size(areg,1),j_wtr2)
     &    + areg_sum(1:size(areg,1),5)
      areg(1:size(areg,1),j_ace2)=areg(1:size(areg,1),j_ace2)
     &    + areg_sum(1:size(areg,1),6)


      end subroutine ground_e

      subroutine conserv_wtg(waterg)
!@sum  conserv_wtg calculates zonal ground water incl snow
!@auth Gavin Schmidt
!@ver  1.0
      use constant, only : rhow
      use model_com, only : fim
      use geom, only : imaxj
      use ghy_com, only : ngm,wbare,wvege,snowbv,fearth
      use veg_com, only : afb
      USE DOMAIN_DECOMP, ONLY : GRID, GET, HERE
      implicit none
!@var waterg zonal ground water (kg/m^2)
      real*8, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO),intent(out)::
     &     waterg

      integer i,j,n
      real*8 wij,fb

C**** define local grid
      integer :: J_0, J_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      do j=J_0,J_1
        waterg(j)=0
        do i=1,imaxj(j)
          if (fearth(i,j).gt.0) then
            fb=afb(i,j)
            wij=fb*snowbv(1,i,j)+(1.-fb)*(wvege(0,i,j)+snowbv(2,i,j))
            do n=1,ngm
              wij=wij+fb*wbare(n,i,j)+(1.-fb)*wvege(n,i,j)
            end do
            waterg(j)=waterg(j)+fearth(i,j)*wij*rhow
          end if
        end do
      end do
      if (HAVE_SOUTH_POLE) waterg(1) =fim*waterg(1)
      if (HAVE_NORTH_POLE) waterg(jm)=fim*waterg(jm)
c****
      end subroutine conserv_wtg

      subroutine conserv_htg(heatg)
!@sum  conserv_htg calculates zonal ground energy incl. snow energy
!@auth Gavin Schmidt
!@ver  1.0
      use model_com, only : fim
      use geom, only : imaxj, dxyp
      use ghy_com, only : ngm,htbare,htvege,fr_snow_ij,nsn_ij,hsn_ij
     *     ,fearth 
      use veg_com, only : afb
      USE DOMAIN_DECOMP, ONLY : GRID, GET, HERE
      implicit none
!@var heatg zonal ground heat (J/m^2)
      real*8, dimension(grid%j_strt_halo:grid%j_stop_halo) :: heatg

      integer i,j
      real*8 hij,fb,fv

C**** define local grid
      integer J_0, J_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      do j=J_0,J_1
        heatg(j)=0
        do i=1,imaxj(j)
          if (fearth(i,j).le.0) cycle
          fb=afb(i,j)
          fv=(1.d0-fb)
          hij=fb*sum( htbare(1:ngm,i,j) )
     &       +  fv*sum( htvege(0:ngm,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( hsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( hsn_ij(1:nsn_ij(2,i,j),2,i,j))
          heatg(j)=heatg(j)+fearth(i,j)*hij
        end do
      end do
      if (HAVE_SOUTH_POLE) heatg(1) =fim*heatg(1)
      if (HAVE_NORTH_POLE) heatg(jm)=fim*heatg(jm)
c****
ccc debugging ...
ccc      print *,'conserv_htg energy ',
ccc     &     sum(heatg(1:jm)*dxyp(1:jm))/(sum(dxyp(1:jm))*im)
      end subroutine conserv_htg


      end module soil_drv

      subroutine check_ghy_conservation( flag )
ccc debugging program: cam be put at the beginning and at the end
ccc of the 'surface' to check water conservation
      use constant, only : rhow
      use geom, only : imaxj
      use model_com, only : im,jm
      use DOMAIN_DECOMP, only : GRID, GET
      use fluxes, only : prec,evapor,runoe
      use ghy_com, only : ngm,wbare,wvege,htbare,htvege,snowbv,dz_ij
     *     ,fearth 
      use veg_com, only : afb
      implicit none
      integer flag
      real*8 total_water(im,jm), error_water
      real*8, save :: old_total_water(im,jm)
!      real*8 total_energy(im,jm), error_energy
!      real*8, save :: old_total_energy(im,jm)
      integer i,j,n
      real*8 fb,fv
ccc enrgy check not implemented yet ...

C**** define local grid
      integer J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      do j=J_0,J_1
        do i=1,imaxj(j)
          if ( fearth(i,j) <= 0.d0 ) cycle

ccc just checking ...
          do n = 1,ngm
            if ( dz_ij(i,j,n) .le. 0.d0 )
     &           call stop_model('incompatible dz',255)
          enddo

          fb = afb(i,j)
          fv = 1.d0 - fb
          total_water(i,j) = fb*sum( wbare(1:ngm,i,j) )
     &         + fv*sum( wvege(0:ngm,i,j) )
     &         + fb*snowbv(1,i,j) + fv*snowbv(2,i,j)
        end do
      end do

      ! call stop_model('just testing...',255)

      if ( flag == 0 ) then
        old_total_water(:,:) = total_water(:,:)
        return
      endif

      do j=J_0,J_1
        do i=1,imaxj(j)

          !print *,'fearth = ', i, j, fearth(i,j)

          if ( fearth(i,j) <= 0.d0 ) cycle
          fb = afb(i,j)
          fv = 1.d0 - fb
          error_water = ( total_water(i,j) - old_total_water(i,j) )*rhow
     &         - prec(i,j) + evapor(i,j,4) + runoe(i,j)

          !print *, 'err H2O: ', i, j, error_water

  !        if ( abs( error_water ) > 1.d-9 ) print *, 'error'
          if ( abs( error_water ) > 1.d-9 ) call stop_model(  ! was -15
     &         'check_ghy_conservation: water conservation problem',255)

        end do
      end do

      end subroutine check_ghy_conservation


      subroutine compute_water_deficit(jday)
      use constant, only : twopi,edpery,rhow
      use ghy_com, only : ngm,imt,dz_ij,q_ij
     &     ,wbare,wvege,fearth
      use veg_com, only : ala,afb
      use model_com, only : focean
      use sle001, only : thm
      use fluxes, only : DMWLDF
      USE DOMAIN_DECOMP, ONLY : GRID, GET

      implicit none
      integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1
      integer k,ibv,m
      real*8 :: w_tot(2),w_stor(2)
      real*8 :: w(0:ngm,2),dz(ngm),q(imt,ngm)
      real*8 :: cosday,sinday,alai
      real*8 :: fb,fv

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

      cosday=cos(twopi/edpery*jday)
      sinday=sin(twopi/edpery*jday)

      DMWLDF(:,:) = 0.d0

      do j=J_0,J_1
        do i=I_0,I_1

          !if( focean(i,j) >= 1.d0 ) then
          ! this condition should be switched to focean(i,j) >= 1.d0
          ! once all ground arrays are properly initialized for focean(i,j)<1
          if( fearth(i,j) <= 0.d0 ) then
            DMWLDF(i,j) = 0.d0
            cycle
          endif

          w(1:ngm,1) =  wbare(1:ngm,i,j)
          w(0:ngm,2) =  wvege(0:ngm,i,j)
          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

          fb = afb(i,j)
          fv=1.-fb

          w_stor(:) = 0.d0
          w_tot(:) = 0.d0
          do ibv=1,2
            do k=1,ngm
              do m=1,imt-1
                w_stor(ibv) = w_stor(ibv) + q(m,k)*thm(0,m)*dz(k)
              end do
              w_tot(ibv) = w_tot(ibv) + w(k,ibv)
            end do
          end do

          ! include canopy water here
          alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
          alai=max(alai,1.d0)
          w_stor(2) = w_stor(2) + .0001d0*alai
          w_tot(2) = w_tot(2) + w(0,2)

          ! total water deficit on kg/m^2
          DMWLDF(i,j) = fb*(w_stor(1) - w_tot(1))
     &         + fv*(w_stor(2) - w_tot(2))
        enddo
      enddo


      !print *,"DMWLDF=",DMWLDF

      end subroutine compute_water_deficit

