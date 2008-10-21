#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

c******************   TRACERS             ******************************
#ifdef TRACERS_ON
      module ghy_tracers

      use sle001, only :
#ifdef TRACERS_WATER
     &    tr_w,tr_wsn,trpr,trdd,tr_surf,ntg,ntgm,atr_evap,atr_rnff,atr_g
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      use sle001,ONLY : aevap
#endif
      use tracer_com, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
     *     ,trname
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,Ntm_dust
#endif
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE
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
      use trdiag_com, only : taijn=>taijn_loc,tij_surf, tij_surfbv,
     *   taijs=>taijs_loc,ijts_isrc,jls_isrc,tajls=>tajls_loc,
     *   itcon_surf
#ifdef TRACERS_WATER
     *     ,tij_evap,tij_grnd,tij_soil,tij_snow
#endif
#ifdef TRACERS_DRYDEP
     *     ,tij_drydep,tij_gsdep,itcon_dd,dtr_dd
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,nDustEmij,nDustEmjl
     &     ,ijts_spec,nDustEv1ij,nDustEv2ij,nDustWthij
     &     ,jls_spec,nDustEv1jl,nDustEv2jl,nDustWthjl
#endif
#ifdef TRACERS_DUST
     &     ,nDustEm2ij,nDustEm2jl
#endif
      use fluxes, only : trsource,trsrfflx,trcsurf
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer,trprec
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,prec,pprec,pevap,dust_flux_glob
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
      use ghy_com, only : tr_w_ij,tr_wsn_ij, ntixw
#endif

ccc extra stuff which was present in "earth" by default
#ifdef TRACERS_WATER
      use ghy_com, only : ngm,nlsn
      use constant, only : rhow
#endif
      use dynamics, only : byam
      use geom, only : byaxyp

      implicit none
      private

      public ghy_tracers_set_step
      public ghy_tracers_set_cell
      public ghy_tracers_save_cell

      real*8 totflux(ntm)
      integer ntx,ntix(ntm)
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
          if ( tr_wd_TYPE(n)==nWATER .or. tr_wd_TYPE(n)==nPART ) then
            ntg = ntg + 1
            if ( ntg > ntgm ) then
               print *,"ntgm too small. maybe set it to ntm=",ntm
               call stop_model("ghy_drv: ntg > ntgm",255)
            endif
            ntixw(ntg) = n
          endif
#endif
        end if
      end do

      end subroutine ghy_tracers_set_step


      subroutine ghy_tracers_set_cell(i,j,ptype,qg,evap_max,pbl_args
     &  ,qm1)
!@sum tracers code to be called before the i,j cell is processed
      use model_com, only : dtsrc
      use pbl_drv, only : t_pbl_args

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE constant,ONLY : sday
      USE model_com,ONLY : jday,jmon,wfcs
      USE geom,ONLY : axyp
      USE ghy_com,ONLY : snowe,wearth,aiearth
      USE tracers_dust,ONLY : d_dust,ers_data,hbaij,ricntd,src_fnct
     &     ,frclay,frsilt,dryhr,vtrsh
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     ,minfr
#endif
#endif
#ifdef TRACERS_WATER
      !!!use sle001, only : qm1
#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: qg,evap_max,ptype
      type (t_pbl_args), intent(inout) :: pbl_args
      real*8, intent(in) :: qm1
      integer n,nx,nsrc

c**** pass indices of tracer arrays to PBL
      pbl_args%ntix(1:ntm) = ntix(1:ntm)
      pbl_args%ntx = ntx
C**** Set up tracers for PBL calculation if required
      do nx=1,ntx
        n = ntix(nx)
C**** Calculate first layer tracer concentration
        pbl_args%trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*byaxyp(i,j)
      end do

ccc tracers variables
#ifdef TRACERS_WATER
      do nx=1,ntg
        n = ntixw(nx)
        ! prognostic vars
        tr_w(nx,0,1) = 0.d0
        tr_w(nx,1:ngm,1) = tr_w_ij(n,1:ngm,1,i,j)
        tr_w(nx,0:ngm,2) = tr_w_ij(n,0:ngm,2,i,j)
        tr_wsn(nx,1:nlsn,1:2) = tr_wsn_ij(n,1:nlsn, 1:2, i, j)
        ! flux in
        trpr(nx) = (trprec(n,i,j)*byaxyp(i,j))/dtsrc ! kg/m^2 s (in precip)
#ifdef TRACERS_DRYDEP
        trdd(nx) = trdrydep(n,itype,i,j)/dtsrc   ! kg/m^2 s (dry dep.)
#else
        trdd(nx) = 0
#endif
        ! concentration of tracers in atm. water at the surface
        if (qm1.gt.0 .and. tr_wd_TYPE(n)==nWATER) then
          tr_surf(nx) = trm(i,j,1,n)*byaxyp(i,j)*rhow/qm1 ! kg/m^3
        else
          tr_surf(nx) = 0.
        end if
      enddo
#endif

      do nx=1,ntx
        n=ntix(nx)
C**** set defaults
        pbl_args%trsfac(nx)=0.
        totflux(nx)=0.
        pbl_args%trconstflx(nx)=0.
#ifdef TRACERS_WATER
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
C**** The select is used to distinguish water from gases or particle
! select removed because of OMP compiler bug
!        select case (tr_wd_TYPE(n))
!        case (nWATER)
        if (tr_wd_TYPE(n) .eq. nWATER) then
C**** no fractionation from ground (yet)
C**** trsfac and trconstflx are multiplied by cq*ws and QG in PBL
          pbl_args%trsfac(nx)=1.
          pbl_args%trconstflx(nx)=gtracer(n,itype,i,j)
!        case (nGAS, nPART)
        elseif (tr_wd_TYPE(n).eq.nGAS .or. tr_wd_TYPE(n).eq.nPART) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) pbl_args%trsfac(nx) = 1.
          !then multiplied by deposition velocity in PBL
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
          do nsrc=1,ntsurfsrc(n)
            totflux(nx) = totflux(nx)+trsource(i,j,nsrc,n)
          end do
          pbl_args%trconstflx(nx)=totflux(nx)*byaxyp(i,j)   ! kg/m^2/s
#ifdef TRACERS_WATER
!        end select
        end if
#endif
      end do

#ifdef TRACERS_WATER
c**** water tracers are also flux limited
      do nx=1,ntx
        n=ntix(nx)
C       pbl_args%tr_evap_max(nx) = evap_max * trsoil_rat(nx)
        pbl_args%tr_evap_max(nx) = evap_max * gtracer(n,itype,i,j)
#ifdef TRACERS_DRYDEP
        if(dodrydep(n)) pbl_args%tr_evap_max(nx) = 1.d30
#endif
      end do
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      pbl_args%snowe=snowe(i,j)
      pbl_args%wearth=wearth(i,j)
      pbl_args%aiearth=aiearth(i,j)
      pbl_args%wfcs=wfcs(i,j)
      pbl_args%ers_data=ers_data(i,j,jmon)
      pbl_args%src_fnct=src_fnct(i,j)
      pbl_args%frclay=frclay(i,j)
      pbl_args%frsilt=frsilt(i,j)
      pbl_args%vtrsh=vtrsh(i,j)
      pbl_args%dryhr=dryhr(i,j)
      pbl_args%hbaij=hbaij(i,j)
      pbl_args%ricntd=ricntd(i,j)
      pbl_args%pprec=pprec(i,j)
      pbl_args%pevap=pevap(i,j,itype)
c**** prescribed dust emission
      pbl_args%d_dust(1:Ntm_dust)=d_dust(i,j,1:Ntm_dust,jday)/Sday
     &     /axyp(i,j)/ptype
#endif
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      pbl_args%minfr(:)=minfr(i,j,:)
#endif

      end subroutine ghy_tracers_set_cell


      subroutine ghy_tracers_save_cell(i,j,ptype,dtsurf,rhosrf,pbl_args
#ifdef INTERACTIVE_WETLANDS_CH4
     & ,ra_temp,ra_sat,ra_gwet
#endif
     & )
!@sum tracers code to be called after the i,j cell is processed
      use model_com, only : itime,qcheck,nisurf
      use pbl_drv, only : t_pbl_args
      use ghy_com, only : tearth
      use somtq_com, only : mz
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7D_acc
#endif
      USE TRACER_COM
 !     use socpbl, only : dtsurf
      use geom, only : axyp
      use sle001, only : nsn,fb,fv
#ifdef TRACERS_GASEXCH_land_CO2
     &     ,agpp,arauto,asoilresp
#endif
#if (defined TRACERS_DUST) && (defined TRACERS_DRYDEP)
      USE trdiag_com,ONLY : rts_save
#endif
#ifdef BIOGENIC_EMISSIONS
      use trdiag_com,ONLY : ijs_isoprene
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      use constant, only : tf
      use tracer_sources, only : n__temp,n__sat,n__gwet
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only: DTR_AMPe
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE tracers_dust,ONLY : hbaij,ricntd
#endif
#ifdef TRACERS_GASEXCH_land_CO2
      USE FLUXES, only : TRGASEX
#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: ptype,dtsurf,rhosrf
      type (t_pbl_args), intent(in) :: pbl_args
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      integer :: n1
#endif
      integer n,nx
#ifdef TRACERS_DRYDEP
      real*8 tdryd,tdd,td1,rtsdt,rts
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      real*8 trc_flux
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      real*8, intent(IN) :: ra_temp,ra_sat,ra_gwet
#endif

ccc tracers

#ifdef TRACERS_GASEXCH_land_CO2
      !!! hack - assume 1 tracer
      n = 1
      TRGASEX(n,4,I,J) =
     &     (arauto+asoilresp-agpp)/dtsurf
      trsrfflx(i,j,n)=trsrfflx(i,j,n) +
     &     (arauto+asoilresp-agpp)/dtsurf *axyp(i,j)*ptype
      taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n))
     &     + (arauto+asoilresp-agpp) * axyp(i,j)*ptype
#endif

#ifdef TRACERS_WATER
      do nx=1,ntg
        n = ntixw(nx)
        tr_w_ij(n,1:ngm,1,i,j) = tr_w(nx,1:ngm,1)
        tr_w_ij(n,0:ngm,2,i,j) = tr_w(nx,0:ngm,2)
        tr_wsn_ij(n,1:nlsn, 1:2, i, j) = tr_wsn(nx,1:nlsn,1:2)
      enddo
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      hbaij(i,j)=pbl_args%hbaij
      ricntd(i,j)=pbl_args%ricntd
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
     &       atr_evap(nx)/dtsurf *axyp(i,j)*ptype
      enddo
#endif
      DO nx=1,ntx
        n=ntix(nx)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
C**** technicallly some of these are ocean emissions, but if
C**** fixed datasets are used, it can happen over land as well.

        select case (trname(n))
        case ('DMS')
          trc_flux=pbl_args%DMS_flux
        case ('seasalt1', 'M_SSA_SS')
          trc_flux=pbl_args%ss1_flux
        case ('seasalt2', 'M_SSC_SS')
          trc_flux=pbl_args%ss2_flux
        case ('M_SSS_SS')
          trc_flux=(pbl_args%ss1_flux+pbl_args%ss2_flux)
#ifdef TRACERS_AMP
        case ('M_DD1_DU')
          trc_flux=sum(pbl_args%dust_flux(1:2))
        case ('M_DD2_DU')
          trc_flux=sum(pbl_args%dust_flux(3:4))
        case ('M_DDD_DU')
          trc_flux=sum(pbl_args%dust_flux(1:4))
#endif
        case default
          trc_flux=0
        end select

        trsrfflx(i,j,n)=trsrfflx(i,j,n)+
     &       trc_flux*axyp(i,j)*ptype
        taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &       trc_flux*axyp(i,j)*ptype*dtsurf

#ifdef TRACERS_AMP
        DTR_AMPe(j,n)=DTR_AMPe(j,n)+trc_flux*axyp(i,j)*ptype*dtsurf
#else
        tajls(j,1,jls_isrc(1,n)) = tajls(j,1,jls_isrc(1,n))+
     *       trc_flux*axyp(i,j)*ptype*dtsurf   ! why not for all aerosols?
#endif

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
          trsrfflx(i,j,n)=trsrfflx(i,j,n)
     &         +pbl_args%dust_flux(n1)*axyp(i,j)*ptype
          taijs(i,j,ijts_isrc(nDustEmij,n))
     &         =taijs(i,j,ijts_isrc(nDustEmij,n))
     &         +pbl_args%dust_flux(n1)
     &         *axyp(i,j)*ptype*dtsurf
          tajls(j,1,jls_isrc(nDustEmjl,n))
     &         =tajls(j,1,jls_isrc(nDustEmjl,n))
     &         +pbl_args%dust_flux(n1)*axyp(i,j)
     &         *ptype*dtsurf
#ifdef TRACERS_DUST
          IF (imDust == 0) THEN
            taijs(i,j,ijts_isrc(nDustEm2ij,n))
     &           =taijs(i,j,ijts_isrc(nDustEm2ij,n))
     &           +pbl_args%dust_flux2(n1)
     &           *axyp(i,j)*ptype*dtsurf
            tajls(j,1,jls_isrc(nDustEm2jl,n))
     &           =tajls(j,1,jls_isrc(nDustEm2jl,n))
     &           +pbl_args%dust_flux2(n1)*axyp(i,j)*ptype*dtsurf
          END IF
#endif

        END SELECT
#endif

#ifdef TRACERS_DRYDEP

ccc accumulate tracer dry deposition
        if(dodrydep(n)) then
          rts=rhosrf*pbl_args%trs(nx)
          rtsdt=rts*dtsurf                             ! kg*s/m^3
          tdryd=-rtsdt*(pbl_args%dep_vel(n)+pbl_args%gs_vel(n))          ! kg/m2
          tdd = tdryd*axyp(i,j)*ptype                    ! kg
          td1 = (trsrfflx(i,j,n)+totflux(nx))*dtsurf   ! kg
          if (trm(i,j,1,n)+td1+tdd.le.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd earth",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,pbl_args%trs(nx),pbl_args%trtop(nx)
            tdd= -max(trm(i,j,1,n)+td1,0d0)
            tdryd= tdd/(axyp(i,j)*ptype)
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
     &         ptype*rtsdt*pbl_args%dep_vel(n)
          taijn(i,j,tij_gsdep ,n)=taijn(i,j,tij_gsdep ,n) +
     &         ptype*rtsdt* pbl_args%gs_vel(n)
#ifdef TRACERS_COSMO
          if (n .eq. n_Be7) BE7D_acc(i,j)=BE7D_acc(i,j)+ptype*rtsdt
     *         *pbl_args%dep_vel(n)+ptype*rtsdt* pbl_args%gs_vel(n)
#endif
          dtr_dd(j,n,1)=dtr_dd(j,n,1)-
     &         ptype*rtsdt*axyp(i,j)*pbl_args%dep_vel(n)
          dtr_dd(j,n,2)=dtr_dd(j,n,2)-
     &         ptype*rtsdt*axyp(i,j)* pbl_args%gs_vel(n)
        end if
#endif
#ifdef BIOGENIC_EMISSIONS
        select case (trname(n))
        case ('Isoprene')
          trsrfflx(i,j,n)=trsrfflx(i,j,n)+pbl_args%emisop
     &         *axyp(i,j)*ptype
          taijs(i,j,ijs_isoprene)=taijs(i,j,ijs_isoprene)+
     &    pbl_args%emisop*axyp(i,j)*ptype*dtsurf
        end select
#endif
      end do

C**** Save surface tracer concentration whether calculated or not
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if (needtrs(n)) then
            nx=nx+1
            taijn(i,j,tij_surf  ,n) = taijn(i,j,tij_surf  ,n)+
     &           pbl_args%trs(nx)*ptype
            taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)+
     &           pbl_args%trs(nx)*ptype*rhosrf
            trcsurf(i,j,n)=trcsurf(i,j,n)+pbl_args%trs(nx)*ptype
          else
            taijn(i,j,tij_surf  ,n) = taijn(i,j,tij_surf  ,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *byaxyp(i,j),0d0)*ptype
            taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *byaxyp(i,j),0d0)*ptype*rhosrf
            trcsurf(i,j,n)=trcsurf(i,j,n)+max((trm(i,j,1,n)-trmom(mz,i,j
     *           ,1,n))*byam(1,i,j)*byaxyp(i,j),0d0)*ptype
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
     &       fb*(sum( tr_w(nx,1:ngm,1) ))+
     &       fv*(sum( tr_w(nx,0:ngm,2) ))
     *       )
        taijn(i,j,tij_snow,n)=taijn(i,j,tij_snow,n) + (
     &       fb*(sum( tr_wsn(nx,1:nsn(1),1) ))+
     &       fv*(sum( tr_wsn(nx,1:nsn(2),2) ))
     *       )
        if (tr_wd_TYPE(n).eq.nWATER) tajls(j,1,jls_isrc(1,n))=
     *       tajls(j,1,jls_isrc(1,n))+trevapor(n,itype,i,j)*ptype
      enddo
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c ..........
c Accumulates dust events. One diagnostic field for all dust tracers
c ..........
      taijs(i,j,ijts_spec(nDustEv1ij))=taijs(i,j,ijts_spec(nDustEv1ij))
     &     +pbl_args%dust_event1
      tajls(j,1,jls_spec(nDustEv1jl))=tajls(j,1,jls_spec(nDustEv1jl))
     &     +pbl_args%dust_event1
      taijs(i,j,ijts_spec(nDustEv2ij))=taijs(i,j,ijts_spec(nDustEv2ij))
     &     +pbl_args%dust_event2
      tajls(j,1,jls_spec(nDustEv2jl))=tajls(j,1,jls_spec(nDustEv2jl))
     &     +pbl_args%dust_event2
      taijs(i,j,ijts_spec(nDustWthij))=taijs(i,j,ijts_spec(nDustWthij))
     &     +pbl_args%wtrsh*ptype
      tajls(j,1,jls_spec(nDustWthjl))=tajls(j,1,jls_spec(nDustWthjl))
     &     +pbl_args%wtrsh*ptype
#endif

c     ..........
c     save global variables for subdaily diagnostics
c     ..........

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      DO n=1,Ntm_dust
        dust_flux_glob(i,j,n)=pbl_args%dust_flux(n)*ptype
#ifdef TRACERS_DUST
        dust_flux2_glob(i,j,n)=pbl_args%dust_flux2(n)*ptype
#endif
      END DO
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
#ifdef TRACERS_DRYDEP
      DO n=1,Ntm
        IF (dodrydep(n)) THEN
          depo_turb_glob(i,j,itype,n)=ptype*rts*pbl_args%dep_vel(n)
          depo_grav_glob(i,j,itype,n)=ptype*rts*pbl_args%gs_vel(n)
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
      use socpbl, only : npbl=>n
#ifndef USE_ENT
      use veg_drv, only : cosday,sinday
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!      USE pbl_drv,ONLY : dust_flux,dust_flux2,z,km,gh,gm,zhat,lmonin,
!     &     wsubtke,wsubwd,wsubwm,dust_event1,dust_event2,wtrsh
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
     *     ,idd_rifl,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,idd_tabl1
     *     ,idd_qabl1,idd_zhat1,idd_e1,idd_km1,idd_ri1,idd_emis
     *     ,idd_emis2,idd_grav,idd_turb
#endif
      implicit none
      private
      save

      public daily_earth
      public ground_e
      public init_gh
      public earth
      !public conserv_wtg
      !public conserv_htg
      !public conserv_wtg_1
      public checke
      public snow_cover

      !real*8 cosday,sinday
      !real*8 cosdaym1, sindaym1               !nyk TEMPORARY for jday-1
#ifndef USE_ENT
      real*8 adlmass          ! accumulator for dleafmass in daily_earth
#endif
      real*8 spgsn !@var spgsn specific gravity of snow
!@dbparam snow_cover_coef coefficient for topography variance in
!@+       snow cover parameterisation for albedo
      real*8 :: snow_cover_coef = .15d0
#ifdef USE_ENT
      integer :: vegCO2X_off = 0
#endif

      integer variable_lk

      ! Indexes used for adiurn and hdiurn
      INTEGER, PARAMETER :: n_idx = 22
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER,PARAMETER :: n_idxd=6
#else
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idxd=12+10*npbl
#endif
#endif

      contains

      subroutine earth (ns,moddsf,moddd)
!@sum EARTH calculates surface fluxes of sensible heat,
!@+   evaporation, thermal radiation, and momentum drag over earth.
!@auth I. Alienov/F. Abramopolous
c****
      use constant, only : grav,rgas,lhe,lhs
     *     ,sha,tf,rhow,deltx,gasc,stbo
      use model_com, only : t,p,q,dtsrc,nisurf,dsig,jdate
     *     ,jday,jhour,nday,itime,u,v
     *     ,im,jm
     &     ,Jyear,Jmon,Jday,Jdate,Jhour
#ifdef SCM
     &     ,I_TARG,J_TARG,NSTEPSCM
#endif
#ifdef SCM
      use SCMCOM , only : SCM_SURFACE_FLAG,ASH,ALH,iu_scm_prt
      use SCMDIAG, only : EVPFLX,SHFLX
#endif
      use DOMAIN_DECOMP, only : GRID, GET
      use DOMAIN_DECOMP, only : HALO_UPDATE, CHECKSUM, NORTH
      use DOMAIN_DECOMP, only : GLOBALSUM, AM_I_ROOT
      use geom, only : imaxj
      use dynamics, only : pmid,pk,pek,pedn,am
      use rad_com, only : trhr,fsf, cosz1,trsurf
#ifdef USE_ENT
     &     ,FSRDIR, SRVISSURF,CO2X, CO2ppm
#endif
      !use surf_albedo, only: albvnh   ! added 5/23/03 from RADIATION.f
      !albvnh(9,6,2)=albvnh(sand+8veg,6bands,2hemi) - only need 1st band
      use sle001, only : advnc,evap_limits,
!!!     &    pr,htpr,prs,htprs,
     &     w,snowd,tp,fice,
     &    ashg,alhg, !fv,fb,
     &    aevap,abetad,
     &    aruns,arunu,aeruns,aerunu,
     &    tbcs,tsns,af0dt,af1dt
!!!     &    qm1,qs,
!!!     &    pres,rho,ts,ch,srht,trht
!!!     &   ,vs,vs0,tprime,qprime
#ifndef USE_ENT
      use veg_drv, only: veg_save_cell,veg_set_cell
#endif
      use fluxes, only : dth1,dq1,uflux1,vflux1,e0,e1,evapor,prec,eprec
     *     ,runoe,erunoe,gtemp,precss,gtempr
      use ghy_com, only : snowbv, fearth,
     &     fr_snow_ij,
     *     canopy_temp_ij,snowe,tearth,tsns_ij,wearth,aiearth,
     &     evap_max_ij, fr_sat_ij, qg_ij, fr_snow_rad_ij,top_dev_ij
#ifndef USE_ENT
      use vegetation, only :
     &    veg_srht=>srht,veg_pres=>pres,veg_ch=>ch,veg_ws=>vsm !ia
#endif
      use pbl_drv, only : pbl, t_pbl_args

      use snow_drvm, only : snow_cover_same_as_rad

      use diag_com , only : HR_IN_DAY,HR_IN_MONTH,NDIUVAR
     &     ,NDIUPT,adiurn_dust,ijdd
#ifdef TRACERS_ON
      use ghy_tracers, only : ghy_tracers_set_step,ghy_tracers_set_cell,
     &     ghy_tracers_save_cell
#endif
#ifdef WATER_PROPORTIONAL
      use tracer_com, only : ntm,trm
      use fluxes, only : gtracer,trevapor,trsrfflx
      use geom, only : axyp
      use dynamics, only : am
      use pblcom, only : qabl,trabl
#endif
#ifdef USE_ENT
      use ent_com, only : entcells
      !use ent_mod, only : ent_prescribe_vegupdate
      use ent_drv, only : update_vegetation_data
#endif


      use ghy_com, only : ngm,imt,nlsn,LS_NFRAC,dz_ij,sl_ij,q_ij,qk_ij
     *     ,top_index_ij,top_dev_ij
     &     ,w_ij,ht_ij,snowbv,nsn_ij,dzsn_ij,wsn_ij
     &     ,hsn_ij,fr_snow_ij,shc_soil_texture
#ifdef USE_ENT
     &     ,Qf_ij
#endif
#ifndef USE_ENT
      use vegetation, only : t_vegcell
#endif
      implicit none

      integer, intent(in) :: ns,moddsf,moddd
      integer i,j,itype,ibv
      real*8 shdt,evhdt,rcdmws,rcdhws
     *     ,cdq,cdm,cdh,elhx,tg,srheat,tg1,ptype,trheat    !,dhgs
     *     ,rhosrf,ma1,tfs,th1,thv1,p1k,psk,ps,pij
     *     ,spring,q1,dlwdt
      real*8 fb,fv
!@var rhosrf0 estimated surface air density
      real*8 rhosrf0


      real*8 qsat
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
#ifndef USE_ENT
      type (t_vegcell) vegcell
#endif
      real*8 qg_sat,ts,qs

      INTEGER :: idx(n_idx)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER :: idxd(n_idxd)
#endif
      INTEGER :: ii, ivar, kr
#ifdef TRACERS_DUST
      INTEGER :: n,n1
#endif
      REAL*8 :: tmp(NDIUVAR)
#ifdef USE_ENT
      real*8 Ca !@Ca concentration of CO2 at surface (mol/m3)
      real*8 vis_rad, direct_vis_rad, cos_zen_angle
      integer hemi(1:IM,grid%J_STRT:grid%J_STOP)
      integer :: JEQUATOR=JM/2
#endif
#ifdef WATER_PROPORTIONAL
      real*8, dimension(ntm) :: conc1
      integer :: lpbl,itr
      real*8 :: trconcflx
#endif

C****   define local grid
      integer J_0, J_1, J_0H, J_1H ,I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

    ! dtsurf=dtsrc/nisurf

      spring=-1.
      if((jday.ge.32).and.(jday.le.212)) spring=1.

c****
c**** outside loop over time steps, executed nisurf times every hour
c****
ccc set i,j - independent stuff for tracers
#ifdef TRACERS_ON
      call ghy_tracers_set_step
#endif

c****
c**** outside loop over j and i, executed once for each grid point
c****
C**** halo update u and v for distributed parallelization
       call halo_update(grid, U, from=NORTH)
       call halo_update(grid, V, from=NORTH)

       idx =
     &     (/ idd_ts,  idd_tg1, idd_qs,  idd_qg,  idd_swg,
     &        idd_lwg, idd_sh,  idd_lh,  idd_hz0, idd_ug,
     &        idd_vg,  idd_wg,  idd_us,  idd_vs,  idd_ws,
     &        idd_cia, idd_cm,  idd_ch,  idd_cq,  idd_eds,
     &        idd_dbl, idd_ev /)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IF (adiurn_dust == 1) THEN
        idxd=(/ idd_wsgcm,idd_wspdf,idd_wtke,idd_wd,idd_wm,idd_wtrsh
#ifdef TRACERS_DUST
     &      ,idd_emis,   idd_emis2,  idd_turb,   idd_grav,   idd_ws2
     &      ,idd_ustar,  idd_us3,    idd_stress, idd_lmon,   idd_rifl,
     &       (idd_zpbl1+i-1,i=1,npbl), (idd_uabl1+i-1,i=1,npbl),
     *       (idd_vabl1+i-1,i=1,npbl), (idd_uvabl1+i-1,i=1,npbl),
     *       (idd_tabl1+i-1,i=1,npbl), (idd_qabl1+i-1,i=1,npbl),
     *       (idd_zhat1+i-1,i=1,npbl-1), (idd_e1+i-1,i=1,npbl-1),
     *       (idd_km1+i-1,i=1,npbl-1), (idd_ri1+i-1,i=1,npbl-1)
#endif
     &      /)
      END IF
#endif

#ifdef USE_ENT
      !--- at the moment update vegetation every time step
      hemi(:,max(JEQUATOR+1,J_0):J_1) = 1
      hemi(:,J_0:min(JEQUATOR,J_1)) = -1
  !    call ent_prescribe_vegupdate(dtsrc/nisurf,entcells,
  !   &     hemi,jday,jyear,.false.,.false.)

      call update_vegetation_data( entcells,
     &     im, jm, 1, im, J_0, J_1, jday, jyear )
#endif

!$OMP  PARALLEL DO PRIVATE
!$OMP*  (ELHX,EVHDT, CDM,CDH,CDQ,
!$OMP*   I,ITYPE,ibv, J, MA1,PIJ,PSK,PS,P1K,PTYPE, QG,
!$OMP*   QG_NSAT, RHOSRF,RHOSRF0,RCDMWS,RCDHWS,SRHEAT,SHDT,dlwdt,
!$OMP*   TRHEAT, TH1,TFS,THV1,TG1,TG,q1,pbl_args,qg_sat,jr,kr,tmp,
!$OMP*   fb,fv,vegcell,ts,qs
#ifdef TRACERS_DUST
!$OMP*   ,n,n1
#endif
!$OMP*   )
!$OMP*   SHARED(ns,moddsf,moddd)
!$OMP*   SCHEDULE(DYNAMIC,2)

      loop_j: do j=J_0,J_1
  !    hemi=1.
  !    if(j.le.jm/2) hemi=-1.

c**** conditions at the south/north pole
  !    pole= ( j.eq.1 .or. j.eq.jm )

      loop_i: do i=I_0,imaxj(j)

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
      !!!qm1=q1*ma1
c     rhosrf=100.*ps/(rgas*tsv)
c     rhosrf=100.*ps/(rgas*tkv)
c****
c**** earth
c****
      if (ptype.le.0.) then
  !      ipbl(4,i,j)=0
        cycle loop_i
      endif
      itype=4

      !!!pr=prec(i,j)/(dtsrc*rhow)
C**** This variable was originally assoicated with super-saturated
C**** large-scale precip, but has not worked for many moons.
C**** If you want to reinstate it, uncomment this calculation.
c      prs=precss(i,j)/(dtsrc*rhow)
      !!!prs=0.
      !!!htpr=eprec(i,j)/dtsrc

ccc tracers variables

      tg1 = tsns_ij(i,j)
      srheat=fsf(itype,i,j)*cosz1(i,j)
c****
c**** boundary layer interaction
c****
      pbl_args%zs1=.5d-2*rgas*pbl_args%tkv*ma1/pmid(1,i,j)
  !    zs1=.5d-2*rgas*tkv*ma1/pmid(1,i,j)

c**** loop over ground time steps
      tg=tg1+tf
      elhx=lhe
      if(tg1.lt.0.)  elhx=lhs
      pbl_args%qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step
      pbl_args%tg=tg
      pbl_args%elhx=elhx
      pbl_args%qsol=srheat   ! solar heating
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
      pbl_args%tr4 = gtempr(4,I,J)**4

      rhosrf0=100.*ps/(rgas*pbl_args%tgv) ! estimated surface density
C**** Obviously there are no ocean currents for earth points, but
C**** variables set for consistency with surfce
      pbl_args%uocean=0 ; pbl_args%vocean=0
      pbl_args%ocean = .FALSE.
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
      call ghy_tracers_set_cell(i,j,ptype,pbl_args%qg_sat
     &     ,pbl_args%evap_max,pbl_args, q1*ma1)
#endif
#ifdef WATER_PROPORTIONAL
#ifdef TRACERS_ON
      call stop_model('GHY: should not get here',255)
#endif
      pbl_args%ntx = 0  ! tracers are activated in PBL, but GHY has none
#endif
      call pbl(i,j,itype,ptype,pbl_args)
c****
      cdm = pbl_args%cm ! cmgs(itype,i,j)
      cdh = pbl_args%ch ! chgs(itype,i,j)
      cdq = pbl_args%cq ! cqgs(itype,i,j)
c***********************************************************************
c**** calculate qs
      qs=pbl_args%qsrf
      ts=pbl_args%tsv/(1.+qs*deltx)
c**** calculate rhosrf*cdm*ws
      rhosrf=100.*ps/(rgas*pbl_args%tsv)
      rcdmws=cdm*pbl_args%ws*rhosrf
      rcdhws=cdh*pbl_args%ws*rhosrf
      trheat=trhr(0,i,j)
c***********************************************************************
c****
c  define extra variables to be passed in surfc:
      !!!!pres   = ps
      !veg_pres = ps
      !!!rho    = rhosrf
      !!!vs     = pbl_args%ws
      !!!vs0    = pbl_args%ws0
      !!!tprime = pbl_args%tprime
      !!!qprime = pbl_args%qprime
      !veg_ws = pbl_args%ws
      !!!ch     = cdh
      !veg_ch = cdh
      !!!srht   = srheat
      !veg_srht = srheat
      !!!trht   = trheat
#ifndef USE_ENT
      veg_pres = ps
      veg_ws = pbl_args%ws
      veg_ch = cdh
      veg_srht = srheat
#endif
c***********************************************************************
c****
c**** calculate ground fluxes
c     call qsbal
!!! insert htprs here ???

#ifdef USE_ENT
ccc stuff needed for dynamic vegetation
      Ca = CO2X*CO2ppm*(1.0D-06)*ps*100.0/gasc/ts
      if (vegCO2X_off==0) Ca = Ca * CO2X
!!      vis_rad        = SRVISSURF(i,j)
      vis_rad        = SRVISSURF(i,j)*cosz1(i,j)*.82
      direct_vis_rad = vis_rad * FSRDIR(i,j)
      cos_zen_angle = cosz1(i,j)
!!!!      call ghinij (i,j)
      !call veg_set_cell(i,j)
!i move the following part inside GHY
!i      call ent_get_value( i, j,
!i     &     canopy_holding_capacity=ws_can,
!i     &     canopy_heat_capacity=shc_can,
!i     &     fraction_of_vegetated_soil=fv )
!i      fb = 1.d0 - fv
!!!!      call advnc( entcells(i,j), Ca,
!!!!     &     cos_zen_angle, vis_rad, direct_vis_rad )
!!!!      call evap_limits( vegcell,
!!!!     &     .false., evap_max_ij(i,j), fr_sat_ij(i,j) )

      !call veg_save_cell(i,j)
!!!!      call ghy_save_cell(i,j)
#else
      !call ghinij (i,j)
      !snowd(1:2) = snowbv(1:2,i,j)
      call veg_set_cell( vegcell, i,j, ps, pbl_args%tsv/(1.+qs*deltx) )
      !call veg_set_cell(  i,j  )
#endif
! snow / var lakes problem at this cell (underwater snow in initial data)
!      if (i==47 .and. j==33) then
!        print *,i,j
!      endif
      call get_fb_fv( fb, fv, i, j )
      call advnc(
#ifdef USE_ENT
     &     entcells(i,j), Ca,
     &     cos_zen_angle, vis_rad, direct_vis_rad,
     &     Qf_ij(i,j),
#else
     &     vegcell,
#endif
     &     w_ij(0:ngm,1:LS_NFRAC,i,j),
     &     ht_ij(0:ngm,1:LS_NFRAC,i,j),
     &     nsn_ij    (1:2, i, j),
     &     dzsn_ij   (1:nlsn, 1:2, i, j),
     &     wsn_ij    (1:nlsn, 1:2, i, j),
     &     hsn_ij    (1:nlsn, 1:2, i, j),
     &     fr_snow_ij(1:2, i, j),
     &     top_index_ij(i, j),
     &     top_dev_ij(i, j),
     &     dz_ij(i,j,1:ngm),
     &     q_ij(i,j,1:imt,1:ngm),
     &     qk_ij(i,j,1:imt,1:ngm),
     &     sl_ij(i,j),
     &     fb,
     &     fv,
     &     prec(i,j)/(dtsrc*rhow),
     &     eprec(i,j)/dtsrc,
     &     0.d0,
     &     0.d0,
     &     srheat,
     &     trheat,
     &     pbl_args%tsv/(1.+qs*deltx),
     &     pbl_args%qsrf,
     &     ps,
     &     rhosrf,
     &     cdh,
     &     q1*ma1,
     &     pbl_args%ws,
     &     pbl_args%ws0,
     &     pbl_args%tprime,
     &     pbl_args%qprime )




      call evap_limits(
#ifndef USE_ENT
     &     vegcell,
#endif
     &     .false., evap_max_ij(i,j), fr_sat_ij(i,j) )

#ifndef USE_ENT
      if ( fv > 0 ) call veg_save_cell(i,j)
#endif
      !call veg_save_cell(i,j)
      !call ghy_save_cell(i,j)
      if ( fb > 0 ) snowbv(1,i,j)   = snowd(1)
      if ( fv > 0 ) snowbv(2,i,j)   = snowd(2)

      tg1=tsns

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
      tearth(i,j) = tbcs
      tsns_ij(i,j) = tsns
c**** wearth+aiearth are used in radiation only
      wearth(i,j)=1000.*( fb*w(1,1)*(1.-fice(1,1)) +
     &     fv*(w(1,2)*(1.-fice(1,2))+w(0,2)*(1.-fice(0,2))) )
      aiearth(i,j)=1000.*( fb*w(1,1)*fice(1,1) +
     &     fv*(w(1,2)*fice(1,2)+w(0,2)*fice(0,2)) )
      gtemp(1,4,i,j)=tsns_ij(i,j)
      gtempr(4,i,j) =tearth(i,j)+tf
c**** calculate fluxes using implicit time step for non-ocean points
      uflux1(i,j)=uflux1(i,j)+ptype*rcdmws*(pbl_args%us) !-uocean)
      vflux1(i,j)=vflux1(i,j)+ptype*rcdmws*(pbl_args%vs) !-vocean)
c**** accumulate surface fluxes and prognostic and diagnostic quantities
      evapor(i,j,4)=evapor(i,j,4)+aevap
      shdt=-ashg
C**** calculate correction for different TG in radiation and surface
      dLWDT = pbl_args%dtsurf*(TRSURF(ITYPE,I,J)-STBO*(tearth(i,j)+TF)
     *     **4)

#ifdef SCM
c     if SCM use sensible and latent heat fluxes provided by ARM
c        values
      if ((i.eq.I_TARG.and.j.eq.J_TARG).and.SCM_SURFACE_FLAG.eq.1) then
           dth1(i,j)=dth1(i,j)
     &             +ash*pbl_args%dtsurf*ptype/(sha*ma1*p1k)
           dq1(i,j) =dq1(i,j)
     &             +alh*pbl_args%dtsurf*ptype/(ma1*lhe)
           EVPFLX = EVPFLX + ALH*ptype
           SHFLX = SHFLX + ASH*ptype
           write(iu_scm_prt,981) i,ptype,dth1(i,j),dq1(i,j),EVPFLX,SHFLX
 981       format(1x,'EARTH ARM   i ptype dth1 dq1 evpflx shflx ',i5,
     &            f9.4,f9.4,f9.5,f9.5,f9.5)
      else    
           dth1(i,j)=dth1(i,j)-(SHDT+dLWDT)*ptype/(sha*ma1*p1k)
           dq1(i,j) =dq1(i,j)+aevap*ptype/ma1
c          write(iu_scm_prt,982) i,ptype,dth1(i,j),dq1(i,j)
c982       format(1x,'EARTH GCM    i ptype dth1 dq1 ',i5,f9.4,f9.4,f9.5)
      endif
#else
      dth1(i,j)=dth1(i,j)-(SHDT+dLWDT)*ptype/(sha*ma1*p1k)
      dq1(i,j) =dq1(i,j)+aevap*ptype/ma1
#endif
  !    qsavg(i,j)=qsavg(i,j)+qs*ptype
c**** save runoff for addition to lake mass/energy resevoirs
      runoe (i,j)=runoe (i,j)+ aruns+ arunu
      erunoe(i,j)=erunoe(i,j)+aeruns+aerunu
c****
      e0(i,j,4)=e0(i,j,4)+af0dt
      e1(i,j,4)=e1(i,j,4)+af1dt

      call ghy_diag(i,j,jr,kr,tmp,ns,moddsf,moddd
     &     ,rcdmws,cdm,cdh,cdq,qg
     &     ,pbl_args, pbl_args%dtsurf
     &     ,idx
#ifdef TRACERS_DUST
     &     ,n,n1
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,idxd
#endif
     &     )

c      if ( moddd == 0 ) then
c        do kr=1,ndiupt
c          if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then
c            ADIURN(idx(:),kr,ih)=ADIURN(idx(:),kr,ih) + tmp(idx(:))
c#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
c    (defined TRACERS_QUARZHEM)
c            ADIURN(idxd(:),kr,ih)=ADIURN(idxd(:),kr,ih) + tmp(idxd(:))
c#endif
c
c#ifndef NO_HDIURN
c            HDIURN(idx(:),kr,ihm)=HDIURN(idx(:),kr,ihm)+tmp(idx(:))
c#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
c    (defined TRACERS_QUARZHEM)
c            HDIURN(idxd(:),kr,ihm)=HDIURN(idxd(:),kr,ihm)+tmp(idxd(:))
c#endif
c#endif
c          endif
c        enddo
c      endif

c**** update tracers
#ifdef TRACERS_ON
      call ghy_tracers_save_cell(i,j,ptype,pbl_args%dtsurf,rhosrf
     &   ,pbl_args
#ifdef INTERACTIVE_WETLANDS_CH4
     & ,tg1,ts-tf,1.d2*abetad/real(nisurf)
#endif
     & )
#endif

#ifdef WATER_PROPORTIONAL
c calculate fluxes of atmosphere-only water tracers.
c assume 1-way until up/down fluxes of vapor are available
c as a PBL diagnostic.
      conc1(1:ntm) = trm(i,j,1,1:ntm)/(q1*axyp(i,j)*am(1,i,j)+1d-20)
      do itr=1,ntm
        if(aevap.ge.0.) then
          trconcflx = gtracer(itr,itype,i,j)
        else
          trconcflx = conc1(itr)
        endif
        trevapor(itr,itype,i,j) = trevapor(itr,itype,i,j) +
     &       aevap*trconcflx
        trsrfflx(i,j,itr)=trsrfflx(i,j,itr)+axyp(i,j)*ptype*
     &       aevap*trconcflx/(pbl_args%dtsurf)
c fill in pbl profile in case it is used to initialize
c another surface type
        do lpbl=1,npbl
          trabl(lpbl,itr,itype,i,j)=conc1(itr)*qabl(lpbl,itype,i,j)
        enddo
      enddo ! itr
#endif

      end do loop_i
      end do loop_j
!$OMP  END PARALLEL DO

      ! land water deficit for changing lake fractions
      !!! not working with Ent
#ifndef USE_ENT
      call compute_water_deficit
#endif

      return
      end subroutine earth

c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine ghy_diag(i,j,jr,kr,tmp,ns,moddsf,moddd
     &     ,rcdmws,cdm,cdh,cdq,qg, pbl_args, dtsurf
     &     ,idx
#ifdef TRACERS_DUST
     &     ,n,n1
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,idxd
#endif
     &     )

      use diag_com , only : aij=>aij_loc,adiurn=>adiurn_loc
#ifndef NO_HDIURN
     &     ,hdiurn=>hdiurn_loc
#endif
     *     ,tsfrez=>tsfrez_loc,tdiurn,jreg
     *     ,ij_rune, ij_arunu, ij_pevap, ij_shdt, ij_beta, ij_trnfp0
     *     ,ij_srtr, ij_neth, ij_ws, ij_ts, ij_us, ij_vs, ij_taus
     *     ,ij_tauus, ij_tauvs, ij_qs, ij_tg1, ij_evap, j_trhdt, j_shdt
     *     ,j_evhdt,j_evap,j_erun,j_run,j_tsrf,j_type,j_tg1,j_tg2,ij_g05
     *     ,ij_g06,ij_g11,ij_g12,ij_g13,ij_g14,ij_g15,ij_g16,ij_g17
     *     ,ij_gpp,ij_rauto,ij_clab,ij_soilresp,ij_soilCpoolsum
     *     ,ij_pblht,ij_g18,ij_g19,ij_g20,ij_g21,ij_g22,ij_g23
     *     ,ij_g24,ij_g25,ij_g26,ij_g27,ijdd,idd_ts,idd_tg1,idd_qs
     *     ,idd_qg,idd_swg,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg
     *     ,idd_wg,idd_us,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq
     *     ,idd_eds,idd_dbl,idd_ev,tf_day1,tf_last,ndiupt
     *     ,HR_IN_DAY,HR_IN_MONTH,NDIUVAR
     &     ,ij_aflmlt,ij_aeruns,ij_aerunu,ij_fveg
     &     ,ij_htsoil,ij_htsnow,ij_aintrcp,ij_trsdn,ij_trsup,adiurn_dust
     &     ,ij_gusti,ij_mccon
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
     *     ,jday,jhour,nday,itime,modrd,itearth
#ifdef SCM
     &     ,I_TARG,J_TARG
#endif
#ifdef SCM
      use SCMDIAG, only : EVPFLX,SHFLX
      use SCMCOM, only : SCM_SURFACE_FLAG,iu_scm_prt
#endif
      use DOMAIN_DECOMP, only : grid
      use geom, only : axyp,lat2d
      use rad_com, only : trhr,fsf, cosz1

      use sle001, only :
     &     tp
     &    ,fv,fb,atrg,ashg,alhg
     &    ,abetad,abetav,abetat
     &    ,abetap,abetab,abeta
     &    ,acna,acnc,agpp,arauto,aclab,asoilresp,asoilCpoolsum
     &    ,aevap,aevapw,aevapd,aevapb
     &    ,aruns,arunu,aeruns,aerunu,aflmlt,aintercep
     &    ,aepc,aepb,aepp,zw,tbcs,tsns
     &    ,qs,ts,ngr=>n,ht,hsn,fr_snow,nsn
     &    ,tg2av,wtr2av,ace2av

      use ghy_com, only : gdeep, fearth
      USE CLOUDS_COM, only : DDMS

      use pbl_drv, only : t_pbl_args

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
      type (t_pbl_args) :: pbl_args
      real*8, intent(in) :: dtsurf
      INTEGER, INTENT(IN) :: idx(:)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER,INTENT(IN) :: idxd(:)
#endif

      real*8 timez
      real*8 trhdt,tg1,shdt,ptype,srheat,srhdt
      real*8 warmer,spring,trheat,evhdt
      integer, parameter :: itype=4
      integer :: kr,jr,ih,ihm
#ifdef TRACERS_DUST
      INTEGER :: n,n1
#endif
      REAL*8 :: tmp(:)

ccc the following values are returned by PBL
      real*8 us,vs,ws,psi,dbl,khs,ug,vg,wg,gusti
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,wsgcm,wspdf
#endif
      us = pbl_args%us
      vs = pbl_args%vs
      ws = pbl_args%ws
      psi = pbl_args%psi
      dbl = pbl_args%dbl
      khs = pbl_args%khs
      ug = pbl_args%ug
      vg = pbl_args%vg
      wg = pbl_args%wg
      gusti = pbl_args%gusti
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      wsgcm=pbl_args%wsgcm
      wspdf=pbl_args%wspdf
#endif

      timez=jday+(mod(itime,nday)+(ns-1.)/nisurf)/nday ! -1 ??
      if(jday.le.31) timez=timez+365.

      ih=1+jhour
      ihm = ih+(jdate-1)*24

      spring=-1.
      if((jday.ge.32).and.(jday.le.212)) spring=1.

      if(lat2d(i,j).lt.0.)  then
         warmer=-spring
       else
         warmer=spring
      end if

      ptype=fearth(i,j)
      jr=jreg(i,j)

      srheat=fsf(itype,i,j)*cosz1(i,j)
      srhdt=srheat*dtsurf
      trheat=trhr(0,i,j)

      !tg1=tbcs
      tg1=tsns
      shdt=-ashg
      evhdt=-alhg

      aij(i,j,ij_fveg)=aij(i,j,ij_fveg)+fv/nisurf
      aij(i,j,ij_g18)=aij(i,j,ij_g18)+aevapb
      aij(i,j,ij_g19)=aij(i,j,ij_g19)+aevapd
      aij(i,j,ij_g20)=aij(i,j,ij_g20)+aevapw
      aij(i,j,ij_g05)=aij(i,j,ij_g05)+abetab/nisurf
      aij(i,j,ij_g06)=aij(i,j,ij_g06)+abetap/nisurf
      aij(i,j,ij_g11)=aij(i,j,ij_g11)+abeta/nisurf
      aij(i,j,ij_g12)=aij(i,j,ij_g12)+acna/nisurf
      aij(i,j,ij_g13)=aij(i,j,ij_g13)+acnc/nisurf
      aij(i,j,ij_gpp)=aij(i,j,ij_gpp)+agpp
      aij(i,j,ij_rauto)=aij(i,j,ij_rauto)+arauto
      aij(i,j,ij_clab)=aij(i,j,ij_clab)+aclab/nisurf
      aij(i,j,ij_soilresp)=aij(i,j,ij_soilresp)+asoilresp
      aij(i,j,ij_soilCpoolsum)=aij(i,j,ij_soilCpoolsum)
     &     + asoilCpoolsum/nisurf
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
      !!!call retp2 (tg2av,wtr2av,ace2av)
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
      call inc_areg(i,j,jr,j_trhdt,trhdt*ptype*axyp(i,j))
      call inc_areg(i,j,jr,j_shdt , shdt*ptype*axyp(i,j))
      call inc_areg(i,j,jr,j_evhdt,evhdt*ptype*axyp(i,j))
      call inc_areg(i,j,jr,j_evap ,aevap*ptype*axyp(i,j))
      call inc_areg(i,j,jr,j_erun ,(aeruns+aerunu)*ptype*axyp(i,j))
      call inc_areg(i,j,jr,j_run  ,  (aruns+arunu)*ptype*axyp(i,j))
      if ( moddsf == 0 ) then
        call inc_areg(i,j,jr,j_tsrf,(ts-tf)*ptype*axyp(i,j))
        call inc_areg(i,j,jr,j_tg1 , tg1   *ptype*axyp(i,j))
        call inc_areg(i,j,jr,j_tg2 , tg2av *ptype*axyp(i,j))
      end if

#ifdef SCM
      if (J.eq.J_TARG.and.I.eq.I_TARG) then
          if (SCM_SURFACE_FLAG.eq.0) then
              EVPFLX = EVPFLX + aevap*PTYPE/DTSURF
              SHFLX  = SHFLX  + SHDT *PTYPE/DTSURF
c             write(iu_scm_prt,*) 'ghy_drv  evpflx shflx ptype ',
c    &                  EVPFLX,SHFLX,ptype
          endif
      endif
#endif

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
        aij(i,j,ij_taus)=aij(i,j,ij_taus)+rcdmws*ws*ptype
        aij(i,j,ij_tauus)=aij(i,j,ij_tauus)+rcdmws*(us)*ptype !-uocean
        aij(i,j,ij_tauvs)=aij(i,j,ij_tauvs)+rcdmws*(vs)*ptype !-vocean
        aij(i,j,ij_qs)=aij(i,j,ij_qs)+qs*ptype
        aij(i,j,ij_tg1)=aij(i,j,ij_tg1)+tg1*ptype
        aij(i,j,ij_pblht)=aij(i,j,ij_pblht)+dbl*ptype
        if(DDMS(I,J).lt.0.)   ! ddms < 0 for down draft
     *     AIJ(I,J,ij_mccon)=AIJ(I,J,ij_mccon)+ptype
        aij(i,j,ij_gusti)=aij(i,j,ij_gusti)+gusti*ptype
chyd       aij(i,j,ij_arunu)=aij(i,j,ij_arunu)
chyd      *  +   (40.6*psoil+.72*(2.*(tss-tfs)-(qsatss-qss)*lhe/sha))

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        aij(i,j,ij_wsgcm)=aij(i,j,ij_wsgcm)+wsgcm*ptype
        aij(i,j,ij_wspdf)=aij(i,j,ij_wspdf)+wspdf*ptype
        aij(i,j,ij_wdry)=aij(i,j,ij_wdry)+pbl_args%wsubwd*ptype
        aij(i,j,ij_wtke)=aij(i,j,ij_wtke)+pbl_args%wsubtke*ptype
        aij(i,j,ij_wmoist)=aij(i,j,ij_wmoist)+pbl_args%wsubwm*ptype
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
              tmp(idd_wtke)=+pbl_args%wsubtke*ptype
              tmp(idd_wd)=+pbl_args%wsubwd*ptype
              tmp(idd_wm)=+pbl_args%wsubwm*ptype
              tmp(idd_wtrsh)=+pbl_args%wtrsh*ptype
#ifdef TRACERS_DUST
              tmp(idd_emis)=0.D0
              tmp(idd_emis2)=0.D0
              tmp(idd_turb)=0.D0
              tmp(idd_grav)=0.D0
              DO n=1,Ntm_dust
                n1=n_clay+n-1
                tmp(idd_emis)=tmp(idd_emis)+pbl_args%dust_flux(n)*ptype
                tmp(idd_emis2)=tmp(idd_emis2)
     &               +pbl_args%dust_flux2(n)*ptype
#ifdef TRACERS_DRYDEP
                IF (dodrydep(n1)) THEN
                  tmp(idd_turb)=tmp(idd_turb)+ptype*rts_save(i,j)
     &                 *pbl_args%dep_vel(n1)
                  tmp(idd_grav)=tmp(idd_grav)+ptype*rts_save(i,j)
     &                 *pbl_args%gs_vel(n1)
                END IF
#endif
              END DO
              tmp(idd_ws2)=+ws*ws*ptype
              tmp(idd_ustar)=+pbl_args%ustar*ptype
              tmp(idd_us3)=+ptype*pbl_args%ustar*pbl_args%ustar
     &             *pbl_args%ustar
              tmp(idd_stress)=+rcdmws*ws*ptype
              tmp(idd_lmon)=+pbl_args%lmonin*ptype
              tmp(idd_rifl)=
     &             +ptype*grav*(ts-(tg1+tf))*pbl_args%zgs
     &             /(ws*ws*(tg1+tf))

              tmp(idd_zpbl1:idd_zpbl1+npbl-1)=ptype*pbl_args%z(1:npbl)
              tmp(idd_uabl1:idd_uabl1+npbl-1)=
     *             ptype*uabl(1:npbl,itype,i,j)
              tmp(idd_vabl1:idd_vabl1+npbl-1)=
     *             ptype*vabl(1:npbl,itype,i,j)
              tmp(idd_uvabl1:idd_uvabl1+npbl-1)=ptype*sqrt(
     *             uabl(1:npbl,itype,i,j)*uabl(1:npbl,itype,i,j)+
     *             vabl(1:npbl,itype,i,j)*vabl(1:npbl,itype,i,j))
              tmp(idd_tabl1:idd_tabl1+npbl-1)=
     *             ptype*tabl(1:npbl,itype,i,j)
              tmp(idd_qabl1:idd_qabl1+npbl-1)=
     *             ptype*qabl(1:npbl,itype,i,j)
              tmp(idd_zhat1:idd_zhat1+npbl-2)=ptype
     *             *pbl_args%zhat(1:npbl-1)
              tmp(idd_e1:idd_e1+npbl-2)=eabl(1:npbl-1,itype,i,j)*ptype
              tmp(idd_km1:idd_km1+npbl-2)=ptype*pbl_args%km(1:npbl-1)
              tmp(idd_ri1:idd_ri1+npbl-2)=ptype*pbl_args%gh(1:npbl-1)
     *             /(pbl_args%gm(1:npbl-1)+1d-20)
#endif

            END IF
#endif

            ADIURN(idx(:),kr,ih)=ADIURN(idx(:),kr,ih) + tmp(idx(:))
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
            ADIURN(idxd(:),kr,ih)=ADIURN(idxd(:),kr,ih) + tmp(idxd(:))
#endif

#ifndef NO_HDIURN
            HDIURN(idx(:),kr,ihm)=HDIURN(idx(:),kr,ihm)+tmp(idx(:))
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
            HDIURN(idxd(:),kr,ihm)=HDIURN(idxd(:),kr,ihm)+tmp(idxd(:))
#endif
#endif
          end if
        end do
      endif
c**** quantities accumulated for surface type tables in diagj
      call inc_aj(i,j,itearth,j_evap ,aevap*ptype)
      call inc_aj(i,j,itearth,j_trhdt,trhdt*ptype)
      call inc_aj(i,j,itearth,j_evhdt,evhdt*ptype)
      call inc_aj(i,j,itearth,j_shdt , shdt*ptype)
      call inc_aj(i,j,itearth,j_erun ,(aeruns+aerunu)*ptype)
      call inc_aj(i,j,itearth,j_run  ,  (aruns+arunu)*ptype)
      if(moddsf.eq.0) then
        call inc_aj(i,j,itearth,j_tsrf,(ts-tf)*ptype)
        call inc_aj(i,j,itearth,j_tg1 ,    tg1*ptype)
        call inc_aj(i,j,itearth,j_tg2 ,  tg2av*ptype)
        call inc_aj(i,j,itearth,j_type,        ptype)
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


      subroutine init_gh(dtsurf,redogh,inisnow,istart,nl_soil)
c**** modifications needed for split of bare soils into 2 types
      use filemanager
      use param
      use constant, only : twopi,rhow,edpery,sha,lhe,tf,shw_kg=>shw
      use DOMAIN_DECOMP, only : GRID, DIST_GRID
      use DOMAIN_DECOMP, only : GET,READT_PARALLEL, DREAD_PARALLEL
      use DOMAIN_DECOMP, only : CHECKSUM, HERE, CHECKSUM_COLUMN
      use DOMAIN_DECOMP, only : GLOBALSUM
      use model_com, only : fearth0,itime,nday,jyear,fland,flice
     &     ,focean
      use lakes_com, only : flake
      use diag_com, only : npts,icon_wtg,icon_htg,conpt0
!      use sle001, only : w,ht,snowd, nsn,dzsn,wsn,hsn,fr_snow, dt,
!     &     hl0, set_snow
      use sle001, only : hl0, dt
#ifdef TRACERS_WATER
      use tracer_com, only : ntm,tr_wd_TYPE,nwater,itime_tr0,needtrs
      use fluxes, only : gtracer
      use veg_com, only:  avh !,afb
#endif
      use fluxes, only : gtemp,gtempr
      use ghy_com
      use dynamics, only : pedn
      use snow_drvm, only : snow_cover_coef2=>snow_cover_coef
     &     ,snow_cover_same_as_rad
#ifndef USE_ENT
      use veg_drv, only : init_vegetation
      use veg_com, only : vdata
      use veg_com, only : ala
#endif
      implicit none

      real*8, intent(in) :: dtsurf
      integer, intent(in) :: istart
      logical, intent(in) :: redogh, inisnow
      integer, intent(out) :: nl_soil
      integer iu_soil,iu_top_index
      integer jday
      real*8 snowdp,wtr1,wtr2,ace1,ace2,tg1,tg2
      logical :: qcon(npts)
      integer i, j, ibv
      logical ghy_data_missing
      character conpt(npts)*10
#ifdef TRACERS_WATER
      real*8 trsoil_tot,wsoil_tot,fm
      integer n
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
      integer I_0, I_1, J_0, J_1
      integer I_0H, I_1H, J_0H, J_1H
      logical present_land
      integer init_flake
      integer kk
      integer :: reset_canopy_ic=0, reset_snow_ic=0
      real*8 :: aa, ht_cap_can, fice_can, fb, fv

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

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
      call sync_param( "wsn_max", wsn_max )
      call sync_param( "ghy_default_data", ghy_default_data )
      call  get_param( "variable_lk", variable_lk )
      call  get_param( "init_flake", init_flake )
#ifdef USE_ENT
      call sync_param( "vegCO2X_off", vegCO2X_off)
#endif
      call sync_param( "reset_canopy_ic", reset_canopy_ic )
      call sync_param( "reset_snow_ic", reset_snow_ic )

c**** set number of layers for vegetation module
      nl_soil = ngm

c**** read land surface parameters or use defaults
      if ( ghy_default_data == 0 ) then ! read from files

        !!!if (istart.le.0) return ! avoid reading unneeded files
c**** read soils parameters
        call openunit("SOIL",iu_SOIL,.true.,.true.)
        ALLOCATE(TEMP_LOCAL(I_0H:I_1H,J_0H:J_1H,11*NGM+1))
        call DREAD_PARALLEL(grid,iu_SOIL,NAMEUNIT(iu_SOIL),TEMP_LOCAL)
        DZ_IJ(:,:,:)   = TEMP_LOCAL(:,:,1:NGM)
         Q_IJ(I_0:I_1,J_0:J_1,:,:) =
     &       RESHAPE( TEMP_LOCAL(I_0:I_1,J_0:J_1,1+NGM:) ,
     *                   (/I_1-I_0+1,J_1-J_0+1,imt,ngm/) )
        QK_IJ(I_0:I_1,J_0:J_1,:,:) =
     *        RESHAPE( TEMP_LOCAL(I_0:I_1,J_0:J_1,1+NGM+NGM*IMT:) ,
     *                   (/I_1-I_0+1,J_1-J_0+1,imt,ngm/) )
        SL_IJ(I_0:I_1,J_0:J_1)  =
     &       TEMP_LOCAL(I_0:I_1,J_0:J_1,1+NGM+NGM*IMT+NGM*IMT)
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
        if ( istart>0 .and. istart<8 ) then ! reset all
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

c**** cosday, sinday should be defined (reset once a day in daily_earth)
      jday=1+mod(itime/nday,365)
#ifndef USE_ENT
      cosday=cos(twopi/edpery*jday)
      sinday=sin(twopi/edpery*jday)
#endif

ccc read and initialize vegetation here
#ifndef USE_ENT
      call init_vegetation(redogh,istart)
#endif

      ! no need to continue computations for postprocessing
      if (istart.le.0) return

c**** check whether ground hydrology data exist at this point.
      ghy_data_missing = .false.
      do j=J_0,J_1
        do i=I_0,I_1
          present_land = .false.
          if (variable_lk==0) then
            if ( fearth(i,j) > 0.d0 ) present_land = .true.
          else
            if ( focean(i,j) < 1.d0 ) present_land = .true.
          endif
          if ( present_land ) then
            if ( top_index_ij(i,j).eq.-1. ) then
              print *,"No top_index data: i,j=",i,j,top_index_ij(i,j)
              ghy_data_missing = .true.
            end if
            if ( sum(dz_ij(i,j,1:ngm)).eq.0 ) then
              print *, "No soil data: i,j=",i,j,dz_ij(i,j,1:ngm)
              ghy_data_missing = .true.
            endif
            if (w_ij(1,1,i,j)<1.d-10 .and. w_ij(1,2,i,j)<1.d-10) then
              print*,"No gh data in restart file: i,j=",i,j,
     &             w_ij(:,1,i,j),w_ij(:,2,i,j)
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

c**** recompute ground hydrology data if necessary (new soils data)
      if (redogh) then

        do j=J_0,J_1
          do i=I_0,I_1
            w_ij(:,:,i,j)=0.d0
            ht_ij(:,:,i,j)=0.d0
            snowbv(:,i,j)=0.d0
            if ( focean(i,j) >= 1.d0 ) cycle
            if ( fearth(i,j) <= 0.d0 .and. variable_lk==0 ) cycle

            call old_gic_2_modele(
     &           w_ij(:,:,i,j), ht_ij(:,:,i,j),snowbv(:,i,j),
     &           wearth(i,j), aiearth(i,j), tearth(i,j), snowe(i,j) )

          end do
        end do
        write (*,*) 'ground hydrology data was made from ground data'
      end if

c**** set gtemp array
      do j=J_0,J_1
        do i=I_0,I_1
          if (fearth(i,j).gt.0) then
            gtemp(1,4,i,j)=tsns_ij(i,j)
            gtempr(4,i,j) =tearth(i,j)+tf
          end if
        end do
      end do

c**** set vegetation temperature to tearth(i,j) if requested
c**** in this case also set canopy water and water tracers to 0
      if ( reset_canopy_ic .ne. 0 .and. istart < 9 ) then
        do j=J_0,J_1
          do i=I_0,I_1
            present_land = .false.
            if (variable_lk==0) then
              if ( fearth(i,j) > 0.d0 ) present_land = .true.
            else
              if ( focean(i,j) < 1.d0 ) present_land = .true.
            endif

            if( .not. present_land ) cycle

            call get_fb_fv( fb, fv, i, j )
            if ( fv <= 0.d0 ) cycle

#ifdef USE_ENT
            !!! probably will not work
         call stop_model("reset_canopy_ic not implemented for Ent",255)
!            call ent_get_exports( entcells(i,j),
!     &           canopy_heat_capacity=ht_cap_can )
#else
            aa=ala(1,i,j)
            ht_cap_can=(.010d0+.002d0*aa+.001d0*aa**2)*shw_kg*rhow
#endif
            ! assume canopy completey dry
            w_ij(0,2,i,j) = 0.d0
            fice_can = 0.d0
            call temperature_to_heat( ht_ij(0,2,i,j),
     &           tsns_ij(i,j), fice_can, w_ij(0,2,i,j), ht_cap_can )
#ifdef TRACERS_WATER
            tr_w_ij(:,0,2,i,j) = 0.d0
#endif
          enddo
        enddo
      endif

C GISS-ESMF EXCEPTIONAL CASE
C-BMP Global sum on evap_max_ij

ccc if not initialized yet, set evap_max_ij, fr_sat_ij, qg_ij
ccc to something more appropriate

      call globalsum(grid, evap_max_ij,evap_max_ij_sum,
     &                all=.true.)
      if ( evap_max_ij_sum > im-1.d0 ) then ! old default
        do j=J_0,J_1
          do i=I_0,I_1
            if ( fearth(i,j) .le. 0.d0 ) cycle
            qg_ij(i,j) = qsat(tsns_ij(i,j)+tf,lhe,pedn(1,i,j))
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
          do i=I_0,I_1
            nsn_ij(:,i,j)     = 1
            dzsn_ij(:,:,i,j)  = 0.d0
            wsn_ij(:,:,i,j)   = 0.d0
            hsn_ij(:,:,i,j)   = 0.d0
            fr_snow_ij(:,i,j) = 0.d0
            if ( focean(i,j) >= 1.d0 ) cycle
            if ( fearth(i,j) <= 0.d0 .and. variable_lk==0 ) cycle

cddd            w(0:ngm,1:LS_NFRAC) = w_ij(0:ngm,1:LS_NFRAC,i,j)
cddd            ht(0:ngm,1:LS_NFRAC) = ht_ij(0:ngm,1:LS_NFRAC,i,j)
cddd            snowd(1:2) =  snowbv(1:2,i,j)
cddd
cddd            call ghinij (i,j)
cddd            call set_snow
cddd
cddd            nsn_ij    (1:2, i, j)         = nsn(1:2)
cddd            !isn_ij    (1:2, i, j)         = isn(1:2)
cddd            dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
cddd            wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
cddd            hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
cddd            fr_snow_ij(1:2, i, j)         = fr_snow(1:2)
cddd
cddd            w_ij (0:ngm,1:LS_NFRAC,i,j) = w (0:ngm,1:LS_NFRAC)
cddd            ht_ij(0:ngm,1:LS_NFRAC,i,j) = ht(0:ngm,1:LS_NFRAC)
cddd            snowbv(1:2,i,j)   = snowd(1:2)

            call set_snow1( w_ij (1:ngm,1:LS_NFRAC,i,j) ,
     &           ht_ij(1:ngm,1:LS_NFRAC,i,j) ,
     &           snowbv(1:2,i,j), q_ij(i,j,:,:), dz_ij(i,j,:),
     &           nsn_ij    (1:2, i, j),
     &           dzsn_ij   (1:nlsn, 1:2, i, j),
     &           wsn_ij    (1:nlsn, 1:2, i, j),
     &           hsn_ij    (1:nlsn, 1:2, i, j),
     &           fr_snow_ij(1:2, i, j) )


          end do
        end do
      end if


c**** remove all land snow from initial conditions
c**** (useful when changing land/vegetation mask)
      if ( reset_snow_ic .ne. 0 .and. istart < 9 ) then
        do j=J_0,J_1
          do i=I_0,I_1
            nsn_ij(:, i, j) = 1
            wsn_ij(:, :, i, j) = 0.d0
            hsn_ij(:, :, i, j) = 0.d0
            dzsn_ij(:, :, i, j) = 0.d0
            fr_snow_ij(:, i, j) = 0.d0
          enddo
        enddo
      endif

!!! hack - remove underwater snow + snow in landice boxes (dealt with separately)
!!! (should not be present in restart file in the first place!)
      do j=J_0,J_1
        do i=I_0,I_1
          if ( fearth(i,j) == 0.d0 ) then
!            if ( maxval(fr_snow_ij(:,i,j)) > 0.d0 ) 
!      &           print *,"removing snow from ",i,j,
!      &           " : cell under lake or land ice" 
            nsn_ij(:, i, j) = 1
            wsn_ij(:, :, i, j) = 0.d0
            hsn_ij(:, :, i, j) = 0.d0
            dzsn_ij(:, :, i, j) = 0.d0
            fr_snow_ij(:, i, j) = 0.d0
          endif
        enddo
      enddo

c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      fr_snow_rad_ij(:,:,:) = 0.d0
      do j=J_0,J_1
        do i=I_0,I_1
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

      ! initialize underwater fraction for variable lakes
      if ( init_flake > 0 .and. variable_lk > 0 .and. istart < 9 )
     &     call init_underwater_soil

      ! land water deficit for changing lake fractions
      !!! not working with Ent
#ifndef USE_ENT
      call compute_water_deficit
#endif

#ifdef TRACERS_WATER
ccc still not quite correct (assumes fw=1)
      do j=J_0,J_1
        do i=I_0,I_1
          gtracer(:,4,i,j)=0.  ! default
          if (fearth(i,j).le.0.d0) cycle
          !fb=afb(i,j) ; fv=1.-fb
          call get_fb_fv( fb, fv, i, j )
          fm=1.d0-exp(-snowbv(2,i,j)/((avh(i,j)*spgsn) + 1d-12))
          if ( fm < 1.d-3 ) fm=0.d0
          wsoil_tot=fb*( w_ij(1,1,i,j)*(1.d0-fr_snow_ij(1,i,j))
     &     + wsn_ij(1,1,i,j)*fr_snow_ij(1,i,j) )
     &     + fv*( w_ij(0,2,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))   !*1.d0
     &     + wsn_ij(1,2,i,j)*fm*fr_snow_ij(2,i,j) )
          do n=1,ntm
            if (itime_tr0(n).gt.itime) cycle
            if ( .not. needtrs(n) ) cycle
            ! should also restrict to TYPE=nWATER ?
            if ( wsoil_tot > 1.d-30 ) then
            gtracer(n,4,i,j) = (
     &           fb*( tr_w_ij(n,1,1,i,j)*(1.d0-fr_snow_ij(1,i,j))
     &           + tr_wsn_ij(n,1,1,i,j) )         !*fr_snow_ij(1,i,j)
     &           + fv*( tr_w_ij(n,0,2,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))
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


      subroutine old_gic_2_modele(
     &           w, ht,snowbv,
     &           wearth, aiearth, tearth, snowe )
      real*8 :: w(:,:), ht(:,:),snowbv(:),
     &           wearth, aiearth, tearth, snowe

      call stop_model(
     &     "conversion old_gic_2_modele not supported yet",255)

      end subroutine old_gic_2_modele


      subroutine set_snow1( w, ht, snowd, q, dz,
     &     nsn, dzsn, wsn, hsn, fr_snow )
!@sum set_snow extracts snow from the first soil layer and initializes
!@+   snow model prognostic variables
!@+   should be called when model restarts from the old restart file
!@+   ( which doesn''t contain new snow model (i.e. 3 layer) data )
c
c input:
c snowd(2) - landsurface snow depth
c w(k,2)   - landsurface water in soil layers
c ht(k,2)  - landsurface heat in soil layers
c fsn      - heat of fusion
c shi      - specific heat of ice
c shc(k,2) - heat capacity of soil layers
c
c output:
c dzsn(lsn,2) - snow layer thicknesses
c wsn(lsn,2)  - snow layer water equivalent depths
c hsn(lsn,2)  - snow layer heat contents
c tsn1(2)     - snow top temperature
c isn(2)      - 0 if no snow, 1 if snow
c nsn(2)      - number of snow layers
c snowd(2)
c w(k,2)
c ht(k,2)
c
c calling sequence:
c
c     assignment of w,ht,snowd
c     call ghinij(i,j,wfc1)
c     call set_snow
c note: only to be called when initializing from landsurface
c       prognostic variables without the snow model.
c
      use constant, only : rhow
     &     ,shi_kg=>shi,lhm
      use snow_model, only: snow_redistr, snow_fraction
      use ghy_com, only : ngm
      use sle001, only : get_soil_properties
      real*8, intent(inout) :: w(:,:), ht(:,:), snowd(:), q(:,:), dz(:)
      integer, intent(out) :: nsn(:)
      real*8, intent(out) :: dzsn(:,:), wsn(:,:), hsn(:,:), fr_snow(:)
      !---
!@var shi heat capacity of pure ice (J/m^3 C)
      real*8, parameter :: shi= shi_kg * rhow
!@var fsn latent heat of melt (J/m^3)
      real*8, parameter :: fsn= lhm * rhow
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2), tsn1(2), fice(2)


      integer ibv

c outer loop over ibv
      do ibv=1,2

        call get_soil_properties( q, dz,
     &     thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )

c initalize all cases to nsn=1
        nsn(ibv)=1


ccc since we don''t know what kind of data we are dealing with,
ccc better check it

        if( snowd(ibv) .gt. w(1,ibv)-dz(1)*thetm(1,ibv)  ) then
          write(99,*) 'snowd corrected: old=', snowd(ibv)
          snowd(ibv) = w(1,ibv)-dz(1)*thetm(1,ibv) - 1.d-10
          write(99,*) '                 new=', snowd(ibv)
          if ( snowd(ibv) .lt. -0.001d0 )
     &         call stop_model('set_snow: neg. snow',255)
          if ( snowd(ibv) .lt. 0.d0 ) snowd(ibv) = 0.d0 ! rounding error
        endif


c if there is no snow, set isn=0.  set snow variables to 0.d0
        if(snowd(ibv).le.0.d0)then
          dzsn(1,ibv)=0.d0
          wsn(1,ibv)=0.d0
          hsn(1,ibv)=0.d0
          tsn1(ibv)=0.d0
          fr_snow(ibv) = 0.d0
        else


c given snow, set isn=1.d0
c!!!        dzsn(1,ibv)=snowd(ibv)/spgsn
c!!!  replacing prev line considering rho_snow = 200
          dzsn(1,ibv)=snowd(ibv) * 5.d0
          wsn(1,ibv)=snowd(ibv)
c!!! actually have to compute fr_snow and modify dzsn ...
          fr_snow(ibv) = 1.d0

c given snow, temperature of first layer can''t be positive.
c the top snow temperature is the temperatre of the first layer.
cddd          if(fsn*w(1,ibv)+ht(1,ibv).lt.0.d0)then
cddd            tsn1(ibv)=(ht(1,ibv)+w(1,ibv)*fsn)/(shc(1,ibv)+w(1,ibv)*shi)
cddd          else
cddd            tsn1(ibv)=0.d0
cddd          endif
cddd          print *,"tsn = ", tsn1(ibv)
          call heat_to_temperature(tsn1(ibv),
     &         fice(ibv), ht(1,ibv), w(1,ibv), shc(1,ibv))
cddd          print *,"      ", tsn1(ibv)
          tsn1(ibv) = min( tsn1(ibv), 0.d0 )

c use snow temperature to get the heat of the snow
cddd          hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
cddd          print *,"hsn = ", hsn(1,ibv)
          call temperature_to_heat(hsn(1,ibv),
     &         tsn1(ibv), 1.d0, wsn(1,ibv), 0.d0)
cddd          print *,"      ", hsn(1,ibv)

          ! hack to get identical results
          hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn

c subtract the snow from the landsurface prognositic variables
          w(1,ibv)=w(1,ibv)-wsn(1,ibv)
          ht(1,ibv)=ht(1,ibv)-hsn(1,ibv)

ccc and now limit all the snow to 5cm water equivalent
          if ( snowd(ibv) .gt. 0.05d0 ) then
            snowd(ibv) = 0.05d0
            dzsn(1,ibv)= snowd(ibv) * 5.d0
            wsn(1,ibv)= snowd(ibv)
            hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
          endif

ccc and now limit all the snow to 50cm water equivalent (for debug)
cddd          if ( snowd(ibv) .gt. 0.0005d0 ) then
cddd            snowd(ibv) = 0.5d0
cddd            dzsn(1,ibv)= snowd(ibv) * 5.d0
cddd            wsn(1,ibv)= snowd(ibv)
cddd            hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
cddd          endif

ccc redistribute snow over the layers and recompute fr_snow
ccc (to make the data compatible with snow model)
          if ( .not. ( hsn(1,ibv) > 0. .or.  hsn(1,ibv) <= 0. ) )
     &        call stop_model("ERR in init_snow: NaN", 255)

          call snow_fraction(dzsn(:,ibv), nsn(ibv), 0.d0, 0.d0,
     &         1.d0, fr_snow(ibv) )
          if ( fr_snow(ibv) > 0.d0 ) then
            call snow_redistr(dzsn(:,ibv), wsn(:,ibv), hsn(:,ibv),
     &           nsn(ibv), 1.d0/fr_snow(ibv) )
            if ( .not. ( hsn(1,ibv) > 0. .or.  hsn(1,ibv) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)
          else
            snowd(ibv) = 0.d0
            dzsn(1,ibv)=0.d0
            wsn(1,ibv)=0.d0
            hsn(1,ibv)=0.d0
            tsn1(ibv)=0.d0
            fr_snow(ibv) = 0.d0
          endif

        endif
!!! debug : check if  snow is redistributed correctly
        if ( dzsn(1,ibv) > 0.d0 .and. dzsn(1,ibv) < .099d0) then
          call stop_model("set_snow: error in dz",255)
        endif

        if ( dzsn(1,ibv) > 0.d0 ) then
          if (  wsn(1,ibv)/dzsn(1,ibv)  < .1d0)
     &         call stop_model("set_snow: error in dz",255)
        endif


      enddo

            if ( .not. ( hsn(1,1) > 0. .or.  hsn(1,1) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)
            if ( .not. ( hsn(1,2) > 0. .or.  hsn(1,2) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)



      return
      end subroutine set_snow1



      subroutine reset_gh_to_defaults( reset_prognostic )
      !use model_com, only: vdata
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use ghy_com
#ifndef USE_ENT
      use veg_drv, only : reset_veg_to_defaults
#endif
      logical, intent(in) :: reset_prognostic
      integer i,j

C**** define local grid
      integer I_0, I_1, J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

ccc ugly, should fix later
#ifdef USE_ENT
      call stop_model(
     &     "reset_gh_to_defaults not implemented yet for Ent",255)
#else
      call reset_veg_to_defaults( reset_prognostic )
#endif

      do j=J_0,J_1
      do i=I_0,I_1

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
      tsns_ij(i,j)= -0.12476520d+00
      wearth(i,j)=  0.29203081d+02
      aiearth(i,j)=  0.93720329d-01
      w_ij(:,1,i,j) = (/0.d0,  0.17837750d-01,  0.40924843d-01,
     &     0.77932012d-01,  0.11919649d+00,  0.57237469d-01,
     &     0.10000000d-11 /)
      w_ij(:,2,i,j) = (/  0.10000000d-11,  0.29362259d-01,
     &     0.50065177d-01,  0.82533140d-01,  0.10383620d+00,
     &     0.31552459d-01,  0.10000000d-11 /)
      ht_ij(:,1,i,j)= (/  0.00000000d+00, -0.15487181d+07,
     &     -0.50720067d+07,  0.18917623d+07,  0.77174974d+07,
     &     0.21716040d+08,  0.44723067d+08 /)
      ht_ij(:,2,i,j)= (/ -0.13991376d+05, -0.53165599d+05,
     &     0.65443775d+06,  0.29276050d+07,  0.81455096d+07,
     &     0.21575081d+08,  0.45952255d+08 /)
      snowbv(1:2,i,j)= (/  0.00000000d+00,  0.65458111d-04 /)

      enddo
      enddo

      end subroutine reset_gh_to_defaults


cddd      subroutine ghinij (i0,j0)
cdddc**** input:
cdddc**** avh(i,j) - array of vegetation heights
cdddc**** spgsn - specific gravity of snow
cdddc**** output:
cdddc**** vh - vegetation height
cdddc**** snowm - snow masking depth
cdddc**** wfcap - water field capacity of top soil layer, m
cdddc****
cddd      use snow_model, only : i_earth,j_earth
cddd      use sle001, only : dz,qk,ng,zb,zc,q,sl,xklh !spgsn,
cddd     *     ,fb,fv,prs,ijdebug,n
cddd     *     ,thets,thetm,ws,thm,nth,shc,shw,htprs,pr !shcap,shtpr,
cddd     *     ,htpr
cddd     *     ,top_index,top_stdev
cddd     &     ,w,ht,snowd,nsn,dzsn,wsn,hsn,fr_snow
cddd#ifdef USE_ENT
cddd     &     ,Ci, Qf, cnc
cddd#endif
cddd      use ghy_com, only : ngm,imt,nlsn,LS_NFRAC,dz_ij,sl_ij,q_ij,qk_ij
cddd     *     ,top_index_ij,top_dev_ij
cddd     &     ,w_ij,ht_ij,snowbv,nsn_ij,dzsn_ij,wsn_ij
cddd     &     ,hsn_ij,fr_snow_ij,shc_soil_texture
cddd#ifdef USE_ENT
cddd     &     ,Ci_ij, Qf_ij, cnc_ij
cddd#endif
cddd#ifndef USE_ENT
cddd      use veg_com, only: afb
cddd#endif
cddd      USE DOMAIN_DECOMP, ONLY : GRID, GET
cddd!      use veg_drv, only : veg_set_cell
cddd
cddd      implicit none
cddd      integer i0,j0
cddd!      real*8 wfcap
cddd      integer k,ibv,i
cddd      real*8 shtpr
cddd!----------------------------------------------------------------------!
cddd      !real*8, parameter :: shcap(imt) = (/2d6,2d6,2d6,2.5d6,2.4d6/)
cddd
cddd
cddd      ijdebug=i0*1000+j0
cddd      i_earth = i0
cddd      j_earth = j0
cddd
cdddccc extracting ghy prognostic vars
cddd      !w(1:ngm,1) =  wbare(1:ngm,i0,j0)
cddd      !w(0:ngm,2) =  wvege(0:ngm,i0,j0)
cddd      !ht(0:ngm,1) = htbare(0:ngm,i0,j0)
cddd      !ht(0:ngm,2) = htvege(0:ngm,i0,j0)
cddd      w(0:ngm,1:LS_NFRAC) = w_ij(0:ngm,1:LS_NFRAC,i0,j0)
cddd      ht(0:ngm,1:LS_NFRAC) = ht_ij(0:ngm,1:LS_NFRAC,i0,j0)
cddd      snowd(1:2)  = snowbv(1:2,i0,j0)
cdddccc extracting snow variables
cddd      nsn(1:2)          = nsn_ij    (1:2, i0, j0)
cddd      !isn(1:2)          = isn_ij    (1:2, i0, j)
cddd      dzsn(1:nlsn, 1:2) = dzsn_ij   (1:nlsn, 1:2, i0, j0)
cddd      wsn(1:nlsn, 1:2)  = wsn_ij    (1:nlsn, 1:2, i0, j0)
cddd      hsn(1:nlsn, 1:2)  = hsn_ij    (1:nlsn, 1:2, i0, j0)
cddd      fr_snow(1:2)      = fr_snow_ij(1:2, i0, j0)
cddd#ifdef USE_ENT
cdddccc extracting vegetation prognostic variables
cddd      cnc = cnc_ij(i0,j0)
cddd      Ci = Ci_ij(i0,j0)
cddd      Qf = Qf_ij(i0,j0)
cddd#endif
cdddccc setting vegetation
cddd !     call veg_set_cell(i0,j0)
cddd
cdddccc passing topmodel parameters
cddd      top_index = top_index_ij(i0, j0)
cddd      top_stdev = top_dev_ij(i0, j0)
cdddc**** set up layers
cddd      dz(1:ngm)=dz_ij(i0,j0,1:ngm)
cddd      q(1:imt,1:ngm)=q_ij(i0,j0,1:imt,1:ngm)
cddd      qk(1:imt,1:ngm)=qk_ij(i0,j0,1:imt,1:ngm)
cddd      sl=sl_ij(i0,j0)
cddd
cddd      n=0
cddd      do k=1,ngm
cddd        if(dz(k).le.0.) exit
cddd        n=k
cddd      end do
cddd
cddd      if(n.le.0) then
cddd         write (99,*) 'ghinij:  n <= 0:  i,j,n=',i0,j0,n,(dz(k),k=1,43)
cddd         call stop_model('stopped in GHY_DRV.f',255)
cddd      end if
cddd
cdddc**** calculate the boundaries, based on the thicknesses.
cddd      zb(1)=0.
cddd      do k=1,n
cddd        zb(k+1)=zb(k)-dz(k)
cddd      end do
cdddc**** calculate the layer centers, based on the boundaries.
cddd      do k=1,n
cddd        zc(k)=.5*(zb(k)+zb(k+1))
cddd      end do
cdddc**** fb,fv: bare, vegetated fraction (1=fb+fv)
cddd!!! I think we still need this
cddd      fb=afb(i0,j0)
cddd      fv=1.-fb
cddd      call get_fb_fv( fb, fv, i0, j0 )
cdddc****
cddd      do ibv=1,2
cddd        do k=1,n
cddd          thets(k,ibv)=0.
cddd          thetm(k,ibv)=0.
cddd          do i=1,imt-1
cddd            thets(k,ibv)=thets(k,ibv)+q(i,k)*thm(0,i)
cddd            thetm(k,ibv)=thetm(k,ibv)+q(i,k)*thm(nth,i)
cddd          end do
cddd          ws(k,ibv)=thets(k,ibv)*dz(k)
cddd        end do
cddd      end do
cddd!veg      ws(0,2)=.0001d0*alai  ! from call veg_set_cell above
cddd  !    wfcap=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
cdddc****
cddd      call xklh(1)
cdddc****
cddd
cddd      do ibv=1,2
cddd        do k=1,n
cddd          shc(k,ibv)=0.
cddd          do i=1,imt
cddd            shc(k,ibv)=shc(k,ibv)+q(i,k)*shc_soil_texture(i)
cddd          end do
cddd          shc(k,ibv)=(1.-thets(k,ibv))*shc(k,ibv)*dz(k)
cddd        end do
cddd      end do
cdddc****
cdddc shc(0,2) is the heat capacity of the canopy
cddd!veg      aa=ala(1,i0,j0)
cddd!veg      shc(0,2)=(.010d0+.002d0*aa+.001d0*aa**2)*shw
cdddc****
cdddc htpr is the heat of precipitation.
cdddc shtpr is the specific heat of precipitation.
cddd      shtpr=0.
cddd      if(pr.gt.0.)shtpr=htpr/pr
cdddc htprs is the heat of large scale precipitation
cddd      htprs=shtpr*prs
cdddc****
cddd      return
cddd      end subroutine ghinij


cddd      subroutine ghy_save_cell(i,j)
cddd      use sle001, only : w,ht,snowd,nsn,dzsn,wsn,hsn,fr_snow
cddd#ifdef USE_ENT
cddd     $     ,cnc,Ci,Qf
cddd#endif
cddd      use ghy_com, only : ngm,nlsn,LS_NFRAC
cddd     &     ,dz_ij,w_ij,ht_ij,snowbv
cddd     &     ,nsn_ij,dzsn_ij,wsn_ij,hsn_ij,fr_snow_ij
cddd#ifdef USE_ENT
cddd     $     ,cnc_ij,Ci_ij,Qf_ij
cddd#endif
cddd      implicit none
cddd      integer, intent(in) :: i,j
cddd
cddd      !wbare(1:ngm,i,j) = w(1:ngm,1)
cddd      !wvege(0:ngm,i,j) = w(0:ngm,2)
cddd      !htbare(0:ngm,i,j) = ht(0:ngm,1)
cddd      !htvege(0:ngm,i,j) = ht(0:ngm,2)
cddd      w_ij (0:ngm,1:LS_NFRAC,i,j) = w (0:ngm,1:LS_NFRAC)
cddd      ht_ij(0:ngm,1:LS_NFRAC,i,j) = ht(0:ngm,1:LS_NFRAC)
cddd      snowbv(1:2,i,j)   = snowd(1:2)
cdddccc copy snow variables back to storage
cddd      nsn_ij    (1:2, i, j)         = nsn(1:2)
cddd      !isn_ij    (1:2, i, j)         = isn(1:2)
cddd      dzsn_ij   (1:nlsn, 1:2, i, j) = dzsn(1:nlsn,1:2)
cddd      wsn_ij    (1:nlsn, 1:2, i, j) = wsn(1:nlsn,1:2)
cddd      hsn_ij    (1:nlsn, 1:2, i, j) = hsn(1:nlsn,1:2)
cddd      fr_snow_ij(1:2, i, j)         = fr_snow(1:2)
cddd#ifdef USE_ENT
cdddccc saving vegetation prognostic variables
cddd      cnc_ij(i,j) = cnc
cddd      Ci_ij(i,j) = Ci
cddd      Qf_ij(i,j) = Qf
cddd#endif
cddd      end subroutine ghy_save_cell


cddd      subroutine ghinht (snowdp,tg1,tg2,wtr1,wtr2,ace1,ace2)
cdddc**** initializes new ground (w,ht,snw) from old (t,w,ice,snw)
cdddc**** evaluates the heat in the soil layers based on the
cdddc**** temperatures.
cdddc**** input:
cdddc**** w - water in soil layers, m
cdddc**** tp - temperature of layers, c
cdddc**** fice - fraction of ice of layers
cdddc**** fsn - heat of fusion of water
cdddc**** shc - specific heat capacity of soil
cdddc**** shi - specific heat capacity of ice
cdddc**** shw - specific heat capcity of water
cdddc**** snowd - snow depth, equivalent water m
cdddc**** output:
cdddc**** ht - heat in soil layers
cdddc**** add calculation of wfc2
cdddc**** based on combination of layers 2-n, as in retp2
cddd      use sle001, only : tp, ht, w, shc, fice, snowd, ws, fb, fv,
cddd     &    n, dz, fsn, thetm, shi, shw, ijdebug
cddd      USE DOMAIN_DECOMP, ONLY : GRID, GET
cddd      implicit none
cddd
cddd      real*8 snowdp,tg1,tg2,wtr1,wtr2,ace1,ace2
cddd      real*8 wfc1, wfc2, wet1, wet2, wmin, fbv
cddd      integer k, ibv, ll
cddd
cddd      wfc1=fb*ws(1,1)+fv*(ws(0,2)+ws(1,2))
cddd      wfc2=0.
cddd      fbv=fb
cddd      do 30 ibv=1,2
cddd      do 20 k=2,n
cddd      wfc2=wfc2+fbv*ws(k,ibv)
cddd   20 continue
cddd      fbv=fv
cddd   30 continue
cddd      wfc1=1000.*wfc1
cddd      wfc2=1000.*wfc2
cddd      fice(0,2)=1.
cddd      fice(1,1)=(ace1+snowdp*1000.)/(wtr1+ace1+snowdp*1000.+1.d-20)
cddd      fice(1,2)=fice(1,1)
cddd      tp(0,2)=tg1
cdddc**** w = snow(if top layer) + wmin + (wmax-wmin)*(wtr+ice)/wfc
cddd      w(0,2)=0.
cddd      do ibv=1,2
cddd        w(1,ibv)=snowdp
cddd        wmin=thetm(1,ibv)*dz(1)
cddd        wet1=(wtr1+ace1)/(wfc1+1.d-20)
cddd        if(wet1.gt.1.) wet1=1.
cddd        w(1,ibv)=w(1,ibv)+wmin+(ws(1,ibv)-wmin)*wet1
cddd        snowd(ibv)=snowdp
cddd        tp(1,ibv)=tg1
cddd        do k=2,n
cddd          fice(k,ibv)=ace2/(wtr2+ace2+1.d-20)
cddd          wmin=thetm(k,ibv)*dz(k)
cddd          wet2=(wtr2+ace2)/(wfc2+1.d-20)
cddd          if(wet2.gt.1.) wet2=1.
cddd          w(k,ibv)=wmin+(ws(k,ibv)-wmin)*wet2
cddd          tp(k,ibv)=tg2
cddd        end do
cddd      end do
cdddc****
cddd!!!      entry ghexht
cdddc****
cdddc**** compute ht (heat w/m+2)
cddd      do ibv=1,2
cddd        ll=2-ibv
cddd        do k=ll,n
cddd          if(tp(k,ibv)) 2,4,6
cddd 2        ht(k,ibv)=tp(k,ibv)*(shc(k,ibv)+w(k,ibv)*shi)-w(k,ibv)*fsn
cddd          cycle
cddd 4        ht(k,ibv)=-fice(k,ibv)*w(k,ibv)*fsn
cddd          cycle
cddd 6        ht(k,ibv)=tp(k,ibv)*(shc(k,ibv)+w(k,ibv)*shw)
cddd        end do
cddd      end do
cddd      if(ijdebug.eq.0)then
cddd       write(99,*)'ghinht id check',ijdebug
cddd       write(99,*)'tg1,tg2',tg1,tg2
cddd       write(99,*)'tp',tp
cddd       write(99,*)'ht',ht
cddd       write(99,*)'w',w
cddd       write(99,*)'wtr1,wtr2',wtr1,wtr2
cddd       write(99,*)'ace1,ace2',ace1,ace2
cddd       write(99,*)'wfc1,wfc2',wfc1,wfc2
cddd       write(99,*)'shc',shc
cddd       write(99,*)'fice',fice
cddd      endif
cddd      return
cddd      end subroutine ghinht

cddd      subroutine retp2 (tg2av,wtr2av,ace2av)
cdddc**** evaluates the mean temperature in the soil layers 2-ngm
cdddc**** as well as the water and ice content.
cdddc**** input:
cdddc**** w - water in soil layers, m
cdddc**** ht - heat in soil layers
cdddc**** fsn - heat of fusion of water
cdddc**** shc - specific heat capacity of soil
cdddc**** shi - specific heat capacity of ice
cdddc**** shw - specific heat capcity of water
cdddc**** output:
cdddc**** tg2av - temperature of layers 2 to ngm, c
cdddc**** ice2av - ice amount in layers 2 to ngm, kg/m+2
cdddc**** wtr2av - water in layers 2 to ngm, kg/m+2
cddd      USE DOMAIN_DECOMP, ONLY : GRID, GET
cddd      use sle001
cddd      implicit none
cddd      real*8 tg2av,wtr2av,ace2av, wc,htc,shcc,tpc,ficec,ftp
cddd      integer k, ibv
cddd      tg2av=0.
cddd      wtr2av=0.
cddd      ace2av=0.
cddd      do 3500 ibv=1,2
cddd      wc=0.
cddd      htc=0.
cddd      shcc=0.
cddd      do k=2,n
cddd        wc=wc+w(k,ibv)
cddd        htc=htc+ht(k,ibv)
cddd        shcc=shcc+shc(k,ibv)
cddd      end do
cddd      tpc=0.
cddd      ficec=0.
cddd      if(wc.ne.0.)  ficec=-htc/(fsn*wc)
cddd      if(fsn*wc+htc.ge.0.)go to 3430
cddd      tpc=(htc+wc*fsn)/(shcc+wc*shi)
cddd      ficec=1.
cddd      go to 3440
cddd 3430 if(htc.le.0.) go to 3440
cddd      tpc=htc/(shcc+wc*shw)
cddd      ficec=0.
cddd 3440 continue
cddd      ftp=fb
cddd      if(ibv.eq.2) ftp=fv
cddd      tg2av=tg2av+tpc*ftp
cddd      wtr2av=wtr2av+wc*ftp*1000.*(1.-ficec)
cddd      ace2av=ace2av+wc*ftp*1000.*ficec
cddd 3500 continue
cddd      return
cddd      end subroutine retp2

      subroutine checke(subr)
!@sum  checke checks whether arrays are reasonable over earth
!@auth original development team
!@ver  1.0
      use model_com, only : itime,wfcs
      use geom, only : imaxj
      use constant, only : rhow
      !use veg_com, only : afb
      use ghy_com, only : tearth,wearth,aiearth,snowe,w_ij,ht_ij
     *     ,snowbv,ngm,fearth,wsn_ij,fr_snow_ij,nsn_ij,LS_NFRAC
#ifdef TRACERS_WATER
     &     ,tr_w_ij,tr_wsn_ij
      USE TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      implicit none
!@var subr identifies where check was called from
      character*6, intent(in) :: subr

      real*8 x,tgl,wtrl,acel
      integer i,j,imax,jmax,n,nsb,nsv
      real*8, parameter :: EPS=1.d-12
      logical QCHECKL
      real*8 relerr, errmax, fb, fv

C**** define local grid
      integer I_0, I_1, J_0, J_1, njpol

      QCHECKL = .false.
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)
      njpol = grid%J_STRT_SKP-grid%J_STRT

c**** check for nan/inf in earth data
      call check3c(w_ij(1:ngm,1,I_0:I_1,J_0:J_1) ,ngm  ,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'wb')
      call check3c(w_ij(0:ngm,2,I_0:I_1,J_0:J_1) ,ngm+1,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'wv')
      call check3c(ht_ij(0:ngm,1,I_0:I_1,J_0:J_1),ngm+1,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'hb')
      call check3c(ht_ij(0:ngm,2,I_0:I_1,J_0:J_1),ngm+1,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'hv')
      call check3c(snowbv(1:LS_nfrac,I_0:I_1,J_0:J_1),2    ,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'sn')

c**** check for reasonable temperatures over earth
      x=1.001
      do j=J_0,J_1
        do i=I_0,imaxj(j)
          if (fearth(i,j).gt.0.) then
            tgl=tearth(i,j)
            wtrl=wearth(i,j)
            acel=aiearth(i,j)
            if ((tgl+60.)*(60.-tgl).le.0.) write (6,901) subr,i,j,itime
     *           ,fearth(i,j),'tg1 off',snowe(i,j),tgl,wtrl,acel
            if (wtrl.lt.0..or.acel.lt.0..or.(wtrl+acel).gt.x*wfcs(i
     *           ,j)) write(6,901) subr,i,j,itime,fearth(i,j),'wtr off'
     *           ,snowe(i,j),tgl,wtrl,acel,wfcs(i,j)
          end if
        end do
      end do

      !print *,"checke: w(51,33) ", w_ij(:,:,51,33)
      !print *,"checke: fearth(51,33) ", fearth(51,33),afb(51,33)


c**** check tracers
#ifdef TRACERS_WATER
      do n=1,ntm
        ! check for neg tracers
        if (t_qlimit(n)) then
          do j=J_0, J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              if ( minval( tr_w_ij(n,:,:,i,j)   ) < -1.d15 .or.
     &             minval( tr_wsn_ij(n,:,:,i,j) ) < -1.d15 ) then
                print*,"Neg tracer in earth after ",SUBR,i,j,trname(n)
     &               , "tr_soil= ", tr_w_ij(n,:,:,i,j)
     &               , "TR_SNOW= ", tr_wsn_ij(n,:,:,i,j)
                QCHECKL=.TRUE.
              end if

            end do
          end do
        end if

        ! check if water == water
        if (trname(n) == 'Water') then
          errmax = 0. ; imax=1 ; jmax=1
          do j=J_0, J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              !fb = afb(i,j)
              !fv = 1.d0 - fb
              call get_fb_fv( fb, fv, i, j )
              relerr= ( fb*sum(abs(
     &             tr_w_ij(n,1:ngm,1,i,j) - w_ij(1:ngm,1,i,j)*rhow))
     &             + fv*sum(abs(
     &             tr_w_ij(n,0:ngm,2,i,j) - w_ij(0:ngm,2,i,j)*rhow)) )
     &             /(  rhow*( fb*sum(w_ij(1:ngm,1,i,j))
     &             + fv*sum(w_ij(0:ngm,2,i,j)) ) + EPS  )
              if (relerr > errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            enddo
          enddo
          !fb = afb(imax,jmax)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, imax, jmax )
          print*,"Relative error in soil after ",trim(subr),":"
     &         ,imax,jmax,errmax
     &         ,( tr_w_ij(n,1:ngm,1,imax,jmax)
     &         - w_ij(1:ngm,1,imax,jmax)*rhow )*fb
     &         ,( tr_w_ij(n,0:ngm,2,imax,jmax)
     &         - w_ij(0:ngm,2,imax,jmax)*rhow )*fv
     &         , rhow*( fb*sum(w_ij(1:ngm,1,imax,jmax))
     &         +        fv*sum(w_ij(0:ngm,2,imax,jmax)) )
cddd     &         ,w_ij(1:ngm,1,imax,jmax)*rhow
cddd     &         ,w_ij(0:ngm,2,imax,jmax)*rhow

          errmax = 0. ; imax=1 ; jmax=1
          do j=J_0, J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              !fb = afb(i,j)
              !fv = 1.d0 - fb
              call get_fb_fv( fb, fv, i, j )
              nsb = nsn_ij(1,i,j)
              nsv = nsn_ij(2,i,j)
              relerr= ( fb*sum(abs( tr_wsn_ij(n,1:nsb,1,i,j)
     &             - wsn_ij(1:nsb,1,i,j)*rhow*fr_snow_ij(1,i,j) ))
     &             + fv*sum(abs( tr_wsn_ij(n,1:nsv,2,i,j)
     &             - wsn_ij(1:nsv,2,i,j)*rhow*fr_snow_ij(2,i,j))) )
     &             /(rhow*(
     &             fb*sum(wsn_ij(1:nsb,1,i,j))*fr_snow_ij(1,i,j) +
     &             fv*sum(wsn_ij(1:nsv,2,i,j))*fr_snow_ij(2,i,j) ) +EPS)
              if (relerr > errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            enddo
          enddo
          !fb = afb(imax,jmax)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, imax, jmax )
          nsb = nsn_ij(1,imax,jmax)
          nsv = nsn_ij(2,imax,jmax)
          print*,"Relative error in snow after ",trim(subr),":"
     &         ,imax,jmax,errmax
     &         ,( tr_wsn_ij(n,1:nsb,1,imax,jmax)
     &         - wsn_ij(1:nsb,1,imax,jmax)*rhow
     &         *fr_snow_ij(1,imax,jmax) )*fb
     &         ,( tr_wsn_ij(n,1:nsv,2,imax,jmax)
     &         - wsn_ij(1:nsv,2,imax,jmax)*rhow
     &         *fr_snow_ij(2,imax,jmax) )*fv
     &         ,rhow*( fb*sum(wsn_ij(1:nsb,1,imax,jmax))
     &         * fr_snow_ij(1,imax,jmax)
     &         +       fv*sum(wsn_ij(1:nsv,2,imax,jmax))
     &         * fr_snow_ij(2,imax,jmax) )

cddd     &         ,wsn_ij(1:nsn_ij(1,imax,jmax),1,imax,jmax)*rhow
cddd     &         *fr_snow_ij(1,imax,jmax)
cddd     &         ,wsn_ij(1:nsn_ij(1,imax,jmax),2,imax,jmax)*rhow
cddd     &         *fr_snow_ij(2,imax,jmax)



        endif
      enddo
#endif

      IF (QCHECKL)
     &     call stop_model('CHECKL: Earth variables out of bounds',255)

      return
 901  format ('0subr,i,j,i-time,pearth,',a7,2i4,i10,f5.2,1x
     *     ,a/' snw,tg1,wtr1,ice1 (wfc1)',5f12.4)

      end subroutine checke

      subroutine daily_earth(end_of_day)
!@sum  daily_earth performs daily tasks for earth related functions
!@auth original development team
!@ver  1.0
!@calls RDLAI
      use constant, only : rhow,twopi,edpery,tf
      use model_com, only : nday,nisurf,jday,jyear,wfcs,focean
#ifndef USE_ENT
      use veg_com, only : vdata                 !nyk
#endif
      use geom, only : imaxj,lat2d
      use diag_com, only : aij=>aij_loc
     *     ,tdiurn,ij_strngts,ij_dtgdts,ij_tmaxe
     *     ,ij_tdsl,ij_tmnmx,ij_tdcomp, ij_dleaf
      use ghy_com, only : snoage, snoage_def,fearth, wsn_max,
     &     q_ij,dz_ij,ngm
#ifdef USE_ENT
     &     ,aalbveg
#else
      use veg_com, only : almass,aalbveg       !nyk
      use vegetation, only: crops_yr,cond_scheme,vegCO2X_off !nyk
#endif
      use surf_albedo, only: albvnh, updsur  !nyk
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      !use sle001, only : fb,fv,ws
      use sle001, only : get_soil_properties
#ifdef USE_ENT
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
#else
      use veg_drv, only : veg_set_cell
      use vegetation, only : t_vegcell
#endif
      !!use ent_com, only : entcells
      !!use ent_mod, only : ent_get_exports

      implicit none
      real*8 tsavg,wfc1
      real*8 aleafmass, aalbveg0, fvp, sfv  !nyk veg ! , aleafmasslast
      integer i,j,itype
      integer northsouth,iv  !nyk
      logical, intent(in) :: end_of_day
      real*8 fb, fv, ws_can
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2)
      integer ibv
#ifdef USE_ENT
!      integer hemi(1:IM,grid%J_STRT:grid%J_STOP)
!      integer :: JEQUATOR=JM/2
#else
      type(t_vegcell) :: vegcell
#endif
C**** define local grid
      integer I_0, I_1, J_0, J_1
      real*8 ws11,ws12

C**** Extract useful local domain parameters from "grid"
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      ! update water and heat in land surface fractions if
      ! lake fraction changed
      if (end_of_day .and. variable_lk > 0) call update_land_fractions

      if (end_of_day .and. wsn_max>0) call remove_extra_snow_to_ocean

!!! testing
!!!      aalbveg(:,:) = 0.08D0
!!!      return

#ifdef USE_ENT
C**** Update vegetation file if necessary  (i.e. if crops_yr=0)
      ! if(crops_yr.eq.0) - should be checked inside Ent
      !call ent_update_crops(jyear)

! the following is probably not needed since we update vegetation
! on every time step
 !     if (jday==1) then ! i guess we need to call it only once per year
 !       hemi(:,JEQUATOR+1:J_1) = 1
 !       hemi(:,J_0:JEQUATOR) = -1
 !       call ent_prescribe_vegupdate(entcells,hemi,jday,jyear, .true.)
 !     endif
      !if(cond_scheme.eq.2) call updsur (0,jday)
      ! we don't use cond_scheme==1 any more, so call it always
      call updsur (0,jday)
c****
#else
      if(crops_yr.eq.0) call updveg(jyear,.true.)
c**** find leaf-area index & water field capacity for ground layer 1
      if(cond_scheme.eq.2) call updsur (0,jday) ! Update vegn albedos
#endif
            !albvnh(9,6,2)=albvnh(1+8veg,6bands,2hemi), band 1 is VIS.
#ifndef USE_ENT
      cosday=cos(twopi/edpery*jday)
      sinday=sin(twopi/edpery*jday)
#endif
      do j=J_0,J_1
        do i=I_0,I_1
          if(lat2d(i,j).lt.0.) then !nyk added northsouth
            northsouth=1            !southern hemisphere
          else
            northsouth=2            !northern hemisphere
          end if
          wfcs(i,j)=24.
          ! if (fearth(i,j).gt.0.) then
          !if (focean(i,j) < 1.d0 ) then
          if (variable_lk==0) then
            if ( fearth(i,j) <= 0.d0 ) cycle
          else
            if ( focean(i,j) >= 1.d0 ) cycle
          endif
#ifndef USE_ENT
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
#endif

!!!            call ghinij(i,j)
            do ibv=1,2
              call get_soil_properties( q_ij(i,j,:,:), dz_ij(i,j,:),
     &             thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )
            enddo

#ifndef USE_ENT
            call veg_set_cell(vegcell, i,j, 1.d0, 1.d0, .true.)
            !call veg_set_cell(i,j, .true.)
!!!            !ws_can = ws(0,2)
!!!            ws_can = vegcell%ws_can
            !ws_can = ws(0,2)
            ws_can = vegcell%ws_can
#else
            call ent_get_exports( entcells(i,j),canopy_max_H2O=ws_can )
#endif
            call get_fb_fv( fb, fv, i, j )
!!!            wfc1=fb*ws(1,1)+fv*(ws_can+ws(1,2))
cddd            wfc1=fb*thets(1,1)*dz_ij(i,j,1) +
cddd     &           fv*( ws_can + thets(1,2)*dz_ij(i,j,1) )
            ws11 = thets(1,1)*dz_ij(i,j,1)
            ws12 = thets(1,2)*dz_ij(i,j,1)
            wfc1=fb*ws11 +
     &           fv*( ws_can + ws12 )
            wfcs(i,j)=rhow*wfc1 ! canopy part changes
cddd            write(934,*) "wfcs", i,j,wfcs(i,j)

         !!! this diag belongs to Ent - commenting out
#ifndef USE_ENT
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
#endif
          !end if
        end do
      end do

      if (end_of_day) then
        do j=J_0,J_1
        do i=I_0,imaxj(j)
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
      use geom, only : imaxj,axyp
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use DOMAIN_DECOMP, only : GLOBALSUM
      use ghy_com, only : snowe, tearth,wearth,aiearth,w_ij
     *     ,snowbv,fr_snow_ij,fr_snow_rad_ij, gdeep, dzsn_ij, nsn_ij,
     *     fearth
      !use veg_com, only : afb
      use diag_com, only : aij=>aij_loc
     *     ,jreg,ij_evap,ij_f0e,ij_evape
     *     ,ij_gwtr,ij_tg1,j_wtr1,j_ace1,j_wtr2,j_ace2
     *     ,j_snow,j_evap,j_type,ij_g01,ij_g07,ij_g04,ij_g10,ij_g28
     *     ,ij_g29,j_rsnow,ij_rsnw,ij_rsit,ij_snow,ij_gice, ij_gwtr1
     &     ,ij_zsnow
      use fluxes, only : e0,e1,evapor,eprec
      implicit none

      real*8 snow,f0dt,f1dt,evap,wtr1,wtr2,ace1,ace2
     *     ,pearth,enrgp,scove,fb,fv
      integer i,j,jr,k

C**** define local grid
      integer J_0, J_1, J_0H, J_1H, I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0      , J_STOP=J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
      do i=I_0,imaxj(j)
      pearth=fearth(i,j)
      jr=jreg(i,j)
      if (pearth.gt.0) then

        snow=snowe(i,j)
        !tg1 = tearth(i,j)
        wtr1= wearth(i,j)
        ace1=aiearth(i,j)
        !tg2=gdeep(i,j,1)
        wtr2=gdeep(i,j,2)
        ace2=gdeep(i,j,3)
        f0dt=e0(i,j,4)
        f1dt=e1(i,j,4)
        evap=evapor(i,j,4)
        enrgp=eprec(i,j)      ! including latent heat
        call get_fb_fv( fb, fv, i, j )

c**** accumulate diagnostics
c**** the following is the actual snow cover of the snow model
c        scove = pearth *
c     *       ( afb(i,j)*fr_snow_ij(1,i,j)
c     *       + (1.-afb(i,j))*fr_snow_ij(2,i,j) )
c**** the following computes the snow cover as it is used in RAD_DRV.f
        scove = pearth *
     *       ( fb*fr_snow_rad_ij(1,i,j)
     *       + fv*fr_snow_rad_ij(2,i,j) )

        !if (snowe(i,j).gt.0.) scove=pearth
        call inc_aj(i,j,itearth,j_rsnow,scove)
        call inc_areg(i,j,jr,j_rsnow,scove*axyp(i,j))
        aij(i,j,ij_rsnw)=aij(i,j,ij_rsnw)+scove
        aij(i,j,ij_snow)=aij(i,j,ij_snow)+snow*pearth
        aij(i,j,ij_rsit)=aij(i,j,ij_rsit)+scove

        call inc_aj(i,j,itearth,j_wtr1,wtr1*pearth)
        call inc_aj(i,j,itearth,j_ace1,ace1*pearth)
        call inc_aj(i,j,itearth,j_wtr2,wtr2*pearth)
        call inc_aj(i,j,itearth,j_ace2,ace2*pearth)
        call inc_aj(i,j,itearth,j_snow,snow*pearth)
        call inc_areg(i,j,jr,j_snow,snow*pearth*axyp(i,j))
        call inc_areg(i,j,jr,j_wtr1,wtr1*pearth*axyp(i,j))
        call inc_areg(i,j,jr,j_ace1,ace1*pearth*axyp(i,j))
        call inc_areg(i,j,jr,j_wtr2,wtr2*pearth*axyp(i,j))
        call inc_areg(i,j,jr,j_ace2,ace2*pearth*axyp(i,j))

        aij(i,j,ij_f0e)  =aij(i,j,ij_f0e)  +f0dt+enrgp
        aij(i,j,ij_gwtr) =aij(i,j,ij_gwtr)+(wtr1+ace1+wtr2+ace2)
        aij(i,j,ij_gwtr1) =aij(i,j,ij_gwtr1)+(wtr1+ace1)
        aij(i,j,ij_gice) =aij(i,j,ij_gice)+(ace1+ace2)
        aij(i,j,ij_evape)=aij(i,j,ij_evape)+evap
        do k=1,3
          aij(i,j,ij_g01+k-1)=aij(i,j,ij_g01+k-1)+w_ij(k,1,i,j)
          aij(i,j,ij_g07+k-1)=aij(i,j,ij_g07+k-1)+w_ij(k-1,2,i,j)
        end do
        aij(i,j,ij_g04)=aij(i,j,ij_g04)+w_ij(6,1,i,j)
        aij(i,j,ij_g10)=aij(i,j,ij_g10)+w_ij(6,2,i,j)
        aij(i,j,ij_g28)=aij(i,j,ij_g28)+snowbv(1,i,j)
        aij(i,j,ij_g29)=aij(i,j,ij_g29)+snowbv(2,i,j)
        aij(i,j,ij_zsnow)=aij(i,j,ij_zsnow) + pearth *
     &       ( fb*fr_snow_ij(1,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(1,i,j),1,i,j) )
     &       + fv*fr_snow_ij(2,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(2,i,j),2,i,j) ) )
      end if
c****
      end do
      end do

      end subroutine ground_e


      end module soil_drv


      subroutine conserv_wtg(waterg)
!@sum  conserv_wtg calculates ground water incl snow
!@auth Gavin Schmidt
!@ver  1.0
      use constant, only : rhow
      use model_com, only : im, fim, focean, jm, flice
      use geom, only : imaxj,AXYP
      use ghy_com, only : ngm,w_ij,wsn_ij,fr_snow_ij,nsn_ij,fearth
      !use veg_com, only : afb
      use LAKES_COM, only : flake
      use LANDICE_COM,only : MDWNIMP
      USE DOMAIN_DECOMP, ONLY : GRID, GET, HERE
      implicit none
!@var waterg ground water (kg/m^2)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO),intent(out)::
     &     waterg

      integer i,j,n
      real*8 wij,fb,fv

C**** define local grid
      integer :: J_0, J_1 ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
      do i=I_0,imaxj(j)
          !if (fearth(i,j).gt.0) then
          if( focean(i,j) + flice(i,j) < 1.d0 ) then
            !fb=afb(i,j)
            !fv=(1.d0-fb)
            call get_fb_fv( fb, fv, i, j )
            wij=fb*sum( w_ij(1:ngm,1,i,j) )
     &       +  fv*sum( w_ij(0:ngm,2,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( wsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( wsn_ij(1:nsn_ij(2,i,j),2,i,j))
            waterg(i,j)=fearth(i,j)*wij*rhow
     &           + flake(i,j)*sum( w_ij(0:ngm,3,i,j) )*rhow
!!! hack to check remove_extra_snow
!!!            waterg(j)=waterg(j)+MDWNIMP(i,j)/axyp(i,j)
          else
            waterg(i,j)=0
          end if
       end do
      end do
      if (HAVE_SOUTH_POLE) waterg(2:im,1) =waterg(1,1)
      if (HAVE_NORTH_POLE) waterg(2:im,jm)=waterg(1,jm)
c****
      end subroutine conserv_wtg


      subroutine conserv_wtg_1(waterg,fearth,flake)
!@sum  conserv_wtg calculates ground water incl snow
!@auth Gavin Schmidt
!@ver  1.0
      use constant, only : rhow
      use model_com, only : im, fim, focean, im, jm
      use geom, only : imaxj
      use ghy_com, only : ngm,w_ij,wsn_ij,fr_snow_ij,nsn_ij
      !use veg_com, only : afb
      !use LAKES_COM, only : flake
      USE DOMAIN_DECOMP, ONLY : GRID, GET, HERE
      implicit none
      real*8,dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                 GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &     intent(in) :: fearth, flake
!@var waterg ground water (kg/m^2)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO),intent(out)::
     &     waterg

      integer i,j,n
      real*8 wij,fb,fv

C**** define local grid
      integer :: J_0, J_1 ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          !if (fearth(i,j).gt.0) then
          if( focean(i,j) < 1.d0 ) then
            !fb=afb(i,j)
            !fv=(1.d0-fb)
            call get_fb_fv( fb, fv, i, j )
            wij=fb*sum( w_ij(1:ngm,1,i,j) )
     &       +  fv*sum( w_ij(0:ngm,2,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( wsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( wsn_ij(1:nsn_ij(2,i,j),2,i,j))
            waterg(i,j)=fearth(i,j)*wij*rhow
     &           + flake(i,j)*sum( w_ij(0:ngm,3,i,j) )*rhow
          else
            waterg(i,j)=0
          end if
       end do
      end do

      if (HAVE_SOUTH_POLE) waterg(2:im,1) =waterg(1,1)
      if (HAVE_NORTH_POLE) waterg(2:im,jm)=waterg(1,jm)
c****
      end subroutine conserv_wtg_1



      subroutine conserv_htg(heatg)
!@sum  conserv_htg calculates ground energy incl. snow energy
!@auth Gavin Schmidt
!@ver  1.0
      use model_com, only : im, fim, focean, jm, flice
      use geom, only : imaxj, dxyp
      use ghy_com, only : ngm,ht_ij,fr_snow_ij,nsn_ij,hsn_ij
     *     ,fearth
      !use veg_com, only : afb
      use LAKES_COM, only : flake
      USE DOMAIN_DECOMP, ONLY : GRID, GET, HERE
      implicit none
!@var heatg ground heat (J/m^2)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: heatg

      integer i,j
      real*8 hij,fb,fv

C**** define local grid
      integer J_0, J_1 ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          heatg(i,j)=0
          !if (fearth(i,j).le.0) cycle
          if ( focean(i,j) + flice(i,j) >= 1.d0 ) cycle
          !fb=afb(i,j)
          !fv=(1.d0-fb)
          call get_fb_fv( fb, fv, i, j )
          hij=fb*sum( ht_ij(1:ngm,1,i,j) )
     &       +  fv*sum( ht_ij(0:ngm,2,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( hsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( hsn_ij(1:nsn_ij(2,i,j),2,i,j))
          heatg(i,j)=fearth(i,j)*hij
     &           + flake(i,j)*sum( ht_ij(0:ngm,3,i,j) )
        end do
      end do
      if (HAVE_SOUTH_POLE) heatg(2:im,1) =heatg(1,1)
      if (HAVE_NORTH_POLE) heatg(2:im,jm)=heatg(1,jm)
c****
ccc debugging ...
ccc      print *,'conserv_htg energy ',
ccc     &     sum(heatg(1:jm)*dxyp(1:jm))/(sum(dxyp(1:jm))*im)
      end subroutine conserv_htg



      subroutine check_ghy_conservation( flag )
ccc debugging program: cam be put at the beginning and at the end
ccc of the 'surface' to check water conservation
      use constant, only : rhow
      use geom, only : imaxj
      use model_com, only : im,jm
      use DOMAIN_DECOMP, only : GRID, GET
      use fluxes, only : prec,evapor,runoe
      use ghy_com, only : ngm,w_ij,ht_ij,snowbv,dz_ij
     *     ,fearth
      !use veg_com, only : afb
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
      integer J_0, J_1 ,I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          if ( fearth(i,j) <= 0.d0 ) cycle

ccc just checking ...
          do n = 1,ngm
            if ( dz_ij(i,j,n) .le. 0.d0 )
     &           call stop_model('incompatible dz',255)
          enddo

          !fb = afb(i,j)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, i, j )
          total_water(i,j) = fb*sum( w_ij(1:ngm,1,i,j) )
     &         + fv*sum( w_ij(0:ngm,2,i,j) )
     &         + fb*snowbv(1,i,j) + fv*snowbv(2,i,j)
        end do
      end do

      ! call stop_model('just testing...',255)

      if ( flag == 0 ) then
        old_total_water(:,:) = total_water(:,:)
        return
      endif

      do j=J_0,J_1
        do i=I_0,imaxj(j)

          !print *,'fearth = ', i, j, fearth(i,j)

          if ( fearth(i,j) <= 0.d0 ) cycle
          !fb = afb(i,j)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, i, j )
          error_water = ( total_water(i,j) - old_total_water(i,j) )*rhow
     &         - prec(i,j) + evapor(i,j,4) + runoe(i,j)

          !print *, 'err H2O: ', i, j, error_water

  !        if ( abs( error_water ) > 1.d-9 ) print *, 'error'
          if ( abs( error_water ) > 1.d-9 ) call stop_model(  ! was -15
     &         'check_ghy_conservation: water conservation problem',255)

        end do
      end do

      end subroutine check_ghy_conservation


      subroutine compute_water_deficit
      use constant, only : twopi,edpery,rhow
      use ghy_com, only : ngm,imt,LS_NFRAC,dz_ij,q_ij
     &     ,w_ij,fearth
      !use veg_com, only : ala !,afb
      use model_com, only : focean
      use sle001, only : thm
      use fluxes, only : DMWLDF
      USE DOMAIN_DECOMP, ONLY : GRID, GET

      implicit none
cddd      integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1
      integer k,ibv,m
      real*8 :: w_tot(2),w_stor(2)
      real*8 :: w(0:ngm,2),dz(ngm),q(imt,ngm)
cddd      real*8 :: cosday,sinday,alai
      real*8 :: fb,fv

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

cddd      cosday=cos(twopi/edpery*jday)
cddd      sinday=sin(twopi/edpery*jday)

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

          w(0:ngm,1:2) = w_ij(0:ngm,1:2,i,j)
          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

          !fb = afb(i,j)
          !fv=1.-fb
          call get_fb_fv( fb, fv, i, j )

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
cddd          alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
cddd          alai=max(alai,1.d0)
cddd          w_stor(2) = w_stor(2) + .0001d0*alai
          w_tot(2) = w_tot(2) + w(0,2)

          ! total water deficit on kg/m^2
          DMWLDF(i,j) = fb*(w_stor(1) - w_tot(1))
     &         + fv*(w_stor(2) - w_tot(2))
          DMWLDF(i,j) = max( DMWLDF(i,j), 0.d0 )
        enddo
      enddo


      !print *,"DMWLDF=",DMWLDF

      end subroutine compute_water_deficit


      subroutine init_underwater_soil

!!!! UNFINISHED
      use constant, only : twopi,edpery,rhow,shw_kg=>shw
      use ghy_com, only : ngm,imt,LS_NFRAC,dz_ij,q_ij
     &     ,w_ij,ht_ij,fearth,shc_soil_texture
#ifdef TRACERS_WATER
     &     ,tr_w_ij
      use TRACER_COM, only : ntm,needtrs,itime_tr0
      use model_com, only : itime
#endif
#ifdef USE_ENT
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
#else
      use veg_com, only : ala !,afb
#endif
      use model_com, only : focean
      use sle001, only : thm
      use fluxes, only : DMWLDF
      USE DOMAIN_DECOMP, ONLY : GRID, GET

      implicit none
cddd      integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1
      integer k,ibv,m
      real*8 :: w_stor(0:ngm), ht_cap(0:ngm)
      real*8 :: w(0:ngm,2),dz(ngm),q(imt,ngm)
      !!real*8 :: cosday,sinday,alai
      real*8 :: fb,fv
      real*8 shc_layer, aa, tp, tpb, tpv, ficeb, ficev, fice
#ifdef TRACERS_WATER
      integer n
      real*8 wsoil_tot
#endif

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)


cddd      cosday=cos(twopi/edpery*jday)
cddd      sinday=sin(twopi/edpery*jday)

      do j=J_0,J_1
        do i=I_0,I_1

          if( focean(i,j) >= 1.d0 ) then
            w_ij (0:ngm,3,i,j) = 0.d0
            ht_ij(0:ngm,3,i,j) = 0.d0
            cycle
          endif

          !w(0:ngm,1:2) = w_ij(0:ngm,1:2,i,j)
          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

#ifdef USE_ENT
          call ent_get_exports( entcells(i,j),
     &         canopy_heat_capacity=ht_cap(0) )
#endif
          !fb = afb(i,j)
          !fv=1.-fb
          call get_fb_fv( fb, fv, i, j )

          ! compute max water storage and heat capacity
          do k=1,ngm
            w_stor(k) = 0.d0
            do m=1,imt-1
              w_stor(k) = w_stor(k) + q(m,k)*thm(0,m)*dz(k)
            enddo
            shc_layer = 0.d0
            do m=1,imt
              shc_layer = shc_layer + q(m,k)*shc_soil_texture(m)
            enddo
            ht_cap(k) = (dz(k)-w_stor(k)) * shc_layer
          enddo

          ! include canopy water here
          !!alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
          !!alai=max(alai,1.d0)
          !!w_stor(0) = .0001d0*alai*fv
          w_stor(0) = 0.d0
!! set above
#ifndef USE_ENT
          aa=ala(1,i,j)
          ht_cap(0)=(.010d0+.002d0*aa+.001d0*aa**2)*shw_kg*rhow
#endif
          ! we will use as a reference average temperature of the
          ! lowest layer
          call heat_to_temperature( tpb, ficeb,
     &         ht_ij(ngm,1,i,j), w_ij(ngm,1,i,j), ht_cap(ngm) )
          call heat_to_temperature( tpv, ficev,
     &         ht_ij(ngm,2,i,j), w_ij(ngm,2,i,j), ht_cap(ngm) )
          tp   = fb*tpb + fv*tpv
          fice = fb*ficeb + fv*ficev

          ! set underground fraction to tp, fice and saturated water
          !! nothing in canopy layer (ie. temp = 0C)
          w_ij(0,3,i,j) = 0.d0
          ht_ij(0,3,i,j) = 0.d0
          do k=1,ngm
            w_ij(k,3,i,j) = w_stor(k)
            call temperature_to_heat( ht_ij(k,3,i,j),
     &           tp, fice, w_ij(k,3,i,j), ht_cap(k) )
          enddo

#ifdef TRACERS_WATER
      ! set underwater tracers to average land tracers (ignore canopy)
          do k=1,ngm
            wsoil_tot = fb*w_ij(k,1,i,j) + fv*w_ij(k,2,i,j)
            do n=1,ntm
              if (itime_tr0(n).gt.itime) cycle
              if ( .not. needtrs(n) ) cycle
              tr_w_ij(n,k,3,i,j) = 0.d0
              if ( wsoil_tot > 1.d-30 ) then
                tr_w_ij(n,k,3,i,j) = (
     &               fb*tr_w_ij(n,k,1,i,j) + fv*tr_w_ij(n,k,2,i,j)
     &               ) / (wsoil_tot) * w_ij(k,3,i,j) !!! removed /rhow
              end if
            enddo
          enddo
#endif


        enddo
      enddo

      end subroutine init_underwater_soil


      subroutine heat_to_temperature(tp, fice, ht, w, ht_cap)
      use constant, only : rhow, lhm, shw_kg=>shw, shi_kg=>shi
      real*8, intent(out) :: tp, fice
      real*8, intent(in) :: ht, w, ht_cap
      ! volumetric quantities
      real*8, parameter :: lhmv=lhm*rhow, shwv=shw_kg*rhow,
     &     shiv=shi_kg*rhow

      fice = 0.d0
      if( lhmv*w+ht < 0.d0 ) then ! all frozen
        tp = ( ht + w*lhmv )/( ht_cap + w*shiv )
        fice = 1.d0
      else if( ht > 0.d0 ) then ! all melted
        tp = ht /( ht_cap + w*shwv )
      else  ! part frozen
        tp = 0.d0
        if( w > 1d-12 ) fice = -ht /( lhmv*w )
      endif

      end subroutine heat_to_temperature


      subroutine temperature_to_heat(ht, tp, fice, w, ht_cap)
      use constant, only : rhow, lhm, shw_kg=>shw, shi_kg=>shi
      real*8, intent(in) :: tp, fice, w, ht_cap
      real*8, intent(out) :: ht
      ! volumetric quantities
      real*8, parameter :: lhmv=lhm*rhow, shwv=shw_kg*rhow,
     &     shiv=shi_kg*rhow
      real*8 lht_ice

      lht_ice = fice*w*lhmv

      if( tp > 0.d0 ) then
        ht = tp*( ht_cap + w*shwv ) - lht_ice
      else
        ht = tp*( ht_cap + w*shiv ) - lht_ice
      endif

      end subroutine temperature_to_heat


      subroutine update_land_fractions !(jday)

!!!! UNFINISHED
      use constant, only : twopi,edpery,rhow,shw_kg=>shw
      use ghy_com, only : ngm,imt,dz_ij,q_ij
     &     ,w_ij,ht_ij,fr_snow_ij,fearth
     &     ,fr_snow_rad_ij,snowbv,top_dev_ij
#ifdef TRACERS_WATER
     &     ,tr_w_ij,tr_wsn_ij
      use TRACER_COM, only : ntm
#endif
      !use veg_com, only : ala !,afb
      use model_com, only : focean,im
      use LAKES_COM, only : flake, svflake
      use sle001, only : thm
      use fluxes, only : DMWLDF, DGML
#ifdef TRACERS_WATER
     &     ,DTRL
#endif
      use GEOM, only : BYAXYP
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      use soil_drv, only : snow_cover ! conserv_wtg_1
      use snow_drvm, only : snow_cover_same_as_rad

      implicit none
      !! integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1,ibv
      integer k,m
      real*8 :: w_stor(0:ngm)
      real*8 :: dz(ngm),q(imt,ngm)
      !!real*8 :: cosday,sinday,alai
      real*8 :: fb,fv
      real*8 :: dw, dw_soil, dw_lake, dht, dht_soil, dht_lake
      real*8 :: sum_water, ht_per_m3

      real*8 dfrac
#ifdef TRACERS_WATER
      real*8 dtr_soil(ntm),dtr_lake(ntm),dtr(ntm),tr_per_m3(ntm)
#endif
      real*8 tmp_before(0:ngm), tmp_after(0:ngm)
c      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
c     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
c     &     w_before_j, w_after_j
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     fearth_old

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)


      fearth_old = fearth + (flake - svflake)

      ! call conserv_wtg_1(w_before_j,fearth_old,svflake)

      do j=J_0,J_1
        do i=I_0,I_1

          if( focean(i,j) >= 1.d0 ) cycle

          if ( svflake(i,j) == flake(i,j) ) cycle

          call get_fb_fv( fb, fv, i, j )

          if ( flake(i,j) < svflake(i,j) ) then ! lake shrunk
            ! no external fluxes, just part of underwater fraction
            ! is transformed into a land fraction

            ! no changes to underwater fraction
            ! just redistribute added quantities over land fractions


            dfrac = svflake(i,j) - flake(i,j)
            ! make sure old fearth >= 0 (i.e. remove round-off errors)
            dfrac = min(dfrac, fearth(i,j))

            tmp_before = ( ht_ij(0:ngm,3,i,j) )*svflake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)-dfrac) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j)-dfrac)

            do k=1,ngm
              dw  = dfrac*w_ij(k,3,i,j)
              dht = dfrac*ht_ij(k,3,i,j)
              w_ij(k,1:2,i,j) =
     &             (w_ij(k,1:2,i,j)*(fearth(i,j)-dfrac) + dw)
     &             / fearth(i,j)
              ht_ij(k,1:2,i,j) =
     &             (ht_ij(k,1:2,i,j)*(fearth(i,j)-dfrac) + dht)
     &             / fearth(i,j)
#ifdef TRACERS_WATER
              dtr(:) = dfrac*tr_w_ij(:,k,3,i,j)
              do ibv=1,2
                tr_w_ij(:,k,ibv,i,j) =
     &               (tr_w_ij(:,k,ibv,i,j)*(fearth(i,j)-dfrac) + dtr(:))
     &               / fearth(i,j)
              enddo
#endif
            enddo
            !vegetation:
cddd            if ( fv > 0.d0 ) then
cddd              w_ij(0,2,i,j) =
cddd     &             (w_ij(0,2,i,j)*(fearth(i,j)-dfrac) + dw/fv)
cddd     &             / fearth(i,j)
cddd              ht_ij(0,2,i,j) =
cddd     &             (ht_ij(0,2,i,j)*(fearth(i,j)-dfrac) + dht/fv)
cddd     &             / fearth(i,j)
cddd            endif
              w_ij(0,2,i,j) =
     &             (w_ij(0,2,i,j)*(fearth(i,j)-dfrac))
     &             / fearth(i,j)
              ht_ij(0,2,i,j) =
     &             (ht_ij(0,2,i,j)*(fearth(i,j)-dfrac))
     &             / fearth(i,j)
#ifdef TRACERS_WATER
              tr_w_ij(:,0,2,i,j) =
     &             (tr_w_ij(:,0,2,i,j)*(fearth(i,j)-dfrac))
     &             / fearth(i,j)
#endif

            ! change snow fraction
            fr_snow_ij(1:2,i,j) =
     &           fr_snow_ij(1:2,i,j)*(1.d0-dfrac/fearth(i,j))
#ifdef TRACERS_WATER
            ! tr_wsn is spread over entire cell (i.e. *fr_snow)
            tr_wsn_ij(:,:,1:2,i,j) = tr_wsn_ij(:,:,1:2,i,j) *
     &           (1.d0 - dfrac/fearth(i,j))
#endif

            tmp_after = ( ht_ij(0:ngm,3,i,j) )*flake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j))
            !print *,"after", (tmp_after), (tmp_after-tmp_before)

          else if ( flake(i,j) > svflake(i,j) ) then ! lake expanded

            ! no need to change land values
            ! for underwater fraction:

            dz(1:ngm) = dz_ij(i,j,1:ngm)
            q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

            ! compute max water storage and heat capacity
            do k=1,ngm
              w_stor(k) = 0.d0
              do m=1,imt-1
                w_stor(k) = w_stor(k) + q(m,k)*thm(0,m)*dz(k)
              enddo
            enddo
            ! include canopy water here
!!!cddd            alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
!!!cddd            alai=max(alai,1.d0)
!!!cddd            w_stor(0) = .0001d0*alai*fv
            ! no uderlake water in canopy
            w_stor(0) = 0.d0
            ! allow any amount of water in upper soil layer
            w_stor(1) = 1.d30

            dfrac = flake(i,j) - svflake(i,j)

            tmp_before = ( ht_ij(0:ngm,3,i,j) )*svflake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)+dfrac) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j)+dfrac)

            sum_water = DMWLDF(i,j)*dfrac/rhow
            if ( sum_water > 1.d-30 ) then
              ht_per_m3 = DGML(i,j)*BYAXYP(I,J)/sum_water
            else
              ht_per_m3 = 0.d0
            endif
#ifdef TRACERS_WATER
            if ( sum_water > 1.d-30 ) then
              tr_per_m3(:) = DTRL(:,i,j)*BYAXYP(I,J)/sum_water
            else
              tr_per_m3(:) = 0.d0
            endif
#endif
            do k=ngm,1,-1  ! do not loop over canopy
              ! dw, dht - total amounts of water and heat added to
              ! underwater fraction
              dw_soil = dfrac*( fb*w_ij(k,1,i,j) + fv*w_ij(k,2,i,j) )
              dw  = min( dfrac*w_stor(k), sum_water + dw_soil )
              dw_lake = dw - dw_soil
              dht_soil = dfrac*( fb*ht_ij(k,1,i,j) + fv*ht_ij(k,2,i,j) )
              dht_lake = dw_lake*ht_per_m3
              dht = dht_soil + dht_lake
              sum_water = sum_water - dw_lake

              w_ij(k,3,i,j) =
     &             (w_ij(k,3,i,j)*svflake(i,j) + dw)/flake(i,j)
              ht_ij(k,3,i,j) =
     &             (ht_ij(k,3,i,j)*svflake(i,j) + dht)/flake(i,j)
#ifdef TRACERS_WATER
              dtr_soil(:) =
     &             dfrac*(fb*tr_w_ij(:,k,1,i,j) + fv*tr_w_ij(:,k,2,i,j))
              dtr_lake(:) = dw_lake*tr_per_m3(:)
              dtr(:) = dtr_soil(:) + dtr_lake(:)
              tr_w_ij(:,k,3,i,j) =
     &             (tr_w_ij(:,k,3,i,j)*svflake(i,j) + dtr(:))/flake(i,j)
#endif
            enddo

            ! dump canopy water into first layer (underwater canopy is 0C)
            dw = dfrac*( fv*w_ij(0,2,i,j) )
            dht = dfrac*( fv*ht_ij(0,2,i,j) )
              w_ij(1,3,i,j) =
     &             w_ij(1,3,i,j) + dw/flake(i,j)
              ht_ij(1,3,i,j) =
     &             ht_ij(1,3,i,j) + dht/flake(i,j)
#ifdef TRACERS_WATER
            dtr(:) = dfrac*( fv*tr_w_ij(:,0,2,i,j) )
            tr_w_ij(:,1,3,i,j) =
     &             tr_w_ij(:,1,3,i,j) + dtr(:)/flake(i,j)
#endif
            tmp_after = ( ht_ij(0:ngm,3,i,j) )*flake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j))

            ! change snow fraction
            if ( fearth(i,j) <= 0.d0 ) then
              print *, "farctions:",i,j,
     &             focean(i,j), fearth(i,j), flake(i,j), svflake(i,j)
              call stop_model("update_land_fractions: fearth<=0",255)
            endif

            fr_snow_ij(1:2,i,j) =
     &           fr_snow_ij(1:2,i,j)*(1.d0+dfrac/fearth(i,j))
#ifdef TRACERS_WATER
            ! tr_wsn is spread over entire cell (i.e. *fr_snow)
            tr_wsn_ij(:,:,1:2,i,j) = tr_wsn_ij(:,:,1:2,i,j) *
     &           (1.d0 + dfrac/fearth(i,j))
#endif


            ! hack to deal with snow in empty fractions (fb, fv)
            if( fb <= 0.d0 )
     &           fr_snow_ij(1,i,j) = min( .95d0, fr_snow_ij(1,i,j) )
            if( fv <= 0.d0 )
     &           fr_snow_ij(2,i,j) = min( .95d0, fr_snow_ij(2,i,j) )

            !print *,"FR_SNOW  after" ,i,j,fr_snow_ij(1:2,i,j)
            if ( fr_snow_ij(1,i,j) > 1.d0 .or.
     &           fr_snow_ij(2,i,j) > 1.d0 ) call stop_model(
     &           "update_land_fractions: fr_snow_ij > 1",255)

          endif

c**** Also reset snow fraction for albedo computation
          if ( snow_cover_same_as_rad == 0 ) then
            ! recompute snow fraction using different formula
            do ibv=1,2
               call snow_cover(fr_snow_rad_ij(ibv,i,j),
     &                snowbv(ibv,i,j), top_dev_ij(i,j) )
               fr_snow_rad_ij(ibv,i,j) = min (
     &            fr_snow_rad_ij(ibv,i,j), fr_snow_ij(ibv, i, j) )
            enddo
          else
            ! snow fraction same as in snow model
            fr_snow_rad_ij(:,i,j) = fr_snow_ij(:, i, j)
          endif

          call set_new_ghy_cells_outputs

          ! reset "FLUXES" arrays
          svflake(i,j) = flake(i,j)
          DMWLDF(i,j) = 0.d0
          DGML(i,j) = 0.d0
        enddo
      enddo

      !call conserv_wtg_1(w_after_j,fearth,flake)
      !print *,"UPDATE_LF CONS_WTRG ", w_after_j-w_before_j


      end subroutine update_land_fractions


      subroutine set_new_ghy_cells_outputs
!@sum set output data for newly created earth cells (when lake shrinks)
      use constant, only : rhow,tf,lhe,lhs
      use ghy_com, only : ngm,imt,dz_ij,q_ij
     &     ,w_ij,ht_ij,fr_snow_ij,fearth,qg_ij,fr_snow_rad_ij
     &     ,shc_soil_texture,canopy_temp_ij,snowe,tearth,wearth,aiearth
     &     ,tsns_ij
#ifdef TRACERS_WATER
     &     ,tr_w_ij
      use TRACER_COM, only : ntm,needtrs,itime_tr0
      use model_com, only : itime
#endif
      use FLUXES, only : gtemp,gtempr
#ifdef TRACERS_WATER
     &     ,gtracer
#endif
      !use veg_com, only : afb
      use sle001, only : thm
      use LAKES_COM, only : flake, svflake
      use dynamics, only : pedn
      USE DOMAIN_DECOMP, ONLY : GRID, GET

      implicit none
      !---
      real*8, external :: qsat
      real*8 :: EPS=1.d-12
      integer i,j,I_0,I_1,J_0,J_1
      real*8 :: dz(ngm), q(imt,ngm), w_stor(ngm), ht_cap(ngm)
      real*8 tg1, fice, fb, fv, tpb, tpv, ficeb, ficev, shc_layer
      real*8 dfrac, elhx, ps
      integer k, m
#ifdef TRACERS_WATER
      integer n
      real*8 wsoil_tot
#endif

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

      do j=J_0,J_1
        do i=I_0,I_1

          dfrac = svflake(i,j) - flake(i,j)
          if ( fearth(i,j) < EPS .or. fearth(i,j)-dfrac > EPS ) cycle

          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

          !fb = afb(i,j)
          !fv=1.-fb
          call get_fb_fv( fb, fv, i, j )

      ! compute max water storage and heat capacity
          do k=1,ngm
            w_stor(k) = 0.d0
            do m=1,imt-1
              w_stor(k) = w_stor(k) + q(m,k)*thm(0,m)*dz(k)
            enddo
            shc_layer = 0.d0
            do m=1,imt
              shc_layer = shc_layer + q(m,k)*shc_soil_texture(m)
            enddo
            ht_cap(k) = (dz(k)-w_stor(k)) * shc_layer
          enddo

          call heat_to_temperature( tpb, ficeb,
     &         ht_ij(1,1,i,j), w_ij(1,1,i,j), ht_cap(ngm) )
          call heat_to_temperature( tpv, ficev,
     &         ht_ij(1,2,i,j), w_ij(1,2,i,j), ht_cap(ngm) )
          tg1   = fb*tpb + fv*tpv
          fice = fb*ficeb + fv*ficev

          elhx=lhe
          if(tg1.lt.0.)  elhx=lhs
          ps=pedn(1,i,j)
      ! ground humidity to be used on next time step
          qg_ij(i,j) = qsat(tg1+tf,elhx,ps) ! all saturated
      ! snow fraction same as in snow model
          fr_snow_rad_ij(:,i,j) = 0.d0 ! no snow in new cell
      ! canopy_temp_ij is not used so far... do we need it?
          canopy_temp_ij(i,j) = 0.d0 ! canopy at 0C

c**** snowe used in RADIATION
          snowe(i,j) = 1000.*0.d0
c**** tearth used only internaly in GHY_DRV
          tearth(i,j) = sqrt(sqrt(fb*(tpb+tf)**4 + fv*(tpv+tf)**4)) - tf
          tsns_ij(i,j) = tg1
c**** wearth+aiearth are used in radiation only
          wearth(i,j)=1000.*( fb*w_ij(1,1,i,j)*(1.-ficeb) +
     &         fv*(w_ij(1,2,i,j)*(1.-ficev)) )
          aiearth(i,j)=1000.*( fb*w_ij(1,1,i,j)*ficeb +
     &         fv*w_ij(1,2,i,j)*ficev )
          gtemp(1,4,i,j)=tsns_ij(i,j)
          gtempr(4,i,j) =tearth(i,j)+tf

#ifdef TRACERS_WATER
      ! I use vegetated ground insted of canopy since canopy has 0 H2O
          wsoil_tot = fb*w_ij(1,1,i,j) + fv*w_ij(1,2,i,j)
          do n=1,ntm
            if (itime_tr0(n).gt.itime) cycle
            if ( .not. needtrs(n) ) cycle
            ! should also restrict to TYPE=nWATER ?
            if ( wsoil_tot > 1.d-30 ) then
              gtracer(n,4,i,j) = (
     &             fb*tr_w_ij(n,1,1,i,j) + fv*tr_w_ij(n,1,2,i,j)
     &             ) / (rhow*wsoil_tot)
            else
              gtracer(n,4,i,j) = 0.
            end if
         !!! test
            gtracer(n,4,i,j) = 1.d0
          enddo
#endif
        enddo
      enddo

      end subroutine set_new_ghy_cells_outputs


      subroutine remove_extra_snow_to_ocean ! bad idea ?
      use constant, only : rhow
      use ghy_com, only : nsn_ij, dzsn_ij, wsn_ij, hsn_ij,
     &     fr_snow_ij,fearth,w_ij,ngm,wsn_max
#ifdef TRACERS_WATER
     &     ,tr_wsn_ij
      use TRACER_COM, only : ntm
#endif
      !use veg_com, only : afb
      use LAKES_COM, only : flake
      use LANDICE_COM, only : MDWNIMP, EDWNIMP
#ifdef TRACERS_WATER
     &     ,TRDWNIMP
#endif
      use GEOM, only : AXYP
      use MODEL_COM, only : ITEARTH
      USE DIAG_COM, only : j_implh,j_implm
     *     ,JREG
      USE DOMAIN_DECOMP, ONLY : GRID, GET, GLOBALSUM

      implicit none
      !! integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1,J_0H,J_1H
      real*8 fbv(2),wsn(3),hsn(3),dzsn(3),fr_snow,wsn_tot,d_wsn,eta
      real*8 dw,dh   ! ,tot0
      integer ibv,nsn,JR

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1
     &        ,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      do j=J_0,J_1
        do i=I_0,I_1

          if( fearth(i,j) <= 0.d0 ) cycle

          JR=JREG(I,J)

          !fbv(1) = afb(i,j)
          !fbv(2) =1.d0 - fbv(1)
          call get_fb_fv( fbv(1), fbv(2), i, j )

c          tot0=(fbv(1)*sum( w_ij(1:ngm,1,i,j) )
c     &         +  fbv(2)*sum( w_ij(0:ngm,2,i,j) )
c     &         +  fbv(1)*fr_snow_ij(1,i,j)*sum( wsn_ij(1:nsn_ij(1,i,j),1
c     *         ,i,j))+  fbv(2)*fr_snow_ij(2,i,j)*sum( wsn_ij(1:nsn_ij(2
c     *         ,i,j),2,i,j)))*fearth(i,j)*rhow+flake(i,j)*sum(
c     *         w_ij(0:ngm,3,i,j) )*rhow


          do ibv=1,2
            if ( fbv(ibv) <= 0.d0 ) cycle

            nsn = nsn_ij(ibv, i, j)
            wsn = 0 ; hsn = 0 ; dzsn = 0
            wsn(1:nsn) = wsn_ij(1:nsn, ibv, i, j)
            hsn(1:nsn) = hsn_ij(1:nsn, ibv, i, j)
            dzsn(1:nsn) = dzsn_ij(1:nsn, ibv, i, j)
            fr_snow = fr_snow_ij(ibv, i, j)
            wsn_tot = sum( wsn(1:nsn) )
            if ( wsn_tot > WSN_MAX ) then
              ! check if snow structure ok for thick snow
              if ( nsn < 3)
     &             call stop_model("remove_extra_snow: nsn<3",255)
              d_wsn = wsn_tot - WSN_MAX
              ! fraction of snow to be removed:
              eta = d_wsn/(wsn(2)+wsn(3))
              wsn_ij(2:3, ibv, i, j)  = (1.d0-eta)*wsn(2:3)
              hsn_ij(2:3, ibv, i, j)  = (1.d0-eta)*hsn(2:3)
              dzsn_ij(2:3, ibv, i, j) = (1.d0-eta)*dzsn(2:3)
              ! extra water and energy
              dw = eta*(wsn(2)+wsn(3))*fr_snow*fbv(ibv)*rhow
              MDWNIMP(i,j) = MDWNIMP(i,j) + dw*fearth(i,j)*axyp(i,j)
              dh = eta*(hsn(2)+hsn(3))*fr_snow*fbv(ibv)
              EDWNIMP(i,j) = EDWNIMP(i,j) + dh*fearth(i,j)*axyp(i,j)
#ifdef TRACERS_WATER
              TRDWNIMP(1:ntm,i,j) = TRDWNIMP(1:ntm,i,j) +
     &             eta*(
     &             tr_wsn_ij(1:ntm,2,ibv,i,j)+tr_wsn_ij(1:ntm,3,ibv,i,j)
     &             )*fbv(ibv)*fearth(i,j)*axyp(i,j)
              tr_wsn_ij(1:ntm,2:3,ibv,i,j) =
     &             (1.d0-eta)*tr_wsn_ij(1:ntm,2:3,ibv,i,j)
#endif
              CALL INC_AJ(I,J,ITEARTH,J_IMPLH,dh*fearth(i,j))
              CALL INC_AJ(I,J,ITEARTH,J_IMPLM,dw*fearth(i,j))
              CALL INC_AREG(I,J,JR,J_IMPLH,dh*fearth(i,j)*axyp(i,j))
              CALL INC_AREG(I,J,JR,J_IMPLM,dw*fearth(i,j)*axyp(i,j))

              !print *,"remove_extra_snow", i,j,ibv,wsn_tot,eta,dw,dh

            endif
          enddo
c          print
c     *         *,i,j,tot0,(fbv(1)*sum( w_ij(1:ngm,1,i,j) )+  fbv(2)*sum(
c     *         w_ij(0:ngm,2,i,j) )+  fbv(1)*fr_snow_ij(1,i,j)*sum(
c     *         wsn_ij(1:nsn_ij(1,i,j),1,i,j))+  fbv(2)*fr_snow_ij(2,i,j)
c     *         *sum( wsn_ij(1:nsn_ij(2,i,j),2,i,j)))*fearth(i,j)*rhow
c     *         +flake(i,j)*sum(w_ij(0:ngm,3,i,j) )*rhow

        enddo
      enddo

      end subroutine remove_extra_snow_to_ocean


      subroutine get_fb_fv( fb, fv, i, j )
!@sum this is a hack to hyde dependence on Ent/non-Ent vegetation
!@+   in most of the code. It returns fb,fv - fractions of bare and
!@+   vegetated soil
#ifdef USE_ENT
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
#else
      use veg_com, only : afb
#endif
      implicit none
      real*8, intent(out) :: fb, fv
      integer, intent(in) :: i, j

#ifdef USE_ENT
      call ent_get_exports( entcells(i,j),
     &     fraction_of_vegetated_soil=fv )
      fb = 1. - fv
#else
      fb = afb(i,j)
      fv=1.d0 - fb
#endif
      end subroutine get_fb_fv

