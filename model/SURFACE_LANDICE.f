! Lakes = itype == 1
! Look in GHY_DRV to see which tracers in SURFACE.f are needed for LANDICE surface type

C****   
C**** SURFACE.f    SURFACE fluxes    2006/12/21
C****
#include "rundeck_opts.h"



      SUBROUTINE SURFACE_LANDICE (ns,moddsf,moddd,TGRND,TGRN2,TGR4,e1)
!      integer, intent(in) :: ns,moddsf,moddd
! TGRND...e1 are passed in from SURFACE.f (it's a local variable there)

!@sum SURFACE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : rgas,lhm,lhe,lhs
     *     ,sha,tf,rhow,shv,shi,stbo,bygrav,by6
     *     ,deltx,teeny,rhows,grav,syr
#ifdef mjo_subdd
     *     ,undef
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     &     ,By3
#endif
      USE ATM_COM, only : u,v,t,p,q
      USE MODEL_COM, only : dtsrc,idacc,nday,itime,jhour,qcheck,jdate
#ifdef mjo_subdd
     *     ,lm
#endif
#ifdef SCM
      USE SCMDIAG, only : EVPFLX,SHFLX
      USE SCMCOM, only : iu_scm_prt, ALH, ASH, SCM_SURFACE_FLAG
     &     ,I_TARG,J_TARG
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, GLOBALSUM
      USE GEOM, only : axyp,imaxj,byaxyp,lat2d
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE ATM_COM, only : pmid,pk,pedn,pek,am,byam
#ifdef mjo_subdd
     *    ,phi,sda
#endif
      USE RAD_COM, only : trhr,fsf,cosz1,trsurf
#ifdef TRACERS_ON
      USE TRACER_COM, only : NTM,itime_tr0,needtrs,trm,trmom,
     *     n_CO2n, n_CFCn, n_Be7, n_Be10, n_clay, trname
#ifndef SKIP_TRACER_SRCS
     *     ,ntsurfsrc
#endif
#ifdef TRACERS_DRYDEP
     *     ,dodrydep
#endif
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE,trw0
#endif
#ifdef TRACERS_DUST
     &     ,Ntm_dust
#endif
#ifdef TRACERS_TOMAS
     &     ,xk,nbins,IDTSO4,IDTECIL
#endif
#endif
      USE PBLCOM, only : tsavg,dclev,eabl,uabl,vabl,tabl,qabl
      USE SOCPBL, only : npbl=>n
      USE PBL_DRV, only : alloc_pbl_args, dealloc_pbl_args
      USE PBL_DRV, only : pbl, t_pbl_args, xdelt
      USE DIAG_COM, only : MODD5S
      USE DIAG_COM, only : ndasf,ia_srf,ia_src,aij=>aij_loc,aijmm
     &     ,itocean,itoice,itlake,itlkice,itlandi
     *     ,tdiurn,adiurn=>adiurn_loc,ndiupt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_popocn
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,ij_tg1,ij_evap,ij_evapo,ij_tgo,ij_f0oc,ij_rhs
     *     ,ij_evapi,ij_f0li,ij_evapli,j_evap,j_evhdt,j_lwcorr
     *     ,j_tsrf,j_shdt,j_trhdt,j_type,j_tg1,j_tg2,ijdd,idd_spr
     *     ,idd_pt5,idd_ts,idd_tg1
     *     ,idd_q5,idd_qs,idd_qg,idd_swg
     *     ,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg,idd_wg,idd_us
     *     ,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq,idd_eds,idd_dbl
     *     ,idd_ev,idd_ldc,idd_dcf,ij_pblht,ndiuvar,NREG,ij_dskin
     *     ,ij_gusti,ij_mccon,ij_sss,ij_trsup,ij_trsdn,ij_fwoc,ij_ssh
     *     ,adiurn_dust, ij_kw, ij_alpha, ij_gasx, ij_silwu, ij_silwd
     *     ,ij_sish ,ij_popwat
     *     ,ij_tsurfmin,ij_tsurfmax
#if (defined mjo_subdd) || (defined etc_subdd)
     *     ,qsen_avg,qlat_avg,pblht_acc
#endif
#ifdef mjo_subdd
     *     ,PW_acc, E_acc,sst_avg,p_avg,lwu_avg
     *     ,u_avg,v_avg,w_avg,t_avg,q_avg,r_avg,z_avg
#endif
#ifndef NO_HDIURN
     *     ,hdiurn=>hdiurn_loc
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
     *     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf
#endif
#ifdef TRACERS_DUST
     *     ,idd_ws2,idd_ustar,idd_us3,idd_stress,idd_lmon
     *     ,idd_rifl,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,idd_tabl1
     *     ,idd_qabl1,idd_zhat1,idd_e1,idd_km1,idd_ri1,idd_grav,idd_turb
#endif
      USE LANDICE, only : z1e,z2li,hc1li,hc2li,ace1li,ace2li,snmin
      USE LANDICE_COM, only : snowli
#ifdef TRACERS_WATER
     *     ,trlndi
#endif
      USE SEAICE, only : xsi,ace1i,alami0,rhoi,byrls,solar_ice_frac
     *     ,tfrez,dEidTi,alami
      USE SEAICE_COM, only : si_atm
      USE LAKES_COM, only : mwl,gml,flake,icelak
      USE LAKES, only : minmld
#ifdef mjo_subdd
      USE GHY_COM, only : FEARTH
#endif
      USE FLUXES, only : dth1,dq1,runoe,erunoe
     *     ,nstype,uflux1,vflux1,tflux1,qflux1
     *     ,UOdrag
     &     ,nisurf,fland,flice,focean
     &     ,atmocn,atmice,atmgla,asflx
#ifdef TRACERS_ON
     *     ,trsrfflx
#ifndef SKIP_TRACER_SRCS
     *     ,trsource
#endif
#ifdef TRACERS_WATER
     *     ,trunoe
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     &     ,pprec,pevap,dust_flux_glob
#ifdef TRACERS_DRYDEP
     &     ,depo_turb_glob,depo_grav_glob
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_AMP) 
     &     ,dust_flux2_glob
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7D_acc
#endif
#ifndef SKIP_TRACER_DIAGS
      USE TRDIAG_COM, only : taijn=>taijn_loc,
     *      taijs=>taijs_loc,ijts_isrc, jls_isrc, tij_surf,
     *      tij_surfbv, tij_evap,ijts_gasex,
     *      tij_grnd, tij_drydep, tij_gsdep, ijts_Sdrydep
#ifdef TRACERS_DRYDEP
     *      , itcon_dd,itcon_surf
#endif
#ifdef BIOGENIC_EMISSIONS
     *     ,  ijs_isoprene
#endif
#endif /*SKIP_TRACER_DIAGS*/
#ifdef TRACERS_ON
      use trdiag_com, only: trcsurf,trcSurfByVol
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      USE tracers_dust, only : hbaij,ricntd,n_soildust
#endif
!#ifdef TRACERS_TOMAS
!      USE TOMAS_AEROSOL, ONLY : TOMAS_EMIS
!#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_COM, only: vol2mass,tr_mm
#ifdef OBIO_ON_GARYocean
      USE MODEL_COM, only: nstep=>itime
#else
      USE HYCOM_SCALARS, only: nstep
#endif
#endif
      USE SOIL_DRV, only: earth,ground_e

!@var DDMS downdraft mass flux in kg/(m^2 s), (i,j)
      USE CLOUDS_COM, only : DDMS
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     &     ,ddml,fss
#endif
      USE Timer_mod, only: Timer_type
      USE TimerList_mod, only: startTimer => start
      USE TimerList_mod, only: stopTimer => stop
      USE itype_enum

      IMPLICIT NONE
! ================ Parameter Declarations =====================
      integer, intent(in) :: ns,moddsf,moddd
      REAL*8, intent(inout), DIMENSION(
     *       NSTYPE,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &       GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *     TGRND,TGRN2,TGR4
      REAL*8, intent(inout), DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSTYPE) ::
     *     E1
! ================ VARIABLE DECLARATIONS ======================
      INTEGER I,J,K,KR,JR,NSTEPS,ITYPE,IH,IHM,IDTYPE
     *     ,ii
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1K
     *     ,ELHX,MSI2,CDTERM,CDENOM,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG
     *     ,dT2,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER
     *     ,PSK,Q1,THV1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,RCDHDWS,RCDQDWS
     *     ,SHEAT,TRHEAT,T2DEN,T2CON,T2MUL,FQEVAP,Z1BY6L,F2
     *     ,FSRI(2),dlwdt,byNIsurf,TGO
#ifdef mjo_subdd
     *     ,PEARTH
#endif

      REAL*8 MA1, MSI1
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: Z2LI3L=Z2LI/(3.*ALAMI0), Z1LIBYL=Z1E/ALAMI0
      REAL*8 QSAT,DQSATDT,TR4
c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 qg_sat,dtsurf,uocean,vocean,qsrf,us,vs,ws,ws0,
     &     dmua_ij,dmva_ij
c      logical pole
c
      logical :: lim_lake_evap,lim_dew ! for tracer convenience
#ifdef TRACERS_ON
      real*8 totflux(ntm)
      integer n,nx,nsrc
      real*8, dimension(ntm) :: trs,trsfac,trconstflx
      integer ntix(ntm), ntx
      real*8, dimension(ntm) :: trgrnd,trgrnd2
      real*8 trc_flux
#ifdef TRACERS_WATER
      real*8  TEV,tevap,dTQS,TDP,TDT1,frac
#ifdef TRACERS_SPECIAL_O18
     *     ,FRACVL,FRACVS
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      real*8 alpha_gas2
#endif      
#ifdef TRACERS_DRYDEP
      real*8 tdryd, tdd, td1, rtsdt, rts
#endif
#ifdef TRACERS_DUST
      INTEGER :: n1
#endif
#endif
C**** some shorthand indices and arrays for diurn diags
      INTEGER, PARAMETER :: n_idx1 = 11
      INTEGER, PARAMETER :: n_idx2 = 22
      INTEGER, PARAMETER :: n_idx3 = 2
      INTEGER, PARAMETER :: n_idx4 = n_idx1+n_idx2
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER,PARAMETER :: n_idxd=5
#else
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idxd=13+6*npbl+4*(npbl-1)
#endif
#endif

      INTEGER :: idx1(n_idx1), idx2(n_idx2), idx3(n_idx3)
      INTEGER :: idx4(n_idx1+n_idx2)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER :: idxd(n_idxd)
#endif
#ifdef TRACERS_TOMAS
      INTEGER ss_bin,num_bin
      real*8 ss_num(nbins),tot_seasalt
      real*8, parameter :: scalesizeSS(nbins)=(/!0.0,0.0,0.0,
     *     6.4614E-08,5.0110E-07,2.7243E-06,1.1172E-05,
     *     3.7192E-05,1.2231E-04,4.4986E-04,1.4821E-03,
     *     3.7403E-03,7.9307E-03,1.8918E-01,7.9705E-01/)
#endif
      REAL*8 :: tmp(NDIUVAR)
C****
      INTEGER :: J_0, J_1, J_0H, J_1H, I_0,I_1
      LOGICAL :: debug

      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,SSS
      REAL*8, DIMENSION(:,:,:), POINTER :: SSI
#ifdef TRACERS_GASEXCH_ocean
      REAL*8, DIMENSION(:,:,:), POINTER :: TRGASEX
#endif

      type (Timer_type), pointer :: aTimer
! ======================= MAIN BODY ======================================

      RSI => SI_ATM%RSI
      MSI => SI_ATM%MSI
      SNOWI => SI_ATM%SNOWI
      SSI => SI_ATM%SSI
      SSS => atmocn%SSS
#ifdef TRACERS_GASEXCH_ocean
      TRGASEX => atmocn%TRGASEX
#endif

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,        J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Initialise constant indices
      idx1 = (/ IDD_SPR, (IDD_PT5+ii-1,ii=1,5), (IDD_Q5+ii-1,ii=1,5) /)
      idx2 = (/ IDD_TS,  IDD_TG1, IDD_QS,  IDD_QG,  IDD_SWG,
     &          IDD_LWG, IDD_SH,  IDD_LH,  IDD_HZ0, IDD_UG,
     &          IDD_VG,  IDD_WG,  IDD_US,  IDD_VS,  IDD_WS,
     &          IDD_CIA, IDD_CM,  IDD_CH,  IDD_CQ,  IDD_EDS,
     &          IDD_DBL, IDD_EV /)
      idx3 = (/ IDD_DCF, IDD_LDC /)
      idx4(:n_idx1)   = idx1
      idx4(n_idx1+1:) = idx2

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IF (adiurn_dust == 1) THEN
        idxd=(/idd_wtke, idd_wd, idd_wm, idd_wsgcm, idd_wspdf
#ifdef TRACERS_DUST
     *       ,idd_grav, idd_turb, idd_ws2, idd_ustar, idd_us3
     *       ,idd_stress, idd_lmon, idd_rifl,
     *       (idd_zpbl1+ii-1,ii=1,npbl), (idd_uabl1+ii-1,ii=1,npbl),
     *       (idd_vabl1+ii-1,ii=1,npbl), (idd_uvabl1+ii-1,ii=1,npbl),
     *       (idd_tabl1+ii-1,ii=1,npbl), (idd_qabl1+ii-1,ii=1,npbl),
     *       (idd_zhat1+ii-1,ii=1,npbl-1), (idd_e1+ii-1,ii=1,npbl-1),
     *       (idd_km1+ii-1,ii=1,npbl-1), (idd_ri1+ii-1,ii=1,npbl-1)
#endif
     &       /)
      END IF
#endif

C****

      call startTimer('SURFACE_LANDICE()')

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      byNIsurf=1.d0/real(NIsurf)
      IH=JHOUR+1
      IHM = IH+(JDATE-1)*24
c avoid uninitialized variable problems when the first gridpoint
c in the domain is ocean
      SNOW = 0.

!! /============================================================\
!! THIS CODE NEEDS to go in SURFACE_LANDICE.f!!!!
!C**** INITIALIZE TGRND: THIS IS USED TO UPDATE T OVER SURFACE STEPS
!! Because we're doing landice, just initialize TGRND(3,...)
!      DO J=J_0,J_1
!      DO I=I_0,I_1
!!        TGRND(2,I,J)=atmice%GTEMP(I,J)
!        TGRND(3,I,J)=atmgla%GTEMP(I,J)
!!        TGRN2(2,I,J)=atmice%GTEMP2(I,J)
!        TGRN2(3,I,J)=atmgla%GTEMP2(I,J)
!!        TGR4(2,I,J)=atmice%GTEMPR(I,J)**4
!        TGR4(3,I,J)=atmgla%GTEMPR(I,J)**4
!      END DO 
!      END DO 


C**** Zero out fluxes summed over type and surface time step
c ---- Don't need this section, it's asflx is a global variable
c and it's already zeroed out before outer loop
!       do itype=1,4
!         asflx(itype)%solar(:,:) = 0.
!         asflx(itype)%dmua(:,:) = 0.
!         asflx(itype)%dmva(:,:) = 0.
!         asflx(itype)%e0(:,:) = 0.
!         asflx(itype)%evapor(:,:) = 0.
! #ifdef TRACERS_WATER
!         asflx(itype)%TREVAPOR(:,:,:) = 0.
! #endif
! #ifdef TRACERS_DRYDEP
!         asflx(itype)%TRDRYDEP(:,:,:) = 0.
! #endif
!       enddo

!      E1=0.    ! I don't think local var E1 is really used.  It is set, but not read.
! /------------ All global variables, zeroed in SURFACE.f ------------\
!      RUNOE=0.   ! From FLUXES, zeroed in SURFACE.f
!      ERUNOE=0.! From FLUXES, zeroed in SURFACE.f
!
!#ifdef TRACERS_WATER
!      TRUNOE = 0.
!#endif
!#ifdef SCM
!      EVPFLX= 0.0d0
!      SHFLX = 0.0d0
!#endif
! #ifdef TRACERS_GASEXCH_ocean
!       TRGASEX = 0.0d0
! #endif
! #if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
!     (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) 
!       dust_flux_glob = 0.d0
!       dust_flux2_glob = 0.d0
!       depo_turb_glob = 0.d0
!       depo_grav_glob = 0.d0
! #endif
! \--------------------------------------------------------------------/
! Next in SURFACE.f comes outside loop over timesteps.

! ==============================================================


#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          nx=nx+1
          ntix(nx) = n
        end if
      end do 
      ntx = nx
#endif

! I don't think this is needed because it's not in GHY_DRV.f
!      call recalc_agrid_uv


      call alloc_pbl_args(pbl_args)
C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

	  ! Model may expect that some things are zero
      ! call stop_model('Please double-check that model might expect some things to be defined for some boxes, even if there is no ice in that box.', -17)
      PLICE=FLICE(I,J)
      PTYPE = PLICE
      IF (PTYPE <= 0) CYCLE


c      debug=i.eq.65.and.j.eq.38

C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
C     RSI   RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
      POICE=RSI(I,J)*PWATER
      POCEAN=PWATER-POICE
C     Not Used.  P == pressure?
      PIJ=P(I,J)
c     ./model/ATMDYN_COM.f:!@var PEDN edge pressure (top of box) (mb)
      PS=PEDN(1,I,J)
c !@var  PEK  PEDN**KAPA
      PSK=PEK(1,I,J)
c !@var  PK   PMID**KAPA
      P1K=PK(1,I,J)
C !@var Q specific humidity (kg water vapor/kg air)
      Q1=Q(I,J,1)
C T=Temperature

      THV1=T(I,J,1)*(1.+Q1*xdelt)
      JR=JREG(I,J)  !@var JREG lat/lon array defining regions for AREG diagnostics
      MA1=AM(1,I,J) !@var MA1 mass of lowest atmospheric layer (kg/m^2)


#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
C**** Calculate first layer tracer concentration
          pbl_args%trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*byaxyp(i,j)
        end if
      end do
#endif

!C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
!         IF(MODDD.EQ.0) THEN
!         DO KR=1,NDIUPT
!           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
!             tmp(IDD_SPR)=PS
!             do ii=1,5
!               tmp(IDD_PT5+ii-1)=PSK*T(I,J,ii)
!               tmp(IDD_Q5+ii-1) =Q(I,J,ii)
!             end do
!             ADIURN(idx1(:),kr,ih)=ADIURN(idx1(:),kr,ih)+tmp(idx1(:))
!#ifndef NO_HDIURN
!             HDIURN(idx1(:),kr,ihm)=HDIURN(idx1(:),kr,ihm)+tmp(idx1(:))
!#endif
!           END IF
!         END DO
!         END IF
!C**** save some ocean diags regardless of PTYPE
!C**** SSH does not work for qflux/fixed SST configurations
!         if(ns.eq.1 .and. focean(i,j).gt.0.)
!     &        aij(i,j,ij_popocn) = aij(i,j,ij_popocn) + pocean
!         IF (FOCEAN(I,J).gt.0. .and. MODDSF.eq.0) THEN
!           AIJ(I,J,IJ_TGO)=AIJ(I,J,IJ_TGO)+atmocn%GTEMP(I,J)*FOCEAN(I,J)
!           AIJ(I,J,IJ_SSS)=AIJ(I,J,IJ_SSS)+SSS(I,J)*FOCEAN(I,J)
!           AIJ(I,J,IJ_SSH)=AIJ(I,J,IJ_SSH)+(atmocn%OGEOZA(I,J)*BYGRAV+
!     *          RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS)*FOCEAN(I,J)
!         END IF
C****
      DO ITYPE=ITYPE_LANDICE, ITYPE_LANDICE    ! no earth type
  !    ipbl(itype,i,j)=0


      ! This is a good way to make sure you dont forget to do something
      ! before you run the model.
      ! call stop_model('Please double-check something or another.', -17)

! BEGIN ---------------------------------------------------------
      PTYPE=PLICE
      IF (PTYPE.gt.0) THEN
          IDTYPE=ITLANDI
          SNOW=SNOWLI(I,J)
          TG1=TGRND(3,I,J)
          TG2=atmgla%GTEMP2(I,J)
          TR4=TGR4(3,I,J)
          SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
          Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
          CDTERM=TG2
          CDENOM=1./(2.*Z1BY6L+Z2LI3L)
          HCG1=HC1LI+SNOW*SHI
          ELHX=LHS
          uocean = 0. ; vocean = 0. ! no land ice velocity
#ifdef TRACERS_WATER
              do nx=1,ntx
                trgrnd2(nx)=TRLNDI(ntix(nx),I,J)/(ACE1LI+ACE2LI)
              end do
#endif
      END IF
! END ---------------------------------------------------------

! ----------------------------------------------------------
! This code is replicated when factorut of SURFACE.f
! Bob should not have to touch this code, it is not LANDICE-specific.

C****
      IF (PTYPE.gt.0) THEN
C****
C**** BOUNDARY LAYER INTERACTION
C****
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.

      TG=TG1+TF
      QG_SAT=QSAT(TG,ELHX,PS)
!      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG_SAT=0.98d0*QG_SAT
      pbl_args%TG=TG   ! actual ground temperature
      pbl_args%TR4=TR4 ! radiative temperature K^4
      pbl_args%ELHX=ELHX   ! relevant latent heat
      pbl_args%QSOL=SRHEAT   ! solar heating
      pbl_args%TGV=TG*(1.+QG_SAT*xdelt)

         write (6,*) 'PBL_DEBUG3',i,j,
     &      TG,TG1,TF
!         write (6,*) 'PBL_DEBUG2',i,j,
!     &      pbl_args%TGV,TG,QG_SAT,xdelt,ELHX,PS,QSAT(TG,ELHX,PS)

#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
C**** set defaults
        trsfac(nx)=0.
        totflux(nx)=0.
        trconstflx(nx)=0.
!#ifdef TRACERS_GASEXCH_ocean
!       IF (ITYPE.EQ.1 .and. focean(i,j).gt.0.) THEN  ! OCEAN
!         trgrnd(nx)=atmocn%gtracer(n,i,j)    
!         trsfac(nx)=1.
!         trconstflx(nx)=trgrnd(nx)
!         if (i.eq.1.and.j.eq.45)
!     .        write(*,'(a,3i5,3e12.4)') 'in SURFACE, gtracer:',
!     .        nstep,I,J,sss(I,J),atmocn%gtracer(n,i,j),trgrnd(nx)
!       END IF
!#endif
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
#ifdef TRACERS_WATER
        pbl_args%tr_evap_max(nx)=1.
C**** This distinguishes water from gases or particle
        if ( tr_wd_TYPE(n) == nWATER ) then
          trgrnd(nx)=asflx(itype)%gtracer(n,i,j)
C**** trsfac and trconstflx are multiplied by cq*ws and QG_SAT in PBL
          trsfac(nx)=1.
          trconstflx(nx)=trgrnd(nx)

        else if ( tr_wd_TYPE(n) == nGAS .or.
     &         tr_wd_TYPE(n) == nPART ) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            trsfac(nx)=1.       !then multiplied by deposition velocity in PBL
#ifdef TRACERS_WATER
            pbl_args%tr_evap_max(nx)=1.d30
            trgrnd(nx)=0.
#endif
          end if
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
C**** Now send kg/m^2/s to PBL, and divided by rho there.
#ifndef SKIP_TRACER_SRCS
          do nsrc=1,ntsurfsrc(n)
            totflux(nx) = totflux(nx)+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux(nx)*byaxyp(i,j)   ! kg/m^2/s
#endif /*SKIP_TRACER_SRCS*/

#ifdef TRACERS_WATER
        endif
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
        !need to redo this here because the previous line has changed trconstflx to zero.
        !because we have no sources. is there a better way to do this?
        trconstflx(nx)=trgrnd(nx) * byaxyp(i,j)   !kg,co2/kg,air/m2
#endif
      end do
#endif
C =====================================================================
      pbl_args%dtsurf = dtsurf
      pbl_args%TKV=THV1*PSK     ! TKV is referenced to the surface pressure
      pbl_args%ZS1=.5d-2*RGAS*pbl_args%TKV*MA1/PMID(1,I,J)
      pbl_args%qg_sat = qg_sat
      pbl_args%qg_aver = qg_sat   ! QG_AVER=QG_SAT
      pbl_args%hemi = sign(1d0,lat2d(i,j))
c      pbl_args%pole = pole
      pbl_args%evap_max = 1.
      pbl_args%fr_sat = 1. ! entire surface is saturated
      pbl_args%uocean = uocean
      pbl_args%vocean = vocean
      pbl_args%psurf = PS
      pbl_args%trhr0 = TRHR(0,I,J)
      pbl_args%ocean = (ITYPE.eq.1 .and. FOCEAN(I,J).gt.0)    !false
      pbl_args%snow = SNOW
#ifdef TRACERS_ON
      pbl_args%trs(1:ntm) = trs(1:ntm)
      pbl_args%trsfac(1:ntm) = trsfac(1:ntm)
      pbl_args%trconstflx(1:ntm) = trconstflx(1:ntm)
      pbl_args%trgrnd2(1:ntm) = trgrnd2(1:ntm)
      pbl_args%ntix(1:ntm) = ntix(1:ntm)
      pbl_args%ntx = ntx
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) 
      pbl_args%hbaij=hbaij(i,j)
      pbl_args%ricntd=ricntd(i,j)
      pbl_args%pprec=pprec(i,j)
      pbl_args%pevap=pevap(i,j,itype)
c**** fractional area of moist convection * fraction of downdrafts
c**** (=1/3), but only if downdrafts reach the lowest atmospheric
c**** layer. It's only needed to constrain dust emission due to
c**** downdrafts for soil type earth, but it's needed here to calculate
c**** wspdf in PBL.f for the other soil types.
      IF (ddml(i,j) == 1) THEN
        pbl_args%mcfrac=(1.D0-fss(1,i,j))*By3
      ELSE
        pbl_args%mcfrac=0.D0
      END IF
#endif

! Calculate drag coefficients, wind speed, air density, etc.
! PBL = "Planetary Boundary Layer"
C**** Call pbl to calculate near surface profile
      CALL PBL(I,J,ITYPE,PTYPE,pbl_args)

#ifdef TRACERS_ON
      trs(1:ntm) = pbl_args%trs(1:ntm)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) 
      hbaij(i,j)=pbl_args%hbaij
      ricntd(i,j)=pbl_args%ricntd
#endif
      us = pbl_args%us
      vs = pbl_args%vs
      ws = pbl_args%ws
      ws0 = pbl_args%ws0
      qsrf = pbl_args%qsrf
      CM = pbl_args%cm
      CH = pbl_args%ch
      CQ = pbl_args%cq
      TS=pbl_args%TSV/(1.+QSRF*xdelt)
C =====================================================================

C**** Adjust ground variables to account for skin effects
      TG = TG + pbl_args%dskin
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (pbl_args%ocean) QG_SAT=0.98d0*QG_SAT
      TG1 = TG - TF
      TR4=(sqrt(sqrt(TR4))+pbl_args%dskin)**4
! %dskin trapped here!!!

C**** CALCULATE RHOSRF*CM*WS AND RHOSRF*CH*WS
      RHOSRF=100.*PS/(RGAS*pbl_args%TSV)
      RCDMWS=CM*WS*RHOSRF
      RCDHWS=CH*WS*RHOSRF
      RCDQWS=CQ*WS*RHOSRF
      RCDHDWS=CH*(WS-WS0)*RHOSRF
      RCDQDWS=CQ*(WS-WS0)*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2) (positive down)
      ! Including gustiness in the sensible heat flux:
      SHEAT=SHA*(RCDHWS*(TS-TG)+RCDHDWS*pbl_args%tprime)
      ! Including gustiness in the latent heat flux:
      EVHEAT=(LHE+TG1*SHV)*(RCDQWS*(QSRF-QG_SAT)+
     *                      RCDQDWS*pbl_args%qprime)
      TRHEAT=TRHR(0,I,J)-STBO*TR4

! BEGIN ------------------------------------------------------------
! Setting up implicit timestep aspects to surface flux calculations
! Setting up derivitives of fluxes w.r.t temperature, etc.
! Might need to change these if I decide to use layering, differnt
! variable names, etc.
C**** CASE (3) ! FLUXES USING IMPLICIT TIME STEP OVER LANDICE
      if ( ITYPE == ITYPE_LANDICE ) then   ! true

        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
        DSHDTG=-RCDHWS*SHA
        DQGDTG=QG_SAT*DQSATDT(TG,ELHX)
        DEVDTG=-RCDQWS*LHE*DQGDTG
        DTRDTG=-4.*STBO*sqrt(sqrt(TR4))**3
        DF0DTG=DSHDTG+DEVDTG+DTRDTG
        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)
        SHDT=DTSURF*(SHEAT+DTG*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+DTG*DEVDTG)
        TRHDT=DTSURF*(TRHEAT+DTG*DTRDTG)
        F1DT=DTSURF*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
        TG1=TG1+DTG
      endif
! END ------------------------------------------------------------

C**** CALCULATE EVAPORATION
      DQ1X =EVHDT/((LHE+TG1*SHV)*MA1)
      EVHDT0=EVHDT
C**** Limit evaporation if lake mass is at minimum
      lim_lake_evap=.false.
      lim_dew=.false.


      IF (DQ1X.GT.Q1+DQ1(I,J)) THEN
          DQ1X=(Q1+DQ1(I,J))
          lim_dew=.true.
      ELSE
          GO TO 3720
      END IF
      EVHDT=DQ1X*(LHE+TG1*SHV)*MA1
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*MA1

! =============== BEGIN Loop over Tracers
#ifdef TRACERS_ON
#ifdef TRACERS_TOMAS
! Aerosol Physics
      ss_bin=0
      num_bin=0
!      TOMAS_emis(I,J,:,1)=0.
#endif
C**** Loop over tracers
      DO NX=1,NTX
        N=NTIX(NX)
#ifdef TRACERS_WATER
        if (tr_wd_TYPE(n).eq.nWATER) THEN
C****
C**** Calculate Water Tracer Evaporation
C****
          IF (ITYPE.EQ.1) THEN  ! OCEAN
          ELSE                  ! ICE AND LAND ICE
C**** tracer flux is set by source tracer concentration
            IF (EVAP.GE.0) THEN ! EVAPORATION
              IF (EVAP.le.SNOW .or. SNOW.lt.SNMIN .or. ITYPE.ne.ITYPE_LANDICE) THEN
                TEVAP=EVAP*trgrnd(nx)
              ELSE ! special treatment for landice when EVAP>SNOW>SNMIN
                TEVAP=SNOW*(trgrnd(nx)-trgrnd2(nx))+EVAP*trgrnd2(nx)
              END IF
            ELSE                ! DEW (fractionates)
#ifdef TRACERS_SPECIAL_O18
              IF (TG1.gt.0) THEN
                frac=FRACVL(TG1,n)
              ELSE
                frac=FRACVS(TG1,n)
              END IF
              TEVAP=EVAP*trs(nx)/(QSRF*frac+teeny)
#else
              TEVAP=EVAP*trs(nx)/(QSRF+teeny)
#endif       // #ifdef TRACERS_SPECIAL_O18
            END IF
          END IF
          TDP = TEVAP*AXYP(I,J)*ptype
          TDT1 = trsrfflx(I,J,n)*DTSURF
#ifdef WATER_PROPORTIONAL
          if(lim_dew) then
#else
          IF (TRM(I,J,1,n)+TDT1+TDP.lt.0..and.tdp.lt.0) THEN
#endif
            IF (QCHECK) WRITE(99,*) "LIMITING TRDEW",I,J,N,TDP,TRM(I,J,1
     *           ,n),TDT1
            TEVAP = -(TRM(I,J,1,n)+TDT1)/(AXYP(I,J)*ptype)
            trsrfflx(I,J,n)= - TRM(I,J,1,n)/DTSURF
          ELSE
            trsrfflx(I,J,n)=trsrfflx(I,J,n)+TDP/DTSURF
          END IF
          !TREVAPOR(n,ITYPE,I,J)=TREVAPOR(n,ITYPE,I,J)+TEVAP
          asflx(itype)%TREVAPOR(n,I,J)=
     *         asflx(itype)%TREVAPOR(n,I,J)+TEVAP
        END IF
#endif // #ifdef TRACERS_WATER
!        end if     ! if (IType.ne.1)


#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
C****
C**** Calculate Aersosol Exchange
C****
        select case (trname(n))
        case ('DMS')
          trc_flux=pbl_args%DMS_flux
        case ('seasalt1', 'M_SSA_SS')
          trc_flux=pbl_args%ss1_flux
        case ('seasalt2', 'M_SSC_SS')
          trc_flux=pbl_args%ss2_flux
        case ('M_SSS_SS')
          trc_flux=(pbl_args%ss1_flux+pbl_args%ss2_flux)
#ifdef TRACERS_AEROSOLS_OCEAN
        case ('OCocean')
          trc_flux=pbl_args%OCocean_flux
#endif  /* TRACERS_AEROSOLS_OCEAN */
        case default
          trc_flux=0
#ifdef TRACERS_TOMAS
        case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04', 
     &         'ANACL_05','ANACL_06','ANACL_07','ANACL_08',
     &         'ANACL_09','ANACL_10','ANACL_11','ANACL_12')

          ss_bin=ss_bin+1
          if(ss_bin.eq.1)
     &     tot_seasalt=(pbl_args%ss2_flux + pbl_args%ss1_flux)

          trc_flux=tot_seasalt*scalesizeSS(ss_bin)
          ss_num(ss_bin)=tot_seasalt*scalesizeSS(ss_bin)
     &         /sqrt(xk(ss_bin)*xk(ss_bin+1))

! No subgrid coagulation for sea-salt
!        TOMAS_EMIS(I,J,ss_bin,1)= trc_flux*axyp(i,j)*ptype

        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04',
     &         'ANUM__05','ANUM__06','ANUM__07','ANUM__08',
     &         'ANUM__09','ANUM__10','ANUM__11','ANUM__12')
           num_bin=num_bin+1
           trc_flux=ss_num(num_bin)
      
#endif
        end select

        trsrfflx(i,j,n)=trsrfflx(i,j,n)+
     &       trc_flux*axyp(i,j)*ptype
        if (ijts_isrc(1,n)>0) then
           taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &       trc_flux*axyp(i,j)*ptype*dtsurf
        end if

#ifdef TRACERS_AMP
            select case (trname(n))
              case ('DMS','M_SSA_SS','M_SSC_SS','M_SSS_SS')
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)
            end select
#else
#ifndef TRACERS_TOMAS

        if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols?
#endif    // #ifndef TRACERS_TOMAS
#endif

      
#ifdef TRACERS_TOMAS

        
            select case (trname(n))

            case ('DMS')              
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *              trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)

        case ('ANACL_01','ANACL_02','ANACL_03','ANACL_04', 
     &       'ANACL_05','ANACL_06','ANACL_07','ANACL_08',
     &       'ANACL_09','ANACL_10','ANACL_11','ANACL_12')
        
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)
        
        if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols?
        
        case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04',
     &       'ANUM__05','ANUM__06','ANUM__07','ANUM__08',
     &       'ANUM__09','ANUM__10','ANUM__11','ANUM__12')
        
!TOMAS - itcon_surf (1,3) is for SO4/EC/OC.
        if (itcon_surf(4,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(4,n),n)

        if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols? 

            end select
#endif
#endif
! We're missing and end if above this point.
!      end do        !**** This is for test

#ifdef TRACERS_DRYDEP
C****
C**** Calculate Tracer Dry Deposition (including gravitational settling)
C****
        if(dodrydep(n))then
          rts=rhosrf*trs(nx)
          rtsdt=rts*dtsurf                             ! kg*s/m^3
          tdryd=-rtsdt*(pbl_args%dep_vel(n)+pbl_args%gs_vel(n))          ! kg/m2
          tdd = tdryd*axyp(i,j)*ptype                    ! kg
          td1 = (trsrfflx(i,j,n)+totflux(nx))*dtsurf   ! kg
          if (trm(i,j,1,n)+td1+tdd.le.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd surface",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,trs(nx),pbl_args%trtop(nx)
            tdd= -max(trm(i,j,1,n)+td1,0d0)
            tdryd=tdd/(axyp(i,j)*ptype)
            trsrfflx(i,j,n)= - trm(i,j,1,n)/dtsurf
          else
            trsrfflx(i,j,n)=trsrfflx(i,j,n)+tdd/dtsurf
          end if
! trdrydep downward flux by surface type (kg/m^2)
          asflx(itype)%trdrydep(n,i,j)=
     &         asflx(itype)%trdrydep(n,i,j) - tdryd
! diagnose turbulent and settling fluxes separately
          taijn(i,j,tij_drydep,n)=taijn(i,j,tij_drydep,n) +
     &         ptype*rtsdt*pbl_args%dep_vel(n)
          taijn(i,j,tij_gsdep ,n)=taijn(i,j,tij_gsdep ,n) +
     &         ptype*rtsdt* pbl_args%gs_vel(n)
#ifdef ACCMIP_LIKE_DIAGS
! estimate stomatal tracer flux:
          if(trname(n)=='Ox')
     &    taijs(i,j,ijts_Sdrydep)=taijs(i,j,ijts_Sdrydep)+ptype*
     &         rtsdt*(pbl_args%stomatal_dep_vel)

#endif
#ifdef TRACERS_COSMO
          if (n .eq. n_Be7) BE7D_acc(i,j)=BE7D_acc(i,j)+ptype*rtsdt
     *         *pbl_args%dep_vel(n)+ptype*rtsdt* pbl_args%gs_vel(n)
#endif

          if (itcon_dd(n,1).gt.0) call inc_diagtcb(i,j,-
     &     ptype*rtsdt*axyp(i,j)*pbl_args%dep_vel(n),itcon_dd(n,1),n)
          if (itcon_dd(n,2).gt.0) call inc_diagtcb(i,j,-
     &     ptype*rtsdt*axyp(i,j)*pbl_args%gs_vel(n),itcon_dd(n,2),n)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
c**** for subdaily diagnostics
          depo_turb_glob( i, j, itype, n ) = depo_turb_glob( i, j, itype
     &         , n ) + ptype * rts * pbl_args%dep_vel( n ) / nisurf
          depo_grav_glob( i, j, itype, n ) = depo_grav_glob( i, j, itype
     &         , n ) + ptype * rts * pbl_args%gs_vel( n ) / nisurf
#endif

        end if     ! if dodrydep(n)

#endif
      END DO       ! do NX=1,NTX
#endif
! =============== END Loop over Tracers



C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
C**** Limit heat fluxes out of lakes if near minimum depth
      !E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      asflx(itype)%E0(I,J)=asflx(itype)%E0(I,J)+F0DT
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT
      !EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP
      asflx(itype)%EVAPOR(I,J)=asflx(itype)%EVAPOR(I,J)+EVAP
#ifdef SCM
      if (J.eq.J_TARG.and.I.eq.I_TARG) then
          if (SCM_SURFACE_FLAG.eq.0.or.SCM_SURFACE_FLAG.eq.2) then
              EVPFLX = EVPFLX -(DQ1X*MA1)*(PTYPE/DTSURF)*LHE
              SHFLX = SHFLX - SHDT*PTYPE/DTSURF
c             write(iu_scm_prt,*) 'srf  evpflx shflx ptype ',
c    *                   EVPFLX,SHFLX,ptype
          endif
      endif
#endif
      TGRND(ITYPE,I,J)=TG1  ! includes skin effects
      TGR4(ITYPE,I,J) =TR4
C**** calculate correction for different TG in radiation and surface
      dLWDT = DTSURF*(TRSURF(ITYPE,I,J)-TRHR(0,I,J))+TRHDT
C**** final fluxes
#ifdef SCM
cccccc for SCM use ARM provided fluxes for designated box
      if ((I.eq.I_TARG.and.J.eq.J_TARG).and.SCM_SURFACE_FLAG.eq.1) then
           DTH1(I,J)=DTH1(I,J)
     &              +ash*DTSURF*ptype/(SHA*MA1*P1K)
           DQ1(I,J)=DQ1(I,J) + ALH*DTSURF*ptype/(MA1*LHE)
           SHFLX = SHFLX + ASH*ptype
           EVPFLX = EVPFLX + ALH*ptype
           write(iu_scm_prt,980) I,PTYPE,DTH1(I,J),DQ1(I,J),
     &           EVPFLX,SHFLX
 980       format(1x,'SURFACE ARM   I PTYPE DTH1 DQ1 evpflx shflx',
     &            i5,f9.4,f9.5,f9.6,f9.5,f9.5)
      else
#endif
      DTH1(I,J)=DTH1(I,J)-(SHDT+dLWDT)*PTYPE/(SHA*MA1*P1K) ! +ve up
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
#ifdef SCM
      if (i.eq.I_TARG.and.j.eq.J_TARG) then
          write(iu_scm_prt,988) I,PTYPE,DTH1(I,J),DQ1(I,J),SHDT,dLWDT
 988      format(1x,'988 SURFACE GCM  I PTYPE DTH1 DQ1 SHDT dLWDT ',
     &           i5,f9.4,f9.5,f9.6,f12.4,f10.4)
      endif
      endif
#endif
      DMUA_IJ=PTYPE*DTSURF*RCDMWS*(US-UOCEAN)
      DMVA_IJ=PTYPE*DTSURF*RCDMWS*(VS-VOCEAN)
      asflx(itype)%DMUA(I,J) = asflx(itype)%DMUA(I,J) + DMUA_IJ
      asflx(itype)%DMVA(I,J) = asflx(itype)%DMVA(I,J) + DMVA_IJ
      uflux1(i,j)=uflux1(i,j)+PTYPE*RCDMWS*(US-UOCEAN)
      vflux1(i,j)=vflux1(i,j)+PTYPE*RCDMWS*(VS-VOCEAN)
C****
C**** ACCUMULATE DIAGNOSTICS FOR EACH SURFACE TIME STEP AND ITYPE
C****
        CALL INC_AJ(I,J,IDTYPE,J_EVAP , EVAP*PTYPE)
        CALL INC_AJ(I,J,IDTYPE,J_EVHDT,EVHDT*PTYPE)
        CALL INC_AJ(I,J,IDTYPE,J_SHDT , SHDT*PTYPE)
        CALL INC_AJ(I,J,IDTYPE,J_TRHDT,TRHDT*PTYPE)
        CALL INC_AJ(I,J,IDTYPE,J_LWCORR,dLWDT*PTYPE)
        IF(MODDSF.EQ.0) THEN
          CALL INC_AJ(I,J,IDTYPE,J_TSRF,(TS-TF)*PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_TYPE,        PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_TG1 ,    TG1*PTYPE)
          CALL INC_AJ(I,J,IDTYPE,J_TG2 ,    TG2*PTYPE)
        END IF
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
        CALL INC_AREG(I,J,JR,J_TRHDT,TRHDT*PTYPE)
        CALL INC_AREG(I,J,JR,J_SHDT,  SHDT*PTYPE)
        CALL INC_AREG(I,J,JR,J_LWCORR,dLWDT*PTYPE)
        CALL INC_AREG(I,J,JR,J_EVHDT,EVHDT*PTYPE)
        CALL INC_AREG(I,J,JR,J_EVAP ,EVAP *PTYPE)
        IF(MODDSF.EQ.0) THEN
          CALL INC_AREG(I,J,JR,J_TSRF,(TS-TF)*PTYPE)
          CALL INC_AREG(I,J,JR,J_TG1 ,    TG1*PTYPE)
          CALL INC_AREG(I,J,JR,J_TG2 ,    TG2*PTYPE)
        END IF

C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
        AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDT*PTYPE
        IF(MODDSF.EQ.0) THEN
          AIJ(I,J,IJ_TRSDN)=AIJ(I,J,IJ_TRSDN)+TRHR(0,I,J)*PTYPE
          AIJ(I,J,IJ_TRSUP)=AIJ(I,J,IJ_TRSUP)+(TRHR(0,I,J)-TRHDT/DTSURF)
     *         *PTYPE
        END IF
        AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHEAT*DTSURF+TRHDT)*PTYPE
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHEAT*DTSURF+TRHDT+SHDT
     *       +EVHDT)*PTYPE
        AIJ(I,J,IJ_EVAP)=AIJ(I,J,IJ_EVAP)+EVAP*PTYPE
#ifdef mjo_subdd
C**** SUBDD E_acc for evaporation *** 
        E_acc(I,J)=E_acc(I,J)+EVAP*PTYPE
#endif
        IF(MODDSF.EQ.0) THEN
          AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WS*PTYPE
          AIJ(I,J,IJ_TS)=AIJ(I,J,IJ_TS)+(TS-TF)*PTYPE
          AIJ(I,J,IJ_US)=AIJ(I,J,IJ_US)+US*PTYPE
          AIJ(I,J,IJ_VS)=AIJ(I,J,IJ_VS)+VS*PTYPE
          AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RCDMWS*WS*PTYPE
          AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RCDMWS*(US-UOCEAN)*PTYPE
          AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RCDMWS*(VS-VOCEAN)*PTYPE
          AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QSRF*PTYPE
          AIJ(I,J,IJ_RHs)=AIJ(I,J,IJ_RHs)+QSRF*PTYPE/qsat(ts,elhx,ps)
          AIJ(I,J,IJ_TG1)=AIJ(I,J,IJ_TG1)+TG1*PTYPE
          AIJ(I,J,IJ_PBLHT)=AIJ(I,J,IJ_PBLHT)+pbl_args%dbl*PTYPE
#if (defined mjo_subdd) || (defined etc_subdd)
C**** SUBDD qblht_acc for PBL height *** YH Chen ***
          pblht_acc(I,J)=pblht_acc(I,J)+pbl_args%dbl*PTYPE
#endif
          if(DDMS(I,J).lt.0.) ! ddms < 0 for down draft
     *         AIJ(I,J,ij_mccon)=AIJ(I,J,ij_mccon)+ptype
          AIJ(I,J,IJ_GUSTI)=AIJ(I,J,IJ_GUSTI)+pbl_args%gusti*PTYPE


#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          aij(i,j,ij_wsgcm)=aij(i,j,ij_wsgcm)+pbl_args%wsgcm*ptype
          aij(i,j,ij_wspdf)=aij(i,j,ij_wspdf)+pbl_args%wspdf*ptype
          aij(i,j,ij_wdry)=aij(i,j,ij_wdry)+pbl_args%wsubwd*ptype
          aij(i,j,ij_wtke)=aij(i,j,ij_wtke)+pbl_args%wsubtke*ptype
          aij(i,j,ij_wmoist)=aij(i,j,ij_wmoist)+pbl_args%wsubwm*ptype
#endif

        END IF

C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
        IF(MODDD.EQ.0) THEN
          DO KR=1,NDIUPT
            IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
              tmp(IDD_TS)=TS*PTYPE
              tmp(IDD_TG1)=(TG1+TF)*PTYPE
              tmp(IDD_QS)=QSRF*PTYPE
              tmp(IDD_QG)=QG_SAT*PTYPE
              tmp(IDD_SWG)=SRHEAT*DTSURF
     *             *PTYPE
              tmp(IDD_LWG)=TRHDT*PTYPE
              tmp(IDD_SH)=SHDT*PTYPE
              tmp(IDD_LH)=EVHDT*PTYPE
              tmp(IDD_HZ0)=
     *             +(SRHEAT*DTSURF+TRHDT+SHDT+EVHDT)*PTYPE
              tmp(IDD_UG)=pbl_args%UG*PTYPE
              tmp(IDD_VG)=pbl_args%VG*PTYPE
              tmp(IDD_WG)=pbl_args%WG*PTYPE
              tmp(IDD_US)=US*PTYPE
              tmp(IDD_VS)=VS*PTYPE
              tmp(IDD_WS)=WS*PTYPE
              tmp(IDD_CIA)=pbl_args%PSI*PTYPE
              tmp(IDD_CM)=CM*PTYPE
              tmp(IDD_CH)=CH*PTYPE
              tmp(IDD_CQ)=CQ*PTYPE
              tmp(IDD_EDS)=pbl_args%KHS*PTYPE
              tmp(IDD_DBL)=pbl_args%DBL*PTYPE
              tmp(IDD_EV)=EVAP*PTYPE

              ADIURN(idx2(:),kr,ih)=ADIURN(idx2(:),kr,ih)+tmp(idx2(:))
#ifndef NO_HDIURN
             HDIURN(idx2(:),kr,ihm)=HDIURN(idx2(:),kr,ihm)+tmp(idx2(:))
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
              IF (adiurn_dust == 1) THEN
                tmp(idd_wsgcm)=pbl_args%wsgcm*ptype
                tmp(idd_wspdf)=pbl_args%wspdf*ptype
                tmp(idd_wtke)=pbl_args%wsubtke*ptype
                tmp(idd_wd)=pbl_args%wsubwd*ptype
                tmp(idd_wm)=pbl_args%wsubwm*ptype
#ifdef TRACERS_DUST
                tmp(idd_turb)=0.D0
                tmp(idd_grav)=0.D0
#ifdef TRACERS_DRYDEP
                nx = 0
                do n = 1,ntm
                  if (itime_tr0( n ) <= itime .and. needtrs( n )) then
                    nx = nx + 1
                    if (dodrydep( n )) then
                      select case(trname( n ))
                      case('Clay', 'Silt1', 'Silt2', 'Silt3')
                        tmp( idd_turb ) = tmp( idd_turb ) + ptype
     &                       * rhosrf * trs( nx ) * pbl_args%dep_vel( n
     &                       )
                        tmp( idd_grav ) = tmp( idd_grav ) + ptype
     &                       * rhosrf * trs( nx ) * pbl_args%gs_vel( n )
                      end select
                    end if
                  end if
                end do
#endif
                tmp(idd_ws2)=ws*ws*ptype
                tmp(idd_ustar)=pbl_args%ustar*ptype
                tmp(idd_us3)=ptype*pbl_args%ustar**3
                tmp(idd_stress)=rcdmws*pbl_args%ws*ptype
                tmp(idd_lmon)=pbl_args%lmonin*ptype
                tmp(idd_rifl)=
     &               +ptype*grav*(ts-tg)*pbl_args%zgs/(ws*ws*tg)

                tmp(idd_zpbl1:idd_zpbl1+npbl-1)=ptype*pbl_args%z(1:npbl)
                tmp(idd_uabl1:idd_uabl1+npbl-1)=
     *               ptype*uabl(1:npbl,itype,i,j)
                tmp(idd_vabl1:idd_vabl1+npbl-1)=
     *               ptype*vabl(1:npbl,itype,i,j)
                tmp(idd_uvabl1:idd_uvabl1+npbl-1)=ptype*sqrt(
     *               uabl(1:npbl,itype,i,j)*uabl(1:npbl,itype,i,j)+
     *               vabl(1:npbl,itype,i,j)*vabl(1:npbl,itype,i,j))
                tmp(idd_tabl1:idd_tabl1+npbl-1)=
     *               ptype*tabl(1:npbl,itype,i,j)
                tmp(idd_qabl1:idd_qabl1+npbl-1)=
     *               ptype*qabl(1:npbl,itype,i,j)
                tmp(idd_zhat1:idd_zhat1+npbl-2)=ptype
     *               *pbl_args%zhat(1:npbl-1)
                tmp(idd_e1:idd_e1+npbl-2)=eabl(1:npbl-1,itype,i,j)*ptype
                tmp(idd_km1:idd_km1+npbl-2)=ptype*pbl_args%km(1:npbl-1)
                tmp(idd_ri1:idd_ri1+npbl-2)=ptype*pbl_args%gh(1:npbl-1)
     *               /(pbl_args%gm(1:npbl-1)+1d-20)
#endif
                ADIURN(idxd(:),kr,ih)=ADIURN(idxd(:),kr,ih)+tmp(idxd(:))
#ifndef NO_HDIURN
                HDIURN(idxd(:),kr,ihm)=HDIURN(idxd(:),kr,ihm)+
     &               tmp(idxd(:))
#endif
              END IF
#endif
            END IF
          END DO
        END IF
C****
#ifdef TRACERS_ON
#ifndef SKIP_TRACER_DIAGS
C**** Save surface tracer concentration whether calculated or not
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          if (needtrs(n)) then
            nx=nx+1
            taijn(i,j,tij_surf  ,n) = taijn(i,j,tij_surf  ,n)
     *           +trs(nx)*ptype
            taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)
     *           +trs(nx)*ptype*rhosrf
            trcsurf(i,j,n)=trcsurf(i,j,n)+trs(nx)*ptype*byNIsurf
            trcSurfByVol(i,j,n)=trcSurfByVol(i,j,n)+trs(nx)*ptype*rhosrf
     &           *byNIsurf
          else
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *byaxyp(i,j),0d0)*ptype
            taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *byaxyp(i,j),0d0)*ptype*rhosrf
            trcsurf(i,j,n)=trcsurf(i,j,n)+max((trm(i,j,1,n)-trmom(mz,i,j
     *           ,1,n))*byam(1,i,j)*byaxyp(i,j),0d0)*ptype*byNIsurf
            trcSurfByVol(i,j,n)=trcSurfByVol(i,j,n)+max((trm(i,j,1,n)
     &           -trmom(mz,i,j,1,n))*byam(1,i,j)*byaxyp(i,j),0d0)*ptype
     &           *rhosrf*byNIsurf
          end if

#ifdef TRACERS_GASEXCH_ocean
          if (n.eq.n_CO2n .or. n.eq.n_CFCn) then
          if (focean(i,j).gt.0) then

! original versions
             if (POCEAN.gt.0) then ! original coding
               if (MODDSF.EQ.0) THEN
                 AIJ(i,j,ij_kw) = AIJ(i,j,ij_kw)  
     .            + pbl_args%Kw_gas * focean(i,j) ! m/s
     .            * (1.d0 - RSI(i,j)) ! only over open water
                 AIJ(i,j,ij_alpha) = AIJ(i,j,ij_alpha) 
     .           + pbl_args%alpha_gas * focean(i,j) ! mol,CO2/m3/uatm
     .           * (1.d0 - RSI(i,j)) ! only over open water
               endif
               if(NS==NIsurf .and. itype==1) then
                 AIJ(i,j,ij_gasx) = AIJ(i,j,ij_gasx) 
     .                + TRGASEX(n,I,J) * focean(i,j)
     .                * 3600.*24.*365.    ! mol,CO2/m2/yr
     .                * (1.d0 - RSI(i,j)) ! only over open water
               endif
            end if
! tracer diag versions
            if (ITYPE.eq.1) then 
! gas exchange
            taijs(i,j,ijts_gasex(3,n)) = taijs(i,j,ijts_gasex(3,n)) 
     .           + pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs(nx) 
     .         - pbl_args%alpha_gas * trgrnd(nx) ) 
     .         * 1d6/vol2mass(n)
     .         * dtsurf/dtsrc   !in order to accumulate properly over time
     .         * ptype * syr          ! mol/m2/yr

! zonal mean diag accumulates kgCO2
            if (jls_isrc(1,n)>0) call inc_tajls(i,j,1,jls_isrc(1,n),
     *           - pbl_args%Kw_gas * ( pbl_args%beta_gas  * trs(nx) 
     .         - pbl_args%alpha_gas * trgrnd(nx) ) 
     .         * 1d6/vol2mass(n) * dtsurf  
     .         * ptype*tr_mm(n)*1d-3*axyp(i,j))

              if (MODDSF.EQ.0) THEN
! piston velocity
                taijs(i,j,ijts_gasex(1,n)) = taijs(i,j,ijts_gasex(1,n)) 
     .                + pbl_args%Kw_gas * pocean ! m/s only over open water
! solubility mol/m3/uatm
                taijs(i,j,ijts_gasex(2,n)) = taijs(i,j,ijts_gasex(2,n)) 
     .               + pbl_args%alpha_gas * focean(i,j) 
              endif

            elseif (POCEAN.eq.0) then  ! ITYPE=2 and all ice covered
! solubility mol/m3/uatm ice covered area
               if (MODDSF.EQ.0) taijs(i,j,ijts_gasex(2,n)) = taijs(i,j
     $              ,ijts_gasex(2,n))+ alpha_gas2(tgo,pbl_args%sss_loc)
     $              * focean(i,j) 
            endif                ! itype
          endif                  !focean
          endif                  !gasexch tracers
#endif
#ifdef TRACERS_WATER
          if (tr_wd_type(n).eq.nWater) then
            taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *           asflx(itype)%trevapor(n,i,j)*ptype
! ==== DEBUGGING: Remove these lines to get exact match in regression tests
            if (jls_isrc(1,n)>0) call inc_tajls2(i,j,1,jls_isrc(1,n),
     *           asflx(itype)%trevapor(n,i,j)*ptype)
            if (focean(i,j)>0 .and. jls_isrc(2,n)>0) call inc_tajls2
     *          (i,j,1,jls_isrc(2,n),asflx(itype)%trevapor(n,i,j)*ptype)
! ===== END DEBUGGING
          end if
          taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     *         asflx(itype)%gtracer(n,i,j)*ptype
#endif
        end if
      end do
#endif /* SKIP_TRACER_DIAGS */
#endif
C****
C**** SAVE SOME TYPE DEPENDENT FLUXES/DIAGNOSTICS
C****
!!!      CASE (1)  ! ocean
C****
!!!      CASE (3) ! land ice
      if ( ITYPE == ITYPE_LANDICE ) then
        IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
        IF (MODDSF.eq.0) AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(TS-TF)*PTYPE
        AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)+ SHDT*PTYPE
        AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)  +EVHDT*PTYPE
        AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)  +TRHDT*PTYPE
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)    + F0DT*PTYPE
        AIJ(I,J,IJ_EVAPLI)=AIJ(I,J,IJ_EVAPLI)+ EVAP*PTYPE
C****
      endif
C****
      END IF
      END DO  ! end of itype loop
      END DO  ! end of I loop
      END DO  ! end of J loop
! ============================================================
! Now we're outside the loop over grid points

      call dealloc_pbl_args(pbl_args)

      call stopTimer('SURFACE_LANDICE()')

      RETURN
C****
      END SUBROUTINE SURFACE_LANDICE
