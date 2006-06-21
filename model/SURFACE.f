#include "rundeck_opts.h"

      SUBROUTINE SURFCE
!@sum SURFCE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : rgas,lhm,lhe,lhs
     *     ,sha,tf,rhow,shv,shi,stbo,bygrav,by6
     *     ,deltx,teeny,rhows     ! ,byrt3
#ifdef TRACERS_DUST
     &     ,grav
#endif
      USE MODEL_COM, only : im,jm,dtsrc,nisurf,u,v,t,p,q
     *     ,idacc,dsig,ndasf,fland,flice,focean
     *     ,nday,modrd,itime,jhour,itocean
     *     ,itoice,itlake,itlkice,itlandi,qcheck,UOdrag,jdate
      USE DOMAIN_DECOMP, only : GRID, GET, CHECKSUM, HALO_UPDATE, SOUTH
      USE DOMAIN_DECOMP, only : NORTH
      USE DOMAIN_DECOMP, only : AM_I_ROOT, GLOBALSUM
      USE GEOM, only : dxyp,imaxj,bydxyp,idjj,idij,rapj,kmaxj,sinip
     *     ,cosip
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE DYNAMICS, only : pmid,pk,pedn,pek,am,byam
      USE RAD_COM, only : trhr,fsf,cosz1
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
#ifdef TRACERS_DRYDEP
     *     ,dodrydep
#endif
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE,trname,trw0
#endif
#ifdef TRACERS_DUST
     &     ,Ntm_dust,n_clay
#endif
#if (defined TRACERS_AEROSOLS_Koch) ||(defined TRACERS_AMP)
      USE AEROSOL_SOURCES, only: SHDTT
#endif
#endif
C**** Interface to PBL
  !    USE SOCPBL, only : ZS1,TGV,TKV,QG_SAT,QG_AVER,HEMI,DTSURF,POLE
  !   &     ,US,VS,WS,WSM,WSH,TSV,QSRF,PSI,DBL,KHS !,KQS !,PPBL ! ,KMS
  !   &     ,UG,VG,WG !,WINT

      USE PBLCOM, only : tsavg,dclev
#ifdef TRACERS_DUST
     &     ,eabl,uabl,vabl,tabl,qabl
#endif
      USE PBL_DRV, only : pbl, t_pbl_args
  !   &     ,evap_max,fr_sat,uocean,vocean,psurf,trhr0
#ifdef TRACERS_ON
     *     ,trtop,trs,trsfac,trconstflx,ntx,ntix
#ifdef TRACERS_WATER
     *     ,tr_evap_max
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,z,km,gh,gm,zhat,lmonin,wsubtke,wsubwd,wsubwm
#endif
#ifdef TRACERS_DRYDEP
     *     ,dep_vel,gs_vel
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     *     ,DMS_flux, ss1_flux, ss2_flux
#endif
#endif
      USE DIAG_COM, only : oa,aij=>aij_loc
     *     ,tdiurn,aj=>aj_loc,areg,adiurn,ndiupt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,ij_tg1,ij_evap,ij_evapo,ij_tgo,ij_f0oc
     *     ,ij_f0oi,ij_evapi,ij_f0li,ij_evapli,j_evap,j_evhdt
     *     ,j_tsrf,j_shdt,j_trhdt,j_type,j_tg1,j_tg2,ijdd,idd_spr
     *     ,idd_pt5,idd_pt4,idd_pt3,idd_pt2,idd_pt1,idd_ts,idd_tg1
     *     ,idd_q5,idd_q4,idd_q3,idd_q2,idd_q1,idd_qs,idd_qg,idd_swg
     *     ,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg,idd_wg,idd_us
     *     ,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq,idd_eds,idd_dbl
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     *     ,idd_ev,idd_ldc,idd_dcf,ij_pblht,ndiuvar,NREG
#else
     *     ,idd_ev,idd_ldc,idd_dcf,hdiurn,ij_pblht,ndiuvar,NREG
#endif
     &     ,ij_sss,ij_trsup,ij_trsdn,ij_fwoc,ij_ssh,adiurn_dust
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
     *     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf
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
     &     ,idd_grav,idd_turb
#endif
      USE LANDICE, only : z1e,z2li,hc1li,hc2li
      USE LANDICE_COM, only : snowli
      USE SEAICE, only : xsi,ace1i,alami,byrli,byrls, ! z1i,
     *     solar_ice_frac
      USE SEAICE_COM, only : rsi,msi,snowi,flag_dsws
      USE LAKES_COM, only : mwl,gml,flake
      USE LAKES, only : minmld
      USE FLUXES, only : dth1,dq1,e0,e1,evapor,runoe,erunoe,sss
     *     ,solar,dmua,dmva,gtemp,nstype,uflux1,vflux1,tflux1,qflux1
     *     ,uosurf,vosurf,uisurf,visurf,ogeoza
#ifdef TRACERS_ON
     *     ,trsrfflx,trsource
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,trs_glob
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
      USE TRDIAG_COM, only : taijn=>taijn_loc , tij_surf
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_isrc
      USE TRDIAG_COM, only : tajls=>tajls_loc,jls_source,jls_isrc
#ifdef TRACERS_WATER
     *     ,tij_evap,tij_grnd
#endif
#ifdef TRACERS_DRYDEP
     *     ,tij_drydep,tij_gsdep,itcon_dd,dtr_dd
#endif
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only : EMIS_SOURCE
#endif
      USE DVEG_COUPLER, only : step_dveg
      USE SOIL_DRV, only: earth

      IMPLICIT NONE
      integer rc

      INTEGER I,J,K,KR,JR,NS,NSTEPS,MODDSF,MODDD,ITYPE,IH,IHM,IDTYPE,IM1
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1K
     *     ,BETA,ELHX,ACE2,CDTERM,CDENOM,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,BETAUP,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG !,HSDEN,HSCON !,dEVdQS
     *     ,dT2,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER !,HSMUL,dHS,dQS,dTS
     *     ,PSK,Q1,THV1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,SHEAT,TRHEAT
     *     ,T2DEN,T2CON,T2MUL,FQEVAP ! ,QSDEN,QSCON,QSMUL,TGDEN
     *     ,Z1BY6L,EVAPLIM,F2,FSRI(2),HTLIM

      REAL*8 MA1, MSI1
      REAL*8, DIMENSION(NSTYPE,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *                                                       TGRND,TGRN2
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: ! S1BYG1 = BYRT3, Z1IBYL=Z1I/ALAMI,
     &     Z2LI3L=Z2LI/(3.*ALAMI), Z1LIBYL=Z1E/ALAMI
      REAL*8 QSAT,DQSATDT
c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 hemi,qg_sat,dtsurf,uocean,vocean,qsrf,us,vs,ws
      logical pole
c
#ifdef TRACERS_ON
      real*8 rhosrf0, totflux(ntm)
      integer n,nx,nsrc
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: tevaplim,trgrnd
      real*8  TEV,dTEVdTQS,tevap,dTQS,TDP,TDT1,frac
#ifdef TRACERS_SPECIAL_O18
     *     ,FRACVL,FRACVS,FRACLK
#endif
#endif
#ifdef TRACERS_DRYDEP
      real*8 tdryd, tdd, td1, rtsdt, rts
#endif
#ifdef TRACERS_DUST
      INTEGER :: n1
#endif
#endif

      INTEGER, PARAMETER :: n_idx1 = 11
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER,PARAMETER :: n_idx2=27
#else
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idx2=111
#else
      INTEGER, PARAMETER :: n_idx2 = 22
#endif
#endif
      INTEGER, PARAMETER :: n_idx3 = 2
      INTEGER, PARAMETER :: n_idx4 = n_idx1+n_idx2

      REAL*8, DIMENSION(n_idx4,grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUPT) :: diurn_part
      REAL*8, 
     &     DIMENSION(n_idx3,grid%J_STRT_HALO:grid%J_STOP_HALO,NDIUPT)::
     &     diurn_partb
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,n_idx3)::
     &     diurn_temp
      INTEGER :: idx1(n_idx1), idx2(n_idx2), idx3(n_idx3)
      INTEGER :: idx4(n_idx1+n_idx2)
      REAL*8 :: tmp(NDIUVAR)
      INTEGER :: ii, ivar
      REAL*8, DIMENSION(n_idx4, NDIUPT) :: DIURNSUM
      REAL*8, DIMENSION(n_idx3, NDIUPT) :: DIURNSUMb
      INTEGER, PARAMETER :: n_areg = 7
      REAL*8, DIMENSION(NREG,GRID%J_STRT_HALO:GRID%J_STOP_HALO,n_areg)::
     *     AREG_part
      REAL*8 :: AREGSUM(NREG,n_areg)
      INTEGER :: idx_areg(n_areg)

      INTEGER :: J_0, J_1, J_0H, J_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H, 
     *               J_STRT=J_0,        J_STOP=J_1)

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      IH=JHOUR+1
      IHM = IH+(JDATE-1)*24

C**** ZERO OUT ENERGY AND EVAPORATION FOR GROUND AND INITIALIZE TGRND
      DO J=J_0,J_1
      DO I=1,IM
        TGRND(2,I,J)=GTEMP(1,2,I,J)
        TGRND(3,I,J)=GTEMP(1,3,I,J)
C       TGRND(4,I,J)=GTEMP(1,4,I,J)
        TGRN2(2,I,J)=GTEMP(2,2,I,J)
        TGRN2(3,I,J)=GTEMP(2,3,I,J)
      END DO
      END DO

C**** Zero out fluxes summed over type and surface time step
      E0=0. ; E1=0. ; EVAPOR=0. ; RUNOE=0. ; ERUNOE=0.
      DMUA=0. ; DMVA=0. ; SOLAR=0.
#ifdef TRACERS_WATER
      TREVAPOR = 0. ; TRUNOE = 0.
#endif
#ifdef TRACERS_DRYDEP
      TRDRYDEP = 0.
      dtr_dd=0.
#endif
C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NIsurf TIMES EVERY HOUR
C****
      DO NS=1,NIsurf
         AREG_part = 0.
         MODDSF=MOD(NSTEPS+NS-1,NDASF*NIsurf+1)
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1
         MODDD=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??
C**** ZERO OUT FLUXES ACCUMULATED OVER SURFACE TYPES
         DTH1=0. ;  DQ1 =0. ;  uflux1=0. ; vflux1=0.
#ifdef TRACERS_ON
         trsrfflx = 0.
#endif

      call loadbl

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

      Call HALO_UPDATE(GRID, uosurf, FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID, vosurf, FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID, uisurf, FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID, visurf, FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID,u,FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID,v,FROM=SOUTH+NORTH)

      diurn_part=0

      idx1 = (/ IDD_SPR, 
     &     IDD_PT5, IDD_PT4, IDD_PT3, IDD_PT2, IDD_PT1,
     &     IDD_Q5,  IDD_Q4,  IDD_Q3,  IDD_Q2,  IDD_Q1 /)
      idx2 = (/ IDD_TS,  IDD_TG1, IDD_QS,  IDD_QG,  IDD_SWG,
     &          IDD_LWG, IDD_SH,  IDD_LH,  IDD_HZ0, IDD_UG,
     &          IDD_VG,  IDD_WG,  IDD_US,  IDD_VS,  IDD_WS,
     &          IDD_CIA, IDD_CM,  IDD_CH,  IDD_CQ,  IDD_EDS,
     &          IDD_DBL, IDD_EV
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     *     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf
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
     &     ,idd_grav,idd_turb
#endif
     &     /)
      idx3 = (/ IDD_DCF, IDD_LDC /)
      idx4(:n_idx1)   = idx1
      idx4(n_idx1+1:) = idx2


C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
!$OMP   PARALLEL DO PRIVATE (ACE2, BETA,BETAUP,   CM,CH,CQ,
!$OMP*  CDTERM,CDENOM,DSHDTG,DQGDTG,DEVDTG,DTRDTG,
!$OMP*  DF0DTG,DFDTG,DTG,DQ1X,DF1DTG,DSNDTG,   ! DEVDQS,
!$OMP*  DT2, EVAP,EVAPLIM,ELHX,EVHDT,EVHEAT,EVHDT0, ! DHS,DQS,DTS,
!$OMP*  F0DT,F1DT,F0,F1,F2,FSRI, HCG1,HCG2,  ! HSDEN,HSCON,
!$OMP*  HTLIM,I,ITYPE,IDTYPE,IM1, J,K, !,HSMUL
!$OMP*  KR, MA1,MSI1, PS,P1K,PLAND,PWATER,
!$OMP*  PLICE,PIJ,POICE,POCEAN,PTYPE,PSK, Q1, ! QSDEN,
!$OMP*  RHOSRF,RCDMWS,RCDHWS,RCDQWS, SHEAT,SRHEAT, ! QSCON,QSMUL,
!$OMP*  SNOW,SHDT, T2DEN,T2CON,T2MUL,TS,  ! TGDEN,
!$OMP*  THV1,TG,TG1,TG2,TRHDT,TRHEAT,Z1BY6L,
!$OMP*  HEMI,POLE,UOCEAN,VOCEAN,QG_SAT,US,VS,WS,QSRF,pbl_args,jr
#if defined(TRACERS_ON)
!$OMP*  ,n,nx,nsrc,rhosrf0,totflux
#if defined(TRACERS_WATER)
!$OMP*  ,tevaplim,tevap,trgrnd,TEV,dTEVdTQS,dTQS,TDP,TDT1,frac
#endif
#if defined(TRACERS_DRYDEP)
!$OMP*  ,tdryd,tdd,td1,rtsdt,rts
#endif
#ifdef TRACERS_DUST
!$OMP*  ,n1
#endif
#endif
!$OMP*  )
!$OMP*  SCHEDULE(DYNAMIC,2)
C
      DO J=J_0,J_1
      HEMI=1.
      IF(J.LE.JM/2) HEMI=-1.
      POLE= (J.EQ.1 .or. J.EQ.JM)

      IM1=IM
      DO I=1,IMAXJ(J)

      EVAPLIM = 0. ; HTLIM=0.  ! need initialisation
#ifdef TRACERS_WATER
      tevaplim = 0.
#endif
C****
C**** DETERMINE SURFACE CONDITIONS
C****
      PLAND=FLAND(I,J)
      PWATER=1.-PLAND
      PLICE=FLICE(I,J)
      POICE=RSI(I,J)*PWATER
      POCEAN=PWATER-POICE
      PIJ=P(I,J)
      PS=PEDN(1,I,J)
      PSK=PEK(1,I,J)
      P1K=PK(1,I,J)
      Q1=Q(I,J,1)
      THV1=T(I,J,1)*(1.+Q1*deltx)
CCC   JR=JREG(I,J)
      MA1=AM(1,I,J) !@var MA1 mass of lowest atmospheric layer (kg/m^2)
c     MSUM = (PS*100.)/GRAV !@var MSUM total mass of atmosphere (kg/m^2)
c     PGK = (PS*100.)**KAPA
c     PKDN = (GRAV*(MSUM-MA1*0.25))**KAPA
#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
C**** Calculate first layer tracer concentration
          trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
        end if
      end do
#endif

C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
         IF(MODDD.EQ.0) THEN
         DO KR=1,NDIUPT
           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             tmp(IDD_SPR)=+PS
             tmp(IDD_PT5)=+PSK*T(I,J,5)
             tmp(IDD_PT4)=+PSK*T(I,J,4)
             tmp(IDD_PT3)=+PSK*T(I,J,3)
             tmp(IDD_PT2)=+PSK*T(I,J,2)
             tmp(IDD_PT1)=+PSK*T(I,J,1)
             tmp(IDD_Q5)=+Q(I,J,5)
             tmp(IDD_Q4)=+Q(I,J,4)
             tmp(IDD_Q3)=+Q(I,J,3)
             tmp(IDD_Q2)=+Q(I,J,2)
             tmp(IDD_Q1)=+Q1
             DIURN_part(1:n_idx1,J,kr)=DIURN_part(1:n_idx1,J,kr)+
     &            tmp(idx1(:))
           END IF
         END DO
         END IF
C**** save some ocean diags regardless of PTYPE
C**** SSH does not work for qflux/fixed SST configurations
         IF (FOCEAN(I,J).gt.0. .and. MODDSF.eq.0) THEN
           AIJ(I,J,IJ_TGO)=AIJ(I,J,IJ_TGO)+GTEMP(1,1,I,J)
           AIJ(I,J,IJ_SSS)=AIJ(I,J,IJ_SSS)+SSS(I,J)
           AIJ(I,J,IJ_SSH)=AIJ(I,J,IJ_SSH)+OGEOZA(I,J)*BYGRAV+
     *          RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS 
         END IF
C****
      DO ITYPE=1,3       ! no earth type
  !    ipbl(i,j,itype)=0
!!!      SELECT CASE (ITYPE)
C****
C**** OPEN OCEAN/LAKE
C****
!!!      CASE (1)
      if ( ITYPE == 1 ) then

      PTYPE=POCEAN
      IF (PTYPE.gt.0) THEN
      TG1=GTEMP(1,1,I,J)
      TG2=GTEMP(2,1,I,J)   ! diagnostic only
      IF (FLAKE(I,J).gt.0) THEN
C**** limit evap/cooling if between MINMLD and 40cm, no evap below 40cm
        IF (MWL(I,J).lt.MINMLD*RHOW*FLAKE(I,J)*DXYP(J)) THEN
          EVAPLIM=MAX(0.5*(MWL(I,J)/(FLAKE(I,J)*DXYP(J))-0.4d0*RHOW),
     *         0d0)
        ELSE
          EVAPLIM=MWL(I,J)/(FLAKE(I,J)*DXYP(J))-(0.5*MINMLD+0.2d0)*RHOW
        END IF
#ifdef TRACERS_WATER
C**** limit on tracer evporation from lake
        TEVAPLIM(1:NTX)=EVAPLIM*GTRACER(NTIX(1:NTX),1,I,J)
#endif
        HTLIM = GML(I,J)/(FLAKE(I,J)*DXYP(J)) + 0.5*LHM*EVAPLIM
        IDTYPE=ITLAKE
        uocean = 0. ; vocean = 0. ! no velocities for lakes
      ELSE
        IDTYPE=ITOCEAN
        if (UOdrag.eq.1) then ! use uocean for drag calculation
C**** Convert UOSURF,VOSURF from C grid to A grid
C**** Note that uosurf,vosurf start with j=1, (not j=2 as in atm winds)
          if (pole) then
            uocean = 0. ; vocean = 0.
            do k=1,kmaxj(j)
              uocean = uocean + rapj(k,j)*(uosurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)-hemi*vosurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
              vocean = vocean + rapj(k,j)*(vosurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)+hemi*uosurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
            end do
          else
            uocean = 0.5*(uosurf(i,j)+uosurf(im1,j))
            vocean = 0.5*(vosurf(i,j)+vosurf(i,j-1))
          end if
        else
          uocean=0. ; vocean=0.
        end if
      END IF
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(1,I,J)=SOLAR(1,I,J)+DTSURF*SRHEAT
            OA(I,J,5)=OA(I,J,5)+SRHEAT*DTSURF
      BETA=1.
      ELHX=LHE
      END IF
C****
C**** OCEAN/LAKE ICE
C****
!!!      CASE (2)
      else if ( ITYPE == 2 ) then

      PTYPE=POICE
      IF (PTYPE.gt.0) THEN
      IF (FLAKE(I,J).gt.0) THEN
        IDTYPE=ITLKICE
        uocean = 0. ; vocean = 0. ! no dynamic ice for lakes
      ELSE
        IDTYPE=ITOICE
        if (UOdrag.eq.1) then ! use ice velcoities in drag calculation
C**** Convert UISURF,VISURF from C grid to A grid
C**** Note that uisurf,visurf start with j=1, (not j=2 as in atm winds)
          if (pole) then
            uocean = 0. ; vocean = 0.
            do k=1,kmaxj(j)
              uocean = uocean + rapj(k,j)*(uisurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)-hemi*visurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
              vocean = vocean + rapj(k,j)*(visurf(idij(k,i,j),idjj(k,j)
     *             -1)*cosip(k)+hemi*uisurf(idij(k,i,j),idjj(k,j)-1)
     *             *sinip(k))
            end do
          else
            uocean = 0.5*(uisurf(i,j)+uisurf(im1,j))
            vocean = 0.5*(visurf(i,j)+visurf(i,j-1))
          end if
        else
          uocean = 0. ; vocean =0.
        end if
      END IF
      SNOW=SNOWI(I,J)
      TG1=TGRND(2,I,J)
      TG2=TGRN2(2,I,J)
      MSI1 = SNOW+ACE1I ! snow and first layer ice mass (kg/m^2)
      ACE2=MSI(I,J) ! second (physical) layer ice mass (kg/m^2)
      dF1dTG = 2./(ACE1I*BYRLI+SNOW*BYRLS)
      HCG1 = SHI*XSI(1)*MSI1 ! heat capacity of top ice layer (J/C*m^2)
      HCG2 = SHI*XSI(2)*MSI1 ! heat capacity of second layer ice
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      SOLAR(2,I,J)=SOLAR(2,I,J)+DTSURF*SRHEAT
C**** fraction of solar radiation leaving layer 1 and 2
      IF (SRHEAT.gt.0) THEN ! only bother if there is sun
        call solar_ice_frac(SNOW,ACE2,FLAG_DSWS(I,J),FSRI,2)
      ELSE
        FSRI(1:2) = 0
      END IF
            OA(I,J,12)=OA(I,J,12)+SRHEAT*DTSURF
      BETA=1.
      ELHX=LHS

      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)

      END IF
C****
C**** LAND ICE
C****
!!!      CASE (3)
      else if ( ITYPE == 3 ) then

      PTYPE=PLICE
      IF (PTYPE.gt.0) THEN
      IDTYPE=ITLANDI
      SNOW=SNOWLI(I,J)
      TG1=TGRND(3,I,J)
      TG2=GTEMP(2,3,I,J)
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
      CDTERM=TG2
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)
      HCG1=HC1LI+SNOW*SHI
      BETA=1.
      ELHX=LHS
      uocean = 0. ; vocean = 0. ! no land ice velocity
      END IF
!!!      END SELECT
      endif
C****
      IF (PTYPE.gt.0) THEN
C****
C**** BOUNDARY LAYER INTERACTION
C****
  !    TKV=THV1*PSK  ! TKV is referenced to the surface pressure
  !    pbl_args%ZS1=.5*DSIG(1)*RGAS*BYGRAV*TKV*PIJ/PMID(1,I,J)
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.

      TG=TG1+TF
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG_SAT=0.98d0*QG_SAT
      pbl_args%TGV=TG*(1.+QG_SAT*deltx)
   !   psurf=PS   ! extra values to pass to PBL, possibly temporary
   !   trhr0 = TRHR(0,I,J)
#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
C**** set defaults
        trsfac(nx)=0.
        totflux(nx)=0.
        trconstflx(nx)=0.
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
#ifdef TRACERS_WATER
        tr_evap_max(nx)=1.
C**** The select is used to distinguish water from gases or particle
!!!        select case (tr_wd_TYPE(n))
!!!        case (nWATER)
        if ( tr_wd_TYPE(n) == nWATER ) then
          trgrnd(nx)=gtracer(n,itype,i,j)*QG_SAT
C**** trsfac and trconstflx are multiplied by cq*wsh in PBL
          trsfac(nx)=1.
          trconstflx(nx)=trgrnd(nx)
!!!        case (nGAS, nPART)
        else if ( tr_wd_TYPE(n) == nGAS .or.
     &         tr_wd_TYPE(n) == nPART ) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
          rhosrf0=100.*ps/(rgas*pbl_args%tgv) ! estimated surface density
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) then
            trsfac(nx)=1. 
     &      !then multiplied by deposition velocity in PBL
#ifdef TRACERS_WATER
            tr_evap_max(nx)=1.d30
            trgrnd(nx)=0.
#endif
          end if
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
C**** Now send kg/m^2/s to PBL, and dived by rho there.
          do nsrc=1,ntsurfsrc(n)
            totflux(nx) = totflux(nx)+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux(nx)*bydxyp(j)   ! kg/m^2/s
#ifdef TRACERS_WATER
!!!        end select
        endif
#endif
      end do
#endif
C =====================================================================
      pbl_args%dtsurf = dtsurf
      pbl_args%TKV=THV1*PSK     ! TKV is referenced to the surface pressure
      pbl_args%ZS1=.5*DSIG(1)*RGAS*BYGRAV* pbl_args%TKV*PIJ/PMID(1,I,J)
      pbl_args%qg_sat = qg_sat
      pbl_args%qg_aver = qg_sat   ! QG_AVER=QG_SAT
      pbl_args%hemi = hemi
      pbl_args%pole = pole
      pbl_args%evap_max = 1.
      pbl_args%fr_sat = 1. ! entire surface is saturated
      pbl_args%uocean = uocean
      pbl_args%vocean = vocean
      pbl_args%psurf = PS
      pbl_args%trhr0 = TRHR(0,I,J)
      
      CALL PBL(I,J,ITYPE,PTYPE,pbl_args)

      us = pbl_args%us
      vs = pbl_args%vs
      ws = pbl_args%ws
      qsrf = pbl_args%qsrf
      CM = pbl_args%cm
      CH = pbl_args%ch
      CQ = pbl_args%cq
C =====================================================================
      TS=pbl_args%TSV/(1.+QSRF*deltx)
C**** CALCULATE RHOSRF*CM*WS AND RHOSRF*CH*WS
      RHOSRF=100.*PS/(RGAS*pbl_args%TSV)
      RCDMWS=CM*pbl_args%WSM*RHOSRF
      RCDHWS=CH*pbl_args%WSH*RHOSRF
      RCDQWS=CQ*pbl_args%WSH*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2)
      SHEAT=SHA*RCDHWS*(TS-TG)
      BETAUP = BETA
      IF (QSRF .GT. QG_SAT) BETAUP = 1.
      EVHEAT=(LHE+TG1*SHV)*BETAUP*RCDQWS*(QSRF-QG_SAT)
      TRHEAT=TRHR(0,I,J)-STBO*(TG*TG)*(TG*TG)
C****
!!!      SELECT CASE (ITYPE)

!!!      CASE (1) ! FLUXES USING IMPLICIT TIME STEP FOR OCEAN POINTS
      if ( ITYPE == 1) then
c       DSHDTG=-RCDHWS*SHA
c       dEVdQS = LHE*RCDQWS
c       dHS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
c    A       ((1.+2.*S1BYG1)*DTSURF*PGK*RCDHWS+MA1*PKDN)
c       dTS = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/
c    A       (MA1*PKDN*SHA-(1.+2.*S1BYG1)*DTSURF*PGK*DSHDTG)
c       dQS = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/
c    A       ((1.+2.*S1BYG1)*DTSURF*dEVdQS+MA1*LHE)
c       SHDT = DTSURF*(SHEAT-dTS*DSHDTG)
c       EVHDT=DTSURF*(EVHEAT+dQS*dEVdQS) ! latent heat flux
        SHDT = DTSURF*SHEAT
        EVHDT=DTSURF*EVHEAT              ! latent heat flux
        TRHDT=DTSURF*TRHEAT
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        SHDTT(I,J)=SHDT
#endif
C****
!!!      CASE (2) ! FLUXES USING IMPLICIT TIME STEP FOR ICE POINTS
      else if ( ITYPE == 2 ) then

! heat flux on first/second/third layers (W/m^2)
        F1 = (TG1-TG2)*dF1dTG + SRHEAT*FSRI(1)
        F2 = SRHEAT*FSRI(2)
        EVHEAT = LHE*RCDQWS*(QSRF-QG_SAT) ! latent heat flux (W/m^2)
        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        dSNdTG=-RCDHWS*SHA
        dQGdTG=QG_SAT*DQSATDT(TG,ELHX) ! d(QG)/dTG
        dEVdTG = -dQGdTG*LHE*RCDQWS ! d(EVHEAT)/dTG
        dTRdTG = -4*STBO*TG*TG*TG ! d(TRHEAT)/dTG
        dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG
cc      dSNdHS = RCDHWS ! d(SHEAT)/dHS - kg/(sec*m^2)
c       dEVdQS = LHE*RCDQWS     ! d(EVHEAT)/dQS
c       HSDEN = -(1.+2.*S1BYG1)*DTSURF*PGK*BETA*dSNdTG+MA1*PKDN*SHA
c       HSCON = -(1.+2.*S1BYG1)*DTSURF*PGK*SHEAT/HSDEN ! (J*sec)/kg
c       HSMUL=-(1.+2.*S1BYG1)*DTSURF*PGK*BETA*dSNdTG/HSDEN ! J/(kg*degC)
c       QSDEN = (1.+2.*S1BYG1)*BETA*DTSURF*dEVdQS+MA1*LHE
c       QSCON = -(1.+2.*S1BYG1)*DTSURF*EVHEAT/QSDEN
c       QSMUL = -(1.+2.*S1BYG1)*DTSURF*BETA*dEVdTG/QSDEN
        T2DEN = HCG2+BETA*DTSURF*dF1dTG
        T2CON = DTSURF*(F1-F2)/T2DEN
        T2MUL = BETA*DTSURF*dF1dTG/T2DEN
c       TGDEN = HCG1-BETA*DTSURF*(dF0dTG-dF1dTG-
c    A       HSMUL*dSNdTG+QSMUL*dEVdQS+T2MUL*dF1dTG) ! W/(m^2*degC)
c       dTG = DTSURF*(F0-F1+BETA*
c    A       (QSCON*dEVdQS-HSCON*dSNdTG+T2CON*dF1dTG))/TGDEN ! degC

        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)

        IF (TG1+dTG .GT. 0.) dTG = -TG1
c       dHS = HSCON+HSMUL*dTG   ! (J*sec)/kg
c       dQS = QSCON+QSMUL*dTG
        dT2 = T2CON+T2MUL*dTG
c       SHDT = DTSURF*(SHEAT+BETA*((dTG-dHS)*dSNdTG)) ! sensible
c       EVHDT = DTSURF*(EVHEAT+BETA*(dTG*dEVdTG+dQS*dEVdQS)) ! latent
        SHDT = DTSURF*(SHEAT+BETA*dTG*dSNdTG) ! sensible
        EVHDT = DTSURF*(EVHEAT+BETA*dTG*dEVdTG) ! latent

        TRHDT = DTSURF*(TRHEAT+BETA*dTG*dTRdTG) ! thermal flux (J/m^2)
        F1DT = DTSURF*(F1+BETA*(dTG*dF1dTG-dT2*dF1dTG))
        TG1 = TG1+dTG          ! first layer sea ice temperature (degC)
        TG2 = TG2+dT2          ! second layer sea ice temperature (degC)
        TGRN2(ITYPE,I,J) = TG2
C****
!!!      CASE (3) ! IMPLICIT TIME STEP OVER LANDICE
      else if ( ITYPE == 3 ) then

        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
        DSHDTG=-RCDHWS*SHA
        DQGDTG=QG_SAT*DQSATDT(TG,ELHX)
        DEVDTG=-RCDQWS*LHE*BETAUP*DQGDTG
        DTRDTG=-4.*STBO*TG*TG*TG
        DF0DTG=DSHDTG+DEVDTG+DTRDTG
        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)
        SHDT=DTSURF*(SHEAT+DTG*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+DTG*DEVDTG)
        TRHDT=DTSURF*(TRHEAT+DTG*DTRDTG)
        F1DT=DTSURF*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
        TG1=TG1+DTG

!!!      END SELECT
      endif

C**** CALCULATE EVAPORATION
      DQ1X =EVHDT/((LHE+TG1*SHV)*MA1)
      EVHDT0=EVHDT
C**** Limit evaporation if lake mass is at minimum
      IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and.
     *     (EVAPOR(I,J,1)-DQ1X*MA1).gt.EVAPLIM) THEN
        if (QCHECK) WRITE(99,*) "Lake EVAP limited: I,J,EVAP,MWL",I,J
     *     ,EVAPOR(I,J,1)-DQ1X*MA1,MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J))
        DQ1X=(EVAPOR(I,J,1)-EVAPLIM)*BYAM(1,I,J)
      ELSEIF (DQ1X.GT.Q1+DQ1(I,J)) THEN
        DQ1X=(Q1+DQ1(I,J))
      ELSE
        GO TO 3720
      END IF
      EVHDT=DQ1X*(LHE+TG1*SHV)*MA1
      IF (ITYPE.NE.1) TG1=TG1+(EVHDT-EVHDT0)/HCG1
 3720 EVAP=-DQ1X*MA1

#ifdef TRACERS_ON
C**** Loop over tracers
      DO NX=1,NTX
        N=NTIX(NX)
#ifdef TRACERS_WATER
        if (tr_wd_TYPE(n).eq.nWATER) THEN
C****
C**** Calculate Water Tracer Evaporation
C****
          IF (ITYPE.EQ.1) THEN  ! OCEAN
C**** do calculation implicitly for TQS
#ifdef TRACERS_SPECIAL_O18
            TEV=-RCDQWS*(trs(nx)-trgrnd(nx)*
     *           fracvl(tg1,trname(n)))*FRACLK(pbl_args%WSM,trname(n))
            dTEVdTQS =-RCDQWS*FRACLK(pbl_args%WSM,trname(n))
#else
            TEV=-RCDQWS*(trs(nx)-trgrnd(nx))
            dTEVdTQS =-RCDQWS
#endif
c           dTQS = -(1.+2.*S1BYG1)*DTSURF*TEV/
c    *           ((1.+2.*S1BYG1)*DTSURF*dTEVdTQS-MA1)
c           TEVAP=DTSURF*(TEV+dTQS*dTEVdTQS)
            TEVAP=DTSURF*TEV
          ELSE                  ! ICE AND LAND ICE
C**** tracer flux is set by source tracer concentration
            IF (EVAP.GE.0) THEN ! EVAPORATION
              TEVAP=EVAP*trgrnd(nx)/QG_SAT
            ELSE                ! DEW (fractionates)
#ifdef TRACERS_SPECIAL_O18
              IF (TG1.gt.0) THEN
                frac=FRACVL(TG1,trname(n))
              ELSE
                frac=FRACVS(TG1,trname(n))
              END IF
              TEVAP=EVAP*trs(nx)/(QSRF*frac+teeny)
#else
              TEVAP=EVAP*trs(nx)/(QSRF+teeny)
#endif
            END IF
          END IF
C**** Limit evaporation if lake mass is at minimum
          IF (ITYPE.EQ.1 .and. FLAKE(I,J).GT.0 .and.
     *         (TREVAPOR(n,1,I,J)+TEVAP.gt.TEVAPLIM(NX))) THEN
            IF (QCHECK) WRITE(99,*) "Lake TEVAP limited: I,J,TEVAP,TMWL"
     *           ,N,TREVAPOR(n,1,I,J)+TEVAP,TEVAPLIM(NX)
            TEVAP= TEVAPLIM(NX)-TREVAPOR(n,1,I,J)
          END IF
          TDP = TEVAP*DXYP(J)*ptype
          TDT1 = trsrfflx(I,J,n)*DTSURF
          IF (TRM(I,J,1,n)+TDT1+TDP.lt.0..and.tdp.lt.0) THEN
            IF (QCHECK) WRITE(99,*) "LIMITING TRDEW",I,J,N,TDP,TRM(I,J,1
     *           ,n),TDT1
            TEVAP = -(TRM(I,J,1,n)+TDT1)/(DXYP(J)*ptype)
            trsrfflx(I,J,n)= - TRM(I,J,1,n)/DTSURF
          ELSE
            trsrfflx(I,J,n)=trsrfflx(I,J,n)+TDP/DTSURF
          END IF
          TREVAPOR(n,ITYPE,I,J)=TREVAPOR(n,ITYPE,I,J)+TEVAP
        END IF
#endif
c#ifdef TRACERS_AEROSOLS_Koch
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
#endif
        end select
c#endif
#ifdef TRACERS_DRYDEP
C****
C**** Calculate Tracer Dry Deposition (including gravitational settling)
C****
        if(dodrydep(n))then
          rts=rhosrf*trs(nx)
          rtsdt=rts*dtsurf                             ! kg*s/m^3
          tdryd=-rtsdt*(dep_vel(n)+gs_vel(n))          ! kg/m2
          tdd = tdryd*dxyp(j)*ptype                    ! kg
          td1 = (trsrfflx(i,j,n)+totflux(nx))*dtsurf   ! kg
          if (trm(i,j,1,n)+td1+tdd.le.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd surfce",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,trs(nx),trtop(nx)
            tdd= -max(trm(i,j,1,n)+td1,0d0)
            tdryd=tdd/(dxyp(j)*ptype)
            trsrfflx(i,j,n)= - trm(i,j,1,n)/dtsurf
          else
            trsrfflx(i,j,n)=trsrfflx(i,j,n)+tdd/dtsurf
          end if
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
      END DO
#endif
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
C**** Limit heat fluxes out of lakes if near minimum depth
      IF (ITYPE.eq.1 .and. FLAKE(I,J).gt.0 .and.
     *     E0(I,J,1)+F0DT+HTLIM.lt.0) THEN
        if (QCHECK) write(6,*) "Limiting heat flux from lake",i,j,SHDT
     *       ,F0DT,E0(I,J,1),DTSURF*SRHEAT,TRHDT,EVHDT,HTLIM
        SHDT = -(HTLIM+E0(I,J,1)+DTSURF*SRHEAT+TRHDT+EVHDT)
        F0DT = -E0(I,J,1)-HTLIM
        if (QCHECK) write(6,*) "New SHDT,F0DT",i,j,SHDT,F0DT
      END IF
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP
      TGRND(ITYPE,I,J)=TG1
      DTH1(I,J)=DTH1(I,J)-SHDT*PTYPE/(SHA*MA1*P1K)
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
      DMUA(I,J,ITYPE)=DMUA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*(US-UOCEAN)
      DMVA(I,J,ITYPE)=DMVA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*(VS-VOCEAN)
      uflux1(i,j)=uflux1(i,j)+PTYPE*RCDMWS*(US-UOCEAN)
      vflux1(i,j)=vflux1(i,j)+PTYPE*RCDMWS*(VS-UOCEAN)
C****
C**** ACCUMULATE DIAGNOSTICS FOR EACH SURFACE TIME STEP AND ITYPE
C****
        AJ(J,J_EVAP ,IDTYPE)=AJ(J,J_EVAP ,IDTYPE)+ EVAP*PTYPE
        AJ(J,J_EVHDT,IDTYPE)=AJ(J,J_EVHDT,IDTYPE)+EVHDT*PTYPE
        AJ(J,J_SHDT ,IDTYPE)=AJ(J,J_SHDT ,IDTYPE)+ SHDT*PTYPE
        AJ(J,J_TRHDT,IDTYPE)=AJ(J,J_TRHDT,IDTYPE)+TRHDT*PTYPE
        IF(MODDSF.EQ.0) THEN
          AJ(J,J_TSRF,IDTYPE)=AJ(J,J_TSRF,IDTYPE)+(TS-TF)*PTYPE
          AJ(J,J_TYPE,IDTYPE)=AJ(J,J_TYPE,IDTYPE)+        PTYPE
          AJ(J,J_TG1 ,IDTYPE)=AJ(J,J_TG1 ,IDTYPE)+    TG1*PTYPE
          AJ(J,J_TG2 ,IDTYPE)=AJ(J,J_TG2 ,IDTYPE)+    TG2*PTYPE
        END IF
C**** QUANTITIES ACCUMULATED FOR REGIONS IN DIAGJ
CCC     AREG(JR,J_TRHDT)=AREG(JR,J_TRHDT)+TRHDT*PTYPE*DXYP(J)
CCC     AREG(JR,J_SHDT )=AREG(JR,J_SHDT )+SHDT *PTYPE*DXYP(J)
CCC     AREG(JR,J_EVHDT)=AREG(JR,J_EVHDT)+EVHDT*PTYPE*DXYP(J)
CCC     AREG(JR,J_EVAP )=AREG(JR,J_EVAP )+EVAP *PTYPE*DXYP(J)
CCC     IF(MODDSF.EQ.0)
CCC  *       AREG(JR,J_TSRF)=AREG(JR,J_TSRF)+(TS-TF)*PTYPE*DXYP(J)
        JR=JREG(I,J)
        AREG_part(JR,J,1)=AREG_part(JR,J,1)+TRHDT*PTYPE*DXYP(J)
        AREG_part(JR,J,2)=AREG_part(JR,J,2)+SHDT *PTYPE*DXYP(J)
        AREG_part(JR,J,3)=AREG_part(JR,J,3)+EVHDT*PTYPE*DXYP(J)
        AREG_part(JR,J,4)=AREG_part(JR,J,4)+EVAP *PTYPE*DXYP(J)
        IF(MODDSF.EQ.0) THEN
          AREG_part(JR,J,5)=AREG_part(JR,J,5)+(TS-TF)*PTYPE*DXYP(J)
          AREG_part(JR,J,6)=AREG_part(JR,J,6)+    TG1*PTYPE*DXYP(J)
          AREG_part(JR,J,7)=AREG_part(JR,J,7)+    TG2*PTYPE*DXYP(J)
        END IF
C
C**** QUANTITIES ACCUMULATED FOR LATITUDE-LONGITUDE MAPS IN DIAGIJ
        AIJ(I,J,IJ_SHDT)=AIJ(I,J,IJ_SHDT)+SHDT*PTYPE
        IF(MODDSF.EQ.0) THEN
          AIJ(I,J,IJ_TRSDN)=AIJ(I,J,IJ_TRSDN)+TRHR(0,I,J)*PTYPE
          AIJ(I,J,IJ_TRSUP)=AIJ(I,J,IJ_TRSUP)+(TRHR(0,I,J)-TRHDT/DTSURF)
     *         *PTYPE
        END IF
        IF(MODRD.EQ.0) AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)+TRHDT
     *       *PTYPE/DTSRC
        AIJ(I,J,IJ_SRTR)=AIJ(I,J,IJ_SRTR)+(SRHEAT*DTSURF+TRHDT)*PTYPE
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+(SRHEAT*DTSURF+TRHDT+SHDT
     *       +EVHDT)*PTYPE
        AIJ(I,J,IJ_EVAP)=AIJ(I,J,IJ_EVAP)+EVAP*PTYPE
        IF(MODDSF.EQ.0) THEN
          AIJ(I,J,IJ_WS)=AIJ(I,J,IJ_WS)+WS*PTYPE
          AIJ(I,J,IJ_TS)=AIJ(I,J,IJ_TS)+(TS-TF)*PTYPE
          AIJ(I,J,IJ_US)=AIJ(I,J,IJ_US)+US*PTYPE
          AIJ(I,J,IJ_VS)=AIJ(I,J,IJ_VS)+VS*PTYPE
          AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RCDMWS*pbl_args%WSM*PTYPE
          AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RCDMWS*(US-UOCEAN)*PTYPE
          AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RCDMWS*(VS-VOCEAN)*PTYPE
          AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QSRF*PTYPE
          AIJ(I,J,IJ_TG1)=AIJ(I,J,IJ_TG1)+TG1*PTYPE
          AIJ(I,J,IJ_PBLHT)=AIJ(I,J,IJ_PBLHT)+pbl_args%dbl*PTYPE

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          aij(i,j,ij_wsgcm)=aij(i,j,ij_wsgcm)+pbl_args%wsgcm*ptype
          aij(i,j,ij_wspdf)=aij(i,j,ij_wspdf)+pbl_args%wspdf*ptype
          aij(i,j,ij_wdry)=aij(i,j,ij_wdry)+wsubwd*ptype
          aij(i,j,ij_wtke)=aij(i,j,ij_wtke)+wsubtke*ptype
          aij(i,j,ij_wmoist)=aij(i,j,ij_wmoist)+wsubwm*ptype
#endif

        END IF

c     ..........
c     save global variables for subdaily diagnostics
c     ..........

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

C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
        IF(MODDD.EQ.0) THEN
          DO KR=1,NDIUPT
            IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
              tmp(IDD_TS)=+TS*PTYPE
              tmp(IDD_TG1)=+(TG1+TF)*PTYPE
              tmp(IDD_QS)=+QSRF*PTYPE
              tmp(IDD_QG)=+QG_SAT*PTYPE
              tmp(IDD_SWG)=+SRHEAT*DTSURF
     *             *PTYPE
              tmp(IDD_LWG)=+TRHDT*PTYPE
              tmp(IDD_SH)=+SHDT*PTYPE
              tmp(IDD_LH)=+EVHDT*PTYPE
              tmp(IDD_HZ0)=
     *             +(SRHEAT*DTSURF+TRHDT+SHDT+EVHDT)*PTYPE
              tmp(IDD_UG)=+pbl_args%UG*PTYPE
              tmp(IDD_VG)=+pbl_args%VG*PTYPE
              tmp(IDD_WG)=+pbl_args%WG*PTYPE
              tmp(IDD_US)=+US*PTYPE
              tmp(IDD_VS)=+VS*PTYPE
              tmp(IDD_WS)=+WS*PTYPE
              tmp(IDD_CIA)=+pbl_args%PSI*PTYPE
              tmp(IDD_CM)=+CM*PTYPE
              tmp(IDD_CH)=+CH*PTYPE
              tmp(IDD_CQ)=+CQ*PTYPE
              tmp(IDD_EDS)=+pbl_args%KHS*PTYPE
              tmp(IDD_DBL)=+pbl_args%DBL*PTYPE
              tmp(IDD_EV)=+EVAP*PTYPE

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
            IF (adiurn_dust == 1) THEN
              tmp(idd_wsgcm)=+pbl_args%wsgcm*ptype
              tmp(idd_wspdf)=+pbl_args%wspdf*ptype
              tmp(idd_wtke)=+wsubtke*ptype
              tmp(idd_wd)=+wsubwd*ptype
              tmp(idd_wm)=+wsubwm*ptype
            END IF
#endif
#ifdef TRACERS_DUST
              IF (adiurn_dust == 1) THEN

                tmp(idd_turb)=0.D0
                tmp(idd_grav)=0.D0
#ifdef TRACERS_DRYDEP
                DO n=1,Ntm_dust
                  n1=n_clay+n-1
                  IF (dodrydep(n1)) THEN
                    tmp(idd_turb)=+ptype*rts*dep_vel(n1)
                    tmp(idd_grav)=+ptype*rts*gs_vel(n1)
                  END IF
                END DO
#endif

                tmp(idd_ws2)=+ws*ws*ptype
                tmp(idd_ustar)=+pbl_args%ustar*ptype
                tmp(idd_us3)=+ptype*pbl_args%ustar**3
                tmp(idd_stress)=+rcdmws*pbl_args%wsm*ptype
                tmp(idd_lmon)=+lmonin*ptype
                tmp(idd_rifl)=
     &               +ptype*grav*(ts-(tg1+tf))*pbl_args%zgs
     &               /(ws*ws*(tg1+tf))

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
     *               +ptype*sqrt( uabl(1,i,j,itype)*uabl(1,i,j,itype)
     *               +            vabl(1,i,j,itype)*vabl(1,i,j,itype)) 
                tmp(idd_uvabl2)=
     *               +ptype*sqrt( uabl(2,i,j,itype)*uabl(2,i,j,itype)
     *               +            vabl(2,i,j,itype)*vabl(2,i,j,itype))
                tmp(idd_uvabl3)=
     *               +ptype*sqrt( uabl(3,i,j,itype)*uabl(3,i,j,itype)
     *               +            vabl(3,i,j,itype)*vabl(3,i,j,itype))
                tmp(idd_uvabl4)=
     *               +ptype*sqrt( uabl(4,i,j,itype)*uabl(4,i,j,itype)
     *               +            vabl(4,i,j,itype)*vabl(4,i,j,itype))
                tmp(idd_uvabl5)=
     *               +ptype*sqrt( uabl(5,i,j,itype)*uabl(5,i,j,itype)
     *               +            vabl(5,i,j,itype)*vabl(5,i,j,itype))
                tmp(idd_uvabl6)=
     *               +ptype*sqrt( uabl(6,i,j,itype)*uabl(6,i,j,itype)
     *               +            vabl(6,i,j,itype)*vabl(6,i,j,itype))
                tmp(idd_uvabl7)=
     *               +ptype*sqrt( uabl(7,i,j,itype)*uabl(7,i,j,itype)
     *               +            vabl(7,i,j,itype)*vabl(7,i,j,itype))
                tmp(idd_uvabl8)=
     *               +ptype*sqrt( uabl(8,i,j,itype)*uabl(8,i,j,itype)
     *               +            vabl(8,i,j,itype)*vabl(8,i,j,itype))

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
              DIURN_part(n_idx1+1:n_idx4,J,kr)=
     &             DIURN_part(n_idx1+1:n_idx4,J,kr)+tmp(idx2(:))
            END IF
          END DO
        END IF
C****
#ifdef TRACERS_ON
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
#ifdef TRACERS_WATER
          if (tr_wd_type(n).eq.nWater) then
            taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *           trevapor(n,itype,i,j)*ptype
            tajls(j,1,jls_source(1,n))=tajls(j,1,jls_source(1,n))
     *           +trevapor(n,itype,i,j)*ptype
            if (focean(i,j).gt.0) tajls(j,1,jls_source(2,n))=tajls(j,1
     *           ,jls_source(2,n))+trevapor(n,itype,i,j)*ptype
          end if
          taijn(i,j,tij_grnd,n)=taijn(i,j,tij_grnd,n)+
     *         gtracer(n,itype,i,j)*ptype
#endif
        end if
      end do
#endif
C****
C**** SAVE SOME TYPE DEPENDENT FLUXES/DIAGNOSTICS
C****
!!!      SELECT CASE (ITYPE)
!!!      CASE (1)  ! ocean
      if ( ITYPE == 1 ) then
        OA(I,J,6)=OA(I,J,6)+TRHDT
        OA(I,J,7)=OA(I,J,7)+SHDT
        OA(I,J,8)=OA(I,J,8)+EVHDT
        AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+EVAP*PTYPE
        IF (FOCEAN(I,J).gt.0) THEN
          AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC) +F0DT*PTYPE
          AIJ(I,J,IJ_FWOC)=AIJ(I,J,IJ_FWOC) -EVAP*PTYPE
        END IF
C****
!!!      CASE (2)  ! seaice
      else if ( ITYPE == 2 ) then
        IF (TG1.GT.TDIURN(I,J,7)) TDIURN(I,J,7) = TG1
        OA(I,J,9)=OA(I,J,9)+TRHDT
        OA(I,J,10)=OA(I,J,10)+SHDT
        OA(I,J,11)=OA(I,J,11)+EVHDT
        AIJ(I,J,IJ_F0OI) =AIJ(I,J,IJ_F0OI) +F0DT*PTYPE
        AIJ(I,J,IJ_EVAPI)=AIJ(I,J,IJ_EVAPI)+EVAP*PTYPE
C****
!!!      CASE (3) ! land ice
      else if ( ITYPE == 3 ) then
        IF (TG1.GT.TDIURN(I,J,8)) TDIURN(I,J,8) = TG1
        IF (MODDSF.eq.0) AIJ(I,J,IJ_TSLI)=AIJ(I,J,IJ_TSLI)+(TS-TF)
        AIJ(I,J,IJ_SHDTLI)=AIJ(I,J,IJ_SHDTLI)+SHDT
        AIJ(I,J,IJ_EVHDT)=AIJ(I,J,IJ_EVHDT)+EVHDT
        AIJ(I,J,IJ_TRHDT)=AIJ(I,J,IJ_TRHDT)+TRHDT
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+F0DT
        AIJ(I,J,IJ_EVAPLI)=AIJ(I,J,IJ_EVAPLI)+EVAP
C****
!!!      END SELECT
      endif
C****
      END IF
      END DO   ! end of itype loop
      IM1=I
      END DO   ! end of I loop

      END DO   ! end of J loop
!$OMP  END PARALLEL DO

      idx_areg = (/ J_TRHDT, J_SHDT, J_EVHDT, J_EVAP, J_TSRF, 
     &     J_TG1, J_TG2 /)
      CALL GLOBALSUM(grid, AREG_PART, AREGSUM)
      AREG(:,idx_areg) = AREG(:,idx_areg) + AREGSUM

      CALL GLOBALSUM(grid, DIURN_part, DIURNSUM)

      IF (AM_I_ROOT()) THEN
         ADIURN(ih,idx4,:)=ADIURN(ih,idx4,:)   + DIURNSUM
#ifndef TRACERS_DUST
#ifndef TRACERS_MINERALS
#ifndef TRACERS_QUARZHEM
         HDIURN(ihm,idx4,:)=HDIURN(ihm,idx4,:) + DIURNSUM
#endif
#endif
#endif
      END IF
      

C****
C**** dycamic vegetation time step
C****
      call step_dveg(dtsurf)
C****
C**** EARTH
C****
      CALL EARTH(NS,MODDSF,MODDD)
C****
C**** UPDATE FIRST LAYER QUANTITIES
C****
!$OMP  PARALLEL DO PRIVATE (I,J,FTEVAP,FQEVAP,P1K)
!$OMP*          SCHEDULE(DYNAMIC,2)
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        FTEVAP=0
        IF (DTH1(I,J)*T(I,J,1).lt.0) FTEVAP=-DTH1(I,J)/T(I,J,1)
        FQEVAP=0
        IF (DQ1(I,J).lt.0.and.Q(I,J,1).gt.0) FQEVAP=-DQ1(I,J)/Q(I,J,1)
! Z-moments should be set from PBL
        TMOM(:,I,J,1) = TMOM(:,I,J,1)*(1.-FTEVAP)
        QMOM(:,I,J,1) = QMOM(:,I,J,1)*(1.-FQEVAP)
        IF ( Q(I,J,1)+DQ1(I,J) .LT. qmin ) THEN
          WRITE(99,*)
     &         ITime,'I,J:',I,J,' Q1:',Q(I,J,1)+DQ1(I,J),'->',qmin
          dq1(i,j)=qmin-q(i,j,1)
          QMOM(:,I,J,1)=0.
        ENDIF
c****   retrieve fluxes
        P1K=PK(1,I,J)
        tflux1(i,j)=dth1(i,j)*(-AM(1,I,J)*P1K)/(dtsurf)
        qflux1(i,j)=dq1(i,j)*(-AM(1,I,J))/(dtsurf)
C**** Diurnal cycle of temperature diagnostics
        tdiurn(i,j,5)=tdiurn(i,j,5)+(tsavg(i,j)-tf)
        if(tsavg(i,j).gt.tdiurn(i,j,6)) tdiurn(i,j,6)=tsavg(i,j)
        if(tsavg(i,j).lt.tdiurn(i,j,9)) tdiurn(i,j,9)=tsavg(i,j)
      END DO
      END DO
!$OMP  END PARALLEL DO
#ifdef TRACERS_ON
C****
C**** Apply tracer surface sources and sinks
C****
      call apply_tracer_2Dsource(dtsurf)
#endif
c****
c**** apply surface fluxes to the first layer of the atmosphere
c****  (replaced with dummy sub when ATURB is used)
c****
      call apply_fluxes_to_atm(dtsurf)

C**** Call dry convection or aturb depending on rundeck
      CALL ATM_DIFFUS(1,1,dtsurf)
C****
C**** ACCUMULATE SOME ADDITIONAL BOUNDARY LAYER DIAGNOSTICS
C****
      IF(MODDD.EQ.0) THEN
        DIURN_partb = 0

        DO KR=1,NDIUPT
C**** CHECK IF DRY CONV HAS HAPPENED FOR THIS DIAGNOSTIC
C**** For distributed implementation - ensure point is on local process.          

          I = IJDD(1,KR)
          J = IJDD(2,KR)
          IF ((J >= J_0) .AND. (J <= J_1)) THEN
            IF(DCLEV(I,J).GT.1.) THEN
              tmp(1)=+1.
              tmp(2)=+DCLEV(I,J)
              DIURN_partb(:,J,KR)=DIURN_partb(:,J,KR)+tmp(1:2)
            END IF
          END IF
       END DO

       CALL GLOBALSUM(grid,  DIURN_partb, DIURNSUMb)

       IF (AM_I_ROOT()) THEN
          ADIURN(ih,idx3,:)=ADIURN(ih,idx3,:)   + DIURNSUMb
#ifndef TRACERS_DUST
#ifndef TRACERS_MINERALS
#ifndef TRACERS_QUARZHEM
          HDIURN(ihm,idx3,:)=HDIURN(ihm,idx3,:) + DIURNSUMb
#endif
#endif
#endif
       END IF

      END IF
C****
      END DO   ! end of surface time step

#ifdef TRACERS_DRYDEP
C**** Save for tracer dry deposition conservation quantity:
      do n=1,ntm
        if(dodrydep(n)) then
          call diagtcb(dtr_dd(:,n,1),itcon_dd(n,1),n)  ! turb dep
          if (itcon_dd(n,2).gt.0)
     *         call diagtcb(dtr_dd(:,n,2),itcon_dd(n,2),n) ! grav sett
        end if
      end do
#endif

      RETURN
C****
      END SUBROUTINE SURFCE
