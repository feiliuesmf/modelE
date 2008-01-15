C****   
C**** SURFACE.f    SURFACE fluxes    2006/12/21
C****
#include "rundeck_opts.h"

      SUBROUTINE SURFCE
!@sum SURFCE calculates the surface fluxes which include
!@+   sensible heat, evaporation, thermal radiation, and momentum
!@+   drag.  It also calculates instantaneous surface temperature,
!@+   surface specific humidity, and surface wind components.
!@auth Nobody will claim responsibilty
      USE CONSTANT, only : rgas,lhm,lhe,lhs
     *     ,sha,tf,rhow,shv,shi,stbo,bygrav,by6
     *     ,deltx,teeny,rhows,grav
      USE MODEL_COM, only : im,jm,dtsrc,nisurf,u,v,t,p,q
     *     ,idacc,ndasf,fland,flice,focean,IVSP,IVNP
     *     ,nday,modrd,itime,jhour,itocean
     *     ,itoice,itlake,itlkice,itlandi,qcheck,UOdrag,jdate
      USE DOMAIN_DECOMP, only : GRID, GET, CHECKSUM, HALO_UPDATE, SOUTH
      USE DOMAIN_DECOMP, only : NORTH
      USE DOMAIN_DECOMP, only : AM_I_ROOT, GLOBALSUM
      USE GEOM, only : dxyp,imaxj,bydxyp,idjj,idij,rapj,kmaxj,sinip
     *     ,cosip
      USE SOMTQ_COM, only : tmom,qmom,mz
      USE DYNAMICS, only : pmid,pk,pedn,pek,am,byam
      USE RAD_COM, only : trhr,fsf,cosz1,trsurf
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,needtrs,trm,trmom,ntsurfsrc
     $     ,tr_mm, n_Be7, n_Be10
#ifdef TRACERS_DRYDEP
     *     ,dodrydep
#endif
#ifdef TRACERS_WATER
     *     ,nWATER,nGAS,nPART,tr_wd_TYPE,trname,trw0
#endif
#ifdef TRACERS_DUST
     &     ,Ntm_dust,n_clay
#endif
#endif
      USE PBLCOM, only : tsavg,dclev,eabl,uabl,vabl,tabl,qabl
      USE SOCPBL, only : npbl=>n
      USE PBL_DRV, only : pbl, t_pbl_args
      USE DIAG_COM, only : oa,aij=>aij_loc
     *     ,tdiurn,aj=>aj_loc,aregj=>aregj_loc,adiurn,ndiupt,jreg
     *     ,ij_tsli,ij_shdtli,ij_evhdt,ij_trhdt,ij_shdt,ij_trnfp0
     *     ,ij_srtr,ij_neth,ij_ws,ij_ts,ij_us,ij_vs,ij_taus,ij_tauus
     *     ,ij_tauvs,ij_qs,ij_tg1,ij_evap,ij_evapo,ij_tgo,ij_f0oc
     *     ,ij_f0oi,ij_evapi,ij_f0li,ij_evapli,j_evap,j_evhdt
     *     ,j_tsrf,j_shdt,j_trhdt,j_type,j_tg1,j_tg2,ijdd,idd_spr
     *     ,idd_pt5,idd_ts,idd_tg1
     *     ,idd_q5,idd_qs,idd_qg,idd_swg
     *     ,idd_lwg,idd_sh,idd_lh,idd_hz0,idd_ug,idd_vg,idd_wg,idd_us
     *     ,idd_vs,idd_ws,idd_cia,idd_cm,idd_ch,idd_cq,idd_eds,idd_dbl
     *     ,idd_ev,idd_ldc,idd_dcf,ij_pblht,ndiuvar,NREG,ij_dskin
     *     ,ij_gusti,ij_mccon,ij_sss,ij_trsup,ij_trsdn,ij_fwoc,ij_ssh
     *     ,adiurn_dust
#ifndef NO_HDIURN
     *     ,hdiurn
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
      USE LANDICE, only : z1e,z2li,hc1li,hc2li
      USE LANDICE_COM, only : snowli
      USE SEAICE, only : xsi,ace1i,alami,byrli,byrls,solar_ice_frac
      USE SEAICE_COM, only : rsi,msi,snowi,flag_dsws
      USE LAKES_COM, only : mwl,gml,flake
      USE LAKES, only : minmld
      USE COSMO_SOURCES, only : BE7D_acc 
      USE FLUXES, only : dth1,dq1,e0,e1,evapor,runoe,erunoe,sss
     *     ,solar,dmua,dmva,gtemp,nstype,uflux1,vflux1,tflux1,qflux1
     *     ,uosurf,vosurf,uisurf,visurf,ogeoza,gtempr
#ifdef TRACERS_ON
     *     ,trsrfflx,trsource
#ifdef TRACERS_GASEXCH_Natassa
     *     ,TRGASEX,GTRACER
#endif
#ifdef TRACERS_WATER
     *     ,trevapor,trunoe,gtracer
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
     &     ,trs_glob,pprec,pevap
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
      USE TRDIAG_COM, only : taijn=>taijn_loc, tajls=>tajls_loc,
     *      taijs=>taijs_loc,ijts_isrc,jls_isrc, jls_isrc, tij_surf,
     *      tij_surfbv, tij_gasx, tij_kw, tij_alpha, tij_evap,
     *      tij_grnd, tij_drydep, tij_gsdep
#ifdef TRACERS_DRYDEP
     *      , itcon_dd, dtr_dd
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only: DTR_AMPe
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE tracers_dust, only : hbaij,ricntd
#endif
#endif
      USE SOIL_DRV, only: earth

!@var DDMS downdraft mass flux in kg/(m^2 s), (i,j)
      USE CLOUDS_COM, only : DDMS
      
      IMPLICIT NONE

      INTEGER I,J,K,KR,JR,NS,NSTEPS,MODDSF,MODDD,ITYPE,IH,IHM,IDTYPE,IM1
     *     ,ii
      REAL*8 PLAND,PLICE,POICE,POCEAN,PIJ,PS,P1K
     *     ,ELHX,ACE2,CDTERM,CDENOM,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CM,CH,CQ,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG
     *     ,dT2,DQ1X,EVHDT0,EVAP,F0DT,FTEVAP,PWATER
     *     ,PSK,Q1,THV1,PTYPE,TG1,SRHEAT,SNOW,TG2
     *     ,SHDT,TRHDT,TG,TS,RHOSRF,RCDMWS,RCDHWS,RCDQWS,RCDHDWS,RCDQDWS
     *     ,SHEAT,TRHEAT,T2DEN,T2CON,T2MUL,FQEVAP,Z1BY6L,EVAPLIM,F2
     *     ,FSRI(2),HTLIM,dlwdt

      REAL*8 MA1, MSI1
      REAL*8, DIMENSION(NSTYPE,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *     TGRND,TGRN2,TGR4
      REAL*8, PARAMETER :: qmin=1.d-12
      REAL*8, PARAMETER :: Z2LI3L=Z2LI/(3.*ALAMI), Z1LIBYL=Z1E/ALAMI
      REAL*8 QSAT,DQSATDT,TR4
c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 hemi,qg_sat,dtsurf,uocean,vocean,qsrf,us,vs,ws,ws0
      logical pole
c
#ifdef TRACERS_ON
      real*8 totflux(ntm)
      integer n,nx,nsrc
      real*8, dimension(ntm) :: trs,trsfac,trconstflx
      integer ntix(ntm), ntx
      real*8, dimension(ntm) :: trgrnd
      real*8 trc_flux
#ifdef TRACERS_WATER
      real*8, dimension(ntm) :: tevaplim
      real*8  TEV,tevap,dTQS,TDP,TDT1,frac
#ifdef TRACERS_SPECIAL_O18
     *     ,FRACVL,FRACVS
#endif
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
      INTEGER,PARAMETER :: n_idxd=9+10*npbl
#endif
#endif

      REAL*8, DIMENSION(n_idx4,grid%J_STRT_HALO:grid%J_STOP_HALO,
     &     NDIUPT) :: diurn_part
      REAL*8,
     &     DIMENSION(n_idx3,grid%J_STRT_HALO:grid%J_STOP_HALO,NDIUPT)::
     &     diurn_partb
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      REAL*8,
     &     DIMENSION(n_idxd,grid%J_STRT_HALO:grid%J_STOP_HALO,NDIUPT) ::
     &     diurn_partd
#endif
      INTEGER :: idx1(n_idx1), idx2(n_idx2), idx3(n_idx3)
      INTEGER :: idx4(n_idx1+n_idx2)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER :: idxd(n_idxd)
#endif
      REAL*8 :: tmp(NDIUVAR)
      REAL*8, DIMENSION(n_idx4, NDIUPT) :: DIURNSUM
      REAL*8, DIMENSION(n_idx3, NDIUPT) :: DIURNSUMb
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      REAL*8,DIMENSION(n_idxd, NDIUPT) :: DIURNSUMd
#endif
C****
      INTEGER :: J_0, J_1, J_0H, J_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,        J_STOP=J_1)

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
     *       ,idd_ws2, idd_ustar, idd_us3, idd_stress, idd_lmon,
     *       idd_rifl,
     *       (idd_zpbl1+ii-1,ii=1,npbl), (idd_uabl1+ii-1,ii=1,npbl),
     *       (idd_vabl1+ii-1,ii=1,npbl), (idd_uvabl1+ii-1,ii=1,npbl),
     *       (idd_tabl1+ii-1,ii=1,npbl), (idd_qabl1+ii-1,ii=1,npbl),  
     *       (idd_zhat1+ii-1,ii=1,npbl-1), (idd_e1+ii-1,ii=1,npbl-1),
     *       (idd_km1+ii-1,ii=1,npbl-1), (idd_ri1+ii-1,ii=1,npbl-1),
     *       idd_grav,   idd_turb
#endif
     &       /)
      END IF
#endif

C****

      NSTEPS=NIsurf*ITime
      DTSURF=DTsrc/NIsurf
      IH=JHOUR+1
      IHM = IH+(JDATE-1)*24

C**** INITIALIZE TGRND: THIS IS USED TO UPDATE T OVER SURFACE STEPS
      DO J=J_0,J_1
      DO I=1,IM
        TGRND(2,I,J)=GTEMP(1,2,I,J)
        TGRND(3,I,J)=GTEMP(1,3,I,J)
        TGRN2(2,I,J)=GTEMP(2,2,I,J)
        TGRN2(3,I,J)=GTEMP(2,3,I,J)
        TGR4(2,I,J)=GTEMPR(2,I,J)**4
        TGR4(3,I,J)=GTEMPR(3,I,J)**4
      END DO
      END DO


C**** Zero out fluxes summed over type and surface time step
      E0=0. ; E1=0. ; EVAPOR=0. ; RUNOE=0. ; ERUNOE=0.
      DMUA=0. ; DMVA=0. ; SOLAR=0.
#ifdef TRACERS_WATER
      TREVAPOR = 0. ; TRUNOE = 0.
#endif
#ifdef TRACERS_DRYDEP
      TRDRYDEP = 0. ; dtr_dd=0.
#endif
#ifdef TRACERS_AMP
      DTR_AMPe(J_0:J_1,:) = 0.d0
#endif
C****
C**** OUTSIDE LOOP OVER TIME STEPS, EXECUTED NIsurf TIMES EVERY HOUR
C****
      DO NS=1,NIsurf
         MODDSF=MOD(NSTEPS+NS-1,NDASF*NIsurf+1)
         IF(MODDSF.EQ.0) IDACC(3)=IDACC(3)+1
         MODDD=MOD(1+ITime/NDAY+NS,NIsurf)   ! 1+ not really needed ??
C**** ZERO OUT FLUXES ACCUMULATED OVER SURFACE TYPES
         DTH1=0. ;  DQ1 =0. ;  uflux1=0. ; vflux1=0.
#ifdef TRACERS_ON
         trsrfflx = 0.
#ifdef TRACERS_GASEXCH_Natassa
         trgasex(:,:,:,:) = 0.
#endif
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

      diurn_part=0
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      if (adiurn_dust == 1) diurn_partd=0.D0
#endif

      Call HALO_UPDATE(GRID, vosurf, FROM=SOUTH)
      Call HALO_UPDATE(GRID, uisurf, FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID, visurf, FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID,u,FROM=SOUTH+NORTH)
      Call HALO_UPDATE(GRID,v,FROM=SOUTH+NORTH)

C****
C**** OUTSIDE LOOP OVER J AND I, EXECUTED ONCE FOR EACH GRID POINT
C****
!$OMP   PARALLEL DO PRIVATE (ACE2, CM,CH,CQ,
!$OMP*  CDTERM,CDENOM,DSHDTG,DQGDTG,DEVDTG,DTRDTG,
!$OMP*  DF0DTG,DFDTG,DTG,DQ1X,DF1DTG,DSNDTG,
!$OMP*  DT2, EVAP,EVAPLIM,ELHX,EVHDT,EVHEAT,EVHDT0,
!$OMP*  F0DT,F1DT,F0,F1,F2,FSRI, HCG1,HCG2,
!$OMP*  HTLIM,I,ITYPE,IDTYPE,IM1, J,K,
!$OMP*  KR, MA1,MSI1, PS,P1K,PLAND,PWATER,
!$OMP*  PLICE,PIJ,POICE,POCEAN,PTYPE,PSK, Q1,
!$OMP*  RHOSRF,RCDMWS,RCDHWS,RCDQWS,RCDHDWS,RCDQDWS, SHEAT,SRHEAT,
!$OMP*  SNOW,SHDT, T2DEN,T2CON,T2MUL,TS,
!$OMP*  THV1,TG,TG1,TG2,TR4,TRHDT,TRHEAT,Z1BY6L,dlwdt,
!$OMP*  HEMI,POLE,UOCEAN,VOCEAN,QG_SAT,US,VS,WS,WS0,QSRF,pbl_args,jr,tmp
#if defined(TRACERS_ON)
!$OMP*  ,n,nx,nsrc,totflux,trs,trsfac,trconstflx,trgrnd
!$OMP*  ,trc_flux
#if defined(TRACERS_WATER)
!$OMP*  ,tevaplim,tevap,TEV,dTQS,TDP,TDT1,frac
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

C**** Start loop over grid points
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
      JR=JREG(I,J)
      MA1=AM(1,I,J) !@var MA1 mass of lowest atmospheric layer (kg/m^2)

#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
C**** Calculate first layer tracer concentration
          pbl_args%trtop(nx)=trm(i,j,1,n)*byam(1,i,j)*bydxyp(j)
        end if
      end do
#endif

C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
         IF(MODDD.EQ.0) THEN
         DO KR=1,NDIUPT
           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             tmp(IDD_SPR)=PS
             do ii=1,5
               tmp(IDD_PT5+ii-1)=PSK*T(I,J,ii)
               tmp(IDD_Q5+ii-1) =Q(I,J,ii)
             end do
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
C****
C**** OPEN OCEAN/LAKE
C****
      if ( ITYPE == 1 ) then

      PTYPE=POCEAN
      IF (PTYPE.gt.0) THEN
      TG1=GTEMP(1,1,I,J)
      TG2=GTEMP(2,1,I,J)   ! diagnostic only
      TR4=GTEMPR(1,I,J)**4
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
          if (j==1) then
            uocean = uosurf(im,1)   ; vocean = uosurf(ivsp,1)
          else if (j==jm) then
            uocean = uosurf(im,jm)  ; vocean = uosurf(ivnp,jm)
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
      ELHX=LHE
      END IF
C****
C**** OCEAN/LAKE ICE
C****
      else if ( ITYPE == 2 ) then

      PTYPE=POICE
      IF (PTYPE.gt.0) THEN
      IF (FLAKE(I,J).gt.0) THEN
        IDTYPE=ITLKICE
        uocean = 0. ; vocean = 0. ! no dynamic ice for lakes
      ELSE
        IDTYPE=ITOICE
        if (UOdrag.eq.1) then ! use ice velocities in drag calculation
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
      TR4=TGR4(2,I,J)
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
      ELHX=LHS

      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)

      END IF
C****
C**** LAND ICE
C****
      else if ( ITYPE == 3 ) then

      PTYPE=PLICE
      IF (PTYPE.gt.0) THEN
      IDTYPE=ITLANDI
      SNOW=SNOWLI(I,J)
      TG1=TGRND(3,I,J)
      TG2=GTEMP(2,3,I,J)
      TR4=TGR4(3,I,J)
      SRHEAT=FSF(ITYPE,I,J)*COSZ1(I,J)
      Z1BY6L=(Z1LIBYL+SNOW*BYRLS)*BY6
      CDTERM=TG2
      CDENOM=1./(2.*Z1BY6L+Z2LI3L)
      HCG1=HC1LI+SNOW*SHI
      ELHX=LHS
      uocean = 0. ; vocean = 0. ! no land ice velocity
      END IF

      endif
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
      IF (ITYPE.eq.1 .and. focean(i,j).gt.0) QG_SAT=0.98d0*QG_SAT
      pbl_args%TG=TG   ! actual ground temperature
      pbl_args%TR4=TR4 ! radiative temperature K^4
      pbl_args%ELHX=ELHX   ! relevant latent heat
      pbl_args%QSOL=SRHEAT   ! solar heating
      pbl_args%TGV=TG*(1.+QG_SAT*deltx)
      IF (ITYPE.EQ.1) pbl_args%sss_loc=sss(I,J)

#ifdef TRACERS_ON
C**** Set up b.c. for tracer PBL calculation if required
      do nx=1,ntx
        n=ntix(nx)
C**** set defaults
        trsfac(nx)=0.
        totflux(nx)=0.
        trconstflx(nx)=0.
#ifdef TRACERS_GASEXCH_Natassa
       IF (ITYPE.EQ.1 .and. focean(i,j).gt.0.) THEN  ! OCEAN
          pbl_args%alati=sss(I,J)
          trgrnd(nx)=gtracer(n,itype,i,j)
          trsfac(nx)=1.
          trconstflx(nx)=trgrnd(nx)
       END IF
#endif
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
#ifdef TRACERS_WATER
        pbl_args%tr_evap_max(nx)=1.
C**** This distinguishes water from gases or particle
        if ( tr_wd_TYPE(n) == nWATER ) then
          trgrnd(nx)=gtracer(n,itype,i,j)
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
          do nsrc=1,ntsurfsrc(n)
            totflux(nx) = totflux(nx)+trsource(i,j,nsrc,n)
          end do
          trconstflx(nx)=totflux(nx)*bydxyp(j)   ! kg/m^2/s

#ifdef TRACERS_WATER
        endif
#endif
#ifdef TRACERS_GASEXCH_CO2_Natassa
        !need to redo this here because the previous line has changed trconstflx to zero.
        !because we have no sources. is there a better way to do this?
        trconstflx(nx)=trgrnd(nx)
#endif
      end do
#endif
C =====================================================================
      pbl_args%dtsurf = dtsurf
      pbl_args%TKV=THV1*PSK     ! TKV is referenced to the surface pressure
      pbl_args%ZS1=.5d-2*RGAS*pbl_args%TKV*MA1/PMID(1,I,J)
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
      pbl_args%ocean = (ITYPE.eq.1 .and. FOCEAN(I,J).gt.0)
#ifdef TRACERS_ON
      pbl_args%trs(:) = trs(:)
      pbl_args%trsfac(:) = trsfac(:)
      pbl_args%trconstflx(:) = trconstflx(:)
      pbl_args%ntix(:) = ntix(:)
      pbl_args%ntx = ntx
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      pbl_args%hbaij=hbaij(i,j)
      pbl_args%ricntd=ricntd(i,j)
      pbl_args%pprec=pprec(i,j)
      pbl_args%pevap=pevap(i,j,itype)
#endif

C**** Call pbl to calculate near surface profile
      CALL PBL(I,J,ITYPE,PTYPE,pbl_args)

#ifdef TRACERS_ON
      trs(:) = pbl_args%trs(:)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
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
      TS=pbl_args%TSV/(1.+QSRF*deltx)
C =====================================================================
C**** Adjust ground variables to account for skin effects
      TG = TG + pbl_args%dskin
      QG_SAT=QSAT(TG,ELHX,PS)
      IF (pbl_args%ocean) QG_SAT=0.98d0*QG_SAT
      TG1 = TG - TF
      TR4=(sqrt(sqrt(TR4))+pbl_args%dskin)**4

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
      TRHEAT=TRHR(0,I,J)-STBO*(TG*TG)*(TG*TG)
c      TRHEAT=TRHR(0,I,J)-STBO*TR4

C**** CASE (1) ! FLUXES USING EXPLICIT TIME STEP FOR OCEAN POINTS
      if ( ITYPE == 1) then
        SHDT = DTSURF*SHEAT
        EVHDT=DTSURF*EVHEAT              ! latent heat flux
        TRHDT=DTSURF*TRHEAT

C**** CASE (2) ! FLUXES USING IMPLICIT TIME STEP FOR ICE POINTS
      else if ( ITYPE == 2 ) then

! heat flux on first/second/third layers (W/m^2)
        F1 = (TG1-TG2)*dF1dTG + SRHEAT*FSRI(1)
        F2 = SRHEAT*FSRI(2)
        ! Including gustiness in the latent heat flux:
        EVHEAT=LHE*(RCDQWS*(QSRF-QG_SAT)+RCDQDWS*pbl_args%qprime) ! why is this different to above?
        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        dSNdTG=-RCDHWS*SHA
        dQGdTG=QG_SAT*DQSATDT(TG,ELHX) ! d(QG)/dTG
        dEVdTG = -dQGdTG*LHE*RCDQWS ! d(EVHEAT)/dTG
        dTRdTG = -4*STBO*TG*TG*TG ! d(TRHEAT)/dTG
c        dTRdTG = -4*STBO*sqrt(sqrt(TR4))**3 ! d(TRHEAT)/dTG
        dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG

        T2DEN = HCG2+DTSURF*dF1dTG
        T2CON = DTSURF*(F1-F2)/T2DEN
        T2MUL = DTSURF*dF1dTG/T2DEN

        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)

        IF (TG1+dTG .GT. 0.) dTG = -TG1
        dT2 = T2CON+T2MUL*dTG

        SHDT  = DTSURF*(SHEAT +dTG*dSNdTG) ! sensible
        EVHDT = DTSURF*(EVHEAT+dTG*dEVdTG) ! latent
        TRHDT = DTSURF*(TRHEAT+dTG*dTRdTG) ! thermal flux (J/m^2)
        F1DT = DTSURF*(F1+(dTG*dF1dTG-dT2*dF1dTG))
        TG1 = TG1+dTG          ! first layer sea ice temperature (degC)
        TG2 = TG2+dT2          ! second layer sea ice temperature (degC)
        TGRN2(ITYPE,I,J) = TG2

C**** CASE (3) ! FLUXES USING IMPLICIT TIME STEP OVER LANDICE
      else if ( ITYPE == 3 ) then

        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        F1=(TG1-CDTERM-F0*Z1BY6L)*CDENOM
        DSHDTG=-RCDHWS*SHA
        DQGDTG=QG_SAT*DQSATDT(TG,ELHX)
        DEVDTG=-RCDQWS*LHE*DQGDTG
        DTRDTG=-4.*STBO*TG*TG*TG
c        DTRDTG=-4.*STBO*sqrt(sqrt(TR4))**3
        DF0DTG=DSHDTG+DEVDTG+DTRDTG
        DFDTG=DF0DTG-(1.-DF0DTG*Z1BY6L)*CDENOM
        DTG=(F0-F1)*DTSURF/(HCG1-DTSURF*DFDTG)
        SHDT=DTSURF*(SHEAT+DTG*DSHDTG)
        EVHDT=DTSURF*(EVHEAT+DTG*DEVDTG)
        TRHDT=DTSURF*(TRHEAT+DTG*DTRDTG)
        F1DT=DTSURF*(TG1-CDTERM-(F0+DTG*DFDTG)*Z1BY6L)*CDENOM
        TG1=TG1+DTG

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
#ifdef TRACERS_SPECIAL_O18
            TEV=-(RCDQWS*(trs(nx)-trgrnd(nx)*QG_SAT*fracvl(tg1,n))
     *           +RCDQDWS*pbl_args%trprime(nx))*pbl_args%frack(nx)
#else
            TEV=-(RCDQWS*(trs(nx)-trgrnd(nx)*QG_SAT)
     *           +RCDQDWS*pbl_args%trprime(nx))
#endif
            TEVAP=DTSURF*TEV
          ELSE                  ! ICE AND LAND ICE
C**** tracer flux is set by source tracer concentration
            IF (EVAP.GE.0) THEN ! EVAPORATION
              TEVAP=EVAP*trgrnd(nx)
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
#ifdef TRACERS_GASEXCH_Natassa
C****
C**** Calculate Tracer Gas Exchange
C****
       IF (ITYPE.EQ.1 .and. focean(i,j).gt.0.) THEN  ! OCEAN
#ifdef TRACERS_GASEXCH_CFC_Natassa
          TRGASEX(n,ITYPE,I,J) =
     .        pbl_args%Kw_gas * (pbl_args%beta_gas*trs(nx)-trgrnd(nx))
          trsrfflx(i,j,n)=trsrfflx(i,j,n)
     .         +pbl_args%Kw_gas * (pbl_args%beta_gas*trs(nx)-trgrnd(nx))
     .               * dxyp(j)*ptype
         taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n))
     .         +pbl_args%Kw_gas * (pbl_args%beta_gas*trs(nx)-trgrnd(nx))
     .               * dxyp(j)*ptype*dtsurf
#endif
#ifdef TRACERS_GASEXCH_CO2_Natassa
          !this is modeled in complete accordance to what Watson is doing
          TRGASEX(n,ITYPE,I,J) =
     .        pbl_args%Kw_gas * (pbl_args%beta_gas*trs(nx)-
     .                           pbl_args%alpha_gas*1.024e-3*trgrnd(nx))
          trsrfflx(i,j,n)=trsrfflx(i,j,n)
     .         +pbl_args%Kw_gas * (pbl_args%beta_gas*trs(nx)-
     .                           pbl_args%alpha_gas*1.024e-3*trgrnd(nx))
     .               * dxyp(j)*ptype
         taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n))
     .         +pbl_args%Kw_gas * (pbl_args%beta_gas*trs(nx)-
     .                           pbl_args%alpha_gas*1.024e-3*trgrnd(nx))
     .               * dxyp(j)*ptype*dtsurf
#endif
       END IF
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
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
        case default
          trc_flux=0
        end select

        trsrfflx(i,j,n)=trsrfflx(i,j,n)+
     &       trc_flux*dxyp(j)*ptype
        taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &       trc_flux*dxyp(j)*ptype*dtsurf

#ifdef TRACERS_AMP
        DTR_AMPe(j,n)=DTR_AMPe(j,n)+trc_flux*dxyp(j)*ptype*dtsurf
#else
        tajls(j,1,jls_isrc(1,n)) = tajls(j,1,jls_isrc(1,n))+
     *       trc_flux*dxyp(j)*ptype*dtsurf   ! why not for all aerosols?
#endif
#endif

#ifdef TRACERS_DRYDEP
C****
C**** Calculate Tracer Dry Deposition (including gravitational settling)
C****
        if(dodrydep(n))then
          rts=rhosrf*trs(nx)
          rtsdt=rts*dtsurf                             ! kg*s/m^3
          tdryd=-rtsdt*(pbl_args%dep_vel(n)+pbl_args%gs_vel(n))          ! kg/m2
          tdd = tdryd*dxyp(j)*ptype                    ! kg
          td1 = (trsrfflx(i,j,n)+totflux(nx))*dtsurf   ! kg
          if (trm(i,j,1,n)+td1+tdd.le.0.and.tdd.lt.0) then
            if (qcheck) write(99,*) "limiting tdryd surfce",i,j,n,tdd
     *           ,trm(i,j,1,n),td1,trs(nx),pbl_args%trtop(nx)
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
     &         ptype*rtsdt*pbl_args%dep_vel(n)
          taijn(i,j,tij_gsdep ,n)=taijn(i,j,tij_gsdep ,n) +
     &         ptype*rtsdt* pbl_args%gs_vel(n)
          if (n .eq. n_Be7) then 
            BE7D_acc(i,j)=BE7D_acc(i,j)+ptype*rtsdt*pbl_args%dep_vel(n)
     *           +ptype*rtsdt* pbl_args%gs_vel(n)
          end if

          dtr_dd(j,n,1)=dtr_dd(j,n,1)-
     &         ptype*rtsdt*dxyp(j)*pbl_args%dep_vel(n)
          dtr_dd(j,n,2)=dtr_dd(j,n,2)-
     &         ptype*rtsdt*dxyp(j)* pbl_args%gs_vel(n)
        end if
#endif
      END DO
#endif
C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSURF*SRHEAT+TRHDT+SHDT+EVHDT
C**** Limit heat fluxes out of lakes if near minimum depth
      IF (ITYPE.eq.1 .and. FLAKE(I,J).gt.0 .and.
     *     E0(I,J,1)+F0DT+HTLIM.lt.0 .and. HTLIM.gt.0) THEN
        if (QCHECK) write(6,*) "Limiting heat flux from lake",i,j,SHDT
     *       ,F0DT,E0(I,J,1),DTSURF*SRHEAT,TRHDT,EVHDT,HTLIM
        SHDT = -(HTLIM+E0(I,J,1)+DTSURF*SRHEAT+TRHDT+EVHDT)
        F0DT = -E0(I,J,1)-HTLIM
        if (QCHECK) write(6,*) "New SHDT,F0DT",i,j,SHDT,F0DT
      END IF
      E0(I,J,ITYPE)=E0(I,J,ITYPE)+F0DT
      E1(I,J,ITYPE)=E1(I,J,ITYPE)+F1DT
      EVAPOR(I,J,ITYPE)=EVAPOR(I,J,ITYPE)+EVAP

      TGRND(ITYPE,I,J)=TG1  ! includes skin effects
      TGR4(ITYPE,I,J) =TR4
C**** calculate correction for different TG in radiation and surface
      dLWDT = DTSURF*(TRSURF(ITYPE,I,J)-STBO*(TG1+TF)**4)
c      dLWDT = DTSURF*(TRSURF(ITYPE,I,J)-STBO*TR4
C**** final fluxes
      DTH1(I,J)=DTH1(I,J)-(SHDT+dLWDT)*PTYPE/(SHA*MA1*P1K)  ! +ve up
      DQ1(I,J) =DQ1(I,J) -DQ1X*PTYPE
      DMUA(I,J,ITYPE)=DMUA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*(US-UOCEAN)
      DMVA(I,J,ITYPE)=DMVA(I,J,ITYPE)+PTYPE*DTSURF*RCDMWS*(VS-VOCEAN)
      uflux1(i,j)=uflux1(i,j)+PTYPE*RCDMWS*(US-UOCEAN)
      vflux1(i,j)=vflux1(i,j)+PTYPE*RCDMWS*(VS-VOCEAN)
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
        AREGJ(JR,J,J_TRHDT)=AREGJ(JR,J,J_TRHDT)+TRHDT*PTYPE*DXYP(J)
        AREGJ(JR,J,J_SHDT )=AREGJ(JR,J,J_SHDT )+SHDT *PTYPE*DXYP(J)
        AREGJ(JR,J,J_EVHDT)=AREGJ(JR,J,J_EVHDT)+EVHDT*PTYPE*DXYP(J)
        AREGJ(JR,J,J_EVAP )=AREGJ(JR,J,J_EVAP )+EVAP *PTYPE*DXYP(J)
        IF(MODDSF.EQ.0) THEN
          AREGJ(JR,J,J_TSRF)=AREGJ(JR,J,J_TSRF)+(TS-TF)*PTYPE*DXYP(J)
          AREGJ(JR,J,J_TG1 )=AREGJ(JR,J,J_TG1 )+    TG1*PTYPE*DXYP(J)
          AREGJ(JR,J,J_TG2 )=AREGJ(JR,J,J_TG2 )+    TG2*PTYPE*DXYP(J)
        END IF

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
          AIJ(I,J,IJ_TAUS)=AIJ(I,J,IJ_TAUS)+RCDMWS*WS*PTYPE
          AIJ(I,J,IJ_TAUUS)=AIJ(I,J,IJ_TAUUS)+RCDMWS*(US-UOCEAN)*PTYPE
          AIJ(I,J,IJ_TAUVS)=AIJ(I,J,IJ_TAUVS)+RCDMWS*(VS-VOCEAN)*PTYPE
          AIJ(I,J,IJ_QS)=AIJ(I,J,IJ_QS)+QSRF*PTYPE
          AIJ(I,J,IJ_TG1)=AIJ(I,J,IJ_TG1)+TG1*PTYPE
          AIJ(I,J,IJ_PBLHT)=AIJ(I,J,IJ_PBLHT)+pbl_args%dbl*PTYPE
          if(DDMS(I,J).lt.0.) ! ddms < 0 for down draft
     *         AIJ(I,J,ij_mccon)=AIJ(I,J,ij_mccon)+ptype
          AIJ(I,J,IJ_GUSTI)=AIJ(I,J,IJ_GUSTI)+pbl_args%gusti*PTYPE

          if (ITYPE==1)
     *         AIJ(I,J,IJ_DSKIN)=AIJ(I,J,IJ_DSKIN)+pbl_args%dskin

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          aij(i,j,ij_wsgcm)=aij(i,j,ij_wsgcm)+pbl_args%wsgcm*ptype
          aij(i,j,ij_wspdf)=aij(i,j,ij_wspdf)+pbl_args%wspdf*ptype
          aij(i,j,ij_wdry)=aij(i,j,ij_wdry)+pbl_args%wsubwd*ptype
          aij(i,j,ij_wtke)=aij(i,j,ij_wtke)+pbl_args%wsubtke*ptype
          aij(i,j,ij_wmoist)=aij(i,j,ij_wmoist)+pbl_args%wsubwm*ptype
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
          depo_turb_glob(i,j,itype,n)=ptype*rts*pbl_args%dep_vel(n)
          depo_grav_glob(i,j,itype,n)=ptype*rts*pbl_args%gs_vel(n)
        END IF
      END DO
#endif
#endif

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

              DIURN_part(n_idx1+1:n_idx4,J,kr)=
     &             DIURN_part(n_idx1+1:n_idx4,J,kr)+tmp(idx2(:))

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
                DO n=1,Ntm_dust
                  n1=n_clay+n-1
                  IF (dodrydep(n1)) THEN
                    tmp(idd_turb)=ptype*rts*pbl_args%dep_vel(n1)
                    tmp(idd_grav)=ptype*rts*pbl_args%gs_vel(n1)
                  END IF
                END DO
#endif
                tmp(idd_ws2)=ws*ws*ptype
                tmp(idd_ustar)=pbl_args%ustar*ptype
                tmp(idd_us3)=ptype*pbl_args%ustar**3
                tmp(idd_stress)=rcdmws*pbl_args%ws*ptype
                tmp(idd_lmon)=pbl_args%lmonin*ptype
                tmp(idd_rifl)=
     &               +ptype*grav*(ts-tg)*pbl_args%zgs/(ws*ws*tg)

                tmp(idd_zpbl1:idd_zpbl1+npbl-1)=ptype*pbl_args%z(1:npbl)
                tmp(idd_uabl1:idd_uabl1+npbl-1)=ptype*uabl(1:npbl,i,j
     *               ,itype)
                tmp(idd_vabl1:idd_vabl1+npbl-1)=ptype*vabl(1:npbl,i,j
     *               ,itype)
                tmp(idd_uvabl1:idd_uvabl1+npbl-1)=ptype*sqrt(
     *               uabl(1:npbl,i,j,itype)*uabl(1:npbl,i,j,itype)+
     *               vabl(1:npbl,i,j,itype)*vabl(1:npbl,i,j,itype))
                tmp(idd_tabl1:idd_tabl1+npbl-1)=ptype*tabl(1:npbl,i,j
     *               ,itype)
                tmp(idd_qabl1:idd_qabl1+npbl-1)=ptype*qabl(1:npbl,i,j
     *               ,itype)
                tmp(idd_zhat1:idd_zhat1+npbl-2)=ptype
     *               *pbl_args%zhat(1:npbl-1)
                tmp(idd_e1:idd_e1+npbl-2)=eabl(1:npbl-1,i,j,itype)*ptype
                tmp(idd_km1:idd_km1+npbl-2)=ptype*pbl_args%km(1:npbl-1)
                tmp(idd_ri1:idd_ri1+npbl-2)=ptype*pbl_args%gh(1:npbl-1)
     *               /(pbl_args%gm(1:npbl-1)+1d-20)
#endif
                DIURN_partd(:,J,kr)=DIURN_partd(:,J,kr)+tmp(idxd(:))
              END IF
#endif
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
            taijn(i,j,tij_surf  ,n) = taijn(i,j,tij_surf  ,n)
     *           +trs(nx)*ptype
            taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)
     *           +trs(nx)*ptype*rhosrf
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
            trs_glob(i,j,itype,n)=trs(nx)*ptype
#endif
          else
            taijn(i,j,tij_surf,n) = taijn(i,j,tij_surf,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *bydxyp(j),0d0)*ptype
            taijn(i,j,tij_surfbv,n) = taijn(i,j,tij_surfbv,n)
     *           +max((trm(i,j,1,n)-trmom(mz,i,j,1,n))*byam(1,i,j)
     *           *bydxyp(j),0d0)*ptype*rhosrf
          end if
#ifdef TRACERS_GASEXCH_Natassa
          if (focean(i,j).gt.0.) then
            taijn(i,j,tij_kw,n) = taijn(i,j,tij_kw,n)
     *                          + pbl_args%Kw_gas*ptype
            taijn(i,j,tij_alpha,n) = taijn(i,j,tij_alpha,n)
     *                             + pbl_args%alpha_gas*ptype
            taijn(i,j,tij_gasx,n) = taijn(i,j,tij_gasx,n)
     *                            + trgrnd(nx)*ptype
          endif
#endif
#ifdef TRACERS_WATER
          if (tr_wd_type(n).eq.nWater) then
            taijn(i,j,tij_evap,n)=taijn(i,j,tij_evap,n)+
     *           trevapor(n,itype,i,j)*ptype
            tajls(j,1,jls_isrc(1,n))=tajls(j,1,jls_isrc(1,n))
     *           +trevapor(n,itype,i,j)*ptype
            if (focean(i,j).gt.0) tajls(j,1,jls_isrc(2,n))=tajls(j,1
     *           ,jls_isrc(2,n))+trevapor(n,itype,i,j)*ptype
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
      END DO   ! end of itype loop
      IM1=I
      END DO   ! end of I loop

      END DO   ! end of J loop
!$OMP  END PARALLEL DO

      CALL GLOBALSUM(grid, DIURN_part, DIURNSUM)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IF (adiurn_dust == 1)
     &   CALL globalsum(grid,diurn_partd,diurnsumd)
#endif

      IF (AM_I_ROOT()) THEN
         ADIURN(ih,idx4,:)=ADIURN(ih,idx4,:)   + DIURNSUM
#ifndef NO_HDIURN
         HDIURN(ihm,idx4,:)=HDIURN(ihm,idx4,:) + DIURNSUM
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
         IF (adiurn_dust == 1) THEN
           adiurn(ih,idxd,:)=adiurn(ih,idxd,:)+diurnsumd
#ifndef NO_HDIURN
           hdiurn(ihm,idxd,:)=hdiurn(ihm,idxd,:)+diurnsumd
#endif
         END IF
#endif
      END IF


C****
C**** dynamic vegetation time step
C****
!!! probably don't need this call unless something can be done
!   separately from ground hydrology on i,j grid
!      call step_dveg(dtsurf)
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
          QMOM(:,I,J,1)=0.
        ENDIF
c****   retrieve fluxes
        P1K=PK(1,I,J)
        tflux1(i,j)=-dth1(i,j)*AM(1,I,J)*P1K/(dtsurf)
        qflux1(i,j)=-dq1(i,j)*AM(1,I,J)/(dtsurf)
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
              tmp(1)=1.
              tmp(2)=DCLEV(I,J)
              DIURN_partb(:,J,KR)=DIURN_partb(:,J,KR)+tmp(1:2)
            END IF
          END IF
       END DO

       CALL GLOBALSUM(grid,  DIURN_partb, DIURNSUMb)

       IF (AM_I_ROOT()) THEN
          ADIURN(ih,idx3,:)=ADIURN(ih,idx3,:)   + DIURNSUMb
#ifndef NO_HDIURN
          HDIURN(ihm,idx3,:)=HDIURN(ihm,idx3,:) + DIURNSUMb
#endif
       END IF

      END IF
C****
      END DO   ! end of surface time step

#ifdef TRACERS_DRYDEP
C**** Save for tracer dry deposition conservation quantity:
      do n=1,ntm
        if(dodrydep(n)) then
          if (itcon_dd(n,1).gt.0)
     *    call diagtcb(dtr_dd(:,n,1),itcon_dd(n,1),n)  ! turb dep
          if (itcon_dd(n,2).gt.0)
     *         call diagtcb(dtr_dd(:,n,2),itcon_dd(n,2),n) ! grav sett
        end if
      end do
#endif

      RETURN
C****
      END SUBROUTINE SURFCE
