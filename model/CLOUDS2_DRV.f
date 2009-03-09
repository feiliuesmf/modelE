#include "rundeck_opts.h"

      SUBROUTINE CONDSE
!@sum   CONDSE driver for moist convection AND large-scale condensation
!@auth  M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver   1.0 (taken from CB265)
!@calls CLOUDS:MSTCNV,CLOUDS:LSCOND
      USE CONSTANT, only : bygrav,lhm,rgas,grav,tf,lhe,lhs,sha,deltx
     *     ,teeny,sday
      USE MODEL_COM, only : im,jm,lm,p,u,v,t,q,wm,JHOUR
     *     ,ls1,psf,ptop,dsig,bydsig,sig,DTsrc,ftype,jdate
     *     ,ntype,itime,focean,fland,flice,jyear,jmon
#ifdef SCM
     &     ,I_TARG,J_TARG,NSTEPSCM
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID,GET,AM_I_ROOT
      USE DOMAIN_DECOMP_ATM, only : GLOBALSUM
      USE QUSDEF, only : nmom
      USE SOMTQ_COM, only : t3mom=>tmom,q3mom=>qmom
      USE GEOM, only : imaxj,axyp,byaxyp, kmaxj
#ifndef CUBE_GRID
      USE GEOM, only : ravj
#endif
      USE RANDOM
      USE RAD_COM, only : cosz1
      USE CLOUDS_COM, only : ttold,qtold,svlhx,svlat,rhsav,cldsav
     &     ,isccp_reg2d,ukm,vkm
#ifdef CLD_AER_CDNC
     *     ,oldnl,oldni,ctem,cd3d,cl3d,ci3d,clwp,cdn3d,cre3d  ! for 3 hrly diag
#endif
     *     ,tauss,taumc,cldss,cldmc,csizmc,csizss,fss,cldsav1
     *     ,tls,qls,tmc,qmc,ddm1,airx,lmc
     *     ,ddms,tdn1,qdn1,ddml
      USE DIAG_COM, only : aij=>aij_loc,
     *     aijl=>aijl_loc,adiurn=>adiurn_loc,jreg,ij_pscld,
     *     ij_pdcld,ij_scnvfrq,ij_dcnvfrq,ij_wmsum,ij_snwf,ij_prec,
     *     ij_neth,ij_f0oc,j_eprcp,j_prcpmc,j_prcpss,ijl_mc,
     *     ijdd,idd_pr,idd_ecnd,idd_mcp,idd_dmc,idd_smc,idd_ssp,
     &     jl_mcmflx,jl_sshr,jl_mchr,jl_dammc,jl_rhe,jl_mchphas,
     *     jl_mcdtotw,jl_mcldht,jl_mcheat,jl_mcdry,ij_ctpi,ij_taui,
     *     ij_lcldi,ij_mcldi,ij_hcldi,ij_tcldi,ij_sstabx,isccp_diags,
     *     ndiupt,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss,ij_scldi,
     *     jl_mcshlw,jl_mcdeep,
#ifndef NO_HDIURN
     *     hdiurn=>hdiurn_loc,
#endif
     *     ntau,npres,aisccp=>aisccp_loc,ij_precmc,ij_cldw,ij_cldi,
     *     ij_fwoc,p_acc,pm_acc,ndiuvar,nisccp,adiurn_dust,jl_mcdflx
     *     ,lh_diags,ijl_llh,ijl_mctlh,ijl_mcdlh,ijl_mcslh
#ifdef CLD_AER_CDNC
     *     ,jl_cnumwm,jl_cnumws,jl_cnumim,jl_cnumis
     *     ,ij_3dnwm,ij_3dnws,ij_3dnim,ij_3dnis
     *     ,ij_3drwm,ij_3drws,ij_3drim,ij_3dris
     *     ,ij_3dlwm,ij_3dlws,ij_3dlim,ij_3dlis
     *     ,ijl_rewm,ijl_rews,ijl_cdwm,ijl_cdws,ijl_cwwm,ijl_cwws
     *     ,ij_wmclwp,ij_wmctwp
     *     ,ijl_reim,ijl_reis,ijl_cdim,ijl_cdis,ijl_cwim,ijl_cwis
#endif
#ifdef TRACERS_DUST
     &     ,idd_wet
#endif
#ifdef HTAP_LIKE_DIAGS
     &     ,IJ_MCamFX
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
      USE AERO_CONFIG, only: NMODES
      USE AMP_AEROSOL, only: NACTV, NACTC
#endif
#endif
#ifdef TRACERS_ON
      USE TRACER_COM, only: itime_tr0,TRM,TRMOM,NTM,trname,trdn1
#ifdef TRACERS_COSMO
     *     ,n_Be10,n_Be7
#endif
#ifdef TRACERS_DUST
     *     ,n_clay,n_clayilli,n_sil1quhe
#endif
#ifdef TRACERS_WATER
     *     ,trwm,trw0,dowetdep
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,Ntm_dust
#endif
#ifdef TRACERS_DUST
     &     ,imDust
#else
#endif
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7W_acc
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE LIGHTNING, only : RNOx_lgt
#endif
#ifndef SKIP_TRACER_DIAGS
      USE TRDIAG_COM,only: jlnt_mc,jlnt_lscond,itcon_mc
     *     ,itcon_ss,taijn=>taijn_loc,taijs=>taijs_loc
#ifdef TRACERS_WATER
     *     ,jls_prec,tij_prec
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     *     ,jls_incloud,ijts_aq
#endif
#ifdef TRDIAG_WETDEPO
     &     ,jls_trdpmc,jls_trdpls,ijts_trdpmc,ijts_trdpls
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,jls_wet,ijts_wet,itcon_wt
#endif
#endif
#endif /*SKIP_TRACER_DIAGS*/
      USE CLOUDS, only : tm,tmom,trdnl ! local  (i,j)
     *     ,ntx,ntix              ! global (same for all i,j)
#ifdef TRACERS_WATER
     *     ,trwml,trsvwml,trprmc,trprss
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
     *     ,dt_sulf_mc,dt_sulf_ss
#endif
#ifdef TRDIAG_WETDEPO
     &     ,trcond_mc,trdvap_mc,trflcw_mc,trprcp_mc,trnvap_mc,trwash_mc
     &     ,trwash_ls,trevap_ls,trclwc_ls,trprcp_ls,trclwe_ls,trcond_ls
     &     ,diag_wetdep
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,tm_dust,tmom_dust,trprc_dust
#endif
#endif
#endif
      USE CLOUDS, only : BYDTsrc,mstcnv,lscond ! glb var & subs
     *     ,airm,byam,etal,sm,smom,qm,qmom,isc,dxypij,lp50,hcndss
     *     ,tl,ris,ri1,ri2,mcflx,sshr,dgdsm,dphase,dtotw,dqcond,dctei
     *     ,wml,sdl,u_0,v_0,um,vm,um1,vm1,qs,us,vs,dcl,airxl,prcpss
     *     ,prcpmc,pearth,ts,taumcl,cldmcl,svwmxl,svlatl,svlhxl,dgdqm
     *     ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum
     *     ,aq,dpdt,th,ql,wmx,ttoldl,rh,taussl,cldssl,cldsavl,rh1
     *     ,kmax,ra,pl,ple,plk,rndssl,lhp,debug,fssl,pland,cldsv1
     *     ,smommc,smomls,qmommc,qmomls,ddmflx,wturb,ncol
     *     ,tvl,w2l,gzl,savwl,savwl1,save1l,save2l
     *     ,dphashlw,dphadeep,dgshlw,dgdeep,tdnl,qdnl,prebar1
#ifdef CLD_AER_CDNC
     *     ,acdnwm,acdnim,acdnws,acdnis,arews,arewm,areis,areim
     *     ,alwim,alwis,alwwm,alwws,nlsw,nlsi,nmcw,nmci
     *     ,oldcdl,oldcdi,sme
     *     ,cteml,cd3dl,cl3dl,ci3dl,cdn3dl,cre3dl,smlwp
     *     ,wmclwp,wmctwp
#endif
#ifdef SCM
      USE SCMCOM , only : SCM_SAVE_Q,SCM_SAVE_T,SCM_DEL_Q,SCM_DEL_T,
     *                    SCM_ATURB_FLAG,iu_scm_prt
      USE SCMDIAG , only : WCUSCM,WCUALL,WCUDEEP,PRCCDEEP,NPRCCDEEP,
     &                    MPLUMESCM,MPLUMEALL,MPLUMEDEEP,
     &                    ENTSCM,ENTALL,ENTDEEP,
     &                    DETRAINDEEP,TPALL,PRCSS,PRCMC
#endif
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg,tgvavg,qgavg,dclev,egcm
     *  ,w2gcm
      USE DYNAMICS, only : pk,pek,pmid,pedn,sd_clouds,gz,ptold,pdsig,sda
     *     ,ltropo
     &     ,ua=>ualij,va=>valij
      USE SEAICE_COM, only : rsi
      USE GHY_COM, only : snoage,fearth
      USE LAKES_COM, only : flake
      USE FLUXES, only : prec,eprec,precss,gtempr
#ifdef TRACERS_WATER
     *     ,trprec
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,trprec_dust
#endif
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only : AQsulfRATE
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      use tracer_sources, only : n__prec
#endif
      USE FILEMANAGER, only: openunit,closeunit
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE tracers_dust,ONLY : prelay
#endif
      IMPLICIT NONE

#ifdef TRACERS_ON
!@var tmsave holds tracer value (for diagnostics)
      REAL*8 tmsave(lm,ntm),tmomsv(nmom,lm,ntm)
      INTEGER NX
#ifdef TRACERS_SPECIAL_Shindell
!@var Lfreeze Lowest level where temperature is below freezing (TF)
      INTEGER Lfreeze
#endif
#endif

!@var TLS,QLS,TMC,QMC temperature and humidity work arrays
!@var FSS fraction of the grid box for large-scale cloud
C      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
c     *           TLS,QLS,TMC,QMC

      REAL*8, DIMENSION(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO)
     &     :: UASV              ! for U tendency diagnostic

!@param ENTCON fractional rate of entrainment (km**-1)
      REAL*8,  PARAMETER :: ENTCON = .2d0

      INTEGER I,J,K,L,N,LL  !@var I,J,K,L,N loop variables
      INTEGER JR,KR,ITYPE,IT,IH,LP850,LP600,IHM,KMAX_NONPOLAR
!@var JR = JREG(I,J)
!@var KR index for regional diagnostics
!@var ITYPE index for snow age
!@var IT index for surface types
!@var LP850 layer near 850 mb
!@var LP600 layer near 600 mb
!@var LERR,IERR error reporting
      INTEGER :: LERR, IERR

      REAL*8 :: HCNDMC,PRCP,TPRCP,EPRCP,ENRGP,WMERR,ALPHA1,ALPHA2,ALPHAS
      REAL*8 :: DTDZ,DTDZS,DUDZ,DVDZ,DUDZS,DVDZS,THSV,THV1,THV2,QG,TGV
      REAL*8 :: DH1S,BYDH1S,DH12,BYDH12,DTDZG,DUDZG,DVDZG,SSTAB,DIFT,CSC
     *     ,TSV,WM1,WMI
cECON*     ,E,E1,W1,ep,ep1,TSV,q0,q1,q2,WM1,WMI
!@var HCNDMC heating due to moist convection
!@var PRCP precipitation
!@var TPRCP temperature of mc. precip  (deg. C)
!@var EPRCP sensible heat of precip
!@var ENRGP total energy of precip
!@var WMERR DH12,BYDH12,DH1S,BYDH1S,SSTAB dummy variable
!@var THSV,THV1,THV2 vertual potential temperatures
!@var QG,TGV ground humidity,virt.temperature from pbl
!@var ALPHA1,ALPHA2,ALPHAS,DIFT,CSC dummy variables
!@var DTDZ,DTDZS,DTDZG vertical potential temperature gradients
!@var DUDZ,DVDZ,DUDZS,DVDZS,DUDZG,DVDZG vertical wind gradients
!@var TSV virtual surface temperature (K)

C**** parameters and variables for isccp diags
      real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
      real*8 skt,conv(lm),qv(lm)
      real*8 pfull(lm),at(lm),cc(lm),dtau_s(lm),dtau_c(lm)
      real*8 dem_s(lm),dem_c(lm),phalf(lm+1)
      real*8 fq_isccp(ntau,npres),ctp,tauopt
      integer itau,itrop,nbox,sunlit,ipres
C****

C
Cred*                       Reduced Arrays 1                 *********
C        not clear yet whether they still speed things up
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,LM) ::
     &     GZIL,SD_CLDIL,WMIL
      real*8, dimension(NMOM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,LM) ::
     &     TMOMIL,QMOMIL
      INTEGER ICKERR, JCKERR, JERR, seed, NR
      REAL*8  RNDSS(3,LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                   GRID%J_STRT_HALO:GRID%J_STOP_HALO),xx
      integer :: nij_before_j0,nij_after_j1,nij_after_i1
      REAL*8  UKMSP(IM,LM), VKMSP(IM,LM), UKMNP(IM,LM),VKMNP(IM,LM)
      REAL*4 WCU500(IM,16),SAVWCU(IM,16,LM),SAVEN1(IM,16,LM),
     *  SAVEN2(IM,16,LM),W500P1(16),ENTJ(16),SAVWC1(IM,16,LM)
      INTEGER :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      INTEGER, PARAMETER :: n_idx1 = 5
      INTEGER, PARAMETER :: n_idx2 = 3
      INTEGER, PARAMETER :: n_idx3 = 6
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idxd=1
#endif
      INTEGER :: idx1(n_idx1), idx2(n_idx2), idx3(n_idx3)
#ifdef TRACERS_DUST
      INTEGER :: idxd(n_idxd)
#endif
      REAL*8 :: tmp(NDIUVAR)
#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER :: n1,n_fidx
#endif
#endif

C**** Initialize
#ifdef TRACERS_SPECIAL_Shindell
      RNOx_lgt(:,:)=0.d0
#endif
      idx1 = (/ IDD_PR, IDD_ECND, IDD_MCP, IDD_DMC, IDD_SMC /)
      idx2 = (/ IDD_PR, IDD_ECND, IDD_SSP /)
      idx3 = (/ IDD_PR, IDD_ECND, IDD_MCP, IDD_DMC, IDD_SMC, IDD_SSP /)
#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) idxd=(/idd_wet/)
#endif


C**** define local grid
      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_HALO=J_0H,    J_STOP_HALO=J_1H,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE        )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C
C     OBTAIN RANDOM NUMBERS FOR PARALLEL REGION
C
C     Burn some random numbers corresponding to latitudes off
C     processor
      CALL BURN_RANDOM(nij_before_j0(J_0)*LP50*3)

      DO J=J_0,J_1
      CALL BURN_RANDOM((I_0-1)*LP50*3)
      DO I=I_0,IMAXJ(J)
        DO L=LP50,1,-1
          DO NR=1,3
            RNDSS(NR,L,I,J) = RANDU(xx)
          END DO
        END DO
C     Do not bother to save random numbers for isccp_clouds
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1)*LP50*3)
      END DO

      CALL BURN_RANDOM(nij_after_j1(J_1)*LP50*3)

C     But save the current seed in case isccp_routine is activated
      if (isccp_diags.eq.1) CALL RFINAL(seed)
      WCU500=0.
      SAVWCU=0.
      SAVWC1=0.
      SAVEN1=0.
      SAVEN2=0.
      W500P1=0.
      ENTJ=0.

      call recalc_agrid_uv ! may not be necessary - check later

c
c collect staggered velocities to be mixed into an A-grid array
c
      kmax_nonpolar = minval(kmaxj(j_0:j_1))
      call replicate_uv_to_agrid(ukm,vkm,kmax_nonpolar,
     &     ukmsp,vkmsp,ukmnp,vkmnp)

C
C**** SAVE UC AND VC, AND ZERO OUT CLDSS AND CLDMC
      TLS=T
      QLS=Q
      TMC=T
      QMC=Q
      FSS=1.
      IH=JHOUR+1
      IHM = IH+(JDATE-1)*24
#ifdef TRACERS_ON
C**** Find the ntx active tracers ntix(1->ntx)
      nx = 0
      do n=1,ntm
        if (itime.lt.itime_tr0(n)) cycle
        nx = nx+1
        ntix(nx) = n
      end do
      ntx = nx

#ifdef TRACERS_AMP
      AQsulfRATE = 0.d0
#endif

#endif

C****
C**** MAIN J LOOP
C****
       ICKERR=0
       JCKERR=0

       ! Burn random numbers for earlier latitudes here.
       ! Actual generation of random numbers is in CLOUDS2.f::ISCCP_CLOUD_TYPES
      if (isccp_diags.eq.1) then
        CALL BURN_RANDOM(nij_before_j0(J_0)*NCOL*(LM+1))
      end if
!xOMP  PARALLEL DO PRIVATE (
#ifdef TRACERS_ON
!xOMP*  NX,tmsave,tmomsv,
#endif
!xOMP*  tmp,ALPHAS,ALPHA1,ALPHA2,AT,BYDH1S,BYDH12, CC,CONV,CTP,
!xOMP*  DH1S,DH12,DTDZ,DTDZG,DTDZS,DUDZ,DUDZG,DUDZS,DVDZ,DVDZG,DVDZS,
!xOMP*  DTAU_S,DTAU_C,DEM_S,DEM_C, FQ_ISCCP, ENRGP,EPRCP,
!xOMP*  HCNDMC, I,ITYPE,IT,ITAU, IPRES,
#ifdef TRACERS_SPECIAL_Shindell
!xOMP*  Lfreeze,
#endif
#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!xOMP*  n1,n_fidx,
#endif
#endif
!xOMP*  ITROP,IERR, J,JERR,JR, K,KR, L,LERR, N,NBOX, PRCP,PFULL,PHALF,
!xOMP*  GZIL, SD_CLDIL, WMIL, TMOMIL, QMOMIL,        ! reduced arrays
!xOMP*  QG,QV, SKT,SSTAB, TGV,TPRCP,THSV,THV1,THV2,TAUOPT,TSV, WMERR,
!xOMP*  LP600,LP850,CSC,DIFT, WM1,WMI,sunlit
cECON !xOMP*  E,E1,W1,ep,ep1,q0,q1,q2,
!xOMP*    ) SCHEDULE(DYNAMIC,2)
!xOMP*    REDUCTION(+:ICKERR,JCKERR)
C
      DO J=J_0,J_1

       ! Burn random numbers for earlier longitudes here.
       ! Actual generation of random numbers is in CLOUDS2.f::ISCCP_CLOUD_TYPES
      if (isccp_diags.eq.1) then
        CALL BURN_RANDOM((I_0-1)*NCOL*(LM+1))
      end if
C
Cred* Reduced Arrays 2
C
      DO L=1,LM
      DO I=I_0,I_1
        GZIL(I,L) = GZ(I,J,L)
#ifdef SCM        
        SD_CLDIL(I,L) = SD_CLOUDS(I,J,L)
#else
        SD_CLDIL(I,L) = SDA(I,J,L)/DTsrc ! averaged SD
#endif
        WMIL(I,L) = WM(I,J,L)
        TMOMIL(:,I,L) = T3MOM(:,I,J,L)
        QMOMIL(:,I,L) = Q3MOM(:,I,J,L)
      END DO
      END DO
Cred* end Reduced Arrays 2
      kmax = kmaxj(j)
C****
C**** MAIN I LOOP
C****
      DO I=I_0,IMAXJ(J)
        DXYPIJ=AXYP(I,J)
        JR=JREG(I,J)
C****
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS FOR MSTCNV
C****
      DEBUG = .FALSE.   ! use for individual box diags in clouds
      PEARTH=FEARTH(I,J)
      PLAND=FLAND(I,J)
      TS=TSAVG(I,J)
      QS=QSAVG(I,J)
      US=USAVG(I,J)
      VS=VSAVG(I,J)
      TGV=TGVAVG(I,J)
      QG=QGAVG(I,J)
      TSV=TS*(1+QS*DELTX)
!!!   DCL=NINT(DCLEV(I,J))   ! prevented by openMP bug
      DCL=INT(DCLEV(I,J)+.5)
#ifdef SCM
      if (I.eq.I_TARG .and. J.eq.J_TARG) then
          do LL=1,LM
             do L=1,LM
                WCUALL(L,1,LL)=0.
                WCUALL(L,2,LL)=0.
                MPLUMEALL(L,1,LL)=0.
                MPLUMEALL(L,2,LL)=0.
                ENTALL(L,1,LL)=0.
                ENTALL(L,2,LL)=0.
                DETRAINDEEP(L,1,LL) = 0.0
                DETRAINDEEP(L,2,LL) = 0.0
                TPALL(L,1,LL)=0.
                TPALL(L,2,LL)=0.
                PRCCDEEP(L,1,LL) = 0.0
                PRCCDEEP(L,2,LL)  = 0.0
                NPRCCDEEP(L,1,LL) = 0.0
                NPRCCDEEP(L,2,LL) = 0.0
             enddo
          enddo
          do L=1,LM
             WCUDEEP(L,1) = 0.0
             WCUDEEP(L,2) = 0.0
             MPLUMEDEEP(L,1) = 0.0
             MPLUMEDEEP(L,2) = 0.0
             ENTDEEP(L,1) = 0.0
             ENTDEEP(L,2) = 0.0
          enddo
      endif
#endif
#ifdef CUBE_GRID
      ra = .5d0
#else
#ifdef ALT_CLDMIX_UV
      ra(1:kmax) = 1d0
#else
      DO K=1,KMAX
        RA(K)=RAVJ(K,J)
      END DO
#endif
#endif
C**** PRESSURES, AND PRESSURE TO THE KAPA
      PL(:) =PMID(:,I,J)
      PLE(:)=PEDN(:,I,J)
      PLK(:)=PK(:,I,J)
      AIRM(:)=PDSIG(:,I,J)
      BYAM(:)=1./AIRM(:)
      WTURB(:)=SQRT(.6666667*EGCM(:,I,J))
#ifdef SCM
      if (SCM_ATURB_FLAG.eq.0) then
c****     for SCM run with DRY convection - zero out WTURB
          WTURB(:) = 0.d0
      else
c****     for SCM run with ATURB
          WTURB(:)=SQRT(.6666667*EGCM(:,I,J))
      endif
#endif

C**** other fields where L is the leading index
      SVLHXL(:)=SVLHX(:,I,J)
      TTOLDL(:)=TTOLD(:,I,J)
      CLDSAVL(:)=CLDSAV(:,I,J)
      CLDSV1(:)=CLDSAV1(:,I,J)
      RH(:)=RHSAV(:,I,J)
#ifdef CLD_AER_CDNC
        OLDCDL(:)=OLDNL(:,I,J)
        OLDCDI(:)=OLDNI(:,I,J)  ! OLDNI is for rsf save
        SME(:)  =EGCM(:,I,J)  !saving 3D TKE value
        CTEML(:) =CTEM(:,I,J)
c       if(l.eq.2)write(6,*)"CTEM_DRV",CTEML(L),SME(L),OLDCDL(L)
        CD3DL(:) =CD3D(:,I,J)
        CL3DL(:) =CL3D(:,I,J)
        CI3DL(:) =CI3D(:,I,J)
        CDN3DL(:)=CDN3D(:,I,J)
        CRE3DL(:)=CRE3D(:,I,J)
        SMLWP=CLWP(I,J)
#ifdef TRACERS_AMP
C**not sure if this needed
#ifdef BLK_2MOM
       do n=1,nmodes
        NACTC(:,n)= NACTV(I,J,:,n)
       enddo
#endif
#endif
#endif
      FSSL(:)=FSS(:,I,J)
      DPDT(1:LS1-1)=SIG(1:LS1-1)*(P(I,J)-PTOLD(I,J))*BYDTsrc
      DPDT(LS1:LM)=0.
      DO L=1,LM
C**** TEMPERATURES
        SM(L)  =T(I,J,L)*AIRM(L)
        SMOM(:,L) =TMOMIL(:,I,L)*AIRM(L)
        SMOMMC(:,L) =SMOM(:,L)
        SMOMLS(:,L) =SMOM(:,L)
        TL(L)=T(I,J,L)*PLK(L)
C**** MOISTURE (SPECIFIC HUMIDITY)
        QM(L)  =Q(I,J,L)*AIRM(L)
        QMOM(:,L) =QMOMIL(:,I,L)*AIRM(L)
        QMOMMC(:,L) =QMOM(:,L)
        QMOMLS(:,L) =QMOM(:,L)
        WML(L)=WMIL(I,L)
        QL(L) =Q(I,J,L)
C**** others
        SDL(L)=SD_CLDIL(I,L)*BYAXYP(I,J)
        TVL(L)=TL(L)*(1.+DELTX*QL(L))
        W2L(L)=W2GCM(L,I,J)
        SAVWL(L)=0.
        SAVWL1(L)=0.
        SAVE1L(L)=0.
        SAVE2L(L)=0.
        IF(L.LE.LM-2)
     *    ETAL(L+1)=.5*ENTCON*(GZIL(I,L+2)-GZIL(I,L))*1.d-3*BYGRAV
        IF(L.LE.LM-2) GZL(L+1)=ETAL(L+1)/ENTCON
      END DO


      ETAL(LM)=ETAL(LM-1)
      ETAL(1)=0.     ! not used
      GZL(LM)=GZL(LM-1)
      GZL(1)=0.
#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
      do l=1,lm
        tm(l,nx) = trm(i,j,l,ntix(nx))
        tmom(:,l,nx) = trmom(:,i,j,l,ntix(nx))
      end do
      end do
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2M
C** Add activated fraction, NACTV
      do nx=1,nmodes
      do l=1,lm
        NACTC(l,nx)=NACTV(i,j,l,nx)
      end do
      end do
#endif
#endif

C**** SURROUNDING WINDS

      if(j.eq.1 .and. have_south_pole) then
        U_0(1:KMAX,:) = UKMSP(1:KMAX,:)
        V_0(1:KMAX,:) = VKMSP(1:KMAX,:)
      elseif(j.eq.jm .and. have_north_pole) then
        U_0(1:KMAX,:) = UKMNP(1:KMAX,:)
        V_0(1:KMAX,:) = VKMNP(1:KMAX,:)
      else
        U_0(1:KMAX,:) = UKM(1:KMAX,:,I,J)
        V_0(1:KMAX,:) = VKM(1:KMAX,:,I,J)
      endif
      DO L=1,LM
        DO K=1,KMAX
          UM(K,L) = U_0(K,L)*AIRM(L)
          VM(K,L) = V_0(K,L)*AIRM(L)
          UM1(K,L) = UM(K,L)
          VM1(K,L) = VM(K,L)
        END DO
      END DO

C**** INITIALISE PRECIPITATION AND LATENT HEAT
      PRCP=0.
      ENRGP=0.
C**** temperature of precip is based on pre-mstcnv profile
      TPRCP=T(I,J,1)*PLK(1)-TF
#ifdef TRACERS_WATER
      TRPREC(:,I,J) = 0.
#endif

C**** SET DEFAULTS FOR AIR MASS FLUX (STRAT MODEL)
      AIRX(I,J)=0.
      DDM1(I,J)=0.
      DDMS(I,J)=0.
      DDML(I,J)=0.
      TDN1(I,J)=0.
      QDN1(I,J)=0.
#ifdef TRACERS_ON
      TRDN1(:,I,J)=0.
#endif
C****
C**** Energy conservation note: For future reference the energy function
C**** for these column calculations (assuming energy reference level
C**** of 0 K for air, and 0 C liquid for water) is:
C****  E = SH + LH_vapour + LH_clw + ENRGP
C****    =  (sum(TL(:)*AIRM(:))*SHA + sum(QM(:))*LHE +
C****        sum(WML(:)*(LHE-SVLHXL(:))*AIRM(:)))*100.*BYGRAV
C**** The LH_clw term is slightly different after MSTCNV:
C****   LH_clw = sum((WML(:)*(LHE-SVLHXL(:))+SVWMXL(:)*(LHE-SVLATL(:)))
C****                *AIRM(:))*100.*BYGRAV
C**** After LSCOND, latent heat changes to:
C****          = sum(WMX(:)*(LHE-SVLHXL(:))*AIRM(:))*100.*BYGRAV
C****
C**** Note that the column changes after MSTCNV only apply to the
C**** moist convective fraction (1-FSSL(:)), and after LSCOND, FSSL(:).
C**** Condensate is always defined over the whole box.
C****
c**** uncomment lines marked ECON to check energy conservation
c**** uncomment lines marked QCON to check water conservation

cQCON q0 = sum(QM(:)+WML(:)*AIRM(:))*100.*BYGRAV
cECON  E = (sum(TL(:)*AIRM(:))*SHA + sum(QM(:))*LHE +sum(WML(:)*(LHE
cECON*     -SVLHXL(:))*AIRM(:)))*100.*BYGRAV

C**** MOIST CONVECTION
      CALL MSTCNV(IERR,LERR,i,j)

cECON E1 = ( sum( ((T(I,J,:)*PLK(:)-TL(:))*AIRM(:)*SHA + (Q(I,J,:)
cECON*     *AIRM(:)-QM(:))*LHE)*(1.-FSSL(:)))-sum(SVWMXL(:)*(LHE
cECON*     -SVLATL(:))*AIRM(:)))*100.*BYGRAV 
cQCON q1 = sum( (Q(I,J,:)*AIRM(:)-QM(:))*(1.-FSSL(:))-SVWMXL(:)*AIRM(:))
cQCON*     *100.*BYGRAV

C**** Error reports
      if (ierr.gt.0) then
        write(6,*) "Error in moist conv: i,j,l=",i,j,lerr
        if (ierr.eq.2) ickerr = ickerr + 1
      end if

#ifdef TRACERS_SPECIAL_Shindell
C**** Calculate NOx from lightning:
C**** first, need the local freezing level:
      IF(LMCMAX.gt.0)THEN
        Lfreeze=1
        DO L=1,LMCMAX
          IF(T(I,J,L)*PLK(L).lt.TF) THEN
            Lfreeze=L
            EXIT
          END IF
        END DO
        CALL calc_lightning(I,J,LMCMAX,Lfreeze)
      END IF
#endif

C**** ACCUMULATE MOIST CONVECTION DIAGNOSTICS
      IF (LMCMIN.GT.0) THEN
        AIJ(I,J,IJ_PSCLD)=AIJ(I,J,IJ_PSCLD)+CLDSLWIJ
        AIJ(I,J,IJ_PDCLD)=AIJ(I,J,IJ_PDCLD)+CLDDEPIJ
        IF(CLDSLWIJ.GT.1e-6) AIJ(I,J,IJ_SCNVFRQ)=AIJ(I,J,IJ_SCNVFRQ)+1.
        IF(CLDDEPIJ.GT.1e-6) AIJ(I,J,IJ_DCNVFRQ)=AIJ(I,J,IJ_DCNVFRQ)+1.
        AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
#ifdef CLD_AER_CDNC
        AIJ(I,J,IJ_WMCLWP)=AIJ(I,J,IJ_WMCLWP)+WMCLWP
        AIJ(I,J,IJ_WMCTWP)=AIJ(I,J,IJ_WMCTWP)+WMCTWP
#endif

        HCNDMC=0.
        DO L=1,LMCMAX
          HCNDMC=HCNDMC+DGDSM(L)+DPHASE(L)
          call inc_ajl(i,j,l,jl_mchr,DGDSM(L)*BYDSIG(L))
          call inc_ajl(i,j,l,jl_mchphas,DPHASE(L)*BYDSIG(L))
          call inc_ajl(i,j,l,jl_mcdtotw,DTOTW(L)*BYDSIG(L))
CCC       IF(J.GE.J5S.AND.J.LE.J5N) AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+
CCC  *         (DGDSM(L)+DPHASE(L))*(AXYP(I,J)*BYDSIG(L))
          AIJL(I,J,L,IJL_MC) = AIJL(I,J,L,IJL_MC) + (DPHASE(L)+DGDSM(L))
          call inc_ajl(i,j,l,jl_mcheat,(DPHASE(L)+DGDSM(L)))
          call inc_ajl(i,j,l,jl_mcdry,(DQCOND(L)-DGDQM(L)))
          call inc_ajl(i,j,l,jl_mcshlw,(DPHASHLW(L)+DGSHLW(L)))
          call inc_ajl(i,j,l,jl_mcdeep,(DPHADEEP(L)+DGDEEP(L)))
C*** Begin Accumulate 3D convective latent heating
         if(lh_diags.eq.1) then
          AIJL(I,J,L,IJL_MCTLH)=AIJL(I,J,L,IJL_MCTLH)+
     &         (DPHASE(L)+DGDSM(L))
          AIJL(I,J,L,IJL_MCDLH)=AIJL(I,J,L,IJL_MCDLH)+
     &         (DPHADEEP(L)+DGDEEP(L))
          AIJL(I,J,L,IJL_MCSLH)=AIJL(I,J,L,IJL_MCSLH)+
     &         (DPHASHLW(L)+DGSHLW(L))
         endif
C*** End Accumulate 3D convective latent heating
         call inc_ajl(i,j,l,jl_mcmflx,MCFLX(L))
         call inc_ajl(i,j,l,jl_cldmc,CLDMCL(L)*AIRM(L))
         call inc_ajl(i,j,l,jl_mcdflx,DDMFLX(L))
         call inc_ajl(i,j,l,jl_csizmc,CSIZEL(L)*CLDMCL(L)*AIRM(L))
#ifdef HTAP_LIKE_DIAGS
          AIJ(I,J,IJ_MCamFX(L))=AIJ(I,J,IJ_MCamFX(L))+MCFLX(L)
#endif
        END DO
        DO IT=1,NTYPE
          CALL INC_AJ(I,J,IT,J_PRCPMC,PRCPMC*FTYPE(IT,I,J))
        END DO
        CALL INC_AREG(I,J,JR,J_PRCPMC,PRCPMC)
        DO KR=1,NDIUPT
        IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
          tmp(IDD_PR)  =+PRCPMC
          tmp(IDD_ECND)=+HCNDMC
          tmp(IDD_MCP) =+PRCPMC
          tmp(IDD_DMC) =+CLDDEPIJ
          tmp(IDD_SMC) =+CLDSLWIJ
          ADIURN(IDX1(:),KR,IH)=ADIURN(IDX1(:),KR,IH)+TMP(IDX1(:))
#ifndef NO_HDIURN
          HDIURN(IDX1(:),KR,IHM)=HDIURN(IDX1(:),KR,IHM)+TMP(IDX1(:))
#endif
        END IF
        END DO
#ifdef CLD_AER_CDNC
        DO L =1,LM
        IF (NMCW.ge.1) then
          AIJ(I,J,IJ_3dNWM)=AIJ(I,J,IJ_3dNWM)+ACDNWM(L)
          AIJ(I,J,IJ_3dRWM)=AIJ(I,J,IJ_3dRWM)+AREWM(L)
          AIJ(I,J,IJ_3dLWM)=AIJ(I,J,IJ_3dLWM)+ALWWM(L)
          AIJL(I,J,L,IJL_REWM)= AIJL(I,J,L,IJL_REWM)+AREWM(L)
          AIJL(I,J,L,IJL_CDWM)= AIJL(I,J,L,IJL_CDWM)+ACDNWM(L)
          AIJL(I,J,L,IJL_CWWM)= AIJL(I,J,L,IJL_CWWM)+ALWWM(L)
          call inc_ajl(i,j,l,JL_CNUMWM,ACDNWM(L)*AIRM(L))
c         write(6,*)"IJL_REWM",AIJL(I,J,L,IJL_REWM),I,J,L,
c     *   AIJL(I,J,L,IJL_CDWM),AIJL(I,J,L,IJL_CWWM),ALWWM(L)
        ENDIF

        IF (NMCI.ge.1) then
          AIJ(I,J,IJ_3dNIM)=AIJ(I,J,IJ_3dNIM)+ACDNIM(L)
          AIJ(I,J,IJ_3dRIM)=AIJ(I,J,IJ_3dRIM)+AREIM(L)
          AIJ(I,J,IJ_3dLIM)=AIJ(I,J,IJ_3dLIM)+ALWIM(L)
          call inc_ajl(i,j,l,JL_CNUMIM,ACDNIM(L)*AIRM(L))
          AIJL(I,J,L,IJL_REIM)= AIJL(I,J,L,IJL_REIM)+AREIM(L)
          AIJL(I,J,L,IJL_CDIM)= AIJL(I,J,L,IJL_CDIM)+ACDNIM(L)
          AIJL(I,J,L,IJL_CWIM)= AIJL(I,J,L,IJL_CWIM)+ALWIM(L)
c         write(6,*)"IJL_REIM",AIJL(I,J,L,IJL_REIM),I,J,L,
c     *   AIJL(I,J,L,IJL_CDIM),AIJL(I,J,L,IJL_CWIM),ALWIM(L)
        ENDIF

        ENDDO
#endif
C**** ACCUMULATE PRECIP
        PRCP=PRCPMC*100.*BYGRAV
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        precnvl(1)=precnvl(1)+prcpmc*bygrav
#endif
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
C**** NEED TO TAKE ACCOUNT OF LATENT HEAT THOUGH
        IF (TPRCP.gt.0) THEN
C         EPRCP=PRCP*TPRCP*SHW
          EPRCP=0.
          ENRGP=ENRGP+EPRCP
cECON     ep=0.
        ELSE
C         EPRCP=PRCP*TPRCP*SHI
          EPRCP=0.
          ENRGP=ENRGP+EPRCP-PRCP*LHM
cECON     ep=-PRCP*LHM
          AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCP
        END IF
        AIJ(I,J,IJ_PRECMC)=AIJ(I,J,IJ_PRECMC)+PRCP

C**** Uncomment next lines for check on conservation
cECON   if (abs(E1-ep).gt.0.01) print*,"energy err0",i,j,(E1-ep)
cECON*       *GRAV/100.,E,E1,ep,prcp,tprcp
cQCON   if (abs(q1-prcp).gt.0.01) print*,"water err0",i,j,(q1-prcp)
cQCON*       *GRAV/100.,q0,q1,prcp

        DO L=1,LMCMAX
          T(I,J,L)=(1.-FSSL(L))*SM(L)*BYAM(L)+FSSL(L)*TLS(I,J,L)
          Q(I,J,L)=(1.-FSSL(L))*QM(L)*BYAM(L)+FSSL(L)*QLS(I,J,L)
          TMC(I,J,L)=SM(L)*BYAM(L)
          QMC(I,J,L)=QM(L)*BYAM(L)
          SMOMMC(:,L)=SMOM(:,L)
          QMOMMC(:,L)=QMOM(:,L)
          DO K=1,KMAX
            UM1(K,L)=UM(K,L)
            VM1(K,L)=VM(K,L)
          END DO
        END DO

        CSIZMC(1:LMCMAX,I,J)=CSIZEL(1:LMCMAX)
        FSS(:,I,J)=FSSL(:)
        AIRX(I,J) = AIRXL*AXYP(I,J)
        DO L=1,DCL
          DDML(I,J)=L                    ! the lowest downdraft layer
          IF(DDMFLX(L).GT.0.d0) EXIT
        END DO
        IF (DDML(I,J).gt.0) THEN
          TDN1(I,J)=TDNL(DDML(I,J))     ! downdraft temperature
          QDN1(I,J)=QDNL(DDML(I,J))     ! downdraft humidity
          DDMS(I,J)=-100.*DDMFLX(DDML(I,J))/(GRAV*DTsrc) ! downdraft mass flux
#ifdef TRACERS_ON
          TRDN1(:,I,J)=1d-2*TRDNL(:,DDML(I,J))*GRAV*BYAXYP(I,J) ! downdraft tracer conc
#endif
        END IF

C**** level 1 downfdraft mass flux/rho (m/s)
        DDM1(I,J) = (1.-FSSL(1))*DDMFLX(1)*RGAS*TSV/(GRAV*PEDN(1,I,J)
     *       *DTSrc)
      END IF
#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
        n = ntix(nx)
        do l=1,lm
#ifndef SKIP_TRACER_DIAGS
           if (itcon_mc(n).gt.0) call inc_diagtcb(i,j,
     *          (tm(l,nx)-trm(i,j,l,n))*(1.-fssl(l))
#ifdef TRACERS_WATER
     *         + trsvwml(nx,l)
#endif
     *          ,itcon_mc(n),n)
          call inc_tajln(i,j,l,jlnt_mc,n,
     &         (tm(l,nx)-trm(i,j,l,n))*(1.-fssl(l))
#ifdef TRACERS_WATER
     *         + trsvwml(nx,l)
#endif
     *         )
#endif /*SKIP_TRACER_DIAGS*/
#ifdef TRACERS_WATER
          trwml(nx,l) = trwm(i,j,l,n)+trsvwml(nx,l)
#endif
          tmsave(l,nx) = tm(l,nx) ! save for tajln(large-scale condense)
          tmomsv(:,l,nx) = tmom(:,l,nx)
          tm(l,nx) = trm(i,j,l,n)*fssl(l)   ! kg in lsc fraction only
          tmom(:,l,nx) = trmom(:,i,j,l,n)*fssl(l)
        end do
#ifdef TRACERS_WATER
        trprec(n,i,j) = trprmc(nx)
#endif
      end do
#endif
      LMC(1,I,J) = LMCMIN
      LMC(2,I,J) = LMCMAX+1
C****
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS FOR LSCOND
C****
      DO L=1,LM
        TL(L)=TLS(I,J,L)*PLK(L)
        TH(L)=TLS(I,J,L)
        QL(L)=QLS(I,J,L)
        SMOM(:,L)=SMOMLS(:,L)
        QMOM(:,L)=QMOMLS(:,L)
      END DO
      WMX(:)=WML(:)+SVWMXL(:)
      AQ(:)=(QL(:)-QTOLD(:,I,J))*BYDTsrc
#ifdef SCM
      if (I.eq.I_TARG .and. J.eq.J_TARG) AQ(:) = ((SCM_SAVE_Q(:)
     *     +SCM_DEL_Q(:))-QTOLD(:,I,J))*BYDTsrc 
#endif
      RNDSSL(:,1:LP50)=RNDSS(:,1:LP50,I,J)
      FSSL(:)=FSS(:,I,J)
      DO L=1,LM
        DO K=1,KMAX
          UM(K,L) = U_0(K,L)*AIRM(L)
          VM(K,L) = V_0(K,L)*AIRM(L)
        END DO
      END DO
C****
C**** COMPUTE STRATOCUMULUS CLOUDS USING PHILANDER'S FORMULA
C****
      IF (ISC.EQ.1.AND.FOCEAN(I,J).GT..5) THEN
        CSC=0.D0
        LP600=LM
        LP850=LM
        DO L=2,LM
          IF(L.GT.LP600) EXIT
          IF(PL(L).LT.600.) THEN
            LP600=L
            IF(600.-PL(L).GT.PL(L-1)-600.) LP600=L-1
          ENDIF
        ENDDO
        DO L=2,LM
          IF(L.GT.LP850) EXIT
          IF(PL(L).LT.850.) THEN
            LP850=L
            IF(850.-PL(L).GT.PL(L-1)-850.) LP850=L-1
          ENDIF
        ENDDO
        IF(SDL(LP600)+SDL(LP600+1).GT.0.) THEN
          DIFT=TL(LP850)-TGV/(1.+DELTX*QG)
          CSC=.031D0*DIFT+.623D0
          IF(CSC.LT.0.) CSC=0.
        ENDIF
        CLDMCL(1)=CLDMCL(1)+CSC
        IF(CSC.GT.0.) TAUMCL(1)=AIRM(1)*.08D0
        IF(CLDMCL(1).GT.1.) CLDMCL(1)=1.
C     IF(CSC.GT.0.) WRITE (6,*) I,J,DCL,TL(LP850),TGV/(1.+DELTX*QG),CSC
      ENDIF

C**** COMPUTE RICHARDSON NUMBER FROM SURFACE CONDITIONS WHEN DEPTH OF
C**** BOUNDARY LAYER IS AT OR BELOW FIRST LAYER (E.G. AT NIGHT)
      IF(DCL.LE.1) THEN
        THSV=TS*(1.+DELTX*QS)/PEK(1,I,J)
        THV1=TH(1)*(1.+DELTX*QL(1))
        THV2=TH(2)*(1.+DELTX*QL(2))
        ALPHAS=2./((TGV/(1.+DELTX*QG)+TS)/PEK(1,I,J))
        ALPHA1=2./(TH(1)+TS/PEK(1,I,J))
        ALPHA2=2./(TH(1)+TH(2))
        DH1S=(PLE(1)-PL(1))*TL(1)*RGAS/(GRAV*PL(1))
        BYDH1S=1./DH1S
        DH12=(GZIL(I,2)-GZIL(I,1))*BYGRAV
        BYDH12=1./DH12
        DTDZS=(THV1-THSV)*BYDH1S
        DTDZ=(THV2-THV1)*BYDH12
        DUDZ=(UA(2,I,J)-UA(1,I,J))*BYDH12
        DVDZ=(VA(2,I,J)-VA(1,I,J))*BYDH12
        DUDZS=(UA(1,I,J)-US)*BYDH1S
        DVDZS=(VA(1,I,J)-VS)*BYDH1S
        DUDZG=.1d0*US
        DVDZG=.1d0*VS
        DTDZG=.1d0*(THSV-TGV/PEK(1,I,J))
        RIS=(GRAV*ALPHAS*DTDZG)/(DUDZG*DUDZG+DVDZG*DVDZG+teeny)
        RI1=(GRAV*ALPHA1*DTDZS)/(DUDZS*DUDZS+DVDZS*DVDZS+teeny)
        RI2=(GRAV*ALPHA2*DTDZ)/(DUDZ*DUDZ+DVDZ*DVDZ+teeny)
C       WRITE (6,*)'I,J,QG,TGV,THSV,RIS,RI1=',I,J,QG,TGV,THSV,RIS,RI1
      ENDIF

C**** Uncomment next few lines for check on conservation
cECON W1 = sum( (WML(1:LP50)*(LHE-SVLHXL(1:LP50))
cECON*     +SVWMXL(1:LP50)*(LHE-SVLATL(1:LP50)))*AIRM(1:LP50))
cECON E = ( sum(TL(1:LP50)*AIRM(1:LP50))*SHA + sum(QL(1:LP50)
cECON*     *AIRM(1:LP50))*LHE )*100.*BYGRAV

C**** LARGE-SCALE CLOUDS AND PRECIPITATION
      CALL LSCOND(IERR,WMERR,LERR,i,j)

cECON E1 = ( sum( ((TLS(I,J,1:LP50)-TH(1:LP50))*PLK(1:LP50)*AIRM(1:LP50)
cECON*     *SHA +(QLS(I,J,1:LP50)-QL(1:LP50))*AIRM(1:LP50)*LHE)
cECON*     *FSSL(1:LP50))+W1-sum(WMX(1:LP50)*(LHE-SVLHXL(1:LP50))
cECON*     *AIRM(1:LP50)) )*100.*BYGRAV

C**** Error reports
      IF (IERR.ne.0) WRITE(99,'(I10,3I4,A,D14.5,A)')
     *       Itime,I,J,LERR,' CONDSE:H2O<0',WMERR,' ->0'

C**** Accumulate diagnostics of LSCOND
         AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
         DO IT=1,NTYPE
           CALL INC_AJ(I,J,IT,J_PRCPSS,PRCPSS*FTYPE(IT,I,J))
         END DO
         CALL INC_AREG(I,J,JR,J_PRCPSS,PRCPSS)
         DO KR=1,NDIUPT
         IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
           tmp(IDD_PR)  =+PRCPSS
           tmp(IDD_ECND)=+HCNDSS
           tmp(IDD_SSP) =+PRCPSS
           ADIURN(IDX2(:),KR,IH)=ADIURN(IDX2(:),KR,IH)+TMP(IDX2(:))
#ifndef NO_HDIURN
           HDIURN(IDX2(:),KR,IHM)=HDIURN(IDX2(:),KR,IHM)+TMP(IDX2(:))
#endif
         END IF
         END DO

C**** TOTAL PRECIPITATION AND AGE OF SNOW
      PRCP=PRCP+PRCPSS*100.*BYGRAV

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      DO l=1,lm
         prelay(i,j,l)=((prebar1(l)*DTsrc*100.+precnvl(l)*100.)+
     &   (prebar1(l+1)*DTsrc*100.+precnvl(l+1)*100.))/2.
      END DO
#endif

C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
C**** NEED TO TAKE ACCOUNT OF LATENT HEAT THOUGH
      IF (LHP(1).ne.LHS) THEN
C       EPRCP=PRCPSS*100.*BYGRAV*TPRCP*SHW
        EPRCP=0.
        ENRGP=ENRGP+EPRCP
cECON   ep1=0.
      ELSE
C       EPRCP=PRCPSS*100.*BYGRAV*TPRCP*SHI
        EPRCP=0.
        ENRGP=ENRGP+EPRCP-PRCPSS*100.*BYGRAV*LHM
cECON   ep1=-PRCPSS*100.*BYGRAV*LHM
        AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCPSS*100.*BYGRAV
      END IF

cECON if (abs(E1-ep1).gt.0.01) print*,"energy err1",i,j,(E1-ep1)
cECON*     *GRAV*1d-2,E1-ep1,E1,ep1,prcpss*100.*BYGRAV,lhp(1)
      
C**** PRECIPITATION DIAGNOSTICS
        DO IT=1,NTYPE
          CALL INC_AJ(I,J,IT,J_EPRCP,ENRGP*FTYPE(IT,I,J))
        END DO
        CALL INC_AREG(I,J,JR,J_EPRCP,ENRGP)
        AIJ(I,J,IJ_PREC)=AIJ(I,J,IJ_PREC)+PRCP
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+ENRGP
        AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+
     *       ENRGP*FOCEAN(I,J)*(1.-RSI(I,J))
        AIJ(I,J,IJ_FWOC)=AIJ(I,J,IJ_FWOC)+
     *       PRCP*FOCEAN(I,J)*(1.-RSI(I,J))

      IF(ENRGP.LT.0.) THEN    ! MODIFY SNOW AGES AFTER SNOW FALL
        DO ITYPE=1,3
          SNOAGE(ITYPE,I,J)=SNOAGE(ITYPE,I,J)*EXP(-PRCP)
        END DO
      END IF

C**** cloud water diagnostics
      WM1=0  ; WMI=0
      DO L=1,LP50
        WM1=WM1+WMX(L)*AIRM(L)
        IF (SVLHXL(L).eq.LHS) WMI=WMI+WMX(L)*AIRM(L)
      END DO
      AIJ(I,J,IJ_CLDW)=AIJ(I,J,IJ_CLDW)+WM1*100.*BYGRAV   ! all condensate
      AIJ(I,J,IJ_CLDI)=AIJ(I,J,IJ_CLDI)+WMI*100.*BYGRAV   ! ice only

C**** Calculate ISCCP cloud diagnostics if required
      if (isccp_diags.eq.1) then
        do l=1,lm
          cc(l)=cldmcl(LM+1-L)+cldssl(LM+1-L)
          if(cc(l) .gt. 1.) then
            cc(l)=1.
          endif
          conv(l)=cldmcl(LM+1-L)
          if(conv(l) .gt. 1.) then
            conv(l)=1.
          endif

          dtau_s(l)=taussl(LM+1-L)
          dtau_c(l)=taumcl(LM+1-L)
          pfull(l)=pl(LM+1-L)*100.
          phalf(l)=ple(LM+2-L)*100.
          at(l)=tl(LM+1-L)  ! in situ temperature

C**** set skt from radiative temperature
          skt=sqrt(sqrt(
     *         (focean(i,j)+flake(i,j))*(1.-rsi(i,j))*gtempr(1,i,j)**4+
     *         (focean(i,j)+flake(i,j))*    rsi(i,j) *gtempr(2,i,j)**4+
     *         flice(i,j) *gtempr(3,i,j)**4+
     *         fearth(i,j)*gtempr(4,i,j)**4))
          dem_s(l)=0.
          dem_c(l)=0.
          if(svlhxl(LM+1-L) .eq. lhe )   ! large-scale water cloud
     *      dem_s(l)=1.-exp(-taussl(LM+1-L)*bywc)
          if(svlatl(LM+1-L) .eq. lhe )   ! convective water cloud
     *      dem_c(l)=1.-exp(-taumcl(LM+1-L)*bywc)
          if(svlhxl(LM+1-L) .eq. lhs )   ! large-scale ice cloud
     *      dem_s(l)=1.-exp(-taussl(LM+1-L)*byic)
          if(svlatl(LM+1-L) .eq. lhs )   ! convective ice cloud
     *      dem_c(l)=1.-exp(-taumcl(LM+1-L)*byic)

          qv(l)=ql(LM+1-L)
        end do
        phalf(lm+1)=ple(1)*100.
        itrop = LM+1-LTROPO(I,J)
        sunlit=0
        if (cosz1(i,j).gt.0) sunlit=1

        call ISCCP_CLOUD_TYPES(sunlit,pfull,phalf,qv,
     &       cc,conv,dtau_s,dtau_c,skt,
     &       at,dem_s,dem_c,itrop,fq_isccp,ctp,tauopt,nbox,jerr)
        if(jerr.ne.0) jckerr = jckerr + 1

C**** set ISCCP diagnostics
        AIJ(I,J,IJ_SCLDI) = AIJ(I,J,IJ_SCLDI) + sunlit
        if (nbox.gt.0.and.sunlit.gt.0) then
          AIJ(I,J,IJ_CTPI) = AIJ(I,J,IJ_CTPI) + ctp
          AIJ(I,J,IJ_TAUI) = AIJ(I,J,IJ_TAUI) + tauopt
          AIJ(I,J,IJ_TCLDI)= AIJ(I,J,IJ_TCLDI)+ 1.
C**** note LOW CLOUDS:       ipres=6,7
C****      MID-LEVEL CLOUDS: ipres=4,5
C****      HIGH CLOUDS:      ipres=1,2,3
C**** Sum over itau=2,ntau (itau=1 is no cloud)
          AIJ(I,J,IJ_LCLDI)=AIJ(I,J,IJ_LCLDI)+sum(fq_isccp(2:ntau,6:7))
          AIJ(I,J,IJ_MCLDI)=AIJ(I,J,IJ_MCLDI)+sum(fq_isccp(2:ntau,4:5))
          AIJ(I,J,IJ_HCLDI)=AIJ(I,J,IJ_HCLDI)+sum(fq_isccp(2:ntau,1:3))
C**** Save area weighted isccp histograms
          n=isccp_reg2d(i,j)
          if (n.gt.0) AISCCP(:,:,n) = AISCCP(:,:,n) +
     &         fq_isccp(:,:)*axyp(i,j)
        end if
      end if

C**** Peak static stability diagnostic
      SSTAB=-1.d30
      DO L=1,DCL
Cred    IF(SSTAB.lt.(TH(L+1)-TH(L))/(GZ(I,J,L+1)-GZ(I,J,L)))
Cred *     SSTAB =  (TH(L+1)-TH(L))/(GZ(I,J,L+1)-GZ(I,J,L))
        IF(SSTAB.lt.(TH(L+1)-TH(L))/(GZIL(I,L+1)-GZIL(I,L)))
     *     SSTAB =  (TH(L+1)-TH(L))/(GZIL(I,L+1)-GZIL(I,L))
      END DO
      AIJ(I,J,ij_sstabx) = AIJ(I,J,ij_sstabx) + SSTAB

C**** WRITE TO GLOBAL ARRAYS
      TAUMC(:,I,J)=TAUMCL(:)
      CLDMC(:,I,J)=CLDMCL(:)
      SVLAT(:,I,J)=SVLATL(:)

      TAUSS(:,I,J)=TAUSSL(:)
      CLDSS(:,I,J)=CLDSSL(:)
      CLDSAV(:,I,J)=CLDSAVL(:)
      CLDSAV1(:,I,J)=CLDSV1(:)
      SVLHX(:,I,J)=SVLHXL(:)
      CSIZSS(:,I,J)=CSIZEL(:)

      RHSAV(:,I,J)=RH(:)
#ifdef CLD_AER_CDNC
         OLDNL(:,I,J)=OLDCDL(:)
         OLDNI(:,I,J)=OLDCDI(:)
         EGCM(:,I,J) =SME(:)
         CTEM(:,I,J) =CTEML(:)
         CD3D(:,I,J) =CD3DL(:)
         CL3D(:,I,J) =CL3DL(:)
         CI3D(:,I,J) =CI3DL(:)
         CDN3D(:,I,J)=CDN3DL(:)
         CRE3D(:,I,J)=CRE3DL(:)
         CLWP(I,J) = SMLWP
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
       do n=1,nmodes
        NACTV(I,J,:,n) =  NACTC(:,n)
       enddo
#endif
#endif
#endif
      TTOLD(:,I,J)=TH(:)
      QTOLD(:,I,J)=QL(:)
 
      PREC(I,J)=PRCP            ! total precip mass (kg/m^2)
      EPREC(I,J)=ENRGP          ! energy of precipitation (J/m^2)
C**** The PRECSS array is only used if a distinction is being made
C**** between kinds of rain in the ground hydrology.
      PRECSS(I,J)=PRCPSS*100.*BYGRAV  ! large scale precip (kg/m^2)
C**** accumulate precip specially for SUBDD
      P_acc(I,J)=P_acc(I,J)+PRCP
      PM_acc(I,J)=PM_acc(I,J)+PRCP-PRECSS(I,J)
#ifdef SCM
c**** save total precip for time step (in mm/hr) for SCM
      if (I.eq.I_TARG .and. J.eq.J_TARG) then
          PRCSS = PRECSS(I,J)*(3600./DTsrc)
          PRCMC = (PREC(I,J)-PRECSS(I,J))*(3600./DTsrc)
      endif
#endif

#ifdef INTERACTIVE_WETLANDS_CH4
C**** update running-average of precipitation (in mm/day):
      call running_average(prcp*sday*byDTsrc,I,J,1.d0,n__prec)
#endif

      DO L=1,LM
        call inc_ajl(i,j,l,JL_SSHR,SSHR(L))
C*** Begin Accumulate 3D heating by large scale condensation --
       if(lh_diags.eq.1) then
        AIJL(I,J,L,IJL_LLH)=AIJL(I,J,L,IJL_LLH)+SSHR(L)
       endif
C*** End Accumulate 3D heating by large scale condensation --
       call inc_ajl(i,j,l,JL_MCLDHT,DCTEI(L))
       call inc_ajl(i,j,l,JL_RHE,RH1(L))
       call inc_ajl(i,j,l,JL_CLDSS,CLDSSL(L)*AIRM(L))
       call inc_ajl(i,j,l,JL_CSIZSS,CSIZEL(L)*CLDSSL(L)*AIRM(L))
c       write(6,*) "CTEM_DRV",CTEML(L),CTEM(I,J,L),L,I,J

        T(I,J,L)=TH(L)*FSSL(L)+TMC(I,J,L)*(1.-FSSL(L))
        Q(I,J,L)=QL(L)*FSSL(L)+QMC(I,J,L)*(1.-FSSL(L))
        SMOM(:,L)=SMOM(:,L)*FSSL(L)+SMOMMC(:,L)*(1.-FSSL(L))
        QMOM(:,L)=QMOM(:,L)*FSSL(L)+QMOMMC(:,L)*(1.-FSSL(L))
C**** update moment changes
        TMOMIL(:,I,L)=SMOM(:,L)*BYAM(L)
        QMOMIL(:,I,L)=QMOM(:,L)*BYAM(L)
        WMIL(I,L)=WMX(L)

C**** CALCULATE WIND TENDENCIES AND STORE IN UKM,VKM
         IF(J.EQ.1 .AND. HAVE_SOUTH_POLE)  THEN
            DO K=1,IM ! KMAX
              UKMSP(K,L)=(UM(K,L)*FSSL(L)+UM1(K,L)*(1.-FSSL(L)))*BYAM(L)
     *             -UKMSP(K,L)
              VKMSP(K,L)=(VM(K,L)*FSSL(L)+VM1(K,L)*(1.-FSSL(L)))*BYAM(L)
     *             -VKMSP(K,L)
            END DO
         ELSE IF(J.EQ.JM .AND. HAVE_NORTH_POLE)  THEN
            DO K=1,IM ! KMAX
              UKMNP(K,L)=(UM(K,L)*FSSL(L)+UM1(K,L)*(1.-FSSL(L)))*BYAM(L)
     *             -UKMNP(K,L)
              VKMNP(K,L)=(VM(K,L)*FSSL(L)+VM1(K,L)*(1.-FSSL(L)))*BYAM(L)
     *             -VKMNP(K,L)
            END DO
         ELSE
            DO K=1,KMAX
            UKM(K,L,I,J)=(UM(K,L)*FSSL(L)+UM1(K,L)*(1.-FSSL(L)))*BYAM(L)
     *             -UKM(K,L,I,J)
            VKM(K,L,I,J)=(VM(K,L)*FSSL(L)+VM1(K,L)*(1.-FSSL(L)))*BYAM(L)
     *             -VKM(K,L,I,J)
            END DO
         END IF
      ENDDO

C**** Uncomment next two lines for check on water conservation
cQCON q2=sum((Q(I,J,:)+WMX(:))*AIRM(:))*100*BYGRAV+PRCP
cQCON if (abs(q2-q0).gt.1d-13) print*,"water err1",i,j,q2-q0,q2,q0,q1
cQCON*     ,prcp

#ifdef CLD_AER_CDNC
        DO L=1,LM
         IF (NLSW.ge.1) then
          AIJ(I,J,IJ_3dNWS)=AIJ(I,J,IJ_3dNWS)+ACDNWS(L)
          AIJ(I,J,IJ_3dRWS)=AIJ(I,J,IJ_3dRWS)+AREWS(L)
          AIJ(I,J,IJ_3dLWS)=AIJ(I,J,IJ_3dLWS)+ALWWS(L)
          AIJL(I,J,L,IJL_REWS)= AIJL(I,J,L,IJL_REWS)+AREWS(L)
          AIJL(I,J,L,IJL_CDWS)= AIJL(I,J,L,IJL_CDWS)+ACDNWS(L)
          AIJL(I,J,L,IJL_CWWS)= AIJL(I,J,L,IJL_CWWS)+ALWWS(L)
          call inc_ajl(i,j,l,JL_CNUMWS,ACDNWS(L)*AIRM(L))
c     if(AIJ(I,J,IJ_3dNWS).gt.0.)write(6,*)"OUTDRV",AIJ(I,J,IJ_3dNWS)
c    * ,ACDNWS(L),NLSW,itime,l
         ENDIF

        IF(NLSI.ge.1) then
         AIJ(I,J,IJ_3dNIS)=AIJ(I,J,IJ_3dNIS)+ACDNIS(L)
         AIJ(I,J,IJ_3dRIS)=AIJ(I,J,IJ_3dRIS)+AREIS(L)
         AIJ(I,J,IJ_3dLIS)=AIJ(I,J,IJ_3dLIS)+ALWIS(L)
         call inc_ajl(i,j,l,JL_CNUMIS,ACDNIS(L)*AIRM(L))
         AIJL(I,J,L,IJL_REIS)= AIJL(I,J,L,IJL_REIS)+AREIS(L)
         AIJL(I,J,L,IJL_CDIS)= AIJL(I,J,L,IJL_CDIS)+ACDNIS(L)
         AIJL(I,J,L,IJL_CWIS)= AIJL(I,J,L,IJL_CWIS)+ALWIS(L)
        ENDIF

        ENDDO

#endif
cQCON q2 = sum((Q(I,J,:)+WMX(:))*AIRM(:))*100.*BYGRAV+PRCP
cQCON if (abs(q0-q2).gt.1d-13) print*,"pr1",i,j,q0,q1,q2,prcp,prcpss*100
cQCON*     *bygrav

#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
        n = ntix(nx)
        do l=1,lp50
#ifndef SKIP_TRACER_DIAGS
           if (itcon_ss(n).gt.0) call inc_diagtcb(i,j,
     *          tm(l,nx)-trm(i,j,l,n)*fssl(l)
#ifdef TRACERS_WATER
     &         + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
#endif
     *          ,itcon_ss(n),n)
          call inc_tajln(i,j,l,jlnt_lscond,n,
     &         tm(l,nx)-trm(i,j,l,n)*fssl(l)
#ifdef TRACERS_WATER
     &         + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
#endif
     &         )
#endif  /*SKIP_TRACER_DIAGS*/
#ifdef TRACERS_WATER
          trwm(i,j,l,n) = trwml(nx,l)
#endif
          trm(i,j,l,n) = tm(l,nx)+tmsave(l,nx)*(1.-fssl(l))
          trmom(:,i,j,l,n) = tmom(:,l,nx)+tmomsv(:,l,nx)*(1.-fssl(l))
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
          if (trname(n).eq."SO2".or.trname(n).eq."SO4".or.trname(n).eq."
     *         H2O2_s") then
            call inc_tajls(i,j,l,jls_incloud(1,n),
     *           dt_sulf_mc(n,l)*(1.-fssl(l)))
            call inc_tajls(i,j,l,jls_incloud(2,n),dt_sulf_ss(n,l))
          if (ijts_aq(n).gt.0) then
            taijs(i,j,ijts_aq(n))=taijs(i,j,ijts_aq(n))+
     *           dt_sulf_mc(n,l)*(1.-fssl(l))+dt_sulf_ss(n,l)
          endif
          end if
#endif
#ifdef TRACERS_AMP
           if (trname(n).eq."M_ACC_SU") then
           AQsulfRATE(i,j,l)=  dt_sulf_mc(n,l)+dt_sulf_ss(n,l)
           endif
#endif
        end do
#ifdef TRACERS_WATER
        trprec(n,i,j) = trprec(n,i,j)+trprss(nx)
C**** diagnostics
        if (dowetdep(n)) then
#ifndef SKIP_TRACER_DIAGS
          if (jls_prec(1,n).gt.0) call inc_tajls(i,j,1,jls_prec(1,n),
     *         trprec(n,i,j)*byaxyp(i,j))
          if (jls_prec(2,n).gt.0) call inc_tajls(i,j,1,jls_prec(2,n),
     *         trprec(n,i,j)*focean(i,j)*byaxyp(i,j))
          taijn(i,j,tij_prec,n) =taijn(i,j,tij_prec,n) +
     *         trprec(n,i,j)*byaxyp(i,j)
#ifdef TRACERS_COSMO
          if (n .eq. n_Be7) BE7W_acc(i,j)=BE7W_acc(i,j)+
     *         trprec(n,i,j)*byaxyp(i,j)
#endif
#endif /*SKIP_TRACER_DIAGS*/
#ifdef TRDIAG_WETDEPO
c     ..........
c     accumulates special wet depo diagnostics
c     ..........
          IF (diag_wetdep == 1) THEN
            DO l=1,lmcmax
              IF (jls_trdpmc(1,n) > 0) call inc_tajls(i,j,l,jls_trdpmc(1
     *             ,n),trcond_mc(l,nx))
              IF (jls_trdpmc(2,n) > 0) call inc_tajls(i,j,l,jls_trdpmc(2
     *             ,n),trdvap_mc(l,nx))
              IF (jls_trdpmc(3,n) > 0) call inc_tajls(i,j,l,jls_trdpmc(3
     *             ,n),trflcw_mc(l,nx))
              IF (jls_trdpmc(4,n) > 0) call inc_tajls(i,j,l,jls_trdpmc(4
     *             ,n),trprcp_mc(l,nx))
              IF (jls_trdpmc(5,n) > 0) call inc_tajls(i,j,l,jls_trdpmc(5
     *             ,n),trnvap_mc(l,nx))
              IF (jls_trdpmc(6,n) > 0) call inc_tajls(i,j,l,jls_trdpmc(6
     *             ,n),trwash_mc(l,nx))
            END DO
            IF (ijts_trdpmc(1,n) > 0) taijs(i,j,ijts_trdpmc(1,n))
     &          =taijs(i,j,ijts_trdpmc(1,n))+SUM(trcond_mc(1:lmcmax,nx))
            IF (ijts_trdpmc(2,n) > 0) taijs(i,j,ijts_trdpmc(2,n))
     &          =taijs(i,j,ijts_trdpmc(2,n))+SUM(trdvap_mc(1:lmcmax,nx))
            IF (ijts_trdpmc(3,n) > 0) taijs(i,j,ijts_trdpmc(3,n))
     &          =taijs(i,j,ijts_trdpmc(3,n))+SUM(trflcw_mc(1:lmcmax,nx))
            IF (ijts_trdpmc(4,n) > 0) taijs(i,j,ijts_trdpmc(4,n))
     &          =taijs(i,j,ijts_trdpmc(4,n))+SUM(trprcp_mc(1:lmcmax,nx))
            IF (ijts_trdpmc(5,n) > 0) taijs(i,j,ijts_trdpmc(5,n))
     &          =taijs(i,j,ijts_trdpmc(5,n))+SUM(trnvap_mc(1:lmcmax,nx))
            IF (ijts_trdpmc(6,n) > 0) taijs(i,j,ijts_trdpmc(6,n))
     &          =taijs(i,j,ijts_trdpmc(6,n))+SUM(trwash_mc(1:lmcmax,nx))
            DO l=1,lp50
              IF (jls_trdpls(1,n) > 0) call inc_tajls(i,j,l,jls_trdpls(1
     *             ,n),trwash_ls(l,nx))
              IF (jls_trdpls(2,n) > 0) call inc_tajls(i,j,l,jls_trdpls(2
     *             ,n),trprcp_ls(l,nx))
              IF (jls_trdpls(3,n) > 0) call inc_tajls(i,j,l,jls_trdpls(3
     *             ,n),trclwc_ls(l,nx))
              IF (jls_trdpls(4,n) > 0) call inc_tajls(i,j,l,jls_trdpls(4
     *             ,n),trevap_ls(l,nx))
              IF (jls_trdpls(5,n) > 0) call inc_tajls(i,j,l,jls_trdpls(5
     *             ,n),trclwe_ls(l,nx))
              IF (jls_trdpls(6,n) > 0) call inc_tajls(i,j,l,jls_trdpls(6
     *             ,n),trcond_ls(l,nx))
            END DO
            IF (ijts_trdpls(1,n) > 0) taijs(i,j,ijts_trdpls(1,n))
     &           =taijs(i,j,ijts_trdpls(1,n))+SUM(trwash_ls(1:lp50,nx))
            IF (ijts_trdpls(2,n) > 0) taijs(i,j,ijts_trdpls(2,n))
     &           =taijs(i,j,ijts_trdpls(2,n))+SUM(trprcp_ls(1:lp50,nx))
            IF (ijts_trdpls(3,n) > 0) taijs(i,j,ijts_trdpls(3,n))
     &           =taijs(i,j,ijts_trdpls(3,n))+SUM(trclwc_ls(1:lp50,nx))
            IF (ijts_trdpls(4,n) > 0) taijs(i,j,ijts_trdpls(4,n))
     &           =taijs(i,j,ijts_trdpls(4,n))+SUM(trevap_ls(1:lp50,nx))
            IF (ijts_trdpls(5,n) > 0) taijs(i,j,ijts_trdpls(5,n))
     &           =taijs(i,j,ijts_trdpls(5,n))+SUM(trclwe_ls(1:lp50,nx))
            IF (ijts_trdpls(6,n) > 0) taijs(i,j,ijts_trdpls(6,n))
     &           =taijs(i,j,ijts_trdpls(6,n))+SUM(trcond_ls(1:lp50,nx))
          END IF
#endif
#ifdef TRACERS_DUST
          IF (adiurn_dust == 1) THEN
            DO kr=1,Ndiupt
              IF(i == ijdd(1,kr) .AND. j == ijdd(2,kr)) THEN
                SELECT CASE (trname(n))
                CASE ('Clay','Silt1','Silt2','Silt3','Silt4')
                  tmp(idd_wet)=+trprec(n,i,j)*byaxyp(i,j)/Dtsrc
                  ADIURN(IDXD(:),KR,IH)=ADIURN(IDXD(:),KR,IH)+
     &                 TMP(IDXD(:))
#ifndef NO_HDIURN
                  HDIURN(IDXD(:),KR,IHM)=HDIURN(IDXD(:),KR,IHM)+
     &                 TMP(IDXD(:))
#endif
                END SELECT
              END IF
            END DO
          END IF
#endif
        end if
#endif
      end do
#endif

#ifndef TRACERS_WATER
c     ..........
c     call simple wet deposition scheme for dust/mineral tracers
c     ..........

#ifdef TRACERS_DUST
      n_fidx=n_clay
#else
#ifdef TRACERS_MINERALS
      n_fidx=n_clayilli
#else
#ifdef TRACERS_QUARZHEM
      n_fidx=n_sil1quhe
#endif
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)

      DO n=1,Ntm_dust
        n1=n_fidx+n-1
        DO l=1,Lm
          tm_dust(l,n)=trm(i,j,l,n1)
          tmom_dust(:,l,n)=trmom(:,i,j,l,n1)
        END DO
      END DO

      CALL dust_wet(i,j)

      DO n=1,Ntm_dust
        n1=n_fidx+n-1
        trprec_dust(n,i,j)=0.D0
        DO l=1,Lm
          if (itcon_wt(n).gt.0) call inc_diagtcb(i,j,
     *          tm_dust(l,n)-trm(i,j,l,n1),itcon_wt(n),n)
          trm(i,j,l,n1)=tm_dust(l,n)
          trmom(:,i,j,l,n1)=tmom_dust(:,l,n)
          trprec_dust(n,i,j)=trprec_dust(n,i,j)+trprc_dust(l,n)
          call inc_tajls(i,j,l,jls_wet(n1),trprc_dust(l,n))
          taijs(i,j,ijts_wet(n1))=taijs(i,j,ijts_wet(n1))
     &         +trprc_dust(l,n)
        END DO
      END DO

#endif

#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN
        DO n=1,Ntm_dust
          DO kr=1,Ndiupt
            IF(i == ijdd(1,kr) .AND. j == ijdd(2,kr)) THEN
              SELECT CASE (trname(n))
              CASE ('Clay','Silt1','Silt2','Silt3','Silt4')
                tmp(idd_wet)=+trprec_dust(n,i,j)*byaxyp(i,j)/Dtsrc
                ADIURN(IDXD(:),KR,IH)=ADIURN(IDXD(:),KR,IH)+
     &               TMP(IDXD(:))
#ifndef NO_HDIURN
                HDIURN(IDXD(:),KR,IHM)=HDIURN(IDXD(:),KR,IHM)+
     &               TMP(IDXD(:))
#endif
              END SELECT
            END IF
          END DO
        END DO
      END IF
#endif
#endif

      END DO
C**** END OF MAIN LOOP FOR INDEX I

       ! Burn random numbers for later longitudes here.
       ! Actual generation of random numbers is in CLOUDS2.f::ISCCP_CLOUD_TYPES
      if (isccp_diags.eq.1) then
        CALL BURN_RANDOM(nij_after_i1(I_1)*NCOL*(LM+1))
      end if

C****
Cred*           Reduced Arrays 3
C****
         DO L=1,LM
         DO I=I_0,I_1
            WM(I,J,L) = WMIL(I,L)
            T3MOM(:,I,J,L) = TMOMIL(:,I,L)
            Q3MOM(:,I,J,L) = QMOMIL(:,I,L)
         END DO
         END DO
Cred*       end Reduced Arrays 3
      END DO
C**** END OF MAIN LOOP FOR INDEX J
!xOMP  END PARALLEL DO
C****
C
C     WAS THERE AN ERROR IN SUBSID ??
C
      IF(ICKERR.NE.0)  THEN
         WRITE(6,*)  'SUBSID ERROR: ABS(C) > 1'
         call stop_model('SUBSID ERROR: ABS(C) > 1',255)
      END IF
C
C     WAS THERE AN ERROR IN ISCCP CLOUD TYPING ??
C
      IF(JCKERR.NE.0)  THEN
         WRITE(6,*)  'ISCCP CLOUD TYPING ERROR'
         call stop_model('ISCCP CLOUD TYPING ERROR',255)
      END IF

#ifdef SKIP_TRACER_DIAGS
#ifdef TRACERS_WATER
      call trac_accum_clouds
#endif
#endif

C
C     NOW UPDATE THE MODEL WINDS
C
#ifndef SCM
      call avg_replicated_duv_to_vgrid(ukm,vkm,kmax_nonpolar,
     &     ukmsp,vkmsp,ukmnp,vkmnp)
#else
      I=I_TARG
      J=J_TARG
      DO L=1,LM
         DO K=1,2 ! KMAXJ(J)
            U(I,J,L)=U(I,J,L)+UKM(K,L,I,J)
            V(I,J,L)=V(I,J,L)+VKM(K,L,I,J)
         END DO
      END DO
#endif

C**** ADD IN CHANGE OF MOMENTUM BY MOIST CONVECTION AND CTEI
      UASV(:,I_0:I_1,J_0S:J_1S) = UA(:,I_0:I_1,J_0S:J_1S)
      call recalc_agrid_uv ! add option for tendency computation?
      DO J=J_0S,J_1S
      DO I=I_0,I_1
      DO L=1,LM
        call inc_ajl(i,j,l,JL_DAMMC,
     &       (UA(L,I,J)-UASV(L,I,J))*PDSIG(L,I,J))
      END DO
      END DO
      END DO

      if (isccp_diags.eq.1) CALL RINIT(seed) ! reset random number sequ.

  415 FORMAT(1X,'W500 AT I=21 L=5 TIME= ',I10/,1X,10F8.3/,1X,10F8.3)
  420 FORMAT(1X,'ENT  AT I=21 L=5'/,1X,10F8.2/,1X,10F8.2)

      RETURN
      END SUBROUTINE CONDSE

      SUBROUTINE init_CLD
!@sum  init_CLD initialises parameters for MSTCNV and LSCOND
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE CONSTANT, only : grav,by3,radian
      USE MODEL_COM, only : jm,lm,dtsrc,ls1,plbot,pednl00
      USE DOMAIN_DECOMP_ATM, only : GRID, AM_I_ROOT
      USE GEOM, only : lat2d, kmaxj

      USE CLOUDS, only : lmcm,bydtsrc,xmass,brcld,bybr,U00wtrX,U00ice
     *  ,U00a,U00b       ! tuning knobs to replace U00ice and U00wtrX
     *  ,HRMAX,ISC,lp50,RICldX,RWCldOX,xRIcld,do_blU00,tautab,invtau
     *  ,funio_denominator,autoconv_multiplier,radius_multiplier
     *  ,entrainment_cont1,entrainment_cont2
     &  ,RA,UM,VM,UM1,VM1,U_0,V_0
      USE CLOUDS_COM, only : llow,lmid,lhi
     &     ,isccp_reg2d,UKM,VKM
      USE DIAG_COM, only : nisccp,isccp_reg,isccp_late
     &     ,isccp_diags,ntau,npres
      USE PARAM
      USE FILEMANAGER, only : openunit, closeunit

      IMPLICIT NONE
      REAL*8 PLE
      INTEGER L,I,J,n,iu_ISCCP
      INTEGER :: I_0,I_1,J_0,J_1, I_0H,I_1H,J_0H,J_1H
      CHARACTER TITLE*80

      I_0 =GRID%I_STRT
      I_1 =GRID%I_STOP
      J_0 =GRID%J_STRT
      J_1 =GRID%J_STOP
      I_0H =GRID%I_STRT_HALO
      I_1H =GRID%I_STOP_HALO
      J_0H =GRID%J_STRT_HALO
      J_1H =GRID%J_STOP_HALO

c
c allocate space for the varying number of staggered
c wind data to be vertically mixed by clouds on the A grid
c
      n = maxval(kmaxj(j_0:j_1))
      allocate(RA(n))
      allocate(UM(n,lm),VM(n,lm),UM1(n,lm),VM1(n,lm))
      allocate(U_0(n,lm),V_0(n,lm))
      n = minval(kmaxj(j_0:j_1))
      allocate(UKM(n,lm,i_0h:i_1h,j_0h:j_1h),
     &         VKM(n,lm,i_0h:i_1h,j_0h:j_1h))

      call sync_param( 'U00wtrX', U00wtrX )
      call sync_param( 'U00ice', U00ice )
      call sync_param( 'U00a', U00a )
      call sync_param( 'U00b', U00b )
      call sync_param( "LMCM", LMCM )
      call sync_param( "HRMAX", HRMAX )
      call sync_param( "RICldX", RICldX )
      xRIcld = .001d0*(RICldX-1.d0)
      call sync_param( "RWCldOX", RWCldOX )
      call sync_param( "ISC", ISC)
      call sync_param( "do_blU00", do_blU00)
      call sync_param( "funio_denominator",funio_denominator)
      call sync_param( "autoconv_multiplier",autoconv_multiplier)
      call sync_param( "radius_multiplier",radius_multiplier)
      call sync_param( "entrainment_cont1",entrainment_cont1)
      call sync_param( "entrainment_cont2",entrainment_cont2)

      IF(LMCM.LT.0) LMCM = LS1-1
      call set_param( "LMCM", LMCM, 'o' )

      BYDTsrc=1./DTsrc
      XMASS=0.1d0*DTsrc*GRAV

      BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3

C**** SEARCH FOR THE 50 MB LEVEL
      LP50=LM
      DO L=LM-1,1,-1
        PLE=.25*(PEDNL00(L)+2.*PEDNL00(L+1)+PEDNL00(L+2))
        IF (PLE.LT.50.) LP50=L
      END DO
      if (AM_I_ROOT())  write(6,*)
     *     "Maximum level for LSCOND calculations (50mb): ",LP50

C**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS (MATCHES ISCCP DEFNs)
      DO L=1,LM
        LLOW=L
        IF (.5*(PLbot(L+1)+PLbot(L+2)).LT.680.) EXIT
      END DO
      DO L=LLOW+1,LM
        LMID=L
        IF (.5*(PLbot(L+1)+PLbot(L+2)).LT.440.) EXIT
      END DO
      LHI=LM
      IF (LMID+1.GT.LHI) LHI=LMID+1
      if (AM_I_ROOT()) WRITE (6,47) LLOW,LLOW+1,LMID,LMID+1,LHI
 47   FORMAT (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN',
     *     ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)

C**** Define regions for ISCCP diagnostics

c allocate/define distributed 2D ISCCP arrays
c      if (isccp_diags.eq.1) then
        allocate(isccp_reg2d(i_0h:i_1h,j_0h:j_1h))
        do j=j_0,j_1
        do i=i_0,i_1
          isccp_reg2d(i,j)=0
          do n=1,nisccp
           if(dble(nint(lat2d(i,j)/radian)).ge.isccp_late(n) .and.
     &        dble(nint(lat2d(i,j)/radian)).lt.isccp_late(n+1)) then
              isccp_reg2d(i,j)=n
              exit
           endif
          enddo
        enddo
        enddo
c      endif

C**** Read in tau/invtau tables for ISCCP calculations
      call openunit("ISCCP",iu_ISCCP,.true.,.true.)
      read(iu_ISCCP) title,tautab,invtau
      if (AM_I_ROOT())  write(6,*) "Read ISCCP:",trim(title)
      call closeunit(iu_ISCCP)
      END SUBROUTINE init_CLD
