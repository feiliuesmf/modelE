#include "rundeck_opts.h"

      SUBROUTINE CONDSE
!@sum   CONDSE driver for moist convection AND large-scale condensation
!@auth  M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver   1.0 (taken from CB265)
!@calls CLOUDS:MSTCNV,CLOUDS:LSCOND
      USE CONSTANT, only : bygrav,lhm,rgas,grav,tf,lhe,lhs,sha,deltx
     *     ,teeny
      USE MODEL_COM, only : im,jm,lm,p,u,v,t,q,wm,JHOUR,fearth
     *     ,ls1,psf,ptop,dsig,bydsig,jeq,sig,DTsrc,ftype,jdate
     *     ,ntype,itime,fim,focean,fland,flice
#ifdef TRACERS_AEROSOLS_Koch
     *     ,jyear,jmon
#endif
      USE DOMAIN_DECOMP, only : HALO_UPDATE,GRID,GET
      USE DOMAIN_DECOMP, only : CHECKSUM, NORTH,SOUTH
      USE DOMAIN_DECOMP, only : HALO_UPDATE_COLUMN,CHECKSUM_COLUMN
      USE DOMAIN_DECOMP, only : GLOBALSUM,AM_I_ROOT
      USE QUSDEF, only : nmom
      USE SOMTQ_COM, only : t3mom=>tmom,q3mom=>qmom
      USE GEOM, only : bydxyp,dxyp,imaxj,kmaxj,ravj,idij,idjj
      USE RANDOM
      USE CLOUDS_COM, only : ttold,qtold,svlhx,svlat,rhsav,cldsav
#ifdef CLD_AER_CDNC
     *     ,oldno,oldnl,smfpm
#endif
     *     ,tauss,taumc,cldss,cldmc,csizmc,csizss,ddm1,airx,lmc
      USE DAGCOM, only : aj,areg,aij,ajl,ail,adiurn,jreg,ij_pscld,
     *     ij_pdcld,ij_scnvfrq,ij_dcnvfrq,ij_wmsum,ij_snwf,ij_prec,
     *     ij_neth,ij_f0oc,j_eprcp,j_prcpmc,j_prcpss,il_mceq,j5s,j5n,
     *     ijdd,idd_pr,idd_ecnd,idd_mcp,idd_dmc,idd_smc,idd_ssp,
     &     jl_mcmflx,jl_sshr,jl_mchr,jl_dammc,jl_rhe,jl_mchphas
     *     ,jl_mcdtotw,jl_mcldht,jl_mcheat,jl_mcdry,ij_ctpi,ij_taui
     *     ,ij_lcldi,ij_mcldi,ij_hcldi,ij_tcldi,ij_sstabx,isccp_diags
     *     ,ndiupt,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss,hdiurn,
     *     ntau,npres,aisccp,isccp_reg
#ifdef CLD_AER_CDNC
     *     ,jl_cnumwm,jl_cnumws,jl_cnumim,jl_cnumis
     *     ,ij_3dnwm,ij_3dnws,ij_3dnim,ij_3dnis
     *     ,ij_3drwm,ij_3drws,ij_3drim,ij_3dris
#endif
#ifdef TRACERS_ON
      USE TRACER_COM, only: itime_tr0,TRM,TRMOM,NTM,trname
#ifdef TRACERS_WATER
     *     ,trwm,trw0,dowetdep
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE LIGHTNING, only : RNOx_lgt
#endif
      USE TRACER_DIAG_COM,only: tajln,jlnt_mc,jlnt_lscond,itcon_mc
     *     ,itcon_ss
#ifdef TRACERS_WATER
     *     ,jls_prec,taijn,tajls,tij_prec
#ifdef TRACERS_AEROSOLS_Koch
     *     ,jls_incloud,ijts_aq,taijs
#endif
#endif
      USE CLOUDS, only : tm,tmom ! local  (i,j)
     *     ,ntx,ntix              ! global (same for all i,j)
#ifdef TRACERS_WATER
     *     ,trwml,trsvwml,trprmc,trprss
#ifdef TRACERS_AEROSOLS_Koch
     *     ,dt_sulf_mc,dt_sulf_ss
#endif
#endif
#endif
      USE CLOUDS, only : BYDTsrc,mstcnv,lscond ! glb var & subroutines
     *     ,airm,byam,etal,sm,smom,qm,qmom,isc,dxypj,lp50
     *     ,tl,ris,ri1,ri2,mcflx,sshr,dgdsm,dphase,dtotw,dqcond,dctei
     *     ,wml,sdl,u_0,v_0,um,vm,qs,us,vs,dcl,airxl,prcpss,hcndss
     *     ,prcpmc,pearth,ts,taumcl,cldmcl,svwmxl,svlatl,svlhxl,dgdqm
     *     ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum
     *     ,aq,dpdt,th,ql,wmx,ttoldl,rh,taussl,cldssl,cldsavl,rh1
     *     ,kmax,ra,pl,ple,plk,rndssl,lhp,pland,debug,ddmflx
#ifdef CLD_AER_CDNC
     *     ,acdnwm,acdnim,acdnws,acdnis,arews,arewm,areis,areim
     *     ,nlsw,nlsi,nmcw,nmci
     *     ,oldcdo,oldcdl,smfpml
#endif
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg,tgvavg,qgavg,dclev
      USE DYNAMICS, only : pk,pek,pmid,pedn,sd_clouds,gz,ptold,pdsig
     *     ,ltropo,dke
      USE SEAICE_COM, only : rsi
      USE GHYCOM, only : snoage
      USE LAKES_COM, only : flake
      USE FLUXES, only : prec,eprec,precss,gtemp
#ifdef TRACERS_WATER
     *     ,trprec
#endif
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE

#ifdef TRACERS_ON
!@var tmsave holds tracer value (for diagnostics)
      REAL*8 tmsave(lm,ntm),
     *       dtr_mc(GRID%J_STRT_HALO:GRID%J_STOP_HALO,ntm),
     *       dtr_ss(GRID%J_STRT_HALO:GRID%J_STOP_HALO,ntm)
      INTEGER NX
#ifdef TRACERS_SPECIAL_Shindell
!@var Lfreeze Lowest level where temperature is below freezing (TF)
      INTEGER Lfreeze
#endif
#endif

!@var UC,VC,UZM,VZM velocity work arrays
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
     *        :: UC,VC
      REAL*8, DIMENSION(2,LM) :: UZM,VZM

!@param ENTCON fractional rate of entrainment (km**-1)
      REAL*8,  PARAMETER :: ENTCON = .2d0

      INTEGER I,J,K,L,N  !@var I,J,K,L,N loop variables
      INTEGER JR,KR,ITYPE,IT,IH,LP850,LP600,IHM
!@var JR = JREG(I,J)
!@var KR index for regional diagnostics
!@var ITYPE index for snow age
!@var IT index for surface types
!@var LP850 layer near 850 mb
!@var LP600 layer near 600 mb
!@var LERR,IERR error reporting
      INTEGER :: LERR, IERR
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID

      REAL*8 :: HCNDMC,PRCP,TPRCP,EPRCP,ENRGP,WMERR,ALPHA1,ALPHA2,ALPHAS
      REAL*8 :: DTDZ,DTDZS,DUDZ,DVDZ,DUDZS,DVDZS,THSV,THV1,THV2,QG,TGV
      REAL*8 :: DH1S,BYDH1S,DH12,BYDH12,DTDZG,DUDZG,DVDZG,SSTAB,DIFT,CSC
     *     ,E,E1,ep,TSV
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
      integer itau,itrop,nbox
C****

C 
Cred*                       Reduced Arrays 1                 *********
C        not clear yet whether they still speed things up
      REAL*8  GZIL(IM,LM), SD_CLDIL(IM,LM), WMIL(IM,LM)
      REAL*8  TMOMIL(NMOM,IM,LM),  QMOMIL(NMOM,IM,LM)
Cred*                   end Reduced Arrays 1
      INTEGER ICKERR, JCKERR, JERR, seed, NR
      REAL*8  RNDSS(3,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),xx
!     REAL*8  AJEQIL(J5N-J5S+1,IM,LM),
      REAL*8  AJEQIL(GRID%J_STRT_HALO:GRID%J_STOP_HALO,IM,LM),
     *        AREGIJ(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,3)
      REAL*8  UKP1(IM,LM), VKP1(IM,LM), UKPJM(IM,LM),VKPJM(IM,LM)
      REAL*8  UKM(4,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *        VKM(4,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
      INTEGER :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,J_0STG,J_1STG
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

C**** Initialize
      AJEQIL(:,:,:)=0.
C**** define local grid
      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_HALO=J_0H,    J_STOP_HALO=J_1H,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE        )

C 
C     OBTAIN RANDOM NUMBERS FOR PARALLEL REGION
C 
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        DO L=LP50,1,-1
          DO NR=1,3
            RNDSS(NR,L,I,J) = RANDU(xx)
          END DO
        END DO
C     Do not bother to save random numbers for isccp_clouds
      END DO
      END DO
C     But save the current seed in case isccp_routine is activated
      if (isccp_diags.eq.1) CALL RFINAL(seed)
C
C**** UDATE HALOS of U and V FOR DISTRIBUTED PARALLELIZATION
      CALL CHECKSUM   (grid, U, __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, U, from= NORTH)
      CALL CHECKSUM   (grid, V, __LINE__, __FILE__)
      CALL HALO_UPDATE(grid, V, from= NORTH)
C 
C**** SAVE UC AND VC, AND ZERO OUT CLDSS AND CLDMC
      UC=U
      VC=V
C**** COMPUTE ZONAL MEAN U AND V AT POLES
      DO L=1,LM
        UZM(1,L)=0.
        UZM(2,L)=0.
        VZM(1,L)=0.
        VZM(2,L)=0.
      ENDDO
      DO L=1,LM
        IF(HAVE_SOUTH_POLE) THEN
          DO I=1,IM
            UZM(1,L)=UZM(1,L)+UC(I,2,L)
            VZM(1,L)=VZM(1,L)+VC(I,2,L)
          ENDDO
        ENDIF
        IF(HAVE_NORTH_POLE) THEN
          DO I=1,IM
            UZM(2,L)=UZM(2,L)+UC(I,JM,L)
            VZM(2,L)=VZM(2,L)+VC(I,JM,L)
          ENDDO
        ENDIF
        UZM(1,L)=UZM(1,L)/FIM
        UZM(2,L)=UZM(2,L)/FIM
        VZM(1,L)=VZM(1,L)/FIM
        VZM(2,L)=VZM(2,L)/FIM
      ENDDO
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
#ifdef TRACERS_SPECIAL_Shindell
C**** Make sure the NOx from lightning is initialized:
      RNOx_lgt=0.
#endif
#endif
C****
C**** MAIN J LOOP
C****
       ICKERR=0
       JCKERR=0
!$OMP  PARALLEL DO PRIVATE (
#ifdef TRACERS_ON
!$OMP*  NX,tmsave,
#endif
!$OMP*  ALPHAS,ALPHA1,ALPHA2,AT,BYDH1S,BYDH12, CC,CONV,CTP,
!$OMP*  DH1S,DH12,DTDZ,DTDZG,DTDZS,DUDZ,DUDZG,DUDZS,DVDZ,DVDZG,DVDZS,
!$OMP*  DTAU_S,DTAU_C,DEM_S,DEM_C, FQ_ISCCP, ENRGP,EPRCP,
!$OMP*  HCNDMC, I,ITYPE,IT,ITAU, IDI,IDJ,
#ifdef TRACERS_SPECIAL_Shindell
!$OMP*  Lfreeze,
#endif
!$OMP*  ITROP,IERR, J,JERR, K,KR, L,LERR, N,NBOX, PRCP,PFULL,PHALF,
!$OMP*  GZIL, SD_CLDIL, WMIL, TMOMIL, QMOMIL,        ! reduced arrays
!$OMP*  QG,QV, SKT,SSTAB, TGV,TPRCP,THSV,THV1,THV2,TAUOPT,TSV, WMERR,
!$OMP*  LP600,LP850,CSC,DIFT, E,E1,ep)
!$OMP*    SCHEDULE(DYNAMIC,2)
!$OMP*    REDUCTION(+:ICKERR,JCKERR)
C 
      DO J=J_0,J_1
C 
Cred* Reduced Arrays 2
C 
        DXYPJ=DXYP(J)
      DO L=1,LM
         GZIL(:,L) = GZ(:,J,L)
      END DO
      DO L=1,LM
         SD_CLDIL(:,L) = SD_CLOUDS(:,J,L)
      END DO
      DO L=1,LM
         WMIL(:,L) = WM(:,J,L)
      END DO
      DO L=1,LM
         TMOMIL(:,:,L) = T3MOM(:,:,J,L)
      END DO
      DO L=1,LM
         QMOMIL(:,:,L) = Q3MOM(:,:,J,L)
      END DO
Cred* end Reduced Arrays 2
#ifdef TRACERS_ON
      dtr_mc(j,:)=0. ; dtr_ss(j,:)=0.
#endif
      kmax = kmaxj(j)
C****
C**** MAIN I LOOP
C****
      DO I=1,IMAXJ(J)
cc       JR=JREG(I,J)  ! summing done outside parallel region
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

      DO K=1,KMAX
        RA(K)=RAVJ(K,J)
        IDI(K)=IDIJ(K,I,J)
        IDJ(K)=IDJJ(K,J)
      END DO
C**** PRESSURES, AND PRESSURE TO THE KAPA
      PL(:) =PMID(:,I,J)
      PLE(:)=PEDN(:,I,J)
      PLK(:)=PK(:,I,J)
      AIRM(:)=PDSIG(:,I,J)
      BYAM(:)=1./AIRM(:)
C**** other fields where L is the leading index
      SVLHXL(:)=SVLHX(:,I,J)
      TTOLDL(:)=TTOLD(:,I,J)
      CLDSAVL(:)=CLDSAV(:,I,J)
      RH(:)=RHSAV(:,I,J)
#ifdef CLD_AER_CDNC
        OLDCDL(:)=OLDNL(:,I,J)
        OLDCDO(:)=OLDNO(:,I,J)  ! OLDN is for rsf save
        SMFPML(:)=SMFPM(:,I,J)
#endif
      DPDT(1:LS1-1)=SIG(1:LS1-1)*(P(I,J)-PTOLD(I,J))*BYDTsrc
      DPDT(LS1:LM)=0.
      DO L=1,LM
C**** TEMPERATURES
        SM(L)  =T(I,J,L)*AIRM(L)
Cred    SMOM(:,L) =T3MOM(:,I,J,L)*AIRM(L)
        SMOM(:,L) =TMOMIL(:,I,L)*AIRM(L)
        TL(L)=T(I,J,L)*PLK(L)
C**** MOISTURE (SPECIFIC HUMIDITY)
        QM(L)  =Q(I,J,L)*AIRM(L)
Cred    QMOM(:,L) =Q3MOM(:,I,J,L)*AIRM(L)
        QMOM(:,L) =QMOMIL(:,I,L)*AIRM(L)
Cred    WML(L)=WM(I,J,L)
        WML(L)=WMIL(I,L)
C**** others
Cred    SDL(L)=SD_CLOUDS(I,J,L)*BYDXYP(J)
        SDL(L)=SD_CLDIL(I,L)*BYDXYP(J)
        IF(L.LE.LM-2)
Cred *    ETAL(L+1)=.5*ENTCON*(GZ(I,J,L+2)-GZ(I,J,L))*1.d-3*BYGRAV
     *    ETAL(L+1)=.5*ENTCON*(GZIL(I,L+2)-GZIL(I,L))*1.d-3*BYGRAV
      END DO
      ETAL(LM)=ETAL(LM-1)
      ETAL(1)=0.     ! not used
#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
      do l=1,lm
        tm(l,nx) = trm(i,j,l,ntix(nx))
        tmom(:,l,nx) = trmom(:,i,j,l,ntix(nx))
      end do
      end do
#endif
C**** SURROUNDING WINDS
      DO L=1,LM
        DO K=1,KMAX
          U_0(K,L) = UC(IDI(K),IDJ(K),L)
          V_0(K,L) = VC(IDI(K),IDJ(K),L)
          UM(K,L) = U_0(K,L)*AIRM(L)
          VM(K,L) = V_0(K,L)*AIRM(L)
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
C**** and again after LSCOND:
C****          = sum(WMX(:)*(LHE-SVLHXL(:))*AIRM(:))*100.*BYGRAV
C****

C**** MOIST CONVECTION
      CALL MSTCNV(IERR,LERR,I,J)

C**** Error reports
      if (ierr.gt.0) then
        write(6,*) "Error in moist conv: i,j,l=",i,j,lerr
ccc     if (ierr.eq.2) call stop_model("Subsid error: abs(c) > 1",255)
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
        HCNDMC=0.
        DO L=1,LMCMAX
          HCNDMC=HCNDMC+DGDSM(L)+DPHASE(L)
          AJL(J,L,JL_MCHR)=AJL(J,L,JL_MCHR)+DGDSM(L)*BYDSIG(L)
          AJL(J,L,JL_MCHPHAS)=AJL(J,L,JL_MCHPHAS)+DPHASE(L)*BYDSIG(L)
          AJL(J,L,JL_MCDTOTW)=AJL(J,L,JL_MCDTOTW)+DTOTW(L)*BYDSIG(L)
CCC       IF(J.GE.J5S.AND.J.LE.J5N) AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+
CCC  *         (DGDSM(L)+DPHASE(L))*(DXYP(J)*BYDSIG(L))
          IF(J.GE.J5S.AND.J.LE.J5N)     ! add in after parallel region
     *      AJEQIL(J,I,L) = (DGDSM(L)+DPHASE(L))*
     *                            (DXYP(J)*BYDSIG(L))
          AJL(J,L,JL_MCHEAT)=AJL(J,L,JL_MCHEAT)+
     &         (DPHASE(L)+DGDSM(L))*BYDSIG(L)
          AJL(J,L,JL_MCDRY)=AJL(J,L,JL_MCDRY)+
     &         (DQCOND(L)-DGDQM(L))*BYDSIG(L)
          AJL(J,L,JL_MCMFLX)=AJL(J,L,JL_MCMFLX)+MCFLX(L)
          AJL(J,L,JL_CLDMC) =AJL(J,L,JL_CLDMC) +CLDMCL(L)
          AJL(J,L,JL_CSIZMC)=AJL(J,L,JL_CSIZMC)+CSIZEL(L)*CLDMCL(L)
        END DO
        DO IT=1,NTYPE
          AJ(J,J_PRCPMC,IT)=AJ(J,J_PRCPMC,IT)+PRCPMC*FTYPE(IT,I,J)
        END DO
CCC     AREG(JR,J_PRCPMC)=AREG(JR,J_PRCPMC)+PRCPMC*DXYP(J)
        AREGIJ(I,J,1)=PRCPMC*DXYP(J)  ! add in after parallel region
        DO KR=1,NDIUPT
          IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
            ADIURN(IH,IDD_PR  ,KR)=ADIURN(IH,IDD_PR  ,KR)+PRCPMC
            ADIURN(IH,IDD_ECND,KR)=ADIURN(IH,IDD_ECND,KR)+HCNDMC
            ADIURN(IH,IDD_MCP ,KR)=ADIURN(IH,IDD_MCP ,KR)+PRCPMC
            ADIURN(IH,IDD_DMC ,KR)=ADIURN(IH,IDD_DMC ,KR)+CLDDEPIJ
            ADIURN(IH,IDD_SMC ,KR)=ADIURN(IH,IDD_SMC ,KR)+CLDSLWIJ
            HDIURN(IHM,IDD_PR  ,KR)=HDIURN(IHM,IDD_PR  ,KR)+PRCPMC
            HDIURN(IHM,IDD_ECND,KR)=HDIURN(IHM,IDD_ECND,KR)+HCNDMC
            HDIURN(IHM,IDD_MCP ,KR)=HDIURN(IHM,IDD_MCP ,KR)+PRCPMC
            HDIURN(IHM,IDD_DMC ,KR)=HDIURN(IHM,IDD_DMC ,KR)+CLDDEPIJ
            HDIURN(IHM,IDD_SMC ,KR)=HDIURN(IHM,IDD_SMC ,KR)+CLDSLWIJ
          END IF
        END DO
#ifdef CLD_AER_CDNC
        DO L =1,LM
        IF (NMCW.ge.1) then
         AIJ(I,J,IJ_3dNWM)=AIJ(I,J,IJ_3dNWM)+ACDNWM(L)   !/NMCW
         AIJ(I,J,IJ_3dRWM)=AIJ(I,J,IJ_3dRWM)+AREWM(L)   !/NMCW
         AJL(J,L,JL_CNUMWM)=AJL(J,L,JL_CNUMWM)+ACDNWM(L)   !/NMCW
        ENDIF
        IF (NMCI.ge.1) then
         AIJ(I,J,IJ_3dNIM)=AIJ(I,J,IJ_3dNIM)+ACDNIM(L)    !/NMCI
         AIJ(I,J,IJ_3dRIM)=AIJ(I,J,IJ_3dRIM)+AREIM(L)   !/NMCI
         AJL(J,L,JL_CNUMIM)=AJL(J,L,JL_CNUMIM)+ACDNIM(L)   !/NMCI
        ENDIF
        ENDDO
#endif
C**** ACCUMULATE PRECIP
        PRCP=PRCPMC*100.*BYGRAV
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
C**** NEED TO TAKE ACCOUNT OF LATENT HEAT THOUGH
        IF (TPRCP.gt.0) THEN
C         EPRCP=PRCP*TPRCP*SHW
          EPRCP=0.
          ENRGP=ENRGP+EPRCP
        ELSE
C         EPRCP=PRCP*TPRCP*SHI
          EPRCP=0.
          ENRGP=ENRGP+EPRCP-PRCP*LHM
          AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCP
        END IF

        DO L=1,LMCMAX
          T(I,J,L)=  SM(L)*BYAM(L)
          Q(I,J,L)=  QM(L)*BYAM(L)
        END DO
        CSIZMC(1:LMCMAX,I,J)=CSIZEL(1:LMCMAX)
        AIRX(I,J) = AIRXL*DXYP(J)
C**** level 1 downfdraft mass flux/rho (m/s)
        DDM1(I,J) = DDMFLX(1)*RGAS*TSV/(GRAV*PEDN(1,I,J)*DTSrc)
      END IF
#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
        n = ntix(nx)
        do l=1,lm
          dtr_mc(j,nx)=dtr_mc(j,nx)+(tm(l,nx)-trm(i,j,l,n))
#ifdef TRACERS_WATER
     *         + trsvwml(nx,l)
#endif
          tajln(j,l,jlnt_mc,n) = tajln(j,l,jlnt_mc,n) +
     &          (tm(l,nx)-trm(i,j,l,n))
#ifdef TRACERS_WATER
     *         + trsvwml(nx,l)
          trwml(nx,l) = trwm(i,j,l,n)+trsvwml(nx,l)
#endif
          tmsave(l,nx) = tm(l,nx) ! save for tajln(large-scale condense)
        end do
#ifdef TRACERS_WATER
        trprec(n,i,j) = trprmc(nx)
#endif
      end do
#endif
      LMC(1,I,J) = LMCMIN
      LMC(2,I,J) = LMCMAX+1
C****
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS
C****
      DO L=1,LM
        TL(L)=T(I,J,L)*PLK(L)
        TH(L)=T(I,J,L)
        QL(L)=Q(I,J,L)
      END DO
      WMX(:)=WML(:)+SVWMXL(:)
      AQ(:)=(QL(:)-QTOLD(:,I,J))*BYDTsrc
      RNDSSL(:,1:LP50)=RNDSS(:,1:LP50,I,J)
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
CAOO        IF (J.EQ.1) THEN
        IF(HAVE_SOUTH_POLE .AND. J.EQ.J_0) THEN
          DUDZ=(UZM(1,2)-UZM(1,1))*BYDH12
          DVDZ=(VZM(1,2)-VZM(1,1))*BYDH12
          DUDZS=(UZM(1,1)-US)*BYDH1S
          DVDZS=(VZM(1,1)-VS)*BYDH1S
        ENDIF
CAOO        IF (J.EQ.JM) THEN
        IF(HAVE_NORTH_POLE .AND. J.EQ.J_1) THEN
          DUDZ=(UZM(2,2)-UZM(2,1))*BYDH12
          DVDZ=(VZM(2,2)-VZM(2,1))*BYDH12
          DUDZS=(UZM(2,1)-US)*BYDH1S
          DVDZS=(VZM(2,1)-VS)*BYDH1S
        ENDIF
        IF(J.GT.1.AND.J.LT.JM) THEN
          DUDZ=(UC(IDI(1),IDJ(1),2)+UC(IDI(2),IDJ(2),2)+
     *         UC(IDI(3),IDJ(3),2)+UC(IDI(4),IDJ(4),2)-
     *         UC(IDI(1),IDJ(1),1)-UC(IDI(2),IDJ(2),1)-
     *         UC(IDI(3),IDJ(3),1)-UC(IDI(4),IDJ(4),1))*.25*BYDH12
          DVDZ=(VC(IDI(1),IDJ(1),2)+VC(IDI(2),IDJ(2),2)+
     *         VC(IDI(3),IDJ(3),2)+VC(IDI(4),IDJ(4),2)-
     *         VC(IDI(1),IDJ(1),1)-VC(IDI(2),IDJ(2),1)-
     *         VC(IDI(3),IDJ(3),1)-VC(IDI(4),IDJ(4),1))*.25*BYDH12
          DUDZS=(UC(IDI(1),IDJ(1),1)+UC(IDI(2),IDJ(2),1)+
     *         UC(IDI(3),IDJ(3),1)+UC(IDI(4),IDJ(4),1)-
     *         4.*US)*.25*BYDH1S
          DVDZS=(VC(IDI(1),IDJ(1),1)+VC(IDI(2),IDJ(2),1)+
     *         VC(IDI(3),IDJ(3),1)+VC(IDI(4),IDJ(4),1)-
     *         4.*VS)*.25*BYDH1S
        ENDIF
        DUDZG=.1d0*US
        DVDZG=.1d0*VS
        DTDZG=.1d0*(THSV-TGV/PEK(1,I,J))
        RIS=(GRAV*ALPHAS*DTDZG)/(DUDZG*DUDZG+DVDZG*DVDZG+teeny)
        RI1=(GRAV*ALPHA1*DTDZS)/(DUDZS*DUDZS+DVDZS*DVDZS)
        RI2=(GRAV*ALPHA2*DTDZ)/(DUDZ*DUDZ+DVDZ*DVDZ)
C       WRITE (6,*)'I,J,QG,TGV,THSV,RIS,RI1=',I,J,QG,TGV,THSV,RIS,RI1
      ENDIF

c**** uncomment lines marked ECON to check energy conservation
cECON  E = ( sum(TL(1:LP50)*AIRM(1:LP50))*SHA + sum(QL(1:LP50)
cECON *     *AIRM(1:LP50))*LHE +sum( (WML(1:LP50)*(LHE-SVLHXL(1:LP50))
cECON *     +SVWMXL(1:LP50)*(LHE-SVLATL(1:LP50)))*AIRM(1:LP50))  )*100.
cECON *     *BYGRAV

C**** LARGE-SCALE CLOUDS AND PRECIPITATION
      CALL LSCOND(IERR,WMERR,LERR,i,j)

cECON  E1 = ( sum(TL(1:LP50)*AIRM(1:LP50))*SHA + sum(QL(1:LP50)
cECON *     *AIRM(1:LP50))*LHE +sum(WMX(1:LP50)*(LHE-SVLHXL(1:LP50))
cECON *     *AIRM(1:LP50)) )*100.*BYGRAV

C**** Error reports
      IF (IERR.ne.0) WRITE(99,'(I10,3I4,A,D14.5,A)')
     *       Itime,I,J,LERR,' CONDSE:H2O<0',WMERR,' ->0'

C**** Accumulate diagnostics of LSCOND
         AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
         DO IT=1,NTYPE
           AJ(J,J_PRCPSS,IT)=AJ(J,J_PRCPSS,IT)+PRCPSS*FTYPE(IT,I,J)
         END DO
CCC      AREG(JR,J_PRCPSS)=AREG(JR,J_PRCPSS)+PRCPSS*DXYP(J)
         AREGIJ(I,J,2)=PRCPSS*DXYP(J)  ! add in after parallel region
         DO KR=1,NDIUPT
           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             ADIURN(IH,IDD_PR  ,KR)=ADIURN(IH,IDD_PR  ,KR)+PRCPSS
             ADIURN(IH,IDD_ECND,KR)=ADIURN(IH,IDD_ECND,KR)+HCNDSS
             ADIURN(IH,IDD_SSP ,KR)=ADIURN(IH,IDD_SSP ,KR)+PRCPSS
             HDIURN(IHM,IDD_PR  ,KR)=HDIURN(IHM,IDD_PR  ,KR)+PRCPSS
             HDIURN(IHM,IDD_ECND,KR)=HDIURN(IHM,IDD_ECND,KR)+HCNDSS
             HDIURN(IHM,IDD_SSP ,KR)=HDIURN(IHM,IDD_SSP ,KR)+PRCPSS
           END IF
         END DO

C**** TOTAL PRECIPITATION AND AGE OF SNOW
      PRCP=PRCP+PRCPSS*100.*BYGRAV
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
C**** NEED TO TAKE ACCOUNT OF LATENT HEAT THOUGH
      IF (LHP(1).ne.LHS) THEN
C       EPRCP=PRCPSS*100.*BYGRAV*TPRCP*SHW
        EPRCP=0.
        ENRGP=ENRGP+EPRCP
cECON    ep=0.
      ELSE
C       EPRCP=PRCPSS*100.*BYGRAV*TPRCP*SHI
        EPRCP=0.
        ENRGP=ENRGP+EPRCP-PRCPSS*100.*BYGRAV*LHM
cECON    ep=-PRCPSS*100.*BYGRAV*LHM
        AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCPSS*100.*BYGRAV
      END IF

cECON  if (abs(E-E1-ep).gt.0.01) print*,"energy err",i,j,E-E1-ep,
cECON *     E,E1,ep,prcpss,lhp(1)

C**** PRECIPITATION DIAGNOSTICS
        DO IT=1,NTYPE
          AJ(J,J_EPRCP,IT)=AJ(J,J_EPRCP,IT)+ENRGP*FTYPE(IT,I,J)
        END DO
CCC     AREG(JR,J_EPRCP)=AREG(JR,J_EPRCP)+ENRGP*DXYP(J)
        AREGIJ(I,J,3)=ENRGP*DXYP(J)  ! add in after parallel region
        AIJ(I,J,IJ_PREC)=AIJ(I,J,IJ_PREC)+PRCP
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+ENRGP
        IF (FOCEAN(I,J).gt.0) AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+
     *       ENRGP*FOCEAN(I,J)*(1.-RSI(I,J))

      IF(ENRGP.LT.0.) THEN ! MODIFY SNOW AGES AFTER SNOW FALL
        DO ITYPE=1,3
          SNOAGE(ITYPE,I,J)=SNOAGE(ITYPE,I,J)*EXP(-PRCP)
        END DO
      END IF

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

C**** set tg1 from GTEMP array (or save in SURFACE?)
c          skt=tf+tg1(i,j)
          skt=tf + (focean(i,j)+flake(i,j))*(1.-rsi(i,j))*gtemp(1,1,i,j)
     *         + (focean(i,j)+flake(i,j))*rsi(i,j)*gtemp(1,2,i,j)
     *         + flice(i,j)*gtemp(1,3,i,j)+fearth(i,j)*gtemp(1,4,i,j)
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

        call ISCCP_CLOUD_TYPES(pfull,phalf,qv,
     &       cc,conv,dtau_s,dtau_c,skt,
     &       at,dem_s,dem_c,itrop,fq_isccp,ctp,tauopt,nbox,jerr)
        if(jerr.ne.0) jckerr = jckerr + 1
C**** set ISCCP diagnostics
        if (nbox.gt.0) then
          AIJ(I,J,IJ_CTPI) = AIJ(I,J,IJ_CTPI) + ctp
          AIJ(I,J,IJ_TAUI) = AIJ(I,J,IJ_TAUI) + tauopt
          AIJ(I,J,IJ_TCLDI)= AIJ(I,J,IJ_TCLDI)+ 1.
        end if
C**** note LOW CLOUDS:       ipres=6,7, MID-LEVEL CLOUDS: ipres=4,5 ,
C****      HIGH CLOUDS:      ipres=1,2,3
C**** Sum over itau=2,ntau (itau=1 is no cloud)
        do itau=2,ntau
          AIJ(I,J,IJ_LCLDI)= AIJ(I,J,IJ_LCLDI) + fq_isccp(itau,6)
     *         + fq_isccp(itau,7)
          AIJ(I,J,IJ_MCLDI)= AIJ(I,J,IJ_MCLDI) + fq_isccp(itau,4)
     *         + fq_isccp(itau,5)
          AIJ(I,J,IJ_HCLDI)= AIJ(I,J,IJ_HCLDI) + fq_isccp(itau,1)
     *         + fq_isccp(itau,2) + fq_isccp(itau,3)
        end do
C**** Save area weighted isccp histograms
        n=isccp_reg(j)
        if (n.gt.0) then
          AISCCP(:,:,n) = AISCCP(:,:,n) + fq_isccp(:,:)*DXYP(j)
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
      SVLHX(:,I,J)=SVLHXL(:)
      CSIZSS(:,I,J)=CSIZEL(:)
      RHSAV(:,I,J)=RH(:)

#ifdef CLD_AER_CDNC
         OLDNL(:,I,J)=OLDCDL(:)
         OLDNO(:,I,J)=OLDCDO(:)
         SMFPM(:,I,J)=SMFPML(:)
#endif
      TTOLD(:,I,J)=TH(:)
      QTOLD(:,I,J)=QL(:)

      PREC(I,J)=PRCP            ! total precip mass (kg/m^2)
      EPREC(I,J)=ENRGP          ! energy of precipitation (J/m^2)
C**** The PRECSS array is only used if a distinction is being made
C**** between kinds of rain in the ground hydrology.
      PRECSS(I,J)=PRCPSS*100.*BYGRAV  ! large scale precip (kg/m^2)

      DO L=1,LM
        AJL(J,L,JL_SSHR)=AJL(J,L,JL_SSHR)+SSHR(L)
        AJL(J,L,JL_MCLDHT)=AJL(J,L,JL_MCLDHT)+DCTEI(L)
        AJL(J,L,JL_RHE)=AJL(J,L,JL_RHE)+RH1(L)
        AJL(J,L,JL_CLDSS) =AJL(J,L,JL_CLDSS) +CLDSSL(L)
        AJL(J,L,JL_CSIZSS)=AJL(J,L,JL_CSIZSS)+CSIZEL(L)*CLDSSL(L)

        T(I,J,L)=TH(L)
        Q(I,J,L)=QL(L)
C**** update moment changes
Cred    T3MOM(:,I,J,L)=SMOM(:,L)*BYAM(L)
Cred    Q3MOM(:,I,J,L)=QMOM(:,L)*BYAM(L)
Cred    WM(I,J,L)=WMX(L)
        TMOMIL(:,I,L)=SMOM(:,L)*BYAM(L)
        QMOMIL(:,I,L)=QMOM(:,L)*BYAM(L)
        WMIL(I,L)=WMX(L)

C**** UPDATE MODEL WINDS
CCC     DO K=1,KMAX       !  add in after parallel region (order)
CCC       U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
CCC  &         +(UM(K,L)*BYAM(L)-UC(IDI(K),IDJ(K),L))
CCC       V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
CCC  &         +(VM(K,L)*BYAM(L)-VC(IDI(K),IDJ(K),L))
CCC     ENDDO
         IF(J.EQ.1)  THEN
            DO K=1,IM ! KMAX
               UKP1(K,L)=(UM(K,L)*BYAM(L)-UC(IDI(K),IDJ(K),L))
               VKP1(K,L)=(VM(K,L)*BYAM(L)-VC(IDI(K),IDJ(K),L))
            END DO
         ELSE IF(J.EQ.JM)  THEN
            DO K=1,IM ! KMAX
               UKPJM(K,L)=(UM(K,L)*BYAM(L)-UC(IDI(K),IDJ(K),L))
               VKPJM(K,L)=(VM(K,L)*BYAM(L)-VC(IDI(K),IDJ(K),L))
            END DO
         ELSE
            DO K=1,4 ! KMAX
               UKM(K,I,J,L)=(UM(K,L)*BYAM(L)-UC(IDI(K),IDJ(K),L))
               VKM(K,I,J,L)=(VM(K,L)*BYAM(L)-VC(IDI(K),IDJ(K),L))
            END DO
         END IF
      ENDDO
#ifdef CLD_AER_CDNC
        DO L=1,LM
         IF (NLSW.ge.1) then
          AIJ(I,J,IJ_3dNWS)=AIJ(I,J,IJ_3dNWS)+ACDNWS(L) !/NLSW
          AIJ(I,J,IJ_3dRWS)=AIJ(I,J,IJ_3dRWS)+AREWS(L)  !/NLSW
          AJL(J,L,JL_CNUMWS)=AJL(J,L,JL_CNUMWS)+ACDNWS(L)  !/NLSW
c     if(AIJ(I,J,IJ_3dNWS).gt.1500.)write(6,*)"OUTDRV",AIJ(I,J,IJ_3dNWS)
c    * ,ACDNWS,NLSW,itime
         ENDIF
        IF(NLSI.ge.1) then
         AIJ(I,J,IJ_3dNIS)=AIJ(I,J,IJ_3dNIS)+ACDNIS(L)!/NLSI
         AIJ(I,J,IJ_3dRIS)=AIJ(I,J,IJ_3dRIS)+AREIS(L)!/NLSI
         AJL(J,L,JL_CNUMIS)=AJL(J,L,JL_CNUMIS)+ACDNIS(L)!/NLSI
        ENDIF
c     if(AIJ(I,J,IJ_3dNWS).gt.1500.)write(6,*)"ODRV",AIJ(I,J,IJ_3dNWS)
c    * ,ACDNWS(L),L
        ENDDO
#endif

#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
        n = ntix(nx)
        do l=1,lp50
          dtr_ss(j,nx)=dtr_ss(j,nx)+(tm(l,nx)-tmsave(l,nx))
#ifdef TRACERS_WATER
     &         + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
#endif
          trm(i,j,l,n) = tm(l,nx)
          trmom(:,i,j,l,n) = tmom(:,l,nx)
          tajln(j,l,jlnt_lscond,n) = tajln(j,l,jlnt_lscond,n) +
     &         (tm(l,nx)-tmsave(l,nx))
#ifdef TRACERS_WATER
     &         + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
          trwm(i,j,l,n) = trwml(nx,l)
#ifdef TRACERS_AEROSOLS_Koch
          if (trname(n).eq."SO2".or.trname(n).eq."SO4".or.trname(n).eq."
     *         H2O2_s") then
            tajls(j,l,jls_incloud(1,n))=tajls(j,l,jls_incloud(1,n))+
     *           dt_sulf_mc(n,l)
            tajls(j,l,jls_incloud(2,n))=tajls(j,l,jls_incloud(2,n))+
     *           dt_sulf_ss(n,l)
          if (ijts_aq(n).gt.0) then
           taijs(i,j,ijts_aq(n))=taijs(i,j,ijts_aq(n))+
     *           dt_sulf_mc(n,l)+dt_sulf_ss(n,l)
          end if
          end if
#endif
#endif
        end do
#ifdef TRACERS_WATER
        trprec(n,i,j) = trprec(n,i,j)+trprss(nx)
C**** diagnostics
        if (dowetdep(n)) then
          if (jls_prec(1,n).gt.0) tajls(j,1,jls_prec(1,n))=tajls(j,1
     *         ,jls_prec(1,n))+trprec(n,i,j)*bydxyp(j)
          if (jls_prec(2,n).gt.0) tajls(j,1,jls_prec(2,n))=tajls(j,1
     *         ,jls_prec(2,n))+trprec(n,i,j)*focean(i,j)*bydxyp(j)
          taijn(i,j,tij_prec,n) =taijn(i,j,tij_prec,n) +
     *         trprec(n,i,j)*bydxyp(j)
        end if
#endif
      end do
#endif

      END DO
C**** END OF MAIN LOOP FOR INDEX I

C****
Cred*           Reduced Arrays 3
C****
         DO L=1,LM
            WM(:,J,L) = WMIL(:,L)
         END DO
         DO L=1,LM
            T3MOM(:,:,J,L) = TMOMIL(:,:,L)
         END DO
         DO L=1,LM
            Q3MOM(:,:,J,L) = QMOMIL(:,:,L)
         END DO
Cred*       end Reduced Arrays 3
      END DO
C**** END OF MAIN LOOP FOR INDEX J
!$OMP  END PARALLEL DO
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

#ifdef TRACERS_ON
C**** Save the conservation quantities for tracers
      do nx=1,ntx
        n=ntix(nx)
        if (itcon_mc(n).gt.0) call diagtcb(dtr_mc(1,nx),itcon_mc(n),n)
        if (itcon_ss(n).gt.0) call diagtcb(dtr_ss(1,nx),itcon_ss(n),n)
      end do
#endif

!ESMF-GLOBAL SUM accumulated into AIL!!!
! note that sum is only over l<lmc(2,i,j) with each addition inside an if block
C**** Delayed summations (to control order of summands)
!     DO J=MAX(J_0,J5S),MIN(J_1,J5N)
      DO J=J_0,J_1
      DO I=1,IM
        IF(LMC(1,I,J).GT.0) THEN
          DO L=1,LMC(2,I,J)-1
            AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+AJEQIL(J,I,L)
          END DO
        END IF
      END DO
      END DO
C****EXCEPTION: partial-global sum above!
C     CALL GLOBAL_SUM(AIL.....)
C 
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
         JR=JREG(I,J)
         IF(LMC(1,I,J).GT.0)
     *     AREG(JR,J_PRCPMC)=AREG(JR,J_PRCPMC)+AREGIJ(I,J,1)
         AREG(JR,J_PRCPSS)=AREG(JR,J_PRCPSS)+AREGIJ(I,J,2)
         AREG(JR,J_EPRCP) =AREG(JR,J_EPRCP) +AREGIJ(I,J,3)
      END DO
      END DO
C 
C     NOW REALLY UPDATE THE MODEL WINDS
C 
CAOO      J=1
      IF(HAVE_SOUTH_POLE) THEN
        DO K=1,IM ! KMAXJ(J)
          IDI(K)=IDIJ(K,1,1)
          IDJ(K)=IDJJ(K,1)
        END DO
        DO L=1,LM
          DO K=1,IM ! KMAXJ(J)
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKP1(K,L)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKP1(K,L)
          END DO
        END DO
      ENDIF
C 
C**** Initialize dummy work arrays

!$OMP  PARALLEL DO PRIVATE(I,J,K,L,IDI,IDJ)
      DO L=1,LM
        DO J=J_0S,J_1-1
           DO K=1,4  !  KMAXJ(J)
             IDJ(K)=IDJJ(K,J)
           END DO
           DO I=1,IM
             DO K=1,4 ! KMAXJ(J)
               IDI(K)=IDIJ(K,I,J)
               U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)
               V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)
             END DO
           END DO
        END DO
C**** First half of loop cycle for j=j_1 for internal blocks
        IF (.not. HAVE_NORTH_POLE) then
          J=J_1
          DO K=1,2  !  KMAXJ(J)
            IDJ(K)=IDJJ(K,J)
          END DO
          DO I=1,IM
            DO K=1,2 ! KMAXJ(J)
              IDI(K)=IDIJ(K,I,J)
              U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)
              V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)
            END DO
          END DO
        ENDIF

C**** Second half of southern neighbor's j=j_1 cycle (equivalent to j=j_0-1
C**** in this block).

        if (.not. HAVE_SOUTH_POLE) then
          CALL CHECKSUM_COLUMN(grid, UKM, __LINE__, __FILE__)
          CALL CHECKSUM_COLUMN(grid, VKM, __LINE__, __FILE__)

          CALL HALO_UPDATE_COLUMN(grid, UKM, from=SOUTH)
          CALL HALO_UPDATE_COLUMN(grid, VKM, from=SOUTH)

C**** ....then, accumulate neighbors contribution to
C**** U,V, at the J=J_0 (B-grid) corners --i.e.do a
C**** K=3,4 iterations on the newly updated J=J_0-1 box.
          J=J_0-1
          DO K=3,4  !  KMAXJ(J)
            IDJ(K)=IDJJ(K,J)
          END DO
          DO I=1,IM
            DO K=3,4 ! KMAXJ(J)
              IDI(K)=IDIJ(K,I,J)
              U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)
              V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)
            END DO
          END DO
        end if           !NOT SOUTH POLE
      END DO       !LM
!$OMP  END PARALLEL DO
C 
CAOO      J=JM
      IF(HAVE_NORTH_POLE) THEN
        DO K=1,IM  !  KMAXJ(J)
          IDI(K)=IDIJ(K,1,JM)
          IDJ(K)=IDJJ(K,JM)
        END DO
        DO L=1,LM
          DO K=1,IM  !  KMAXJ(J)
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKPJM(K,L)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKPJM(K,L)
          END DO
        END DO
      END IF

C
C**** ADD IN CHANGE OF MOMENTUM BY MOIST CONVECTION AND CTEI
C**** and save changes in KE for addition as heat later
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM
      DO J=J_0STG,J_1STG
      DO I=1,IM
         AJL(J,L,JL_DAMMC)=AJL(J,L,JL_DAMMC)+
     &         (U(I,J,L)-UC(I,J,L))*PDSIG(L,I,J)
         DKE(I,J,L)=0.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L)
     *       -UC(I,J,L)*UC(I,J,L)-VC(I,J,L)*VC(I,J,L))
      END DO
      END DO
      END DO
!$OMP  END PARALLEL DO

      if (isccp_diags.eq.1) CALL RINIT(seed) ! reset random number sequ.

      RETURN
      END SUBROUTINE CONDSE

      SUBROUTINE init_CLD
!@sum  init_CLD initialises parameters for MSTCNV and LSCOND
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE CONSTANT, only : grav,by3
      USE MODEL_COM, only : dtsrc,ls1,sige,lm,psfmpt,ptop,plbot,jm
      USE DOMAIN_DECOMP, only : GRID
      USE GEOM, only : lat_dg
      USE CLOUDS, only : lmcm,bydtsrc,xmass,brcld,bybr,U00wtrX,U00ice
     *  ,HRMAX,ISC,lp50,RICldX,RWCldOX,xRIcld,do_blU00
      USE CLOUDS_COM, only : llow,lmid,lhi
      USE DAGCOM, only : nisccp,isccp_reg,isccp_late
      USE PARAM

      IMPLICIT NONE
      REAL*8 PLE
      INTEGER L,J,n
      INTEGER :: J_0,J_1

      J_0 =GRID%J_STRT
      J_1 =GRID%J_STOP

      call sync_param( 'U00wtrX', U00wtrX )
      call sync_param( 'U00ice', U00ice )
      call sync_param( "LMCM", LMCM )
      call sync_param( "HRMAX", HRMAX )
      call sync_param( "RICldX", RICldX )
      xRIcld = .001d0*(RICldX-1.d0)
      call sync_param( "RWCldOX", RWCldOX )
      call sync_param( "ISC", ISC)
      call sync_param( "do_blU00", do_blU00)

      IF(LMCM.LT.0) LMCM = LS1-1
      call set_param( "LMCM", LMCM, 'o' )

      BYDTsrc=1./DTsrc
      XMASS=0.1d0*DTsrc*GRAV

      BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3

C**** SEARCH FOR THE 50 MB LEVEL
      LP50=LM
      DO L=LM-1,1,-1
        PLE=.25*(SIGE(L)+2.*SIGE(L+1)+SIGE(L+2))*PSFMPT+PTOP
        IF (PLE.LT.50.) LP50=L
      END DO
      write(6,*) "Maximum level for LSCOND calculations (50mb): ",LP50

C**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS
      DO L=1,LM
        LLOW=L
        IF (.5*(PLbot(L+1)+PLbot(L+2)).LT.750.) EXIT ! was 786. 4/16/97
      END DO
      DO L=LLOW+1,LM
        LMID=L
        IF (.5*(PLbot(L+1)+PLbot(L+2)).LT.430.) EXIT
      END DO
      LHI=LM
      IF (LMID+1.GT.LHI) LHI=LMID+1
      WRITE (6,47) LLOW,LLOW+1,LMID,LMID+1,LHI
 47   FORMAT (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN',
     *     ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)

C**** Define regions for ISCCP diagnostics
      do j=J_0,J_1
        isccp_reg(j)=0
        do n=1,nisccp
           if(lat_dg(j,1).ge.isccp_late(n) .and.
     &        lat_dg(j,1).lt.isccp_late(n+1)) then
              isccp_reg(j)=n
              exit
           endif
        enddo
      end do

      END SUBROUTINE init_CLD
