#include "rundeck_opts.h"

      SUBROUTINE CONDSE
!@sum   CONDSE driver for moist convection AND large-scale condensation
!@auth  M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver   1.0 (taken from CB265)
!@calls CLOUDS:MSTCNV,CLOUDS:LSCOND

      USE CONSTANT, only : bygrav,lhm,rgas,grav,tf,lhe,lhs
      USE MODEL_COM, only : im,jm,lm,p,u,v,t,q,wm,JHOUR,fearth
     *     ,ls1,psf,ptop,dsig,bydsig,jeq,fland,sig,DTsrc,ftype
     *     ,ntype,itime,fim,airx,lmc,focean,fland,flice
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : bydxyp,dxyp,imaxj,kmaxj,ravj,idij,idjj
      USE CLOUDS_COM, only : ttold,qtold,svlhx,svlat,rhsav,cldsav
     *     ,pbltop,tauss,taumc,cldss,cldmc,csizmc,csizss
      USE CLOUDS, only : BYDTsrc,mstcnv,lscond ! glb var & subroutines
     *     ,airm,byam,etal,sm,smomij=>smom,qm,qmomij=>qmom
     *     ,tl,ris,ri1,ri2,aj8,aj11,aj13,aj50,aj51,aj52,aj53,aj57
     *     ,wml,sdl,u_0,v_0,um,vm,tf,qs,us,vs,dcl,airxl,prcpss,hcndss
     *     ,prcpmc,pearth,ts,taumcl,cldmcl,svwmxl,svlatl,svlhxl
     *     ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum
     *     ,aq,dpdt,th,ql,wmx,ttoldl,rh,taussl,cldssl,cldsavl,rh1
     *     ,kmax,ra,pl,ple,plk,rndss1l,rndss2l
#ifdef TRACERS_ON
#ifdef TRACERS_WATER
     *     ,trwml,trsvwml,trprmc,trprss
#endif
     *     ,tm,trmomij=>tmom      ! local  (i,j)
     *     ,ntx,ntix              ! global (same for all i,j)
      USE TRACER_DIAG_COM,only: tajln,jlnt_mc,jlnt_lscond,itcon_mc
     *     ,itcon_ss
#ifdef TRACERS_WATER
     *     ,jls_source,taijn,tajls,tij_prec
#endif
      USE TRACER_COM, only: itime_tr0,TRM,TRMOM,NTM
#ifdef TRACERS_WATER
     *     ,trwm
#endif
#endif
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg,tgvavg,qgavg,dclev
      USE DAGCOM, only : aj,areg,aij,ajl,ail,adiurn,jreg,ij_pscld,
     *     ij_pdcld,ij_scnvfrq,ij_dcnvfrq,ij_wmsum,ij_snwf,ij_prec,
     *     ij_neth,ij_f0oc,j_eprcp,j_prcpmc,j_prcpss,il_mceq,j5s,j5n, 
     *     ijdd,idd_pr,idd_ecnd,idd_mcp,idd_dmc,idd_smc,idd_ssp,
     &     jl_mcmflx,jl_sshr,jl_mchr,jl_dammc,jl_rhe,
     &     jl_mchphas,jl_mcdtotw,jl_mcldht,jl_mcheat,jl_mcdry,
     *     ij_ctpi,ij_taui,ij_lcldi,ij_mcldi,ij_hcldi,ij_tcldi,
     *     isccp_diags
      USE DYNAMICS, only : pk,pek,pmid,pedn,sd_clouds,gz,ptold,pdsig
     *     ,ltropo
      USE SEAICE_COM, only : rsi
      USE GHYCOM, only : snoage
      USE LAKES_COM, only : flake
      USE FLUXES, only : prec,eprec,precss,gtemp
#ifdef TRACERS_WATER
     *     ,trprec
#endif
      USE RANDOM
      USE QUSDEF, only : nmom
      IMPLICIT NONE

#ifdef TRACERS_ON
!@var tmsave holds tracer value (for diagnostics)
      REAL*8 tmsave(lm,ntm),dtr_mc(jm,ntm),dtr_ss(jm,ntm)
      INTEGER NX,N
#endif

!@var UC,VC,UZM,VZM velocity work arrays
      REAL*8, DIMENSION(IM,JM,LM) :: UC,VC
      REAL*8, DIMENSION(2,LM) :: UZM,VZM

!@param ENTCON fractional rate of entrainment (km**-1)
      REAL*8,  PARAMETER :: ENTCON = .2d0
      REAL*8,  PARAMETER :: DELTX=.608d0

      INTEGER I,J,K,L  !@var I,J,K,L loop variables
      INTEGER JR,KR,ITYPE,IT,IH
!@var JR = JREG(I,J)
!@var KR index for regional diagnostics
!@var ITYPE index for snow age
!@var IT index for surface types
!@var LERR,IERR error reporting
      INTEGER :: LERR, IERR
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID

      REAL*8 :: HCNDMC,PRCP,TPRCP,EPRCP,ENRGP,WMERR,ALPHA1,ALPHA2,ALPHAS
      REAL*8 :: DTDZ,DTDZS,DUDZ,DVDZ,DUDZS,DVDZS,THSV,THV1,THV2,QG,TGV
      REAL*8 :: DH1S,BYDH1S,DH12,BYDH12,DTDZG,DUDZG,DVDZG
!@var HCNDMC heating due to moist convection
!@var PRCP precipipation
!@var TPRCP temperature of precip  (deg. C)
!@var EPRCP sensible heat of precip
!@var ENRGP energy of precip
!@var WMERR DH12,BYDH12,DH1S,BYDH1S dummy variable
!@var THSV,THV1,THV2 vertual potential temperatures
!@var QG,TGV ground humidity,virt.temperature from pbl
!@var ALPHA1,ALPHA2,ALPHAS dummy variables
!@var DTDZ,DTDZS,DTDZG vertical potential temperature gradients
!@var DUDZ,DVDZ,DUDZS,DVDZS,DUDZG,DVDZG vertical wind gradients

C**** parameters and variables for isccp diags
      integer, parameter :: ntau=7,npres=7
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
      INTEGER ICKERR, JCKERR, JERR, seed
      REAL*8  RNDSS1(LM,IM,JM), RNDSS2(LM-1,IM,JM),xx
      REAL*8  AJEQIL(J5N-J5S+1,IM,JM), AREGIJ(IM,JM,3)
      REAL*8  UKP1(IM,LM), VKP1(IM,LM), UKPJM(IM,LM),VKPJM(IM,LM)
      REAL*8  UKM(4,IM,2:JM-1,LM), VKM(4,IM,2:JM-1,LM)

C
C     OBTAIN RANDOM NUMBERS FOR PARALLEL REGION
C
      DO J=1,JM
      DO I=1,IMAXJ(J)
        RNDSS1(LM,I,J)  = RANDU(xx)
        DO L=LM-1,1,-1
          RNDSS1(L,I,J) = RANDU(xx)
          RNDSS2(L,I,J) = RANDU(xx)
        END DO
C     Do not bother to save random numbers for isccp_clouds
      END DO
      END DO
C     But save the current seed in case isccp_routine is activated
      if (isccp_diags.eq.1) CALL RFINAL(seed)
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
        DO I=1,IM
          UZM(1,L)=UZM(1,L)+UC(I,2,L)
          UZM(2,L)=UZM(2,L)+UC(I,JM,L)
          VZM(1,L)=VZM(1,L)+VC(I,2,L)
          VZM(2,L)=VZM(2,L)+VC(I,JM,L)
        ENDDO
        UZM(1,L)=UZM(1,L)/FIM
        UZM(2,L)=UZM(2,L)/FIM
        VZM(1,L)=VZM(1,L)/FIM
        VZM(2,L)=VZM(2,L)/FIM
      ENDDO
      IH=JHOUR+1
#ifdef TRACERS_ON
C**** Find the ntx active tracers ntix(1->ntx)
      nx = 0
      do n=1,ntm
        if (itime.lt.itime_tr0(n)) cycle
        nx = nx+1
        ntix(nx) = n
      end do
      ntx = nx
#endif
C****
C**** MAIN J LOOP
C****
       ICKERR=0
       JCKERR=0
C$OMP  PARALLEL DO PRIVATE (
#ifdef TRACERS_ON
C$OMP*  NX,N,tmsave,
#endif
C$OMP*  ALPHAS,ALPHA1,ALPHA2,AT,BYDH1S,BYDH12, CC,CONV,CTP,
C$OMP*  DH1S,DH12,DTDZ,DTDZG,DTDZS,DUDZ,DUDZG,DUDZS,DVDZ,DVDZG,DVDZS,
C$OMP*  DTAU_S,DTAU_C,DEM_S,DEM_C, FQ_ISCCP, ENRGP,EPRCP,
C$OMP*  HCNDMC, I,ITYPE,IT,ITAU, IDI,IDJ,
C$OMP*  ITROP,IERR, J,JERR, K,KR, L,LERR, NBOX, PRCP,PFULL,PHALF,
C$OMP*  GZIL, SD_CLDIL, WMIL, TMOMIL, QMOMIL,        ! reduced arrays
C$OMP*  QG,QV, SKT, TGV,TPRCP,THSV,THV1,THV2,TAUOPT, WMERR)
C$OMP*    SCHEDULE(DYNAMIC,2)
C
      DO J=1,JM
C
Cred* Reduced Arrays 2
C
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
         TMOMIL(:,:,L) = TMOM(:,:,J,L)
      END DO
      DO L=1,LM
         QMOMIL(:,:,L) = QMOM(:,:,J,L)
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
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS
C****
      PEARTH=FEARTH(I,J)
      TS=TSAVG(I,J)
      QS=QSAVG(I,J)
      US=USAVG(I,J)
      VS=VSAVG(I,J)
      TGV=TGVAVG(I,J)
      QG=QGAVG(I,J)
      DCL=NINT(DCLEV(I,J))

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
      DPDT(1:LS1-1)=SIG(1:LS1-1)*(P(I,J)-PTOLD(I,J))*BYDTsrc
      DPDT(LS1:LM)=0.
      DO L=1,LM
C**** TEMPERATURES
        SM(L)  =T(I,J,L)*AIRM(L)
Cred    SMOMIJ(:,L) =TMOM(:,I,J,L)*AIRM(L)
        SMOMIJ(:,L) =TMOMIL(:,I,L)*AIRM(L)
        TL(L)=T(I,J,L)*PLK(L)
C**** MOISTURE (SPECIFIC HUMIDITY)
        QM(L)  =Q(I,J,L)*AIRM(L)
Cred    QMOMIJ(:,L) =QMOM(:,I,J,L)*AIRM(L)
        QMOMIJ(:,L) =QMOMIL(:,I,L)*AIRM(L)
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
        trmomij(:,l,nx) = trmom(:,i,j,l,ntix(nx))
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

C**** SET PRECIPITATION AND LATENT HEAT
      PREC(I,J)=0.
      TPRCP=T(I,J,1)*PK(1,I,J)-TF
#ifdef TRACERS_WATER
      TRPREC(:,I,J) = 0.
#endif

C**** SET DEFAULT FOR AIR MASS FLUX (STRAT MODEL)
      AIRX(I,J)=0.

C**** MOIST CONVECTION
      CALL MSTCNV(IERR,LERR)

C**** Error reports
      if (ierr.gt.0) then
        write(6,*) "Error in moist conv: i,j,l=",i,j,lerr
ccc     if (ierr.eq.2) stop "Subsid error: abs(c) > 1"
        if (ierr.eq.2) ickerr = 1
      end if

C**** ACCUMULATE MOIST CONVECTION DIAGNOSTICS
      IF (LMCMIN.GT.0) THEN
        AIJ(I,J,IJ_PSCLD)=AIJ(I,J,IJ_PSCLD)+CLDSLWIJ
        AIJ(I,J,IJ_PDCLD)=AIJ(I,J,IJ_PDCLD)+CLDDEPIJ
        IF(CLDSLWIJ.GT.1e-6) AIJ(I,J,IJ_SCNVFRQ)=AIJ(I,J,IJ_SCNVFRQ)+1.
        IF(CLDDEPIJ.GT.1e-6) AIJ(I,J,IJ_DCNVFRQ)=AIJ(I,J,IJ_DCNVFRQ)+1.
        AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
        HCNDMC=0.
        DO L=1,LMCMAX
          HCNDMC=HCNDMC+AJ13(L)+AJ50(L)
          AJL(J,L,JL_MCHR)=AJL(J,L,JL_MCHR)+AJ13(L)*BYDSIG(L)
          AJL(J,L,JL_MCHPHAS)=AJL(J,L,JL_MCHPHAS)+AJ50(L)*BYDSIG(L)
          AJL(J,L,JL_MCDTOTW)=AJL(J,L,JL_MCDTOTW)+AJ51(L)*BYDSIG(L)
CCC       IF(J.GE.J5S.AND.J.LE.J5N) AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+
CCC  *         (AJ13(L)+AJ50(L))*(DXYP(J)*BYDSIG(L))
          IF(J.GE.J5S.AND.J.LE.J5N)     ! add in after parallel region
     *      AJEQIL(J-J5S+1,I,L) = (AJ13(L)+AJ50(L))*(DXYP(J)*BYDSIG(L))
          AJL(J,L,JL_MCHEAT)=AJL(J,L,JL_MCHEAT)+
     &         (AJ50(L)+AJ13(L))*BYDSIG(L)
          AJL(J,L,JL_MCDRY)=AJL(J,L,JL_MCDRY)+
     &         (AJ52(L)-AJ57(L))*BYDSIG(L)
          AJL(J,L,JL_MCMFLX)=AJL(J,L,JL_MCMFLX)+AJ8(L)
        END DO
        DO IT=1,NTYPE
          AJ(J,J_PRCPMC,IT)=AJ(J,J_PRCPMC,IT)+PRCPMC*FTYPE(IT,I,J)
        END DO
CCC     AREG(JR,J_PRCPMC)=AREG(JR,J_PRCPMC)+PRCPMC*DXYP(J)
        AREGIJ(I,J,1)=PRCPMC*DXYP(J)  ! add in after parallel region
        DO KR=1,4
          IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
            ADIURN(IH,IDD_PR  ,KR)=ADIURN(IH,IDD_PR  ,KR)+PRCPMC
            ADIURN(IH,IDD_ECND,KR)=ADIURN(IH,IDD_ECND,KR)+HCNDMC
            ADIURN(IH,IDD_MCP ,KR)=ADIURN(IH,IDD_MCP ,KR)+PRCPMC
            ADIURN(IH,IDD_DMC ,KR)=ADIURN(IH,IDD_DMC ,KR)+CLDDEPIJ
            ADIURN(IH,IDD_SMC ,KR)=ADIURN(IH,IDD_SMC ,KR)+CLDSLWIJ
          END IF
        END DO

C**** WRITE TO SOME GLOBAL ARRAYS
        PREC(I,J)=PRCPMC*100.*BYGRAV
        DO L=1,LMCMAX
          T(I,J,L)=  SM(L)*BYAM(L)
          Q(I,J,L)=  QM(L)*BYAM(L)
        END DO
        CSIZMC(1:LMCMAX,I,J)=CSIZEL(1:LMCMAX)
        AIRX(I,J) = AIRXL*DXYP(J)
      END IF                    ! should this be after tracers....????
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
      RNDSS1L(:)=RNDSS1(:,I,J)
      RNDSS2L(:)=RNDSS2(:,I,J)

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
        DH12=(GZ(I,J,2)-GZ(I,J,1))*BYGRAV
        BYDH12=1./DH12
        DTDZS=(THV1-THSV)*BYDH1S
        DTDZ=(THV2-THV1)*BYDH12
        IF (J.EQ.1) THEN
          DUDZ=(UZM(1,2)-UZM(1,1))*BYDH12
          DVDZ=(VZM(1,2)-VZM(1,1))*BYDH12
          DUDZS=(UZM(1,1)-US)*BYDH1S
          DVDZS=(VZM(1,1)-VS)*BYDH1S
        ENDIF
        IF (J.EQ.JM) THEN
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
        RIS=(GRAV*ALPHAS*DTDZG)/(DUDZG*DUDZG+DVDZG*DVDZG)
        RI1=(GRAV*ALPHA1*DTDZS)/(DUDZS*DUDZS+DVDZS*DVDZS)
        RI2=(GRAV*ALPHA2*DTDZ)/(DUDZ*DUDZ+DVDZ*DVDZ)
C       WRITE (6,*)'I,J,QG,TGV,THSV,RIS,RI1=',I,J,QG,TGV,THSV,RIS,RI1
      ENDIF
C**** LARGE-SCALE CLOUDS AND PRECIPITATION

      CALL LSCOND(IERR,WMERR,LERR)

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
         DO KR=1,4
           IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             ADIURN(IH,IDD_PR  ,KR)=ADIURN(IH,IDD_PR  ,KR)+PRCPSS
             ADIURN(IH,IDD_ECND,KR)=ADIURN(IH,IDD_ECND,KR)+HCNDSS
             ADIURN(IH,IDD_SSP ,KR)=ADIURN(IH,IDD_SSP ,KR)+PRCPSS
           END IF
         END DO

C**** TOTAL PRECIPITATION AND AGE OF SNOW
      PREC(I,J)=PREC(I,J)+PRCPSS*100.*BYGRAV
      PRCP=PREC(I,J)
      PRECSS(I,J)=PRCPSS*100.*BYGRAV
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
      IF (TPRCP.GE.0.) THEN
C       EPRCP=PRCP*TPRCP*SHW
        EPRCP=0.
        ENRGP=EPRCP
      ELSE
C       EPRCP=PRCP*TPRCP*SHI
        EPRCP=0.
        ENRGP=EPRCP-PRCP*LHM
        AIJ(I,J,IJ_SNWF)=AIJ(I,J,IJ_SNWF)+PRCP
      END IF
      EPREC(I,J)=ENRGP  ! energy of precipitation
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

      IF(TPRCP.LT.0.) THEN ! MODIFY SNOW AGES AFTER SNOW FALL
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
        if(jerr.ne.0) jckerr = 1
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
      end if

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
      TTOLD(:,I,J)=TH(:)
      QTOLD(:,I,J)=QL(:)

      DO L=1,LM
        AJL(J,L,JL_SSHR)=AJL(J,L,JL_SSHR)+AJ11(L)
        AJL(J,L,JL_MCLDHT)=AJL(J,L,JL_MCLDHT)+AJ53(L)
        AJL(J,L,JL_RHE)=AJL(J,L,JL_RHE)+RH1(L)

        T(I,J,L)=TH(L)
        Q(I,J,L)=QL(L)
C**** update moment changes
Cred    TMOM(:,I,J,L)=SMOMIJ(:,L)*BYAM(L)
Cred    QMOM(:,I,J,L)=QMOMIJ(:,L)*BYAM(L)
Cred    WM(I,J,L)=WMX(L)
        TMOMIL(:,I,L)=SMOMIJ(:,L)*BYAM(L)
        QMOMIL(:,I,L)=QMOMIJ(:,L)*BYAM(L)
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

#ifdef TRACERS_ON
C**** TRACERS: Use only the active ones
      do nx=1,ntx
        n = ntix(nx)
        do l=1,lm
          dtr_ss(j,nx)=dtr_ss(j,nx)+(tm(l,nx)-tmsave(l,nx))
#ifdef TRACERS_WATER
     &         + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
#endif
          trm(i,j,l,n) = tm(l,nx)
          trmom(:,i,j,l,n) = trmomij(:,l,nx)
          tajln(j,l,jlnt_lscond,n) = tajln(j,l,jlnt_lscond,n) +
     &          (tm(l,nx)-tmsave(l,nx))
#ifdef TRACERS_WATER
     &         + (trwml(nx,l)-trwm(i,j,l,n)-trsvwml(nx,l))
          trwm(i,j,l,n) = trwml(nx,l)
#endif
        end do
#ifdef TRACERS_WATER
        trprec(n,i,j) = trprec(n,i,j)+trprss(nx)
C**** diagnostics
        tajls(j,1,jls_source(3,n))=tajls(j,1,jls_source(3,n))
     *       +trprec(n,i,j)*bydxyp(j)
        tajls(j,1,jls_source(4,n))=tajls(j,1,jls_source(4,n))
     *       +trprec(n,i,j)*focean(i,j)*bydxyp(j)
        taijn(i,j,tij_prec,n) =taijn(i,j,tij_prec,n)+trprec(n,i,j)
     *       *bydxyp(j)
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
            TMOM(:,:,J,L) = TMOMIL(:,:,L)
         END DO
         DO L=1,LM
            QMOM(:,:,J,L) = QMOMIL(:,:,L)
         END DO
Cred*       end Reduced Arrays 3
      END DO
C**** END OF MAIN LOOP FOR INDEX J
C$OMP  END PARALLEL DO
C****
C
C     WAS THERE AN ERROR IN SUBSID ??
C
      IF(ICKERR.NE.0)  THEN
         WRITE(6,*)  'SUBSID ERROR: ABS(C) > 1'
         STOP 800
      END IF
C
C     WAS THERE AN ERROR IN ISCCP CLOUD TYPING ??
C
      IF(JCKERR.NE.0)  THEN
         WRITE(6,*)  'ISCCP CLOUD TYPING ERROR'
         STOP 900
      END IF

#ifdef TRACERS_ON
C**** Save the conservation quantities for tracers
      do nx=1,ntx
        n=ntix(nx)
        call diagtcb(dtr_mc(1,nx),itcon_mc(n),n)
        call diagtcb(dtr_ss(1,nx),itcon_ss(n),n)
      end do
#endif

C**** Delayed summations (to control order of summands)
      DO J=J5S,J5N
      DO I=1,IM
        IF(LMC(1,I,J).GT.0) THEN
          DO L=1,LMC(2,I,J)-1
            AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+AJEQIL(J-J5S+1,I,L)
          END DO
        END IF
      END DO
      END DO
C
      DO J=1,JM
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
      J=1
      DO K=1,IM ! KMAXJ(J)
         IDI(K)=IDIJ(K,1,J)
         IDJ(K)=IDJJ(K,J)
      END DO
      DO L=1,LM
      DO K=1,IM ! KMAXJ(J)
         U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKP1(K,L)
         V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKP1(K,L)
      END DO
      END DO
C
C$OMP  PARALLEL DO PRIVATE(I,J,K,L,IDI,IDJ)
      DO L=1,LM
      DO J=2,JM-1
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
      END DO
C$OMP  END PARALLEL DO
C
      J=JM
      DO K=1,IM  !  KMAXJ(J)
         IDI(K)=IDIJ(K,1,J)
         IDJ(K)=IDJJ(K,J)
      END DO
      DO L=1,LM
      DO K=1,IM  !  KMAXJ(J)
         U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKPJM(K,L)
         V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKPJM(K,L)
      END DO
      END DO
C
C**** ADD IN CHANGE OF MOMENTUM BY MOIST CONVECTION AND CTEI
C$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO L=1,LM
      DO J=2,JM
      DO I=1,IM
         AJL(J,L,JL_DAMMC)=AJL(J,L,JL_DAMMC)+
     &         (U(I,J,L)-UC(I,J,L))*PDSIG(L,I,J)
      END DO
      END DO
      END DO
C$OMP  END PARALLEL DO

      if (isccp_diags.eq.1) CALL RINIT(seed) ! reset random number sequ.

      RETURN
      END SUBROUTINE CONDSE

      SUBROUTINE init_CLD
!@sum  init_CLD initialises parameters for MSTCNV and LSCOND
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE CONSTANT, only : grav,by3
      USE MODEL_COM, only : dtsrc,ls1
      USE CLOUDS, only : lmcm,bydtsrc,xmass,brcld,bybr,U00wtr,U00ice
     *  ,HRMAX
      USE PARAM

      IMPLICIT NONE

      call sync_param( 'U00wtr', U00wtr )
      call sync_param( 'U00ice', U00ice )
      call sync_param( "LMCM", LMCM )
      call sync_param( "HRMAX", HRMAX )

      IF(LMCM.LT.0) LMCM = LS1-1
      call set_param( "LMCM", LMCM, 'o' )

      BYDTsrc=1./DTsrc
      XMASS=0.1d0*DTsrc*GRAV

      BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3

      END SUBROUTINE init_CLD
