      SUBROUTINE CONDSE
!@sum   CONDSE driver for moist convection AND large-scale condensation
!@auth  M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver   1.0 (taken from CB265)
!@calls CLOUDS:MSTCNV,CLOUDS:LSCOND

      USE CONSTANT, only : bygrav,lhm,rgas,grav,tf,lhe
      USE MODEL_COM, only : im,jm,lm,p,u,v,t,q,wm,JHOUR,fearth
     *     ,ls1,psf,ptop,dsig,bydsig,jeq,fland,sig,DTsrc,ftype
     *     ,ntype,itime,fim,airx,lmc,focean,fland,flice
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : bydxyp,dxyp,imaxj,kmaxj,ravj,idij,idjj
      USE CLOUDS_COM, only : ttold,qtold,svlhx,svlat,rhsav,cldsav
     *     ,pbltop,tauss,taumc,cldss,cldmc,csizmc,csizss
      USE CLOUDS, only : kmax,ra,pl,ple,plk
     *     ,airm,byam,etal,sm,smomij=>smom,qm,qmomij=>qmom
     *     ,tl,ri1,ri2,aj13
     *     ,aj50,aj51,aj52,aj57,aj8,aj11,wml,sdl,u_0,v_0,um,vm,tf
     *     ,prcpmc,pearth,ts,taumcl,cldmcl,svwmxl,svlatl,svlhxl
     *     ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum
     *     ,mstcnv,qs,us,vs,dcl
     *     ,aq,dpdt,th,ql,wmx,ttoldl,rh,lpbl,taussl,cldssl,cldsavl,
     *     prcpss,hcndss,aj53,BYDTsrc,lscond,airxl
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg,dclev
      USE DAGCOM, only : aj,areg,aij,ajl,ail,adiurn,jreg,ij_pscld,
     *     ij_pdcld,ij_scnvfrq,ij_dcnvfrq,ij_wmsum,ij_snwf,ij_prec,
     *     ij_neth,j_eprcp,j_prcpmc,j_prcpss,il_mceq,j5s,j5n,
     *     ijdd,idd_pr,idd_ecnd,idd_mcp,idd_dmc,idd_smc,idd_ssp,
     &     jl_mcmflx,jl_sshr,jl_mchr,jl_dammc,
     &     jl_mchphas,jl_mcdtotw,jl_mcdlht,jl_mcheat,jl_mcdry,
     *     ij_ctpi,ij_taui,ij_lcldi,ij_mcldi,ij_hcldi,ij_tcldi,
     *     isccp_diags
      USE DYNAMICS, only : pk,pek,pmid,pedn,sd_clouds,gz,ptold,pdsig
     *     ,ltropo
      USE SEAICE_COM, only : rsi
      USE GHYCOM, only : snoage
      USE LAKES_COM, only : flake
      USE FLUXES, only : prec,eprec,precss,gtemp

      IMPLICIT NONE

!@var UC,VC,UZM,VZM velocity work arrays
      REAL*8, DIMENSION(IM,JM,LM) :: UC,VC
      REAL*8, DIMENSION(2,LM) :: UZM,VZM

!@param ENTCON fractional rate of entrainment (km**-1)
      REAL*8,  PARAMETER :: ENTCON = .2d0
      REAL*8,  PARAMETER :: DELTX=.608d0

      INTEGER I,J,K,L  !@var I,J,K,L loop variables
      INTEGER IMAX,JR,KR,ITYPE,IT,IM1,IH
!@var IMAX maximum number of zonal grid points used
!@var IM1 IM-1
!@var JR = JREG(I,J)
!@var KR index for regional diagnostics
!@var ITYPE index for snow age
!@var IT index for surface types
!@var LERR,IERR error reporting 
      INTEGER :: LERR=0, IERR=0
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID

      REAL*8 :: HCNDMC,PRCP,TPRCP,EPRCP,ENRGP,WMERR,ALPHA1,ALPHA2,THV1
      REAL*8 :: DH12,DTDZ,DTDZS,DUDZ,DVDZ,DUDZS,DVDZS,THSV,THV2
      REAL*8 :: DH1S,BYDH1S,BYDH12
!@var HCNDMC heating due to moist convection
!@var PRCP precipipation
!@var TPRCP temperature of precip  (deg. C)
!@var EPRCP sensible heat of precip
!@var ENRGP energy of precip
!@var WMERR DH12,BYDH12,DH1S,BYDH1S,THV1,THV2 dummy variables
!@var ALPHA1,ALPHA2,DTDZ,DTDZS,DUDZ,DVDZ,DUDZS,DVDZS dummy variables

C**** parameters and variables for isccp diags
      integer, parameter :: ntau=7,npres=7
      real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
      real*8 skt,conv(lm),qv(lm)
      real*8 pfull(lm),at(lm),cc(lm),dtau_s(lm),dtau_c(lm)
      real*8 dem_s(lm),dem_c(lm),phalf(lm+1)
      real*8 fq_isccp(ntau,npres),ctp,tauopt,sumcld
      integer itau,ipres,itrop,nbox
C**** 

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
        UZM(1,L)=UZM(1,L)+U(I,2,L)
        UZM(2,L)=UZM(2,L)+U(I,JM,L)
        VZM(1,L)=VZM(1,L)+V(I,2,L)
        VZM(2,L)=VZM(2,L)+V(I,JM,L)
       ENDDO
       UZM(1,L)=UZM(1,L)/FIM
       UZM(2,L)=UZM(2,L)/FIM
       VZM(1,L)=VZM(1,L)/FIM
       VZM(2,L)=VZM(2,L)/FIM
      ENDDO
      IH=JHOUR+1
C****
C**** MAIN J LOOP
C****
      DO J=1,JM

      KMAX=KMAXJ(J)
      IMAX=IMAXJ(J)
C****
C**** MAIN I LOOP
C****
      IM1=IM
      DO I=1,IMAX
         JR=JREG(I,J)
C****
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS
C****
      PEARTH=FEARTH(I,J)
      TS=TSAVG(I,J)
      QS=QSAVG(I,J)
      US=USAVG(I,J)
      VS=VSAVG(I,J)
      DCL=DCLEV(I,J)
      LPBL=1

      DO K=1,KMAX
         RA(K)=RAVJ(K,J)
         IDI(K)=IDIJ(K,I,J)
         IDJ(K)=IDJJ(K,J)
      END DO
C**** PRESSURES, AND PRESSURE TO THE KAPA
      DO 150 L=1,LM
      PL(L) =PMID(L,I,J)
      PLE(L)=PEDN(L,I,J)
      PLK(L)=PK(L,I,J)
      AIRM(L)=PDSIG(L,I,J)
      BYAM(L)=1./AIRM(L)
      IF(L.LE.LM-2) ETAL(L+1)=.5*ENTCON*(GZ(I,J,L+2)-GZ(I,J,L))*
     *     1.d-3*BYGRAV
C**** TEMPERATURES
      SM(L)  =T(I,J,L)*AIRM(L)
      SMOMIJ(:,L) =TMOM(:,I,J,L)*AIRM(L)
      TL(L)=T(I,J,L)*PLK(L)
C**** MOISTURE (SPECIFIC HUMIDITY)
      QM(L)  =Q(I,J,L)*AIRM(L)
      QMOMIJ(:,L) =QMOM(:,I,J,L)*AIRM(L)
c      QL(L)=Q(I,J,L)
      WML(L)=WM(I,J,L)
      SDL(L)=SD_CLOUDS(I,J,L)*BYDXYP(J)
      SVLHXL(L)=SVLHX(L,I,J)
         TTOLDL(L)=TTOLD(L,I,J)
         CLDSAVL(L)=CLDSAV(L,I,J)
         RH(L)=RHSAV(L,I,J)
         DPDT(L)=SIG(L)*(P(I,J)-PTOLD(I,J))*BYDTsrc
         IF(L.GE.LS1) DPDT(L)=0.
  150 CONTINUE
      ETAL(LM)=ETAL(LM-1)
      PLE(LM+1)=PEDN(LM+1,I,J)
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

C**** SET DEFAULT FOR AIR MASS FLUX (STRAT MODEL)
      AIRX(I,J)=0.

C**** MOIST CONVECTION
      CALL MSTCNV(IERR,LERR)

C**** Error reports
      if (ierr.gt.0) then  
        write(6,*) "Error in moist conv: i,j,l=",i,j,lerr
        if (ierr.eq.2) stop "Subsid error: abs(c) > 1"
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
            IF(J.GE.J5S.AND.J.LE.J5N) AIL(I,L,IL_MCEQ)=AIL(I,L,IL_MCEQ)+
     *           (AJ13(L)+AJ50(L))*(DXYP(J)*BYDSIG(L))
            AJL(J,L,JL_MCHEAT)=AJL(J,L,JL_MCHEAT)+
     &           (AJ50(L)+AJ13(L))*BYDSIG(L)
            AJL(J,L,JL_MCDRY)=AJL(J,L,JL_MCDRY)+
     &           (AJ52(L)-AJ57(L))*BYDSIG(L)
            AJL(J,L,JL_MCMFLX)=AJL(J,L,JL_MCMFLX)+AJ8(L)
         END DO
         DO IT=1,NTYPE
           AJ(J,J_PRCPMC,IT)=AJ(J,J_PRCPMC,IT)+PRCPMC*FTYPE(IT,I,J)
         END DO
         AREG(JR,J_PRCPMC)=AREG(JR,J_PRCPMC)+PRCPMC*DXYP(J)
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
            CSIZMC(L,I,J)=CSIZEL(L)
         END DO
         AIRX(I,J) = AIRXL*DXYP(J)
         LMC(1,I,J) = LMCMIN
         LMC(2,I,J) = LMCMAX+1
      END IF
C****
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS
C****
      DO L=1,LM
         TL(L)=T(I,J,L)*PLK(L)
         TH(L)=T(I,J,L)
         QL(L)=Q(I,J,L)
         WMX(L)=WML(L)+SVWMXL(L)
         AQ(L)=(QL(L)-QTOLD(L,I,J))*BYDTsrc
      END DO

C**** COMPUTE RICHARDSON NUMBER
      IF(DCL.LE.1) THEN
        THSV=TS*(1.+DELTX*QS)/PEK(1,I,J)
        THV1=TH(1)*(1.+DELTX*QL(1))
        THV2=TH(2)*(1.+DELTX*QL(2))
        ALPHA1=2./(TH(1)+TS/PEK(1,I,J))
        ALPHA2=2./(TH(1)+TH(2))
        DH1S=(PLE(1)-PL(1))*TL(1)*RGAS/(GRAV*PL(1))
        BYDH1S=1./DH1S
        DH12=(GZ(I,J,2)-GZ(I,J,1))/GRAV
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
          DUDZ=(U(IDI(1),IDJ(1),2)+U(IDI(2),IDJ(2),2)+
     *         U(IDI(3),IDJ(3),2)+U(IDI(4),IDJ(4),2)-
     *         U(IDI(1),IDJ(1),1)-U(IDI(2),IDJ(2),1)-
     *         U(IDI(3),IDJ(3),1)-U(IDI(4),IDJ(4),1))*.25*BYDH12
          DVDZ=(V(IDI(1),IDJ(1),2)+V(IDI(2),IDJ(2),2)+
     *         V(IDI(3),IDJ(3),2)+V(IDI(4),IDJ(4),2)-
     *         V(IDI(1),IDJ(1),1)-V(IDI(2),IDJ(2),1)-
     *         V(IDI(3),IDJ(3),1)-V(IDI(4),IDJ(4),1))*.25*BYDH12
          DUDZS=(U(IDI(1),IDJ(1),1)+U(IDI(2),IDJ(2),1)+
     *         U(IDI(3),IDJ(3),1)+U(IDI(4),IDJ(4),1)-
     *         4.*US)*.25*BYDH1S
          DVDZS=(V(IDI(1),IDJ(1),1)+V(IDI(2),IDJ(2),1)+
     *         V(IDI(3),IDJ(3),1)+V(IDI(4),IDJ(4),1)-
     *         4.*VS)*.25*BYDH1S
        ENDIF
        RI1=(GRAV*ALPHA1*DTDZS)/(DUDZS*DUDZS+DVDZS*DVDZS)
        RI2=(GRAV*ALPHA2*DTDZ)/(DUDZ*DUDZ+DVDZ*DVDZ)
C       WRITE (6,*)'I,J,TS,THSV,THV1,RI1,RI2=',I,J,TS,THSV,THV1,RI1,RI2
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
         AREG(JR,J_PRCPSS)=AREG(JR,J_PRCPSS)+PRCPSS*DXYP(J)
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
        AREG(JR,J_EPRCP)=AREG(JR,J_EPRCP)+ENRGP*DXYP(J)
        AIJ(I,J,IJ_PREC)=AIJ(I,J,IJ_PREC)+PRCP
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+ENRGP

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
          if(svlhxl(LM+1-L) .eq. lhe ) then   ! water cloud
            dem_s(l)=1.-exp(-taussl(LM+1-L)*bywc)  
            dem_c(l)=1.-exp(-taumcl(LM+1-L)*bywc)  
          else                                     ! ice cloud
            dem_s(l)=1.-exp(-taussl(LM+1-L)*byic)  
            dem_c(l)=1.-exp(-taumcl(LM+1-L)*byic)  
          endif
          
          qv(l)=ql(LM+1-L)  
        end do
        phalf(lm+1)=ple(1)*100.
        itrop = LM+1-LTROPO(I,J)

        call ISCCP_CLOUD_TYPES(pfull,phalf,qv,
     &       cc,conv,dtau_s,dtau_c,skt,
     &       at,dem_s,dem_c,itrop,fq_isccp,ctp,tauopt,nbox)
          
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
      DO L=1,LM
         TAUMC(L,I,J)=TAUMCL(L)
         CLDMC(L,I,J)=CLDMCL(L)
         SVLAT(L,I,J)=SVLATL(L)

         TAUSS(L,I,J)=TAUSSL(L)
         CLDSS(L,I,J)=CLDSSL(L)
         CLDSAV(L,I,J)=CLDSAVL(L)
         SVLHX(L,I,J)=SVLHXL(L)
         CSIZSS(L,I,J)=CSIZEL(L)
         AJL(J,L,JL_SSHR)=AJL(J,L,JL_SSHR)+AJ11(L)
         AJL(J,L,JL_MCDLHT)=AJL(J,L,JL_MCDLHT)+AJ53(L)

         T(I,J,L)=TH(L)
         Q(I,J,L)=QL(L)
C**** update moment changes
         TMOM(:,I,J,L)=SMOMIJ(:,L)*BYAM(L)
         QMOM(:,I,J,L)=QMOMIJ(:,L)*BYAM(L)
         RHSAV(L,I,J)=RH(L)
         TTOLD(L,I,J)=TH(L)
         QTOLD(L,I,J)=QL(L)
         WM(I,J,L)=WMX(L)

C**** UPDATE MODEL WINDS
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UM(K,L)*BYAM(L)-UC(IDI(K),IDJ(K),L))
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VM(K,L)*BYAM(L)-VC(IDI(K),IDJ(K),L))
         ENDDO
      ENDDO

      IM1=I
      END DO
C**** END OF MAIN LOOP FOR INDEX I
      END DO
C**** END OF MAIN LOOP FOR INDEX J
C****

C**** ADD IN CHANGE OF MOMENTUM BY MOIST CONVECTION AND CTEI
      DO L=1,LM
         DO J=2,JM
            DO I=1,IM
               AJL(J,L,JL_DAMMC)=AJL(J,L,JL_DAMMC)+
     &              (U(I,J,L)-UC(I,J,L))*PDSIG(L,I,J)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE CONDSE

      SUBROUTINE init_CLD
!@sum  init_CLD initialises parameters for MSTCNV and LSCOND
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE CONSTANT, only : grav,by3
      USE MODEL_COM, only : dtsrc,ls1
      USE CLOUDS, only : lmcm,bydtsrc,xmass,brcld,bybr,U00wtr,U00ice
      USE PARAM

      IMPLICIT NONE

      call sync_param( 'U00wtr', U00wtr )
      call sync_param( 'U00ice', U00ice )
      call sync_param( "LMCM", LMCM )

      IF(LMCM.LT.0) LMCM = LS1-1
      call set_param( "LMCM", LMCM, 'o' )

      BYDTsrc=1./DTsrc
      XMASS=0.1d0*DTsrc*GRAV

      BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3

      END SUBROUTINE init_CLD
