      SUBROUTINE CONDSE
!@sum   CONDSE driver for moist convection AND large-scale condensation
!@auth  M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver   1.0 (taken from CB265)
!@calls MSTCNV and LSCOND

      USE CONSTANT, only : bygrav,lhm
      USE E001M12_COM, only : im,jm,lm,p,u,v,t,q,wm,JHOUR,fearth
     *     ,ls1,psf,ptop,dsig,bydsig,jeq,fland,ijd6,gdata,sig,DTsrc
      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : bydxyp,dxyp,imaxj,kmaxj,raj,idij,idjj
      USE CLD01_COM_E001
      USE CLD01, only : kmax,ra,pl,ple,plk
     *     ,airm,byam,etal,sm,smomij=>smom,qm,qmomij=>qmom
     *     ,tl,aj13
     *     ,aj50,aj51,aj52,aj57,aj8,aj11,wml,sdl,u_0,v_0,um,vm,tf
     *     ,prcpmc,pearth,ts,taumcl,cldmcl,svwmxl,svlatl,svlhxl
     *     ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum
     *     ,mstcnv
     *     ,aq,dpdt,th,ql,wmx,ttoldl,rh,lpbl,taussl,cldssl,cldsavl,
     *     prcpss,hcndss,aj55,BYDTsrc,lscond
      USE PBLCOM, only : tsavg
      USE DAGCOM, only : aj,bj,cj,areg,aij,ajl,ail,adaily,jreg,ij_pscld
     *     ,ij_pdcld,ij_scnvfrq,ij_dcnvfrq,ij_wmsum,ij_snwf,ij_prec
     *     ,ij_neth,j_eprcp,j_prcpmc,j_prcpss
      USE DYNAMICS, only : pk,pmid,pedn,sd_clouds,gz,ptold,pdsig
      USE SEAICE_COM, only : rsi

      IMPLICIT NONE

!@var UC,VC velocity work arrays
      REAL*8, DIMENSION(IM,JM,LM) :: UC,VC

      REAL*8,  PARAMETER :: ENTCON = .2d0  !@param ENTCON  ???

      INTEGER I,J,K,L  !@var I,J,K,L loop variables
      INTEGER IMAX,JR,KR
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID

      REAL*8 :: HCNDMC,PRCP,TPRCP,EPRCP,ENRGP

C**** SAVE UC AND VC, AND ZERO OUT CLDSS AND CLDMC
      UC=U
      VC=V
C****
C**** MAIN J LOOP
C****
      DO J=1,JM

      KMAX=KMAXJ(J)
      IMAX=IMAXJ(J)
C****
C**** MAIN I LOOP
C****
      DO I=1,IMAX
         JR=JREG(I,J)
C****
C**** SET UP VERTICAL ARRAYS, OMITTING THE J AND I SUBSCRIPTS
C****
      PEARTH=FEARTH(I,J)
      TS=TSAVG(I,J)
      LPBL=1

      DO K=1,KMAX
         RA(K)=RAJ(K,J)
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
      TPREC(I,J)=T(I,J,1)*PK(1,I,J)-TF

C**** MOIST CONVECTION
      CALL MSTCNV

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
            AJL(J,L,13)=AJL(J,L,13)+AJ13(L)*BYDSIG(L)
            AJL(J,L,50)=AJL(J,L,50)+AJ50(L)*BYDSIG(L)
            AJL(J,L,51)=AJL(J,L,51)+AJ51(L)*BYDSIG(L)
            IF(J.GE.JEQ-2.AND.J.LE.JEQ) AIL(I,L,6)=AIL(I,L,6)+
     *           (AJ13(L)+AJ50(L))*(DXYP(J)*BYDSIG(L))
            AJL(J,L,56)=AJL(J,L,56)+(AJ50(L)+AJ13(L))*BYDSIG(L)
            AJL(J,L,57)=AJL(J,L,57)+(AJ52(L)-AJ57(L))*BYDSIG(L)
            AJL(J,L,8)=AJL(J,L,8)+AJ8(L)
         END DO
         AJ(J,J_PRCPMC)=AJ(J,J_PRCPMC)+PRCPMC*(1.-FLAND(I,J))*
     *        (1.-RSI(I,J))
         BJ(J,J_PRCPMC)=BJ(J,J_PRCPMC)+PRCPMC*FLAND(I,J)
         CJ(J,J_PRCPMC)=CJ(J,J_PRCPMC)+PRCPMC*RSI(I,J)*
     *        (1.-FLAND(I,J))
         AREG(JR,J_PRCPMC)=AREG(JR,J_PRCPMC)+PRCPMC*DXYP(J)
         DO KR=1,4
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(JHOUR+1,49,KR)=ADAILY(JHOUR+1,49,KR)+PRCPMC
               ADAILY(JHOUR+1,5,KR)=ADAILY(JHOUR+1,5,KR)+HCNDMC
               ADAILY(JHOUR+1,63,KR)=ADAILY(JHOUR+1,63,KR)+PRCPMC
               ADAILY(JHOUR+1,51,KR)=ADAILY(JHOUR+1,51,KR)+CLDDEPIJ
               ADAILY(JHOUR+1,52,KR)=ADAILY(JHOUR+1,52,KR)+CLDSLWIJ
            END IF
         END DO

C**** WRITE TO SOME GLOBAL ARRAYS
         PREC(I,J)=PRCPMC*100.*BYGRAV
         DO L=1,LMCMAX
            T(I,J,L)=  SM(L)*BYAM(L)
            Q(I,J,L)=  QM(L)*BYAM(L)
            CSIZMC(L,I,J)=CSIZEL(L)
         END DO
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

C**** LARGE-SCALE CLOUDS AND PRECIPITATION

      CALL LSCOND(I,J)

C**** Accumulate diagnostics of LSCOND
         AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
         AJ(J,J_PRCPSS)=AJ(J,J_PRCPSS)+PRCPSS*(1.-FLAND(I,J))*
     *        (1.-RSI(I,J))
         BJ(J,J_PRCPSS)=BJ(J,J_PRCPSS)+PRCPSS*FLAND(I,J)
         CJ(J,J_PRCPSS)=CJ(J,J_PRCPSS)+PRCPSS*RSI(I,J)*
     *        (1.-FLAND(I,J))
         AREG(JR,J_PRCPSS)=AREG(JR,J_PRCPSS)+PRCPSS*DXYP(J)
         DO KR=1,4
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(JHOUR+1,49,KR)=ADAILY(JHOUR+1,49,KR)+PRCPSS
               ADAILY(JHOUR+1,5,KR)=ADAILY(JHOUR+1,5,KR)+HCNDSS
               ADAILY(JHOUR+1,62,KR)=ADAILY(JHOUR+1,62,KR)+PRCPSS
            END IF
         END DO

C**** TOTAL PRECIPITATION AND AGE OF SNOW
      PREC(I,J)=PREC(I,J)+PRCPSS*100.*BYGRAV
      PRCP=PREC(I,J)
      PRECSS(I,J)=PRCPSS*100.*BYGRAV
C**** CALCULATE PRECIPITATION HEAT FLUX (FALLS AT 0 DEGREES CENTIGRADE)
      TPRCP=TPREC(I,J)
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
      EPREC(1,I,J)=EPRCP  ! assuming liquid water
      EPREC(2,I,J)=ENRGP  ! including latent heat
C**** PRECIPITATION DIAGNOSTICS
        AREG(JR,J_EPRCP)=AREG(JR,J_EPRCP)+ENRGP*DXYP(J)
        AIJ(I,J,IJ_PREC)=AIJ(I,J,IJ_PREC)+PRCP
        AIJ(I,J,IJ_NETH)=AIJ(I,J,IJ_NETH)+ENRGP

      IF(TPREC(I,J).GE.0.) PRCP=0.
C**** MODIFY SNOW AGES AFTER SNOW FALL
      GDATA(I,J,9)=GDATA(I,J,9)*EXP(-PRCP)
      GDATA(I,J,10)=GDATA(I,J,10)*EXP(-PRCP)
      GDATA(I,J,11)=GDATA(I,J,11)*EXP(-PRCP)

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
         AJL(J,L,11)=AJL(J,L,11)+AJ11(L)
         AJL(J,L,55)=AJL(J,L,55)+AJ55(L)

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

      END DO
C**** END OF MAIN LOOP FOR INDEX I
      END DO
C**** END OF MAIN LOOP FOR INDEX J
C****

C**** ADD IN CHANGE OF MOMENTUM BY MOIST CONVECTION AND CTEI
      DO L=1,LM
         DO J=2,JM
            DO I=1,IM
               AJL(J,L,39)=AJL(J,L,39)+(U(I,J,L)-UC(I,J,L))*PDSIG(L,I,J)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE CONDSE

      SUBROUTINE init_CLD
!@sum  init_CLD initialises parameters for MSTCNV and LSCOND
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE CONSTANT, only : grav,by3
      USE E001M12_COM, only : dtsrc,ls1
      USE CLD01, only : lmcm,bydtsrc,xmass,brcld,bybr

      IMPLICIT NONE

      LMCM = LS1-1
      BYDTsrc=1./DTsrc
      XMASS=0.1d0*DTsrc*GRAV

      BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3

      END SUBROUTINE init_CLD
