      SUBROUTINE MC_COND
!@sum   MSTCNV driver for moist convction
!@sum   CONDSE driver for large-scale condensation
!@auth  M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver   1.0 (taken from CB265)
!@calls MSTCNV_loc

      USE CONSTANT, only : rgas,grav,lhe,lhs,lhm,kapa,sha,bysha
      USE E001M12_COM, only : im,jm,lm,p,u,v,t,q,wm,tofday,fearth
     *     ,ls1,psf,ptop,dsig,bydsig,jeq,fland,ijd6,gdata,sig
      USE SOMTQ_COM
      USE GEOM, only : bydxyp,dxyp,imaxj,kmaxj,raj,idij,idjj
      USE CLD01_COM_E001
      USE CLD01, only : kmax,ra,pl,ple,plk
     *     ,airm,byam,etal,sm,sxm,sym,szm,sxxm,sxym,syym,syzm,szzm,szxm
     *     ,qm,qxm,qym,qzm,qxxm,qxym,qyym,qyzm,qzzm,qzxm,tl,aj13
     *     ,aj50,aj51,aj52,aj57,aj8,aj11,wml,sdl,u_0,v_0,um,vm,tf
     *     ,prcpmc,pearth,ts,bygrav,taumcl,cldmcl,svwmxl,svlatl,svlhxl
     *     ,cldslwij,clddepij,csizel,precnvl,vsubl,lmcmax,lmcmin,wmsum
     *     ,mstcnv_loc
     *     ,aq,dpdt,th,ql,wmx,ttoldl,rh,lpbl,taussl,cldssl,cldsavl,
     *     prcpss,hcndss,aj55,bydtcn,condse_loc
      USE PBLCOM, only : tsavg
      USE DAGCOM  !, only : aj,bj,cj,areg,aij,ajl,ail,adaily,jreg
      USE DYNAMICS, only : pk,pmid,pedn,sd_clouds,gz,ptold,pdsig
      USE OCEAN, only : odata

      IMPLICIT NONE

!@var UC,VC velocity work arrays
      REAL*8, DIMENSION(IM,JM,LM) :: UC,VC

      REAL*8,  PARAMETER :: ENTCON = .2d0  !@param ENTCON  ???

      INTEGER I,J,K,L  !@var I,J,K,L loop variables
      INTEGER IHOUR,IMAX,JR,KR
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID 

      REAL*8 :: HCNDMC,PRCP

C**** SAVE UC AND VC, AND ZERO OUT CLDSS AND CLDMC
      UC=U
      VC=V
         IHOUR=1.5+TOFDAY
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
      SXM(L) =TX(I,J,L)*AIRM(L)
      SYM(L) =TY(I,J,L)*AIRM(L)
      SZM(L) =TZ(I,J,L)*AIRM(L)
      SXXM(L)=TXX(I,J,L)*AIRM(L)
      SYYM(L)=TYY(I,J,L)*AIRM(L)
      SZZM(L)=TZZ(I,J,L)*AIRM(L)
      SXYM(L)=TXY(I,J,L)*AIRM(L)
      SYZM(L)=TYZ(I,J,L)*AIRM(L)
      SZXM(L)=TZX(I,J,L)*AIRM(L)
      TL(L)=T(I,J,L)*PLK(L)
C**** MOISTURE (SPECIFIC HUMIDITY)
      QM(L)  =Q(I,J,L)*AIRM(L)
      QXM(L) =QX(I,J,L)*AIRM(L)
      QYM(L) =QY(I,J,L)*AIRM(L)
      QZM(L) =QZ(I,J,L)*AIRM(L)
      QXXM(L)=QXX(I,J,L)*AIRM(L)
      QYYM(L)=QYY(I,J,L)*AIRM(L)
      QZZM(L)=QZZ(I,J,L)*AIRM(L)
      QXYM(L)=QXY(I,J,L)*AIRM(L)
      QYZM(L)=QYZ(I,J,L)*AIRM(L)
      QZXM(L)=QZX(I,J,L)*AIRM(L)
c      QL(L)=Q(I,J,L)
      WML(L)=WM(I,J,L)
      SDL(L)=SD_CLOUDS(I,J,L)*BYDXYP(J)
      SVLHXL(L)=SVLHX(I,J,L)
         TTOLDL(L)=TTOLD(I,J,L)
         CLDSAVL(L)=CLDSAV(I,J,L)
         RH(L)=RHSAV(I,J,L)
         DPDT(L)=SIG(L)*(P(I,J)-PTOLD(I,J))*BYDTCN
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

      CALL MSTCNV_LOC

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
     *        (1.-ODATA(I,J,2))
         BJ(J,J_PRCPMC)=BJ(J,J_PRCPMC)+PRCPMC*FLAND(I,J)
         CJ(J,J_PRCPMC)=CJ(J,J_PRCPMC)+PRCPMC*ODATA(I,J,2)*
     *        (1.-FLAND(I,J))
         AREG(JR,J_PRCPMC)=AREG(JR,J_PRCPMC)+PRCPMC*DXYP(J)
         DO KR=1,4
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(IHOUR,49,KR)=ADAILY(IHOUR,49,KR)+PRCPMC
               ADAILY(IHOUR,5,KR)=ADAILY(IHOUR,5,KR)+HCNDMC
               ADAILY(IHOUR,63,KR)=ADAILY(IHOUR,63,KR)+PRCPMC
               ADAILY(IHOUR,51,KR)=ADAILY(IHOUR,51,KR)+CLDDEPIJ
               ADAILY(IHOUR,52,KR)=ADAILY(IHOUR,52,KR)+CLDSLWIJ
            END IF
         END DO

C**** UPDATE MODEL TEMPERATURE, SPECIFIC HUMIDITY AND PRECIPITATION
         PREC(I,J)=PRCPMC*100.*BYGRAV
         DO L=1,LMCMAX
            T(I,J,L)=  SM(L)*BYAM(L)
            Q(I,J,L)=  QM(L)*BYAM(L)
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
         AQ(L)=(QL(L)-QTOLD(I,J,L))*BYDTCN
      END DO

C**** LARGE-SCALE CLOUDS AND PRECIPITATION

      CALL CONDSE_LOC(I,J)

C**** Accumulate diagnostics of CONDSE
         AIJ(I,J,IJ_WMSUM)=AIJ(I,J,IJ_WMSUM)+WMSUM
         AJ(J,J_PRCPSS)=AJ(J,J_PRCPSS)+PRCPSS*(1.-FLAND(I,J))*
     *        (1.-ODATA(I,J,2))
         BJ(J,J_PRCPSS)=BJ(J,J_PRCPSS)+PRCPSS*FLAND(I,J)
         CJ(J,J_PRCPSS)=CJ(J,J_PRCPSS)+PRCPSS*ODATA(I,J,2)*
     *        (1.-FLAND(I,J))
         AREG(JR,J_PRCPSS)=AREG(JR,J_PRCPSS)+PRCPSS*DXYP(J)
         DO KR=1,4
            IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
               ADAILY(IHOUR,49,KR)=ADAILY(IHOUR,49,KR)+PRCPSS
               ADAILY(IHOUR,5,KR)=ADAILY(IHOUR,5,KR)+HCNDSS
               ADAILY(IHOUR,62,KR)=ADAILY(IHOUR,62,KR)+PRCPSS
            END IF
         END DO
C**** TOTAL PRECIPITATION AND AGE OF SNOW
      PREC(I,J)=PREC(I,J)+PRCPSS*100.*BYGRAV
      PRCP=PREC(I,J)
      PRECSS(I,J)=PRCPSS*100.*BYGRAV
      IF(TPREC(I,J).GE.0.) PRCP=0.
C**** MODIFY SNOW AGES AFTER SNOW FALL
      GDATA(I,J,9)=GDATA(I,J,9)*EXP(-PRCP)
      GDATA(I,J,10)=GDATA(I,J,10)*EXP(-PRCP)
      GDATA(I,J,11)=GDATA(I,J,11)*EXP(-PRCP)

C**** WRITE TO GLOBAL ARRAYS
      DO L=1,LM
         TAUMC(L,I,J)=TAUMCL(L)
         CLDMC(L,I,J)=CLDMCL(L)
         SVLAT(I,J,L)=SVLATL(L)
         CSIZE(1,L,I,J)=CSIZEL(L,1)

         TAUSS(L,I,J)=TAUSSL(L)
         CLDSS(L,I,J)=CLDSSL(L)
         CLDSAV(I,J,L)=CLDSAVL(L)
         SVLHX(I,J,L)=SVLHXL(L)
         CSIZE(2,L,I,J)=CSIZEL(L,2)
         AJL(J,L,11)=AJL(J,L,11)+AJ11(L)
         AJL(J,L,55)=AJL(J,L,55)+AJ55(L)

         T(I,J,L)=TH(L)
         Q(I,J,L)=QL(L)
C**** update moment changes
         TX(I,J,L)= SXM(L)*BYAM(L)
         TY(I,J,L)= SYM(L)*BYAM(L)
         TZ(I,J,L)= SZM(L)*BYAM(L)
         TXX(I,J,L)=SXXM(L)*BYAM(L)
         TXY(I,J,L)=SXYM(L)*BYAM(L)
         TYY(I,J,L)=SYYM(L)*BYAM(L)
         TYZ(I,J,L)=SYZM(L)*BYAM(L)
         TZZ(I,J,L)=SZZM(L)*BYAM(L)
         TZX(I,J,L)=SZXM(L)*BYAM(L)
         QX(I,J,L)= QXM(L)*BYAM(L)
         QY(I,J,L)= QYM(L)*BYAM(L)
         QZ(I,J,L)= QZM(L)*BYAM(L)
         QXX(I,J,L)=QXXM(L)*BYAM(L)
         QXY(I,J,L)=QXYM(L)*BYAM(L)
         QYY(I,J,L)=QYYM(L)*BYAM(L)
         QYZ(I,J,L)=QYZM(L)*BYAM(L)
         QZZ(I,J,L)=QZZM(L)*BYAM(L)
         QZX(I,J,L)=QZXM(L)*BYAM(L)
         RHSAV(I,J,L)=RH(L)
         TTOLD(I,J,L)=TH(L)
         QTOLD(I,J,L)=QL(L)
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
      END SUBROUTINE MC_COND

      SUBROUTINE init_CLD
!@sum  init_CLD initialises parameters for MSTCNV and CONDSE
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE CONSTANT, only : rgas,grav,lhe,lhs,kapa,bysha,sday  !,by3
      USE E001M12_COM, only : dt,ncnds,LS1
      USE CLD01

      IMPLICIT NONE

      LMCM = LS1-1
      DTCNDS=NCNDS*DT
      BYDTCN=1./DTCNDS

      BXCONS=.622d0/RGAS
      AXCONS=LOG(6.1071D0)
      XMASS=0.1d0*DTCNDS*GRAV
      PK1000=1000.**KAPA
      BYBR=((1.-BRCLD)*(1.-2.*BRCLD))**BY3
      SLHE=LHE*BYSHA
      SLHS=LHS*BYSHA
      DQDTX=.622d0*LHE/RGAS
      DTPERD=DTCNDS/SDAY
      AGESNX=1.-DTPERD/50.

      END SUBROUTINE init_CLD
