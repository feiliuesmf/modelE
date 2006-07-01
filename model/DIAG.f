#include "rundeck_opts.h"
#define JJ(J) (J)-J_0H+1

!@sum  DIAG ModelE diagnostic calculations
!@auth G. Schmidt/J. Lerner/R. Ruedy/M. Kelley
!@ver  1.0
C**** AJ(J,N)  (ZONAL SUM OVER LONGITUDE AND TIME)
C****   See j_defs for contents
C****                                                             IDACC
C**** CONTENTS OF APJ(J,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  P (100 PA)                                              4 DA
C****   2  4*P4I (100 PA)  (UV GRID)                               4 DA
C****
C**** CONTENTS OF AJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   See jl_defs for contents
C****
C**** CONTENTS OF ASJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   See jls_defs for contents
C****
C**** CONTENTS OF AIJ(I,J,N)  (SUM OVER TIME OF)
C****   See ij_defs for contents
C****
C**** CONTENTS OF AIL(I,L,N)  (SUM OVER TIME OF)
C****   See il_defs for contents
C****
C**** CONTENTS OF IDACC(N), NUMBER OF ACCUMULATION TIMES OF
C****   1  SOURCE TERMS  (dt: DTSRC)
C****   2  RADIATION SOURCE TERMS  (dt: NRAD*DTsrc)
C****   3  SURFACE INTERACTION SOURCE TERMS  (dt: NDASF*DTsrc+DTsurf)
C****   4  QUANTITIES IN DIAGA  (dt: NDAA*DTsrc+2*DTdyn)
C****   5  ENERGY NUMBERS IN DIAG4  (DEYERMINED BY NDA4)
C****   6  KINETIC ENERGY IN DIAG5 FROM DYN'CS (dt: NDA5K*DTsrc+2*DTdyn)
C****   7  ENERGY IN DIAG5 FROM DYNAMICS  (dt: NDA5D*DTsrc)
C****   8  ENERGY IN DIAG5 FROM SOURCES  (DETERMINED BY NDA5S)
C****   9  WAVE ENERGY IN DIAG7  (dt: 12 HOURS)
C****  10  ENERGY IN DIAG5 FROM FILTER  (DT: NFILTR*DTsrc)
C****  11  NOT USED
C****  12  ALWAYS =1 (UNLESS SEVERAL RESTART FILES WERE ACCUMULATED)
C****

      MODULE DIAG_LOC
!@sum DIAG_LOC is a local module for some saved diagnostic calculations
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,imh,jm,lm
      IMPLICIT NONE
      SAVE
C**** Variables passed from DIAGA to DIAGB
!@var W,TX vertical velocity and in-situ temperature calculations
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: W
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TX
!@var TJL0 zonal mean temperatures prior to advection
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TJL0

C**** Variables used in DIAG5 calculations
!@var FCUVA,FCUVB fourier coefficients for velocities
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: FCUVA,FCUVB

C**** Some local constants
!@var JET, LDEX model levels for various pressures
!@var LUPA,LDNA shorthand for above/below levels
!@var PMO,PLO,PM,PL some shorthand pressure level
      INTEGER :: JET
      INTEGER, DIMENSION(3) :: LDEX
      REAL*8, DIMENSION(LM) :: LUPA,LDNA
      REAL*8, DIMENSION(LM) :: PMO,PLO
      REAL*8, DIMENSION(LM+1) :: PM,PL

      END MODULE DIAG_LOC

      SUBROUTINE ALLOC_DIAG_LOC(grid)
      USE DOMAIN_DECOMP, only : GET
      USE DOMAIN_DECOMP, only : DIST_GRID
      USE MODEL_COM, only : im,imh,lm
      USE DIAG_LOC, only  : W,TX,TJL0,FCUVA,FCUVB
      IMPLICIT NONE
      LOGICAL, SAVE :: init=.false.
      INTEGER :: J_1H    , J_0H
      INTEGER :: IER
      TYPE(DIST_GRID) :: grid

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE( W(IM, J_0H:J_1H, LM), TX(IM, J_0H:J_1H, LM),
     &     STAT = IER)
      ALLOCATE( TJL0(J_0H:J_1H, LM),
     &     STAT = IER)
      ALLOCATE( FCUVA(0:IMH, J_0H:J_1H, LM, 2),
     &          FCUVB(0:IMH, J_0H:J_1H, LM, 2),
     &     STAT = IER)

      !hack hack hack!
      TX(:,:,:) = 0.d0

      RETURN
      END SUBROUTINE ALLOC_DIAG_LOC

      SUBROUTINE DIAGA
!@sum  DIAGA accumulate various diagnostics during dynamics
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,lhe,sha,bygrav,bbyg,gbyrb,tf
     *     ,rvap,gamd,teeny,undef
      USE MODEL_COM, only : im,imh,fim,byim,jm,jeq,lm,ls1,idacc,ptop
     *     ,pmtop,psfmpt,mdyn,mdiag,sig,sige,dsig,zatmo,WM,ntype,ftype
     *     ,u,v,t,p,q
      USE GEOM, only : areag,cosp,dlat,dxv,dxyn,dxyp,dxys,dxyv,dyp,fcor
     *     ,imaxj,ravpn,ravps,sinp,bydxyv
      USE RAD_COM, only : rqt,lm_req
      USE DIAG_COM, only : aj=>aj_loc,areg,jreg,apj=>apj_loc
     *     ,ajl=>ajl_loc,asjl=>asjl_loc,ail,j50n,j70n,j5nuv
     *     ,j5suv,j5s,j5n,aij=>aij_loc
     *     ,ij_dtdp,ij_dsev,ij_phi1k,ij_pres,ij_puq
     *     ,ij_pvq,ij_slp,ij_t850,ij_t500,ij_t300,ij_q850,ij_q500
     *     ,ij_RH1,ij_RH850,ij_RH500,ij_RH300,ij_qm
     *     ,ij_q300,ij_ujet,ij_vjet,j_tx1,j_tx,j_qp,j_dtdjt,j_dtdjs
     *     ,j_dtdgtr,j_dtsgst,j_rictr,j_rostr,j_ltro,j_ricst,j_rosst
     *     ,j_lstr,j_gamm,j_gam,j_gamc,lstr,il_ueq,il_veq,il_weq,il_teq
     *     ,il_qeq,il_w50n,il_t50n,il_u50n,il_w70n,il_t70n,il_u70n
     *     ,kgz_max,pmb,ght,jl_dtdyn,jl_zmfntmom,jl_totntmom,jl_ape
     *     ,jl_uepac,jl_vepac,jl_uwpac,jl_vwpac,jl_wepac,jl_wwpac
     *     ,jl_epflxn,jl_epflxv,ij_p850,z_inst,rh_inst,t_inst,plm
     *     ,ij_p1000,ij_p925,ij_p700,ij_p600,ij_p500
      USE DYNAMICS, only : pk,phi,pmid,plij, pit,sd,pedn,am
      USE PBLCOM, only : tsavg
      USE DIAG_LOC, only : w,tx,lupa,ldna,jet,tjl0
      USE DOMAIN_DECOMP, only : GET, CHECKSUM, HALO_UPDATE,
     &                          CHECKSUM_COLUMN, HALO_UPDATE_COLUMN,
     &                          GRID, SOUTH, NORTH, GLOBALSUM
      USE DOMAIN_DECOMP, only : AM_I_ROOT
      IMPLICIT NONE
      REAL*8, DIMENSION(LM) :: GMEAN
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     GMEAN_part
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &        TIL,UI,UMAX,PI,EL,RI,DUDVSQ
      REAL*8, DIMENSION(NTYPE,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &        SPTYPE
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        THJL,THSQJL,SPI,PHIPI,TPI
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM-1) ::
     &        SDMEAN
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &        PUV
      REAL*8, DIMENSION(LM_REQ) :: TRI
      REAL*8, DIMENSION(IM) :: THSEC,PSEC,SQRTP,PDA
      REAL*8,
     & DIMENSION(size(areg,1),GRID%J_STRT_HALO:GRID%J_STOP_HALO,2)
     & :: AREG_part
      REAL*8 :: AREGSUM(size(areg,1),2)
      CHARACTER*16 TITLE
      REAL*8, PARAMETER :: ONE=1.,P1000=1000.
      INTEGER :: I,IM1,J,K,L,JR,LDN,LUP,
     &     IP1,LM1,LP1,LR,MBEGIN,IT
      REAL*8 THBAR ! external
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        THGM_part
      REAL*8, DIMENSION(LM):: THGM
      REAL*8 ::
     &     BBYGV,BDTDL,BYSDSG,CDTDL,DLNP,DLNP12,DLNP23,DBYSD,
     &     DLNS,DP,DS,DT2,DTHDP,DU,DUDP,DUDX,DV,DXYPJ,ELX,
     *     ESEPS,FPHI,GAMC,GAMM,GAMX,GMEANL,P1,P4,P4I,
     &     PDN,PE,PEQ,PEQM1,PEQM2,PHIRI,PIBYIM,PIJ,PITIJ,PITMN,
     *     PKE,PL,PRT,PU4I,PUP,PUV4I,PV4I,PVTHP,
     *     QLH,ROSSX,SDMN,SDPU,SMALL,SP,SP2,SS,T4,THETA,THMN,TPIL,
     *     TZL,UAMAX,UMN,UPE,VPE,X,Z4,THI,TIJK,QIJK
      LOGICAL qpress,qabove
      INTEGER nT,nQ,nRH
      REAL*8, PARAMETER :: EPSLON=1.
      INTEGER, PARAMETER ::
     *     I150E = IM*(180+150)/360+1, ! WEST EDGE OF 150 EAST
     *     I110W = IM*(180-110)/360+1, ! WEST EDGE OF 110 WEST
     *     I135W = IM*(180-135)/360+1  ! WEST EDGE OF 135 WEST

      REAL*8 QSAT
      REAL*8 :: gsum(IM, LM)
      REAL*8, dimension(:,:,:), allocatable :: TMP
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(MBEGIN)

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      IDACC(4)=IDACC(4)+1

      BYSDSG=1./(1.-SIGE(LM+1))
      DLNP12=LOG(.75/.35)
      DLNP23=LOG(.35/.1)
C****
C**** FILL IN HUMIDITY AND SIGMA DOT ARRAYS AT THE POLES
C****
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          DO I=2,IM
            Q(I,1,L)=Q(1,1,L)
          END DO
        END DO
      ENDIF        ! HAVE_SOUTH_POLE
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          DO I=2,IM
            Q(I,JM,L)=Q(1,JM,L)
          END DO
        END DO
      ENDIF        ! HAVE_NORTH_POLE
C****
C**** CALCULATE PK AND TX, THE REAL TEMPERATURE
C****
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          TX(1,1,L)=T(1,1,L)*PK(L,1,1)
          DO I=2,IM
            T(I,1,L)=T(1,1,L)
            TX(I,1,L)=TX(1,1,L)
          END DO
        END DO
      ENDIF        ! HAVE_SOUTH_POLE
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
          DO I=2,IM
            T(I,JM,L)=T(1,JM,L)
            TX(I,JM,L)=TX(1,JM,L)
          END DO
        END DO
      ENDIF          ! HAVE_NORTH_POLE
      DO L=1,LM
        DO J=J_0S,J_1S
          DO I=1,IM
            TX(I,J,L)=T(I,J,L)*PK(L,I,J)
          END DO
        END DO
      END DO
C****
C**** CALCULATE PUV, THE MASS WEIGHTED PRESSURE
C****
      CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      DO J=J_0STG,J_1STG
        I=IM
        DO IP1=1,IM
          PUV(I,J)=RAVPN(J-1)*(P(I,J-1)+P(IP1,J-1))+
     *             RAVPS(  J)*(P(I,  J)+P(IP1,  J))
          I=IP1
        END DO
      END DO
C****
C**** J LOOPS FOR ALL PRIMARY GRID ROWS
C****
      AREG_part=0.
      DO J=J_0,J_1
        DXYPJ=DXYP(J)
C**** NUMBERS ACCUMULATED FOR A SINGLE LEVEL
        PI(J)=0.
        SPTYPE(:,J)=0
        DO I=1,IMAXJ(J)
          JR=JREG(I,J)
          DO IT=1,NTYPE
            SPTYPE(IT,J)=SPTYPE(IT,J)+FTYPE(IT,I,J)
            AJ(J,J_TX1,IT)=AJ(J,J_TX1,IT)+(TX(I,J,1)-TF)*FTYPE(IT,I,J)
          END DO
C****         AREG(JR,J_TX1)=AREG(JR,J_TX1)+(TX(I,J,1)-TF)*DXYPJ
          AREG_part(JR,J,1)=AREG_part(JR,J,1)+(TX(I,J,1)-TF)*DXYPJ
          PI(J)=PI(J)+P(I,J)
          AIJ(I,J,IJ_PRES)=AIJ(I,J,IJ_PRES)+ P(I,J)
          AIJ(I,J,IJ_SLP)=AIJ(I,J,IJ_SLP)+((P(I,J)+PTOP)*(1.+BBYG
     *         *ZATMO(I,J)/TSAVG(I,J))**GBYRB-P1000)
          AIJ(I,J,IJ_RH1)=AIJ(I,J,IJ_RH1)+Q(I,J,1)/QSAT(TX(I,J,1),LHE,
     *        PMID(1,I,J))
        END DO
        APJ(J,1)=APJ(J,1)+PI(J)
C**** CALCULATE GEOPOTENTIAL HEIGHTS AT SPECIFIC MILLIBAR LEVELS
        DO I=1,IMAXJ(J)
          K=1
          L=1
          rh_inst(:,i,j) = undef ; t_inst(:,i,j) = undef
          z_inst(:,i,j) = undef
 172      L=L+1
          PDN=PMID(L-1,I,J)
          PL=PMID(L,I,J)
          IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
C**** Select pressure levels on which to save temperature and humidity
C**** Use masking for 850 mb temp/humidity
 174      qpress = .false.
          qabove = pmb(k).le.pedn(l-1,i,j)
          SELECT CASE (NINT(PMB(K)))
          CASE (850)            ! 850 mb
            nT = IJ_T850 ; nQ = IJ_Q850 ; nRH = IJ_RH850 ; qpress=.true.
            if (.not. qabove) qpress = .false.
            if (qpress) aij(i,j,ij_p850) = aij(i,j,ij_p850) + 1.
          CASE (500)            ! 500 mb
            nT = IJ_T500 ; nQ = IJ_Q500 ; nRH = IJ_RH500 ; qpress=.true.
          CASE (300)            ! 300 mb
            nT = IJ_T300 ; nQ = IJ_Q300 ; nRH = IJ_RH300 ; qpress=.true.
          END SELECT
C**** calculate geopotential heights + temperatures
          IF (ABS(TX(I,J,L)-TX(I,J,L-1)).GE.EPSLON) THEN
            BBYGV=(TX(I,J,L-1)-TX(I,J,L))/(PHI(I,J,L)-PHI(I,J,L-1))
            AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *           -TX(I,J,L)*((PMB(K)/PL)**(RGAS*BBYGV)-1.)/BBYGV-GHT(K)
     *           *GRAV)
            IF (qabove) then
              TIJK=(TX(I,J,L)-TF
     *           +(TX(I,J,L-1)-TX(I,J,L))*LOG(PMB(K)/PL)/LOG(PDN/PL))
              Z_inst(K,I,J)=(PHI(I,J,L)
     *           -TX(I,J,L)*((PMB(K)/PL)**(RGAS*BBYGV)-1.)/BBYGV-GHT(K)
     *             *GRAV)
            END IF
          ELSE
            AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *           -RGAS*TX(I,J,L)*LOG(PMB(K)/PL)-GHT(K)*GRAV)
            IF (qabove) then
              TIJK=TX(I,J,L)-TF
              Z_inst(K,I,J)=(PHI(I,J,L)
     *             -RGAS*TX(I,J,L)*LOG(PMB(K)/PL)-GHT(K)*GRAV)
            END IF
          END IF
          if (qabove) then
            QIJK=Q(I,J,L)+(Q(I,J,L-1)-Q(I,J,L))*(PMB(K)-PL)/(PDN-PL)
            RH_inst(K,I,J)=QIJK/qsat(TIJK+TF,LHE,PMB(K))
            T_inst(K,I,J) =TIJK
            if (qpress) then
              AIJ(I,J,nT)=AIJ(I,J,nT)+TIJK
              AIJ(I,J,nQ)=AIJ(I,J,nQ)+QIJK
              AIJ(I,J,nRH)=AIJ(I,J,nRH)+QIJK/qsat(TIJK+TF,LHE,PMB(K))
            end if
          end if
C****
          IF (K.LT.KGZ_max) THEN
            K=K+1
            IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
            GO TO 174
          END IF
C**** BEGIN AMIP
          IF((P(I,J)+PTOP).LT.1000.)AIJ(I,J,IJ_P1000)=
     *      AIJ(I,J,IJ_P1000)+1.
          IF((P(I,J)+PTOP).LT.925.)AIJ(I,J,IJ_P925)=AIJ(I,J,IJ_P925)+1.
          IF((P(I,J)+PTOP).LT.700.)AIJ(I,J,IJ_P700)=AIJ(I,J,IJ_P700)+1.
          IF((P(I,J)+PTOP).LT.600.)AIJ(I,J,IJ_P600)=AIJ(I,J,IJ_P600)+1.
          IF((P(I,J)+PTOP).LT.500.)AIJ(I,J,IJ_P500)=AIJ(I,J,IJ_P500)+1.
C**** END AMIP
        END DO
      END DO
C****          AREG(JR,J_TX1)=AREG(JR,J_TX1)+(TX(I,J,1)-TF)*DXYPJ

      CALL GLOBALSUM(grid,AREG_part(1:size(AREG,1),:,1:1),
     &    AREGSUM(1:size(AREG,1),1:1),ALL=.TRUE.)
      AREG(1:size(AREG,1),J_TX1)=AREG(1:size(AREG,1),J_TX1)
     &  +AREGSUM(1:size(AREG,1),1)

C**** ACCUMULATION OF TEMP., POTENTIAL TEMP., Q, AND RH
      AREG_part = 0.
      DO J=J_0,J_1
        DXYPJ=DXYP(J)
        DO L=1,LM
          TPI(J,L)=0.
          PHIPI(J,L)=0.
          SPI(J,L)=0.
          THI=0.
          DBYSD=DSIG(L)*BYSDSG
          DO I=1,IMAXJ(J)
            JR=JREG(I,J)
            PIJ=PLIJ(L,I,J)
            AIJ(I,J,IJ_QM)=AIJ(I,J,IJ_QM)+Q(I,J,L)*AM(L,I,J)
            DO IT=1,NTYPE
              AJ(J,J_TX,IT)=AJ(J,J_TX,IT)+(TX(I,J,L)-TF)*FTYPE(IT,I,J)
     *             *DBYSD
              AJ(J,J_QP,IT)=AJ(J,J_QP,IT)+(Q(I,J,L)+WM(I,J,L))*PIJ
     *             *DSIG(L)*FTYPE(IT,I,J)
            END DO
C****            AREG(JR,J_QP)=AREG(JR,J_QP)+(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)
C****     *           *DXYPJ
C****            AREG(JR,J_TX)=AREG(JR,J_TX)+(TX(I,J,L)-TF)*DBYSD*DXYPJ
            AREG_part(JR,J,1)=AREG_part(JR,J,1)+(Q(I,J,L)+WM(I,J,L))*
     *                        PIJ*DSIG(L)*DXYPJ
            AREG_part(JR,J,2)=AREG_part(JR,J,2)+
     *                        (TX(I,J,L)-TF)*DBYSD*DXYPJ
            TPI(J,L)=TPI(J,L)+(TX(I,J,L)-TF)*PIJ
            PHIPI(J,L)=PHIPI(J,L)+PHI(I,J,L)*PIJ
            SPI(J,L)=SPI(J,L)+T(I,J,L)*PIJ
            THI=THI+T(I,J,L)
          END DO
          AJL(J,L,JL_DTDYN)=AJL(J,L,JL_DTDYN)+THI-TJL0(J,L)
        END DO
      END DO

      CALL GLOBALSUM(grid,AREG_part(1:size(AREG,1),:,1:2),
     &   AREGSUM(1:size(AREG,1),1:2),ALL=.TRUE.)
      AREG(1:size(AREG,1),J_QP)=AREG(1:size(AREG,1),J_QP)
     &   +AREGSUM(1:size(AREG,1),1)
      AREG(1:size(AREG,1),J_TX)=AREG(1:size(AREG,1),J_TX)
     &   +AREGSUM(1:size(AREG,1),2)

C****
C**** NORTHWARD GRADIENT OF TEMPERATURE: TROPOSPHERIC AND STRATOSPHERIC
C****
      CALL HALO_UPDATE(grid, TX, FROM=NORTH+SOUTH)

      DO J=J_0S,J_1S
C**** MEAN TROPOSPHERIC NORTHWARD TEMPERATURE GRADIENT
        DO L=1,LS1-1
        DO I=1,IM
        DO IT=1,NTYPE
          AJ(J,J_DTDJT,IT)=AJ(J,J_DTDJT,IT)+(TX(I,J+1,L)-TX(I,J-1,L))
     *         *FTYPE(IT,I,J)*DSIG(L)
        END DO
        END DO
        END DO
C**** MEAN STRATOSPHERIC NORTHWARD TEMPERATURE GRADIENT
        DO L=LS1,LSTR
        DO I=1,IM
        DO IT=1,NTYPE
          AJ(J,J_DTDJS,IT)=AJ(J,J_DTDJS,IT)+(TX(I,J+1,L)-TX(I,J-1,L))
     *         *FTYPE(IT,I,J)*DSIG(L)
        END DO
        END DO
        END DO
      END DO
C****
C**** STATIC STABILITIES: TROPOSPHERIC AND STRATOSPHERIC
C****
      AREG_part=0.
      DO J=J_0,J_1
      DXYPJ=DXYP(J)
C**** OLD TROPOSPHERIC STATIC STABILITY
      DO I=1,IMAXJ(J)
        JR=JREG(I,J)
        SS=(T(I,J,LS1-1)-T(I,J,1))/(PHI(I,J,LS1-1)-PHI(I,J,1)+teeny)
        DO IT=1,NTYPE
          AJ(J,J_DTDGTR,IT)=AJ(J,J_DTDGTR,IT)+SS*FTYPE(IT,I,J)
        END DO
C****        AREG(JR,J_DTDGTR)=AREG(JR,J_DTDGTR)+SS*DXYPJ
        AREG_part(JR,J,1)=AREG_part(JR,J,1)+SS*DXYPJ
        AIJ(I,J,IJ_DTDP)=AIJ(I,J,IJ_DTDP)+SS
      END DO
C**** OLD STRATOSPHERIC STATIC STABILITY (USE LSTR as approx 10mb)
      DO I=1,IMAXJ(J)
        JR=JREG(I,J)
        SS=(T(I,J,LSTR)-T(I,J,LS1-1))/((PHI(I,J,LSTR)-PHI(I,J,LS1-1))
     *       +teeny)
        DO IT=1,NTYPE
          AJ(J,J_DTSGST,IT)=AJ(J,J_DTSGST,IT)+SS*FTYPE(IT,I,J)
        END DO
C****        AREG(JR,J_DTSGST)=AREG(JR,J_DTSGST)+SS*DXYPJ
        AREG_part(JR,J,2)=AREG_part(JR,J,2)+SS*DXYPJ
      END DO
C****
C**** NUMBERS ACCUMULATED FOR THE RADIATION EQUILIBRIUM LAYERS
C****
      DO LR=1,LM_REQ
        TRI(LR)=0.
        DO I=1,IMAXJ(J)
          TRI(LR)=TRI(LR)+RQT(LR,I,J)
        END DO
        ASJL(J,LR,1)=ASJL(J,LR,1)+(TRI(LR)-TF*IMAXJ(J))
      END DO
      PHIRI=0.
      DO I=1,IMAXJ(J)
        PHIRI=PHIRI+(PHI(I,J,LM)+RGAS*.5*(TX(I,J,LM)+RQT(1,I,J))
     *       *LOG((SIG(LM)*PSFMPT+PTOP)/PLM(LM+1)))
      END DO
      ASJL(J,1,2)=ASJL(J,1,2)+PHIRI
      PHIRI=PHIRI+RGAS*.5*(TRI(1)+TRI(2))*DLNP12
      ASJL(J,2,2)=ASJL(J,2,2)+PHIRI
      PHIRI=PHIRI+RGAS*.5*(TRI(2)+TRI(3))*DLNP23
      ASJL(J,3,2)=ASJL(J,3,2)+PHIRI
      END DO

C****        AREG(JR,J_DTDGTR)=AREG(JR,J_DTDGTR)+SS*DXYPJ
C****        AREG(JR,J_DTSGST)=AREG(JR,J_DTSGST)+SS*DXYPJ
      CALL GLOBALSUM(grid,AREG_part(1:size(AREG,1),:,1:2),
     &   AREGSUM(1:size(AREG,1),1:2),ALL=.TRUE.)
      AREG(1:size(AREG,1),J_DTDGTR)=AREG(1:size(AREG,1),J_DTDGTR)
     &    +AREGSUM(1:size(AREG,1),1)
      AREG(1:size(AREG,1),J_DTSGST)=AREG(1:size(AREG,1),J_DTSGST)
     &    +AREGSUM(1:size(AREG,1),2)

C****
C**** RICHARDSON NUMBER , ROSSBY NUMBER , RADIUS OF DEFORMATION
C****
C**** NUMBERS ACCUMULATED OVER THE TROPOSPHERE
      DO J=J_0STG,J_1STG
        DUDVSQ(J)=0.
        UMAX(J)=0.
        DO I=1,IM
          DU=U(I,J,LS1-1)-U(I,J,1)
          DV=V(I,J,LS1-1)-V(I,J,1)
          DUDVSQ(J)=DUDVSQ(J)+(DU*DU+DV*DV)*PUV(I,J)
        END DO
      END DO

      CALL HALO_UPDATE(grid, DUDVSQ, FROM=NORTH)

      DO J=J_0S,J_1S
        PIBYIM=PI(J)*BYIM
        DLNP=LOG((SIG(1)*PIBYIM+PTOP)/(SIG(LS1-1)*PIBYIM+PTOP))
        DLNS=LOG(SPI(J,LS1-1)/SPI(J,1))
        DS=SPI(J,LS1-1)-SPI(J,1)
        EL(J)=SQRT(DLNS/DLNP)
        RI(J)=DS*DLNP/(.5*(DUDVSQ(J)+DUDVSQ(J+1)))
      END DO
      DO L=1,LS1-1
        DO J=J_0STG,J_1STG
          UI(J)=0.
          DO I=1,IM
            UI(J)=UI(J)+U(I,J,L)
          END DO
        END DO

        CALL HALO_UPDATE(grid, UI, FROM=NORTH)

        DO J=J_0S,J_1S
          UAMAX=ABS(UI(J)+UI(J+1))
          IF (UAMAX.GT.UMAX(J)) UMAX(J)=UAMAX
        END DO
      END DO
      DO J=J_0S,J_1S
        ROSSX=DYP(J)/(DXYP(J)*SINP(J))
        ELX=1./SINP(J)
        DO IT=1,NTYPE
          AJ(J,J_RICTR,IT)=AJ(J,J_RICTR,IT)+RI(J)  *SPTYPE(IT,J)
          AJ(J,J_ROSTR,IT)=AJ(J,J_ROSTR,IT)+UMAX(J)*SPTYPE(IT,J)*ROSSX
          AJ(J,J_LTRO ,IT)=AJ(J,J_LTRO ,IT)+EL(J)  *SPTYPE(IT,J)*ELX
        END DO
      END DO
C**** NUMBERS ACCUMULATED OVER THE LOWER STRATOSPHERE
C**** LSTR is approx. 10mb level. This maintains consistency over
C**** the different model tops
      DO J=J_0STG,J_1STG
        DUDVSQ(J)=0.
        UMAX(J)=0.
        DO I=1,IM
          DU=U(I,J,LSTR)-U(I,J,LS1-1)
          DV=V(I,J,LSTR)-V(I,J,LS1-1)
          DUDVSQ(J)=DUDVSQ(J)+(DU*DU+DV*DV)*PUV(I,J)
        END DO
      END DO

      CALL HALO_UPDATE(grid, DUDVSQ, FROM=NORTH)

      DO J=J_0S,J_1S
        PIBYIM=PI(J)*BYIM
        DLNP=LOG((SIG(LS1-1)*PIBYIM+PTOP)/(SIG(LSTR)*PSFMPT+PTOP))
        DLNS=LOG(SPI(J,LSTR)/SPI(J,LS1-1))
        DS=SPI(J,LSTR)-SPI(J,LS1-1)
        EL(J)=SQRT(DLNS/DLNP)
        RI(J)=DS*DLNP/(.5*(DUDVSQ(J)+DUDVSQ(J+1)))
      END DO
      DO L=LS1,LSTR
        DO J=J_0STG,J_1STG
          UI(J)=0.
          DO I=1,IM
            UI(J)=UI(J)+U(I,J,L)
          END DO
        END DO

        CALL HALO_UPDATE(grid, UI, FROM=NORTH)

        DO J=J_0S,J_1S
          UAMAX=ABS(UI(J)+UI(J+1))
          IF (UAMAX.GT.UMAX(J)) UMAX(J)=UAMAX
        END DO
      END DO
      DO J=J_0S,J_1S
        ROSSX=DYP(J)/(DXYP(J)*SINP(J))
        ELX=1./SINP(J)
        DO IT=1,NTYPE
          AJ(J,J_RICST,IT)=AJ(J,J_RICST,IT)+RI(J)  *SPTYPE(IT,J)
          AJ(J,J_ROSST,IT)=AJ(J,J_ROSST,IT)+UMAX(J)*SPTYPE(IT,J)*ROSSX
          AJ(J,J_LSTR ,IT)=AJ(J,J_LSTR ,IT)+EL(J)  *SPTYPE(IT,J)*ELX
        END DO
      END DO
C****
C**** MEAN TROPOSPHERIC LAPSE RATES:  MOIST CONVECTIVE, ACTUAL,
C****    DRY ADIABATIC
C****
      X=RGAS*LHE*LHE/(SHA*RVAP)
      DO J=J_0,J_1
        GAMM=0.
        DO L=1,LS1-1
          TZL=TPI(J,L)/PI(J)+TF
          PRT=(SIG(L)*PI(J)*BYIM+PTOP)*RGAS*TZL
          ESEPS=QSAT(TZL,LHE,ONE)
          GAMM=GAMM+(PRT+LHE*ESEPS)/(PRT+X*ESEPS/TZL)*DSIG(L)
        END DO
        GAMX=(TPI(J,1)-TPI(J,LS1-1))/(PHIPI(J,LS1-1)-PHIPI(J,1))
        DO IT=1,NTYPE
          AJ(J,J_GAMM,IT)=AJ(J,J_GAMM,IT)+GAMM*SPTYPE(IT,J)
          AJ(J,J_GAM ,IT)=AJ(J,J_GAM ,IT)+GAMX*SPTYPE(IT,J)
        END DO
      END DO
C**** DRY ADIABATIC LAPSE RATE
      DO J=J_0,J_1
        TPIL=0.
        DO L=1,LS1-1
          TPIL=TPIL+TPI(J,L)*DSIG(L)
        END DO
        TIL(J)=TPIL/(PI(J)*(SIGE(1)-SIGE(LS1)))
      END DO

      CALL HALO_UPDATE(grid, TIL, FROM=NORTH+SOUTH)

      DO J=J_0S,J_1S
        X=SINP(J)*GRAV/(COSP(J)*RGAS*2.*DLAT)
        DT2=TIL(J+1)-TIL(J-1)
        GAMC=GAMD+X*DT2/(TIL(J)+TF)
        DO IT=1,NTYPE
          AJ(J,J_GAMC,IT)=AJ(J,J_GAMC,IT)+GAMC*SPTYPE(IT,J)
        END DO
      END DO
C****
C**** EASTWARD TRANSPORTS
C****

      CALL HALO_UPDATE(grid, U, FROM=NORTH)

      DO L=1,LM
      DO J=J_0S,J_1S
      I=IM
      DO IP1=1,IM
        AIJ(I,J,IJ_PUQ)=AIJ(I,J,IJ_PUQ)+(PLIJ(L,I,J)+PLIJ(L,IP1,J))*
     *       (U(I,J,L)+U(I,J+1,L))*(Q(I,J,L)+Q(IP1,J,L))*DSIG(L)
        I=IP1
      END DO
      END DO
      END DO
C****
C**** MOMENTUM, KINETIC ENERGY, NORTHWARD TRANSPORTS, ANGULAR MOMENTUM
C****

!Not necessary here, done above      CALL CHECKSUM(grid, P, __LINE__, __FILE__)
!Not necessary here, done above      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, PLIJ, FROM=SOUTH)
!Not necessary here, done above      CALL CHECKSUM(grid, TX, __LINE__, __FILE__)
!Not necessary here, done above      CALL HALO_UPDATE(grid, TX, FROM=SOUTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, Q, FROM=SOUTH)

      DO J=J_0STG,J_1STG
      P4I=0.
      I=IM
      DO IP1=1,IM
        P4=P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J)
        P4I=P4I+P4
        AIJ(I,J,IJ_UJET)=AIJ(I,J,IJ_UJET)+U(I,J,JET)
        AIJ(I,J,IJ_VJET)=AIJ(I,J,IJ_VJET)+V(I,J,JET)
        I=IP1
      END DO
      APJ(J,2)=APJ(J,2)+P4I
      DO L=1,LM
        PU4I=0.
        PV4I=0.
        PUV4I=0.
        I=IM
        DO IP1=1,IM
          P4=PLIJ(L,I,J-1)+PLIJ(L,IP1,J-1)+PLIJ(L,I,J)+PLIJ(L,IP1,J)
          IF(L.EQ.LS1) P4I=FIM*P4
          PU4I=PU4I+P4*U(I,J,L)
          PV4I=PV4I+P4*V(I,J,L)
          PUV4I=PUV4I+P4*U(I,J,L)*V(I,J,L)
          T4=TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L)
          Z4=PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L)
          AIJ(I,J,IJ_DSEV)=AIJ(I,J,IJ_DSEV)+P4*(SHA*T4+Z4)*V(I,J,L)
     *         *DSIG(L)*DXV(J)
          SP2=PLIJ(L,IP1,J-1)+PLIJ(L,IP1,J)
          AIJ(IP1,J,IJ_PVQ)=AIJ(IP1,J,IJ_PVQ)+SP2
     *         *(V(I,J,L)+V(IP1,J,L))*(Q(IP1,J-1,L)+Q(IP1,J,L))*DSIG(L)
          I=IP1
        END DO
        AJL(J,L,JL_ZMFNTMOM)=AJL(J,L,JL_ZMFNTMOM)+PU4I*PV4I/P4I
        AJL(J,L,JL_TOTNTMOM)=AJL(J,L,JL_TOTNTMOM)+PUV4I
      END DO
      END DO
C****
C**** EVEN LEVEL GEOPOTENTIALS, VERTICAL WINDS AND VERTICAL TRANSPORTS
C****
      DO L=1,LM-1
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        PIJ=PLIJ(L,I,J)
        PE=SIGE(L+1)*PIJ+PTOP
        PKE=PE**KAPA
        THETA=THBAR(T(I,J,L+1),T(I,J,L))
        W(I,J,L)=SD(I,JJ(J),L)*THETA*PKE/PE
      END DO
      END DO
      END DO
C****
C**** AVAILABLE POTENTIAL ENERGY
C****
C**** SET UP FOR CALCULATION
      DO 710 L=1,LM
  710 GMEAN_part(:,L)=0.
      DO 740 J=J_0,J_1
      DO 720 I=1,IMAXJ(J)
  720 SQRTP(I)=SQRT(P(I,J))
C**** GMEAN CALCULATED FOR EACH LAYER, THJL, THSQJL ARRAYS FILLED
      DO 730 L=1,LM
      LDN=LDNA(L)
      LUP=LUPA(L)
      THJL(J,L)=0.
      THSQJL(J,L)=0.
      DO 730 I=1,IMAXJ(J)
      THJL(J,L)=THJL(J,L)+T(I,J,L)*SQRTP(I)
      THSQJL(J,L)=THSQJL(J,L)+T(I,J,L)*T(I,J,L)*P(I,J)
  730 GMEAN_part(J,L)=GMEAN_part(J,L)+
     *     (SIG(L)*P(I,J)+PTOP)*(T(I,J,LUP)-T(I,J,LDN))*
     *     DXYP(J)/(P(I,J)*PK(L,I,J))
  740 CONTINUE
C**** CALCULATE APE
      DO 760 L=1,LM
      IF(HAVE_SOUTH_POLE) THEN
        THJL(1,L)=THJL(1,L)*FIM
        THSQJL(1,L)=THSQJL(1,L)*FIM
      ENDIF
      IF(HAVE_NORTH_POLE) THEN
        THJL(JM,L)=THJL(JM,L)*FIM
        THSQJL(JM,L)=THSQJL(JM,L)*FIM
      ENDIF
      THGM_part(:,L)=0.
      DO J=J_0,J_1
        THGM_part(J,L)=THJL(J,L)*DXYP(J)
      ENDDO
  760 continue

      CALL GLOBALSUM(grid,THGM_part(:,1:LM),THGM(1:LM),ALL=.TRUE.)
      THGM=THGM/AREAG
      CALL GLOBALSUM(grid,GMEAN_part(:,1:LM),GMEAN(1:LM),ALL=.TRUE.)

      DO 761 L=1,LM
      LP1=LUPA(L)
      LM1=LDNA(L)
      GMEANL=GMEAN(L)/((SIG(LM1)-SIG(LP1))*AREAG)
      DO 761 J=J_0,J_1
  761 AJL(J,L,JL_APE)=AJL(J,L,JL_APE)+
     &        (THSQJL(J,L)-2.*THJL(J,L)*THGM(L)+THGM(L)*THGM(L)*
     *  FIM)/GMEANL
C****
C**** CERTAIN HORIZONTAL WIND AVERAGES
C****
      DO L=1,LM
      DO J=J_0STG,J_1STG
      DO I=I135W,I110W      ! EAST PACIFIC
        AJL(J,L,JL_UEPAC)=AJL(J,L,JL_UEPAC)+U(I,J,L)
        AJL(J,L,JL_VEPAC)=AJL(J,L,JL_VEPAC)+V(I,J,L)
      END DO
      DO I=I150E,IM             ! WEST PACIFIC
        AJL(J,L,JL_UWPAC)=AJL(J,L,JL_UWPAC)+U(I,J,L)
        AJL(J,L,JL_VWPAC)=AJL(J,L,JL_VWPAC)+V(I,J,L)
      END DO
      END DO
      END DO


      CALL GLOBALSUM(grid, U(1:IM,:,1:LM), gsum,
     &      jband=(/J5SUV,J5NUV/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_UEQ) =
     &       AIL(1:IM,1:LM,IL_UEQ)+gsum
      CALL GLOBALSUM(grid, V(1:IM,:,1:LM), gsum,
     &       jband=(/J5SUV,J5NUV/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_VEQ) =
     &       AIL(1:IM,1:LM,IL_VEQ)+gsum

      CALL GLOBALSUM(grid, TX(1:IM,:,1:LM)-TF, gsum,
     &       jband=(/J5S,J5N/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_TEQ) =
     &       AIL(1:IM,1:LM,IL_TEQ)+gsum
      allocate(TMP(1:IM,
     &  grid%J_STRT_HALO:grid%J_STOP_HALO,1:LM))
      DO L=1,LM
      DO J = MAX(J_0,J5S),MIN(J_1,J5N)
      DO I=1,IM
         TMP(I,J,L) = Q(I,J,L)/QSAT(TX(I,J,L),LHE,PMID(L,I,J))
      END DO
      END DO
      END DO
      CALL GLOBALSUM(grid, TMP(1:IM,:,1:LM), gsum,
     &       jband=(/J5S,J5N/))
      deallocate(TMP)
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_QEQ) =
     &       AIL(1:IM,1:LM,IL_QEQ)+gsum

      CALL GLOBALSUM(grid, TX(1:IM,:,1:LM)-TF, gsum,
     &       jband=(/J50N,J50N/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_T50N)=
     &       AIL(1:IM,1:LM,IL_T50N)+gsum
      CALL GLOBALSUM(grid, U(1:IM,:,1:LM), gsum,
     &       jband =(/J50N,J50N+1/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_U50N)=
     &       AIL(1:IM,1:LM,IL_U50N)+gsum

      CALL GLOBALSUM(grid, TX(1:IM,:,1:LM)-TF, gsum,
     &       jband=(/J70N,J70N/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_T70N)=
     &       AIL(1:IM,1:LM,IL_T70N)+gsum
      CALL GLOBALSUM(grid, U(1:IM,:,1:LM), gsum,
     &       jband =(/J70N,J70N+1/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM,IL_U70N)=
     &       AIL(1:IM,1:LM,IL_U70N)+gsum

C****
C**** CERTAIN VERTICAL WIND AVERAGES
C****
      DO L=1,LM-1
        DO J=J_0S,J_1S
          DO I=I135W,I110W      ! EAST PACIFIC
            AJL(J,L,JL_WEPAC)=AJL(J,L,JL_WEPAC)+W(I,J,L)
          END DO
          DO I=I150E,IM         ! WEST PACIFIC
            AJL(J,L,JL_WWPAC)=AJL(J,L,JL_WWPAC)+W(I,J,L)
          END DO
        END DO
      END DO
      CALL GLOBALSUM(grid, W(1:IM,:,1:LM-1), gsum(1:IM,1:LM-1),
     &     jband=(/J5S,J5N/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM-1,IL_WEQ) =
     &    AIL(1:IM,1:LM-1,IL_WEQ)+gsum(1:IM,1:LM-1)
      CALL GLOBALSUM(grid, W(1:IM,:,1:LM-1), gsum(1:IM,1:LM-1),
     &     jband=(/J50N,J50N/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM-1,IL_W50N) =
     &     AIL(1:IM,1:LM-1,IL_W50N)+gsum(1:IM,1:LM-1)
      CALL GLOBALSUM(grid, W(1:IM,:,1:LM-1), gsum(1:IM,1:LM-1),
     &     jband=(/J70N,J70N/))
      IF (AM_I_ROOT()) AIL(1:IM,1:LM-1,IL_W70N) =
     &     AIL(1:IM,1:LM-1,IL_W70N)+gsum(1:IM,1:LM-1)
C****
C**** ELIASSEN PALM FLUX
C****
C**** NORTHWARD TRANSPORT
!Not necessary here, done above      CALL CHECKSUM(grid, P, __LINE__, __FILE__)
!Not necessary here, done above      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      CALL HALO_UPDATE(grid, T, FROM=SOUTH)

      DO 868 J=J_0STG,J_1STG
      I=IM
      DO 862 IP1=1,IM
      PDA(I)=.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+P(IP1,J-1))*
     *  DXYN(J-1))
      PSEC(I)=PDA(I)*BYDXYV(J)
  862 I=IP1
      DO 868 L=1,LM
      DUDP=0.
      DTHDP=0.
      UMN=0.
      THMN=0.
      LDN=LDNA(L)
      LUP=LUPA(L)
      I=IM
      DO 864 IP1=1,IM
      DUDP=DUDP+U(I,J,LUP)-U(I,J,LDN)
      DTHDP=DTHDP+T(I,J,LUP)+T(I,J-1,LUP)-T(I,J,LDN)-T(I,J-1,LDN)
      UMN=UMN+U(I,J,L)
      THMN=THMN+T(I,J,L)+T(I,J-1,L)
      THSEC(I)=T(I,J,L)+T(IP1,J,L)+T(I,J-1,L)+T(IP1,J-1,L)
  864 I=IP1
      UMN=UMN*BYIM
      THMN=2.*THMN/FIM
      FPHI=0.
      SMALL=.0002*FIM*T(1,J,L)
c      IF (DTHDP.LT.SMALL) WRITE (6,999) J,L,DTHDP,SMALL
      IF (DTHDP.LT.SMALL) DTHDP=SMALL
      DO 866 I=1,IM
      SP=PSEC(I)
      IF(L.GE.LS1) SP=PSFMPT
  866 FPHI=FPHI+SP*V(I,J,L)*(.5*(THSEC(I)-THMN)*DUDP/DTHDP
     *   -U(I,J,L)+UMN)
  868 AJL(J,L,JL_EPFLXN)=AJL(J,L,JL_EPFLXN)+FPHI
C**** VERTICAL TRANSPORT
!Not necessary here, done above      CALL CHECKSUM(grid, U, __LINE__, __FILE__)
!Not necessary here, done above      CALL HALO_UPDATE(grid, U, FROM=NORTH)
      CALL HALO_UPDATE(grid, V, FROM=NORTH)

      DO 878 J=J_0S,J_1S
      PITMN=0.
      DO 870 I=1,IM
  870 PITMN=PITMN+PIT(I,JJ(J))
      PITMN=PITMN/FIM
      DO 878 L=1,LM-1
      IF(L.GE.LS1-1) PITMN=0.
      THMN=0.
      SDMN=0.
      DTHDP=0.
      DO 872 I=1,IM
      DTHDP=DTHDP+T(I,J,L+1)-T(I,J,L)
      THMN=THMN+T(I,J,L+1)+T(I,J,L)
  872 SDMN=SDMN+SD(I,JJ(J),L)
      SMALL=.0001*FIM*T(1,J,L+1)
c      IF (DTHDP.LT.SMALL) WRITE (6,999) J,L,DTHDP,SMALL
      IF (DTHDP.LT.SMALL) DTHDP=SMALL
      THMN=THMN/FIM
      SDMN=SDMN/FIM
      DUDX=0.
      PVTHP=0.
      SDPU=0.
      IM1=IM
      DO 874 I=1,IM
      DUDX=DUDX+DXV(J+1)*(U(I,J+1,L)+U(I,J+1,L+1))-DXV(J)*
     *   (U(I,J,L)+U(I,J,L+1))
      UPE=U(IM1,J,L)+U(IM1,J+1,L)+U(I,J,L)+U(I,J+1,L)+
     *    U(IM1,J,L+1)+U(IM1,J+1,L+1)+U(I,J,L+1)+U(I,J+1,L+1)
      VPE=V(IM1,J,L)+V(IM1,J+1,L)+V(I,J,L)+V(I,J+1,L)+
     *    V(IM1,J,L+1)+V(IM1,J+1,L+1)+V(I,J,L+1)+V(I,J+1,L+1)
      DP=(SIG(L)-SIG(L+1))*P(I,J)
      IF(L.GE.LS1) DP=(SIG(L)-SIG(L+1))*PSFMPT
      IF(L.EQ.LS1-1) DP=P(I,J)*SIG(L)-PSFMPT*SIG(LS1)
      PVTHP=PVTHP+DP*VPE*(T(I,J,L)+T(I,J,L+1)-THMN)
      PITIJ=PIT(I,JJ(J))
      IF(L.GE.LS1-1) PITIJ=0.
      SDPU=SDPU+(SD(I,JJ(J),L)-SDMN+(PITIJ-PITMN)*SIGE(L+1))*UPE
  874 IM1=I
      AJL(J,L,JL_EPFLXV)=AJL(J,L,JL_EPFLXV)+
     &     (.5*FIM*FCOR(J)-.25*DUDX)*PVTHP/DTHDP + SDPU
  878 CONTINUE

C**** ACCUMULATE TIME USED IN DIAGA
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
  999 FORMAT (' DTHETA/DP IS TOO SMALL AT J=',I4,' L=',I4,2F15.6)

      ENTRY DIAGA0
C****
C**** INITIALIZE TJL0 ARRAY (FROM PRIOR TO ADVECTION)
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      DO L=1,LM
      DO J=J_0,J_1
        THI=0.
        DO I=1,IMAXJ(J)
          THI=THI+T(I,J,L)
        END DO
        TJL0(J,L)=THI
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE DIAGA


      SUBROUTINE DIAGB
!@sum DIAGB calculate constant pressure diagnostics from within DYNAM
C****
C**** CONTENTS OF AJK(J,K,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   See jks_defs for contents
C****
C**** CONTENTS OF AIJK(I,J,K,N)   (SUM OVER TIME OF)
C****   See ijks_defs for contents
C****
      USE CONSTANT, only : lhe,omega,sha,tf,teeny
      USE MODEL_COM, only :
     &     im,imh,fim,byim,jm,jeq,lm,ls1,idacc,ptop,psfmpt,jdate,
     &     mdyn,mdiag, ndaa,sig,sige,dsig,Jhour,u,v,t,p,q,wm
      USE GEOM, only : bydxyp,
     &     COSV,DXV,DXYN,DXYP,DXYS,DXYV,DYP,DYV,FCOR,IMAXJ,RADIUS
      USE DIAG_COM, only : ajk=>ajk_loc,aijk=>aijk_loc,speca,nspher,  ! adiurn,hdiurn
     &     nwav_dag,ndiupt,hr_in_day,ijk_u,ijk_v,ijk_t,ijk_q,ijk_dp
     *     ,ijk_dse,klayer,idd_w,ijdd,ijk_r,ijk_w,ijk_pf,
     *      ijk_uv,ijk_vt,ijk_vq,ijk_vv,ijk_uu,ijk_tt,
     &      JK_DPA,JK_DPB,JK_TEMP,JK_HGHT,JK_Q,JK_THETA,
     &      JK_RH,JK_U,JK_V,JK_ZMFKE,JK_TOTKE,JK_ZMFNTSH,
     &      JK_TOTNTSH,JK_ZMFNTGEO,JK_TOTNTGEO,JK_ZMFNTLH,
     &      JK_TOTNTLH,JK_ZMFNTKE,JK_TOTNTKE,JK_ZMFNTMOM,JK_TOTNTMOM,
     &      JK_P2KEDPGF,JK_DPSQR,JK_NPTSAVG,
     &      JK_VVEL,JK_ZMFVTDSE,JK_TOTVTDSE,JK_ZMFVTLH,JK_TOTVTLH,
     &      JK_VTGEOEDDY,JK_BAREKEGEN,JK_POTVORT,JK_VTPV,
     &      JK_VTPVEDDY,JK_NPTSAVG1,JK_TOTVTKE,
     &      JK_VTAMEDDY,JK_TOTVTAM,JK_SHETH,JK_DUDTMADV,JK_DTDTMADV,
     &      JK_DUDTTEM,JK_DTDTTEM,JK_EPFLXNCP,JK_EPFLXVCP,
     &      JK_UINST,JK_TOTDUDT,JK_TINST,
     &      JK_TOTDTDT,JK_EDDVTPT,JK_CLDH2O
      USE DYNAMICS, only : phi,dut,dvt,plij,sd
      USE DIAG_LOC, only : w,tx,pm,pl,pmo,plo
      USE DOMAIN_DECOMP, only : GET, CHECKSUM, HALO_UPDATE, GRID
      USE DOMAIN_DECOMP, only : HALO_UPDATEj
      USE DOMAIN_DECOMP, only : SOUTH, NORTH, GLOBALSUM
      IMPLICIT NONE
      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE
      REAL*8, DIMENSION
     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &     ZX,STB,UDX
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &     STJK,DPJK,UJK,VJK,WJK,TJK,
     &     PSIJK,UP,TY,PSIP,WTJK,UVJK,WUJK
      REAL*8, DIMENSION(IM) :: PSEC,X1
      REAL*8, DIMENSION(LM) :: SHETH,DPM,DTH
c      REAL*8, DIMENSION(size(adiurn,1),size(adiurn,3),
c     &        GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: ADIURN_part
c      REAL*8, DIMENSION(size(hdiurn,1),size(hdiurn,3),
c     &        GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: HDIURN_part
c      REAL*8 :: ADIURNSUM,HDIURNSUM

      INTEGER ::
     &     I,IH,IHM,IM1,INCH,INCHM,IP1,IZERO,J,J45N,
     &     JHEMI,K,KDN,KR,KS,KS1,KSPHER,KUP,KX,L,
     &     LUP,MBEGIN,N,NM,KM

      REAL*8 ::
     &     BYDP,BYFIM,DP,DPDN,DPDX,DPDY,
     &     DPE,DPI,DPK,DPSQI,DPUP,DPUV,DUTI,DUTK,DVTI,DVTK,FIMI,
     &     PAI,PAK,PDN,PHIPI,PMK,PQ4I,PQ4K,PQV4I,PS,PS4I,
     &     PS4K,PSIY,PSV4I,PT4I,PT4K,PTK,PTV4I,PUI,PUK,PUP,
     &     PUVI,PV2,PV2I,PVI,PVK,PWWI,PWWVI,PY,PZ4I,PZ4K,
     &     PZV4I,QK,QKI,QLH,QPI,QSATL,RHPI,SDK,
     &     SMALL,SP,SP2,SQRTDP,THK,THKI,THPI,TK,TKI,TPI,
     &     UDUTI,    UEARTH,UK,UKI,UY,VDVTI,VK,VSTAR,W2,W2I,W4,
     &     W4I,WI,WKE4I,WMPI,WNP,WPA2I,WPV4I,WQI,WSP,WSTAR,WTHI,
     &     WTI,WU4I,WUP,WZI,ZK,ZKI
     &     ,AMRHT,AMRHQ,AMUV,AMVQ,AMVT,AMUU,AMVV,AMTT

      REAL*8, PARAMETER :: BIG=1.E20
      REAL*8 :: QSAT
      REAL*8 :: pm_ge_ps(im,grid%j_strt_halo:grid%j_stop_halo,lm)
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(MBEGIN)

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      pm_ge_ps(:,:,:)=-1.
C****
C**** INTERNAL QUANTITIES T,TH,Q,RH
C****
      KM=LM
      QLH=LHE
      DO 170 J=J_0,J_1
      DO 170 K=1,KM
      DPI=0.
      TPI=0.
      PHIPI=0.
      QPI=0.
      WMPI=0.
      THPI=0.
      RHPI=0.
      FIMI=0.
      DO 160 I=1,IMAXJ(J)
C**** FIND L=L(K) AND LUP=L(K+1) S.T. P(LUP).GT.P(K+1)
      SP=PLIJ(K,I,J)
      PS=SP+PTOP
      IF (PM(K+1).GE.PS) GO TO 160
      L=1
      PDN=PS
      IF (PM(K).GE.PS) GO TO 120
      PDN=PM(K)
  110 IF (PM(K).GT.SP*SIGE(L+1)+PTOP) GO TO 120
      L=L+1
      GO TO 110
  120 LUP=L
  130 IF (PM(K+1).GE.SP*SIGE(LUP+1)+PTOP) GO TO 140
      LUP=LUP+1
      GO TO 130
  140 CONTINUE
C**** ACCUMULATE HERE
      DPI=DPI+PDN-PM(K+1)
      FIMI=FIMI+1.
  150 PUP=SP*SIGE(L+1)+PTOP
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      TPI=TPI+(TX(I,J,L)-TF)*DP
      PHIPI=PHIPI+PHI(I,J,L)*DP
      QPI=QPI+Q(I,J,L)*DP
      WMPI=WMPI+WM(I,J,L)*DP
CW       IF(WMPI.GT.1.E-3) WRITE(6,169) I,J,L,DP,WM(I,J,L),WMPI
CW169 FORMAT(1X,'1616--',3I5,3E15.2)
      THPI=THPI+T(I,J,L)*DP
      QSATL=QSAT(TX(I,J,L),QLH,SIG(L)*SP+PTOP)
      IF (QSATL.GT.1.) QSATL=1.
      RHPI=RHPI+Q(I,J,L)*DP/QSATL
      IF (L.EQ.LUP) GO TO 160
      L=L+1
      PDN=SP*SIGE(L)+PTOP
      GO TO 150
  160 CONTINUE
      AJK(J,K,JK_NPTSAVG1)=AJK(J,K,JK_NPTSAVG1)+FIMI
      AJK(J,K,JK_DPA)=AJK(J,K,JK_DPA)+DPI
      AJK(J,K,JK_TEMP)=AJK(J,K,JK_TEMP)+TPI
      AJK(J,K,JK_HGHT)=AJK(J,K,JK_HGHT)+PHIPI
      AJK(J,K,JK_Q)=AJK(J,K,JK_Q)+QPI
      AJK(J,K,JK_THETA)=AJK(J,K,JK_THETA)+THPI
      AJK(J,K,JK_RH)=AJK(J,K,JK_RH)+RHPI
      AJK(J,K,JK_CLDH2O)=AJK(J,K,JK_CLDH2O)+WMPI
         TJK(J,K)=THPI/(DPI+teeny)
         IF (IDACC(4).EQ.1) AJK(J,K,JK_TINST)=TJK(J,K)
         AJK(J,K,JK_TOTDTDT)=TJK(J,K)-AJK(J,K,JK_TINST)
  170 CONTINUE
C****
C**** CALCULATE STABILITY AT ODD LEVELS ON PU GRID
C****
      DO 230 J=J_0,J_1
      I=IMAXJ(J)
      DO 230 IP1=1,IMAXJ(J)
      SP2=P(I,J)+P(IP1,J)
      SP=.5*SP2
      DO 175 L=1,LS1-1
      PLO(L)=SP*SIG(L)+PTOP
  175 PL(L)=SP*SIGE(L)+PTOP
      DO 180 L=1,LM-1
      DTH(L)=(T(I,J,L)+T(IP1,J,L)-T(I,J,L+1)-T(IP1,J,L+1))/
     *  (2.*(PLO(L)-PLO(L+1)))
  180 CONTINUE
      DO 220 K=1,KM
      STB(I,J,K)=0.
      IF (PM(K+1).GE.PL(1)) GO TO 220
      PMK=PMO(K)
      IF (PM(K).GT.PL(1)) PMK=.5*(SP+PTOP+PM(K+1))
      L=2
      IF (PMK.GE.PL(2)) GO TO 210
  190 LUP=L+1
      IF (L.EQ.LM) GO TO 210
      IF (PMK.GE.PL(LUP)) GO TO 200
      L=LUP
      GO TO 190
  200 DPUP=PMK-PL(LUP)
      DPDN=PL(L)-PMK
      STB(I,J,K)=(DTH(L-1)*DPUP+DTH(L)*DPDN)/(DPUP+DPDN+teeny)
      GO TO 220
C**** SPECIAL CASES,  L=2, L=LM
  210 STB(I,J,K)=DTH(L-1)
  220 CONTINUE
  230 I=IP1
C**** CALCULATE STJK; THE MEAN STATIC STABILITY
      DO 260 J=J_0,J_1
      DO 260 K=1,KM
      STJK(J,K)=0.
      DPJK(J,K)=0.
      I=IMAXJ(J)
      DO 250 IP1=1,IMAXJ(J)
      PS=.5*(P(I,J)+P(IP1,J))+PTOP
      IF (PM(K+1).GT.PS) GO TO 250
      STJK(J,K)=STJK(J,K)+STB(I,J,K)
      DPJK(J,K)=DPJK(J,K)+1.
  250 I=IP1
      STJK(J,K)=STJK(J,K)/(DPJK(J,K)+teeny)
      SMALL=.0001
      IF (ABS(STJK(J,K)).LT.SMALL) STJK(J,K)=-SMALL
  260 CONTINUE
C****
C**** CONSTANT PRESSURE DIAGNOSTICS:  FLUX, ENERGY, ANGULAR MOMENTUM
C****
      ZX(:,:,:)=0.
      IF(HAVE_SOUTH_POLE) THEN
        UDX(:,1,:)=0.
      ENDIF

C Needs to check later if all these halo calls are necessary.
C DIAGA may have contained relevant halo calls
C and since DIAGB is called immediately after DIAGA
C there may not be a need for these calls if
C the concerned arrays have not been updated
C from the previous halo call.
      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      CALL HALO_UPDATE(grid, TX, FROM=SOUTH)
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, Q, FROM=SOUTH)
      CALL HALO_UPDATE(grid, T, FROM=SOUTH)
      CALL HALO_UPDATEj(grid, STJK, FROM=SOUTH)
c***      DO L=1,LM
c***         CALL HALO_UPDATE(grid, STJK(:,L), FROM=SOUTH)
c***      END DO

      DO 390 J=J_0STG,J_1STG
      I=IM
      DO 280 IP1=1,IM
      PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
      DO  K=1,KM
        UDX(I,J,K)=0.
      END DO
      DO 275 L=1,LS1-1
      DUT(I,J,L)=DUT(I,J,L)/(PSEC(I)*DXYV(J)*DSIG(L))
  275 DVT(I,J,L)=DVT(I,J,L)/(PSEC(I)*DXYV(J)*DSIG(L))
      DO 276 L=LS1,LM
      DUT(I,J,L)=DUT(I,J,L)/(PSFMPT*DXYV(J)*DSIG(L))
  276 DVT(I,J,L)=DVT(I,J,L)/(PSFMPT*DXYV(J)*DSIG(L))
  280 I=IP1
      DO 350 K=1,KM
      DPI=0.
      DPSQI=0.
      FIMI=0.
      PUI=0.
      PVI=0.
      PWWI=0.
      PT4I=0.
      PTV4I=0.
      PZ4I=0.
      PZV4I=0.
      PQ4I=0.
      PQV4I=0.
      PWWVI=0.
      PUVI=0.
      DVTI=0.
      VDVTI=0.
      DUTI=0.
      UDUTI=0.
      PS4I=0.
      PSV4I=0.
      I=IM
      DO 340 IP1=1,IM
      SP=PSEC(I)
      PS=SP+PTOP
      DO 286 L=1,LS1-1
  286 PL(L)=SP*SIGE(L)+PTOP
      IF (PM(K+1).GE.PS) THEN
        pm_ge_ps(i,j,k) = 1.
        UDX(I,J,K)=BIG
      ELSE
        L=1
        PDN=PS
        IF (PM(K).GE.PS) GO TO 300
        PDN=PM(K)
  290   IF (PM(K).GT.PL(L+1)) GO TO 300
        L=L+1
        GO TO 290
  300   LUP=L
  310   IF (PM(K+1).GE.PL(LUP+1)) GO TO 320
        LUP=LUP+1
        GO TO 310
  320   CONTINUE
        DPK=PDN-PM(K+1)
        PUK=0.
        PVK=0.
        PT4K=0.
        PZ4K=0.
        PQ4K=0.
        DUTK=0.
        DVTK=0.
        PS4K=0.
C**** FOR AMIP
        AMVQ=0.
        AMVT=0.
        AMUU=0.
        AMVV=0.
        AMUV=0.
        AMTT=0.
C**** END AMIP
C**** INTERPOLATE HERE
  330 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      PUK=PUK+DP*U(I,J,L)
      PVK=PVK+DP*V(I,J,L)
      PT4K=PT4K+DP*(TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L))
      PZ4K=PZ4K+DP*(PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L))
      PQ4K=PQ4K+DP*(Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L))
      DUTK=DUTK+DP*DUT(I,J,L)
      DVTK=DVTK+DP*DVT(I,J,L)
      PS4K=PS4K+DP*(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))
C**** FOR AMIP 2
      AMRHT=.25*(TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L))
      AMRHQ=.25*(Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L))
      AMVQ=AMVQ+DP*V(I,J,L)*AMRHQ
      AMVT=AMVT+DP*V(I,J,L)*AMRHT
      AMUU=AMUU+DP*U(I,J,L)*U(I,J,L)
      AMVV=AMVV+DP*V(I,J,L)*V(I,J,L)
      AMUV=AMUV+DP*U(I,J,L)*V(I,J,L)
      AMTT=AMTT+DP*AMRHT*AMRHT
C**** END AMIP
      IF (LUP.EQ.L) GO TO 332
      L=L+1
      PDN=PL(L)
      GO TO 330
C**** ACCUMULATE HERE
  332 FIMI=FIMI+1.
      DPI=DPI+DPK
      DPSQI=DPSQI+DPK*DPK
      IF (DPK.LT.teeny) DPK=teeny
      BYDP=1./DPK
      PUI=PUI+PUK
      PVI=PVI+PVK
      PWWI=PWWI+BYDP*(PUK*PUK+PVK*PVK)
      PWWVI=PWWVI+BYDP*BYDP*(PUK*PUK+PVK*PVK)*PVK
      PUVI=PUVI+BYDP*PUK*PVK
      PT4I=PT4I+PT4K
      PTV4I=PTV4I+BYDP*PT4K*PVK
      PZ4I=PZ4I+PZ4K
      PZV4I=PZV4I+BYDP*PZ4K*PVK
      PQ4I=PQ4I+PQ4K
      PQV4I=PQV4I+BYDP*PQ4K*PVK
      DVTI=DVTI+DVTK
      VDVTI=VDVTI+BYDP*PVK*DVTK
      DUTI=DUTI+DUTK
      UDUTI=UDUTI+BYDP*PUK*DUTK
!!    IF(SKIPSE.EQ.1.) GO TO 334
      AIJK(I,J,K,IJK_U)  =AIJK(I,J,K,IJK_U)  +PUK
      AIJK(I,J,K,IJK_V)  =AIJK(I,J,K,IJK_V)  +PVK
      AIJK(I,J,K,IJK_DSE)=AIJK(I,J,K,IJK_DSE)+SHA*PT4K+PZ4K
      AIJK(I,J,K,IJK_DP) =AIJK(I,J,K,IJK_DP) +DPK
      AIJK(I,J,K,IJK_T)  =AIJK(I,J,K,IJK_T)  +PT4K
      AIJK(I,J,K,IJK_Q)  =AIJK(I,J,K,IJK_Q)  +PQ4K
      AIJK(I,J,K,IJK_R)  =AIJK(I,J,K,IJK_R)  +
     *     PQ4K/QSAT(0.25*PT4K/DPK,LHE,PMO(K))
      AIJK(I,J,K,IJK_PF)  =AIJK(I,J,K,IJK_PF)+1.
C     *  *  *  FOR AMIP 2  *  *  *
      AIJK(I,J,K,IJK_UV)=AIJK(I,J,K,IJK_UV)+AMUV
      AIJK(I,J,K,IJK_VQ)=AIJK(I,J,K,IJK_VQ)+AMVQ
      AIJK(I,J,K,IJK_VT)=AIJK(I,J,K,IJK_VT)+AMVT
      AIJK(I,J,K,IJK_UU)=AIJK(I,J,K,IJK_UU)+AMUU
      AIJK(I,J,K,IJK_VV)=AIJK(I,J,K,IJK_VV)+AMVV
      AIJK(I,J,K,IJK_TT)=AIJK(I,J,K,IJK_TT)+AMTT
C**** END AMIP
C**** EDDY TRANSPORT OF THETA;  VORTICITY
  334   PS4I=PS4I+PS4K
        PSV4I=PSV4I+BYDP*PVK*PS4K
        UDX(I,J,K)=BYDP*PUK*DXV(J)
!ESMF   IF (UDX(I,J-1,K).LT.BIG) ZX(I,J-1,K)=UDX(I,J,K)-UDX(I,J-1,K)
!ESMF   IF (UDX(I,J-1,K).GE.BIG) ZX(I,J-1,K)=0.
!ESMF   IF (ZX(I,J-1,K).GE.BIG) ZX(I,J-1,K)=0.
      END IF                            !---> (PM(K+1).GE.PS)

  340 I=IP1     !-->END I Loop (IP1) from IM to IM-1 (1-IM).
      DPM(K)=DPI/(FIMI+teeny)
      DPJK(J,K)=DPI
      AJK(J,K,JK_DPB)=AJK(J,K,JK_DPB)+DPI
      AJK(J,K,JK_DPSQR)=AJK(J,K,JK_DPSQR)+DPSQI
      AJK(J,K,JK_NPTSAVG)=AJK(J,K,JK_NPTSAVG)+FIMI
      IF (DPI.LT.teeny) DPI=teeny
      AJK(J,K,JK_U)=AJK(J,K,JK_U)+PUI
      AJK(J,K,JK_V)=AJK(J,K,JK_V)+PVI
      AJK(J,K,JK_ZMFKE)=AJK(J,K,JK_ZMFKE)+(PUI*PUI+PVI*PVI)/DPI
      AJK(J,K,JK_TOTKE)=AJK(J,K,JK_TOTKE)+PWWI
      AJK(J,K,JK_ZMFNTSH)=AJK(J,K,JK_ZMFNTSH)+PT4I*PVI/DPI
      AJK(J,K,JK_TOTNTSH)=AJK(J,K,JK_TOTNTSH)+PTV4I
      AJK(J,K,JK_ZMFNTGEO)=AJK(J,K,JK_ZMFNTGEO)+PZ4I*PVI/DPI
      AJK(J,K,JK_TOTNTGEO)=AJK(J,K,JK_TOTNTGEO)+PZV4I
      AJK(J,K,JK_ZMFNTLH)=AJK(J,K,JK_ZMFNTLH)+PQ4I*PVI/DPI
      AJK(J,K,JK_TOTNTLH)=AJK(J,K,JK_TOTNTLH)+PQV4I
      AJK(J,K,JK_ZMFNTKE)=AJK(J,K,JK_ZMFNTKE)+PWWI*PVI/DPI
      AJK(J,K,JK_TOTNTKE)=AJK(J,K,JK_TOTNTKE)+PWWVI
      AJK(J,K,JK_ZMFNTMOM)=AJK(J,K,JK_ZMFNTMOM)+PUI*PVI/DPI
      AJK(J,K,JK_TOTNTMOM)=AJK(J,K,JK_TOTNTMOM)+PUVI
      AJK(J,K,JK_P2KEDPGF)=AJK(J,K,JK_P2KEDPGF)+VDVTI+UDUTI-
     *   (PUI*DUTI+PVI*DVTI)/DPI
      SHETH(K)=(PSV4I-PS4I*PVI/DPI)*DXYV(J)/(STJK(J-1,K)*DXYN(J-1)+
     *   STJK(J,K)*DXYS(J))
         UJK(J,K)=PUI/DPI
         VJK(J,K)=PVI/DPI
         PSIJK(J,K)=.25*SHETH(K)/DPI
         UVJK(J,K)=(PUVI-PUI*PVI/DPI)/DPI
         IF (IDACC(4).EQ.1) AJK(J,K,JK_UINST)=UJK(J,K)
         AJK(J,K,JK_TOTDUDT)=UJK(J,K)-AJK(J,K,JK_UINST)
  350 AJK(J,K,JK_SHETH)=AJK(J,K,JK_SHETH)+SHETH(K)
  390 CONTINUE

C**** ZX for distributed parallelization
c****
      CALL HALO_UPDATE( grid, UDX, from=NORTH )
      CALL HALO_UPDATE( grid, pm_ge_ps, from=NORTH)

      DO J=J_0,J_1S
        DO K=1,KM
          DO I=1,IM
            if (pm_ge_ps(i,j+1,k) < 0) then
            IF (UDX(I,J,K).LT.BIG ) ZX(I,J,K)=-UDX(I,J,K)+UDX(I,J+1,K)
            IF (UDX(I,J,K).GE.BIG)  ZX(I,J,K)=0.
            IF (ZX(I,J,K).GE.BIG)   ZX(I,J,K)=0
            end if
          END DO
        END DO
      END DO
C****
C**** alternate vertical mass flux diagnostic (from SD)
C****
      DO J=J_0,J_1
        W(:,J,:)=0.
      END DO
C**** interpolate SD to constant pressure
      DO J=J_0,J_1
        I=IM
        DO IP1=1,IM
          DO K=1,KM-1
            DPK=0.
            SDK=0.
            SP=P(I,J)
            DO L=1,LS1-1
              PL(L)=SP*SIGE(L)+PTOP
            END DO
            IF (PM(K+1).GE.SP+PTOP) GO TO 860
            L=1
            PDN=SP+PTOP
            IF (PM(K).GE.SP+PTOP) GO TO 820
            PDN=PM(K)
 810        IF (PM(K).GT.PL(L+1)) GO TO 820
            L=L+1
            GO TO 810
 820        LUP=L
 830        IF (PM(K+1).GE.PL(LUP+1)) GO TO 840
            LUP=LUP+1
            GO TO 830
 840        CONTINUE
C**** INTERPOLATE HERE
 850        PUP=PL(L+1)
            IF (LUP.EQ.L) PUP=PM(K+1)
            DPK=DPK+(PDN-PUP)
            SDK=SDK+(PDN-PUP)*SD(I,JJ(J),L)
            IF (LUP.EQ.L) GO TO 860
            L=L+1
            PDN=PL(L)
            GO TO 850
 860        CONTINUE
C**** ACCUMULATE HERE (SHOULD I ACCUMULATE A WEIGHTING FUNCTION?)
            W(I,J,K)=0.
            IF (DPK.gt.0) THEN
              AIJK(I,J,K,IJK_W)=AIJK(I,J,K,IJK_W)+SDK*BYDXYP(J)/DPK
              W(I,J,K)=SDK/DPK
            END IF
          END DO
          I=IP1
        END DO
      END DO

C**** ACCUMULATE ALL VERTICAL WINDS
!!    DO 558 J=J_0,J_1
!!    DO 558 I=1,IM
!!    DO KR=1,NDIUPT
!!       IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
!!*** Warning:     This diagnostic has 3 flaws   (?)
!!***          1 - It assumes that DTsrc=1hr, (DTsrc=3600.)
!!***          2 - since DTdaa-Ndaa*DTsrc=2*DTdyn rather than 0,
!!***              some hours are skipped once in a while
!!***          3 - Some of the first Ndaa hours are skipped at the
!!***              beginning of a month and overcounted at the end;
!!***              this happens to balance out, if and only if
!!***              mod(days_in_month,ndaa)=0  (i.e. February if Ndaa=7)
!!***          In addition, IHM occasionally is out-of-bounds.
!!          IH=JHOUR+1
!!          IHM = IH+(JDATE-1)*24
!!          DO INCH=1,NDAA
!!            IF(IH.GT.HR_IN_DAY) IH=IH-HR_IN_DAY
!!            ADIURN(IH,IDD_W,KR)=ADIURN(IH,IDD_W,KR)+1.E5*W(I,J,3)
!!   *             /DXYP(J)
!!            HDIURN(IHM,IDD_W,KR)=HDIURN(IHM,IDD_W,KR)+1.E5*W(I,J,3)
!!   *             /DXYP(J)
!!            IH=IH+1
!!            IHM=IHM+1
!!          END DO
!!       END IF
!!    END DO
!!558 CONTINUE

      DO 565 J=J_0,J_1
      DO 565 K=1,KM
      WI=0.
      DO I=1,IMAXJ(J)
        WI=WI+W(I,J,K)
      END DO
  565 AJK(J,K,JK_VVEL)=AJK(J,K,JK_VVEL)+WI
C****
C**** ACCUMULATE T,Z,Q VERTICAL TRANSPORTS
C****
      DO 610 J=J_0,J_1
      DO 610 K=2,KM
      WI=0.
      TKI=0.
      QKI=0.
      ZKI=0.
      WTI=0.
      WQI=0.
      WZI=0.
         THKI=0.
         WTHI=0.
      FIMI=0.
      DO 600 I=1,IMAXJ(J)
      SP=P(I,J)
      DO 569 L=1,LS1-1
  569 PLO(L)=SP*SIG(L)+PTOP
      IF (PM(K).GE.SP+PTOP) GO TO 600
      L=1
      IF (PM(K).GE.PLO(1)) GO TO 580
  570 LUP=L+1
      IF (L.EQ.LM) GO TO 580
      IF (PM(K).GE.PLO(LUP)) GO TO 575
      L=LUP
      GO TO 570
  575 DPUP=PM(K)-PLO(LUP)
      DPDN=PLO(L)-PM(K)
      BYDP=1./(DPDN+DPUP)
      TK=BYDP*(TX(I,J,L)*DPUP+TX(I,J,LUP)*DPDN)
      QK=Q(I,J,L)*Q(I,J,LUP)/(BYDP*(Q(I,J,L)*DPDN+Q(I,J,LUP)*DPUP)+
     *  teeny)
      ZK=BYDP*(PHI(I,J,L)*DPUP+PHI(I,J,LUP)*DPDN)
         THK=BYDP*(T(I,J,L)*DPUP+T(I,J,LUP)*DPDN)
      GO TO 590
C**** SPECIAL CASES;  L=1, L=LM
  580 TK=TX(I,J,L)
      QK=Q(I,J,L)
      ZK=PHI(I,J,L)
         THK=T(I,J,L)
C**** MERIDIONAL AVERAGING
  590 WI=WI+W(I,J,K)
      TKI=TKI+TK
      QKI=QKI+QK
      ZKI=ZKI+ZK
      WTI=WTI+W(I,J,K)*TK
      WQI=WQI+W(I,J,K)*QK
      WZI=WZI+W(I,J,K)*ZK
         THKI=THKI+THK
         WTHI=WTHI+W(I,J,K)*THK
      FIMI=FIMI+1.
  600 CONTINUE
      BYFIM=teeny
      IF (FIMI.GT.teeny) BYFIM=1./FIMI
      AJK(J,K-1,JK_ZMFVTDSE)=AJK(J,K-1,JK_ZMFVTDSE)+
     &     BYFIM*(SHA*TKI+ZKI)*WI
      AJK(J,K-1,JK_TOTVTDSE)=AJK(J,K-1,JK_TOTVTDSE)+SHA*WTI+WZI
      AJK(J,K-1,JK_ZMFVTLH)=AJK(J,K-1,JK_ZMFVTLH)+BYFIM*QKI*WI
      AJK(J,K-1,JK_TOTVTLH)=AJK(J,K-1,JK_TOTVTLH)+WQI
      AJK(J,K-1,JK_VTGEOEDDY)=AJK(J,K-1,JK_VTGEOEDDY)+WZI-BYFIM*WI*ZKI
C     AJK(J,K-1,JK_BAREKEGEN)=AJK(J,K-1,JK_BAREKEGEN)+WTI-BYFIM*WI*TKI
         WJK(J,K)=BYFIM*WI/DXYP(J)
         WTJK(J,K)=BYFIM*(WTHI-BYFIM*WI*THKI)/DXYP(J)
         AJK(J,K-1,JK_EDDVTPT)=AJK(J,K-1,JK_EDDVTPT)+WTJK(J,K)
  610 CONTINUE
C****
C**** BAROCLINIC EDDY KINETIC ENERGY GENERATION
C****
      DO 630 J=J_0,J_1
      DO 630 K=1,KM
      FIMI=0.
      W2I=0.
      PAI=0.
      WPA2I=0.
      DO 626 I=1,IMAXJ(J)
      SP=P(I,J)
      DO 611 L=1,LS1-1
  611 PL(L)=SP*SIGE(L)+PTOP
      PS=SP+PTOP
      IF (PM(K+1).GE.PS) GO TO 626
      L=1
      PDN=PS
      IF (PM(K).GE.PS) GO TO 614
      PDN=PM(K)
  612 IF (PM(K).GT.PL(L+1)) GO TO 614
      L=L+1
      GO TO 612
  614 LUP=L
  616 IF (PM(K+1).GE.PL(LUP+1)) GO TO 618
      LUP=LUP+1
      GO TO 616
  618 CONTINUE
      PTK=0.
C**** INTERPOLATE HERE
  620 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      PTK=PTK+DP*TX(I,J,L)
      IF (LUP.EQ.L) GO TO 622
      L=L+1
      PDN=PL(L)
      GO TO 620
C**** ACCUMULATE HERE
  622 FIMI=FIMI+1.
      WUP=0.
      IF (K.LT.KM) WUP=W(I,J,K+1)
      W2I=W2I+W(I,J,K)+WUP
      PY=PMO(K)
      IF (PM(K).GE.PS) PY=.5*(PS+PM(K+1))
      PAK=PTK/PY
      PAI=PAI+PAK
      WPA2I=WPA2I+(W(I,J,K)+WUP)*PAK
  626 CONTINUE
  630 AJK(J,K,JK_BAREKEGEN)=AJK(J,K,JK_BAREKEGEN)-
     &     (WPA2I-W2I*PAI/(FIMI+teeny))
C****
C**** ACCUMULATE UV VERTICAL TRANSPORTS
C****
C**** DOUBLE POLAR WINDS
      DO 640 K=1,KM
        IF(HAVE_SOUTH_POLE) THEN
          WSP=2.*W(1,1,K)/FIM
          DO I=1,IM
            W(I,1,K)=WSP
          ENDDO
        ENDIF
        IF(HAVE_NORTH_POLE) THEN
          WNP=2.*W(1,JM,K)/FIM
          DO I=1,IM
            W(I,JM,K)=WNP
          ENDDO
        ENDIF
  640 CONTINUE

C P already halo'ed; no need      CALL CHECKSUM(grid, P, __LINE__, __FILE__)
C P already halo'ed; no need     CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      CALL HALO_UPDATE(grid, W, FROM=SOUTH)

      DO 710 J=J_0STG,J_1STG
      UEARTH=RADIUS*OMEGA*COSV(J)
      I=IM
      DO 650 IP1=1,IM
      PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
  650 I=IP1
      DO 710 K=2,KM
      W4I=0.
      UKI=0.
      WU4I=0.
      WKE4I=0.
      FIMI=0.
      I=IM
      DO 700 IP1=1,IM
      SP=PSEC(I)
      DO 660 L=1,LS1-1
  660 PLO(L)=SP*SIG(L)+PTOP
      IF (PM(K).GE.SP+PTOP) GO TO 700
      L=1
      IF (PM(K).GE.PLO(1)) GO TO 680
  670 LUP=L+1
      IF (L.EQ.LM) GO TO 680
      IF (PM(K).GE.PLO(LUP)) GO TO 675
      L=LUP
      GO TO 670
  675 DPUP=PM(K)-PLO(LUP)
      DPDN=PLO(L)-PM(K)
      BYDP=1./(DPDN+DPUP)
      UK=BYDP*(U(I,J,L)*DPUP+U(I,J,LUP)*DPDN)
      VK=BYDP*(V(I,J,L)*DPUP+V(I,J,LUP)*DPDN)
      GO TO 690
C**** SPECIAL CASES;  L=1,L=LM
  680 UK=U(I,J,L)
      VK=V(I,J,L)
C**** MERIDIONAL AVERAGING
  690 W4=W(I,J-1,K)+W(IP1,J-1,K)+W(I,J,K)+W(IP1,J,K)
      W4I=W4I+W4
      UKI=UKI+UK
      WU4I=WU4I+W4*UK
      WKE4I=WKE4I+W4*(UK*UK+VK*VK)
      FIMI=FIMI+1.
  700 I=IP1
      BYFIM=1./(FIMI+teeny)
         WUJK(J,K)=.25*(WU4I-W4I*UKI*BYFIM)*BYFIM/DXYV(J)
      AJK(J,K-1,JK_TOTVTKE)=AJK(J,K-1,JK_TOTVTKE)+WKE4I
      AJK(J,K-1,JK_VTAMEDDY)=AJK(J,K-1,JK_VTAMEDDY)+WU4I-BYFIM*W4I*UKI
  710 AJK(J,K-1,JK_TOTVTAM)=AJK(J,K-1,JK_TOTVTAM)+WU4I   !+W4I*UEARTH
C****
C**** POTENTIAL VORTICITY AND VERTICAL TRANSPORT OF POT. VORT.
C****
      DO 760 J=J_0S,J_1S
      JHEMI=1
      IF (J.LT.1+JM/2) JHEMI=-1
      DO 730 K=1,KM
      PVI=0.
      DO 720 I=1,IM
      DUT(I,J,K)=JHEMI*STB(I,J,K)*(ZX(I,J,K)-FCOR(J))
  720 PVI=PVI+DUT(I,J,K)
  730 AJK(J,K,JK_POTVORT)=AJK(J,K,JK_POTVORT)+PVI
      DO 760 K=2,KM
      W2I=0.
      PV2I=0.
      WPV4I=0.
      FIMI=0.
      I=IM
      DO 740 IP1=1,IM
      PS=.5*(P(I,J)+P(IP1,J))+PTOP
      IF (PM(K).GE.PS) GO TO 740
      W2=W(I,J,K)+W(IP1,J,K)
      W2I=W2I+W2
      PV2=DUT(I,J,K-1)+DUT(I,J,K)
      PV2I=PV2I+PV2
      WPV4I=WPV4I+W2*PV2
      FIMI=FIMI+1.
  740 I=IP1
      AJK(J,K-1,JK_VTPV)=AJK(J,K-1,JK_VTPV)+WPV4I
  760 AJK(J,K-1,JK_VTPVEDDY)=AJK(J,K-1,JK_VTPVEDDY)+
     &     WPV4I-W2I*PV2I/(FIMI+teeny)
C****
C**** SPECIAL MEAN/EDDY DIAGNOSTICS ARE CALCULATED
C****
      DO 770 J=J_0STG,J_1STG
      DO 765 K=2,KM
      DPE=PMO(K)-PMO(K-1)
      UP(J,K)=(UJK(J,K)-UJK(J,K-1))/DPE
  765 PSIP(J,K)=(PSIJK(J,K)-PSIJK(J,K-1))/DPE
      UP(J,1)=UP(J,2)
      PSIP(J,1)=PSIP(J,2)
  770 CONTINUE
      DO 780 K=1,KM
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      KDN=K-1
      IF (K.EQ.1) KDN=1

      CALL HALO_UPDATEj(grid, TJK, FROM=SOUTH)
c***      DO L=1,LM
c***         CALL HALO_UPDATE(grid, TJK(:,L), FROM=SOUTH)
c***      END DO

      DO 780 J=J_0STG,J_1STG
      TY(J,K)=(TJK(J,K)-TJK(J-1,K))/DYV(J)
C**** E-P FLUX NORTHWARD COMPONENT
      AJK(J,K,JK_EPFLXNCP)=AJK(J,K,JK_EPFLXNCP)+
     &     PSIJK(J,K)*(UJK(J,KUP)-UJK(J,KDN))/
     *  (PMO(KUP)-PMO(KDN))-UVJK(J,K)
  780 CONTINUE

      CALL HALO_UPDATEj(grid, PSIJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, UJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, VJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, WJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, UP, FROM=NORTH)
      CALL HALO_UPDATEj(grid, TY, FROM=NORTH)
      CALL HALO_UPDATEj(grid, PSIP, FROM=NORTH)

c***      DO L=1,LM
c***         CALL HALO_UPDATE(grid, PSIJK(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, UJK(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, VJK(:,L), FROM=NORTH)
c***         If (L > 1) THEN
c***           CALL HALO_UPDATE(grid, WJK(:,L), FROM=NORTH)
c***         END IF
c***         CALL HALO_UPDATE(grid, UP(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, TY(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, PSIP(:,L), FROM=NORTH)
c***      END DO

      DO 800 J=J_0S,J_1S
      DO 800 K=2,KM-1
      UY=(UJK(J+1,K)*DXV(J+1)-UJK(J,K)*DXV(J)-FCOR(J))/DXYP(J)
      PSIY=(PSIJK(J+1,K)*DXV(J+1)-PSIJK(J,K)*DXV(J))/DXYP(J)
C**** ZONAL MEAN MOMENTUM EQUATION   (MEAN ADVECTION)
      AJK(J,K,JK_DUDTMADV)=AJK(J,K,JK_DUDTMADV)-
     &     .5*UY*(VJK(J,K)+VJK(J+1,K))-
     *  .25*((UP(J+1,K+1)+UP(J,K+1))*WJK(J,K+1)+(UP(J+1,K)+UP(J,K))*
     *   WJK(J,K))
C**** ZONAL MEAN HEAT EQUATION   (MEAN ADVECTION)
      AJK(J,K,JK_DTDTMADV)=AJK(J,K,JK_DTDTMADV)-
     &     .5*(TY(J,K)*VJK(J,K)+TY(J+1,K)*VJK(J+1,K))
     *  -.5*STJK(J,K)*(WJK(J,K+1)+WJK(J,K))
C**** LAGRANGIAN MEAN MOMENTUM EQUATION  (MEAN ADVECTION)
      VSTAR=.5*(VJK(J,K)+VJK(J+1,K)-.5*(PSIP(J,K)+PSIP(J,K+1)
     *  +PSIP(J+1,K)+PSIP(J+1,K+1)))
      WSTAR=.5*(WJK(J,K)+WJK(J,K+1))+PSIY
      AJK(J,K,JK_DUDTTEM)=AJK(J,K,JK_DUDTTEM)-
     &     UY*VSTAR-.25*(UP(J,K)+UP(J+1,K)+
     *  UP(J,K+1)+UP(J+1,K+1))*WSTAR
      AJK(J,K,JK_DTDTTEM)=AJK(J,K,JK_DTDTTEM)-
     &     .5*(TY(J+1,K)+TY(J,K))*VSTAR-
     *  STJK(J,K)*WSTAR
C**** VERTICAL E-P FLUX
      AJK(J,K-1,JK_EPFLXVCP)=AJK(J,K-1,JK_EPFLXVCP)-
     &     WUJK(J,K)-.5*PSIJK(J,K)*UY
      AJK(J,K,JK_EPFLXVCP)=AJK(J,K,JK_EPFLXVCP)-.5*PSIJK(J,K)*UY
  800 CONTINUE
C****
C**** SPECTRAL ANALYSIS OF KINETIC ENERGIES AT CONSTANT PRESSURE
C****
      IZERO=0
      NM=1+IM/2
      J45N=2+.75*(JM-1.)
c      KS1=LS1
C**** TOTAL THE KINETIC ENERGIES
      KE(:,:)=0.
      KE_part(:,:,:)=0.

C P already halo'ed; no need      CALL CHECKSUM(grid, P, __LINE__, __FILE__)
C P already halo'ed; no need     CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      DO J=J_0STG,J_1STG
        I=IM
        DO IP1=1,IM
          PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
          I=IP1
        ENDDO
        DO K=1,KM
          KSPHER=KLAYER(K)
          IF (J.GT.JEQ) KSPHER=KSPHER+1
          DO KX=IZERO,LM,LM
            DO I=1,IM
              DPUV=0.
              SP=PSEC(I)
              DO 2025 L=1,LS1-1
              PLO(L)=SP*SIG(L)+PTOP                       ! PL or PLO ??
 2025         PL(L)=SP*SIGE(L)+PTOP                       ! PLE or PL ??
              PS=SP+PTOP
              IF (PM(K+1).GE.PLO(1)) GO TO 2090           ! really ?? not PL?
              L=1
              PDN=PS
              IF (PM(K).GE.PLO(1)) GO TO 2040             ! really ?? not PL?
              PDN=PM(K)
 2030         IF (PM(K).GT.PL(L+1)) GO TO 2040
              L=L+1
              GO TO 2030
 2040         LUP=L
 2050         IF (PM(K+1).GE.PL(LUP+1)) GO TO 2060
              LUP=LUP+1
              GO TO 2050
 2060         CONTINUE
C**** ACCUMULATE HERE
              SQRTDP=SQRT(PDN-PM(K+1))
 2070         PUP=PL(L+1)
              IF (LUP.EQ.L) PUP=PM(K+1)
              DP=PDN-PUP
              IF(KX.EQ.IZERO) DPUV=DPUV+DP*U(I,J,L)
              IF(KX.EQ.LM)    DPUV=DPUV+DP*V(I,J,L)
              IF (LUP.EQ.L) GO TO 2080
              L=L+1
              PDN=PL(L)
              GO TO 2070
 2080         IF (SQRTDP.EQ.0.) SQRTDP=teeny
              DPUV=DPUV/SQRTDP
 2090         X1(I)=DPUV
            ENDDO
            CALL FFTE (X1,X1)
            IF (J.NE.JEQ) THEN
              DO N=1,NM
                KE_part(N,J,KSPHER)=KE_part(N,J,KSPHER)+X1(N)*DXYV(J)
              ENDDO
              IF (J.EQ.J45N) THEN
                DO N=1,NM
                  KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                  X1(N)*DXYV(J)
                ENDDO
              ENDIF
            ELSE
              DO N=1,NM
                KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                X1(N)*DXYV(J)
                KE_part(N,J,KSPHER  )=KE_part(N,J,KSPHER)+
     &                                .5D0*X1(N)*DXYV(J)
                KE_part(N,J,KSPHER+1)=KE_part(N,J,KSPHER+1)+
     &                                .5D0*X1(N)*DXYV(J)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL GLOBALSUM(grid, KE_part, KE, ALL=.true.)

      DO 2150 KS=1,NSPHER
      DO 2150 N=1,NM
 2150 SPECA(N,18,KS)=SPECA(N,18,KS)+KE(N,KS)
C**** ACCUMULATE TIME USED IN DIAGA
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAGB


      SUBROUTINE DIAG7A
C****
C**** THIS ROUTINE ACCUMULATES A TIME SEQUENCE FOR SELECTED
C**** QUANTITIES AND FROM THAT PRINTS A TABLE OF WAVE FREQUENCIES.
C****
      USE CONSTANT, only : grav,bygrav
      USE MODEL_COM, only : im,imh,jm,lm,
     &     IDACC,JEQ,LS1,MDIAG,P,PTOP,PSFMPT,SIG,SIGE,U,V
      USE DYNAMICS, only : PHI
      USE DIAG_COM, only : nwav_dag,wave,max12hr_sequ,j50n,kwp,re_and_im
      USE DIAG_LOC, only : ldex
      USE DOMAIN_DECOMP, only : GRID,GET,GLOBALSUM
      IMPLICIT NONE

      REAL*8, DIMENSION(0:IMH) :: AN,BN
      REAL*8, DIMENSION(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP,
     & grid%J_STRT_HALO:grid%J_STOP_HALO) :: WAVE_part
      INTEGER, PARAMETER :: KM=6,KQMAX=12
      INTEGER :: NMAX=nwav_dag
      REAL*8, DIMENSION(IM,KM) :: HTRD
      REAL*8, DIMENSION(KM), PARAMETER ::
     &     PMB=(/922.,700.,500.,300.,100.,10./),
     &     GHT=(/500.,2600.,5100.,8500.,15400.,30000./)
      REAL*8 :: PIJ50N,PL,PLE,PLM1,SLOPE
      INTEGER I,IDACC9,JLK,K,KQ,L,LX,MNOW,N
      INTEGER :: J_0, J_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

      IDACC9=IDACC(9)+1
      IDACC(9)=IDACC9
      IF (IDACC9.GT.Max12HR_sequ) RETURN

      WAVE_part(:,:,:,1:6,:)=0.
      IF(J_0 <= JEQ .and. JEQ <= J_1) THEN
        DO KQ=1,3
          CALL FFT (U(1,JEQ,LDEX(KQ)),AN,BN)
          DO N=1,NMAX
            WAVE_part(1,IDACC9,N,2*KQ-1,JEQ)=AN(N)
            WAVE_part(2,IDACC9,N,2*KQ-1,JEQ)=BN(N)
          ENDDO
          CALL FFT (V(1,JEQ,LDEX(KQ)),AN,BN)
          DO N=1,NMAX
            WAVE_part(1,IDACC9,N,2*KQ,JEQ)=AN(N)
            WAVE_part(2,IDACC9,N,2*KQ,JEQ)=BN(N)
          ENDDO
        ENDDO
      ENDIF

      CALL GLOBALSUM(grid, WAVE_part(1,IDACC9:IDACC9,1:NMAX,1:6,:),
     &    WAVE(1,IDACC9:IDACC9,1:NMAX,1:6), ALL=.TRUE.)
      CALL GLOBALSUM(grid, WAVE_part(2,IDACC9:IDACC9,1:NMAX,1:6,:),
     &    WAVE(2,IDACC9:IDACC9,1:NMAX,1:6), ALL=.TRUE.)

      WAVE_part(:,:,:,7:,:)=0.
      IF(J_0 <= J50N .and. J50N <= J_1) THEN
        DO 150 I=1,IM
          PIJ50N=P(I,J50N)
          K=1
          L=1
          PL=SIG(1)*P(I,J50N)+PTOP
 130      L=L+1
          IF(L.GE.LS1) PIJ50N=PSFMPT
          PLM1=PL
          PL=SIG(L)*PIJ50N+PTOP
          IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
C**** ASSUME THAT PHI IS LINEAR IN LOG P
          SLOPE=(PHI(I,J50N,L-1)-PHI(I,J50N,L))/LOG(PLM1/PL)
 140      HTRD(I,K)=(PHI(I,J50N,L)+SLOPE*LOG(PMB(K)/PL))*BYGRAV-GHT(K)
          IF (K.GE.KM) GO TO 150
          K=K+1
          IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
          GO TO 140
 150    CONTINUE
        DO KQ=7,KQMAX
          CALL FFT(HTRD(1,KQ-6),AN,BN)
          DO N=1,NMAX
            WAVE_part(1,IDACC9,N,KQ,J50N)=AN(N)
            WAVE_part(2,IDACC9,N,KQ,J50N)=BN(N)
          END DO
        END DO
      ENDIF
      CALL GLOBALSUM(grid, WAVE_part(1,IDACC9:IDACC9,1:NMAX,7:KQMAX,:),
     &    WAVE(1,IDACC9:IDACC9,1:NMAX,7:KQMAX), ALL=.TRUE.)
      CALL GLOBALSUM(grid, WAVE_part(2,IDACC9:IDACC9,1:NMAX,7:KQMAX,:),
     &    WAVE(2,IDACC9:IDACC9,1:NMAX,7:KQMAX), ALL=.TRUE.)

      CALL TIMER (MNOW,MDIAG)
      RETURN
      END SUBROUTINE DIAG7A


      SUBROUTINE DIAGCA (M)
!@sum  DIAGCA Keeps track of the conservation properties of angular
!@+    momentum, kinetic energy, mass, total potential energy and water
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : mdiag,itime
#ifdef TRACERS_ON
      USE TRACER_COM, only: itime_tr0,ntm  !xcon
#endif
      USE DIAG_COM, only : icon_AM,icon_KE,icon_MS,icon_TPE
     *     ,icon_WM,icon_LKM,icon_LKE,icon_EWM,icon_WTG,icon_HTG
     *     ,icon_OMSI,icon_OHSI,icon_OSSI,icon_LMSI,icon_LHSI,icon_MLI
     *     ,icon_HLI,title_con
      USE SOIL_DRV, only: conserv_WTG,conserv_HTG
      IMPLICIT NONE
!@var M index denoting from where DIAGCA is called
      INTEGER, INTENT(IN) :: M
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1  INITIALIZE CURRENT QUANTITY
C****   2  AFTER DYNAMICS
C****   3  AFTER CONDENSATION
C****   4  AFTER RADIATION
C****   5  AFTER PRECIPITATION
C****   6  AFTER LAND SURFACE (INCL. RIVER RUNOFF)
C****   7  AFTER FULL SURFACE INTERACTION
C****   8  AFTER FILTER
C****   9  AFTER OCEAN DYNAMICS (from MAIN)
C****  10  AFTER DAILY
C****  11  AFTER OCEAN DYNAMICS (from ODYNAM)
C****  12  AFTER OCEAN SUB-GRIDSCALE PHYS
C****
      EXTERNAL conserv_AM,conserv_KE,conserv_MS,conserv_PE
     *     ,conserv_WM,conserv_EWM,conserv_LKM,conserv_LKE,conserv_OMSI
     *     ,conserv_OHSI,conserv_OSSI,conserv_LMSI,conserv_LHSI
     *     ,conserv_MLI,conserv_HLI
      INTEGER MNOW
      INTEGER NT

C**** ATMOSPHERIC ANGULAR MOMENTUM
      CALL conserv_DIAG(M,conserv_AM,icon_AM)

C**** ATMOSPHERIC KINETIC ENERGY
      CALL conserv_DIAG(M,conserv_KE,icon_KE)

C**** ATMOSPHERIC MASS
      CALL conserv_DIAG(M,conserv_MS,icon_MS)

C**** ATMOSPHERIC TOTAL POTENTIAL ENERGY
      CALL conserv_DIAG(M,conserv_PE,icon_TPE)

C**** ATMOSPHERIC TOTAL WATER MASS
      CALL conserv_DIAG(M,conserv_WM,icon_WM)

C**** ATMOSPHERIC TOTAL WATER ENERGY
      CALL conserv_DIAG(M,conserv_EWM,icon_EWM)

C**** LAKE MASS AND ENERGY
      CALL conserv_DIAG(M,conserv_LKM,icon_LKM)
      CALL conserv_DIAG(M,conserv_LKE,icon_LKE)

C**** OCEAN ICE MASS, ENERGY, SALT
      CALL conserv_DIAG(M,conserv_OMSI,icon_OMSI)
      CALL conserv_DIAG(M,conserv_OHSI,icon_OHSI)
      CALL conserv_DIAG(M,conserv_OSSI,icon_OSSI)

C**** LAKE ICE MASS, ENERGY
      CALL conserv_DIAG(M,conserv_LMSI,icon_LMSI)
      CALL conserv_DIAG(M,conserv_LHSI,icon_LHSI)

C**** GROUND WATER AND ENERGY
      CALL conserv_DIAG(M,conserv_WTG,icon_WTG)
      CALL conserv_DIAG(M,conserv_HTG,icon_HTG)

C**** LAND ICE MASS AND ENERGY
      CALL conserv_DIAG(M,conserv_MLI,icon_MLI)
      CALL conserv_DIAG(M,conserv_HLI,icon_HLI)

C**** OCEAN CALLS ARE DEALT WITH SEPARATELY
      CALL DIAGCO (M)

#ifdef TRACERS_ON
C**** Tracer calls are dealt with separately
      do nt=1,ntm
        CALL DIAGTCA(M,NT)
      end do
#endif
C****
      CALL TIMER (MNOW,MDIAG)
      RETURN
      END SUBROUTINE DIAGCA


      module DIAG
      contains

      SUBROUTINE DIAGCD (grid,M,UX,VX,DUT,DVT,DT1,PIT)
!@sum  DIAGCD Keeps track of the conservation properties of angular
!@+    momentum and kinetic energy inside dynamics routines
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : omega,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,mdiag,mdyn
      USE GEOM, only : cosv,radius,ravpn,ravps
      USE DIAG_COM, only : consrv=>consrv_loc
      USE DOMAIN_DECOMP, only : GET, CHECKSUM, HALO_UPDATE, DIST_GRID
      USE DOMAIN_DECOMP, only : SOUTH
      IMPLICIT NONE
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCD IS BEING CALLED
C**** M=1  AFTER ADVECTION IN DYNAMICS
C****   2  AFTER CORIOLIS FORCE IN DYNAMICS
C****   3  AFTER PRESSURE GRADIENT FORCE IN DYNAMICS
C****   4  AFTER STRATOS DRAG IN DYNAMICS
C****   5  AFTER FLTRUV IN DYNAMICS
C****   6  AFTER GRAVITY WAVE DRAG IN DYNAMICS
C****
      TYPE (DIST_GRID), INTENT(IN) :: grid
!@var M index denoting from where DIAGCD is called
      INTEGER, INTENT(IN) :: M
!@var DT1 current time step
      REAL*8, INTENT(IN) :: DT1
!@var UX,VX current velocities
      REAL*8, INTENT(IN),
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        UX,VX
!@var DUT,DVT current momentum changes
      REAL*8, INTENT(IN),
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        DUT,DVT
!@var PIT current pressure tendency
      REAL*8, INTENT(IN), OPTIONAL,
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PIT
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PI
      INTEGER :: I,J,L,MBEGIN,N,IP1
      LOGICAL dopit
      REAL*8 :: DUTI,DUTIL,RKEI,RKEIL
      INTEGER, DIMENSION(6) ::
     *     NAMOFM=(/2,3,4,5,6,7/), NKEOFM=(/14,15,16,17,18,19/)

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(MBEGIN)

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

C****
C**** PRESSURE TENDENCY FOR CHANGE BY ADVECTION
C****
      IF (M.eq.1) THEN
        dopit=.true.
        IF(HAVE_SOUTH_POLE) PI(1)=FIM*PIT(1,1)
        IF(HAVE_NORTH_POLE) PI(JM)=FIM*PIT(1,JM)
        DO J=J_0S,J_1S
          PI(J)=SUM(PIT(:,J))
        END DO
      ELSE
        PI=0.
        dopit=.false.
      END IF
C****
C**** CHANGE OF ANGULAR MOMENTUM AND KINETIC ENERGY BY VARIOUS
C**** PROCESSES IN DYNAMICS
C****
C****

      CALL HALO_UPDATE(grid, PI, FROM=SOUTH)

!$OMP PARALLEL DO PRIVATE (J,L,I,DUTIL,RKEIL,DUTI,RKEI,N)
      DO J=J_0STG,J_1STG
        DUTIL=0.
        RKEIL=0.
        DO L=1,LM
          DUTI=0.
          RKEI=0.
          DO I=1,IM
            DUTI=DUTI+DUT(I,J,L)
            RKEI=RKEI+(UX(I,J,L)*DUT(I,J,L)+VX(I,J,L)*DVT(I,J,L))
          END DO
          DUTIL=DUTIL+DUTI
          RKEIL=RKEIL+RKEI
        END DO
        N=NAMOFM(M)
        if (dopit) DUTIL=DUTIL+2.*DT1*RADIUS*OMEGA*COSV(J)*
     *       (PI(J-1)*RAVPN(J-1)+PI(J)*RAVPS(J))
        CONSRV(J,N)=CONSRV(J,N)+DUTIL*COSV(J)*RADIUS*mb2kg
        N=NKEOFM(M)
        CONSRV(J,N)=CONSRV(J,N)+RKEIL*mb2kg
      END DO
!$OMP END PARALLEL DO
C****
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAGCD

      end module DIAG


      SUBROUTINE conserv_DIAG (M,CONSFN,ICON)
!@sum  conserv_DIAG generic routine keeps track of conserved properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : jm
      USE DOMAIN_DECOMP, only : GET, GRID, CHECKSUMj
      USE DIAG_COM, only : consrv=>consrv_loc,nofm
      IMPLICIT NONE
!@var M index denoting from where routine is called
      INTEGER, INTENT(IN) :: M
!@var ICON index for the quantity concerned
      INTEGER, INTENT(IN) :: ICON
!@var CONSFN external routine that calculates total conserved quantity
      EXTERNAL CONSFN
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TOTAL
      INTEGER :: I,J,NM,NI
      INTEGER :: J_0,J_1

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1)

C**** NOFM contains the indexes of the CONSRV array where each
C**** change is to be stored for each quantity. If NOFM(M,ICON)=0,
C**** no calculation is done.
C**** NOFM(1,ICON) is the index for the instantaneous value.
      IF (NOFM(M,ICON).gt.0) THEN
C**** Calculate current value TOTAL
        CALL CONSFN(TOTAL)
        NM=NOFM(M,ICON)
        NI=NOFM(1,ICON)
C**** Accumulate difference from last time in CONSRV(NM)
        IF (M.GT.1) THEN
          DO J=J_0,J_1
            CONSRV(J,NM)=CONSRV(J,NM)+(TOTAL(J)-CONSRV(J,NI))
          END DO
        END IF
C**** Save current value in CONSRV(NI)
        DO J=J_0,J_1
          CONSRV(J,NI)=TOTAL(J)
        END DO
      END IF
      RETURN
C****
      END SUBROUTINE conserv_DIAG


      SUBROUTINE conserv_AM(AM)
!@sum  conserv_AM calculates total atmospheric angular momentum
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : omega,radius,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,ls1,dsig,p,u,psfmpt,pstrat
      USE GEOM, only : cosv,dxyn,dxys,dxyv
      USE DOMAIN_DECOMP, only : GET, SOUTH, HALO_UPDATE, GRID
      USE DOMAIN_DECOMP, only : CHECKSUM
      IMPLICIT NONE

      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM,PI
      INTEGER :: I,IP1,J,L
      REAL*8 :: UMI,UMIL

      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

C****
C**** ANGULAR MOMENTUM
C****
      IF(HAVE_SOUTH_POLE) PI(1)=FIM*P(1,1)
      IF(HAVE_SOUTH_POLE) AM(1)=0.
      IF(HAVE_NORTH_POLE) PI(JM)=FIM*P(1,JM)
      DO J=J_0S,J_1S
        PI(J)=0.
        DO I=1,IM
          PI(J)=PI(J)+P(I,J)
        END DO
      END DO

      CALL HALO_UPDATE(grid, PI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      DO J=J_0STG,J_1STG
        UMIL=0.
        DO L=1,LM
          UMI=0.
          I=IM
          DO IP1=1,IM
            IF(L.LT.LS1) THEN
              UMI=UMI+U(I,J,L)*((P(I,J-1)+P(IP1,J-1))*DXYN(J-1)
     *             +(P(I,J)+P(IP1,J))*DXYS(J))
            ELSE
              UMI=UMI+U(I,J,L)*(2.*PSFMPT*DXYV(J))
            END IF
            I=IP1
          END DO
          UMIL=UMIL+UMI*DSIG(L)
        END DO
        AM(J)=(RADIUS*OMEGA*COSV(J)*((PI(J-1)*DXYN(J-1)+PI(J)*DXYS(J))
     *       +FIM*PSTRAT*DXYV(J))+.5*UMIL)*COSV(J)*RADIUS*mb2kg
      END DO
      RETURN
C****
      END SUBROUTINE conserv_AM


      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates total atmospheric kinetic energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,dsig,ls1,p,u,v,psfmpt
      USE GEOM, only : dxyn,dxys,dxyv
      USE DOMAIN_DECOMP, only : GET, CHECKSUM, HALO_UPDATE, GRID
      USE DOMAIN_DECOMP, only : SOUTH
      IMPLICIT NONE

      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RKE
      INTEGER :: I,IP1,J,L
      INTEGER :: J_0STG,J_1STG
      REAL*8 :: RKEI,RKEIL
      LOGICAL :: HAVE_SOUTH_POLE

      CALL GET(grid, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

C****
C**** KINETIC ENERGY
C****

      CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      IF (HAVE_SOUTH_POLE) RKE(1)=0.
      DO J=J_0STG,J_1STG
        RKEIL=0.
        DO L=1,LM
          RKEI=0.
          I=IM
          DO IP1=1,IM
            IF(L.LT.LS1) THEN
              RKEI=RKEI+(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
     *             *((P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J
     *             ))*DXYS(J))
            ELSE
              RKEI=RKEI+(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*
     *             (2.*PSFMPT*DXYV(J))
            END IF
            I=IP1
          END DO
          RKEIL=RKEIL+RKEI*DSIG(L)
        END DO
        RKE(J)=0.25*RKEIL*mb2kg
      END DO
      RETURN
C****
      END SUBROUTINE conserv_KE


      SUBROUTINE conserv_MS(RMASS)
!@sum  conserv_MA calculates total atmospheric mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,fim,p,pstrat
      USE DOMAIN_DECOMP, only : GET, GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RMASS
      INTEGER :: I,J
      INTEGER :: J_0S,J_1S
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grid, J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)


C****
C**** MASS
C****
      IF(HAVE_SOUTH_POLE) RMASS(1) =FIM*(P(1,1) +PSTRAT)*mb2kg
      IF(HAVE_NORTH_POLE) RMASS(JM)=FIM*(P(1,JM)+PSTRAT)*mb2kg
      DO J=J_0S,J_1S
        RMASS(J)=FIM*PSTRAT
        DO I=1,IM
          RMASS(J)=RMASS(J)+P(I,J)
        END DO
        RMASS(J)=RMASS(J)*mb2kg
      END DO
      RETURN
C****
      END SUBROUTINE conserv_MS


      SUBROUTINE conserv_PE(TPE)
!@sum  conserv_TPE calculates total atmospheric potential energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,t,p,ptop,zatmo
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pk,pdsig
      USE DOMAIN_DECOMP, only : GET,GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TPE
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      REAL*8 :: TPEI,TPEIL,SGEOI

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C****
C**** TOTAL POTENTIAL ENERGY (J/m^2)
C****
      DO J=J_0,J_1
        TPEIL=0.
        DO L=1,LM
          TPEI=0.
          DO I=1,IMAXJ(J)
            TPEI=TPEI+T(I,J,L)*PK(L,I,J)*PDSIG(L,I,J)
          END DO
          TPEIL=TPEIL+TPEI
        END DO
        SGEOI=0.
        DO I=1,IMAXJ(J)
          SGEOI=SGEOI+ZATMO(I,J)*(P(I,J)+PTOP)
        END DO
        TPE(J)=(SGEOI+TPEIL*SHA)*mb2kg
      END DO
      IF(HAVE_SOUTH_POLE) TPE(1)=FIM*TPE(1)
      IF(HAVE_NORTH_POLE) TPE(JM)=FIM*TPE(JM)
      RETURN
C****
      END SUBROUTINE conserv_PE

      SUBROUTINE conserv_WM(WATER)
!@sum  conserv_WM calculates total atmospheric water mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,wm,q
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig
      USE DOMAIN_DECOMP, only : GET, GRID
      IMPLICIT NONE

      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: WATER
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C****
C**** TOTAL WATER MASS (kg/m^2)
C****
      DO J=J_0,J_1
        WATER(J) = 0.
        DO L=1,LM
          DO I=1,IMAXJ(J)
            WATER(J)=WATER(J)+(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)*mb2kg
          END DO
        END DO
      END DO
      IF (HAVE_SOUTH_POLE) WATER(1) = FIM*WATER(1)
      IF (HAVE_NORTH_POLE) WATER(JM)= FIM*WATER(JM)
      RETURN
C****
      END SUBROUTINE conserv_WM


      SUBROUTINE conserv_EWM(EWATER)
!@sum  conserv_EWM calculates total atmospheric water energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg,shv,grav,lhe
      USE MODEL_COM, only : im,jm,lm,fim,wm,t,q,p
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig, pmid, pk
      USE CLOUDS_COM, only : svlhx
      USE DOMAIN_DECOMP, only : GET, GRID
      IMPLICIT NONE
      REAL*8, PARAMETER :: HSCALE = 7.8d0 ! km
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: EWATER
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      REAL*8 W,EL

      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

C****
C**** TOTAL WATER ENERGY (J/m^2)
C****
      DO J=J_0,J_1
        EWATER(J) = 0.
        DO L=1,LM
          DO I=1,IMAXJ(J)
c this calculation currently only calculates latent heat
            W =(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)*mb2kg
            EL=(Q(I,J,L)*LHE+WM(I,J,L)*(LHE-SVLHX(L,I,J)))*PDSIG(L,I,J)
     *           *mb2kg
            EWATER(J)=EWATER(J)+EL !+W*(SHV*T(I,J,L)*PK(L,I,J)+GRAV
!     *           *HSCALE*LOG(P(I,J)/PMID(L,I,J)))
          END DO
        END DO
      END DO
      IF(HAVE_SOUTH_POLE) EWATER(1) = FIM*EWATER(1)
      IF(HAVE_NORTH_POLE) EWATER(JM)= FIM*EWATER(JM)
      RETURN
C****
      END SUBROUTINE conserv_EWM


      SUBROUTINE DIAG5D (M5,NDT,DUT,DVT)
      USE MODEL_COM, only : im,imh,jm,lm,fim,
     &     DSIG,JEQ,LS1,MDIAG,MDYN
      USE DIAG_COM, only : speca,nspher,klayer
      USE DIAG_LOC, only : FCUVA,FCUVB
      USE DOMAIN_DECOMP, only : GRID,GET,GLOBALSUM, WRITE_PARALLEL
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        DUT,DVT

c      REAL*8, DIMENSION(0:IMH,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,2) :: FCUVA,FCUVB
c      COMMON/WORK7/FCUVA,FCUVB

      INTEGER :: M5,NDT

      REAL*8, DIMENSION(IMH+1) :: X
      REAL*8, DIMENSION(0:IMH) :: FA,FB
      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE
      REAL*8, DIMENSION
     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part

      INTEGER :: J,J45N,KUV,KSPHER,L,MBEGIN,MKE,N,NM
      INTEGER :: J_0STG,J_1STG

      CALL GET(GRID,J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG)

      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      MKE=M5

      GO TO (810,810,810,100,100,  100,810),M5
C****  810 WRITE (6,910) M5
  810 CALL WRITE_PARALLEL(M5, UNIT=6, format=
     & "('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)")
C****  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)
      call stop_model('INCORRECT VALUE OF M5 WHEN CALLING DIAG5D',255)
C****
C**** KINETIC ENERGY
C****
C**** TRANSFER RATES FOR KINETIC ENERGY IN THE DYNAMICS
  100 CALL GETTIME(MBEGIN)
      KE(:,:)=0.
      KE_part(:,:,:)=0.

      DO L=1,LM
        DO J=J_0STG,J_1STG
          KSPHER=KLAYER(L)
          IF (J > JEQ) KSPHER= KSPHER+1
          DO KUV=1,2 ! loop over u,v
            IF(KUV.EQ.1) CALL FFT(DUT(1,J,L),FA,FB)
            IF(KUV.EQ.2) CALL FFT(DVT(1,J,L),FA,FB)
            DO N=1,NM
              X(N)=.5*FIM*
     &          (FA(N-1)*FCUVA(N-1,J,L,KUV)+FB(N-1)*FCUVB(N-1,J,L,KUV))
            ENDDO
            X(1)=X(1)+X(1)
            X(NM)=X(NM)+X(NM)
            IF (J.NE.JEQ) KE_part(:,J,KSPHER)=KE_part(:,J,KSPHER)+
     &                                        X(:)*DSIG(L)
            IF (J.EQ.J45N) THEN     ! 45 N
               KE_part(:,J,KSPHER+2)=KE_part(:,J,KSPHER+2)+X(:)*DSIG(L)
            ELSE IF (J.EQ.JEQ) THEN ! EQUATOR
              DO N=1,NM
                KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                X(N)*DSIG(L)
                KE_part(N,J,KSPHER  )=KE_part(N,J,KSPHER  )+
     &                                .5D0*X(N)*DSIG(L)       ! CONTRIB TO SH
                KE_part(N,J,KSPHER+1)=KE_part(N,J,KSPHER+1)+
     &                                .5D0*X(N)*DSIG(L)       ! CONTRIB TO NH
              ENDDO
              IF (KUV.EQ.2) KSPHER=KSPHER+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL GLOBALSUM(grid, KE_part(1:NM,:,1:NSPHER), KE(1:NM,1:NSPHER),
     &   ALL=.TRUE.)

      DO 180 KSPHER=1,NSPHER
      DO 180 N=1,NM
  180 SPECA(N,MKE,KSPHER)=SPECA(N,MKE,KSPHER)+KE(N,KSPHER)/NDT
C****
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAG5D


      SUBROUTINE DIAG5A (M5,NDT)
C****
C**** THIS DIAGNOSTICS ROUTINE PRODUCES A SPECTRAL ANALYSIS OF KINETIC
C**** AND AVAILABLE POTENTIAL ENERGIES AND THEIR TRANSFER RATES BY
C**** VARIOUS ATMOSPHERIC PROCESSES.
C****
C**** THE PARAMETER M INDICATES WHAT IS STORED IN SPECA(N,M,KSPHER),
C**** IT ALSO INDICATES WHEN DIAG5A IS BEING CALLED.
C**** M=1  MEAN STANDING KINETIC ENERGY            BEFORE SOURCES
C****   2  MEAN KINETIC ENERGY                     BEFORE DYNAMICS
C****   3  MEAN POTENTIAL ENERGY
C****   4  CONVERSION OF K.E. BY ADVECTION         AFTER ADVECTION
C****   5  CONVERSION OF K.E. BY CORIOLIS FORCE    AFTER CORIOLIS TERM
C****   6  CONVERSION FROM P.E. INTO K.E.          AFTER PRESS GRAD FORC
C****   7  CHANGE OF K.E. BY DYNAMICS              AFTER DYNAMICS
C****   8  CHANGE OF P.E. BY DYNAMICS
C****   9  CHANGE OF K.E. BY CONDENSATION          AFTER CONDENSATION
C****  10  CHANGE OF P.E. BY CONDENSATION
C****  11  CHANGE OF P.E. BY RADIATION             AFTER RADIATION
C****  12  CHANGE OF K.E. BY SURFACE               AFTER SURFACE
C****  13  CHANGE OF P.E. BY SURFACE
C****  14  CHANGE OF K.E. BY FILTER                AFTER FILTER
C****  15  CHANGE OF P.E. BY FILTER
C****  16  CHANGE OF K.E. BY DAILY                 AFTER DAILY
C****  17  CHANGE OF P.E. BY DAILY
C****  18  UNUSED
C****  19  LAST KINETIC ENERGY
C****  20  LAST POTENTIAL ENERGY
C****
      USE CONSTANT, only : sha
      USE MODEL_COM, only : im,imh,jm,lm,fim,
     &     DSIG,IDACC,JEQ,LS1,MDIAG,
     &     P,PTOP,PSFMPT,SIG,T,U,V,ZATMO
      USE GEOM, only : AREAG,DXYN,DXYP,DXYS
      USE DIAG_COM, only : speca,atpe,nspher,kspeca,klayer
      USE DIAG_COM, only : SQRTM
      USE DYNAMICS, only : sqrtp,pk
      USE DOMAIN_DECOMP, only : GRID,GET,CHECKSUM,HALO_UPDATE, AM_I_ROOT
      USE DOMAIN_DECOMP, only : GLOBALSUM, SOUTH, WRITE_PARALLEL

      IMPLICIT NONE
      INTEGER :: M5,NDT
      REAL*8, DIMENSION(IM) :: X
      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE,APE
      REAL*8, DIMENSION
     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part
      REAL*8, DIMENSION(IMH+1,4,LM) :: VAR
      REAL*8, DIMENSION(IMH+1,4,LM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
     &     :: VAR_part
      REAL*8 :: VARSUM
      REAL*8, DIMENSION(2) :: TPE
      REAL*8               :: TPE_sum
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: TPE_psum
CMoved to DAGCOM so it could be declared allocatable      REAL*8, SAVE, DIMENSION(IM,JM) :: SQRTM

      REAL*8, DIMENSION(LM) :: THJSP,THJNP,THGM

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KSPECA), PARAMETER ::
     &     MTPEOF=(/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/)

      INTEGER :: I,IJL2,IP1,J,J45N,JH,JHEMI,JP,K,KS,KSPHER,L,LDN,
     &     LUP,MAPE,MKE,MNOW,MTPE,N,NM

      REAL*8 :: GMEAN(LM),GMTMP,SQRTPG,SUMI,SUMT,THGSUM,THJSUM
      REAL*8 :: THGSUM_part(grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 :: GMSUM(grid%J_STRT_HALO:grid%J_STOP_HALO,LM)

      INTEGER :: J_0S,J_1S,J_0STG,J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(GRID, J_STRT_SKP=J_0S   , J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      SQRTPG = SQRT(PSFMPT)
      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      IJL2=IM*JM*LM*2

      MKE=M5
      MAPE=M5
C****
C**** Note: KSPHER has been re-arranged from previous models to better
C****       deal with optionally resolved stratosphere. The higher
C****       values are only used if the model top is high enough.
C****
C**** KSPHER=1 SOUTHERN TROPOSPHERE         2 NORTHERN TROPOSPHERE
C****        3 EQUATORIAL TROPOSPHERE       4 45 DEG NORTH TROPOSPHERE
C****
C****        5 SOUTHERN LOW STRATOSPHERE    6 NORTHERN LOW STRATOSPHERE
C****        7 EQUATORIAL LOW STRATOSPHERE  8 45 DEG NORTH LOW STRATOSPH
C****
C****        9 SOUTHERN MID STRATOSPHERE   10 NORTHERN MID STRATOSPHERE
C****       11 EQUATORIAL MID STRATOSPHERE 12 45 DEG NORTH MID STRATOSPH
C****
C****       13 SOUTHERN UPP STRATOSPHERE   14 NORTHERN UPP STRATOSPHERE
C****       15 EQUATORIAL UPP STRATOSPHERE 16 45 DEG NORTH UPP STRATOSPH
C****
      GO TO (200,200,810,810,810,  810,200,810,205,810,
     *       296,205,810,205,810,  205,810,810,810,810),M5

C***  810 WRITE (6,910) M5
  810 CALL WRITE_PARALLEL(M5, UNIT=6, format=
     & "('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.  M5=',I5)")
C****  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.  M5=',I5)
      call stop_model('INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.',255)
C**** MASS FOR KINETIC ENERGY
  200 CONTINUE

      I=IM
      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      DO J=J_0STG,J_1STG
        DO  IP1=1,IM
          SQRTM(I,J)=SQRT(.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+
     *         P(IP1,J-1))*DXYN(J-1)))
          I=IP1
        END DO
      END DO

C****
  205 CONTINUE

      MAPE=MKE+1
      KE(:,:)=0.
      KE_part(:,:,:)=0.
C**** CURRENT KINETIC ENERGY
      DO L=1,LM
        DO J=J_0STG,J_1STG
          IF (J <= JEQ) THEN
            KSPHER=KLAYER(L)
          ELSE
            KSPHER=KLAYER(L)+1
          END IF
          DO K = IZERO,LM,LM
            IF(K.EQ.IZERO)X(1:IM)=U(1:IM,J,L)*SQRTM(1:IM,J)
            IF(K.EQ.LM)   X(1:IM)=V(1:IM,J,L)*SQRTM(1:IM,J)
            CALL FFTE (X,X)
            IF (J.EQ.JEQ) THEN
              DO N=1,NM
                KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+X(N)*DSIG(L)
                KE_part(N,J,KSPHER  )=KE_part(N,J,KSPHER  )+
     &                                .5D0*X(N)*DSIG(L)
                KE_part(N,J,KSPHER+1)=KE_part(N,J,KSPHER+1)+
     &                                .5D0*X(N)*DSIG(L)
              ENDDO
cgsfc              IF(K.EQ.LM)KSPHER=KSPHER+1
            ELSE
              DO N=1,NM
                KE_part(N,J,KSPHER)=KE_part(N,J,KSPHER)+X(N)*DSIG(L)
              ENDDO
              IF (J.EQ.J45N) THEN
                DO N=1,NM
                  KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                X(N)*DSIG(L)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL GLOBALSUM(grid, KE_part, KE, ALL=.true.)

      IF (NDT /= 0) THEN
C**** TRANSFER RATES AS DIFFERENCES OF KINETIC ENERGY
        DO KS=1,NSPHER
          DO N=1,NM
           SPECA(N,MKE,KS)=SPECA(N,MKE,KS)+(KE(N,KS)-SPECA(N,19,KS))/NDT
          END DO
        END DO
      END IF
      DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,19,KS)=KE(N,KS)
        END DO
      END DO

C****
C**** POTENTIAL ENERGY
C****
  296 CONTINUE

      APE(:,:)=0.
C**** CURRENT AVAILABLE POTENTIAL ENERGY

      DO L = 0, LM
        LDN = MAX(L - 1, 1)
        LUP = MIN(L + 1,LM)

        IF (L < LM) THEN
          IF(LUP.GE.LS1) THEN
            IF(HAVE_SOUTH_POLE) THJSP(LUP) = T(1,1,LUP)*SQRTPG
            IF(HAVE_NORTH_POLE) THJNP(LUP) = T(1,JM,LUP)*SQRTPG
          ELSE
            IF(HAVE_SOUTH_POLE) THJSP(LUP)=T(1,1,LUP)*SQRTP(1,1)
            IF(HAVE_NORTH_POLE) THJNP(LUP)=T(1,JM,LUP)*SQRTP(1,JM)
          ENDIF
          IF(HAVE_SOUTH_POLE) THGSUM_part(1)  = FIM*THJSP(LUP)*DXYP(1)
          IF(HAVE_NORTH_POLE) THGSUM_part(JM) = FIM*THJNP(LUP)*DXYP(JM)

          DO J=J_0S,J_1S
            THJSUM=0.
            DO I=1,IM
              THJSUM=THJSUM+T(I,J,LUP)*SQRTP(I,J)
            ENDDO
            THGSUM_part(J) = THJSUM*DXYP(J)
          ENDDO

          CALL GLOBALSUM(grid, THGSUM_part, THGSUM, ALL=.TRUE.)
          THGM(LUP)=THGSUM/AREAG

        End IF

        IF (LUP < 2) CYCLE

        VAR(2:NM,1:2,L)=0.
        VAR_part(:,:,L,:)=0
        IF(HAVE_SOUTH_POLE) THEN
          VAR_part(1,1,L,1)=.5*(THJSP(L)-THGM(L))**2*DXYP(1)*FIM
         GMSUM(1,L)=((THJSP(LUP)-THJSP(LDN))*DXYP(1)*
     *         (SIG(L)*P(1,1)+PTOP)/
     *         (SQRTP(1,1)*P(1,1)*PK(L,1,1)))*FIM
        END IF
        IF(HAVE_NORTH_POLE) THEN
          VAR_part(1,2,L,JM)=.5*(THJNP(L)-THGM(L))**2*DXYP(JM)*FIM
          GMSUM(JM,L)=((THJNP(LUP)-THJNP(LDN))*DXYP(JM)*
     *       (SIG(L)*P(1,JM)+PTOP)/(SQRTP(1,JM)*P(1,JM)*PK(L,1,JM)))*FIM
        END IF

        DO J=J_0S,J_1S

          IF (J < JEQ) THEN
            JHEMI = 1
          ELSE
            JHEMI = 2
          END IF

          GMTMP = 0
          DO I=1,IM
            X(I)=T(I,J,L)*SQRTP(I,J)-THGM(L)
            GMTMP=GMTMP+(T(I,J,LUP)-T(I,J,LDN))*(SIG(L)*P(I,J)+PTOP)/
     *           (P(I,J)*PK(L,I,J))
          END DO
          GMSUM(J,L) = GMTMP * DXYP(J)
          CALL FFTE (X,X)
          DO N=1,NM
            VAR_part(N,JHEMI,L,J)=VAR_part(N,JHEMI,L,J)+X(N)*DXYP(J)
          END DO
          IF (J == JEQ-1) THEN
            DO N=1,NM
              VAR_part(N,3,L,J)=X(N)*DXYP(J)
            END DO
          END IF
          IF (J == J45N-1) THEN
            DO N=1,NM
              VAR_part(N,4,L,J)=X(N)*DXYP(J)
            END DO
          END IF
        END DO

      END DO

      CALL GLOBALSUM(grid, GMSUM, GMEAN)
      CALL GLOBALSUM(grid, VAR_part, VAR)

      IF (AM_I_ROOT()) THEN
        DO L = 1, LM
          LDN = MAX(L - 1, 1)
          LUP = MIN(L + 1,LM)
          GMEAN(L)=DSIG(L)*AREAG*(SIG(LDN)-SIG(LUP))/GMEAN(L)
          KS=KLAYER(L)
          DO JHEMI=1,4
            DO N=1,NM
              APE(N,KS)=APE(N,KS)+VAR(N,JHEMI,L)*GMEAN(L)
            END DO
            KS=KS+1
          END DO
        END DO
      END IF
C**** CURRENT TOTAL POTENTIAL ENERGY
 450  CONTINUE

      IF (HAVE_SOUTH_POLE) THEN
        J=1
        SUMT=0
        DO L=1, LM
          SUMT=SUMT + T(1,J,L)*PK(L,1,J)*DSIG(L)
        END DO
        TPE_psum(J)=FIM*DXYP(J)*(ZATMO(1,J)*(P(1,J)+PTOP)+
     *       SUMT*SHA*P(1,J))
      END IF
      IF (HAVE_NORTH_POLE) THEN
        J=JM
        SUMT=0
        DO L=1, LM
          SUMT=SUMT + T(1,J,L)*PK(L,1,J)*DSIG(L)
        END DO
        TPE_psum(J)=FIM*DXYP(J)*(ZATMO(1,J)*(P(1,J)+PTOP)+
     *       SUMT*SHA*P(1,J))
      END IF
      DO J=J_0S, J_1S
        SUMI=0
        DO I=1,IM
          SUMT=0
          DO L=1,LM
            SUMT=SUMT + T(I,J,L)*PK(L,I,J)*DSIG(L)
          END DO
          SUMI=SUMI+ZATMO(I,J)*(P(I,J)+PTOP)+SUMT*SHA*P(I,J)
        END DO
        TPE_psum(J) = SUMI*DXYP(J)
      END DO
      CALL GLOBALSUM(grid, TPE_psum, TPE_sum, TPE)

      IF (NDT /= 0) THEN
        MTPE=MTPEOF(MAPE)

C**** TRANSFER RATES AS DIFFERENCES FOR POTENTIAL ENERGY
        DO KS=1,NSPHER
          DO N=1,NM
        SPECA(N,MAPE,KS)=SPECA(N,MAPE,KS)+(APE(N,KS)-SPECA(N,20,KS))/NDT
          END DO
        END DO

        ATPE(MTPE,1)=ATPE(MTPE,1)+(TPE(1)-ATPE(8,1))/NDT
        ATPE(MTPE,2)=ATPE(MTPE,2)+(TPE(2)-ATPE(8,2))/NDT
      END IF
      DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,20,KS)=APE(N,KS)
        END DO
      END DO
      ATPE(8,1)=TPE(1)
      ATPE(8,2)=TPE(2)

      IF (M5.EQ.2) THEN
C**** ACCUMULATE MEAN KINETIC ENERGY AND MEAN POTENTIAL ENERGY
        DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,2,KS)=SPECA(N,2,KS)+KE(N,KS)
          SPECA(N,3,KS)=SPECA(N,3,KS)+APE(N,KS)
        END DO
        END DO
        ATPE(1,1)=ATPE(1,1)+TPE(1)
        ATPE(1,2)=ATPE(1,2)+TPE(2)
      END IF

      CALL TIMER (MNOW,MDIAG)
      RETURN
      END SUBROUTINE DIAG5A


      SUBROUTINE DIAG5F(UX,VX)
C**** FOURIER COEFFICIENTS FOR CURRENT WIND FIELD
C****
      USE MODEL_COM, only : im,imh,jm,lm,
     &     IDACC,MDIAG,MDYN
      USE DIAG_LOC, only : FCUVA,FCUVB
      USE DOMAIN_DECOMP, only : GRID,GET
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        UX,VX
c      REAL*8, DIMENSION(0:IMH,JM,LM,2) :: FCUVA,FCUVB
c      COMMON/WORK7/FCUVA,FCUVB
      INTEGER :: J,L,MBEGIN
      INTEGER :: J_0STG, J_1STG

      CALL GET(GRID, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)
      CALL GETTIME(MBEGIN)
      IDACC(6)=IDACC(6)+1
      DO L=1,LM
         DO J=J_0STG,J_1STG
            CALL FFT(UX(1,J,L),FCUVA(0,J,L,1),FCUVB(0,J,L,1))
            CALL FFT(VX(1,J,L),FCUVA(0,J,L,2),FCUVB(0,J,L,2))
         ENDDO
      ENDDO
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)

      RETURN
      END SUBROUTINE DIAG5F


      SUBROUTINE DIAG4A
C****
C**** THIS ROUTINE PRODUCES A TIME HISTORY OF ENERGIES
C****
      USE MODEL_COM, only : im,jm,lm,
     &     IDACC,JEQ,LS1,ISTRAT          !! ,SKIPSE
      USE GEOM, only : DXYV
      USE DIAG_COM, only : energy,speca,ajk=>ajk_loc
     &     ,aijk=>aijk_loc,ijk_u,ijk_v,ijk_dp,ned
      IMPLICIT NONE

      INTEGER ::
     &     I,IDACC5,J,KS,L,N,NM
      REAL*8 ::
     &     BYIADA,PU4TI,PV4TI,SEKE,SKE4I

      IF (IDACC(4).LE.0.OR.IDACC(7).LE.0) RETURN
C      JEQ=2.+.5*(JM-1.)
      NM=1+IM/2
C****
C**** LOAD ENERGIES INTO TIME HISTORY ARRAY
C****
      IDACC5=IDACC(5)+1
      IF (IDACC5.GT.100) RETURN
      DO I=0,1+ISTRAT  ! loop over number of 'spheres'
        ENERGY(1+NED*I,IDACC5)=SPECA(1,19,1+4*I)   ! SH
        ENERGY(2+NED*I,IDACC5)=SPECA(1,19,2+4*I)   ! NH
        ENERGY(5+NED*I,IDACC5)=SPECA(2,19,2+4*I)   ! NH wave 1
        ENERGY(6+NED*I,IDACC5)=SPECA(3,19,2+4*I)   ! NH wave 2
        ENERGY(7+NED*I,IDACC5)=SPECA(1,20,1+4*I)
        ENERGY(8+NED*I,IDACC5)=SPECA(1,20,2+4*I)
        DO N=2,NM
        ENERGY( 3+NED*I,IDACC5)=ENERGY( 3+10*I,IDACC5)+SPECA(N,19,1+4*I)
        ENERGY( 4+NED*I,IDACC5)=ENERGY( 4+10*I,IDACC5)+SPECA(N,19,2+4*I)
        ENERGY( 9+NED*I,IDACC5)=ENERGY( 9+10*I,IDACC5)+SPECA(N,20,1+4*I)
        ENERGY(10+NED*I,IDACC5)=ENERGY(10+10*I,IDACC5)+SPECA(N,20,2+4*I)
        END DO
      END DO
      IDACC(5)=IDACC5
      RETURN
C****
      END SUBROUTINE DIAG4A

      module subdaily
!@sum SUBDAILY defines variables associated with the sub-daily diags
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,lm,itime
      USE FILEMANAGER, only : openunit, closeunits, nameunit
      USE DIAG_COM, only : kgz_max,pmname,P_acc
      USE PARAM
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE TRACER_COM, only : trm
#endif
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch)
     &       ,tr_mm
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRACER_COM, only : n_Ox
#endif
#ifdef TRACERS_AEROSOLS_Koch
      USE TRACER_COM, only : n_SO4
#ifdef TRACERS_HETCHEM
     *       ,n_SO4_d1,n_SO4_d2, n_SO4_d3
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE TRACER_COM,ONLY : Ntm,Ntm_dust,trname
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#ifdef TRACERS_WATER
     &     ,dowetdep
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_AEROSOLS_Koch)
     &     ,n_clay
#else
#ifdef TRACERS_MINERALS
     &     ,n_clayilli
#else
#ifdef TRACERS_QUARZHEM
     &     ,n_sil1quhe
#endif
#endif
#endif
#ifdef TRACERS_SPECIAL_O18
      USE TRACER_COM, only : n_water, n_HDO, trm, trw0
#endif
      IMPLICIT NONE
      SAVE
!@var kddmax maximum number of sub-daily diags output files
      INTEGER, PARAMETER :: kddmax = 42
!@var kdd total number of sub-daily diags
      INTEGER :: kdd
!@var kddunit total number of sub-daily files
      INTEGER :: kddunit
!@var namedd array of names of sub-daily diags
      CHARACTER*10, DIMENSION(kddmax) :: namedd
!@var iu_subdd array of unit numbers for sub-daily diags output
      INTEGER, DIMENSION(kddmax) :: iu_subdd
!@var subddt = subdd + subdd1 = all variables for sub-daily diags
      CHARACTER*128 :: subddt = " "
!@dbparam subdd string contains variables to save for sub-daily diags
!@dbparam subdd1 additional string of variables for sub-daily diags
C**** Note: for longer string increase MAX_CHAR_LENGTH in PARAM
      CHARACTER*64  :: subdd  = "SLP"
      CHARACTER*64  :: subdd1 = " "
!@dbparam Nsubdd: DT_save_SUBDD =  Nsubdd*DTsrc sub-daily diag freq.
      INTEGER :: Nsubdd = 0
!@dbparam LmaxSUBDD: the max L when writing "ALL" levels
      INTEGER :: LmaxSUBDD = LM

      contains

      subroutine init_subdd(aDATE)
!@sum init_subdd initialise sub daily diags and position files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate
      character*12 name
      integer :: i,j,k,kunit,kk

      call sync_param( "subdd" ,subdd)
      call sync_param( "subdd1" ,subdd1)
      call sync_param( "Nsubdd",Nsubdd)
      call sync_param( "LmaxSUBDD",LmaxSUBDD)

      if (nsubdd.ne.0) then
C**** combine strings subdd1 and subdd2:
        subddt=trim(subdd)//' '//subdd1
C**** calculate how many names
        k=0
        i=1
 10     j=index(subddt(i:len(subddt))," ")
        if (j.gt.1) then
          k=k+1
          i=i+j
        else
          i=i+1
        end if
        if (i.lt.len(subddt)) goto 10
        kdd=k
        if (kdd.gt.kddmax) call stop_model
     *       ("Increase kddmax: No. of sub-daily diags too big",255)

C**** make array of names
        read(subddt,*) namedd(1:kdd)

C**** open units and position
C**** Some names have more than one unit associated (i.e. "ZALL")
        kunit=0
        do k=1,kdd
          if (namedd(k)(len_trim(namedd(k))-2:len_trim(namedd(k))).eq.
     *         "ALL") then
            select case (namedd(k)(1:1))
            case ("U", "V", "W", "C")! velocities/cloud opt depth on layers
              kunit=kunit+1
              write(name,'(A1,A3,A7)') namedd(k)(1:1),'ALL',aDATE(1:7)
              call openunit(name,iu_SUBDD(kunit),.true.,.false.)
              call io_POS(iu_SUBDD(kunit),Itime,im*jm,Nsubdd)
            case ("Z", "T", "R", "Q")! heights, temps, rel/spec hum,PMB levels
              do kk=1,kgz_max
                kunit=kunit+1
                call openunit(namedd(k)(1:1)//trim(PMNAME(kk))//
     *               aDATE(1:7),iu_SUBDD(kunit),.true.,.false.)
                call io_POS(iu_SUBDD(kunit),Itime,im*jm,Nsubdd)
              end do
#ifdef TRACERS_SPECIAL_Shindell
            case ("O")  ! Ox tracer
              kunit=kunit+1
              write(name,'(A2,A3,A7)') namedd(k)(1:2),'ALL',aDATE(1:7)
              call openunit(name,iu_SUBDD(kunit),.true.,.false.)
              call io_POS(iu_SUBDD(kunit),Itime,im*jm,Nsubdd)
#endif
#ifdef TRACERS_SPECIAL_O18
            case ("D")  ! HDO tracer
              kunit=kunit+1
              write(name,'(A1,A3,A7)') namedd(k)(1:1),'ALL',aDATE(1:7)
              call openunit(name,iu_SUBDD(kunit),.true.,.false.)
              call io_POS(iu_SUBDD(kunit),Itime,im*jm,Nsubdd)
#endif
            end select
          else                  ! single file per name
            kunit=kunit+1
            call openunit(trim(namedd(k))//aDATE(1:7),iu_SUBDD(kunit),
     *           .true.,.false.)
            call io_POS(iu_SUBDD(kunit),Itime,im*jm,Nsubdd)
          endif
        end do
        kddunit=kunit
      end if
C**** initialise special subdd accumulation
      P_acc=0.

      return
      end subroutine init_subdd

      subroutine reset_subdd(aDATE)
!@sum reset_subdd resets sub daily diag files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate
      character*12 name
      integer :: k,kunit,kk

      if (nsubdd.ne.0) then
C**** close and re-open units
        call closeunits ( iu_SUBDD, kddunit )
        kunit=0
        do k=1,kdd
          if (namedd(k)(len_trim(namedd(k))-2:len_trim(namedd(k))).eq.
     *         "ALL") then
            select case (namedd(k)(1:1))
            case ("U", "V", "W", "C")! velocities on model layers
              kunit=kunit+1
              write(name,'(A1,A3,A7)') namedd(k)(1:1),'ALL',aDATE(1:7)
              call openunit(name,iu_SUBDD(kunit),.true.,.false.)
            case ("Z", "T", "R", "Q")! heights, temps, rel/spec hum PMB levels
              do kk=1,kgz_max
                kunit=kunit+1
                call openunit(namedd(k)(1:1)//trim(PMNAME(kk))//
     *               aDATE(1:7),iu_SUBDD(kunit),.true.,.false.)
              end do
#ifdef TRACERS_SPECIAL_Shindell
            case ("O")  ! Ox tracer
              kunit=kunit+1
              write(name,'(A2,A3,A7)') namedd(k)(1:2),'ALL',aDATE(1:7)
              call openunit(name,iu_SUBDD(kunit),.true.,.false.)
#endif
#ifdef TRACERS_SPECIAL_O18
            case ("D")  ! HDO tracer
              kunit=kunit+1
              write(name,'(A1,A3,A7)') namedd(k)(1:1),'ALL',aDATE(1:7)
              call openunit(name,iu_SUBDD(kunit),.true.,.false.)
#endif
            end select
          else                  ! single file per name
            kunit=kunit+1
            call openunit(trim(namedd(k))//aDATE(1:7),iu_SUBDD(kunit),
     *           .true.,.false.)
          endif
        end do
      end if
C****
      return
      end subroutine reset_subdd

      subroutine get_subdd
!@sum get_SUBDD saves instantaneous variables at sub-daily frequency
!@+   every ABS(NSUBDD)
!@+   Note that TMIN, TMAX can only be output once a day
!@+   Current options: SLP, PS, SAT, PREC, QS, LCLD, MCLD, HCLD, PTRO
!@+                    QLAT, QSEN, SWD, SWU, LWD, LWU, LWT, STX, STY,
!@+                    ICEF, SNOWD, TCLD, SST, SIT, US, VS, TMIN, TMAX
!@+                    Z*, R*, T*, Q*  (on any fixed pressure level)
!@+                    U*, V*, W*, C*  (on any model level)
!@+                    Ox*         (on any model level with chemistry)
!@+                    D*          (HDO on any model level)
!@+                    SO4
#ifdef CLD_AER_CDNC
!@+                    CTEM,CD3D,CI3D,CL3D,CDN3D,CRE3D,CLWP
!@+                    TAUSS,TAUMC,CLDSS,CLDMC
#endif
#ifdef TRACERS_HETCHEM
!@+                    SO4_d1,SO4_d2,SO4_d3,
!@+                    Clay, Silt1, Silt2, Silt3
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@+                    DUEMIS,DUDEPTURB,DUDEPGRAV,DUDEPWET,DUTRS,DULOAD
#endif
#ifdef TRACERS_DUST
!@+                    DUEMIS2
#endif
!@+   More options can be added as extra cases in this routine
!@auth Gavin Schmidt/Reto Ruedy
      USE CONSTANT, only : grav,rgas,bygrav,bbyg,gbyrb,sday,tf,mair,sha
     *     ,lhe,rhow,undef,stbo
      USE MODEL_COM, only : lm,p,ptop,zatmo,dtsrc,u,v,focean,fearth
     *     ,flice,nday
      USE GEOM, only : imaxj,dxyp
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,bydxyp
#endif
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg
      USE CLOUDS_COM, only : llow,lmid,lhi,cldss,cldmc,taumc,tauss,fss
#ifdef CLD_AER_CDNC
     *           ,ctem,cd3d,ci3d,cl3d,cdn3d,cre3d,clwp
#endif
      USE DYNAMICS, only : ptropo,am,wsave
      USE FLUXES, only : prec,dmua,dmva,tflux1,qflux1,uflux1,vflux1
     *     ,gtemp
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,dust_flux_glob,trs_glob
#ifdef TRACERS_DRYDEP
     &     ,depo_turb_glob,depo_grav_glob
#endif
#ifdef TRACERS_WATER
     &     ,trprec
#else
     &     ,trprec_dust
#endif
#endif
#ifdef TRACERS_DUST
     &     ,dust_flux2_glob
#endif
#ifdef TRACERS_ON
      USE TRACER_COM
#endif
      USE SEAICE_COM, only : rsi,snowi
      USE LANDICE_COM, only : snowli
      USE LAKES_COM, only : flake
      USE GHY_COM, only : snowe
      USE RAD_COM, only : trhr,srdn,salb,cfrac,cosz1
      USE DIAG_COM, only : z_inst,rh_inst,t_inst,kgz_max,pmname,tdiurn
     *     ,p_acc,pmb
      USE DOMAIN_DECOMP, only : GRID,GET,WRITEI_PARALLEL
      IMPLICIT NONE
      REAL*4, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: DATA
      INTEGER :: I,J,K,L,kp,kunit
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,n,n1,n_fidx
#endif
      CHARACTER namel*3
      REAL*8 POICE,PEARTH,PLANDI,POCEAN,QSAT
      INTEGER :: J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      kunit=0
C**** depending on namedd string choose what variables to output
      do k=1,kdd

C**** simple diags
        select case (namedd(k))
        case ("SLP")            ! sea level pressure (mb)
          do j=J_0,J_1
          do i=1,imaxj(j)
            data(i,j)=(p(i,j)+ptop)*(1.+bbyg*zatmo(i,j)/tsavg(i,j))
     *           **gbyrb
          end do
          end do
        case ("PS")             ! surface pressure (mb)
          do j=J_0,J_1
          do i=1,imaxj(j)
            data(i,j)=p(i,j)+ptop
          end do
          end do
        case ("SAT")            ! surf. air temp (C)
          data=tsavg-tf
        case ("US")             ! surf. u wind (m/s)
          data=usavg
        case ("VS")             ! surf. v wind (m/s)
          data=vsavg
        case ("SST")            ! sea surface temp (C)
          do j=J_0,J_1
            do i=1,imaxj(j)
              if (FOCEAN(I,J)+FLAKE(I,J).gt.0) then
                data(i,j)=GTEMP(1,1,i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("SIT")            ! surface ice temp (C)
          do j=J_0,J_1
            do i=1,imaxj(j)
              if (RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J)).gt.0) then
                data(i,j)=GTEMP(1,2,i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("QS")             ! surf humidity (kg/kg)
          data=qsavg
        case ("PREC")           ! precip (mm/day)
c          data=sday*prec/dtsrc
          data=sday*P_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          P_acc=0.
        case ("SNOWD")     ! snow depth (w.e. mm)
          do j=J_0,J_1
            do i=1,imaxj(j)
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              data(i,j)=1d3*(SNOWI(I,J)*POICE+SNOWLI(I,J)*PLANDI+SNOWE(I
     *             ,J)*PEARTH)/RHOW
            end do
          end do
        case ("QLAT")           ! latent heat (W/m^2)
          data=qflux1*lhe
        case ("QSEN")           ! sensible heat flux (W/m^2)
          data=tflux1*sha
        case ("SWD")            ! solar downward flux at surface (W/m^2)
          data=srdn*cosz1       ! multiply by instant cos zenith angle
        case ("SWU")            ! solar upward flux at surface (W/m^2)
! estimating this from the downward x albedo, since that's already saved
          data=srdn*(1.-salb)*cosz1
        case ("LWD")            ! LW downward flux at surface (W/m^2)
          data=TRHR(0,:,:)
        case ("LWU")            ! LW upward flux at surface (W/m^2)
          do j=J_0,J_1
            do i=1,imaxj(j)
              POCEAN=(1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              data(i,j)=STBO*(POCEAN*(GTEMP(1,1,I,J)+TF)**4+
     *             POICE *(GTEMP(1,2,I,J)+TF)**4+
     *             PLANDI*(GTEMP(1,3,I,J)+TF)**4+
     *             PEARTH*(GTEMP(1,4,I,J)+TF)**4)
            end do
          end do
        case ("LWT")            ! LW upward flux at TOA (P1) (W/m^2)
          do j=J_0,J_1     ! sum up all cooling rates + net surface emission
            do i=1,imaxj(j)
              POCEAN=(1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              data(i,j)=-SUM(TRHR(0:LM,I,J))+
     *             STBO*(POCEAN*(GTEMP(1,1,I,J)+TF)**4+
     *             POICE *(GTEMP(1,2,I,J)+TF)**4+
     *             PLANDI*(GTEMP(1,3,I,J)+TF)**4+
     *             PEARTH*(GTEMP(1,4,I,J)+TF)**4)
            end do
          end do
        case ("ICEF")           ! ice fraction over open water (%)
          data=RSI*100.
        case ("STX")            ! E-W surface stress (N/m^2)
          data=uflux1
        case ("STY")            ! N-S surface stress (N/m^2)
          data=vflux1
        case ("LCLD")           ! low level cloud cover (%)
          data=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=1,imaxj(j)
              do l=1,llow
                data(i,j)=data(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              data(i,j)=data(i,j)*100./real(llow,kind=8)
            end do
          end do
        case ("MCLD")           ! mid level cloud cover (%)
          data=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=1,imaxj(j)
              do l=llow+1,lmid
                data(i,j)=data(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              data(i,j)=data(i,j)*100./real(lmid-llow,kind=8)
            end do
          end do
        case ("HCLD")           ! high level cloud cover (%)
          data=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=1,imaxj(j)
              do l=lmid+1,lhi
                data(i,j)=data(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              data(i,j)=data(i,j)*100./real(lhi-lmid,kind=8)
            end do
          end do
        case ("TCLD")           ! total cloud cover (%) (As seen by rad)
          data=cfrac*100.
        case ("PTRO")           ! tropopause pressure (mb)
          data = ptropo
        case ("TMIN")           ! min daily temp (C)
          if (mod(itime+1,Nday).ne.0) then ! only at end of day
            kunit=kunit+1
            cycle
          end if
          data=tdiurn(:,:,9)
        case ("TMAX")           ! max daily temp (C)
          if (mod(itime+1,Nday).ne.0) then ! only at end of day
            kunit=kunit+1
            cycle
          end if
          data=tdiurn(:,:,6)
#ifdef TRACERS_AEROSOLS_Koch
        case ("SO4")      ! sulfate in L=1
          do j=J_0,J_1
          do i=1,imaxj(j)
            data(i,j)=trm(i,j,1,n_SO4)
#ifdef TRACERS_HETCHEM
     *               +trm(i,j,1,n_SO4_d1)
     *               +trm(i,j,1,n_SO4_d2)
     *               +trm(i,j,1,n_SO4_d3)
#endif

          end do
          end do
#endif
        case default
          goto 10
        end select
        kunit=kunit+1
C**** fix polar values
        IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
        IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
C**** write out
        !call writei(iu_subdd(kunit),itime,data,im*jm)
        call writei_parallel(grid,
     &       iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
        cycle

C**** diags on fixed pressure levels or velocity
 10     select case (namedd(k)(1:1))
        case ("Z","R","T","Q")  ! heights, rel/spec humidity or temp
C**** get pressure level
          do kp=1,kgz_max
            if (namedd(k)(2:5) .eq. PMNAME(kp)) then
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("Z")        ! geopotential heights
                data=z_inst(kp,:,:)
              case ("R")        ! relative humidity (wrt water)
                data=rh_inst(kp,:,:)
              case ("Q")        ! specific humidity
                do j=J_0,J_1
                do i=1,imaxj(j)
                  data(i,j)=rh_inst(kp,i,j)*qsat(t_inst(kp,i,j),lhe
     *                 ,PMB(kp))
                end do
                end do
              case ("T")        ! temperature (C)
                data=t_inst(kp,:,:)
              end select
C**** fix polar values
              IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
              IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
              cycle
            end if
          end do
          if (namedd(k)(2:4) .eq. "ALL") then
            do kp=1,kgz_max
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("Z")        ! geopotential heights
                data=z_inst(kp,:,:)
              case ("R")        ! relative humidity (wrt water)
                data=rh_inst(kp,:,:)
              case ("Q")        ! specific humidity
                do j=J_0,J_1
                do i=1,imaxj(j)
                  data(i,j)=rh_inst(kp,i,j)*qsat(t_inst(kp,i,j),lhe
     *                 ,PMB(kp))
                end do
                end do
              case ("T")        ! temperature (C)
                data=t_inst(kp,:,:)
              end select
C**** fix polar values
              IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
              IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
            end do
            cycle
          end if
        case ("U","V","W","C")    ! velocity/clouds on levels
          if (namedd(k)(2:4) .eq. "ALL") then
            kunit=kunit+1
            do kp=1,LmaxSUBDD
              select case (namedd(k)(1:1))
              case ("U")        ! E-W velocity
                data=u(:,:,kp)
              case ("V")        ! N-S velocity
                data=v(:,:,kp)
              case ("W")        ! vertical velocity
                data=wsave(:,:,kp)
C**** fix polar values for W only (calculated on tracer points)
                IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
                IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
              case ("C")        ! estimate of cloud optical depth
                do j=J_0,J_1
                  do i=1,imaxj(j)
                    data(i,j)=(1.-fss(kp,i,j))*taumc(kp,i,j)+fss(kp,i,j)
     *                   *tauss(kp,i,j)
                  end do
                end do
C**** fix polar values
                IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
                IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
              end select
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
            end do
            cycle
          end if
C**** get model level
          do l=1,lm
            if (l.lt.10) then
              write(namel,'(I1)') l
            else
              write(namel,'(I2)') l
            end if
            if (trim(namedd(k)(2:5)) .eq. trim(namel)) then
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("U")        ! U velocity
                data=u(:,:,l)
              case ("V")        ! V velocity
                data=v(:,:,l)
              case ("W")        ! W velocity
                data=wsave(:,:,l)
C**** fix polar values for W only (calculated on tracer points)
                data(2:im,1) =data(1,1)
                data(2:im,jm)=data(1,jm)
              end select
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
              cycle
            end if
          end do
#ifdef TRACERS_SPECIAL_Shindell
        case ("O")        ! Ox ozone tracer (ppmv)
          if (namedd(k)(3:5) .eq. "ALL") then
            kunit=kunit+1
            do kp=1,LmaxSUBDD
              do j=J_0,J_1
                do i=1,imaxj(j)
                  data(i,j)=1.d6*trm(i,j,kp,n_Ox)*mair/
     *                 (tr_mm(n_Ox)*am(kp,i,j)*dxyp(j))
                end do
              end do
C**** fix polar values
              IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
              IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
            end do
            cycle
          end if
C**** get model level
          do l=1,lm
            if (l.lt.10) then
              write(namel,'(I1)') l
            else
              write(namel,'(I2)') l
            end if
            if (trim(namedd(k)(3:6)) .eq. trim(namel)) then
              kunit=kunit+1
              do j=J_0,J_1
                do i=1,imaxj(j)
                  data(i,j)=1d6*trm(i,j,l,n_Ox)*mair/
     *                 (tr_mm(n_Ox)*am(l,i,j)*dxyp(j))
                end do
              end do
C**** fix polar values
              IF(HAVE_SOUTH_POLE) data(2:im,1) =data(1,1)
              IF(HAVE_NORTH_POLE) data(2:im,jm)=data(1,jm)
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
              cycle
            end if
          end do
#endif
#ifdef TRACERS_SPECIAL_O18
        case ("D")        ! HDO tracer (permil)
          if (namedd(k)(2:4) .eq. "ALL") then
            kunit=kunit+1
            do kp=1,LmaxSUBDD
              do j=J_0,J_1
                do i=1,imaxj(j)
                  data(i,j)=1d3*(trm(i,j,kp,n_HDO)/(trm(i,j,kp,n_water)
     *                 *trw0(n_HDO))-1.)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
            end do
            cycle
          end if
C**** get model level
          do l=1,lm
            if (l.lt.10) then
              write(namel,'(I1)') l
            else
              write(namel,'(I2)') l
            end if
            if (trim(namedd(k)(2:5)) .eq. trim(namel)) then
              kunit=kunit+1
              do j=J_0,J_1
                do i=1,imaxj(j)
                  data(i,j)=1d3*(trm(i,j,kp,n_HDO)/(trm(i,j,kp,n_water)
     *                 *trw0(n_HDO))-1.)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
C**** write out
              !call writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
              cycle
            end if
          end do
#endif
        end select
#ifdef TRACERS_HETCHEM
            select case (namedd(k))
        case ("SO2")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_SO2)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
        case ("SO4")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_SO4)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
        case ("SO4_d1")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_SO4_d1)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do        
        case ("SO4_d2")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)= trm(i,j,l,n_SO4_d2)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
        case ("SO4_d3")  
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_SO4_d3)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do        
        case ("Clay")  
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_Clay)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do        
        case ("Silt1")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_Silt1)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do  
        case ("Silt2")  
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_Silt2)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
        case ("Silt3")  
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=trm(i,j,l,n_Silt3)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
             end select
#endif
#ifdef CLD_AER_CDNC   !for 3 hrly diagnostics
            select case (namedd(k))
              case ("CTEM")
              kunit=kunit+1
              do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=ctem(l,i,j)     ! cld temp (K) at cld top
c                 write(6,*)"CTEM1",ctem(l,i,j),i,j,l,k
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
            end do
             end select

             select case (namedd(k))
              case ("CL3D")
              kunit=kunit+1
              do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=cl3d(l,i,j)     ! cld LWC (kg m-3)
c                write(6,*)"CL3D",cl3d(l,i,j),i,j,l,k
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
            end do
             end select

             select case (namedd(k))
              case ("CI3D")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=ci3d(l,i,j)     ! cld IWC (kg m-3)
c                write(6,*)"CI3D",ci3d(l,i,j),i,j,l,k
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
            end do
            end select

            select case (namedd(k))
             case ("CD3D")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                  data(i,j)=cd3d(l,i,j)     ! cld thickness (m)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
            end do
             end select

             select case (namedd(k))
              case ("CLDSS")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                  data(i,j)=100.d0*cldss(l,i,j)    ! Cld cover LS(%)
c                 if(cldss(l,i,j).gt.0.d0)
c    &            write(6,*)"CLDSS",cldss(l,i,j),i,j,l,k
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
             end select

             select case (namedd(k))
              case ("CLDMC")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=100.d0*cldmc(l,i,j)   ! Cld cover MC(%)
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
             end select

             select case (namedd(k))
              case ("CDN3D")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=cdn3d(l,i,j)    ! cld CDNC (cm^-3)
c        if(cdn3d(l,i,j).gt.20.d0)write(6,*)"CDN",cdn3d(l,i,j),i,j,l,k
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
            end do
             end select

             select case (namedd(k))
              case ("CRE3D")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=1.d-6*cre3d(l,i,j)    ! cld Reff (m)
c        if(cre3d(l,i,j).gt.2.d0)write(6,*)"CRe",cre3d(l,i,j),i,j,l,k
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
            end do
             end select

             select case (namedd(k))
              case ("TAUSS")
            kunit=kunit+1
            do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=tauss(l,i,j)       ! LS cld tau
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
           end do
             end select

            select case (namedd(k))
              case ("TAUMC")
            kunit=kunit+1
           do l=1,lm
              do j=1,jm
                do i=1,imaxj(j)
                 data(i,j)=taumc(l,i,j)       ! MC cld tau
                end do
              end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
             call writei(iu_subdd(kunit),itime,data,im*jm)
          end do
             end select

           select case (namedd(k))  !(1:5))
             case ("CLWP")         !LWP (kg m-2)
            kunit=kunit+1
             do j=1,jm
                do i=1,imaxj(j)
                  data(i,j)=clwp(i,j)
c       if(clwp(i,j).ne.0.) write(6,*)"CLWP",data(i,j),i,j,k,kunit
                end do
             end do
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
C**** write out
            call writei(iu_subdd(kunit),itime,data,im*jm)
            end select
#endif
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
        SELECT CASE (namedd(k))
        CASE ('DUEMIS')       ! Dust emission flux [kg/m^2/s]
          kunit=kunit+1
          DO n=1,Ntm_dust
            data(:,:)=dust_flux_glob(:,:,n)
C**** fix polar values
            data(2:im,1) =data(1,1)
            data(2:im,jm)=data(1,jm)
            !CALL writei(iu_subdd(kunit),itime,data,im*jm)
            call writei_parallel(grid,
     &           iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
          END DO
#ifdef TRACERS_DUST
        CASE ('DUEMIS2')   ! Dust emission flux 2 (diag. var. only) [kg/m^2/s]
          kunit=kunit+1
          DO n=1,Ntm_dust
            data(:,:)=dust_flux2_glob(:,:,n)
C**** fix polar values
            data(2:im,1) =data(1,1)
            data(2:im,jm)=data(1,jm)
            !CALL writei(iu_subdd(kunit),itime,data,im*jm)
            call writei_parallel(grid,
     &           iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
          END DO
#endif
#ifdef TRACERS_DRYDEP
        CASE ('DUDEPTURB')  ! Turb. deposition flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          DO n=1,Ntm_dust
            n1=n_fidx+n-1
            IF (dodrydep(n1)) THEN
              data(:,:)=depo_turb_glob(:,:,1,n1)+
     &             depo_turb_glob(:,:,2,n1)+depo_turb_glob(:,:,3,n1)
     &             +depo_turb_glob(:,:,4,n1)
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
              !CALL writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
            END IF
          END DO
        CASE ('DUDEPGRAV')  ! Gravit. settling flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          DO n=1,Ntm_dust
            n1=n_fidx+n-1
            IF (dodrydep(n1)) THEN
              data(:,:)=depo_grav_glob(:,:,1,n1)
     &             +depo_grav_glob(:,:,2,n1)+depo_grav_glob(:,:,3,n1)
     &             +depo_grav_glob(:,:,4,n1)
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
              CALL writei(iu_subdd(kunit),itime,data,im*jm)
            END IF
          END DO
#endif
          CASE ('DUDEPWET')   ! Wet deposition flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          DO n=1,Ntm_dust
#ifdef TRACERS_WATER
            n1=n_fidx+n-1
            IF (dowetdep(n1)) THEN
              DO j=1,Jm
                data(:,j)=trprec(n1,:,j)*bydxyp(j)/Dtsrc
              END DO
#else
              DO j=1,Jm
                data(:,j)=trprec_dust(n,:,j)*bydxyp(j)/Dtsrc
              END DO
#endif
C**** fix polar values
              data(2:im,1) =data(1,1)
              data(2:im,jm)=data(1,jm)
              !CALL writei(iu_subdd(kunit),itime,data,im*jm)
              call writei_parallel(grid,
     &             iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
#ifdef TRACERS_WATER
            END IF
#endif
          END DO
        CASE ('DUTRS')       ! Mixing ratio of dust tracers at surface [kg/kg]
          kunit=kunit+1
          DO n=1,Ntm_dust
            n1=n_fidx+n-1
            data(:,:)=trs_glob(:,:,1,n1)+trs_glob(:,:,2,n1)
     &           +trs_glob(:,:,3,n1)+trs_glob(:,:,4,n1)
C**** fix polar values
            data(2:im,1) =data(1,1)
            data(2:im,jm)=data(1,jm)
            !CALL writei(iu_subdd(kunit),itime,data,im*jm)
            call writei_parallel(grid,
     &           iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
          END DO
          CASE ('DULOAD')       ! Dust load [kg/m^2]
          kunit=kunit+1
          DO n=1,Ntm_dust
            n1=n_fidx+n-1
            data(:,:)=0.D0
            DO j=1,Jm
              DO l=1,LmaxSUBDD
                data(:,j)=data(:,j)+trm(:,j,l,n1)
              END DO
              data(:,j)=data(:,j)*bydxyp(j)
            END DO
C**** fix polar values
            data(2:im,1) =data(1,1)
            data(2:im,jm)=data(1,jm)
            !CALL writei(iu_subdd(kunit),itime,data,im*jm)
            call writei_parallel(grid,
     &           iu_subdd(kunit),nameunit(iu_subdd(kunit)),data,itime)
          END DO
        END SELECT
#endif
      end do ! end of  k=1,kdd loop
c****
      return
      end subroutine get_subdd

      end module subdaily

      subroutine ahourly
!@sum ahourly saves instantaneous variables at sub-daily frequency
!@+   for diurnal cycle diagnostics
!@auth Reha Cakmur/Jan Perlwitz

      USE MODEL_COM, only : u,v,t,p,q,jhour,ptop,sig
      USE CONSTANT, only : bygrav
      USE domain_decomp,ONLY : grid,get
      USE GEOM, only : imaxj,dxyp,bydxyp
      USE DYNAMICS, only : phi,wsave,pek,byam
      USE rad_com,ONLY : cosz1,srnflb_save,trnflb_save,ttausv_save,
     &     ttausv_cs_save
      USE diag_com,ONLY : adiurn_dust,ndiupt,ijdd,adiurn
#ifdef TRACERS_DUST
     *     ,idd_u1,idd_u2,idd_u3,idd_u4,idd_u5,idd_u6,idd_u7
     *     ,idd_u8,idd_u9,idd_u10,idd_u11,idd_v1,idd_v2,idd_v3,idd_v4
     *     ,idd_v5,idd_v6,idd_v7,idd_v8,idd_v9,idd_v10,idd_v11
     *     ,idd_uv1,idd_uv2,idd_uv3,idd_uv4,idd_uv5,idd_uv6,idd_uv7
     *     ,idd_uv8,idd_uv9,idd_uv10,idd_uv11,idd_t1,idd_t2,idd_t3
     *     ,idd_t4,idd_t5,idd_t6,idd_t7,idd_t8,idd_t9,idd_t10,idd_t11
     *     ,idd_qq1,idd_qq2,idd_qq3,idd_qq4,idd_qq5,idd_qq6,idd_qq7
     *     ,idd_qq8,idd_qq9,idd_qq10,idd_qq11,idd_p1,idd_p2,idd_p3
     *     ,idd_p4,idd_p5,idd_p6,idd_p7,idd_p8,idd_p9,idd_p10,idd_p11
     *     ,idd_w1,idd_w2,idd_w3,idd_w4,idd_w5,idd_w6,idd_w7,idd_w8
     *     ,idd_w9,idd_w10,idd_w11,idd_phi1,idd_phi2,idd_phi3,idd_phi4
     *     ,idd_phi5,idd_phi6,idd_phi7,idd_phi8,idd_phi9,idd_phi10
     *     ,idd_phi11,idd_sr1,idd_sr2
     *     ,idd_sr3,idd_sr4,idd_sr5,idd_sr6,idd_sr7,idd_sr8,idd_sr9
     *     ,idd_sr10,idd_sr11,idd_tr1,idd_tr2,idd_tr3,idd_tr4
     *     ,idd_tr5,idd_tr6,idd_tr7,idd_tr8,idd_tr9,idd_tr10,idd_tr11
     *     ,idd_load1,idd_load2,idd_load3,idd_load4,idd_load5
     *     ,idd_load6,idd_load7,idd_load8,idd_load9,idd_load10
     *     ,idd_load11,idd_conc1,idd_conc2,idd_conc3,idd_conc4
     *     ,idd_conc5,idd_conc6,idd_conc7,idd_conc8,idd_conc9
     *     ,idd_conc10,idd_conc11
     *     ,idd_tau1,idd_tau2,idd_tau3,idd_tau4,idd_tau5,idd_tau6
     *     ,idd_tau7,idd_tau8,idd_tau9,idd_tau10,idd_tau11
     *     ,idd_tau_cs1,idd_tau_cs2,idd_tau_cs3,idd_tau_cs4
     *     ,idd_tau_cs5,idd_tau_cs6,idd_tau_cs7,idd_tau_cs8
     *     ,idd_tau_cs9,idd_tau_cs10,idd_tau_cs11
#endif
#ifdef TRACERS_ON
      USE TRACER_COM, only : trm
#ifdef TRACERS_DUST
     &     ,Ntm_dust,n_clay
#endif
#endif
      INTEGER :: i,j,ih,kr,n,n1,nx
      REAL*8 :: psk

C****   define local grid
      INTEGER J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL get(grid, J_STRT=J_0, J_STOP=J_1)

      ih=jhour+1
!$OMP PARALLEL DO PRIVATE(i,j,kr,n,psk,n1)
!$OMP*   SCHEDULE(DYNAMIC,2)
      do j=j_0,j_1
      do i=1,imaxj(j)
      psk=pek(1,i,j)
      do kr=1,ndiupt
        if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then
#ifdef TRACERS_DUST
          IF (adiurn_dust == 1) THEN

            adiurn(ih,idd_u1,kr)=adiurn(ih,idd_u1,kr)+u(i,j,1)
            adiurn(ih,idd_u2,kr)=adiurn(ih,idd_u2,kr)+u(i,j,2)
            adiurn(ih,idd_u3,kr)=adiurn(ih,idd_u3,kr)+u(i,j,3)
            adiurn(ih,idd_u4,kr)=adiurn(ih,idd_u4,kr)+u(i,j,4)
            adiurn(ih,idd_u5,kr)=adiurn(ih,idd_u5,kr)+u(i,j,5)
            adiurn(ih,idd_u6,kr)=adiurn(ih,idd_u6,kr)+u(i,j,6)
            adiurn(ih,idd_u7,kr)=adiurn(ih,idd_u7,kr)+u(i,j,7)
            adiurn(ih,idd_u8,kr)=adiurn(ih,idd_u8,kr)+u(i,j,8)
            adiurn(ih,idd_u9,kr)=adiurn(ih,idd_u9,kr)+u(i,j,9)
            adiurn(ih,idd_u10,kr)=adiurn(ih,idd_u10,kr)+u(i,j,10)
            adiurn(ih,idd_u11,kr)=adiurn(ih,idd_u11,kr)+u(i,j,11)

            adiurn(ih,idd_v1,kr)=adiurn(ih,idd_v1,kr)+v(i,j,1)
            adiurn(ih,idd_v2,kr)=adiurn(ih,idd_v2,kr)+v(i,j,2)
            adiurn(ih,idd_v3,kr)=adiurn(ih,idd_v3,kr)+v(i,j,3)
            adiurn(ih,idd_v4,kr)=adiurn(ih,idd_v4,kr)+v(i,j,4)
            adiurn(ih,idd_v5,kr)=adiurn(ih,idd_v5,kr)+v(i,j,5)
            adiurn(ih,idd_v6,kr)=adiurn(ih,idd_v6,kr)+v(i,j,6)
            adiurn(ih,idd_v7,kr)=adiurn(ih,idd_v7,kr)+v(i,j,7)
            adiurn(ih,idd_v8,kr)=adiurn(ih,idd_v8,kr)+v(i,j,8)
            adiurn(ih,idd_v9,kr)=adiurn(ih,idd_v9,kr)+v(i,j,9)
            adiurn(ih,idd_v10,kr)=adiurn(ih,idd_v10,kr)+v(i,j,10)
            adiurn(ih,idd_v11,kr)=adiurn(ih,idd_v11,kr)+v(i,j,11)

            adiurn(ih,idd_uv1,kr)=adiurn(ih,idd_uv1,kr)
     *           +sqrt( u(i,j,1)*u(i,j,1) + v(i,j,1)*v(i,j,1))
            adiurn(ih,idd_uv2,kr)=adiurn(ih,idd_uv2,kr)
     *           +sqrt( u(i,j,2)*u(i,j,2) + v(i,j,2)*v(i,j,2))
            adiurn(ih,idd_uv3,kr)=adiurn(ih,idd_uv3,kr)
     *           +sqrt( u(i,j,3)*u(i,j,3) + v(i,j,3)*v(i,j,3))
            adiurn(ih,idd_uv4,kr)=adiurn(ih,idd_uv4,kr)
     *           +sqrt( u(i,j,4)*u(i,j,4) + v(i,j,4)*v(i,j,4))
            adiurn(ih,idd_uv5,kr)=adiurn(ih,idd_uv5,kr)
     *           +sqrt( u(i,j,5)*u(i,j,5) + v(i,j,5)*v(i,j,5))
            adiurn(ih,idd_uv6,kr)=adiurn(ih,idd_uv6,kr)
     *           +sqrt( u(i,j,6)*u(i,j,6) + v(i,j,6)*v(i,j,6))
            adiurn(ih,idd_uv7,kr)=adiurn(ih,idd_uv7,kr)
     *           +sqrt( u(i,j,7)*u(i,j,7) + v(i,j,7)*v(i,j,7))
            adiurn(ih,idd_uv8,kr)=adiurn(ih,idd_uv8,kr)
     *           +sqrt( u(i,j,8)*u(i,j,8) + v(i,j,8)*v(i,j,8))
            adiurn(ih,idd_uv9,kr)=adiurn(ih,idd_uv9,kr)
     *           +sqrt( u(i,j,9)*u(i,j,9) + v(i,j,9)*v(i,j,9))
            adiurn(ih,idd_uv10,kr)=adiurn(ih,idd_uv10,kr)
     *           +sqrt( u(i,j,10)*u(i,j,10) + v(i,j,10)*v(i,j,10))
            adiurn(ih,idd_uv11,kr)=adiurn(ih,idd_uv11,kr)
     *           +sqrt( u(i,j,11)*u(i,j,11) + v(i,j,11)*v(i,j,11))

            adiurn(ih,idd_t1,kr)=adiurn(ih,idd_t1,kr)+t(i,j,1)*psk
            adiurn(ih,idd_t2,kr)=adiurn(ih,idd_t2,kr)+t(i,j,2)*psk
            adiurn(ih,idd_t3,kr)=adiurn(ih,idd_t3,kr)+t(i,j,3)*psk
            adiurn(ih,idd_t4,kr)=adiurn(ih,idd_t4,kr)+t(i,j,4)*psk
            adiurn(ih,idd_t5,kr)=adiurn(ih,idd_t5,kr)+t(i,j,5)*psk
            adiurn(ih,idd_t6,kr)=adiurn(ih,idd_t6,kr)+t(i,j,6)*psk
            adiurn(ih,idd_t7,kr)=adiurn(ih,idd_t7,kr)+t(i,j,7)*psk
            adiurn(ih,idd_t8,kr)=adiurn(ih,idd_t8,kr)+t(i,j,8)*psk
            adiurn(ih,idd_t9,kr)=adiurn(ih,idd_t9,kr)+t(i,j,9)*psk
            adiurn(ih,idd_t10,kr)=adiurn(ih,idd_t10,kr)+t(i,j,10)*psk
            adiurn(ih,idd_t11,kr)=adiurn(ih,idd_t11,kr)+t(i,j,11)*psk

            adiurn(ih,idd_qq1,kr)=adiurn(ih,idd_qq1,kr)+q(i,j,1)
            adiurn(ih,idd_qq2,kr)=adiurn(ih,idd_qq2,kr)+q(i,j,2)
            adiurn(ih,idd_qq3,kr)=adiurn(ih,idd_qq3,kr)+q(i,j,3)
            adiurn(ih,idd_qq4,kr)=adiurn(ih,idd_qq4,kr)+q(i,j,4)
            adiurn(ih,idd_qq5,kr)=adiurn(ih,idd_qq5,kr)+q(i,j,5)
            adiurn(ih,idd_qq6,kr)=adiurn(ih,idd_qq6,kr)+q(i,j,6)
            adiurn(ih,idd_qq7,kr)=adiurn(ih,idd_qq7,kr)+q(i,j,7)
            adiurn(ih,idd_qq8,kr)=adiurn(ih,idd_qq8,kr)+q(i,j,8)
            adiurn(ih,idd_qq9,kr)=adiurn(ih,idd_qq9,kr)+q(i,j,9)
            adiurn(ih,idd_qq10,kr)=adiurn(ih,idd_qq10,kr)+q(i,j,10)
            adiurn(ih,idd_qq11,kr)=adiurn(ih,idd_qq11,kr)+q(i,j,11)

            adiurn(ih,idd_p1,kr)=adiurn(ih,idd_p1,kr)
     *           +p(i,j)*sig(1)+ptop
            adiurn(ih,idd_p2,kr)=adiurn(ih,idd_p2,kr)
     *           +p(i,j)*sig(2)+ptop
            adiurn(ih,idd_p3,kr)=adiurn(ih,idd_p3,kr)
     *           +p(i,j)*sig(3)+ptop
            adiurn(ih,idd_p4,kr)=adiurn(ih,idd_p4,kr)
     *           +p(i,j)*sig(4)+ptop
            adiurn(ih,idd_p5,kr)=adiurn(ih,idd_p5,kr)
     *           +p(i,j)*sig(5)+ptop
            adiurn(ih,idd_p6,kr)=adiurn(ih,idd_p6,kr)
     *           +p(i,j)*sig(6)+ptop
            adiurn(ih,idd_p7,kr)=adiurn(ih,idd_p7,kr)
     *           +p(i,j)*sig(7)+ptop
            adiurn(ih,idd_p8,kr)=adiurn(ih,idd_p8,kr)
     *           +p(i,j)*sig(8)+ptop
            adiurn(ih,idd_p9,kr)=adiurn(ih,idd_p9,kr)
     *           +p(i,j)*sig(9)+ptop
            adiurn(ih,idd_p10,kr)=adiurn(ih,idd_p10,kr)
     *           +p(i,j)*sig(10)+ptop
            adiurn(ih,idd_p11,kr)=adiurn(ih,idd_p11,kr)
     *           +p(i,j)*sig(11)+ptop

            adiurn(ih,idd_w1,kr)=adiurn(ih,idd_w1,kr)+wsave(i,j,1)
            adiurn(ih,idd_w2,kr)=adiurn(ih,idd_w2,kr)+wsave(i,j,2)
            adiurn(ih,idd_w3,kr)=adiurn(ih,idd_w3,kr)+wsave(i,j,3)
            adiurn(ih,idd_w4,kr)=adiurn(ih,idd_w4,kr)+wsave(i,j,4)
            adiurn(ih,idd_w5,kr)=adiurn(ih,idd_w5,kr)+wsave(i,j,5)
            adiurn(ih,idd_w6,kr)=adiurn(ih,idd_w6,kr)+wsave(i,j,6)
            adiurn(ih,idd_w7,kr)=adiurn(ih,idd_w7,kr)+wsave(i,j,7)
            adiurn(ih,idd_w8,kr)=adiurn(ih,idd_w8,kr)+wsave(i,j,8)
            adiurn(ih,idd_w9,kr)=adiurn(ih,idd_w9,kr)+wsave(i,j,9)
            adiurn(ih,idd_w10,kr)=adiurn(ih,idd_w10,kr)+wsave(i,j,10)
            adiurn(ih,idd_w11,kr)=adiurn(ih,idd_w11,kr)+wsave(i,j,11)

            adiurn(ih,idd_phi1,kr)=adiurn(ih,idd_phi1,kr)
     *           +phi(i,j,1)*bygrav
            adiurn(ih,idd_phi2,kr)=adiurn(ih,idd_phi2,kr)
     *           +phi(i,j,2)*bygrav
            adiurn(ih,idd_phi3,kr)=adiurn(ih,idd_phi3,kr)
     *           +phi(i,j,3)*bygrav
            adiurn(ih,idd_phi4,kr)=adiurn(ih,idd_phi4,kr)
     *           +phi(i,j,4)*bygrav
            adiurn(ih,idd_phi5,kr)=adiurn(ih,idd_phi5,kr)
     *           +phi(i,j,5)*bygrav
            adiurn(ih,idd_phi6,kr)=adiurn(ih,idd_phi6,kr)
     *           +phi(i,j,6)*bygrav
            adiurn(ih,idd_phi7,kr)=adiurn(ih,idd_phi7,kr)
     *           +phi(i,j,7)*bygrav
            adiurn(ih,idd_phi8,kr)=adiurn(ih,idd_phi8,kr)
     *           +phi(i,j,8)*bygrav
            adiurn(ih,idd_phi9,kr)=adiurn(ih,idd_phi9,kr)
     *           +phi(i,j,9)*bygrav
            adiurn(ih,idd_phi10,kr)=adiurn(ih,idd_phi10,kr)
     *           +phi(i,j,10)*bygrav
            adiurn(ih,idd_phi11,kr)=adiurn(ih,idd_phi11,kr)
     *           +phi(i,j,11)*bygrav

            adiurn(ih,idd_sr1,kr)=adiurn(ih,idd_sr1,kr)
     *           +srnflb_save(i,j,1)*cosz1(i,j)
            adiurn(ih,idd_sr2,kr)=adiurn(ih,idd_sr2,kr)
     *           +srnflb_save(i,j,2)*cosz1(i,j)
            adiurn(ih,idd_sr3,kr)=adiurn(ih,idd_sr3,kr)
     *           +srnflb_save(i,j,3)*cosz1(i,j)
            adiurn(ih,idd_sr4,kr)=adiurn(ih,idd_sr4,kr)
     *           +srnflb_save(i,j,4)*cosz1(i,j)
            adiurn(ih,idd_sr5,kr)=adiurn(ih,idd_sr5,kr)
     *           +srnflb_save(i,j,5)*cosz1(i,j)
            adiurn(ih,idd_sr6,kr)=adiurn(ih,idd_sr6,kr)
     *           +srnflb_save(i,j,6)*cosz1(i,j)
            adiurn(ih,idd_sr7,kr)=adiurn(ih,idd_sr7,kr)
     *           +srnflb_save(i,j,7)*cosz1(i,j)
            adiurn(ih,idd_sr8,kr)=adiurn(ih,idd_sr8,kr)
     *           +srnflb_save(i,j,8)*cosz1(i,j)
            adiurn(ih,idd_sr9,kr)=adiurn(ih,idd_sr9,kr)
     *           +srnflb_save(i,j,9)*cosz1(i,j)
            adiurn(ih,idd_sr10,kr)=adiurn(ih,idd_sr10,kr)
     *           +srnflb_save(i,j,10)*cosz1(i,j)
            adiurn(ih,idd_sr11,kr)=adiurn(ih,idd_sr11,kr)
     *           +srnflb_save(i,j,11)*cosz1(i,j)

            adiurn(ih,idd_tr1,kr)=adiurn(ih,idd_tr1,kr)
     *           +trnflb_save(i,j,1)
            adiurn(ih,idd_tr2,kr)=adiurn(ih,idd_tr2,kr)
     *           +trnflb_save(i,j,2)
            adiurn(ih,idd_tr3,kr)=adiurn(ih,idd_tr3,kr)
     *           +trnflb_save(i,j,3)
            adiurn(ih,idd_tr4,kr)=adiurn(ih,idd_tr4,kr)
     *           +trnflb_save(i,j,4)
            adiurn(ih,idd_tr5,kr)=adiurn(ih,idd_tr5,kr)
     *           +trnflb_save(i,j,5)
            adiurn(ih,idd_tr6,kr)=adiurn(ih,idd_tr6,kr)
     *           +trnflb_save(i,j,6)
            adiurn(ih,idd_tr7,kr)=adiurn(ih,idd_tr7,kr)
     *           +trnflb_save(i,j,7)
            adiurn(ih,idd_tr8,kr)=adiurn(ih,idd_tr8,kr)
     *           +trnflb_save(i,j,8)
            adiurn(ih,idd_tr9,kr)=adiurn(ih,idd_tr9,kr)
     *           +trnflb_save(i,j,9)
            adiurn(ih,idd_tr10,kr)=adiurn(ih,idd_tr10,kr)
     *           +trnflb_save(i,j,10)
            adiurn(ih,idd_tr11,kr)=adiurn(ih,idd_tr11,kr)
     *           +trnflb_save(i,j,11)

            DO n=1,Ntm_dust
              n1=n_clay+n-1

              adiurn(ih,idd_load1,kr)=adiurn(ih,idd_load1,kr)
     *             +trm(i,j,1,n1)/dxyp(j)
              adiurn(ih,idd_load2,kr)=adiurn(ih,idd_load2,kr)
     *             +trm(i,j,2,n1)/dxyp(j)
              adiurn(ih,idd_load3,kr)=adiurn(ih,idd_load3,kr)
     *             +trm(i,j,3,n1)/dxyp(j)
              adiurn(ih,idd_load4,kr)=adiurn(ih,idd_load4,kr)
     *             +trm(i,j,4,n1)/dxyp(j)
              adiurn(ih,idd_load5,kr)=adiurn(ih,idd_load5,kr)
     *             +trm(i,j,5,n1)/dxyp(j)
              adiurn(ih,idd_load6,kr)=adiurn(ih,idd_load6,kr)
     *             +trm(i,j,6,n1)/dxyp(j)
              adiurn(ih,idd_load7,kr)=adiurn(ih,idd_load7,kr)
     *             +trm(i,j,7,n1)/dxyp(j)
              adiurn(ih,idd_load8,kr)=adiurn(ih,idd_load8,kr)
     *             +trm(i,j,8,n1)/dxyp(j)
              adiurn(ih,idd_load9,kr)=adiurn(ih,idd_load9,kr)
     *             +trm(i,j,9,n1)/dxyp(j)
              adiurn(ih,idd_load10,kr)=adiurn(ih,idd_load10,kr)
     *             +trm(i,j,10,n1)/dxyp(j)
              adiurn(ih,idd_load11,kr)=adiurn(ih,idd_load11,kr)
     *             +trm(i,j,11,n1)/dxyp(j)

              adiurn(ih,idd_conc1,kr)=adiurn(ih,idd_conc1,kr)
     *             +trm(i,j,1,n1)*byam(1,i,j)*bydxyp(j)
              adiurn(ih,idd_conc2,kr)=adiurn(ih,idd_conc2,kr)
     *             +trm(i,j,2,n1)*byam(2,i,j)*bydxyp(j)
              adiurn(ih,idd_conc3,kr)=adiurn(ih,idd_conc3,kr)
     *             +trm(i,j,3,n1)*byam(3,i,j)*bydxyp(j)
              adiurn(ih,idd_conc4,kr)=adiurn(ih,idd_conc4,kr)
     *             +trm(i,j,4,n1)*byam(4,i,j)*bydxyp(j)
              adiurn(ih,idd_conc5,kr)=adiurn(ih,idd_conc5,kr)
     *             +trm(i,j,5,n1)*byam(5,i,j)*bydxyp(j)
              adiurn(ih,idd_conc6,kr)=adiurn(ih,idd_conc6,kr)
     *             +trm(i,j,6,n1)*byam(6,i,j)*bydxyp(j)
              adiurn(ih,idd_conc7,kr)=adiurn(ih,idd_conc7,kr)
     *             +trm(i,j,7,n1)*byam(7,i,j)*bydxyp(j)
              adiurn(ih,idd_conc8,kr)=adiurn(ih,idd_conc8,kr)
     *             +trm(i,j,8,n1)*byam(8,i,j)*bydxyp(j)
              adiurn(ih,idd_conc9,kr)=adiurn(ih,idd_conc9,kr)
     *             +trm(i,j,9,n1)*byam(9,i,j)*bydxyp(j)
              adiurn(ih,idd_conc10,kr)=adiurn(ih,idd_conc10,kr)
     *             +trm(i,j,10,n1)*byam(10,i,j)*bydxyp(j)
              adiurn(ih,idd_conc11,kr)=adiurn(ih,idd_conc11,kr)
     *             +trm(i,j,11,n1)*byam(11,i,j)*bydxyp(j)

              adiurn(ih,idd_tau1,kr)=adiurn(ih,idd_tau1,kr)
     *             +ttausv_save(i,j,n1,1)
              adiurn(ih,idd_tau2,kr)=adiurn(ih,idd_tau2,kr)
     *             +ttausv_save(i,j,n1,2)
              adiurn(ih,idd_tau3,kr)=adiurn(ih,idd_tau3,kr)
     *             +ttausv_save(i,j,n1,3)
              adiurn(ih,idd_tau4,kr)=adiurn(ih,idd_tau4,kr)
     *             +ttausv_save(i,j,n1,4)
              adiurn(ih,idd_tau5,kr)=adiurn(ih,idd_tau5,kr)
     *             +ttausv_save(i,j,n1,5)
              adiurn(ih,idd_tau6,kr)=adiurn(ih,idd_tau6,kr)
     *             +ttausv_save(i,j,n1,6)
              adiurn(ih,idd_tau7,kr)=adiurn(ih,idd_tau7,kr)
     *             +ttausv_save(i,j,n1,7)
              adiurn(ih,idd_tau8,kr)=adiurn(ih,idd_tau8,kr)
     *             +ttausv_save(i,j,n1,8)
              adiurn(ih,idd_tau9,kr)=adiurn(ih,idd_tau9,kr)
     *             +ttausv_save(i,j,n1,9)
              adiurn(ih,idd_tau10,kr)=adiurn(ih,idd_tau10,kr)
     *             +ttausv_save(i,j,n1,10)
              adiurn(ih,idd_tau11,kr)=adiurn(ih,idd_tau11,kr)
     *             +ttausv_save(i,j,n1,11)

              adiurn(ih,idd_tau_cs1,kr)=adiurn(ih,idd_tau_cs1,kr)
     *             +ttausv_cs_save(i,j,n1,1)
              adiurn(ih,idd_tau_cs2,kr)=adiurn(ih,idd_tau_cs2,kr)
     *             +ttausv_cs_save(i,j,n1,2)
              adiurn(ih,idd_tau_cs3,kr)=adiurn(ih,idd_tau_cs3,kr)
     *             +ttausv_cs_save(i,j,n1,3)
              adiurn(ih,idd_tau_cs4,kr)=adiurn(ih,idd_tau_cs4,kr)
     *             +ttausv_cs_save(i,j,n1,4)
              adiurn(ih,idd_tau_cs5,kr)=adiurn(ih,idd_tau_cs5,kr)
     *             +ttausv_cs_save(i,j,n1,5)
              adiurn(ih,idd_tau_cs6,kr)=adiurn(ih,idd_tau_cs6,kr)
     *             +ttausv_cs_save(i,j,n1,6)
              adiurn(ih,idd_tau_cs7,kr)=adiurn(ih,idd_tau_cs7,kr)
     *             +ttausv_cs_save(i,j,n1,7)
              adiurn(ih,idd_tau_cs8,kr)=adiurn(ih,idd_tau_cs8,kr)
     *             +ttausv_cs_save(i,j,n1,8)
              adiurn(ih,idd_tau_cs9,kr)=adiurn(ih,idd_tau_cs9,kr)
     *             +ttausv_cs_save(i,j,n1,9)
              adiurn(ih,idd_tau_cs10,kr)=adiurn(ih,idd_tau_cs10,kr)
     *             +ttausv_cs_save(i,j,n1,10)
              adiurn(ih,idd_tau_cs11,kr)=adiurn(ih,idd_tau_cs11,kr)
     *             +ttausv_cs_save(i,j,n1,11)

            END DO

          END IF
#endif
        endif
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine ahourly

      SUBROUTINE init_DIAG(istart,num_acc_files)
!@sum  init_DIAG initializes the diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday,kapa,undef
      USE MODEL_COM, only : lm,Itime,ItimeI,Itime0,sige,sig,ptop
     *     ,pmtop,psfmpt,nfiltr,jhour,jdate,jmon,amon,jyear
     *     ,jhour0,jdate0,jmon0,amon0,jyear0,idacc,ioread_single
     *     ,xlabel,iowrite_single,iyear1,nday,dtsrc,dt,nmonav
     *     ,ItimeE,lrunid,focean
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE DIAG_COM, only : TSFREZ => TSFREZ_loc
      USE DIAG_COM, only : NPTS, NAMDD, NDIUPT, IJDD, ISCCP_DIAGS
      USE DIAG_COM, only : monacc, acc_period, keyct, KEYNR, PLE
      USE DIAG_COM, only : PLM, p1000k, icon_AM, NOFM
      USE DIAG_COM, only : PLE_DN, icon_KE, NSUM_CON, IA_CON, SCALE_CON
      USE DIAG_COM, only : TITLE_CON, PSPEC, LSTR, NSPHER, KLAYER
      USE DIAG_COM, only : ISTRAT, kgz, pmb, kgz_max
      USE DIAG_COM, only : TF_DAY1, TF_LAST, TF_LKON, TF_LKOFF
      USE DIAG_COM, only : name_consrv, units_consrv, lname_consrv
      USE DIAG_COM, only : CONPT0, icon_MS, icon_TPE, icon_WM, icon_EWM
      USE diag_com,ONLY : adiurn_dust
      USE DIAG_LOC
      USE PARAM
      USE FILEMANAGER
      USE DOMAIN_DECOMP, only: GRID,GET,WRITE_PARALLEL
      IMPLICIT NONE
      integer, intent(in) :: istart,num_acc_files
      INTEGER I,J,L,K,KL,ioerr,months,years,mswitch,ldate,iu_AIC
     *     ,jday0,jday,moff,kb,iu_ACC,l850,l300,l50
      REAL*8 PLE_tmp
      CHARACTER FILENM*100
      CHARACTER CONPT(NPTS)*10
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER :: J_0,J_1
!@var out_line local variable to hold mixed-type output for parallel I/O
      character(len=300) :: out_line

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      call sync_param( "NAMDD", NAMDD, NDIUPT )
      call sync_param( "IJDD", IJDD(1:2,1), 2*NDIUPT )
      call sync_param( "isccp_diags",isccp_diags)
      call sync_param( "adiurn_dust",adiurn_dust)

      IF(ISTART.LT.1) THEN  ! initialize for post-processing
        call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jday0,Jdate0,Jhour0
     *       ,amon0)
        call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour
     *       ,amon)
        months=1 ; years=monacc(jmon0) ; mswitch=0 ; moff=0 ; kb=jmon0
        do kl=jmon0+1,jmon0+11
          k = kl
          if (k.gt.12) k=k-12
          if (monacc(k).eq.years) then
            months=months+1
          else if (monacc(k).ne.0) then
C****            write(6,*) 'uneven period:',monacc
            CALL WRITE_PARALLEL(monacc, UNIT=6, format=
     &                          "('uneven period:',12I5)")
            call stop_model( 'uneven period', 255 )
          end if
          if(monacc(k).ne.monacc(kb)) mswitch = mswitch+1
          if(mswitch.eq.2) moff = moff+1
          kb = k
        end do
        if (mswitch.gt.2) then
C****          write(6,*) 'non-consecutive period:',monacc
            CALL WRITE_PARALLEL(monacc, UNIT=6, format=
     &                          "('non-consecutive period:',12I5)")
          call stop_model( 'non-consecutive period', 255 )
        end if
        call aPERIOD (JMON0,JYEAR0,months,years,moff, acc_period,Ldate)
        if (num_acc_files.gt.1) then  ! save the summed acc-file
          write(out_line,*) num_acc_files,' files are summed up'
          CALL WRITE_PARALLEL(TRIM(out_line), UNIT=6)
          keyct=1 ; KEYNR=0
          XLABEL(128:132)='     '
          XLABEL(120:132)=acc_period(1:3)//' '//acc_period(4:Ldate)
C****          write(6,*) XLABEL
          CALL WRITE_PARALLEL(XLABEL, UNIT=6)
          call openunit(acc_period(1:Ldate)//'.acc'//XLABEL(1:LRUNID)
     *         ,iu_ACC,.true.,.false.)
          call io_rsf (iu_ACC,Itime,iowrite_single,ioerr)
          call closeunit(iu_ACC)
        end if
        ItimeE = -1
        close (6)
        open(6,file=acc_period(1:Ldate)//'.'//XLABEL(1:LRUNID)//'.PRT',
     *       FORM='FORMATTED')
      END IF

C**** Initialize certain arrays used by more than one print routine
      DO L=1,LM
        PLE(L)=SIGE(L+1)*PSFMPT+PTOP
        PLE_DN(L)=SIGE(L)*PSFMPT+PTOP
        PLM(L)=SIG(L)*PSFMPT+PTOP
      END DO
      PLM(LM+1)=.75d0*PMTOP
      PLM(LM+2)=.35d0*PMTOP
      PLM(LM+3)=.1d0*PMTOP
      p1000k=1000.0**kapa

C**** Initialise some local constants (replaces IFIRST constructions)
C**** From DIAGA:
      DO L=1,LM
        LUPA(L)=L+1
        LDNA(L)=L-1
      END DO
      LDNA(1)=1
      LUPA(LM)=LM

C**** From DIAGB
      PM(1)=1200.
      DO L=2,LM+1
        PL(L)=PSFMPT*SIGE(L)+PTOP
        PM(L)=PSFMPT*SIGE(L)+PTOP
      END DO
      DO L=1,LM
        PLO(L)=PSFMPT*SIG(L)+PTOP
        PMO(L)=.5*(PM(L)+PM(L+1))
      END DO

C**** From DIAG7A
      L850=LM
      L300=LM
      L50=LM
      DO L=LM-1,1,-1
        PLE_tmp=.25*(SIGE(L)+2.*SIGE(L+1)+SIGE(L+2))*PSFMPT+PTOP
        IF (PLE_tmp.LT.850.) L850=L
        IF (PLE_tmp.LT.300.) L300=L
        IF (PLE_tmp.LT.250.) JET=L
        IF (PLE_tmp.LT.50.) L50=L
      END DO
C      WRITE (6,888) JET
      CALL WRITE_PARALLEL(JET, UNIT=6, format=
     & "(' JET WIND LEVEL FOR DIAG',I3)")
C 888  FORMAT (' JET WIND LEVEL FOR DIAG',I3)
C****      WRITE (6,889) L850,L300,L50
      WRITE (out_line,889) L850,L300,L50
      CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
 889  FORMAT (' LEVELS FOR WIND WAVE POWER DIAG  L850=',I3,
     *     ' L300=',I3,' L50=',I3)
      LDEX(1)=L850
      LDEX(2)=L300
      LDEX(3)=L50

C**** Initialize conservation diagnostics
C**** NCON=1:25 are special cases: Angular momentum and kinetic energy
      icon_AM=1
      NOFM(:,icon_AM) = (/  1, 8, 0, 0, 0, 0, 9,10, 0,11, 0, 0/)
      icon_KE=2
      NOFM(:,icon_KE) = (/ 13,20,21, 0, 0, 0,22,23, 0,24, 0, 0/)
      NSUM_CON(1:25) = (/-1,-1,-1,-1,-1,-1,-1,12,12,12,12, 0,
     *                   -1,-1,-1,-1,-1,-1,-1,25,25,25,25,25, 0/)
      IA_CON(1:25) =   (/12, 1, 1, 1, 1, 1, 1, 7, 8,10, 9,12,
     *                   12, 1, 1, 1, 1, 1, 1, 7, 8, 8,10, 9,12/)
      SCALE_CON(1)              = 1d-9
      SCALE_CON((/2,3,4,5,6,7,8,9/))= 1d-2/DTSRC
      SCALE_CON(10)              = 1d-2/(NFILTR*DTSRC)
      SCALE_CON(11)             = 2d-2/SDAY
      SCALE_CON((/12,25/))      = 1.
      SCALE_CON(13)             = 1d-3
      SCALE_CON((/14,15,16,17,18,19,20,21,22/)) = 1d3/DTSRC
      SCALE_CON(23)             = 1d3/(NFILTR*DTSRC)
      SCALE_CON(24)             = 2d3/SDAY
      TITLE_CON(1:25) = (/
     *  ' INSTANTANE AM (10**9 J*S/M^2)  ',
     *  '     DELTA AM BY ADVECTION      ',
     *  '     DELTA AM BY CORIOLIS FORCE ',
     *  '     DELTA AM BY PRESSURE GRAD  ',
     *  '     DELTA AM BY STRATOS DRAG   ',
     *  '     DELTA AM BY UV FILTER      ',
     *  '     DELTA AM BY GW DRAG        ',
     *  ' CHANGE OF AM BY DYNAMICS       ',
     *  ' CHANGE OF AM BY SURF FRIC+TURB ',
     *  ' CHANGE OF AM BY FILTER         ',
     *  ' CHANGE OF AM BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**2 J/M^2)   ',
     *  '0INSTANTANEOUS KE (10**3 J/M^2) ',
     *  '     DELTA KE BY ADVECTION      ',
     *  '     DELTA KE BY CORIOLIS FORCE ',
     *  '     DELTA KE BY PRESSURE GRAD  ',
     *  '     DELTA KE BY STRATOS DRAG   ',
     *  '     DELTA KE BY UV FILTER      ',
     *  '     DELTA KE BY GW DRAG        ',
     *  ' CHANGE OF KE BY DYNAMICS       ',
     *  ' CHANGE OF KE BY MOIST CONVEC   ',
     *  ' CHANGE OF KE BY SURF + DC/TURB ',
     *  ' CHANGE OF KE BY FILTER         ',
     *  ' CHANGE OF KE BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**-3 W/M^2)  '/)
      name_consrv(1:25) = (/
     *     'inst_AM   ','del_AM_ADV','del_AM_COR','del_AM_PRE',
     *     'del_AM_STR','del_AM_UVF','del_AM_GWD','chg_AM_DYN'
     *     ,'chg_AM_SUR','chg_AM_FIL','chg_AM_DAI','sum_chg_AM'
     *     ,'inst_KE   ','del_KE_ADV','del_KE_COR','del_KE_PRE'
     *     ,'del_KE_STR','del_KE_UVF','del_KE_GWD','chg_KE_DYN'
     *     ,'chg_KE_MOI','chg_KE_SUR','del_KE_FIL','chg_KE_DAI'
     *     ,'sum_chg_KE'/)
      units_consrv(1)    ="10**9 J*S/M^2"
      units_consrv(2:12) ="10**2 J/M^2"
      units_consrv(13)   ="10**3 J/M^2"
      units_consrv(14:24)="10**-3 W/M^2"
      lname_consrv(1:25)=TITLE_CON(1:25)
C**** To add a new conservation diagnostic:
C****    i) Add 1 to NQUANT, and increase KCON in DIAG_COM.f
C****   ii) Set up a QCON, and call SET_CON to allocate array numbers,
C****       set up scales, titles, etc. The icon_XX index must be
C****       declared in DIAG_COM.f for the time being
C**** QCON denotes when the conservation diags should be done
C**** 1:NPTS ==> DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C****  iii) Write a conserv_XYZ routine that returns the zonal average
C****       of your quantity
C****   iv) Add a line to DIAGCA that calls conserv_DIAG (declared
C****       as external)
C****    v) Note that the conserv_XYZ routine, and call to SET_CON
C****       should be in the driver module for the relevant physics

C**** Set up atmospheric component conservation diagnostics
      CONPT=CONPT0
C**** Atmospheric mass
      QCON=(/ T, F, F, F, F, F, T, F, T, F, F/)
      CALL SET_CON(QCON,CONPT,"MASS    ","(KG/M^2)        ",
     *     "(10**-8 KG/SM^2)",1d0,1d8,icon_MS)
C**** Atmospheric total potential energy
      CONPT(8)="SURF+TURB" ; CONPT(6)="KE DISSIP"
      QCON=(/ T, T, T, F, F, T, T, T, F, F, F/)
      CALL SET_CON(QCON,CONPT,"TPE     ","(10**5 J/M^2)   ",
     *     "(10**-2 W/M^2)  ",1d-5,1d2,icon_TPE)
C**** Atmospheric water mass
      CONPT(6)="SURF+TURB"
      QCON=(/ T, T, F, F, F, T, F, F, F, F, F/)
      CALL SET_CON(QCON,CONPT,"ATM WAT ","(10**-2 KG/M^2) ",
     *     "(10**-8 KG/SM^2)",1d2,1d8,icon_WM)
C**** Atmospheric water latent heat (at some point should include
C**** sensible + potential energy associated with water mass as well)
      QCON=(/ T, T, F, F, F, T, F, F, F, F, F/)
      CALL SET_CON(QCON,CONPT,"ENRG WAT","(10**3 J/M^2)   ",
     *     "(10**-2 W/M^2)  ",1d-3,1d2,icon_EWM)

C**** Initialize layering for spectral diagnostics
C**** add in epsilon=1d-5 to avoid roundoff mistakes
      KL=1
      DO L=1,LM
        IF (PTOP+PSFMPT*SIGE(L+1)+1d-5.lt.PSPEC(KL) .and.
     *      PTOP+PSFMPT*SIGE(L)+1d-5.gt.PSPEC(KL)) THEN
          IF (KL.eq.2) LSTR = L  ! approx. 10mb height
          KL=KL+1
        END IF
        KLAYER(L)=4*(KL-1)+1
      END DO
      IF (KL*4 .gt. NSPHER) THEN
C****        WRITE(6,*) "Inconsistent definitions of stratosphere:"
        CALL WRITE_PARALLEL("Inconsistent definitions of stratosphere:"
     &                       ,UNIT=6)
C****        WRITE(6,*) "Adjust PSPEC, ISTRAT so that KL*4 = NSPHER"
        CALL WRITE_PARALLEL("Adjust PSPEC, ISTRAT so that KL*4 = NSPHER"
     &                       ,UNIT=6)
C****        WRITE(6,*) "ISTRAT,PSPEC,NSPHER,KL=",ISTRAT,PSPEC,NSPHER,KL
        WRITE(out_line,*) "ISTRAT,PSPEC,NSPHER,KL=",
     &                     ISTRAT,PSPEC,NSPHER,KL
        CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        call stop_model(
     *    "Stratospheric definition problem for spectral diags.",255)
      END IF

C**** Calculate the max number of geopotential heights
      do k=1,kgz
        if (pmb(k).le.pmtop) exit
        kgz_max = k
      end do
C****      write(6,'(a)') " Geopotential height diagnostics at (mb): "
      CALL WRITE_PARALLEL(" Geopotential height diagnostics at (mb): ",
     &                      UNIT=6)
C*****      write(6,'(20F9.3)') PMB(1:kgz_max)
      CALL WRITE_PARALLEL(PMB(1:kgz_max), UNIT=6, format="(20F9.3)")

c**** Initialize acc-array names, units, idacc-indices
      call def_acc

C**** Ensure that diagnostics are reset at the beginning of the run
      IF (Itime.le.ItimeI .and. ISTART.gt.0) THEN
        call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour
     *       ,amon)
        CALL reset_DIAG(0)
C**** Initiallise ice freeze diagnostics at beginning of run
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            TSFREZ(I,J,TF_DAY1)=365.
            TSFREZ(I,J,TF_LAST)=365.
            IF (FOCEAN(I,J)+FLAKE(I,J).gt.0) then
              IF (RSI(I,J).gt.0) then
                TSFREZ(I,J,TF_LKON) = JDAY-1
                TSFREZ(I,J,TF_LKOFF) = JDAY
              ELSE
                TSFREZ(I,J,TF_LKON) = JDAY
                TSFREZ(I,J,TF_LKOFF) = undef
              END IF
            ELSE
              TSFREZ(I,J,TF_LKON) = undef
              TSFREZ(I,J,TF_LKOFF) = undef
            END IF
          END DO
        END DO
        CALL daily_DIAG
      END IF

      RETURN
      END SUBROUTINE init_DIAG


      SUBROUTINE reset_DIAG(isum)
!@sum  reset_DIAG resets/initializes diagnostics
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : Itime,iyear1,nday,kradia,
     *     Itime0,jhour0,jdate0,jmon0,amon0,jyear0,idacc,u
      USE DIAG_COM
      USE PARAM
#ifdef TRACERS_ON
      USE TRDIAG_COM, only: taijln, TAIJLN_loc
      USE TRDIAG_COM, only: taijn,  TAIJN_loc
      USE TRDIAG_COM, only: taijs,  TAIJS_loc
      USE TRDIAG_COM, only: TAJLN,  TAJLN_loc
      USE TRDIAG_COM, only: TAJLS,  TAJLS_loc
      USE TRDIAG_COM, only: TCONSRV,TCONSRV_loc
#endif
      USE DOMAIN_DECOMP, only: grid, CHECKSUM
      IMPLICIT NONE
      INTEGER :: isum !@var isum if =1 preparation to add up acc-files
      INTEGER jd0

      IDACC(1:12)=0
      if (kradia.gt.0) then
        AFLX_ST = 0.
        if (isum.eq.1) return
        go to 100
      end if
      AJ_loc=0    ; AREG=0

      APJ_loc=0   ; AJL_loc=0  ; ASJL_loc=0   ; AIJ_loc=0
      AIL=0   ; ENERGY=0 ; CONSRV_loc=0
      SPECA=0 ; ATPE=0 ; WAVE=0
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      AJK_loc=0   ; AIJK_loc=0
#else
      AJK_loc=0   ; AIJK_loc=0 ; HDIURN=0
#endif
      AISCCP=0
#ifdef TRACERS_ON
       TAJLN=0. ; TAJLS=0. ; TCONSRV=0.
       TAJLN_loc=0. ; TAJLS_loc=0. ; TCONSRV_loc=0.
       TAIJLN=0.    ; TAIJN=0.     ; TAIJS=0.
      TAIJLN_loc=0. ; TAIJN_loc=0. ; TAIJS_loc=0.
#endif
      call reset_ODIAG(isum)  ! ocean diags if required
      call reset_icdiag       ! ice dynamic diags if required

      if (isum.eq.1) then ! prepare to add up acc-files
         AJ=0 ; APJ=0 ; AJL=0 ; ASJL=0 ; AIJ=0 ; AJK=0 ; AIJK=0
         return
      end if

      AIJ_loc(:,:,IJ_TMNMX)=1000. ; IDACC(12)=1

      CALL EPFLXI (U)  ! strat

  100 Itime0=Itime
      call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jd0,
     *     Jdate0,Jhour0,amon0)

      RETURN
      END SUBROUTINE reset_DIAG


      SUBROUTINE daily_DIAG
!@sum  daily_DIAG resets diagnostics at beginning of each day
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : undef
      USE MODEL_COM, only : im,jm,jday,fearth,focean
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE DIAG_COM, only : aij=>aij_loc
     *     ,ij_lkon,ij_lkoff,ij_lkice,tsfrez=>tsfrez_loc,tdiurn
     *     ,tf_lkon,tf_lkoff,tf_day1,tf_last
      USE DOMAIN_DECOMP, only : GRID,GET
      IMPLICIT NONE
      INTEGER I,J
      INTEGER :: J_0, J_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      IF (JDAY.EQ.32) THEN
         DO J=MAX(J_0,1+JM/2),MIN(J_1,JM)
            DO I=1,IM
               TSFREZ(I,J,TF_DAY1)=JDAY
            END DO
         END DO
         DO J=MAX(J_0,1),MIN(J_1,JM/2)
            DO I=1,IM
               TSFREZ(I,J,TF_LAST)=JDAY
            END DO
         END DO
      ELSEIF (JDAY.EQ.213) THEN
         DO J=MAX(J_0,1),MIN(J_1,JM/2)
            DO I=1,IM
              TSFREZ(I,J,TF_DAY1)=JDAY
            END DO
         END DO
      END IF

C**** set and initiallise freezing diagnostics
C**** Note that TSFREZ saves the last day of no-ice and some-ice.
C**** The AIJ diagnostics are set once a year (zero otherwise)
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF (J.le.JM/2) THEN
C**** initiallise/save South. Hemi. on Feb 28
            IF (JDAY.eq.59 .and. TSFREZ(I,J,TF_LKOFF).ne.undef) THEN
              AIJ(I,J,IJ_LKICE)=1.
              AIJ(I,J,IJ_LKON) =MOD(NINT(TSFREZ(I,J,TF_LKON)) +307,365)
              AIJ(I,J,IJ_LKOFF)=MOD(NINT(TSFREZ(I,J,TF_LKOFF))+306,365)
     *             +1
              IF (RSI(I,J).gt.0) THEN
                TSFREZ(I,J,TF_LKON) = JDAY-1
              ELSE
                TSFREZ(I,J,TF_LKOFF) = undef
              END IF
            END IF
          ELSE
C**** initiallise/save North. Hemi. on Aug 31
C**** Note that for continuity across the new year, the julian days
C**** are counted from Sep 1 (NH only).
            IF (JDAY.eq.243 .and. TSFREZ(I,J,TF_LKOFF).ne.undef) THEN
              AIJ(I,J,IJ_LKICE)=1.
              AIJ(I,J,IJ_LKON) =MOD(NINT(TSFREZ(I,J,TF_LKON)) +123,365)
              AIJ(I,J,IJ_LKOFF)=MOD(NINT(TSFREZ(I,J,TF_LKOFF))+122,365)
     *             +1
              IF (RSI(I,J).gt.0) THEN
                TSFREZ(I,J,TF_LKON) = JDAY-1
              ELSE
                TSFREZ(I,J,TF_LKOFF) = undef
              END IF
            END IF
          END IF
C**** set ice on/off days
          IF (FOCEAN(I,J)+FLAKE(I,J).gt.0) THEN
            IF (RSI(I,J).eq.0.and.TSFREZ(I,J,TF_LKOFF).eq.undef)
     *           TSFREZ(I,J,TF_LKON)=JDAY
            IF (RSI(I,J).gt.0) TSFREZ(I,J,TF_LKOFF)=JDAY
          END IF
        END DO
      END DO

C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF EACH DAY
      DO J=J_0,J_1
         DO I=1,IM
            TDIURN(I,J,1)= 1000.
            TDIURN(I,J,2)=-1000.
            TDIURN(I,J,3)= 1000.
            TDIURN(I,J,4)=-1000.
            TDIURN(I,J,5)=    0.
            TDIURN(I,J,6)=-1000.
            TDIURN(I,J,7)=-1000.
            TDIURN(I,J,8)=-1000.
            TDIURN(I,J,9)= 1000.
            IF (FEARTH(I,J).LE.0.) THEN
               TSFREZ(I,J,TF_DAY1)=365.
               TSFREZ(I,J,TF_LAST)=365.
            END IF
         END DO
      END DO

      END SUBROUTINE daily_DIAG


      SUBROUTINE SET_CON(QCON,CONPT,NAME_CON,INST_UNIT,SUM_UNIT,INST_SC
     *     ,CHNG_SC,ICON)
!@sum  SET_CON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday
      USE MODEL_COM, only : dtsrc,nfiltr
      USE DIAG_COM, only : kcon,nquant,npts,title_con,scale_con,nsum_con
     *     ,nofm,ia_con,kcmx,ia_d5d,ia_d5s,ia_filt,ia_12hr,name_consrv
     *     ,lname_consrv,units_consrv
      USE DOMAIN_DECOMP, only : WRITE_PARALLEL
      IMPLICIT NONE
!@var QCON logical variable sets where conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(NPTS) :: QCON
!@var CONPT names for points where conservation diags are saved
      CHARACTER*10, INTENT(IN),DIMENSION(NPTS) :: CONPT
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*16, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*16, INTENT(IN) :: SUM_UNIT
!@var ICON index for the conserved quantity
      INTEGER, INTENT(OUT) :: ICON
!@var out_line local variable to hold mixed-type output for parallel I/O
      character(len=300) :: out_line

      INTEGER NI,NM,NS,N,k
      INTEGER, SAVE :: NQ = 2   ! first 2 special cases AM + KE

      NQ=NQ+1
      IF (NQ.gt.NQUANT) THEN
C****        WRITE(6,*) "Number of conserved quantities larger than NQUANT"
C****     *       ,NQUANT,NQ
        WRITE(out_line,*)
     *       "Number of conserved quantities larger than NQUANT"
     *       ,NQUANT,NQ
        CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        call stop_model("Change NQUANT in diagnostic common block",255)
      END IF
C**** remove spaces in NAME_CON for netcdf names
      sname=TRIM(NAME_CON)
      do k=1,len_trim(NAME_CON)
        if (sname(k:k).eq." ") sname(k:k)="_"
      end do
C****
      NI=KCMX+1
      NOFM(1,NQ) = NI
      TITLE_CON(NI) = "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_CON(NI) = INST_SC
      name_consrv(NI) ="inst_"//sname
      lname_consrv(NI) = "INSTANT "//TRIM(NAME_CON)
      units_consrv(NI) = INST_UNIT
      IA_CON(NI) = 12
      NM=NI
      DO N=1,NPTS
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFM(N+1,NQ) = NM
          TITLE_CON(NM) = " CHANGE OF "//TRIM(NAME_CON)//" BY "//
     *         CONPT(N)
          name_consrv(NM) ="chg_"//trim(sname)//"_"//TRIM(CONPT(N)(1:3))
          lname_consrv(NM) = TITLE_CON(NM)
          units_consrv(NM) = SUM_UNIT
          SELECT CASE (N)
          CASE (1)
            SCALE_CON(NM) = CHNG_SC/DTSRC
            IA_CON(NM) = ia_d5d
          CASE (2,3,4,5,6,8,10,11)
            SCALE_CON(NM) = CHNG_SC/DTSRC
            IA_CON(NM) = ia_d5s
          CASE (7)
            SCALE_CON(NM) = CHNG_SC/(NFILTR*DTSRC)
            IA_CON(NM) = ia_filt
          CASE (9)
            SCALE_CON(NM) = CHNG_SC*2./SDAY
            IA_CON(NM) = ia_12hr
          END SELECT
        ELSE
          NOFM(N+1,NQ) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KCON) THEN
C****        WRITE(6,*) "KCON not large enough for extra conserv diags",
C****     *       KCON,NI,NM,NQ,NS,NAME_CON
        WRITE(out_line,*)
     *      "KCON not large enough for extra conserv diags",
     *       KCON,NI,NM,NQ,NS,NAME_CON
        CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        call stop_model("Change KCON in diagnostic common block",255)
      END IF
      TITLE_CON(NS) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_consrv(NS) ="sum_chg_"//trim(sname)
      lname_consrv(NS) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_consrv(NS) = SUM_UNIT
      SCALE_CON(NS) = 1.
      IA_CON(NS) = 12
      NSUM_CON(NI) = -1
      NSUM_CON(NI+1:NS-1) = NS
      NSUM_CON(NS) = 0
      KCMX=NS
      ICON=NQ
      RETURN
C****
      END SUBROUTINE set_con

      SUBROUTINE UPDTYPE
!@sum UPDTYPE updates FTYPE array to ensure correct budget diagnostics
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,focean,fearth,flice,itocean
     *     ,itoice,itlandi,itearth,itlake,itlkice,ftype
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE DOMAIN_DECOMP, only : GRID,GET
      IMPLICIT NONE
      INTEGER I,J
      INTEGER :: J_0,J_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          FTYPE(ITOICE ,I,J)=FOCEAN(I,J)*RSI(I,J)
          FTYPE(ITOCEAN,I,J)=FOCEAN(I,J)-FTYPE(ITOICE,I,J)
          FTYPE(ITLKICE,I,J)=FLAKE(I,J)*RSI(I,J)
          FTYPE(ITLAKE ,I,J)=FLAKE(I,J)-FTYPE(ITLKICE,I,J)
C**** set land components of FTYPE array. Summation is necessary for
C**** cases where Earth and Land Ice are lumped together
          FTYPE(ITLANDI,I,J)=0.
          FTYPE(ITEARTH,I,J)=FEARTH(I,J)
          FTYPE(ITLANDI,I,J)=FTYPE(ITLANDI,I,J)+FLICE(I,J)
        END DO
      END DO
      RETURN
C****
      END SUBROUTINE UPDTYPE
