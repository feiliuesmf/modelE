#include "rundeck_opts.h"
#ifdef SKIP_TRACER_DIAGS
#undef TRACERS_SPECIAL_O18
#endif

!@sum  DIAG ModelE diagnostic calculations
!@auth G. Schmidt/J. Lerner/R. Ruedy/M. Kelley
!@ver  1.0
C**** AJ(J,N)  (ZONAL SUM OVER LONGITUDE AND TIME)
C****   See j_defs for contents
C****                                                             IDACC
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
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
C**** Variables passed from DIAGA to DIAGB
!@var W,TX vertical velocity and in-situ temperature calculations
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: W
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TX

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
      USE DOMAIN_DECOMP_ATM, only : GET
      USE DOMAIN_DECOMP_ATM, only : DIST_GRID
      USE MODEL_COM, only : lm
      USE DIAG_LOC, only  : W,TX
      IMPLICIT NONE
      LOGICAL, SAVE :: init=.false.
      INTEGER :: I_0H,I_1H,J_0H,J_1H
      INTEGER :: IER
      TYPE(DIST_GRID) :: grid

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = GRID%I_STRT_HALO
      I_1H = GRID%I_STOP_HALO

      ALLOCATE( W(I_0H:I_1H, J_0H:J_1H, LM),
     &         TX(I_0H:I_1H, J_0H:J_1H, LM),
     &     STAT = IER)

      !hack hack hack!
      TX(:,:,:) = 0.d0

      RETURN
      END SUBROUTINE ALLOC_DIAG_LOC

      SUBROUTINE DIAGA
!@sum  DIAGA accumulate various diagnostics during dynamics
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,lhe,sha,bygrav,tf
     *     ,rvap,gamd,teeny,undef,radius,omega,kg2mb
      USE MODEL_COM, only : im,jm,lm,ls1,idacc,ptop
     *     ,pmtop,psfmpt,mdyn,mdiag,sig,sige,dsig,zatmo,WM,ntype,ftype
     *     ,u,v,t,p,q,lm_req,req_fac_m,pmidl00
      USE GEOM, only : sinlat2d,coslat2d,axyp,imaxj,ddy_ci,ddy_cj,
     &     lon2d_dg,byaxyp
      USE RAD_COM, only : rqt
      USE DIAG_COM, only : ia_dga,jreg,
     *     aijl=>aijl_loc
     *     ,aij=>aij_loc,ij_dtdp,ij_phi1k,ij_pres,ij_slpq,ij_presq
     *     ,ij_slp,ij_t850,ij_t500,ij_t300,ij_q850,ij_q500
     *     ,ij_RH1,ij_RH850,ij_RH500,ij_RH300,ij_qm,ij_q300,ij_ujet
     *     ,ij_vjet,j_tx1,j_tx,j_qp,j_dtdjt,j_dtdjs,j_dtdgtr,j_dtsgst
     &     ,ijl_dp,ijk_dp,ijl_u,ijl_v,ijl_w,ijk_tx,ijk_q,ijk_rh
     *     ,j_rictr,j_rostr,j_ltro,j_ricst,j_rosst,j_lstr,j_gamm,j_gam
     *     ,j_gamc,lstr,kgz_max,pmb,ght,ple
     *     ,jl_dtdyn,jl_dpa
     *     ,jl_epacwt,jl_uepac,jl_vepac,jl_wepac
     *     ,jl_wpacwt,jl_uwpac,jl_vwpac,jl_wwpac
     *     ,jk_dpwt,jk_tx,jk_hght,jk_q,jk_rh,jk_cldh2o
     *     ,ij_p850,z_inst,rh_inst,t_inst,plm,ij_p1000,ij_p925,ij_p700
     *     ,ij_p600,ij_p500
#ifdef HTAP_LIKE_DIAGS
     *     ,ij_templ,ij_gridh,ij_husl
#endif
      USE DYNAMICS, only : pk,phi,pmid,pdsig,plij, SD,pedn,am
     &     ,ua=>ualij,va=>valij
      USE PBLCOM, only : tsavg
      USE DIAG_LOC, only : w,tx,jet
      USE DOMAIN_DECOMP_ATM, only : GET, HALO_UPDATE, GRID
      USE GETTIME_MOD
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     TX_TROP,TX_STRAT
      REAL*8, DIMENSION(LM_REQ) :: TRI

      REAL*8, PARAMETER :: ONE=1.,P1000=1000.
      INTEGER :: I,IM1,J,K,L,JR,
     &     IP1,LM1,LP1,LR,MBEGIN,IT
      REAL*8 THBAR ! external
      REAL*8 ::
     &     BBYGV,BYSDSG,DLNP,DLNP01,DLNP12,DLNP23,DBYSD,
     &     DXYPJ,
     *     ESEPS,GAMC,GAMM,GAMX,
     &     PDN,PE,PHI_REQ,PIJ,
     *     PKE,PL,PRT,W2MAX,RICHN,
     *     ROSSN,ROSSL,BYFCOR,BYBETA,BYBETAFAC,NH,SS,THETA,
     *     TZL,X,TIJK,QIJK,DTXDY
      LOGICAL qpress,qabove
      INTEGER nT,nQ,nRH
      REAL*8, PARAMETER :: EPSLON=1.

      REAL*8 QSAT, SLP, PS, ZS, QLH

      real*8, dimension(lm+1) :: pecp,pedge
      real*8, dimension(lm) :: dpwt,txdp,phidp,qdp,rhdp,wmdp,rh
      integer, parameter :: lmxmax=2*lm
      integer :: lmx
      real*8, dimension(lmxmax) :: dpx
      integer, dimension(lmxmax) :: lmod,lcp

      INTEGER :: J_0, J_1, J_0S, J_1S, I_0,I_1,I_0H,I_1H
     &     ,IM1S,IP1S,IP1E
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(MBEGIN)

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      I_0H = GRID%I_STRT_HALO
      I_1H = GRID%I_STOP_HALO

c
c get winds on the atmospheric primary grid
c
      call recalc_agrid_uv

      IDACC(ia_dga)=IDACC(ia_dga)+1

      BYSDSG=1./(1.-SIGE(LM+1))
      DLNP01=LOG(pmidl00(lm)/PLM(LM+1))
      DLNP12=LOG(REQ_FAC_M(1)/REQ_FAC_M(2))  ! LOG(.75/.35)
      DLNP23=LOG(REQ_FAC_M(2)/REQ_FAC_M(3))  ! LOG(.35/.1)
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
        DO I=I_0,I_1
          TX(I,J,L)=T(I,J,L)*PK(L,I,J)
        END DO
      END DO
      END DO


C****
C**** J LOOPS FOR ALL PRIMARY GRID ROWS
C****
      DO J=J_0,J_1
C**** NUMBERS ACCUMULATED FOR A SINGLE LEVEL
        DO I=I_0,IMAXJ(J)
          DXYPJ=AXYP(I,J)
          JR=JREG(I,J)
          DO IT=1,NTYPE
            CALL INC_AJ(I,J,IT,J_TX1,(TX(I,J,1)-TF)*FTYPE(IT,I,J))
          END DO
          CALL INC_AREG(I,J,JR,J_TX1,(TX(I,J,1)-TF))
          PS=P(I,J)+PTOP
          ZS=BYGRAV*ZATMO(I,J)
          AIJ(I,J,IJ_PRES)=AIJ(I,J,IJ_PRES)+ PS
          AIJ(I,J,IJ_SLP)=AIJ(I,J,IJ_SLP)+SLP(PS,TSAVG(I,J),ZS)-P1000
C**** calculate pressure diags including water 
          PS=PS+SUM((Q(I,J,:)+WM(I,J,:))*AM(:,I,J))*kg2mb
          AIJ(I,J,IJ_PRESQ)=AIJ(I,J,IJ_PRESQ)+ PS
          AIJ(I,J,IJ_SLPQ)=AIJ(I,J,IJ_SLPQ)+SLP(PS,TSAVG(I,J),ZS)-P1000

          AIJ(I,J,IJ_RH1)=AIJ(I,J,IJ_RH1)+Q(I,J,1)/QSAT(TX(I,J,1),LHE,
     *        PMID(1,I,J))
        END DO
c        APJ(J,1)=APJ(J,1)+PI(J)
C**** CALCULATE GEOPOTENTIAL HEIGHTS AT SPECIFIC MILLIBAR LEVELS
        DO I=I_0,IMAXJ(J)
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

C**** ACCUMULATION OF TEMP., POTENTIAL TEMP., Q, AND RH
      DO J=J_0,J_1
        DO L=1,LM
          DBYSD=DSIG(L)*BYSDSG
          DO I=I_0,IMAXJ(J)
            DXYPJ=AXYP(I,J)
            JR=JREG(I,J)
            PIJ=PLIJ(L,I,J)
            aijl(i,j,l,ijl_dp) = aijl(i,j,l,ijl_dp) + pdsig(l,i,j)
            call inc_ajl(i,j,l,jl_dpa,pdsig(l,i,j))
c ajl(jl_dtdyn) was incremented by -t(i,j,l) before dynamics
            call inc_ajl(i,j,l,jl_dtdyn,tx(i,j,l)*pdsig(l,i,j))
            AIJ(I,J,IJ_QM)=AIJ(I,J,IJ_QM)+Q(I,J,L)*AM(L,I,J)
#ifdef HTAP_LIKE_DIAGS
            AIJ(I,J,IJ_TEMPL(L))=AIJ(I,J,IJ_TEMPL(L))+TX(I,J,L)
            AIJ(I,J,IJ_HUSL(L))=AIJ(I,J,IJ_HUSL(L))+Q(I,J,L)
            AIJ(I,J,IJ_GRIDH(L))=AIJ(I,J,IJ_GRIDH(L))+
     &      rgas/grav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
#endif
            DO IT=1,NTYPE
              CALL INC_AJ(I,J,IT,J_TX,(TX(I,J,L)-TF)*FTYPE(IT,I,J)*
     *             DBYSD)
              CALL INC_AJ(I,J,IT,J_QP,(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)
     *             *FTYPE(IT,I,J)) 
            END DO
            CALL INC_AREG(I,J,JR,J_QP,(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)) 
            CALL INC_AREG(I,J,JR,J_TX,(TX(I,J,L)-TF)*DBYSD)
          END DO
        END DO
      END DO


C****
C**** STATIC STABILITIES: TROPOSPHERIC AND STRATOSPHERIC
C****
      DO J=J_0,J_1
C**** OLD TROPOSPHERIC STATIC STABILITY
        DO I=I_0,IMAXJ(J)
          DXYPJ=AXYP(I,J)
          JR=JREG(I,J)
          SS=(T(I,J,LS1-1)-T(I,J,1))/(PHI(I,J,LS1-1)-PHI(I,J,1)+teeny)
          DO IT=1,NTYPE
            CALL INC_AJ(I,J,IT,J_DTDGTR,SS*FTYPE(IT,I,J))
          END DO
          CALL INC_AREG(I,J,JR,J_DTDGTR,SS)
          AIJ(I,J,IJ_DTDP)=AIJ(I,J,IJ_DTDP)+SS
        END DO
C**** OLD STRATOSPHERIC STATIC STABILITY (USE LSTR as approx 10mb)
        DO I=I_0,IMAXJ(J)
          JR=JREG(I,J)
          SS=(T(I,J,LSTR)-T(I,J,LS1-1))/((PHI(I,J,LSTR)-PHI(I,J,LS1-1))
     *         +teeny)
          DO IT=1,NTYPE
            CALL INC_AJ(I,J,IT,J_DTSGST,SS*FTYPE(IT,I,J))
          END DO
          CALL INC_AREG(I,J,JR,J_DTSGST,SS)
        END DO

C****
C**** NUMBERS ACCUMULATED FOR THE RADIATION EQUILIBRIUM LAYERS
C****
        DO I=I_0,IMAXJ(J)
          DO LR=1,LM_REQ
            TRI(LR)=RQT(LR,I,J)
            call inc_asjl(i,j,lr,1,RQT(LR,I,J)-TF)
          END DO
          PHI_REQ=PHI(I,J,LM)
          PHI_REQ=PHI_REQ+RGAS*.5*(TX(I,J,LM)+TRI(1))*DLNP01
          call inc_asjl(i,j,1,2,PHI_REQ)
          PHI_REQ=PHI_REQ+RGAS*.5*(TRI(1)+TRI(2))*DLNP12
          call inc_asjl(i,j,2,2,PHI_REQ)
          PHI_REQ=PHI_REQ+RGAS*.5*(TRI(2)+TRI(3))*DLNP23
          call inc_asjl(i,j,3,2,PHI_REQ)
        END DO
      END DO

C****
C**** RICHARDSON NUMBER , ROSSBY NUMBER , RADIUS OF DEFORMATION
C****

c
c Spatial mean Rossby number is computed using the spatial mean of
c column-maximum wind speeds.
c
c Spatial mean Rossby radius is computed using the spatial mean of
c min(N*H/f,sqrt(N*H/beta)) where H is the either the depth of
c the troposphere or the lower stratosphere. (N*H)**2 is taken as
c log(theta2/theta1)*grav*H where theta2,theta1 are the potential
c temperatures at the top/bottom of the troposphere or lower stratosphere.
c f is the coriolis parameter and beta its latitudinal derivative.
c
c The Richardson number is not being computed.  It is set to 99 or 999
c to flag that it is missing.  For the large vertical ranges of these
c calculations, it should be replaced by something like the ratio of
c available mechanical energy to the work required to homogenize
c potential temperature.
c

      X=RGAS*LHE*LHE/(SHA*RVAP)
      bybetafac = radius/(2.*omega)

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

        byfcor = 1d0/(2.*omega*abs(sinlat2d(i,j))+teeny)
        bybeta = bybetafac/(coslat2d(i,j)+teeny)
c troposphere
        w2max = maxval(ua(1:ls1-1,i,j)**2+va(1:ls1-1,i,j)**2)+teeny
        rossn = sqrt(w2max)*byfcor
        nh = sqrt((phi(i,j,ls1-1)-phi(i,j,1))*
     &       log(t(i,j,ls1-1)/t(i,j,1)))
        rossl = min(nh*byfcor,sqrt(nh*bybeta))
        richn = 99d0 ! for now

c lapse rates
        gamx = (tx(i,j,1)-tx(i,j,ls1-1))/
     &       (phi(i,j,ls1-1)-phi(i,j,1))
        gamm=0.
        do l=1,ls1-1
          tzl=tx(i,j,l)
          prt=(sig(l)*p(i,j)+ptop)*rgas*tzl
          eseps=qsat(tzl,lhe,one)
          gamm=gamm+dsig(l)*(prt+lhe*eseps)/(prt+x*eseps/tzl)
        end do

        do it=1,ntype
          call inc_aj(i,j,it,j_rostr,rossn*ftype(it,i,j))
          call inc_aj(i,j,it,j_ltro,rossl*ftype(it,i,j))
          call inc_aj(i,j,it,j_rictr,richn*ftype(it,i,j))

          call inc_aj(i,j,it,j_gam ,gamx*ftype(it,i,j))
          call inc_aj(i,j,it,j_gamm,gamm*ftype(it,i,j))

        end do

c stratosphere
        w2max = maxval(ua(ls1:lstr,i,j)**2+va(ls1:lstr,i,j)**2)
        rossn = sqrt(w2max)*byfcor
        nh = sqrt((phi(i,j,lstr)-phi(i,j,ls1))*
     &       log(t(i,j,lstr)/t(i,j,ls1)))
        rossl = min(nh*byfcor,sqrt(nh*bybeta))
        richn = 999d0 ! for now
        do it=1,ntype
          call inc_aj(i,j,it,j_rosst,rossn*ftype(it,i,j))
          call inc_aj(i,j,it,j_lstr,rossl*ftype(it,i,j))
          call inc_aj(i,j,it,j_ricst,richn*ftype(it,i,j))
        end do

      ENDDO
      ENDDO

C****
C**** NORTHWARD GRADIENT OF TEMPERATURE: TROPOSPHERIC AND STRATOSPHERIC
C**** GAMC, THE DYNAMICALLY DETERMINED LAPSE RATE IN THE EXTRATROPICS
C****
      do j=j_0,j_1
      do i=i_0,i_1
        tx_trop(i,j) = 0.
        tx_strat(i,j) = 0.
      enddo
      enddo
      do l=1,ls1-1
      do j=j_0,j_1
      do i=i_0,i_1
        tx_trop(i,j) = tx_trop(i,j) + tx(i,j,l)*dsig(l)
      enddo
      enddo
      enddo
      do l=ls1,lstr
      do j=j_0,j_1
      do i=i_0,i_1
        tx_strat(i,j) = tx_strat(i,j) + tx(i,j,l)*dsig(l)
      enddo
      enddo
      enddo
      do j=j_0,j_1
      do i=i_0,i_1
        tx_trop(i,j) = tx_trop(i,j)/(SIGE(1)-SIGE(LS1))
        tx_strat(i,j) = tx_strat(i,j)/(SIGE(LS1)-SIGE(LSTR+1)+1d-12)
      enddo
      enddo

      call halo_update(grid, tx_trop)
      call halo_update(grid, tx_strat)

      if(i_0h.lt.i_0) then      ! halo cells exist in i direction
        im1s=i_0-1; ip1s=i_0+1; ip1e=i_1+1
      else                      ! periodic
        im1s=im-1; ip1s=1; ip1e=im
      endif
      do j=j_0s,j_1s
        im1=im1s
        i=im1s+1
        do ip1=ip1s,ip1e
          dtxdy = (tx_trop(ip1,j)-tx_trop(im1,j))*ddy_ci(i,j)
     &          + (tx_trop(i,j+1)-tx_trop(i,j-1))*ddy_cj(i,j)
          gamc = gamd+grav*radius*dtxdy*sinlat2d(i,j)/
     &         (rgas*tx_trop(i,j)*(coslat2d(i,j)+.001))
          do it=1,ntype
            call inc_aj(i,j,it,j_dtdjt,dtxdy*ftype(it,i,j))
            call inc_aj(i,j,it,j_gamc,gamc*ftype(it,i,j))
          enddo
          dtxdy = (tx_strat(ip1,j)-tx_strat(im1,j))*ddy_ci(i,j)
     &          + (tx_strat(i,j+1)-tx_strat(i,j-1))*ddy_cj(i,j)
          do it=1,ntype
            call inc_aj(i,j,it,j_dtdjs,dtxdy*ftype(it,i,j))
          enddo
          im1=i
          i=ip1
        enddo
      enddo

      DO J=J_0S,J_1S
      DO I=I_0,I_1
        AIJ(I,J,IJ_UJET)=AIJ(I,J,IJ_UJET)+UA(JET,I,J)
        AIJ(I,J,IJ_VJET)=AIJ(I,J,IJ_VJET)+VA(JET,I,J)
      ENDDO
      ENDDO

C****
C**** CONVERT VERTICAL WINDS TO UNITS PROPORTIONAL TO M/S
C****
      DO L=1,LM-1
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        PIJ=PLIJ(L,I,J)
        PE=SIGE(L+1)*PIJ+PTOP
        PKE=PE**KAPA
        THETA=THBAR(T(I,J,L+1),T(I,J,L))
        W(I,J,L)=SD(I,J,L)*THETA*PKE/PE
      END DO
       if(have_south_pole .and. J==1) W(2:IM,J,L) = W(1,J,L)
       if(have_north_pole .and. J==JM) W(2:IM,J,L) = W(1,J,L)
      END DO
      END DO

c
c accumulate AIJL: U,V,W on model layers
c

      DO L=1,LM
      DO J=J_0S,J_1S
      DO I=I_0,I_1
        AIJL(I,J,L,IJL_U) = AIJL(I,J,L,IJL_U) + UA(L,I,J)
        AIJL(I,J,L,IJL_V) = AIJL(I,J,L,IJL_V) + VA(L,I,J)
      ENDDO
      ENDDO
      ENDDO ! L
      DO L=1,LM-1
      DO J=J_0,J_1
      DO I=I_0,I_1
        AIJL(I,J,L,IJL_W) = AIJL(I,J,L,IJL_W) + W(I,J,L)
      ENDDO
      ENDDO
      ENDDO

C****
C**** CERTAIN HORIZONTAL WIND AVERAGES
C****
      do j=j_0,j_1
      do i=i_0,imaxj(j)

        if(lon2d_dg(i,j).ge.-135. .and. lon2d_dg(i,j).le.-110.) then
c east pacific
          do l=1,lm
            call inc_ajl(i,j,l,jl_epacwt,pdsig(l,i,j))
            call inc_ajl(i,j,l,jl_uepac, pdsig(l,i,j)*ua(l,i,j))
            call inc_ajl(i,j,l,jl_vepac, pdsig(l,i,j)*va(l,i,j))
          enddo
          do l=1,lm-1 ! using pdsig for vertical velocity weight
            call inc_ajl(i,j,l,jl_wepac,
     &           pdsig(l,i,j)*w(i,j,l)*byaxyp(i,j))
          enddo
        elseif(lon2d_dg(i,j).ge.150.) then
c west pacific
          do l=1,lm
            call inc_ajl(i,j,l,jl_wpacwt,pdsig(l,i,j))
            call inc_ajl(i,j,l,jl_uwpac, pdsig(l,i,j)*ua(l,i,j))
            call inc_ajl(i,j,l,jl_vwpac, pdsig(l,i,j)*va(l,i,j))
          enddo
          do l=1,lm-1 ! using pdsig for vertical velocity weight
            call inc_ajl(i,j,l,jl_wwpac,
     &           pdsig(l,i,j)*w(i,j,l)*byaxyp(i,j))
          enddo
        endif

      enddo
      enddo

c
c constant-pressure diagnostics
c
      pecp(2:lm+1) = ple(1:lm)
      pecp(1) = 1d30 ! ensure that all column mass is included
      qlh=lhe
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        pedge(:) = pedn(:,i,j)
        call get_dx_intervals(pedge,lm,pecp,lm,dpx,lmod,lcp,lmx,lmxmax)
        do l=1,lm
          dpwt(l) = 0d0
          txdp(l) = 0d0
          phidp(l) = 0d0
          qdp(l) = 0d0
          rhdp(l) = 0d0
          wmdp(l) = 0d0
          rh(l) = q(i,j,l)/min(1.,QSAT(TX(I,J,L),QLH,pmid(l,i,j)))
        enddo
        do l=1,lmx
          dpwt(lcp(l))  = dpwt(lcp(l))  + dpx(l)
          txdp(lcp(l))  = txdp(lcp(l))  + dpx(l)*(tx(i,j,lmod(l))-tf)
          phidp(lcp(l)) = phidp(lcp(l)) + dpx(l)*phi(i,j,lmod(l))
          qdp(lcp(l))   = qdp(lcp(l))   + dpx(l)*q(i,j,lmod(l))
          rhdp(lcp(l))  = rhdp(lcp(l))  + dpx(l)*rh(lmod(l))
          wmdp(lcp(l))  = wmdp(lcp(l))  + dpx(l)*wm(i,j,lmod(l))
        enddo
        do l=1,lm
          aijl(i,j,l,ijk_dp) = aijl(i,j,l,ijk_dp) + dpwt(l)
          aijl(i,j,l,ijk_tx) = aijl(i,j,l,ijk_tx) + txdp(l)
          aijl(i,j,l,ijk_q)  = aijl(i,j,l,ijk_q)  + qdp(l)
          aijl(i,j,l,ijk_rh) = aijl(i,j,l,ijk_rh) + rhdp(l)
          call inc_ajl(i,j,l,jk_dpwt,  dpwt(l))
          call inc_ajl(i,j,l,jk_tx,    txdp(l))
          call inc_ajl(i,j,l,jk_hght,  phidp(l))
          call inc_ajl(i,j,l,jk_q,     qdp(l))
          call inc_ajl(i,j,l,jk_rh,    rhdp(l))
          call inc_ajl(i,j,l,jk_cldh2o,wmdp(l))
        enddo
      enddo
      enddo

C**** ACCUMULATE TIME USED IN DIAGA
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN

      ENTRY DIAGA0
c increment ajl(jl_dtdyn) by -t before dynamics.
c ajl(jl_dtdyn) will be incremented by +t after the dynamics, giving
c the tendency.

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

      DO L=1,LM
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        call inc_ajl(i,j,l,jl_dtdyn,-t(i,j,l)*pk(l,i,j)*pdsig(l,i,j))
      END DO
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE DIAGA

      subroutine get_dx_intervals(
     &     xesrc,nsrc,xedst,ndst,dx,lsrc,ldst,nxchng,nmax)
c
c This routine returns information needed for conservative
c remapping in 1 dimension.
c Given two lists of points xesrc(1:nsrc+1) and xedst(1:ndst+1)
c along the "x" axis, merge the two lists and calculate the
c distance increments dx(1:nxchng) separating the points
c in the merged list.  The i_th increment lies in the intervals
c xesrc(lsrc(i):lsrc(i)+1) and xedst(ldst(i):ldst(i)+1).
c Points in xedst lying outside the range xesrc(1:nsrc+1) are
c not included.
c
      implicit none
      integer :: nsrc,ndst,nxchng,nmax
      real*8, dimension(nsrc+1) :: xesrc
      real*8, dimension(ndst+1) :: xedst
      real*8, dimension(nmax) :: dx
      integer, dimension(nmax) :: lsrc,ldst
      integer :: ls,ld
      real*8 :: xlast,s
      s = sign(1d0,xesrc(2)-xesrc(1))
      ld = 1
      xlast = max(s*xesrc(1),s*xedst(1))
      do while(s*xedst(ld).le.xlast)
        ld = ld + 1
      enddo
      ls = 2
      do nxchng=1,nmax
        lsrc(nxchng) = ls-1
        ldst(nxchng) = ld-1
        if(s*xedst(ld).le.s*xesrc(ls)) then
          if(xedst(ld).eq.xesrc(ls)) ls = ls + 1
          dx(nxchng) = s*xedst(ld)-xlast
          xlast = s*xedst(ld)
          ld = ld + 1
        else
          dx(nxchng) = s*xesrc(ls)-xlast
          xlast = s*xesrc(ls)
          ls = ls + 1
        endif
        if(ld.gt.ndst+1) exit
        if(ls.gt.nsrc+1) exit
      enddo
      return
      end subroutine get_dx_intervals

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
      !USE SOIL_DRV, only: conserv_WTG,conserv_HTG
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
     *     ,conserv_MLI,conserv_HLI,conserv_WTG,conserv_HTG
      INTEGER MNOW
      INTEGER NT

#ifdef SCM
c     write(0,*) 'SCM - no diags   DIAGCA'
      return
#endif

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

      SUBROUTINE conserv_DIAG (M,CONSFN,ICON)
!@sum  conserv_DIAG generic routine keeps track of conserved properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE GEOM, only : j_budg, j_0b, j_1b, imaxj
      USE DIAG_COM, only : consrv=>consrv_loc,nofm, jm_budg,wtbudg
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE
!@var M index denoting from where routine is called
      INTEGER, INTENT(IN) :: M
!@var ICON index for the quantity concerned
      INTEGER, INTENT(IN) :: ICON
!@var CONSFN external routine that calculates total conserved quantity
      EXTERNAL CONSFN
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TOTAL
      REAL*8, DIMENSION(JM_BUDG) :: TOTALJ
      INTEGER :: I,J,NM,NI
      INTEGER :: I_0,I_1,J_0,J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C**** NOFM contains the indexes of the CONSRV array where each
C**** change is to be stored for each quantity. If NOFM(M,ICON)=0,
C**** no calculation is done.
C**** NOFM(1,ICON) is the index for the instantaneous value.
      IF (NOFM(M,ICON).gt.0) THEN
C**** Calculate current value TOTAL
        CALL CONSFN(TOTAL)
        NM=NOFM(M,ICON)
        NI=NOFM(1,ICON)
C**** Calculate zonal sums
        TOTALJ(J_0B:J_1B)=0.
        DO J=J_0,J_1
          DO I=I_0,IMAXJ(J) ! not I_1 b/c latlon wtbudg differs at poles
            TOTALJ(J_BUDG(I,J)) = TOTALJ(J_BUDG(I,J)) + TOTAL(I,J)
     &           *WTBUDG(I,J)
          END DO
        END DO
C**** Accumulate difference from last time in CONSRV(NM)
        IF (M.GT.1) THEN
          DO J=J_0B,J_1B
            CONSRV(J,NM)=CONSRV(J,NM)+(TOTALJ(J)-CONSRV(J,NI))
          END DO
        END IF
C**** Save current value in CONSRV(NI)
        DO J=J_0B,J_1B
          CONSRV(J,NI)=TOTALJ(J)
        END DO
      END IF
      RETURN
C****
      END SUBROUTINE conserv_DIAG


      SUBROUTINE conserv_MS(RMASS)
!@sum  conserv_MA calculates total atmospheric mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,p,pstrat
      USE GEOM, only : imaxj
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RMASS
      INTEGER :: I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0,    J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** MASS
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        RMASS(I,J)=(P(I,J)+PSTRAT)*mb2kg
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) RMASS(2:im,1) =RMASS(1,1)
      IF(HAVE_NORTH_POLE) RMASS(2:im,JM)=RMASS(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_MS


      SUBROUTINE conserv_PE(TPE)
!@sum  conserv_TPE calculates total atmospheric potential energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,t,p,ptop,zatmo
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pk,pdsig
      USE DOMAIN_DECOMP_ATM, only : GET,GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TPE
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL POTENTIAL ENERGY (J/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        TPE(I,J)=0.
        DO L=1,LM
          TPE(I,J)=TPE(I,J)+T(I,J,L)*PK(L,I,J)*PDSIG(L,I,J)
        ENDDO
        TPE(I,J)=(TPE(I,J)*SHA+ZATMO(I,J)*(P(I,J)+PTOP))*mb2kg
      ENDDO
      ENDDO
      IF(HAVE_SOUTH_POLE) TPE(2:im,1) =TPE(1,1)
      IF(HAVE_NORTH_POLE) TPE(2:im,JM)=TPE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_PE

      SUBROUTINE conserv_WM(WATER)
!@sum  conserv_WM calculates total atmospheric water mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,wm,q
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE

      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: WATER
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL WATER MASS (kg/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        WATER(I,J) = 0.
        DO L=1,LM
          WATER(I,J)=WATER(I,J)+(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)
        ENDDO
        WATER(I,J)=WATER(I,J)*mb2kg
      ENDDO
      ENDDO
      IF (HAVE_SOUTH_POLE) WATER(2:im,1) = WATER(1,1)
      IF (HAVE_NORTH_POLE) WATER(2:im,JM)= WATER(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_WM


      SUBROUTINE conserv_EWM(EWATER)
!@sum  conserv_EWM calculates total atmospheric water energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg,shv,grav,lhe
      USE MODEL_COM, only : im,jm,lm,wm,t,q,p
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig, pmid, pk
      USE CLOUDS_COM, only : svlhx
      USE DOMAIN_DECOMP_ATM, only : GET, GRID
      IMPLICIT NONE
      REAL*8, PARAMETER :: HSCALE = 7.8d0 ! km
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: EWATER
      INTEGER :: I,J,L
      INTEGER :: J_0,J_1,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      REAL*8 EL!,W

      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** TOTAL WATER ENERGY (J/m^2)
C****
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EWATER(I,J) = 0.
        DO L=1,LM
c this calculation currently only calculates latent heat
c          W =(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)*mb2kg
          EL=(Q(I,J,L)*LHE+WM(I,J,L)*(LHE-SVLHX(L,I,J)))*PDSIG(L,I,J)
          EWATER(I,J)=EWATER(I,J)+EL !+W*(SHV*T(I,J,L)*PK(L,I,J)+GRAV
!     *           *HSCALE*LOG(P(I,J)/PMID(L,I,J)))
        ENDDO
        EWATER(I,J)=EWATER(I,J)*mb2kg
      ENDDO
      ENDDO
      IF(HAVE_SOUTH_POLE) EWATER(2:im,1) = EWATER(1,1)
      IF(HAVE_NORTH_POLE) EWATER(2:im,JM)= EWATER(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_EWM

      SUBROUTINE DIAG4A
C****
C**** THIS ROUTINE PRODUCES A TIME HISTORY OF ENERGIES
C****
      USE MODEL_COM, only : im,istrat,IDACC
      USE DIAG_COM, only : energy,speca,ned
      IMPLICIT NONE

      INTEGER :: I,IDACC5,N,NM

#ifdef SCM
c     write(0,*) 'SCM no diags    DIAG4A'
      return
#endif

      IF (IDACC(4).LE.0.OR.IDACC(7).LE.0) RETURN
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
      USE FILEMANAGER, only : openunit, closeunit, nameunit
      USE DIAG_COM, only : kgz_max,pmname,P_acc,PM_acc
      USE PARAM
#ifdef CALCULATE_FLAMMABILITY
      use flammability_com, only : raP_acc
#endif
#ifdef TRACERS_ON
#ifndef SKIP_TRACER_DIAGS
      USE TRACER_COM, only : ntm, trm, trname
     *     , mass2vol, n_Ox, n_SO4,
     *     n_SO4_d1,n_SO4_d2,n_SO4_d3,n_clay,n_clayilli,n_sil1quhe,
     *     n_water, n_HDO, n_Be7, n_NOx, n_CO
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     *     ,Ntm_dust
#endif
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#ifdef TRACERS_WATER
     &     ,dowetdep, trw0
#endif
#endif /* SKIP_TRACER_DIAGS */
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7D_acc,BE7W_acc
#endif
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
!@var lst level strings
      character*2, dimension(lm) :: lst

      contains

      subroutine init_subdd(aDATE)
!@sum init_subdd initialise sub daily diags and position files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate
      integer :: i,j,k,l,kunit

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
        call open_subdd(aDATE)

C**** position correctly
        do kunit=1,kddunit
          call io_POS(iu_SUBDD(kunit),Itime-1,im*jm,Nsubdd)
        end do

      end if

C**** define lst
      do l=1,lm
        if (l.lt.10) write(lst(l)(1:2),'(I1,1X)') l
        if (l.ge.10) write(lst(l)(1:2),'(I2)') l
      end do

C**** initialise special subdd accumulation
      P_acc=0.
      PM_acc=0.
#ifdef TRACERS_COSMO
      BE7W_acc=0.
      BE7D_acc=0.
#endif

      return
      end subroutine init_subdd

      subroutine open_subdd(aDATE)
!@sum open_subdd opens sub daily diag files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate
      character*12 name
      integer :: k,kunit,kk

      kunit=0
      do k=1,kdd
C**** Some names have more than one unit associated (i.e. "ZALL")
        if (namedd(k)(len_trim(namedd(k))-2:len_trim(namedd(k))).eq.
     *       "ALL") then
          select case (namedd(k)(1:1))
          case ("U","V","W","C","D","O","B","N","t","q","z","r","n","c")
            ! velocities/tracers on model layers
            kunit=kunit+1
            write(name,'(A1,A3,A7)') namedd(k)(1:1),'ALL',aDATE(1:7)
            call openunit(name,iu_SUBDD(kunit),.true.,.false.)
          case ("Z", "T", "R", "Q") ! heights, temps, rel/spec hum PMB levels
            do kk=1,kgz_max
              kunit=kunit+1
              call openunit(namedd(k)(1:1)//trim(PMNAME(kk))//
     *             aDATE(1:7),iu_SUBDD(kunit),.true.,.false.)
            end do
          end select
          select case (namedd(k)(1:2))
          case ("GT","GW","GI")
            ! soil variables on soil levels
            kunit=kunit+1
            write(name,'(A2,A3,A7)') namedd(k)(1:2),'ALL',aDATE(1:7)
            call openunit(name,iu_SUBDD(kunit),.true.,.false.)
          end select
        else                    ! single file per name
          kunit=kunit+1
          call openunit(trim(namedd(k))//aDATE(1:7),iu_SUBDD(kunit),
     *         .true.,.false.)
        endif
      end do
      kddunit=kunit
C****
      return
      end subroutine open_subdd

      subroutine reset_subdd(aDATE)
!@sum reset_subdd resets sub daily diag files
!@auth Gavin Schmidt
      implicit none
      character*14, intent(in) :: adate

      if (nsubdd.ne.0) then
C**** close and re-open units
        call closeunit ( iu_SUBDD(1:kddunit) )
        call open_subdd( aDATE )
      end if
C****
      return
      end subroutine reset_subdd

      subroutine get_subdd
!@sum get_SUBDD saves instantaneous variables at sub-daily frequency
!@+   every ABS(NSUBDD)
!@+   Note that TMIN, TMAX, AOD can only be output once a day
!@+   Current options: SLP, PS, SAT, PREC, QS, LCLD, MCLD, HCLD, PTRO
!@+                    QLAT, QSEN, SWD, SWU, LWD, LWU, LWT, STX, STY,
!@+                    ICEF, SNOWD, TCLD, SST, SIT, US, VS, TMIN, TMAX
!@+                    MCP, SNOWC, RS, GT1, GTD, GW0, GWD, GI0, GID,
!@+                    GTALL, GWALL, GIALL (on soil levels)
!@+                    Z*, R*, T*, Q*  (on any fixed pressure level)
!@+                    z*, r*, t*, q*  (on any model level, note lowercase)
!@+                    U*, V*, W*, C*  (on any model level)
!@+                    O*, N*, c*, n*  (Ox, NOx, CO, NO2 on any model level)
!@+                    D*          (HDO on any model level)
!@+                    B*          (BE7 on any model level)
!@+                    SO4, RAPR
!@+                    7BEW, 7BED, BE7ATM
!@+                    CTEM,CD3D,CI3D,CL3D,CDN3D,CRE3D,CLWP  ! aerosol
!@+                    TAUSS,TAUMC,CLDSS,CLDMC
!@+                    SO4_d1,SO4_d2,SO4_d3,   ! het. chem
!@+                    Clay, Silt1, Silt2, Silt3  ! dust
!@+                    DUEMIS,DUDEPTURB,DUDEPGRAV,DUDEPWET,DUTRS,DULOAD
!@+                    DUEMIS2
!@+                    AOD aer opt dep (1,NTRACE in rad code) daily avg
!@+
!@+   More options can be added as extra cases in this routine
!@auth Gavin Schmidt/Reto Ruedy
      USE CONSTANT, only : grav,rgas,bygrav,bbyg,gbyrb,sday,tf,mair,sha
     *     ,lhe,rhow,undef,stbo
      USE MODEL_COM, only : lm,p,ptop,zatmo,dtsrc,u,v,focean
     *     ,flice,nday,t,q
      USE GEOM, only : imaxj,axyp,byaxyp
      USE PBLCOM, only : tsavg,qsavg,usavg,vsavg
      USE GHY_COM, only : gdeep,gsaveL,ngm
      USE CLOUDS_COM, only : llow,lmid,lhi,cldss,cldmc,taumc,tauss,fss
#ifdef CLD_AER_CDNC
     *           ,ctem,cd3d,ci3d,cl3d,cdn3d,cre3d,clwp
#endif
      USE DYNAMICS, only : ptropo,am,wsave,pk,phi,pmid
      USE FLUXES, only : prec,dmua,dmva,tflux1,qflux1,uflux1,vflux1
     *     ,gtemp,gtempr
#ifdef TRACERS_ON
     *     ,trcsurf
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,dust_flux_glob
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
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only : mNO2
#endif
      USE SEAICE_COM, only : rsi,snowi
      USE LANDICE_COM, only : snowli
      USE LAKES_COM, only : flake
      USE GHY_COM, only : snowe,fearth,wearth,aiearth
      USE RAD_COM, only : trhr,srdn,salb,cfrac,cosz1
#ifdef HTAP_LIKE_DIAGS
     & ,ttausv_sum,ttausv_count,ntrix
#endif
#ifdef HTAP_LIKE_DIAGS
      USE RADPAR, only : NTRACE
#endif
      USE DIAG_COM, only : z_inst,rh_inst,t_inst,kgz_max,pmname,tdiurn
     *     ,p_acc,pm_acc,pmb
      USE DOMAIN_DECOMP_ATM, only : GRID,GET,am_i_root
#ifdef TRACERS_ON
      USE TRACER_COM
#endif
      IMPLICIT NONE
      REAL*4, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: DATA
      INTEGER :: I,J,K,L,kp,ks,kunit,n,n1,n_fidx
      REAL*8 POICE,PEARTH,PLANDI,POCEAN,QSAT,PS,SLP, ZS
      INTEGER :: J_0,J_1,J_0S,J_1S,I_0,I_1
      LOGICAL :: polefix,have_south_pole,have_north_pole

      CALL GET(GRID,J_STRT=J_0, J_STOP=J_1,
     &              J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=have_south_pole,
     &               HAVE_NORTH_POLE=have_north_pole)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

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

      kunit=0
C**** depending on namedd string choose what variables to output
      nameloop: do k=1,kdd

C**** simple diags (one record per file)
        select case (namedd(k))
        case ("SLP")            ! sea level pressure (mb)
          do j=J_0,J_1
          do i=I_0,imaxj(j)
            ps=(p(i,j)+ptop)
            zs=bygrav*zatmo(i,j)
            data(i,j)=slp(ps,tsavg(i,j),zs)
          end do
          end do
        case ("PS")             ! surface pressure (mb)
          data=p+ptop
        case ("SAT")            ! surf. air temp (C)
          data=tsavg-tf
        case ("US")             ! surf. u wind (m/s)
          data=usavg
        case ("VS")             ! surf. v wind (m/s)
          data=vsavg
        case ("SST")            ! sea surface temp (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (FOCEAN(I,J)+FLAKE(I,J).gt.0) then
                data(i,j)=GTEMP(1,1,i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("SIT")            ! surface ice temp (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J)).gt.0) then
                data(i,j)=GTEMP(1,2,i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("GT1")      ! level 1 ground temp (LAND) (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                data(i,j)=gtemp(1,4,i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("GTD")   ! avg levels 2-6 ground temp (LAND) (C)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                data(i,j)=gdeep(i,j,1)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("GWD")  ! avg levels 2-6 ground liq water (m)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                data(i,j)=gdeep(i,j,2)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("GID")  ! avg levels 2-6 ground ice (m liq. equiv.)
          do j=J_0,J_1          
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                data(i,j)=gdeep(i,j,3)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("GW0")  ! ground lev 1 + canopy liq water (m)
          do j=J_0,J_1          
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                data(i,j)=wearth(i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("GI0")  ! ground lev 1 + canopy ice (m liq. equiv.)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if (fearth(i,j).gt.0) then
                data(i,j)=aiearth(i,j)
              else
                data(i,j)=undef
              end if
            end do
          end do
        case ("QS")             ! surf spec humidity (kg/kg)
          data=qsavg
        case ("RS")             ! surf rel humidity
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              data(i,j)=qsavg(i,j)/qsat(tsavg(i,j),lhe,p(i,j)+ptop)
            enddo
          enddo
        case ("PREC")           ! precip (mm/day)
c          data=sday*prec/dtsrc
          data=sday*P_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          P_acc=0.
#ifdef CALCULATE_FLAMMABILITY
        case ("RAPR")   !running avg precip (mm/day)
          data=sday*raP_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          raP_acc=0.
#endif
        case ("MCP")       ! moist conv precip (mm/day)
          data=sday*PM_acc/(Nsubdd*dtsrc) ! accum over Nsubdd steps
          PM_acc=0.
        case ("SNOWD")     ! snow depth (w.e. mm)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              data(i,j)=1d3*(SNOWI(I,J)*POICE+SNOWLI(I,J)*PLANDI+SNOWE(I
     *             ,J)*PEARTH)/RHOW
            end do
          end do
        case ("SNOWC")     ! snow cover (fraction of grid)
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              data(i,j)=0.d0
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              if(SNOWI(I,J) > 0.)data(i,j)=data(i,j)+POICE
              PEARTH=FEARTH(I,J)
              if(SNOWE(I,J) > 0.)data(i,j)=data(i,j)+PEARTH
              PLANDI=FLICE(I,J)
              if(SNOWLI(I,J) > 0.)data(i,j)=data(i,j)+PLANDI
              data(i,j)=min(1.0,data(i,j))
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
            do i=I_0,imaxj(j)
              POCEAN=(1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              data(i,j)=STBO*(POCEAN*GTEMPR(1,I,J)**4+
     *             POICE *GTEMPR(2,I,J)**4+PLANDI*GTEMPR(3,I,J)**4+
     *             PEARTH*GTEMPR(4,I,J)**4)
            end do
          end do
        case ("LWT")            ! LW upward flux at TOA (P1) (W/m^2)
          do j=J_0,J_1     ! sum up all cooling rates + net surface emission
            do i=I_0,imaxj(j)
              POCEAN=(1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))
              POICE=RSI(I,J)*(FOCEAN(I,J)+FLAKE(I,J))
              PEARTH=FEARTH(I,J)
              PLANDI=FLICE(I,J)
              data(i,j)=-SUM(TRHR(0:LM,I,J))+
     *             STBO*(POCEAN*GTEMPR(1,I,J)**4+
     *             POICE *GTEMPR(2,I,J)**4+PLANDI*GTEMPR(3,I,J)**4+
     *             PEARTH*GTEMPR(4,I,J)**4)
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
            do i=I_0,imaxj(j)
              do l=1,llow
                data(i,j)=data(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              data(i,j)=data(i,j)*100./real(llow,kind=8)
            end do
          end do
        case ("MCLD")           ! mid level cloud cover (%)
          data=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              do l=llow+1,lmid
                data(i,j)=data(i,j)+(cldss(l,i,j)+cldmc(l,i,j))
              end do
              data(i,j)=data(i,j)*100./real(lmid-llow,kind=8)
            end do
          end do
        case ("HCLD")           ! high level cloud cover (%)
          data=0.               ! Warning: these can be greater >100!
          do j=J_0,J_1
            do i=I_0,imaxj(j)
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
          data=trm(:,:,1,n_SO4)
#ifdef TRACERS_HETCHEM
     *     +trm(:,:,1,n_SO4_d1)+trm(:,:,1,n_SO4_d2)+trm(:,:,1,n_SO4_d3)
#endif
#endif
#ifdef CLD_AER_CDNC
        case ("CLWP")             !LWP (kg m-2)
          data=clwp
#endif
#ifdef TRACERS_COSMO
        case ("7BEW")
          data=Be7w_acc
          Be7w_acc=0.

        case ("7BED")
          data=Be7d_acc
          Be7d_acc=0.

        case ("7BES")
          data=1.d6*trcsurf(:,:,n_Be7)   ! 10^-6 kg/kg
#endif
        case default
          goto 10
        end select
        kunit=kunit+1
        polefix=.true.
        call write_data(data,kunit,polefix)
        cycle
 10     continue

C**** diags on soil levels
        select case (namedd(k)(1:2))
        case ("GT","GW","GI") ! soil temperature, water, ice
          if (namedd(k)(3:5) .eq. "ALL") then
            kunit=kunit+1
            do ks=1,ngm
              select case (namedd(k)(1:2))
              case ("GT")        ! soil temperature lvls 1-6, land (C)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    if (fearth(i,j).gt.0) then
                      data(i,j)=gsaveL(i,j,ks,1)
                    else
                      data(i,j)=undef
                    end if
                  end do
                end do
              case ("GW")        ! ground wetness lvls 1-6, land (m)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    if (fearth(i,j).gt.0) then
                      data(i,j)=gsaveL(i,j,ks,2)
                    else
                      data(i,j)=undef
                    end if
                  end do
                end do
              case ("GI")  ! ground ice lvls 1-6, land (m, liq equiv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    if (fearth(i,j).gt.0) then
                      data(i,j)=gsaveL(i,j,ks,3)
                    else
                      data(i,j)=undef
                    end if
                  end do
                end do
              end select
              polefix=.true.
              call write_data(data,kunit,polefix)
            end do
            cycle
          end if
        end select

C**** diags on fixed pressure levels or velocity
        select case (namedd(k)(1:1))
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
                do i=I_0,imaxj(j)
                  data(i,j)=rh_inst(kp,i,j)*qsat(t_inst(kp,i,j),lhe
     *                 ,PMB(kp))
                end do
                end do
              case ("T")        ! temperature (C)
                data=t_inst(kp,:,:)
              end select
              polefix=.true.
              call write_data(data,kunit,polefix)
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
                do i=I_0,imaxj(j)
                  data(i,j)=rh_inst(kp,i,j)*qsat(t_inst(kp,i,j),lhe
     *                 ,PMB(kp))
                end do
                end do
              case ("T")        ! temperature (C)
                data=t_inst(kp,:,:)
              end select
              polefix=.true.
              call write_data(data,kunit,polefix)
            end do
            cycle
          end if

C**** diagnostics on model levels
        case ("U","V","W","C","O","B","D","N","t","q","z","r","c","n")
             ! velocity/clouds/tracers, temp,spec.hum.,geo.ht
          if (namedd(k)(2:4) .eq. "ALL") then
            kunit=kunit+1
            do kp=1,LmaxSUBDD
              select case (namedd(k)(1:1))
              case ("t")        ! temperature (C)
                if(have_south_pole) data(1:im,1)=
     &          t(1,1,kp)*pk(kp,1, 1)-tf
                if(have_north_pole) data(1:im,jm)=
     &          t(1,jm,kp)*pk(kp,1,jm)-tf
                data(:,J_0S:J_1S)=
     &          t(:,J_0S:J_1S,kp)*pk(kp,:,J_0S:J_1S)-tf
              case ("r")        ! relative humidity
                if(have_south_pole) data(1:im, 1)=q(1,1,kp)/
     &          qsat(t(1,1,kp)*pk(kp,1,1),lhe,pmid(kp,1,1))
                if(have_north_pole) data(1:im,jm)=q(1,jm,kp)/
     &          qsat(t(1,jm,kp)*pk(kp,1,jm),lhe,pmid(kp,1,jm))
                do j=J_0S,J_1S; do i=i_0,i_1
                  data(i,j)=q(i,j,kp)/qsat(t(i,j,kp)*pk(kp,i,j),
     &            lhe,pmid(kp,i,j))
                enddo         ; enddo
              case ("q")        ! specific humidity
                if(have_south_pole) data(1:im, 1)=q(1, 1,kp)
                if(have_north_pole) data(1:im,jm)=q(1,jm,kp)
                data(:,J_0S:J_1S)=q(:,J_0S:J_1S,kp)
              case ("z")        ! geopotential height
                if(have_south_pole) data(1:im, 1)=phi(1, 1,kp)
                if(have_north_pole) data(1:im,jm)=phi(1,jm,kp)
                data(:,J_0S:J_1S)=phi(:,J_0S:J_1S,kp)
              case ("U")        ! E-W velocity
                data=u(:,:,kp)
              case ("V")        ! N-S velocity
                data=v(:,:,kp)
              case ("W")        ! vertical velocity
                data=wsave(:,:,kp)
              case ("C")        ! estimate of cloud optical depth
                data=(1.-fss(kp,:,:))*taumc(kp,:,:)+fss(kp,:,:)
     *               *tauss(kp,:,:)
#ifdef TRACERS_SPECIAL_Shindell
              case ("O")                ! Ox ozone tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*trm(i,j,kp,n_Ox)*mass2vol(n_Ox)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
              case ("N")                ! NOx tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*trm(i,j,kp,n_NOx)*mass2vol(n_NOx)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
              case ("c")                ! CO tracer (ppbv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d9*trm(i,j,kp,n_CO)*mass2vol(n_CO)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
              case ("n")                ! NO2 (not a tracer) (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*mNO2(i,j,kp)*mair/(46.0055
     *                   *am(kp,i,j)*axyp(i,j))
                  end do
                end do
#endif
#ifdef TRACERS_COSMO
              case ("B")                ! Be7 tracer
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*trm(i,j,kp,n_Be7)* mass2vol(n_Be7)/
     *                   (am(kp,i,j)*axyp(i,j))
                  end do
                end do
#endif
#ifdef TRACERS_SPECIAL_O18
              case ("D")                ! HDO tracer (permil)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1d3*(trm(i,j,kp,n_HDO)/(trm(i,j,kp,n_water
     *                   )*trw0(n_HDO))-1.)
                  end do
                end do
#endif
              end select
              polefix=(namedd(k)(1:1).ne."U".and.namedd(k)(1:1).ne."V")
              call write_data(data,kunit,polefix)
            end do
            cycle
          end if
C**** get model level
          do l=1,lm
            if (trim(namedd(k)(2:5)) .eq. lst(l)) then
              kunit=kunit+1
              select case (namedd(k)(1:1))
              case ("t")        ! temperature (C)
                if(have_south_pole) data(1:im,1)=
     &               t(1,1,l)*pk(l,1, 1)-tf
                if(have_north_pole) data(1:im,jm)=
     &               t(1,jm,l)*pk(l,1,jm)-tf
                data(:,J_0S:J_1S)=
     &               t(:,J_0S:J_1S,l)*pk(l,:,J_0S:J_1S)-tf
              case ("r")        ! relative humidity
                if(have_south_pole) data(1:im, 1)=q(1, 1,l)/
     &               qsat(t(1,1,l)*pk(l,1,1),lhe,pmid(l,1,1))
                if(have_north_pole) data(1:im,jm)=q(1,jm,l)/
     &               qsat(t(1,jm,l)*pk(l,1,jm),lhe,pmid(l,1,jm))
                do j=J_0S,J_1S; do i=i_0,i_1
                  data(i,j)=q(i,j,l)/qsat(t(i,j,l)*pk(l,i,j),
     &                 lhe,pmid(l,i,j))
                enddo
              enddo
              case ("q")        ! specific humidity
                if(have_south_pole) data(1:im, 1)=q(1, 1,l)
                if(have_north_pole) data(1:im,jm)=q(1,jm,l)
                data(:,J_0S:J_1S)=q(:,J_0S:J_1S,l)
              case ("z")        ! geopotential height
                if(have_south_pole) data(1:im, 1)=phi(1, 1,l)
                if(have_north_pole) data(1:im,jm)=phi(1,jm,l)
                data(:,J_0S:J_1S)=phi(:,J_0S:J_1S,l)
              case ("U")        ! U velocity
                data=u(:,:,l)
              case ("V")        ! V velocity
                data=v(:,:,l)
              case ("W")        ! W velocity
                data=wsave(:,:,l)
              case ("C")        ! estimate of cloud optical depth
                data=(1.-fss(l,:,:))*taumc(l,:,:)+fss(l,:,:)
     *               *tauss(l,:,:)
#ifdef TRACERS_SPECIAL_Shindell
              case ("O")                ! Ox ozone tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*trm(i,j,l,n_Ox)*mass2vol(n_Ox)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
              case ("N")                ! NOx tracer (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*trm(i,j,l,n_NOx)*mass2vol(n_NOx)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
              case ("c")                ! CO tracer (ppbv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d9*trm(i,j,l,n_CO)*mass2vol(n_CO)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
              case ("n")                ! NO2 (not a tracer) (ppmv)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*mNO2(i,j,l)*mair/(46.0055
     *                   *am(l,i,j)*axyp(i,j))
                  end do
                end do
#endif
#ifdef TRACERS_COSMO
              case ("B")                ! Be7 tracer
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1.d6*trm(i,j,l,n_Be7)* mass2vol(n_Be7)/
     *                   (am(l,i,j)*axyp(i,j))
                  end do
                end do
#endif
#ifdef TRACERS_SPECIAL_O18
              case ("D")                ! HDO tracer (permil)
                do j=J_0,J_1
                  do i=I_0,imaxj(j)
                    data(i,j)=1d3*(trm(i,j,l,n_HDO)/(trm(i,j,l,n_water
     *                   )*trw0(n_HDO))-1.)
                  end do
                end do
#endif
              end select
              polefix=(namedd(k)(1:1).ne."U".and.namedd(k)(1:1).ne."V")
              call write_data(data,kunit,polefix)
              cycle nameloop
            end if
          end do
        end select

C**** Additional diags - multiple records per file
        select case (namedd(k))

C**** cases using all levels up to LmaxSUBDD
          case ("SO2", "SO4", "SO4_d1", "SO4_d2", "SO4_d3", "Clay",
     *         "Silt1", "Silt2", "Silt3", "CTEM", "CL3D", "CI3D", "CD3D"
     *         , "CLDSS", "CLDMC", "CDN3D", "CRE3D", "TAUSS", "TAUMC",
     *         "CLWP")
          kunit=kunit+1
          do l=1,LmaxSUBDD
            select case(namedd(k))
#ifdef TRACERS_HETCHEM
            case ("SO2")
              data=trm(:,:,l,n_SO2)
            case ("SO4")
              data=trm(:,:,l,n_SO4)
            case ("SO4_d1")
              data=trm(:,:,l,n_SO4_d1)
            case ("SO4_d2")
              data= trm(:,:,l,n_SO4_d2)
            case ("SO4_d3")
              data=trm(:,:,l,n_SO4_d3)
            case ("Clay")
              data=trm(:,:,l,n_Clay)
            case ("Silt1")
              data=trm(:,:,l,n_Silt1)
            case ("Silt2")
              data=trm(:,:,l,n_Silt2)
            case ("Silt3")
              data=trm(:,:,l,n_Silt3)
#endif
#ifdef CLD_AER_CDNC
            case ("CTEM")
              data=ctem(l,:,:) ! cld temp (K) at cld top
            case ("CL3D")
              data=cl3d(l,:,:) ! cld LWC (kg m-3)
            case ("CI3D")
              data=ci3d(l,:,:) ! cld IWC (kg m-3)
            case ("CD3D")
              data=cd3d(l,:,:) ! cld thickness (m)
            case ("CLDSS")
              data=100.d0*cldss(l,:,:) ! Cld cover LS(%)
            case ("CLDMC")
              data=100.d0*cldmc(l,:,:) ! Cld cover MC(%)
            case ("CDN3D")
              data=cdn3d(l,:,:) ! cld CDNC (cm^-3)
            case ("CRE3D")
              data=1.d-6*cre3d(l,:,:) ! cld Reff (m)
            case ("TAUSS")
              data=tauss(l,:,:) ! LS cld tau
            case ("TAUMC")
              data=taumc(l,:,:) ! MC cld tau
#endif
            end select
            polefix=.true.
            call write_data(data,kunit,polefix)
          end do
          cycle

#ifdef HTAP_LIKE_DIAGS
C**** for AOD multiple tracers are written to one file 
          case ('AOD') ! aerosol optical depths daily average
            kunit=kunit+1
            if(mod(itime+1,Nday).ne.0) cycle ! only at end of day
            if(ttausv_count==0.)call stop_model('ttausv_count=0',255)
            do n=1,NTRACE
              if(ntrix(n) > 0)then
                data=ttausv_sum(:,:,n)/ttausv_count
                polefix=.true.
                call write_data(data,kunit,polefix)
              endif
            enddo
            cycle
#endif

C**** cases where multiple records go to one file for dust

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          case ('DUEMIS','DUEMIS2','DUTRS', 'DULOAD')
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_fidx+n-1
C**** first set: no 'if' tests
            select case (namedd(k))

            CASE ('DUEMIS')     ! Dust emission flux [kg/m^2/s]
              data=dust_flux_glob(:,:,n)
#ifdef TRACERS_DUST
            CASE ('DUEMIS2')    ! Dust emission flux 2 (diag. var. only) [kg/m^2/s]
              data=dust_flux2_glob(:,:,n)
#endif
            CASE ('DUTRS')      ! Mixing ratio of dust tracers at surface [kg/kg]
              data=trcsurf(:,:,n1)
            CASE ('DULOAD')     ! Dust load [kg/m^2]
              data=0.D0
              DO j=J_0,J_1
                DO l=1,LmaxSUBDD
                  data(:,j)=data(:,j)+trm(:,j,l,n1)
                END DO
                data(:,j)=data(:,j)*byaxyp(:,j)
              END DO

            end select
            polefix=.true.
            call write_data(data,kunit,polefix)
          end do
          cycle

C**** other dust special cases

#ifdef TRACERS_DRYDEP
          CASE ('DUDEPTURB')        ! Turb. deposition flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_fidx+n-1
            IF (dodrydep(n1)) THEN
              data=depo_turb_glob(:,:,1,n1)+depo_turb_glob(:,:,2,n1)
     *             +depo_turb_glob(:,:,3,n1)+depo_turb_glob(:,:,4,n1)
              polefix=.true.
              call write_data(data,kunit,polefix)
            END IF
          end do
          cycle

          CASE ('DUDEPGRAV')      ! Gravit. settling flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_fidx+n-1
            IF (dodrydep(n1)) THEN
              data=depo_grav_glob(:,:,1,n1)+depo_grav_glob(:,:,2,n1)
     *             +depo_grav_glob(:,:,3,n1)+depo_grav_glob(:,:,4,n1)
              polefix=.true.
              call write_data(data,kunit,polefix)
            END IF
          end do
          cycle
#endif
          CASE ('DUDEPWET')         ! Wet deposition flux of dust tracers [kg/m^2/s]
          kunit=kunit+1
          do n=1,Ntm_dust
            n1=n_fidx+n-1
#ifdef TRACERS_WATER
            IF (dowetdep(n1)) THEN
              DO j=J_0,J_1
                data(:,j)=trprec(n1,:,j)*byaxyp(:,j)/Dtsrc
              END DO
#else
              DO j=J_0,J_1
                data(:,j)=trprec_dust(n,:,j)*byaxyp(:,j)/Dtsrc
              END DO
#endif
              polefix=.true.
              call write_data(data,kunit,polefix)
#ifdef TRACERS_WATER
            END IF
#endif
          end do
          cycle
#endif
C**** this prevents tokens that are not caught from messing up the file data
        case default
          kunit=kunit+1
        end select

      end do nameloop ! end of  k=1,kdd loop
c****
      return
      end subroutine get_subdd

      subroutine write_data(data,kunit,polefix)
!@sum write out subdd data array with optional pole fix
      use domain_decomp_1d, only : grid,get,writei_parallel
! need to add writei_parallel to domain_decomp_atm
      real*4, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: data
      integer kunit
      logical :: polefix

c**** fix polar values
      if (polefix) then
        if(grid%have_south_pole) data(2:im,1) =data(1,1)
        if(grid%have_north_pole) data(2:im,jm)=data(1,jm)
      end if
      call writei_parallel(grid,iu_subdd(kunit),
     *     nameunit(iu_subdd(kunit)),data,itime)

      end subroutine write_data

      end module subdaily

      subroutine ahourly
!@sum ahourly saves instantaneous variables at sub-daily frequency
!@+   for diurnal cycle diagnostics
!@auth Reha Cakmur/Jan Perlwitz

      USE MODEL_COM, only : u,v,t,p,q,jdate,jhour,ptop,sig
      USE CONSTANT, only : bygrav
      USE domain_decomp_atm, ONLY : am_i_root,get,globalsum,grid
      USE GEOM, only : imaxj,axyp,byaxyp
      USE DYNAMICS, only : phi,wsave,pek,byam
      USE rad_com,ONLY : cosz1,srnflb_save,trnflb_save,ttausv_save,
     &     ttausv_cs_save
      USE diag_com,ONLY : adiurn_dust,ndiupt,ndiuvar,ijdd
     &     ,adiurn=>adiurn_loc
#ifndef NO_HDIURN
     &     ,hdiurn=>hdiurn_loc
#endif
#ifdef TRACERS_DUST
     *     ,idd_u1,idd_v1,idd_uv1,idd_t1,idd_qq1,idd_p1,idd_w1,idd_phi1
     *     ,idd_sr1,idd_tr1,idd_load1,idd_conc1,idd_tau1,idd_tau_cs1
#endif
#ifdef TRACERS_ON
      USE TRACER_COM, only : trm
#ifdef TRACERS_DUST
     &     ,Ntm_dust,n_clay
#endif
#endif

      IMPLICIT NONE

      INTEGER :: i,j,ih,ihm,kr,n,n1
      REAL*8 :: psk
      INTEGER,PARAMETER :: lmax_dd2=11, n_idxd=14*lmax_dd2
      INTEGER :: idxd(n_idxd)
      REAL*8 :: tmp(NDIUVAR)

C****   define local grid
      INTEGER J_0, J_1, I_0, I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL get(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN

        idxd=(/
     *       (idd_u1+i-1,i=1,lmax_dd2), (idd_v1+i-1,i=1,lmax_dd2),
     *       (idd_uv1+i-1,i=1,lmax_dd2), (idd_t1+i-1,i=1,lmax_dd2),
     *       (idd_qq1+i-1,i=1,lmax_dd2), (idd_p1+i-1,i=1,lmax_dd2),
     *       (idd_w1+i-1,i=1,lmax_dd2), (idd_phi1+i-1,i=1,lmax_dd2),
     *       (idd_sr1+i-1,i=1,lmax_dd2), (idd_tr1+i-1,i=1,lmax_dd2),
     *       (idd_load1+i-1,i=1,lmax_dd2), (idd_conc1+i-1,i=1,lmax_dd2),
     *       (idd_tau1+i-1,i=1,lmax_dd2), (idd_tau_cs1+i-1,i=1,lmax_dd2)
     *       /)

      END IF
#endif

      ih=jhour+1
      ihm=ih+(jdate-1)*24
!$OMP PARALLEL DO PRIVATE(i,j,kr,n,psk,n1,tmp)
!$OMP*   SCHEDULE(DYNAMIC,2)
      do j=j_0,j_1
      do i=I_0,imaxj(j)
      psk=pek(1,i,j)
      do kr=1,ndiupt
        if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then
#ifdef TRACERS_DUST
          IF (adiurn_dust == 1) THEN
            tmp=0.D0

            tmp(idd_u1:idd_u1+lmax_dd2-1)=u(i,j,1:lmax_dd2)
            tmp(idd_v1:idd_v1+lmax_dd2-1)=v(i,j,1:lmax_dd2)
            tmp(idd_uv1:idd_uv1+lmax_dd2-1)=sqrt(u(i,j,1:lmax_dd2)*u(i
     *           ,j,1:lmax_dd2)+v(i,j,1:lmax_dd2)*v(i,j,1:lmax_dd2))
            tmp(idd_t1:idd_t1+lmax_dd2-1)=t(i,j,1:lmax_dd2)*psk
            tmp(idd_qq1:idd_qq1+lmax_dd2-1)=q(i,j,1:lmax_dd2)
            tmp(idd_p1:idd_p1+lmax_dd2-1)=p(i,j)*sig(1:lmax_dd2)+ptop
            tmp(idd_w1:idd_w1+lmax_dd2-1)=wsave(i,j,1:lmax_dd2)
            tmp(idd_phi1:idd_phi1+lmax_dd2-1)=phi(i,j,1:lmax_dd2)*bygrav
            tmp(idd_sr1:idd_sr1+lmax_dd2-1)=srnflb_save(i,j,1:lmax_dd2)
     *           *cosz1(i,j)
            tmp(idd_tr1:idd_tr1+lmax_dd2-1)=trnflb_save(i,j,1:lmax_dd2)

            DO n=1,Ntm_dust
              n1=n_clay+n-1

              tmp(idd_load1:idd_load1+lmax_dd2-1)
     *             =tmp(idd_load1:idd_load1+lmax_dd2-1)+trm(i,j
     *             ,1:lmax_dd2,n1)*byaxyp(i,j)
              tmp(idd_conc1:idd_conc1+lmax_dd2-1)
     *             =tmp(idd_conc1:idd_conc1+lmax_dd2-1)+trm(i,j
     *             ,1:lmax_dd2,n1)*byam(1,i,j)*byaxyp(i,j)
              tmp(idd_tau1:idd_tau1+lmax_dd2-1)=tmp(idd_tau1:idd_tau1
     *             +lmax_dd2-1)+ttausv_save(i,j,n1,1:lmax_dd2)
              tmp(idd_tau_cs1:idd_tau_cs1+lmax_dd2-1)
     *             =tmp(idd_tau_cs1:idd_tau_cs1+lmax_dd2-1)
     *             +ttausv_cs_save(i,j,n1,1:lmax_dd2)

            END DO

            ADIURN(idxd(:),kr,ih)=ADIURN(idxd(:),kr,ih)+tmp(idxd(:))
#ifndef NO_HDIURN
            HDIURN(idxd(:),kr,ihm)=HDIURN(idxd(:),kr,ihm)+tmp(idxd(:))
#endif

          END IF
#endif
        endif
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine ahourly

      module msu_wts_mod
      implicit none
      save
      integer, parameter :: nmsu=200 , ncols=4
      real*8 plbmsu(nmsu),wmsu(ncols,nmsu)
      contains
      subroutine read_msu_wts
      use filemanager
      integer n,l,iu_msu
c**** read in the MSU weights file
      call openunit('MSU_wts',iu_msu,.false.,.true.)
      do n=1,4
        read(iu_msu,*)
      end do
      do l=1,nmsu
        read(iu_msu,*) plbmsu(l),(wmsu(n,l),n=1,ncols)
      end do
      call closeunit(iu_msu)
      end subroutine read_msu_wts
      end module msu_wts_mod

      subroutine diag_msu(pland,ts,tlm,ple,tmsu2,tmsu3,tmsu4)
!@sum diag_msu computes MSU channel 2,3,4 temperatures as weighted means
!@auth Reto A Ruedy (input file created by Makiko Sato)
      USE MODEL_COM, only : lm
      use msu_wts_mod
      implicit none
      real*8, intent(in) :: pland,ts,tlm(lm),ple(lm+1)
      real*8, intent(out) :: tmsu2,tmsu3,tmsu4

      real*8 tlmsu(nmsu),tmsu(ncols)
      real*8 plb(0:lm+2),tlb(0:lm+2)
      integer l

c**** find edge temperatures (assume continuity and given means)
      tlb(0)=ts ; plb(0)=plbmsu(1) ; tlb(1)=ts
      plb(1:lm+1) = ple(1:lm+1)
      do l=1,lm
        tlb(l+1)=2*tlm(l)-tlb(l)
      end do
      tlb(lm+2)=tlb(lm+1) ; plb(lm+2)=0.
      call vntrp1 (lm+2,plb,tlb, nmsu-1,plbmsu,tlmsu)
c**** find MSU channel 2,3,4 temperatures
      tmsu(:)=0.
      do l=1,nmsu-1
        tmsu(:)=tmsu(:)+tlmsu(l)*wmsu(:,l)
      end do
      tmsu2 = (1-pland)*tmsu(1)+pland*tmsu(2)
      tmsu3 = tmsu(3)
      tmsu4 = tmsu(4)

      return
      end subroutine diag_msu

      SUBROUTINE init_DIAG(istart,num_acc_files)
!@sum  init_DIAG initializes the diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday,kapa,undef
      USE MODEL_COM, only : lm,Itime,ItimeI,Itime0,pmtop,nfiltr,jhour
     *     ,jdate,jmon,amon,jyear,jhour0,jdate0,jmon0,amon0,jyear0,idacc
     *     ,ioread_single,xlabel,iowrite_single,iyear1,nday,dtsrc,dt
     *     ,nmonav,ItimeE,lrunid,focean,pednl00,pmidl00,lm_req
      USE GEOM, only : axyp,imaxj,lon2d_dg,lat2d_dg
      USE GEOM, only : lonlat_to_ij
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE DIAG_COM, only : TSFREZ => TSFREZ_loc
      USE DIAG_COM, only : NPTS, NAMDD, NDIUPT, IJDD,LLDD, ISCCP_DIAGS
      USE DIAG_COM, only : monacc, acc_period, keyct, KEYNR, PLE
      USE DIAG_COM, only : PLM, p1000k, icon_AM, NOFM
      USE DIAG_COM, only : PLE_DN, icon_KE, NSUM_CON, IA_CON, SCALE_CON
      USE DIAG_COM, only : TITLE_CON, PSPEC, LSTR, NSPHER, KLAYER
      USE DIAG_COM, only : ISTRAT, kgz, pmb, kgz_max
      USE DIAG_COM, only : TF_DAY1, TF_LAST, TF_LKON, TF_LKOFF
      USE DIAG_COM, only : name_consrv, units_consrv, lname_consrv
      USE DIAG_COM, only : CONPT0, icon_MS, icon_TPE, icon_WM, icon_EWM
      USE DIAG_COM, only : nreg,jreg,titreg,namreg,sarea_reg
      USE diag_com,ONLY : adiurn_dust,adiurn_loc,areg_loc,aisccp_loc
#ifndef NO_HDIURN
     &     ,hdiurn_loc
#endif
      USE diag_com,only : lh_diags
      USE DIAG_LOC
      USE PARAM
      USE FILEMANAGER
      USE DOMAIN_DECOMP_ATM, only: GRID,GET,WRITE_PARALLEL,
     &     AM_I_ROOT,GLOBALSUM
      use msu_wts_mod
      IMPLICIT NONE
      integer, intent(in) :: istart,num_acc_files
      INTEGER I,J,L,K,KL,n,ioerr,months,years,mswitch,ldate
     *     ,jday0,jday,moff,kb,l850,l300,l50
      REAL*8 PLE_tmp
      CHARACTER CONPT(NPTS)*10
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER :: J_0,J_1, I_0,I_1
!@var out_line local variable to hold mixed-type output for parallel I/O
      character(len=300) :: out_line
      character(len=80) :: filenm
!@var iu_REG unit number for regions file
      INTEGER iu_REG
#ifdef CUBE_GRID
C***  regions defined as rectangles 
      integer, dimension(23) :: NRECT
      character*4, dimension(23,6) :: CORLON,CORLAT 
      real*8, dimension(23,6) :: DCORLON,DCORLAT   !lat-lon coordinates of rect. corners 
      integer :: ireg,irect,icorlon,icorlat
      real*8::lon,lat
#else
c dummy global array to read special-format regions file lacking a
c a parallelized i/o routine that understands it
      integer :: jreg_glob(im,jm)
#endif
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     area_part

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP


C****   READ SPECIAL REGIONS
#ifndef CUBE_GRID
      call openunit("REG",iu_REG,.true.,.true.)
      READ(iu_REG) TITREG,JREG_glob,NAMREG
      jreg(i_0:i_1,j_0:j_1) = jreg_glob(i_0:i_1,j_0:j_1)
#else
c**** Regions are defined in reg.txt input file as union of rectangles
c**** independant of resolution & grid type
        
      call openunit("REG",iu_REG,.false.,.true.)
      READ(iu_REG,'(A80)') TITREG !read title
      READ (iu_REG,'(I2)') (NRECT(I), I=1,23 ) !#of rectangles per region
      READ (iu_REG,'(A4,1X,A4)') (NAMREG(1,I),NAMREG(2,I),I=1,23) !Read region name
        
c**** Read cordinates of rectangles
c**** (NWcorner long, NW corner lat, SE corner long, SE corner lat)(1:23)
c**** 0555 = no rectangle
      READ (iu_REG,'(11(A4,1X),A4)') (      
     &       CORLON(I,1),CORLAT(I,1),
     &       CORLON(I,2),CORLAT(I,2),
     &       CORLON(I,3),CORLAT(I,3),
     &       CORLON(I,4),CORLAT(I,4),
     &       CORLON(I,5),CORLAT(I,5),
     &       CORLON(I,6),CORLAT(I,6), I=1,23)
 
c**** Convert to integer
      do i=1,23
        do j=1,6
          read(corlon(i,j),'(I4)') icorlon
          dcorlon(i,j) =icorlon
          read(corlat(i,j),'(I4)') icorlat
          dcorlat(i,j) =icorlat
        enddo
      enddo

c**** determine the region to which each cell belongs
      JREG(:,:)=24
      do j=j_0,j_1
        do i=i_0,i_1
          lon=lon2d_dg(i,j)
          lat=lat2d_dg(i,j)
          do ireg=1,23
            do irect=1,NRECT(ireg)
              if ( lat > dcorlat(ireg,2*irect-1) .or.
     &             lat < dcorlat(ireg,2*irect  ) )  cycle
              if(dcorlon(ireg,2*irect-1)<dcorlon(ireg,2*irect)) then
                if( lon < dcorlon(ireg,2*irect-1) .or.
     &              lon > dcorlon(ireg,2*irect  ) ) cycle
              else ! wraparound
                if( lon < dcorlon(ireg,2*irect-1) .and.
     &              lon > dcorlon(ireg,2*irect  ) ) cycle
              endif
              JREG(i,j)=ireg
            enddo
          enddo
        enddo
      enddo
#endif
      IF (AM_I_ROOT()) then
        WRITE(6,*) ' read REGIONS from unit ',iu_REG,': ',TITREG
      endif 
      call closeunit(iu_REG)
c
c calculate the areas of the special regions
c
      do n=1,nreg
        where(jreg(i_0:i_1,j_0:j_1).eq.n)
          area_part(i_0:i_1,j_0:j_1) = axyp(i_0:i_1,j_0:j_1)
        elsewhere
          area_part(i_0:i_1,j_0:j_1) = 0d0
        end where
        call globalsum(grid,area_part,sarea_reg(n),all=.true.)
      enddo


C**** Initialse diurnal diagnostic locations (taken from the 4x5 res)
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      NAMDD =
     &   (/'AUSD', 'MWST', 'SAHL', 'EPAC', 'AF01',
     &     'AF02', 'AF03', 'AF04', 'AF05', 'ASA1',
     &     'AF06', 'AME1', 'AF07', 'AF08', 'AF09',
     &     'ARAB', 'AF10', 'AF11', 'ASA2', 'AUS1',
     &     'AUS2', 'AUS3', 'AF12', 'AF13', 'ASA3',
     &     'AF14', 'ASA4', 'AUS4', 'AME2', 'AF15',
     &     'AME3', 'AF16', 'AF17', 'AF18' /)
      LLDD = RESHAPE( (/
     &   132.5,-26.,  -97.5, 42.,    2.5, 14.,
     &  -117.5, -2.,   32.5, 18.,   22.5, 30.,
     &    -2.5, 26.,  -12.5, 26.,   12.5, 18.,
     &    62.5, 38.,   17.5, 18.,  -67.5,-42.,
     &    -7.5, 26.,   22.5, 26.,   27.5, 30.,
     &    42.5, 34.,   12.5, 14.,   -2.5, 22.,
     &    62.5, 42.,  137.5,-30.,  137.5,-26.,
     &   127.5,-30.,   22.5, 18.,   -2.5, 18.,
     &   102.5, 42.,    7.5, 30.,   72.5, 26.,
     &   132.5,-30.,  -67.5,-38.,    2.5, 26.,
     &   -67.5,-46.,   22.5, 22.,    2.5, 22.,
     &    17.5, 30.  /),(/2,NDIUPT/))
#else
c defaults for diurnal diagnostics
      NAMDD = (/ 'AUSD', 'MWST', 'SAHL', 'EPAC' /)      
      LLDD = RESHAPE( (/  
     &      132.5, -26.,
     &      -97.5,  42.,
     &        2.5,  14.,
     &     -117.5,  -2.
     &     /),(/2,4/))
#endif
      call sync_param( "LLDD", LLDD(1:2,1), 2*NDIUPT )
      do n=1,ndiupt
        call lonlat_to_ij(lldd(1,n),ijdd(1,n))
      enddo

      call sync_param( "NAMDD", NAMDD, NDIUPT )
c if people still want to specify dd points as ij, let them
      call sync_param( "IJDD", IJDD(1:2,1), 2*NDIUPT )

      call sync_param( "isccp_diags",isccp_diags)
      call sync_param( "adiurn_dust",adiurn_dust)
      call sync_param( "lh_diags",lh_diags)

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
          filenm=acc_period(1:Ldate)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
        end if
        ItimeE = -1
        close (6)
        open(6,file=acc_period(1:Ldate)//'.'//XLABEL(1:LRUNID)//'.PRT',
     *       FORM='FORMATTED')
      END IF

C**** Initialize certain arrays used by more than one print routine
      DO L=1,LM
        PLE(L)   =pednl00(l+1)
        PLE_DN(L)=pednl00(l)
        PLM(L)   =pmidl00(l)
      END DO
      PLM(LM+1:LM+LM_REQ)=pmidl00(lm+1:lm+lm_req)

      p1000k=1000.0**kapa

C**** Initialise some local constants (replaces IFIRST constructions)
C**** From DIAGA:
      DO L=1,LM
        LUPA(L)=L+1
        LDNA(L)=L-1
      END DO
      LDNA(1)=1
      LUPA(LM)=LM

C**** From DIAGB (PM, PMO are fixed, PL,PLO will vary)
      PM(1)=1200.   ! ensures below surface for extrapolation
      DO L=2,LM+1
        PL(L)=pednl00(l)
        PM(L)=pednl00(l)
      END DO
      DO L=1,LM
        PLO(L)=pmidl00(l)
        PMO(L)=.5*(PM(L)+PM(L+1))
      END DO

C**** From DIAG7A
      L850=LM
      L300=LM
      L50=LM
      DO L=LM-1,1,-1
        PLE_tmp=.25*(PEDNL00(L)+2.*PEDNL00(L+1)+PEDNL00(L+2))
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
        IF (PEDNL00(L+1)+1d-5.lt.PSPEC(KL) .and.
     *      PEDNL00(L)  +1d-5.gt.PSPEC(KL)) THEN
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
          DO I=I_0,IMAXJ(J)
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

c
c zero out certain non-distributed arrays
c
      areg_loc = 0
      aisccp_loc = 0
      adiurn_loc = 0
#ifndef NO_HDIURN
      hdiurn_loc = 0
#endif

c
c read MSU weighting functions for diagnostics
c
      call read_msu_wts

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
      USE DOMAIN_DECOMP_ATM, only: grid,am_i_root
      IMPLICIT NONE
      INTEGER :: isum !@var isum if =1 preparation to add up acc-files
      INTEGER jd0

      IDACC(1:12)=0
      if (kradia.gt.0) then
        AFLX_ST = 0.
        if (isum.eq.1) return
        go to 100
      end if

      if(am_i_root()) then
        aj = 0; ajl = 0; asjl = 0; agc = 0; consrv = 0
      endif
      AJ_loc=0    ; AREG_loc=0; AREG=0
      AJL_loc=0  ; ASJL_loc=0   ; AIJ_loc=0
      AIJL_loc=0   ; ENERGY=0 ; CONSRV_loc=0
      SPECA=0 ; ATPE=0 ; WAVE=0 ; AGC_loc=0   ; AIJK_loc=0
#ifndef NO_HDIURN
      HDIURN=0; HDIURN_loc=0
#endif
      ADIURN=0 ; ADIURN_loc=0; AISCCP=0; AISCCP_loc=0
#ifdef TRACERS_ON
      call reset_trdiag
#endif
      call reset_ODIAG(isum)  ! ocean diags if required
      call reset_icdiag       ! ice dynamic diags if required

      if (isum.eq.1) return ! just adding up acc-files

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
      USE MODEL_COM, only : im,jm,jday,focean
      USE GEOM, only : imaxj,lat2d
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE GHY_COM, only : fearth
      USE DIAG_COM, only : aij=>aij_loc
     *     ,ij_lkon,ij_lkoff,ij_lkice,tsfrez=>tsfrez_loc,tdiurn
     *     ,tf_lkon,tf_lkoff,tf_day1,tf_last
      USE DOMAIN_DECOMP_ATM, only : GRID,GET,am_i_root
#ifdef HTAP_LIKE_DIAGS
      USE RAD_COM, only:  ttausv_sum,ttausv_count
#endif
      IMPLICIT NONE
      INTEGER I,J
      INTEGER :: J_0, J_1, I_0,I_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      IF (JDAY.EQ.32) THEN
        DO J=J_0,J_1
        DO I=I_0,I_1
          if(lat2d(i,j).gt.0.) then
            TSFREZ(I,J,TF_DAY1)=JDAY
          else
            TSFREZ(I,J,TF_LAST)=JDAY
          endif
        ENDDO
        ENDDO
      ELSEIF (JDAY.EQ.213) THEN
        DO J=J_0,J_1
        DO I=I_0,I_1
          if(lat2d(i,j).lt.0.) then
            TSFREZ(I,J,TF_DAY1)=JDAY
          endif
        ENDDO
        ENDDO
      END IF
C**** set and initiallise freezing diagnostics
C**** Note that TSFREZ saves the last day of no-ice and some-ice.
C**** The AIJ diagnostics are set once a year (zero otherwise)
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          if(lat2d(i,j).lt.0.) then
C**** initialize/save South. Hemi. on Feb 28
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
         DO I=I_0,I_1
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
#ifdef HTAP_LIKE_DIAGS
            ttausv_sum(I,J,:)=0.d0 
#endif
         END DO
      END DO
#ifdef HTAP_LIKE_DIAGS
      ttausv_count=0.d0
#endif
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
      USE DOMAIN_DECOMP_ATM, only : WRITE_PARALLEL
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
      USE MODEL_COM, only : im,jm,focean,flice,itocean
     *     ,itoice,itlandi,itearth,itlake,itlkice,ftype
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : flake
      USE GHY_COM, only : fearth
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
      INTEGER I,J
      INTEGER :: J_0,J_1,I_0,I_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
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

      subroutine calc_derived_aij
      USE CONSTANT, only : grav,rgas,bygrav,tf,teeny
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE MODEL_COM, only : idacc,zatmo,fearth0,flice,focean,lm,pmtop
      USE GEOM, only : imaxj
      USE DIAG_COM, only : aij=>aij_loc,tsfrez=>tsfrez_loc,
     &  aijl=>aijl_loc,ia_ij,ia_src,ia_inst,ia_dga,tf_last,tf_day1,
     *  ij_topo, ij_wsmn, ij_wsdir, ij_jet, ij_jetdir, ij_grow,
     *  ij_netrdp, ij_albp, ij_albg, ij_albv,
     *  ij_fland, ij_dzt1, ij_albgv, ij_clrsky, ij_pocean, ij_ts,
     *  ij_RTSE, ij_HWV, ij_PVS,
     &  IJ_TRNFP0,IJ_SRNFP0,IJ_TRSUP,IJ_TRSDN,IJ_EVAP,IJ_QS,IJ_PRES,
     &  IJ_SRREF,IJ_SRVIS,IJ_SRINCP0,IJ_SRINCG,IJ_SRNFG,IJ_PHI1K,
     &  IJ_US,IJ_VS,IJ_UJET,IJ_VJET,IJ_CLDCV,IJ_TATM,IJK_DP,IJK_TX,
     &  IJ_MSU2,IJ_MSU3,IJ_MSU4,
     &  KGZ_MAX,GHT,PMB
      IMPLICIT NONE
      INTEGER :: I,J,L,K,K1,K2
      INTEGER :: J_0,J_1,I_0,I_1
      REAL*8 :: SCALEK
      real*8 :: ts,pland,tlm(lm),ple(lm+1),dp,tmsu2,tmsu3,tmsu4

      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      J_0 = GRID%J_STRT
      J_1 = GRID%J_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)

        k = ij_topo
        aij(i,j,k) = zatmo(i,j)*bygrav*idacc(ia_inst)

        k = ij_fland
        aij(i,j,k) = (fearth0(i,j)+flice(i,j))*idacc(ia_src)

        k = ij_netrdp
        aij(i,j,k) = aij(i,j,ij_trnfp0)+aij(i,j,ij_srnfp0)

        k = ij_grow
        aij(i,j,k) = (tsfrez(i,j,tf_last)-tsfrez(i,j,tf_day1))
     &       *idacc(ia_inst)

        k = ij_rtse
        aij(i,j,k) = aij(i,j,ij_trsup) - aij(i,j,ij_trsdn)

        k = ij_hwv
        aij(i,j,k) = aij(i,j,ij_evap)

        k = ij_pvs
        aij(i,j,k) = aij(i,j,ij_qs)*
     &       aij(i,j,ij_pres)/idacc(ia_ij(ij_pres))

        k = ij_albv
        aij(i,j,k) = aij(i,j,ij_srref)

        k = ij_albgv
        aij(i,j,k) = aij(i,j,ij_srvis)

        k = ij_albp
        aij(i,j,k) = aij(i,j,ij_srincp0)-aij(i,j,ij_srnfp0)*
     &       idacc(ia_ij(ij_srincp0))/idacc(ia_ij(ij_srnfp0))

        k = ij_albg
        aij(i,j,k) = aij(i,j,ij_srincg)-aij(i,j,ij_srnfg)*
     &       idacc(ia_ij(ij_srincg))/idacc(ia_ij(ij_srnfg))

        do k=ij_dzt1,ij_dzt1+kgz_max-2
          k1 = k-ij_dzt1+1  ; k2 = ij_phi1k + k1
          scalek = 1./(rgas*log(pmb(k1)/pmb(k1+1)))
          aij(i,j,k) = (scalek*(ght(k1+1)-ght(k1))*grav-tf)*
     &         idacc(ia_ij(ij_phi1k))
     &         +scalek*(aij(i,j,k2)-aij(i,j,k2-1))
        enddo

        k = ij_jet
        aij(i,j,k) = sqrt(aij(i,j,ij_ujet)**2+aij(i,j,ij_vjet)**2)

        k = ij_wsmn
        aij(i,j,k) = sqrt(aij(i,j,ij_us)**2+aij(i,j,ij_vs)**2)

        k = ij_wsdir
        aij(i,j,k) = atan2(aij(i,j,ij_us)+teeny,aij(i,j,ij_vs)+teeny)
     &       *idacc(ia_inst)

        k = ij_jetdir
        aij(i,j,k)=atan2(aij(i,j,ij_ujet)+teeny,aij(i,j,ij_vjet)+teeny)
     &       *idacc(ia_inst)

        k = ij_clrsky
        aij(i,j,k) = idacc(ia_ij(ij_cldcv))-aij(i,j,ij_cldcv)

        k = ij_pocean
        aij(i,j,k) = idacc(ia_ij(k))*focean(i,j)

        k = ij_tatm
        aij(i,j,k) = sum(aijl(i,j,:,ijk_tx))

C**** Find MSU channel 2,3,4 temperatures (simple lin.comb. of Temps)
        pland = fearth0(i,j)+flice(i,j)
        ts = aij(i,j,ij_ts)/idacc(ia_ij(ij_ts))
        ple(lm+1) = pmtop
        do l=lm,1,-1
          dp = aijl(i,j,l,ijk_dp)
          ple(l) = ple(l+1)+dp/idacc(ia_dga)
          tlm(l) = ts
          if(dp.gt.0.) tlm(l)=aijl(i,j,l,ijk_tx)/dp
        enddo
        call diag_msu(pland,ts,tlm,ple,tmsu2,tmsu3,tmsu4)
        aij(i,j,ij_msu2) = tmsu2*idacc(ia_inst)
        aij(i,j,ij_msu3) = tmsu3*idacc(ia_inst)
        aij(i,j,ij_msu4) = tmsu4*idacc(ia_inst)

      ENDDO
      ENDDO

      return
      end subroutine calc_derived_aij

      SUBROUTINE DIAGJ_PREP
      USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : dtsrc,idacc,ntype
      USE DIAG_COM, only : jm_budg,
     &     aj,ntype_out,wt=>wtj_comp,aj_out,areg,areg_out,
     &     nreg,kaj,j_albp0,j_srincp0,j_albg,j_srincg,
     &     j_srabs,j_srnfp0,j_srnfg,j_trnfp0,j_hsurf,j_trhdt,j_trnfp1,
     *     j_hatm,j_rnfp0,j_rnfp1,j_srnfp1,j_rhdt,j_hz1,j_prcp,j_prcpss,
     *     j_prcpmc,j_hz0,j_implh,j_shdt,j_evhdt,j_eprcp,j_erun,
     *     j_hz2,j_ervr,
     *     ia_src,ia_rad,
     &     sarea=>sarea_reg,
     &     hemis_j,dxyp_budg
      IMPLICIT NONE
      REAL*8 :: A1BYA2,hemfac
      INTEGER :: J,JR,J1,J2,K,M,IT

      IF(.NOT.AM_I_ROOT()) RETURN

C**** CALCULATE THE DERIVED QUANTTIES
      A1BYA2=IDACC(ia_src)/(IDACC(ia_rad)+teeny)
      DO JR=1,23
        AREG(JR,J_SRABS) =AREG(JR,J_SRNFP0)-AREG(JR,J_SRNFG)
        AREG(JR,J_RNFP0) =AREG(JR,J_SRNFP0)+AREG(JR,J_TRNFP0)
        AREG(JR,J_RNFP1) =AREG(JR,J_SRNFP1)+AREG(JR,J_TRNFP1)
        AREG(JR,J_ALBP0)=AREG(JR,J_SRINCP0)-AREG(JR,J_SRNFP0)
        AREG(JR,J_ALBG)=AREG(JR,J_SRINCG)-AREG(JR,J_SRNFG)
        AREG(JR,J_RHDT)  =A1BYA2*AREG(JR,J_SRNFG)*DTSRC+AREG(JR,J_TRHDT)
        AREG(JR,J_PRCP)  =AREG(JR,J_PRCPSS)+AREG(JR,J_PRCPMC)
        AREG(JR,J_HZ0)=AREG(JR,J_RHDT)+AREG(JR,J_SHDT)+
     *                 AREG(JR,J_EVHDT)+AREG(JR,J_EPRCP)
        AREG(JR,J_HZ1)=AREG(JR,J_HZ0)+AREG(JR,J_ERVR)
        AREG(JR,J_HZ2)=AREG(JR,J_HZ1)-AREG(JR,J_ERUN)-AREG(JR,J_IMPLH)
      END DO
      DO J=1,JM_BUDG
      DO IT=1,NTYPE
        AJ(J,J_SRABS ,IT)=AJ(J,J_SRNFP0,IT)-AJ(J,J_SRNFG,IT)
        AJ(J,J_RNFP0 ,IT)=AJ(J,J_SRNFP0,IT)+AJ(J,J_TRNFP0,IT)
        AJ(J,J_RNFP1 ,IT)=AJ(J,J_SRNFP1,IT)+AJ(J,J_TRNFP1,IT)
        AJ(J,J_ALBP0,IT)=AJ(J,J_SRINCP0,IT)-AJ(J,J_SRNFP0,IT)
        AJ(J,J_ALBG,IT)=AJ(J,J_SRINCG,IT)-AJ(J,J_SRNFG,IT)
        AJ(J,J_RHDT  ,IT)=A1BYA2*AJ(J,J_SRNFG,IT)*DTSRC+AJ(J,J_TRHDT,IT)
        AJ(J,J_PRCP  ,IT)=AJ(J,J_PRCPSS,IT)+AJ(J,J_PRCPMC,IT)
        AJ(J,J_HZ0,IT)=AJ(J,J_RHDT,IT)+AJ(J,J_SHDT,IT)+
     *                 AJ(J,J_EVHDT,IT)+AJ(J,J_EPRCP,IT)
        AJ(J,J_HZ1,IT)=AJ(J,J_HZ0,IT)+AJ(J,J_ERVR,IT)
        AJ(J,J_HZ2,IT)=AJ(J,J_HZ1,IT)-AJ(J,J_ERUN,IT)-AJ(J,J_IMPLH,IT)
      END DO
      END DO

c
c construct the composites of surface types
c
      aj_out = 0d0
      do m=1,ntype_out
        do it=1,ntype
          aj_out(:,:,m) = aj_out(:,:,m) + aj(:,:,it)*wt(m,it)
        enddo
      enddo

c
c compute hemispheric and global means
c
      hemfac = 2./sum(dxyp_budg)
      do m=1,ntype_out
      do k=1,kaj
        j1 = 1; j2 = jm_budg/2
        hemis_j(1,k,m) = hemfac*sum(aj_out(j1:j2,k,m)*dxyp_budg(j1:j2))
        j1 = jm_budg/2+1; j2 = jm_budg
        hemis_j(2,k,m) = hemfac*sum(aj_out(j1:j2,k,m)*dxyp_budg(j1:j2))
        hemis_j(3,k,m) = .5*(hemis_j(1,k,m)+hemis_j(2,k,m))
      enddo
      enddo

c
c scale areg by area
c
      do jr=1,nreg
        areg_out(jr,:) = areg(jr,:)/sarea(jr)
      enddo

      RETURN
      END SUBROUTINE DIAGJ_PREP

      subroutine diagjl_prep
      use model_com, only : lm,lm_req
      use domain_decomp_atm, only : am_i_root
      use diag_com, only : kajl,jm_budg,
     &     ajl,asjl,jl_srhr,jl_trcr,jl_rad_cool,
     &     dxyp_budg,hemis_jl,vmean_jl
      implicit none
      integer :: j,j1,j2,l,k,lr
      real*8 :: hemfac

      if(.not.am_i_root()) return

      do j=1,jm_budg
        do lr=1,lm_req
          asjl(j,lr,5)=asjl(j,lr,3)+asjl(j,lr,4)
        enddo
        do l=1,lm
          ajl(j,l,jl_rad_cool)=ajl(j,l,jl_srhr)+ajl(j,l,jl_trcr)
        enddo
      enddo

c
c compute hemispheric/global means and vertical sums
c
      hemfac = 2./sum(dxyp_budg)
      do k=1,kajl
        do l=1,lm
          j1 = 1; j2 = jm_budg/2
          hemis_jl(1,l,k) = hemfac*sum(ajl(j1:j2,l,k)*dxyp_budg(j1:j2))
          j1 = jm_budg/2+1; j2 = jm_budg
          hemis_jl(2,l,k) = hemfac*sum(ajl(j1:j2,l,k)*dxyp_budg(j1:j2))
          hemis_jl(3,l,k) = .5*(hemis_jl(1,l,k)+hemis_jl(2,l,k))
        enddo
        vmean_jl(1:jm_budg,1,k) = sum(ajl(:,:,k),dim=2)
        vmean_jl(jm_budg+1:jm_budg+3,1,k) = sum(hemis_jl(:,:,k),dim=2)
      enddo

      return
      end subroutine diagjl_prep

      subroutine diag_isccp_prep
c calculates the denominator array for ISCCP histograms
      use diag_com, only : aij=>aij_loc,ij_scldi,nisccp,wisccp
      use clouds_com, only : isccp_reg2d
      use domain_decomp_atm, only : grid,sumxpe
      use geom, only : axyp
      implicit none
      real*8 wisccp_loc(nisccp)
      integer :: i,j,n,j_0,j_1,i_0,i_1
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      wisccp_loc(:) = 0.
      do j=j_0,j_1
      do i=i_0,i_1
        n=isccp_reg2d(i,j)
        if(n.gt.0)
     &       wisccp_loc(n)=wisccp_loc(n)+aij(i,j,ij_scldi)*axyp(i,j)
      enddo
      enddo
      call sumxpe(wisccp_loc,wisccp)
      return
      end subroutine diag_isccp_prep

      SUBROUTINE VNTRP1 (KM,P,AIN,  LMA,PE,AOUT)
C**** Vertically interpolates a 1-D array
C**** Input:       KM = number of input pressure levels
C****            P(K) = input pressure levels (mb)
C****          AIN(K) = input quantity at level P(K)
C****             LMA = number of vertical layers of output grid
C****           PE(L) = output pressure levels (mb) (edges of layers)
C**** Output: AOUT(L) = output quantity: mean between PE(L-1) & PE(L)
C****
      implicit none
      integer, intent(in) :: km,lma
      REAL*8, intent(in)  :: P(0:KM),AIN(0:KM),    PE(0:LMA)
      REAL*8, intent(out) :: AOUT(LMA)

      integer k,k1,l
      real*8 pdn,adn,pup,aup,psum,asum

C****
      PDN = PE(0)
      ADN = AIN(0)
      K=1
C**** Ignore input levels below ground level pe(0)=p(0)
      IF(P(1).GT.PE(0)) THEN
         DO K1=2,KM
         K=K1
         IF(P(K).LT.PE(0)) THEN  ! interpolate to ground level
           ADN=AIN(K)+(AIN(K-1)-AIN(K))*(PDN-P(K))/(P(K-1)-P(K))
           GO TO 300
         END IF
         END DO
         STOP 'VNTRP1 - error - should not get here'
      END IF
C**** Integrate - connecting input data by straight lines
  300 DO 330 L=1,LMA
      ASUM = 0.
      PSUM = 0.
      PUP = PE(L)
  310 IF(P(K).le.PUP)  GO TO 320
      PSUM = PSUM + (PDN-P(K))
      ASUM = ASUM + (PDN-P(K))*(ADN+AIN(K))/2.
      PDN  = P(K)
      ADN  = AIN(K)
      K=K+1
      IF(K.LE.KM) GO TO 310
      stop 'VNTRP1 - should not happen'
C****
  320 AUP  = AIN(K) + (ADN-AIN(K))*(PUP-P(K))/(PDN-P(K))
      PSUM = PSUM + (PDN-PUP)
      ASUM = ASUM + (PDN-PUP)*(ADN+AUP)/2.
      AOUT(L) = ASUM/PSUM
      PDN = PUP
  330 ADN = AUP
C****
      RETURN
      END subroutine vntrp1

      subroutine calc_derived_acc_atm
      use diag_com, only : isccp_diags
      implicit none
      call gather_zonal_diags
      call collect_scalars
      call calc_derived_aij
      call diagj_prep
      call diagjl_prep
      call diaggc_prep
      call diag_river_prep
      if(isccp_diags.eq.1) call diag_isccp_prep
      return
      end subroutine calc_derived_acc_atm
