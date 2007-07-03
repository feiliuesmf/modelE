#include "rundeck_opts.h"

!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont init_LI,PRECIP_LI,GROUND_LI,daily_LI

      SUBROUTINE init_LI(istart)
!@sum  init_ice initialises landice arrays
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : edpery,sday,lhm
      USE MODEL_COM, only : im,jm,flice,focean,dtsrc
      USE GEOM, only : dxyp,imaxj
      USE LANDICE, only: ace1li,ace2li,glmelt_on,glmelt_fac_nh
     *     ,glmelt_fac_sh,fwarea_sh,fwarea_nh,accpda,accpdg,eaccpda
     *     ,eaccpdg
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,traccpda,traccpdg
#endif
#endif               /* TNL: inserted */
      USE LANDICE_COM, only : tlandi,snowli,loc_gla,loc_glg,mdwnimp
     *     ,edwnimp
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trli0,trdwnimp
#endif
      USE FLUXES, only : gtemp,gmelt,egmelt
#ifdef TRACERS_WATER
     *     ,gtracer
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif
#endif
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only :  trglac
#endif
      USE DIAG_COM, only : npts,icon_MLI,icon_HLI,title_con,conpt0
      USE PARAM
      USE DOMAIN_DECOMP, only : GRID,GET, GLOBALSUM
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER, INTENT(IN) :: istart

      INTEGER, PARAMETER :: NBOXMAX=18
      INTEGER IFW(NBOXMAX),JFW(NBOXMAX)
      INTEGER :: JML, JMU, IML1, IMU1, IML2, IMU2, NBOX
      REAL*8 FAC_SH,FAC_NH
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO)::FWAREA_part
      LOGICAL :: do_glmelt
      INTEGER I,J,N
      INTEGER :: J_0,J_1

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
C**** set GTEMP array for landice
      DO J=J_0,J_1
        DO I=1,IM
          IF (FLICE(I,J).gt.0) THEN
            GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
#ifdef TRACERS_WATER
            if (istart.ge.9) then
            IF (SNOWLI(I,J).gt.1d-5) THEN
              GTRACER(:,3,I,J)=TRSNOWLI(:,I,J)/SNOWLI(I,J)
            ELSE
              GTRACER(:,3,I,J)=TRLNDI(:,I,J)/(ACE1LI+ACE2LI)
            END IF
            end if
#endif
          END IF
        END DO
      END DO

C**** Calculate (fixed) iceberg melt terms from Antarctica and Greenland
      call sync_param("glmelt_on",glmelt_on)
      call sync_param("glmelt_fac_nh",glmelt_fac_nh)
      call sync_param("glmelt_fac_sh",glmelt_fac_sh)

C**** Note these parameters are highly resolution dependent!
C**** Around Antarctica, fresh water is added from 78S to 62S
C**** and from 0-60W and 135W to 165E
C**** Around Greenland the freshwater is added to both the east and
C**** west coasts in the one grid box closest to the coast which is
C**** determined by the direction in the river flow file.
C**** This information could be read in from a file.
      IF (JM.eq.46) THEN        ! 4x5
        JML=4 ; JMU=8 ; IML1=24 ; IMU1=36
        IML2=69 ; IMU2=11 ; NBOX=10
        IFW(1:NBOX) = (/26,25,25,25,29,29,30,31,32,33/)
        JFW(1:NBOX) = (/39,40,41,42,39,40,40,40,41,42/)
        do_glmelt=.true.
      ELSEIF (JM.eq.90) THEN    ! 2x2.5
        JML=7 ; JMU=11 ; IML1=48 ; IMU1=69
        IML2=136 ; IMU2=18 ; NBOX=18
        IFW(1:NBOX) = (/50,50,51,51,51,52,53,56,56,57,58,59,60,61,62,63
     *       ,64,64/)
        JFW(1:NBOX) = (/82,81,80,79,78,77,76,76,77,78,78,78,79,79,79,80
     *       ,81,82/)
        do_glmelt=.true.
      ELSEIF (JM.eq.24) THEN    ! 8x10 (do nothing)
        do_glmelt=.false.
      ELSE  ! unknown resolution
        WRITE(*,*) "Glacial Melt flux: unknown resolution JM=",JM
        WRITE(*,*) "Please edit init_LI if glacial melt is required."
        do_glmelt=.false.
      END IF

C**** integrate area (which will depend on resolution/landmask)
C**** Antarctica melt area
      FWAREA_part=0.
      DO J=MAX(J_0,JML),MIN(J_1,JMU)
        DO I=1,IM
          IF (FOCEAN(I,J).GT.0.) THEN
            IF ((I.GE.IML1.AND.I.LE.IMU1) .or. I.GE.IML2 .or. I.LE
     *           .IMU2) THEN
              LOC_GLA(I,J)=.TRUE.
              FWAREA_part(J)=FWAREA_part(J)+DXYP(J)*FOCEAN(I,J)
            END IF
          END IF
        END DO
      END DO
      CALL GLOBALSUM(grid, FWAREA_part, FWAREA_SH, all=.true.)

C**** Greenland melt area
      FWAREA_part=0.
      DO N=1,NBOX
        I=IFW(N)
        J=JFW(N)
        If (J >= J_0 .and. J <= J_1) THEN
          LOC_GLG(I,J)=.TRUE.
          FWAREA_part(J)=FWAREA_part(J)+DXYP(J)*FOCEAN(I,J)
        END IF
      END DO
      CALL GLOBALSUM(grid, FWAREA_part, FWAREA_NH, all=.true.)

C**** Intialise gmelt fluxes
      GMELT = 0. ; EGMELT = 0.
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
      TRGMELT = 0.
#endif
#endif /* TNL: inserted */

      if (do_glmelt .and. glmelt_on > 0) then
C****  Note that water goes in with as if it is ice at 0 deg.
C****  Possibly this should be a function of in-situ freezing temp?
C****
C**** gmelt_fac_nh/sh should be used to match accumulation in control
C**** run so that global mean sea level is nearly constant. We can't
C**** get perfect values, but the drift should hopefully be less
C**** than 1 mm/year. This will minimise shock once interactive values
C**** are calculated.
C****
C**** Initial average mass fluxes for Antarctica/Greenland

      if (istart.lt.8 .or. ACCPDA.eq.0) then
C**** Initiallise total mass/energy fluxes (only at start of run)
C**** The net accumulation from IPCC2 report is 2016x10**12 kg/year
C**** for Antarctica and for Greenland it is 316x10**12 kg/year
        ACCPDA = glmelt_fac_sh*2016d12 ! kg/year
        ACCPDG = glmelt_fac_nh*316d12 ! kg/year
        EACCPDA = -LHM*ACCPDA ; EACCPDG = -LHM*ACCPDG ! J/year
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
        TRACCPDA(:) = trglac(:)*ACCPDA ! kg/year
        TRACCPDG(:) = trglac(:)*ACCPDG ! kg/year
#endif
#endif /* TNL: inserted */

C**** initiallise implicit accumulators (note that these arrays are
C**** not used until at least one full year has passed)
        MDWNIMP=0.
        EDWNIMP=0.
#ifdef TRACERS_WATER
        TRDWNIMP=0.
#endif
      end if

! accumulation (kg per source time step) per water column
      FAC_SH=DTsrc/(EDPERY*SDAY*FWAREA_SH)
      FAC_NH=DTsrc/(EDPERY*SDAY*FWAREA_NH)

C**** Set GL MELT arrays
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF (LOC_GLA(I,J)) THEN
             GMELT(I,J) =  ACCPDA*FAC_SH*DXYP(J)*FOCEAN(I,J) ! kg
            EGMELT(I,J) = EACCPDA*FAC_SH*DXYP(J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            TRGMELT(:,I,J)= TRACCPDA(:)*FAC_SH*DXYP(J)*FOCEAN(I,J)  ! kg
#endif
#endif /* TNL: inserted */
          END IF
          IF (LOC_GLG(I,J)) THEN
             GMELT(I,J) =  ACCPDG*FAC_NH*DXYP(J)*FOCEAN(I,J) ! kg
            EGMELT(I,J) = EACCPDG*FAC_NH*DXYP(J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            TRGMELT(:,I,J) = TRACCPDG(:)*FAC_NH*DXYP(J)*FOCEAN(I,J) ! kg
#endif
#endif /* TNL: inserted */
          END IF
        END DO
      END DO

      end if

C**** Set conservation diagnostics for land ice mass, energy
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LAND ICE","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_MLI)
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LNDI ENR","(10**6 J/M^2)   ",
     *     "(10**-3 W/M^2)  ",1d-6,1d3,icon_HLI)
C****
      END SUBROUTINE init_LI

      SUBROUTINE PRECIP_LI
!@sum  PRECIP_LI driver for applying precipitation to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:PRECLI
      USE MODEL_COM, only : im,jm,flice,itlandi
      USE GEOM, only : imaxj,dxyp,bydxyp
      USE FLUXES, only : runoli,prec,eprec,gtemp
#ifdef TRACERS_WATER
     *     ,trunoli,trprec,gtracer
#endif
      USE LANDICE, only: ace1li,ace2li,precli
      USE LANDICE_COM, only : snowli,tlandi,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,ntm,trli0,trdwnimp
#endif
      USE DIAG_COM, only : aj=>aj_loc,aregj=>aregj_loc,aij=>aij_loc
     *     ,jreg,ij_f0li,ij_f1li,ij_erun2
     *     ,ij_runli,j_run,j_implh,j_implm
      USE DOMAIN_DECOMP, only : GRID,GET
      USE DOMAIN_DECOMP, only : GLOBALSUM, CHECKSUM, CHECKSUM_COLUMN
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0,PLICE,DXYPJ
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
C**** Get useful grid parameters
      INTEGER :: I, J, JR
      INTEGER :: J_0, J_1, J_0H, J_1H

      CALL GET(GRID,J_STRT=J_0      , J_STOP=J_1      ,
     &              J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )

      DO J=J_0,J_1
      DXYPJ=DXYP(J)
      DO I=1,IMAXJ(J)
      PLICE=FLICE(I,J)
      PRCP=PREC(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0
#ifdef TRACERS_WATER
      TRUNOLI(:,I,J)=0.
#endif
      IF (PLICE.gt.0 .and. PRCP.gt.0) THEN

        ENRGP=EPREC(I,J)      ! energy of precipitation
        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J)
        TRSNOW(:)=TRSNOWLI(:,I,J)
        TRPRCP(:)=TRPREC(:,I,J)*BYDXYP(J)
#endif
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP*PLICE

        CALL PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *       TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *       EDIFS,DIFS,ERUN2,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J)  =RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
C**** accumulate implicit fluxes for setting ocean balance
        MDWNIMP(I,J)=MDWNIMP(I,J)+DIFS *PLICE*DXYPJ
        EDWNIMP(I,J)=EDWNIMP(I,J)+ERUN2*PLICE*DXYPJ
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
        TRDWNIMP(:,I,J)=TRDWNIMP(:,I,J)+TRDIFS(:)*PLICE*DXYPJ
        IF (SNOW.gt.1d-5) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
c       AJ(J,J_TYPE, ITLANDI)=AJ(J,J_TYPE, ITLANDI)+      PLICE
        AJ(J,J_RUN,  ITLANDI)=AJ(J,J_RUN,  ITLANDI)+RUN0 *PLICE
C       AJ(J,J_ERUN ,ITLANDI)=AJ(J,J_ERUN ,ITLANDI)+ERUN0*PLICE ! (Tg=0)
        AJ(J,J_IMPLM,ITLANDI)=AJ(J,J_IMPLM,ITLANDI)+DIFS *PLICE
        AJ(J,J_IMPLH,ITLANDI)=AJ(J,J_IMPLH,ITLANDI)+ERUN2*PLICE
        AREGJ(JR,J,J_RUN)  =AREGJ(JR,J,J_RUN)  +RUN0 *PLICE*DXYPJ
c       AREGJ(JR,J,J_ERUN) =AREGJ(JR,J,J_ERUN) +ERUN0*PLICE*DXYPJ ! (Tg=0)
        AREGJ(JR,J,J_IMPLM)=AREGJ(JR,J,J_IMPLM)+DIFS *PLICE*DXYPJ
        AREGJ(JR,J,J_IMPLH)=AREGJ(JR,J,J_IMPLH)+ERUN2*PLICE*DXYPJ
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS*PLICE
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+ERUN2*PLICE
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0 *PLICE
      END IF
      END DO
      END DO

      END SUBROUTINE PRECIP_LI

      SUBROUTINE GROUND_LI
!@sum  GROUND_LI driver for applying surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:LNDICE
      USE MODEL_COM, only : im,jm,flice,itlandi,itocean,itoice
      USE GEOM, only : imaxj,dxyp,bydxyp
      USE LANDICE, only : lndice,ace1li,ace2li
      USE SEAICE_COM, only : rsi
      USE SEAICE, only : rhos
      USE DIAG_COM, only : aj=>aj_loc,aregj=>aregj_loc,aij=>aij_loc
     *     ,jreg,ij_runli,ij_f1li,ij_erun2
     *     ,j_wtr1,j_ace1,j_wtr2,j_ace2,j_snow,j_run
     *     ,j_implh,j_implm,j_rsnow,ij_rsnw,ij_rsit,ij_snow,ij_f0oc
     *     ,j_rvrd,j_ervr,ij_mrvr,ij_ervr,ij_zsnow,ij_fwoc,ij_li
      USE LANDICE_COM, only : snowli,tlandi,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,ntm,trsnowli,trlndi,trli0,trdwnimp
#endif
      USE FLUXES, only : e0,e1,evapor,gtemp,runoli,gmelt,egmelt
#ifdef TRACERS_WATER
     *     ,trunoli,trevapor,gtracer
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
#endif
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif   /* TNL: inserted */
      USE TRDIAG_COM, only : taijn=>taijn_loc , tij_rvr
#endif
#endif
      USE DOMAIN_DECOMP, only : GRID,GET
      USE DOMAIN_DECOMP, only : GLOBALSUM
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0,PLICE,DXYPJ
     *     ,SCOVLI
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
      INTEGER I,J,JR
      INTEGER :: J_0,J_1, J_0H, J_1H

      CALL GET(GRID,J_STRT=J_0      ,J_STOP=J_1
     &             ,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      DO J=J_0,J_1
      DXYPJ=DXYP(J)
      DO I=1,IMAXJ(J)
      PLICE=FLICE(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0.
      IF (PLICE.gt.0) THEN

        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)
        F0DT=E0(I,J,3)
        F1DT=E1(I,J,3)
        EVAP=EVAPOR(I,J,3)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J)
        TRSNOW(:)=TRSNOWLI(:,I,J)
        TREVAP(:)=TREVAPOR(:,3,I,J)
#ifdef TRACERS_DRYDEP
     *       -trdrydep(:,3,i,j)
#endif
#endif
        CALL LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)
C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J) = RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
C**** accumulate implicit fluxes for setting ocean balance
        MDWNIMP(I,J)=MDWNIMP(I,J)+DIFS *PLICE*DXYPJ
        EDWNIMP(I,J)=EDWNIMP(I,J)+EDIFS*PLICE*DXYPJ
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
        TRDWNIMP(:,I,J)=TRDWNIMP(:,I,J)+TRDIFS(:)*PLICE*DXYPJ
        IF (SNOW.gt.1d-5) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
        SCOVLI=0
        IF (SNOWLI(I,J).GT.0.) SCOVLI=PLICE
        AJ(J,J_RSNOW,ITLANDI)=AJ(J,J_RSNOW,ITLANDI)+SCOVLI
        AREGJ(JR,J,J_RSNOW)=AREGJ(JR,J,J_RSNOW)+SCOVLI*DXYPJ
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVLI
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*PLICE
        AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+PLICE
        AIJ(I,J,IJ_LI)  =AIJ(I,J,IJ_LI)  +PLICE
        AIJ(I,J,IJ_ZSNOW)=AIJ(I,J,IJ_ZSNOW)+PLICE*SNOW/RHOS

        AJ(J,J_RUN,ITLANDI)  =AJ(J,J_RUN,ITLANDI)  +RUN0  *PLICE
        AJ(J,J_SNOW,ITLANDI) =AJ(J,J_SNOW,ITLANDI) +SNOW  *PLICE
        AJ(J,J_ACE1,ITLANDI) =AJ(J,J_ACE1,ITLANDI) +ACE1LI*PLICE
        AJ(J,J_ACE2,ITLANDI) =AJ(J,J_ACE2,ITLANDI) +ACE2LI*PLICE
        AJ(J,J_IMPLH,ITLANDI)=AJ(J,J_IMPLH,ITLANDI)+EDIFS *PLICE
        AJ(J,J_IMPLM,ITLANDI)=AJ(J,J_IMPLM,ITLANDI)+DIFS  *PLICE

        AREGJ(JR,J,J_RUN)   =AREGJ(JR,J,J_RUN)   +RUN0  *PLICE*DXYPJ
        AREGJ(JR,J,J_SNOW)  =AREGJ(JR,J,J_SNOW)  +SNOW  *PLICE*DXYPJ
        AREGJ(JR,J,J_ACE1)  =AREGJ(JR,J,J_ACE1)  +ACE1LI*PLICE*DXYPJ
        AREGJ(JR,J,J_ACE2)  =AREGJ(JR,J,J_ACE2)  +ACE2LI*PLICE*DXYPJ
        AREGJ(JR,J,J_IMPLH) =AREGJ(JR,J,J_IMPLH) +EDIFS *PLICE*DXYPJ
        AREGJ(JR,J,J_IMPLM) =AREGJ(JR,J,J_IMPLM) +DIFS  *PLICE*DXYPJ

        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +(EDIFS+F1DT)*PLICE
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0 *PLICE
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+EDIFS*PLICE
      END IF

C**** Accumulate diagnostics related to iceberg flux here also
      AJ(J,J_RVRD,ITOCEAN) = AJ(J,J_RVRD,ITOCEAN)+(1.-RSI(I,J))
     *     * GMELT(I,J)*BYDXYP(J)
      AJ(J,J_ERVR,ITOCEAN) = AJ(J,J_ERVR,ITOCEAN)+(1.-RSI(I,J))
     *     *EGMELT(I,J)*BYDXYP(J)
      AJ(J,J_RVRD,ITOICE)  = AJ(J,J_RVRD,ITOICE) +    RSI(I,J)
     *     * GMELT(I,J)*BYDXYP(J)
      AJ(J,J_ERVR,ITOICE)  = AJ(J,J_ERVR,ITOICE) +    RSI(I,J)
     *     *EGMELT(I,J)*BYDXYP(J)
      AIJ(I,J,IJ_F0OC) = AIJ(I,J,IJ_F0OC)+EGMELT(I,J)*BYDXYP(J)
      AIJ(I,J,IJ_FWOC) = AIJ(I,J,IJ_FWOC)+ GMELT(I,J)*BYDXYP(J)

      AREGJ(JR,J,J_RVRD) = AREGJ(JR,J,J_RVRD) +  GMELT(I,J)
      AREGJ(JR,J,J_ERVR) = AREGJ(JR,J,J_ERVR) + EGMELT(I,J)
      AIJ(I,J,IJ_MRVR)=AIJ(I,J,IJ_MRVR) +  GMELT(I,J)
      AIJ(I,J,IJ_ERVR)=AIJ(I,J,IJ_ERVR) + EGMELT(I,J)
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
      TAIJN(I,J,TIJ_RVR,:)=TAIJN(I,J,TIJ_RVR,:)+ TRGMELT(:,I,J)
     *     *BYDXYP(J)
#endif
#endif  /* TNL: inserted */
C****
      END DO
      END DO

      END SUBROUTINE GROUND_LI

      SUBROUTINE conserv_MLI(ICE)
!@sum  conserv_MLI calculates total amount of snow and ice in land ice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : snowli
      USE LANDICE, only : lndice,ace1li,ace2li
      USE DOMAIN_DECOMP, only : GRID,GET
      IMPLICIT NONE
!@var ICE total land ice snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: ICE
      INTEGER I,J
      INTEGER :: J_0,J_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        ICE(J)=0
        DO I=1,IMAXJ(J)
          ICE(J)=ICE(J)+FLICE(I,J)*(ACE1LI+ACE2LI+SNOWLI(I,J))
        END DO
      END DO
      IF(HAVE_SOUTH_POLE) ICE(1) =FIM*ICE(1)
      IF(HAVE_NORTH_POLE) ICE(JM)=FIM*ICE(JM)
      RETURN
C****
      END SUBROUTINE conserv_MLI

      SUBROUTINE conserv_HLI(EICE)
!@sum  conserv_HLI calculates total land ice energy
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shi,lhm
      USE MODEL_COM, only : im,jm,fim,flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : ace1li,ace2li
      USE DOMAIN_DECOMP, only : GRID,GET
      IMPLICIT NONE
!@var EICE total land ice energy (J/m^2)
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: EICE
      INTEGER I,J
      INTEGER :: J_0,J_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)


      DO J=J_0,J_1
        EICE(J)=0
        DO I=1,IMAXJ(J)
          EICE(J)=EICE(J)+FLICE(I,J)*((TLANDI(1,I,J)*SHI-LHM)*(ACE1LI
     *         +SNOWLI(I,J))+(TLANDI(2,I,J)*SHI-LHM)*ACE2LI)
        END DO
      END DO
      IF(HAVE_SOUTH_POLE) EICE(1) =FIM*EICE(1)
      IF(HAVE_NORTH_POLE) EICE(JM)=FIM*EICE(JM)
      RETURN
C****
      END SUBROUTINE conserv_HLI

      SUBROUTINE daily_LI
!@sum  daily_ice does daily landice things
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : edpery,sday,lhm,shi
      USE MODEL_COM, only : im,jm,flice,focean,dtsrc,jday,jyear
     *     ,itime,itimei,nday,JDperY
      USE GEOM, only : dxyp,imaxj
      USE LANDICE, only: ace1li,ace2li,glmelt_on,glmelt_fac_nh
     *     ,glmelt_fac_sh,fwarea_sh,fwarea_nh,accpda,accpdg,eaccpda
     *     ,eaccpdg
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,traccpda,traccpdg
#endif
#endif    /* TNL: inserted */
      USE LANDICE_COM, only : tlandi,snowli,mdwnimp,edwnimp,loc_gla
     *     ,loc_glg
#ifdef TRACERS_WATER
     *     ,ntm,trsnowli,trlndi,trli0,trdwnimp
#endif
      USE FLUXES, only : gtemp,gmelt,egmelt
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif
#endif    /* TNL: inserted */
#ifdef TRACERS_WATER
     *     ,gtracer
      USE TRACER_COM, only :  trw0
#endif
      USE PARAM
      USE DOMAIN_DECOMP, only : GRID, GET, GLOBALSUM, AM_I_ROOT,
     *     ESMF_BCAST
      IMPLICIT NONE
!@var gm_relax Glacial Melt relaxation parameter (1/year)
      REAL*8, PARAMETER :: gm_relax = 0.1d0  ! 10 year relaxation

      REAL*8 mdwnimp_SH,mdwnimp_NH,edwnimp_SH,edwnimp_NH,FAC_SH,FAC_NH
#ifdef TRACERS_WATER
      REAL*8 trdwnimp_SH(NTM),trdwnimp_NH(NTM)
#endif
      REAL*8 :: HSUM(2),GSUM
      INTEGER :: J_0,J_1,I,J,ITM

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)

C**** Every year update gmelt factors in order to balance downward
C**** implicit fluxes. If this is used, then ice sheets/snow are FORCED to
C**** be in balance. This may not be appropriate for transient runs but
C**** we aren't getting that right anyway.

      IF (JDAY.eq.1) THEN   ! Jan 1. only

! only adjust after at least one full year
        IF (itime.ge.itimei+JDperY*nday .and. glmelt_on==1) THEN

! mass and energy (kg, J)
          CALL GLOBALSUM(grid, MDWNIMP, gsum, hsum ,ALL=.TRUE.)
          mdwnimp_SH=hsum(1) ;  mdwnimp_NH=hsum(2)
          CALL GLOBALSUM(grid, EDWNIMP, gsum, hsum ,ALL=.TRUE.)
          edwnimp_SH=hsum(1) ;  edwnimp_NH=hsum(2)

#ifdef TRACERS_WATER
          DO ITM=1,NTM
            CALL GLOBALSUM(grid, TRDWNIMP(ITM,:,:), gsum, hsum ,ALL=
     *           .TRUE.)
            trdwnimp_SH(ITM)=hsum(1) ;  trdwnimp_NH(ITM)=hsum(2)
            print*,'Tracers:',(trdwnimp_NH(ITM)/mdwnimp_NH/trw0(itm)-1.
     *           )*1000,(trdwnimp_SH(ITM)/mdwnimp_SH/trw0(itm)-1.)*1000
     *           ,trdwnimp_NH(ITM),mdwnimp_NH,trdwnimp_SH(ITM)
     *           ,mdwnimp_SH
          END DO
#endif

C**** adjust hemispheric mean glacial melt amounts (only on root processor)
      if (AM_I_ROOT()) THEN
          write(6,*) "Adjusting glacial melt: ", jday,jyear
          write(6,*) "Mass (before): ",accpda,accpdg,mdwnimp_SH
     *         ,mdwnimp_NH
          write(6,*) "Temp (before): ",(eaccpda/accpda+lhm)/shi
     *         ,(eaccpdg/accpdg+lhm)/shi
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
          write(6,*) "Tracers (before)",1000*(traccpda(:)/accpda/trw0(:)
     *         -1.)
#endif
#endif   /* TNL: inserted */
           accpda =  accpda + gm_relax*(mdwnimp_SH -  accpda)
           accpdg =  accpdg + gm_relax*(mdwnimp_NH -  accpdg)
          eaccpda = eaccpda + gm_relax*(edwnimp_SH - eaccpda)
          eaccpdg = eaccpdg + gm_relax*(edwnimp_NH - eaccpdg)
          write(6,*) "Mass (after): ",accpda,accpdg
          write(6,*) "Temp (after): ",(eaccpda/accpda+lhm)/shi,
     *         (eaccpdg/accpdg+lhm)/shi
        END IF
        call ESMF_BCAST(grid,  ACCPDA)
        call ESMF_BCAST(grid,  ACCPDG)
        call ESMF_BCAST(grid, EACCPDA)
        call ESMF_BCAST(grid, EACCPDG)

#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
        if (AM_I_ROOT()) THEN
          traccpda(:)=traccpda(:)+gm_relax*(trdwnimp_SH(:)-traccpda(:))
          traccpdg(:)=traccpdg(:)+gm_relax*(trdwnimp_NH(:)-traccpdg(:))

          write(6,*) "Tracers (after)",1000*(traccpda(:)/accpda/trw0(:)
     *         -1.)
        ENDIF
        call ESMF_BCAST(grid, TRACCPDA)
        call ESMF_BCAST(grid, TRACCPDG)
#endif
#endif   /* TNL: inserted */
! accumulation (kg per source time step) per water column
      FAC_SH=DTsrc/(EDPERY*SDAY*FWAREA_SH)
      FAC_NH=DTsrc/(EDPERY*SDAY*FWAREA_NH)

C**** Set GL MELT arrays
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF (LOC_GLA(I,J)) THEN
              GMELT(I,J)  =  ACCPDA*FAC_SH*DXYP(J)*FOCEAN(I,J)  ! kg
              EGMELT(I,J) = EACCPDA*FAC_SH*DXYP(J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
              TRGMELT(:,I,J)=TRACCPDA(:)*FAC_SH*DXYP(J)*FOCEAN(I,J)  ! kg
#endif
#endif    /* TNL: inserted */
          END IF
          IF (LOC_GLG(I,J)) THEN
             GMELT(I,J) =  ACCPDG*FAC_NH*DXYP(J)*FOCEAN(I,J) ! kg
            EGMELT(I,J) = EACCPDG*FAC_NH*DXYP(J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            TRGMELT(:,I,J) = TRACCPDG(:)*FAC_NH*DXYP(J)*FOCEAN(I,J)  ! kg
#endif
#endif    /* TNL: inserted */
          END IF
        END DO
      END DO

      END IF  ! run not in its first year

C**** reset implicit accumulators
      MDWNIMP=0.
      EDWNIMP=0.
#ifdef TRACERS_WATER
      TRDWNIMP=0.
#endif

      END IF  ! start of new year

      END SUBROUTINE daily_LI


