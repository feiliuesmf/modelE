!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont init_LI,PRECIP_LI,GROUND_LI

      SUBROUTINE init_LI
!@sum  init_ice initialises landice arrays
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE LANDICE_COM, only : tlandi
      USE FLUXES, only : gtemp
      IMPLICIT NONE
      INTEGER I,J

C**** set GTEMP array for landice
      DO J=1,JM
        DO I=1,IM
          GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
        END DO
      END DO
C****
      END SUBROUTINE init_LI

      SUBROUTINE PRECIP_LI
!@sum  PRECIP_LI driver for applying precipitation to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:PRECLI
      USE MODEL_COM, only : im,jm,flice,itlandi
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runoli,prec,eprec
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : precli
      USE DAGCOM, only : aj,areg,aij,jreg,
     *     ij_f0li,ij_f1li,ij_erun2,ij_runli,j_difs
     *     ,j_run1,j_edifs,j_erun2,j_imelt,j_type
      USE FLUXES, only : gtemp
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0,PLICE,DXYPJ
      INTEGER I,J,IMAX,JR

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      PLICE=FLICE(I,J)
      PRCP=PREC(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0
      IF (PLICE.gt.0) THEN

        ENRGP=EPREC(I,J)      ! energy of precipitation
        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)

        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP

        CALL PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J)  =RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
C**** ACCUMULATE DIAGNOSTICS
        AJ(J,J_DIFS, ITLANDI)=AJ(J,J_DIFS, ITLANDI)+DIFS *PLICE
        AJ(J,J_RUN1, ITLANDI)=AJ(J,J_RUN1, ITLANDI)+RUN0 *PLICE
        AJ(J,J_TYPE, ITLANDI)=AJ(J,J_TYPE, ITLANDI)+      PLICE
        AJ(J,J_IMELT,ITLANDI)=AJ(J,J_IMELT,ITLANDI)+DIFS *PLICE
        AJ(J,J_EDIFS,ITLANDI)=AJ(J,J_EDIFS,ITLANDI)+EDIFS*PLICE
        AJ(J,J_ERUN2,ITLANDI)=AJ(J,J_ERUN2,ITLANDI)+ERUN2*PLICE
C       AJ(J,J_ERUN1,ITLANDI)=AJ(J,J_ERUN1,ITLANDI)+ERUN0*PLICE ! (Tg=0)
        AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS*PLICE*DXYPJ
        AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0*PLICE*DXYPJ
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+ERUN2
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
      END IF
      END DO
      END DO
      END SUBROUTINE PRECIP_LI

      SUBROUTINE GROUND_LI
!@sum  GROUND_LI driver for applying surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:LNDICE
      USE MODEL_COM, only : im,jm,flice,itlandi
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runoli
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : lndice,ace1li,ace2li
      USE DAGCOM, only : aj,areg,aij,jreg,ij_evap,ij_f0li,ij_evapli
     *     ,ij_runli,ij_f1li,ij_erun2,ij_tg1,j_tg2,j_tg1,j_difs,j_wtr1
     *     ,j_ace1,j_wtr2,j_ace2,j_snow,j_run1,j_f2dt,j_edifs,j_f1dt
     *     ,j_erun2,j_imelt,j_run2,j_evap,j_rsnow,ij_rsnw
     *     ,ij_rsit,ij_snow
      USE FLUXES, only : e0,e1,evapor,gtemp
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0,PLICE,DXYPJ
     *     ,SCOVLI
      INTEGER I,J,IMAX,JR

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
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
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+F0DT
        AIJ(I,J,IJ_EVAPLI)=AIJ(I,J,IJ_EVAPLI)+EVAP

        CALL LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J) = RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
C**** ACCUMULATE DIAGNOSTICS
        SCOVLI=0
        IF (SNOWLI(I,J).GT.0.) SCOVLI=PLICE
        AJ(J,J_RSNOW,ITLANDI)=AJ(J,J_RSNOW,ITLANDI)+SCOVLI
        AREG(JR,J_RSNOW)=AREG(JR,J_RSNOW)+SCOVLI*DXYPJ
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVLI
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*PLICE
        AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+PLICE

        AJ(J,J_TG1,ITLANDI)  =AJ(J,J_TG1,ITLANDI)  +TG1  *PLICE
        AJ(J,J_TG2,ITLANDI)  =AJ(J,J_TG2,ITLANDI)  +TG2  *PLICE
        AJ(J,J_SNOW,ITLANDI) =AJ(J,J_SNOW,ITLANDI) +SNOW *PLICE
        AJ(J,J_F1DT,ITLANDI) =AJ(J,J_F1DT,ITLANDI) +F1DT *PLICE
        AJ(J,J_EVAP,ITLANDI) =AJ(J,J_EVAP,ITLANDI) +EVAP *PLICE
        AJ(J,J_RUN1,ITLANDI) =AJ(J,J_RUN1,ITLANDI) +RUN0 *PLICE
        AJ(J,J_DIFS,ITLANDI) =AJ(J,J_DIFS,ITLANDI) +DIFS *PLICE
        AJ(J,J_IMELT,ITLANDI)=AJ(J,J_IMELT,ITLANDI)+DIFS *PLICE
        AJ(J,J_EDIFS,ITLANDI)=AJ(J,J_EDIFS,ITLANDI)+EDIFS*PLICE
        AJ(J,J_ERUN2,ITLANDI)=AJ(J,J_ERUN2,ITLANDI)+EDIFS*PLICE
        IF (JR.ne.24) THEN
        AREG(JR,J_TG1) =AREG(JR,J_TG1) +TG1   *PLICE*DXYPJ
        AREG(JR,J_TG2) =AREG(JR,J_TG2) +TG2   *PLICE*DXYPJ
        AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0  *PLICE*DXYPJ
        AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS  *PLICE*DXYPJ
        AREG(JR,J_SNOW)=AREG(JR,J_SNOW)+SNOW  *PLICE*DXYPJ
        AREG(JR,J_ACE1)=AREG(JR,J_ACE1)+ACE1LI*PLICE*DXYPJ
        AREG(JR,J_ACE2)=AREG(JR,J_ACE2)+ACE2LI*PLICE*DXYPJ
        END IF
        AIJ(I,J,IJ_TG1)  =AIJ(I,J,IJ_TG1)  +TG1 *PLICE
        AIJ(I,J,IJ_EVAP) =AIJ(I,J,IJ_EVAP) +EVAP*PLICE
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS+F1DT
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+EDIFS
C****
      END IF
      END DO
      END DO
      END SUBROUTINE GROUND_LI
