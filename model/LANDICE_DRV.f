!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont PRECIP_LI,GROUND_LI

      SUBROUTINE PRECIP_LI
!@sum  PRECIP_LI driver for applying precipitation to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls PRECLI
      USE E001M12_COM, only : im,jm,flice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runoli,prec,eprec
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : precli
      USE DAGCOM, only : bj,areg,aij,jreg,
     *     ij_f0li,ij_f1li,ij_erun2,ij_runli,j_eprcp,j_difs
     *     ,j_run1,j_edifs,j_erun2,j_imelt
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
        BJ(J,J_EPRCP)=BJ(J,J_EPRCP)+ENRGP*PLICE
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP
        
        CALL PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0)
        
C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J)  =RUN0
C**** ACCUMULATE DIAGNOSTICS
        BJ(J,J_DIFS) =BJ(J,J_DIFS) +DIFS *PLICE
        BJ(J,J_RUN1) =BJ(J,J_RUN1) +RUN0 *PLICE
        BJ(J,J_IMELT)=BJ(J,J_IMELT)+DIFS *PLICE
        BJ(J,J_EDIFS)=BJ(J,J_EDIFS)+EDIFS*PLICE
        BJ(J,J_ERUN2)=BJ(J,J_ERUN2)+ERUN2*PLICE
C       BJ(J,J_ERUN1)=BJ(J,J_ERUN1)+ERUN0*PLICE ! land ice (Tg=0)
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
!@calls LNDICE
      USE E001M12_COM, only : im,jm,flice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runoli
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : lndice,ace1li,ace2li
      USE DAGCOM, only : bj,areg,aij,jreg,ij_evap,ij_f0li,ij_evapli
     *     ,ij_runli,ij_f1li,ij_erun2,ij_tg1,j_tg2,j_tg1,j_difs,j_wtr1
     *     ,j_ace1,j_wtr2,j_ace2,j_snow,j_run1,j_f2dt,j_edifs,j_f1dt
     *     ,j_erun2,j_imelt,j_run2,j_evap
      USE FLUXES, only : e0,e1,evapor
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0,PLICE,DXYPJ
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
C**** ACCUMULATE DIAGNOSTICS
        BJ(J,J_TG1)  =BJ(J,J_TG1)  +TG1  *PLICE
        BJ(J,J_TG2)  =BJ(J,J_TG2)  +TG2  *PLICE
        BJ(J,J_SNOW) =BJ(J,J_SNOW) +SNOW *PLICE
        BJ(J,J_F1DT) =BJ(J,J_F1DT) +F1DT *PLICE
        BJ(J,J_EVAP) =BJ(J,J_EVAP) +EVAP *PLICE
        BJ(J,J_RUN1) =BJ(J,J_RUN1) +RUN0 *PLICE
        BJ(J,J_DIFS) =BJ(J,J_DIFS) +DIFS *PLICE
        BJ(J,J_IMELT)=BJ(J,J_IMELT)+DIFS *PLICE
        BJ(J,J_EDIFS)=BJ(J,J_EDIFS)+EDIFS*PLICE
        BJ(J,J_ERUN2)=BJ(J,J_ERUN2)+EDIFS*PLICE
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
