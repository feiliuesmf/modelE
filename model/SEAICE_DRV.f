!@sum  SEAICE_DRV contains drivers for SEAICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont PRECIP_SI,GROUND_SI

      SUBROUTINE PRECIP_SI
!@sum  PRECIP_SI driver for applying precipitation to sea ice fraction
!@auth Original Development team
!@ver  1.0
!@calls PRECSI
      USE E001M12_COM, only : im,jm,fland,kocean
      USE GEOM, only : imaxj,dxyp
      USE CLD01_COM_E001, only : prec,tprec,eprec
      USE FLUXES, only : runosi
      USE SEAICE_COM, only : rsi,msi,snowi,tsi   !,hsi soon
      USE SEAICE, only : prec_si, ace1i
      USE DAGCOM, only : cj,areg,aij,jreg,
     *     ij_f0oi,ij_erun2,j_eprcp,j_difs
     *     ,j_run1,j_edifs,j_erun2,j_imelt
      IMPLICIT NONE

      REAL*8 SNOW,MSI2,TG1,TG2,TG3,TG4,PRCP,TPRCP,ENRGP
     *     ,RUN0,DIFS,EDIFS,ERUN2,DXYPJ,POICE
      LOGICAL QFIXR
      INTEGER I,J,IMAX,JR

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      POICE=RSI(I,J)*(1.-FLAND(I,J))
      PRCP=PREC(I,J)
        JR=JREG(I,J)
      RUNOSI(I,J)=0  
      IF (POICE.gt.0) THEN

        TPRCP=TPREC(I,J)
        ENRGP=EPREC(2,I,J)      ! including latent heat
        SNOW=SNOWI(I,J)
        MSI2=MSI(I,J)
        CJ(J,J_EPRCP)=CJ(J,J_EPRCP)+ENRGP*POICE
        AIJ(I,J,IJ_F0OI)=AIJ(I,J,IJ_F0OI)+ENRGP*POICE

        TG1 = TSI(1,I,J)      ! first layer sea ice temperature
        TG2 = TSI(2,I,J)      ! second layer sea ice temperature
        TG3 = TSI(3,I,J)      ! third layer sea ice temperature
        TG4 = TSI(4,I,J)      ! fourth layer sea ice temperature
        QFIXR = .TRUE.
        IF (KOCEAN.eq.1) QFIXR=.FALSE.
c       IF (FLAKE(I,J).gt.0) QFIXR=.FALSE.   ! will soon be implemented
        
C**** CALL SUBROUTINE FOR CALCULATION OF PRECIPITATION OVER SEA ICE
        
        CALL PREC_SI(SNOW,MSI2,TG1,TG2,TG3,TG4,PRCP,ENRGP
     *       ,RUN0,DIFS,EDIFS,ERUN2,QFIXR)
        
        SNOWI(I,J)=SNOW
        RUNOSI(I,J) =RUN0
        TSI(1,I,J)=TG1
        TSI(2,I,J)=TG2
        IF (.not. QFIXR) THEN
          TSI(3,I,J)=TG3
          TSI(4,I,J)=TG4
          MSI(I,J)=MSI2
        END IF
        
C**** ACCUMULATE DIAGNOSTICS
        IF (QFIXR) THEN
          CJ(J,J_IMELT)=CJ(J,J_IMELT)+DIFS *POICE
          CJ(J,J_ERUN2)=CJ(J,J_ERUN2)+ERUN2*POICE
          AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS*POICE*DXYPJ
        END IF
        CJ(J,J_DIFS) =CJ(J,J_DIFS) +DIFS *POICE
        CJ(J,J_RUN1) =CJ(J,J_RUN1) +RUN0 *POICE
        CJ(J,J_EDIFS)=CJ(J,J_EDIFS)+EDIFS*POICE
        AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0*POICE*DXYPJ
        
      END IF
      END DO
      END DO
C****
      END SUBROUTINE PRECIP_SI

      SUBROUTINE GROUND_SI
!@sum  GROUND_SI driver for applying surface fluxes to sea ice fraction
!@auth Original Development team
!@ver  1.0
!@calls SEA_ICE
      USE CONSTANT, only : lhm,byshi
      USE E001M12_COM, only : im,jm,dtsrc,fland,kocean,focean
     *     ,flake
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : e0,e1,evapor,runosi,erunosi
      USE SEAICE_COM, only : rsi,msi,snowi,tsi
      USE SEAICE, only : sea_ice,ace1i,xsi1,xsi2,xsi3,xsi4
      USE OCEAN, only : tocean
      USE LAKES_COM, only : tlake
      USE DAGCOM, only : aj,cj,areg,aij,jreg,
     *     ij_f0oi,ij_erun2,ij_rsoi,ij_msi2,ij_evapi,j_eprcp,j_difs
     *     ,j_run1,j_edifs,j_erun2,j_imelt,j_f1dt,j_f2dt,j_evap,ij_evap
      IMPLICIT NONE

      REAL*8 SNOW,ROICE,TG1,TG2,TG3,TG4,MSI2
     *     ,F0DT,F1DT,EVAP,HSI1,HSI2,HSI3,HSI4,TGW,RUN0
     *     ,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,DXYPJ,POICE,PWATER
      LOGICAL QFIXR
      INTEGER I,J,IMAX,JR

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      PWATER=FOCEAN(I,J)+FLAKE(I,J)   ! 1.-FLAND(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      JR=JREG(I,J)
      RUNOSI(I,J)=0  
      ERUNOSI(I,J)=0  
      IF (PWATER.gt.0) THEN
        
        F0DT=E0(I,J,2) ! heat flux to the top ice surface (J/m^2)
        F1DT=E1(I,J,2) ! heat flux between 1st and 2nd ice layers (J/m^2)
        EVAP=EVAPOR(I,J,2)  ! evaporation/dew at the ice surface (kg/m^2)
        SNOW= SNOWI(I,J)  ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        TG1 = TSI(1,I,J)  ! first layer sea ice temperature
        TG2 = TSI(2,I,J)  ! second layer sea ice temperature
        TG3 = TSI(3,I,J)  ! third layer sea ice temperature
        TG4 = TSI(4,I,J)  ! fourth layer sea ice temperature
        TGW = TOCEAN(1,I,J) ! ocean temperature
c        TGW = TLAKE(I,J) ! lake temperature

        AIJ(I,J,IJ_RSOI) =AIJ(I,J,IJ_RSOI) +POICE
        AIJ(I,J,IJ_MSI2) =AIJ(I,J,IJ_MSI2) +MSI2*POICE
        AIJ(I,J,IJ_F0OI) =AIJ(I,J,IJ_F0OI) +F0DT*POICE
        AIJ(I,J,IJ_EVAPI)=AIJ(I,J,IJ_EVAPI)+EVAP*POICE
        
        QFIXR = .TRUE.
        IF (KOCEAN.eq.1) QFIXR=.FALSE.
c       IF (FLAKE(I,J).gt.0) QFIXR=.FALSE.   ! will soon be implemented
        
        CALL SEA_ICE(DTSRC,SNOW,ROICE,TG1,TG2,TG3,TG4,MSI2
     *       ,F0DT,F1DT,EVAP,HSI1,HSI2,HSI3,HSI4,TGW,RUN0
     *       ,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,QFIXR)
        
C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) =SNOW
        TSI(1,I,J) =HSI1  !TG1   ! this is a temporary reassignment for 
        TSI(2,I,J) =HSI2  !TG2   ! bytewise compatible-ness in GROUND routines
        TSI(3,I,J) =HSI3  !TG3
        TSI(4,I,J) =HSI4  !TG4
        IF (.not. QFIXR) THEN
          MSI(I,J) = MSI2
        END IF

        RUNOSI(I,J) = RUN0+ACE2M
        ERUNOSI(I,J)= F2DT
C**** ACCUMULATE DIAGNOSTICS
        IF (KOCEAN .EQ. 1) THEN
          CJ(J,J_DIFS) =CJ(J,J_DIFS) +DIFSI *PWATER
          CJ(J,J_EDIFS)=CJ(J,J_EDIFS)+EDIFSI*PWATER
c     AJ(J,J_IMELT)=AJ(J,J_IMELT)+ACE2M *POICE  ???
          IF (JR.ne.24) AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFSI*PWATER
     *         *DXYPJ
        ELSE
          CJ(J,J_DIFS) =CJ(J,J_DIFS) +DIFS *POICE
          CJ(J,J_IMELT)=CJ(J,J_IMELT)+DIFS *POICE
          CJ(J,J_EDIFS)=CJ(J,J_EDIFS)+EDIFS*POICE
          CJ(J,J_ERUN2)=CJ(J,J_ERUN2)+EDIFS*POICE
          IF (JR.ne.24) AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS*POICE*DXYPJ
        END IF
        CJ(J,J_RUN1)=CJ(J,J_RUN1)+RUN0*POICE
        CJ(J,J_F1DT)=CJ(J,J_F1DT)+F1DT*POICE
        CJ(J,J_F2DT)=CJ(J,J_F2DT)+F2DT*POICE
        CJ(J,J_EVAP)=CJ(J,J_EVAP)+EVAP*POICE
        IF (JR.ne.24) AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0*POICE*DXYPJ
        AIJ(I,J,IJ_EVAP)=AIJ(I,J,IJ_EVAP)+EVAP*POICE
        
      END IF
      END DO
      END DO
C****       
      END SUBROUTINE GROUND_SI

      SUBROUTINE FORM_SI
!@sum  FORM_SI driver for adding new sea ice 
!@auth Original Development team
!@ver  1.0
!@calls SEA_ICE
      USE E001M12_COM, only : im,jm,focean,flake,kocean
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi,erunosi
      USE SEAICE_COM, only : rsi,msi,snowi,tsi
      USE SEAICE, only : ace1i,addice
      USE DAGCOM, only : cj,areg,aij,jreg,j_difs,j_tg1,j_tg2,j_rsi
     *     ,j_ace1,j_ace2,j_snow,j_edifs,ij_tg1,j_rsi
      USE OCEAN, only : tfo, fleadoc
c      USE LAKES, only : tfl, fleadlk
      USE FLUXES, only : dmsi,dhsi
      IMPLICIT NONE

      REAL*8 SNOW,ROICE,TG1,TG2,TG3,TG4,MSI2,HSI1,HSI2
     *     ,HSI3,HSI4,DIFSI,EDIFSI,ENRGFO,ACEFO,ACE2F,ENRGFI
     *     ,DXYPJ,POICE,PWATER
      LOGICAL QFIXR,QCMPR
      INTEGER I,J,IMAX,JR

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      PWATER=FOCEAN(I,J)+FLAKE(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      JR=JREG(I,J)
      IF (PWATER.gt.0) THEN
        
        SNOW= SNOWI(I,J)      ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
c        update HSI... (temporary assignment)
        HSI1 = TSI(1,I,J)      ! first layer sea ice temperature
        HSI2 = TSI(2,I,J)      ! second layer sea ice temperature
        HSI3 = TSI(3,I,J)      ! third layer sea ice temperature
        HSI4 = TSI(4,I,J)      ! fourth layer sea ice temperature

        ACEFO=DMSI(1,I,J)
        ACE2F=DMSI(2,I,J)
        ENRGFO=DHSI(1,I,J)
        ENRGFI=DHSI(2,I,J)

        QFIXR = .TRUE.
        QCMPR = .FALSE.
        IF (KOCEAN.eq.1) QFIXR=.FALSE.
        IF (KOCEAN.eq.1.and.FLAKE(I,J).le.0) QCMPR=.TRUE.

        CALL ADDICE (SNOW,ROICE,TG1,TG2,TG3,TG4,MSI2,HSI1,HSI2
     *       ,HSI3,HSI4,DIFSI,EDIFSI,ENRGFO,ACEFO,ACE2F,ENRGFI,TFO
     *       ,FLEADOC,QFIXR,QCMPR)

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) =SNOW
        TSI(1,I,J) =TG1
        TSI(2,I,J) =TG2
        TSI(3,I,J) =TG3
        TSI(4,I,J) =TG4
        IF (.not. QFIXR) THEN
          RSI(I,J)=ROICE
          MSI(I,J)=MSI2
        END IF

C**** ACCUMULATE DIAGNOSTICS
        IF (KOCEAN .EQ. 1) THEN
          CJ(J,J_DIFS) =CJ(J,J_DIFS) +DIFSI*PWATER
          CJ(J,J_EDIFS)=CJ(J,J_EDIFS)+EDIFSI*PWATER
          IF (JR.ne.24) AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFSI*PWATER
     *         *DXYPJ
        END IF
        CJ(J,J_TG1) =CJ(J,J_TG1) +TG1*POICE
        CJ(J,J_TG2) =CJ(J,J_TG2) +TG2*POICE
        CJ(J,J_RSI) =CJ(J,J_RSI) +POICE
        CJ(J,J_ACE2)=CJ(J,J_ACE2)+MSI2*POICE
        CJ(J,J_SNOW)=CJ(J,J_SNOW)+SNOW*POICE
        IF (JR.ne.24) THEN
        AREG(JR,J_RSI) =AREG(JR,J_RSI) +      POICE*DXYPJ
        AREG(JR,J_TG1) =AREG(JR,J_TG1) +TG1  *POICE*DXYPJ
        AREG(JR,J_TG2) =AREG(JR,J_TG2) +TG2  *POICE*DXYPJ
        AREG(JR,J_SNOW)=AREG(JR,J_SNOW)+SNOW *POICE*DXYPJ
        AREG(JR,J_ACE1)=AREG(JR,J_ACE1)+ACE1I*POICE*DXYPJ
        AREG(JR,J_ACE2)=AREG(JR,J_ACE2)+MSI2 *POICE*DXYPJ
        END IF
        AIJ(I,J,IJ_TG1)=AIJ(I,J,IJ_TG1)+TG1  *POICE

      END IF
      END DO
      END DO
C****       
      END SUBROUTINE FORM_SI
