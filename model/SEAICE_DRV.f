!@sum  SEAICE_DRV contains drivers for SEAICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont PRECIP_SI,GROUND_SI

      SUBROUTINE PRECIP_SI
!@sum  PRECIP_SI driver for applying precipitation to sea ice fraction
!@auth Original Development team
!@ver  1.0
!@calls PRECSI
      USE CONSTANT, only : byshi,lhm
      USE E001M12_COM, only : im,jm,fland,kocean,itoice,itlkice,focean
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi,prec,eprec
      USE SEAICE_COM, only : rsi,msi,snowi,hsi
      USE SEAICE, only : prec_si, ace1i, lmi,xsi
      USE DAGCOM, only : aj,areg,aij,jreg,ij_f0oi,ij_erun2
     *     ,j_difs,j_run1,j_edifs,j_erun2,j_imelt
      USE FLUXES, only : gtemp
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,TSIL
      REAL*8 SNOW,MSI2,PRCP,ENRGP,RUN0,DIFS,EDIFS,ERUN2,DXYPJ,POICE,th
     *     ,thn,tm,tmn
      LOGICAL QFIXR
      INTEGER I,J,IMAX,JR,ITYPE,N

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
        JR=JREG(I,J)
      POICE=RSI(I,J)*(1.-FLAND(I,J))
      RUNOSI(I,J)=0
      IF (POICE.gt.0) THEN

        IF (FOCEAN(I,J).gt.0) THEN
          ITYPE=ITOICE
          QFIXR = .TRUE.
          IF (KOCEAN.eq.1) QFIXR=.FALSE.
        ELSE
          ITYPE=ITLKICE
          QFIXR=.FALSE.
        END IF
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)      ! energy of precip
        SNOW=SNOWI(I,J)
        MSI2=MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice temperatures

        AIJ(I,J,IJ_F0OI)=AIJ(I,J,IJ_F0OI)+ENRGP*POICE
C**** CALL SUBROUTINE FOR CALCULATION OF PRECIPITATION OVER SEA ICE

        CALL PREC_SI(SNOW,MSI2,HSIL,TSIL,PRCP,ENRGP,RUN0,DIFS,EDIFS
     *       ,ERUN2,QFIXR)

        SNOWI(I,J)=SNOW
        RUNOSI(I,J) =RUN0
        HSI(1:2,I,J)=HSIL(1:2)
C**** set gtemp array
        GTEMP(1:2,2,I,J)=TSIL(1:2)
        IF (.not. QFIXR) THEN
          HSI(3:4,I,J)=HSIL(3:4)
          MSI(I,J)=MSI2
        END IF

C**** ACCUMULATE DIAGNOSTICS
        IF (QFIXR) THEN
          AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)+DIFS *POICE
          AJ(J,J_ERUN2,ITYPE)=AJ(J,J_ERUN2,ITYPE)+ERUN2*POICE
          AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFS*POICE*DXYPJ
        END IF
        AJ(J,J_DIFS, ITYPE)=AJ(J,J_DIFS, ITYPE)+DIFS *POICE
        AJ(J,J_RUN1, ITYPE)=AJ(J,J_RUN1, ITYPE)+RUN0 *POICE
        AJ(J,J_EDIFS,ITYPE)=AJ(J,J_EDIFS,ITYPE)+EDIFS*POICE
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
      USE CONSTANT, only : lhm,byshi,rhow
      USE E001M12_COM, only : im,jm,dtsrc,fland,kocean,focean
     *     ,itoice,itlkice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : e0,e1,evapor,runosi,erunosi,gtemp
      USE SEAICE_COM, only : rsi,msi,snowi,hsi
      USE SEAICE, only : sea_ice,lmi,xsi,ace1i
      USE LAKES_COM, only : mwl,gml,flake
      USE DAGCOM, only : aj,areg,aij,jreg,
     *     ij_f0oi,ij_erun2,ij_rsoi,ij_msi2,ij_evapi,j_difs
     *     ,j_run1,j_edifs,j_erun2,j_imelt,j_f1dt,j_f2dt,j_evap,ij_evap
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL
      REAL*8 SNOW,ROICE,MSI2,F0DT,F1DT,EVAP,TGW,RUN0
     *     ,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,DXYPJ,POICE,PWATER
      LOGICAL QFIXR, QFLUXLIM
      REAL*8 FLUXLIM
      INTEGER I,J,IMAX,JR,ITYPE

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
        F1DT=E1(I,J,2) ! heat flux between 1st and 2nd ice layer (J/m^2)
        EVAP=EVAPOR(I,J,2) ! evaporation/dew at the ice surface (kg/m^2)
        SNOW= SNOWI(I,J)  ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)  ! first layer sea ice enthalpy
        TGW=GTEMP(1,1,I,J)    ! underneath water temperature
        QFLUXLIM=.FALSE.
        IF (FOCEAN(I,J).gt.0) THEN
          ITYPE=ITOICE
          QFIXR = .TRUE.
          IF (KOCEAN.eq.1) QFIXR=.FALSE.
        ELSE
          ITYPE=ITLKICE
          QFIXR=.FALSE.
C**** Limit lake-to-ice flux if lake is too shallow (< 20cm)
          IF (MWL(I,J).lt.0.2d0*RHOW*FLAKE(I,J)*DXYP(J)) THEN
            QFLUXLIM=.TRUE.
            FLUXLIM= -GML(I,J)/(FLAKE(I,J)*DXYP(J))
          END IF
        END IF
        TGW=MIN(TGW,0d0) ! don't ever allow TGW>0 (should be tfrez?)

        AIJ(I,J,IJ_RSOI) =AIJ(I,J,IJ_RSOI) +POICE
        AIJ(I,J,IJ_MSI2) =AIJ(I,J,IJ_MSI2) +MSI2*POICE
        AIJ(I,J,IJ_F0OI) =AIJ(I,J,IJ_F0OI) +F0DT*POICE
        AIJ(I,J,IJ_EVAPI)=AIJ(I,J,IJ_EVAPI)+EVAP*POICE

        
        CALL SEA_ICE(DTSRC,SNOW,ROICE,HSIL,MSI2,F0DT,F1DT,EVAP,TGW,
     *       RUN0,DIFSI,EDIFSI,DIFS,EDIFS,ACE2M,F2DT,QFIXR,
     *       QFLUXLIM,FLUXLIM)
        
C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) =SNOW
        HSI(:,I,J) =HSIL(:)
        IF (.not. QFIXR) MSI(I,J) = MSI2

        RUNOSI(I,J) = RUN0+ACE2M
        ERUNOSI(I,J)= F2DT
C**** ACCUMULATE DIAGNOSTICS
        IF (.not. QFIXR) THEN
          AJ(J,J_DIFS ,ITYPE)=AJ(J,J_DIFS ,ITYPE)+DIFSI *PWATER
          AJ(J,J_EDIFS,ITYPE)=AJ(J,J_EDIFS,ITYPE)+EDIFSI*PWATER
c         AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)+ACE2M *POICE  ???
          IF (JR.ne.24) AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFSI*PWATER
     *         *DXYPJ
        ELSE
          AJ(J,J_DIFS ,ITYPE)=AJ(J,J_DIFS ,ITYPE)+DIFS *POICE
          AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)+DIFS *POICE
          AJ(J,J_EDIFS,ITYPE)=AJ(J,J_EDIFS,ITYPE)+EDIFS*POICE
          AJ(J,J_ERUN2,ITYPE)=AJ(J,J_ERUN2,ITYPE)+EDIFS*POICE
          IF (JR.ne.24) THEN
            AREG(JR,J_DIFS) =AREG(JR,J_DIFS) +DIFS *POICE*DXYPJ
            AREG(JR,J_ERUN2)=AREG(JR,J_ERUN2)+EDIFS*POICE*DXYPJ
          END IF
        END IF
        AJ(J,J_RUN1,ITYPE)=AJ(J,J_RUN1,ITYPE)+RUN0*POICE
        AJ(J,J_F1DT,ITYPE)=AJ(J,J_F1DT,ITYPE)+F1DT*POICE
        AJ(J,J_F2DT,ITYPE)=AJ(J,J_F2DT,ITYPE)+F2DT*POICE
        AJ(J,J_EVAP,ITYPE)=AJ(J,J_EVAP,ITYPE)+EVAP*POICE
        IF (JR.ne.24) THEN
          AREG(JR,J_RUN1)=AREG(JR,J_RUN1)+RUN0*POICE*DXYPJ
          AREG(JR,J_F2DT)=AREG(JR,J_F2DT)+F2DT*POICE*DXYPJ
        END IF
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
      USE E001M12_COM, only : im,jm,focean,kocean,ftype,fland
     *     ,itocean,itoice,itlake,itlkice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi,erunosi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi
      USE SEAICE, only : ace1i,addice,lmi,fleadoc,fleadlk
      USE DAGCOM, only : aj,areg,aij,jreg,j_difs,j_tg1,j_tg2,j_rsi
     *     ,j_ace1,j_ace2,j_snow,j_edifs,ij_tg1,j_type
      USE FLUXES, only : dmsi,dhsi,gtemp
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,TSIL
      REAL*8 SNOW,ROICE,MSI2,DIFSI,EDIFSI,ENRGFO,ACEFO,ACE2F,ENRGFI
     *     ,DXYPJ,POICE,PWATER,FLEAD
      LOGICAL QFIXR,QCMPR
      INTEGER I,J,IMAX,JR,ITYPE,ITYPEO

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      PWATER=1.-FLAND(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      JR=JREG(I,J)
      IF (PWATER.gt.0) THEN

        SNOW= SNOWI(I,J)      ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice enthalpy

        IF (FOCEAN(I,J).gt.0) THEN
          FLEAD=FLEADOC
          ITYPE=ITOICE
          ITYPEO=ITOCEAN
          IF (KOCEAN.eq.1) THEN
            QFIXR=.FALSE.
            QCMPR=.TRUE.
          ELSE
            QFIXR=.TRUE.
            QCMPR=.FALSE.
          END IF
        ELSE
          FLEAD=FLEADLK
          ITYPE=ITLKICE
          ITYPEO=ITLAKE
          QFIXR=.FALSE.
c          QCMPR=.FALSE.
          QCMPR=.TRUE.   ! allow compression now that minimum is thinner
        END IF

        ACEFO=DMSI(1,I,J)
        ACE2F=DMSI(2,I,J)
        ENRGFO=DHSI(1,I,J)
        ENRGFI=DHSI(2,I,J)

        CALL ADDICE (SNOW,ROICE,HSIL,MSI2,TSIL,DIFSI,EDIFSI
     *       ,ENRGFO,ACEFO,ACE2F,ENRGFI,FLEAD,QFIXR,QCMPR)

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) =SNOW
        HSI(:,I,J) =HSIL(:)
        IF (.not. QFIXR) THEN
          RSI(I,J)=ROICE
          MSI(I,J)=MSI2
C**** set ftype arrays
          FTYPE(ITYPE ,I,J)=(1.-FLAND(I,J))*RSI(I,J)
          FTYPE(ITYPEO,I,J)=(1.-FLAND(I,J))-FTYPE(ITYPE,I,J)
        END IF
C**** set gtemp array
        GTEMP(1:2,2,I,J)=TSIL(1:2)

C**** ACCUMULATE DIAGNOSTICS
        IF (QCMPR) THEN
          AJ(J,J_DIFS, ITYPE)=AJ(J,J_DIFS, ITYPE)+DIFSI *PWATER
          AJ(J,J_EDIFS,ITYPE)=AJ(J,J_EDIFS,ITYPE)+EDIFSI*PWATER
          IF (JR.ne.24) AREG(JR,J_DIFS)=AREG(JR,J_DIFS)+DIFSI*PWATER
     *         *DXYPJ
        END IF
        AJ(J,J_TG1, ITYPE)=AJ(J,J_TG1, ITYPE)+TSIL(1)*POICE
        AJ(J,J_TG2, ITYPE)=AJ(J,J_TG2, ITYPE)+TSIL(2)*POICE
        AJ(J,J_RSI, ITYPE)=AJ(J,J_RSI, ITYPE)+        POICE
        AJ(J,J_ACE2,ITYPE)=AJ(J,J_ACE2,ITYPE)+MSI2   *POICE
        AJ(J,J_SNOW,ITYPE)=AJ(J,J_SNOW,ITYPE)+SNOW   *POICE
        IF (JR.ne.24) THEN
        AREG(JR,J_TG1) =AREG(JR,J_TG1) +TSIL(1)*POICE*DXYPJ
        AREG(JR,J_TG2) =AREG(JR,J_TG2) +TSIL(2)*POICE*DXYPJ
        AREG(JR,J_RSI) =AREG(JR,J_RSI) +        POICE*DXYPJ
        AREG(JR,J_SNOW)=AREG(JR,J_SNOW)+SNOW   *POICE*DXYPJ
        AREG(JR,J_ACE1)=AREG(JR,J_ACE1)+ACE1I  *POICE*DXYPJ
        AREG(JR,J_ACE2)=AREG(JR,J_ACE2)+MSI2   *POICE*DXYPJ
        END IF
        AIJ(I,J,IJ_TG1)=AIJ(I,J,IJ_TG1)+TSIL(1)*POICE

      END IF
      END DO
      END DO
C****
      END SUBROUTINE FORM_SI

      SUBROUTINE vflx_OCEAN
!@sum  vflx_OCEAN saves quantities for OHT calculations
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : im,jm,focean
      USE CONSTANT, only : byshi,lhm
      USE DAGCOM, only : oa
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : msi,hsi,snowi
      IMPLICIT NONE
      INTEGER I,J
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  TG1OI  (INSTANTANEOUS AT NOON GMT)
C****       3  TG2OI  (INSTANTANEOUS AT NOON GMT)
C****
      DO J=1,JM
        DO I=1,IM
          IF (FOCEAN(I,J).gt.0) THEN
            OA(I,J,1)=ACE1I+SNOWI(I,J)
            OA(I,J,2)=((HSI(1,I,J)+HSI(2,I,J))/(ACE1I+SNOWI(I,J)) + LHM)
     *           *BYSHI
            OA(I,J,3)=((HSI(3,I,J)+HSI(4,I,J))/MSI(I,J) + LHM)*BYSHI
c            OA(I,J,2)=TSI(1,I,J)*XSI(1)+TSI(2,I,J)*XSI(2)
c            OA(I,J,3)=TSI(3,I,J)*XSI(3)+TSI(4,I,J)*XSI(4)
          END IF
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE vflx_OCEAN

      SUBROUTINE init_ice
!@sum  init_ice initialises ice arrays
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : byshi,lhm
      USE E001M12_COM, only : im,jm,kocean
      USE SEAICE_COM, only : rsi,msi,hsi,snowi
      USE SEAICE, only : xsi,ace1i,ac2oim
      USE FLUXES, only : gtemp
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 MSI1

      IF (KOCEAN.EQ.0) THEN
C****   set defaults for no ice case
        DO J=1,JM
        DO I=1,IM
           IF (RSI(I,J).eq.0) THEN
             MSI1        =ACE1I
             MSI(I,J)    =AC2OIM
             SNOWI(I,J)  =0.
             HSI(1:2,I,J)=-LHM*XSI(1:2)*MSI1
             HSI(3:4,I,J)=-LHM*XSI(3:4)*AC2OIM
           END IF
        END DO
        END DO
      END IF
C**** set GTEMP array for ice
      DO J=1,JM
      DO I=1,IM
        MSI1=SNOWI(I,J)+ACE1I
        GTEMP(1:2,2,I,J)=(HSI(1:2,I,J)/(XSI(1:2)*MSI1)+LHM)*BYSHI
      END DO
      END DO

C**** Set conservation diagnostics for ice mass and energy
c      QCON=(/ F, F, F, T, T, T, F, F, T/)
c      CALL SET_CON(QCON,"ICE MASS","(10**10 KG/M^2)  ",
c     *     "(10**3 KG/S/M^2)",1d-10,1d-3)
c      QCON=(/ F, F, F, T, T, T, F, F, T/)
c      CALL SET_CON(QCON,"ICE ENRG","(10**14 J/M^2)   ",
c     *     "(10**8 J/S/M^2) ",1d-14,1d-8)
C****
      END SUBROUTINE init_ice

      SUBROUTINE conserv_ICE(ICE)
!@sum  conserv_ICE calculates total amount of snow and ice over water
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : im,jm,fim,fland
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : ace1i
      IMPLICIT NONE
!@var ICE total snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(JM) :: ICE
      INTEGER I,J

      DO J=1,JM
        ICE(J)=0
        DO I=1,IMAXJ(J)
          ICE(J)=ICE(J)+RSI(I,J)*(MSI(I,J)+ACE1I+SNOWI(I,J))
     *         *(1.-FLAND(I,J))
        END DO
      END DO
      ICE(1) =FIM*ICE(1)
      ICE(JM)=FIM*ICE(JM)
      RETURN
C****
      END SUBROUTINE conserv_ICE

      SUBROUTINE conserv_EIC(EICE)
!@sum  conserv_EIC calculates total ice energy over water
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : im,jm,fim,fland
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi,hsi
      IMPLICIT NONE
!@var EICE total snow and ice energy (J/m^2)
      REAL*8, DIMENSION(JM) :: EICE
      INTEGER I,J

      DO J=1,JM
        EICE(J)=0
        DO I=1,IMAXJ(J)
          EICE(J)=EICE(J)+RSI(I,J)*(1.-FLAND(I,J))*
     *         (HSI(1,I,J)+HSI(2,I,J)+HSI(3,I,J)+HSI(4,I,J))
        END DO
      END DO
      EICE(1) =FIM*EICE(1)
      EICE(JM)=FIM*EICE(JM)
      RETURN
C****
      END SUBROUTINE conserv_EIC



