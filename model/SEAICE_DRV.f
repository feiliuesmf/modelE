#include "rundeck_opts.h"

!@sum  SEAICE_DRV contains drivers for SEAICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont PRECIP_SI,GROUND_SI

      SUBROUTINE PRECIP_SI
!@sum  PRECIP_SI driver for applying precipitation to sea ice fraction
!@auth Original Development team
!@ver  1.0
!@calls seaice:prec_si
      USE CONSTANT, only : byshi,lhm,teeny,rhoi
      USE MODEL_COM, only : im,jm,fland,kocean,itoice,itlkice,focean
     *     ,jday,itocean,itlake
      USE GEOM, only : imaxj,dxyp,bydxyp
      USE FLUXES, only : runpsi,prec,eprec,srunpsi,gtemp
#ifdef TRACERS_WATER
     *     ,trprec,trunpsi,gtracer
#endif
      USE SEAICE, only : prec_si, ace1i, lmi,xsi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,flag_dsws,pond_melt
#ifdef TRACERS_WATER
     *     ,ntm,trsi
#endif
      USE DAGCOM, only : aj,areg,aij,jreg,ij_f0oi,ij_erun2,j_imelt
     *     ,j_smelt
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 SNOW,MSI2,PRCP,ENRGP,RUN0,POICE,SRUN0
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, DIMENSION(NTM) :: TRUN0,TRPRCP
#endif
      INTEGER I,J,JR,ITYPE,N,ITYPEO
      LOGICAL WETSNOW

      DO J=1,JM
      DO I=1,IMAXJ(J)
        JR=JREG(I,J)
      POICE=    RSI(I,J) *(1.-FLAND(I,J))
      RUNPSI(I,J)=0
      SRUNPSI(I,J)=0
      IF (POICE.gt.0) THEN

        IF (FOCEAN(I,J).gt.0) THEN
          ITYPE=ITOICE
          ITYPEO=ITOCEAN
        ELSE
          ITYPE=ITLKICE
          ITYPEO=ITLAKE
        END IF
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)      ! energy of precip
        SNOW=SNOWI(I,J)
        MSI2=MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice temperatures
        SSIL(:) = SSI(:,I,J)      ! sea ice salt
#ifdef TRACERS_WATER
        TRSIL(:,:)=TRSI(:,:,I,J)  ! sea ice tracers
        TRPRCP(:)=TRPREC(:,I,J)*BYDXYP(J)   ! tracer in precip
#endif

        AIJ(I,J,IJ_F0OI)=AIJ(I,J,IJ_F0OI)+ENRGP*POICE
C**** CALL SUBROUTINE FOR CALCULATION OF PRECIPITATION OVER SEA ICE

        CALL PREC_SI(SNOW,MSI2,HSIL,TSIL,SSIL,PRCP,ENRGP,RUN0,SRUN0,
#ifdef TRACERS_WATER
     *       TRSIL,TRPRCP,TRUN0,
#endif
     *       WETSNOW)

C**** check for snowice formation
c        if (RHOI*(SNOW+ACE1I+MSI2).gt.(ACE1I+MSI2)*1030.) then
c          print*,"snowice possibility",i,j,SNOW,MSI2+ACE1I,(SNOW+ACE1I
c     *         +MSI2)/1030.-(ACE1I+MSI2)/RHOI
c        end if

        SNOWI(I,J)  =SNOW
        RUNPSI(I,J) =RUN0
        SRUNPSI(I,J)=SRUN0
        HSI(:,I,J)=HSIL(:)
        SSI(:,I,J)=SSIL(:)
#ifdef TRACERS_WATER
        TRSI(:,:,I,J)=TRSIL(:,:)
        TRUNPSI(:,I,J)=TRUN0(:) ! tracer in runoff
#endif
        FLAG_DSWS(I,J)=FLAG_DSWS(I,J).or.WETSNOW
C**** reset flag if there was fresh snow (i.e. prcp but no rain!)
        IF (.not. WETSNOW .and. PRCP.gt.0.) FLAG_DSWS(I,J)=.FALSE.
C**** pond_melt accumulates in melt season only
        if ((J.gt.JM/2 .and. (jday.ge.152 .and. jday.lt.212)) .or.
     *       (J.lt.JM/2 .and. (jday.ge.334 .or. jday.lt.31))) then
          pond_melt(i,j)=pond_melt(i,j)+0.3d0*RUN0
        end if

C**** set gtemp array
        MSI(I,J)=MSI2
        GTEMP(1:2,2,I,J)=TSIL(1:2)
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J) = TRSIL(:,1)/(XSI(1)*(SNOW+ACE1I)-SSIL(1))
#endif

C**** Accumulate diagnostics for ice fraction
        AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)-RUN0 *POICE
        AJ(J,J_SMELT,ITYPE)=AJ(J,J_SMELT,ITYPE)-SRUN0*POICE
c       AJ(J,J_HMELT,ITYPE)=AJ(J,J_HMELT,ITYPE)-ERUN *POICE  ! ==0
C**** Accumulate diagnostics for water fraction
        AJ(J,J_IMELT,ITYPEO)=AJ(J,J_IMELT,ITYPEO)+RUN0 *POICE
        AJ(J,J_SMELT,ITYPEO)=AJ(J,J_SMELT,ITYPEO)+SRUN0*POICE
c       AJ(J,J_HMELT,ITYPEO)=AJ(J,J_HMELT,ITYPEO)+ERUN *POICE ! ==0
C**** Accumulate regional diagnostics
        AREG(JR,J_IMELT)=AREG(JR,J_IMELT)-RUN0 *POICE*DXYP(J)
        AREG(JR,J_SMELT)=AREG(JR,J_SMELT)-SRUN0*POICE*DXYP(J)
c       AREG(JR,J_HMELT)=AREG(JR,J_HMELT)-ERUN *POICE*DXYP(J) ! ==0

      END IF
      END DO
      END DO
C****
      END SUBROUTINE PRECIP_SI

      SUBROUTINE UNDERICE
!@sum  underice calculates basal fluxes under sea and lake ice
!@+    saves the resulting fluxes
!@auth Gavin Schmidt
!@ver  1.0
!@calls iceocean,icelake
      USE CONSTANT, only : rhow,lhm,omega,rhoi,byshi,shw
      USE MODEL_COM, only : im,jm,focean,dtsrc,qcheck,kocean
      USE GEOM, only : sinp,imaxj,dxyp
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm, trname
#endif
      USE SEAICE, only : lmi,xsi,icelake,iceocean,ac2oim,alami,alpha
     *     ,tfrez
      USE SEAICE_COM, only : msi,hsi,ssi,rsi
#ifdef TRACERS_WATER
     *     ,trsi
#endif
      USE FLUXES, only : fmsi_io,fhsi_io,fssi_io,ui2rho,gtemp,sss,mlhc
#ifdef TRACERS_WATER
     *     ,ftrsi_io,gtracer
#endif
      USE LAKES_COM, only : tlake,mwl,flake,gml,mldlk
#ifdef TRACERS_WATER
     *     ,trlake
#endif
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 coriol,ustar,Tm,Sm,Si,Ti,dh,mflux,hflux,sflux,fluxlim
     *     ,mlsh   !,mfluxmax
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: Trm,Tri,trflux,tralpha
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
#endif
#endif

      DO J=1,JM
C**** Coriolis parameter (on tracer grid)
        coriol = ABS(2.*OMEGA*SINP(J))
        DO I=1,IMAXJ(J)
          IF ((FOCEAN(I,J)+FLAKE(I,J))*RSI(I,J).gt.0) THEN
C**** Set mixed layer conditions
            Tm = GTEMP(1,1,I,J)
#ifdef TRACERS_WATER
            Trm(:)=GTRACER(:,1,I,J)
            Tri(:)=TRSI(:,LMI,I,J)/(XSI(LMI)*MSI(I,J)-SSI(LMI,I,J))
#ifdef TRACERS_SPECIAL_O18
            tralpha(:)=fracls(trname(:))
#else
            tralpha(:)=1.
#endif
#endif
            dh = 0.5*(XSI(LMI)*MSI(I,J))/RHOI
c            mfluxmax = (MSI(I,J)-AC2OIM)/dtsrc
            IF (FOCEAN(I,J).gt.0) THEN
C**** Ice lowest layer conditions
              Si = 1d3*SSI(LMI,I,J)/(XSI(LMI)*MSI(I,J))
              Ti = (HSI(LMI,I,J)/(XSI(LMI)*MSI(I,J))+LHM)*BYSHI
C**** for the Salinity thermodynamics case (or something similar)
c             Ti = TICE(HSI(LMI,I,J),SSI(LMI,I,J),XSI(LMI)*MSI(I,J))
              IF (KOCEAN.eq.1) THEN
C**** should we calculate ocean rho(Tm,Sm) here?
                Ustar = MAX(5d-4,SQRT(UI2rho(I,J)/RHOW))
                Sm = SSS(I,J)
                mlsh = MLHC(I,J)
                call iceocean(Ti,Si,Tm,Sm,dh,Ustar,Coriol,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *               Tri,Trm,trflux,tralpha,
#endif
     *               mflux,sflux,hflux)
              ELSE ! for fixed SST assume freezing temp at base,implicit
                hflux=alami*(Ti-tfrez(sss(i,j)))/(dh+alpha*dtsrc*alami
     *               *byshi/(XSI(LMI)*MSI(I,J)))
                mflux=0.
                sflux=0.
#ifdef TRACERS_WATER
                trflux= 0.
#endif
              END IF
            ELSE   ! for lakes (no salinity so solve directly)
              Ti = (HSI(LMI,I,J)/(XSI(LMI)*MSI(I,J))+LHM)*BYSHI
              mlsh=SHW*MLDLK(I,J)*RHOW
              sflux = 0.
              call icelake(Ti,Tm,dh,dtsrc,mlsh,
#ifdef TRACERS_WATER
     *             Tri,Trm,trflux,tralpha,
#endif
     *             mflux,hflux)
C**** Limit lake-to-ice flux if lake is too shallow (< 40cm)
              IF (MWL(I,J).lt.0.4d0*RHOW*FLAKE(I,J)*DXYP(J)) THEN
                FLUXLIM=-GML(I,J)/(DTSRC*FLAKE(I,J)*DXYP(J))
                IF (hflux.lt.FLUXLIM) hflux = FLUXLIM
                if (mflux.lt.0) then
                  mflux = 0.
#ifdef TRACERS_WATER
                  trflux= 0.
#endif
                end if
                if (qcheck) print*,"Flux limiting",I,J,MWL(I,J)/
     *               (RHOW*FLAKE(I,J)*DXYP(J)),FLUXLIM*DTSRC
              END IF
            END IF
            FMSI_IO(I,J) = mflux*dtsrc   ! positive down
            FHSI_IO(I,J) = hflux*dtsrc
            FSSI_IO(I,J) = sflux*dtsrc
#ifdef TRACERS_WATER
            FTRSI_IO(:,I,J)=trflux(:)*dtsrc
#endif
          ELSE
            FMSI_IO(I,J) = 0.
            FHSI_IO(I,J) = 0.
            FSSI_IO(I,J) = 0.
#ifdef TRACERS_WATER
            FTRSI_IO(:,I,J)=0.
#endif
          END IF
        END DO
      END DO

C****
      RETURN
      END SUBROUTINE UNDERICE

      SUBROUTINE MELT_SI
!@sum  MELT_SI driver for lateral melt of sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
!@calls SEAICE:SIMELT
      USE CONSTANT, only : sday
      USE MODEL_COM, only : im,jm,kocean,focean,itoice,itlkice,jhour
      USE GEOM, only : dxyp,imaxj
      USE DAGCOM, only : aj,j_imelt,j_hmelt,j_smelt,areg,jreg
      USE SEAICE, only : simelt,tfrez
      USE SEAICE_COM, only : rsi,hsi,msi,lmi,snowi,ssi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE LAKES_COM, only : flake
      USE FLUXES, only : sss,flowo,eflowo,sflowo,gtemp
#ifdef TRACERS_WATER
     *     ,trflowo
#endif
      IMPLICIT NONE
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 MSI2,ROICE,SNOW,ENRGW,ENRGUSED,ANGLE,RUN0,SALT,POCEAN,TFO
     *     ,PWATER,Tm,DT
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: TRSIL
      REAL*8, DIMENSION(NTM) :: TRUN0
#endif
      INTEGER I,J,ITYPE,JR

C**** CALCULATE LATERAL MELT ONCE A DAY (ALSO ELIMINATE SMALL AMOUNTS)
C**** We could put this in daily but it then we need an extra routine to
C**** add fluxes to oceans/lakes. This is a little more confusing, but
C**** we can use the S/E/FLOWO arrays to pass the information.
      IF (Jhour.eq.0) THEN
        DT=SDAY    ! if called more frequently this should change
        DO J=1,JM
        DO I=1,IMAXJ(J)
          PWATER=FOCEAN(I,J)+FLAKE(I,J)
          POCEAN=FOCEAN(I,J)
          RUN0=0. ; ENRGUSED=0. ; SALT=0.
#ifdef TRACERS_WATER
          TRUN0(:) = 0.
#endif
C**** Call simelt if (lake and v. small ice) or (q-flux ocean, some ice)
          IF ( (RSI(I,J).lt.1d-4 .and. FLAKE(I,J)*RSI(I,J).gt.0)
     *         .or. (KOCEAN.eq.1.and.POCEAN*RSI(I,J).gt.0) ) THEN
            JR=JREG(I,J)
            IF (POCEAN.gt.0) THEN
              ITYPE =ITOICE
              TFO = tfrez(sss(i,j))
            ELSE
              ITYPE =ITLKICE
              TFO = 0.
            END IF
            Tm=GTEMP(1,1,I,J)
            ROICE=RSI(I,J)
            MSI2=MSI(I,J)
            SNOW=SNOWI(I,J)     ! snow mass
            HSIL(:)= HSI(:,I,J) ! sea ice enthalpy
            SSIL(:)= SSI(:,I,J) ! sea ice salt
#ifdef TRACERS_WATER
            TRSIL(:,:)=TRSI(:,:,I,J) ! tracer content of sea ice
#endif
            CALL SIMELT(DT,ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,Tm,TFO,TSIL
#ifdef TRACERS_WATER
     *           ,TRSIL,TRUN0
#endif
     *           ,ENRGUSED,RUN0,SALT)

C**** accumulate diagnostics
            AJ(J,J_HMELT,ITYPE)=AJ(J,J_HMELT,ITYPE)-ENRGUSED
            AJ(J,J_SMELT,ITYPE)=AJ(J,J_SMELT,ITYPE)+    SALT
            AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)+    RUN0
            AREG(JR,J_HMELT)=AREG(JR,J_HMELT)-ENRGUSED*DXYP(J)
            AREG(JR,J_SMELT)=AREG(JR,J_SMELT)+    SALT*DXYP(J)
            AREG(JR,J_IMELT)=AREG(JR,J_IMELT)+    RUN0*DXYP(J)
C**** Update prognostic sea ice variables
            RSI(I,J)=ROICE
            MSI(I,J)=MSI2
            SNOWI(I,J)=SNOW
            HSI(:,I,J)=HSIL(:)
            SSI(:,I,J)=SSIL(:)
#ifdef TRACERS_WATER
            TRSI(:,:,I,J)=TRSIL(:,:)
#endif
          END IF
C**** Save fluxes (in kg, J etc.), positive into ocean
          FLOWO(I,J) = RUN0*PWATER*DXYP(J)
          EFLOWO(I,J)=-ENRGUSED*PWATER*DXYP(J)
          SFLOWO(I,J)= SALT*PWATER*DXYP(J)
#ifdef TRACERS_WATER
          TRFLOWO(:,I,J)=TRUN0(:)*PWATER*DXYP(J)
#endif
        END DO
        END DO
      ELSE
        FLOWO=0. ; EFLOWO=0. ; SFLOWO=0.
#ifdef TRACERS_WATER
        TRFLOWO=0.
#endif
      END IF
C****
      RETURN
      END SUBROUTINE MELT_SI

      SUBROUTINE GROUND_SI
!@sum  GROUND_SI driver for applying surface + base fluxes to sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
!@calls SEAICE:SEA_ICE
      USE CONSTANT, only : lhm,byshi,rhow
      USE MODEL_COM, only : im,jm,dtsrc,fland,kocean,focean
     *     ,itoice,itlkice,jday
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : e0,e1,evapor,runosi,erunosi,srunosi,solar
     *     ,fmsi_io,fhsi_io,fssi_io
#ifdef TRACERS_WATER
     *     ,ftrsi_io,trevapor,trunosi
#endif
      USE SEAICE, only : sea_ice,ssidec,lmi,xsi,ace1i,qsfix
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,pond_melt,flag_dsws
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE LAKES_COM, only : mwl,gml,flake
      USE DAGCOM, only : aj,areg,aij,jreg,ij_erun2,ij_rsoi,ij_msi2
     *     ,j_imelt,j_hmelt,j_smelt,j_rsnow,ij_rsit,ij_rsnw,ij_snow
     *     ,ij_mltp
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,SSIL
      REAL*8 SNOW,ROICE,MSI2,F0DT,F1DT,EVAP,SROX(2)
     *     ,FMOC,FHOC,FSOC,POICE,PWATER,SCOVI
      REAL*8 MSI1,MFLUX,HFLUX,SFLUX,RUN,ERUN,SRUN,MELT12
      INTEGER I,J,JR,ITYPE
      LOGICAL WETSNOW
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: trsil
      REAL*8, DIMENSION(NTM) :: trflux,ftroc,trevap,trrun
#endif

      DO J=1,JM
      DO I=1,IMAXJ(J)
      PWATER=FOCEAN(I,J)+FLAKE(I,J)   ! 1.-FLAND(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      JR=JREG(I,J)
      SOLAR(3,I,J)=0
      RUNOSI(I,J)=0
      ERUNOSI(I,J)=0
      SRUNOSI(I,J)=0
#ifdef TRACERS_WATER
      TRUNOSI(:,I,J) = 0.
#endif
      IF (POICE.gt.0) THEN

        F0DT=E0(I,J,2) ! heat flux to the top ice surface (J/m^2)
        F1DT=E1(I,J,2) ! heat flux between 1st and 2nd ice layer (J/m^2)
        EVAP=EVAPOR(I,J,2) ! evaporation/dew at the ice surface (kg/m^2)
        SROX(1)=SOLAR(2,I,J) ! solar radiation absrbd by sea ice (J/m^2)
        FMOC=fmsi_io(i,j)  ! mass flux at base (kg/m^2)
        FSOC=fssi_io(i,j)  ! salt flux at base (kg/m^2)
        FHOC=fhsi_io(i,j)  ! heat flux at base (J/M^2)
        SNOW= SNOWI(I,J)  ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)  ! sea ice enthalpy
        SSIL(:) = SSI(:,I,J)  ! sea ice salt
        WETSNOW=FLAG_DSWS(I,J)  ! wetness of snow
#ifdef TRACERS_WATER
        TREVAP(:) = TREVAPOR(:,2,I,J)
        FTROC(:)  = ftrsi_io(:,i,j)
        TRSIL(:,:)= TRSI(:,:,I,J)
#endif
        IF (FOCEAN(I,J).gt.0) THEN
          ITYPE=ITOICE
        ELSE
          ITYPE=ITLKICE
        END IF

        AIJ(I,J,IJ_RSOI) =AIJ(I,J,IJ_RSOI) +POICE
        AIJ(I,J,IJ_MSI2) =AIJ(I,J,IJ_MSI2) +MSI2*POICE

        CALL SEA_ICE(DTSRC,SNOW,ROICE,HSIL,SSIL,MSI2,F0DT,F1DT,EVAP,SROX
#ifdef TRACERS_WATER
     *       ,TRSIL,TREVAP,FTROC,TRRUN
#endif
     *       ,FMOC,FHOC,FSOC,RUN,ERUN,SRUN,WETSNOW,MELT12)

C**** Decay sea ice salinity
        MSI1 = ACE1I + SNOW
        if (.not. qsfix .and. FOCEAN(I,J).gt.0) then
          CALL SSIDEC(I,J,MSI1,MSI2,HSIL,SSIL,DTsrc,
#ifdef TRACERS_WATER
     *         TRSIL,TRFLUX,
#endif
     *         MFLUX,HFLUX,SFLUX)
        else
          MFLUX=0. ; SFLUX=0. ; HFLUX=0.
#ifdef TRACERS_WATER
          TRFLUX = 0.
#endif
        end if

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J)=SNOW
        HSI(:,I,J)=HSIL(:)
        SSI(:,I,J)=SSIL(:)
        MSI(I,J) = MSI2
#ifdef TRACERS_WATER
        TRSI(:,:,I,J) = TRSIL(:,:)
#endif
        FLAG_DSWS(I,J)=WETSNOW
C**** pond_melt accumulates in melt season only
        if ((J.gt.JM/2 .and. (jday.ge.152 .and. jday.lt.212)) .or.
     *       (J.lt.JM/2 .and. (jday.ge.334 .or. jday.lt.31))) then
          pond_melt(i,j)=pond_melt(i,j)+0.3d0*MELT12
          pond_melt(i,j)=MIN(pond_melt(i,j),0.5*(MSI2+SNOW+ACE1I))
        end if

C**** Net fluxes to ocean
        RUNOSI(I,J) = FMOC + RUN  + MFLUX
        ERUNOSI(I,J)= FHOC + ERUN + HFLUX
        SRUNOSI(I,J)= FSOC + SRUN + SFLUX
        SOLAR(3,I,J)= SROX(2)
#ifdef TRACERS_WATER
        TRUNOSI(:,I,J) = FTROC(:) + TRRUN(:) + TRFLUX(:)
#endif

C**** ACCUMULATE DIAGNOSTICS
          SCOVI=0.
          IF (SNOWI(I,J).GT.0) SCOVI=POICE

          AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVI
          AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*POICE
          AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+POICE
          AIJ(I,J,IJ_MLTP)=AIJ(I,J,IJ_MLTP)+pond_melt(i,j)*POICE

          AJ(J,J_RSNOW,ITYPE)=AJ(J,J_RSNOW,ITYPE)+SCOVI
          AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)+(FMOC+RUN+MFLUX)*POICE
          AJ(J,J_HMELT,ITYPE)=AJ(J,J_HMELT,ITYPE)+(FHOC+ERUN+HFLUX)
     *         *POICE
          AJ(J,J_SMELT,ITYPE)=AJ(J,J_SMELT,ITYPE)+(FSOC+SRUN+SFLUX)
     *         *POICE

          AREG(JR,J_RSNOW)=AREG(JR,J_RSNOW)+SCOVI*DXYP(J)
          AREG(JR,J_IMELT)=AREG(JR,J_IMELT)+(FMOC+RUN+MFLUX)*POICE
     *         *DXYP(J)
          AREG(JR,J_HMELT)=AREG(JR,J_HMELT)+(FHOC+ERUN+HFLUX)*POICE
     *         *DXYP(J)
          AREG(JR,J_SMELT)=AREG(JR,J_SMELT)+(FSOC+SRUN+SFLUX)*POICE
     *         *DXYP(J)

      END IF
      END DO
      END DO
C****
      END SUBROUTINE GROUND_SI

      SUBROUTINE FORM_SI
!@sum  FORM_SI driver for adding new sea ice
!@auth Original Development team
!@ver  1.0
!@calls seaice:addice
      USE CONSTANT, only : lhm,byshi
      USE MODEL_COM, only : im,jm,focean,kocean,fland
     *     ,itocean,itoice,itlake,itlkice,itime
      USE GEOM, only : imaxj,dxyp
#ifdef TRACERS_WATER
      USE TRACER_COM, only : itime_tr0,tr_wd_type,nWater
#endif
      USE DAGCOM, only : aj,areg,aij,jreg,j_rsi,j_ace1,j_ace2,j_snow
     *     ,j_smelt,j_imelt,j_hmelt,ij_tsi,ij_ssi1,ij_ssi2
      USE SEAICE, only : ace1i,addice,lmi,fleadoc,fleadlk,xsi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
      USE TRACER_DIAG_COM, only : tij_seaice,taijn
#endif
      USE FLUXES, only : dmsi,dhsi,dssi,gtemp
#ifdef TRACERS_WATER
     *     ,dtrsi,gtracer
#endif
      USE LAKES_COM, only : flake
      IMPLICIT NONE

      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 SNOW,ROICE,MSI2,ENRGFO,ACEFO,ACE2F,ENRGFI,SALTO,SALTI
     *     ,POICE,PWATER,FLEAD,POCEAN
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,LMI) :: trsil
      REAL*8, DIMENSION(NTM) :: tro,tri
#endif
      LOGICAL QFIXR
      INTEGER I,J,JR,ITYPE,ITYPEO,N

      DO J=1,JM
      DO I=1,IMAXJ(J)
      PWATER=FOCEAN(I,J)+FLAKE(I,J)
      ROICE=RSI(I,J)
      POICE=ROICE*PWATER
      POCEAN=(1.-ROICE)*PWATER
      JR=JREG(I,J)
      IF (PWATER.gt.0) THEN

        SNOW= SNOWI(I,J)      ! snow mass (kg/m^2)
        MSI2= MSI(I,J)
        HSIL(:) = HSI(:,I,J)      ! sea ice enthalpy
        SSIL(:) = SSI(:,I,J)      ! sea ice salt
#ifdef TRACERS_WATER
        TRSIL(:,:)= TRSI(:,:,I,J)
#endif

        IF (FOCEAN(I,J).gt.0) THEN
          FLEAD=FLEADOC
          ITYPE=ITOICE
          ITYPEO=ITOCEAN
          IF (KOCEAN.eq.1) THEN
            QFIXR=.FALSE.
          ELSE
            QFIXR=.TRUE.
          END IF
        ELSE
          FLEAD=FLEADLK
          ITYPE=ITLKICE
          ITYPEO=ITLAKE
          QFIXR=.FALSE.
        END IF

        ACEFO=DMSI(1,I,J)
        ACE2F=DMSI(2,I,J)
        ENRGFO=DHSI(1,I,J)
        ENRGFI=DHSI(2,I,J)
        SALTO=DSSI(1,I,J)
        SALTI=DSSI(2,I,J)
#ifdef TRACERS_WATER
        TRO(:) = DTRSI(:,1,I,J)
        TRI(:) = DTRSI(:,2,I,J)
#endif
C**** ice formation diagnostics on the atmospheric grid
C**** open ocean diagnostics
        AJ(J,J_SMELT,ITYPEO)=AJ(J,J_SMELT,ITYPEO)-SALTO *POCEAN
        AJ(J,J_HMELT,ITYPEO)=AJ(J,J_HMELT,ITYPEO)-ENRGFO*POCEAN
        AJ(J,J_IMELT,ITYPEO)=AJ(J,J_IMELT,ITYPEO)-ACEFO *POCEAN
C**** Ice-covered ocean diagnostics
        AJ(J,J_SMELT,ITYPE)=AJ(J,J_SMELT,ITYPE)-SALTI *POICE
        AJ(J,J_HMELT,ITYPE)=AJ(J,J_HMELT,ITYPE)-ENRGFI*POICE
        AJ(J,J_IMELT,ITYPE)=AJ(J,J_IMELT,ITYPE)-ACE2F *POICE
C**** regional diagnostics
        AREG(JR,J_IMELT)=AREG(JR,J_IMELT)-
     *       (ACEFO *POCEAN+ACE2F *POICE)*DXYP(J)
        AREG(JR,J_HMELT)=AREG(JR,J_HMELT)-
     *       (ENRGFO*POCEAN+ENRGFI*POICE)*DXYP(J)
        AREG(JR,J_SMELT)=AREG(JR,J_SMELT)-
     *       (SALTO *POCEAN+SALTI *POICE)*DXYP(J)

        CALL ADDICE (SNOW,ROICE,HSIL,SSIL,MSI2,TSIL,ENRGFO,ACEFO,ACE2F,
     *       ENRGFI,SALTO,SALTI,
#ifdef TRACERS_WATER
     *       TRSIL,TRO,TRI,
#endif
     *       FLEAD,QFIXR)

C**** RESAVE PROGNOSTIC QUANTITIES
        SNOWI(I,J) = SNOW
        HSI(:,I,J) = HSIL(:)
        SSI(:,I,J) = SSIL(:)
#ifdef TRACERS_WATER
        TRSI(:,:,I,J) = TRSIL(:,:)
#endif
        IF (.not. QFIXR) THEN
          RSI(I,J)=ROICE
          MSI(I,J)=MSI2
        END IF
C**** set gtemp array
        GTEMP(1:2,2,I,J)=TSIL(1:2)
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J)=TRSIL(:,1)/(XSI(1)*(SNOW+ACE1I)-SSIL(1))
#endif

C**** ACCUMULATE DIAGNOSTICS
        AJ(J,J_RSI, ITYPE)=AJ(J,J_RSI, ITYPE)+        POICE
        AJ(J,J_ACE1,ITYPE)=AJ(J,J_ACE1,ITYPE)+ACE1I  *POICE
        AJ(J,J_ACE2,ITYPE)=AJ(J,J_ACE2,ITYPE)+MSI2   *POICE
        AJ(J,J_SNOW,ITYPE)=AJ(J,J_SNOW,ITYPE)+SNOW   *POICE
        IF (JR.ne.24) THEN
          AREG(JR,J_RSI) =AREG(JR,J_RSI) +        POICE*DXYP(J)
          AREG(JR,J_SNOW)=AREG(JR,J_SNOW)+SNOW   *POICE*DXYP(J)
          AREG(JR,J_ACE1)=AREG(JR,J_ACE1)+ACE1I  *POICE*DXYP(J)
          AREG(JR,J_ACE2)=AREG(JR,J_ACE2)+MSI2   *POICE*DXYP(J)
        END IF
        AIJ(I,J,IJ_TSI)=AIJ(I,J,IJ_TSI)+
     *       POICE*(XSI(3)*TSIL(3)+XSI(4)*TSIL(4))
        AIJ(I,J,IJ_SSI1)=AIJ(I,J,IJ_SSI1)+POICE*(SSIL(1)+SSIL(2))/ACE1I
        AIJ(I,J,IJ_SSI2)=AIJ(I,J,IJ_SSI2)+POICE*(SSIL(3)+SSIL(4))
     *       /MSI(I,J)

#ifdef TRACERS_WATER
C**** Save sea ice tracer amount
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. tr_wd_TYPE(n).eq.nWater) then
          taijn(i,j,tij_seaice,n)=taijn(i,j,tij_seaice,n)+
     *         rsi(i,j)*sum(trsil(n,:))
        end if
      end do
#endif
      END IF
      END DO
      END DO
C**** replicate ice values at the north pole
      DO I=2,IM
        RSI(I,JM)=RSI(1,JM)
        MSI(I,JM)=MSI(1,JM)
        HSI(:,I,JM)=HSI(:,1,JM)
        SSI(:,I,JM)=SSI(:,1,JM)
        SNOWI(I,JM)=SNOWI(1,JM)
        GTEMP(1:2,2,I,JM)=GTEMP(1:2,2,1,JM)
#ifdef TRACERS_WATER
        TRSI(:,:,I,JM) = TRSI(:,:,1,JM)
        GTRACER(:,2,I,JM) = GTRACER(:,2,1,JM)
#endif
      END DO
C****
      END SUBROUTINE FORM_SI

      SUBROUTINE vflx_OCEAN
!@sum  vflx_OCEAN saves quantities for OHT calculations
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,focean
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
C****       2  MSI2   (INSTANTANEOUS AT NOON GMT)
C****       3  HSIT   (INSTANTANEOUS AT NOON GMT)
C****
      DO J=1,JM
        DO I=1,IM
          IF (FOCEAN(I,J).gt.0) THEN
            OA(I,J,1)=ACE1I+SNOWI(I,J)
            OA(I,J,2)=MSI(I,J)
            OA(I,J,3)=SUM(HSI(:,I,J))
          END IF
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE vflx_OCEAN

      SUBROUTINE init_ice(iniOCEAN)
!@sum  init_ice initialises ice arrays
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : byshi,lhm,shi,rhow
      USE MODEL_COM, only : im,jm,kocean,focean
      USE SEAICE, only : xsi,ace1i,ac2oim,ssi0,tfrez,oi_ustar0,silmfac
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi,pond_melt,flag_dsws
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,sss,ui2rho
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DAGCOM, only : npts,icon_MSI,icon_HSI,icon_SSI,title_con
     *     ,conpt0
      USE PARAM
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE. , iniOCEAN
      CHARACTER CONPT(NPTS)*10
      INTEGER I,J
      REAL*8 MSI1,TFO

C**** set up a default ice-ocean stress field. This can be changed by
C**** adjusting oi_ustar0 in the parameter list. If ice dynamics
C**** is used, this is overwritten.
      call sync_param("oi_ustar0",oi_ustar0)
      UI2rho = rhow*(oi_ustar0)**2

C**** Adjust degree of lateral melt by changing silmfac
C**** Default is 2.5d-8, but could be changed by a factor of 2.
      call sync_param("silmfac",silmfac)

      IF (KOCEAN.EQ.0.and.iniOCEAN) THEN
C****   set defaults for no ice case
        DO J=1,JM
        DO I=1,IM
          IF (RSI(I,J).le.0) THEN
            MSI1        =ACE1I
            MSI(I,J)    =AC2OIM
            SNOWI(I,J)  =0.
            IF (FOCEAN(I,J).gt.0) THEN
              SSI(1:2,I,J)=SSI0*XSI(1:2)*ACE1I
              SSI(3:4,I,J)=SSI0*XSI(3:4)*AC2OIM
              TFO = TFREZ(SSS(I,J))
            ELSE
              SSI(:,I,J)  = 0.
              TFO = 0.
            END IF
            HSI(1:2,I,J)=(SHI*TFO-LHM)*XSI(1:2)*ACE1I
            HSI(3:4,I,J)=(SHI*TFO-LHM)*XSI(3:4)*AC2OIM
#ifdef TRACERS_WATER
            TRSI(:,:,I,J)=0.
#endif
            pond_melt(i,j) = 0.
            flag_dsws(i,j) = .FALSE.
          END IF
        END DO
        END DO
      END IF
C**** set GTEMP array for ice
      DO J=1,JM
      DO I=1,IM
        MSI1=SNOWI(I,J)+ACE1I
        GTEMP(1:2,2,I,J)=(HSI(1:2,I,J)/(XSI(1:2)*MSI1)+LHM)*BYSHI
#ifdef TRACERS_WATER
        GTRACER(:,2,I,J) = TRSI(:,1,I,J)/(XSI(1)*MSI1-SSI(1,I,J))
#endif
      END DO
      END DO

C**** Set conservation diagnostics for ice mass, energy, salt
      CONPT=CONPT0
      CONPT(3)="LAT. MELT" ; CONPT(4)="PRECIP"
      CONPT(5)="THERMO"    ; CONPT(6)="ADVECT"
      CONPT(8)="OCN FORM"
      QCON=(/ F, F, T, T, T, T, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"ICE MASS","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_MSI)
      QCON=(/ F, F, T, T, T, T, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"ICE ENRG","(10**6 J/M^2)   ",
     *     "(10**-3 J/S M^2)",1d-6,1d3,icon_HSI)
      QCON=(/ F, F, T, T, T, T, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT,"ICE SALT","(10**-3 KG/M^2) ",
     *     "(10**-9 KG/SM^2)",1d3,1d9,icon_SSI)
C****
      END SUBROUTINE init_ice

      SUBROUTINE conserv_MSI(ICE)
!@sum  conserv_MSI calculates total amount of snow and ice over water
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,fland
      USE GEOM, only : imaxj
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
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
      END SUBROUTINE conserv_MSI

      SUBROUTINE conserv_HSI(EICE)
!@sum  conserv_HSI calculates total ice energy over water
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,fland
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
      END SUBROUTINE conserv_HSI

      SUBROUTINE conserv_SSI(SALT)
!@sum  conserv_SSI calculates total amount of salt in ocean ice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,focean
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : rsi,ssi,lmi
      IMPLICIT NONE
!@var SALT total salt in ocean ice (kg/m^2)
      REAL*8, DIMENSION(JM) :: SALT
      INTEGER I,J,L

      DO J=1,JM
        SALT(J)=0
        DO I=1,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
            DO L=1,LMI
              SALT(J)=SALT(J)+FOCEAN(I,J)*RSI(I,J)*SSI(L,I,J)
            END DO
          END IF
        END DO
      END DO
      SALT(1) =FIM*SALT(1)
      SALT(JM)=FIM*SALT(JM)
      RETURN
C****
      END SUBROUTINE conserv_SSI

      SUBROUTINE daily_ice
!@sum daily_ice performs ice processes that are needed everyday
!@auth Gavin Schmidt
      USE MODEL_COM, only : jm,jday
      USE GEOM, only : imaxj
      USE SEAICE_COM, only : flag_dsws,pond_melt
      IMPLICIT NONE
      INTEGER I,J

C**** Every day adjust pond-melt
      DO J=1,JM
      DO I=1,IMAXJ(J)

C**** pond_melt decreases linearly in shoulder season
        if (J.gt.JM/2 .and. (jday.ge.212 .and. jday.lt.244)) then
          pond_melt(i,j)=pond_melt(i,j)*(1.-1./(244.-jday))
        elseif (J.lt.JM/2 .and. (jday.ge.31 .and. jday.lt.60)) then
          pond_melt(i,j)=pond_melt(i,j)*(1.-1./(60.-jday))
C**** allow pond_melt to accumulate in melt season only
C**** otherwise pond_melt is zero
        elseif (.not.((J.gt.JM/2 .and. (jday.ge.152 .and. jday.lt.212))
     *         .or. (J.lt.JM/2 .and. (jday.ge.334 .or. jday.lt.31))) )
     *         then
          pond_melt(i,j)=0.
        end if
      END DO
      END DO

      RETURN
      END SUBROUTINE daily_ice
