#include "rundeck_opts.h"

!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont init_LI,PRECIP_LI,GROUND_LI

      SUBROUTINE init_LI
!@sum  init_ice initialises landice arrays
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : edpery,sday,lhm
      USE MODEL_COM, only : im,jm,flice,focean,dtsrc
      USE GEOM, only : dxyp
      USE LANDICE, only: ace1li,ace2li,glmelt_on
      USE LANDICE_COM, only : tlandi,snowli
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trli0
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
C**** The net accumulation from IPCC2 report is 2016x10**12 kg/year
C**** for Antarctica and for Greenland it is 316x10**12 kg/year
C****
!@var ACCPDA total accumulation per year for Antarctica (kg/yr)
!@var ACCPDG total accumulation per year for Greenland (kg/yr)
      REAL*8, PARAMETER :: ACCPDA = 2016d12, ACCPDG = 316d12

      INTEGER, PARAMETER :: NBOXMAX=18
      INTEGER IFW(NBOXMAX),JFW(NBOXMAX)
      INTEGER :: JML, JMU, IML1, IMU1, IML2, IMU2, NBOX
      REAL*8 ACCPCA,ACCPCG,FWAREA
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
            IF (SNOWLI(I,J).gt.1d-5) THEN
              GTRACER(:,3,I,J)=TRSNOWLI(:,I,J)/SNOWLI(I,J)
            ELSE
              GTRACER(:,3,I,J)=trli0(:) !TRLNDI(:,I,J)/(ACE1LI+ACE2LI)
            END IF
#endif
          END IF
        END DO
      END DO

C**** Calculate (fixed) iceberg melt terms from Antarctica and Greenland
      call sync_param("glmelt_on",glmelt_on)

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

      GMELT = 0. ; EGMELT = 0.
#ifdef TRACERS_OCEAN
      TRGMELT = 0.
#endif

      if (do_glmelt .and. glmelt_on.eq.1) then
C****  Note that water goes in with as if it is ice at 0 deg.
C****  Possibly this should be a function of in-situ freezing temp?
C****
C**** Antarctica
! accumulation (kg per source time step) per water column
C**** integrate area (which will depend on resolution/landmask)
      FWAREA_part=0.
      DO J=MAX(J_0,JML),MIN(J_1,JMU)
        DO I=1,IM
          IF (FOCEAN(I,J).GT.0.) THEN
            IF ((I.GE.IML1.AND.I.LE.IMU1) .or. I.GE.IML2 .or. I.LE
     *           .IMU2) THEN
              FWAREA_part(J)=FWAREA_part(J)+DXYP(J)*FOCEAN(I,J)
            END IF
          END IF
        END DO
      END DO
      CALL GLOBALSUM(grid, FWAREA_part, FWAREA, all=.true.)

      ACCPCA = ACCPDA*DTsrc/(EDPERY*SDAY*FWAREA)      ! kg/m^2
      DO J=MAX(J_0,JML),MIN(J_1,JMU)
        DO I=1,IM
          IF (FOCEAN(I,J).GT.0.) THEN
            IF ((I.GE.IML1.AND.I.LE.IMU1) .or. I.GE.IML2 .or. I.LE
     *           .IMU2) THEN
              GMELT(I,J)  =  ACCPCA*DXYP(J)*FOCEAN(I,J)  ! kg
              EGMELT(I,J) = -LHM*ACCPCA *DXYP(J)*FOCEAN(I,J) ! J
#ifdef TRACERS_OCEAN
              TRGMELT(:,I,J)= trglac(:)*ACCPCA*DXYP(J)*FOCEAN(I,J)  ! kg
#endif
            END IF
          END IF
        END DO
      END DO

C**** Greenland
! accumulation (kg per source time step) per water column
C**** integrate area (which will depend on resolution/landmask)
      FWAREA_part=0.
      DO N=1,NBOX
        I=IFW(N)
        J=JFW(N)
        If (J >= J_0 .and. J <= J_1) 
     &       FWAREA_part(J)=FWAREA_part(J)+DXYP(J)*FOCEAN(I,J)      
      END DO
      CALL GLOBALSUM(grid, FWAREA_part, FWAREA, all=.true.)

      ACCPCG = ACCPDG*DTsrc/(EDPERY*SDAY*FWAREA)  ! kg/m^2 
      DO N=1,NBOX
        I=IFW(N)
        J=JFW(N)
        If (J >= J_0 .and. J <= J_1) THEN
          GMELT(I,J)  =  ACCPCG*DXYP(J)*FOCEAN(I,J) ! kg
          EGMELT(I,J) = -LHM*ACCPCG*DXYP(J)*FOCEAN(I,J) ! J
#ifdef TRACERS_OCEAN
        TRGMELT(:,I,J) = trglac(:)*ACCPCG*DXYP(J)*FOCEAN(I,J)  ! kg
#endif
        End If
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
      USE LANDICE_COM, only : snowli,tlandi
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,ntm,trli0
#endif
      USE DIAG_COM, only : aj=>aj_loc,areg,aij=>aij_loc
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
      REAL*8  :: AREG_SUM
      REAL*8, DIMENSION(
     &        size(AREG,1),grid%j_strt_halo:grid%j_stop_halo,4 )
     &        :: AREG_PART
C**** Get useful grid parameters
      INTEGER :: I, J, JR
      INTEGER :: J_0, J_1, J_0H, J_1H
      CALL GET(GRID,J_STRT=J_0      , J_STOP=J_1      ,
     &              J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )

C**** Initialize work array
      AREG_PART(:,J_0H:J_1H,1:4)=0.
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
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP

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
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
c       TRDIFS(:)     !  diagnostic?
        IF (SNOW.gt.1d-5) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=trli0(:)  !TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
c        AJ(J,J_TYPE, ITLANDI)=AJ(J,J_TYPE, ITLANDI)+      PLICE
        AJ(J,J_RUN,  ITLANDI)=AJ(J,J_RUN,  ITLANDI)+RUN0 *PLICE
C       AJ(J,J_ERUN ,ITLANDI)=AJ(J,J_ERUN ,ITLANDI)+ERUN0*PLICE ! (Tg=0)
        AJ(J,J_IMPLM,ITLANDI)=AJ(J,J_IMPLM,ITLANDI)+DIFS *PLICE
        AJ(J,J_IMPLH,ITLANDI)=AJ(J,J_IMPLH,ITLANDI)+ERUN2*PLICE
        AREG_PART(JR,J,1)=AREG_PART(JR,J,1) + RUN0 *PLICE*DXYPJ
c       AREG_PART(JR,J,2)=AREG_PART(JR,J,2)+ERUN0*PLICE*DXYPJ ! (Tg=0)
        AREG_PART(JR,J,3)=AREG_PART(JR,J,3)+DIFS *PLICE*DXYPJ
        AREG_PART(JR,J,4)=AREG_PART(JR,J,4)+ERUN2*PLICE*DXYPJ
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+ERUN2
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
      END IF
      END DO
      END DO
      CALL CHECKSUM(grid, SNOWLI,__LINE__,__FILE__)
      CALL CHECKSUM_COLUMN(grid, TLANDI,__LINE__,__FILE__)

C**** Finish summing and store total accumulations into AREG.
      DO JR=1,SIZE(AREG,1)
        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,1), AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_RUN  )=AREG(JR,J_RUN) + AREG_SUM

c       AREG_SUM=0.
c       CALL GLOBALSUM(GRID,AREG_PART(JR,:,2), AREG_SUM, ALL=.TRUE.)
c       AREG(JR,J_ERUN )=AREG(JR,J_ERUN) + AREG_SUM         ! (Tg=0)

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,3), AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM) + AREG_SUM 

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,4), AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH) + AREG_SUM 
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
      USE DIAG_COM, only : aj=>aj_loc,areg,aij=>aij_loc
     *     ,jreg,ij_runli,ij_f1li,ij_erun2
     *     ,j_wtr1,j_ace1,j_wtr2,j_ace2,j_snow,j_run
     *     ,j_implh,j_implm,j_rsnow,ij_rsnw,ij_rsit,ij_snow,ij_f0oc
     *     ,j_rvrd,j_ervr,ij_mrvr,ij_ervr
      USE LANDICE_COM, only : snowli,tlandi
#ifdef TRACERS_WATER
     *     ,ntm,trsnowli,trlndi,trli0
#endif
      USE FLUXES, only : e0,e1,evapor,gtemp,runoli,gmelt,egmelt
#ifdef TRACERS_WATER
     *     ,trunoli,trevapor,gtracer
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
#endif
#ifdef TRACERS_OCEAN
     *     ,trgmelt
      USE TRDIAG_COM, only : taijn,tij_rvr
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
      REAL*8  :: AREG_SUM
      REAL*8, DIMENSION(
     &        size(AREG,1),grid%j_strt_halo:grid%j_stop_halo,9 )
     &        :: AREG_PART

      INTEGER I,J,JR
      INTEGER :: J_0,J_1, J_0H, J_1H
      CALL GET(GRID,J_STRT=J_0      ,J_STOP=J_1 
     &             ,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

C**** Initialize work array
      AREG_PART(:,J_0H:J_1H,1:9) = 0.

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
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
c       TRDIFS(:)     !  diagnostic?
        IF (SNOW.gt.1d-5) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=trli0(:) !TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
        SCOVLI=0
        IF (SNOWLI(I,J).GT.0.) SCOVLI=PLICE
        AJ(J,J_RSNOW,ITLANDI)=AJ(J,J_RSNOW,ITLANDI)+SCOVLI
        AREG_PART(JR,J,1)=AREG_PART(JR,J,1)+SCOVLI*DXYPJ
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVLI
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*PLICE
        AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+PLICE

        AJ(J,J_RUN,ITLANDI)  =AJ(J,J_RUN,ITLANDI)  +RUN0 *PLICE
        AJ(J,J_SNOW,ITLANDI) =AJ(J,J_SNOW,ITLANDI) +SNOW *PLICE
        AJ(J,J_ACE1,ITLANDI) =AJ(J,J_ACE1,ITLANDI)+ACE1LI*PLICE
        AJ(J,J_ACE2,ITLANDI) =AJ(J,J_ACE2,ITLANDI)+ACE2LI*PLICE
        AJ(J,J_IMPLH,ITLANDI)=AJ(J,J_IMPLH,ITLANDI)+EDIFS*PLICE
        AJ(J,J_IMPLM,ITLANDI)=AJ(J,J_IMPLM,ITLANDI)+DIFS *PLICE

        AREG_PART(JR,J,2) =AREG_PART(JR,J,2) +RUN0  *PLICE*DXYPJ
        AREG_PART(JR,J,3) =AREG_PART(JR,J,3) +SNOW  *PLICE*DXYPJ
        AREG_PART(JR,J,4) =AREG_PART(JR,J,4) +ACE1LI*PLICE*DXYPJ
        AREG_PART(JR,J,5) =AREG_PART(JR,J,5) +ACE2LI*PLICE*DXYPJ
        AREG_PART(JR,J,6) =AREG_PART(JR,J,6) +EDIFS *PLICE*DXYPJ
        AREG_PART(JR,J,7) =AREG_PART(JR,J,7) +DIFS  *PLICE*DXYPJ

        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS+F1DT
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+EDIFS
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
      AIJ(I,J,IJ_F0OC) = AIJ(I,J,IJ_F0OC)        +(1.-RSI(I,J))
     *     *EGMELT(I,J)*BYDXYP(J)
      
      AREG_PART(JR,J,8)   = AREG_PART(JR,J,8) +  GMELT(I,J)
      AREG_PART(JR,J,9)   = AREG_PART(JR,J,9) + EGMELT(I,J)
      AIJ(I,J,IJ_MRVR)=AIJ(I,J,IJ_MRVR) +  GMELT(I,J)
      AIJ(I,J,IJ_ERVR)=AIJ(I,J,IJ_ERVR) + EGMELT(I,J)
#ifdef TRACERS_OCEAN
      TAIJN(I,J,TIJ_RVR,:)=TAIJN(I,J,TIJ_RVR,:)+ TRGMELT(:,I,J)
     *     *BYDXYP(J)
#endif
C****
      END DO
      END DO

      DO JR=1,SIZE(AREG,1)
        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,1),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_RSNOW)=AREG(JR,J_RSNOW)+AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,2),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_RUN)  =AREG(JR,J_RUN)  +AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,3),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_SNOW) =AREG(JR,J_SNOW) +AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,4),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_ACE1) =AREG(JR,J_ACE1) +AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,5),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_ACE2) =AREG(JR,J_ACE2) +AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,6),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_IMPLH)=AREG(JR,J_IMPLH)+AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,7),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_IMPLM)=AREG(JR,J_IMPLM)+AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,8),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_RVRD) = AREG(JR,J_RVRD)+AREG_SUM          

        AREG_SUM=0.
        CALL GLOBALSUM(GRID,AREG_PART(JR,:,9),AREG_SUM, ALL=.TRUE.)
        AREG(JR,J_ERVR) = AREG(JR,J_ERVR)+AREG_SUM          
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
