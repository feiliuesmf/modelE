#include "rundeck_opts.h"

      SUBROUTINE AG2OG_precip

!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid

      USE OCEAN, only : oDXYPO=>DXYPO,oDLATM=>DLATM, OXYP
      Use GEOM,  only : AXYP

      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : NTM
#ifdef TRACERS_WATER
      USE FLUXES,  only : aTRPREC=>TRPREC, aTRUNPSI=>TRUNPSI
      USE OFLUXES, only : oTRPREC, oTRUNPSI
#endif
#endif
      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aPREC=>PREC, aEPREC=>EPREC
     *     , aRUNPSI=>RUNPSI, aSRUNPSI=>SRUNPSI, aERUNPSI=>ERUNPSI

      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI

      USE INT_AG2OG_MOD, only : INT_AG2OG

      IMPLICIT NONE

      INTEGER N

      REAL*8, allocatable :: aWEIGHT(:,:)

      character*80 :: name

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      aWEIGHT(:,:) = 1.- aRSI(:,:)  !!  open ocean fraction
      CALL INT_AG2OG(aPREC,oPREC, aWEIGHT)

      CALL INT_AG2OG(aEPREC,oEPREC, aWEIGHT)
      oEPREC(:,:) = oEPREC(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aRUNPSI,oRUNPSI, aWEIGHT)

      CALL INT_AG2OG(aSRUNPSI,oSRUNPSI, aWEIGHT)
      oSRUNPSI(:,:) = oSRUNPSI(:,:)*OXYP(:,:)

      CALL INT_AG2OG(aERUNPSI,oERUNPSI, aWEIGHT)
      oERUNPSI(:,:) = oERUNPSI(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)

      aWEIGHT(:,:) = 1.- aRSI(:,:)
      DO N=1,NTM
        aTRPREC(N,:,:) = aTRPREC(N,:,:)/AXYP(:,:)
      END DO
      CALL INT_AG2OG(aTRPREC,oTRPREC, aWEIGHT, NTM)
      DO N=1,NTM
        oTRPREC(N,:,:) = oTRPREC(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNPSI,oTRUNPSI, aWEIGHT, NTM)
      DO N=1,NTM
        oTRUNPSI(N,:,:) = oTRUNPSI(N,:,:)*OXYP(:,:)
      END DO
#endif

      RETURN
      END SUBROUTINE AG2OG_precip


      SUBROUTINE OG2AG_TOC2SST
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oLM=>LMO, oFOCEAN=>FOCEAN
     *                , oDXYPO=>DXYPO, OXYP, oIMAXJ=>IMAXJ
     *                , oCOSU=>COSU,oSINU=>SINU
     *                , IVSPO=>IVSP,IVNPO=>IVNP

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : MO, UO,VO, G0M
     *     , S0M, OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , TRMO
#endif
      USE AFLUXES, only : aMO, aUO1,aVO1, aG0
     *     , aS0, aOGEOZ,aOGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , aTRAC
#endif
      USE OFLUXES, only : oRSI

#ifdef TRACERS_OCEAN
      USE TRACER_COM, only: NTM
      USE OCN_TRACER_COM, only: conc_from_fw
#endif
#ifdef TRACERS_OceanBiology
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field
      USE obio_com, only: tot_chlo
      USE FLUXES, only : CHL
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE obio_com, only: pCO2
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      INTEGER IER, I,J, L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1
      REAL*8, allocatable ::
     * oWEIGHT(:,:), oFOCEAN_loc(:,:)

      REAL*8, ALLOCATABLE :: oG0(:,:,:), oS0(:,:,:)
     *                     , oUO1(:,:), oVO1(:,:), oTRAC(:,:,:)
     *                     , oTOT_CHLO_loc(:,:),opCO2_loc(:,:)
     *                     ,atest(:,:)
      character*80 :: name
      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP

      ALLOCATE
     *  (oG0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oS0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
      ALLOCATE
     *  (oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE
     *  (oTRAC(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,NTM), STAT = IER)
#endif
      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)
      ALLOCATE
     *  (opCO2_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) ,STAT = IER)

      allocate(oFOCEAN_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &         oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

      CALL OCN_UNPACK (oGRID,oFOCEAN,oFOCEAN_loc)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(MO,aMO,oWEIGHT,oLM,2,.FALSE.)

      oG0(:,:,:) = 0.d0
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oG0(I,J,L) = G0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      CALL INT_OG2AG(oG0,aG0, oWEIGHT, 2,2,.TRUE.)

      oS0(:,:,:) = 0.d0
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oS0(I,J,L) = S0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      CALL INT_OG2AG(oS0,aS0, oWEIGHT, 2,2,.TRUE.)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(OGEOZ,aOGEOZ, oWEIGHT, .FALSE.)

      CALL INT_OG2AG(OGEOZ_SV,aOGEOZ_SV, oWEIGHT, .FALSE.)

      oUO1(:,:) = UO(:,:,1)
      oVO1(:,:) = VO(:,:,1)
      oWEIGHT(:,:) = 1.d0
      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
     *             , IVSPO,IVNPO)

#ifdef TRACERS_OCEAN
C**** surface tracer concentration
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO NT = 1,NTM
        if (conc_from_fw(nt)) then  ! define conc from fresh water
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J)
     *             -S0M(I,J,1))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        else  ! define conc from total sea water mass
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        end if
      END DO
      CALL INT_OG2AG(oTRAC,aTRAC, oWEIGHT, NTM,NTM,.TRUE.)

#ifdef TRACERS_OceanBiology
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation
      oWEIGHT(:,:) = oFOCEAN_loc(:,:) * (1.d0 - oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            oTOT_CHLO_loc(I,J) = tot_chlo(I,J)
          ELSE
            oTOT_CHLO_loc(I,J)=0.
          END IF
        END DO
      END DO
      CALL INT_OG2AG(oTOT_CHLO_loc,CHL, oWEIGHT, .FALSE.)
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      oWEIGHT(:,:) = oFOCEAN_loc(:,:)*(1.d0-oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            opCO2_loc(I,J) = pCO2(I,J)
          ELSE
            opCO2_loc(I,J)=0.
          END IF
        END DO
      END DO
      DO NT = 1,NTM
      CALL INT_OG2AG(opCO2_loc,aTRAC(:,:,NT), oWEIGHT, .FALSE.)
      END DO
#endif
#endif

      DEALLOCATE(oG0, oS0, oUO1,oVO1, oTOT_CHLO_loc, opCO2_loc)
#ifdef TRACERS_OCEAN
      DEALLOCATE(oTRAC)
#endif

      RETURN
      END SUBROUTINE OG2AG_TOC2SST

      SUBROUTINE AG2OG_oceans
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oFOCEAN=>FOCEAN
     *                , oDXYPO=>DXYPO, OXYP
     *                , IVSPO=>IVSP,IVNPO=>IVNP

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#if defined (TRACERS_OceanBiology) && !defined (TRACERS_GASEXCH_ocean)
      USE OCN_TRACER_COM, only: NTM
#else
      USE TRACER_COM, only: NTM
#endif
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only : pack_data
      USE OCEANR_DIM, only : ogrid

      USE SEAICE_COM, only : aRSI=>RSI

      USE FLUXES, only : aSOLAR=>SOLAR, aE0=>E0, aEVAPOR=>EVAPOR
     *     , aRUNOSI=>RUNOSI,aERUNOSI=>ERUNOSI,aSRUNOSI=>SRUNOSI
     *     , aFLOWO=>FLOWO,aEFLOWO=>EFLOWO, aAPRESS=>APRESS
     *     , aMELTI=>MELTI,aEMELTI=>EMELTI,aSMELTI=>SMELTI
     *     , aDMUA=>DMUA, aDMVA=>DMVA, aDMUI=>DMUI, aDMVI=>DMVI
     *     , aGMELT=>GMELT, aEGMELT=>EGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO=>TRFLOWO, aTREVAPOR=>TREVAPOR
     *     , aTRUNOSI=>TRUNOSI, aTRMELTI=>TRMELTI
     *     , aTRGMELT=>TRGMELT
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP=>TRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE FLUXES, only : aTRGASEX=>TRGASEX
      USE TRACER_GASEXCH_COM, only: tracflx_glob
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only : avisdir    => FSRDIR
     *                   ,asrvissurf => SRVISSURF
     *                   ,avisdif    => FSRDIF
     *                   ,anirdir    => DIRNIR
     *                   ,anirdif    => DIFNIR
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : aCOSZ1=>COSZ1
      USE PBLCOM, only : awind=>wsavg
#endif

      USE OFLUXES, only : oRSI, oSOLAR, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMUA,oDMVA, oDMUI,oDMVI
     *     , oGMELT, oEGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
     *     , oTRGMELT
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE OFLUXES, only : oTRGASEX
#endif
#ifdef OBIO_RAD_coupling
      USE obio_forc, only : ovisdir,ovisdif,onirdir,onirdif
#endif
#ifdef TRACERS_OceanBiology
      USE obio_forc, only : osolz
      USE obio_forc, only : owind
#endif

      Use GEOM,  only : AXYP,aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN
      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob

      USE INT_AG2OG_MOD, only : INT_AG2OG

      IMPLICIT NONE

      INTEGER, PARAMETER :: NSTYPE=4
      INTEGER I,J, N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER oJ_0,oJ_1

      REAL*8,
     * DIMENSION(aGRID%I_STRT:aGRID%I_STOP,aGRID%J_STRT:aGRID%J_STOP)::
     * aFact

      REAL*8, allocatable :: aweight(:,:)
      real*8, dimension(oIM,oJM) :: otest_glob
      real*4, dimension(oIM,oJM) :: tr4

      character*80 :: title,name

      allocate (aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP

      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)
      
      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
            aFLOWO(I,J) = aFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aFLOWO,oFLOWO, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aEFLOWO(I,J) = aEFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aEFLOWO,oEFLOWO, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aMELTI(I,J) = aMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aMELTI,oMELTI, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aEMELTI(I,J) = aEMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aEMELTI,oEMELTI, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aSMELTI(I,J) = aSMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aSMELTI,oSMELTI, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aGMELT(I,J) = aGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aGMELT,oGMELT, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/AXYP(I,J)
            aEGMELT(I,J) = aEGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aEGMELT,oEGMELT, aWEIGHT)
      oEGMELT(:,:) = oEGMELT(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aRUNOSI,oRUNOSI, aWEIGHT)

      CALL INT_AG2OG(aERUNOSI,oERUNOSI, aWEIGHT)

      CALL INT_AG2OG(aSRUNOSI,oSRUNOSI, aWEIGHT)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aAPRESS,oAPRESS, aWEIGHT)

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aE0,oE0, aWEIGHT, NSTYPE,1)

      CALL INT_AG2OG(aEVAPOR,oEVAPOR, aWEIGHT, NSTYPE,1)

      CALL INT_AG2OG(aSOLAR,oSOLAR, aWEIGHT, 3,3,1)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aSOLAR,oSOLAR, aWEIGHT, 3,3,3)

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              aTRFLOWO(N,I,J) = aTRFLOWO(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRFLOWO,oTRFLOWO, aWEIGHT, NTM)

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aTRMELTI(N,I,J) = aTRMELTI(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRMELTI,oTRMELTI, aWEIGHT, NTM)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNOSI,oTRUNOSI, aWEIGHT, NTM)

      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/AXYP(I,J)
              aTRGMELT(N,I,J) = aTRGMELT(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRGMELT,oTRGMELT, aWEIGHT, NTM)
      DO N=1,NTM
        oTRGMELT(N,:,:) = oTRGMELT(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aTREVAPOR,oTREVAPOR, aWEIGHT, NTM, NSTYPE,1)

#ifdef TRACERS_DRYDEP
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRDRYDEP,oTRDRYDEP, aWEIGHT, NTM, NSTYPE,1)
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRGASEX,oTRGASEX, tracflx_glob, aWEIGHT
     *             , NTM, NSTYPE,1)
#endif

#ifdef OBIO_RAD_coupling
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aVISDIR,aSRVISSURF,oVISDIR, aWEIGHT)

      CALL INT_AG2OG(aVISDIF,oVISDIF, aWEIGHT)

      CALL INT_AG2OG(aNIRDIR,oNIRDIR, aWEIGHT)

      CALL INT_AG2OG(aNIRDIF,oNIRDIF, aWEIGHT)
#endif

#ifdef TRACERS_OceanBiology
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aCOSZ1,oSOLZ, aWEIGHT)

      CALL INT_AG2OG(aWIND,oWIND, aWEIGHT)
#endif
      aWEIGHT(:,:) = 1.d0
      
      CALL INT_AG2OG(aDMUA,aDMVA,oDMUA,oDMVA, aWEIGHT,aFOCEAN_loc
     *              ,NSTYPE,1)

      CALL INT_AG2OG(aDMUI,aDMVI,oDMUI,oDMVI, aWEIGHT, IVSPO,IVNPO)

      END SUBROUTINE AG2OG_oceans

c      subroutine debug_acoupling(name,ar8)
c!@auth D. Gueyffier
c      USE RESOLUTION, only : aIM=>IM, aJM=>JM
c      use domain_decomp_atm, only : agrid=>grid,pack_data,am_i_root
c      real*8, intent(in) :: ar8(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
c     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
c      real*8, allocatable :: ar8glob(:,:,:)
c      real*4, allocatable :: ar4(:,:,:)
c      character*80 :: name,title,fname
c
c      allocate(ar8glob(aIM,aJM,6),ar4(aIM,aJM,6))
c
c      call pack_data(agrid,ar8,ar8glob)
c
c      if (am_i_root()) then
c         ar4=ar8glob
c         title="testa-"//trim(name)
c         fname=title
c         open(20,FILE=fname,FORM='unformatted', STATUS='unknown')
c         write(20) title,ar4
c         write(*,*) "max "//trim(name),"=",maxval(ar8glob),
c     &        "min "//trim(name),"=",minval(ar8glob),
c     &        "maxloc "//trim(name),"=",maxloc(ar8glob),
c     &        "minloc "//trim(name),"=",minloc(ar8glob)
c         close(20)
c         name=""
c      endif
c
c      deallocate(ar8glob,ar4)
c
c      end subroutine debug_acoupling
c

c      subroutine debug_ocoupling(name,or8)
c!@auth D. Gueyffier
c      USE OCEAN, only : oIM=>IM, oJM=>JM, oFOCEAN=>FOCEAN
c      use domain_decomp_1d, only : pack_data,am_i_root
c      USE OCEANR_DIM, only : ogrid
c      real*8, intent(in) :: or8(oIM,
c     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
c      real*8, allocatable :: or8glob(:,:)
c      real*4, allocatable :: or4(:,:)
c      character*80 :: title,name,fname
c
c      allocate(or8glob(oIM,oJM),or4(oIM,oJM))
c
cc      write(*,*) "max o dist=",maxval(
cc     &     or8(:,oGRID%J_STRT:oGRID%J_STOP)
cc     &                                ),
cc     &           "min odist=",minval(
cc     &     or8(:,oGRID%J_STRT:oGRID%J_STOP)
cc     &                               )
c
c      call pack_data(ogrid,or8,or8glob)
c
c      if (am_i_root()) then
c         or4=or8glob
c         title="testo-"//trim(name)
c         fname=title
c         open(30,FILE=fname,FORM='unformatted', STATUS='unknown')
c         write(30) title,or4
c         write(*,*) "max "//trim(name),"=",maxval(or8glob),
c     &        "min "//trim(name),"=",minval(or8glob),
c     &        "maxloc "//trim(name),"=",maxloc(or8glob),
c     &        "minloc "//trim(name),"=",minloc(or8glob)
c 
c         close(30)
c         name=""
c      endif
c
c      deallocate(or8glob,or4)
c
c      end subroutine debug_ocoupling
c
      
      SUBROUTINE OG2AG_oceans
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM,oJM=>JM, oFOCEAN=>FOCEAN

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : UNPACK_DATA
      USE OCEANR_DIM, only : ogrid

      USE MODEL_COM, ONLY : aFOCEAN_loc=>FOCEAN
      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob

      USE FLUXES, only : aDMSI=>DMSI,aDHSI=>DHSI,aDSSI=>DSSI
#ifdef TRACERS_ON
     *     , aDTRSI=>DTRSI
#endif

      USE OFLUXES, only : oDMSI, oDHSI, oDSSI
#ifdef TRACERS_OCEAN
     *     , oDTRSI
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      REAL*8,
     * DIMENSION(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)::
     * oWEIGHT, oFOCEAN_loc

      CALL UNPACK_DATA (oGRID,oFOCEAN,oFOCEAN_loc)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(oDMSI,aDMSI, oWEIGHT, 2)

      CALL INT_OG2AG(oDHSI,aDHSI, oWEIGHT, 2)

      CALL INT_OG2AG(oDSSI,aDSSI, oWEIGHT, 2)

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (NTM == NTM_ATM) THEN

        oWEIGHT(:,:) = oFOCEAN_loc(:,:)
        CALL INT_OG2AG(oDTRSI,aDTRSI, oWEIGHT, 2, NTM)

      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif
      RETURN
      END SUBROUTINE OG2AG_oceans


