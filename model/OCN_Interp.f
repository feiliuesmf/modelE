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

      deallocate(aweight)

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
     *                , oCOSI=>COSIC,oSINI=>SINIC
     *                , IVSPO=>IVSP,IVNPO=>IVNP
     *                , sinpo, sinvo
#ifndef CUBE_GRID
      Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH
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
      USE MODEL_COM, only: nstep=>itime
      USE obio_com, only: pCO2
      USE TRACER_COM, only : vol2mass
      USE FLUXES, only : gtracer
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      INTEGER IER, I,J, L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S
      REAL*8, allocatable ::
     * oWEIGHT(:,:), oFOCEAN_loc(:,:)
      REAL*8 :: UNP,VNP,AWT1,AWT2
      REAL*8, ALLOCATABLE :: oG0(:,:,:), oS0(:,:,:)
     *                     , oUO1(:,:), oVO1(:,:), oTRAC(:,:,:)
     *                     , oTOT_CHLO_loc(:,:),opCO2_loc(:,:)
     *                     ,atest(:,:)
      character*80 :: name
      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP

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

c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      call halo_update(ogrid,vo(:,:,1),from=south)
      do j=oJ_0S,oJ_1S
c area weights that would have been used by HNTRP for ocean C -> ocean A
        awt1 = (sinpo(j)-sinvo(j-1))/(sinvo(j)-sinvo(j-1))
        awt2 = 1.-awt1
        i=1
          oUO1(i,j) = .5*(UO(i,j,1)+UO(oIM,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        do i=2,oIM
          oUO1(i,j) = .5*(UO(i,j,1)+UO(i-1,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        enddo
      enddo
      if(oGRID%have_south_pole) then
        oUO1(:,1) = 0.; oVO1(:,1) = 0.
      endif
      if(oGRID%have_north_pole) then ! NP U,V from prognostic polar U,V
        oUO1(:,oJM) = UO(oIM,oJM,1)*oCOSI(:) + UO(IVNPO,oJM,1)*oSINI(:)
! oVO1 currently has no effect when atm is lat-lon
        oVO1(:,oJM) = UO(IVNPO,oJM,1)*oCOSI(:) - UO(oIM,oJM,1)*oSINI(:)
      endif
      oWEIGHT(:,:) = 1.d0
      CALL INT_OG2AG(oUO1, aUO1, oWEIGHT, .FALSE., AvgPole=.FALSE.)
      CALL INT_OG2AG(oVO1, aVO1, oWEIGHT, .FALSE., AvgPole=.FALSE.)
#ifndef CUBE_GRID
      if(aGRID%have_north_pole) then ! latlon atm needs single polar vector
        UNP = SUM(aUO1(:,aJM)*aCOSI(:))*2/aIM
        VNP = SUM(aUO1(:,aJM)*aSINI(:))*2/aIM
        aUO1(1,aJM) = UNP
        aVO1(1,aJM) = VNP
      endif
#endif

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

         !opco2 is in uatm, convert to kg,CO2/kg,air
         aTRAC(:,:,NT) = aTRAC(:,:,NT) * vol2mass(nt)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air

         !gtracer is first set in TRACER_DRV, then atrac is interpolated
         !here from pco2 in the ocean and later OCNDYN sets gtracer=atrac
         !Therefore in timesetep 0 (cold start) pco2 has not yet been defined
         !and atrac has to be hard coded in here, so that we do not have
         !urealistic tracer flux at the air-sea interface during step0.

         if (nstep.eq.0) aTRAC(:,:,NT) = gtracer(NT,1,:,:)

      END DO
#endif
#endif

      DEALLOCATE(oG0, oS0, oUO1,oVO1, oTOT_CHLO_loc, opCO2_loc)
#ifdef TRACERS_OCEAN
      DEALLOCATE(oTRAC)
#endif

      deallocate(oweight,ofocean_loc)

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
     *     , aDMUA=>DMUA, aDMVA=>DMVA
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
     *     , oDMUA,oDMVA
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
      USE MODEL_COM, only: nstep=>itime
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

      deallocate(aweight)

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


      MODULE IGOG_regrid_info
!@sum IGOG_info saves instances of hntrp_type for use in
!@+   ice <-> ocean regrids
      use hntrp_mod, only : hntrp_type
      implicit none
      save

      type(hntrp_type) :: hntrp_i2o_u ! ice u C -> ocn u C
      type(hntrp_type) :: hntrp_i2o_v ! ice v C -> ocn v C
      logical :: hntrp_i2o_uv_need_init = .true.

      type(hntrp_type) :: hntrp_o2i_u ! ocn u C -> ice u A
      type(hntrp_type) :: hntrp_o2i_v ! ocn v C -> ice v A
      logical :: hntrp_o2i_uv_need_init = .true.

      END MODULE IGOG_regrid_info

      SUBROUTINE IG2OG_oceans
!@sum IG2OG_oceans interpolates relevant DYNSI outputs to the
!@+   ocean grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      use icedyn_com, only : iDMUI=>DMUI, iDMVI=>DMVI
      use ofluxes,    only : oDMUI,oDMVI
      USE ICEDYN, only : iIM=>imicdyn,iJM=>jmicdyn
      USE OCEAN,  only : oIM=>im,oJM=>jm
      USE ICEDYN, only : iDLATM=>DLATM
      USE OCEAN,  only : oDLATM=>DLATM
      USE ICEDYN,     only : iGRID=>grid_icdyn
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK
      implicit none
      real*8, dimension(:,:), allocatable ::
     &     ones_band,idmui_band,idmvi_band
      integer :: jmin,jmax

      if(iIM.eq.oIM .and. iJM.eq.oJM) then
        oDMUI(:,:) = iDMUI(:,:)
        oDMVI(:,:) = iDMVI(:,:)
        return
      endif

      if(hntrp_i2o_uv_need_init) then
        call Init_Hntrp_Type(hntrp_i2o_u,
     &     iGRID, .5d0,iDLATM,
     &     oGRID, .5d0,oDLATM,
     &     0.d0)
        call Init_Hntrp_Type(hntrp_i2o_v,
     &     iGRID, 0.d0,iDLATM,
     &     oGRID, 0.d0,oDLATM,
     &     0.d0,
     &     JMA_4interp=iJM-1,JMB_4interp=oJM-1) ! for secondary lats
        hntrp_i2o_uv_need_init = .false.
      endif

c regrid DMUI from ice C to ocn C
      jmin = hntrp_i2o_u%bpack%jband_strt
      jmax = hntrp_i2o_u%bpack%jband_stop
      ALLOCATE(iDMUI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_u%bpack, iDMUI, iDMUI_band)
      if(iGRID%have_north_pole) then
        iDMUI_band(:,iJM) = 0.
      endif
      call HNTR8_band (ones_band, iDMUI_band, hntrp_i2o_u, oDMUI)
      deallocate(iDMUI_band,ones_band)

c regrid DMVI from ice C to ocn C
      jmin = hntrp_i2o_v%bpack%jband_strt
      jmax = hntrp_i2o_v%bpack%jband_stop
      ALLOCATE(iDMVI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_v%bpack, iDMVI, iDMVI_band)
      if(iGRID%have_north_pole) then
        iDMVI_band(:,iJM) = 0.
      endif
      call HNTR8_band (ones_band, iDMVI_band, hntrp_i2o_v, oDMVI)
      deallocate(iDMVI_band,ones_band)

      if(oGRID%have_north_pole) then ! INT_AG2OG_Vector2 set them to zero
        oDMUI(:,oJM) = 0.0
        oDMVI(:,oJM) = 0.0
      endif

      return
      END SUBROUTINE IG2OG_oceans

      SUBROUTINE OG2IG_uvsurf
!@sum OG2IG_uvsurf interpolates ocean surface velocity to the
!@+   DYNSI A-grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE ICEDYN, only     : iIM=>imicdyn,iJM=>jmicdyn
      USE OCEAN,  only     : oIM=>im,oJM=>jm
      USE FLUXES, only : atm_uosurf=>uosurf,atm_vosurf=>vosurf
      USE ICEDYN_COM, only : uosurf,vosurf
      USE ICEDYN, only : iDLATM=>DLATM
      USE OCEAN,  only : oDLATM=>DLATM
      USE ICEDYN,     only : iGRID=>grid_icdyn
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK
      USE OCEAN, only : UO,VO
      implicit none
      real*8, dimension(:,:), allocatable ::
     &     ones_band,ocnu_band,ocnv_band
      integer :: jmin,jmax

c If DYNSI grid == ATM grid, simply replicate ATM copy
c (Note this requires that TOC2SST has been called first)
      if(aIM.ne.aJM .and. ! extra check if the atm is lat-lon
     &     aIM.eq.iIM .and. aJM.eq.iJM) then
        uosurf(:,:) = atm_uosurf(:,:)
        vosurf(:,:) = atm_vosurf(:,:)
        return
      endif

c call HNTRP.  this section has not been tested yet.
c      if(hntrp_o2i_uv_need_init) then
c        call Init_Hntrp_Type(hntrp_o2i_u,
c     &       oGRID, .5d0,oDLATM,
c     &       iGRID, 0.d0,iDLATM,
c     &       0.d0)
c        call Init_Hntrp_Type(hntrp_o2i_v,
c     &       oGRID, 0.d0,oDLATM,
c     &       iGRID, 0.d0,iDLATM,
c     &       0.d0,
c     &       JMA_4interp=oJM-1) ! secondary ocean lats
c        hntrp_o2i_uv_need_init = .false.
c      endif
cc regrid uosurf
c      jmin = hntrp_o2i_u%bpack%jband_strt
c      jmax = hntrp_o2i_u%bpack%jband_stop
c      ALLOCATE(ocnu_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
c      ones_band(:,:) = 1d0
c      call BAND_PACK (hntrp_o2i_u%bpack, uo(:,:,1), ocnu_band)
c      call HNTR8_band (ones_band, ocnu_band, hntrp_o2i_u, uosurf)
c      deallocate(ocnu_band,ones_band)
cc regrid vosurf
c      jmin = hntrp_o2i_v%bpack%jband_strt
c      jmax = hntrp_o2i_v%bpack%jband_stop
c      ALLOCATE(ocnv_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
c      ones_band(:,:) = 1d0
c      call BAND_PACK (hntrp_o2i_v%bpack, vo(:,:,1), ocnv_band)
c      call HNTR8_band (ones_band, ocnv_band, hntrp_o2i_v, vosurf)
c      deallocate(ocnv_band,ones_band)
cc fix up the north pole
      return
      END SUBROUTINE OG2IG_uvsurf



!------------------------------------------------------------
!!!! Unfinished code ! Unfinished code ! Unfinished code !

!#define TEST_ARRAY_BUNDLE
#ifdef TEST_ARRAY_BUNDLE

      module bundle_arrays
      implicit none
      private

      public lookup_str
      public ba_init, ba_add, ba_bundle, ba_unbundle

      interface ba_add
         module procedure ba_add_2
         module procedure ba_add_3
      end interface

      integer, parameter :: N_LOOKUP_MAX=256

      type lookup_record
        integer :: rank,km
        real*8, pointer :: src2(:,:), dest2(:,:)
        real*8, pointer :: src_w2(:,:), dest_w2(:,:)
        real*8, pointer :: src3(:,:,:), dest3(:,:,:)
        real*8, pointer :: src_w3(:,:,:), dest_w3(:,:,:)
      end type lookup_record

      type lookup_str
        integer :: si_0, si_1, sj_0, sj_1  ! source bounds
        integer :: di_0, di_1, dj_0, dj_1  ! destination bounds
        integer  :: n_lookup=0
        type (lookup_record) :: lr(N_LOOKUP_MAX)
      end type lookup_str


      contains

      subroutine ba_init( lstr, si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 )
      type (lookup_str) :: lstr
      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      
      lstr%si_0 = si_0
      lstr%si_1 = si_1
      lstr%sj_0 = sj_0
      lstr%sj_1 = sj_1
      lstr%di_0 = di_0
      lstr%di_1 = di_1
      lstr%dj_0 = dj_0
      lstr%dj_1 = dj_1

      lstr%n_lookup = 0
      end subroutine ba_init


      subroutine ba_add_2( lstr, src, dest, src_w, dest_w )
      type (lookup_str) :: lstr
      real*8, dimension(:,:), target :: src, dest
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      lstr%n_lookup = lstr%n_lookup + 1
      if ( lstr%n_lookup > N_LOOKUP_MAX )
     &     call stop_model("ba_add: increase N_LOOKUP_MAX", 255)

      lstr%lr(lstr%n_lookup)%src2 => src(:,:)
      lstr%lr(lstr%n_lookup)%dest2 => dest(:,:)
      lstr%lr(lstr%n_lookup)%rank = 2
      lstr%lr(lstr%n_lookup)%km = 1

      if ( present(src_w) .and. present(dest_w) ) then
        lstr%lr(lstr%n_lookup)%src_w2 => src_w(:,:)
        lstr%lr(lstr%n_lookup)%dest_w2 => dest_w(:,:)
      else
        nullify( lstr%lr(lstr%n_lookup)%src_w2 )
        nullify( lstr%lr(lstr%n_lookup)%dest_w2 )
      endif

      end subroutine ba_add_2


      subroutine ba_add_3( lstr, src, dest, src_w, dest_w )
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), target :: src, dest
      real*8, dimension(:,:,:), target, optional :: src_w, dest_w

      lstr%n_lookup = lstr%n_lookup + 1
      if ( lstr%n_lookup > N_LOOKUP_MAX )
     &     call stop_model("ba_add: increase N_LOOKUP_MAX", 255)

      lstr%lr(lstr%n_lookup)%src3 => src(:,:,:)
      lstr%lr(lstr%n_lookup)%dest3 => dest(:,:,:)
      lstr%lr(lstr%n_lookup)%rank = 3
      lstr%lr(lstr%n_lookup)%km = size(src,1)

      if ( present(src_w) .and. present(dest_w) ) then
        lstr%lr(lstr%n_lookup)%src_w3 => src_w(:,:,:)
        lstr%lr(lstr%n_lookup)%dest_w3 => dest_w(:,:,:)
      else
        nullify( lstr%lr(lstr%n_lookup)%src_w3 )
        nullify( lstr%lr(lstr%n_lookup)%dest_w3 )
      endif

      end subroutine ba_add_3


      subroutine ba_bundle( lstr, buf_s, buf_d )
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer k, n

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      n = 1
      do k=1,lstr%n_lookup
        n = n+lstr%lr(k)%km
      enddo

      allocate( buf_s(n, si_0:si_1, sj_0:sj_1) )
      allocate( buf_d(n, di_0:di_1, dj_0:dj_1) )

      n = 1
      do k=1,lstr%n_lookup
        select case(lstr%lr(k)%rank)
        case(2)
          buf_s(n,:,:) = lstr%lr(k)%src2(:,:)
          if ( associated(lstr%lr(k)%src_w2 ) )
     &         buf_s(n,:,:) = buf_s(n,:,:) * lstr%lr(k)%src_w2(:,:)
        case(3)
          buf_s(n:n+lstr%lr(k)%km-1,:,:) = lstr%lr(k)%src3(:,:,:)
          if ( associated(lstr%lr(k)%src_w3 ) )
     &         buf_s(n:n+lstr%lr(k)%km-1,:,:) =
     &         buf_s(n:n+lstr%lr(k)%km-1,:,:) * lstr%lr(k)%src_w3(:,:,:)
        end select
        n = n+lstr%lr(k)%km
      enddo

      end subroutine ba_bundle


      subroutine ba_unbundle( lstr, buf_s, buf_d )
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer k, n

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      n = 1
      do k=1,lstr%n_lookup
        select case(lstr%lr(k)%rank)
        case(2)
          lstr%lr(k)%dest2(:,:) = buf_d(n,:,:)
          if ( associated(lstr%lr(k)%dest_w2 ) )
     &         lstr%lr(k)%dest2(:,:) =
     &         lstr%lr(k)%dest2(:,:) / lstr%lr(k)%dest_w2(:,:)
        case(3)
          lstr%lr(k)%dest3(:,:,:) = buf_d(n:n+lstr%lr(k)%km-1,:,:)
          if ( associated(lstr%lr(k)%dest_w3 ) )
     &         lstr%lr(k)%dest3(:,:,:) =
     &         lstr%lr(k)%dest3(:,:,:) * lstr%lr(k)%dest_w3(:,:,:)
        end select
        n = n+lstr%lr(k)%km
      enddo

      deallocate( buf_s )
      deallocate( buf_d )     

      end subroutine ba_unbundle



      end module bundle_arrays


      subroutine do_interpalation(lstr)
      USE bundle_arrays

      implicit none
      type (lookup_str) :: lstr
      real*8, pointer :: buf_s(:,:,:), buf_d(:,:,:)

      call ba_bundle( lstr, buf_s, buf_d )

      !! do interpotation here from bundle buf_s(:,:,:)
      !! to bundle buf_d(:,:,:)

      call ba_unbundle( lstr, buf_s, buf_d )

      end subroutine do_interpalation


      SUBROUTINE AG2OG_precip_example

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

      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aPREC=>PREC, aEPREC=>EPREC
     *     , aRUNPSI=>RUNPSI, aSRUNPSI=>SRUNPSI, aERUNPSI=>ERUNPSI

      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI

      USE INT_AG2OG_MOD, only : INT_AG2OG
      USE OCEANR_DIM, only : ogrid
      USE bundle_arrays

      IMPLICIT NONE

      INTEGER N

      REAL*8, allocatable :: aWEIGHT(:,:), oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:), oWEIGHT1(:,:)

      character*80 :: name

      type (lookup_str) :: lstr

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- need to allocate separate array for different weights
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))


      !-- initializing "lstr"
      call ba_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aWEIGHT(:,:) = 1.- aRSI(:,:)  !!  open ocean fraction

      !!CALL INT_AG2OG(aEPREC,oEPREC, aWEIGHT)
      !-- example of edding an array with weights
      call ba_add( lstr, aEPREC, oEPREC, aWEIGHT, oWEIGHT )
      !-- example of edding an array without weights
      call ba_add( lstr, aWEIGHT, oWEIGHT )

      oEPREC(:,:) = oEPREC(:,:)*OXYP(:,:)

      !-- use aWEIGHT1 sine aWEIGHT is already busy
      aWEIGHT1(:,:) = aRSI(:,:)
      !!CALL INT_AG2OG(aRUNPSI,oRUNPSI, aWEIGHT)
      CALL ba_add( lstr,aRUNPSI,oRUNPSI, aWEIGHT1, oWEIGHT1)
      CALL ba_add( lstr,aWEIGHT1, oWEIGHT1)

      !CALL INT_AG2OG(aSRUNPSI,oSRUNPSI, aWEIGHT)
      CALL ba_add( lstr,aSRUNPSI,oSRUNPSI, aWEIGHT1, oWEIGHT1)
      oSRUNPSI(:,:) = oSRUNPSI(:,:)*OXYP(:,:)

      !CALL INT_AG2OG(aERUNPSI,oERUNPSI, aWEIGHT)
      CALL ba_add( lstr,aERUNPSI,oERUNPSI, aWEIGHT1, oWEIGHT1)
      oERUNPSI(:,:) = oERUNPSI(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = 1.d0
      !CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)
      !-- since weights==1 no need to register them
      CALL ba_add( lstr, aRSI, oRSI)

      call do_interpalation(lstr)

      !-- deallocate only after interpolation is done
      deallocate(aweight, oweight)
      deallocate(aweight1, oweight1)

      RETURN
      END SUBROUTINE AG2OG_precip_example


#endif

