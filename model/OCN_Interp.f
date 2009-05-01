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

      REAL*8, 
     * DIMENSION(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)::
     * aWEIGHT

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


      SUBROUTINE OG2AG
!@sum  OG2AG gathers all necessary arrays on the ocean grid, interpolates
!!      them to the atmospheric grid, scatters them on the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT

      IMPLICIT NONE

      call gather_ocean1 ! mo,uo,vo,g0m,s0m,ogz's,tr
!
!!!  Do interpolation to the atmospheric grid
!
#ifndef CUBE_GRID
      if(AM_I_ROOT()) then
        call INT_OG2AG
      end if
#endif

      call scatter_ocean1

      RETURN
      END SUBROUTINE OG2AG

      SUBROUTINE INT_OG2AG
!@sum  INT_OG2AG is for interpolation of arrays from ocean grid
!!      to the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0
      USE MODEL_COM, only :
#if defined(TRACERS_GASEXCH_ocean) || defined(TRACERS_OceanBiology)
      USE MODEL_COM, only: nstep=>itime
#endif

      USE RESOLUTION, only : IMA=>IM, JMA=>JM
      USE OCEANRES,   only : IMO, JMO, LMO

#ifdef TRACERS_OCEAN
      USE TRACER_COM, only: ntm
      USE OCN_TRACER_COM, only: conc_from_fw
#endif

      USE OCEAN, only : MO_glob, UO_glob,VO_glob, G0M_glob
     *     ,S0M_glob, OGEOZ_glob, OGEOZ_SV_glob, oFOCEAN=>FOCEAN
#ifdef TRACERS_OCEAN
     *     , TRMO_glob
#endif

      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob
      USE AFLUXES, only : aMO_glob, aUO1_glob,aVO1_glob, aG0_glob
     *     , aS0_glob, aOGEOZ_glob, aOGEOZ_SV_glob
#ifdef TRACERS_OCEAN
     *     , aTRAC_glob
#endif
#ifdef TRACERS_OceanBiology
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field
     *     , CHL_glob
      USE obio_com, only: tot_chlo_glob
      USE FLUXES, only : CHL
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE obio_com, only: pCO2_glob
#endif

      USE OCEAN, only : oDXYPO=>DXYPO, oIMAXJ=>IMAXJ,oDLATM=>DLATM
     *     , IVSPO=>IVSP, IVNPO=>IVNP, oCOSU=>COSU,oSINU=>SINU
#ifndef CUBE_GRID
      Use GEOM,  only : aIMAXJ=>IMAXJ,aDLATM=>DLATM
     *     , aCOSI=>COSIP,aSINI=>SINIP
#endif
      IMPLICIT NONE

      integer I,J,L, NT ,im1

      REAL*8, DIMENSION(IMA,JMA) :: aFtemp
      REAL*8, DIMENSION(IMO,JMO) :: oFtemp, oONES, oFweight
      REAL*8  oVOsp, oVOnp
      REAL*8  aUO1sp, aVO1sp, aUO1np, aVO1np
      REAL*8  SUM_oG0M, SUM_oFtemp, SUM_aG0, diff


      oONES(:,:) = 1.d0
#ifndef CUBE_GRID
      call HNTR80 (IMO,JMO,0.d0,oDLATM, IMA,JMA,0.d0,aDLATM, 0.d0)

!!!  Ocean mass for the 1st two layers

      DO L = 1,2
        oFtemp(:,:) = MO_glob(:,:,L)
        oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
        call HNTR8P (oFOCEAN, oFtemp, aFtemp)
        aMO_glob(:,:,L) = aFtemp(:,:)
      END DO

!!!  Enthalpy for the 1st two ocean layers

      oFweight(:,:) = MO_glob(:,:,1)*oFOCEAN(:,:)
      oFweight(2:IMO,JMO) = oFweight(1,JMO)
      DO L = 1,2
        SUM_oG0M = 0.
        SUM_oFtemp = 0.
        DO J=1,JMO
          DO I=1,oIMAXJ(J)
            IF (oFOCEAN(I,J).gt.0.) THEN
              oFtemp(I,J) = G0M_glob(I,J,L)/(MO_glob(I,J,L)*oDXYPO(J))

              SUM_oG0M = SUM_oG0M + G0M_glob(I,J,L)
              SUM_oFtemp = SUM_oFtemp
     *          + (oFtemp(I,J)*MO_glob(I,J,L)*oDXYPO(J)*oFOCEAN(I,J))

            END IF
          END DO
        END DO

        oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
        call HNTR8P (oFweight, oFtemp, aFtemp)
        aG0_glob(:,:,L) = aFtemp(:,:)

      END DO

!!!  Salinity for the 1st two ocean layers

      DO L = 1,2
        DO J=1,JMO
          DO I=1,oIMAXJ(J)
            IF (oFOCEAN(I,J).gt.0.) THEN
              oFtemp(I,J) = S0M_glob(I,J,L)/(MO_glob(I,J,L)*oDXYPO(J))
            END IF
          END DO
        END DO

        oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
        call HNTR8P (oFweight, oFtemp, aFtemp)
        aS0_glob(:,:,L) = aFtemp(:,:)
      END DO

      OGEOZ_glob(2:IMO,JMO) = OGEOZ_glob(1,JMO)
      call HNTR8P (oFOCEAN, OGEOZ_glob, aOGEOZ_glob)
      OGEOZ_SV_glob(2:IMO,JMO) = OGEOZ_SV_glob(1,JMO)
      call HNTR8P (oFOCEAN, OGEOZ_SV_glob, aOGEOZ_SV_glob)

#ifdef TRACERS_OCEAN
C**** surface tracer concentration
      DO NT = 1,NTM
        if (conc_from_fw(nt)) then  ! define conc from fresh water
        DO J=1,JMO
          DO I=1,oIMAXJ(J)
            IF (oFOCEAN(I,J).gt.0.) THEN
              oFtemp(I,J)=TRMO_glob(I,J,1,NT)/(MO_glob(I,J,1)*oDXYPO(J)
     *             -S0M_glob(I,J,1))
            ELSE
              oFtemp(I,J)=0.
            END IF
          END DO
        END DO
        else  ! define conc from total sea water mass
        DO J=1,JMO
          DO I=1,oIMAXJ(J)
            IF (oFOCEAN(I,J).gt.0.) THEN
              oFtemp(I,J)=TRMO_glob(I,J,1,NT)/(MO_glob(I,J,1)*oDXYPO(J))
            ELSE
              oFtemp(I,J)=0.
            END IF
          END DO
        END DO
        end if
        oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
        call HNTR8P (oFweight, oFtemp, aFtemp)
        aTRAC_glob(:,:,NT)=aFtemp(:,:)
      END DO

#ifdef TRACERS_OceanBiology   
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field 
      DO J=1,JMO
        DO I=1,oIMAXJ(J)
          IF (oFOCEAN(I,J).gt.0.) THEN
            oFtemp(I,J) = tot_chlo_glob(I,J)
!           write(*,'(a,3i5,e12.4)')'OCN_Interp:',
!    .       nstep,i,j,tot_chlo_glob(i,j)
            ELSE
              oFtemp(I,J)=0.
          END IF
        END DO
      END DO
      oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
      call HNTR8P (oFweight, oFtemp, aFtemp)
      CHL_glob(:,:) = aFtemp(:,:)
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
!in the CO2 gas exchange experiments we do not use trmo_glob but rather pco2
!pco2 is in uatm
      DO NT = 1,NTM
        DO J=1,JMO
          DO I=1,oIMAXJ(J)
            IF (oFOCEAN(I,J).gt.0.) THEN
              !!oFtemp(I,J)=pCO2_glob(i,j)/(MO_glob(I,J,1)*oDXYPO(J))
              !atrac is defined here in uatm (ppmv) 
              !for consistensy with the units of trs=tr=atmco2 in PBL.f
              oFtemp(I,J)=pCO2_glob(i,j)
            ELSE
              oFtemp(I,J)=0.
            END IF
          END DO
        END DO
        oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
        call HNTR8P (oFweight, oFtemp, aFtemp)
        aTRAC_glob(:,:,NT)=aFtemp(:,:)
      END DO
#endif
#endif

!!!  U velocity for the 1st ocean layer.

      oVOsp = UO_glob(IVSPO,  1,1)
      oVOnp = UO_glob(IVNPO,JMO,1)

      UO_glob(:,  1,1) = UO_glob(IMO,  1,1)*oCOSU(:) - oVOsp*oSINU(:)
      UO_glob(:,JMO,1) = UO_glob(IMO,JMO,1)*oCOSU(:) + oVOnp*oSINU(:)

      call HNTR80 (IMO,JMO,0.5d0,oDLATM, IMA,JMA,0.d0,aDLATM, 0.d0)
      call HNTR8  (oONES, UO_glob(1,1,1), aUO1_glob)  !!  U-grid => A-grid

      aUO1sp = SUM(aUO1_glob(:,  1)*aCOSI(:))*2/IMA
      aVO1sp = SUM(aUO1_glob(:,  1)*aSINI(:))*2/IMA
      aUO1np = SUM(aUO1_glob(:,JMA)*aCOSI(:))*2/IMA
      aVO1np = SUM(aUO1_glob(:,JMA)*aSINI(:))*2/IMA

      aUO1_glob(1,  1) = aUO1sp
      aUO1_glob(1,JMA) = aUO1np

!!!  V velocity for the 1st ocean layer

      call HNTR80 (IMO,JMO-1,0.d0,oDLATM, IMA,JMA,0.d0,aDLATM, 0.d0)
      call HNTR8  (oONES, VO_glob(1,1,1), aVO1_glob)  !!  V-grid => A-grid

      aVO1_glob(1,  1) = aVO1sp
      aVO1_glob(1,JMA) = aVO1np
#endif
      RETURN
      END SUBROUTINE INT_OG2AG



      SUBROUTINE gather_ocean1

!@sum  gather_ocean1  gathers necessary arrays on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      use domain_decomp_atm, only: agrid=>grid, atm_pack=>pack_data
      use domain_decomp_1d, only: ocn_pack=>pack_data
      use OCEANR_DIM, only : ogrid

      USE MODEL_COM, ONLY : aFOCEAN_loc=>FOCEAN
      USE OCEAN, only : MO, UO,VO, G0M
     *     , S0M, OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , TRMO
#endif
      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob
      USE OCEAN, only : MO_glob, UO_glob,VO_glob, G0M_glob
     *     ,S0M_glob, OGEOZ_glob, OGEOZ_SV_glob
#ifdef TRACERS_OCEAN
     *     , TRMO_glob
#endif

      CALL ATM_PACK(agrid, aFOCEAN_loc, aFOCEAN)

      CALL OCN_PACK(ogrid,   MO   ,    MO_glob)
      CALL OCN_PACK(ogrid,   UO   ,    UO_glob)
      CALL OCN_PACK(ogrid,   VO   ,    VO_glob)

      CALL OCN_PACK(ogrid,OGEOZ   , OGEOZ_glob)
      CALL OCN_PACK(ogrid,OGEOZ_SV,OGEOZ_SV_glob)

      CALL OCN_PACK(ogrid,  G0M   ,   G0M_glob)
      CALL OCN_PACK(ogrid,  S0M   ,   S0M_glob)
#ifdef TRACERS_OCEAN
      CALL OCN_PACK(ogrid,  TRMO  ,  TRMO_glob)
#endif

      RETURN
      END SUBROUTINE gather_ocean1

      SUBROUTINE scatter_ocean1

!@sum  scatter_ocean1  scatters necessary arrays on the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0

      use domain_decomp_atm, only: agrid=>grid,unpack_data

      USE AFLUXES, only : aMO, aUO1,aVO1, aG0
     *     , aS0, aOGEOZ,aOGEOZ_SV
     *     , aMO_glob, aUO1_glob,aVO1_glob, aG0_glob
     *     , aS0_glob, aOGEOZ_glob, aOGEOZ_SV_glob
#ifdef TRACERS_ON
     *     , aTRAC, aTRAC_glob
#endif
#ifdef TRACERS_OceanBiology   
!do not do this in seawifs case
     *     , CHL_glob
      USE FLUXES, only : CHL 
#endif

      CALL UNPACK_DATA(agrid,       aMO_glob, aMO  )
      CALL UNPACK_DATA(agrid,      aUO1_glob, aUO1 )
      CALL UNPACK_DATA(agrid,      aVO1_glob, aVO1 )

      CALL UNPACK_DATA(agrid,    aOGEOZ_glob, aOGEOZ )
      CALL UNPACK_DATA(agrid, aOGEOZ_SV_glob, aOGEOZ_SV )

      CALL UNPACK_DATA(agrid,       aG0_glob, aG0 )
      CALL UNPACK_DATA(agrid,       aS0_glob, aS0 )
#ifdef TRACERS_ON
      CALL UNPACK_DATA(agrid,     aTRAC_glob, aTRAC )
#endif
#ifdef TRACERS_OceanBiology   
!do not do this in seawifs case
!when chl is computed in the ocean biology module
!it needs to be passed back into the atmosphere
      CALL UNPACK_DATA(agrid,     CHL_glob, CHL )
#endif

      RETURN
      END SUBROUTINE scatter_ocean1

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

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid, PACK_DATA
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

      REAL*8, 
     * DIMENSION(aGRID%I_STRT:aGRID%I_STOP,
     * aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)::
     * aWEIGHT

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


