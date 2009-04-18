#include "rundeck_opts.h"

      MODULE INT_AG2OG_MOD

!@sum INT_AG2OG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from atm. to the ocean grid 
!@auth Larissa Nazarenko
!@ver  1.0
      PRIVATE
      PUBLIC INT_AG2OG

      Interface INT_AG2OG
      Module Procedure INT_AG2OG_2D
      End Interface

      contains

      SUBROUTINE INT_AG2OG_2D(aA,oA,aWEIGHT,aFACTOR,oFACTOR)

!@sum INT_AG2OG_2D is for conversion 2D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_DATA, UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER

      REAL*8, INTENT(IN)  :: aA(aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, DIMENSION(aIM,aJM) :: aFtemp, aFactor, aWEIGHT
      REAL*8, DIMENSION(oIM,oJM) :: oFtemp, oFactor

      REAL*8, ALLOCATABLE :: aA_glob(:,:), oA_glob(:,:)

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(oIM,oJM), STAT = IER)
      end if

C***  Gather 2D array on atmospheric grid into the global array

      CALL PACK_DATA (aGRID, aA, aA_glob)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        aFtemp(:,:) = aA_glob(:,:)*aFactor(:,:)
        aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
        call HNTR8P (aWEIGHT, aFtemp, oFtemp)
        oA_glob(:,:) = oFtemp(:,:)*oFactor(:,:)
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_DATA (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob)
      end if

      RETURN
      END SUBROUTINE INT_AG2OG_2D

      SUBROUTINE INT_AG2OG_3D(aA,oA,aWEIGHT,aFACTOR,oFACTOR,NT)

!@sum INT_AG2OG_2D is for conversion 2D arrays from atm. to the ocean grid 

!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm

      USE OCEAN, only : oDLATM=>DLATM
      Use GEOM,  only : aDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only : aGRID=>GRID, AM_I_ROOT
     *                           , PACK_COLUMN, UNPACK_COLUMN
      Use OCEANR_DIM,       only : oGRID

      IMPLICIT NONE

      INTEGER :: IER, NT,N

      REAL*8, INTENT(IN)  :: 
     *  aA(NT,aIM,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) :: 
     *  oA(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)

      REAL*8, DIMENSION(aIM,aJM) :: aFtemp, aFactor, aWEIGHT
      REAL*8, DIMENSION(oIM,oJM) :: oFtemp, oFactor

      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), oA_glob(:,:,:)

      if(AM_I_ROOT()) then
        ALLOCATE(aA_glob(NT,aIM,aJM), STAT = IER)
        ALLOCATE(oA_glob(NT,oIM,oJM), STAT = IER)
      end if

C***  Gather 2D array on atmospheric grid into the global array

      CALL PACK_COLUMN (aGRID, aA, aA_glob)

C***  Interpolate aA_glob from atmospheric grid to ocean grid 

      if(AM_I_ROOT()) then

        call HNTR80 (aIM,aJM,0.d0,aDLATM, oIM,oJM,0.d0,oDLATM, 0.d0)

        DO N=1,NT
          aFtemp(:,:) = aA_glob(N,:,:)*aFactor(:,:)
          aFtemp(2:aIM,aJM) = aFtemp(1,aJM)
          call HNTR8P (aWEIGHT, aFtemp, oFtemp)
          oA_glob(N,:,:) = oFtemp(:,:)*oFactor(:,:)
        END DO
      end if

C***  Scatter global array oA_glob to the ocean grid

      CALL UNPACK_COLUMN (ogrid, oA_glob, oA)

      if(AM_I_ROOT()) then
        DEALLOCATE(aA_glob, oA_glob)
      end if

      RETURN
      END SUBROUTINE INT_AG2OG_3D

      END MODULE INT_AG2OG_MOD

#ifndef CUBE_GRID
      SUBROUTINE INT_AG2OG_precip(aRSI_glob)

!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE OCEAN, only : oDXYPO=>DXYPO,oDLATM=>DLATM
      Use GEOM,  only : aDXYP, aDLATM=>DLATM

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

      INTEGER J

      REAL*8, INTENT(IN), DIMENSION(aIM,aJM) :: aRSI_glob

      REAL*8, DIMENSION(aIM,aJM) :: aFactor, aWEIGHT
      REAL*8, DIMENSION(oIM,oJM) :: oFactor

      aWEIGHT(:,:) = 1.- aRSI_glob(:,:)  !!  open ocean fraction
      aFactor(:,:) = 1.d0
      oFactor(:,:) = 1.d0
      CALL INT_AG2OG(aPREC,oPREC, aWEIGHT, aFactor,oFactor)

      DO J = 1,oJM
        oFactor(:,J) = oDXYPO(J)           
      END DO
      CALL INT_AG2OG(aEPREC,oEPREC, aWEIGHT, aFactor,oFactor)
      
      aWEIGHT(:,:) = aRSI_glob(:,:)  
      aFactor(:,:) = 1.d0
      oFactor(:,:) = 1.d0
      CALL INT_AG2OG(aRUNPSI,oRUNPSI, aWEIGHT, aFactor,oFactor)

      DO J = 1,oJM
        oFactor(:,J) = oDXYPO(J)           
      END DO
      CALL INT_AG2OG(aSRUNPSI,oSRUNPSI, aWEIGHT, aFactor,oFactor)

      CALL INT_AG2OG(aERUNPSI,oERUNPSI, aWEIGHT, aFactor,oFactor)

      aWEIGHT(:,:) = 1.d0
      aFactor(:,:) = 1.d0
      oFactor(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT, aFactor,oFactor)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)

      aWEIGHT(:,:) = 1.- aRSI_glob(:,:) 
      DO J = 1,aJM
        aFactor(:,J) = 1./aDXYP(J)           
      END DO
      DO J = 1,oJM
        oFactor(:,J) = oDXYPO(J)           
      END DO
      CALL INT_AG2OG(aTRPREC,oTRPREC, aWEIGHT, aFactor,oFactor,NTM)
      
      aWEIGHT(:,:) = aRSI_glob(:,:)  
      aFactor(:,:) = 1.d0
      CALL INT_AG2OG(aTRUNPSI,oTRUNPSI, aWEIGHT, aFactor,oFactor,NTM)
#endif

      RETURN
      END SUBROUTINE INT_AG2OG_precip
#endif  
!!  CUBE_GRID

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

#ifndef CUBE_GRID
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
      Use GEOM,  only : aIMAXJ=>IMAXJ,aDLATM=>DLATM
     *     , aCOSI=>COSIP,aSINI=>SINIP

      IMPLICIT NONE

      integer I,J,L, NT ,im1

      REAL*8, DIMENSION(IMA,JMA) :: aFtemp
      REAL*8, DIMENSION(IMO,JMO) :: oFtemp, oONES, oFweight
      REAL*8  oVOsp, oVOnp
      REAL*8  aUO1sp, aVO1sp, aUO1np, aVO1np
      REAL*8  SUM_oG0M, SUM_oFtemp, SUM_aG0, diff

!      write (555,*) ' INT_OG2AG#1: dLATMO,aDLATM= ',oDLATM,aDLATM

      oONES(:,:) = 1.d0

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

      RETURN
      END SUBROUTINE INT_OG2AG
#endif


      SUBROUTINE gather_ocean1

!@sum  gather_ocean1  gathers necessary arrays on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      use domain_decomp_1d, only: agrid=>grid, pack_data
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

      CALL PACK_DATA(agrid, aFOCEAN_loc, aFOCEAN)

      CALL PACK_DATA(ogrid,   MO   ,    MO_glob)
      CALL PACK_DATA(ogrid,   UO   ,    UO_glob)
      CALL PACK_DATA(ogrid,   VO   ,    VO_glob)

      CALL PACK_DATA(ogrid,OGEOZ   , OGEOZ_glob)
      CALL PACK_DATA(ogrid,OGEOZ_SV,OGEOZ_SV_glob)

      CALL PACK_DATA(ogrid,  G0M   ,   G0M_glob)
      CALL PACK_DATA(ogrid,  S0M   ,   S0M_glob)
#ifdef TRACERS_OCEAN
      CALL PACK_DATA(ogrid,  TRMO  ,  TRMO_glob)
#endif

      RETURN
      END SUBROUTINE gather_ocean1

      SUBROUTINE scatter_ocean1

!@sum  scatter_ocean1  scatters necessary arrays on the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0

      use domain_decomp_1d, only: agrid=>grid,unpack_data

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

#ifndef CUBE_GRID
      SUBROUTINE INT_AG2OG_oceans(
     *       aRSI_glob, aFLOWO_glob,aEFLOWO_glob
     *     , aRUNOSI_glob,aSRUNOSI_glob,aERUNOSI_glob
     *     , aMELTI_glob,aEMELTI_glob,aSMELTI_glob
     *     , aAPRESS_glob,aSOLAR_glob,aE0_glob,aEVAPOR_glob
     *     , aDMUA_glob, aDMVA_glob, aDMUI_glob, aDMVI_glob
     *     , aGMELT_glob, aEGMELT_glob
     *     , oRSI_glob, oFLOWO_glob,oEFLOWO_glob
     *     , oRUNOSI_glob,oSRUNOSI_glob,oERUNOSI_glob
     *     , oMELTI_glob,oEMELTI_glob,oSMELTI_glob
     *     , oAPRESS_glob,oSOLAR_glob,oE0_glob,oEVAPOR_glob
     *     , oDMUA_glob, oDMVA_glob, oDMUI_glob, oDMVI_glob
     *     , oGMELT_glob, oEGMELT_glob
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO_glob, aTRUNOSI_glob
     *     , aTRMELTI_glob, aTREVAPOR_glob
     *     , oTRFLOWO_glob, oTRUNOSI_glob
     *     , oTRMELTI_glob, oTREVAPOR_glob
     *     , aTRGMELT_glob, oTRGMELT_glob
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP_glob, oTRDRYDEP_glob
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
     *     , aTRGASEX_glob, oTRGASEX_glob
#endif
#ifdef OBIO_RAD_coupling
     *     , avisdir_glob,ovisdir_glob
     *     , asrvissurf_glob
     *     , avisdif_glob,ovisdif_glob
     *     , anirdir_glob,onirdir_glob
     *     , anirdif_glob,onirdif_glob
#endif
#ifdef TRACERS_OceanBiology
     *     , aCOSZ1_glob, osolz_glob
     *     , awind_glob,  owind_glob
#endif
     *     )
!@sum  INT_AG2OG is for interpolation of arrays from subr.GROUND_OC
!!      from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : IMA=>IM,JMA=>JM
      USE OCEAN, only : IMO=>IM,JMO=>JM, IVSPO=>IVSP, IVNPO=>IVNP
      USE MODEL_COM, only : IVSPA=>IVSP, IVNPA=>IVNP

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#if defined (TRACERS_OceanBiology) && !defined (TRACERS_GASEXCH_ocean)
      USE OCN_TRACER_COM, only: ntm
#else
      USE TRACER_COM, only: ntm
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE MODEL_COM, only: nstep=>itime
      USE TRACER_GASEXCH_COM, only: tracflx_glob
#endif

      USE OCEAN, only : oDXYPO=>DXYPO, oFOCEAN=>FOCEAN, oIMAXJ=>IMAXJ,
     *     oDLATM=>DLATM, oSINI=>SINIC, oCOSI=>COSIC,
     *     oCOSU=>COSU,oSINU=>SINU
      Use GEOM,  only : aDXYP, aIMAXJ=>IMAXJ,
     *     aDLATM=>DLATM, aSINI=>SINIP, aCOSI=>COSIP,
     *     aSINU=>SINU, aCOSU=>COSU

      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob

      IMPLICIT NONE

      INTEGER I,J, N

!!!   Global arrays on atmospheric grid

      INTEGER, PARAMETER :: NSTYPE=4

      REAL*8, INTENT(IN), DIMENSION(IMA,JMA) ::
     *        aRSI_glob, aFLOWO_glob,aEFLOWO_glob
     *      , aRUNOSI_glob, aSRUNOSI_glob, aERUNOSI_glob
     *      , aMELTI_glob,aEMELTI_glob,aSMELTI_glob, aAPRESS_glob
     *      , aDMUI_glob, aDMVI_glob, aGMELT_glob, aEGMELT_glob
      REAL*8, INTENT(IN), DIMENSION(3,IMA,JMA) :: aSOLAR_glob
      REAL*8, INTENT(IN), DIMENSION(IMA,JMA,NSTYPE) ::
     *        aE0_glob, aEVAPOR_glob
     *      , aDMUA_glob, aDMVA_glob
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      REAL*8, INTENT(IN), DIMENSION(NTM,IMA,JMA) :: aTRFLOWO_glob
     *      , aTRUNOSI_glob, aTRMELTI_glob, aTRGMELT_glob
      REAL*8, INTENT(IN), DIMENSION(NTM,NSTYPE,IMA,JMA) ::
     *        aTREVAPOR_glob
#ifdef TRACERS_DRYDEP
      REAL*8, INTENT(IN), DIMENSION(NTM,NSTYPE,IMA,JMA) ::
     *        aTRDRYDEP_glob
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      REAL*8, INTENT(IN), DIMENSION(NTM,NSTYPE,IMA,JMA) :: aTRGASEX_glob
#endif
#ifdef TRACERS_OceanBiology
      REAL*8, INTENT(IN), DIMENSION(IMA,JMA) :: avisdir_glob
     *    ,asrvissurf_glob
     *    ,avisdif_glob,anirdir_glob,anirdif_glob
#endif
#ifdef TRACERS_OceanBiology
      REAL*8, INTENT(IN), DIMENSION(IMA,JMA) :: aCOSZ1_glob,awind_glob
#endif

!!!   Global arrays on ocean grid

      REAL*8, INTENT(INOUT), DIMENSION(IMO,JMO) ::
     *        oRSI_glob, oFLOWO_glob,oEFLOWO_glob
     *      , oRUNOSI_glob, oSRUNOSI_glob, oERUNOSI_glob
     *      , oMELTI_glob,oEMELTI_glob,oSMELTI_glob, oAPRESS_glob
     *      , oDMUI_glob, oDMVI_glob, oGMELT_glob, oEGMELT_glob
      REAL*8, INTENT(INOUT), DIMENSION(3,IMO,JMO) :: oSOLAR_glob
      REAL*8, INTENT(INOUT), DIMENSION(IMO,JMO,1) ::
     *        oE0_glob, oEVAPOR_glob
     *      , oDMUA_glob, oDMVA_glob
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT), DIMENSION(NTM,IMO,JMO) :: oTRFLOWO_glob
     *      , oTRUNOSI_glob, oTRMELTI_glob, oTRGMELT_glob
      REAL*8, INTENT(INOUT), DIMENSION(NTM,1,IMO,JMO) ::
     *        oTREVAPOR_glob
#ifdef TRACERS_DRYDEP
      REAL*8, INTENT(INOUT), DIMENSION(NTM,1,IMO,JMO) ::
     *        oTRDRYDEP_glob
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      REAL*8, INTENT(INOUT), DIMENSION(NTM,1,IMO,JMO) :: oTRGASEX_glob
#endif
#ifdef OBIO_RAD_coupling
      REAL*8, INTENT(INOUT), DIMENSION(IMO,JMO) :: ovisdir_glob,
     *   ovisdif_glob,onirdir_glob,onirdif_glob
#endif
#ifdef TRACERS_OceanBiology
      REAL*8, INTENT(INOUT), DIMENSION(IMO,JMO) :: osolz_glob,owind_glob
#endif
      REAL*8, DIMENSION(IMA,JMA) :: aFtemp, aONES, aOOCN
      REAL*8, DIMENSION(IMO,JMO) :: oFtemp
      REAL*8  aDMUAsp, aDMVAsp, aDMUAnp, aDMVAnp
      REAL*8  oDMUAsp, oDMVAsp, oDMUAnp, oDMVAnp
      REAL*8  aDMUIsp, aDMVIsp, aDMUInp, aDMVInp
      REAL*8  oDMUIsp, oDMVIsp, oDMUInp, oDMVInp
      REAL*8  diff, SUM_oFtemp, SUM_aFtemp, SUM_aORIG

      aONES(:,:) = 1.d0
      aOOCN(:,:) = 1.d0 - aRSI_glob(:,:)

      call HNTR80 (IMA,JMA,0.d0,aDLATM, IMO,JMO,0.d0,oDLATM, 0.d0)

      aFtemp(:,:) = aRSI_glob(:,:)
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oRSI_glob(:,:) = oFtemp(:,:)

      aFtemp(:,:) = 0.
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aFLOWO_glob(I,J)/(aDXYP(I,J)*aFOCEAN(I,J))  !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oFLOWO_glob(:,:) = oFtemp(:,:)                                  !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aEFLOWO_glob(I,J)/(aDXYP(I,J)*aFOCEAN(I,J)) !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oEFLOWO_glob(:,:) = oFtemp(:,:)                                 !  J/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aMELTI_glob(I,J)/(aDXYP(I,J)*aFOCEAN(I,J))  !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oMELTI_glob(:,:) = oFtemp(:,:)                                  !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aEMELTI_glob(I,J)/(aDXYP(I,J)*aFOCEAN(I,J)) !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oEMELTI_glob(:,:) = oFtemp(:,:)                                 !  J/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aSMELTI_glob(I,J)/(aDXYP(I,J)*aFOCEAN(I,J)) !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oSMELTI_glob(:,:) = oFtemp(:,:)                                 !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aRUNOSI_glob(I,J)                           !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aRSI_glob, aFtemp, oFtemp)
      oRUNOSI_glob(:,:) = oFtemp(:,:)                                 !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aERUNOSI_glob(I,J)                          !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aRSI_glob, aFtemp, oFtemp)
      oERUNOSI_glob(:,:) = oFtemp(:,:)                                !  J/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aSRUNOSI_glob(I,J)                          !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aRSI_glob, aFtemp, oFtemp)
      oSRUNOSI_glob(:,:) = oFtemp(:,:)                                !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aE0_glob(I,J,1)                             !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aOOCN, aFtemp, oFtemp)
      oE0_glob(:,:,1) = oFtemp(:,:)                                   !  J/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aEVAPOR_glob(I,J,1)                         !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aOOCN, aFtemp, oFtemp)
      oEVAPOR_glob(:,:,1) = oFtemp(:,:)                               !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aSOLAR_glob(1,I,J)                          !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aOOCN, aFtemp, oFtemp)
      oSOLAR_glob(1,:,:) = oFtemp(:,:)                                !  J/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aSOLAR_glob(3,I,J)                          !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aRSI_glob, aFtemp, oFtemp)
      oSOLAR_glob(3,:,:) = oFtemp(:,:)                                !  J/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aAPRESS_glob(I,J)                           !  Pa
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oAPRESS_glob(:,:) = oFtemp(:,:)                                 !  Pa

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aGMELT_glob(I,J)/(aDXYP(I,J)*aFOCEAN(I,J))  !  kg/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oGMELT_glob(:,:) = oFtemp(:,:)                                  !  kg/m^2

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aEGMELT_glob(I,J)/aDXYP(I,J)                !  J/m^2
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      DO J = 1,JMO
        oEGMELT_glob(:,J) = oFtemp(:,J)*oDXYPO(J)                     !  J
      END DO

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTRFLOWO_glob(N,I,J)/(aDXYP(I,J)*aFOCEAN(I,J)) !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oTRFLOWO_glob(N,:,:) = oFtemp(:,:)                                 !  kg/m^2
      END DO

      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTRMELTI_glob(N,I,J)/(aDXYP(I,J)*aFOCEAN(I,J)) !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oTRMELTI_glob(N,:,:) = oFtemp(:,:)                                 !  kg/m^2
      END DO

      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTRUNOSI_glob(N,I,J)                           !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aRSI_glob, aFtemp, oFtemp)
      oTRUNOSI_glob(N,:,:) = oFtemp(:,:)                                 !  kg/m^2
      END DO

      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTREVAPOR_glob(N,1,I,J)                        !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aOOCN, aFtemp, oFtemp)
      oTREVAPOR_glob(N,1,:,:) = oFtemp(:,:)                              !  kg/m^2
      END DO

      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTRGMELT_glob(N,I,J)/aDXYP(I,J)                !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      DO J = 1,JMO
        oTRGMELT_glob(N,:,J) = oFtemp(:,J)*oDXYPO(J)                     !  kg
      END DO
      END DO

#ifdef TRACERS_DRYDEP
      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTRDRYDEP_glob(N,1,I,J)                        !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oTRDRYDEP_glob(N,1,:,:) = oFtemp(:,:)                              !  kg/m^2
      END DO
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      DO N=1,NTM
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aTRGASEX_glob(N,1,I,J)                         !  kg/m^2
          ELSE
            aFtemp(I,J) = 0.
          END IF
        END DO
      END DO

      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aONES, aFtemp, oFtemp)
      oTRGASEX_glob(N,1,:,:) = oFtemp(:,:)                               !  kg/m^2
      tracflx_glob(:,:,N)=oTRGASEX_glob(N,1,:,:)
      END DO
  
#endif

#ifdef OBIO_RAD_coupling
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = avisdir_glob(I,J)
     .                  * asrvissurf_glob(I,J)
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aFOCEAN, aFtemp, oFtemp)
      ovisdir_glob(:,:) = oFtemp(:,:)

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = avisdif_glob(I,J)
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aFOCEAN, aFtemp, oFtemp)
      ovisdif_glob(:,:) = oFtemp(:,:)

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = anirdir_glob(I,J)
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aFOCEAN, aFtemp, oFtemp)
      onirdir_glob(:,:) = oFtemp(:,:)

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = anirdif_glob(I,J)
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aFOCEAN, aFtemp, oFtemp)
      onirdif_glob(:,:) = oFtemp(:,:)
#endif

#ifdef TRACERS_OceanBiology
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aCOSZ1_glob(I,J)
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aFOCEAN, aFtemp, oFtemp)
      osolz_glob(:,:) = oFtemp(:,:)
      write(*,*)'OCNDYN, solz=',osolz_glob(16,45)

      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = awind_glob(I,J)
          END IF
        END DO
      END DO
      aFtemp(2:IMA,JMA) = aFtemp(1,JMA)
      call HNTR8P (aFOCEAN, aFtemp, oFtemp)
      owind_glob(:,:) = oFtemp(:,:)
#endif

      aDMUAsp = aDMUA_glob(1,1,1)
      aDMVAsp = aDMVA_glob(1,1,1)
      aDMUAnp = aDMUA_glob(1,JMA,1)
      aDMVAnp = aDMVA_glob(1,JMA,1)

      aFtemp(:,1  ) = aDMUAsp*aCOSI(:) - aDMVAsp*aSINI(:)
      aFtemp(:,JMA) = aDMUAnp*aCOSI(:) + aDMVAnp*aSINI(:)

      DO J=2,JMA-1
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aDMUA_glob(I,J,1)                           !  kg/m s
          END IF
        END DO
      END DO
      call HNTR80 (IMA,JMA,0.d0,aDLATM, IMO,JMO,0.d0,oDLATM, 0.d0)
      call HNTR8  (aONES, aFtemp, oDMUA_glob)          !!  A-grid => A-grid

      oDMUAsp =  SUM(oDMUA_glob(:,  1,1)*oCOSI(:))*2/IMO
      oDMVAsp = -SUM(oDMUA_glob(:,  1,1)*oSINI(:))*2/IMO
      oDMUAnp =  SUM(oDMUA_glob(:,JMO,1)*oCOSI(:))*2/IMO
      oDMVAnp =  SUM(oDMUA_glob(:,JMO,1)*oSINI(:))*2/IMO

      oDMUA_glob(1,  1,1) = oDMUAsp
      oDMUA_glob(1,JMO,1) = oDMUAnp

      aFtemp(:,  1) = aDMVAsp*aCOSI(:) + aDMUAsp*aSINI(:)
      aFtemp(:,JMA) = aDMVAnp*aCOSI(:) - aDMUAnp*aSINI(:)

      DO J=2,JMA-1
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aFtemp(I,J) = aDMVA_glob(I,J,1)                           !  kg/m s
          END IF
        END DO
      END DO
      call HNTR80 (IMA,JMA,0.d0,aDLATM, IMO,JMO,0.d0,oDLATM, 0.d0)
      call HNTR8  (aONES, aFtemp, oDMVA_glob)          !!  A-grid => A-grid

      oDMVA_glob(1,  1,1) = oDMVAsp
      oDMVA_glob(1,JMO,1) = oDMVAnp

      aDMUIsp = 0.0  !!  aDMUI_glob(  IMA,  1)
      aDMVIsp = 0.0  !!  aDMVI_glob(IVSPA,  1)
      aDMUInp = 0.0  !!  aDMUI_glob(  IMA,JMA)
      aDMVInp = 0.0  !!  aDMVI_glob(IVNPA,JMA)

      aFtemp(:,:)   = aDMUI_glob(:,:)
      aFtemp(:,JMA) = 0.0
      call HNTR80 (IMA,JMA,.5d0,aDLATM, IMO,JMO,.5d0,oDLATM, 0.d0)
      call HNTR8  (aONES, aFtemp, oDMUI_glob)          !!  U-grid => U-grid

      oDMUIsp = 0.0  !!   SUM(oDMUI_glob(:,  1)*oCOSU(:))*2/IMO
      oDMVIsp = 0.0  !!  -SUM(oDMUI_glob(:,  1)*oSINU(:))*2/IMO
      oDMUInp = 0.0  !!   SUM(oDMUI_glob(:,JMO)*oCOSU(:))*2/IMO
      oDMVInp = 0.0  !!   SUM(oDMUI_glob(:,JMO)*oSINU(:))*2/IMO

      oDMUI_glob(  IMO,  1) = oDMUIsp
      oDMUI_glob(IVSPO,  1) = oDMUIsp
      oDMUI_glob(  IMO,JMO) = oDMUInp
      oDMUI_glob(IVNPO,  1) = oDMUInp

      call HNTR80 (IMA,JMA-1,0.d0,aDLATM, IMO,JMO-1,0.d0,oDLATM, 0.d0)
      call HNTR8  (aONES, aDMVI_glob, oDMVI_glob)      !!  V-grid => V-grid
      oDMVI_glob(:,JMO) = 0.0                          !!  should not be used

      RETURN
      END SUBROUTINE INT_AG2OG_oceans
#endif

#ifndef CUBE_GRID
      SUBROUTINE INT_OG2AG_oceans(
     *       oDMSI_glob, oDHSI_glob, oDSSI_glob
     *     , aDMSI_glob, aDHSI_glob, aDSSI_glob
#ifdef TRACERS_OCEAN
     *     , oDTRSI_glob
#endif
#ifdef TRACERS_ON
     *     , aDTRSI_glob
#endif
     *     )
!@sum  INT_OG2AG is for interpolation of arrays from subr.GROUND_OC
!!      from ocean grid to the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : IMA=>IM,JMA=>JM
      USE OCEAN, only : IMO=>IM,JMO=>JM, oFOCEAN=>FOCEAN,
     *     oDXYPO=>DXYPO, oIMAXJ=>IMAXJ, oDLATM=>DLATM

      USE AFLUXES, only : aFOCEAN=>aFOCEAN_glob

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE GEOM,  only : aIMAXJ=>IMAXJ, aDLATM=>DLATM

      IMPLICIT NONE

      INTEGER I,J, N, NT

!!!   Global arrays on ocean grid

      REAL*8, INTENT(IN), DIMENSION(2,IMO,JMO) ::
     *     oDMSI_glob, oDHSI_glob, oDSSI_glob
#ifdef TRACERS_OCEAN
      REAL*8, INTENT(IN), DIMENSION(NTM,2,IMO,JMO) :: oDTRSI_glob
#endif

!!!   Global arrays on atmospheric grid

      REAL*8, INTENT(INOUT), DIMENSION(2,IMA,JMA) ::
     *     aDMSI_glob, aDHSI_glob, aDSSI_glob
#ifdef TRACERS_ON
      REAL*8, INTENT(INOUT), DIMENSION(NTM_ATM,2,IMA,JMA) :: aDTRSI_glob
#endif

      REAL*8, DIMENSION(IMA,JMA) :: aFtemp
      REAL*8, DIMENSION(IMO,JMO) :: oFtemp

      call HNTR80 (IMO,JMO,0.d0,oDLATM, IMA,JMA,0.d0,aDLATM, 0.d0)

      DO N=1,2
      DO J=1,JMO
        DO I=1,oIMAXJ(J)
          IF (oFOCEAN(I,J).gt.0.) THEN
            oFtemp(I,J) = oDMSI_glob(N,I,J)                 !  kg/m^2
          END IF
        END DO
      END DO
      oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
      call HNTR8P (oFOCEAN, oFtemp, aFtemp)
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aDMSI_glob(N,I,J) = aFtemp(I,J)                 !  kg/m^2
          END IF
        END DO
      END DO

      DO J=1,JMO
        DO I=1,oIMAXJ(J)
          IF (oFOCEAN(I,J).gt.0.) THEN
            oFtemp(I,J) = oDHSI_glob(N,I,J)                 !  J/m^2
          END IF
        END DO
      END DO
      oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
      call HNTR8P (oFOCEAN, oFtemp, aFtemp)
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aDHSI_glob(N,I,J) = aFtemp(I,J)                 !  J/m^2
          END IF
        END DO
      END DO

      DO J=1,JMO
        DO I=1,oIMAXJ(J)
          IF (oFOCEAN(I,J).gt.0.) THEN
            oFtemp(I,J) = oDSSI_glob(N,I,J)                 !  kg/m^2
          END IF
        END DO
      END DO
      oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
      call HNTR8P (oFOCEAN, oFtemp, aFtemp)
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aDSSI_glob(N,I,J) = aFtemp(I,J)                 !  kg/m^2
          END IF
        END DO
      END DO

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (NTM == NTM_ATM) THEN
      DO NT=1,NTM
      DO J=1,JMO
        DO I=1,oIMAXJ(J)
          IF (oFOCEAN(I,J).gt.0.) THEN
            oFtemp(I,J) = oDTRSI_glob(NT,N,I,J)             !  kg/m^2
          END IF
        END DO
      END DO
      oFtemp(2:IMO,JMO) = oFtemp(1,JMO)
      call HNTR8P (oFOCEAN, oFtemp, aFtemp)
      DO J=1,JMA
        DO I=1,aIMAXJ(J)
          IF (aFOCEAN(I,J).gt.0.) THEN
            aDTRSI_glob(NT,N,I,J) = aFtemp(I,J)             !  kg/m^2
          END IF
        END DO
      END DO
      END DO
      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif
      END DO   !N-loop

      END SUBROUTINE INT_OG2AG_oceans
#endif

      SUBROUTINE AG2OG_oceans
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : IMA=>IM, JMA=>JM

      USE OCEAN, only : IMO=>IM,JMO=>JM

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#ifdef TRACERS_OceanBiology
      USE OCN_TRACER_COM, only: ntm
#else
      USE TRACER_COM, only: ntm
#endif
#endif

      USE DOMAIN_DECOMP_1D, only : agrid=>grid, AM_I_ROOT
     *     , PACK_DATA,UNPACK_DATA, PACK_COLUMN, UNPACK_COLUMN
     *     , PACK_BLOCK, UNPACK_BLOCK
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

      IMPLICIT NONE

!!!   Global arrays on atmospheric grid

      INTEGER, PARAMETER :: NSTYPE=4
      REAL*8, DIMENSION(IMA,JMA) :: aRSI_glob, aFLOWO_glob,aEFLOWO_glob
     *      , aRUNOSI_glob, aSRUNOSI_glob, aERUNOSI_glob, aAPRESS_glob
     *      , aMELTI_glob, aEMELTI_glob, aSMELTI_glob
     *      , aDMUI_glob, aDMVI_glob, aGMELT_glob, aEGMELT_glob
      REAL*8, DIMENSION(3,IMA,JMA) :: aSOLAR_glob
      REAL*8, DIMENSION(IMA,JMA,NSTYPE) :: aE0_glob, aEVAPOR_glob
     *      , aDMUA_glob, aDMVA_glob
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,IMA,JMA) :: aTRFLOWO_glob
     *      , aTRUNOSI_glob, aTRMELTI_glob, aTRGMELT_glob
      REAL*8, DIMENSION(NTM,NSTYPE,IMA,JMA) :: aTREVAPOR_glob
#ifdef TRACERS_DRYDEP
      REAL*8, DIMENSION(NTM,NSTYPE,IMA,JMA) :: aTRDRYDEP_glob
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      REAL*8, DIMENSION(NTM,NSTYPE,IMA,JMA) :: aTRGASEX_glob
#endif
#ifdef OBIO_RAD_coupling
      REAL*8, DIMENSION(IMA,JMA) :: avisdir_glob
     *                             ,asrvissurf_glob,avisdif_glob
     *                             ,anirdir_glob,anirdif_glob
#endif
#ifdef TRACERS_OceanBiology
      REAL*8, DIMENSION(IMA,JMA) :: aCOSZ1_glob, awind_glob
#endif

!!!   Global arrays on ocean grid

      REAL*8, DIMENSION(IMO,JMO) :: oRSI_glob, oFLOWO_glob,oEFLOWO_glob
     *      , oRUNOSI_glob, oSRUNOSI_glob, oERUNOSI_glob, oAPRESS_glob
     *      , oMELTI_glob, oEMELTI_glob, oSMELTI_glob
     *      , oDMUI_glob, oDMVI_glob, oGMELT_glob, oEGMELT_glob
      REAL*8, DIMENSION(3,IMO,JMO) :: oSOLAR_glob
      REAL*8, DIMENSION(IMO,JMO,1) :: oE0_glob, oEVAPOR_glob
     *      , oDMUA_glob, oDMVA_glob
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM,IMO,JMO) :: oTRFLOWO_glob
     *      , oTRUNOSI_glob, oTRMELTI_glob, oTRGMELT_glob
      REAL*8, DIMENSION(NTM,1,IMO,JMO) :: oTREVAPOR_glob
#ifdef TRACERS_DRYDEP
      REAL*8, DIMENSION(NTM,1,IMO,JMO) :: oTRDRYDEP_glob
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      REAL*8, DIMENSION(NTM,1,IMO,JMO) :: oTRGASEX_glob
#endif
#ifdef OBIO_RAD_coupling
      REAL*8, DIMENSION(IMO,JMO) :: ovisdir_glob,
     *  ovisdif_glob,onirdir_glob,onirdif_glob
#endif
#ifdef TRACERS_OceanBiology
      REAL*8, DIMENSION(IMO,JMO) :: osolz_glob, owind_glob
#endif


C***  Gather arrays on atmospheric grid into the global arrays

      CALL PACK_DATA  (agrid,     aRSI,     aRSI_glob)
      CALL PACK_DATA  (agrid,   aFLOWO,   aFLOWO_glob)
      CALL PACK_DATA  (agrid,  aEFLOWO,  aEFLOWO_glob)
      CALL PACK_DATA  (agrid,  aRUNOSI,  aRUNOSI_glob)
      CALL PACK_DATA  (agrid, aSRUNOSI, aSRUNOSI_glob)
      CALL PACK_DATA  (agrid, aERUNOSI, aERUNOSI_glob)
      CALL PACK_DATA  (agrid,  aAPRESS,  aAPRESS_glob)
      CALL PACK_DATA  (agrid,   aMELTI,   aMELTI_glob)
      CALL PACK_DATA  (agrid,  aEMELTI,  aEMELTI_glob)
      CALL PACK_DATA  (agrid,  aSMELTI,  aSMELTI_glob)
      CALL PACK_COLUMN(agrid,   aSOLAR,   aSOLAR_glob)
      CALL PACK_DATA  (agrid,      aE0,      aE0_glob)
      CALL PACK_DATA  (agrid,  aEVAPOR,  aEVAPOR_glob)
      CALL PACK_DATA  (agrid,    aDMUA,    aDMUA_glob)
      CALL PACK_DATA  (agrid,    aDMVA,    aDMVA_glob)
      CALL PACK_DATA  (agrid,    aDMUI,    aDMUI_glob)
      CALL PACK_DATA  (agrid,    aDMVI,    aDMVI_glob)
      CALL PACK_DATA  (agrid,   aGMELT,   aGMELT_glob)
      CALL PACK_DATA  (agrid,  aEGMELT,  aEGMELT_glob)
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      CALL PACK_COLUMN(agrid,  aTRGMELT,  aTRGMELT_glob)
      CALL PACK_COLUMN(agrid,  aTRFLOWO,  aTRFLOWO_glob)
      CALL PACK_COLUMN(agrid,  aTRUNOSI,  aTRUNOSI_glob)
      CALL PACK_COLUMN(agrid,  aTRMELTI,  aTRMELTI_glob)
      CALL PACK_BLOCK (agrid, aTREVAPOR, aTREVAPOR_glob)
#ifdef TRACERS_DRYDEP
      CALL PACK_BLOCK (agrid, aTRDRYDEP, aTRDRYDEP_glob)
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      CALL PACK_BLOCK (agrid,  aTRGASEX,  aTRGASEX_glob)
#endif
#ifdef OBIO_RAD_coupling
      CALL PACK_DATA  (agrid, avisdir, avisdir_glob)
      CALL PACK_DATA  (agrid, asrvissurf, asrvissurf_glob)
      CALL PACK_DATA  (agrid, avisdif, avisdif_glob)
      CALL PACK_DATA  (agrid, anirdir, anirdir_glob)
      CALL PACK_DATA  (agrid, anirdif, anirdif_glob)
#endif
#ifdef TRACERS_OceanBiology
      CALL PACK_DATA  (agrid,  aCOSZ1,  aCOSZ1_glob)
      CALL PACK_DATA  (agrid,   awind,   awind_glob)
#endif
C*
C***  Do interpolation to the ocean grid
C*
      if(AM_I_ROOT()) then
#ifndef CUBE_GRID
        call INT_AG2OG_oceans(
     *       aRSI_glob, aFLOWO_glob,aEFLOWO_glob
     *     , aRUNOSI_glob,aSRUNOSI_glob,aERUNOSI_glob
     *     , aMELTI_glob,aEMELTI_glob,aSMELTI_glob
     *     , aAPRESS_glob,aSOLAR_glob,aE0_glob,aEVAPOR_glob
     *     , aDMUA_glob, aDMVA_glob, aDMUI_glob, aDMVI_glob
     *     , aGMELT_glob, aEGMELT_glob
     *     , oRSI_glob, oFLOWO_glob,oEFLOWO_glob
     *     , oRUNOSI_glob,oSRUNOSI_glob,oERUNOSI_glob
     *     , oMELTI_glob,oEMELTI_glob,oSMELTI_glob
     *     , oAPRESS_glob,oSOLAR_glob,oE0_glob,oEVAPOR_glob
     *     , oDMUA_glob, oDMVA_glob, oDMUI_glob, oDMVI_glob
     *     , oGMELT_glob, oEGMELT_glob
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO_glob, aTRUNOSI_glob
     *     , aTRMELTI_glob, aTREVAPOR_glob
     *     , oTRFLOWO_glob, oTRUNOSI_glob
     *     , oTRMELTI_glob, oTREVAPOR_glob
     *     , aTRGMELT_glob, oTRGMELT_glob
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP_glob, oTRDRYDEP_glob
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
     *     , aTRGASEX_glob, oTRGASEX_glob
#endif
#ifdef OBIO_RAD_coupling
     *     , avisdir_glob,ovisdir_glob
     *     , asrvissurf_glob
     *     , avisdif_glob,ovisdif_glob
     *     , anirdir_glob,onirdir_glob
     *     , anirdif_glob,onirdif_glob
#endif
#ifdef TRACERS_OceanBiology
     *     , aCOSZ1_glob, osolz_glob
     *     ,  awind_glob, owind_glob
#endif
     *     )
#endif
      end if

C***  Scatter the arrays to the ocean grid

      CALL UNPACK_DATA  (ogrid,     oRSI_glob,     oRSI)
      CALL UNPACK_DATA  (ogrid,   oFLOWO_glob,   oFLOWO)
      CALL UNPACK_DATA  (ogrid,  oEFLOWO_glob,  oEFLOWO)
      CALL UNPACK_DATA  (ogrid,  oRUNOSI_glob,  oRUNOSI)
      CALL UNPACK_DATA  (ogrid, oSRUNOSI_glob, oSRUNOSI)
      CALL UNPACK_DATA  (ogrid, oERUNOSI_glob, oERUNOSI)
      CALL UNPACK_DATA  (ogrid,  oAPRESS_glob,  oAPRESS)
      CALL UNPACK_DATA  (ogrid,   oMELTI_glob,   oMELTI)
      CALL UNPACK_DATA  (ogrid,  oEMELTI_glob,  oEMELTI)
      CALL UNPACK_DATA  (ogrid,  oSMELTI_glob,  oSMELTI)
      CALL UNPACK_COLUMN(ogrid,   oSOLAR_glob,   oSOLAR)
      CALL UNPACK_DATA  (ogrid,      oE0_glob,      oE0)
      CALL UNPACK_DATA  (ogrid,  oEVAPOR_glob,  oEVAPOR)
      CALL UNPACK_DATA  (ogrid,    oDMUA_glob,    oDMUA)
      CALL UNPACK_DATA  (ogrid,    oDMVA_glob,    oDMVA)
      CALL UNPACK_DATA  (ogrid,    oDMUI_glob,    oDMUI)
      CALL UNPACK_DATA  (ogrid,    oDMVI_glob,    oDMVI)
      CALL UNPACK_DATA  (ogrid,   oGMELT_glob,   oGMELT)
      CALL UNPACK_DATA  (ogrid,  oEGMELT_glob,  oEGMELT)
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      CALL UNPACK_COLUMN(ogrid,  oTRGMELT_glob,  oTRGMELT)
      CALL UNPACK_COLUMN(ogrid,  oTRFLOWO_glob,  oTRFLOWO)
      CALL UNPACK_COLUMN(ogrid,  oTRUNOSI_glob,  oTRUNOSI)
      CALL UNPACK_COLUMN(ogrid,  oTRMELTI_glob,  oTRMELTI)
      CALL UNPACK_BLOCK (ogrid, oTREVAPOR_glob, oTREVAPOR)
#ifdef TRACERS_DRYDEP
      CALL UNPACK_BLOCK (ogrid, oTRDRYDEP_glob, oTRDRYDEP)
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      CALL UNPACK_BLOCK (ogrid,  oTRGASEX_glob,  oTRGASEX)
#endif
#ifdef OBIO_RAD_coupling
      CALL UNPACK_DATA  (ogrid,ovisdir_glob,ovisdir)
      CALL UNPACK_DATA  (ogrid,ovisdif_glob,ovisdif)
      CALL UNPACK_DATA  (ogrid,onirdir_glob,onirdir)
      CALL UNPACK_DATA  (ogrid,onirdif_glob,onirdif)
#endif
#ifdef TRACERS_OceanBiology
      CALL UNPACK_DATA  (ogrid,  osolz_glob,  osolz)
      CALL UNPACK_DATA  (ogrid,  owind_glob,  owind)
#endif

      RETURN
      END SUBROUTINE AG2OG_oceans

      SUBROUTINE OG2AG_oceans
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : IMA=>IM, JMA=>JM

      USE OCEAN, only : IMO=>IM,JMO=>JM, FOCEAN, IMAXJ

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE DOMAIN_DECOMP_1D, only : agrid=>grid, AM_I_ROOT
     *     , PACK_COLUMN, UNPACK_COLUMN, PACK_DATA
     *     , PACK_BLOCK, UNPACK_BLOCK
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

      IMPLICIT NONE

!!!   Global arrays on atmospheric grid

      INTEGER, PARAMETER :: NSTYPE=4
      REAL*8, DIMENSION(2,IMA,JMA) :: aDMSI_glob, aDHSI_glob, aDSSI_glob
#ifdef TRACERS_ON
      REAL*8, DIMENSION(NTM_ATM,2,IMA,JMA) :: aDTRSI_glob
#endif

!!!   Global arrays on ocean grid

      REAL*8, DIMENSION(2,IMO,JMO) :: oDMSI_glob, oDHSI_glob, oDSSI_glob
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM,2,IMO,JMO) :: oDTRSI_glob
#endif

C***  Gather the arrays for the sea ice on atm. grid
C***    since they were already filled in for LAKES

      CALL PACK_DATA  (agrid, aFOCEAN_loc, aFOCEAN)

      CALL PACK_COLUMN(agrid,    aDMSI,    aDMSI_glob)
      CALL PACK_COLUMN(agrid,    aDHSI,    aDHSI_glob)
      CALL PACK_COLUMN(agrid,    aDSSI,    aDSSI_glob)
#ifdef TRACERS_ON
      CALL PACK_BLOCK  (agrid, aDTRSI, aDTRSI_glob)
#endif

C***  Gather arrays on ocean grid into the global arrays

      CALL PACK_COLUMN (ogrid, oDMSI, oDMSI_glob)
      CALL PACK_COLUMN (ogrid, oDHSI, oDHSI_glob)
      CALL PACK_COLUMN (ogrid, oDSSI, oDSSI_glob)
#ifdef TRACERS_OCEAN
      CALL PACK_BLOCK  (ogrid, oDTRSI, oDTRSI_glob)
#endif
C*
C***  Do interpolation to the atmospheric grid
C*
#ifndef CUBE_GRID
      if(AM_I_ROOT()) then
        call INT_OG2AG_oceans(
     *       oDMSI_glob, oDHSI_glob, oDSSI_glob
     *     , aDMSI_glob, aDHSI_glob, aDSSI_glob
#ifdef TRACERS_OCEAN
     *     , oDTRSI_glob
#endif
#ifdef TRACERS_ON
     *     , aDTRSI_glob
#endif
     *     )
      end if
#endif

C***  Scatter the arrays to the atmospheric grid

      CALL UNPACK_COLUMN (agrid, aDMSI_glob, aDMSI)
      CALL UNPACK_COLUMN (agrid, aDHSI_glob, aDHSI)
      CALL UNPACK_COLUMN (agrid, aDSSI_glob, aDSSI)
#ifdef TRACERS_ON
      CALL UNPACK_BLOCK  (agrid, aDTRSI_glob, aDTRSI)
#endif

      RETURN
      END SUBROUTINE OG2AG_oceans

      Subroutine HNTR80 (IMA,JMA,OFFIA,DLATA,
     *                   IMB,JMB,OFFIB,DLATB, DATMIS)
C****
C**** HNTR80 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR8.
C**** The 5 Real input values are expected to be Real*8.
C****
C**** Input: IMA = number of cells in east-west direction of grid A
C****        JMA = number of cells in north-south direction of grid A
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit Real*8 (A-H,O-Z)
      Parameter (TWOPI=6.283185307179586477d0)
      Real*8 OFFIA,DLATA, OFFIB,DLATB, DATMIS,DATMCB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMCB, INA,JNA, INB,JNB
C****
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or.
     *    IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 0
      GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 1-FMAX(IB)
  150 IB = IBp1
      IMAX(IMB) = IMAX(IMB) + IMA
C       WRITE (0,915) 'IMIN=',IMIN(1:IMB)
C       WRITE (0,915) 'IMAX=',IMAX(1:IMB)
C       WRITE (0,916) 'FMIN=',FMIN(1:IMB)
C       WRITE (0,916) 'FMAX=',FMAX(1:IMB)
C****
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 SINA(JA) = Sin (RJA*TWOPI/(360*60))
      SINA(0)  = -1
      SINA(JMA)=  1
C****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 SINB(JB) = Sin (RJB*TWOPI/(360*60))
      SINB(0)  = -1
      SINB(JMB)=  1
C****
      JMIN(1) = 1
      GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0
      GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA) - SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB) - SINA(JA-1)
  350 Continue
      JMAX(JMB) = JMA
      GMAX(JMB) = 0
C       WRITE (0,915) 'JMIN=',JMIN(1:JMB)
C       WRITE (0,915) 'JMAX=',JMAX(1:JMB)
C       WRITE (0,916) 'GMIN=',GMIN(1:JMB)
C       WRITE (0,916) 'GMAX=',GMAX(1:JMB)
      Return
C****
C**** Invalid parameters or dimensions out of range
C****
  400 Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS
      Stop 400
C****
C 915 Format (/ 1X,A5 / (20I6))
C 916 Format (/ 1X,A5 / (20F6.2))
  940 Format ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     *                    ' when integrated WTA = 0'/
     *  '0These arguments are invalid or out of range.')
      End Subroutine HNTR80

      Subroutine HNTR8 (WTA,A,B)
C****
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*8.
C****
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
C**** Interpolate the A grid onto the B grid
C****
      Do 20 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      Do 20 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0
      VALUE = 0
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      Do 10 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      If (JA==JAMIN)  G = G - GMIN(JB)
      If (JA==JAMAX)  G = G - GMAX(JB)
      Do 10 IAREV=IAMIN,IAMAX
      IA  = 1 + Mod(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1
      If (IAREV==IAMIN)  F = F - FMIN(IB)
      If (IAREV==IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
   10 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = DATMIS
      If (WEIGHT.ne.0)  B(IJB) = VALUE/WEIGHT
   20 Continue
      Return
      End Subroutine HNTR8

      Subroutine HNTR8P (WTA,A,B)
C****
C**** HNTR8P is similar to HNTR8 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*8.
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
      Call HNTR8 (WTA,A,B)
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      Do 40 JB=1,JMB,JMB-1
      BMEAN  = DATMIS
      WEIGHT = 0
      VALUE  = 0
      Do 10 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      If (B(IJB) == DATMIS)  GoTo 20
      WEIGHT = WEIGHT + 1
      VALUE  = VALUE  + B(IJB)
   10 Continue
      If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
   20 Do 30 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
   30 B(IJB) = BMEAN
   40 Continue
      Return
      End Subroutine HNTR8P

