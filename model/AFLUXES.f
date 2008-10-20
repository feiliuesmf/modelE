#include "rundeck_opts.h"

      MODULE AFLUXES
!@sum  AFLUXES contains some ocean arrays defined on the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, ONLY : IMA=>IM, JMA=>JM

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: NTM
#endif

! Arrays  aMO,aUO1,aVO1,aG0M,aS0M,aOGEOZ,aOGEOZ_SV,aTRMO are ocean quantities
!   after interpolation from ocean to atmospheric grid  
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: aMO, aG0, aS0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: aOGEOZ,aOGEOZ_SV
     *      , aUO1,aVO1
#ifdef TRACERS_OCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: aTRMO
#endif

! Global arrays needed for interpolation from ocean to atmospheric grid  
      REAL*8, DIMENSION(IMA,JMA,2) :: aMO_glob, aG0_glob, aS0_glob
      REAL*8, DIMENSION(IMA,JMA)   :: aUO1_glob, aVO1_glob
      REAL*8, DIMENSION(IMA,JMA)   :: aOGEOZ_glob, aOGEOZ_SV_glob 
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(IMA,JMA,1,NTM) :: aTRMO_glob
#endif

      END MODULE AFLUXES

      SUBROUTINE ALLOC_AFLUXES(grd_dum)
!@sum   Initializes AFLUXES''s arrays
!@auth  Larissa Nazarenko 
!@ver  1.0

      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE AFLUXES

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      ALLOCATE(  aMO (I_0H:I_1H,J_0H:J_1H,2), STAT = IER)
      ALLOCATE( aG0  (I_0H:I_1H,J_0H:J_1H,2), STAT = IER)
      ALLOCATE( aS0  (I_0H:I_1H,J_0H:J_1H,2), STAT = IER)

      ALLOCATE( aOGEOZ     (I_0H:I_1H, J_0H:J_1H),
     &          aOGEOZ_SV  (I_0H:I_1H, J_0H:J_1H),
     &          aUO1       (I_0H:I_1H, J_0H:J_1H), 
     &          aVO1       (I_0H:I_1H, J_0H:J_1H), 
     &   STAT=IER )

#ifdef TRACERS_OCEAN
      ALLOCATE(aTRMO(I_0H:I_1H,J_0H:J_1H,1,NTM), STAT = IER)
#endif

      END SUBROUTINE ALLOC_AFLUXES

