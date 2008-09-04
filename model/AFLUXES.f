#include "rundeck_opts.h"

      MODULE AFLUXES
!@sum  AFLUXES contains some ocean arrays defined on the atmospheric grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, ONLY : IMA=>IM, JMA=>JM

! Arrays  MOA,UOA,VOA,G0M,S0M,OGEOZA,OGEOZ_SVA,TRMOA are ocean quantities
!   after interpolation from ocean to atmospheric grid  
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: MOA, G0MA, S0MA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OGEOZ_A,OGEOZ_SV_A
     *      , UOA1,VOA1
#ifdef TRACERS_OCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRMOA
#endif

! Global arrays needed for interpolation from ocean to atmospheric grid  
      REAL*8, DIMENSION(IMA,JMA,2) :: MOA_glob, G0MA_glob, S0MA_glob
      REAL*8, DIMENSION(IMA,JMA) :: UOA1_glob, VOA1_glob
      REAL*8, DIMENSION(IMA,JMA) :: OGEOZ_A_glob, OGEOZ_SV_A_glob 
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(IMA,JMA,1,NTM) :: TRMOA_glob
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

      ALLOCATE(   MOA(I_0H:I_1H,J_0H:J_1H,2), STAT = IER)
      ALLOCATE( G0MA (I_0H:I_1H,J_0H:J_1H,2), STAT = IER)
      ALLOCATE( S0MA (I_0H:I_1H,J_0H:J_1H,2), STAT = IER)

      ALLOCATE( OGEOZ_A    (I_0H:I_1H, J_0H:J_1H),
     &          OGEOZ_SV_A (I_0H:I_1H, J_0H:J_1H),
     &          UOA1       (I_0H:I_1H, J_0H:J_1H), 
     &          VOA1       (I_0H:I_1H, J_0H:J_1H), 
     &   STAT=IER )

#ifdef TRACERS_OCEAN
      ALLOCATE(TRMOA(I_0H:I_1H,J_0H:J_1H,1,NTM), STAT = IER)
#endif

      END SUBROUTINE ALLOC_AFLUXES

