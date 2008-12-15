#include "rundeck_opts.h"
      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth T. Clune
!@ver  1.0 (CS grid version)
!@cont GEOM_CS

      use fv_grid_tools_mod, only: 


      IMPLICIT NONE
      SAVE

!@var AREAG global integral of area (m^2)
      REAL*8 :: AREAG
!@var  LAT latitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: LAT(:,:)
!@var  LON longitude of mid points of primary grid box (radians)
      REAL*8, ALLOCATABLE :: LON(:,:)
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DXYP,BYDXYP

!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!@+     (+inverse)
      REAL*8, DIMENSION(JM) :: DXP,DYP,BYDXP,BYDYP

      SUBROUTINE GEOM_CS
!@sum  GEOM_CS Calculate spherical geometry for CS grid
!@auth T. Clune
!@ver  1.0 (CS grid version)

      use fv_grid_tools_mod, only: CSgrid => grid
      use DOMAIN_DECOMP, only: grid, get
      implicit none

      integer :: i0h, i1h
      integer :: j0h, j1h
      
      call get(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h, I_STRT_HALO=i0h, I_STOP_HALO=i1h)

      allocate(LAT(i0h:i1h, j0h:j1h)
      allocate(LON(i0h:i1h, j0h:j1h)
      allocate(DXYP(i0h:i1h, j0h:j1h)
      allocate(BYDXYP(i0h:i1h, j0h:j1h)

      AREAG = 2 * TWOPI * RADIUS*RADIUS

      ! From FV CS grid
      LON = CSgrid(:,:,1)
      LAT = CSgrid(:,:,2)
      DXYP = area(:,:)
      BYDXYP = 1/DXYP

      END SUBROUTINE GEOM_CS

      END MODULE GEOM


