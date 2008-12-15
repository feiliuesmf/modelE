#include "rundeck_opts.h"
      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth T. Clune
!@ver  1.0 (CS grid version)
!@cont GEOM_CS


      IMPLICIT NONE
      SAVE

!@var AREAG global integral of area (m^2)
      REAL*8 :: AREAG
!@var  LAT latitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lat2d(:,:)
      REAL*8, ALLOCATABLE :: lat2d_dg(:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sinlat2d, coslat2d
!@var  LON longitude of mid points of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lon2d(:,:)
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DXYP,BYDXYP


      SUBROUTINE GEOM_CS
!@sum  GEOM_CS Calculate spherical geometry for CS grid
!@auth T. Clune
!@ver  1.0 (CS grid version)
      USE CONSTANT, only : RADIUS,TWOPI,radian
      use fv_grid_tools_mod, only: CSgrid => grid, area
      use DOMAIN_DECOMP, only: grid, get
      implicit none

      integer :: i0h, i1h
      integer :: j0h, j1h
      
      call get(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h, I_STRT_HALO=i0h, I_STOP_HALO=i1h)

      allocate(lat2d(i0h:i1h, j0h:j1h)
      allocate(lat2d_deg(i0h:i1h, j0h:j1h)
      allocate(lon2d(i0h:i1h, j0h:j1h)
      allocate(sinlat2d(i0h:i1h, j0h:j1h)
      allocate(coslon2d(i0h:i1h, j0h:j1h)

      allocate(DXYP(i0h:i1h, j0h:j1h)
      allocate(BYDXYP(i0h:i1h, j0h:j1h)

      AREAG = 2 * TWOPI * RADIUS*RADIUS

      ! From FV CS grid
      lon2d = CSgrid(:,:,1)
      lat2d = CSgrid(:,:,2)
      lat2d_dg = lat2d/radian
      sinlat2d = sin(lat2d)
      coslat2d = cos(lat2d)

      DXYP = area(:,:)
      BYDXYP = 1/DXYP

      END SUBROUTINE GEOM_CS

      END MODULE GEOM


