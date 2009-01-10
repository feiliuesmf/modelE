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
!@var  lat2d latitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lat2d(:,:)
!@var  lat2d_dg latitude of mid point of primary grid box (degrees)
      REAL*8, ALLOCATABLE :: lat2d_dg(:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sinlat2d, coslat2d
!@var  lon2d longitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lon2d(:,:)
!@var  lon2d_dg longitude of mid point of primary grid box (degrees)
      REAL*8, ALLOCATABLE :: lon2d_dg(:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sinlon2d, coslon2d

!@var  lat2d_corner latitude of corner point (radians)
      REAL*8, ALLOCATABLE :: lat2d_corner(:,:)
!@var  lon2d_corner longitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE :: lon2d_corner(:,:)

!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: axyp, byaxyp, dxyp, bydxyp,
     &     dxp, dyp

      integer, allocatable, dimension(:) :: imaxj

!@var J_BUDG a mapping array that takes every grid point to the
!@+   zonal mean budget array
      integer, allocatable, dimension(:,:) :: J_BUDG
!@var j_0b, j_1b are the min/max zonal budget latitudes for this processor
      integer :: j_0b, j_1b

      CONTAINS

      SUBROUTINE GEOM_CS
!@sum  GEOM_CS Calculate spherical geometry for CS grid
!@auth T. Clune
!@ver  1.0 (CS grid version)
      USE CONSTANT, only : RADIUS,TWOPI,radian
      use fv_grid_tools_mod, only: agrid, corner_grid => grid, area
      use DOMAIN_DECOMP, only: grid, get
      implicit none

      integer :: i0h, i1h, i1
      integer :: j0h, j1h, j0, j1
      
      call get(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h,
     &     J_STRT=j0, J_STOP=j1,
     &     I_STRT_HALO=i0h, I_STOP_HALO=i1h, 
     &     I_STOP=i1)

      allocate(lat2d(i0h:i1h, j0h:j1h))
      allocate(lon2d(i0h:i1h, j0h:j1h))

      allocate(lat2d_dg(i0h:i1h, j0h:j1h))
      allocate(lon2d_dg(i0h:i1h, j0h:j1h))

      allocate(sinlat2d(i0h:i1h, j0h:j1h))
      allocate(coslat2d(i0h:i1h, j0h:j1h))

      allocate(lat2d_corner(i0h:i1h+1, j0h:j1h+1))
      allocate(lon2d_corner(i0h:i1h+1, j0h:j1h+1))

      allocate(axyp(i0h:i1h, j0h:j1h))
      allocate(byaxyp(i0h:i1h, j0h:j1h))
      allocate(imaxj(j0:j1))

      AREAG = 2 * TWOPI * RADIUS*RADIUS

      ! From FV CS grid
      lon2d = agrid(:,:,1)
      lat2d = agrid(:,:,2)

      lon2d_corner = corner_grid(:,:,1)
      lat2d_corner = corner_grid(:,:,2)

      lat2d_dg = lat2d/radian
      lon2d_dg = lon2d/radian
      sinlat2d = sin(lat2d)
      coslat2d = cos(lat2d)

      DXYP = area(:,:)
      BYDXYP = 1/DXYP

      call set_j_budg   !called after lat2d_dg is initialized
      call set_wtbudg !sets area weights

      imaxj(:)=i1

c     values below are wrong
      DXP(:,:)  = 1.0 ! this needs to be calculated 
      DYP(:,:)  = 1.0 ! this needs to be calculated

      END SUBROUTINE GEOM_CS

      END MODULE GEOM


