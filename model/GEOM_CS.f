#include "rundeck_opts.h"
      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth T. Clune
!@ver  1.0 (CS grid version)
!@cont GEOM_CS

      USE MODEL_COM, only : IM,JM,LM,FIM,BYIM
      USE CONSTANT, only : radius
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
!@var ddx_ci, ddx_cj,ddy_ci ddy_cj 
      REAL*8, ALLOCATABLE :: ddx_ci(:,:),ddx_cj(:,:)
      REAL*8, ALLOCATABLE :: ddy_ci(:,:),ddy_cj(:,:)
      real*8 :: dloni,dlonj,dlati,dlatj,bydet

!@var quantities below are temporary fix - quantity may not be used on CS
      REAL*8, parameter :: dlon =1.0 
      REAL*8, parameter :: dlat_dg = 1.0
      REAL*8, dimension(JM) :: lat
      REAL*8, dimension(IM) :: lon
      REAL*8, dimension(JM,2) :: lat_dg
      REAL*8, dimension(IM,2) :: lon_dg
      REAL*8, DIMENSION(JM) :: COSP,COSV
      REAL*8, DIMENSION(JM) :: SINP
      REAL*8 :: SINV
      REAL*8, DIMENSION(JM,2,2) :: WTJ
      INTEGER, parameter :: J1U = 2
      INTEGER, parameter, dimension(2,2,2) :: JRANGE_HEMI = reshape(
     *  (/1,JM/2,  1+JM/2,JM,  J1U,J1U-1+JM/2, J1U-1+JM/2,JM+J1U-2/),
     *  (/2,2,2/))
      REAL*8, DIMENSION(JM) :: DXYV,BYDXYV,DYV
      REAL*8, DIMENSION(IM) :: SINIV,COSIV,SINIP,COSIP
!@var


!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: axyp, byaxyp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: dxyp, bydxyp

      integer, allocatable, dimension(:) :: imaxj

!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM) :: KMAXJ

      integer, parameter :: jg_u=2, jg_ke=2

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
      use DOMAIN_DECOMP_ATM, only: grid, get
      implicit none

      integer :: i0h, i1h, i0, i1
      integer :: j0h, j1h, j0, j1, j0s, j1s
      integer :: i,j
      call get(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h,
     &     J_STRT=j0, J_STOP=j1,
     &     I_STRT_HALO=i0h, I_STOP_HALO=i1h, 
     &     I_STRT=i0,I_STOP=i1, 
     &     J_STRT_SKP=j0s, J_STOP_SKP=j1s)

      write(*,*) "geom  i0h,i1h,i1,j0h,j1h,j0,j1=",
     &     i0h,i1h,i1,j0h,j1h,j0,j1

      allocate(lat2d(i0h:i1h, j0h:j1h))
      allocate(lon2d(i0h:i1h, j0h:j1h))

      allocate(lat2d_dg(i0h:i1h, j0h:j1h))
      allocate(lon2d_dg(i0h:i1h, j0h:j1h))

      allocate(sinlat2d(i0h:i1h, j0h:j1h))
      allocate(coslat2d(i0h:i1h, j0h:j1h))

      allocate(sinlon2d(i0h:i1h, j0h:j1h))
      allocate(coslon2d(i0h:i1h, j0h:j1h))

      allocate(lat2d_corner(i0h:i1h+1, j0h:j1h+1))
      allocate(lon2d_corner(i0h:i1h+1, j0h:j1h+1))

      allocate(
     &    ddx_ci(i0h:i1h,j0h:j1h)
     &    ,ddx_cj(i0h:i1h,j0h:j1h)
     &    ,ddy_ci(i0h:i1h,j0h:j1h)
     &    ,ddy_cj(i0h:i1h,j0h:j1h))

      allocate(axyp(i0h:i1h, j0h:j1h))
      allocate(byaxyp(i0h:i1h, j0h:j1h))
      allocate(dxyp(j0h:j1h))
      allocate(bydxyp(j0h:j1h))

      allocate(imaxj(j0:j1))

      AREAG = 2 * TWOPI * RADIUS*RADIUS

      ! From FV CS grid
c      lon2d = agrid(:,:,1)
c      lat2d = agrid(:,:,2)
       lon2d(:,:) = 1.0  ! until problem with Initialize(fv) is fixed
       lat2d(:,:) = 1.0  ! idem

c      lon2d_corner = corner_grid(:,:,1)
c      lat2d_corner = corner_grid(:,:,2)

      lat2d_dg = lat2d/radian
      lon2d_dg = lon2d/radian
      sinlat2d = sin(lat2d)
      coslat2d = cos(lat2d)

c      AXYP = area(:,:)
c      BYAXYP = 1/AXYP
      axyp(:,:) =1.0  ! until problem with Initialize(fv) is fixed
      byaxyp(:,:) =1.0

      call set_j_budg   !called after lat2d_dg is initialized
      call set_wtbudg !sets area weights

      imaxj(:)=i1

      do j=j0s,j1s
      do i=i0,i1
      dloni = lon2d(i+1,j)-lon2d(i-1,j)
      dlati = lat2d(i+1,j)-lat2d(i-1,j)
      dlonj = lon2d(i,j+1)-lon2d(i,j-1)
      dlatj = lat2d(i,j+1)-lat2d(i,j-1)
      if(dloni.lt.0.) dloni=dloni+twopi
      if(dlonj.lt.0.) dlonj=dlonj+twopi
      bydet = 1d0/(dloni*dlatj-dlonj*dlati)
      ddx_ci(i,j) =  dlatj*bydet/(radius*coslat2d(i,j))
      ddx_cj(i,j) = -dlonj*bydet/(radius*coslat2d(i,j))
      ddy_ci(i,j) = -dlati*bydet/radius
      ddy_cj(i,j) =  dloni*bydet/radius
      enddo
      enddo

c     some values below may be inexact or deprecated - temporary fix for compilation
      kmaxj(:)= 2
      DYV(:) = 1.0
      lat(:) = 1.0 
      lon(:) = 1.0
      lat_dg(:,:) = 1.0
      lon_dg(:,:) = 1.0
      cosp(:) = 1.0
      cosv(:) = 1.0
      sinp(:) = 1.0
      sinv = 1.0
      wtj(:,:,:) = 1.0
      cosip(:)= 1.0
      sinip(:) = 1.0
c     end inexact

      write(*,*) "WARNING : kmaxj, DYP, lat, lon, lat_dg, 
     &     lon_dg, ... may have inexact values"

      END SUBROUTINE GEOM_CS

      END MODULE GEOM


