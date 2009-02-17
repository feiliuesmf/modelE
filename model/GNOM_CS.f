#include "rundeck_opts.h"
      MODULE GEOM
!@sum  GEOM contains geometric variables and arrays
!@auth M. Kelley
!@ver  1.0 (Gnomonic CS grid version)
!@cont GEOM_CS

      USE MODEL_COM, only : IM,JM
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

!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: axyp, byaxyp

      integer, allocatable, dimension(:) :: imaxj

!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM) :: KMAXJ


!@var J_BUDG a mapping array that takes every grid point to the
!@+   zonal mean budget array
      integer, allocatable, dimension(:,:) :: J_BUDG
!@var j_0b, j_1b are the min/max zonal budget latitudes for this processor
      integer :: j_0b, j_1b


      CONTAINS

      SUBROUTINE GEOM_CS
!@sum  GEOM_CS Calculate geometry for CS grid
!@auth M. Kelley
!@ver  1.0 (Gnomonic CS grid version)
      USE CONSTANT, only : RADIUS,PI,TWOPI,radian
      use DOMAIN_DECOMP_ATM, only: grid, get, halo_update
      implicit none
      real*8 :: x,y
      integer :: i0h, i1h, i0, i1
      integer :: j0h, j1h, j0, j1
      integer :: i,j
      real*8 :: dloni,dlonj,dlati,dlatj,bydet,shiftwest

      real*8, dimension(grid%i_strt:grid%i_stop+1,
     &                  grid%j_strt:grid%j_stop+1) :: axyp_int
      call get(grid, J_STRT_HALO=j0h, J_STOP_HALO=j1h,
     &     J_STRT=j0, J_STOP=j1,
     &     I_STRT_HALO=i0h, I_STOP_HALO=i1h, 
     &     I_STRT=i0,I_STOP=i1)

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
     &     ddx_ci(i0h:i1h,j0h:j1h)
     &    ,ddx_cj(i0h:i1h,j0h:j1h)
     &    ,ddy_ci(i0h:i1h,j0h:j1h)
     &    ,ddy_cj(i0h:i1h,j0h:j1h))

      allocate(axyp(i0h:i1h, j0h:j1h))
      allocate(byaxyp(i0h:i1h, j0h:j1h))

      allocate(imaxj(j0:j1))

      AREAG = 2 * TWOPI * RADIUS*RADIUS

c     shift grid 10 degrees West to avoid corner over Japan
      shiftwest = twopi/36.
c
c calculate corner lons/lats and areas integrated from the center
c of this cube face
c
      do j=j0,j1+1
      do i=i0,i1+1
        x = -1d0 + 2d0*(i-1)/im
        y = -1d0 + 2d0*(j-1)/im
        call csxy2ll(x,y,grid%tile,lon2d_corner(i,j),lat2d_corner(i,j))
        lon2d_corner(i,j)=lon2d_corner(i,j)-shiftwest
        if ( lon2d_corner(i,j) .lt. -pi) lon2d_corner(i,j)=
     &       lon2d_corner(i,j) + twopi
        axyp_int(i,j) = aint(x,y)
      enddo
      enddo
c
c cell areas, lons/lats at cell centers
c
      do j=j0,j1
      do i=i0,i1
        axyp(i,j) = axyp_int(i+1,j+1)-axyp_int(i,j+1)
     &             -axyp_int(i+1,j  )+axyp_int(i,j  )
        axyp(i,j) = axyp(i,j)*radius*radius
        byaxyp(i,j) = 1d0/axyp(i,j)
        x = -1d0 + 2d0*(dble(i)-.5d0)/im
        y = -1d0 + 2d0*(dble(j)-.5d0)/im
        call csxy2ll(x,y,grid%tile,lon2d(i,j),lat2d(i,j))
        lon2d(i,j)=lon2d(i,j)-shiftwest
        if ( lon2d(i,j) .lt. -pi) lon2d(i,j)= lon2d(i,j) + twopi
        lat2d_dg(i,j) = lat2d(i,j)/radian
        lon2d_dg(i,j) = lon2d(i,j)/radian
        sinlat2d(i,j) = sin(lat2d(i,j))
        coslat2d(i,j) = cos(lat2d(i,j))
        lon2d(i,j) = lon2d(i,j) + pi ! IDL has a value of zero
      enddo
      enddo

c halo everything just in case
      call halo_update(grid,axyp)
      call halo_update(grid,byaxyp)
      call halo_update(grid,lon2d)
      call halo_update(grid,lat2d)
      call halo_update(grid,lon2d_dg)
      call halo_update(grid,lat2d_dg)
      call halo_update(grid,sinlat2d)
      call halo_update(grid,coslat2d)

      call set_j_budg   !called after lat2d_dg is initialized
      call set_wtbudg !sets area weights

      imaxj(:)=i1

      do j=j0,j1
      do i=i0,i1
        dloni = lon2d(i+1,j)-lon2d(i-1,j)
        dlati = lat2d(i+1,j)-lat2d(i-1,j)
        dlonj = lon2d(i,j+1)-lon2d(i,j-1)
        dlatj = lat2d(i,j+1)-lat2d(i,j-1)
c account for discontinuity of lon2d at the international date line
        if(abs(dloni).gt.pi) dloni=dloni-twopi*sign(1d0,dloni)
        if(abs(dlonj).gt.pi) dlonj=dlonj-twopi*sign(1d0,dlonj)
        bydet = 1d0/(dloni*dlatj-dlonj*dlati)
        ddx_ci(i,j) =  dlatj*bydet/(radius*coslat2d(i,j))
        ddx_cj(i,j) = -dlati*bydet/(radius*coslat2d(i,j))
        ddy_ci(i,j) = -dlonj*bydet/radius
        ddy_cj(i,j) =  dloni*bydet/radius
      enddo
      enddo

      kmaxj(:)= 2

      return
      END SUBROUTINE GEOM_CS

      subroutine csxy2ll(x,y,tile,lon,lat)
c This routine places the center of tile 1 at the IDL.
c Gnomonic (great circle) grid generation rules:
c 1. y is proportional to latitude when lon=pi/4 (x=const)
c 2. tan(lat) proportional to cos(lon) along a great circle y=const
c x and y are nondimensional cube coordinates ranging from -1 to +1.
c On a polar face, the transformation X=tan(g*x), Y=tan(g*y) makes
c latitude lines into X*X+Y*Y=const and longitude lines into
c Y/X=const, where g=.5*acos(1/3).
      USE CONSTANT, only : PI
      implicit none
      real*8 :: x,y ! input
      integer :: tile ! input
      real*8 :: lon,lat ! output
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: gx,gy,tmpx,tmpy,cosgx,tangx,tangy,coslon
      gx = g*x
      gy = g*y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        tmpx = gx
        tmpy = gy
        gx = +tmpy
        gy = -tmpx
      elseif(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
        tmpx = gx
        tmpy = gy
        gx = -tmpy
        gy = -tmpx
      endif
      if(tile.eq.3 .or. tile.eq.6) then
        tangx = tan(gx)
        tangy = tan(gy)
        lat = atan(1d0/sqrt(2d0*(tangx**2 +tangy**2)+1d-40))
        lon = atan2(tangy,tangx)
        if(tile.eq.6) lat = -lat
      else
        cosgx = cos(gx)
        coslon = cosgx/sqrt(2d0-cosgx*cosgx)
        lat = atan(coslon*sqrt(2d0)*tan(gy))
        lon = sign(acos(coslon),gx)
! add longitude offset to tiles 1,2,4,5. integer arithmetic.
        lon = lon + (mod(tile,3)-1)*pi/2. -pi*(1-(tile/3))
        if(lon.lt.-pi) lon=lon+2.*pi
      endif
      return
      end subroutine csxy2ll

      function aint(x,y)
c calculates the area integral from the center of a cube face
c to the point x,y.
c x and y are nondimensional cube coordinates ranging from -1 to +1.
      implicit none
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: x,y
      real*8 :: aint
      real*8 :: tangx,tangy
      tangx = tan(g*x)
      tangy = tan(g*y)
      aint = atan(2.*tangx*tangy/sqrt(1.+2.*(tangx*tangx+tangy*tangy)))
      return
      end function aint

      subroutine lonlat_to_ij(ll,ij)
c converts lon,lat=ll(1:2) into model i,j=ij(1:2)
c this is a dummy version to allow compilation of new river/lakes routine
c Will be fixed soon.
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      ij(1) = 1
      ij(2) = 1
      return
      end subroutine lonlat_to_ij
      END MODULE GEOM

      subroutine e1e2(x,y,tile,e1,e2)
c
c this routine computes latlon-oriented basis vectors e1,e2 that
c are parallel to cubed-sphere gridlines at a position x,y on tile.
c e1(1:2) is parallel to the constant-y great circles and is
c  proportional to (/ cos(lat)*dlon/dx, dlat/dx /)
c e2(1:2) is parallel to the constant-x great circles and is
c  proportional to (/ cos(lat)*dlon/dy, dlat/dy /)
c
      implicit none
      real*8 :: x,y ! input
      integer :: tile ! input
      real*8, dimension(2) :: e1,e2
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: gx,gy,tangx,tangy,r2,tmpx,tmpy,tmpe(2)
      gx = g*x
      gy = g*y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        tmpx = gx
        tmpy = gy
        gx = +tmpy
        gy = -tmpx
      elseif(tile.eq.6) then    ! tile 6 = tile 3 flipped around the axis x+y=0
        tmpx = gx
        tmpy = gy
        gx = -tmpy
        gy = -tmpx
      endif
      tangx = tan(gx)
      tangy = tan(gy)
      r2 = 1d0/sqrt(1d0+2d0*(tangx**2+tangy**2))
      if(tile.eq.3 .or. tile.eq.6) then
        e1(1) = -tangy
        e1(2) = -tangx*r2
        e2(1) = tangx
        e2(2) = -tangy*r2
        if(tile.eq.6) then
          tmpe(:) = -e1
          e1(:) = -e2(:)
          e2(:) = tmpe(:)
          e1(2) = -e1(2)
          e2(2) = -e2(2)
        endif
      else
        e1(1)=1d0
        e1(2)=-2d0*tangx*tangy*r2
        e2(1)=0d0
        e2(2)=1d0
        if(tile.eq.4 .or. tile.eq.5) then
          e2(1)=1d0
          e2(2)=e1(2)
          e1(1)=0d0
          e1(2)=-1d0
        endif
      endif
      e1=e1/sqrt(sum(e1*e1))
      e2=e2/sqrt(sum(e2*e2))
      return
      end subroutine e1e2
