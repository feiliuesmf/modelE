#include "rundeck_opts.h"

c
c Definitions of non-dimensional cube coordinates x,y:
c
c Except where indicated otherwise, all routines in this file assume
c that input/output x,y are defined to vary linearly with grid indices
c i,j (im=jm):
c x_edge(i) = -1+2*(i-1)/im    x_midpoint(i) = -1+2*(i-1/2)/im
c y_edge(j) = -1+2*(j-1)/jm    y_midpoint(j) = -1+2*(j-1/2)/jm
c
c Grid spacing is taken to give equal lengths along cube edges, i.e.
c y = lat/g when x=1; g is the latitude of cube corners.
c g = acos(1/3)/2 = asin(1/sqrt(3)) = atan(1/sqrt(2))
c
c For convenience, many routines internally define an alternative
c space denoted here as "capital" X,Y:
c X = sqrt(2)*tan(g*x)  Y = sqrt(2)*tan(g*y)
c X = tan(rotated_lon)  Y = tan(rotated_lat)/cos(rotated_lon)
c In addition to allowing many formulas to be written more compactly,
c this transformation has the following properties:
c - great circles are straight lines in X-Y space
c - latitude circles on polar tiles are X*X + Y*Y = const
c - latitude circles on equatorial tiles are Y*Y = const*(1+X*X)
c

      MODULE GEOM
!@sum  GEOM contains geometric variables and arrays
!@auth M. Kelley
!@ver  1.0 (Gnomonic CS grid version)
!@cont GEOM_CS

      USE MODEL_COM, only : IM,JM
      USE CONSTANT, only : radius,twopi
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

c     shift grid 10 degrees West to avoid corner over Japan
      real*8, parameter :: shiftwest = twopi/36.

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
      real*8 :: dloni,dlonj,dlati,dlatj,bydet

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

c
c calculate corner lons/lats, and areas integrated from the center
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
        lon2d(i,j) = lon2d(i,j)-shiftwest
        lat2d_dg(i,j) = lat2d(i,j)/radian
        lon2d_dg(i,j) = lon2d(i,j)/radian
        if(lon2d_dg(i,j) .lt. -180.) lon2d_dg(i,j)=lon2d_dg(i,j)+360.
        sinlat2d(i,j) = sin(lat2d(i,j))
        coslat2d(i,j) = cos(lat2d(i,j))
        lon2d(i,j) = lon2d(i,j) + pi ! IDL has a value of zero
        if(lon2d(i,j) .lt. 0.) lon2d(i,j)= lon2d(i,j) + twopi
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
c converts x,y,tile to lon,lat (radians)
c This routine places the center of tile 1 at the IDL.
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
c this version is for the gnomonic grid.  ll are in degrees east.
c if lon,lat do not lie on this tile, i,j are set to -99
      use model_com, only : im,jm
      use constant, only : radian,pi,twopi
      use domain_decomp_atm, only : grid
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      real*8 :: x,y,lonshift
      integer :: tile
      lonshift = ll(1)*radian+shiftwest
      if(lonshift.gt.pi) lonshift=lonshift-twopi
      call ll2csxy(lonshift,ll(2)*radian,x,y,tile)
      if(tile.eq.grid%tile) then
        ij(1) = min(1+int(.5*(x+1d0)*im),im)
        ij(2) = min(1+int(.5*(y+1d0)*jm),jm)
      else
        ij = -99
      endif
      return
      end subroutine lonlat_to_ij

      subroutine ll2csxy(lon,lat,x,y,tile)
c converts lon,lat (radians) to x,y,tile.
c This routine places the center of face 1 at the IDL.
      USE CONSTANT, only : PI
      implicit none
      real*8, parameter :: byg=1.62474893308877d0 ! byg=2/acos(1/3)
      real*8 :: lon,lat ! input
      real*8 :: x,y ! output
      integer :: tile ! output
      real*8 :: modlon,coslon,tanlat,tmpx,tmpy,bytanfac
      modlon = lon
      do while(modlon.lt.-pi/4d0)
        modlon = modlon + pi/2d0
      enddo
      do while(modlon.gt.+pi/4d0)
        modlon = modlon - pi/2d0
      enddo
      coslon = cos(modlon)
      tanlat = tan(lat)
      y = byg*atan(tanlat/(coslon*sqrt(2d0)))
      if(abs(y).le.1d0) then
c equatorial face
        x = byg*atan(tan(modlon)/sqrt(2d0))
c determine which face (1,2,4,5) we are on.  integer arithmetic
        tile = 1 + int((lon+1.25*pi)*2./pi)
        tile = tile +(tile/3) -5*(tile/5)
        if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
          tmpx = x
          tmpy = y
          x = -tmpy
          y = +tmpx
        endif
      else
c polar face
        bytanfac = sqrt(.5d0)/abs(tanlat)
        x = byg*atan(bytanfac*cos(lon))
        y = byg*atan(bytanfac*sin(lon))
        if(lat.gt.0.) then
          tile = 3
        else ! tile 6 x,y = tile 3 x,y flipped around the axis x+y=0
          tile = 6
          tmpx = x
          tmpy = y
          x = -tmpy
          y = -tmpx
        endif
      endif
      return
      end subroutine ll2csxy

      subroutine ll2csxy_vec(x_in,y_in,tile,ull,vll,uxy,vxy)
c
c Given latlon-oriented vector components ull,vll located at a
c position x_in,y_in on tile, return "vector" components uxy,vxy
c defining the directed great circle parallel to ull,vll.
c This great circle is a straight line in the X-Y space of this tile:
c         vxy*dX = uxy*dY
c See the beginning of this file for "capital" X,Y definitions.
c This routine could be made independent of the grid variant if
c position were specified as lon,lat or "capital" X,Y.
c
      implicit none
      real*8 :: x_in,y_in,ull,vll ! input
      integer :: tile ! input
      real*8 :: uxy,vxy ! output
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: x,y,rtxy,mag,bymag,utmp,vtmp
      x = x_in
      y = y_in
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        x = -y_in
        y = +x_in
      endif
      x = tan(g*x)*sqrt(2d0)
      y = tan(g*y)*sqrt(2d0)
      rtxy = sqrt(1d0+x*x+y*y)
      if(tile.eq.3 .or. tile.eq.6) then
        uxy = -y*ull - x*rtxy*vll
        vxy = +x*ull - y*rtxy*vll
        if(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
          uxy = -uxy
          vxy = -vxy
        endif
      else
        uxy = (1.+x*x)*ull
        vxy = x*y*ull + rtxy*vll
        if(tile.eq.4 .or. tile.eq.5) then
          utmp = uxy; vtmp = vxy
          uxy = -vtmp
          vxy = +utmp
        endif
      endif
      mag = sqrt(uxy**2 + vxy**2)
      if (mag .gt. 1.d-8) then
         bymag=1./mag
         uxy = uxy*bymag
         vxy = vxy*bymag
      endif
      return
      end subroutine ll2csxy_vec

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
      real*8 :: gx,gy,tangx,tangy,r2,tmpx,tmpy
      gx = g*x
      gy = g*y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        tmpx = gx
        tmpy = gy
        gx = +tmpy
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
          e1 = -e1
          e2 = -e2
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

      subroutine trigint(tile,x_in,y_in,a,cltcln,cltsln,slt)
c computes the integrals from 0 to x_in, 0 to y_in of
c area:                   a
c cos(lat)*cos(lon)*area: cltcln
c cos(lat)*sin(lon)*area: cltsln
c sin(lat)*area:          slt
c 
      implicit none
      integer :: tile
      real*8 :: x_in,y_in
      real*8 :: a,cltcln,cltsln,slt
      real*8 :: x,y,tmpx,tmpy,rtx,rty,atanx,atany
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      x = tan(g*x_in)*sqrt(2d0)
      y = tan(g*y_in)*sqrt(2d0)
      a = atan(x*y/sqrt(1.+x*x+y*y))
      tmpx = x; tmpy = y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        x = +tmpy
        y = -tmpx
      elseif(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
        x = -tmpy
        y = -tmpx
      endif
      rtx = sqrt(1.+x*x)
      rty = sqrt(1.+y*y)
      atanx = -.5*atan(x/rty)/rty
      atany = -.5*atan(y/rtx)/rtx
      if(tile.eq.3 .or. tile.eq.6) then ! polar tile
        cltsln = atanx
        cltcln = atany
        slt    = -(x*atany + y*atanx)
      else
        slt    = atanx
c longitude offsets; rotation on tiles 4/5
        if(tile.eq.1 .or. tile.eq.4) then
          cltsln = -atany
          cltcln = x*atany + y*atanx
        else
          cltsln = x*atany + y*atanx
          cltcln = atany
        endif
      endif
      if(tile.ge.4) slt = -slt
      return
      end subroutine trigint

      END MODULE GEOM
