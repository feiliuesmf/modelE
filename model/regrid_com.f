      module REGRID_COM
!@sum  the regrid module contains variables and procedures 
!@+    to remap data between two grids.
!@+ 
!@+    We define the x_grid container which stores information and data about the 
!@+    exchange grid. This data is common to the source and target grids
!@+    We also define the x_2grids container which stores an exchange grid
!@+    plus info about the direction of the remapping and about the resolution 
!@+    of the source and target grids
!@+ 
!@auth Denis Gueyffier (dgueyffier@nasa.giss.gov)

#ifndef DECOMP
      use DOMAIN_DECOMP, only: grid,dist_grid,AM_I_ROOT,SUMXPE
#endif
#ifdef CUBE_GRID
ccc Remove lines below when cs is instanciated inside DOMAIN_DECOMP.f
      use mpp_mod,only : mpp_sum
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain, mytile=>tile
ccc end remove
#else
      integer ::  is,ie,js,je,isd,ied,jsd,jed,gid,mytile
#endif

      integer,parameter:: ilonm=288,jlatm=180,ntiles=6,ic=48,jc=48
c      integer,parameter:: ilonm=144,jlatm=90,ntiles=6,ic=48,jc=48
c      integer,parameter:: ilonm=72,jlatm=46,ntiles=6,ic=48,jc=48
      integer,parameter:: dom_per_tile=4   ! #domains per cube face

ccc   routine members
      public :: init_regrid
      public :: init_csgrid_debug

ccc   variable members
      private :: xgrid
      public :: regrid

      type x_grid   ! stores x-grid information common to source and target grids
      integer  :: ncells
      real*8,  allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:) :: tile
      integer, allocatable, dimension(:,:) :: ijcub
      integer, allocatable, dimension(:,:) :: ijlatlon
      
      integer, allocatable, dimension(:) :: index_from_key
      integer, allocatable, dimension(:) :: jlat_from_key
      real*8,  allocatable, dimension(:) :: xarea_from_key
      integer :: maxkey      
      end type x_grid


      type x_2grids   ! stores x-grid info plus source and target grid info 
      type (x_grid) :: xgrid
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      end type x_2grids

      end module REGRID_COM
