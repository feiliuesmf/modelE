      module REGRID_COM
!@sum  the regrid module contains variables and procedures 
!@+    to remap data between two grids.
!@+ 
!@+    We define the x_grid container which stores information and data about the 
!@+    exchange grid. This data is common to the source and target grids
!@+    We also define the x_2grids type which contains an x_grid (exchange grid)
!@+    plus info about the direction of the remapping and about the resolution 
!@+    of the source and target grids
!@+ 
!@auth Denis Gueyffier (dgueyffier@nasa.giss.gov)

#ifndef DECOMP
      use DOMAIN_DECOMP, only: grid,dist_grid,AM_I_ROOT,SUMXPE
#endif
#ifdef CUBE_GRID
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain, mytile=>tile
ccc end remove
#else
      integer ::  is,ie,js,je,isd,ied,jsd,jed,gid,mytile
#endif

      integer,parameter:: dom_per_tile=4   ! #domains per cube face
      integer, parameter :: nrecmax=100 ! max #records in input files

ccc   routine members
c      public :: init_regrid
c      public :: init_csgrid_debug

ccc   variable members
      private :: xgrid
      public :: x_2grids

ccc   derived types
      type x_grid   ! stores x-grid information common to source and target grids
                    ! x_grid is distributed and contains data local to each domain
      integer, allocatable, dimension(:) :: index_key
      integer, allocatable, dimension(:) :: icub_key
      integer, allocatable, dimension(:) :: jcub_key
      integer, allocatable, dimension(:) :: ilon_key
      integer, allocatable, dimension(:) :: jlat_key
      real*8, allocatable, dimension(:) :: xarea_key
      integer, allocatable, dimension(:) :: itile_key
      integer :: maxkey      
      integer :: ncells
      end type x_grid

      type x_gridroot   ! stores x-grid information common to source and target grids
                        ! x_grid data is global and stays on root proc
      integer  :: ncells
      integer, allocatable, dimension(:,:) :: ijcub
      integer, allocatable, dimension(:,:) :: ijlatlon
      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:) :: tile
      integer :: maxkey      
      end type x_gridroot

      type x_2grids   ! stores x-grid info plus source and target grid info 
      type (x_grid) :: xgrid
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      end type x_2grids

      type x_2gridsroot   ! stores x-grid info plus source and target grid info 
      type (x_gridroot) :: xgridroot
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      end type x_2gridsroot

      end module REGRID_COM
