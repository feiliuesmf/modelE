      module REGRID_COM
!@sum Contains conservative regridding, interpolation and zonal mean routines
!@+   to exchange data between grids (latlon and cubed sphere) 
!@+ 
!@+    We define the x_grid container which stores information and data about the 
!@+    exchange grid. This data is common to the source and target grids
!@+    see http://www.gfdl.noaa.gov/~vb/gridstd/gridstdse2.html#x4-160002.7
!@     for precise definition of the concept of exchange grid
!@+    We also define the x_2grids type which contains an x_grid (exchange grid)
!@+    plus info about the direction of the remapping and about the resolution 
!@+    of the source and target grids
!@+ 
!@auth Denis Gueyffier


      use DOMAIN_DECOMP_1D, only: AM_I_ROOT,SUMXPE
      use dd2d_utils, only : dist_grid

      integer, parameter :: nrecmax=200 ! max #records in input files

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
      integer :: is,ie,js,je,isd,ied,jsd,jed ! domain decomp. for source grid
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

c***  create specializations of above types, may be moved to a different layer later
      type (x_2grids) :: xO2A    ! ocean to atmosphere
      type (x_2grids) :: xA2O    ! atmosphere to ocean
      type (x_2grids) :: xI2A    ! ice grid to atmosphere
      type (x_2grids) :: xA2I    ! atmosphere to ice grid

c*** interface for functions
      interface
          function vminus(p,q)
          real*8 :: vminus(2),p(2),q(2)
          end function vminus
      end interface

      interface 
          function cross_product(p,q)
          real*8 :: cross_product,p(2),q(2)
          end function cross_product
      end interface

      end module REGRID_COM
