      module REGRID_COM
!@sum  REGRID  contains variables and procedures to regrid between cubed-sphere and latlon grids 
!@auth Denis Gueyffier (dgueyffier@nasa.giss.gov)

      use DOMAIN_DECOMP, only: grid,AM_I_ROOT,SUMXPE
      integer ::  is,ie,js,je,isd,ied,jsd,jed,gid,mytile
      integer  :: ncells
      integer,parameter:: ilonm=288,jlatm=180,ntiles=6,ic=48,jc=48
c      integer,parameter:: ilonm=144,jlatm=90,ntiles=6,ic=48,jc=48
      integer,parameter:: dom_per_tile=4   ! #domains per cube face

      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:,:) :: ijcub,ijlatlon
      integer, allocatable, dimension(:) :: tile
      
      integer, allocatable, dimension(:) :: az12,az22
      real*8, allocatable, dimension(:) :: az32
      integer :: maxkey

      interface init_regrid
      subroutine init_regrid()
      end subroutine init_regrid
      end interface 

      interface init_xgrid_zonal
      subroutine init_xgrid_zonal()
      end subroutine init_xgrid_zonal
      end interface 
      
      end module REGRID_COM
