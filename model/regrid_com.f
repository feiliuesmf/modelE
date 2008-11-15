      module REGRID_COM
!@sum  REGRID  contains variables and procedures to regrid between cubed-sphere and latlon grids 
!@auth Denis Gueyffier (dgueyffier@nasa.giss.gov)

#ifndef DECOMP
      use DOMAIN_DECOMP, only: grid,AM_I_ROOT,SUMXPE
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

      integer  :: ncells
      integer,parameter:: ilonm=288,jlatm=180,ntiles=6,ic=48,jc=48
c      integer,parameter:: ilonm=144,jlatm=90,ntiles=6,ic=48,jc=48
c      integer,parameter:: ilonm=72,jlatm=46,ntiles=6,ic=48,jc=48
      integer,parameter:: dom_per_tile=4   ! #domains per cube face

      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:,:) :: ijcub,ijlatlon
      integer, allocatable, dimension(:) :: tile
      
      integer, allocatable, dimension(:) :: az12,az22
      real*8, allocatable, dimension(:) :: az32
      integer :: maxkey


      interface init_regrid
      subroutine init_regrid
      end subroutine
      end interface 

      interface init_csgrid_debug
      subroutine init_csgrid_debug
      end subroutine
      end interface

      interface readt_regrid_parallel
      subroutine readt_regrid_parallel(iunit,name,nskip,tcubloc,ipos)
      integer, intent(in) :: iunit
      character*16, intent(in) :: name
      integer, intent(in) :: nskip
      real*8 :: tcubloc(:,:) !output regridded and scattered data 
      integer, intent(in) :: ipos
      end subroutine
      end interface

      interface test_globalsum_exact
      subroutine test_globalsum_exact
      end subroutine
      end interface

      interface test_regrid_exact
      subroutine test_regrid_exact
      end subroutine
      end interface

      end module REGRID_COM
