      module REGRID_COM

      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed

      integer  :: ncells
c      integer,parameter:: im=288,jm=180,ntiles=6,ic=48,jc=48
      integer,parameter:: im=144,jm=90,ntiles=6,ic=48,jc=48
      integer,parameter:: dom_per_tile=4   ! #domains per cube face

      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:,:) :: ijcub,ijlatlon
      integer, allocatable, dimension(:) :: tile
      
      integer, allocatable, dimension(:) :: az11,az21
      integer, allocatable, dimension(:) :: az31,az41,az12,az22
      real*8, allocatable, dimension(:) :: az51,az32
      integer :: maxkey

      interface init_regrid
      subroutine init_regrid(pelist)
      integer, dimension(:), allocatable :: pelist
      end subroutine init_regrid
      end interface 

      interface init_regrid_rootpe
      subroutine init_regrid_rootpe
      end subroutine init_regrid_rootpe
      end interface 

      interface init_grid_temp
      subroutine init_grid_temp()
      end subroutine init_grid_temp
      end interface 

      interface init_xgrid_unrolled
      subroutine init_xgrid_unrolled()
      end subroutine init_xgrid_unrolled
      end interface 

      interface init_xgrid_loop
      subroutine init_xgrid_loop()
      end subroutine init_xgrid_loop
      end interface 
      
      end module REGRID_COM
