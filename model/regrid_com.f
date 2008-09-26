      module REGRID_COM
      
      integer  :: ncells
      integer,parameter:: im=288,jm=180,ntiles=6,icc=48,jcc=48
      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:,:) :: ijcub,ijlatlon
      integer, allocatable, dimension(:) :: tile
      
      interface init_regrid
      subroutine init_regrid(pelist)
      integer, dimension(:), allocatable :: pelist
      end subroutine init_regrid
      end interface 
      
      end module REGRID_COM

