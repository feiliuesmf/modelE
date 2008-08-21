c****
      module ZONAL_COM

      USE MODEL_COM, only : im,jm
      integer  :: ncells
c      integer,parameter:: im=144,jm=90,ic=48,jc=48,ndomains=6
      integer,parameter:: ic=48,jc=48,ndomains=6

      integer, pointer, dimension(:,:) :: azonal1,azonal2
      integer, pointer, dimension(:,:) :: azonal3,azonal4
      real*8, pointer, dimension(:,:) :: azonal5
      integer, dimension(6) :: keymax 

      interface init_xgrid
      subroutine init_xgrid(im,jm,ic,jc,ncells,ndomains,
     *     azonal1,azonal2,azonal3,azonal4,azonal5,keymax)
      integer  :: ncells
      integer:: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: azonal1,azonal2
      integer, pointer, dimension(:,:) :: azonal3,azonal4
      real*8, pointer, dimension(:,:) :: azonal5
      integer, dimension(6) :: keymax 
      end subroutine init_xgrid
      end interface 

      end module ZONAL_COM
c**** 
