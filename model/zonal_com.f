c****
      module ZONAL_COM

      USE MODEL_COM, only : im,jm 
      integer  :: ncells
c      integer,parameter:: im=144,jm=90,ic=48,jc=48,ndomains=6
      integer,parameter:: ic=48,jc=48,ndomains=6

      integer, pointer, dimension(:,:) :: az11,az21
      integer, pointer, dimension(:,:) :: az31,az41,az12,az22
      real*8, pointer, dimension(:,:) :: az51,az32
      integer, dimension(6) :: keymax

      interface init_xgrid_unrolled
      subroutine init_xgrid_unrolled(im,jm,ic,jc,ncells,ndomains,
     *     az11,az21,az31,az41,az51,keymax)
      integer  :: ncells
      integer:: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: az11,az21
      integer, pointer, dimension(:,:) :: az31,az41
      real*8, pointer, dimension(:,:) :: az51
      integer, dimension(6) :: keymax 
      end subroutine init_xgrid_unrolled
      end interface 


      interface init_xgrid_loop
      subroutine init_xgrid_loop(im,jm,ic,jc,ncells,ndomains,
     *     az12,az22,az32,keymax)
      integer  :: ncells
      integer:: im,jm,ic,jc,ndomains
      integer, pointer, dimension(:,:) :: az12,az22
      real*8, pointer, dimension(:,:) :: az32
      integer, dimension(6) :: keymax 
      end subroutine init_xgrid_loop
      end interface 

      end module ZONAL_COM
c**** 
