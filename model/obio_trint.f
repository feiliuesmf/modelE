#include "rundeck_opts.h"

      subroutine obio_trint

#ifdef OBIO_ON_GARYocean
      use oceanres,  only: idm=>imo,jdm=>jmo,kdm=>lmo
      use oceanr_dim, only : ogrid
      use obio_com, only: tracer       !use global array
      use ocn_tracer_com, only : ntrcr=>ntm
      use model_com, only : nstep=>itime
#else
      use hycom_dim_glob, only : ntrcr,idm,jdm,kdm
      use hycom_dim, only: ogrid, ifu, ilu, isu
      use hycom_arrays, only : tracer,dpinit,scp2
      use hycom_scalars, only : nstep,huge
#endif
      use domain_decomp_1d, only: am_i_root, globalsum, get

      implicit none

      integer i,j,k,l,kn,nn
      integer ntr
      real sumo

      integer :: iTracer
      integer :: j_0, j_1, j_0h, j_1h

      real*8, allocatable :: summ(:)
      real*8, allocatable :: partialTracerSums(:,:)

      call get(ogrid, j_strt = j_0, j_stop = j_1,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h)

      sumo = getTotalOceanVolume()

      ! Form partial sum for each latitude for each tracer
      allocate(partialTracerSums(j_0h:j_1h, ntrcr))
      do iTracer = 1, ntrcr
         partialTracerSums(:,iTracer) = 
     &        partialIntegration(tracer(:,:,:,iTracer))
      end do

      allocate(summ(ntrcr))
      call globalSum(ogrid, partialTracerSums, summ)

      if (am_i_root()) then
         do iTracer = 1, ntrcr
            write(*,*) 'total intgrl tracer:', iTracer, nstep, 
     &           summ(iTracer), summ(iTracer)/sumo
         end do
      end if

      deallocate(partialTracerSums, summ)

      contains

#ifdef OBIO_ON_GARYocean
      function partialIntegration(quantity)
      use geom, only : dxyp
      use oceanres, only: dzo
      real*8, intent(in) :: quantity(:,j_0h:,:)
      real*8 :: partialIntegration(j_0h:j_1h)
      
      real*8 :: gridCellVolume

      partialIntegration = 0
      do k = 1, kdm
         do j = j_0, j_1
            gridCellVolume = dzo(k) * dxyp(j)
            do i= 1, size(quantity, 1)
               partialIntegration(j) = partialIntegration(j) + 
     &              quantity(i,j,k) * gridCellVolume
            end do
         end do
      end do

      end function partialIntegration
#else
      function partialIntegration(quantity)
      real*8, intent(in) :: quantity(:,j_0h:,:)
      real*8 :: partialIntegration(j_0h:j_1h)
      
      integer :: i, j, l, k

      partialIntegration = 0
      do k = 1, kdm
         do j = j_0, j_1
            do l = 1, isu(j)
               do i= ifu(j,l), ilu(j,l)
                  partialIntegration(j) = partialIntegration(j) + 
     &                 quantity(i,j,k)) * dpinit(i,j,j)*scp2(i,j)
               end do
            end do
         end do
      end do

      end function partialIntegration
#endif

      real*8 function getTotalOceanVolume() result(oceanVolume)
      real*8, allocatable :: partialOceanVolumes(:)
      real*8, allocatable :: ones(:,:,:)

      allocate(partialOceanVolumes(j_0h:j_1h))
      allocate(ones(idm, J_0H:J_1H, kdm))
      
      ones = 1
      partialOceanVolumes(:) = partialIntegration(ones)
      call globalSum(ogrid, partialOceanVolumes, oceanVolume)

      deallocate(ones)

      end function getTotalOceanVolume

      end subroutine obio_trint
