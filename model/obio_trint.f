#include "rundeck_opts.h"

      subroutine obio_trint

#ifdef OBIO_ON_GARYocean
      use oceanres,  only: idm=>imo,jdm=>jmo,kdm=>lmo
      use oceanr_dim, only : ogrid
      use obio_com, only: tracer => tracer_loc ! rename local
      use ocn_tracer_com, only : ntrcr=>ntm
      use model_com, only : nstep=>itime
#else
      use hycom_dim_glob, only : ntrcr,idm,jdm,kdm
      use hycom_dim, only: ogrid
      use hycom_arrays, only : tracer,dpinit,scp2
      use hycom_scalars, only : nstep
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_GASEXCH_COM, only : tracflx
#endif
      use domain_decomp_1d, only: am_i_root, globalsum, get

      implicit none

      integer i,j,k,l
      integer ntr
      real sumo, areao

      integer :: iTracer
      integer :: j_0, j_1, j_0h, j_1h

      real*8, allocatable :: summ(:)
      real*8, allocatable :: partialTracerSums(:,:)
      real*8, allocatable :: partialfluxsum(:)
      real*8 :: sumflux

      call get(ogrid, j_strt = j_0, j_stop = j_1,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h)

      call getOceanVolumeAndArea(sumo, areao)

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

      !integrate flux
      allocate(partialfluxsum(j_0h:j_1h))
      partialfluxsum(:) = partialIntegration(tracflx(:,:,1:1))

      call globalSum(ogrid, partialfluxsum, sumflux)

      if (am_i_root()) then
       write(*,*)'global averaged flux=',sumflux/areao
      endif

      deallocate(partialTracerSums, summ)
      deallocate(partialfluxsum)

      contains

#ifdef OBIO_ON_GARYocean
      function partialIntegration(quantity)
      use ocean, only : dxypo, focean
      use oceanres, only: dzo
      real*8, intent(in) :: quantity(:,j_0h:,:)
      real*8 :: partialIntegration(j_0h:j_1h)
      
      real*8 :: gridCellVolume
      integer :: numLevels

      partialIntegration = 0
      numLevels = size(quantity,3)
      do k = 1, numLevels
         do j = j_0, j_1
            gridCellVolume = dzo(k) * dxypo(j)
            do i= 1, idm
               if (focean(i,j) > 0) then
                  partialIntegration(j) = partialIntegration(j) + 
     &                 quantity(i,j,k) * gridCellVolume
               end if
            end do
         end do
      end do

      end function partialIntegration
#else
      function partialIntegration(quantity)
      use hycom_scalars, only : huge
      use hycom_dim, only: ifp, ilp, isp
      real*8, intent(in) :: quantity(:,j_0h:,:)
      real*8 :: partialIntegration(j_0h:j_1h)
      
      integer :: i, j, l, k
      integer :: numLevels

      partialIntegration = 0
      numLevels = size(quantity,3)
      do k = 1, numLevels
         do j = j_0, j_1
            do i = 1, idm
               if (dpinit(i,j,k) < huge) then
                  partialIntegration(j) = partialIntegration(j) + 
     &                 quantity(i,j,k) * dpinit(i,j,k)*scp2(i,j)
               end if
            end do
         end do
      end do

      end function partialIntegration
#endif

      subroutine getOceanVolumeAndArea(oceanVolume, oceanArea)
      real*8, intent(out) :: oceanVolume
      real*8, intent(out) :: oceanArea

      real*8, allocatable :: partialOceanVolumes(:)
      real*8, allocatable :: ones(:,:,:)

      allocate(partialOceanVolumes(j_0h:j_1h))
      allocate(ones(idm, J_0H:J_1H, kdm))
      
      ones = 1
      partialOceanVolumes(:) = partialIntegration(ones)
      call globalSum(ogrid, partialOceanVolumes, oceanVolume)

      partialOceanVolumes(:) = partialIntegration(ones(:,:,1:1))
      call globalSum(ogrid, partialOceanVolumes, oceanArea)

      deallocate(ones)
      deallocate(partialOceanVolumes)

      end subroutine getOceanVolumeAndArea

      end subroutine obio_trint
