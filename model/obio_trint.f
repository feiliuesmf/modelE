#include "rundeck_opts.h"

      subroutine obio_trint(iflg)
!iflg - 0 compute  trac_old  (before obio)
!iflg - 1 zero-out trac_old  (after  obio)

#ifdef OBIO_ON_GARYocean
      use oceanres,  only: idm=>imo,jdm=>jmo,kdm=>lmo
      use oceanr_dim, only : ogrid
      use obio_com, only: tracer => tracer_loc ! rename local
      use ocn_tracer_com, only : ntrcr=>ntm
      use model_com, only : nstep=>itime
      use ocean, only: trmo
#else
      use hycom_dim_glob, only : ntrcr,idm,jdm,kdm
      use hycom_dim, only: ogrid
      use hycom_arrays, only : tracer,dpinit,scp2
      use hycom_scalars, only : nstep,onem
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_GASEXCH_COM, only : tracflx
#endif
      USE obio_incom, only: mgchltouMC
      USE obio_com,   only: trac_old,obio_deltath
      use domain_decomp_1d, only: am_i_root, globalsum, get

      implicit none

      integer i,j,k,l
      integer ntr
      real sumo, areao, glb_carbon_invntry

      integer :: iTracer, iflg
      integer :: j_0, j_1, j_0h, j_1h

      real*8, allocatable :: summ(:)
      real*8 :: sumFlux(1), fluxNorm1, fluxNorm2
      real*8 :: sumFlux2(1)

      call get(ogrid, j_strt = j_0, j_stop = j_1,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h)

      sumo =  getTotalOceanVolume()
      areao =  getTotalOceanArea()

      ! Form partial sum for each latitude for each tracer
      allocate(summ(ntrcr))
      summ = volumeIntegration(tracer)

      !integrate flux
      ! using resize to force tracflx to act as 4D array with
      ! size 1 in the uninteresting directions to match
      ! expected interface.
      sumFlux = volumeIntegration(reshape(tracflx(:,:,1),
     &     (/size(tracflx,1), size(tracflx,2), 1, 1/) ))
      fluxNorm1 = avgDepthOfTopLayer()
      fluxNorm2 = sumDepthOfTopLayer()
!     sumFlux = sumFlux / fluxNorm

      sumFlux2= areaIntegration(tracflx(:,:,1))

      if (am_i_root()) then
         do iTracer = 1, ntrcr
            write(*,*) 'total intgrl tracer:', iTracer, nstep, 
     &           summ(iTracer), summ(iTracer)/sumo
         end do

       write(*,*)'---------- tracer conservation--------------------'
       print*, 'mgchltouMC=',mgchltouMC
       print*, 'area, volume ocean=', areao, sumo
       write(*,*)'global averaged flux=',nstep,sumFlux, sumFlux/areao
       write(*,*)'global averaged flux2=',
     .            nstep,sumFlux/fluxNorm1,sumFlux/fluxNorm2
       write(*,*)'global averaged flux3=',
     .            nstep,sumFlux2
       !volume integrated carbon inventory:
        glb_carbon_invntry = 
     .                     (summ(5)+summ(6)
     .                    + summ(7)+summ(8)
     .                    + summ(9)) * mgchltouMC * 1e-3     !mgm3*m3 -> uM*m3=mili,molC -> mol,C
     .                    + summ(11) *1e-3 /12               !micro-grC/lt*m3 -> mili,grC-> mol,C
     .                    +(summ(14)+summ(15))*1e-3          !mol,C (or mol,CO2)
       write(*,*)'glb carbon inventory bfre=',nstep,trac_old
       write(*,*)'glb carbon inventory aftr=',nstep,glb_carbon_invntry
       if (iflg.eq.0) then
        trac_old = glb_carbon_invntry
       else
        write(*,'(a,i5,1x,2(e18.11,1x))')'carbon conserv.',nstep,
     .       (glb_carbon_invntry-trac_old)/(obio_deltath*3600.d0),
     .       sumFlux2
        write(*,'(a,i5,1x,e18.11)')'carbon residual',nstep,
     .       (glb_carbon_invntry-trac_old)/(obio_deltath*3600.d0)
     .       -sumFlux2 
        trac_old = 0.d0
       endif

       write(*,*)'----------------------------------------------------'
      endif

#ifdef OBIO_ON_GARYocean
      summ = volumeIntegration(trmo)

      if (am_i_root()) then
         do iTracer = 1, ntrcr
            write(*,*) 'total intgrl trmo:', iTracer, nstep, 
     &           summ(iTracer), summ(iTracer)/sumo
         end do
       write(*,*)'----------   trmo conservation--------------------'
       write(*,*)'global averaged flux=',nstep,sumFlux, sumFlux/areao
       !volume integrated carbon inventory:
       !trmo units are Kg
       write(*,*)'global carbon inventory=',nstep,
     .                                     ((summ(5)+summ(6)
     .                                       +summ(7)+summ(8)
     .                                       +summ(9))* mgchltouMC   !mol,C
     .                                       +summ(11)
     .                                       +summ(14)+summ(15))     !mol,C
       write(*,*)'----------------------------------------------------'
      endif
#endif

#ifdef OBIO_ON_GARYocean
      summ = volumeIntegration(trmo)

      if (am_i_root()) then
         do iTracer = 1, ntrcr
            write(*,*) 'total intgrl trmo:', iTracer, nstep, 
     &           summ(iTracer), summ(iTracer)/sumo
         end do
       write(*,*)'----------   trmo conservation--------------------'
       write(*,*)'global averaged flux=',nstep,sumFlux, sumFlux/areao
       !volume integrated carbon inventory:
       !trmo units are Kg
       write(*,*)'global carbon inventory=',nstep,
     .                                      ((summ(5)+summ(6)
     .                                       +summ(7)+summ(8)
     .                                       +summ(9))* mgchltouMC   !mol,C
     .                                       +summ(11)
     .                                       +summ(14)+summ(15))     !mol,C
       write(*,*)'----------------------------------------------------'
      endif
#endif

      deallocate(summ)

      contains

      real*8 function avgDepthOfTopLayer()
      real*8, allocatable :: ones(:,:,:,:)
      real*8 :: volumeOfTopLayer(1), area

      allocate(ones(idm, J_0H:J_1H, 1, 1))
      ones = 1
      volumeOfTopLayer = volumeIntegration(ones)
      area = areaIntegration(ones(:,:,1,1))
      deallocate(ones)

      avgDepthOfTopLayer = volumeOfTopLayer(1)/area
      end function avgDepthOfTopLayer

#ifdef OBIO_ON_GARYocean
      function volumeIntegration(quantity)
      use ocean, only : dxypo, focean
      use oceanres, only: dzo
      real*8, intent(in) :: quantity(:,j_0h:,:,:)
      real*8 :: volumeIntegration(size(quantity,4))
      
      real*8 :: partialIntegration(j_0h:j_1h,size(quantity,4))
      real*8 :: gridCellVolume
      integer :: numLevels, numTracers
      integer :: i,j,k,n

      partialIntegration = 0
      numLevels = size(quantity,3)
      numTracers = size(quantity,4)

      do k = 1, numLevels
        do j = j_0, j_1
          gridCellVolume = dzo(k) * dxypo(j)
          do i= 1, idm
            do n = 1, numTracers
               if (focean(i,j) > 0) then
                 partialIntegration(j,n) = partialIntegration(j,n) + 
     &                 quantity(i,j,k,n) * gridCellVolume
              end if
            end do
          end do
        end do
      end do
      
      call globalSum(ogrid, partialIntegration, volumeIntegration,
     &     all=.true.)

      end function volumeIntegration

      real*8 function areaIntegration(quantity)
      use ocean, only : dxypo, focean
      real*8, intent(in) :: quantity(:,j_0h:)
      real*8 :: partialAreaIntegration(j_0h:j_1h)
      
      real*8 :: gridCellArea
      integer :: numLevels

      integer :: i, j

      partialAreaIntegration = 0
      do j = j_0, j_1
        gridCellArea = dxypo(j)
        do i= 1, idm
          if (focean(i,j) > 0) then
            partialAreaIntegration(j) = partialAreaIntegration(j) + 
     &            quantity(i,j) * gridCellArea
          end if
        end do
      end do

      call globalSum(ogrid, partialAreaIntegration, areaIntegration,
     &     all=.true.)

      end function areaIntegration

      real*8 function sumDepthOfTopLayer()
      use oceanres, only: dzo
      use ocean, only : focean
      real*8 :: partialIntegration(j_0h:j_1h)
      
      integer :: i, j

      partialIntegration = 0
      do j = j_0, j_1
        do i= 1, idm
          if (focean(i,j) > 0) then
            partialIntegration(j) = partialIntegration(j) + dzo(1)
          end if
        end do
      end do

      call globalSum(ogrid, partialIntegration, sumDepthOfTopLayer, 
     &     all=.true.)

      end function sumDepthOfTopLayer
#else
      function volumeIntegration(quantity)
      use hycom_scalars, only : huge
      real*8, intent(in) :: quantity(:,j_0h:,:,:)
      real*8 :: volumeIntegration(size(quantity,4))
      
      real*8 :: partialVolumeIntegration(j_0h:j_1h, size(quantity,4))
      integer :: i, j, l, k, n
      integer :: numLevels, numTracers

      partialVolumeIntegration = 0
      numLevels = size(quantity,3)
      numTracers = size(quantity,4)

      do k = 1, numLevels
        do j = j_0, j_1
          do i = 1, idm
            if (dpinit(i,j,k) < huge) then
              do n = 1, numTracers
                 partialVolumeIntegration(j,n) = 
     &                partialVolumeIntegration(j,n) + 
     &                quantity(i,j,k,n)
     &                    *dpinit(i,j,k)/onem*scp2(i,j)
              end do
            end if
          end do
        end do
      end do

      call globalSum(ogrid, partialVolumeIntegration, volumeIntegration,
     &     all=.true.)

      end function volumeIntegration

      real*8 function areaIntegration(quantity)
      use hycom_scalars, only : huge
      real*8, intent(in) :: quantity(:,j_0h:)
      real*8 :: partialAreaIntegration(j_0h:j_1h)
      
      integer :: i, j, l

      partialAreaIntegration = 0
      do j = j_0, j_1
         do i = 1, idm
            if (dpinit(i,j,1) < huge) then
               partialAreaIntegration(j) = partialAreaIntegration(j) + 
     &              quantity(i,j) * scp2(i,j)
            end if
         end do
      end do

      call globalSum(ogrid, partialAreaIntegration, areaIntegration,
     &     all=.true.)

      end function areaIntegration

      real*8 function sumDepthOfTopLayer()
      use hycom_scalars, only : huge
      
      real*8 :: partialIntegration(j_0h:j_1h)
      integer :: i, j

      partialIntegration = 0
      do j = j_0, j_1
        do i = 1, idm
          if (dpinit(i,j,1) < huge) then
            partialIntegration(j) = partialIntegration(j)
     .                            + dpinit(i,j,1)/onem
          end if
        end do
      end do

      call globalSum(ogrid, partialIntegration, sumDepthOfTopLayer,
     &     all=.true.)

      end function sumDepthOfTopLayer

#endif

      real*8 function getTotalOceanVolume() result(oceanVolume)
      real*8, allocatable :: ones(:,:,:,:)

      real*8 :: volumeAsArray(1)
      allocate(ones(idm, J_0H:J_1H, kdm, 1))

      ones = 1
      volumeAsArray = volumeIntegration(ones)
      oceanVolume = volumeAsArray(1)
      deallocate(ones)

      end function getTotalOceanVolume

      real*8 function getTotalOceanArea() result(oceanArea)
      real*8, allocatable :: ones(:,:)

      allocate(ones(idm, J_0H:J_1H))

      ones = 1
      oceanArea = areaIntegration(ones)
      deallocate(ones)

      end function getTotalOceanArea

      end subroutine obio_trint
