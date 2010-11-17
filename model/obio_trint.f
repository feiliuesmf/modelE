#include "rundeck_opts.h"

      subroutine obio_trint(iflg)
!iflg - 0 compute  carb_old  (before obio)
!iflg - 1 zero-out carb_old  (after  obio)

#ifdef OBIO_ON_GARYocean
      use oceanres,  only: idm=>imo,jdm=>jmo,kdm=>lmo
      use oceanr_dim, only : ogrid
      use obio_com, only: tracer => tracer_loc ! rename local
      use ocn_tracer_com, only : ntrcr=>ntm, obio_tr_mm
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
      USE MODEL_COM, only: JMON
      USE obio_incom, only: mgchltouMC,solFe
      USE obio_forc, only: atmFe
      USE obio_com,   only: carb_old,obio_deltath,iron_old,p1d
      use domain_decomp_1d, only: am_i_root, globalsum, get

      implicit none

      integer i,j,k,l
      integer ntr
      real sumo,areao,glb_carbon_invntry,ironFlux,glb_iron_invntry
      real h1

      integer :: iTracer, iflg
      integer :: j_0, j_1, j_0h, j_1h

      real*8, allocatable :: summ(:)
      real*8 :: sumFlux(1)

      call get(ogrid, j_strt = j_0, j_stop = j_1,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h)

      sumo  =  getTotalOceanVolume()
      areao =  getTotalOceanArea()
        h1  = avgDepthOfTopLayer()

      ! Form partial sum for each latitude for each tracer
      allocate(summ(ntrcr))
      summ = volumeIntegration(tracer)

#ifdef TRACERS_GASEXCH_ocean
      !integrate flux
      ! using resize to force tracflx to act as 4D array with
      ! size 1 in the uninteresting directions to match
      ! expected interface.
#ifdef OBIO_ON_GARYocean
      sumFlux= areaIntegration(tracflx(:,:,1))
#else
      sumFlux= areaIntegration(tracflx(:,:,1))
#endif
#else
      sumFlux=0.   ! no surface flux
#endif
      ironFlux= areaIntegration(atmFe(:,:,JMON))
#ifdef zero_ironflux
      ironFlux=0.d0
#endif
      ironFlux= ironFlux*solFe*1.d-3/max(h1,1.d-3)

      if (am_i_root()) then
         do iTracer = 1, ntrcr
            write(*,*) 'total intgrl tracer:', iTracer, nstep, 
     &           summ(iTracer), summ(iTracer)/sumo
         end do

       write(*,*)'---------- carbon conservation--------------------'
       print*, 'area, volume ocean=', areao, sumo
       write(*,*)'global averaged carbon flux=',
     .            nstep,sumFlux
       !volume integrated carbon inventory:
        glb_carbon_invntry = 
     .                     (summ(5)+summ(6)
     .                    + summ(7)+summ(8)
     .                    + summ(9)) * mgchltouMC * 1e-3     !mgm3*m3 -> uM*m3=mili,molC -> mol,C
     .                    + summ(11) *1e-3 /12               !micro-grC/lt*m3 -> mili,grC-> mol,C
     .                    +(summ(14)+summ(15))*1e-3          !mol,C (or mol,CO2)
       write(*,*)'glb carbon inventory bfre=',nstep,carb_old
       write(*,*)'glb carbon inventory aftr=',nstep,glb_carbon_invntry
       if (iflg.eq.0) then
        carb_old = glb_carbon_invntry
       else
        write(*,'(a,i5,1x,2(e18.11,1x))')'carbon conserv.',nstep,
     .       (glb_carbon_invntry-carb_old)/(obio_deltath*3600.d0),
     .       sumFlux
        write(*,'(a,i5,1x,e18.11)')'carbon residual (no units)',nstep,
     .       ((glb_carbon_invntry-carb_old)/(obio_deltath*3600.d0)
     .       -sumFlux)/glb_carbon_invntry
        carb_old = 0.d0
       endif
       write(*,*)'----------------------------------------------------'

       write(*,*)'-------------- iron conservation--------------------'
       print*, 'area, volume ocean=', areao, sumo
       write(*,*)'global averaged iron deposition flux=',
     .            nstep,ironFlux
       !volume integrated iron inventory:
        glb_iron_invntry = (summ(4) + summ(13))    !nano-mol,Iron/m3
       write(*,*)'glb iron inventory bfre=',nstep,iron_old
       write(*,*)'glb iron inventory aftr=',nstep,glb_iron_invntry
       if (iflg.eq.0) then
        iron_old = glb_iron_invntry
       else
        write(*,'(a,i5,1x,2(e18.11,1x))')'iron conserv.',nstep,
     .       (glb_iron_invntry-iron_old)/(obio_deltath*3600.d0),
     .       ironFlux
        write(*,'(a,i5,1x,e18.11)')'iron residual (no units)',nstep,
     .       ((glb_iron_invntry-iron_old)/(obio_deltath*3600.d0)
     .       -ironFlux)/glb_iron_invntry
        iron_old = 0.d0
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
       write(*,*)'global averaged flux (kg,C per s) =',
     .            nstep,sumFlux*44.d0*1e-3/3.664d0     !convert mole,CO2/s->kg,C/s
       !volume integrated carbon inventory:
       !trmo units are Kg, summ is in Kg*m3
       write(*,*)'global carbon inventory (kg,C)=',nstep,
     .   (summ(5)+summ(6)+summ(7)+summ(8)+summ(9)  ! kg,C*m3
     . +  summ(11)                                 ! kg,C*m3
     . +  summ(14)+summ(15)                        ! kg,C*m3
     .                     )/sumo                  ! kg,C
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
      use ocean, only : dxypo, focean,imaxj
      use oceanres,  only: dzo
      real*8, intent(inout) :: quantity(:,j_0h:,:,:)
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
          if (j.eq.jdm) gridCellVolume = gridCellVolume * idm
          do i= 1, imaxj(j)
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
      use ocean, only : dxypo, focean,imaxj
      real*8, intent(in) :: quantity(:,j_0h:)
      real*8 :: partialAreaIntegration(j_0h:j_1h)
      
      real*8 :: gridCellArea
      integer :: numLevels

      integer :: i, j

      partialAreaIntegration = 0
      do j = j_0, j_1
        gridCellArea = dxypo(j)
        if (j.eq.jdm) gridCellArea = gridCellArea *idm
        do i= 1, imaxj(j)
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
      use oceanres,  only: dzo
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
