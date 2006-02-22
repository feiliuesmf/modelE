#include "rundeck_opts.h"

      module ent_drv
!@sum veg_drv contains variables and routines for vegetation driver
!@auth I. Alienov, N. Kiang

      use model_com, only : im,jm
      implicit none
      private
      save

      public init_vegetation,reset_veg_to_defaults
     &     ,veg_set_cell, veg_save_cell,upd_gh

      real*8,public :: cosday,sinday

      contains



      subroutine init_module_ent(iniENT, Jday, Jyear)
!@sum initializes vegetation
      use param
      use ent_GISSveg, only :: ent_GISS_init
      use ent_com, only :: entcells
      use DOMAIN_DECOMP, only : GRID, GET
      integer, intent(in) :: Jday, Jyear
      logical, intent(in) :: iniENT
      !---
      integer I_0, I_1, J_0, J_1



      if (iniENT ) then
        call ent_GISS_init( entcells(I_0:I_1,J_0:J_1),
     &       IM, JM, I_0, I_1, J_0, J_1, jday, year )
      endif

      call read_veg_data(redogh,istart)
      end subroutine init_vegetation


      subroutine read_veg_data(redogh,istart)
!@sum reads vegetation arrays and rundeck parameters
      use filemanager
      use param
      use DOMAIN_DECOMP, only : GRID, GET, READT_PARALLEL
      use vegetation, only : cond_scheme,vegCO2X_off,crops_yr
      use veg_com
      use model_com, only : fearth,jyear

      implicit none

      integer, intent(in) :: istart
      logical, intent(in) :: redogh

      INTEGER :: J_1, J_0, J_1H, J_0H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG

      integer iu_veg
      integer i, j, k
      logical veg_data_missing

      real*8, allocatable :: veg_c4(:,:)
      real*8 :: vc4
      integer :: read_c4_grass = 0

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

c**** read rundeck parameters
      call sync_param( "cond_scheme", cond_scheme)  !nyk 5/1/03
      call sync_param( "vegCO2X_off", vegCO2X_off)  !nyk 3/2/04
      call sync_param( "crops_yr", crops_yr)
      call sync_param( "read_c4_grass", read_c4_grass)

c**** read land surface parameters or use defaults
      call openunit("VEG",iu_VEG,.true.,.true.)
      do k=1,10                 !  11 ????
        CALL READT_PARALLEL
     *    (grid,iu_VEG,NAMEUNIT(iu_VEG),0,vdata(:,:,K),1)
      end do
c**** zero-out vdata(11) until it is properly read in
      vdata(:,:,11) = 0.
      vdata(:,:,12) = 0.
      call closeunit(iu_VEG)
c**** add data on c4 grass
      if ( read_c4_grass >= 1 ) then
        print *,"Adding c4 grass data"
        call openunit("VEG_C4",iu_VEG,.true.,.true.)
        allocate( veg_c4(im,J_0H:J_1H) )
        CALL READT_PARALLEL
     *       (grid,iu_VEG,NAMEUNIT(iu_VEG),0,veg_c4(:,:),1)
        do j=J_0,J_1
          do i=1,im
            if (fearth(i,j).gt.0) then
              ! normalize to earth fraction
              vc4 = veg_c4(i,j) / fearth(i,j)
              ! add c4 grass up to max amount of current grass
              vdata(i,j,12) = min( vdata(i,j,3), vc4 )
              vdata(i,j,3) = vdata(i,j,3) - vdata(i,j,12)
            end if
          enddo
        enddo      
        deallocate( veg_c4 )
        call closeunit(iu_VEG)
      endif
C**** Update vegetation file if necessary (i.e. crops_yr =0 or >0)
      if(crops_yr.eq.0) call updveg(jyear,.false.)
      if(crops_yr.gt.0) call updveg(crops_yr,.false.)

      if (istart.le.2 .or. redogh) then ! initial. foliage arrays (adf)
        Cint(:,:)=0.0127D0      ! internal CO2
        Qfol(:,:)=3.D-6         ! surface mixing ratio
        cnc_ij(:,:) = 0.d0
      end if
      if (istart.le.0) return   ! avoid reading unneeded files

     !!! spgsn=.1d0

c**** check whether ground hydrology data exist at this point.
      veg_data_missing = .false.
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).gt.0) then
            if ( sum(vdata(i,j,1:12)).eq.0 ) then
              print *,"No vegetation data: i,j=",i,j,vdata(i,j,1:12)
              veg_data_missing = .true.
            end if
          end if
        enddo
      enddo
      if ( veg_data_missing ) then
        write(6,*) 'Vegetation data is missing at some pts'
        write(6,*) 'If you have a non-standard land mask, please'
        write(6,*) 'consider using extended GH data and rfs file.'
        call stop_model(
     &       'Vegetation data is missing at some cells',255)
      endif

      end subroutine read_veg_data




      subroutine reset_veg_to_defaults( reset_prognostic )
      use veg_com, only: vdata
      use DOMAIN_DECOMP, only : GRID, GET
      logical, intent(in) :: reset_prognostic
      integer i,j

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)

      do j=J_0,J_1
      do i=1,im

      vdata(i,j,1:12)= (/  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.62451148d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.37548852d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00 /)
      enddo
      enddo

      end subroutine reset_veg_to_defaults



