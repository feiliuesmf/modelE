#include "rundeck_opts.h"

      module ent_drv
!@sum ent_drv contains variables and routines for vegetation driver
!@auth I. Alienov, N. Kiang

      use model_com, only : im,jm
      use ent_mod
      implicit none
      private
      save

      public init_module_ent

      contains



      subroutine init_module_ent(iniENT, Jday, Jyear)
!@sum initializes vegetation
      use param
      use ent_com, only : entcells
      use DOMAIN_DECOMP, only : GRID, GET
      integer, intent(in) :: Jday, Jyear
      logical, intent(in) :: iniENT
      !---
      integer I_0, I_1, J_0, J_1
      ! the following are rundeck parameters which need to be
      ! passed to ent (and used there)
      integer cond_scheme, vegCO2X_off, crops_yr, read_c4_grass

      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               I_STRT     =I_0,    I_STOP     =I_1)

      !--- read rundeck parameters
      call sync_param( "cond_scheme", cond_scheme)  !nyk 5/1/03
      call sync_param( "vegCO2X_off", vegCO2X_off)  !nyk 3/2/04
      call sync_param( "crops_yr", crops_yr)
      call sync_param( "read_c4_grass", read_c4_grass)


      ! maybe call "ent_initialize" ? , i.e.
      ! call ent_initialize(cond_scheme,vegCO2X_off,crops_yr,nl_soil, etc)
      ! ask Max if OK

      if (iniENT ) then
        call set_vegetation_data( entcells(I_0:I_1,J_0:J_1),
     &       IM, JM, I_0, I_1, J_0, J_1, jday, jyear )
      endif

      ! the following parts of the code should be implemented
      ! somewhere in dynamic vegetation
#ifdef UNFINISHED_CODE
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

#endif
      end subroutine init_module_ent


      subroutine set_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_GISSveg, only : GISS_vegdata
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: frootdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TYPES,I0:I1,J0:J1) :: soil_texture
      !-----Local---------

      !Read land surface parameters or use defaults
      !GISS data sets:
      call GISS_vegdata(jday, year, 
     &     im,jm,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     frootdata,popdata,soildata,soil_texture)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, nmdata, frootdata, soildata, albedodata, soil_texture)

      end subroutine set_vegetation_data


      end module ent_drv
