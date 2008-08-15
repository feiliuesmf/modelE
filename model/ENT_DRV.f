#include "rundeck_opts.h"

      module ent_drv
!@sum ent_drv contains variables and routines for vegetation driver
!@auth I. Alienov, N. Kiang

      use model_com, only : im,jm
      use ent_mod
      implicit none
      private
      save

      public init_module_ent, update_vegetation_data

      contains



      subroutine init_module_ent(iniENT, Jday, Jyear, focean1)
!@sum initializes vegetation
      use param
      use ent_com, only : entcells,Cint,Qfol,cnc_ij
      use model_com, only : focean
      use DOMAIN_DECOMP, only : GRID, GET
      integer, intent(in) :: Jday, Jyear
      logical, intent(inout) :: iniENT
      real*8, intent(in) :: focean1(:,:)
      !---
      integer I_0, I_1, J_0, J_1, i, j
      ! the following are rundeck parameters which need to be
      ! passed to ent (and used there)
      integer cond_scheme, vegCO2X_off, crops_yr, read_c4_grass
      integer :: force_init_ent=0

      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               I_STRT     =I_0,    I_STOP     =I_1)

!!! hack
      call sync_param( "init_ent", force_init_ent)
      if ( force_init_ent == 1 ) iniENT = .true.

      !--- read rundeck parameters
      call sync_param( "cond_scheme", cond_scheme)  !nyk 5/1/03
      call sync_param( "vegCO2X_off", vegCO2X_off)  !nyk 3/2/04
      call sync_param( "crops_yr", crops_yr)
      call sync_param( "read_c4_grass", read_c4_grass)


      ! maybe call "ent_initialize" ? , i.e.
      ! call ent_initialize(cond_scheme,vegCO2X_off,crops_yr,nl_soil, etc)
      ! ask Max if OK
      call ent_initialize(
     &     do_soilresp=.true.
     &     ,do_phenology=.false.
     &     ,do_frost_hardiness=.false.
     &     ,do_patchdynamics=.false.
     &     )

      if (iniENT ) then
      ! initialize ent cells to something meaningful
        do j=J_0,J_1
          do i=I_0,I_1
            if (focean(i,j) <= 0) then
              !print *,"EARTH",i,j,focean(i,j)
              call ent_cell_construct( entcells(i,j) )
            else
              !print *,"OCEAN",i,j,focean(i,j)
            end if
          enddo
        enddo

        !print *,'init_module_ent printing cells 1'
        !call ent_cell_print(6, entcells)
        !print *,'init_module_ent end printing cells 1'

        call set_vegetation_data( entcells, ! (I_0:I_1,J_0:J_1),
     &       IM, JM, I_0, I_1, J_0, J_1, jday, jyear )

        ! probably it is ok to initialize these here since 
        ! if .not. iniENT we should read everything from restart file
        Cint(:,:)=0.0127D0      ! internal CO2
        Qfol(:,:)=3.D-6         ! surface mixing ratio
        cnc_ij(:,:) = 0.d0

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
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
      use ent_prescr_veg, only : prescr_calc_shc,prescr_calcconst
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpool_ini
      !-----Local---------
      integer i,j
      real*8 heat_capacity


cddd      do j=j0,j1
cddd        do i=i0,i1
cddd          print *,'set_vegetation_data i,j = ',i,j
cddd          call ent_cell_print( 6, entcells(i,j) )
cddd        enddo
cddd      enddo

!      print *,'set_vegetation_data printing cells 1'
!      call ent_cell_print(entcells)
!      print *,'set_vegetation_data end printing cells 1'


      !Read land surface parameters or use defaults
      !GISS data sets:
!      call GISS_vegdata(jday, year, 
!     &     im,jm,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
!     &     frootdata,popdata,soildata,soil_texture)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
!      call ent_cell_set(entcells, vegdata, popdata, laidata,
!     &     hdata, nmdata, frootdata, soildata, albedodata, soil_texture)

      call prescr_calcconst()

      call prescr_vegdata(jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpool_ini,.true.)
      !print *,"popdata in ent_GISS_init: ",popdata
      !vegdata(1:2,:,:) = 0.0
      !vegdata(3,:,:) = 1.0
      !vegdata(4:N_COVERTYPES,:,:) = 0.0
      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      !!! hack
      !!!Tpool_ini = 0.d0

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soildata, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini)

      !!! hack to set constant heat capacity
cddd      do j=j0,j1
cddd        do i=i0,i1
cddd          heat_capacity = GISS_calc_shc( vegdata(1:N_COVERTYPES,i,j ) )
cddd          call ent_cell_update(entcells(i,j),
cddd     &         heat_capacity=heat_capacity
cddd     &         )
cddd        enddo
cddd      enddo
      
      ! just in case, do nothing, just set heat capacities
      call ent_prescribe_vegupdate(entcells)


      end subroutine set_vegetation_data


      subroutine update_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata
      !use ent_prescr_veg, only: prescr_get_laidata,prescr_veg_albedodata
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      !-----Local---------
      integer i,j


      call prescr_get_laidata(jday,JM,I0,I1,J0,J1,laidata)
      call prescr_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)

cddd      do j=j0,j1
cddd        do i=i0,i1
cddd          call ent_cell_update(entcells(i,j),
cddd     &         pft_leaf_area_index=laidata(2:2+N_PFT-1,i,j),
cddd     &         pft_vegalbedo=albedodata(:,2:2+N_PFT-1,i,j)
cddd     &         )
cddd        enddo
cddd      enddo

      call ent_prescribe_vegupdate(entcells
     &     ,laidata=laidata(2:2+N_PFT-1,:,:)
     &     ,albedodata=albedodata(:,2:2+N_PFT-1,:,:)
     &     )

      end subroutine update_vegetation_data


      end module ent_drv
