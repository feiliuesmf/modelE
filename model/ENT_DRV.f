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

      logical :: initialized = .false.
      integer :: crops_yr = 0

      contains



      subroutine init_module_ent(iniENT_in, Jday, Jyear, focean1)
!@sum initializes vegetation
      use param
      use ent_com, only : entcells,Cint,Qfol,cnc_ij
      use ent_prescr_veg, only : prescr_calc_shc,prescr_calcconst
      use model_com, only : focean, FLICE
      use DOMAIN_DECOMP_ATM, only : GRID, GET
      integer, intent(in) :: Jday, Jyear
      logical, intent(in) :: iniENT_in
      real*8, intent(in) :: focean1(:,:)
      !---
      integer I_0, I_1, J_0, J_1, i, j
      ! the following are rundeck parameters which need to be
      ! passed to ent (and used there)
      integer cond_scheme, vegCO2X_off, read_c4_grass, year
      integer :: force_init_ent=0
      logical iniENT

      if ( initialized ) return
      initialized = .true.

      iniENT = iniENT_in

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

      if ( crops_yr .ne. 0 ) then
        year = crops_yr
      else
        year = Jyear
      endif

      ! maybe call "ent_initialize" ? , i.e.
      ! call ent_initialize(cond_scheme,vegCO2X_off,crops_yr,nl_soil, etc)
      ! ask Max if OK
      call ent_initialize(
     &     do_soilresp=.true.
     &     ,do_phenology_activegrowth=.false.
     &     ,do_frost_hardiness=.true. ! .false.
     &     ,do_patchdynamics=.false.
     &     )

      call prescr_calcconst() ! moved above

      if (iniENT ) then
      ! initialize ent cells to something meaningful
        do j=J_0,J_1
          do i=I_0,I_1
            !!if (focean(i,j) <= 0) then
            if ( focean(i,j) < 1.d0 ) then
              !print *,"EARTH",i,j,focean(i,j)
              !print *,"ent_construct",i,j
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
     &       IM, JM, I_0, I_1, J_0, J_1, jday, year )

        ! probably it is ok to initialize these here since 
        ! if .not. iniENT we should read everything from restart file
        Cint(:,:)=0.0127D0      ! internal CO2
        Qfol(:,:)=3.D-6         ! surface mixing ratio
        cnc_ij(:,:) = 0.d0

      else ! i.e. *not* iniENT

        ! just in case, do nothing, just set heat capacities
        call ent_prescribe_vegupdate(entcells)        

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
        do i=I_0,I_1
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
      use DOMAIN_DECOMP_ATM, only : GRID
      use geom, only : lat2d
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata

      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata,
     &     prescr_get_hdata,prescr_get_woodydiameter,prescr_get_pop,
     &     prescr_get_crownrad,prescr_get_carbonplant,prescr_get_initnm,
     &     prescr_get_rootprof,prescr_get_soilcolor,
     &     prescr_get_soilpools,prescr_get_soil_C_total


      use ent_prescr_veg, only : prescr_calc_shc,prescr_calcconst
      use ghy_com, only : q_ij, qk_ij, dz_ij
      use ghy_com, only : fearth
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soil_color ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      !real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture1
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpool_ini
      !-----Local---------
      integer i,j,k
      real*8 heat_capacity
      real*8 :: vdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      real*8 :: cropdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 :: soil_C_total_H(N_CASA_LAYERS,
     &     grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      integer hemi(I0:I1,J0:J1)

            !* Set hemisphere flags.
      where(lat2d(I0:I1,J0:J1) <= 0.)
        hemi(I0:I1,J0:J1) = -1  ! S
      elsewhere
        hemi(I0:I1,J0:J1) = +1  ! N
      end where

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

      !call prescr_calcconst() ! moved above

cddd      call prescr_vegdata(jday, year, 
cddd     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
cddd     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
cddd     &     soil_color,soil_texture,Tpool_ini,.true.,.false.,.false.)
      !print *,"popdata in ent_GISS_init: ",popdata
      !vegdata(1:2,:,:) = 0.0
      !vegdata(3,:,:) = 1.0
      !vegdata(4:N_COVERTYPES,:,:) = 0.0
      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      !!! hack
      !!!Tpool_ini = 0.d0

      vegdata=0
      albedodata=0
      laidata=0
      hdata=0
      nmdata=0

      popdata=0
      dbhdata=0
      craddata=0
      cpooldata=0
      rootprofdata=0

      soil_color=0
      soil_texture=0
      Tpool_ini=0

      !Read vegdata
      call get_vdata(vdata_H)
      call get_cropdata(year, cropdata_H)
      do k=1,N_COVERTYPES
        vegdata(k,I0:I1,J0:J1) = vdata_H(I0:I1,J0:J1,k)
     &       *(1.d0-cropdata_H(I0:I1,J0:J1))
      enddo
      vegdata(9,I0:I1,J0:J1) = cropdata_H(I0:I1,J0:J1)



      call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)


         do j=J0,J1
           do i=I0,I1
             call prescr_get_hdata(hdata(:,i,j)) !height
             call prescr_get_woodydiameter(hdata(:,i,j), dbhdata(:,i,j))
             call prescr_get_pop(dbhdata(:,i,j), popdata(:,i,j))
             call prescr_get_crownrad(popdata(:,i,j), craddata(:,i,j))
           enddo
         enddo
         call prescr_get_carbonplant(I0,I1,J0,J1,
     &        laidata,hdata,dbhdata,popdata,cpooldata)
      call prescr_get_initnm(nmdata) !nm ! mean canopy nitrogen
      call prescr_get_rootprof(rootprofdata)
      call prescr_get_soilcolor(soil_color)

#ifdef SET_SOILCARBON_GLOBAL_TO_ZERO
      Tpool_ini = 0.d0
#else
      !call prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,soil_C_total)
      call get_soil_C_total(N_CASA_LAYERS, soil_C_total_H)
      call prescr_get_soilpools(I0,I1,J0,J1,
     &      soil_C_total_H(:,I0:I1,J0:J1), Tpool_ini)
#endif


      ! Carbon pools not ready yet
      ! Tpool_ini = 0.d0
      

      !Compute soil textures for upper 30cm
      do j=j0,j1
        do i=i0,i1
          !if ( fearth(i,j) < .01 ) cycle

          soil_texture(:,i,j) = (
     &            q_ij(i,j,:,1)*dz_ij(i,j,1)
     &         +  q_ij(i,j,:,2)*dz_ij(i,j,2)
     &         +  q_ij(i,j,:,3)*(.3d0 - dz_ij(i,j,1) - dz_ij(i,j,2))
     &         ) / .3d0

cddd          write(887,'(a1,5f10.5)') '1', q_ij(i,j,:,1)
cddd          write(887,'(a1,5f10.5)') '2', q_ij(i,j,:,2)
cddd          write(887,'(a1,5f10.5)') '3', q_ij(i,j,:,3)
cddd          write(887,'(a1,5f10.5)') 'o', soil_texture(:,i,j)
cddd          write(887,'(a1,5f10.5)') 'n', soil_texture1(:,i,j)
cddd          write(887,'(a1,5f10.5)') 'd', soil_texture1(:,i,j)
cddd     &         -soil_texture(:,i,j)
cddd          write(887,'(a1,5f10.5)') 's',sum( soil_texture1(:,i,j) )
        enddo
      enddo

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soil_color, albedodata, soil_texture,
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
     &     im, jm, i0, i1, j0, j1, jday, jyear )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use DOMAIN_DECOMP_ATM, only : GRID
      use geom, only : lat2d
      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata,prescr_get_cropdata
      !use ent_prescr_veg, only: prescr_get_laidata,prescr_veg_albedodata
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, jyear
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      !real*8, dimension(I0:I1,J0:J1) :: cropsdata
      !-----Local---------
      integer hemi(I0:I1,J0:J1)
      integer i,j
      integer year
      integer, save :: year_old = -1
      real*8 :: cropdata_H(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)


! update crops here 
!          year1 = 1965
      if ( crops_yr .ne. 0 ) then
        year = crops_yr
      else
        year = Jyear
      endif

      if( year .ne. year_old ) then
        !call prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropsdata)
        call get_cropdata(year, cropdata_H)
        call ent_prescribe_vegupdate(entcells,
     &       do_giss_lai=.false.,
     &       cropsdata=cropdata_H(I0:I1,J0:J1) )
        year_old = year
      endif

!!!! HACK : trying to update Ent exactly like in ent_prog

            !* Set hemisphere flags.
      where(lat2d(I0:I1,J0:J1) <= 0.)
        hemi(I0:I1,J0:J1) = -1  ! S
      elsewhere
        hemi(I0:I1,J0:J1) = +1  ! N
      end where
 
          call ent_prescribe_vegupdate(entcells,hemi,jday,year,
     &         do_giss_phenology=.true.,
     &         do_giss_albedo=.true.,
     &         do_giss_lai=.true.,
     &         update_crops=.false. )


      return

      call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)

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
