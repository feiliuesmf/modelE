      module ent_mod
!@sum this module contains the interface to Ent dynamic vegetation model
!@+   this is the only module that should be visible to GCM
 
      ! need to know about internal structure of Ent types
      use ent_types, only : entcelltype, patch, cohort, timestruct,
     &     MAX_PATCHES, MAX_COHORTS, ent_config
      use ent_const, only : N_BANDS, N_COVERTYPES, N_DEPTH, N_SOIL_TYPES
     &     , N_BPOOLS, N_PFT, N_CASA_LAYERS,NPOOLS,NLIVE,CARBON,PTRACE  !added last 5 for soil bgc -PK 
      !use ent_prescrveg
      use entcells

      !use prescr_veg ! just for compilation purposes
      implicit none

      private

      !--- public constants ---
      public N_BANDS, N_COVERTYPES, N_DEPTH, N_SOIL_TYPES, N_BPOOLS
      public N_PFT, N_CASA_LAYERS  !added last one -PK

      public entcelltype_public, ent_cell_pack, ent_cell_unpack
      public ent_get_exports, ent_set_forcings
      public ent_cell_construct, ent_cell_destruct, ent_cell_nullify
      public ent_fast_processes,ent_seasonal_update,ent_vegcover_update
      public ent_cell_set !, ent_cell_update
      public ent_prescribe_vegupdate
      public ent_cell_print
      public ent_initialize

      type entcelltype_public
        private
        type(entcelltype), pointer :: entcell
      end type entcelltype_public

      !---- public interfaces ---

      !--- consttructor/destructor ---
      interface ent_cell_construct
        module procedure ent_cell_construct_single
        module procedure ent_cell_construct_array_1d
        module procedure ent_cell_construct_array_2d
      end interface

      interface ent_cell_destruct
        module procedure ent_cell_destruct_single
        module procedure ent_cell_destruct_array_1d
        module procedure ent_cell_destruct_array_2d
      end interface

      interface ent_cell_nullify
        module procedure ent_cell_nullify_single
        module procedure ent_cell_nullify_array_1d
        module procedure ent_cell_nullify_array_2d
      end interface

      interface ent_cell_pack
      module procedure ent_cell_pack
      module procedure ent_cell_pack_2d
      end interface

      !--- passing initial data to ent cells ---
      interface ent_cell_set
        module procedure ent_cell_set_single
        module procedure ent_cell_set_array_1d
        module procedure ent_cell_set_array_2d
      end interface

      !--- passing updated prescribed data to ent cells ---
cddd      interface ent_cell_update
cddd        module procedure ent_cell_update_single
cddd      end interface ent_cell_update

      !--- set forcings / get exports ---
      interface ent_set_forcings
        module procedure ent_set_forcings_single
        module procedure ent_set_forcings_array_1d
        module procedure ent_set_forcings_array_2d
      end interface

      interface ent_get_exports
        module procedure ent_get_exports_single
        module procedure ent_get_exports_array_1d
        module procedure ent_get_exports_array_2d
      end interface

      !--- run model for fast/medium/slow physics ---
      interface ent_fast_processes
        module procedure ent_fast_processes_single
        module procedure ent_fast_processes_array_1d
        module procedure ent_fast_processes_array_2d
      end interface

      interface ent_seasonal_update
        module procedure ent_seasonal_update_single
        module procedure ent_seasonal_update_array_1d
        module procedure ent_seasonal_update_array_2d
      end interface

      interface ent_vegcover_update
        module procedure ent_vegcover_update_single
        module procedure ent_vegcover_update_array_1d
        module procedure ent_vegcover_update_array_2d
      end interface

      interface ent_cell_print
        module procedure ent_cell_print_single
        module procedure ent_cell_print_array_1d
        module procedure ent_cell_print_array_2d
      end interface


!!! do we need 1d and 2d array interfaces for pack/unpack ?

      !---- private interfaces ----
      interface copy_vars
        module procedure copy_vars_single
        module procedure copy_vars_array
        module procedure copy_vars_i_single
        module procedure copy_vars_i_array
      end interface


      !---- global data ----
      type(ent_config) config

      contains

!*************************************************************************
      subroutine ent_initialize(
     &     do_soilresp
     &     )
!@sum initializes Ent module. This subroutine should set all the flags
!@+   and all the variables that are constant during the run.
      logical, optional :: do_soilresp

      ! first set some defaults:
      config%do_soilresp = .false.

      ! now overwrite defaults with explicitly passed values
      if ( present(do_soilresp) ) config%do_soilresp = do_soilresp

      end subroutine ent_initialize

!*************************************************************************
!---- interfaces to run the model one time step --------------------------
      subroutine ent_prescribe_vegupdate(entcell,hemi,jday,year,
     &     update_crops, do_giss_phenology,
     &     laidata, albedodata, cropsdata)
!@sum updates prescribed vegatation parameters. This parameters can
!@+   be passed directly in form of arrays like laidata or one can
!@+   set a flag requesting certain action like do_giss_phenology.
!@+   All arguments except entcell are optional.
      use ent_prescribed_updates, only:  entcell_vegupdate
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer, intent(in), optional, target :: hemi(:,:)
      integer,intent(in), optional :: jday,year
      logical, intent(in), optional :: update_crops
      logical, intent(in), optional :: do_giss_phenology
      real*8, intent(in), optional, target :: laidata(:,:,:)
      real*8, intent(in), optional, target :: albedodata(:,:,:,:)
      real*8, intent(in), optional, target :: cropsdata(:,:)
      !---
      real*8, allocatable :: cropsdata_loc(:,:)
      real*8, pointer :: laidata_1(:), albedodata_1(:,:), cropsdata_1
      integer, pointer :: hemi_1
      integer i, ic, j, jc

!      write(780,*) __FILE__,__LINE__,present(hemi)

      ic = size(entcell,1)
      jc = size(entcell,2)

      if ( present(update_crops) ) then
        if ( update_crops ) then
          if ( .not. present(year) )
     &         call stop_model("ent_prescribe_vegupdate: need year",255)
          allocate ( cropsdata_loc(ic,jc) )
          ! insert call to get_crops here
          ! maybe we can avoid it ? I mean could we always
          ! pass cropdata from the driver?
        endif
      endif

      nullify( laidata_1, albedodata_1, cropsdata_1, hemi_1 )

      do j=1,jc
        do i=1,ic
          ! skip uninitialized cells (no land)
          if ( .not. associated(entcell(i,j)%entcell) ) cycle

          if ( present(laidata) ) laidata_1 => laidata(:,i,j)
          if ( present(albedodata) ) albedodata_1 => albedodata(:,:,i,j)
          if ( present(cropsdata) ) cropsdata_1 => cropsdata(i,j)
          if ( present(hemi) ) hemi_1 => hemi(i,j)

!          write(780,*) __FILE__,__LINE__,hemi_1
          
          call entcell_vegupdate(entcell(i,j)%entcell, hemi_1,
     &         jday, do_giss_phenology, 
     &         laidata_1, albedodata_1, cropsdata_1)
        enddo
      enddo

      end subroutine ent_prescribe_vegupdate
      
!*************************************************************************

      subroutine ent_fast_processes_single(entcell, dt)
      use ent, only : ent_biophysics
      type(entcelltype_public), intent(inout) :: entcell
      real*8, intent(in) :: dt
      !---
      
      call ent_biophysics(dt, entcell%entcell, config)

      end subroutine ent_fast_processes_single


      subroutine ent_fast_processes_array_1d(entcell, dt)
      type(entcelltype_public), intent(inout) :: entcell(:)
      real*8, intent(in) :: dt
      !---
      integer n, nc

      nc = size(entcell)
      do n=1,nc
        call ent_fast_processes_single( entcell(n), dt )
      enddo

      end subroutine ent_fast_processes_array_1d


      subroutine ent_fast_processes_array_2d(entcell, dt )
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      real*8, intent(in) :: dt
      !---
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_fast_processes_single( entcell(i,j),dt )
        enddo
      enddo

      end subroutine ent_fast_processes_array_2d

!*************************************************************************

      subroutine ent_seasonal_update_single(entcell,
     &     dt,
! insert any needed input parameters here
     &     time) !KIM - for phenology:
                 !KIM - consiquently time is added in subroutines *_1d & *_2d 
      use ent, only : ent_integrate_GISS, ent_integrate
!!! it is not clear yet for me how this call will be implemented ...
!@sum this call updates variable that change on a long time scale.
!@+   Right now (before real dynamic vegetation is implemented)
!@+   it should perform prescribed seasonal update of vegatation
!@+   parameters (LAI, root fraction etc.)
!@+   I think extra input parameters needed here should be passed 
!@+   as formal parameters and not be packed into entcell structure.
!@+   It seems that for prescribed variation of vegeatation
!@+   parameters we need only "jday"
!@+   Is it OK from ESMF point of view?
      !use ent_driver, only : ent_update_veg_structure
      type(entcelltype_public), intent(inout) :: entcell !KIM - changed in->inout
      real*8, intent(in) :: time !KIM - for phenology
      real*8, intent(in) :: dt !Time step (s)
      !---
      
      !call ent_update_veg_structure( entcell%entcell, jday )

#ifdef PFT_MODEL_ENT
      call ent_integrate(dt, entcell%entcell, time) 
#else
      call ent_integrate_GISS(entcell%entcell,dt)
#endif

      end subroutine ent_seasonal_update_single


      subroutine ent_seasonal_update_array_1d(entcell, dt,time)
      type(entcelltype_public), intent(inout) :: entcell(:)
      real*8, intent(in) :: time 
      real*8, intent(in) :: dt !Time step (s)
      !---
      integer n, nc

      nc = size(entcell)
      do n=1,nc
        call ent_seasonal_update_single( entcell(n), dt, time)
      enddo

      end subroutine ent_seasonal_update_array_1d


      subroutine ent_seasonal_update_array_2d(entcell,dt,time)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      real*8, intent(in) :: time 
      real*8, intent(in) :: dt !Time step (s)
!      integer, intent(in) :: jday
      !---
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_seasonal_update_single( entcell(i,j), dt, time ) 
        enddo
      enddo

      end subroutine ent_seasonal_update_array_2d

!*************************************************************************

      subroutine ent_vegcover_update_single(entcell,
     $     jday,
     $     jyear
     $     )
      type(entcelltype_public), intent(in) :: entcell
      integer, intent(in) :: jday, jyear
       

      ! no code for vegcover_update yet ...

      end subroutine ent_vegcover_update_single


      subroutine ent_vegcover_update_array_1d(entcell, jday, jyear)
      type(entcelltype_public), intent(inout) :: entcell(:)
      integer, intent(in) :: jday, jyear
      !---
      integer n, nc

      nc = size(entcell)
      do n=1,nc
        call ent_vegcover_update_single( entcell(n), jday, jyear )
      enddo

      end subroutine ent_vegcover_update_array_1d


      subroutine ent_vegcover_update_array_2d(entcell, jday, jyear)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer, intent(in) :: jday, jyear
      !---      
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_vegcover_update_single( entcell(i,j), jday, jyear )
        enddo
      enddo

      end subroutine ent_vegcover_update_array_2d

!---- END interfaces to run the model one time step ----
!*************************************************************************


!---- Constructor / Destructor -------------------------------------------

      subroutine ent_cell_construct_single(entcell)
      use entcells, only : entcell_construct
      type(entcelltype_public), intent(inout) :: entcell

cddd      allocate( entcell%entcell )
cddd
cddd      ! allocate internal arrays
cddd      allocate( entcell%entcell%fracroot(N_DEPTH) )
cddd      allocate( entcell%entcell%betadl(N_DEPTH) )
cddd      allocate( entcell%entcell%Soilmoist(N_DEPTH) )
cddd      allocate( entcell%entcell%Soilmp(N_DEPTH) )
cddd      allocate( entcell%entcell%fice(N_DEPTH) )
cddd
cddd      ! don't allocate patches, just set all pointers to NULL
cddd      nullify( entcell%entcell%youngest )
cddd      nullify( entcell%entcell%oldest   )
cddd      ! for now set all values o zero or defaults
cddd      call zero_entcell(entcell%entcell)

      !print *,"ent_cell_constr"
      call entcell_construct(entcell%entcell)
      !call entcell_print(6, entcell%entcell)

      end subroutine ent_cell_construct_single


      subroutine ent_cell_construct_array_1d(entcell)
      type(entcelltype_public), intent(inout) :: entcell(:)
      integer n, nc

      nc = size(entcell)

      do n=1,nc
        call ent_cell_construct_single( entcell(n) )
      enddo

      end subroutine ent_cell_construct_array_1d


      subroutine ent_cell_construct_array_2d(entcell)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer i, ic, j, jc

      !print *,"ent_cell_construct_array_2d:"

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          !print *,"i,j=",i,j
          call ent_cell_construct_single( entcell(i,j) )
        enddo
      enddo

      end subroutine ent_cell_construct_array_2d

!*************************************************************************

      subroutine ent_cell_destruct_single(entcell)
      use entcells, only : entcell_destruct
      type(entcelltype_public), intent(inout) :: entcell

      call entcell_destruct( entcell%entcell )

      end subroutine ent_cell_destruct_single


      subroutine ent_cell_destruct_array_1d(entcell)
      type(entcelltype_public), intent(inout) :: entcell(:)
      integer n, nc

      nc = size(entcell)

      do n=1,nc
        call ent_cell_destruct_single( entcell(n) )
      enddo

      end subroutine ent_cell_destruct_array_1d


      subroutine ent_cell_destruct_array_2d(entcell)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_cell_destruct_single( entcell(i,j) )
        enddo
      enddo

      end subroutine ent_cell_destruct_array_2d


      subroutine ent_cell_nullify_single(entcell)
      use entcells, only : entcell_destruct
      type(entcelltype_public), intent(inout) :: entcell

      nullify( entcell%entcell )

      end subroutine ent_cell_nullify_single


      subroutine ent_cell_nullify_array_1d(entcell)
      type(entcelltype_public), intent(inout) :: entcell(:)
      integer n, nc

      nc = size(entcell)

      do n=1,nc
        nullify( entcell(n)%entcell )
      enddo

      end subroutine ent_cell_nullify_array_1d


      subroutine ent_cell_nullify_array_2d(entcell)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          nullify( entcell(i,j)%entcell )
        enddo
      enddo

      end subroutine ent_cell_nullify_array_2d


!---- END of  Constructor / Destructor -----
!*************************************************************************

      
      subroutine ent_cell_set_single(entcell,
     &     veg_fraction,
     &     pft_population_density,
     &     leaf_area_index,
     &     pft_heights,
     &     pft_dbh,
     &     pft_crad,
     &     pft_cpool,
     &     pft_nmdata,
     &     pft_froots,
     &     pft_soil_type,
     &     vegalbedo,
     &     soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini
     &     )
      type(entcelltype_public), intent(out) :: entcell
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     veg_fraction,
     &     leaf_area_index
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     pft_heights,
     &     pft_dbh,
     &     pft_crad,
     &     pft_nmdata,
     &     pft_population_density
      real*8, dimension(:,:) :: pft_cpool !Carbon pools in individuals
      real*8, dimension(:,:)  :: pft_froots
      integer, dimension(:)  :: pft_soil_type
      real*8, dimension(:,:)  ::  vegalbedo ! dim=N_BANDS,N_COVERTYPES
      real*8, dimension(:)  ::  soil_texture ! dim=N_SOIL_TYPES
      real*8 :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini 

      call init_simple_entcell( entcell%entcell,
     &     veg_fraction, pft_population_density, leaf_area_index,
     &     pft_heights, pft_dbh,pft_crad,pft_cpool, pft_nmdata, 
     &     pft_froots, pft_soil_type, vegalbedo, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
      
      end subroutine ent_cell_set_single


      subroutine ent_cell_set_array_1d(entcell,
     &     veg_fraction,
     &     pft_population_density,
     &     leaf_area_index,
     &     pft_heights,
     &     pft_dbh,
     &     pft_crad,
     &     pft_cpool,
     &     pft_nmdata,
     &     pft_froots,
     &     pft_soil_type,
     &     vegalbedo,
     &     soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
      type(entcelltype_public), intent(out) :: entcell(:)
      real*8, dimension(:,:)  ::   ! dim=N_COVERTYPES, n
     &     veg_fraction,
     &     leaf_area_index
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     pft_heights,
     &     pft_dbh,
     &     pft_crad,
     &     pft_nmdata,
     &     pft_population_density
      real*8, dimension(:,:,:) :: pft_cpool !Carbon pools in individuals
      real*8, dimension(:,:)  :: pft_froots
      integer, dimension(:)  :: pft_soil_type
      real*8, dimension(:,:,:)  ::  vegalbedo ! dim=N_BANDS,N_COVERTYPES
      real*8, dimension(:,:)  ::    soil_texture ! dim=N_SOIL_TYPES
      real*8, dimension(:) :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini
      !---
      integer n, nc

      nc = size(entcell, 1)

      do n=1,nc
        call init_simple_entcell( entcell(n)%entcell,
     &       veg_fraction(:,n), pft_population_density,
     &       leaf_area_index(:,n), pft_heights, pft_dbh, pft_crad, 
     &       pft_cpool(:,:,n), pft_nmdata, pft_froots, 
     &       pft_soil_type, vegalbedo(:,:,n), soil_texture(:,n),
     &       Ci_ini(n), CNC_ini(n), Tcan_ini(n), Qf_ini(n))
      enddo
      
      end subroutine ent_cell_set_array_1d


      subroutine ent_cell_set_array_2d(entcell,
     &     veg_fraction,
     &     pft_population_density,
     &     leaf_area_index,
     &     pft_heights,
     &     pft_dbh,
     &     pft_crad,
     &     pft_cpool,
     &     pft_nmdata,
     &     pft_froots,
     &     pft_soil_type,
     &     vegalbedo,
     &     soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
      type(entcelltype_public), intent(out) :: entcell(:,:)
      real*8, dimension(:,:,:)  ::   ! dim=N_COVERTYPES, n
     &     veg_fraction,
     &     leaf_area_index
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     pft_heights,
     &     pft_dbh,
     &     pft_crad,
     &     pft_nmdata,
     &     pft_population_density
      real*8, dimension(:,:,:,:) :: pft_cpool !Carbon pools in individuals
      real*8, dimension(:,:)  :: pft_froots
      integer, dimension(:)  :: pft_soil_type
      real*8, dimension(:,:,:,:)  ::  vegalbedo ! dim=N_COVERTYPES, n
      real*8, dimension(:,:,:)  ::  soil_texture ! dim=N_SOIL_TYPES
      real*8, dimension(:,:) :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini
      !---
      integer i, j, ic, jc

      ic = size(entcell, 1)
      jc = size(entcell, 2)

      do j=1,jc
        do i=1,ic
          !print *,"ent_cell_set_array_2d i,j=",i,j
          if ( .not. associated(entcell(i,j)%entcell) ) cycle
!      if ( .not. associated(ecp) ) 
!     &      call stop_model("init_simple_entcell 1",255)
          !call entcell_print(6,entcell(i,j)%entcell)

          call init_simple_entcell( entcell(i,j)%entcell,
     &         veg_fraction(:,i,j),pft_population_density,
     &         leaf_area_index(:,i,j),
     &         pft_heights,pft_dbh,pft_crad,pft_cpool(:,:,i,j),
     &         pft_nmdata,
     &         pft_froots,
     &         pft_soil_type,vegalbedo(:,:,i,j), soil_texture(:,i,j),
     &         Ci_ini(i,j), CNC_ini(i,j), Tcan_ini(i,j), Qf_ini(i,j))
        enddo
      enddo
      
      end subroutine ent_cell_set_array_2d

!*************************************************************************

cddd      subroutine ent_cell_update_single(entcell,
cddd     &     pft_population_density,
cddd     &     pft_leaf_area_index,
cddd     &     pft_heights,
cddd     &     pft_dbh,
cddd     &     pft_crad,
cddd     &     pft_cpool,
cddd     &     pft_nmdata,
cddd     &     pft_froots,
cddd     &     pft_vegalbedo,
cddd     &     heat_capacity
cddd     &     )
cddd      use ent_prescribed_updates, only : entcell_update_lai,
cddd     &     entcell_update_albedo, entcell_update_shc
cddd      type(entcelltype_public), intent(out) :: entcell
cddd      real*8, dimension(:), optional  ::   ! dim=N_PFT
cddd     &     pft_leaf_area_index,
cddd     &     pft_heights,
cddd     &     pft_dbh,
cddd     &     pft_crad,
cddd     &     pft_nmdata,
cddd     &     pft_population_density
cddd      real*8, dimension(:,:), optional :: pft_cpool !Carbon pools in individuals
cddd      real*8, dimension(:,:), optional :: pft_froots
cddd      real*8, dimension(:,:), optional ::  pft_vegalbedo ! dim=N_BANDS,N_PFT
cddd      ! the following is needed for a hack to set GISS canopy heat capacity
cddd      real*8, optional :: heat_capacity
cddd
cddd      ! if cell is not initialized (no land) do nothing
cddd      if ( .not. associated(entcell%entcell) ) return
cddd
cddd      if ( present(pft_leaf_area_index) ) then
cddd        call entcell_update_lai( entcell%entcell, pft_leaf_area_index )
cddd      endif
cddd
cddd      if ( present(pft_vegalbedo) ) then
cddd        call entcell_update_albedo( entcell%entcell, pft_vegalbedo )
cddd      endif
cddd
cddd      !!! hack, should be replaced with computation of heat capacity
cddd      !!! inside ent
cddd      if ( present(heat_capacity) ) then
cddd        entcell%entcell%heat_capacity = heat_capacity
cddd      endif
cddd
cddd      !!!! subroutine litter(dtsec, pp) is not called here !!!
cddd      !!!! should be called somewhere else
cddd     
cddd
cddd      if ( present(pft_population_density)
cddd     &     .or. present(pft_heights)
cddd     &     .or. present(pft_dbh)
cddd     &     .or. present(pft_crad)
cddd     &     .or. present(pft_cpool)
cddd     &     .or. present(pft_nmdata)
cddd     &     .or. present(pft_froots)
cddd     &      ) then
cddd        call stop_model("ent_cell_update: var not supported yet", 255)
cddd      endif
cddd      
cddd      ! just in case, maybe we don't need to summarize if all updates
cddd      ! also update summarized values
cddd      call summarize_entcell(entcell%entcell)
cddd
cddd      ! make sure heat capacity of canopy is up-to-date
cddd      if ( present(pft_leaf_area_index) ) then ! i.e. LAI changed
cddd        call entcell_update_shc( entcell%entcell )
cddd      endif
cddd
cddd
cddd      end subroutine ent_cell_update_single


!*************************************************************************

      subroutine ent_cell_pack_2d(dbuf, entcell)
!@sum allocate single linear arrays dbuf and pack contents of
!@+   entcells(i,j) into it
      real*8, pointer :: dbuf(:)
      type(entcelltype_public), intent(in) :: entcell(:,:)
      !---
      type real8_ptr
        real*8, pointer :: ptr(:)
      end type real8_ptr
      type(real8_ptr), dimension(:,:), allocatable :: buf2d
      integer i, j, ic, jc, dc, dcc, lsize

      ic = size(entcell, 1)
      jc = size(entcell, 2)

      allocate( buf2d(ic,jc) )
      
      dc = 0
      do j=1,jc
        do i=1,ic
          !print *,"ent_cell_pack_2d i,j=",i,j
          nullify( buf2d(i,j)%ptr )
          if ( .not. associated(entcell(i,j)%entcell) ) cycle

          call ent_cell_pack(buf2d(i,j)%ptr, entcell(i,j))
          dc = dc + size(buf2d(i,j)%ptr, 1) + 3 ! 3 =: i,j,size
          
        enddo
      enddo

      allocate( dbuf(dc) )
      dcc = 1
      do j=1,jc
        do i=1,ic
          if ( .not. associated( buf2d(i,j)%ptr ) ) cycle

          !print *,"ent_cell_pack_2d i,j,dcc,lsize=",i,j,dcc,lsize
          !print *,buf2d(i,j)%ptr(1:lsize)
          lsize = size(buf2d(i,j)%ptr, 1)
          dbuf(dcc) = i; dcc = dcc+1
          dbuf(dcc) = j; dcc = dcc+1
          dbuf(dcc) = lsize; dcc = dcc+1
          dbuf(dcc:dcc+lsize-1) = buf2d(i,j)%ptr(1:lsize)
          dcc = dcc+lsize
          deallocate( buf2d(i,j)%ptr )

        enddo
      enddo

      if ( dcc-1 .ne. dc ) call stop_model("ent_cell_pack_2d: dcc",255)
     
      deallocate( buf2d )

      end subroutine ent_cell_pack_2d


      subroutine ent_cell_pack(dbuf, entcell)
!@sum allocate single linear arrays dbuf and pack contents of
!@+   entcell into it
      real*8, pointer :: dbuf(:)
      type(entcelltype_public), intent(in) :: entcell ! pointer ?
      !---
      type(patch), pointer :: p  !@var p current patch
      type(cohort), pointer :: c !@var current cohort
      integer :: np              !@var np number of patches in the cell
      integer :: nc(MAX_PATCHES) !@var nc number of cohorts in the patch
      integer :: dc, ndbuf, nn
      real*8, pointer :: NULL(:) !@var NULL dummy pointer

      nullify(NULL)

      ! return "-1" for not associated cells
      if ( .not. associated(entcell%entcell) ) then
        allocate( dbuf(1) )
        dbuf(1) = -1.d0;
        return
      endif

      ! first compute number of patches and cohorts in the cell
      ! this actually can be save in the cell structure 
      ! for optimization ...
      ! also count the number of real*8 values to be saved
      ndbuf = 0
      np = 0
      p => entcell%entcell%oldest      
      do while ( associated(p) )
        np = np + 1
        if ( np > MAX_PATCHES )
     &       call stop_model("ent_cell_pack: too many patches",255)
        call copy_patch_vars(NULL, nn, p, 0); ndbuf = ndbuf + nn
        nc(np) = 0
        c => p%tallest
        do while ( associated(c) )
          nc(np) = nc(np) + 1
          if ( nc(np) > MAX_COHORTS )
     &         call stop_model("ent_cell_pack: too many cohorts",255)
          !save cohort
          !dbuf(dc) = c%_any_value_ ; dc = dc + 1
          call copy_cohort_vars(NULL, nn, c, 0); ndbuf = ndbuf + nn
          c => c%shorter
        enddo
        p => p%younger
      enddo

      allocate( dbuf(ndbuf+1+np) ) !i.e. num reals + num int's
      dc = 0
      dbuf(dc+1) = real( np, kind(0d0) );               dc = dc + 1
      !print *,"pack ", np, dbuf(1)
      dbuf(dc+1:dc+np) = real( nc(1:np), kind(0d0) ); dc = dc + np
      !print *,"pack1 ", nc(1:np), dbuf(2:dc) 

      ! now do the real saving
      ! no need to count patches and cohorts again, but leaving it here
      ! for a while for debugging
      ! save cell vars here (if there are any...), i.e. 
      ! call copy_cell_vars(dbuf, nn, p, -1);
      np = 0
      p => entcell%entcell%oldest      
      do while ( associated(p) )
        np = np + 1
        if ( np > MAX_PATCHES )
     &       call stop_model("ent_cell_pack: too many patches",255)
        !save patch
        call copy_patch_vars(dbuf(dc+1:), nn, p, -1); dc = dc + nn
        nc(np) = 0
        c => p%tallest
        do while ( associated(c) )
          nc(np) = nc(np) + 1
          if ( nc(np) > MAX_COHORTS )
     &         call stop_model("ent_cell_pack: too many cohorts",255)
          !save cohort
          call copy_cohort_vars(dbuf(dc+1:), nn, c, -1); dc = dc + nn
         c => c%shorter
        enddo
        p => p%younger
      enddo

      if ( dbuf(1) .ne. np ) then
        print *,"GGGGGGGGGGG", np, nc(1:np), "XX", dbuf
      endif

      end subroutine ent_cell_pack

!*************************************************************************

      subroutine ent_cell_unpack(dbuf, entcell)
! this program is not finished yet: have to assign all the pointers
      use cohorts, only : cohort_construct
      use patches, only : patch_construct
      real*8, intent(inout) :: dbuf(0:)
      type(entcelltype_public), intent(out) :: entcell ! pointer ?
      !---
      type(patch), pointer :: p, pprev  !@var p current patch
      type(cohort), pointer :: c, cprev !@var current cohort
      integer :: np              !@var np number of patches in the cell
      integer, allocatable :: nc(:) !@var nc number of cohorts in the patch
      integer dc, nn
      integer i, j
      integer npdebug, ncdebug ! these are for debuging

      dc = 0

      ! doesn't seem that we need to restore anything for the cell
      np = nint( dbuf(dc) ); dc = dc + 1
      if ( np == -1 ) return  ! no data for this cell

      allocate( nc(np) )
      nc(1:np) = nint( dbuf(dc:dc+np-1) ); dc = dc + np

      if ( np <= 0 ) return  ! nothing to restore...

      nullify( pprev )
      do i=1,np
        !allocate( p )
        call patch_construct(p, entcell%entcell, 0.d0, -1)
        call copy_patch_vars(dbuf(dc:), nn, p, 1); dc = dc + nn
        p%older => pprev
        nullify( cprev)
        do j=1,nc(i)
          !allocate( c )
          call cohort_construct(c, p)
          call copy_cohort_vars(dbuf(dc:), nn, c, 1); dc = dc + nn
          c%taller => cprev
          cprev => c
        enddo
        p%shortest => cprev
        pprev => p
      enddo
      entcell%entcell%youngest => p

      ! now restore pointer lists in opposite direction
      npdebug = 0
      nullify( pprev )
      p => entcell%entcell%youngest
      do while ( associated(p) )
        p%younger => pprev
        p%cellptr => entcell%entcell
        npdebug = npdebug + 1
        if ( npdebug > np )
     &       call stop_model("ent_cell_unpack: broken struct: np",255)
        ncdebug = 0
        nullify( cprev)
        c => p%shortest
        do while ( associated(c) )
          c%shorter => cprev
          c%pptr => p
          c%cellptr => entcell%entcell
          ncdebug = ncdebug + 1
          if ( ncdebug > nc(np-npdebug+1) )
     &         call stop_model("ent_cell_unpack: broken struct: nc",255)
          cprev => c
          c => c%taller
        enddo
        p%tallest => cprev
        pprev => p
        p => p%older
      enddo
      entcell%entcell%oldest => pprev

      deallocate( nc )

      call summarize_entcell(entcell%entcell)

      end subroutine ent_cell_unpack

!*************************************************************************

      subroutine copy_vars_single( buf, n, var, flag )
!@copy variable to/from buffer
!@+   !!! may need to write similar for arrays and create an interface
!@+   !!! in that case "n" will have non-triial value
      real*8, intent(inout) :: buf(:)
      integer, intent(out) :: n
      real*8, intent(inout):: var
!@var flag defines the actual action:
!@+     -1 copy from var to buffer
!@+      1 copy from buffer to var
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      
      n = 1
      if ( flag == 0 ) return

      if ( flag == -1 ) then
        buf(1) = var
      else if ( flag == 1 ) then
        var = buf(1)
      else
        call stop_model("ent_mod:copy_vars: flag .ne. 0,-1,1",255)
      endif

      end subroutine copy_vars_single

      subroutine copy_vars_array( buf, n, var, flag )
!@copy variable to/from buffer
!@+   !!! may need to write similar for arrays and create an interface
!@+   !!! in that case "n" will have non-triial value
      real*8, intent(inout) :: buf(:)
      integer, intent(out) :: n
      real*8, intent(inout):: var(:)
!@var flag defines the actual action:
!@+     -1 copy from var to buffer
!@+      1 copy from buffer to var
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      
      n = size(var)
      if ( flag == 0 ) return

      if ( flag == -1 ) then
        buf(1:n) = var(1:n)
      else if ( flag == 1 ) then
        var(1:n) = buf(1:n)
      else
        call stop_model("ent_mod:copy_vars: flag .ne. 0,-1,1",255)
      endif

      end subroutine copy_vars_array

      subroutine copy_vars_i_single( buf, n, var, flag )
!@copy variable to/from buffer
!@+   !!! may need to write similar for arrays and create an interface
!@+   !!! in that case "n" will have non-triial value
      real*8, intent(inout) :: buf(:)
      integer, intent(out) :: n
      integer, intent(inout):: var
!@var flag defines the actual action:
!@+     -1 copy from var to buffer
!@+      1 copy from buffer to var
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      
      n = 1
      if ( flag == 0 ) return

      if ( flag == -1 ) then
        buf(1) = real( var, kind(0d0) )
      else if ( flag == 1 ) then
        var = nint( buf(1) )
      else
        call stop_model("ent_mod:copy_vars: flag .ne. 0,-1,1",255)
      endif

      end subroutine copy_vars_i_single

      subroutine copy_vars_i_array( buf, n, var, flag )
!@copy variable to/from buffer
!@+   !!! may need to write similar for arrays and create an interface
!@+   !!! in that case "n" will have non-triial value
      real*8, intent(inout) :: buf(:)
      integer, intent(out) :: n
      integer, intent(inout):: var(:)
!@var flag defines the actual action:
!@+     -1 copy from var to buffer
!@+      1 copy from buffer to var
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      
      n = size(var)
      if ( flag == 0 ) return

      if ( flag == -1 ) then
        buf(1:n) = real( var(1:n), kind(0d0) )
      else if ( flag == 1 ) then
        var(1:n) = nint( buf(1:n) )
      else
        call stop_model("ent_mod:copy_vars: flag .ne. 0,-1,1",255)
      endif

      end subroutine copy_vars_i_array

!*************************************************************************

!**************************************************************
!   the following two functions are all that user has to modify
!   when the list of i/o variables is changed
!   I wrote it in such a complicated way so that the list of
!   i/o variable appears only once (and is used both for input
!   and output). This prevents possible confusion due to
!   non-synchronized input and output lists.

!   i didn't include any i/o sub for cell since it looks like 
!   patch will not have any i/o vars

      subroutine copy_patch_vars(buf, n, p, flag)
      real*8, intent(inout) :: buf(0:)
      integer, intent(out) :: n
      type(patch), intent(inout):: p
!@var flag defines the actual action:
!@+     -1 copy from patch to buffer
!@+      1 copy from buffer to patch
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      integer dc, nn, i

      dc = 0

      ! include all patch variables that need i/o
      call copy_vars( buf(dc:), nn,  p%age,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  p%area, flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  p%Ci,   flag ); dc = dc + nn
      do i=1,N_CASA_LAYERS      !need b/c Tpool now rank 3  -PK  
       call copy_vars( buf(dc:), nn,  p%Tpool(1,:,i),flag );dc = dc + nn
       call copy_vars( buf(dc:), nn,  p%Tpool(2,:,i),flag );dc = dc + nn
      end do
      ! not sure about the following, probably can be restored from 
      ! other data...
      call copy_vars( buf(dc:), nn,  p%soil_type, flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  p%GCANOPY, flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  p%albedo, flag ); dc = dc + nn

      n = dc

      end subroutine copy_patch_vars


      subroutine copy_cohort_vars(buf, n, c, flag)
      real*8, intent(inout) :: buf(0:)
      integer, intent(out) :: n
      type(cohort), intent(inout):: c
!@var flag defines the actual action:
!@+     -1 copy from patch to buffer
!@+      1 copy from buffer to patch
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      integer dc, nn

      dc = 0

      ! include all cohort variables that need i/o
      call copy_vars( buf(dc:), nn,  c%pft,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%n,    flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%nm,   flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%lai,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%h,    flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%dbh,  flag ); dc = dc + nn
!      call copy_vars( buf(dc:), nn,  c%_any_var2_, flag ); dc = dc + nn
      ! data for Tpool, do we need these?
      call copy_vars( buf(dc:), nn,  c%C_fol,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%C_froot,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%C_hw,  flag ); dc = dc + nn
      ! I guess fracroot is also needed ...
      call copy_vars( buf(dc:), nn,  c%fracroot,  flag ); dc = dc + nn

      ! added new data to restore checkpoint after sumcohort was removed...
      call copy_vars( buf(dc:), nn,  c%Ci,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%gcanopy,  flag ); dc = dc + nn

      n = dc

      end subroutine copy_cohort_vars

!******************************************************************

      subroutine ent_set_forcings_single( entcell,
     &     air_temperature, !KIM - for phenology
     &     canopy_temperature,
     &     canopy_air_humidity,      
     &     surf_pressure,            
     &     surf_CO2,                 
     &     heat_transfer_coef,       
     &     wind_speed,               
     &     total_visible_rad,
     &     direct_visible_rad,
     &     cos_solar_zenith_angle,
!     &     soil_temp30cm,       !added soil T, volum moist (avg top 30 cm) -PK 6/28/06
!     &     soil_moist30cm,
     &     soil_temp,   !now want explicit depth-structured soiltemp,soilmoist--see ent_types -PK 7/07
     &     soil_moist,
!     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
     &     ) ! need to pass Ci, Qf ??
      type(entcelltype_public), intent(out) :: entcell
      ! forcings probably should not be optional ...
      real*8, intent(in) ::
     &     air_temperature, !KIM - for phenology
     &     canopy_temperature,
     &     canopy_air_humidity,
     &     surf_pressure,
     &     surf_CO2,
     &     heat_transfer_coef,
     &     wind_speed,
     &     total_visible_rad,
     &     direct_visible_rad,
     &     cos_solar_zenith_angle
!     &     soil_temp30cm,        
!     &     soil_moist30cm
      real*8, dimension(:), intent(in) ::
     &     soil_temp,
     &     soil_moist,
!     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
      !----------
      integer n

      entcell%entcell%TairC = air_temperature !KIM - for phenology
      entcell%entcell%TcanopyC = canopy_temperature
      entcell%entcell%Qf = canopy_air_humidity
      entcell%entcell%P_mbar = surf_pressure
      entcell%entcell%Ca = surf_CO2
      entcell%entcell%Ch = heat_transfer_coef
      entcell%entcell%U = wind_speed
      entcell%entcell%IPARdif = total_visible_rad-direct_visible_rad
      entcell%entcell%IPARdir = direct_visible_rad
      entcell%entcell%CosZen = cos_solar_zenith_angle
!      entcell%entcell%Soiltemp = soil_temp30cm
!      entcell%entcell%Soilmoist = soil_moist30cm
      do n=1,N_CASA_LAYERS
        entcell%entcell%Soiltemp(n) = soil_temp(n)
        entcell%entcell%Soilmoist(n) = soil_moist(n)
      end do
      do n=1,N_DEPTH
        entcell%entcell%Soilmp(n) = soil_matric_pot(n)
        entcell%entcell%fice(n) = soil_ice_fraction(n)
      enddo

      end subroutine ent_set_forcings_single


      subroutine ent_set_forcings_array_1d( entcell,
     &     air_temperature, !KIM - for phenology
     &     canopy_temperature,
     &     canopy_air_humidity,      
     &     surf_pressure,            
     &     surf_CO2,                 
     &     heat_transfer_coef,       
     &     wind_speed,               
     &     total_visible_rad,
     &     direct_visible_rad,
     &     cos_solar_zenith_angle,
!     &     soil_temp30cm,
!     &     soil_moist30cm,
     &     soil_temp,
     &     soil_moist,
!     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
     &     ) ! need to pass Ci, Qf ??
      type(entcelltype_public), dimension(:), intent(out) :: entcell
      ! forcings probably should not be optional ...
      real*8, dimension(:)  ::
     &     air_temperature, !KIM - for phenology
     &     canopy_temperature,
     &     canopy_air_humidity,
     &     surf_pressure,
     &     surf_CO2,
     &     heat_transfer_coef,
     &     wind_speed,
     &     total_visible_rad,
     &     direct_visible_rad,
     &     cos_solar_zenith_angle
!     &     soil_temp30cm,       
!     &     soil_moist30cm
      real*8, dimension(:,:), intent(in) ::
     &     soil_temp,
     &     soil_moist,
!     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
      !----------
      integer n
      integer i, ic

      ic = size(entcell, 1)

      do i=1,ic
        entcell(i)%entcell%TairC = air_temperature(i) !KIM - for phenology
        entcell(i)%entcell%TcanopyC = canopy_temperature(i)
        entcell(i)%entcell%Qf = canopy_air_humidity(i)
        entcell(i)%entcell%P_mbar = surf_pressure(i)
        entcell(i)%entcell%Ca = surf_CO2(i)
        entcell(i)%entcell%Ch = heat_transfer_coef(i)
        entcell(i)%entcell%U = wind_speed(i)
        entcell(i)%entcell%IPARdif = total_visible_rad(i)-
     &       direct_visible_rad(i)
        entcell(i)%entcell%IPARdir = direct_visible_rad(i)
        entcell(i)%entcell%CosZen = cos_solar_zenith_angle(i)
!        entcell(i)%entcell%Soiltemp = soil_temp30cm(i)
!        entcell(i)%entcell%Soilmoist = soil_moist30cm(i)
        do n=1,N_CASA_LAYERS
          entcell(i)%entcell%Soiltemp(n) = soil_temp(n,i)
          entcell(i)%entcell%Soilmoist(n) = soil_moist(n,i)
        end do
        do n=1,N_DEPTH
          entcell(i)%entcell%Soilmp(n) = soil_matric_pot(n,i)
          entcell(i)%entcell%fice(n) = soil_ice_fraction(n,i)
        enddo
      enddo

      end subroutine ent_set_forcings_array_1d


      subroutine ent_set_forcings_array_2d( entcell,
     &     air_temperature, !KIM - for phenology
     &     canopy_temperature,
     &     canopy_air_humidity,      
     &     surf_pressure,            
     &     surf_CO2,                 
     &     heat_transfer_coef,       
     &     wind_speed,               
     &     total_visible_rad,
     &     direct_visible_rad,
     &     cos_solar_zenith_angle,
!     &     soil_temp30cm,
!     &     soil_moist30cm,
     &     soil_temp,
     &     soil_moist,
!     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
     &     ) ! need to pass Ci, Qf ??
      type(entcelltype_public), dimension(:,:), intent(out) :: entcell
      ! forcings probably should not be optional ...
      real*8, dimension(:,:)  ::
     &     air_temperature, !KIM - for phenology
     &     canopy_temperature,
     &     canopy_air_humidity,
     &     surf_pressure,
     &     surf_CO2,
     &     heat_transfer_coef,
     &     wind_speed,
     &     total_visible_rad,
     &     direct_visible_rad,
     &     cos_solar_zenith_angle
!     &     soil_temp30cm,
!     &     soil_moist30cm
      real*8, dimension(:,:,:), intent(in) ::
     &     soil_temp,
     &     soil_moist,
!     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
      !----------
      integer n
      integer i, j, ic, jc

      ic = size(entcell, 1)
      jc = size(entcell, 2)
      
      do j=1,jc
        do i=1,ic
          entcell(i,j)%entcell%TairC = air_temperature(i,j) !KIM - for phenoloygy
          entcell(i,j)%entcell%TcanopyC = canopy_temperature(i,j)
          entcell(i,j)%entcell%Qf = canopy_air_humidity(i,j)
          entcell(i,j)%entcell%P_mbar = surf_pressure(i,j)
          entcell(i,j)%entcell%Ca = surf_CO2(i,j)
          entcell(i,j)%entcell%Ch = heat_transfer_coef(i,j)
          entcell(i,j)%entcell%U = wind_speed(i,j)
          entcell(i,j)%entcell%IPARdif = total_visible_rad(i,j)-
     &       direct_visible_rad(i,j)
          entcell(i,j)%entcell%IPARdir = direct_visible_rad(i,j)
          entcell(i,j)%entcell%CosZen = cos_solar_zenith_angle(i,j)
!          entcell(i,j)%entcell%Soiltemp = soil_temp30cm(i,j)
!          entcell(i,j)%entcell%Soilmoist = soil_moist30cm(i,j)
          do n=1,N_CASA_LAYERS
            entcell(i,j)%entcell%Soiltemp(n) = soil_temp(n,i,j)
            entcell(i,j)%entcell%Soilmoist(n) = soil_moist(n,i,j)
          end do
          do n=1,N_DEPTH
            entcell(i,j)%entcell%Soilmp(n) = soil_matric_pot(n,i,j)
            entcell(i,j)%entcell%fice(n) = soil_ice_fraction(n,i,j)
          enddo
        enddo
      enddo

      end subroutine ent_set_forcings_array_2d

!******************************************************************

      subroutine ent_get_exports_single( entcell,
     &     canopy_conductance,
     &     beta_soil_layers,
     &     shortwave_transmit,
     &     foliage_CO2,
     &     foliage_humidity,
     &     canopy_gpp,
     &     roughness_length,
     &     flux_CO2,
     &     albedo,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     vegetation_fractions,
           !added next 2 -PK 7/07
     &     soilresp,
     &     soilcpools
     &     )
      type(entcelltype_public), intent(in) :: entcell
      real*8, optional, intent(out) ::
     &     canopy_conductance,
     &     shortwave_transmit,
     &     foliage_CO2,
     &     foliage_humidity,
     &     canopy_gpp,
     &     roughness_length,
     &     flux_CO2,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     soilresp  
      real*8, dimension(:), optional, intent(out) ::
     &     beta_soil_layers,
     &     albedo,
     &     vegetation_fractions
      real*8, dimension(:,:,:), optional, intent(out) ::
     &     soilcpools 
      !----------
      integer n,p,ii

      if ( present(canopy_conductance) )
     &     canopy_conductance = entcell%entcell%GCANOPY

      if ( present(shortwave_transmit) )
     &     shortwave_transmit = entcell%entcell%TRANS_SW

      if ( present(foliage_CO2) )
     &     foliage_CO2 = entcell%entcell%Ci

      if ( present(foliage_humidity) )
     &     foliage_humidity = entcell%entcell%Qf

      if ( present(canopy_gpp) )
     &     canopy_gpp = entcell%entcell%GPP

      if ( present(roughness_length) )
     &     roughness_length = entcell%entcell%z0

      if ( present(flux_CO2) )
     &     flux_CO2 = entcell%entcell%CO2flux

      if ( present(canopy_max_H2O) )
     &     canopy_max_H2O = entcell%entcell%LAI * .0001d0 !!! GISS setting

      if ( present(canopy_heat_capacity) ) then
        !aa=ala(1,i0,j0)
        canopy_heat_capacity=entcell%entcell%heat_capacity
        !call stop_model("not implemmented yet",255)
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        fraction_of_vegetated_soil = entcell%entcell%fv
      endif

      if ( present(soilresp) )  !PK
     &     soilresp = entcell%entcell%Soil_resp

      if ( present(beta_soil_layers) ) then
        do n=1,N_DEPTH
          beta_soil_layers(n) = entcell%entcell%betadl(n)
        enddo
      endif

      if ( present(albedo) ) then
        do n=1,N_BANDS
          albedo(n) = entcell%entcell%albedo(n)
        enddo
      endif

      if ( present(vegetation_fractions) ) then
        !do n=1,N_COVERTYPES
        ! extract those here ?
        !vegetation_fractions(:) = vdata(i,j,:)
        !  call stop_model("not implemmented yet",255)
        !enddo
        call entcell_extract_pfts(entcell%entcell, vegetation_fractions)
      endif

      if ( present(soilcpools) ) then  !PK
        do n=1,N_CASA_LAYERS
         do p=1,PTRACE
          do ii=1,NPOOLS
            soilcpools(p,ii,n) = entcell%entcell%Tpool(p,ii,n)
          end do
         end do
        enddo
      endif

      end subroutine ent_get_exports_single


      subroutine ent_get_exports_array_1d( entcell,
     &     canopy_conductance,
     &     beta_soil_layers,
     &     shortwave_transmit,
     &     foliage_CO2,
     &     foliage_humidity,
     &     canopy_gpp,
     &     roughness_length,
     &     flux_CO2,
     &     albedo,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     vegetation_fractions,
     &     soilresp,
     &     soilcpools
     &     )
      type(entcelltype_public), dimension(:), intent(in) :: entcell
      real*8, dimension(:), optional, intent(out) ::
     &     canopy_conductance,
     &     shortwave_transmit,
     &     foliage_CO2,
     &     foliage_humidity,
     &     canopy_gpp,
     &     roughness_length,
     &     flux_CO2,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     soilresp
      real*8, dimension(:,:), optional, intent(out) ::
     &     beta_soil_layers,
     &     albedo,
     &     vegetation_fractions
      real*8, dimension(:,:,:,:), optional, intent(out) ::
     &     soilcpools 
      !----------
      integer n, i, ic, p, ii

      ic = size(entcell, 1)

      do i=1,ic
      if ( present(canopy_conductance) )
     &     canopy_conductance(i) = entcell(i)%entcell%gcanopy

      if ( present(shortwave_transmit) )
     &     shortwave_transmit(i) = entcell(i)%entcell%TRANS_SW

      if ( present(foliage_CO2) )
     &     foliage_CO2(i) = entcell(i)%entcell%Ci

      if ( present(foliage_humidity) )
     &     foliage_humidity(i) = entcell(i)%entcell%Qf

      if ( present(canopy_gpp) )
     &     canopy_gpp(i) = entcell(i)%entcell%GPP

      if ( present(roughness_length) )
     &     roughness_length(i) = entcell(i)%entcell%z0

      if ( present(flux_CO2) )
     &     flux_CO2(i) = entcell(i)%entcell%CO2flux

      if ( present(canopy_max_H2O) )
     &     canopy_max_H2O(i) = entcell(i)%entcell%LAI * .0001d0 !!! GISS setting

      if ( present(canopy_heat_capacity) ) then
        !aa=ala(1,i0,j0)
        !canopy_heat_capacity=(.010d0+.002d0*aa+.001d0*aa**2)*shw
        call stop_model("not implemmented yet",255)
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        ! compute it here ?
        call stop_model("not implemmented yet",255)
      endif

      if ( present(soilresp) )  !PK
     &     soilresp(i) = entcell(i)%entcell%Soil_resp

      if ( present(beta_soil_layers) ) then
        do n=1,N_DEPTH
          beta_soil_layers(n,i) = entcell(i)%entcell%betadl(n)
        enddo
      endif

      if ( present(albedo) ) then
        do n=1,N_BANDS
          albedo(n,i) = entcell(i)%entcell%albedo(n)
        enddo
      endif

      if ( present(vegetation_fractions) ) then
        do n=1,N_COVERTYPES
        ! extract those here ?
        !vegetation_fractions(:) = vdata(i,j,:)
          call stop_model("not implemmented yet",255)
        enddo
      endif

      if ( present(soilcpools) ) then  !PK
        do n=1,N_CASA_LAYERS
         do p=1,PTRACE
          do ii=1,NPOOLS
            soilcpools(p,ii,n,i) = 
     &      entcell(i)%entcell%Tpool(p,ii,n)
          end do
         end do
        enddo
      endif

      enddo

      end subroutine ent_get_exports_array_1d


      subroutine ent_get_exports_array_2d( entcell,
     &     canopy_conductance,
     &     beta_soil_layers,
     &     shortwave_transmit,
     &     foliage_CO2,
     &     foliage_humidity,
     &     canopy_gpp,
     &     roughness_length,
     &     flux_CO2,
     &     albedo,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     vegetation_fractions,
     &     soilresp,
     &     soilcpools
     &     )
      type(entcelltype_public), dimension(:,:), intent(in) :: entcell
      real*8, dimension(:,:), optional, intent(out) ::
     &     canopy_conductance,
     &     shortwave_transmit,
     &     foliage_CO2,
     &     foliage_humidity,
     &     canopy_gpp,
     &     roughness_length,
     &     flux_CO2,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     soilresp
      real*8, dimension(:,:,:), optional, intent(out) ::
     &     beta_soil_layers,
     &     albedo,
     &     vegetation_fractions
      real*8, dimension(:,:,:,:,:), optional, intent(out) ::
     &     soilcpools
      !----------
      integer n, i, j, ic, jc, p,ii

      ic = size(entcell, 1)
      jc = size(entcell, 2)

      do j=1,jc
      do i=1,ic
      if ( present(canopy_conductance) )
     &     canopy_conductance(i,j) = 
     &       entcell(i,j)%entcell%GCANOPY

      if( j > 32768 ) then
        print*,'Nancys compiler needs print stmt in ent_mod here'
      endif


      if ( present(shortwave_transmit) )
     &     shortwave_transmit(i,j) = 
     &     entcell(i,j)%entcell%TRANS_SW

      if ( present(foliage_CO2) )
     &     foliage_CO2(i,j) = entcell(i,j)%entcell%Ci

      if ( present(foliage_humidity) )
     &     foliage_humidity(i,j) = entcell(i,j)%entcell%Qf

      if ( present(canopy_gpp) )
     &     canopy_gpp(i,j) = entcell(i,j)%entcell%GPP

      if ( present(roughness_length) )
     &     roughness_length(i,j) = entcell(i,j)%entcell%z0

      if ( present(flux_CO2) )
     &     flux_CO2(i,j) = entcell(i,j)%entcell%CO2flux

      if ( present(canopy_max_H2O) )
     &     canopy_max_H2O(i,j) = 
     &     entcell(i,j)%entcell%LAI * .0001d0 !!! GISS setting

      if ( present(canopy_heat_capacity) ) then
        !aa=ala(1,i0,j0)
        !canopy_heat_capacity=(.010d0+.002d0*aa+.001d0*aa**2)*shw
        call stop_model("not implemmented yet",255)
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        ! compute it here ?
        !call stop_model("not implemmented yet",255)
        if ( associated(entcell(i,j)%entcell) ) then
          fraction_of_vegetated_soil(i,j) =
     &         entcell(i,j)%entcell%fv
        else
          fraction_of_vegetated_soil(i,j) = 0.d0
        endif
      endif

      if ( present(soilresp) )
     &     soilresp(i,j) = entcell(i,j)%entcell%Soil_resp

      if ( present(beta_soil_layers) ) then
        do n=1,N_DEPTH
          beta_soil_layers(n,i,j) = 
     &         entcell(i,j)%entcell%betadl(n)
        enddo
      endif

      if ( present(albedo) ) then
        do n=1,N_BANDS
          albedo(n,i,j) = entcell(i,j)%entcell%albedo(n)
        enddo
      endif

      if ( present(vegetation_fractions) ) then
        do n=1,N_COVERTYPES
        ! extract those here ?
        !vegetation_fractions(:) = vdata(i,j,:)
          call stop_model("not implemmented yet",255)
        enddo
      endif

      if ( present(soilcpools) ) then
        do n=1,N_CASA_LAYERS
         do p=1,PTRACE
          do ii=1,NPOOLS
            soilcpools(p,ii,n,i,j) = 
     &      entcell(i,j)%entcell%Tpool(p,ii,n)
          end do
         end do
        enddo
      endif

      enddo  !i
      enddo  !j

      end subroutine ent_get_exports_array_2d

!******************************************************************

      subroutine ent_cell_print_single(iu, entcell)
      use entcells, only : entcell_destruct
      integer, intent(in) :: iu
      type(entcelltype_public), intent(in) :: entcell

      if ( .not. associated(entcell%entcell) ) then
        print *, "ent_cell_print_single: Empty entcell"
        return
      endif

      call entcell_print( iu, entcell%entcell )

      end subroutine ent_cell_print_single


      subroutine ent_cell_print_array_1d(iu, entcell)
      integer, intent(in) :: iu
      type(entcelltype_public), intent(in) :: entcell(:)
      integer n, nc

      nc = size(entcell)

      do n=1,nc
        call ent_cell_print_single( iu, entcell(n) )
      enddo

      end subroutine ent_cell_print_array_1d


      subroutine ent_cell_print_array_2d(iu, entcell)
      integer, intent(in) :: iu
      type(entcelltype_public), intent(in) :: entcell(:,:)
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_cell_print_single( iu, entcell(i,j) )
        enddo
      enddo

      end subroutine ent_cell_print_array_2d


      end module ent_mod


