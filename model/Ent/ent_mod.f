      module ent_mod
!@sum this module contains the interface to Ent dynamic vegetation model
!@+   this is the only module that should be visible to GCM
 
      ! need to know about internal structure of Ent types
      use ent_types, only : entcelltype, patch, cohort, timestruct,
     &     MAX_PATCHES, MAX_COHORTS
      use ent_const, only : N_BANDS, N_COVERTYPES, N_DEPTH
      !use ent_GISSveg
      use entcells
      implicit none

      private

      !--- public constants ---
      public N_BANDS, N_COVERTYPES, N_DEPTH

      public entcelltype_public, ent_cell_pack, ent_cell_unpack
      public ent_get_exports, ent_set_forcings
      public ent_cell_construct, ent_cell_destruct
      public ent_fast_processes,ent_seasonal_update,ent_vegcover_update
      public ent_cell_set
      public ent_prescribe_vegupdate

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

      !--- passing initial data to ent cells ---
      interface ent_cell_set
        module procedure ent_cell_set_single
        module procedure ent_cell_set_array_1d
        module procedure ent_cell_set_array_2d
      end interface

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


!!! do we need 1d and 2d array interfaces for pack/unpack ?

      !---- private interfaces ----
      interface copy_vars
        module procedure copy_vars_single
        module procedure copy_vars_array
      end interface

      contains

!---- interfaces to run the model one time step ----
      subroutine ent_prescribe_vegupdate(entcell,hemi,jday,year,
     &     update_crops)
      use ent_GISSveg, only:  ent_GISS_vegupdate
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer, intent(in) :: hemi(:,:)
      integer,intent(in) :: jday,year
      logical, intent(in) :: update_crops
      !---
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic      
          call ent_GISS_vegupdate(entcell(i,j)%entcell, hemi(i,j),
     &         jday, year, update_crops)
        enddo
      enddo

      end subroutine ent_prescribe_vegupdate
      


      subroutine ent_fast_processes_single(entcell, dt)
      use ent, only : ent_biophysics
      type(entcelltype_public), intent(inout) :: entcell
      real*8, intent(in) :: dt
      !---
      
      call ent_biophysics(dt, entcell%entcell)

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


      subroutine ent_fast_processes_array_2d(entcell, dt)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      real*8, intent(in) :: dt
      !---
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_fast_processes_single( entcell(i,j), dt )
        enddo
      enddo

      end subroutine ent_fast_processes_array_2d


      subroutine ent_seasonal_update_single(entcell,
     &     jday
! insert any needed input parameters here
     &     )
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
      type(entcelltype_public), intent(in) :: entcell
      integer, intent(in) :: jday !@var jday Julian day of the year
      !---
      
      !call ent_update_veg_structure( entcell%entcell, jday )

      end subroutine ent_seasonal_update_single


      subroutine ent_seasonal_update_array_1d(entcell, jday)
      type(entcelltype_public), intent(inout) :: entcell(:)
      integer, intent(in) :: jday
      !---
      integer n, nc

      nc = size(entcell)
      do n=1,nc
        call ent_seasonal_update_single( entcell(n), jday )
      enddo

      end subroutine ent_seasonal_update_array_1d


      subroutine ent_seasonal_update_array_2d(entcell, jday)
      type(entcelltype_public), intent(inout) :: entcell(:,:)
      integer, intent(in) :: jday
      !---
      integer i, ic, j, jc

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_seasonal_update_single( entcell(i,j), jday )
        enddo
      enddo

      end subroutine ent_seasonal_update_array_2d


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


!---- Constructor / Destructor -----

      subroutine ent_cell_construct_single(entcell)
      use entcells, only : zero_entcell
      type(entcelltype_public), intent(inout) :: entcell

      allocate( entcell%entcell )

      ! don't allocate patches, just set all pointers to NULL
      nullify( entcell%entcell%youngest )
      nullify( entcell%entcell%oldest   )
      nullify( entcell%entcell%sumpatch )
      ! for now set all values o zero or defaults
      call zero_entcell(entcell%entcell)

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

      ic = size(entcell,1)
      jc = size(entcell,2)

      do j=1,jc
        do i=1,ic
          call ent_cell_construct_single( entcell(i,j) )
        enddo
      enddo

      end subroutine ent_cell_construct_array_2d


      subroutine ent_cell_destruct_single(entcell)
      use entcells, only : zero_entcell
      type(entcelltype_public), intent(inout) :: entcell
      type(patch), pointer :: p  !@var p current patch
      type(cohort), pointer :: c !@var current cohort
      integer np, nc


      np = 0
      p => entcell%entcell%oldest      
      do while ( associated(p) )
        np = np + 1
        if ( np > MAX_PATCHES )
     &       call stop_model("ent_cell_destruct: too many patches",255)
        nc = 0
        c => p%tallest
        do while ( associated(c) )
          nc = nc + 1
          if ( nc > MAX_COHORTS )
     &         call stop_model("ent_cell_destruct: too many chrts",255)
          ! destroy cohort here
          deallocate( c )
          c => c%shorter
        enddo
        ! destroy patch here
        if (associated(p%sumcohort) ) deallocate( p%sumcohort )
        deallocate( p )
        p => p%younger
      enddo

      ! destroy cell here
      if ( associated(entcell%entcell%sumpatch) )
     &     deallocate( entcell%entcell%sumpatch )

      deallocate( entcell%entcell )
      nullify( entcell%entcell )

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


!---- END of  Constructor / Destructor -----

      
      subroutine ent_cell_set_single(entcell,
     &     veg_fraction,
     &     leaf_area_index,
     &     pft_hights,
     &     pft_nmdata,
     &     pft_froots,
     &     pft_population_density
     &     )
      type(entcelltype_public), intent(out) :: entcell
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     veg_fraction,
     &     leaf_area_index
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     pft_hights,
     &     pft_nmdata,
     &     pft_population_density
      real*8, dimension(:,:)  :: pft_froots

      call init_simple_entcell( entcell%entcell,
     &     veg_fraction,  leaf_area_index,
     &     pft_hights, pft_nmdata, pft_froots, pft_population_density)
      
      end subroutine ent_cell_set_single


      subroutine ent_cell_set_array_1d(entcell,
     &     veg_fraction,
     &     leaf_area_index,
     &     pft_hights,
     &     pft_nmdata,
     &     pft_froots,
     &     pft_population_density
     &     )
      type(entcelltype_public), intent(out) :: entcell(:)
      real*8, dimension(:,:)  ::   ! dim=N_COVERTYPES, n
     &     veg_fraction,
     &     leaf_area_index
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     pft_hights,
     &     pft_nmdata,
     &     pft_population_density
      real*8, dimension(:,:)  :: pft_froots
      !---
      integer n, nc

      nc = size(entcell, 1)

      do n=1,nc
        call init_simple_entcell( entcell(n)%entcell,
     &       veg_fraction(:,n),
     &       leaf_area_index(:,n),
     &       pft_hights, pft_nmdata, pft_froots, pft_population_density)
      enddo
      
      end subroutine ent_cell_set_array_1d


      subroutine ent_cell_set_array_2d(entcell,
     &     veg_fraction,
     &     leaf_area_index,
     &     pft_hights,
     &     pft_nmdata,
     &     pft_froots,
     &     pft_population_density
     &     )
      type(entcelltype_public), intent(out) :: entcell(:,:)
      real*8, dimension(:,:,:)  ::   ! dim=N_COVERTYPES, n
     &     veg_fraction,
     &     leaf_area_index
      real*8, dimension(:)  ::   ! dim=N_COVERTYPES
     &     pft_hights,
     &     pft_nmdata,
     &     pft_population_density
      real*8, dimension(:,:)  :: pft_froots
      !---
      integer i, j, ic, jc

      ic = size(entcell, 1)
      jc = size(entcell, 2)

      do j=1,jc
        do i=1,ic
          call init_simple_entcell( entcell(i,j)%entcell,
     &         veg_fraction(:,i,j),
     &         leaf_area_index(:,i,j),
     &         pft_hights,pft_nmdata,pft_froots,pft_population_density)
        enddo
      enddo
      
      end subroutine ent_cell_set_array_2d



      subroutine ent_cell_pack(dbuf, entcell)
!@sum allocate two linear arrays ibuf, dbuf and pack contents of
!@+   entcell into them
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

      allocate( dbuf(0:ndbuf-1+1+np) ) !i.e. num reals + num int's
      dc = 0
      dbuf(dc) = real( np, kind(0d0) );               dc = dc + 1
      dbuf(dc:dc+np-1) = real( nc(1:np), kind(0d0) ); dc = dc + 1

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
        call copy_patch_vars(dbuf(dc:), nn, p, -1); dc = dc + nn
        nc(np) = 0
        c => p%tallest
        do while ( associated(c) )
          nc(np) = nc(np) + 1
          if ( nc(np) > MAX_COHORTS )
     &         call stop_model("ent_cell_pack: too many cohorts",255)
          !save cohort
          call copy_cohort_vars(dbuf, nn, c, -1); dc = dc + nn
         c => c%shorter
        enddo
        p => p%younger
      enddo

      end subroutine ent_cell_pack


      subroutine ent_cell_unpack(dbuf, entcell)
! this program is not finished yet: have to assign all the pointers
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
      allocate( nc(np) )
      nc(1:np) = nint( dbuf(dc:dc+np-1) ); dc = dc + np

      if ( np <= 0 ) return  ! nothing to restore...

      nullify( pprev )
      do i=1,np
        allocate( p )
        call copy_patch_vars(dbuf(dc:), nn, p, 1); dc = dc + nn
        p%older => pprev
        nullify( cprev)
        do j=1,nc(i)
          allocate( c )
          call copy_cohort_vars(dbuf, nn, c, 1); dc = dc + nn
          c%taller => cprev
          cprev => c
        enddo
        p%shortest => c
        pprev => p
      enddo
      entcell%entcell%youngest => p

      ! now restore pointer lists in opposite direction
      npdebug = 0
      nullify( pprev )
      p => entcell%entcell%youngest
      do while ( associated(p) )
        p%younger => pprev
        npdebug = npdebug + 1
        if ( npdebug > np )
     &       call stop_model("ent_cell_unpack: broken struct: np",255)
        ncdebug = 0
        nullify( cprev)
        c => p%shortest
        do while ( associated(c) )
          c%shorter => cprev
          ncdebug = ncdebug + 1
          if ( ncdebug > nc(npdebug) )
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

      end subroutine ent_cell_unpack


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
      integer dc, nn

      dc = 0

      ! include all patch variables that need i/o
      call copy_vars( buf(dc:), nn,  p%age,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  p%area, flag ); dc = dc + nn

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
      call copy_vars( buf(dc:), nn,  c%lai,  flag ); dc = dc + nn
!      call copy_vars( buf(dc:), nn,  c%_any_var2_, flag ); dc = dc + nn

      n = dc

      end subroutine copy_cohort_vars


!******************************************************************


      subroutine ent_set_forcings_single( entcell,
     &     canopy_temperature,
     &     canopy_air_humidity,      
     &     surf_pressure,            
     &     surf_CO2,                 
     &     precip,                   
     &     heat_transfer_coef,       
     &     wind_speed,               
     &     total_visible_rad,
     &     direct_visible_rad,
     &     solar_zenith_angle,
     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
     &     ) ! need to pass Ci, Qf ??
      type(entcelltype_public), intent(out) :: entcell
      ! forcings probably should not be optional ...
      real*8, intent(in) ::
     &     canopy_temperature,
     &     canopy_air_humidity,
     &     surf_pressure,
     &     surf_CO2,
     &     precip,
     &     heat_transfer_coef,
     &     wind_speed,
     &     total_visible_rad,
     &     direct_visible_rad,
     &     solar_zenith_angle
      real*8, dimension(:), intent(in) ::
     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
      !----------
      integer n

      entcell%entcell%TcanopyC = canopy_temperature
      entcell%entcell%Qv = canopy_air_humidity
      entcell%entcell%P_mbar = surf_pressure
      entcell%entcell%Ca = surf_CO2
      entcell%entcell%Precip = precip
      entcell%entcell%Ch = heat_transfer_coef
      entcell%entcell%U = wind_speed
      entcell%entcell%Idir = total_visible_rad
      entcell%entcell%Ivis = direct_visible_rad
      entcell%entcell%Solarzen = solar_zenith_angle
      do n=1,N_DEPTH
        entcell%entcell%Soilmoist(n) = soil_water(n)
        entcell%entcell%Soilmp(n) = soil_matric_pot(n)
        entcell%entcell%fice(n) = soil_ice_fraction(n)
      enddo

      end subroutine ent_set_forcings_single


      subroutine ent_set_forcings_array_1d( entcell,
     &     canopy_temperature,
     &     canopy_air_humidity,      
     &     surf_pressure,            
     &     surf_CO2,                 
     &     precip,                   
     &     heat_transfer_coef,       
     &     wind_speed,               
     &     total_visible_rad,
     &     direct_visible_rad,
     &     solar_zenith_angle,
     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
     &     ) ! need to pass Ci, Qf ??
      type(entcelltype_public), dimension(:), intent(out) :: entcell
      ! forcings probably should not be optional ...
      real*8, dimension(:)  ::
     &     canopy_temperature,
     &     canopy_air_humidity,
     &     surf_pressure,
     &     surf_CO2,
     &     precip,
     &     heat_transfer_coef,
     &     wind_speed,
     &     total_visible_rad,
     &     direct_visible_rad,
     &     solar_zenith_angle
      real*8, dimension(:,:), intent(in) ::
     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
      !----------
      integer n
      integer i, ic

      ic = size(entcell, 1)

      do i=1,ic
        entcell(i)%entcell%TcanopyC = canopy_temperature(i)
        entcell(i)%entcell%Qv = canopy_air_humidity(i)
        entcell(i)%entcell%P_mbar = surf_pressure(i)
        entcell(i)%entcell%Ca = surf_CO2(i)
        entcell(i)%entcell%Precip = precip(i)
        entcell(i)%entcell%Ch = heat_transfer_coef(i)
        entcell(i)%entcell%U = wind_speed(i)
        entcell(i)%entcell%Idir = total_visible_rad(i)
        entcell(i)%entcell%Ivis = direct_visible_rad(i)
        entcell(i)%entcell%Solarzen = solar_zenith_angle(i)
        do n=1,N_DEPTH
          entcell(i)%entcell%Soilmoist(n) = soil_water(n,i)
          entcell(i)%entcell%Soilmp(n) = soil_matric_pot(n,i)
          entcell(i)%entcell%fice(n) = soil_ice_fraction(n,i)
        enddo
      enddo

      end subroutine ent_set_forcings_array_1d


      subroutine ent_set_forcings_array_2d( entcell,
     &     canopy_temperature,
     &     canopy_air_humidity,      
     &     surf_pressure,            
     &     surf_CO2,                 
     &     precip,                   
     &     heat_transfer_coef,       
     &     wind_speed,               
     &     total_visible_rad,
     &     direct_visible_rad,
     &     solar_zenith_angle,
     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
     &     ) ! need to pass Ci, Qf ??
      type(entcelltype_public), dimension(:,:), intent(out) :: entcell
      ! forcings probably should not be optional ...
      real*8, dimension(:,:)  ::
     &     canopy_temperature,
     &     canopy_air_humidity,
     &     surf_pressure,
     &     surf_CO2,
     &     precip,
     &     heat_transfer_coef,
     &     wind_speed,
     &     total_visible_rad,
     &     direct_visible_rad,
     &     solar_zenith_angle
      real*8, dimension(:,:,:), intent(in) ::
     &     soil_water,
     &     soil_matric_pot,
     &     soil_ice_fraction
      !----------
      integer n
      integer i, j, ic, jc

      ic = size(entcell, 1)
      jc = size(entcell, 2)
      
      do j=1,jc
        do i=1,ic
          entcell(i,j)%entcell%TcanopyC = canopy_temperature(i,j)
          entcell(i,j)%entcell%Qv = canopy_air_humidity(i,j)
          entcell(i,j)%entcell%P_mbar = surf_pressure(i,j)
          entcell(i,j)%entcell%Ca = surf_CO2(i,j)
          entcell(i,j)%entcell%Precip = precip(i,j)
          entcell(i,j)%entcell%Ch = heat_transfer_coef(i,j)
          entcell(i,j)%entcell%U = wind_speed(i,j)
          entcell(i,j)%entcell%Idir = total_visible_rad(i,j)
          entcell(i,j)%entcell%Ivis = direct_visible_rad(i,j)
          entcell(i,j)%entcell%Solarzen = solar_zenith_angle(i,j)
          do n=1,N_DEPTH
            entcell(i,j)%entcell%Soilmoist(n) = soil_water(n,i,j)
            entcell(i,j)%entcell%Soilmp(n) = soil_matric_pot(n,i,j)
            entcell(i,j)%entcell%fice(n) = soil_ice_fraction(n,i,j)
          enddo
        enddo
      enddo

      end subroutine ent_set_forcings_array_2d




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
     &     vegetation_fractions
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
     &     fraction_of_vegetated_soil
      real*8, dimension(:), optional, intent(out) ::
     &     beta_soil_layers,
     &     albedo,
     &     vegetation_fractions
      !----------
      integer n

      if ( present(canopy_conductance) )
     &     canopy_conductance = entcell%entcell%gcanopy

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
        !canopy_heat_capacity=(.010d0+.002d0*aa+.001d0*aa**2)*shw
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        ! compute it here ?
      endif

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
        do n=1,N_COVERTYPES
        ! extract those here ?
        !vegetation_fractions(:) = vdata(i,j,:)
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
     &     vegetation_fractions
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
     &     fraction_of_vegetated_soil
      real*8, dimension(:,:), optional, intent(out) ::
     &     beta_soil_layers,
     &     albedo,
     &     vegetation_fractions
      !----------
      integer n, i, ic

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
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        ! compute it here ?
      endif

      if ( present(beta_soil_layers) ) then
        do n=1,N_DEPTH
          beta_soil_layers(n,i) = entcell(i)%entcell%betadl(n)
        enddo
      endif

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
     &     vegetation_fractions
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
     &     fraction_of_vegetated_soil
      real*8, dimension(:,:,:), optional, intent(out) ::
     &     beta_soil_layers,
     &     albedo,
     &     vegetation_fractions
      !----------
      integer n, i, j, ic, jc

      ic = size(entcell, 1)
      jc = size(entcell, 2)

      do j=1,jc
      do i=1,ic
      if ( present(canopy_conductance) )
     &     canopy_conductance(i,j) = entcell(i,j)%entcell%gcanopy

      if ( present(shortwave_transmit) )
     &     shortwave_transmit(i,j) = entcell(i,j)%entcell%TRANS_SW

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
     &     canopy_max_H2O(i,j) = entcell(i,j)%entcell%LAI * .0001d0 !!! GISS setting

      if ( present(canopy_heat_capacity) ) then
        !aa=ala(1,i0,j0)
        !canopy_heat_capacity=(.010d0+.002d0*aa+.001d0*aa**2)*shw
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        ! compute it here ?
      endif

      if ( present(beta_soil_layers) ) then
        do n=1,N_DEPTH
          beta_soil_layers(n,i,j) = entcell(i,j)%entcell%betadl(n)
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
        enddo
      endif
      enddo
      enddo

      end subroutine ent_get_exports_array_2d





      end module ent_mod


#ifdef COMMENTED_AREA
      Questions:

 1    How are we going to deal with bare soil? Will it be an empty
      cell?, i.e. np = 0
      Or should it be treated as a type of vegetation?


#endif
