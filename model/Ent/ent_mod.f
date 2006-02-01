      module ent_mod
!@sum this module contains the interface to Ent dynamic vegetation model
!@+   this is the only module that should be visible to GCM
 
      ! need to know about internal structure of Ent types
      use ent_types, only : entcelltype, patch, cohort
      implicit none

      private

      public entcelltype_public, ent_cell_pack, ent_cell_unpack
      public ent_get_results, ent_set_forcings
      public ent_run, ent_run_slow_pysics

      type entcelltype_public
        private
        type(entcelltype) :: entcell
      end type entcelltype_public

      interface ent_set_forcings
       module procedure ent_set_forcings_single
       module procedure ent_set_forcings_array
      end interface

      interface ent_get_results
       module procedure ent_get_results_single
       module procedure ent_get_results_array
      end interface

      interface ent_run
        module procedure ent_run_single
        module procedure ent_run_array
      end interface

      interface copy_vars
        copy_vars_single
        copy_vars_array
      end interface

      contains

      subroutine ent_run_single(entcell)
      use ent_driver, only : ent_model
      type(entcelltype_public), intent(inout) :: entcell
      !---

      call ent_model( entcell%entcell )

      end subroutine ent_run_single


      subroutine ent_run_array(entcell)
      use ent_driver, only : ent_model
      type(entcelltype_public), intent(inout) :: entcell(:)
      !---
      integer n, nc

      nc = size(entcell)
      do n=1,nc
        call ent_model( entcell(n)%entcell )
      enddo

      end subroutine ent_run_array



      subroutine ent_run_slow_pysics(entcell,
     &     jday
! insert any needed input parameters here
     &     )
!@sum this call updates variable that change on a long time scale.
!@+   Right now (before real dynamic vegetation is implemented)
!@+   it should perform prescribed seasonal update of vegatation
!@+   parameters (LAI, root fraction etc.)
!@+   I think extra input parameters needed here should be passed 
!@+   as formal parameters and not be packed into entcell structure.
!@+   Is it OK from ESMF point of view?
      use ent_driver, only : ent_update_veg_structure
      type(entcelltype_public), intent(in) :: entcell
      integer jday
      !---
      
      ent_update_veg_structure( entcell%entcell, jday )

      end subroutine ent_run_slow_pysics


      subroutine ent_cell_pack(ibuf, dbuf, entcell)
!@sum allocate two linear arrays ibuf, dbuf and pack contents of
!@+   entcell into them
      integer, pointer, intent(out) :: ibuf(:)
      real*8, pointer, intent(out) :: dbuf(:)
      type(entcelltype_public), intent(in) :: entcell ! pointer ?
      !---
      ! allocate tmp arrays of a fixed size MAX_NIBIF+1, MAX_NDBUF+1
! the MAX_* params may need to be relocated to better palce (ent_types?)
!@var MAX_PATCHES maximal number of patches per cell
!@var MAX_COHORTS maximal number of cohorts per patch
      integer, parameter :: MAX_PATCHES=32, MAX_COHORTS=64
      type(patch), pointer :: p  !@var p current patch
      type(cohort), pointer :: c !@var current cohort
      integer :: np              !@var np number of patches in the cell
      integer :: nc(MAX_PATCHES) !@var nc number of cohorts in the patch
      integer :: ic, dc


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
        call copy_patch_vars(dbuf, nn, p, 0); ndbuf = ndbuf + nn
        nc(np) = 0
        c => p%tallest
        do while ( associated(c) )
          nc(np) = nc(np) + 1
          if ( nc(np) > MAX_COHORTS )
     &         call stop_model("ent_cell_pack: too many cohorts",255)
          !save cohort
          !dbuf(dc) = c%_any_value_ ; dc = dc + 1
          call copy_cohort_vars(dbuf, nn, c, 0); ndbuf = ndbuf + nn
          c => c%shorter
        enddo
        p => p%younger
      enddo

      ! ibuf contains np and nc(np)
      allocate( ibuf(0:np) )
      ic = 0
      ibuf(ic) = np; ic = ic + 1
      ibuf(ic:ic+np-1) = nc(1:np); ic = ic + 1

      allocate( dbuf(0:ndbuf-1) )
      dc = 0

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


      subroutine ent_cell_unpack(ibuf, dbuf, entcell)
! this program is not finished yet: have to assign all the pointers
      integer, intent(in) :: ibuf(0:)
      real*8, intent(in) :: dbuf(0:)
      type(entcelltype_public), intent(out) :: entcell ! pointer ?
      !---
      type(patch), pointer :: p, pprev  !@var p current patch
      type(cohort), pointer :: c, cprev !@var current cohort
      integer :: np              !@var np number of patches in the cell
      integer, allocatable :: nc(:) !@var nc number of cohorts in the patch
      integer ic, dc
      integer i, j
      integer npdebug, ncdebug ! these are for debuging

      ic = 0; dc = 0

      ! doesn't seem that we need to restore anything for the cell
      np = ibuf(ic); ic = ic + 1
      allocate( nc(np) )
      nc(1:np) = ibuf(ic:ic+np-1); ic = ic + np

      if ( np <= 0 ) return  ! nothing to restore...

      nullify( pprev )
      do i=1,np
        allocate( p )
        call copy_patch_vars(dbuf(dc:), nn, p, 1); dc = dc + nn
        p%older => pprev
        nullify( cprev)
        do j=1,nc(i)
          allocate( c )
          call copy_cohort_vars(dbuf, nn, c, -1); dc = dc + nn
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
      real*8, intent(out) :: buf(:)
      integer, intent(out) :: n
      real*8, intent(in):: var
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
      real*8, intent(out) :: buf(:)
      integer, intent(out) :: n
      real*8, intent(in):: var(:)
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
      real*8, intent(out) :: buf(0:)
      integer, intent(out) :: n
      type(patch), intent(in):: p
!@var flag defines the actual action:
!@+     -1 copy from patch to buffer
!@+      1 copy from buffer to patch
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      integer dc

      dc = 0

      ! include all patch variables that need i/o
      call copy_vars( buf(dc:), nn,  p%age,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  p%area, flag ); dc = dc + nn

      n = dc

      end subroutine copy_patch_vars


      subroutine copy_cohort_vars(buf, n, c, flag)
      real*8, intent(out) :: buf(0:)
      integer, intent(out) :: n
      type(cohort), intent(in):: c
!@var flag defines the actual action:
!@+     -1 copy from patch to buffer
!@+      1 copy from buffer to patch
!@+      0 do nothing - just return the number of fields
      integer, intent(in) :: flag
      !---
      integer dc

      dc = 0

      ! include all cohort variables that need i/o
      call copy_vars( buf(dc:), nn,  c%_any_var_,  flag ); dc = dc + nn
      call copy_vars( buf(dc:), nn,  c%_any_var2_, flag ); dc = dc + nn

      n = dc

      end subroutine copy_cohort_vars


!******************************************************************


      subroutine ent_set_forcings_single( entcell,
     &     canopy_temperature,
     &     direct_visible_rad,
     &     total_visible_rad,
     &     )
      type(entcelltype_public), intent(in) :: entcell
      ! forcings probably should not be optional ...
      real*8 :: canopy_temperature
      real*8 :: direct_short_wave
      !----------
      real*8 alai

      entcell%entcell%TcanopyC = canopy_temperature
      entcell%entcell%Idir     = total_visible_rad
      entcell%entcell%Ivis     = direct_visible_rad

      end subroutine ent_set_forcings_single

      subroutine ent_set_forcings_array( entcell,
     &     canopy_temperature,
     &     direct_visible_rad,
     &     total_visible_rad,
     &     )
      type(entcelltype_public), intent(in) :: entcell(:)
      ! forcings probably should not be optional ...
      real*8 :: canopy_temperature(:)
      real*8 :: direct_visible_rad(:)
      real*8 :: total_visible_rad(:)
     !----------
      real*8 alai

      entcell(:)%entcell%TcanopyC = canopy_temperature(:)
      entcell(:)%entcell%Idir     = total_visible_rad(:)
      entcell(:)%entcell%Ivis     = direct_visible_rad(:)

      end subroutine ent_set_forcings_array



      subroutine ent_get_results_single( entcell,
     &     canopy_conductance,
     &     canopy_gpp,
     &     canopy_max_H2O,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil,
     &     vegetation_fractions
     &     )
      type(entcelltype_public), intent(in) :: entcell
      real*8, optional :: canopy_conductance
      real*8, optional :: canopy_gpp
      real*8, optional :: canopy_holding_capacity
      real*8, optional :: canopy_heat_capacity
      real*8, optional :: fraction_of_vegetated_soil
      real*8, optional, dimension(:) :: vegetation_fractions
      !----------
      real*8 alai

      if ( present(canopy_conductance) ) then
        canopy_conductance = entcell%entcell%GCANOPY
      endif

      if ( present(canopy_gpp) ) then
        canopy_gpp = entcell%entcell%GPP
      endif

      ! the values below are needed by GISS GCM
      ! haven't decided yet what is the best way to extract them

      if ( present(canopy_holding_capacity) ) then
        alai = entcell%entcell%alai
        canopy_max_H2O = .0001d0 * alai
      endif

      if ( present(canopy_heat_capacity) ) then
        !aa=ala(1,i0,j0)
        !canopy_heat_capacity=(.010d0+.002d0*aa+.001d0*aa**2)*shw
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        !fraction_of_vegetated_soil = sum( vdata(i,j,2:9 )
      endif

      if ( present(vegetation_fractions) ) then
        !vegetation_fractions(:) = vdata(i,j,:)
      endif


      end subroutine ent_get_results_single


      subroutine ent_get_results_array( entcell,
     &     canopy_conductance,
     &     canopy_gpp
     &     )
      type(entcelltype_public), intent(in) :: entcell(:)
      real*8, optional :: canopy_conductance(:)
      real*8, optional :: canopy_gpp(:)
      !----------
      real*8 alai

!!!! Can I rely on correctly passed dimensions here, i.e.
!!!  is (:) always the correct extent ?

      if ( present(canopy_conductance) ) then
        canopy_conductance(:) = entcell(:)%entcell%GCANOPY
      endif

      if ( present(canopy_gpp) ) then
        canopy_gpp(:) = entcell(:)%entcell%GPP
      endif

      end subroutine ent_get_results_array


      end module ent_mod



#ifdef DEBUG_PACK_UNPACK
      program foo
      use ent_mod
      type ent_storage
        real*8, pointer :: dbuf
        iteger, pointer :: ibuf
      end type
      type(ent_storage) :: buffer(NUMCELLS)
      type(entcelltype_public) entcells(NUMCELLS)
      integer i
      real*8 Tcan(NUMCELLS),Idir(NUMCELLS),Ivis(NUMCELLS) ! forcings
      real*8 cond(NUMCELLS)                               ! results

      ! to pack ent cells entcells(1:NUMCELLS) into arrays dbuf, ibuf :
      do i=1,NUMCELLS
        ! this call also allocates ibuf, dbuf
        ! entcells is left untouched
        call ent_cell_pack( buffer(i)%ibif, buffer(i)%dbuf, entcells(i) )
      enddo
      ! ic and dc contain total lengths of arrays used

      !to unpack:
      do i=1,NUMCELLS
        ! this call also allocates structures inside entcells
        ! buffer(i)%ibif, buffer(i)%dbuf are left untouched
        call ent_cell_unpack( buffer(i)%ibif, buffer(i)%dbuf, entcells(i) )
      enddo

!!! One should be carefull not unpacking into allocated cells,
!!! since the old memory will not be de-allocated and this will
!!! lead to a memor leak

! and here is an example of how one performs one time step

      call ent_set_forcings( entcells,
     &     canopy_temperature=Tcan,
     &     direct_visible_rad=Idir,
     &     total_visible_rad=Ivis,
     &     )

      call ent_run( entcells )

      call ent_get_results( entcells,
     &     canopy_conductance=cond
     &     )

      end

#endif


#ifdef COMMENTED_AREA
      Questions:

 1    How are we going to deal with bare soil? Will it be an empty
      cell?, i.e. np = 0
      Or should it be treated as a type of vegetation?

 2    Who allocates entcell_public (i.e. the the cell itself, not
      its contents) ? (in above example it is allocated by GCM).


#end
