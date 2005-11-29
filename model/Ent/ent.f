      module ent
      
!@sum Contains main routines to perform single time step Ent model simulation
!@sum on a single grid cell, entcell, whose parameters are set previously by
!@sum ent_driver routines and the subroutine ent_model called by a 
!@sum main program.

!@auth N. Kiang

      !Ent MODULES TO USE
      use ent_const
      use ent_types

      implicit none
      private
      save

      public ent_integrate

      contains
      !*********************************************************************
      subroutine ent_integrate(dtsec, time, entcell)
      use phenology
      use disturbance
      use canopyrad
      use disturbance
      use reproduction
      use cohorts
      use patches
      use util

      implicit none
      real*8 :: dtsec  !dt in seconds
      type(timestruct) :: time !Time in year.fraction, Greenwich Mean Time
      type(entcelltype) :: entcell

      !dtsec = dtyr*YEARSEC(entdata%tt%year)
      !Convert to local time here?
      !---------

      !Flag when it's time to update vegetation structure
      !May be at set time intervals, or function of biomass accumulation, etc.

      write(*,*) 'Ecosystem dynamics for (long,lat)=(',
     & entcell%long,entcell%lat,'),time=',time

      if (STRUCT_FLAG(time,entcell)) then
        !* Update phenology and disturbance
        call phenology_update (dtsec,time, entcell)
        call fire_frequency (dtsec,time, entcell)
        call recalc_radpar (entcell)
      end if

      call calc_cell_disturbance_rates(dtsec,time,entcell)
      call ent_bgc(dtsec,time,entcell)

      if (STRUCT_FLAG(time,entcell)) then
        call reproduction_calc(dtsec, time, entcell)
        call reorganize_cohorts(entcell)
      end if

      if (STRUCT_FLAG(time,entcell)) then
        call reorganize_patches(entcell)
      end if

      end subroutine ent_integrate

      !*********************************************************************

      function STRUCT_FLAG(time, entcell) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple beginning-of-the-month flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: time
        type(entcelltype) :: entcell
        logical :: update_struct
        
        real*8 :: hourfrac
        
        hourfrac = time%hour + time%minute/60.0 + time%seconds/3600.0
        if ((time%day.eq.1).and.
     &       (hourfrac.eq.0.0)) then
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG

      !*********************************************************************
      subroutine ent_bgc(dtsec, time, entcell)
!@sum Calculate canopy conductance and fluxes of CO2 and N, update pools.

      use biophysics
      use growthallometry
      use phenology
      use soilbgc

      implicit none

      real*8 :: dtsec
      type(timestruct) :: time
      type(entcelltype) :: entcell
      type(patch),pointer :: pp

      pp = entcell%oldest
      do while (allocated(pp))
        call photosynth_cond(dtsec, pp)
        call uptake_N(dtsec, pp)
        call litter(dtsec, time, pp)
        call soil_bgc(dtsec, pp)
        pp = pp%younger
      end do

      end subroutine ent_bgc

      subroutine ent_bgc_GISShack(dtsec, time, entcell)
!@sum Calculate canopy conductance and fluxes of CO2 and N, update pools.
!@sum Temporary hack to replicate the grid-cell-level averaging of vegetation
!@sum properties done in the GISS GCM. 

      use biophysics
      use growthallometry
      use phenology
      use soilbgc
      use patches

      implicit none

      real*8 :: dtsec
      type(timestruct) :: time
      type(entcelltype) :: entcell
      !--------Local vars--------
      type(patch),pointer :: pp
      type(patch),pointer :: tempp

      real*8 :: sfv, salai, svh, snm, snf
      real*8 :: tfv

      call allocate(tempp)
      call init_patch(tempp,entcell,entcell%area)
      sfv = 0.0
      salai = 0.0
      svh = 0.0
      snm = 0.0
      snf = 0.0
      almass = 0.0
      tfv = 0.0

      !Average vegetation properties to grid cell level and put into a
      !dummy patch structure to pass to biophysics.
      pp = entcell%oldest
      do while (allocated(pp))
        sfv = pp%area/pp%cellptr%area  !vegfraction
        tfv = tfv + sfv                !Make sure vegfraction adds up to 1.
        salai = salai + sfv*pp%tallest%LAI !Weighted average by vegfraction
        svh = svh + sfv*pp%tallest%h
        snm = snm + sfv*pp%tallest%nm
        snf = snf + sfv*pfpar(pp%tallest%pft)%nf
        
        pp = pp%younger
      end do
      salai = salai/tfv                !Account for bare soil not covered.
      svh = svh/tfv
      snm = snm/tfv
      snf = snf/tfv

      call sum_roots_patches2cell(entcell)

      !* Put entcell grid-average values into a hack patch data structure.
      !* One cohort with grid-averaged values.
      !* Dummy 1.0 values are used for parameters not relevant to GISS ModelE
      call insert_cohort(tempp,entcell%oldest%tallest%pft,
     &     1.0, svh, 1.0, 1.0, 1.0, 1.0, salai, 1.0, entcell%froot,
     &     1.0, 1.0,1.0,1.0,1.0,1.0,1.0,
     &     
       
      call photosynth_cond(dtsec, tempp)

      end subroutine ent_bgc_GISShack

!*****************************************************************************
 
      end module ent
