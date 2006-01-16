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
      subroutine ent_ecosystem_dynamics(dtsec,time,ecp)
!@sum Ent ecosystem dynamics.
      use phenology
      use disturbance
      use canopyrad
      use disturbance
      
!      write(*,*) 'Ecosystem dynamics for (long,lat)=(',
!     & ecp%long,ecp%lat,'),time=',time

      call ent_integrate(dtsec,time,ecp) !Biophysics, growth/allom, reproduction
      
      !Flag when it's time to update disturbance.
      !May be at set time intervals, or function of biomass accumulation, etc.
      !For now, monthly update as a place holder.
      if (STRUCT_FLAG_MONTH(time,ecp)) then
        !* Update phenology and disturbance
        !call phenology_update (dtsec,time, pp) !UPDATE LAI - put in ent_integrate
        call fire_frequency_cell (dtsec,time, ecp) !DUMMY
        call recalc_radpar_cell (ecp) !
        call reorganize_patches(ecp)
        call calc_cell_disturbance_rates(dtsec,time,ecp)
      else
        call calc_cell_disturbance_rates(dtsec,time,ecp)
      end if

      end subroutine ent_ecosystem_dynamics

      !*********************************************************************
      subroutine ent_integrate(dtsec, time, ecp)
!@sum Ent biophysics/biogeochemistry
      use reproduction
      use cohorts
      use patches
      use util

      implicit none
      real*8 :: dtsec  !dt in seconds
      type(timestruct) :: time !Time in year.fraction, Greenwich Mean Time
      type(entcell),pointer :: ecp
      !-----local--------
      type(patch),pointer :: pp

      pp = ecp%youngest
      do while (ASSOCIATED(pp)) 
        call photosynth_cond(dtsec, pp)
        call uptake_N(dtsec, pp)
        call litter(dtsec, time, pp)
        call soil_bgc(dtsec, pp)
        if (STRUCT_FLAG(time,ecp)) then
          call reproduction_calc(dtsec, time, pp)
          call reorganize_cohorts(pp)
        end if
      end do

      if (STRUCT_FLAG(time,pp%cellptr)) then
        call phenology_update (dtsec,time, pp) !UPDATE LAI
        call recalc_radpar (pp) !UPDATE canopy radiative transfer
      end if

      call summarize_patches(time,ecp)
      end subroutine ent_integrate

      !*********************************************************************

      function STRUCT_FLAG(time, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
        type(timestruct) :: time
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !------local------
        update_struct = STRUCT_FLAG_DAY(time,ecp)

      end function STRUCT_FLAG
      !*********************************************************************

      function STRUCT_FLAG_DAY(time, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple end-of-day flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: time
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !-----local----------
        real*8 :: hourfrac
        
        hourfrac = time%hour + time%minute/60.0 + time%seconds/3600.0
!        if (hourfrac.le.dtsec) then !Midnight
        if (hourfrac.eq.0.0)) then  !Midnight
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_DAY

      !*********************************************************************

      function STRUCT_FLAG_MONTH(time, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple beginning-of-the-month flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: time
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        
        real*8 :: hourfrac
        
        hourfrac = time%hour + time%minute/60.0 + time%seconds/3600.0
        if ((time%day.eq.1).and.
     &       (hourfrac.eq.0.0)) then
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_MONTH

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
