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

      public ent_ecosystem_dynamics

      contains
      !*********************************************************************
      subroutine ent_ecosystem_dynamics(dtsec,tt,ecp)
!@sum Ent ecosystem dynamics.
      use phenology
      use disturbance
      use canopyrad
      use disturbance

      real*8,intent(in) :: dtsec
      type(timestruct),pointer :: tt
      type(entcelltype) :: ecp
!      write(*,*) 'Ecosystem dynamics for (long,lat)=(',
!     & ecp%long,ecp%lat,'),tt=',tt

      call ent_integrate(dtsec,tt,ecp) !Biophysics, growth/allom, reproduction
      
      !Flag when it's time to update disturbance.
      !May be at set time intervals, or function of biomass accumulation, etc.
      !For now, monthly update as a place holder.
      if (STRUCT_FLAG_MONTH(tt,ecp)) then
        !* Update phenology and disturbance
        !call phenology_update (dtsec,tt, pp) !UPDATE LAI - put in ent_integrate
        call fire_frequency_cell (dtsec,tt, ecp) !DUMMY
        call recalc_radpar_cell (ecp) !
        call reorganize_patches(ecp)
        call calc_cell_disturbance_rates(dtsec,tt,ecp)
      else
        call calc_cell_disturbance_rates(dtsec,tt,ecp)
      end if

      call summarize_entcell(ecp)

      end subroutine ent_ecosystem_dynamics

      !*********************************************************************
      subroutine ent_integrate(dtsec, tt, ecp)
!@sum Ent biophysics/biogeochemistry
      use reproduction
      use cohorts
      use patches
      use util

      implicit none
      real*8 :: dtsec  !dt in seconds
      type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      type(entcelltype) :: ecp
      !-----local--------
      type(patch),pointer :: pp

      !****Ent final version:  loop through patches
      !pp = ecp%youngest
      !do while (ASSOCIATED(pp)) 
      !***GISS version:  pass in ecp%sumpatch
      pp = ecp%sumpatch
        call photosynth_cond(dtsec, pp)
        call uptake_N(dtsec, pp)
        call litter(dtsec, tt, pp)
        call soil_bgc(dtsec, pp)
        if (STRUCT_FLAG(tt,ecp)) then
          call reproduction_calc(dtsec, tt, pp)
          call reorganize_cohorts(pp)
        end if
      !end do !****Ent final version

      if (STRUCT_FLAG(tt,pp%cellptr)) then
        call phenology_update (dtsec,tt, pp) !UPDATE LAI
        call recalc_radpar (tt,pp) !UPDATE canopy radiative transfer
      end if

      call summarize_patch(tt,pp)

      end subroutine ent_integrate

      !*********************************************************************

      function STRUCT_FLAG(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !------local------
        update_struct = STRUCT_FLAG_DAY(tt,ecp)

      end function STRUCT_FLAG
      !*********************************************************************

      function STRUCT_FLAG_DAY(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple end-of-day flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !-----local----------
        real*8 :: hourfrac
        
        hourfrac = tt%hour + tt%minute/60.0 + tt%seconds/3600.0
!        if (hourfrac.le.dtsec) then !Midnight
        if (hourfrac.eq.0.0) then  !Midnight
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_DAY

      !*********************************************************************

      function STRUCT_FLAG_MONTH(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple beginning-of-the-month flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        
        real*8 :: hourfrac
        
        hourfrac = tt%hour + tt%minute/60.0 + tt%seconds/3600.0
        if ((tt%day.eq.1).and.
     &       (hourfrac.eq.0.0)) then
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_MONTH


!*****************************************************************************
 
      end module ent
