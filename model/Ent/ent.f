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

      public ent_ecosystem_dynamics, ent_integrate_GISS, ent_biophysics
      public ent_integrate !Added by KIM

      contains
      !*********************************************************************
      subroutine ent_ecosystem_dynamics(dtsec,tt,ecp,ALBEDO_FLAG)
!@sum Ent ecosystem dynamics.
      use phenology
      use disturbance
      use canopyrad
      use disturbance
      use reproduction
      use cohorts, only : reorganize_cohorts
      use patches, only : reorganize_patches, summarize_patch
      use entcells, only : summarize_entcell

      real*8,intent(in) :: dtsec
      type(timestruct),pointer :: tt
      type(entcelltype) :: ecp
      logical, optional, intent(in) :: ALBEDO_FLAG
      !---
      type(patch), pointer :: pp
!      write(*,*) 'Ecosystem dynamics for (long,lat)=(',

      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
      !#### THIS LOOP: NEED TO REPLACE ALL CALLS WITH ecp TO CALLS WITH pp ##

        if (ALBEDO_FLAG) then
          call get_patchalbedo(pp)
        end if

!        call ent_integrate(dtsec,ecp,0.0) !Biophysics, respiration

        if (STRUCT_FLAG(tt,ecp)) then
          call reproduction_calc(dtsec, tt, pp)
          call reorganize_cohorts(pp)
!          call phenology_update (dtsec,tt, pp) !UPDATE LAI
          call recalc_radpar (pp) !UPDATE canopy radiative transfer
        end if
        call summarize_patch(pp)
      
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
        pp => pp%younger
      end do

      call summarize_entcell(ecp)

        end subroutine ent_ecosystem_dynamics

      !*********************************************************************
      subroutine ent_integrate_GISS(ecp, dtsec)
!@sum Ent biophysics/biogeochemistry - THIS SUBROUTINE WILL NOT BE NEEDED.
      use reproduction
      use cohorts
      use patches
      use biophysics, only : photosynth_cond
      use growthallometry, only : uptake_N
      use phenology, only : litter
      use soilbgc, only : soil_bgc
      use phenology, only : phenology_update
      use canopyrad, only : recalc_radpar
      use entcells, only : summarize_entcell, entcell_print

      implicit none
      type(entcelltype) :: ecp
      real*8 :: dtsec  !dt in seconds
!      type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      !-----local--------
      type(patch),pointer :: pp

      pp => ecp%oldest  !changed to => (!) -PK 7/11/06
      do while (ASSOCIATED(pp)) 
        call photosynth_cond(dtsec, pp)
        call uptake_N(dtsec, pp) !Dummy
        call litter(dtsec, pp)  !Update litter pools
        call soil_bgc(dtsec, pp)
        pp%age = pp%age + dtsec
        call summarize_patch(pp)
        pp => pp%younger  !changed to => (!) -PK 7/11/06
      end do
      call summarize_entcell(ecp)


      end subroutine ent_integrate_GISS

      !*********************************************************************

      subroutine ent_integrate(dtsec, ecp, time)
!@sum Ent biophysics/biogeochemistry
      use reproduction
      use cohorts
      use patches
      use biophysics, only : photosynth_cond
      use growthallometry, only : uptake_N
      use phenology, only : litter
      use soilbgc, only : soil_bgc
      use phenology, only : phenology_update, phenology_stats
      use canopyrad, only : recalc_radpar
      use entcells, only : summarize_entcell, entcell_print

      implicit none
      real*8 :: dtsec  !dt in seconds
      !type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      type(entcelltype) :: ecp
      !-----local--------
      type(patch),pointer :: pp
      real*8 :: time            !temp. for Phenology
      logical :: dailyupdate    !temp. for Phenology
       
      if (mod(time,86400.d0) .EQ. 0.d0) then !temporarily, later use timestruct
         dailyupdate=.true.
      else 
         dailyupdate=.false.
      end if


      !* Loop through patches
      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 
        call photosynth_cond(dtsec, pp)
        call uptake_N(dtsec, pp) !?
        !call growth(...)
        call phenology_stats(dtsec,pp,dailyupdate)
        if (dailyupdate) then
           call phenology_update(pp)
           call litter(dtsec, pp) 
        end if
!        call litter(dtsec, pp) !Moved to call after phenology_update on daily time step.
        call soil_bgc(dtsec, pp)
        pp%CO2flux = -pp%NPP + pp%Soil_resp
        pp%age = pp%age + dtsec
        call summarize_patch(pp)
        pp => pp%younger 
      end do 

      call summarize_entcell(ecp)

#ifdef DEBUG      !# DEBUG
      print *,"End of ent_biophysics"
      call entcell_print(6, ecp)
      print *,"*"
      !write(90,*) ecp%GCANOPY
      !write(90,*) ecp%GCANOPY
      !write(91,*) ecp%Ci
#endif

      end subroutine ent_integrate


      !*********************************************************************
      subroutine ent_biophysics(dtsec, ecp, config)
!@sum  Photosynthesis CO2 uptake.
!@+    If do_soilresp, then also  soil respiration for net CO2 fluxes.
      use biophysics, only : photosynth_cond
      use soilbgc, only : soil_bgc
      use patches, only : summarize_patch
      use entcells, only : summarize_entcell, entcell_print
      use phenology,only : litter !### Igor won't like this here.
      implicit none
      real*8 :: dtsec  !dt in seconds
      type(entcelltype) :: ecp
      type(ent_config) :: config
      integer :: patchnum
      !---Local--------
      type(patch),pointer :: pp

      patchnum = 0
      pp => ecp%oldest
      do while(ASSOCIATED(pp))
        patchnum = patchnum + 1
        !print*,'NEXT PATCH'
        !print*,'Calling photosynth_cond'
        call photosynth_cond(dtsec, pp)
        if (config%do_soilresp) then 
          !print*,'Calling soil_bgc'
          call soil_bgc(dtsec,pp)
          pp%CO2flux = -pp%NPP + pp%Soil_resp
          
          !*********** DIAGNOSTICS FOR PLOTTING ********************!
          write(995,*)  !Fluxes are positive up.
     &         patchnum,pp%cellptr%IPARdir + pp%cellptr%IPARdif, 
     &         pp%cellptr%coszen,
     &         pp%tallest%pft,pp%lai, pp%Tpool(CARBON,:), 
     &         pp%GPP,pp%R_auto,pp%Soil_resp,
     &         pp%NPP,pp%CO2flux,pp%GCANOPY, pp%tallest%C_lab

          !* Nancy's diagnostics *!
          write(996,*) pp%tallest%pft, pp%lai, pp%Tpool(CARBON,:),
     &         pp%GPP, pp%R_auto, pp%Soil_resp, pp%NPP, pp%CO2flux,
     &         pp%GCANOPY
          !*********************************************************!
        else 
          pp%CO2flux = UNDEF
        endif
        ! litter is updated here since it is an integration over time
        ! is do_soilresp flag ok or different flag is needed ?
        if ( config%do_soilresp ) call litter(dtsec,pp)
        pp%age = pp%age + dtsec
        !call summarize_patch(pp)

        pp => pp%younger

      end do
      call summarize_entcell(ecp)

#ifdef DEBUG      !# DEBUG
      print *,"End of ent_biophysics"
      call entcell_print(6, ecp)
      print *,"*"
      !write(90,*) ecp%GCANOPY
      !write(90,*) ecp%GCANOPY
      !write(91,*) ecp%Ci
#endif
      end subroutine ent_biophysics
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
