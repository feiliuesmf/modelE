      module patches
!@sum Routines to organize patches in an entcell.

      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_patch(gp, area)
      !* Insert patch at youngest end of patch list. *!
      !* Blank patch with no cohorts.
      type(entcelltype), pointer :: gp
      type(patch), pointer :: pp
      real*8 :: area

      !* If !allocated(gp) ERROR *!

      !* If there exist patches, then insert.  If no patches, then allocate.
      if (allocated(gp%youngest)) then
        call allocate(gp%youngest%younger)
        gp%youngest%older = gp%youngest
        gp%youngest = gp%youngest%younger
        call init_patch_defaults(gp%youngest,gp,area)
      else !first patch on entcell
        call allocate(gp%oldest)
        gp%youngest = gp%oldest
        call nullify(gp%oldest%older)
        call nullify(gp%oldest%younger)
        call init_patch_defaults(gp%oldest, gp, area)
      end if

      end subroutine insert_patch

      !*********************************************************************
      subroutine delete_patch(gp, pp)
      !* Delete patch pointed to by pp
      !* NOTE:  THIS DOES NOT AUTOMATICALLY UPDATE ENTCELL SUMMMARY VALUES
      type(entcelltype), pointer :: gp
      type(patch), pointer :: pp

      if (pp.eq.gp%oldest) then
        gp%oldest = gp%oldest%younger
        call nullify(gp%oldest%older)
        call deallocate(pp)
      else if(pp.eq.gp%youngest) then
        gp%youngest = gp%youngest%older
        call nullify(gp%youngest%younger)
        call deallocate(pp)
      else !pp points somewhere in the middle of patch list
        pp%older%younger = pp%younger
        pp%%younger%older = pp%older
        call nullify(pp)
        call deallocate(pp)
      end if

      end subroutine delete_patch

      !*********************************************************************
      subroutine summarize_patch(pp)
      !* Recalculates patch-level summary values, parsing through cohorts.
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      
      !If !allocated(pp) ERROR
      
      pp%age = pp%age + !*TIME SINCE LAST UPDATE *!
      !pp%area=area
      !pp%cellptr = gp
                                !* Flux variables for GCM/EWB - patch total
      !pp%albedo =               !## Get GISS albveg ##!
      pp%z0 =                   !## Get GISS z0 ######!
      pp%GCANOPY = 0.d0         !Will be updated in biophysics.f
      pp%CO2flux = 0.d0 

      cop = pp%tallest
      do while(allocated(cop))
        !Sum up biomass, carbon, nitrogen pools
        !---------------------------------------------------------------
        !              FILL IN CODE                                    
        !---------------------------------------------------------------

        !* Variables calculated by GCM/EWB - downscaled from grid cell
        pp%Soilmoist =          !## Get GISS soil moisture layers ##!
        pp%N_deposit = 0.d0

        !* Variables for biophysics and biogeochemistry
        !pp%crad%##### = !## Get GORT canopy radiation params ##!

        !* Disturbance values
        pp%fuel = 0.d0          !## Dummy ##!
        pp%ignition_rate =0.d0  !## Dummy ##!
        pp%lambda1 = 0.d0       !## Dummy ##!
        pp%disturbance_rate = 0.d0 !## Dummy ##!

        !* DIAGNOSTIC SUMMARIES
        !* Biomass pools - patch total
        pp%LAI_p = pp%LAI_p + cop%LAI
        pp%plant_ag_Cp = 0.d0   !## Dummy ##!
        pp%plant_bg_Cp = 0.d0   !## Dummy ##!
        pp%plant_ag_Np = 0.d0   !## Dummy ##!
        pp%plant_bg_Np = 0.d0   !## Dummy ##!
        
        !* Biomass pools - by pft
        pp%LAI(pnum) = alai(longi,latj,pnum)
        pp%plant_ag_C(pnum) = 0.d0 !## Dummy ##!
        pp%plant_bg_C(pnum) = 0.d0 !## Dummy ##!
        pp%plant_ag_N(pnum) = 0.d0 !## Dummy ##!
        pp%plant_bg_N(pnum) = 0.d0 !## Dummy ##!
        
        !* Soil pools - patch total
        pp%REW =                !## GET FROM GISS SOIL MOISTURE
        pp%soil_labile_C = 0.d0 !## Dummy ##!
        pp%soil_slow_C = 0.d0   !## Dummy ##!
        pp%soil_labile_N = 0.d0 !## Dummy ##!
        pp%soil_slow_N = 0.d0   !## Dummy ##!
        pp%mineral_N = 0.d0     !## Dummy ##!

        !* Rates of change - patch total
        pp%dadt = 0.d0          !## Dummy ##!
        pp%dpdt = 0.d0          !## Dummy ##!
        pp%dwdt = 0.d0          !## Dummy ##!
        
        !* Activity diagnostics - can be summed by month, year, etc.
        pp%GPP = 0.d0           !## Dummy ##!
        pp%NPP = 0.d0           !## Dummy ##!
        pp%Soil_resp = 0.d0     !## Dummy ##!
        
        !Next
        cop = cop%shorter
      end do

      end subroutine summarize_patch

      !*********************************************************************

      subroutine init_patch_defaults(pp,gp,area)
      type(patch),pointer :: pp
      type(entcelltype),pointer :: gp
      real*8 :: area

            pp%age=0.d0
            pp%area=area
            pp%cellptr = gp

            !* Flux variables for GCM/EWB - patch total
            pp%albedo = !## Get GISS albveg ##!
            pp%z0 =     !## Get GISS z0 ######!
            pp%GCANOPY = 0.d0 !Will be updated in biophysics.f
            pp%CO2flux = 0.d0 

            !* Variables calculated by GCM/EWB - downscaled from grid cell
            pp%Soilmoist = !## Get GISS soil moisture layers ##!
            pp%N_deposit = 0.d0

            !* Variables for biophysics and biogeochemistry
            !pp%crad%##### = !## Get GORT canopy radiation params ##!

            !* Disturbance values
            pp%fuel = 0.d0      !## Dummy ##!
            pp%ignition_rate =0.d0 !## Dummy ##!
            pp%lambda1 = 0.d0 !## Dummy ##!
            pp%disturbance_rate = 0.d0 !## Dummy ##!

            !* DIAGNOSTIC SUMMARIES
            !* Biomass pools - patch total
            pp%LAI_p = !## Get GISS alai ##!
            pp%plant_ag_Cp = 0.d0 !## Dummy ##!
            pp%plant_bg_Cp = 0.d0 !## Dummy ##!
            pp%plant_ag_Np = 0.d0 !## Dummy ##!
            pp%plant_bg_Np = 0.d0 !## Dummy ##!

            !* Biomass pools - by pft
            pp%LAI(pnum) = alai(longi,latj,pnum)
            pp%plant_ag_C(pnum) = 0.d0 !## Dummy ##!
            pp%plant_bg_C(pnum) = 0.d0 !## Dummy ##!
            pp%plant_ag_N(pnum) = 0.d0 !## Dummy ##!
            pp%plant_bg_N(pnum) = 0.d0 !## Dummy ##!

            !* Soil pools - patch total
            pp%REW = !## CALCULATE FROM GISS SOIL MOISTURE
            pp%soil_labile_C = 0.d0 !## Dummy ##!
            pp%soil_slow_C = 0.d0 !## Dummy ##!
            pp%soil_labile_N = 0.d0 !## Dummy ##!
            pp%soil_slow_N = 0.d0 !## Dummy ##!
            pp%mineral_N = 0.d0 !## Dummy ##!

            !* Rates of change - patch total
            pp%dadt = 0.d0 !## Dummy ##!
            pp%dpdt = 0.d0 !## Dummy ##!
            pp%dwdt = 0.d0 !## Dummy ##!

            !* Activity diagnostics - can be summed by month, year, etc.
            pp%GPP = 0.d0 !## Dummy ##!
            pp%NPP = 0.d0 !## Dummy ##!
            pp%Soil_resp = 0.d0 !## Dummy ##!

      end subroutine init_patch_defaults

      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      
      subroutine reorganize_patches(entcell)
      type(entcelltype) :: entcell


      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not needed in GISS replication test.
      
      end subroutine reorganize_patches


      !*********************************************************************


      end module patches
