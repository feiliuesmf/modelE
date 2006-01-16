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
        call nullify(gp%youngest%tallest)
        call nullify(gp%youngeset%shortest)
        call allocate(gp%youngest%sumcohort)
        call allocate(gp%youngest%avgcohort)
        call init_patch(gp%youngest,gp,area)
      else !first patch on entcell
        call allocate(gp%oldest)
        gp%youngest = gp%oldest
        call nullify(gp%youngest%older)
        call nullify(gp%youngest%younger)
        call nullify(gp%youngest%tallest)
        call nullify(gp%youngeset%shortest)
        call allocate(gp%youngest%sumcohort)
        call allocate(gp%youngest%avgcohort)
        call init_patch(gp%youngest,gp,area)
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
      subroutine summarize_patch_OLD(pp)
      !* Recalculates patch-level summary values, parsing through cohorts.
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      
      !If !allocated(pp) ERROR
      
      pp%age = pp%age + !*TIME SINCE LAST UPDATE *!
      !pp%area=area
      !pp%cellptr = gp

      !* Flux variables for GCM/EWB - patch total
      !pp%albedo =               !## Get GISS albveg ##!
      !pp%z0 =                   !## Get GISS z0 ######!
      !pp%GCANOPY = 0.d0         !Will be updated in biophysics.f
      !pp%CO2flux = 0.d0         !Will be updated in biophysics.f

      !* Zero out summary variables *!
      do n=1,N_PFT
        pp%albedo(N_BANDS) = 0.0
        pp%LAI(n) = 0.0
        pp%froot(n) = 0.0
        pp%plant_ag_Cp(n) = 0.0
        pp%plant_bg_C(n) = 0.d0 !## Dummy ##!
        pp%plant_ag_N(n) = 0.d0 !## Dummy ##!
        pp%plant_bg_N(n) = 0.d0 !## Dummy ##!

      end do

      cop = pp%tallest
      do while(allocated(cop))
        !Sum up biomass, carbon, nitrogen pools
        !---------------------------------------------------------------
        !              FILL IN CODE                                    
        !---------------------------------------------------------------

        !* Variables for biophysics and biogeochemistry
        !pp%crad%##### = !## Get GORT canopy radiation params ##!
        !pp%albedo(n) =

        !* Disturbance values
        !pp%fuel = 0.d0          !## Dummy ##!
        !pp%ignition_rate =0.d0  !## Dummy ##!
        !pp%lambda1 = 0.d0       !## Dummy ##!
        !pp%disturbance_rate = 0.d0 !## Dummy ##!

        !* DIAGNOSTIC SUMMARIES
        !* Biomass pools - patch total *!
        pp%LAI_p = pp%LAI_p + cop%LAI
        call sum_roots_cohorts2patch(pp) !froot and C_froot

        sC_froot = sC_froot + cop%C_froot
        do n=1,N_DEPTH
          pp%C_froot(n) = froot(n) + cop%C_froot*cop%froot(n)
        end do

        !* NOTE:  Labile stored carbon can be both above- and below-ground.
        !*        May want to separate this pool, but for now is added to ag.
        pp%plant_ag_Cp = pp%plant_ag_Cp + 
                         cop%C_fol + cop%C_sw + cop%C_hw + cop%C_lab
        pp%plant_bg_Cp = pp%plant_bg_Cp + cop%C_froot + cop%C_croot
        pp%plant_ag_Np = pp%plant_ag_Np + 
                         cop%N_fol + cop%N_sw + cop%N_hw + cop%N_lab
        pp%plant_bg_Np = pp%plant_bg_Np + cop%N_froot + cop%N_croot
        
        !* Biomass pools - by pft *!
        pp%LAI(cop%pft) = pp%LAI(cop%pft) + cop%LAI
        
        if (cop%pft == pnum) then 
!        pp%plant_ag_C(pnum) = pp%plant_ag_C(pnum) + 
!     &       cop%C_fol + cop%C_sw + cop%C_hw
!        pp%plant_bg_C(pnum) = pp%plant_bg_C(pnum) +
!     &       cop%C_froot + cop %C_croot
!        pp%plant_ag_N(pnum) = pp%plant_ag_N(pnum) +
!     &       cop%N_fol + cop%N_sw + cop%N_hw
!        pp%plant_bg_N(pnum) = pp%plant_bg_N(pnum) +
!     &       cop%N_froot + cop%N_croot
        end if

        !* Soil pools - patch total *!
!        pp%REW =                !## GET FROM GISS SOIL MOISTURE
!        pp%soil_labile_C = pp%soil_labile_C +
!        pp%soil_slow_C = pp%soil_slow_C +
!        pp%soil_labile_N = pp%soil_labile_N +
!        pp%soil_slow_N = pp%soil_slow_N +
!        pp%mineral_N = pp%mineral_N +

        !* Rates of change - patch total
!        pp%dadt = 0.d0          !## Dummy ##!
!        pp%dpdt = 0.d0          !## Dummy ##!
!        pp%dwdt = 0.d0          !## Dummy ##!
        
        !* Activity diagnostics - can be summed by month, year, etc.
        pp%GPP = pp%GPP + cop%GPP
!        pp%NPP = pp%NPP + cop%NPP
!        pp%Soil_resp = 0.d0     !## Dummy ##!
        
        !Next
        cop = cop%shorter
      end do

      end subroutine summarize_patch_OLD

      !*********************************************************************
      subroutine summarize_patch(pp)
      !* Calculates patch-level summary values of cohort pools.
      ! * Intensive properties (e.g. geometry, LMA) are averages weighted by
      ! total number of individuals (may want to do by biomass of cohort)
      ! * Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: scop, cop
      integer :: nc  !#individuals

      scop = pp%sumcohort

      !* Zero out summary variables *!
      call zero_cohort(scop)
      scop%pnum = pp%tallest%pft
      nc = 0
      scop%n = nc

      cop = pp%tallest
      do while(allocated(cop))
        nc = nc + cop%n  !Number of individuals summed for wtd avg.
        !* PFT PARAMETERS
        ! Only need to index array of pftypes.
         
        !* NITROGEN status - CANOPY*/
        scop%nm = scop%nm + cop%nm*cop%n  !wtd avg
        scop%Ntot = scop%Ntot + cop%Ntot  !Total
        scop%LAI = scop%LAI + cop%LAI  !Total

         !* ALL QUANTITIES BELOW ARE FOR AN INDIVIDUAL *!

         !* GEOMETRY - WEIGHTED AVERAGES
        scop%h = scop%h + cop%h*cop%n  !wtd avg
        scop%crown_dx = scop%crown_dx + cop%crown_dx*cop%n !wtd avg
        scop%crown_dy = scop%crown_dy + cop%crown_dy*cop%n !wtd avg
        scop%dbh = scop%dbh + cop%dbh*cop%n !wtd avg
        scop%root_d = scop%root_d + cop%root_d*cop%n !wtd avg
        scop%clump = scop%clump + cop%clump*cop%n !wtd avg
        !* Do froot(:) outside this loop, below.

         !* BIOMASS POOLS - TOTALS
        scop%LMA = scop%LMA + cop%LMA*cop%n !wtd avg
        scop%C_fol = scop%C_fol + cop%C_fol !Total
        scop%N_fol = scop%N_fol + cop%N_fol
        scop%C_sw = scop%C_sw + cop%C_sw    
        scop%N_sw = scop%N_sw + cop%N_sw    
        scop%C_hw = scop%C_hw + cop%C_hw    
        scop%N_hw = scop%N_hw + cop%N_hw    
        scop%C_lab = scop%C_lab + cop%C_lab   
        scop%N_lab = scop%N_lab + cop%N_lab   
        scop%C_froot = scop%C_froot + cop%C_froot 
        scop%N_froot = scop%N_froot + cop%N_froot 
        scop%C_croot = scop%C_croot + cop%C_croot 
        scop%N_croot = scop%N_croot + cop%N_croot 

         !* FLUXES - TOTALS
        scop%gcanopy = scop%gcanopy + cop%gcanopy 
        scop%GPP = scop%GPP + cop%GPP     
        scop%NPP = scop%NPP + cop%NPP
        scop%R_growth = scop%R_growth + cop%R_growth
        scop%R_maint = scop%R_maint + cop%R_maint 
        scop%N_up = scop%N_up + cop%N_up
        scop%C_litter = scop%C_litter + cop%C_litter
        scop%N_litter = scop%N_litter + cop%N_litter  
        scop%C_to_Nfix = scop%C_to_Nfix + cop%C_to_Nfix 

         !* REPRODUCTION
         !scop%------        ! ASK PAUL
        
        cop = cop%shorter
      end do

      !* ------- DO AVERAGES ----------------------------------------------*!
      !* NITROGEN status */
      scop%nm = scop%nm/nc

      !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
      scop%h = sumh/nc
      scop%crown_dx = acop%crown_dx/nc
      scop%crown_dy = acop%crown_dy/nc
      scop%dbh = cop%dbh/nc
      scop%root_d = scop%root_d/nc
      scop%clump = scop%clump/nc
      CALL  sum_roots_cohorts2patch(pp) !froot and C_froot
      
      !* BIOMASS POOLS - TOTALS
      scop%LMA = scop%LMA/nc

      !* Total individuals *!
      scop%n = nc
      end subroutine summarize_patch
      !*********************************************************************

      subroutine summarize_patch_OLD(pp)
      !* Calculates patch-level summary values of cohort pools.
      !* - Put sums in pp%sumcohort. Where sum is not meaningful, 
      !* such as height, the max value is taken.
      !* - Put averages in pp%avgcohort.

      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: scop, acop, cop
      integer :: n

      !If !allocated(pp) ERROR
      
      !* Flux variables for GCM/EWB - patch total
      !pp%albedo =               !## Get GISS albveg ##!
      !pp%z0 =                   !## Get GISS z0 ######!
      !pp%GCANOPY = 0.d0         !Will be updated in biophysics.f
      !pp%CO2flux = 0.d0         !Will be updated in biophysics.f

      scop = pp%avgcohort
      acop = pp%sumcohort

      !* Zero out summary variables *!
      call zero_cohort(scop)
      call zero_cohort(acop)

      cop = pp%tallest
      do while(allocated(cop))
        !Sum up biomass, carbon, nitrogen pools
        !---------------------------------------------------------------
        !              FILL IN CODE                                    
        !---------------------------------------------------------------

        !* Variables for biophysics and biogeochemistry
        !pp%crad%##### = !## Get GORT canopy radiation params ##!
        !pp%albedo(n) =

        !* Disturbance values
        !pp%fuel = 0.d0          !## Dummy ##!
        !pp%ignition_rate =0.d0  !## Dummy ##!
        !pp%lambda1 = 0.d0       !## Dummy ##!
        !pp%disturbance_rate = 0.d0 !## Dummy ##!

        !* DIAGNOSTIC SUMMARIES
        !* Biomass pools - patch total *!
        scop%LAI = scop%LAI + cop%LAI
        call sum_roots_cohorts2patch(pp) !froot and C_froot

        sC_froot = sC_froot + cop%C_froot
        do n=1,N_DEPTH
          pp%C_froot(n) = froot(n) + cop%C_froot*cop%froot(n)
        end do

        !* NOTE:  Labile stored carbon can be both above- and below-ground.
        !*        May want to separate this pool, but for now is added to ag.
        pp%plant_ag_Cp = pp%plant_ag_Cp + cop%n*(
     &       cop%C_fol + cop%C_sw + cop%C_hw + cop%C_lab)
        pp%plant_bg_Cp = pp%plant_bg_Cp + cop%n*(
     &       cop%C_froot + cop%C_croot
        pp%plant_ag_Np = pp%plant_ag_Np + cop%n*(
     &       cop%N_fol + cop%N_sw + cop%N_hw + cop%N_lab
        pp%plant_bg_Np = pp%plant_bg_Np + cop%n*(
     &       cop%N_froot + cop%N_croot
        
        !* Biomass pools - by pft *!
        pp%LAI(cop%pft) = pp%LAI(cop%pft) + cop%LAI
        
        if (cop%pft == pnum) then 
!        pp%plant_ag_C(pnum) = pp%plant_ag_C(pnum) + 
!     &       cop%C_fol + cop%C_sw + cop%C_hw
!        pp%plant_bg_C(pnum) = pp%plant_bg_C(pnum) +
!     &       cop%C_froot + cop %C_croot
!        pp%plant_ag_N(pnum) = pp%plant_ag_N(pnum) +
!     &       cop%N_fol + cop%N_sw + cop%N_hw
!        pp%plant_bg_N(pnum) = pp%plant_bg_N(pnum) +
!     &       cop%N_froot + cop%N_croot
        end if

        !* Soil pools - patch total *!
!        pp%REW =                !## GET FROM GISS SOIL MOISTURE
!        pp%soil_labile_C = pp%soil_labile_C +
!        pp%soil_slow_C = pp%soil_slow_C +
!        pp%soil_labile_N = pp%soil_labile_N +
!        pp%soil_slow_N = pp%soil_slow_N +
!        pp%mineral_N = pp%mineral_N +

        !* Rates of change - patch total
!        pp%dadt = 0.d0          !## Dummy ##!
!        pp%dpdt = 0.d0          !## Dummy ##!
!        pp%dwdt = 0.d0          !## Dummy ##!
        
        !* Activity diagnostics - can be summed by month, year, etc.
        pp%GPP = pp%GPP + cop%GPP
!        pp%NPP = pp%NPP + cop%NPP
!        pp%Soil_resp = 0.d0     !## Dummy ##!
        
        !Next
        cop = cop%shorter
      end do

      end subroutine summarize_patch_OLD

      !*********************************************************************

      subroutine init_patch(pp,gp,area)
      type(patch),pointer :: pp
      type(entcelltype),pointer :: gp
      real*8 :: area
      integer :: n

            pp%age=0.d0
            pp%area=area
            pp%cellptr = gp

            !* sumcohort *!
            call zero_cohort(pp%sumcohort)
            cop = pp%sumcohort
            cop%pptr = pp
            cop%pft = 0
            cop%n = 0
            cop%nm = 0.0
            cop%Ntot = 0.0

            !* avgcohort *!
            call zero_cohort(pp%avgcohort)
            cop = pp%avgcohort
            cop%pptr = pp
            cop%pft = 0
            cop%n = 0
            cop%nm = 0.0
            cop%Ntot = 0.0

            !* Flux variables for GCM/EWB - patch total
            pp%albedo = 0.0!## Get GISS albveg ##!
            pp%z0 = 0.0    !## Get GISS z0 ######!
            pp%GCANOPY = 0.d0 !Will be updated in biophysics.f
            pp%CO2flux = 0.d0 

            !* Variables calculated by GCM/EWB - downscaled from grid cell
            do n=1,N_DEPTH
              pp%Soilmoist = 0.0   !## Get GISS soil moisture layers ##!
            end do
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
            do n=1,N_DEPTH
              pp%froot(n) = 0.0
            end do

            pp%plant_ag_Cp = 0.d0 !## Dummy ##!
            pp%plant_bg_Cp = 0.d0 !## Dummy ##!
            pp%plant_ag_Np = 0.d0 !## Dummy ##!
            pp%plant_bg_Np = 0.d0 !## Dummy ##!

            !* Biomass pools - by pft
            pp%LAI(cop%pft) = pp%LAI(cop%pft) + cop%LAI
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

      end subroutine init_patch

      !*********************************************************************
      !*********************************************************************
      subroutine sum_roots_cohorts2patch(pp)
      !@sum Calculate patch-level depth- and mass-weighted average
      !@sum of fine roots.
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      integer :: n
      real*8 :: froot(N_DEPTH)
      real*8 :: frootC_total

      do n=1,N_DEPTH  !Initialize
        froot(n) = 0.0
      end do 
      frootC_total = 0.0

      cop = pp%tallest
      do while(allocated(cop)) 
        frootC_total = frootC_total + cop%n*cop%C_froot
        do n=1,N_DEPTH
          froot(n) = froot(n) + cop%n*cop%C_froot*cop%froot(n)
        end do
        cop = cop%shorter
      end do
      if (frootC_total > 0.0) then
        pp%froot = froot/frootC_total
      end if
      pp%C_froot = frootC_total

      end subroutine sum_roots_cohorts2patch
      !*********************************************************************

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
