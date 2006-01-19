      module patches
!@sum Routines to organize patches in an entcell.

      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_patch(gp, area)
      !* Insert patch at youngest end of patch list. *!
      !* Blank patch with no cohorts.
      implicit none
      type(entcelltype) :: gp
      !type(patch), pointer :: pp
      real*8 :: area

      !* If !ASSOCIATED(gp) ERROR *!

      !* If there exist patches, then insert.  If no patches, then allocate.
      if (ASSOCIATED(gp%youngest)) then
        allocate(gp%youngest%younger)
        gp%youngest%older = gp%youngest
        gp%youngest = gp%youngest%younger
        nullify(gp%youngest%tallest)
        nullify(gp%youngest%shortest)
        allocate(gp%youngest%sumcohort)
        call init_patch(gp%youngest,gp,area)
      else !first patch on entcell
        allocate(gp%oldest)
        gp%youngest = gp%oldest
        nullify(gp%youngest%older)
        nullify(gp%youngest%younger)
        nullify(gp%youngest%tallest)
        nullify(gp%youngest%shortest)
        allocate(gp%youngest%sumcohort)
        call init_patch(gp%youngest,gp,area)
      end if

      end subroutine insert_patch

      !*********************************************************************
      subroutine delete_patch(gp, pp)
      !* Delete patch pointed to by pp
      !* NOTE:  THIS DOES NOT AUTOMATICALLY UPDATE ENTCELL SUMMMARY VALUES
      implicit none
      type(entcelltype) :: gp
      type(patch), pointer :: pp

      if (.NOT.ASSOCIATED(pp%older)) then !is oldest
        gp%oldest = gp%oldest%younger
        nullify(gp%oldest%older)
        deallocate(pp)
      else if(.NOT.ASSOCIATED(pp%younger)) then !is youngest
        gp%youngest = gp%youngest%older
        nullify(gp%youngest%younger)
        deallocate(pp)
      else !pp points somewhere in the middle of patch list
        pp%older%younger = pp%younger
        pp%younger%older = pp%older
        nullify(pp)
        deallocate(pp)
      end if

      end subroutine delete_patch

      !*********************************************************************

      subroutine summarize_patch(ttime,pp)
      !* Calculates patch-level summary values of cohort pools.
      ! * Intensive properties (e.g. geometry, LMA) are averages weighted by
      ! total number of individuals (may want to do by biomass of cohort)
      ! * Extensive properties (e.g. biomass, Ntot) are totals per m2 ground
      use cohorts, only: zero_cohort
      use canopyrad
      implicit none
      type(timestruct),pointer :: ttime
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: scop, cop
      integer :: nc  !#individuals

      scop = pp%sumcohort

      !* Zero out summary variables *!
      call zero_cohort(scop)
      scop%pft = pp%tallest%pft
      nc = 0
      scop%n = nc

      cop = pp%tallest
      do while(ASSOCIATED(cop))
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
      scop%h = scop%h/nc
      scop%crown_dx = scop%crown_dx/nc
      scop%crown_dy = scop%crown_dy/nc
      scop%dbh = scop%dbh/nc
      scop%root_d = scop%root_d/nc
      scop%clump = scop%clump/nc

      CALL  sum_roots_cohorts2patch(pp) !froot and C_froot
      
      !* BIOMASS POOLS - TOTALS
      scop%LMA = scop%LMA/nc

      !* Total individuals *!
      scop%n = nc

      !* ------- DO PATCH-LEVEL SUMMARIES OF scop ------------------------*!
      !* Structural variables
      pp%pft = pp%tallest%pft  !Patch tallest pft is considered dominant.
      pp%LAI = scop%LAI
      pp%nm = scop%nm
      pp%h = scop%h

      !* Flux variables for GCM/EWB - patch total
      call get_patchalbedo(ttime,pp)
      !pp%z0 =0.d0               !## Dummy ##!
      !pp%GCANOPY      !Calculated by biophysics
      !pp%CO2flux      !Calculate by biophysics

      !* Variables calculated by GCM/EWB - downscaled from grid cell
      !Soilmoist(:)        !Calculated by GCM/EWB
      !pp%N_deposit = 0.0  !Dummy until we have N cycle
      
      !* Variables for biophysics and biogeochemistry
      !nullify(pp%crad)     !Dummy until we have GORT clumping.

      !* Disturbance values - UPDATE WHEN WE INCLUDE DISTURBANCE
      !pp%fuel = 0.d0          !## Dummy ##!
      !pp%ignition_rate =0.d0  !## Dummy ##!
      !pp%lambda1(T_SUB) = 0.d0!## Dummy ##!
      !pp%disturbance_rate(N_DIST_TYPES)=0.d0 !## Dummy ##!

      end subroutine summarize_patch

      !*********************************************************************

      subroutine init_patch(pp,gp,area)
!@sum Initialize patch, zeroing variables, and nulling or setting pointers.
      use cohorts, only: zero_cohort
      implicit none
      type(patch),pointer :: pp
      type(entcelltype) :: gp
      real*8 :: area
      !----Local----
      integer :: n
      type(cohort),pointer :: cop

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
            pp%LAI = 0.0 !## Get GISS alai ##!
            do n=1,N_DEPTH
              pp%froot(n) = 0.0
            end do
            pp%C_froot = 0.d0

#ifdef NEWDIAG
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
            pp%REW = 0.0!## CALCULATE FROM GISS SOIL MOISTURE
            pp%soil_labile_C = 0.d0 !## Dummy ##!
            pp%soil_slow_C = 0.d0 !## Dummy ##!
            pp%soil_labile_N = 0.d0 !## Dummy ##!
            pp%soil_slow_N = 0.d0 !## Dummy ##!
            pp%mineral_N = 0.d0 !## Dummy ##!

            !* Rates of change - patch total
            pp%dadt = 0.d0 !## Dummy ##!
            pp%dpdt = 0.d0 !## Dummy ##!
            pp%dwdt = 0.d0 !## Dummy ##!
#endif

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
      implicit none
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
      do while(ASSOCIATED(cop)) 
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
      implicit none
      type(entcelltype) :: entcell


      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not needed in GISS replication test.
      
      end subroutine reorganize_patches


      !*********************************************************************


      end module patches
