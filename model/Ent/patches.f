      module patches
!@sum Routines to organize patches in an entcell.

      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_patch(gp, area, soil_type)
      !* Insert patch at youngest end of patch list. *!
      !* Blank patch with no cohorts.
      implicit none
      type(entcelltype) :: gp
      real*8, intent(in) :: area
      integer, intent(in) :: soil_type
      !---
      type(patch), pointer :: pp

      !* If !ASSOCIATED(gp) ERROR *!

      ! create patch
      call patch_construct(pp, gp, area, soil_type)

      !* If there exist patches, then insert.  If no patches, then allocate.
      if (ASSOCIATED(gp%youngest)) then
        !allocate(gp%youngest%younger)
        gp%youngest%younger => pp
        pp%older => gp%youngest
        gp%youngest => pp
      else !first patch on entcell
        gp%oldest => pp
        gp%youngest => pp
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
        gp%oldest => gp%oldest%younger
        nullify(gp%oldest%older)
      else if(.NOT.ASSOCIATED(pp%younger)) then !is youngest
        gp%youngest => gp%youngest%older
        nullify(gp%youngest%younger)
      else !pp points somewhere in the middle of patch list
        pp%older%younger => pp%younger
        pp%younger%older => pp%older
        nullify(pp)
      end if
      call patch_destruct(pp)

      end subroutine delete_patch

      !*********************************************************************

      subroutine summarize_patch(pp)
      !* Calculates patch-level summary values of cohort pools.
      ! * Intensive properties (e.g. geometry, LMA) are averages weighted by
      ! total number of individuals (may want to do by biomass of cohort)
      ! * Extensive properties (e.g. biomass, Ntot) are totals per m2 ground
      use cohorts, only: zero_cohort
      use canopyrad
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: scop, cop
      real*8 :: nc  !density sum
      integer :: ia  !array index
      integer :: pft

      scop => pp%sumcohort

      !* Zero out summary variables *!
      call zero_cohort(scop)
      nc = 0.d0
      scop%n = nc
      pp%Tpool(:,1:NLIVE) = 0.d0 !###

      if ( .not. associated(pp%tallest) ) return ! no cohorts in this patch

      do ia=1,N_COVERTYPES
        pp%LAI(ia) = 0.d0
      end do

      scop%pft = pp%tallest%pft

      cop => pp%tallest
      do while(ASSOCIATED(cop))
        pft = cop%pft
        nc = nc + cop%n  !Number of individuals summed for wtd avg.
        !* PFT PARAMETERS
        ! Only need to index array of pftypes.

        !* NITROGEN status - CANOPY*/
        scop%nm = scop%nm + cop%nm*cop%n  !wtd avg
        scop%Ntot = scop%Ntot + cop%Ntot  !Total

        scop%LAI = scop%LAI + cop%LAI  !Total
        pp%LAI(pft) = pp%LAI(pft) + cop%LAI
        pp%Tpool(CARBON,LEAF) = pp%Tpool(CARBON,LEAF)
     &       + cop%C_fol * cop%n    !kg-C/m^2-ground
        pp%Tpool(CARBON,FROOT) = pp%Tpool(CARBON,FROOT)
     &       + cop%C_froot * cop%n  
        pp%Tpool(CARBON,WOOD) = pp%Tpool(CARBON,WOOD)
     &       + (cop%C_hw + cop%C_croot) * cop%n
        !* Tpool(NITROGEN,:,:) gets updated in casa_bgfluxes.f
        pp%LAI(cop%pft) = pp%LAI(cop%pft) + cop%LAI

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
        !scop%LMA = scop%LMA + cop%LMA*cop%n !wtd avg
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
!        scop%C_litter = scop%C_litter + cop%C_litter
!        scop%N_litter = scop%N_litter + cop%N_litter  
        scop%C_to_Nfix = scop%C_to_Nfix + cop%C_to_Nfix

         !* REPRODUCTION
         !scop%------        ! ASK PAUL

        cop => cop%shorter
      end do

      !* ------- DO AVERAGES ----------------------------------------------*!

      !** ## CHECK IF NC = 0.D0 ## **!

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
      !scop%LMA = scop%LMA/nc
      !pp%Tpool(CARBON,LEAF) = pp%Tpool(CARBON,LEAF)/pp%area !m-2 basis - NYK
      !pp%Tpool(CARBON,FROOT) = pp%Tpool(CARBON,FROOT)/pp%area
      !pp%Tpool(CARBON,WOOD) = pp%Tpool(CARBON,WOOD)/pp%area
      
      !* Total individuals *!
      scop%n = nc


      !* ------- DO PATCH-LEVEL SUMMARIES OF scop ------------------------*!
      !* Structural variables
      pp%nm = scop%nm

      !* Flux variables for GCM/EWB - patch total
      !pp%z0 =0.d0               !## Dummy ##!
      !##### call get_patchalbedo(jday,pp)######
      !##### Calculate TRANS_SW ################
      !pp%GCANOPY      !Calculated by biophysics

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

      subroutine zero_patch(pp)
!@sum Initialize patch, zeroing variables.
      implicit none
      type(patch),pointer :: pp
      !----Local----

            pp%age=0.d0
            pp%area = 0.d0

            !* Structural variables *!
            pp%nm = 0.0d0

            !* Flux variables for GCM/EWB - patch total
            pp%z0 = 0.0d0    !## Get GISS z0 ######!
            pp%albedo = 0.0d0!## Get GISS albveg ##!
            pp%TRANS_SW = 0.0d0 !## Calculate ##!
            pp%GCANOPY = 0.d0 !Will be updated in biophysics.f
            pp%CO2flux = 0.d0 
            pp%Ci = 0.0127D0  !Initial value not zero.
            pp%betad = 0.d0
            pp%betadl = 0.d0
            pp%Soil_resp = 0.d0
            pp%R_can = 0.d0  !##Temporary hack. Belongs in sumcohort.
            pp%R_root = 0.d0  !##Temporary hack. Belongs in sumcohort.

            !* Variables calculated by GCM/EWB - downscaled from grid cell
            pp%Soilmoist = 0.0d0 !## Get GISS soil moisture layers ##!
            pp%N_deposit = 0.d0

            !* Variables for biophysics and biogeochemistry
            !pp%crad%##### = !## Get GORT canopy radiation params ##!
            ! initializing it at least to something 
            nullify(pp%crad%heights)
            
            nullify(pp%crad%lai)
            pp%crad%gortclump = 0.d0

            !* Disturbance values
            pp%fuel = 0.d0      !## Dummy ##!
            pp%ignition_rate =0.d0 !## Dummy ##!
            pp%lambda1 = 0.d0 !## Dummy ##!
            pp%disturbance_rate = 0.d0 !## Dummy ##!

            !* DIAGNOSTIC SUMMARIES
            !* Biomass pools - patch total
            pp%LAI(:) = 0.d0
            pp%C_froot = 0.d0
            pp%Tpool(:,:) = 0.d0

            !* Soil type
            pp%soil_type = UNDEF    ! set to undefined soil type (maybe use -1?)

            !* Soil vars for CASA -PK
            pp%Soilmoist = 0.d0

#ifdef NEWDIAG
            pp%plant_ag_Cp = 0.d0 !## Dummy ##!
            pp%plant_bg_Cp = 0.d0 !## Dummy ##!
            pp%plant_ag_Np = 0.d0 !## Dummy ##!
            pp%plant_bg_Np = 0.d0 !## Dummy ##!

            !* Biomass pools - by pft
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
!            pp%Soil_resp = 0.d0 !## Dummy ##! !moved to fluxes -PK 6/14/06

      end subroutine zero_patch

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

      cop => pp%tallest
      do while(ASSOCIATED(cop)) 
        !frootC_total = frootC_total + cop%n*cop%C_froot !Mass wtd avg.
        frootC_total = frootC_total + cop%n !##HACK UNTIL CAN DO frootC wtd avg
        do n=1,N_DEPTH
          !froot(n) = froot(n) + cop%n*cop%C_froot*cop%froot(n) !Mass wtd avg.
          froot(n) = froot(n) + cop%n*cop%froot(n) !Population wtd avg.##HACK
        end do
        cop => cop%shorter
      end do
      if (frootC_total > 0.0) then
        pp%sumcohort%froot = froot/frootC_total
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

      subroutine patch_construct(pp, parent_entcell, area, soil_type)
      use cohorts, only : cohort_construct
      implicit none
      type(patch), pointer :: pp
      type(entcelltype), intent(in), target :: parent_entcell
      real*8, intent(in) :: area
      integer, intent(in) :: soil_type
      !---

      if ( area > 1.d0 ) call stop_model("patch_construct: area >1",255)
      ! allocate memory
      allocate( pp )
!      allocate( pp%Soilmoist(N_DEPTH) )  !remove -- have made into scalar -PK 6/28/06
      allocate( pp%betadl(N_DEPTH) )
!      allocate( pp%LAI(N_COVERTYPES))

      ! set pointers
      pp%cellptr => parent_entcell
      nullify( pp%older )
      nullify( pp%younger )
      nullify( pp%tallest )
      nullify( pp%shortest )

      ! set variables
      call zero_patch(pp)

      pp%area = area
      pp%soil_type = soil_type
      
      ! create sumcohort
      call cohort_construct(pp%sumcohort, pp)

      end subroutine patch_construct

      !*********************************************************************

      subroutine patch_destruct(pp)
      use cohorts, only : cohort_destruct
      implicit none
      type(patch), pointer :: pp
      !---
      type(cohort),pointer :: cop, cop_tmp

      ! destroy sumcohort
      call cohort_destruct(pp%sumcohort)

      ! destroy other cohorts
      cop => pp%tallest
      do while( associated(cop) )
        cop_tmp => cop
        cop => cop%shorter
        call cohort_destruct(cop_tmp)
      enddo
      !nullify(pp%tallest)
      !nullify(pp%shortest)

      ! deallocate memory
!      deallocate( pp%Soilmoist )
      deallocate( pp )
      nullify( pp )

      end subroutine patch_destruct


      subroutine patch_print(pp,prefix)
      use cohorts, only : cohort_print
      type(patch) :: pp
      character*(*), optional :: prefix
      integer n, nc
      type(cohort),pointer :: cop
      character*8 prefix_c

      prefix_c = "        "

      print '(a,a," = ",f10.7)',prefix,"area",pp%area
      print '(a,a," = ",f10.2)',prefix,"age ",pp%age
      print '(a,a," = ",f10.7)',prefix,"GCANOPY ",pp%GCANOPY
      print '(a,a," = ",f10.7)',prefix,"Ci ",pp%Ci
      print '(a,a," = ",f10.7)',prefix,"nm ",pp%nm
      print '(a,"soil moisture:")',prefix
      do n=1,N_BANDS
        print '(a,"      ",f10.7)',prefix,pp%albedo(n)
      enddo
      print '(a,"      ",f10.7)',prefix,pp%Soilmoist
!      do n=1,N_DEPTH
!        print '(a,"      ",f10.7)',prefix,pp%Soilmoist(n)
!      enddo
      print '(a,a," = ",i7)',prefix,"soil type",pp%soil_type
      print '(a,"cohorts:")',prefix

      cop => pp%tallest
      nc = 0
      do while( associated(cop) )
        nc = nc + 1
        write( prefix_c, '(i2,"      ")' ) nc
        call cohort_print(cop,prefix//prefix_c)
        cop => cop%shorter
      enddo

      print '(a,"sumcohort:")',prefix
      call cohort_print(pp%sumcohort,prefix//"s       ")

      end subroutine patch_print


      subroutine patch_extract_pfts(pp, vdata)
      use cohorts, only : cohort_print
      type(patch) :: pp
      real*8 :: vdata(:)
      !---
      type(cohort),pointer :: cop
      real*8 :: density

      vdata(:) = 0.d0
      density = 0.d0
      cop => pp%tallest
      do while( associated(cop) )
        vdata( cop%pft ) = vdata( cop%pft ) + cop%n
        density = density + cop%n
        cop => cop%shorter
      enddo

      if ( density>0.d0 ) vdata(:) = vdata(:)/density

      end subroutine patch_extract_pfts
     

!**************************************************************************

      subroutine print_Tpool(Tpool)
      real*8 :: Tpool(PTRACE,NPOOLS)
      print*,'Tpool(CARBON,LEAF)',Tpool(CARBON,LEAF)
      print*,'Tpool(CARBON,FROOT)',Tpool(CARBON,FROOT)
      print*,'Tpool(CARBON,WOOD)',Tpool(CARBON,WOOD)
      print*,'Tpool(CARBON,SURFMET)',Tpool(CARBON,SURFMET)
      print*,'Tpool(CARBON,SURFSTR)',Tpool(CARBON,SURFSTR)
      print*,'Tpool(CARBON,SOILMET)',Tpool(CARBON,SOILMET)
      print*,'Tpool(CARBON,SOILSTR)',Tpool(CARBON,SOILSTR)
      print*,'Tpool(CARBON,CWD)',Tpool(CARBON,CWD)
      print*,'Tpool(CARBON,SOILMIC)',Tpool(CARBON,SOILMIC)
      print*,'Tpool(CARBON,SLOW)',Tpool(CARBON,SLOW)
      print*,'Tpool(CARBON,PASSIVE)',Tpool(CARBON,PASSIVE)
      end subroutine print_Tpool

!**************************************************************************
  
      end module patches
