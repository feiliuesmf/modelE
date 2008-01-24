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
      subroutine assign_patch(pp,
     &     Ci_ini, CNC_ini, pft, Tpool_ini)
     
      use ent_const  !for assigning Tpool -PK 12/07
      
      !Eventually may want to include all patch variables as optional vars.
      type(patch),pointer :: pp
      integer, intent(in) :: pft
      real*8 :: Ci_ini, CNC_ini
      
      !for prescribed soil C_org (or N) pools  -PK
      integer :: i, n
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS) ::
     &         Tpool_ini  !in g/m2 -PK

      pp%Ci = Ci_ini
      pp%GCANOPY = CNC_ini

#ifdef PRESCR_SOILCARB  !PK 
      do n=1,N_CASA_LAYERS 
       do i=NLIVE+1,NPOOLS
        pp%Tpool(CARBON,i,n) = Tpool_ini(pft,CARBON,i-NLIVE,n)  
       end do
      end do
#else
      pp%Tpool(:,:,:) = 0.d0  
#endif

#ifdef OFFLINE  
!      print*, 'soil C pools: ', pp%Tpool(:,CARBON,:,:)  !optional test -PK
#endif
      end subroutine assign_patch
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
      use cohorts, only: calc_CASArootfrac  !PK 7/07
      use canopyrad
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      real*8 :: nc, nsum  !density, sum
      integer :: ia  !array index
      integer :: pft
      real*8 :: fracrootCASA(N_CASA_LAYERS)  !to map fracroot to fracrootCASA -PK 7/07

      !* Zero out summary variables *!
      call zero_patch_cohortsum(pp)
      pp%Tpool(:,1:NLIVE,:) = 0.d0 !###

      if ( .not. associated(pp%tallest) ) return ! no cohorts in this patch

      do ia=1,N_COVERTYPES
        pp%LAIpft(ia) = 0.d0
      end do

      nsum = 0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
        pft = cop%pft
        
      !assign root fractions for CASA layers -PK 11/06 
      call calc_CASArootfrac(cop,fracrootCASA)
!      print *, 'from patches: fracrootCASA(:) =', fracrootCASA !***test*** -PK 11/27/06 

        !*- - - - - - COHORT SUMMARY VARIABLES - - - - - - - - - -
        nc = cop%n  
        nsum = nsum + nc !Number of individuals for summing.

        !* PFT PARAMETERS
        ! Only need to index array of pftypes.

        !* NITROGEN status and LEAF */
        pp%nm = pp%nm + cop%nm*cop%LAI !wtd average
        pp%Ntot = pp%Ntot + cop%Ntot  !Total
        pp%LAI = pp%LAI + cop%LAI  !Total
        pp%LAIpft(pft) = pp%LAIpft(pft) + cop%LAI  !Total
        pp%LMA = pp%LMA + cop%LMA * cop%LAI  !wtd avg

         !* GEOMETRY - WEIGHTED AVERAGES
        pp%h = pp%h + cop%h * nc  !wtd avg
        pp%crown_dx = pp%crown_dx + cop%crown_dx * nc !wtd avg
        pp%crown_dy = pp%crown_dy + cop%crown_dy * nc !wtd avg
        pp%clump = pp%clump + cop%clump * nc !wtd avg - ##SHOULD COME FROM RADIATIVE TRANSFER
        !* Do fracroot(:) outside this loop, below.

         !* BIOMASS POOLS - TOTALS
        pp%C_fol = pp%C_fol + cop%C_fol * nc !Total
        pp%N_fol = pp%N_fol + cop%N_fol * nc
        pp%C_w = pp%C_w + (cop%C_sw  + cop%C_hw) * nc
        pp%N_w = pp%N_w + (cop%N_sw  + cop%N_hw) * nc
        pp%C_lab = pp%C_lab + cop%C_lab * nc 
        pp%N_lab = pp%N_lab + cop%N_lab * nc 
        pp%C_froot = pp%C_froot + (cop%C_froot) * nc
        pp%N_froot = pp%N_froot + (cop%N_froot) * nc
        pp%C_root = pp%C_root + (cop%C_froot + cop%C_croot) * nc
        pp%N_root = pp%N_root + (cop%N_froot + cop%N_croot) * nc

         !* FLUXES - TOTALS
        pp%Ci = pp%Ci + cop%Ci * cop%LAI !wtd average
        pp%GCANOPY = pp%GCANOPY + cop%GCANOPY !total
        pp%GPP = pp%GPP + cop%GPP     !total
        pp%NPP = pp%NPP + cop%NPP     !total
        pp%R_auto = pp%R_auto         !total
        pp%R_root = pp%R_root         !total  !PK 5/15/07
        pp%N_up = pp%N_up + cop%N_up  !total
!        pp%C_litter = pp%C_litter + cop%C_litter
!        pp%N_litter = pp%N_litter + cop%N_litter  
!        pp%C_to_Nfix = pp%C_to_Nfix + cop%C_to_Nfix

         !* REPRODUCTION
        !*- - - - - end cohort summary variables - - - - - - - - - - - - -

        !* CASA pools * These are redundant but not 1-1 with the cohort pools.
       do ia=1,N_CASA_LAYERS  !loop over all CASA layers -PK 7/07
        if (ia.eq.1) then  !only top CASA layer has leaves & wood
         pp%Tpool(CARBON,LEAF,ia) = pp%Tpool(CARBON,LEAF,ia)
     &        + cop%C_fol * cop%n    !kg-C/m^2-ground
         pp%Tpool(CARBON,WOOD,ia) = pp%Tpool(CARBON,WOOD,ia)  !note: phenology.f uses C_croot instead of C_sw -PK 7/07
     &       + (cop%C_sw + cop%C_hw) * cop%n
        else    
         pp%Tpool(CARBON,LEAF,ia) = 0.d0  
         pp%Tpool(CARBON,WOOD,ia) = 0.d0
        end if 
         pp%Tpool(CARBON,FROOT,ia) = pp%Tpool(CARBON,FROOT,ia)
     &        + fracrootCASA(ia)*cop%C_froot * cop%n
       end do  !N_CASA_LAYERS
        !* Tpool(NITROGEN,:,:) gets updated in casa_bgfluxes.f

        cop => cop%shorter
      end do    !loop through cohorts

      !* ------- DO AVERAGES ----------------------------------------------*!

      !* CHECK IF NC = 0. If ASSOCIATE(pp%tallest) then should not be a problem.*!
!      pp%n = nsum                          !* Total individuals *!

      !* NITROGEN status and LEAF */
      pp%nm = pp%nm / pp%LAI    !wtd average
      pp%LMA = pp%LMA / pp%LAI  !wtd avg

      !* GEOMETRY - WEIGHTED AVERAGES
      pp%h = pp%h / nsum        !wtd avg
      pp%crown_dx = pp%crown_dx / nsum !wtd avg
      pp%crown_dy = pp%crown_dy / nsum !wtd avg
      pp%clump = pp%clump / nsum !wtd avg - ##SHOULD COME FROM RADIATIVE TRANSFER
      CALL  sum_roots_cohorts2patch(pp) !froot and C_froot

      !* FLUXES 
      pp%Ci = pp%Ci / pp%LAI !wtd average
!
      !* Flux variables for GCM/EWB - patch total
      !pp%z0 =0.d0               !## Dummy ##!
      !##### call get_patchalbedo(jday,pp)######
      !##### Calculate TRANS_SW ################

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

      subroutine zero_patch_cohortsum(pp)
!@sum Zero only the patch's cohort summary variables
      implicit none
      type(patch),pointer :: pp

      pp%nm = 0.d0
      pp%Ntot = 0.d0
      pp%LAI = 0.d0
      pp%LMA = 0.d0
      pp%h = 0.d0
      pp%crown_dx = 0.d0
      pp%crown_dy = 0.d0
      pp%clump = 0.d0
      pp%fracroot = 0.d0
      pp%C_fol = 0.d0
      pp%N_fol = 0.d0
      pp%C_w = 0.d0
      pp%N_w = 0.d0
      pp%C_lab = 0.d0
      pp%N_lab = 0.d0
      pp%C_froot = 0.d0
      pp%N_froot = 0.d0
      pp%C_root = 0.d0
      pp%N_root = 0.d0
      pp%Ci = 0.d0
      pp%GCANOPY = 0.d0
      pp%GPP = 0.d0
      pp%NPP = 0.d0
      pp%R_auto = 0.d0
      pp%R_root = 0.d0  !PK 5/15/07
      pp%N_up = 0.d0

      end subroutine zero_patch_cohortsum
      !*********************************************************************

      subroutine zero_patch(pp)
!@sum Initialize patch, zeroing variables.
      implicit none
      type(patch),pointer :: pp
      !----Local----

            pp%age=0.d0
            pp%area = 0.d0
            
            !* Structural variables *!
            call zero_patch_cohortsum(pp)

            !* Flux variables for GCM/EWB - patch total
            pp%z0 = 0.0d0    !## Get GISS z0 ######!
            pp%albedo = 0.0d0!## Get GISS albveg ##!
            pp%TRANS_SW = 1.0d0 !## Calculate for zero veg.##!
            pp%GCANOPY = 0.d0 !
            pp%CO2flux = 0.d0 
            pp%Ci = 0.0127D0  !Initial value not zero.
            pp%betad = 0.d0
            pp%betadl = 0.d0
            pp%Soil_resp = 0.d0

            !* Variables calculated by GCM/EWB - downscaled from grid cell
            pp%Soilmoist(:) = 0.0d0 !## Get GISS soil moisture layers ##!
!            pp%N_deposit = 0.d0

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
            pp%LAIpft(:) = 0.d0
            pp%Tpool(:,:,:) = 0.d0

            !* Soil type
            pp%soil_type = -1    ! set to undefined soil type (maybe use -1?)

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

      end subroutine zero_patch

      !*********************************************************************
      !*********************************************************************
      subroutine sum_roots_cohorts2patch(pp)
      !@sum Calculate patch-level depth- and mass-weighted average
      !@sum of fine roots depth fractions, fracroot.
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      integer :: n
      real*8 :: fracroot(N_DEPTH)
      real*8 :: frootC_total

      !Initialize
      do n=1,N_DEPTH  
        fracroot(n) = 0.0
      end do 
      frootC_total = 0.0

      cop => pp%tallest
      do while(ASSOCIATED(cop)) 
        frootC_total = frootC_total + cop%n*cop%C_froot !Mass sum per area
        !* NYK - haven't set up fracroot yet.
        do n=1,N_DEPTH
          !fracroot(n) = fracroot(n) + cop%n*cop%C_froot*cop%fracroot(n) !Mass wtd avg.
          fracroot(n) = fracroot(n) + cop%n*cop%fracroot(n) !Population wtd avg.##HACK
        end do
        cop => cop%shorter
      end do
      if (frootC_total > 0.0) then
        pp%fracroot = fracroot/frootC_total
      end if
!      pp%C_froot = frootC_total !Already calculated outside of this routine.

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

      if ( area > 1.d0 ) then 
        print *, "patch area > 1, area=",area
        call stop_model("patch_construct: area >1",255)
      endif
      ! allocate memory
      allocate( pp )
      allocate( pp%LAIpft(N_COVERTYPES) )  !remove -- have made into scalar -PK 6/28/06
      allocate( pp%betadl(N_DEPTH) )
      allocate( pp%fracroot(N_DEPTH) )
!      allocate( pp%LAIpft(N_COVERTYPES))

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
      
      end subroutine patch_construct

      !*********************************************************************

      subroutine patch_destruct(pp)
      use cohorts, only : cohort_destruct
      implicit none
      type(patch), pointer :: pp
      !---
      type(cohort),pointer :: cop, cop_tmp

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


      !*********************************************************************
      subroutine patch_print(iu,pp,prefix)
      use cohorts, only : cohort_print
      integer, intent(in) :: iu
      type(patch), intent(in) :: pp
      character*(*), optional, intent(in) :: prefix
      !---
      integer n, nc, m, i
      type(cohort),pointer :: cop
      character*8 prefix_c

      prefix_c = "        "

      write(iu,'(a,a," = ",f10.7)') prefix,"area",pp%area
      write(iu,'(a,a," = ",f10.2)') prefix,"age ",pp%age
      write(iu,'(a,a," = ",f10.7)') prefix,"GCANOPY ",
     &     pp%GCANOPY
      write(iu,'(a,a," = ",f10.7)') prefix,"Ci ",pp%Ci
      write(iu,'(a,a," = ",f10.7)') prefix,"nm ",pp%nm
      write(iu,'(a,a," = ",f10.7)') prefix,"soil H2O ",pp%Soilmoist
      write(iu,'(a,"albedo:")') prefix
      do n=1,N_BANDS
        write(iu,'(a,"      ",f10.7)') prefix,pp%albedo(n)
      enddo
      write(iu,'(a,"Tpool:")') prefix
      do m=1,PTRACE
        do n=1,NPOOLS,3
         do i=1,N_CASA_LAYERS
          write(iu,'(a,"      ",i1,"  ",e12.4,e12.4,e12.4)') prefix,m
     &         ,pp%Tpool(m,n,i),pp%Tpool(m,n+1,i),pp%Tpool(m,n+2,i)
         end do
        enddo
      enddo
      write(iu,'(a,a," = ",i7)') prefix,"soil type",pp%soil_type
      write(iu,'(a,"cohorts:")') prefix

      cop => pp%tallest
      nc = 0
      do while( associated(cop) )
        nc = nc + 1
        write( prefix_c, '(i2,"      ")' ) nc
        call cohort_print(iu, cop, prefix//prefix_c)
        cop => cop%shorter
      enddo

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
      real*8 :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)
      print*,'Tpool(CARBON,LEAF,:)',Tpool(CARBON,LEAF,:)
      print*,'Tpool(CARBON,FROOT,:)',Tpool(CARBON,FROOT,:)
      print*,'Tpool(CARBON,WOOD,:)',Tpool(CARBON,WOOD,:)
      print*,'Tpool(CARBON,SURFMET,:)',Tpool(CARBON,SURFMET,:)
      print*,'Tpool(CARBON,SURFSTR,:)',Tpool(CARBON,SURFSTR,:)
      print*,'Tpool(CARBON,SOILMET,:)',Tpool(CARBON,SOILMET,:)
      print*,'Tpool(CARBON,SOILSTR,:)',Tpool(CARBON,SOILSTR,:)
      print*,'Tpool(CARBON,CWD,:)',Tpool(CARBON,CWD,:)
      print*,'Tpool(CARBON,SOILMIC,:)',Tpool(CARBON,SURFMIC,:)
      print*,'Tpool(CARBON,SURFMIC,:)',Tpool(CARBON,SOILMIC,:)
      print*,'Tpool(CARBON,SLOW,:)',Tpool(CARBON,SLOW,:)
      print*,'Tpool(CARBON,PASSIVE,:)',Tpool(CARBON,PASSIVE,:)
      end subroutine print_Tpool

!**************************************************************************
  
      end module patches
