      module entcells
!@sum entcells contains routines to perform summaries at the ent grid cell
!@sum level.


      !Ent MODULES TO USE
      use ent_const
      use ent_types
      use patches
      use cohorts

      implicit none
      private
      save

      public zero_entcell, summarize_entcell
      public init_simple_entcell, entcell_construct, entcell_destruct


      contains
!**************************************************************************

      subroutine zero_entcell(ecp)
      !@sum Zeros import/export variables of entdata type
      type(entcelltype) :: ecp
      integer :: n

      ecp%area = 0.0

      !Cell-level summary values - PHYSICAL
      !EXPORT - from radiative transfer
      do n = 1,N_BANDS
        ecp%albedo(n) = 0.0      !Albedo may be in bands or hyperspectral
      end do
      ecp%TRANS_SW = 0.0
      
      !SOIL - IMPORT
      ecp%soil_Phi = 0.0         !Soil porosity (m3/m3)
      ecp%soildepth = 0.0        !Soil depth (m)
      ecp%theta_max = 0.0        !Saturated soil water volume (m/m)
      ecp%k_sat = 0.0            !Saturated hydraulic conductivity
      ecp%root_Phi = 0.0         !Infiltration factor promoted by roots (units?)

      !VEGETATION - EXPORT STATE
      ecp%z0 = 0.0               !Roughness length (m)
      ecp%GCANOPY = 0.0          !Canopy cond. of water vap (mm s-1) = CNC
      ecp%CO2flux = 0.0          !CO2 flux (umol m-2 s-1)
      ecp%GPP = 0.0              !GPP
      !ecp%NPP = 0.0          !NPP
      !ecp%VOCflux = 0.0     !Other kind of fluxes, aerosols from fire, etc.
      !Cell-level diagnostic values - BIOLOGICAL
      !e.g. LAI, biomass pools, nitrogen pools, PFT fractions, GDD, GPP, etc
      ecp%LAI = 0.0
      ecp%betad = 0.0          !Water stress #CALC FROM Soilmoist & SSTAR by PFT
      do n=1,N_DEPTH           !Water stress in layers
        ecp%betadl(n) = 0.0
      end do
      do n=1,N_DEPTH
        ecp%froot(n) = 0.0 !Fraction of roots in soil layer
      end do
      ecp%C_froot = 0.0

      !-----
         
      !VEGETATION - PRIVATE - Initial values not zero.
      ecp%Ci = 0.0127D0         !Internal foliage CO2 (mol/m3) !!Cohort level or patch??
      ecp%Qf = 3.D-6            !Foliage surface vapor mixing ratio (kg/kg)
      
      !METEOROLOGICAL - IMPORT STATE VARIABLES
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      ecp%TcanopyC = 0.0         !Canopy temperatue (Celsius)
      ecp%Qv = 0.0               !Canopy air specif humidity (kg vapor/ kg air)
      ecp%P_mbar = 0.0           !Atmospheric pressure (mb)
      ecp%Ca = 0.0               !@Atmos CO2 conc at surface height (mol/m3).
      ecp%Soilmoist(:) = 0.0 !May be an array by depth (units TBA)
      ecp%fice = 0.0             !Fraction of soil layer that is ice
      ecp%Precip = 0.0           !Precipitation (mm)
      ecp%Ch = 0.0               !Ground to surface heat transfer coefficient 
      ecp%U = 0.0                !Surface layer wind speed (m s-1)

      !Radiation - IMPORT STATE VARIABLE
      !may later be broken down into hyperspectral increments.
      ! in an array
      ecp%Ivis = 0.0              !Incident shortwave 100-2000 nm (W m-2)
      ecp%Idir = 0.0             !Incident PAR 400-700 nm (W m-2)
      ecp%Solarzen = 0.0         !Solar zenith angle
!!! probably shoud remove next line and also remove fdir from structure
!!      ecp%fdir = 0.0             !Fraction of surface vis rad that is direct 

      end subroutine zero_entcell
!**************************************************************************
      subroutine summarize_entcell(ecp)
!@sum Sum patch properties within a grid cell.
      type(entcelltype) :: ecp
      !-----Local variables------------
      type(patch),pointer :: spp, pp
      !real*8 :: albedo(N_BANDS)
      !real*8 :: z0              !Roughness length, average over patch
      !real*8 :: GCANOPY         !Canopy conductance of water vapor (mm/s)
      !real*8 :: CO2flux         !Net CO2 flux up (umol-CO2/m2-gnd/s)
      !real*8 :: Soilmoist(N_DEPTH) !Available soil moisture by depth (mm)
      !real*8 :: N_deposit       !N deposition (kgN/m2)
      !real*8 :: LAI
      !real*8 :: froot(N_DEPTH)
      integer :: ip             !#patches
      integer :: ia             !counter variable

      ecp%LAI = 0.0  !Re-zero
      !call init_patch(ecp%sumpatch,ecp,0.d0) !Summary patch reset to zero area.
      call zero_patch(ecp%sumpatch)
      spp => ecp%sumpatch

      ip = 0
      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
        ip = ip + 1
        spp%age = spp%age + pp%age
        spp%area = spp%area + pp%area 

        ecp%LAI = ecp%LAI + pp%sumcohort%LAI*pp%area

        do ia=1,N_BANDS  !Area-weighted average
          spp%albedo(ia) = spp%albedo(ia) + pp%albedo(ia)*pp%area
        end do
        spp%z0 = spp%z0 + pp%z0*pp%area !Area-weighted average
          
        !* Flux variables for GCM/EWB - patch total
        spp%GCANOPY = spp%GCANOPY + pp%GCANOPY*pp%area
        spp%CO2flux = spp%CO2flux + pp%CO2flux*pp%area

        !* Variables calculated by GCM/EWB - downscaled from grid cell
        do ia=1,N_DEPTH
          spp%Soilmoist(ia) = spp%Soilmoist(ia)+pp%Soilmoist(ia)*pp%area
        end do
        !spp%N_deposit     !N deposition (kgN/m2)

        !* Variables for biophysics and biogeochemistry
        !type(canradtype) :: crad !Data structure for light profile

        !* Disturbance values
        spp%fuel = spp%fuel + spp%fuel*pp%area
        spp%ignition_rate = spp%ignition_rate+spp%ignition_rate*pp%area
        !spp%lambda1(T_SUB)      !Site-averaged fire dist. rate during year
        !spp%disturbance_rate(N_DIST_TYPES)

        !* DIAGNOSTIC SUMMARIES - TO BE COMPLETED LATER
        !* Biomass pools - patch total
        !* SEE avgcohort and sumcohort

        !* Soil pools - patch total
        !spp%REW                 !Relative extractable water (REW)
        !spp%soil_labile_C       !Labile soil carbon (kgC/m2)
        !spp%soil_slow_C         !Slow soil carbon pool (kgC/m2)
        !spp%soil_labile_N       !Labile soil nitrogen (kgN/m2)
        !spp%soil_slow_N         !Slow soil nitrogen (kgN/m2)
        !spp%mineral_N           !Mineralized N (kgN/m2)


        !* Rates of change - patch total
        !spp%dadt                !Rate of change of patch age = 1
        !spp%dpdt          !Rate of change of patch area
        !spp%dwdt                !Rate of change of available soil water

        !* Activity diagnostics - can be summed by month, year, etc.
        spp%GPP = spp%GPP + pp%GPP*pp%area
        !spp%NPP                 !Net primary productivity (kgC/m2/day)
        !spp%Soil_resp           !Soil heterotrophic respiration (kgC/m2/day)

        pp => pp%younger
      end do

      !------DO AVERAGES-----------------------------------------
      !!!CHECK IF SPP%AREA IS ZERO!
      if (ASSOCIATED(ecp%oldest)) then
        spp%age = spp%age/spp%area
        ecp%LAI = ecp%LAI/spp%area
        
        do ia=1,N_BANDS         !Area-weighted average
          spp%albedo(ia) = spp%albedo(ia)/spp%area
        end do
        spp%z0 = spp%z0/spp%area !Area-weighted average
        
        !* Flux variables for GCM/EWB - patch total
        spp%GCANOPY = spp%GCANOPY/spp%area
        spp%CO2flux = spp%CO2flux/spp%area

        !* Variables calculated by GCM/EWB - downscaled from grid cell
        do ia=1,N_DEPTH
          spp%Soilmoist(ia) = spp%Soilmoist(ia)/spp%area
        end do
        !spp%N_deposit     !N deposition (kgN/m2)

        !* Variables for biophysics and biogeochemistry
        !type(canradtype) :: crad !Data structure for light profile

        !* Disturbance values
        spp%fuel = spp%fuel/spp%area
        spp%ignition_rate = spp%ignition_rate/spp%area
      end if

      end subroutine summarize_entcell
!**************************************************************************
!**************************************************************************


      subroutine sum_roots_patches2cell(ecp)
      !@sum Calculate grid-averaged depth-, mass-, and cover-weighted average
      !@sum of fine roots.
      type(entcelltype),pointer :: ecp
      !-----Local variables-------
      type(patch),pointer :: pp
      integer :: n
      real*8 :: froot(N_DEPTH)
      real*8 :: frootC_total
      real*8 :: cf, tcf !cover fraction, total cover fraction

      !* Re-zero summary variable.
      do n=1,N_DEPTH
        froot(n) = 0.0
      end do
      frootC_total = 0.0
      cf = 0.0
      tcf = 0.0

      pp = ecp%oldest
      do while (ASSOCIATED(pp))
        cf = pp%area/ecp%area
        tcf = tcf + cf
        frootC_total = frootC_total + pp%C_froot
        do n=1,N_DEPTH
          froot(n) = froot(n) + cf*pp%sumcohort%froot(n)*pp%C_froot  
        end do
        pp = pp%younger
      end do

      ecp%froot = froot/(tcf*frootC_total)

      end subroutine sum_roots_patches2cell

!**************************************************************************

      subroutine init_simple_entcell( ecp,
     i vegdata,popdens,laidata,hdata,nmdata,frootdata,soildata )
      !@sum Initializes an entcell assuming one cohort per patch.
      use patches, only : summarize_patch
      type(entcelltype) :: ecp
      real*8,intent(in) :: vegdata(N_COVERTYPES) !Veg cover fractions.
      real*8,intent(in) :: popdens(N_COVERTYPES) !Veg population density
      real*8,intent(in) :: laidata(N_COVERTYPES) !LAI
      real*8,intent(in) :: hdata(N_COVERTYPES) !Height
      real*8,intent(in) :: nmdata(N_COVERTYPES) !Nitrogen parameter
      real*8,intent(in) :: frootdata(N_COVERTYPES,N_DEPTH) !Root profile.
      integer,intent(in) :: soildata(N_COVERTYPES)
      !-----Local---------
      integer :: pnum
      type(patch),pointer :: pp, pp_tmp

      ! destroy all existing patches since we are goin to 
      ! re-initialize the cell
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        call patch_destruct(pp_tmp)
      enddo
      nullify(ecp%oldest)
      nullify(ecp%youngest)

      do pnum=1,N_PFT           !One patch with one cohort per pft
      !Get from GISS GCM ## vfraction of grid cell and area.
        if (vegdata(pnum)>0.0) then
          !call insert_patch(ecp,GCMgridareas(j)*vegdata(pnum))
          call insert_patch(ecp,vegdata(pnum),soildata(pnum))
          pp => ecp%youngest
          !## Supply also geometry, clumping index
          ! insert cohort only if population density > 0 (i.e. skip bare soil)
          if ( popdens(pnum) > EPS ) then 
            call insert_cohort(pp,pnum,0.d0,hdata(pnum),
     &           nmdata(pnum),
     &           0.d0,0.d0,0.d0,0.d0,laidata(pnum),0.d0,
     &           frootdata(pnum,:),
     &           0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &           0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &           0.d0,0.d0,0.d0,0.d0,0.d0,
     &           0.d0,0.d0,0.d0,0.d0)
          endif
          call summarize_patch(pp)
        end if
      end do
      call summarize_entcell(ecp)

      call entcell_print(ecp)

      end subroutine init_simple_entcell

 !*********************************************************************

      subroutine entcell_construct(ecp)
      implicit none
      type(entcelltype), pointer :: ecp

      ! allocate memory
      allocate( ecp )
      allocate( ecp%froot(N_DEPTH) )
      allocate( ecp%betadl(N_DEPTH) )
      allocate( ecp%Soilmoist(N_DEPTH) )
      allocate( ecp%Soilmp(N_DEPTH) )
      allocate( ecp%fice(N_DEPTH) )

      ! set pointers
      nullify( ecp%youngest )
      nullify( ecp%oldest   )

      ! for now set all values to zero or defaults
      call zero_entcell(ecp)

      ! maybe the following is not needed, but I guess we wanted
      ! to initialize the cells by default to bare soil (sand)
      call insert_patch(ecp,1.d0,1)

      ! create sumpatch
      call patch_construct(ecp%sumpatch, ecp, 0.d0, 0)

      end subroutine entcell_construct

 !*********************************************************************

      subroutine entcell_destruct(ecp)
      implicit none
      type(entcelltype), pointer :: ecp
      !---
      type(patch), pointer :: pp, pp_tmp

      ! destroy sumpatch
      call patch_destruct(ecp%sumpatch)

      ! destroy other patches
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        call patch_destruct(pp_tmp)
      enddo

      ! deallocate memory
      deallocate( ecp%fice )
      deallocate( ecp%Soilmp )
      deallocate( ecp%Soilmoist )
      deallocate( ecp%betadl )
      deallocate( ecp%froot )
      deallocate( ecp )
      nullify( ecp )

      end subroutine entcell_destruct

 !*********************************************************************
      
      subroutine entcell_print(ecp)
      use patches, only : patch_print
      type(entcelltype) :: ecp
      type(patch), pointer :: pp
      character*1 :: prefix=" "
      character*8 prefix_p
      integer np

      print '(a,"patches:")',prefix

      pp => ecp%oldest
      np = 0
      do while( associated(pp) )
        np = np + 1
        write( prefix_p, '(i2,"      ")' ) np
        call patch_print(pp,prefix//prefix_p)
        pp => pp%younger
      enddo

      print '(a,"sumpatch:")',prefix
      call patch_print(ecp%sumpatch,prefix//"s       ")

      end subroutine entcell_print

 !*********************************************************************
      end module entcells
