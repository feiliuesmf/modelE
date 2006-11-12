      module entcells
!@sum entcells contains routines to perform summaries at the ent grid cell
!@sum level.


      !Ent MODULES TO USE
      use ent_const
      use ent_types
      use ent_pfts
      use patches
      use cohorts

      implicit none
      private
      save

      public zero_entcell, summarize_entcell, entcell_print
      public init_simple_entcell, entcell_construct, entcell_destruct
      public entcell_extract_pfts


      contains
!**************************************************************************

      subroutine zero_entcell(ecp)
      !@sum Zeros import/export variables of entdata type
      type(entcelltype) :: ecp

      ecp%area = 0.0
      
      !Cell-level summary values - PHYSICAL
      !EXPORT - from radiative transfer
!      do n = 1,N_BANDS
!        ecp%albedo(n) = 0.0      !Albedo may be in bands or hyperspectral
!      end do
!      ecp%TRANS_SW = 0.0
      
      !SOIL - IMPORT
      ecp%soil_Phi = 0.0         !Soil porosity (m3/m3)
      ecp%soildepth = 0.0        !Soil depth (m)
      ecp%theta_max = 0.0        !Saturated soil water volume (m/m)
      ecp%k_sat = 0.0            !Saturated hydraulic conductivity
      ecp%root_Phi = 0.0         !Infiltration factor promoted by roots (units?)

      !VEGETATION - EXPORT STATE
      ecp%z0 = 0.0               !Roughness length (m)
!      ecp%GCANOPY = 0.0          !Canopy cond. of water vap (mm s-1) = CNC
!      ecp%CO2flux = 0.0          !CO2 flux (umol m-2 s-1)
!      ecp%GPP = 0.0              !GPP
      !ecp%NPP = 0.0          !NPP
      !ecp%VOCflux = 0.0     !Other kind of fluxes, aerosols from fire, etc.
      !Cell-level diagnostic values - BIOLOGICAL
      !e.g. LAI, biomass pools, nitrogen pools, PFT fractions, GDD, GPP, etc
      ecp%fv = 0.d0
      ecp%LAI= 0.d0
      ecp%heat_capacity = 0.d0
!      ecp%betad = 0.0          !Water stress #CALC FROM Soilmoist & SSTAR by PFT
!      do n=1,N_DEPTH           !Water stress in layers
!        ecp%betadl(n) = 0.0
!      end do
!      do n=1,N_DEPTH
!        ecp%froot(n) = 0.0 !Fraction of roots in soil layer
!      end do
!      ecp%C_froot = 0.0

      !-----
         
      !VEGETATION - PRIVATE - Initial values not zero.
!      ecp%Ci = 0.0127D0         !Internal foliage CO2 (mol/m3) !!Cohort level or patch??
      ecp%Qf = 3.D-6            !Foliage surface vapor mixing ratio (kg/kg)
      
      !METEOROLOGICAL - IMPORT STATE VARIABLES
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      ecp%TcanopyC = 0.0         !Canopy temperatue (Celsius)
!      ecp%Qv = 0.0               !Canopy air specif humidity (kg vapor/ kg air)
      ecp%P_mbar = 0.0           !Atmospheric pressure (mb)
      ecp%Ca = 0.0               !@Atmos CO2 conc at surface height (mol/m3).
      ecp%Soilmoist = 0.0 !Soil moisture avg top 30 cm (volumetric fraction) -PK 6/28/06
      ecp%Soiltemp = 0.0  !Soil temp avg top 30 cm (Celsius)
      ecp%fice = 0.0             !Fraction of soil layer that is ice
      ecp%Ch = 0.0               !Ground to surface heat transfer coefficient 
      ecp%U = 0.0                !Surface layer wind speed (m s-1)

      !Radiation - IMPORT STATE VARIABLE
      !may later be broken down into hyperspectral increments.
      ! in an array
      ecp%IPARdif = 0.0         !Incident diffuse PAR 400-700 nm (W m-2)
      ecp%IPARdir = 0.0         !Incident direct PAR 400-700 nm (W m-2)
      ecp%CosZen = 0.0         !Solar zenith angle
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
      real*8 lai
      real*8 vdata(N_COVERTYPES) ! needed for a hack to compute canopy
                                 ! heat capacity exactly as in GISS GCM

      !call init_patch(ecp%sumpatch,ecp,0.d0) !Summary patch reset to zero area.
      ecp%LAI = 0.d0 !Re-zero
      ecp%fv = 0.d0
      call zero_patch(ecp%sumpatch)
      spp => ecp%sumpatch

      ip = 0
      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
        ip = ip + 1
        spp%age = spp%age + pp%age
        spp%area = spp%area + pp%area

        if ( associated( pp%tallest ) ) ecp%fv = ecp%fv + pp%area
        
        spp%nm = spp%nm + pp%nm*pp%area

        do ia=1,N_COVERTYPES
          spp%LAI(ia) = spp%LAI(ia) + pp%LAI(ia)*pp%area
          ecp%LAI = ecp%LAI + pp%LAI(ia)*pp%area !Area-weighted average
        end do

        do ia=1,N_BANDS  !Area-weighted average
          spp%albedo(ia) = spp%albedo(ia) + pp%albedo(ia)*pp%area
        end do
        spp%z0 = spp%z0 + pp%z0*pp%area !Area-weighted average
        spp%TRANS_SW = spp%TRANS_SW + pp%TRANS_SW*pp%area !Area-weighted average

        !* Flux variables for GCM/EWB - patch total
        spp%GCANOPY = spp%GCANOPY + pp%GCANOPY*pp%area
        spp%CO2flux = spp%CO2flux + pp%CO2flux*pp%area
        spp%Ci = spp%Ci + pp%Ci*pp%area
        spp%betad = spp%betad + pp%betad*pp%area
        spp%soil_resp = spp%soil_resp + pp%soil_resp*pp%area     !added soil resp (umolC/m2/s) -PK 6/14/06
        spp%Tpool = spp%Tpool + pp%Tpool*pp%area                                                
        !* Variables calculated by GCM/EWB - downscaled from grid cell
        spp%Soilmoist = spp%Soilmoist + pp%Soilmoist*pp%area
        do ia=1,N_DEPTH
          spp%betadl(ia) = spp%betadl(ia) + pp%betadl(ia)*pp%area
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
        !spp%Soil_resp           !Soil heterotrophic respiration (kgC/m2/day) !moved to fluxes -PK 6/14/06

        pp => pp%younger
      end do

      !------DO AVERAGES-----------------------------------------
      !!!CHECK IF SPP%AREA IS ZERO!
      if (ASSOCIATED(ecp%oldest)) then
        spp%age = spp%age/spp%area
        spp%nm = spp%nm/spp%area
        do ia=1,N_COVERTYPES
          spp%LAI(ia) = spp%LAI(ia)/spp%area
        end do
        ecp%LAI = ecp%LAI/spp%area

        do ia=1,N_BANDS         !Area-weighted average
          spp%albedo(ia) = spp%albedo(ia)/spp%area
        end do
        spp%z0 = spp%z0/spp%area !Area-weighted average
        spp%TRANS_SW = spp%TRANS_SW/spp%area !Area-weighted average
        
        !* Flux variables for GCM/EWB - patch total
        spp%GCANOPY = spp%GCANOPY/spp%area
        spp%CO2flux = spp%CO2flux/spp%area
        spp%Ci = spp%Ci/spp%area
        spp%betad = spp%betad/spp%area
        spp%soil_resp = spp%soil_resp/spp%area     !added soil resp (umolC/m2/s) -PK 6/14/06
        spp%Tpool = spp%Tpool/spp%area       

        !* Variables calculated by GCM/EWB - downscaled from grid cell
        spp%Soilmoist = spp%Soilmoist/spp%area
        do ia=1,N_DEPTH
          spp%betadl(ia) = spp%betadl(ia)/spp%area
        end do
        !spp%N_deposit     !N deposition (kgN/m2)

        !* Variables for biophysics and biogeochemistry
        !type(canradtype) :: crad !Data structure for light profile

        !* Disturbance values
        spp%fuel = spp%fuel/spp%area
        spp%ignition_rate = spp%ignition_rate/spp%area
      end if

      ! use original formula for canopy heat cpacity
      !lai = ecp%LAI
      !ecp%heat_capacity=(.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow
      !extract vdata exactly like in GISS GCM
      !!vdata(:) = 0.d0
      !!call entcell_extract_pfts(ecp, vdata(2:) )
      !!ecp%heat_capacity=GISS_calc_shc(vdata)

      end subroutine summarize_entcell
!**************************************************************************
!**************************************************************************

#ifdef SUMROOTSCELL
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
#endif

!**************************************************************************

      subroutine init_simple_entcell( ecp,
     i     vegdata,popdens,laidata,hdata,dbhdata,craddata,
     i     cpooldata,nmdata,
     i     frootdata,soildata,albedodata,soil_texture)
      !@sum Initializes an entcell assuming one cohort per patch.
      use patches, only : summarize_patch
      type(entcelltype) :: ecp
      real*8,intent(in) :: vegdata(N_COVERTYPES) !Veg cover fractions.
      real*8,intent(in) :: popdens(N_COVERTYPES) !Veg population density
      real*8,intent(in) :: laidata(N_COVERTYPES) !LAI
      real*8,intent(in) :: hdata(N_COVERTYPES) !Height
      real*8,intent(in) :: dbhdata(N_COVERTYPES) !Woody plant diameter (cm)
      real*8,intent(in) :: craddata(N_COVERTYPES)
      real*8,intent(in) :: cpooldata(N_COVERTYPES,N_BPOOLS)
      real*8,intent(in) :: nmdata(N_COVERTYPES) !Nitrogen parameter
      real*8,intent(in) :: frootdata(N_COVERTYPES,N_DEPTH) !Root profile.
      integer,intent(in) :: soildata(N_COVERTYPES)
      real*8,intent(in) :: albedodata(N_BANDS,N_COVERTYPES) !patch, NOTE:snow
      real*8,intent(in) :: soil_texture(N_SOIL_TYPES) !Veg cover fractions.
      !-----Local---------
      integer :: ncov, pft
      type(patch),pointer :: pp, pp_tmp

      print *,"entered init_simple_entcell"
!      if ( .not. associated(ecp) ) 
!     &      call stop_model("init_simple_entcell 1",255)
      call entcell_print(6,ecp)

      ! destroy all existing patches since we are going to 
      ! re-initialize the cell
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        call patch_destruct(pp_tmp)
      enddo
      nullify(ecp%oldest)
      nullify(ecp%youngest)

!      do ncov=1,N_PFT           !One patch with one cohort per pft or bare
      do ncov=1,N_COVERTYPES           !One patch with one cohort per pft
        pft = ncov - COVEROFFSET
        !### Get from GISS GCM ## vfraction of grid cell and area.

        if (vegdata(ncov)>0.0) then
         !call insert_patch(ecp,GCMgridareas(j)*vegdata(pnum))
          call insert_patch(ecp,vegdata(ncov),soildata(ncov))
          pp => ecp%youngest
          !## Supply also geometry, clumping index
          ! insert cohort only if population density > 0 (i.e. skip bare soil)
          if ( popdens(ncov) > EPS ) then 
            if ( pft < 1 .or. pft > N_PFT ) then
              print *,"init_simple_entcell: wrong pft:", pft
              call patch_print(6,pp,"ERROR ")
              call stop_model("init_simple_entcell: wrong pft",255)
            endif
            call insert_cohort(pp,pft,popdens(ncov),hdata(ncov),
     &           nmdata(ncov),
     &           craddata(ncov),0.d0,dbhdata(ncov),0.d0,laidata(ncov),
     &           0.d0, frootdata(ncov,:),
     &           0.d0,cpooldata(ncov,FOL),0.d0,cpooldata(ncov,SW),0.d0,
     &           cpooldata(ncov,HW),0.d0,
     &           cpooldata(ncov,LABILE),0.d0,
     &           cpooldata(ncov,FR),0.d0,0.d0,0.d0,
     &           0.d0,0.d0,0.d0,0.d0,0.d0,
     &           0.d0,0.d0)
          endif
          call summarize_patch(pp)

          !CALL CALC_ALBEDO HERE
          pp%albedo = albedodata(:,ncov) !##GISS HACK
        end if
      end do

      ! soil textures for CASA
      ecp%soil_texture(:) = soil_texture(:)

      call summarize_entcell(ecp)

      print *,"leaving init_simple_entcell:"
      call entcell_print(6,ecp)

      end subroutine init_simple_entcell

 !*********************************************************************

      subroutine entcell_construct(ecp)
      implicit none
      type(entcelltype), pointer :: ecp

      ! allocate memory
      allocate( ecp )
      allocate( ecp%Soilmp(N_DEPTH) )
      allocate( ecp%fice(N_DEPTH) )

      ! set pointers
      nullify( ecp%youngest )
      nullify( ecp%oldest   )
      nullify( ecp%sumpatch )

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
      deallocate( ecp )
      nullify( ecp )

      end subroutine entcell_destruct

 !*********************************************************************
      
      subroutine entcell_print(iu, ecp)
      use patches, only : patch_print
      integer, intent(in) :: iu
      type(entcelltype), intent(in) :: ecp
      !---
      type(patch), pointer :: pp
      character*1 :: prefix=" "
      character*8 prefix_p
      integer np

      write(iu, '(a,"entcell:")') prefix
      !print '(a," = ",f10.7)',"GCANOPY ",ecp%GCANOPY

      write(iu, '(a,"patches:")') prefix
      pp => ecp%oldest
      np = 0
      do while( associated(pp) )
        np = np + 1
        write( prefix_p, '(i2,"      ")' ) np
        call patch_print(iu, pp, prefix//prefix_p)
        pp => pp%younger
      enddo

      write(iu, '(a,"sumpatch:")') prefix
      call patch_print(iu, ecp%sumpatch, prefix//"s       ")

      end subroutine entcell_print

 !*********************************************************************

      subroutine entcell_extract_pfts(ecp, vdata)
      type(entcelltype) :: ecp
      real*8 :: vdata(:)
      !---
      type(patch), pointer :: pp
      real*8 :: vdata_patch(size(vdata))

      vdata(:) = 0.d0
      pp => ecp%oldest
      do while( associated(pp) )
        call patch_extract_pfts(pp, vdata_patch)
        vdata(:) = vdata(:) + vdata_patch(:)*pp%area
        pp => pp%younger
      enddo

      end subroutine entcell_extract_pfts

 !*********************************************************************
      end module entcells
