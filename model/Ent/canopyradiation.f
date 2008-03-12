      module canopyrad
!@sum Routines for calculating canopy radiation and albedo.
!@auth W.Ni-Meister

      !Ent MODULES TO USE
      use ent_const
      use ent_types

      implicit none
      private
      save

      public recalc_radpar_cell, get_patchalbedo, recalc_radpar
      public get_canopy_rad, GORT_clumping, TwoStream       
      public :: gort_input
      public :: pft_pstate_type

      ! type definition
      type gort_input
        real*8 :: dens_tree         ! Crown count density (m-2)
        real*8 :: dens_foliage      ! foliage area volume density of a single crown(m-1)
        real*8 :: h1, h2            ! lower and upper bound of crown centers (m)
        real*8 :: delta_z           ! height increase for each height level
        real*8 :: horz_radius       ! horizontal crown radius (m)
        real*8 :: vert_radius       ! vertical crown radius (m)
        real*8 :: zenith            ! solar zenith angle (rad)
        real*8 :: LAI
        integer :: pft
      end type gort_input

      type profile_params
        real*8, dimension(:), pointer :: height_levels ! height level for each profile
        real*8, dimension(:), pointer :: fp            ! foliage profile, no clump involved
        real*8, dimension(:), pointer :: rdfp           ! clumping factor * foliage profile
        real*8, dimension(:), pointer :: rifp           ! clumping factor * foliage profile
        real*8, dimension(:), pointer :: efp           ! effective foliage profile for each profile
      end type profile_params

      !----------------------------------------------------
      ! pft physical state variables structure
      !----------------------------------------------------
      type pft_pstate_type
        real*8, pointer :: albd(:,:)    !surface albedo (direct)                       (numrad)
        real*8, pointer :: albi(:,:)    !surface albedo (indirect)                     (numrad)
        real*8, pointer :: fabd(:,:)    !flux absorbed by veg per unit direct flux     (numrad)
        real*8, pointer :: fabi(:,:)    !flux absorbed by veg per unit diffuse flux    (numrad)
        real*8, pointer :: ftdd(:,:)    !down direct flux below veg per unit dir flx   (numrad)
        real*8, pointer :: ftid(:,:)    !down diffuse flux below veg per unit dir flx  (numrad)
        real*8, pointer :: ftii(:,:)    !down diffuse flux below veg per unit dif flx  (numrad)
        real*8, pointer :: frid(:,:)    !up diffuse flux below veg per unit dir flx    (numrad)
        real*8, pointer :: frii(:,:)    !up diffuse flux below veg per unit dif flx    (numrad) 
        real*8, pointer :: gdir(:)      !leaf projection in solar direction (0 to 1)
        real*8, pointer :: omega(:,:)   !fraction of intercepted radiation that is scattered (0 to 1)
        real*8, pointer :: isha(:,:)    !shaded flux per unit incident flux            (numrad)         
        real*8, pointer :: isun(:,:)    !sunlit flux per unit incident flux            (numrad)         
      end type pft_pstate_type

      real*8 :: K
      !real*8, parameter :: PI=3.14159265
      real*8, parameter :: ang_dif = 1.0071   ! 57.7 deg in rad
      integer, parameter :: numrad      =   2   ! number of solar radiation bands: vis, nir
      integer, parameter :: lbp = 1
      integer, parameter :: ubp = 1
      integer, parameter :: num_vegsol = ubp - lbp + 1
      integer, parameter :: num_pft = 16
 
      contains
      !*********************************************************************
      
      subroutine get_patchalbedo(pp)
!@sum !* Return albedo for patch.
      !* This version reads in vegetation structure from GISS data set.
      !use ent_GISSveg, only : GISS_veg_albedo
      implicit none
!!!       integer :: jday  !Day of year. - not needed for real simulation
      !type(timestruct) :: tt
      type(patch),pointer :: pp
      !------

      !---------------------------------------------------------------
      !* GISS version for calculating seasonal aalbveg.
      !* Assumes each GISS grid cell vegetation fraction corresponds to
      !* one patch, which contains one cohort, which gives the albedo.
      !  Should be same as assigning
      !  pp%albedo = aalbveg(i,j)
      !---------------------------------------------------------------

      !* Prescribed albedo is passed in.
 !     call GISS_veg_albedo(pp%cellptr%latj,pp%tallest%pft,
 !    &     jday, pp%albedo)


      !* Prognostic albedo is calculated here.
      !---------------------------------------------------------------
      !* Ent template for GORT clumping index canopy radiative transfer.
      !---------------------------------------------------------------
!      do cop=pp%tallest
        !Get albedo of cop%pft for given tt
        !   - Summarize foliage density in layers
        !* CALCULATE PATCH ALBEDO USING GORT CLUMPING INDEX
!        cop = cop%shorter
!        if (.NOT.ASSOCIATED(cop)) then exit
!      end do

      end subroutine get_patchalbedo


      !*********************************************************************


      subroutine recalc_radpar_cell(pptr)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.  At patch level.

      type(entcelltype) :: pptr

!!! passing actual structure instead of pointer, change if necessary
!!!      if (ASSOCIATED(pptr)) then 
        !call get_patchalbedo(tt,pptr) !*This is called by summarize_patch
        !call GORT_clumping(pptr)
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
!!!      end if

      end subroutine recalc_radpar_cell


      subroutine recalc_radpar(pptr)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.  At patch level.

      type(patch) :: pptr

!!! passing actual structure instead of pointer, change if necessary
!!!      if (ASSOCIATED(pptr)) then 
        !call get_patchalbedo(tt,pptr) !*This is called by summarize_patch
        !call GORT_clumping(pptr)
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
!!!      end if

      end subroutine recalc_radpar


      !*********************************************************************
      subroutine recalc_radpar_OLD(entcell)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.  At entcell level.

      type(entcelltype) :: entcell
      type(patch),pointer :: pptr

      pptr = entcell%youngest

      do while (ASSOCIATED(pptr)) 
        !call GORT_clumping(pptr)
        pptr = pptr%older
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      end do

      end subroutine recalc_radpar_OLD

      !*********************************************************************

      subroutine calc_canopy_rad(pptr, h)
!@sum Get incident light profiles in canopy and return in pptr%crad
      type(patch),pointer :: pptr
      real*8 :: h               !Height in canopy
      !--------------------------


      !Get incident light profiles in canopy and return in crad
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine calc_canopy_rad
      
      !*********************************************************************
      subroutine get_canopy_rad(pptr)

      type(patch),pointer :: pptr
      ! real*8 :: h               !Height in canopy
      
      real*8 :: Id, Ii
      real*8, dimension(:), pointer :: ffp, rdfp, rifp, sunlit, shaded
      real*8, dimension(:), pointer :: height_levels
      integer :: i, n_height_level
      
      integer :: filter_vegsol(num_vegsol), ivt(num_vegsol) 
      integer :: pcolumn(num_vegsol)
      real*8 :: vai(num_vegsol) 
      real*8 :: rho(num_vegsol,numrad), tau(num_vegsol,numrad)
      real*8 :: coszen(num_vegsol), esai(num_vegsol), elai(num_vegsol)
      real*8 :: rhol(num_pft,numrad), taul(num_pft,numrad) 
      real*8 :: rhos(num_pft,numrad), taus(num_pft,numrad)

      type(pft_pstate_type) :: sout
    
      integer  :: fp,g,c,p                ! indices
      integer  :: ib                      ! band index
      integer  :: len
      real*8, parameter :: mpe = 1.e-06     ! prevents overflow for division by zero
      real*8 :: wl(lbp:ubp)                 ! fraction of LAI+SAI that is LAI
      real*8 :: ws(lbp:ubp)                 ! fraction of LAI+SAI that is SAI
      real*8 :: itot(2)                     ! red and nir
      real*8 :: reflect(2), tran(2), summ(2)
      real*8 :: dz
    
      data (ivt(i),i=1,num_vegsol) /1/
      data (pcolumn(i),i=1,num_vegsol) /num_vegsol*1/    ! 1-soil 2-snow
      data (elai(i),i=1,num_vegsol) /num_vegsol*1.2/
      data (esai(i),i=1,num_vegsol) /num_vegsol*1/
      data (rhol(i,1),i=1,16) /0.07, 0.07, 0.07, 0.10, 0.10, 0.10, 
     &   0.10, 0.10, 0.07, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11/
      data (rhol(i,2),i=1,16) /0.35, 0.35, 0.35, 0.45, 0.45, 0.45,
     &   0.45, 0.45, 0.35, 0.45, 0.45, 0.58, 0.58, 0.58, 0.58, 0.58/
      data (rhos(i,1),i=1,16) /0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 
     &    0.16, 0.16, 0.16, 0.16, 0.16, 0.36, 0.36, 0.36, 0.36, 0.36/
      data (rhos(i,2),i=1,16) /0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 
     &    0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58, 0.58, 0.58, 0.58/
      data (taul(i,1),i=1,16) /0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
     &   0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.07, 0.07, 0.07, 0.07/
      data (taul(i,2),i=1,16) /0.10, 0.10, 0.10, 0.25, 0.25, 0.25,  
     &   0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25/
      data (taus(i,1),i=1,16) /0.001, 0.001, 0.001, 0.001, 0.001, 
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.220, 0.220, 
     &   0.220, 0.220, 0.220/
      data (taus(i,2),i=1,16) /0.001, 0.001, 0.001, 0.001, 0.001, 
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380, 
     &   0.380, 0.380, 0.380/

      Id = 1.0
      Ii = 1 - Id
      call GORT_clumping(pptr, height_levels, ffp, rdfp, rifp, sunlit, 
     &    shaded)
    
      ! Weight reflectance/transmittance by lai and sai
      ! Only perform on vegetated pfts where coszen > 0

      do i=1, num_vegsol
         filter_vegsol(i) = i
         coszen(i) = pptr%cellptr%CosZen
      end do

      do fp = 1,num_vegsol
         p = filter_vegsol(fp)
         vai(p) = elai(p) + esai(p)
         wl(p) = elai(p) / max( vai(p), mpe )
         ws(p) = esai(p) / max( vai(p), mpe )
      end do

      do ib = 1, numrad
         do fp = 1,num_vegsol
            p = filter_vegsol(fp)
            rho(p,ib) = max( rhol(ivt(p),ib)*wl(p) + 
     &        rhos(ivt(p),ib)*ws(p), mpe )
            tau(p,ib) = max( taul(ivt(p),ib)*wl(p) + 
     &        taus(ivt(p),ib)*ws(p), mpe )
         end do
      end do

      call TwoStream(lbp, ubp, filter_vegsol, num_vegsol, pcolumn, 
     &    coszen, Id, height_levels, rdfp, rifp, rho, tau, sout)

      pptr%albedo(1) = Id * sout%albd(1,1) + Ii * sout%albi(1,1)
      pptr%albedo(2) = Id * sout%albd(1,2) + Ii * sout%albi(1,2)
 
      end subroutine get_canopy_rad

      !*********************************************************************
      subroutine GORT_clumping(pptr, height_levels, ffp, rdfp, rifp, 
     &   sunlit, shaded)
!@sum Calculate the GORT clumping index in canopy layers and save into
!@sum variable ppt%crad
      type(patch), pointer :: pptr
      type(gort_input) ::  gin1  ! input structure for gort model
      type(gort_input), dimension(:), pointer :: gin
      integer :: num_profiles
      real*8, dimension(:), pointer :: height_levels, ffp, rdfp, rifp 
      real*8, dimension(:), pointer :: sunlit, shaded
      
      call geo_to_gin(pptr, gin)
      num_profiles = size(gin)
      if (num_profiles .eq. 1.d0) then
          gin1%dens_tree = gin(1)%dens_tree
          gin1%h1 = gin(1)%h1
          gin1%h2 = gin(1)%h2
          gin1%horz_radius = gin(1)%horz_radius
          gin1%vert_radius = gin(1)%vert_radius
          gin1%zenith = gin(1)%zenith
          gin1%delta_z = gin(1)%delta_z
          gin1%dens_foliage = gin(1)%dens_foliage
          gin1%pft = gin(1)%pft
          gin1%LAI = gin(1)%LAI
          call run_single_gort(gin1, pptr, height_levels, ffp, rdfp, 
     &        rifp, sunlit, shaded)
      else
          call run_convolute_gort(gin, pptr, height_levels, ffp, rdfp, 
     &        rifp, sunlit, shaded)
      endif
      
      end subroutine GORT_clumping

      !*********************************************************************

      subroutine run_single_gort(gin, pptr, height_levels, fp, rdfp, 
     &    rifp, sunlit, shaded)
!      subroutine run_single_gort(gin, direct_light_ratio, ivt, height_levels,        &
!            fp, rdfp, rifp, transmit, sunlit, shaded)

      ! Input
      type(gort_input) ::  gin  ! input structure for gort model
      type(patch), pointer :: pptr

      integer :: N_height_level
      real*8 :: clumpd, clumpi
      real*8, dimension(:), pointer :: fp, rdfp, rifp, efp 
      real*8, dimension(:), pointer :: height_levels 
      real*8, dimension(:), pointer :: sunlit     ! sunlit leaf area fraction
      real*8, dimension(:), pointer :: shaded     ! shaded leaf area fraction
      real*8 :: tmp_zenith

      ! ----------------------------------------------------------
      ! DEBUG
      ! print *, 'I am in run_single_gort ..........'
   
     
      K = get_K(cos(gin%zenith),gin%pft)

      ! Check the validity of the inputs
      call check_inputs(gin)

      ! define the vertical layers
      call get_height_level(gin, height_levels)
      N_height_level = size(height_levels)

      ! get intermediate variables, clump and fp
      allocate(fp(N_height_level))
      allocate(rdfp(N_height_level))
      allocate(rifp(N_height_level))
      allocate(efp(N_height_level))
      call get_foliage_profile(gin,height_levels,fp)
      clumpd = get_analytical_clump(gin)
      tmp_zenith=gin%zenith
      gin%zenith = ang_dif
      clumpi = get_analytical_clump(gin)
      gin%zenith = tmp_zenith

      ! add clump factor into foliage profile data
      efp = K * clumpd * fp
      rdfp = clumpd * fp
      rifp = clumpi * fp
      
      allocate(sunlit(N_height_level))
      allocate(shaded(N_height_level))
      sunlit = T(efp, gin%delta_z, gin%zenith)
      shaded = 1 - sunlit
            
      end subroutine run_single_gort

      !*********************************************************************

      subroutine run_convolute_gort (gin, pptr, height_levels, fp, rdfp, 
     &    rifp, sunlit, shaded)            ! convoluted profile values
!      subroutine run_convolute_gort (gin, direct_light_ratio, ivt,     &   ! All inputs for gort
!          height_levels, fp, rdfp, rifp, transmit, sunlit, shaded)            ! convoluted profile values

      ! Input
      type(gort_input), dimension(:), pointer :: gin
      type(patch), pointer :: pptr


      ! Output
      real*8, dimension(:), pointer   :: rdfp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: rifp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: sunlit     ! sunlit leaf area fraction
      real*8, dimension(:), pointer   :: shaded     ! shaded leaf area fraction
      real*8, dimension(:), pointer   :: height_levels  ! Convolute profile height levels



      ! local variables
      integer :: num_profiles, iprofile, N_height_level, i
      integer :: N_convolute_height_level, ilevel, iszn
      real*8    :: szn, delta_z, dszn
      real*8    :: tmp_zenith
      real*8    :: clumpd, clumpi
      type(profile_params), dimension(:), pointer :: all_convolute_input
      real*8, dimension(:), pointer   :: fp, efp   ! convoluted effective foliage profile
      real*8, dimension(:), pointer   :: tmp_fp, tmp_rdfp, tmp_rifp  
      real*8, dimension(:), pointer   :: tmp_height_levels, tmp_efp

      real*8, dimension(:), pointer   :: d_fp, d_rdfp, d_rifp, d_efp 
      real*8, dimension(:), pointer   :: d_height_levels
      real*8, dimension(:,:), pointer :: convoluted_szn_transmit
      !-----------------------------------------------------------------

      ! DEBUG
      print *, 'I am in run_convolute_gort.........'
      K = get_K2(gin)
      num_profiles = size(gin)


      ! check to verify that all profiles in all_gort_in have the same zenith angle
      szn = gin(1)%zenith
      do iprofile=2, num_profiles
        if (gin(iprofile)%zenith .ne. szn) then
          print *, 'The sun zenith angle of all the participating 
     &	     profiles should be the same. Exit.....'
          exit
        end if
      end do

      !print *, 'The sun zenith angle for the convoluted profile is ', szn


      ! Get effective efp for each profile

      nullify(all_convolute_input)
      allocate(all_convolute_input(num_profiles))

      do iprofile=1, num_profiles

        call check_inputs(gin(iprofile))

        ! define the vertical layers
        call get_height_level(gin(iprofile), tmp_height_levels) !tmp_height_levels is allocated in the subroutine
        N_height_level = size(tmp_height_levels)

        ! get intermediate variables, clump and fp
        allocate(tmp_fp(N_height_level))
        allocate(tmp_rdfp(N_height_level))
        allocate(tmp_rifp(N_height_level))
        allocate(tmp_efp(N_height_level))
   
        ! consider putting the following two subroutines as function instead
        call get_foliage_profile(gin(iprofile), tmp_height_levels, 
     &      tmp_fp) 
        clumpd = get_analytical_clump(gin(iprofile))
        tmp_zenith=gin(iprofile)%zenith
        gin(iprofile)%zenith = ang_dif
        clumpi = get_analytical_clump(gin(iprofile))
        gin(iprofile)%zenith = tmp_zenith

        ! add clump factor into foliage profile data 
        tmp_efp = K * clumpd * tmp_fp 
        tmp_rdfp = clumpd * tmp_fp
        tmp_rifp = clumpi * tmp_fp

        !DEBUG
        !print *, tmp_fp(10), clumpd, tmp_efp(10)
        !print *, minval(tmp_height_levels), maxval(tmp_height_levels)
        !print *, 'Assinging single_gort output to convolute_input'

        ! fill in convolute_input
        all_convolute_input(iprofile)%fp => tmp_fp
        all_convolute_input(iprofile)%rdfp => tmp_rdfp
        all_convolute_input(iprofile)%rifp => tmp_rifp
        all_convolute_input(iprofile)%efp => tmp_efp
        all_convolute_input(iprofile)%height_levels => tmp_height_levels
     
        !DEBUG
        !print *, all_convolute_input(iprofile)%efp(10)

        ! deallocate and null the tmp pointers
        nullify(tmp_fp)
        nullify(tmp_rdfp)
        nullify(tmp_rifp)
        nullify(tmp_efp)
        nullify(tmp_height_levels)

      end do


      ! Convolute fp, efp by calling convolute, but here only fp is final output for run_gort
      ! efp will be used for getting transmittance for direct light

      ! print *, 'Calling convolute ......' 
      call convolute(all_convolute_input,       ! input
     &   fp, rdfp, rifp, efp, height_levels)    ! convoluted output, also run_gort output


      ! Compute the transmit, sunlit, shaded from efp
      ! Sunlit and Shaded will be the run_gort output, but transmit will be depending on 
      ! if transmit from diffuse light needs to be considered
      N_convolute_height_level = size(height_levels)
      delta_z = height_levels(2)-height_levels(1)
      allocate(sunlit(N_convolute_height_level))
      allocate(shaded(N_convolute_height_level))
      sunlit = T(efp, delta_z, szn)
      shaded = 1 - sunlit

      do iprofile=1, num_profiles
         deallocate(all_convolute_input(iprofile)%fp)
         deallocate(all_convolute_input(iprofile)%rdfp)
         deallocate(all_convolute_input(iprofile)%rifp)
         deallocate(all_convolute_input(iprofile)%efp)
         deallocate(all_convolute_input(iprofile)%height_levels)
      end do
      deallocate(all_convolute_input)

      end subroutine run_convolute_gort

      !*********************************************************************

      subroutine TwoStream(lbp, ubp, filter_vegsol, num_vegsol, pcolumn, 
     &    coszen, Id, height_levels, rdfp, rifp, rho, tau, sout)
!      subroutine TwoStream (lbp, ubp, filter_vegsol, num_vegsol, ivt, pcolumn, &
!                         coszen, Id, height_levels, rdfp, rifp, rho, tau, sout)
!
! !DESCRIPTION:
! Two-stream fluxes for canopy radiative transfer
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
!
! !USES:
!    use clmtype
!    use clm_varpar, only : numrad
!    use clm_varcon, only : omegas, tfrz, betads, betais
!
! !ARGUMENTS:
!    implicit none
      type(patch),pointer :: pptr
      integer , intent(in)  :: lbp, ubp                 ! pft bounds
      integer , intent(in)  :: filter_vegsol(ubp-lbp+1) ! filter for vegetated pfts with coszen>0
      integer , intent(in)  :: num_vegsol               ! number of vegetated pfts where coszen>0

!      integer , intent(in)  :: ivt(lbp:ubp)             ! pft vegetation type
      integer :: pcolumn(num_vegsol)
      ! integer , intent(in)  :: pcolumn(lbp:ubp)         ! column of corresponding pft

      real*8 :: coszen(num_vegsol)
      ! real*8, intent(in)  :: coszen(lbp:ubp)        ! cosine solar zenith angle for next time step
      real*8, intent(in)  :: Id                       ! direct light ratio, used for isha and isun
      real*8, intent(in), dimension(:)  :: height_levels
      real*8, intent(in), dimension(:)  :: rdfp, rifp
      real*8, intent(in)  :: rho(lbp:ubp,numrad)    ! leaf/stem refl weighted by fraction LAI and SAI
      real*8, intent(in)  :: tau(lbp:ubp,numrad)    ! leaf/stem tran weighted by fraction LAI and SAI

      type(pft_pstate_type) :: sout    ! output structure for 2stream model

! !Constant from CLM code to avoid using 'use' model    
      real*8, parameter :: betads  = 0.5            ! two-stream parameter betad for snow
      real*8, parameter :: betais  = 0.5            ! two-stream parameter betai for snow
      real*8 :: omegas(numrad)           ! two-stream parameter omega for snow by band
      data (omegas(i),i=1,numrad) /0.8, 0.4/
      real*8, parameter :: SHR_CONST_TKFRZ   = 273.15       ! freezing T of fresh water          ~ K 
      real*8, parameter :: tfrz   = SHR_CONST_TKFRZ !freezing temperature [K]

!      integer, parameter :: numpft      =   16   ! number of plant function types
!      real :: xl(numpft)         ! ecophys const - leaf/stem orientation index
!      data (xl(i),i=1,numpft) /0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, &
!                 0.25, 0.01, 0.25, 0.25, -0.30, -0.30, -0.30, -0.30, -0.30/
      real*8 :: albgrd(2,numrad)   ! ground albedo (direct) (column-level)
      real*8 :: albgri(2,numrad)   ! ground albedo (diffuse)(column-level)
      data (albgrd(1,i),i=1,numrad) /0.075, 0.314/
      data (albgri(1,i),i=1,numrad) /0.075, 0.314/
      data (albgrd(2,i),i=1,numrad) /0.939, 0.787/
      data (albgri(2,i),i=1,numrad) /0.939, 0.787/
      real*8 :: Ii                 !diffuse light ratio
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! Modified for speedup: Mariana Vertenstein, 8/26/02
! Vectorized routine: Mariana Vertenstein:  8/20/03
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
!    integer , pointer :: pcolumn(:)    ! column of corresponding pft
!    real, pointer :: albgrd(:,:)   ! ground albedo (direct) (column-level)
!    real, pointer :: albgri(:,:)   ! ground albedo (diffuse)(column-level)
!    real, pointer :: t_veg(:)      ! vegetation temperature (Kelvin)
!    real, pointer :: fwet(:)       ! fraction of canopy that is wet (0 to 1)
!    integer , pointer :: ivt(:)        ! pft vegetation type
!    real, pointer :: xl(:)         ! ecophys const - leaf/stem orientation index
!
! local pointers to implicit out scalars
!
      real*8, pointer :: albd(:,:)     ! surface albedo (direct)
      real*8, pointer :: albi(:,:)     ! surface albedo (diffuse)
      real*8, pointer :: fabd(:,:)     ! flux absorbed by veg per unit direct flux
      real*8, pointer :: fabi(:,:)     ! flux absorbed by veg per unit diffuse flux
      real*8, pointer :: ftdd(:,:)     ! down direct flux below veg per unit dir flx
      real*8, pointer :: ftid(:,:)     ! down diffuse flux below veg per unit dir flx
      real*8, pointer :: ftii(:,:)     ! down diffuse flux below veg per unit dif flx
      real*8, pointer :: frid(:,:)     ! up diffuse flux below veg per unit dir flx
      real*8, pointer :: frii(:,:)     ! up diffuse flux below veg per unit dif flx
      real*8, pointer :: gdir(:)       ! leaf projection in solar direction (0 to 1)
      real*8, pointer :: omega(:,:)    ! fraction of intercepted radiation that is scattered (0 to 1)
      real*8, pointer :: isha(:,:)     ! shaded flux per unit incident flux                     
      real*8, pointer :: isun(:,:)     ! sunlit flux per unit incident flux                     
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
      integer  :: fp,p,c         ! array indices
      !integer  :: ic               ! 0=unit incoming direct; 1=unit incoming diffuse
      integer  :: ib, i          ! waveband number
      real*8 :: cosz             ! 0.001 <= coszen <= 1.000
      real*8 :: asu              ! single scattering albedo
      real*8 :: chil(lbp:ubp)    ! -0.4 <= xl <= 0.6
      real*8 :: twostext(lbp:ubp)! optical depth of direct beam per unit leaf area
      real*8 :: avmu(lbp:ubp)    ! average diffuse optical depth
      real*8 :: omegal           ! omega for leaves
      real*8 :: betai            ! upscatter parameter for diffuse radiation
      real*8 :: betail           ! betai for leaves
      real*8 :: betad            ! upscatter parameter for direct beam radiation
      real*8 :: betadl           ! betad for leaves
      real*8 :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 ! temporary
      real*8 :: p1,p2,p3,p4,s1,s2,u1,u2,u3                        ! temporary
      real*8 :: b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10   ! temporary
      real*8 :: phi1,phi2,sigma                                   ! temporary
      real*8 :: temp0(lbp:ubp),temp1,temp2(lbp:ubp)               ! temporary
      real*8 :: t1
      real*8, pointer :: evai(:)
      real*8 :: delta_z, vai(lbp:ubp)

      integer :: N_height_level, ilevel
!      integer , intent(in)  :: ivt(lbp:ubp)             ! pft vegetation type

!-----------------------------------------------------------------------
      N_height_level = size(height_levels)
      Ii = 1 - Id

!    if (ubp-lbp > 0) then 
!       allocate(vai(lbp:ubp))
!    end if
      allocate(evai(N_height_level))
      allocate(sout%albd(lbp:ubp,numrad))
      allocate(sout%albi(lbp:ubp,numrad))
      allocate(sout%fabd(N_height_level,numrad))
      allocate(sout%fabi(N_height_level,numrad))
      allocate(sout%ftdd(N_height_level,numrad))
      allocate(sout%ftid(N_height_level,numrad))
      allocate(sout%ftii(N_height_level,numrad))
      allocate(sout%frid(N_height_level,numrad))
      allocate(sout%frii(N_height_level,numrad))
      allocate(sout%gdir(lbp:ubp))
      allocate(sout%omega(lbp:ubp,numrad))
      allocate(sout%isha(N_height_level,numrad))
      allocate(sout%isun(N_height_level,numrad))

    ! Assign local pointers to derived subtypes components (column-level)

!    albgrd  => clm3%g%l%c%cps%albgrd
!    albgri  => clm3%g%l%c%cps%albgri

    ! Assign local pointers to derived subtypes components (pft-level)

!      pcolumn => clm3%g%l%c%p%column
!      fwet    => sout%fwet
!      t_veg   => clm3%g%l%c%p%pes%t_veg
!      ivt     => clm3%g%l%c%p%itype
      albd    => sout%albd
      albd    => sout%albd
      albi    => sout%albi
      fabd    => sout%fabd
      fabi    => sout%fabi
      ftdd    => sout%ftdd
      ftid    => sout%ftid
      ftii    => sout%ftii
      frid    => sout%frid
      frii    => sout%frii
      gdir    => sout%gdir
      omega   => sout%omega
      isha    => sout%isha
      isun    => sout%isun
!      xl      => pftcon%xl

      ! Calculate two-stream parameters omega, betad, betai, avmu, gdir, twostext.
      ! Omega, betad, betai are adjusted for snow. Values for omega*betad
      ! and omega*betai are calculated and then divided by the new omega
      ! because the product omega*betai, omega*betad is used in solution.
      ! Also, the transmittances and reflectances (tau, rho) are linear
      ! weights of leaf and stem values.

      do fp = 1,num_vegsol
         p = filter_vegsol(fp)
       
         ! note that the following limit only acts on cosz values > 0 and less than 
         ! 0.001, not on values cosz = 0, since these zero have already been filtered
         ! out in filter_vegsol
         ! cosz = max(0.001, coszen(p))
       
         ! chil(p) = min( max(xl(ivt(p)), -0.4), 0.6 )
         ! if (abs(chil(p)) <= 0.01) chil(p) = 0.01
         ! phi1 = 0.5 - 0.633*chil(p) - 0.330*chil(p)*chil(p)
         ! phi2 = 0.877 * (1.-2.*phi1)
         ! gdir(p) = phi1 + phi2*cosz
         ! twostext(p) = gdir(p)/cosz
         gdir(p) = K
         twostext(p) = gdir(p)/cosz
	 
         ! print *, 'coszen=', coszen(p)
         avmu(p) = (1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
         temp0(p) = gdir(p) + phi2*cosz
         temp1 = phi1*cosz
         temp2(p) = (1. - temp1/temp0(p) * log((temp1+temp0(p))/temp1))
      end do

      do ib = 1, numrad
         do fp = 1,num_vegsol
            p = filter_vegsol(fp)
            c = pcolumn(p)

            omegal = rho(p,ib) + tau(p,ib)
            asu = 0.5*omegal*K/temp0(p) *temp2(p)
            betadl = (1.+avmu(p)*twostext(p))/(omegal*avmu(p)*
     &         twostext(p))*asu
            betail = 0.5 *((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) 
     &         * ((1.+chil(p))/2.)**2) / omegal

            ! Adjust omega, betad, and betai for intercepted snow

            !if (t_veg(p) > tfrz) then                             !no snow
               tmp0 = omegal
               tmp1 = betadl
               tmp2 = betail
            !else
            !   tmp0 =   (1.-fwet(p))*omegal        + fwet(p)*omegas(ib)
            !   tmp1 = ( (1.-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0
            !   tmp2 = ( (1.-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0
            !end if
            omega(p,ib) = tmp0           
            betad = tmp1 
            betai = tmp2  

            ! Absorbed, reflected, transmitted fluxes per unit incoming radiation

            b = 1. - omega(p,ib) + omega(p,ib)*betai
            c1 = omega(p,ib)*betai
            tmp0 = avmu(p)*twostext(p)
            d = tmp0 * omega(p,ib)*betad
            f = tmp0 * omega(p,ib)*(1.-betad)
            tmp1 = b*b - c1*c1
            h = sqrt(tmp1) / avmu(p)
            sigma = tmp0*tmp0 - tmp1
            p1 = b + avmu(p)*h
            p2 = b - avmu(p)*h
            p3 = b + tmp0
            p4 = b - tmp0
          
            ! Determine fluxes for vegetated pft for unit incoming direct 
            ! Loop over incoming direct and incoming diffuse
            ! 0=unit incoming direct; 1=unit incoming diffuse

            ! ic = 0 unit incoming direct flux
            ! ========================================

            u1 = b - c1/albgrd(c,ib)
            u2 = b - c1*albgrd(c,ib)
            u3 = f + c1*albgrd(c,ib)

            tmp2 = u1 - avmu(p)*h
            tmp3 = u1 + avmu(p)*h
            tmp4 = u2 + avmu(p)*h
            tmp5 = u2 - avmu(p)*h
            h1 = -d*p4 - c1*f
            tmp6 = d - h1*p3/sigma
            h4 = -f*p3 - c1*d
            tmp8 = h4/sigma

            delta_z = height_levels(2)-height_levels(1)
            vai(p) = sum(rdfp(1:N_height_level))*delta_z
            ! PET, 3/1/04: added this test to avoid floating point errors in exp()
            t1 = min(h*vai(p), 40.d0)
            s1 = exp(-t1)
            t1 = min(twostext(p)*vai(p), 40.d0)
            s2 = exp(-t1)
          
            d1 = p1*tmp2/s1 - p2*tmp3*s1
            d2 = tmp4/s1 - tmp5*s1          
            tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
            h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
            h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
            tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
            h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
            h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
            h7 = (c1*tmp2) / (d1*s1)
            h8 = (-c1*tmp3*s1) / d1
            h9 = tmp4 / (d2*s1)
            h10 = (-tmp5*s1) / d2


            ! Flux reflected by vegetation (ic = 0)
            albd(p,ib) = h1/sigma + h2 + h3

            do ilevel=1, N_height_level
               evai(ilevel) = sum(rdfp(ilevel:N_height_level))*delta_z   ! cumulative clumping factor * foliage profile

               t1 = min(h*evai(ilevel), 40.d0)
               s1 = exp(-t1)
               t1 = min(twostext(p)*evai(ilevel), 40.d0)
               s2 = exp(-t1)

               ! Downward direct and diffuse fluxes below vegetation (ic = 0)

               ftdd(ilevel,ib) = s2
               ftid(ilevel,ib) = h4*s2/sigma + h5*s1 + h6/s1
               frid(ilevel,ib) = h1*s2/sigma + h2*s1 + h3/s1

               ! Flux absorbed by vegetation (ic = 0)

               fabd(ilevel,ib) = 1. - albd(p,ib) 
     &             - (1.-albgrd(c,ib))*ftdd(ilevel,ib) 
     &             - (1.-albgri(c,ib))*ftid(ilevel,ib)
          
            end do
          

            ! ic = 1 unit incoming diffuse
            ! ========================================

            u1 = b - c1/albgri(c,ib)
            u2 = b - c1*albgri(c,ib)
            u3 = f + c1*albgri(c,ib)

            tmp2 = u1 - avmu(p)*h
            tmp3 = u1 + avmu(p)*h
            tmp4 = u2 + avmu(p)*h
            tmp5 = u2 - avmu(p)*h
            h1 = -d*p4 - c1*f
            tmp6 = d - h1*p3/sigma
            h4 = -f*p3 - c1*d
            tmp8 = h4/sigma

            vai(p) = sum(rifp(1:N_height_level))*delta_z
            t1 = min(h*vai(p), 40.d0)
            s1 = exp(-t1)
            t1 = min(twostext(p)*vai(p), 40.d0)
            s2 = exp(-t1)

            d1 = p1*tmp2/s1 - p2*tmp3*s1
            d2 = tmp4/s1 - tmp5*s1
            tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
            h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
            h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
            tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
            h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
            h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
            h7 = (c1*tmp2) / (d1*s1)
            h8 = (-c1*tmp3*s1) / d1
            h9 = tmp4 / (d2*s1)
            h10 = (-tmp5*s1) / d2

            ! Flux reflected by vegetation
  
            albi(p,ib) = h7 + h8

            do ilevel=1, N_height_level
               evai(ilevel) = sum(rifp(ilevel:N_height_level))*delta_z   ! cumulative clumping factor * foliage profile

               t1 = min(h*evai(ilevel), 40.d0)
               s1 = exp(-t1)
               t1 = min(twostext(p)*evai(ilevel), 40.d0)
               s2 = exp(-t1)

               ! Downward direct and diffuse fluxes below vegetation

               ftii(ilevel,ib) = h9*s1 + h10/s1
               frii(ilevel,ib) = h7*s1 + h8/s1

               ! Flux absorbed by vegetation

               fabi(ilevel,ib) = 1. - albi(p,ib) 
     &            - (1.-albgri(c,ib))*ftii(ilevel,ib)
	     
	       ! Sunlit/shaded intensity
	     
               isha(ilevel,ib) = (1. - omegal) * (Id * (ftid(ilevel,ib) 
     &            + frid(ilevel,ib)) + Ii * (ftii(ilevel,ib) 
     &            + frii(ilevel,ib)))
               isun(ilevel,ib) = (1. - omegal) * twostext(p) * Id 
     &            + isha(ilevel,ib)
             
            end do   ! end of foliage level loop
         end do   ! end of pft loop
      end do   ! end of radiation band loop
    
      end subroutine TwoStream

! ----------------------------------------------------------

      function T(fp, delta_z, szn)

      ! return a transmittance array based on the size of input fp

      real*8, dimension(:), pointer :: fp
      real*8 :: delta_z, szn
      real*8, dimension(size(fp)) :: T
   
      ! local vars
      integer :: N_height_level, ilevel

      N_height_level = size(fp)
      do ilevel=1, N_height_level
        T(ilevel) = exp((-sum(fp(ilevel:N_height_level))*delta_z) 
     &     /cos(szn))
      end do

      end function 
! ----------------------------------------------------------
      subroutine check_inputs(gin)

      ! check the validity of inputs
      type(gort_input) :: gin                  ! gort input

      real*8 :: dens_tree            ! Crown count density (m-2)
      real*8 :: dens_foliage         ! foliage area volume density of a single crown(m-1)
      real*8 :: h1, h2               ! lower and upper bound of crown centers (m)
      real*8 :: delta_z              ! height increase for each height level
      real*8 :: horz_radius          ! horizontal crown radius (m)
      real*8 :: vert_radius          ! vertical crown radius (m)
      real*8 :: zenith               ! solar zenith angle (degree)

      dens_tree = gin%dens_tree
      dens_foliage = gin%dens_foliage
      h1 = gin%h1
      h2 = gin%h2
      delta_z = gin%delta_z
      horz_radius = gin%horz_radius
      vert_radius = gin%vert_radius
      zenith = gin%zenith

      ! check density
      if ( (dens_tree <= 0.) .or. (dens_foliage <= 0.)) then
        write(*,*) 'Density must be greater than 0. STOPPING ...'
        stop
      end if

      ! check the height
      if ((h1 <= 0.) .or. (h2 <= 0.)) then
        write(*,*) 'Center of crown height must be greater than 0. 
     & 	     STOPPING ...'
        stop
      end if

      if (h1 >= h2) then
        write(*,*) 'h1 must be less than h2. STOPPING  ...'
        stop
      end if

      if ( h1 < vert_radius ) then
        write(*,*) 'h1 must be greater than vert_radius. STOPPING...'
        stop
      end if

      ! check radius
      if ((horz_radius <= 0.) .or. (vert_radius <= 0.)) then
        write(*,*) 'Radius must be greater than 0. STOPPING ...'
        stop
      end if

      ! check angle
      if ((zenith < 0.) .or. (zenith > 1.)) then
        write (*,*) 'Solar zenith angle, zenith, is in radius, and must 
     &	     be within [0.,1.]. STOPPING ...'
        stop
      end if

      end subroutine check_inputs

! ----------------------------------------------------------
      subroutine get_height_level(gin, height_levels)

        ! compute the height levels between h1 and h2
        ! NOTE : here last layer could be less than delta_z.
        !        It may be necessary to return the dz at the last layer to
        !        ensure the accuracy of calculation, otherwise, we can make
        !        sure (h2-h1) can be integer divided by delta_z.

      type(gort_input), intent(in) :: gin
      real*8, dimension(:), pointer :: height_levels

      ! local
      real*8 :: h1, h2, vert_radius, delta_z
      integer :: n_levels, i, temp1
      real*8 :: jh1, jh2, temp2, dif

      ! -------------------------------------------


      h1 = gin%h1
      h2 = gin%h2
      vert_radius = gin%vert_radius
      delta_z = gin%delta_z

      jh1 = h1 - vert_radius
      jh2 = h2 + vert_radius

      ! Added the rounding checking for different machines.

      temp1 = ceiling((jh2-jh1)/delta_z)
      temp2 = (jh2-jh1)/delta_z

      if ((temp1-temp2) > 0.999) then
        n_levels = ceiling((jh2-jh1)/delta_z)
      else
        n_levels = ceiling((jh2-jh1)/delta_z) + 1
      end if


      allocate(height_levels(n_levels))

      do i =1, n_levels-1
        height_levels(i) = jh1 + delta_z * (i-1)
      end do

      height_levels(n_levels) = jh2

      end subroutine get_height_level

! ----------------------------------------------------------

      subroutine get_foliage_profile(gin, height_levels, fp)

      type(gort_input) :: gin
      real*8, dimension(:), intent(in) :: height_levels
      real*8, dimension(:), intent(inout) :: fp

      ! local
      real*8 :: h1, h2, vert_radius, tmp
      real*8 :: z, jh11, jh12, jh21, jh22, hdif
      integer :: n_levels, i

      real*8, dimension(:), allocatable :: xx

      ! ---------------------------------

      h1 = gin%h1
      h2 = gin%h2
      vert_radius = gin%vert_radius

      tmp=gin%dens_tree * gin%dens_foliage * PI * (gin%horz_radius ** 2)

      jh11 = h1 - vert_radius
      jh12 = h1 + vert_radius
      jh21 = h2 - vert_radius
      jh22 = h2 + vert_radius
      hdif = h2 - h1

      n_levels = size(height_levels)
      allocate(xx(n_levels))

      ! define xx based on the location of each height level
      do i = 1, n_levels

        z = height_levels(i)
        if (jh12 < jh21) then
           if ((z >= jh11) .and. (z <= jh12))   
     &         xx(i) = height_function1(z, h1, vert_radius, hdif)
           if ((z > jh12) .and. (z <= jh21))   
     &         xx(i) = height_function3(hdif, vert_radius)
           if ((z > jh21) .and. (z <= jh22))   
     &         xx(i) = height_function2(z, h2, vert_radius, hdif)
        else ! (jh12 >= jh21)
           if ((z >= jh11) .and. (z <= jh21))   
     &         xx(i) = height_function1(z, h1, vert_radius, hdif)
           if ((z > jh21) .and. (z <= jh12))   
     &         xx(i) = height_function4(z, h1, h2, vert_radius)
           if ((z > jh12) .and. (z <= jh22))   
     &         xx(i) = height_function2(z, h2, vert_radius, hdif)
        endif

      fp(i) = tmp * xx(i)

      end do

      deallocate(xx)

      end subroutine get_foliage_profile

! ----------------------------------------------------
 
      subroutine convolute(all_convolute_input,      ! input
     &   fp, rdfp, rifp, efp, height_levelC)         ! convoluted output


      ! NOTE: This subroutine convolute multiple profiles into one.
      !       The algorithm used here is mapping each profile to the desired height levels
      !       and then convolute (either sum or multiply) profiles together.
      !       The height_levelC will be based upon the input for desired height, however,
      !       it will be adjusted to cover all the participating profiles.
      !
      ! note the following params are going to be decided by the input profiles
      !  real :: hc_low       ! Low end for the convoluted profile
      !  real :: hc_high      ! High end for the convoluted profile
      !  real :: delta_z      ! Delta z for the convoluted profile
      ! --------------------------------------------------------------------------


      ! Input
      type(profile_params), dimension(:), pointer :: all_convolute_input


      ! Output
      real*8, dimension(:), pointer   :: fp         ! foliage profile
      real*8, dimension(:), pointer   :: rdfp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: rifp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: efp        ! effective foliage profile
      real*8, dimension(:), pointer   :: height_levelC  ! Convolute profile height levels


      ! Local variables
      real*8, dimension(:), allocatable :: arrthl, arrthh, arrdz
      real*8 :: hc_low, hc_high, delta_z
      integer :: n_profiles, iprofile, n_convolute_levels
      type(profile_params), pointer :: interpolate_out

      !---------------------------------------------------------------------

      print *, 'I am in convolute now....'

      n_profiles = size(all_convolute_input)

      print *, all_convolute_input(1)%fp(10)
      print *, all_convolute_input(2)%fp(10)


      ! Get the height levels for the convoluted profile
      !----------------------------------------------------------------------------------------------
      allocate(arrthl(n_profiles))
      allocate(arrthh(n_profiles))
      allocate(arrdz(n_profiles))

      ! compare with all profiles and adjust the low end and high end

      do iprofile = 1, n_profiles
        arrthl(iprofile) = minval(
     &	    all_convolute_input(iprofile)%height_levels)
        arrthh(iprofile) = maxval(
     &	    all_convolute_input(iprofile)%height_levels)
        arrdz(iprofile) = all_convolute_input(iprofile)%height_levels(2) 
     &	    - all_convolute_input(iprofile)%height_levels(1)
      end do

      print *, arrthl
      print *, arrthh
      print *, arrdz

      hc_low=minval(arrthl)
      hc_high=maxval(arrthh)
      delta_z=minval(arrdz)
      print *,hc_low, hc_high, delta_z

      deallocate(arrthl, arrthh, arrdz)

      print *, 'Getting convolute height levels ...'
      call get_convolute_height_level(hc_low, hc_high, delta_z,
     &   height_levelC)
      n_convolute_levels=size(height_levelC)


      ! allocate output variables
      print *, 'Allocating for the convoluted outputs ....'

      allocate(fp(n_convolute_levels))
      allocate(rdfp(n_convolute_levels))
      allocate(rifp(n_convolute_levels))
      allocate(efp(n_convolute_levels))

      ! allocate temporary interpolate_out
      print *, 'Allocating for the interpolation output ....'
      allocate(interpolate_out)
      allocate(interpolate_out%fp(n_convolute_levels))
      allocate(interpolate_out%rdfp(n_convolute_levels))
      allocate(interpolate_out%rifp(n_convolute_levels))
      allocate(interpolate_out%efp(n_convolute_levels))
      allocate(interpolate_out%height_levels(n_convolute_levels))

      interpolate_out%height_levels = height_levelC


      ! interpolate profiles to the convoluted height level
      !---------------------------------------------------------

      ! interpolate profile1 to the convoluted height levels
      do iprofile = 1, n_profiles


         print *, 'Interpolating for profile No. ', iprofile
         call interpolate_profile (all_convolute_input(iprofile),   
     &      interpolate_out)
         print *, interpolate_out%fp(10)

         print *, 'Reassign or convolute the output'

         if (iprofile == 1) then
            fp = interpolate_out%fp
            rdfp = interpolate_out%rdfp
            rifp = interpolate_out%rifp
            efp = interpolate_out%efp
         else
            fp = fp + interpolate_out%fp
            rdfp = rdfp + interpolate_out%rdfp
            rifp = rifp + interpolate_out%rifp
            efp = efp + interpolate_out%efp
         end if
            print *, fp(10)

      end do

      deallocate(interpolate_out%fp)
      deallocate(interpolate_out%rdfp)
      deallocate(interpolate_out%rifp)
      deallocate(interpolate_out%efp)
      deallocate(interpolate_out%height_levels)
      deallocate(interpolate_out)


      return


      end subroutine convolute


 ! ----------------------------------------------------------------------------------

      subroutine interpolate_profile (interpolate_in, interpolate_out)

      ! Input
      type(profile_params), intent(in) :: interpolate_in          ! data for input profile

      ! Output
      type(profile_params), pointer :: interpolate_out            ! interpolated profile

      ! local variables
      integer :: i, j
      integer :: nlevel_in, nlevel_out
      real*8 :: h, hin_low, hin_high
      integer :: id1, id2
      integer, dimension(1) :: id
      real*8, dimension(:), pointer :: height_in, height_out      ! height levels for input, output profile
      real*8, dimension(:), pointer :: diff

      ! ---------------------------------------
      print *, interpolate_in%fp(10)

      height_in => interpolate_in%height_levels
      height_out => interpolate_out%height_levels

      nlevel_in = size(height_in)
      nlevel_out = size (height_out)
      hin_low = height_in(1)
      hin_high = height_in(nlevel_in)
      allocate(diff(nlevel_in))


      ! loop through each height level for the convoluted height levels.
      ! Based upon the height comparison to get the interpolated profile

      do i=1, nlevel_out

        h = height_out(i)

        if ( (h > hin_low) .and. (h < hin_high)) then

          ! find the closest two height level
          do j=1, nlevel_in
            diff(j) = (height_in(j) - h)
          end do
          id = minloc ( abs(diff) ) ! find the location that has the height closes to h
          if (diff(id(1)) > 0) then
            id1 = id(1) - 1
            id2 = id(1)
          else
            id1 = id(1)
            id2 = id(1) + 1
          end if

          ! use the data from the above two points to interpolate for the desired height
       

          interpolate_out%fp(i) = interpolate_in%fp(id1) + (h - 
     & 	     height_in(id1)) * (interpolate_in%fp(id2) - 
     &       interpolate_in%fp(id1)) / (height_in(id2) - height_in(id1))
          interpolate_out%rdfp(i) = interpolate_in%rdfp(id1) + (h - 
     &	     height_in(id1)) * (interpolate_in%rdfp(id2) -  
     &       interpolate_in%rdfp(id1)) / (height_in(id2) - 
     &       height_in(id1))
          interpolate_out%rifp(i) = interpolate_in%rifp(id1) + (h - 
     &	     height_in(id1)) * (interpolate_in%rifp(id2) -
     &       interpolate_in%rifp(id1)) / (height_in(id2) - 
     &       height_in(id1))
          interpolate_out%efp(i) = interpolate_in%efp(id1) + (h - 
     &	     height_in(id1)) * (interpolate_in%efp(id2) -   
     &       interpolate_in%efp(id1)) /(height_in(id2) - height_in(id1))
       

        else

          ! just fill in the output profile with the ceiling or bottom value

          if (h <= hin_low) then
            interpolate_out%fp(i) = interpolate_in%fp(1)
            interpolate_out%rdfp(i) = interpolate_in%rdfp(1)
            interpolate_out%rifp(i) = interpolate_in%rifp(1)
            interpolate_out%efp(i) = interpolate_in%efp(1)
          end if

          if (h >= hin_high) then
            interpolate_out%fp(i) = interpolate_in%fp(nlevel_in)
            interpolate_out%rdfp(i) = interpolate_in%rdfp(nlevel_in)
            interpolate_out%rifp(i) = interpolate_in%rifp(nlevel_in)
            interpolate_out%efp(i) = interpolate_in%efp(nlevel_in)
          end if

        end if

      end do

      nullify(height_in)
      nullify(height_out)

      deallocate(diff)
      return

      end subroutine interpolate_profile


 !-------------------------------------------------------------------------------

      subroutine get_convolute_height_level(h1, h2, delta_z, 
     &   height_levels)

      ! compute the height levels between h1 and h2
      ! NOTE : here last layer could be less than delta_z.
      !        It may be necessary to return the dz at the last layer to
      !        ensure the accuracy of calculation, otherwise, we can make
      !        sure (h2-h1) can be integer divided by delta_z.

      real*8, intent(in) :: h1, h2, delta_z
      real*8, dimension(:), pointer :: height_levels

      ! local
      integer :: n_levels, i, temp1
      real*8 :: temp2, dif

      ! -------------------------------------------

      ! Added the rounding checking for different machines.

      temp1 = ceiling((h2-h1)/delta_z)
      temp2 = (h2-h1)/delta_z

      if ((temp1-temp2) > 0.999) then
        n_levels = ceiling((h2-h1)/delta_z)
      else
        n_levels = ceiling((h2-h1)/delta_z) + 1
      end if


      allocate(height_levels(n_levels))

      do i =1, n_levels-1
        height_levels(i) = h1 + delta_z * (i-1)
      end do

        height_levels(n_levels) = h2

      end subroutine get_convolute_height_level

! ----------------------------------------------------

      subroutine geo_to_gin(pp, gin)
!@sum calculate gin from Ent input, get mean patch geometry from cohort value,
!@sum then get the gort input at patch level, note the whole Entcell share the 
!@sum same input radiation  
!@sum note in calculating clumping indices, dbh is not needed
      
      type(patch),pointer :: pp
      type(gort_input),dimension(:),pointer :: gin

      !----Local----------------!
      type(cohort),pointer :: cop
      integer :: count, i
      real*8 :: vc
      
#ifdef DEBUG
      print *,"Started geo_to_gin" ! with patch:"
#endif

      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
        pp%crad%heights = 0.d0
        pp%crad%LAI = 0.d0
        pp%crad%GORTclump = 0.d0
        return
      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *,"GORT_clumping: wrong pft = ", pp%tallest%pft
        call patch_print(6,pp,"ERROR ")
        call stop_model("GORT_clumping: wrong pft",255)
      endif
      
      !* LOOP THROUGH COHORTS *!
      cop => pp%tallest
      do while (ASSOCIATED(cop))
        !* Assign vegpar

        if (cop%LAI.gt.0.d0) then
	  count = count + 1
        endif

        cop => cop%shorter
      end do
      
      if (count.gt.0.d0) then
	if (count.eq.1.d0) then ! single cohort
          allocate(gin(1))   
          cop => pp%tallest
      
          gin(1)%dens_tree = cop%n
          gin(1)%h1 = cop%h - cop%crown_dy
          gin(1)%h2 = gin(1)%h1
          gin(1)%horz_radius = cop%crown_dx
          gin(1)%vert_radius = cop%crown_dy
          gin(1)%zenith = acos(pp%cellptr%CosZen)
          gin(1)%delta_z = 0.1
          vc = 1.33333 * PI * gin(1)%horz_radius ** 2 * 
     &	     gin(1)%vert_radius
          gin(1)%dens_foliage = cop%LAI / (vc * gin(1)%dens_tree)
          gin(1)%pft = cop%pft
          gin(1)%LAI = cop%LAI
	else ! multi cohort
          allocate(gin(count)) 
	    
          !* LOOP THROUGH COHORTS *!
          cop => pp%tallest
          i = 0;
          do while (ASSOCIATED(cop))
      
            if (cop%LAI.gt.0.d0) then
	      i = i + 1
              gin(i)%dens_tree = cop%n
              gin(i)%h1 = cop%h - cop%crown_dy
              gin(i)%h2 = gin(i)%h1
              gin(i)%horz_radius = cop%crown_dx
              gin(i)%vert_radius = cop%crown_dy
              gin(i)%zenith = acos(pp%cellptr%CosZen)
              gin(i)%delta_z = 0.1
              vc = 1.33333 * PI * gin(i)%horz_radius ** 2 * 
     &	         gin(i)%vert_radius
              gin(i)%dens_foliage = cop%LAI / (vc * gin(i)%dens_tree)
              gin(i)%pft = cop%pft
	      gin(i)%LAI = cop%LAI
            endif 
	
            cop => cop%shorter
          end do
	endif
	
      else
        print *,"GORT_clumping: no vagetation "
        call patch_print(6,pp,"ERROR ")
        call stop_model("GORT_clumping: no vagetation",255)
      endif
      
      end subroutine geo_to_gin

! ----------------------------------------------------
      function get_analytical_clump(gin)
!@sum calculate clumping index from canopy geometry 

      type(gort_input),intent(in) :: gin
      real*8 :: get_analytical_clump

      ! local variables:
      real*8 :: tr, a, b, c, d, G

      a = (tan(gin%zenith)) ** 2
      b = (gin%vert_radius/gin%horz_radius) ** 2
      c = sqrt((1+a)/(1+a*b))
      G = get_K(cos(gin%zenith),gin%pft)
      tr = G * gin%dens_foliage * gin%horz_radius * c

      d = (1 - (2*tr+1) * exp(-2*tr)) / (2*tr*tr)

      get_analytical_clump = 3*(1-d)/(4*tr)

      end function get_analytical_clump

! ------------------------------------------------------

      function height_function1(z, h, b, hdif)
!@sum the four height functions describes 3D GORT geometry, 
!@sum from Ni-Meister et al., 2001

      real*8,intent(in) :: z, h, b, hdif
      real*8 :: height_function1

      height_function1 = ((b+z-h)**2) * (2*b + h - z) / (3*hdif*b**2)

      end function height_function1

! ------------------------------------------------------

      function height_function2(z, h, b, hdif)

      real*8,intent(in) :: z, h, b, hdif
      real*8 :: height_function2

      height_function2 = ((b+h-z)**2) * (2*b + z - h) / (3*hdif*b**2)

      end function height_function2

! -------------------------------------------------------

      function height_function3(hdif, b)

      real*8,intent(in) :: b, hdif
      real*8 :: height_function3

      height_function3 = (4./3.) * b / hdif

      end function height_function3

! -------------------------------------------------------

      function height_function4(z, h1, h2, b)

      real*8,intent(in) :: z, h1, h2, b
      real*8 :: height_function4

      height_function4 = 1 - ((h2-h1)**2 + 3*(z-h1)*(z-h2)) / (3*b**2)

      end function height_function4

! --------------------------------------------------------

      function get_K(CosZen, pft)
!@sum G is the relative projected area of leaf elements in the direction zenith
!@sum get_K here is actually G, following the definition of Li et al., 1995
!@sum G function comes from Sellers, 1985 

      real*8,intent(in) :: CosZen
      integer,intent(in) :: pft
      real*8 :: get_K

      integer, parameter :: numpft = 16   ! number of plant function types
      integer :: i 
      real*8 :: xl(numpft)       ! ecophys const - leaf/stem orientation index
      data (xl(i),i=1,numpft) /0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, 
     &    0.25, 0.01, 0.25, 0.25, -0.30, -0.30, -0.30, -0.30, -0.30/
      real*8 :: phi1,phi2,chil   
                                    
      chil = min( max(xl(pft), -0.4d0), 0.6d0 )
      if (abs(chil) <= 0.01) chil = 0.01
      phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
      phi2 = 0.877 * (1.-2.*phi1)
      get_K = phi1 + phi2*CosZen

      end function get_K

! --------------------------------------------------------

      function get_K2(gin)
!@sum G is the relative projected area of leaf elements in the direction zenith
!@sum get_K2 here is actually G for single or multiple layers

      type(gort_input), dimension(:), pointer :: gin  ! input structure for gort model
      real*8 :: get_K2
      real*8 :: szn

      ! local variables
      integer :: num_profiles, iprofile
      real*8, dimension(:), pointer :: clump, G
      real*8 :: eLAI, rLAI
      
      num_profiles = size(gin)
      allocate(clump(num_profiles))
      allocate(G(num_profiles))
      eLAI = 0
      rLAI = 0
      
      do iprofile = 1, num_profiles
         szn = gin(iprofile)%zenith
         G(iprofile) = get_K(cos(szn),gin(iprofile)%pft)
         clump(iprofile) = get_analytical_clump(gin(iprofile))
         eLAI = eLAI + clump(iprofile) * G(iprofile) * gin(iprofile)%LAI
         rLAI = rLAI + clump(iprofile) * gin(iprofile)%LAI
      end do
      get_K2 = eLAI / rLAI

      end function get_K2

      end module canopyrad
