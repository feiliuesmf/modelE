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
        call GORT_clumping(pptr)
        pptr = pptr%older
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      end do

      end subroutine recalc_radpar_OLD


      !*********************************************************************
      subroutine GORT_clumping(pptr)
!@sum Calculate the GORT clumping index in canopy layers and save into
!@sum variable ppt%crad
      type(patch),pointer :: pptr
      !--------------------------
!      type(canradtype) :: crad

      !DEFINE LAYERS FOR LAI AND CLUMPING INDICES
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------


        !Allocate leaf area/crown volume within canopy layers
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      !Calculate GORT clumping index for each layer 
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      !* Return (pptr%crad)
      end subroutine GORT_clumping

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
      subroutine get_canopy_rad(pptr,h)
!@sum Retrieve incident light at height h from pptr%crad
      type(patch),pointer :: pptr
      real*8 :: h               !Height in canopy
      !--------------------------
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine get_canopy_rad
      !*********************************************************************

      end module canopyrad
