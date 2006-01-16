      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types
      use ent_const, only : JEQUATOR

      implicit none


      public phenology_update, litter

      contains
      !*********************************************************************
      subroutine phenology_update_cell(dtsec,time,ecp)
!@sum Updates phenology for all patches in an entcell.
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcell),pointer :: ecp
      !------local--------
      type(patch),pointer :: pp

      pp = ecp%youngest
      do while (ASSOCIATED(pp))
        call phenology_update(dtsec, time, pp)
      end do
      end subroutine phenology_update_cell
      !*********************************************************************
      
      subroutine phenology_update(dtsec, time, pp)
!@sum Update phenology for a patch.
      use ent_GCM_coupler, only : GISS_calc_lai
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(patch),pointer :: pp
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: laip  !patch-level summary of LAI
      real*8 :: laig  !entcell grid-level summary of LAI
      integer :: hemi !-1: S.hemisphere, 1: N.hemisphere
      !---------------------------------------------------------------
      !* DUMMY TEMPORARY GISS 
      !* JUST UPDATES LAI USING MATTHEWS PRESCRIPTION
      if (allocated(pp)) then
        if (pp%cellptr%latj < JEQUATOR) then 
          hemi = -1
        else 
          hemi = 1
        end if
        laip = 0.0
        cop = pp%tallest
        do while (allocated(cop))
          ccount = ccount + 1
          cop%lai = GISS_calc_lai(cop%pft, time%jday, hemi)
          laip = laip + cop%LAI
          cop = cop%shorter
        end do
        pp%LAI_p = laip
      endif

      end subroutine phenology_update


      !*********************************************************************
      subroutine litter(dtsec, time, entcell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: entcell

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine litter
      !*********************************************************************

      end module phenology
