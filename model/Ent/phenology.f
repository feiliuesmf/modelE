      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types

      implicit none


      public phenology_update, litter

      contains
      !*********************************************************************
      
      subroutine phenology_update(dtsec, time, entcell)
      use ent_GCM_coupler, only : GISS_calc_lai
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: entcell
      !-------local-----
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      real*8 :: laip  !patch-level summary of LAI
      real*8 :: laig  !entcell grid-level summary of LAI
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !* DUMMY TEMPORARY GISS 
      !* JUST UPDATES LAI USING MATTHEWS PRESCRIPTION
      if (entcell%latj < JEQUATOR) then 
        hemi = -1
      else 
        hemi = 1
      end if
      laig = 0.0
      pp = entcell%oldest
      do while (allocated(pp))
        laip = 0.0
        cop = pp%tallest
        do while (allocated(cop))
          call GISS_calc_lai(cop%LAI,cop%pft, time%jday, hemi)
          laip = laip + cop%LAI
          cop = cop%shorter
        end do
        pp%LAI_p = laip
        laig = laig + pp%LAI_p
        pp = pp%younger
      end do
      entcell%LAI = laig
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
