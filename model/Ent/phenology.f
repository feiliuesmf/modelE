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
      
      subroutine phenology_update(dtsec, tt, pp)
!@sum Update phenology for a patch.
      use ent_GCM_coupler, only : GISS_phenology_update
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: tt  !Greenwich Mean Time
      type(patch),pointer :: pp

      !---------------------------------------------------------------
      !* GISS VERSION:  JUST UPDATES LAI USING MATTHEWS PRESCRIPTION
      call GISS_phenology_update(tt%jday, pp)

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
