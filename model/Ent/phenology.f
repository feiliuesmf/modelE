      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types

      implicit none


      public phenology_update, litter

      contains
      !*********************************************************************
      
      subroutine phenology_update(dtsec, time, entcell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: entcell

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      
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
