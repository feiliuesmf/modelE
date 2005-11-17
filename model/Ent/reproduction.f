      module reproduction
!@sum Routines to calculate reproduction rates in an entcell.

      use ent_types

      implicit none


      contains
      !*********************************************************************

      subroutine reproduction_calc(dtsec, time, entcell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: entcell

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not needed in GISS replication test.

      end subroutine reproduction_calc

      !*********************************************************************
      
      end module reproduction
