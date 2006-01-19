      module growthallometry
!@sum Routines to calculate uptake and allocation of nutrients within an
!@sum individual plant, growth, and allometry.

      use ent_types

      implicit none


      contains
      !*********************************************************************
      
      subroutine uptake_N(dtsec, pp)
      real*8 :: dtsec           !dt in seconds
      type(patch),pointer :: pp

      !*****************************************************************
      !*              FILL IN CODE                                     *
      !*****************************************************************

      
      end subroutine uptake_N
      !*********************************************************************

      subroutine init_rootdistr(froot, pft)
      use ent_GCM_coupler, only : GISS_calc_froot
      real*8 :: froot(N_DEPTH)
      integer :: pft

      !GISS TEMPORARY
      call GISS_calc_froot(froot, pft)

      end subroutine init_rootdistr
      !*********************************************************************
      end module growthallometry
