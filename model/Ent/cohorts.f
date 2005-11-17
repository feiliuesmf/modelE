      module cohorts
!@sum Routines to organize cohorts within an entcell.


      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_cohort(pp,pft,n, h,
     &     crown_dx, crown_dy,dbh, root_d,LAI,clump,
     &     LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     GCANOPY, GPP, NPP, R_growth,R_maint,
     &     N_up, C_litter,N_litter,C_to_Nfix

      type(patch),pointer :: pp
      integer :: pft, n
      real*8 :: h,
     &     crown_dx, crown_dy,dbh, root_d,LAI,clump,
     &     LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     GCANOPY, GPP, NPP, R_growth,R_maint,
     &     N_up, C_litter,N_litter,C_to_Nfix
      !------------------
      type(cohort),pointer :: cop, csp, newc

      !If !allocated(pp) then ERROR

      !* If pp has cohorts, then insert, else allocate.
      !* Insert at correct height in list, starting from shortest,
      !* since most new cohorts will be shortest.
        call allocate(newc)
        newc%pptr = pp
        call assign_cohort(newc,pft,n, h,
     &       crown_dx, crown_dy,dbh, root_d,LAI,clump,
     &       LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &       C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &       GCANOPY, GPP, NPP, R_growth,R_maint,
     &       N_up, C_litter,N_litter,C_to_Nfix)

        if (allocated(pp%shortest)) then !A. There are other cohorts.
          if (pp%shortest%h.ge.newc%h) then !newc is shortest
            pp%shortest%shorter = newc
            newc%taller = pp%shortest
            pp%shortest = newc
          else if (pp%tallest%h.lt.newc%h) then !newc is tallest
            pp%tallest%taller = newc
            newc%shorter = pp%tallest
            pp%tallest = newc
          else !newc is neither tallest nor shortest
            cop = pp%shortest
            do while (cop%h.lt.newc%h) !find next taller cohort
              cop = cop%taller
            end do
            newc%taller = cop
            newc%shorter = cop%shorter
            newc%shorter%taller = newc
            cop%shorter = new
          end if
          !Now parse through csp's
          csp = newc%taller
          if (allocated(csp)) then
            do while (csp%pft.ne.newc%pft)
              if (allocated(csp%taller).and.(csp%pft.ne.newc%pft)) then
                csp = csp%taller
              else
                exit !exit loop
              end if
            end do
            if (csp%pft.eq.newc%pft) then !1.newc is not tallest csp
              newc%csptaller = csp
              newc%cspshorter = csp%cspshorter
              csp%cspshorter = newc
            else !2. no taller con-specifics
              call nullify(newc%csptaller)
            end if
          else  !3. no taller con-specifics
            call nullify(newc%csptaller)
          end if
          if (unallocated(newc%cspshorter)) then !Case 1 did not hold
            csp = newc%shorter
            if (allocated(csp)) then
              do while (csp%pft.ne.newc%pft))
                if (allocated(csp%shorter).and.
     &               (csp%pft.ne.newc%pft)) then
                  csp = csp%shorter
                else
                  exit
                end if
              end do
              if (csp%pft.eq.newc%pft) then !4. newc is not shortest csp
                newc%cspshorter = csp
                newc%csptaller = csp%csptaller
                csp%csptaller = newc
              else !5. no shorter con-specifics
                call nullify(newc%cspshorter)
              end if
            else !6. no shorter con-specifics
              call nullify(newc%cspshorter)
            end if
          end if
        else !B. newc is the only cohort
          pp%tallest = newc
          pp%shortest = newc
          call nullify(newc%taller) 
          call nullify(newc%shorter)
          call nullify(newc%csptaller)
          call nullify(newc%cspshorter)
        end if
      end subroutine insert_cohort
      !*********************************************************************
      
      subroutine assign_cohort(cop,pft,n, h,
     &     crown_dx, crown_dy,dbh, root_d,LAI,clump,
     &     LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     GCANOPY, GPP, NPP, R_growth,R_maint,
     &     N_up, C_litter,N_litter,C_to_Nfix)
      !Given cohort's characteristics, assign to cohort data variable.

      type(cohort),pointer :: cop
      integer :: pft, n
      real*8 :: h,
     &     crown_dx, crown_dy,dbh, root_d,LAI,clump,
     &     LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     GCANOPY, GPP, NPP, R_growth,R_maint,
     &     N_up, C_litter,N_litter,C_to_Nfix

      cop%pft = pft
      cop%n = n
      cop%h = h
      cop%crown_dx =  crown_dx 
      cop%crown_dy =  crown_dy
      cop%dbh =  dbh 
      cop%root_d = root_d
      cop%LAI = LAI
      cop%clump = clump
      cop%LMA =  LMA 
      cop%C_fol =  C_fol 
      cop%N_fol = N_fol
      cop%C_sw =  C_sw 
      cop%N_sw =  N_sw 
      cop%C_hw =  C_hw 
      cop%N_hw = N_hw
      cop%C_lab =  C_lab 
      cop%N_lab =  N_lab 
      cop%C_froot = C_froot
      cop%N_froot =  N_froot 
      cop%C_croot =  C_croot 
      cop%N_croot = N_croot
      cop%GCANOPY =  GCANOPY 
      cop%GPP =  GPP 
      cop%NPP =  NPP 
      cop%R_growth = R_growth
      cop%R_maint = R_maint
      cop%N_up =  N_up 
      cop%N_litter = N_litter
      cop%C_to_N = fixC_to_Nfix

      end subroutine assign_cohort
      !*********************************************************************
      subroutine init_cohort_defaults(cop,pnum)
      !Initialize a cohort with default values
      type(cohort),pointer :: cop
      integer :: pnum

            cop%pft = pnum
            cop%n = 1   !## Dummy ##!
!            cop%pptr = pp
!            call nullify(cop%taller) !Only one cohort
!            call nullify(cop%shorter)
!            call nullify(cop%csptaller)
!            call nullify(cop%cspshorter)            

            cop%nm = pfpar(pnum,5)  !## This is nf.  Need to calc nm ##!
            cop%Ntot = !## Get from GISS GCM ##!

            !* Individual plant properties *!
            !* GEOMETRY *!
            cop%h
            cop%crown_dx
            cop%crown_dy
            cop%dbh
            cop%root_d
            cop%LAI
            cop%clump

            !* BIOMASS POOLS *!
            cop%LMA
            cop%C_fol
            cop%N_fol
            cop%C_sw
            cop%N_sw
            cop%C_lab
            cop%N_lab
            cop%C_froot
            cop%N_froot
            cop%C_croot
            cop%N_croot

            !* FLUXES *!
            cop%GCANOPY
            cop%GPP
            cop%NPP
            cop%R_growth
            cop%R_maint
            cop%N_up
            cop%C_litter
            cop%N_litter
            cop%C_to_Nfix

            !* REPRODUCTION *!
            !cop%
      end subroutine init_cohort_defaults


      !*********************************************************************
      !*********************************************************************
       
      subroutine reorganize_cohorts(entcell)
      type(entcelltype) :: entcell

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not need in GISS replication test.
      
      end subroutine reorganize_cohorts


      !*********************************************************************
  

      end module cohorts
