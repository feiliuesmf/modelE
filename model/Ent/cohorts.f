      module cohorts
!@sum Routines to organize cohorts within an entcell.

      use ent_const
      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_cohort(pp,pft,n, h,
     &     nm,LAI,
     &     crown_dx, crown_dy,dbh, clump,LMA, root_d,fracroot,
     &     C_fol, N_fol, C_sw, N_sw, C_hw, N_hw, C_lab, N_lab,
     &     C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix, 
     &     phenofactor_c, phenofactor_d, phenofactor, dphdt, CB_d) !KIM - 5 vars for phenology
      !NYK - stressH2O and stressH2Ol depend on soil moisture, are calculated in biophysics.f.

      type(patch),pointer :: pp
      integer :: pft
      real*8 :: n, h, nm, LAI,
     &     crown_dx, crown_dy,dbh, clump
      real*8 :: root_d,fracroot(N_DEPTH)
      real*8 :: LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix,
     &     phenofactor_c, phenofactor_d, phenofactor, dphdt, CB_d
!     &     stressH2O, stressH2O(N_DEPTH) !No need to assign biophysical values initialized in cohort_construct.
      !------------------
      type(cohort),pointer :: cop, csp, newc

      !If !ASSOCIATED(pp) then ERROR

      !* If pp has cohorts, then insert, else allocate.
      !* Insert at correct height in list, starting from shortest,
      !* since most new cohorts will be shortest.
        !!ALLOCATE(newc)
        !!newc%pptr = pp
      call cohort_construct(newc, pp, pft)
        call assign_cohort(newc,pft,n, h, nm, LAI,
     &       crown_dx, crown_dy,dbh, clump, LMA, root_d,fracroot,
     &       C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &       C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &       Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &       N_up, C_to_Nfix,
     &       phenofactor_c, phenofactor_d, phenofactor, dphdt, CB_d)

        newc%Ntot = nm*LAI

        ! >>>>
        ! Nancy, make sure that you distinguish between = and =>
        ! in the following code
        !changed '=' to '=>' to avoid memory faults -PK 9/28/07
        if (ASSOCIATED(pp%shortest)) then !A. There are other cohorts.
          if (pp%shortest%h.ge.newc%h) then !newc is shortest
            pp%shortest%shorter => newc  !changed = to => -PK 9/28/07 
            newc%taller => pp%shortest
            pp%shortest => newc
          else if (pp%tallest%h.lt.newc%h) then !newc is tallest
            pp%tallest%taller => newc  
            newc%shorter => pp%tallest
            pp%tallest => newc
          else !newc is neither tallest nor shortest
            cop => pp%shortest
            do while (cop%h.lt.newc%h) !find next taller cohort
              cop => cop%taller
            end do
            newc%taller => cop
            newc%shorter => cop%shorter
            newc%shorter%taller => newc
            cop%shorter => newc
          end if
          !Now parse through csp's
          csp => newc%taller
          if (ASSOCIATED(csp)) then
            do while (csp%pft.ne.newc%pft)
              if (ASSOCIATED(csp%taller).and.(csp%pft.ne.newc%pft)) then
                csp => csp%taller
              else
                exit !exit loop
              end if
            end do
            if (csp%pft.eq.newc%pft) then !1.newc is not tallest csp
              newc%csptaller => csp
              newc%cspshorter => csp%cspshorter
              csp%cspshorter => newc
            else !2. no taller con-specifics
              nullify(newc%csptaller)
            end if
          else  !3. no taller con-specifics
            nullify(newc%csptaller)
          end if
          if (.NOT.ASSOCIATED(newc%cspshorter)) then !Case 1 did not hold
            csp => newc%shorter
            if (ASSOCIATED(csp)) then
              do while (csp%pft.ne.newc%pft)
                if (ASSOCIATED(csp%shorter).and.
     &               (csp%pft.ne.newc%pft)) then
                  csp => csp%shorter
                else
                  exit
                end if
              end do
              if (csp%pft.eq.newc%pft) then !4. newc is not shortest csp
                newc%cspshorter => csp
                newc%csptaller => csp%csptaller
                csp%csptaller => newc
              else !5. no shorter con-specifics
                nullify(newc%cspshorter)
              end if
            else !6. no shorter con-specifics
              nullify(newc%cspshorter)
            end if
          end if
        else !B. newc is the only cohort
          pp%tallest => newc
          pp%shortest => newc
          ! actually following are already nullified in cohort_construct
          !nullify(newc%taller) 
          !nullify(newc%shorter)
          !nullify(newc%csptaller)
          !nullify(newc%cspshorter)
        end if
      end subroutine insert_cohort
      !*********************************************************************
      
      subroutine assign_cohort(cop,pft,n, h, nm, LAI,
     &     crown_dx, crown_dy,dbh, clump,LMA,root_d,fracroot,
     &     C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix,
     &     phenofactor_c, phenofactor_d, phenofactor, dphdt, CB_d)
!     &     stressH2O, stressH2Ol)

      !Given cohort's characteristics, assign to cohort data variable.

      type(cohort) :: cop
      integer :: pft
      real*8,optional :: n, h, nm, LAI,
     &     crown_dx, crown_dy,dbh, root_d, fracroot(:),clump,
     &     LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix,
     &     phenofactor_c, phenofactor_d, phenofactor, dphdt, CB_d
!     &     stressH2O, stressH2Ol(:)
      cop%pft = pft
      cop%n = n
      cop%nm = nm
      cop%LAI = LAI
!      write(777,*) __FILE__,__LINE__,cop%LAI
      cop%h = h
      cop%crown_dx =  crown_dx 
      cop%crown_dy =  crown_dy
      cop%dbh =  dbh 
      cop%root_d = root_d
      cop%fracroot(:) = fracroot(:)
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
      cop%Ci = Ci
      cop%GCANOPY =  GCANOPY 
      cop%GPP =  GPP 
      cop%NPP =  NPP 
      cop%R_auto = 0.0
      cop%R_root = 0.0 !PK 5/15/07
      cop%N_up =  N_up 
!      cop%C_litter = C_litter
!      cop%N_litter = N_litter
      cop%C_to_Nfix = C_to_Nfix
      cop%phenofactor_c = phenofactor_c
      cop%phenofactor_d = phenofactor_d
      cop%phenofactor = phenofactor
      cop%dphdt = dphdt
      cop%CB_d = CB_d
!      cop%stressH2O = stressH2O      !Calculated in biophysics
!      cop%stressH2Ol = stressH2Ol    !Calculated in biophysics
      end subroutine assign_cohort
      !*********************************************************************
!!! basically replaced with cohort_construct()
cddd      subroutine init_cohort_defaults(cop,pnum)
cddd!     @sum Initialize a cohort with default values
cddd      type(cohort),pointer :: cop
cddd      integer :: pnum
cddd
cddd      cop%pft = pnum
cddd      cop%n = 1                 !## Dummy ##!
cddd!     cop%pptr = pp
cddd!     call nullify(cop%taller) !Only one cohort
cddd!     call nullify(cop%shorter)
cddd!     call nullify(cop%csptaller)
cddd!     call nullify(cop%cspshorter)            
cddd
cddd      call zero_cohort(cop)
cddd      
cddd      end subroutine init_cohort_defaults

      subroutine zero_cohort(cop)
!@sum Zero all real variables in cohort record.      
      use growthallometry,only : init_rootdistr
      type(cohort),pointer :: cop

      cop%nm = 0.0
      cop%Ntot = 0.0

      !* Individual plant properties *!
      !* GEOMETRY *!
      cop%h = 0.0
      cop%crown_dx = 0.0
      cop%crown_dy = 0.0
      cop%dbh = 0.0
      cop%root_d = 0.0
      cop%LAI = 0.0
      cop%clump = 0.0
      cop%fracroot(:) = 0.0

      !* BIOMASS POOLS *!
      cop%LMA = 0.0
      cop%C_fol = 0.0
      cop%N_fol = 0.0
      cop%C_sw = 0.0
      cop%N_sw = 0.0
      cop%C_hw = 0.0
      cop%N_hw = 0.0
      cop%C_lab = 0.0
      cop%N_lab = 0.0
      cop%C_froot = 1.0         !Dummy
      cop%N_froot = 0.0
      cop%C_croot = 0.0
      cop%N_croot = 0.0
      
      !* FLUXES *!
      cop%Ci =  0.0127D0        !Initial value not zero.
      cop%GCANOPY = 0.0
      cop%GPP = 0.0
      cop%NPP = 0.0
      cop%R_auto = 0.0
      cop%R_root = 0.0 !PK 5/15/07
      cop%N_up = 0.0
!      cop%C_litter = 0.0
!      cop%N_litter = 0.0
      cop%C_to_Nfix = 0.0

      !* PHENOLOGY/GROWTH *!
      cop%phenofactor_c = 0.d0
      cop%phenofactor_d = 0.d0
      cop%phenofactor = 0.d0
      cop%dphdt = 0.d0
      cop%CB_d = 0.d0

      !* PHYSIOLOGICAL STATUS *!
      cop%stressH2O = 1.d0 !Default no stress.
      cop%stressH2Ol(:) = 1.d0 !Default no stress.
      cop%senescefrac = 0.d0
      end subroutine zero_cohort


      !*********************************************************************
      !*********************************************************************
       
      subroutine reorganize_cohorts(pp)
      type(patch),pointer :: pp

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not need in GISS replication test.
      
      end subroutine reorganize_cohorts
      !*********************************************************************


      subroutine cohort_construct(cop, parent_patch, pnum)
      !@sum create a cohort with default values. if optional values
      !@+ are provided - set them
      ! this function may eventually be combined with assign_cohort
      ! for better performance
      type(cohort),pointer :: cop
      integer, optional :: pnum
      type(patch), optional, target :: parent_patch

      ! allocate memory
      allocate( cop )
      allocate( cop%fracroot(N_DEPTH) )
      allocate( cop%stressH2Ol(N_DEPTH) )

      ! set pointers if any
      nullify(cop%cellptr )
      nullify(cop%pptr )
      nullify(cop%taller )
      nullify(cop%shorter )
      nullify(cop%csptaller )
      nullify(cop%cspshorter )
      if ( present(parent_patch) ) then
        cop%pptr => parent_patch
        if ( associated( cop%pptr%cellptr ) )
     &       cop%cellptr => cop%pptr%cellptr
      endif

      ! set variables
      cop%pft = -1              ! = -1 if pft not set
      if ( present(pnum) ) cop%pft = pnum
      cop%n = 0.0

      call zero_cohort(cop)

      end subroutine cohort_construct
      !*********************************************************************


      subroutine cohort_destruct(cop)
      !@sum deallocate memory used by cohort
      type(cohort),pointer :: cop

      ! we may want ot collapse hole between "taller" and "shorter"
      ! here if this functionality is needed

      ! deallocate all memory
      deallocate( cop%stressH2Ol )
      deallocate( cop%fracroot )
      deallocate( cop )
      nullify( cop )

      end subroutine cohort_destruct


      subroutine cohort_print(iu, cop, prefix)
      integer, intent(in) :: iu
      type(cohort), intent(in) :: cop
      character*(*), optional, intent(in) :: prefix
      !---
      integer n

      write(iu, '(a,a," = ",i7)') prefix,"pft ",cop%pft
      write(iu, '(a,a," = ",f10.7)') prefix,"n   ",cop%n
      write(iu, '(a,a," = ",f10.7)') prefix,"nm  ",cop%nm
! Ntot doesn't seem terribly important ( always LAI*nm? )
      !write(iu, '(a,a," = ",f10.7)') prefix,"Ntot",cop%Ntot
      write(iu, '(a,a," = ",f10.7)') prefix,"LAI ",cop%LAI
      write(iu, '(a,a," = ",f10.7)') prefix,"h   ",cop%h
      write(iu, '(a,a," = ",f10.7)') prefix,"dbh ",cop%dbh
      write(iu, '(a,"fracroot = " )') prefix
      do n=1,N_DEPTH
        write(iu, '(a,"      ",f10.7)') prefix,cop%fracroot(n)
      enddo
      write(iu, '(a,a," = ",f10.7)') prefix,"Gcan",cop%gcanopy
      write(iu, '(a,a," = ",f10.7)') prefix,"GPP ",cop%GPP


      end subroutine cohort_print
      
      
      !*********************************************************************
      subroutine calc_CASArootfrac(cop,fracrootCASA)  !PK 11/06
      !maps fracroot(N_DEPTH) to fracrootCASA(N_CASA_LAYERS)
      !needs to be customized based on thicknesses of CASA layers and GCM layers 
      type(cohort),intent(in) :: cop
      real*8,intent(out) :: fracrootCASA(N_CASA_LAYERS)
      
      if (N_CASA_LAYERS == 1) then
         fracrootCASA = 1.d0  !if there is no explicit depth structure
      else
      !***scheme for N_CASA_LAYERS=2 (layers: 0-30, 30-100 cm)*** 
         fracrootCASA(1) = cop%fracroot(1) + cop%fracroot(2)  !CASA layer 1 --> GISS GCM layers 1,2
         fracrootCASA(2) = cop%fracroot(3) + cop%fracroot(4)  !CASA layer 2 --> GISS layers 3,4
     &                + cop%fracroot(5)                    !need to add 5th GISS layer (mainly for trees) -PK 6/26/07
      end if
                               
      end subroutine calc_CASArootfrac
      !*********************************************************************


      end module cohorts
