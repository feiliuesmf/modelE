      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types
      use ent_const
 
      implicit none


      public phenology_update, litter

      contains
      !*********************************************************************
      subroutine phenology_update_cell(dtsec,time,ecp)
!@sum Updates phenology for all patches in an entcell.
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: ecp
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
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: tt  !Greenwich Mean Time
      !integer,intent(in) :: latj !j index for latitude of entcell
      type(patch),pointer :: pp

      !---------------------------------------------------------------
      !* Prognostic phenology is calculated here.

      end subroutine phenology_update


      !*********************************************************************
      subroutine litter(dtsec, pp)
      !* Determine litter from live carbon pools and update Tpool.
      !* After CASA. - NYK 7/27/06

      real*8 :: dtsec           !dt in seconds
      !type(timestruct) :: tt    !Greenwich Mean Time
      type(patch),pointer :: pp
      !--Local-----------------
      type(cohort),pointer :: cop
      real*8 :: Closs(PTRACE,NPOOLS) !Litter per cohort.
      real*8 :: Clossacc(PTRACE,NPOOLS) !Litter accumulator.
      integer :: pft

      Closs(:,:) = 0.d0
      Clossacc(:,:) = 0.d0

      !* Calculate fresh litter from each cohort *!
      cop => pp%tallest  !changed to => (!) -PK 7/11/06
      do while(ASSOCIATED(cop)) 
        !*## NOTE:  betad should eventually be defined at cohort level!!##
        pft = cop%pft

        !* NLIVE POOLS *! !* betad/SECPY to replace stressCD is temporary HACK ##
        Closs(CARBON,LEAF) = 
     &       cop%C_fol * cop%n * (annK(pft,LEAF)+pp%betad/SECPY)*dtsec !* x tune factor
        Closs(CARBON,FROOT) = 
     &       cop%C_froot * cop%n * annK(pft,FROOT)*dtsec
        Closs(CARBON,WOOD) = 
     &       cop%C_hw * cop%n
     &       *(1.d0-exp(-annK(pft,WOOD)*dtsec)) !* expr is kdt; x tune factor

!        write(98,*) 'In litter: ',dtsec
!        write(98,*) cop%pft, cop%C_fol, cop%n, annK(pft,LEAF),pp%betad
!        write(98,*) cop%pft, cop%C_froot, cop%n, annK(pft,FROOT)
!        write(98,*) cop%pft, cop%C_hw, cop%n, annK(pft,WOOD)
!        write(98,*) 'solubfract(pft)', solubfract(pft)
!        write(98,*) 'Closs(CARBON,LEAF)',Closs(CARBON,LEAF)
!        write(98,*) 'Closs(CARBON,FROOT)',Closs(CARBON,FROOT)
        
        Clossacc(CARBON,LEAF) = Clossacc(CARBON,LEAF)
     &       + Closs(CARBON,LEAF)
        Clossacc(CARBON,FROOT) = Clossacc(CARBON,FROOT) 
     &       + Closs(CARBON,FROOT)
        Clossacc(CARBON,WOOD) = Clossacc(CARBON,WOOD) 
     &       + Closs(CARBON,WOOD)

        !* NDEAD POOLS *!
        Clossacc(CARBON,SURFMET) = Clossacc(CARBON,SURFMET) 
     &       + Closs(CARBON,LEAF) * solubfract(pft)
        Clossacc(CARBON,SOILMET) = Clossacc(CARBON,SOILMET) 
     &       + Closs(CARBON,FROOT) * solubfract(pft)
        Clossacc(CARBON,SURFSTR) = Clossacc(CARBON,SURFSTR)
     &       + Closs(CARBON,LEAF) * (1-solubfract(pft))
        Clossacc(CARBON,SOILSTR) = Clossacc(CARBON,SOILSTR) 
     &       + Closs(CARBON,FROOT) * (1-solubfract(pft))
        Clossacc(CARBON,CWD) = Clossacc(CARBON,CWD) 
     &       + Closs(CARBON,WOOD)
     
        cop => cop%shorter  !added -PK 7/12/06
      end do

      !* NDEAD POOLS *!
      pp%Tpool(CARBON,SURFMET) = pp%Tpool(CARBON,SURFMET) 
     &     + Clossacc(CARBON,SURFMET)
      pp%Tpool(CARBON,SOILMET) = pp%Tpool(CARBON,SOILMET) 
     &     + Clossacc(CARBON,SOILMET)
      pp%Tpool(CARBON,SURFSTR) = pp%Tpool(CARBON,SURFSTR)
     &     + Clossacc(CARBON,SURFSTR)
      pp%Tpool(CARBON,SOILSTR) = pp%Tpool(CARBON,SOILSTR) 
     &     + Clossacc(CARBON,SOILSTR)
      pp%Tpool(CARBON,CWD) = pp%Tpool(CARBON,CWD) 
     &     + Clossacc(CARBON,CWD)
      
      end subroutine litter
      !*********************************************************************

      end module phenology
