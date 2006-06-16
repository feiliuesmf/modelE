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
      !use ent_GISSveg, only : GISS_phenology
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: tt  !Greenwich Mean Time
      !integer,intent(in) :: latj !j index for latitude of entcell
      type(patch),pointer :: pp

      !---------------------------------------------------------------
      !* GISS VERSION:  JUST UPDATES LAI USING MATTHEWS PRESCRIPTION
      !* Prescribed phenology is passed in.
      !call GISS_phenology(tt%jday,latj, pp)
      !* Prognostic phenology is calculated here.

      end subroutine phenology_update


      !*********************************************************************
      subroutine litter(dtsec, pp)
      !* Update pp%Tpool
      !* This is imitating CASA casa_litterfall.F.

      real*8 :: dtsec           !dt in seconds
      !type(timestruct) :: tt    !Greenwich Mean Time
      type(patch),pointer :: pp
      !--Local-----------------
      type(cohort),pointer :: cop
      real*8 :: Closs(PTRACE,NPOOLS) !Litter per cohort.
      real*8 :: Clossacc(PTRACE,NPOOLS) !Litter accumulator.
      integer :: n,pft

      Closs(:,:) = 0.d0
      Clossacc(:,:) = 0.d0

      !* Calculate fresh litter from each cohort *!
      cop = pp%tallest
      do while(ASSOCIATED(cop)) 
        !*## NOTE:  betad should eventually be defined at cohort level!!##
        pft = cop%pft

        !* NLIVE POOLS *!
        Closs(CARBON,LEAF) = 
     &       pp%Tpool(CARBON,LEAF) * (annK(pft,LEAF)-pp%betad)*dtsec !* x tune factor
        Closs(CARBON,FROOT) = 
     &       pp%Tpool(CARBON,FROOT) * annK(pft,FROOT)*dtsec
        Closs(CARBON,WOOD) = 
     &       pp%Tpool(CARBON,WOOD) * 
     &       (1.d0-exp(-annK(pft,WOOD)*dtsec)) !* expr is kdt; x tune factor

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
      end do

      !* Patch summary of Tpool *!
      !* Don't let NLIVE Tpools go below zero. In CASA, bound is LeafMin *!
      do n = 1,NLIVE
        pp%Tpool(CARBON,n) = 
     &       max(0.d0,pp%Tpool(CARBON,n)-Clossacc(CARBON,n))
      end do
      !* NDEAD POOLS *!
      pp%Tpool(CARBON,SURFMET) = pp%Tpool(CARBON,SURFMET) 
     &     + Clossacc(CARBON,LEAF)
      pp%Tpool(CARBON,SOILMET) = pp%Tpool(CARBON,SOILMET) 
     &     + Clossacc(CARBON,FROOT)
      pp%Tpool(CARBON,SURFSTR) = pp%Tpool(CARBON,SURFSTR)
     &     + Clossacc(CARBON,LEAF)
      pp%Tpool(CARBON,SOILSTR) = pp%Tpool(CARBON,SOILSTR) 
     &     + Clossacc(CARBON,FROOT)
      pp%Tpool(CARBON,CWD) = pp%Tpool(CARBON,CWD) 
     &     + Clossacc(CARBON,WOOD)
      
      end subroutine litter
      !*********************************************************************

      end module phenology
