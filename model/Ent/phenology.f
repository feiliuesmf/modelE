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
      !* This is imitating CASA casa_littefall.F.

      real*8 :: dtsec           !dt in seconds
      !type(timestruct) :: tt    !Greenwich Mean Time
      type(patch),pointer :: pp
      !--Local-----------------
      type(cohort),pointer :: cop
      real*8 :: Closs(PTRACE,NPOOLS,N_PFT)
      integer :: n,p,pft

      Closs(:,:,:) = 0.d0

      !* Calculate fresh litter *!
      cop = pp%tallest
      do while(ASSOCIATED(cop)) 
        !*## NOTE:  betad should eventually be defined at cohort level!!##
        pft = cop%pft
        Closs(CARBON,LEAF,pft) = Closs(CARBON,LEAF,pft) + 
     &       pp%Tpool(CARBON,LEAF,pft) * (annK(pft,LEAF)-pp%betad)*dtsec !* x tune factor
        Closs(CARBON,FROOT,pft) = Closs(CARBON,FROOT,pft) + 
     &       pp%Tpool(CARBON,FROOT,pft) * annK(pft,FROOT)*dtsec
        Closs(CARBON,WOOD,pft) = Closs(CARBON,WOOD,pft) + 
     &       pp%Tpool(CARBON,WOOD,pft) * 
     &       (1.d0-exp(-annK(pft,WOOD)*dtsec)) !* expr is kdt; x tune factor
      end do

      !* Don't let Tpool go below zero. In CASA, bound is LeafMin *!
      do p = 1,N_PFT
        do n = 1,NLIVE
          pp%Tpool(CARBON,n,p) = 
     &        max(0.d0,pp%Tpool(CARBON,n,p)-Closs(CARBON,n,p))
        end do
      end do

      do n = 1,N_PFT
        pp%Tpool(CARBON,SURFMET,n) = pp%Tpool(CARBON,SURFMET,n) 
     &       + Closs(CARBON,LEAF,n) * solubfract(n)
        pp%Tpool(CARBON,SOILMET,n) = pp%Tpool(CARBON,SOILMET,n) 
     &       + Closs(CARBON,FROOT,n) * solubfract(n)
        pp%Tpool(CARBON,SURFSTR,n) = pp%Tpool(CARBON,SURFSTR,n)
     &       + Closs(CARBON,LEAF,n) * (1-solubfract(n))
        pp%Tpool(CARBON,SOILSTR,n) = pp%Tpool(CARBON,SOILSTR,n) 
     &       + Closs(CARBON,FROOT,n) * (1-solubfract(n))
        pp%Tpool(CARBON,CWD,n) = pp%Tpool(CARBON,CWD,n) 
     &       + Closs(CARBON,WOOD,n)
      end do

      end subroutine litter
      !*********************************************************************

      end module phenology
