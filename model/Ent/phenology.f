      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types
      use ent_const
      use ent_pfts
 
      implicit none


      public phenology_update, litter
      public phenology_stats

      contains
      !*********************************************************************
!KIM      subroutine phenology_update_cell(dtsec,time,ecp)
!KIM!@sum Updates phenology for all patches in an entcell.
!KIM      real*8 :: dtsec           !dt in seconds
!KIM      type(timestruct) :: time  !Greenwich Mean Time
!KIM      type(entcelltype) :: ecp
!KIM      !------local--------
!KIM      type(patch),pointer :: pp
!KIM
!KIM      pp = ecp%youngest
!KIM      do while (ASSOCIATED(pp))
!KIM        call phenology_update(dtsec, time, pp)
!KIM      end do
!KIM      end subroutine phenology_update_cell
      !*********************************************************************
      subroutine phenology_update(pp)
!@sum Update phenology for a patch.

      use ent_prescr_veg !Temoprary

      type(patch),pointer :: pp 
      type(cohort), pointer :: cop
      integer :: pft
      real*8, parameter :: laimax(N_PFT) = !temp. from ent_ENTveg.f
     &     alamax((COVEROFFSET+1):(COVEROFFSET+N_PFT))
!     $     (/ 8.0d0, 8.0d0, 10.0d0, 10.0d0, 6.0d0 ,6.0d0, 4.0d0
!     &     ,2.5d0, 2.5d0, 2.0d0, 2.0d0, 2.0d0, 4.5d0/)

      !---------------------------------------------------------------
      !* Prognostic phenology is calculated here.
      
      cop => pp%tallest

      do while(ASSOCIATED(cop))
        pft = cop%pft   

        !Temp. until the allocation scheme available:
        !(Leaf Area Index) 
        ! = (Maximum Leaf Area Index) * (Phenology Factor)
        !No carbon pool update for now.

        cop%LAI= laimax(pft) * cop%phenofactor
        
        cop => cop%shorter 
 
      end do

!      print *, 'phenology_update ok'

      end subroutine phenology_update

      !*********************************************************************
      subroutine phenology_stats(dtsec, pp, dailyupdate)
!@sum Update statstics for phneology_update    

      real*8 :: dtsec           !dt in seconds
      type(patch) :: pp
      logical, intent(in) :: dailyupdate

      !----Local----------
      type(cohort), pointer :: cop

      real*8 :: airtemp_par, gdd_par1, gdd_par2, gdd_par3 
      real*8 :: gdd_threshold, ddfacu
      real*8 :: ld_threshold, tsoil_threshold1, tsoil_threshold2
      real*8 :: wat_threshold
      real*8 :: sandfrac, clayfrac
      real*8 :: watsat, smpsat, bch, watdry 

      real*8 :: airtemp
      real*8 :: soilmoist    
      real*8 :: soiltemp     
      real*8 :: coszen

      integer :: pft
      integer :: phenotype

      real*8 :: zweight

      real*8 :: soiltemp_10d      
      real*8 :: airtemp_10d
      real*8 :: soilmoist_10d

      real*8 :: gdd
      real*8 :: ncd
      real*8 :: ld
      real*8 :: wat
  
      real*8 :: phenofactor_c
      real*8 :: phenofactor_d
      real*8 :: phenofactor

      real*8 :: phenofactor_c_old

      airtemp_par = 5.d0
      gdd_par1 = -68.d0   
      gdd_par2 = 638.d0
      gdd_par3 = -0.01d0
      ld_threshold = 655.d0
      tsoil_threshold1 = 11.15d0
      tsoil_threshold2 = 2.d0
      ddfacu = 1.d0/15.d0      
      wat_threshold = 0.5d0

      sandfrac = pp%cellptr%soil_texture(1)
      clayfrac = pp%cellptr%soil_texture(3)
      
      airtemp =  pp%cellptr%TairC
!      soilmoist = pp%Soilmoist  !vary with patches, but currently doesn't work at the patch level
      soilmoist = pp%cellptr%Soilmoist
      soiltemp = pp%cellptr%Soiltemp  
      coszen = pp%cellptr%CosZen
    
  
!*****************************************     
!* Update Climate/Environmental Conditions
!*****************************************

      zweight=exp(-1.d0/(10.d0*86400.d0/dtsec)) !weighting factor for 10-day running average

      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      soilmoist_10d = pp%cellptr%soilmoist_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd
      ld =  pp%cellptr%ld
      wat = pp%cellptr%wat

      cop => pp%tallest
      do while(ASSOCIATED(cop))

         !10-day running average of Soil Temperature
         soiltemp_10d=zweight*soiltemp_10d+(1.0-zweight)*soiltemp
      
         !10-day running average of Air Temperature
         airtemp_10d=zweight*airtemp_10d+(1.0-zweight)*airtemp
     
         !10-day running average of Soil Moisture
         soilmoist_10d=zweight*soilmoist_10d+(1.0-zweight)*soilmoist

         !GDD & NCD - Once a day
         if (dailyupdate) then
            !Growing degree days
            if (airtemp_10d .lt. airtemp_par) then
               gdd = 0.d0
            else
               gdd = gdd + airtemp_10d 
            end if

            !Number of chilling days
            if (airtemp_10d .lt. airtemp_par) then
               ncd = ncd +  1.d0
            end if

         end if
  
         !Photoperiod (Day length) in minute
         if ( coszen .gt. 0.d0 ) then
            ld = ld + dtsec/60.d0
         end if

         !Soil Saturation - following the calc. in soilbgc.f
         watsat = 0.489d0 - 0.00126d0*sandfrac 
         smpsat = -10.d0 * ( 10.d0**(1.88d0-0.0131d0*sandfrac) )
         bch = 2.91d0 + 0.159d0*clayfrac
         watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
         wat = min( max(soilmoist_10d-watdry,0.d0) 
     &         /(watsat-watdry), 1.d0)      

!**************************    
!* Update Phenology Factor
!**************************

         if (dailyupdate) then
            phenofactor_c=cop%phenofactor_c
            phenofactor_d=cop%phenofactor_d
            phenofactor=cop%phenofactor
            pft=cop%pft
            phenotype=pfpar(pft)%phenotype
            phenofactor_c_old=cop%phenofactor_c

            ! Cold Deciduousness

            if ((phenotype .eq. 2) .OR. (phenotype .eq. 4)) then 
               gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
               if ( ((ld .le. ld_threshold) .AND. 
     &           (soiltemp_10d .le. tsoil_threshold1)) 
     &          .OR. (soiltemp_10d .le. tsoil_threshold2) ) then
                  phenofactor_c = max (0.0d0, phenofactor_c - ddfacu)
               else if (gdd .gt. gdd_threshold .AND.
     &                  gdd .le. 1000.0d0 ) then

                  phenofactor_c = min (1.0d0, phenofactor_c + ddfacu)
!option for 300 growing days, instead of 15 days
!                  phenofactor_c = min (1.0d0, 
!     &                            (gdd-gdd_threshold)/300.0d0)
               end if
            else 
               phenofactor_c = 1.0d0
            end if

            ! Drought Deciduousness
            if ((phenotype .eq. 3) .OR. (phenotype .eq. 4)) then
               phenofactor_d =  min (1.0d0, wat/wat_threshold)                    
            else
               phenofactor_d = 1.0d0
            end if

            phenofactor = phenofactor_c * phenofactor_d

            cop%phenofactor_c=phenofactor_c
            cop%phenofactor_d=phenofactor_d
            cop%phenofactor=phenofactor
            !* cop%senescefrac = INSERT CALCULATION IF <> phenofactor

            ld =  0.d0 !zero-out every day
            
            if ((phenofactor_c_old .gt. phenofactor_c) .AND.
     &         phenofactor_c .eq. 0.d0 ) then 
               ncd = 0.d0 !zero-out every summer
            end if

         end if   

         cop => cop%shorter 
 
      end do        

      pp%cellptr%soiltemp_10d = soiltemp_10d
      pp%cellptr%airtemp_10d = airtemp_10d
      pp%cellptr%soilmoist_10d = soilmoist_10d
      pp%cellptr%gdd = gdd
      pp%cellptr%ncd = ncd
      pp%cellptr%wat = wat
      pp%cellptr%ld =  ld

      end subroutine phenology_stats

      !*********************************************************************
      subroutine litter( pp)
      !* Determine litter from live carbon pools and update Tpool.
      !* After CASA, but called at daily time step. - NYK 7/27/06

      !real*8 :: dtsec           !dt in seconds
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

        !* NLIVE POOLS *! !* betad to replace stressCD is temporary HACK ##
!        Closs(CARBON,LEAF) =
!     &       cop%C_fol * cop%n * (annK(pft,LEAF)+pp%betad/SECPY)*SDAY !* x tune factor
!        Closs(CARBON,FROOT) =
!     &       cop%C_froot * cop%n * annK(pft,FROOT)*SDAY
!        Closs(CARBON,WOOD) =
!     &       cop%C_hw * cop%n
!     &       *(1.d0-exp(-annK(pft,WOOD)*SDAY)) !* expr is kdt; x tune factor


!*   PUSHKER - COMMENT OUT THE ABOVE LINES AND UNCOMMENT THESE LINES FOR NEW
!*             LITTER FLUX - NANCY
!*   Later will replace above with senescefrac factors, which can be calculated by
!*   either prescribed or prognostic phenology: ****** NYK!
        Closs(CARBON,LEAF) = 
     &       cop%C_fol * cop%n *
     &       (annK(pft,LEAF)*SDAY +cop%senescefrac) !* x tune factor
        Closs(CARBON,FROOT) = 
     &       cop%C_froot * cop%n * 
     &       (annK(pft,FROOT)*SDAY + cop%senescefrac) !* x tune factor
        Closs(CARBON,WOOD) = 
     &       (cop%C_hw + cop%C_croot) * cop%n
     &       *(1.d0-exp(-annK(pft,WOOD)*SDAY)) !* expr is kdt; x tune factor


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
