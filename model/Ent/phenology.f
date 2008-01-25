      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types
      use ent_const
      use ent_pfts
 
      implicit none

      public phenology_update !may change the name into veg_update 
      public phenology_stats
      public litter   
      private growth_cpools
      private phenology_cpools

      !l_fract: fraction of leaves retained after leaf fall (unitless)
      real*8, parameter :: l_fract = 0.50d0 
      !q: ratio of root to leaf biomass (unitless)
      real*8, parameter :: q=1.0d0 
      !iqsw: sapwood biomass per (leaf area x wood height) (kgC/m2/m)
      !3900.0: leaf area per sapwood area (m2/m2) 
      !1000.0: sapwood density (kg/m3)
      !2.0:  biomass per carbon (kg/kgC)
      !(qsw)=(iqsw*sla) (1/m) & (qsw*h): ratio of sapwood to leaf biomass (unitless)
      real*8, parameter :: iqsw=1000.0d0/3900.0d0/2.0d0
      !hw_fract: ratio of above ground stem to total stem (stem plus structural roots)
      real*8, parameter :: hw_fract = 0.70d0 
      !C2B: ratio of biomass to carbon (kg-Biomass/kg-Carbon)
      real*8, parameter :: C2B = 2.0d0

      contains

      !*********************************************************************
      subroutine phenology_stats(dtsec, pp, dailyupdate,time)
!@sum Update statstics for phneology_update    
      real*8,intent(in) :: dtsec           !dt in seconds
      type(patch) :: pp
      logical, intent(in) :: dailyupdate
      real*8,intent(in) :: time
      !--Local-----
      type(cohort), pointer :: cop
 
      !temperature constrain for cold-deciduous PFTs (Botta et al. 1997)
      !airtemp_par: base temperature to estimate growing degree day (gdd) (degree C)
      real*8, parameter :: airtemp_par = 5.d0 
      ! gdd_par1/2/3: paramters to estimate the threshold for gdd
      ! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
      real*8, parameter :: gdd_par1 = -68.d0 
      real*8, parameter :: gdd_par2 = 638.d0
      real*8, parameter :: gdd_par3 = -0.01d0 
      !ld_threshold (minute): light length constraint for cold-deciduous woody PFTs (White et al. 1997)
      real*8, parameter :: ld_threshold = 655.d0
      !tsoil_threshold1, tsoil_threshold2 : soil temperature constraint for cold-deciduous woody PFTs (White et al. 1997)
      real*8, parameter :: tsoil_threshold1 = 11.15d0
      real*8, parameter :: tsoil_threshold2 = 2.d0
      !ddfacu: the rate of leaf fall (1/day) (IBIS)
      real*8, parameter :: ddfacu = 1.d0/15.d0
      !wat_threshold: water threshold (unitless)       
      real*8, parameter :: wat_threshold = 0.5d0

      real*8 :: gdd_threshold

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
      

      sandfrac = pp%cellptr%soil_texture(1)
      clayfrac = pp%cellptr%soil_texture(3)
      
      airtemp =  pp%cellptr%TairC
!      soilmoist = pp%Soilmoist  !vary with patches, but currently doesn't work at the patch level
      soilmoist = pp%cellptr%Soilmoist(1)  !assign these only to top 30 cm for now -PK 3/16/07
      soiltemp = pp%cellptr%Soiltemp(1)
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
         
         !daily carbon balance for growth
         cop%CB_d =  cop%CB_d + cop%NPP*dtsec/cop%n

         !10-day running average of Soil Temperature
         soiltemp_10d=zweight*soiltemp_10d+(1.0d0-zweight)*soiltemp
      
         !10-day running average of Air Temperature
         airtemp_10d=zweight*airtemp_10d+(1.0d0-zweight)*airtemp
   
         !10-day running average of Soil Moisture
         soilmoist_10d=zweight*soilmoist_10d+(1.0d0-zweight)*soilmoist

         !GDD & NCD - Once a day
         if (dailyupdate) then
            !Growing degree days
            if (airtemp_10d .lt. airtemp_par) then
               gdd = 0.d0
            else
               gdd = gdd + ( airtemp_10d - airtemp_par )
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
            cop%dphdt = phenofactor - cop%phenofactor
            cop%phenofactor_c=phenofactor_c
            cop%phenofactor_d=phenofactor_d
            cop%phenofactor=phenofactor

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
      subroutine phenology_update(dtsec,pp,dailyupdate,monthlyupdate)
!@sum Update the vegetation state and carbon pools.
!@sum i.e., LAI, senescefrac, DBH, height &
!@sum carbon pools of foliage, sapwood, fineroot, hardwood, coarseroot

      use ent_prescr_veg
  
      type(patch),pointer :: pp 
      type(cohort), pointer :: cop
      real*8 :: dtsec
      logical :: dailyupdate
      logical :: monthlyupdate
      integer :: pft
      ! phenofactor phenological elongation factor [0,1] (unitless)
      real*8 :: phenofactor
      ! dphdt difference between today's (unpdated) and yesterday's phenofact (unitless)
      real*8 :: dphdt
      real*8 :: C_fol_old
      real*8 :: C_fol
      real*8 :: C_lab
      ! Cactive active carbon pool, including foliage, sapwood and fine root (gC/pool/individual)
      real*8 :: Cactive
      ! Cactive_max maximum active carbon pool allowed by the allometric constraint 
      real*8 :: Cactive_max
      real*8 :: C_hw
      real*8 :: C_croot
      real*8 :: C_froot
      ! Cdead dead carbon pool, including hardwood and coarse root (gC/pool/individual)
      real*8 :: Cdead
      real*8 :: dCdead
      ! dC_lab_growth change in the labile carbon due to carbon allocation to active and dead pool
      real*8 :: dC_lab_growth
      ! dC_lab_growth change in the labile carbon due to carbon allocation to active pool
      real*8 :: dC_lab_growth_active
      ! dC_lab_litter change in the labile carbon due to litter
      real*8 :: dC_lab_litter
      ! dC_fol_growth change in the foliage carbon due to carbon allocation to active pool
      real*8 :: dC_fol_growth
      ! dC_fol_litter change in the foliage carbon due to litter
      real*8 :: dC_fol_litter
      ! dC_fol sum of dC_fol_growth+dC_fol_litter
      real*8 :: dC_fol
      real*8 :: Closs(PTRACE,NPOOLS) !Litter per cohort.
      real*8 :: Clossacc(PTRACE,NPOOLS) !Litter accumulator.
      ! qsw 
      real*8 :: qsw
      ! dbh diameter at the breast height (cm)
      real*8 :: dbh
      ! h plant height (m)
      real*8 :: h
      ! plant population (#-individual/m2-ground)
      real*8 :: nplant
      real*8 :: ialloc
      real*8 :: senescefrac
      ! CB_d daily carbon balance (gC/individual)
      real*8 :: CB_d
      real*8 :: dC_lab
      real*8 laipatch

      laipatch = 0.d0
      cop => pp%tallest

      do while(ASSOCIATED(cop))
    
         if (dailyupdate) then 
            pft = cop%pft
            phenofactor = cop%phenofactor
            dphdt = cop%dphdt
            dbh = cop%dbh
            h = cop%h
            nplant = cop%n

            if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops
               qsw=0.0d0  !no allocation to the wood
            else
               qsw = pfpar(pft)%sla*iqsw
            end if
            ialloc = 1.0d0+q+h*qsw 
            Cactive = cop%C_froot * ialloc / q 
            if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops    
               Cactive_max = height2Cfol(pft,5.0d0,nplant)*ialloc 
               !no allometric constraints in the carbon allocation, 
               !then 5.0d0 is arbitrary number (must be tall enough, then 5m)
            else 
               Cactive_max = dbh2Cfol(pft,dbh)* ialloc
            end if
c$$$            C_lab = cop%C_lab*1000.0 
c$$$     &           - Cactive * (1.0-(phenofactor-dphdt)) / ialloc
            C_lab = cop%C_lab*1000.d0 
            C_fol_old = cop%C_fol
 
            C_froot = cop%C_froot
            C_hw = cop%C_hw
            C_croot = cop%C_croot
            Cdead = C_hw + C_croot   
            
            CB_d = cop%CB_d*1000.d0
           
            !allocate the labile carbon to the active/dead pools or 
            !store the carbon for the reproduction use
            call growth_cpools (
     i           pft,
     i           dbh,
     i           phenofactor, 
     i           dphdt,
     i           Cactive,
     i           Cactive_max,
     i           C_lab,
     i           CB_d,
     i           Cdead, 
     i           monthlyupdate, 
     o           dC_lab_growth, 
     o           dC_lab_growth_active, 
     o           dCdead)           

            !change in carbon allocation due to the phenological constraints
            call phenology_cpools(
     i           phenofactor, 
     i           dphdt, 
     i           ialloc,
     i           Cactive,
     i           C_fol_old,
     i           dC_lab_growth_active,
     o           dC_lab_litter, 
     o           dC_fol_litter, 
     o           dC_fol_growth,
     o           senescefrac)   
           
            !update the active pool 
            dC_fol = dC_fol_growth + dC_fol_litter !dC_fol is dummy.
            dC_lab = dC_lab_growth + dC_lab_litter
            Cactive = Cactive - (dC_lab_growth_active + dC_lab_litter) 
            cop%C_fol = Cactive 
     &            * phenofactor * (1.0d0/ialloc)
            cop%C_sw = Cactive *(h*qsw/ialloc)
            cop%C_froot = Cactive *(q/ialloc)  
            cop%C_lab = 0.001d0*( C_lab + dC_lab )  !g/individual -> kg/ind.
   
            !update the dead pool
            Cdead= Cdead+dCdead
            cop%C_hw = Cdead * hw_fract
            cop%C_croot = Cdead * (1-hw_fract)

            !update the plant size 
            if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops   
               cop%dbh = 0.0d0
               cop%h = Cfol2height(pft,cop%C_fol,nplant)
            else        
               cop%dbh = Cdead2dbh(pft,Cdead)
               cop%h = dbh2height(pft,cop%dbh)
            end if
 
            !update the senescefrac & LAI
            cop%senescefrac = senescefrac
            cop%LAI=cop%C_fol/1000.0d0*pfpar(pft)%sla*cop%n
            if (cop%LAI .lt. EPS) cop%LAI=EPS
         
            cop%Ntot = cop%nm * cop%LAI 
            laipatch = laipatch + cop%lai 

            !zero-out the daily accumulated carbon 
            cop%CB_d = 0.d0    

         end if
         cop => cop%shorter 
      end do

      pp%LAI = laipatch

      end subroutine phenology_update
      !*********************************************************************
      subroutine growth_cpools(
     i           pft,
     i           dbh,
     i           phenofactor, 
     i           dphdt, 
     i           Cactive,
     i           Cactive_max,
     i           C_lab,
     i           CB_d,
     i           Cdead, 
     i           monthlyupdate, 
     o           dC_lab_growth,  
     o           dC_lab_growth_active,
     o           dCdead)

      use ent_prescr_veg

      real*8 :: dtsec        

      type(cohort),pointer :: cop
      integer, intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: phenofactor
      real*8, intent(in) :: dphdt
      real*8, intent(in) :: Cactive
      real*8, intent(in) :: Cactive_max
      real*8, intent(in) :: C_lab
      real*8, intent(in) :: CB_d
      real*8, intent(in) :: Cdead
      logical, intent(in) :: monthlyupdate
      
      real*8 :: dC_lab_growth
      real*8 :: dC_lab_growth_active
      real*8 :: dCdead
      
      real*8 :: phehofactor_old
      real*8 :: Cactive_pot
      real*8 :: Cactive_inc

      !r_fract: fraction of excess c going to seed reproduction 
      real*8 :: r_fract 
      !c_fract: fraction of excess c going to clonal reproduction
      real*8 :: c_fract 
      !c_fract: fraction of excess c going to structural growth
      real*8 :: gr_fract  
      
      dC_lab_growth = 0.d0
      
      !* fast growth - Once a day   
      if ((phenofactor .eq. 0.d0) .AND. (CB_d .gt. 0.d0)) then
         dC_lab_growth = 0.d0
      else
         if (phenofactor .gt. 0.d0) then
            Cactive_pot = Cactive + CB_d
            !compare the max numer (max. allowed pool size according to the DBH)
            !with the potential number (current size + daily accumulated carbon) 
            if (Cactive_pot .gt. Cactive_max) then
               if (Cactive .lt. Cactive_max) then
                  dC_lab_growth = Cactive -Cactive_max
               else
                  dC_lab_growth = Cactive_max-Cactive
               end if
             else 
               dC_lab_growth = -CB_d
             end if
         else
            dC_lab_growth = -CB_d !negative CB_d    
         end if
      end if
   
      !* fast growth - Once a day only during growing season
      if ((dphdt .ge. 0.d0) .AND. (phenofactor .gt. 0.d0)) then
         if (Cactive_max .gt. (Cactive-dC_lab_growth)) then
            dC_lab_growth = dC_lab_growth -
     &                      min((Cactive_max - (Cactive-dC_lab_growth))
     &                      , (C_lab+dC_lab_growth))
         end if
      end if

      dC_lab_growth_active = dC_lab_growth 
                                
      !* slow growth - Once a month   
      if (monthlyupdate) then
         if (C_lab+dC_lab_growth .gt. 0.d0) then
            if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops(?) 
               r_fract = 0.3
               c_fract = 0.7
               gr_fract = 0.0 !no structural pools
            else   
               r_fract = 0.3
               c_fract = 0.0
               if (dbh .le. 0.99*maxdbh(pft))then
                  gr_fract = 1.0 - r_fract - c_fract
               else
                  gr_fract = 0.0 
               end if
            end if

            dCdead = gr_fract * (C_lab + dC_lab_growth)
            dC_lab_growth = dc_lab_growth-
     &             (r_fract+c_fract+gr_fract)*(C_lab+dC_lab_growth) 
                   !(r_fract+c_fract)*C_lab to be saved for reproduction!
           
         end if
      end if

      end subroutine growth_cpools
      !*********************************************************************
      subroutine phenology_cpools(
     i           phenofactor, 
     i           dphdt,
     i           ialloc,
     i           Cactive,
     i           C_fol_old,
     i           dC_lab_growth_active, 
     o           dC_lab_litter, 
     o           dC_fol_litter, 
     o           dC_fol_growth,
     o           senescefrac)      
 
      real*8 :: phenofactor
      real*8 :: dphdt
      real*8 :: ialloc
      real*8 :: Cactive
      real*8 :: dC_lab_growth_active
      real*8 :: dC_lab_litter
      real*8 :: dC_fol_litter 
      real*8 :: dC_fol_growth
      real*8 :: C_fol
      real*8 :: C_fol_old
      real*8 :: dC_fol
      real*8 :: senescefrac

      C_fol = phenofactor * (Cactive - dC_lab_growth_active) / ialloc  !foliage update from drought/cold stress
      dC_fol = C_fol - C_fol_old
      
      dC_fol_litter = 0.0
      dC_fol_growth = 0.0
      dC_lab_litter = 0.0
      senescefrac = 0.0
      if (dphdt .lt. 0.0 .AND. phenofactor .lt. 1.0) then
         if (dC_fol .lt. 0.0 ) then 
           dC_fol_litter = dC_fol
           dC_lab_litter = - dC_fol_litter * l_fract
           if (C_fol .ne. 0.0 ) then 
              senescefrac = - dC_fol_litter * (1.0 - l_fract) /C_fol_old
           endif
         else
           dC_fol_growth = dC_fol
         endif
      else
         dC_fol_growth = dC_fol
      end if
      
      end subroutine phenology_cpools

      !*********************************************************************
      subroutine litter( pp)
      !* Determine litter from live carbon pools and update Tpool.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      
      use cohorts, only : calc_CASArootfrac 

      !real*8 :: dtsec           !dt in seconds
      !type(timestruct) :: tt    !Greenwich Mean Time
      type(patch),pointer :: pp
      !--Local-----------------
      type(cohort),pointer :: cop
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdttotal, adj !Total, adjustment factor if larger than C_lab.

      Closs(:,:,:) = 0.d0
      Clossacc(:,:,:) = 0.d0

      !* Calculate fresh litter from each cohort *!
      cop => pp%tallest  !changed to => (!) -PK 7/11/06
      do while(ASSOCIATED(cop)) 
        pft = cop%pft
        
      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(cop,fracrootCASA)
!      print *, 'from litter(pheno*.f): fracrootCASA(:) =', fracrootCASA !***test*** -PK 11/27/06  

       do i=1,N_CASA_LAYERS  !do this over all CASA layers -PK
        !* NLIVE POOLS *! 
        turnoverdtleaf = annK(pft,LEAF)*SDAY
        turnoverdtfroot = annK(pft,FROOT)*SDAY
        turnoverdtwood = 1.d0-exp(-annK(pft,WOOD)*SDAY) !Sapwood not hardwood

        !* UPDATE C_LAB: Turnover should draw down C_lab. *!
        ! Check that amount not too large 
        turnoverdttotal = turnoverdtleaf+turnoverdtfroot+turnoverdtwood
        !if (turnoverdttotal.lt.cop%C_lab) then
          adj = 1.d0  !No adjustment, enough C_lab
        !else
        !  adj = (cop%C_lab - EPS)/turnoverdttotal !Turnover can reduce C_lab only to EPS.
        !endif
        !cop%C_lab = cop%C_lab - adj*turnoverdttotal
        !* NEED TO PUT RETRANSLOCATION IN HERE, TOO *!

        !* Calculate litter *!
        ! Senescefrac factor can be calculated by either prescribed or prognostic phenology: ****** NYK!
        if (i.eq.1) then  !only top CASA layer has leaf and wood litter -PK   
         Closs(CARBON,LEAF,i) = 
     &       cop%C_fol * cop%n *
     &       (adj*turnoverdtleaf + cop%senescefrac) !* x tune factor
         Closs(CARBON,WOOD,i) = 
     &       (cop%C_hw + cop%C_croot) * cop%n
     &       *(adj*turnoverdtwood) !* expr is kdt; x tune factor
        else    
         Closs(CARBON,LEAF,i) = 0.d0 
         Closs(CARBON,WOOD,i) = 0.d0
        end if
         Closs(CARBON,FROOT,i) =  !both layers have root litter -PK 
     &       fracrootCASA(i)*cop%C_froot * cop%n * 
     &       (adj*turnoverdtfroot + cop%senescefrac) !* x tune factor

!        write(98,*) 'In litter: ',dtsec
!        write(98,*) cop%pft, cop%C_fol, cop%n, annK(pft,LEAF),pp%betad
!        write(98,*) cop%pft, cop%C_froot, cop%n, annK(pft,FROOT)
!        write(98,*) cop%pft, cop%C_hw, cop%n, annK(pft,WOOD)
!        write(98,*) 'solubfract(pft)', solubfract(pft)
!        write(98,*) 'Closs(CARBON,LEAF)',Closs(CARBON,LEAF)
!        write(98,*) 'Closs(CARBON,FROOT)',Closs(CARBON,FROOT)
         Clossacc(CARBON,LEAF,i) = Clossacc(CARBON,LEAF,i)
     &        + Closs(CARBON,LEAF,i)
         Clossacc(CARBON,FROOT,i) = Clossacc(CARBON,FROOT,i) 
     &        + Closs(CARBON,FROOT,i)
         Clossacc(CARBON,WOOD,i) = Clossacc(CARBON,WOOD,i) 
     &        + Closs(CARBON,WOOD,i)

        !* NDEAD POOLS *!
         Clossacc(CARBON,SURFMET,i) = Clossacc(CARBON,SURFMET,i) 
     &        + Closs(CARBON,LEAF,i) * solubfract(pft)
         Clossacc(CARBON,SOILMET,i) = Clossacc(CARBON,SOILMET,i) 
     &        + Closs(CARBON,FROOT,i) * solubfract(pft)
         Clossacc(CARBON,SURFSTR,i) = Clossacc(CARBON,SURFSTR,i)
     &        + Closs(CARBON,LEAF,i) * (1-solubfract(pft))
         Clossacc(CARBON,SOILSTR,i) = Clossacc(CARBON,SOILSTR,i) 
     &        + Closs(CARBON,FROOT,i) * (1-solubfract(pft))
         Clossacc(CARBON,CWD,i) = Clossacc(CARBON,CWD,i) 
     &        + Closs(CARBON,WOOD,i)
        end do  !loop through CASA layers-->cumul litter per pool per layer -PK
     
        cop => cop%shorter  !added -PK 7/12/06
      end do  !loop through cohorts

      !* NDEAD POOLS *!
       do i=1,N_CASA_LAYERS
        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
     &     + Clossacc(CARBON,SURFMET,i)
        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
     &     + Clossacc(CARBON,SOILMET,i)
        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
     &     + Clossacc(CARBON,SURFSTR,i)
        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
     &     + Clossacc(CARBON,SOILSTR,i)
        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
     &     + Clossacc(CARBON,CWD,i)
       end do   !loop through CASA layers-->total C per pool per layer -PK
!       print *, __FILE__,__LINE__,'pp%Tpool=',pp%Tpool(CARBON,:,:) !***test*** -PK 7/24/07  

      end subroutine litter
 
!*********************************************************************
      real*8 function dbh2Cfol(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8 :: maxdbh

      if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops  
         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      else !woody
         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      end if
      
      dbh2Cfol=1000.0d0*(1.0/C2B) *pfpar(pft)%b1Cf 
     &        * min(dbh, maxdbh)**pfpar(pft)%b2Cf

      end function dbh2Cfol
!*************************************************************************
      real*8 function maxdbh(pft)
      integer,intent(in) :: pft       

      if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops 
         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      else !woody   
         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      end if
 
      end function maxdbh
!*************************************************************************
      real*8 function dbh2Cdead(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dbh2Cdead = (1.0/C2B) * pfpar(pft)%b1Cd * 
     &             dbh**pfpar(pft)%b2Cd * 1000.0d0

      end function  dbh2Cdead
!*************************************************************************
      real*8 function dbh2height(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dbh2height = 1.3 + pfpar(pft)%b1Ht * 
     &             (1.0-exp(pfpar(pft)%b2Ht*dbh))

      end function dbh2height
!*************************************************************************
      real*8 function height2dbh(pft,h)
      integer,intent(in) :: pft
      real*8, intent(in) :: h

      height2dbh = log(1.0-(h-1.3)/pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht

      end function height2dbh
!*************************************************************************
      real*8 function Cdead2dbh(pft,Cdead)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

      Cdead2dbh = (Cdead/1000.0d0*C2B/pfpar(pft)%b1Cd)
     &            **(1.0d0/pfpar(pft)%b2Cd)

      end function Cdead2dbh
!*************************************************************************
      real*8 function height2Cfol(pft,height,nplant)
      integer,intent(in) :: pft
      real*8, intent(in) :: height
      real*8, intent(in) :: nplant
      real*8,parameter :: h1Cf = 5.0d0  !1.660d0
      real*8,parameter :: h2Cf = 1.20d0 !1.500d0

      height2Cfol=(1.0d0/C2B) *h1Cf 
     &        *((height*100.0d0)**h2Cf)/nplant

      end function height2Cfol
!*************************************************************************
      real*8 function Cfol2height(pft,Cfol,nplant)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol !gC/pool/plant
      real*8, intent(in) :: nplant !# of individual/m2-ground
      real*8,parameter :: h1Cf = 5.0d0  !1.660d0
      real*8,parameter :: h2Cf = 1.20d0 !1.500d0

      Cfol2height=exp(log(Cfol*C2B/h1Cf*nplant)/h2Cf)/100.0d0

      end function Cfol2height
!*************************************************************************

      end module phenology
