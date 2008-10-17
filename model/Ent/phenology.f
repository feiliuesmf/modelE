      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence

      use ent_types
      use ent_const
      use ent_pfts
 
      implicit none
      public veg_init
      public clim_stats
      public pheno_update
      public veg_update !may change the name into veg_update 
      public litter   
      private growth_cpools_active
      private growth_cpools_structural
      private senesce_cpools
      private photosyn_acclim

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
      subroutine veg_init(pp)

      use ent_const
      use ent_pfts

      implicit none
      save

      type(patch) :: pp
      type(cohort), pointer :: cop 
      integer :: pft
      real*8 :: dbh
      real*8 :: h
      real*8 :: qsw
      real*8 :: ialloc

   
      pp%cellptr%soiltemp_10d = 0.0d0
      pp%cellptr%airtemp_10d = 0.0d0
      pp%cellptr%paw_10d = 0.50d0

      cop => pp%tallest   

      do while(ASSOCIATED(cop))
         pft=cop%pft
         h=cop%h
         if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops 
            dbh = 0.0d0  
            cop%n=cop%LAI/pfpar(pft)%sla/(height2Cfol(pft,h)/1000.0d0) 
            qsw = 0.d0
            cop%C_hw = 0.0d0
            cop%C_croot =  0.0d0
         else
            dbh=height2dbh(pft,h)
            cop%n=cop%LAI/pfpar(pft)%sla/(dbh2Cfol(pft,dbh)/1000.0d0)
            qsw = pfpar(pft)%sla*iqsw
            cop%C_hw = dbh2Cdead(pft,dbh) * hw_fract
            cop%C_croot = dbh2Cdead(pft,dbh) * (1.0d0 - hw_fract)
         end if
         cop%dbh=dbh
         cop%C_fol = 1000.0d0*cop%LAI/pfpar(pft)%sla/cop%n
         cop%C_froot =  q*1000.0d0*cop%LAI/pfpar(pft)%sla/cop%n
         ialloc = (1.0d0+q+h*qsw)
         cop%C_sw = cop%C_fol * h * qsw  
   
         cop => cop%shorter  

      end do
      end subroutine veg_init
      !*********************************************************************
      subroutine clim_stats(dtsec, pp, config,dailyupdate)
!@sum Calculate climate statistics such as 10 day running average   
      real*8,intent(in) :: dtsec           !dt in seconds
      type(patch) :: pp
      type(ent_config) :: config
      logical, intent(in) :: dailyupdate

      !--Local-----
      type(cohort), pointer :: cop
 
      !temperature constrain for cold-deciduous PFTs (Botta et al. 1997)
      !airtemp_par: base temperature to estimate growing degree day (gdd) (degree C)
      real*8, parameter :: airtemp_par = 5.d0 

      real*8 :: gdd_threshold

      real*8 :: sandfrac, clayfrac
      real*8 :: smpsat, bch
      real*8 :: watsat, watdry 
      real*8 :: paw

      real*8 :: airtemp
      real*8 :: wat   
      real*8 :: soiltemp     
      real*8 :: coszen

      real*8 :: soiltemp_10d      
      real*8 :: airtemp_10d
      real*8 :: paw_10d

      real*8 :: gdd
      real*8 :: ncd
      real*8 :: ld

      real*8 :: zweight
      
      sandfrac = pp%cellptr%soil_texture(1)
      clayfrac = pp%cellptr%soil_texture(3)
      smpsat = -10.d0 * ( 10.d0**(1.88d0-0.0131d0*sandfrac) )
      bch = 2.91d0 + 0.159d0*clayfrac
      watsat = 0.489d0 - 0.00126d0*sandfrac 
      watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch) 
  
c$$$      smpsc = -200000.d0 !grass

      airtemp =  pp%cellptr%TairC
      wat = pp%cellptr%Soilmoist(1)  !assign these only to top 30 cm for now -PK 3/16/07
      soiltemp = pp%cellptr%Soiltemp(1)
      coszen = pp%cellptr%CosZen
      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      paw_10d = pp%cellptr%paw_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd
      ld =  pp%cellptr%ld

      zweight=exp(-1.d0/(10.d0*86400.d0/dtsec)) !weighting factor for 10-day running average

      cop => pp%tallest
      do while(ASSOCIATED(cop))
         
         !daily carbon balance for growth
         cop%CB_d =  cop%CB_d + cop%NPP*dtsec/cop%n

         !10-day running average of Soil Temperature
         soiltemp_10d=zweight*soiltemp_10d+(1.0d0-zweight)*soiltemp
      
         !10-day running average of Air Temperature
         airtemp_10d=zweight*airtemp_10d+(1.0d0-zweight)*airtemp
   
         !10-day running average of Plant Available Water  
         paw = min( max(wat-watdry,0.d0) 
     &         /(watsat-watdry), 1.d0)   
         paw_10d=zweight*paw_10d+(1.0d0-zweight)*paw

c$$$         soilmetric = min (smpsc, smpsat*paw**(-bch))      

         !GDD & NCD - Once a day
         if (dailyupdate) then
            !Growing degree days
            if (airtemp_10d .ge. airtemp_par) then
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


!**************************************************************
!* Update photosynthetic acclimation factor for evergreen veg
!**************************************************************

         if (config%do_frost_hardiness) then
            if (pfpar(cop%pft)%phenotype==1) then !evergreen = 1
               call photosyn_acclim(dtsec,airtemp_10d,cop%Sacclim) 
            else
               cop%Sacclim = UNDEF
            endif
	 else
             cop%Sacclim = UNDEF
         endif


         cop => cop%shorter 
 
      end do        

      pp%cellptr%soiltemp_10d = soiltemp_10d
      pp%cellptr%airtemp_10d = airtemp_10d
      pp%cellptr%paw_10d = paw_10d
      pp%cellptr%gdd = gdd
      pp%cellptr%ncd = ncd
      pp%cellptr%ld =  ld

      end subroutine clim_stats
      !*********************************************************************   
      subroutine pheno_update(dtsec, pp)
!@sum Update statstics for phneology_update    
      real*8,intent(in) :: dtsec           !dt in seconds
      type(patch) :: pp
      !--Local-----
      type(cohort), pointer :: cop
 
      !temperature constrain for cold-deciduous PFTs (Botta et al. 1997)
      ! gdd_par1/2/3: paramters to estimate the threshold for gdd
      ! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
      real*8, parameter :: gdd_par1 = -68.d0 
      real*8, parameter :: gdd_par2 = 638.d0
      real*8, parameter :: gdd_par3 = -0.01d0 
      !gdd_length
      real*8, parameter :: gdd_length = 50.d0
      !airt_threshold
      real*8, parameter :: airt_max_w = 15.d0
      real*8, parameter :: airt_min_w = 5.d0
      !ld_threshold (minute): light length constraint for cold-deciduous woody PFTs (White et al. 1997)
      real*8, parameter :: ld_threshold = 655.d0
      !tsoil_threshold1, tsoil_threshold2 : soil temperature constraint for cold-deciduous woody PFTs (White et al. 1997)
c$$$      real*8, parameter :: tsoil_threshold1 = 11.15d0
c$$$      real*8, parameter :: tsoil_threshold2 = 2.d0
      !ddfacu: the rate of leaf fall (1/day) (IBIS)
      real*8, parameter :: ddfacu = 1.d0/15.d0
      !paw_threshold: water threshold (unitless)       
      real*8, parameter :: paw_max_w = 0.3d0
      real*8, parameter :: paw_min_w = 0.0d0
      real*8, parameter :: paw_res_w = 0.25d0
      real*8, parameter :: paw_max_h = 0.3d0
      real*8, parameter :: paw_min_h = 0.2d0
      real*8, parameter :: paw_res_h = 1.0d0
  
      real*8 :: gdd_threshold

      integer :: pft
      integer :: phenotype

      real*8 :: soiltemp_10d      
      real*8 :: airtemp_10d
      real*8 :: paw_10d

      real*8 :: gdd
      real*8 :: ncd
  
      real*8 :: phenofactor_c
      real*8 :: phenofactor_d
      real*8 :: phenofactor

      integer :: phenostatus
      
      real*8 :: pawmin,pawmax,pawres

      logical :: temp_limit, light_limit, water_limit
      logical :: woody
      logical :: fall

      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      paw_10d = pp%cellptr%paw_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd

      gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)  

      pp%cellptr%light_prev=pp%cellptr%light
      pp%cellptr%light=pp%cellptr%CosZen
      if (pp%cellptr%light .le. pp%cellptr%light_prev ) then
         fall = .true.
      else
         fall = .false.
      end if
      
      cop => pp%tallest
      do while(ASSOCIATED(cop))       
      
         phenofactor_c=cop%phenofactor_c
         phenofactor_d=cop%phenofactor_d
         phenofactor=cop%phenofactor
         phenostatus=cop%phenostatus
         pft=cop%pft
         phenotype=pfpar(pft)%phenotype

         if (phenotype .eq. 1) then !Cold-deciduous woody plants
            temp_limit = .true.
            water_limit = .false.
         else if (phenotype .eq. 2) then !Drought-deciduous woody plants
            temp_limit = .true.
            water_limit = .true.
         else if (phenotype .ge. 3) then !Herbaceous plants
            temp_limit = .true.
            water_limit = .true.            
         end if

         if ((pft .ge. 11) .and. (pft .le. 15)) then !Herbaceous plants 
            woody = .false.
         else
            woody = .true.
         end if
      
         if (temp_limit .and. woody) then
            if ((.not. fall) .and.
     &         (phenostatus.le.2).and.(gdd.gt.gdd_threshold)) then
               phenofactor_c = min (1.d0,(gdd-gdd_threshold)/gdd_length)
               if (phenofactor_c .lt. 1.d0) then
                  phenostatus = 2
               else
                  phenostatus = 3
               end if
            end if
            if (fall .and. 
     &         (phenostatus.ge.3).and.(airtemp_10d.lt.airt_max_w)) then
               phenofactor_c = min(phenofactor_c,max(0.d0,
     &            (airtemp_10d-airt_min_w)/(airt_max_w-airt_min_w)))
               if (phenofactor_c .eq. 0.d0) then
                  phenostatus = 1
                  ncd = 0.d0             !zero-out
                  gdd = 0.d0
               else  
                  phenostatus = 4
               end if 
            end if
         end if
         
         if (temp_limit .and. (.not. woody) )then
            phenofactor_c = 1.d0 !not yet implemented
         end if            
                 
         if (water_limit .and. woody .and. (phenostatus.ge.2)) then
            phenofactor_d = min(1.d0,max(0.d0,
     &         ((paw_10d-paw_min_w)/(paw_max_w-paw_min_w))**paw_res_w))
         end if
  
         if (water_limit .and. (.not. woody)) then
            if ((phenostatus.le.2).and.(paw_10d.gt.paw_min_h))then
               phenofactor_d = min(1.d0,
     &          ((paw_10d-paw_min_h)/(paw_max_h-paw_min_h))**paw_res_h)
               if (phenofactor_d .lt. 1.d0) then
                  phenostatus = 2
               else
                  phenostatus = 3
               end if
            else if ((phenostatus.ge.3).and.(paw_10d.lt.paw_max_h))then
               phenofactor_d = max(0.d0,
     &          ((paw_10d-paw_min_h)/(paw_max_h-paw_min_h))**paw_res_h)
               if (phenofactor_d .eq. 0.d0) then
                  phenostatus = 1
               else
                  phenostatus = 4
               end if 
            end if
         end if    
 
         if (.not.temp_limit) phenofactor_c = 1.d0
         if (.not.water_limit) phenofactor_d = 1.d0
      
         phenofactor = phenofactor_c * phenofactor_d
         cop%phenofactor_c=phenofactor_c
         cop%phenofactor_d=phenofactor_d
         cop%phenofactor=phenofactor
         cop%phenostatus=phenostatus

         
    
         cop => cop%shorter 
      
      end do   

!      pp%cellptr%ld =  0.d0     !zero-out every day
      pp%cellptr%gdd = gdd      !it should be cohort level!
      pp%cellptr%ncd = ncd


      end subroutine pheno_update
      !*********************************************************************   
      subroutine veg_update(dtsec,pp)
!@sum Update the vegetation state and carbon pools.
!@sum i.e., LAI, senescefrac, DBH, height &
!@sum carbon pools of foliage, sapwood, fineroot, hardwood, coarseroot
      
      use ent_prescr_veg
      implicit none
      type(patch),pointer :: pp 
      type(cohort), pointer :: cop
      real*8 :: dtsec
      integer :: pft
      ! phenofactor phenological elongation factor [0,1] (unitless)
      real*8 :: phenofactor
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
      real*8 :: laipatch
      real*8 :: C_fol_max
      real*8 :: qf
      real*8 :: dCrepro

      logical :: woody
      logical :: annual

      laipatch = 0.d0
      cop => pp%tallest

      do while(ASSOCIATED(cop))
               
         pft = cop%pft
         phenofactor = cop%phenofactor
         dbh = cop%dbh
         h = cop%h
         nplant = cop%n

         if ((pft .ge. 11) .and. (pft .le. 15)) then !Herbaceous plants 
            woody = .false.
         else !Woody
            woody = .true.
         end if

         if (pft .eq. 13) then !C3 annual grass
            annual = .true.
         else 
            annual = .false. 
         end if
       
         if (woody) then 
            qsw = pfpar(pft)%sla*iqsw
         else 
            qsw=0.0d0           !no allocation to the wood
         end if
         
         if (annual) then
            qf = q * phenofactor 
         else
            qf = q 
         endif 
         ialloc = phenofactor+qf+h*qsw
         
         Cactive = cop%C_froot + cop%C_fol + cop%C_sw
         if (woody) then  
            Cactive_max=dbh2Cfol(pft,dbh)*ialloc
         else
            Cactive_max=height2Cfol(pft,5.0d0)*ialloc
           !no allometric constraints in the carbon allocation, 
           !then 5.0d0 is arbitrary number (must be tall enough, then 5m)
         end if
         C_lab = cop%C_lab*1000.d0 
         C_fol_old = cop%C_fol
         C_fol = cop%C_fol
         C_froot = cop%C_froot
         C_hw = cop%C_hw
         C_croot = cop%C_croot
         Cdead = C_hw + C_croot   
         
         CB_d = cop%CB_d*1000.d0
         
        !allocate the labile carbon to the active/dead pools or 
        !store the carbon for the reproduction use
         
        !active growth
         call  growth_cpools_active(phenofactor,Cactive_max,
     &        CB_d,Cactive,C_lab)
         if (phenofactor .ne. 0.d0) then    
            C_fol = phenofactor* Cactive/ialloc !C_fol should be updated for structural growth!
         else
            C_fol = 0.d0
         endif

         !structural/active/reproductive
         call growth_cpools_structural(pft,dbh,h,qsw,qf,
     &        C_fol, Cactive_max, 
     &        Cactive,C_lab,Cdead,dCrepro) 
         if (phenofactor .ne. 0.d0) then    
            C_fol = phenofactor* Cactive/ialloc
         else
            C_fol = 0.d0
         endif

         !senescening
         call senesce_cpools(annual,C_fol_old,C_fol,C_lab,Cactive,
     &        senescefrac)
         
         !update the active and structural pools   
         if (phenofactor .ne. 0.d0) then    
            cop%C_fol = phenofactor* Cactive/ialloc
         else
            cop%C_fol = 0.d0
         endif
         cop%C_froot = Cactive *(qf/ialloc)
         cop%C_sw = Cactive *(h*qsw/ialloc)
         cop%C_hw = Cdead * hw_fract
         cop%C_croot = Cdead * (1-hw_fract)                  

         !update the labile and reproduction
         cop%C_lab = 0.001d0*C_lab !g/individual -> kg/ind.
         cop%pptr%Reproduction(cop%pft) = 
     &        cop%pptr%Reproduction(cop%pft)+ dCrepro*cop%n

         !update the plant size 
         if (woody) then
            cop%dbh = Cdead2dbh(pft,Cdead)
            cop%h = dbh2height(pft,cop%dbh)
         else
            cop%dbh = 0.0d0
            cop%h = Cfol2height(pft,cop%C_fol)
         end if
         
         !update the senescefrac & LAI & Nitrogen
         cop%senescefrac = senescefrac
         cop%LAI=cop%C_fol/1000.0d0*pfpar(pft)%sla*cop%n
         if (cop%LAI .lt. EPS) cop%LAI=EPS
         cop%Ntot = cop%nm * cop%LAI 
         laipatch = laipatch + cop%lai  

         !reproduction
c$$$         if (annual) then
c$$$            if (phenofactor .gt. 0.d0 .AND. cop%h .lt. 0.05d0) then
c$$$               cop%h = 0.05d0
c$$$               cop%C_fol = height2Cfol(pft,cop%h)
c$$$               cop%LAI=cop%n*pfpar(pft)%sla*
c$$$     &              (height2Cfol(pft,cop%h)/1000.0d0) 
c$$$               cop%C_froot = q*cop%C_fol
c$$$            else if(phenofactor .eq. 0.d0) then
c$$$               cop%C_lab =0.d0
c$$$            end if
c$$$         end if

                
         !zero-out the daily accumulated carbon 
         cop%CB_d = 0.d0    

         
!         write(990,'(i4,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g,g)')
!     &        cop%pft,  
!     &        cop%phenofactor_c,
!     &        cop%phenofactor_d,
!     &        cop%phenostatus,
!     &        cop%LAI,
!     &        cop%C_fol,
!     &        cop%C_lab,
!     &        cop%C_sw,
!     &        cop%C_hw,
!     &        cop%C_froot,
!     &        cop%C_croot,
!     &        cop%NPP,
!     &        cop%dbh,
!     &        cop%h,
!     &        cop%CB_d,
!     &        cop%senescefrac,
!     &        cop%cellptr%airtemp_10d,
!     &        cop%cellptr%soiltemp_10d,
!     &        cop%cellptr%paw_10d,
!     &        cop%cellptr%gdd,
!     &        cop%cellptr%ncd,
!     &        cop%cellptr%CosZen

         cop => cop%shorter 
      end do
      pp%LAI = laipatch
      pp%cellptr%ld =  0.d0  
      end subroutine veg_update

      !*********************************************************************
      subroutine growth_cpools_active(phenofactor, 
     &     Cactive_max,CB_d,Cactive,C_lab)

      real*8, intent(in) :: phenofactor
      real*8, intent(in) :: Cactive_max
      real*8, intent(in) :: CB_d
      real*8, intent(inout) :: Cactive
      real*8, intent(inout) :: C_lab
    
      real*8 :: Cactive_pot
      real*8 ::dC_lab

                              
      if ((phenofactor .eq. 0.d0) .AND. (CB_d .gt. 0.d0)) then
         dC_lab = 0.d0 !store the carbon in the labile
      else
         if (phenofactor .gt. 0.d0) then
            !Cactive_max (max. allowed pool size according to the DBH)
            !Cactive_pot (current size + daily accumulated carbon)  
            !Cactive (cuurent size)
            Cactive_pot = Cactive + CB_d          
            if (Cactive .gt. Cactive_max) then 
              dC_lab = Cactive - Cactive_max
            else if (Cactive_pot .gt. Cactive_max) then
              dC_lab= Cactive - Cactive_max  
            else
              dC_lab = -CB_d
            end if
         else
           !negative CB_d 
           !compensate it with the active carbon rather than the labile carbon.
!           dC_lab_growth_a1 = -CB_d 
           !compensate it with the labile carbon rather than the active carbon.
           dC_lab = 0.d0
         end if
      end if
      C_lab = C_lab + dC_lab
      Cactive = Cactive - dC_lab 


      end subroutine growth_cpools_active


      !*********************************************************************
      subroutine growth_cpools_structural(pft,dbh,h,qsw,qf,
     &      C_fol,Cactive_max,Cactive,
     &      C_lab, Cdead,dCrepro)

      use ent_prescr_veg

      integer, intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: h
      real*8, intent(in) :: qsw
      real*8, intent(in) :: qf
      real*8, intent(in) :: C_fol
      real*8, intent(in) :: Cactive_max
      real*8, intent(inout) :: Cactive
      real*8, intent(inout) :: C_lab
      real*8, intent(inout) :: Cdead
      
      real*8, intent(out) :: dCrepro
      real*8 :: dCdead 

      real*8 :: dCactive
      !r_fract: fraction of excess c going to seed reproduction 
      real*8 :: r_fract 
      !c_fract: fraction of excess c going to clonal reproduction
      real*8 :: c_fract 
      !gr_fract: fraction of excess c going to structural growth
      real*8 :: gr_fract  
      real*8 :: qs
      real*8 :: dCfoldCdead
      real*8 :: dCfrootdCdead
      real*8 :: dHdCdead
      real*8 :: dCswdCdead

      !herbaceous
      r_fract = 0.3d0
      c_fract = 0.7d0 
      if (pft .ge. 11 .and. pft .le. 15) then !grasses/crops(?) 
         if (C_lab .gt. 0.d0 )then
            qs = 0.d0  !no structural pools
            gr_fract = 1.d0 - (r_fract + c_fract)
         else
            qs = 0.d0 
            gr_fract = 1.d0
         end if
      !woody
      r_fract = 0.3d0
      c_fract = 0.0d0 
      else
         if (C_fol .gt. 0.d0 .and. 
     &       Cactive .ge. Cactive_max .and. C_lab .gt. 0.d0) then
            if (dbh .le. maxdbh(pft))then
               dCfoldCdead = dDBHdCdead(pft,Cdead)/dDBHdCfol(pft,C_fol)
               dCfrootdCdead = qf *dCfoldCdead
               dHdCdead = dHdDBH(pft, dbh)  * dDBHdCdead(pft,Cdead) 
               dCswdCdead = qsw*
     &                      (h*dCfoldCdead + C_fol*dHdCdead)
               qs = 1.d0 / (dCfoldCdead + dCfrootdCdead + dCswdCdead)
               gr_fract = 1.d0 - r_fract
             else
               qs = 1.d0
               gr_fract = 1.d0 - (r_fract + c_fract)
             end if
         else
            qs = 0.d0
            gr_fract = 1.d0
         end if
      end if
       
      dCdead = gr_fract * qs / (1.d0 + qs) * C_lab
      dCactive = gr_fract / (1.d0 + qs) * C_lab
      dCrepro =  ( 1.d0 - gr_fract ) * C_lab 

      C_lab = 0.d0
      Cdead = Cdead + dCdead
      Cactive = Cactive + dCactive

      end subroutine growth_cpools_structural
      !*********************************************************************
      subroutine senesce_cpools(annual, C_fol_old,C_fol,C_lab,Cactive,
     &                          senescefrac)
      logical, intent(in) :: annual
      real*8, intent(in) :: C_fol_old
      real*8, intent(in) :: C_fol
      real*8, intent(inout) :: C_lab
      real*8, intent(inout) :: Cactive
      real*8, intent(out) :: senescefrac
   
      real*8 :: dC_fol
      real*8 :: dCactive
      real*8 :: dC_lab

      senescefrac = 0.d0
      if (C_fol_old .gt. C_fol) then 
         !with senescening, the part of foliage carbon goes to the litter pool
         ! and the remained goes to the labile.
         dC_fol = C_fol - C_fol_old !negative
         dCactive = dC_fol 
         dC_lab = - dCactive * l_fract
         if (C_fol_old .ne. 0.d0 ) then 
            senescefrac = -dC_fol * (1.d0 - l_fract) 
     &           /C_fol_old
         endif    
      else
         dCactive = 0.d0
         dC_lab = 0.d0         
      endif

      C_lab = C_lab + dC_lab  
      Cactive= Cactive + dCactive
      
      end subroutine senesce_cpools

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

      subroutine photosyn_acclim(dtsec,Ta,Sacc)
!@sum Model for state of acclimation/frost hardiness based 
!@sum for boreal coniferous forests based on Repo et al (1990), 
!@sum Hanninen & Kramer (2007),and  Makela et al (2006)
      implicit none
      real*8,intent(in) :: dtsec ! time step size [sec]
      real*8,intent(in) :: Ta ! air temperature [deg C]
      real*8,intent(inout) :: Sacc ! state of acclimation [deg C]

      !----Local-----
      real*8,parameter :: tau_inv = 2.22222e-6 
                          ! inverse of time constant of delayed
                          ! response to ambient temperature [sec] = 125 hr 
                          ! Makela et al (2004) for Scots pine

!      Use a first-order Euler scheme
       Sacc = Sacc + dtsec*(tau_inv)*(Ta - Sacc) 

!!     Predictor-corrector method requires temperature from next timestep
!       Sacc_old = Sacc
!       Sacc = Sacc_old + dtsec*(1/tau_acclim)*(Ta - Sacc ) 
!       Sacc = Sacc_old + ((1/tau_acclim)*(Ta - Sacc_old)+
!     &                     (1/tau_acclim)*(Ta_next - Sacc))*0.5d*dtsec

      end subroutine photosyn_acclim

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
      real*8 function height2Cfol(pft,height)
      integer,intent(in) :: pft
      real*8, intent(in) :: height
      real*8,parameter :: h1Cf = 5.0d0  !1.660d0
      real*8,parameter :: h2Cf = 1.20d0 !1.500d0
      real*8,parameter :: nplant = 3000.0d0
      height2Cfol=(1.0d0/C2B) *h1Cf 
     &        *((height*100.0d0)**h2Cf)/nplant

      end function height2Cfol
!*************************************************************************
      real*8 function Cfol2height(pft,Cfol)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol !gC/pool/plant
      real*8,parameter :: h1Cf = 5.0d0  !1.660d0
      real*8,parameter :: h2Cf = 1.20d0 !1.500d0
      real*8,parameter :: nplant = 3000.0d0

      Cfol2height=exp(log(Cfol*C2B/h1Cf*nplant)/h2Cf)/100.0d0

      end function Cfol2height
!*************************************************************************
      real*8 function dDBHdCdead(pft,Cdead)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

      dDBHdCdead=(C2B/1000.0d0/pfpar(pft)%b1Cd)**(1.0d0/pfpar(pft)%b2Cd)
     &          *Cdead**((1.0d0/pfpar(pft)%b2Cd)-1.0d0)
     &          /pfpar(pft)%b2Cd

      end function dDBHdCdead
!*************************************************************************
      real*8 function dDBHdCfol(pft,Cfol)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol

      dDBHdCfol=(C2B/1000.0d0/pfpar(pft)%b1Cf)**(1.0d0/pfpar(pft)%b2Cf)
     &          *Cfol**((1.0d0/pfpar(pft)%b2Cf)-1.0d0)
     &          /pfpar(pft)%b2Cf

      end function dDBHdCfol
!*************************************************************************
      real*8 function dHdDBH(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dHdDBH = - pfpar(pft)%b1Ht*pfpar(pft)%b2Ht
     &         * exp(pfpar(pft)%b2Ht * dbh)
      end function dHdDBH
!*************************************************************************
      end module phenology
