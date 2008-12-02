      module ent_prescribed_updates
!@sum Routines for updating prescribed vegetation. These routines
!@+   work on entcell level or lower.

!#define DEBUG TRUE  !NYK

      use ent_types

      implicit none
      private

      !public entcell_update_lai, entcell_update_albedo
      !public entcell_update_shc
      public entcell_vegupdate

      contains

      subroutine entcell_update_lai_poolslitter( ecp,
     i    laidata)
!@sum sets prescribed LAI over the cell and update  C_fol, C_froot, senescefrac
      !use ent_prescrveg, only : prescr_plant_cpools
      use phenology, only : litter_patch
      type(entcelltype) :: ecp
      real*8,intent(in) :: laidata(N_PFT) !@var LAI for all PFT's 
      !-----Local---------
      real*8, parameter :: dt = 1800.d0 !seconds, time step
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 laipatch, lai_old,lai_new
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.

      laipatch = 0.d0
      Clossacc(:,:,:) = 0.d0
      pp => ecp%oldest      
      do while ( associated(pp) )
        laipatch = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
          lai_old = cop%LAI
          lai_new = laidata(cop%pft)
          !* Update biomass pools, senescefrac, accumulate litter
          call prescr_veglitterupdate_cohort(cop,lai_new,Clossacc)

          !* Calculate senescence factor for next time step litterfall routine.
          !* OLD SCHEME - replaced by presc_veglitterupdate_cohort - NYK
!          if (cop%LAI.gt.0d0) then !Prescribed senescence fraction.
!            cop%senescefrac = max(0.d0,(lai_old-cop%LAI)/lai_old)
!          else
!            cop%senescefrac = 0.d0
!          endif

          laipatch = laipatch + cop%lai
          cop => cop%shorter
        enddo
        call litter_patch(pp,Clossacc) 
        pp%LAI = laipatch
        pp => pp%younger
      enddo

      end subroutine entcell_update_lai_poolslitter


      subroutine entcell_update_height( ecp,
     i    hdata,init)
!@sum sets prescribed LAI over the cell
      use ent_prescr_veg, only : prescr_plant_cpools, popdensity,
     &     ED_woodydiameter
      use ent_pfts
      type(entcelltype) :: ecp
      real*8,intent(in) :: hdata(N_PFT) !@var LAI for all PFT's 
      logical,intent(in) :: init
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 :: cpool(N_BPOOLS)
      real*8 :: C_sw_old, C_hw_old,C_croot_old

      cpool(:) = 0.d0
      pp => ecp%oldest      
      do while ( associated(pp) )
        cop => pp%tallest
        do while ( associated(cop) )
          C_sw_old = cop%C_sw
          C_hw_old = cop%C_hw
          C_croot_old = cop%C_croot
          cop%h = hdata(cop%pft)

          if (pfpar(cop%pft)%woody) then !update dbhuse
            cop%dbh = ED_woodydiameter(cop%pft,cop%h)
            if (init) then !Set population density
              cop%n = popdensity(cop%pft,cop%dbh) 
            endif
          endif

          !* Update biomass pools except for C_lab.
          call prescr_plant_cpools(cop%pft, cop%lai, cop%h, 
     &         cop%dbh, cop%n, cpool )
          cop%C_fol = cpool(FOL)
          cop%C_sw = cpool(SW)
          cop%C_hw = cpool(HW)
          cop%C_froot = cpool(FR)
          cop%C_croot = cpool(CR)
          !* Update C_lab
          cop%C_lab = cop%C_lab - max(0.d0,cop%C_hw - C_hw_old)
     &         - max(0.d0,cop%C_sw - C_sw_old) 
     &         - max(0.d0,cop%C_croot-C_croot_old)
          cop => cop%shorter
        enddo
        pp => pp%younger
      enddo
      !summarize_patch - called outside this routine
      end subroutine entcell_update_height


      subroutine entcell_update_albedo( ecp,
     i    albedodata)
!@sum sets prescribed albedo in vegetated patches of the cell (skips bare soil)
!@+   This subroutine assumes one cohort per patch !!!
      type(entcelltype) :: ecp
      real*8,intent(in) :: albedodata(N_BANDS,N_PFT) !@var albedo for all PFTs 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch

      pp => ecp%oldest      
      do while ( associated(pp) )
        ! update albedo of vegetated patches only
        if ( associated(pp%tallest) ) then ! have vegetation
          if ( pp%tallest%pft > N_PFT .or.  pp%tallest%pft < 1 )
     &         call stop_model("entcell_update_albedo: bad pft",255)
          pp%albedo(1:N_BANDS) = albedodata(1:N_BANDS, pp%tallest%pft)
        endif
        pp => pp%younger
      enddo

      end subroutine entcell_update_albedo

      subroutine entcell_update_crops( ecp,
     i    cropsdata)
!@sum sets prescribed albedo in vegetated patches of the cell (skips bare soil)
!@+   This subroutine assumes one cohort per patch !!!
      type(entcelltype) :: ecp
      real*8,intent(in) :: cropsdata !@var albedo for all PFTs 
      !-----Local---------

      !!! nothing here yet ...

      end subroutine entcell_update_crops


      subroutine entcell_update_shc( ecp )
      use ent_prescr_veg, only : prescr_calc_shc
      use entcells, only : entcell_extract_pfts
      type(entcelltype) :: ecp
      !-----Local---------
      real*8 vdata(N_COVERTYPES) ! needed for a hack to compute canopy

      vdata(:) = 0.d0
      call entcell_extract_pfts( ecp, vdata(2:) )
      ecp%heat_capacity=prescr_calc_shc(vdata)

      end subroutine entcell_update_shc
     

      subroutine entcell_vegupdate(entcell, hemi, jday
     &     ,do_giss_phenology, do_giss_lai, do_giss_albedo
     &     ,laidata, hdata, albedodata, cropsdata, init )
!@sum updates corresponding data on entcell level (and down). 
!@+   everything except entcell is optional. coordinate-dependent
!@+   is given pointer attribute to provide a way to tell the 
!@+   program that an argument is actually optional and missing
!@+   (see how it is used in ent_prescribe_vegupdate)
      use patches, only : summarize_patch
      use entcells,only : summarize_entcell!,entcell_extract_pfts
      use ent_prescr_veg, only : prescr_calc_shc, prescr_veg_albedo
      implicit none
      type(entcelltype) :: entcell
      integer,intent(in) :: jday
      integer,intent(in) :: hemi
      logical, intent(in) :: do_giss_phenology
      logical, intent(in) :: do_giss_lai
      logical, intent(in) :: do_giss_albedo
      real*8,  pointer :: laidata(:)  !Array of length N_PFT
      real*8,  pointer :: hdata(:)  !Array of length N_PFT
      real*8,  pointer :: albedodata(:,:)
      real*8,  pointer :: cropsdata
      logical, intent(in) :: init
      !----Local------
      type(patch),pointer :: pp

      if (.not.do_giss_lai) then  !* Prescribed non-GISS veg structure.
      ! update with external data first
        if ( associated(laidata) )
     &       call entcell_update_lai_poolslitter(entcell, laidata)

        if ( associated(hdata) )
     &       call entcell_update_height(entcell, hdata, init)
        
        if ( associated(albedodata) )
     &       call entcell_update_albedo(entcell, albedodata)
        
        if ( associated(cropsdata) )
     &       call entcell_update_crops(entcell, cropsdata)
      endif
      ! and then do GISS phenology if required
      if ( do_giss_phenology ) then  !do_giss_phenology is redundant with do_giss_lai.
        if ( hemi<-2 .or. jday <-2 )
     &       call stop_model("entcell_vegupdate: needs hemi,jday",255)
        pp => entcell%oldest
        do while (ASSOCIATED(pp))
          !* LAI, SENESCEFRAC *!
          call prescr_phenology(jday,hemi, pp, do_giss_lai)
          !* ALBEDO *!
          if ( ASSOCIATED(pp%tallest).and.do_giss_albedo ) then ! update if have vegetation or not prognostic albedo
            call prescr_veg_albedo(hemi, pp%tallest%pft, 
     &           jday, pp%albedo)
          endif
          !call summarize_patch(pp) !* Redundant because summarize_entcell is called.
          pp => pp%younger
        end do
      endif

      call summarize_entcell(entcell)

      call entcell_update_shc(entcell)
cddd      vdata(:) = 0.d0
cddd      call entcell_extract_pfts(entcell, vdata(2:) )
cddd      entcell%heat_capacity=GISS_calc_shc(vdata)

      !if (YEAR_FLAG.eq.0) call ent_GISS_init(entcellarray,im,jm,jday,year)
      !!!### REORGANIZE WTIH ent_prog.f ####!!!
      
      end subroutine entcell_vegupdate


      subroutine prescr_phenology(jday,hemi,pp,do_giss_lai)
      !* DAILY TIME STEP *!
      !* Calculate new LAI, biomass poos, and senescefrac 
      !* for given jday, for prescr vegetation. *!
      use ent_pfts
      use ent_prescr_veg, only : prescr_calc_lai,prescr_plant_cpools,
     &     prescr_veg_albedo
      use phenology, only : litter_patch
      implicit none
      integer,intent(in) :: jday !Day of year.
      integer,intent(in) :: hemi !@var hemi -1: S.hemisphere, 1: N.hemisphere
      type(patch),pointer :: pp
      logical, intent(in) :: do_giss_lai
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: laipatch
      real*8 :: cpool(N_BPOOLS)
      real*8 :: lai_new
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      real*8 :: resp_growth_root !g-C/individ/s
      real*8 :: resp_growth_patch, resp_growth_root_patch !kg-C/m/s

      if (ASSOCIATED(pp)) then
        !* LAI *AND* BIOMASS - carbon pools *!
        laipatch = 0.d0         !Initialize for summing
        cpool(:) = 0.d0
        Clossacc(:,:,:) = 0.d0
        resp_growth_patch = 0.d0
        resp_growth_root = 0.d0
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          if ( do_giss_lai ) then
            lai_new = prescr_calc_lai(cop%pft+COVEROFFSET, jday, hemi)
            call prescr_veglitterupdate_cohort(cop,
     &           lai_new,Clossacc)
            laipatch = laipatch + cop%LAI
          endif
         !* Summarize for patch level. - NOTE:This is total flux, not just the growth increment.
         resp_growth_patch = resp_growth_patch + cop%R_auto
         resp_growth_root_patch = resp_growth_root_patch + cop%R_root

          cop => cop%shorter
        end do
        pp%R_auto = resp_growth_patch !Total flux includes growth increment.
        pp%R_root = resp_growth_root !Total flux includes growth increment.
        pp%NPP = pp%GPP - resp_growth_patch
        call litter_patch(pp,Clossacc) !Update Tpools following litter accumulation. Daily time step
        pp%LAI = laipatch

        !* ALBEDO *! - Moved these lines to entcell_vegupdate.
!        if ( ASSOCIATED(pp%tallest) ) then ! update if have vegetation
!          call prescr_veg_albedo(hemi, pp%tallest%pft, 
!     &         jday, pp%albedo)
!        endif
      endif
      end subroutine prescr_phenology

!******************************************************************************
      subroutine prescr_veglitterupdate_cohort(cop,lai_new,Clossacc)
!@sum prescr_veglitterupdate_cohort Given new LAI, update cohort biomass
!@sum pools, litter, senescefrac. - NYK
      use phenology, only : litter_cohort
      use ent_prescr_veg, only :prescr_plant_cpools
      implicit none
      type(cohort),pointer :: cop
      real*8,intent(in) :: lai_new
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
       !-------local-----
      real*8 :: cpool(N_BPOOLS)
      real*8 :: lai_old
      real*8 :: C_fol_old,C_froot_old,C_hw_old,C_croot_old
      real*8 :: resp_growth, resp_growth_root
      real*8 :: i2a  !Convert g-C/individual to kg-C/m^2
#ifdef DEBUG
      integer :: i
      real*8 :: Csum
#endif

      lai_old = cop%LAI
      C_fol_old = cop%C_fol
      C_froot_old = cop%C_froot
      C_hw_old = cop%C_hw
      C_croot_old = cop%C_croot
      cop%LAI = lai_new
 
      !* Update biomass pools except for C_lab.
      call prescr_plant_cpools(cop%pft, cop%lai, cop%h, 
     &     cop%dbh, cop%n, cpool )
      cop%C_fol = cpool(FOL)
      cop%C_sw = cpool(SW)
      cop%C_hw = cpool(HW)
      cop%C_froot = cpool(FR)
      cop%C_croot = cpool(CR)

      !* Update litter, C_lab, senescefrac
      call litter_cohort(SDAY,
     i     C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     cop,Clossacc)
        !*## DEBUG ##*!
#ifdef DEBUG
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo
      write(999,*) Csum
#endif
        
      
      cop%Ntot = cop%nm * cop%LAI !This should eventually go into N allocation routine if dynamic nm.

      i2a = 1d-3*cop%n          !Convert g-C/individual to kg-C/m^2
      cop%R_auto = cop%R_auto + i2a*resp_growth
      cop%R_root = cop%R_root + i2a*resp_growth_root
      cop%NPP = cop%NPP - i2a*resp_growth
 
      end subroutine prescr_veglitterupdate_cohort
!******************************************************************************

      end module ent_prescribed_updates
