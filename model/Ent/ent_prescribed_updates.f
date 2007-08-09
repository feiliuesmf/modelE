      module ent_prescribed_updates
!@sum Routines for updating prescribed vegetation. These routines
!@+   work on entcell level or lower.
      use ent_types

      implicit none
      private

      !public entcell_update_lai, entcell_update_albedo
      !public entcell_update_shc
      public entcell_vegupdate

      contains

      subroutine entcell_update_lai( ecp,
     i    laidata)
!@sum sets prescribed LAI over the cell
      !use ent_prescrveg, only : prescr_plant_cpools
      type(entcelltype) :: ecp
      real*8,intent(in) :: laidata(N_PFT) !@var LAI for all PFT's 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 laipatch, lai_old

      pp => ecp%oldest      
      do while ( associated(pp) )
        laipatch = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
          lai_old = cop%LAI
          cop%lai = laidata(cop%pft)
!          write(777,*) __FILE__,__LINE__,cop%LAI
          laipatch = laipatch + cop%lai

          !* Calculate senescence factor for next time step litterfall routine.
          if (cop%LAI.gt.0d0) then !Prescribed senescence fraction.
            cop%senescefrac = max(0.d0,(lai_old-cop%LAI)/lai_old)
          else
            cop%senescefrac = 0.d0
          endif

          cop => cop%shorter
        enddo
        pp%LAI = laipatch
        pp => pp%younger
      enddo

      end subroutine entcell_update_lai


      subroutine entcell_update_height( ecp,
     i    hdata)
!@sum sets prescribed LAI over the cell
      !use ent_prescrveg, only : prescr_plant_cpools
      type(entcelltype) :: ecp
      real*8,intent(in) :: hdata(N_PFT) !@var LAI for all PFT's 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort

      pp => ecp%oldest      
      do while ( associated(pp) )
        cop => pp%tallest
        do while ( associated(cop) )
          cop%h = hdata(cop%pft)
          cop => cop%shorter
        enddo
        pp => pp%younger
      enddo

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
     &     ,do_giss_phenology, do_giss_lai
     &     ,laidata, hdata, albedodata, cropsdata )
!@sum updates corresponding data on entcell level (and down). 
!@+   everything except entcell is optional. coordinate-dependent
!@+   is given pointer attribute to provide a way to tell the 
!@+   program that an argument is actually optional and missing
!@+   (see how it is used in ent_prescribe_vegupdate)
      use patches, only : summarize_patch
      use entcells,only : summarize_entcell!,entcell_extract_pfts
      use ent_prescr_veg, only : prescr_calc_shc

      type(entcelltype) :: entcell
      integer,intent(in) :: jday
      integer,intent(in) :: hemi
      logical, intent(in) :: do_giss_phenology
      logical, intent(in) :: do_giss_lai
      real*8,  pointer :: laidata(:)  !Array of length N_PFT
      real*8,  pointer :: hdata(:)  !Array of length N_PFT
      real*8,  pointer :: albedodata(:,:)
      real*8,  pointer :: cropsdata
      !----Local------
      type(patch),pointer :: pp

      ! update with external data first
      if ( associated(laidata) )
     &     call entcell_update_lai(entcell, laidata)

      if ( associated(hdata) )
     &     call entcell_update_height(entcell, hdata)

      if ( associated(albedodata) )
     &     call entcell_update_albedo(entcell, albedodata)

      if ( associated(cropsdata) )
     &     call entcell_update_crops(entcell, cropsdata)

      ! and then do GISS phenology if required
      if ( do_giss_phenology ) then
        if ( hemi<-2 .or. jday <-2 )
     &       call stop_model("entcell_vegupdate: needs hemi,jday",255)
        pp => entcell%oldest
        do while (ASSOCIATED(pp))
          !* LAI, ALBEDO *!
          call prescr_phenology(jday,hemi, pp, do_giss_lai)
          call summarize_patch(pp)
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
      !* Calculate new LAI and albedo for given jday, for prescr vegetation. *!
      !* TBA:  THIS ROUTINE WILL ALSO UPDATE LIVE BIOMASS POOLS.           *!
      use ent_pfts
      use ent_prescr_veg, only : prescr_calc_lai,prescr_plant_cpools,
     &     prescr_veg_albedo
      use phenology, only : litter
      implicit none
      integer,intent(in) :: jday !Day of year.
      integer,intent(in) :: hemi !@var hemi -1: S.hemisphere, 1: N.hemisphere
      type(patch),pointer :: pp
      logical, intent(in) :: do_giss_lai
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: laipatch
      real*8 :: cpool(N_BPOOLS)
      real*8 :: lai_old

!      write(779,*) __FILE__,__LINE__
!     &   ,hemi,present(laidata)

      if (ASSOCIATED(pp)) then

        !* LAI *AND* BIOMASS - carbon pools *!
        laipatch = 0.d0         !Initialize for summing
        cpool(:) = 0.d0
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          if ( do_giss_lai ) then
            lai_old = cop%LAI
            cop%LAI = prescr_calc_lai(cop%pft+COVEROFFSET, jday, hemi)
            !* Calculate senescence factor for next time step litterfall routine.
            if (cop%LAI.gt.0d0) then !Prescribed senescence fraction.
              cop%senescefrac = max(0.d0,(lai_old-cop%LAI)/lai_old)
            else
              cop%senescefrac = 0.d0
            endif
          endif

          laipatch = laipatch + cop%LAI

          call prescr_plant_cpools(cop%pft, cop%lai, cop%h, 
     &         cop%dbh, cop%n, cpool )
          cop%C_fol = cpool(FOL)
          cop%C_sw = cpool(SW)
          cop%C_hw = cpool(HW)
!          cop%C_lab = cpool(LABILE)
          cop%C_froot = cpool(FR)
          cop%C_croot = cpool(CR)

          cop%Ntot = cop%nm * cop%LAI  !This should eventually go into N allocation routine if dynamic nm.

          cop => cop%shorter
        end do
        call litter(pp) !Daily time step
        pp%LAI = laipatch

        !* ALBEDO *!
        call prescr_veg_albedo(hemi, pp%tallest%pft, 
     &       jday, pp%albedo)

      endif
      end subroutine prescr_phenology


      end module ent_prescribed_updates
