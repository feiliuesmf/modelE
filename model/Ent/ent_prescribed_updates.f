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
      !use ent_GISSveg, only : GISS_plant_cpools
      type(entcelltype) :: ecp
      real*8,intent(in) :: laidata(N_PFT) !@var LAI for all PFT's 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 laipatch
      real*8 :: cpool(N_BPOOLS)

      pp => ecp%oldest      
      do while ( associated(pp) )
        laipatch = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
          cop%lai = laidata(cop%pft)
          laipatch = laipatch + cop%lai
          !!! this is hack, but don't know what to do with it at the moment...
cddd          call GISS_plant_cpools(cop%pft, cop%lai, cop%h, 
cddd     &         cop%dbh, cop%n, cpool )
cddd          cop%C_fol = cpool(FOL)
cddd          cop%C_sw = cpool(SW)
cddd          cop%C_hw = cpool(HW)
cddd          cop%C_lab = cpool(LABILE)
cddd          cop%C_froot = cpool(FR)
cddd          cop%C_croot = cpool(CR)

          cop => cop%shorter
        enddo
        pp%sumcohort%LAI = laipatch
        pp => pp%younger
      enddo

      end subroutine entcell_update_lai
 

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
          if ( pp%sumcohort%pft > N_PFT .or.  pp%sumcohort%pft < 1 )
     &         call stop_model("entcell_update_albedo: bad pft",255)
          pp%albedo(1:N_BANDS) = albedodata(1:N_BANDS, pp%sumcohort%pft)
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
      use ent_GISSveg, only : GISS_calc_shc
      use entcells, only : entcell_extract_pfts
      type(entcelltype) :: ecp
      !-----Local---------
      real*8 vdata(N_COVERTYPES) ! needed for a hack to compute canopy

      vdata(:) = 0.d0
      call entcell_extract_pfts( ecp, vdata(2:) )
      ecp%heat_capacity=GISS_calc_shc(vdata)

      end subroutine entcell_update_shc
     

      subroutine entcell_vegupdate(entcell, hemi, jday
     &     ,do_giss_phenology
     &     ,laidata, albedodata, cropsdata )
!@sum updates corresponding data on entcell level (and down). 
!@+   everything except entcell is optional. coordinate-dependent
!@+   is given pointer attribute to provide a way to tell the 
!@+   program that an argument is actually optional and missing
!@+   (see how it is used in ent_prescribe_vegupdate)
      use patches, only : summarize_patch
      use entcells,only : summarize_entcell!,entcell_extract_pfts
      use ent_GISSveg, only : GISS_calc_shc

      type(entcelltype) :: entcell
      integer,intent(in), optional :: jday
      integer, intent(in), pointer, optional :: hemi
      logical, intent(in), optional :: do_giss_phenology
      real*8, intent(in), pointer, optional :: laidata(:)
      real*8, intent(in), pointer, optional :: albedodata(:,:)
      real*8, intent(in), pointer, optional :: cropsdata
      !----Local------
      type(patch),pointer :: pp
      type(cohort),pointer :: cop

      if ( present(do_giss_phenology) ) then
        if ( do_giss_phenology ) then
          if ( .not. (present(hemi).and.present(jday)) )
     &         call stop_model("entcell_vegupdate: needs hemi,jday",255)
          if ( .not. associated(hemi) )
     &         call stop_model("entcell_vegupdate: needs hemi",255)
          pp => entcell%oldest
          do while (ASSOCIATED(pp))
          !* LAI, ALBEDO *!
            call GISS_phenology(jday,hemi, pp)
            call summarize_patch(pp)
            pp => pp%younger
          end do
        endif
      endif

      if ( present(laidata) ) then
        if ( associated(laidata) )
     &       call entcell_update_lai(entcell, laidata)
      endif

      if ( present(albedodata) ) then
        if ( associated(albedodata) )
     &       call entcell_update_albedo(entcell, albedodata)
      endif

      if ( present(cropsdata) ) then
        if ( associated(cropsdata) )
     &       call entcell_update_crops(entcell, cropsdata)
      endif

      call summarize_entcell(entcell)

      call entcell_update_shc(entcell)
cddd      vdata(:) = 0.d0
cddd      call entcell_extract_pfts(entcell, vdata(2:) )
cddd      entcell%heat_capacity=GISS_calc_shc(vdata)

      !if (YEAR_FLAG.eq.0) call ent_GISS_init(entcellarray,im,jm,jday,year)
      !!!### REORGANIZE WTIH ent_prog.f ####!!!
      
      end subroutine entcell_vegupdate


      subroutine GISS_phenology(jday,hemi, pp)
      !* Calculate new LAI and albedo for given jday, for GISS vegetation. *!
      !* TBA:  THIS ROUTINE WILL ALSO UPDATE LIVE BIOMASS POOLS.           *!
      use ent_pfts
      use ent_GISSveg, only : GISS_calc_lai,GISS_plant_cpools,
     &     GISS_veg_albedo
      implicit none
      integer,intent(in) :: jday !Day of year.
      integer,intent(in) :: hemi !@var hemi -1: S.hemisphere, 1: N.hemisphere
      type(patch),pointer :: pp
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: laipatch
      real*8 :: cpool(N_BPOOLS)

      if (ASSOCIATED(pp)) then

        !* LAI *AND* BIOMASS - carbon pools *!
        laipatch = 0.d0 !Initialize for summing
        cpool(:) = 0.d0
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          cop%lai = GISS_calc_lai(cop%pft+COVEROFFSET, jday, hemi)
          laipatch = laipatch + cop%lai

          call GISS_plant_cpools(cop%pft, cop%lai, cop%h, 
     &         cop%dbh, cop%n, cpool )
          cop%C_fol = cpool(FOL)
          cop%C_sw = cpool(SW)
          cop%C_hw = cpool(HW)
          cop%C_lab = cpool(LABILE)
          cop%C_froot = cpool(FR)
          cop%C_croot = cpool(CR)

          cop => cop%shorter
        end do
        pp%sumcohort%LAI = laipatch

        !* ALBEDO *!
        call GISS_veg_albedo(hemi, pp%sumcohort%pft, 
     &       jday, pp%albedo)

      endif
      end subroutine GISS_phenology


      end module ent_prescribed_updates
