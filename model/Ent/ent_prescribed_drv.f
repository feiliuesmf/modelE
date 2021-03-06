#include "rundeck_opts.h"

      module ent_prescribed_drv

      use ent_const
      use ent_pfts
      use ent_prescr_veg

      implicit none
      private
      save


      public 
     &     init_canopy_physical,
     &     prescr_vegdata,
     &     init_vdata, prescr_get_laidata, 
     &     prescr_update_vegcrops, prescr_veg_albedodata,
     &     prescr_get_height, !YKIM
     &     prescr_get_soilpools,  !for prescribing soil C, N pools -PK 12/07
     &     prescr_get_cropdata,
     &     prescr_get_soil_C_total

      public prescr_get_hdata,prescr_get_woodydiameter,prescr_get_pop,
     &     prescr_get_crownrad,prescr_get_carbonplant,prescr_get_initnm,
     &     prescr_get_rootprof,prescr_get_soilcolor

      public init_ent_laidata, init_ent_hdata,  prescr_get_ent_plant

      contains

      !*********************************************************************
      !*    SUBROUTINES TO READ IN prescr VEGETATION DATA SETS 
      !*********************************************************************

!***************************************************************************
      subroutine init_canopy_physical(
     & I0,I1,J0,J1,Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
      integer,intent(in) :: I0,I1,J0,J1
      real*8, DIMENSION(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini

      Ci_ini(:,:) = 0.0127d0
      CNC_ini(:,:) = 0.d0
      Tcan_ini(:,:) = 0.d0           !Should be a forcing from land surface model.
      Qf_ini(:,:) = 0.d0             !Should be a forcing from land surface model.

      end subroutine init_canopy_physical
      
!***************************************************************************
      subroutine prescr_get_soilpools(I0,I1,J0,J1,
     &     soil_C_total, Tpool_ini)
      !this routine reads in total soil pool amounts (measured), 
      !and individual soil pool fractions (modeled pft-dependent values from spinup runs),
      !then prescribes individual amounts **all carbon amounts should be in g/m2** -PK 12/07   
      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in) ::
     &     soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      real*8,intent(out) :: 
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,  !prescribed soil pools, g/m2
     &                I0:I1,J0:J1)
      !-----Local------
!      first 3 for eventually reading in globally gridded dataset, e.g. ISRIC-WISE
      integer :: iu_SOILCARB
      integer :: n,p,nn
      real*8, dimension(N_CASA_LAYERS) :: total_Cpool  !site-specific total measured soil C_org
      real*8, dimension(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs  !modeled soil C_org pool fractions
      real*8, dimension(NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_tmp !YK

      Tpool_ini(:,:,:,:,:,:) = 0.d0  !initialize all pools to zero

#ifdef PFT_MODEL_ENT
!YK - temp. values, modified from 8 GISS pfts below
        Cpool_fracs(1,:,1) = (/ !ever_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(2,:,1) = (/ !ever_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(3,:,1) = (/ !ever_ES_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(4,:,1) = (/ !ever_LS_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(5,:,1) = (/ !cold_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(6,:,1) = (/ !cold_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(7,:,1) = (/ !drought_broad
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(8,:,1) = (/ !decid_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(9,:,1) = (/ !shrub_cold
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(10,:,1) = (/ !shrub_arid
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(11,:,1) = (/ !c3grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(12,:,1) = (/ !c4grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(13,:,1) = (/ !c3grass_ann
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(14,:,1) = (/ !c3grass_arctic
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(15,:,1) = (/ !cropsc4
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(16,:,1) = (/ !cropstree
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
#else
!***  for now define for 8 GISS pfts, one-layer only -PK 1/23/08***
        Cpool_fracs(1,:,1) = (/ !tundra (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(2,:,1) = (/ !C3 grass (Vaira)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(3,:,1) = (/ !shrub (for now=savanna)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(4,:,1) = (/ !savanna (Tonzi)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(5,:,1) = (/ !decid broadl (MMSF)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(6,:,1) = (/ !evergr needl (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(7,:,1) = (/ !trop rainf (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(8,:,1) = (/ !crops (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
#endif

#ifdef SOILCARB_SITE
!####temporary hack (for site-specific runs)####
!for now, external file should be named as below and should be organized as follows:
!(1) there should be 1 or 2 columns (corresponding to each soil bgc layer);
!(2) first non-header row should have total site-measured pool (in g/m2);
!(3) 9 subsequent rows correspond to modeled 9 soil pool fractions
!###YK- hack but better way...at least, do not need to change according to the sites.
      call openunit("SOILCARB_site",iu_SOILCARB,.false.,.true.)  !formatted dataset
      read(iu_SOILCARB,*)  !skip optional header row(s)
      read(iu_SOILCARB,*) total_Cpool(:)
      do nn=1,NPOOLS-NLIVE
        read(iu_SOILCARB,*) Cpool_tmp(nn,:)
      end do

      do p=1,N_PFT      
       do n=1,N_CASA_LAYERS 
        do nn=NLIVE+1,NPOOLS
         Cpool_fracs(p,nn-NLIVE,n) = Cpool_tmp(nn-NLIVE,n)
         Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &          Cpool_fracs(p,nn-NLIVE,n)*total_Cpool(n)
        end do
       end do
      end do
#else
        !assign Tpool_ini values (pft-specific)
        do p=1,N_PFT
          do n=1,N_CASA_LAYERS 
            do nn=NLIVE+1,NPOOLS
              Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &             Cpool_fracs(p,nn-NLIVE,n)
     &             * soil_C_total(n,I0:I1,J0:J1)*1d3  
            end do
          end do
        end do
#endif

      end subroutine prescr_get_soilpools


      subroutine prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,
     &     soil_C_total)
      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      real*8,intent(out) ::
     &     soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      !---
      real*4 :: buf(N_CASA_LAYERS,IM,JM)
      character*80 :: title
      integer :: iu_SOILCARB

      call openunit("SOILCARB_global",iu_SOILCARB,.true.,.true.)
      read (iu_SOILCARB) title, buf !data in kg/m2 (converted to g/m2 below)
      call closeunit(iu_SOILCARB)

      soil_C_total(1:N_CASA_LAYERS,I0:I1,J0:J1) =
     &     buf(1:N_CASA_LAYERS,I0:I1,J0:J1)

      end subroutine prescr_get_soil_C_total

      
!***************************************************************************
      subroutine prescr_vegdata(jday, year, IM,JM,I0,I1,J0,J1,
     &     vegdata,albedodata,laidata,hdata,nmdata,popdata,dbhdata,
     &     craddata,cpooldata,rootprofdata,soil_color,soil_texture,
     &     Tpooldata, 
     &     do_soilinit,do_phenology_activegrowth,do_read_from_files)
      implicit none
      integer,intent(in) :: jday, year
      integer,intent(in) :: IM,JM,I0,I1,J0,J1 !long/lat grid number range
      real*8,intent(out) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: laidata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: nmdata(N_COVERTYPES)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)
      real*8,intent(out) :: popdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: dbhdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: craddata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      integer,intent(out) :: soil_color(N_COVERTYPES)
      real*8,intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpooldata !in g/m2 -PK
      logical,intent(in) :: do_soilinit
      logical,intent(in) :: do_phenology_activegrowth
      logical,intent(in) :: do_read_from_files

      !-----Local------
      integer :: i,j, jeq
      integer hemi(I0:I1,J0:J1)
      REAL*8 :: soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)

      jeq = JM/2
      do j=J0,J1
        hemi(:,j) = 1
        if (j <= jeq) hemi(:,j) = -1
      enddo


!YKIM
cddd      call init_vdata(IM,JM,I0,I1,J0,J1,vegdata)   !veg fractions
cddd      call prescr_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)
cddd      call prescr_get_laidata(jday,JM,I0,I1,J0,J1,laidata) !lai
cddd      call prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1,vegdata)
cddd      call prescr_get_hdata(hdata) !height
cddd      !print *,'hdata',hdata
cddd      call prescr_get_initnm(nmdata) !nm
cddd      call prescr_get_rootprof(rootprofdata)
cddd      call prescr_get_woodydiameter(hdata,dbhdata)
cddd      !print *,"dbhdata",dbhdata
cddd      call prescr_get_pop(dbhdata,popdata)
cddd      !print *,"popdata",popdata
cddd      call prescr_get_crownrad(popdata,craddata)
cddd      call prescr_get_carbonplant(IM,JM,I0,I1,J0,J1,
cddd     &     laidata,hdata,dbhdata,popdata,cpooldata)
cddd      !do i=1,N_COVERTYPES
cddd      !  print*,"cpooldata(ncov)",i,cpooldata(i,:,I0:I1,J0:J1)
cddd      !end do
cddd      call prescr_get_soilcolor(soil_color)
cddd      call prescr_get_soiltexture(IM,JM,I0,I1,J0,J1,
cddd     &     soil_texture)
cddd      call prescr_soilpools(IM,JM,I0,I1,J0,J1,Tpooldata,do_soilinit)

!YKIM
!change sequences of calls
!to have options according to do_phenology_activegrowth

      if ( do_read_from_files )
     &     call init_vdata(IM,JM,I0,I1,J0,J1,vegdata)   !veg fractions
      if ( do_read_from_files )
     &     call prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1,vegdata)
      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)
      if (.not.do_phenology_activegrowth) then
         call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata) !lai
         do j=J0,J1
           do i=I0,I1
             call prescr_get_hdata(hdata(:,i,j)) !height
             call prescr_get_woodydiameter(hdata(:,i,j), dbhdata(:,i,j))
             call prescr_get_pop(dbhdata(:,i,j), popdata(:,i,j))
             call prescr_get_crownrad(popdata(:,i,j), craddata(:,i,j))
           enddo
         enddo
         call prescr_get_carbonplant(I0,I1,J0,J1,
     &        laidata,hdata,dbhdata,popdata,cpooldata)
      else !if do_phenology_activegrowth=true
         call init_ent_laidata(IM,JM,I0,I1,J0,J1,laidata) !lai
         call init_ent_hdata(IM,JM,I0,I1,J0,J1,hdata) !height
         !update diameter, population density, carbon plant &  crown rad
         !can be more modular like the above - Should I???  -YKIM
         call prescr_get_ent_plant(I0,I1,J0,J1, 
     &        laidata,hdata,dbhdata,popdata,craddata,cpooldata)
      endif
      call prescr_get_initnm(nmdata) !nm ! mean canopy nitrogen
      call prescr_get_rootprof(rootprofdata)
      call prescr_get_soilcolor(soil_color)
      if ( do_read_from_files )
     &     call prescr_get_soiltexture(IM,JM,I0,I1,J0,J1,
     &     soil_texture)
#ifdef SET_SOILCARBON_GLOBAL_TO_ZERO
      Tpooldata = 0.d0
#else
      if ( do_soilinit ) then
#ifndef SOILCARB_SITE
        call prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,soil_C_total)
#endif
        call prescr_get_soilpools(I0,I1,J0,J1,soil_C_total,Tpooldata)
      else
        Tpooldata = 0.d0
      endif
#endif
      !print*,'vegdata(:,I1,J1)',vegdata(:,I1,J1)
      !print*,'hdata',hdata
      !print*,'nmdata',nmdata
      !print*,'dbhdata',dbhdata
      !print*,'popdata',popdata
      !print*,'craddata',craddata
!      print*,'cpooldata(GRASSC3+COVEROFFSET,:,I1,J1)',
!     &     cpooldata(GRASSC3+COVEROFFSET,:,I1,J1)
!      print*,'cpooldata(SHRUB+COVEROFFSET,:,I1,J1)',
!     &     cpooldata(SHRUB+COVEROFFSET,:,I1,J1)
!      print*,'cpooldata(SAVANNA+COVEROFFSET,:,I1,J1)',
!     &     cpooldata(SAVANNA+COVEROFFSET,:,I1,J1)
      !print*,'soil_color',soil_color

      end subroutine prescr_vegdata


!***************************************************************************
      subroutine init_vdata(im,jm,I0,I1,J0,J1,vdata)
      !* This version reads in vegetation structure from prescr data set.
      use FILEMANAGER, only : openunit,closeunit,nameunit
      integer, intent(in) :: im,jm,I0,I1,J0,J1
      real*8, intent(out) :: vdata(N_COVERTYPES,I0:I1,J0:J1) 
      !------Local---------------------
      !1    2    3    4    5    6    7    8    9   10   11    12
      !BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE GRAC4
      character*80 :: title
      real*4 :: buf(im,jm)
      integer :: iu_VEG
      integer :: k,i,j
      real*8 :: s

      ! Make sure that unused fractions are set to 0
      vdata(:,:,:) = 0.d0
      call openunit("VEG",iu_VEG,.true.,.true.)

      do k=1,N_COVERTYPES-N_OTHER !## Skip algae and grac4 #HACK
        !print *,k
        read(iu_VEG) title , buf
        vdata(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        !print *,"read VEG:", title
      end do
      !print *,"vdata", vdata(:,I0:I1,J0:J1) !#DEBUG
      call closeunit(iu_VEG)

      ! make sure that veg fractions are reasonable
      do j=J0,J1
        do i=I0,I1
          do k=1,N_COVERTYPES
            ! get rid of unreasonably small fractions
            if ( vdata(k,i,j) < 1.d-4 ) vdata(k,i,j) = 0.d0
          enddo
          s = sum( vdata(:,i,j) )
          if ( s > .9d0 ) then
            vdata(:,i,j) = vdata(:,i,j)/s
          else if ( s < .1d0 ) then
            print *, "missing veg data at ",i,j,"assume bare soil"
            vdata(:,i,j) = 0.d0
            vdata(COVER_SAND,i,j) = 1.d0
          else
            call stop_model("Incorrect data in VEG file",255)
          endif
        enddo
      enddo
          
      end subroutine init_vdata

!**************************************************************************
      subroutine prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata)
      !* This version reads in crop distribution from prescr data set.
      !* And calculates crop fraction for given year.
      use FILEMANAGER, only : openunit,closeunit,nameunit
      integer, intent(in) :: year
      integer, intent(in) :: IM, JM, I0, I1, J0, J1
      real*8, intent(out) :: cropdata(I0:I1,J0:J1)
      integer i
      !----------
      integer :: iu_CROPS
      integer :: year1, year2
      real*4 crop4(im,jm)
      real*8 wt, crop1(I0:I1,J0:J1), crop2(I0:I1,J0:J1)
      character*80 title

      !* Calculate fraction for given gcmtime:  interpolate between years*/

      year1 = -32768 ; crop1(:,:) = 0.d0
      year2 = -32767 ; crop2(:,:) = 0.d0
      wt = 1.d0

      call openunit("CROPS",iu_CROPS,.true.,.true.)
      do while( year2 < year )
        year1 = year2
        crop1(:,:) = crop2(:,:)
        read (iu_CROPS,end=10) title , crop4
        read(title,*) year2 !Read year integer out of character array title
        crop2(I0:I1,J0:J1) = crop4(I0:I1,J0:J1)
      enddo
      wt = (year-year1)/(real(year2-year1,kind=8))
 10   continue
      call closeunit(iu_CROPS)

      cropdata(:,:) = max(0.d0, crop1(:,:)
     &     + wt * (crop2(:,:) - crop1(:,:)))  !Set min to zero, since no land mask yet -nyk 1/22/08

      end subroutine prescr_get_cropdata

!**************************************************************************

      subroutine prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1,
     &     vegdata)
      !* Modify vegdata given new cropdata.
      integer,intent(in) :: year
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8, intent(inout) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      !--------
      real*8,ALLOCATABLE,dimension(:,:) :: cropdata !grid array
      integer :: i,j
      real*8 crops_old

      ALLOCATE(cropdata(I0:I1,J0:J1))

      !* Loop *!
!!! hack --> remove crops
!!!      cropdata(:,:) = 0.d0
      call prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata) !crop fraction

      !* If cropdata was prepared somewhere else, then cover is as simple as
      !* modifying the vegetation fractions.  Need to update cohort and
      !* patch summary variables.
      do j=J0,J1
        do i=I0,I1
          if ( cropdata(i,j) == 1.d0 ) then
            vegdata(:,i,j) = 0.d0
            vegdata(CROPS+COVEROFFSET,i,j) = 1.d0
          else
            crops_old = vegdata(CROPS+COVEROFFSET,i,j)
            if ( crops_old == 1.d0 ) then
              call stop_model("incompatible crops: old=1, new<1",255)
            endif
            vegdata(:,i,j) = vegdata(:,i,j)
     $           * (1.d0-cropdata(i,j))/(1.d0-crops_old)
            vegdata(CROPS+COVEROFFSET,i,j) = cropdata(i,j)
          endif
        end do
        !write(*,*) vegdata(CROPS+COVEROFFSET,:,j)
      end do

      DEALLOCATE(cropdata)
 
      end subroutine prescr_update_vegcrops

!**************************************************************************

      subroutine prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
!@sum Returns prescr GCM leaf area index for entire grid and given jday.
      use FILEMANAGER, only : openunit,closeunit,nameunit !for VEG_PROGNOSTIC
      use ent_const,only : N_COVERTYPES
      integer,intent(in) :: jday
      integer, intent(in) :: I0,I1,J0,J1
      real*8 :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      integer :: hemi(I0:I1,J0:J1) !@var hemi =1 in N. hemisphere, =-1 in South
      !----------
      integer :: n !@var cover type
      integer i,j,jeq

      !jeq = JM/2

      do j=J0,J1
        !hemi = 1
        !if (j <= jeq) hemi = -1
        do i=I0,I1
          do n=1,N_COVERTYPES
            laidata(n,i,j) = prescr_calc_lai(n,jday,hemi(i,j))
          enddo
        enddo
      enddo

      !* Return lai for each vegetation type.
      end subroutine prescr_get_laidata

!**************************************************************************

      subroutine init_ent_laidata(IM,JM,I0,I1,J0,J1,laidata)
!@sum YKIM - read the initial LAI for the prognostic vegetation 
!@sum (i.e., prognostic phenology/growth with Ent PFTs)

      use FILEMANAGER, only : openunit,closeunit,nameunit !for VEG_PROGNOSTIC
      use ent_const,only : N_COVERTYPES
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8 :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      !----------
      character*80 :: title
      real*4 :: buf(IM,JM)
      integer :: iu_LAI
      integer :: k !@var cover type

      call openunit("LAIent",iu_LAI,.true.,.true.)

      do k=1,N_COVERTYPES
        read(iu_LAI) title , buf
        laidata(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      end do
      !print *,"laidata", laidata(:,I0:I1,J0:J1) !#DEBUG

      call closeunit(iu_LAI)

      !* Return lai for each vegetation type.
      end subroutine init_ent_laidata


!**************************************************************************
      subroutine prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)
      integer,intent(in) :: jday
      integer, intent(in) :: I0,I1,J0,J1
      real*8 :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      integer :: hemi(I0:I1,J0:J1)
      !----------
      integer :: pft
      integer i,j,jeq
      
      !jeq = JM/2

      do j=J0,J1
        !hemi = 1
        !if (j <= jeq) hemi = -1
        do i=I0,I1
          do pft = 1, N_COVERTYPES
            call prescr_veg_albedo(hemi(i,j),pft,jday,
     &           albedodata(:,pft,i,j))
          end do
        enddo
      enddo

      end subroutine prescr_veg_albedodata

!**************************************************************************

      subroutine prescr_get_height(hdata)
!@sum Returns prescr GCM leaf area index for entire grid and given jday.
      use ent_const,only : N_COVERTYPES
      real*8 :: hdata(N_COVERTYPES) 
      !----------
      
      call prescr_get_hdata(hdata)


      end subroutine prescr_get_height

!**************************************************************************

      subroutine init_ent_hdata(IM,JM,I0,I1,J0,J1,hdata3d)
!@sum YKIM - read the initial height for the prognostic vegetation 
!@sum (i.e., prognostic phenology/growth with Ent PFTs)
      use FILEMANAGER, only : openunit,closeunit,nameunit !for VEG_PROGNOSTIC
      use ent_const,only : N_COVERTYPES
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8 :: hdata3d(N_COVERTYPES,I0:I1,J0:J1) 
      !----------
      character*80 :: title
      real*4 :: buf(IM,JM)
      integer :: iu_HITE
      integer :: k !@var cover type
     
      call openunit("HITEent",iu_HITE,.true.,.true.)

      do k=1,N_COVERTYPES
        read(iu_HITE) title , buf
        hdata3d(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      end do
      !print *,"hdata", hdata(:,I0:I1,J0:J1) !#DEBUG

      call closeunit(iu_HITE)

      end subroutine init_ent_hdata

!**************************************************************************
      subroutine prescr_get_carbonplant(I0,I1,J0,J1,
     &     laidata, hdata, dbhdata, popdata, cpooldata)
      !*  Calculate per plant carbon pools (g-C/plant).
      !*  After Moorcroft, et al. (2001).
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in) :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(in) :: hdata(N_COVERTYPES)
      real*8,intent(in) :: dbhdata(N_COVERTYPES)
      real*8,intent(in) :: popdata(N_COVERTYPES) 
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      !Array: 1-foliage, 2-sapwood, 3-hardwood, 4-labile,5-fine root, 6-coarse root
      !----Local----
      integer :: p,pft !@var pft vegetation type
      integer i,j,jeq

      cpooldata(:,:,:,:) = 0.d0  !Zero initialize
      do j=J0,J1
        do i=I0,I1
          do pft=1,N_PFT
            p = pft + COVEROFFSET
            call prescr_plant_cpools(pft,laidata(p,i,j),hdata(p),
     &           dbhdata(p), popdata(p),cpooldata(p,:,i,j))
            call prescr_init_Clab(pft,popdata(p),cpooldata(p,:,i,j))
          enddo
        enddo
      enddo
      
      end subroutine prescr_get_carbonplant
!**************************************************************************
      subroutine prescr_get_ent_plant(I0,I1,J0,J1,
     &     laidata, hdata3d, dbhdata3d, popdata3d, craddata3d,cpooldata)
!@sum YKIM- calculate woody diameter, population denisty, crown radiation
!@sum & carbon pools for the Ent prognostic vegetation
      use phenology, only: update_plant_cpools, height2dbh,nplant
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in) :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(in) :: hdata3d(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: dbhdata3d(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: popdata3d(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(out) :: craddata3d(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      !Array: 1-foliage, 2-sapwood, 3-hardwood, 4-labile,5-fine root, 6-coarse root
      !----Local----
      integer :: pft !@var pft vegetation type
      integer :: i,j, n

!Zero initialize.
      dbhdata3d(:,:,:) = 0.0 
      popdata3d(:,:,:) = 0.0 
      craddata3d(:,:,:) = 0.0
      cpooldata(:,:,:,:) = 0.d0 
      do j=J0,J1
        do i=I0,I1
          do pft = 1,N_PFT
            n = pft + COVEROFFSET
            !update woody diameter
            if (pfpar(pft)%woody)  dbhdata3d(n,i,j) 
     &                             =height2dbh(pft,hdata3d(n,i,j))  
            !update population denisty
            popdata3d(n,i,j) = nplant(pft,dbhdata3d(n,i,j), 
     &                       hdata3d(n,i,j), laidata(n,i,j))
            !update crown rad
            craddata3d(n,i,j) = 0.5*sqrt(1/popdata3d(n,i,j))
            !update cabon pools
            call update_plant_cpools(pft,laidata(n,i,j),hdata3d(n,i,j),
     &           dbhdata3d(n,i,j), popdata3d(n,i,j),cpooldata(n,:,i,j))
            call prescr_init_Clab(pft,popdata3d(n,i,j),
     &           cpooldata(n,:,i,j))
          enddo
        enddo
      enddo
      
      end subroutine prescr_get_ent_plant

!*************************************************************************
      subroutine prescr_get_soiltexture(im,jm,I0,I1,J0,J1,
     &     soil_texture)
      !* Return arrays of GISS soil color and texture.
      use FILEMANAGER, only : openunit,closeunit
      integer, intent(in) :: im,jm,I0,I1,J0,J1
      real*8, intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      !------
      real*8 :: buf(im,jm,N_SOIL_TEXTURES)
      integer :: iu_SOIL
      integer k

      call openunit("soil_textures",iu_SOIL,.true.,.true.)
      !print *,IM,JM,N_COVERTYPES !#DEBUG
      read(iu_SOIL) buf
      call closeunit(iu_SOIL)
      !print *,"soil fractions:",buf(I0,J0,:)!#DEBUG

      do k=1,N_SOIL_TEXTURES
        soil_texture(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1,k)
      enddo

!      do j=J0,J1
!        do i=I0,I1
!          print *,'from GISS_get_soiltypes (in ent_GISSveg):' 
!          print *, soil_texture(:,i,j) , sum(soil_texture(:,i,j))
!          print *,'soiltexture=',soil_texture(:,i,j)
!     &           ,'sum(soiltextures)=',sum(soil_texture(:,i,j))
!        enddo
!      enddo
      end subroutine prescr_get_soiltexture
!*************************************************************************
      end module ent_prescribed_drv

