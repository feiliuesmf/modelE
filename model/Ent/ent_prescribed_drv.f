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
     &     prescr_get_vdata, prescr_get_laidata, 
     &     prescr_update_vegcrops, prescr_veg_albedodata,
     &     prescr_soilpools  !for prescribing soil C, N pools -PK 12/07

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
      subroutine prescr_soilpools(IM,JM,I0,I1,J0,J1,Tpool_ini)
      !this routine reads in total soil pool amounts (measured), 
      !and individual soil pool fractions (modeled pft-dependent values from spinup runs),
      !then prescribes individual amounts **all amounts in g/m2** -PK 12/07   
      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      real*8,intent(out) :: 
     &      Tpool_ini(PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,I0:I1,J0:J1)  !prescribed soil pools
      !-----Local------
!      first 3 for eventually reading in globally gridded dataset, e.g. ISRIC-WISE
!      real*4 :: soilC_data(N_CASA_LAYERS,IM,JM)
      integer :: iu_SOILCARB  !and pft? (for reading in pool fractions)
!      character*80 :: title
      integer :: i,n
      real*8, dimension(N_CASA_LAYERS) :: total_Cpool  !measured total C_org pool
      real*8, dimension(N_CASA_LAYERS,NPOOLS-NLIVE) :: Cpool_fracs  !modeled soil C_org pool fractions

      Tpool_ini(:,:,:,:,:) = 0.d0  !initialize all pools to zero
      
!***for eventual reading of global data***      
!      call openunit("SOILCARB_global",iu_SOILCARB,.true.,.true.)  !globally gridded binary dataset
!      read (iu_SOILCARB) title, soilC_data
!      .....add code to also read in pft-specific pool fractions (e.g., array of arrays?) 
!      do n=1,N_CASA_LAYERS 
!       do i=NLIVE+1,NPOOLS
!        Tpool_ini(CARBON,i,n,...,...) = Cpool_fracs(n,i-NLIVE,pft?)*soilC_data(n,...,...)
!       end do
!      end do

!####temporary hack (for site-specific runs)####
!for now, external file should be named as below and should be organized as follows:
!(1) there should be 1 or 2 columns (corresponding to each soil bgc layer);
!(2) first non-header row should have total site-measured pool (in g/m2);
!(3) 9 subsequent rows correspond to modeled 9 soil pool fractions
      call openunit("SOILCARB_site",iu_SOILCARB,.false.,.true.)  !unformatted dataset
      read(iu_SOILCARB,*)  !skip optional header row(s)
      read(iu_SOILCARB,*) total_Cpool(:)
      do i=1,NPOOLS-NLIVE
        read(iu_SOILCARB,*) Cpool_fracs(:,i)
      end do
!####
      
      do n=1,N_CASA_LAYERS 
       do i=NLIVE+1,NPOOLS
        Tpool_ini(CARBON,i-NLIVE,n,I0:I1,J0:J1) =
     &          Cpool_fracs(n,i-NLIVE)*total_Cpool(n)
       end do
      end do

      call closeunit(iu_SOILCARB)

      end subroutine prescr_soilpools
      
!***************************************************************************
      subroutine prescr_vegdata(jday, year, IM,JM,I0,I1,J0,J1,
     &     vegdata,albedodata,laidata,hdata,nmdata,popdata,dbhdata,
     &     craddata,cpooldata,rootprofdata,soil_color,soil_texture)
      integer,intent(in) :: jday, year
      integer,intent(in) :: IM,JM,I0,I1,J0,J1 !long/lat grid number range
      real*8,intent(out) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: laidata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: hdata(N_COVERTYPES)
      real*8,intent(out) :: nmdata(N_COVERTYPES)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)
      real*8,intent(out) :: popdata(N_COVERTYPES)
      real*8,intent(out) :: dbhdata(N_COVERTYPES)
      real*8,intent(out) :: craddata(N_COVERTYPES)
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      integer,intent(out) :: soil_color(N_COVERTYPES)
      real*8,intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      !-----Local------

      call prescr_get_vdata(IM,JM,I0,I1,J0,J1,vegdata)   !veg fractions
      call prescr_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)
      call prescr_get_laidata(jday,JM,I0,I1,J0,J1,laidata) !lai
      call prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1,vegdata)

      call prescr_get_hdata(hdata) !height
      call prescr_get_initnm(nmdata) !nm
      call prescr_get_rootprof(rootprofdata)
      call prescr_get_woodydiameter(hdata,dbhdata)
      call prescr_get_pop(dbhdata,popdata)
      call prescr_get_crownrad(popdata,craddata)
      call prescr_get_carbonplant(IM,JM,I0,I1,J0,J1,
     &     laidata,hdata,dbhdata,popdata,cpooldata)
      call prescr_get_soilcolor(soil_color)
      call prescr_get_soiltexture(IM,JM,I0,I1,J0,J1,
     &     soil_texture)

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
      subroutine prescr_get_vdata(im,jm,I0,I1,J0,J1,vdata)
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
      integer :: k

      ! Make sure that unused fractions are set to 0
      vdata(:,:,:) = 0.d0
      call openunit("VEG",iu_VEG,.true.,.true.)

      do k=1,N_COVERTYPES-2  !## Skip algae and grac4 #HACK
        read(iu_VEG) title , buf
        vdata(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        print *,"read VEG:", title
      end do
      print *,"vdata", vdata(:,I0:I1,J0:J1) !#DEBUG
      call closeunit(iu_VEG)
      end subroutine prescr_get_vdata

!**************************************************************************
      subroutine prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata)
      !* This version reads in crop distribution from prescr data set.
      !* And calculates crop fraction for given year.
#define tempdebug
#ifdef tempdebug
      use FILEMANAGER, only : openunit,closeunit,nameunit
#endif
      integer, intent(in) :: year
      integer, intent(in) :: IM, JM, I0, I1, J0, J1
      real*8, intent(out) :: cropdata(I0:I1,J0:J1)
      !----------
#ifdef tempdebug
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
        print *, "Got here in prescr_get_cropdata"
        year1 = year2
        crop1(:,:) = crop2(:,:)
        read (iu_CROPS,end=10) title , crop4
        read(title,*) year2 !Read year integer out of character array title
        print *,"read CROPS:",title,year2
        crop2(I0:I1,J0:J1) = crop4(I0:I1,J0:J1)
      enddo
      wt = (year-year1)/(real(year2-year1,kind=8))
 10   continue
      call closeunit(iu_CROPS)

      cropdata(:,:) = crop1(:,:)
     &     + wt * (crop2(:,:) - crop1(:,:))
#else
      !*TEMPORARY ZERO OUT CROPDATA IFDEF *!
      do i=I0,I1
          cropdata(i,J0:J1) = 0.0
      end do
#endif
      !* Return cropdata single layer crop fraction.
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
      call prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata) !crop fraction

      !* If cropdata was prepared somewhere else, then cover is as simple as
      !* modifying the vegetation fractions.  Need to update cohort and
      !* patch summary variables.
      do j=J0,J1
        do i=I0,I1
          if ( cropdata(i,j) == 1.d0 ) then
            vegdata(:,i,j) = 0.d0
            vegdata(:,i,j) = 1.d0
          else
            crops_old = vegdata(9,i,j)
            if ( crops_old == 1.d0 ) then
              call stop_model("incompatible crops: old=1, new<1",255)
            endif
            vegdata(:,i,j) = vegdata(:,i,j)
     $           * (1.d0-cropdata(i,j))/(1.d0-crops_old)
            vegdata(9,i,j) = cropdata(i,j)
          endif
        end do
      end do

      DEALLOCATE(cropdata)
 
      end subroutine prescr_update_vegcrops

!**************************************************************************

      subroutine prescr_get_laidata(jday,JM,I0,I1,J0,J1,laidata)
!@sum Returns prescr GCM leaf area index for entire grid and given jday.
      use ent_const,only : N_COVERTYPES
      integer,intent(in) :: jday
      integer, intent(in) :: JM,I0,I1,J0,J1
      real*8 :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      !----------
      integer :: n !@var cover type
      integer :: hemi !@var hemi =1 in N. hemisphere, =-1 in South
      integer i,j,jeq

      jeq = JM/2

      do j=J0,J1
        hemi = 1
        if (j <= jeq) hemi = -1
        do i=I0,I1
          do n=1,N_COVERTYPES
            laidata(n,i,j) = prescr_calc_lai(n,jday,hemi)
          enddo
        enddo
      enddo

      !* Return lai for each vegetation type.
      end subroutine prescr_get_laidata


!**************************************************************************
      subroutine prescr_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)
      integer,intent(in) :: jday
      integer, intent(in) :: JM,I0,I1,J0,J1
      real*8 :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      !----------
      integer :: hemi, pft
      integer i,j,jeq
      
      jeq = JM/2

      do j=J0,J1
        hemi = 1
        if (j <= jeq) hemi = -1
        do i=I0,I1
          do pft = 1, N_COVERTYPES
            call prescr_veg_albedo(hemi,pft,jday,
     &           albedodata(:,pft,i,j))
          end do
        enddo
      enddo

      end subroutine prescr_veg_albedodata

!**************************************************************************

      subroutine prescr_get_carbonplant(IM,JM,I0,I1,J0,J1,
     &     laidata, hdata, dbhdata, popdata, cpooldata)
      !*  Calculate per plant carbon pools (g-C/plant).
      !*  After Moorcroft, et al. (2001).

      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      real*8,intent(in) :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(in) :: hdata(N_COVERTYPES)
      real*8,intent(in) :: dbhdata(N_COVERTYPES)
      real*8,intent(in) :: popdata(N_COVERTYPES) 
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      !Array: 1-foliage, 2-sapwood, 3-hardwood, 4-labile,5-fine root, 6-coarse root
      !----Local----
      integer :: n,pft !@var pft vegetation type
      integer :: hemi !@var hemi =1 in N. hemisphere, =-1 in South
      integer i,j,jeq

      jeq = JM/2

      cpooldata(:,:,:,:) = 0.0  !Zero initialize
      do j=J0,J1
        hemi = 1
        if (j <= jeq) hemi = -1
        do i=I0,I1
          do pft=1,N_PFT
            n = pft + COVEROFFSET
            call prescr_plant_cpools(pft,laidata(n,i,j),hdata(n),
     &           dbhdata(n), popdata(n),cpooldata(n,:,i,j))
          enddo
        enddo
      enddo
      
      end subroutine prescr_get_carbonplant

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
      print *,IM,JM,N_COVERTYPES !#DEBUG
      read(iu_SOIL) buf
      call closeunit(iu_SOIL)
      print *,"soil fractions:",buf(I0,J0,:)!#DEBUG

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

