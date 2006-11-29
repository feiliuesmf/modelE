      module ent_GISSveg

      use ent_const
      use ent_types
      use ent_pfts
      !use GCM_module, only:  GCMi, GCMj !Fix to names from GCM

      implicit none
      private
      save


      public 
     &     GISS_vegdata,ent_GISS_vegupdate,
     &     GISS_get_vdata, GISS_get_laidata, 
     &     GISS_get_hdata, GISS_get_initnm,
     &     GISS_update_vegcrops, GISS_phenology,
     &     GISS_veg_albedo,GISS_calc_rootprof,
     &     GISS_calcconst
!     &     GCM_get_grid, GCM_get_time, GCM_getdrv_cell, GCM_EWB,
      public GISS_calc_shc,GISS_veg_albedodata,GISS_plant_cpools


!*********************************************************************
!--- sand tundr  grass  shrub  trees  decid evrgr  rainf crops bdirt algae  c4grass
      real*8, parameter :: alamax(N_COVERTYPES) =
     $     (/ 0.d0, 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,10.0d0,8.0d0,4.5d0
     &     ,0.d0, 0.d0, 2.d0 /)
      real*8, parameter :: alamin(N_COVERTYPES) =
     $     (/ 0.d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 8.0d0,6.0d0,1.0d0
     &     ,0.d0, 0.d0, 1.d0 /)
      integer, parameter :: laday(N_COVERTYPES) =
     $     (/ 0, 196,  196,  196,  196,  196,  196,  196,  196
     &     ,0, 0, 196 /)

      real*8,parameter :: EDPERY=365. !GISS CONST.f

      !*********************************************************************
      contains
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************

!***************************************************************************
      subroutine GISS_calcconst() 
      !*    SUBROUTINE TO CALCULATE CONSTANT ARRAYS                     
      use ent_const
      !--Local------
      integer :: n
      real*8 :: lnscl

      !* annK - Turnover time of litter and soil carbon *!
      do n = 1, N_PFT
        if(pfpar(n)%lrage.gt.0.0d0)then
          annK(n,LEAF)  = (1.0d0/(pfpar(n)%lrage*secpy))
          annK(n,FROOT) = (1.0d0/(pfpar(n)%lrage*secpy))
        else
          annK(n,LEAF)  = 1.0d-40  !CASA originally 1.0e-40
          annK(n,FROOT) = 1.0d-40
        end if

        if(pfpar(n)%woodage.gt.0.0d0)then
          annK(n,WOOD)  = 1.0d0/(pfpar(n)%woodage*secpy)
        else
          annK(n,WOOD)  = 1.0d-40
        end if
       !* iyf: 1/(turnover times) for dead pools.  Want annK in sec-1.
        annK(n,SURFMET)    = 14.8d0  /secpy
        annK(n,SURFMIC)    = 6.0d0   /secpy
        annK(n,SURFSTR)    = 3.9d0   /secpy
        annK(n,SOILMET)    = 18.5d0  /secpy
        annK(n,SOILMIC)    = 7.3d0   /secpy
        annK(n,SOILSTR)    = 4.9d0   /secpy        ! 4.8 in casa v3.0
        annK(n,CWD)        = 0.2424d0/secpy
        annK(n,SLOW)       = 0.2d0   /secpy
        annK(n,PASSIVE)    = 0.1d0 * 0.02d0  /secpy
      enddo

#ifdef DEBUG
      do n = 1,N_PFT
        write(98,*) 'pft',n
        write(98,*) 'annK(n,LEAF)',annK(n,LEAF)
        write(98,*) 'annK(n,FROOT)',annK(n,FROOT)
        write(98,*) 'annK(n,WOOD)' ,annK(n,WOOD) 
        write(98,*) 'annK(n,SURFMET)',annK(n,SURFMET)
        write(98,*) 'annK(n,SURFMIC)',annK(n,SURFMIC)
        write(98,*) 'annK(n,SURFSTR)',annK(n,SURFSTR)
        write(98,*) 'annK(n,SOILMET)',annK(n,SOILMET)
        write(98,*) 'annK(n,SOILMIC)',annK(n,SOILMIC)
        write(98,*) 'annK(n,SOILSTR)',annK(n,SOILSTR)
        write(98,*) 'annK(n,CWD)',annK(n,CWD)    
        write(98,*) 'annK(n,SLOW)',annK(n,SLOW)   
        write(98,*) 'annK(n,PASSIVE)',annK(n,PASSIVE)
      enddo
#endif

      !* solubfrac - Soluble fraction of litter*!
      !structurallignin, lignineffect - frac of structural C from lignin, effect of lignin on decomp -PK 6/29/06 
      do n = 1,N_PFT
        lnscl = pfpar(n)%lit_C2N * pfpar(n)%lignin * 2.22 !lignin:nitrogen scalar        
        solubfract(n) = 0.85 - (0.018 * lnscl)
        structurallignin(n) = (pfpar(n)%lignin * 0.65 * 2.22) 
     &                      / (1. - solubfract(n))
        lignineffect(n) = exp(-3.0 * structuralLignin(n))
      end do
      print *,'from GISS_calcconst (ent_GISSveg): lignineffect[n_pft]='
     &       ,lignineffect

      end subroutine GISS_calcconst

      !*********************************************************************
      !*********************************************************************
      !*    SUBROUTINES TO READ IN GISS VEGETATION DATA SETS 
      !*********************************************************************

!***************************************************************************
      subroutine GISS_vegdata(jday, year, IM,JM,I0,I1,J0,J1,
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
      real*8,intent(out) :: soil_texture(N_SOIL_TYPES,I0:I1,J0:J1)
      !-----Local------

      call GISS_get_vdata(IM,JM,I0,I1,J0,J1,vegdata)   !veg fractions
      call GISS_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)
      call GISS_get_laidata(jday,JM,I0,I1,J0,J1,laidata) !lai
      call GISS_update_vegcrops(year,IM,JM,I0,I1,J0,J1,vegdata)

      call GISS_get_hdata(hdata) !height
      call GISS_get_initnm(nmdata) !nm
      call GISS_get_rootprof(rootprofdata)
      call GISS_get_woodydiameter(hdata,dbhdata)
      call GISS_get_pop(dbhdata,popdata)
      call GISS_get_crownrad(popdata,craddata)
      call GISS_get_carbonplant(IM,JM,I0,I1,J0,J1,
     &     laidata,hdata,dbhdata,popdata,cpooldata)
      call GISS_get_soil_types(IM,JM,I0,I1,J0,J1,
     &     soil_color,soil_texture)

      print*,'vegdata(:,I1,J1)',vegdata(:,I1,J1)
      print*,'hdata',hdata
      print*,'nmdata',nmdata
      print*,'dbhdata',dbhdata
      print*,'popdata',popdata
      print*,'craddata',craddata
!      print*,'cpooldata(GRASSC3+COVEROFFSET,:,I1,J1)',
!     &     cpooldata(GRASSC3+COVEROFFSET,:,I1,J1)
!      print*,'cpooldata(SHRUB+COVEROFFSET,:,I1,J1)',
!     &     cpooldata(SHRUB+COVEROFFSET,:,I1,J1)
!      print*,'cpooldata(SAVANNA+COVEROFFSET,:,I1,J1)',
!     &     cpooldata(SAVANNA+COVEROFFSET,:,I1,J1)
      print*,'soil_color',soil_color

      end subroutine GISS_vegdata

!***************************************************************************
      subroutine ent_GISS_vegupdate(dt,entcell,hemi,jday,year,
     &     update_crops, update_soil)
      use patches, only : summarize_patch
      use entcells,only : summarize_entcell,entcell_extract_pfts
      use phenology,only : litter !### Igor won't like this here.

      real*8,intent(in) :: dt
      type(entcelltype) :: entcell
      integer,intent(in) :: jday,year,hemi
      logical,intent(in) :: update_crops
      logical,intent(in) :: update_soil
      !----Local------
      type(patch),pointer :: pp
      type(cohort),pointer :: cop
      real*8 vdata(N_COVERTYPES) ! needed for a hack to compute canopy

      pp => entcell%oldest
      cop => pp%tallest
      do while (ASSOCIATED(pp))
        !* Soil *!
        if (update_soil)  call litter(dt,pp) !###Dependency tree?

        !* LAI, ALBEDO *!
        call GISS_phenology(jday,hemi, pp)  !## SHOULD HAVE update_crops here
        if (update_crops) then
        ! re-read crops data and update the vegetation
      ! this function is located up in the dependency tree
      ! can't be called here ... IA
        endif  

        call summarize_patch(pp)
        pp => pp%younger
      end do

      call summarize_entcell(entcell)

      vdata(:) = 0.d0
      call entcell_extract_pfts(entcell, vdata(2:) )
      entcell%heat_capacity=GISS_calc_shc(vdata)

      !if (YEAR_FLAG.eq.0) call ent_GISS_init(entcellarray,im,jm,jday,year)
      !!!### REORGANIZE WTIH ent_prog.f ####!!!
      
      end subroutine ent_GISS_vegupdate

!***************************************************************************
      subroutine GISS_get_vdata(im,jm,I0,I1,J0,J1,vdata)
      !* This version reads in vegetation structure from GISS data set.
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

      do k=1,10  !## HACK - this should be related to N_COVERTYPES

        read(iu_VEG) title, buf
        vdata(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        print *,"read VEG:", title
      end do

      call closeunit(iu_VEG)
      end subroutine GISS_get_vdata

!**************************************************************************
      subroutine GISS_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata)
      !* This version reads in crop distribution from GISS data set.
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
        year1 = year2
        crop1(:,:) = crop2(:,:)
        read (iu_CROPS,end=10) title, crop4
        read(title,*) year2
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
      end subroutine GISS_get_cropdata

!**************************************************************************

      subroutine GISS_update_vegcrops(year,IM,JM,I0,I1,J0,J1,
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
      call GISS_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata) !crop fraction

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
 
      end subroutine GISS_update_vegcrops

!**************************************************************************

      subroutine GISS_get_laidata(jday,JM,I0,I1,J0,J1,laidata)
!@sum Returns GISS GCM leaf area index for entire grid and given jday.
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
            laidata(n,i,j) = GISS_calc_lai(n,jday,hemi)
          enddo
        enddo
      enddo

      !* Return lai for each vegetation type.
      end subroutine GISS_get_laidata

!**************************************************************************
      
      real*8 function GISS_calc_lai(pnum,jday,hemi ) RESULT(lai)
!@sum Returns GISS GCM leaf area index for given vegetation type, julian day
!@+   and hemisphere
      use ent_const
      !real*8, intent(out) :: lai !@var lai leaf area index - returned
      integer, intent(in) :: pnum !@var pnum cover type
      integer, intent(in) :: jday !@var jday julian day
      integer, intent(in) :: hemi !@var hemi =1 in N. hemisphere, =-1 S.hemi
      !-----Local variables------
      real*8 dphi

      dphi = 0
      if ( hemi < 0 ) dphi = 2d0*pi*.5d0

      !* Return lai *!
      lai =  .5d0 * (alamax(pnum) + alamin(pnum))
     $     + .5d0 * (alamax(pnum) - alamin(pnum))
     $     * cos( 2d0*pi*(laday(pnum)-jday)/dble(EDPERY) + dphi )

      end function GISS_calc_lai


!**************************************************************************
      
      real*8 function GISS_calc_shc(vdata) RESULT(shc)
!@sum Returns GISS GCM specific heat capacity for canopy
      use ent_const
      real*8, intent(in) :: vdata(:) !@var pnum cover type
      !-----Local variables------
      real*8 lai, fsum
      integer pnum

      lai = 0.d0
      fsum = 0.d0

      do pnum=COVEROFFSET+1,COVEROFFSET+N_PFT
        lai = lai + .5d0 * (alamax(pnum) + alamin(pnum))*vdata(pnum)
        fsum = fsum + vdata(pnum)
      enddo

      if ( fsum > EPS ) lai = lai/fsum

      shc = (.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow

      end function GISS_calc_shc


!**************************************************************************
      real*8 function GISS_calc_shoot(pnum,hdata,dbhdata) Result(Bshoot)
!@sum Returns GISS GCM veg shoot kg-C per plant for given vegetation type.
!@+   From Moorcroft, et al. (2001), who takes allometry data from
!@+   Saldarriaga et al. (1998).
      integer,intent(in) :: pnum !@var pnum vegetation type
      real*8,intent(in) :: hdata(N_COVERTYPES), dbhdata(N_COVERTYPES)
      !-----Local-------
      real*8 :: wooddens
      integer :: n !covertypes index

      n = pnum + COVEROFFSET
      wooddens = wooddensity_gcm3(pnum)
      Bshoot = 0.069 * (hdata(n))**0.572
     &     * (dbhdata(n))**1.94 * (wooddens**0.931)

      end function GISS_calc_shoot

!**************************************************************************
      subroutine GISS_phenology(jday,hemi, pp)
      !* Calculate new LAI and albedo for given jday, for GISS vegetation. *!
      !* TBA:  THIS ROUTINE WILL ALSO UPDATE LIVE BIOMASS POOLS.           *!
      use ent_types
      use ent_pfts
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

!**************************************************************************
      subroutine GISS_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)
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
            call GISS_veg_albedo(hemi,pft,jday,
     &           albedodata(:,pft,i,j))
          end do
        enddo
      enddo

      end subroutine GISS_veg_albedodata

!**************************************************************************

      subroutine GISS_veg_albedo(hemi, pft, jday, albedo)
!@sum returns albedo for vegetation of type pft 
!@+   as it is computed in GISS modelE
      integer, intent(in) :: hemi !@hemi hemisphere (-1 south, +1 north)
      integer, intent(in) :: pft !@var pftlike iv, plant functional type
      integer, intent(in) :: jday !@jday julian day
      real*8, intent(out) :: albedo(N_BANDS) !@albedo returned albedo
      !----------Local----------
      integer, parameter :: NV=12
      !@var SEASON julian day for start of season (used for veg albedo calc)
C                      1       2       3       4
C                    WINTER  SPRING  SUMMER  AUTUMN
      real*8, parameter, dimension(4)::
     *     SEASON=(/ 15.00,  105.0,  196.0,  288.0/)
C**** parameters used for vegetation albedo
!@var albvnd veg alb by veg type, season and band
      real*8, parameter :: ALBVND(NV,4,6) = RESHAPE( (/
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.067,.089,.089,.078,.100,.067,.061,.089,.000,.200,.089,
     2 .500,.062,.100,.100,.073,.055,.067,.061,.100,.000,.200,.100,
     3 .500,.085,.091,.139,.085,.058,.083,.061,.091,.000,.200,.091,
     4 .500,.080,.090,.111,.064,.055,.061,.061,.090,.000,.200,.090,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.080,.107,.107,.093,.120,.080,.073,.107,.000,.200,.107,
     2 .500,.082,.140,.120,.096,.083,.080,.073,.140,.000,.200,.140,
     3 .500,.119,.145,.167,.119,.115,.100,.073,.145,.000,.200,.145,
     4 .500,.102,.126,.132,.081,.087,.073,.073,.126,.000,.200,.126,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.100,.133,.133,.116,.150,.100,.091,.133,.000,.200,.133,
     2 .500,.103,.175,.150,.120,.109,.100,.091,.175,.000,.200,.175,
     3 .500,.148,.182,.208,.148,.144,.125,.091,.182,.000,.200,.182,
     4 .500,.127,.157,.166,.102,.109,.091,.091,.157,.000,.200,.157,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.020,.027,.027,.023,.030,.020,.018,.027,.000,.200,.027,
     2 .500,.021,.035,.030,.024,.022,.020,.018,.035,.000,.200,.035,
     3 .500,.030,.036,.042,.030,.029,.025,.018,.036,.000,.200,.036,
     4 .500,.026,.032,.033,.020,.022,.018,.018,.032,.000,.200,.032
     *     /),(/NV,4,6/) )
C
ccc or pass k-vegetation type, L-band and 1 or 2 for Hemisphere
      integer k,kh1,kh2,l
      real*8 seasn1,seasn2,wt2,wt1
c
c                      define seasonal albedo dependence
c                      ---------------------------------
c
      seasn1=-77.0d0
      do k=1,4
        seasn2=SEASON(k)
        if(jday.le.seasn2) go to 120
        seasn1=seasn2
      end do
      k=1
      seasn2=380.0d0
  120 continue
      wt2=(jday-seasn1)/(seasn2-seasn1)
      wt1=1.d0-wt2
      if ( hemi == -1 ) then    ! southern hemisphere
        kh1=1+mod(k,4)
        kh2=1+mod(k+1,4)
      else                      ! northern hemisphere
        kh1=1+mod(k+2,4)
        kh2=k
      endif

      do l=1,6
        albedo(l)=wt1*ALBVND(pft,kh1,l)+wt2*ALBVND(pft,kh2,l)
      enddo

      end subroutine GISS_veg_albedo

!**************************************************************************
      subroutine GISS_calc_rootprof(rootprof, pnum)
      !Return array rootprof of fractions of roots in soil layer
      !Cohort/patch level.
      real*8 :: rootprof(:)
      integer :: pnum !plant functional type
      !-----Local variables------------------
      real*8,parameter :: dz_soil(1:6)=  !N_DEPTH
     &     (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
 !--- tundr  grass  shrub  trees  decid evrgr  rainf crops bdirt algae  c4grass
      real*8, parameter :: aroot(N_COVERTYPES) = 
     $     (/ 0.d0,12.5d0, 0.9d0, 0.8d0,0.25d0,0.25d0,0.25d0,1.1d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: broot(N_COVERTYPES) = 
     $     (/ 0.d0, 1.0d0, 0.9d0, 0.4d0,2.00d0,2.00d0,2.00d0,0.4d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      integer :: n,l
      real*8 :: z, frup,frdn

c**** calculate root fraction afr averaged over vegetation types
      !Initialize zero
      do l=1,N_DEPTH
        rootprof(l) = 0.0
      end do
      do n=1,N_DEPTH
        if (dz_soil(n) <= 0.0) exit !Get last layer w/roots in it.
      end do
      n=n-1
      z=0.
      frup=0.
      do l=1,n
        z=z+dz_soil(l)
        frdn=aroot(pnum)*z**broot(pnum) !cumulative root distrib.
        !frdn=min(frdn,one)
        frdn=min(frdn,1d0)
        if(l.eq.n)frdn=1.
        rootprof(l) = frdn-frup
        frup=frdn
      end do
      !Return rootprof(:)
      end subroutine GISS_calc_rootprof
!**************************************************************************
      subroutine GISS_get_rootprof(rootprofdata)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH) 
      !---Local--------
      integer :: pnum !plant functional type      

      do pnum=1,N_COVERTYPES
        call GISS_calc_rootprof(rootprofdata(pnum,:), pnum)
        !Return array rootprof of fractions of roots in soil layer
        !by vegetation type.
      end do
      end subroutine GISS_get_rootprof

!**************************************************************************

      subroutine GISS_get_hdata(hdata)
      !* Return array parameter of GISS vegetation heights.
      real*8 :: hdata(N_COVERTYPES) 
      !------
      real*8, parameter :: vhght(N_COVERTYPES) =
      !* bsand tundrv  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
     $     (/0.d0, 0.1d0, 1.5d0,   5d0,  15d0,  20d0,  30d0, 25d0,1.75d0
     &     ,0.d0, 0.d0, 1.5d0 /)

      !* Copy GISS code for calculating seasonal canopy height here.
      ! For GISS Model E replication, don't need to fill in an
      ! i,j array of vegetation height, but just can use
      ! constant arry for each vegetation pft type.
      ! For full-fledged Ent, will need to read in a file entdata
      ! containing vegetation heights.

      !* Return hdata heights for all vegetation types
      hdata = vhght
      end subroutine GISS_get_hdata

!**************************************************************************
      subroutine GISS_get_initnm(nmdata)
!@sum  Mean canopy nitrogen (nmv; g/m2[leaf])
      real*8 :: nmdata(N_COVERTYPES)
      !-------
      real*8, parameter :: nmv(N_COVERTYPES) =
     $     (/0.d0,1.6d0,0.82d0,2.38d0,1.03d0,1.25d0,2.9d0,2.7d0,2.50d0
     &     ,0.d0, 0.d0, 0.82d0 /)

      !* Return intial nm for all vegetation and cover types
      nmdata = nmv
      end subroutine GISS_get_initnm

!*************************************************************************
      subroutine GISS_get_pop(dbhdata,popdata)
      !* Return array of GISS vegetation population density (#/m2)
      !* Derived from Moorcroft, et al. (2001)
      real*8,intent(in) :: dbhdata(N_COVERTYPES)
      real*8,intent(out) :: popdata(N_COVERTYPES)
      !---Local-----------
      integer :: n,pft
      real*8 :: Blmax, wooddens

      popdata(:) = 0.0 !Zero initialize, and zero bare soil.
      do pft=1,N_PFT
        n = pft + COVEROFFSET
        if (pft.eq.GRASSC3) then
          popdata(n) = 1.0      !Grass is just a large ensemble
        else
          wooddens = wooddensity_gcm3(pft)
          Blmax = 0.0419 * ((dbhdata(n))**1.56) * (wooddens**0.55)
          popdata(n) = (alamax(n)/pfpar(pft)%sla)/Blmax
          !print*,'pft,wooddens,Blmax,popd',pft,wooddens,Blmax,popdata(n)
        endif
      enddo

      end subroutine GISS_get_pop

!*************************************************************************
      subroutine GISS_get_woodydiameter(hdata, wddata)
      !* Return array of woody plant diameters at breast height (dbh, cm)
      real*8,intent(in) :: hdata(N_COVERTYPES)
      real*8,intent(out) :: wddata(N_COVERTYPES)
      !----Local---------
      integer :: n,pft

      wddata(:) = 0.0 !Zero initialize.
      do pft = 1,N_PFT
        n = pft + COVEROFFSET
        if (pft.ne.GRASSC3) then !Woody
          if (pft.eq.TUNDRA) then
            wddata(n) = ED_woodydiameter(pft,hdata(n)) * 20 !FUDGE UNTIL HAVE MIXED CANOPIES
          else                  !Most trees
            wddata(n) = ED_woodydiameter(pft,hdata(n))
          end if
        endif
      enddo
      end subroutine GISS_get_woodydiameter
!*************************************************************************
      real*8 function ED_woodydiameter(pft,h) Result(dbh)
      !* Return woody plant diameter (m).
      !* From Moorcroft, et al. (2001)
      integer,intent(in) :: pft !plant functional type
      real*8,intent(in) ::  h !height (m)
      !real*8,intent(out) :: dbh !(cm)

      if (pft.eq.SAVANNA) then
        dbh = 30.0 !Estimate from Tonzi Ranch, NYK
      else
        dbh = ((1/2.34)*h)**(1/0.64)
      endif
      
      end function ED_woodydiameter
!*************************************************************************
      subroutine GISS_get_crownrad(popdata,craddata)
      real*8,intent(in) :: popdata(N_COVERTYPES)
      real*8,intent(out) :: craddata(N_COVERTYPES)
      !---Local----
      integer :: n, pft

      craddata(:) = 0.0 !Zero initialize.
      do pft=1,N_PFT
        n = pft + COVEROFFSET
        craddata(n) = 0.5*sqrt(1/popdata(n))
      end do

      end subroutine GISS_get_crownrad

!*************************************************************************
      subroutine GISS_get_carbonplant(IM,JM,I0,I1,J0,J1,
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
            call GISS_plant_cpools(pft,laidata(n,i,j),hdata(n),
     &           dbhdata(n), popdata(n),cpooldata(n,:,i,j))
          enddo
        enddo
      enddo
      
      end subroutine GISS_get_carbonplant

!*************************************************************************

      subroutine GISS_plant_cpools(pft, lai, h, dbh, popdens, cpool )
      !* Calculate plant carbon pools for single plant (g-C/plant)
      !* After Moorcroft, et al. (2001).
      integer,intent(in) :: pft !plant functional type
      real*8, intent(in) :: lai,h,dbh,popdens  !lai, h(m), dbh(cm),popd(#/m2)
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      !----Local------
      
      cpool(FOL) = lai/pfpar(pft)%sla/popdens *1e3!Bl
      cpool(FR) = cpool(FOL)   !Br
      cpool(LABILE) = 0.d0      !dummy
      if (pft.ne.GRASSC3) then  !Woody
        cpool(SW) = 0.00128 * pfpar(pft)%sla * cpool(FR) * h *1e3 !Bsw
        cpool(HW) = 0.069*(h**0.572)*(dbh**1.94) * wooddensity_gcm3(pft)
     &       *1e3
        cpool(CR) = 0.d0        !dummy
      else
        cpool(SW) = 0.d0
        cpool(HW) = 0.d0
        cpool(CR) = 0.d0
      endif
      !write(97,*) pft, lai, h, dbh, popdens, cpool
      end subroutine GISS_plant_cpools
!*************************************************************************
      real*8 function wooddensity_gcm3(pft) Result(wooddens)
      integer,intent(in) :: pft
      !* Wood density (g cm-3). Moorcroft et al. (2001).

      wooddens = max(0.5d0, 0.5d0 + 0.2d0*(pfpar(pft)%lrage-1.d0))

      end function wooddensity_gcm3
!*************************************************************************

      subroutine GISS_get_soil_types(im,jm,I0,I1,J0,J1,
     &     soil_color,soil_texture)
      !* Return arrays of GISS soil color and texture.
      use FILEMANAGER, only : openunit,closeunit
      integer, intent(in) :: im,jm,I0,I1,J0,J1
      integer, intent(out) :: soil_color(N_COVERTYPES)
      real*8, intent(out) :: soil_texture(N_SOIL_TYPES,I0:I1,J0:J1)
      !------
      integer, parameter :: soil_color_prescribed(N_COVERTYPES) =
      !* bsand tundr  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
     $     (/1, 2, 2,  2, 2, 2, 2, 2
     &     ,2, 2, 2, 2 /)
      real*8 :: buf(im,jm,N_SOIL_TYPES)
      integer :: iu_SOIL
      integer k

      soil_color(:) = soil_color_prescribed(:)

      call openunit("soil_textures",iu_SOIL,.true.,.true.)
      read(iu_SOIL) buf
      call closeunit(iu_SOIL)

      do k=1,N_SOIL_TYPES
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
      end subroutine GISS_get_soil_types
!*************************************************************************
      end module ent_GISSveg

