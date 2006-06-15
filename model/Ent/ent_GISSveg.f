      module ent_GISSveg

      use ent_const, only : pi, N_COVERTYPES, N_DEPTH
     &     
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
     &     GISS_veg_albedo,GISS_calc_froot,
     &     GISS_calcconst
!     &     GCM_get_grid, GCM_get_time, GCM_getdrv_cell, GCM_EWB,


      !*********************************************************************

      real*8,parameter :: EDPERY=365. !GISS CONST.f

      !*********************************************************************
      contains
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*    SUBROUTINE TO CALCULATE CONSTANT ARRAYS                     
      subroutine GISS_calcconst() 
      use ent_const
      !--Local------
      integer :: n
      real*8, DIMENSION(N_PFT) :: lrage, woodage
      real*8 :: lnscl

      !* annK - Turnover time of litter and soil carbon *!
      do n = 1, N_PFT
        if(lrage(n).gt.0.0)then
          annK(n,LEAF)  = (1.0/(lrage(n)*secpy))
          annK(n,FROOT) = (1.0/(lrage(n)*secpy))
        else
          annK(n,LEAF)  = 1.0d-40  !CASA originally 1.0e-40
          annK(n,FROOT) = 1.0d-40
        end if

        if(woodage(n).gt.0.0)then
          annK(n,WOOD)  = 1.0/(woodage(n)*secpy)
        else
          annK(n,WOOD)  = 1.0d-40
        end if
       !* iyf: 1/(turnover times) for dead pools.  Want annK in sec-1.
        annK(n,SURFMET)    = 14.8  /secpy
        annK(n,SURFMIC)    = 6.0   /secpy
        annK(n,SURFSTR)    = 3.9   /secpy
        annK(n,SOILMET)    = 18.5  /secpy
        annK(n,SOILMIC)    = 7.3   /secpy
        annK(n,SOILSTR)    = 4.9   /secpy        ! 4.8 in casa v3.0
        annK(n,CWD)        = 0.2424/secpy
        annK(n,SLOW)       = 0.2   /secpy
        annK(n,PASSIVE)    = 0.1 * 0.02  /secpy
      enddo

      !* solubfrac - Soluble fraction of litter*!
      do n = 1,N_PFT
        lnscl = pfpar(n)%lit_C2N * pfpar(n)%lignin * 2.22 !lignin:nitrogen scalar        
        solubfract(n) = 0.85 - (0.018 * lnscl)
      end do


      !*******************************************************************
      !* INSERT CALCULATION OF SOIL TEXTURE PARAMETER FOR MICROBE TURNOVER
      !*******************************************************************

      end subroutine GISS_calcconst

      !*********************************************************************
      !*********************************************************************
      !*    SUBROUTINES TO READ IN GISS VEGETATION DATA SETS 
      !*********************************************************************

      subroutine GISS_vegdata(jday, year, IM,JM,I0,I1,J0,J1,
     &     vegdata,albedodata,laidata,hdata,nmdata,frootdata,
     &     popdata,soil_color,soil_texture)
      integer,intent(in) :: jday, year
      integer,intent(in) :: IM,JM,I0,I1,J0,J1 !long/lat grid number range
      real*8,intent(out) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: laidata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: hdata(N_COVERTYPES)
      real*8,intent(out) :: nmdata(N_COVERTYPES)
      real*8,intent(out) :: frootdata(N_COVERTYPES,N_DEPTH)
      real*8,intent(out) :: popdata(N_COVERTYPES)
      integer,intent(out) :: soil_color(N_COVERTYPES)
      real*8,intent(out) :: soil_texture(N_SOIL_TYPES,I0:I1,J0:J1)
 
      !-----Local------

      call GISS_get_vdata(IM,JM,I0,I1,J0,J1,vegdata)   !veg fractions
      call GISS_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)
      call GISS_get_laidata(jday,JM,I0,I1,J0,J1,laidata) !lai
      call GISS_update_vegcrops(year,IM,JM,I0,I1,J0,J1,vegdata)

      call GISS_get_hdata(hdata) !height
      call GISS_get_initnm(nmdata) !nm
      call GISS_get_froot(frootdata)
      call GISS_get_pop(popdata)
      call GISS_get_soil_types(IM,JM,I0,I1,J0,J1,
     &     soil_color,soil_texture)

      end subroutine GISS_vegdata

      !*********************************************************************
      subroutine ent_GISS_vegupdate(entcell,hemi,jday,year,update_crops)
      use patches, only : summarize_patch
      use entcells,only : summarize_entcell
      type(entcelltype) :: entcell
      integer,intent(in) :: jday,year,hemi
      logical,intent(in) :: update_crops
      !----Local------
      type(patch),pointer :: pp

      pp => entcell%oldest
      do while (ASSOCIATED(pp))
        call GISS_phenology(jday,hemi, pp)
        call summarize_patch(pp)
        pp => pp%younger
      end do

      call summarize_entcell(entcell)

      if (update_crops) then
        ! re-read crops data and update the vegetation
      endif

      ! this function is located up in the dependency tree
      ! can't be called here ... IA
      !if (YEAR_FLAG.eq.0) call ent_GISS_init(entcellarray,im,jm,jday,year)
      !!!### REORGANIZE WTIH ent_prog.f ####!!!
      
      end subroutine ent_GISS_vegupdate
      !*********************************************************************

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

      ! make sure that unused fractions are set to 0
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
      integer :: pft !@var pft vegetation type
      integer :: hemi !@var hemi =1 in N. hemisphere, =-1 in South
      integer i,j,jeq

      jeq = JM/2

      do j=J0,J1
        hemi = 1
        if (j <= jeq) hemi = -1
        do i=I0,I1
          do pft=1,N_COVERTYPES
            laidata(pft,i,j) = GISS_calc_lai(pft,jday,hemi)
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
      integer, intent(in) :: pnum !@var pnum vegetation type
      integer, intent(in) :: jday !@var jday julian day
      integer, intent(in) :: hemi !@var hemi =1 in N. hemisphere, =-1 S.hemi
      !-----Local variables------
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
      real*8 dphi

      dphi = 0
      if ( hemi < 0 ) dphi = 2d0*pi*.5d0

      !* Return lai *!
      lai =  .5d0 * (alamax(pnum) + alamin(pnum))
     $     + .5d0 * (alamax(pnum) - alamin(pnum))
     $     * cos( 2d0*pi*(laday(pnum)-jday)/dble(EDPERY) + dphi )

      end function GISS_calc_lai


!**************************************************************************
      subroutine GISS_phenology(jday,hemi, pp)
      use ent_types
      use ent_pfts
      implicit none
      integer,intent(in) :: jday !Day of year.
      integer,intent(in) :: hemi !@var hemi -1: S.hemisphere, 1: N.hemisphere
      type(patch),pointer :: pp
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: coplaiprev

      real*8 :: laip  !patch-level summary of LAI
!      real*8 :: laig  !entcell grid-level summary of LAI
      if (ASSOCIATED(pp)) then
        laip = 0.0
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          coplaiprev = cop%lai
          cop%lai = GISS_calc_lai(cop%pft, jday, hemi)
          !call GISS_calc_litter(cop, OUT POOLS HERE)
          laip = laip + cop%LAI
          cop => cop%shorter
        end do
        pp%sumcohort%LAI = laip
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
      subroutine GISS_calc_froot(froot, pnum)
      !Return array froot of fractions of roots in soil layer
      !Cohort/patch level.
      real*8 :: froot(:)
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
        froot(l) = 0.0
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
        froot(l) = frdn-frup
        frup=frdn
      end do
      !Return froot(:)
      end subroutine GISS_calc_froot
!**************************************************************************
      subroutine GISS_get_froot(frootdata)
      real*8,intent(out) :: frootdata(N_COVERTYPES,N_DEPTH) 
      !---Local--------
      integer :: pnum !plant functional type      

      do pnum=1,N_COVERTYPES
        call GISS_calc_froot(frootdata(pnum,:), pnum)
        !Return array froot of fractions of roots in soil layer
        !by vegetation type.
      end do
      end subroutine GISS_get_froot

!**************************************************************************

      subroutine GISS_get_hdata(hdata)
      !* Return array parameter of GISS vegetation heights.
      real*8 :: hdata(N_COVERTYPES) 
      !------
      real*8, parameter :: vhght(N_COVERTYPES) =
      !* bsand tundr  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
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
      subroutine GISS_get_pop(popdata)
      !* Return array parameter of GISS vegetation heights.
      real*8 :: popdata(N_COVERTYPES) 
      !------
      real*8, parameter :: popdens(N_COVERTYPES) =
      !* xx tundr  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
     $     (/0.d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0
     &     ,1.0d0 ,0.d0, 0.d0, 1.0d0 /)

      ! For GISS Model E replication, don't need to fill in an
      ! i,j array of vegetation density, but just can use
      ! constant arry for each vegetation pft type.
      ! For full-fledged Ent, will need to read in a file entdata
      ! containing vegetation density.

      !* Return popdata population density for all vegetation types
      popdata = popdens

      end subroutine GISS_get_pop
!*************************************************************************

      subroutine GISS_get_soil_types(im,jm,I0,I1,J0,J1,
     &     soil_color,soil_texture)
      !* Return array parameter of GISS vegetation heights.
      use FILEMANAGER, only : openunit,closeunit
      integer, intent(in) :: im,jm,I0,I1,J0,J1
      integer, intent(out) :: soil_color(N_COVERTYPES)
      real*8, intent(out) :: soil_texture(N_SOIL_TYPES,I0:I1,J0:J1)
      !------
      integer, parameter :: soil_color_prescribed(N_COVERTYPES) =
      !* tundr  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
     $     (/1, 2, 2,  2, 2, 2, 2, 2
     &     ,2, 2, 2, 2 /)
      real*8 :: buf(im,jm,N_SOIL_TYPES)
      integer :: iu_SOIL
      integer k,i,j
      ! For GISS Model E replication, don't need to fill in an
      ! i,j array of vegetation density, but just can use
      ! constant arry for each vegetation pft type.
      ! For full-fledged Ent, will need to read in a file entdata
      ! containing vegetation density.

      !* Return popdata population density for all vegetation types
      soil_color(:) = soil_color_prescribed(:)

      call openunit("soil_textures",iu_SOIL,.true.,.true.)
      read(iu_SOIL) buf
      call closeunit(iu_SOIL)

      do k=1,N_SOIL_TYPES
        soil_texture(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1,k)
      enddo

      do j=J0,J1
        do i=I0,I1
          print *, soil_texture(:,i,j) , sum(soil_texture(:,i,j))
        enddo
      enddo
      end subroutine GISS_get_soil_types
!*************************************************************************
      end module ent_GISSveg

