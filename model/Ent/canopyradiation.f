      module canopyrad
!@sum Routines for calculating canopy radiation and albedo.
!@auth W.Ni-Meister

      !Ent MODULES TO USE
      use ent_const
      use ent_types

      implicit none
      private
      save

      public recalc_radpar, get_patchalbedo

      contains
      !*********************************************************************
      
      subroutine get_patchalbedo(gcmtime,pp)
      !* Parse through cohorts to calculate patch albedo via GORT clumping.
      !* This version reads in vegetation structure from GISS data set.
      type(timestruct) :: gcmtime
      type(patch),pointer :: pp
      !------
      type(cohort),pointer :: cop

      do cop=pp%tallest
      !* Copy GISS code for calculating seasonal aalbveg here.
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
        !Get albedo of cop%pft for given gcmtime
        !   - Summarize foliage density in layers
        !Calculate patch albedo using GORT clumping.
        cop = cop%shorter
        if DEALLOCATED(cop) then exit
      end do
      
      !* Return seasonal albedo for patch
      !   - To reproduce GISS albedoes, just assign the one cohort's albedo
      !     to the patch.
      !  should be same as assigning
      !  pp%albedo = aalbveg(i,j)
      end subroutine get_patchalbedo


      subroutine GISS_veg_albedo(pft, jday, hemi, albedo)
!@sum returns albedo for vegetation of type pft 
!@+   as it is computed in GISS modelE
      integer, intent(in) :: pft !@var pftlike iv, plant functional type
      integer, intent(in) :: jday !@jday julian day
      integer, intent(in) :: hemi !@hemi hemisphere (-1 south, +1 north)
      real*8, intent(out) :: albedo(6) !@albedo returned albedo

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
      integer ks1,ks2,kn1,kn2,l,kh1,kh2
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

      end subroutine get_gcm_albedo



      !*********************************************************************
      subroutine recalc_radpar(entcell)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.

      type(entcelltype) :: entcell
      type(patch),pointer :: pptr

      pptr = entcell%youngest

      do while (ASSOCIATED(pptr)) 
        call GORT_clumping(pptr)
        pptr = pptr%older
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      end do

      end subroutine recalc_radpar


      !*********************************************************************
      subroutine GORT_clumping(pptr)
!@sum Calculate the GORT clumping index in canopy layers and save into
!@sum variable ppt%crad
      type(patch),pointer :: pptr
      !--------------------------
      type(cohort),pointer :: cptrhi, cptrlo, cptr
      type(canradtype) :: crad

      cptrhi = pptr%tallest
      cptrlo = pptr%shortest
      cptr = cptrhi

      !DEFINE LAYERS FOR LAI AND CLUMPING INDICES
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------


      do while (ASSOCIATED(cptr)) 
        !Allocate leaf area/crown volume within canopy layers
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

        
        cptr = cptr%shorter  !next
      end do

      !Calculate GORT clumping index for each layer 
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      !* Return (pptr%crad)
      end subroutine GORT_clumping

      !*********************************************************************

      subroutine calc_canopy_rad(pptr, h)
!@sum Get incident light profiles in canopy and return in pptr%crad
      type(patch),pointer :: pptr
      real*8 :: h               !Height in canopy
      !--------------------------

      type(canradtype) :: crad

      real*8 :: Solarzen, Isw, IPAR, Ibeam, Idiff 

      Solarzen = pptr%cellptr%Solarzen
      Isw = pptr%cellptr%Isw
      IPAR = pptr%cellptr%IPAR
      Ibeam = pptr%cellptr%Ibeam
      Idiff = pptr%cellptr%Idiff

      !Get incident light profiles in canopy and return in crad
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------


      pptr%crad = crad
      end subroutine calc_canopy_rad
      
      !*********************************************************************
      subroutine get_canopy_rad(pptr,h)
!@sum Retrieve incident light at height h from pptr%crad
      type(patch),pointer :: pptr
      real*8 :: h               !Height in canopy
      !--------------------------
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine get_canopy_rad
      !*********************************************************************

      end module canopyrad
