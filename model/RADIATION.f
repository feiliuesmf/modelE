        
      MODULE SURF_ALBEDO
!@sum SURF_ALBEDO contains parameters/variables needed for albedo calc
!@auth A. Lacis/V. Oinas (modifications by I. Alienov/G. Schmidt)
      implicit none
      save
!@var NKBAND number of K-bands
      integer, parameter :: NKBAND=33

!@var NVEG number of real vegetation types (not including bare soil)
!@var NV total number of vegetation types
      integer, parameter :: NVEG = 9, NV=11

!@var SRFOAM look up table for ocean foam as a function of wind speed
      real*8, parameter, dimension(25) :: SRFOAM = (/
     *     0.000,0.000,0.000,0.000,0.001,0.002,0.003,0.005,0.007,0.010,
     *     0.014,0.019,0.025,0.032,0.041,0.051,0.063,0.077,0.094,0.112,
     *     0.138,0.164,0.191,0.218,0.246/)

!@var SEASON julian day for start of season (used for veg albedo calc)
C                      1       2       3       4
C                    WINTER  SPRING  SUMMER  AUTUMN
      real*8, parameter, dimension(4)::
     *     SEASON=(/ 15.00,  105.0,  196.0,  288.0/)
C**** parameters used for vegetation albedo
!@var albvnd veg alb by veg type, season and band
      real*8, parameter :: ALBVND(NV,4,6) = RESHAPE( (/
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .067, .089, .089, .078, .100, .067, .061, .089, .000, .200,
     2 .500, .062, .100, .100, .073, .055, .067, .061, .100, .000, .200,
     3 .500, .085, .091, .139, .085, .058, .083, .061, .091, .000, .200,
     4 .500, .080, .090, .111, .064, .055, .061, .061, .090, .000, .200,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .200, .267, .267, .233, .300, .200, .183, .267, .000, .200,
     2 .500, .206, .350, .300, .241, .218, .200, .183, .350, .000, .200,
     3 .500, .297, .364, .417, .297, .288, .250, .183, .364, .000, .200,
     4 .500, .255, .315, .333, .204, .218, .183, .183, .315, .000, .200,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .200, .267, .267, .233, .300, .200, .183, .267, .000, .200,
     2 .500, .206, .350, .300, .241, .218, .200, .183, .350, .000, .200,
     3 .500, .297, .364, .417, .297, .288, .250, .183, .364, .000, .200,
     4 .500, .255, .315, .333, .204, .218, .183, .183, .315, .000, .200,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .080, .107, .107, .093, .120, .080, .073, .107, .000, .200,
     2 .500, .082, .140, .120, .096, .083, .080, .073, .140, .000, .200,
     3 .500, .119, .145, .167, .119, .115, .100, .073, .145, .000, .200,
     4 .500, .102, .126, .132, .081, .087, .073, .073, .126, .000, .200,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .100, .133, .133, .116, .150, .100, .091, .133, .000, .200,
     2 .500, .103, .175, .150, .120, .109, .100, .091, .175, .000, .200,
     3 .500, .148, .182, .208, .148, .144, .125, .091, .182, .000, .200,
     4 .500, .127, .157, .166, .102, .109, .091, .091, .157, .000, .200,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)
C        1     2     3     4     5     6     7     8     9    10    11
C      BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     1 .500, .020, .027, .027, .023, .030, .020, .018, .027, .000, .200,
     2 .500, .021, .035, .030, .024, .022, .020, .018, .035, .000, .200,
     3 .500, .030, .036, .042, .030, .029, .025, .018, .036, .000, .200,
     4 .500, .026, .032, .033, .020, .022, .018, .018, .032, .000, .200
     *     /), (/11,4,6/) )
C
!@var VTMASK vegetation depth mask by type (kg/m^2)
      real*8, parameter :: VTMASK(NV) = (/
C        1     2     3     4     5     6     7     8     9    10    11
C     BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE
     * 1d1,  2d1,  2d1,  5d1,  2d2,  5d2,  1d3, 25d2,  2d1,  1d1, .001d0
     *     /)

!@var ASHZOI,ANHZOI hemisph.Ice Albedo half-max depth (m) (orig.version)
      real*8, parameter :: ASHZOI=.1d0, ANHZOI=.1d0
!@var ASNALB snow albedo (original version)
!@var AOIALB seaice albedo (original version)
!@var ALIALB land ice albedo (original version)
C****  albedo = asnalb+snowage*snowage_fac so these numbers are the min
      real*8 ::
C                        VIS  NIR1  NIR2  NIR3  NIR4  NIR5    NIR
     *     ASNALB(7)=(/.60d0,.55d0,.55d0,.30d0,.10d0,.05d0, .35d0/),
     *     AOIALB(7)=(/.55d0,.50d0,.45d0,.25d0,.10d0,.05d0, .30d0/),
     *     ALIALB(7)=(/.60d0,.55d0,.50d0,.30d0,.10d0,.05d0, .35d0/)
C**** shorthand for the 2 band version
      real*8 ASNVIS,ASNNIR, AOIVIS,AOINIR, ALIVIS,ALINIR
      equivalence (ASNVIS,ASNALB(1)),(ASNNIR,ASNALB(7))
      equivalence (AOIVIS,AOIALB(1)),(AOINIR,AOIALB(7))
      equivalence (ALIVIS,ALIALB(1)),(ALINIR,ALIALB(7))

C**** variables that control original snow aging calculation
!@var AGEXPF exponent in snowage calculation depends on hemi/surf type
!@var ALBDIF difference in albedo as function of snowage
      REAL*8, PARAMETER, DIMENSION(3,2) ::
     *     AGEXPF = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0 /), (/3,2/) ),
     *     ALBDIF = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0/), (/3,2/) )

!@var DMOICE, DMLICE masking depth for snow on ice and land ice
      real*8, parameter :: DMOICE = 10., DMLICE = 10.

!@var AOImin seaice albedo (Hansen)
!@var AOImax seaice albedo (Hansen)
!@var ASNwet wet snow albedo (Hansen)
!@var ASNdry dry snow albedo (Hansen)
!@var AMPmin min melt pond albedo (Hansen)
      real*8, parameter ::
C                         VIS   NIR1   NIR2   NIR3   NIR4   NIR5
     *     AOImin(6)=(/ .05d0, .05d0, .05d0, .050d0, .05d0, .03d0/),
     *     AOImax(6)=(/ .62d0, .42d0, .30d0, .120d0, .05d0, .03d0/),
     *     ASNwet(6)=(/ .85d0, .75d0, .50d0, .175d0, .03d0, .01d0/),
     *     ASNdry(6)=(/ .90d0, .85d0, .65d0, .450d0, .10d0, .10d0/),
     *     AMPmin(6)=(/ .10d0, .05d0, .05d0, .050d0, .05d0, .03d0/)

!@var AOCEAN K-band dependent Thermal radiation characteristics for ocn
      real*8, parameter, dimension(NKBAND) :: AOCEAN = (/
     +        0.04000,0.09566,0.10273,0.10389,0.10464,0.10555,0.10637,
     +        0.10666,0.10697,0.10665,0.10719,0.10728,0.11007,0.04009,
     +        0.04553,0.05554,0.08178,0.09012,0.09464,0.09548,0.09532,
     +        0.09558,0.09558,0.09568,0.09565,0.05771,0.04985,0.04670,
     +        0.04630,0.04575,0.04474,0.04468,0.04500/)

!@var AGSIDV K-band dependent Thermal radiation for other types
C                     AGSNOW  AGLICE  AGROCK  AGVEG
      real*8, parameter, dimension(NKBAND,4) :: AGSIDV = RESHAPE( (/
     +        0.01400,0.09262,0.09170,0.07767,0.07130,0.06603,0.06540,
     +        0.06397,0.06358,0.06361,0.06365,0.06386,0.06564,0.01354,
     +        0.01537,0.02320,0.04156,0.03702,0.03633,0.03417,0.03346,
     +        0.03342,0.03322,0.03350,0.03170,0.01967,0.01845,0.01977,
     +        0.01986,0.01994,0.02013,0.02041,0.02100,
     +        0.01400,0.09262,0.09170,0.07767,0.07130,0.06603,0.06540,
     +        0.06397,0.06358,0.06361,0.06365,0.06386,0.06564,0.01354,
     +        0.01537,0.02320,0.04156,0.03702,0.03633,0.03417,0.03346,
     +        0.03342,0.03322,0.03350,0.03170,0.01967,0.01845,0.01977,
     +        0.01986,0.01994,0.02013,0.02041,0.02100,
     +        0.04500,0.10209,0.08806,0.05856,0.04835,0.04052,0.04001,
     +        0.03775,0.03687,0.03740,0.03637,0.03692,0.03570,0.07001,
     +        0.05665,0.05326,0.05349,0.04356,0.03845,0.03589,0.03615,
     +        0.03610,0.03602,0.03613,0.03471,0.13687,0.14927,0.16484,
     +        0.16649,0.16820,0.17199,0.17484,0.18000,
     +        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     *        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /),
     *        (/ NKBAND,4 /) )

!@var AVSCAT,ANSCAT,AVFOAM,ANFOAM for ocean albedo calc
      real*8, parameter ::
     *     AVSCAT=0.0156d0, ANSCAT=0d0, AVFOAM=0.2197d0, ANFOAM=0.1514d0

C**** miscellaneous tuning parameters
      real*8 ::
!@var WETTRA,WETSRA adjustment factors for wet earth albedo calc
     *     WETTRA=1.0, WETSRA=1.0,
!@var ZOCSRA,ZSNSRA,ZICSRA,ZDSSRA,ZVGSRA adjustment factors for
!@+   solar zenith angle effects (not fully enabled)
     *     ZOCSRA=1.0, ZSNSRA=1.0, ZICSRA=1.0, ZDSSRA=1.0, ZVGSRA=1.0,
!@var EOCTRA,ESNTRA,EICTRA,EDSTRA,EVGTRA adjustment factors for
!@+   thermal radiation effects  (not fully enabled)
     *     EOCTRA=1.0, EDSTRA=1.0, ESNTRA=1.0, EICTRA=1.0, EVGTRA=1.0

C**** ALBVNH is set only once a day then saved
!@var ALBVNH hemispherically varying vegetation albedo
      real*8, dimension(NV,6,2) :: ALBVNH

!@var GZSNOW asymmetry parameter for snow over three types
!@+   from Wiscombe and Warren (1980) JAS
!@+   Note this is used for ice + melt-ponds as well.
      REAL*8, PARAMETER, DIMENSION(7,3,2) :: GZSNOW = RESHAPE( (/
C       VIS     NIR1    NIR2     NIR3     NIR4     NIR5    NIRT
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! SH EA
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    OI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    LI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! NH EA
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,!    OI
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0 !    LI
     *     /), (/7,3,2/) )

      END MODULE SURF_ALBEDO

      MODULE RADPAR
!@sum radiation module based originally on rad00b.radcode1.F
!@auth A. Lacis/V. Oinas/R. Ruedy
      IMPLICIT NONE

C-----------------------------------------
C     Grid parameters: Vertical resolution
C-----------------------------------------

!@var LX max.number of vertical layer edges of the radiation (1D)-model
!@+
!@+   The Radiation Model can accomodate  arbitrary vertical resolution,
!@+              the number of layers may be time or location dependent,
!@+              but it cannot exceed LX-1.
!@+   Since the GCM uses 3 radiative equilibrium layers on top of the
!@+   model atmosphere, the number LM of GCM layers may be at most LX-4.
      integer, parameter :: LX = 53+4

!@var PTOPTR top of sigma layer part of vertical grid (mb), i.e.
!@+          for purposes of repartitioning the prescribed constituents,
!@+          the pressure levels above PTOPTR mb are assumed to be
!@+          fixed in time and space, below PTOPTR mb the RATIOS of
!@+          the layer thicknesses are fixed in time and space.
!@+          (altern.: Use REPART to go from mean to current P-levels)
      real*8 :: PTOPTR = 150.d0
!@var MRELAY if not 0, gases/aerosols are repartitioned to new layering
!@+   KEEP10, RO3COL, NO3COL may be used to modify this repartitioning:
!@+           if NO3COL=1, column amount of O3 is reset to RO3COL
!@+   (for OFF-line use if number of layers is time-dependent only)
      integer :: MRELAY=0, KEEP10=0, NO3COL=0 ; real*8 :: RO3COL=1.

C-------------------------------------------     This should eventually
C     Grid parameters: Horizontal resolution     be moved outside RADPAR
C-------------------------------------------     into the driver

!@var MLAT46,MLON72 horizontal grid dimensions referred to in this model
!@+   The Radiation Model utilizes Data with 72x46 (lon,lat) resolution.
!@+               For GCM resolution other than 72x46, set JLAT and ILON
!@+               to appropriately Sample  (rather than interpolate) the
!@+               72x46 aerosol, ozone, cloud heterogeneity data sets
      integer, parameter ::  MLAT46=46,MLON72=72

!@var JNORTH latitude index defining northern hemisphere : jlat>jnorth
      integer, parameter ::  JNORTH=MLAT46/2

!@var DLAT46 latitudes of box centers (degrees)
      real*8, parameter, dimension(46) :: DLAT46=(/
     A    -90.,-86.,-82.,-78.,-74.,-70.,-66.,-62.,-58.,-54.,-50.,-46.,
     B    -42.,-38.,-34.,-30.,-26.,-22.,-18.,-14.,-10., -6., -2.,  2.,
     C      6., 10., 14., 18., 22., 26., 30., 34., 38., 42., 46., 50.,
     D     54., 58., 62., 66., 70., 74., 78., 82., 86., 90./)

!@var DLON72 longitude difference to date line (degrees)
      real*8, parameter, dimension(72) :: DLON72=(/
     A      0.,  5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55.,
     B     60., 65., 70., 75., 80., 85., 90., 95.,100.,105.,110.,115.,
     C    120.,125.,130.,135.,140.,145.,150.,155.,160.,165.,170.,175.,
     D    180.,185.,190.,195.,200.,205.,210.,215.,220.,225.,230.,235.,
     E    240.,245.,250.,255.,260.,265.,270.,275.,280.,285.,290.,295.,
     F    300.,305.,310.,315.,320.,325.,330.,335.,340.,345.,350.,355./)

C----------------
C     Input data               for the 1-d radiation
C----------------

!@var LASTVC if not < zero, atmosph. and ground data are being set
      integer :: LASTVC=-123456  ! for OFF-line use only (set_io)

!@var COSZ          cosine of zenith angle  (1)
      real*8 cosz
!@var JLAT,ILON     lat,lon index  w.r.to 72x46 lon-lat grid
!@var NL,NLP        number of rad. model layers, layer edges (NLP=NL+1)
!@var LS1_loc       local tropopause level, used to limit H2O-scaling
      integer   :: JLAT,ILON,NL,NLP, LS1_loc
!@var JYEAR,JDAY    current year, Julian date
      integer :: JYEAR=1980, JDAY=1

!@var PLB           layer pressure (mb) at bottom of layer
!@var HLB           height (km) at bottom of layer - currently NOT Used
!@var TLm           mean layer temperature (K)
!@var TLb,TLt       bottom,top layer temperature (K) - computed from TLm
!@+                 except if TLGRAD<0 ; TLGRAD,PTLISO control T-profile
      real*8 :: PLB(LX),HLB(LX),TLB(LX),TLT(LX),TLM(LX)
!@var       TLGRAD     0<TLGRAD<1 controls T-profile within a layer, but
!@var       PTLISO     tlt=tlb=tlm above PTLISO mb independent of TLGRAD
      real*8 :: TLGRAD=1.,  PTLISO=2.5d0 ! control param

!@var ULGAS,U0GAS   gas amounts: curr.,ref.values (cm atm) (get/setgas)
!@var TAUWC,TAUIC   opt.depth of water,ice cloud layer (1)
!@var SIZEWC,SIZEIC particle size of water,ice clouds (micron)
!@var CLDEPS        cloud heterogeneity; is computed using KCLDEP,EPSCON
      real*8 :: ULGAS(LX,13),TAUWC(LX),TAUIC(LX),SIZEWC(LX),SIZEIC(LX)
     *     ,CLDEPS(LX)
!@var       EPSCON  cldeps=EPSCON if KCLDEP=1
!@var       KCLDEP  KCLDEP=0->CLDEPS=0, 1->=EPSCON, 2->as is, 3,4->isccp
      real*8 :: EPSCON=0. ; integer :: KCLDEP=4 ! control param

!@var SHL,RHL       layer specific,relative humidity (1)
      real*8 :: SHL(LX),RHL(LX)
!@var       KEEPRH  if 0, RHL is computed from SHL, else SHL from RHL
      integer :: KEEPRH=0       ! control param
!@var KDELIQ Flag for dry(0) or wet(1) air deliquescence
      integer :: KDELIQ(LX,4)

!@var SRBALB,SRXALB diffuse,direct surface albedo (1); see KEEPAL
      real*8 :: SRBALB(6),SRXALB(6)
!@var       KEEPAL  if 0, SRBALB,SRXALB are computed in SET/GETSUR
      integer :: KEEPAL=0       ! control param
!@dbparm    KSIALB  sea ice albedo computation flag: 0=Hansen 1=Lacis
      integer :: KSIALB=0
!@var PVT           frac. of surf.type (bareWhite+veg*8+bareDark+ocn)(1)
!@var AGESN 1-3     age of snow    (over soil,oice,land ice) (days)
!@var SNOWE,SNOWLI  amount of snow (over soil,land ice)   (kg/m^2)
!@var SNOWOI        amount of snow (over ocean/lake ice)  (kg/m^2)
!@var WEARTH        soil wetness (1)
!@var WMAG          wind speed (m/s)
!@var POCEAN        fraction of box covered by ocean or lake  (1)
!@var PLAKE         fraction of box covered by lake           (1)
!@var PEARTH        fraction of box covered by soil           (1)
!@var POICE         fraction of box covered by ocean/lakeice  (1)
!@var PLICE         fraction of box covered by glacial ice    (1)
!@var TGO           top layer water temperature (K) of ocean/lake
!@var TGE,TGOI,TGLI top layer ground temperature (K) soil,seaice,landice
!@var TSL           surface air temperature (K)
      real*8 PVT(11),AGESN(3),SNOWE,SNOWOI,SNOWLI,WEARTH,WMAG,POCEAN
     *     ,PEARTH,POICE,PLICE,PLAKE,TGO,TGE,TGOI,TGLI,TSL
!@var KZSNOW        =1 for snow/ice albedo zenith angle dependence
      integer :: KZSNOW=1
!     Additional info for Schramm/Schmidt/Hansen sea ice albedo KSIALB=0
!@var ZSNWOI        depth of snow over ocean ice (m)
!@var zoice         depth of ocean ice (m)
!@var zmp           depth of melt pond (m)
!@var fmp           fraction of melt pond area (1)
!@var zlake         lake depth (m)
!@var flags         true if snow is wet
!@var snow_frac(2)  fraction of snow over bare(1),vegetated(2) soil (1)
!@var snoage_fac_max  max snow age reducing-factor for sea ice albedo
      REAL*8 :: zsnwoi,zoice,zmp,fmp,zlake,snow_frac(2)
      REAL*8 :: snoage_fac_max=.5d0

!@var TRACER array to add up to 8 additional aerosol species
      real*8    :: TRACER(LX,8)
!@var FSTOPX,FTTOPX scales optional aerosols (solar,thermal component)
      real*8    :: FSTOPX(8),FTTOPX(8)
!@var O3_IN column variable for importing ozone field from rest of model
!@var use_tracer_ozone =0 normal case, =1 means
!@+   RCOMPX will use O3_IN(L) for U0GAS(L,3) for GCM levels
      real*8, dimension(lx) :: O3_IN
      integer use_tracer_ozone
      LOGICAL*4 :: flags

      COMMON/RADPAR_INPUT_IJDATA/    !              Input data to RCOMPX
     A              PLB,HLB,TLB,TLT,TLM,ULGAS
     B             ,TAUWC,TAUIC,SIZEWC,SIZEIC,CLDEPS
     C             ,SHL,RHL,TRACER,SRBALB,SRXALB
     D             ,PVT,AGESN,SNOWE,SNOWOI,SNOWLI,WEARTH,WMAG
     E             ,POCEAN,PEARTH,POICE,PLICE,PLAKE
     F             ,TGO,TGE,TGOI,TGLI,TSL,COSZ,FSTOPX,FTTOPX,O3_IN
     X             ,zsnwoi,zoice,zmp,fmp,snow_frac,zlake
C      integer variables start here, followed by logicals
     Y             ,JLAT,ILON,NL,NLP, LS1_loc,flags,use_tracer_ozone
     Z             ,KDELIQ                ! is updated by rad. after use
!$OMP  THREADPRIVATE(/RADPAR_INPUT_IJDATA/)

C     array with local and global entries: repeat this section in driver
      REAL*8 U0GAS(LX,13)
      COMMON/RADPAR_hybrid/U0GAS
!$OMP THREADPRIVATE(/RADPAR_hybrid/)
C     end of section to be repeated in driver (needed for 'copyin')

C--------------------------------------------------------
C     Output data     (from RCOMPX)  grid point dependent
C--------------------------------------------------------

!@var TRDFLB,TRUFLB,TRNFLB  thrml down,up,net flux at layr bottom (W/m2)
!@var SRDFLB,SRUFLB,SRNFLB  solar down,up,net flux at layr bottom (W/m2)
!@var TRFCRL,SRFHRL         layer LW cooling rate,SW heating rate (W/m2)
!@var etc etc
!@var etc etc
!sl!@var FTAUSL,TAUSL,...  surface layer computations commented out: !sl
!@var LBOTCL,LTOPCL  bottom and top cloud level (lbot < ltop)
!@var O3_OUT column variable for exporting ozone field to rest of model

      real*8, dimension(lx) :: TRDFLB,TRUFLB,TRNFLB,TRFCRL,
     *     SRDFLB,SRUFLB,SRNFLB,SRFHRL,O3_OUT
      real*8 :: SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,
     *     SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,
     *     SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,
     *     TRDFGW,TRUFGW,TRUFTW,BTEMPW,SRXVIS,SRXNIR
      real*8, dimension(4) :: FSRNFG,FTRUFG,DTRUFG !,SRXATM
      real*8 :: WINDZF(3),WINDZT(3),TOTLZF(3),TOTLZT(3),SRKINC(16),
     *     SRKALB(16),SRKGAX(16,4),SRKGAD(16,4),SKDFLB(LX,17),
     *     SKUFLB(LX,17),SKNFLB(LX,17),SKFHRL(LX,17)
!sl  K             ,FTAUSL(33),TAUSL(33)    ! input rather than output ?
!nu  K             ,TRDFSL,TRUFSL,TRSLCR,SRSLHR,TRSLWV   !nu = not used
!sl  K             ,TRSLTS,TRSLTG,TRSLBS

      integer :: LBOTCL,LTOPCL

      COMMON/RADPAR_OUTPUT_IJDATA/
     A              TRDFLB,TRUFLB,TRNFLB,TRFCRL
     B             ,SRDFLB,SRUFLB,SRNFLB,SRFHRL
     C             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR  !,SRXATM
     D             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG
     E             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG
     F             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,DTRUFG
     G             ,WINDZF,WINDZT,TOTLZF,TOTLZT,SRKINC
     I             ,SRKALB,SRKGAX,SRKGAD,SKDFLB
     J             ,SKUFLB,SKNFLB,SKFHRL,SRXVIS,SRXNIR
!sl  K             ,FTAUSL,TAUSL    ! input rather than output ?
!nu  K             ,TRDFSL,TRUFSL,TRSLCR,SRSLHR,TRSLWV   !nu = not used
!sl  K             ,TRSLTS,TRSLTG,TRSLBS
     K             ,O3_OUT
     L             ,LBOTCL,LTOPCL   ! integers last for alignment
!$OMP THREADPRIVATE(/RADPAR_OUTPUT_IJDATA/)
!nu   EQUIVALENCE (SRXATM(1),SRXVIS),(SRXATM(2),SRXNIR)
!nu   EQUIVALENCE (SRXATM(3),XXAVIS),(SRXATM(4),XXANIR)  !nu = not used

C----------------   scratch pad for temporary arrays that are passed to
C     Work arrays   other routines while working on a lat/lon point; but
C----------------   in multi-cpu mode, each cpu needs its own copy !!

      real*8, dimension(LX,6) ::
     *     SRAEXT,SRASCT,SRAGCB,SRBEXT,SRBSCT,SRBGCB,
     *     SRDEXT,SRDSCT,SRDGCB,SRVEXT,SRVSCT,SRVGCB,
     *     SRCEXT,SRCSCT,SRCGCB,DBLEXT,DBLSCT,DBLGCB,DBLPI0,SRCPI0,
     *     QDUST,SDUST,CDUST                    !nu ,QAERO,SAERO,CAERO
      real*8, dimension(LX,33) :: TRCALK,TRAALK,TRBALK,TRTAUK,TRDALK
     *     ,TRVALK,TRGXLK,DFLB,UFLB,WFLB,ADUST  !nu ,AAERO
      real*8, dimension(33) :: TRCTCA,DFSL,UFSL,WFSL,CLPI0,TXCTPG,TSCTPG
     *     ,TGCTPG,AVH2S,TRGALB,BGFEMT,BGFEMD
      real*8, dimension(LX) :: PL,DPL,WTLB,WTLT
     *     ,ENA,ENB,ENC,TRA,TRB,TRC,AERX1,AERS1,AERG1,TRAXNL,AERX2,AERS2
     *     ,AERG2,UGAS0,UGASR,RNB,RNX,TNB,TNX,XNB,XNX,SRB,SRX,VRU,VRD
     *     ,FAC,DNA,DNB,DNC,O2FHRL,SRAXNL,SRASNL,SRAGNL,AO3X,O2FHRB,AO3D
     *     ,AO3U,HTPROF
      real*8 :: UXGAS(LX,9),TAUN(33*LX),BXA(7),PRNB(6,4),PRNX(6,4)
     *     ,Q55H2S,RIJTCK(6,33),FDXTCK(3,33),ALBTCK(3,33),FEMTCK(3,33)
     *     ,QVH2S(6),SVH2S(6),GVH2S(6),XTRU(LX,4),XTRD(LX,4)
      integer, dimension(LX) :: ITLB,ITLT
C**** local except for special radiative aerosol diagnostics aadiag
      real*8  ATAULX(LX,6)
      integer NRHNAN(LX,8)

      COMMON/WORKDATA/          !          Temp data generated by RCOMPX
     A              SRAEXT,SRASCT,SRAGCB,TRCALK
     B             ,SRBEXT,SRBSCT,SRBGCB,TRAALK
     C             ,SRDEXT,SRDSCT,SRDGCB,TRBALK
     D             ,SRVEXT,SRVSCT,SRVGCB,TRTAUK
     E             ,DBLEXT,DBLSCT,DBLGCB,DBLPI0
     F             ,SRCEXT,SRCSCT,SRCGCB,SRCPI0
     G             ,TRDALK,TRVALK,TRGXLK,TRCTCA
     H             ,DFLB,UFLB,WFLB,PL,DPL
     I             ,TRGALB,BGFEMT,BGFEMD,WTLB,WTLT
     J             ,ENA,ENB,ENC,TRA,TRB,TRC
     K             ,DFSL,UFSL,WFSL
     L             ,AERX1,AERS1,AERG1,TRAXNL
     M             ,AERX2,AERS2,AERG2,UGAS0,UGASR
     N             ,RNB,RNX,TNB,TNX,XNB,XNX
     O             ,SRB,SRX,VRU,VRD,FAC
     P             ,UXGAS,TAUN,BXA,PRNB,PRNX
     Q             ,DNA,DNB,DNC,Q55H2S
     R             ,RIJTCK,FDXTCK,ALBTCK,CLPI0
     S             ,FEMTCK,TXCTPG,TSCTPG,TGCTPG
!nu  T             ,QAERO,SAERO,CAERO,AAERO
     U             ,QDUST,SDUST,CDUST,ADUST
     V             ,O2FHRL,SRAXNL,SRASNL,SRAGNL,AO3X
     W             ,O2FHRB,AO3D,AO3U
     X             ,HTPROF,QVH2S,SVH2S,GVH2S,AVH2S
     F             ,XTRU,XTRD,                  ATAULX
     I             ,ITLB,ITLT,                  NRHNAN   ! integers last
!$OMP  THREADPRIVATE(/WORKDATA/)

      real*8 ::                      !  Temp data used by WRITER, WRITET
     A           SRCQPI(6,15),TRCQPI(33,15),
     B           TRAQAB(33,11),TRBQAB(33,10),TRCQAB(33,15),TRDQAB(33,25)
      integer :: NORDER(16),NMWAVA(16),NMWAVB(16)

C------------------------------------------
C     Reference data, Tables, Climatologies
C------------------------------------------

      real*8, parameter, dimension(16) :: DKS0=(/
     *           .010, .030, .040, .040, .040, .002, .004, .013,
     +           .002, .003, .003, .072, .200, .480, .050, .011/)

      integer ::  NKSLAM=14
      integer,parameter,dimension(16) ::
     *     KSLAM=(/1,1,2,2,5,5,5,5,1,1,1,3,4,6,6,1/)

      real*8 ::                 !   Model parameters generated by RCOMP1
     H              HLB0(LX),PLB0(LX),TLM0(LX),U0GAS3(LX)
     A             ,TKPFW(630),TKPFT(900),AO3(460)
     D             ,FPXCO2(LX),FPXOZO(LX) !nu ,PIAERO(10)
     C             ,SRAX(LX,6,5),SRAS(LX,6,5),SRAC(LX,6,5),ZTABLE(LX,11)
     E             ,QXDUST(6,8),QSDUST(6,8),QCDUST(6,8),ATDUST(33,8)
     D             ,QDST55(8),TRAX(LX,33,5),DBLN(30),TCLMIN

C            RADDAT_TR_SGP_TABLES          read from  radfile1, radfile2
      real*8 ::
     A              TAUTBL(154000),PLANCK(8250),XKCFC(12,8,4)
     B             ,TAUWV0(154000),H2OCN8(33,8,14),H2OCF8(33,8,5)
     D             ,ULOX(285),DUX(285),GTAU(51,11,143),TGDATA(122,13)
     E           ,SALBTG(768,14),TAUGSA(1001,14),TAUTGD(122),TAUTGS(768)
     F            ,XUCH4(9,15),XUN2O(9,15),XTRUP(24,3,15),XTRDN(24,3,15)
     G             ,CXUO3(7,15),CXUCO2(7,15),XTU0(24,3)
     H             ,XTD0(24,3),XUCH40(9),XUN2O0(9)
      integer, parameter :: ITPFT0=123 ,ITNEXT=250 ! offsets for lookups

C            RADDAT_AERCLD_MIEPAR          read from            radfile3
      real*8 ::
     A              SRAQEX( 6,11),SRAQSC( 6,11),SRAQCB( 6,11),Q55A11(11)
     B             ,TRAQEX(33,11),TRAQSC(33,11),TRAQCB(33,11),REFA11(11)
     C             ,SRBQEX( 6,10),SRBQSC( 6,10),SRBQCB( 6,10),Q55B10(10)
     D             ,TRBQEX(33,10),TRBQSC(33,10),TRBQCB(33,10),REFB10(10)
     E             ,SRCQEX( 6,15),SRCQSC( 6,15),SRCQCB( 6,15),Q55C15(15)
     F             ,TRCQEX(33,15),TRCQSC(33,15),TRCQCB(33,15),REFC15(15)
     G             ,TRCQAL(33,15),VEFC15(15)   ,VEFA11(   11),VEFB10(10)
     H             ,SRDQEX( 6,25),SRDQSC( 6,25),SRDQCB( 6,25),Q55D25(25)
     I             ,TRDQEX(33,25),TRDQSC(33,25),TRDQCB(33,25),REFD25(25)
     J             ,TRDQAL(33,25),VEFD25(25)
     K         ,SRVQEX( 6,20,6),SRVQSC( 6,20,6),SRVQCB( 6,20,6)
     L         ,TRVQEX(33,20,6),TRVQSC(33,20,6),TRVQCB(33,20,6)
     M         ,TRVQAL(33,20,6),Q55V20(20,6),REFV20(20,6),VEFV20(20,6)
     N         ,SRUQEX( 6,120),SRUQSC( 6,120),SRUQCB( 6,120),Q55U22(120)
     O         ,TRUQEX(33,120),TRUQSC(33,120),TRUQCB(33,120),REFU22(120)
     P         ,TRUQAL(33,120),VEFU22(120),TRSQAL(33,25),VEFS25(25)
     Q             ,SRSQEX( 6,25),SRSQSC( 6,25),SRSQCB( 6,25),Q55S25(25)
     R             ,TRSQEX(33,25),TRSQSC(33,25),TRSQCB(33,25),REFS25(25)

      real*8    SRQV( 6,20),SRSV( 6,20),SRGV( 6,20),Q55V(   20),REFV(20)
      real*8    TRQV(33,20),TRSV(33,20),TRGV(33,20),TRAV(33,20),VEFV(20)
      EQUIVALENCE (SRVQEX(1,1,6),SRQV(1,1)), (SRVQSC(1,1,6),SRSV(1,1))
      EQUIVALENCE (SRVQCB(1,1,6),SRGV(1,1)),   (Q55V20(1,6),Q55V(1))
      EQUIVALENCE (TRVQEX(1,1,6),TRQV(1,1)), (TRVQSC(1,1,6),TRSV(1,1))
      EQUIVALENCE (TRVQCB(1,1,6),TRGV(1,1)), (TRVQAL(1,1,6),TRAV(1,1))
      EQUIVALENCE   (REFV20(1,6),REFV(1)),     (VEFV20(1,6),VEFV(1))

C            RADDAT_CLDCOR_TRSCAT           read from           radfileE
      real*8 :: RIJTPG(6,49,17,21),FDXTPG(3,49,17,21),FEMTPG(3,49,17,21)



C--------------------------------------    This also should be moved out
C     History files (+ control options)    of RADPAR, which should just
C--------------------------------------    have to deal  1 point in time

!     -------------------------------------------------------i/o control
!@var MADxxx  Model Add-on Data of Extended Climatology Enable Parameter
!@+   ------   if 0   input process is skipped
!@+ 2 MADAER   =  1   Reads  Aerosol tropospheric climatology
!@+ 3 MADDST   =  1   Reads  Dust-windblown mineral climatology   RFILE6
!@+ 4 MADVOL   =  1   Reads  Volcanic 1950-00 aerosol climatology RFILE7
!@+ 5 MADEPS   =  1   Reads  Epsilon cloud heterogeneity data     RFILE8
!@+ 6 MADLUV   =  1   Reads  Lean's SolarUV 1882-1998 variability RFILE9
!@+   MADGHG   =  1          Enables UPDGHG update. MADGHG=0: no update
!@+   MADSUR   =  1   Reads  Vegetation,Topography data    RFILEC,RFILED
!@+   MADBAK   if 1          Adds background aerosols
!     ------------------------------------------------------------------
      integer :: MADO3M=1,MADAER=1,MADDST=1,MADVOL=1,MADEPS=1,MADLUV=1
      integer :: MADGHG=1,MADSUR=0,MADBAK=0 ! MADSUR=1 for OFF-line only

!     ------------------------------------------------------time control
!@var KYEARx,KJDAYx if both are 0   : data are updated to current yr/day
!@+   -------------    only KJDAYx=0: data cycle through year KYEARx
!@+                    neither is 0 : yr/day=KYEARx/KJDAYx data are used
!@+   KYEARS,KJDAYS: Solar Trend
!@+   KYEARO,KJDAYO: Ozone Trend
!@+   KYEARD,KJDAYD: Dust Trend
!@+   KYEARE,KJDAYE: CldEps Trend
!@+   KYEARG,KJDAYG: GHG  Trend
!@+   KYEARR,KJDAYR: RVegeTrend (Ground Albedo)
!@+   KYEARV,KJDAYV: Volc.Aerosol Trend
!@+   KYEARA,KJDAYA: trop.Aerosol Trend
!     ------------------------------------------------------------------
      integer ::                    KYEARS=0,KJDAYS=0, KYEARG=0,KJDAYG=0
     *          ,KYEARO=0,KJDAYO=0, KYEARA=0,KJDAYA=0, KYEARD=0,KJDAYD=0
     *          ,KYEARV=0,KJDAYV=0, KYEARE=0,KJDAYE=0, KYEARR=0,KJDAYR=0

      INTEGER, PARAMETER :: NLO3=49 !  # of layers in ozone data files
      real*8 :: O3JDAY(NLO3,MLON72,MLAT46)
      COMMON/O3JCOM/O3JDAY
C**** PLBO3(NLO3+1) could be read off the titles of the decadal files
      REAL*8, PARAMETER, DIMENSION(NLO3+1) :: PLBO3 = (/
     *       984d0, 934d0, 854d0, 720d0, 550d0, 390d0, 285d0, 210d0, 
     *       150d0, 125d0, 100d0,  80d0,  60d0,  55d0,  50d0, 
     *        45d0,  40d0,  35d0,  30d0,  25d0,  20d0,  15d0, 
     *       10.d0,  7.d0,  5.d0,  4.d0,  3.d0,  2.d0,  1.5d0, 
     *        1.d0,  7d-1,  5d-1,  4d-1,  3d-1,  2d-1,  1.5d-1, 
     *        1d-1,  7d-2,  5d-2,  4d-2,  3d-2,  2d-2,  1.5d-2, 
     *        1d-2,  7d-3,  5d-3,  4d-3,  3d-3,  1d-3,  1d-7/)

!@var PLBA09 Vert. Layering for tropospheric aerosols/dust (reference)
      real*8, parameter, dimension(10) :: PLBA09=(/
     *   984.,934.,854.,720.,550.,390.,255.,150., 70., 10./)
C     Layer  1    2    3    4    5    6    7    8    9

C            RADMAD3_DUST_SEASONAL            (user SETDST)     radfile6
      real*4 TDUST(72,46,9,8,13)
      real*8 DDJDAY(9,8,72,46)

C            RADMAD4_VOLCAER_DECADAL          (user SETVOL)     radfile7
      real*8         V4TAUR(1800,24,5),FDATA(80),GDATA(80)
     C              ,HTFLAT(49,4),TAULAT(49),SIZLAT(49)

C            RADMAD5_CLDEPS_3D_SEASONAL       (user SETCLD)     radfile8
      real*4 EPLMHC(72,46,12,4)
      real*8 EPLOW(72,46),EPMID(72,46),EPHIG(72,46),EPCOL(72,46)

C            RADMAD6_SOLARUV_DECADAL          (user SETSOL)     radfile9
!@var iy1S0,MS0X first year, max.number of months for S0 history
!@var icycs0,mcycs0 solar cycle in yrs,months used to extend S0 history
!@var KSOLAR controls which data are used: <0 Thekaekara, else Lean:
!@+          1: use monthly data, 2: use annual data, 0: constant data
      integer :: KSOLAR=1       ! MADLUV=KSOLAR=0 only possible OFF-line

      integer, parameter :: iy1S0=1882, MS0X=12*(1998-iy1S0+1)
      integer, parameter :: icycs0=11,  mcycs0=icycs0*12
      real*4 UVLEAN(Ms0X,190),yr1S0,yr2S0
      real*8 TSI1(Ms0X),TSI2(Ms0X),FSLEAN(190),W1LEAN(190)

      real*8 :: S00WM2=1366.2911d0, S0=1366.d0, RATLS0=1.

      real*8 :: WSOLAR(190),FSOLAR(190)

C***  alternate sources to get WSOLAR,FSOLAR:
      real*8, dimension(190) :: WSLEAN,DSLEAN,FRLEAN
      common/LEAN1950/   WSLEAN,DSLEAN,FRLEAN              ! if MADLUV=0

      real*8, parameter, dimension(190) :: WTHEK=(/        ! if KSOLAR<0
     *           .115,.120,.125,.130,.140,.150,.160,.170,.180,.190,.200,
     1 .210,.220,.225,.230,.235,.240,.245,.250,.255,.260,.265,.270,.275,
     2      .280,.285,.290,.295,.300,.305,.310,.315,.320,.325,.330,.335,
     3           .340,.345,.350,.355,.360,.365,.370,.375,.380,.385,.390,
     4           .395,.400,.405,.410,.415,.420,.425,.430,.435,.440,.445,
     5           .450,.455,.460,.465,.470,.475,.480,.485,.490,.495,.500,
     6           .505,.510,.515,.520,.525,.530,.535,.540,.545,.550,.555,
     7           .560,.565,.570,.575,.580,.585,.590,.595,.600,.605,.610,
     8           .620,.630,.640,.650,.660,.670,.680,.690,.700,.710,.720,
     9           .730,.740,.750,.760,.770,.780,.790,.800,.810,.820,.830,
     A .840,.850,.860,.870,.880,.890,.900,.910,.920,.930,.940,.950,.960,
     B 0.97,0.98,0.99,1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,
     C 1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00,2.10,2.20,
     D 2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,3.50,
     E 3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,4.50,4.60,4.70,4.80,
     F  4.9, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,12.0,13.0,14.0,15.00/)

      real*8, parameter, dimension(190) :: FTHEK=(/
     *         .007,.900,.007,.007,.030,.070,.230,.630,1.25,2.71,10.7,
     1 22.9,57.5,64.9,66.7,59.3,63.0,72.3,70.4,104.,130.,185.,232.,204.,
     2    222.,315.,482.,584.,514.,603.,689.,764.,830.,975.,1059.,1081.,
     31074.,1069.,1093.,1083.,1068.,1132.,1181.,1157.,1120.,1098.,1098.,
     41189.,1429.,1644.,1751.,1774.,1747.,1693.,1639.,1663.,1810.,1922.,
     52006.,2057.,2066.,2048.,2033.,2044.,2074.,1976.,1950.,1960.,1942.,
     61920.,1882.,1833.,1833.,1852.,1842.,1818.,1783.,1754.,1725.,1720.,
     71695.,1705.,1712.,1719.,1715.,1712.,1700.,1682.,1666.,1647.,1635.,
     81602.,1570.,1544.,1511.,1486.,1456.,1427.,1402.,1389.,1344.,1314.,
     91290.,1260.,1235.,1211.,1185.,1159.,1134.,1109.,1085.,1060.,1036.,
     A1013.,990.,968.,947.,926.,908.,891.,880.,869.,858.,847.,837.,820.,
     B 803.,785.,767.,748.,668.,593.,535.,485.,438.,397.,358.,337.,312.,
     C 288.,267.,245.,223.,202.,180.,159.,142.,126.,114.,103., 90., 79.,
     D 69.0,62.0,55.0,48.0,43.0,39.0,35.0,31.0,26.0,22.6,19.2,16.6,14.6,
     E 13.5,12.3,11.1,10.3, 9.5,8.70,7.80,7.10,6.50,5.92,5.35,4.86,4.47,
     F  4.11,3.79,1.82,0.99,.585,.367,.241,.165,.117,.0851,.0634,.0481/)

!icb         RADMAD7_VEG_TOPOG          (user SETSUR)  radfileC,radfileD
!icb                 FVEG11(72,46,11),FOLGIZ(72,46,9)

C            RADMAD8_RELHUM_AERDATA     (user SETAER,SETREL)    radfileH
!@var KRHAER(4) 0/1=off/on flag to make aeros.sizes humidity dependent
      integer, dimension(4) :: KRHAER=(/1,1,1,1/) ! SO4,SeaSalt,Nitr,org
!@var KRHTRA(8) 0/1=off/on flag to make tracer aeros.sizes humidity dep.
      integer, dimension(8) :: KRHTRA=(/1,1,1,1,1,1,1,1/)
      real*8 ::
     A               SRHQEX(6,190,4),SRHQSC(6,190,4),SRHQCB( 6,190,4)
     B              ,TRHQAB(33,190,4),RHINFO(190,9,4),A6JDAY(9,6,72,46)
     C              ,SRTQEX(6,190,8),SRTQSC(6,190,8),SRTQCB( 6,190,8)
     D              ,TRTQAB(33,190,8),RTINFO(190,9,8)
!new
!new  save TSOIL,TVEGE                  (not implemented)
!nu   DIMENSION PI0TRA(11)
!new  save FTRUFS,FTRUFV,DTRUFS,DTRUFV  (not implemented)

C     -----------------------
C     Ozone absorption tables
C     -----------------------
      real*8, parameter, dimension(226) ::        XWAVO3=(/
     *            .2002,.2012,.2022,.2032,.2042,.2052,.2062,.2072,.2082,
     A.2092,.2102,.2112,.2122,.2132,.2142,.2152,.2162,.2172,.2182,.2192,
     B.2202,.2212,.2222,.2232,.2242,.2252,.2262,.2272,.2282,.2292,.2302,
     C.2312,.2322,.2332,.2342,.2352,.2362,.2372,.2382,.2392,.2400,.2402,
     D.2412,.2422,.2432,.2438,.2444,.2452,.2458,.2463,.2472,.2478,.2482,
     E.2490,.2492,.2500,.2508,.2519,.2527,.2539,.2543,.2553,.2562,.2566,
     F.2571,.2575,.2579,.2587,.2597,.2604,.2617,.2624,.2635,.2643,.2650,
     G.2654,.2662,.2669,.2675,.2682,.2692,.2695,.2702,.2712,.2718,.2722,
     H.2732,.2742,.2746,.2752,.2762,.2772,.2782,.2792,.2802,.2812,.2822,
     I.2830,.2842,.2852,.2862,.2872,.2882,.2892,.2902,.2912,.2922,.2932,
     J.2942,.2952,.2962,.2972,.2982,.2992,.2998,
     &            .3004,.3016,.3021,.3029,.3036,.3037,.3051,.3053,.3059,
     A.3061,.3066,.3075,.3077,.3083,.3085,.3092,.3098,.3100,.3104,.3106,
     B.3109,.3112,.3130,.3135,.3146,.3148,.3151,.3154,.3167,.3170,.3173,
     C.3176,.3190,.3194,.3199,.3200,.3209,.3210,.3216,.3220,.3223,.3226,
     D.3239,.3242,.3245,.3248,.3253,.3255,.3269,.3272,.3275,.3279,.3292,
     E.3295,.3299,.3303,.3309,.3312,.3328,.3332,.3334,.3338,.3357,.3365,
     F.3369,.3372,.3391,.3395,.3398,.3401,.3417,.3421,.3426,.3430,.3437,
     G.3439,.3451,.3455,.3460,.3463,.3466,.3472,.3481,.3485,.3489,.3493,
     H.3499,.3501,.3506,.3514,.3521,.3523,.3546,.3550,.3554,.3556,.3561,
     I.3567,.3572,.3573,.3588,.3594,.3599,.3600,.3604,.3606,.3639,.3647,
     J.3650,.3654,.3660/)
      real*8, dimension(226) ::  UVA
      real*8, parameter, dimension(226) ::  FUVKO3=(/
     *             8.3,  8.3,  8.1,  8.3,  8.6,  9.0,  9.7, 10.8, 11.7,
     A 13.0, 14.3, 16.0, 18.0, 20.6, 23.0, 26.1, 29.3, 32.6, 36.9, 40.8,
     B 46.9, 51.4, 56.7, 63.4, 69.1, 76.6, 84.0, 91.4, 99.9,110.0,118.0,
     C126.0,136.0,145.0,154.0,164.0,175.0,186.0,192.0,201.0,210.0,212.0,
     D221.0,230.0,239.0,248.0,250.0,259.0,264.0,264.0,273.0,277.0,275.0,
     E283.0,283.0,290.0,283.0,297.0,290.0,300.0,290.0,302.0,295.0,283.0,
     F293.0,290.0,286.0,297.0,281.0,280.0,271.0,275.0,254.0,264.0,250.0,
     G248.0,242.0,228.0,230.0,216.0,213.0,211.0,199.0,188.0,188.0,178.0,
     H169.0,153.0,155.0,148.0,136.0,127.0,117.0,108.0, 97.0, 88.7, 81.3,
     I 78.7, 67.9, 61.4, 54.3, 49.6, 43.1, 38.9, 34.6, 30.2, 27.5, 23.9,
     J 21.0, 18.6, 16.2, 14.2, 12.3, 10.7,  9.5,
     &            8.880,7.520,6.960,6.160,5.810,5.910,4.310,4.430,4.130,
     A4.310,4.020,3.330,3.390,3.060,3.100,2.830,2.400,2.490,2.330,2.320,
     B2.120,2.200,1.436,1.595,1.074,1.138,1.068,1.262,0.818,0.948,0.860,
     C1.001,0.543,0.763,0.665,0.781,0.382,0.406,0.373,0.608,0.484,0.601,
     D0.209,0.276,0.259,0.470,0.319,0.354,0.131,0.223,0.185,0.339,0.080,
     E0.093,0.079,0.184,0.139,0.214,0.053,0.074,0.068,0.152,0.038,0.070,
     F.0540000,.1030000,.0240000,.0382500,.0292500,.0550000,.0135000,
     G.0155250,.0127500,.0188250,.0167250,.0262500,.0115500,.0140250,
     H.0099750,.0115500,.0081000,.0104250,.0050100,.0057000,.0046650,
     I.0073425,.0051825,.0055275,.0040575,.0077700,.0048900,.0054600,
     J.0015375,.0017775,.0013275,.0014100,.0011550,.0023325,.0018825,
     K.0019650,.0009600,.0013650,.0011925,.0013200,.0008925,.0009825,
     L.0001350,.0006300,.0004500,.0006225,0.0/)

C     ------------------------------------------------------------------
C          NO2 Trace Gas Vertical Distribution and Concentration Profile
C     ------------------------------------------------------------------

      real*8, parameter, dimension(42) ::
     *     CMANO2=(/            ! every 2 km starting at 0km
     1  8.66E-06,5.15E-06,2.85E-06,1.50E-06,9.89E-07,6.91E-07,7.17E-07,
     2  8.96E-07,3.67E-06,4.85E-06,5.82E-06,6.72E-06,7.77E-06,8.63E-06,
     3  8.77E-06,8.14E-06,6.91E-06,5.45E-06,4.00E-06,2.67E-06,1.60E-06,
     4  8.36E-07,3.81E-07,1.58E-07,6.35E-08,2.57E-08,1.03E-08,4.18E-09,
     5  1.66E-09,6.57E-10,2.58E-10,1.02E-10,4.11E-11,1.71E-11,7.73E-12,
     6  9.07E-12,4.63E-12,2.66E-12,1.73E-12,1.28E-12,1.02E-12,1.00E-30/)

C     ------------------------------------------------------------------
C     TRACE GAS REFERENCE AMOUNTS & DISTRIBUTIONS ARE DEFINED IN  SETGAS
C     ------------------------------------------------------------------

C-------------------------
C     Scaling/kill factors
C-------------------------

!@var FULGAS scales the various atmospheric constituents:
!@+         H2O CO2 O3 O2 NO2 N2O CH4 F11 F12 N2C CFC11 CFC12 SO2
!@+   Note: FULGAS(1) only acts in the stratosphere (unless LS1_loc=1)
      real*8, dimension(13) :: FULGAS = (/    ! scales ULGAS

C      H2O CO2  O3  O2 NO2 N2O CH4 F11 F12 N2C CFC11+ CFC12+ SO2
C        1   2   3   4   5   6   7   8   9  10    11     12   13
     +   1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,   1.,    1.,  0./)

!@var FGOLDH scales background aerosols for Glb Ocn Land Desert Haze
      real*8, dimension(5) ::       ! used by setbak/getbak only
C               GLOBAL  OCEAN   LAND  DESERT    HAZE
C                    1      2      3       4       5
     +   FGOLDH=(/ 1d0, .68d0, .32d0, 1.d-20, 1.d-20 /)

!@var FSXAER,FTXAER scales solar,thermal opt.depth for var. aerosols:
!@+          1: total 2:background 3: AClim 4:dust 5:volcanic
      real*8, dimension(5) ::  FSXAER=(/1.,1.,1.,1.,1./)
      real*8, dimension(5) ::  FTXAER=(/1.,1.,1.,1.,1./)
      real*8 FSTAER,FSBAER,FSAAER,FSDAER,FSVAER
     *      ,FTTAER,FTBAER,FTAAER,FTDAER,FTVAER
      EQUIVALENCE (FSXAER(1),FSTAER),  (FTXAER(1),FTTAER)
      EQUIVALENCE (FSXAER(2),FSBAER),  (FTXAER(2),FTBAER)
      EQUIVALENCE (FSXAER(3),FSAAER),  (FTXAER(3),FTAAER)
      EQUIVALENCE (FSXAER(4),FSDAER),  (FTXAER(4),FTDAER)
      EQUIVALENCE (FSXAER(5),FSVAER),  (FTXAER(5),FTVAER)

!@var PIVMAX limits PI0 of volcanic aerosols
      real*8 :: PIVMAX=1.0
!@var ECLTRA,KCLDEM scales,enables full cloud scattering correction
      real*8 :: ECLTRA=1. ; integer :: KCLDEM=1
!@var FCLDTR,FCLDSR scales opt.depth of clouds - not used (yet)
!@var FRAYLE        scales Rayleigh parameter
      real*8 ::   FCLDTR=1.,  FCLDSR=1.,  FRAYLE=1.

!@var KUVFAC,UVFACT,UVWAVL,KSNORM rescale UV spectral flux distribution
      integer :: KUVFAC=0,  KSNORM=0  ! no rescaling
      real*8, dimension(3)  :: UVWAVL=(/0.0010,0.0020,0.0030/)
      real*8, dimension(3)  :: UVFACT=(/1.0000,1.0000,1.0000/)

!@var SRCGSF Scaling Factors for Cloud Asymmetry Parameter for
!@+                                      Water    Ice    MieIce
      real*8, dimension(3) ::  SRCGSF=(/ 1.000,  1.000,  1.000/)

!@var TAUWC0,TAUIC0 lower limits for water/ice cloud opt.depths
      real*8 ::  TAUWC0=1d-3, TAUIC0=1d-3

!@var KPFCO2,KPFOZO if > 0 scale CO2,O3 to stand. vertical profile
      integer :: KPFCO2=0,  KPFOZO=0

!@var KANORM,KCNORM if > 0 renormalize aerosols,cloud albedos
      integer :: KANORM=0, KCNORM=0

!@var KWVCON        ON/OFF flag for water vapor continuum absorption
!@var KUFH2O,KUFCO2 H2O,CO2 column absorb.scaling
!@var KCSELF,KCFORN H2O_ContSelf-Broadening,CO2_ContForeign-Broadening
      integer :: KWVCON=1, KUFH2O=1,  KUFCO2=1,  KCSELF=1,  KCFORN=1

!@var ICE012 pick ice droplet type: 0 liquid, 1 ice non-spher, 2 ice Mie
      integer :: ICE012=1

!@var VEFF0 effective volc. aerosol size distribution variance
      real*8  :: VEFF0=0.35d0,  REFF0=0.30d0      ! REFF0 not used

!@var NORMS0 if =1, Incident (TOA) Solar flux is normalized to equal S0
      integer :: NORMS0=1

!@var KORDER,KWTRAB controls WRITER-output (Mie-scattering info)
      integer :: KWTRAB=0, KORDER=0

C-----------------------------------------------------------------------
C      COMPOSITION & VERTICAL DISTRIBUTION FOR 5 SPECIFIED AEROSOL TYPES
C-----------------------------------------------------------------------
C TYPE
C    1   STRATOSPHERIC GLOBAL AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    2    TROPOSPHERIC  OCEAN AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    3    TROPOSPHERIC   LAND AEROSOL  A,B,C ARE GLOBAL AVERAGE VALUES
C    4    TROPOSPHERIC DESERT AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES
C    5    TROPOSPHERIC   HAZE AEROSOL  A,B,C ARE  LOCAL AVERAGE VALUES

C        1     2     3     4     5     6     7     8     9    10    11
C      ACID1 SSALT SLFT1 SLFT2 BSLT1 BSLT2 DUST1 DUST2 DUST3 CARB1 CARB2
      real*8, dimension(11,5) :: AGOLDH=reshape( (/
     1 .005,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, .020, .010, .010, .005,   .0, .010,   .0,   .0, .005,   .0,
     3   .0,   .0,   .0, .020, .005,   .0, .010, .010,   .0,   .0, .015,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, .020, .010,   .0,   .0,
     5   .0,   .0,   .0, .010,   .0,   .0,   .0,   .0,   .0,   .0, .005/
     *  ),(/11,5/) )
      real*8, dimension(11,5) :: BGOLDH=reshape( (/
     1 20.0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, 1.00, 4.00, 1.00, 4.00, 1.00, 4.00,   .0,   .0, 1.00,   .0,
     3   .0,   .0,   .0, 0.00, 2.00,   .0, 4.00, 2.00,   .0,   .0, 0.00,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 2.00, 0.00,   .0,   .0,
     5   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0, 0.00/
     *  ),(/11,5/) )
      real*8, dimension(11,5) :: CGOLDH=reshape( (/
     1 3.00,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,   .0,
     2   .0, 1.00, 3.00, 2.00, 3.00, 1.00, 2.00,   .0,   .0, 1.00,   .0,
     3   .0,   .0,   .0, 1.00, 3.00,   .0, 1.00, 1.00,   .0,   .0, 1.00,
     4   .0,   .0,   .0,   .0,   .0,   .0,   .0, 1.00, 1.00,   .0,   .0,
     5   .0,   .0,   .0, 1.00,   .0,   .0,   .0,   .0,   .0,   .0, 1.00/
     *  ),(/11,5/) )

!nu   real*8, dimension(11) :: PI0VIS=(/
!nu         1          2          3          4          5          6
!nu       ACID1      SSALT      SLFT1      SLFT2      BSLT1      BSLT2
!nu  1   1.00000,   1.00000,   1.00000,   1.00000,   0.98929,   0.95609,
!nu
!nu         7          8          9         10         11
!nu       DUST1      DUST2      DUST3      CARB1      CARB2
!nu  2   0.91995,   0.78495,   0.63594,   0.31482,   0.47513/)

      real*8, dimension(8) ::
C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                  SO4    SEA    ANT    OCX    BCI    BCB    DST   VOL
     *  REFDRY=(/0.200, 1.000, 0.300, 0.300, 0.100, 0.100, 1.000,1.000/)

     * ,REFWET=(/0.272, 1.808, 0.398, 0.318, 0.100, 0.100, 1.000,1.000/)

     * ,DRYM2G=(/4.667, 0.866, 4.448, 5.018, 9.000, 9.000, 1.000,1.000/)

CKoch   DRYM2G=(/5.000, 2.866, 8.000, 8.000, 9.000, 9.000, 1.000,1.000/)

!nu     RHTMAG=(/1.788, 3.310, 1.756, 1.163, 1.000, 1.000, 1.000,1.000/)
!nu alt RHTMAG=(/1.982, 3.042, 1.708, 1.033, 1.000, 1.000, 1.000,1.000/)
!nu  *  WETM2G=(/8.345, 2.866, 7.811, 5.836, 9.000, 9.000, 1.000,1.000/)
!nu alt WETM2G=(/9.250, 2.634, 7.598, 5.180, 9.000, 9.000, 1.000,1.000/)
     * ,Q55DRY=(/2.191, 2.499, 3.069, 3.010, 1.560, 1.560, 1.000,1.000/)

     * ,DENAER=(/1.760, 2.165, 1.725, 1.500, 1.300, 1.300, 2.000,2.000/)

C     TROP AEROSOL 1850 BACKGROUND, INDUSTRIAL & BIO-BURNING PARAMETERS
      real*8, dimension(13) :: AERMIX=(/
C      Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C      ---------------------------------  ------------------  --------
C       1    2    3    4    5    6    7    8    9   10   11   12   13
C      SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
     + 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9/)

      real*8, dimension(8) ::
C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                  SO4    SEA    ANT    OCX    BCI    BCB    DST   VOL
     *  FS8OPX=(/1.000, 1.000, 1.000, 1.000, 2.000, 2.000, 1.000, 1.00/)

     * ,FT8OPX=(/1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.300, 1.00/)

     * ,FRSULF=(/0.000, 0.000, 0.000, 0.330, 0.000, 0.000, 0.000, 1.00/)

     * ,PI0MAX=(/1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.00/)

!nu  * ,A8VEFF=(/ .200,  .200,  .200,  .200,  .200,  .200,  .200, .200/)

      real*8, dimension(8) ::
C                          MINERAL DUST PARAMETERS
C                         CLAY                  SILT
     *   REDUST=(/ 0.1, 0.2, 0.4, 0.8,   1.0, 2.0, 4.0, 8.0/)
!nu  *  ,VEDUST=(/ 0.2, 0.2, 0.2, 0.2,   0.2, 0.2, 0.2, 0.2/)
!nu  *  ,RODUST=(/ 2.5, 2.5, 2.5, 2.5,   2.6, 2.6, 2.6, 2.6/)
!nu  *  ,FSDUST=(/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/)
!nu  *  ,FTDUST=(/ 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0/)

C-----------------------------------------------------------------------
C     GHG 1980 Reference Concentrations and Vertical Profile Definitions
C-----------------------------------------------------------------------

!@var KTREND if > 0 table GHG concentrations (Trend G) are used for
!@+             yr/day KYEARG/KJDAYG; if KTREND=0, GHG are set to PPMVK0
      integer :: KTREND=1

!@var PPMV80  reference GHG concentrations (ppm)
      real*8, dimension(13) ::
C     GAS NUMBER    1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
     *   PPMV80=(/0d0, 337.90d0, 0d0,  21d4, 0d0,  .3012d0,   1.5470d0
     *     ,.1666d-03,.3003d-03, 0d0,   .978D-04,  .0010D-10,  .0420d0/)
C              CCL3F1    CCL2F2   N2     CFC-Y       CFC-Z         SO2
C     GAS NUMBER    8         9   10        11          12          13

!@var PPMVK0  user set  GHG concentrations (ppm), used if KTREND=0
      real*8, dimension(12) ::
C     GAS  NUMBER   1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
     *   PPMVK0=(/0d0, 337.90d0, 0d0, 21.d4, 0d0,  .3012d0,   1.5470d0
     *               ,.1666d-03,  .3003d-03, 0d0, .978D-04, 0.0010D-10/)
C                        CCL3F1      CCL2F2   N2     CFC-Y       CFC-Z
C     GAS  NUMBER             8           9   10        11          12

C     Makiko's GHG Trend Compilation  GHG.1850-2050.Dec1999 in GTREND
C     ---------------------------------------------------------------
!@var nghg nr. of well-mixed GHgases: CO2 N2O CH4 CFC-11 CFC-12 others
!@var nyrsghg max.number of years of prescr. greenhouse gas history
      integer, parameter :: nghg=6, nyrsghg=2050-1850+1

!@var ghgyr1,ghgyr2 first and last year of GHG history
      integer ghgyr1,ghgyr2
!@var ghgam,xref,xnow     GHG-mixing ratios in ppm,ppm,ppm,ppb,ppb,ppb
      real*8 GHGAM(nghg,nyrsghg),XREF(nghg+1),XNOW(nghg+1)
      common/ghgcom/ghgyr1,ghgyr2,ghgam,xref,xnow

C     GTREND:  1980.,  337.9,  .3012,  1.547,  .1666,  .3003,  .0978,
C     ---------------------------------------------------------------

!@var KGGVDF,KPGRAD,KLATZ0 control parameters for vertical GHG profiles
!@+   -----------------------------------------------------------------
!@+   Minschwaner et al JGR (1998) CH4, N2O, CFC-12 Vertical profiles
!@+   IF(KGGVDF.GT.0) Then:
!@+      Gas decreases are linear with pressure, from unity at ground to
!@+      the fractional value PPMVDF(NGAS) at the top of the atmosphere.
!@+   Exponential decrease by EXP(-(Z-Z0)/H) is superimposed on this.
!@+   IF(KLATZ0.GT.0) Then: Z0 depends on latitude, KGGVDF not used
!@+   KPGRAD>0: Pole-to-Pole lat. gradient (PPGRAD) is also superimposed
!@+   ------------------------------------------------------------------
!@var Z0,ZH   scale heights used for vertical profile (km)
!@var PPMVDF  frac. value at top of atmosphere (used if KGGVDF > 0)
!@var PPGRAD  Pole-to-Pole latitud.gradient for GHG (used if KPGRAD > 0)
      integer :: KGGVDF=0, KPGRAD=1, KLATZ0=1

      real*8, dimension(12) ::
C     NUMBER   1    2    3    4  5    6    7    8     9   10   11  12
C             H2O  CO2  O3   O2 NO2  N2O  CH4 CFC11 CFC12 N2 CF-Y  CF-Z
     *   Z0=(/0.0, 0.0,0.0, 0.0,0.0, 16., 16., 16., 16., 0.0, 16., 16./)
     *  ,ZH=(/8.0, 8.0,8.0, 8.0,8.0, 30., 50., 30., 30., 0.0, 30., 30./)

C     GAS NUMBER    1     2    3    4    5         6         7
C                 H2O   CO2   O3   O2  NO2       N2O       CH4
     *  ,PPMVDF=(/1.0,  1.0, 1.0, 1.0, 1.0,  0.88888,  0.88888,
     *              0.88888,  0.88888, 1.0,  0.88888,  0.88888/)
C                    CCL3F1    CCL2F2   N2     CFC-Y     CFC-Z
C     GAS NUMBER          8         9   10        11        12

C     GAS  NUMBER   1     2    3    4    5         6         7
C                 H2O   CO2   O3   O2  NO2       N2O       CH4
     *  ,PPGRAD=(/0.0,  0.0, 0.0, 0.0, 0.0,   0.0100,   0.0900,
     *               0.0600,   0.0600, 0.0,   0.0600,   0.0600/)
C                    CCL3F1    CCL2F2   N2     CFC-Y     CFC-Z
C     GAS  NUMBER         8         9   10        11        12

C---------------------
C     Optional Tracers    used via setbak/getbak
C---------------------
      integer, dimension(8) :: ITR=(/0,0,0,0, 0,0,0,0/)
      integer :: NTRACE=0

      real*8, dimension(8) ::
C                TRACER AEROSOL COMPOSITIONAL/TYPE PARAMETERS
     *  TRRDRY=(/ .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0/)
!nu  * ,TRVEFF=(/ .2d0, .2d0, .2d0, .2d0, .2d0, .2d0, .2d0, .2d0/)
!nu  * ,TRADEN=(/ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
!loc * ,FSTOPX=(/ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
!loc * ,FTTOPX=(/ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)

      SAVE

      CONTAINS

      SUBROUTINE RCOMP1(NRFUN)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     Solar,GHG Trend, VolcAer Size Selection Parameters:    Defaults
C                                           Process       KYEARX  KJDAYX
c                                         SolarCon, UV       0       0
c                                         GH Gas Trend       0       0
c                                                         REFF0= 0.3
c                                                         VEFF0= 0.35
C     ------------------------------------------------------------------

c     NRFUN is now set as an argument from calling routine so that unit
c     numbers can be set automatically
      INTEGER, DIMENSION(14) :: NRFUN
C          radfile1   2   3   4   5   6   7   8   9   A   B   C   D   E
!?    DATA NRFN0/71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84/

      INTEGER, SAVE :: IFIRST=1 ! ,NRFN0
      CHARACTER*80 EPSTAG,TITLE

      REAL*4 OZONLJ(44,46),R72X46(72,46),VTAUR4(1800,24)
      REAL*8 :: EJMLAT(47),E20LAT(20)
      INTEGER :: I,J,K,L,M,N,N1,N2,NRFU,KK,NN,IYEAR,IMONTH,JJDAYS,JYEARS
     *     ,JJDAYG,JYEARG
      REAL*8 :: WAVNA,WAVNB,PFWI,TKOFPF,SUM,EPK,EPL,DEP,SFNORM,D,O,Q,S
     *     ,OCM,WCM

!?    IF(LASTVC.GT.0) NRFUN=NRFN0
      IF(IFIRST.LT.1) GO TO 9999

C     ------------------------------------------------------------------
C     Input data are read as specified in the first CALL RCOMP1 (NRFUN).
C     Subsequent calls to RCOMP1 can be used to re-initialize parameters
C     in SETXXX subroutines to different values, but no new data is read
C     ------------------------------------------------------------------

C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips input process)
C
C     MADO3M   =  1   Reads  Decadal Ozone files and Ozone trend file
C     MADAER   =  2   Reads  Aerosol 50y tropospheric climatology RFILE5
C     MADDST   =  3   Reads  Dust-windblown mineral climatology   RFILE6
C     MADVOL   =  4   Reads  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS   =  5   Reads  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV   =  6   Reads  Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in RADPAR
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1   Reads  V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------


C              Initialize variables that might not otherwise get defined
C              ---------------------------------------------------------

      DO 110 L=1,LX
      TAUWC(L)=0.D0
      TAUIC(L)=0.D0
      SIZEWC(L)=0.D0
      SIZEIC(L)=0.D0
      CLDEPS(L)=0.D0
      FPXCO2(L)=1.D0
      FPXOZO(L)=1.D0
      TLB(L)=250.D0
      TLT(L)=250.D0
      TLM(L)=250.D0
      SHL(L)=0.D0
      RHL(L)=0.D0
      DO 101 I=1,6
      SRAEXT(L,I)=0.D0
      SRASCT(L,I)=0.D0
      SRAGCB(L,I)=0.D0
      SRBEXT(L,I)=0.D0
      SRBSCT(L,I)=0.D0
      SRBGCB(L,I)=0.D0
      SRDEXT(L,I)=0.D0
      SRDSCT(L,I)=0.D0
      SRDGCB(L,I)=0.D0
      SRVEXT(L,I)=0.D0
      SRVSCT(L,I)=0.D0
      SRVGCB(L,I)=0.D0
      DBLEXT(L,I)=0.D0
      DBLSCT(L,I)=0.D0
      DBLGCB(L,I)=0.D0
      DBLPI0(L,I)=0.D0
      SRCEXT(L,I)=0.D0
      SRCSCT(L,I)=0.D0
      SRCGCB(L,I)=0.D0
      SRCPI0(L,I)=0.D0
  101 CONTINUE
      DO 102 I=1,33
      TRAALK(L,I)=0.D0
      TRBALK(L,I)=0.D0
      TRDALK(L,I)=0.D0
      TRVALK(L,I)=0.D0
      TRCALK(L,I)=0.D0
      TRGXLK(L,I)=0.D0
  102 CONTINUE
      DO 103 I=1,13
      U0GAS(L,I)=0.D0
      ULGAS(L,I)=0.D0
  103 CONTINUE
      DO 104 I=1,8
      TRACER(L,I)=0.D0
  104 CONTINUE
  110 CONTINUE

      IF(LASTVC.GT.0) CALL SETATM
      IF(NLP.GT.LX)   call stop_model('rcomp1: increase LX',255)

C**** Use (global mean) pressures to get standard mid-latitude summer
C**** values for height, density, temperature, ozone, water vapor
      DO 120 L=1,NLP
      PLB0(L)=PLB(L)
      CALL PHATMO(PLB0(L),HLB0(L),D,TLB(L),O,Q,S,OCM,WCM,1,2)
  120 CONTINUE
      DO 121 L=1,NL
      TLT(L)=TLB(L+1)
      TLM(L)=0.5D0*(TLB(L)+TLT(L))
  121 CONTINUE

!sl   De-activate surface layer computations
!sl   DO I=1,33
!sl      TAUSL(I)=0.0
!sl      FTAUSL(I)=0.0
!sl   END DO

C-----------------------------------------------------------------------
CR(1) Reads GTAU Asymmetry Parameter Conversion Table used within SGPGXG
C
C       (SGPGXG does Multiple Scattering Parameterization used in SOLAR)
C       ----------------------------------------------------------------

      NRFU=NRFUN(1)
      READ (NRFU) GTAU,TGDATA
      CALL SETGTS


C-----------------------------------------------------------------------
CR(2)    Reads in Merged k-Distribution Tau Tables for Thermal Radiation
C        CFCs, H2O Continuum Tau Table, Merged k-Distr Planck Flux Table
C
C        (Reads: TAUTBL,TAUWV0,PLANCK,XKCFC,H2OCN8,H2OCF8
c                DUCH4,SDUCH4,DUN2O,SDUN2O,ULOX,DUX      used in TAUGAS)
C       ----------------------------------------------------------------

      NRFU=NRFUN(2)
      READ(NRFU) TAUTBL,TAUWV0,PLANCK,XKCFC,H2OCN8,H2OCF8, XUN2O,XUN2O0,
     *  XUCH4,XUCH40, XTRUP,XTU0, XTRDN,XTD0, CXUCO2,CXUO3,ULOX,DUX


C        Define Window Flux to Brightness Temperature Conversion Factors
C        ---------------------------------------------------------------

      WAVNA=850.D0
      WAVNB=900.D0
      DO I=1,630
      PFWI=0.001D0*I
      IF(I.GT.100) PFWI=(0.1D0+0.01D0*(I-100))
      IF(I.GT.190) PFWI=(1.0D0+0.10D0*(I-190))
      TKPFW(I)=TKOFPF(WAVNA,WAVNB,PFWI)
      END DO
      DO I=1,900
      PFWI=I
      TKPFT(I)=TKOFPF(0.D0,10000.D0,PFWI)
      END DO

C                            PLANCK Table interpolation limit parameters
C                            -------------------------------------------
C     ITPFT0=123
C     ITNEXT=250
C-----------------------------------------------------------------------
CR(3)        Read Mie Scattering Parameters [Qext, Qscat, AsymParameter]
C            (1) Tropospheric Aerosols [11 Background, 8 Trop8 Aerosols]
C            (2) Clouds [5 Water, 5 non-spherical Ice, 5 Mie Ice Clouds]
C            (3) Desert Dust Aerosols  [25 particle sizes - to select 8]
C            (4) Volcanic Aerosols [20 particle sizes, 5 size variances]
C            (5) Sulfate  Aerosols [22 particle sizes, 0.1 - 10. micron]
C            (6) Soot   Aerosols [25 particle sizes, 0.001 - 5.0 micron]
C            -----------------------------------------------------------

      NRFU=NRFUN(3)

C                               GCM 11 background aerosol Mie parameters
C                               ----------------------------------------
      DO 301 N=1,11
      READ (NRFU,3000) TITLE
 3000 FORMAT(A80)
      READ (NRFU,3001) (SRAQEX(K,N),K=1,6)
 3001 FORMAT( 18X,6(F7.5,1X))
      READ (NRFU,3001) (SRAQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRAQCB(K,N),K=1,6)
  301 CONTINUE
      READ (NRFU,3002) (Q55A11(N),N=1,11)
 3002 FORMAT( 18X,6(F7.5,1X)/18X,6(F7.5,1X))
      READ (NRFU,3003) (REFA11(N),N=1,11)
 3003 FORMAT( 18X,6(F7.3,1X)/18X,5(F7.3,1X))
      READ (NRFU,3003) (VEFA11(N),N=1,11)
      DO 302 N=1,11
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRAQEX(K,N),K=1,33)
 3004 FORMAT( 14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
 3005 FORMAT(/14X,7(F7.5,1X),4(/14X,7(F7.5,1X)))
      READ (NRFU,3005) (TRAQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRAQCB(K,N),K=1,33)
  302 CONTINUE

C                       GCM 9 (of 10) climatology aerosol Mie parameters
C                       ------------------------------------------------
      DO 303 N=1,10
      IF(N.EQ.6) GO TO 303
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRBQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRBQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRBQCB(K,N),K=1,6)
  303 CONTINUE
      READ (NRFU,3002) (Q55B10(N),N=1,5),(Q55B10(N),N=7,10)
      READ (NRFU,3003) (REFB10(N),N=1,5),(REFB10(N),N=7,10)
      READ (NRFU,3003) (VEFB10(N),N=1,5),(VEFB10(N),N=7,10)
      DO 304 N=1,10
      IF(N.EQ.6) GO TO 304
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRBQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRBQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRBQCB(K,N),K=1,33)
  304 CONTINUE


C                               Cloud Water, Ice-non, Ice-Mie parameters
C                               ----------------------------------------
      DO 305 N=1,15
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRCQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRCQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRCQCB(K,N),K=1,6)
  305 CONTINUE
      READ (NRFU,3006) (Q55C15(N),N=1,15)
 3006 FORMAT( 18X,6(F7.5,1X)/18X,6(F7.5,1X)/18X,6(F7.5,1X))
      READ (NRFU,3007) (REFC15(N),N=1,15)
 3007 FORMAT( 18X,6(F7.3,1X)/18X,6(F7.3,1X)/18X,6(F7.3,1X))
      READ (NRFU,3007) (VEFC15(N),N=1,15)
      DO 306 N=1,15
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRCQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRCQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRCQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRCQAL(K,N),K=1,33)
  306 CONTINUE

C                               Desert Dust 25 sizes, Mie parameter data
C                               ----------------------------------------
      DO 307 N=1,25
      READ (NRFU,3001) (SRDQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRDQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRDQCB(K,N),K=1,6)
  307 CONTINUE
      READ (NRFU,3008) (Q55D25(N),N=1,25)
 3008 FORMAT( 18X,5(F7.5,1X),4(/18X,5(F7.5,1X)))
      READ (NRFU,3009) (REFD25(N),N=1,25)
 3009 FORMAT( 18X,12(F3.1,1X)/18X,12(F3.1,1X)/18X,F3.0)
      READ (NRFU,3010) (VEFD25(N),N=1,25)
 3010 FORMAT( 18X,12(F3.1,1X)/18X,12(F3.1,1X)/18X,F3.1)
      DO 308 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRDQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRDQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRDQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRDQAL(K,N),K=1,33)
  308 CONTINUE

      TRDQAB(:,:)=TRDQEX(:,:)-TRDQSC(:,:)  !  used in writer only

C                               Volcanic aerosol Mie size, variance data
C                               ----------------------------------------
      DO 313 M=1,5
      IF(M.EQ.4) GO TO 313
      DO 311 N=1,20
      READ (NRFU,3001) (SRVQEX(K,N,M),K=1,6)
      READ (NRFU,3001) (SRVQSC(K,N,M),K=1,6)
      READ (NRFU,3001) (SRVQCB(K,N,M),K=1,6)
  311 CONTINUE
      READ (NRFU,3011) (Q55V20(N,M),N=1,20)
 3011 FORMAT( 18X,5(F7.5,1X),3(/18X,5(F7.5,1X)))
      READ (NRFU,3012) (REFV20(N,M),N=1,20)
 3012 FORMAT( 18X,12(F3.1,1X)/18X,8(F3.1,1X))
      READ (NRFU,3012) (VEFV20(N,M),N=1,20)
      DO 312 N=1,20
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRVQEX(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQSC(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQCB(K,N,M),K=1,33)
      READ (NRFU,3005) (TRVQAL(K,N,M),K=1,33)
  312 CONTINUE
  313 CONTINUE
      DO 316 N=1,20
      DO 314 K=1,6
      SRVQEX(K,N,4)=(SRVQEX(K,N,3)+SRVQEX(K,N,5))/2.D0
      SRVQSC(K,N,4)=(SRVQSC(K,N,3)+SRVQSC(K,N,5))/2.D0
      SRVQCB(K,N,4)=(SRVQCB(K,N,3)+SRVQCB(K,N,5))/2.D0
  314 CONTINUE
      Q55V20(N,4)=(Q55V20(N,3)+Q55V20(N,5))/2.D0
      REFV20(N,4)=(REFV20(N,3)+REFV20(N,5))/2.D0
      VEFV20(N,4)=(VEFV20(N,3)+VEFV20(N,5))/2.D0
      DO 315 K=1,33
      TRVQEX(K,N,4)=(TRVQEX(K,N,3)+TRVQEX(K,N,5))/2.D0
      TRVQSC(K,N,4)=(TRVQSC(K,N,3)+TRVQSC(K,N,5))/2.D0
      TRVQCB(K,N,4)=(TRVQCB(K,N,3)+TRVQCB(K,N,5))/2.D0
      TRVQAL(K,N,4)=(TRVQAL(K,N,3)+TRVQAL(K,N,5))/2.D0
  315 CONTINUE
  316 CONTINUE

C                            Sulfate aerosol, Mie parameter 22-size data
C                            -------------------------------------------
      DO 321 N=1,22
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRUQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRUQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRUQCB(K,N),K=1,6)
  321 CONTINUE
      READ (NRFU,3008) (Q55U22(N),N=1,22)
      READ (NRFU,3013) (REFU22(N),N=1,22)
 3013 FORMAT( 18X,5(F7.3,1X),4(/18X,5(F7.3,1X)))
      READ (NRFU,3013) (VEFU22(N),N=1,22)
      DO 322 N=1,22
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRUQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRUQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRUQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRUQAL(K,N),K=1,33)
  322 CONTINUE

C                               Soot aerosol, Mie parameter 25-size data
C                               ----------------------------------------
      DO 323 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRSQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRSQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRSQCB(K,N),K=1,6)
  323 CONTINUE
      READ (NRFU,3008) (Q55S25(N),N=1,25)
      READ (NRFU,3013) (REFS25(N),N=1,25)
      READ (NRFU,3013) (VEFS25(N),N=1,25)
      DO 324 N=1,25
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRSQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRSQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRSQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRSQAL(K,N),K=1,33)
  324 CONTINUE

C                            Seasalt aerosol, Mie parameter 22-size data
C                            Nitrate aerosol, Mie parameter 22-size data
C                            (Water) aerosol, Mie parameter 22-size data
C                            Organic aerosol, Mie parameter 22-size data
C                            -------------------------------------------
      N1=23
      DO 326 KK=1,4
      N2=N1+21
      DO 325 N=N1,N2
      READ (NRFU,3000) TITLE
      READ (NRFU,3001) (SRUQEX(K,N),K=1,6)
      READ (NRFU,3001) (SRUQSC(K,N),K=1,6)
      READ (NRFU,3001) (SRUQCB(K,N),K=1,6)
  325 CONTINUE
      READ (NRFU,3008) (Q55U22(N),N=N1,N2)
      READ (NRFU,3013) (REFU22(N),N=N1,N2)
      READ (NRFU,3013) (VEFU22(N),N=N1,N2)
      N1=N2+1
  326 CONTINUE
      N1=23
      DO 328 KK=1,4
      N2=N1+21
      DO 327 N=N1,N2
      READ (NRFU,3000) TITLE
      READ (NRFU,3004) (TRUQEX(K,N),K=1,33)
      READ (NRFU,3005) (TRUQSC(K,N),K=1,33)
      READ (NRFU,3005) (TRUQCB(K,N),K=1,33)
      READ (NRFU,3005) (TRUQAL(K,N),K=1,33)
  327 CONTINUE
      N1=N2+1
  328 CONTINUE


C-----------------------------------------------------------------------
CR(6) DUST:   Monthly-Mean Desert Dust (Clay,Silt) 8-Size Optical Depths
C                  Map: IJ=72x46,  Lay: L=1-9, Siz: S=1-8, Month: M=1-12
C                  -----------------------------------------------------

      IF(MADDST.LT.1) GO TO 699
      NRFU=NRFUN(6)
      READ (NRFU) TDUST

  699 CONTINUE

C-----------------------------------------------------------------------
CR(7)        Read Makiko's Stratospheric binary data made in April, 2002
C                               (1800 months (1850-1999) x 24 latitudes)
C                               ----------------------------------------

      IF(MADVOL.LT.1) GO TO 799
      NRFU=NRFUN(7)
      READ (NRFU) TITLE
      IF(TITLE(1:13).eq.'Optical Depth')
     &     call stop_model('rcomp1: use new RADN7',255)
      REWIND (NRFU)
      DO K=1,5
        READ (NRFU) TITLE,VTAUR4
        DO J=1,24
        DO I=1,1800
          V4TAUR(I,J,K)=VTAUR4(I,J)
        END DO
        END DO
      END DO

  799 CONTINUE

C-----------------------------------------------------------------------
CR(8)  ISCCP Derived Cloud Variance (EPSILON) Cloud Optical Depth Factor
C      Low, Mid, High  Cloud Optical Depths are Reduced by (1 - EPSILON)
C
C               INPUT DATA FILE:  UNIT = INFILE
C                                 TAG  = EPSTAG  (CHARACTER*80)
C                                 DATA = EPLMHC  (72,46,12,4) REAL*4
C
C      Data are 72X46 Monthly Mean Low, Mid, High, Column EPSILON Values
C      Cloud Heterogeneity selections used in UPDEPS, GETEPS (in SETCLD)
C
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C               --------------------------------------------------------

      IF(MADEPS.LT.1) GO TO 899
      NRFU=NRFUN(8)
      READ (NRFU) EPSTAG,EPLMHC

      DO 811 N=1,4
      DO 810 M=1,12
      DO 809 I=1,72
      L=0
  801 CONTINUE
      L=L+1
      IF(EPLMHC(I,L,M,N).LT.0.0) GO TO 801
      IF(L.GT.1) THEN
      DO 802 J=1,L
      EPLMHC(I,J,M,N)=EPLMHC(I,L,M,N)
  802 CONTINUE
      ENDIF
      L=47
  803 CONTINUE
      L=L-1
      IF(EPLMHC(I,L,M,N).LT.0.0) GO TO 803
      IF(L.LT.46) THEN
      DO 804 J=L,46
      EPLMHC(I,J,M,N)=EPLMHC(I,L,M,N)
  804 CONTINUE
      ENDIF
      J=0
  805 CONTINUE
      J=J+1
      IF(J.GT.46) GO TO 808
      IF(EPLMHC(I,J,M,N).GE.0.0) GO TO 805
      K=J-1
      EPK=EPLMHC(I,K,M,N)
  806 CONTINUE
      J=J+1
      IF(J.GT.46) GO TO 808
      IF(EPLMHC(I,J,M,N).LT.0.0) GO TO 806
      L=J
      EPL=EPLMHC(I,L,M,N)
      NN=L-K-1
      DEP=(EPL-EPK)/(L-K)
      DO 807 L=1,NN
      EPLMHC(I,K+L,M,N)=EPK+L*DEP
  807 CONTINUE
      GO TO 805
  808 CONTINUE
  809 CONTINUE
  810 CONTINUE
  811 CONTINUE

  899 CONTINUE


C-----------------------------------------------------------------------
CR(E)
C             KCLDEM  Selects: Top-Cloud (Thermal) Scattering Correction
C                     KCLDEM =  0  Utilizes Non-scattering approximation
C                     KCLDEM =  1  Modifies emission and transmission by
C                                  top cloud (over-rides old correction)
C             ----------------------------------------------------------

      NRFU=NRFUN(14)
      READ (NRFU) RIJTPG,FDXTPG,FEMTPG


C-----------------------------------------------------------------------
CR(9)         Read Judith Lean's Solar UV and Solar Constant Variability
C                                      Monthly-Mean Solar UV (1882-1998)
C                                      ---------------------------------

      IF(KSOLAR.LT.0) GO TO 949
      IF(MADLUV.LT.1) GO TO 909
      NRFU=NRFUN(9)

      READ(NRFU,'(a80)') TITLE
      if(ksolar.ge.2 .and. TITLE(1:3).ne.'ANN')
     &     call stop_model('rcomp1: change RADN9 to ann.file',255)
      if(ksolar.lt.2 .and. TITLE(1:3).eq.'ANN')
     &     call stop_model('rcomp1: change RADN9 to monthly file',255)
      READ(NRFU,'(5F14.2)') (WSLEAN(I),I=1,190)
      READ(NRFU,'(a80)') TITLE
      READ(NRFU,'(5E14.3)') (DSLEAN(I),I=1,190)
      DO I=1,190
        WSLEAN(I)=WSLEAN(I)/1000.D0
        DSLEAN(I)=DSLEAN(I)/1000.D0
        W1LEAN(I)=WSLEAN(I)-0.5D0*DSLEAN(I)
      END DO
      READ(NRFU,'(a80)') TITLE
      READ(NRFU,'(a80)') TITLE
      IF(KSOLAR.LT.2) THEN
C****   Read in monthly-mean data
        DO I=1,Ms0X
          READ(NRFU,'(2I6,3F17.6)') IYEAR,IMONTH,TSI1(I),TSI2(I)
          READ(NRFU,'(5E14.6)')     (FSLEAN(K),K=1,190)
          SUM=0.D0
          DO K=1,190
            SUM=SUM+FSLEAN(K)*DSLEAN(K)
          END DO
          SFNORM=TSI1(I)/SUM
          DO K=1,190
            UVLEAN(I,K)=FSLEAN(K)*SFNORM
          END DO
        END DO
      ELSE
C****   Read in annual-mean data
        DO I=1,Ms0X
          READ(NRFU,'(F12.1,2F15.4)',end=908) yr2S0,TSI1(I),TSI2(I)
          if(I.eq.1) yr1S0 = yr2S0
          READ(NRFU,'(5E14.6)')               (UVLEAN(I,K),K=1,190)
        END DO
  908   write(6,*) 'read S0-history: ',yr1S0,' - ',yr2S0
      END IF
      GO TO 949

  909 CONTINUE
      DO I=1,190
        WSLEAN(I)=WSLEAN(I)/1000.D0
        DSLEAN(I)=DSLEAN(I)/1000.D0
        W1LEAN(I)=WSLEAN(I)-0.5D0*DSLEAN(I)
      END DO

  949 CONTINUE


C-----------------------------------------------------------------------
CR(C)     Read:    Elaine Mathews 10 Fractional Vegetation Distributions
C         10 global maps (72x46) depict fractional vegetation/soil types
C         Map-1 (bright sand) + Map-10 (black dirt) define desert albedo
C         (sum of Maps 1-10 over land-area (ILON,JLAT) grid boxes = 1.0)
C
C         Map-11 refers to plankton concentrations over ocean areas that
C         are yet to be implemented.
C         --------------------------------------------------------------





C-----------------------------------------------------------------------
CR(D)      Read:   1   FOCEAN   72x46 ocean fraction   (FOCEAN = 0 or 1)
C                  2   FLAKE    72x46 lake  fraction
C                  3   FGRND    72x46 lake  fraction
C                  4   FGICE    72x46 glacial ice fraction
C                               (FLAKE + FGRND + FGICE + FOCEAN = 1.000)
C


C                  5   ZATMO    72x46 topography (ocean = 0.0)
C                  6   HOCEAN   72x46 ocean depth
C                  7   HLAKE    72x46 lake  depth
C                  8   HGICE    72x46 glice depth
C                  9   ZSOLID   72x46 topography of solid ground surface
C                  -----------------------------------------------------
C
C     FOLGIZ is for off-line use only, and is not used in GCM radiation.
C     GCM supplies dynamically changing POCEAN,POICE,PEARTH,PLICE values
C     ------------------------------------------------------------------



  999 CONTINUE

      IFIRST=0
 9999 CONTINUE


C     ---------------------------------------------------------------
C     LASTVC     Initialize:  Default Atmospheric Layering, Structure
C                (for Off-Line use)  as Specified by LASTVC Parameter
C                If LASTVC < 0, GCM defines all Radiation Model Input
C     otherwise:
C                Each LASTVC digit(6) specifies a model configuration
C          e.g.:    LASTVC= 123456
C                L=0,1,..9  Layers NL=  Any,GCM12,GCM23,Pset,Hset,etc
C                A=0,1,..6  Atmosphere  Any,Trop,MLS,MLW,SAS,SAW,Std
C                S=0,1,..9  Surf Types  POCEAN=1,PEARTH=1,POICE=1,etc
C                T=0,1,..9  Tracer Aer  Tau=0,  Tau=0.1 Aer Comp(1-9)
C                V=0,1,..9  Vegetation  Sand,Tundra,Grass,Shrubs, etc
C                C=0,1,..9  Cloud,R=10  Clim Cloud Tau in Layer(1,-9)
C                ----------------------------------------------------

      IF(LASTVC.GE.0) CALL SETATM

C             -------------------------------------------------------
C             Set Solar Constant for Default Reference Time: Jan 1950
C             Default used for KSOLAR(=1) is that specified in RADPAR
C             -------------------------------------------------------

      JJDAYS=1
      JYEARS=1950
      IF(KJDAYS.GT.0)             JJDAYS=KJDAYS
      IF(KYEARS.GT.0)             JYEARS=KYEARS
C----------------------------------------------
                      CALL SETSOL(JYEARS,JJDAYS)
C----------------------------------------------


C             -------------------------------------------------------
C             Set Default Greenhouse Gas Reference Year to:  Mid 1980
C             Default used for KTREND(=1) is that specified in RADPAR
C             -------------------------------------------------------

      JJDAYG=184
      JYEARG=1980
C----------------------------------------------
                      CALL SETGHG(JYEARG,JJDAYG)
C----------------------------------------------
      IF(KJDAYG.GT.0)             JJDAYG=KJDAYG
      IF(KYEARG.GT.0)             JYEARG=KYEARG
C----------------------------------------------
                      CALL UPDGHG(JYEARG,JJDAYG)
C----------------------------------------------

C--------------------------------
                      CALL SETGAS
C
                                     CALL SETBAK
      IF(MADAER.GT.0.or.NTRACE.gt.0) CALL SETAER
      IF(MADDST.GT.0) CALL SETDST
C--------------------------------


C               -----------------------------------------------------
C               Set Volcanic Aerosol Effective Variance Default Value
C                   Particle Size(REFF0=0.3) when not known from data
C                   (VEFF0=0.35 is value based on thermal ISAMS data)
C                   -------------------------------------------------

C----------------------------------------------
      IF(MADVOL.GT.0) CALL SETVOL
C----------------------------------------------

C--------------------------------
                      CALL SETCLD

C--------------------------------

      CALL SOLAR0

      RETURN
      END SUBROUTINE RCOMP1

      SUBROUTINE RCOMPT
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C     Time Trend Selection Parameters and Options:
C     -------------------------------------------
C
C             The Nominal Default Values are KYEARX = 0, and KJDAYX = 0,
C             in which case RADPAR supplied Time JYEAR and JDAY are used
C
C             When Non-Zero Values are specified for  KYEARX and KJDAYX,
C             the JYEAR,JDAY Time Dependence of the Specified Process is
C             over-ridden by the Non-Zero KYEARX and KJDAYX Value.
C             ----------------------------------------------------------
C                                  Process       KYEARX  KJDAYX
c             KYEARS,KJDAYS        SolarCon, UV       0       0
c             KYEARG,KJDAYG        GH Gas Trend       0       0
c             KYEARO,KJDAYO        Ozone Distr        0       0
c             KYEARA,KJDAYA        AerClimtolgy       0       0
c             KYEARD,KJDAYD        Desert Dust        0       0
c             KYEARV,KJDAYV        Volcanic Aer       0       0
c             KYEARE,KJDAYE        Epsilon Clds       0       0
c             KYEARR,KJDAYR        Refl Surface       0       0

C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips input process)
C
C     MADAER  =  2  Updates  Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  3  Updates  Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  4  Updates  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  5  Updates  Epsilon cloud heterogeneity data     RFILE8
C     MADLUV  =  6  Updates  Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in RADPAR
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1          V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------
      INTEGER JJDAYS,JYEARS,JJDAYG,JYEARG,JJDAYO,JYEARO,JJDAYA,JYEARA
     *     ,JJDAYD,JYEARD,JJDAYV,JYEARV,JJDAYE,JYEARE,JJDAYR,JYEARR

C                      -------------------------------------------------
C                      Set Seasonal and Time (JDAY) Dependent Quantities
C                      -------------------------------------------------

      JJDAYS=JDAY
      JYEARS=JYEAR
      IF(KJDAYS.GT.0)             JJDAYS=KJDAYS
      IF(KYEARS.GT.0)             JYEARS=KYEARS
C----------------------------------------------
      IF(MADLUV.GT.0) CALL UPDSOL(JYEARS,JJDAYS)
C----------------------------------------------

      JJDAYG=JDAY
      JYEARG=JYEAR
      IF(KJDAYG.GT.0)             JJDAYG=KJDAYG
      IF(KYEARG.GT.0)             JYEARG=KYEARG
C----------------------------------------------
      IF(MADGHG.GT.0) CALL UPDGHG(JYEARG,JJDAYG)
C----------------------------------------------


      JJDAYO=JDAY
      JYEARO=JYEAR
      IF(KJDAYO.NE.0)             JJDAYO=KJDAYO
      IF(KYEARO.NE.0)             JYEARO=KYEARO
C----------------------------------------------
                      CALL UPDO3D(JYEARO,JJDAYO)
C----------------------------------------------

      JJDAYA=JDAY
      JYEARA=JYEAR
      IF(KJDAYA.GT.0)             JJDAYA=KJDAYA
      IF(KYEARA.GT.0)             JYEARA=KYEARA
C----------------------------------------------
      IF(MADAER.NE.0) CALL UPDAER(JYEARA,JJDAYA)
C----------------------------------------------

      JJDAYD=JDAY
      JYEARD=JYEAR
      IF(KJDAYD.GT.0)             JJDAYD=KJDAYD
      IF(KYEARD.GT.0)             JYEARD=KYEARD
C----------------------------------------------
      IF(MADDST.GT.0) CALL UPDDST(JYEARD,JJDAYD)
C----------------------------------------------

      JJDAYV=JDAY
      JYEARV=JYEAR
      IF(KJDAYV.GT.0)             JJDAYV=KJDAYV
      IF(KYEARV.GT.0)             JYEARV=KYEARV
C----------------------------------------------
      IF(MADVOL.GT.0) CALL UPDVOL(JYEARV,JJDAYV)
C----------------------------------------------

      JJDAYE=JDAY
      JYEARE=JYEAR
      IF(KJDAYE.GT.0)             JJDAYE=KJDAYE
      IF(KYEARE.GT.0)             JYEARE=KYEARE
C----------------------------------------------
      IF(MADEPS.GT.0) CALL UPDEPS(JYEARE,JJDAYE)
C----------------------------------------------

      JJDAYR=JDAY
      JYEARR=JYEAR
      IF(KJDAYR.GT.0)             JJDAYR=KJDAYR
      IF(KYEARR.GT.0)             JYEARR=KYEARR
C----------------------------------------------
                      CALL UPDSUR(JYEARR,JJDAYR)
C----------------------------------------------

      RETURN
      END SUBROUTINE RCOMPT

      SUBROUTINE RCOMPX
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     MADVEL  Model Add-on Data of Extended Climatology Enable Parameter
C             Each MADVEL digit is ON/OFF switch for corresponding input
C             e.g. MADVEL=123456   (zero digit skips process)
C
C     MADO3M  =  1           Makiko's 1951-1997 Ozone climatology RFILEA
C     MADAER  =  2  Updates  Aerosol 50y tropospheric climatology RFILE5
C     MADDST  =  3  Updates  Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  4  Updates  Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  5           Epsilon cloud heterogeneity data     RFILE8
C     MADLUV  =  6           Lean's SolarUV 1882-1998 variability RFILE9
C
C                 Related Model Add-on Data Parameters set in RADPAR
C
C     MADGHG   =  1  Default Enables UPDGHG update. (MADGHG=0),no update
C     MADSUR   =  1          V72X46N.1.cor Vegetation type data   RFILEC
C                            Z72X46N Ocean fraction, topography   RFILED
C     ------------------------------------------------------------------
C
C      -----------------------------------------------------------------
C      Get Surface, Atmosphere, Sun Angle, Radiative Forcing, etc. Input
C      to compute Solar/Thermal Radiation for given (JLAT,ILON) Grid-box
C
C      The Radiation Model utilizes Data with 72x46 (lon,lat) resolution
C                 for GCM resolution other than 72x46, set JLAT and ILON
C                 to appropriately Sample  (rather than interpolate) the
C                 72x46 aerosol, ozone, cloud heterogeneity data sets
C
C      The Radiation Model can accommodate arbitrary vertical resolution
C      -----------------------------------------------------------------


C--------------------------------
!!!                   CALL GETO3D(ILON,JLAT)
      CALL REPART(O3JDAY(1,ILON,JLAT),PLBO3,NLO3+1,U0GAS(1,3),PLB0,NL+1)
      O3_OUT(:)=U0GAS(:,3) ! to save 3D ozone field in SUBR. RADIA
      if(use_tracer_ozone .eq. 1) U0GAS(1:NL-3,3)=O3_IN(1:NL-3)
      ! The -3 in the line above is just a fudge for the 23-layer model.
      ! Gavin said he'd think about how to do this properly.
                      CALL GETGAS
C--------------------------------


C--------------------------------
      SRBEXT=1.d-20 ; SRBSCT=0. ; SRBGCB=0. ; TRBALK=0.
      IF(MADBAK.GT.0) CALL GETBAK

      IF(MADAER.NE.0.OR.NTRACE.GT.0) THEN ; CALL GETAER
       ELSE ; SRAEXT=0.     ; SRASCT=0. ; SRAGCB=0. ; TRAALK=0. ; END IF
      IF(MADDST.GT.0) THEN ; CALL GETDST
       ELSE ; SRDEXT=0.     ; SRDSCT=0. ; SRDGCB=0. ; TRDALK=0. ; END IF
      IF(MADVOL.GT.0) THEN ; CALL GETVOL
       ELSE ; SRVEXT=0.     ; SRVSCT=0. ; SRVGCB=0. ; TRVALK=0. ; END IF
C--------------------------------


C--------------------------------  (GETSUR sets albedo needed by GETCLD)
                      CALL GETSUR
                      CALL GETEPS
                      CALL GETCLD
C--------------------------------

C--------------------------------
                      CALL THERML

                      CALL SOLARM
C--------------------------------

      RETURN
      END SUBROUTINE RCOMPX

      SUBROUTINE SETSOL(JYEARS,JJDAYS)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C     SETSOL Parameters:
C----------------------
C            KSOLAR    Selects Solar Spectrum, (Lean vs Thekaekara Flux)
C            JYEARS    JYEAR Proxy:  Sets: Solar Constant Reference Year
C            JJDAYS    JDAY  Proxy:  Sets Reference Year Month JDAY/30.5
C                      (Nominal Reference: JYEARS= 1950 JJDAYS= January)
C
C-----------------------------------------------------------------------
C            KSOLAR   SOLSPEC     UVWAVLs       UVFACTs         KUVFAC
C-----------------------------------------------------------------------
C              -1     THEK      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C               0     LEAN      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C               1     LEAN      Can be set    Can be set   (if KUVFAC=1)
C-----------------------------------------------------------------------
C
C                               (Option to Modify Solar UV Fluxes)
C            UVWAVL    Specified Edges of UV Flux Variation SubIntervals
C            UVFACT    Factors to Change the Amplitude of UV Variability
C
C            KUVFAC    ON/OFF switch for activating UV Flux Modification
C            KSNORM    Re-Normalize S0 (VIS) (after UV Amplitude Change)
C                         (Nominal UVWAVLs are: 0.295,0.310,0.366)
C
C-----------------------------------------------------------------------
C     SETSOL Output:
C------------------
C
C        AO3 =   Ozone Absorption Table AO3(460)
C                (Solar UV Flux Weighted Absorption Table is used by the
C                FUNCTION AO3ABS(OCM) in SOLAR to compute Ozone Heating)
C                AO3 is the fraction of total Solar Flux absorbed by O3.
C
C     S00WM2 =   Solar Constant Reference Value for Time = JYEARS,JJDAYS
C                (Thekaekara, if KSOLAR=-1, Reference = 1367 WATTS/M**2)
C
C
C     SETSOL  is Generally Called once at Model Initialization to Select
C                Solar Flux (LEAN,THEK), and to Define S00WM2 (RATLS0=1)
C
C-----------------------------------------------------------------------
C NOTE:
C-----
C     S00WM2 = Nominal Reference Solar Constant 1366.448785D0 WATTS/M**2
C                (Spectral Integral: Lean99 Solar Flux for January 1950)
C
C     KSOLAR=-1  Reproduces Thekaekhara Ozone Absorption, e.g., XRAD83XX
C     KSOLAR= 0  Uses Lean99 Solar Flux as set for Time= (JYEARS,JJDAYS)
C     KSOLAR= 1  Sets Lean99 Solar Flux to Current Time= (JYEARS,JJDAYS)
C     KSOLAR= 2  same as 1 but based on annual (not monthly) data
C                (JJDAYS used to select the specified Monthly-Mean Flux)
C
C-----------------------------------------------------------------------
C
C     UPDSOL Parameters:
C----------------------
C            JYEARS    JYEAR Proxy:  Selects Solar Constant Current Year
C            JJDAYS    JDAY  Proxy:  Selects Lean Data Month JJDAYS/30.5
C
C     UPDSOL Output:
C------------------
C
C        AO3 =   Ozone Absorption Table AO3(460)
C                (Solar UV Flux Weighted Absorption Table is used by the
C                FUNCTION AO3ABS(OCM) in SOLAR to compute Ozone Heating)
C                AO3 is the fraction of total Solar Flux absorbed by O3.
C
C     RATLS0 =   Ratio:  Current-Time Solar Constant to Reference S00WM2
C
C-----------------------------------------------------------------------
C     Remark:
C
C     UPDSOL  is Called in RCOMPT to Update Solar Constant and Ozone AO3
C                Solar UV Absorption Dependence.  (Monthly-Mean Data are
C                NOT Interpolated in Time, but get Updated with Changing
C                Month, i.e., whenever JDAY/30.5 Reaches Integer Value.)
C
C-----------------------------------------------------------------------
      REAL*8, PARAMETER :: CORFAC=1366.2911D0/1366.4487855D0
      INTEGER, INTENT(IN) :: JYEARS,JJDAYS
      INTEGER, SAVE :: LMOREF=0
      INTEGER JMO,LMO,Is0x,K,I,NWSUV,II,J
      REAL*8 FLXSUM,F1FLUX,F2FLUX,F3FLUX,UVNORM,XX,OCM
     *     ,TAUK,UVWAVA,UVWAVB,AO33

      IF(KSOLAR.GE.0) THEN
C                                           Lean99 Solar Flux, UV Option
C                                           ----------------------------
      if(Ksolar.lt.2) then    ! monthly data
        JMO=1+JJDAYS/30.5D0
        IF(JMO.GT.12) JMO=12
        LMO=(JYEARS-iy1S0)*12+JMO
        IF(LMO.GT.Ms0X) LMO=LMO-mcycs0*((LMO-Ms0X+mcycs0-1)/mcycs0)
        IF(LMO.LT.1) LMO=LMO+mcycs0*((mcycs0-lmo)/mcycs0)
      else                    ! annual data
        Is0x = nint( yr2s0-yr1s0+1 )
        lmo = nint( jyears - yr1s0 + 1.5 )
        IF(LMO.GT.Is0X) LMO=LMO-icycs0*((LMO-Is0X+icycs0-1)/icycs0)
        IF(LMO.LT.1) LMO=LMO+icycs0*((icycs0-lmo)/icycs0)
      end if
      LMOREF=LMO

C                        IF(MADLUV.EQ.0) Default Option is then in force
C                        Default (FRLEAN) = Lean 1950 Jan Solar, UV flux
C                        CORFAC accounts for DSLEAN units in BLOCK DATA,
C                        and TSI1/TSI2 normalization of Lean input data.
C                        -----------------------------------------------

c      CORFAC=1366.2911D0/1366.4487855D0
      FLXSUM=0.D0
      DO 110 K=1,190
      IF(MADLUV.EQ.0) FLXSUM=FLXSUM+FRLEAN(K)*DSLEAN(K)*CORFAC
      IF(MADLUV.GT.0) FLXSUM=FLXSUM+UVLEAN(LMO,K)*DSLEAN(K)
  110 CONTINUE
      S00WM2=FLXSUM

      I=0
      DO 120 K=1,50
      I=I+1
      WSOLAR(I)=W1LEAN(K)
      IF(MADLUV.EQ.0) FSOLAR(I)=FRLEAN(K)
      IF(MADLUV.GT.0) FSOLAR(I)=UVLEAN(LMO,K)
      I=I+1
      WSOLAR(I)=W1LEAN(K+1)
      FSOLAR(I)=FSOLAR(I-1)
  120 CONTINUE
      NWSUV=100
C                                          Thekaekhara Solar Flux Option
C                                          -----------------------------
      ELSE
      DO 130 I=1,190
      WSOLAR(I)=WTHEK(I)
      FSOLAR(I)=FTHEK(I)
  130 CONTINUE
      S00WM2=1367.D0
      LMOREF=-111
      NWSUV=190
      ENDIF

C                                         Option to Modify Solar UV Flux
C                                         ------------------------------
      IF(KUVFAC.EQ.1) THEN
      F1FLUX=0.D0
      F2FLUX=0.D0
      F3FLUX=0.D0
      DO 160 I=1,NWSUV
      IF(WSOLAR(I).GT.UVWAVL(1)) GO TO 140
      F1FLUX=F1FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(1)
      GO TO 160
  140 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(2)) GO TO 150
      F2FLUX=F2FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(2)
      GO TO 160
  150 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(3)) GO TO 170
      F3FLUX=F3FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(3)
  160 CONTINUE
  170 CONTINUE
      UVNORM=F1FLUX*(1.D0-UVFACT(1))+F2FLUX*(1.D0-UVFACT(2))
     +      +F3FLUX*(1.D0-UVFACT(3))
      IF(KSNORM.EQ.0) S00WM2=S00WM2-UVNORM
      ENDIF
C                  -----------------------------------------------------
C                  When KUVFAC=1 option multiplicative factors UVFACT(I)
C                  are used to change the UV spectral flux distribution,
C                  KSNORM=1 provides the option to keep S00WM2 constant.
C                  -----------------------------------------------------

      RATLS0=1.D0

      DO 190 I=1,460
      II=(I-10)/90-4
      XX=I-((I-10)/90)*90
      OCM=XX*10.D0**II
      DO 180 J=1,226
      TAUK=FUVKO3(J)*OCM
      IF(TAUK.GT.35.D0) TAUK=35.D0
      UVA(J)=1.D0-EXP(-TAUK)
  180 CONTINUE
      UVWAVA=0.100D0
      UVWAVB=0.400D0
      CALL FXGINT(UVA,XWAVO3,226,FSOLAR,WSOLAR,NWSUV,UVWAVA,UVWAVB,AO33)
      AO3(I)=AO33/S00WM2
  190 CONTINUE

C                       ------------------------------------------------
C                       NOTE:  AO3 is the Ozone-path Absorption Function
C                              AO3 convolves O3 asborption with solar UV
C                              spectral variations by FXGINT integration
C                              AO3 is expressed as the absorbed fraction
C                              of the total solar flux (S00WM2=1366W/m2)
C                              -----------------------------------------

      RETURN


C--------------------------------
      ENTRY UPDSOL(JYEARS,JJDAYS)
C--------------------------------


      IF(KSOLAR.LT.1) GO TO 300
      IF(JYEARS.LT.1) GO TO 300

      if(Ksolar.lt.2) then    ! monthly data       ksolar=1
        JMO=1+JJDAYS/30.5D0
        IF(JMO.GT.12) JMO=12
        LMO=(JYEARS-iy1S0)*12+JMO
        IF(LMO.GT.Ms0X) LMO=LMO-mcycs0*((LMO-Ms0X+mcycs0-1)/mcycs0)
        IF(LMO.LT.1) LMO=LMO+mcycs0*((mcycs0-lmo)/mcycs0)
      else                    ! annual data        ksolar=2
        Is0x = nint( yr2s0-yr1s0+1 )
        lmo = nint( jyears - yr1s0 + 1.5 )
        IF(LMO.GT.Is0X) LMO=LMO-icycs0*((LMO-Is0X+icycs0-1)/icycs0)
        IF(LMO.LT.1) LMO=LMO+icycs0*((icycs0-lmo)/icycs0)
      end if

      IF(LMO.EQ.LMOREF) GO TO 300
      LMOREF=LMO

      IF(KSOLAR.GE.0) THEN
C                                               Select Lean99 Solar Flux
C                                               ------------------------
      FLXSUM=0.D0
      DO 210 K=1,190
      FLXSUM=FLXSUM+UVLEAN(LMO,K)*DSLEAN(K)
  210 CONTINUE

      I=0
      DO 220 K=1,50
      I=I+1
      WSOLAR(I)=W1LEAN(K)
      FSOLAR(I)=UVLEAN(LMO,K)
      I=I+1
      WSOLAR(I)=W1LEAN(K+1)
      FSOLAR(I)=FSOLAR(I-1)
  220 CONTINUE
      NWSUV=100
      ELSE
C                                          Select Thekaekhara Solar Flux
C                                          -----------------------------
      DO 230 I=1,190
      WSOLAR(I)=WTHEK(I)
      FSOLAR(I)=FTHEK(I)
  230 CONTINUE
      FLXSUM=1367.D0
      LMOREF=-111
      NWSUV=190
      ENDIF
C                                         Option to Modify Solar UV Flux
C                                         ------------------------------
      IF(KUVFAC.EQ.1) THEN
      F1FLUX=0.D0
      F2FLUX=0.D0
      F3FLUX=0.D0
      DO 260 I=1,NWSUV
      IF(WSOLAR(I).GT.UVWAVL(1)) GO TO 240
      F1FLUX=F1FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(1)
      GO TO 260
  240 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(2)) GO TO 250
      F2FLUX=F2FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(2)
      GO TO 260
  250 CONTINUE
      IF(WSOLAR(I).GT.UVWAVL(3)) GO TO 270
      F3FLUX=F3FLUX+FSOLAR(I)
      FSOLAR(I)=FSOLAR(I)*UVFACT(3)
  260 CONTINUE
  270 CONTINUE
      UVNORM=F1FLUX*(1.D0-UVFACT(1))+F2FLUX*(1.D0-UVFACT(2))
     +      +F3FLUX*(1.D0-UVFACT(3))
      IF(KSNORM.EQ.0) FLXSUM=FLXSUM-UVNORM
      ENDIF

      RATLS0=FLXSUM/S00WM2

      DO 290 I=1,460
      II=(I-10)/90-4
      XX=I-((I-10)/90)*90
      OCM=XX*10.D0**II
      DO 280 J=1,226
      TAUK=FUVKO3(J)*OCM
      IF(TAUK.GT.35.D0) TAUK=35.D0
      UVA(J)=1.D0-EXP(-TAUK)
  280 CONTINUE
      UVWAVA=0.100D0
      UVWAVB=0.400D0
      CALL FXGINT(UVA,XWAVO3,226,FSOLAR,WSOLAR,NWSUV,UVWAVA,UVWAVB,AO33)
      AO3(I)=AO33/FLXSUM
  290 CONTINUE

  300 CONTINUE

      RETURN
      END SUBROUTINE SETSOL

      SUBROUTINE SETGHG(JYEARG,JJDAYG)
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C     SETGHG  Sets Default Greenhouse Gas Reference Year (for FULGAS)
C
C     Control Parameter:
C                     KTREND (specified in RADPAR) activates GH Trend
C               Default
C     KTREND  =   1
C     Selects   GTREND
C     ---------------------------------------------------------------
      INTEGER, INTENT(IN) :: JYEARG,JJDAYG
      REAL*8 TREF,TNOW
      INTEGER I
C
      TREF=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND.EQ.0) THEN
        XREF(1)=PPMV80(2)
        XREF(2)=PPMV80(6)
        XREF(3)=PPMV80(7)
        XREF(4)=PPMV80(8)*1000.D0
        XREF(5)=PPMV80(9)*1000.D0
        XREF(6)=PPMV80(11)*1000.D0  ! YREF11=PPMV80(11)*1000.D0
        XREF(7)=PPMV80(12)*1000.D0  ! ZREF12=PPMV80(12)*1000.D0
        RETURN
      END IF

      CALL GTREND(XREF,TREF)     ! finds xref 1-6 (yref11=xx6=xref(6))
      XREF(7)=1.D-13             ! ZREF12=1.D-13
      DO 120 I=1,NGHG
      IF(XREF(I).LT.1.D-06) XREF(I)=1.D-06
  120 CONTINUE
      PPMV80(2)=XREF(1)
      PPMV80(6)=XREF(2)
      PPMV80(7)=XREF(3)
      PPMV80(8)=XREF(4)/1000.D0
      PPMV80(9)=XREF(5)/1000.D0
      PPMV80(11)=XREF(6)/1000.d0   ! YREF11/1000.D0
      PPMV80(12)=XREF(7)/1000.D0   ! ZREF12/1000.D0
      RETURN
C
C--------------------------------
      ENTRY UPDGHG(JYEARG,JJDAYG)
C--------------------------------
C
      TNOW=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND.EQ.0) THEN
        FULGAS(2)=PPMVK0(2)/XREF(1)
        FULGAS(6)=PPMVK0(6)/XREF(2)
        FULGAS(7)=PPMVK0(7)/XREF(3)
        FULGAS(8)=PPMVK0(8)/XREF(4)
        FULGAS(9)=PPMVK0(9)/XREF(5)
        FULGAS(11)=PPMVK0(11)/XREF(6) ! YREF11
        FULGAS(12)=PPMVK0(12)/XREF(7) ! .../ZREF12
        RETURN
      END IF

      CALL GTREND(XNOW,TNOW) ! finds xnow 1-6 (ynow11=xx6=xnow(6))
      XNOW(7)=1.D-20         ! ZNOW12=1.D-20
      FULGAS(2)=XNOW(1)/XREF(1)
      FULGAS(6)=XNOW(2)/XREF(2)
      FULGAS(7)=XNOW(3)/XREF(3)
      FULGAS(8)=XNOW(4)/XREF(4)
      FULGAS(9)=XNOW(5)/XREF(5)
      FULGAS(11)=XNOW(6)/XREF(6) ! YNOW11/YREF11
      FULGAS(12)=XNOW(7)/XREF(7) ! ZNOW12/ZREF12
C
      RETURN
      END SUBROUTINE SETGHG

!     new Ozone routines

      SUBROUTINE UPDO3D(JYEARO,JJDAYO)


!!!   use RADPAR, only : IM=>MLON72,JM=>MLAT46,NL,PLB0,U0GAS,MADO3M
      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT NONE
      integer, parameter :: im=72,jm=46  ! ??? temporarily

!     In aug2003, 9 decadal files and an ozone trend data file
!     have been defined using 18-layer PLB pressure levels

      integer, PARAMETER ::
     *     NFO3x=9,      !  max. number of decadal ozone files used
!!!  *     NLO3=49,      !  number of layers of ozone data files
     *     IYIO3=1850, IYEO3=2050, ! beg & end year of O3 trend file
     *     LMONTR=12*(IYEO3-IYIO3+1) ! length of O3 trend file

      REAL*4 O3YEAR(MLON72,MLAT46,NLO3,0:12),OTREND(MLAT46,NLO3,LMONTR)
      REAL*4 O3ICMA(im,jm,NLO3,12),O3JCMA(im,jm,NLO3,12),A(im,jm)
      INTEGER LJTTRO(jm)
!!!   COMMON/O3JCOM/O3JDAY(NLO3,im,jm) ! to 'save' & for offline testing

C     UPDO3D CALLs GTREND to get CH4 to interpolate tropospheric O3
C     -----------------------------------------------------------------

      CHARACTER*80 TITLE
      logical, parameter :: qbinary=.true.   ; logical qexist

C**** The data statements below are only used if  MADO3M > -1
      CHARACTER*40, DIMENSION(NFO3X) :: DDFILE = (/
     1            'aug2003_o3_shindelltrop_72x46x49x12_1850'
     2           ,'aug2003_o3_shindelltrop_72x46x49x12_1890'
     3           ,'aug2003_o3_shindelltrop_72x46x49x12_1910'
     4           ,'aug2003_o3_shindelltrop_72x46x49x12_1930'
     5           ,'aug2003_o3_shindelltrop_72x46x49x12_1950'
     6           ,'aug2003_o3_shindelltrop_72x46x49x12_1960'
     7           ,'aug2003_o3_shindelltrop_72x46x49x12_1970'
     8           ,'aug2003_o3_shindelltrop_72x46x49x12_1980'
     9           ,'aug2003_o3_shindelltrop_72x46x49x12_1990'/)
      INTEGER, DIMENSION(NFO3X) :: IYEAR =
     *     (/1850,1890,1910,1930,1950,1960,1970,1980,1990/)
      INTEGER :: NFO3 = NFO3X
      CHARACTER*40 :: OTFILE ='aug2003_o3timetrend_46x49x2412_1850_2050'
      INTEGER :: IFILE=11            ! not used in GCM runs

!!!   REAL*8 :: PLBO3(NLO3+1) = (/ ! could be read off the titles
!!!  *       984d0, 934d0, 854d0, 720d0, 550d0, 390d0, 285d0, 210d0, 
!!!  *       150d0, 125d0, 100d0,  80d0,  60d0,  55d0,  50d0, 
!!!  *        45d0,  40d0,  35d0,  30d0,  25d0,  20d0,  15d0, 
!!!  *       10.d0,  7.d0,  5.d0,  4.d0,  3.d0,  2.d0,  1.5d0, 
!!!  *        1.d0,  7d-1,  5d-1,  4d-1,  3d-1,  2d-1,  1.5d-1, 
!!!  *        1d-1,  7d-2,  5d-2,  4d-2,  3d-2,  2d-2,  1.5d-2, 
!!!  *        1d-2,  7d-3,  5d-3,  4d-3,  3d-3,  1d-3,  1d-7/)
C**** LJTTRO(jm)    could be computed from PLBO3
      DATA LJTTRO/6*6,7*7,20*8,7*7,6*6/           !   Top of troposphere
      INTEGER, SAVE :: IYR=0, JYRNOW=0, IYRDEC=0, IFIRST=1, JYR

      save nfo3,iyear,ljttro,otrend,o3year
!!!   save plbo3
      save ddfile,ifile

      INTEGER :: JYEARO,JJDAYO
      INTEGER I,J,L,M,N,IY,JYEARX,jy,MI,MJ,MK,MN,NLT,ILON,JLAT
      REAL*8 WTTI,WTTJ, WTSI,WTSJ,WTMJ,WTMI, XMI,DSO3

      IF(IFIRST.EQ.1) THEN

      IF(MADO3M.lt.0) then
C****   Find O3 data files and fill array IYEAR from title(1:4)
        nfo3=0
        write(6,'(/a)') ' List of O3 files used:'
        do n=1,nfo3X    !  files have the generic names O3file_01,....
          ddfile(n)=' '
          write (ddfile(n),'(a7,i2.2)') 'O3file_',n
          inquire (file=ddfile(n),exist=qexist)
          if(.not.qexist) go to 10 !  exit
          call openunit (ddfile(n),ifile,qbinary)
          read(ifile) title
          call closeunit (ifile)
          write(6,'(a,a)') ' read O3 file, record 1: ',title
          read(title(1:4),*) IYEAR(n)
          nfo3=nfo3+1
        end do

   10   continue
        if(nfo3.eq.1) JYEARO=-IYEAR(1)
        if(nfo3.eq.0) call stop_model('updo3d: no Ozone files',255)
        OTFILE='O3trend '
      END IF

      IY=0
      IF(JYEARO.LT.0) THEN ! 1 year of O3 data is used in a cyclical way
        do n=1,nfo3        ! check whether we need the O3 trend array
           if(IYEAR(n).eq.-JYEARO) IY=n
        end do
      end if

      if(IY.EQ.0) then ! READ strat O3 time trend for strat O3 interpol.
        call openunit (OTFILE,ifile,qbinary)
        READ (IFILE) OTREND    ! strat O3 time-trend OTREND(JM,18,2412)
        call closeunit (ifile)
        if(MADO3M.lt.0) write(6,'(a,a)') ' read ',OTfile
      else
        call openunit (ddfile(IY),ifile,qbinary)
        DO 40 M=1,12
        DO 30 L=1,NLO3
        READ (IFILE) TITLE,A
        DO 20 J=1,JM         ! (decadal climatology) O3YEAR(IM,JM,18,12)
        DO 20 I=1,IM
           O3YEAR(I,J,L,M)=A(I,J)
   20   CONTINUE
   30   CONTINUE
   40   CONTINUE
        DO 50 L=1,NLO3
        DO 50 J=1,JM
        DO 50 I=1,IM
        O3YEAR(I,J,L,0)=O3YEAR(I,J,L,12)
   50   CONTINUE
        call closeunit (ifile)
        JYRNOW=-JYEARO  ! insures that O3YEAR is no longer changed
      end if

      IFIRST=0
      ENDIF

C     To time input data READs, JYEARX is set ahead of JYEARO by 15 days
C     ------------------------------------------------------------------
      IF(JYEARO.LT.0) THEN    !              ... except in cyclical case
        JYEARX=-JYEARO        ! Use fixed-year decadal climatology
        IYRDEC=JYEARX-1       ! force continuous cycles
      ELSE
        JYEARX=MIN(JYEARO+(JJDAYO+15)/366,IYEO3+1) ! +1 for continuity
      END IF                                       !         at Dec 15

      IF(JYEARX.EQ.JYRNOW) GO TO 500    ! Get O3JDAY from current O3YEAR

C****
C**** Get 13 months of O3 data O3YEAR starting with the leading December
C****
      do jy=1,nfo3                  ! nfo3 is at least 2, if we get here
        if(iyear(jy).gt.JYEARx) go to 100
      end do
      jy=nfo3
  100 if(jy.le.1) jy=2
      iy=jy-1
      IYR=IYEAR(IY)
      JYR=IYEAR(JY)

C**** Get first decadal file
      call openunit (ddfile(IY),ifile,qbinary)                     ! IYR
      DO 140 M=1,12
      DO 130 L=1,NLO3
      READ (IFILE) TITLE,A
      DO 120 J=1,JM
      DO 110 I=1,IM
      O3ICMA(I,J,L,M)=A(I,J)    !   IYR decadal data O3ICMA(IM,JM,18,12)
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
      call closeunit (ifile)

      IF(JYEARX.EQ.IYR.AND.IYRDEC.NE.JYEARX-1.AND.IY.GT.1) THEN
C        READ and use prior decadal file to define prior year December
C        (only when starting up with JYEARO=1890,1910,1930,...1980)

      call openunit (ddfile(IY-1),ifile,qbinary)                   ! KYR
      DO 240 M=1,12
      DO 230 L=1,NLO3
      READ(IFILE) TITLE,A
      DO 220 J=1,JM
      DO 210 I=1,IM
      O3JCMA(I,J,L,M)=A(I,J)   !  (prior decade) KYR temporary data file
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
      call closeunit (ifile)

C     Tropospheric & stratospheric ozone timetrend interpolation weights
C       Tropospheric ozone time variability is proportional to CH4 trend
C          Stratospheric ozone (above level LJTTRO(J)) is linear in time
C     ------------------------------------------------------------------

      CALL O3_WTS (IYIO3,LMONTR, IYR,IYEAR(IY-1), JYEARX-1,12,     ! in
     *             WTTI,WTTJ, WTSI,WTSJ, MI,MK,MN)                 ! out

      DO 290 J=1,JM
      NLT=LJTTRO(J)     !      NLT=LJTTRO(J) is top layer of troposphere
      DO 260 L=1,NLT
      DO 250 I=1,IM
      O3YEAR(I,J,L,0)=WTTI*O3ICMA(I,J,L,12)+WTTJ*O3JCMA(I,J,L,12)
      IF(O3YEAR(I,J,L,0).LT.0.) O3YEAR(I,J,L,0)=0.
  250 CONTINUE
  260 CONTINUE

C     DSO3 = add-on residual intra-decadal stratospheric O3 variability
C     ------------------------------------------------------------------
      DO 280 L=NLT+1,NLO3
      DSO3=OTREND(J,L,MN)-WTSI*OTREND(J,L,MI)-WTSJ*OTREND(J,L,MK)
      DO 270 I=1,IM
      O3YEAR(I,J,L,0)=WTSI*O3ICMA(I,J,L,12)+WTSJ*O3JCMA(I,J,L,12)+ DSO3
      IF(O3YEAR(I,J,L,0).LT.0.) O3YEAR(I,J,L,0)=0.
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
      IYRDEC=JYEARX    !   Set flag to indicate December data is current
      ENDIF

C**** Get next  decadal file
      call openunit (ddfile(JY),ifile,qbinary)                    !  JYR
      DO 340 M=1,12
      DO 330 L=1,NLO3
      READ(IFILE) TITLE,A
      DO 320 J=1,JM
      DO 310 I=1,IM
      O3JCMA(I,J,L,M)=A(I,J)    !   JYR decadal data O3JCMA(IM,JM,18,12)
  310 CONTINUE
  320 CONTINUE
  330 CONTINUE
  340 CONTINUE
      call closeunit (ifile)

      IF(JYEARX.eq.IYRDEC) GO TO 410    ! done with prior December

      IF(JYEARX.eq.IYRDEC+1) THEN      ! copy data from M=12 -> M=0
        DO L=1,NLO3                    ! normal non-start-up case
        DO J=1,JM
        DO I=1,IM
        O3YEAR(I,J,L,0)=O3YEAR(I,J,L,12)      !     DEC from prior year
        END DO
        END DO
        END DO
        IYRDEC=JYEARX   !  Set flag to indicate December data is current
      ELSE
C       Interpolate prior December from the decadal files - start-up
        CALL O3_WTS (IYIO3,LMONTR, IYR,JYR, JYEARX-1,12,        ! in
     *               WTTI,WTTJ, WTSI,WTSJ, MI,MJ,MN)            ! out

        DO 400 J=1,JM
        NLT=LJTTRO(J)     !    NLT=LJTTRO(J) is top layer of troposphere
        DO L=1,NLT
        DO I=1,IM
        O3YEAR(I,J,L,0)=WTTI*O3ICMA(I,J,L,12)+WTTJ*O3JCMA(I,J,L,12)
        IF(O3YEAR(I,J,L,0).LT.0.) O3YEAR(I,J,L,0)=0.
        END DO
        END DO

        DO L=NLT+1,NLO3
        DSO3=OTREND(J,L,MN)-WTSI*OTREND(J,L,MI)-WTSJ*OTREND(J,L,MJ)
        DO I=1,IM
        O3YEAR(I,J,L,0)=WTSI*O3ICMA(I,J,L,12)+WTSJ*O3JCMA(I,J,L,12)+DSO3
        IF(O3YEAR(I,J,L,0).LT.0.) O3YEAR(I,J,L,0)=0.
        END DO
        END DO
  400   CONTINUE
        IYRDEC=JYEARX  !   Set flag to indicate December data is current
      END IF

C            Fill in a full year of O3 data by interpolation
C            -----------------------------------------------
  410 CONTINUE
      CALL O3_WTS (IYIO3,LMONTR, IYR,JYR, JYEARX,0,           ! in
     *             WTTI,WTTJ, WTSI,WTSJ, MI,MJ,MN)            ! out

C            Tropospheric O3 interpolation is in proportion to CH4 trend
C                         ----------------------------------------------
      DO 490 M=1,12
      DO 480 J=1,JM
      NLT=LJTTRO(J)     !      NLT=LJTTRO(J) is top layer of troposphere
      DO 450 L=1,NLT
      DO 440 I=1,IM
      O3YEAR(I,J,L,M)=WTTI*O3ICMA(I,J,L,M)+WTTJ*O3JCMA(I,J,L,M)
      IF(O3YEAR(I,J,L,M).LT.0.) O3YEAR(I,J,L,M)=0.
  440 CONTINUE
  450 CONTINUE

C     DSO3 = add-on residual intra-decadal stratospheric O3 variability
C     ------------------------------------------------------------------
      DO 470 L=NLT+1,NLO3
      DSO3=OTREND(J,L,M+MN)-WTSI*OTREND(J,L,M+MI)-WTSJ*OTREND(J,L,M+MJ)
      DO 460 I=1,IM
      O3YEAR(I,J,L,M)=WTSI*O3ICMA(I,J,L,M)+WTSJ*O3JCMA(I,J,L,M) + DSO3
      IF(O3YEAR(I,J,L,M).LT.0.) O3YEAR(I,J,L,M)=0.
  460 CONTINUE
  470 CONTINUE
  480 CONTINUE
  490 CONTINUE
      JYRNOW=JYEARX

C****
C**** O3JDAY is interpolated daily from O3YEAR seasonal data via JJDAYO
C****

  500 CONTINUE
C     the formula below yields M near the middle of month M
      XMI=(JJDAYO+JJDAYO+31-(JJDAYO+15)/61+(JJDAYO+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYO
      WTMI=1.D0-WTMJ
      IF(MI.GT.11) MI=0
      MJ=MI+1
      DO 530 J=1,JM
      DO 520 I=1,IM
      DO 510 L=1,NLO3
      O3JDAY(L,I,J)=WTMI*O3YEAR(I,J,L,MI)+WTMJ*O3YEAR(I,J,L,MJ)
  510 CONTINUE
  520 CONTINUE
  530 CONTINUE
      RETURN

!!    ENTRY GETO3D (ILON,JLAT)
C

!!    CALL REPART(O3JDAY(1,ILON,JLAT),PLBO3,NLO3+1,U0GAS(1,3),PLB0,NL+1)

!!    RETURN
      END SUBROUTINE UPDO3D

      SUBROUTINE O3_WTS (IYI,MONTHS, IY1,IY2, IYX,MON,        !  input
     *                   WTTI,WTTJ, WTSI,WTSJ, MI,MJ,MN)      ! output
!@sum O3_WTS finds the weights needed for Ozone interpolation
!@auth A. Lacis/R. Ruedy
      implicit none
      integer IYI,MONTHS  ! first year, length of O3-trend data - input
      integer IY1,IY2     ! 2 distinct years with O3 data       - input
      integer IYX,MON     ! current year and month              - input
      integer MI, MJ, MN  ! indices for ozone trend array       - output
      real*8  WTTI,WTTJ   ! tropospheric weights                - output
      real*8  WTSI,WTSJ   ! stratospheric weights               - output

      real*8  GHGAS(6)     ! greenhouse gas conentrations (3=CH4)
      real*8  dYEAR,CH4IY1,CH4IY2,CH4NOW  ! dummies

C     Tropospheric O3 interpolation is in proportion to CH4 trend
C                    GTREND returns mid-year (annual mean) values
C     -----------------------------------------------------------
      CALL GTREND(GHGAS,IY1+.5d0)
      CH4IY1=GHGAS(3)

      CALL GTREND(GHGAS,IY2+.5d0)
      CH4IY2=GHGAS(3)

      CALL GTREND(GHGAS,IYX+.5d0)
      CH4NOW=GHGAS(3)

      WTTI=(CH4IY2-CH4NOW)/(CH4IY2-CH4IY1)      !  Trop O3 varies as CH4
      WTTJ=1.D0-WTTI

C     Strat O3 interpolation uses relative monthly variability in OTREND
C     ------------------------------------------------------------------
      DYEAR=IY2-IY1
      WTSI=(IY2-IYX)/DYEAR                      ! Strat O3 = time linear
      IF(WTSI.LT.0.D0) WTSI=0.D0
      IF(WTSI.GT.1.D0) WTSI=1.D0
      WTSJ=1.d0-WTSI

      MI=    (IY1-IYI)*12 + MON
      MJ=    (IY2-IYI)*12 + MON
      MN=MAX((IYX-IYI)*12,0)
      IF(MN.GT.MONTHS-12) MN=MONTHS-12
      MN=MN+MON
c     write(0,*) 'IYI,MON,IY1,IY2,IYX,MON',IYI,MONTHS,IY1,IY2,IYX,MON
c     write(0,*) 'MI,MJ,MN',MI,MJ,MN
c     write(0,*) 'WTTI,WTTJ, WTSI,WTSJ',WTTI,WTTJ, WTSI,WTSJ


      RETURN
      END SUBROUTINE O3_WTS

      SUBROUTINE SETGAS
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Global   U.S. (1976) Standard Atmosphere  P, T, Geo Ht  Parameters
C-----------------------------------------------------------------------

      REAL*8, PARAMETER, DIMENSION(36) ::
     *     P36 = (/
     $1.2000D+03, .9720D+03, .9445D+03, .9065D+03, .8515D+03, .7645D+03,
     $ .6400D+03, .4975D+03, .3695D+03, .2795D+03, .2185D+03, .1710D+03,
     $ .1250D+03, .8500D+02, .6000D+02, .4000D+02, .2500D+02, .1500D+02,
     $ .7500D+01, .4000D+01, .2500D+01, .1500D+01, .7500D+00, .4000D+00,
     $ .2500D+00, .1500D+00, .7810D-01, .4390D-01, .2470D-01, .1390D-01,
     $ .7594D-02, .3623D-02, .1529D-02, .7030D-03, .2059D-03, .0D0/),
     *     UFAC36 =(/
     $     0.800d0,0.800d0,0.800d0,0.750d0,0.750d0,0.750d0,0.700d0,
     *     0.750d0,0.846d0,0.779d0,0.892d0,0.886d0,0.881d0,0.875d0,
     *     0.870d0,0.846d0,0.840d0,0.902d0,0.880d0,0.775d0,0.796d0,
     *     0.842d0,0.866d0,0.861d0,0.821d0,0.903d0,1.264d0,1.732d0,
     *     2.000d0,1.701d0,1.609d0,1.478d0,1.253d0,1.372d0,1.571d0,
     *     1.571d0/)

      REAL*8, PARAMETER :: HPCON=34.16319d0,P0=1013.25d0,
     *     PI=3.141592653589793D0
      REAL*8, SAVE :: SINLAT(46)
      INTEGER, SAVE :: IFIRST=1, NL0,NL1
      INTEGER I,NLAY,NATM,L,J,K,N
      REAL*8 RADLAT,RHP,EST,FWB,FWT,PLT,DP,EQ,ES,ACM,HI,FI,HL,HJ,FJ,DH
     *     ,FF,GGVDF,ZT,ZB,EXPZT,EXPZB,PARTTR,PARTTG,PTRO,DL,DLS,DLN
     *     ,Z0LAT,SUMCOL,ULGASL

      IF(IFIRST.EQ.1) THEN
      DO 90 I=1,MLAT46
      RADLAT=DLAT46(I)*PI/180.D0
      SINLAT(I)=SIN(RADLAT)
   90 CONTINUE
      NL0=NL
      IFIRST=0
      ENDIF
C                  -----------------------------------------------------
C                  Use PLB to fix Standard Heights for Gas Distributions
C                  -----------------------------------------------------

      NL1=NL0+1
!nu   PS0=PLB0(1)

      DO 100 L=1,NL0
      DPL(L)=PLB0(L)-PLB0(L+1)
      PL(L)=(PLB0(L)+PLB0(L+1))*0.5D0
!nu   HLB(L)=HLB0(L)
  100 CONTINUE
!nu   HLB(NL1)=HLB0(NL1)
      CALL RETERP(UFAC36,P36,36,FPXCO2,PL,NL0)
cc    IUFAC=1
cc    IF(IUFAC.EQ.0) THEN
cc      DO I=1,NL
cc        FPXCO2(I)=1.
cc      END DO
cc    ENDIF

      NLAY=LASTVC/100000
      NATM=(LASTVC-NLAY*100000)/10000
      IF(NATM.GT.0) GO TO 112

C     ----------------------------------------------------------------
C     Define Default Global Mean Gas Amounts for Off-Line Use Purposes
C
C                                         Global Mean H2O Distribution
C                                         ----------------------------
      RHP=0.77D0
      EST=10.D0**(9.4051D0-2353.D0/TLB(1))
      FWB=0.662D0*RHP*EST/(PLB0(1)-RHP*EST)
      DO 111 L=1,NL0
      PLT=PLB0(L+1)
      DP=PLB0(L)-PLT
      RHP=0.77D0*(PLT/P0-0.02D0)/.98D0
      EST=10.D0**(9.4051D0-2353.D0/TLT(L))
      FWT=0.662D0*RHP*EST/(PLT-RHP*EST)
      IF(FWT.GT.3.D-06) GO TO 110
      FWT=3.D-06
      RHP=FWT*PLT/(EST*(FWT+0.662D0))
  110 CONTINUE
      ULGASL=0.5D0*(FWB+FWT)*DP*1268.75D0
      U0GAS(L,1)=ULGASL
      SHL(L)=ULGASL/(ULGASL+1268.75D0*DP)
      EQ=0.5D0*(PLB0(L)+PLT)*SHL(L)/(0.662D0+0.378D0*SHL(L))
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      RHL(L)=EQ/ES
      FWB=FWT
  111 CONTINUE
  112 CONTINUE

C                                         ----------------------------
C                                         Global Mean NO2 Distribution
C                                         ----------------------------
      ACM=0.D0
      HI=0.D0
      FI=CMANO2(1)
      HL=HLB0(2)
      L=1
      J=1
  130 CONTINUE
      J=J+1
      IF(J.GT.42) GO TO 133
      HJ=HI+2.D0
      FJ=CMANO2(J)
  131 CONTINUE
      DH=HJ-HI
      IF(HJ.GT.HL) GO TO 132
      ACM=ACM+(FI+FJ)*DH*0.5D0
      HI=HJ
      FI=FJ
      GO TO 130
  132 CONTINUE
      FF=FI+(FJ-FI)*(HL-HI)/DH
      DH=HL-HI
      ACM=ACM+(FI+FJ)*DH*0.5D0
      U0GAS(L,5)=ACM
      ACM=0.D0
      HI=HL
      FI=FF
      IF(L.EQ.NL0) GO TO 133
      L=L+1
      HL=HLB0(L+1)
      GO TO 131
  133 CONTINUE
      U0GAS(L,5)=ACM
      ACM=0.D0
      L=L+1
      IF(L.LT.NL1) GO TO 133
C                            -----------------------------------------
C                            (CO2,O2) Uniformly Mixed Gas Distribution
C                            -----------------------------------------
      DO 141 K=2,4,2
      DO 140 N=1,NL0
      U0GAS(N,K)=PPMV80(K)*0.8D0*DPL(N)/P0
  140 CONTINUE
  141 CONTINUE
C                -----------------------------------------------------
C                (N20,CH4,F11,F12) Specified Vertical Gas Distribution
C                -----------------------------------------------------
      DO 151 K=6,12
      IF(K.EQ.10) GO TO 151
      DO 150 N=1,NL0
      GGVDF=1.D0-(1.D0-PPMVDF(K))*(1.D0-PLB0(N)/PLB0(1))
      IF(KGGVDF.LT.1) GGVDF=1.D0
      U0GAS(N,K)=PPMV80(K)*0.8D0*DPL(N)/P0*GGVDF
      ZT=(HLB0(N+1)-Z0(K))/ZH(K)
      IF(ZT.LE.0.D0) GO TO 150
      ZB=(HLB0(N)-Z0(K))/ZH(K)
      EXPZT=EXP(-ZT)
      EXPZB=EXP(-ZB)
      IF(ZB.LT.0.D0) EXPZB=1.D0-ZB
      U0GAS(N,K)=U0GAS(N,K)*(EXPZB-EXPZT)/(ZT-ZB)
  150 CONTINUE
  151 CONTINUE
C                         --------------------------------------------
C                         Specification of  FULGAS  Scaled Gas Amounts
C                         --------------------------------------------

      DO 230 L=1,NL0
Cc*** Adjust all water levels
cc    ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
c**** Only adjust stratospheric water levels (above LS1_loc)
      IF (L.lt.LS1_loc) THEN
        ULGAS(L,1)=U0GAS(L,1)
      ELSE
        ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
      END IF
C****
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3)
!obso ULGAS(L,5)=U0GAS(L,5)*FULGAS(5) ! redefined below
  230 CONTINUE

      IF(KPFOZO.EQ.1) THEN
      DO 235 L=1,NL0
      ULGAS(L,3)=ULGAS(L,3)*FPXOZO(L)
  235 CONTINUE
      ENDIF

      PARTTR=(PLB(1)-PTOPTR)/(PLB0(1)-PTOPTR)
      DO 240 L=1,NL0
      IF(PLB(L).LE.PTOPTR) PARTTR=1.D0
      DO 239 K=2,12
      IF(K.EQ.3) GO TO 239
      PARTTG=PARTTR
!nojl IF(KPGRAD.GT.0) PARTTG=PARTTG*(1.D0+0.5D0*PPGRAD(K)*SINLAT(JLAT))
      ULGAS(L,K)=U0GAS(L,K)*FULGAS(K)*PARTTG
  239 CONTINUE
      ULGAS(L,13)=U0GAS(L,13)*FULGAS(13)
  240 CONTINUE

      IF(KPFCO2.EQ.1) THEN
      DO 245 L=1,NL0
      ULGAS(L,2)=ULGAS(L,2)*FPXCO2(L)
  245 CONTINUE
      ENDIF

      RETURN


C-----------------
      ENTRY GETGAS
C-----------------
C                        ---------------------------------------------
C                        Specify ULGAS: Get Gas Absorption from TAUGAS
C                        ---------------------------------------------

C                -----------------------------------------------------
C                N20,CH4,F11,F12 Specified Latitudinal Z0 Distribution
C                -----------------------------------------------------

      IF(KLATZ0.GT.0) THEN
      PTRO=100.D0
      DL=DLAT46(JLAT)
      DLS=-40.D0
      DLN= 40.D0
      IF(DL.LT.DLS) PTRO=189.D0-(DL+40.D0)*2.22D0
      IF(DL.GT.DLN) PTRO=189.D0+(DL-40.D0)*2.22D0
      DO 249 N=1,NL0
      IF(PLB0(N).GE.PTRO) Z0LAT=HLB0(N)  ! hlb not hlb0
  249 CONTINUE
      DO 251 K=6,12
      IF(K.EQ.10) GO TO 251
      DO 250 N=1,NL0
      U0GAS(N,K)=PPMV80(K)*0.8D0*(PLB0(N)-PLB0(N+1))/P0
      ZT=(HLB0(N+1)-Z0LAT)/ZH(K)         ! hlb not hlb0
      IF(ZT.LE.0.D0) GO TO 250
      ZB=(HLB0(N)-Z0LAT)/ZH(K)           ! hlb not hlb0
      EXPZT=EXP(-ZT)
      EXPZB=EXP(-ZB)
      IF(ZB.LT.0.D0) EXPZB=1.D0-ZB
      U0GAS(N,K)=U0GAS(N,K)*(EXPZB-EXPZT)/(ZT-ZB)
  250 CONTINUE
  251 CONTINUE
      ENDIF

      DO 300 L=1,NL
      DPL(L)=PLB(L)-PLB(L+1)
      PL(L)=(PLB(L)+PLB(L+1))*0.5D0
  300 CONTINUE

      IF(KEEPRH.EQ.1) GO TO 311
      DO 310 L=1,NL
      U0GAS(L,1)=1268.75D0*DPL(L)*SHL(L)/(1.D0-SHL(L))
      EQ=PL(L)*SHL(L)/(0.662D0+0.378D0*SHL(L))
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      RHL(L)=EQ/ES
  310 CONTINUE
      GO TO 313
  311 CONTINUE
      DO 312 L=1,NL
      ES=10.D0**(9.4051D0-2353.D0/TLM(L))
      SHL(L)=0.622D0*(RHL(L)*ES)/(PL(L)-0.378D0*(RHL(L)*ES))
      U0GAS(L,1)=1268.75D0*DPL(L)*SHL(L)/(1.D0-SHL(L))
  312 CONTINUE
  313 CONTINUE
      DO 314 L=1,NL
Cc*** Adjust ALL water levels
cc    ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
c**** Only adjust stratospheric water levels (above LS1_loc)
      IF (L.lt.LS1_loc) THEN
        ULGAS(L,1)=U0GAS(L,1)
      ELSE
        ULGAS(L,1)=U0GAS(L,1)*FULGAS(1)
      END IF
C****
  314 CONTINUE

      DO 330 L=1,NL0
      ULGAS(L,3)=U0GAS(L,3)*FULGAS(3)
!obso ULGAS(L,5)=U0GAS(L,5)*FULGAS(5) ! redefined below
  330 CONTINUE

      IF(KPFOZO.EQ.1) THEN
      DO 335 L=1,NL0
      ULGAS(L,3)=ULGAS(L,3)*FPXOZO(L)
  335 CONTINUE
      ENDIF

      PARTTR=(PLB(1)-PTOPTR)/(PLB0(1)-PTOPTR)
      DO 340 L=1,NL0
      IF(PLB(L).LE.PTOPTR) PARTTR=1.D0
      DO 339 K=2,12
      IF(K.EQ.3) GO TO 339
      PARTTG=PARTTR
      IF(KPGRAD.GT.0) PARTTG=PARTTG*(1.D0+0.5D0*PPGRAD(K)*SINLAT(JLAT))
      ULGAS(L,K)=U0GAS(L,K)*FULGAS(K)*PARTTG
  339 CONTINUE
      ULGAS(L,13)=U0GAS(L,13)*FULGAS(13)
  340 CONTINUE

      IF(KPFCO2.EQ.1) THEN
      DO 345 L=1,NL0
      ULGAS(L,2)=ULGAS(L,2)*FPXCO2(L)
  345 CONTINUE
      ENDIF

      IF(MRELAY.GT.0) THEN          ! for offline use only if NL changes
        IF(NO3COL.GT.0) THEN        ! rescale ozone to col.amount RO3COL
          SUMCOL=0.D0
          DO L=1,NL0
          SUMCOL=SUMCOL+U0GAS(L,3)
          END DO
          DO L=1,NL0
          ULGAS(L,3)=U0GAS(L,3)*RO3COL/SUMCOL
          END DO
        ENDIF
        DO 450 K=2,12                      ! repartition to new layering
          IF(K.EQ.10.AND.KEEP10.GT.0) GO TO 450
          DO L=1,NL0
          UGAS0(L)=ULGAS(L,K)
          END DO
          CALL REPART(UGAS0,PLB0,NL1,UGASR,PLB,NLP)
          DO L=1,NL
          ULGAS(L,K)=UGASR(L)
          END DO
  450   CONTINUE
        IF(KEEP10.GT.0) THEN
          IF(KEEP10.LT.10) THEN
            DO L=1,NL
            ULGAS(L,KEEP10)=ULGAS(L,10)
            END DO
          ENDIF
          IF(KEEP10.GT.10) THEN
            DO L=1,NL
            ULGAS(L,KEEP10-10)=ULGAS(L,KEEP10-10)+ULGAS(L,10)
            END DO
          ENDIF
        ENDIF
      ENDIF

C-----------------
      CALL  TAUGAS
C-----------------

      RETURN
      END SUBROUTINE SETGAS

      SUBROUTINE SETO2A
      IMPLICIT NONE

      REAL*8, PARAMETER ::
     *     SFWM2(18) = (/
     A 2.196E-03, 0.817E-03, 1.163E-03, 1.331E-03, 1.735E-03, 1.310E-03,
     B 1.311E-03, 2.584E-03, 2.864E-03, 4.162E-03, 5.044E-03, 6.922E-03,
     C 6.906E-03,10.454E-03, 5.710E-03, 6.910E-03,14.130E-03,18.080E-03
     *     /),
     *     SIGMA(18,6) = RESHAPE( (/
     A     2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19, 2.74E-19,
     B     4.33E-21, 4.89E-21, 6.63E-21, 1.60E-20, 7.20E-20, 1.59E-18,
     C     2.10E-21, 2.32E-21, 3.02E-21, 6.30E-21, 3.46E-20, 7.52E-19,
     D     5.95E-22, 9.72E-22, 2.53E-21, 7.57E-21, 7.38E-20, 7.44E-19,
     E     3.33E-22, 1.02E-22, 4.09E-21, 1.63E-20, 8.79E-20, 3.81E-19,
     F     1.09E-21, 1.16E-21, 1.45E-21, 3.32E-21, 2.00E-20, 4.04E-19,
     G     1.15E-21, 1.30E-21, 1.90E-21, 4.89E-21, 2.62E-20, 4.08E-19,
     H     3.90E-22, 4.90E-22, 9.49E-22, 3.33E-21, 2.14E-20, 2.39E-19,
     I     1.29E-22, 2.18E-22, 8.28E-22, 3.46E-21, 1.94E-20, 1.06E-19,
     J     6.26E-23, 7.80E-23, 2.62E-22, 1.83E-21, 1.25E-20, 3.95E-20,
     K     2.74E-23, 3.58E-23, 8.64E-23, 4.03E-22, 2.13E-21, 1.95E-20,
     L     1.95E-23, 2.44E-23, 4.89E-23, 2.87E-22, 1.95E-21, 1.36E-20,
     M     1.84E-23, 1.96E-23, 2.71E-23, 8.52E-23, 6.48E-22, 3.89E-21,
     N     1.80E-23, 1.81E-23, 1.87E-23, 2.69E-23, 1.34E-22, 1.52E-21,
     O     1.80E-23, 1.80E-23, 1.82E-23, 2.40E-23, 5.71E-23, 5.70E-22,
     P     1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 1.76E-23, 3.50E-23,
     Q     1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 1.71E-23, 2.68E-23,
     R     1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23, 1.00E-23/)
     *     , (/ 18, 6 /) )

      REAL*8, PARAMETER, DIMENSION(6) :: WTKO2 =
     *     (/0.05,0.20,0.25,0.25,0.20,0.05/)
      INTEGER, SAVE :: NL0,NL1

      REAL*8, PARAMETER :: STPMOL=2.68714D+19
      INTEGER, PARAMETER :: NW=18, NZ=11, NKO2=6
      INTEGER, SAVE :: IFIRST=1
      REAL*8 FSUM,SUMMOL,ZCOS,WSUM,TAU,DLFLUX,WTI,WTJ
      INTEGER I,J,K,L,JI,JJ,N,LL

      IF(IFIRST.EQ.1) THEN
      NL0=NL
      DO 10 N=1,NL0
      ULGAS(N,4)=PPMV80(4)*0.8D0*(PLB0(N)-PLB0(N+1))/PLB0(1)
   10 CONTINUE
      IFIRST=0
      ENDIF

      NL1=NL0+1
      FSUM=0.0
      DO 100 I=1,NW
      FSUM=FSUM+SFWM2(I)
  100 CONTINUE
      DO 110 J=1,NZ
      ZTABLE(NL1,J)=FSUM
  110 CONTINUE
      SUMMOL=0.D0
      DO 150 N=1,NL0
      L=NL1-N
      SUMMOL=SUMMOL+ULGAS(L,4)*STPMOL
      DO 140 J=1,NZ
      ZCOS=0.01D0*(1/J)+0.1D0*(J-1)
      FSUM=0.D0
      DO 130 I=1,NW
      WSUM=0.D0
      DO 120 K=1,NKO2
      TAU=SIGMA(I,K)*SUMMOL/ZCOS
      IF(TAU.GT.30.D0) TAU=30.D0
      WSUM=WSUM+WTKO2(K)*EXP(-TAU)
  120 CONTINUE
      FSUM=FSUM+WSUM*SFWM2(I)
  130 CONTINUE
      ZTABLE(L,J)=FSUM
  140 CONTINUE
  150 CONTINUE
      DO 170 J=1,NZ
      DO 160 L=1,NL0
      DLFLUX=ZTABLE(L+1,J)-ZTABLE(L,J)
      ZTABLE(L,J)=DLFLUX/1366.D0
  160 CONTINUE
  170 CONTINUE
      RETURN

C-----------------
      ENTRY GETO2A
C-----------------

C              ---------------------------------------------------------
C              UV absorption by Oxygen is expressed as a fraction of the
C              total solar flux S0. Hence, O2FHRL(L)=ZTABLE(L,J) must be
C              normalized within SOLARM, dividing the GETO2A absorptions
C              O2FHRL(L) and O2FHRB(L) by the fraction of the solar flux
C              within the spectral interval DKS0(15), nominally by 0.05.
C              ---------------------------------------------------------

      ZCOS=1.D0+10.D0*COSZ
      JI=ZCOS
      IF(JI.GT.10) JI=10
      JJ=JI+1
      WTJ=ZCOS-JI
      WTI=1.0-WTJ
      DO 200 LL=1,NL
      L=NL1-LL
      IF(L.LT.1) GO TO 210
      O2FHRL(L)=(WTI*ZTABLE(L,JI)+WTJ*ZTABLE(L,JJ))
      O2FHRB(L)=ZTABLE(L,6)
  200 CONTINUE
  210 CONTINUE
      RETURN
      END SUBROUTINE SETO2A

      SUBROUTINE SETBAK
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SETBAK,GETBAK  Initializes Background Aerosol Specification, i.e.,
C                    Aerosol Composition and Distribution that is set in
C                    RADPAR by AGOLDH, BGOLDH, CGOLDH Factors
C                    and controlled by FGOLDH ON/OFF Scaling Parameters.
C                    Optional tracers may be added in SETAER/GETAER
C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSBAER   FTBAER  ! Bgrnd Aerosol
C
C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSBAER    SW   SETBAKonly Aerosol Optical Depth
C                        FTBAER    LW   SETBAKonly Aerosol Optical Depth
C                        -----------------------------------------------
      REAL*8, DIMENSION(5) :: SGOLDH,TGOLDH

      INTEGER, SAVE :: IFIRST=1
      INTEGER, SAVE :: NL0=0

      REAL*8 C,BC,ABC,HXPB,HXPT,ABCD,SRAQX,SRAQS
     *      ,SUMABS,EXTSUM,SCTSUM,COSSUM
      INTEGER I,J,K,L,JJ

C**** Background aerosols
C     ------------------------------------------------------------------
C     Thermal: Set (5) Aerosol Type Compositions & Vertical Distribution
C     ------------------------------------------------------------------
      IF(IFIRST.EQ.1) THEN
      NL0=NL
      IFIRST=0
      ENDIF

      DO 100 J=1,5
      DO 100 K=1,33
      DO 100 L=1,NL0
      TRAX(L,K,J)=0.D0
  100 CONTINUE

      DO 105 I=1,11
      DO 103 J=1,5
      IF(AGOLDH(I,J).LT.1.D-06) GO TO 103
      C=CGOLDH(I,J)
      BC=EXP(-BGOLDH(I,J)/C)
      ABC=AGOLDH(I,J)*(1.D0+BC)

      HXPB=1.D0
      DO 102 L=1,NL0
      HXPT=HLB0(L+1)/C
      IF(HXPT.GT.80.D0) GO TO 102
      HXPT=EXP(HXPT)
      ABCD=ABC/(1.D0+BC*HXPB)
     +    -ABC/(1.D0+BC*HXPT)
      HXPB=HXPT
      DO 101 K=1,33
      TRAX(L,K,J)=TRAX(L,K,J)+ABCD*(TRAQEX(K,I)-TRAQSC(K,I))
  101 CONTINUE
  102 CONTINUE
  103 CONTINUE
      DO 104 K=1,33
      TRAQAB(K,I)=TRAQEX(K,I)-TRAQSC(K,I)
  104 CONTINUE
  105 CONTINUE

      DO 107 K=1,33
      DO 106 L=1,NL0
      TRBALK(L,K)=0.D0
  106 CONTINUE
  107 CONTINUE

C-----------------------------------------------------------------------
C     SOLAR:   Set (5) Aerosol Type Compositions & Vertical Distribution
C-----------------------------------------------------------------------

      DO 110 J=1,5
      DO 110 K=1,6
      DO 110 L=1,NL0
      SRAX(L,K,J)=1.D-20
      SRAS(L,K,J)=1.D-30
      SRAC(L,K,J)=0.D0
  110 CONTINUE

      DO 114 I=1,11
      DO 113 J=1,5
      IF(AGOLDH(I,J).LT.1.D-06) GO TO 113
      C=CGOLDH(I,J)
      BC=EXP(-BGOLDH(I,J)/C)
      ABC=AGOLDH(I,J)*(1.D0+BC)

      HXPB=1.D0
      DO 112 L=1,NL0
      HXPT=HLB0(L+1)/C
      IF(HXPT.GT.80.D0) GO TO 112
      HXPT=EXP(HXPT)
      ABCD=ABC/(1.D0+BC*HXPB)
     +    -ABC/(1.D0+BC*HXPT)
      HXPB=HXPT
      DO 111 K=1,6
      SRAQX=SRAQEX(K,I)
      SRAQS=SRAQSC(K,I)
      SRAX(L,K,J)=SRAX(L,K,J)+ABCD*SRAQX
      SRAS(L,K,J)=SRAS(L,K,J)+ABCD*SRAQS
      SRAC(L,K,J)=SRAC(L,K,J)+ABCD*SRAQCB(K,I)*SRAQS
  111 CONTINUE
  112 CONTINUE
  113 CONTINUE
  114 CONTINUE

      DO 117 J=1,5
      DO 116 K=1,6
      DO 115 L=1,NL0
      SRAC(L,K,J)=SRAC(L,K,J)/SRAS(L,K,J)
  115 CONTINUE
  116 CONTINUE
  117 CONTINUE

      DO 119 K=1,6
      DO 118 L=1,NL0
      SRBEXT(L,K)=1.D-20
      SRBSCT(L,K)=0.D0
      SRBGCB(L,K)=0.D0
!nu   SRBPI0(L,K)=0.D0
  118 CONTINUE
  119 CONTINUE

      RETURN

C-----------------
      ENTRY GETBAK
C-----------------

C     ------------------------------------------------------------------
C     GETBAK   Specifies Background Aerosol Contribution and Initializes
C                    (1) Thermal Radiation Aerosol Coefficient Table:
C                        TRAALK(L,K), for (L=1,NL), (K=1,33)
C
C                    (2) Solar Radiation Coefficient Tables:
C                        SRAEXT(L,K),SRASCT(L,K),SRAGCB(L,K) for (K=1,6)
C                    ---------------------------------------------------

C                                                              (Thermal)
C                                                              ---------
      DO 200 J=1,5
      TGOLDH(J)=FTTAER*FTBAER*FGOLDH(J)
  200 CONTINUE
      DO 203 L=1,NL0
      DO 202 K=1,33
      SUMABS=1.D-20
      DO 201 J=1,5
      SUMABS=SUMABS+TGOLDH(J)*TRAX(L,K,J)
  201 CONTINUE
      TRBALK(L,K)=SUMABS
  202 CONTINUE
  203 CONTINUE

C                                                                (Solar)
C                                                                -------

      DO 210 J=1,5
      SGOLDH(J)=FSTAER*FSBAER*FGOLDH(J)
  210 CONTINUE
      DO 212 K=1,6
      DO 212 L=1,NL0
      EXTSUM=1.D-20
      SCTSUM=1.D-30
      COSSUM=0.D0
      DO 211 J=1,5
      EXTSUM=EXTSUM+SGOLDH(J)*SRAX(L,K,J)
      SCTSUM=SCTSUM+SGOLDH(J)*SRAS(L,K,J)
      COSSUM=COSSUM+SGOLDH(J)*SRAS(L,K,J)*SRAC(L,K,J)
  211 CONTINUE
      SRBEXT(L,K)=EXTSUM
      SRBSCT(L,K)=SCTSUM
      SRBGCB(L,K)=COSSUM/SCTSUM
  212 CONTINUE


      RETURN
      END SUBROUTINE SETBAK


      SUBROUTINE SETAER

cc    INCLUDE  'rad00def.radCOMMON.f'
      IMPLICIT NONE
C     ---------------------------------------------------------------
C     GISS MONTHLY-MEAN (1850-2050)  TROPOSPHERIC AEROSOL CLIMATOLOGY
C     ---------------------------------------------------------------

C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSAAER   FTAAER  ! AClim Aerosol

C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSAAER    SW   AClim Aer  Aerosol Optical Depth
C                        FTAAER    LW   AClim Aer  Aerosol Optical Depth
C                        -----------------------------------------------

!nu   DIMENSION ATAU09(9)
cc    DIMENSION PLBA09(10)          !       Aerosol data pressure levels
cc    DATA PLBA09/984.,934.,854.,720.,550.,390.,255.,150.,70.,10./
      real*8, parameter, dimension(4) ::
C              Crystallization RH               Deliquescence RH
     *  RHC=(/.38d0,.47d0,.28d0,.38d0/), RHD=(/.80d0,.75d0,.62d0,.80d0/)

C     ------------------------------------------------------------------
C                  Define aerosol size according to REFDRY specification
C                                      (if KRHAER(NA)=0, REFWET is used)
C                  FRSULF= Sulfate fraction of basic aerosol composition
C
C          Set size SO4 (NA=1) = Sulfate aerosol  (Nominal dry Reff=0.2)
C          Set size SEA (NA=2) = SeaSalt aerosol  (Nominal dry Reff=1.0)
C          Set size ANT (NA=3) = Nitrate aerosol  (Nominal dry Reff=0.3)
C          Set size OCX (NA=4) = Organic aerosol  (Nominal dry Reff=0.3)
C     ------------------------------------------------------------------
      REAL*8 AREFF, XRH,FSXTAU,FTXTAU,SRAGQL,RHFTAU,q55
      REAL*8          TTAULX(LX,8),   SRBGQL
      INTEGER K,L,NA,N,NRH,M,KDREAD,NT

      IF(MADAER.LE.0) GO TO 150
      DO 110 NA=1,4
      AREFF=REFDRY(NA)
      IF(KRHAER(NA).EQ.0) AREFF=REFWET(NA)
      CALL GETMIE(NA,AREFF,SRHQEX(1,1,NA),SRHQSC(1,1,NA),SRHQCB(1,1,NA)
     +                    ,TRHQAB(1,1,NA),Q55DRY(NA))
      DRYM2G(NA)=0.75D0/DENAER(NA)*Q55DRY(NA)/AREFF
      RHINFO(1,1,NA)=0.D0                                     !  Rel Hum
      RHINFO(1,2,NA)=1.D0                                     !  TAUFAC
      RHINFO(1,3,NA)=AREFF                                    !  AerSize
      RHINFO(1,4,NA)=0.D0                                     !  LW g/m2
      RHINFO(1,5,NA)=1.33333333D0*AREFF*DENAER(NA)/Q55DRY(NA) !  Dryg/m2
      RHINFO(1,6,NA)=1.33333333D0*AREFF*DENAER(NA)/Q55DRY(NA) !  Totg/m2
      RHINFO(1,7,NA)=1.D0                                     !  Xmas fr
      RHINFO(1,8,NA)=DENAER(NA)                               !  Density
      RHINFO(1,9,NA)=Q55DRY(NA)                               !  Q55 Ext
  110 CONTINUE

C     Set size BCI (NA=5) = Black Carbon (Industrial) (Nominal Reff=0.1)
C     Set size BCB (NA=6) = Black Carbon (BioBurning) (Nominal Reff=0.1)
C     ------------------------------------------------------------------
      DO 120 NA=5,6
      AREFF=REFDRY(NA)
      CALL GETMIE(NA,AREFF,SRBQEX(1,NA),SRBQSC(1,NA),SRBQCB(1,NA)
     +                    ,TRBQAB(1,NA),Q55DRY(NA))
      DRYM2G(NA)=0.75D0/DENAER(NA)*Q55DRY(NA)/AREFF
  120 CONTINUE

              !      Extend default dry aerosol coefficients for N=2,190
      DO 135 N=2,190
      DO 134 NA=1,4
      DO 131 K=1,6
      SRHQEX(K,N,NA)=SRHQEX(K,1,NA)
      SRHQSC(K,N,NA)=SRHQSC(K,1,NA)
      SRHQCB(K,N,NA)=SRHQCB(K,1,NA)
  131 CONTINUE
      DO 132 K=1,33
      TRHQAB(K,N,NA)=TRHQAB(K,1,NA)
  132 CONTINUE
      DO 133 M=1,9
      RHINFO(N,M,NA)=RHINFO(N,1,NA)
  133 CONTINUE
  134 CONTINUE
  135 CONTINUE

      KDREAD=71           !  Over-write dry coefficients if KRHAER(NA)=1
      DO 140 NA=1,4
      IF(KRHAER(NA).GT.0) THEN
      CALL SETREL(REFDRY(NA),NA ,kdread
     A           ,SRUQEX,SRUQSC,SRUQCB
     B           ,TRUQEX,TRUQSC,TRUQCB
     C           ,REFU22,Q55U22
     D      ,SRHQEX(1,1,NA),SRHQSC(1,1,NA),SRHQCB(1,1,NA)
     E      ,TRHQAB(1,1,NA)
     F      ,RHINFO(1,1,NA))
      ENDIF
  140 CONTINUE

  150 CONTINUE
      IF(NTRACE.LE.0) RETURN

C**** Optional Tracer aerosols initializations
      DO NT=1,NTRACE
      NA=ITR(NT)
      AREFF=TRRDRY(NT)
      CALL GETMIE(NA,AREFF,SRTQEX(1,1,NT),SRTQSC(1,1,NT),SRTQCB(1,1,NT)
     +                    ,TRTQAB(1,1,NT),Q55)
      RTINFO(1,1,NT)=0.0
      RTINFO(1,2,NT)=1.0
      RTINFO(1,3,NT)=AREFF
      RTINFO(1,4,NT)=0.0
      RTINFO(1,5,NT)=1.33333333D0*AREFF*DENAER(NT)/Q55
      RTINFO(1,6,NT)=1.33333333D0*AREFF*DENAER(NT)/Q55
      RTINFO(1,7,NT)=1.0
      RTINFO(1,8,NT)=DENAER(NT)
      RTINFO(1,9,NT)=Q55
      END DO
            !      Define default dry aerosol coefficients for N=2,190
      DO N=2,190
      DO NT=1,NTRACE
        DO K=1,6
        SRTQEX(K,N,NT)=SRTQEX(K,1,NT)
        SRTQSC(K,N,NT)=SRTQSC(K,1,NT)
        SRTQCB(K,N,NT)=SRTQCB(K,1,NT)
        END DO
        DO K=1,33
        TRTQAB(K,N,NT)=TRTQAB(K,1,NT)
        END DO
        DO M=1,9
        RTINFO(N,M,NT)=RTINFO(N,1,NT)
        END DO
      END DO
      END DO

      KDREAD=71           !  Over-write dry coefficients if KRHTRA(NT)=1
      DO NT=1,NTRACE
      NA=ITR(NT)
      IF(KRHTRA(NT).GT.0) THEN
      CALL SETREL(REFDRY(NT),NA,KDREAD
     A           ,SRUQEX,SRUQSC,SRUQCB
     B           ,TRUQEX,TRUQSC,TRUQCB
     C           ,REFU22,Q55U22
     D      ,SRTQEX(1,1,NT),SRTQSC(1,1,NT),SRTQCB(1,1,NT)
     E      ,TRTQAB(1,1,NT)
     F      ,RTINFO(1,1,NT))
      ENDIF
      END DO

      RETURN


C-----------------
      ENTRY GETAER
C-----------------

      NRHNAN(:,:) = 1
      DO 230 L=1,NL
      IF(RHL(L).GT.0.9005D0) THEN
      XRH=(RHL(L)-0.899499D0)*1000.D0
      NRH=XRH+64
      IF(NRH.GT.164) NRH=164
      ELSE
      XRH=RHL(L)*100.D0+0.5D0
      NRH=XRH-26
      IF(NRH.LT.0) NRH=0
      ENDIF
      DO 220 NA=1,4
      IF(KDELIQ(L,NA).EQ.0) THEN
      IF(RHL(L).GT.RHD(NA)) KDELIQ(L,NA)=1
      ELSE
      IF(RHL(L).LT.RHC(NA)) KDELIQ(L,NA)=0
      ENDIF
      NRHNAN(L,NA)=NRH*KDELIQ(L,NA)+1
  220 CONTINUE
  230 CONTINUE

      IF(MADAER.LE.0) GO TO 500

      DO NA=1,6
      CALL REPART(A6JDAY(1,NA,ILON,JLAT),PLBA09,10,ATAULX(1,NA),PLB,NLP)
      END DO

      FSXTAU=FSTAER*FSAAER+1.D-10
      FTXTAU=FTTAER*FTAAER
                           !            (Solar BCI,BCB components)
      DO 250 L=1,NL
      DO 240 K=1,6
      SRAEXT(L,K)=SRBQEX(K,5)*ATAULX(L,5)*FSXTAU*FS8OPX(5)
     +           +SRBQEX(K,6)*ATAULX(L,6)*FSXTAU*FS8OPX(6)
      SRASCT(L,K)=SRBQSC(K,5)*ATAULX(L,5)*FSXTAU*FS8OPX(5)
     +           +SRBQSC(K,6)*ATAULX(L,6)*FSXTAU*FS8OPX(6)
      SRAGQL     =SRBQSC(K,5)*ATAULX(L,5)*SRBQCB(K,5)
     +           +SRBQSC(K,6)*ATAULX(L,6)*SRBQCB(K,6)
      SRAGCB(L,K)=SRAGQL/(SRASCT(L,K)+1.D-10)
  240 CONTINUE
  250 CONTINUE
                          !           (Thermal BCI,BCB components)
      DO 270 L=1,NL
      DO 260 K=1,33
      TRAALK(L,K)=TRBQAB(K,5)*ATAULX(L,5)*FTXTAU*FT8OPX(5)
     +           +TRBQAB(K,6)*ATAULX(L,6)*FTXTAU*FT8OPX(6)
      IF(PLB(L).GT.10.D0) GO TO 260
      TRAALK(L,K)=0.D0
  260 CONTINUE
  270 CONTINUE

      DO 330 L=1,NL
      DO 320 K=1,6
      DO 310 NA=1,4
      RHFTAU=RHINFO(NRHNAN(L,NA),2,NA)*ATAULX(L,NA)*FSXTAU*FS8OPX(NA)
      SRAEXT(L,K)=SRAEXT(L,K)
     +        +SRHQEX(K,NRHNAN(L,NA),NA)*RHFTAU
      SRASCT(L,K)=SRASCT(L,K)
     +        +SRHQSC(K,NRHNAN(L,NA),NA)*RHFTAU
      SRAGQL     =SRAGCB(L,K)*SRASCT(L,K)+SRHQCB(K,NRHNAN(L,NA),NA)
     +           *SRHQSC(K,NRHNAN(L,NA),NA)*RHFTAU
      SRAGCB(L,K)=SRAGQL/(SRASCT(L,K)+1.D-10)
  310 CONTINUE
  320 CONTINUE
  330 CONTINUE

      DO 360 L=1,NL
      DO 350 K=1,33
      DO 340 NA=1,4
      RHFTAU=RHINFO(NRHNAN(L,NA),2,NA)*ATAULX(L,NA)*FTXTAU*FT8OPX(NA)
      TRAALK(L,K)=TRAALK(L,K)
     +        +TRHQAB(K,NRHNAN(L,NA),NA)*RHFTAU
  340 CONTINUE
  350 CONTINUE
  360 CONTINUE

  500 CONTINUE
      IF(NTRACE.LE.0) RETURN

C     ------------------------------------------------------------------
C     Option to add on Tracer Type aerosol thermal & solar contributions
C
C     NOTE:  Aerosol carried as a tracer is assumed to be in kg/m2 units
C     ------------------------------------------------------------------

      FSXTAU=FSTAER*FSBAER+1.D-10
      FTXTAU=FTTAER*FTBAER

      DO L=1,NL
      DO K=1,6
      DO NT=1,NTRACE
      NA=ITR(NT)
      RHFTAU=RTINFO(NRHNAN(L,NA),2,NT)*TTAULX(L,NT)*FSXTAU*FSTOPX(NT)
      SRBEXT(L,K)=SRBEXT(L,K)+SRTQEX(K,NRHNAN(L,NA),NT)*RHFTAU
      SRBSCT(L,K)=SRBSCT(L,K)+SRTQSC(K,NRHNAN(L,NA),NT)*RHFTAU
      SRBGQL     =SRBGCB(L,K)*SRBSCT(L,K)+SRTQCB(K,NRHNAN(L,NA),NT)
     +           *SRTQSC(K,NRHNAN(L,NA),NT)*RHFTAU
      SRBGCB(L,K)=SRBGQL/(SRBSCT(L,K)+1.D-10)
      END DO
      END DO
      END DO

      DO L=1,NL
      DO K=1,33
      DO NT=1,NTRACE
      NA=ITR(NT)
      RHFTAU=RTINFO(NRHNAN(L,NA),2,NT)*TTAULX(L,NT)*FTXTAU*FTTOPX(NT)
      TRBALK(L,K)=TRBALK(L,K)+TRTQAB(K,NRHNAN(L,NA),NT)*RHFTAU
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE SETAER


      SUBROUTINE UPDAER(JYEARA,JJDAYA)

cc    INCLUDE 'rad00def.radCOMMON.f'
C     ------------------------------------------------------------------
C     Reads: sep2003_XXX_Koch_kg_m2_72x46x9_1850-1990 aerosol kg/m2 data
C     for SUI,OCI,BCI, and PRE (PRE=SNP,SBP,SSP,ANP,ONP,OBP,ANI,OCB,BCB)
C
C     Makes: A6YEAR(72,46,9,0:12,6), A6JDAY(9,6,72,46) (dry aerosol Tau)
C     ------------------------------------------------------------------

      USE FILEMANAGER, only : openunit,closeunit
      implicit none

      integer, intent(in) :: jyeara,jjdaya
      REAL*4 A6YEAR(72,46,9,0:12,6)
      REAL*4  PREDD(72,46,9,12,10),SUIDD(72,46,9,12,8)
      REAL*4  OCIDD(72,46,9,12, 8),BCIDD(72,46,9,12,8),AMON(72,46,9)
      save A6YEAR,PREDD,SUIDD,OCIDD,BCIDD

      CHARACTER*80 XTITLE
      CHARACTER*40, dimension(4) :: RDFILE = (/            !  Input data
     1            'sep2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850'
     2           ,'sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990'
     3           ,'sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990'
     4           ,'sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990'/)

      CHARACTER*40, dimension(4) :: RDFGEN = (/          ! generic names
     * 'TAero_PRE','TAero_SUI','TAero_OCI','TAero_BCI'/)

C                TROPOSPHERIC AEROSOL COMPOSITIONAL/TYPE PARAMETERS
C                   SO4    SEA    ANT    OCX    BCI    BCB   *BCB  *BCB
C     DATA REFDRY/0.200, 1.000, 0.300, 0.300, 0.100, 0.100, 0.200,0.050/
C
C     DATA REFWET/0.272, 1.808, 0.398, 0.318, 0.100, 0.100, 0.200,0.050/
C
C     DATA DRYM2G/4.667, 0.866, 4.448, 5.018, 9.000, 9.000, 5.521,8.169/
C
CKoch DATA DRYM2G/5.000, 2.866, 8.000, 8.000, 9.000, 9.000, 5.521,8.169/
C
C     DATA RHTMAG/1.788, 3.310, 1.756, 1.163, 1.000, 1.000, 1.000,1.000/
C
CRH70 DATA WETM2G/8.345, 2.866, 7.811, 5.836, 9.000, 9.000, 5.521,8.169/
C
C     DATA Q55DRY/2.191, 2.499, 3.069, 3.010, 1.560, 1.560, 1.914,0.708/
C
C     DATA DENAER/1.760, 2.165, 1.725, 1.500, 1.300, 1.300, 1.300,1.300/
C
C     ------------------------------------------------------------------
C          DRYM2G(I) = 0.75/DENAER(I)*Q55DRY(I)/REFDRY(I)
C          WETM2G(I) = DRYM2G(I)*RHTMAG(I)
C          RHTMAG(I) = Rel Humidity TAU Magnification factor  at RH=0.70
C          REFWET(I) = Rel Humidity REFDRY Magnification      at RH=0.70
C     ------------------------------------------------------------------

C     TROP AEROSOL 1850 BACKGROUND, INDUSTRIAL & BIO-BURNING PARAMETERS
C     DATA AERMIX/
C       Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C       ---------------------------------  ------------------  --------
C        1    2    3    4    5    6    7    8    9   10   11   12   13
C       SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
C    +  1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9/

C      A6YEAR          PRE                  SUI         OCI        BCI
C     ------------------------------------------------------------------
C     NAER=1=SO4 = SNP*1+SBP*2         +  SUI*I,J
C          2=SEA = SSP*3
C          3=ANT = ANP*4+ANI*0,8
C          4=OCX = ONP*5+OBP*6+OCB*0,9             +  OCI*I,J
C          5=BCI =                                               BCI*I,J
C          6=BCB = BBP*7,BCB*0,10
C     ------------------------------------------------------------------
C     Aerosol input data is from designated source files PRE SUI OCI BCI
C     Aerosol output is accumulated for 6-A6YEAR designated compositions
C           SNP*1 represents AERMIX(1)*PRE(I,J,L,M,1) = 1850 Natural SO4
C           SBP*2 represents AERMIX(2)*PRE(I,J,L,M,2) = 1850 BioBurn SO4
C         SUI*I,J represents AERMIX(8)*(SUI(I,J,L,M,I)intSUI(I,J,L,M,J))
C           SSP*1 represents AERMIX(3)*PRE(I,J,L,M,3) = 1850 SeaSalt
C         BCB*0,10represents AERMIX(11)*(0_interpol_PRE(I,J,L,M,10)) BCB
C         (which is interpolated linearly in time from 0 amount in 1850)
C
C         1850 Background   Sulfate  SO4 = 0.870 Natural + 0.130 BioBurn
C         1850 Background  Sea Salt  SEA =  all Natural-Mean SSP SeaSalt
C         1850 Background AmmNitrate ANT = ANP=(1.26/5.27)*1990(ANP+ANI)
C         1850 Background Org Carbon OCX = 0.162 Natural + 0.838 BioBurn
C         1850 Background Blk Carbon BCI = 0  (No Industrial BC in 1850)
C         1850 Background Blk Carbon BCB = BBP,all of 1850 BC is BioBurn
C     ------------------------------------------------------------------
      logical, parameter :: qbinary=.true.   ; logical qexist
      integer, save :: IFILE=11, IFIRST=1, JYRNOW=0

      integer ia,idd,ndd,m,mi,mj,i,j,l,n,jyearx,iys,jys,iyc,jyc
      real*8 WTANI,WTOCB,WTBCB,wt75,swti,swtj,cwti,cwtj,xmi,wtmi,wtmj
      IF(IFIRST.EQ.1) THEN
C                                       READ Input PRE,SUI,OCI,BCI Files
C                                       --------------------------------
      inquire (file=RDFGEN(1),exist=qexist) ! decide whether specific or
      if(qexist) RDFILE=RDFGEN              !     generic names are used
      inquire (file=RDFILE(1),exist=qexist) !     stop if neither exist
      if(.not.qexist) call stop_model('updaer: no TropAero files',255)

      DO 106 IA=1,4
      call openunit (RDFILE(IA),ifile,qbinary)   ! unform.:qbinary=true
      NDD=8
      IF(IA.EQ.1) NDD=10
      DO 105 IDD=1,NDD
      DO 104 M=1,12
      READ (IFILE) XTITLE,AMON
      DO 103 L=1,9
      DO 102 J=1,46
      DO 101 I=1,72
      IF(IA.EQ.1) PREDD(I,J,L,M,IDD)=AMON(I,J,L)
      IF(IA.EQ.2) SUIDD(I,J,L,M,IDD)=AMON(I,J,L)
      IF(IA.EQ.3) OCIDD(I,J,L,M,IDD)=AMON(I,J,L)
      IF(IA.EQ.4) BCIDD(I,J,L,M,IDD)=AMON(I,J,L)
  101 CONTINUE
  102 CONTINUE
  103 CONTINUE
  104 CONTINUE
  105 CONTINUE
      call closeunit (ifile)
  106 CONTINUE
      IFIRST=0
      ENDIF


C     To time input data READs, JYEARX is set ahead of JYEARA by 15 days
C     ------------------------------------------------------------------
      JYEARX=MIN(JYEARA+(JJDAYA+15)/366,2050)

      IF(JYEARX.EQ.JYRNOW) GO TO 500    ! Get A6JDAY from current A6YEAR

C     Begin current A6YEAR  with 1850 Background SO4,SEA,ANT,OCX,BCI,BCB
      DO 114 M=1,12
      DO 113 L=1,9
      DO 112 J=1,46
      DO 111 I=1,72
      A6YEAR(I,J,L,M,1)=AERMIX(1)*PREDD(I,J,L,M,1)*1000.D0*DRYM2G(1)
     +                 +AERMIX(2)*PREDD(I,J,L,M,2)*1000.D0*DRYM2G(1)
      A6YEAR(I,J,L,M,2)=AERMIX(3)*PREDD(I,J,L,M,3)*1000.D0*DRYM2G(2)
      A6YEAR(I,J,L,M,3)=AERMIX(4)*PREDD(I,J,L,M,4)*1000.D0*DRYM2G(3)
      A6YEAR(I,J,L,M,4)=AERMIX(5)*PREDD(I,J,L,M,5)*1000.D0*DRYM2G(4)
     +                 +AERMIX(6)*PREDD(I,J,L,M,6)*1000.D0*DRYM2G(4)
      A6YEAR(I,J,L,M,5)=0.D0
      A6YEAR(I,J,L,M,6)=AERMIX(7)*PREDD(I,J,L,M,7)*1000.D0*DRYM2G(6)
  111 CONTINUE
  112 CONTINUE
  113 CONTINUE
  114 CONTINUE
      DO 117 L=1,9                !     Define 1849 Background  Dec data
      DO 116 J=1,46
      DO 115 I=1,72
      A6YEAR(I,J,L,0,1)=A6YEAR(I,J,L,12,1)
      A6YEAR(I,J,L,0,2)=A6YEAR(I,J,L,12,2)
      A6YEAR(I,J,L,0,3)=A6YEAR(I,J,L,12,3)
      A6YEAR(I,J,L,0,4)=A6YEAR(I,J,L,12,4)
      A6YEAR(I,J,L,0,5)=A6YEAR(I,J,L,12,5)
      A6YEAR(I,J,L,0,6)=A6YEAR(I,J,L,12,6)
  115 CONTINUE
  116 CONTINUE
  117 CONTINUE

      IF(JYEARA.GT.1850) THEN                           !   (JYEAR>1850)
      WTANI=GLOPOP(JYEARA)
      WTOCB=(JYEARA-1850)/140.D0
      WTBCB=(JYEARA-1850)/140.D0
      DO 124 M=1,12              !  Add time dependent JYEAR ANI,OCB,BCB
      DO 123 L=1,9
      DO 122 J=1,46
      DO 121 I=1,72
      A6YEAR(I,J,L,M,3)=A6YEAR(I,J,L,M,3)+
     +             AERMIX( 9)*WTANI*PREDD(I,J,L,M, 8)*1000.D0*DRYM2G(3)
      A6YEAR(I,J,L,M,4)=A6YEAR(I,J,L,M,4)+
     +             AERMIX(12)*WTOCB*PREDD(I,J,L,M, 9)*1000.D0*DRYM2G(4)
      A6YEAR(I,J,L,M,6)=A6YEAR(I,J,L,M,6)+
     +             AERMIX(13)*WTBCB*PREDD(I,J,L,M,10)*1000.D0*DRYM2G(6)
  121 CONTINUE
  122 CONTINUE
  123 CONTINUE
  124 CONTINUE
      WTANI=GLOPOP(JYEARA-1)
      WTOCB=(JYEARA-1851)/140.D0
      WTBCB=(JYEARA-1851)/140.D0
      M=12
      DO 127 L=1,9    !  Add time dependent JYEAR-1 ANI,OCB,BCB Dec data
      DO 126 J=1,46
      DO 125 I=1,72
      A6YEAR(I,J,L,0,3)=A6YEAR(I,J,L,0,3)+
     +             AERMIX( 9)*WTANI*PREDD(I,J,L,M, 8)*1000.D0*DRYM2G(3)
      A6YEAR(I,J,L,0,4)=A6YEAR(I,J,L,0,4)+
     +             AERMIX(12)*WTOCB*PREDD(I,J,L,M, 9)*1000.D0*DRYM2G(4)
      A6YEAR(I,J,L,0,6)=A6YEAR(I,J,L,0,6)+
     +             AERMIX(13)*WTBCB*PREDD(I,J,L,M,10)*1000.D0*DRYM2G(6)
  125 CONTINUE
  126 CONTINUE
  127 CONTINUE
      ENDIF

      IF(JYEARA.GT.1850.AND.JYEARA.LT.1876) THEN   !   (1850<JYEAR<1876)
      WT75=(JYEARA-1850)/25.D0
      DO 134 M=1,12            !    Add time dependent JYEAR SUI,OCI,BCI
      DO 133 L=1,9
      DO 132 J=1,46
      DO 131 I=1,72
      A6YEAR(I,J,L,M,1)=A6YEAR(I,J,L,M,1)+
     +              WT75*SUIDD(I,J,L,M,1)*AERMIX( 8)*1000.D0*DRYM2G(1)
      A6YEAR(I,J,L,M,4)=A6YEAR(I,J,L,M,4)+
     +              WT75*OCIDD(I,J,L,M,1)*AERMIX(10)*1000.D0*DRYM2G(4)
      A6YEAR(I,J,L,M,5)=A6YEAR(I,J,L,M,5)+
     +              WT75*BCIDD(I,J,L,M,1)*AERMIX(11)*1000.D0*DRYM2G(5)

  131 CONTINUE
  132 CONTINUE
  133 CONTINUE
  134 CONTINUE
      WT75=(JYEARA-1851)/25.D0
      M=12            !  Add time dependent JYEAR-1 SUI,OCI,BCI Dec data
      DO 137 L=1,9
      DO 136 J=1,46
      DO 135 I=1,72
      A6YEAR(I,J,L,0,1)=A6YEAR(I,J,L,0,1)+
     +              WT75*SUIDD(I,J,L,M,1)*AERMIX( 8)*1000.D0*DRYM2G(1)
      A6YEAR(I,J,L,0,4)=A6YEAR(I,J,L,0,4)+
     +              WT75*OCIDD(I,J,L,M,1)*AERMIX(10)*1000.D0*DRYM2G(4)
      A6YEAR(I,J,L,0,5)=A6YEAR(I,J,L,0,5)+
     +              WT75*BCIDD(I,J,L,M,1)*AERMIX(11)*1000.D0*DRYM2G(5)
  135 CONTINUE
  136 CONTINUE
  137 CONTINUE
      ENDIF


      IF(JYEARA.GT.1875) THEN                         !     (JYEAR>1875)
      CALL STREND(JYEARA,IYS,JYS,SWTI,SWTJ)
      CALL CTREND(JYEARA,IYC,JYC,CWTI,CWTJ)
      DO 144 M=1,12            !    Add time dependent JYEAR SUI,OCI,BCI
      DO 143 L=1,9
      DO 142 J=1,46
      DO 141 I=1,72
      A6YEAR(I,J,L,M,1)=A6YEAR(I,J,L,M,1)+AERMIX( 8)*1000.D0*DRYM2G(1)*
     +                 (SWTI*SUIDD(I,J,L,M,IYS)+SWTJ*SUIDD(I,J,L,M,JYS))
      A6YEAR(I,J,L,M,4)=A6YEAR(I,J,L,M,4)+AERMIX(10)*1000.D0*DRYM2G(4)*
     +                 (CWTI*OCIDD(I,J,L,M,IYC)+CWTJ*OCIDD(I,J,L,M,JYC))
      A6YEAR(I,J,L,M,5)=A6YEAR(I,J,L,M,5)+AERMIX(11)*1000.D0*DRYM2G(5)*
     +                 (CWTI*BCIDD(I,J,L,M,IYC)+CWTJ*BCIDD(I,J,L,M,JYC))
  141 CONTINUE
  142 CONTINUE
  143 CONTINUE
  144 CONTINUE
      CALL STREND(JYEARA-1,IYS,JYS,SWTI,SWTJ)
      CALL CTREND(JYEARA-1,IYC,JYC,CWTI,CWTJ)
      M=12            !  Add time dependent JYEAR-1 SUI,OCI,BCI Dec data
      DO 147 L=1,9
      DO 146 J=1,46
      DO 145 I=1,72
      A6YEAR(I,J,L,0,1)=A6YEAR(I,J,L,0,1)+AERMIX( 8)*1000.D0*DRYM2G(1)*
     +                 (SWTI*SUIDD(I,J,L,M,IYS)+SWTJ*SUIDD(I,J,L,M,JYS))
      A6YEAR(I,J,L,0,4)=A6YEAR(I,J,L,0,4)+AERMIX(10)*1000.D0*DRYM2G(4)*
     +                 (CWTI*OCIDD(I,J,L,M,IYC)+CWTJ*OCIDD(I,J,L,M,JYC))
      A6YEAR(I,J,L,0,5)=A6YEAR(I,J,L,0,5)+AERMIX(11)*1000.D0*DRYM2G(5)*
     +                 (CWTI*BCIDD(I,J,L,M,IYC)+CWTJ*BCIDD(I,J,L,M,JYC))
  145 CONTINUE
  146 CONTINUE
  147 CONTINUE
      ENDIF
      JYRNOW=JYEARX


C      A6JDAY is interpolated daily from A6YEAR seasonal data via JJDAYA
C      -----------------------------------------------------------------

  500 CONTINUE
      XMI=(JJDAYA+JJDAYA+31-(JJDAYA+15)/61+(JJDAYA+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYA
      WTMI=1.D0-WTMJ
      IF(MI.GT.11) MI=0
      MJ=MI+1
      DO 540 J=1,46
      DO 530 I=1,72
      DO 520 N=1,6
      DO 510 L=1,9
      A6JDAY(L,N,I,J)=WTMI*A6YEAR(I,J,L,MI,N)+WTMJ*A6YEAR(I,J,L,MJ,N)
  510 CONTINUE
  520 CONTINUE
  530 CONTINUE
  540 CONTINUE
      RETURN        !  A6JDAY(9,6,72,46) is used in GETAER via ILON,JLAT
      END SUBROUTINE UPDAER


      REAL*8 FUNCTION GLOPOP(JYEAR)
      IMPLICIT none

C     ----------------------------------------------------------------
C     GLOPOP = normalized global population trend set to unity in 1990
C              based on UN statistics & population projections to 2050
C
C     GLOPOP = 0.000 for 1850 and earlier
C            = 1.000 for 1990
C            = 1.658 for 2050 and later
C     ----------------------------------------------------------------

      integer, intent(in) :: jyear
      real*8 :: GPNORM = 5.27-1.26 ,DNGPOP(21), GPOP(21) = (/
C               1850                     1900                     1950
     A          1.26,1.33,1.41,1.49,1.57,1.65,1.75,1.86,2.07,2.30,2.52
C                                        2000                     2050
     B              ,3.02,3.70,4.44,5.27,6.06,6.79,7.50,8.11,8.58,8.91/)
      integer i,iy
      real*8 xy,dy

      DO 110 I=1,21
      DNGPOP(I)=(GPOP(I)-GPOP(1))/GPNORM
  110 CONTINUE
      XY=(JYEAR-1840)/10.D0
      IY=XY
      DY=XY-IY
      IF(IY.LT.1) THEN
      IY=1
      DY=0.D0
      ENDIF
      IF(IY.GT.20) THEN
      IY=20
      DY=1.D0
      ENDIF
      GLOPOP=DNGPOP(IY)+DY*(DNGPOP(IY+1)-DNGPOP(IY))
      RETURN
      END FUNCTION GLOPOP

      SUBROUTINE SETDST
      IMPLICIT NONE

C     ---------------------------------------------------------------
C     MONTHLY-MEAN DESERT DUST CLIMATOLOGY
C     ---------------------------------------------------------------

C           OUTPUT: via SRDEXT(L,K)   D Dust Extinction Optical Depth
C                       SRDSCT(L,K)   D Dust Scattering Optical Depth
C                       SRDGCB(L,K)   D Dust Asymmetry Parameter  g
C                       TRDALK(L,K)  Thermal Absorption Optical Depth

C     Tau Scaling Factors:    Solar    Thermal    apply to:
C                             FSTAER   FTTAER  ! Total Aerosol
C                             FSDAER   FTDAER  ! Dust  Aerosol
C
C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW   (All-type) Aerosol Optical Depth
C                        FTTAER    LW   (All-type) Aerosol Optical Depth
C                        FSDAER    SW   Dust Aer   Aerosol Optical Depth
C                        FTDAER    LW   Dust Aer   Aerosol Optical Depth
C                        -----------------------------------------------

C     Select Desert Dust (NA=7) Mie scattering parameters for REDUST(N)
!nu   REAL*8 TAUCON(8),pidust(8)
      INTEGER, INTENT(IN) :: JYEARD,JJDAYD
      REAL*8 XMI,WTMI,WTMJ,SRDGQL,FSXTAU,FTXTAU,DTAULX(LX,8)
      INTEGER I,J,K,L,N,MI,MJ

      DO 110 N=1,8
      CALL GETMIE(7,REDUST(N),QXDUST(1,N),QSDUST(1,N),QCDUST(1,N)
     +                    ,ATDUST(1,N),QDST55(N))
!nu   TAUCON(N)=0.75E+03*QDST55(N)/(RODUST(N)*REDUST(N))
!nu   PIDUST(N)=QSDUST(6,N)/(QXDUST(6,N)+1.D-10)
  110 CONTINUE

      RETURN


C--------------------------------
      ENTRY UPDDST(JYEARD,JJDAYD)
C--------------------------------

C     ------------------------------------------------------------------
C     Makes DDJDAY(9,8,72,46) from TDUST(72,46,9,8,13) read in in RCOMP1
C
C      DDJDAY is interpolated daily from  TDUST seasonal data via JJDAYD
C      -----------------------------------------------------------------
!nu   JYEARX=MIN(JYEARD,(JJDAYD+15)/366,2050)

  500 CONTINUE
      XMI=(JJDAYD+JJDAYD+31-(JJDAYD+15)/61+(JJDAYD+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYD
      WTMI=1.D0-WTMJ
      IF(MI.LT.1) MI=12
      IF(MI.GT.12) MI=1
      MJ=MI+1
      IF(MJ.GT.12) MJ=1
      DO 540 J=1,46
      DO 530 I=1,72
      DO 520 N=1,8
      DO 510 L=1,9
      DDJDAY(L,N,I,J)=WTMI*TDUST(I,J,L,N,MI)+WTMJ*TDUST(I,J,L,N,MJ)
  510 CONTINUE
  520 CONTINUE
  530 CONTINUE
  540 CONTINUE
      RETURN        !  DDJDAY(9,8,72,46) is used in GETDST via ILON,JLAT


C-----------------
      ENTRY GETDST
C-----------------

      DO 200 N=1,8
      CALL REPART(DDJDAY(1,N,ILON,JLAT),PLBA09,10,DTAULX(1,N),PLB,NLP)
  200 CONTINUE

C                     Apply Solar/Thermal Optical Depth Scaling Factors
C                              Dust Aerosol  Solar   FSXD=FSTAER*FSDAER
C                              Dust Aerosol Thermal  FTXD=FSTAER*FTDAER
C                              ----------------------------------------

      FSXTAU=FSTAER*FSDAER+1.D-10
      FTXTAU=FTTAER*FTDAER

      DO 220 K=1,6
      DO 210 L=1,NL
      SRDEXT(L,K)=2.D-10
      SRDSCT(L,K)=1.D-10
      SRDGCB(L,K)=0.D0
  210 CONTINUE
  220 CONTINUE

      DO 240 K=1,33
      DO 230 L=1,NL
      TRDALK(L,K)=0.D0
  230 CONTINUE
  240 CONTINUE

      DO 270 L=1,NL
      DO 260 K=1,6
      DO 250 N=1,8
      SRDEXT(L,K)=SRDEXT(L,K)+QDUST(L,K)*DTAULX(L,N)*FSXTAU*FS8OPX(7)
      SRDSCT(L,K)=SRDSCT(L,K)+SDUST(L,K)*DTAULX(L,N)*FSXTAU*FS8OPX(7)
      SRDGQL     = CDUST(L,K)*SDUST(L,K)*DTAULX(L,N)*FSXTAU*FS8OPX(7)
      SRDGCB(L,K)=SRDGQL/(SRDSCT(L,K)+1.D-10)
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE

      DO 300 L=1,NL
      DO 290 K=1,33
      DO 280 N=1,8
      TRDALK(L,K)=TRDALK(L,K)+ADUST(L,K)*DTAULX(L,N)*FTXTAU*FT8OPX(7)
  280 CONTINUE
  290 CONTINUE
  300 CONTINUE

      RETURN
      END SUBROUTINE SETDST

      SUBROUTINE SETVOL
      IMPLICIT NONE

      REAL*8, SAVE :: E24LAT(25),EJMLAT(47)
      REAL*8  HLATTF(4)
      REAL*8, PARAMETER :: HLATKM(5) = (/15.0, 20.0, 25.0, 30.0, 35.0/)
      INTEGER, SAVE :: LATVOL = 0

      real*8, parameter :: htplim=1.d-3
      REAL*8, SAVE :: FSXTAU,FTXTAU
      INTEGER, SAVE :: NJ25,NJJM
      INTEGER, INTENT(IN) :: JYEARV,JDAYVA
      INTEGER J,L,MI,MJ,K
      REAL*8 XYYEAR,XYI,WMI,WMJ,SUM,SUMHTF,SIZVOL

C     ------------------------------------------------------------------
C     Tau Scaling Factors:    Solar    Thermal    apply to:
c                             FSTAER   FTTAER  ! Total Aerosol
c                             FSVAER   FTVAER  ! SETVOL Aer

C     Control Parameters/Aerosol Scaling (kill) Factors
C                        FSTAER    SW  (All-type) Aerosol Optical Depth
C                        FTTAER    LW  (All-type) Aerosol Optical Depth
C                        FSVAER    SW  SETVOLonly Aerosol Optical Depth
C                        FTVAER    LW  SETVOLonly Aerosol Optical Depth
C                        -----------------------------------------------

C     -----------------------------------------------------------------
C     VEFF0   Selects Size Distribution Variance (this affects Thermal)
C     REFF0   Selects Effective Particle Size for Archive Volcanic Data
C     -----------------------------------------------------------------

      FSXTAU=FSTAER*FSVAER
      FTXTAU=FTTAER*FTVAER

C                   Set Grid-Box Edge Latitudes for Data Repartitioning
C                   ---------------------------------------------------
      NJ25=25
      DO 110 J=2,24
      E24LAT(J)=-90.D0+(J-1.5D0)*180.D0/23.D0
  110 CONTINUE
      E24LAT( 1)=-90.D0
      E24LAT(25)= 90.D0
      NJJM=46+1
      DO 120 J=2,46
      EJMLAT(J)=-90.D0+(J-1.5D0)*180.D0/(MLAT46-1)
  120 CONTINUE
      EJMLAT(   1)=-90.D0
      EJMLAT(NJJM)= 90.D0

      DO 150 L=1,NL
      HTPROF(L)=0.D0
  150 CONTINUE

C                       -----------------------------------------------
C                       Initialize H2SO4 Q,S,C,A Tables for Input VEFF0
C                       -----------------------------------------------
C     ------------------
      CALL SETQVA(VEFF0)
C     ------------------

      RETURN


C--------------------------------
      ENTRY UPDVOL(JYEARV,JDAYVA)
C--------------------------------

C                                          (Makiko's 1850-1999 data)
C                                          -------------------------
      XYYEAR=JYEARV+JDAYVA/366.D0
      IF(XYYEAR.LT.1850.D0) XYYEAR=1850.D0
      XYI=(XYYEAR-1850.D0)*12.D0+1.D0
      IF(XYI.GT.1799.999D0) XYI=1799.999D0
      MI=XYI
      WMJ=XYI-MI
      WMI=1.D0-WMJ
      MJ=MI+1
      DO 250 J=1,24
      GDATA(J)=WMI*V4TAUR(MI,J,5)+WMJ*V4TAUR(MJ,J,5)
  250 CONTINUE
      CALL RETERP(GDATA,E24LAT,NJ25,SIZLAT,EJMLAT,NJJM)
      DO 270 K=1,4
      DO 260 J=1,24
      FDATA(J)=WMI*V4TAUR(MI,J,K)+WMJ*V4TAUR(MJ,J,K)
  260 CONTINUE
      CALL RETERP(FDATA,E24LAT,NJ25,HTFLAT(1,K),EJMLAT,NJJM)
  270 CONTINUE
      DO 290 J=1,24
      SUM=0
      DO 280 K=1,4
      SUM=SUM+HTFLAT(J,K)
  280 CONTINUE
      TAULAT(J)=SUM
  290 CONTINUE

      RETURN


C-----------------
      ENTRY GETVOL
C-----------------

      IF(MRELAY.GT.0)    GO TO 300
      IF(JLAT.EQ.LATVOL) GO TO 350

C                      Set JLAT Dependent Aerosol Distribution and Size
C                      ------------------------------------------------
  300 CONTINUE

      DO 310 K=1,4
      HLATTF(K)=HTFLAT(JLAT,K)
  310 CONTINUE
      CALL REPART(HLATTF,HLATKM,5,HTPROF,HLB0,NLP)
!nu   LHPMAX=0       ! not used
!nu   LHPMIN=NL      ! not used
!nu   DO L=1,NL
!nu   N=NLP-L
!nu   IF(HTPROF(L).GE.HTPLIM) LHPMAX=L
!nu   IF(HTPROF(N).GE.HTPLIM) LHPMIN=N
!nu   END DO
      SUMHTF=1.D-10
      DO 330 L=1,NL
      IF(HTPROF(L).LT.HTPLIM) HTPROF(L)=0.D0
      SUMHTF=SUMHTF+HTPROF(L)
  330 CONTINUE

      SIZVOL=SIZLAT(JLAT)

C                        Select H2SO4 Q,S,C,A Tables for  Size = SIZVOL
C                        ----------------------------------------------

C------------------------
      CALL GETQVA(SIZVOL)
C------------------------

      LATVOL=JLAT
  350 CONTINUE
C                                  ------------------------------------
C                                  H2SO4 Thermal Contribution in TRVALK
C                                  ------------------------------------
      DO 420 K=1,33
      DO 410 L=1,NL
      TRVALK(L,K)=HTPROF(L)*AVH2S(K)*FTXTAU
  410 CONTINUE
  420 CONTINUE

C                      H2SO4 Solar Contribution in SRVEXT,SRVSCT,SRVGCB
C                      ------------------------------------------------

      DO 440 K=1,6
      DO 430 L=1,NL
      SRVEXT(L,K)=QVH2S(K)*HTPROF(L)*FSXTAU
      SRVSCT(L,K)=SVH2S(K)*HTPROF(L)*FSXTAU*PIVMAX
      SRVGCB(L,K)=GVH2S(K)
  430 CONTINUE
  440 CONTINUE

      RETURN
      END SUBROUTINE SETVOL

      SUBROUTINE SETQVA(VEFF)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SETQVA   Selects (interpolates) H2SO4 Mie Parameters for specified
C              Variance VEFF for subsequent Size interpolation by GETQVA
C     ------------------------------------------------------------------

ceq   REAL*8 SRQV( 6,20),SRSV( 6,20),SRGV( 6,20),Q55V(   20),REFV(20)
ceq   REAL*8 TRQV(33,20),TRSV(33,20),TRGV(33,20),TRAV(33,20),VEFV(20)
      REAL*8 TRAB(33,20),Q5(5),RV20(20),QV20(20)
      REAL*8, PARAMETER, DIMENSION(5) ::
     *     V5=(/0.1D0,0.2D0,0.3D0,0.4D0,0.5D0/)
      SAVE TRAB

C     ------------------------------------------------------------------
C     SRVQEX Volcanic Aerosol sizes (Reff) range from 0.1 to 5.0 microns
C     To utilize equal interval interpolation, Reff N=9,20 are redefined
C     so Volcanic Aerosol sizes have effective range of 0.1-2.0 microns.
C     ------------------------------------------------------------------
      REAL*8, INTENT(IN) :: SIZVOL,VEFF
      REAL*8 REFN,RADX,WTJHI,WTJLO
      INTEGER I,K,N,JRXLO,JRXHI

      DO 130 N=1,20
      RV20(N)=REFV20(N,1)
      VEFV(N)=VEFF
      REFV(N)=N/10.D0
      DO 110 I=1,5
      Q5(I)=Q55V20(N,I)
  110 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,Q55V(N),1.D0,1.D0,1)
      DO 115 K=1,6
      DO 111 I=1,5
      Q5(I)=SRVQEX(K,N,I)
  111 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,SRQV(K,N),1.D0,1.D0,1)
      DO 112 I=1,5
      Q5(I)=SRVQSC(K,N,I)
  112 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,SRSV(K,N),1.D0,1.D0,1)
      DO 113 I=1,5
      Q5(I)=SRVQCB(K,N,I)
  113 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,SRGV(K,N),1.D0,1.D0,1)
  115 CONTINUE
      DO 120 K=1,33
      DO 116 I=1,5
      Q5(I)=TRVQEX(K,N,I)
  116 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRQV(K,N),1.D0,1.D0,1)
      DO 117 I=1,5
      Q5(I)=TRVQSC(K,N,I)
  117 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRSV(K,N),1.D0,1.D0,1)
      DO 118 I=1,5
      Q5(I)=TRVQCB(K,N,I)
  118 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRGV(K,N),1.D0,1.D0,1)
      DO 119 I=1,5
      Q5(I)=TRVQAL(K,N,I)
  119 CONTINUE
      CALL SPLINE(V5,Q5,5,VEFF,TRAV(K,N),1.D0,1.D0,1)
  120 CONTINUE
      DO 125 K=1,33
      TRAB(K,N)=TRQV(K,N)-TRSV(K,N)
  125 CONTINUE
  130 CONTINUE
      DO 131 N=1,20
      QV20(N)=Q55V(N)
  131 CONTINUE
      DO 132 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,Q55V(N),1.D0,1.D0,1)
  132 CONTINUE
      DO 140 K=1,6
      DO 133 N=1,20
      QV20(N)=SRQV(K,N)
  133 CONTINUE
      DO 134 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRQV(K,N),1.D0,1.D0,1)
  134 CONTINUE
      DO 135 N=1,20
      QV20(N)=SRSV(K,N)
  135 CONTINUE
      DO 136 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRSV(K,N),1.D0,1.D0,1)
  136 CONTINUE
      DO 137 N=1,20
      QV20(N)=SRGV(K,N)
  137 CONTINUE
      DO 138 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,SRGV(K,N),1.D0,1.D0,1)
  138 CONTINUE
  140 CONTINUE
      DO 150 K=1,33
      DO 141 N=1,20
      QV20(N)=TRQV(K,N)
  141 CONTINUE
      DO 142 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRQV(K,N),1.D0,1.D0,1)
  142 CONTINUE
      DO 143 N=1,20
      QV20(N)=TRSV(K,N)
  143 CONTINUE
      DO 144 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRSV(K,N),1.D0,1.D0,1)
  144 CONTINUE
      DO 145 N=1,20
      QV20(N)=TRGV(K,N)
  145 CONTINUE
      DO 146 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRGV(K,N),1.D0,1.D0,1)
  146 CONTINUE
      DO 147 N=1,20
      QV20(N)=TRAV(K,N)
  147 CONTINUE
      DO 148 N=9,20
      REFN=REFV(N)
      CALL SPLINE(RV20,QV20,20,REFN,TRAV(K,N),1.D0,1.D0,1)
  148 CONTINUE
      DO 149 N=9,20
      TRAB(K,N)=TRQV(K,N)-TRSV(K,N)
  149 CONTINUE
  150 CONTINUE

      RETURN

C-------------------------
      ENTRY GETQVA(SIZVOL)
C-------------------------

C     ------------------------------------------------------------------
C     Volcanic Aerosol sizes have effective range of  0.1 - 2.0 microns.
C     ------------------------------------------------------------------

      RADX=SIZVOL*10.D0
      IF(RADX.LT.1.000001D0) RADX=1.000001D0
      IF(RADX.GT.19.99999D0) RADX=19.99999D0
      JRXLO=RADX
      WTJHI=RADX-JRXLO
      WTJLO=1.D0-WTJHI
      JRXHI=JRXLO+1

      DO 210 I=1,6
      QVH2S(I)=WTJLO*SRQV(I,JRXLO)+WTJHI*SRQV(I,JRXHI)
      SVH2S(I)=WTJLO*SRSV(I,JRXLO)+WTJHI*SRSV(I,JRXHI)
      GVH2S(I)=WTJLO*SRGV(I,JRXLO)+WTJHI*SRGV(I,JRXHI)
  210 CONTINUE

      Q55H2S=WTJLO*Q55V(JRXLO)+WTJHI*Q55V(JRXHI)

      DO 220 I=1,33
      AVH2S(I)=WTJLO*TRAB(I,JRXLO)+WTJHI*TRAB(I,JRXHI)
  220 CONTINUE

      RETURN
      END SUBROUTINE SETQVA

      SUBROUTINE SETCLD
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Control Parameters used in SETCLD,GETCLD,GETEPS: defined in RADPAR
C
C             ICE012  Selects Water, Non-Mie, Mie  Ice Cloud Qex,Qsc,Pi0
C             TAUWC0  Minimum Optical Depth for Water Clouds
C             TAUIC0  Minimum Optical Depth for   Ice Clouds
C             FCLDTR  Scaling Factor for Thermal Cloud Optical Depth
C             FCLDSR  Scaling Factor for  Solar  Cloud Optical Depth
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C
C-----------------------------------------------------------------------
C                    Define Cloud Absorption Cross-Sections
C
C     Selected by:   ICE012 = 0    Liquid Water Droplets   (N =  1 -  5)
C                    ICE012 = 1    Ice - Non-Spherical     (N =  6 - 10)
C                    ICE012 = 2    Ice - Mie (Spherical)   (N = 11 - 15)
C
C     Define Solar,Thermal Cloud Single Scattering Albedo: SRCQPI( 6,15)
C                                                          TRCQPI(33,15)
C-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: JYEARE,JJDAYE
      REAL*8 SIZWCL,SIZICL,XRW,XMW,XPW,EPS,VEP,VEP1,VEP2,VEPP,TAUWCL
     *     ,TAUICL,QAWATK,QPWATK,SRCGFW,QXWATK,QSWATK,QGWATK,XRI,XMI,XPI
     *     ,QAICEK,QPICEK,SRCGFC,QXICEK,QSICEK,QGICEK,SCTTAU,GCBICE
     *     ,SCTGCB,TCTAUW,TCTAUC,ALWATK,WTI,WTW,ALICEK,TRCTCI,XJDAY,XMO
     *     ,WTMJ,WTMI
      INTEGER I,J,N,K,L,LBOTCW,LTOPCW,LBOTCI,LTOPCI,IRWAT,IRICE,MI,MJ

      DO 120 N=1,15
      DO 110 K=1,33
      TRCQAB(K,N)=TRCQEX(K,N)-TRCQSC(K,N)
      TRCQPI(K,N)=TRCQSC(K,N)/TRCQEX(K,N)
  110 CONTINUE
  120 CONTINUE

      DO 140 N=1,15
      DO 130 K=1,6
      SRCQPI(K,N)=SRCQSC(K,N)/SRCQEX(K,N)
  130 CONTINUE
  140 CONTINUE

C                          Initialize  GETCLD Output Parameters to Zero
C                          --------------------------------------------
      DO 150 K=1,33
      TRCTCA(K)=0.D0
  150 CONTINUE
      DO 180 L=1,NL
      DO 160 K=1,33
      TRCALK(L,K)=0.D0
  160 CONTINUE
      DO 170 K=1,6
      SRCEXT(L,K)=1.D-20
      SRCSCT(L,K)=0.D0
      SRCGCB(L,K)=0.D0
  170 CONTINUE
  180 CONTINUE

      RETURN
C-----------------
      ENTRY GETCLD
C-----------------

C-----------------------------------------------------------------------
C           Define:    TRCALK(LX,33)  Thermal Radiation Cloud Absorption
C                      TRCTCA(33)     Thermal Radiation Top Cloud Albedo
C
C                      SRCEXT(LX,6)   Solar Radiation Cloud Ext Op Depth
C                      SRCSCT(LX,6)   Solar Radiation Cloud Sct Op Depth
C                      SRCGCB(LX,6)   Solar Radiation Cloud Asym Param g
C
C                         LTOPCL      Top Cloud Layer Location
C                         LBOTCL      Bot Cloud Layer Location
C
C                         LTOPCW      Top Water Cloud Layer Location
C                         LBOTCW      Bot Water Cloud Layer Location
C
C                         LTOPCI      Top Ice Cloud Layer Location
C                         LBOTCI      Bot Ice Cloud Layer Location
C
C-----------------------------------------------------------------------

      LBOTCW=0
      LTOPCW=0
      LBOTCI=0
      LTOPCI=0
      DO 210 K=1,33
      TRCTCA(K)=0.D0
  210 CONTINUE
      DO 280 L=1,NL
      DO 220 K=1,33
      TRCALK(L,K)=0.D0
  220 CONTINUE
      DO 230 K=1,6
      SRCEXT(L,K)=1.D-20
      SRCSCT(L,K)=1.D-30
      SRCGCB(L,K)=0.D0
      SRCPI0(L,K)=0.D0
  230 CONTINUE
C                                         Water Cloud Size Interpolation
C                                         ------------------------------

      IF(TAUWC(L).GT.TAUWC0) THEN
      SIZWCL=SIZEWC(L)
      LTOPCW=L
      IF(LBOTCW.EQ.0) LBOTCW=L
      IF(SIZWCL.LT.15.D0) THEN
      IF(SIZWCL.LT.3.0D0) SIZWCL=3.0D0
      IRWAT=2
      XRW=SIZWCL/10.0D0-1.00D0
      ELSE
      IF(SIZWCL.GT.25.D0) SIZWCL=25.D0
      IRWAT=4
      XRW=SIZWCL/10.0D0-2.00D0
      ENDIF
      XMW=1.D0-XRW-XRW
      XPW=1.D0+XRW+XRW
      EPS=CLDEPS(L)
      VEP=EPS/(1.D0-EPS)
      VEP1=1.D0+VEP
      TAUWCL=TAUWC(L)
      DO 240 K=1,33
      QAWATK=XMW*XPW*TRCQAB(K,IRWAT)
     +      -XMW*XRW*TRCQAB(K,IRWAT-1)+XPW*XRW*TRCQAB(K,IRWAT+1)
      QPWATK=XMW*XPW*TRCQPI(K,IRWAT)
     +      -XMW*XRW*TRCQPI(K,IRWAT-1)+XPW*XRW*TRCQPI(K,IRWAT+1)
      VEPP=VEP*QPWATK
      TRCALK(L,K)=TRCALK(L,K)+TAUWCL*QAWATK/(VEP1-VEPP)
  240 CONTINUE
      SRCGFW=SRCGSF(1)
      DO 250 K=1,6
      QXWATK=XMW*XPW*SRCQEX(K,IRWAT)
     +      -XMW*XRW*SRCQEX(K,IRWAT-1)+XPW*XRW*SRCQEX(K,IRWAT+1)
      QSWATK=XMW*XPW*SRCQSC(K,IRWAT)
     +      -XMW*XRW*SRCQSC(K,IRWAT-1)+XPW*XRW*SRCQSC(K,IRWAT+1)
      QGWATK=XMW*XPW*SRCQCB(K,IRWAT)
     +      -XMW*XRW*SRCQCB(K,IRWAT-1)+XPW*XRW*SRCQCB(K,IRWAT+1)
      QPWATK=XMW*XPW*SRCQPI(K,IRWAT)
     +      -XMW*XRW*SRCQPI(K,IRWAT-1)+XPW*XRW*SRCQPI(K,IRWAT+1)
      QGWATK=QGWATK*SRCGFW
      VEPP=VEP*QPWATK
      VEP2=VEP1-VEPP
      SRCEXT(L,K)=SRCEXT(L,K)+TAUWCL*QXWATK/VEP1
      SRCSCT(L,K)=TAUWCL*QSWATK/(VEP1*VEP2)
      SRCGCB(L,K)=QGWATK*VEP2/(VEP1-VEPP*QGWATK)
  250 CONTINUE
      ENDIF
C                                           Ice Cloud Size Interpolation
C                                           ----------------------------
      IF(TAUIC(L).GT.TAUIC0) THEN
      SIZICL=SIZEIC(L)
      LTOPCI=L
      IF(LBOTCI.EQ.0) LBOTCI=L
      IF(SIZICL.LT.25.D0) THEN
      IF(SIZICL.LT.3.0D0) SIZICL=3.0D0
      IRICE=2+ICE012*5
      XRI=SIZICL/20.D0-0.75D0
      ELSE
      IF(SIZICL.GT.75.D0) SIZICL=75.D0
      IRICE=4+ICE012*5
      XRI=SIZICL/50.D0-1.00D0
      ENDIF
      XMI=1.D0-XRI-XRI
      XPI=1.D0+XRI+XRI
      EPS=CLDEPS(L)
      VEP=EPS/(1.D0-EPS)
      VEP1=1.D0+VEP
      TAUICL=TAUIC(L)
      DO 260 K=1,33
      QAICEK=XMI*XPI*TRCQAB(K,IRICE)
     +      -XMI*XRI*TRCQAB(K,IRICE-1)+XPI*XRI*TRCQAB(K,IRICE+1)
      QPICEK=XMI*XPI*TRCQPI(K,IRICE)
     +      -XMI*XRI*TRCQPI(K,IRICE-1)+XPI*XRI*TRCQPI(K,IRICE+1)
      VEPP=VEP*QPICEK
      TRCALK(L,K)=TRCALK(L,K)+TAUICL*QAICEK/(VEP1-VEPP)
  260 CONTINUE

      SRCGFC=SRCGSF(2)
      IF(ICE012.EQ.2) SRCGFC=SRCGSF(3)
      DO 270 K=1,6
      QXICEK=XMI*XPI*SRCQEX(K,IRICE)
     +      -XMI*XRI*SRCQEX(K,IRICE-1)+XPI*XRI*SRCQEX(K,IRICE+1)
      QSICEK=XMI*XPI*SRCQSC(K,IRICE)
     +      -XMI*XRI*SRCQSC(K,IRICE-1)+XPI*XRI*SRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*SRCQCB(K,IRICE)
     +      -XMI*XRI*SRCQCB(K,IRICE-1)+XPI*XRI*SRCQCB(K,IRICE+1)
      QPICEK=XMI*XPI*SRCQPI(K,IRICE)
     +      -XMI*XRI*SRCQPI(K,IRICE-1)+XPI*XRI*SRCQPI(K,IRICE+1)
      QGICEK=QGICEK*SRCGFC
      VEPP=VEP*QPICEK
      VEP2=VEP1-VEPP
      SRCEXT(L,K)=SRCEXT(L,K)+TAUICL*QXICEK/VEP1
      SCTTAU=TAUICL*QSICEK/(VEP1*VEP2)
      GCBICE=QGICEK*VEP2/(VEP1-VEPP*QGICEK)
      SCTGCB=SRCSCT(L,K)*SRCGCB(L,K)+SCTTAU*GCBICE
      SRCSCT(L,K)=SRCSCT(L,K)+SCTTAU
      SRCGCB(L,K)=SCTGCB/SRCSCT(L,K)
  270 CONTINUE
      ENDIF
  280 CONTINUE

C     ------------------------------------------------------------------
C     Identify Top Cloud (LTOPCL) and define top cloud albedo correction
C
C     Full Scattering Correction:      KCLDEM=1   ECLTRA=1.0   (default)
C     Partial(rad99a) Correction:      KCLDEM=0   ECLTRA=1.0
C       No Scattering Correction:      KCLDEM=0   ECLTRA=0.0
C
C     KCLDEM=1 Top-cloud scattering correction uses TXCTPG,TSCTPG,TGCTPG
C              to generate correction (over-rides old ECLTRA correction)
C              (KCLDEM correction is computed in THRMAL at LTOPCL level)
C     ------------------------------------------------------------------

      LTOPCL=LTOPCI
      IF(LTOPCI.GT.LTOPCW) GO TO 330
      IF(LTOPCW.LT.1) GO TO 350
      LTOPCL=LTOPCW
      TCTAUW=TAUWC(LTOPCL)
      DO 310 K=1,33
      ALWATK=XMW*XPW*TRCQAL(K,IRWAT)
     +      -XMW*XRW*TRCQAL(K,IRWAT-1)+XPW*XRW*TRCQAL(K,IRWAT+1)
      QXWATK=XMW*XPW*TRCQEX(K,IRWAT)
     +      -XMW*XRW*TRCQEX(K,IRWAT-1)+XPW*XRW*TRCQEX(K,IRWAT+1)
      TRCTCA(K)=(1.D0-EXP(-TAUWC(LTOPCL)*QXWATK))*ALWATK*ECLTRA
      QSWATK=XMW*XPW*TRCQSC(K,IRWAT)
     +      -XMW*XRW*TRCQSC(K,IRWAT-1)+XPW*XRW*TRCQSC(K,IRWAT+1)
      QGWATK=XMW*XPW*TRCQCB(K,IRWAT)
     +      -XMW*XRW*TRCQCB(K,IRWAT-1)+XPW*XRW*TRCQCB(K,IRWAT+1)
      TXCTPG(K)=QXWATK*TCTAUW
      TSCTPG(K)=QSWATK*TCTAUW
      TGCTPG(K)=QGWATK
  310 CONTINUE
      LBOTCL=LBOTCW
      IF(LBOTCI.LT.1) GO TO 360
      IF(LBOTCI.LE.LBOTCW) LBOTCL=LBOTCI
      IF(LTOPCI.EQ.LTOPCW) THEN
      TCTAUW=TAUWC(LTOPCL)
      TCTAUC=TAUIC(LTOPCL)
      WTI=TAUIC(LTOPCL)/(TAUIC(LTOPCL)+TAUWC(LTOPCL))
      WTW=TAUWC(LTOPCL)/(TAUIC(LTOPCL)+TAUWC(LTOPCL))
      DO 320 K=1,33
      ALICEK=XMI*XPI*TRCQAL(K,IRICE)
     +      -XMI*XRI*TRCQAL(K,IRICE-1)+XPI*XRI*TRCQAL(K,IRICE+1)
      QXICEK=XMI*XPI*TRCQEX(K,IRICE)
     +      -XMI*XRI*TRCQEX(K,IRICE-1)+XPI*XRI*TRCQEX(K,IRICE+1)
      TRCTCI=(1.D0-EXP(-TAUIC(LTOPCL)*QXICEK))*ALICEK*ECLTRA
      TRCTCA(K)=WTW*TRCTCA(K)+WTI*TRCTCI
      QSICEK=XMI*XPI*TRCQSC(K,IRICE)
     +      -XMI*XRI*TRCQSC(K,IRICE-1)+XPI*XRI*TRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*TRCQCB(K,IRICE)
     +      -XMI*XRI*TRCQCB(K,IRICE-1)+XPI*XRI*TRCQCB(K,IRICE+1)
      TXCTPG(K)=TXCTPG(K)+QXICEK*TCTAUC
      SCTGCB=TSCTPG(K)*TGCTPG(K)+QSICEK*TCTAUC*QGICEK
      TSCTPG(K)=TSCTPG(K)+QSICEK*TCTAUC
      TGCTPG(K)=SCTGCB/(1.D-10+TSCTPG(K))
  320 CONTINUE
      ENDIF
      GO TO 360
  330 CONTINUE
      LTOPCL=LTOPCI
      TCTAUC=TAUIC(LTOPCL)
      DO 340 K=1,33
      ALICEK=XMI*XPI*TRCQAL(K,IRICE)
     +      -XMI*XRI*TRCQAL(K,IRICE-1)+XPI*XRI*TRCQAL(K,IRICE+1)
      QXICEK=XMI*XPI*TRCQEX(K,IRICE)
     +      -XMI*XRI*TRCQEX(K,IRICE-1)+XPI*XRI*TRCQEX(K,IRICE+1)
      TRCTCA(K)=(1.D0-EXP(-TAUIC(LTOPCL)*QXICEK))*ALICEK*ECLTRA
      QSICEK=XMI*XPI*TRCQSC(K,IRICE)
     +      -XMI*XRI*TRCQSC(K,IRICE-1)+XPI*XRI*TRCQSC(K,IRICE+1)
      QGICEK=XMI*XPI*TRCQCB(K,IRICE)
     +      -XMI*XRI*TRCQCB(K,IRICE-1)+XPI*XRI*TRCQCB(K,IRICE+1)
      TXCTPG(K)=QXICEK*TCTAUC
      TSCTPG(K)=QSICEK*TCTAUC
      TGCTPG(K)=QGICEK
  340 CONTINUE
      LBOTCL=LBOTCI
      IF(LBOTCW.EQ.0) GO TO 360
      LBOTCL=LBOTCW
      GO TO 360
  350 CONTINUE
      LBOTCL=0
      LTOPCL=0
  360 CONTINUE

      RETURN


C--------------------------------
      ENTRY UPDEPS(JYEARE,JJDAYE)
C--------------------------------

C                 Select ISCCP-Based Cloud Heterogeneity Time Dependence
C                 ------------------------------------------------------

      XJDAY=JJDAYE-0.999D0
      XMO=XJDAY/30.5D0+.5D0
      MI=XMO
      WTMJ=XMO-MI
      WTMI=1.D0-WTMJ
      IF(MI.LT.1) MI=12
      MJ=MI+1
      IF(MJ.GT.12) MJ=1

      DO 420 J=1,46
      DO 410 I=1,72
      EPLOW(I,J)=WTMI*EPLMHC(I,J,MI,1)+WTMJ*EPLMHC(I,J,MJ,1)
      EPMID(I,J)=WTMI*EPLMHC(I,J,MI,2)+WTMJ*EPLMHC(I,J,MJ,2)
      EPHIG(I,J)=WTMI*EPLMHC(I,J,MI,3)+WTMJ*EPLMHC(I,J,MJ,3)
      EPCOL(I,J)=WTMI*EPLMHC(I,J,MI,4)+WTMJ*EPLMHC(I,J,MJ,4)
  410 CONTINUE
  420 CONTINUE

      RETURN


C-----------------
      ENTRY GETEPS
C-----------------

C             ----------------------------------------------------------
C                     Select Cloud Heterogeneity CLDEPS Options
C             EPSCON  Column Cloud Inhomogeneity EPSILON (when KCLDEP=1)
C             KCLDEP  Selects Cloud Inhomogeneity Option (0-4):
C                     KCLDEP =  0  Sets Column CLDEPS to Zero
C                     KCLDEP =  1  Sets Column CLDEPS to EPSCON
C                     KCLDEP =  2  Keeps whatever is specified in CLDEPS
C                     KCLDEP =  3  Uses: Column EPCOL(72,46) Climatology
C                     KCLDEP =  4  Uses: Ht Dep EPLOW, EPMID, EPHIG Data
C                     --------------------------------------------------

      IF(KCLDEP.EQ.0) THEN
      DO 510 L=1,NL
      CLDEPS(L)=0.D0
  510 CONTINUE
      ENDIF
      IF(KCLDEP.EQ.1) THEN
      DO 520 L=1,NL
      CLDEPS(L)=EPSCON
  520 CONTINUE
      ENDIF
      I=ILON
      J=JLAT
      IF(KCLDEP.EQ.3) THEN
      DO 530 L=1,NL
      CLDEPS(L)=EPCOL(I,J)
  530 CONTINUE
      ENDIF
      IF(KCLDEP.EQ.4) THEN
      DO 540 L=1,NL
      CLDEPS(L)=EPMID(I,J)
      IF(PLB(L).GT.750.D0) CLDEPS(L)=EPLOW(I,J)
      IF(PLB(L).LT.430.D0) CLDEPS(L)=EPHIG(I,J)
  540 CONTINUE
      ENDIF

      RETURN
      END SUBROUTINE SETCLD


      SUBROUTINE TAUGAS
      IMPLICIT NONE
C     ----------------------------------------------------------
C     TAUGAS INPUT REQUIRES:  NL,PL,DPL,TLM,ULGAS
C                             TAUTBL,TAUWV0,XKCFC,H2OCN8,H2OCF8
C                             XUCH4,XUCH40,XUN2O,XUN2O0,ULOX,DUX
C                             XTRUP,XTU0,XTRDN,XTD0,CXUCO2,CXUO3
C     TAUGAS OUTPUT DATA IS:  TRGXLK,XTRU,XTRD
C     ----------------------------------------------------------

      INTEGER, PARAMETER :: NTX=8, NPX=19, NGUX=1008, NPUX=19, NPU2=14,
     *     NPU=5
      REAL*8, PARAMETER :: TLOX=181.d0, DTX=23.d0, P0=1013.25d0

      REAL*8, PARAMETER :: PX(19)= (/1d3,750d0,5d2,3d2,2d2,1d2,5d1
     *     ,2d1,1d1,5d0,2d0,1d0,.5d0,.2d0,.1d0,.03d0,.01d0,.003d0,.001d0
     *     /)

      INTEGER, PARAMETER :: NGX(4) = (/12,12,08,33/),
     *     IG1X(4) = (/2,14,26,1/)
      REAL*8, PARAMETER :: PDPU2(14) = (/1.D4,1.D5,2.D5,5.D5,1.D6,2.D6,
     *     5.D6,1.D7,2.D7,5.D7,1.D8,2.D8,5.D8,1.D9/)
      REAL*8, PARAMETER ::  PU(5) = (/  50.,200.,800.,3200.,12800./)
      INTEGER, PARAMETER :: IGASX(20) = (/ 1, 2, 3, 1, 1, 2, 2, 3, 3, 6,
     *     6, 6, 7,13,13, 8, 8, 9, 9, 1/)
      INTEGER, PARAMETER :: KGX(20) = (/ 1, 2, 3, 2, 3, 1, 3, 1, 2, 1, 2
     *     , 3, 1, 1,3, 2, 3, 2, 3, 4/)
      INTEGER, PARAMETER :: NUX(15) = (/25, 9, 9, 9, 9, 5, 5, 5, 5, 2, 2
     *     , 2, 2, 2,2/)
      INTEGER, PARAMETER :: IGUX(15) = (/ 0,300,408,480,588,660,720,760
     *     ,820,880,904,928,944,968,992/)

      REAL*8, PARAMETER, DIMENSION(8,2) ::  XKCFCW = RESHAPE( (/
     + 12.4414,11.7842,11.3630,10.8109,10.3200, 9.8900, 9.3916, 8.8933,
     +  5.3994, 5.6429, 5.8793, 6.1687, 6.2300, 6.5200, 6.8650, 7.2100/)
     *     , (/8,2/) )

      REAL*8, PARAMETER ::  P24(24) = (/
     $ .100D+04,.973D+03,.934D+03,.865D+03,.752D+03,.603D+03,
     $ .439D+03,.283D+03,.156D+03,.754D+02,.350D+02,.162D+02,
     $ .754D+01,.350D+01,.162D+01,.743D+00,.340D+00,.152D+00,
     $ .701D-01,.347D-01,.159D-01,.750D-02,.350D-02,.100D-02/)
      REAL*8, PARAMETER ::  DP24(24) = (/
     $     24.4,32.0,46.6,89.8,136.8,162.0,165.4,146.9,106.5,55.2,25.6
     *     ,11.9,5.52,2.56,1.20,.551,.256,.119,.0452,.0256,.0119,.005,
     *     .003,.002/)
      REAL*8, PARAMETER ::  DLSQ2 =.1505d0,  ULMNH2=2.060d0,
     *     ULMNCH=-1.028d0, ULMNN2=-1.530d0, DLOG2 =.3010d0,
     *     ULMNO3=-1.393d0, ULMNCO= 1.529d0
      REAL*8 XTU(24,3),XTD(24,3),XUCH(9),XUN2(9),CXUO(7),CXUC(7)
      INTEGER MLGAS(20)
      INTEGER I,K,L,IP,IM,IULOW,IU1,IU2,IU,II,IPX,IAA,IBA,ITX,ITX1
     *     ,ITX2,IPXM1,IGAS,NG,KK,IK0,IKF,IPU,IPUI,IPU1,IK,IGCFC,NU,IUA
     *     ,IUB,IAAA,IAAB,IABA,IABB,IBAA,IBAB,IBBA,IBBB,IH2O0,IG
      REAL*8 UH2O,UCO2,UO3,UCH4,UN2O,UH2OL,UCO2L,UO3LL,UCH4L,UN2OL,CXCO2
     *     ,CXO3,DUH2,DU1,DU2,DUCO,DUO3,DUCH,XCH4,DUN2,XN2O,PLI,DENO
     *     ,ANUM,DELP24,PRAT,XTR0,WPB,WTB,WBB,WBA,WAB,WAA,UGAS,XF
     *     ,FPL,U,PU2,WTPU,TAUT1,TAUT2,TAUHCN,UP,TAUHFB,XA,XB,XK,TAUCF
     *     ,XUA,XUB,QAA,QAB,QBA,QBB,FNU1,UAB,UBB,UAA,UBA,WAAA,WAAB,WABA
     *     ,WABB,WBAA,WBAB,WBBA,WBBB,TAUIPG,TAUSUM,TAU11,TAU12

      DO 100 K=1,33
      DO 100 L=1,NL
      TRGXLK(L,K)=0.D0
  100 CONTINUE

      DO 110 I=1,20
      MLGAS(I)=1
  110 CONTINUE

C              KWVCON = ON/OFF flag for water vapor continuum absorption
C              ---------------------------------------------------------
      IF(KWVCON.LT.1) MLGAS(18)=0

      UH2O=1.D-10
      UCO2=1.D-10
      UO3=1.D-10
      UCH4=1.D-10
      UN2O=1.D-10
      DO 120 IP=1,NL
      UH2O=UH2O+ULGAS(IP,1)
      UCO2=UCO2+ULGAS(IP,2)
      UO3=UO3+ULGAS(IP,3)
      UCH4=UCH4+ULGAS(IP,7)
      UN2O=UN2O+ULGAS(IP,6)
  120 CONTINUE
      UH2OL=LOG10(UH2O)
      UCO2L=LOG10(UCO2)
      UO3LL=LOG10(UO3)
      UCH4L=LOG10(UCH4)
      UN2OL=LOG10(UN2O)

      IULOW=0
      IF(UH2O.LT.1.1D-10) THEN
      IULOW=1
      DO 130 IU=1,9
      XUCH(IU)=XUCH40(IU)
      XUN2(IU)=XUN2O0(IU)
  130 CONTINUE
      DO 140 IM=1,3
      DO 140 I=1,24
      XTU(I,IM)=XTU0(I,IM)
      XTD(I,IM)=XTD0(I,IM)
  140 CONTINUE
      CXCO2=0.D0
      CXO3=0.D0
      GO TO 180
      ENDIF

      DUH2=UH2OL-ULMNH2
      IF(DUH2.LT.0.) DUH2=0.D0
      IU1=DUH2/DLSQ2+1.D0
      IF(IU1.LT.1) IU1=1
      IF(IU1.GT.14) IU1=14
      IU2=IU1+1
      DU1=DUH2-(IU1-1)*DLSQ2
      DU2=DLSQ2-DU1

      DO 160 IM=1,3
      DO 150 I=1,24
      XTU(I,IM)=(XTRUP(I,IM,IU2)*DU1+XTRUP(I,IM,IU1)*DU2)/DLSQ2
      XTD(I,IM)=(XTRDN(I,IM,IU2)*DU1+XTRDN(I,IM,IU1)*DU2)/DLSQ2
  150 CONTINUE
  160 CONTINUE
      DO 170 IU=1,9
      XUCH(IU)=(XUCH4(IU,IU2)*DU1+XUCH4(IU,IU1)*DU2)/DLSQ2
      XUN2(IU)=(XUN2O(IU,IU2)*DU1+XUN2O(IU,IU1)*DU2)/DLSQ2
      IF(IU.GT.7) GO TO 170
      CXUO(IU)=(CXUO3(IU,IU2)*DU1+CXUO3(IU,IU1)*DU2)/DLSQ2
      CXUC(IU)=(CXUCO2(IU,IU2)*DU1+CXUCO2(IU,IU1)*DU2)/DLSQ2
  170 CONTINUE

      DUCO=UCO2L-ULMNCO
      IF(DUCO.LT.0.) DUCO=0.
      IU1=DUCO/DLOG2+1
      IF(IU1.LT.1) IU1=1
      IF(IU1.GT.6) IU1=6
      IU2=IU1+1
      DU1=DUCO-(IU1-1)*DLOG2
      DU2=DLOG2-DU1
      CXCO2=(CXUC(IU2)*DU1+CXUC(IU1)*DU2)/DLOG2

      DUO3=UO3LL-ULMNO3
      IF(DUO3.LT.0.) DUO3=0.D0
      IU1=DUO3/DLOG2+1
      IF(IU1.LT.1) IU1=1
      IF(IU1.GT.6) IU1=6
      IU2=IU1+1
      DU1=DUO3-(IU1-1)*DLOG2
      DU2=DLOG2-DU1
      CXO3=(CXUO(IU2)*DU1+CXUO(IU1)*DU2)/DLOG2

  180 CONTINUE

      DUCH=UCH4L-ULMNCH
      IF(DUCH.LT.0.) DUCH=0.D0
      IU1=DUCH/DLOG2+1
      IF(IU1.LT.1) IU1=1
      IF(IU1.GT.8) IU1=8
      IU2=IU1+1
      DU1=DUCH-(IU1-1)*DLOG2
      DU2=DLOG2-DU1
      XCH4=(XUCH(IU2)*DU1+XUCH(IU1)*DU2)/DLOG2

      DUN2=UN2OL-ULMNN2
      IF(DUN2.LT.0.) DUN2=0.D0
      IU1=DUN2/DLOG2+1
      IF(IU1.LT.1) IU1=1
      IF(IU1.GT.8) IU1=8
      IU2=IU1+1
      DU1=DUN2-(IU1-1)*DLOG2
      DU2=DLOG2-DU1
      XN2O=(XUN2(IU2)*DU1+XUN2(IU1)*DU2)/DLOG2

      DO 190 I=1,NL
      XTRU(I,1)=1.D0
      XTRD(I,1)=1.D0
  190 CONTINUE

      IP=2
      DO 230 I=1,NL
      PLI=PL(I)
      IF(PLI.GE.P24(1)) THEN
      DO 200 IM=1,3
      XTRU(I,IM+1)=1.D0-(1.D0-XTU(1,IM))*DPL(I)/DP24(1)
      XTRD(I,IM+1)=1.D0-(1.D0-XTD(1,IM))*DPL(I)/DP24(1)
  200 CONTINUE
      GO TO 230
      ENDIF
  210 IF(PLI.GE.P24(IP)) THEN
      DENO=P24(IP)-P24(IP-1)
      ANUM=DP24(IP-1)*(P24(IP)-PLI)+DP24(IP)*(PLI-P24(IP-1))
      DELP24=ANUM/DENO
      PRAT=DPL(I)/DELP24
      DO 220 IM=1,3
      ANUM=XTU(IP-1,IM)*(P24(IP)-PLI)+XTU(IP,IM)*(PLI-P24(IP-1))
      XTR0=ANUM/DENO
      XTRU(I,IM+1)=1.D0-(1.D0-XTR0)*PRAT
      ANUM=XTD(IP-1,IM)*(P24(IP)-PLI)+XTD(IP,IM)*(PLI-P24(IP-1))
      XTR0=(ANUM/DENO)
      XTRD(I,IM+1)=1.D0-(1.D0-XTR0)*PRAT
  220 CONTINUE
      ELSE
      IP=IP+1
      IF(IP.GT.24) GO TO 240
      GO TO 210
      ENDIF
  230 CONTINUE
      GO TO 260
  240 CONTINUE
      DO 250 IM=1,3
      DO 250 II=I,NL
      XTRU(II,IM+1)=XTU(24,IM)
      XTRD(II,IM+1)=XTD(24,IM)
  250 CONTINUE
  260 CONTINUE
      DO 270 IM=2,4
      XTRD(NL,IM)=1.D0
  270 CONTINUE

      IPX=2
      DO 600 IP=1,NL
  280 CONTINUE
      WPB = (PL(IP)-PX(IPX))/(PX(IPX-1)-PX(IPX))
      IF(WPB.GE.0.D0.OR.IPX.GE.NPX) GO TO 290
      IPX = IPX+1
      GO TO 280
  290 CONTINUE
      WTB = (TLM(IP)-TLOX)/DTX
      ITX = MIN(MAX(INT(WTB),0),NTX-2)
      WTB = WTB-FLOAT(ITX)

      WBB = WPB*WTB
      WBA = WPB-WBB
      WAB = WTB-WBB
      WAA = 1.D0-(WBB+WBA+WAB)

      IAA = NGUX*(ITX+NTX*(IPX-1))
      IBA = IAA-NGUX*NTX

      ITX1=ITX+1
      ITX2=ITX+2
      IPXM1=IPX-1

      DO 500 IGAS=1,20
      IF(MLGAS(IGAS).LT.1) GO TO 500
      NG = NGX(KGX(IGAS))
      UGAS = ULGAS(IP,IGASX(IGAS))
      KK=IG1X(KGX(IGAS))
C                       Apply absorber scaling for H2O in CO2 & O3 bands
C                       ------------------------------------------------
      IF(IGAS.EQ.4) THEN
      IF(PL(IP).GT.350.D0) THEN
      XF=1.D0+CXCO2*(PL(IP)-350.D0)
      IF(XF.LT.0.D0) XF=0.D0
      UGAS=UGAS*XF
      ENDIF
      ENDIF

      IF(IGAS.EQ.5) THEN
      IF(PL(IP).GT.350.D0) THEN
      XF=1.D0+CXO3*(PL(IP)-350.D0)
      IF(XF.LT.0.D0) XF=0.D0
      UGAS=UGAS*XF
      ENDIF
      ENDIF

C                                         Modified scaling for CH4 & N2O
C                                         ------------------------------
      IF(PL(IP).GT.100.D0) THEN
      FPL=MIN((PL(IP)-100.D0)/300.D0,1.D0)
      IF(IGAS.EQ.13) UGAS=UGAS*(1.D0+FPL*(XCH4-1.D0))
      IF(IGAS.GT.9.AND.IGAS.LT.13)
     +               UGAS=UGAS*(1.D0+FPL*(XN2O-1.D0))
      ENDIF

C                                               Absorber scaling for SO2
C                                               ------------------------
      IF(IGAS.EQ.14.AND.IULOW.EQ.0) UGAS=UGAS*3.
      IF(IGAS.EQ.14.AND.IULOW.EQ.1) UGAS=UGAS*2.
      IF(IGAS.EQ.15) UGAS=UGAS*1.3

C                                 Apply water vapor continuum absorption
C                                 --------------------------------------
      IF(IGAS.EQ.20) THEN

C                  KWSELF = ON/FF flag for H2O self broadening continuum
C                  -----------------------------------------------------
      IF(KCSELF.GT.0) THEN
      U=UGAS*1.15D0
      IK0=1
      IKF=1
      DO 330 I=1,2
      IF(I.EQ.2) THEN
      U=UGAS
      IK0=2
      IKF=33
      ENDIF
      PU2=PL(IP)/DPL(IP)*U**2
      IF(PU2.GT.PDPU2(1)) THEN
      DO 300 IPU=2,NPU2
      IPUI=IPU
      IF(PU2.LE.PDPU2(IPU)) GO TO 310
  300 CONTINUE
  310 CONTINUE
      IPU=IPUI
      IPU1=IPU-1
      WTPU=(PU2-PDPU2(IPU1))/(PDPU2(IPU)-PDPU2(IPU1))
      ELSE
      WTPU=PU2/PDPU2(1)
      ENDIF

      DO 320 IK=IK0,IKF
      IF(PU2.GT.PDPU2(1)) THEN
      TAUT1=WTPU*(H2OCN8(IK,ITX1,IPU)-H2OCN8(IK,ITX1,IPU1))+
     $      H2OCN8(IK,ITX1,IPU1)
      TAUT2=WTPU*(H2OCN8(IK,ITX2,IPU)-H2OCN8(IK,ITX2,IPU1))+
     $      H2OCN8(IK,ITX2,IPU1)
      ELSE
      TAUT1=WTPU*H2OCN8(IK,ITX1,1)
      TAUT2=WTPU*H2OCN8(IK,ITX2,1)
      ENDIF

      TAUHCN=WTB*(TAUT2-TAUT1)+TAUT1
      TRGXLK(IP,KK)=TRGXLK(IP,KK)+TAUHCN
      KK=KK+1
  320 CONTINUE
  330 CONTINUE
      ENDIF
C               KWFORN = ON/FF flag for H2O foreign broadening continuum
C               --------------------------------------------------------

      IF(KCFORN.LT.1) GO TO 500
      KK=IG1X(KGX(IGAS))
      U=UGAS*1.15D0
      IK0=1
      IKF=1
      DO 370 I=1,2
      IF(I.EQ.2) THEN
      U=UGAS
      IK0=2
      IKF=33
      ENDIF
      UP=PL(IP)/P0*U
      IF(UP.GT.PU(1)) THEN
      DO 340 IPU=2,NPU
      IPUI=IPU
      IF(UP.LE.PU(IPU)) GO TO 350
  340 CONTINUE
  350 CONTINUE
      IPU=IPUI
      IPU1=IPU-1
      WTPU=(UP-PU(IPU1))/(PU(IPU)-PU(IPU1))
      ELSE
      WTPU=UP/PU(1)
      ENDIF
      DO 360 IK=IK0,IKF
      IF(UP.GT.PU(1)) THEN
      TAUT1=WTPU*(H2OCF8(IK,ITX1,IPU)-H2OCF8(IK,ITX1,IPU1))+
     +      H2OCF8(IK,ITX1,IPU1)
      TAUT2=WTPU*(H2OCF8(IK,ITX2,IPU)-H2OCF8(IK,ITX2,IPU1))+
     +      H2OCF8(IK,ITX2,IPU1)
      ELSE
      TAUT1=WTPU*H2OCF8(IK,ITX1,1)
      TAUT2=WTPU*H2OCF8(IK,ITX2,1)
      ENDIF
      TAUHFB=WTB*(TAUT2-TAUT1)+TAUT1
      TRGXLK(IP,KK)=TRGXLK(IP,KK)+TAUHFB
      KK=KK+1
  360 CONTINUE
  370 CONTINUE
      GO TO 500
      ENDIF

      IF(IGAS.GT.15) THEN
      IGCFC=IGAS-15
      UGAS=UGAS*1.85
      DO 380 IK=1,NG
      XA=WTB*(XKCFC(IK,ITX2,IGCFC)-XKCFC(IK,ITX1,IGCFC))+
     $ XKCFC(IK,ITX1,IGCFC)
      XB=WTB*(XKCFC(IK,ITX2,IGCFC)-XKCFC(IK,ITX1,IGCFC))+
     $ XKCFC(IK,ITX1,IGCFC)
      XK=WPB*(XA-XB)+XB
      TAUCF=XK*UGAS
      TRGXLK(IP,KK)=TRGXLK(IP,KK)+TAUCF
      KK=KK+1
  380 CONTINUE
      GO TO 500
      ENDIF
      IU = IPX + NPUX*(IGAS-1)
      NU = NUX(IGAS)
      IF(NU.GT.1) GO TO 390
      XUA = 0.D0
      XUB = 0.D0
      GO TO 400
  390 CONTINUE
      XUA = (UGAS-ULOX(IU))/DUX(IU)
      XUB = (UGAS-ULOX(IU-1))/DUX(IU-1)
  400 CONTINUE
      IUA = INT(XUA)
      IUB = INT(XUB)

      QAA = 1.D0
      QAB = 1.D0
      IF(XUA.GT.0.D0.AND.IUA.LT.NU-1) GO TO 410
      FNU1=NU-1
      XUA = MIN(MAX(XUA,0.D0),FNU1)
      IUA = MIN(INT(XUA),NU-2)
      QAA = UGAS/(ULOX(IU)+DUX(IU)*FLOAT(IUA))
      QAB = UGAS/(ULOX(IU)+DUX(IU)*FLOAT(IUA+1))
  410 CONTINUE
      QBA = 1.D0
      QBB = 1.D0
      IF(XUB.GT.0.D0.AND.IUB.LT.NU-1) GO TO 420
      FNU1=NU-1
      XUB = MIN(MAX(XUB,0.D0),FNU1)
      IUB = MIN(INT(XUB),NU-2)
      QBA = UGAS/(ULOX(IU-1)+DUX(IU-1)*FLOAT(IUB))
      QBB = UGAS/(ULOX(IU-1)+DUX(IU-1)*FLOAT(IUB+1))
  420 CONTINUE
      UAB = XUA-FLOAT(IUA)
      UBB = XUB-FLOAT(IUB)
      UAA = 1.D0-UAB
      UBA = 1.D0-UBB

      WAAA = WAA*UAA*QAA
      WAAB = WAA*UAB*QAB
      WABA = WAB*UAA*QAA
      WABB = WAB*UAB*QAB
      WBAA = WBA*UBA*QBA
      WBAB = WBA*UBB*QBB
      WBBA = WBB*UBA*QBA
      WBBB = WBB*UBB*QBB

      IAAA = IAA+IGUX(IGAS) + NG*IUA
      IAAB = IAAA+NG
      IABA = IAAA+NGUX
      IABB = IABA+NG
      IBAA = IBA+IGUX(IGAS) + NG*IUB
      IBAB = IBAA+NG
      IBBA = IBAA+NGUX
      IBBB = IBBA+NG

      IH2O0=0
      IF(IGAS.EQ.6.OR.IGAS.EQ.8.OR.IGAS.EQ.10.OR.IGAS.EQ.13.OR.
     +   IGAS.EQ.14) THEN
      IF(IULOW.EQ.1) IH2O0=1
      ENDIF

      DO 430 IG=1,NG
      IF(IH2O0.EQ.0) THEN
      TAUIPG=
     +       WAAA*TAUTBL(IAAA+IG)
     +     + WAAB*TAUTBL(IAAB+IG)
     +     + WABA*TAUTBL(IABA+IG)
     +     + WABB*TAUTBL(IABB+IG)
     +     + WBAA*TAUTBL(IBAA+IG)
     +     + WBAB*TAUTBL(IBAB+IG)
     +     + WBBA*TAUTBL(IBBA+IG)
     +     + WBBB*TAUTBL(IBBB+IG)
      ELSE
      TAUIPG=
     +       WAAA*TAUWV0(IAAA+IG)
     +     + WAAB*TAUWV0(IAAB+IG)
     +     + WABA*TAUWV0(IABA+IG)
     +     + WABB*TAUWV0(IABB+IG)
     +     + WBAA*TAUWV0(IBAA+IG)
     +     + WBAB*TAUWV0(IBAB+IG)
     +     + WBBA*TAUWV0(IBBA+IG)
     +     + WBBB*TAUWV0(IBBB+IG)
      ENDIF

      TAUSUM=TRGXLK(IP,KK)+TAUIPG
      IF(TAUSUM.GT.0.D0) TRGXLK(IP,KK)=TAUSUM
      KK=KK+1
  430 CONTINUE
  500 CONTINUE
C                               CFC11 and CFC12 Window Absorption (1997)
C                               ----------------------------------------

      IF(MLGAS(16).EQ.1.OR.MLGAS(17).EQ.1) THEN
      XK=WTB*(XKCFCW(ITX2,1)-XKCFCW(ITX1,1))+XKCFCW(ITX1,1)
      TAU11=XK*ULGAS(IP,8)
      TRGXLK(IP,1)=TRGXLK(IP,1)+TAU11
      ENDIF
      IF(MLGAS(18).EQ.1.OR.MLGAS(19).EQ.1) THEN
      XK=WTB*(XKCFCW(ITX2,2)-XKCFCW(ITX1,2))+XKCFCW(ITX1,2)
      TAU12=XK*ULGAS(IP,9)
      TRGXLK(IP,1)=TRGXLK(IP,1)+TAU12
      ENDIF
  600 CONTINUE

      RETURN
      END SUBROUTINE TAUGAS



      SUBROUTINE THERML
      IMPLICIT NONE
C     ------------------------------------------------------------------
C             Top-cloud Thermal Scattering Correction Control Parameters
C             ----------------------------------------------------------
C
C             ECLTRA = 1.0  Scattering correction is enabled
C             with KCLDEM = 1, Rigorous scattering correction is applied
C             with KCLDEM = 0, Approximate scattering correction is used
C
C             ECLTRA = 0.0  No scattering correction is used
C                                          (Independent of KCLDEM value)
C
C     ------------------------------------------------------------------
C                                   Lower Edge Temperature Interpolation
C                                   ------------------------------------
C     TLGRAD=1.0  (Default)
C                 Layer-mean temperatures (TLM) supplied by GCM are used
C                 to define the layer edge temperature TLT (top) and TLB
C                 (bottom) using overall atmospheric temperature profile
C                 to establish temperature gradient within each layer so
C                 as to minimize the temperature discontinuities between
C                 layer edges and to conserve layer thermal energy.
C
C     TLGRAD=0.0  This results in isothermal layers with TLT = TLB = TLM
C
C     TLGRAD<0.0  TLT and TLB are used as specified, without any further
C                 adjustments.  This is mainly for off-line use when the
C                 temperature profile (TLM,TLT,TLB) can be fully defined
C                 from a continuous temperature profile.
C
C     NOTE:       TLGRAD can also accommodate values between 0.0 and 1.0
C
C     PTLISO      (Default PTLISO=2.5mb)
C                 Pressure level above which model layers are defined to
C                 be isothermal.  This is appropriate for optically thin
C                 layers where emitted flux depends on mean temperature.
C     ------------------------------------------------------------------
      REAL*8, PARAMETER :: R6=.16666667D0, R24=4.1666667D-02
      REAL*8, PARAMETER :: A=0.3825D0,B=0.5742D0,C=0.0433D0

      REAL*8 TA,TB,TC,P1,P2,P3,P4,DT1CPT,DTHALF,CLTAUX,CLTAUS,CLCOSB
     *     ,CLPI0K,CTX,DT2,DT1,CTG,DG2,DG1,WT1,WT2,WT3,WT4,WT5,WT6,WT7
     *     ,WT8,BG,DNACUM,DNBCUM,DNCCUM,TAUAG,TAUAP,TAUBP,TAUCP,TAUAX
     *     ,TAUBX,TAUCX,XTRDL,BTOP,BBOT,BBAR,TX,PLBN,F,TAUA,TAUB,TAUC
     *     ,BDIF,BBTA,BBTB,BBTC,TAUA2,TRANA,TAUB2,TRANB,TAUC2,TRANC,DEC
     *     ,DEB,DEA,COALB1,COALB2,COALB3,FDNABC,UNA,UNB,UNC,FUNABC
     *     ,PFW,DPF,CTP,DP1,DP2,TAUBG,TAUCG,DDFLUX,XTRUL
      INTEGER K,L,N,II,ITL,ICT,IT1,IT2,IP1,IP2,ICG,IG1,IG2,ITK0,IMOL
     *     ,IPF,ICP

      IF(TLGRAD.LT.0.D0) GO TO 130
      TA=TLM(1)
      TB=TLM(2)
      P1=PLB(1)
      P2=PLB(2)
      P3=PLB(3)
      DT1CPT=0.5D0*TA*(PLB(1)**0.286D0-PLB(2)**0.286D0)/PL(1)**0.286D0
      DTHALF=(TA-TB)*(P1-P2)/(P1-P3)
      IF(DTHALF.GT.DT1CPT) DTHALF=DT1CPT
      TLB(1)=TA+DTHALF*TLGRAD
      TLT(1)=TA-DTHALF*TLGRAD
      DO 110 L=3,NL
      TC=TLM(L)
      P4=PLB(L+1)
      DTHALF=0.5D0*((TA-TB)/(P1-P3)+(TB-TC)/(P2-P4))*(P2-P3)*TLGRAD
      TLB(L-1)=TB+DTHALF
      TLT(L-1)=TB-DTHALF
      TA=TB
      TB=TC
      P1=P2
      P2=P3
      P3=P4
  110 CONTINUE
      DTHALF=(TA-TB)*(P2-P3)/(P1-P3)*TLGRAD
      TLB(NL)=TC+DTHALF
      TLT(NL)=TC-DTHALF
      DO 120 L=NL,1,-1
      IF(PLB(L).GT.PTLISO) GO TO 130
      TLT(L)=TLM(L)
      TLB(L)=TLM(L)
  120 CONTINUE
  130 CONTINUE
      TLB(NLP)=TLT(NL)

C     ------------------------------------------------------------------
C                   WEIGHT ASSIGNMENTS FOR PLANCK FUNCTION INTERPOLATION
C                    (Effective range is from TK = 123 K to TK = 373 K)
C     ------------------------------------------------------------------

      DO 140 L=1,NL
      ITL=TLB(L)
      WTLB(L)=TLB(L)-ITL
      ITLB(L)=ITL-ITPFT0
      IF(ITLB(L).LT.1) ITLB(L)=1
      IF(ITLB(L).GT.249) ITLB(L)=249
      ITL=TLT(L)
      WTLT(L)=TLT(L)-ITL
      ITLT(L)=ITL-ITPFT0
      IF(ITLT(L).LT.1) ITLT(L)=1
      IF(ITLT(L).GT.249) ITLT(L)=249
  140 CONTINUE

      DO 180 L=1,NL
      IF(L.EQ.LTOPCL) THEN
      DO 170 K=1,33
      CLTAUX=TXCTPG(K)+TRGXLK(L,K)+1.D-10
      CLTAUS=TSCTPG(K)
      CLCOSB=TGCTPG(K)
      CLPI0K=CLTAUS*ECLTRA/CLTAUX
      CLPI0(K)=CLPI0K
      CTX=CLTAUX*10.D0
      IF(CLTAUX.GE.3.D0) THEN
      CTX=CLTAUX*2.D0+24.D0
      IF(CTX.GT.47.999999D0) CTX=47.999999D0
      ENDIF
      ICT=CTX
      DT2=CTX-ICT
      DT1=1.D0-DT2
      IT1=ICT+1
      IT2=ICT+2
      CTP=CLPI0K*20.D0
      ICP=CTP
      DP2=CTP-ICP
      DP1=1.D0-DP2
      IP1=ICP+1
      IP2=ICP+2
      CTG=CLCOSB*20.D0
      ICG=CTG
      DG2=CTG-ICG
      DG1=1.D0-DG2
      IG1=ICG+1
      IG2=ICG+2
      WT1=DT1*DP1*DG1
      WT2=DT2*DP1*DG1
      WT3=DT2*DP2*DG1
      WT4=DT1*DP2*DG1
      WT5=DT1*DP1*DG2
      WT6=DT2*DP1*DG2
      WT7=DT2*DP2*DG2
      WT8=DT1*DP2*DG2
      DO 150 II=1,6
      RIJTCK(II,K)=WT1*RIJTPG(II,IT1,IP1,IG1)+WT2*RIJTPG(II,IT2,IP1,IG1)
     +            +WT3*RIJTPG(II,IT2,IP2,IG1)+WT4*RIJTPG(II,IT1,IP2,IG1)
     +            +WT5*RIJTPG(II,IT1,IP1,IG2)+WT6*RIJTPG(II,IT2,IP1,IG2)
     +            +WT7*RIJTPG(II,IT2,IP2,IG2)+WT8*RIJTPG(II,IT1,IP2,IG2)
  150 CONTINUE
      DO 160 II=1,3
      FEMTCK(II,K)=WT1*FEMTPG(II,IT1,IP1,IG1)+WT2*FEMTPG(II,IT2,IP1,IG1)
     +            +WT3*FEMTPG(II,IT2,IP2,IG1)+WT4*FEMTPG(II,IT1,IP2,IG1)
     +            +WT5*FEMTPG(II,IT1,IP1,IG2)+WT6*FEMTPG(II,IT2,IP1,IG2)
     +            +WT7*FEMTPG(II,IT2,IP2,IG2)+WT8*FEMTPG(II,IT1,IP2,IG2)
      FDXTCK(II,K)=WT1*FDXTPG(II,IT1,IP1,IG1)+WT2*FDXTPG(II,IT2,IP1,IG1)
     +            +WT3*FDXTPG(II,IT2,IP2,IG1)+WT4*FDXTPG(II,IT1,IP2,IG1)
     +            +WT5*FDXTPG(II,IT1,IP1,IG2)+WT6*FDXTPG(II,IT2,IP1,IG2)
     +            +WT7*FDXTPG(II,IT2,IP2,IG2)+WT8*FDXTPG(II,IT1,IP2,IG2)
  160 CONTINUE
  170 CONTINUE
      ENDIF
  180 CONTINUE
      DO 190 L=1,NLP
      TRDFLB(L)=0.D0
      TRUFLB(L)=0.D0
  190 CONTINUE

      BG=BGFEMT(1)
      TOTLZF(1)=0.D0
      TOTLZF(2)=0.D0
      TOTLZF(3)=0.D0
!sl   TRSLTS=0.D0
!sl   TRSLTG=0.D0
!sl   TRSLBS=0.D0

      K=0
      IMOL=0
  200 CONTINUE
      ITK0=K*ITNEXT
      K=K+1
      IF(K.GT.33) GO TO 300
      BG=BGFEMT(K)
      IF(K.GT.1.AND.K.LT.14) IMOL=1
      IF(K.GT.13.AND.K.LT.26) IMOL=2
      IF(K.GT.25) IMOL=3
      DFLB(NLP,K)=0.D0
      DNACUM=0.D0
      DNBCUM=0.D0
      DNCCUM=0.D0
      L=NL
  205 CONTINUE
      TAUAG=TRGXLK(L,K)
      TAUAP=TRCALK(L,K)+TRAALK(L,K)+TRBALK(L,K)+TRDALK(L,K)+TRVALK(L,K)
      TAUAX=TAUAG+TAUAP
      IF(TAUAX.GT.1.D-06) GO TO 215
      DFLB(L,K)=0.D0
      ENA(L)=0.D0
      DNA(L)=0.D0
      TRA(L)=1.D0
      ENB(L)=0.D0
      DNB(L)=0.D0
      TRB(L)=1.D0
      ENC(L)=0.D0
      DNC(L)=0.D0
      TRC(L)=1.D0
      L=L-1
      IF(L.GT.0) GO TO 205
  210 CONTINUE
      L=L+1
      UFLB(L,K)=BG
      TRUFLB(L)=TRUFLB(L)+BG
      IF(L.LT.NLP) GO TO 210
      TOTLZF(1)=TOTLZF(1)+BG
      TOTLZF(2)=TOTLZF(2)+BG
      TOTLZF(3)=TOTLZF(3)+BG
      GO TO 200
  215 CONTINUE
  220 CONTINUE
      XTRDL=XTRD(L,IMOL+1)

      ITL=ITLT(L)+ITK0
      BTOP=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLT(L)
      ITL=ITLB(L)+ITK0
      BBOT=PLANCK(ITL)-(PLANCK(ITL)-PLANCK(ITL+1))*WTLB(L)
      TAUAG=TRGXLK(L,K)
      TAUAP=TRCALK(L,K)+TRAALK(L,K)+TRBALK(L,K)+TRDALK(L,K)+TRVALK(L,K)
      TAUAX=TAUAG+TAUAP

C               Optically thin limit emission/transmission approximation
C               --------------------------------------------------------

      IF(TAUAX.LT.1.D-04) THEN
      TAUBX=TAUAX+TAUAX
      TAUCX=10.D0*TAUAX
      BBAR=0.5D0*(BTOP+BBOT)
      TRA(L)=1.D0-TAUAX
      ENA(L)=BBAR*TAUAX
      DNA(L)=ENA(L)
      TX=TRA(L)*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNACUM=DNACUM*TX+DNA(L)
      TRB(L)=1.D0-TAUBX
      ENB(L)=BBAR*TAUBX
      DNB(L)=ENB(L)
      TX=TRB(L)*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNBCUM=DNBCUM*TX+DNB(L)
      TRC(L)=1.D0-TAUCX
      ENC(L)=BBAR*TAUCX
      DNC(L)=ENC(L)
      TX=TRC(L)*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNCCUM=DNCCUM*TX+DNC(L)
      GO TO 230
      ENDIF

C                     TAUB absorber-dependent extinction path adjustment
C                     --------------------------------------------------

      PLBN=PLB(L)
      TAUBG=TAUAG+TAUAG
      TAUCG=10.D0*TAUAG

       IF(IMOL.EQ.3.AND.PLBN.GT.500.D0) THEN
       IF(TAUAG.GT.0.05D0.AND.TAUAG.LT.0.25D0) THEN
      F=23.71D0*TAUAG**2-7.06D0*TAUAG+1.266D0
      TAUBG=TAUBG*F
      GO TO 221
       ENDIF
       ENDIF

      IF(TAUAG.GT.0.1) THEN
         IF(IMOL.EQ.1) THEN
         IF(PLBN.GT.250.D0) THEN
      F=0.75D0
      IF(TAUAG.LT.3.D0) F=0.92D0-0.053D0*TAUAG
      TAUBG=TAUBG*F
         ELSE
      F=0.70D0
      IF(TAUAG.LT.2.5D0) F=0.90D0-0.073D0*TAUAG
      TAUBG=TAUBG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.2) THEN
         IF(PLBN.GT.250.D0) THEN
      F=0.58D0
      IF(TAUAG.LT.3.5D0) F=0.93D0-0.097D0*TAUAG
      TAUBG=TAUBG*F
         ELSE
      F=0.70D0
      IF(TAUAG.LT.3.5D0) F=0.92D0-0.062D0*TAUAG
      TAUBG=TAUBG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.3) THEN
         IF(PLBN.GT.250.D0) THEN
      F=0.95D0
      IF(TAUAG.LT.0.5D0) F=0.99D0-0.016D0*TAUAG
      TAUBG=TAUBG*F
         ELSE
      F=0.75D0
      IF(TAUAG.LT.3.7D0) F=0.97D0-0.060D0*TAUAG
      TAUBG=TAUBG*F
      ENDIF
      ENDIF
      ENDIF

C                     TAUC absorber-dependent extinction path adjustment
C                     --------------------------------------------------
  221 CONTINUE

       IF(IMOL.EQ.3.AND.PLBN.GT.500.D0) THEN
       IF(TAUAG.GT.0.01D0.AND.TAUAG.LT.0.25D0) THEN
      F=26.14D0*TAUAG**2-6.93D0*TAUAG+1.0567D0
      TAUCG=TAUCG*F
      GO TO 222
       ENDIF
       ENDIF

      IF(TAUAG.GT.0.01D0) THEN
         IF(IMOL.EQ.1) THEN
         IF(PLBN.GT.250.D0) THEN
      F=0.65D0
      IF(TAUAG.LT.0.37D0) F=0.96D0-0.67D0*TAUAG
      TAUCG=TAUCG*F
         ELSE
      F=0.50D0
      IF(TAUAG.LT.0.47D0) F=0.87D0-0.71D0*TAUAG
      TAUCG=TAUCG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.2) THEN
         IF(PLBN.GT.250.D0) THEN
      F=0.50D0
      IF(TAUAG.LT.0.75D0) F=0.95D0-0.32D0*TAUAG
      TAUCG=TAUCG*F
         ELSE
      F=0.50D0
      IF(TAUAG.LT.0.70D0) F=0.90D0-0.59D0*TAUAG
      TAUCG=TAUCG*F
      ENDIF
      ENDIF
         IF(IMOL.EQ.3) THEN
         IF(PLBN.GT.250.D0) THEN
      F=0.95D0
      IF(TAUAG.LT.0.5D0) F=0.98D0-0.039D0*TAUAG
      TAUCG=TAUCG*F
         ELSE
      F=0.75
      IF(TAUAG.LT.0.70D0) F=0.98D0-0.29D0*TAUAG
      TAUCG=TAUCG*F
      ENDIF
      ENDIF
      ENDIF
  222 CONTINUE

      TAUBP=TAUAP+TAUAP
      TAUCP=10.D0*TAUAP
      TAUA=TAUAG+TAUAP
      TAUB=TAUBG+TAUBP
      TAUC=TAUCG+TAUCP

      IF(L.EQ.LTOPCL.AND.KCLDEM.EQ.1) GO TO 225

      BDIF=BBOT-BTOP
      BBTA=BDIF/TAUA
      BBTB=BDIF/TAUB
      BBTC=BDIF/TAUC

C            Optically thick limit non-scattering emission approximation
C            -----------------------------------------------------------

      IF(TAUA.GT.9.D0) THEN
      TRA(L)=0.D0
      TRB(L)=0.D0
      TRC(L)=0.D0
      ENA(L)=BTOP+BBTA
      ENB(L)=BTOP+BBTB
      ENC(L)=BTOP+BBTC
      DNA(L)=BBOT-BBTA
      DNB(L)=BBOT-BBTB
      DNC(L)=BBOT-BBTC
      DNACUM=BBOT-BBTA
      DNBCUM=BBOT-BBTB
      DNCCUM=BBOT-BBTC
      GO TO 230
      ENDIF

      IF(TAUA.LT.0.5D0) THEN
      TAUA2=TAUA*TAUA
      TRANA=1.D0-TAUA+(0.5D0-R6*TAUA+R24*TAUA2)*TAUA2
      ELSE
      TRANA=EXP(-TAUA)
      ENDIF
      IF(TAUB.LT.0.5D0) THEN
      TAUB2=TAUB*TAUB
      TRANB=1.D0-TAUB+(0.5D0-R6*TAUB+R24*TAUB2)*TAUB2
      ELSE
      TRANB=EXP(-TAUB)
      ENDIF
      IF(TAUC.LT.0.5D0) THEN
      TAUC2=TAUC*TAUC
      TRANC=1.D0-TAUC+(0.5D0-R6*TAUC+R24*TAUC2)*TAUC2
      ELSE
      TRANC=EXP(-TAUC)
      ENDIF

      TRA(L)=TRANA
      ENA(L)=BTOP+BBTA-(BBOT+BBTA)*TRANA
      DNA(L)=BBOT-BBTA-(BTOP-BBTA)*TRANA
      TX=TRANA*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNACUM=DNACUM*TX+DNA(L)
      TRB(L)=TRANB
      ENB(L)=BTOP+BBTB-(BBOT+BBTB)*TRANB
      DNB(L)=BBOT-BBTB-(BTOP-BBTB)*TRANB
      TX=TRANB*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNBCUM=DNBCUM*TX+DNB(L)
      TRC(L)=TRANC
      ENC(L)=BTOP+BBTC-(BBOT+BBTC)*TRANC
      DNC(L)=BBOT-BBTC-(BTOP-BBTC)*TRANC
      TX=TRANC*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNCCUM=DNCCUM*TX+DNC(L)
      GO TO 230

C                          ---------------------------------------------
C                          Top-cloud multiple scattering corrections for
C                          emitted, transmitted, and reflected radiances
C                          and fluxes at the top-cloud (L=LTOPCL) level.
C                          ---------------------------------------------

  225 CONTINUE
      TRA(L)=EXP(-TAUAG-TAUAP*FDXTCK(3,K))
      TRB(L)=EXP(-TAUBG-TAUBP*FDXTCK(2,K))
      TRC(L)=EXP(-TAUCG-TAUCP*FDXTCK(1,K))
      DEC=C*DNCCUM*RIJTCK(1,K)+B*DNBCUM*RIJTCK(2,K)+A*DNACUM*RIJTCK(3,K)
      DEB=C*DNCCUM*RIJTCK(2,K)+B*DNBCUM*RIJTCK(4,K)+A*DNACUM*RIJTCK(5,K)
      DEA=C*DNCCUM*RIJTCK(3,K)+B*DNBCUM*RIJTCK(5,K)+A*DNACUM*RIJTCK(6,K)
      ALBTCK(1,K)=C*RIJTCK(1,K)+B*RIJTCK(2,K)+A*RIJTCK(3,K)
      ALBTCK(2,K)=C*RIJTCK(2,K)+B*RIJTCK(4,K)+A*RIJTCK(5,K)
      ALBTCK(3,K)=C*RIJTCK(3,K)+B*RIJTCK(5,K)+A*RIJTCK(6,K)
      COALB1=1.D0-ALBTCK(1,K)
      COALB2=1.D0-ALBTCK(2,K)
      COALB3=1.D0-ALBTCK(3,K)
      TAUA=TAUAG+TAUAP*FEMTCK(3,K)
      TAUB=TAUBG+TAUBP*FEMTCK(2,K)
      TAUC=TAUCG+TAUCP*FEMTCK(1,K)
      TRANA=EXP(-TAUA)
      TRANB=EXP(-TAUB)
      TRANC=EXP(-TAUC)
      BDIF=BBOT-BTOP
      BBTA=BDIF/TAUA
      BBTB=BDIF/TAUB
      BBTC=BDIF/TAUC
      ENA(L)=(BTOP+BBTA-(BBOT+BBTA)*TRANA)*COALB3
      DNA(L)=(BBOT-BBTA-(BTOP-BBTA)*TRANA)*COALB3
      TX=TRA(L)*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNACUM=DNACUM*TX+DNA(L)
      ENB(L)=(BTOP+BBTB-(BBOT+BBTB)*TRANB)*COALB2
      DNB(L)=(BBOT-BBTB-(BTOP-BBTB)*TRANB)*COALB2
      TX=TRB(L)*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNBCUM=DNBCUM*TX+DNB(L)
      ENC(L)=(BTOP+BBTC-(BBOT+BBTC)*TRANC)*COALB1
      DNC(L)=(BBOT-BBTC-(BTOP-BBTC)*TRANC)*COALB1
      TX=TRC(L)*XTRDL
      IF(TX.GT.1.D0) TX=1.D0
      DNCCUM=DNCCUM*TX+DNC(L)
      ENC(L)=ENC(L)+DEC
      ENB(L)=ENB(L)+DEB
      ENA(L)=ENA(L)+DEA
  230 CONTINUE
      FDNABC=A*DNACUM+B*DNBCUM+C*DNCCUM
      TRDFLB(L)=TRDFLB(L)+FDNABC
      DFLB(L,K)=FDNABC
      L=L-1
      IF(L.GT.0) GO TO 220

C             Old form of scattering correction is skipped when KCLDEM=1
C             ----------------------------------------------------------

      IF(KCLDEM.EQ.1) GO TO 235
      IF(LTOPCL.LT.1) GO TO 235
      ENA(LTOPCL)=ENA(LTOPCL)*(1.0-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
      ENB(LTOPCL)=ENB(LTOPCL)*(1.0-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
      ENC(LTOPCL)=ENC(LTOPCL)*(1.0-TRCTCA(K))+TRCTCA(K)*DFLB(LTOPCL+1,K)
  235 CONTINUE

!sl   ------------------------------------------------------------------
!sl                                       SURFACE LAYER FLUX COMPUTATION
!sl   with TAUSL,FTAUSL=0 defaults, surface layer calculation is skipped
!sl   ------------------------------------------------------------------

      L=1
      DFSL(K)=FDNABC
!sl   TAUA=TAUSL(K)+FTAUSL(K)
!sl   IF(TAUA.GT.1.D-06) GO TO 240
      BG=BG+FDNABC*TRGALB(K)
      UNA=BG
      UNB=BG
      UNC=BG
      FUNABC=BG
!sl   GO TO 245
  240 CONTINUE
!sl   ITS=TSL
!sl   WTS=TSL-ITS
!sl   WTS1=1.D0-WTS
!sl   ITS=ITS-ITPFT0
!sl   ITS=ITS+(K-1)*ITNEXT
!sl   BS=PLANCK(ITS)*WTS1+PLANCK(ITS+1)*WTS
!sl   TA=EXP(-TAUA)
!sl   TB=TA*TA
!sl   TC=(TB*TB*TA)**2
!sl   DNA(1)=(DNA(1)-BS)*TA+BS
!sl   DNB(1)=(DNB(1)-BS)*TB+BS
!sl   DNC(1)=(DNC(1)-BS)*TC+BS
!sl   FDNABC=A*DNA(1)+B*DNB(1)+C*DNC(1)
!sl   BG=BGFEMT(K)+FDNABC*TRGALB(K)
!sl   UNA=(BG-BS)*TA+BS
!sl   UNB=(BG-BS)*TB+BS
!sl   UNC=(BG-BS)*TC+BS
!sl   FUNABC=A*UNA+B*UNB+C*UNC
!sl   BSP=PLANCK(ITS+1)*WTS1+PLANCK(ITS+2)*WTS
!sl   BSM=PLANCK(ITS-1)*WTS1+PLANCK(ITS  )*WTS
!sl   SLABS=1.D0-A*TA-B*TB-C*TC
!sl   TRSLTS=TRSLTS+(BSP-BSM)*SLABS
!sl   TRSLTG=TRSLTG+BGFEMD(K)*SLABS
!sl   TRSLBS=TRSLBS+BS*SLABS

C     ------------------------------------------------------------------
C                                                UPWARD FLUX COMPUTATION
C     ------------------------------------------------------------------

  245 CONTINUE
      TRUFLB(L)=TRUFLB(L)+FUNABC
      UFLB(L,K)=FUNABC

C       ----------------------------------------------------------------
C       At top-cloud level, compute component of upwelling flux relected
C       downward by cloud bottom and add to downwelling flux below cloud
C       ----------------------------------------------------------------

      IF(L.EQ.LTOPCL.AND.KCLDEM.EQ.1) THEN
      DEC=C*UNC*RIJTCK(1,K)+B*UNB*RIJTCK(2,K)+A*UNA*RIJTCK(3,K)
      DEB=C*UNC*RIJTCK(2,K)+B*UNB*RIJTCK(4,K)+A*UNA*RIJTCK(5,K)
      DEA=C*UNC*RIJTCK(3,K)+B*UNB*RIJTCK(5,K)+A*UNA*RIJTCK(6,K)
      N=L
  250 CONTINUE
      DNA(N)=DNA(N)+DEA
      DNB(N)=DNB(N)+DEB
      DNC(N)=DNC(N)+DEC
      DDFLUX=A*DEA+B*DEB+C*DEC
      TRDFLB(N)=TRDFLB(N)+DDFLUX
      DFLB(N,K)=DFLB(N,K)+DDFLUX
      IF(N.LT.2) GO TO 255
      N=N-1
      DEA=DEA*TRA(N)
      DEB=DEB*TRB(N)
      DEC=DEC*TRC(N)
      GO TO 250
  255 CONTINUE
      ENDIF
      XTRUL=XTRU(L,IMOL+1)
      TX=TRA(L)*XTRUL
      IF(TX.GT.1.D0) TX=1.D0
      UNA=UNA*TX+ENA(L)
      TX=TRB(L)*XTRUL
      IF(TX.GT.1.D0) TX=1.D0
      UNB=UNB*TX+ENB(L)
      TX=TRC(L)*XTRUL
      IF(TX.GT.1.D0) TX=1.D0
      UNC=UNC*TX+ENC(L)
      FUNABC=A*UNA+B*UNB+C*UNC

      L=L+1
      IF(L.LT.NLP) GO TO 245
      IF(K.EQ.1) THEN
      TRUFTW=FUNABC
      TRDFGW=TRDFLB(1)
      TRUFGW=BG
      WINDZF(1)=UNA
      WINDZF(2)=UNB
      WINDZF(3)=UNC
      ENDIF

      TRUFLB(NLP)=TRUFLB(NLP)+FUNABC
      UFLB(NLP,K)=FUNABC
      UFSL(K)=UFLB(1,K)
      TOTLZF(1)=TOTLZF(1)+UNA
      TOTLZF(2)=TOTLZF(2)+UNB
      TOTLZF(3)=TOTLZF(3)+UNC

      GO TO 200
  300 CONTINUE

      DO 310 L=1,NLP
      TRNFLB(L)=TRUFLB(L)-TRDFLB(L)
  310 CONTINUE
      DO 320 L=1,NL
      TRFCRL(L)=TRNFLB(L+1)-TRNFLB(L)
  320 CONTINUE

C**** Window region and spectr. integrated total flux diagnostics
      PFW=max( 1.0001d-2, 10.D0*TRUFTW )
      IF(PFW.GT.719.999d0) PFW=719.999d0
      IPF=PFW
      DPF=0.  ! initialise?
      IF(IPF.LT.10) GO TO 330
      DPF=PFW-IPF
      IPF=IPF+180
      GO TO 350
  330 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
      IF(IPF.LT.10) GO TO 340
      DPF=PFW-IPF
      IPF=IPF+90
      GO TO 350
  340 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
C     IF(IPF.LT.1) IPF=1
  350 CONTINUE
      BTEMPW=TKPFW(IPF)+DPF*(TKPFW(IPF+1)-TKPFW(IPF))

      DO 390 II=1,3
      PFW=max( 1.d0, TOTLZF(II) )
      IF(PFW.GT.899.999d0) PFW=899.999d0
      IPF=PFW
      DPF=PFW-IPF
C     IF(IPF.LT.1) IPF=1
C     IF(IPF.GT.899) IPF=899
      TOTLZT(II)=TKPFT(IPF)+DPF*(TKPFT(IPF+1)-TKPFT(IPF))

      PFW=max( 1.0001d-2, 10.D0*WINDZF(II) )
      IF(PFW.GT.719.999d0) PFW=719.999d0
      IPF=PFW
      DPF=0.    ! initialise?
      IF(IPF.LT.10) GO TO 360
      DPF=PFW-IPF
      IPF=IPF+180
      GO TO 380
  360 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
      IF(IPF.LT.10) GO TO 370
      DPF=PFW-IPF
      IPF=IPF+90
      GO TO 380
  370 CONTINUE
      PFW=10.D0*PFW
      IPF=PFW
  380 CONTINUE
      WINDZT(II)=TKPFW(IPF)+DPF*(TKPFW(IPF+1)-TKPFW(IPF))
  390 CONTINUE
      RETURN
      END SUBROUTINE THERML

      SUBROUTINE SOLAR0
      IMPLICIT NONE

      INTEGER, PARAMETER, DIMENSION(17) :: NMKWAV = (/ 200, 360, 770,
     *     795, 805, 810, 860,1250,1500,1740,2200,3000,3400,3600,3800
     *     ,4000,9999/)
      INTEGER, PARAMETER, DIMENSION(16) :: LORDER = (/15,14, 8, 7, 6, 5,
     *     13, 12, 4, 3, 2, 1, 11, 10, 9,16/)
      INTEGER N2,I

      N2=1
      DO 1 I=1,30
      N2=N2*2
      DBLN(I)=N2
    1 CONTINUE
      DO 2 I=1,16
      NORDER(I)=LORDER(I)
      NMWAVA(I)=NMKWAV(I)
      NMWAVB(I)=NMKWAV(I+1)
    2 CONTINUE

      TCLMIN=MIN(TAUIC0,TAUWC0)

      RETURN
      END SUBROUTINE SOLAR0

      SUBROUTINE SOLARM
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SOLARM Returns:
C                      SRDFLB   Solar downward flux at layer bottom edge
C                      SRUFLB   Solar  upward  flux at layer bottom edge
C                      SRNFLB   Solar net downward flux  (in Watts/m**2)
C                      SRFHRL   Solar heating rate/layer (in Watts/m**2)
C                      FSRNFG   Solar flux abs at ground by surface-type
C                               (see explanatory note at end of SOLARM)
C               Also:
C                      TOA:   SRIVIS SROVIS PLAVIS  SRINIR SRONIR PLANIR
C                      BOA:   SRDVIS SRUVIS ALBVIS  SRDNIR SRUNIR ALBNIR
C                      ATM:   SRTVIS SRRVIS SRAVIS  SRTNIR SRRNIR SRANIR
C                             SRXVIS SRXNIR (Direct beam only at ground)
C
C           Spectral:  (by k-distribution/pseudo-spectral) breakdown:
C                      SKDFLB   Solar downward flux at layer bottom edge
C                      SKUFLB   Solar  upward  flux at layer bottom edge
C                      SKNFLB   Solar net downward flux  (in Watts/m**2)
C                      SKFHRL   Solar heating rate/layer (in Watts/m**2)
C
C                      SRKALB   Planetary albedo (by spectral breakdown)
C                      SRKINC   Incident fluxedo (by spectral breakdown)
C                      SRKGAX   Direct  k-d flux absorbed by ground-type
C                      SRKGAD   Diffuse k-d flux absorbed by ground-type
C     ------------------------------------------------------------------
C     Remarks:
C              NORMS0=1 Incident (TOA) Solar flux normalized to equal S0
C                       (COSZ dependence included in calculated results)
C                        The returned solar fluxes have to be multiplied
C                       by COSZ to yield actual atmospheric heating rate
C
C              NMKWAV   Spectral/k-distribution subdivisions are nominal
C                       (due to spectral trading of absorption features)
C
C                 VIS   Designates solar visible wavelengths (   <770nm)
C                 NIR   Designates solar near-IR wavelengths (770>   nm)
C                       VIS comprises .53 of S0, NIR comprises .47 of S0
C     ------------------------------------------------------------------
C
C     ------------------------------------------------------------------
C         Fractional solar flux k-distribution/pseudo-spectral intervals
C
C         KSLAM=    1     1     2     2     5     5     5     5
C             K=    1     2     3     4     5     6     7     8
C     DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
C         KSLAM=    1     1     1     3     4     6     6     1
C             K=    9    10    11    12    13    14    15    16
C    +           .002, .003, .003, .072, .200, .480, .050, .011/
C
C     ------------------------------------------------------------------
C     The nominal spectral order for k-dist/pseudo-spectral intervals is
C     (WavA and WavB designate approximate spectral interval boundaries)
C
C             L=   12    11    10     9     6     5     4     3
C     WavA (nm)= 3000  2200  1740  1500   810   805   795   770
C     WavB (nm)= 3400  3000  2200  1740   860   810   805   795
C             K=    1     2     3     4     5     6     7     8
C     DATA DKS0/ .010, .030, .040, .040, .040, .002, .004, .013,
C
C             L=   15    14    13     8     7     2     1    16
C     WavA (nm)= 3800  3500  3400  1250   860   360   200  4000
C     WavB (nm)= 4000  3800  3600  1500  1250   770   360  9999
C             K=    9    10    11    12    13    14    15    16
C    +           .002, .003, .003, .072, .200, .480, .050, .011/
C
C     ------------------------------------------------------------------
C     6 spectral intervals overlap the 16 solar k-distribution intervals
C
C     Cloud and aerosol Mie scattering parameters (also surface albedos)
C     are averaged over these spectral intervals. These intervals are in
C     reverse spectral order. Thus spectral interval 6 refers to visible
C     (VIS) wavelengths, intervals 1-5 refer to nearIR (NIR) wavelengths
C     KSLAM designates the spectral interval of first 14 k-distributions
C     (K=15 for UV ozone absorption refers to (VIS) spectral interval 6)
C     (K=16 represents strong absorbing spectral regions via interval 1)
C
C     The nominal Mie scattering spectral band subdivisions are:
C
C                         -------------NIR------------      VIS
C                     L=    1     2     3     4     5        6
C             WavA (nm)=  2200  1500  1250   860   770      300
C             WavB (nm)=  4000  2200  1500  1250   860      770
C
C     ------------------------------------------------------------------

      REAL*8 COLEXT(6),COLSCT(6),COLGCB(6)   ! ,ALLGCB(6)

C                            -------------------------------------------
C                            NO2, O3 Chappuis Band, Rayleigh, parameters
C                            -------------------------------------------
      REAL*8, PARAMETER :: XCMNO2=5.465d0, XCMO3=.0399623d0,
     *     TOTRAY=0.000155d0
      REAL*8 S0COSZ,COSMAG,SECZ,TAURAY,RTAU,SUMEXT,SUMCST,SUMCGB,COLPFG
     *     ,SURFBB,TAUSBB,ALLTAU,TAULAY,GCBLAY,RTAUL,DKS0X,RBNB,RBNX
     *     ,RCNB,RCNX,TLN,PLN,ULN,TERMA,TERMB,TAU1,TAU,PIZERO,PR,PT,DBLS
     *     ,XANB,XANX,TANB,TANX,XXT,RASB,RASX,BNORM,XNORM,RARB,RARX,XATB
     *     ,DENOM,DB,DX,UB,UX,RBXTOA,ATOPX,ATOPD,O3CMX,O3CMD
     *     ,SUMSCT,SUMGCB,XXG,SURX,PFF,XATC,XBNB,XBNX,TBNB,TBNX,XBTB
     *     ,ABOTX,ABOTD,AO3UXN,AO3UDN,SRKA16,DKS0XX,TRNC,CLX,TRNU,TRN1
     *     ,TRN2,TRN3,TAUG,TAU2,TAU3,S0VIS,S0NIR,SUMX,SUMD,SUMU,SUMN
     *     ,SUMH,SGPG
      INTEGER I,K,KK,L,M,N,NN,KLAM,NDBLS

      S0COSZ=S0
      IF(NORMS0.EQ.0) S0COSZ=S0*COSZ

      DO 10 N=1,NLP
      SRNFLB(N)=0.D0
      SRDFLB(N)=0.D0
      SRUFLB(N)=0.D0
      SRFHRL(N)=0.D0
      SKDFLB(N,16)=0.D0
      SKUFLB(N,16)=0.D0
   10 CONTINUE
      DO N=1,16
        SRKALB(N)=0.D0          ! for WRITER only
      END DO
      dblext=0. ; dblsct=0. ; dblgcb=0. ; dblpi0=0. ! for writer only
      skdflb=0. ; sknflb=0. ; skuflb=0.             ! for writer only
      skfhrl=0. ; srkgax=0. ; srkgad=0.             ! for writer only
C                     TOA solar flux VIS/NIR subdivision
C                     (incident, outgoing, plane albedo)
C                     ----------------------------------
      SRIVIS=0.D0
      SROVIS=0.D0
      PLAVIS=1.D0
      SRINIR=0.D0
      SRONIR=0.D0
      PLANIR=1.D0
C                     BOA solar flux VIS/NIR subdivision
C                     (incident, upward, surface albedo)
C                     ----------------------------------
      SRDVIS=0.D0
      SRUVIS=0.D0
      ALBVIS=1.D0
      SRDNIR=0.D0
      SRUNIR=0.D0
      ALBNIR=1.D0
C                     Fractional atmos only flux VIS/NIR subdivision
C                     (fractions reflected, transmitted, & absorbed)
C                     ----------------------------------------------
      SRRVIS=1.D0
      SRTVIS=0.D0
      SRAVIS=0.D0
      SRRNIR=0.D0
      SRTNIR=0.D0
      SRANIR=0.D0
C                     Direct beam, fractional S0 VIS/NIR subdivision
C                     ----------------------------------------------
      SRXVIS=0.D0
      SRXNIR=0.D0
C                     Ground surface absorbed solar flux subdivision
C                     according to 4 fractional surface-type albedos
C                     ----------------------------------------------
      FSRNFG(1)=0.D0
      FSRNFG(2)=0.D0
      FSRNFG(3)=0.D0
      FSRNFG(4)=0.D0

      IF(COSZ.LT.0.001D0) RETURN
      COSMAG=35.D0/SQRT(1224.D0*COSZ*COSZ+1.D0)
      SECZ=1.D0/COSZ

      TAURAY=TOTRAY*FRAYLE

      DO 90 K=1,6
      RTAU=1.D-10
      IF(K.EQ.6) RTAU=TAURAY
      COLEXT(K)=0.D0
      COLSCT(K)=0.D0
      COLGCB(K)=0.D0
      DO 30 L=1,NL
      RTAUL=RTAU*(PLB(L)-PLB(L+1))
      SUMEXT=RTAUL+SRCEXT(L,K)+SRAEXT(L,K)+SRBEXT(L,K)
     +                        +SRDEXT(L,K)+SRVEXT(L,K)
      SUMSCT=RTAUL+SRCSCT(L,K)+SRASCT(L,K)+SRBSCT(L,K)
     +                        +SRDSCT(L,K)+SRVSCT(L,K)
      SUMGCB=      SRCSCT(L,K)*SRCGCB(L,K)+SRASCT(L,K)*SRAGCB(L,K)
     +            +SRBSCT(L,K)*SRBGCB(L,K)+SRDSCT(L,K)*SRDGCB(L,K)
     +            +SRVSCT(L,K)*SRVGCB(L,K)
      DBLEXT(L,K)=SUMEXT
      DBLSCT(L,K)=SUMSCT
      DBLGCB(L,K)=SUMGCB/(SUMSCT+1.D-10)
      DBLPI0(L,K)=SUMSCT/(SUMEXT+1.D-10)
      COLEXT(K)=COLEXT(K)+DBLEXT(L,K)
      COLSCT(K)=COLSCT(K)+DBLSCT(L,K)
      COLGCB(K)=COLGCB(K)+DBLSCT(L,K)*DBLGCB(L,K)
   30 CONTINUE
      COLGCB(K)=COLGCB(K)/(COLSCT(K)+1.D-10)

      IF(KANORM.GT.0) THEN
C     -----------------------------------------------------------------
C     KANORM (default = 0)  Option to renormalize aerosol column albedo
C                           to make column albedo less dependent on the
C                           number of model layers due to SGP treatment
C
C                           KANORM=1  aerosol column only is normalized
C
C                           KANORM=2  aerosol plus ground is normalized
C                                     with Tau equivalent ground albedo
C                                     ---------------------------------
      COLPFG=COLGCB(K)
      SURFBB=SRBALB(K)
      TAUSBB=0.D0
      IF(KANORM.GT.1) CALL GTSALB(XXG,XXT,SURX,SURFBB,COLPFG,TAUSBB,2)
      DBLEXT(NLP,K)=TAUSBB
      ALLTAU=TAUSBB+COLEXT(K)
      CALL SGPGXG(COSZ,ALLTAU,COLPFG,SGPG)
cc    ALLGCB(K)=SGPG
      DO 40 L=1,NL
      DBLGCB(L,K)=SGPG
   40 CONTINUE
      ELSE

      DO 50 L=1,NL
      TAULAY=DBLEXT(L,K)
      GCBLAY=DBLGCB(L,K)
      CALL SGPGXG(COSZ,TAULAY,GCBLAY,SGPG)
      DBLGCB(L,K)=SGPG
   50 CONTINUE
      ENDIF

      IF(LTOPCL.LT.1) GO TO 90
      RTAU=1.D-10
      IF(K.EQ.6) RTAU=TAURAY
      COLEXT(K)=0.D0
      COLSCT(K)=0.D0
      COLGCB(K)=0.D0
      DO 60 L=1,NL
      IF(SRCEXT(L,K).LT.TCLMIN) GO TO 60
      RTAUL=RTAU*(PLB(L)-PLB(L+1))
      SUMEXT=RTAUL+SRCEXT(L,K)+SRAEXT(L,K)+SRBEXT(L,K)
     +                        +SRDEXT(L,K)+SRVEXT(L,K)
      SUMSCT=RTAUL+SRCSCT(L,K)+SRASCT(L,K)+SRBSCT(L,K)
     +                        +SRDSCT(L,K)+SRVSCT(L,K)
      SUMGCB=      SRCSCT(L,K)*SRCGCB(L,K)+SRASCT(L,K)*SRAGCB(L,K)
     +            +SRBSCT(L,K)*SRBGCB(L,K)+SRDSCT(L,K)*SRDGCB(L,K)
     +            +SRVSCT(L,K)*SRVGCB(L,K)
      DBLEXT(L,K)=SUMEXT
      DBLSCT(L,K)=SUMSCT
      DBLGCB(L,K)=SUMGCB/(SUMSCT+1.D-10)
      DBLPI0(L,K)=SUMSCT/(SUMEXT+1.D-10)
      COLEXT(K)=COLEXT(K)+DBLEXT(L,K)
      COLSCT(K)=COLSCT(K)+DBLSCT(L,K)
      COLGCB(K)=COLGCB(K)+DBLSCT(L,K)*DBLGCB(L,K)
   60 CONTINUE
      COLGCB(K)=COLGCB(K)/(COLSCT(K)+1.D-10)

      IF(KCNORM.GT.0) THEN
C     -----------------------------------------------------------------
C     KCNORM (default = 0)  Option to renormalize  cloud  column albedo
C                           to make column albedo less dependent on the
C                           number of model layers due to SGP treatment
C
C                           KCNORM=1    cloud column only is normalized
C
C                           KCNORM=2    cloud plus ground is normalized
C                                     with Tau equivalent ground albedo
C                                     ---------------------------------
      COLPFG=COLGCB(K)
      SURFBB=SRBALB(K)
      TAUSBB=0.D0
      IF(KCNORM.GT.1) CALL GTSALB(XXG,XXT,SURX,SURFBB,COLPFG,TAUSBB,2)
      DBLEXT(NLP,K)=TAUSBB
      ALLTAU=TAUSBB+COLEXT(K)
      CALL SGPGXG(COSZ,ALLTAU,COLPFG,SGPG)
cc    ALLGCB(K)=SGPG
      DO 70 L=1,NL
      IF(SRCEXT(L,K).LT.TCLMIN) GO TO 70
      DBLGCB(L,K)=SGPG
   70 CONTINUE
      ELSE
      DO 80 L=1,NL
      IF(SRCEXT(L,K).LT.TCLMIN) GO TO 80
      TAULAY=DBLEXT(L,K)
      GCBLAY=DBLGCB(L,K)
      CALL SGPGXG(COSZ,TAULAY,GCBLAY,SGPG)
      DBLGCB(L,K)=SGPG
   80 CONTINUE
      ENDIF
   90 CONTINUE

      K = 0
  300 CONTINUE
      K = K+1

      KLAM=KSLAM(K)
      DKS0X=DKS0(K)*S0COSZ
      RBNB=SRBALB(KLAM)
      RBNX=SRXALB(KLAM)
      RCNB=0.D0
      RCNX=0.D0
      SRKINC(K)=DKS0X

      N = 0
  200 CONTINUE
      N = N+1

      SRB(N)=RBNB
      SRX(N)=RBNX
      TLN=TLM(N)
      PLN=PL(N)
      ULN=ULGAS(N,1)

C     Select parameterized k-distribution gas absorption by H2O, O2, CO2
C     ------------------------------------------------------------------

      GO TO (101,102,103,104,105,106,107,108,109,110,111,112,113,114),K
 101  CONTINUE
C--------K=6-------H2O       DS0=.01
      TERMA=(35.66+TLN*(.0416-.0004622*TLN+.001057*PLN))*(1.+.04286*PLN)
      TERMB=(1.+.00171*ULN)*(1.+PLN*(189.088+.1316*PLN))
      TAU1 =TERMA/TERMB
      IF(TAU1.GT.0.02343) TAU1=0.02343
      TAU=TAU1*ULN
      GO TO 120
 102  CONTINUE
C--------K=5-------H2O       DS0=.03
      TERMA=(2.792+TLN*(.0914-.0002848*TLN+.0003395*PLN))
     +     *(1.+.02964*PLN)
      TERMB=(1.0+.000657*ULN)*(1.+PLN*(240.70+.13847*PLN))
      TAU1 =TERMA/TERMB
      IF(TAU1.GT.0.00520) TAU1=0.00520
      TAU=TAU1*ULN
      GO TO 120
 103  CONTINUE
C--------K=4-------H2O       DS0=.04
      TERMA=(.4768+.467E-04*PLN*TLN)*(1.+TLN*(.00191-.719E-05*TLN))
      TERMB=(1.+.717E-04*ULN)*(1.+PLN*(130.56+.0876*PLN))/(1.+.0266*PLN)
      TAU1 =TERMA/TERMB
      IF(TAU1.GT.0.00150) TAU1=0.0015
      TAU=TAU1*ULN
      GO TO 120
 104  CONTINUE
C--------K=3-------H2O       DS0=.04
      TERMA=(.000247*TLN-.091+PLN*(.00035+.78E-06*TLN))*(1.+.2847*PLN)
      TERMB=(1.+.2066E-04*ULN)*(1.+PLN*(137.17+.16132*PLN))
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 105  CONTINUE
C--------K=2-------H2O       DS0=.04
      TERMA=(PLN*(1.974/TLN+.0001117*TLN)-10.713)*(1.+.005788*TLN)
     +     *(1.+.001517*PLN)
      TERMB=(1.+.3218E-04*ULN)*(1.+PLN*(863.44+.2048*PLN))
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 106  CONTINUE
C--------K=4-------O2        DS0=.002
      ULN=ULGAS(N,4)
      TERMA=(.2236E-05-.1181E-09*TLN)*(1.+PLN*(.6364E-05*PLN+.001168))
      TERMB=1.+.1521E-05*ULN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 107  CONTINUE
C--------K=3-------O2        DS0=.004
      ULN=ULGAS(N,4)
      TERMA=(.3179E-06-.9263E-11*TLN)*(1.+PLN*(.8832E-05*PLN+.0005292))
      TERMB=1.+.1968E-06*ULN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 108  CONTINUE
C--------K=2-------O2        DS0=.013
      ULN=ULGAS(N,4)
      TERMA=(.2801E-07-.1638E-12*TLN)*(1.+PLN*(.1683E-04*PLN-.001721))
      TERMB=1.+.8097E-07*ULN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 109  CONTINUE
C--------K=4-------CO2       DS0=.002
      ULN=ULGAS(N,2)
      TERMA=(50.73-.03155*TLN-PLN*(.5543+.00091*TLN))*(1.-.1004*PLN)
      TERMB=(1.+.006468*ULN)*(1.+PLN*(49.51+.8285*PLN))
      TAU  =(TERMA/TERMB)*ULN
      IF(PLN.LT.175.0) TAU=(.00018*PLN+0.00001)*ULN
      GO TO 120
 110  CONTINUE
C--------K=3-------CO2       DS0=.003
      ULN=ULGAS(N,2)
      TERMA=(1.+.01319*TLN)*(PLN*(.008001*ULN+.4589E-03)-.8396*ULN)
      TERMB=ULN*(PLN+295.7+1.967*ULN)+.15126*PLN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 111  CONTINUE
C--------K=2-------CO2       DS0=.003
      ULN=ULGAS(N,2)
      TERMA=(1.+.02257*TLN)*(PLN*(.002295*ULN-.5489E-04)-.7571*ULN)
      TERMB=ULN*(PLN+803.9+2.477*ULN)-.09899*PLN
      TAU  =(TERMA/TERMB)*ULN
      GO TO 120
 112  CONTINUE
      TAU=0.D0
      GO TO 120
 113  CONTINUE
      TAU=0.D0
      GO TO 120
 114  CONTINUE
      TAU=XCMNO2*ULGAS(N,5)+XCMO3*ULGAS(N,3)
 120  CONTINUE

C     With 10 doublings to get to Tau=1.0, maximum seed tau is < 1/1024.
C     ------------------------------------------------------------------

      IF(TAU.LT.0.D0) TAU=0.D0

      TAU=TAU+DBLEXT(N,KLAM)
      IF(TAU.LT.1.D-06) GO TO 180
      PIZERO=DBLSCT(N,KLAM)/TAU
      IF(PIZERO.LT.0.001D0) GO TO 180

      PFF=DBLGCB(N,KLAM)

      NDBLS=0
      PR=1.D0-PFF
      PT=1.D0+PFF
      IF(TAU.GT.0.0019531D0) THEN
      DBLS=10.D0+1.44269D0*LOG(TAU)
      NDBLS=DBLS
      TAU=TAU/DBLN(NDBLS)
      ENDIF

C     Set optically thin limit values of R,T,X using PI0 renormalization
C     ------------------------------------------------------------------

      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      TANB=PT*XANB
      XXT=(SECZ-2.D0)*TAU
      TANX=PT*SECZ
     +    *(.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))*XANX
      RASB=PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT=(SECZ+2.D0)*TAU
      RASX=PR*SECZ
     +    *(.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM=(1.D0-XANB)/(RASB+TANB)*PIZERO
      XNORM=(1.D0-XANX)/(RASX+TANX)*PIZERO
      RASB=RASB*BNORM
      RASX=RASX*XNORM
      TANB=TANB*BNORM
      TANX=TANX*XNORM

C     Compute and record R,T,X atmospheric layer doubling/adding results
C     ------------------------------------------------------------------

      IF(NDBLS.LT.1) GO TO 170
      DO 160 NN=1,NDBLS
      RARB=RASB*RASB
      RARX=XANX*RASX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+RARX*RASB)/DENOM
      UB=RASB*(XANB+DB)
      UX=RARX+RASB*DX
      RASB=RASB+XATB*UB
      RASX=RASX+XATB*UX
      TANB=XANB*TANB+XATB*DB
      TANX=XANX*TANX+XATB*DX
      XANB=XANB*XANB
      XANX=XANX*XANX
  160 CONTINUE
  170 CONTINUE
      RARB=RASB*RBNB
      RARX=RASB*RBNX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+XANX*RARX)/DENOM
      UB=RBNB*(XANB+DB)
      UX=RBNX*XANX+RBNB*DX
      RBNB=RASB+XATB*UB
      RBNX=RASX+XATB*UX
      XATC=XATB/(1.D0-RASB*RCNB)
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC
      RCNB=RASB+RCNB*XATB*XATC
      GO TO 190
  180 CONTINUE
      RASB=0.D0
      RASX=0.D0
      TANB=0.D0
      TANX=0.D0
      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      DX=0.D0
      UX=RBNX*XANX
      RBNB=RBNB*XANB*XANB
      RBNX=UX*XANB
      RCNB=RCNB*XANB*XANB
      RCNX=RCNX*XANX*XANB
  190 CONTINUE
      RNB(N)=RASB
      RNX(N)=RASX
      TNB(N)=TANB
      TNX(N)=TANX
      XNB(N)=XANB
      XNX(N)=XANX
      IF(N.LT.NL) GO TO 200

C     Record fluxes, spectral components at TOA, & top-layer bottom edge
C     ------------------------------------------------------------------

      SRDFLB(NLP)=SRDFLB(NLP)+DKS0X
      SRUFLB(NLP)=SRUFLB(NLP)+DKS0X*RBNX
      SRDFLB(NL)=SRDFLB(NL)+DKS0X*(XANX+DX)
      SRUFLB(NL)=SRUFLB(NL)+DKS0X*UX
      SKDFLB(NLP,K)=DKS0X
      SKUFLB(NLP,K)=DKS0X*RBNX
      SKDFLB(NL,K)=DKS0X*(XANX+DX)
      SKUFLB(NL,K)=DKS0X*UX
      RBXTOA=RBNX
      SRKALB(K)=RBNX

C     Add successively layer N (at bottom) to form upper composite layer
C     ------------------------------------------------------------------

      DO 230 M=2,NL
      N=NLP-M
      XBNB=XNB(N)
      XBNX=XNX(N)
      RBNX=RNX(N)
      IF(RBNX.GT.1.D-05) GO TO 210
      RASB=RASB*XBNB*XBNB
      TANX=TANX*XBNB
      GO TO 220
  210 CONTINUE
      RBNB=RNB(N)
      TBNB=TNB(N)
      TBNX=TNX(N)
      RARB=RASB*RBNB
      XBTB=XBNB+TBNB
      DENOM=1.D0-RARB
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)/DENOM
      RASB=RBNB+XBTB*XBTB*RASB/DENOM
  220 CONTINUE
      XANX=XANX*XBNX
      RBNB=SRB(N)
      RBNX=SRX(N)
      DX=(TANX+XANX*RBNX*RASB)/(1.D0-RASB*RBNB)
      UX=RBNX*XANX+RBNB*DX
      SRUFLB(N)=SRUFLB(N)+DKS0X*UX
      SRDFLB(N)=SRDFLB(N)+DKS0X*(XANX+DX)
      SKUFLB(N,K)=DKS0X*UX
      SKDFLB(N,K)=DKS0X*(XANX+DX)
  230 CONTINUE

C     Record absorbed spectral flux at ground for surface type fractions
C     ------------------------------------------------------------------

      DO 231 I=1,4
      SRKGAX(K,I)=DKS0X*XANX*(1.D0-PRNX(KLAM,I))
      SRKGAD(K,I)=DKS0X*  DX*(1.D0-PRNB(KLAM,I))
  231 CONTINUE

      IF(K.EQ.NKSLAM) GO TO 301
      SRINIR=SRINIR+DKS0X
      SRONIR=SRONIR+DKS0X*RBXTOA
      SRDNIR=SRDNIR+SKDFLB(1,K)
      SRUNIR=SRUNIR+SKUFLB(1,K)
      SRRNIR=SRRNIR+DKS0X*RCNX
      SRTNIR=SRTNIR+DKS0X*(TANX+XANX)
      SRXNIR=SRXNIR+DKS0X*XANX
      GO TO 300

  301 CONTINUE
      SRIVIS=DKS0X
      SROVIS=DKS0X*RBXTOA
      SRDVIS=SKDFLB(1,K)
      SRUVIS=SKUFLB(1,K)
      SRRVIS=DKS0X*RCNX
      SRTVIS=DKS0X*(TANX+XANX)
      SRXVIS=SRXVIS+DKS0X*XANX

C     ------------------------------------------------------------------
C     UV absorption by O3 and O2 within solar spectral band DKS0(15)=.05
C     ------------------------------------------------------------------

      K=15
      DKS0X=DKS0(K)*S0COSZ
      SRKINC(K)=DKS0X

      N=NLP
      ATOPX=0.D0
      ATOPD=0.D0
      O3CMX=0.D0
      O3CMD=0.D0
  302 CONTINUE
      N=N-1
      O3CMX=O3CMX+COSMAG*ULGAS(N,3)
      O3CMD=O3CMD+1.90D0*ULGAS(N,3)
      CALL AO3ABS(O3CMX,ABOTX)
      CALL AO3ABS(O3CMD,ABOTD)
      AO3X(N)=(ABOTX-ATOPX)/DKS0(15)
      AO3D(N)=(ABOTD-ATOPD)/DKS0(15)
      ATOPX=ABOTX
      ATOPD=ABOTD
      IF(N.GT.1) GO TO 302
  303 CONTINUE
      O3CMX=O3CMX+1.90D0*ULGAS(N,3)
      O3CMD=O3CMD+1.90D0*ULGAS(N,3)
      CALL AO3ABS(O3CMX,ATOPX)
      CALL AO3ABS(O3CMD,ATOPD)
      AO3UXN=(ATOPX-ABOTX)/DKS0(15)
      AO3UDN=(ATOPD-ABOTD)/DKS0(15)
      AO3U(N)=XNX(N)*AO3UXN+(1.D0-XNX(N))*AO3UDN
      ABOTX=ATOPX
      ABOTD=ATOPD
      N=N+1
      IF(N.LT.NLP) GO TO 303
      RBNB=SRBALB(KLAM)
      RBNX=SRXALB(KLAM)
      RCNB=0.D0
      RCNX=0.D0
C                                  -------------------------------------
C                                  Get Oxygen UV absorption contribution
C                                  -------------------------------------
C----------------
      CALL GETO2A
C----------------
C                 ------------------------------------------------------
C                 Add Layers from Ground up. Retain Composite RBNB, RBNX
C             R,T,X of "A" (above) layer are corrected for O3 absorption
C             ----------------------------------------------------------

      DO 304 N=1,NL
      O2FHRL(N)=O2FHRL(N)/DKS0(15)*FULGAS(4)
      O2FHRB(N)=O2FHRB(N)/DKS0(15)*FULGAS(4)
      SRB(N)=RBNB
      SRX(N)=RBNX
      XANX=XNX(N)*(1.D0-AO3X(N)-O2FHRL(N))
      XANB=XNB(N)*(1.D0-AO3D(N)-O2FHRB(N))
      RASX=RNX(N)*(1.D0-AO3U(N))
      RASB=RNB(N)*(1.D0-AO3U(N))
      TANX=TNX(N)*(1.D0-AO3D(N))
      TANB=TNB(N)*(1.D0-AO3D(N))
!nu      ABSRTX=1.D0-XANX-TANX-RASX
!nu      ABSRTB=1.D0-XANB-TANB-RASB
      RARB=RASB*RBNB
      RARX=RASB*RBNX
      XATB=XANB+TANB
      DENOM=1.D0-RARB
      DB=(TANB+XANB*RARB)/DENOM
      DX=(TANX+XANX*RARX)/DENOM
      UB=RBNB*(XANB+DB)
      UX=RBNX*XANX+RBNB*DX
      RBNB=RASB+XATB*UB
      RBNX=RASX+XATB*UX
      XATC=XATB/(1.D0-RASB*RCNB)
      RCNX=RASX+(XANX*RCNX+TANX*RCNB)*XATC
      RCNB=RASB+RCNB*XATB*XATC
  304 CONTINUE
      VRD(NLP)=1.D0
      VRU(NLP)=RBNX
      SRKALB(15)=RBNX
      N=NL
      VRD(N)=XANX+DX
      VRU(N)=UX
  310 CONTINUE
      N=N-1
      XBNX=XNX(N)*(1.D0-AO3X(N)-O2FHRL(N))
      XBNB=XNB(N)*(1.D0-AO3D(N)-O2FHRB(N))
      RBNX=RNX(N)*(1.D0-AO3U(N))
      RBNB=RNB(N)*(1.D0-AO3U(N))
      TBNX=TNX(N)*(1.D0-AO3D(N))
      TBNB=TNB(N)*(1.D0-AO3D(N))

C     Add successively layer N (at bottom) to form upper composite layer
C     ------------------------------------------------------------------

      RARB=RASB*RBNB
      XBTB=XBNB+TBNB
      DENOM=1.D0/(1.D0-RARB)
      TANX=TBNX*XANX+XBTB*(TANX+XANX*RBNX*RASB)*DENOM
      RASB=RBNB+XBTB*XBTB*RASB*DENOM
      XANX=XANX*XBNX

C     Add upper & bottom composite layers to get flux at layer interface
C     ------------------------------------------------------------------

      RBNB=SRB(N)
      RBNX=SRX(N)
      DX=(TANX+XANX*RBNX*RASB)/(1.D0-RASB*RBNB)
      UX=RBNX*XANX+RBNB*DX
      VRD(N)=XANX+DX
      VRU(N)=UX
      IF(N.GT.1) GO TO 310
      DO 321 I=1,4
      SRKGAX(15,I)=DKS0X*XANX*(1.D0-PRNX(6,I))
      SRKGAD(15,I)=DKS0X*  DX*(1.D0-PRNB(6,I))
  321 CONTINUE

      DO 325 N=1,NLP
      VRD(N)=VRD(N)*DKS0X
      VRU(N)=VRU(N)*DKS0X
      SKDFLB(N,K)=VRD(N)
      SKUFLB(N,K)=VRU(N)
  325 CONTINUE
      SRIVIS=SRIVIS+VRD(NLP)
      SROVIS=SROVIS+VRU(NLP)
      PLAVIS=SROVIS/SRIVIS
      SRDVIS=SRDVIS+VRD(1)
      SRUVIS=SRUVIS+VRU(1)
      ALBVIS=SRUVIS/(SRDVIS+1.D-10)
      SRRVIS=SRRVIS+DKS0X*RCNX
      SRTVIS=SRTVIS+DKS0X*(TANX+XANX)
      SRXVIS=SRXVIS+DKS0X*XANX
      SRAVIS=1.D0-SRRVIS-SRTVIS

C     K16 strong absorbing contributions are computed without scattering
C     ------------------------------------------------------------------

      K=16
      DKS0X=DKS0(16)*S0COSZ
      SRKINC(16)=DKS0X
      SRKA16=0.D0
      DO 344 I=1,4
      SRKGAX(16,I)=0.D0
      SRKGAD(16,I)=0.D0
  344 CONTINUE
      DO 345 KK=1,3
      IF(KK.EQ.1) DKS0XX=DKS0X*0.002D0/0.011D0
      IF(KK.EQ.2) DKS0XX=DKS0X*0.008D0/0.011D0
      IF(KK.EQ.3) DKS0XX=DKS0X*0.001D0/0.011D0
      N=NLP
      TRNC=1.D0
      DO 330 M=1,NL
      N=N-1
      PLN=PL(N)
      CLX=DBLEXT(N,1)-DBLSCT(N,1)

C--------K=5-------CO2       DS0=.002
      IF(KK.EQ.1) THEN
      TRN1=0.D0
      ULN=ULGAS(N,2)*SECZ
      IF(ULN.GT.7.D0) ULN=7.D0
      TERMA=.003488*PLN*(1.+39.59*EXP(-8.769*ULN/(1.+4.419*ULN)))
     +     *(1.+ULN*(.001938*PLN-.00503*ULN))
      TERMB=(1.+.04712*PLN*(1.+.4877*ULN))
      TAUG=TERMA/TERMB*ULN
      TAU1=TAUG+CLX*SECZ
      IF(TAU1.LT.10.0) TRN1=EXP(-TAU1)
      FAC(N)=TRN1
      ENDIF

C--------K=7-------H2O       DS0=.008
      IF(KK.EQ.2) THEN
      TRN2=0.D0
      ULN=ULGAS(N,1)*SECZ
      TERMA=.001582*PLN*(1.+6.769*EXP(-9.59*ULN/(1.+5.026*ULN)))
     +     *(1.+ULN*(.2757E-03*PLN+.001429*ULN))
      TERMB=(1.+.003683*PLN*(1.+1.187*ULN))
      TAUG=TERMA/TERMB*ULN
      TAU2=TAUG+CLX*SECZ
      IF(TAU2.LT.10.0) TRN2=EXP(-TAU2)
      FAC(N)=TRN2
      ENDIF

C--------K=5-------O2        DS0=.001
      IF(KK.EQ.3) THEN
      TRN3=0.D0
      ULN=ULGAS(N,4)*SECZ
      TERMA=(.1366E-03-.2203E-07*TLN)*(1.+PLN*(.1497E-06*ULN+.001261))
      TERMB=(1.+.3867E-03*ULN)/(1.+.2075E-04*ULN)
      TAUG=TERMA/TERMB*ULN
      TAU3=TAUG+CLX*SECZ
      IF(TAU3.LT.10.0) TRN3=EXP(-TAU3)
      FAC(N)=TRN3
      ENDIF

      TRNC=TRNC*FAC(N)
      SRDFLB(N)=SRDFLB(N)+DKS0XX*TRNC
      SKDFLB(N,K)=SKDFLB(N,K)+DKS0XX*TRNC
  330 CONTINUE
      SRDFLB(NLP)=SRDFLB(NLP)+DKS0XX
      SRUFLB(1)=SRUFLB(1)+DKS0XX*TRNC*SRXALB(1)
      SKDFLB(NLP,K)=SKDFLB(NLP,K)+DKS0XX
      SKUFLB(1,K)=SKUFLB(1,K)+DKS0XX*TRNC*SRXALB(1)

C     For completeness, any incident flux at ground is relflected upward
C     ------------------------------------------------------------------

      TRNU=TRNC
      DO 340 N=2,NLP
      TRNU=TRNU*FAC(N-1)
      SRUFLB(N)=SRUFLB(N)+DKS0XX*TRNC*SRXALB(1)*TRNU
      SKUFLB(N,K)=SKUFLB(N,K)+DKS0XX*TRNC*SRXALB(1)*TRNU
  340 CONTINUE
      DO 341 I=1,4
      SRKGAX(16,I)=SRKGAX(16,I)+DKS0XX*TRNC*(1.D0-PRNX(1,I))
  341 CONTINUE
      SRKA16=SRKA16+TRNU*SRXALB(1)

      SRINIR=SRINIR+DKS0XX
      SRONIR=SRONIR+DKS0XX*TRNU*SRXALB(1)
      SRDNIR=SRDNIR+SKDFLB(1,K)
      SRUNIR=SRUNIR+SKUFLB(1,K)
  345 CONTINUE
      PLANIR=SRONIR/SRINIR
      ALBNIR=SRUNIR/(SRDNIR+1.D-10)
      SRKALB(16)=SRKA16/DKS0X

      DO 350 N=1,NLP
      SRDFLB(N)=SRDFLB(N)+VRD(N)
      SRUFLB(N)=SRUFLB(N)+VRU(N)
      SRNFLB(N)=SRDFLB(N)-SRUFLB(N)
  350 CONTINUE
      DO 360 N=1,NL
      SRFHRL(N)=SRNFLB(N+1)-SRNFLB(N)
  360 CONTINUE
      SRRNIR=SRRNIR+DKS0X*RCNX
      SRTNIR=SRTNIR+DKS0X*(TANX+XANX)
      SRXNIR=SRXNIR+DKS0X*XANX

      S0VIS=0.53D0*S0
      SRTVIS=SRTVIS/S0VIS
      SRRVIS=SRRVIS/S0VIS
      SRXVIS=SRXVIS/S0VIS
      SRAVIS=1.D0-SRTVIS-SRRVIS

      S0NIR=0.47D0*S0
      SRTNIR=SRTNIR/S0NIR
      SRRNIR=SRRNIR/S0NIR
      SRXNIR=SRXNIR/S0NIR
      SRANIR=1.D0-SRTNIR-SRRNIR


C     ------------------------------------------------------------------
C     FSRNFG defines the total solar flux absorbed at the ground surface
C              taking into account the albedo of different surface types
C     Thus:
C              SRNFLB(1)=POCEAN*FSRNFG(1)+PEARTH*FSRNFG(2)
C                       + POICE*FSRNFG(3)+ PLICE*FSRNFG(4)
C
C     NOTE:    If any surface type POCEAN, PEARTH, POICE, PLICE are Zero
C              the corresponding FSRNFG(I) absorbed solar flux at ground
C              is computed with that surface-type albedo set equal to 0.
C              ---------------------------------------------------------

      DO 420 I=1,4
      SUMX=0.D0
      SUMD=0.D0
      DO 410 K=1,16
      SUMX=SUMX+SRKGAX(K,I)
      SUMD=SUMD+SRKGAD(K,I)
  410 CONTINUE
      FSRNFG(I)=SUMX+SUMD
  420 CONTINUE


      DO 520 L=1,NLP
      DO 510 K=1,16
      SKNFLB(L,K)=SKDFLB(L,K)-SKUFLB(L,K)
  510 CONTINUE
  520 CONTINUE

      DO 540 L=1,NL
      DO 530 K=1,16
      SKFHRL(L,K)=SKNFLB(L+1,K)-SKNFLB(L,K)
  530 CONTINUE
  540 CONTINUE

      DO 560 L=1,NLP
      SUMD=0.D0
      SUMU=0.D0
      SUMN=0.D0
      SUMH=0.D0
      DO 550 K=1,16
      SUMD=SUMD+SKDFLB(L,K)
      SUMU=SUMU+SKUFLB(L,K)
      SUMN=SUMN+SKNFLB(L,K)
      SUMH=SUMH+SKFHRL(L,K)
  550 CONTINUE
      SKDFLB(L,17)=SUMD
      SKUFLB(L,17)=SUMU
      SKNFLB(L,17)=SUMN
      SKFHRL(L,17)=SUMH
  560 CONTINUE

      RETURN
      END SUBROUTINE SOLARM

      SUBROUTINE SETGTS
CCC   SUBROUTINE GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: GIN,TAUIN,RBBIN,EGIN
      INTEGER, INTENT(IN) :: KGTAUR
      REAL*8, INTENT(OUT) :: RBBOUT,TAUOUT
      REAL*8 FFKG(4,3),RBBK(3)
      REAL*8, PARAMETER, DIMENSION(14) :: GVALUE = (/.0,.25,.45,.50,.55,
     *     .60,.65,.70,.75,.80,.85,.90,.95,1./)
      REAL*8 CWM,CWE,TIJ,RBBI,RBB,BTAU,CUSPWM,CUSPWE,G,TAU,EG,DELTAU,TI
     *     ,WTJ,WTI,GI,WGI,WGJ,F1,F2,F3,F4,F21,F32,F43,F3221,F4332,A,B,C
     *     ,D,XF,FFCUSP,XE,XEXM,CUSPWT,FFLINR,RB2,RB3,TBB,TB2,TB3,XG,XM
     *     ,XP,RBBB,RI,WRJ,WRI,EI,WEI,WEJ,DELALB,X1,X2,X3,X4,XX,BB,DTAU
      INTEGER I,J,K,KTERPL,IT,JT,IG,JG,ITERPL,IGM,JGP,KG,IR,JR,IE,JE,IEM
     *     ,JEP

      DO 100 I=1,122
      TAUTGD(I)=(I-1)*0.1D0
      IF(I.GT.24) TAUTGD(I)=(I-24)*0.2D0+2.2D0
      IF(I.GT.48) TAUTGD(I)=(I-48)*0.5D0+7.0D0
      IF(I.GT.72) TAUTGD(I)=(I-72)+19.0D0
      IF(I.GT.96) TAUTGD(I)=(I-96)*5.0D0+40.0D0
      IF(I.GT.112) TAUTGD(I)=(I-112)*100.0D0+100.0D0
      IF(I.EQ.121) TAUTGD(I)=9999.99D0
      IF(I.EQ.122) TAUTGD(I)=12000.0D0
  100 CONTINUE

      DO 110 I=1,768
      IF(I.LT.602) TAUTGS(I)=(I-1)*0.05D0
      IF(I.GT.601) TAUTGS(I)=(I-601)*0.50D0+30.0D0
      IF(I.GT.741) TAUTGS(I)=(I-741)*50.0D0+100.D0
      IF(I.GT.758) TAUTGS(I)=(I-758)*1000.D0
  110 CONTINUE

      DO 130 J=1,13
      DO 120 I=1,768
      CWM=0.5
      CWE=0.5
      IF(I.GT.759) CWM=0.0
      IF(I.GT.759) CWE=0.0
      TIJ=TAUTGS(I)
      CALL SPLINE(TAUTGD,TGDATA(1,J),122,TIJ,RBBI,CWM,CWE,0)
      SALBTG(I,J)=RBBI
  120 CONTINUE
  130 CONTINUE
      DO 150 J=1,13
      DO 140 I=2,1000
      RBB=(I-1)*0.001D0
      CWM=0.5
      CWE=0.5
      CALL SPLINE(SALBTG(1,J),TAUTGS,768,RBB,BTAU,CWM,CWE,0)
      TAUGSA(I,J)=BTAU
  140 CONTINUE
  150 CONTINUE
      DO 160 J=1,14
      SALBTG(1,J)=0.D0
      TAUGSA(1,J)=0.D0
      TAUGSA(1001,J)=10000.D0
  160 CONTINUE

      DO 170 I=1,768
      SALBTG(I,14)=SALBTG(I,13)*2.D0-SALBTG(I,12)
  170 CONTINUE
      DO 180 I=1,1001
!nu   SALBI=(I-1)*0.001D0
      TAUGSA(I,14)=TAUGSA(I,13)*2.D0-TAUGSA(I,12)
  180 CONTINUE
      RETURN

      ENTRY GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)

      KTERPL=0
      CUSPWM=0.5
      CUSPWE=0.5

      G=GIN
      TAU=TAUIN
      RBB=RBBIN
      EG=EGIN

      RBBOUT=0.0
      TAUOUT=0.0
C                                           ---------------------------
C                                           OPTICAL DEPTH INTERPOLATION
C                                           0.05 ON (0.00 < TAU < 30.0)
C                                           0.50 ON (30.0 < TAU < 100.)
C                                           50.0 ON (100. < TAU < 1000)
C                                           ---------------------------

      IF(KGTAUR.EQ.2) GO TO 300

  200 CONTINUE

      DELTAU=0.05D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT.GT.599) GO TO 210
      WTJ=TI-IT
      WTI=1.D0-WTJ
      IT=IT+1
      GO TO 240
  210 CONTINUE
      DELTAU=0.50D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT.GT.199) GO TO 220
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+541
      GO TO 240
  220 CONTINUE
      DELTAU=50.0D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT.GT.19) GO TO 230
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+649
      GO TO 240
  230 CONTINUE
      DELTAU=1000.0D0
      TI=TAU/DELTAU
      IT=TI
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+758
  240 CONTINUE
      JT=IT+1

C                                    ---------------------------------
C                                    ASYMMETRY PARAMETER INTERPOLATION
C                                    0.05 CUBIC SPLINE (0.5 < G < 0.9)
C                                    0.25 QUADRATIC ON (0.0 < G < 0.5)
C                                    LINEAR EXTRAP FOR (.95 < G < 1.0)
C                                    ---------------------------------

      GI=G*20.D0
      IF(GI.GT.10.0) GO TO 250
      IG=2
      JG=3
      ITERPL=1
      GO TO 260

  250 CONTINUE
      ITERPL=4
      IG=GI
      WGJ=GI-IG
      WGI=1.D0-WGJ
      IG=IG-6
      IF(IG.GT.12) THEN
      ITERPL=2
      IG=12
      ENDIF
      JG=IG+1

  260 CONTINUE

      IGM=IG-1
      JGP=JG+1

      K=0
      DO 270 KG=IGM,JGP
      K=K+1
      F1=SALBTG(IT-1,KG)
      F2=SALBTG(IT  ,KG)
      F3=SALBTG(JT  ,KG)
      F4=SALBTG(JT+1,KG)
      IF(IT.EQ.1) F1=-F3
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WTJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      FFKG(K,1)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      FFKG(K,2)=F2
      FFKG(K,3)=F3
  270 CONTINUE

      IF(ITERPL.LT.4) GO TO 290

      DO 280 K=1,3
      F1=FFKG(1,K)
      F2=FFKG(2,K)
      F3=FFKG(3,K)
      F4=FFKG(4,K)
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WGJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      RBBK(K)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
  280 CONTINUE
      RBB=RBBK(1)
      RB2=RBBK(2)
      RB3=RBBK(3)
      TBB=TAU
      TB2=TAUTGS(IT)
      TB3=TAUTGS(JT)
      IF(KGTAUR.EQ.1) RETURN
      IF(KTERPL.EQ.1) GO TO 400
      GO TO 300

  290 CONTINUE
      XG=G*2.D0-0.5D0
      IF(ITERPL.EQ.2) XG=G*10.D0-9.D0
      XM=1.D0-XG-XG
      XP=1.D0+XG+XG
      RBB=XM*XP*FFKG(ITERPL+1,1)-XG*XM*FFKG(ITERPL,1)+XG*XP*FFKG(4,1)
      RB2=XM*XP*FFKG(ITERPL+1,2)-XG*XM*FFKG(ITERPL,2)+XG*XP*FFKG(4,2)
      RB3=XM*XP*FFKG(ITERPL+1,3)-XG*XM*FFKG(ITERPL,3)+XG*XP*FFKG(4,3)

      IF(KGTAUR.EQ.1) RETURN
      IF(KTERPL.EQ.1) GO TO 400

  300 CONTINUE
      RBBB=RBB

      RI=RBB*1000.D0
      IR=RI
      WRJ=RI-IR
      WRI=1.D0-WRJ
      IR=IR+1
      JR=IR+1

      EI=EG*20.D0
      IF(EI.GT.10.0) GO TO 310
      IE=2
      JE=3
      ITERPL=1
      GO TO 320
  310 CONTINUE

      ITERPL=4
      IE=EI
      WEJ=EI-IE
      WEI=1.D0-WEJ
      IE=IE-6
      IF(IE.GT.12) THEN
      ITERPL=2
      IE=12
      ENDIF
      JE=IE+1
  320 CONTINUE

      DELALB=0.001D0
      IEM=IE-1
      JEP=JE+1
      K=0
      DO 330 KG=IEM,JEP
      K=K+1
      F1=TAUGSA(IR-1,KG)
      F2=TAUGSA(IR  ,KG)
      F3=TAUGSA(JR  ,KG)
      F4=TAUGSA(JR+1,KG)
      IF(IR.EQ.1) F1=-F3
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WRJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      FFKG(K,1)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      FFKG(K,2)=F2
      FFKG(K,3)=F3
  330 CONTINUE
      X1=GVALUE(IE-1)
      X2=GVALUE(IE  )
      X3=GVALUE(JE  )
      X4=GVALUE(JE+1)
      XX=WEJ
      IF(ITERPL.LT.4) GO TO 350

      DO 340 K=1,3
      F1=FFKG(1,K)
      F2=FFKG(2,K)
      F3=FFKG(3,K)
      F4=FFKG(4,K)
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WEJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      RBBK(K)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
  340 CONTINUE
      TBB=RBBK(1)
      TB2=RBBK(2)
      TB3=RBBK(3)

      IF(KTERPL.EQ.1) GO TO 400
      GO TO 390

  350 CONTINUE
      XG=EG*2.D0-0.5D0
      IF(ITERPL.EQ.2) XG=G*10.D0-9.D0
      XM=1.D0-XG-XG
      XP=1.D0+XG+XG
      TBB=XM*XP*FFKG(ITERPL+1,1)-XG*XM*FFKG(ITERPL,1)+XG*XP*FFKG(4,1)
      TB2=XM*XP*FFKG(ITERPL+1,2)-XG*XM*FFKG(ITERPL,2)+XG*XP*FFKG(4,2)
      TB3=XM*XP*FFKG(ITERPL+1,3)-XG*XM*FFKG(ITERPL,3)+XG*XP*FFKG(4,3)
      IF(KTERPL.EQ.1) GO TO 400
  390 CONTINUE
      KTERPL=1
      TAU=TBB
      G=EGIN
      GO TO 200
  400 CONTINUE
      IF(ABS(WTI*WTJ).LT.0.1D0) DTAU=(RBBB-RB2)/(RB3-RB2)
      IF(ABS(WTI*WTJ).GE.0.1D0) THEN
      C=(RB3-RBB)/WTI-(RBB-RB2)/WTJ
      B=(RBB-RB2)/WTJ-WTJ*C
      A=RB2
      BB=B*B+4.D0*C*(RBBB-A)
      IF(BB.GT.0.D0) DTAU=(SQRT(BB)-B)/(C+C)
      ENDIF
      TAUOUT=(IT-1+DTAU)*DELTAU
      RBBOUT=RBBB

      RETURN
CCC   END SUBROUTINE GTSALB
      END SUBROUTINE SETGTS

      SUBROUTINE SGPGXG(XMU,TAU,G,GG)
      IMPLICIT NONE
C     ----------------------------------------------------------------
C     COSBAR ADJUSTMENT TO REPRODUCE THE SOLAR ZENITH ANGLE DEPENDENCE
C     FOR AEROSOL ALBEDO FOR OPTICAL THICKNESSES [0.0 < TAU < 10000.0]
C     ----------------------------------------------------------------
      REAL*8, INTENT(IN) :: XMU,TAU,G
      REAL*8, INTENT(OUT) :: GG
      REAL*8 XI,WXI,WXJ,GI,WGI,WGJ,TI,WTJ,WTI
      INTEGER IX,JX,IG,JG,IT,IT0,JT
C                          -------------------------------------------
C                          XMU (COSZ) SOLAR ZENITH ANGLE INTERPOLATION
C                          DATA INTERVAL:  0.02  ON  [0.0 < XMU < 1.0]
C                          -------------------------------------------

      XI=XMU*50.D0+0.999999D0  ! >1 since XMU=COSZ>.001
      IX=XI
      JX=IX+1
      WXJ=XI-IX
      WXI=1.D0-WXJ

C                                      -------------------------------
C                                      COSBAR DEPENDENCE INTERPOLATION
C                                         0.10 ON [0.0 < COSBAR < 1.0]
C                                      -------------------------------

      GI=G*10.D0
      IG=GI
      WGJ=GI-IG
      WGI=1.D0-WGJ
      IG=IG+1
      JG=IG+1

C                            -----------------------------------------
C                               AEROSOL TAU INTERPOLATION INTERVALS
C                            -----------------------------------------
C                            dTau      1      1  (Lin Int)  61     62
C                            0.10 ON [0.00 , 0.00 < TAU < 6.00 , 6.10]
C                                     63     64            92     93
C                            0.50 ON [5.50 , 6.00 < TAU < 20.0 , 20.5]
C                                     94     95           111    112
C                            5.00 ON [15.0 , 20.0 < TAU < 100. , 105.]
C                                     113    114          132    133
C                            50.0 ON [50.0 , 100. < TAU < 1000 , 1050]
C                                            134          143
C                            1000 ON [     , 1000 < TAU < 10000,     ]
C                            -----------------------------------------

      IF(TAU.LT.6.D0) THEN
      TI=TAU*10.D0+1.
      IT=TI
      WTJ=TI-IT
      IT0=0

      ELSEIF(TAU.LT.20.D0) THEN
      TI=(TAU-6.D0)*2.00D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=62

      ELSEIF(TAU.LT.100.D0) THEN
      TI=(TAU-20.D0)*0.20D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=93

      ELSEIF(TAU.LT.1000.D0) THEN
      TI=(TAU-100.D0)*0.02D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=112

      ELSE
      TI=TAU*0.001D0+1.D-6
      IT=TI
      WTJ=TI-IT
      IF(IT.GT.9) IT=9
      IT0=133
      ENDIF

      WTI=1.D0-WTJ
      IT=IT+IT0
      JT=IT+1
      GG=WGI*(WTI*(WXI*GTAU(IX,IG,IT)+WXJ*GTAU(JX,IG,IT))
     +      + WTJ*(WXI*GTAU(IX,IG,JT)+WXJ*GTAU(JX,IG,JT)))
     +  +WGJ*(WTI*(WXI*GTAU(IX,JG,IT)+WXJ*GTAU(JX,JG,IT))
     +      + WTJ*(WXI*GTAU(IX,JG,JT)+WXJ*GTAU(JX,JG,JT)))

      RETURN
      END SUBROUTINE SGPGXG

      SUBROUTINE BCTAUW(XJYEAR,IDEC,JDEC,BCWTID,BCWTJD)
      IMPLICIT NONE
C-------------------------------------------------------------------
C     Black Carbon interdecadal TAU interpolation is based on linear
C     TAU trend (between decadal global TAUmaps) with a superimposed
C     intra-decadal time dependence scaled to the Black Carbon Total
C     emission rate.
C
C        INPUT:  XJYEAR  (Fractional Julian year)
C                                                 by 25     by 10
C       OUTPUT:  IDEC    (Map Index: I= -3,4-7  (0) -1925,1950-1980)
C                JDEC    (Map Index: J=1-4,5-8  1875-1950,1960-1990)
C
C              BCWTID    (Multiplicative Weight for BC TAU-Map IDEC)
C              BCWTJD    (Multiplicative Weight for BC TAU-Map JDEC)
C
C-------------------------------------------------------------------

      REAL*8, DIMENSION(50) :: EYEAR,BCEHC,BCEBC,BCEDI,BCTOT

      REAL*8 BCE(5,45)

C    Global Annual Emissions of BC U   Emission (Mt/yr)

      DATA BCE/
C      Year    Hard_Coal    Brown_Coal      Diesel        Total
     A 50.0, 2.280581713, 0.4449132979, 0.1599090248, 2.885536671,
     B 51.0, 2.443193913, 0.4855868816, 0.1884280443, 3.117194653,
     C 52.0, 2.473641872, 0.5115299225, 0.2027695477, 3.187930107,
     D 53.0, 2.481340885, 0.5448409319, 0.2149295360, 3.241089582,
     E 54.0, 2.505670071, 0.5780177116, 0.2343477309, 3.317960978,
     F 55.0, 2.698692560, 0.6238067150, 0.2733324766, 3.595800638,
     G 56.0, 2.855226278, 0.6531309485, 0.3043369055, 3.812692404,
     H 57.0, 2.975781679, 0.6821750998, 0.3207367063, 3.978575468,
     I 58.0, 3.341105223, 0.7035279870, 0.3370627165, 4.381746292,
     J 59.0, 3.638528824, 0.7075053453, 0.3695519567, 4.715488434,
     K 60.0, 3.770926714, 0.7416650057, 0.3832504749, 4.896034241,
     L 61.0, 3.392980337, 0.7805693150, 0.4217525721, 4.595387459,
     M 62.0, 3.288835049, 0.8179932237, 0.4603823125, 4.567360401,
     N 63.0, 3.359177589, 0.8604368567, 0.5090782642, 4.728550911,
     O 64.0, 3.432664871, 0.8952696323, 0.5388473868, 4.866865158,
     P 65.0, 3.529418945, 0.8819132447, 0.5785927773, 4.989773750,
     Q 66.0, 3.577459812, 0.8817394972, 0.6323299408, 5.091631413,
     R 67.0, 3.418204546, 0.8635972142, 0.6592246890, 4.941041946,
     S 68.0, 3.452457905, 0.8943673372, 0.7338049412, 5.080585003,

ceq   DATA BCE2/
     A 69.0, 3.626069546, 0.9298774004, 0.7889106274, 5.344810009,
     B 70.0, 3.264039755, 0.9229136109, 0.8880128860, 5.074741840,
     C 71.0, 3.437611580, 0.9374827743, 0.9531223178, 5.328329086,
     D 72.0, 3.473345757, 0.7836616039, 1.0180075170, 5.274850368,
     E 54.0, 2.505670071, 0.5780177116, 0.2343477309, 3.317960978,
     F 74.0, 3.506143808, 0.8251076341, 1.0828053950, 5.413989067,
     G 75.0, 3.906814098, 0.8527192473, 1.0454736950, 5.804963112,
     H 76.0, 4.005736828, 0.8900613785, 1.1400985720, 6.035901546,
     I 77.0, 4.236912251, 0.9103702307, 1.2190728190, 6.366260529,
     J 78.0, 4.459666252, 0.9303293228, 1.2408012150, 6.630728722,
     K 79.0, 4.697422504, 0.9856286645, 1.3019220830, 6.984815121,
     L 80.0, 4.796229839, 0.9959300756, 1.2336660620, 7.026207924,
     M 81.0, 4.789204121, 1.0459070210, 1.1664049630, 7.001126766,
     N 82.0, 4.872739315, 1.0975246430, 1.1601715090, 7.130136490,
     O 83.0, 4.983223438, 1.1424025300, 1.1732926370, 7.298912525,
     P 84.0, 5.265352249, 1.2178678510, 1.2251536850, 7.708741188,
     Q 85.0, 5.763637543, 1.2965050940, 1.2428865430, 8.303324699,
     R 86.0, 5.924767494, 1.3386499880, 1.2930148840, 8.556744576,
     S 87.0, 6.155550480, 1.3738890890, 1.3162037130, 8.845513344,

ceq   DATA BCE3/
     A 88.0, 6.379704475, 1.3670797350, 1.3813229800, 9.127896309,
     B 89.0, 6.594299316, 1.4169263840, 1.4029121400, 9.414231300,
     C 90.0, 6.566919804, 1.4685817960, 1.4224120380, 9.458042145,
     D 91.0, 6.661097050, 1.2067918780, 1.4163945910, 9.284657478,
     E 92.0, 7.737902641, 1.3509917260, 1.4471185210, 10.53625107,
     F 93.0, 7.393332005, 1.2448183300, 1.4543261530, 10.09271908,
     G 94.0, 7.515841007, 1.2333894970, 1.4780857560, 10.22745800/

C     REAL*8, PARAMETER :: POST90=50.D0
C     REAL*8, PARAMETER :: PRE50 =50.D0
      REAL*8, PARAMETER :: POST90=-250.D0
      INTEGER, SAVE :: IFIRST=1
      INTEGER, INTENT(OUT) :: IDEC,JDEC
      REAL*8, INTENT(OUT) :: BCWTID,BCWTJD
      REAL*8, INTENT(IN) :: XJYEAR
      INTEGER I,IYDI,IYDJ,IYEAR,IYYI,IYYJ
      REAL*8 XDEC,DYEAR,DDEC,DELTAY,DELDEC,BCED,BCEY,RATYD
      SAVE EYEAR,BCEHC,BCEBC,BCEDI,BCTOT,BCE

      IF(IFIRST.EQ.1) THEN
      DO 110 I=1,45
      EYEAR(I)=BCE(1,I)+1900.D0
      BCEHC(I)=BCE(2,I)
      BCEBC(I)=BCE(3,I)
      BCEDI(I)=BCE(4,I)
      BCTOT(I)=BCE(5,I)
  110 CONTINUE
      IFIRST=0
      ENDIF

      IF(XJYEAR.LT.1950.D0) THEN
C**** 0 until 1850, then lin. interpolate obs. data (every 25 years)
        XDEC=(XJYEAR-1850.d0)/25.d0
        IF (XDEC.lt.0.) Xdec=0.
        IDEC=XDEC
        DDEC=XDEC-IDEC
        BCWTID=1.-DDEC
        JDEC=IDEC+1
        BCWTJD=DDEC
        IF(IDEC.LT.1) THEN
          IDEC=1
          JDEC=1
          BCWTID=XDEC
          BCWTJD=0.
        END IF
      ELSE IF(XJYEAR.GE.1990.D0) THEN
C**** Slow reduction after 1990 (POST90=-250 years e-folding time)
C**** Actually we will use no reduction after 1990
        DYEAR=XJYEAR-1990.D0
        BCWTID=0.D0
        BCWTJD=1.D0 !  BCWTJD=1.D0*EXP(DYEAR/POST90)
        IDEC=8
        JDEC=8
      ELSE
C**** lin. interpolate obs. data (every 10 years) 1950-1990
        DYEAR=XJYEAR-1950.D0
        IYEAR=DYEAR
        DELTAY=DYEAR-IYEAR
        DELDEC=DYEAR/10.D0-IYEAR/10
        IDEC=(IYEAR+10)/10
        JDEC=IDEC+1
        IYDI=1+(IDEC-1)*10
        IYDJ=IYDI+10
        IYYI=IYEAR+1
        IYYJ=IYYI+1
        BCED=BCTOT(IYDI)+DELDEC*(BCTOT(IYDJ)-BCTOT(IYDI))
        BCEY=BCTOT(IYYI)+DELTAY*(BCTOT(IYYJ)-BCTOT(IYYI))
        RATYD=BCEY/BCED
        BCWTID=RATYD*(1.D0-DELDEC)
        BCWTJD=RATYD*DELDEC
        IDEC=IDEC+3
        JDEC=JDEC+3
      END IF

      RETURN
      END SUBROUTINE BCTAUW

      SUBROUTINE SUTAUW(XJYEAR,IDEC,JDEC,SUWTID,SUWTJD)
      IMPLICIT NONE
C-------------------------------------------------------------------
C     Anthropogenic Sulfate inter-decadal TAU interpolation is based
C     on a linear TAU trend (between decadal global TAU-maps) with a
C     superimposed intradecadal time dependence scaled in proportion
C     to the Anthropogenic Sulfate global emission rate.
C
C        INPUT:  XJYEAR  (Fractional Julian year)
C                                                 by 25     by 10
C       OUTPUT:  IDEC    (Map Index: I= -3,4-7  (0) -1925,1950-1980)
C                JDEC    (Map Index: J=1-4,5-8  1875-1950,1960-1990)
C
C              SUWTID    (Multiplicative Weight for SU TAU-Map IDEC)
C              SUWTJD    (Multiplicative Weight for SU TAU-Map JDEC)
C
C-------------------------------------------------------------------

      REAL*8, DIMENSION(50) :: EYEAR,SUANT,SUNAT

      REAL*8 SUE(3,41)

C     Global Emission of Sulfate

C     Emission (Mt/yr)
C               year      Anthropogenic_Sulfate Natural_Sulfate
      DATA SUE/
     A          1950.0,     30.46669769,           14.4,
     B          1951.0,     32.38347244,           14.4,
     C          1952.0,     32.18632889,           14.4,
     D          1953.0,     32.83379745,           14.4,
     E          1954.0,     32.79270935,           14.4,
     F          1955.0,     35.79611969,           14.4,
     G          1956.0,     39.93603897,           14.4,
     H          1957.0,     38.68806839,           14.4,
     I          1958.0,     39.35904312,           14.4,
     J          1959.0,     41.06065369,           14.4,
     K          1960.0,     42.67050934,           14.4,
     L          1961.0,     41.32410431,           14.4,
     M          1962.0,     41.80470276,           14.4,
     N          1963.0,     43.26312637,           14.4,
     O          1964.0,     44.68368530,           14.4,
     P          1965.0,     45.81701660,           14.4,
     Q          1966.0,     46.61584091,           14.4,
     R          1967.0,     46.42276001,           14.4,
     S          1968.0,     47.77438354,           14.4,
ceq   DATA SUE2/
     A          1969.0,     49.30817032,           14.4,
     B          1970.0,     52.81050873,           14.4,
     C          1971.0,     52.95043945,           14.4,
     D          1972.0,     54.10167694,           14.4,
     E          1973.0,     55.93037415,           14.4,
     F          1974.0,     57.31056213,           14.4,
     G          1975.0,     58.52788162,           14.4,
     H          1976.0,     59.71361542,           14.4,
     I          1977.0,     62.59599304,           14.4,
     J          1978.0,     61.98198318,           14.4,
     K          1979.0,     64.71042633,           14.4,
     L          1980.0,     65.28986359,           14.4,
     M          1981.0,     63.23768234,           14.4,
     N          1982.0,     62.88000488,           14.4,
     O          1983.0,     61.45023346,           14.4,
     P          1984.0,     63.85008621,           14.4,
     Q          1985.0,     66.47412872,           14.4,
     R          1986.0,     68.00902557,           14.4,
     S          1987.0,     69.87956238,           14.4,
ceq   DATA SUE3/
     A          1988.0,     70.52937317,           14.4,
     B          1989.0,     72.06355286,           14.4,
     C          1990.0,     71.29174805,           14.4/

C     REAL*8, PARAMETER :: POST90=100.D0
C     REAL*8, PARAMETER :: PRE50=50.D0
      REAL*8, PARAMETER :: POST90=-250.D0
      INTEGER, SAVE :: IFIRST=1
      SAVE  EYEAR,SUANT,SUNAT,SUE
      REAL*8, INTENT(IN) :: XJYEAR
      REAL*8, INTENT(OUT) :: SUWTID,SUWTJD
      INTEGER, INTENT(OUT) :: IDEC,JDEC
      REAL*8 XDEC,DDEC,DYEAR,SUED,SUEY,RATYD,DELDEC,DELTAY
      INTEGER I,IYEAR,IYDI,IYDJ,IYYI,IYYJ

      IF(IFIRST.EQ.1) THEN
      DO 110 I=1,41
      EYEAR(I)=SUE(1,I)
      SUANT(I)=SUE(2,I)
      SUNAT(I)=SUE(3,I)
  110 CONTINUE
      IFIRST=0
      ENDIF

      IF(XJYEAR.LT.1950.D0) THEN
C**** 0 until 1850, then lin. interpolate obs. data (every 25 years)
        XDEC=(XJYEAR-1850.d0)/25.d0
        IF (XDEC.lt.0.) Xdec=0.
        IDEC=XDEC
        DDEC=XDEC-IDEC
        SUWTID=1.-DDEC
        JDEC=IDEC+1
        SUWTJD=DDEC
        IF(IDEC.LT.1) THEN
          IDEC=1
          JDEC=1
          SUWTID=XDEC
          SUWTJD=0.
        END IF
      ELSE IF(XJYEAR.GE.1990.D0) THEN
C**** Slow reduction after 1990 (POST90=-250 years e-folding time)
C**** Actually we will use no reduction after 1990
        DYEAR=XJYEAR-1990.D0
        SUWTID=0.D0
        SUWTJD=1.D0  !  SUWTJD=1.D0*EXP(DYEAR/POST90)
        IDEC=8
        JDEC=8
      ELSE
C**** lin. interpolate obs. data (every 10 years) 1950-1990
        DYEAR=XJYEAR-1950.D0
        IYEAR=DYEAR
        DELTAY=DYEAR-IYEAR
        DELDEC=DYEAR/10.D0-IYEAR/10
        IDEC=(IYEAR+10)/10
        JDEC=IDEC+1
        IYDI=1+(IDEC-1)*10
        IYDJ=IYDI+10
        IYYI=IYEAR+1
        IYYJ=IYYI+1
        SUED=SUANT(IYDI)+DELDEC*(SUANT(IYDJ)-SUANT(IYDI))
        SUEY=SUANT(IYYI)+DELTAY*(SUANT(IYYJ)-SUANT(IYYI))
        RATYD=SUEY/SUED
        SUWTID=RATYD*(1.D0-DELDEC)
        SUWTJD=RATYD*DELDEC
        IDEC=IDEC+3
        JDEC=JDEC+3
      END IF

      RETURN
      END SUBROUTINE SUTAUW

      SUBROUTINE GETMIE(NA,AREFF,SQEX,SQSC,SQCB,TQAB,Q55)

c     INCLUDE 'rad00def.radCOMMON.f'

      INTEGER, INTENT(IN) :: NA
      real*8,  intent(in) :: areff
      real*8   SQEX(6),SQSC(6),SQCB(6),TQEX(33),TQSC(33),TQAB(33),Q55
      real*8   QXAERN(25),QSAERN(25),QGAERN(25),Q55AER(25)

      real*8 wts,wta,QGAERX,pi,vreff
      integer n0,k,n,nn
                          !                               1   2   3   4
      IF(NA.LT.5) THEN    !    NA : Aerosol compositions SO4,SEA,ANT,OCX
      N0=0
      IF(NA.EQ.2) N0=22
      IF(NA.EQ.3) N0=44
      IF(NA.EQ.4) N0=88
      DO 112 K=1,6
      DO 111 N=1,22
      NN=N0+N
      WTS=FRSULF(NA)
      WTA=1.D0-WTS
      QXAERN(N)=SRUQEX(K,NN)*WTA+SRUQEX(K,N)*WTS
      QSAERN(N)=SRUQSC(K,NN)*WTA+SRUQSC(K,N)*WTS
      QGAERX=SRUQCB(K,NN)*SRUQSC(K,NN)*WTA+SRUQCB(K,N)*SRUQSC(K,N)*WTS
      QGAERN(N)=QGAERX/QSAERN(N)
  111 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,AREFF,SQEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,AREFF,SQSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,AREFF,SQCB(K),1.D0,1.D0,1)

      PI=SQSC(K)/SQEX(K)
      IF(PI.GT.PI0MAX(NA)) SQSC(K)=SQSC(K)*PI0MAX(NA)/PI
  112 CONTINUE
      DO 114 K=1,33
      DO 113 N=1,22
      NN=N0+N
      WTS=FRSULF(NA)
      WTA=1.D0-WTS
      QXAERN(N)=TRUQEX(K,NN)*WTA+TRUQEX(K,N)*WTS
      QSAERN(N)=TRUQSC(K,NN)*WTA+TRUQSC(K,N)*WTS
      QGAERX=TRUQCB(K,NN)*TRUQSC(K,NN)*WTA+TRUQCB(K,N)*TRUQSC(K,N)*WTS
      QGAERN(N)=QGAERX/(QSAERN(N)+1.d-20)
  113 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC(K),1.D0,1.D0,1)
      TQAB(K)=TQEX(K)-TQSC(K)
  114 CONTINUE
      DO 115 N=1,22
      NN=N0+N
      WTS=FRSULF(NA)
      WTA=1.D0-WTS
      Q55AER(N)=Q55U22(NN)*WTA+Q55U22(N)*WTS
  115 CONTINUE
      CALL SPLINE(REFU22,Q55U22,22,AREFF,Q55,1.D0,1.D0,1)
      ENDIF

                                   !                              5   6
      IF(NA.EQ.5.OR.NA.EQ.6) THEN  !   NA : Aerosol compositions BIC,BCB
cc    AREFF=REFDRY(NA)
      DO 122 K=1,6
      DO 121 N=1,25
      QXAERN(N)=SRSQEX(K,N)
      QSAERN(N)=SRSQSC(K,N)
      QGAERN(N)=SRSQCB(K,N)
  121 CONTINUE
      CALL SPLINE(REFS25,QXAERN,25,AREFF,SQEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QSAERN,25,AREFF,SQSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QGAERN,25,AREFF,SQCB(K),1.D0,1.D0,1)
  122 CONTINUE
      DO 124 K=1,33
      DO 123 N=1,25
      QXAERN(N)=TRSQEX(K,N)
      QSAERN(N)=TRSQSC(K,N)
      QGAERN(N)=TRSQCB(K,N)
  123 CONTINUE
      CALL SPLINE(REFS25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
      TQAB(K)=TQEX(K)-TQSC(K)
  124 CONTINUE
      CALL SPLINE(REFS25,Q55S25,25,AREFF,Q55,1.D0,1.D0,1)
      ENDIF

                                        !                             7
      IF(NA.EQ.7) THEN                  !   NA : Aerosol composition DST
cc    AREFF=REFDRY(NA)
      DO 132 K=1,6
      DO 131 N=1,25
      QXAERN(N)=SRDQEX(K,N)
      QSAERN(N)=SRDQSC(K,N)
      QGAERN(N)=SRDQCB(K,N)
  131 CONTINUE
      CALL SPLINE(REFS25,QXAERN,25,AREFF,SQEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QSAERN,25,AREFF,SQSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QGAERN,25,AREFF,SQCB(K),1.D0,1.D0,1)
  132 CONTINUE
      DO 134 K=1,33
      DO 133 N=1,25
      QXAERN(N)=TRDQEX(K,N)
      QSAERN(N)=TRDQSC(K,N)
      QGAERN(N)=TRDQCB(K,N)
  133 CONTINUE
      CALL SPLINE(REFS25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFS25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
      TQAB(K)=TQEX(K)-TQSC(K)
  134 CONTINUE
      CALL SPLINE(REFD25,Q55D25,25,AREFF,Q55,1.D0,1.D0,1)
      ENDIF

                               !                                      8
      IF(NA.EQ.8) THEN         !     NA : Aerosol composition(H2SO4) VOL
      VREFF=AREFF
      IF(VREFF.LT.0.1D0) VREFF=0.1D0
      IF(VREFF.GT.2.0D0) VREFF=2.0D0
      CALL GETQVA(VREFF)
      DO 141 K=1,6
      SQEX(K)=QVH2S(K)
      SQSC(K)=SVH2S(K)
      SQCB(K)=GVH2S(K)
  141 CONTINUE
      DO 142 K=1,33
      TQAB(K)=AVH2S(K)
  142 CONTINUE
      Q55=Q55H2S
      ENDIF
      RETURN
      END SUBROUTINE GETMIE

      SUBROUTINE AO3ABS(OCM,O3ABS)
      IMPLICIT NONE
C              ---------------------------------------------------------
C              UV absorption by Ozone  is expressed as a fraction of the
C              total solar flux S0. Hence O3ABS (fraction of total solar
C              flux absored by OCM cm ofozone) must be normalized within
C              SOLARM by dividing O3ABS by the corresponding fraction of
C              the solar flux within the spectral interval DKS0(15)=0.05
C              ---------------------------------------------------------
      REAL*8, INTENT(IN) :: OCM
      REAL*8, INTENT(OUT) :: O3ABS
      REAL*8 XX,DX
      INTEGER IP,IX

      O3ABS=AO3(460)
      IP=0
      XX=OCM*1.D+04
      IX=XX
      IF(IX.GT.99) GO TO 110
      IF(IX.LT.1 ) GO TO 130
      GO TO 120
  110 CONTINUE
      IP=IP+90
      XX=XX*0.1D0
      IX=XX
      IF(IX.GT.99) GO TO 110
  120 CONTINUE
      DX=XX-IX
      IX=IX+IP
      IF(IX.GT.459) GO TO 140
      O3ABS=AO3(IX)+DX*(AO3(IX+1)-AO3(IX))
      GO TO 140
  130 CONTINUE
      O3ABS=XX*AO3(1)
  140 CONTINUE

      RETURN
      END SUBROUTINE AO3ABS

      SUBROUTINE WRITER(KWRU,INDEX)
C
      USE SURF_ALBEDO, only : AVSCAT, ANSCAT, AVFOAM, ANFOAM,
     *     WETTRA, WETSRA, ZOCSRA, ZSNSRA, ZICSRA, ZDSSRA, ZVGSRA,
     *     EOCTRA, ESNTRA, EICTRA, EDSTRA, EVGTRA, AGEXPF, ALBDIF
      IMPLICIT NONE
C
C     ------------------------------------------------------------------
C     WRITER Radiative Input/Output Cloud/Aerosol Data/Conrol Parameters
C
C         INDEX
C           0     control parameter defaults in RADPAR
C           1     RADPAR Radiative control/scaling param's; GHG defaults
C           2     RADPAR Atmospheric composition P,H,T,Cld,Aer  profiles
C           3     RADPAR Computed LW SW fluxes cooling and heating rates
C           4     Aerosol and Cloud: Mie scattering radiative parameters
C                 A  SW aerosol Mie scattering Qx,Qs,g in use parameters
C                 B  SW  cloud  Mie scattering Qx,Qs,g in use parameters
C                 C  SW cld+aer Mie scattering Qx,Qs,g in use parameters
C                 D  SW LW aerosol 11-compositon  Mie Qx,Qs,g parameters
C                 E  SW LW aerosol  6-compositon  Mie Qx,Qs,g parameters
C                 F  SW LW aerosol 8-size D dust  Mie Qx,Qs,g parameters
C                 G  SW LW  cloud  15-size/phase  Mie Qx,Qs,g parameters
C           5     LW cld,aer,gas total optical k-distribution extinction
C           6     LW gas absorb: total optical k-distribution extinction
C           7     A  LW cloud   TRCALK optical k-distribution extinction
C                 B  LW aerosol TRAALK optical k-distribution extinction
C           8     SW Spectral/k-dist flux, albedo, absorption components
C                 A  Spectral components of downward & upward solar flux
C                 B  Spectral components of net solar flux, heating rate
C           9     LW flux contribution from each k-distribution interval
C                 1  Downward LW flux  from each k-distribution interval
C                 2  Upward   LW flux  from each k-distribution interval
C                 3  Net (Up) LW flux  from each k-distribution interval
C                 4  Flux cooling rate from each k-distribution interval
C                 5  Fraction coolrate from each k-distribution interval
C        NOTE:
C                 KWTRAB sets LW Mie parameters in 4-D,E,F,G
C                        KWTRAB=0 (default) sets LW output to be Mie Qab
C                        KWTRAB=1           sets LW output to be Mie Qex
C                        KWTRAB=2           sets LW output to be Mie Qsc
C                        KWTRAB=3           sets LW output to be Mie Qcb
C                        KWTRAB=4           sets LW output to be Mie Pi0
C
C                 INDEX  0-9 : show item 'INDEX' only
C                 INDEX 11-19: show items 1->last digit of 'INDEX'
C                 INDEX 21-29: show items 0->last digit of 'INDEX'
C                 KWRU directs the output to selected (KWRU) file number
C     ------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: INDEX
      real*8, dimension(10) ::       ! no longer needed except in writer
C!nu               TROPOSPHERIC AEROSOL effective radius
C!nu               BCI  OCI  SUI  SEA  SUN    ANT  OCN  OCB  BCB  SSB
     *   REAERO=(/ 0.1, 0.3, 0.3, 2.0, 0.3,   1.0, 0.3, 0.3, 0.2, 0.5/)
      CHARACTER*8, PARAMETER :: FTYPE(5) =
     *     (/'DOWNWARD','  UPWARD','UPWD NET','COOLRATE','FRACTION'/)
      CHARACTER*6, PARAMETER :: GHG(12) =
     *     (/'   H2O','   CO2','    O3','    O2','   NO2','   N2O'
     +      ,'   CH4','CCL3P1','CCL2P2','    N2',' CFC-Y',' CFC-Z'/)
      CHARACTER*3 TRABCD(5),TRAXSG(5),snotyp
      DATA TRABCD/'TRA','TRB','TRC','TRD','TRE'/
      DATA TRAXSG/'QAB','QEX','QSC','QCB','PI0'/

      REAL*8 TKEFF(3),TRPI0K(25)

      REAL*8, DIMENSION(33) :: BGFLUX,BGFRAC,TAUSUM
      REAL*8 SUM0(20),SUM1(LX),SUM2(LX),SUM3(LX)
      REAL*8, DIMENSION(LX,6) :: WSREXT,WSRSCT,WSRGCB,WSRPI0
      REAL*8, DIMENSION(17) :: FSR1,FSR2
      INTEGER :: ISR1(16), KWRU
      INTEGER, PARAMETER, DIMENSION(16) ::
     *     KSLAMW=(/1,1,2,2,5,5,5,5,1,1,1,3,4,6,6,1/),
     *     IORDER=(/12,11,10, 9, 6, 5, 4, 3,15,14,13, 8, 7, 2, 1,16/)

      character*1,parameter,dimension(4) :: AUXGAS = (/'0','L','X','X'/)
      REAL*8, PARAMETER :: P0=1013.25,SIGMA=5.6697D-08,
     *     xxxxxx=0. ! dummy for obsolete variables
      REAL*8 ACOLX,BCOLX,DCOLX,VCOLX,TCOLX,FACTOR,PPMCO2,PPMO2
     *     ,PPMN2O,PPMCH4,PPMF11,PPMF12,PPMY11,PPMZ12,EPS,TAER,HLM,TLAPS
     *     ,TAU55,TGMEAN,PSUM,SRALB,STNFLB,CRHRF,STFHR,TRDCR,SRDHR,STDHR
     *     ,PFW,DPF,FRACSL,SUM,SIGT4,WTG,SUMK,SUMT,SUMK1,SUMK2
     *     ,ASUM1,BSUM1,CSUM1,DSUM1,ESUM1,FSUM1,ASUM2,BSUM2,CSUM2,DSUM2
     *     ,ESUM2,FSUM2,ASUM3,BSUM3,CSUM3,DSUM3,ESUM3,FSUM3,SUML,SUMA
     *     ,SUMB,SUMC,SUMD,SUME,SUMF
      INTEGER I,J,K,L,KW,INDJ,INDI,INDX,KPAGE,NPAGE,LUXGAS,LGS,IPI0,IRHL
     *     ,N,II,IPF,ITG,LK,KK,NW,LINFIL

      DO 20 K=1,6
      DO 10 L=1,NL
      WSREXT(L,K)=SRAEXT(L,K)+SRBEXT(L,K)
     +           +SRDEXT(L,K)+SRVEXT(L,K)
      WSRSCT(L,K)=SRASCT(L,K)+SRBSCT(L,K)
     +           +SRDSCT(L,K)+SRVSCT(L,K)
      WSRGCB(L,K)=SRASCT(L,K)*SRAGCB(L,K)+SRBSCT(L,K)*SRBGCB(L,K)
     +           +SRDSCT(L,K)*SRDGCB(L,K)+SRVSCT(L,K)*SRVGCB(L,K)
      WSRPI0(L,K)=WSRSCT(L,K)/(WSREXT(L,K)+1.E-10)
      WSRGCB(L,K)=WSRGCB(L,K)/(WSRSCT(L,K)+1.D-10)
   10 CONTINUE
   20 CONTINUE
C
      ACOLX=0.0
      BCOLX=0.0
      DCOLX=0.0
      VCOLX=0.0
      DO 30 L=1,NL
      ACOLX=ACOLX+SRAEXT(L,6)
      BCOLX=BCOLX+SRBEXT(L,6)
      DCOLX=DCOLX+SRDEXT(L,6)
      VCOLX=VCOLX+SRVEXT(L,6)
   30 CONTINUE
      TCOLX=ACOLX+BCOLX+DCOLX+VCOLX
C
      KW=KWRU
      INDJ=MOD(INDEX,10)
      IF(INDJ.LT.1.and.INDEX.gt.0) INDJ=10
      INDI=1
      IF(INDEX.GT.20.or.INDEX.eq.0) INDI=0
      IF(INDEX.LT.11) INDI=INDJ
      DO 9999 INDX=INDI,INDJ
C
      KPAGE=1
      IF(INDX.EQ.0) GO TO 90
      GO TO (100,200,300,400,500,600,700,800,900,1000),INDX
C
C-------------
   90 CONTINUE
C-------------
C
      WRITE(KW,6000)
 6000 FORMAT(' CALL WRITER(KW,0) :',2X,'PAGE 1/2  '
     +          ,'CONTROL PARAMS   DEFINITIONS'/
     +      /' CONTROL PARAMTER      DEFAULT  PARAMETER DESCRIPTION')
      WRITE(KW,6001)                              KUVFAC,KSNORM
     + ,KWTRAB,KGGVDF,KPGRAD,KLATZ0,KCLDEM,KANORM,KPFCO2,KPFOZO,KSIALB
     + ,KORDER,KUFH2O,KUFCO2,KCSELF,KCFORN
 6001 FORMAT( ! 7X,'   KVRAER = ',I1,'     1      Repartition Aer VDist'
!nu  2    ! /7X,'   MEANAC = ',I1,'     0      Use Ann-Mean Aer Clim'
!nu  3    ! /7X,'   MEANDD = ',I1,'     0      Use Ann-Mean Des Dust'
!nu  4    ! /7X,'   MEANVA = ',I1,'     0      Use Ann-Mean Volc Aer'
!nu  5      /7X,'   NCARO3 = ',I1,'     0      NCAR London 1976 Ozon'/
     6       7X,'   KUVFAC = ',I1,'     0      ON/OFF UV Mult Factor'
     7      /7X,'   KSNORM = ',I1,'     0      Norm S0 when KUVFAC=1'
     8      /7X,'   KWTRAB = ',I1,'     0      WRITER: Qab,Qex,Qsc,g'
     9      /7X,'   KGGVDF = ',I1,'     0      Use GHG VertProf Grad'
     A      /7X,'   KPGRAD = ',I1,'     1      Pole-to-Pole GHG Grad'
     1      /7X,'   KLATZ0 = ',I1,'     1      Use GHG VDist Lat Dep'
     2      /7X,'   KCLDEM = ',I1,'     1      Use TopCloud Scat Cor'
     3      /7X,'   KANORM = ',I1,'     0      Use SGP Atmo Col Norm'
     4      /7X,'   KPFCO2 = ',I1,'     0      1=MOD CO2PROF: FPXCO2'
     5      /7X,'   KPFOZO = ',I1,'     0      1=MOD O3 PROF: FPXOZO'
     6      /7X,'   KSIALB = ',I1,'     0      Schramm"s ocn ice alb'
     7      /7X,'   KORDER = ',I1,'     0      WRITER k-d spec order'
     8      /7X,'   KUFH2O = ',I1,'     1      Col Absorber Scal H2O'
     +      /7X,'   KUFCO2 = ',I1,'     1      Col Absorber Scal CO2'
     1      /7X,'   KCSELF = ',I1,'     1      H2O Cont Self-Broaden'
     2      /7X,'   KCFORN = ',I1,'     1      H2O Con Foreign-Broad'
     R      )
C
      WRITE(KW,6004)
 6004 FORMAT(/' CONTROL PARAMTER    DEFAULT       SNOW/ICE FACTORS')

      WRITE(KW,6005) agexpf,albdif
 6005 FORMAT(7X   ,'   AGEXPF = ',F7.3,'      SNOWAGE XPFACTOR SH EARTH'
     A      /7X   ,'    SH O  = ',F7.3,'         "       "     "  OCICE'
     B      /7X   ,'    SH L  = ',F7.3,'         "       "     "  LDICE'
     C      /7X   ,'    NH E  = ',F7.3,'         "       "     NH EARTH'
     D      /7X   ,'    NH O  = ',F7.3,'         "       "     "  OCICE'
     E      /7X   ,'    NH L  = ',F7.3,'         "       "     "  LDICE'
     F      /7X   ,'   ALBDIF = ',F7.3,'      SNOW/ICE ALBDIF  SH EARTH'
     G      /7X   ,'    SH O  = ',F7.3,'         "       "     "  OCICE'
     H      /7X   ,'    SH L  = ',F7.3,'         "       "     "  LDICE'
     I      /7X   ,'    NH E  = ',F7.3,'         "       "     NH EARTH'
     J      /7X   ,'    NH O  = ',F7.3,'         "       "     "  OCICE'
     K      /7X   ,'    NH L  = ',F7.3,'         "       "     "  LDICE'
     L      )
C
      WRITE(KW,6006)
 6006 FORMAT('0CONTROL PARAMTER    VALUE',16X,' DEFAULT')
      WRITE(KW,6007) PTOPTR,REFF0 ,VEFF0
     +              ,AVSCAT,ANSCAT,AVFOAM,ANFOAM
 6007 FORMAT(7X,'   PTOPTR = ',F7.1,6X,'TropTop (SIGMA lev) Pressure'
     A      /7X,'   REFF0  = ',F7.3,'                 0.300         '
     B      /7X,'   VEFF0  = ',F7.3,'                 0.350         '
     H      /7X,'   AVSCAT = ',F7.5,'                 0.01560       '
     I      /7X,'   ANSCAT = ',F7.5,'                 0.00020       '
     J      /7X,'   AVFOAM = ',F7.5,'                 0.21970       '
     K      /7X,'   ANFOAM = ',F7.5,'                 0.15140       '
     X      )
      WRITE(KW,6008)
 6008 FORMAT(/10X,'UV Solar Flux Spectral Partitions and Factors')
      WRITE(KW,6009) UVWAVL,UVFACT
 6009 FORMAT(10X,'UVWAVL = ',F7.5,2F8.5/10X,'UVFACT = ',F7.5,2F8.5)
C
!nu   WRITE(KW,6013)
 6013 FORMAT(/' CONTROL PARAMETER  PI0VIS     PI0TRA      DEFAULT')
!nu   WRITE(KW,6014) PI0VIS,PI0TRA
 6014 FORMAT(7X,'   ACID1 = ',F8.6,F11.6,'     1.0                 '
     A      /7X,'   SSALT = ',F8.6,F11.6,'     1.0                 '
     B      /7X,'   SLFT1 = ',F8.6,F11.6,'     1.0                 '
     C      /7X,'   SLFT2 = ',F8.6,F11.6,'     1.0                 '
     D      /7X,'   BSLT1 = ',F8.6,F11.6,'      .98929             '
     E      /7X,'   BSLT2 = ',F8.6,F11.6,'      .95609             '
     F      /7X,'   DUST1 = ',F8.6,F11.6,'      .91995             '
     G      /7X,'   DUST2 = ',F8.6,F11.6,'      .78495             '
     H      /7X,'   DUST3 = ',F8.6,F11.6,'      .63576             '
     I      /7X,'   CARB1 = ',F8.6,F11.6,'      .31482             '
     J      /7X,'   CARB2 = ',F8.6,F11.6,'      .47513             '
     K      )

      WRITE(KW,6019)
 6019 FORMAT(/'  GHGAS',9X,'PPMVK0    PPMVDF    PPGRAD')
      WRITE(KW,6020) (ghg(I),PPMVK0(I),PPMVDF(I),PPGRAD(I),I=1,12)
 6020 FORMAT(1X,a6,' ',F15.7,F10.5,F10.5)
      GO TO 9999
C
C-------------
  100 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6101)
      WRITE(KW,6102)
      FACTOR=P0/(PLB(1)-PLB(2))*1.25
      PPMCO2=ULGAS(1,2)*FACTOR
      PPMO2 =ULGAS(1,4)*FACTOR
      PPMN2O=ULGAS(1,6)*FACTOR
      PPMCH4=ULGAS(1,7)*FACTOR
      PPMF11=ULGAS(1,8)*FACTOR
      PPMF12=ULGAS(1,9)*FACTOR
      PPMY11=ULGAS(1,11)*FACTOR
      PPMZ12=ULGAS(1,12)*FACTOR
      WRITE(KW,6103) (FULGAS(I),I=1,3),(FULGAS(I),I=6,9)
     +              ,FULGAS(11),FULGAS(12),    (FGOLDH(I),I=1,5)
C     IF(KGASSR.GT.0)
C    +WRITE(KW,6104) (FULGAS(I+9),I=1,2),(FULGAS(I+9),I=4,9)
C    +              ,FULGAS(11),FULGAS(12),   (FGOLDH(I+9),I=1,5)
      WRITE(KW,6105) PPMCO2,PPMN2O,PPMCH4,PPMF11,PPMF12,PPMY11,PPMZ12
     +             ,(FSTOPX(I),I=1,4),PPMV80(2),(PPMV80(I),I=6,9)
     +             ,(PPMV80(I),I=11,12),KTREND,JYEAR,JDAY,LASTVC
      WRITE(KW,6106) TAUWC0,FCLDTR,EOCTRA,ZOCSRA,KZSNOW,KCLDEM,NTRACE
     +             ,FSAAER,FTTAER,MADO3M,KCLDEP,NL
      WRITE(KW,6107) TAUIC0,FCLDSR,ESNTRA,ZSNSRA,WETTRA,KSIALB,ITR(1)
     +             ,ITR(5),FSBAER,FTBAER,MADO3M,KEEPAL,NLP
      WRITE(KW,6108)       FRAYLE,EICTRA,ZICSRA,WETSRA,KCNORM,ITR(2)
     +             ,ITR(6),FSAAER,FTAAER,       KEEP10,MLAT46
      WRITE(KW,6109) TLGRAD,ECLTRA,EDSTRA,ZDSSRA,KANORM,KPGRAD,ITR(3)
     +             ,ITR(7),FSDAER,FTDAER,KWVCON,ICE012,MLON72
      WRITE(KW,6110) PTLISO       ,EVGTRA,ZVGSRA,KEEPRH,KLATZ0,ITR(4)
     +             ,ITR(8),FSVAER,FTVAER,KSOLAR,NORMS0
C
 6101 FORMAT(' (1)FUL:  1',7X,'2',8X,'3',7X,'6',7X,'7',8X,'8',8X,'9'
     +      ,8X,'11',7X,'12',4X,'RADPAR 1/F: (Control/Default'
     +      ,'/Scaling Parameters)')
 6102 FORMAT(4X,'GAS: ','H2O',5X,'CO2',7X,'O3'
     +      ,5X,'N2O',5X,'CH4',5X,'CFC-11',3X,'CFC-12'
     +      ,3X,'CFY-11',3X,'CFZ-12'
     +      ,2X,'Aerosol Global    Ocean     Land  Desert    Haze')
 6103 FORMAT(1X,'FULGAS=',F5.3,F10.5,F7.3,F9.5,F8.5,4F9.5
     +      ,2X,'FGOLDH=',F7.5,2F9.6,2F8.5)
C6104 FORMAT('+',T84,'T'
C    +      /1X,'FULGAS=',1P,1E7.1,1P,2E8.1,1P,2E8.1,1P,4E9.1
C    +      ,' S','FGOLDH=',1P,1E7.1,1P,2E9.2,1P,2E8.1)
 6105 FORMAT(1X,'PPM(1)=(now)',2X,F8.3,8X,F8.5,F8.5,4(1X,F8.7)
     +      ,2X,'TRACER=',F7.5,2F9.6,F8.5
     +      /' PPMV80=(ref)=',0P,F9.3,8X,2F8.5,4(1X,F8.7),2X
     +      ,'KTREND=',I1,2X,'JYEAR=',I4,' JDAY=',I3,5X,'LASTVC=',I7)
 6106 FORMAT(1X,'TAUWC0=',1P,E6.0,' FCLDTR=',0P,F4.2,' EOCTRA=',F3.1
     +      ,1X,'ZOCSRA=',   F3.1,' KZSNOW=',     I4,' KCLDEM=',  I3
     +      ,1X,'NTRACE=',    I3,2X,'FSTAER=',  F3.1,' FTTAER=',F3.1
     +      ,1X,'MADO3M=',    I1,1X,'KCLDEP=',    I1,'     NL=',  I2)
 6107 FORMAT(1X,'TAUIC0=',1P,E6.0,' FCLDSR=',0P,F4.2,' ESNTRA=',F3.1
     +      ,1X,'ZSNSRA=',  F3.1,1X,'WETTRA=',  F4.2,' KSIALB=',  I3
     +      ,1X,'ITR(1)=',    2I2,1X,'FSBAER=', F3.1,' FTBAER=',F3.1
     +      ,1X,'K03LON=',    I1,1X,'KEEPAL=',    I1,'    NLP=',  I2)
 6108 FORMAT(1X,'       ',  6X ,  ' FRAYLE=',0P,F4.1,' EICTRA=',F3.1
     +      ,1X,'ZICSRA=',  F3.1,1X,'WETSRA=',  F4.2,' KCNORM=',  I3
     +      ,1X,'ITR(2)=',    2I2,1X,'FSAAER=', F3.1,' FTAAER=',F3.1
     +      ,1X,'       ',   ' ',1X,'KEEP10=',    I1,' MLAT46=',  I2)
 6109 FORMAT(1X,'TLGRAD=',   F6.2,' ECLTRA=',0P,F4.2,' EDSTRA=',F3.1
     +      ,1X,'ZDSSRA=',  F3.1,1X,'KANORM=',  I4  ,' KPGRAD=',  I3
     +      ,1X,'ITR(3)=',    2I2,1X,'FSDAER=', F3.1,' FTDAER=',F3.1
     +      ,1X,'KWVCON=',    I1,1X,'ICE012=',    I1,' MLON72=',  I2)
 6110 FORMAT(1X,'PTLISO=',  F6.1,1X,'           '   ,' EVGTRA=',F3.1
     +      ,1X,'ZVGSRA=',  F3.1,1X,'KEEPRH=',  I4  ,' KLATZ0=',  I3
     +      ,1X,'ITR(4)=',    2I2,1X,'FSVAER=', F3.1,' FTVAER=',F3.1
     +      ,1X,'KSOLAR=',    I1,1X,'NORMS0=',    I1,'        ')
      GO TO 9999
C
C-------------
  200 CONTINUE
C-------------
C
      NPAGE=0
      LUXGAS=0
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6201) AUXGAS(LUXGAS+1),S00WM2,S0,COSZ
      DO 202 K=1,9
      DO 201 L=1,NL
      UXGAS(L,K)=ULGAS(L,K)
  201 CONTINUE
  202 CONTINUE
      IF(LUXGAS.LT.2) GO TO 205
      LGS=(LUXGAS-2)*9
      DO 203 L=1,NL
      UXGAS(L,1)=U0GAS(L,1)*FULGAS(1+LGS)
      UXGAS(L,3)=U0GAS(L,3)*FULGAS(3+LGS)
  203 UXGAS(L,5)=U0GAS(L,5)*FULGAS(5+LGS)
C
      DO 204 L=1,NL
      UXGAS(L,2)=U0GAS(L,2)*FULGAS(2+LGS)
      UXGAS(L,4)=U0GAS(L,4)*FULGAS(4+LGS)
      UXGAS(L,6)=U0GAS(L,6)*FULGAS(6+LGS)
      UXGAS(L,7)=U0GAS(L,7)*FULGAS(7+LGS)
      UXGAS(L,8)=U0GAS(L,8)*FULGAS(8+LGS)
  204 UXGAS(L,9)=U0GAS(L,9)*FULGAS(9+LGS)
  205 CONTINUE
      DO 206 L=NL,1,-1
      EPS=CLDEPS(L)
      TAER=WSREXT(L,6)
      IPI0=WSRPI0(L,6)*1000.D0+1.D-05
      HLM=0.5D0*(HLB0(L+1)+HLB0(L))
      TLAPS=(TLT(L)-TLB(L))/(HLB0(L+1)-HLB0(L))
      IRHL=RHL(L)*100.0
      IF(PL(L).LT.1.D0) THEN
      WRITE(KW,6212) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +     ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +     ,SIZEWC(L),SIZEIC(L),TAUWC(L),TAUIC(L),EPS,TAER,IPI0
      ELSE
      IF(UXGAS(L,1).GE.1.D0) THEN
      WRITE(KW,6202) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +     ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +     ,SIZEWC(L),SIZEIC(L),TAUWC(L),TAUIC(L),EPS,TAER,IPI0
      ELSE
      WRITE(KW,6211) L,PL(L),HLM,TLM(L),TLAPS,SHL(L),IRHL
     +     ,(UXGAS(L,K),K=1,3),(UXGAS(L,K),K=6,9),UXGAS(L,5)
     +     ,SIZEWC(L),SIZEIC(L),TAUWC(L),TAUIC(L),EPS,TAER,IPI0
      ENDIF
      ENDIF
  206 CONTINUE
      DO 207 I=1,16
  207 SUM0(I)=0.
      DO 210 L=1,NL
      DO 208 I=1,9
      SUM0(I)=SUM0(I)+ULGAS(L,I)
  208 CONTINUE
      DO 209 I=1,4
      SUM0(12+I)=SUM0(12+I)+TRACER(L,I)
  209 CONTINUE
      SUM0(10)=SUM0(10)+TAUWC(L)
      SUM0(11)=SUM0(11)+TAUIC(L)
  210 CONTINUE
      TAU55=0.0
      DO 211 L=1,NL
      TAU55=TAU55+WSREXT(L,6)
  211 CONTINUE
      SUM0(12)=TAU55
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      WRITE(KW,6203) (SUM0(I),I=1,3),(SUM0(I),I=6,9),SUM0(5)
     +               ,SUM0(10),SUM0(11),SUM0(12)
      WRITE(KW,6204) POCEAN,TGO,PLAKE,zlake,SUM0(13),JYEAR
     +             ,BXA(4:5),LASTVC
      WRITE(KW,6205) PEARTH,TGE,SNOWE,ZSNWOI,SUM0(14),JDAY,BXA(6:7)
      WRITE(KW,6206) POICE,TGOI,SNOWOI,ZOICE,SUM0(15),JLAT
     +             ,(SRBALB(I),I=1,6)
      WRITE(KW,6207) PLICE,TGLI,SNOWLI,zmp,SUM0(16),ILON
     +             ,(SRXALB(I),I=1,6)
      PSUM=POCEAN+PEARTH+POICE+PLICE
      snotyp='DRY' ; if(flags) snotyp='WET'
      WRITE(KW,6208) TGMEAN,snotyp,fmp
     +               ,PSUM,TSL,WMAG,LS1_loc,(PVT(I),I=1,11)
      write(kw,6213) snow_frac(1),snow_frac(2),agesn(1),
     +  agesn(2),agesn(3),wearth,fulgas(4),fulgas(5),fulgas(10)
 6213 FORMAT(1X,'FSNWds=',F6.4,' FSNWvg=',F6.4,'  AGESN=[EA:',F6.3,
     +      ' OI:',F6.3,' LI:',F6.3,'] WEARTH=',F6.4,1X,
     +      ' FULGAS[ 4=O2:',F3.1,' 5=NO2:',F3.1,' 10=N2C:',F3.1,']')
      WRITE(KW,6209) (PRNB(1:2,I),PRNX(1:2,I),I=1,4),BXA(1:3)
      WRITE(KW,6210)
 6201 FORMAT(' (2) RADPAR G/L: (Input Data)'
     +      ,2X,'Absorber Amount per Layer:'
     +      ,'  U',1A1,'GAS(L,K) in cm**3(STP)/cm**2',2X,'S00WM2=',F9.4
     +      ,1X,'S0=',F9.4,2X,'COSZ=',F6.4
     +      /' LN     PL  HLM   TLM  TLAP   SHL  .RH     '
     +      ,'H2O   CO2    O3    N2O    CH4  CFC-11'
     +      ,'  CFC-12   NO2   WC.SIZ.IC  WC.TAU.IC CLEP  A TAU PI0')
 6202 FORMAT(1X,I2,F7.2,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.2,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
 6203 FORMAT(24X,' Column Amount',F7.1,F7.2,1X,F6.5
     +       ,1X,F5.4,F7.4,1P,3E8.1,0P,10X,F6.2,F5.2,5X,F6.3)
 6204 FORMAT( 1X,'PWATER=',F6.4,'    TGO=' ,F6.2,1X,' PLAKE=',F6.3
     +      , 1X,' ZLAKE=',F6.3,' TRACER 1=',F5.3,' JYEAR=',I4
     +      , 3X,'BSNVIS=',F6.4,' BSNNIR=' ,F6.4,7X,'LASTVC=',I7)
 6205 FORMAT(    ' PEARTH=',F6.4,'    TGE=',F6.2,'  SNOWE=',F6.3
     +      ,    '  ZSNOW=',F6.3,'  Sums: 2=',F5.3
     +      ,     '  JDAY=',I4  ,2X,' XSNVIS=',F6.4,' XSNNIR=',F6.4
     +      , 8X,'NIRALB VISALB')
 6206 FORMAT(    '  POICE=',F6.4,'   TGOI=',F6.2,' SNOWOI=',F6.3
     +      ,    '  ZOICE=',F6.3,'        3=',F5.3
     +      ,     '  JLAT=',I4,  2X,' SRBALB=',F6.4
     +      ,4F7.4,F7.4)
 6207 FORMAT(    '  PLICE=',F6.4,'   TGLI=',F6.2,' SNOWLI=',F6.3
     +      ,    '  ZMLTP=',F6.3,'        4=',F5.3
     +      ,     '  ILON=',I4,  2X,' SRXALB=',F6.4
     +      ,4F7.4,F7.4)
 6208 FORMAT(8X,6('-'),' TGMEAN=',F6.2,'    SNOW : ',a3,'  FMLTP='
     +      ,F6.3,'  BSAND TUNDRA GRASSL SHRUBS  TREES DECIDF'
     +       ,' EVERGF','  RAINF','  CROPS','  BDIRT','  ALGAE'
     +      /    '   PSUM=',F6.4,'    TSL=',F6.2,' WINDSP=',F6.3
     +       ,'  LS1L=',I2,T54,'PVT=',F6.4,10F7.4)
 6209 FORMAT(' BOCVIS BOCNIR XOCVIS XOCNIR BEAVIS BEANIR XEAVIS XEANIR'
     +      ,' BOIVIS BOINIR XOIVIS XOINIR BLIVIS BLINIR XLIVIS XLINIR'
     +      ,' EXPSNE EXPSNO EXPSNL'/1X,F6.4,18F7.4)
 6210 FORMAT(' ')
 6211 FORMAT(1X,I2,F7.2,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.5,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
 6212 FORMAT(1X,I2,F7.4,F5.1,F7.2,F5.1,1X,F7.6,I3
     +      ,F8.5,F6.2,1X,F6.5,1X,F5.4
     +      ,F7.4,1P,3E8.1,0P,2F5.1,F6.2,F5.2,1X,F4.3,F6.3,I5)
C
      GO TO 9999
C
C-------------
  300 CONTINUE
C-------------
C
      NPAGE=0
      IF(INDEX.LT.11) NPAGE=KPAGE
      IF(NL.GT.13) NPAGE=1
      L=NLP
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)
      STNFLB=SRNFLB(L)-TRNFLB(L)
      WRITE(KW,6301) NORMS0
      WRITE(KW,6302) L,PLB(L),HLB0(L),TLT(L-1) ! TLB(LN+1) unused/set
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),STNFLB,SRALB
      DO 301 N=1,NL
      L=NLP-N
      CRHRF=8.4167/(PLB(L)-PLB(L+1))
      STNFLB=SRNFLB(L)-TRNFLB(L)
      STFHR =SRFHRL(L)-TRFCRL(L)
      TRDCR =TRFCRL(L)*CRHRF
      SRDHR =SRFHRL(L)*CRHRF
      STDHR=STFHR*CRHRF
      SRALB =SRUFLB(L)/(SRDFLB(L)+1.E-10)
!eq   SRXVIS=SRXATM(1)
!eq   SRXNIR=SRXATM(2)
      IF(PLB(L).LT.1.D0) THEN
      WRITE(KW,6313) L,PLB(L),HLB0(L),TLB(L),TLT(L)
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB
      ELSE
      WRITE(KW,6303) L,PLB(L),HLB0(L),TLB(L),TLT(L)
     +             ,TRDFLB(L),TRUFLB(L),TRNFLB(L),TRFCRL(L)
     +             ,SRDFLB(L),SRUFLB(L),SRNFLB(L),SRFHRL(L)
     +             ,STNFLB,STFHR,STDHR,TRDCR,SRDHR,SRALB
      ENDIF
  301 CONTINUE
C
      DO 302 II=1,3
      PFW=TRDFLB(1)
      IF(II.EQ.2) PFW=TRUFLB(1)
      IF(II.EQ.3) PFW=TRUFLB(NLP)
      IPF=PFW
      DPF=PFW-IPF
      IF(IPF.LT.1) IPF=1
      IF(IPF.GT.899) IPF=899
  302 TKEFF(II)=TKPFT(IPF)+DPF*(TKPFT(IPF+1)-TKPFT(IPF))
C
      WRITE(KW,6304) WINDZF(1),WINDZT(1),TOTLZF(1),TOTLZT(1),
     +              (FSRNFG(I),I=1,4),LTOPCL,JLAT,JYEAR
      WRITE(KW,6305) WINDZF(2),WINDZT(2),TOTLZF(2),TOTLZT(2),
     +              (FTRUFG(I),I=1,4),LBOTCL,ILON,JDAY
      IF(KORDER.EQ.0) WRITE(KW,6306) WINDZF(3),WINDZT(3),TOTLZF(3)
     +              ,TOTLZT(3),(I,I=1,16)
      IF(KORDER.EQ.1) WRITE(KW,6307) WINDZF(3),WINDZT(3),TOTLZF(3)
     +              ,TOTLZT(3),(I,I=1,16)
      FRACSL=0.D0
      IF(KORDER.EQ.0) WRITE(KW,6308) TKEFF(1),TKEFF(2),TKEFF(3)
     +               ,(SRKALB(NORDER(I)),I=1,16),BTEMPW,TRUFTW,SRIVIS
     +               ,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      IF(KORDER.EQ.1) WRITE(KW,6308) TKEFF(1),TKEFF(2),TKEFF(3)
     +               ,(SRKALB(I),I=1,16),BTEMPW,TRUFTW,SRIVIS
     +               ,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      WRITE(KW,6309) TRDFGW,TRUFGW,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR
     +             ,ALBNIR
      WRITE(KW,6310) SRXVIS,SRXNIR,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR
     +             ,SRANIR
C
C
 6301 FORMAT(/' (3) RADPAR M/S: (Output Data)'
     +      ,T37,'Thermal Fluxes (W/M**2)',4X,'Solar Fluxes (W/M**2)'
     +      ,1X,'NORMS0=',I1,'  Energy  Input  Heat/Cool Deg/Day Alb'
     +      ,'do'/' LN     PLB   HLB    TLB    TLT '
     +      ,'  TRDFLB TRUFLB TRNFLB TRFCRL   SRDFLB  SRUFLB  SRNFLB'
     +      ,' SRFHRL  STNFLB  STFHR  SR-TR TR=CR SR=HR SRALB')
 6302 FORMAT(1X,I2,F9.3,F6.2,1X,F6.2,8X,3F7.2,8X,3F8.2,7X,F8.2,26X,F6.4)
 6303 FORMAT(1X,I2,F9.3,F6.2,2F7.2,1X,3F7.2,F7.2,1X,3F8.2,F7.2,1X,F7.2
     +      ,1X,F6.2,1X,3F6.2,1X,F5.4)
 6304 FORMAT( 1X,'XMU  WINDZF WINDZT   TOTLZF TOTLZT'/
     +        1X,'1.0',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'FR.SRNLB1'
     +      ,' OCEAN=',F7.2,' EARTH=',F7.2,'  OICE=',F7.2,'   LICE='
     +      ,F7.2,1X,' LTOPCL=',I2,' JLAT=',I2,' JYEAR=',I4)
 6305 FORMAT( 1X,'0.5',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'FR.TRULB1'
     +      ,' OCEAN=',F7.4,' EARTH=',F7.4,'  OICE=',F7.4,'   LICE='
     +      ,F7.4,1X,' LBOTCL=',I2,' ILON=',I2,'  JDAY=',I4)
 6306 FORMAT( 1X,'0.1',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'L=',I3,15I6)
 6307 FORMAT( 1X,'0.1',1X,F7.3,F7.2,2X,F7.3,F7.2,2X,'K=',I3,15I6)
 6308 FORMAT(' TKeff= ',F6.2,2F7.2,'  SRKALB=',16F6.4/
     +        1X,'At Top of Atm: ',' BTEMPW=',F6.2,1X,' TRUFTW=',F6.3
     +      , 2X,' SRIVIS=',F6.2,' SROVIS=',F6.2,   ' PLAVIS=',F6.4
     +      , 2X,' SRINIR=',F6.2,' SRONIR=',F6.2,   ' PLANIR=',F6.4)
 6309 FORMAT( 1X,'At Bot of Atm: ',' TRDFGW=',F6.3,1X,' TRUFGW=',F6.3
     +      , 2X,' SRDVIS=',F6.2,' SRUVIS=',F6.2,   ' ALBVIS=',F6.4
     +      , 2X,' SRDNIR=',F6.2,' SRUNIR=',F6.2,   ' ALBNIR=',F6.4)
 6310 FORMAT( 1X,'In Atmosphere: ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4
     +      , 2X,' SRTVIS=',F6.4,' SRRVIS=',F6.4,   ' SRAVIS=',F6.4
     +      , 2X,' SRTNIR=',F6.4,' SRRNIR=',F6.4,   ' SRANIR=',F6.4)
 6311 FORMAT(' ')
 6313 FORMAT(1X,I2,F9.5,F6.2,2F7.2,1X,F7.4,2F7.2,F7.4,1X,3F8.2,F7.4
     +      ,1X,F7.2,F7.4,1X,3F6.2,1X,F5.4)
      GO TO 9999
C
C-------------
  400 CONTINUE
C-------------
C
C                                (4A)  Total Aerosol Qx, Qs, g, Pi0
C                                ----------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6401)
      DO 402 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 401 L=1,NL
      SUM1(K)=SUM1(K)+WSREXT(L,K)
      SUM2(K)=SUM2(K)+WSRSCT(L,K)
      SUM3(K)=SUM3(K)+WSRSCT(L,K)*WSRGCB(L,K)
  401 CONTINUE
      SUM3(K)=SUM3(K)/(SUM2(K)+1.D-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.D-10)
  402 CONTINUE
      WRITE(KW,6402) (K,K=1,6),(K,K=1,6)
      DO 403 N=1,NL
      L=NLP-N
      WRITE(KW,6403) L,PLB(L),HLB0(L)
     +              ,(WSREXT(L,J),J=1,6),(WSRSCT(L,J),J=1,6)
  403 CONTINUE
      WRITE(KW,6404) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL.GT.13) NPAGE=1
      WRITE(KW,6405) KANORM
      WRITE(KW,6406) (K,K=1,6),(K,K=1,6)
      DO 404 N=1,NL
      L=NLP-N
      WRITE(KW,6407) L,PL(L),DPL(L)
     +              ,(WSRGCB(L,J),J=1,6),(WSRPI0(L,J),J=1,6)
  404 CONTINUE
      WRITE(KW,6408) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
C     WRITE(KW,6420) (SRBALB(K),K=1,6)
C     WRITE(KW,6421) (SRXALB(K),K=1,6)
C     WRITE(KW,6422)
      SUM=0.
      DO 406 J=1,5
      TAU55=0.
      DO 405 I=1,11  !  NAERO
  405 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6423) J,FGOLDH(J),TAU55
  406 SUM=SUM+TAU55
      WRITE(KW,6424) SUM
      WRITE(KW,6424) BCOLX,ACOLX,DCOLX,VCOLX,TCOLX
      DO 407 I=1,8
      WRITE(KW,6425)
  407 CONTINUE
C
C                                (4B)  Water/Ice Cloud Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6411)
      DO 412 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 411 L=1,NL
      SUM1(K)=SUM1(K)+SRCEXT(L,K)
      SUM2(K)=SUM2(K)+SRCSCT(L,K)
      SUM3(K)=SUM3(K)+SRCSCT(L,K)*SRCGCB(L,K)
  411 SRCPI0(L,K)=SRCSCT(L,K)/(SRCEXT(L,K)+1.D-10)
      SUM3(K)=SUM3(K)/(SUM2(K)+1.D-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.D-10)
  412 CONTINUE
      WRITE(KW,6412) (K,K=1,6),(K,K=1,6)
      DO 413 N=1,NL
      L=NLP-N
      WRITE(KW,6413) L,PLB(L),HLB0(L)
     +              ,(SRCEXT(L,J),J=1,6),(SRCSCT(L,J),J=1,6)
  413 CONTINUE
      WRITE(KW,6414) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL.GT.13) NPAGE=1
      WRITE(KW,6415) KANORM
      WRITE(KW,6416) (K,K=1,6),(K,K=1,6)
      DO 414 N=1,NL
      L=NLP-N
      WRITE(KW,6417) L,PL(L),DPL(L)
     +              ,(SRCGCB(L,J),J=1,6),(SRCPI0(L,J),J=1,6)
  414 CONTINUE
      WRITE(KW,6418) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
      WRITE(KW,6420) (SRBALB(K),K=1,6)
      WRITE(KW,6421) (SRXALB(K),K=1,6)
      WRITE(KW,6422)
      SUM=0.
      DO 416 J=1,5
      TAU55=0.
      DO 415 I=1,11  !  NAERO
  415 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6423) J,FGOLDH(J),TAU55
  416 SUM=SUM+TAU55
      WRITE(KW,6424) SUM
      DO 417 I=1,2
      WRITE(KW,6425)
  417 CONTINUE
C
C                                (4C)  Aerosol + Cloud Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6426)
      DO 419 K=1,6
      SUM1(K)=0.
      SUM2(K)=0.
      SUM3(K)=0.
      DO 418 L=1,NL
      SUM1(K)=SUM1(K)+DBLEXT(L,K)
      SUM2(K)=SUM2(K)+DBLSCT(L,K)
      SUM3(K)=SUM3(K)+DBLSCT(L,K)*DBLGCB(L,K)
  418 DBLPI0(L,K)=DBLSCT(L,K)/(DBLEXT(L,K)+1.E-10)
      SUM3(K)=SUM3(K)/(SUM2(K)+1.E-10)
      SUM0(K)=SUM2(K)/(SUM1(K)+1.E-10)
  419 CONTINUE
      WRITE(KW,6427) (K,K=1,6),(K,K=1,6)
      DO 420 N=1,NL
      L=NLP-N
      WRITE(KW,6428) L,PLB(L),HLB0(L)
     +              ,(DBLEXT(L,J),J=1,6),(DBLSCT(L,J),J=1,6)
  420 CONTINUE
      WRITE(KW,6429) (SUM1(K),K=1,6),(SUM2(K),K=1,6)
      NPAGE=0
      IF(NL.GT.13) NPAGE=1
      WRITE(KW,6430) KANORM
      WRITE(KW,6431) (K,K=1,6),(K,K=1,6)
      DO 421 N=1,NL
      L=NLP-N
      WRITE(KW,6432) L,PL(L),DPL(L)
     +              ,(DBLGCB(L,J),J=1,6),(DBLPI0(L,J),J=1,6)
  421 CONTINUE
      WRITE(KW,6433) (SUM3(K),K=1,6),(SUM0(K),K=1,6)
      WRITE(KW,6434) (SRBALB(K),K=1,6)
      WRITE(KW,6435) (SRXALB(K),K=1,6)
      WRITE(KW,6436)
      SUM=0.
      DO 423 J=1,5
      TAU55=0.
      DO 422 I=1,11  !  NAERO
  422 TAU55=TAU55+AGOLDH(I,J)*FGOLDH(J)
      WRITE(KW,6437) J,FGOLDH(J),TAU55
  423 SUM=SUM+TAU55
      WRITE(KW,6438) SUM
      DO 424 I=1,2
      WRITE(KW,6439)
  424 CONTINUE
C
C                                (4D)  11-Comp Aerosol Qx, Qs, g, Pi0
C                                ------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6440)  KWTRAB,(N,N=1,11)
      WRITE(KW,6441)
      DO 425 K=1,6
      WRITE(KW,6442) K,(SRAQEX(K,N),N=1,11)
  425 CONTINUE
      WRITE(KW,6443)
      DO 426 K=1,6
      WRITE(KW,6442) K,(SRAQSC(K,N),N=1,11)
  426 CONTINUE
      WRITE(KW,6444)
      DO 427 K=1,6
      WRITE(KW,6442) K,(SRAQCB(K,N),N=1,11)
  427 CONTINUE
      WRITE(KW,6445) TRABCD(1),TRAXSG(KWTRAB+1)
      DO 428 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6442) K,(TRAQAB(K,N),N=1,11)
      IF(KWTRAB.EQ.1) WRITE(KW,6442) K,(TRAQEX(K,N),N=1,11)
      IF(KWTRAB.EQ.2) WRITE(KW,6442) K,(TRAQSC(K,N),N=1,11)
      IF(KWTRAB.EQ.3) WRITE(KW,6442) K,(TRAQCB(K,N),N=1,11)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,11
      TRPI0K(N)=TRAQSC(K,N)/(1.D-10+TRAQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,11)
      ENDIF
  428 CONTINUE
      DO 429 I=1,1
      WRITE(KW,6446)
  429 CONTINUE
C
C
C                                (4E)  10-Comp Aerosol Qx, Qs, g, Pi0
C                                ------------------------------------
      WRITE(KW,6450) KWTRAB,(N,N=1, 6),(REFDRY(N),N=1, 6)
      WRITE(KW,6451)
      DO 435 K=1,6
      WRITE(KW,6452) K,(SRHQEX(K,1,N),N=1, 4),(SRBQEX(K,N),N=5, 6)
  435 CONTINUE
      WRITE(KW,6453)
      DO 436 K=1,6
      WRITE(KW,6452) K,(SRHQSC(K,1,N),N=1, 4),(SRBQSC(K,N),N=5, 6)
  436 CONTINUE
      WRITE(KW,6454)
      DO 437 K=1,6
      WRITE(KW,6452) K,(SRHQCB(K,1,N),N=1, 4),(SRBQCB(K,N),N=5, 6)
  437 CONTINUE
      WRITE(KW,6455) TRABCD(2),TRAXSG(1) !obs TRAXSG(KWTRAB+1)
      DO 438 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6442) K,(TRHQAB(K,1,N),N=1, 4),
     *                                 (TRBQAB(K,N),N=5, 6)
!obs  IF(KWTRAB.EQ.1) WRITE(KW,6442) K,(TRHQEX(K,1,N),N=1, 4),
!obs *                                 (TRBQEX(K,N),N=5, 6)
!obs  IF(KWTRAB.EQ.2) WRITE(KW,6442) K,(TRHQSC(K,1,N),N=1, 4),
!obs *                                 (TRBQSC(K,N),N=5, 6)
!obs  IF(KWTRAB.EQ.3) WRITE(KW,6442) K,(TRHQCB(K,1,N),N=1, 4),
!obs *                                 (TRBQCB(K,N),N=5, 6)
!obs  IF(KWTRAB.EQ.4) THEN
!obs  DO N=1,4
!obs  TRPI0K(N)=TRHQCB(K,1,N)/(1.D-10+TRHQEX(K,1,N))
!obs  END DO
!obs  DO N=5,6 ! 10
!obs  TRPI0K(N)=TRBQSC(K,N)/(1.D-10+TRBQEX(K,N))
!obs  END DO
!obs  WRITE(KW,6442) K,(TRPI0K(N),N=1, 6)
!obs  ENDIF
  438 CONTINUE
      DO 439 I=1,1
      WRITE(KW,6456)
  439 CONTINUE
C
C                             (4F  8-size Dust Aerosol Qx, Qs, g, Pi0
C                             ---------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6460) KWTRAB,(N,N=1,8),(REDUST(N),N=1,8)
      WRITE(KW,6461)
      DO 445 K=1,6
      WRITE(KW,6462) K,(SRAQEX(K,N),N=1,8)
  445 CONTINUE
      WRITE(KW,6463)
      DO 446 K=1,6
      WRITE(KW,6462) K,(SRAQSC(K,N),N=1,8)
  446 CONTINUE
      WRITE(KW,6464)
      DO 447 K=1,6
      WRITE(KW,6462) K,(SRAQCB(K,N),N=1,8)
  447 CONTINUE
      WRITE(KW,6465) TRABCD(4),TRAXSG(KWTRAB+1)
      DO 448 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6442) K,(TRDQAB(K,N),N=1,8)
      IF(KWTRAB.EQ.1) WRITE(KW,6442) K,(TRDQEX(K,N),N=1,8)
      IF(KWTRAB.EQ.2) WRITE(KW,6442) K,(TRDQSC(K,N),N=1,8)
      IF(KWTRAB.EQ.3) WRITE(KW,6442) K,(TRDQCB(K,N),N=1,8)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,8
      TRPI0K(N)=TRDQSC(K,N)/(1.D-10+TRDQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,8)
      ENDIF
  448 CONTINUE
      DO 449 I=1,1
      WRITE(KW,6466)
  449 CONTINUE
C
C                           (4G  15-Size/phase Cloud Qx, Qs, g, Pi0
C                           ---------------------------------------
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
      WRITE(KW,6470) KWTRAB,(N,N=1,15)
      WRITE(KW,6471)
      DO 430 K=1,6
      WRITE(KW,6472) K,(SRCQEX(K,N),N=1,15)
  430 CONTINUE
      WRITE(KW,6473)
      DO 431 K=1,6
      WRITE(KW,6472) K,(SRCQSC(K,N),N=1,15)
  431 CONTINUE
      WRITE(KW,6474)
      DO 432 K=1,6
      WRITE(KW,6472) K,(SRCQCB(K,N),N=1,15)
  432 CONTINUE
      WRITE(KW,6475) TRABCD(3),TRAXSG(KWTRAB+1)
      DO 433 K=1,33
      IF(KWTRAB.EQ.0) WRITE(KW,6472) K,(TRCQAB(K,N),N=1,15)
      IF(KWTRAB.EQ.1) WRITE(KW,6472) K,(TRCQEX(K,N),N=1,15)
      IF(KWTRAB.EQ.2) WRITE(KW,6472) K,(TRCQSC(K,N),N=1,15)
      IF(KWTRAB.EQ.3) WRITE(KW,6472) K,(TRCQCB(K,N),N=1,15)
      IF(KWTRAB.EQ.4) THEN
      DO N=1,15
      TRPI0K(N)=TRCQSC(K,N)/(1.D-10+TRCQEX(K,N))
      END DO
      WRITE(KW,6442) K,(TRPI0K(N),N=1,15)
      ENDIF
  433 CONTINUE
      DO 434 I=1,2
      WRITE(KW,6476)
  434 CONTINUE
C
 6401 FORMAT(' (4A) Aerosol Input for Solar Radiation:'
     +      ,' Aerosol Radiative Parameters'
     +      ,T81,'LIST: SRAEXT(L,K),SRASCT(L,K),SRAGCB(L,K),SRAPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6402 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6403 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6404 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6405 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6406 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6407 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6408 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
C
 6411 FORMAT(' (4B) Cloud Input for Solar Radiation:'
     +      ,'   Cloud Radiative Parameters'
     +      ,T81,'LIST: SRCEXT(L,K),SRCSCT(L,K),SRCGCB(L,K),SRCPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6412 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6413 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6414 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6415 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6416 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6417 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6418 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
C
 6420 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)
 6421 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)
 6422 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)
 6423 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)
 6424 FORMAT(/T11,'SUM COLUMN TAU(0.55) =    BkGrnd    ClimAer   D Dust'
     +       ,'    VolAer    TotAer'/T33,5F10.5)
 6425 FORMAT(' ')
C
 6426 FORMAT(' (4C) Cloud+Aerosol  Output from SOLARM/SGPGXG:'
     +      ,'   Cloud+Aerosol Rad Parameters'
     +      ,T81,'LIST: DBLEXT(L,K),DBLSCT(L,K),DBLGCB(L,K),DBLPI0(L,K)'
     +      //T42,'TAU -- EXTINCTION',T99,'TAU -- SCATTERING'
     +      ,/T24,53('-'),4X,53('-'))
 6427 FORMAT(' LN    PLB     HLB     K=',I3,5I9,7X,'K=',I3,5I9)
 6428 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6429 FORMAT(/1X,T7,'COLUMN AMOUNT=',2X,6F9.6,3X,6F9.6)
 6430 FORMAT(6X,'KANORM=',1I1/T48,'COSBAR',T105,'PIZERO'
     +      ,/T24,53('-'),4X,53('-'))
 6431 FORMAT(' LN     PL     DPL     K=',I3,5I9,7X,'K=',I3,5I9)
 6432 FORMAT(1X,I2,2F8.3,3X,6F9.6,3X,6F9.6)
 6433 FORMAT(/1X,T7,'COLUMN   MEAN=',2X,6F9.6,3X,6F9.6)
 6434 FORMAT(/1X,T7,'ALBEDO RSURFB=',2X,6F9.6,3X,6F9.6)
 6435 FORMAT( 1X,T7,'ALBEDO RSURFX=',2X,6F9.6,3X,6F9.6)
 6436 FORMAT(///T44,'AEROSOL COMPOSITION AND TYPE MIX:'
     +      ,T81,'FACTOR',6X,'VALUE',T107,'TAU(0.55)'/)
 6437 FORMAT(T81,'FGOLDH(',I1,') =',1P,E9.2,5X,0P,F7.4)
 6438 FORMAT(/T81,'SUM COLUMN TAU(0.55) =',F10.4)
 6439 FORMAT(' ')
C
 6440 FORMAT(' (4D) Background Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRAQEX(L,K),SRAQST(L,K),SRAQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,11I8/
     +        '   AEROSOL  ACID1   SSALT   SLFT1   SLFT2   BSLT1'
     +        ,'   BSLT2   DUST1   DUST2   DUST3   CARB1   CARB2'/
     +        '   SIZE      0.5     2.0     0.3     1.0     0.5 '
     +        ,'    2.0     0.5     2.0     8.0     0.1     0.5 ')
 6441 FORMAT('  K  SRAQEX')
 6442 FORMAT(I3,6X,15F8.5)
 6443 FORMAT('  K  SRAQSC')
 6444 FORMAT('  K  SRAQCB')
 6445 FORMAT('  K  ',2A3)
 6446 FORMAT(' ')
C
 6450 FORMAT(' (4E) Climatology Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRBQEX(L,K),SRBQST(L,K),SRBQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X, 6I8/
     +        '   AEROSOL   SO4     SEA     ANT     OCX     BCI '
     +        ,'    BCB'/ ! OCN     OCB     BCB     SSB         '/
     +        '   SIZE ', 6F8.1)
 6451 FORMAT('  K  SRBQEX - DRY')
 6452 FORMAT(I3,6X,15F8.5)
 6453 FORMAT('  K  SRBQSC - DRY')
 6454 FORMAT('  K  SRBQCB - DRY')
 6455 FORMAT('  K  ',2A3,' - DRY')
 6456 FORMAT(' ')
C
 6460 FORMAT(' (4F) Desert Dust Aerosol Solar and Thermal Mie '
     +      ,'Scattering Parameters:'
     +      ,T81,'List: SRDQEX(L,K),SRDQST(L,K),SRDQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,8I8/
     +        '   AEROSOL  CLAY1   CLAY2   CLAY3   CLAY4   SILT1'
     +        ,'   SILT2   SILT3   SILT4                        '/
     +        '   SIZE ',8F8.1)
 6461 FORMAT('  K  SRDQEX')
 6462 FORMAT(I3,6X,15F8.5)
 6463 FORMAT('  K  SRDQSC')
 6464 FORMAT('  K  SRDQCB')
 6465 FORMAT('  K  ',2A3)
 6466 FORMAT(' ')
C
 6470 FORMAT(' (4G) Cloud Input for Solar, Thermal Radiation:'
     +      ,'  Mie Cloud Radiative Properties'
     +      ,T81,'List: SRCQEX(L,K),SRCQST(L,K),SRCQCB(L,K), TRAB Q S G'
     +      /'      KWTRAB=',I1/7X,15I8/
     +        ' WIM CLOUD  WAT05   WAT10   WAT15   WAT20   WAT25'
     +                ,'   ICE05   ICE15   ICE25   ICE50   ICE75'
     +                ,'   MIC05   MIC15   MIC25   MIC50   MIC75')
 6471 FORMAT('  K  SRCQEX')
 6472 FORMAT(I3,6X,15F8.5)
 6473 FORMAT('  K SRCQSC')
 6474 FORMAT('  K SRCQCB')
 6475 FORMAT('  K ',2A3)
 6476 FORMAT(' ')
      GO TO 9999
C
C-------------
  500 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 501 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  501 CONTINUE
      LK=0
      DO 503 K=1,33
      TAUSUM(K)=0. !!sl TAUSL(K)
      DO 502 L=1,NL
      TRTAUK(L,K)=TRGXLK(L,K)+TRCALK(L,K)+TRAALK(L,K)
      TAUSUM(K)=TAUSUM(K)+TRGXLK(L,K)+TRCALK(L,K)+TRAALK(L,K)
  502 CONTINUE
  503 CONTINUE
      WRITE(KW,6501)
      WRITE(KW,6502) (K,K=1,13)
      DO 504 N=1,NL
      L=NLP-N
      WRITE(KW,6503) L,PL(L),TLM(L),(TRTAUK(L,K),K=1,13)
  504 CONTINUE
!sl   WRITE(KW,6504) (TAUSL(K),K=1,13)
      WRITE(KW,6505) (TAUSUM(K),K=1,13)
      WRITE(KW,6506) SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6507) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL.GT.13)  NPAGE=1
      WRITE(KW,6508) NPAGE
      WRITE(KW,6509) (K,K=14,33)
      DO 505 N=1,NL
      L=NLP-N
      WRITE(KW,6510) L,(TRTAUK(L,K),K=14,33)
  505 CONTINUE
!sl   WRITE(KW,6511) ( TAUSL(K),K=14,33)
      WRITE(KW,6512) (TAUSUM(K),K=14,33)
      WRITE(KW,6513) (BGFLUX(K),K=14,33)
      WRITE(KW,6514) (BGFRAC(K),K=14,33)
      DO 506 I=1,10
      WRITE(KW,6515)
  506 CONTINUE
C
 6501 FORMAT(' (5) TAU TABLE FOR THERMAL RADIATION: CONTAINS'
     +      ,' TOTAL SPECIFIED GAS, CLOUD & AEROSOL ABSORPTION'
     +      ,T99,'TRGXLK(L,K),TRCALK(L,K),TRCAAK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6502 FORMAT(' LN     PL     TLM      K=',I1,4X,'K=',I2,9I9,3I8)
 6503 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.4,3F8.3)
!sl6504 FORMAT(/4X,'SURFACE LAYER= ',10F9.4,3F8.3)
 6505 FORMAT( 4X,'COLUMN AMOUNT= ',10F9.4,3F8.3)
 6506 FORMAT(/1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6507 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6508 FORMAT(1I1/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6509 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6510 FORMAT( 1X,  I2,6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
!sl6511 FORMAT(/1X,'SL',6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
 6512 FORMAT( 1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2)
 6513 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6514 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6515 FORMAT(' ')
      GO TO 9999
C
C-------------
  600 CONTINUE
C-------------
C
      NPAGE=1
      IF(INDEX.LT.11) NPAGE=KPAGE
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 601 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  601 CONTINUE
      WRITE(KW,6601)
      WRITE(KW,6602) (K,K=1,13)
      DO 602 N=1,NL
      L=NLP-N
      WRITE(KW,6603) L,PL(L),TLM(L),(TRGXLK(L,K),K=1,13)
  602 CONTINUE
      LK=0
      DO 604 K=1,33
      TAUSUM(K)=0. !!sl  TAUSL(K)
      DO 603 L=1,NL
  603 TAUSUM(K)=TAUSUM(K)+TRGXLK(L,K)
  604 CONTINUE
!sl   WRITE(KW,6604) (TAUSL(K),K=1,13)
      WRITE(KW,6605) (TAUSUM(K),K=1,13)
      WRITE(KW,6606) SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6607) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL.GT.13)  NPAGE=1
      WRITE(KW,6608)
      WRITE(KW,6609) (K,K=14,33)
      DO 605 N=1,NL
      L=NLP-N
      WRITE(KW,6610) L,(TRGXLK(L,K),K=14,33)
  605 CONTINUE
!sl   WRITE(KW,6611) ( TAUSL(K),K=14,33)
      WRITE(KW,6612) (TAUSUM(K),K=14,33)
      WRITE(KW,6613) (BGFLUX(K),K=14,33)
      WRITE(KW,6614) (BGFRAC(K),K=14,33)
      DO 606 I=1,10
      WRITE(KW,6615)
  606 CONTINUE
C
 6601 FORMAT(' (6) TAU TABLE FOR THERMAL RADIATION: INCLUDES ANY'
     +      ,' SPECIFIED OVERLAP, CLOUD & AEROSOL ABSORPTION'
     +      ,T114,'TRGXLK(L,K),TAUSL(L)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6602 FORMAT(' LN     PL     TLM      K=',I1,4X,'K=',I2,9I9,3I8)
 6603 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.4,3F8.3)
 6604 FORMAT(/4X,'SURFACE LAYER= ',10F9.4,3F8.3)
 6605 FORMAT( 4X,'COLUMN AMOUNT= ',10F9.4,3F8.3)
 6606 FORMAT(/1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6607 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6608 FORMAT(/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6609 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6610 FORMAT( 1X,  I2,6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
!sl6611 FORMAT(/1X,'SL',6F7.4,2F6.3,3F6.2,1F6.1,4F7.4,3F6.3,F6.2)
 6612 FORMAT( 1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2)
 6613 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6614 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6615 FORMAT(' ')
      GO TO 9999
C
C-------------
  700 CONTINUE
C-------------
C
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 701 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  701 CONTINUE
      WRITE(KW,6701)
      WRITE(KW,6702) (K,K=1,13)
      DO 702 N=1,NL
      L=NLP-N
      WRITE(KW,6703) L,PL(L),TLM(L),(TRCALK(L,K),K=1,13)
  702 CONTINUE
      LK=0
      DO 704 K=1,33
      TAUSUM(K)=0.0
      DO 703 L=1,NL
      LK=LK+1
  703 TAUSUM(K)=TAUSUM(K)+TRCALK(L,K)
  704 CONTINUE
      WRITE(KW,6704) (TAUSUM(K),K=1,13),(TRCTCA(K),K=1,13)
      WRITE(KW,6705)
      WRITE(KW,6706)  SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6707) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
C
      WRITE(KW,6708)
      WRITE(KW,6709) (K,K=14,33)
      DO 705 N=1,NL
      L=NLP-N
      WRITE(KW,6710) L,(TRCALK(L,K),K=14,33)
  705 CONTINUE
      WRITE(KW,6711) (TAUSUM(K),K=14,33),(TRCTCA(K),K=14,33)
      WRITE(KW,6712) (BGFLUX(K),K=14,33)
      WRITE(KW,6713) (BGFRAC(K),K=14,33)
      DO 706 I=1,8
      WRITE(KW,6714)
  706 CONTINUE
C
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      SUMK=0.0
      DO 711 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      SUMK=SUMK+BGFLUX(K)
      ITG=ITG+ITNEXT
  711 CONTINUE
      WRITE(KW,6721)
      WRITE(KW,6722) (K,K=1,13)
      DO 712 N=1,NL
      L=NLP-N
      WRITE(KW,6723) L,PL(L),TLM(L),(TRAALK(L,K),K=1,13)
  712 CONTINUE
      DO 714 K=1,33
      TAUSUM(K)=0.0
      DO 713 L=1,NL
  713 TAUSUM(K)=TAUSUM(K)+TRAALK(L,K)
  714 CONTINUE
      WRITE(KW,6724) (TAUSUM(K),K=1,13)
      WRITE(KW,6725)
      WRITE(KW,6726)         SUMK,(BGFLUX(K),K=1,13)
      WRITE(KW,6727) TGMEAN,SIGT4,(BGFRAC(K),K=1,13)
      NPAGE=0
      IF(NL.GT.13)  NPAGE=1
      WRITE(KW,6728) NPAGE
      WRITE(KW,6729) (K,K=14,33)
      DO 715 N=1,NL
      L=NLP-N
      WRITE(KW,6730) L,(TRAALK(L,K),K=14,33)
  715 CONTINUE
      WRITE(KW,6731) (TAUSUM(K),K=14,33)
      WRITE(KW,6732) (BGFLUX(K),K=14,33)
      WRITE(KW,6733) (BGFRAC(K),K=14,33)
      DO 716 I=1,12
      WRITE(KW,6734)
  716 CONTINUE
C
 6701 FORMAT(' (7A) TRCALK TABLE FOR THERMAL RADIATION: CONTAINS'
     +      ,' 33 KD CLOUD ABSORPTION OPTICAL DEPTHS AT'
     +      ,' THERMAL WAVELENGTHS ',T117,'LIST: TRCALK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6702 FORMAT(' LN     PL     TLM      K=',I1,6X,I2,9I9,3I8)
 6703 FORMAT(1X,I2,F8.3,F7.2,1X,9F9.5,4F8.5)
 6704 FORMAT(/4X,'COLUMN AMOUNT= ',9F9.4,4F8.5
     +       /4X,'TOPCLD ALBEDO= ',9F9.4,4F8.5)
 6705 FORMAT(/' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:')
 6706 FORMAT( 1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6707 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6708 FORMAT(/T25,'CARBON DIOXIDE:   PRINCIPAL ABSORBER REGION'
     +      ,T93,'OZONE:   PRINCIPAL ABSORBER REGION'
     +      /4X,77('-'),1X,51('-'))
 6709 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6710 FORMAT( 1X,  I2,6F7.4,5F6.4,F6.4,4F7.4,3F6.4,F6.4)
 6711 FORMAT(/1X,'CA',5F7.4,1F7.3,3F6.2,2F6.1,1F6.0,4F7.4,2F6.3,2F6.2
     +       /1X,'TA',5F7.4,1F7.4,3F6.3,2F6.3,1F6.3,4F7.4,2F6.3,2F6.3)
 6712 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6713 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6714 FORMAT(' ')
C
 6721 FORMAT(' (7B) AEROSOL TAU TABLE FOR THERMAL RADIATION:'
     +      ,'  AEROSOL ABSORPTION OPTICAL DEPTH AT THERMAL WAVELENGTHS'
     +      ,T116,'LIST:  TRAALK(L,K)'/
     +      ,/1X,'K-DIST BREAKDOWN:',T23,'WINDOW'
     +      ,3X,'WATER VAPOR:',T71,'PRINCIPAL ABSORBER REGION'
     +      ,/T23,6('-'),3X,101('-'))
 6722 FORMAT(' LN     PL     TLM      K=',I1,6X,I2,9I9,3I8)
 6723 FORMAT(1X,I2,F8.3,F7.2,1X,10F9.5,3F8.5)
 6724 FORMAT(/4X,'COLUMN AMOUNT= ',10F9.5,3F8.5)
 6725 FORMAT(' K-INTERVAL CONTRIBUTIONS:'/' COMPARE WITH GROUND FLUX:')
 6726 FORMAT( 1X,'PF W/M**2= '  ,F6.2,1X,10F9.3,3F8.3)
 6727 FORMAT( 1X,'TG=',F6.2,'= ',F6.2,1X,10F9.4,3F8.3)
 6728 FORMAT(1I1/4X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +      ,T83,'OZONE:',T100,'PRINCIPAL ABSORBER REGION'
     +      /4X,76('-'),2X,50('-'))
 6729 FORMAT(1X,'LN  K=',I2,5I7,6I6,3X,'K=',I2,3I7,6I6)
 6730 FORMAT( 1X,  I2,6F7.5,2F6.4,3F6.4,F6.4,4F7.4,3F6.4,F6.4)
 6731 FORMAT( 1X,'CA',5F7.5,1F7.5,3F6.4,2F6.4,1F6.4,4F7.4,2F6.4,2F6.4)
 6732 FORMAT(/1X,'PF',1F7.4,5F7.3,1F6.2,3F6.3,2F6.3,2F7.3,2F7.4,4F6.3)
 6733 FORMAT( 1X,'FR',6F7.4,2F6.3,3F6.3,1F6.3,4F7.4,3F6.3,F6.3)
 6734 FORMAT(' ')
      GO TO 9999
C
C-------------
  800 CONTINUE
C-------------
C
      WRITE(KW,6800)
      DO 801 K=1,16
      ISR1(K)=NORDER(K)
      IF(KORDER.EQ.1) ISR1(K)=K
  801 CONTINUE
      WRITE(KW,6801) (ISR1(K),K=1,16)
      SUMK=0.0
      DO 802 K=1,16
      FSR1(K)=DKS0(NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=DKS0(K)
      SUMK=SUMK+FSR1(K)
  802 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6802) (FSR1(K),K=1,17)
      DO 803 K=1,16
      ISR1(K)=NMWAVA(K)
      IF(KORDER.EQ.1) ISR1(K)=NMWAVA(IORDER(K))
  803 CONTINUE
      WRITE(KW,6803) (ISR1(K),K=1,16)
      DO 804 K=1,16
      ISR1(K)=NMWAVB(K)
      IF(KORDER.EQ.1) ISR1(K)=NMWAVB(IORDER(K))
  804 CONTINUE
      WRITE(KW,6804) (ISR1(K),K=1,16)
      IF(KORDER.EQ.0) WRITE(KW,6805)
      IF(KORDER.EQ.1) WRITE(KW,6806)
      DO 805 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  805 CONTINUE
      WRITE(KW,6807) (ISR1(K),K=1,16)
      DO 807 N=1,NLP
      L=NLP+1-N
      SUMK=0.0
      DO 806 K=1,16
      FSR1(K)=SKDFLB(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKDFLB(L,K)
      SUMK=SUMK+FSR1(K)
  806 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6808) L,(FSR1(K),K=1,17)
  807 CONTINUE
      DO 808 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1) ISR1(K)=IORDER(K)
  808 CONTINUE
      WRITE(KW,6809) (ISR1(K),K=1,16)
      DO 810 N=1,NLP
      L=NLP+1-N
      SUMK=0.0
      DO 809 K=1,16
      FSR1(K)=SKUFLB(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKUFLB(L,K)
      SUMK=SUMK+FSR1(K)
  809 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6810) L,(FSR1(K),K=1,17)
  810 CONTINUE
      DO 811 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1) ISR1(K)=IORDER(K)
  811 CONTINUE
      WRITE(KW,6811) (ISR1(K),K=1,16)
      SUMT=0.D0
      SUMK=0.D0
      DO 812 K=1,16
      FSR1(K)=SRKALB(NORDER(K))
      FSR2(K)=DKS0(NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SRKALB(K)
      IF(KORDER.EQ.1) FSR2(K)=DKS0(K)
      SUMK=SUMK+FSR1(K)*FSR2(K)
  812 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6812) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      SUMK1=0.D0
      SUMK2=0.D0
      DO 813 K=1,16
      FSR1(K)=SKNFLB(NLP,NORDER(K))-SKNFLB(1,NORDER(K))
      FSR2(K)=SKDFLB(NLP,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKNFLB(NLP,K)-SKNFLB(1,K)
      IF(KORDER.EQ.1) FSR2(K)=SKDFLB(NLP,K)
      SUMK1=SUMK1+FSR1(K)
      SUMK2=SUMK2+FSR2(K)
      FSR1(K)=FSR1(K)/(FSR2(K)+1.d-20)
  813 CONTINUE
      FSR1(17)=SUMK1/(SUMK2+1.d-20)
      WRITE(KW,6813) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      SUMK1=0.D0
      SUMK2=0.D0
      DO 814 K=1,16
      FSR1(K)=SKNFLB(1,NORDER(K))
      FSR2(K)=SKDFLB(NLP,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKNFLB(1,K)
      IF(KORDER.EQ.1) FSR2(K)=SKDFLB(NLP,K)
      SUMK1=SUMK1+FSR1(K)
      SUMK2=SUMK2+FSR2(K)
      FSR1(K)=FSR1(K)/(FSR2(K)+1.d-20)
  814 CONTINUE
      FSR1(17)=SUMK1/(SUMK2+1.d-20)
      WRITE(KW,6814) (FSR1(K),K=1,17)
      SUMT=SUMT+FSR1(17)
      DO 815 K=1,16
      ISR1(K)=KSLAMW(NORDER(K))
      IF(KORDER.EQ.1) ISR1(K)=KSLAMW(K)
  815 CONTINUE
      WRITE(KW,6815) SUMT,(ISR1(K),K=1,16)
      SUMK=0.D0
      DO 816 K=1,16
      KK=KSLAMW(NORDER(K))
      IF(KORDER.EQ.1) KK=KSLAMW(K)
      FSR1(K)=SRBALB(KK)
      FSR2(K)=SRXALB(KK)
  816 CONTINUE
      WRITE(KW,6816) (FSR1(K),K=1,16)
      WRITE(KW,6817) (FSR2(K),K=1,16)
      WRITE(KW,6818)   COSZ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR
      WRITE(KW,6819) SRXVIS,SRXNIR,SRDVIS,SRUVIS,ALBVIS,SRDNIR
     +                            ,SRUNIR,ALBNIR
      WRITE(KW,6820) SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR
      DO 817 I=1,1
      IF(KORDER.EQ.1) WRITE(KW,6821)
  817 CONTINUE
C
      WRITE(KW,6840)
      DO 821 K=1,16
      ISR1(K)=NORDER(K)
      IF(KORDER.EQ.1) ISR1(K)=K
  821 CONTINUE
      WRITE(KW,6841) (ISR1(K),K=1,16)
      SUMK=0.0
      DO 822 K=1,16
      FSR1(K)=DKS0(NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=DKS0(K)
      SUMK=SUMK+FSR1(K)
  822 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6842) (FSR1(K),K=1,17)
      IF(KORDER.EQ.0) WRITE(KW,6843)
      IF(KORDER.EQ.1) WRITE(KW,6844)
      DO 825 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  825 CONTINUE
      WRITE(KW,6845) (ISR1(K),K=1,16)
      DO 827 N=1,NLP
      L=NLP+1-N
      SUMK=0.0
      DO 826 K=1,16
      FSR1(K)=SKNFLB(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKNFLB(L,K)
      SUMK=SUMK+FSR1(K)
  826 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6846) L,(FSR1(K),K=1,17)
  827 CONTINUE
      DO 828 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1) ISR1(K)=IORDER(K)
  828 CONTINUE
      WRITE(KW,6847) (ISR1(K),K=1,16)
      DO 830 N=1,NL
      L=NL+1-N
      SUMK=0.0
      DO 829 K=1,16
      FSR1(K)=SKFHRL(L,NORDER(K))
      IF(KORDER.EQ.1) FSR1(K)=SKFHRL(L,K)
      SUMK=SUMK+FSR1(K)
  829 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6848) L,(FSR1(K),K=1,17)
  830 CONTINUE
      DO 831 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  831 CONTINUE
      WRITE(KW,6849) (ISR1(K),K=1,16)
      DO 833 N=1,4
      SUMK=0.0
      DO 832 K=1,16
      FSR1(K)=SRKGAX(NORDER(K),N)
      IF(KORDER.EQ.1) FSR1(K)=SRKGAX(K,N)
      SUMK=SUMK+FSR1(K)
  832 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6850) L,(FSR1(K),K=1,17)
  833 CONTINUE
      DO 834 K=1,16
      ISR1(K)=K
      IF(KORDER.EQ.1)ISR1(K)=IORDER(K)
  834 CONTINUE
      WRITE(KW,6851)
      DO 836 N=1,4
      SUMK=0.0
      DO 835 K=1,16
      FSR1(K)=SRKGAD(NORDER(K),N)
      IF(KORDER.EQ.1) FSR1(K)=SRKGAD(K,N)
      SUMK=SUMK+FSR1(K)
  835 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6852) N,(FSR1(K),K=1,17)
  836 CONTINUE
      WRITE(KW,6853)
      DO 838 N=1,4
      SUMK=0.0
      DO 837 K=1,16
      FSR1(K)=SRKGAX(NORDER(K),N)
      FSR2(K)=SRKGAD(NORDER(K),N)
      IF(KORDER.EQ.1) FSR1(K)=SRKGAX(K,N)
      IF(KORDER.EQ.1) FSR2(K)=SRKGAD(K,N)
      FSR1(K)=FSR1(K)+FSR2(K)
      SUMK=SUMK+FSR1(K)
  837 CONTINUE
      FSR1(17)=SUMK
      WRITE(KW,6854) N,(FSR1(K),K=1,17)
  838 CONTINUE
      WRITE(KW,6855) SRNFLB(1),POCEAN,FSRNFG(1),PEARTH,FSRNFG(2)
     +                        ,POICE ,FSRNFG(3),PLICE ,FSRNFG(4)
C
 6800 FORMAT(' (8A)     SPECTRAL/k-DISTRIBUTION COMPONENT BREAKDOWN'
     +      ,' FOR DOWNWARD AND UPWARD SOLAR RADIATIVE FLUXES'
     +      ,T108,'SKDFLB(L,K)  SKUFLB(L,K)  SRKALB(K)'/)
 6801 FORMAT('     K=',I5,2I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6802 FORMAT('  DKS0=',F6.3,2F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6803 FORMAT(' NMWAVA=',I6,2I8,2I7,I8,I9,7I8,I7,I8)
 6804 FORMAT(' NMWAVB=',I6,2I8,2I7,I8,I9,7I8,I7,I8)
 6805 FORMAT(' ABSORB'/'   GAS= O3,O2  O3,NO2      O2     O2     O2'
     +        ,'     H2O',22X,'H2O     H2O     H2O     H2O     CO2'
     +        ,'     CO2    CO2  CO2,H2O,O2'
     +        /' SKDFLB (Downward Spectral Flux)'
     +        /6X,6('-'),'VIS',6('-'),2X,46('-'),'NIR',59('-'))
 6806 FORMAT(' ABSORB'/'   GAS=   H2O     H2O     H2O    H2O    H2O'
     +        ,'      O2       O2      O2     CO2     CO2     CO2',18X
     +        ,'O3,NO2  O3,O2 CO2,H2O,O2'/'SKDFLB  (Downard Spectral'
     +        ,' Flux)',T110,6('-'),'VIS',5('-'))
 6807 FORMAT('  N  L=',I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6808 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6809 FORMAT(/' SKUFLB (Upward Spectral Flux)'/'  N  L='
     +         ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6810 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6811 FORMAT(/' SRKALB '
     +         ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6812 FORMAT('   TOA='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6813 FORMAT(' ABSORB'/'  ATMO='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6814 FORMAT(' ABSORB'/'  SURF='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6815 FORMAT(' ALSURF',T133,'Sum=',F6.4
     +       /' KSLAM= ',I3,I9,I8,2I7,I8,I9,7I8,I7,I8)
 6816 FORMAT('   SRX='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6817 FORMAT('   SRB='
     +       ,F5.4,F9.4,F8.4,2F7.4,F8.4,F9.4,7F8.4,F7.4,F8.4,F11.4)
 6818 FORMAT(/' At Top of Atm:  ',' COSZ  =',F6.4,14X
     +      , 2X,' SRIVIS=',F7.3,'  SROVIS=',F7.3,   '   PLAVIS=',F6.4
     +      , 2X,' SRINIR=',F7.3,'  SRONIR=',F7.3,   '   PLANIR=',F6.4)
 6819 FORMAT( ' At Bot of Atm:  ',' SRXVIS=',F6.4,1X,' SRXNIR=',F6.4
     +      , 1X,' SRDVIS=',F7.3,'  SRUVIS=',F7.3,   '   ALBVIS=',F6.4
     +      , 2X,' SRDNIR=',F7.3,'  SRUNIR=',F7.3,   '   ALBNIR=',F6.4)
 6820 FORMAT( ' In Atmosphere:  ',' (VIS=0.53*S0)',2X,'(NIR=0.47*S0)'
     +      , 1X,' SRTVIS=',F7.5,'  SRRVIS=',F7.5,   '   SRAVIS=',F6.4
     +      , 2X,' SRTNIR=',F7.5,'  SRRNIR=',F7.5,   '   SRANIR=',F6.4)
 6821 FORMAT(' ')

 6840 FORMAT(' (8B)  SPECTRAL/k-DISTRIBUTION COMPONENT BREAKDOWN'
     +      ,' FOR NET DOWNWARD SOLAR FLUX & HEATING RATE'
     +      ,T106,'SKNFLB(L,K)  SKFHRL(L,K)  SRKGAX(L,I)'/)
 6841 FORMAT('     K=',I5,2I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6842 FORMAT('  DKS0=',F6.3,2F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6843 FORMAT('   GAS= O3,O2  O3,NO2      O2     O2     O2'
     +        ,'     H2O',22X,'H2O     H2O     H2O     H2O     CO2',5X
     +        ,'CO2    CO2  CO3,H2O,O2'/' SKNFLB (Spectral Net Flux)')
 6844 FORMAT('   GAS=   H2O     H2O     H2O    H2O    H2O'
     +        ,'      O2       O2      O2     CO2     CO2     CO2',18X
     +        ,'O3,NO2  O3,O2 CO2,H2O,O2'
     +        /' SKDFLB (Spectral Net Flux)',T110,6('-'),'VIS',5('-'))
 6845 FORMAT('  N  L=',I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6846 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6847 FORMAT(/' SKFHRL (Spectral Heating Rate)'/'  N  L='
     +       ,I4,I9,I8,2I7,I8,I9,7I8,I7,I8,'       Total')
 6848 FORMAT(I3,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6849 FORMAT(/' SRKGAX (Direct Beam Spectral Absorption at Ground)'
     +       /' N   L=',I4,8I8,I9,4I8,I7,I8,'      Total')
 6850 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6851 FORMAT( ' SRKGAD (Diffuse Spectral Absorption at Ground)')
 6852 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6853 FORMAT( ' SRKGAD (Total Spectral Absorption at Ground)')
 6854 FORMAT(I2,1X,2F9.3,F8.3,2F7.3,F8.3,F9.3,7F8.3,F7.3,F8.3,F11.3)
 6855 FORMAT(/' Absorption at Ground by Surface-type'
     +       ,T39,'SRNFLB(1) = POCEAN * FSRNFG(1) + PEARTH * FSRNFG(2) '
     +                     ,'+  POICE * FSRNFG(3) +  PLICE * FSRNFG(4) '
     +       /T39,F7.3,'   = ',F6.4,' *',F8.3,'   + ',F6.4,' *',F8.3
     +                ,'   + ',F6.4,' *',F8.3,'   + ',F6.4,' *',F8.3)
      GO TO 9999
C-------------
  900 CONTINUE
C-------------
C
c      SIGMA=5.6697D-08
      TGMEAN=POCEAN*TGO**4+PEARTH*TGE**4+PLICE*TGLI**4+POICE*TGOI**4
      TGMEAN=SQRT(TGMEAN)
      TGMEAN=SQRT(TGMEAN)
      SIGT4=SIGMA*TGMEAN**4
      ITG=TGMEAN
      WTG=TGMEAN-ITG
      ITG=ITG-ITPFT0
      DO 901 K=1,33
      BGFLUX(K)=PLANCK(ITG)-(PLANCK(ITG)-PLANCK(ITG+1))*WTG
      BGFRAC(K)=BGFLUX(K)/SIGT4
      ITG=ITG+ITNEXT
  901 CONTINUE
      DO 910 NW=1,5
      DO 903 K=1,33
      DO 902 L=1,NLP
      IF(NW.EQ.1) WFLB(L,K)=DFLB(L,K)
      IF(NW.EQ.2) WFLB(L,K)=UFLB(L,K)
      IF(NW.EQ.3) WFLB(L,K)=UFLB(L,K)-DFLB(L,K)
      IF(NW.GT.3.AND.L.GT.NL) GO TO 902
      IF(NW.EQ.4) WFLB(L,K)=WFLB(L+1,K)-WFLB(L,K)
      IF(NW.EQ.5.AND.ABS(TRFCRL(L)).LT.1.E-10) WFLB(L,K)=1.E-30
      IF(NW.EQ.5) WFLB(L,K)=WFLB(L,K)/(ABS(TRFCRL(L))+1.E-10)
  902 CONTINUE
      IF(NW.EQ.1) WFSL(K)=DFSL(K)
      IF(NW.EQ.2) WFSL(K)=UFSL(K)
      IF(NW.EQ.3) WFSL(K)=UFSL(K)-DFSL(K)
      IF(NW.EQ.4) WFSL(K)=WFSL(K)-UFLB(1,K)+DFLB(1,K)
!sl   IF(NW.EQ.5.AND.ABS(TRSLCR).LT.1.E-10) WFSL(K)=1.E-30
      IF(NW.EQ.5) WFSL(K)=0.   !nu =WFSL(K)/(ABS(TRSLCR)+1.E-10)
  903 CONTINUE
      DO 907 L=1,NLP
      IF(L.GT.NL.AND.NW.GT.3) GO TO 907
      ASUM1=0.
      BSUM1=0.
      CSUM1=0.
      DSUM1=0.
      ESUM1=0.
      FSUM1=0.
      SUM=0.
      DO 904 K=2,13
      ASUM1=ASUM1+  WFSL(K)
      BSUM1=BSUM1+BGFEMT(K)
      CSUM1=CSUM1+BGFLUX(K)
      DSUM1=DSUM1+BGFRAC(K)
      ESUM1=ESUM1+TRCTCA(K)
      FSUM1=FSUM1+TRGALB(K)
  904 SUM=SUM+WFLB(L,K)
      SUM1(L)=SUM
      ASUM2=0.
      BSUM2=0.
      CSUM2=0.
      DSUM2=0.
      ESUM2=0.
      FSUM2=0.
      SUM=0.
      DO 905 K=14,25
      ASUM2=ASUM2+  WFSL(K)
      BSUM2=BSUM2+BGFEMT(K)
      CSUM2=CSUM2+BGFLUX(K)
      DSUM2=DSUM2+BGFRAC(K)
      ESUM2=ESUM2+TRCTCA(K)
      FSUM2=FSUM2+TRGALB(K)
  905 SUM=SUM+WFLB(L,K)
      SUM2(L)=SUM
      ASUM3=0.
      BSUM3=0.
      CSUM3=0.
      DSUM3=0.
      ESUM3=0.
      FSUM3=0.
      SUM=0.
      DO 906 K=26,33
      ASUM3=ASUM3+  WFSL(K)
      BSUM3=BSUM3+BGFEMT(K)
      CSUM3=CSUM3+BGFLUX(K)
      DSUM3=DSUM3+BGFRAC(K)
      ESUM3=ESUM3+TRCTCA(K)
      FSUM3=FSUM3+TRGALB(K)
  906 SUM=SUM+WFLB(L,K)
      SUM3(L)=SUM
  907 CONTINUE
C
      NPAGE=1
      WRITE(KW,6901) NW,FTYPE(NW)
      WRITE(KW,6902) (K,K=1,13)
      DO 908 N=1,NLP
      L=NLP+1-N
      IF(L.GT.NL.AND.NW.GT.3) GO TO 908
      SUML=SUM1(L)+SUM2(L)+SUM3(L)+WFLB(L,1)
      WRITE(KW,6903) L,PL(L),SUML,SUM1(L),SUM2(L),SUM3(L)
     +               ,(WFLB(L,K),K=1,13)
  908 CONTINUE
      SUMA=ASUM1+ASUM2+ASUM3+  WFSL(1)
      SUMB=BSUM1+BSUM2+BSUM3+BGFEMT(1)
      SUMC=CSUM1+CSUM2+CSUM3+BGFLUX(1)
      SUMD=DSUM1+DSUM2+DSUM3+BGFRAC(1)
      SUME=ESUM1+ESUM2+ESUM3+TRCTCA(1)
      SUMF=FSUM1+FSUM2+FSUM3+TRGALB(1)
      WRITE(KW,6904) SUMA,ASUM1,ASUM2,ASUM3,(  WFSL(K),K=1,13)
      WRITE(KW,6905) SUMB,BSUM1,BSUM2,BSUM3,(BGFEMT(K),K=1,13)
      WRITE(KW,6906) SUMC,CSUM1,CSUM2,CSUM3,(BGFLUX(K),K=1,13)
      WRITE(KW,6907) SUMD,DSUM1,DSUM2,DSUM3,(BGFRAC(K),K=1,13)
      WRITE(KW,6908) SUME,ESUM1,ESUM2,ESUM3,(TRCTCA(K),K=1,13)
      WRITE(KW,6909) SUMF,FSUM1,FSUM2,FSUM3,(TRGALB(K),K=1,13)
      NPAGE=0
      WRITE(KW,6910) NPAGE
      WRITE(KW,6911) (K,K=14,33)
      DO 909 N=1,NLP
      L=NLP+1-N
      IF(L.GT.NL.AND.NW.GT.3) GO TO 909
      WRITE(KW,6912) L,(WFLB(L,K),K=14,33)
  909 CONTINUE
      WRITE(KW,6913) (  WFSL(K),K=14,33)
      WRITE(KW,6914) (BGFEMT(K),K=14,33)
      WRITE(KW,6915) (BGFLUX(K),K=14,33)
      WRITE(KW,6916) (BGFRAC(K),K=14,33)
      WRITE(KW,6917) (TRCTCA(K),K=14,33)
      WRITE(KW,6918) (TRGALB(K),K=14,33)
      LINFIL=2
      IF(NW.GT.3) LINFIL=4
      DO 911 I=1,LINFIL
      WRITE(KW,6919)
  911 CONTINUE
  910 CONTINUE
C
 6901 FORMAT(' (9.',I1,') THERMAL RADIATION: K-DISTRIBUTION'
     +       ,' BREAKDOWN FOR  ',1A8,' FLUX'/
     +       /T21,'PRINCIPAL REGION SUM',2X,'WINDOW'
     +       ,T52,'WATER VAPOR:',T76,'PRINCIPAL ABSORBER REGION'
     +       /20X,20('-'),2X,6('-'),3X,81('-'))
 6902 FORMAT(1X,'LN     PL    TOTAL    H2O   CO2    O3     K='
     +       ,I2,4X,'K=',I2,12I7)
 6903 FORMAT(1X,I2,2F8.2,3F7.2,F8.3,1X,12F7.3)
 6904 FORMAT(/' SL',  9X,4F7.2,F8.3,1X,12F7.3)
 6905 FORMAT(/' BG',  9X,4F7.2,F8.3,1X,12F7.3)
 6906 FORMAT( ' PF',  9X,4F7.2,F8.3,1X,12F7.3)
 6907 FORMAT( ' FR',  9X,4F7.2,F8.3,1X,12F7.3)
 6908 FORMAT(/' AC',  9X,4F7.2,F8.3,1X,12F7.3)
 6909 FORMAT( ' AG',  9X,4F7.2,F8.3,1X,12F7.3)
 6910 FORMAT(1I1/5X,'CARBON DIOXIDE:',T36,'PRINCIPAL ABSORBER REGION'
     +       ,T85,'OZONE:',T101,'PRINCIPAL ABSORBER REGION'
     +       /5X,76('-'),3X,48('-'))
 6911 FORMAT(1X,'LN  K=',I2,6I7,5I6,4X,'K=',I2,1I7,6I6)
 6912 FORMAT( 1X,I2,7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6913 FORMAT(/' SL',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6914 FORMAT(/' BG',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6915 FORMAT( ' PF',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6916 FORMAT( ' FR',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6917 FORMAT(/' AC',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6918 FORMAT( ' AG',7F7.3,5F6.3,1X,2F7.3,6F6.3)
 6919 FORMAT(' ')
      RETURN
C-------------
 1000 CONTINUE
C-------------
C
 9999 CONTINUE
      RETURN
      END SUBROUTINE WRITER

      SUBROUTINE WRITET(KWRU,INDEX,JYRREF,JYRNOW,JMONTH,KLIMIT)
      IMPLICIT NONE
C
C
C     ------------------------------------------------------------------
C     WRITET  GHG, Solar UV, Ozone, Aerosol Trend Diagnostic Information
C
C         INDEX
C           1     GHG DT0 Trends / FULGAS Ratios for CO2,NO2,CH4,F11,F12
C           2     GHG DF  Change / Ann Increase Rate CO2,NO2,CH4,F11,F12
C           3     Lean Solar Constant, UV Spectral Variation Time Trends
C           4     Ozone Zonal-mean (Latitude and Vertical) Distributions
C           5     Ozone Surface-150mb, 150mb-TOA, Column Longitude Distr
C                 A  O3 (Wang-Jacobs) Relative Longitudinal Distribution
C                 B  O3 (London-NCAR) Relative Longitudinal Distribution
C                 C  O3 (W-J, London) Relative Longitudinal Distribution
C           6     Tropospheric Climatology Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           7     Tropospheric Desert Dust Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           8     Stratospheric (Volcanic) Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C           9     Total Column Atmospheric Aerosol Latitude/Height Distr
C                 A  Zonal-mean Extinction Optical Depth
C                 B  Zonal-mean Single Scattering Albedo
C                 C  Zonal-mean Asymmetry Parameter
C        NOTE:
C                 Time Trend (year) Specification is by JYRREF to JYRNOW
C                 Time Specification (O3,Aerosol) is by JYRREF to JMONTH
C                      (If JMONTH = 0, JDAY is used)
C
C                 INDEX < 10 is selective, INDEX > 10 is digit inclusive
C                 KLIMIT = 0 full output,  KLIMIT > 0 abbreviated output
C                 KWRU directs the output to selected (KWRU) file number
C     ------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: KWRU,INDEX,JYRREF,JYRNOW,JMONTH,KLIMIT

      REAL*8 WREF(7),WDAT(7),WPPM(7),DWM2(7),XDT0(5),XRAT(5)
      REAL*8, DIMENSION(49,LX) :: QX,QS,QG,QP,O3
      REAL*8, DIMENSION(49) :: QXCOL,QSCOL,QGCOL,QPCOL,O3COL
      REAL*8 SFL0(5),SFLX(5),DFLX(5),RFLX(5),O3L(46,72)
      INTEGER :: LO3(36)
C
      INTEGER, PARAMETER :: NSW1=24, NSW2=32, NSW3=40, NSW4=48
C
      CHARACTER*32, DIMENSION(4), PARAMETER :: CHAER =(/
     *     'Tropospheric Climatology Aerosol',
     +     'Tropospheric Desert Dust Aerosol',
     +     'Stratospheric (Volcanic) Aerosol',
     +     'Total Column Atmospheric Aerosol'/)

      REAL*8 YREF11,ZREF12,YNOW11,ZNOW12,SDT0,DFL0,SUMO3,QOSH,QONH,QOGL
     *     ,SUMXL,SUMGL,SUMSL,QXSH,QXNH,QXGL,QSSH,QSNH,QSGL,QPSH,QPNH
     *     ,QPGL,QGSH,QGNH,QGGL
      INTEGER KW,INDJ,INDI,INDX,KINDEX,I,JJDAYG,JYEARG,KWSKIP,J,IYEAR
     *     ,NSPACE,LMO,K,mavg,iyr1,lmax,icyc,JYEARS,M,JJDAYO,L,JJ,N,N1
     *     ,N2,II,KAEROS,L1,KA,JJDAY
C
      KW=KWRU
      INDJ=MOD(INDEX,10)
      IF(INDJ.LT.1) INDJ=10
      INDI=1
      IF(INDEX.EQ.0)  INDJ=1
      IF(INDEX.LT.11) INDI=INDJ
      DO 9999 INDX=INDI,INDJ
C
      GO TO (100,100,300,400,500,600,600,600,600,1000),INDX
C
C-------------
  100 CONTINUE
C-------------
C
      KINDEX=INDX
      DO 110 I=1,5
      WREF(I)=XREF(I)
      WPPM(I)=PPMV80(I+4)
  110 CONTINUE
      WREF(6)=PPMV80(11)*1000.D0
      WREF(7)=PPMV80(12)*1000.D0
      WPPM(1)=PPMV80(2)
      WPPM(6)=PPMV80(11)
      WPPM(7)=PPMV80(12)
C
      YREF11=PPMV80(11)
      ZREF12=PPMV80(12)
C
      JJDAYG=184
C
      IF(KINDEX.EQ.1) THEN
      WRITE(KW,6101) JJDAYG
 6101 FORMAT(/1X,'(1)=INDEX'
     +      ,T12,'JDAY=',I3,'   RCM RAD EQUIL NO-FEEDBACK DT0'
     +      ,T55,'PRESENT TREND UPDGHG INPUT DATA TO GCM'
     +      ,T96,'FULGAS FACTOR RELATIVE TO 1980 AMOUNTS')
       WRITE(KW,6102) KTREND
 6102 FORMAT(1X,'KTREND=',I2,1X,40('-'),3X,38('-'),3X,38('-')
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPBF11  PPBF12'
     +         ,         '   FULCO2  FULN2O  FULCH4  FULF11  FULF12')
      ENDIF
C
      IF(KINDEX.EQ.2) THEN
      WRITE(KW,6201) JJDAYG
 6201 FORMAT(/1X,'(2)=INDEX'
     +      ,T12,'JDAY=',I3,'   RCM EQ NO-FEEDBACK DFLUX W/M2'
     +      ,T55,'PRESENT TREND UPDGHG INPUT DATA TO GCM'
     +      ,T96,'ANNUAL CHANGE RATE OF TRACE GAS AMOUNT')
      WRITE(KW,6202) KTREND
 6202 FORMAT(1X,'KTREND=',I2,1X,40('-'),3X,38('-'),3X,38('-')
     +      /1X,'YEAR DTSUM  *DTCO2   DTN2O   DTCH4   DTF11   DTF12'
     +         ,         '   PPMCO2  PPMN20  PPMCH4  PPBF11  PPBF12'
     +         ,         '   RATCO2  RATN2O  RATCH4  RATF11  RATF12')
       ENDIF
C
      JYEARG=JYRREF-1
      CALL UPDGHG(JYEARG,JJDAYG)
      YNOW11=FULGAS(11)*YREF11
      ZNOW12=FULGAS(12)*ZREF12
      XDT0(4)=XDT0(4)+0.066*(YNOW11-YREF11)
      XDT0(5)=XDT0(5)+0.084*(ZNOW12-ZREF12)
C
      DO 120 I=1,5
      WDAT(I)=XNOW(I)
  120 CONTINUE
C
      DO 230 J=JYRREF,JYRNOW
      KWSKIP=0
      IF(J.GT.JYRREF) KWSKIP=KLIMIT
      IF(J.EQ.1980)   KWSKIP=0
      IF(J.EQ.JYRNOW) KWSKIP=0
      JYEARG=J
      CALL UPDGHG(JYEARG,JJDAYG)
      YNOW11=FULGAS(11)*YREF11
      ZNOW12=FULGAS(12)*ZREF12
      XDT0(4)=XDT0(4)+0.066*(YNOW11-YREF11)
      XDT0(5)=XDT0(5)+0.084*(ZNOW12-ZREF12)
      SDT0=0.D0
      DFL0=0.D0
      DO 220 I=1,5
      SDT0=SDT0+XDT0(I)
      XRAT(I)=(XNOW(I)-WDAT(I))/(1.D-10+WDAT(I))
      IF(XRAT(I).GT.9.9999) XRAT(I)=9.9999
      WDAT(I)=XNOW(I)
      DWM2(I)=XDT0(I)*3.300D0
      DFL0=DFL0+DWM2(I)
  220 CONTINUE
      IYEAR=JYEARG
      IF(KINDEX.EQ.1) THEN
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6103) IYEAR,SDT0,(XDT0(I),I=1,5),(XNOW(I),I=1,5)
     +              ,FULGAS(2),(FULGAS(I),I=6,9)
 6103 FORMAT(1X,I4,F6.3,5F8.4,1X,F8.2,4F8.4,1X,5F8.4)
      ENDIF
      IF(KINDEX.EQ.2) THEN
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6203) IYEAR,DFL0,(DWM2(I),I=1,5),(XNOW(I),I=1,5)
     +              ,(XRAT(I),I=1,5)
 6203 FORMAT(1X,I4,F6.3,5F8.4,1X,F8.2,4F8.4,1X,5F8.4)
      ENDIF
      NSPACE=IYEAR-(IYEAR/10)*10
      IF(KLIMIT.GT.0) GO TO 230
      IF(NSPACE.EQ.0) WRITE(KW,6104)
 6104 FORMAT(' ')
  230 CONTINUE
      GO TO 9999
C
C-------------
  300 CONTINUE
C-------------
C
      if(ksolar.lt.0) go to 9999
      LMO=(1950-iy1S0)*12+1
      if(ksolar.gt.1) LMO=nint(1950 - yr1s0 + 1.5)
      DO 310 I=1,5
      SFL0(I)=0.D0
  310 CONTINUE
      DO 320 K=1,190
      IF(K.LE.NSW1)              SFL0(1)=SFL0(1)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW1.AND.K.LE.NSW2)SFL0(2)=SFL0(2)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW2.AND.K.LE.NSW3)SFL0(3)=SFL0(3)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW3.AND.K.LE.NSW4)SFL0(4)=SFL0(4)+UVLEAN(LMO,K)*DSLEAN(K)
                                 SFL0(5)=SFL0(5)+UVLEAN(LMO,K)*DSLEAN(K)
  320 CONTINUE
C
      if(ksolar.eq.2)
     *   WRITE(KW,6299) int(yr1s0),int(yr2s0),JYRREF,JYRNOW,SFL0(5)
      if(ksolar.lt.2) WRITE(KW,6300) JYRREF,JYRNOW,SFL0(5)
 6299 FORMAT(/' (3)=INDEX  Annual-mean Solar flux (from J.Lean annual'
     +      ,I6,'-',I4,' data) for JYRREF=',I4,' to JYRNOW=',I4,'  mid'
     +      ,' 1950 Ref S00WM2=',F9.4/12X,'Solar UV Spectral Flux W/m2'
     +      ,T57,'Delta Solar UV Spectral Flux W/m2'
     +      ,T97,'Solar UV Spectral Flux Ratios'
     +      /'  YEAR    0-280 280-320 320-360 360-400   Total '
     +      ,6X,'0-280 280-320 320-360 360-400   Total '
     +      ,4X,'0-280 280-320 320-360 360-400   Total ')
 6300 FORMAT(/' (3)=INDEX  Annual-mean Solar flux (from J.Lean monthly'
     +      ,' 1882-1998 data) for JYRREF=',I4,' to JYRNOW=',I4,'  Jan'
     +      ,' 1950 Ref S00WM2=',F9.4/12X,'Solar UV Spectral Flux W/m2'
     +      ,T57,'Delta Solar UV Spectral Flux W/m2'
     +      ,T97,'Solar UV Spectral Flux Ratios'
     +      /'  YEAR    0-280 280-320 320-360 360-400   Total '
     +      ,6X,'0-280 280-320 320-360 360-400   Total '
     +      ,4X,'0-280 280-320 320-360 360-400   Total ')
C
      if(ksolar.lt.2) then
        mavg = 12
        iyr1 = iy1s0
        lmax = ms0x
        icyc = mcycs0
      else
        mavg = 1
        iyr1 = yr1s0
        lmax = nint( yr2s0-yr1s0+1 )
        icyc = icycs0
      end if
      DO 370 J=JYRREF,JYRNOW
      KWSKIP=0
      IF(J.GT.JYRREF) KWSKIP=KLIMIT
      IF(J.EQ.JYRNOW) KWSKIP=0
      JYEARS=J
      DO 330 I=1,5
      SFLX(I)=0.D0
  330 CONTINUE
      LMO=(JYEARS-iyr1)*mavg
      DO 350 M=1,mavg
      LMO=LMO+1
      IF(LMO.GT.lmax) LMO=LMO-icyc*((LMO-lmax+icyc-1)/icyc)
      IF(LMO.LT.1) LMO=LMO+icyc*((icyc-LMO)/icyc)
      DO 340 K=1,190
      IF(K.LE.NSW1)              SFLX(1)=SFLX(1)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW1.AND.K.LE.NSW2)SFLX(2)=SFLX(2)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW2.AND.K.LE.NSW3)SFLX(3)=SFLX(3)+UVLEAN(LMO,K)*DSLEAN(K)
      IF(K.GT.NSW3.AND.K.LE.NSW4)SFLX(4)=SFLX(4)+UVLEAN(LMO,K)*DSLEAN(K)
                                 SFLX(5)=SFLX(5)+UVLEAN(LMO,K)*DSLEAN(K)
  340 CONTINUE
  350 CONTINUE
      DO 360 I=1,5
      SFLX(I)=SFLX(I)/mavg
      DFLX(I)=SFLX(I)-SFL0(I)
      RFLX(I)=SFLX(I)/SFL0(I)
  360 CONTINUE
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6301) JYEARS,(SFLX(I),I=1,5),(DFLX(I),I=1,5)
     +                     ,(RFLX(I),I=1,5)
 6301 FORMAT(2X,I4,1X,4F8.4,F10.4,2X,5F8.4,2X,5F8.5)
      NSPACE=JYEARS-(JYEARS/10)*10
      IF(KLIMIT.GT.0) GO TO 370
      IF(NSPACE.EQ.0) WRITE(KW,6302)
 6302 FORMAT(' ')
  370 CONTINUE
      GO TO 9999
C
C-------------
  400 CONTINUE
C-------------
C
      JJDAYO=JMONTH*30-15
      IF(JMONTH.LT.1) JJDAYO=JDAY
      CALL UPDO3D(JYRREF,JJDAYO)
      DO 450 J=1,46
      DO 410 L=1,NL
      O3(J,L)=0.D0
  410 CONTINUE
      JLAT=J
      DO 430 I=1,72
      ILON=I
!!!   CALL GETO3D(ILON,JLAT)
      CALL REPART(O3JDAY(1,ILON,JLAT),PLBO3,NLO3+1,U0GAS(1,3),PLB0,NL+1)
      DO 420 L=1,NL
      O3(J,L)=O3(J,L)+U0GAS(L,3)/72.D0
  420 CONTINUE
  430 CONTINUE
      SUMO3=0.D0
      DO 440 L=1,NL
      SUMO3=SUMO3+O3(J,L)
  440 CONTINUE
      O3COL(J)=SUMO3
  450 CONTINUE
      CALL BOXAV1(DLAT46,O3COL,46, 1,23,QOSH)
      CALL BOXAV1(DLAT46,O3COL,46,24,46,QONH)
      CALL BOXAV1(DLAT46,O3COL,46, 1,46,QOGL)
      O3COL(47)=QOSH
      O3COL(48)=QONH
      O3COL(49)=QOGL
      DO 460 L=1,NL
      CALL BOXAV1(DLAT46,O3(1,L),46, 1,23,QOSH)
      CALL BOXAV1(DLAT46,O3(1,L),46,24,46,QONH)
      CALL BOXAV1(DLAT46,O3(1,L),46, 1,46,QOGL)
      O3(47,L)=QOSH
      O3(48,L)=QONH
      O3(49,L)=QOGL
  460 CONTINUE
C
      IF(KLIMIT.GT.0)
     +WRITE(KW,6400) JYRREF,JJDAYO,JMONTH,MADO3M,(L,L=2,15)
 6400 FORMAT(/' (4)=INDEX  JYRREF=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone: Zonal-mean Vertical Distribution (cmSTP)'
     +      ,T126,'MADO3=',I2/'  JLAT DLAT46   COLUMN  L =   1',14I7)
      IF(KLIMIT.LT.1)
     +WRITE(KW,7400) JYRREF,JJDAYO,JMONTH,MADO3M
     +              ,(PLB0(I),I=1,15),(L,L=2,15)
 7400 FORMAT(/' (4)=INDEX  JYRREF=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone: Zonal-mean Vertical Distribution (cmSTP)'
     +      ,T126,'MADO3=',I2//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =   1',14I7)
C
      DO 470 JJ=1,46
      J=47-JJ
      IF(KLIMIT.GT.0) GO TO 470
      WRITE(KW,6401) J,DLAT46(J),O3COL(J),(O3(J,L),L=1,15)
 6401 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  470 CONTINUE
      IF(KLIMIT.LT.1) WRITE(KW,6402)
 6402 FORMAT(' ')
      WRITE(KW,6403) O3COL(48),(O3(48,L),L=1,15)
 6403 FORMAT(11X,'NH',F9.5,4X,15F7.5)
      IF(KLIMIT.LT.1) WRITE(KW,6402)
      WRITE(KW,6404) O3COL(47),(O3(47,L),L=1,15)
 6404 FORMAT(11X,'SH',F9.5,4X,15F7.5)
      IF(KLIMIT.LT.1) WRITE(KW,6402)
      WRITE(KW,6405) O3COL(49),(O3(49,L),L=1,15)
 6405 FORMAT( 7X,'GLOBAL',F9.5,4X,15F7.5)
      GO TO 9999
C
C
C-------------
  500 CONTINUE
C-------------
C
      JJDAYO=JMONTH*30-15
      IF(JMONTH.LT.1) JJDAYO=JDAY
      CALL UPDO3D(JYRREF,JJDAYO)
      DO 590 N=1,3
      N1=1
      N2=8
      IF(N.EQ.2) N1=9
      IF(N.GT.1) N2=NL
      DO 530 J=1,46
      JLAT=J
      DO 520 I=1,72
      ILON=I
!!!   CALL GETO3D(ILON,JLAT)
      CALL REPART(O3JDAY(1,ILON,JLAT),PLBO3,NLO3+1,U0GAS(1,3),PLB0,NL+1)
      SUMO3=0.D0
      DO 510 L=N1,N2
      SUMO3=SUMO3+U0GAS(L,3)
  510 CONTINUE
      O3L(J,I)=SUMO3
  520 CONTINUE
  530 CONTINUE
      DO 560 J=1,46
      SUMO3=0.D0
      DO 540 I=1,72
      SUMO3=SUMO3+O3L(J,I)/72.D0
  540 CONTINUE
      DO 550 I=1,72
      O3L(J,I)=O3L(J,I)/SUMO3
  550 CONTINUE
  560 CONTINUE
C
      IF(N.EQ.1)
     +WRITE(KW,6510) JYRREF,JJDAYO,JMONTH,MADO3M,(I,I=10,310,10)
 6510 FORMAT(/'  5A=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Troposphere'
     +      ,' (Wang-Jacobs) Surf to 150 mb',T126,'MADO3=',I2
     +      /'  J LON=0',31I4)
      IF(N.EQ.2)
     +WRITE(KW,6520) JYRREF,JJDAYO,JMONTH,MADO3M,(I,I=10,310,10)
 6520 FORMAT(/'  5B=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Stratosphere'
     +      ,' (London-NCAR) 150 mb to TOA',T126,'MADO3=',I2
     +      /'  J LON=0',31I4)
      IF(N.EQ.3.AND.KLIMIT.LT.1)
     +WRITE(KW,6530) JYRREF,JJDAYO,JMONTH,MADO3M,(I,I=10,310,10)
 6530 FORMAT(/'  5C=INDEX  JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' Ozone Longitudinal Variation:  Total Column'
     +      ,' (W-J/London) Surface to TOA',T126,'MADO3=',I2
     +      /'  J LON=0',31I4)
      IF(KLIMIT.LT.1) WRITE(KW,6540)
 6540 FORMAT(' ')
C
      IF(N.EQ.3.AND.KLIMIT.GT.0) GO TO 590
      DO 580 JJ=1,46
      J=47-JJ
      KWSKIP=KLIMIT
      IF(J.EQ.36) KWSKIP=0
      IF(J.EQ.24) KWSKIP=0
      IF(J.EQ.12) KWSKIP=0
      DO 570 I=1,36
      II=I*2-1
      LO3(I)=O3L(J,II)*100.D0+0.5D0
  570 CONTINUE
      IF(KWSKIP.EQ.0)
     +WRITE(KW,6501) J,(LO3(I),I=1,32)
 6501 FORMAT(I4,1X,36I4)
  580 CONTINUE
  590 CONTINUE
C
      GO TO 9999
C
C-------------
  600 CONTINUE
C-------------
C
      KAEROS=4
      IF(INDX.EQ.6) KAEROS=1
      IF(INDX.EQ.7) KAEROS=2
      IF(INDX.EQ.8) KAEROS=3
      L1=1
      IF(INDX.EQ.8.and.NL.gt.15) L1=NL-14
      JJDAY=JMONTH*30-15
      IF(JMONTH.LT.1) JJDAY=JDAY
      K=6
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) CALL UPDAER(JYRREF,JJDAY)
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) CALL UPDDST(JYRREF,JJDAY)
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) CALL UPDVOL(JYRREF,JJDAY)
C
      DO 650 J=1,46
      DO 610 L=1,NL
      QX(J,L)=0.D0
      QS(J,L)=0.D0
      QG(J,L)=0.D0
  610 CONTINUE
      JLAT=J
      DO 630 I=1,72
      ILON=I
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) CALL GETAER
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) CALL GETDST
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) CALL GETVOL
      DO 620 L=1,NL
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) QX(J,L)=QX(J,L)+SRAEXT(L,K)/72.D0
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) QX(J,L)=QX(J,L)+SRDEXT(L,K)/72.D0
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) QX(J,L)=QX(J,L)+SRVEXT(L,K)/72.D0
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) QS(J,L)=QS(J,L)+SRASCT(L,K)/72.D0
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) QS(J,L)=QS(J,L)+SRDSCT(L,K)/72.D0
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) QS(J,L)=QS(J,L)+SRVSCT(L,K)/72.D0
      IF(KAEROS.EQ.1.OR.KAEROS.GT.3) QG(J,L)=QG(J,L)+SRAGCB(L,K)
     +                                              *SRASCT(L,K)/72.D0
      IF(KAEROS.EQ.2.OR.KAEROS.GT.3) QG(J,L)=QG(J,L)+SRDGCB(L,K)
     +                                              *SRDSCT(L,K)/72.D0
      IF(KAEROS.EQ.3.OR.KAEROS.GT.3) QG(J,L)=QG(J,L)+SRVGCB(L,K)
     +                                              *SRVSCT(L,K)/72.D0
  620 CONTINUE
  630 CONTINUE
      SUMXL=1.D-10
      SUMSL=1.D-20
      SUMGL=1.D-20
      DO 640 L=1,NL
      SUMXL=SUMXL+QX(J,L)
      SUMSL=SUMSL+QS(J,L)
      SUMGL=SUMGL+QG(J,L)
      QG(J,L)=(1.D-20+QG(J,L))/(1.D-10+QS(J,L))
      QP(J,L)=(1.D-20+QS(J,L))/(1.D-10+QX(J,L))
      IF(QP(J,L).GT.0.99999D0) QP(J,L)=0.99999D0
  640 CONTINUE
      QXCOL(J)=SUMXL
      QSCOL(J)=SUMSL
      QGCOL(J)=(1.D-15+SUMGL)/(1.D-05+SUMSL)
      QPCOL(J)=(1.D-20+SUMSL)/(1.D-10+SUMXL)
  650 CONTINUE
      CALL BOXAV1(DLAT46,QXCOL,46, 1,23,QXSH)
      CALL BOXAV1(DLAT46,QXCOL,46,24,46,QXNH)
      CALL BOXAV1(DLAT46,QXCOL,46, 1,46,QXGL)
      QXCOL(47)=QXSH
      QXCOL(48)=QXNH
      QXCOL(49)=QXGL
      CALL BOXAV1(DLAT46,QSCOL,46, 1,23,QSSH)
      CALL BOXAV1(DLAT46,QSCOL,46,24,46,QSNH)
      CALL BOXAV1(DLAT46,QSCOL,46, 1,46,QSGL)
      QSCOL(47)=QSSH
      QSCOL(48)=QSNH
      QSCOL(49)=QSGL
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46, 1,23,QPSH)
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46,24,46,QPNH)
      CALL BOXAV2(DLAT46,QXCOL,QPCOL,46, 1,46,QPGL)
      QPCOL(47)=QPSH
      QPCOL(48)=QPNH
      QPCOL(49)=QPGL
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46, 1,23,QGSH)
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46,24,46,QGNH)
      CALL BOXAV2(DLAT46,QSCOL,QGCOL,46, 1,46,QGGL)
      QGCOL(47)=QGSH
      QGCOL(48)=QGNH
      QGCOL(49)=QGGL
      DO 660 L=1,NL
      CALL BOXAV1(DLAT46,QX(1,L),46, 1,23,QXSH)
      CALL BOXAV1(DLAT46,QX(1,L),46,24,46,QXNH)
      CALL BOXAV1(DLAT46,QX(1,L),46, 1,46,QXGL)
      QX(47,L)=QXSH
      QX(48,L)=QXNH
      QX(49,L)=QXGL
      CALL BOXAV1(DLAT46,QS(1,L),46, 1,23,QSSH)
      CALL BOXAV1(DLAT46,QS(1,L),46,24,46,QSNH)
      CALL BOXAV1(DLAT46,QS(1,L),46, 1,46,QSGL)
      QS(47,L)=QSSH
      QS(48,L)=QSNH
      QS(49,L)=QSGL
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46, 1,23,QPSH)
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46,24,46,QPNH)
      CALL BOXAV2(DLAT46,QX(1,L),QP(1,L),46, 1,46,QPGL)
      QP(47,L)=QPSH
      QP(48,L)=QPNH
      QP(49,L)=QPGL
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46, 1,23,QGSH)
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46,24,46,QGNH)
      CALL BOXAV2(DLAT46,QS(1,L),QG(1,L),46, 1,46,QGGL)
      QG(47,L)=QGSH
      QG(48,L)=QGNH
      QG(49,L)=QGGL
  660 CONTINUE
C
      KA=KAEROS
      IF(KLIMIT.GT.0)
     +WRITE(KW,6600) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA),(L,L=L1,L1+14)
 6600 FORMAT(/I3,'A=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL OPTICAL DEPTH',T100,A32
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
      IF(KLIMIT.LT.1)
     +WRITE(KW,7600) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=L1,L1+14),(L,L=L1,L1+14)
 7600 FORMAT(/I3,'A=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL OPTICAL DEPTH',T100,A32
     +      //21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
C
      IF(KLIMIT.LT.1) THEN
      DO 670 JJ=1,46
      J=47-JJ
      WRITE(KW,6601) J,DLAT46(J),QXCOL(J),(QX(J,L),L=L1,L1+14)
 6601 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  670 CONTINUE
      WRITE(KW,6602) QXCOL(48),(QX(48,L),L=L1,L1+14)
 6602 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6603) QXCOL(47),(QX(47,L),L=L1,L1+14)
 6603 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6604) QXCOL(49),(QX(49,L),L=L1,L1+14)
 6604 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
      ENDIF
      IF(KLIMIT.GT.0) THEN
      WRITE(KW,6605) QXCOL(48),(QX(48,L),L=L1,L1+14)
 6605 FORMAT( 11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6606) QXCOL(47),(QX(47,L),L=L1,L1+14)
 6606 FORMAT( 11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6607) QXCOL(49),(QX(49,L),L=L1,L1+14)
 6607 FORMAT( 7X,'GLOBAL',F9.5,4X,15F7.5)
      ENDIF
C
      IF(KLIMIT.GT.0) GO TO 699
      WRITE(KW,6610) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=L1,L1+14),(L,L=L1,L1+14)
 6610 FORMAT(/I3,'B=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL SINGLE SCATTERING ALBEDO'
     +      ,T100,A32//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
C
      DO 680 JJ=1,46
      J=47-JJ
      WRITE(KW,6611) J,DLAT46(J),QPCOL(J),(QP(J,L),L=L1,L1+14)
 6611 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  680 CONTINUE
      WRITE(KW,6612) QPCOL(48),(QP(48,L),L=L1,L1+14)
 6612 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6613) QPCOL(47),(QP(47,L),L=L1,L1+14)
 6613 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6614) QPCOL(49),(QP(49,L),L=L1,L1+14)
 6614 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
C
      WRITE(KW,6620) INDX,JYRREF,JJDAY,JMONTH,CHAER(KA)
     +              ,(PLB0(I),I=L1,L1+14),(L,L=L1,L1+14)
 6620 FORMAT(/I3,'C=INDEX JYEAR=',I5,'  JDAY=',I3,'   JMONTH=',I2
     +      ,T50,' ZONAL MEAN AEROSOL ASYMMETRY PARAMETER'
     +      ,T100,A32//21X,'PLB0 =',F6.1,9F7.1,5F7.2
     +      /'  JLAT DLAT46   COLUMN  L =',I4,14I7)
C
      DO 690 JJ=1,46
      J=47-JJ
      WRITE(KW,6621) J,DLAT46(J),QGCOL(J),(QG(J,L),L=L1,L1+14)
 6621 FORMAT(I5,F8.2,F9.5,4X,15F7.5)
  690 CONTINUE
      WRITE(KW,6622) QGCOL(48),(QG(48,L),L=L1,L1+14)
 6622 FORMAT(/11X,'NH',F9.5,4X,15F7.5)
      WRITE(KW,6623) QGCOL(47),(QG(47,L),L=L1,L1+14)
 6623 FORMAT(/11X,'SH',F9.5,4X,15F7.5)
      WRITE(KW,6624) QGCOL(49),(QG(49,L),L=L1,L1+14)
 6624 FORMAT(/7X,'GLOBAL',F9.5,4X,15F7.5)
  699 CONTINUE
      GO TO 9999
C
 1000 CONTINUE
C
 9999 CONTINUE
C
      RETURN
      END SUBROUTINE WRITET
      END MODULE RADPAR

      SUBROUTINE BOXAV1(DEGLAT,TAULAT,NLAT,JALIM,JBLIM,TAU)
      IMPLICIT NONE
C
C--------------------------------------------------------------------
C     BOXAV1  Performs:
C                       Latitudinal average (area-weighted) of TAULAT
C
C              DEGLAT   Center latitude of grid-box variable (TAULAT)
C                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
C
C              TAULAT   Zonal average value is constant over grid-bos
C
C        JALIM, JBLIM   Latitude boxes for which (TAULAT) is averaged
C
C                 TAU   Area-weighted (TAULAT) latitude average value
C--------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: NLAT, JALIM, JBLIM
      REAL*8, DIMENSION(NLAT), INTENT(IN) :: DEGLAT, TAULAT
      REAL*8, INTENT(OUT) :: TAU
      REAL*8 ASUM,TSUM,PI,RADIAN,RLAT1,ALAT1,RLAT2,ALAT2,ALATJ
      INTEGER J1,J,J2

      ASUM=0.D0
      TSUM=0.D0
      PI=ACOS(-1.D0)
      RADIAN=180.D0/PI
      J1=JALIM-1
      IF(J1.LT.1) J1=1
      RLAT1=(0.5D0*(DEGLAT(J1)+DEGLAT(JALIM))+90.D0)/RADIAN
      ALAT1=SIN(RLAT1)
      DO 110 J=JALIM,JBLIM
      J2=J+1
      IF(J2.GT.NLAT) J2=NLAT
      RLAT2=(0.5D0*(DEGLAT(J)+DEGLAT(J2))+90.D0)/RADIAN
      ALAT2=SIN(RLAT2)
      ALATJ=0.5D0*(ALAT1+ALAT2)/(RLAT2-RLAT1)
      ASUM=ASUM+ALATJ
      TSUM=TSUM+ALATJ*TAULAT(J)
      RLAT1=RLAT2
      ALAT1=ALAT2
  110 CONTINUE
      TAU=TSUM/ASUM
      RETURN
      END SUBROUTINE BOXAV1

      SUBROUTINE BOXAV2(DEGLAT,TAULAT,SIZLAT,NLAT,JALIM,JBLIM,SIZ)
      IMPLICIT NONE
C
C--------------------------------------------------------------------
C     BOXAV2  Performs:
C                       TAULAT-weighted latitudinal average of SIZLAT
C
C              DEGLAT   Center latitude of grid-box variable (TAULAT)
C                       of the form:  DEGLAT = -90+(J-1)*180/(NLAT-1)
C
C              TAULAT   Zonal average value is constant over grid-box
C              SIZLAT   Zonal average value is constant over grid-box
C
C        JALIM, JBLIM   Latitude boxes for which variable is averaged
C
C                 SIZ   TAULAT-weighted latitudinal average of SIZLAT
C--------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: NLAT,JALIM,JBLIM
      REAL*8, DIMENSION(NLAT), INTENT(IN) :: DEGLAT,TAULAT,SIZLAT
      REAL*8, INTENT(OUT) :: SIZ
      REAL*8 ASUM,TSUM,PI,RADIAN,RLAT1,RLAT2,ALAT1,ALAT2,ALATJ
      INTEGER J,J1,J2

      ASUM=0.D0
      TSUM=0.D0
      PI=ACOS(-1.D0)
      RADIAN=180.D0/PI
      J1=JALIM-1
      IF(J1.LT.1) J1=1
      RLAT1=(0.5D0*(DEGLAT(J1)+DEGLAT(JALIM))+90.D0)/RADIAN
      ALAT1=SIN(RLAT1)
      DO 110 J=JALIM,JBLIM
      J2=J+1
      IF(J2.GT.NLAT) J2=NLAT
      RLAT2=(0.5D0*(DEGLAT(J)+DEGLAT(J2))+90.D0)/RADIAN
      ALAT2=SIN(RLAT2)
      ALATJ=0.5D0*(ALAT1+ALAT2)/(RLAT2-RLAT1)
      ASUM=ASUM+ALATJ*TAULAT(J)
      TSUM=TSUM+ALATJ*TAULAT(J)*SIZLAT(J)
      RLAT1=RLAT2
      ALAT1=ALAT2
  110 CONTINUE
      SIZ=(1.D-20+TSUM)/(1.D-10+ASUM)
      RETURN
      END SUBROUTINE BOXAV2


      SUBROUTINE GTREND(XNOW,TNOW)
C
      USE RADPAR, only: nghg,ghgyr1,ghgyr2,ghgam
      IMPLICIT NONE
      real*8 xnow(nghg),tnow,year,dy,frac
      integer iy,n
C
C-------------------------------------------------------------
C        Makiko's GHG Trend Compilation  GHG.1850-2050.Dec1999
C
C        Annual-Mean      Greenhouse Gas Mixing Ratios
C-------------------------------------------------------------
C                 CO2     N2O     CH4   CFC-11  CFC-12  others
C        Year     ppm     ppm     ppm     ppb     ppb     ppb
C-------------------------------------------------------------
C     Read from external file - outside table: use value from
C                                      years ghgyr1 or ghgyr2
      YEAR=TNOW
      IF(TNOW.LE.ghgyr1+.5D0) YEAR=ghgyr1+.5D0
      IF(TNOW.gE.ghgyr2+.49999D0) YEAR=ghgyr2+.49999D0
      DY=YEAR-(ghgyr1+.5D0)
      IY=DY
      frac=DY-IY
      IY=IY+1
C
C     CO2 N2O CH4 CFC-11 CFC-12 other_GHG  SCENARIO
C--------------------------------------------------
C
      do n=1,nghg
        XNOW(N)=GHGAM(N,IY)+frac*(GHGAM(N,IY+1)-GHGAM(N,IY))
      end do
C
      RETURN
      END SUBROUTINE GTREND

      SUBROUTINE PHATMO(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     -------------     MCCLATCHY (1972) ATMOSPHERE DATA     -----------
C     ------------------------------------------------------------------
C
C        INPUT DATA
C------------------
C                  NATM=0  GIVES ABREVIATED DATA FOR  STANDARD ATMOSPHER
C                                (INPUT: P OR H) (RETURNS: H OR P & D,T)
C
C                  NATM=1  GIVES ATMOSPHERE DATA FOR  TROPICAL LATITUDES
C                  NATM=2  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE SUMMER
C                  NATM=3  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE WINTER
C                  NATM=4  GIVES ATMOSPHERE DATA FOR  SUBARCTIC SUMMER
C                  NATM=5  GIVES ATMOSPHERE DATA FOR  SUBARCTIC WINTER
C                  NATM=6  GIVES ATMOSPHERE DATA FOR  STANDARD ATMOSPHER
C
C                  NPHD=1  RETURNS H,D,T,O,Q,S DATA FOR GIVEN PRESSURE P
C                  NPHD=2  RETURNS P,D,T,O,Q,S DATA FOR GIVEN   HEIGHT H
C                  NPHD=3  RETURNS P,H,T,O,Q,S DATA FOR GIVEN  DENSITY D
C
C       OUTPUT DATA
C------------------
C                  P = PRESSURE IN MILLIBARS
C                  H = HEIGHT IN KILOMETERS
C                  D = DENSITY IN GRAMS/METER**3
C                  T = TEMPERATURE (ABSOLUTE)
C                  O = OZONE MIXING RATIO (GRAMS OZONE)/(GRAMS AIR)
C                  Q = SPECIFIC HUMIDITY (GRAMS WATER VAPOR)/(GRAMS AIR)
C                  S = SATURATION RATIO (GRAMS WATER VAPOR)/(GRAMS AIR)
C                  OCM = OZONE (CM-STP) ABOVE GIVEN HEIGHT
C                  WCM = WATER VAPOR (CM-STP) ABOVE GIVEN HEIGHT
C
C           REMARKS
C------------------
C                  INPUT P,H,D PARAMETERS ARE NOT ALTERED
C                  P,D INTERPOLATION IS EXPONENTIAL WITH HEIGHT
C                  NO EXTRAPOLATION IS MADE OUTSIDE 0-100 KM INTERVAL
C                  S  IS NOT COMPUTED ABOVE 40 KM (FORMULA NOT ACCURATE)
C
C                  R = Q/S          GIVES RELATIVE HUMIDITY
C                  W = Q/(1-Q)      GIVES WATER VAPOR MIXING RATIO
C                  N = D*2.079E 16  GIVES NUMBER DENSITY PER CM**3
C
      REAL*8, DIMENSION(33) ::
     *             PRS1,PRS2,PRS3,PRS4,PRS5,PRS6
     1            ,DNS1,DNS2,DNS3,DNS4,DNS5,DNS6
     2            ,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6
     3            ,WVP1,WVP2,WVP3,WVP4,WVP5,WVP6
     4            ,OZO1,OZO2,OZO3,OZO4,OZO5,OZO6
      REAL*8, DIMENSION(33,6) :: PRES,DENS,TEMP,WVAP,OZON

      EQUIVALENCE
     +     (PRES(1,1),PRS1(1)),(DENS(1,1),DNS1(1)),(TEMP(1,1),TMP1(1))
     +    ,(PRES(1,2),PRS2(1)),(DENS(1,2),DNS2(1)),(TEMP(1,2),TMP2(1))
     +    ,(PRES(1,3),PRS3(1)),(DENS(1,3),DNS3(1)),(TEMP(1,3),TMP3(1))
     +    ,(PRES(1,4),PRS4(1)),(DENS(1,4),DNS4(1)),(TEMP(1,4),TMP4(1))
     +    ,(PRES(1,5),PRS5(1)),(DENS(1,5),DNS5(1)),(TEMP(1,5),TMP5(1))
     +    ,(PRES(1,6),PRS6(1)),(DENS(1,6),DNS6(1)),(TEMP(1,6),TMP6(1))
      EQUIVALENCE (WVAP(1,1),WVP1(1)),(OZON(1,1),OZO1(1))
      EQUIVALENCE (WVAP(1,2),WVP2(1)),(OZON(1,2),OZO2(1))
      EQUIVALENCE (WVAP(1,3),WVP3(1)),(OZON(1,3),OZO3(1))
      EQUIVALENCE (WVAP(1,4),WVP4(1)),(OZON(1,4),OZO4(1))
      EQUIVALENCE (WVAP(1,5),WVP5(1)),(OZON(1,5),OZO5(1))
      EQUIVALENCE (WVAP(1,6),WVP6(1)),(OZON(1,6),OZO6(1))

      REAL*8, PARAMETER, DIMENSION(33) :: HTKM =
     *     (/1d-9, 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0,10d0,11d0
     1     ,12d0,13d0,14d0,15d0,16d0,17d0,18d0,19d0,20d0,21d0,22d0,23d0
     *     ,24d0,25d0,30d0,35d0,40d0,45d0,50d0,70d0,99.9d0/)

C----------------------------------------------------------------------
C0000 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS
C----------------------------------------------------------------------

      REAL*8, PARAMETER, DIMENSION(8) :: SPLB=(/
     *     1013.25d0, 226.32d0, 54.748d0 , 8.6801d0,
     *     1.109d0  , .66938d0, .039564d0, 3.7338d-03/),
     *     STLB=(/288.15d0, 216.65d0, 216.65d0, 228.65d0,
     *            270.65d0, 270.65d0, 214.65d0, 186.87d0/),
     *     SHLB=(/0d0,11d0,20d0,32d0,47d0,51d0,71d0,84.852d0/),
     *     SDLB=(/-6.5d0,0d0,1d0,2.8d0,0d0,-2.8d0,-2d0,0d0/)
      REAL*8, PARAMETER :: HPCON = 34.16319d0

C-----------------------------------------------------------------------
C1111 TROPICAL LATITUDES      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS1/      1.013D 03,9.040D 02,8.050D 02,7.150D 02,6.330D 02,
     1      5.590D 02,4.920D 02,4.320D 02,3.780D 02,3.290D 02,2.860D 02,
     2      2.470D 02,2.130D 02,1.820D 02,1.560D 02,1.320D 02,1.110D 02,
     3      9.370D 01,7.890D 01,6.660D 01,5.650D 01,4.800D 01,4.090D 01,
     4      3.500D 01,3.000D 01,2.570D 01,1.220D 01,6.000D 00,3.050D 00,
     5      1.590D 00,8.540D-01,5.790D-02,3.000D-04/
      DATA DNS1/      1.167D 03,1.064D 03,9.689D 02,8.756D 02,7.951D 02,
     1      7.199D 02,6.501D 02,5.855D 02,5.258D 02,4.708D 02,4.202D 02,
     2      3.740D 02,3.316D 02,2.929D 02,2.578D 02,2.260D 02,1.972D 02,
     3      1.676D 02,1.382D 02,1.145D 02,9.515D 01,7.938D 01,6.645D 01,
     4      5.618D 01,4.763D 01,4.045D 01,1.831D 01,8.600D 00,4.181D 00,
     5      2.097D 00,1.101D 00,9.210D-02,5.000D-04/
      DATA TMP1/  300.0,294.0,288.0,284.0,277.0,270.0,264.0,257.0,250.0,
     1244.0,237.0,230.0,224.0,217.0,210.0,204.0,197.0,195.0,199.0,203.0,
     2207.0,211.0,215.0,217.0,219.0,221.0,232.0,243.0,254.0,265.0,270.0,
     3  219.0,210.0/
      DATA WVP1/1.9D 01,1.3D 01,9.3D 00,4.7D 00,2.2D 00,1.5D 00,8.5D-01,
     1  4.7D-01,2.5D-01,1.2D-01,5.0D-02,1.7D-02,6.0D-03,1.8D-03,1.0D-03,
     2  7.6D-04,6.4D-04,5.6D-04,5.0D-04,4.9D-04,4.5D-04,5.1D-04,5.1D-04,
     3  5.4D-04,6.0D-04,6.7D-04,3.6D-04,1.1D-04,4.3D-05,1.9D-05,6.3D-06,
     4  1.4D-07,1.0D-09/
      DATA OZO1/5.6D-05,5.6D-05,5.4D-05,5.1D-05,4.7D-05,4.5D-05,4.3D-05,
     1  4.1D-05,3.9D-05,3.9D-05,3.9D-05,4.1D-05,4.3D-05,4.5D-05,4.5D-05,
     2  4.7D-05,4.7D-05,6.9D-05,9.0D-05,1.4D-04,1.9D-04,2.4D-04,2.8D-04,
     3  3.2D-04,3.4D-04,3.4D-04,2.4D-04,9.2D-05,4.1D-05,1.3D-05,4.3D-06,
     4  8.6D-08,4.3D-11/

C-----------------------------------------------------------------------
C2222 MIDLATITUDE SUMMER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS2/      1.013D 03,9.020D 02,8.020D 02,7.100D 02,6.280D 02,
     1      5.540D 02,4.870D 02,4.260D 02,3.720D 02,3.240D 02,2.810D 02,
     2      2.430D 02,2.090D 02,1.790D 02,1.530D 02,1.300D 02,1.110D 02,
     3      9.500D 01,8.120D 01,6.950D 01,5.950D 01,5.100D 01,4.370D 01,
     4      3.760D 01,3.220D 01,2.770D 01,1.320D 01,6.520D 00,3.330D 00,
     5      1.760D 00,9.510D-01,6.710D-02,3.000D-04/
      DATA DNS2/      1.191D 03,1.080D 03,9.757D 02,8.846D 02,7.998D 02,
     1      7.211D 02,6.487D 02,5.830D 02,5.225D 02,4.669D 02,4.159D 02,
     2      3.693D 02,3.269D 02,2.882D 02,2.464D 02,2.104D 02,1.797D 02,
     3      1.535D 02,1.305D 02,1.110D 02,9.453D 01,8.056D 01,6.872D 01,
     4      5.867D 01,5.014D 01,4.288D 01,1.322D 01,6.519D 00,3.330D 00,
     5      1.757D 00,9.512D-01,6.706D-02,5.000D-04/
      DATA TMP2/  294.0,290.0,285.0,279.0,273.0,267.0,261.0,255.0,248.0,
     1242.0,235.0,229.0,222.0,216.0,216.0,216.0,216.0,216.0,216.0,217.0,
     2218.0,219.0,220.0,222.0,223.0,224.0,234.0,245.0,258.0,270.0,276.0,
     3  218.0,210.0/
      DATA WVP2/1.4D 01,9.3D 00,5.9D 00,3.3D 00,1.9D 00,1.0D 00,6.1D-01,
     1  3.7D-01,2.1D-01,1.2D-01,6.4D-02,2.2D-02,6.0D-03,1.8D-03,1.0D-03,
     2  7.6D-04,6.4D-04,5.6D-04,5.0D-04,4.9D-04,4.5D-04,5.1D-04,5.1D-04,
     3  5.4D-04,6.0D-04,6.7D-04,3.6D-04,1.1D-04,4.3D-05,1.9D-05,6.3D-06,
     4  1.4D-07,1.0D-09/
      DATA OZO2/6.0D-05,6.0D-05,6.0D-05,6.2D-05,6.4D-05,6.6D-05,6.9D-05,
     1  7.5D-05,7.9D-05,8.6D-05,9.0D-05,1.1D-04,1.2D-04,1.5D-04,1.8D-04,
     2  1.9D-04,2.1D-04,2.4D-04,2.8D-04,3.2D-04,3.4D-04,3.6D-04,3.6D-04,
     3  3.4D-04,3.2D-04,3.0D-04,2.0D-04,9.2D-05,4.1D-05,1.3D-05,4.3D-06,
     4  8.6D-08,4.3D-11/

C-----------------------------------------------------------------------
C3333 MIDLATITUDE WINTER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS3/      1.018D 03,8.973D 02,7.897D 02,6.938D 02,6.081D 02,
     1      5.313D 02,4.627D 02,4.016D 02,3.473D 02,2.992D 02,2.568D 02,
     2      2.199D 02,1.882D 02,1.610D 02,1.378D 02,1.178D 02,1.007D 02,
     3      8.610D 01,7.350D 01,6.280D 01,5.370D 01,4.580D 01,3.910D 01,
     4      3.340D 01,2.860D 01,2.430D 01,1.110D 01,5.180D 00,2.530D 00,
     5      1.290D 00,6.820D-01,4.670D-02,3.000D-04/
      DATA DNS3/      1.301D 03,1.162D 03,1.037D 03,9.230D 02,8.282D 02,
     1      7.411D 02,6.614D 02,5.886D 02,5.222D 02,4.619D 02,4.072D 02,
     2      3.496D 02,2.999D 02,2.572D 02,2.206D 02,1.890D 02,1.620D 02,
     3      1.388D 02,1.188D 02,1.017D 02,8.690D 01,7.421D 01,6.338D 01,
     4      5.415D 01,4.624D 01,3.950D 01,1.783D 01,7.924D 00,3.625D 00,
     5      1.741D 00,8.954D-01,7.051D-02,5.000D-04/
      DATA TMP3/  272.2,268.7,265.2,261.7,255.7,249.7,243.7,237.7,231.7,
     1225.7,219.7,219.2,218.7,218.2,217.7,217.2,216.7,216.2,215.7,215.2,
     2215.2,215.2,215.2,215.2,215.2,215.2,217.4,227.8,243.2,258.5,265.7,
     3  230.7,210.2/
      DATA WVP3/3.5D 00,2.5D 00,1.8D 00,1.2D 00,6.6D-01,3.8D-01,2.1D-01,
     1  8.5D-02,3.5D-02,1.6D-02,7.5D-03,6.9D-03,6.0D-03,1.8D-03,1.0D-03,
     2  7.6D-04,6.4D-04,5.6D-04,5.0D-04,4.9D-04,4.5D-04,5.1D-04,5.1D-04,
     3  5.4D-04,6.0D-04,6.7D-04,3.6D-04,1.1D-04,4.3D-05,1.9D-05,6.3D-06,
     4  1.4D-07,1.0D-09/
      DATA OZO3/6.0D-05,5.4D-05,4.9D-05,4.9D-05,4.9D-05,5.8D-05,6.4D-05,
     1  7.7D-05,9.0D-05,1.2D-04,1.6D-04,2.1D-04,2.6D-04,3.0D-04,3.2D-04,
     2  3.4D-04,3.6D-04,3.9D-04,4.1D-04,4.3D-04,4.5D-04,4.3D-04,4.3D-04,
     3  3.9D-04,3.6D-04,3.4D-04,1.9D-04,9.2D-05,4.1D-05,1.3D-05,4.3D-06,
     4  8.6D-08,4.3D-11/

C-----------------------------------------------------------------------
C4444 SUBARCTIC SUMMER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS4/      1.010D 03,8.960D 02,7.929D 02,7.000D 02,6.160D 02,
     1      5.410D 02,4.730D 02,4.130D 02,3.590D 02,3.107D 02,2.677D 02,
     2      2.300D 02,1.977D 02,1.700D 02,1.460D 02,1.250D 02,1.080D 02,
     3      9.280D 01,7.980D 01,6.860D 01,5.890D 01,5.070D 01,4.360D 01,
     4      3.750D 01,3.227D 01,2.780D 01,1.340D 01,6.610D 00,3.400D 00,
     5      1.810D 00,9.870D-01,7.070D-02,3.000D-04/
      DATA DNS4/      1.220D 03,1.110D 03,9.971D 02,8.985D 02,8.077D 02,
     1      7.244D 02,6.519D 02,5.849D 02,5.231D 02,4.663D 02,4.142D 02,
     2      3.559D 02,3.059D 02,2.630D 02,2.260D 02,1.943D 02,1.671D 02,
     3      1.436D 02,1.235D 02,1.062D 02,9.128D 01,7.849D 01,6.750D 01,
     4      5.805D 01,4.963D 01,4.247D 01,1.338D 01,6.614D 00,3.404D 00,
     5      1.817D 00,9.868D-01,7.071D-02,5.000D-04/
      DATA TMP4/  287.0,282.0,276.0,271.0,266.0,260.0,253.0,246.0,239.0,
     1232.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,
     2225.0,225.0,225.0,225.0,226.0,228.0,235.0,247.0,262.0,274.0,277.0,
     3  216.0,210.0/
      DATA WVP4/9.1D 00,6.0D 00,4.2D 00,2.7D 00,1.7D 00,1.0D 00,5.4D-01,
     1  2.9D-01,1.3D-02,4.2D-02,1.5D-02,9.4D-03,6.0D-03,1.8D-03,1.0D-03,
     2  7.6D-04,6.4D-04,5.6D-04,5.0D-04,4.9D-04,4.5D-04,5.1D-04,5.1D-04,
     3  5.4D-04,6.0D-04,6.7D-04,3.6D-04,1.1D-04,4.3D-05,1.9D-05,6.3D-06,
     4  1.4D-07,1.0D-09/
      DATA OZO4/4.9D-05,5.4D-05,5.6D-05,5.8D-05,6.0D-05,6.4D-05,7.1D-05,
     1  7.5D-05,7.9D-05,1.1D-04,1.3D-04,1.8D-04,2.1D-04,2.6D-04,2.8D-04,
     2  3.2D-04,3.4D-04,3.9D-04,4.1D-04,4.1D-04,3.9D-04,3.6D-04,3.2D-04,
     3  3.0D-04,2.8D-04,2.6D-04,1.4D-04,9.2D-05,4.1D-05,1.3D-05,4.3D-06,
     4  8.6D-08,4.3D-11/

C-----------------------------------------------------------------------
C5555 SUBARCTIC WINTER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT
C-----------------------------------------------------------------------

      DATA PRS5/      1.013D 03,8.878D 02,7.775D 02,6.798D 02,5.932D 02,
     1      5.158D 02,4.467D 02,3.853D 02,3.308D 02,2.829D 02,2.418D 02,
     2      2.067D 02,1.766D 02,1.510D 02,1.291D 02,1.103D 02,9.431D 01,
     3      8.058D 01,6.882D 01,5.875D 01,5.014D 01,4.277D 01,3.647D 01,
     4      3.109D 01,2.649D 01,2.256D 01,1.020D 01,4.701D 00,2.243D 00,
     5      1.113D 00,5.719D-01,4.016D-02,3.000D-04/
      DATA DNS5/      1.372D 03,1.193D 03,1.058D 03,9.366D 02,8.339D 02,
     1      7.457D 02,6.646D 02,5.904D 02,5.226D 02,4.538D 02,3.879D 02,
     2      3.315D 02,2.834D 02,2.422D 02,2.071D 02,1.770D 02,1.517D 02,
     3      1.300D 02,1.113D 02,9.529D 01,8.155D 01,6.976D 01,5.966D 01,
     4      5.100D 01,4.358D 01,3.722D 01,1.645D 01,7.368D 00,3.330D 00,
     5      1.569D 00,7.682D-01,5.695D-02,5.000D-04/
      DATA TMP5/  257.1,259.1,255.9,252.7,247.7,240.9,234.1,227.3,220.6,
     1217.2,217.2,217.2,217.2,217.2,217.2,217.2,216.6,216.0,215.4,214.8,
     2214.1,213.6,213.0,212.4,211.8,211.2,216.0,222.2,234.7,247.0,259.3,
     3  245.7,210.0/
      DATA WVP5/1.2D 00,1.2D 00,9.4D-01,6.8D-01,4.1D-01,2.0D-01,9.8D-02,
     1  5.4D-02,1.1D-02,8.4D-03,5.5D-03,3.8D-03,2.6D-03,1.8D-03,1.0D-03,
     2  7.6D-04,6.4D-04,5.6D-04,5.0D-04,4.9D-04,4.5D-04,5.1D-04,5.1D-04,
     3  5.4D-04,6.0D-04,6.7D-04,3.6D-04,1.1D-04,4.3D-05,1.9D-05,6.3D-06,
     4  1.4D-07,1.0D-09/
      DATA OZO5/4.1D-05,4.1D-05,4.1D-05,4.3D-05,4.5D-05,4.7D-05,4.9D-05,
     1  7.1D-05,9.0D-05,1.6D-04,2.4D-04,3.2D-04,4.3D-04,4.7D-04,4.9D-04,
     2  5.6D-04,6.2D-04,6.2D-04,6.2D-04,6.0D-04,5.6D-04,5.1D-04,4.7D-04,
     3  4.3D-04,3.6D-04,3.2D-04,1.5D-04,9.2D-05,4.1D-05,1.3D-05,4.3D-06,
     4  8.6D-08,4.3D-11/

C----------------------------------------------------------------------
C6666 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS
C----------------------------------------------------------------------

      DATA PRS6/    1.01325D+03,8.987D+02,7.950D+02,7.011D+02,6.164D+02,
     1      5.402D+02,4.718D+02,4.106D+02,3.560D+02,3.074D+02,2.644D+02,
     2      2.263D+02,1.933D+02,1.651D+02,1.410D+02,1.204D+02,1.029D+02,
     3      8.787D+01,7.505D+01,6.410D+01,5.475D+01,4.678D+01,4.000D+01,
     4      3.422D+01,2.931D+01,2.511D+01,1.172D+01,5.589D+00,2.775D+00,
     5      1.431D+00,7.594D-01,4.634D-02,2.384D-04/
      DATA DNS6/      1.225D+03,1.112D+03,1.006D+03,9.091D+02,8.191D+02,
     1      7.361D+02,6.597D+02,5.895D+02,5.252D+02,4.663D+02,4.127D+02,
     2      3.639D+02,3.108D+02,2.655D+02,2.268D+02,1.937D+02,1.654D+02,
     3      1.413D+02,1.207D+02,1.031D+02,8.803D+01,7.487D+01,6.373D+01,
     4      5.428D+01,4.627D+01,3.947D+01,1.801D+01,8.214D+00,3.851D+00,
     5      1.881D+00,9.775D-01,7.424D-02,4.445D-04/
      DATA TMP6/
     1         288.150,281.650,275.150,268.650,262.150,255.650,249.150,
     2         242.650,236.150,229.650,223.150,216.650,216.650,216.650,
     3         216.650,216.650,216.650,216.650,216.650,216.650,216.650,
     4         217.650,218.650,219.650,220.650,221.650,226.650,237.050,
     5         251.050,265.050,270.650,217.450,186.870/
      DATA WVP6/      1.083D+01,6.323D+00,3.612D+00,2.015D+00,1.095D+00,
     1      5.786D-01,2.965D-01,1.469D-01,7.021D-02,3.226D-02,1.419D-02,
     2      5.956D-03,5.002D-03,4.186D-03,3.490D-03,2.896D-03,2.388D-03,
     3      1.954D-03,1.583D-03,1.267D-03,9.967D-04,8.557D-04,7.104D-04,
     4      5.600D-04,4.037D-04,2.406D-04,5.404D-05,2.464D-05,1.155D-05,
     5      5.644D-06,2.932D-06,2.227D-07,1.334D-09/
      DATA OZO6/      7.526D-05,3.781D-05,6.203D-05,3.417D-05,5.694D-05,
     1      3.759D-05,5.970D-05,4.841D-05,7.102D-05,6.784D-05,9.237D-05,
     2      9.768D-05,1.251D-04,1.399D-04,1.715D-04,1.946D-04,2.300D-04,
     3      2.585D-04,2.943D-04,3.224D-04,3.519D-04,3.714D-04,3.868D-04,
     4      3.904D-04,3.872D-04,3.728D-04,2.344D-04,9.932D-05,3.677D-05,
     5      1.227D-05,4.324D-06,5.294D-08,1.262D-10/

      REAL*8, INTENT(INOUT) :: H,P,D
      INTEGER, INTENT(IN) :: NATM,NPHD
      REAL*8, INTENT(OUT) :: O,Q,S,OCM,WCM,T
      REAL*8 :: XX,XI,XJ,DELTA,RAT,PI,PJ,DI,DJ,DP,ES,RS,OI,OJ,QI,QJ
      INTEGER :: I,J,K,N

      IF(NATM.GT.0) GO TO 200
      O=1.E-10
      Q=1.E-10
      S=1.E-10
      OCM=1.E-10
      WCM=1.E-10
      IF(NPHD.LT.2) GO TO 150
      DO 110 N=2,8
      IF(H.LT.SHLB(N)) GO TO 120
  110 CONTINUE
      N=9
  120 CONTINUE
      N=N-1
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 130
      P=SPLB(N)*(1.+SDLB(N)/STLB(N)*(H-SHLB(N)))**(-HPCON/SDLB(N))
      GO TO 140
  130 CONTINUE
      P=SPLB(N)*EXP(-HPCON/STLB(N)*(H-SHLB(N)))
  140 CONTINUE
      T=STLB(N)+SDLB(N)*(H-SHLB(N))
      D=P/T*28.9644E 05/8.31432E 03    ! P/T*mair/gasc
      RETURN

  150 CONTINUE
      DO 160 N=2,8
      IF(P.GT.SPLB(N)) GO TO 170
  160 CONTINUE
      N=9
  170 CONTINUE
      N=N-1
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 180
      H=SHLB(N)+STLB(N)/SDLB(N)*((SPLB(N)/P)**(SDLB(N)/HPCON)-1.)
      GO TO 190
  180 CONTINUE
      H=SHLB(N)+STLB(N)/HPCON*LOG(SPLB(N)/P)
  190 CONTINUE
      T=STLB(N)+SDLB(N)*(H-SHLB(N))
      D=P/T*28.9644E 05/8.31432E 03
      RETURN

 200  CONTINUE
      IF(NPHD.EQ.1) GO TO 240
      IF(NPHD.EQ.2) GO TO 220
      XX=D
      XI=DENS(1,NATM)
      IF(D.GT.XI) XX=XI
      IF(D.LT.5.0E-04) GO TO 280
      DO 210 J=2,33
      XJ=DENS(J,NATM)
      IF(XX.GT.XJ) GO TO 260
      XI=XJ
  210 CONTINUE
  220 CONTINUE
      XX=H
      XI=HTKM(1)
      IF(H.LT.XI) XX=XI
      IF(H.GT.99.9) GO TO 280
      DO 230 J=2,33
      XJ=HTKM(J)
      IF(XX.LT.XJ) GO TO 260
      XI=XJ
  230 CONTINUE
  240 CONTINUE
      XX=P
      XI=PRES(1,NATM)
      IF(P.GT.XI) XX=XI
      IF(P.LT.3.0E-04) GO TO 280
      DO 250 J=2,33
      XJ=PRES(J,NATM)
      IF(XX.GT.XJ) GO TO 260
      XI=XJ
  250 CONTINUE
  260 CONTINUE
      DELTA=(XX-XI)/(XJ-XI)
      I=J-1
      IF(NPHD.NE.2) THEN
      RAT=LOG(XX/XI)/LOG(XJ/XI)
      H=HTKM(I)+(HTKM(J)-HTKM(I))*RAT
      T=TEMP(I,NATM)+(TEMP(J,NATM)-TEMP(I,NATM))*RAT
      ENDIF
      PI=PRES(I,NATM)
      PJ=PRES(J,NATM)
      DI=DENS(I,NATM)
      DJ=DENS(J,NATM)
      IF(NPHD.EQ.1) D=DI+DELTA*(DJ-DI)
      IF(NPHD.EQ.3) P=PI+DELTA*(PJ-PI)
      IF(NPHD.EQ.2) THEN
      P=PI*(PJ/PI)**DELTA
      D=DI*(DJ/DI)**DELTA
      T=TEMP(I,NATM)+DELTA*(TEMP(J,NATM)-TEMP(I,NATM))
      ENDIF
      O=OZON(I,NATM)/DI+DELTA*(OZON(J,NATM)/DJ-OZON(I,NATM)/DI)
      Q=WVAP(I,NATM)/DI+DELTA*(WVAP(J,NATM)/DJ-WVAP(I,NATM)/DI)
      ES=10.D0**(9.4051D0-2353.D0/T)
      IF(P.LT.PI) PI=P
      S=1.D+06
      RS=(PI-ES+0.622*ES)/(0.622*ES)
      IF(RS.GT.1.E-06) S=1./RS
      OI=O
      QI=Q
      OCM=0.D0
      WCM=0.D0
      DO 270 K=J,33
      PJ=PRES(K,NATM)
      DJ=DENS(K,NATM)
      OJ=OZON(K,NATM)/DJ
      QJ=WVAP(K,NATM)/DJ
      DP=PI-PJ
      OCM=OCM+0.5D0*(OI+OJ)*DP
      WCM=WCM+0.5D0*(QI+QJ)*DP
      OI=OJ
      QI=QJ
      PI=PJ
  270 CONTINUE
      WCM=WCM/0.980D0*22420.7D0/18.D0
      OCM=OCM/0.980D0*22420.7D0/48.D0
      RETURN
  280 CONTINUE
      T=210.D0
      IF(NATM.EQ.6) T=186.87
      O=1.D-10
      Q=1.D-10
      S=1.D-10
      OCM=1.D-10
      WCM=1.D-10
      IF(NPHD.NE.1) P=1.D-05
      IF(NPHD.NE.2) H=99.99
      IF(NPHD.NE.3) D=2.D-05
      RETURN
      END SUBROUTINE PHATMO

      REAL*8 FUNCTION PFOFTK(WAVNA,WAVNB,TK)
C     ------------------------------------------------------------------
C
C        INPUT DATA
C                  WAVNA,WAVNB  SPECTRAL INTERVAL IN WAVENUMBERS
C                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)
C
C                  TK           ABSOLUTE TEMPERATURE IN DEGREES KELVIN
C
C       OUTPUT DATA
C                  PFOFTK       PLANCK FLUX (W/M**2)
C
C
C           REMARKS
C                   PLANCK INTENSITY (W/M**2/STER) IS GIVEN BY PFOFTK/PI
C
C     ------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, PARAMETER, DIMENSION(21) ::
     *     BN = (/1D0, -1D0, 1D0, -1D0, 1D0, -1D0, 5D0, -691D0, 7D0,
     1     -3617D0, 43867D0, -174611D0, 854513D0, -236364091D0,
     2     8553103D0, -23749461029D0, 8615841276005D0, -7709321041217D0,
     *     2577687858367D0, -2631527155305348D4, 2929993913841559D0 /),
     *     BD = (/1D0, 2D0, 6D0, 30D0, 42D0, 30D0, 66D0, 2730D0, 6D0
     1     ,510D0, 798D0, 330D0, 138D0, 2730D0, 6D0, 870D0, 14322D0
     2     ,510D0, 6D0, 1919190D0, 6D0/)
      REAL*8, PARAMETER :: PI4 =97.40909103400244D0
C     REAL*8, PARAMETER :: PI =3.141592653589793D0
      REAL*8, PARAMETER :: HCK=1.43879D0
      REAL*8, PARAMETER :: DGXLIM=1D-06

      REAL*8, INTENT(IN) :: WAVNA, WAVNB, TK
      REAL*8 GSUM,B,DG,DGB,GX,PNORM,GTERM,GXA,GXB,X,XX,XN,XN3,XNN,XNM
     *     ,XNF,XNX
      INTEGER II,NB,NNB,N

      PFOFTK=0D0
      IF(TK.LT.1D-06) RETURN
      DO 160 II=1,2
      IF(II.EQ.1) X=HCK*WAVNA/TK
      IF(II.EQ.2) X=HCK*WAVNB/TK
      IF(X.GT.2.3D0) GO TO 120
      XX=X*X
      GSUM=1D0/3D0-X/8D0+XX/60D0
      NB=3
      XNF=XX/2D0
      DO 100 N=4,38,2
      NB=NB+1
      NNB=NB
      B=BN(NB)/BD(NB)
      XN3=N+3
      XNM=N*(N-1)
      XNF=XNF*(XX/XNM)
      DG=B/XN3*XNF
      GSUM=GSUM+DG
      DGB=DG
      IF(ABS(DG).LT.DGXLIM) GO TO 110
  100 CONTINUE
  110 CONTINUE
      GX=GSUM*XX*X
      GO TO 150
  120 CONTINUE
      GSUM=PI4/15.D0
      DO 130 N=1,20
      NNB=N
      XN=N
      XNN=XN*XN
      XNX=XN*X
      IF(XNX.GT.100.D0) GO TO 140
      GTERM=(X*X*(3.D0+XNX)+6.D0*(1.D0+XNX)/XNN)/XNN
      DG=GTERM*EXP(-XNX)
      GSUM=GSUM-DG
      DGB=DG
      IF(DG.LT.DGXLIM) GO TO 140
  130 CONTINUE
  140 CONTINUE
      GX=GSUM
  150 CONTINUE
      IF(II.EQ.1) GXA=GX
      IF(II.EQ.2) GXB=GX
  160 CONTINUE
      PNORM=15.D0/PI4
      PFOFTK=ABS(GXB-GXA)*PNORM
      PFOFTK=PFOFTK*5.6692D-08*TK**4
      RETURN
      END FUNCTION PFOFTK

      REAL*8 FUNCTION TKOFPF(WAVNA,WAVNB,FLUXAB)
C     ------------------------------------------------------------------
C
C        INPUT DATA
C------------------
C                  WAVNA,WAVNB  SPECTRAL INTERVAL IN WAVENUMBERS
C                               (ORDER OF WAVNA,WAVNB NOT IMPORTANT)
C                  FLUXAB       PLANCK FLUX (W/M**2) IN INTERVAL
C                                                       (WAVNA,WAVNB)
C
C       OUTPUT DATA
C------------------
C                  TK           BRIGHTNESS TEMPERATURE IN DEGREES KELVIN
C
C
C           REMARKS
C------------------
C                   TKOFPF IS INVERSE FUNCTION OF PFOFTK(WAVNA,WAVNB,TK)
C                   THE OUTPUT OF TKOFPF SATISFIES THE IDENTITY
C                                 FLUXAB=PFOFTK(WAVNA,WAVNB,TK)
C                   (UNITS FOR FLUXAB AND PFOFTK MUST BE IDENTICAL)
C
C     ------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL LOGFIT
      REAL*8, PARAMETER :: DELFIT=1.D-06
      INTEGER, PARAMETER :: NMAX=20
      REAL*8, INTENT(IN) :: WAVNA,WAVNB,FLUXAB
      REAL*8 XA,YA,XB,YB,XX,YY,XC,YC,XBA,XCA,XBC,YBA,YCA,YBC,YXBA,YXCA,A
     *     ,B,C,ROOT,PF,PFOFTK
      INTEGER NFIT

      IF(FLUXAB.LE.0.D0) RETURN
      LOGFIT=.FALSE.
      NFIT=0
      PF=FLUXAB
      XA=0.D0
      YA=0.D0
      XB=250.D0
      YB=PFOFTK(WAVNA,WAVNB,XB)
      XX=PF*XB/YB
      YY=PFOFTK(WAVNA,WAVNB,XX)
      IF(ABS(YY-PF).LT.DELFIT) GO TO 200
      IF((YY/PF).LT.0.5D0) GO TO 150
      IF((YY/PF).GT.2.0D0) GO TO 170
      IF(XX.GT.XB) GO TO 110
      XC=XB
      YC=YB
      XB=XX
      YB=YY
      GO TO 120
  110 CONTINUE
      XC=XX
      YC=YY
  120 CONTINUE
      XBA=XB-XA
      XCA=XC-XA
      XBC=XB-XC
      YBA=YB-YA
      YCA=YC-YA
      YBC=YB-YC
      NFIT=NFIT+1
      IF(NFIT.GT.NMAX) GO TO 200
      YXBA=YBA/XBA
      YXCA=YCA/XCA
      C=(YXBA-YXCA)/XBC
      B=YXBA-(XB+XA)*C
      A=YA-XA*(B+XA*C)
      ROOT=SQRT(B*B+4.D0*C*(PF-A))
      XX=0.5D0*(ROOT-B)/C
      IF(XX.LT.XA.OR.XX.GT.XC) XX=-0.5D0*(ROOT+B)/C
      YY=PFOFTK(WAVNA,WAVNB,XX)
      IF(LOGFIT) YY=LOG(YY)
      IF(ABS(YY-PF).LT.DELFIT) GO TO 200
      IF(XX.GT.XB) GO TO 130
      XC=XB
      YC=YB
      GO TO 140
  130 CONTINUE
      XA=XB
      YA=YB
  140 CONTINUE
      XB=XX
      YB=YY
      GO TO 120
  150 CONTINUE
      XA=XX
      YA=YY
  160 CONTINUE
      XC=XB
      YC=YB
      XB=XB/2.D0
      YB=PFOFTK(WAVNA,WAVNB,XB)
      IF(YB.LT.YA) GO TO 190
      IF(YB.GT.PF) GO TO 160
      XA=XB
      YA=YB
      GO TO 190
  170 CONTINUE
      XC=XX
      YC=YY
  180 CONTINUE
      XA=XB
      YA=YB
      XB=XB*2.D0
      YB=PFOFTK(WAVNA,WAVNB,XB)
      IF(YB.GT.YC) GO TO 190
      IF(YB.LT.PF) GO TO 180
      XC=XB
      YC=YB
  190 CONTINUE
      XB=XA+(PF-YA)*(XC-XA)/(YC-YA)
      YB=PFOFTK(WAVNA,WAVNB,XB)
      XX=XB
      IF(ABS(YB-PF).LT.DELFIT) GO TO 200
      PF=LOG(PF)
      YA=LOG(YA)
      YB=LOG(YB)
      YC=LOG(YC)
      LOGFIT=.TRUE.
      GO TO 120
  200 CONTINUE
      TKOFPF=XX
      RETURN
      END FUNCTION TKOFPF

      SUBROUTINE SPLINE(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT NONE

      integer, intent(in)  :: NXF,              KXTRAP
      real*8 , intent(in)  :: X(NXF),F(NXF),XX, CUSPWM,CUSPWE
      real*8 , intent(out) :: FF

C---------------------------------------------------------------------
C
C    SPLINE locates XX between points (F2,X2)(F3,X3) on 4-point spread
C       and returns 4-point Cubic Spline interpolated value FF = F(XX)
C
C    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
C    (X-Coordinate may be specified in increasing or decreasing order)
C
C---------------------------------------------------------------------
C
C    Custom Control Parameters:  CUSPWM,CUSPWE,KXTRAP
C------------------------------
C
C    In cases where data points are unevenly spaced and/or data points
C    exhibit abrupt changes in value, Spline Interpolation may produce
C    undesirable bulging of interpolated values. In more extreme cases
C    Linear Interpolation may be less problematic to use.
C
C    Interpolation can be weighted between: Cubic Spline and Linear by
C    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
C
C    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
C    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
C
C    For example, with:
C
C    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
C    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
C
C---------------------------------------------------------------------
C
C     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
C
C               KXTRAP = 0    No Extrapolation  (i.e., sets F(XX)=0.0)
C                        1    Fixed Extrapolation (F(XX) = edge value)
C                        2    Linear Extrapolation using 2 edge points
C
C---------------------------------------------------------------------

      REAL*8 x1,x2,x3,x4, x21,x32,x43,x31,x42, betw,FFCUSP,FFLINR,CUSPWT
      REAL*8 f1,f2,f3,f4, f21,f32,f43,f3221,f4332,  a,b,c,d, xf,xe,xexm
      INTEGER K

      K=2
      X2=X(K)
      X3=X(NXF-1)
      BETW=(XX-X2)*(X3-XX)
      IF(BETW.LE.0.D0) GO TO 120

  100 CONTINUE
      K=K+1
      X3=X(K)
      BETW=(XX-X2)*(X3-XX)
      IF(BETW.GE.0.D0) GO TO 110
      X2=X3
      GO TO 100

  110 CONTINUE
      F3=F(K)
      F4=F(K+1)
      X4=X(K+1)
      F2=F(K-1)
      X2=X(K-1)
      F1=F(K-2)
      X1=X(K-2)
      X21=X2-X1
      X31=X3-X1
      X32=X3-X2
      X43=X4-X3
      X42=X4-X2
      F21=(F2-F1)/(X21*X21)
      F32=(F3-F2)/(X32*X32)
      F43=(F4-F3)/(X43*X43)
      F3221=(F32+F21)/X31*X21
      F4332=(F43+F32)/X42*X43
      A=F2
      B=X32*F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)/X32
      XF=XX-X2

C                             FFCUSP= Cubic Spline Interpolation Result
C                             -----------------------------------------

      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=(X3+X2-XX-XX)/X32
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE

C                                   FFLINR= Linear Interpolation Result
C                                   -----------------------------------
      FFLINR=A+XF*F32*X32
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

C                Edge Point Interval Interpolation and/or Extrapolation
C                ------------------------------------------------------
  120 CONTINUE
      BETW=(X2-XX)*(X3-X2)
      IF(BETW.LT.0.D0) GO TO 140

C                          X(1),X(2)  Edge Point Interval Interpolation
C                          --------------------------------------------
      X1=X(1)
      F1=F(1)
      F2=F(2)
      X21=X2-X1
      F21=(F2-F1)/X21
      XF=XX-X1
      BETW=(X2-XX)*XF
      IF(BETW.LT.0.D0) GO TO 130
      F3=F(3)
      X3=X(3)
      X32=X3-X2
      X31=X3-X1
      C=((F3-F2)/X32-F21)/X31
      B=F21-X21*C
      A=F1
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F21
      XE=1.D0-2.D0*XF/X21
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

  130 CONTINUE
C                  Extrapolation for XX Outside of Interval X(1) - X(2)
C                  ----------------------------------------------------
C                  IF(KXTRAP.EQ.0)  (No Extrapolation:  sets F(XX)=0.0)
C                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP.EQ.0) FF=0.D0
      IF(KXTRAP.EQ.1) FF=F1
      IF(KXTRAP.EQ.2) FF=F1+XF*F21
      GO TO 160

  140 CONTINUE
C                    X(NXF-1),X(NXF)  Edge Point Interval Interpolation
C                    --------------------------------------------------
      F3=F(NXF)
      X3=X(NXF)
      F2=F(NXF-1)
      X2=X(NXF-1)
      X32=X3-X2
      F32=(F3-F2)/X32
      XF=XX-X3
      BETW=(X2-XX)*(XX-X3)
      IF(BETW.LT.0.D0) GO TO 150
      F1=F(NXF-2)
      X1=X(NXF-2)
      X21=X2-X1
      X31=X3-X1
      F21=(F2-F1)/X21
      XF=XX-X2

C                    3-Point Quadratic Interpolation for Edge Intervals
C                    --------------------------------------------------
C
C      (Edge Option)     ----------------------------------------------
C                        For Linear Interpolation within Edge Intervals
C                        between X(1),X(2), and between X(NXF-1),X(NXF)
C                        set the value of coefficient C below, to C=0.0
C                        ----------------------------------------------

      C=(F32-F21)/X31
      B=F21+X21*C
      A=F2
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F32
      XE=1.D0-2.D0*XF/X32
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

  150 CONTINUE
C              Extrapolation for X Outside of Interval  X(NXF-1)-X(NXF)
C              --------------------------------------------------------
C                  IF(KXTRAP.EQ.0)  (No Extrapolation:  sets F(XX)=0.0)
C                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP.EQ.0) FF=0.D0
      IF(KXTRAP.EQ.1) FF=F3
      IF(KXTRAP.EQ.2) FF=F3+XF*(F3-F2)/(X3-X2)

  160 CONTINUE
      RETURN
      END SUBROUTINE SPLINE

      SUBROUTINE REPART(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT NONE

C     ------------------------------------------------------------------
C
C     REPART   Repartitions FXL (a histogram-tyoe distribution function)
C              where XLB depicts the NXB partitions that define FXL data
C              FXL is assumed to be constant between XLB(N) AND XLB(N+1)
C
C              GYL(N) is the new histogram distribution function defined
C              by the (input) NYB partitions YLB(N).  The YLB partitions
C              can differ from XLB in number and/or spacing, or in upper
C              and lower limits.  XLB and YLB coordinates are assumed to
C              be linear both increasing or decreasing in the same sense
C
C     RETURNS  GYL as a histogram-type distribution function, with NYB-1
C              values assumed to be constant between YLB(N) and YLB(N+1)
C
C       NOTE:  The column amount of a vertically distributed quantity is
C              conserved, within the Repartition Interval that is common
C              to to XLB and YLB (layer bottom edge and top edge) limits
C
C     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NXB,NYB
      REAL*8, INTENT(IN) :: FXL(NXB),XLB(NXB),YLB(NYB)
      REAL*8, INTENT(OUT) :: GYL(NYB)
      INTEGER NXF,NYG
      REAL*8 SUMG,SUMY,XA,YA,XB,YB,XAYA,PART
      INTEGER I,J

      NXF=NXB-1
      NYG=NYB-1
      SUMG=0.D0
      DO 50 I=1,NYG
      GYL(I)=0.D0
   50 CONTINUE
      SUMY=0.D0
      I=1
      XA=XLB(I)
      J=1
      YA=YLB(J)
      XB=XLB(I+1)
      IF(XB.LT.XA) GO TO 200
  100 CONTINUE
      YB=YLB(J+1)
      IF(YB.GT.XA) GO TO 110
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 160
      YA=YB
      GO TO 100
  110 CONTINUE
      XB=XLB(I+1)
      IF(XB.GT.YA) GO TO 120
      I=I+1
      IF(I.GT.NXF) GO TO 160
      XA=XB
      GO TO 110
  120 CONTINUE
      XAYA=XA
      IF(YA.GT.XA) XAYA=YA
      IF(YB.GT.XB) GO TO 130
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG
      J=J+1
      IF(J.GT.NYG) GO TO 160
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 120
  130 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 140
      XA=XB
      XB=XLB(I+1)
      GO TO 120
  140 CONTINUE
      GYL(J)=SUMG
  150 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 160
      GYL(J)=0.D0
      GO TO 150
  160 CONTINUE
      GO TO 300

  200 CONTINUE
      YB=YLB(J+1)
      IF(YB.LT.XA) GO TO 210
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 260
      YA=YB
      GO TO 200
  210 CONTINUE
      XB=XLB(I+1)
      IF(XB.LT.YA) GO TO 220
      I=I+1
      IF(I.GT.NXF) GO TO 260
      XA=XB
      GO TO 210
  220 CONTINUE
      XAYA=XA
      IF(YA.LT.XA) XAYA=YA
      IF(YB.LT.XB) GO TO 230
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG
      J=J+1
      IF(J.GT.NYG) GO TO 260
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 220
  230 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 240
      XA=XB
      XB=XLB(I+1)
      GO TO 220
  240 CONTINUE
      GYL(J)=SUMG
  250 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 260
      GYL(J)=0.D0
      GO TO 250
  260 CONTINUE

  300 CONTINUE
      RETURN
      END SUBROUTINE REPART

      SUBROUTINE RETERP(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     RETERP
C              Interpolates FXL (a histogram-type distribution function)
C              where XLB depicts the NXB partitions that define FXL data
C              FXL is assumed to be constant between XLB(N) AND XLB(N+1)
C
C              GYL(N) is the new histogram distribution function defined
C              by the (input) NYB partitions YLB(N).  The YLB partitions
C              can differ from XLB in number and/or spacing, or in upper
C              and lower limits.  XLB and YLB coordinates are assumed to
C              be linear both increasing or decreasing in the same sense
C
C     RETURNS  GYL as a histogram-type distribution function, with NYB-1
C              values assumed to be constant between YLB(N) and YLB(N+1)
C
C     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NXB,NYB
      REAL*8, INTENT(IN) :: FXL(NXB),XLB(NXB),YLB(NYB)
      REAL*8, INTENT(OUT) :: GYL(NYB)
      INTEGER NXF,NYG
      REAL*8 SUMG,SUMY,XA,YA,XB,YB,XAYA,PART
      INTEGER I,J

      NXF=NXB-1
      NYG=NYB-1
      SUMG=0.D0
      DO 50 I=1,NYG
      GYL(I)=0.D0
   50 CONTINUE
      SUMY=0.D0
      I=1
      XA=XLB(I)
      J=1
      YA=YLB(J)
      XB=XLB(I+1)
      IF(XB.LT.XA) GO TO 200
  100 CONTINUE
      YB=YLB(J+1)
      IF(YB.GT.XA) GO TO 110
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 160
      YA=YB
      GO TO 100
  110 CONTINUE
      XB=XLB(I+1)
      IF(XB.GT.YA) GO TO 120
      I=I+1
      IF(I.GT.NXF) GO TO 160
      XA=XB
      GO TO 110
  120 CONTINUE
      XAYA=XA
      IF(YA.GT.XA) XAYA=YA
      IF(YB.GT.XB) GO TO 130
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG/SUMY
      J=J+1
      IF(J.GT.NYG) GO TO 160
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 120
  130 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 140
      XA=XB
      XB=XLB(I+1)
      GO TO 120
  140 CONTINUE
      GYL(J)=SUMG/SUMY
  150 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 160
      GYL(J)=0.D0
      GO TO 150
  160 CONTINUE
      GO TO 300

  200 CONTINUE
      YB=YLB(J+1)
      IF(YB.LT.XA) GO TO 210
      GYL(J)=0.D0
      J=J+1
      IF(J.GT.NYG) GO TO 260
      YA=YB
      GO TO 200
  210 CONTINUE
      XB=XLB(I+1)
      IF(XB.LT.YA) GO TO 220
      I=I+1
      IF(I.GT.NXF) GO TO 260
      XA=XB
      GO TO 210
  220 CONTINUE
      XAYA=XA
      IF(YA.LT.XA) XAYA=YA
      IF(YB.LT.XB) GO TO 230
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG/SUMY
      J=J+1
      IF(J.GT.NYG) GO TO 260
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 220
  230 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I.GT.NXF) GO TO 240
      XA=XB
      XB=XLB(I+1)
      GO TO 220
  240 CONTINUE
      GYL(J)=SUMG/SUMY
  250 CONTINUE
      J=J+1
      IF(J.GT.NYG) GO TO 260
      GYL(J)=0.D0
      GO TO 250
  260 CONTINUE

  300 CONTINUE
      RETURN
      END SUBROUTINE RETERP

      SUBROUTINE FABINT(F,X,NX,ALIM,BLIM,ABINT)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     FABINT  PERFORMS NUMERICAL INTEGRATION (AREA UNDER CURVE) OF F(X)
C             BETWEEN THE LIMITS X=ALIM AND X=BLIM  (WITH BLIM GT ALIM)
C
C       F(X)  IS DEFINED BY CONNECTING SUCCESSIVE F(X) DATA POINTS USING
C             STRAIGHT-LINE SEGMENTS, I.E. F(X) IS PIECE-WISE CONTINUOUS
C             THE  X  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
C
C             (F(X) IS ZERO OUTSIDE THE INTERVAL BETWEEN X(1) AND X(NX))
C     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NX
      REAL*8, INTENT(IN) :: F(NX),X(NX),ALIM,BLIM
      REAL*8, INTENT(OUT) :: ABINT
      REAL*8, PARAMETER :: DELTA=1.D-07
      REAL*8 XA,XB,XX,XMIN,XMAX,XJ,XI,FI,FJ,BF,AF,DINT,X2,X1
      INTEGER JX,KX,IX

      ABINT=0.D0
      JX=1
      KX=1
      XA=X(JX)
      XB=X(NX)
      XX=XA
      IF(XB.GT.XA) GO TO 120
      XA=XB
      XB=XX
      JX=NX
      KX=-1
 120  XMIN=XA
      XMAX=XB
      IF(XMIN.GE.BLIM) RETURN
      IF(XMAX.LE.ALIM) RETURN
      IF(XMIN.LT.ALIM) XMIN=ALIM
      IF(XMAX.GT.BLIM) XMAX=BLIM
 130  JX=JX+KX
      XJ=X(JX)
      IF(XJ.LE.XMIN) GO TO 130
      IX=JX-KX
      XI=X(IX)
      IF((XJ-XI).LT.DELTA) GO TO 130
      FI=F(IX)
      FJ=F(JX)
      BF=(FJ-FI)/(XJ-XI)
      AF=FJ-BF*XJ
      X2=XMIN
 160  X1=X2
      X2=XJ
      IF(X2.GT.XMAX) X2=XMAX
      DINT=AF*(X2-X1)+BF*(X2**2-X1**2)/2.D0
      ABINT=ABINT+DINT
      IF(DABS(X2-XMAX).LT.DELTA) RETURN
      IF((XJ-X2).GT.DELTA) GO TO 160
 170  XI=XJ
      FI=FJ
      IX=JX
      JX=JX+KX
      XJ=X(JX)
      FJ=F(JX)
      IF(DABS(XJ-XI).LT.DELTA) GO TO 170
      BF=(FJ-FI)/(XJ-XI)
      AF=FJ-BF*XJ
      GO TO 160
      END SUBROUTINE FABINT

      SUBROUTINE FXGINT(F,X,NX,G,Y,NY,ALIM,BLIM,ABINT)
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     FXGINT  PERFORMS NUMERICAL INTEGRATION (AREA UNDER CURVE) OF  F*G
C             BETWEEN THE LIMITS X=ALIM AND X=BLIM  (WITH BLIM GT ALIM)
C
C       F(X)  IS DEFINED BY CONNECTING SUCCESSIVE F(X) DATA POINTS USING
C             STRAIGHT-LINE SEGMENTS, I.E. F(X) IS PIECE-WISE CONTINUOUS
C             THE  X  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
C
C       G(Y)  IS DEFINED BY CONNECTING SUCCESSIVE G(Y) DATA POINTS USING
C             STRAIGHT-LINE SEGMENTS, I.E. G(Y) IS PIECE-WISE CONTINUOUS
C             THE  Y  COORDINATE CAN BE IN ASCENDING OR DESCENDING ORDER
C
C             (X,Y ARE THE SAME LINEAR COORDINATE INDEPENDENTLY DEFINED)
C
C             (F(X) IS ZERO OUTSIDE THE INTERVAL BETWEEN X(1) AND X(NX))
C             (G(Y) IS ZERO OUTSIDE THE INTERVAL BETWEEN Y(1) AND Y(NY))
C     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NX,NY
      REAL*8, INTENT(IN) :: F(NX),X(NX),G(NY),Y(NY),ALIM,BLIM
      REAL*8, INTENT(OUT) :: ABINT
      REAL*8, PARAMETER :: DELTA=1.D-07
      REAL*8 XA,YA,XB,YB,XX,XMIN,XMAX,XJ,XI,FI,FJ,BF,AF,YI,YJ,GI,GJ,AG
     *     ,BG,DINT,X2,X1
      INTEGER JX,JY,KX,KY,IX,IY

      ABINT=0.D0
      JX=1
      JY=1
      KX=1
      KY=1
      XA=X(JX)
      YA=Y(JY)
      XB=X(NX)
      YB=Y(NY)
      XX=XA
      IF(XB.GT.XA) GO TO 100
      XA=XB
      XB=XX
      JX=NX
      KX=-1
 100  XX=YA
      IF(YB.GT.YA) GO TO 120
      YA=YB
      YB=XX
      JY=NY
      KY=-1
 120  XMIN=MAX(XA,YA)
      XMAX=MIN(XB,YB)
      IF(XMIN.GE.BLIM) RETURN
      IF(XMAX.LE.ALIM) RETURN
      IF(XMIN.LT.ALIM) XMIN=ALIM
      IF(XMAX.GT.BLIM) XMAX=BLIM
 130  JX=JX+KX
      XJ=X(JX)
      IF(XJ.LE.XMIN) GO TO 130
      IX=JX-KX
      XI=X(IX)
      IF((XJ-XI).LT.DELTA) GO TO 130
      FI=F(IX)
      FJ=F(JX)
      BF=(FJ-FI)/(XJ-XI)
      AF=FJ-BF*XJ
 140  JY=JY+KY
      YJ=Y(JY)
      IF(YJ.LE.XMIN) GO TO 140
      IY=JY-KY
      YI=Y(IY)
      IF((YJ-YI).LT.DELTA) GO TO 140
      GI=G(IY)
      GJ=G(JY)
      BG=(GJ-GI)/(YJ-YI)
      AG=GJ-BG*YJ
      X2=XMIN
 160  X1=X2
      X2=MIN(XJ,YJ)
      IF(X2.GT.XMAX) X2=XMAX
      DINT=(AF*AG)*(X2-X1)
     *    +(AF*BG+BF*AG)*(X2**2-X1**2)/2.D0
     *    +(BF*BG)*(X2**3-X1**3)/3.D0
      ABINT=ABINT+DINT
      IF(DABS(X2-XMAX).LT.DELTA) RETURN
      IF((XJ-X2).GT.DELTA) GO TO 180
 170  XI=XJ
      FI=FJ
      IX=JX
      JX=JX+KX
      XJ=X(JX)
      FJ=F(JX)
      IF(DABS(XJ-XI).LT.DELTA) GO TO 170
      BF=(FJ-FI)/(XJ-XI)
      AF=FJ-BF*XJ
 180  CONTINUE
      IF(YJ.GT.X2) GO TO 160
      YI=YJ
      GI=GJ
      IY=JY
      JY=JY+KY
      YJ=Y(JY)
      GJ=G(JY)
      IF(DABS(YJ-YI).LT.DELTA) GO TO 180
      BG=(GJ-GI)/(YJ-YI)
      AG=GJ-BG*YJ
      GO TO 160
      END SUBROUTINE FXGINT


c*********************************************************************
c  The following part computes the albedo. It has been taken out
c  from the radiation module with the intention to convert it later
c  into a separate module.
c*********************************************************************

      SUBROUTINE UPDSUR(JYEARR,JJDAYR)
!@sum UPDSUR updates variables for surface albedo once a day
!@auth A. Lacis/V. Oinas (modifications by I. Alienov/G. Schmidt)
      use SURF_ALBEDO, only : albvnh,albvnd,season,nv
      implicit none
!@var jyearr radiation year (not used for anything yet)
!@var jjdayr julian day (used for seasonality of veg albedo)
      integer, intent(in) :: jyearr, jjdayr
      integer K,KS1,KS2,KN1,KN2,L
      real*8 XJDAY,SEASN1,SEASN2,WT2,WT1
C
C                      Define Seasonal Albedo Dependence
C                      ---------------------------------
C
      XJDAY=JJDAYR
      SEASN1=-77.0D0
      DO K=1,4
        SEASN2=SEASON(K)
        IF(XJDAY.LE.SEASN2) GO TO 120
        SEASN1=SEASN2
      END DO
      K=1
      SEASN2=380.0D0
  120 CONTINUE
      WT2=(XJDAY-SEASN1)/(SEASN2-SEASN1)
      WT1=1.D0-WT2
      KS1=1+MOD(K,4)
      KS2=1+MOD(K+1,4)
      KN1=1+MOD(K+2,4)
      KN2=K

      DO K=1,NV
      DO L=1,6
C     -------------------
C     Southern Hemisphere
C     -------------------
        ALBVNH(K,L,1)=WT1*ALBVND(K,KS1,L)+WT2*ALBVND(K,KS2,L)
C     -------------------
C     Northern Hemisphere
C     -------------------
        ALBVNH(K,L,2)=WT1*ALBVND(K,KN1,L)+WT2*ALBVND(K,KN2,L)
      END DO
      END DO
      RETURN
      END SUBROUTINE UPDSUR

      SUBROUTINE GETSUR
!@sum GETSUR computes surface albedo for each grid box
!@auth A. Lacis/V. Oinas (modifications by I. Alienov/G. Schmidt)
      use SURF_ALBEDO
      use RADPAR, only:
C**** config data
     *     MLAT46,jnorth,KEEPAL,KSIALB,snoage_fac_max,KZSNOW,MADSUR,
C**** input from radiation
     *     COSZ,PLANCK,ITNEXT,ITPFT0,
C**** input from driver
     *     AGESN,ILON,JLAT,POCEAN,POICE,PEARTH,PLICE,PLAKE,zlake,
     *     TGO,TGOI,TGE,TGLI,TSL,ZOICE,FLAGS,FMP,ZSNWOI,zmp,
     *     SNOWOI,SNOWE,SNOWLI,SNOW_FRAC,WEARTH,WMAG,PVT,
C**** output
     *     BXA,PRNB,PRNX,SRBALB,SRXALB,TRGALB,BGFEMD,BGFEMT,
     *     DTRUFG,FTRUFG
      implicit none
      integer K,J,I,L,JH,IWM,JWM,ITOC,ITEA,ITOI,ITLI
      real*8 ASNAGE,FSNAGE,BOCSUM,BEASUM,BOISUM,BLISUM,AVSCUM,ANSCUM,WMJ
     *   ,WMI,FRFOAM,AV,BV,WTOC,BOCM,BOCP,TRAPOC,BOCM1,BOCP1,BOC
     *   ,DSFRAC,VGFRAC,SEAVIS,SEANIR,VTFRAC,WTEA,BEAM,BEAP,TRAPEA
     *   ,BEAM1,BEAP1,BEA,WTOI,BOIM,BOIP,TRAPOI,BOIM1,BOIP1,BOI,AMEAN
     *   ,WTLI,BLIM,BLIP,BGF,TRAPLI,BLIM1,BLIP1,BLI,KKZSNO,FDZICE,AHMZOI

C**** local variables for albedos for each surface type
      real*8 BOCVIS,BEAVIS,BOIVIS,BLIVIS,BOCNIR,BEANIR,BOINIR,BLINIR,
     +       XOCVIS,XEAVIS,XOIVIS,XLIVIS,XOCNIR,XEANIR,XOINIR,XLINIR,
     +       EXPSNE,EXPSNO,EXPSNL,BSNVIS,BSNNIR,XSNVIS,XSNNIR,
     +       BVSOIL,BNSOIL,BVVEGE,BNVEGE,XVSOIL,XNSOIL,XVVEGE,XNVEGE,
     +       BVSURF,BNSURF,XVSURF,XNSURF

C**** arrays needed if 6 band albedo is used
      real*8, dimension(6) :: BOCVN,BEAVN,BOIVN,BLIVN,BSNVN,BVNSUR,
     *                        XOCVN,XEAVN,XOIVN,XLIVN,XSNVN,XVNSUR

C**** Equivalence 2 band variables to 6 band array for easier passing
      EQUIVALENCE
     *     (BOCVN(1),BOCVIS),(BOCVN(2),BOCNIR),
     *     (BEAVN(1),BEAVIS),(BEAVN(2),BEANIR),
     *     (BOIVN(1),BOIVIS),(BOIVN(2),BOINIR),
     *     (BLIVN(1),BLIVIS),(BLIVN(2),BLINIR),
     *     (BSNVN(1),BSNVIS),(BSNVN(2),BSNNIR),
     *     (XOCVN(1),XOCVIS),(XOCVN(2),XOCNIR),
     *     (XEAVN(1),XEAVIS),(XEAVN(2),XEANIR),
     *     (XOIVN(1),XOIVIS),(XOIVN(2),XOINIR),
     *     (XLIVN(1),XLIVIS),(XLIVN(2),XLINIR),
     *     (XSNVN(1),XSNVIS),(XSNVN(2),XSNNIR),
     *     (BVNSUR(1),BVSURF),(BVNSUR(2),BNSURF),
     *     (XVNSUR(1),XVSURF),(XVNSUR(2),XNSURF)

C**** variables used for sea ice albedo calculation (6 bands, Hansen)
      real*8, dimension(6) :: almp6,alsf6
      real*8 :: patchy,snagfac

C     -----------------------------------------------------------------
C     Ocean Albedo Dependence on Zenith Angle and Wind Speed
C
      real*8 BVH2O, XVH2O, WMAG1, X
      BVH2O(WMAG1)=.0488D0+.0974D0/(5.679D0+WMAG1)+
     &     .0004D0/(.3333D0+WMAG1)
      XVH2O(WMAG1,X)=.021D0+X*X*(.0421D0+X*(.1283D0+X*(-.04D0+X*(3.117D0
     +              /(5.679D0+WMAG1)+X*.025D0/(.3333D0+WMAG1)))))
C     -----------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Select albedo computation using KSIALB
C     KSIALB= 0  Schramm oi.alb, Antarc/Greenl alb=.8 (J.Hansen)
C     KSIALB= 1  6-band original albedo - no 'fixups' (Andy Lacis)

C     For offline use: if MADSUR=1 get vegetation fractions from ij-map
      if(MADSUR.eq.1) call getveg(ilon,jlat)

C           Get Albedo, Thermal Flux, Flux Derivative for each Surf Type
C           ------------------------------------------------------------

      JH=1
      IF(JLAT.GT.JNORTH) JH=2
      KKZSNO=KZSNOW
      IF(COSZ.LT.0.001) KKZSNO=0
C
      EXPSNE=1.D0 ; EXPSNO=1.D0 ;  EXPSNL=1.D0
C
      DO K=1,NKBAND
        TRGALB(K)=0.D0
        BGFEMD(K)=0.D0
        BGFEMT(K)=0.D0
      END DO
C
      BOCSUM=0.D0
      BEASUM=0.D0
      BOISUM=0.D0
      BLISUM=0.D0
      DO K=1,4
        DTRUFG(K)=0.D0
      END DO
C
      AVSCUM=0.D0
      ANSCUM=0.D0
C
      BOCVN=0. ; BEAVN=0. ; BOIVN=0. ; BLIVN=0. ; BSNVN=0.
      XOCVN=0. ; XEAVN=0. ; XOIVN=0. ; XLIVN=0. ; XSNVN=0. ; BXA=0.
C
C                                             --------------------------
C                                             Ocean Albedo Specification
C                                             --------------------------
C
      IF(POCEAN.LT.1.D-04) GO TO 400
      X=0.5D0+(0.5D0-COSZ)*ZOCSRA
      BOCVIS=BVH2O(WMAG)  +AVSCAT+AVSCUM
      XOCVIS=XVH2O(WMAG,X)+AVSCAT+AVSCUM
      BOCNIR=BVH2O(WMAG)  +ANSCAT+ANSCUM
      XOCNIR=XVH2O(WMAG,X)+ANSCAT+ANSCUM
C
      IWM=WMAG
      IF(IWM.LT.1) IWM=1
      IF(IWM.GT.24) IWM=24
      JWM=IWM+1
      WMJ=WMAG-IWM
      WMI=1.D0-WMJ
      FRFOAM=WMI*SRFOAM(IWM)+WMJ*SRFOAM(JWM)
C
      BOCVIS=BOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
      XOCVIS=XOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
      BOCNIR=BOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
      XOCNIR=XOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
      DO L=3,6                  ! fill in higher bands
        BOCVN(L)=BOCNIR         ! 1/2 already equivalenced
        XOCVN(L)=XOCNIR
      END DO

C**** For lakes increase albedo if lakes are very shallow
C**** This is a fix to prevent lakes from overheating when they
C**** are not allowed to evaporate away. We assume that departing
C**** lakes would leave bare soil behind.
C**** Reduce effect by a half...
      if (PLAKE.gt.0 .and. zlake.lt.1.) then  ! < 1m
        DO L=1,6
          BOCVN(L)=(BOCVN(L)*(zlake-0.4d0)+
     *         (1.-zlake)*0.5*ALBVNH(1,L,1))/0.6d0
          XOCVN(L)=(XOCVN(L)*(zlake-0.4d0)+
     *         (1.-zlake)*0.5*ALBVNH(1,L,1))/0.6d0
        END DO
      end if
C
      X=1.D0/(1.D0+WMAG)
      AV=(-.0147087D0*X*X+.0292266D0*X-.0081079D0)*EOCTRA
      BV=(1.01673D0-0.0083652D0*WMAG)*EOCTRA
C
      ITOC=TGO
      WTOC=TGO-ITOC
      ITOC=ITOC-ITPFT0
      IF(ITOC.LT.124) ITOC=124
      BOCSUM=0.D0
      BOCM=0.D0
      BOCP=0.D0
C
      DO K=1,NKBAND
        TRAPOC=AV+BV*AOCEAN(K)
        BOCM1 =(PLANCK(ITOC-1)-(PLANCK(ITOC-1)-PLANCK(ITOC  ))*WTOC)
     +       *(1.D0-TRAPOC)
        BOCM  =BOCM+BOCM1
        BOCP1 =(PLANCK(ITOC+1)-(PLANCK(ITOC+1)-PLANCK(ITOC+2))*WTOC)
     +       *(1.D0-TRAPOC)
        BOCP  =BOCP+BOCP1
        BOC   =(PLANCK(ITOC  )-(PLANCK(ITOC  )-PLANCK(ITOC+1))*WTOC)
     +       *(1.D0-TRAPOC)
        BOCSUM=BOCSUM+BOC
        ITOC=ITOC+ITNEXT
C
        TRGALB(K)=TRGALB(K)+POCEAN*TRAPOC
        BGFEMD(K)=BGFEMD(K)+POCEAN*(BOCP1-BOCM1)
        BGFEMT(K)=BGFEMT(K)+POCEAN*BOC
      END DO
      DTRUFG(1)=0.5D0*(BOCP-BOCM)
C
  400 CONTINUE
      IF(PEARTH.LT.1.D-04) GO TO 500
C
C                                         ------------------------------
C                                         Land Snow Albedo Specification
C                                         ------------------------------
      ASNAGE=ALBDIF(1,JH)*EXP(-AGEXPF(1,JH)*AGESN(1))
      DO L=1,6
        FSNAGE=1.D0
        IF(L.GT.2) FSNAGE=2.0D0/L
        BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
C**** Set zenith angle dependence if required
        IF (KKZSNO.GT.0) THEN
          CALL RXSNOW(BSNVN(L),COSZ,GZSNOW(L,1,JH),XSNVN(L))
        ELSE
          XSNVN(L)=BSNVN(L)
        END IF
      END DO

C                                          -----------------------------
C                                          Soil/Veg Albedo Specification
C                                          -----------------------------
c**** In the following code when computing albedo for snow covered soil
c**** we use snow_frac(1:2) which is the snow fraction cover for
c**** bare/vegetated soil. It is computed in GHY_DRV.f in accordance
c**** with the surface topography.
c**** The final snow cover is minimum of snow_frac and the snow fraction
c**** obtained using the vegetation masking.
      DSFRAC=PVT(1)+PVT(10)
      VGFRAC=1.D0-DSFRAC
      IF(SNOWE .LE.1.D-04) THEN
        DO L=1,6
          BEAVN(L)=PVT(1)*ALBVNH(1,L,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        END DO
        BVSOIL=BEAVN(1)
        BNSOIL=BEAVN(2)
        DO K=2,NVEG
          DO L=1,6
            BEAVN(L)=BEAVN(L)+PVT(K)*ALBVNH(K,L,JH)
          END DO
        END DO
        SEAVIS=BEAVN(1)
        SEANIR=BEAVN(2)
        BVVEGE=BVSOIL
        BNVEGE=BNSOIL
        IF(VGFRAC.GT.0.001D0) THEN
          BVVEGE=(BEAVN(1)-BVSOIL*DSFRAC)/VGFRAC
          BNVEGE=(BEAVN(2)-BNSOIL*DSFRAC)/VGFRAC
        ENDIF
      ELSE
        VTFRAC=PVT(1)*MAX((1.d0-snow_frac(1)),EXP(-SNOWE/VTMASK(1)))
        EXPSNE=VTFRAC +
     &       PVT(10)*MAX((1.d0-snow_frac(1)),EXP(-SNOWE/VTMASK(10)))
        DSFRAC=EXPSNE
        DO L=1,6
          BEAVN(L)=VTFRAC*ALBVNH(1,L,JH)*(1.D0-0.5D0*WEARTH*WETSRA)
        END DO
        DO K=2,NVEG
          VTFRAC=PVT(K)*MAX((1.d0-snow_frac(2)),EXP(-SNOWE/VTMASK(K)))
          DO L=1,6
            BEAVN(L)=BEAVN(L)+VTFRAC*ALBVNH(K,L,JH)
          END DO
          EXPSNE=EXPSNE+VTFRAC
        END DO
      END IF
      DO L=1,6
        XEAVN(L)=BEAVN(L)
      END DO
      DO L=1,6
        BEAVN(L)=BEAVN(L)+BSNVN(L)*(1.D0-EXPSNE)
        XEAVN(L)=XEAVN(L)+XSNVN(L)*(1.D0-EXPSNE)
      END DO
      VGFRAC=EXPSNE-DSFRAC
      XVSOIL=BVSOIL
      XNSOIL=BNSOIL
      XVVEGE=BVVEGE
      XNVEGE=BNVEGE

      ITEA=TGE
      WTEA=TGE-ITEA
      ITEA=ITEA-ITPFT0
      BEASUM=0.D0
      BEAM=0.D0
      BEAP=0.D0

      DO K=1,NKBAND
        TRAPEA=AGSIDV(K,1)*(1.D0-EXPSNE)
     +        +AGSIDV(K,3)*DSFRAC*EDSTRA*(1.D0-WETTRA*WEARTH)
     +        +AGSIDV(K,4)*VGFRAC*EVGTRA
        BEAM1 =(PLANCK(ITEA-1)-(PLANCK(ITEA-1)-PLANCK(ITEA  ))*WTEA)
     +       *(1.D0-TRAPEA)
        BEAM  =BEAM+BEAM1
        BEAP1 =(PLANCK(ITEA+1)-(PLANCK(ITEA+1)-PLANCK(ITEA+2))*WTEA)
     +       *(1.D0-TRAPEA)
        BEAP  =BEAP+BEAP1
        BEA   =(PLANCK(ITEA  )-(PLANCK(ITEA  )-PLANCK(ITEA+1))*WTEA)
     +       *(1.D0-TRAPEA)
        BEASUM=BEASUM+BEA
        ITEA=ITEA+ITNEXT
C
        TRGALB(K)=TRGALB(K)+PEARTH*TRAPEA
        BGFEMD(K)=BGFEMD(K)+PEARTH*(BEAP1-BEAM1)
        BGFEMT(K)=BGFEMT(K)+PEARTH*BEA
      END DO
      DTRUFG(2)=0.5D0*(BEAP-BEAM)
C
  500 CONTINUE
      IF(POICE.LT.1.D-04) GO TO 600
C
C                                         ------------------------------
C                                         Ocean Ice Albedo Specification
C                                         ------------------------------
      IF(KSIALB.eq.1) THEN
C****                      original version (6-band)
        EXPSNO=EXP(-SNOWOI/DMOICE)
C**** Set snow albedo over sea ice
        ASNAGE=ALBDIF(2,JH)*EXP(-AGEXPF(2,JH)*AGESN(2))
        DO L=1,6
          FSNAGE=1.D0
          IF(L.GT.2) FSNAGE=2.0D0/L
          BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
        END DO

C**** set ice albedo
        AHMZOI=ASHZOI
        IF(JLAT.GT.JNORTH) AHMZOI=ANHZOI
        FDZICE=zoice/(zoice+AHMZOI)   ! ZOICE = ice depth (m)
        DO L=1,6
          BOIVN(L)=FDZICE*AOIALB(L)*EXPSNO+BSNVN(L)*(1.D0-EXPSNO)
C**** Set zenith angle dependence if required
          IF (KKZSNO.GT.0) THEN
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,2,JH),XOIVN(L))
          ELSE
            XOIVN(L)=BOIVN(L)
          END IF
        END DO
C**** end of original version (6-band)

      else        ! KSIALB=0
C**** Schramm/J. Hansen's sea ice albedo formulas (6 spectral bands)
C**** Bare ice:
        BOIVN(1:6) = aoimax(1:6)
        if(zoice.lt.1.)then       ! zoice: ice depth (at least Z1I=.1m)
          BOIVN(1:4)=aoimin(1:4)+(aoimax(1:4)-aoimin(1:4))*sqrt(zoice)
        endif
C**** Snow:    patchy: snow_cover_fraction (<1 if snow depth < .1m)
        patchy = 0.       !  max(0.d0 , min(1.d0, 10.d0*zsnwoi) )
        if(zsnwoi.gt.0.)then
          if(zsnwoi.ge.0.1d0)then      ! snow deeper than .1m
            patchy=1d0
          else
            patchy=zsnwoi/0.1d0
          endif
          if(flags)then         ! wet snow
            alsf6(1:6)=asnwet(1:6)
          else                  ! dry snow
            alsf6(1:6)=asndry(1:6)
          endif
C       snow aging based on Loth and Graf (1998)
C       Dry, Wet(thick), Wet(thin) snow decreases by
C       0.006,  0.015 and 0.071 per day, respectively (for mean)
C       assume decrease for each band is proportional
          if (flags) then
            if (zsnwoi.gt.0.25) then
              snagfac = 0.015d0/0.7d0 * AGESN(2)
            else
              snagfac = 0.071d0/0.7d0 * AGESN(2)
            end if
          else
            snagfac = 0.006d0/0.82d0  * AGESN(2)
          end if
C       make sure snow albedo doesn't get too low!
          snagfac=min(snoage_fac_max,snagfac)
          alsf6(1:6)=alsf6(1:6)*(1.-snagfac)
C       combine bare ice and snow albedos
          BOIVN(1:6)=BOIVN(1:6)*(1.-patchy)+alsf6(1:6)*patchy
        endif
C**** Melt ponds:
        almp6(1:6)=ampmin(1:6)    ! used if melt pond deeper than .5m
        if(zmp.gt.0. .and. zmp.lt.0.5)then
          almp6(1:6)=almp6(1:6)+(aoimax(1:6)-almp6(1:6))*(1.-2.*zmp)**2
        end if
c**** combined sea ice albedo
        BOIVN(1:6)=BOIVN(1:6)*(1.-fmp)+almp6(1:6)*fmp
C**** set zenith angle dependence
        IF (KKZSNO.GT.0) THEN     ! for all surface types
          DO L=1,6
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,2,JH),XOIVN(L))
          END DO
        ELSE
          XOIVN(1:6)=BOIVN(1:6)
        END IF
        EXPSNO=1.-patchy
      end if  !  KSIALB.ne.1:  Schramm/Hansen
C*
      ITOI=TGOI
      WTOI=TGOI-ITOI
      ITOI=ITOI-ITPFT0
      BOISUM=0.D0
      BOIM=0.D0
      BOIP=0.D0
C
      DO K=1,NKBAND
        TRAPOI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNO)
     +        +AGSIDV(K,2)*EICTRA*EXPSNO
        BOIM1 =(PLANCK(ITOI-1)-(PLANCK(ITOI-1)-PLANCK(ITOI  ))*WTOI)
     +       *(1.D0-TRAPOI)
        BOIM  =BOIM+BOIM1
        BOIP1 =(PLANCK(ITOI+1)-(PLANCK(ITOI+1)-PLANCK(ITOI+2))*WTOI)
     +       *(1.D0-TRAPOI)
        BOIP  =BOIP+BOIP1
        BOI   =(PLANCK(ITOI  )-(PLANCK(ITOI  )-PLANCK(ITOI+1))*WTOI)
     +       *(1.D0-TRAPOI)
        BOISUM=BOISUM+BOI
        ITOI=ITOI+ITNEXT
C
        TRGALB(K)=TRGALB(K)+POICE*TRAPOI
        BGFEMD(K)=BGFEMD(K)+POICE*(BOIP1-BOIM1)
        BGFEMT(K)=BGFEMT(K)+POICE*BOI
      END DO
      DTRUFG(3)=0.5D0*(BOIP-BOIM)
C
  600 CONTINUE
      IF(PLICE.LT.1.E-04) GO TO 700
C                                          -----------------------------
C                                          Land Ice Albedo Specification
C                                          -----------------------------
C**** Set snow albedo over land ice
      ASNAGE=ALBDIF(3,JH)*EXP(-AGEXPF(3,JH)*AGESN(3))
      DO L=1,6
        FSNAGE=1.D0
        IF(L.GT.2) FSNAGE=2.0D0/L
        BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
      END DO
      XSNVN(1:6)=BSNVN(1:6)

      EXPSNL=EXP(-SNOWLI/DMLICE)
      DO L=1,6
        BLIVN(L)=ALIALB(L)*EXPSNL+BSNVN(L)*(1.D0-EXPSNL)
      END DO
C**** For KSIALB != 1
C**** Specify the Albedo for Antarctica and Greenland: vis.alb = 95%
C**** and mean albedo=80%, i.e. AMEAN = .585*BLIVIS+.415*BLINIR = .80
C****
      if (ksialb.ne.1) then
      IF( JLAT.LT.NINT(MLAT46/6.) .OR.
     *   (JLAT.LT.45.AND.JLAT.GT.38.AND.ILON.LT.33.AND.ILON.GT.23)) THEN
        IF (KKZSNO.GT.0) THEN
          AMEAN=.78d0  ! compensate for zenith angle effects in the mean
        ELSE
          AMEAN=.8d0
        END IF
        BLIVIS=.95d0
        BLINIR=(AMEAN-.585d0*BLIVIS)/.415d0
        DO L=3,6  ! fill in higher bands
          BLIVN(L)=BLINIR
        END DO
        END IF
      end if

C**** zenith angle dependence if required
      DO L=1,6
        IF (KKZSNO.GT.0) THEN
          CALL RXSNOW(BLIVN(L),COSZ,GZSNOW(L,3,JH),XLIVN(L))
        ELSE
          XLIVN(L)=BLIVN(L)
        END IF
      END DO
C
      ITLI=TGLI
      WTLI=TGLI-ITLI
      ITLI=ITLI-ITPFT0
      BLISUM=0.D0
      BLIM=0.D0
      BLIP=0.D0
      BGF=0.D0
      DO K=1,NKBAND
        TRAPLI=AGSIDV(K,1)*ESNTRA*(1.-EXPSNL)
     +        +AGSIDV(K,2)*EICTRA*EXPSNL
        BLIM1 =(PLANCK(ITLI-1)-(PLANCK(ITLI-1)-PLANCK(ITLI  ))*WTLI)
     +       *(1.D0-TRAPLI)
        BLIM  =BLIM+BLIM1
        BLIP1 =(PLANCK(ITLI+1)-(PLANCK(ITLI+1)-PLANCK(ITLI+2))*WTLI)
     +       *(1.D0-TRAPLI)
        BLIP  =BLIP+BLIP1
        BLI   =(PLANCK(ITLI  )-(PLANCK(ITLI  )-PLANCK(ITLI+1))*WTLI)
     +       *(1.D0-TRAPLI)
        BLISUM=BLISUM+BLI
        ITLI=ITLI+ITNEXT
        TRGALB(K)=TRGALB(K)+PLICE*TRAPLI
        BGFEMD(K)=BGFEMD(K)+PLICE*(BLIP1-BLIM1)
        BGFEMT(K)=BGFEMT(K)+PLICE*BLI
      END DO
      DTRUFG(4)=0.5D0*(BLIP-BLIM)
  700 CONTINUE

C**** write some BXA for diagnostic output (in WRITER) (replaces equiv.)
      BXA(1)=EXPSNE ; BXA(2)=EXPSNO ; BXA(3)=EXPSNL
      BXA(4)=BSNVIS ; BXA(5)=BSNNIR ; BXA(6)=XSNVIS ; BXA(7)=XSNNIR

C**** calculate final variables always over 6-bands
      DO L=1,6
        BVNSUR(L)=POCEAN*BOCVN(L)+PEARTH*BEAVN(L)
     +           + POICE*BOIVN(L)+ PLICE*BLIVN(L)
        XVNSUR(L)=POCEAN*XOCVN(L)+PEARTH*XEAVN(L)
     +           + POICE*XOIVN(L)+ PLICE*XLIVN(L)
      END DO
      DO L=1,6
        J=7-L
        PRNB(J,1)=BOCVN(L)
        PRNB(J,2)=BEAVN(L)
        PRNB(J,3)=BOIVN(L)
        PRNB(J,4)=BLIVN(L)
        PRNX(J,1)=XOCVN(L)
        PRNX(J,2)=XEAVN(L)
        PRNX(J,3)=XOIVN(L)
        PRNX(J,4)=XLIVN(L)
      END DO
      IF(KEEPAL.NE.1) THEN
        DO J=1,6
          L=7-J
          SRBALB(J)=BVNSUR(L)
          SRXALB(J)=XVNSUR(L)
        END DO
      ENDIF
C
C                     --------------------------------------------------
C                     Define each Surface Flux Factors, Flux Derivatives
C                     --------------------------------------------------
      BGF=0.D0
      DO K=1,NKBAND
        BGFEMD(K)=BGFEMD(K)*0.5D0
        BGF=BGF+BGFEMT(K)
      END DO
C
!nu   BGM=BOCM*POCEAN+BEAM*PEARTH+BOIM*POICE+BLIM*PLICE
!nu   BGP=BOCP*POCEAN+BEAP*PEARTH+BOIP*POICE+BLIP*PLICE
!nu   TTRUFG=0.5D0*(BGP-BGM)
      FTRUFG(1)=BOCSUM/BGF
      FTRUFG(2)=BEASUM/BGF
      FTRUFG(3)=BOISUM/BGF
      FTRUFG(4)=BLISUM/BGF

      RETURN
      END SUBROUTINE GETSUR

      SUBROUTINE RXSNOW(RBSNO,XCOSZ,GGSNO,RXSNO)
!@sum RXSNOW calculate zenith angle dependence for snow/ice albedo
!@auth A. Lacis (modified by G. Schmidt)
      USE RADPAR, only : gtsalb,sgpgxg
      IMPLICIT NONE
!@var RBSNO diffuse albedo
      REAL*8, INTENT(IN) :: RBSNO
!@var XCOSZ zenith angle
      REAL*8, INTENT(IN) :: XCOSZ
!@var GGSNO Asymmetry parameter for snow
      REAL*8, INTENT(IN) :: GGSNO
!@var RXSNO direct albedo
      REAL*8, INTENT(OUT) :: RXSNO
      INTEGER NDBLS,NN
      REAL*8 XXG,XXT,GGSN,RBSN,FRTOP,TAU,TAUSN,GPFF,PR,PT,DBLS,SECZ,XANB
     *     ,XANX,TANB,TANX,RASB,RASX,BNORM,XNORM,RARB,RARX,XATB,DENOM,DB
     *     ,DX,UB,UX,DRBRAT,RBBOUT

      IF(RBSNO.LT.0.05D0) THEN
        RXSNO=RBSNO
        RETURN
      ENDIF
      XXG=0.D0
      XXT=0.D0
      GGSN=GGSNO
      IF(GGSNO.GT.0.9D0) GGSN=0.9D0
      RBSN=RBSNO
      FRTOP=1.D0
      IF(RBSNO.GT.0.5D0) THEN
        RBSN=0.5D0
        FRTOP=((1.D0-RBSNO)/0.5D0)**2
      ENDIF

      CALL GTSALB(XXG,XXT,RBBOUT,RBSN,GGSN,TAUSN,2)
      CALL SGPGXG(XCOSZ,TAUSN,GGSN,GPFF)
      PR=1.D0-GPFF
      PT=1.D0+GPFF
      DBLS=10.D0+1.44269D0*LOG(TAUSN)
      NDBLS=DBLS
      TAU=TAUSN/2**NDBLS
C     Set optically thin limit values of R,T,X using PI0 renormalization
C     ------------------------------------------------------------------
C
      SECZ=1.D0/XCOSZ
      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      TANB=PT*XANB
      XXT=(SECZ-2.D0)*TAU
      TANX=PT*SECZ
     +    *(.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))*XANX
      RASB=PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT=(SECZ+2.D0)*TAU
      RASX=PR*SECZ
     +    *(.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM=(1.D0-XANB)/(RASB+TANB)
      XNORM=(1.D0-XANX)/(RASX+TANX)
      RASB=RASB*BNORM
      RASX=RASX*XNORM
      TANB=TANB*BNORM
      TANX=TANX*XNORM
      DO NN=1,NDBLS
        RARB=RASB*RASB
        RARX=XANX*RASX
        XATB=XANB+TANB
        DENOM=1.D0-RARB
        DB=(TANB+XANB*RARB)/DENOM
        DX=(TANX+RARX*RASB)/DENOM
        UB=RASB*(XANB+DB)
        UX=RARX+RASB*DX
        RASB=RASB+XATB*UB
        RASX=RASX+XATB*UX
        TANB=XANB*TANB+XATB*DB
        TANX=XANX*TANX+XATB*DX
        XANB=XANB*XANB
        XANX=XANX*XANX
      END DO
      DRBRAT=RASX/RBSN-1.D0
      RXSNO=RBSNO*(1.D0+DRBRAT*FRTOP)
      RETURN
      END SUBROUTINE RXSNOW

      SUBROUTINE CTREND(JYEAR,IDEC,JDEC,CWTI,CWTJ)
      IMPLICIT NONE

C-------------------------------------------------------------------
C     Black Carbon interdecadal TAU interpolation is based on linear
C     TAU trend (between decadal global TAUmaps) with a superimposed
C     intra-decadal time dependence scaled to the Black Carbon Total
C     emission rate.
C
C        INPUT: JYEAR   (Julian year)
C
C                 CTREND coefficients refer to sep2003_OCI_Koch maps
C                 CTREND coefficients refer to sep2003_BCI_Koch maps
C                 --------------------------------------------------
C
C                 Map=  1850 1875 1900 1925 1950 1960 1970 1980 1990
C       OUTPUT:  IDEC=   (0)   1    2    3    4    5    6    7    8
C                JDEC=  IDEC + 1    (returned IDEC,JDEC are (1 to 8)
C
C                CWTI=   (Multiplicative Weight for BC DataMap IDEC)
C                CWTJ=   (Multiplicative Weight for BC DataMap JDEC)
C
C        NOTE:  Time dependence is linear before 1950. Industrial BC
C               is assumed 0 in 1850 so CWTI=0, and IDEC is set to 1
C-------------------------------------------------------------------

      INTEGER, INTENT(IN)  :: JYEAR
      INTEGER, INTENT(OUT) :: IDEC,JDEC
      REAL*8,  INTENT(OUT) :: CWTI,CWTJ

C    Global Annual Emissions of BC U   Emission (Mt/yr)

      real*8, parameter, dimension(5,45) :: BCE = RESHAPE( (/
C      Year    Hard_Coal    Brown_Coal      Diesel        Total
     A 50.0, 2.280581713, 0.4449132979, 0.1599090248, 2.885536671,
     B 51.0, 2.443193913, 0.4855868816, 0.1884280443, 3.117194653,
     C 52.0, 2.473641872, 0.5115299225, 0.2027695477, 3.187930107,
     D 53.0, 2.481340885, 0.5448409319, 0.2149295360, 3.241089582,
     E 54.0, 2.505670071, 0.5780177116, 0.2343477309, 3.317960978,
     F 55.0, 2.698692560, 0.6238067150, 0.2733324766, 3.595800638,
     G 56.0, 2.855226278, 0.6531309485, 0.3043369055, 3.812692404,
     H 57.0, 2.975781679, 0.6821750998, 0.3207367063, 3.978575468,
     I 58.0, 3.341105223, 0.7035279870, 0.3370627165, 4.381746292,
     J 59.0, 3.638528824, 0.7075053453, 0.3695519567, 4.715488434,
     K 60.0, 3.770926714, 0.7416650057, 0.3832504749, 4.896034241,
     L 61.0, 3.392980337, 0.7805693150, 0.4217525721, 4.595387459,
     M 62.0, 3.288835049, 0.8179932237, 0.4603823125, 4.567360401,
     N 63.0, 3.359177589, 0.8604368567, 0.5090782642, 4.728550911,
     O 64.0, 3.432664871, 0.8952696323, 0.5388473868, 4.866865158,
     P 65.0, 3.529418945, 0.8819132447, 0.5785927773, 4.989773750,
     Q 66.0, 3.577459812, 0.8817394972, 0.6323299408, 5.091631413,
     R 67.0, 3.418204546, 0.8635972142, 0.6592246890, 4.941041946,
     S 68.0, 3.452457905, 0.8943673372, 0.7338049412, 5.080585003,
     A 69.0, 3.626069546, 0.9298774004, 0.7889106274, 5.344810009,
     B 70.0, 3.264039755, 0.9229136109, 0.8880128860, 5.074741840,
     C 71.0, 3.437611580, 0.9374827743, 0.9531223178, 5.328329086,
     D 72.0, 3.473345757, 0.7836616039, 1.0180075170, 5.274850368,
     E 73.0, 3.495583296, 0.8056778908, 1.1174367670, 5.418928623,
     F 74.0, 3.506143808, 0.8251076341, 1.0828053950, 5.413989067,
     G 75.0, 3.906814098, 0.8527192473, 1.0454736950, 5.804963112,
     H 76.0, 4.005736828, 0.8900613785, 1.1400985720, 6.035901546,
     I 77.0, 4.236912251, 0.9103702307, 1.2190728190, 6.366260529,
     J 78.0, 4.459666252, 0.9303293228, 1.2408012150, 6.630728722,
     K 79.0, 4.697422504, 0.9856286645, 1.3019220830, 6.984815121,
     L 80.0, 4.796229839, 0.9959300756, 1.2336660620, 7.026207924,
     M 81.0, 4.789204121, 1.0459070210, 1.1664049630, 7.001126766,
     N 82.0, 4.872739315, 1.0975246430, 1.1601715090, 7.130136490,
     O 83.0, 4.983223438, 1.1424025300, 1.1732926370, 7.298912525,
     P 84.0, 5.265352249, 1.2178678510, 1.2251536850, 7.708741188,
     Q 85.0, 5.763637543, 1.2965050940, 1.2428865430, 8.303324699,
     R 86.0, 5.924767494, 1.3386499880, 1.2930148840, 8.556744576,
     S 87.0, 6.155550480, 1.3738890890, 1.3162037130, 8.845513344,
     A 88.0, 6.379704475, 1.3670797350, 1.3813229800, 9.127896309,
     B 89.0, 6.594299316, 1.4169263840, 1.4029121400, 9.414231300,
     C 90.0, 6.566919804, 1.4685817960, 1.4224120380, 9.458042145,
     D 91.0, 6.661097050, 1.2067918780, 1.4163945910, 9.284657478,
     E 92.0, 7.737902641, 1.3509917260, 1.4471185210, 10.53625107,
     F 93.0, 7.393332005, 1.2448183300, 1.4543261530, 10.09271908,
     G 94.0, 7.515841007, 1.2333894970, 1.4780857560, 10.22745800
     *   /), (/5,45/) )

      real*8 XDEC
      integer IBCDEC,JBCDEC,IJYEAR

      IF(JYEAR.LT.1876) THEN
      CWTJ=(JYEAR-1850)/25.D0
      IF(CWTJ.LT.0.D0) CWTJ=0.D0
      CWTI=0.D0
      IDEC=1
      JDEC=1
      GO TO 100
      ENDIF

      IF(JYEAR.LT.1950) THEN
      XDEC=(JYEAR-1850)/25.D0
      IDEC=XDEC
      JDEC=IDEC+1
      CWTJ=XDEC-IDEC
      CWTI=1.D0-CWTJ
      GO TO 100
      ENDIF

      IF(JYEAR.LT.1990) THEN
      IDEC=(JYEAR-1910)/10
      JDEC=IDEC+1
      IBCDEC=1+(IDEC-4)*10
      JBCDEC=IBCDEC+10
      IJYEAR=JYEAR-1949
      CWTJ=(BCE(5,IJYEAR)-BCE(5,IBCDEC))
     +    /(BCE(5,JBCDEC)-BCE(5,IBCDEC))
      CWTI=1.D0-CWTJ
      GO TO 100
      ENDIF

      IF(JYEAR.GT.1989) THEN
      IDEC=7
      JDEC=8
      IJYEAR=JYEAR-1949
      IF(IJYEAR.GT.45) IJYEAR=45
      CWTJ=BCE(5,IJYEAR)/BCE(5,41)
      CWTI=0.D0
      ENDIF
  100 CONTINUE

      RETURN
      END SUBROUTINE CTREND

      SUBROUTINE STREND(JYEAR,IDEC,JDEC,SWTI,SWTJ)
      IMPLICIT NONE

C-------------------------------------------------------------------
C     Anthropogenic Sulfate inter-decadal TAU interpolation is based
C     on a linear TAU trend (between decadal global TAU-maps) with a
C     superimposed intradecadal time dependence scaled in proportion
C     to the Anthropogenic Sulfate global emission rate.
C
C        INPUT: JYEAR   (Julian year)
C
C                 CTREND coefficients refer to sep2003_SUI_Koch maps
C                 --------------------------------------------------
C
C                 Map=  1850 1875 1900 1925 1950 1960 1970 1980 1990
C       OUTPUT:  IDEC=   (0)   1    2    3    4    5    6    7    8
C                JDEC=  IDEC + 1    (returned IDEC,JDEC are (1 to 8)
C
C                SWTI=  (Multiplicative Weight for SUI DataMap IDEC)
C                SWTJ=  (Multiplicative Weight for SUI DataMap JDEC)
C
C        NOTE:  Time dependence linear before 1950.   Industrial SUI
C               is assumed 0 in 1850 so SWTI=0, and IDEC is set to 1
C-------------------------------------------------------------------

      INTEGER, INTENT(IN)  :: JYEAR
      INTEGER, INTENT(OUT) :: IDEC,JDEC
      REAL*8,  INTENT(OUT) :: SWTI,SWTJ

C     Global Emission of Sulfate

C     Emission (Mt/yr)
C               year      Anthropogenic_Sulfate Natural_Sulfate
      real*8, parameter, dimension(3,41) :: SUE = reshape( (/
     A          1950.0,     30.46669769,           14.4,
     B          1951.0,     32.38347244,           14.4,
     C          1952.0,     32.18632889,           14.4,
     D          1953.0,     32.83379745,           14.4,
     E          1954.0,     32.79270935,           14.4,
     F          1955.0,     35.79611969,           14.4,
     G          1956.0,     39.93603897,           14.4,
     H          1957.0,     38.68806839,           14.4,
     I          1958.0,     39.35904312,           14.4,
     J          1959.0,     41.06065369,           14.4,
     K          1960.0,     42.67050934,           14.4,
     L          1961.0,     41.32410431,           14.4,
     M          1962.0,     41.80470276,           14.4,
     N          1963.0,     43.26312637,           14.4,
     O          1964.0,     44.68368530,           14.4,
     P          1965.0,     45.81701660,           14.4,
     Q          1966.0,     46.61584091,           14.4,
     R          1967.0,     46.42276001,           14.4,
     S          1968.0,     47.77438354,           14.4,
     A          1969.0,     49.30817032,           14.4,
     B          1970.0,     52.81050873,           14.4,
     C          1971.0,     52.95043945,           14.4,
     D          1972.0,     54.10167694,           14.4,
     E          1973.0,     55.93037415,           14.4,
     F          1974.0,     57.31056213,           14.4,
     G          1975.0,     58.52788162,           14.4,
     H          1976.0,     59.71361542,           14.4,
     I          1977.0,     62.59599304,           14.4,
     J          1978.0,     61.98198318,           14.4,
     K          1979.0,     64.71042633,           14.4,
     L          1980.0,     65.28986359,           14.4,
     M          1981.0,     63.23768234,           14.4,
     N          1982.0,     62.88000488,           14.4,
     O          1983.0,     61.45023346,           14.4,
     P          1984.0,     63.85008621,           14.4,
     Q          1985.0,     66.47412872,           14.4,
     R          1986.0,     68.00902557,           14.4,
     S          1987.0,     69.87956238,           14.4,
     A          1988.0,     70.52937317,           14.4,
     B          1989.0,     72.06355286,           14.4,
     C          1990.0,     71.29174805,           14.4
     *              /), (/3,41/) )

      real*8 xdec
      integer ISUDEC,JSUDEC,IJYEAR

      IF(JYEAR.LT.1876) THEN
      SWTJ=(JYEAR-1850)/25.D0
      IF(SWTJ.LT.0.D0) SWTJ=0.D0
      SWTI=0.D0
      IDEC=1
      JDEC=1
      GO TO 100
      ENDIF

      IF(JYEAR.LT.1950) THEN
      XDEC=(JYEAR-1850)/25.D0
      IDEC=XDEC
      JDEC=IDEC+1
      SWTJ=XDEC-IDEC
      SWTI=1.D0-SWTJ
      GO TO 100
      ENDIF

      IF(JYEAR.LT.1990) THEN
      IDEC=(JYEAR-1910)/10
      JDEC=IDEC+1
      ISUDEC=1+(IDEC-4)*10
      JSUDEC=ISUDEC+10
      IJYEAR=JYEAR-1949
      SWTJ=(SUE(2,IJYEAR)-SUE(2,ISUDEC))
     +    /(SUE(2,JSUDEC)-SUE(2,ISUDEC))
      SWTI=1.D0-SWTJ
      GO TO 100
      ENDIF

      IF(JYEAR.GT.1989) THEN
      IDEC=7
      JDEC=8
      IJYEAR=JYEAR-1949
      IF(IJYEAR.GT.41) IJYEAR=41
      SWTJ=SUE(2,IJYEAR)/SUE(2,41)
      SWTI=0.D0
      ENDIF
  100 CONTINUE

      RETURN
      END SUBROUTINE STREND

      SUBROUTINE SPLINV(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT NONE

      integer, intent(in)  :: NXF,              KXTRAP
      real*8 , intent(in)  :: X(NXF),F(NXF),FF, CUSPWM,CUSPWE
      real*8 , intent(out) :: XX

C---------------------------------------------------------------------
C    Inverse spline:
C    SPLINV locates FF between points (F2,X2)(F3,X3) on 4-point spread
C    and returns 4-point Cubic Spline value of XX such that FF = F(XX)
C
C    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
C    (X-Coordinate may be specified in increasing or decreasing order)
C
C---------------------------------------------------------------------
C
C    Custom Control Parameters:  CUSPWM,CUSPWE
C------------------------------
C
C    In cases where data points are unevenly spaced and/or data points
C    exhibit abrupt changes in value, Spline Interpolation may produce
C    undesirable bulging of interpolated values. In more extreme cases
C    Linear Interpolation may be less problematic to use.
C
C    Interpolation can be weighted between: Cubic Spline and Linear by
C    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
C
C    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
C    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
C
C    For example, with:
C
C    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
C    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
C
C---------------------------------------------------------------------
C
C     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
C
C               KXTRAP = 0    No Extrapolation   (i.e., sets XX = 0.0)
C                        1    Fixed Extrapolation (sets XX=edge value)
C                        2    Linear Extrapolation using 2 edge points
C
C---------------------------------------------------------------------
C
C
C    NOTE:  F(X) is assumed to be monotonic between F(1) and F(NXF)
C
C------------------------------------------------------------------

      REAL*8 x1,x2,x3,x4, x21,x32,x43,x31,x42, BETW,FFCUSP,FFLINR,CUSPWT
      REAL*8 f1,f2,f3,f4, f21,f32,f43,f3221,f4332, a,b,c,d, xf,xe,xexm
      REAL*8 DX,gg,xg,xy,deltx,slopec,slopel,slopes
      integer k,kk

      BETW=(F(2)-FF)*(F(NXF)-F(1))
      IF(BETW.GT.0.D0) GO TO 120
      BETW=(FF-F(NXF-1))*(F(NXF)-F(1))
      IF(BETW.GT.0.D0) GO TO 140

      DO 100 K=3,NXF-1
      BETW=(FF-F(K-1))*(F(K)-FF)
      DX=(FF-F(K-1))/(F(K)-F(K-1))
      XX=X(K-1)+DX*(X(K)-X(K-1))
      IF(BETW.GE.0.D0) GO TO 110
  100 CONTINUE

  110 CONTINUE
      DO 115 KK=1,5
      X1=X(K-2)
      X2=X(K-1)
      X3=X(K)
      X4=X(K+1)
      F1=F(K-2)
      F2=F(K-1)
      F3=F(K)
      F4=F(K+1)
      X21=X2-X1
      X31=X3-X1
      X32=X3-X2
      X43=X4-X3
      X42=X4-X2
      F21=(F2-F1)/(X21*X21)
      F32=(F3-F2)/(X32*X32)
      F43=(F4-F3)/(X43*X43)
      F3221=(F32+F21)/X31*X21
      F4332=(F43+F32)/X42*X43
      A=F2
      B=X32*F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)/X32
      XF=XX-X2

C                             FFCUSP= Cubic Spline Interpolation Result
C                             -----------------------------------------

      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=(X3+X2-XX-XX)/X32
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE

C                                   FFLINR= Linear Interpolation Result
C                                   -----------------------------------
      FFLINR=A+XF*F32*X32
      GG=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      SLOPEC=B+2.D0*C*XF+3.D0*D*XF**2
      SLOPEL=F32*X32
      SLOPES=SLOPEC*CUSPWT+SLOPEL*(1.D0-CUSPWT)
      XG=XF
      XY=XX
      DELTX=(GG-FF)/SLOPES
      XX=XF-(GG-FF)/SLOPES+X2
  115 CONTINUE
      GO TO 160

C                Edge Point Interval Interpolation and/or Extrapolation
C                ------------------------------------------------------
  120 CONTINUE
      BETW=(F(1)-FF)*(F(NXF)-F(1))
      IF(BETW.GT.0.D0) GO TO 130

C                          F(1),F(2)  Edge Point Interval Interpolation
C                          --------------------------------------------
      DO 125 KK=2,6
      X1=X(1)
      X2=X(2)
      X3=X(3)
      F1=F(1)
      F2=F(2)
      F3=F(3)
      XX=X1+(FF-F(1))/(F(2)-F(1))*(X2-X1)
      XF=XX-X1
      X21=X2-X1
      F21=(F2-F1)/X21
      X32=X3-X2
      X31=X3-X1
      C=((F3-F2)/X32-F21)/X31
      B=F21-X21*C
      A=F1
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F21
      XE=1.D0-2.D0*XF/X21
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      GG=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      SLOPEC=B+2.D0*C*XF
      SLOPEL=F21
      SLOPES=SLOPEC*CUSPWT+SLOPEL*(1.D0-CUSPWT)
      XG=XF
      DELTX=(GG-FF)/SLOPES
      XX=XF-(GG-FF)/SLOPES+X1
  125 CONTINUE
      GO TO 160

  130 CONTINUE
C                  Extrapolation for FF Outside of Interval F(1) - F(2)
C                  ----------------------------------------------------
C                  IF(KXTRAP.EQ.0)  (No Extrapolation:   sets XX = 0.0)
C                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP.EQ.0) XX=0.D0
      IF(KXTRAP.EQ.1) XX=X(1)
      IF(KXTRAP.EQ.2) XX=X(1)-(F(1)-FF)/(F(2)-F(1))*(X(2)-X(1))
      GO TO 160

  140 CONTINUE
      BETW=(FF-F(NXF))*(F(NXF)-F(1))
      IF(BETW.GT.0.D0) GO TO 150

C                    F(NXF-1),F(NXF)  Edge Point Interval Interpolation
C                    --------------------------------------------------
      DO 145 KK=3,7
      X1=X(NXF-2)
      X2=X(NXF-1)
      X3=X(NXF)
      F1=F(NXF-2)
      F2=F(NXF-1)
      F3=F(NXF)
      XX=X2+(FF-F2)/(F3-F2)*(X3-X2)
      XF=XX-X2
      X32=X3-X2
      F32=(F3-F2)/X32
      X21=X2-X1
      X31=X3-X1
      F21=(F2-F1)/X21

C                    3-Point Quadratic Interpolation for Edge Intervals
C                    --------------------------------------------------
C
C      (Edge Option)     ----------------------------------------------
C                        For Linear Interpolation within Edge Intervals
C                        between F(1),F(2), and between F(NXF-1),F(NXF)
C                        set the value of coefficient C below, to C=0.0
C                        ----------------------------------------------

      C=(F32-F21)/X31
      B=F21+X21*C
      A=F2
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F32
      XE=1.D0-2.D0*XF/X32
      IF(XE.LT.0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      GG=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      SLOPEC=B+2.D0*C*XF
      SLOPEL=F21
      SLOPES=SLOPEC*CUSPWT+SLOPEL*(1.D0-CUSPWT)
      XG=XF
      DELTX=(GG-FF)/SLOPES
      XX=XF-(GG-FF)/SLOPES+X2
  145 CONTINUE
      GO TO 160

  150 CONTINUE
C              Extrapolation for F Outside of Interval  F(NXF-1)-F(NXF)
C              --------------------------------------------------------
C                  IF(KXTRAP.EQ.0)  (No Extrapolation:   sets XX = 0.0)
C                  IF(KXTRAP.EQ.1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP.EQ.2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP.EQ.0) XX=0.D0
      IF(KXTRAP.EQ.1) XX=X(NXF)
      IF(KXTRAP.EQ.2) XX=X(NXF)
     +                  -(F(NXF)-FF)/(F(NXF-1)-F(NXF))*(X(NXF-1)-X(NXF))

  160 CONTINUE
      RETURN
      END SUBROUTINE SPLINV

      SUBROUTINE SPLN44(Q,NI,NJ,IR,DR,JN,DN,QQ)
      IMPLICIT NONE

      integer, intent(in)  ::   NI,NJ,  IR, JN
      real*8,  intent(in)  :: Q(NI,NJ), DR, DN
      real*8,  intent(out) :: QQ

!nu   real*8,save :: CUSPWM=1., CUSPWE=1. ,CUSPWT,fflinr
      real*8 QK(4),  ffcusp,a,b,c,d,xf,xe,xexm
      REAL*8 f1,f2,f3,f4, f21,f32,f43,f3221,f4332
      integer k,kr,irm,irp

      K=0
      IRM=IR-1
      IRP=IR+2
      DO 110 KR=IRM,IRP
      K=K+1
      F1=Q(KR,JN-1)
      F2=Q(KR,JN)
      F3=Q(KR,JN+1)
      F4=Q(KR,JN+2)
      F21=F2-F1
      F32=F3-F2
      F43=F4-F3
      F3221=0.5D0*(F32+F21)
      F4332=0.5D0*(F43+F32)
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=F3221+F4332-F32-F32
      XF=DN
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
!=1   CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
!nu   FFLINR=A+XF*F32
      QK(K)=FFCUSP ! *CUSPWT+FFLINR*(1.D0-CUSPWT)
  110 CONTINUE
      F1=QK(1)
      F2=QK(2)
      F3=QK(3)
      F4=QK(4)
      F21=F2-F1
      F32=F3-F2
      F43=F4-F3
      F3221=0.5D0*(F32+F21)
      F4332=0.5D0*(F43+F32)
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=F3221+F4332-F32-F32
      XF=DR
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
!=1   CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
!nu   FFLINR=A+XF*F32
      QQ=FFCUSP ! *CUSPWT+FFLINR*(1.D0-CUSPWT)
      RETURN
      END SUBROUTINE SPLN44

      SUBROUTINE SPLNI4(Q,NI,NJ,IR,JN,DN,QQ)
      IMPLICIT NONE

      integer, intent(in)  ::   NI,NJ,  IR, JN
      real*8,  intent(in)  :: Q(NI,NJ),     DN
      real*8,  intent(out) :: QQ

!nu   real*8,save :: CUSPWM=1., CUSPWE=1. ,CUSPWT,fflinr
      real*8 ffcusp,a,b,c,d,xf,xe,xexm
      REAL*8 f1,f2,f3,f4, f21,f32,f43,f3221,f4332

      F1=Q(IR,JN-1)
      F2=Q(IR,JN)
      F3=Q(IR,JN+1)
      F4=Q(IR,JN+2)
      F21=F2-F1
      F32=F3-F2
      F43=F4-F3
      F3221=0.5D0*(F32+F21)
      F4332=0.5D0*(F43+F32)
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=F3221+F4332-F32-F32
      XF=DN
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=1.D0-XF-XF
      IF(XE.LT.0.0) XE=-XE
      XEXM=XE**2
!=1   CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
!nu   FFLINR=A+XF*F32
      QQ=FFCUSP ! *CUSPWT+FFLINR*(1.D0-CUSPWT)
      RETURN
      END SUBROUTINE SPLNI4


      SUBROUTINE SETREL(REFF0,NAER,KDREAD           !  input parameters
     A                 ,SRUQEX,SRUQSC,SRUQCB        !  SW input ( 6,110)
     B                 ,TRUQEX,TRUQSC,TRUQCB        !  LW input (33,110)
     C                 ,REFU22,Q55U22               !  RQ input    (110)
                                                    !     SETREL  output
     D                 ,SRHQEX,SRHQSC,SRHQCB,TRHQAB !  3(6,190),(33,190)
     E                 ,RHDATA)                     !  RH info   (190,9)

      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT REAL*8(A-H,O-Z)

      real*8 REFF0,SRUQEX( 6,110),SRUQSC( 6,110),SRUQCB( 6,110)
     A            ,TRUQEX(33,110),TRUQSC(33,110),TRUQCB(33,110)
     B            ,REFU22(110),Q55U22(110)
      real*8 SRHQEX(6,190),SRHQSC(6,190),SRHQCB( 6,190)
     A            ,TRHQAB(33,190),RHDATA(190,9)

C     ------------------------------------------------------------------
C     REFF0  = Effective radius for dry aerosol seed size (in microns)
C     NAER   = Aerosol composition index
C     KDREAD = IO READ unit number for Q(m,r),g(m,r) data used by SETREL
C     ------------------------------------------------------------------
C     Aerosol index = NAER    Composition         Input data order = NNA
C                      1      SO4  Sulfate                            1
C                      2      SEA  Sea Salt                           2
C                      3      NO3  Nitrate                            3
C                                  Pure Water                         4
C                      4      ORG  Organinc                           5
C     ------------------------------------------------------------------

      character*40, save :: dtfile='oct2003.relhum.nr.Q633G633.table'
      logical, parameter :: qbinary=.false.  ; logical qexist

C     Local variables

      real*8    RHRHRH(190),RHTAUF(190),RHREFF(190)
     A         ,RHWGM2(190),RHDGM2(190),RHTGM2(190)
     B         ,RHXMFX(190),RHDENS(190),RHQ550(190),RHINFO(190,9)
      EQUIVALENCE (RHINFO(1,1),RHRHRH(1)),(RHINFO(1,2),RHTAUF(1))
      EQUIVALENCE (RHINFO(1,3),RHREFF(1)),(RHINFO(1,4),RHWGM2(1))
      EQUIVALENCE (RHINFO(1,5),RHDGM2(1)),(RHINFO(1,6),RHTGM2(1))
      EQUIVALENCE (RHINFO(1,7),RHXMFX(1)),(RHINFO(1,8),RHDENS(1))
      EQUIVALENCE (RHINFO(1,9),RHQ550(1))

      real*8    R633NR(890),XNR(31),Q633NR(890,31),G633NR(890,31)
      real*8    Q880M1(890),G880M1(890),Q880M0(890),G880M0(890)
      real*8    Q880N1(890),Q880N0(890),R550NR(890)
      real*8    RR0RHX(190),XNRRHX(190),QRH633(190),GRH633(190)
      real*8    XXMF(190),XXRH(190),XRR0(190)
      real*8    XXNR(190),XROA(190),DNRX(190)

      real*8    QXAERN(33),QSAERN(33),QGAERN(33)
     A         ,SR1QEX( 6),SR1QSC( 6),SR1QCB( 6)
     B         ,SR2QEX( 6),SR2QSC( 6),SR2QCB( 6)
     C         ,SR3QEX( 6),SR3QSC( 6),SR3QCB( 6)
     D         ,SR4QEX( 6),SR4QSC( 6),SR4QCB( 6)
     E         ,TR1QEX(33),TR1QSC(33),TR1QCB(33)
     F         ,TR2QEX(33),TR2QSC(33),TR2QCB(33)
     G         ,TR3QEX(33),TR3QSC(33),TR3QCB(33)
     H         ,TR4QEX(33),TR4QSC(33),TR4QCB(33)
     I         ,TRHQEX(33),TRHQSC(33),TRHQCB(33)

      integer, parameter, dimension(4) :: NRHCRY=(/38,47,28,38/)

!nu   CHARACTER*8 AERTYP(4)
!nu   DATA AERTYP/'Sulfate ','SeaSalt ','Nitrate ','Organic '/

C     ------------------------------------------------------------------
C     Hygroscopic aerosols (Sulfate,SeaSalt,Nitrate) physical properties
C     formulas from Tang and Munkelwitz (1994, 1996) in JGR 99, JGR 101.
C
C     AW=water activity RO=density  BX=growth factor RX=refractive index
C     SO4 = ammonium sulfate;   SEA = sea salt;   NO3 = ammonium nitrate
C     ------------------------------------------------------------------

C     functions

      real*8 AWSO4,ROSO4,BXSO4,RXSO4,DRWSO4,DRDSO4
      real*8 AWSEA,ROSEA,BXSEA,RXSEA,DRWSEA,DRDSEA,RRSEA,VVSEA,GXSEA
      real*8 AWNO3,RONO3,BXNO3,R1NO3,R2NO3, DRXNO3

      AWSO4(X)=1.D0-0.2715*X+0.3113*X**2-2.336*X**3+1.412*X**4    ! TM94
      ROSO4(X)=0.9971D0+5.92D-01*X-5.036D-02*X**2+1.024D-02*X**3  ! TM94
      BXSO4(X)=(1.D0/X*1.760D0/ROSO4(X))**(1.D0/3.D0)             ! TM96
      RXSO4(X)=1.3330+0.16730*X-0.0395*X**2                       ! TM91
      DRWSO4(RH)=1.002146-0.00149*RH+0.001*RH/(1.0+0.911*RH**10)
      DRDSO4(RH)=1.002503     ! ratio of wet & dry nr(0.550) / nr(0.633)

      AWSEA(X)=1.0D0-0.6366*X+0.8624*X**2-11.58*X**3+15.18*X**4   ! TM96
      ROSEA(X)=0.9971+0.741*X-0.3741*X**2+2.252*X**3-2.060*X**4   ! TM96
      BXSEA(X)=(1.D0/X*2.165D0/ROSEA(X))**(1.D0/3.D0)
      RRSEA(X)=3.70958+(8.95-3.70958)/(1.D0+(1.0-X)/X*58.448/18.0)
      VVSEA(X)=(18.0+(58.448-18.0)/(1.0+(1.0-X)/X*58.448/18.0))/ROSEA(X)
      GXSEA(X)=SQRT((2.D0*RRSEA(X)+VVSEA(X))/(VVSEA(X)-RRSEA(X))) ! TM96
      RXSEA(X)=1.333+(GXSEA(X)-1.333)*(1.490-1.333)/(1.544-1.333)
      DRWSEA(RH)=1.00212-0.001625*RH+0.00131*RH/(1.0+0.928*RH**3)
      DRDSEA(RH)=1.003007     ! ratio of wet & dry nr(0.550) / nr(0.633)

      AWNO3(X)=1.D0-3.65D-01*X-9.155D-02*X**2-2.826D-01*X**3      ! TM96
      RONO3(X)=0.9971D0+4.05D-01*X+9.0D-02*X**2                   ! TM96
      BXNO3(X)=(1.D0/X*1.725D0/RONO3(X))**(1.D0/3.D0)             ! TM96
      R1NO3(X)=1.3330+0.119D0*X          !  (X<0.205)              TWM81
      R2NO3(X)=1.3285+0.145D0*X          !  (X>0.205)              TWM81
      DRXNO3(RH)=1.001179     ! ratio of wet & dry nr(0.550) / nr(0.633)


C     ------------------------------------------------------------------
C     Q,G Mie data (879x31) at 0.633 microns, use 31 points to cover the
C     refractive index from 1.30 to 1.60 with equal spacing of 0.01
C
C     Q,G data effective radius spans the range from 0.0 to 20.4 microns
C     in (3) segments of equally spaced data for optimized 4-point Cubic
C     Spline interpolation.  The equally spaced segments are as follows:
C
C     Index:    1 - 303   304 - 603   604 -  879   881 - 885   886 - 890
C     Reff:   0.00-3.02   3.04-9.02   9.04-20.04   2.98-3.04   8.96-9.08
C     Delta:     0.01        0.02        0.04         0.02        0.04
C
C     The last two intervals are constructed to accommodate transitions
C     between the (3) segments using 4-point Cubic Spline interpolation
C     ------------------------------------------------------------------

      inquire (file=dtfile,exist=qexist)
      if(.not.qexist) dtfile='RH_QG_Mie  ' ! generic name used by GCM
      inquire (file=dtfile,exist=qexist)
      if(.not.qexist) call stop_model('setrel: no RH_QG files',255)
      call openunit(dtfile,kdread,qbinary) ! formatted:qbinary=.false.

      READ (KDREAD,7000) (XNR(J),J=1,31)
 7000 FORMAT(12X,F5.3,30F8.3)
      DO 101 I=1,880
      READ (KDREAD,7001) R633NR(I),(Q633NR(I,J),J=1,31)
 7001 FORMAT(3X,F6.2,31F8.5)
  101 CONTINUE
      READ (KDREAD,7000) (XNR(J),J=1,31)
      DO 102 I=1,880
      READ (KDREAD,7001) R633NR(I),(G633NR(I,J),J=1,31)
  102 CONTINUE
      call CLOSEunit (KDREAD)

      J=880
      DO 104 K=299,305
      IF(K.EQ.300) GO TO 104
      IF(K.EQ.302) GO TO 104
      J=J+1
      R633NR(J)=R633NR(K)
      DO 103 I=1,31
      Q633NR(J,I)=Q633NR(K,I)
      G633NR(J,I)=G633NR(K,I)
  103 CONTINUE
  104 CONTINUE
      DO 106 K=600,606
      IF(K.EQ.601) GO TO 106
      IF(K.EQ.603) GO TO 106
      J=J+1
      R633NR(J)=R633NR(K)
      DO 105 I=1,31
      Q633NR(J,I)=Q633NR(K,I)
      G633NR(J,I)=G633NR(K,I)
  105 CONTINUE
  106 CONTINUE

C         Set dry mass fraction XXMF and relative humidity RHRHRH scales
C         --------------------------------------------------------------
      DO 110 I=1,190
      XXMF(I)=1.D0-(I-1)/100.D0
      IF(I.GT.91)  XXMF(I)=0.10D0-(I-91)/1000.D0
      RHRHRH(I)=(I-1)/100.D0
      IF(I.GT.91) RHRHRH(I)=0.90D0+(I-91)/1000.D0
  110 CONTINUE

C         Define RH (=AW), RO, BX, RX as functions of X for NAER aerosol
C         --------------------------------------------------------------
      NRHN1=NRHCRY(NAER)+1
      DO 111 I=1,190
      RHI=RHRHRH(I)
      RR0RHX(I)=1.D0
      RHXMFX(I)=1.D0
      IF(NAER.EQ.1) THEN       !    Dry Sulfate refrac index and density
      XNRRHX(I)=1.526
      RHDENS(I)=1.760
      IF(I.LT.NRHN1) DNRX(I)=DRDSO4(RHI)
      IF(I.GE.NRHN1) DNRX(I)=DRWSO4(RHI)
      ENDIF
      IF(NAER.EQ.2) THEN       !    Dry SeaSalt refrac index and density
      XNRRHX(I)=1.490
      RHDENS(I)=2.165
      IF(I.LT.NRHN1) DNRX(I)=DRDSEA(RHI)
      IF(I.GE.NRHN1) DNRX(I)=DRWSEA(RHI)
      ENDIF
      IF(NAER.EQ.3) THEN       !    Dry Nitrate refrac index and density
      XNRRHX(I)=1.554
      RHDENS(I)=1.725
      DNRX(I)=DRXNO3(RHRHRH(I))
      ENDIF
      IF(NAER.EQ.4) THEN       !    Dry Organic refrac index and density
      XNRRHX(I)=1.526          !                  (representative value)
      RHDENS(I)=1.5            !                  (representative value)
      IF(I.LT.NRHN1) DNRX(I)=DRDSO4(RHI)
      IF(I.GE.NRHN1) DNRX(I)=DRWSO4(RHI)
      ENDIF
  111 CONTINUE
      DO 112 I=1,190
      X=XXMF(I)
      IF(NAER.EQ.1) THEN          !       RH dependent Sulfate X,R,NR,RO
      XXRH(I)=AWSO4(X)
      XRR0(I)=BXSO4(X)
      XXNR(I)=RXSO4(X)
      XROA(I)=ROSO4(X)
      ENDIF
      IF(NAER.EQ.2) THEN          !       RH dependent SeaSalt X,R,NR,RO
      XXRH(I)=AWSEA(X)
      XRR0(I)=BXSEA(X)
      XXNR(I)=RXSEA(X)
      XROA(I)=ROSEA(X)
      IF(I.LT.NRHN1) XXRH(I)=I/1000.0
      ENDIF
      IF(NAER.EQ.3) THEN          !       RH dependent Nitrate X,R,NR,RO
      XXRH(I)=AWNO3(X)
      XRR0(I)=BXNO3(X)
      XXNR(I)=R1NO3(X)
      XROA(I)=RONO3(X)
      IF(X.GT.0.205) XXNR(I)=R2NO3(X)
      ENDIF
      IF(NAER.EQ.4) THEN          !       RH dependent Organic X,R,NR,RO
      XXRH(I)=AWSO4(X)**0.50      !       (to yield growth factor G=1.1)
      XRR0(I)=BXSO4(X)**0.30      !       (at RH=84 Virkkula et al 1999)
      XXNR(I)=RXSO4(X)
      XROA(I)=ROSO4(X)*(1.0-X*(1.760-1.650)/1.760)
      ENDIF
  112 CONTINUE

C            Invert X, RO, BX, RX functions of (X) to be functions of RH
C            -----------------------------------------------------------
      DO 113 I=NRHN1,190
      RHI=RHRHRH(I)
      CALL SPLINE(XXRH,XRR0,190,RHI,RR0RHX(I),1.D0,1.D0,1)
      CALL SPLINE(XXRH,XXNR,190,RHI,XNRRHX(I),1.D0,1.D0,1)
      CALL SPLINE(XXRH,XXMF,190,RHI,RHXMFX(I),1.D0,1.D0,1)
      CALL SPLINE(XXRH,XROA,190,RHI,RHDENS(I),1.D0,1.D0,1)
      RHDENS(I)=MAX(RHDENS(I),1.D0)
  113 CONTINUE

C     ------------------------------------------------------------------
C     Find Qdry(r),gdry(r) from Q(m,r),g(m,r) maps for each aerosol type
C     Find Qwet(r),gwet(r) from Q(m,r),g(m,r) maps for each aerosol type
C     also locate MAXDRY,MAXWET pts where Qdry(r),Qwet(r) are at maximum
C          (M1 refers to mass fraction X of 1.0, i.e., "dry" aerosol)
C          (M0 refers to mass fraction X of 0.0, i.e., "wet" aerosol)
C     ------------------------------------------------------------------
      MAXDRY=1
      MAXWET=1
      QQDMAX=0.D0
      QQWMAX=0.D0
      XDRY=XNRRHX(1)
C     IF(MCRYON.EQ.1) XDRY=XNRRHX(NRHN1) ! If "dry" = RHC reference line
      SDRY=XDRY*100.D0-129
      JDRY=SDRY
      DDRY=SDRY-JDRY
      XWET=1.3330D0                     !  Pure water Nr = "wet" aerosol
      SWET=XWET*100.D0-129
      JWET=SWET
      DWET=SWET-JWET
      DO 114 I=1,880
      CALL SPLNI4(Q633NR,890,31,I,JDRY,DDRY,Q880M1(I))
      CALL SPLNI4(G633NR,890,31,I,JDRY,DDRY,G880M1(I))
      CALL SPLNI4(Q633NR,890,31,I,JWET,DWET,Q880M0(I))
      CALL SPLNI4(G633NR,890,31,I,JWET,DWET,G880M0(I))
      IF(Q880M1(I).GT.QQDMAX) THEN
      QQDMAX=Q880M1(I)
      MAXDRY=I
      ENDIF
      IF(Q880M0(I).GT.QQWMAX) THEN
      QQWMAX=Q880M0(I)
      MAXWET=I
      ENDIF
  114 CONTINUE
      RQDMAX=R633NR(MAXDRY)
      RQWMAX=R633NR(MAXWET)

C     Define:  Qdry(r) and Qwet(r) at the reference wavelength of 550 nm
C              using refractive index off-set and size parameter scaling
C     ------------------------------------------------------------------
      XDRY=XNRRHX(1)*DNRX(1)             !      Dry aerosol Nr at 550 nm
C     IF(MCRYON.EQ.1) XDRY=XNRRHX(NRHN1) ! If "dry" = RHC reference line
      SDRY=XDRY*100.D0-129
      JDRY=SDRY
      DDRY=SDRY-JDRY
      XWET=1.3330D0*1.001179     !       Pure water aerosol Nr at 550 nm
      SWET=XWET*100.D0-129
      JWET=SWET
      DWET=SWET-JWET
      DO 115 I=1,880
      CALL SPLNI4(Q633NR,890,31,I,JDRY,DDRY,Q880N1(I))
      CALL SPLNI4(Q633NR,890,31,I,JWET,DWET,Q880N0(I))
      R550NR(I)=R633NR(I)*(0.550/0.633) !  Size shift refers Q to 550 nm
  115 CONTINUE
      CALL SPLINE(R550NR,Q880N1,880,REFF0,Q55DRY,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M1,880,REFF0,Q63DRY,1.D0,1.D0,1)

C     Find Q(RH),g(RH) paths in Q(m,r),g(m,r) maps for seed size = REFF0
C     2-coordinate paths defined via XN0=XNRRHX(I) & RR0=REFF0*RR0RHX(I)
C     ------------------------------------------------------------------
      DO 116 I=1,190
      XN0=XNRRHX(I)
      XN1=XN0*100.D0-129
      IN1=XN1
      DWN=XN1-IN1
      RR0=REFF0*RR0RHX(I)
      IF(RR0.LT.0.01) RR0=0.01
      IF(RR0.LE.3.00D0) XR1=RR0*100.D0+1
      IF(RR0.GT.3.00D0.AND.RR0.LT.3.04D0) XR1=RR0*50.0D0+732
      IF(RR0.GE.3.04D0.AND.RR0.LE.9.00D0) XR1=RR0*50.0D0+152
      IF(RR0.GT.9.00D0.AND.RR0.LT.9.08D0) XR1=RR0*25.0D0+662
      IF(RR0.GE.9.08D0) THEN
      XR1=RR0*25.0D0+378
      IF(XR1.GT.877.9999D0) XR1=877.9999D0
      ENDIF
      IR1=XR1
      DWR=XR1-IR1
      CALL SPLN44(Q633NR,890,31,IR1,DWR,IN1,DWN,QRH633(I))
      CALL SPLN44(G633NR,890,31,IR1,DWR,IN1,DWN,GRH633(I))
  116 CONTINUE

C     Define Q55(RH) by tracing path in Q(m,r) map for RH dependent size
C     via 2-coordinate path XN0=XNRRHX(I)*DNRX(I), RR0=RRH(I)*(.633/.55)
C     ------------------------------------------------------------------
      DO 117 I=1,190
      XN0=XNRRHX(I)*DNRX(I)
      XN1=XN0*100.D0-129
      IN1=XN1
      DWN=XN1-IN1
      RR0=REFF0*RR0RHX(I)*(0.633D0/0.550D0)
      IF(RR0.LT.0.01) RR0=0.01
      IF(RR0.LE.3.00D0) XR1=RR0*100.D0+1
      IF(RR0.GT.3.00D0.AND.RR0.LT.3.04D0) XR1=RR0*50.0D0+732
      IF(RR0.GE.3.04D0.AND.RR0.LE.9.00D0) XR1=RR0*50.0D0+152
      IF(RR0.GT.9.00D0.AND.RR0.LT.9.08D0) XR1=RR0*25.0D0+662
      IF(RR0.GE.9.08D0) THEN
      XR1=RR0*25.0D0+378
      IF(XR1.GT.877.9999D0) XR1=877.9999D0
      ENDIF
      IR1=XR1
      DWR=XR1-IR1
      CALL SPLN44(Q633NR,890,31,IR1,DWR,IN1,DWN,RHQ550(I))
      RHREFF(I)=RR0RHX(I)*REFF0
  117 CONTINUE

C        Aerosol liquid water content is in kg/m2 per unit optical depth
C                with the aerosol effective radius expressed in microns.
C                -------------------------------------------------------
      DO 118 I=1,190
      RHTAUF(I)=(RHQ550(I)/Q55DRY)*RR0RHX(I)**2 ! if input data are adj.
      AERMAS=1.33333333D0*RHREFF(I)*RHDENS(I)/RHQ550(I)*RHTAUF(I)
      RHTGM2(I)=AERMAS
      RHDGM2(I)=RHTGM2(1)
      RHWGM2(I)=RHTGM2(I)-RHDGM2(I)
  118 CONTINUE

      DO 120 N=1,9
      DO 119 I=1,190
      RHDATA(I,N)=RHINFO(I,N)
  119 CONTINUE
  120 CONTINUE

C     Determination of RH dependent Mie scattering tables for GCM input.
C     Find equivalent aersol dry sizes (RD1,RD2) and wet sizes (RW1,RW2)
C     and corresponding weights to match the RH dependent Q(r) and g(r).
C     Fits made to form: QRH=X*[Y*QD1+(1-Y)*QD2]+(1-X)*[Z*WD1+(1-Z)*WD2]
C     ------------------------------------------------------------------
      RFAC=1.25                                !   Spectral range factor
      J1=MAXWET
      JHIMAX=881-MAXWET
      K1=MAXDRY
      KHIMAX=881-MAXDRY
      NP=190-NRHN1+1
      DO 140 I=1,190
      RHRHI=RHRHRH(I)
      XFDRY=RHXMFX(I)
      REFFI=RHREFF(I)
      RRH=RR0RHX(I)*REFF0
      GRH=GRH633(I)
      QRH=QRH633(I)
      QD1=QRH
      QD2=QRH
      QW1=QRH
      QW2=QRH
      IF(QW1.GT.QQWMAX) QW1=QQWMAX
      IF(QW2.GT.QQWMAX) QW2=QQWMAX
      CALL SPLINV(R633NR    ,Q880M0    ,MAXWET,RW1,QW1,1.D0,1.D0,1)
      CALL SPLINV(R633NR(J1),Q880M0(J1),JHIMAX,RW2,QW2,1.D0,1.D0,1)
      CALL SPLINE(R633NR,G880M0,880,RW1,GW1,1.D0,1.D0,1)
      CALL SPLINE(R633NR,G880M0,880,RW2,GW2,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M0,880,RW1*RFAC,PW1,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M0,880,RW2*RFAC,PW2,1.D0,1.D0,1)
      IF(I.GE.NRHN1.AND.QRH.GT.QQWMAX) THEN
      QD1=QQWMAX+(QRH-QQWMAX)/XFDRY ! QD1 such that  QRH=X*QD1+(1-X)*QW1
      QD2=2.3D0                     ! 2 dry sizes are used if QD1>QQWMAX
      ENDIF
      CALL SPLINV(R633NR    ,Q880M1    ,MAXDRY,RD1,QD1,1.D0,1.D0,1)
      CALL SPLINV(R633NR(K1),Q880M1(K1),KHIMAX,RD2,QD2,1.D0,1.D0,1)
      CALL SPLINE(R633NR,G880M1,880,RD1,GD1,1.D0,1.D0,1)
      CALL SPLINE(R633NR,G880M1,880,RD2,GD2,1.D0,1.D0,1)
      RRHFAC=RRH*RFAC
      CALL SPLINE(RHREFF(NRHN1),QRH633(NRHN1),NP,RRHFAC,PRH,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M1,880,RD1*RFAC,PD1,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M1,880,RD2*RFAC,PD2,1.D0,1.D0,1)

      IF(I.LT.NRHN1) THEN         !          Pure dry aerosol region (1)
      WTX=1.D0
      WTY=1.D0
      IF(REFF0.GT.RQDMAX) WTY=0.D0
      WTZ=1.D0
      ELSE            !         Dry/wet weighted average regions (2)-(6)
      IF(QRH.LE.QQWMAX.AND.REFFI.LT.RW1) THEN
      IF(QRH.LT.2.2D0) THEN     !          Small-size aerosol region (2)
      WTZ=1.D0
      WTY=1.D0
      WTX=(GRH-GW1)/(GD1-GW1)
      ELSE                       !      Moderate-size aerosol region (3)
      WTZ=1.D0
      WTY=((PD2-PW1)*(GRH-GW1)-(PRH-PW1)*(GD2-GW1))
     +   /((PRH-PW1)*(GD1-GD2)-(PD1-PD2)*(GRH-GW1))
      WTX=(GRH-GW1)/(WTY*(GD1-GD2)+(GD2-GW1))
      ENDIF
      ENDIF
      IF(QRH.GT.QQWMAX) THEN       !      Medium-size aerosol region (4)
C     Fit form: QRH=X*(Y*QD1+(1-Y)*QD2)+(1-X)*QWmax & QRH=/QD1=/QD2=/QW1
      WTZ=1.D0
      WTY=((GRH-GD2)*QRH*QD2+(GD2-GW1)*QD2*QW1+(GW1-GRH)*QRH*QW1)
     +   /((GD1-GRH)*QRH*QD1+(GRH-GD2)*QRH*QD2+(GW1-GD1)*QD1*QW1
     +    +(GD2-GW1)*QD2*QW1)
      WTX=(QRH-QW1)/(WTY*(QD1-QD2)+(QD2-QW1))
      ENDIF
      IF(QRH.LE.QQWMAX.AND.REFFI.GT.RW1) THEN   !  Large size region (5)
      WTY=0.D0
      WTZ=((PD2-PW2)*(GRH-GW2)-(PRH-PW2)*(GD2-GW2))
     +   /((PW2-PW1)*(GD2-GRH)+(PRH-PD2)*(GW2-GW1))
      IF(WTZ.GT.1.D0) WTZ=1.D0
      IF(WTZ.LT.0.D0) WTZ=0.D0
      WTX=((GRH-GW2)+WTZ*(GW2-GW1))/((GD2-GW2)+WTZ*(GW2-GW1))
      ENDIF
      ENDIF
      IF(REFFI.GT.RQWMAX.AND.RHRHI.GT.0.995) THEN   ! High RH region (6)
      WTY=0.D0
      WTX=XFDRY
      WTZ=((GRH-GW2)-(GD2-GW2)*WTX)/((1.D0-WTX)*(GW1-GW2))
      ENDIF

      VD1=WTX*WTY
      VD2=WTX*(1.D0-WTY)
      VW1=WTZ*(1.D0-WTX)
      VW2=(1.D0-WTZ)*(1.D0-WTX)
      RD1=MIN(RD1,10.D0)
      RD2=MIN(RD2,10.D0)
      RW1=MIN(RW1,10.D0)
      RW2=MIN(RW2,10.D0)

C     Computed weight factors are for Lab reference wavelength of 633nm.
C     Rescale spectral extinction to 550 nm & renormalize weight factors
C     ------------------------------------------------------------------
      CALL SPLINE(R550NR,Q880N1,880,RD1,Q550,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M1,880,RD1,Q633,1.D0,1.D0,1)
      WD1=VD1*(Q550/Q633)
      CALL SPLINE(R550NR,Q880N1,880,RD2,Q550,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M1,880,RD2,Q633,1.D0,1.D0,1)
      WD2=VD2*(Q550/Q633)
      CALL SPLINE(R550NR,Q880N0,880,RW1,Q550,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M0,880,RW1,Q633,1.D0,1.D0,1)
      WW1=VW1*(Q550/Q633)
      CALL SPLINE(R550NR,Q880N0,880,RW2,Q550,1.D0,1.D0,1)
      CALL SPLINE(R633NR,Q880M0,880,RW2,Q633,1.D0,1.D0,1)
      WW2=VW2*(Q550/Q633)
      SUMW=WD1+WD2+WW1+WW2
      W1=WD1/SUMW
      W2=WD2/SUMW
      W3=WW1/SUMW
      W4=WW2/SUMW

C     ------------------------------------------------------------------
C     Tabulate relative humidity dependent solar, thermal Mie scattering
C     parameters SRHQEX,SRHQSC,SRHQCS, TRHQAB for each aerosol type NAER
C     These are mass weighted averages of equivalent dry and wet aerosol
C     parameters for sizes matching the relative humidity dependent Q(r)
C     ------------------------------------------------------------------

      N0=0                      !      Select Mie parameters for Sulfate
      IF(NAER.EQ.2) N0=22       !      Select Mie parameters for SeaSalt
      IF(NAER.EQ.3) N0=44       !      Select Mie parameters for Nitrate
      IF(NAER.EQ.4) N0=88       !      Select Mie parameters for Organic
      N1=N0+1
      DO 122 K=1,6                            !   SW dry sizes RD1 & RD2
      DO 121 N=1,22
      NN=N0+N
      QXAERN(N)=SRUQEX(K,NN)
      QSAERN(N)=SRUQSC(K,NN)
      QGAERN(N)=SRUQCB(K,NN)
  121 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,RD1,SR1QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RD1,SR1QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RD1,SR1QCB(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QXAERN,22,RD2,SR2QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RD2,SR2QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RD2,SR2QCB(K),1.D0,1.D0,1)
  122 CONTINUE

      DO 124 K=1,33                           !   LW dry sizes RD1 & RD2
      DO 123 N=1,22
      NN=N0+N
      QXAERN(N)=TRUQEX(K,NN)
      QSAERN(N)=TRUQSC(K,NN)
      QGAERN(N)=TRUQCB(K,NN)
  123 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,RD1,TR1QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RD1,TR1QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RD1,TR1QCB(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QXAERN,22,RD2,TR2QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RD2,TR2QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RD2,TR2QCB(K),1.D0,1.D0,1)
  124 CONTINUE
      CALL SPLINE(REFU22,Q55U22(N1),22,RD1,Q55RH1,1.D0,1.D0,1)
      CALL SPLINE(REFU22,Q55U22(N1),22,RD2,Q55RH2,1.D0,1.D0,1)

      N0=66                    !   Select Mie parameters for pure water
      N1=N0+1
      DO 126 K=1,6                           !   SW wet sizes RW1 & RW2
      DO 125 N=1,22
      NN=N0+N
      QXAERN(N)=SRUQEX(K,NN)
      QSAERN(N)=SRUQSC(K,NN)
      QGAERN(N)=SRUQCB(K,NN)
  125 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,RW1,SR3QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RW1,SR3QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RW1,SR3QCB(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QXAERN,22,RW2,SR4QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RW2,SR4QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RW2,SR4QCB(K),1.D0,1.D0,1)
  126 CONTINUE

      DO 128 K=1,33                          !   LW wet sizes RW1 & RW2
      DO 127 N=1,22
      NN=N0+N
      QXAERN(N)=TRUQEX(K,NN)
      QSAERN(N)=TRUQSC(K,NN)
      QGAERN(N)=TRUQCB(K,NN)
  127 CONTINUE
      CALL SPLINE(REFU22,QXAERN,22,RW1,TR3QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RW1,TR3QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RW1,TR3QCB(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QXAERN,22,RW2,TR4QEX(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QSAERN,22,RW2,TR4QSC(K),1.D0,1.D0,1)
      CALL SPLINE(REFU22,QGAERN,22,RW2,TR4QCB(K),1.D0,1.D0,1)
  128 CONTINUE
      CALL SPLINE(REFU22,Q55U22(N1),22,RW1,Q55RH3,1.D0,1.D0,1)
      CALL SPLINE(REFU22,Q55U22(N1),22,RW2,Q55RH4,1.D0,1.D0,1)

                        !       Weighted GCM SW Mie cattering parameters
      DO 129 K=1,6
      SRHQEX(K,I)=W1*SR1QEX(K)+W2*SR2QEX(K)+W3*SR3QEX(K)+W4*SR4QEX(K)
      SRHQSC(K,I)=W1*SR1QSC(K)+W2*SR2QSC(K)+W3*SR3QSC(K)+W4*SR4QSC(K)
      QSCQCB     =W1*SR1QCB(K)*SR1QSC(K)+W2*SR2QCB(K)*SR2QSC(K)
     +           +W3*SR3QCB(K)*SR3QSC(K)+W4*SR4QCB(K)*SR4QSC(K)
      SRHQCB(K,I)=QSCQCB/SRHQSC(K,I)
  129 CONTINUE
                        !       Weighted GCM LW Mie cattering parameters
      DO 130 K=1,33
      TRHQEX(K)=W1*TR1QEX(K)+W2*TR2QEX(K)+W3*TR3QEX(K)+W4*TR4QEX(K)
      TRHQSC(K)=W1*TR1QSC(K)+W2*TR2QSC(K)+W3*TR3QSC(K)+W4*TR4QSC(K)
      QSCQCB   =W1*TR1QCB(K)*TR1QSC(K)+W2*TR2QCB(K)*TR2QSC(K)
     +         +W3*TR3QCB(K)*TR3QSC(K)+W4*TR4QCB(K)*TR4QSC(K)
      TRHQCB(K)=QSCQCB/TRHQSC(K)
      TRHQAB(K,I)=TRHQEX(K)-TRHQSC(K)
  130 CONTINUE

!dbug Diagnostic variables  (not required for SETREL output)
!dbug QGFIT4=VD1*GD1*QD1+VD2*GD2*QD2+VW1*GW1*QW1+VW2*GW2*QW2
!dbug Q63FIT=VD1*QD1+VD2*QD2+VW1*QW1+VW2*QW2
!dbug G63FIT=QGFIT4/Q63FIT
!dbug Q63RAT=QRH633(I)/QRH633(1)
!dbug TAUF63=Q63RAT*RR0RHX(I)**2   !TAUF63 close to RHTAUF but not equal
!dbug GRATIO=GRH633(I)/G63FIT      !633nm weighted G-fit  equal to unity
!dbug Q55RH4=W1*Q55RH1+W2*Q55RH2+W3*Q55RH3+W4*Q55RH4   !  Weighted Q-fit
!dbug QRATIO=Q55RH4/RHQ550(I) ! only close to unity due to coarse interp

  140 CONTINUE

      RETURN
      END SUBROUTINE SETREL
