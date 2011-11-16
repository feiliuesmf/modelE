#include "rundeck_opts.h"


	  !@var  itype  1, ocean; 2, ocean ice; 3, land ice; 4, land
	  module itype_enum
	  	integer, parameter :: ITYPE_MIN = 1
	  	integer, parameter :: ITYPE_OCEAN = 1
	  	integer, parameter :: ITYPE_OCEANICE = 2
	  	integer, parameter :: ITYPE_LANDICE = 3
	  	integer, parameter :: ITYPE_LAND = 4
	  	integer, parameter :: ITYPE_MAX = 4
	  end module itype_enum


      MODULE EXCHANGE_TYPES ! todo: move to another file.
      use dist_grid_mod, only : dist_grid
#ifdef CUBED_SPHERE
      !USE cs2ll_utils, only : cs2llint_type,ll2csint_type
#else
      use domain_decomp_1d, only : band_pack_type
#endif
      IMPLICIT NONE

      type simple_bounds_type ! todo: move to another module
         INTEGER :: I_0,I_1, J_0,J_1  ! bounds of domain
         INTEGER :: I_0H,I_1H, J_0H,J_1H ! bounds of arrays
         LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
         INTEGER, DIMENSION(:), POINTER :: IMAXJ
      end type simple_bounds_type

      type, extends(simple_bounds_type) :: atmsrf_xchng_vars
!@var grid a pointer to the grid object whose domain bounds were
!@+   used to allocate an instance of this type
         type(dist_grid), pointer :: grid

         REAL*8, DIMENSION(:,:), POINTER ::
!@var E0 net energy flux at surface (J/m^2)
!@var SOLAR absorbed solar radiation (J/m^2)
     &      E0,SOLAR
C**** Momemtum stresses are calculated as if they were over whole box
!@var DMUA,DMVA momentum flux from atmosphere (kg/m s)
!@+   On atmospheric A grid (tracer point)
     &     ,DMUA, DMVA
!@var EVAPOR evaporation (kg/m^2) 
     &     ,PREC,EPREC,EVAPOR
!@var GTEMP temperature of surface (C)
!@var GTEMP2 "ground" temperature of "second" layer (C)
!@var GTEMPR radiative ground temperature over surface type (K)
     &     ,GTEMP,GTEMP2,GTEMPR
!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
!@var LAT latitude of gridbox (radians)
     &     ,COSZ1,LAT
!@var WSAVG     COMPOSITE SURFACE WIND MAGNITUDE (M/S) (over all types)
     &     ,WSAVG
#ifdef STANDALONE_OCEAN
!@var USAVG,VSAVG,TSAVG,QSAVG surface air wind components, temp, humidity
     &     ,USAVG,VSAVG,TSAVG,QSAVG ! in coupled model, PBL still owns
!@var FLONG, FSHORT downwelling longwave, shortwave radiation at surface
     &     ,FLONG,FSHORT
#endif
!@var SRFP actual surface pressure (hecto-Pascals)
     &     ,SRFP
#ifdef TRACERS_ON
!@var GTRACER ground concentration of tracer on atmospheric grid (kg/kg)
         REAL*8, DIMENSION(:,:,:), POINTER :: GTRACER
#ifdef TRACERS_WATER
         REAL*8, DIMENSION(:,:,:), POINTER ::
!@var TRPREC tracers in precip (kg/m^2)
     &      TRPREC
!@var TREVAPOR tracer evaporation over each type (kg/m^2) 
     &     ,TREVAPOR
#endif
#ifdef TRACERS_DRYDEP
!@var TRDRYDEP tracer dry deposition by type (kg/m^2) (positive down)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRDRYDEP
#endif
#endif
!@var WORK1,WORK2 temporary workspace
         REAL*8, DIMENSION(:,:), POINTER ::
     &      WORK1,WORK2

!@var AIJ a pointer to modelE diag accumulation arrays
         REAL*8, DIMENSION(:,:,:), POINTER :: AIJ

#ifdef TRACERS_ON
!@var TAIJN a pointer to modelE diag accumulation arrays
         REAL*8, DIMENSION(:,:,:,:), POINTER :: TAIJN
#endif

!@var JREG lat/lon array defining regions for modelE AREG diagnostics
         INTEGER, DIMENSION(:,:), POINTER :: JREG

#ifndef CUBED_SPHERE
!@var dlatm latitudinal gridsize in minutes, for certain regrid routines
         real*8 :: dlatm
!@var sini,cosi sin(lon),cos(lon) for vector regrids
         real*8, dimension(:), pointer :: sini,cosi
#endif

#ifdef TRACERS_ON
         integer :: ntm=0
#endif

      end type atmsrf_xchng_vars

      type, extends(atmsrf_xchng_vars) :: atmocn_xchng_vars
         REAL*8, DIMENSION(:,:), POINTER ::
!@var FOCEAN ocean fraction
     &      FOCEAN
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg/m^2, J/m^2)
     &     ,FLOWO, EFLOWO
!@var GMELT,EGMELT mass,energy from glacial melt into ocean (kg/m^2, J/m^2)
     &     ,GMELT, EGMELT
!@var UOSURF, VOSURF ocean surface velocity (cell center) (m/s)
!@+   components defined along true N/S and E/W directions
!@+   At the NP, U points from 90E to 90W, V from IDL to GM
!@+   At the SP, U points from 9OW to 90E, V from GM to IDL 
     &     ,UOSURF,VOSURF
!@var OGEOZA ocean surface height geopotential (m^2/s^2)
     &     ,OGEOZA
!@var MLHC ocean mixed layer heat capacity (J/m^2 C) 
     &     ,MLHC
!@var SSS sea surface salinity (ppt)
     &     ,SSS
#ifdef STANDALONE_OCEAN
!@var SSSOBS,SSSRESFAC observed salinity and restoration factor toward it
!@var RSIOBS observed sea ice fraction
     &     ,SSSOBS,SSSRESFAC
     &     ,RSIOBS
#endif
#ifdef TRACERS_WATER
!@var TRFLOWO tracer in river runoff into ocean (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRFLOWO
#ifdef TRACERS_OCEAN
!@var TRGMELT tracer from glacial melt into ocean (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRGMELT
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
!@var TRGASEX  tracer gas exchange (mol,CO2/m^2/s)
         REAL*8, DIMENSION(:,:,:), POINTER :: TRGASEX
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
         REAL*8, DIMENSION(:,:), POINTER :: pCO2
#endif
#ifdef OBIO_RAD_coupling
         REAL*8, DIMENSION(:,:), POINTER ::
     &     DIRVIS,DIFVIS,DIRNIR,DIFNIR
#endif
#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
C**** array of Chlorophyll data for use in ocean albedo calculation
!@var CHL Chlorophyll concentration data (mgr/m**3)
         REAL*8, DIMENSION(:,:), POINTER :: CHL
#endif

!@var eflow_gl global integral of eflowo
         real*8 :: eflow_gl
         logical :: need_eflow_gl=.false.

!@var modd5s,jm_budg,area_of_zone,conserv,nofm
!@+   permit an ocean model to accumulate modelE-style
!@+   conservation diagnostics.  see diag_com.
         integer :: modd5s,jm_budg
         real*8, dimension(:), pointer :: area_of_zone
         real*8, dimension(:,:), pointer :: consrv
         integer, dimension(:,:), pointer :: nofm

#ifdef TRACERS_ON
! Some atmosphere-declared tracer info for uses within ocean codes.
! See TRACER_COM.f
         real*8, dimension(:), pointer :: trw0
#ifdef TRACERS_GASEXCH_ocean
         integer :: ntm_gasexch=0
         real*8, dimension(:), pointer :: vol2mass
#endif
#endif

#ifdef TRACERS_OCEAN
!@var natmtrcons,tconsrv,nofmt
!@+   permit an ocean model to accumulate modelE-style
!@+   tracer conservation diagnostics.  see diag_com.
         integer :: natmtrcons
         real*8, dimension(:,:,:), pointer :: tconsrv
         integer, dimension(:,:), pointer :: nofmt
#endif

      end type atmocn_xchng_vars

      type, extends(atmsrf_xchng_vars) :: atmice_xchng_vars
         REAL*8, DIMENSION(:,:), POINTER ::
!@var FOCEAN ocean fraction
     &      FOCEAN
!@var E1 net energy flux at layer 1
     &     ,E1
!@var UISURF, VISURF dynamic ice surface velocity (Atm A grid) (m/s)
!@+   directions as for UOSURF/VOSURF
     &     ,UISURF,VISURF
!@var FWSIM fresh water sea ice mass (kg/m^2) (used for qflux model)
!@var HSICNV,MSICNV heat and fresh water sea ice mass convergence
!@+   after advsi (kg/m^2) (used for qflux model)
     &     ,FWSIM,MSICNV,HSICNV
!@var RSI fraction of water area covered in ice
!@var SNOWI snow amount on sea ice (kg/m^2)
     &     ,RSI,SNOWI
!@var USI,VSI ice velocities (m/s)
     &     ,USI,VSI ! temporary while ice still on atm grid
     &     ,SNOAGE ! really belongs to icestate

! Some arrays and indices for diagnostic purposes.  they are
! placed in this atm-ice interface type because the diagnostics
! are administrated by the atm model.  Many of them were
! introduced to maintain identicality of diagnostics through
! code reorganizations and can be eliminated if small changes
! in diagnostics are acceptable.
         REAL*8, DIMENSION(:,:), POINTER ::
     &     RSIstart,MSIsave,SNOWsave,TICEsave,TI1save,SSI1save,SSI2save
     &     ,SIHC,SNTOSI,SITOPMLT,MSNFLOOD,HSNFLOOD, SNOWsave2,MSIsave2
     &     ,MUSI,MVSI ,HUSI,HVSI ,SUSI,SVSI
         INTEGER :: IJ_RSNW,IJ_SNOW,IJ_RSIT,IJ_ZSNOW
         INTEGER :: IJ_MLTP,IJ_TSICE
     &     ,IJ_SIHC,IJ_SNTOSI,IJ_MSNFLOOD,IJ_HSNFLOOD,IJ_SIBOTMLT
     &     ,IJ_SMFX,IJ_SIGRLT,IJ_FWIO,IJ_HTIO,IJ_STIO,IJ_SIGRFR
     &     ,IJ_SIGRCG,IJ_SSI1,IJ_SSI2,IJ_TSI,IJ_F0OI,IJ_SISNWF
     &     ,IJ_RSOI,IJ_MSI,IJ_SITOPMLT
!@var IJ_[MHS][UV]SI indices for sea ice mass/heat/salt transport diags
     &     ,IJ_MUSI,IJ_MVSI,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI
         INTEGER :: J_IMELT,J_HMELT,J_SMELT
     &     ,j_implm,j_implh
     &     ,j_rsnow,j_rsi,j_ace1,j_ace2,j_snow
         INTEGER :: itoice,itlkice,itocean,itlake
#ifdef TRACERS_WATER
         LOGICAL, DIMENSION(:), ALLOCATABLE :: DO_ACCUM
         REAL*8, DIMENSION(:,:,:), POINTER :: TUSI,TVSI,TRSIsum
         integer :: tij_icocflx, tij_seaice, tij_tusi, tij_tvsi
#endif

      end type atmice_xchng_vars

      ! -----------------------------------------------------
      ! For coupling of atmosphere with glacial ice.
      type, extends(atmsrf_xchng_vars) :: atmgla_xchng_vars
         !@var E1 net energy flux at layer 1
         REAL*8, DIMENSION(:,:), POINTER :: E1
      end type atmgla_xchng_vars
      ! -----------------------------------------------------

      type, extends(atmsrf_xchng_vars) :: atmlnd_xchng_vars
!@var bare_soil_wetness bare_soil_wetness (1)
         REAL*8, DIMENSION(:,:), POINTER :: bare_soil_wetness
      end type atmlnd_xchng_vars

      type, extends(simple_bounds_type) :: iceocn_xchng_vars
!@var grid a pointer to the grid object whose domain bounds were
!@+   used to allocate an instance of this type
         type(dist_grid), pointer :: grid

         REAL*8, DIMENSION(:,:), POINTER ::
!@var FWATER water fraction of gridbox
     &      FWATER
!@var RSI fraction of water area covered in ice
     &     ,RSI
!@var CORIOL coriolis parameter (1/s)
     &     ,CORIOL
!@var SOLAR solar radiation penetrating the ice absorbed by ocean (J/m^2)
     &     ,SOLAR
!@var APRESS total atmos + sea ice pressure (at base of sea ice) (Pa)
     &     ,APRESS
!@var RUNOSI run off from sea/lake ice after surface (kg/m^2)
!@var ERUNOSI energy of run off from sea/lake ice after surface (J/m^2)
!@var SRUNOSI salt in run off from sea/lake ice after surface (kg/m^2)
     &     ,RUNOSI, ERUNOSI, SRUNOSI
!@var MELTI,EMELTI,SMELTI mass,energy,salt from simelt into ocn (kg/m^2,J/m^2)
     &     ,MELTI, EMELTI, SMELTI
!@var RUNPSI run off from sea/lake ice after precip (kg/m^2)
!@var ERUNPSI energy of run off from sea/lake ice after precip (J/m^2)
!@var SRUNPSI salt in run off from sea/lake ice after precip (kg/m^2)
     &     ,RUNPSI, ERUNPSI, SRUNPSI
!@var DMUI,DMVI momentum flux from sea ice to ocean (kg/m s)
!@+   On C grid for now
     &     ,DMUI, DMVI
!@var fmsi_io,fhsi_io,fssi_io basal ice-ocean fluxes (kg or J/m^2)
     &     ,fmsi_io,fhsi_io,fssi_io
!@var UI2rho Ustar*2*rho ice-ocean friction velocity on atmospheric grid
     &     ,UI2rho
!@var UOSURF,VOSURF ocean surface velocity (m/s)
     &     ,UOSURF,VOSURF
!@var OGEOZA ocean surface height geopotential (m^2/s^2)
     &     ,OGEOZA
!@var MLDLK mixed layer depth in lake (m)
!@var DLAKE depth of lake (m)
!@var GLAKE lake heat content (J/m2)
     &     ,mldlk,dlake,glake ! lakes only

C**** DMSI,DHSI,DSSI are fluxes for ice formation within water column
!@var DMSI mass flux of sea ice 1) open water and 2) under ice (kg/m^2)
!@var DHSI energy flux of sea ice 1) open water and 2) under ice (J/m^2)
!@var DSSI salt flux in sea ice 1) open water and 2) under ice (kg/m^2)
         REAL*8, DIMENSION(:,:,:), POINTER ::
     &      DMSI, DHSI, DSSI

#ifdef TRACERS_WATER
         REAL*8, DIMENSION(:,:,:), POINTER ::
!@var TRUNPSI tracer in run off from sea/lake ice after precip (kg/m^2)
!@var TRUNOSI tracer in run off from sea/lake ice after surface (kg/m^2)
!@var TRMELTI tracer from simelt into ocean (kg/m^2)
!@var ftrsi_io ice-ocean tracer fluxes under ice (kg/m^2)
     &     TRUNPSI, TRUNOSI, TRMELTI, ftrsi_io
#endif
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER) /* huh? */
!@var DTRSI tracer flux in sea ice under ice and on open water (kg/m^2)
         REAL*8, DIMENSION(:,:,:,:), POINTER :: DTRSI
#endif

!@var dlatm latitudinal gridsize in minutes, for certain regrid routines
         real*8 :: dlatm
#ifdef CUBED_SPHERE
#else
!@var pack_a2i,pack_i2a contain info for redistributing data from
!@+   atmos. domains to icedyn domains and vice versa.
         type(band_pack_type), pointer :: pack_a2i,pack_i2a
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
!@var ntm number of tracers participating in sea ice formation.
         integer :: ntm=0
#endif

      end type iceocn_xchng_vars

      interface alloc_xchng_vars
        module procedure alloc_atmsrf_xchng_vars
        module procedure alloc_atmocn_xchng_vars
        module procedure alloc_atmice_xchng_vars
        module procedure alloc_atmgla_xchng_vars
        module procedure alloc_atmlnd_xchng_vars
        module procedure alloc_iceocn_xchng_vars
      end interface alloc_xchng_vars

      CONTAINS

      subroutine set_simple_bounds_type(grd_dum,bds)
      USE DOMAIN_DECOMP_1D, ONLY : GET,DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(simple_bounds_type) :: BDS

      INTEGER :: I_0H,I_1H, J_0H,J_1H, I_0,I_1, J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(grd_dum,
     &     I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1,
     &     I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      bds % I_0H = I_0H
      bds % I_1H = I_1H
      bds % J_0H = J_0H
      bds % J_1H = J_1H

      bds % I_0 = I_0
      bds % I_1 = I_1
      bds % J_0 = J_0
      bds % J_1 = J_1

      bds % HAVE_SOUTH_POLE = HAVE_SOUTH_POLE
      bds % HAVE_NORTH_POLE = HAVE_NORTH_POLE

      ALLOCATE( bds % IMAXJ(J_0H:J_1H) )

      bds % IMAXJ(:) = I_1
      IF(HAVE_SOUTH_POLE) bds%IMAXJ(J_0) = 1
      IF(HAVE_NORTH_POLE) bds%IMAXJ(J_1) = 1

      return
      end subroutine set_simple_bounds_type

      subroutine alloc_atmsrf_xchng_vars(grd_dum,this)
      USE CONSTANT, only : tf
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmsrf_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_ON
      integer :: ntm
#endif

#ifdef TRACERS_ON
      ntm = this%ntm
#endif

      call set_simple_bounds_type(grd_dum,this%simple_bounds_type)

#ifndef CUBED_SPHERE
      allocate(this%sini(grd_dum%im_world),this%cosi(grd_dum%im_world))
#endif

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % E0      ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EVAPOR  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DMUA    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DMVA    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % PREC    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EPREC   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SOLAR   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % GTEMP   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % GTEMP2  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % GTEMPR  ( I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_ON
     &          this % GTRACER ( NTM, I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_WATER
     &          this % TRPREC  ( NTM, I_0H:I_1H , J_0H:J_1H ),
     &          this % TREVAPOR( NTM, I_0H:I_1H , J_0H:J_1H ),
#endif
#ifdef TRACERS_DRYDEP
     &          this % TRDRYDEP( NTM, I_0H:I_1H , J_0H:J_1H ),
#endif
#endif
     &          this % LAT     ( I_0H:I_1H , J_0H:J_1H ),
     &          this % COSZ1   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % WSAVG   ( I_0H:I_1H , J_0H:J_1H ),
#ifdef STANDALONE_OCEAN
     &          this % USAVG   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VSAVG   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % TSAVG   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % QSAVG   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FLONG   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FSHORT  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SRFP    ( I_0H:I_1H , J_0H:J_1H ),
#endif
     &          this % WORK1   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % WORK2   ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)

      this % GTEMP = 0.    ! initialize at 0 C
      this % GTEMP2 = 0.   ! initialize at 0 C
      this % GTEMPR = TF   ! initialize at 273 K
      this % EVAPOR = 0.

#ifdef TRACERS_ON
      this % GTRACER = 0.
#ifdef TRACERS_DRYDEP
      this % TRDRYDEP = 0.
#endif
#endif
      return
      end subroutine alloc_atmsrf_xchng_vars

      subroutine setptr_atmsrf_loopable(src,dest)
      IMPLICIT NONE
      TYPE(atmsrf_xchng_vars) :: src,dest
      dest%e0     => src%e0
      dest%evapor => src%evapor
      dest%solar  => src%solar
      dest%dmua   => src%dmua
      dest%dmva   => src%dmva
      dest%gtempr => src%gtempr
      dest%gtemp  => src%gtemp
      dest%gtemp2 => src%gtemp2
#ifdef TRACERS_ON
      dest%gtracer  => src%gtracer
#endif
#ifdef TRACERS_WATER
      dest%trevapor => src%trevapor
#endif
#ifdef TRACERS_DRYDEP
      dest%trdrydep => src%trdrydep
#endif
      return
      end subroutine setptr_atmsrf_loopable

      subroutine alloc_atmocn_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmocn_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_ON
      integer :: ntm
#endif
#ifdef TRACERS_GASEXCH_ocean
      integer :: ntm_gasexch
#endif

#ifdef TRACERS_ON
      ntm = this%ntm
#endif

#ifdef TRACERS_GASEXCH_ocean
      ntm_gasexch = this%ntm_gasexch
#endif

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % FOCEAN  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % GMELT   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EGMELT  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FLOWO   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EFLOWO  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % OGEOZA  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SSS     ( I_0H:I_1H , J_0H:J_1H ),
#ifdef STANDALONE_OCEAN
     &          this % SSSOBS  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SSSRESFAC ( I_0H:I_1H , J_0H:J_1H ),
     &          this % RSIOBS  ( I_0H:I_1H , J_0H:J_1H ),
#endif
     &          this % MLHC    ( I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_WATER
     &          this % TRFLOWO ( NTM , I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_OCEAN
     &          this % TRGMELT ( NTM , I_0H:I_1H , J_0H:J_1H ),
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
     &          this % TRGASEX ( NTM_gasexch , I_0H:I_1H , J_0H:J_1H ),
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
     &          this % pCO2    ( I_0H:I_1H , J_0H:J_1H ),
#endif
#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
     &          this % CHL     ( I_0H:I_1H , J_0H:J_1H ),
#endif
     &   STAT = IER)

      this % UOSURF = 0.
      this % VOSURF = 0.
      this % OGEOZA = 0.

#ifdef TRACERS_GASEXCH_ocean
      this % TRGASEX = 0.
! I dimensioned trgasex by ntm_gasexch to avoid confusion elsewhere
! in the code.  If ntm differs from ntm_gasexch, fluxes need to be
! stored in the appropriate positions in trgasex or this%trgasex.  - MK
      if(ntm /= ntm_gasexch) call stop_model(
     &     'alloc_atmocn_xchng_vars: ntm /= ntm_gasexch',255)
#endif

#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
      this % CHL = 0.
#endif

#ifdef OBIO_RAD_coupling
      allocate(
     &          this % DIRVIS  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DIFVIS  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DIRNIR  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DIFNIR  ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif

      this % modd5s = -999

      return
      end subroutine alloc_atmocn_xchng_vars

      subroutine alloc_atmice_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmice_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#ifdef TRACERS_WATER
      integer :: ntm
#endif

#ifdef TRACERS_WATER
      ntm = this%ntm
#endif

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % E1 ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UISURF ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VISURF ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FWSIM  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MSICNV ( I_0H:I_1H , J_0H:J_1H ),
     &          this % HSICNV ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FOCEAN ( I_0H:I_1H , J_0H:J_1H ),
     &          this % USI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VSI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MUSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MVSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % HUSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % HVSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SUSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SVSI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SNOAGE ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      this % UISURF = 0.
      this % VISURF = 0.
      this % MSICNV = 0.
      this % HSICNV = 0.
      this % SNOAGE = 0.

#ifdef TRACERS_WATER
      ALLOCATE( this % TUSI   (I_0H:I_1H, J_0H:J_1H, NTM),
     &          this % TVSI   (I_0H:I_1H, J_0H:J_1H, NTM),
     &          this % TRSIsum(NTM, I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)
#endif

      ALLOCATE( this % MSIsave(I_0H:I_1H, J_0H:J_1H),
     &          this % SNOWsave(I_0H:I_1H, J_0H:J_1H),
     &          this % TICEsave(I_0H:I_1H, J_0H:J_1H),
     &          this % TI1save(I_0H:I_1H, J_0H:J_1H),
     &          this % SSI1save(I_0H:I_1H, J_0H:J_1H),
     &          this % SSI2save(I_0H:I_1H, J_0H:J_1H),
     &          this % SIHC(I_0H:I_1H, J_0H:J_1H),
     &          this % RSIstart(I_0H:I_1H, J_0H:J_1H),
     &          this % SNTOSI(I_0H:I_1H, J_0H:J_1H),
     &          this % SITOPMLT(I_0H:I_1H, J_0H:J_1H),
     &          this % MSNFLOOD(I_0H:I_1H, J_0H:J_1H),
     &          this % HSNFLOOD(I_0H:I_1H, J_0H:J_1H),
     &          this % MSIsave2(I_0H:I_1H, J_0H:J_1H),
     &          this % SNOWsave2(I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)

      this % MSIsave = 0.
      this % SNOWsave = 0.
      this % TICEsave = 0.
      this % TI1save = 0.
      this % SIHC = 0.
      this % RSIstart = 0.
      this % SNTOSI = 0.
      this % SITOPMLT = 0.
      this % MSNFLOOD = 0.
      this % HSNFLOOD = 0.

      return
      end subroutine alloc_atmice_xchng_vars

      subroutine alloc_atmgla_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmgla_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % E1 ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      return
      end subroutine alloc_atmgla_xchng_vars

      subroutine alloc_atmlnd_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(atmlnd_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      call alloc_atmsrf_xchng_vars(grd_dum,this%atmsrf_xchng_vars)
      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO
      ALLOCATE( this % bare_soil_wetness ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
      return
      end subroutine alloc_atmlnd_xchng_vars

      subroutine alloc_iceocn_xchng_vars(grd_dum,this)
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      TYPE(iceocn_xchng_vars) :: THIS
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      integer :: ntm
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      ntm = this%ntm
#endif


      call set_simple_bounds_type(grd_dum,this%simple_bounds_type)

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      ALLOCATE( this % RUNOSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          this % ERUNOSI ( I_0H:I_1H , J_0H:J_1H ), 
     &          this % SRUNOSI ( I_0H:I_1H , J_0H:J_1H ),
     &          this % RUNPSI  ( I_0H:I_1H , J_0H:J_1H ), 
     &          this % SRUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          this % ERUNPSI ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DMUI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % DMVI    ( I_0H:I_1H , J_0H:J_1H ),
     &          this % MELTI   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % EMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SMELTI  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % APRESS  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % RSI     ( I_0H:I_1H , J_0H:J_1H ),
     &          this % fmsi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          this % fhsi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          this % fssi_io ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UI2rho  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % SOLAR   ( I_0H:I_1H , J_0H:J_1H ),
     &          this % FWATER  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % CORIOL  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % UOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % VOSURF  ( I_0H:I_1H , J_0H:J_1H ),
     &          this % OGEOZA  ( I_0H:I_1H , J_0H:J_1H ),
!     &          this % MLHC    ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER)

       ALLOCATE( this % DMSI  (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           this % DHSI  (  2  , I_0H:I_1H , J_0H:J_1H ), 
     &           this % DSSI  (  2  , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
       this % DMSI = 0.
       this % DHSI = 0.
       this % DSSI = 0.

#ifdef TRACERS_WATER
      ALLOCATE( this % TRUNPSI (NTM, I_0H:I_1H, J_0H:J_1H),
     &          this % TRUNOSI (NTM, I_0H:I_1H, J_0H:J_1H),
     &          this % TRMELTI (NTM, I_0H:I_1H, J_0H:J_1H),
     &          this % FTRSI_IO(NTM, I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)
      this % TRUNPSI = 0.
      this % TRUNOSI = 0.
      this % TRMELTI = 0.
      this % FTRSI_IO = 0.
#endif
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER) /* huh? */
      ALLOCATE( this % DTRSI(NTM, 2, I_0H:I_1H, J_0H:J_1H),
     &   STAT = IER)
      this % DTRSI = 0.
#endif

      return
      end subroutine alloc_iceocn_xchng_vars

      END MODULE EXCHANGE_TYPES

#ifndef STANDALONE_OCEAN

      MODULE FLUXES
!@sum  FLUXES contains the fluxes between various atm-grid components
!@auth Gavin Schmidt
      USE RESOLUTION, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,atmice_xchng_vars,
     &     atmsrf_xchng_vars,atmlnd_xchng_vars,atmgla_xchng_vars
#ifdef TRACERS_ON
      USE TRACER_COM, only: NTM
#ifndef SKIP_TRACER_SRCS
     *     ,ntsurfsrcmax,nt3Dsrcmax
#endif
#endif
      IMPLICIT NONE

!@dbparam NIsurf: DT_Surface  =  DTsrc/NIsurf
      INTEGER :: NIsurf = 2

!@var Fxx fraction of gridbox of type xx (land,ocean,...)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAND
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FOCEAN
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLICE
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FEARTH0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FLAKE0


!@var RUNOE run off from earth (kg/m^2)
!@var ERUNOE energy of run off from earth (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNOE, ERUNOE

!@var RUNOLI run off from land ice (kg/m^2) (Energy always=0)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNOLI

!@param NSTYPE number of surface types for radiation purposes
      INTEGER, PARAMETER :: NSTYPE=4

!@dbparam UOdrag parameter that decides whether ocean.ice velocities
!@+   feed into drag calculation in surface (default = 0)
      INTEGER :: UOdrag = 0

C**** currently saved - should be replaced by fluxed quantities
!@var DTH1,DQ1 heat/water flux from atmos. summed over type (C, kg/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DTH1,DQ1

!@var uflux1 surface turbulent u-flux (=-<uw>) 
!@var vflux1 surface turbulent v-flux (=-<vw>)
!@var tflux1 surface turbulent t-flux (=-<tw>)
!@var qflux1 surface turbulent q-flux (=-<qw>)
      real*8, allocatable, dimension(:,:) :: 
     &        uflux1,vflux1,tflux1,qflux1

C**** The E/FLOWO, E/S/MELTI, E/GMELT arrays are used to flux quantities 
C**** to the ocean that are not tied to the open water/ice covered 
C**** fractions. This is done separately for river flow, complete
C**** sea ice melt and iceberg/glacial melt.
!@var FLOWO,EFLOWO mass, energy from rivers into ocean (kg/m^2, J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FLOWO,EFLOWO
!@var GMELT,EGMELT mass,energy from glacial melt into ocean (kg/m^2,J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: GMELT,EGMELT

!@var PREC precipitation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PREC
!@var EPREC energy of preciptiation (J/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PRECSS

#ifdef IRRIGATION_ON
!@var Actual irrigation rate (& energy) [m/s] [W/m2]
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: irrig_water_act
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:) :: irrig_energy_act
#ifdef TRACERS_WATER
!@var Actual irrigation tracer rate [kg/s] 
      REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: irrig_tracer_act
#endif
#endif

C**** fluxes associated with variable lake fractions
!@var DMWLDF  water deficit over land surface (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DMWLDF
!@var DGML energy associated with DMWLDF (J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DGML

#ifdef TRACERS_ON
!@var TRSOURCE non-interactive surface sources/sinks for tracers (kg/s)
#ifndef SKIP_TRACER_SRCS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trsource
#endif
!@var TRSRFFLX interactive surface sources/sinks for tracers (kg/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trsrfflx
!@var TRFLUX1 total surface flux for each tracer (kg/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trflux1
!@var TR3DSOURCE 3D sources/sinks for tracers (kg/s)
#ifndef SKIP_TRACER_SRCS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:):: tr3Dsource
#endif

#ifdef TRACERS_WATER
!@var TRPREC tracers in precip (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: TRPREC
!@var TRUNOE tracer runoff from earth (kg/m^2)
!@var TRUNOLI tracer runoff from land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: 
     &        TRUNOE, TRUNOLI

!@var TRFLOWO tracer in river runoff into ocean (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRFLOWO

C**** fluxes associated with variable lake fractions
!@var DTRL tracers associate with DMWLDF (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DTRL

#ifdef TRACERS_OCEAN
!@var TRGMELT tracer from glacial melt into ocean (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRGMELT
#endif
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
!@var trprec_dust dust/mineral tracers in precip [kg]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:):: trprec_dust
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
!@var pprec precipitation at previous time step [kg/m^2]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: pprec
!@var pevap evaporation at previous time step [kg/m^2]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: pevap
!@var dust_flux_glob global array of dust emission flux [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dust_flux_glob
#endif
!@var dust_flux2_glob global array of cubic dust emission flux (for diags only)
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: dust_flux2_glob

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
#ifdef TRACERS_DRYDEP
!@var depo_turb_glob global array of flux due to dry turb. dep. of tracers
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: depo_turb_glob
!@var depo_grav_glob global array of flux due to gravit. settling of tracers
!@+   [kg/m^2/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: depo_grav_glob
#endif
#endif

#endif

!@var atmocn,atmic,atmgla,atmlnd derived-type strucures containing
!@+   variables (or pointers thereto) needed for atmospheric
!@+   interactions with water, floating ice, glacial ice, and the
!@+   land surface.
      type(atmocn_xchng_vars) :: atmocn ! ocean and lakes
      type(atmice_xchng_vars) :: atmice ! ocean and lakes
      type(atmgla_xchng_vars) :: atmgla ! glacial ice
      type(atmlnd_xchng_vars) :: atmlnd ! land surface

!@var asflx an array for looping over atmocn,atmice,atmgla,atmlnd
      type(atmsrf_xchng_vars), dimension(4) :: asflx

      TARGET :: 
     &      FOCEAN, PREC, EPREC
     &     ,FLOWO, EFLOWO, GMELT, EGMELT
#ifdef TRACERS_WATER
      TARGET :: TRPREC,TRFLOWO
#ifdef TRACERS_OCEAN
      TARGET :: TRGMELT
#endif
#endif

      END MODULE FLUXES

      SUBROUTINE ALLOC_FLUXES !(grd_dum)
!@sum   Initializes FLUXES''s arrays
!@auth  Rosalinda de Fainchtein
      USE FILEMANAGER
      USE EXCHANGE_TYPES, only : alloc_xchng_vars,setptr_atmsrf_loopable
      USE DOMAIN_DECOMP_ATM, ONLY : GRD_DUM=>GRID,READT_PARALLEL,
     &     HALO_UPDATE,HASSOUTHPOLE,HASNORTHPOLE
      USE FLUXES
      USE GEOM, only : lat2d
#ifndef CUBED_SPHERE
      USE GEOM, only : dlatm,sinip,cosip
#endif
#ifdef TRACERS_ON
      USE tracer_com,ONLY : Ntm=> NTM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     &     ,Ntm_dust
#endif
#ifdef TRACERS_GASEXCH_ocean
      use tracer_com, only : ntm_gasexch
#endif
#endif
      USE ATM_COM, only : srfp
      USE Dictionary_mod
      IMPLICIT NONE
      !TYPE (DIST_GRID), INTENT(IN) :: grd_dum
      INTEGER :: iu_TOPO
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: I, J, I_0, I_1, J_1, J_0
      INTEGER :: IER

      call sync_param( "NIsurf", NIsurf )
      call sync_param( "UOdrag", UOdrag )

      I_0H = grd_dum%I_STRT_HALO
      I_1H = grd_dum%I_STOP_HALO
      J_0H = grd_dum%J_STRT_HALO
      J_1H = grd_dum%J_STOP_HALO

      I_0 = grd_dum%I_STRT
      I_1 = grd_dum%I_STOP
      J_0 = grd_dum%J_STRT
      J_1 = grd_dum%J_STOP

      ALLOCATE(FLAND(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FOCEAN(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FLICE(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FEARTH0(I_0H:I_1H,J_0H:J_1H), STAT = IER)
      ALLOCATE(FLAKE0(I_0H:I_1H,J_0H:J_1H), STAT = IER)

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
C**** Note that FLAKE0 is read in only to provide initial values
C**** Actual array is set from restart file.
      call openunit("TOPO",iu_TOPO,.true.,.true.)

      CALL READT_PARALLEL(grd_dum,iu_TOPO,NAMEUNIT(iu_TOPO),FOCEAN,1) ! Ocean fraction
      CALL HALO_UPDATE(GRD_DUM, FOCEAN)

      CALL READT_PARALLEL(grd_dum,iu_TOPO,NAMEUNIT(iu_TOPO),FLAKE0,1) ! Orig. Lake fraction
      CALL READT_PARALLEL(grd_dum,iu_TOPO,NAMEUNIT(iu_TOPO),FEARTH0,1) ! Earth frac. (no LI)
      CALL HALO_UPDATE(GRD_DUM, FEARTH0)

      CALL READT_PARALLEL(grd_dum,iu_TOPO,NAMEUNIT(iu_TOPO),FLICE ,1) ! Land ice fraction

C**** Deal with single -> double precision problems and potential
C**** ocean/lake inconsistency. Adjust FLAKE0 and FLICE if necessary.
      DO J=J_0,J_1
      DO I=I_0,I_1 !IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          FLAND(I,J)=1.-FOCEAN(I,J) ! Land fraction if focean>0
          IF (FLAKE0(I,J).gt.0) THEN
            WRITE(6,*) "Ocean and lake cannot co-exist in same grid box"
     *       ,i,j,FOCEAN(I,J),FLAKE0(I,J)
            FLAKE0(I,J)=0
          END IF
        ELSEIF (FLAKE0(I,J).gt.0) THEN
          FLAND(I,J)=1.-FLAKE0(I,J)  ! for initialization only
        ELSE
          FLAND(I,J)=1.              ! for initialization only
        END IF
C**** Ensure that no round off error effects land with ice and earth
        IF (FLICE(I,J)-FLAND(I,J).gt.-1d-4 .and. FLICE(I,J).gt.0) THEN
          FLICE(I,J)=FLAND(I,J)
        END IF
      END DO
      END DO
      CALL HALO_UPDATE(GRD_DUM, FLAKE0)
      CALL HALO_UPDATE(GRD_DUM, FLICE)
      call closeunit(iu_TOPO)
      If (HASSOUTHPOLE(GRD_DUM)) Then
         FLAND(2:IM,1)=FLAND(1,1)
         FLICE(2:IM,1)=FLICE(1,1)
      End If
      If (HASNORTHPOLE(GRD_DUM)) Then
         FLAND(2:IM,JM)=FLAND(1,JM)
         FLICE(2:IM,JM)=FLICE(1,JM)
      End If

      !I-J arrays
      ALLOCATE(
     &          RUNOE   ( I_0H:I_1H , J_0H:J_1H ), 
     &          ERUNOE  ( I_0H:I_1H , J_0H:J_1H ),
     &          RUNOLI  ( I_0H:I_1H , J_0H:J_1H ),
     &          DTH1    ( I_0H:I_1H , J_0H:J_1H ),
     &          DQ1     ( I_0H:I_1H , J_0H:J_1H ),
     &   STAT=IER )
      ALLOCATE( uflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          vflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          tflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          qflux1  ( I_0H:I_1H , J_0H:J_1H ),
     &          FLOWO   ( I_0H:I_1H , J_0H:J_1H ),
     &          EFLOWO  ( I_0H:I_1H , J_0H:J_1H ),
     &          GMELT   ( I_0H:I_1H , J_0H:J_1H ),
     &          EGMELT  ( I_0H:I_1H , J_0H:J_1H ),
     &          PREC    ( I_0H:I_1H , J_0H:J_1H ),
     &          EPREC   ( I_0H:I_1H , J_0H:J_1H ),
     &          PRECSS  ( I_0H:I_1H , J_0H:J_1H ),
     &          DMWLDF  ( I_0H:I_1H , J_0H:J_1H ),
     &          DGML    ( I_0H:I_1H , J_0H:J_1H ),
#ifdef IRRIGATION_ON
     &          irrig_water_act ( I_0H:I_1H , J_0H:J_1H ),
     &          irrig_energy_act( I_0H:I_1H , J_0H:J_1H ),
#ifdef TRACERS_WATER
     &          irrig_tracer_act(NTM, I_0H:I_1H , J_0H:J_1H),
#endif
#endif
     &   STAT=IER)

!TRACERS_ON**********************

#ifdef TRACERS_ON
      !(I,J,:,:)  array
#ifndef SKIP_TRACER_SRCS
      ALLOCATE(trsource (I_0H:I_1H,J_0H:J_1H,ntsurfsrcmax,NTM)
     &  ,STAT = IER)
      trsource = 0.
#endif

      !(I,J,:) arrays
      ALLOCATE( trsrfflx( I_0H:I_1H , J_0H:J_1H , NTM    ),
     &          trflux1 ( I_0H:I_1H , J_0H:J_1H , NTM    ),
     &   STAT = IER)

      !I-J-L-:-: array
#ifndef SKIP_TRACER_SRCS
      ALLOCATE( tr3Dsource(I_0H:I_1H,J_0H:J_1H,LM,nt3Dsrcmax,NTM)
     &  ,STAT = IER)
#endif

#ifdef TRACERS_WATER
                                                    !(:)-(:)-I-J arrays
       !:-I-J arrays
       ALLOCATE( TRPREC  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRUNOE  ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRUNOLI ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           TRFLOWO ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &           DTRL    ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#ifdef TRACERS_OCEAN
       ALLOCATE( TRGMELT ( NTM , I_0H:I_1H , J_0H:J_1H ),
     &   STAT = IER)
#endif
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      ALLOCATE(trprec_dust(Ntm_dust,I_0H:I_1H ,J_0H:J_1H),STAT=ier)
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      ALLOCATE(pprec(I_0H:I_1H,J_0H:J_1H),STAT = IER)
      ALLOCATE(pevap(I_0H:I_1H,J_0H:J_1H,NSTYPE),STAT = IER)
      ALLOCATE(dust_flux_glob(I_0H:I_1H,J_0H:J_1H,Ntm_dust),STAT = IER)
#ifdef TRACERS_DRYDEP
      ALLOCATE(depo_turb_glob(I_0H:I_1H,J_0H:J_1H,Nstype,Ntm)
     &     ,STAT = IER)
      ALLOCATE(depo_grav_glob(I_0H:I_1H,J_0H:J_1H,Nstype,Ntm)
     &     ,STAT = IER)
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      ALLOCATE(dust_flux2_glob(I_0H:I_1H,J_0H:J_1H,Ntm_dust),STAT = IER)
#endif

#endif

#ifdef TRACERS_ON
      atmocn%ntm = ntm
      atmice%ntm = ntm
      atmgla%ntm = ntm
      atmlnd%ntm = ntm
#endif
#ifdef TRACERS_GASEXCH_ocean
      atmocn%ntm_gasexch = ntm_gasexch
#endif
      call alloc_xchng_vars(grd_dum,atmocn)
      atmocn%grid => grd_dum
      atmocn%focean(:,:) = focean(:,:)
      atmocn%lat(:,:) = lat2d(:,:)
#ifndef CUBED_SPHERE
      atmocn%dlatm = dlatm
      atmocn%sini(:) = sinip(:)
      atmocn%cosi(:) = cosip(:)
#endif
      deallocate(atmocn%prec);  atmocn%prec => prec
      deallocate(atmocn%eprec); atmocn%eprec => eprec
      deallocate(atmocn%flowo);  atmocn%flowo => flowo
      deallocate(atmocn%eflowo); atmocn%eflowo => eflowo
      deallocate(atmocn%gmelt);  atmocn%gmelt => gmelt
      deallocate(atmocn%egmelt); atmocn%egmelt => egmelt
#ifdef TRACERS_WATER
      deallocate(atmocn%trprec);  atmocn%trprec => trprec
      deallocate(atmocn%trflowo);  atmocn%trflowo => trflowo
#ifdef TRACERS_OCEAN
      deallocate(atmocn%trgmelt);  atmocn%trgmelt => trgmelt
#endif
#endif
      atmocn%srfp => srfp

      call alloc_xchng_vars(grid,atmice)
      atmice%grid => grd_dum
      atmice%focean(:,:) = focean(:,:)
      atmice%lat(:,:) = lat2d(:,:)
#ifndef CUBED_SPHERE
      atmice%dlatm = dlatm
      atmice%sini(:) = sinip(:)
      atmice%cosi(:) = cosip(:)
#endif
      deallocate(atmice%prec);  atmice%prec => prec
      deallocate(atmice%eprec); atmice%eprec => eprec
#ifdef TRACERS_WATER
      deallocate(atmice%trprec);  atmice%trprec => trprec
#endif
      atmice%srfp => srfp

      call alloc_xchng_vars(grid,atmgla)
      atmgla%grid => grd_dum
      atmgla%srfp => srfp
      atmgla%lat(:,:) = lat2d(:,:)

      call alloc_xchng_vars(grid,atmlnd)
      atmlnd%grid => grd_dum
      atmlnd%srfp => srfp
      atmlnd%lat(:,:) = lat2d(:,:)

! set pointers for looping over surface types
      call setptr_atmsrf_loopable(atmocn%atmsrf_xchng_vars,asflx(1))
      call setptr_atmsrf_loopable(atmice%atmsrf_xchng_vars,asflx(2))
      call setptr_atmsrf_loopable(atmgla%atmsrf_xchng_vars,asflx(3))
      call setptr_atmsrf_loopable(atmlnd%atmsrf_xchng_vars,asflx(4))

      END SUBROUTINE ALLOC_FLUXES

      subroutine def_rsf_fluxes(fid)
!@sum  def_rsf_fluxes defines structure of coupler arrays in restart files
!@auth M. Kelley
!@ver  beta
      use fluxes
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      character(len=17) :: ijstr
      ijstr='(dist_im,dist_jm)'
      !call defvar(agrid,fid,atmocn%gtemp,'gtemp'//ijstr)
      !call defvar(agrid,fid,atmocn%gtempr,'gtempr'//ijstr)
      call defvar(grid,fid,atmocn%gtemp,'asst'//ijstr)
      call defvar(grid,fid,atmocn%gtempr,'atempr'//ijstr)
      call defvar(grid,fid,atmocn%sss,'sss'//ijstr)
      call defvar(grid,fid,atmocn%ogeoza,'ogeoza'//ijstr)
      call defvar(grid,fid,atmocn%uosurf,'uosurf'//ijstr)
      call defvar(grid,fid,atmocn%vosurf,'vosurf'//ijstr)
      return
      end subroutine def_rsf_fluxes

      subroutine new_io_fluxes(fid,iaction)
!@sum  new_io_fluxes read/write coupler arrays from/to restart files
!@auth M. Kelley
!@ver  beta
      use model_com, only : iowrite,ioread
      use fluxes
      use domain_decomp_atm, only: grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        !call write_dist_data(grid,fid,'gtemp',atmocn%gtemp)
        !call write_dist_data(grid,fid,'gtempr',atmocn%gtempr)
        call write_dist_data(grid,fid,'asst',atmocn%gtemp)
        call write_dist_data(grid,fid,'atempr',atmocn%gtempr)
        call write_dist_data(grid,fid,'sss',atmocn%sss)
        call write_dist_data(grid,fid,'ogeoza',atmocn%ogeoza)
        call write_dist_data(grid,fid,'uosurf',atmocn%uosurf)
        call write_dist_data(grid,fid,'vosurf',atmocn%vosurf)
      case (ioread)             ! input from restart file
        !call read_dist_data(grid,fid,'gtemp',atmocn%gtemp)
        !call read_dist_data(grid,fid,'gtempr',atmocn%gtempr)
        call read_dist_data(grid,fid,'asst',atmocn%gtemp)
        call read_dist_data(grid,fid,'atempr',atmocn%gtempr)
        call read_dist_data(grid,fid,'sss',atmocn%sss)
        call read_dist_data(grid,fid,'ogeoza',atmocn%ogeoza)
        call read_dist_data(grid,fid,'uosurf',atmocn%uosurf)
        call read_dist_data(grid,fid,'vosurf',atmocn%vosurf)
      end select
      return
      end subroutine new_io_fluxes

#endif /* not STANDALONE_OCEAN */
