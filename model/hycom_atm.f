#include "rundeck_opts.h"

#if defined(CUBED_SPHERE) || defined(NEW_IO)
#else
#define USE_ATM_GLOBAL_ARRAYS
#endif

      module hycom_atm
!@sum module for atmospheric variables to be passed to/from hylom.
!@+   hycom will see them as global arrays
      USE MODEL_COM, only : im,jm
      USE HYCOM_DIM, only : aI_0, aI_1, aI_0H, aI_1H
      USE HYCOM_DIM, only : aJ_0, aJ_1, aJ_0H, aJ_1H
      use FLUXES, only: NSTYPE

      USE FLUXES, only : PREC_loc => PREC
      USE FLUXES, only : EVAPOR_loc => EVAPOR
      USE FLUXES, only : FLOWO_loc => FLOWO
      USE FLUXES, only : GMELT_loc => GMELT
      USE FLUXES, only : MELTI_loc => MELTI
      USE FLUXES, only : RUNOSI_loc => RUNOSI
      USE FLUXES, only : RUNPSI_loc => RUNPSI
      USE FLUXES, only : E0_loc => E0
      USE FLUXES, only : EPREC_loc => EPREC
      USE FLUXES, only : EFLOWO_loc => EFLOWO
      USE FLUXES, only : EGMELT_loc => EGMELT
      USE FLUXES, only : EMELTI_loc => EMELTI
      USE FLUXES, only : ERUNOSI_loc => ERUNOSI
      USE FLUXES, only : ERUNPSI_loc => ERUNPSI
      USE FLUXES, only : SRUNOSI_loc => SRUNOSI
      USE FLUXES, only : SRUNPSI_loc => SRUNPSI
      USE FLUXES, only : SMELTI_loc => SMELTI
      USE FLUXES, only : DMUA_loc => DMUA
      USE FLUXES, only : DMVA_loc => DMVA
      USE ICEDYN, only : grid_icdyn
      USE ICEDYN_COM, only : pack_i2a ! to send dmui,dmvi to the atm grid
      USE ICEDYN_COM, only : pack_a2i
      USE ICEDYN_COM, only : DMUI_loc => DMUI
      USE ICEDYN_COM, only : DMVI_loc => DMVI
      USE FLUXES, only : SOLAR_loc => SOLAR
      USE FLUXES, only : SSS_loc => SSS
      USE FLUXES, only : UOSURF_loc => UOSURF
      USE FLUXES, only : VOSURF_loc => VOSURF
      USE ICEDYN_COM, only : UOSURF_4DYNSI_loc => UOSURF
      USE ICEDYN_COM, only : VOSURF_4DYNSI_loc => VOSURF
      USE FLUXES, only : OGEOZA_loc => OGEOZA
      USE FLUXES, only : GTEMP_loc => GTEMP
      USE FLUXES, only : GTEMPR_loc => GTEMPR
      USE FLUXES, only : DMSI_loc => DMSI
      USE FLUXES, only : DHSI_loc => DHSI
      USE FLUXES, only : DSSI_loc => DSSI

      USE SEAICE_COM, only : RSI_loc => RSI ! seems to be used for diags only?
      USE MODEL_COM, only : FOCEAN_loc => FOCEAN

#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_COM, only : ntm
      USE FLUXES, only : GTRACER_loc => GTRACER, TRGASEX_loc => TRGASEX
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : COSZ1_loc => COSZ1
      USE PBLCOM, only : wsavg_loc => wsavg
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only: FSRDIR_loc => FSRDIR
     .                  ,FSRDIF_loc => FSRDIF
     .                  ,DIRNIR_loc => DIRNIR
     .                  ,DIFNIR_loc => DIFNIR
     .                  ,SRVISSURF_loc => SRVISSURF
#endif
#ifdef pCO2_ONLINE
      USE PBL_DRV, only : t_pbl_args
#endif


      implicit none
      private

      public PREC_loc
      public EVAPOR_loc
      public FLOWO_loc
      public GMELT_loc
      public MELTI_loc
      public RUNOSI_loc
      public RUNPSI_loc
      public E0_loc
      public EPREC_loc
      public EFLOWO_loc
      public EGMELT_loc
      public EMELTI_loc
      public ERUNOSI_loc
      public ERUNPSI_loc
      public SRUNOSI_loc
      public SRUNPSI_loc
      public SMELTI_loc
      public DMUA_loc
      public grid_icdyn
      public pack_i2a
      public pack_a2i
      public DMUI_loc
      public DMVA_loc
      public DMVI_loc
      public SOLAR_loc
      public SSS_loc
      public UOSURF_loc
      public VOSURF_loc
      public UOSURF_4DYNSI_loc
      public VOSURF_4DYNSI_loc
      public OGEOZA_loc
      public GTEMP_loc
      public GTEMPR_loc
      public DMSI_loc
      public DHSI_loc
      public DSSI_loc

      public RSI_loc
      public FOCEAN_loc

      public alloc_hycom_atm,
     &     gather_atm_before_checkpoint, scatter_atm_after_checkpoint

      public im,jm

#ifdef USE_ATM_GLOBAL_ARRAYS
      public SSS
      public UOSURF
      public VOSURF
      public OGEOZA
      public GTEMP
      public GTEMPR
      public DMSI
      public DHSI
      public DSSI
      public FOCEAN
      public asst
      public atempr
#endif

      public ataux_loc,atauy_loc,aflxa2o_loc
     .     ,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc
     .     ,admui_loc,admvi_loc
     .     ,asst_loc,atempr_loc

#ifdef TRACERS_GASEXCH_ocean
      public atracflx_loc
      public GTRACER_glob, GTRACER_loc
      public TRGASEX_loc
#endif
#ifdef TRACERS_OceanBiology
      public asolz_loc
      public asolz
      public awind_loc
      public awind
      public wsavg, wsavg_loc
      public COSZ1, COSZ1_loc
#endif
#ifdef OBIO_RAD_coupling
      public avisdir_loc
      public avisdif_loc
      public anirdir_loc
      public anirdif_loc
      public avisdir
      public avisdif
      public anirdir
      public anirdif
      public FSRDIR, FSRDIR_loc
      public FSRDIF, FSRDIF_loc
      public DIRNIR, DIRNIR_loc
      public DIFNIR, DIFNIR_loc
      public SRVISSURF, SRVISSURF_loc
      public aCHL
#endif

#ifdef USE_ATM_GLOBAL_ARRAYS
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SSS
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: UOSURF
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: VOSURF
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OGEOZA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GTEMP
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GTEMPR
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DHSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DSSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FOCEAN
      ! arrays from cpl.h
      ! so far they should stay global since they are needed only 
      ! for coupling
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: asst
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: atempr
#endif

      ! accumulators for output on atmospheric grid
      ! (shouldn't these actually be accumulated on ocean grid
      !  and then interpolated to atmospheric grid?)
      real*8, allocatable, dimension(:,:) :: ataux_loc,atauy_loc
     .     ,aflxa2o_loc,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc
     .     ,admui_loc,admvi_loc ! == dmui_loc,dmvi_loc on atm. domain
     .     ,asst_loc,atempr_loc

#ifdef TRACERS_GASEXCH_ocean
      real, ALLOCATABLE, DIMENSION(:,:,:) :: atracflx_loc
      real, ALLOCATABLE, DIMENSION(:,:,:,:) :: GTRACER_glob
#endif
#ifdef TRACERS_OceanBiology
      real, ALLOCATABLE, DIMENSION(:,:)    :: asolz
      real, ALLOCATABLE, DIMENSION(:,:)    :: asolz_loc
      real, ALLOCATABLE, DIMENSION(:,:)    :: awind           !wind speed from modelE (see hycom2.f)
      real, ALLOCATABLE, DIMENSION(:,:)    :: awind_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: wsavg
      real, ALLOCATABLE, DIMENSION(:,:) :: COSZ1
#endif
#ifdef OBIO_RAD_coupling
      real*8, ALLOCATABLE, DIMENSION(:,:)    :: avisdir,avisdif
     .                                         ,anirdir,anirdif
      real*8, ALLOCATABLE, DIMENSION(:,:) :: avisdir_loc,avisdif_loc
     .     ,anirdir_loc,anirdif_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: FSRDIR,FSRDIF
      real, ALLOCATABLE, DIMENSION(:,:) :: DIRNIR,DIFNIR
      real, ALLOCATABLE, DIMENSION(:,:) :: SRVISSURF    
      real, ALLOCATABLE, DIMENSION(:,:) :: aCHL
#endif


      contains

      subroutine alloc_hycom_atm
      USE HYCOM_DIM, only : iia,jja ! actually should use IM,JM !!

#ifdef USE_ATM_GLOBAL_ARRAYS
      ALLOCATE( SSS( im, jm ) )
      ALLOCATE( UOSURF( im, jm ) )
      ALLOCATE( VOSURF( im, jm ) )
      ALLOCATE( OGEOZA( im, jm ) )
      ALLOCATE( GTEMP( 2 , NSTYPE, im, jm ) )
      ALLOCATE( GTEMPR( NSTYPE, im, jm ) )
      ALLOCATE( DMSI(  2  , im, jm ) )
      ALLOCATE( DHSI(  2  , im, jm ) )
      ALLOCATE( DSSI(  2  , im, jm ) )
      ALLOCATE( FOCEAN( im, jm ) )
      ALLOCATE( asst( im, jm ) )
      ALLOCATE( atempr( im, jm ) )
#endif

      ALLOCATE(
     &     ataux_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     atauy_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aflxa2o_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aemnp_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aice_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     asalt_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     austar_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aswflx_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     asst_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     atempr_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     admui_loc(im,aJ_0H:aJ_1H),
     &     admvi_loc(im,aJ_0H:aJ_1H)
     &     )

#ifdef TRACERS_GASEXCH_ocean
      ALLOCATE( atracflx_loc(iia,aJ_0H:aJ_1H,ntm) )
      ALLOCATE( GTRACER_glob ( NTM , NSTYPE , im , jm ) )
      GTRACER_glob = 0
#endif
#ifdef TRACERS_OceanBiology
      ALLOCATE(asolz(iia,jja))
      ALLOCATE(asolz_loc(iia,aJ_0H:aJ_1H))
      ALLOCATE(awind(iia,jja))
      ALLOCATE(awind_loc(iia,aJ_0H:aJ_1H))
      ALLOCATE( wsavg(im,jm) )
      ALLOCATE( COSZ1(im,jm) )
#endif
#ifdef OBIO_RAD_coupling
      ALLOCATE(avisdir(iia,jja),avisdif(iia,jja)
     .        ,anirdir(iia,jja),anirdif(iia,jja))
      ALLOCATE(avisdir_loc(iia,aJ_0H:aJ_1H),avisdif_loc(iia,aJ_0H:aJ_1H)
     &     ,anirdir_loc(iia,aJ_0H:aJ_1H),anirdif_loc(iia,aJ_0H:aJ_1H) )
      ALLOCATE(FSRDIR(im,jm) )
      ALLOCATE(FSRDIF(im,jm) )
      ALLOCATE(DIRNIR(im,jm) )
      ALLOCATE(DIFNIR(im,jm) )
      ALLOCATE(SRVISSURF(im,jm) )
      ALLOCATE(aCHL(im,jm) )
#endif

      end subroutine alloc_hycom_atm

cddd      subroutine alloc_locals
cddd
cddd      ALLOCATE( PREC_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EVAPOR_loc( I_0H:I_1H , J_0H:J_1H , NSTYPE ) )
cddd      ALLOCATE( FLOWO_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( GMELT_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( MELTI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( RUNOSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( RUNPSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( E0_loc( I_0H:I_1H , J_0H:J_1H , NSTYPE ) )
cddd      ALLOCATE( EPREC_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EFLOWO_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EGMELT_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( EMELTI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( ERUNOSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( ERUNPSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SRUNOSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SRUNPSI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SMELTI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DMUA_loc( I_0H:I_1H , J_0H:J_1H , NSTYPE ) )
cddd      ALLOCATE( DMUI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DMVA_loc( I_0H:I_1H , J_0H:J_1H, NSTYPE ) )
cddd      ALLOCATE( DMVI_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SOLAR_loc(  3  , I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( SSS_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( UOSURF_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( VOSURF_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( OGEOZA_loc( I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( GTEMP_loc( 2 , NSTYPE, I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( GTEMPR_loc( NSTYPE, I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DMSI_loc(  2  , I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DHSI_loc(  2  , I_0H:I_1H , J_0H:J_1H ) )
cddd      ALLOCATE( DSSI_loc(  2  , I_0H:I_1H , J_0H:J_1H ) )
cddd
cddd
cddd      end subroutine alloc_locals


      subroutine gather_atm_before_checkpoint
#ifdef USE_ATM_GLOBAL_ARRAYS
      USE DOMAIN_DECOMP_1D, ONLY: GRID,PACK_DATA,PACK_COLUMN,PACK_BLOCK

      call pack_data( grid,  ASST_loc, ASST )
      call pack_data( grid,  ATEMPR_loc, ATEMPR )
      call pack_data( grid,  SSS_loc, SSS )
      call pack_data( grid,  UOSURF_loc, UOSURF )
      call pack_data( grid,  VOSURF_loc, VOSURF )
      call pack_data( grid,  OGEOZA_loc, OGEOZA )
       call pack_column( grid,  DMSI_loc, DMSI )
       call pack_column( grid,  DHSI_loc, DHSI )
       call pack_column( grid,  DSSI_loc, DSSI )

c these arrays are not referenced by io_ocean.  why gather?

#ifdef TRACERS_OceanBiology
      call pack_data(  grid,wsavg_loc,wsavg)
      call pack_data(  grid,COSZ1_loc,COSZ1)
#endif
#ifdef OBIO_RAD_coupling
      call pack_data(grid,FSRDIR_loc,FSRDIR)
      call pack_data(grid,FSRDIF_loc,FSRDIF)
      call pack_data(grid,DIRNIR_loc,DIRNIR)
      call pack_data(grid,DIFNIR_loc,DIFNIR)
      call pack_data(grid,SRVISSURF_loc,SRVISSURF)
#endif

#endif /* USE_ATM_GLOBAL_ARRAYS */
      end subroutine gather_atm_before_checkpoint

      
      subroutine scatter_atm_after_checkpoint
#ifdef USE_ATM_GLOBAL_ARRAYS
      USE DOMAIN_DECOMP_1D, ONLY: GRID, UNPACK_DATA, UNPACK_COLUMN,
     &     UNPACK_BLOCK

      call unpack_data( grid,  ASST, ASST_loc )
      call unpack_data( grid,  ATEMPR, ATEMPR_loc )
      call unpack_data( grid,  SSS, SSS_loc )
      call unpack_data( grid,  UOSURF, UOSURF_loc )
      call unpack_data( grid,  VOSURF, VOSURF_loc )
c UOSURF and VOSURF are also needed on the ice dynamics A-grid.
c For the moment, HYCOM only runs with modelE configurations having
c identical atmosphere and ice dynamics grids, so the atmospheric
c copy of UOSURF,VOSURF can be used.
      if(grid_icdyn%have_domain) then ! ice dyn may run on subset of PEs
        call unpack_data( grid_icdyn,  UOSURF, UOSURF_4DYNSI_loc)
        call unpack_data( grid_icdyn,  VOSURF, VOSURF_4DYNSI_loc) 
      endif
      call unpack_data( grid,  OGEOZA, OGEOZA_loc )
       call unpack_column( grid,  DMSI, DMSI_loc )
       call unpack_column( grid,  DHSI, DHSI_loc )
       call unpack_column( grid,  DSSI, DSSI_loc )

#endif /* USE_ATM_GLOBAL_ARRAYS */
      end subroutine scatter_atm_after_checkpoint

      end module hycom_atm
