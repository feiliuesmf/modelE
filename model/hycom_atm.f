#include "rundeck_opts.h"

      module hycom_atm
!@sum module for atmospheric variables to be passed to/from hylom.
!@+   hycom will see them as global arrays
      USE HYCOM_DIM, only : aI_0, aI_1, aI_0H, aI_1H
      USE HYCOM_DIM, only : aJ_0, aJ_1, aJ_0H, aJ_1H
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      implicit none

      private

      type(atmocn_xchng_vars), public :: ocnatm

      public alloc_hycom_atm

      public ataux_loc,atauy_loc,aflxa2o_loc
     .     ,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc
     .     ,admui_loc,admvi_loc

#ifdef TRACERS_GASEXCH_ocean
      public atracflx_loc
#endif
#ifdef TRACERS_OceanBiology
      public asolz_loc
      public awind_loc
#endif
#ifdef OBIO_RAD_coupling
      public avisdir_loc
      public avisdif_loc
      public anirdir_loc
      public anirdif_loc
#endif

      ! accumulators for output on atmospheric grid
      ! (shouldn't these actually be accumulated on ocean grid
      !  and then interpolated to atmospheric grid?)
      real*8, allocatable, dimension(:,:) :: ataux_loc,atauy_loc
     .     ,aflxa2o_loc,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc
     .     ,admui_loc,admvi_loc ! == dmui_loc,dmvi_loc on atm. domain

#ifdef TRACERS_GASEXCH_ocean
      real, ALLOCATABLE, DIMENSION(:,:,:) :: atracflx_loc
#endif
#ifdef TRACERS_OceanBiology
      real, ALLOCATABLE, DIMENSION(:,:)    :: asolz_loc
!wind speed from modelE (see hycom2.f)
      real, ALLOCATABLE, DIMENSION(:,:)    :: awind_loc
#endif
#ifdef OBIO_RAD_coupling
      real*8, ALLOCATABLE, DIMENSION(:,:) :: avisdir_loc,avisdif_loc
     .     ,anirdir_loc,anirdif_loc
#endif


      contains

      subroutine alloc_hycom_atm(atmocn)
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      USE EXCHANGE_TYPES, only : alloc_xchng_vars
      USE HYCOM_DIM, only : ogrid

      type(atmocn_xchng_vars) :: atmocn

#ifdef TRACERS_GASEXCH_ocean
      ocnatm % ntm = atmocn % ntm
      ocnatm % ntm_gasexch = atmocn % ntm_gasexch
#endif
      call alloc_xchng_vars(ogrid,ocnatm)

      ALLOCATE(
     &     ataux_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     atauy_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aflxa2o_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aemnp_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aice_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     asalt_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     austar_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aswflx_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     admui_loc(aI_0H:aI_1H,aJ_0H:aJ_1H), ! temporary
     &     admvi_loc(aI_0H:aI_1H,aJ_0H:aJ_1H)  ! temporary
     &     )

#ifdef TRACERS_GASEXCH_ocean
      ALLOCATE(atracflx_loc(aI_0H:aI_1H,aJ_0H:aJ_1H,atmocn%ntm_gasexch))
#endif
#ifdef TRACERS_OceanBiology
      ALLOCATE(asolz_loc(aI_0H:aI_1H,aJ_0H:aJ_1H))
      ALLOCATE(awind_loc(aI_0H:aI_1H,aJ_0H:aJ_1H))
#endif
#ifdef OBIO_RAD_coupling
      ALLOCATE(avisdir_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &         avisdif_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &         anirdir_loc(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &         anirdif_loc(aI_0H:aI_1H,aJ_0H:aJ_1H) )
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

      end module hycom_atm
