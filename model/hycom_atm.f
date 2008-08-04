#include "rundeck_opts.h"

      module hycom_atm
!@sum module for atmospheric variables to be passed to/from hylom.
!@+   hycom will see them as global arrays
      USE MODEL_COM, only : im,jm
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
      USE FLUXES, only : DMUI_loc => DMUI
      USE FLUXES, only : DMVA_loc => DMVA
      USE FLUXES, only : DMVI_loc => DMVI
      USE FLUXES, only : SOLAR_loc => SOLAR
      USE FLUXES, only : SSS_loc => SSS
      USE FLUXES, only : UOSURF_loc => UOSURF
      USE FLUXES, only : VOSURF_loc => VOSURF
      USE FLUXES, only : OGEOZA_loc => OGEOZA
      USE FLUXES, only : GTEMP_loc => GTEMP
      USE FLUXES, only : GTEMPR_loc => GTEMPR
      USE FLUXES, only : DMSI_loc => DMSI
      USE FLUXES, only : DHSI_loc => DHSI
      USE FLUXES, only : DSSI_loc => DSSI

      USE SEAICE_COM, only : RSI_loc => RSI ! seems to be used for diags only?
      USE MODEL_COM, only : FOCEAN_loc => FOCEAN

#ifdef TRACERS_GASEXCH_Natassa
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
      public DMUI_loc
      public DMVA_loc
      public DMVI_loc
      public SOLAR_loc
      public SSS_loc
      public UOSURF_loc
      public VOSURF_loc
      public OGEOZA_loc
      public GTEMP_loc
      public GTEMPR_loc
      public DMSI_loc
      public DHSI_loc
      public DSSI_loc

      public RSI_loc
      public FOCEAN_loc

      public alloc_hycom_atm, gather_atm, scatter_atm

      public im,jm

      public PREC
      public EVAPOR
      public FLOWO
      public GMELT
      public MELTI
      public RUNOSI
      public RUNPSI
      public E0
      public EPREC
      public EFLOWO
      public EGMELT
      public EMELTI
      public ERUNOSI
      public SRUNOSI
      public SRUNPSI
      public ERUNPSI
      public SMELTI
      public DMUA
      public DMUI
      public DMVA
      public DMVI
      public SOLAR
      public SSS
      public UOSURF
      public VOSURF
      public OGEOZA
      public GTEMP
      public GTEMPR
      public DMSI
      public DHSI
      public DSSI

      public RSI
      public FOCEAN

      public asst
      public atempr

      public ataux,atauy,aflxa2o
     .     ,aemnp,aice,asalt
     .     ,austar,aswflx
      public ataux_loc,atauy_loc,aflxa2o_loc
     .     ,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc

#ifdef TRACERS_GASEXCH_Natassa
      public GTRACER, GTRACER_loc
      public TRGASEX, TRGASEX_loc
#endif
#ifdef TRACERS_OceanBiology
      public wsavg, wsavg_loc
      public COSZ1, COSZ1_loc
#endif
#ifdef OBIO_RAD_coupling
      public FSRDIR, FSRDIR_loc
      public FSRDIF, FSRDIF_loc
      public DIRNIR, DIRNIR_loc
      public DIFNIR, DIFNIR_loc
      public SRVISSURF, SRVISSURF_loc
#endif

      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PREC
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: EVAPOR
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FLOWO
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: GMELT
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MELTI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNOSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RUNPSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: E0
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EPREC
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EFLOWO
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EGMELT
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EMELTI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ERUNOSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ERUNPSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRUNOSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRUNPSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SMELTI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMUA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DMUI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMVA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DMVI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SOLAR
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SSS
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: UOSURF
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: VOSURF
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OGEOZA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GTEMP
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GTEMPR
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DHSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DSSI

      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: FOCEAN

      ! arrays from cpl.h
      ! so far they should stay global since they are needed only 
      ! for coupling
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: asst
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: atempr

      ! accumulators for output on atmospheric grid
      ! (shouldn't these actually be accumulated on ocean grid
      !  and then interpolated to atmospheric grid?)
      real*8, allocatable, dimension(:,:) :: ataux,atauy,aflxa2o
     .     ,aemnp,aice,asalt
     .     ,austar,aswflx
      real*8, allocatable, dimension(:,:) :: ataux_loc,atauy_loc
     .     ,aflxa2o_loc,aemnp_loc,aice_loc,asalt_loc
     .     ,austar_loc,aswflx_loc

#ifdef TRACERS_GASEXCH_Natassa
      real, ALLOCATABLE, DIMENSION(:,:,:,:) :: GTRACER
      real, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRGASEX
#endif
#ifdef TRACERS_OceanBiology
      real, ALLOCATABLE, DIMENSION(:,:) :: wsavg
      real, ALLOCATABLE, DIMENSION(:,:) :: COSZ1

#endif
#ifdef OBIO_RAD_coupling
      real, ALLOCATABLE, DIMENSION(:,:) :: FSRDIR,FSRDIF
      real, ALLOCATABLE, DIMENSION(:,:) :: DIRNIR,DIFNIR
      real, ALLOCATABLE, DIMENSION(:,:) :: SRVISSURF    
#endif


      contains

      subroutine alloc_hycom_atm
      USE HYCOM_DIM, only : iia,jja ! actually should use IM,JM !!

      ALLOCATE( PREC( im, jm ) )
      ALLOCATE( EVAPOR( im, jm , NSTYPE ) )
      ALLOCATE( FLOWO( im, jm ) )
      ALLOCATE( GMELT( im, jm ) )
      ALLOCATE( MELTI( im, jm ) )
      ALLOCATE( RUNOSI( im, jm ) )
      ALLOCATE( RUNPSI( im, jm ) )
      ALLOCATE( E0( im, jm , NSTYPE ) )
      ALLOCATE( EPREC( im, jm ) )
      ALLOCATE( EFLOWO( im, jm ) )
      ALLOCATE( EGMELT( im, jm ) )
      ALLOCATE( EMELTI( im, jm ) )
      ALLOCATE( ERUNOSI( im, jm ) )
      ALLOCATE( ERUNPSI( im, jm ) )
      ALLOCATE( SRUNOSI( im, jm ) )
      ALLOCATE( SRUNPSI( im, jm ) )
      ALLOCATE( SMELTI( im, jm ) )
      ALLOCATE( DMUA( im, jm , NSTYPE ) )
      ALLOCATE( DMUI( im, jm ) )
      ALLOCATE( DMVA( im, jm, NSTYPE ) )
      ALLOCATE( DMVI( im, jm ) )
      ALLOCATE( SOLAR(  3  , im, jm ) )
      ALLOCATE( SSS( im, jm ) )
      ALLOCATE( UOSURF( im, jm ) )
      ALLOCATE( VOSURF( im, jm ) )
      ALLOCATE( OGEOZA( im, jm ) )
      ALLOCATE( GTEMP( 2 , NSTYPE, im, jm ) )
      ALLOCATE( GTEMPR( NSTYPE, im, jm ) )
      ALLOCATE( DMSI(  2  , im, jm ) )
      ALLOCATE( DHSI(  2  , im, jm ) )
      ALLOCATE( DSSI(  2  , im, jm ) )

      ALLOCATE( RSI( im, jm ) )
      ALLOCATE( FOCEAN( im, jm ) )

      ALLOCATE( asst( im, jm ) )
      ALLOCATE( atempr( im, jm ) )

      ALLOCATE( ataux(iia,jja),atauy(iia,jja),aflxa2o(iia,jja)
     .     ,aemnp(iia,jja),aice(iia,jja),asalt(iia,jja)
     .     ,austar(iia,jja),aswflx(iia,jja) )
      ALLOCATE( ataux_loc(iia,aJ_0H:aJ_1H),atauy_loc(iia,aJ_0H:aJ_1H)
     .     ,aflxa2o_loc(iia,aJ_0H:aJ_1H)
     .     ,aemnp_loc(iia,aJ_0H:aJ_1H),aice_loc(iia,aJ_0H:aJ_1H),
     .      asalt_loc(iia,aJ_0H:aJ_1H)
     .     ,austar_loc(iia,aJ_0H:aJ_1H),aswflx_loc(iia,aJ_0H:aJ_1H) )

#ifdef TRACERS_GASEXCH_Natassa
      ALLOCATE( GTRACER ( NTM , NSTYPE , im , jm ) )
      GTRACER = 0
      ALLOCATE( TRGASEX(NTM , NSTYPE , im , jm ) )
#endif
#ifdef TRACERS_OceanBiology
      ALLOCATE( wsavg(im,jm) )
      ALLOCATE( COSZ1(im,jm) )
#endif
#ifdef OBIO_RAD_coupling
      ALLOCATE(FSRDIR(im,jm) )
      ALLOCATE(FSRDIF(im,jm) )
      ALLOCATE(DIRNIR(im,jm) )
      ALLOCATE(DIFNIR(im,jm) )
      ALLOCATE(SRVISSURF(im,jm) )
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


!!!!!!! probably over-kill: need to check what actually needs 
!!!!!!! packing and what needs unpacking

      subroutine gather_atm
      USE DOMAIN_DECOMP, ONLY: GRID, PACK_DATA, PACK_COLUMN, PACK_BLOCK

      call pack_data( grid,  PREC_loc, PREC )
      call pack_data( grid,  EVAPOR_loc, EVAPOR )
      call pack_data( grid,  FLOWO_loc, FLOWO )
      call pack_data( grid,  GMELT_loc, GMELT )
      call pack_data( grid,  MELTI_loc, MELTI )
      call pack_data( grid,  RUNOSI_loc, RUNOSI )
      call pack_data( grid,  RUNPSI_loc, RUNPSI )
      call pack_data( grid,  E0_loc, E0 )
      call pack_data( grid,  EPREC_loc, EPREC )
      call pack_data( grid,  EFLOWO_loc, EFLOWO )
      call pack_data( grid,  EGMELT_loc, EGMELT )
      call pack_data( grid,  EMELTI_loc, EMELTI )
      call pack_data( grid,  ERUNOSI_loc, ERUNOSI )
      call pack_data( grid,  ERUNPSI_loc, ERUNPSI )
      call pack_data( grid,  SRUNOSI_loc, SRUNOSI )
      call pack_data( grid,  SRUNPSI_loc, SRUNPSI )
      call pack_data( grid,  SMELTI_loc, SMELTI )
      call pack_data( grid,  DMUA_loc, DMUA )
      call pack_data( grid,  DMUI_loc, DMUI )
      call pack_data( grid,  DMVA_loc, DMVA )
      call pack_data( grid,  DMVI_loc, DMVI )
       call pack_column( grid,  SOLAR_loc, SOLAR )
      call pack_data( grid,  SSS_loc, SSS )
      call pack_data( grid,  UOSURF_loc, UOSURF )
      call pack_data( grid,  VOSURF_loc, VOSURF )
      call pack_data( grid,  OGEOZA_loc, OGEOZA )
       call pack_block( grid,  GTEMP_loc, GTEMP )
       call pack_column( grid,  GTEMPR_loc, GTEMPR )
       call pack_column( grid,  DMSI_loc, DMSI )
       call pack_column( grid,  DHSI_loc, DHSI )
       call pack_column( grid,  DSSI_loc, DSSI )

#ifdef TRACERS_GASEXCH_Natassa
      call pack_block( grid,GTRACER_loc,GTRACER)
      call pack_block( grid,TRGASEX_loc,TRGASEX)
#endif
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

      ! these are not changed by hycom - no need to scatter
      call pack_data( grid,  RSI_loc, RSI )
      call pack_data( grid,  FOCEAN_loc, FOCEAN )


      end subroutine gather_atm

      
      subroutine scatter_atm
      USE DOMAIN_DECOMP, ONLY: GRID, UNPACK_DATA, UNPACK_COLUMN,
     &     UNPACK_BLOCK

      call unpack_data( grid,  PREC, PREC_loc )
      call unpack_data( grid,  EVAPOR, EVAPOR_loc )
      call unpack_data( grid,  FLOWO, FLOWO_loc )
      call unpack_data( grid,  GMELT, GMELT_loc )
      call unpack_data( grid,  MELTI, MELTI_loc )
      call unpack_data( grid,  RUNOSI, RUNOSI_loc )
      call unpack_data( grid,  RUNPSI, RUNPSI_loc )
      call unpack_data( grid,  E0, E0_loc )
      call unpack_data( grid,  EPREC, EPREC_loc )
      call unpack_data( grid,  EFLOWO, EFLOWO_loc )
      call unpack_data( grid,  EGMELT, EGMELT_loc )
      call unpack_data( grid,  EMELTI, EMELTI_loc )
      call unpack_data( grid,  ERUNOSI, ERUNOSI_loc )
      call unpack_data( grid,  ERUNPSI, ERUNPSI_loc )
      call unpack_data( grid,  SRUNOSI, SRUNOSI_loc )
      call unpack_data( grid,  SRUNPSI, SRUNPSI_loc )
      call unpack_data( grid,  SMELTI, SMELTI_loc )
      call unpack_data( grid,  DMUA, DMUA_loc )
      call unpack_data( grid,  DMUI, DMUI_loc )
      call unpack_data( grid,  DMVA, DMVA_loc )
      call unpack_data( grid,  DMVI, DMVI_loc )
       call unpack_column( grid,  SOLAR, SOLAR_loc )
      call unpack_data( grid,  SSS, SSS_loc )
      call unpack_data( grid,  UOSURF, UOSURF_loc )
      call unpack_data( grid,  VOSURF, VOSURF_loc )
      call unpack_data( grid,  OGEOZA, OGEOZA_loc )
       call unpack_block( grid,  GTEMP, GTEMP_loc )
       call unpack_column( grid,  GTEMPR, GTEMPR_loc )
       call unpack_column( grid,  DMSI, DMSI_loc )
       call unpack_column( grid,  DHSI, DHSI_loc )
       call unpack_column( grid,  DSSI, DSSI_loc )

#ifdef TRACERS_GASEXCH_Natassa
      call unpack_block( grid,GTRACER,GTRACER_loc)
#endif

      end subroutine scatter_atm

      end module hycom_atm
