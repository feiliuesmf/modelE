C****   
C**** OCEAN_COM.f    Module Variables for Ocean    2007/01/12
C****
#include "rundeck_opts.h"

      Module OCEAN
!@sum  OCEAN dynamic ocean related variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
C**** Note that we currently use the same horizontal grid as for the
C**** atmosphere. However, we can redefine im,jm if necessary.
      Use CONSTANT,  Only: TWOPI
      Use OCEANRES,  Only: IM=>IMO,JM=>JMO, LMO, LMO_MIN, LSRPD, dZO
#ifdef TRACERS_OCEAN
      Use OCN_TRACER_COM, Only : ntm
#endif
      Use SparseCommunicator_mod
#if defined(CUBED_SPHERE) || defined(CUBE_GRID)
      use cs2ll_utils, only : aoremap_type=>xgridremap_type
#else
      use hntrp_mod, only :   aoremap_type=>hntrp_type
#endif
      Implicit None
      Integer*4,Parameter ::
     *  IVSP = 3*IM/4,      !  V at south pole is stored in U(IVSP,1)
     *  IVNP =   IM/4       !  V at north pole is stored in U(IVNP,JM)
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM

!@dbparam OBottom_drag use ocean bottom drag routine (default=1)
!@dbparam OCoastal_drag use ocean coastal drag routine (default=1)
      Integer*4 ::
     *     OBottom_drag = 1,    
     *     OCoastal_drag = 1

!@dbparam oc_mean_salinity define mean salinity of ocean (if set)
!@dbparam oc_tracer_mean mean tracer ratio of ocean 
!@+       (in permil for water isotopes)
      REAL*8 :: 
     *     oc_salt_mean = -999. 
#ifdef TRACERS_OCEAN
     *     , oc_tracer_mean(ntm) = -999. 
#endif 
!@var MO mass of ocean (kg/m^2)
!@var UO E-W velocity on C-grid (m/s)
!@var VO N-S velocity on C-grid (m/s)
!@var UOD E-W velocity on D-grid (m/s)
!@var VOD N-S velocity on D-grid (m/s)
!@var G0M,GXMO,GYMO,GZMO pot. enthalpy of ocean (+moments) (J)
!@var S0M,SXMO,SYMO,SZMO salinity of ocean (+moments) (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: MO,UO,VO,UOD,VOD,
     *     G0M,GXMO,GYMO,GZMO, S0M,SXMO,SYMO,SZMO
C**** Global arrays needed for i/o, GM,straits,odiff ?
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: 
     *     MO_glob,UO_glob,VO_glob,UOD_glob,VOD_glob,
     *     G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob,
     *     S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob

      Real*8 UONP(LMO), !  U component at north pole, points down 90W
     *       VONP(LMO)  !  V component at north pole, points down 0 (GM)

C**** ocean geometry (should this be in a separate module?)
      Real*8 
     *     DLON,    !@var DLON longitudinal spacing in radians
     *     DLAT,    !@var DLAT latitudinal spacing in radians
     *     DLATM,   !@var DLATM latitudinal spacing in minutes
     *     FJEQ     !@var FJEQ location of equator in grid units 
     *   , oDLAT_DG ! grid spacing in latitude (deg) 
     *   , oDLON_DG ! grid spacing in longitude (deg) 

      REAL*8, ALLOCATABLE, DIMENSION(:,:):: OXYP

      REAL*8, DIMENSION(JM) :: DXYPO,DXPO,DYPO,DXVO,DYVO
     *     ,COSPO,SINPO,DXYVO,DXYSO,DXYNO,RAMVS,RAMVN,RLAT,BYDXYPO
      REAL*8, DIMENSION(0:JM) :: SINVO,COSVO
      REAL*8, DIMENSION(0:LMO) :: ZE
      Real*8
     *  DXPGF(0:JM),! DXYV/dYPGF is north-south distance used in PGF
     *  DYPGF(JM), !  DXYP/dXPGF is east-west distance used in PGF
     *  COSM(JM),  !  = .5*COSV(J-1) + .5*COSV(J)
     *  COSQ(JM),  !  sQuare of COSine = .5*COSV(J-1)^2 + .5*COSV(J)^2
     *  SINIC(IM), !  SINe of longitude of grid cell center from IDL
     *  COSIC(IM), !  COSine of longitude of grid cell center from IDL
     *  SINU(IM),  !  SINe of longitude of eastern edge of grid cell
     *  COSU(IM),  !  COSine of longitude of eastern edge of grid cell
     *  SINxY(JM), !  DLAT * SINe of latitude used by Coriolis force
     *  TANxY(JM)  !  DLAT * TANgent of latitude used by metric term
!@var  oLAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JM,2) :: oLAT_DG
!@var  oLON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, DIMENSION(IM,2) :: oLON_DG
!@var  IMAXJ varying number of used longitudes
      INTEGER, DIMENSION(JM) :: IMAXJ
      INTEGER, DIMENSION(IM,JM) :: LMM,LMU,LMV
      Integer*4 J40S, !  maximum grid cell below 40S (used in OPFIL)
     *           J1O  !  most southern latitude (J) where ocean exists
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: oLAT2D_DG !distributed latitute array (in degrees)
      REAL*8, DIMENSION(IM,JM) :: HATMO,HOCEAN,FOCEAN
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: FOCEAN_loc
C**** ocean related parameters
      INTEGER NDYNO,MDYNO,MSGSO
!@dbparam DTO timestep for ocean dynamics (s)
      REAL*8 :: DTO=450.        ! default. setable parameter
      REAL*8 DTOFS,DTOLF,DTS,BYDTS
!@var budget grid quantities (defined locally on each proc.)
      REAL*8, ALLOCATABLE, DIMENSION(:,:):: owtbudg
      INTEGER, ALLOCATABLE, DIMENSION(:,:):: oJ_BUDG
      INTEGER :: oJ_0B,oJ_1B
!@var OPRESS Anomalous pressure at surface of ocean (under ice) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OPRESS ! (IM,JM)
      REAL*8, DIMENSION(IM,JM) :: OPRESS_glob ! needed for straits ?
!@var OGEOZ ocean geopotential at surface (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OGEOZ,OGEOZ_SV
      REAL*8, DIMENSION(IM,JM) :: OGEOZ_glob, OGEOZ_SV_glob ! for i/o
!@var OPBOT ocean bottom pressure (diagnostic only) (Pa)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: OPBOT

#ifdef TRACERS_OCEAN
!@var TRMO,TXMO,TYMO,TZMO tracer amount (+moments) in ocean (kg)
! for i/o, straits, GM ?
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: 
     *       TRMO_glob,TXMO_glob,TYMO_glob,TZMO_glob
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRMO,TXMO,TYMO,TZMO
#endif
     
      type (SparseCommunicator_type), save :: mySparseComm_type

!@var nbyz[muvc]: # of basins at each lat/depth
!@var i[12]yz[muvc]: start/end i-indices for each basin
! m: cell center ! u: east edge ! v: north edge ! c: northeast corner
      integer, parameter :: nbyzmax=20 ! suffices up to 1x1.25 deg res
      integer, dimension(:,:), allocatable :: nbyzm,nbyzu,nbyzv,nbyzc
      integer, dimension(:,:,:), allocatable ::
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc

!@var remap_a2o,remap_o2a atm->ocn,ocn->atm interpolation info
      type(aoremap_type) ::
     &     remap_a2o            ! atm A -> ocn A
     &    ,remap_o2a            ! ocn A -> atm A

      REAL*8, SAVE, ALLOCATABLE, DIMENSION(:) ::
     *     BYDXYV, KHP,KHV,TANP,TANV,BYDXV,BYDXP,BYDYV,BYDYP
      REAL*8, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::
     *      UXA,UXB,UXC,UYA,UYB,UYC,VXA,VXB,VXC,VYA,VYB,VYC
      REAL*8, SAVE :: BYDXYPJM
      REAL*8, SAVE, DIMENSION(LMO) :: UYPB
      REAL*8, SAVE, DIMENSION(IM,LMO) :: UYPA
      REAL*8, PARAMETER :: FSLIP=0.

      contains

      subroutine gather_ocean (icase)

      use domain_decomp_1d, only: pack_data
      use OCEANR_DIM, only : grid=>ogrid

      integer, intent(in) :: icase

c**** icase=0:full i/o, 1:ini_straits, 2:serialized ocn dynamics
      CALL PACK_DATA(grid,   MO   ,    MO_glob)
      if (icase.ne.1) then      ! needed for i/o and ODIFF only
        CALL PACK_DATA(grid,   UO   ,    UO_glob)
        CALL PACK_DATA(grid,   VO   ,    VO_glob)
        CALL PACK_DATA(grid,   UOD   ,    UOD_glob)
        CALL PACK_DATA(grid,   VOD   ,    VOD_glob)
      end if
      if (icase.lt.1) then      ! needed for i/o only
        CALL PACK_DATA(grid,OGEOZ   , OGEOZ_glob)
        CALL PACK_DATA(grid,OGEOZ_SV,OGEOZ_SV_glob)
      end if

      CALL PACK_DATA(grid,  G0M   ,   G0M_glob)
      CALL PACK_DATA(grid,  GXMO  ,  GXMO_glob)
      CALL PACK_DATA(grid,  GYMO  ,  GYMO_glob)
      CALL PACK_DATA(grid,  GZMO  ,  GZMO_glob)
      CALL PACK_DATA(grid,  S0M   ,   S0M_glob)
      CALL PACK_DATA(grid,  SXMO  ,  SXMO_glob)
      CALL PACK_DATA(grid,  SYMO  ,  SYMO_glob)
      CALL PACK_DATA(grid,  SZMO  ,  SZMO_glob)
#ifdef TRACERS_OCEAN
      CALL PACK_DATA(grid,  TRMO  ,  TRMO_glob)
      CALL PACK_DATA(grid,  TXMO  ,  TXMO_glob)
      CALL PACK_DATA(grid,  TYMO  ,  TYMO_glob)
      CALL PACK_DATA(grid,  TZMO  ,  TZMO_glob)
#endif
      if (icase.lt.2) return

c**** icase=2: still serialized non-i/o parts of ocn dynamics
               ! for straits:  mo,G0M,...,S0M,...,TRMO,...,opress
      CALL PACK_DATA(grid,  OPRESS,OPRESS_glob)
               ! for OCNGM:    mo,G0M,...,S0M,...,TRMO,...
               ! for ODIFF:    mo,uo,vo

      RETURN
      end subroutine gather_ocean

      subroutine scatter_ocean (icase)

      use domain_decomp_1d, only: unpack_data
      use OCEANR_DIM, only : grid=>ogrid

      integer, intent(in) :: icase

c**** icase=-1: i/o no_trc 0:full i/o, 1:ini_straits, 2:serial ocn dyn
      CALL UNPACK_DATA(grid,       MO_glob,   MO )
      if (icase.lt.1) then            ! needed for i/o only
        CALL UNPACK_DATA(grid,       UO_glob,   UO )
        CALL UNPACK_DATA(grid,       VO_glob,   VO )
        CALL UNPACK_DATA(grid,       UOD_glob,   UOD )
        CALL UNPACK_DATA(grid,       VOD_glob,   VOD )
        CALL UNPACK_DATA(grid,    OGEOZ_glob, OGEOZ   )
        CALL UNPACK_DATA(grid, OGEOZ_SV_glob, OGEOZ_SV)
      end if

      CALL UNPACK_DATA(grid,      G0M_glob,  G0M )
      CALL UNPACK_DATA(grid,     GXMO_glob,  GXMO)
      CALL UNPACK_DATA(grid,     GYMO_glob,  GYMO)
      CALL UNPACK_DATA(grid,     GZMO_glob,  GZMO)
      CALL UNPACK_DATA(grid,      S0M_glob,  S0M )
      CALL UNPACK_DATA(grid,     SXMO_glob,  SXMO)
      CALL UNPACK_DATA(grid,     SYMO_glob,  SYMO)
      CALL UNPACK_DATA(grid,     SZMO_glob,  SZMO)
#ifdef TRACERS_OCEAN
      if (icase.lt.0) return                   ! IC w/o tracers
      CALL UNPACK_DATA(grid,     TRMO_glob,  TRMO)
      CALL UNPACK_DATA(grid,     TXMO_glob,  TXMO)
      CALL UNPACK_DATA(grid,     TYMO_glob,  TYMO)
      CALL UNPACK_DATA(grid,     TZMO_glob,  TZMO)
#endif
      if (icase.lt.2) return

c**** icase=2: still serialized non-i/o parts of ocn dynamics
               ! for straits: mo,G0M,...,S0M,...,TRMO,...,opress
      CALL UNPACK_DATA(grid,   OPRESS_glob,OPRESS)

      RETURN
      end subroutine scatter_ocean

      subroutine gather_straits_to_global (icase)
      use SparseCommunicator_mod, only: gatherIJ
      integer, intent(in) :: icase

c**** icase=0:full i/o, 1:ini_straits, 2:serialized ocn dynamics
      CALL gatherIJ(mySparseComm_type,   MO   ,    MO_glob)
      if (icase.ne.1) then      ! needed for i/o and ODIFF only
        CALL gatherIJ(mySparseComm_type,   UO   ,    UO_glob)
        CALL gatherIJ(mySparseComm_type,   VO   ,    VO_glob)
      end if
      if (icase.lt.1) then      ! needed for i/o only
        CALL gatherIJ(mySparseComm_type,OGEOZ   , OGEOZ_glob)
        CALL gatherIJ(mySparseComm_type,OGEOZ_SV,OGEOZ_SV_glob)
      end if

      CALL gatherIJ(mySparseComm_type,  G0M   ,   G0M_glob)
      CALL gatherIJ(mySparseComm_type,  GXMO  ,  GXMO_glob)
      CALL gatherIJ(mySparseComm_type,  GYMO  ,  GYMO_glob)
      CALL gatherIJ(mySparseComm_type,  GZMO  ,  GZMO_glob)
      CALL gatherIJ(mySparseComm_type,  S0M   ,   S0M_glob)
      CALL gatherIJ(mySparseComm_type,  SXMO  ,  SXMO_glob)
      CALL gatherIJ(mySparseComm_type,  SYMO  ,  SYMO_glob)
      CALL gatherIJ(mySparseComm_type,  SZMO  ,  SZMO_glob)
#ifdef TRACERS_OCEAN
      CALL gatherIJ(mySparseComm_type,  TRMO  ,  TRMO_glob)
      CALL gatherIJ(mySparseComm_type,  TXMO  ,  TXMO_glob)
      CALL gatherIJ(mySparseComm_type,  TYMO  ,  TYMO_glob)
      CALL gatherIJ(mySparseComm_type,  TZMO  ,  TZMO_glob)
#endif
      if (icase.lt.2) return

c**** icase=2: still serialized non-i/o parts of ocn dynamics
               ! for straits:  mo,G0M,...,S0M,...,TRMO,...,opress
      CALL gatherIJ(mySparseComm_type,  OPRESS,OPRESS_glob)
               ! for OCNGM:    mo,G0M,...,S0M,...,TRMO,...
               ! for ODIFF:    mo,uo,vo

      RETURN
      end subroutine gather_straits_to_global

      subroutine scatter_straits_from_global (icase)
      use SparseCommunicator_mod, only: scatterIJ
      integer, intent(in) :: icase

c**** icase=-1: i/o no_trc 0:full i/o, 1:ini_straits, 2:serial ocn dyn
      CALL scatterIJ(mySparseComm_type,       MO_glob,   MO )
      if (icase.lt.1) then            ! needed for i/o only
        CALL scatterIJ(mySparseComm_type,       UO_glob,   UO )
        CALL scatterIJ(mySparseComm_type,       VO_glob,   VO )
        CALL scatterIJ(mySparseComm_type,    OGEOZ_glob, OGEOZ   )
        CALL scatterIJ(mySparseComm_type, OGEOZ_SV_glob, OGEOZ_SV)
      end if

      CALL scatterIJ(mySparseComm_type,      G0M_glob,  G0M )
      CALL scatterIJ(mySparseComm_type,     GXMO_glob,  GXMO)
      CALL scatterIJ(mySparseComm_type,     GYMO_glob,  GYMO)
      CALL scatterIJ(mySparseComm_type,     GZMO_glob,  GZMO)
      CALL scatterIJ(mySparseComm_type,      S0M_glob,  S0M )
      CALL scatterIJ(mySparseComm_type,     SXMO_glob,  SXMO)
      CALL scatterIJ(mySparseComm_type,     SYMO_glob,  SYMO)
      CALL scatterIJ(mySparseComm_type,     SZMO_glob,  SZMO)
#ifdef TRACERS_OCEAN
      if (icase.lt.0) return                   ! IC w/o tracers
      CALL scatterIJ(mySparseComm_type,     TRMO_glob,  TRMO)
      CALL scatterIJ(mySparseComm_type,     TXMO_glob,  TXMO)
      CALL scatterIJ(mySparseComm_type,     TYMO_glob,  TYMO)
      CALL scatterIJ(mySparseComm_type,     TZMO_glob,  TZMO)
#endif
      if (icase.lt.2) return

c**** icase=2: still serialized non-i/o parts of ocn dynamics
               ! for straits: mo,G0M,...,S0M,...,TRMO,...,opress
      CALL scatterIJ(mySparseComm_type,   OPRESS_glob,OPRESS)

      RETURN
      end subroutine scatter_straits_from_global

      subroutine alloc_odiff(grid)
      use DOMAIN_DECOMP_1D, only: dist_grid, get
      type (dist_grid) :: grid

      integer :: J_0H, J_1H

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      allocate( BYDXYV(grid%j_strt_halo:grid%j_stop_halo) )
      allocate( KHP   (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( KHV   (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( TANP  (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( TANV  (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDXV (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDXP (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDYV (grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BYDYP (grid%j_strt_halo:grid%j_stop_halo) )
      
      allocate( UXA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UXB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UXC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UYA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UYB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( UYC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VXA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VXB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VXC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VYA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VYB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
      allocate( VYC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

      end subroutine alloc_odiff

      END Module OCEAN

      Module OCEAN_DYN
!@sum  OCEAN_DYN contains variables used in ocean dynamics
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      Use OCEAN, Only : im,jm,lmo
!@var DH height of each ocean layer
!@var VBAR mean specific volume of each layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DH,VBAR !  (IM,JM,LMO)
     &     ,dZGdP

!@var GUP,GDN specific pot enthropy upper,lower part of layer (J/kg)
!@var SUP,SDN salinity at           upper,lower part of layer (1)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GUP,GDN,SUP,SDN

C**** momentum and mass fluxes
!@var MMI initial mass field (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: MMI  !  (IM,JM,LMO)
!@var SMU,SMV,SMW integrated mass fluxes
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SMU,SMV,SMW ! (IM,JM,LMO)
!@var CONV mass flux convergence
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CONV      !   (IM,JM,LMO)
!@var MU,MV,MW instantaneous mass fluxes
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: MU,MV,MW  !   (IM,JM,LMO)
C****
      END Module OCEAN_DYN

      Module SW2OCEAN
!@sum  SW2OCEAN variables for putting solar radiation into ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      Use OCEAN, Only : ze,lsrpd
      Implicit None
      REAL*8, DIMENSION(LSRPD) :: FSR,FSRZ,dFSRdZ,dFSRdZB
      REAL*8, PARAMETER :: RFRAC=.62d0, ZETA1=1.5d0, ZETA2=2d1

      CONTAINS

      SUBROUTINE init_solar
!@sum  init_solar calculates penetration of solar radiation
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      REAL*8 EF,EFZ,Z
      INTEGER L
      EF(Z) = RFRAC*EXP(-Z/ZETA1) + (1d0-RFRAC)*EXP(-Z/ZETA2)
      EFZ(Z)=ZETA1*RFRAC*EXP(-Z/ZETA1)+ZETA2*(1d0-RFRAC)*EXP(-Z/ZETA2)
C****
C**** Calculate the fraction of solar energy absorbed in each layer
C****
      do l=1,LSRPD
         FSR(l) = EF(ZE(l-1))
         FSRZ(l) = -3d0*(EF(ZE(l-1))+EF(ZE(l))) +
     *        6d0*(EFZ(ZE(l-1))-EFZ(ZE(l)))/(ZE(l)-ZE(l-1))
         dFSRdZB(l) = FSR(l)/(ZE(l)-ZE(l-1))
      end do
      do l=1,LSRPD-1
         dFSRdZ(l) = (FSR(l)-FSR(l+1))/(ZE(l)-ZE(l-1))
      end do
      dFSRdZ (LSRPD) = 0.
C****
      END SUBROUTINE init_solar
C****
      END Module SW2OCEAN

      SUBROUTINE alloc_ocean
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel
!@ver  1.0
      USE DOMAIN_DECOMP_ATM, only : dist_grid,agrid=>grid
      USE DOMAIN_DECOMP_1D, only : get, am_i_root
      USE OCEANR_DIM, only : ogrid,J_0H,J_1H,init_oceanr_grid  

      USE OCEANRES, only : IM=>IMO, JM=>JMO, LMO 

      USE OCEAN, only : MO,UO,VO,G0M,GXMO,GYMO,GZMO, OGEOZ,OGEOZ_SV
      USE OCEAN, only : UOD,VOD
      USE OCEAN, only :          S0M,SXMO,SYMO,SZMO, OPRESS,OPBOT
      USE OCEAN, only :
     *     MO_glob,UO_glob,VO_glob,UOD_glob,VOD_glob,
     *     G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob,
     *     S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob
      USE OCEAN, only : OXYP,OLAT2D_DG,OJ_BUDG,OWTBUDG,FOCEAN_loc
#ifdef TRACERS_OCEAN
      USE OCEAN, only : TRMO,TXMO,TYMO,TZMO
     *       ,TRMO_glob,TXMO_glob,TYMO_glob,TZMO_glob
      USE OCN_TRACER_COM, only : ntm
#endif
      USE OCEAN, only : nbyzmax,
     &     nbyzm,nbyzu,nbyzv,nbyzc,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc

      USE OCEAN, only: alloc_odiff
#ifdef TRACERS_OceanBiology
      USE obio_forc, only: alloc_obio_forc
      USE obio_com,  only: alloc_obio_com
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE TRACER_GASEXCH_COM, only: alloc_gasexch_com
#endif

      USE OCEAN_DYN, only : DH,VBAR, dZGdP, GUP,GDN, SUP,SDN
      USE OCEAN_DYN, only : MMI,SMU,SMV,SMW,CONV,MU,MV,MW

      IMPLICIT NONE

      INTEGER :: IER
      integer :: img, jmg, lmg

C*
C**** Define the ocean (Russell) grid 
C*
      call init_oceanr_grid  
C****

      CALL GET(ogrid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
 
      ALLOCATE(   MO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   UO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   VO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   UOD(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   VOD(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( G0M (IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( S0M (IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SZMO(IM,J_0H:J_1H,LMO), STAT = IER)

      if (am_i_root()) then
        img = im
        jmg = jm
        lmg = lmo
      else
        img = 1
        jmg = 1
        lmg = 1
      end if

      ALLOCATE(   MO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   UO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   VO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   UOD_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   VOD_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   G0M_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   GXMO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   GYMO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   GZMO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   S0M_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   SXMO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   SYMO_glob(IMG,JMG,LMG), STAT = IER)
      ALLOCATE(   SZMO_glob(IMG,JMG,LMG), STAT = IER)
      
      ALLOCATE( FOCEAN_loc(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OPRESS(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OPBOT (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OGEOZ (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OGEOZ_SV (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OXYP  (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OLAT2D_DG  (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OJ_BUDG  (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OWTBUDG  (IM,J_0H:J_1H), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE( TRMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TXMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TYMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TZMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)

      ALLOCATE( TRMO_glob(IMG,JMG,LMG,NTM), STAT = IER)
      ALLOCATE( TXMO_glob(IMG,JMG,LMG,NTM), STAT = IER)
      ALLOCATE( TYMO_glob(IMG,JMG,LMG,NTM), STAT = IER)
      ALLOCATE( TZMO_glob(IMG,JMG,LMG,NTM), STAT = IER)
#endif
!!!   ALLOCATE(   PO(IM,J_0H:J_1H,LMO), STAT = IER)
!!!   ALLOCATE(  PHI(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   DH(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( VBAR(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(dZGdP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  GUP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  GDN(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SUP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SDN(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  MMI(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SMU(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SMV(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(  SMW(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( CONV(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   MU(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   MV(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   MW(IM,J_0H:J_1H,LMO), STAT = IER)

C**** Necessary initiallisation?
      MU=0. ; MV=0. ; MW=0. ; CONV=0. ; MMI=0.
      UO=0. ; VO=0.
      SMU=0.; SMV=0.; SMW=0.

      ALLOCATE(NBYZM(J_0H:J_1H,LMO))
      ALLOCATE(NBYZU(J_0H:J_1H,LMO))
      ALLOCATE(NBYZV(J_0H:J_1H,LMO))
      ALLOCATE(NBYZC(J_0H:J_1H,LMO))
      ALLOCATE(I1YZM(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZM(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I1YZU(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZU(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I1YZV(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZV(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I1YZC(NBYZMAX,J_0H:J_1H,LMO))
      ALLOCATE(I2YZC(NBYZMAX,J_0H:J_1H,LMO))

c??   call ALLOC_GM_COM(agrid)
      call ALLOC_KPP_COM(ogrid)
      call alloc_odiag(ogrid)
      call alloc_afluxes(agrid)
      call ALLOC_OFLUXES(ogrid)

#ifdef TRACERS_OceanBiology
      call alloc_obio_forc
      call alloc_obio_com
#endif
#ifdef TRACERS_GASEXCH_ocean
      call alloc_gasexch_com
#endif
      call alloc_odiff(ogrid)

      return
      end subroutine alloc_ocean
