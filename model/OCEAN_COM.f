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
      Use GEOM,      Only: imaxj
      Use OCEANRES,  Only: IM=>IMO,JM=>JMO, LMO, LMO_MIN, LSRPD, dZO        
#ifdef TRACERS_OCEAN
      Use OCN_TRACER_COM, Only : ntm
#endif
      Use SparseCommunicator_mod
      Implicit None
      Integer*4,Parameter ::
     *  IVSP = 3*IM/4,      !  V at south pole is stored in U(IVSP,1)
     *  IVNP =   IM/4,      !  V at north pole is stored in U(IVNP,JM)
     *  JLATD = 180/(JM-1)  !  LATitudinal spacing in degrees
      Real*8,Parameter ::
     *  DLON  = TWOPI/IM,        !  LONgitudinal spacing in radians
     *  DLAT  = TWOPI*JLATD/360, !  LATitudinal spacing in radians
     *  FJEQ  = .5*(1+JM)        !  J coordinate of the EQuator
      REAL*8, PARAMETER :: FIM=IM, BYIM=1./FIM

      Integer*4 ::
     *  OBottom_drag = 1,  !  @dbparam use ocean bottom drag routine
     *  OCoastal_drag = 1  !  @dbparam use ocean coastal drag routine

      REAL*8 ::
     *     oc_salt_mean = -999. ! @dbparam mean salinity of ocean (if set)
#ifdef TRACERS_OCEAN
     *     , oc_tracer_mean(ntm) = -999. ! @dbparam mean tracer ratio of ocean
                                         ! Note: in permil for water isotopes
#endif 

!@var MO mass of ocean (kg/m^2)
!@var UO E-W velocity on C-grid (m/s)
!@var VO N-S velocity on C-grid (m/s)
!@var G0M,GXMO,GYMO,GZMO pot. enthalpy of ocean (+moments) (J)
!@var S0M,SXMO,SYMO,SZMO salinity of ocean (+moments) (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: MO,UO,VO,
     *     G0M,GXMO,GYMO,GZMO, S0M,SXMO,SYMO,SZMO
C**** Global arrays needed for i/o, GM,straits,odiff ?
      REAL*8, DIMENSION(IM,JM,LMO) :: MO_glob, UO_glob, VO_glob,
     *     G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob,
     *     S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob

      Real*8 UONP(LMO), !  U component at north pole, points down 90W
     *       VONP(LMO)  !  V component at north pole, points down 0 (GM)

C**** ocean geometry (should this be in a separate module?)
      REAL*8, DIMENSION(JM) :: DXYPO,DXPO,DYPO,DXVO,DYVO
     *     ,COSPO,SINPO,DXYVO,DXYSO,DXYNO,RAMVS,RAMVN,RLAT,BYDXYPO
      REAL*8, DIMENSION(0:JM) :: COSVO
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
      INTEGER, DIMENSION(IM,JM) :: LMM,LMU,LMV
!@var RATOC,ROCAT Ratio of areas for converting atm. fluxes to ocean
      REAL*8, DIMENSION(JM) :: RATOC,ROCAT
      Integer*4 J40S, !  maximum grid cell below 40S (used in OPFIL)
     *           J1O  !  most southern latitude (J) where ocean exists

      REAL*8, DIMENSION(IM,JM) :: HATMO,HOCEAN,FOCEAN
C**** ocean related parameters
      INTEGER NDYNO,MDYNO,MSGSO
!@dbparam DTO timestep for ocean dynamics (s)
      REAL*8 :: DTO=450.        ! default. setable parameter
      REAL*8 DTOFS,DTOLF,DTS,BYDTS

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
      REAL*8, DIMENSION(IM,JM,LMO,NTM) :: ! for i/o, straits, GM ?
     *       TRMO_glob,TXMO_glob,TYMO_glob,TZMO_glob
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRMO,TXMO,TYMO,TZMO
#endif
     
      type (SparseCommunicator_type), save :: mySparseComm_type

      contains

      subroutine gather_ocean (icase)

      use domain_decomp, only: pack_data
      use OCEANR_DIM, only : grid=>ogrid

      integer, intent(in) :: icase

c**** icase=0:full i/o, 1:ini_straits, 2:serialized ocn dynamics
      CALL PACK_DATA(grid,   MO   ,    MO_glob)
      if (icase.ne.1) then      ! needed for i/o and ODIFF only
        CALL PACK_DATA(grid,   UO   ,    UO_glob)
        CALL PACK_DATA(grid,   VO   ,    VO_glob)
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

      use domain_decomp, only: unpack_data
      use OCEANR_DIM, only : grid=>ogrid

      integer, intent(in) :: icase

c**** icase=-1: i/o no_trc 0:full i/o, 1:ini_straits, 2:serial ocn dyn
      CALL UNPACK_DATA(grid,       MO_glob,   MO )
      if (icase.lt.1) then            ! needed for i/o only
        CALL UNPACK_DATA(grid,       UO_glob,   UO )
        CALL UNPACK_DATA(grid,       VO_glob,   VO )
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

      subroutine gather_ocean_straits (icase)
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
      end subroutine gather_ocean_straits

      subroutine scatter_ocean_straits (icase)
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
      end subroutine scatter_ocean_straits

      END Module OCEAN

      Module OCEAN_DYN
!@sum  OCEAN_DYN contains variables used in ocean dynamics
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0
      Use OCEAN, Only : im,jm,lmo
!@var DH height of each ocean layer
!@var VBAR mean specific volume of each layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DH,VBAR !  (IM,JM,LMO)
      REAL*8, DIMENSION(IM,JM,LMO) :: DH_glob,VBAR_glob ! for serial ocnGM ???

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
      Use OCEAN, Only : ze
      Implicit None
      INTEGER, PARAMETER :: LSRPD = 3
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

      USE DOMAIN_DECOMP, only : dist_grid,get

      USE OCEANR_DIM

      USE OCEANRES, only : IM=>IMO, JM=>JMO, LMO 

      USE OCEAN, only : MO,UO,VO,G0M,GXMO,GYMO,GZMO, OGEOZ,OGEOZ_SV
      USE OCEAN, only :          S0M,SXMO,SYMO,SZMO, OPRESS,OPBOT
#ifdef TRACERS_OCEAN
      USE OCEAN, only : TRMO,TXMO,TYMO,TZMO
      USE OCN_TRACER_COM, only : ntm
#endif
      USE OCEAN_DYN, only : DH,VBAR, GUP,GDN, SUP,SDN
      USE OCEAN_DYN, only : MMI,SMU,SMV,SMW,CONV,MU,MV,MW

      IMPLICIT NONE

      INTEGER :: IER
C*
C**** Define the ocean (Russell) grid 
C*
      call init_oceanr_grid  
C****

      CALL GET(ogrid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
 
!      write (555,*) 
!      write (555,*) ' alloc_ocean: halos', J_0H,J_1H

      ALLOCATE(   MO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   UO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   VO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( G0M (IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( GZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( S0M (IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SXMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SYMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( SZMO(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( OPRESS(IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OPBOT (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OGEOZ (IM,J_0H:J_1H), STAT = IER)
      ALLOCATE( OGEOZ_SV (IM,J_0H:J_1H), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE( TRMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TXMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TYMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TZMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
#endif
!!!   ALLOCATE(   PO(IM,J_0H:J_1H,LMO), STAT = IER)
!!!   ALLOCATE(  PHI(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   DH(IM,J_0H:J_1H,LMO), STAT = IER)
!!!   ALLOCATE(DZGDP(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE( VBAR(IM,J_0H:J_1H,LMO), STAT = IER)
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
      MU=0. ; MV=0. ; MW=0. ; CONV=0.

c??   call ALLOC_GM_COM(grid)
      call ALLOC_KPP_COM(ogrid)
      call alloc_odiag(ogrid)

      return
      end subroutine alloc_ocean
