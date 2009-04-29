#include "rundeck_opts.h"

      MODULE INT_AG2OG_MOD
!@sum INT_AG2OG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from atm. to the ocean grid 
!@auth Larissa Nazarenko

      PRIVATE
      PUBLIC INT_AG2OG

      Interface INT_AG2OG
      Module Procedure INT_AG2OG_2Da
      Module Procedure INT_AG2OG_2Db
      Module Procedure INT_AG2OG_3Da
      Module Procedure INT_AG2OG_3Db
      Module Procedure INT_AG2OG_3Dc
      Module Procedure INT_AG2OG_4Da
      Module Procedure INT_AG2OG_4Db
      Module Procedure INT_AG2OG_Vector1
      Module Procedure INT_AG2OG_Vector2
      End Interface

      contains

      SUBROUTINE INT_AG2OG_2Da(aA,oA,aWEIGHT)
!@sum INT_AG2OG_2D regridding 2D arrays from CS atm. to ocean grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: aA(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: oA_glob(:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing

      missing=-1.e30
    
      allocate(oA_glob(oIM,oJM),oArea(oIM,oJM))

C***  regridding from atmospheric grid to ocean grid 
      call repr_regrid_wt(xA2O,aWEIGHT,missing,aA,oA_glob,oArea)   

C***  Scatter global array oA_glob to the ocean grid

      call OCN_UNPACK (ogrid, oA_glob, oA)

      deallocate(oA_glob,oArea)
      
      END SUBROUTINE INT_AG2OG_2Da
c*


      SUBROUTINE INT_AG2OG_2Db(aA1,aA2,oA,aWEIGHT) 

!@sum INT_AG2OG_2D regridding 2D arrays from CS atm. to ocean grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: 
     &     aA1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aA2(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: oA_glob(:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:),atmp(:,:)
      real*8 missing

      missing=-1.e30
    
      allocate(oA_glob(oIM,oJM),oArea(oIM,oJM),
     &     atmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

C***  regridding from atmospheric grid to ocean grid 
      atmp(:,:) = aA1(:,:)*aA2(:,:)
      call repr_regrid_wt(xA2O,aWEIGHT,missing,atmp,oA_glob,oArea)   

C***  Scatter global array oA_glob to the ocean grid

      CALL OCN_UNPACK (ogrid, oA_glob, oA)

      deallocate(oA_glob, oArea)
      
      END SUBROUTINE INT_AG2OG_2Db


      SUBROUTINE INT_AG2OG_3Da(aA,oA,aWEIGHT,NT)
!@sum INT_AG2OG_3D regridding 3D arrays from CS atm. to ocean grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK_COL=>UNPACK_COLUMN
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: 
     &     aA(NT,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(NT,oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      integer, intent(in) :: NT 
      REAL*8, ALLOCATABLE :: oA_glob(:,:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing
      integer :: N

      missing=-1.e30
    
      allocate(oA_glob(NT,oIM,oJM),oArea(oIM,oJM))

C***  regridding from atmospheric grid to ocean grid 
      do N=1,NT
         call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(N,:,:),
     &        oA_glob(N,:,:),oArea)   
      enddo

C***  Scatter global array oA_glob to the ocean grid
      call OCN_UNPACK_COL (ogrid, oA_glob, oA)

      deallocate(oA_glob,oArea)

      END SUBROUTINE INT_AG2OG_3Da
c*


      SUBROUTINE INT_AG2OG_3Db(aA,oA,aWEIGHT, aN, oN)
!@sum INT_AG2OG_3D regridding 3D arrays from CS atm. to ocean grid 
!@+   extract oN-th component of original aA array
!@+   this subroutine and its latlon couterpart may experience issues if oN > 1. 
!@    the routine should probably return only the oN-th component of oA (i.e. a 2D array)
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_COLUMN
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: aN, oN
      REAL*8, INTENT(IN)  :: 
     &     aA(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)
      REAL*8, ALLOCATABLE :: oA_glob(:,:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing

      missing=-1.e30

c*** make sure that oN <= aN
      if (oN .gt. aN) 
     &     call stop_model("out of bounds",255)     
      
      allocate(oA_glob(oIM,oJM,oN),oArea(oIM,oJM))

C***  regridding from atmospheric grid to ocean grid 
      call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(:,:,oN),
     &        oA_glob(:,:,oN),oArea)   

      CALL OCN_UNPACK(ogrid, oA_glob, oA)

      DEALLOCATE(oA_glob,oArea)
      
      END SUBROUTINE INT_AG2OG_3Db
c*


      SUBROUTINE INT_AG2OG_3Dc(aA,oA,aWEIGHT, aN,oN,oNin)
!@sum INT_AG2OG_3D regridding 3D arrays from CS atm. to ocean grid 
!@+   extract oNin-th component of original aA array
!@+   this subroutine and its latlon couterpart should be modified 
!@+   so that it takes as argument only the oNin-th 
!@+   component of aA and returns the oNin-th component of oA, 
!@+   i.e. we actually pass 2D arrays instead of 3D 
!@auth Larissa Nazarenko, Denis Gueyffier
      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only :  OCN_UNPACK_COL=>UNPACK_COLUMN
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: aN, oN, oNin
      REAL*8, INTENT(IN)  :: 
     &     aA(aN,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(oN,oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: oA_glob(:,:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing

      missing=-1.e30

c*** make sure that oNin <= aN
      if (oNin .gt. aN .OR. oNin .gt. oN) 
     &     call stop_model("out of bounds",255)
    
      allocate(oA_glob(oN,oIM,oJM),oArea(oIM,oJM))

      call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(oNin,:,:),
     &        oA_glob(oNin,:,:),oArea)

C***  Scatter global array oA_glob to the ocean grid

      CALL OCN_UNPACK_COL (ogrid, oA_glob, oA)

      DEALLOCATE(oA_glob,oArea)
      
      END SUBROUTINE INT_AG2OG_3Dc
c*


      SUBROUTINE INT_AG2OG_4Da(aA,oA,aWEIGHT, NT, aN,oN)
!@+   extract oN-th component of original aA array
!@auth Larissa Nazarenko, Denis Gueyffier
      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK_BLK=>UNPACK_BLOCK
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: aN, oN, NT
      REAL*8, INTENT(IN)  :: 
     &     aA(NT,aN,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(NT,aN,oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: oA_glob(:,:,:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing
      integer :: N

      missing=-1.e30

c*** make sure that oN <= aN
      if (oN .gt. aN) 
     &     call stop_model("out of bounds",255)     
      
      allocate(oA_glob(NT,oN,oIM,oJM),oArea(oIM,oJM))

C***  regridding from atmospheric grid to ocean grid 
      do N=1,NT
         call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(N,oN,:,:),
     &        oA_glob(N,oN,:,:),oArea)   
      enddo

C***  Scatter global array oA_glob to the ocean grid
      CALL OCN_UNPACK_BLK (ogrid, oA_glob, oA)

      DEALLOCATE(oA_glob,oArea)

      END SUBROUTINE INT_AG2OG_4Da
c*

      SUBROUTINE INT_AG2OG_4Db(aA,oA,oB,aWEIGHT, NT, aN,oN)
!@auth Larissa Nazarenko, Denis Gueyffier
      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK_BLK=>UNPACK_BLOCK
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: aN, oN, NT
      REAL*8, INTENT(IN)  :: 
     &     aA(NT,aN,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(NT,aN,oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, DIMENSION(NT,oIM,oJM), INTENT(OUT) :: oB
      REAL*8, ALLOCATABLE :: oA_glob(:,:,:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing
      integer :: N

      missing=-1.e30

c*** make sure that oN <= aN
      if (oN .gt. aN) 
     &     call stop_model("out of bounds",255)     
      
      allocate(oA_glob(NT,oN,oIM,oJM),oArea(oIM,oJM))

C***  regridding from atmospheric grid to ocean grid 
      do N=1,NT
         call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(N,oN,:,:),
     &        oA_glob(N,oN,:,:),oArea)   
         oB(:,:,N)=oA_glob(N,oN,:,:)
      enddo

C***  Scatter global array oA_glob to the ocean grid

      CALL OCN_UNPACK_BLK (ogrid, oA_glob, oA)

      DEALLOCATE(oA_glob,oArea)

      END SUBROUTINE INT_AG2OG_4Db
c*


      SUBROUTINE INT_AG2OG_Vector1(aU,aV, oU,oV, aWEIGHT,aFOCEAN,aIMAXJ, 
     *           aSINI,aCOSI, oSINI,oCOSI, aN,oN) 
!@sum INT_AG2OG_Vector1 is for conversion vector from atm. CS grid to ocean A grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: aN, oN
      INTEGER, DIMENSION(aJM), intent(in) :: aIMAXJ  !not used. Kept only to satisfy interface contract
      REAL*8, INTENT(IN)  ::  
     &     aU(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN), 
     &     aV(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) ::
     &     oU(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)
     &     ,oV(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)
      REAL*8, DIMENSION(aIM,aJM), intent(in) :: aFOCEAN
      REAL*8, DIMENSION(aIM), intent(in) :: aSINI,aCOSI
      REAL*8, DIMENSION(oIM), intent(in) :: oSINI,oCOSI
      REAL*8, DIMENSION(oIM,oJM) :: oFtemp 
      REAL*8  oUsp, oVsp, oUnp, oVnp
      REAL*8, ALLOCATABLE :: oU_glob(:,:,:),oV_glob(:,:,:)
      REAL*8, allocatable :: oArea(:,:),aftemp(:,:)
      integer :: i,j
      real*8 :: missing

      missing=-1.e30

      ALLOCATE(
     &     oU_glob(oIM,oJM,oN), 
     &     oV_glob(oIM,oJM,oN),
     &     aftemp(aGRID%I_STRT:aGRID%I_STOP,
     &     aGRID%J_STRT:aGRID%J_STOP),
     &     oArea(oIM,oJM)
     &     )
      
      do j=oGRID%J_STRT,oGRID%J_STOP
         do i=oGRID%I_STRT,oGRID%I_STOP
            if (aFOCEAN(i,j).gt.0.) then
               aftemp(i,j) = aU(i,j,oN) !  kg/m s
            endif
        enddo
      enddo

      call repr_regrid_wt(xA2O,aWEIGHT,missing,aftemp,
     &     oU_glob(:,:,oN),oArea)  
      
      oUsp =  SUM(oU_glob(:,  1,oN)*oCOSI(:))*2/oIM
      oVsp = -SUM(oU_glob(:,  1,oN)*oSINI(:))*2/oIM
      oUnp =  SUM(oU_glob(:,oJM,oN)*oCOSI(:))*2/oIM
      oVnp =  SUM(oU_glob(:,oJM,oN)*oSINI(:))*2/oIM
      oU_glob(1,  1,oN) = oUsp
      oU_glob(1,oJM,oN) = oUnp
      
      do j=oGRID%J_STRT,oGRID%J_STOP
         do i=oGRID%I_STRT,oGRID%I_STOP
            if (aFOCEAN(i,j).gt.0.) then
               aftemp(i,j) = aV(i,j,oN) !  kg/m s
            endif
         enddo
      enddo
      
      call repr_regrid_wt(xA2O,aWEIGHT,missing,aftemp,
     &     oV_glob(:,:,oN),oArea)  
      
      oV_glob(1,  1,oN) = oVsp
      oV_glob(1,oJM,oN) = oVnp
     
      
C***  Scatter global array oA_glob to the ocean grid
      CALL OCN_UNPACK (ogrid, oU_glob, oU)
      CALL OCN_UNPACK (ogrid, oV_glob, oV)
      
      DEALLOCATE(oU_glob,oV_glob,aftemp,oArea)
      
      END SUBROUTINE INT_AG2OG_Vector1
c*


      SUBROUTINE INT_AG2OG_Vector2(aU,aV, oU,oV, aWEIGHT, IVSPO,IVNPO)
!@sum INT_AG2OG_Vector2 is for conversion vector from atm. (U,V) grid to ocean (U,V) grid 
!@auth Larissa Nazarenko, Denis Gueyffier
      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xA2O

      IMPLICIT NONE
      REAL*8, INTENT(IN)  ::  
     &     aU(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO), 
     &     aV(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8, INTENT(OUT) ::
     &     oU(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     &     ,oV(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      integer, intent(in) :: IVSPO,IVNPO
      REAL*8, DIMENSION(oIM,oJM) :: oFtemp 
      REAL*8  oUsp, oUnp
      REAL*8, ALLOCATABLE :: oU_glob(:,:),oV_glob(:,:),oArea(:,:)
      real*8 :: missing

      missing=-1.e30

      allocate(oU_glob(oIM,oJM),oV_glob(oIM,oJM),oArea(oIM,oJM))

C***  Interpolate aU from atmospheric grid to ocean grid 

      call repr_regrid_wt(xA2O,aWEIGHT,missing,aU,
     &     oU_glob,oArea)  

        oUsp = 0.0  
        oUnp = 0.0   

        oU_glob(IVSPO,  1) = oUsp
        oU_glob(  oIM,oJM) = oUnp
        oU_glob(IVNPO,  1) = oUnp

      call repr_regrid_wt(xA2O,aWEIGHT,missing,aV,
     &     oV_glob,oArea)  

        oV_glob(:,oJM) = 0.0                       !!  should not be used

C***  Scatter global array oA_glob to the ocean grid

      CALL OCN_UNPACK(ogrid, oU_glob, oU)
      CALL OCN_UNPACK(ogrid, oV_glob, oV)

      DEALLOCATE(oU_glob,oV_glob,oArea)
      
      END SUBROUTINE INT_AG2OG_Vector2
c*

      END MODULE INT_AG2OG_MOD
c*
