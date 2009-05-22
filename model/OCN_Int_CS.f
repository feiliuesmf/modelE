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
      REAL*8, INTENT(INOUT)  :: 
     &     aA(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,aN),
     &     aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: oA(oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)
      REAL*8, ALLOCATABLE :: oA_glob(:,:),oAtmp(:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing

      missing=-1.e30

c*** make sure that oN <= aN
      if (oN .gt. aN) 
     &     call stop_model("out of bounds",255)     
      
      allocate(oA_glob(oIM,oJM),oArea(oIM,oJM), oAtmp(oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

C***  regridding from atmospheric grid to ocean grid 
      call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(:,:,oN),
     &        oA_glob,oArea)   

      CALL OCN_UNPACK(ogrid, oA_glob, oAtmp)

      aA(:,:,oN)=oAtmp(:,:)

      DEALLOCATE(oA_glob,oArea,oAtmp)
      
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
      REAL*8, INTENT(INOUT) :: oA(oN,oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: oA_glob(:,:),oAtmp(:,:)
      REAL*8, ALLOCATABLE :: oArea(:,:)
      real*8 missing

      missing=-1.e30

c*** make sure that oNin <= aN
      if (oNin .gt. aN .OR. oNin .gt. oN) 
     &     call stop_model("out of bounds",255)
    
      allocate(oA_glob(oIM,oJM),oArea(oIM,oJM),oAtmp(oIM,
     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

      call repr_regrid_wt(xA2O,aWEIGHT,missing,aA(oNin,:,:),
     &        oA_glob(:,:),oArea)

C***  Scatter global array oA_glob to the ocean grid

      CALL OCN_UNPACK_COL (ogrid, oA_glob, oAtmp)

      oA(oNin,:,:)=oAtmp(:,:)

      DEALLOCATE(oA_glob,oArea,oAtmp)
      
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


      SUBROUTINE INT_AG2OG_Vector1(aU,aV,oU,oV,aWEIGHT,aFOCEAN,aN,oN) 
!@sum INT_AG2OG_Vector1 is for conversion vector from atm. CS grid to ocean A grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA,AM_I_ROOT
      Use OCEANR_DIM,       only : oGRID
      USE OCEAN,            only : oSINI=>SINIC, oCOSI=>COSIC   
      use regrid_com, only : xA2O

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: aN, oN
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
      REAL*8, intent(in) :: 
     &     aFOCEAN(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
      REAL*8, DIMENSION(oIM,oJM) :: oFtemp 
      REAL*8  oUsp, oVsp, oUnp, oVnp
      REAL*8, ALLOCATABLE :: oUtmp(:,:),oVtmp(:,:)
      REAL*8, ALLOCATABLE :: oUt(:,:),oVt(:,:)
      REAL*8, allocatable :: oArea(:,:),aftemp(:,:)
      integer :: i,j
      real*8 :: missing
      

      missing=-1.e30

      ALLOCATE(
     &     aftemp(aGRID%I_STRT:aGRID%I_STOP,
     &     aGRID%J_STRT:aGRID%J_STOP),
     &     oArea(oIM,oJM),
     &     oUt(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     oVt(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     oUtmp(oIM,oJM),oVtmp(oIM,oJM) )

      do j=aGRID%J_STRT,aGRID%J_STOP
         do i=aGRID%I_STRT,aGRID%I_STOP
            if (aFOCEAN(i,j).gt.0.) then
               aftemp(i,j) = aU(i,j,oN) !  kg/m s
            endif
        enddo
      enddo

      call repr_regrid_wt(xA2O,aWEIGHT,missing,aftemp,
     &     oUtmp,oArea)  
      
      oUsp =  SUM(oUtmp(:,  1)*oCOSI(:))*2/oIM
      oVsp = -SUM(oUtmp(:,  1)*oSINI(:))*2/oIM
      oUtmp(1,  1) = oUsp
     
      oUnp =  SUM(oUtmp(:,oJM)*oCOSI(:))*2/oIM
      oVnp =  SUM(oUtmp(:,oJM)*oSINI(:))*2/oIM
      oUtmp(1,oJM) = oUnp
     
      
      do j=aGRID%J_STRT,aGRID%J_STOP
         do i=aGRID%I_STRT,aGRID%I_STOP
            if (aFOCEAN(i,j).gt.0.) then
               aftemp(i,j) = aV(i,j,oN) !  kg/m s
            endif
         enddo
      enddo
      
      call repr_regrid_wt(xA2O,aWEIGHT,missing,aftemp,
     &     oVtmp,oArea)  

      oVtmp(1,  1) = oVsp
      oVtmp(1,oJM) = oVnp
     
C***  Scatter global array oA_glob to the ocean grid
      CALL OCN_UNPACK (ogrid, oUtmp, oUt)
      CALL OCN_UNPACK (ogrid, oVtmp, oVt)

      oU(:,:,oN)=oUt
      oV(:,:,oN)=oVt
      
      DEALLOCATE(aftemp,oArea,oUt,oVt,oUtmp,oVtmp)
      
      END SUBROUTINE INT_AG2OG_Vector1
c*


      SUBROUTINE INT_AG2OG_Vector2(aU,aV, oU,oV, aWEIGHT, IVSPO,IVNPO)
!@sum INT_AG2OG_Vector2 is for conversion vector from atm. (U,V) grid to ocean (U,V) grid 
!@auth Larissa Nazarenko, Denis Gueyffier
      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA,get
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

      oUnp = 0.0   
      oU_glob(  oIM,oJM) = oUnp
      oU_glob(IVNPO,  1) = oUnp
      
      oUsp = 0.0  
      oU_glob(IVSPO,  1) = oUsp

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


      MODULE INT_OG2AG_MOD

!@sum INT_OG2AG_MOD contains subroutines for conversion 2D, 3D, etc. 
!!    arrays from ocean to the atm. grid 
!@auth Larissa Nazarenko
!@ver  1.0
      PRIVATE
      PUBLIC INT_OG2AG

      Interface INT_OG2AG
      Module Procedure INT_OG2AG_test
      Module Procedure INT_OG2AG_2Da
      Module Procedure INT_OG2AG_3Da
      Module Procedure INT_OG2AG_3Db
      Module Procedure INT_OG2AG_4Da
      Module Procedure INT_OG2AG_Vector1
      End Interface

      contains

      SUBROUTINE INT_OG2AG_2Da(oA,aA,oWEIGHT,CopyPole)

!@sum regridding 2D arrays from ocean to the CS atm. grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE OCEAN,      only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,ATM_UNPACK=>UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : get
      Use OCEANR_DIM,       only : oGRID
      use regrid_com, only : xO2A

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: CopyPole
      REAL*8 :: oWEIGHT(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8 :: aA(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: aA_glob(:,:,:),aArea(:,:,:),oFtemp(:,:)
      real*8 :: missing
      logical :: HAVE_NORTH_POLE
      missing=-1.e30

      ALLOCATE(aA_glob(aIM,aJM,6),aArea(aIM,aJM,6),
     &     oFtemp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

C***  Interpolate aA_glob from ocean grid to atmospheric grid 
      call get(ogrid, HAVE_NORTH_POLE=HAVE_NORTH_POLE) 

      if (CopyPole .and. HAVE_NORTH_POLE) 
     &     oWEIGHT(2:oIM,oJM) = oWEIGHT(1,oJM)

      oFtemp(:,:) = oA(:,:)
      if (HAVE_NORTH_POLE) oFtemp(2:oIM,oJM) = oFtemp(1,oJM)

      call repr_regrid_wt(xO2A,oWEIGHT,missing,oFtemp,aA_glob,aArea)

C***  Scatter global array aA_glob to the atmospheric grid
      CALL ATM_UNPACK (agrid, aA_glob, aA)

      DEALLOCATE(aA_glob,aArea,oFtemp)
      
      END SUBROUTINE INT_OG2AG_2Da
c*


      SUBROUTINE INT_OG2AG_3Da(oA,aA,oWEIGHT,NT)

!@sum regridding 3D arrays from ocean to CS atm. grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,ATM_PACK=>PACK_DATA,
     &     ATM_UNPACK=>UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : get
      Use OCEANR_DIM,       only : oGRID
      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN
      use regrid_com, only : xO2A

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: 
     &        oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      integer, intent(in) :: NT 
      REAL*8  :: 
     &  aA(NT,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: 
     *  oA(NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: aA_glob(:,:,:,:), aFOCEAN(:,:,:),
     &     aArea(:,:,:),aFtemp(:,:,:),oFtemp(:,:)
      integer :: i,j,k,N
      logical :: HAVE_NORTH_POLE
      real*8 :: missing

      missing=-1.e30

      call get(ogrid, HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      ALLOCATE(aFOCEAN(aIM,aJM,6),aA_glob(NT,aIM,aJM,6),
     &     aArea(aIM,aJM,6),aFtemp(aIM,aJM,6),
     &     oFtemp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

      call ATM_PACK (aGRID,aFOCEAN_loc,aFOCEAN)

      do N=1,NT
          oFtemp(:,:) = oA(N,:,:)
          if (HAVE_NORTH_POLE) 
     &         oFtemp(2:oIM,oJM) = oFtemp(1,oJM)

         call repr_regrid_wt(xO2A,oWEIGHT,missing,oFtemp,
     &        aFtemp,aArea)
         do k=1,6
            do J=1,aJM
               do I=1,aIM
                  IF (aFOCEAN(I,J,K).gt.0.) THEN
                     aA_glob(N,I,J,K) = aFtemp(I,J,K)  
                  END IF
               enddo
            enddo
         enddo
      enddo

C***  Scatter global array oA_glob to the ocean grid

      CALL ATM_UNPACK(agrid, aA_glob, aA)

      DEALLOCATE(aFocean, aA_glob, aArea, aFtemp, oFtemp)
      
      END SUBROUTINE INT_OG2AG_3Da
c*
      subroutine int_og2ag_test(aA,aN)
      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,UNPACK_DATA,get
      IMPLICIT NONE
      integer, intent(in) :: aN
      REAL*8, intent(inout):: 
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO,aN)
      REAL*8, ALLOCATABLE :: aA_glob(:,:,:), aAtmp(:,:)
      integer :: N

      ALLOCATE(
     &     aAtmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &           aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     aA_glob(aIM,aJM,6) )

      do N=1,aN
         aA_glob=1.
         call UNPACK_DATA(agrid, aA_glob, aAtmp)
         aA(:,:,N)=aAtmp
      enddo

      DEALLOCATE(aA_glob,aAtmp)

      end subroutine int_og2ag_test

      SUBROUTINE INT_OG2AG_3Db(oA,aA,oWEIGHT,oN,aN,CopyPole)

!@sum regridding 3D arrays from ocean to CS atm. grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : get
      Use OCEANR_DIM,       only : ogrid
      use regrid_com, only : xO2A
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: CopyPole
      REAL*8 :: 
     &     oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      integer, intent(in) :: aN,oN
      REAL*8, intent(inout):: 
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO,aN) 
      REAL*8, intent(in) :: 
     *     oA(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,oN)
      REAL*8, ALLOCATABLE :: aA_glob(:,:,:),aArea(:,:,:),oFtemp(:,:),
     &     aAtmp(:,:)
      integer :: N
      real*8 :: missing
      logical :: HAVE_NORTH_POLE

      missing=-1.e30

      ALLOCATE(aArea(aIM,aJM,6),
     &     aAtmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &           aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     oFtemp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     aA_glob(aIM,aJM,6)
     &        )

      call get(ogrid, HAVE_NORTH_POLE=HAVE_NORTH_POLE) 

      if (CopyPole .and. HAVE_NORTH_POLE) 
     &     oWEIGHT(2:oIM,oJM) = oWEIGHT(1,oJM)
      
      do N=1,aN
         oFtemp(:,:) = oA(:,:,N)
         if (HAVE_NORTH_POLE) oFtemp(2:oIM,oJM) = oFtemp(1,oJM)
         aA_glob=1.
         call UNPACK_DATA(agrid, aA_glob, aAtmp)
         aA(:,:,N)=aAtmp

      enddo

      DEALLOCATE(aArea,oFtemp,aA_glob,aAtmp)
      
      END SUBROUTINE INT_OG2AG_3Db
c*



      SUBROUTINE INT_OG2AG_4Da(oA,aA,oWEIGHT,NT,NTM)

!@sum regridding 4D arrays from ocean to the CS atm. grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN, only : oIM=>im,oJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,ATM_PACK=>PACK_DATA,
     &     ATM_UNPACK=>UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : OCN_PACK_COL=>PACK_COLUMN,get
      Use OCEANR_DIM,       only : oGRID
      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN
      use regrid_com, only : xO2A

      IMPLICIT NONE
      INTEGER :: IER, NTM,NT,N1,N2, I,J,K,N
      REAL*8, INTENT(IN)  :: 
     &     oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) 
      REAL*8  :: 
     &  aA(NTM,NT,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) 
      REAL*8 :: 
     &  oA(NTM,NT,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8, ALLOCATABLE :: aA_glob(:,:,:,:,:), aArea(:,:,:),
     &                       aFOCEAN(:,:,:), aFtemp(:,:,:),
     &                       oFtemp(:,:)
      real*8 :: missing

      logical :: HAVE_NORTH_POLE

      call get(ogrid, HAVE_NORTH_POLE=HAVE_NORTH_POLE) 

      missing=-1.e30

      ALLOCATE(
     &     aFOCEAN(aIM,aJM,6),
     &     aA_glob(NTM,NT,aIM,aJM,6),
     &     aFtemp(aIM,aJM,6),
     &     aArea(aIM,aJM,6),
     &     oFtemp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

C***  Gather 3D array on atmospheric grid into the global array

      CALL ATM_PACK(agrid, aFOCEAN_loc, aFOCEAN)
      
      DO N1=1,NTM
         DO N2=1,NT
            oFtemp(:,:) = oA(N1,N2,:,:)
            oFtemp(2:oIM,oJM) = oFtemp(1,oJM)
            
            call repr_regrid_wt(xO2A,oWEIGHT,missing,oFtemp,
     &           aFtemp,aArea)

            do K=1,6
               DO J=1,aJM
                  DO I=1,aIM
                     IF (aFOCEAN(I,J,K).gt.0.) THEN
                        aA_glob(N1,N2,I,J,K) = aFtemp(I,J,K)  
                     END IF
                  enddo
               enddo
            enddo
         enddo
      enddo

C***  Scatter global array oA_glob to the ocean grid

      call ATM_UNPACK (agrid, aA_glob, aA)

      deallocate(aA_glob,aArea,aFOCEAN,aFtemp,oFtemp)

      END SUBROUTINE INT_OG2AG_4Da


      SUBROUTINE INT_OG2AG_Vector1(oUO1,oVO1,aUO1,aVO1,oWEIGHT,
     &     IVSPO,IVNPO) 
!@sum INT_AG2OG_Vector1 is for conversion vector from ocean grid to atm. CS grid 
!@auth Larissa Nazarenko, Denis Gueyffier

      USE OCEAN, only : oIM=>im,oJM=>jm
      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,get,
     &     ATM_UNPACK=>UNPACK_DATA
      Use OCEANR_DIM,       only : oGRID
      USE OCEAN,            only : oSINU=>SINU, oCOSU=>COSU   
      use regrid_com, only : xO2A

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IVSPO,IVNPO
      REAL*8 ::  
     &     aUO1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO), 
     &     aVO1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO),
     &     oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8 ::
     &     oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      REAL*8  oUsp, oVsp, oUnp, oVnp
      REAL*8, ALLOCATABLE :: aUO1_glob(:,:,:),aVO1_glob(:,:,:)
      REAL*8, allocatable :: aArea(:,:,:)
      logical :: HAVE_NORTH_POLE,HAVE_SOUTH_POLE
      real*8 :: missing

      missing=-1.e30
      
      ALLOCATE(
     &     aUO1_glob(aIM,aJM,6), 
     &     aVO1_glob(aIM,aJM,6),
     &     aArea(oIM,oJM,6)
     &     )
      
      call get(ogrid, HAVE_NORTH_POLE=HAVE_NORTH_POLE,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE) 

!!!   U velocity for the 1st ocean layer.
    
      if (HAVE_SOUTH_POLE) then
         oVsp = oUO1(IVSPO,  1)
         oUO1(:,  1) = oUO1(oIM,  1)*oCOSU(:) - oVsp*oSINU(:)
      endif

      if (HAVE_NORTH_POLE) then
         oVnp = oUO1(IVNPO,oJM)
         oUO1(:,oJM) = oUO1(oIM,oJM)*oCOSU(:) + oVnp*oSINU(:)
      endif

      call repr_regrid_wt(xO2A,oWEIGHT,missing,oUO1,
     &     aUO1_glob,aArea)  
         
      call repr_regrid_wt(xO2A,oWEIGHT,missing,oVO1,
     &     aVO1_glob,aArea)  
      
      
C***  Scatter global array oA_glob to the ocean grid
      CALL ATM_UNPACK (agrid, aUO1_glob, aUO1)
      CALL ATM_UNPACK (agrid, aVO1_glob, aVO1)
      
      DEALLOCATE(aUO1_glob,aVO1_glob,aArea)
      
      END SUBROUTINE INT_OG2AG_Vector1
c*

      END MODULE INT_OG2AG_MOD
