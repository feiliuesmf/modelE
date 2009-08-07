c**** 
C**** ICEDYN_DRV.f    Sea ICE DYNamics    2006/12/21
C****
#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
      MODULE ICEDYN_COM
!@sum  ICEDYN_COM holds global variables for dynamic sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP_ATM, only : DIST_GRID
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      USE DIAG_COM, only : lname_strlen,sname_strlen,units_strlen
      IMPLICIT NONE
      SAVE

C**** Dimensions of ice advection grid (EDIT FOR ADVSI GRID CHANGE) are the same as 
C**** atmospheric grid
      INTEGER, parameter :: IMIC=IM, JMIC=JM

C**** Ice advection grid, same as atmospheric grid (CS or latlon)
C**** dimensions IMIC = IM, JMIC = JM
      TYPE(DIST_GRID) :: grid_MIC

C**** variables used in advection (on ICE grid)
!@var RSIX,RSIY first order moments for seaice concentration
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RSIX,RSIY

!@var USIDT,VSIDT sea ice fluxes, saved for advection (m)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: USIDT,VSIDT

C**** Needed for ADVSI (on ATM grid)
!@var RSISAVE saved value of sea ice concentration before DYNSI
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RSISAVE

C**** Ice advection diagnostics
      INTEGER, PARAMETER :: KICIJ=12
!@var IJ_xxx Names for ICIJ diagnostics
      INTEGER IJ_USI,IJ_VSI,IJ_DMUI,IJ_DMVI,IJ_PICE,IJ_MUSI,IJ_MVSI
     *     ,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI,IJ_RSI
!@var ICIJ lat-lon ice dynamic diagnostics (on atm grid)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: ICIJ
!@var lname_icij Long names for ICIJ diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KICIJ) :: LNAME_ICIJ
!@var sname_icij Short names for ICIJ diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KICIJ) :: SNAME_ICIJ
!@var units_icij Units for ICIJ diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KICIJ) :: UNITS_ICIJ
!@var ia_icij IDACC numbers for ICIJ diagnostics
      INTEGER, DIMENSION(KICIJ) :: IA_ICIJ
!@var scale_icij scales for ICIJ diagnostics
      REAL*8, DIMENSION(KICIJ) :: SCALE_ICIJ
!@var [ij]grid_icij Grid descriptors for ICIJ diagnostics
       INTEGER, DIMENSION(KICIJ) :: IGRID_ICIJ,JGRID_ICIJ
!@var denom_icij denominators for ICIJ diagnostics
       INTEGER, DIMENSION(KICIJ) :: DENOM_ICIJ
!@var cdl_icij consolidated metadata for ICIJ output fields in CDL notation
       character(len=100), dimension(kicij*6) :: cdl_icij
#ifdef TRACERS_WATER
!@var KTICIJ number of lat/lon ice dynamic tracer diagnostics
      INTEGER, PARAMETER :: KTICIJ=2
!@var TICIJ  lat/lon ice dynamic tracer diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)  :: TICIJ
!@var ticij_xxx indices for TICIJ diags
      INTEGER :: tICIJ_tusi,tICIJ_tvsi
!@var lname_ticij Long names for TICIJ diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KTICIJ) :: LNAME_TICIJ
!@var sname_ticij Short names for TICIJ diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KTICIJ) :: SNAME_TICIJ
!@var units_ticij Units for TICIJ diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KTICIJ) :: UNITS_TICIJ
!@var ia_ticij IDACC numbers for TICIJ diagnostics
      INTEGER, DIMENSION(KTICIJ) :: IA_TICIJ
!@var scale_ticij scales for TICIJ diagnostics
      REAL*8, DIMENSION(KTICIJ) :: SCALE_TICIJ
!@var ijgrid_ticij Grid descriptor for TICIJ diagnostics
      INTEGER, DIMENSION(KTICIJ) :: IJGRID_TICIJ
#endif

      END MODULE ICEDYN_COM

      SUBROUTINE ALLOC_ICEDYN_COM(grid)
!@sum ALLOC_ICEDYN_COM allocates arrays defined in the ICEDYN_COM module
!@auth Rosalinda de Fainchtein
!     grid = atm grid

      USE DOMAIN_DECOMP_ATM, only : GET,DIST_GRID
      USE MODEL_COM, only : im
      USE ICEDYN_COM, only : grid_MIC,imic
      USE ICEDYN_COM, only : KICIJ
      USE ICEDYN_COM, only : RSIX,RSIY,USIDT,VSIDT,RSISAVE,ICIJ
#ifdef TRACERS_WATER
      USE ICEDYN_COM, only : TICIJ,KTICIJ,NTM
#endif
      IMPLICIT NONE

      LOGICAL, SAVE :: init=.false.
      INTEGER :: I_1H    , I_0H
      INTEGER :: J_1H    , J_0H
      INTEGER :: I_1H_MIC, I_0H_MIC
      INTEGER :: J_1H_MIC, J_0H_MIC
      INTEGER :: IER
      TYPE(DIST_GRID) :: grid

      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

!*** For now set grid_MIC to be the same as grid
!    This is consistent with the current status of the code for
!    parallelization along latitude (j)
      grid_MIC=grid   

C**** Get dimensioning parameters for arrays defined in the grid  and
C**** grid_MIC stencils.

      CALL GET(grid    , I_STRT_HALO=I_0H    , I_STOP_HALO=I_1H    ,
     &     J_STRT_HALO=J_0H    , J_STOP_HALO=J_1H    )
      CALL GET(grid_MIC, I_STRT_HALO=I_0H_MIC, I_STOP_HALO=I_1H_MIC,
     &     J_STRT_HALO=J_0H_MIC, J_STOP_HALO=J_1H_MIC)

      ALLOCATE( RSIX(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC),
     &          RSIY(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC),
     &     STAT = IER)

      ALLOCATE( USIDT(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC),
     &          VSIDT(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC),
     &     STAT = IER)

      ALLOCATE( RSISAVE(I_0H:I_1H, J_0H:J_1H),
     &     STAT = IER)

      ALLOCATE(  ICIJ(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC,KICIJ),
     &     STAT = IER)

#ifdef TRACERS_WATER
      ALLOCATE( TICIJ(I_0H_MIC:I_1H_MIC, J_0H_MIC:J_1H_MIC,KTICIJ, NTM),
     &     STAT = IER)
#endif

      return
      END SUBROUTINE ALLOC_ICEDYN_COM


      SUBROUTINE io_icedyn(kunit,iaction,ioerr)
!@sum  io_icedyn reads and writes dynamic ice arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irsficno,irsficnt
     *     ,irerun,lhead
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT, PACK_DATA, UNPACK_DATA
      USE ICEDYN_COM
      USE ICEDYN, only : grid_ICDYN,imicdyn,jmicdyn,USI,VSI
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICEDYN01"
      REAL*8, DIMENSION(IMIC, JMIC) :: RSIX_glob,RSIY_glob
      REAL*8, DIMENSION(IMICDYN,JMICDYN) :: USI_glob,VSI_glob

      write(MODULE_HEADER(lhead+1:80),'(a7,i3,a1,i3,a)')
     *     'R8 dim(',imic,',',jmic,'):RSIX,RSIY,USI,VSI'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        CALL PACK_DATA(grid_MIC, RSIX, RSIX_GLOB)
        CALL PACK_DATA(grid_MIC, RSIY, RSIY_GLOB)
        CALL PACK_DATA(grid_ICDYN,  USI,  USI_GLOB)    ! TODO: usi/vsi not on atm grid anymore
        CALL PACK_DATA(grid_ICDYN,  VSI,  VSI_GLOB)    ! idem
        IF (AM_I_ROOT())
     &   WRITE (kunit,err=10) MODULE_HEADER,RSIX_glob,RSIY_glob
     &                                     , USI_glob, VSI_glob
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFICNO)           ! initial conditions (no ocean)
        CASE (ioread,irerun,irsfic,irsficnt)    ! restarts
          if (AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,       RSIX_glob,RSIY_glob
     &                                         , USI_glob, VSI_glob
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER,
     &             MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_DATA(grid_MIC, RSIX_GLOB, RSIX)
          CALL UNPACK_DATA(grid_MIC, RSIY_GLOB, RSIY)
          CALL UNPACK_DATA(grid_MIC,  USI_GLOB,  USI)
          CALL UNPACK_DATA(grid_MIC,  VSI_GLOB,  VSI)
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icedyn

      SUBROUTINE io_icdiag(kunit,it,iaction,ioerr)
!@sum  io_icdiag reads and writes ice dynamic diagnostic arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irsfic,irsficnt,irerun,ioread_single,lhead
      USE DOMAIN_DECOMP_ATM, only : GET, AM_I_ROOT
      USE DOMAIN_DECOMP_ATM, only : PACK_DATA, UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : ESMF_BCAST
      USE ICEDYN_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICDIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var ICIJ4 dummy arrays for reading diag. files
      REAL*8, DIMENSION(:,:,:), allocatable  :: ICIJ4
      REAL*4, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ4_GLOB
      REAL*8, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ4_GLOB8
      REAL*8, DIMENSION(IMIC,JMIC,KICIJ)  :: ICIJ_GLOB
      INTEGER :: J_0H_MIC, J_1H_MIC

#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), allocatable  :: TICIJ4
      REAL*4, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ4_GLOB
      REAL*8, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ4_GLOB8
      REAL*8, DIMENSION(IMIC,JMIC,KTICIJ,NTM)  :: TICIJ_GLOB
!@var TR_HEADER Character string label for individual tracer records
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TRICDIAG01"

      write(TR_MODULE_HEADER(lhead+1:80),'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')
     *     'R8 Ticij(',imic,',',jmic,',',kticij,',',ntm,'),it'
#endif
      write(MODULE_HEADER(lhead+1:80),'(a8,i3,a1,i3,a1,i2,a4)')
     *     'R8 ICij(',imic,',',jmic,',',kicij,'),it'

      CALL GET(grid_MIC, J_STRT_HALO=J_0H_MIC, J_STOP_HALO=J_1H_MIC)

      SELECT CASE (IACTION)
      CASE (IOWRITE)  ! output to standard restart file
        CALL PACK_DATA(grid_mic, icij, icij_glob)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) MODULE_HEADER,ICIJ_glob,it
#ifdef TRACERS_WATER
        CALL PACK_DATA(grid_mic, TICIJ, TICIJ_GLOB)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) TR_MODULE_HEADER,TICIJ_GLOB,it
#endif
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        CALL PACK_DATA(grid_mic, icij, icij_glob)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) MODULE_HEADER,REAL(ICIJ_GLOB,KIND=4),it
#ifdef TRACERS_WATER
        TR_MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        CALL PACK_DATA(grid_mic, TICIJ, TICIJ_GLOB)
        IF (AM_I_ROOT())
     &     WRITE (kunit,err=10) TR_MODULE_HEADER
     &          ,REAL(TICIJ_GLOB,KIND=4),it
#endif
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,ICIJ4_GLOB,it
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
          endif
C**** accumulate diagnostics
          allocate( ICIJ4(IMIC, J_0H_MIC:J_1H_MIC, KICIJ) )
          ICIJ4 = 0.d0 ! should do "halo_update" instead?
          ICIJ4_GLOB8 = ICIJ4_GLOB ! convert to real*8
          CALL UNPACK_DATA(grid_MIC, ICIJ4_GLOB8, ICIJ4)
          call ESMF_BCAST(grid_MIC, it)   !MPI_BCAST instead
          ICIJ(:,J_0H_MIC:J_1H_MIC,:)=ICIJ(:,J_0H_MIC:J_1H_MIC,:)
     &                            +ICIJ4(:,J_0H_MIC:J_1H_MIC,:)
          deallocate( ICIJ4 )
#ifdef TRACERS_WATER
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TR_HEADER,TICIJ4_GLOB,it
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TR_HEADER
     *             ,TR_MODULE_HEADER
              GO TO 10
            END IF
          endif
C**** accumulate diagnostics
          allocate( TICIJ4(IMIC, J_0H_MIC:J_1H_MIC, KTICIJ, NTM) )
          TICIJ4 = 0.d0 ! should do "halo_ipdate" instead?
          TICIJ4_GLOB8 = TICIJ4_GLOB ! convert to real*8
          CALL UNPACK_DATA(grid_MIC, TICIJ4_GLOB8, TICIJ4)
          call ESMF_BCAST(grid_MIC, it)
          TICIJ(:,J_0H_MIC:J_1H_MIC,:,:)=TICIJ(:,J_0H_MIC:J_1H_MIC,:,:)
     &                                 +TICIJ4(:,J_0H_MIC:J_1H_MIC,:,:)
          deallocate( TICIJ4 )
#endif
        CASE (ioread)    ! restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,ICIJ_GLOB,it
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_DATA(grid_MIC, ICIJ_GLOB, ICIJ)
          call ESMF_BCAST(grid_MIC, it)
#ifdef TRACERS_WATER
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TR_HEADER,TICIJ_GLOB,it
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TR_HEADER
     *             ,TR_MODULE_HEADER
              GO TO 10
            END IF
          end if
          CALL UNPACK_DATA(grid_MIC, TICIJ_GLOB, TICIJ)
          call ESMF_BCAST(grid_MIC, it)
#endif
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icdiag

#ifdef NEW_IO
      subroutine def_rsf_icedyn(fid)
!@sum  def_rsf_icedyn defines ice dynam array structure in restart files
!@auth M. Kelley
!@ver  beta
      use icedyn, only : usi,vsi
      use icedyn_com, only : grid=>grid_MIC,
     &     rsix,rsiy
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,rsix,'rsix(dist_imic,dist_jmic)')
      call defvar(grid,fid,rsiy,'rsiy(dist_imic,dist_jmic)')
      call defvar(grid,fid,usi,'usi(dist_imic,dist_jmic)')
      call defvar(grid,fid,vsi,'vsi(dist_imic,dist_jmic)')
      return
      end subroutine def_rsf_icedyn

      subroutine new_io_icedyn(fid,iaction)
!@sum  new_io_icedyn read/write ice dynam arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use icedyn, only : usi,vsi
      use icedyn_com, only : grid=>grid_MIC,
     &     rsix,rsiy
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart file
        call write_dist_data(grid, fid, 'rsix', rsix)
        call write_dist_data(grid, fid, 'rsiy', rsiy)
        call write_dist_data(grid, fid, 'usi', usi)
        call write_dist_data(grid, fid, 'vsi', vsi)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'rsix', rsix)
        call read_dist_data(grid, fid, 'rsiy', rsiy)
        call read_dist_data(grid, fid, 'usi', usi)
        call read_dist_data(grid, fid, 'vsi', vsi)
      end select
      return      
      end subroutine new_io_icedyn

      subroutine def_rsf_icdiag(fid,r4_on_disk)
!@sum  def_rsf_icdiag defines ice diag array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use icedyn_com, only : grid=>grid_MIC,
     &     icij
      use icedyn_com
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,icij,'icij(dist_imic,dist_jmic,kicij)',
     &     r4_on_disk=r4_on_disk)
      return
      end subroutine def_rsf_icdiag

      subroutine new_io_icdiag(fid,iaction)
!@sum  new_io_icdiag read/write ice diag arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use icedyn_com, only : grid=>grid_MIC,
     &     icij
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)           ! output to restart or acc file
        call write_dist_data(grid, fid, 'icij', icij)
      case (ioread)            ! input from restart or acc file
        call read_dist_data(grid, fid, 'icij', icij)
      end select
      return      
      end subroutine new_io_icdiag

      subroutine def_meta_icdiag(fid)
!@sum  def_meta_icdiag defines icedyn metadata in acc files
!@auth M. Kelley
!@ver  beta
      use icedyn_com, only : grid=>grid_MIC,
     &     ia_icij,denom_icij,scale_icij,sname_icij,cdl_icij
      use pario, only : defvar,write_attr
      use geom, only : lon2d_dg,lat2d_dg ! TEMPORARY
      implicit none
      integer :: fid         !@var fid file id

      call defvar(grid,fid,lat2d_dg(:,1),'latic(jmic)')
      call defvar(grid,fid,lat2d_dg(:,2),'latic2(jmic)')
      call defvar(grid,fid,lon2d_dg(:,1),'lonic(imic)')
      call defvar(grid,fid,lon2d_dg(:,2),'lonic2(imic)')

      call write_attr(grid,fid,'icij','reduction','sum')
      call write_attr(grid,fid,'icij','split_dim',3)
      call defvar(grid,fid,ia_icij,'ia_icij(kicij)')
      call defvar(grid,fid,denom_icij,'denom_icij(kicij)')
      call defvar(grid,fid,scale_icij,'scale_icij(kicij)')
      call defvar(grid,fid,sname_icij,'sname_icij(sname_strlen,kicij)')
      call defvar(grid,fid,cdl_icij,'cdl_icij(cdl_strlen,kcdl_icij)')

      return
      end subroutine def_meta_icdiag

      subroutine write_meta_icdiag(fid)
!@sum  write_meta_icdiag write icedyn accumulation metadata to file
!@auth M. Kelley
!@ver  beta
      use icedyn_com, only : grid=>grid_MIC,
     &     ia_icij,denom_icij,scale_icij,sname_icij,cdl_icij
      use pario, only : write_dist_data,write_data
      use geom, only : lon2d_dg,lat2d_dg ! TEMPORARY
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'latic',lat2d_dg(:,1))
      call write_data(grid,fid,'latic2',lat2d_dg(:,2))
      call write_data(grid,fid,'lonic',lon2d_dg(:,1))
      call write_data(grid,fid,'lonic2',lon2d_dg(:,2))

      call write_data(grid,fid,'ia_icij',ia_icij)
      call write_data(grid,fid,'denom_icij',denom_icij)
      call write_data(grid,fid,'scale_icij',scale_icij)
      call write_data(grid,fid,'sname_icij',sname_icij)
      call write_data(grid,fid,'cdl_icij',cdl_icij)

      return
      end subroutine write_meta_icdiag

      subroutine set_ioptrs_iceacc_default
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays used during normal operation. 
c temporarily empty.
      return
      end subroutine set_ioptrs_iceacc_default

#endif /* NEW_IO */


      SUBROUTINE reset_icdiag
!@sum reset_icdiag resets ice dynamic diagnostic arrays
!@auth Gavin Schmidt
      USE ICEDYN_COM
      IMPLICIT NONE

      ICIJ=0.
#ifdef TRACERS_WATER
      TICIJ=0.
#endif
      RETURN
      END SUBROUTINE reset_icdiag

      SUBROUTINE DYNSI
!@sum  DYNSI calculates ice velocites
!@+    Note that the ice velocities are calculated on the ice grid
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang), D. Gueyffier (cubed sphere)
!@ver  1.0
      USE CONSTANT, only : rhoi,grav,omega,rhows
      USE MODEL_COM, only : dts=>dtsrc,focean
      USE RESOLUTION, only : aIM=>IM, aJM=>JM
      USE ICEDYN, only : imicdyn,jmicdyn,  !dimensions of icedyn grid
     &     nx1,ny1
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,aGET=>GET
      USE DOMAIN_DECOMP_1D, only : DIST_GRID, iGET=>GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, NORTH, SOUTH
      USE GEOM, only : abyaxyp=>byaxyp,imaxj
#ifndef CUBE_GRID
     &     ,adxyn=>dxyn,adxys=>dxys,adxyv=>dxyv
#endif
      USE ICEDYN, only : dxp,dyv
      USE ICEDYN, only : grid_ICDYN
      USE ICEDYN, only : press,heffm,uvm,dwatn,cor
     *     ,sinwat,coswat,bydts,sinen,uice,vice,heff,area,gairx,gairy
     *     ,gwatx,gwaty,pgfub,pgfvb,amass,uicec,vicec,uib,vib,dmu,dmv
     *     ,usi,vsi
      USE ICEDYN_COM, only : usidt,vsidt,rsisave,icij,ij_usi
     *     ,ij_vsi,ij_dmui,ij_dmvi,ij_pice,ij_rsi
      USE FLUXES, only : dmua,dmva,dmui,dmvi,UI2rho,ogeoza,uosurf,vosurf
     *     ,apress,uisurf,visurf
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      IMPLICIT NONE
      SAVE
      include 'mpif.h'
C**** intermediate calculation for pressure gradient terms
      REAL*8, DIMENSION(IMICDYN, 
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO) ::
     &                            PGFU,PGFV
C****
      real*8, allocatable, dimension(:,:) ::
     &     aPtmp,aheff,aarea,iPtmp,iRSI,iMSI,iFOCEAN,
     &     iUOSURF,iVOSURF,iDMUA,iDMVA,
     &     iDMUI,iDMVI,aDMU,aDMV,aUSI,aVSI

      REAL*8, PARAMETER :: BYRHOI=1D0/RHOI
      REAL*8 :: hemi
      INTEGER I,J,ip1,im1,ierr,mype
      REAL*8 USINP,DMUINP,duA,dvA
      INTEGER :: iJ_1   , iJ_0
      INTEGER :: iJ_1S  , iJ_0S
      INTEGER :: iJ_1H  , iJ_0H
      INTEGER :: iJ_1STG, iJ_0STG
      INTEGER :: aI_1   , aI_0
      INTEGER :: aJ_1   , aJ_0
      INTEGER :: aI_1H  , aI_0H
      INTEGER :: aJ_1H  , aJ_0H
      INTEGER :: aJ_1S  , aJ_0S
      INTEGER :: aJ_1STG, aJ_0STG

C**** Get loop indices  corresponding to grid_ICDYN and atm. grid structures
      CALL iGET(grid_ICDYN, J_STRT=iJ_0       , J_STOP=iJ_1
     &                   , J_STRT_SKP=iJ_0S   , J_STOP_SKP=iJ_1S
     &                   , J_STRT_HALO=iJ_0H  , J_STOP_HALO=iJ_1H )
      call iGET(grid_ICDYN, J_STRT_STGR=iJ_0STG, J_STOP_STGR=iJ_1STG)
      call aGET(agrid    , I_STRT=aI_0        , I_STOP=aI_1     
     &                   , J_STRT=aJ_0        , J_STOP=aJ_1     )
      call aGET(agrid    , I_STRT_HALO=aI_0H  , I_STOP_HALO=aI_1H     
     &                   , J_STRT_HALO=aJ_0H  , J_STOP_HALO=aJ_1H )
      call aGET(agrid    , J_STRT_SKP=aJ_0S   , J_STOP_SKP=aJ_1S)
      call aGET(agrid    , J_STRT_STGR=aJ_0STG, J_STOP_STGR=aJ_1STG)

      allocate(
     &     aPtmp(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     iPtmp(1:IMICDYN,iJ_0H:iJ_1H),
     &     iRSI(1:IMICDYN,iJ_0H:iJ_1H),
     &     iMSI(1:IMICDYN,iJ_0H:iJ_1H),
     &     aHEFF(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aAREA(aI_0H:aI_1H,aJ_0H:aJ_1H),
#ifndef CUBE_GRID 
     &     aDMU(NX1,aJ_0H:aJ_1H),    
     &     aDMV(NX1,aJ_0H:aJ_1H),    
#else
     &     aDMU(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aDMV(aI_0H:aI_1H,aJ_0H:aJ_1H),
#endif
     &     aUSI(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     aVSI(aI_0H:aI_1H,aJ_0H:aJ_1H),
     &     iFOCEAN(1:IMICDYN,iJ_0H:iJ_1H),
     &     iUOSURF(1:IMICDYN,iJ_0H:iJ_1H),
     &     iVOSURF(1:IMICDYN,iJ_0H:iJ_1H),
     &     iDMUA(1:IMICDYN,iJ_0H:iJ_1H),
     &     iDMVA(1:IMICDYN,iJ_0H:iJ_1H),
     &     iDMUI(1:IMICDYN,iJ_0H:iJ_1H),
     &     iDMVI(1:IMICDYN,iJ_0H:iJ_1H)
     &     )

C**** Start main loop
C**** Replicate polar boxes
      if (agrid%HAVE_NORTH_POLE) then
        RSI(2:aIM,aJM)=RSI(1,aJM)
        MSI(2:aIM,aJM)=MSI(1,aJM)
        DMUA(2:aIM,aJM,2) = DMUA(1,aJM,2)
        DMVA(2:aIM,aJM,2) = DMVA(1,aJM,2)
      end if

C**** save current value of sea ice concentration for ADVSI
C**** RSISAVE is on atmospheric grid
      RSISAVE(:,:)=RSI(:,:)

C**** Pressure anomaly at surface APRESS: calculated by sea ice routines
C**** APRESS is on atmospheric grid. We are now no longer using this as
C**** a forcing in the sea ice dynamics because it is partially
C**** included already in the internal ice pressure gradient. The
C**** atmospheric pressure gradient term does not in general produce a
C**** horizontal force in a solid (such as ice).

C**** calculate sea surface tilt on atmospheric C grid
C**** (using OGEOZA on atmospheric grid plus displacement of free
C**** surface due to presence of ice). This is ignored in favour of
C**** geostrophy if osurf_tilt=0.
C**** PGF is an accelaration

      if (grid_ICDYN%HAVE_NORTH_POLE) PGFU(1:IMICDYN,JMICDYN)=0
      if (grid_ICDYN%HAVE_SOUTH_POLE) PGFU(1:IMICDYN, 1)=0  !RKF
C****  define scalar pressure on atm grid then regrid it to the icedyn grid  
      DO J=aJ_0,aJ_1            
         DO I=aI_0,aI_1           
             aPtmp(I,J)=(OGEOZA(I,J)
     *            +(RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I))*GRAV/RHOWS) 
        END DO
      END DO

      CALL HALO_UPDATE(agrid, FOCEAN, from=NORTH )
      CALL HALO_UPDATE(agrid, RSI   , from=NORTH )
      CALL HALO_UPDATE(agrid, MSI   , from=NORTH )

      call INT_AtmA2IceA(aPtmp,iPtmp)   
      call INT_AtmA2IceA(Focean,iFocean)   
      call INT_AtmA2IceA(RSI,iRSI)

      CALL HALO_UPDATE(grid_ICDYN, iPtmp , from=NORTH )
      CALL HALO_UPDATE(grid_ICDYN, iFOCEAN, from=NORTH )
      CALL HALO_UPDATE(grid_ICDYN, iRSI   , from=NORTH )

c*** Calculate gradient on ice dyn. grid
      DO J=iJ_0S,iJ_1S
        I=IMICDYN
        DO IP1=1,IMICDYN
          IF(iFOCEAN(I,J).gt.0 .and. iFOCEAN(IP1,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(IP1,J).gt.0.) THEN
            PGFU(I,J)=-(iPtmp(IP1,J)-iPtmp(I,J))/DXP(J)
          ELSE
            PGFU(I,J)=0.
          END IF
          I=IP1
        END DO
      END DO

      DO J=iJ_0,iJ_1S
        DO I=1,IMICDYN
          IF(iFOCEAN(I,J+1).gt.0 .and. iFOCEAN(I,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(I,J+1).gt.0.) THEN
            PGFV(I,J)=-(iPtmp(I,J+1)-iPtmp(I,J))/DYV(J+1)
          ELSE
            PGFV(I,J)=0.
          END IF
        END DO
      END DO
c*

C**** DMUA is defined over the whole box (not just over ptype)
C**** Convert to stress over ice fraction only (on atmospheric grid)
      DO J=aJ_0,aJ_1
        do i=aI_0,aI_1 
          IF (FOCEAN(I,J)*RSI(I,J).gt.0) THEN
            DMUA(I,J,2) = DMUA(I,J,2)/(FOCEAN(I,J)*RSI(I,J))
            DMVA(I,J,2) = DMVA(I,J,2)/(FOCEAN(I,J)*RSI(I,J))
          ELSE
            DMUA(I,J,2) = 0.
            DMVA(I,J,2) = 0.
          END IF
        END DO
      END DO

      call INT_AtmA2IceA(MSI,iMSI)   

C**** Set up ice grid variables
C**** HEFF,AREA on primary (tracer) grid for ice
      do j=iJ_0,iJ_1
         do i=2,NX1-1
            HEFF(I,J)=iRSI(I-1,J)*(ACE1I+iMSI(I-1,J))*BYRHOI
            HEFF(I,J)=HEFF(I,J)*HEFFM(I,J)
            AREA(I,J)=iRSI(I-1,J)
        enddo
      enddo

C**** fill in overlap regions
      DO J=iJ_0,iJ_1
        HEFF(1,J)=HEFF(NX1-1,J)
        AREA(1,J)=AREA(NX1-1,J)
        HEFF(NX1,J)=HEFF(2,J)
        AREA(NX1,J)=AREA(2,J)
      END DO

C****
C**** Set up mass per unit area and coriolis term (on ice grid B)
C****
C**** Update halo for HEFF
      CALL HALO_UPDATE(grid_ICDYN, HEFF, from=NORTH    )

      DO J=iJ_0,iJ_1S
      DO I=1,NX1-1
        AMASS(I,J)=RHOI*0.25*(HEFF(I,J)
     *       +HEFF(I+1,J)+HEFF(I,J+1)+HEFF(I+1,J+1))
        COR(I,J)=AMASS(I,J)*2.0*OMEGA*SINEN(I,J)
      END DO
      END DO
c**** set north pole
      if (grid_ICDYN%HAVE_NORTH_POLE) then
        do i=1,nx1
          AMASS(i,jmicdyn)= RHOI*HEFF(1,JMICDYN)
          COR  (i,jmicdyn)= AMASS(i,jmicdyn)
     *         *2.0*OMEGA*SINEN(1,JMICDYN)
        end do
      end if                    !end NORTH_POLE block if

c**** interpolate air, current and ice velocity from C grid to B grid
C**** This should be more generally from ocean grid to ice grid
C**** NOTE: UOSURF, VOSURF are expected to be on the A-grid

      CALL HALO_UPDATE(agrid, UOSURF )
      CALL HALO_UPDATE(agrid, VOSURF )
c**** getting instance of (UOSURF, VOSURF) on the icedyn grid
      call INT_AtmA2IceA(UOSURF,iUOSURF)   !already in latlon basis
      call INT_AtmA2IceA(VOSURF,iVOSURF)   !already in latlon basis

C**** Update halo for USI,UOSURF,VOSURF,PGFU
      CALL HALO_UPDATE(grid_ICDYN, USI   , from=NORTH    )
      CALL HALO_UPDATE(grid_ICDYN, iUOSURF , from=NORTH    )
      CALL HALO_UPDATE(grid_ICDYN, iVOSURF , from=NORTH    )
      CALL HALO_UPDATE(grid_ICDYN, PGFU  , from=NORTH    )

      do j=iJ_0,iJ_1S
        im1=imicdyn
        do i=1,imicdyn
          UIB  (i,j)=0.5*(USI (im1,j)  +USI (im1,j+1))   ! iceC--> iceB
          GWATX(i,j)=0.25*(iUOSURF(im1,j)  +iUOSURF(im1,j+1)
     &         +iUOSURF(i,j)+iUOSURF(i,j+1))                     ! ocean -> iceB  
          PGFUB(i,j)=0.5*(PGFU(im1,j)  +PGFU(im1,j+1))   ! iceC--> iceB 
          VIB  (i,j)=0.5*(VSI (im1,j)  +VSI (i,j))       ! y component
          GWATY(i,j)=0.25*(iVOSURF(im1,j)  +iVOSURF(im1,j+1)
     &         +iVOSURF(i,j)+iVOSURF(i,j+1))                     ! y component
          PGFVB(i,j)=0.5*(PGFV(im1,j)  +PGFV(i,j))       ! y component
          im1=i
        enddo
      enddo
c**** set north pole
      if (grid_ICDYN%HAVE_NORTH_POLE) then
        do i=1,imicdyn
          UIB  (i,jmicdyn)=USI(1,jmicdyn)
          GWATX(i,jmicdyn)=0.
          PGFUB(i,jmicdyn)=0.
          VIB  (i,jmicdyn)=0.
          GWATY(i,jmicdyn)=0.
          PGFVB(i,jmicdyn)=0.
        enddo
      end if                    !end NORTH_POLE block if
      DO J=iJ_0,iJ_1
        UIB(nx1-1,J)=UIB(1,J)
        VIB(nx1-1,J)=VIB(1,J)
        GWATX(nx1-1,J)=GWATX(1,J)
        GWATY(nx1-1,J)=GWATY(1,J)
        PGFUB(nx1-1,J)=PGFUB(1,J)
        PGFVB(nx1-1,J)=PGFVB(1,J)
        UIB(nx1,J)=UIB(2,J)
        VIB(nx1,J)=VIB(2,J)
        GWATX(nx1,J)=GWATX(2,J)
        GWATY(nx1,J)=GWATY(2,J)
        PGFUB(nx1,J)=PGFUB(2,J)
        PGFVB(nx1,J)=PGFVB(2,J)
      END DO

c**** interpolate air stress from A grid in atmos, to B grid in ice
C**** change of unit from change of momentum, to flux
C**** Update halo for USI,UOSURF,PGFU
      CALL HALO_UPDATE(agrid, DMUA  )	 
      CALL HALO_UPDATE(agrid, DMVA  )

c**** getting instance of (DMUA, DVMA) on the icedyn grid
      call INT_AtmA2IceA(DMUA(:,:,2),iDMUA)   !stays on latlon basis
      call INT_AtmA2IceA(DMVA(:,:,2),iDMVA)   !stays on latlon basis

      CALL HALO_UPDATE(grid_ICDYN, iDMUA  , from=NORTH    )  
      CALL HALO_UPDATE(grid_ICDYN, iDMVA  , from=NORTH    )

      do j=iJ_0,iJ_1S
        im1=imicdyn
        do i=1,imicdyn
          GAIRX(i,j)=0.25*(idmua(i,j)+idmua(im1,j)+idmua(im1,j+1)
     &         +idmua(i,j+1))*bydts  
          GAIRY(i,j)=0.25*(idmva(i,j)+idmva(im1,j)+idmva(im1,j+1)
     &         +idmva(i,j+1))*bydts  
          im1=i
        enddo
      enddo
      IF (grid_ICDYN%HAVE_NORTH_POLE) THEN
        GAIRX(1:nx1,jmicdyn)=idmua(1,jmicdyn)*bydts
        GAIRY(1:nx1,jmicdyn)=idmva(1,jmicdyn)*bydts
      END IF
      do j=iJ_0,iJ_1S
       GAIRX(nx1-1,j)=GAIRX(1,j)
       GAIRY(nx1-1,j)=GAIRY(1,j)
       GAIRX(nx1,j)=GAIRX(2,j)
       GAIRY(nx1,j)=GAIRY(2,j)
      enddo

c**** read in sea ice velocity
      DO J=iJ_0,iJ_1
      DO I=1,NX1
       UICE(I,J,1)=UIB(I,J)
       VICE(I,J,1)=VIB(I,J)
       UICE(I,J,2)=0.
       VICE(I,J,2)=0.
       UICE(I,J,3)=0.
       VICE(I,J,3)=0.
      END DO
      END DO


c      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
c      write(220+myPE,*) HEFFM
c      write(230+myPE,*) COR
c      write(240+myPE,*) pgfub
c      write(250+myPE,*) pgfvb
c      write(260+myPE,*) gairx 
c      write(270+myPE,*) gairy 
c      write(280+myPE,*) gwatx 
c      write(290+myPE,*) gwaty 
c      write(300+myPE,*) uib 
c      write(310+myPE,*) vib
c      write(320+myPE,*) HEFF
c      write(330+myPE,*) amass

C**** do the looping over pseudo-timesteps
      CALL VPICEDYN

C**** Calculate stress on ice velocity grid (B grid)
      DO J=iJ_0,iJ_1
        hemi=-1.
        if (J.gt.NY1/2) hemi=1.
        DO I=1,NX1-1
          DMU(i,j)=DTS*dwatn(i,j)*(COSWAT*(UICE(i,j,1)-GWATX(i,j))-
     *         HEMI*SINWAT*(VICE(i,j,1)-GWATY(i,j)))
          DMV(i,j)=DTS*dwatn(i,j)*(HEMI*SINWAT*(UICE(i,j,1)-GWATX(i,j))
     *         +COSWAT*(VICE(i,j,1)-GWATY(i,j)))
        END DO
      END DO
      DO J=iJ_0,iJ_1
        DMU(1,J)=DMU(NX1-1,J)
        DMU(NX1,J)=DMU(2,J)
      END DO

C**** interpolate ice velocity and stress from B grid to C grid in atm
C**** Update halos for UICE and DMU
      CALL HALO_UPDATE(grid_ICDYN,  UICE, from=SOUTH     )
      CALL HALO_UPDATE(grid_ICDYN,   DMU, from=SOUTH     )
 
      do j=iJ_0STG,iJ_1STG
        i=imicdyn
        do ip1=1,imicdyn
          usi(i,j)=0.5*(uice(i+1,j-1,1)+uice(i+1,j,1))
          IF (abs(USI(I,J)).lt.1d-10) USI(I,J)=0
          iDMUI(I,J) = 0.5*(dmu(i+1,j-1)+dmu(i+1,j))
          i=ip1
        enddo
      enddo

C**** Update halos for FOCEAN, DXYS, DXYV
      CALL HALO_UPDATE(agrid,  FOCEAN, from=NORTH     )

      do j=iJ_0,iJ_1s
        do i=1,imicdyn
          vsi(i,j)=0.5*(vice(i,j,1)+vice(i+1,j,1))
          IF (abs(VSI(I,J)).lt.1d-10) VSI(I,J)=0
          iDMVI(I,J) = 0.5*(dmv(i,j)+dmv(i+1,j))
       enddo
      enddo

      call INT_IceC2AtmC_U(iDMUI,DMUI) !resampled on atm grid but still expressed in latlon basis
      call INT_IceC2AtmC_V(iDMVI,DMVI) !resampled on atm grid but still expressed in latlon basis

C**** Update halos for FOCEAN, and RSI
      CALL HALO_UPDATE(agrid,  FOCEAN)
      CALL HALO_UPDATE(agrid,  RSI )

      do j=aJ_0STG,aJ_1STG
        do i=aI_0,aI_1
           ip1=i+1
           if (ip1 .eq. aIM+1) ip1=1
C**** Rescale DMUI to be net momentum into ocean
          DMUI(I,J) = 0.5*DMUI(I,J)*(FOCEAN(I,J)*RSI(I,J)+FOCEAN(ip1,J)
     *         *RSI(ip1,J))
        enddo
      enddo

      do j=aJ_0,aJ_1s
        do i=aI_0,aI_1
C**** Rescale DMVI to be net momentum into ocean
          IF (J.lt.aJM-1) DMVI(I,J) = 
     &         0.5*DMVI(I,J)
#ifdef CUBE_GRID
     &          *0.5*(FOCEAN(I,J)*RSI(I,J)
     &          +FOCEAN(I,J+1)*RSI(I,J+1))
#else
     &          *(FOCEAN(I,J)*RSI(I,J)*aDXYN(J)
     &          +FOCEAN(I,J+1)*RSI(I,J+1)*aDXYS(J+1))
     &          /aDXYV(J+1)
#endif
          IF (J.eq.aJM-1) DMVI(I,aJM-1) = 
     &         0.5*DMVI(I,aJM-1)
#ifdef CUBE_GRID
     &         *0.5*(FOCEAN(I,aJM-1)*RSI(I,aJM-1)
     &         +FOCEAN(1,aJM)*RSI(1,aJM))
#else
     &         *(FOCEAN(I,aJM-1)*RSI(I,aJM-1)*aDXYN(aJM-1)
     &         +FOCEAN(1,aJM)*RSI(1,aJM)*aDXYS(aJM))
     &         /aDXYV(aJM)
#endif
        enddo
      enddo

      if (grid_ICDYN%HAVE_SOUTH_POLE) then
        usi(1:imicdyn,1)=0.
      endif
      if (grid_ICDYN%HAVE_NORTH_POLE) then
        vsi(1:imicdyn,jmicdyn)=0.
      endif
      IF (agrid%HAVE_SOUTH_POLE) then
        dmui(1:aim,1)=0.
      END IF
      IF (agrid%HAVE_NORTH_POLE) then
        dmvi(1:aim,ajm)=0.
      END IF

C**** Calculate ustar*2*rho for ice-ocean fluxes on atmosphere grid
C**** UI2rho = | tau |
C**** Update halos for DMU, and DMV
      CALL HALO_UPDATE(grid_ICDYN,  DMU   , from=SOUTH     )
      CALL HALO_UPDATE(grid_ICDYN,  DMV   , from=SOUTH     )

c**** resampled on atm grid (CS or latlon) but still expressed in latlon basis
      call INT_IceB2AtmB_NX1(DMU,aDMU) 
      call INT_IceB2AtmB_NX1(DMV,aDMV)

      CALL HALO_UPDATE(agrid,  aDMU  )
      CALL HALO_UPDATE(agrid,  aDMV  )
      
      do j=aJ_0,aJ_1
         do i=aI_0,imaxj(j)
            UI2rho(i,j)=0
            if (FOCEAN(I,J)*RSI(i,j).gt.0) THEN
#ifdef CUBE_GRID
C**** 4 points average on cubed sphere
               duA = 0.25*(admu(i+1,j)+admu(i,j)
     *              +admu(i+1,j-1)+admu(i,j-1) )
               dvA = 0.25*(admv(i+1,j)+admv(i,j)
     *              +admv(i+1,j-1)+admv(i,j-1) )
#else 
C**** calculate 4 point average of B grid values of stresses
               duA = 0.5*(aDXYN(J)*(admu(i+1,j)+admu(i,j))
     *              +aDXYS(j)*(admu(i+1
     *              ,j-1)+admu(i,j-1)))*aBYAXYP(I,J)
               dvA = 0.5*(aDXYN(J)*(admv(i+1,j)+admv(i,j))
     *              +aDXYS(j)*(admv(i+1
     *              ,j-1)+admv(i,j-1)))*aBYAXYP(I,J)
#endif
               UI2rho(i,j)= sqrt (duA**2 + dvA**2) * bydts
            end if
         end do
      end do

C**** set north pole 
      IF (grid_ICDYN%HAVE_NORTH_POLE) THEN
        USINP=0.
        do i=1,imicdyn
          USINP = USINP + USI(i,jmicdyn)
        enddo
        USINP=USINP/IMICDYN  
        USI(1:imicdyn,jmicdyn)=USINP
        VSI(1:imicdyn,jmicdyn)=0.
      END IF
C**** dmui,dmvi on atm C grid
      IF (agrid%HAVE_NORTH_POLE) THEN
        DMUINP=0.
        do i=1,aim
          DMUINP = DMUINP + DMUI(i,ajm)
        enddo
        DMUINP=DMUINP/aIM
        DMUI(1:aim,ajm)=DMUINP
        DMVI(1:aim,ajm)=0.
      END IF


C**** calculate mass fluxes for the ice advection
C**** Update halos for FOCEAN, and RSI
      CALL HALO_UPDATE(agrid,  FOCEAN, from=NORTH     )
      CALL HALO_UPDATE(agrid,  RSI   , from=NORTH     )

      CALL HALO_UPDATE(grid_ICDYN,  iFOCEAN, from=NORTH     )
      CALL HALO_UPDATE(grid_ICDYN,  iRSI   , from=NORTH     )

c*** interpolate ice velocities to the atm C-grid but 
c*** keep latlon orientation 
      call Int_IceC2AtmC_U(USI,aUSI)
      call Int_IceC2AtmC_V(VSI,aVSI)

c*** Computation of usidt/vsidt: note that usidt/vsidt is seen only by the latlon version of ADVSI 
      DO J=aJ_0,aJ_1S
        DO I=aI_0,aI_1
         ip1=i+1
         if (ip1 .eq. aIM+1) ip1=1
          USIDT(I,J)=0.
          VSIDT(I,J)=0.
          IF (FOCEAN(I,J).gt.0 .and. FOCEAN(IP1,J).gt.0. .and.
     &         RSI(I,J)+RSI(IP1,J).gt.1d-4) 
     &       USIDT(I,J)=aUSI(I,J)*DTS
          IF (FOCEAN(I,J+1).gt.0 .and. FOCEAN(I,J).gt.0. .and.
     &         RSI(I,J)+RSI(I,J+1).gt.1d-4) 
     &       VSIDT(I,J)=aVSI(I,J)*DTS
        END DO
      END DO
      IF (agrid%HAVE_NORTH_POLE) THEN
        VSIDT(1:aIM,aJM)=0.
        USIDT(1:aIM,aJM)=aUSI(1,aJM)*DTS
      END IF


c*** diagnostics
      DO J=iJ_0,iJ_1S
        DO I=1,imicdyn
         ip1=i+1
         if (ip1 .eq. IMICDYN+1) ip1=1
          IF (iFOCEAN(I,J).gt.0 .and. iFOCEAN(IP1,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(IP1,J).gt.1d-4) THEN
            ICIJ(I,J,IJ_USI) =ICIJ(I,J,IJ_USI) +(iRSI(I,J)+iRSI(IP1,J))
     *           *USI(i,j)
            ICIJ(I,J,IJ_DMUI)=ICIJ(I,J,IJ_DMUI)+DMUI(i,j)
          END IF
          IF (iFOCEAN(I,J+1).gt.0 .and. iFOCEAN(I,J).gt.0. .and.
     *         iRSI(I,J)+iRSI(I,J+1).gt.1d-4) THEN
            ICIJ(I,J,IJ_VSI) =ICIJ(I,J,IJ_VSI) +(iRSI(I,J)+iRSI(I,J+1))
     *           *VSI(i,j)
            ICIJ(I,J,IJ_DMVI)=ICIJ(I,J,IJ_DMVI)+DMVI(i,j)
          END IF
          ICIJ(I,J,IJ_PICE)=ICIJ(I,J,IJ_PICE)+ iRSI(I,J)*press(i+1,j)
          ICIJ(I,J,IJ_RSI) =ICIJ(I,J,IJ_RSI) + iRSI(I,J)
        END DO
      END DO
      IF (agrid%HAVE_NORTH_POLE) THEN
        ICIJ(1,JMICDYN,IJ_USI) =ICIJ(1,JMICDYN,IJ_USI) 
     &       +iRSI(1,JMICDYN)*USI(1,JMICDYN)
        ICIJ(1,JMICDYN,IJ_DMUI)=ICIJ(1,JMICDYN,IJ_DMUI)+DMUI(1,JMICDYN)
        ICIJ(1,JMICDYN,IJ_RSI) =ICIJ(1,JMICDYN,IJ_RSI) +iRSI(1,JMICDYN)
        ICIJ(1,JMICDYN,IJ_PICE)=ICIJ(1,JMICDYN,IJ_PICE)
     &       +iRSI(1,JMICDYN)*press(1,JMICDYN)
      END IF


C**** Set uisurf,visurf (on atm A grid) for use in atmos. drag calc.
      call get_uisurf(aUSI,aVSI,uisurf,visurf) !uisurf/visurf are on atm grid but are latlon oriented
C****

c      do i=1,aim
c         do j=aj_0,aj_1s
c            write(340+myPE,*) i,j,DMUI(i,j)
c         enddo
c      enddo
c      do i=1,aim
c         do j=aj_0,aj_1s
c            write(350+myPE,*) i,j,DMVI(i,j)
c         enddo
c      enddo
c      write(360+myPE,*) DMU
c      write(370+myPE,*) DMV
c      write(380+mype,*) usi
c      write(390+mype,*) vsi
c      write(400+mype,*) usidt
c      write(410+mype,*) vsidt
c      write(420+mype,*) uisurf
c      write(430+mype,*) visurf
c      write(440+mype,*) ui2rho

c      write(800+mype,*) rsi
c      write(900+mype,*) focean

      deallocate(aPtmp,aheff,aarea,iPtmp,iRSI,iMSI,iFOCEAN,
     &     iUOSURF,iVOSURF,iDMUA,iDMVA,aDMU,aDMV,aUSI,aVSI,
     &     iDMUI,iDMVI)

      RETURN
      END SUBROUTINE DYNSI

#ifndef CUBE_GRID
      SUBROUTINE ADVSI
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
      USE CONSTANT, only : grav,tf
      USE MODEL_COM, only :  im,jm,focean,p,ptop,kocean
      USE DOMAIN_DECOMP_ATM, only : grid, GET
      USE DOMAIN_DECOMP_ATM, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : SOUTH, NORTH
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE_COLUMN
      USE GEOM, only : dxyp,dyp,dxp,dxv,bydxyp,imaxj   !atmosphere grid geom
      USE ICEDYN_COM, only : usidt,vsidt,rsix,rsiy,rsisave,icij,ij_musi
     *     ,ij_mvsi,ij_husi,ij_hvsi,ij_susi,ij_svsi
#ifdef TRACERS_WATER
     *     ,ticij,ticij_tusi,ticij_tvsi
#endif
      USE ICEDYN_COM, only : grid_MIC
      USE SEAICE, only : ace1i,xsi,Ti,Ei
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,lmi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,apress,msicnv,fwsim,gtempr
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DIAG_COM, only : oa
      IMPLICIT NONE
!      REAL*8, DIMENSION(IM) :: FAW,FASI,FXSI,FYSI
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     FASI, FXSI, FYSI, FAW
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=2+2*LMI
#else
      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
      INTEGER ITR
      REAL*8 TRSNOW(NTM), TRICE(NTM)
#endif
      REAL*8 FMSI(NTRICE,IM),SFMSI(NTRICE),AMSI(NTRICE)
      REAL*8 FMSJ(IM,NTRICE,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 BYFOA(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cmsi1(grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cmsi2(grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cssi1(grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,cssi2(grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER I,J,L,IM1,IP1,K
      REAL*8 SFASI,DMHSI,ASI,YRSI,XRSI,FRSI,SICE,TMP,TICE,ENRG
!@var MHS mass/heat/salt content of sea ice
      REAL*8 MHS(NTRICE,IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
C****
C**** FLUXCB  USIDT  U compon of time integrated sea ice velocity (m)
C****         VSIDT  V compon of time integrated sea ice velocity (m)
C****
C****         FAW    flux of surface water area (m^2) = USIDT*DYP
C****         FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C****         FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)

      INTEGER J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
C**** Get grid parameters
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S ,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO =J_1H ,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Regularise ice concentration gradients to prevent advection errors
      DO J=J_0S, J_1S
      DO I=1,IM
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)

          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
      END DO
      END DO

C**** update RSISAVE for diagnostics
      RSISAVE(:,:)=RSI(:,:)

C**** set up local MHS array to contain all advected quantities
C**** MHS(1:2) = MASS, MHS(3:2+LMI) = HEAT, MHS(3+LMI:2+2*LMI)=SALT
C**** Currently this is on atmospheric grid
      MHS(1,:,J_0:J_1) = ACE1I + SNOWI(:,J_0:J_1)
      MHS(2,:,J_0:J_1) = MSI(:,J_0:J_1)
      DO L=1,LMI
        MHS(L+2,:,J_0:J_1) = HSI(L,:,J_0:J_1)
        MHS(L+2+LMI,:,J_0:J_1) = SSI(L,:,J_0:J_1)
      END DO
#ifdef TRACERS_WATER
C**** add tracers to advected arrays
      DO J=J_0, J_1
        DO I=1,IM
          DO ITR=1,NTM
          IF (SNOWI(I,J)*XSI(2).gt.XSI(1)*ACE1I) THEN ! layer 1:all snow
            SICE=SSI(1,I,J)+SSI(2,I,J)
            TRSNOW(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J)*MAX(1.
     *           -(ACE1I-SICE)/(XSI(2)*(ACE1I+SNOWI(I,J))-SICE),0d0)
          ELSE                  ! first layer is snow and some ice
            TRSNOW(ITR) = TRSI(ITR,1,I,J)*MIN(SNOWI(I,J)/(XSI(1)*(ACE1I
     *           +SNOWI(I,J))-SSI(1,I,J)),1d0)
          END IF
          TRICE(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J) - TRSNOW(ITR)
          MHS(1+2+(1+ITR)*LMI,I,J)=TRSNOW(ITR)
          MHS(2+2+(1+ITR)*LMI,I,J)=TRICE(ITR)
          DO L=3,LMI
            MHS(L+2+(1+ITR)*LMI,I,J)=TRSI(ITR,L,I,J)
          END DO
          END DO
        END DO
      END DO
#endif
C**** define inverse area array
      DO J=J_0, J_1
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0) THEN
          BYFOA(I,J)=BYDXYP(J)/FOCEAN(I,J)
        ELSE
          BYFOA(I,J)=0.
        END IF
      END DO
      END DO
C****
C**** North-South Advection of Sea Ice
C****
      SFASI  = 0.
      SFMSI(1:NTRICE) = 0.
C**** Update halo of DXV,RSIY,RSI,RSIX,FOCEAN,BYDXYP,and MHS
      CALL HALO_UPDATE(grid, RSI  , FROM=NORTH)
      CALL HALO_UPDATE(grid,FOCEAN, FROM=NORTH+SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, MHS  , FROM=NORTH)

      CALL HALO_UPDATE(grid, VSIDT, FROM=SOUTH)

      CALL HALO_UPDATE(grid_MIC, RSIY, FROM=NORTH)
      CALL HALO_UPDATE(grid_MIC, RSIX, FROM=NORTH)

      DO I=1,IM
C****
C**** Calculate south-north sea ice fluxes at grid box edges
C****
      DO 120 J=J_0S,MIN(JM-2,J_1)
      IF(VSIDT(I,J).eq.0.)  GO TO 120
      FAW(I,J) = VSIDT(I,J)*DXV(J+1) ! be careful with atm. grid index
      IF(VSIDT(I,J).le.0.) THEN
C**** Sea ice velocity is southward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J+1)-
     *       (1d0+FAW(I,J)*BYDXYP(J+1))*RSIY(I,J+1))*FOCEAN(I,J+1)
        FXSI(I,J)=FAW(I,J)*RSIX(I,J+1)*FOCEAN(I,J+1)
        FYSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J+1)*
     *       FAW(I,J)*RSIY(I,J+1)*FOCEAN(I,J+1) - 3d0*FASI(I,J))
        FMSJ(I,1:NTRICE,J) = FASI(I,J)*MHS(1:NTRICE,I,J+1)
      ELSE
C**** Sea ice velocity is northward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J)+(1d0-FAW(I,J)*BYDXYP(J))*RSIY(I,J))
     *       *FOCEAN(I,J)
        FXSI(I,J)=FAW(I,J)*RSIX(I,J)*FOCEAN(I,J)
        FYSI(I,J)=FAW(I,J)*
     *       (FAW(I,J)*BYDXYP(J)*FAW(I,J)*RSIY(I,J)*FOCEAN(I,J)
     *       -3d0*FASI(I,J))
        FMSJ(I,1:NTRICE,J) = FASI(I,J)*MHS(1:NTRICE,I,J)
      END IF
        ICIJ(I,J,IJ_MVSI)=ICIJ(I,J,IJ_MVSI)+SUM(FMSJ(I,1:2,J))
        ICIJ(I,J,IJ_HVSI)=ICIJ(I,J,IJ_HVSI)+SUM(FMSJ(I,3:2+LMI,J))
        ICIJ(I,J,IJ_SVSI)=ICIJ(I,J,IJ_SVSI)+SUM(FMSJ(I,3+LMI:2+2*LMI,J))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,J,TICIJ_TVSI,ITR)=TICIJ(I,J,TICIJ_TVSI,ITR)+
     *          SUM(FMSJ(I,3+(1+ITR)*LMI:2+(2+ITR)*LMI,J))
         END DO
#endif
  120 CONTINUE
C****
C**** Calculate south-north sea ice fluxes near North Pole
C****
      IF (HAVE_NORTH_POLE) THEN
        IF(VSIDT(I,JM-1).eq.0.) cycle
        FAW(I,JM-1) = VSIDT(I,JM-1)*DXV(JM) ! careful with atm.grid index!
        IF(VSIDT(I,JM-1).le.0.) THEN
C**** Sea ice velocity is southward from North Pole box
          FASI(I,JM-1) = FAW(I,JM-1)*RSI(1,JM)*FOCEAN(1,JM)
          FXSI(I,JM-1) = 0.
          FYSI(I,JM-1) = -FAW(I,JM-1)*FASI(I,JM-1)
          FMSJ(I,1:NTRICE,JM-1) = FASI(I,JM-1)*MHS(1:NTRICE,1,JM)
        ELSE
C**** Sea ice velocity is northward into North Pole box
          FASI(I,JM-1) = FAW(I,JM-1)*FOCEAN(I,JM-1)*
     *       (RSI(I,JM-1)+(1d0-FAW(I,JM-1)*BYDXYP(JM-1))*RSIY(I,JM-1))
          FXSI(I,JM-1) = FAW(I,JM-1)*RSIX(I,JM-1)*FOCEAN(I,JM-1)
          FYSI(I,JM-1) = FAW(I,JM-1)*(FAW(I,JM-1)*BYDXYP(JM-1)*
     *       FAW(I,JM-1)*RSIY(I,JM-1)*FOCEAN(I,JM-1)-3d0*FASI(I,JM-1))
          FMSJ(I,1:NTRICE,JM-1) = FASI(I,JM-1)*MHS(1:NTRICE,I,JM-1)
        END IF
C**** Accumulate sea ice leaving and entering North Pole box
        SFASI = SFASI + FASI(I,JM-1)
        SFMSI(1:NTRICE) = SFMSI(1:NTRICE) + FMSJ(I,1:NTRICE,JM-1)
         ICIJ(I,JM-1,IJ_MVSI)=ICIJ(I,JM-1,IJ_MVSI)+SUM(FMSJ(I,1:2,JM-1))
         ICIJ(I,JM-1,IJ_HVSI)=ICIJ(I,JM-1,IJ_HVSI)
     &                        +SUM(FMSJ(I,3:2+LMI,JM-1))
         ICIJ(I,JM-1,IJ_SVSI)=ICIJ(I,JM-1,IJ_SVSI)+
     *      SUM(FMSJ(I,3+LMI:2+2*LMI,JM-1))
#ifdef TRACERS_WATER
           DO ITR=1,NTM
             TICIJ(I,JM-1,TICIJ_TVSI,ITR)=TICIJ(I,JM-1,TICIJ_TVSI,ITR)+
     *            SUM(FMSJ(I,3+(1+ITR)*LMI:2+(2+ITR)*LMI,JM-1))
           END DO
#endif
      ENDIF

      END DO ! I loop
C****
C**** Update sea ice variables due to south-north fluxes
C****

C****Update halo of VSIDT, FASI, FAW, FOCEAN, FMSI, and FXSI
      CALL HALO_UPDATE_COLUMN(grid, FMSJ, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FASI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FXSI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FYSI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, FAW,  FROM=SOUTH)

      DO I=1,IM
      DO 330 J=J_0S, J_1S
      IF(VSIDT(I,J-1)) 240,210,280
C**** VSIDT(J-1)=0.
  210 IF(VSIDT(I,J))  220,330,230
C**** VSIDT(J-1)=0, VSIDT(J)<0.
  220 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I,J)
      DO 225 K=1,NTRICE
  225 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSJ(I,K,J)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FYSI(I,J)
     *  + 3d0*(FAW(I,J)*ASI-DXYP(J)*FASI(I,J))) / (DXYP(J)-FAW(I,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) - FXSI(I,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)=0, VSIDT(J)>0.
  230 RSI(I,J)  =  RSI(I,J) -  FASI(I,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I,J)*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I,J)*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)<0.
  240 IF(VSIDT(I,J))  260,250,270
C**** VSIDT(J-1)<0, VSIDT(J)=0.
  250 RSI(I,J)  =  RSI(I,J) +  FASI(I,J-1)*BYFOA(I,J)
      TMP = (1d0+FAW(I,J-1)*FOCEAN(I,J-1)*BYFOA(I,J))
      RSIX(I,J) = RSIX(I,J)*TMP
      RSIY(I,J) = RSIY(I,J)*TMP**2
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)<0  or  VSIDT(J-1)>0, VSIDT(J) not 0.
  260 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + ( FASI(I,J-1)- FASI(I,J))
      DO 265 K=1,NTRICE
  265 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSJ(I,K,J-1)-FMSJ(I,K,J))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+
     *     (FYSI(I,J-1)-FYSI(I,J)) + 3d0*((FAW(I,J-1)+
     *     FAW(I,J))*ASI-DXYP(J)*(FASI(I,J-1)+FASI(I,J))))
     *    / (DXYP(J) + (FAW(I,J-1)-FAW(I,J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + (FXSI(I,J-1)-FXSI(I,J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)>0.
  270 RSI(I,J)  =  RSI(I,J) + (FASI(I,J-1)-FASI(I,J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(I,J-1)*FOCEAN(I,J-1)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(I,J-1)*FOCEAN(I,J-1)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))**2
      GO TO 310
C**** VSIDT(J-1)>0.
  280 IF(VSIDT(I,J).ne.0.)  GO TO 260
C**** VSIDT(J-1)>0, VSIDT(J)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(I,J-1)
      DO 285 K=1,NTRICE
  285 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSJ(I,K,J-1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 320
      YRSI = (RSIY(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FYSI(I,J-1)
     *    + 3d0*(FAW(I,J-1)*ASI-DXYP(J)*FASI(I,J-1))) /
     *     (DXYP(J)+FAW(I,J-1))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIY(I,J) = YRSI*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J) + FXSI(I,J-1)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  310 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 330
C**** Sea ice crunches into itself and completely covers grid box
  320 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYFOA(I,J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(I,J) -1d0 / ASI )
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C**** End of loop over J
  330 CONTINUE
C**** End of loop over I
      END DO
C****
C**** Advection of Sea Ice leaving and entering North Pole box
C****
      IF (HAVE_NORTH_POLE) THEN
        ASI = RSI(1,JM)*DXYP(JM)*FOCEAN(1,JM) + SFASI/IM
        DO 345 K=1,NTRICE
  345   AMSI(K) = RSI(1,JM)*DXYP(JM)*MHS(K,1,JM)*FOCEAN(1,JM)+
     &            SFMSI(K)/IM
        IF(ASI.gt.DXYP(JM)*FOCEAN(1,JM))  GO TO 350
        RSI(1,JM)   = ASI*BYFOA(1,JM)
        IF (ASI.gt.0) MHS(1:NTRICE,1,JM) = AMSI(1:NTRICE)/ASI
        GO TO 400
C**** Sea ice crunches into itself at North Pole box
  350   RSI(1,JM)   = 1d0
        MHS(1,1,JM) = AMSI(1)/ASI
        MHS(2,1,JM) =(AMSI(1)+AMSI(2))*BYFOA(1,JM)-MHS(1,1,JM)
        DO K=1,(NTRICE-2)/LMI
          MHS(3+LMI*(K-1),1,JM) = AMSI(3+LMI*(K-1)) / ASI
          MHS(4+LMI*(K-1),1,JM) = AMSI(4+LMI*(K-1)) / ASI
          DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+
     &              AMSI(5+LMI*(K-1))
     *         +AMSI(6+LMI*(K-1)))*(BYFOA(1,JM) -1d0/ ASI)
          MHS(5+LMI*(K-1),1,JM) = AMSI(5+LMI*(K-1)) / ASI +
     *         XSI(3)*DMHSI
          MHS(6+LMI*(K-1),1,JM) = AMSI(6+LMI*(K-1)) / ASI +
     *         XSI(4)*DMHSI
        END DO
      END IF   !HAVE_NORTH_POLE
C****
C**** East-West Advection of Sea Ice
C****
  400 DO 640 J=J_0S, J_1S
C****
C**** Calculate west-east sea ice fluxes at grid box edges
C****
      I=IM
      DO IP1=1,IM
      IF(USIDT(I,J).eq.0.)  GO TO 420
      FAW(I,J) = USIDT(I,J)*DYP(J)
      IF(USIDT(I,J).le.0.) THEN
C**** Sea ice velocity is westward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(IP1,J)-
     *       (1d0+FAW(I,J)*BYDXYP(J))*RSIX(IP1,J))
     *       *FOCEAN(IP1,J)
        FXSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J)*
     *       FAW(I,J)*RSIX(IP1,J)*FOCEAN(IP1,J)-3d0*FASI(I,J))
        FYSI(I,J)=FAW(I,J)*RSIY(IP1,J)*FOCEAN(IP1,J)
        FMSI(1:NTRICE,I) = FASI(I,J)*MHS(1:NTRICE,IP1,J)
      ELSE
C**** Sea ice velocity is eastward at grid box edge
        FASI(I,J)=FAW(I,J)*(RSI(I,J)+(1d0-FAW(I,J)*BYDXYP(J))*RSIX(I,J))
     *       *FOCEAN(I,J)
        FXSI(I,J)=FAW(I,J)*(FAW(I,J)*BYDXYP(J)*FAW(I,J)*
     *       RSIX(I,J)*FOCEAN(I,J)-3d0*FASI(I,J))
        FYSI(I,J)=FAW(I,J)*RSIY(I,J)*FOCEAN(I,J)
        FMSI(1:NTRICE,I) = FASI(I,J)*MHS(1:NTRICE,I,J)
      END IF
         ICIJ(I,J,IJ_MUSI)=ICIJ(I,J,IJ_MUSI)+SUM(FMSI(1:2,I))
         ICIJ(I,J,IJ_HUSI)=ICIJ(I,J,IJ_HUSI)+SUM(FMSI(3:2+LMI,I))
         ICIJ(I,J,IJ_SUSI)=ICIJ(I,J,IJ_SUSI)+SUM(FMSI(3+LMI:2+2*LMI,I))
#ifdef TRACERS_WATER
         DO ITR=1,NTM
           TICIJ(I,J,TICIJ_TUSI,ITR)=TICIJ(I,J,TICIJ_TUSI,ITR)+
     *          SUM(FMSI(3+(1+ITR)*LMI:2+(2+ITR)*LMI,I))
         END DO
#endif
  420 I=IP1
      END DO
C****
C**** Update sea ice variables due to west-east fluxes
C****
      IM1=IM
      DO 630 I=1,IM
      IF(USIDT(IM1,J)) 540,510,580
C**** USIDT(IM1)=0.
  510 IF(USIDT(I,J))  520,630,530
C**** USIDT(IM1)=0, USIDT(I)<0.
  520 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) -  FASI(I,J)
      DO 525 K=1,NTRICE
  525 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) - FMSI(K,I)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) - FXSI(I,J)
     *    + 3d0*(FAW(I,J)*ASI-DXYP(J)*FASI(I,J))) / (DXYP(J)-FAW(I,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) - FYSI(I,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)=0, USIDT(I)>0.
  530 RSI(I,J)  =  RSI(I,J) -  FASI(I,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I,J)*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I,J)*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)<0.
  540 IF(USIDT(I,J))  560,550,570
C**** USIDT(IM1)<0, USIDT(I)=0.
  550 RSI(I,J)  =  RSI(I,J) +  FASI(IM1,J)*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(IM1,J)*FOCEAN(IM1,J)*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(IM1,J)*FOCEAN(IM1,J)*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)<0  or  USIDT(IM1)>0, USIDT(I) not 0.
  560 ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + (FASI(IM1,J)- FASI(I,J))
      DO 565 K=1,NTRICE
  565 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) +
     *       (FMSI(K,IM1)-FMSI(K,I))
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J)+
     *     (FXSI(IM1,J)-FXSI(I,J)) + 3d0*((FAW(IM1,J)+FAW(I,J))*ASI-
     *     DXYP(J)*(FASI(IM1,J)+FASI(I,J))))
     *    / (DXYP(J) + (FAW(IM1,J)-FAW(I,J)))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + (FYSI(IM1,J)-FYSI(I,J))*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)>0.
  570 RSI(I,J)  =  RSI(I,J) + (FASI(IM1,J)-FASI(I,J))*BYFOA(I,J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(IM1,J)*FOCEAN(IM1,J)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(IM1,J)*FOCEAN(IM1,J)
     *     -FAW(I,J)*FOCEAN(I,J))*BYFOA(I,J))
      GO TO 610
C**** USIDT(IM1)>0.
  580 IF(USIDT(I,J).ne.0.)  GO TO 560
C**** USIDT(IM1)>0, USIDT(I)=0.
      ASI = RSI(I,J)*DXYP(J)*FOCEAN(I,J) + FASI(IM1,J)
      DO 585 K=1,NTRICE
  585 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J)*FOCEAN(I,J) + FMSI(K,IM1)
      IF(ASI.gt.DXYP(J)*FOCEAN(I,J))  GO TO 620
      XRSI = (RSIX(I,J)*DXYP(J)*DXYP(J)*FOCEAN(I,J) + FXSI(IM1,J)
     *    + 3d0*(FAW(IM1,J)*ASI-DXYP(J)*FASI(IM1,J))) /
     *     (DXYP(J)+FAW(IM1,J))
      RSI(I,J)  = ASI*BYFOA(I,J)
      RSIX(I,J) = XRSI*BYFOA(I,J)
      RSIY(I,J) = RSIY(I,J) + FYSI(IM1,J)*BYFOA(I,J)
      IF (ASI.gt.0) MHS(1:NTRICE,I,J) = AMSI(1:NTRICE)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  610 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 630
C**** Sea ice crunches into itself and completely covers grid box
  620 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYFOA(I,J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYFOA(I,J) -1d0/ ASI)
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI(4)*DMHSI
      END DO
C**** End of loop over I
  630 IM1=I
C**** End of loop over J
  640 CONTINUE

      IF (KOCEAN.ge.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=J_0, J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
C**** Fresh water sea ice mass convergence (needed for qflux model)
            MSICNV(I,J) = RSI(I,J)*(MHS(1,I,J)+MHS(2,I,J)-SUM(MHS(3
     *           +LMI:2*LMI+2,I,J))) - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     *           +MSI(I,J)-SUM(SSI(1:LMI,I,J)))
C**** sea ice prognostic variables
            SNOWI(I,J)= MAX(0d0,MHS(1,I,J) - ACE1I)
            MSI(I,J)  = MHS(2,I,J)
            DO L=1,LMI
              HSI(L,I,J) = MHS(L+2,I,J)
            END DO
C**** ensure that salinity is only associated with ice
            SICE=MHS(1+2+LMI,I,J)+MHS(2+2+LMI,I,J)
            IF (SNOWI(I,J).gt.XSI(2)*(ACE1I+SNOWI(I,J))) THEN
              SSI(1,I,J)=0.
            ELSE
              SSI(1,I,J)=SICE*(XSI(1)*ACE1I-XSI(2)*SNOWI(I,J))/ACE1I
            END IF
            SSI(2,I,J)=SICE-SSI(1,I,J)
C**** correction of heat energy to compensate for salinity fix
            TICE=Ti(HSI(1,I,J)/(XSI(1)*(ACE1I+SNOWI(I,J))),1d3*MHS(1+2
     *           +LMI,I,J)/(XSI(1)*(ACE1I+SNOWI(I,J))))
            ENRG=XSI(1)*(ACE1I+SNOWI(I,J))*(
     *         Ei(TICE,1d3*MHS(1+2+LMI,I,J)/(XSI(1)*(ACE1I+SNOWI(I,J))))
     *        -Ei(TICE,1d3*SSI(1,I,J)/(XSI(1)*(ACE1I+SNOWI(I,J)))) )
            HSI(1,I,J)=HSI(1,I,J)-ENRG
            HSI(2,I,J)=HSI(2,I,J)+ENRG
C****
            DO L=3,LMI
               SSI(L,I,J) = MHS(L+2+LMI,I,J)
            END DO
#ifdef TRACERS_WATER
C**** reconstruct tracer arrays
            DO ITR=1,NTM
              IF (ACE1I.gt.XSI(2)*(SNOWI(I,J)+ACE1I)) THEN
                TRSI(ITR,1,I,J)= MHS(2+2+(1+ITR)*LMI,I,J) *(ACE1I
     *               -XSI(2)*(SNOWI(I,J)+ACE1I))/ACE1I +MHS(1+2+(1+ITR)
     *               *LMI,I,J)
              ELSE
                TRSI(ITR,1,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)*XSI(1)*(ACE1I
     *               +SNOWI(I,J))/SNOWI(I,J)
              END IF
              TRSI(ITR,2,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)+MHS(2+2+(1+ITR)
     *             *LMI,I,J)-TRSI(ITR,1,I,J)
              DO L=3,LMI
                TRSI(ITR,L,I,J)=MHS(L+2+(1+ITR)*LMI,I,J)
              END DO
            END DO
#endif
            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     *           SUM(SSI(1:LMI,I,J)))
            END IF
          END DO
        END DO
C**** Set atmospheric arrays
        DO J=J_0, J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
C**** set total atmopsheric pressure anomaly in case needed by ocean
              APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+RSI(I,J)
     *             *(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV
              GTEMP(1,2,I,J)=Ti(HSI(1,I,J)/(XSI(1)*(SNOWI(I,J)+ACE1I))
     *             ,1d3*SSI(1,I,J)/(XSI(1)*(SNOWI(I,J)+ACE1I)))
              GTEMP(2,2,I,J)=Ti(HSI(2,I,J)/(XSI(2)*(SNOWI(I,J)+ACE1I))
     *             ,1d3*SSI(2,I,J)/(XSI(2)*(SNOWI(I,J)+ACE1I)))
              GTEMPR(2,I,J) = GTEMP(1,2,I,J)+TF
#ifdef TRACERS_WATER
              GTRACER(:,2,I,J)=TRSI(:,1,I,J)/(XSI(1)*MHS(1,I,J)
     *             -SSI(1,I,J))
#endif
C**** adjust rad fluxes for change in ice fraction
              if (rsi(i,j).gt.rsisave(i,j)) ! ice from ocean
     *       call RESET_SURF_FLUXES(I,J,1,2,RSISAVE(I,J),RSI(I,J))
              if (rsi(i,j).lt.rsisave(i,j)) ! ocean from ice
     *       call RESET_SURF_FLUXES(I,J,2,1,1.-RSISAVE(I,J),1.-RSI(I,J))
C****
            END IF
          END DO
        END DO
        IF (HAVE_NORTH_POLE) THEN
          IF (FOCEAN(1,JM).gt.0) THEN
            DO I=2,IM           ! North pole
              APRESS(I,JM)=APRESS(1,JM)
              GTEMP(1:2,2,I,JM)= GTEMP(1:2,2,1,JM)
              GTEMPR(2,I,JM)   = GTEMPR(2,1,JM)
#ifdef TRACERS_WATER
              GTRACER(:,2,I,JM)=GTRACER(:,2,1,JM)
#endif
            END DO
          END IF
        END IF
        IF (HAVE_SOUTH_POLE) THEN
          IF (FOCEAN(1,1).gt.0) THEN
            DO I=2,IM           ! North pole
              APRESS(I,1)=APRESS(1,1)
              GTEMP(1:2,2,I,1)= GTEMP(1:2,2,1,1)
              GTEMPR(2,I,1)   = GTEMPR(2,1,1)
#ifdef TRACERS_WATER
              GTRACER(:,2,I,1)=GTRACER(:,2,1,1)
#endif
            END DO
          END IF
        END IF
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=J_0, J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
              OA(I,J,13)=OA(I,J,13)+(RSI(I,J)*SUM(MHS(3:2+LMI,I,J))
     *             -RSISAVE(I,J)*SUM(HSI(1:LMI,I,J)))
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****
      RETURN
      END SUBROUTINE ADVSI
#endif

      subroutine INT_AtmA2IceA(aA,iA)
!@sum interpolate from Atm A-grid to Ice A-grid for either scalars
!@+   or U or V component of vector  
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,aGET=>GET
      USE DOMAIN_DECOMP_1D, only : iGET=>GET,ICE_UNPACK=>UNPACK_DATA
      USE ICEDYN, only : IMICDYN,JMICDYN,grid_ICDYN
      IMPLICIT NONE
      real*8 ::
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO),     
     &     iA(1:IMICDYN,
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO)     
#ifdef CUBE_GRID
      real*8, allocatable :: iA_glob(:,:)
      allocate (iA_glob(IMICDYN,JMICDYN))
!interpolate then sumxpe global array on target latlon icedyn grid then scatter
      call parallel_bilin_CS_A_2_latlon_A(agrid,grid_ICDYN,aA,iA_glob,
     &     IMICDYN,JMICDYN)
      call ICE_UNPACK(grid_ICDYN,iA_glob,iA)
      deallocate (iA_glob)
#else 
c***  for the moment me assume that the atm and icedyn grids are 
c***  both latlon with equal resolution
      iA=aA
#endif
      end subroutine INT_AtmA2IceA
c*

      subroutine INT_IceC2AtmC_U(iA,aA)
!@sum interpolate from Ice C-grid to Atm C-grid for U component of vector
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,aGET=>GET
      USE DOMAIN_DECOMP_1D, only : iGET=>GET
      USE ICEDYN, only : IMICDYN,grid_ICDYN
      IMPLICIT NONE
      real *8 ::
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO),     
     &     iA(1:IMICDYN,grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO)     

#ifdef CUBE_GRID
!sumxpe global array on target latlon icedyn grid then scatter

#else 
c***  for the moment me assume that the atm and icedyn grids are 
c***  both latlon with equal resolution
      aA=iA
#endif
      end subroutine INT_IceC2AtmC_U
c*

      subroutine INT_IceC2AtmC_V(iA,aA)
!@sum interpolate from ICe C-grid to Atm C-grid for V component of vector
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,aGET=>GET
      USE DOMAIN_DECOMP_1D, only : iGET=>GET
      USE ICEDYN, only : IMICDYN,grid_ICDYN
      IMPLICIT NONE
      real *8 ::
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO),     
     &     iA(1:IMICDYN,grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO)     

#ifdef CUBE_GRID
!sumxpe global array on target latlon icedyn grid then scatter

#else 
c***  for the moment me assume that the atm and icedyn grids are 
c***  both latlon with equal resolution
      aA=iA
#endif
      end subroutine INT_IceC2AtmC_V
c*

      subroutine INT_IceB2AtmB_NX1(iA,aA)
!@sum interpolate from Ice B-grid to Atm B-grid for either U or V component of vector 
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,aGET=>GET
      USE DOMAIN_DECOMP_1D, only : iGET=>GET
      USE ICEDYN, only : NX1,IMICDYN,grid_ICDYN
      IMPLICIT NONE
      real *8 ::
     &     iA(NX1,
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO),
#ifdef CUBE_GRID
     &     aA(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO)  
#else
     &     aA(NX1,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO)
#endif

#ifdef CUBE_GRID
!sumxpe global array on target latlon icedyn grid then scatter

#else 
c***  for the moment me assume that the atm and icedyn grids are 
c***  both latlon with equal resolution
      aA=iA
#endif
      end subroutine INT_IceB2AtmB_NX1
c*

      subroutine INT_IceB2AtmA(iAb,aAa)
!@sum  interpolation from Ice B-grid to Atm A-grid for either U or V component of vector 
      USE RESOLUTION, only : aIM=>IM,aJM=>JM
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid,aGET=>GET
      USE ICEDYN, only : IMICDYN,grid_ICDYN
      IMPLICIT NONE
      real *8 ::
     &     aAa(agrid%I_STRT_HALO:agrid%I_STOP_HALO,   ! on atm A grid
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO),     
     &     iAb(1:IMICDYN,                             ! on ice B grid
     &     grid_ICDYN%J_STRT_HALO:grid_ICDYN%J_STOP_HALO)
      INTEGER :: aI_1,aI_0,aI_1H,aI_0H,aJ_1,aJ_0,aJ_1H,aJ_0H
 
#ifdef CUBE_GRID
      call aGET(agrid    , I_STRT=aI_0        , I_STOP=aI_1     
     &                   , I_STRT_HALO=aI_0H  , I_STOP_HALO=aI_1H
     &                   , J_STRT=aJ_0        , J_STOP=aJ_1    
     &                   , J_STRT_HALO=aJ_0H  , J_STOP_HALO=aJ_1H )
      
c**** TODO
#endif
      end subroutine INT_IceB2AtmA


      SUBROUTINE init_icedyn(iniOCEAN)
!@sum  init_icedyn initializes ice dynamics variables
!@auth Gavin Schmidt
      USE MODEL_COM, only : im,jm,dtsrc,foceanA=>focean
      USE DIAG_COM, only : ia_src
      USE ICEDYN_COM
      USE ICEDYN, only : focean,osurf_tilt,bydts,usi,vsi
      USE ICEDYN, only : grid_ICDYN,GEOMICDYN
      USE FLUXES, only : uisurf,visurf
      USE PARAM
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER k,kk
      character(len=10) :: xstr,ystr

C**** setup ice dynamics grid
C**** Currently using ATM grid:
C**** -------------------------
C**** If a different grid (but still lat/lon) is required, edit
C**** definition of focean here and the definition of imic,jmic in
C**** ICEDYN.f, and obviously the code in DYNSI.
C**** -------------------------
C**** Set land masks for ice dynamics
      focean(:,:)=foceanA(:,:)         ! EDIT FOR GRID CHANGE

C**** Set up ice momentum grid geometry
      call GEOMICDYN()

      bydts = 1./dtsrc

C**** Initialise ice dynamics if ocean model needs initialising
      if (iniOCEAN) THEN
        RSIX=0.
        RSIY=0.
        USI=0.
        VSI=0.
      end if

C**** set uisurf,visurf for atmospheric drag calculations
      call get_uisurf(usi,vsi,uisurf,visurf)

C**** set properties for ICIJ diagnostics
      do k=1,kicij
        ia_icij(k)=ia_src
        denom_icij(k)=0
        igrid_icij(k)=1
        jgrid_icij(k)=1
      enddo
      k=0

      k=k+1
      IJ_RSI=k
      lname_icij(k)="Ocean ice fraction (ice dynamic grid)"
      sname_icij(k)="icij_rsi"
      units_icij(k)="%"
      scale_icij(k)=100.

      k=k+1
      IJ_USI=k
c      denom_icij(k)=IJ_RSIU ! need to add IJ_RSIU
      lname_icij(k)="Sea ice EW velocity x POICEU"
      sname_icij(k)="icij_usi"
      units_icij(k)="m/s"
      scale_icij(k)=1.
      igrid_icij(k)=2

      k=k+1
      IJ_VSI=k
c      denom_icij(k)=IJ_RSIV ! need to add IJ_RSIV
      lname_icij(k)="Sea ice NS velocity x POICEV"
      sname_icij(k)="icij_vsi"
      units_icij(k)="m/s"
      scale_icij(k)=1.
      jgrid_icij(k)=2

      k=k+1
      IJ_DMUI=k
      lname_icij(k)="Ice-ocean EW stress"
      sname_icij(k)="icij_dmui"
      units_icij(k)="kg/m s^2"
      scale_icij(k)=1./dtsrc

      k=k+1
      IJ_DMVI=k
      lname_icij(k)="Ice-ocean NS stress"
      sname_icij(k)="icij_dmvi"
      units_icij(k)="kg/m s^2"
      scale_icij(k)=1./dtsrc

      k=k+1
      IJ_PICE=k
      denom_icij(k)=IJ_RSI
      lname_icij(k)="Sea ice internal pressure"
      sname_icij(k)="icij_psi"
      units_icij(k)="10^3 kg/m s^2"
      scale_icij(k)=1d-3

      k=k+1
      IJ_MUSI=k
      lname_icij(k)="Sea ice EW mass flux"
      sname_icij(k)="icij_musi"
      units_icij(k)="10^7 kg/s"
      scale_icij(k)=1d-7/dtsrc
      igrid_icij(k)=2

      k=k+1
      IJ_MVSI=k
      lname_icij(k)="Sea ice NS mass flux"
      sname_icij(k)="icij_mvsi"
      units_icij(k)="10^7 kg/s"
      scale_icij(k)=1d-7/dtsrc
      jgrid_icij(k)=2

      k=k+1
      IJ_HUSI=k
      lname_icij(k)="Sea ice EW heat flux"
      sname_icij(k)="icij_husi"
      units_icij(k)="10^12 W"
      scale_icij(k)=1d-12/dtsrc
      igrid_icij(k)=2

      k=k+1
      IJ_HVSI=k
      lname_icij(k)="Sea ice NS heat flux"
      sname_icij(k)="icij_hvsi"
      units_icij(k)="10^12 W"
      scale_icij(k)=1d-12/dtsrc
      jgrid_icij(k)=2

      k=k+1
      IJ_SUSI=k
      lname_icij(k)="Sea ice EW salt flux"
      sname_icij(k)="icij_susi"
      units_icij(k)="10^3 kg/s"
      scale_icij(k)=1d-3/dtsrc
      igrid_icij(k)=2

      k=k+1
      IJ_SVSI=k
      lname_icij(k)="Sea ice NS salt flux"
      sname_icij(k)="icij_svsi"
      units_icij(k)="10^3 kg/s"
      scale_icij(k)=1d-3/dtsrc
      jgrid_icij(k)=2

      if (k.gt.KICIJ) then
        write(6,*) "Too many ICIJ diags: increase KICIJ to at least"
     *       ,k
        call stop_model("ICIJ diagnostic error",255)
      end if

#ifdef TRACERS_WATER
C**** simple tracer diags same description for all tracers
C**** set properties for TICIJ diagnostics
      k=0

      k=k+1
      TICIJ_TUSI=k
      lname_ticij(k)="Sea ice NS tracer flux"
      sname_ticij(k)="ticij_tusi"
      units_ticij(k)="kg/s"
      ia_ticij(k)=ia_src
      scale_ticij(k)=1./dtsrc
      ijgrid_ticij(k)=2

      k=k+1
      TICIJ_TVSI=k
      lname_ticij(k)="Sea ice EW tracer flux"
      sname_ticij(k)="ticij_tvsi"
      units_ticij(k)="kg/s"
      ia_ticij(k)=ia_src
      scale_ticij(k)=1./dtsrc
      ijgrid_ticij(k)=2

      if (k.gt.KTICIJ) then
        write(6,*) "Too many TICIJ diags: increase KTICIJ to at least"
     *       ,k
        call stop_model("TICIJ diagnostic error",255)
      end if
#endif

#ifdef NEW_IO
c
c Declare the dimensions and metadata of output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c

      cdl_icij = ''
      cdl_icij(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_icij(3),'(a,i3,a)') '   lonic = ',imic,' ;'
      write(cdl_icij(4),'(a,i3,a)') '   latic = ',jmic,' ;'
      write(cdl_icij(5),'(a,i3,a)') '   lonic2 = ',imic,' ;'
      write(cdl_icij(6),'(a,i3,a)') '   latic2 = ',jmic,' ;'
      cdl_icij(7:15)(:) = (/
     &     'variables:                       ',
     &     'float lonic(lonic) ;               ',
     &     '   lonic:units = "degrees_east" ; ',
     &     'float latic(latic) ;               ',
     &     '   latic:units = "degrees_north" ;',
     &     'float lonic2(lonic2) ;             ',
     &     '   lonic2:units = "degrees_east" ;',
     &     'float latic2(latic2) ;             ',
     &     '   latic2:units = "degrees_north" ;'
     &     /)
      kk = count(len_trim(cdl_icij).gt.0)
      do k=1,kicij
        if(trim(sname_icij(k)).eq.'unused') cycle
        xstr='lonic) ;'
        if(igrid_icij(k).eq.2) xstr='lonic2) ;'
        ystr='(latic,'
        if(jgrid_icij(k).eq.2) ystr='(latic2,'
        kk = kk + 1
        cdl_icij(kk) = 'float '//trim(sname_icij(k))//
     &       trim(ystr)//trim(xstr)
        kk = kk + 1
        cdl_icij(kk) = '   '//trim(sname_icij(k))//':long_name = "'//
     &       trim(lname_icij(k))//'" ;'
        kk = kk + 1
        cdl_icij(kk) = '   '//trim(sname_icij(k))//':units = "'//
     &       trim(units_icij(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_icij(kk) = '}'
#endif

      RETURN
      END SUBROUTINE init_icedyn

      SUBROUTINE diag_ICEDYN
!@sum  diag_ICEDYN prints out diagnostics for ice dynamics
!@&    ESMF: It should only be called from a serial region.
!@$          It is NOT parallelized
!@auth Gavin Schmidt
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm,trname,ntrocn
#endif
      USE DIAG_COM, only : qdiag,acc_period,
     &     lname_strlen,sname_strlen,units_strlen
      USE ICEDYN_COM
      USE ICEDYN, only : focean
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE
      REAL*8, DIMENSION(IMIC,JMIC) :: Q,ADENOM
      INTEGER I,J,L,N,KXLB,IP1,k1,k
      CHARACTER XLB*30
      CHARACTER TITLE*80
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units
      character*50 :: unit_string
      REAL*8 QJ(JM),QSUM
      REAL*8 byiacc

      IF (.not. QDIAG) RETURN
C**** determine label to be added to all titles
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)

C**** Open output files
      call open_ij(trim(acc_period)//'.icij'//
     *     XLABEL(1:LRUNID),imic,jmic)
C****
C**** Simple scaled ICIJ diagnostics
C****
      DO K=1,KICIJ
        byiacc=1./(IDACC(IA_ICIJ(K))+teeny)
        lname=lname_icij(k)
        if(denom_icij(k).gt.0) then
          adenom=icij(1:imic,1:jmic,denom_icij(k)) * byiacc
        else
          adenom=1.
        endif
        k1 = index(lname,' x ')
        if (k1 .gt. 0) then
          if (index(lname,' x POICEU') .gt. 0) then
            do j=1,jmic
            i=imic
            do ip1=1,imic
              adenom(i,j)=0.5*(icij(i,j,ij_rsi)+icij(ip1,j,ij_rsi))
     *             * byiacc
              i=ip1
            end do
            end do
          else if (index(lname,' x POICEV') .gt. 0) then
            do j=1,jmic-1
            do i=1,imic
              adenom(i,j)=0.5*(icij(i,j,ij_rsi)+icij(i,j+1,ij_rsi))
     *             * byiacc
            end do
            end do
          end if
          lname(k1:50) = ' '
        end if

        Q=UNDEF ; QJ=UNDEF ; QSUM=UNDEF
        DO J=1,JMIC
          DO I=1,IMIC
            IF (ADENOM(I,J).gt.0 .and. FOCEAN(I,J).gt.0.5)
     *           Q(I,J)=SCALE_ICIJ(K)*ICIJ(I,J,K)*byiacc/adenom(i,j)
          END DO
        END DO
        Q(2:IMIC,JMIC)=Q(1,JMIC)
        Q(2:IMIC,1)=Q(1,1)
        TITLE=trim(LNAME)//" ("//trim(UNITS_ICIJ(K))//") "
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME_ICIJ(K),LNAME_ICIJ(K),UNITS_ICIJ(K),Q
     *       ,QJ,QSUM,IGRID_ICIJ(K),JGRID_ICIJ(K))

      END DO

#ifdef TRACERS_WATER
C**** simple tracer diags (no need for weighting)
C**** Name and scale are tracer dependent
      DO K=1,KTICIJ
        byiacc=1./(IDACC(IA_TICIJ(K))+teeny)

        DO N=1,NTM
        lname=trim(trname(n))//" "//lname_ticij(k)
        sname=trim(trname(n))//"_"//sname_ticij(k)
        Q=UNDEF
        DO J=1,JMIC
          DO I=1,IMIC
            IF (FOCEAN(I,J).gt.0.5) Q(I,J)=10**(-ntrocn(n))*
     *           SCALE_TICIJ(K)*TICIJ(I,J,K,N)*byiacc
          END DO
        END DO
        Q(2:IMIC,JMIC)=Q(1,JMIC)
        Q(2:IMIC,1)=Q(1,1)
        UNITS=unit_string(ntrocn(n),UNITS_TICIJ(K))
        TITLE=trim(LNAME)//" ("//trim(UNITS)//")"
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,
     &       IJGRID_TICIJ(K),IJGRID_TICIJ(K)) ! assume igrid=jgrid
      END DO
      END DO
#endif
      call close_ij
C****
      RETURN
      END SUBROUTINE diag_ICEDYN

      SUBROUTINE GET_UISURF(aUSI,aVSI,UISURF,VISURF)
!@sum calculate atmos. A grid winds from C grid
!@auth Gavin Schmidt
C**** temporary assumption that ice velocities are on atmos B grid 
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP_ATM, only : get,agrid=>grid,HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : SOUTH,NORTH
      USE GEOM, only : imaxj
#ifndef CUBE_GRID
      use GEOM,only : idjj,idij,rapj,kmaxj,sinip,cosip
#endif
      IMPLICIT NONE
      REAL*8,INTENT(IN),
     &     DIMENSION(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO) :: aUSI,aVSI
      REAL*8,INTENT(OUT), 
     &     DIMENSION(agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &     agrid%J_STRT_HALO:agrid%J_STOP_HALO) :: UISURF,VISURF
      INTEGER :: I_0,I_1,J_0,J_1,I,J,K,IM1
      LOGICAL :: pole
      REAL*8 :: hemi

      CALL GET(agrid, J_STRT=J_0, J_STOP=J_1)
      I_0 = agrid%I_STRT
      I_1 = agrid%I_STOP

      Call HALO_UPDATE(aGRID, uisurf)
      Call HALO_UPDATE(aGRID, visurf)

      POLE=.FALSE.

C**** go from atm C to atm A
      do j=J_0,J_1
        do i=I_0,IMAXJ(J)
         im1=i-1
#ifndef CUBE_GRID
         if (im1 .le. 0) im1=IM
         POLE= (J.EQ.1 .or. J.EQ.JM)

          HEMI=1.
          IF(J.LE.JM/2) HEMI=-1.
C**** Note that usi,vsi start with j=1, (not j=2 as in atm winds)
          if (pole) then
            uisurf(i,j) = 0. ; visurf(i,j) = 0.
            do k=1,kmaxj(j)
              uisurf(i,j) = uisurf(i,j) + rapj(k,j)*(ausi(idij(k,i,j)
     *             ,idjj(k,j)-1)*cosip(k)-hemi*avsi(idij(k,i,j),idjj(k
     *             ,j)-1)*sinip(k))
              visurf(i,j) = visurf(i,j) + rapj(k,j)*(avsi(idij(k,i,j)
     *             ,idjj(k,j)-1)*cosip(k)+hemi*ausi(idij(k,i,j),idjj(k
     *             ,j)-1)*sinip(k))
            end do
          else
#endif
            uisurf(i,j) = 0.5*(ausi(i,j)+ausi(im1,j))
            visurf(i,j) = 0.5*(avsi(i,j)+avsi(i,j-1))
#ifndef CUBE_GRID
          end if
#endif
        end do
      end do

      RETURN
      END SUBROUTINE GET_UISURF
