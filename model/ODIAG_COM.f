#include "rundeck_opts.h"

      MODULE ODIAG
!@sum  ODIAG ocean diagnostic arrays (incl. dynamic sea ice)
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm
#endif
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE DIAG_COM, only : npts  ! needed for conservation diags
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: KOIJ=6,KOIJL=22,KOL=6,KOLNST=8
!@var OIJ   lat-lon ocean diagnostics (on ocean grid)
!@var OIJL  3-dimensional ocean diagnostics
!@var OL    vertical ocean diagnostics
!@var OLNST strait diagnostics
!ny?  logical :: allocated_odiag_glob = .false.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: OIJ_loc   !ny? ,OIJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: OIJL_loc !ny? ,OIJL
      REAL*8, DIMENSION(IM,JM,KOIJ)      :: OIJ
      REAL*8, DIMENSION(IM,JM,LMO,KOIJL) :: OIJL
      REAL*8, DIMENSION(LMO,KOL)   :: OL
      REAL*8, DIMENSION(LMO,NMST,KOLNST):: OLNST
!@var IJ_xxx Names for OIJ diagnostics
      INTEGER IJ_HBL,IJ_BO,IJ_BOSOL,IJ_USTAR,IJ_SSH,IJ_PB
!@var lname_oij Long names for OIJ diagnostics
      CHARACTER*50, DIMENSION(KOIJ) :: LNAME_OIJ
!@var sname_oij Short names for OIJ diagnostics
      CHARACTER*30, DIMENSION(KOIJ) :: SNAME_OIJ
!@var units_oij Units for OIJ diagnostics
      CHARACTER*50, DIMENSION(KOIJ) :: UNITS_OIJ
!@var ia_oij IDACC numbers for OIJ diagnostics
      INTEGER, DIMENSION(KOIJ) :: IA_OIJ
!@var scale_oij scales for OIJ diagnostics
      REAL*8, DIMENSION(KOIJ) :: SCALE_OIJ
!@var ijgrid_oij Grid descriptor for OIJ diagnostics
       INTEGER, DIMENSION(KOIJ) :: IJGRID_OIJ

!@var IJL_xxx Names for OIJL diagnostics
      INTEGER IJL_MO,IJL_G0M,IJL_S0M,IJL_GFLX,IJL_SFLX,IJL_MFU,IJL_MFV
     *     ,IJL_MFW,IJL_GGMFL,IJL_SGMFL,IJL_KVM,IJL_KVG,IJL_WGFL
     *     ,IJL_WSFL
!@var LN_xxx Names for OLNST diagnostics
      INTEGER LN_KVM,LN_KVG,LN_WGFL,LN_WSFL,LN_MFLX,LN_GFLX,LN_SFLX
     *     ,LN_ICFL
!@var L_xxx Names for OL diagnostics
      INTEGER L_RHO,L_TEMP,L_SALT
!@var icon_xx indexes for conservation quantities
      INTEGER icon_OCE,icon_OKE,icon_OAM,icon_OMS,icon_OSL
!@var kbasin integer index of which basin a particular ocean point is in
      INTEGER, DIMENSION(IM,JM) :: KBASIN
!@var XLB label for diagnostic titles
      CHARACTER XLB*30
!@var FLAT latitude values on primary and secondary ocean grids
      REAL*8, DIMENSION(JM,2) :: FLAT
!@var FLON longitude values on primary and secondary ocean grids
      REAL*8, DIMENSION(IM,2) :: FLON
!@var iu_otj unit number for ascii output of ocean transports
      INTEGER iu_otj
!@var BASIN names of ocean basins for diag output
      CHARACTER*16, DIMENSION(4) :: BASIN=
     *     (/"Atlantic","Pacific ","Indian  ","Global  "/)
C****

#ifdef TRACERS_OCEAN
!@var KTOIJL number of 3-dimensional ocean tracer diagnostics
      INTEGER, PARAMETER :: KTOIJL=10
!@var TOIJL  3-dimensional ocean tracer diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TOIJL_loc !ny? ,TOIJL
      REAL*8, DIMENSION(IM,JM,LMO,KTOIJL,NTM)    :: TOIJL
!@var toijl_xxx indices for TOIJL diags
      INTEGER, PARAMETER :: toijl_conc=1,toijl_tflx=2,toijl_gmfl=6
     *     ,toijl_wtfl=10
!@var TLNST strait diagnostics
      REAL*8, DIMENSION(LMO,NMST,KOLNST,NTM):: TLNST
#endif

#ifdef NEW_IO
c declarations that facilitate summation of acc-files when using
c the new i/o system
      target :: oij_loc,oijl_loc,ol,olnst
      real*8, dimension(:,:,:), target, allocatable ::
     &     oij_fromdisk
      real*8, dimension(:,:,:,:), target, allocatable ::
     &     oijl_fromdisk
      real*8, DIMENSION(LMO,KOL), target :: OL_fromdisk
      real*8, DIMENSION(LMO,NMST,KOLNST), target :: OLNST_fromdisk
      real*8, dimension(:,:), pointer :: ol_ioptr
      real*8, dimension(:,:,:), pointer :: oij_ioptr,olnst_ioptr
      real*8, dimension(:,:,:,:), pointer :: oijl_ioptr
#ifdef TRACERS_OCEAN
      target :: toijl_loc,tlnst
      real*8, dimension(:,:,:,:,:), target, allocatable ::
     &     toijl_fromdisk
      real*8, DIMENSION(LMO,NMST,KOLNST,NTM), target:: TLNST_fromdisk
      real*8, dimension(:,:,:,:,:), pointer :: toijl_ioptr
      real*8, dimension(:,:,:,:), pointer :: tlnst_ioptr
#endif
#endif

      END MODULE ODIAG

      SUBROUTINE io_ocdiag(kunit,it,iaction,ioerr)
!@sum  io_ocdiag reads and writes ocean diagnostic arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irsfic,irerun,irsficno,ioread_single,lhead
      USE ODIAG

!      USE DOMAIN_DECOMP_1D, only: grid, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var OIJ4,OIJL4,OLNST4,OL4 dummy arrays for reading diag. files
      logical , save :: de_alloc = .true.
      REAL*4, DIMENSION(IM,JM,KOIJ)  :: OIJ4
      REAL*4, DIMENSION(IM,JM,LMO,KOIJL) :: OIJL4
      REAL*4, DIMENSION(LMO,KOL)   :: OL4
      REAL*4, DIMENSION(LMO,NMST,KOLNST):: OLNST4
#ifdef TRACERS_OCEAN
!@var TOIJL_glob work array for read/wrote operations
      REAL*4, DIMENSION(IM,JM,LMO,KTOIJL,NTM) :: TOIJL4
      REAL*4, DIMENSION(LMO,NMST,KOLNST,NTM)  :: TLNST4
!@var TR_HEADER Character string label for individual tracer records
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TROCDIAG01"

      write(TR_MODULE_HEADER(lhead+1:80),'(a19,i2,a1,i2,a9,i4,a4)')
     *     'R8 Toijl(im,jm,lmo,',ktoijl,',',ntm,'), Tlnst(',
     *     LMO*NMST*KOLNST*NTM,'),it'
#endif
      write(MODULE_HEADER(lhead+1:80),'(a13,i2,a13,i2,a1,  i2,a5,i2,
     *  a1,i2,a8,i4,a)') 'R8 Oij(im,jm,',koij,'),Oijl(im,jm,',lmo,',',
     *  koijl,'),Ol(',lmo,   ',',kol,'),OLNST(',LMO*NMST*KOLNST,'),it'

!ny?  if(.not.allocated_odiag_glob) then
!ny?    call alloc_odiag_glob
!ny?    allocated_odiag_glob = .true.
!ny?  end if

      SELECT CASE (IACTION)
      CASE (IOWRITE)  ! output to standard restart file
        call gather_odiags ()
        IF (AM_I_ROOT())
     *    WRITE (kunit,err=10) MODULE_HEADER,
     *      OIJ,OIJL,OL,OLNST,it
#ifdef TRACERS_OCEAN
        IF (AM_I_ROOT())
     *    WRITE (kunit,err=10) TR_MODULE_HEADER,TOIJL,TLNST,it
#endif
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        call gather_odiags ()
        IF (AM_I_ROOT())
     *    WRITE (kunit,err=10) MODULE_HEADER,REAL(OIJ,KIND=4),
     *      REAL(OIJL,KIND=4),REAL(OL,KIND=4),REAL(OLNST,KIND=4),it
#ifdef TRACERS_OCEAN
        TR_MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        IF (AM_I_ROOT())
     *    WRITE (kunit,err=10) TR_MODULE_HEADER,
     *      REAL(TOIJL,KIND=4),REAL(TLNST,KIND=4),it
#endif
      CASE (IOREAD:)            ! input from acc/restart file
        SELECT CASE (IACTION)
        CASE (IRSFICNO)         ! initial conditions
        CASE (ioread_single)    ! input from acc files (post-processing)
         if  (AM_I_ROOT()) then
            READ (kunit,err=10) HEADER,OIJ4,OIJL4,OL4,OLNST4,it
C**** accumulate diagnostics
            if (de_alloc) then   !  1st time only
                de_alloc = .false. ; OIJ=0. ; OIJL=0.
#ifdef TRACERS_OCEAN
                                   TOIJL=0.
#endif
            end if
            OIJ=OIJ+OIJ4
            OIJL=OIJL+OIJL4
            OL=OL+OL4
            OLNST=OLNST+OLNST4
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
#ifdef TRACERS_OCEAN
            READ (kunit,err=10) TR_HEADER,TOIJL4,TLNST4,it
C**** accumulate diagnostics
            TOIJL=TOIJL+TOIJL4
            TLNST=TLNST+TLNST4
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
              GO TO 10
            END IF
#endif
         end if
         call scatter_odiags ()
        CASE (ioread)    ! restarts
          IF (AM_I_ROOT()) then
            READ (kunit,err=10) HEADER,OIJ,OIJL,OL,OLNST,it
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
#ifdef TRACERS_OCEAN
            READ (kunit,err=10) TR_HEADER,TOIJL,TLNST,it
            IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
              GO TO 10
            END IF
#endif
          END IF
          call scatter_odiags ()
        END SELECT
      END SELECT

!ny?  if(de_alloc) then
!ny?    allocated_odiag_glob = .false.
!ny?    call de_alloc_odiag_glob
!ny?  end if

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocdiag

#ifdef NEW_IO
      subroutine def_rsf_ocdiag(fid,r4_on_disk)
!@sum  def_rsf_ocdiag defines ocean diag array structure in restart+acc files
!@auth M. Kelley
!@ver  beta
      use odiag, only : ol,olnst,oij=>oij_loc,oijl=>oijl_loc
#ifdef TRACERS_OCEAN
      use odiag, only : tlnst,toijl=>toijl_loc
#endif
      USE OCEANR_DIM, only : grid=>ogrid
      use pario, only : defvar
      implicit none
      integer fid            !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,oij,'oij(dist_imo,dist_jmo,koij)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,oijl,
     &     'oijl(dist_imo,dist_jmo,lmo,koijl)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,ol,'ol(lmo,kol)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,olnst,'olnst(lmo,nmst,kolnst)',
     &     r4_on_disk=r4_on_disk)
#ifdef TRACERS_OCEAN
      call defvar(grid,fid,toijl,
     &     'toijl(dist_imo,dist_jmo,lmo,ktoijl,ntm)',
     &     r4_on_disk=r4_on_disk)
      call defvar(grid,fid,tlnst,'tlnst(lmo,nmst,kolnst,ntm)',
     &     r4_on_disk=r4_on_disk)
#endif
      return
      end subroutine def_rsf_ocdiag

      subroutine new_io_ocdiag(fid,iaction)
!@sum  new_io_ocdiag read/write ocean arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      USE OCEANR_DIM, only : grid=>ogrid
c in the postprocessing case where arrays are read from disk and summed,
c these i/o pointers point to temporary arrays.  Otherwise, they point to
c the instances of the arrays used during normal operation. 
      use odiag, only : ol=>ol_ioptr,olnst=>olnst_ioptr,
     &     oij=>oij_ioptr,oijl=>oijl_ioptr
#ifdef TRACERS_OCEAN
      use odiag, only : tlnst=>tlnst_ioptr,toijl=>toijl_ioptr
#endif
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart or acc file
        call write_dist_data(grid,fid,'oij',oij)
        call write_dist_data(grid,fid,'oijl',oijl)
        call write_data(grid,fid,'ol',ol)
c straits arrays
        call write_data(grid,fid,'olnst',olnst)
#ifdef TRACERS_OCEAN
        call write_dist_data(grid,fid,'toijl',toijl)
        call write_data(grid,fid,'tlnst',tlnst)
#endif
      case (ioread)            ! input from restart or acc file
        call read_dist_data(grid,fid,'oij',oij)
        call read_dist_data(grid,fid,'oijl',oijl)
        call read_data(grid,fid,'ol',ol,bcast_all=.true.)
c straits arrays
        call read_data(grid,fid,'olnst',olnst,bcast_all=.true.)
#ifdef TRACERS_OCEAN
        call read_dist_data(grid,fid,'toijl',toijl)
        call read_data(grid,fid,'tlnst',tlnst,bcast_all=.true.)
#endif
      end select
      return
      end subroutine new_io_ocdiag

      subroutine set_ioptrs_ocnacc_default
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays used during normal operation. 
      use odiag
      implicit none
      oij_ioptr    => oij_loc
      oijl_ioptr   => oijl_loc
      ol_ioptr     => ol
      olnst_ioptr  => olnst
#ifdef TRACERS_OCEAN
      toijl_ioptr  => toijl_loc
      tlnst_ioptr  => tlnst
#endif
      return
      end subroutine set_ioptrs_ocnacc_default

      subroutine set_ioptrs_ocnacc_sumfiles
c point i/o pointers for diagnostic accumlations to temporary
c arrays that hold data read from disk
      use odiag
      USE OCEANR_DIM, only : grid=>ogrid
      implicit none
      integer :: j_0h,j_1h
c
c allocate arrays (local size) and set i/o pointers
c
      if(allocated(oij_fromdisk)) return
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      allocate(oij_fromdisk(im,j_0h:j_1h,koij))
      oij_ioptr => oij_fromdisk
      allocate(oijl_fromdisk(im,j_0h:j_1h,lmo,koijl))
      oijl_ioptr => oijl_fromdisk
      ol_ioptr     => ol_fromdisk
      olnst_ioptr  => olnst_fromdisk
#ifdef TRACERS_OCEAN
      allocate(toijl_fromdisk(im,j_0h:j_1h,lmo,ktoijl,ntm))
      toijl_ioptr  => toijl_fromdisk
      tlnst_ioptr  => tlnst_fromdisk
#endif
      return
      end subroutine set_ioptrs_ocnacc_sumfiles

      subroutine sumfiles_ocnacc
c increment diagnostic accumlations with the data that was
c read from disk and stored in the _fromdisk arrays.
      use odiag
      implicit none
      oij_loc    = oij_loc    + oij_fromdisk
      oijl_loc   = oijl_loc   + oijl_fromdisk
      ol         = ol         + ol_fromdisk
      olnst      = olnst      + olnst_fromdisk
#ifdef TRACERS_OCEAN
      toijl_loc  = toijl_loc  + toijl_fromdisk
      tlnst      = tlnst      + tlnst_fromdisk
#endif
      return
      end subroutine sumfiles_ocnacc

#endif /* NEW_IO */

      SUBROUTINE DIAGCO (M)
!@sum  DIAGCO Keeps track of the ocean conservation properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE ODIAG, only : icon_OCE,icon_OKE,icon_OMS,icon_OSL,icon_OAM
      IMPLICIT NONE
!@var M index denoting from where DIAGCO is called (see DIAGCA)
      INTEGER, INTENT(IN) :: M
      REAL*8, EXTERNAL :: conserv_OCE,conserv_OKE,conserv_OMS
     *     ,conserv_OSL,conserv_OAM

C**** OCEAN MASS
      CALL conserv_ODIAG(M,conserv_OMS,icon_OMS)

C**** OCEAN ANGULAR MOMENTUM
      CALL conserv_ODIAG(M,conserv_OAM,icon_OAM)

C**** OCEAN KINETIC ENERGY
      CALL conserv_ODIAG(M,conserv_OKE,icon_OKE)

C**** OCEAN POTENTIAL ENTHALPY
      CALL conserv_ODIAG(M,conserv_OCE,icon_OCE)

C**** OCEAN SALT
      CALL conserv_ODIAG(M,conserv_OSL,icon_OSL)

C****
      RETURN
      END SUBROUTINE DIAGCO

      SUBROUTINE conserv_ODIAG (M,CONSFN,ICON)
!@sum  conserv_ODIAG generic routine keeps track of conserved properties
!@+    (cloned version from AGCM for better regridding)
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : jm
      USE DOMAIN_DECOMP_1D, only : GET, GRID, CHECKSUMj
      USE DIAG_COM, only : consrv=>consrv_loc,nofm
      IMPLICIT NONE
!@var M index denoting from where routine is called
      INTEGER, INTENT(IN) :: M
!@var ICON index for the quantity concerned
      INTEGER, INTENT(IN) :: ICON
!@var CONSFN external routine that calculates total conserved quantity
      EXTERNAL CONSFN
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TOTAL
      INTEGER :: I,J,NM,NI
      INTEGER :: J_0,J_1

      CALL GET(grid, J_STRT=J_0,         J_STOP=J_1)

C**** NOFM contains the indexes of the CONSRV array where each
C**** change is to be stored for each quantity. If NOFM(M,ICON)=0,
C**** no calculation is done.
C**** NOFM(1,ICON) is the index for the instantaneous value.
      IF (NOFM(M,ICON).gt.0) THEN
C**** Calculate current value TOTAL
        CALL CONSFN(TOTAL)
        NM=NOFM(M,ICON)
        NI=NOFM(1,ICON)
C**** Accumulate difference from last time in CONSRV(NM)
        IF (M.GT.1) THEN
          DO J=J_0,J_1
            CONSRV(J,NM)=CONSRV(J,NM)+(TOTAL(J)-CONSRV(J,NI))
          END DO
        END IF
C**** Save current value in CONSRV(NI)
        DO J=J_0,J_1
          CONSRV(J,NI)=TOTAL(J)
        END DO
      END IF
      RETURN
C****
      END SUBROUTINE conserv_ODIAG

      SUBROUTINE init_ODIAG
!@sum  init_ODIAG initialises ocean diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav
      USE MODEL_COM, only : dtsrc
      USE OCEAN, only : ze
      USE DIAG_COM, only : ia_src,conpt0,zoc,zoc1
      USE ODIAG
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T = .TRUE. , F = .FALSE.
      CHARACTER CONPT(NPTS)*10
      INTEGER k

C**** Set names for OIJL diagnostics
      IJL_MO=1  ; IJL_G0M=2  ; IJL_S0M=3   ; IJL_MFU=4    ; IJL_MFV=5
      IJL_MFW=6 ; IJL_KVM=13; IJL_KVG=14 ; IJL_WGFL=15 ; IJL_WSFL=16
C**** These flux diagnostics need 3 spots each
      IJL_GFLX=7 ; IJL_SFLX=10 ; IJL_GGMFL=17 ; IJL_SGMFL=20
C**** Set names for OLNST diagnostics
      LN_KVM=1  ; LN_KVG=2  ; LN_WGFL=3 ; LN_WSFL=4
      LN_MFLX=5 ; LN_GFLX=6 ; LN_SFLX=7 ; LN_ICFL=8
C**** Set names for OL diagnostics
      L_RHO=1   ; L_TEMP=2  ; L_SALT=3

C**** set properties for OIJ diagnostics
      k=0
c
      k=k+1
      IJ_HBL=k
      lname_oij(k)="Ocean Boundary layer depth (KPP) x PO4"
      sname_oij(k)="oij_hbl"
      units_oij(k)="m"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25
      ijgrid_oij(k)=1

      k=k+1
      IJ_BO=k
      lname_oij(k)="Surface buoyancy forcing (KPP) x PO4"
      sname_oij(k)="oij_bo"
      units_oij(k)="10^-7 m^2/s^3"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25*1d7
      ijgrid_oij(k)=1

      k=k+1
      IJ_BOSOL=k
      lname_oij(k)="Surface solar buoyancy flux x PO4"
      sname_oij(k)="oij_bosol"
      units_oij(k)="10^-7 m^2/s^3"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25*1d7
      ijgrid_oij(k)=1

      k=k+1
      IJ_USTAR=k
      lname_oij(k)="Surface friction speed x PO4"
      sname_oij(k)="oij_ustar"
      units_oij(k)="m/s"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25
      ijgrid_oij(k)=1

      k=k+1
      IJ_SSH=k
      lname_oij(k)="Ocean surface height"
      sname_oij(k)="oij_ssh"
      units_oij(k)="m"
      ia_oij(k)=ia_src
      scale_oij(k)=bygrav
      ijgrid_oij(k)=1

      k=k+1
      IJ_PB=k
      lname_oij(k)="Ocean bottom pressure anomaly"
      sname_oij(k)="oij_pb"
      units_oij(k)="Pa"
      ia_oij(k)=ia_src
      scale_oij(k)=1.
      ijgrid_oij(k)=1

      if (k.gt.KOIJ) then
        write(6,*) "Too many OIJ diagnostics: increase KOIJ to at least"
     *       ,k
        call stop_model("OIJ diagnostic error",255)
      end if

C**** Set up oceanic component conservation diagnostics
C**** Oceanic mass
      CONPT=CONPT0
      CONPT(8)="OCN PHYS"
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN MASS","(10**2 KG/M^2) ",
     *     "(10**-8 KG/SM^2)",1d-2,1d8,icon_OMS)
C**** Oceanic angular momentum
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN AM  ","(10**12 JS/M^2)",
     *     "(10**2 J/M^2)   ",1d-12,1d-2,icon_OAM)
C**** Oceanic kinetic energy
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCEAN KE","(J/M^2)        ",
     *     "(10**-6 W/M^2)  ",1d0,1d6,icon_OKE)
C**** Oceanic potential enthalpy (heat)
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN HEAT","(10**6 J/M^2)  ",
     *     "(10**-2 W/M^2)  ",1d-6,1d2,icon_OCE)
C**** Oceanic salt mass
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN SALT","(10 KG/M^2)    ",
     *     "(10**-9 KG/SM^2)",1d-1,1d9,icon_OSL)
C**** Initialise ocean basins
      CALL OBASIN

C**** Define ocean depths for diagnostic output
      ZOC1(1:LMO+1) = ZE(0:LMO)
      ZOC(1:LMO) = 0.5*(ZE(1:LMO)+ZE(0:LMO-1))

C****
      RETURN
      END SUBROUTINE init_ODIAG

      SUBROUTINE reset_odiag(isum)
!@sum  reset_odiag zeros out ocean diagnostics if needed
!@auth G. Schmidt
!@ver  1.0
      USE ODIAG, only : oij=>oij_loc,oijl=>oijl_loc,ol,olnst
#ifdef TRACERS_OCEAN
     *     ,toijl=>toijl_loc,tlnst
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: isum  ! needed for plug-play compatibility

      OIJ=0. ; OIJL=0. ; OL=0. ; OLNST=0.

#ifdef TRACERS_OCEAN
      TOIJL=0. ; TLNST = 0.
#endif

      return
      END SUBROUTINE reset_odiag

      SUBROUTINE alloc_odiag(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Reto Ruedy
!@ver  1.0

      USE DOMAIN_DECOMP_1D, only : dist_grid,get

      USE ODIAG

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(        OIJ_loc (IM,J_0H:J_1H,KOIJ), STAT=IER )
      ALLOCATE(       OIJL_loc (IM,J_0H:J_1H,LMO,KOIJL), STAT=IER )
#ifdef TRACERS_OCEAN
      ALLOCATE(      TOIJL_loc (IM,J_0H:J_1H,LMO,KTOIJL,NTM), STAT=IER )
#endif

      END SUBROUTINE alloc_odiag

!ny?  SUBROUTINE alloc_odiag_glob    ! not yet in use
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Reto Ruedy
!@ver  1.0

!ny?  USE ODIAG

!ny?  IMPLICIT NONE
!ny?  integer ier

!ny?  ALLOCATE(       OIJ  (IM,JM,KOIJ), STAT=IER )
!ny?  ALLOCATE(       OIJL (IM,JM,LMO,KOIJL), STAT=IER )
#ifdef TRACERS_OCEAN
!ny?  ALLOCATE(      TOIJL (IM,JM,LMO,KTOIJL,NTM), STAT=IER )
#endif

!ny?  END SUBROUTINE alloc_odiag_glob

!ny?  SUBROUTINE de_alloc_odiag_glob
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Reto Ruedy
!@ver  1.0

!ny?  USE ODIAG

!ny?  IMPLICIT NONE

!ny?  DEALLOCATE( OIJ, OIJL )
#ifdef TRACERS_OCEAN
!ny?  DEALLOCATE( TOIJL )
#endif

!ny?  END SUBROUTINE de_alloc_odiag_glob

      SUBROUTINE gather_odiags ()
!@sum  collect the local acc-arrays into global arrays
!@+    run-time
!@auth Reto Ruedy
!@ver  1.0

      USE ODIAG

!      use domain_decomp_1d, only : grid, pack_data
      use domain_decomp_1d, only : pack_data
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE

      call pack_data (grid, OIJ_loc  , OIJ)
      call pack_data (grid, OIJL_loc , OIJL)
#ifdef TRACERS_OCEAN
      call pack_data (grid, TOIJL_loc, TOIJL)
#endif

      END SUBROUTINE gather_odiags

      SUBROUTINE scatter_odiags ()
!@sum  To distribute the global acc-arrays to the local pieces
!@auth Reto Ruedy
!@ver  1.0

      USE ODIAG

!      use domain_decomp_1d, only : grid, unpack_data, ESMF_BCAST
      use domain_decomp_1d, only : unpack_data, ESMF_BCAST
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE

      call unpack_data (grid, OIJ  , OIJ_loc)
      call unpack_data (grid, OIJL , OIJL_loc)
      CALL ESMF_BCAST(grid, OL)
      CALL ESMF_BCAST(grid, OLNST)
#ifdef TRACERS_OCEAN
      call unpack_data (grid, TOIJL, TOIJL_loc)
      CALL ESMF_BCAST(grid, TLNST)
#endif

      END SUBROUTINE scatter_odiags

